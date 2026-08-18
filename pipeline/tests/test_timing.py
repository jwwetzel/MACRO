"""Unit tests for macro_core.timing — the S3 shared time axis.

Layout mirrors the module:

* DATE-OBS parsing on every format the archive actually contains;
* the mid-exposure policy per readout family, including the StackPro
  worst-case bound;
* the BJD_TDB wrapper, anchored three independent ways:
    1. the TDB-UTC scale offset against its literature value
       (32.184 s + 37 leap seconds in the 2017+ era),
    2. a hand-computed heliocentric light-travel time at the 2024
       perihelion (almanac distance x light time per au),
    3. an external-software cross-check against a MaxIm-written JD-HELIO
       header recorded in the immutable archive (coe M16 frame);
* the pure eclipse-fitting helpers on synthetic data.

Every astropy-touching test uses the 'builtin' ephemeris so the suite
stays green offline; one extra test pins builtin against de440s when the
JPL kernel is available locally (skipped, with a reason, when not).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_core import manifest as m           # noqa: E402
from macro_core import timing as tm            # noqa: E402


# ---------------------------------------------------------------------------
# DATE-OBS parsing
# ---------------------------------------------------------------------------
class TestParseDateObs:
    def test_three_digit_fraction(self):
        # 2024-05-10T03:07:54.259 -> JD checked by hand:
        # 03:07:54.259 = 0.130489... day past midnight; JD midnight = .5
        jd = tm.parse_date_obs("2024-05-10T03:07:54.259")
        assert abs(jd - 2460440.6304891088) < 1e-8   # < 1 ms

    def test_two_digit_fraction(self):
        # The Fast/Mode0 era writes 2-digit fractions ('.60'), which
        # Python 3.10 fromisoformat rejects raw — the parser must pad.
        jd = tm.parse_date_obs("2026-03-22T03:53:57.60")
        assert jd is not None
        # .60 means 600 ms, NOT 60 ms: check against the whole-second JD.
        # Tolerance 1 ms: a full JD eats ~47 of the double's 52 mantissa
        # bits, so differences of two JDs carry ~40 us of float noise.
        base = tm.parse_date_obs("2026-03-22T03:53:57")
        assert abs((jd - base) * 86400.0 - 0.60) < 1e-3

    def test_no_fraction(self):
        jd = tm.parse_date_obs("2023-06-01T08:18:08")
        # Midnight JD fraction is 0.5; 08:18:08 = 29888 s.
        assert abs(jd - (2460096.5 + 29888.0 / 86400.0)) < 1e-8

    def test_garbage_and_missing(self):
        assert tm.parse_date_obs(None) is None
        assert tm.parse_date_obs("") is None
        assert tm.parse_date_obs("  ") is None
        assert tm.parse_date_obs("10/05/24") is None
        assert tm.parse_date_obs("not a date") is None

    def test_epoch_constant_matches_s0(self):
        # S0 and S3 must share one epoch or every comparison is garbage.
        assert tm.UNIX_EPOCH_JD == m.UNIX_EPOCH_JD


# ---------------------------------------------------------------------------
# Mid-exposure policy
# ---------------------------------------------------------------------------
class TestMidPolicy:
    def test_plain_frame_half_exptime(self):
        mid, method = tm.jd_utc_mid(2460000.0, 120.0, "Fast")
        assert method == tm.MID_PLAIN
        assert abs((mid - 2460000.0) * 86400.0 - 60.0) < 1e-3

    def test_stackpro_same_arithmetic_distinct_label(self):
        mid, method = tm.jd_utc_mid(2460000.0, 1024.0, "High Gain StackPro")
        assert method == tm.MID_STACKPRO
        assert abs((mid - 2460000.0) * 86400.0 - 512.0) < 1e-3

    def test_stackpro_detection_is_substring_and_case_blind(self):
        assert tm.is_stackpro("High Gain StackPro")
        assert tm.is_stackpro("low gain stackpro")   # future-proofing
        assert not tm.is_stackpro("High Gain")
        assert not tm.is_stackpro(None)
        assert not tm.is_stackpro("")

    def test_no_jd(self):
        mid, method = tm.jd_utc_mid(None, 60.0, "Fast")
        assert mid is None and method == tm.MID_NO_JD
        mid, method = tm.jd_utc_mid(float("nan"), 60.0, "Fast")
        assert mid is None and method == tm.MID_NO_JD

    def test_exptime_nonpos_uses_start(self):
        for bad in (None, 0.0, -1.0, float("nan")):
            mid, method = tm.jd_utc_mid(2460000.25, bad, "Mode0")
            assert method == tm.MID_EXPTIME_NONPOS
            assert mid == 2460000.25

    def test_stackpro_worst_case_bound(self):
        # The policy's stated worst case: half the empirical dead-time
        # bound, and comfortably sub-second (the S3 precision target).
        assert tm.worst_case_mid_error_s("High Gain StackPro") == \
            pytest.approx(tm.STACKPRO_DEADTIME_BOUND_S / 2.0)
        assert tm.worst_case_mid_error_s("High Gain StackPro") < 0.5
        assert tm.worst_case_mid_error_s("Mode0") == 0.0


# ---------------------------------------------------------------------------
# BJD_TDB — three independent anchors
# ---------------------------------------------------------------------------
class TestBjdTdb:
    def test_tdb_minus_utc_literature_value(self):
        """Anchor 1: the UTC->TDB scale offset.

        Literature: TDB - UTC = 32.184 s (TT-TAI) + 37 s (leap seconds
        since 2017-01-01, IERS Bulletin C) + periodic terms < 2 ms.
        A wrong leap-second table or scale confusion shows up here as a
        whole-second error — the classic timing-pipeline failure mode.
        """
        _, _, tdb_minus_utc = tm.bjd_tdb_from_utc(
            2460310.5, 90.0, 66.56, ephemeris="builtin")  # 2024-01-01, NEP
        assert abs(float(tdb_minus_utc) - 69.184) < 0.005

    def test_perihelion_hand_computed_heliocentric(self):
        """Anchor 2: hand-computed light-travel time at 2024 perihelion.

        Almanac facts (USNO/Astronomical Almanac 2024): Earth perihelion
        2024 Jan 3 ~00:38 UT at r = 0.98330 au; light time for 1 au =
        499.00478 s.  A target AT the Sun's geocentric position (RA
        282.9 deg, Dec -22.87 deg on that date) puts the whole Sun-Earth
        distance along the line of sight, so the light-travel term must
        be -r * 499.00478 s = -490.67 s (negative: Earth trails the Sun
        along that line of sight by exactly r, so the Sun-centered clock
        records the photons EARLIER).  Tolerance 2.5 s covers the ~1.7
        s/degree sensitivity to the hand-read solar coordinates plus the
        observatory term (<0.03 s).
        """
        _, ltt_s = tm.hjd_utc_from_utc(
            2460312.5266, 282.9, -22.87, ephemeris="builtin")
        assert abs(float(ltt_s) - (-490.67)) < 2.5

    def test_maxim_jd_helio_header_cross_check(self):
        """Anchor 3: an independent software's own heliocentric stamp.

        Frame coe/cas/2023/152/cas15208.fts.fz (immutable archive), a
        1024 s High Gain StackPro exposure of M16, header cards recorded
        here verbatim on 2026-08-18:

            JD       = 2460096.8459347221   (UTC exposure start)
            EXPTIME  = 1024.0
            JD-HELIO = 2460096.8571224902
            CRVAL1/2 = 274.665588238, -13.8092470258

        MaxIm computed JD-HELIO at MID-exposure (start + 512 s).  Our own
        HJD_UTC at the same instant and position must agree to well under
        a second (MaxIm's solar-position approximation is the ~0.1 s
        class).  This one test cross-checks the mid-exposure policy AND
        the heliocentric geometry against software we did not write.
        """
        jd_mid, method = tm.jd_utc_mid(2460096.8459347221, 1024.0,
                                       "High Gain StackPro")
        assert method == tm.MID_STACKPRO
        hjd, _ = tm.hjd_utc_from_utc(jd_mid, 274.665588238, -13.8092470258,
                                     ephemeris="builtin")
        assert abs((float(hjd) - 2460096.8571224902) * 86400.0) < 0.5

    def test_bjd_exceeds_hjd_by_scale_offset(self):
        """BJD_TDB and HJD_UTC toward the same target differ by the scale
        offset (~69.2 s in 2023+) plus the helio-vs-bary geometry (< 4 s)
        — the sanity relation every downstream consumer relies on."""
        bjd, _, _ = tm.bjd_tdb_from_utc(2460100.5, 161.2, 33.35,
                                        ephemeris="builtin")
        hjd, _ = tm.hjd_utc_from_utc(2460100.5, 161.2, 33.35,
                                     ephemeris="builtin")
        diff_s = (float(bjd) - float(hjd)) * 86400.0
        assert 69.184 - 4.5 < diff_s < 69.184 + 4.5

    def test_vectorized_matches_scalar(self):
        jds = np.array([2460100.5, 2460200.75])
        ras = np.array([161.2, 274.7])
        decs = np.array([33.35, -13.81])
        vec, vec_ltt, _ = tm.bjd_tdb_from_utc(jds, ras, decs,
                                              ephemeris="builtin")
        for i in range(2):
            one, one_ltt, _ = tm.bjd_tdb_from_utc(jds[i], ras[i], decs[i],
                                                  ephemeris="builtin")
            assert abs(float(one) - vec[i]) < 1e-9      # < 0.1 ms
            assert abs(float(one_ltt) - vec_ltt[i]) < 1e-4

    def test_de440s_agrees_with_builtin_when_available(self):
        """The precision claim: builtin (analytic) vs DE440s (JPL) differ
        by the analytic ephemeris error — literature says km-class Earth
        position, i.e. < 20 ms of light time.  Skipped cleanly when the
        kernel is not cached on this host (offline CI)."""
        if tm.resolve_ephemeris() != "de440s":
            pytest.skip("de440s kernel not available on this host")
        b1, _, _ = tm.bjd_tdb_from_utc(2460440.63, 200.0, 10.0,
                                       ephemeris="de440s")
        b2, _, _ = tm.bjd_tdb_from_utc(2460440.63, 200.0, 10.0,
                                       ephemeris="builtin")
        assert abs(float(b1) - float(b2)) * 86400.0 < 0.020


# ---------------------------------------------------------------------------
# Eclipse helpers (pure)
# ---------------------------------------------------------------------------
class TestEclipseHelpers:
    def test_fold_phase_basic(self):
        ph = tm.fold_phase(np.array([100.0, 100.25, 100.5, 100.75]),
                           epoch=100.0, period=1.0)
        assert np.allclose(ph, [0.0, 0.25, -0.5, -0.25])

    def test_fold_phase_many_cycles(self):
        # 1330 cycles later, a frame exactly at conjunction folds to 0.
        p, e = 1.3590176, 2458211.194
        assert abs(tm.fold_phase(e + 1330 * p, e, p)) < 1e-9

    def test_fit_recovers_synthetic_offset(self):
        """A Gaussian dip injected at ph0 = +0.013 with noise must come
        back within its stated error (deterministic seed)."""
        rng = np.random.default_rng(42)
        ph = np.concatenate([rng.uniform(-0.1, 0.1, 60),
                             rng.uniform(-0.02, 0.04, 30)])
        true = 0.55 * np.exp(-((ph - 0.013) ** 2) / (2 * 0.018 ** 2))
        dm = true + rng.normal(0, 0.03, ph.size)
        fit = tm.fit_eclipse_offset(ph, dm, errs=np.full(ph.size, 0.03))
        assert fit["ph0"] is not None
        assert abs(fit["ph0"] - 0.013) < max(3 * fit["ph0_err"], 0.004)
        assert 0.3 < fit["depth"] < 0.8

    def test_fit_rejects_pure_noise_or_brightening(self):
        rng = np.random.default_rng(7)
        ph = rng.uniform(-0.1, 0.1, 80)
        # A BRIGHTENING (negative dmag dip) is not an eclipse: the fitter
        # must not report a positive-depth solution deeper than noise.
        dm = -0.5 * np.exp(-(ph ** 2) / (2 * 0.02 ** 2)) \
            + rng.normal(0, 0.02, ph.size)
        fit = tm.fit_eclipse_offset(ph, dm)
        assert fit["ph0"] is None or fit["depth"] < 0.1
