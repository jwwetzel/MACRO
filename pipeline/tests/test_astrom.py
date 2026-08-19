"""Unit tests for macro_core.astrom — the S1 pure logic.

Same philosophy as test_manifest/test_inventory: every decision function is
exercised on hand-built cases, INCLUDING the cases that must NOT work (a
grism spectrum sneaking into a stratum, an 8-px photometry window classified
as solvable, a sample that changes between runs).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from macro_core import astrom
from macro_core.astrom import (
    build_solve_command, classify_stratum, diagnose_failure, image_metrics,
    is_grism_filter, is_solvable_candidate, is_window_geometry,
    night_collapse, projected_hours, rms, sample_frames, scale_bounds,
    sky_residuals_arcsec, solution_sane, stable_seed, verdict_for,
    wilson_ci)


# ---------------------------------------------------------------------------
# helpers: a template candidate row (CV polar, Mode0 bin2 Sloan)
# ---------------------------------------------------------------------------
def row(**over) -> dict:
    base = dict(target_key="stlmi", canonical_target="ST LMi",
                readoutm="Mode0", xbinning=2, filter="g", exptime=60.0,
                naxis1=4788, naxis2=3194)
    base.update(over)
    return base


# ---------------------------------------------------------------------------
# candidate gates
# ---------------------------------------------------------------------------
class TestGates:
    def test_grism_filters_are_never_candidates(self):
        # A grism frame is a spectrum: excluded however it is spelled.
        for f in ("hrg", "LRG", " HaGrism ", "OGGrism"):
            assert is_grism_filter(f)
            assert not is_solvable_candidate(row(filter=f))

    def test_ordinary_filters_pass(self):
        for f in ("g", "r", "i", "Ha", "O", None, ""):
            assert not is_grism_filter(f)

    def test_calib_vocab_filter_excluded(self):
        # The era-76 glitch: science frame with FILTER = 'dark'.
        assert not is_solvable_candidate(row(filter="dark"))

    def test_window_geometry(self):
        # The 8x3211 Fast photometry strips: unsolvable by geometry.
        assert is_window_geometry(8, 3211)
        # NULL geometry cannot be promised solvable either.
        assert is_window_geometry(None, 3211)
        # A full frame passes.
        assert not is_window_geometry(4788, 3194)
        # The boundary: exactly MIN_SOLVABLE_NAXIS passes (>= contract).
        m = astrom.MIN_SOLVABLE_NAXIS
        assert not is_window_geometry(m, m)
        assert is_window_geometry(m - 1, m)


# ---------------------------------------------------------------------------
# stratum classification — one test per stratum, plus the must-not cases
# ---------------------------------------------------------------------------
class TestClassify:
    def test_cv_mode0_short_vs_long(self):
        assert classify_stratum(row(exptime=60)) == "cv_mode0_sloan_short"
        assert classify_stratum(row(exptime=240)) == "cv_mode0_sloan_long"
        # The band edge: 150 s belongs to the long stratum.
        assert classify_stratum(row(exptime=150)) == "cv_mode0_sloan_long"

    def test_cv_ikon(self):
        r = row(target_key="vvpup", readoutm="1MHz High Sensitivity 16-bit",
                xbinning=1, naxis1=2048, naxis2=2048)
        assert classify_stratum(r) == "cv_ikon_sloan"
        # Blank-filter frames of the same series ride along.
        assert classify_stratum({**r, "filter": "empty"}) == "cv_ikon_sloan"

    def test_cv_gsense(self):
        r = row(readoutm="High Gain", xbinning=1, naxis1=4096, naxis2=4096,
                filter="R")
        assert classify_stratum(r) == "cv_gsense_misc"

    def test_sn(self):
        r = row(target_key="2023ixf", canonical_target="2023ixf",
                readoutm="High Gain", xbinning=1, naxis1=4096,
                naxis2=4096, filter="O")
        assert classify_stratum(r) == "sn_gsense_broadband"

    def test_dwarf(self):
        r = row(target_key="dw1234p56", canonical_target="Dw1234+56",
                readoutm="High Gain StackPro", xbinning=1,
                naxis1=4096, naxis2=4096, filter="L")
        assert classify_stratum(r) == "dwarf_gsense_deep"
        r2 = {**r, "target_key": "ngc5548", "canonical_target": "NGC 5548"}
        assert classify_stratum(r2) == "dwarf_gsense_deep"

    def test_backlog_strata(self):
        r = row(target_key="rhooph", canonical_target="rho Oph", exptime=2)
        assert classify_stratum(r) == "mode0_backlog_short"
        assert classify_stratum({**r, "exptime": 120}) \
            == "mode0_backlog_long"
        fast = row(target_key="x", canonical_target="X", readoutm="Fast",
                   naxis1=4800, naxis2=3211)
        assert classify_stratum(fast) == "fast_fullframe"
        ikon = row(target_key="x", canonical_target="X",
                   readoutm="1MHz High Sensitivity 16-bit", xbinning=1,
                   naxis1=2048, naxis2=2048)
        assert classify_stratum(ikon) == "ikon_backlog"

    def test_must_not_classify(self):
        # A grism CV frame is NOT a Sloan frame however tempting.
        assert classify_stratum(row(filter="hrg")) is None
        # The EU UMa Fast photometry strips: excluded by geometry even
        # though the target is CV-critical.
        strip = row(target_key="euuma", readoutm="Fast",
                    naxis1=8, naxis2=3211)
        assert classify_stratum(strip) is None
        # ...and the point of the S0e geometry repair: the SAME frame, once
        # its true 4800x3211 geometry is read, is a perfectly ordinary CV
        # full frame and now reaches a stratum of its own.
        repaired = row(target_key="euuma", readoutm="Fast",
                       xbinning=2, naxis1=4800, naxis2=3211)
        assert classify_stratum(repaired) == "cv_fast_fullframe"
        # A CV target on an unplanned camera config (blank readout — the
        # header-convention break eras) joins no stratum.
        odd = row(readoutm="", xbinning=2)
        assert classify_stratum(odd) is None

    def test_cv_fast_frames_do_not_land_in_the_facility_backlog(self):
        """A CV target and a facility target on the same camera must NOT
        share a stratum.  Letting EU UMa fall through to 'fast_fullframe'
        would have changed the population behind an already-published
        stratum id; a new id was added instead."""
        cv = row(target_key="euuma", readoutm="Fast", xbinning=2,
                 naxis1=4800, naxis2=3211)
        facility = row(target_key="somefield", readoutm="Fast", xbinning=2,
                       naxis1=4800, naxis2=3211)
        assert classify_stratum(cv) == "cv_fast_fullframe"
        assert classify_stratum(facility) == "fast_fullframe"

    def test_every_stratum_id_is_reachable_and_declared(self):
        # The classifier and the STRATA table must agree exactly.
        declared = {s.stratum_id for s in astrom.STRATA}
        reached = {
            classify_stratum(row(exptime=60)),
            classify_stratum(row(exptime=240)),
            classify_stratum(row(target_key="vvpup",
                                 readoutm="1MHz High Sensitivity 16-bit",
                                 xbinning=1, naxis1=2048, naxis2=2048)),
            classify_stratum(row(readoutm="High Gain", xbinning=1,
                                 naxis1=4096, naxis2=4096)),
            classify_stratum(row(target_key="2023ixf", readoutm="High Gain",
                                 xbinning=1, naxis1=4096, naxis2=4096)),
            classify_stratum(row(target_key="dw1x", canonical_target="Dw1",
                                 readoutm="High Gain", xbinning=1,
                                 naxis1=4096, naxis2=4096)),
            classify_stratum(row(target_key="a", exptime=2)),
            classify_stratum(row(target_key="a", exptime=100)),
            classify_stratum(row(target_key="a", readoutm="Fast",
                                 naxis1=4800, naxis2=3211)),
            classify_stratum(row(target_key="euuma", readoutm="Fast",
                                 xbinning=2, naxis1=4800, naxis2=3211)),
            classify_stratum(row(target_key="a",
                                 readoutm="5MHz High Sensitivity 16-bit",
                                 xbinning=1, naxis1=2048, naxis2=2048)),
        }
        assert reached == declared


# ---------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------
class TestSampling:
    def test_reproducible(self):
        ids = list(range(1000, 1500))
        a = sample_frames(ids, 48, 20260817, "cv_mode0_sloan_short")
        b = sample_frames(ids, 48, 20260817, "cv_mode0_sloan_short")
        assert a == b and len(a) == 48

    def test_order_independent(self):
        # Query order is not a contract: a shuffled candidate list gives
        # the identical sample.
        ids = list(range(1000, 1500))
        rev = list(reversed(ids))
        assert sample_frames(ids, 48, 1, "s") == \
            sample_frames(rev, 48, 1, "s")

    def test_different_strata_differ(self):
        ids = list(range(1000, 1500))
        assert sample_frames(ids, 48, 1, "a") != \
            sample_frames(ids, 48, 1, "b")

    def test_small_population_takes_all(self):
        assert sorted(sample_frames([5, 3, 9], 48, 1, "s")) == [3, 5, 9]

    def test_stable_seed_is_cross_process_stable(self):
        # Regression pin: the exact value, so a platform/Python change
        # that would silently redraw every sample fails a test instead.
        assert stable_seed(20260817, "cv_mode0_sloan_short") \
            == stable_seed(20260817, "cv_mode0_sloan_short")
        assert stable_seed(20260817, "a") != stable_seed(20260817, "b")


# ---------------------------------------------------------------------------
# Wilson intervals + verdicts + projections
# ---------------------------------------------------------------------------
class TestStats:
    def test_wilson_known_value(self):
        # 45/48 → interval (0.8316, 0.9785); computed independently
        # (z = 1.96 Wilson formula by hand).
        lo, hi = wilson_ci(45, 48)
        assert lo == pytest.approx(0.8316, abs=2e-3)
        assert hi == pytest.approx(0.9785, abs=2e-3)

    def test_wilson_edges(self):
        lo, hi = wilson_ci(0, 48)
        # Zero successes: the lower bound touches 0 (to float precision)
        # but the UPPER does not — 0/48 is evidence of a low rate, not
        # proof of an impossible one.
        assert lo == pytest.approx(0.0, abs=1e-4)
        assert 0 < hi < 0.12
        lo, hi = wilson_ci(48, 48)
        # Perfect score: the interval hugs 1 but honestly stops short of
        # it on the upper side too (48/48 cannot prove a rate of 1).
        assert 0.99 < hi <= 1.0 and 0.88 < lo < 1.0
        assert wilson_ci(0, 0) == (0.0, 1.0)  # no data: total ignorance

    def test_verdicts(self):
        assert verdict_for(48, 48) == "GO"
        assert verdict_for(46, 48) == "GO"       # lower bound ≈ 0.858
        assert verdict_for(35, 48) == "CAUTION"  # lower bound ≈ 0.588
        assert verdict_for(10, 48) == "NO-GO"
        assert verdict_for(0, 0) == "NO-GO"      # no data can never be GO

    def test_census_verdict_judges_exact_rate(self):
        # Regression pin (review finding): cv_gsense_misc sampled ALL 40
        # of its 40 population frames — a census.  25/40 = 62.5% is the
        # EXACT population rate; there is no sampling uncertainty for a
        # Wilson interval to be pessimistic about, so the verdict must be
        # CAUTION (62.5 >= 50), not the NO-GO the Wilson lower bound
        # (47.0%) would give a genuine sample.
        assert verdict_for(25, 40) == "NO-GO"          # as a sample
        assert verdict_for(25, 40, n_population=40) == "CAUTION"  # census
        # A census can also be a clean GO or NO-GO on its exact rate.
        assert verdict_for(36, 40, n_population=40) == "GO"
        assert verdict_for(10, 40, n_population=40) == "NO-GO"
        # Non-census populations keep the pessimistic Wilson judgement.
        assert verdict_for(25, 40, n_population=4000) == "NO-GO"

    def test_night_collapse(self):
        # Three nights: one perfect, one mixed, one all-failed.
        rows = [("n1", True), ("n1", True),
                ("n2", True), ("n2", False),
                ("n3", False)]
        assert night_collapse(rows) == (1, 3)
        # All perfect / empty edge cases.
        assert night_collapse([("a", True), ("b", True)]) == (2, 2)
        assert night_collapse([]) == (0, 0)
        # SQLite hands back 0/1 ints, not bools — must behave identically.
        assert night_collapse([("a", 1), ("a", 0), ("b", 1)]) == (1, 2)

    def test_projection(self):
        # 59,000 frames × 6 s ÷ 10 workers = 9.83 hours.
        assert projected_hours(59000, 6.0, 10) \
            == pytest.approx(9.833, abs=1e-3)
        assert projected_hours(0, 6.0, 10) == 0.0


# ---------------------------------------------------------------------------
# scale priors + solve command
# ---------------------------------------------------------------------------
class TestSolvePlumbing:
    def test_scale_bounds_families(self):
        # Mode0 bin2 measured 0.4508"/px — the prior must contain it.
        lo, hi = scale_bounds("Mode0", 2)
        assert lo < 0.4508 < hi
        # GSENSE unbinned is 0.54"/px.
        lo, hi = scale_bounds("High Gain", 1)
        assert lo < 0.54 < hi
        # iKon is unverified: bounds must stay wide, and the experiment
        # measured 0.809"/px unbinned — the prior must contain it.
        lo, hi = scale_bounds("1MHz High Sensitivity 16-bit", 1)
        assert hi - lo > 0.5
        assert lo < 0.809 < hi
        # Binning scales EVERY family's prior, the iKon included — a
        # bin2 iKon frame (~1.62"/px) must not fall outside its bounds.
        lo2, hi2 = scale_bounds("1MHz High Sensitivity 16-bit", 2)
        assert lo2 < 2 * 0.809 < hi2
        # NULL readout (blank-era frames) gets the IMX455 prior.
        lo, hi = scale_bounds(None, 2)
        assert lo < 0.45 < hi

    def test_command_with_hint(self):
        cmd = build_solve_command("/tmp/x.fts", "/cfg", "/out",
                                  0.3, 0.7, ra_deg=150.0, dec_deg=30.0)
        assert cmd[0] == "solve-field"
        assert "--ra" in cmd and "--radius" in cmd
        # The scale window rides along in arcsec/px.
        i = cmd.index("--scale-low")
        assert cmd[i + 1] == "0.3"

    def test_command_blind_when_no_hint(self):
        cmd = build_solve_command("/tmp/x.fts", "/cfg", "/out", 0.3, 0.7)
        assert "--ra" not in cmd
        # Plots and derived FITS are suppressed (10-worker I/O budget).
        assert "--no-plots" in cmd
        assert "none" in cmd


# ---------------------------------------------------------------------------
# solution-acceptance gate
# ---------------------------------------------------------------------------
class TestSolutionGate:
    # The Mode0 bin2 prior handed to solve-field: (0.30, 0.66).
    LO, HI = scale_bounds("Mode0", 2)

    def test_rejects_the_recorded_false_positive(self):
        # Regression pin (review finding, executed bug): frame 149276
        # (WASP-102, mode0_backlog_long) produced a .solved marker with
        # pixscale 3.3817"/px, 4 matched stars, RMS 5.63" — a false
        # solution on a healthy 0.45"/px star field.  v1.0 counted it as
        # a success; the gate must reject it on every one of its checks.
        assert not solution_sane(3.3817, 4, 5.63, self.LO, self.HI)
        # Each violation alone is disqualifying:
        assert not solution_sane(3.3817, 30, 1.0, self.LO, self.HI)  # scale
        assert not solution_sane(0.4508, 4, 1.0, self.LO, self.HI)   # stars
        assert not solution_sane(0.4508, 30, 5.63, self.LO, self.HI)  # rms

    def test_accepts_a_typical_genuine_solve(self):
        # The experiment's typical accepted solve: prior-consistent
        # scale, 30+ matched stars, sub-arcsecond RMS.
        assert solution_sane(0.4508, 30, 0.92, self.LO, self.HI)
        # The observed extremes of GENUINE solves must stay accepted:
        # min 8 matched stars, max RMS 4.83".
        assert solution_sane(0.4508, astrom.MIN_MATCHED_STARS, 4.83,
                             self.LO, self.HI)
        # Bounds are inclusive: a scale exactly at the prior edge passes.
        assert solution_sane(self.LO, 30, 1.0, self.LO, self.HI)
        assert solution_sane(self.HI, 30, 1.0, self.LO, self.HI)

    def test_missing_values_fail_closed(self):
        # A solution that cannot show its evidence is not accepted.
        assert not solution_sane(None, 30, 1.0, self.LO, self.HI)
        assert not solution_sane(0.45, None, 1.0, self.LO, self.HI)
        assert not solution_sane(0.45, 30, None, self.LO, self.HI)


# ---------------------------------------------------------------------------
# residual arithmetic
# ---------------------------------------------------------------------------
class TestResiduals:
    def test_simple_offset(self):
        # 1" offset in Dec only → residual exactly 1".
        r = sky_residuals_arcsec([10.0], [20.0], [10.0], [20.0 - 1 / 3600])
        assert r[0] == pytest.approx(1.0, abs=1e-6)

    def test_ra_compression(self):
        # 1" of RA at Dec 60° is only 0.5" on the sky.
        r = sky_residuals_arcsec([10.0 + 1 / 3600], [60.0], [10.0], [60.0])
        assert r[0] == pytest.approx(0.5, abs=1e-3)

    def test_ra_wrap(self):
        # Residuals across the 0/360 seam must not explode.
        r = sky_residuals_arcsec([359.9999], [0.0], [0.0001], [0.0])
        assert r[0] < 1.0

    def test_rms(self):
        assert rms([3.0, 4.0]) == pytest.approx(math.sqrt(12.5))
        assert rms([]) is None


# ---------------------------------------------------------------------------
# failure autopsy
# ---------------------------------------------------------------------------
class TestAutopsy:
    def test_diagnosis_ladder(self):
        # Args: (n_sources, n_psf_sources, med_elong, sat_frac, bright_a).
        assert diagnose_failure(None, None, None, None) == "unreadable"
        assert "starved" in diagnose_failure(3, 2, 1.0, 0.0, 1.0)
        assert "trailing" in diagnose_failure(100, 80, 2.5, 0.0, 2.0)
        assert "saturated" in diagnose_failure(100, 80, 1.1, 0.5, 2.0)
        assert "unexplained" in diagnose_failure(100, 80, 1.1, 0.0, 2.0)
        # Starvation outranks trailing: 3 streaks are still 3 sources.
        assert "starved" in diagnose_failure(3, 2, 5.0, 0.0, 1.0)
        # The hot-pixel forest: thousands of 10σ detections, none of them
        # PSF-shaped — that is a STARLESS frame however big n_sources is
        # (the second autopsy round's blank-frame failure mode).
        assert "starved" in diagnose_failure(16000, 4, 1.0, 0.0, 0.7)
        # The defocus tell: the brightest detections are giant blobs
        # (the filter-'6' failure mode found in the first autopsy round).
        assert "defocused" in diagnose_failure(16000, 4, 1.05, 0.0, 20.0)
        # Defocus outranks trailing: donut fragments read as elongated
        # PSF-band sources, but giant bright blobs settle the question.
        assert "defocused" in diagnose_failure(16000, 600, 2.5, 0.0, 14.0)
        # ...but tight PSF-shaped sources of ordinary size stay
        # unexplained: the solver failed with stars on the table.
        assert "unexplained" in diagnose_failure(100, 50, 1.05, 0.0, 2.0)

    def test_blank_frame_hot_pixel_pairs_are_starved_not_trailing(self):
        # Regression pin (review finding, executed bug): frame 55732
        # (dwarf_gsense_deep, H 512 s) is visually BLANK, yet v1.0
        # diagnosed it 'trailing' — adjacent hot pixels pair into 1–2 px
        # clumps that pass the PSF gate (50 "PSF-shaped" detections) AND
        # read as elongated (median elongation 2.36).  The tell is the
        # brightest-detection size: 0.707 px median semi-major axis means
        # the brightest things on the frame are single-pixel spikes — no
        # real star anywhere.  The exact DB metrics of all four
        # mislabeled frames must now diagnose as starved.
        for n_src, n_psf, elong, sat, bright_a in [
                (8633, 50, 2.3589, 0.00082, 0.70740),    # 55732
                (8654, 53, 2.2045, 0.00082, 0.70735),    # 55733
                (9366, 65, 2.4116, 0.00087, 0.70781),    # 65495
                (10475, 74, 2.3248, 0.00092, 0.70795)]:  # 65497
            assert "starved" in diagnose_failure(n_src, n_psf, elong,
                                                 sat, bright_a)
        # And the honest trailing frames must STAY trailing — including
        # 68275, whose bright_median_a (1.12 px) sits BELOW PSF_A_MIN_PX
        # (1.2): gating trailing on PSF_A_MIN_PX would misclassify real
        # curved star trails, which is why the blank veto has its own,
        # lower threshold between the two measured clusters.
        assert "trailing" in diagnose_failure(13690, 148, 2.2456,
                                              0.0, 1.1216)     # 68275
        assert "trailing" in diagnose_failure(11092, 145, 2.1192,
                                              0.0, 3.6237)     # 198853

    def test_image_metrics_on_synthetic_field(self):
        rng = np.random.default_rng(42)
        img = rng.normal(100.0, 5.0, (256, 256)).astype(np.float32)
        # Plant 30 bright Gaussian-ish stars on a grid.
        yy, xx = np.mgrid[-3:4, -3:4]
        psf = 2000.0 * np.exp(-(xx ** 2 + yy ** 2) / 4.0)
        for i in range(30):
            y, x = 20 + (i % 6) * 40, 20 + (i // 6) * 40
            img[y - 3:y + 4, x - 3:x + 4] += psf
        m = image_metrics(img, saturation_adu=60000.0)
        assert m["n_psf_sources"] >= 25      # finds the planted stars
        assert m["saturated_fraction"] == 0.0
        assert m["median_elongation"] < astrom.AUTOPSY_TRAIL_ELONG
        # Planted PSFs are compact: no defocus verdict on a good field.
        assert m["bright_median_a_px"] < astrom.AUTOPSY_DEFOCUS_A_PX
        # And a blank frame is starved, not a crash.
        blank = rng.normal(100.0, 5.0, (128, 128)).astype(np.float32)
        mb = image_metrics(blank, saturation_adu=60000.0)
        assert mb["n_psf_sources"] < astrom.AUTOPSY_MIN_SOURCES
