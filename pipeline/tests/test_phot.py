"""Unit tests for macro_phot — every pure S4 function.

The centerpiece is the REQUIRED ensemble-recovery test: the Honeycutt
solver must recover a synthetic injected zero-point pattern within
tolerance, drop an injected variable star from the comparison set, and
leave a constant star's statistics consistent with its noise.  The
must-NOT cases (a flipped field failing to match under the wrong parity,
an aperture blowing past its guard rails, matching across the tolerance)
matter as much as the must cases.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Make the package importable regardless of pytest's working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import photometry as ph
from macro_phot import ensemble as ens
from macro_phot import errors as err
from macro_phot import gaia as gg


# ---------------------------------------------------------------------------
# Plate scale and apertures
# ---------------------------------------------------------------------------
class TestPlateScale:
    def test_cmos_era_matches_wcs(self):
        # Era 76 headers: XPIXSZ 7.52 um (binned), FOCALLEN 3454 mm; the
        # plate-solved frames' CD matrix says 0.451 "/px.
        s = ph.plate_scale_arcsec_per_px(7.52, 3454.0)
        assert abs(s - 0.449) < 0.002

    def test_ikon_era(self):
        # Era 72: Andor iKon 13.5 um pixels, same telescope.
        s = ph.plate_scale_arcsec_per_px(13.5, 3454.0)
        assert abs(s - 0.806) < 0.002

    def test_missing_cards_give_none(self):
        assert ph.plate_scale_arcsec_per_px(None, 3454.0) is None
        assert ph.plate_scale_arcsec_per_px(7.52, None) is None
        assert ph.plate_scale_arcsec_per_px(0.0, 3454.0) is None

    def test_aperture_scales_inversely_with_plate_scale(self):
        # 4" at 0.449 "/px = 8.9 px; at 0.806 "/px = 5.0 px — the SAME sky
        # aperture across eras is the whole point.
        a_cmos = ph.aperture_radius_px(0.449)
        a_ikon = ph.aperture_radius_px(0.806)
        assert abs(a_cmos - 4.0 / 0.449) < 1e-9
        assert abs(a_ikon - 4.0 / 0.806) < 1e-9
        assert a_cmos > a_ikon

    def test_aperture_guard_rails(self):
        # Absurd plate scales clamp instead of silently producing 1-px or
        # 60-px apertures.
        assert ph.aperture_radius_px(10.0) == ph.APERTURE_MIN_PX
        assert ph.aperture_radius_px(0.01) == ph.APERTURE_MAX_PX
        assert ph.aperture_radius_px(None) is None

    def test_fwhm_from_ab_gaussian(self):
        # A circular Gaussian with sigma = 2 px has FWHM 4.71 px.
        f = ph.fwhm_from_ab(np.array([2.0]), np.array([2.0]))
        assert abs(f[0] - 2.3548 * 2.0) < 1e-3


# ---------------------------------------------------------------------------
# Reference-frame choice
# ---------------------------------------------------------------------------
class TestChooseReference:
    def test_most_detections_among_sharp(self):
        stats = [(1, 100, 5.0, 0), (2, 150, 5.5, 0), (3, 200, 9.0, 0),
                 (4, 120, 4.0, 0)]
        # FWHM 9.0 is beyond the 75th percentile — frame 3 is barred even
        # though it has the most detections.
        assert ph.choose_reference(stats) == 2

    def test_solved_breaks_ties(self):
        stats = [(1, 100, 5.0, 0), (2, 100, 5.0, 1)]
        assert ph.choose_reference(stats) == 2

    def test_empty_series(self):
        assert ph.choose_reference([]) is None
        assert ph.choose_reference([(1, 0, None, 0)]) is None

    def test_ranking_orders_all_candidates(self):
        # The ranking exists so a doubled top candidate can be REJECTED
        # and the next-best adopted; it must agree with choose_reference
        # at rank 0 and order the rest by the same rule.
        stats = [(1, 100, 5.0, 0), (2, 150, 5.5, 0), (3, 200, 9.0, 0),
                 (4, 120, 4.0, 0)]
        ranking = ph.rank_references(stats)
        assert ranking[0] == ph.choose_reference(stats)
        assert ranking == [2, 4, 1]       # frame 3 barred by FWHM cut
        assert ph.rank_references([]) == []


class TestReferenceDoubleQC:
    """The double-image detector (executed bug: both VV Pup references of
    the first build were double-imaged exposures — 83-89% of their stars
    were equal-brightness pairs at one constant offset — and the
    most-detections selection rule actively rewarded them)."""

    def _clean_field(self, n=400, seed=11):
        rng = np.random.default_rng(seed)
        xy = rng.uniform(0, 4000, size=(n, 2))
        flux = 10 ** rng.uniform(3, 6, size=n)
        return xy, flux

    def test_clean_field_is_not_flagged(self):
        xy, flux = self._clean_field()
        frac = ph.paired_fraction(xy, flux, radius_px=12.0)
        assert frac < ph.REF_DOUBLED_MAX_FRAC

    def test_doubled_field_is_flagged(self):
        # Double image: every star + a companion at ONE constant offset
        # with near-equal flux (the guiding-jump signature; the real bad
        # references paired at (-0.5, -9.9) px with flux ratios ~0.9).
        xy, flux = self._clean_field()
        xy2 = np.vstack([xy, xy + np.array([-0.5, -9.9])])
        flux2 = np.concatenate([flux, 0.9 * flux])
        frac = ph.paired_fraction(xy2, flux2, radius_px=12.0)
        assert frac > 0.8
        assert frac > ph.REF_DOUBLED_MAX_FRAC

    def test_faint_companions_do_not_count(self):
        # A near neighbour 10x fainter is a chance blend, not a double
        # image — the flux-ratio arm must refuse it.
        xy, flux = self._clean_field(n=200)
        xy2 = np.vstack([xy, xy + np.array([0.0, 8.0])])
        flux2 = np.concatenate([flux, flux / 10.0])
        frac = ph.paired_fraction(xy2, flux2, radius_px=12.0)
        assert frac < ph.REF_DOUBLED_MAX_FRAC

    def test_tiny_catalogs(self):
        assert ph.paired_fraction(np.zeros((0, 2)), np.zeros(0), 10.0) == 0.0
        assert ph.paired_fraction(np.array([[1.0, 1.0]]),
                                  np.array([5.0]), 10.0) == 0.0


# ---------------------------------------------------------------------------
# One-to-one matching
# ---------------------------------------------------------------------------
class TestMatching:
    def test_exact_match(self):
        ref = np.array([[0.0, 0.0], [10.0, 0.0], [0.0, 10.0]])
        xy = ref + 0.1
        idx = ph.match_one_to_one(ref, xy, tol_px=1.0)
        assert list(idx) == [0, 1, 2]

    def test_tolerance_is_respected(self):
        ref = np.array([[0.0, 0.0]])
        xy = np.array([[2.5, 0.0]])
        assert ph.match_one_to_one(ref, xy, tol_px=2.0)[0] == -1
        assert ph.match_one_to_one(ref, xy, tol_px=3.0)[0] == 0

    def test_one_to_one_greedy_by_distance(self):
        # Two detections compete for one reference star: the closer wins,
        # the farther goes unmatched (must NOT double-assign).
        ref = np.array([[0.0, 0.0]])
        xy = np.array([[0.5, 0.0], [0.2, 0.0]])
        idx = ph.match_one_to_one(ref, xy, tol_px=1.0)
        assert idx[1] == 0 and idx[0] == -1

    def test_empty_inputs(self):
        assert len(ph.match_one_to_one(np.zeros((0, 2)),
                                       np.zeros((0, 2)), 1.0)) == 0

    def test_tolerance_scaling(self):
        assert ph.match_tolerance_px(None) == ph.MATCH_TOL_MIN_PX
        assert ph.match_tolerance_px(10.0) == ph.MATCH_TOL_FWHM * 10.0


# ---------------------------------------------------------------------------
# Magnitudes
# ---------------------------------------------------------------------------
class TestMagnitudes:
    def test_exposure_normalization(self):
        # Twice the flux in twice the exposure = the same magnitude.
        m1 = ph.instrumental_mag(np.array([1000.0]), 10.0)
        m2 = ph.instrumental_mag(np.array([2000.0]), 20.0)
        assert abs(m1[0] - m2[0]) < 1e-12

    def test_nonpositive_flux_is_nan(self):
        m = ph.instrumental_mag(np.array([0.0, -5.0, 100.0]), 10.0)
        assert np.isnan(m[0]) and np.isnan(m[1]) and np.isfinite(m[2])

    def test_mag_error_snr(self):
        # SNR 100 -> ~10.86 mmag.
        e = ph.mag_error(np.array([100.0]), np.array([1.0]))
        assert abs(e[0] - 0.010857) < 1e-5
        assert np.isnan(ph.mag_error(np.array([0.0]), np.array([1.0]))[0])


# ---------------------------------------------------------------------------
# The Honeycutt ensemble — REQUIRED synthetic-recovery tests
# ---------------------------------------------------------------------------
def _synthetic(seed=42, S=30, F=80, noise=0.01, missing_frac=0.15):
    """Build a synthetic series: constant stars + a known ZP pattern.

    The injected zero-point pattern mixes a smooth night-long drift with a
    discrete 0.3-mag cloud event — both must come back out.
    """
    rng = np.random.default_rng(seed)
    M_true = rng.uniform(12.0, 17.0, size=S)
    t = np.linspace(0, 1, F)
    ZP_true = 0.15 * np.sin(2 * np.pi * t) + np.where(
        (t > 0.4) & (t < 0.5), 0.3, 0.0)     # drift + cloud step
    ZP_true -= ZP_true.mean()                # match the solver's gauge
    sig = np.full((S, F), noise)
    mag = (M_true[:, None] + ZP_true[None, :]
           + rng.normal(0, noise, size=(S, F)))
    # Miss some observations at random (an inhomogeneous ensemble by
    # construction — Honeycutt's whole reason for existing).
    mask = rng.random((S, F)) < missing_frac
    mag[mask] = np.nan
    return M_true, ZP_true, mag, sig


class TestEnsembleSolver:
    def test_recovers_injected_zp_pattern(self):
        # THE required test: the solver must recover the injected pattern
        # within tolerance.  With ~25 stars voting at 10 mmag noise, each
        # frame's ZP is determined to sigma ~ 2 mmag.  Assert on the RMS
        # of the recovery error (expectation ~2 mmag, tolerance 6 = a real
        # 3x margin) — asserting on the MAX over 80 frames, as an earlier
        # version did, sets the statistic's own expectation at ~2.9 sigma
        # ~ 5.8 mmag, so a 6 mmag ceiling had ~1.0x margin and neighboring
        # seeds crossed it (seed 3: 6.10 mmag).  A loose 10 mmag cap on
        # the max (~3.5 sigma per frame) still catches any gross failure.
        M_true, ZP_true, mag, sig = _synthetic()
        sol = ens.solve_ensemble(mag, sig)
        assert sol.converged
        err_zp = sol.zp - ZP_true
        assert np.sqrt(np.mean(err_zp ** 2)) < 0.006
        assert np.max(np.abs(err_zp)) < 0.010
        # Star magnitudes come back too (same gauge).
        assert np.nanmax(np.abs(sol.mean_mag - M_true)) < 0.01

    def test_frame_with_no_voting_stars_gets_nan_zp(self):
        # A frame whose comparison stars are ALL unmeasured must come back
        # with ZP = NaN — never a fabricated 0.0 that downstream code
        # would mistake for a valid zero point (executed bug: 34 real
        # frames carried zp = 0.0 with n_star_used = 0, and their target
        # measurements entered light curves uncorrected).
        M_true, ZP_true, mag, sig = _synthetic()
        mag[:, 17] = np.nan               # frame 17: nobody measured
        sol = ens.solve_ensemble(mag, sig)
        assert np.isnan(sol.zp[17])
        assert np.isnan(sol.zp_err[17])
        assert sol.n_star_used[17] == 0
        # The gauge still holds exactly over the frames that HAVE a ZP.
        assert abs(np.nanmean(sol.zp)) < 1e-9
        # And the other frames' recovery is unharmed.
        ok = np.isfinite(sol.zp)
        resid = sol.zp[ok] - ZP_true[ok]
        assert np.sqrt(np.mean((resid - resid.mean()) ** 2)) < 0.006

    def test_masked_out_frame_gets_nan_zp(self):
        # Same guarantee through the comp-selection path: when the MASK
        # removes every voting star of one frame (all comps clipped, say),
        # that frame must not receive an invented zero point.
        _, _, mag, sig = _synthetic()
        mask = np.isfinite(mag)
        mask[:, 3] = False
        sol = ens.solve_ensemble(mag, sig, mask=mask)
        assert np.isnan(sol.zp[3])
        assert sol.n_star_used[3] == 0

    def test_gauge_mean_zp_zero(self):
        _, _, mag, sig = _synthetic()
        sol = ens.solve_ensemble(mag, sig)
        assert abs(np.mean(sol.zp)) < 1e-9

    def test_outlier_is_clipped_not_absorbed(self):
        # One cosmic-ray hit must be masked, not smeared into the ZP.
        M_true, ZP_true, mag, sig = _synthetic()
        mag[3, 10] += 1.0                 # 100-sigma outlier
        sol = ens.solve_ensemble(mag, sig)
        assert sol.clipped[3, 10]
        assert abs(sol.zp[10] - ZP_true[10]) < 0.01

    def test_empty_and_single_star(self):
        # Degenerate inputs must not hang or crash.
        sol = ens.solve_ensemble(np.full((2, 3), np.nan),
                                 np.full((2, 3), np.nan))
        assert np.all(np.isnan(sol.mean_mag))
        sol1 = ens.solve_ensemble(np.array([[15.0, 15.1]]),
                                  np.array([[0.01, 0.01]]))
        assert np.isfinite(sol1.mean_mag[0])


class TestCompSelection:
    def test_variable_star_is_dropped(self):
        # Inject a 0.4-mag sinusoidal variable: it must land in
        # 'dropped_unstable', never in the comp set.
        M_true, ZP_true, mag, sig = _synthetic()
        F = mag.shape[1]
        mag[7] += 0.4 * np.sin(np.linspace(0, 6 * np.pi, F))
        sel = ens.select_comps(mag, sig, target_row=0)
        assert sel.role[7] == "dropped_unstable"
        assert sel.role[0] == "target"
        assert (sel.role == "comp").sum() >= ens.MIN_COMPS
        # The ZP still comes back clean despite the variable.
        assert np.max(np.abs(sel.solution.zp - ZP_true)) < 0.01

    def test_target_never_comps_even_if_constant(self):
        _, _, mag, sig = _synthetic()
        sel = ens.select_comps(mag, sig, target_row=5)
        assert sel.role[5] == "target"

    def test_check_stars_held_out(self):
        _, _, mag, sig = _synthetic()
        sel = ens.select_comps(mag, sig, target_row=0)
        assert (sel.role == "check").sum() == ens.N_CHECK_STARS
        # Checks must not overlap comps by construction.
        assert not np.any((sel.role == "check") & (sel.role == "comp"))

    def test_ephemeral_star_stays_field(self):
        # A star seen on 3 of 80 frames can never anchor a zero point.
        _, _, mag, sig = _synthetic()
        mag[11, :] = np.nan
        mag[11, :3] = 15.0
        sel = ens.select_comps(mag, sig, target_row=0)
        assert sel.role[11] == "field"


class TestStarStats:
    def test_constant_star_chi2_near_one(self):
        M_true, ZP_true, mag, sig = _synthetic(noise=0.02)
        sol = ens.solve_ensemble(mag, sig)
        _, rms, nobs, chi2nu = ens.star_stats(mag, sig, sol.zp)
        # With correct errors, reduced chi2 of constant stars ~ 1 (the
        # weight floor biases it slightly low; allow a generous band).
        assert 0.5 < np.nanmedian(chi2nu) < 1.5
        assert np.all(nobs > 0)


# ---------------------------------------------------------------------------
# Error model
# ---------------------------------------------------------------------------
class TestErrorModel:
    def test_inflation_white_noise_is_one(self):
        rng = np.random.default_rng(1)
        chi = rng.chisquare(50, size=200) / 50   # true chi2nu draws, nu=50
        f = err.inflation_factor(chi)
        assert 0.9 < f < 1.1

    def test_inflation_underestimated_errors(self):
        # Errors underestimated 2x -> chi2nu ~ 4 -> inflation ~ 2.
        assert abs(err.inflation_factor(np.full(10, 4.0)) - 2.0) < 1e-12
        assert np.isnan(err.inflation_factor(np.array([])))

    def test_allan_white_noise_slope(self):
        # White noise integrates down as tau^-1/2: on a doubling ladder
        # each rung falls by sqrt(2).  Fit the log-log slope.
        rng = np.random.default_rng(7)
        y = rng.normal(0, 0.01, size=4096)
        taus, adevs, n = err.allan_deviation(y, dt_s=60.0)
        assert len(taus) >= 8
        slope = np.polyfit(np.log(taus[:8]), np.log(adevs[:8]), 1)[0]
        assert abs(slope + 0.5) < 0.1
        # First rung = the per-point sigma (up to sampling noise).
        assert abs(adevs[0] - 0.01) < 0.002

    def test_allan_too_short(self):
        taus, adevs, n = err.allan_deviation(np.ones(3), 60.0)
        assert len(taus) == 0

    def test_longest_run_splits_on_gaps(self):
        # Two runs separated by a 2-hour gap: the longer one wins.
        jd = np.concatenate([np.arange(5) * 60 / 86400.0,
                             0.5 + np.arange(20) * 60 / 86400.0])
        a, b = err.longest_run(jd, max_gap_s=1800.0)
        assert (a, b) == (5, 25)

    def test_rms_vs_mag_binning(self):
        mags = np.array([15.1, 15.2, 16.1, 16.2, 16.3])
        vals = np.array([0.01, 0.03, 0.1, 0.2, 0.3])
        c, med, n = err.rms_vs_mag_curve(mags, vals, bin_width=1.0)
        assert list(n) == [2, 3]
        assert med[0] == pytest.approx(0.02)
        assert med[1] == pytest.approx(0.2)


# ---------------------------------------------------------------------------
# Gaia geometry
# ---------------------------------------------------------------------------
class TestGaiaGeometry:
    def test_projection_round_trip(self):
        ra = np.array([166.0, 166.2, 165.9])
        dec = np.array([45.0, 45.1, 44.9])
        xi, eta = gg.tangent_project(ra, dec, 166.1, 45.05)
        ra2, dec2 = gg.tangent_deproject(xi, eta, 166.1, 45.05)
        assert np.allclose(ra, ra2, atol=1e-9)
        assert np.allclose(dec, dec2, atol=1e-9)

    def test_projection_scale_is_arcsec(self):
        # 0.1 deg east at the equator of the tangent point = 360 arcsec.
        xi, eta = gg.tangent_project(np.array([10.1]), np.array([0.0]),
                                     10.0, 0.0)
        assert abs(xi[0] - 360.0) < 0.1
        assert abs(eta[0]) < 0.01

    def test_parity_candidates_flip_xi_only(self):
        xi, eta = np.array([1.0, -2.0]), np.array([3.0, 4.0])
        cands = dict(gg.parity_candidates(xi, eta))
        assert np.allclose(cands["direct"][:, 0], xi)
        assert np.allclose(cands["flipped"][:, 0], -xi)
        assert np.allclose(cands["direct"][:, 1], cands["flipped"][:, 1])

    def test_median_offset_robust(self):
        g = np.array([15.0, 16.0, 17.0, 18.0, 30.0])   # one junk match
        m = np.array([10.0, 11.0, 12.0, 13.0, 13.5])
        off, mad, n = gg.median_offset(g, m)
        assert off == pytest.approx(5.0)
        assert n == 5
        off2, _, n2 = gg.median_offset(np.array([np.nan]),
                                       np.array([1.0]))
        assert np.isnan(off2) and n2 == 0


# ---------------------------------------------------------------------------
# Deterministic alignment (seeded RANSAC)
# ---------------------------------------------------------------------------
class TestSeededAlign:
    """Regression for an executed bug: astroalign 2.6.2 shuffles its RANSAC
    with ``np.random.default_rng()`` — fresh OS entropy per call, beyond
    the reach of ``np.random.seed`` — so re-running the match stage flipped
    1-13 star identities per frame and even two back-to-back calls in one
    process disagreed.  ``seeded_ransac`` routes those calls to
    deterministic streams; these tests pin the contract."""

    def _field(self, seed=3):
        # A reference field and the 'same' frame under a small similarity
        # transform + centroid noise — realistic astroalign input.
        rng = np.random.default_rng(seed)
        ref = rng.uniform(0, 4000, size=(200, 2))
        th = 0.003
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        frame = (ref @ R.T) * 1.0015 + [8.0, -5.0]
        frame = frame + rng.normal(0, 0.25, size=frame.shape)
        return frame, ref

    def test_same_seed_same_transform(self):
        from macro_phot import extract as ext
        frame, ref = self._field()
        params = []
        for _ in range(3):
            tf = ext.find_series_transform(frame, ref, seed=12345)
            params.append((tf.scale, tf.rotation,
                           tf.translation[0], tf.translation[1]))
        # Bit-identical across repeated calls: the whole point.
        assert params[0] == params[1] == params[2]
        # A recovered transform, not a degenerate one.
        assert params[0][0] == pytest.approx(1.0015, abs=5e-3)

    def test_different_seeds_still_agree_physically(self):
        # Different seeds explore RANSAC differently but must land on the
        # same physical transform to well below the match tolerance.
        from macro_phot import extract as ext
        frame, ref = self._field()
        t1 = ext.find_series_transform(frame, ref, seed=1)
        t2 = ext.find_series_transform(frame, ref, seed=2)
        assert t1.scale == pytest.approx(t2.scale, abs=1e-3)
        assert t1.rotation == pytest.approx(t2.rotation, abs=1e-3)

    def test_rng_is_restored_even_on_failure(self):
        # The patch must never leak: after the context exits (success OR
        # the MaxIterError path), numpy's default_rng is the original.
        from macro_phot import extract as ext
        orig = np.random.default_rng
        frame, ref = self._field()
        ext.find_series_transform(frame, ref, seed=7)
        assert np.random.default_rng is orig
        with pytest.raises(Exception):
            # Two unrelated point clouds cannot be aligned.
            rng = np.random.default_rng(0)
            ext.find_series_transform(rng.uniform(0, 100, (30, 2)),
                                      rng.uniform(0, 100, (30, 2)) + 5000,
                                      seed=7)
        assert np.random.default_rng is orig

    def test_seed_none_leaves_library_untouched(self):
        # seed=None must not patch anything — library default behavior.
        from macro_phot.extract import seeded_ransac
        orig = np.random.default_rng
        with seeded_ransac(None):
            assert np.random.default_rng is orig
