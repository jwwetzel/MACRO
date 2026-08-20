"""Unit tests for the CV Phase-2 completion arithmetic, plus the invariants
the built product has to satisfy.

Two halves, and they catch different things.

The FIRST half tests ``macro_phot.phase2`` against synthetic data whose right
answer is known by construction: a night with a cloud in the middle of it, a
star field with a known second-order colour term, a similarity transform with
a known rotation, an aperture with a known flux.  These protect the
arithmetic.

The SECOND half tests the PRODUCT.  Every defect this stage found in its own
first run was of a kind no pure-function test could catch: an AIRMASS card of
6,877 that made a fit meaningless while every function did exactly what its
docstring said, and a block that published 66 forced "detections" through a
transform that closed to 1,650 pixels.  Those tests skip when the product is
absent, so a fresh checkout still runs green.
"""

from __future__ import annotations

import math
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import phase2 as p2                            # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PHOT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"


# ===========================================================================
# 1.  Robust statistics
# ===========================================================================
class TestRobustStats:
    def test_mad_of_a_gaussian_recovers_sigma(self):
        rng = np.random.default_rng(7)
        x = rng.normal(10.0, 2.0, 20_000)
        assert p2.median_abs_deviation(x) == pytest.approx(2.0, abs=0.05)

    def test_mad_of_one_value_is_nan_not_zero(self):
        """A spread of one number is UNKNOWN, not zero.  Returning 0.0 would
        make a one-frame night look like the steadiest night of the
        campaign, and the exemplar-night picker sorts on exactly this."""
        assert math.isnan(p2.median_abs_deviation([3.0]))
        assert math.isnan(p2.median_abs_deviation([]))

    def test_running_median_window_shrinks_at_the_ends(self):
        y = [1.0, 2.0, 3.0, 100.0, 5.0, 6.0, 7.0]
        out = p2.running_median(y, 2)
        assert out[0] == pytest.approx(2.0)      # median of [1,2,3]
        assert out[3] == pytest.approx(5.0)      # median of [2,3,100,5,6]
        assert out[-1] == pytest.approx(6.0)     # median of [5,6,7]

    def test_running_median_ignores_nan_and_reports_all_nan_windows(self):
        out = p2.running_median([np.nan, np.nan, np.nan], 0)
        assert np.isnan(out).all()
        out2 = p2.running_median([1.0, np.nan, 3.0], 1)
        assert out2[1] == pytest.approx(2.0)

    def test_mann_whitney_finds_a_shift_and_clears_a_null(self):
        rng = np.random.default_rng(11)
        a = rng.normal(0.0, 1.0, 400)
        b = rng.normal(0.0, 1.0, 400)
        _u, _z, p_null = p2.mann_whitney_u(a, b)
        assert p_null > 0.05
        _u, _z, p_shift = p2.mann_whitney_u(a + 1.0, b)
        assert p_shift < 1e-6

    def test_mann_whitney_handles_ties_without_dividing_by_zero(self):
        u, z, p = p2.mann_whitney_u([1, 1, 1, 1], [1, 1, 1, 1])
        assert math.isnan(z) and math.isnan(p)
        assert u == pytest.approx(8.0)

    def test_mann_whitney_refuses_a_sample_of_one(self):
        """Returning p = 1 for an untestable comparison would read as
        'tested and cleared', which is the one answer it must not give."""
        assert all(math.isnan(v) for v in p2.mann_whitney_u([1.0], [2, 3, 4]))

    def test_two_proportion_z(self):
        d, z, p = p2.two_proportion_z(50, 100, 10, 100)
        assert d == pytest.approx(0.40)
        assert p < 1e-6
        d2, _z2, p2_ = p2.two_proportion_z(10, 100, 10, 100)
        assert d2 == pytest.approx(0.0)
        assert p2_ > 0.9


# ===========================================================================
# 2.  The cloud veto
# ===========================================================================
def _synthetic_night(n_frames=80, n_stars=25, cloud_at=(40, 43),
                     depth=0.55, seed=3):
    """A night of ensemble photometry with a cloud of known depth in it."""
    rng = np.random.default_rng(seed)
    trans = np.ones(n_frames)
    trans[cloud_at[0]:cloud_at[1]] = depth
    star_ids = np.repeat(np.arange(n_stars), n_frames)
    frame_idx = np.tile(np.arange(n_frames), n_stars)
    base = rng.uniform(500.0, 5000.0, n_stars)
    flux = base[star_ids] * trans[frame_idx] * rng.normal(1.0, 0.008,
                                                          star_ids.size)
    return star_ids, frame_idx, flux, trans


class TestCoreEnsemble:
    def test_a_star_present_on_every_frame_qualifies(self):
        sid, fj, _f, _t = _synthetic_night()
        assert len(p2.core_ensemble(sid, fj)) == 25

    def test_a_star_that_comes_and_goes_is_excluded(self):
        """The rule that makes the statistic a measurement of the SKY rather
        than of the detection threshold: if membership can change, a cloud
        shrinks the sum twice — once by dimming, once by bookkeeping."""
        sid = np.array([1, 1, 1, 1, 2, 2, 3])
        fj = np.array([0, 1, 2, 3, 0, 1, 0])
        core = p2.core_ensemble(sid, fj, min_frac=0.8)
        assert core == [1]

    def test_empty_input(self):
        assert p2.core_ensemble([], []) == []


class TestEnsembleFluxRatio:
    def test_ratio_tracks_the_injected_transparency(self):
        sid, fj, flux, trans = _synthetic_night()
        core = p2.core_ensemble(sid, fj)
        ratio, cnt = p2.ensemble_flux_ratio(sid, fj, flux, 80, core)
        assert np.allclose(ratio, trans, atol=0.01)
        assert (cnt == 25).all()

    def test_a_lost_star_does_not_move_the_ratio(self):
        """The subset normalisation is the whole point: a frame that lost
        one core star must still report 'what fraction of its own normal did
        this ensemble deliver', not 'the sum went down by one star'."""
        sid, fj, flux, _t = _synthetic_night(seed=5)
        core = p2.core_ensemble(sid, fj)
        full, _c = p2.ensemble_flux_ratio(sid, fj, flux, 80, core)
        # Drop the BRIGHTEST core star from frame 10 only.
        keep = ~((fj == 10) & (sid == np.argmax(np.bincount(
            sid, weights=flux)[:25])))
        part, cnt = p2.ensemble_flux_ratio(sid[keep], fj[keep], flux[keep],
                                           80, core)
        assert cnt[10] == 24
        assert abs(part[10] - full[10]) < 0.01

    def test_a_frame_with_no_core_star_is_nan_not_zero(self):
        sid = np.array([1, 1, 2, 2])
        fj = np.array([0, 2, 0, 2])          # frame 1 has nothing
        flux = np.array([100.0, 100.0, 200.0, 200.0])
        ratio, cnt = p2.ensemble_flux_ratio(sid, fj, flux, 3, [1, 2])
        assert math.isnan(ratio[1])
        assert cnt[1] == 0


class TestVeto:
    def test_it_finds_the_cloud_and_nothing_else(self):
        sid, fj, flux, _t = _synthetic_night(cloud_at=(40, 43), depth=0.6)
        core = p2.core_ensemble(sid, fj)
        ratio, _c = p2.ensemble_flux_ratio(sid, fj, flux, 80, core)
        vetoed, local = p2.veto_from_ratio(ratio, 0.90)
        assert set(np.nonzero(vetoed)[0]) == {40, 41, 42}
        assert np.isfinite(local).all()

    def test_a_smooth_airmass_decline_is_not_vetoed(self):
        """Extinction across a night is already absorbed exactly by the
        ensemble zero point.  A veto that fired on it would delete the end
        of every night."""
        n = 120
        sid = np.repeat(np.arange(20), n)
        fj = np.tile(np.arange(n), 20)
        trend = np.linspace(1.0, 0.55, n)     # 45% loss across the night
        rng = np.random.default_rng(4)
        base = rng.uniform(800.0, 4000.0, 20)
        flux = base[sid] * trend[fj] * rng.normal(1.0, 0.006, sid.size)
        core = p2.core_ensemble(sid, fj)
        ratio, _c = p2.ensemble_flux_ratio(sid, fj, flux, n, core)
        vetoed, _l = p2.veto_from_ratio(ratio, 0.90)
        assert not vetoed.any()

    def test_a_nan_local_normal_never_vetoes(self):
        vetoed, local = p2.veto_from_ratio([np.nan] * 5, 0.90)
        assert not vetoed.any()
        assert np.isnan(local).all()


class TestZmagCalibration:
    def test_transmission_from_a_dimmed_zero_point(self):
        z = np.full(41, 22.0)
        z[20] = 22.0 - 0.75                  # lost half the light
        t = p2.zmag_transmission(z, half_width=10)
        assert t[20] == pytest.approx(0.5, abs=0.01)
        assert t[0] == pytest.approx(1.0, abs=1e-6)

    def test_zero_zmag_is_absent_not_faint(self):
        """The Sloan-era polar frames write ZMAG = 0.  Reading that as a
        zero point would label every one of them catastrophically
        attenuated."""
        t = p2.zmag_transmission(np.array([22.0, 0.0, 22.0]), half_width=1)
        assert math.isnan(t[1])

    def test_labels_leave_a_deliberate_gap(self):
        lab = p2.label_by_zmag([1.0, 10 ** (-0.4 * 0.08),
                                10 ** (-0.4 * 0.30)])
        assert list(lab) == ["clear", "", "attenuated"]

    def test_roc_and_threshold_choice(self):
        # 100 clear frames at ratio 1.0, 100 attenuated at 0.7.
        rel = [1.0] * 100 + [0.7] * 100
        lab = ["clear"] * 100 + ["attenuated"] * 100
        rows = p2.roc_table(rel, lab, [0.6, 0.8, 0.95, 1.05])
        by_t = {r["threshold"]: r for r in rows}
        assert by_t[0.6]["recall"] == 0.0
        assert by_t[0.8]["recall"] == 1.0
        assert by_t[0.8]["false_veto_rate"] == 0.0
        assert by_t[1.05]["false_veto_rate"] == 1.0
        chosen, why = p2.choose_threshold(rows, max_false_veto=0.01)
        assert chosen == 0.95          # highest inside the budget
        assert "false-veto" in why

    def test_choose_threshold_reports_failure_instead_of_guessing(self):
        rows = p2.roc_table([1.0] * 10, ["clear"] * 10, [1.5])
        chosen, why = p2.choose_threshold(rows, max_false_veto=0.01)
        assert chosen is None
        assert "no threshold" in why


class TestSculptingTest:
    def test_a_veto_uncorrelated_with_brightness_is_cleared(self):
        rng = np.random.default_rng(2)
        mag = rng.normal(19.0, 0.4, 500)
        vetoed = rng.random(500) < 0.1
        out = p2.sculpting_test(mag, vetoed)
        assert out["verdict"] == "NO SCULPTING DETECTED"

    def test_a_veto_that_eats_the_faint_end_is_named(self):
        mag = np.concatenate([np.full(200, 18.0), np.full(200, 21.0)])
        vetoed = np.concatenate([np.zeros(200, bool), np.ones(200, bool)])
        out = p2.sculpting_test(mag, vetoed)
        assert out["verdict"] == "FAINT-PHASE VETO EXCESS"

    def test_a_bright_side_excess_is_not_called_sculpting(self):
        """Two of the 23 real series trip the significance bar in this
        direction.  Calling a bright-side excess 'sculpting' would report
        two alarms that cannot carve a low state out of anything, and would
        train the reader to ignore the word."""
        mag = np.concatenate([np.full(200, 18.0), np.full(200, 21.0)])
        vetoed = np.concatenate([np.ones(200, bool), np.zeros(200, bool)])
        out = p2.sculpting_test(mag, vetoed)
        assert out["verdict"] == "BRIGHT-PHASE VETO EXCESS"

    def test_untestable_cases_say_so(self):
        out = p2.sculpting_test([19.0] * 4, [False] * 4)
        assert out["verdict"].startswith("NOT TESTABLE")


# ===========================================================================
# 3.  Second-order colour extinction
# ===========================================================================
class TestTwoWayCenter:
    def test_row_and_column_means_are_removed(self):
        rng = np.random.default_rng(9)
        n_r, n_c = 12, 15
        row = np.repeat(np.arange(n_r), n_c)
        col = np.tile(np.arange(n_c), n_r)
        v = (rng.normal(0, 1, n_r)[row] + rng.normal(0, 1, n_c)[col]
             + rng.normal(0, 0.01, n_r * n_c))
        out = p2.two_way_center(v, row, col)
        for i in range(n_r):
            assert abs(out[row == i].mean()) < 1e-6
        for j in range(n_c):
            assert abs(out[col == j].mean()) < 1e-6

    def test_it_survives_a_sparse_table(self):
        row = np.array([0, 0, 1, 2, 2])
        col = np.array([0, 1, 0, 1, 2])
        out = p2.two_way_center(np.array([1.0, 2, 3, 4, 5]), row, col)
        assert np.isfinite(out).all()

    def test_empty(self):
        assert p2.two_way_center([], [], []).size == 0


def _synthetic_extinction(kpp=0.030, n_stars=40, n_frames=90, noise=0.010,
                          seed=17):
    rng = np.random.default_rng(seed)
    star = rng.integers(0, n_stars, n_stars * n_frames)
    frame = rng.integers(0, n_frames, n_stars * n_frames)
    colour = rng.uniform(0.2, 1.4, n_stars)[star]
    airmass = rng.uniform(1.0, 2.2, n_frames)[frame]
    truth = kpp * (colour - np.median(colour)) * (airmass -
                                                  np.median(airmass))
    resid = (truth + rng.normal(0, 0.02, n_stars)[star]
             + rng.normal(0, 0.02, n_frames)[frame]
             + rng.normal(0, noise, star.size))
    sigma = np.full(star.size, noise)
    return resid, sigma, colour, airmass, star, frame


class TestFitKpp:
    def test_it_recovers_an_injected_coefficient(self):
        args = _synthetic_extinction(kpp=0.030)
        fit = p2.fit_kpp(*args)
        assert fit.kpp == pytest.approx(0.030, abs=4 * fit.kpp_err)
        assert fit.significant

    def test_it_returns_zero_when_there_is_nothing_there(self):
        args = _synthetic_extinction(kpp=0.0)
        fit = p2.fit_kpp(*args)
        assert abs(fit.t_stat) < p2.KPP_SIGNIFICANCE_T
        assert not fit.significant

    def test_the_star_and_frame_nuisance_terms_do_not_leak_in(self):
        """The injected data carry a random offset per star AND per frame —
        exactly what the Honeycutt solver has already removed.  If the
        design column were not two-way centred, those offsets would bias
        the coefficient."""
        args = _synthetic_extinction(kpp=0.0, seed=31)
        fit = p2.fit_kpp(*args)
        assert abs(fit.kpp) < 5 * fit.kpp_err

    def test_too_few_points_is_reported_not_fitted(self):
        fit = p2.fit_kpp([1.0, 2.0], [0.1, 0.1], [0.5, 0.6], [1.0, 1.1],
                         [0, 1], [0, 1])
        assert math.isnan(fit.kpp)
        assert not fit.significant

    def test_error_is_inflated_when_the_fit_is_over_dispersed(self):
        tight = _synthetic_extinction(kpp=0.03, noise=0.002, seed=41)
        # Same data, but claim a 10x smaller sigma -> chi2nu ~ 100.
        loose = (tight[0], tight[1] / 10.0) + tight[2:]
        f_t, f_l = p2.fit_kpp(*tight), p2.fit_kpp(*loose)
        assert f_l.chi2nu > 10 * f_t.chi2nu
        assert f_l.kpp_err > 0.5 * f_t.kpp_err

    def test_bound_and_sentence(self):
        args = _synthetic_extinction(kpp=0.0)
        fit = p2.fit_kpp(*args)
        b = p2.bound_mmag(fit, 1.2, 1.2)
        assert b > 0
        s = p2.budget_sentence(fit, 1.2, 1.2)
        assert "consistent with zero" in s and "bound" in s


class TestBootstrapKppError:
    def test_the_star_bootstrap_exceeds_the_formal_error(self):
        """The formal error counts POINTS; the real uncertainty is set by
        the number of independent STARS.  On the production data the two
        differ by 3-35x, which is what moved 8 of 10 groups from
        'significant' to 'consistent with zero'."""
        args = _synthetic_extinction(kpp=0.03, n_stars=30, n_frames=120,
                                     seed=23)
        fit = p2.fit_kpp(*args)
        boot = p2.bootstrap_kpp_error(*args, n_boot=16)
        assert math.isfinite(boot)
        assert boot > 0.4 * fit.kpp_err

    def test_it_declines_rather_than_guess_with_too_few_stars(self):
        args = _synthetic_extinction(kpp=0.03, n_stars=6, n_frames=40)
        assert math.isnan(p2.bootstrap_kpp_error(*args, n_boot=16))

    def test_it_is_deterministic(self):
        args = _synthetic_extinction(kpp=0.02, seed=13)
        a = p2.bootstrap_kpp_error(*args, n_boot=8, seed=99)
        b = p2.bootstrap_kpp_error(*args, n_boot=8, seed=99)
        assert a == b


# ===========================================================================
# 4.  Cross-era transformation
# ===========================================================================
class TestOlsLine:
    def test_it_recovers_an_injected_line(self):
        rng = np.random.default_rng(5)
        x = rng.uniform(0.1, 1.3, 300)
        y = 0.02 + 0.05 * (x - 0.7) + rng.normal(0, 0.01, 300)
        fit = p2.ols_line(x, y, np.full(300, 0.01))
        assert fit.b == pytest.approx(0.05, abs=4 * fit.b_err)
        assert fit.a == pytest.approx(0.02 + 0.05 * (fit.x_ref - 0.7),
                                      abs=4 * fit.a_err)

    def test_centring_decorrelates_the_offset_from_the_slope(self):
        """`a` must be the offset AT THE COLOUR THE STARS HAVE, not an
        extrapolation to zero colour that no star in the campaign
        occupies."""
        rng = np.random.default_rng(6)
        x = rng.uniform(1.0, 1.2, 200)      # a narrow, far-from-zero range
        y = 0.5 * x + rng.normal(0, 0.01, 200)
        fit = p2.ols_line(x, y)
        assert fit.x_ref == pytest.approx(np.median(x), abs=1e-9)
        assert fit.a == pytest.approx(0.5 * fit.x_ref, abs=0.01)

    def test_outliers_are_clipped(self):
        x = np.linspace(0, 1, 120)
        y = 0.1 * x.copy()
        y[:4] += 5.0                        # four wild points
        fit = p2.ols_line(x, y, np.full(120, 0.01))
        assert fit.b == pytest.approx(0.1, abs=0.02)
        assert fit.n < 120

    def test_too_few_points(self):
        fit = p2.ols_line([1.0, 2.0], [1.0, 2.0])
        assert math.isnan(fit.b)


class TestFitTransform:
    def test_it_measures_the_bandpass_difference(self):
        rng = np.random.default_rng(8)
        colour = rng.uniform(0.2, 1.5, 250)
        m_from = rng.uniform(14.0, 18.0, 250)
        m_to = m_from + 0.01 + 0.06 * (colour - 0.8) + rng.normal(0, 0.02, 250)
        fit = p2.fit_transform(m_from, m_to, colour, np.full(250, 0.02))
        assert fit.b == pytest.approx(0.06, abs=4 * fit.b_err)


class TestDisciplineChecks:
    def test_era_mixing_is_caught(self):
        n, ex = p2.verify_no_era_mixing([("a|e7|G", 7, 7), ("a|e7|G", 7, 76)])
        assert n == 1
        assert "era 76" in ex[0]

    def test_a_clean_product_reports_zero(self):
        n, ex = p2.verify_no_era_mixing([("a|e7|G", 7, 7)] * 50)
        assert n == 0 and ex == []

    def test_a_transformed_target_magnitude_is_caught(self):
        rows = [("s", 15.0, 17.5, 2.5),              # exactly mag - zp
                ("s", 15.03, 17.5, 2.5)]             # a colour term crept in
        n, ex = p2.verify_no_target_transform(rows)
        assert n == 1
        assert "cal_mag - (mag - zp)" in ex[0]

    def test_nulls_are_skipped_not_failed(self):
        n, _ex = p2.verify_no_target_transform([("s", None, 17.5, 2.5)])
        assert n == 0


# ===========================================================================
# 5.  Forced photometry and upper limits
# ===========================================================================
class TestSimilarity:
    def test_it_recovers_a_known_rotation_scale_and_shift(self):
        rng = np.random.default_rng(1)
        src = rng.uniform(0, 2000, (60, 2))
        th, sc = math.radians(2.5), 1.0008
        m = np.array([[sc * math.cos(th), -sc * math.sin(th)],
                      [sc * math.sin(th), sc * math.cos(th)]])
        dst = src @ m.T + np.array([13.0, -7.5])
        t = p2.similarity_from_pairs(src, dst)
        assert t.scale == pytest.approx(sc, abs=1e-6)
        assert t.rotation_deg == pytest.approx(2.5, abs=1e-4)
        assert t.rms_px < 1e-6
        x, y = p2.apply_similarity(t, 0.0, 0.0)
        assert (x, y) == pytest.approx((13.0, -7.5), abs=1e-6)

    def test_the_residual_reports_a_bad_match_set(self):
        """The per-frame RMS is what the production gate reads.  It has to
        blow up when the pairs disagree, because that is the only signal
        that separates a good forced position from an aperture on blank
        sky."""
        rng = np.random.default_rng(2)
        src = rng.uniform(0, 2000, (40, 2))
        dst = rng.uniform(0, 2000, (40, 2))       # unrelated
        assert p2.similarity_from_pairs(src, dst).rms_px > 100.0

    def test_too_few_pairs(self):
        t = p2.similarity_from_pairs(np.zeros((1, 2)), np.zeros((1, 2)))
        assert math.isnan(t.rms_px)


class TestPlateModel:
    def test_gnomonic_round_trips_through_an_affine(self):
        rng = np.random.default_rng(12)
        ra0, dec0 = 166.4, 25.1
        ra = ra0 + rng.uniform(-0.3, 0.3, 200) / math.cos(math.radians(dec0))
        dec = dec0 + rng.uniform(-0.3, 0.3, 200)
        xi, eta = p2.gnomonic_project(ra, dec, ra0, dec0)
        mat_true = np.array([[7900.0, -60.0, 2394.0],
                             [55.0, 7900.0, 1597.0]])
        xy = (mat_true @ np.vstack([xi, eta, np.ones_like(xi)])).T
        mat, rms = p2.affine_from_pairs(np.column_stack([xi, eta]), xy)
        assert rms < 1e-6
        assert np.allclose(mat, mat_true, atol=1e-6)

    def test_affine_declines_with_too_few_pairs(self):
        mat, rms = p2.affine_from_pairs(np.zeros((2, 2)), np.zeros((2, 2)))
        assert math.isnan(rms) and np.isnan(mat).all()


def _synthetic_image(flux=0.0, sky=100.0, sky_rms=5.0, fwhm=4.0, size=200,
                     seed=21):
    rng = np.random.default_rng(seed)
    img = np.full((size, size), sky) + rng.normal(0, sky_rms, (size, size))
    if flux > 0:
        yy, xx = np.mgrid[0:size, 0:size]
        s = fwhm / 2.3548
        img += (flux * np.exp(-((xx - size / 2) ** 2 + (yy - size / 2) ** 2)
                              / (2 * s ** 2)) / (2 * math.pi * s ** 2))
    return img


class TestForcedAperture:
    def test_it_recovers_an_injected_flux(self):
        img = _synthetic_image(flux=20000.0)
        m = p2.forced_aperture(img, 100.0, 100.0, 10.0, 20.0, 30.0,
                               gain_e_per_adu=1.0)
        assert m.flux == pytest.approx(20000.0, rel=0.05)
        assert m.sky == pytest.approx(100.0, abs=1.0)
        assert m.snr > 50

    def test_blank_sky_gives_a_small_flux_and_an_honest_noise(self):
        img = _synthetic_image(flux=0.0)
        m = p2.forced_aperture(img, 60.0, 60.0, 10.0, 20.0, 30.0)
        assert abs(m.flux) < 3 * m.flux_err
        # Sky-limited aperture noise: sqrt(n_pix) * sky_rms, plus the
        # sky-level term.  Within 25% is the right tolerance for a
        # soft-edged aperture on one noise realisation.
        assert m.flux_err == pytest.approx(math.sqrt(m.n_pix) * m.sky_rms,
                                           rel=0.25)

    def test_the_soft_edge_makes_flux_continuous_in_position(self):
        """A hard pixel mask makes the measured flux jump as the forced
        position moves by a tenth of a pixel — which is the size of the
        production closure residual, so the jump would be real."""
        img = _synthetic_image(flux=20000.0, sky_rms=0.0)
        a = p2.forced_aperture(img, 100.0, 100.0, 10.0, 20.0, 30.0)
        b = p2.forced_aperture(img, 100.1, 100.0, 10.0, 20.0, 30.0)
        assert abs(a.n_pix - b.n_pix) < 1.0

    def test_an_aperture_off_the_edge_is_refused(self):
        img = _synthetic_image()
        m = p2.forced_aperture(img, 1.0, 1.0, 10.0, 20.0, 30.0)
        assert math.isnan(m.flux)


class TestLimits:
    def test_limit_flux_is_k_sigma_not_flux_plus_k_sigma(self):
        """The convention is a statement about the FRAME's sensitivity, not
        about this realisation of the noise.  Both are in use in the
        literature and they differ by tenths of a magnitude."""
        assert p2.limit_flux(50.0, k=3.0) == pytest.approx(150.0)
        assert math.isnan(p2.limit_flux(0.0))
        assert math.isnan(p2.limit_flux(float("nan")))

    def test_limit_magnitude_matches_the_ensemble_convention(self):
        # inst_mag = -2.5 log10(flux/exptime) + 25 ; mag = inst_mag - zp
        m = p2.limit_magnitude(1000.0, 100.0, -0.25)
        assert m == pytest.approx(-2.5 * math.log10(10.0) + 25.0 + 0.25)

    def test_a_non_positive_flux_has_no_magnitude(self):
        assert math.isnan(p2.limit_magnitude(-5.0, 100.0, 0.0))
        assert math.isnan(p2.limit_magnitude(100.0, 0.0, 0.0))
        assert math.isnan(p2.limit_magnitude(100.0, 10.0, None))


class TestSurvival:
    def test_km_with_no_censoring_is_the_empirical_survivor(self):
        m, s = p2.km_survival([1.0, 2.0, 3.0, 4.0], [False] * 4)
        assert list(s) == pytest.approx([0.75, 0.5, 0.25, 0.0])
        assert p2.km_median([1.0, 2.0, 3.0, 4.0], [False] * 4) == 2.0

    def test_censored_points_keep_the_curve_up(self):
        """A right-censored epoch (an upper limit) says 'fainter than this',
        so it removes a subject from the risk set without being an event.
        The curve therefore stays FLAT across a censored magnitude instead
        of stepping down, which is precisely why the limit-aware duty cycle
        differs from the censored one."""
        vals = [1.0, 2.0, 3.0, 4.0]
        _m, s_plain = p2.km_survival(vals, [False] * 4)
        _m, s_cens = p2.km_survival(vals, [False, True, True, False])
        assert list(s_plain) == pytest.approx([0.75, 0.5, 0.25, 0.0])
        # Two censored epochs in the middle: no event, so no step.
        assert list(s_cens) == pytest.approx([0.75, 0.75, 0.75, 0.0])
        assert s_cens[2] > s_plain[2]

    def test_no_estimable_median_returns_nan(self):
        """More than half the epochs being limits is a real state of this
        data set, and inventing a median from the detections alone is
        exactly the bias this task removes."""
        assert math.isnan(p2.km_median([1.0, 2.0, 3.0],
                                       [False, True, True]))

    def test_state_statistics_report_both_versions(self):
        st = p2.state_statistics([15.0, 15.2, 15.5],
                                 [16.0, 16.2, 15.8, 16.5])
        assert st["n_detected"] == 3
        assert st["n_limit"] == 4
        assert st["detected_fraction_censored"] == 1.0
        assert st["detected_fraction_true"] == pytest.approx(3 / 7)
        assert st["faint_state_fraction"] == pytest.approx(4 / 7)
        assert st["median_censored"] == 15.2

    def test_state_statistics_on_an_empty_series(self):
        st = p2.state_statistics([], [])
        assert st["n_epochs"] == 0
        assert math.isnan(st["detected_fraction_true"])


# ===========================================================================
# 6.  Constants that are claims
# ===========================================================================
class TestConstants:
    def test_the_airmass_window_excludes_the_card_that_broke_the_first_run(self):
        """62 matched CV frames carry AIRMASS between 5 and 6,877.  VV Pup
        cannot exceed 2.1 from this site.  The window has to refuse them."""
        assert p2.AIRMASS_MIN == 1.0
        assert p2.AIRMASS_MAX <= 3.0
        assert not (p2.AIRMASS_MIN <= 6877.0 <= p2.AIRMASS_MAX)

    def test_the_significance_bar_is_three_sigma(self):
        """Two sigma over ~10 (era, filter) fits manufactures a significant
        result by arithmetic alone."""
        assert p2.KPP_SIGNIFICANCE_T >= 3.0

    def test_the_detection_and_limit_thresholds_are_the_same_number(self):
        """'We can see it' and 'we can bound it' must be one threshold, or
        the two categories overlap."""
        assert p2.FORCED_DETECT_SNR == p2.LIMIT_SIGMA

    def test_inst_mag_offset_matches_the_photometry_module(self):
        from macro_phot import photometry as ph
        assert p2.INST_MAG_OFFSET == ph.INST_MAG_OFFSET


# ===========================================================================
# 7.  The built product
# ===========================================================================
def _con(path: Path) -> sqlite3.Connection:
    if not path.exists():
        pytest.skip(f"{path.name} not built in this checkout")
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


@pytest.fixture(scope="module")
def phot():
    con = _con(PHOT_DB)
    yield con
    con.close()


def _has(con, table: str) -> bool:
    return bool(con.execute("SELECT count(*) FROM sqlite_master WHERE "
                            "type='table' AND name=?", (table,)).fetchone()[0])


class TestCloudProduct:
    def test_the_veto_flag_agrees_with_the_recorded_threshold(self, phot):
        """A stored flag that no longer follows from the stored threshold is
        the defect class this repo has already been bitten by once (the
        saturation veto).  Recompute it."""
        if not _has(phot, "p2_cloud_frame"):
            pytest.skip("cloud stage not run")
        thr = phot.execute("SELECT value FROM p2_meta WHERE "
                           "key='cloud_threshold_used'").fetchone()
        if thr is None:
            pytest.skip("no threshold recorded")
        bad = phot.execute(
            "SELECT count(*) FROM p2_cloud_frame WHERE vetoed <> "
            "(CASE WHEN rel_ratio IS NOT NULL AND rel_ratio < ? THEN 1 "
            "ELSE 0 END)", (float(thr[0]),)).fetchone()[0]
        assert bad == 0

    def test_no_series_shows_a_faint_side_veto_excess(self, phot):
        """The one result that would invalidate the veto: a cut that eats
        the target's faint phases sculpts every light curve behind it."""
        if not _has(phot, "p2_cloud_bias"):
            pytest.skip("cloud stage not run")
        bad = phot.execute(
            "SELECT group_concat(series_key) FROM p2_cloud_bias "
            "WHERE verdict='FAINT-PHASE VETO EXCESS'").fetchone()[0]
        assert not bad, f"faint-phase veto excess in: {bad}"

    def test_the_calibration_rested_on_independent_evidence(self, phot):
        if not _has(phot, "p2_cloud_roc"):
            pytest.skip("cloud stage not run")
        n_clear, n_att = phot.execute(
            "SELECT max(n_clear), max(n_attenuated) FROM "
            "p2_cloud_roc").fetchone()
        assert (n_clear or 0) >= 100
        assert (n_att or 0) >= 20


class TestExtinctionProduct:
    def test_no_fit_used_an_impossible_airmass(self, phot):
        if not _has(phot, "p2_extinction"):
            pytest.skip("extinction stage not run")
        bad = phot.execute(
            "SELECT group_concat(era_id || ' ' || filter) FROM p2_extinction "
            "WHERE airmass_max > ? OR airmass_min < ?",
            (p2.AIRMASS_MAX, p2.AIRMASS_MIN)).fetchone()[0]
        assert not bad, f"impossible airmass entered the fit for: {bad}"

    def test_the_published_error_is_never_smaller_than_the_formal_one(self,
                                                                     phot):
        if not _has(phot, "p2_extinction"):
            pytest.skip("extinction stage not run")
        bad = phot.execute(
            "SELECT count(*) FROM p2_extinction WHERE kpp_err IS NOT NULL "
            "AND kpp_err_formal IS NOT NULL AND kpp_err < kpp_err_formal "
            "- 1e-12").fetchone()[0]
        assert bad == 0

    def test_significance_follows_from_the_published_error(self, phot):
        if not _has(phot, "p2_extinction"):
            pytest.skip("extinction stage not run")
        for era, filt, k, e, sig in phot.execute(
                "SELECT era_id, filter, kpp, kpp_err, significant "
                "FROM p2_extinction WHERE kpp IS NOT NULL"):
            expect = abs(k / e) >= p2.KPP_SIGNIFICANCE_T
            assert bool(sig) == expect, f"era {era} {filt}"


class TestCrossEraProduct:
    def test_every_discipline_assertion_holds(self, phot):
        if not _has(phot, "p2_discipline"):
            pytest.skip("crossera stage not run")
        bad = phot.execute(
            "SELECT group_concat(check_id) FROM p2_discipline "
            "WHERE verdict NOT IN ('HOLDS','NOT APPLICABLE')").fetchone()[0]
        assert not bad, f"discipline violated: {bad}"

    def test_a_verdict_of_holds_never_carries_violations(self, phot):
        if not _has(phot, "p2_discipline"):
            pytest.skip("crossera stage not run")
        bad = phot.execute("SELECT count(*) FROM p2_discipline WHERE "
                           "verdict='HOLDS' AND n_violation > 0").fetchone()[0]
        assert bad == 0

    def test_no_transformation_was_applied_to_a_target(self, phot):
        if not _has(phot, "p2_transform"):
            pytest.skip("crossera stage not run")
        bad = phot.execute("SELECT count(*) FROM p2_transform WHERE "
                           "applied_to_targets <> 0").fetchone()[0]
        assert bad == 0

    def test_each_slope_agrees_with_the_tie_colour_terms(self, phot):
        """Two independent routes to the same bandpass difference.  A
        tension above 4 sigma would mean one of them is measuring something
        else."""
        if not _has(phot, "p2_transform"):
            pytest.skip("crossera stage not run")
        bad = phot.execute(
            "SELECT group_concat(target_key || ' ' || band_from) "
            "FROM p2_transform WHERE b_tension_sigma IS NOT NULL "
            "AND abs(b_tension_sigma) > 4.0").fetchone()[0]
        assert not bad, f"slope disagrees with the tie for: {bad}"


class TestLimitsProduct:
    def test_every_published_limit_came_from_a_validated_position(self, phot):
        """The gate that refused EU UMa's era-78 block.  A limit measured at
        an unverifiable position is not a limit."""
        if not _has(phot, "p2_limit_series"):
            pytest.skip("forced stage not run")
        bad = phot.execute(
            "SELECT group_concat(series_key) FROM p2_limit_series "
            "WHERE n_limits > 0 AND (n_closure < ? OR closure_median_px "
            "IS NULL OR closure_median_px > ?)",
            (p2.CLOSURE_MIN_FRAMES, p2.CLOSURE_MAX_MEDIAN_PX)).fetchone()[0]
        assert not bad, f"limits published without closure for: {bad}"

    def test_no_limit_row_carries_a_magnitude_without_a_zero_point(self,
                                                                  phot):
        if not _has(phot, "p2_limits"):
            pytest.skip("forced stage not run")
        bad = phot.execute(
            "SELECT count(*) FROM p2_limits WHERE limit_mag IS NOT NULL "
            "AND zp IS NULL").fetchone()[0]
        assert bad == 0

    def test_the_limit_magnitude_is_reproducible_from_its_own_columns(self,
                                                                     phot):
        if not _has(phot, "p2_limits"):
            pytest.skip("forced stage not run")
        rows = phot.execute(
            "SELECT limit_flux_adu, exptime, zp, limit_mag FROM p2_limits "
            "WHERE limit_mag IS NOT NULL LIMIT 400").fetchall()
        if not rows:
            pytest.skip("no limits recorded")
        for f, t, zp, m in rows:
            assert p2.limit_magnitude(f, t, zp) == pytest.approx(m, abs=1e-9)

    def test_every_forced_outcome_is_one_of_the_declared_four(self, phot):
        if not _has(phot, "p2_limits"):
            pytest.skip("forced stage not run")
        kinds = {r[0] for r in phot.execute(
            "SELECT DISTINCT outcome FROM p2_limits")}
        assert kinds <= {"limit", "detection", "failed", "no_zeropoint"}

    def test_the_stage_wrote_nothing_into_the_light_curve(self, phot):
        """The design commitment that makes both the sculpting test and the
        'transformation is metadata' claim falsifiable: this stage adds no
        column to cv_lightcurve."""
        cols = {r[1] for r in phot.execute(
            "PRAGMA table_info(cv_lightcurve)")}
        assert not (cols & {"vetoed", "cloud_flag", "limit_mag",
                            "forced_mag", "p2_mag"})
