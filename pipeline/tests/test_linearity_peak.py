"""Unit tests for the peak-resolved differential linearity estimators
(rlmt_diagnostics.linearity, second half) and the star photometry and
saturation logic that feed them (rlmt_diagnostics.starphot, .saturation).

The estimator's whole claim is that it measures response deviation against
a star's own peak WITHOUT being fooled by what changes from frame to frame.
So the tests are built around known truth: a linear detector must come out
linear (the null test, with its SIGNED bias checked against zero), a known
roll-off must come out at its injected size, and the two defects found
while the estimator was written — a selection bias and an unconverged
normalisation — each have a test that fails if they return.
"""

from __future__ import annotations

import numpy as np
import pytest

from rlmt_diagnostics import linearity as lin
from rlmt_diagnostics import saturation as sat
from rlmt_diagnostics import starphot as sp

CEIL = 3400.0                     # usable scale in ADU (ceiling - bias)


def ladder_geometry(n=60, seed=1):
    """An exposure ladder on one field: 60 stars, 8/32/64/128 s, with
    frame factors (transparency) that differ by up to 12%."""
    rng = np.random.default_rng(seed)
    t = np.array([8.0, 32.0, 64.0, 128.0])
    z = np.array([1.0, 0.97, 1.02, 0.9])
    z = z / np.median(z)
    peak_rate = 10 ** rng.uniform(0.3, 2.2, n)       # peak ADU per second
    base = peak_rate * 40.0                          # flux ADU per second
    true = base[:, None] * t[None, :] * z[None, :]
    err = np.sqrt(true / 1.06 + 6400.0)              # Poisson + sky floor
    ppf = np.full(true.shape, 1.0 / 40.0 / CEIL)     # peak fraction per flux
    return {"t": t, "z": z, "base": base, "true": true, "err": err,
            "ppf": ppf}


def rolloff(pf):
    return np.where(pf > 0.7, -0.05 * (pf - 0.7) / 0.3, 0.0)


class TestPeakFraction:
    def test_definition(self):
        assert lin.peak_fraction(3496.0, 94.0, 3496.0) == pytest.approx(1.0)
        assert lin.peak_fraction(94.0, 94.0, 3496.0) == pytest.approx(0.0)
        assert lin.peak_fraction(1795.0, 94.0, 3496.0) == pytest.approx(0.5)


class TestEnsembleDeviation:
    def test_recovers_frame_factors_and_base_rates(self):
        g = ladder_geometry()
        pf = g["true"] * g["ppf"]
        res = lin.ensemble_deviation(g["true"], g["err"], g["t"], pf)
        assert res["z"] == pytest.approx(g["z"], rel=1e-6)
        ok = np.isfinite(res["base"])
        assert res["base"][ok] == pytest.approx(g["base"][ok], rel=1e-6)
        assert np.nanmax(np.abs(res["dev"])) < 1e-6

    def test_noiseless_rolloff_is_read_back_exactly(self):
        g = ladder_geometry()
        pf = g["true"] * g["ppf"]
        flux = g["true"] * (1.0 + rolloff(pf))
        res = lin.ensemble_deviation(flux, g["err"], g["t"], pf)
        test = np.isfinite(res["dev"]) & (pf > 0.75)
        assert test.sum() > 5
        assert res["dev"][test] == pytest.approx(rolloff(pf)[test], abs=1e-4)

    def test_a_pure_scale_change_is_invisible(self):
        """Honesty limit, pinned: a response that is wrong by the SAME
        factor at every level is absorbed in the frame factors."""
        g = ladder_geometry()
        pf = g["true"] * g["ppf"]
        res = lin.ensemble_deviation(0.9 * g["true"], g["err"], g["t"], pf)
        assert np.nanmax(np.abs(res["dev"])) < 1e-6

    def test_star_never_in_reference_regime_is_not_anchored(self):
        g = ladder_geometry()
        pf = g["true"] * g["ppf"]
        flux = g["true"].copy()
        bright = np.argmax(g["base"])
        pf[bright, :] = np.maximum(pf[bright, :], 0.5)
        res = lin.ensemble_deviation(flux, g["err"], g["t"], pf)
        assert np.isnan(res["base"][bright])
        assert np.all(np.isnan(res["dev"][bright]))

    def test_reference_frames_can_be_restricted(self):
        """StackPro-as-reference: only flagged frames may anchor."""
        g = ladder_geometry()
        pf = g["true"] * g["ppf"]
        only_first = np.array([True, False, False, False])
        res = lin.ensemble_deviation(g["true"], g["err"], g["t"], pf,
                                     frame_is_reference=only_first)
        assert not res["is_ref"][:, 1:].any()

    def test_missing_measurements_are_tolerated(self):
        g = ladder_geometry()
        pf = g["true"] * g["ppf"]
        flux = g["true"].copy()
        flux[::5, 2] = np.nan
        res = lin.ensemble_deviation(flux, g["err"], g["t"], pf)
        assert res["z"] == pytest.approx(g["z"], rel=1e-6)
        assert np.all(np.isnan(res["dev"][::5, 2]))


class TestInjectionBias:
    """Standing rule 3: the SIGNED matched-cell bias, never |bias|."""

    def test_null_injection_is_unbiased(self):
        g = ladder_geometry()
        rows = lin.injection_bias(g["err"], g["t"], g["base"], g["z"],
                                  g["ppf"], lambda pf: 0.0 * pf,
                                  n_trials=60)
        well = [r for r in rows if r["n_trials"] >= 50 and r["lo"] >= 0.3]
        assert len(well) >= 4
        for r in well:
            assert abs(r["bias_pct"]) < 0.25, r
            assert abs(r["bias_pct"]) < 4 * r["bias_err_pct"] + 0.05, r

    def test_rolloff_is_recovered_at_its_injected_size(self):
        g = ladder_geometry()
        rows = lin.injection_bias(g["err"], g["t"], g["base"], g["z"],
                                  g["ppf"], rolloff, n_trials=60)
        top = [r for r in rows if r["lo"] >= 0.95 and r["n_trials"] >= 40]
        assert top, "the ladder must populate the top of the scale"
        for r in top:
            assert r["injected_pct"] < -4.0
            assert abs(r["bias_pct"]) < 0.4, r

    def test_selection_on_measured_snr_would_be_biased(self):
        """Regression test for the first defect: choosing reference
        measurements by their OWN measured signal-to-noise lets through the
        ones noise pushed up.  Emulate that selection by hand and show the
        reference fluxes it admits are biased high — the mechanism the
        expected-S/N cut in ensemble_deviation exists to avoid."""
        rng = np.random.default_rng(5)
        true = np.full(200_000, 30.0)          # exactly at the threshold
        meas = true + rng.normal(0.0, 1.0, true.size)
        selected = meas[meas >= 30.0]
        assert selected.mean() > true[0] + 0.5

    def test_normalisation_is_iterated_to_convergence(self):
        """Regression test for the second defect: six alternating
        iterations left the frame factors ~1% from their fixed point when
        the frames' reference sets barely overlap."""
        g = ladder_geometry()
        rng = np.random.default_rng(6)
        flux = g["true"] + rng.normal(0, 1, g["true"].shape) * g["err"]
        pf = g["true"] * g["ppf"]
        res = lin.ensemble_deviation(flux, g["err"], g["t"], pf)
        assert res["z"] == pytest.approx(g["z"], rel=0.01)
        assert lin.NORM_ITERS >= 500


class TestDeviationCurve:
    def test_groups_are_averaged_before_binning(self):
        """Forty points from one star are one independent unit, not forty."""
        pf = np.concatenate([np.full(40, 0.55), np.full(3, 0.55)])
        dev = np.concatenate([np.full(40, -0.02), np.array([0.0, 0.0, 0.0])])
        grp = np.concatenate([np.zeros(40), np.array([1, 2, 3])])
        curve = lin.deviation_curve(pf, dev, grp)
        b = next(c for c in curve if c["lo"] == 0.5)
        assert b["n_groups"] == 4 and b["n_points"] == 43
        assert b["dev_pct"] == pytest.approx(0.0)      # median of 4 groups

    def test_thin_bin_is_reported_but_not_measured(self):
        curve = lin.deviation_curve([0.85, 0.86], [-0.01, -0.02], [0, 1])
        assert len(curve) == 1 and curve[0]["measured"] is False

    def test_errors_come_from_group_scatter(self):
        rng = np.random.default_rng(7)
        dev = rng.normal(0.0, 0.01, 400)
        curve = lin.deviation_curve(np.full(400, 0.45), dev, np.arange(400))
        assert curve[0]["dev_err_pct"] == pytest.approx(
            100 * 1.2533 * 0.01 / 20.0, rel=0.25)


def flat_curve(devs, errs=None, measured=None):
    edges = lin.PEAK_BIN_EDGES
    errs = errs or [0.1] * len(devs)
    measured = measured or [True] * len(devs)
    return [{"lo": edges[i], "hi": edges[i + 1],
             "peak_frac": 0.5 * (edges[i] + edges[i + 1]),
             "dev_pct": d, "dev_err_pct": e, "n_points": 50,
             "n_groups": 10, "measured": m}
            for i, (d, e, m) in enumerate(zip(devs, errs, measured))]


class TestRecommendCap:
    def test_cap_stops_at_first_violation(self):
        #        0-.1 .1-.2 .2-.3 .3-.4 .4-.5 .5-.6 .6-.7 .7-.8 .8-.85
        curve = flat_curve([0, 0, 0, -0.1, -0.2, -0.3, -0.6, -0.9, -1.6])
        rec = lin.recommend_cap(curve)
        assert rec["cap_fraction"] == pytest.approx(0.80)
        assert rec["limited_by"] == "nonlinearity"
        assert rec["first_bad"]["lo"] == pytest.approx(0.80)
        assert rec["worst_dev_pct"] == pytest.approx(-0.9)

    def test_data_limited_cap_says_so(self):
        curve = flat_curve([0, 0, 0, 0.1, -0.1, 0.2])
        rec = lin.recommend_cap(curve)
        assert rec["cap_fraction"] == pytest.approx(0.60)
        assert rec["limited_by"] == "data"
        assert rec["first_bad"] is None

    def test_a_gap_is_not_assumed_linear(self):
        m = [True, True, True, True, False, True, True]
        curve = flat_curve([0, 0, 0, 0.1, 0.0, 0.1, 0.1], measured=m)
        rec = lin.recommend_cap(curve)
        assert rec["cap_fraction"] == pytest.approx(0.40)
        assert rec["limited_by"] == "data"

    def test_nonlinear_from_the_first_tested_bin(self):
        rec = lin.recommend_cap(flat_curve([0, 0, 0, -2.5]))
        assert rec["cap_fraction"] is None
        assert rec["first_bad"]["dev_pct"] == pytest.approx(-2.5)

    def test_slope_is_fitted_over_the_bins_below_the_cap(self):
        devs = [0, 0, 0, -0.2, -0.4, -0.6, -0.8]
        rec = lin.recommend_cap(flat_curve(devs))
        assert rec["slope_pct_per_scale"] < 0
        assert rec["slope_err"] > 0

    def test_no_curve_no_cap(self):
        rec = lin.recommend_cap([])
        assert rec["cap_fraction"] is None and rec["limited_by"] == "none"


# ---------------------------------------------------------------------------
# starphot
# ---------------------------------------------------------------------------
def star_field(seed=0, n=50, fwhm=4.7, shift=(0.0, 0.0), size=768):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:size, 0:size]
    img = np.zeros((size, size))
    stars = [(rng.uniform(90, size - 90), rng.uniform(90, size - 90),
              10 ** rng.uniform(3.2, 5.3)) for _ in range(n)]
    sig = fwhm / 2.3548
    for x, y, f in stars:
        img += f / (2 * np.pi * sig * sig) * np.exp(
            -((xx - x - shift[0]) ** 2 + (yy - y - shift[1]) ** 2)
            / (2 * sig * sig))
    noise_rng = np.random.default_rng(seed + 1000)
    return img + noise_rng.normal(100.0, 4.0, img.shape), stars


class TestStarphot:
    def test_fluxes_and_fwhm_are_recovered(self):
        img, stars = star_field()
        m = sp.measure_stars(img)
        assert m["fwhm_med"] == pytest.approx(4.7, rel=0.08)
        truth = np.array([(x, y) for x, y, _f in stars])
        tflux = np.array([f for _x, _y, f in stars])
        meas_xy = np.c_[m["x"], m["y"]]
        idx = sp.match_stars(truth, meas_xy, (0.0, 0.0))
        ok = idx >= 0
        assert ok.sum() >= 40
        k = m["apertures"].index(sp.choose_aperture(m["fwhm_med"]))
        ratio = m["flux"][k][idx[ok]] / tflux[ok]
        assert np.median(ratio) == pytest.approx(1.0, abs=0.01)

    def test_peak_is_the_raw_peak_including_background(self):
        img, stars = star_field(n=5)
        m = sp.measure_stars(img)
        j = int(np.argmax(m["flux"][0]))
        x, y = int(round(m["x"][j])), int(round(m["y"][j]))
        assert m["peak_raw"][j] == pytest.approx(
            img[y - 3:y + 4, x - 3:x + 4].max())
        assert m["peak_raw"][j] > 100.0          # background is in there

    def test_shift_is_found_and_stars_matched(self):
        img_a, _ = star_field(seed=3)
        img_b, _ = star_field(seed=3, shift=(13.0, -7.0))
        a, b = sp.measure_stars(img_a), sp.measure_stars(img_b)
        shift = sp.match_shift(np.c_[a["x"], a["y"]], np.c_[b["x"], b["y"]])
        assert shift == pytest.approx((-13.0, 7.0), abs=0.5)
        idx = sp.match_stars(np.c_[a["x"], a["y"]], np.c_[b["x"], b["y"]],
                             shift)
        assert (idx >= 0).mean() > 0.85
        assert len(set(idx[idx >= 0])) == (idx >= 0).sum()   # one-to-one

    def test_no_common_stars_no_shift(self):
        assert sp.match_shift(np.zeros((2, 2)), np.zeros((2, 2))) is None

    def test_empty_frame(self):
        m = sp.measure_stars(np.full((300, 300), 100.0)
                             + np.random.default_rng(0).normal(0, 3, (300, 300)))
        assert len(m["x"]) == 0 and m["flux"].shape[1] == 0

    def test_choose_aperture(self):
        assert sp.choose_aperture(2.0) == 6.0
        assert sp.choose_aperture(4.0) == 12.0
        assert sp.choose_aperture(50.0) == max(sp.APERTURES_PX)

    def test_native_peak_factor_matches_the_engineer_memo(self):
        """DE.F9: an average-binned peak understates the native peak by
        7-16% for FWHM 6-4 native pixels (3-2 binned pixels)."""
        assert sp.native_peak_factor(3.0) == pytest.approx(1.07, abs=0.02)
        assert sp.native_peak_factor(2.0) == pytest.approx(1.17, abs=0.03)
        # Broader stars hide less; an unresolved one hides the most.
        f = [sp.native_peak_factor(w) for w in (1.5, 2.0, 3.0, 5.0)]
        assert f == sorted(f, reverse=True) and f[-1] > 1.0


# ---------------------------------------------------------------------------
# saturation
# ---------------------------------------------------------------------------
def stamp_with_star(peak=2000.0, fwhm=4.0, sky=100.0, clip=None, seed=0):
    rng = np.random.default_rng(seed)
    n = 2 * sat.STAMP_HALF + 1
    yy, xx = np.mgrid[0:n, 0:n]
    sig = fwhm / 2.3548
    img = sky + peak * np.exp(-((xx - 20) ** 2 + (yy - 20) ** 2)
                              / (2 * sig * sig))
    img = img + rng.normal(0, 3.0, img.shape)
    return np.minimum(img, clip) if clip else img


class TestSaturation:
    def test_unclipped_star(self):
        st = sat.stamp_peak_stats(stamp_with_star())
        assert st["peak_raw"] == pytest.approx(2100.0, abs=15.0)
        assert st["n_at_peak"] <= 2
        assert st["fwhm_px"] == pytest.approx(4.0, rel=0.15)

    def test_clipped_star_is_a_plateau(self):
        st = sat.stamp_peak_stats(stamp_with_star(peak=9000.0, clip=3500.0))
        assert st["peak_raw"] == pytest.approx(3500.0)
        assert st["n_at_peak"] >= sat.FLAT_TOP_MIN_PIXELS

    def test_verdict_ladder(self):
        v = sat.saturation_verdict
        assert v(3490.0, 94.0, 3496.0, 0.8, 1) == "clipped"
        assert v(3100.0, 94.0, 3496.0, 0.8, 6) == "flat_topped"
        assert v(3100.0, 94.0, 3496.0, 0.8, 1) == "above_cap"
        assert v(1500.0, 94.0, 3496.0, 0.8, 1) == "clean"
        # No cap for the mode: only clipping can be judged.
        assert v(3100.0, 94.0, 3496.0, None, 1) == "clean"

    def test_native_equivalent_peak(self):
        assert sat.native_equivalent_peak(3000.0, 94.0, 4.0, 1) == (3000.0, 1.0)
        native, factor = sat.native_equivalent_peak(40000.0, 303.0, 2.0, 2)
        assert factor == pytest.approx(sp.native_peak_factor(2.0))
        assert native == pytest.approx(303.0 + (40000.0 - 303.0) * factor)
        # A missing width falls back to the sharp (conservative) case.
        assert sat.native_equivalent_peak(40000.0, 303.0, float("nan"), 2)[1] \
            == pytest.approx(sp.native_peak_factor(2.0))

    def test_average_binning_hides_a_clipped_native_pixel(self):
        """The point of rule 4: a binned peak of 58,000 ADU passes a
        60,200 veto while the native pixel inside it is at the rail."""
        native, _ = sat.native_equivalent_peak(58000.0, 303.0, 2.0, 2)
        assert native > 65535.0
        assert sat.saturation_verdict(native, 303.0, 65535.0, 0.9, 1) \
            == "clipped"

    def test_fallback_locator(self):
        x = np.array([10.0, 500.0, 520.0])
        y = np.array([10.0, 510.0, 480.0])
        flux = np.array([9e6, 2e4, 5e4])
        # The corner star is the brightest but is not near the centre.
        assert sat.brightest_central_source(x, y, flux, (1024, 1024)) == 2
        assert sat.brightest_central_source(x[:1], y[:1], flux[:1],
                                            (1024, 1024)) is None
        assert sat.brightest_central_source([], [], [], (10, 10)) is None
