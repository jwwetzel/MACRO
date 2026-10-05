"""Unit tests for S3b — the absolute clock check from archived transits.

Layout mirrors ``macro_core.clock_transits`` and the database-facing parts
of ``scripts/build_s3b_clock_transits.py``:

* the occultation model, pinned three independent ways (the analytic
  two-circle overlap for a uniform disc at every radius ratio the stage
  uses; the small-planet limit under limb darkening; the total-eclipse
  limit), plus its symmetry — the one property a mid-time relies on;
* ephemeris propagation and the clock-era map;
* the light-curve assembly rules (native-pixel saturation ratio, the
  model-free scatter estimate, comparison pruning, the cloud cut);
* the mid-time fit: unbiased on synthetic data, BLIND to where the event
  is (a 1,065 s displacement is found without being told), refuses a
  one-sided event, and has a bootstrap error that matches the scatter of
  repeated experiments;
* injection-recovery: signed bias, and an honest "not recoverable" when
  the shift pushes the event out of the window;
* the combination rules (no chi-square rescaling) and the verdict;
* the build script's census arithmetic, series keys, saturation handling
  and the twin (raw vs reduced) rule of the summary, on tiny synthetic
  databases.

No test touches the archive, the manifest or the network.
"""

from __future__ import annotations

import math
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_core import clock_transits as ct          # noqa: E402
from scripts import build_s3b_clock_transits as b    # noqa: E402

DAY = 86400.0


def _transit_shape(**kw) -> ct.EventShape:
    """A WASP-43 b-like transit (P 0.81 d, a/R 4.87, k 0.159, b 0.66)."""
    base = dict(period=0.8135, a_over_r=4.87, k=0.159, b=0.66, u1=0.45,
                u2=0.25)
    base.update(kw)
    return ct.EventShape(**base)


def _synthetic(shape, t0, n=300, span=(-0.08, 0.09), noise=0.003, seed=1,
               exptime_s=30.0, slope=0.02):
    """One noisy synthetic light curve with a tilted baseline."""
    rng = np.random.default_rng(seed)
    t = np.linspace(*span, n) + 2460000.0
    ex = exptime_s / DAY
    f = ct.event_model(t, ex, t0, shape.period, shape.a_over_r, shape.k,
                       shape.b, shape.u1, shape.u2, shape.f1)
    f = f * (1.0 + slope * (t - t.mean())) + rng.normal(0, noise, n)
    return t, f, np.full(n, noise), ex


# ---------------------------------------------------------------------------
# The occultation model
# ---------------------------------------------------------------------------
class TestOccultation:
    @pytest.mark.parametrize("k", [0.05, 0.159, 0.5, 0.916, 1.0, 2.0, 9.1])
    def test_uniform_disc_matches_analytic_overlap(self, k):
        # Independent truth: the lens-area formula.  Covers a planet, an
        # HW Vir secondary (k ~ 1) and an M dwarf over a white dwarf.
        z = np.linspace(0.0, 1.5 + k, 3001)
        got = ct.occultation_depth(z, k)
        want = ct.uniform_overlap_depth(z, k)
        assert np.max(np.abs(got - want)) < 2e-6

    def test_no_overlap_is_zero_and_total_is_one(self):
        assert ct.occultation_depth(np.array([1.2]), 0.159)[0] == 0.0
        # An occulter larger than the star, centred: everything hidden,
        # limb darkening or not.
        assert ct.occultation_depth(np.array([0.0]), 9.1, 0.3, 0.2)[0] \
            == pytest.approx(1.0, abs=1e-9)

    def test_small_planet_limit_under_limb_darkening(self):
        # A tiny planet at disc centre hides k^2 x I(centre)/<I>, and
        # I(centre) = 1 while <I> = 1 - u1/3 - u2/6.
        k, u1, u2 = 0.01, 0.4, 0.2
        want = k * k / (1.0 - u1 / 3.0 - u2 / 6.0)
        got = ct.occultation_depth(np.array([0.0]), k, u1, u2)[0]
        assert got == pytest.approx(want, rel=2e-4)

    def test_limb_darkening_deepens_centre_and_shallows_limb(self):
        k = 0.1
        flat = ct.occultation_depth(np.array([0.0, 0.95]), k)
        dark = ct.occultation_depth(np.array([0.0, 0.95]), k, 0.5, 0.2)
        assert dark[0] > flat[0] and dark[1] < flat[1]

    def test_depth_is_monotonic_in_separation(self):
        z = np.linspace(0, 1.3, 400)
        d = ct.occultation_depth(z, 0.159, 0.45, 0.25)
        assert np.all(np.diff(d) <= 1e-12)

    def test_projected_separation_geometry(self):
        z = ct.projected_separation(np.array([10.0, 10.1, 9.9, 10.5]),
                                    10.0, 1.0, 5.0, 0.3)
        assert z[0] == pytest.approx(0.3)            # conjunction -> b
        assert z[1] == pytest.approx(z[2])           # even about t0
        assert math.isinf(z[3])                      # far side: no event

    def test_duration_formula_against_brute_force(self):
        s = _transit_shape()
        t = np.linspace(-0.05, 0.05, 200001)
        z = ct.projected_separation(t, 0.0, s.period, s.a_over_r, s.b)
        inside = t[z < 1.0 + s.k]
        assert s.duration_d() == pytest.approx(inside.max() - inside.min(),
                                               abs=2e-6)
        # A geometry that never touches has no duration.
        assert ct.event_duration_d(1.0, 5.0, 0.1, 1.2) == 0.0

    def test_model_is_symmetric_about_mid_time(self):
        s = _transit_shape()
        dt = np.linspace(0, 0.04, 50)
        left = ct.event_model(5.0 - dt, 60 / DAY, 5.0, s.period, s.a_over_r,
                              s.k, s.b, s.u1, s.u2)
        right = ct.event_model(5.0 + dt, 60 / DAY, 5.0, s.period,
                               s.a_over_r, s.k, s.b, s.u1, s.u2)
        assert np.max(np.abs(left - right)) < 1e-9

    def test_exposure_smearing_softens_a_sharp_ingress(self):
        # White-dwarf eclipse: 120 s exposures must round the ingress
        # that 1 s exposures resolve.
        args = (0.0, 0.3443, 107.0, 9.1, 0.0, 0.3, 0.2, 0.8)
        t = np.array([-0.0050])                      # mid-ingress region
        t = np.linspace(-0.0060, -0.0040, 41)
        sharp = ct.event_model(t, 1.0 / DAY, *args)
        smooth = ct.event_model(t, 120.0 / DAY, *args)
        assert np.max(np.abs(np.diff(smooth))) < np.max(np.abs(np.diff(sharp)))


# ---------------------------------------------------------------------------
# Ephemerides and eras
# ---------------------------------------------------------------------------
class TestEphemeris:
    def test_prediction_and_error_propagation(self):
        t, sig = ct.predict_event(2455000.0, 1e-4, 2.0, 1e-6, 1000)
        assert t == pytest.approx(2457000.0)
        assert sig == pytest.approx(math.hypot(1e-4, 1e-3))

    def test_quadratic_term(self):
        t, sig = ct.predict_event(0.0, 0.0, 1.0, 0.0, 1000, quad=1e-9,
                                  sig_quad=1e-10)
        assert t == pytest.approx(1000.001)
        assert sig == pytest.approx(1e-4)

    def test_events_in_window_and_nearest(self):
        assert ct.events_in_window(10.2, 13.1, 0.0, 1.0) == [11, 12, 13]
        assert ct.events_in_window(10.2, 10.8, 0.0, 1.0) == []
        assert ct.nearest_epoch(10.8, 0.0, 1.0) == 11

    def test_clock_era_map_and_december_2024_overlap(self):
        assert ct.clock_era("2023-03-08") == "A"
        assert ct.clock_era("2024-03-04", 7) == "B"
        assert ct.clock_era("2024-05-29", 72) == "C"
        # Both cameras took frames in December 2024: the era id decides.
        assert ct.clock_era("2024-12-14", 72) == "C"
        assert ct.clock_era("2024-12-14", 76) == "D"
        assert ct.clock_era("2025-05-23", 76) == "D"
        assert ct.clock_era("2026-01-16", 76) == "E"
        assert ct.clock_era("2026-04-08", 78) == "F"
        assert ct.clock_era("2026-06-29", 81) == "G"
        assert ct.clock_era("2019-01-01") == "?"
        assert "iKon" in ct.clock_era_label("C")

    def test_eb_ephemerides_agree_with_their_alternates(self):
        # A typing error in a literature constant would show up as two
        # published ephemerides of one star disagreeing by minutes in
        # 2026.  They must agree to well under the acceptance criterion.
        for name, t in b.EB_TARGETS.items():
            p, a = t["primary"], t["alt"]
            e = ct.nearest_epoch(2461100.0, p["t0"], p["period"])
            tp, _ = ct.predict_event(p["t0"], p["sig_t0"], p["period"],
                                     p["sig_p"], e)
            ea = ct.nearest_epoch(tp, a["t0"], a["period"])
            ta, _ = ct.predict_event(a["t0"], a["sig_t0"], a["period"],
                                     a["sig_p"], ea, a.get("quad", 0.0))
            assert abs(ta - tp) * DAY < ct.ACCEPT_OC_S / 2, name


# ---------------------------------------------------------------------------
# Light-curve assembly
# ---------------------------------------------------------------------------
class TestAssembly:
    def test_native_peak_ratio(self):
        assert ct.native_peak_ratio(3.0, 1) == 1.0
        r_sharp = ct.native_peak_ratio(1.5, 2)
        r_soft = ct.native_peak_ratio(4.0, 2)
        assert r_sharp > r_soft > 1.0
        # Hand check: FWHM 3.4 binned px -> sigma_native 2.888 px; the
        # four native pixels are at offsets (0,0),(1,0),(0,1),(1,1).
        s2 = 2.0 * (3.4 * 2 / 2.3548) ** 2
        want = 4.0 / (1 + 2 * math.exp(-1 / s2) + math.exp(-2 / s2))
        assert ct.native_peak_ratio(3.4, 2) == pytest.approx(want)

    def test_successive_difference_rms_ignores_trend_and_event(self):
        rng = np.random.default_rng(4)
        y = rng.normal(0, 0.01, 4000)
        assert ct.successive_difference_rms(y) == pytest.approx(0.01,
                                                                rel=0.08)
        ramp = y + np.linspace(0, 1.0, 4000)         # a 100-sigma trend
        assert ct.successive_difference_rms(ramp) == pytest.approx(
            0.01, rel=0.08)
        assert math.isnan(ct.successive_difference_rms([1.0, 2.0]))

    def test_comparison_pruning_drops_the_variable_star(self):
        rng = np.random.default_rng(2)
        flux = 1e5 * (1 + rng.normal(0, 0.003, (200, 5)))
        flux[:, 3] *= 1 + rng.normal(0, 0.05, 200)   # one noisy star
        keep = ct.select_comparisons(flux, np.ones((200, 5), bool))
        assert list(keep) == [True, True, True, False, True]

    def test_comparison_needs_frame_coverage_and_a_minimum(self):
        flux = np.ones((100, 3)) * 1e4
        usable = np.ones((100, 3), bool)
        usable[:50, 0] = False                       # saturated half the run
        keep = ct.select_comparisons(flux, usable)
        assert not keep[0] and keep[1] and keep[2]
        usable[:50, 1] = False                       # now only one remains
        assert not ct.select_comparisons(flux, usable).any()

    def test_cloud_frames_are_dropped(self):
        tf = np.full(10, 100.0)
        cf = np.full((10, 2), 1000.0)
        tf[4], cf[4] = 20.0, 200.0                   # an 80 % cloud
        lc, use = ct.differential_lightcurve(
            tf, np.ones(10, bool), cf, np.ones((10, 2), bool),
            np.array([True, True]))
        assert not use[4] and use.sum() == 9
        assert np.nanmedian(lc[use]) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# The mid-time fit
# ---------------------------------------------------------------------------
class TestFit:
    def test_recovers_mid_time_without_bias(self):
        shape = _transit_shape()
        true = 2460000.003
        errs = []
        for seed in range(12):
            t, f, e, ex = _synthetic(shape, true, seed=seed)
            r = ct.fit_event(t, f, e, ex, shape)
            assert r["status"] == ct.STATUS_OK
            errs.append((r["t0"] - true) * DAY)
        errs = np.array(errs)
        sem = errs.std(ddof=1) / math.sqrt(len(errs))
        assert abs(errs.mean()) < 3 * sem + 2.0      # signed, not |.|
        assert errs.std() < 45.0

    def test_fit_is_blind_to_the_prediction(self):
        # The event sits 1,065 s from where an ephemeris would put it;
        # nothing tells the fit, and it must still land on it.
        shape = _transit_shape()
        true = 2460000.003 + ct.CV_OFFSET_S / DAY
        t, f, e, ex = _synthetic(shape, true, seed=7)
        r = ct.fit_event(t, f, e, ex, shape)
        assert r["status"] == ct.STATUS_OK
        assert abs(r["t0"] - true) * DAY < 90.0

    def test_formal_error_is_in_days_not_scaled_units(self):
        shape = _transit_shape()
        t, f, e, ex = _synthetic(shape, 2460000.003)
        r = ct.fit_event(t, f, e, ex, shape)
        assert 3.0 < r["t0_err_formal"] * DAY < 120.0

    def test_one_sided_event_is_refused(self):
        shape = _transit_shape()
        # Data start in mid-transit: only the egress is seen.
        t, f, e, ex = _synthetic(shape, 2460000.0 - 0.08, span=(-0.08, 0.09))
        r = ct.fit_event(t, f, e, ex, shape)
        assert r["status"] in (ct.STATUS_ONE_SIDED, ct.STATUS_NO_DIP)

    def test_too_few_points(self):
        shape = _transit_shape()
        t, f, e, ex = _synthetic(shape, 2460000.0, n=8)
        assert ct.fit_event(t, f, e, ex, shape)["status"] == \
            ct.STATUS_TOO_FEW

    def test_deep_binary_eclipse_with_free_duration_and_depth(self):
        shape = ct.EventShape(0.34433, 107.0, 9.1, 0.0, 0.3, 0.2, f1=0.8,
                              free=("a_over_r", "f1"), baseline_order=2)
        rng = np.random.default_rng(3)
        t = np.arange(-0.03, 0.03, 130 / DAY) + 2461139.9
        ex = 120 / DAY
        true = 2461139.9004
        f = ct.event_model(t, ex, true, shape.period, 100.0, 9.1, 0.0, 0.3,
                           0.2, 0.7) + rng.normal(0, 0.04, len(t))
        r = ct.fit_event(t, f, np.full(len(t), 0.04), ex, shape)
        assert r["status"] == ct.STATUS_OK
        assert abs(r["t0"] - true) * DAY < 30.0
        assert r["f1"] == pytest.approx(0.7, abs=0.05)
        assert r["depth"] == pytest.approx(0.7, abs=0.05)

    def test_outlier_clip(self):
        r = np.zeros(100)
        r[::2] = 0.01
        r[1::2] = -0.01
        r[10] = 1.0
        keep = ct.clip_outliers(r)
        assert not keep[10] and keep.sum() == 99
        assert ct.clip_outliers(np.zeros(5)).all()   # zero MAD: keep all


class TestBaselineRegressors:
    def test_jump_steps_finds_a_repointing_and_ignores_drift(self):
        n = 100
        x = np.linspace(0, 6, n)                     # slow drift: no jump
        y = np.zeros(n)
        assert ct.jump_steps(x, y).shape == (n, 0)
        x[40:] += 12.0                               # one re-centring
        st = ct.jump_steps(x, y)
        assert st.shape == (n, 1)
        assert st[:40, 0].sum() == 0 and st[40:, 0].all()

    def test_jump_steps_refuses_short_segments_and_caps_the_count(self):
        n = 100
        x = np.zeros(n)
        x[3:] += 20.0                                # 3-point segment
        assert ct.jump_steps(x, np.zeros(n)).shape[1] == 0
        x = np.zeros(200)
        for i in range(10, 200, 10):                 # 19 jumps, growing
            x[i:] += 6.0 + i / 100.0
        st = ct.jump_steps(x, np.zeros(200))
        assert st.shape[1] == ct.JUMP_MAX_STEPS

    def test_step_in_the_light_curve_is_absorbed_not_timed(self):
        # A 1 % step in mid-transit (the star re-centred onto different
        # pixels).  Without the offset the mid-time is dragged; with it,
        # the bias goes away.  Averaged over seeds, signed.
        shape = _transit_shape()
        true = 2460000.003
        plain, stepped = [], []
        for seed in range(8):
            t, f, e, ex = _synthetic(shape, true, seed=40 + seed,
                                     noise=0.002)
            step = (t > true + 0.004).astype(float)
            f = f * (1.0 - 0.010 * step)
            plain.append((ct.fit_event(t, f, e, ex, shape)["t0"] - true)
                         * DAY)
            stepped.append((ct.fit_event(t, f, e, ex, shape,
                                         regressors=step[:, None])["t0"]
                            - true) * DAY)
        assert abs(np.mean(plain)) > 60.0            # the defect is real
        assert abs(np.mean(stepped)) < 30.0          # and the cure works
        assert abs(np.mean(stepped)) < abs(np.mean(plain)) / 3

    def test_regressor_fit_counts_its_parameters(self):
        shape = _transit_shape()
        t, f, e, ex = _synthetic(shape, 2460000.003)
        r0 = ct.fit_event(t, f, e, ex, shape)
        r1 = ct.fit_event(t, f, e, ex, shape,
                          regressors=np.sin(np.arange(len(t)))[:, None])
        assert r1["dof"] == r0["dof"] - 1


class TestErrors:
    def test_block_resample_keeps_length_and_runs(self):
        rng = np.random.default_rng(0)
        res = np.arange(100, dtype=float)
        out = ct.resample_blocks(res, 10, rng)
        assert len(out) == 100
        # Inside a block consecutive values still differ by exactly 1.
        steps = np.diff(out)
        assert (steps == 1.0).sum() >= 85

    def test_block_length_in_points(self):
        t = np.arange(100) * (30.0 / DAY)            # 30 s cadence
        assert ct.block_length_points(t, 10.0) == 20

    def test_bootstrap_error_matches_repeated_experiments(self):
        shape = _transit_shape()
        true = 2460000.003
        errs = []
        for seed in range(16):
            t, f, e, ex = _synthetic(shape, true, seed=100 + seed)
            errs.append((ct.fit_event(t, f, e, ex, shape)["t0"] - true)
                        * DAY)
        empirical = float(np.std(errs))
        t, f, e, ex = _synthetic(shape, true, seed=100)
        fit = ct.fit_event(t, f, e, ex, shape)
        boot = ct.bootstrap_t0(t, f, e, ex, shape, fit, n_boot=120)
        assert 0.5 * empirical < boot["sigma"] * DAY < 2.0 * empirical
        bead = ct.prayer_bead_t0(t, f, e, ex, shape, fit)
        assert 0.3 * empirical < bead["sigma"] * DAY < 3.0 * empirical

    def test_injection_recovers_cv_sized_offset_signed(self):
        shape = _transit_shape()
        t, f, e, ex = _synthetic(shape, 2460000.003, seed=5)
        fit = ct.fit_event(t, f, e, ex, shape)
        for shift in (ct.CV_OFFSET_S, -ct.CV_OFFSET_S):
            inj = ct.inject_and_recover(t, f, e, ex, shape, fit, shift, 4, 1)
            assert inj["n_ok"] == 4
            assert abs(inj["bias_s"]) < 60.0

    def test_injection_out_of_window_is_reported_unrecoverable(self):
        shape = _transit_shape()
        t, f, e, ex = _synthetic(shape, 2460000.003, seed=5)
        fit = ct.fit_event(t, f, e, ex, shape)
        inj = ct.inject_and_recover(t, f, e, ex, shape, fit, 9000.0, 3, 1)
        assert inj["n_ok"] == 0 and math.isnan(inj["bias_s"])


# ---------------------------------------------------------------------------
# Combination and verdict
# ---------------------------------------------------------------------------
class TestCombine:
    def test_weighted_mean_and_both_errors(self):
        c = ct.combine_oc([10.0, 20.0], [10.0, 10.0])
        assert c["wmean"] == pytest.approx(15.0)
        assert c["wmean_err"] == pytest.approx(10 / math.sqrt(2))
        assert c["chi2"] == pytest.approx(0.5) and c["dof"] == 1
        assert c["scatter_err"] == pytest.approx(5.0)
        assert c["sigma_adopted"] == pytest.approx(10 / math.sqrt(2))

    def test_overdispersed_set_adopts_scatter_not_a_rescale(self):
        # Error bars of 1 s on values 100 s apart: the formal error
        # (0.7 s) must NOT be what is quoted.
        c = ct.combine_oc([0.0, 100.0], [1.0, 1.0])
        assert c["sigma_adopted"] == pytest.approx(50.0)
        assert c["chi2"] / c["dof"] > 1000           # reported, not hidden

    def test_single_value_has_no_scatter_error(self):
        c = ct.combine_oc([5.0], [3.0])
        assert c["scatter_err"] is None and c["sigma_adopted"] == 3.0
        assert ct.combine_oc([], [])["n"] == 0

    def test_verdicts(self):
        assert ct.era_verdict(10.0, 20.0) == "PASS"
        assert ct.era_verdict(1065.0, 30.0) == "OFFSET MEASURED"
        assert ct.era_verdict(90.0, 40.0).startswith("PASS (central")
        assert ct.era_verdict(150.0, 100.0) == "UNRESOLVED"


# ---------------------------------------------------------------------------
# Build-script logic
# ---------------------------------------------------------------------------
class TestScriptHelpers:
    def test_series_ids_distinguish_two_requests_of_one_star(self):
        a = b.series_id_of("WASP-52b", "2023-11-20", "R", "High Gain", 7,
                           "rawimage", "WASP 52")
        c = b.series_id_of("WASP-52b", "2023-11-20", "R", "High Gain", 7,
                           "rawimage", "WASP-52")
        assert a != c and "HighGain" in a

    def test_sexagesimal_and_separation(self):
        ra, dec = b.sexagesimal_to_deg("10:19:38.0089", "-09:48:22.603")
        assert ra == pytest.approx(154.908370, abs=1e-5)
        assert dec == pytest.approx(-9.806279, abs=1e-5)
        assert b.ang_sep_arcmin(10.0, 0.0, 10.0, 1.0) == pytest.approx(60.0)

    def test_eclipse_window_holds_a_cv_sized_displacement(self):
        t = np.linspace(-0.2, 0.2, 2001)
        shape = b.shape_of(b.EB_TARGETS["NSVS 07826147"])
        win = b.event_window("eclipse", t, 0.0, shape)
        half = t[win].max()
        # Wide enough that a +-1,065 s shifted eclipse is still inside...
        assert half > ct.CV_OFFSET_S / DAY + shape.duration_d() / 2
        # ...and narrow enough to exclude the secondary eclipse at P/2.
        assert half < shape.period / 2 - shape.duration_d() / 2
        assert b.event_window("transit", t, 0.0, shape).all()

    def test_saturated_target_frames_are_not_used(self):
        n = 40
        rng = np.random.default_rng(1)
        flux = np.ones((n, 4, 3)) * 1e5 * (1 + rng.normal(0, 0.002,
                                                          (n, 4, 3)))
        arr = {"fwhm": np.full(n, 3.4), "flux": flux,
               "err": np.full((n, 4, 3), 300.0),
               "peak": np.full((n, 4), 20000.0),
               "good": np.ones((n, 4), bool)}
        arr["peak"][:10, 0] = 59000.0     # x1.06 native ratio -> over veto
        out = b.build_lightcurves(arr, veto_adu=60200.0, binning=2)
        assert out["n_target_saturated"] == 10
        ap = out["apertures"][out["adopted"]]
        assert not ap["use"][:10].any() and ap["use"][10:].all()
        # The same peaks are fine on an unbinned sensor.
        out1 = b.build_lightcurves(arr, veto_adu=60200.0, binning=1)
        assert out1["n_target_saturated"] == 0

    def test_veto_lookup_reads_the_manifest_never_a_constant(self):
        con = sqlite3.connect(":memory:")
        con.execute("CREATE TABLE detector_params (era_group TEXT, "
                    "quantity TEXT, value REAL)")
        con.executemany("INSERT INTO detector_params VALUES (?,?,?)", [
            ("High Gain", "saturation_veto_adu", 3200.0),
            ("Mode0", "saturation_veto_adu", 60200.0),
            ("1MHz High Sensitivity 16-bit", "saturation_veto_adu", 59500.0),
            ("(blank 2026)", "saturation_veto_adu", 60100.0)])
        assert b.veto_for(con, "High Gain") == 3200.0
        assert b.veto_for(con, "5MHz High Sensitivity 16-bit") == 59500.0
        assert b.veto_for(con, "") == 60100.0
        assert b.veto_for(con, "Never Seen") == 3200.0   # most conservative


def _summary_fixture(tmp_path):
    """A products DB with two raw/reduced twin events and one lone event,
    plus the smallest manifest the summary reads."""
    db = b.open_db(tmp_path / "clock.sqlite")
    db.execute("""CREATE TABLE s3b_series (series_id TEXT, target TEXT,
        kind TEXT, manifest_target TEXT, night TEXT, filter TEXT,
        readoutm TEXT, era_id INTEGER, tree TEXT, clock_era TEXT,
        n_frames INTEGER, jd_start REAL, jd_end REAL, exptime_s REAL,
        offset_arcmin REAL, binning INTEGER, admitted INTEGER)""")
    db.execute("""CREATE TABLE s3b_census_events (event_id TEXT,
        series_id TEXT, target TEXT, epoch INTEGER, t_pred_bjd REAL,
        sig_pred_s REAL, pre_h REAL, post_h REAL, t14_min REAL,
        admitted INTEGER, reason TEXT)""")
    db.execute("""CREATE TABLE s3b_fits (event_id TEXT, status TEXT,
        t0_bjd REAL, sig_stat_s REAL, sig_model_s REAL)""")
    db.execute("""CREATE TABLE s3b_injection (event_id TEXT, shift_s REAL,
        n_ok INTEGER, n_draws INTEGER, bias_s REAL, scatter_s REAL)""")
    t_pred = 2461139.9
    rows = [  # (sid, tree, era, epoch, oc_s)
        ("s_raw", "rawimage", 78, 100, +20.0),
        ("s_red", "reduced", 79, 100, +26.0),
        ("s_lone", "rawimage", 78, 103, -10.0)]
    for sid, tree, era, epoch, oc in rows:
        db.execute("INSERT INTO s3b_series VALUES (?,?,?,?,?,?,?,?,?,?,"
                   "?,?,?,?,?,?,?)",
                   (sid, "GK Vir", "eclipse", "GK Vir",
                    "2026-04-08" if epoch == 100 else "2026-04-09", "g",
                    "Fast", era, tree, "F", 100, 0, 0, 120, 0, 2, 1))
        db.execute("INSERT INTO s3b_census_events VALUES (?,?,?,?,?,?,?,?,"
                   "?,?,?)", (sid + "#e", sid, "GK Vir", epoch,
                              t_pred + (epoch - 100) * 0.34433, 4.5, 1, 1,
                              15, 1, ""))
        db.execute("INSERT INTO s3b_fits VALUES (?,?,?,?,?)",
                   (sid + "#e", ct.STATUS_OK,
                    t_pred + (epoch - 100) * 0.34433 + oc / DAY, 12.0, 5.0))
        db.execute("INSERT INTO s3b_injection VALUES (?,?,?,?,?,?)",
                   (sid + "#e", ct.CV_OFFSET_S, 8, 8, 1.5, 9.0))
    db.commit()
    man = sqlite3.connect(tmp_path / "manifest.sqlite")
    man.execute("""CREATE TABLE s3_clock_eclipses (tag TEXT,
        n_points INTEGER, o_minus_c_s REAL, o_minus_c_err_s REAL,
        clock_bound_s REAL, status TEXT)""")
    man.execute("CREATE TABLE s3_build_meta (key TEXT, value TEXT)")
    man.commit()
    man.close()
    return db


class TestSummary:
    def test_twins_enter_the_mean_once_and_pairs_are_recorded(self, tmp_path):
        db = _summary_fixture(tmp_path)
        b.stage_summary(db, tmp_path / "manifest.sqlite",
                        tmp_path / "no_cv.sqlite", dict(b.EB_TARGETS))
        prim = dict(db.execute("SELECT series_id, is_primary FROM s3b_oc"))
        assert prim == {"s_raw": 1, "s_red": 0, "s_lone": 1}
        pair = db.execute("SELECT kind, delta_s FROM s3b_pairs").fetchall()
        assert pair == [("raw_vs_reduced", pytest.approx(6.0, abs=0.01))]
        era = db.execute("""SELECT n_targets, n_series, oc_s, verdict,
                                   inj1065_n_ok
                            FROM s3b_era WHERE clock_era = 'F'""").fetchone()
        # Mean of the two PRIMARY values (+20, -10), equal weights.
        assert era[0] == 1 and era[1] == 2
        assert era[2] == pytest.approx(5.0, abs=0.05)
        assert era[4] == 16                          # twin not counted
        none = db.execute("SELECT verdict FROM s3b_era WHERE clock_era='D'"
                          ).fetchone()[0]
        assert none == "NO CLOCK TARGET IN ERA"

    def test_standard_systematics_are_carried_not_dropped(self, tmp_path):
        db = _summary_fixture(tmp_path)
        b.stage_summary(db, tmp_path / "manifest.sqlite",
                        tmp_path / "no_cv.sqlite", dict(b.EB_TARGETS))
        sig_meas, sig_eph, sig_sys, sig_tot, alt = db.execute(
            """SELECT sig_meas_s, sig_eph_s, sig_sys_s, sig_total_s,
                      alt_minus_primary_s FROM s3b_oc
               WHERE series_id = 's_raw'""").fetchone()
        assert sig_meas == pytest.approx(13.0)       # hypot(12, 5)
        # GK Vir's systematic IS the disagreement of its two ephemerides.
        assert sig_sys == pytest.approx(abs(alt))
        assert sig_tot == pytest.approx(
            math.sqrt(sig_meas ** 2 + sig_eph ** 2 + sig_sys ** 2))
        tgt_tot = db.execute("SELECT sig_total_s FROM s3b_target_oc"
                             ).fetchone()[0]
        assert tgt_tot > sig_sys                     # never smaller


class TestGradesAndAnchoring:
    def test_low_precision_and_anchored_rows_never_enter_a_mean(
            self, tmp_path):
        db = _summary_fixture(tmp_path)
        # s_lone becomes a 500 s-error "measurement" sitting at +2,000 s;
        # s_red becomes an anchored one-sided fit.
        db.execute("UPDATE s3b_fits SET sig_stat_s = 500, "
                   "t0_bjd = t0_bjd + 2010.0 / 86400 "
                   "WHERE event_id = 's_lone#e'")
        db.execute("UPDATE s3b_fits SET status = ? WHERE event_id = "
                   "'s_red#e'", (ct.STATUS_ANCHORED,))
        db.commit()
        b.stage_summary(db, tmp_path / "manifest.sqlite",
                        tmp_path / "no_cv.sqlite", dict(b.EB_TARGETS))
        grade = dict(db.execute("SELECT series_id, grade FROM s3b_oc"))
        assert grade == {"s_raw": ct.GRADE_TIMING,
                         "s_red": ct.GRADE_ANCHORED,
                         "s_lone": ct.GRADE_LOW}
        oc, n_series, n_low, oc_all = db.execute(
            """SELECT oc_s, n_series, n_low_precision, oc_all_s
               FROM s3b_era WHERE clock_era = 'F'""").fetchone()
        assert n_series == 1 and n_low == 1
        assert oc == pytest.approx(20.0, abs=0.05)   # the outlier is out
        assert oc_all > oc                           # ...but still shown

    def test_era_with_only_low_precision_events_gets_no_verdict(
            self, tmp_path):
        db = _summary_fixture(tmp_path)
        db.execute("UPDATE s3b_fits SET sig_stat_s = 500")
        db.commit()
        b.stage_summary(db, tmp_path / "manifest.sqlite",
                        tmp_path / "no_cv.sqlite", dict(b.EB_TARGETS))
        verdict, oc = db.execute("SELECT verdict, oc_s FROM s3b_era "
                                 "WHERE clock_era = 'F'").fetchone()
        assert verdict == "NO TIMING-GRADE EVENT IN ERA" and oc is None

    def test_anchored_refit_times_an_egress_only_eclipse(self):
        tgt = b.EB_TARGETS["NSVS 07826147"]
        shape = b.shape_of(tgt)
        sd = dict(period=shape.period, a_over_r=shape.a_over_r, k=shape.k,
                  b=shape.b, u1=shape.u1, u2=shape.u2, f1=shape.f1,
                  free=shape.free, baseline_order=shape.baseline_order)
        rng = np.random.default_rng(11)
        true_ar, true_f1 = 6.2, 0.85
        jobs, results = [], []

        def make(eid, t0, lo, hi):
            t = np.arange(lo, hi, 35 / DAY) + t0
            ex = np.full(len(t), 30 / DAY)
            f = ct.event_model(t, ex, t0, shape.period, true_ar, shape.k,
                               shape.b, shape.u1, shape.u2, true_f1)
            f = f + rng.normal(0, 0.005, len(t))
            ap = {"t": t, "flux": f, "err": np.full(len(t), 0.005),
                  "ex": ex}
            return {"event_id": eid, "apertures": [ap], "adopted": 0,
                    "shape": sd, "seed": 1}, ap

        for i in range(3):                           # complete eclipses
            job, ap = make(f"NSVS 07826147|full{i}", 2461000.0 + i,
                           -0.02, 0.02)
            fit = ct.fit_event(ap["t"], ap["flux"], ap["err"], ap["ex"],
                               shape)
            assert fit["status"] == ct.STATUS_OK
            jobs.append(job)
            results.append({"event_id": job["event_id"],
                            "status": ct.STATUS_OK,
                            "a_over_r": fit["a_over_r"], "f1": fit["f1"]})
        t0 = 2461010.0                               # egress only
        job, ap = make("NSVS 07826147|egress", t0, 0.004, 0.035)
        jobs.append(job)
        results.append({"event_id": job["event_id"],
                        "status": ct.STATUS_ONE_SIDED,
                        "depth_expected": 0.6})
        out = b.anchored_refits(jobs, results)
        assert len(out) == 1
        assert out[0]["status"] == ct.STATUS_ANCHORED
        assert abs(out[0]["t0"] - t0) * DAY < 15.0
        # The refused row was replaced, not duplicated.
        assert all(r["event_id"] != job["event_id"] for r in results)

    def test_anchoring_needs_enough_complete_eclipses(self):
        tgt = b.EB_TARGETS["NSVS 07826147"]
        s = b.shape_of(tgt)
        sd = dict(period=s.period, a_over_r=s.a_over_r, k=s.k, b=s.b,
                  u1=s.u1, u2=s.u2, f1=s.f1, free=s.free,
                  baseline_order=s.baseline_order)
        jobs = [{"event_id": "NSVS 07826147|x", "shape": sd}]
        results = [{"event_id": "NSVS 07826147|x",
                    "status": ct.STATUS_ONE_SIDED}]
        assert b.anchored_refits(jobs, results) == []
        assert len(results) == 1
