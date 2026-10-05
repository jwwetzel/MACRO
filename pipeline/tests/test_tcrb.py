"""Tests for macro_tcrb — the T CrB project's decision rules.

Synthetic-recovery tests: every rule is shown to return what was put in, so
a number in the paper is never the product of an untested rule.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from macro_tcrb import detect as dt  # noqa: E402


# ---------------------------------------------------------------------------
# A5a — the pre-registered detection rule
# ---------------------------------------------------------------------------
def _season(n=40, cadence=2.0):
    t = np.arange(n) * cadence
    return t, np.full(n, 0.5)


def test_constant_series_calls_nothing():
    t, s = _season()
    e = np.zeros_like(t)
    assert dt.step_candidates(t, e, s) == []


def test_persistent_step_is_called_and_transient_is_not():
    t, s = _season()
    e = np.where(np.arange(t.size) >= 20, 5.0, 0.0)     # 10 sigma step
    ev = dt.step_candidates(t, e, s)
    assert [(x.i, x.j) for x in ev] == [(19, 20)]
    spike = np.zeros_like(t)
    spike[20] = 5.0                                      # one-night outlier
    assert dt.step_candidates(t, spike, s) == []


def test_gap_longer_than_four_days_breaks_consecutiveness():
    t = np.array([0.0, 1.0, 6.0, 7.0])
    s = np.full(4, 0.5)
    e = np.array([0.0, 0.0, 5.0, 5.0])
    assert dt.step_candidates(t, e, s) == []


def test_corroboration_by_other_grism():
    t, s = _season()
    e = np.where(np.arange(t.size) >= 20, 5.0, 0.0)
    ev = dt.step_candidates(t, e, s)[0]
    assert dt.corroborate(ev, t + 0.1, e, s) == "corroborated"
    assert dt.corroborate(ev, t + 0.1, -e, s) == "contradicted"
    assert dt.corroborate(ev, t + 500.0, e, s) == "single-grism"
    assert dt.corroborate(ev, None, None, None) == "single-grism"


def test_constancy_chi2_and_dof():
    e = np.array([1.0, -1.0, 1.0, -1.0])
    chi2, dof, p = dt.constancy(e, np.ones(4))
    assert chi2 == pytest.approx(4.0) and dof == 3
    assert 0.2 < p < 0.3


def test_false_alarm_rate_is_small_at_three_sigma():
    t, s = _season(60)
    ef = dt.expected_false_steps(t, s, 3.0, n_mc=2000)
    assert ef < dt.MAX_FALSE_PER_SEASON
    thr, ef2 = dt.calibrated_threshold(t, s, n_mc=500)
    assert thr == dt.STEP_SIGMA and ef2 <= dt.MAX_FALSE_PER_SEASON


def test_min_recoverable_step_is_a_few_sigma():
    t, s = _season(30)
    amp = dt.min_recoverable_step(t, s, 3.0, n_mc=200)
    # a 3-sigma step test on sqrt(2)*0.5 needs > 2.1 and is reached by ~3.5
    assert 2.0 < amp < 5.0


def test_flux_verdict():
    assert dt.flux_verdict(1.0, 0.1, 2.0, 0.1) == "accretion"
    assert dt.flux_verdict(1.0, 0.5, 1.2, 0.5) == "continuum-driven"
    assert dt.flux_verdict(1.0, np.nan, 2.0, 0.1) == "untested"


# ---------------------------------------------------------------------------
# P0 — mechanical epochs
# ---------------------------------------------------------------------------
from macro_tcrb import p0  # noqa: E402


def test_crosswalk_flags_hardware_crossing_only():
    cause = {"A:1": "camera", "A:2": "wheel_map", "B:1": "flipstat"}
    v = {x["grism_epoch"]: x["verdict"] for x in p0.crosswalk_verdicts(
        [("g1", "A:1"), ("g1", "A:2"), ("g2", "A:2"), ("g2", "B:1"),
         ("g3", "B:1")], cause)}
    assert v == {"g1": "merged_soft", "g2": "CROSSES_HARDWARE",
                 "g3": "one_to_one"}


def test_master_allowed_requires_equal_epoch():
    assert p0.master_allowed("ASI:2024-12-16", "ASI:2024-12-16")
    assert not p0.master_allowed("ASI:2024-12-16", "ASI:2025-10-11")
    assert not p0.master_allowed(None, None)


def test_temp_group_boundaries():
    assert p0.temp_group(-10.0) == "cold"
    assert p0.temp_group(0.3) == "warm"
    assert p0.temp_group(-6.0) == "intermediate"
    assert p0.temp_group(None) == "unknown"


def test_flank_residual_removes_linear_gradient():
    rows = np.arange(200, dtype=float)[:, None] * np.ones((1, 50))
    strip = 100 + 0.5 * rows
    mean, err = p0.flank_residual(strip, half=12, gap=18, width=30,
                                  centre=100)
    assert abs(mean) < 1e-9


def test_lambda_eff_interpolation():
    known = [(-0.2, 8000.0), (0.0, 6400.0), (0.3, 4400.0)]
    # slope decreasing with wavelength -> sort by slope ascending
    assert p0.lambda_eff_from_slope(0.0, known) == pytest.approx(6400.0)
    assert p0.lambda_eff_from_slope(0.15, known) == pytest.approx(5400.0)


def test_colour_slope_recovers_injected_term():
    rng = np.random.default_rng(3)
    c = rng.uniform(0.3, 2.5, 300)
    m_ref = rng.uniform(10, 15, 300)
    m = m_ref + 0.12 * c + 25.0 + rng.normal(0, 0.01, 300)
    m[:5] += 1.0                                     # five outliers
    sl, err, icp, n = p0.colour_slope(m, m_ref, c)
    assert sl == pytest.approx(0.12, abs=0.005) and n >= 280


def test_zmag_band_match_picks_the_right_band():
    rng = np.random.default_rng(4)
    r = rng.uniform(10, 15, 100)
    g = r + rng.uniform(0.2, 1.2, 100)
    m_inst = r - 21.0 + rng.normal(0, 0.01, 100)
    out = p0.zmag_band_match(m_inst, 21.0, {"r": r, "g": g})
    assert out["r"]["scatter"] < out["g"]["scatter"]
    assert abs(out["r"]["offset"]) < 0.01


def test_exposure_offset_fit_recovers_delta():
    t = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 6.0])
    c = 1e6 * (t - 0.05)
    r = p0.exposure_offset_fit(t, c, 0.001 * c)
    assert r["delta_s"] == pytest.approx(-0.05, abs=1e-6)


# ---------------------------------------------------------------------------
# Phase B / C1
# ---------------------------------------------------------------------------
from macro_tcrb import phot as ph  # noqa: E402


def test_pick_master_is_mode_matched():
    cands = [
        {"kind": "dark", "mode": "High Gain", "code": "", "exptime": 16,
         "night": "2023-10-16", "path": "hg16"},
        {"kind": "dark", "mode": "High Gain StackPro", "code": "",
         "exptime": 16, "night": "2023-10-16", "path": "sp16"},
        {"kind": "dark", "mode": "High Gain StackPro", "code": "",
         "exptime": 128, "night": "2024-02-01", "path": "sp128"},
        {"kind": "flat", "mode": "High Gain", "code": "B", "exptime": 0.7,
         "night": "2023-10-18", "path": "hgB"}]
    assert ph.pick_master(cands, "dark", "High Gain StackPro", 16.0,
                          "2024-02-12")["path"] == "sp16"
    assert ph.pick_master(cands, "flat", "High Gain StackPro", 0, "2024",
                          "B") is None
    assert ph.pick_master(cands, "flat", "High Gain", 0, "2024-02-12",
                          "B")["path"] == "hgB"


def test_huber_line_no_chi2_inflation():
    rng = np.random.default_rng(5)
    x = rng.uniform(0, 2, 200)
    y = 1.0 + 0.2 * (x - 1.0) + rng.normal(0, 0.05, 200)
    s = np.full(200, 0.01)                 # errors 5x too small
    f = ph.huber_line(x, y, s, c_ref=1.0)
    assert f.zp == pytest.approx(1.0, abs=0.02)
    assert f.chi2 / f.dof > 5              # reported, not hidden
    assert f.n >= 190                      # and no good point clipped


def test_peak_verdict():
    assert ph.peak_verdict(100.0, 2000.0, 3500.0, False) == "clean"
    assert ph.peak_verdict(2500.0, 2000.0, 3500.0, False) == "above_cap"
    assert ph.peak_verdict(3500.0, 2000.0, 3500.0, False) == "clipped"
    assert ph.peak_verdict(100.0, 2000.0, 3500.0, True) == "dispersed"


def test_red_noise_has_requested_rms():
    rng = np.random.default_rng(6)
    t = np.sort(rng.uniform(0, 600, 12))
    y = ph.red_noise(12, t, 0.05, 2.0, rng)
    assert np.std(y) == pytest.approx(0.05, rel=1e-6)


def test_detrend_rms_removes_line():
    t = np.linspace(0, 600, 15)
    rms, dof = ph.detrend_rms(t, 0.001 * t + 3.0)
    assert rms < 1e-10 and dof == 13
    assert ph.detrend_model(5, 2) == (0, 0)
    assert ph.detrend_model(7, 2) == (1, 0)
    assert ph.detrend_model(12, 2) == (1, 2)


def test_flicker_limit_exceeds_noise_and_is_monotone():
    rng = np.random.default_rng(8)
    t = np.arange(10) * 60.0
    noise = [rng.normal(0, 0.01, 10) for _ in range(5)]
    obs = rng.normal(0, 0.01, 10)
    r = ph.flicker_upper_limit(t, obs, noise, n_mc=100)
    assert 0.0 < r["limit"][1.0] < 0.1


# ---------------------------------------------------------------------------
# A5 — EW definitions
# ---------------------------------------------------------------------------
from macro_tcrb import spec as sp  # noqa: E402


def _gauss_spec(ew_true=20.0, fwhm=8.0, slope=0.0, dl=0.5):
    w = np.arange(6300, 6800, dl)
    cont = 100.0 + slope * (w - 6562.8)
    sig = fwhm / 2.3548
    prof = np.exp(-0.5 * ((w - 6562.8) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    return w, cont * (1 + ew_true * prof)


def test_ew_recovers_injected_value_on_sloped_continuum():
    w, f = _gauss_spec(20.0, 8.0, slope=0.05)
    r = sp.equivalent_width(w, f, np.full(w.size, 1.0), h=30)
    assert r["ew"] == pytest.approx(20.0, rel=0.01)
    assert r["ew_err"] > 0


def test_degrade_preserves_ew_for_resolved_window():
    w, f = _gauss_spec(20.0, 4.0)
    fd = sp.degrade(w, f, r_native=10000, fwhm_target_a=14.0)
    e0 = sp.equivalent_width(w, f, h=sp.half_width(14.0))["ew"]
    e1 = sp.equivalent_width(w, fd, h=sp.half_width(14.0))["ew"]
    assert e1 == pytest.approx(e0, rel=0.01)


def test_half_width_rule():
    assert sp.half_width(None) == 30.0
    assert sp.half_width(10.0) == 30.0
    assert sp.half_width(30.0) == 45.0


def test_orbital_phase_wraps():
    assert sp.orbital_phase(100.0 + 2.5 * 10.0, 100.0, 10.0) == pytest.approx(0.5)


def test_wavelength_scale_anchor_and_dispersion():
    x = np.arange(4000.0)
    lam = sp.wavelength_scale(x, [-460.0], 2000.0, 2390.0)
    assert lam[2390] == pytest.approx(6562.8)
    assert lam[2391] - lam[2390] == pytest.approx(-0.46)


def test_red_noise_beta_white_is_one():
    rng = np.random.default_rng(9)
    b = ph.red_noise_beta([rng.normal(0, 0.01, 60) for _ in range(20)])
    assert 0.9 < b < 1.2


def test_mean_with_jitter_recovers_mean_and_scatter():
    rng = np.random.default_rng(10)
    y = 5.0 + rng.normal(0, 0.05, 80)
    r = ph.mean_with_jitter(y, np.full(80, 0.01), n_steps=1500)
    assert abs(r["mu"][1] - 5.0) < 0.02
    assert 0.035 < r["jitter"][1] < 0.07


def test_xcorr_lag_and_pair_zero_point():
    x = np.arange(2000.0)
    shape = 1 - 0.5 * np.exp(-0.5 * ((x - 800) / 5) ** 2) \
        - 0.3 * np.exp(-0.5 * ((x - 1460) / 4) ** 2)
    shifted = np.interp(x - 7.3, x, shape)
    assert sp.xcorr_lag(shape, shifted) == pytest.approx(7.3, abs=0.2)
    mins = sp.absorption_minima(shape)
    # 6563 at px 800, 6867 at px 1460 => disp = 304/660 A/px
    assert sp.pair_zero_point(mins, 304 / 660, 6563, 6867) == 800
