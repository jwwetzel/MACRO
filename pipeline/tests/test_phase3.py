"""Unit tests for :mod:`macro_phot.phase3` — the CV Phase-3 arithmetic.

The tests are grouped in the order of the module, and each one asserts a
PROPERTY rather than a remembered number wherever a property is available:
a recovered period is checked against the injected one, an alias fraction
against an analytic window, a Kaplan-Meier duty cycle against a
hand-countable example.  The handful of tests that DO pin a number are the
ones where the number is the deliverable — the celerite2/dense agreement,
and the fact that a running-median detrend suppresses an injected signal
that a joint fit recovers.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import phase3 as p3            # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures: a synthetic multi-night series with a known period, built to look
# like the real ones — nine 0.35 d nights, one per day, 219 s cadence.
# ---------------------------------------------------------------------------
PERIOD_D = 0.07908912
FREQ_CD = 1.0 / PERIOD_D


def _nightly_times(n_nights: int = 9, hours: float = 8.4,
                   cadence_s: float = 219.0, gap_d: float = 1.0) -> np.ndarray:
    """Times for ``n_nights`` runs of ``hours`` each, one night apart."""
    per_night = int(hours * 3600.0 / cadence_s)
    step = cadence_s / 86400.0
    return np.concatenate([
        n * gap_d + step * np.arange(per_night) for n in range(n_nights)])


@pytest.fixture
def series():
    """(times, mags, errors, blocks) with a clean sinusoid at PERIOD_D."""
    t = _nightly_times()
    rng = np.random.default_rng(11)
    e = np.full(t.size, 0.015)
    y = 16.0 + 0.25 * np.sin(2 * np.pi * FREQ_CD * t) + rng.normal(0.0, e)
    return t, y, e, p3.block_index(t)


# ===========================================================================
# 1.  Ephemerides and cycle counts
# ===========================================================================
def test_quoted_precision_counts_decimals():
    # Eight decimals -> half of the eighth decimal place.
    assert p3.quoted_precision_sigma("0.07908912") == pytest.approx(0.5e-8)
    assert p3.quoted_precision_sigma("0.0868") == pytest.approx(0.5e-4)
    # No decimal point is not a period this archive has; say NaN, not 0.
    assert math.isnan(p3.quoted_precision_sigma("3"))
    assert math.isnan(p3.quoted_precision_sigma(""))


def test_quoted_precision_ignores_exponent():
    assert p3.quoted_precision_sigma("1.234e-2") == pytest.approx(0.5e-3)


def test_phase_and_cycle_are_consistent():
    e0 = 2459298.4236
    t = e0 + np.array([0.0, 0.5, 1.0, 2.7]) * PERIOD_D
    ph = p3.phase_of(t, PERIOD_D, e0)
    cyc = p3.cycle_number(t, PERIOD_D, e0)
    assert ph == pytest.approx([0.0, 0.5, 0.0, 0.7], abs=1e-8)
    assert list(cyc) == [0, 0, 1, 2]
    # cycle + phase reconstructs the time exactly.
    assert (e0 + (cyc + ph) * PERIOD_D) == pytest.approx(t, abs=1e-8)


def test_cycle_number_is_floor_not_round():
    e0 = 100.0
    # 0.9 of a cycle in belongs to cycle 0, not cycle 1.
    assert p3.cycle_number([e0 + 0.9 * PERIOD_D], PERIOD_D, e0)[0] == 0


def test_cycle_ambiguity_unique_when_period_is_sharp():
    e0, t_obs = 2459298.4236, 2460734.7
    amb = p3.cycle_ambiguity(t_obs, PERIOD_D, e0, 0.5e-8)
    assert amb.n_cycles == pytest.approx((t_obs - e0) / PERIOD_D)
    assert amb.unique
    assert amb.drift_cycles < 0.5


def test_cycle_ambiguity_fails_when_period_is_loose():
    e0, t_obs = 2459298.4236, 2460734.7
    # A period good only to the fourth decimal cannot count 18,000 cycles.
    amb = p3.cycle_ambiguity(t_obs, PERIOD_D, e0, 0.5e-4)
    assert not amb.unique
    assert amb.drift_cycles > 0.5


def test_sigma_period_max_is_exactly_half_a_cycle():
    e0, t_obs = 0.0, 1000.0 * PERIOD_D
    amb = p3.cycle_ambiguity(t_obs, PERIOD_D, e0, 1e-12)
    # n * sigma_max / P must be exactly 0.5 by construction.
    assert (amb.n_cycles * amb.sigma_period_max_d / PERIOD_D
            == pytest.approx(0.5))


def test_oc_is_zero_on_a_perfect_ephemeris():
    e0 = 2459298.4236
    cyc = np.arange(5)
    t = e0 + cyc * PERIOD_D
    assert p3.oc_seconds(t, cyc, PERIOD_D, e0) == pytest.approx(0.0, abs=1e-6)


def test_linear_ephemeris_recovers_injected_period():
    e0 = 2459298.4236
    cyc = np.arange(0, 60, 3)
    rng = np.random.default_rng(3)
    sig_d = 30.0 / 86400.0
    t = e0 + cyc * PERIOD_D + rng.normal(0.0, sig_d, cyc.size)
    E, sE, P, sP, chi2nu = p3.fit_linear_ephemeris(cyc, t, np.full(cyc.size, sig_d))
    assert P == pytest.approx(PERIOD_D, abs=5 * sP)
    assert E == pytest.approx(e0, abs=5 * sE)
    assert 0.3 < chi2nu < 3.0


def test_linear_ephemeris_inflates_error_when_scatter_exceeds_errors():
    cyc = np.arange(20)
    rng = np.random.default_rng(4)
    claimed = 1.0 / 86400.0
    truth = 20.0 / 86400.0                       # 20x the claimed error
    t = cyc * PERIOD_D + rng.normal(0.0, truth, cyc.size)
    _E, _sE, _P, sP, chi2nu = p3.fit_linear_ephemeris(
        cyc, t, np.full(cyc.size, claimed))
    assert chi2nu > 10.0
    # The rescaling must have grown the bar by roughly sqrt(chi2nu).
    naive = claimed / math.sqrt(np.var(cyc) * cyc.size)
    assert sP > 5 * naive


def test_linear_ephemeris_refuses_degenerate_input():
    out = p3.fit_linear_ephemeris([1, 1, 1], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
    assert all(math.isnan(v) for v in out)


# ===========================================================================
# 2.  Periodograms
# ===========================================================================
def test_frequency_grid_step_is_one_over_oversample_baseline():
    g = p3.frequency_grid(10.0, 8.0, 20.0, oversample=10)
    assert g[0] == pytest.approx(8.0)
    assert (g[1] - g[0]) == pytest.approx(1.0 / 100.0)
    assert g[-1] >= 20.0


def test_frequency_grid_refuses_a_zero_baseline():
    assert p3.frequency_grid(0.0).size == 0


def test_block_index_splits_on_the_gap():
    t = np.array([0.0, 0.01, 0.02, 1.0, 1.01, 5.0])
    b = p3.block_index(t, gap_h=3.0)
    assert list(b) == [0, 0, 0, 1, 1, 2]


def test_block_index_handles_empty():
    assert p3.block_index([]).size == 0


def test_project_out_blocks_removes_per_block_offsets_exactly():
    b = np.array([0, 0, 0, 1, 1])
    w = np.ones(5)
    x = np.array([[1.0, 2.0, 3.0, 10.0, 12.0]])
    out = p3._project_out_blocks(x, b, w)
    assert out[0][:3] == pytest.approx([-1.0, 0.0, 1.0])
    assert out[0][3:] == pytest.approx([-1.0, 1.0])


def test_gls_recovers_the_injected_frequency(series):
    t, y, e, b = series
    grid = p3.frequency_grid(t.max() - t.min(), 8.0, 20.0)
    power = p3.gls_block_power(t, y, e, grid, b)
    f_hat, _ = p3.refine_peak(grid, power, FREQ_CD, 0.05)
    assert 1.0 / f_hat == pytest.approx(PERIOD_D, rel=1e-5)


def test_gls_power_is_immune_to_per_night_offsets(series):
    """The whole justification of the nuisance model, as an assertion."""
    t, y, e, b = series
    grid = p3.frequency_grid(t.max() - t.min(), 12.0, 13.5)
    clean = p3.gls_block_power(t, y, e, grid, b)
    rng = np.random.default_rng(7)
    shifted = y + rng.normal(0.0, 0.4, b.max() + 1)[b]
    dirty = p3.gls_block_power(t, shifted, e, grid, b)
    # Adding an arbitrary constant to each night must not move the power.
    assert dirty == pytest.approx(clean, abs=1e-10)


def test_gls_power_is_bounded_and_normalised(series):
    t, y, e, b = series
    grid = p3.frequency_grid(t.max() - t.min(), 8.0, 20.0)
    power = p3.gls_block_power(t, y, e, grid, b)
    assert np.nanmin(power) >= 0.0
    assert np.nanmax(power) <= 1.0


def test_gls_returns_nan_for_too_few_points():
    out = p3.gls_block_power([1.0, 2.0], [1.0, 2.0], [1.0, 1.0],
                             np.array([1.0, 2.0]))
    assert np.isnan(out).all()


def test_gls_amplitude_recovers_the_injected_amplitude(series):
    t, y, e, b = series
    amp = p3.gls_amplitude(t, y, e, FREQ_CD, b)
    assert amp == pytest.approx(0.25, rel=0.05)


def test_pdm_dips_at_the_true_period(series):
    t, y, _e, b = series
    grid = np.linspace(FREQ_CD - 0.4, FREQ_CD + 0.4, 401)
    theta = p3.pdm_theta(t, y, grid, blocks=b)
    f_hat, th = p3.refine_trough(grid, theta, FREQ_CD, 0.05)
    assert 1.0 / f_hat == pytest.approx(PERIOD_D, rel=1e-3)
    assert th < 0.5                       # a real dip, not noise
    # Far from the period, theta should sit near 1.
    assert theta[0] > 0.8


def test_pdm_refuses_a_series_shorter_than_its_bins():
    out = p3.pdm_theta([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], np.array([1.0]))
    assert np.isnan(out).all()


def test_spectral_window_is_one_at_zero_frequency(series):
    t = series[0]
    assert p3.spectral_window_power(t, [0.0])[0] == pytest.approx(1.0)


def test_spectral_window_of_daily_sampling_peaks_at_one_per_day():
    """One point per night for 20 nights: W(1 c/d) must be ~1."""
    t = np.arange(20.0)
    w = p3.spectral_window_power(t, [1.0, 0.5])
    assert w[0] == pytest.approx(1.0, abs=1e-9)
    assert w[1] < 0.05


def test_alias_family_folds_negatives_and_drops_zero():
    fam = dict(p3.alias_family(2.0, 1.0, orders=(-3, -2, -1, 1)))
    assert fam[-1] == pytest.approx(1.0)
    assert fam[-3] == pytest.approx(1.0)      # |2 - 3| = 1
    assert -2 not in fam                      # 2 - 2 = 0, dropped


def test_alias_fractions_are_high_for_nightly_sampling(series):
    t = series[0]
    fr = p3.alias_window_fractions(t)
    # Nine nights of 8.4 h: the +/-1 c/d sidelobes are the whole problem.
    assert fr[1] == pytest.approx(fr[-1], abs=1e-9)
    assert fr[1] > 0.3


def test_refine_peak_interpolates_off_grid():
    f = np.linspace(10.0, 11.0, 101)
    true_f = 10.5037
    power = np.exp(-((f - true_f) / 0.02) ** 2)
    f_hat, _ = p3.refine_peak(f, power, 10.5, 0.2)
    on_grid = f[int(np.argmax(power))]
    assert abs(f_hat - true_f) < abs(on_grid - true_f)


def test_refine_peak_returns_nan_outside_the_window():
    f = np.linspace(10.0, 11.0, 11)
    assert math.isnan(p3.refine_peak(f, np.ones(11), 50.0, 0.1)[0])


def test_peak_halfwidth_matches_a_known_gaussian():
    f = np.linspace(10.0, 11.0, 2001)
    sigma = 0.02
    power = np.exp(-((f - 10.5) / sigma) ** 2 / 2)
    hw = p3.peak_halfwidth(f, power, 10.5, frac=0.5)
    assert hw == pytest.approx(sigma * math.sqrt(2 * math.log(2)), rel=0.05)


def test_mhb_sigma_scales_as_expected():
    base = p3.mhb_period_sigma(400, 40.0, 0.2, 0.02, PERIOD_D)
    # Four times the points halves the error.
    assert p3.mhb_period_sigma(1600, 40.0, 0.2, 0.02, PERIOD_D) == pytest.approx(base / 2, rel=1e-9)
    # Twice the baseline halves it too.
    assert p3.mhb_period_sigma(400, 80.0, 0.2, 0.02, PERIOD_D) == pytest.approx(base / 2, rel=1e-9)
    # Twice the amplitude halves it.
    assert p3.mhb_period_sigma(400, 40.0, 0.4, 0.02, PERIOD_D) == pytest.approx(base / 2, rel=1e-9)


def test_mhb_sigma_refuses_zero_amplitude():
    assert math.isnan(p3.mhb_period_sigma(400, 40.0, 0.0, 0.02, PERIOD_D))


def test_night_bootstrap_brackets_the_truth(series):
    t, y, e, b = series
    sig, med, n = p3.night_bootstrap_period(t, y, e, b, FREQ_CD, 0.03,
                                            n_boot=40)
    assert n >= 30
    assert sig > 0
    assert med == pytest.approx(PERIOD_D, abs=5 * sig)


def test_night_bootstrap_refuses_fewer_than_three_nights():
    t = np.linspace(0.0, 0.3, 40)
    out = p3.night_bootstrap_period(t, np.zeros(40), np.ones(40),
                                    np.zeros(40, dtype=int), 12.6, 0.1)
    assert math.isnan(out[0])


def test_family_choice_prior_when_sidelobes_are_tall():
    code, sentence = p3.classify_family_choice(12, 0.92, 0.03)
    assert code == "PRIOR"
    assert "published ephemeris" in sentence


def test_family_choice_single_night():
    code, _s = p3.classify_family_choice(1, 0.97, 3.0)
    assert code == "SINGLE-NIGHT"


def test_family_choice_data_only_below_the_constant():
    code, _s = p3.classify_family_choice(12, p3.ALIAS_DECIDABLE_MAX / 2, 0.03)
    assert code == "DATA"


def test_family_choice_unresolved_without_a_prior():
    code, _s = p3.classify_family_choice(12, 0.9, 0.03, prior_available=False)
    assert code == "UNRESOLVED"


def test_agreement_passes_and_fails_at_the_bar():
    dev, ok, _s = p3.agreement(0.0790900, 1e-6, 0.07908912, 0.0)
    assert ok and abs(dev) < p3.AGREE_SIGMA
    dev, ok, _s = p3.agreement(0.0791900, 1e-6, 0.07908912, 0.0)
    assert not ok and abs(dev) > p3.AGREE_SIGMA


def test_agreement_refuses_when_neither_value_has_an_error():
    dev, ok, sentence = p3.agreement(0.079, float("nan"), 0.0791, float("nan"))
    assert not ok and math.isnan(dev)
    assert "no agreement test is possible" in sentence


# ===========================================================================
# 3.  Edge timing
# ===========================================================================
def test_ramp_is_clipped_and_centred():
    t = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    r = p3.ramp(t, 0.0, 1.0)
    assert r == pytest.approx([0.0, 0.0, 0.5, 1.0, 1.0])


def test_fit_edge_recovers_a_clean_injected_edge():
    t = np.arange(0.0, 0.02, 219.0 / 86400.0)
    t_true = 0.010
    width = 120.0 / 86400.0
    y = 15.0 + 1.5 * p3.ramp(t, t_true, width)
    e = np.full(t.size, 0.01)
    fit = p3.fit_edge(t, y, e, p3.edge_time_grid(t_true, 0.004, 401),
                      np.array([width]), 219.0)
    assert fit.accepted
    assert fit.t_edge_d == pytest.approx(t_true, abs=30.0 / 86400.0)
    assert fit.depth_mag == pytest.approx(1.5, abs=0.05)


def test_fit_edge_rejects_a_step_too_small_to_see():
    rng = np.random.default_rng(5)
    t = np.arange(0.0, 0.02, 219.0 / 86400.0)
    y = 15.0 + 0.005 * p3.ramp(t, 0.01, 120 / 86400.0) + rng.normal(0, 0.05, t.size)
    fit = p3.fit_edge(t, y, np.full(t.size, 0.05),
                      p3.edge_time_grid(0.01, 0.004, 201),
                      np.array([120 / 86400.0]), 219.0)
    assert not fit.accepted
    assert "SNR" in fit.reason


def test_fit_edge_rejects_an_edge_that_fell_in_a_gap():
    # Points either side of a 40-minute hole, edge in the middle of it.
    t = np.concatenate([np.arange(0.0, 0.006, 219 / 86400.0),
                        np.arange(0.034, 0.040, 219 / 86400.0)])
    y = 15.0 + 1.5 * p3.ramp(t, 0.020, 120 / 86400.0)
    fit = p3.fit_edge(t, y, np.full(t.size, 0.01),
                      p3.edge_time_grid(0.020, 0.008, 201),
                      np.array([120 / 86400.0]), 219.0)
    assert not fit.accepted
    assert "gap" in fit.reason


def test_fit_edge_rejects_too_few_points():
    fit = p3.fit_edge([0.0, 1.0], [1.0, 2.0], [0.1, 0.1],
                      np.array([0.5]), np.array([0.1]), 219.0)
    assert not fit.accepted
    assert "usable points" in fit.reason


def test_fit_edge_error_bar_grows_when_the_model_misfits():
    """chi2_nu rescaling: the same edge with real scatter must not report the
    same error bar as a noiseless one."""
    t = np.arange(0.0, 0.02, 219.0 / 86400.0)
    width = 120.0 / 86400.0
    clean = 15.0 + 1.5 * p3.ramp(t, 0.010, width)
    rng = np.random.default_rng(9)
    messy = clean + rng.normal(0.0, 0.08, t.size)       # 8x the claimed error
    grid = p3.edge_time_grid(0.010, 0.004, 401)
    e = np.full(t.size, 0.01)
    a = p3.fit_edge(t, clean, e, grid, np.array([width]), 219.0)
    b = p3.fit_edge(t, messy, e, grid, np.array([width]), 219.0)
    assert b.chi2nu > a.chi2nu
    assert b.sigma_t_s > a.sigma_t_s


def test_fit_edge_fixed_depth_uses_one_free_level():
    t = np.arange(0.0, 0.02, 219.0 / 86400.0)
    width = 120.0 / 86400.0
    y = 15.0 + 1.5 * p3.ramp(t, 0.010, width)
    fit = p3.fit_edge(t, y, np.full(t.size, 0.01),
                      p3.edge_time_grid(0.010, 0.004, 401),
                      np.array([width]), 219.0, fixed_depth=1.5)
    assert fit.accepted
    assert fit.depth_mag == pytest.approx(1.5)
    assert fit.t_edge_d == pytest.approx(0.010, abs=30.0 / 86400.0)


def test_band_difference_is_zero_for_identical_edges():
    t = np.array([1.0, 2.0, 3.0])
    d, s, _c = p3.band_difference(t, [10.0] * 3, t, [10.0] * 3)
    assert d == pytest.approx(0.0)
    assert s > 0


def test_band_difference_recovers_a_known_offset():
    t = np.array([1.0, 2.0, 3.0, 4.0])
    offset_s = 45.0
    d, s, _c = p3.band_difference(t + offset_s / 86400.0, [10.0] * 4,
                                  t, [10.0] * 4)
    assert d == pytest.approx(offset_s, abs=1e-6)
    assert s == pytest.approx(math.hypot(10.0, 10.0) / 2.0, rel=0.01)


def test_band_difference_inflates_when_the_offset_is_not_constant():
    t = np.array([1.0, 2.0, 3.0, 4.0])
    wobble = np.array([0.0, 200.0, -200.0, 0.0]) / 86400.0
    d, s, chi2nu = p3.band_difference(t + wobble, [10.0] * 4, t, [10.0] * 4)
    assert chi2nu > 10.0
    assert s > math.hypot(10.0, 10.0) / 2.0


def test_band_difference_returns_nan_with_no_pairs():
    d, s, c = p3.band_difference([], [], [], [])
    assert math.isnan(d) and math.isnan(s) and math.isnan(c)


# ===========================================================================
# 4.  The sigma_t injection test
# ===========================================================================
def test_bright_phase_template_is_brighter_inside_the_bright_phase():
    p = 0.08
    t = np.linspace(0.0, p, 400)
    m = p3.bright_phase_template(t, p, p / 2, 1.5, 0.33, 60 / 86400.0)
    assert m.min() == pytest.approx(-1.5, abs=1e-6)     # brightest at centre
    assert m[0] == pytest.approx(0.0, abs=1e-6)          # faint at the edges
    assert (m <= 1e-9).all()


def test_bright_phase_template_repeats_every_cycle():
    p = 0.08
    t = np.linspace(0.0, p, 97)
    a = p3.bright_phase_template(t, p, p / 2, 1.5, 0.33, 60 / 86400.0)
    b = p3.bright_phase_template(t + 5 * p, p, p / 2, 1.5, 0.33, 60 / 86400.0)
    assert a == pytest.approx(b, abs=1e-9)


def test_sigma_t_injection_is_unbiased_with_the_right_template():
    t = np.arange(0.0, 0.39, 219.0 / 86400.0)
    out = p3.sigma_t_injection(t, np.full(t.size, 0.016), PERIOD_D,
                               depth_mag=1.6, edge_width_d=120 / 86400.0,
                               bright_width_phase=0.33, n_real=60)
    assert out["n_ok"] > 30
    assert abs(out["bias_s"]) < 15.0
    assert out["sigma_t_s"] > 0.0
    assert out["total_error_s"] >= out["sigma_t_s"]


def test_sigma_t_injection_degrades_with_a_wrong_template():
    """The measurement the whole task turns on: a wrong shape must cost."""
    t = np.arange(0.0, 0.39, 219.0 / 86400.0)
    kw = dict(period_d=PERIOD_D, depth_mag=1.6,
              edge_width_d=120 / 86400.0, bright_width_phase=0.33,
              n_real=60)
    good = p3.sigma_t_injection(t, np.full(t.size, 0.016), **kw)
    bad = p3.sigma_t_injection(t, np.full(t.size, 0.016), shape_error=5.0,
                               depth_error=0.2, **kw)
    assert bad["total_error_s"] > good["total_error_s"]


def test_sigma_t_injection_returns_nan_for_a_stub_series():
    out = p3.sigma_t_injection([0.0, 1.0], [0.01, 0.01], PERIOD_D,
                               1.0, 1e-3, 0.3, n_real=3)
    assert math.isnan(out["sigma_t_s"])


def test_sigma_t_injection_is_reproducible():
    t = np.arange(0.0, 0.39, 219.0 / 86400.0)
    kw = dict(period_d=PERIOD_D, depth_mag=1.6, edge_width_d=120 / 86400.0,
              bright_width_phase=0.33, n_real=25, seed=1234)
    a = p3.sigma_t_injection(t, np.full(t.size, 0.016), **kw)
    b = p3.sigma_t_injection(t, np.full(t.size, 0.016), **kw)
    assert a["sigma_t_s"] == b["sigma_t_s"]


def test_contour_verdict_three_branches():
    ok = [{"shape_error": 1.0, "depth_error": 0.0, "total_error_s": 10.0},
          {"shape_error": 5.0, "depth_error": 0.2, "total_error_s": 40.0}]
    assert p3.contour_verdict(ok)[0] == "PUBLISHABLE"
    cond = [{"shape_error": 1.0, "depth_error": 0.0, "total_error_s": 20.0},
            {"shape_error": 5.0, "depth_error": 0.2, "total_error_s": 140.0}]
    assert p3.contour_verdict(cond)[0] == "CONDITIONAL"
    bad = [{"shape_error": 1.0, "depth_error": 0.0, "total_error_s": 90.0},
           {"shape_error": 5.0, "depth_error": 0.2, "total_error_s": 300.0}]
    assert p3.contour_verdict(bad)[0] == "NOT PUBLISHABLE"


def test_contour_verdict_with_no_finite_cells():
    assert p3.contour_verdict([{"total_error_s": float("nan")}])[0] == "NO RESULT"


# ===========================================================================
# 5.  Accretion states
# ===========================================================================
def test_phase_coverage_full_and_partial():
    p, e0 = 0.08, 0.0
    full = np.linspace(0.0, p, 200)
    assert p3.phase_coverage(full, p, e0) == pytest.approx(1.0)
    quarter = np.linspace(0.0, p / 4, 50)
    assert p3.phase_coverage(quarter, p, e0) == pytest.approx(0.25, abs=0.06)


def test_phase_coverage_of_nothing_is_zero():
    assert p3.phase_coverage([], 0.08, 0.0) == 0.0


def test_otsu_splits_two_clean_populations():
    rng = np.random.default_rng(2)
    v = np.concatenate([rng.normal(15.0, 0.1, 60), rng.normal(18.0, 0.1, 40)])
    thr, sep = p3.otsu_threshold(v)
    assert 15.5 < thr < 17.5
    assert sep > 0.8                     # genuinely bimodal


def test_otsu_separability_is_low_for_one_population():
    rng = np.random.default_rng(6)
    thr, sep = p3.otsu_threshold(rng.normal(16.0, 0.3, 400))
    assert np.isfinite(thr)
    assert sep < 0.75                    # a cut through a unimodal blob


def test_otsu_refuses_a_constant_series():
    thr, sep = p3.otsu_threshold(np.full(20, 16.0))
    assert math.isnan(thr) and math.isnan(sep)


def test_bootstrap_threshold_is_small_for_a_clean_split():
    rng = np.random.default_rng(8)
    v = np.concatenate([rng.normal(15.0, 0.1, 60), rng.normal(18.0, 0.1, 40)])
    assert p3.bootstrap_threshold(v, n_boot=80) < 0.6


def test_classify_state_uses_the_uncertainty_band():
    assert p3.classify_state(15.0, 16.5, 0.2) == "HIGH"
    assert p3.classify_state(18.0, 16.5, 0.2) == "LOW"
    assert p3.classify_state(16.55, 16.5, 0.2) == "INTERMEDIATE"
    assert p3.classify_state(float("nan"), 16.5, 0.2) == "UNKNOWN"


def test_duty_cycle_counts_faint_limits_as_low():
    #  3 bright detections, 1 faint detection, 2 limits fainter than the cut.
    mags = [15.0, 15.1, 15.2, 18.0, 19.0, 19.5]
    cens = [False, False, False, False, True, True]
    out = p3.duty_cycle(mags, cens, 16.5)
    assert out["naive"] == pytest.approx(3 / 4)
    assert out["with_limits"] == pytest.approx(3 / 6)
    assert out["n_informative_limits"] == 2
    assert out["bias"] == pytest.approx(0.25)


def test_duty_cycle_will_not_use_a_limit_brighter_than_the_threshold():
    """A limit at 15.0 says 'fainter than 15.0', which does not decide a
    16.5 cut, and must not be counted as if it did."""
    out = p3.duty_cycle([15.0, 15.1, 15.0], [False, False, True], 16.5)
    assert out["n_uninformative"] == 1
    assert out["n_informative_limits"] == 0
    assert out["with_limits"] == pytest.approx(1.0)


def test_duty_cycle_with_no_detections():
    out = p3.duty_cycle([19.0, 19.5], [True, True], 16.5)
    assert math.isnan(out["naive"])
    assert out["with_limits"] == pytest.approx(0.0)


# ===========================================================================
# 6.  Detrending discipline
# ===========================================================================
def test_matern32_is_a_valid_covariance():
    t = np.linspace(0.0, 1.0, 30)
    k = p3.matern32_cov(t, 0.3, 0.2)
    assert k == pytest.approx(k.T)
    assert np.diag(k) == pytest.approx(0.09)
    assert np.linalg.eigvalsh(k).min() > -1e-10


def test_matern32_matches_celerite2_at_the_eps_we_pin():
    """celerite2's ``Matern32Term`` is a two-exponential APPROXIMATION whose
    accuracy is set by ``eps``; the module pins it so the fast path and the
    dense reference are the same kernel."""
    pytest.importorskip("celerite2")
    from celerite2 import terms
    t = np.sort(np.random.default_rng(1).uniform(0.0, 1.0, 40))
    sigma, rho = 0.3, 0.15
    tau = np.abs(t[:, None] - t[None, :])
    ours = p3.matern32_cov(t, sigma, rho)
    pinned = terms.Matern32Term(sigma=sigma, rho=rho,
                                eps=p3.CELERITE_MATERN_EPS).get_value(tau)
    assert ours == pytest.approx(pinned, rel=1e-10, abs=1e-12)


def test_celerite2_default_eps_is_a_different_kernel():
    """The reason the constant exists: celerite2's DEFAULT eps = 0.01 is
    wrong by ~2e-5 relative, and the error shrinks as eps does.  If this
    test ever starts passing at the default, the pin can be reconsidered."""
    pytest.importorskip("celerite2")
    from celerite2 import terms
    t = np.sort(np.random.default_rng(1).uniform(0.0, 1.0, 40))
    tau = np.abs(t[:, None] - t[None, :])
    ours = p3.matern32_cov(t, 0.3, 0.15)

    def rel(eps):
        k = terms.Matern32Term(sigma=0.3, rho=0.15, eps=eps).get_value(tau)
        return float(np.max(np.abs(ours - k) / np.maximum(np.abs(k), 1e-30)))

    assert rel(0.01) > 1e-6                 # the default is NOT our kernel
    assert rel(1e-4) < rel(0.01)            # and the error is eps-controlled
    assert rel(p3.CELERITE_MATERN_EPS) < 1e-10


def test_dense_log_likelihood_matches_celerite2():
    celerite2 = pytest.importorskip("celerite2")
    from celerite2 import GaussianProcess, terms
    rng = np.random.default_rng(2)
    t = np.sort(rng.uniform(0.0, 1.0, 50))
    r = rng.normal(0.0, 0.1, 50)
    e = np.full(50, 0.02)
    ours = p3.gp_log_likelihood_dense(t, r, e, 0.2, 0.1)
    gp = GaussianProcess(
        terms.Matern32Term(sigma=0.2, rho=0.1, eps=p3.CELERITE_MATERN_EPS),
        mean=0.0)
    gp.compute(t, yerr=e)
    assert ours == pytest.approx(gp.log_likelihood(r), rel=1e-9)


def test_joint_gp_fit_backends_agree():
    pytest.importorskip("celerite2")
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0.0, 0.4, 60))
    y = 16.0 + 0.1 * np.sin(2 * np.pi * FREQ_CD * t) + rng.normal(0, 0.01, 60)
    e = np.full(60, 0.01)
    fast = p3.joint_gp_fit(t, y, e, FREQ_CD, use_celerite=True)
    slow = p3.joint_gp_fit(t, y, e, FREQ_CD, use_celerite=False)
    assert fast["backend"] == "celerite2"
    assert slow["backend"] == "dense"
    assert fast["amplitude"] == pytest.approx(slow["amplitude"], rel=1e-6)
    assert fast["loglike"] == pytest.approx(slow["loglike"], rel=1e-8)


def test_joint_gp_fit_recovers_a_signal_under_a_strong_trend():
    rng = np.random.default_rng(4)
    t = np.sort(rng.uniform(0.0, 0.4, 120))
    trend = 0.6 * np.sin(2 * np.pi * t / 0.8)          # slow, big
    y = 16.0 + trend + 0.08 * np.sin(2 * np.pi * FREQ_CD * t) + rng.normal(0, 0.01, 120)
    out = p3.joint_gp_fit(t, y, np.full(120, 0.01), FREQ_CD)
    assert out["amplitude"] == pytest.approx(0.08, rel=0.25)


def test_joint_gp_fit_refuses_a_stub():
    out = p3.joint_gp_fit([0.0, 1.0], [1.0, 2.0], [0.1, 0.1], 12.6)
    assert math.isnan(out["amplitude"])
    assert out["backend"] == "none"


def test_running_median_detrend_removes_a_constant():
    t = np.linspace(0.0, 1.0, 50)
    out = p3.running_median_detrend(t, np.full(50, 7.0), 0.2)
    assert out == pytest.approx(0.0, abs=1e-12)


def test_detrend_then_search_destroys_what_a_joint_fit_keeps():
    """THE demonstration, as a unit test.  A smoothing window shorter than
    the signal period must eat the signal; the joint fit must not."""
    t = np.sort(np.random.default_rng(5).uniform(0.0, 0.4, 100))
    out = p3.detrend_suppression(t, np.full(t.size, 0.01), FREQ_CD,
                                 amplitude_mag=0.10,
                                 window_d=0.5 * PERIOD_D, n_real=4)
    assert out["frac_detrend"] < 0.5           # most of the signal gone
    assert out["frac_joint"] > 0.85            # signal preserved
    assert out["frac_joint"] > out["frac_detrend"]


def test_detrend_can_also_FABRICATE_amplitude():
    """Not just attenuation.  At a window of about 1.5 periods the running
    median leaves MORE amplitude than was injected — the analyst who does
    not know the true amplitude cannot tell which way it went, which is the
    actual argument against detrend-then-search."""
    t = np.sort(np.random.default_rng(5).uniform(0.0, 0.4, 100))
    out = p3.detrend_suppression(t, np.full(t.size, 0.01), FREQ_CD,
                                 amplitude_mag=0.10,
                                 window_d=1.5 * PERIOD_D, n_real=4)
    assert out["frac_detrend"] > 1.05
    assert out["frac_joint"] == pytest.approx(1.0, abs=0.12)


def test_joint_fit_is_insensitive_to_the_window_the_detrend_needs():
    """The joint fit has no window to choose, so it returns the same answer
    where the detrend swings from 0.27 to 1.25."""
    t = np.sort(np.random.default_rng(5).uniform(0.0, 0.4, 100))
    runs = [p3.detrend_suppression(t, np.full(t.size, 0.01), FREQ_CD, 0.10,
                                   window_d=w * PERIOD_D, n_real=3)
            for w in (0.5, 1.5, 10.0)]
    detrend = [r["frac_detrend"] for r in runs]
    joint = [r["frac_joint"] for r in runs]
    assert max(detrend) - min(detrend) > 0.5
    assert max(joint) - min(joint) < 0.02


def test_detrend_suppression_reports_the_window_in_periods():
    t = np.sort(np.random.default_rng(6).uniform(0.0, 0.4, 60))
    out = p3.detrend_suppression(t, np.full(t.size, 0.02), FREQ_CD,
                                 0.1, window_d=2 * PERIOD_D, n_real=2)
    assert out["window_periods"] == pytest.approx(2.0, rel=1e-9)


def test_lsq_amplitude_recovers_a_noiseless_sinusoid():
    t = np.linspace(0.0, 0.4, 200)
    y = 16.0 + 0.13 * np.sin(2 * np.pi * FREQ_CD * t + 0.7)
    assert p3._lsq_amplitude(t, y, np.full(200, 0.01), FREQ_CD) == pytest.approx(0.13, rel=1e-6)


def test_joint_gp_fit_is_immune_to_the_time_origin():
    """Regression: BJD-scale absolute times broke the celerite2 backend.

    These are BJD values near 2.46e6 while the correlation lengths fitted
    are ~0.1 d, so the structure lives in the 8th significant figure.  The
    dense path only ever forms differences and survived; celerite2 keeps the
    absolute times and its factorisation failed at some hyper-parameters —
    silently, because ``quiet=True`` leaves a half-built object instead of
    raising.  Those grid points vanished from one backend's search and not
    the other's, and the two disagreed by 3.6% in recovered amplitude.
    """
    pytest.importorskip("celerite2")
    rng = np.random.default_rng(12)
    t0 = np.sort(rng.uniform(0.0, 0.4, 80))
    y = 16.0 + 0.09 * np.sin(2 * np.pi * FREQ_CD * t0) + rng.normal(0, 0.01, 80)
    e = np.full(80, 0.01)
    based = p3.joint_gp_fit(t0, y, e, FREQ_CD)
    shifted = p3.joint_gp_fit(t0 + 2460734.0, y, e, FREQ_CD)
    # Not bit-identical, and cannot be: storing t + 2460734.0 in a float64
    # costs about seven decimal digits of the fractional part, so the two
    # runs genuinely see slightly different timestamps.  1e-8 is comfortably
    # inside that and comfortably outside the 3.6% the bug produced.
    assert shifted["amplitude"] == pytest.approx(based["amplitude"], rel=1e-8)
    assert shifted["loglike"] == pytest.approx(based["loglike"], rel=1e-8)


def test_joint_gp_backends_agree_on_bjd_scale_times():
    """The same regression, stated as the property that actually matters."""
    pytest.importorskip("celerite2")
    rng = np.random.default_rng(13)
    t = 2460734.0 + np.sort(rng.uniform(0.0, 0.39, 90))
    y = 16.0 + 0.10 * np.sin(2 * np.pi * FREQ_CD * t) + rng.normal(0, 0.02, 90)
    e = np.full(90, 0.02)
    fast = p3.joint_gp_fit(t, y, e, FREQ_CD, use_celerite=True)
    slow = p3.joint_gp_fit(t, y, e, FREQ_CD, use_celerite=False)
    assert fast["backend"] == "celerite2"
    assert fast["amplitude"] == pytest.approx(slow["amplitude"], rel=1e-8)
    assert fast["sigma"] == pytest.approx(slow["sigma"], rel=1e-12)
    assert fast["rho"] == pytest.approx(slow["rho"], rel=1e-12)


def test_joint_gp_fit_survives_a_hyperparameter_that_breaks_celerite2():
    """A grid point celerite2 cannot factorise must be dropped, not crash.

    The guard has to wrap the whole celerite2 block: with ``quiet=True`` the
    failure surfaces out of ``apply_inverse`` as an AttributeError on a
    private array, so a try/except around ``compute`` alone catches nothing.
    """
    pytest.importorskip("celerite2")
    rng = np.random.default_rng(14)
    t = np.sort(rng.uniform(0.0, 0.4, 60))
    y = 16.0 + rng.normal(0, 0.01, 60)
    e = np.full(60, 0.01)
    # A grid deliberately containing degenerate corners (zero-ish and huge
    # correlation lengths) alongside one usable point.
    out = p3.joint_gp_fit(t, y, e, FREQ_CD,
                          rho_grid_d=[1e-12, 0.05, 1e9],
                          sigma_grid=[1e-12, 0.02, 1e9])
    assert np.isfinite(out["amplitude"])
    assert out["backend"] == "celerite2"
