"""Unit tests for :mod:`macro_phot.final_science` — the CV-S10 arithmetic.

Each test asserts a PROPERTY rather than a remembered number wherever one is
available: a hump injected at a known amplitude and phase comes back at that
amplitude and phase; a structure function built from white noise of known
sigma returns that sigma at every lag; a flickering amplitude added in
quadrature to a known floor is recovered by subtracting the floor in
quadrature.  The handful of tests that pin a number are the ones where the
number IS the rule — the coverage gate, the drift bar, the verdict words.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import final_science as fs        # noqa: E402

PERIOD_D = 0.0868
EPOCH = 2_460_000.0


def _run_times(hours: float = 5.0, cadence_s: float = 190.0,
               t0: float = 2_460_360.7) -> np.ndarray:
    """Timestamps that look like a real YZ Cnc run: one filter of a
    three-filter cycle, so ~190 s between points on the target."""
    n = int(hours * 3600.0 / cadence_s)
    return t0 + np.arange(n) * cadence_s / 86400.0


# ===========================================================================
# 1.  Run selection
# ===========================================================================
def test_select_dense_runs_keeps_only_dense_and_sorts_by_local_night():
    rows = [
        {"local_night": "2024-05-01", "is_dense": 1, "state": "QUIESCENT"},
        {"local_night": "2024-02-20", "is_dense": 1, "state": "QUIESCENT"},
        {"local_night": "2024-03-09", "is_dense": 0, "state": "OUTBURST"},
        {"local_night": "2024-02-21", "is_dense": 1, "state": "OUTBURST"},
    ]
    got = fs.select_dense_runs(rows)
    assert [r["local_night"] for r in got] == ["2024-02-20", "2024-02-21",
                                               "2024-05-01"]


def test_select_dense_runs_filters_on_state():
    rows = [
        {"local_night": "2024-05-01", "is_dense": 1, "state": "QUIESCENT"},
        {"local_night": "2024-02-21", "is_dense": 1, "state": "OUTBURST"},
    ]
    assert len(fs.select_dense_runs(rows, "QUIESCENT")) == 1
    assert len(fs.select_dense_runs(rows, "OUTBURST")) == 1
    assert fs.select_dense_runs(rows, "ELEVATED") == []


# ===========================================================================
# 2.  Folding and the orbital hump
# ===========================================================================
def test_orbital_phase_is_in_unit_interval_and_periodic():
    t = _run_times()
    ph = fs.orbital_phase(t, PERIOD_D, EPOCH)
    assert ph.min() >= 0.0 and ph.max() < 1.0
    # One whole period later is the same phase.
    ph2 = fs.orbital_phase(t + PERIOD_D, PERIOD_D, EPOCH)
    assert np.allclose(ph, ph2, atol=1e-9)


def test_harmonic_design_shape_with_and_without_nightly_constants():
    ph = np.linspace(0, 1, 40, endpoint=False)
    assert fs.harmonic_design(ph, 2).shape == (40, 5)
    nights = np.array(["a"] * 20 + ["b"] * 20)
    # Two nightly constants replace the single one: 2 + 2*2 = 6 columns.
    assert fs.harmonic_design(ph, 2, nights).shape == (40, 6)


def test_amplitude_and_phase_round_trips():
    for amp, phase in [(0.05, 0.0), (0.12, 0.25), (0.3, 0.6), (0.02, 0.99)]:
        a = amp * math.cos(2 * math.pi * phase)
        b = amp * math.sin(2 * math.pi * phase)
        got_a, got_p = fs.amplitude_and_phase(a, b)
        assert got_a == pytest.approx(amp)
        assert got_p == pytest.approx(phase, abs=1e-9)


def test_fold_fit_recovers_an_injected_hump():
    """The headline property: inject a known hump, get it back."""
    t = _run_times()
    rng = np.random.default_rng(7)
    amp, phase = 0.08, 0.35
    # Magnitudes: BRIGHTEST at `phase`, so the magnitude model is -A cos().
    m = 14.3 - amp * np.cos(2 * np.pi * (fs.orbital_phase(t, PERIOD_D, EPOCH)
                                         - phase))
    e = np.full(t.size, 0.01)
    m = m + rng.normal(0.0, e)
    fit = fs.fold_fit(t, m, e, PERIOD_D, EPOCH)
    assert fit["amp"] == pytest.approx(amp, rel=0.10)
    assert fit["phase_max"] == pytest.approx(phase, abs=0.02)
    assert fit["mean_mag"] == pytest.approx(14.3, abs=0.01)
    # And chi2/dof is near 1 when the model IS the truth.
    assert 0.5 < fit["chi2nu"] < 2.0


def test_fold_fit_nightly_constants_absorb_a_night_to_night_step():
    """Two nights with a 0.3 mag offset and the SAME hump.

    Without one constant per night the step leaks into the fundamental.
    With them, the amplitude comes back.  This is the whole reason blocks
    are fitted jointly rather than mean-subtracted first.
    """
    t1 = _run_times(2.5, 145.0, 2_460_432.7)
    t2 = _run_times(2.0, 145.0, 2_460_433.7)
    t = np.concatenate([t1, t2])
    nights = np.array(["n1"] * t1.size + ["n2"] * t2.size)
    amp, phase = 0.06, 0.2
    step = np.where(nights == "n2", 0.30, 0.0)
    m = (14.3 + step
         - amp * np.cos(2 * np.pi * (fs.orbital_phase(t, PERIOD_D, EPOCH)
                                     - phase)))
    e = np.full(t.size, 0.01)
    joint = fs.fold_fit(t, m, e, PERIOD_D, EPOCH, night_index=nights)
    shared = fs.fold_fit(t, m, e, PERIOD_D, EPOCH)
    assert joint["amp"] == pytest.approx(amp, rel=0.05)
    assert joint["phase_max"] == pytest.approx(phase, abs=0.02)
    # The shared-constant fit is measurably wrong about the amplitude.
    assert abs(shared["amp"] - amp) > abs(joint["amp"] - amp)


def test_fold_fit_returns_residuals_aligned_to_the_input():
    t = _run_times()
    m = 14.3 + 0.05 * np.sin(2 * np.pi * t / PERIOD_D)
    fit = fs.fold_fit(t, m, np.full(t.size, 0.01), PERIOD_D, EPOCH)
    assert fit["resid"].shape == t.shape
    assert np.std(fit["resid"]) < 1e-6      # the model is exact here


def test_weighted_lstsq_returns_nan_when_underdetermined():
    X = np.ones((3, 5))
    coef, cov, chi2 = fs.weighted_lstsq(X, np.ones(3), np.ones(3))
    assert not np.isfinite(coef).any()
    assert not np.isfinite(chi2)


def test_amplitude_sigma_matches_a_direct_monte_carlo():
    """Linear propagation on A = hypot(a, b), checked by sampling."""
    a, b = 0.08, 0.03
    va, vb, vab = 4e-6, 9e-6, 1e-6
    cov = np.array([[va, vab], [vab, vb]])
    rng = np.random.default_rng(3)
    draws = rng.multivariate_normal([a, b], cov, size=40_000)
    mc = np.std(np.hypot(draws[:, 0], draws[:, 1]))
    assert fs.amplitude_sigma(a, b, va, vb, vab) == pytest.approx(mc,
                                                                  rel=0.05)


def test_circular_mean_handles_the_wrap_at_zero():
    """The ordinary mean of 0.98 and 0.02 is 0.5 — the point furthest from
    both.  The circular mean is 0.0, which is the whole reason this function
    exists."""
    mean, spread = fs.circular_mean_and_spread([0.98, 0.02])
    # Compared ON THE CIRCLE: the answer may legitimately come back as
    # 0.0 or as 0.999999, and those are the same point.
    assert abs(fs.phase_difference(mean, 0.0)) < 1e-9
    assert spread < 0.05


def test_circular_spread_is_zero_for_identical_phases_and_grows():
    _m, s0 = fs.circular_mean_and_spread([0.3, 0.3, 0.3])
    _m, s1 = fs.circular_mean_and_spread([0.28, 0.30, 0.32])
    _m, s2 = fs.circular_mean_and_spread([0.1, 0.4, 0.8])
    assert s0 == pytest.approx(0.0, abs=1e-12)
    assert s0 < s1 < s2


def test_phase_difference_wraps_into_half_a_cycle():
    assert fs.phase_difference(0.98, 0.02) == pytest.approx(-0.04)
    assert fs.phase_difference(0.02, 0.98) == pytest.approx(0.04)
    assert abs(fs.phase_difference(0.77, 0.40)) == pytest.approx(0.37)
    assert -0.5 <= fs.phase_difference(0.9, 0.1) < 0.5


def test_phase_drift_grows_linearly_with_baseline():
    d1 = fs.phase_drift_cycles(1.0, PERIOD_D, 5e-5)
    d71 = fs.phase_drift_cycles(71.0, PERIOD_D, 5e-5)
    assert d71 == pytest.approx(71.0 * d1)
    # The two real cases the pipeline turns on: one night apart is foldable,
    # 71 days apart is not.
    assert d1 < fs.PHASE_DRIFT_BAR_CYCLES
    assert d71 > fs.PHASE_DRIFT_BAR_CYCLES


# ===========================================================================
# 3.  Flickering
# ===========================================================================
def test_structure_function_of_white_noise_returns_its_sigma():
    """For uncorrelated noise the structure function is flat at sigma."""
    t = _run_times(8.0, 100.0)
    rng = np.random.default_rng(11)
    sigma = 0.02
    y = rng.normal(0.0, sigma, t.size)
    edges = fs.log_tau_edges(150.0, 10_000.0, 6)
    tau, sf, npair = fs.structure_function(t, y, edges)
    ok = np.isfinite(sf)
    assert ok.sum() >= 4
    assert np.allclose(sf[ok], sigma, rtol=0.25)


def test_structure_function_returns_nan_for_thin_bins():
    t = _run_times(1.0, 200.0)
    y = np.zeros(t.size)
    edges = fs.log_tau_edges(1.0, 10.0, 3)     # no pairs live down here
    _tau, sf, npair = fs.structure_function(t, y, edges)
    assert npair.sum() == 0
    assert not np.isfinite(sf).any()


def test_structure_function_rises_for_a_random_walk():
    """A red process must give a RISING structure function; that is the
    only thing that separates flickering from white noise here."""
    t = _run_times(8.0, 100.0)
    rng = np.random.default_rng(5)
    y = np.cumsum(rng.normal(0.0, 0.01, t.size))
    edges = fs.log_tau_edges(150.0, 10_000.0, 6)
    _tau, sf, _n = fs.structure_function(t, y, edges)
    ok = np.isfinite(sf)
    assert sf[ok][-1] > 2.0 * sf[ok][0]


def test_quadrature_excess_recovers_an_injected_flicker_amplitude():
    total = np.hypot(np.array([0.05, 0.08]), np.array([0.02, 0.02]))
    got = fs.quadrature_excess(total, np.array([0.02, 0.02]))
    assert got == pytest.approx([0.05, 0.08])


def test_quadrature_excess_is_nan_not_zero_below_the_floor():
    got = fs.quadrature_excess(np.array([0.01, 0.03]), np.array([0.02, 0.02]))
    assert not np.isfinite(got[0])          # below the floor: NOT MEASURED
    assert np.isfinite(got[1])


def test_excess_significance_is_zero_when_target_equals_floor():
    z = fs.excess_significance(np.array([0.02]), np.array([0.002]),
                               np.array([0.02]), np.array([0.002]))
    assert z[0] == pytest.approx(0.0)


def test_excess_significance_grows_with_the_excess():
    z = fs.excess_significance(np.array([0.02, 0.06]),
                               np.array([0.002, 0.002]),
                               np.array([0.02, 0.02]),
                               np.array([0.002, 0.002]))
    assert z[1] > z[0]


def test_sf_sigma_is_capped_by_the_point_count_not_the_pair_count():
    """Pairs are not independent; the conservative n_eff is the cap."""
    s = fs.sf_sigma(np.array([0.1]), np.array([5000]), n_points=50)
    assert s[0] == pytest.approx(0.1 / math.sqrt(2 * 50))


def test_median_floor_ignores_one_variable_interloper():
    """A single variable star in the matched sample must not move the
    floor: that is why the floor is a median and not a mean."""
    quiet = [np.array([0.010, 0.011, 0.012]) for _ in range(9)]
    loud = [np.array([0.20, 0.22, 0.24])]
    med, sig = fs.median_floor(quiet + loud)
    assert med == pytest.approx([0.010, 0.011, 0.012])
    assert np.isfinite(sig).all()


def test_roll_within_nights_never_moves_a_point_between_nights():
    v = np.array([1.0, 2.0, 3.0, 10.0, 20.0, 30.0])
    nights = np.array(["a", "a", "a", "b", "b", "b"])
    rng = np.random.default_rng(2)
    for _ in range(20):
        out = fs.roll_within_nights(v, nights, rng)
        assert set(out[:3]) == {1.0, 2.0, 3.0}
        assert set(out[3:]) == {10.0, 20.0, 30.0}


def test_roll_within_nights_on_one_night_is_a_plain_cyclic_roll():
    v = np.arange(6, dtype=float)
    nights = np.array(["a"] * 6)
    out = fs.roll_within_nights(v, nights, np.random.default_rng(1))
    assert sorted(out) == sorted(v)
    # A cyclic roll preserves every adjacent difference except at the seam.
    diffs = np.diff(np.concatenate([out, out[:1]]))
    assert (diffs == 1.0).sum() == 5


# ===========================================================================
# 4.  Outburst runs
# ===========================================================================
def test_linear_rate_recovers_an_injected_slope():
    t = _run_times(5.0, 190.0)
    rate = -0.05                                   # mag per hour, rising
    m = 13.0 + rate * (t - t.mean()) * 24.0
    got, sig = fs.linear_rate_per_hour(t, m)
    assert got == pytest.approx(rate, rel=1e-6)
    assert sig == pytest.approx(0.0, abs=1e-9)


def test_linear_rate_error_bar_covers_the_truth_on_noisy_data():
    t = _run_times(5.0, 190.0)
    rng = np.random.default_rng(19)
    rate = 0.03
    m = 13.0 + rate * (t - t.mean()) * 24.0 + rng.normal(0, 0.05, t.size)
    got, sig = fs.linear_rate_per_hour(t, m)
    assert abs(got - rate) < 3.0 * sig


def test_percentile_amplitude_ignores_two_outliers():
    m = np.concatenate([np.full(100, 14.0), [10.0, 18.0]])
    # max-min would say 8 mag; p5-p95 says essentially nothing varies.
    assert fs.percentile_amplitude(m) < 0.5


# ===========================================================================
# 5.  Gates and verdicts
# ===========================================================================
def test_snr_gate_passes_and_fails_around_its_ratio():
    assert fs.snr_gate(0.010, 0.060, 5.0)["passes"]
    assert not fs.snr_gate(0.010, 0.040, 5.0)["passes"]
    assert fs.snr_gate(0.010, 0.050, 5.0)["margin_mag"] == pytest.approx(0.0)


def test_snr_gate_refuses_nonsense_inputs():
    assert not fs.snr_gate(0.0, 0.05)["passes"]
    assert not fs.snr_gate(float("nan"), 0.05)["passes"]


def test_detection_call_needs_both_tests():
    assert fs.detection_call(0.10, 0.05, 0.4, 0.2) == "DETECTED"
    assert fs.detection_call(0.10, 0.05, 0.1, 0.2) == "MARGINAL"
    assert fs.detection_call(0.02, 0.05, 0.4, 0.2) == "MARGINAL"
    assert fs.detection_call(0.02, 0.05, 0.1, 0.2) == "NOT DETECTED"
    assert fs.detection_call(0.10, float("nan"), 0.4, 0.2) == "UNTESTED"


def test_detection_call_coverage_gate_overrides_everything():
    """Under 1.5 orbits, a modulation at P_orb and a trend are the same
    statement, so no score may promote the answer to DETECTED."""
    assert fs.detection_call(0.10, 0.05, 0.9, 0.2, cycles=0.9) \
        == "AMPLITUDE ONLY"
    assert fs.detection_call(0.10, 0.05, 0.9, 0.2, cycles=2.0) == "DETECTED"


def test_capability_verdict_both_directions():
    assert fs.capability_verdict(9, 3) == "SUPPORTED"
    assert fs.capability_verdict(2, 3) == "NOT SUPPORTED"
    # Smaller-is-better, e.g. a duty-cycle interval half-width.
    assert fs.capability_verdict(12.0, 15.0, higher_is_better=False) \
        == "SUPPORTED"
    assert fs.capability_verdict(18.0, 15.0, higher_is_better=False) \
        == "NOT SUPPORTED"
    assert fs.capability_verdict(float("nan"), 3) == "UNTESTED"


def test_nights_needed_reports_zero_when_the_bar_is_already_met():
    got = fs.nights_needed(9, 3, 4.5)
    assert got["shortfall"] == 0 and got["seasons"] == 0.0


def test_nights_needed_converts_a_shortfall_into_seasons():
    got = fs.nights_needed(4, 8, 2.0)
    assert got["shortfall"] == 4
    assert got["seasons"] == pytest.approx(2.0)


def test_log_tau_edges_are_log_spaced_and_bracket_the_range():
    e = fs.log_tau_edges(60.0, 14_400.0, 9)
    assert e.size == 10
    assert e[0] == pytest.approx(60.0) and e[-1] == pytest.approx(14_400.0)
    ratios = e[1:] / e[:-1]
    assert np.allclose(ratios, ratios[0])


# ---------------------------------------------------------------------------
# run_night_label -- one run, one name, in Table 4 and in Figure 11
# ---------------------------------------------------------------------------
def test_run_night_label_prefers_the_utc_night():
    """p4_run.scope keys a run by its LOCAL observing night and utc_nights
    records the UTC one.  Table 4 named the uninformative YZ Cnc run
    2024-02-20 and Figure 11 named the same run 2024-02-21."""
    assert fs.run_night_label("2024-02-21", "2024-02-20") == "2024-02-21"


def test_run_night_label_falls_back_to_the_observing_night():
    """A row written before utc_nights was populated must still get a name
    rather than an empty string in a caption."""
    assert fs.run_night_label(None, "2024-02-20") == "2024-02-20"
    assert fs.run_night_label("", "2024-02-20") == "2024-02-20"


def test_run_night_label_names_a_multi_night_block_by_its_first_night():
    assert fs.run_night_label("2024-05-02+2024-05-03", None) == "2024-05-02"


def test_run_night_label_on_nothing_is_empty_not_none():
    """A caption concatenates this; None would print as the word 'None'."""
    assert fs.run_night_label(None, None) == ""
