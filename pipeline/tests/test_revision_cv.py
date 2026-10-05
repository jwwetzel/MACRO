"""Unit tests for :mod:`macro_phot.revision_cv` — the committee revision's
arithmetic (CV-R1...R8), and invariants of the ``rv_`` products it wrote.

The module under test exists because four inferences in the CV paper were
arithmetically right and inferentially wrong.  So these tests are built
the way the review was: each one sets up the situation in which the OLD
method gave the wrong answer — over-stated error bars, a magnitude-space
ramp on bands of different depth, a pairing window on a steep light curve,
a recovery contour read upside down — and demands the right answer back.

Grouped in the order of the module.  The last section checks the built
products and skips when they are absent, like ``test_cv_products.py``.
"""

from __future__ import annotations

import math
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import phase3 as p3            # noqa: E402
from macro_phot import revision_cv as rv       # noqa: E402

PERIOD_D = 0.07908912
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PHOT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"


# ===========================================================================
# 1.  Small-sample distributions
# ===========================================================================
@pytest.mark.parametrize("t, dof, p", [
    (2.228, 10, 0.05), (12.706, 1, 0.05), (2.086, 20, 0.05),
    (3.169, 10, 0.01), (0.0, 7, 1.0)])
def test_student_t_matches_the_tables(t, dof, p):
    assert rv.student_t_two_sided_p(t, dof) == pytest.approx(p, abs=2e-4)


@pytest.mark.parametrize("chi2, dof, p", [
    (3.841, 1, 0.05), (18.307, 10, 0.05), (6.635, 1, 0.01),
    (43.773, 30, 0.05), (0.0, 4, 1.0)])
def test_chi2_tail_matches_the_tables(chi2, dof, p):
    assert rv.chi2_sf(chi2, dof) == pytest.approx(p, abs=2e-4)


def test_sign_test_is_the_exact_binomial():
    # 0 of 5 negative: two-sided 2 * (1/32).
    assert rv.sign_test_p(0, 5) == pytest.approx(2.0 / 32.0)
    # 9 of 9 negative, the data scientist's own number: 2/512 = 0.0039.
    assert rv.sign_test_p(9, 9) == pytest.approx(2.0 / 512.0)
    # 1 of 10: 2 * (1 + 10) / 1024.
    assert rv.sign_test_p(1, 10) == pytest.approx(22.0 / 1024.0)
    assert rv.sign_test_p(5, 10) == 1.0
    assert math.isnan(rv.sign_test_p(0, 0))


def test_wilcoxon_is_exact_and_handles_ties():
    # All five positive and distinct: only the two extreme patterns reach
    # the observed rank sum, so p = 2 / 2^5.
    assert rv.wilcoxon_signed_rank_p([1, 2, 3, 4, 5]) == pytest.approx(
        2.0 / 32.0)
    assert rv.wilcoxon_signed_rank_p([-1, -2, -3, -4, -5, -6, -7, -8, -9,
                                      -10]) == pytest.approx(2.0 / 1024.0)
    # Symmetric about zero: no evidence at all.
    assert rv.wilcoxon_signed_rank_p([-3, -2, -1, 1, 2, 3]) == 1.0
    # Ties get average ranks and the enumeration is over THOSE ranks:
    # |d| = 1,1,2 -> ranks 1.5,1.5,3; all positive -> 2/8.
    assert rv.wilcoxon_signed_rank_p([1, 1, 2]) == pytest.approx(0.25)
    # Zeros are dropped, not ranked.
    assert rv.wilcoxon_signed_rank_p([0, 0, 1, 2, 3, 4, 5]) == pytest.approx(
        2.0 / 32.0)


def test_wilcoxon_agrees_with_scipy_where_scipy_is_exact():
    stats = pytest.importorskip("scipy.stats")
    rng = np.random.default_rng(2)
    for n in (6, 9, 14):
        d = rng.normal(0.4, 1.0, n)
        want = stats.wilcoxon(d, mode="exact").pvalue
        assert rv.wilcoxon_signed_rank_p(d) == pytest.approx(want, rel=1e-9)


def test_signflip_is_exact_and_two_sided():
    # Nine units all of one sign: only the identity and its mirror are as
    # extreme, so p = 2 / 2^9.
    p, basis = rv.signflip_p([-10, -20, -30, -40, -50, -60, -70, -80, -90])
    assert p == pytest.approx(2.0 / 512.0)
    assert basis.startswith("exact") and "512" in basis
    p_pos, _ = rv.signflip_p([10, 20, 30, 40, 50, 60, 70, 80, 90])
    assert p_pos == p
    # Perfectly symmetric data: every pattern is at least as extreme as a
    # mean of zero.
    assert rv.signflip_p([-2, -1, 1, 2])[0] == 1.0


def test_signflip_by_cluster_uses_nights_not_cycles():
    """Twelve cycles on three nights are three independent units.  The
    cycle-level test credits itself with twelve and reports p = 2/4096;
    the night-level one can do no better than 2/8, and says so."""
    d = np.full(12, -100.0) + np.arange(12)
    nights = np.repeat(["a", "b", "c"], 4)
    p_cycle, _ = rv.signflip_p(d)
    p_night, _ = rv.signflip_p(d, clusters=nights)
    assert p_cycle == pytest.approx(2.0 / 4096.0)
    assert p_night == pytest.approx(2.0 / 8.0)


def test_signflip_monte_carlo_branch_agrees_with_the_normal_limit():
    rng = np.random.default_rng(8)
    d = rng.normal(0.5, 1.0, rv.EXACT_PERM_MAX_N + 15)
    p, basis = rv.signflip_p(d, seed=1)
    assert basis.startswith("Monte Carlo")
    t = d.mean() / (d.std(ddof=1) / math.sqrt(d.size))
    assert p == pytest.approx(rv.student_t_two_sided_p(t, d.size - 1),
                              abs=0.01)


def test_paired_tests_find_an_offset_the_budget_hides():
    """The finding in miniature (DS.F1 / RF.B1).  Twelve pairs offset by
    -100 s with 30 s of real scatter and assigned errors of 180 s per pair.
    The budget-weighted statistic is 'under two sigma'; the scatter-based
    one is more than eight; chi2nu says the budget is the one to disbelieve;
    and every distribution-free test agrees with the scatter."""
    rng = np.random.default_rng(1)
    d = -100.0 + rng.normal(0.0, 30.0, 12)
    out = rv.paired_tests(d, clusters=np.repeat(np.arange(6), 2),
                          sigma_budget=np.full(12, 180.0), n_trials=3)
    assert out["n"] == 12 and out["n_clusters"] == 6
    assert out["mean"] == pytest.approx(d.mean())
    assert out["se"] == pytest.approx(d.std(ddof=1) / math.sqrt(12))
    assert abs(out["mean"]) / out["se"] > 8.0
    assert abs(out["mean_budget"]) / out["sigma_budget"] < 2.5
    assert out["sigma_budget"] == pytest.approx(180.0 / math.sqrt(12))
    assert out["dof"] == 11 and out["chi2nu"] < 0.1       # reported, unclipped
    assert out["n_negative"] == 12
    assert out["p_sign"] == pytest.approx(2.0 / 4096.0)
    assert out["p_perm"] == pytest.approx(2.0 / 4096.0)
    assert out["p_wilcoxon"] == pytest.approx(2.0 / 4096.0)
    assert out["p_perm_cluster"] == pytest.approx(2.0 / 64.0)
    assert out["p_perm_bonf"] == pytest.approx(3 * out["p_perm"])
    assert out["p_perm_cluster_bonf"] == pytest.approx(3 * 2.0 / 64.0)


def test_paired_tests_do_not_invent_an_offset():
    rng = np.random.default_rng(4)
    d = rng.normal(0.0, 50.0, 14)
    out = rv.paired_tests(d, n_trials=3)
    assert out["p_perm"] > 0.05 and out["p_t"] > 0.05
    assert out["p_perm_bonf"] == pytest.approx(min(1.0, 3 * out["p_perm"]))


def test_paired_tests_survive_one_pair_and_none():
    one = rv.paired_tests([-80.0])
    assert one["n"] == 1 and math.isnan(one["se"]) and math.isnan(one["t"])
    none = rv.paired_tests([])
    assert none["n"] == 0 and math.isnan(none["mean"])


def test_inverse_variance_combination_and_its_heterogeneity():
    c = rv.combine_inverse_variance([-100.0, -120.0], [20.0, 40.0])
    w = np.array([1 / 400.0, 1 / 1600.0])
    assert c["mean"] == pytest.approx(np.sum(w * [-100, -120]) / w.sum())
    assert c["se"] == pytest.approx(math.sqrt(1 / w.sum()))
    assert c["dof_het"] == 1 and c["p_het"] > 0.5        # the eras agree
    far = rv.combine_inverse_variance([-100.0, +100.0], [20.0, 20.0])
    assert far["p_het"] < 1e-6                           # ... and when not


# ===========================================================================
# 2.  The edge estimators
# ===========================================================================
def _window(t_edge=0.010, width_s=300.0, depth=1.2, cadence_s=219.0,
            noise=0.0, seed=0, n=11):
    """A window of ``n`` points at the real cadence, edge mid-window, in
    MAGNITUDES whose FLUX falls linearly (the achromatic-limb template)."""
    step = cadence_s / 86400.0
    t = t_edge + step * (np.arange(n) - (n - 1) / 2.0 + 0.31)
    flux = rv.flux_edge_template(t, t_edge, width_s / 86400.0, depth)
    mag = 16.0 - 2.5 * np.log10(flux)
    if noise:
        mag = mag + np.random.default_rng(seed).normal(0.0, noise, n)
    return t, mag, np.full(n, max(noise, 0.01))


def test_profile_fit_has_the_published_fits_chi2_surface():
    """``fit_edge_profile`` is the vectorised form of
    ``phase3.fit_edge``: same model, same minimum, same epoch, same
    acceptance, on a problem with a unique minimum."""
    t, mag, e = _window(noise=0.03, seed=3)
    tg = p3.edge_time_grid(0.010, 3 * 219.0 / 86400.0, 361)
    wg = np.array([60.0, 120.0, 240.0, 480.0]) / 86400.0
    old = p3.fit_edge(t, mag, e, tg, wg, 219.0)
    new = rv.fit_edge_profile(t, mag, e, tg, wg, 219.0, tie="first")
    assert new["chi2nu"] == pytest.approx(old.chi2nu, rel=1e-8)
    assert new["dof"] == t.size - 4
    assert new["chi2"] == pytest.approx(old.chi2nu * new["dof"], rel=1e-8)
    assert new["width_d"] == old.width_d
    assert new["accepted"] == old.accepted
    if new["valley_s"] == 0.0:
        assert new["t_edge_d"] == pytest.approx(old.t_edge_d, abs=1e-12)
        assert new["level_bright"] == pytest.approx(old.level_bright,
                                                    abs=1e-8)
        assert new["step"] == pytest.approx(old.depth_mag, abs=1e-8)


def test_a_flat_valley_is_reported_and_its_midpoint_returned():
    """A step sharper than the cadence, fitted with a narrow ramp, leaves
    chi-squared flat across the gap between the two exposures that bracket
    it.  The data say 'somewhere in this gap': the symmetric answer is the
    middle, and the width of the indifference is part of the result."""
    cad = 219.0 / 86400.0
    t = cad * np.arange(12)
    t_true = 5.4 * cad
    y = 16.0 + 1.0 * (t > t_true)
    e = np.full(t.size, 0.02)
    tg = np.linspace(4.0 * cad, 7.0 * cad, 601)
    wg = np.array([20.0 / 86400.0])
    mid = rv.fit_edge_profile(t, y, e, tg, wg, 219.0, tie="centre")
    first = rv.fit_edge_profile(t, y, e, tg, wg, 219.0, tie="first")
    # The valley spans the gap less the ramp width, to the grid step.
    assert mid["valley_s"] == pytest.approx(219.0 - 20.0, abs=2.5)
    assert mid["t_edge_d"] == pytest.approx(5.5 * cad, abs=1.5 / 86400.0)
    assert first["t_edge_d"] < mid["t_edge_d"]
    assert (mid["t_edge_d"] - first["t_edge_d"]) * 86400.0 == pytest.approx(
        mid["valley_s"] / 2.0, abs=1.0)
    with pytest.raises(ValueError):
        rv.fit_edge_profile(t, y, e, tg, wg, 219.0, tie="last")


def test_profile_fit_refuses_what_the_published_fit_refuses():
    few = rv.fit_edge_profile([0.0, 1.0], [1.0, 2.0], [0.1, 0.1],
                              np.array([0.5]), np.array([0.1]), 219.0)
    assert not few["accepted"] and "usable points" in few["reason"]
    # A step lost in the noise fails the S/N gate ...
    t, mag, e = _window(depth=0.05, noise=0.05, seed=5)
    tg = rv.v2_time_grid(0.010, 219.0, PERIOD_D)
    wg = np.array(rv.V2_WIDTH_GRID_S) / 86400.0
    weak = rv.fit_edge_profile(t, mag, e, tg, wg, 219.0)
    assert not weak["accepted"] and "SNR" in weak["reason"]
    # ... unless the caller says the edge's existence is already settled.
    assert rv.fit_edge_profile(t, mag, e, tg, wg, 219.0,
                               min_snr=0.0)["reason"] != weak["reason"]


def test_v2_time_grid_is_centred_stepped_and_inside_the_window():
    tg = rv.v2_time_grid(100.0, 219.0, PERIOD_D)
    assert tg[tg.size // 2] == pytest.approx(100.0)
    assert np.diff(tg) * 86400.0 == pytest.approx(rv.V2_TIME_STEP_S,
                                                  abs=1e-3)
    assert (tg[-1] - tg[0]) * 86400.0 / 2 <= 3 * 219.0 + 1e-6
    # A long cadence must not push the grid outside the fitting window.
    wide = rv.v2_time_grid(100.0, 5000.0, PERIOD_D)
    assert (wide[-1] - 100.0) <= p3.EDGE_WINDOW_PHASE * PERIOD_D + 1e-9


def test_magnitude_midpoint_shift_is_the_geometric_mean():
    assert rv.magnitude_midpoint_shift(1e-9) == pytest.approx(0.5, abs=1e-6)
    for dm in (0.54, 1.32, 1.48):
        x = rv.magnitude_midpoint_shift(dm)
        q = 10 ** (-0.4 * dm)
        # At fraction x of a linear flux ramp the flux is sqrt(q).
        assert 1.0 - (1.0 - q) * x == pytest.approx(math.sqrt(q))
    # Deeper band -> later magnitude midpoint: the sign of the bias.
    assert rv.magnitude_midpoint_shift(1.48) > rv.magnitude_midpoint_shift(
        0.54)


def test_a_magnitude_fit_times_the_deeper_band_later_and_a_flux_fit_does_not():
    """THE ESTIMATOR BIAS, DEMONSTRATED (RF.B2, D3).  One achromatic flux
    ramp, 450 s wide, seen in a shallow band (0.54 mag, like g) and a deep
    one (1.48 mag, like i), noiseless and finely sampled so nothing but the
    estimator can differ.  The magnitude-space ramp puts the deep band's
    edge later by the analytic amount; the flux-space ramp puts both at the
    true instant."""
    w_s = 450.0
    t_true = 0.010
    t = t_true + np.arange(-900.0, 901.0, 15.0) / 86400.0
    tg = t_true + np.arange(-300.0, 300.5, 0.5) / 86400.0
    wg = np.array([w_s]) / 86400.0
    got = {}
    for name, depth in (("shallow", 0.54), ("deep", 1.48)):
        flux = rv.flux_edge_template(t, t_true, w_s / 86400.0, depth)
        mag = 16.0 - 2.5 * np.log10(flux)
        e = np.full(t.size, 0.01)
        fm = rv.fit_edge_profile(t, mag, e, tg, wg, 15.0)
        f, fe, _ = rv.to_relative_flux(mag, e)
        ff = rv.fit_edge_profile(t, f, fe, tg, wg, 15.0)
        got[name] = ((fm["t_edge_d"] - t_true) * 86400.0,
                     (ff["t_edge_d"] - t_true) * 86400.0)
    assert abs(got["shallow"][1]) < 1.0 and abs(got["deep"][1]) < 1.0
    assert got["deep"][0] > got["shallow"][0] > 0.0
    # The differential bias has the sign and roughly the size of the
    # analytic midpoint shift (the least-squares ramp is not exactly the
    # midpoint crossing, hence the generous tolerance).
    want = (rv.magnitude_midpoint_shift(1.48)
            - rv.magnitude_midpoint_shift(0.54)) * w_s
    assert got["deep"][0] - got["shallow"][0] == pytest.approx(want,
                                                              rel=0.6)


def test_flux_edge_time_is_immune_to_constant_dilution():
    """Why flux is the right space: adding a constant, unocculted light
    source (more of it in one band than another) changes the depth in
    magnitudes but not where the flux crosses half-way."""
    t_true = 0.010
    t = t_true + np.arange(-900.0, 901.0, 30.0) / 86400.0
    tg = t_true + np.arange(-200.0, 200.5, 0.5) / 86400.0
    wg = np.array([300.0]) / 86400.0
    flux = rv.flux_edge_template(t, t_true, 300.0 / 86400.0, 1.5)
    out = []
    for dilution in (0.0, 0.5, 2.0):
        mag = 16.0 - 2.5 * np.log10(flux + dilution)
        f, fe, _ = rv.to_relative_flux(mag, np.full(t.size, 0.01))
        out.append(rv.fit_edge_profile(t, f, fe, tg, wg, 30.0)["t_edge_d"])
    assert max(out) - min(out) < 1.0 / 86400.0
    assert out[0] == pytest.approx(t_true, abs=1.0 / 86400.0)


def test_edge_window_selects_and_clips_like_the_published_stage():
    t = 10.0 + np.arange(40) * 219.0 / 86400.0
    y = np.full(40, 16.0)
    y[20] = 21.0                                    # a 5 mag outlier
    e = np.full(40, 0.02)
    t_guess = t[20]
    tt, yy, ee = rv.edge_window(t, y, e, t_guess, PERIOD_D)
    inside = np.abs(t - t_guess) <= p3.EDGE_WINDOW_PHASE * PERIOD_D
    assert tt.size == inside.sum() - 1              # the outlier is gone
    assert 21.0 not in yy
    assert np.all(np.abs(tt - t_guess) <= p3.EDGE_WINDOW_PHASE * PERIOD_D)


def test_quantisation_sigma_is_step_over_root_twelve():
    assert rv.quantisation_sigma_s(3.65) == pytest.approx(3.65 / math.sqrt(12))
    assert math.isnan(rv.quantisation_sigma_s(float("nan")))


# ===========================================================================
# 3.  The per-edge bootstrap
# ===========================================================================
def _boot(noise, seed):
    t, mag, e = _window(noise=noise, seed=seed, width_s=250.0)
    tg = rv.v2_time_grid(0.010, 219.0, PERIOD_D)
    wg = np.array(rv.V2_WIDTH_GRID_S) / 86400.0
    fit = rv.fit_edge_profile(t, mag, np.full(t.size, 0.02), tg, wg, 219.0)

    def refit(ystar):
        return rv.fit_edge_profile(t, ystar, np.full(t.size, 0.02), tg, wg,
                                   219.0)["t_edge_d"]
    return fit, rv.bootstrap_edge(fit, refit, n_boot=150, seed=seed)


def test_bootstrap_error_grows_with_the_misfit():
    _, quiet = _boot(0.02, 1)
    _, loud = _boot(0.15, 1)
    assert quiet["n_ok"] == 150 and loud["n_ok"] == 150
    assert np.isfinite(quiet["sigma_s"]) and loud["sigma_s"] > quiet["sigma_s"]


def test_bootstrap_is_deterministic_and_needs_a_fit():
    _, a = _boot(0.08, 7)
    _, b = _boot(0.08, 7)
    assert a["sigma_s"] == b["sigma_s"]
    empty = rv.bootstrap_edge({"accepted": False}, lambda y: 0.0)
    assert empty["n_ok"] == 0 and math.isnan(empty["sigma_s"])


# ===========================================================================
# 4.  Injection with a known estimator
# ===========================================================================
def _night(seed=0, n=130, cadence_s=219.0, offset_s=0.0):
    t = 100.0 + (np.arange(n) * cadence_s + offset_s) / 86400.0
    return {"t": t, "e": np.full(n, 0.02), "depth_mag": 1.0,
            "pool": np.random.default_rng(seed).normal(0.0, 0.01, 400)}


def test_injection_recovers_a_known_signed_bias_per_band_and_per_pair():
    """Give the injection an 'estimator' with a KNOWN flaw — it reports
    the truth in i and the truth plus 50 s in g — and the injection must
    hand back +50 s for g, 0 for i and +50 s for g-i, signed, with errors
    that are honest about how well it knows them."""
    truth = {}

    def fit_all(band, t, mag, e, t_guess):
        # Recover the true edge as the flux half-crossing of the noiseless
        # part of the template: the midpoint of the faint/bright levels.
        f = 10 ** (-0.4 * (mag - mag.min()))
        half = 0.5 * (f.max() + f.min())
        i = int(np.flatnonzero(f < half)[0])
        frac = (f[i - 1] - half) / (f[i - 1] - f[i])
        te = t[i - 1] + frac * (t[i] - t[i - 1])
        truth[band] = te
        return {"est": te + (50.0 / 86400.0 if band == "g" else 0.0)}

    bands = {"g": _night(1), "i": _night(2, offset_s=146.0)}
    for b in bands.values():
        b["pool"] = np.zeros(400)                  # noiseless: exact answer
    out = rv.inject_recover(bands, PERIOD_D, (480.0,), fit_all, ("est",),
                            n_real=60, seed=5)
    by = {r["band"]: r for r in out["per_band"]}
    assert by["g"]["bias_mean_s"] == pytest.approx(50.0, abs=3.0)
    assert by["i"]["bias_mean_s"] == pytest.approx(0.0, abs=3.0)
    pair = [r for r in out["per_pair"] if r["band_a"] == "g"
            and r["band_b"] == "i"][0]
    assert pair["dbias_mean_s"] == pytest.approx(50.0, abs=4.0)
    assert pair["dbias_mean_s"] > 0                # SIGNED, never |bias|
    assert pair["n_ok"] > 40 and pair["dbias_se_s"] < 3.0


def test_injection_is_reproducible_and_counts_its_misses():
    def never(band, t, mag, e, t_guess):
        return {"est": float("nan")}

    def always(band, t, mag, e, t_guess):
        return {"est": t_guess}

    bands = {"g": _night(3), "i": _night(4)}
    miss = rv.inject_recover(bands, PERIOD_D, (60.0, 240.0), never, ("est",),
                             n_real=20)
    assert all(r["n_ok"] == 0 and r["bias_mean_s"] is None
               for r in miss["per_band"])
    assert all(r["n_try"] == 20 for r in miss["per_band"])
    a = rv.inject_recover(bands, PERIOD_D, (240.0,), always, ("est",),
                          n_real=40, seed=9)
    b = rv.inject_recover(bands, PERIOD_D, (240.0,), always, ("est",),
                          n_real=40, seed=9)
    assert a == b
    # An estimator that just returns its starting guess has the guess's
    # scatter and no bias: the guess scatter is really being applied.
    g = [r for r in a["per_band"] if r["band"] == "g"][0]
    assert g["rms_s"] == pytest.approx(rv.GUESS_SCATTER_S, rel=0.35)
    with pytest.raises(ValueError):
        rv.inject_recover(bands, PERIOD_D, (240.0,), always, ("est",),
                          noise="pink")


def test_per_pair_rows_cover_unequal_widths():
    """'Per-band ramp widths' (CV-R2): every combination of width in band
    a and width in band b is tabulated, not only the diagonal."""
    def always(band, t, mag, e, t_guess):
        return {"est": t_guess}
    out = rv.inject_recover({"g": _night(1), "r": _night(2), "i": _night(3)},
                            PERIOD_D, (60.0, 240.0, 480.0), always, ("est",),
                            n_real=15)
    assert len(out["per_band"]) == 3 * 3
    assert len(out["per_pair"]) == len(rv.BAND_PAIRS) * 3 * 3
    assert {(r["true_width_a_s"], r["true_width_b_s"])
            for r in out["per_pair"]} == {(a, b) for a in (60.0, 240.0, 480.0)
                                          for b in (60.0, 240.0, 480.0)}


def test_residual_pool_is_the_stars_own_noise_without_the_edges():
    """A periodic top-hat light curve plus 3% flickering: the pool must
    return residuals of about 3% and must have dropped the points next to
    the two edges, where a binned profile cannot follow the star."""
    rng = np.random.default_rng(6)
    t = 50.0 + np.arange(400) * 60.0 / 86400.0
    ph = p3.phase_of(t, PERIOD_D, 50.0)
    flux = np.where(ph < 0.3, 1.0, 0.4) * (1.0 + rng.normal(0, 0.03, t.size))
    pool = rv.residual_pool(t, -2.5 * np.log10(flux), PERIOD_D, 50.0)
    assert 0 < pool.size < t.size                   # edge points dropped
    assert np.std(pool) == pytest.approx(0.03 * 0.6, rel=0.6)
    assert rv.residual_pool(t[:30], -2.5 * np.log10(flux[:30]), PERIOD_D,
                            50.0).size == 0         # too sparse: says so


# ===========================================================================
# 5.  The O-C model
# ===========================================================================
def _oc_data(gamma=0.0, era_step=0.0, noise=20.0, seed=0):
    rng = np.random.default_rng(seed)
    cyc = np.concatenate([np.linspace(13100, 13600, 9),
                          np.linspace(18080, 19050, 24),
                          np.linspace(21860, 21900, 3)])
    band = np.array((["g", "r", "i"] * 12))
    era = np.where(cyc < 15000, "7", "76")
    const = {"g": -60.0, "r": -10.0, "i": 40.0}
    x = cyc - cyc.mean()
    oc = (np.array([const[b] for b in band]) + 0.004 * x + gamma * x * x
          + np.where(era == "76", era_step, 0.0)
          + rng.normal(0.0, noise, cyc.size))
    return cyc, oc, np.full(cyc.size, noise), band, era, const


def test_oc_fit_recovers_band_constants_and_the_quadratic():
    gamma = 3.0e-6
    cyc, oc, sig, band, era, const = _oc_data(gamma=gamma, noise=15.0)
    fit = rv.fit_oc_model(cyc, oc, sig, band=band)
    assert fit["n"] == 36 and fit["n_params"] == 5 and fit["dof"] == 31
    for b in "gri":
        i = fit["names"].index(f"band:{b}")
        assert fit["coef"][i] == pytest.approx(
            const[b], abs=4 * math.sqrt(fit["cov"][i, i]))
    assert fit["gamma"] == pytest.approx(gamma,
                                         abs=4 * fit["gamma_sigma_budget"])
    assert abs(fit["gamma"]) > 5 * fit["gamma_sigma_budget"]
    assert 0.4 < fit["chi2nu"] < 2.0
    assert fit["beta"] == pytest.approx(0.004, abs=5 * fit["beta_sigma_budget"])


def test_pooling_bands_leaks_a_band_offset_into_the_quadratic():
    """DS.F2: when the band mix changes with time a band offset aliases
    into the curvature.  Build that case — no true curvature, a 100 s band
    offset, g-only at one end and i-only at the other — and the pooled fit
    'detects' a quadratic that the per-band fit correctly does not."""
    cyc = np.linspace(-5000, 5000, 60)
    band = np.where(np.abs(cyc) > 3000, "i", "g")
    oc = np.where(band == "i", 60.0, -60.0) + np.random.default_rng(3).normal(
        0.0, 10.0, cyc.size)
    sig = np.full(cyc.size, 10.0)
    pooled = rv.fit_oc_model(cyc, oc, sig)
    split = rv.fit_oc_model(cyc, oc, sig, band=band)
    assert abs(pooled["gamma"]) > 5 * pooled["gamma_sigma_scatter"]
    assert abs(split["gamma"]) < 3 * split["gamma_sigma_scatter"]
    assert pooled["chi2nu"] > 5 * split["chi2nu"]


def test_era_offset_is_fitted_and_costs_the_quadratic_its_lever_arm():
    cyc, oc, sig, band, era, _ = _oc_data(era_step=150.0, noise=15.0, seed=2)
    no_era = rv.fit_oc_model(cyc, oc, sig, band=band)
    with_era = rv.fit_oc_model(cyc, oc, sig, band=band, era=era)
    i = with_era["names"].index("era:76")
    assert with_era["coef"][i] == pytest.approx(
        150.0, abs=4 * math.sqrt(with_era["cov"][i, i]))
    # Without the nuisance term the step is forced into the polynomial.
    assert no_era["chi2nu"] > with_era["chi2nu"]
    # With it, the error on the quadratic is larger: honesty has a price.
    assert with_era["gamma_sigma_budget"] > 1.5 * no_era["gamma_sigma_budget"]


def test_scatter_and_budget_errors_differ_by_root_chi2nu_in_both_directions():
    cyc, oc, sig, band, era, _ = _oc_data(noise=10.0, seed=5)
    over = rv.fit_oc_model(cyc, oc, 4.0 * sig, band=band)   # errors 4x big
    under = rv.fit_oc_model(cyc, oc, 0.25 * sig, band=band)  # 4x small
    for f in (over, under):
        assert f["gamma_sigma_scatter"] == pytest.approx(
            f["gamma_sigma_budget"] * math.sqrt(f["chi2nu"]))
    assert over["chi2nu"] < 0.2 and under["chi2nu"] > 5.0
    assert over["gamma_sigma_scatter"] < over["gamma_sigma_budget"]
    assert under["gamma_sigma_scatter"] > under["gamma_sigma_budget"]
    # The scatter-based error does not care what the budget claimed.
    assert over["gamma_sigma_scatter"] == pytest.approx(
        under["gamma_sigma_scatter"], rel=1e-6)


def test_oc_fit_is_well_conditioned_in_raw_cycles_and_names_singletons():
    """Cycle numbers of 2e4 squared are 4e8; the fit must not call a
    well-posed problem singular (it once did), and a band with one epoch
    must be reported as absorbed by its own constant."""
    cyc, oc, sig, band, era, _ = _oc_data(seed=8)
    cyc = np.append(cyc, 13700.0)
    oc = np.append(oc, 55.0)
    sig = np.append(sig, 20.0)
    band = np.append(band, "z")
    era = np.append(era, "47")
    fit = rv.fit_oc_model(cyc, oc, sig, band=band, era=era)
    assert not fit["rank_deficient"]
    assert "band:z" in fit["singleton_levels"]
    assert any(s.startswith("era:47") for s in fit["singleton_levels"])
    assert fit["n_informative"] == fit["n"] - 1
    iz = list(band).index("z")
    assert abs(fit["resid"][iz]) < 1e-6           # fitted exactly: no info
    tiny = rv.fit_oc_model([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], [1.0] * 3)
    assert tiny["rank_deficient"] and tiny["coef"] is None


def test_pdot_from_quadratic_is_the_published_conversion():
    # 0.5 * Pdot * P * E^2 seconds  <=>  gamma = 0.5 * Pdot * P[s].
    pdot = 1.14e-9
    gamma = 0.5 * pdot * PERIOD_D * 86400.0
    assert rv.pdot_from_quadratic(gamma, PERIOD_D) == pytest.approx(pdot)


def test_night_level_epochs_combine_bands_after_removing_their_constants():
    rows = [
        {"night": "n1", "band": "g", "cycle": 10.0, "oc_s": -70.0,
         "sigma_s": 20.0, "n_cycles": 2, "era": 76},
        {"night": "n1", "band": "i", "cycle": 12.0, "oc_s": 50.0,
         "sigma_s": 10.0, "n_cycles": 3, "era": 76},
        {"night": "n2", "band": "g", "cycle": 30.0, "oc_s": -40.0,
         "sigma_s": 20.0, "n_cycles": 1, "era": 76},
    ]
    out = rv.night_level_epochs(rows, {"g": -60.0, "i": 40.0})
    assert [r["night"] for r in out] == ["n1", "n2"]
    n1, n2 = out
    assert n1["n_bands"] == 2 and n1["n_cycles"] == 5 and n1["bands"] == "gi"
    w = np.array([1 / 400.0, 1 / 100.0])
    assert n1["oc_s"] == pytest.approx(np.sum(w * [-10.0, 10.0]) / w.sum())
    assert n1["sigma_s"] == pytest.approx(math.sqrt(1 / w.sum()))
    assert n1["within_rms_s"] == pytest.approx(np.std([-10.0, 10.0], ddof=1))
    assert n2["oc_s"] == pytest.approx(20.0) and n2["within_rms_s"] is None


def test_per_group_chi2_exposes_two_miscalibrations_behind_one_average():
    """Standing rule 1.  One band with errors 2x too small, another 2x too
    large: together chi2 per epoch is about 2, which hides that one is at
    4 and the other at 0.25."""
    rng = np.random.default_rng(10)
    resid = np.concatenate([rng.normal(0, 100, 400), rng.normal(0, 50, 400)])
    sig = np.concatenate([np.full(400, 50.0), np.full(400, 100.0)])
    grp = ["g"] * 400 + ["i"] * 400
    rows = {r["group"]: r for r in rv.per_group_chi2(resid, sig, grp)}
    assert rows["g"]["chi2_per_n"] == pytest.approx(4.0, rel=0.2)
    assert rows["i"]["chi2_per_n"] == pytest.approx(0.25, rel=0.2)
    assert rows["g"]["n"] == 400 and rows["i"]["sigma_median_s"] == 100.0


def test_scatter_sigma_by_band_recovers_the_single_cycle_scatter():
    rng = np.random.default_rng(11)
    n = rng.integers(1, 5, 600).astype(float)
    band = np.array(["g"] * 300 + ["i"] * 300)
    s_true = np.where(band == "g", 120.0, 45.0)
    resid = rng.normal(0.0, s_true / np.sqrt(n))
    got = rv.scatter_sigma_by_band(resid, band, n, dof_total=600)
    assert got["g"] == pytest.approx(120.0, rel=0.12)
    assert got["i"] == pytest.approx(45.0, rel=0.12)
    assert "z" not in rv.scatter_sigma_by_band([1.0, 2.0, 3.0],
                                               ["g", "g", "z"], [1, 1, 1], 2)


# ===========================================================================
# 6.  Longitude and physical scales
# ===========================================================================
def test_seconds_to_degrees_is_360_dt_over_p():
    p_s = PERIOD_D * 86400.0
    assert float(rv.seconds_to_degrees(p_s, PERIOD_D)) == pytest.approx(360.0)
    # The physicist's check: 84 s is 4.4 degrees.
    assert float(rv.seconds_to_degrees(84.0, PERIOD_D)) == pytest.approx(
        4.4, abs=0.05)


def test_gr_scale_is_orders_below_any_bound_this_paper_can_set():
    pd = rv.gr_orbital_pdot(PERIOD_D)
    assert 1e-15 < pd < 1e-12
    # Scaling: Pdot goes as P^(-5/3) and as chirp mass^(5/3).
    assert rv.gr_orbital_pdot(2 * PERIOD_D) == pytest.approx(
        pd * 2 ** (-5.0 / 3.0))
    assert rv.gr_orbital_pdot(PERIOD_D, 1.5, 0.34) == pytest.approx(
        pd * 2 ** (5.0 / 3.0))


def test_cyclotron_harmonic_rises_to_the_blue():
    g = rv.cyclotron_harmonic(4810.0, 12.1)
    i = rv.cyclotron_harmonic(7520.0, 12.1)
    assert g > i > 5
    assert rv.cyclotron_harmonic(10710.0, 100.0) == pytest.approx(1.0)


def test_literature_scales_agree_with_the_literature_package():
    """The predicted scales are arithmetic on constants the cv-literature
    package verified (``committee/work/cv-literature/literature_scales.md``).
    Both packages must print the same numbers, so the values that file
    records are pinned here."""
    assert rv.historical_edge_shift_s(PERIOD_D) == pytest.approx(170.83,
                                                                 abs=0.01)
    assert float(rv.seconds_to_degrees(
        rv.historical_edge_shift_s(PERIOD_D), PERIOD_D)) == pytest.approx(9.0)
    # Redder band ends later: the bluer-minus-redder prediction is negative.
    gi = rv.predicted_chromatic_shift_s(PERIOD_D, "g", "i")
    assert gi == pytest.approx(-37.19, abs=0.05)
    assert rv.predicted_chromatic_shift_s(PERIOD_D, "g", "r") < 0
    assert rv.predicted_chromatic_shift_s(PERIOD_D, "r", "i") < 0
    assert gi == pytest.approx(
        rv.predicted_chromatic_shift_s(PERIOD_D, "g", "r")
        + rv.predicted_chromatic_shift_s(PERIOD_D, "r", "i"))
    assert rv.libration_pdot(PERIOD_D) == pytest.approx(3.5707e-11, rel=1e-3)
    assert rv.secular_pdot(PERIOD_D, 1.0) == pytest.approx(2.1653e-13,
                                                           rel=1e-3)
    assert rv.secular_pdot(PERIOD_D, 5.0) == pytest.approx(4.3307e-14,
                                                           rel=1e-3)
    # Cycle-count drift under Cropper's published sigma_P: 0.022 cycle
    # over the paper's 21,869 cycles, a margin of 23 on half a cycle.
    drift = rv.phase_drift_cycles(21869 * PERIOD_D, PERIOD_D,
                                  rv.CROPPER_PERIOD_SIGMA_D)
    assert drift == pytest.approx(0.0221, abs=1e-4)
    assert 0.5 / drift == pytest.approx(22.6, abs=0.1)
    # YZ Cnc: 0.027 cycle per day under the VSX period's published error.
    assert rv.phase_drift_cycles(1.0, 0.0868, rv.YZCNC_SH_PERIOD_SIGMA_D
                                 ) == pytest.approx(0.0265, abs=1e-4)
    assert round(rv.cyclotron_harmonic(rv.BAND_WAVELENGTH_A["g"],
                                       rv.STLMI_FIELD_MG)) == 18
    assert round(rv.cyclotron_harmonic(rv.BAND_WAVELENGTH_A["i"],
                                       rv.STLMI_FIELD_MG)) == 12


def test_slope_permutation_finds_a_dependence_and_no_more():
    rng = np.random.default_rng(13)
    x = rng.uniform(-1, 1, 40)
    slope, icpt, p = rv.slope_permutation_p(x, 5.0 * x + rng.normal(0, 1, 40),
                                            n_perm=2000)
    assert slope == pytest.approx(5.0, abs=1.0) and p < 0.01
    _, _, p0 = rv.slope_permutation_p(x, rng.normal(0, 1, 40), n_perm=2000)
    assert p0 > 0.05
    assert math.isnan(rv.slope_permutation_p([1, 1, 1, 1], [1, 2, 3, 4])[0])


# ===========================================================================
# 7.  The colour curve
# ===========================================================================
def _two_band_night(night, rng, amp=0.75, lag_s=73.0, cadence_s=219.0,
                    n=160, noise=0.01):
    """Bands a and b on one night: b is exposed ``lag_s`` after a.  Band a
    varies with phase (a steep-sided bright phase); band b is constant, so
    the true colour a-b is band a's light curve."""
    t0 = 200.0 + night
    ta = t0 + np.arange(n) * cadence_s / 86400.0
    tb = ta + lag_s / 86400.0

    def curve(t):
        ph = p3.phase_of(t, PERIOD_D, 200.0)
        return amp * np.clip((0.15 - np.abs(ph - 0.5)) / 0.05, 0.0, 1.0)
    return (ta, curve(ta) + rng.normal(0, noise, n),
            tb, rng.normal(0, noise, n), curve)


def test_interpolated_colour_has_no_first_order_timing_bias():
    """RF.M1.  On a linear slope a fixed pairing lag biases the colour by
    slope x lag; interpolation to the same instant does not."""
    ta = np.arange(0.0, 100.0)
    tb = ta[:-1] + 0.33                 # always a third of a step later
    slope = 0.02
    ma = slope * ta
    mb = np.zeros(tb.size)
    t, col, gap = rv.interpolated_colour(ta, ma, tb, mb, max_gap_s=2 * 86400)
    assert np.allclose(col, slope * tb)          # the colour AT tb
    assert np.allclose(gap, 86400.0)
    # Nearest-neighbour pairing is off by slope * 0.33 everywhere.
    near = ma[:-1] - mb
    assert np.allclose(near - slope * tb, -slope * 0.33)
    # Nothing is interpolated across a gap wider than the gate.
    ta2 = np.array([0.0, 1.0, 50.0, 51.0])
    t2, c2, g2 = rv.interpolated_colour(ta2, ta2 * 0, np.array([0.5, 25.0,
                                                               50.5]),
                                        np.zeros(3), max_gap_s=2 * 86400)
    assert t2.size == 2 and 25.0 not in t2


def test_colour_summary_recovers_amplitude_and_phases_with_night_errors():
    rng = np.random.default_rng(14)
    ph_all, col_all, night_all = [], [], []
    for k in range(6):
        ta, ma, tb, mb, _ = _two_band_night(k, rng)
        t, col, _ = rv.interpolated_colour(ta, ma, tb, mb, 2 * 219.0)
        ph_all.append(p3.phase_of(t, PERIOD_D, 200.0))
        col_all.append(col)
        night_all.append(np.full(t.size, f"n{k}"))
    s = rv.colour_curve_summary(np.concatenate(ph_all),
                                np.concatenate(col_all),
                                np.concatenate(night_all), n_boot=120)
    assert s["n_nights"] == 6
    assert s["amplitude"] == pytest.approx(0.75, abs=0.06)
    assert s["amplitude_err"] is not None and s["amplitude_err"] < 0.08
    # Red (large) at the centre of the bump, blue well away from it.
    assert abs(rv.circ_diff(s["phase_red"], 0.5)) < 0.08
    assert abs(rv.circ_diff(s["phase_blue"], 0.5)) > 0.2
    assert np.isfinite(s["bin_err"]).sum() >= 10


def test_colour_summary_refuses_an_error_bar_from_two_nights():
    rng = np.random.default_rng(15)
    ph = rng.uniform(0, 1, 400)
    col = 0.5 * np.sin(2 * np.pi * ph) + rng.normal(0, 0.01, 400)
    s = rv.colour_curve_summary(ph, col, ["a"] * 200 + ["b"] * 200,
                                n_boot=50)
    assert s["amplitude"] == pytest.approx(1.0, abs=0.12)
    assert s["amplitude_err"] is None and s["n_boot_ok"] == 0
    thin = rv.colour_curve_summary(ph[:20], col[:20], ["a"] * 20, n_boot=10)
    assert thin["amplitude"] is None              # under-filled bins: no curve


def test_binned_curve_leaves_underfilled_bins_empty():
    c, m, n = rv.binned_curve([0.01, 0.02, 0.03, 0.51], [1.0, 2.0, 3.0, 9.0],
                              n_bins=4, min_count=2)
    assert n.tolist() == [3, 0, 1, 0]
    assert m[0] == 2.0 and np.isnan(m[2])


def test_curve_repeatability_is_about_shape_not_zero_point():
    x = np.linspace(0, 1, 20, endpoint=False)
    a = 0.4 * np.sin(2 * np.pi * x)
    err = np.full(20, 0.02)
    same = rv.curve_repeatability(a, err, a + 0.3, err)    # shifted: same
    assert same["rms_diff"] == pytest.approx(0.0, abs=1e-12)
    assert same["r"] == pytest.approx(1.0) and same["scale"] == pytest.approx(1.0)
    assert same["chi2nu"] == pytest.approx(0.0, abs=1e-9) and same["dof"] == 19
    half = rv.curve_repeatability(a, err, 0.5 * a, err)
    assert half["scale"] == pytest.approx(0.5) and half["chi2nu"] > 10
    other = rv.curve_repeatability(a, err, 0.4 * np.cos(2 * np.pi * x), err)
    assert abs(other["r"]) < 0.1
    few = rv.curve_repeatability(a[:4], err[:4], a[:4], err[:4])
    assert few["n_bins"] == 4 and math.isnan(few["rms_diff"])


# ===========================================================================
# 8.  What a recovery contour excludes
# ===========================================================================
def test_a_contour_above_the_floor_excludes_nothing_at_the_floor():
    """PH.P1 / RF.B4.  The first draft read 'the contour sits above the
    50 mmag floor' as 'a superhump would have been seen'.  It is the
    reverse: with contours of 82-237 mmag, a 50 mmag signal is excluded in
    NO run, a 150 mmag one in some, a 250 mmag one in all that have a
    contour — and the runs without a contour exclude nothing at all."""
    contours = [0.130, None, None, None, None, 0.100, None, None, 0.0815,
                0.1266, None, None, None, 0.2368, None, None, None, None]
    s = rv.excluded_amplitude_summary(contours, (0.05, 0.10, 0.15, 0.25))
    assert s["n_run_filters"] == 18 and s["n_with_contour"] == 5
    assert s["contour_min"] == pytest.approx(0.0815)
    assert s["contour_max"] == pytest.approx(0.2368)
    assert s["n_excluding"][0.05] == 0
    assert s["n_excluding"][0.10] == 2
    assert s["n_excluding"][0.15] == 4
    assert s["n_excluding"][0.25] == 5
    # Monotone in the threshold, and never more than the runs with a contour.
    counts = [s["n_excluding"][k] for k in sorted(s["n_excluding"])]
    assert counts == sorted(counts) and counts[-1] <= s["n_with_contour"]
    empty = rv.excluded_amplitude_summary([None, None], (0.05,))
    assert empty["n_with_contour"] == 0 and empty["contour_min"] is None
    assert empty["n_excluding"][0.05] == 0


def test_half_level_crossing_finds_a_known_shift_without_a_model():
    rng = np.random.default_rng(16)
    dt = rng.uniform(-1100, 1100, 1500)

    def profile(x, shift):
        return 1.0 - 0.7 * np.clip((x - shift) / 300.0 + 0.5, 0.0, 1.0)
    early = rv.half_level_crossing(dt, profile(dt, -90.0)
                                   + rng.normal(0, 0.02, dt.size))
    late = rv.half_level_crossing(dt, profile(dt, 0.0)
                                  + rng.normal(0, 0.02, dt.size))
    assert early["t_half_s"] == pytest.approx(-90.0, abs=20.0)
    assert late["t_half_s"] == pytest.approx(0.0, abs=20.0)
    assert early["t_half_s"] - late["t_half_s"] == pytest.approx(-90.0,
                                                               abs=25.0)
    assert early["depth_frac"] == pytest.approx(0.7, abs=0.05)
    flat = rv.half_level_crossing(dt, np.ones(dt.size))
    assert flat["t_half_s"] is None               # no edge: says so


# ===========================================================================
# 9.  The built products (skipped when absent)
# ===========================================================================
@pytest.fixture(scope="module")
def phot():
    if not PHOT_DB.exists():
        pytest.skip("cv_timeseries.sqlite not built in this checkout")
    con = sqlite3.connect(f"file:{PHOT_DB}?mode=ro", uri=True, timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    if "rv_result" not in have or "rv_edge" not in have:
        con.close()
        pytest.skip("revision stage not built in this checkout")
    yield con
    con.close()


class TestRevisionProducts:
    """Invariants of the ``rv_`` tables.  Each names what would be wrong
    with the report if it failed."""

    def test_the_revision_re_measured_exactly_the_published_edges(self, phot):
        """The revision changes the ESTIMATOR, never the sample.  Every
        published accepted ST LMi edge has a v1 row whose epoch is the
        published one, and the refit on the reconstructed window reproduces
        it — which is the proof the reconstruction is right."""
        n_pub = phot.execute("SELECT count(*) FROM p3_edge WHERE "
                             "target_key='stlmi' AND accepted=1").fetchone()[0]
        n_rv, worst = phot.execute(
            "SELECT count(*), max(abs(v1_refit_minus_stored_s)) FROM "
            "rv_edge WHERE estimator='v1'").fetchone()
        assert n_rv == n_pub
        assert worst < 1e-3, (
            f"the revision's v1 refit differs from the published epoch by "
            f"{worst} s: it is not fitting the published windows")
        bad = phot.execute("""
            SELECT count(*) FROM rv_edge e JOIN p3_edge p
              ON p.series_key = e.series_key AND p.cycle = e.cycle
            WHERE e.estimator='v1' AND abs(e.t_edge_bjd - p.t_edge_bjd)
                  > 1e-9""").fetchone()[0]
        assert bad == 0

    def test_every_edge_has_its_own_error(self, phot):
        """RF.B3 / DS.F3: 'a Monte-Carlo error on every edge' was three
        constants.  The bootstrap's must be per edge."""
        n, n_distinct, n_null = phot.execute(
            "SELECT count(*), count(DISTINCT round(sigma_edge_s, 3)), "
            "sum(sigma_edge_s IS NULL) FROM rv_edge WHERE estimator='v1'"
        ).fetchone()
        assert n_null == 0
        # Revised 2026-10-04.  The guard is "not one constant per series".
        # The first form, "> 80 per cent distinct", fails at 50 of 63
        # because the bootstrap error is a percentile range of a time GRID
        # (3.7 s in 2025, up to 8 s in 2024), so edges of one series can
        # share it to the grid step.  The equivalent invariant: no series
        # with three or more edges has a single error, and every tie is
        # between edges on the same grid step.
        for sk, n_e, n_d in phot.execute(
                "SELECT series_key, count(*), count(DISTINCT "
                "round(sigma_edge_s, 3)) FROM rv_edge WHERE estimator='v1' "
                "GROUP BY series_key").fetchall():
            if n_e >= 3:
                assert n_d > 1, f"{sk}: one error for all {n_e} edges"
        for v, n_steps in phot.execute(
                "SELECT round(sigma_edge_s, 3), count(DISTINCT "
                "round(grid_step_s, 3)) FROM rv_edge WHERE estimator='v1' "
                "GROUP BY 1 HAVING count(*) > 1").fetchall():
            assert n_steps == 1, (
                f"edges on different grids share the error {v}: not "
                f"grid quantisation")
        assert n_distinct > 0.5 * n
        q = phot.execute("SELECT min(sigma_quant_s), max(sigma_quant_s) "
                         "FROM rv_edge WHERE estimator='v1'").fetchone()
        assert q[0] > 0, "grid quantisation is missing from the budget"

    def test_injection_ran_on_every_night_that_carries_an_epoch(self, phot):
        """RF.B2: no transported budget.  Every (night, band) with a
        published g/r/i epoch has its own injection rows, including the
        2024 High Gain nights."""
        want = {(r[0], r[1].lower()) for r in phot.execute(
            "SELECT night, filter FROM p3_oc_night WHERE target_key='stlmi' "
            "AND lower(filter) IN ('g','r','i')")}
        have = {(r[0], r[1]) for r in phot.execute(
            "SELECT DISTINCT night, band FROM rv_inject_band WHERE "
            "noise='rolled' AND estimator='v1'")}
        assert want <= have, f"no injection for {sorted(want - have)}"
        assert phot.execute("SELECT count(DISTINCT night) FROM "
                            "rv_inject_band WHERE era_id=7").fetchone()[0] > 0

    def test_injection_biases_are_signed(self, phot):
        """A table of biases that are all non-negative is a table of
        absolute values.  The grid must contain both signs."""
        lo, hi = phot.execute("SELECT min(bias_mean_s), max(bias_mean_s) "
                              "FROM rv_inject_band").fetchone()
        assert lo < 0 < hi

    def test_every_oc_fit_reports_chi2_with_its_dof(self, phot):
        """Standing rule 1."""
        rows = phot.execute("SELECT variant, model, chi2, dof, chi2nu, "
                            "n_epochs, n_params FROM rv_oc_fit WHERE chi2 "
                            "IS NOT NULL").fetchall()
        assert rows
        for v, m, chi2, dof, chi2nu, n, p in rows:
            assert dof == n - p, f"{v} {m}"
            assert chi2nu == pytest.approx(chi2 / dof), f"{v} {m}"

    def test_the_first_drafts_bound_is_reproduced_before_it_is_revised(
            self, phot):
        """The revision's 'one constant, budget errors' fit must return the
        number the first draft published, or the comparison of every other
        model against it compares against nothing."""
        mine = phot.execute("SELECT pdot, pdot_limit3_budget FROM rv_oc_fit "
                            "WHERE variant='v1/pub' AND model='pooled'"
                            ).fetchone()
        pub = phot.execute("SELECT pdot, pdot_limit3 FROM p3_cycle_count "
                           "WHERE target_key='stlmi'").fetchone()
        assert mine[0] == pytest.approx(pub[0], rel=1e-6)
        assert mine[1] == pytest.approx(pub[1], rel=1e-6)

    def test_night_level_epochs_are_one_per_night(self, phot):
        n, n_distinct = phot.execute(
            "SELECT count(*), count(DISTINCT night) FROM rv_oc_epoch WHERE "
            "variant='v1/scatter'").fetchone()
        n_nights = phot.execute(
            "SELECT count(DISTINCT night) FROM p3_oc_night WHERE "
            "target_key='stlmi' AND lower(filter) IN ('g','r','i')"
        ).fetchone()[0]
        assert n == n_distinct == n_nights

    def test_every_scalar_is_explained_and_typed(self, phot):
        rows = phot.execute("SELECT key, fmt, origin, note FROM rv_result"
                            ).fetchall()
        assert len(rows) > 100
        for key, fmt, origin, note in rows:
            assert note and len(note) > 15, key
            assert origin in ("measured", "constant", "literature"), key
            assert fmt in ("int", "f0", "f1", "f2", "f3", "f4", "f6", "f8",
                           "sci1",
                           "sci2", "p", "pct", "text"), key

    def test_external_scales_are_not_recorded_as_measurements(self, phot):
        for key in ("rv pdot gr scale", "rv pdot async scale",
                    "rv off predicted gi s", "rv long historical shift deg",
                    "rv sh peak semi mmag", "rv yz vp period d"):
            r = phot.execute("SELECT origin FROM rv_result WHERE key=?",
                             (key,)).fetchone()
            assert r is not None and r[0] == "literature", key

    def test_revision_macros_are_emitted_and_letters_only(self, phot):
        from macro_phot import numbers_cv as nx
        phot.row_factory = None
        nums = nx.revision_numbers(phot)
        names = [n.macro for n in nums]
        assert len(names) == len(set(names)), "duplicate revision macros"
        assert all(nm.isalpha() and nm.startswith("NumRv") for nm in names)
        assert all(n.value is not None for n in nums), (
            "a revision macro has no value: " + ", ".join(
                n.macro for n in nums if n.value is None))
        kinds = {n.source.split(" ")[0] for n in nums}
        assert kinds <= {"rv_result", "CV-R", "literature:"}
