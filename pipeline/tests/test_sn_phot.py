"""Tests for macro_sn.snphot — the decision rules of the SN 2023ixf release.

Synthetic-recovery tests: every estimator is shown to return what was put
in, so a number in the paper is never the product of an untested rule.
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from macro_sn import snphot as sp  # noqa: E402


def test_wcs_acceptance_fails_closed():
    assert sp.wcs_acceptable(50, 0.3, 0.54)
    assert not sp.wcs_acceptable(sp.WCS_MIN_MATCH - 1, 0.3, 0.54)
    assert not sp.wcs_acceptable(50, 2.0, 0.54)
    assert not sp.wcs_acceptable(None, 0.3, 0.54)
    assert not sp.wcs_acceptable(50, float("nan"), 0.54)
    assert sp.shift_acceptable(sp.WCS_MIN_MATCH_SHIFT, 0.3, 0.54)
    assert not sp.shift_acceptable(sp.WCS_MIN_MATCH_SHIFT - 1, 0.3, 0.54)


def test_check_star_split_is_deterministic_and_about_one_in_four():
    rng = np.random.default_rng(1)
    ra, dec = rng.uniform(210, 211, 4000), rng.uniform(54, 55, 4000)
    a = np.array([sp.is_check_star(r, d) for r, d in zip(ra, dec)])
    b = np.array([sp.is_check_star(r, d) for r, d in zip(ra, dec)])
    assert (a == b).all()
    assert abs(a.mean() - 1 / sp.CHECK_MODULUS) < 0.03


def test_young_scintillation_scaling():
    s1 = sp.young_scintillation_mag(1.0, 1.0)
    assert sp.young_scintillation_mag(4.0, 1.0) == pytest.approx(s1 / 2)
    assert sp.young_scintillation_mag(1.0, 2.0) == pytest.approx(s1 * 2 ** 1.75)
    # 0.5 m, 1515 m, X = 1.2, 0.5 s: a few mmag (strategy §5: ~5-10 mmag)
    assert 0.004 < sp.young_scintillation_mag(0.5, 1.2) < 0.012


def test_chi2nu_and_verdict():
    rng = np.random.default_rng(2)
    g = np.repeat(np.arange(200), 10)
    v = rng.normal(0, 0.01, g.size)
    c2, dof = sp.chi2nu(v, np.full(g.size, 0.01), g)
    assert dof == 200 * 9
    assert c2 == pytest.approx(1.0, abs=0.1)
    assert sp.chi2_verdict(0.4).startswith("DEFECT")
    assert sp.chi2_verdict(2.5).startswith("DEFECT")
    assert sp.chi2_verdict(1.1) == "pass"


def test_fit_floor_recovers_injected_floor():
    rng = np.random.default_rng(3)
    n = 20000
    sph = rng.uniform(0.005, 0.03, n)
    ssc = np.full(n, 0.004)
    r = rng.normal(0, np.sqrt(sph ** 2 + ssc ** 2 + 0.015 ** 2))
    assert sp.fit_floor(r, sph, ssc) == pytest.approx(0.015, abs=0.001)


def test_pooled_scintillation_recovers_C_and_floors():
    rng = np.random.default_rng(4)
    n = 30000
    t = rng.choice([0.5, 2.0, 8.0, 32.0], n)
    su = sp.young_scintillation_mag(t, 1.2)
    band = rng.integers(0, 3, n)
    floors = np.array([0.010, 0.015, 0.020])
    sph = np.full(n, 0.003)
    r = rng.normal(0, np.sqrt(sph ** 2 + (2.0 * su) ** 2 + floors[band] ** 2))
    C, fl = sp.fit_scint_pooled(r, sph, su, band)
    assert C == pytest.approx(2.0, rel=0.1)
    for k in range(3):
        assert fl[k] == pytest.approx(floors[k], abs=0.0015)


def test_solve_band_recovers_colour_term_zero_points_and_k2():
    rng = np.random.default_rng(5)
    S, F = 120, 60
    col = rng.uniform(0.3, 2.0, S)
    cat = rng.uniform(12, 16, S)
    X = rng.uniform(1.0, 2.0, F)
    zp = rng.normal(0, 0.3, F)
    cterm, k2, Z = 0.10, -0.03, 1.7
    sig = np.full((S, F), 0.01)
    mag = (cat[:, None] + Z + cterm * (col[:, None] - sp.COLOUR_PIVOT) + zp[None, :]
           + k2 * (col[:, None] - sp.COLOUR_PIVOT) * (X[None, :] - 1)
           + rng.normal(0, 0.01, (S, F)))
    is_check = np.zeros(S, bool)
    is_check[::4] = True
    sol = sp.solve_band(mag, sig, np.full(F, 0.001), np.ones(S, bool), is_check,
                        col, cat, X, fixed_C=1.0)
    assert sol["cterm"] == pytest.approx(cterm, abs=0.005)
    assert sol["k2"] == pytest.approx(k2, abs=0.005)
    # zero points are recovered up to the gauge constant
    d = (sol["zp"] + sol["Z"]) - (zp + Z)
    assert np.std(d) < 0.003
    c2, dof, *_ = sp.check_star_test(sol["mcorr"], sol["sig"], sol["zp"],
                                     sol["zp_err"], is_check)
    assert 0.8 < c2 < 1.2 and dof > 1000


def test_delta_flat_round_trip():
    rng = np.random.default_rng(6)
    x, y = rng.uniform(0, 4096, 5000), rng.uniform(0, 4096, 5000)
    coef = np.array([0.0, 0.01, -0.005, 0.008, 0.0, -0.004])
    res = sp.eval_delta_flat(coef, x, y) + rng.normal(0, 0.002, x.size)
    fit = sp.fit_delta_flat(x, y, res, np.ones_like(x))
    assert np.allclose(fit[1:], coef[1:], atol=5e-4)


def test_colour_solution_inverts_the_transformation():
    cg, ci = 0.10, -0.01
    for col in (-0.5, 0.0, 0.8, 1.5):
        g_nat = 12.0 + cg * (col - sp.COLOUR_PIVOT)
        i_nat = 12.0 - col + ci * (col - sp.COLOUR_PIVOT)
        assert sp.colour_from_natural(g_nat, i_nat, cg, ci) == pytest.approx(col)
        assert sp.natural_to_ps1(g_nat, col, cg) == pytest.approx(12.0)


def test_weighted_mean():
    mu, err, c2, n = sp.weighted_mean([1.0, 3.0], [1.0, 1.0])
    assert (mu, n) == (2.0, 2)
    assert err == pytest.approx(1 / math.sqrt(2))
    assert c2 == pytest.approx(2.0)


def test_joint_search_finds_a_strong_bump_and_not_noise():
    rng = np.random.default_rng(7)
    t = np.sort(rng.uniform(5, 50, 32))
    s = np.full(t.size, 0.01)
    trend = 11 + 0.03 * (t - 5)
    J = sp.JointSearch(t, s, 8.0)
    y0 = trend + rng.normal(0, 0.01, t.size)
    y1 = y0 - 0.2 * np.exp(-0.5 * ((t - t[15]) / 1.0) ** 2)
    assert J.bump(y1)[0] > 3 * J.bump(y0)[0]
    a = J.bump_at(y1, t[15], 1.0) - J.bump_at(y0, t[15], 1.0)
    assert a == pytest.approx(-0.2, abs=0.03)
    # a pure trend leaves nothing: residual about the trend projector is ~0
    assert np.allclose(J.trend(trend), trend, atol=1e-6)


def test_knot_choice_respects_floor():
    t = np.linspace(5, 50, 30)
    y = 11 + 0.03 * (t - 5)
    ks, bic = sp.choose_knot_spacing(t, y, np.full(t.size, 0.01))
    assert ks >= sp.TREND_KNOT_SPACING_D
    assert all(k >= sp.TREND_KNOT_SPACING_D for k in bic)


def test_calib_star_mask_rejects_disc_knots_without_motion():
    n = 3
    cat = {"rmag": np.full(n, 14.0), "gmag": np.full(n, 14.6),
           "imag": np.full(n, 13.7), "rp1": np.full(n, 99.9),
           "dupvar": np.full(n, 2.0),
           "RA_ICRS": np.array([sp.M101_RA, sp.M101_RA, 211.6]),
           "DE_ICRS": np.array([sp.M101_DEC, sp.M101_DEC, 54.0]),
           "pmRA": np.array([0.1, 20.0, 0.1]), "e_pmRA": np.ones(n),
           "pmDE": np.array([0.1, 0.0, 0.1]), "e_pmDE": np.ones(n),
           "Plx": np.zeros(n), "e_Plx": np.ones(n)}
    m = sp.calib_star_mask(cat)
    assert list(m) == [False, True, True]
