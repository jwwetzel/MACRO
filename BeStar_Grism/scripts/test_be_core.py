"""Unit tests for the Be-star core numerics (truth known, no database).

Run:  /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest -q BeStar_Grism/scripts/test_be_core.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_tscore as T  # noqa: E402
import be_measure as M  # noqa: E402


def test_power_equals_floating_mean_gls():
    """With one offset column the joint-fit power is the floating-mean GLS
    (astropy LombScargle, standard normalisation, equal weights)."""
    from astropy.timeseries import LombScargle
    rng = np.random.default_rng(1)
    t = np.sort(rng.uniform(0, 100, 40))
    y = 0.5 * np.sin(2 * np.pi * t / 7.3) + rng.normal(0, 1, 40)
    N, _ = T.nuisance_design(np.zeros(40, int))
    f = np.linspace(0.02, 1.0, 300)
    P, _, _ = T.power(t, y, T.projector(N), f)
    ref = LombScargle(t, y, fit_mean=True, center_data=True).power(f, normalization="standard", method="cython")
    assert np.allclose(P, ref, atol=1e-8)


def test_degenerate_regressor_is_dropped():
    states = np.array([0] * 10 + [1] * 10)
    reg = {"temp": np.where(states == 0, -10.0, -20.0)}      # constant within each state
    N, names = T.nuisance_design(states, reg)
    assert names == ["offset[0]", "offset[1]"] and N.shape[1] == 2


def test_injected_sinusoid_recovered_at_its_frequency():
    rng = np.random.default_rng(2)
    t = np.sort(rng.uniform(0, 120, 25))
    N, _ = T.nuisance_design(np.zeros(25, int))
    f = T.freq_grid(t, 1.2, p_max=40)
    y = 6 * np.sin(2 * np.pi * t / 9.0) + rng.normal(0, 1, 25)
    P, a, b = T.power(t, y, T.projector(N), f)
    k = np.nanargmax(P)
    assert abs(f[k] - 1 / 9.0) / (1 / 9.0) < 0.01
    assert 4.5 < np.hypot(a[k], b[k]) < 7.5


def test_ew_of_gaussian_line_signed():
    """EW sign convention: absorption positive, emission negative; value to 1%."""
    w = np.arange(6400.0, 6700.0, 0.46)
    for amp, sign in ((-0.5, 1), (2.0, -1)):
        prof = amp * np.exp(-0.5 * ((w - 6563) / 3.0) ** 2)
        f = 1000 * (1 + prof)
        m = M.ew_measure(w, f, np.full_like(w, 1.0))
        truth = -amp * 3.0 * np.sqrt(2 * np.pi)
        assert np.sign(m["ew"]) == sign
        assert abs(m["ew"] - truth) < 0.01 * abs(truth)
