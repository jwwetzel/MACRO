"""Unit tests for :mod:`macro_phot.revision_checks` (CV-R9, CV-R11) and for
the linearity-cap rule in :mod:`macro_phot.series`.

Each test builds a situation with a known answer — a known ramp, a known
step, a known clock offset, a reduction that differs only by a constant —
and demands it back.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import revision_checks as rc   # noqa: E402
from macro_phot import series as sr            # noqa: E402


def test_differential_mags_cancel_a_common_transparency_change():
    rng = np.random.default_rng(1)
    trans = 10 ** (-0.4 * rng.normal(0, 0.3, 50))   # clouds, frame by frame
    star = 1000.0 * trans
    comps = np.outer(trans, [500.0, 2000.0, 800.0])
    d = rc.differential_mags(star, comps)
    assert np.allclose(d, -2.5 * np.log10(1000.0 / 3300.0))


def test_differential_mags_refuse_nonpositive_flux():
    d = rc.differential_mags(np.array([1.0, -1.0, np.nan]),
                             np.array([[1.0], [1.0], [1.0]]))
    assert np.isfinite(d[0]) and np.isnan(d[1]) and np.isnan(d[2])


def test_compare_light_curves_ignores_a_zero_point():
    rng = np.random.default_rng(2)
    a = rng.normal(0, 0.05, 200)
    b = a + 0.123 + rng.normal(0, 0.004, 200)
    s = rc.compare_light_curves(a, b)
    assert s["offset"] == pytest.approx(-0.123, abs=0.002)
    assert s["rms"] == pytest.approx(0.004, rel=0.15)
    assert s["mad_sigma"] == pytest.approx(0.004, rel=0.2)


def test_ramp_fit_recovers_a_known_ramp_through_a_heavy_tail():
    rng = np.random.default_rng(3)
    x = rng.uniform(0, 1, 400)
    slope = 0.02                      # 2 per cent of flux edge to edge-ish
    y = slope * x + rng.normal(0, 0.02, 400)
    bad = rng.choice(400, 40, replace=False)
    y[bad] += rng.normal(0, 0.5, 40)  # blends and variables
    f = rc.ramp_fit(x, y)
    assert f["slope_mag"] == pytest.approx(slope, abs=3 * f["slope_err_mag"])
    assert f["n_clipped"] >= 30
    assert f["ramp_pct"] == pytest.approx(100 * rc.MAG_TO_FRAC * slope,
                                          abs=1.0)


def test_ramp_sign_positive_means_right_side_too_faint():
    x = np.linspace(0, 1, 50)
    f = rc.ramp_fit(x, 0.01 * x)       # residual = ens - cat grows to right
    assert f["ramp_pct"] > 0


def test_state_offset_recovers_a_step():
    rng = np.random.default_rng(4)
    a = 15.0 + rng.normal(0, 0.02, 300)
    b = 15.03 + rng.normal(0, 0.02, 100)
    s = rc.state_offset(a, b)
    assert s["step"] == pytest.approx(0.03, abs=4 * s["step_err"])


def test_interpolated_offset_finds_a_clock_step_in_the_middle_era():
    cyc = np.array([0, 10, 20, 30, 1000, 1010, 1020, 1030, 2000, 2010.0])
    oc = 0.05 * cyc                                   # a period error
    verified = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1], dtype=bool)
    oc = oc + np.where(verified, 0.0, -70.0)          # era-D clock offset
    sig = np.full(cyc.size, 10.0)
    r = rc.interpolated_offset(cyc, oc, sig, verified)
    assert r["offset_s"] == pytest.approx(-70.0, abs=1e-6)
    assert r["chi2_verified"] == pytest.approx(0.0, abs=1e-9)
    assert r["offset_err_s"] > 0


def test_nights_in_ranges_labels_inclusively():
    lab = rc.nights_in_ranges(["2025-01-01", "2025-06-23", "2025-07-01"],
                              [("A", "2024-12-16", "2025-06-23")])
    assert lab == ["A", "A", None]


def test_photometry_veto_is_the_lower_of_veto_and_cap():
    # High Gain era 7 (EGAIN 1.057): cap 2,950 < veto 3,200.
    assert sr.photometry_veto_adu("High Gain", 7) == 2950
    # Mode0: veto 60,200 < cap 62,200.
    assert sr.photometry_veto_adu("Mode0", 76) == 60200
    # iKon 1 MHz: cap 38,900 < veto 59,500.
    assert sr.photometry_veto_adu("1MHz High Sensitivity 16-bit", 72) \
        == 38900
    # An unmeasured mode stays unmeasured.
    assert sr.photometry_veto_adu("Low Gain") is None
