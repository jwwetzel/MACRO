"""Unit tests for macro_dw.dwcore — the Dwarf-Galaxy Hα paper's rules.

No archive, manifest or network: every case is built by hand.
"""
import math

import numpy as np
import pytest

from macro_dw import dwcore as core


# --- 1. frames ---------------------------------------------------------------
def test_night_label_noon_to_noon():
    # 2023-06-25 09:13 UTC (02:13 local) belongs to the night of 06-24.
    assert core.night_label(2460120.8846) == "2023-06-24"
    # 2023-06-25 20:00 UTC (13:00 local) starts the night of 06-25.
    assert core.night_label(2460121.3333) == "2023-06-25"


def test_disposition_precedence():
    assert core.disposition(0, "L", 9.0, False, ["fwhm"]) == "duplicate"
    # a mispointed SPECTRUM is reported as mispointed (2023-03-25 NGC 5548)
    assert core.disposition(1, "6", 9.0, None, []) == "mispointed"
    assert core.disposition(1, "6", 0.01, None, []) == "spectrum"
    assert core.disposition(1, "H", 0.01, False, []) == "unsolved"
    assert core.disposition(1, "H", 0.01, True, ["fwhm"]) == "qc_rejected"
    assert core.disposition(1, "H", 0.01, True, []) == "science"


def test_qc_reasons_thresholds():
    assert core.qc_reasons(5.0, 18.5, 18.6, 18.7) == []
    assert core.qc_reasons(12.5, 18.5, 18.6, 18.7) == ["fwhm"]
    assert core.qc_reasons(5.0, 18.2, 18.6, 18.7) == ["zp_below_night"]
    assert core.qc_reasons(5.0, 17.6, 17.6, 18.7) == ["transparency"]
    assert core.qc_reasons(None, None, None, None) == ["no_fwhm", "no_zp"]


def test_wcs_ok_fail_closed():
    assert core.wcs_ok(50, 0.5, 0.54)
    assert not core.wcs_ok(5, 0.5, 0.54)
    assert not core.wcs_ok(50, 2.0, 0.54)
    assert not core.wcs_ok(None, None, 0.54)


# --- 2. flats ----------------------------------------------------------------
def test_fit_sky_on_template_recovers_pedestal():
    rng = np.random.default_rng(1)
    yy, xx = np.mgrid[0:32, 0:32]
    tmpl = 1.0 - 0.5 * ((xx - 16) ** 2 + (yy - 16) ** 2) / 512.0   # vignetting
    block = 40.0 * tmpl + 12.0 + rng.normal(0, 0.05, tmpl.shape)
    block[3, 3] = 500.0                                              # an outlier
    s0, c = core.fit_sky_on_template(block, tmpl)
    assert s0 == pytest.approx(40.0, rel=0.01)
    assert c == pytest.approx(12.0, abs=0.2)


def test_residual_structure_plane_vs_bump():
    yy, xx = np.mgrid[0:64, 0:64]
    plane = 1.0 + 0.01 * (xx - 32) / 32
    st = core.residual_structure(plane)
    assert st["rms_noplane"] < 1e-6 < st["rms"]
    bump = 1.0 + 0.003 * np.exp(-((xx - 20) ** 2 + (yy - 40) ** 2) / 800)
    assert core.residual_structure(bump)["p95_p05"] > 0.002


# --- 3. depth and structure --------------------------------------------------
def test_roman_limit_scaling():
    a = core.roman_mu_limit(1.0, 25.0, 0.54)
    b = core.roman_mu_limit(0.1, 25.0, 0.54)
    assert b - a == pytest.approx(2.5)
    assert math.isnan(core.roman_mu_limit(0.0, 25.0, 0.54))


def test_sersic_n1_total_mag_matches_numerical_disc():
    zp, pix, re, mu0 = 25.0, 0.54, 10.0, 22.0
    img = core.exp_disc_image(801, 400, 400, mu0, re, zp, pix)
    m_num = zp - 2.5 * np.log10(img.sum())
    mu_e = mu0 + 2.5 * core.sersic_bn(1.0) / math.log(10)
    assert core.sersic_mu0(mu_e, 1.0) == pytest.approx(mu0)
    assert core.sersic_total_mag(mu_e, re, 1.0) == pytest.approx(m_num, abs=0.01)


def test_recovery_contour_interpolates():
    mu = np.array([22.0, 23.0, 24.0, 25.0])
    F = np.array([[1.0], [1.0], [0.6], [0.2]])
    assert core.recovery_contour(F, mu, 0.9)[0] == pytest.approx(23.25)
    assert core.recovery_contour(F, mu, 0.5)[0] == pytest.approx(24.25)


# --- 4. Hα ------------------------------------------------------------------
def test_effective_width_of_a_tophat():
    # A top-hat of width W: a line gives rate F*T, a continuum f_lam*W*T.
    zp, W = 25.0, 50.0
    flam1 = core.flam_per_rate(zp)            # f_lam giving 1 unit of rate
    scale_line = W * flam1                    # line flux giving 1 unit
    assert core.effective_width(scale_line, zp) == pytest.approx(W)


def test_band_offset_and_verdicts():
    assert core.ha_offset_A(229.0) == pytest.approx(5.01, abs=0.01)
    assert core.ha_offset_A(None) is None
    assert core.ha_verdict(5.1, "in") == "detected"
    assert core.ha_verdict(4.0, "in") == "marginal"
    assert core.ha_verdict(-1.0, "out") == "not_detected"
    assert core.upper_limit(-2.0, 1.0) == pytest.approx(3.0)
    assert core.upper_limit(1.0, 1.0) == pytest.approx(4.0)


# --- 5. time series ----------------------------------------------------------
def _series(seed=3, n=300):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(0, 60, n))
    band = rng.integers(0, 3, n)
    cov = rng.normal(0, 1, (n, 2))
    y = (np.array([0.0, 0.3, -0.2])[band] + 0.01 * cov[:, 0]
         + rng.normal(0, 0.005, n))
    return t, band, cov, y


def test_periodogram_recovers_injection_with_offsets_and_covariates():
    t, band, cov, y = _series()
    X = core.nuisance_basis(band, cov)
    f = core.freq_grid(t, 5.0)
    pg = core.Periodogram(t, np.full(len(t), 1 / 0.005 ** 2), X, f)
    fi = 1.37
    Y = y[:, None] + 0.02 * np.sin(2 * np.pi * fi * t[:, None] + 0.4)
    P = pg.power(Y)
    k = np.argmax(P[:, 0])
    assert abs(f[k] - fi) < 1 / np.ptp(t)
    assert pg.amplitude(Y, np.array([k]))[0] == pytest.approx(0.02, rel=0.15)
    # Without the signal the offsets and covariate leave nothing to find.
    assert pg.power(y[:, None]).max() < 0.1


def test_recovery_threshold_and_bias_are_signed():
    rng = np.random.default_rng(5)
    t, band, cov, y = _series(seed=5)
    pg = core.Periodogram(t, np.ones(len(t)), core.nuisance_basis(band, cov),
                          core.freq_grid(t, 5.0))
    thr = core.shuffle_threshold(pg, y - np.mean(y), rng, n=50)
    assert 0 < thr < 1
    a_inj = np.full(4, 10.0)
    a_rec = np.array([8.0, 9.0, 10.0, 50.0])
    ok = np.array([True, True, True, False])
    assert core.signed_cell_bias(a_rec, a_inj, ok) == pytest.approx(-0.1)


def test_amp90():
    amps = np.array([2.0, 5.0, 10.0, 20.0])
    assert core.amp90(amps, np.array([0.0, 0.5, 0.9, 1.0])) == pytest.approx(10.0)
    assert math.isnan(core.amp90(amps, np.array([0.0, 0.1, 0.2, 0.3])))


def test_sysrem_finds_common_trend():
    rng = np.random.default_rng(7)
    trend = np.sin(np.linspace(0, 6, 80))
    amp = rng.uniform(0.5, 2.0, 40)
    R = np.outer(amp, trend) * 0.01 + rng.normal(0, 0.001, (40, 80))
    comps = core.sysrem(R, np.full(R.shape, 0.001), n_comp=1)
    c = np.corrcoef(comps[0], trend)[0, 1]
    assert abs(c) > 0.99
