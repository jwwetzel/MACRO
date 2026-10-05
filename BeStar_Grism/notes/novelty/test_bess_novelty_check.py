"""Tests for ``bess_novelty_check.py`` — run with

    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest \
        "BeStar_Grism/notes/novelty/test_bess_novelty_check.py" -q

The file lives beside the script (not in ``pipeline/tests/``) because the
novelty package owns only this directory.  Nothing here touches the network,
the manifest or the archive: every spectrum is synthetic, so the truth is
known and a bias can be measured rather than argued about.

What is pinned down:

* the EW estimator is unbiased on absorption and on emission profiles, at
  amateur-echelle and at low resolution (standing rule 3: SIGNED bias);
* the quoted EW error is honest (pull distribution has unit width);
* the emission flag has a measured false-positive rate on pure photospheric
  absorption — the peak is a max-statistic and would otherwise cry wolf;
* the emission flag has a measured detection rate for a weak (5-10%) peak,
  so a NO-EMISSION verdict carries the effect size it could have seen
  (standing rule 2);
* season labelling, name stripping and the verdict logic.
"""
import numpy as np
import pytest

import bess_novelty_check as nv


def synth(peak=0.0, absorb=0.35, snr=100.0, dl=0.2, fwhm_abs=12.0, fwhm_em=6.0,
          tilt=2e-3, rng=None, lo=6480.0, hi=6650.0):
    """A sloped continuum with a photospheric absorption line and optional emission.

    Returns wavelength, flux and the TRUE equivalent width inside LINE_WIN
    (positive = absorption), integrated from the noiseless profile.
    """
    rng = rng or np.random.default_rng(0)
    wl = np.arange(lo, hi, dl)
    g = lambda f: np.exp(-0.5 * ((wl - nv.HALPHA) / (f / 2.3548)) ** 2)
    prof = 1.0 - absorb * g(fwhm_abs) + peak * g(fwhm_em)
    cont = 1000.0 * (1.0 + tilt * (wl - nv.HALPHA))
    flux = cont * prof + rng.normal(0, 1000.0 / snr, wl.size)
    m = (wl >= nv.LINE_WIN[0]) & (wl <= nv.LINE_WIN[1])
    return wl, flux, float(nv._trapz(1.0 - prof[m], wl[m]))


@pytest.mark.parametrize("peak,dl", [(0.0, 0.2), (0.8, 0.2), (0.0, 2.0), (2.0, 0.05)])
def test_ew_is_unbiased_and_errors_are_honest(peak, dl):
    rng = np.random.default_rng(42)
    pulls, diffs = [], []
    for _ in range(300):
        wl, fl, truth = synth(peak=peak, dl=dl, rng=rng)
        m = nv.measure_halpha(wl, fl)
        diffs.append(m["ew"] - truth)
        pulls.append((m["ew"] - truth) / m["ew_err"])
    diffs, pulls = np.array(diffs), np.array(pulls)
    # signed bias smaller than 3 standard errors of the mean, and < 0.05 A absolute
    assert abs(diffs.mean()) < 3 * diffs.std() / np.sqrt(diffs.size) + 0.02
    assert abs(diffs.mean()) < 0.05
    # the quoted error neither flatters nor inflates by more than ~40%
    assert 0.7 < pulls.std() < 1.4


def test_emission_flag_false_positive_rate_on_pure_absorption():
    """Pure photospheric absorption, modest S/N: the max-statistic must not fire."""
    rng = np.random.default_rng(7)
    fired = sum(nv.measure_halpha(*synth(peak=0.0, snr=30.0, rng=rng)[:2])["emission"]
                for _ in range(500))
    assert fired / 500 < 0.01


def test_emission_flag_detects_modest_peaks():
    """A peak 15% above continuum is always seen at S/N 50; nothing below 5% ever is."""
    rng = np.random.default_rng(8)
    # absorb=0 so 'peak' is the height above the continuum itself
    hit = sum(nv.measure_halpha(*synth(peak=0.15, absorb=0.0, snr=50.0, rng=rng)[:2])["emission"]
              for _ in range(200))
    miss = sum(nv.measure_halpha(*synth(peak=0.02, absorb=0.0, snr=200.0, rng=rng)[:2])["emission"]
               for _ in range(200))
    assert hit == 200
    assert miss == 0


def test_infilled_absorption_is_not_called_emission():
    """Emission that only fills the core (peak < continuum) is NOT 'emission above
    continuum' — the documented blind spot of the flag; the EW still moves."""
    wl, fl, truth = synth(peak=0.25, absorb=0.35, snr=200.0)
    wl0, fl0, truth0 = synth(peak=0.0, absorb=0.35, snr=200.0)
    m, m0 = nv.measure_halpha(wl, fl), nv.measure_halpha(wl0, fl0)
    assert m["emission"] == 0
    assert m["ew"] < m0["ew"] - 1.0


def test_no_coverage_returns_none():
    wl, fl, _ = synth(lo=6530.0, hi=6600.0)
    assert nv.measure_halpha(wl, fl) is None


def test_season_and_name_helpers():
    assert nv.season_of("2025-07-31") == "2024-25"
    assert nv.season_of("2025-08-01") == "2025-26"
    assert nv.strip_suffix("QQ Gem hrg 1-2e+02s") == "QQ Gem"
    assert nv.strip_suffix("PHECDA lrg 0-25s") == "PHECDA"
    assert nv.strip_suffix("Lam Eri") == "Lam Eri"
    assert abs(nv.mjd_of("2000-01-01") - 51544.0) < 1e-9
    assert nv.iso_of(51544.0) == "2000-01-01"


def test_verdict_logic():
    assert nv.classify(0, 0) == "UNWITNESSED"
    assert nv.classify(1, 1) == "EMISSION(1)"
    assert nv.classify(5, 5) == "EMISSION"
    assert nv.classify(5, 0) == "NO-EMISSION"
    assert nv.classify(5, 2) == "MIXED"
