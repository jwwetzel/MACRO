"""Tests for build_dwarf_novelty.py — run from the repo root with

    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest DwarfGalaxy_AGN_Survey/notes/novelty/test_dwarf_novelty.py -q

They touch neither the archive nor the manifest: the parser is tested on the
stored LVGDB snapshots, the estimators on synthetic data, and the relations
against numbers the LVGDB pages themselves carry.
"""
import csv
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("bdn", HERE / "build_dwarf_novelty.py")
bdn = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bdn)


def lv(i):
    return bdn.parse_lvgdb((HERE / "sources" / "lvgdb" / f"lvgdb_{i}.html").read_text())


def test_lvgdb_page_with_halpha_flux():
    """NGC 5238 carries a published H-alpha flux; the parser must find it."""
    p = lv(370)
    assert p["F_Ha"] == pytest.approx(5.5e-13) and p["e_F_Ha"] == pytest.approx(8.1e-14)
    assert p["m_FUV"] == "15.30" and p["D_method"] == "TRGB"


@pytest.mark.parametrize("i", [1160, 1532, 1585, 1587, 1588, 1589, 1590, 1591, 1671])
def test_lvgdb_candidates_have_no_halpha_flux(i):
    """The fact DW-N1 rests on: no candidate page holds an H-alpha flux.
    If LVGDB adds one and ``fetch`` is re-run, this test fails — on purpose."""
    p = lv(i)
    assert p["F_Ha"] is None and p["m_Ha"] is None
    assert "V_h" in p and "A_B" in p


def test_lvgdb_blank_uv_row_is_not_misread():
    """Dw1559+46 has no FUV SFR; the H-alpha/HI cells must not slide into it."""
    assert lv(1588)["logSFR_FUV"] is None and lv(1591)["logSFR_FUV"] == "-3.06"


@pytest.mark.parametrize("i", [370, 837, 838, 1585, 1587, 1591, 1671])
def test_fuv_relation_reproduces_lvgdb_sfr(i):
    """UNGC eq. 17-18 as coded must give LVGDB's own log SFR(FUV) to 0.1 dex
    (LVGDB rounds m_FUV to 0.1 mag, hence the tolerance)."""
    p = lv(i)
    m_c = float(p["m_FUV"]) - bdn.UNGC_FUV_AB * float(p["A_B"])
    log_sfr = bdn.UNGC_FUV_CONST - 0.4 * m_c + 2 * math.log10(float(p["D_Mpc"]))
    assert log_sfr == pytest.approx(float(p["logSFR_FUV"]), abs=0.1)


def test_predicted_flux_is_sfr_equality():
    """F_pred is defined by SFR(Ha) = SFR(FUV); check the algebra numerically."""
    m, ab, d = 18.4, 0.22, 4.6
    f = bdn.predicted_flux_fuv(m, ab)
    sfr_ha = math.log10(f) + 2 * math.log10(d) + bdn.KK_HA_CONST
    sfr_fuv = bdn.UNGC_FUV_CONST - 0.4 * (m - bdn.UNGC_FUV_AB * ab) + 2 * math.log10(d)
    assert sfr_ha == pytest.approx(sfr_fuv, abs=1e-9)


def test_robust_sigma_on_integer_data_is_not_quantised():
    """Regression for the 7.413/8.896 defect: integer ADU, sigma of a few ADU,
    plus 1 % contaminating bright pixels.  Must recover sigma to 2 %."""
    rng = np.random.default_rng(42)
    for sg in (3.0, 5.3, 7.7):
        x = np.rint(rng.normal(130, sg, (512, 512)))
        x.ravel()[rng.integers(0, x.size, 2600)] += rng.exponential(300, 2600)
        _, s = bdn.robust_sigma(x)
        assert s == pytest.approx(math.sqrt(sg ** 2 + 1 / 12), rel=0.02)


def test_flux_for_one_adu_per_second_matches_memo():
    """Physicist: ZP 18.3, W 65 A -> ~8e-15 erg/s/cm2 per ADU/s."""
    assert bdn.flam_ab(18.3) * 65.0 == pytest.approx(8e-15, rel=0.05)


def test_limit_scales_as_expected():
    one = bdn.limit_3sigma([5.0], 512.0, 18.3, 10.0)
    assert bdn.limit_3sigma([5.0] * 4, 512.0, 18.3, 10.0) == pytest.approx(one / 2)
    assert bdn.limit_3sigma([5.0], 512.0, 18.3, 20.0) == pytest.approx(one * 2)
    assert bdn.limit_3sigma([5.0], 512.0, 18.3, 10.0, width_a=130.0) == pytest.approx(one * 2)


def test_band_status():
    assert bdn.band_status(None)[0] == "velocity unknown"
    assert bdn.band_status(77.0)[0] == "in band"
    assert bdn.band_status(1300.0)[0] == "in band (edge)"
    assert bdn.band_status(2033.0)[0] == "out of band (within one width)"
    assert bdn.band_status(3697.0)[0] == "out of band"


def test_verdict_rules():
    v = bdn.ha_verdict
    assert v(0, "in band", None, 1e-14, None) == "no_ha_data"
    assert v(5, "out of band", None, 1e-13, 1e-15) == "out_of_band"
    assert v(5, "in band", 3e-13, 1e-13, 1e-15) == "already_measured"
    assert v(5, "in band", None, None, 1e-15) == "no_prediction"
    assert v(5, "in band", None, 3.1e-14, 1e-14) == "informative"
    assert v(5, "velocity unknown", None, 2e-14, 1e-14) == "marginal"
    assert v(5, "in band", None, 9e-15, 1e-14) == "below_depth"


def test_literature_table_is_well_formed():
    """Every velocity has a reference key, and every key resolves."""
    rows = list(csv.DictReader(open(HERE / "literature_status.csv")))
    keys = {r["key"] for r in csv.DictReader(open(HERE / "references.csv"))}
    assert len(rows) == 19 and len({r["field"] for r in rows}) == 19
    for r in rows:
        assert r["first_listing"] in keys
        if r["v_hel_kms"]:
            assert r["v_ref"] in keys and r["D_Mpc"]
        if r["m_fuv_pub"]:
            assert r["m_fuv_ref"] in keys
        bdn.sexa_to_deg(r["ra_hms"], r["dec_dms"])


def test_block_sigma_detects_correlated_pixels():
    """White noise must give a correlation factor of 1; pixels smoothed 2x2
    (the kind of thing a stacking readout mode does) must give ~2, i.e. the
    white-noise limit would have been optimistic by that factor."""
    rng = np.random.default_rng(3)
    a, b = rng.normal(0, 7, (2, 512, 512))
    assert bdn.block_sigma(a - b) / (bdn.BLOCK * 7) == pytest.approx(1.0, abs=0.12)
    sm = [(x[:-1, :-1] + x[1:, :-1] + x[:-1, 1:] + x[1:, 1:]) / 4 for x in rng.normal(0, 7, (2, 513, 513))]
    pix = bdn.robust_sigma(sm[0] - sm[1])[1] / math.sqrt(2)
    assert bdn.block_sigma(sm[0] - sm[1]) / (bdn.BLOCK * pix) == pytest.approx(2.0, abs=0.3)
