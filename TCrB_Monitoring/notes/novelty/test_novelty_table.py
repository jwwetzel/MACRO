"""Tests for ``build_novelty_table.py`` (TCRB-N1).

Two kinds of test, as elsewhere in this repository.

The first kind pins the PURE helpers the count rests on: the JD-to-UT-date
conversion that puts RLMT and ARAS on one calendar, the "covers H-alpha"
predicate, the resolution classes, the gap statistics, and the ARAS row
parser — including its refusal to silently drop a row whose layout changed.

The second kind guards the two ways this particular count could go quietly
wrong on real data: a truncated ARAS download (which parses cleanly and
undercounts — it happened during development, 933 rows of 2,214), and a
manifest query that forgets ``is_canonical`` and double-counts every RLMT
frame through its reduced twin.  Those tests skip when the cache or the
manifest is absent, so the file runs anywhere.

Run from the repository root::

    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest \
        TCrB_Monitoring/notes/novelty/test_novelty_table.py -q
"""

from __future__ import annotations

import importlib.util
import sqlite3
from datetime import date
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    "build_novelty_table", HERE / "build_novelty_table.py")
nt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nt)


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------
def test_jd_to_utc_date_matches_a_known_frame():
    # First RLMT hrg frame: DATE-OBS 2025-02-27T09:07:41 UT, local night
    # 2025-02-26.  The UT date, not the local night, is the unit.
    assert nt.jd_to_utc_date(2460733.88033646) == date(2025, 2, 27)
    # JD x.5 is 00:00 UT: the boundary belongs to the new day.
    assert nt.jd_to_utc_date(2460733.5) == date(2025, 2, 27)
    assert nt.jd_to_utc_date(2460733.4999) == date(2025, 2, 26)


def test_aras_listing_jd_convention():
    # The listing prints JD - 2400000; its own first row is the check:
    # 2026-09-25 19:35 UT <-> 61309.316.
    assert nt.jd_to_utc_date(2400000.0 + 61309.316) == date(2026, 9, 25)


def test_covers_requires_both_sides_of_the_line():
    assert nt.covers(3800, 7500, nt.HA_LOOSE)
    assert nt.covers(6500, 6700, nt.HA_LOOSE)          # narrow echelle window
    assert not nt.covers(6500, 6700, nt.HA_STRICT)     # ...fails the strict one
    assert not nt.covers(3800, 6500, nt.HA_LOOSE)      # stops short of the line
    assert not nt.covers(6600, 9000, nt.HA_LOOSE)      # starts past it


def test_r_classes_are_exhaustive_and_ordered():
    assert [nt.r_class(r) for r in (410, 1999.9, 2000, 6999, 7000, 86000)] == \
        ["low", "low", "mid", "mid", "high", "high"]
    with pytest.raises(ValueError):
        nt.r_class(-1)


def test_gap_stats_reports_the_longest_gap_not_just_the_median():
    d = [date(2025, 3, 1), date(2025, 3, 2), date(2025, 3, 3), date(2025, 4, 2)]
    st = nt.gap_stats(d + d)                           # duplicates collapse
    assert st["n_dates"] == 4
    assert st["median_gap_d"] == 1
    assert st["max_gap_d"] == 30
    assert nt.gap_stats([date(2025, 3, 1)])["n_dates"] == 1


def test_month_range_keeps_empty_months():
    assert nt.month_range(date(2024, 11, 20), date(2025, 2, 1)) == \
        ["2024-11", "2024-12", "2025-01", "2025-02"]


def test_orbital_phase_convention():
    # At T0 the phase is 0; half a period later it is 0.5; it wraps.
    assert nt.orbital_phase(2459978.37, 227.5528, 2459978.37) == 0.0
    assert nt.orbital_phase(2459978.37 + 113.7764, 227.5528, 2459978.37) == \
        pytest.approx(0.5)
    assert nt.orbital_phase(2459978.37 - 1.0, 227.5528, 2459978.37) == \
        pytest.approx(1.0 - 1.0 / 227.5528)
    # date_to_jd against a known frame: 2025-02-27 00:00 UT = JD 2460733.5.
    assert nt.date_to_jd(date(2025, 2, 27)) == 2460733.5


_ROW = ("<tr><td>2025-03-09</td><td>23:10</td><td>60744.465</td><td>XDu</td>"
        "<td>sro-fr</td><td>15337</td><td><svg/></td><td>6500</td><td>6700</td>"
        "<td><svg/></td><td><a href=\"spectra/x.fit\">f</a></td>"
        "<td><a href=\"figures/x.png\">p</a></td><td></td></tr>")


def test_parse_aras_row_and_case_normalisation():
    (row,) = nt.parse_aras("<table><tr><th>Date</th></tr>" + _ROW + "</table>")
    assert row["observer"] == "XDU" and row["site"] == "SRO-FR"
    assert row["resolving_power"] == 15337
    assert (row["lambda_min"], row["lambda_max"]) == (6500, 6700)
    assert row["jd"] == pytest.approx(2460744.465)
    assert row["file"] == "spectra/x.fit"


def test_parse_aras_refuses_a_changed_layout():
    broken = _ROW.replace("<td><svg/></td>", "", 1)    # twelve cells
    with pytest.raises(ValueError):
        nt.parse_aras(broken)


def test_verdict_does_not_fire_when_rlmt_is_denser():
    # Synthetic: RLMT observes more dates than ARAS in every editor month.
    monthly = [{"month": m, "aras_dates_ha": 5, "rlmt_dates_dispersed": 20,
                "rlmt_dates_labelled": 20} for m in nt.EDITOR_MONTHS]
    facts = {"aras_dates_ha": 25, "rlmt_dates_labelled": 100,
             "rlmt_dates_dispersed": 100,
             "aras_vs_rlmt_R": [{"branch": "b", "rlmt_R": 300.0,
                                 "aras_frac_spectra_above": 1.0,
                                 "aras_dates_above": 25,
                                 "aras_median_above": True}],
             "best_series": {"n_dates": 10}, "best_observer": {"n_dates": 12}}
    v = nt.verdict(monthly, facts)
    assert v["aras_denser"] is False
    assert v["rule_fires_on_any_branch"] is False
    assert v["rlmt_margin_over_best_single_series"] == 90


# ---------------------------------------------------------------------------
# Guards on the real inputs
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not nt.ARAS_RAW.exists(), reason="no cached ARAS page")
def test_cached_aras_page_is_complete():
    import gzip
    page = gzip.open(nt.ARAS_RAW, "rb").read().decode("utf-8", "replace")
    assert len(nt.parse_aras(page)) == nt.declared_count(page)


@pytest.mark.skipif(not nt.MANIFEST.exists(), reason="no manifest")
def test_rlmt_query_counts_each_exposure_once_and_cannot_write():
    rows = nt.load_rlmt()
    ids = [r["obs_rowid"] for r in rows]
    assert len(ids) == len(set(ids))
    # 240 s exposures cannot start within 240 s of each other in one filter:
    # a pair that does is the same exposure counted through two trees.
    for filt in nt.RLMT_GRISM_2025:
        jd = sorted(r["jd"] for r in rows if r["filter"] == filt)
        assert all((b - a) * 86400.0 > 200.0 for a, b in zip(jd, jd[1:]))
    con = sqlite3.connect(f"file:{nt.MANIFEST}?mode=ro", uri=True)
    with pytest.raises(sqlite3.OperationalError):
        con.execute("CREATE TABLE _novelty_should_fail (x)")
    con.close()
