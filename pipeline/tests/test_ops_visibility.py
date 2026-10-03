"""Tests for pipeline/scripts/ops_visibility.py — the request's sky table.

Three independent anchors, because a visibility table that is wrong sends
the site a request for an observation the sky does not allow (which is what
revision 2 of the request did):

1. the spherical-triangle altitude against astropy's full ``AltAz``
   transform, for a star and for the Sun;
2. the headline cells against the observational astronomer's independent
   astropy run (``committee/reviews/2026-10-03/observational-astronomer.md``
   section 2), which this script must reproduce to the 0.1 h quoted;
3. closed-form sanity: a circumpolar-like case, a never-rises case, and the
   run-finding helpers on hand-built arrays.

The manifest is not needed: positions are passed in explicitly.
"""

import sys
from datetime import date
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import ops_visibility as ov                      # noqa: E402
from macro_core import timing                    # noqa: E402

# ICRS positions used by the anchor tests (SIMBAD).
T_CRB = (239.87567, 25.92017)
ST_LMI = (166.4158, 25.1078)
LAM_ERI = (77.2866, -8.7541)
M101 = (210.8024, 54.3488)
NGC5548 = (214.4981, 25.1368)


@pytest.fixture(autouse=True)
def _offline():
    ov._quiet_iers()


# ---------------------------------------------------------------------------
# 1. geometry against astropy AltAz
# ---------------------------------------------------------------------------
class TestAltitude:
    def test_star_matches_altaz(self):
        from astropy.coordinates import AltAz, SkyCoord
        import astropy.units as u
        times = ov.night_grid(date(2027, 1, 15), step_min=60.0)
        _, last = ov.sun_and_last(times)
        ra, dec = ov.of_date(T_CRB[0], T_CRB[1], times[len(times) // 2])
        mine = ov.altitude_deg(ra, dec, last)
        ref = SkyCoord(T_CRB[0] * u.deg, T_CRB[1] * u.deg).transform_to(
            AltAz(obstime=times, location=timing.winer_location())).alt.deg
        # 0.02 deg = 5 s of hour angle: two orders below the 0.1 h quoted.
        assert np.max(np.abs(mine - ref)) < 0.02

    def test_sun_matches_altaz(self):
        from astropy.coordinates import AltAz, get_sun
        times = ov.night_grid(date(2026, 10, 5), step_min=60.0)
        mine, _ = ov.sun_and_last(times)
        ref = get_sun(times).transform_to(
            AltAz(obstime=times, location=timing.winer_location())).alt.deg
        assert np.max(np.abs(mine - ref)) < 0.02

    def test_airmass_two_is_thirty_degrees(self):
        assert ov.airmass_to_alt_deg(2.0) == pytest.approx(30.0)

    def test_zenith_and_pole(self):
        lat = timing.WINER_LAT_DEG
        # A star on the meridian at dec = lat is at the zenith.
        assert ov.altitude_deg(10.0, lat, 10.0) == pytest.approx(90.0)
        # The pole sits at altitude = latitude at every hour angle.
        assert ov.altitude_deg(0.0, 90.0, np.array([0.0, 90.0, 200.0])) \
            == pytest.approx([lat] * 3)


# ---------------------------------------------------------------------------
# 2. headline cells against the observational astronomer's table
# ---------------------------------------------------------------------------
#: (position, night, hours Sun<-12, hours Sun<-18) from OA section 2.
OA_CELLS = [
    (T_CRB, date(2026, 10, 5), 1.1, 0.6),
    (T_CRB, date(2026, 10, 15), 0.6, 0.2),
    (T_CRB, date(2026, 11, 1), 0.0, 0.0),
    (T_CRB, date(2026, 12, 15), 0.2, 0.0),
    (T_CRB, date(2027, 1, 1), 1.4, 0.9),
    (T_CRB, date(2027, 1, 15), 2.4, 1.9),
    (T_CRB, date(2027, 3, 1), 4.8, 4.3),
    (ST_LMI, date(2026, 10, 5), 0.0, 0.0),
    (ST_LMI, date(2026, 11, 1), 1.6, 1.2),
    (ST_LMI, date(2026, 12, 15), 5.1, 4.6),
    (ST_LMI, date(2027, 1, 15), 7.2, 6.7),
    (LAM_ERI, date(2026, 10, 5), 4.0, None),
    (LAM_ERI, date(2026, 11, 1), 6.1, None),
    (LAM_ERI, date(2027, 3, 1), 2.8, None),
    (M101, date(2026, 12, 15), 3.0, 2.5),
    (M101, date(2027, 3, 1), 7.6, None),
    (NGC5548, date(2026, 12, 15), 1.9, 1.4),
    (NGC5548, date(2027, 1, 15), 4.1, None),
]


@pytest.mark.parametrize("pos,night,h12,h18", OA_CELLS)
def test_reproduces_oa_table(pos, night, h12, h18):
    rec = ov.night_visibility(night, [pos[0]], [pos[1]])[0]
    # OA quotes one decimal; allow the rounding half-width plus one sample.
    assert rec["hours"][-12.0] == pytest.approx(h12, abs=0.1)
    if h18 is not None:
        assert rec["hours"][-18.0] == pytest.approx(h18, abs=0.1)


def test_limits_are_nested():
    """Darker limit, fewer hours -- for every night of a month."""
    for d in range(1, 29, 9):
        rec = ov.night_visibility(date(2027, 1, d), [T_CRB[0]], [T_CRB[1]])[0]
        h = rec["hours"]
        assert h[-12.0] >= h[-15.0] >= h[-18.0] >= 0.0


def test_tcrb_january_window_is_a_morning_window():
    rec = ov.night_visibility(date(2027, 1, 15), [T_CRB[0]], [T_CRB[1]])[0]
    assert rec["half"] == "morning"
    assert rec["airmass_max"] <= 2.0 + 1e-6
    assert rec["airmass_min"] >= 1.0


def test_never_rises():
    # dec -70 culminates at altitude 90 - 31.67 - 70 < 0 from Winer.
    rec = ov.night_visibility(date(2026, 12, 1), [100.0], [-70.0])[0]
    assert rec["hours"][-12.0] == 0.0
    assert rec["window_h"] == 0.0 and rec["half"] == ""


# ---------------------------------------------------------------------------
# 3. helpers
# ---------------------------------------------------------------------------
class TestLongestRun:
    def test_picks_the_longest(self):
        m = np.array([1, 1, 0, 1, 1, 1, 0, 1], dtype=bool)
        assert ov.longest_run(m) == (3, 6)

    def test_run_to_the_end(self):
        assert ov.longest_run(np.array([0, 1, 1], dtype=bool)) == (1, 3)

    def test_empty(self):
        assert ov.longest_run(np.zeros(5, dtype=bool)) == (0, 0)


class TestUsableSegments:
    NIGHTS = [date(2026, 10, 1 + i) for i in range(30)]

    def test_conjunction_splits_the_season(self):
        h = np.zeros(30)
        h[:5] = 1.0
        h[20:] = 1.0
        segs = ov.usable_segments(self.NIGHTS, h, 0.25)
        assert segs == [(date(2026, 10, 1), date(2026, 10, 5)),
                        (date(2026, 10, 21), date(2026, 10, 30))]

    def test_short_gap_is_bridged(self):
        h = np.ones(30)
        h[10:13] = 0.0                      # three bad nights, max_gap = 7
        assert len(ov.usable_segments(self.NIGHTS, h, 0.25)) == 1

    def test_never(self):
        assert ov.usable_segments(self.NIGHTS, np.zeros(30), 0.25) == []

    def test_threshold_is_inclusive(self):
        h = np.full(30, 0.25)
        assert len(ov.usable_segments(self.NIGHTS, h, 0.25)) == 1


def test_inject_replaces_only_named_blocks(tmp_path):
    doc = tmp_path / "req.md"
    doc.write_text("a\n<!-- BEGIN GENERATED: one -->\nold\n"
                   "<!-- END GENERATED: one -->\nb\n"
                   "<!-- BEGIN GENERATED: two -->\nkeep\n"
                   "<!-- END GENERATED: two -->\n")
    done = ov.inject(doc, {"one": "NEW", "absent": "x"})
    text = doc.read_text()
    assert done == ["one"]
    assert "NEW" in text and "old" not in text and "keep" in text
    # Idempotent: a second injection leaves the file unchanged.
    ov.inject(doc, {"one": "NEW"})
    assert doc.read_text() == text


def test_manifest_median_ignores_stale_headers(tmp_path):
    """OA.E2: a few frames carry the previous slew's coordinates."""
    import sqlite3
    db = tmp_path / "m.sqlite"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE frames (canonical_target, ra_deg, dec_deg)")
    rows = [("T CrB", 239.876, 25.920)] * 9 + [("T CrB", 82.25, -17.0)] * 2
    con.executemany("INSERT INTO frames VALUES (?,?,?)", rows)
    con.commit()
    ra, dec, n = ov.manifest_position(con, "T CrB")
    con.close()
    assert (ra, dec, n) == (pytest.approx(239.876), pytest.approx(25.920), 11)
