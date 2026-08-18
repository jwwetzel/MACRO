"""Regression tests for the S3 build script's database-facing logic.

Everything here was written in response to an adversarial review of the
S3 track (2026-08-18).  Each test builds a TINY in-memory manifest that
reproduces one defect the review found in the real one, then asserts the
build no longer produces it:

* calibration frames wearing ``IMAGETYP = 'Light Frame'`` reaching the
  shared time axis (232 of them did, including 164 master flats);
* a readout family with NO start-vs-mid evidence being certified anyway
  (3,901 pyscope frames were);
* raw/reduced copies of one exposure sitting on the axis with
  contradictory BJDs and nothing marking either (46 did);
* a one-sided eclipse night publishing a converged O-C (two did), a
  configured night vanishing with no row at all (one did), and a clock
  bound quoted tighter than the ephemeris it was derived from.

The unit-testable PURE parts of the same fixes live in
``test_timing.py``; this file covers the parts that need SQL.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_core import timing as tm                  # noqa: E402
from scripts import build_s3_timing as b             # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures: the smallest manifest each test needs
# ---------------------------------------------------------------------------
def _frames_db() -> sqlite3.Connection:
    """An in-memory manifest with the two tables SCIENCE_WHERE reads."""
    con = sqlite3.connect(":memory:")
    con.execute("""CREATE TABLE frames (
        path TEXT PRIMARY KEY, obs_rowid INTEGER, era_id INTEGER,
        readoutm TEXT, jd REAL, exptime REAL, ra_deg REAL, dec_deg REAL,
        imagetyp TEXT, is_canonical INTEGER, night TEXT)""")
    con.execute("""CREATE TABLE calib_frames (
        path TEXT PRIMARY KEY, kind TEXT, is_master INTEGER)""")
    return con


def _add_frame(con, path, imagetyp="Light Frame", readoutm="Mode0",
               is_canonical=1, jd=2460000.0, exptime=60.0):
    con.execute("INSERT INTO frames VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (path, abs(hash(path)) % 10 ** 6, 1, readoutm, jd, exptime,
                 161.2, 33.35, imagetyp, is_canonical, "2024-02-22"))


# ---------------------------------------------------------------------------
# SCIENCE_WHERE: what is allowed onto the shared time axis
# ---------------------------------------------------------------------------
class TestScienceWhere:
    def test_header_mislabelled_calibrations_are_excluded(self):
        """A master flat with IMAGETYP='Light Frame' must not get a BJD.

        A master frame is a STACK: its header JD is not an exposure
        instant at all, so a barycentric time for it is meaningless and —
        worse — indistinguishable from a real science time on the axis.
        S0b's calib_frames classification is the authority; the header
        card is not.
        """
        con = _frames_db()
        _add_frame(con, "rawimage/2024-02-22/real_science.fts.fz")
        _add_frame(con, "Calibrations/2025-01/Flat/master_flat_g.fts")
        con.execute("INSERT INTO calib_frames VALUES (?,?,?)",
                    ("Calibrations/2025-01/Flat/master_flat_g.fts",
                     "flat", 1))
        got = {r[0] for r in con.execute(
            f"SELECT path FROM frames WHERE {b.SCIENCE_WHERE}")}
        assert got == {"rawimage/2024-02-22/real_science.fts.fz"}

    def test_real_science_with_null_imagetyp_still_qualifies(self):
        # The blank-header science frames must NOT be collateral damage:
        # they belong on the axis (even if only to record why they have
        # no BJD), so the calibration clause must not swallow them.
        con = _frames_db()
        _add_frame(con, "rawimage/2024-02-22/blank_header.fts.fz",
                   imagetyp=None)
        got = [r[0] for r in con.execute(
            f"SELECT path FROM frames WHERE {b.SCIENCE_WHERE}")]
        assert got == ["rawimage/2024-02-22/blank_header.fts.fz"]


# ---------------------------------------------------------------------------
# start_evidence: which families the two probes actually cover
# ---------------------------------------------------------------------------
class TestStartEvidence:
    @staticmethod
    def _audit_db(rows):
        con = sqlite3.connect(":memory:")
        con.execute("""CREATE TABLE s3_header_audit (
            path TEXT PRIMARY KEY, family TEXT, exptime_s REAL,
            jd_helio_header REAL, telut_minus_dateobs_s REAL,
            helio_resid_mid_s REAL)""")
        con.executemany("INSERT INTO s3_header_audit VALUES (?,?,?,?,?,?)",
                        rows)
        return con

    def test_pyscope_family_is_not_certified(self):
        """No JD-HELIO card and a TELUT that is a copy of DATE-OBS is not
        evidence of anything — those two rows are exactly what the 2026
        pyscope eras carry, and the old code called the semantics
        'proven in every era family' on their behalf."""
        con = self._audit_db([
            ("p1", "(blank)", 300.0, None, 0.0, None),
            ("p2", "(blank)", 0.25, None, 0.0, None),
            ("m1", "Mode0", 60.0, 2460000.5, 61.5, -0.02),
        ])
        proven = b.families_with_start_evidence(con)
        assert "Mode0" in proven
        assert "(blank)" not in proven

    def test_a_millisecond_exposure_proves_nothing(self):
        """The probes separate start from mid by EXPTIME/2, so a 1 ms
        frame cannot discriminate.  This is what keeps the HDR family —
        whose only long exposures give 836 s residuals — out of the
        proven set."""
        con = self._audit_db([
            ("h1", "HDR", 0.00099, 2460000.5, None, -0.016),
            ("h2", "HDR", 1024.0, 2460000.5, None, -836.0),
        ])
        assert b.families_with_start_evidence(con) == set()

    def test_telut_alone_is_enough_when_it_is_a_real_read(self):
        con = self._audit_db([
            ("f1", "Fast", 120.0, None, 121.4, None),
        ])
        assert b.families_with_start_evidence(con) == {"Fast"}


# ---------------------------------------------------------------------------
# sibling_jd_drift: the raw/reduced disagreement the JD audit cannot see
# ---------------------------------------------------------------------------
class TestSiblingDrift:
    @staticmethod
    def _links_db(rows):
        con = sqlite3.connect(":memory:")
        con.execute("""CREATE TABLE raw_reduced_links (
            raw_path TEXT, reduced_path TEXT, jd_drift_s REAL)""")
        con.executemany("INSERT INTO raw_reduced_links VALUES (?,?,?)", rows)
        return con

    def test_both_sides_of_an_on_axis_pair_are_flagged(self):
        con = self._links_db([("raw/a.fts", "red/a.fts", 271.17)])
        drifts, n_pairs = sibling = b.sibling_jd_drift(
            con, {"raw/a.fts", "red/a.fts"})
        assert n_pairs == 1
        assert drifts["raw/a.fts"] == pytest.approx(271.17)
        assert drifts["red/a.fts"] == pytest.approx(271.17)
        assert sibling  # tuple is truthy; keeps flake8 quiet about unused

    def test_links_whose_sibling_is_off_axis_are_not_flagged(self):
        # Most reduced copies are non-canonical duplicates that never
        # reach frame_times; flagging their raw parent would be noise,
        # since that raw row is then the sole authoritative stamp.
        con = self._links_db([("raw/b.fts", "red/b.fts", 0.4)])
        drifts, n_pairs = b.sibling_jd_drift(con, {"raw/b.fts"})
        assert drifts == {} and n_pairs == 0

    def test_worst_disagreement_wins_when_a_path_has_several_links(self):
        con = self._links_db([("raw/c.fts", "red/c1.fts", 5.0),
                              ("raw/c.fts", "red/c2.fts", -4199.6)])
        drifts, _ = b.sibling_jd_drift(
            con, {"raw/c.fts", "red/c1.fts", "red/c2.fts"})
        assert drifts["raw/c.fts"] == pytest.approx(4199.6)


# ---------------------------------------------------------------------------
# fit_clock: the coverage gate, the missing row, the ephemeris term
# ---------------------------------------------------------------------------
def _clock_db(points) -> sqlite3.Connection:
    """In-memory DB holding synthetic AG LMi photometry.

    ``points`` is a list of (night, phase, dmag) triples; the rest of the
    columns fit_clock reads are filled with one config so every point
    lands in a single baseline group.
    """
    con = sqlite3.connect(":memory:")
    con.execute("""CREATE TABLE s3_clock_points (
        path TEXT PRIMARY KEY, night TEXT, readoutm TEXT, filter TEXT,
        hjd_utc_mid REAL, phase REAL, dmag REAL, dmag_err REAL)""")
    eph = b.AGLMI
    for i, (night, phase, dmag) in enumerate(points):
        # An HJD ~1500 cycles past the epoch, consistent with the phase:
        # fit_clock uses it only for the cycle count in the error term.
        hjd = eph["epoch_hjd"] + (1500 + phase) * eph["period_d"]
        con.execute("INSERT INTO s3_clock_points VALUES (?,?,?,?,?,?,?,?)",
                    (f"p{i}", night, "High Gain", "g", hjd, phase, dmag,
                     0.01))
    return con


def _dip(phase, ph0=0.0, depth=0.6, width=0.025):
    """A symmetric eclipse of the shape fit_eclipse_offset models."""
    return depth * float(np.exp(-((phase - ph0) ** 2) / (2 * width ** 2)))


def _synthetic_nights():
    """One well-covered night, one one-sided arc, one starved night.

    Mirrors the real AG LMi sample: 2024-02-22 brackets the minimum,
    2023-03-18 sits entirely on one flank, 2024-03-01 has 5 points.
    """
    pts = []
    for ph in np.linspace(-0.11, 0.11, 40):          # two-sided night
        pts.append(("good", float(ph), _dip(ph)))
    for ph in np.linspace(0.0191, 0.0357, 19):       # one-sided arc
        pts.append(("onesided", float(ph), _dip(ph)))
    for ph in np.linspace(-0.085, -0.078, 5):        # starved night
        pts.append(("starved", float(ph), _dip(ph)))
    return pts


@pytest.fixture()
def fitted(monkeypatch):
    """Run fit_clock over the synthetic nights and return its rows."""
    monkeypatch.setattr(b, "CLOCK_ECLIPSE_NIGHTS",
                        ("good", "onesided", "starved"))
    con = _clock_db(_synthetic_nights())
    b.fit_clock(con, b.AGLMI)
    rows = {r[0]: r for r in con.execute(
        "SELECT tag, n_points, o_minus_c_s, o_minus_c_err_s, "
        "clock_bound_s, status FROM s3_clock_eclipses")}
    meta = dict(con.execute("SELECT key, value FROM s3_build_meta"))
    return rows, meta


class TestClockCoverageGate:
    def test_one_sided_night_is_not_a_measurement(self, fitted):
        """Regression: a symmetric template fitted to a 0.017-wide arc
        that never reaches phase 0 converged and published
        O-C = -2,231 +- 3,053 s as status 'ok'."""
        rows, _ = fitted
        assert rows["onesided"][5] == tm.CLOCK_STATUS_ONE_SIDED
        assert rows["onesided"][2] is None      # no O-C at all
        assert rows["onesided"][4] is None      # and no bound

    def test_starved_night_still_gets_a_row(self, fitted):
        """Regression: summarize() returned early on < 8 points, so the
        fourth configured night vanished from the table with no trace
        while the page's text still claimed four."""
        rows, _ = fitted
        assert "starved" in rows
        assert rows["starved"][5] == tm.CLOCK_STATUS_TOO_FEW
        assert rows["starved"][1] == 5

    def test_global_fit_uses_only_gated_nights(self, fitted):
        """Regression: the global fit folded all three nights together,
        pulling the answer away from the one trustworthy night."""
        rows, meta = fitted
        assert rows["global"][5] == tm.CLOCK_STATUS_OK
        assert rows["global"][1] == rows["good"][1]
        assert rows["global"][2] == pytest.approx(rows["good"][2], abs=1.0)
        assert meta["clock_nights_gated"] == "good"

    def test_well_covered_night_recovers_the_injected_centre(self, fitted):
        # The gate must not be so strict that a real fit stops working:
        # the dip was injected at phase 0, so O-C must come back ~0.
        rows, _ = fitted
        assert abs(rows["good"][2]) < 60.0


class TestClockBoundHonesty:
    def test_bound_includes_the_real_ephemeris_uncertainty(self, fitted):
        """Regression: the bound propagated only VSX's quoted last digit
        (0.5e-7 d), giving 519 s — while the same build had already
        recorded a 3.0e-5 d Gaia period uncertainty worth ~3,947 s.  The
        page asserted a 519 s ceiling and admitted a +-3,947 s envelope
        two sentences later.
        """
        rows, meta = fitted
        cycles = float(meta["clock_mean_cycle"])
        gaia_term = b.GAIA_EB["period_err_d"] * cycles * 86400.0
        vsx_term = b.PERIOD_QUANT_D * cycles * 86400.0
        # The credible term is the one that dominates ...
        assert gaia_term > 100 * vsx_term
        # ... and the published bound cannot be tighter than it.
        assert rows["global"][4] >= gaia_term
        assert float(meta["clock_eph_term_s"]) >= gaia_term

    def test_the_tight_wrong_number_is_still_recorded(self, fitted):
        # Kept so the report can show the difference rather than quietly
        # swapping one number for another.
        _, meta = fitted
        assert float(meta["clock_vsx_quant_term_s"]) < \
            float(meta["clock_eph_term_s"])
