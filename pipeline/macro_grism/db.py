"""The grism track's manifest tables (NEW tables only — S0/S0b/S1 tables
are never modified; the manifest stays append-only across stages).

* ``g_gate_calib``   — the era's camera-orientation calibration: one row
                       per solved imaging frame, plus the adopted CD.
* ``g_extractions``  — one row per (frame, background method): the gate
                       verdict, trace statistics, wavelength anchors, the
                       Halpha snippet, and the flanking-vs-dark debt.
* ``g_build_meta``   — code version, timestamps, sample definition.

Writes go through one transaction per batch; a crashed run leaves the
last complete batch, and the runner resumes by skipping rows that exist
(the (obs_rowid, method) uniqueness is the resume key).
"""

from __future__ import annotations

import sqlite3

G_SCHEMA = """
CREATE TABLE IF NOT EXISTS g_gate_calib (
    calib_id     INTEGER PRIMARY KEY,
    era_id       INTEGER,
    frame_path   TEXT,           -- the imaging frame that was solved
    night        TEXT,
    status       TEXT,           -- solved / unsolved / bad_solve
    cd1_1 REAL, cd1_2 REAL, cd2_1 REAL, cd2_2 REAL,   -- deg/px
    pixscale_arcsec REAL,
    rotation_deg REAL,
    rms_arcsec   REAL,
    n_matched    INTEGER,
    adopted      INTEGER DEFAULT 0    -- 1 = the CD the gate uses
);
CREATE TABLE IF NOT EXISTS g_extractions (
    obs_rowid    INTEGER NOT NULL,   -- frames.obs_rowid
    method       TEXT NOT NULL,      -- 'flanking' | 'masterdark'
    path         TEXT,
    target       TEXT,
    filter       TEXT,               -- hrg | lrg
    night        TEXT,
    jd           REAL,               -- header JD, UTC exposure START (S3
                                     -- owns BJD; nothing here converts)
    exptime      REAL,
    era_id       INTEGER,
    role         TEXT,               -- 'tcrb_sample' | 'gate_bad' |
                                     -- 'gate_good' | 'calibrator'
    layout       TEXT,               -- FITS packaging that was resolved
    -- identity gate ------------------------------------------------------
    gate_verdict TEXT,               -- ACCEPT | REJECT
    gate_reason  TEXT,
    pointing_offset_deg REAL,
    u_obs REAL, u_pred REAL, u_resid_px REAL,
    gate_parity  TEXT,
    n_gaia       INTEGER,
    brightest_g  REAL,
    -- trace ---------------------------------------------------------------
    trace_height REAL,
    trace_slope  REAL,
    trace_c0 REAL, trace_c1 REAL, trace_c2 REAL,   -- centers poly (deg 2)
    trace_rms_px REAL,
    trace_n_centroids INTEGER,
    -- extraction ----------------------------------------------------------
    bg_method    TEXT,               -- 'flanking' | 'masterdark+flanking'
    dark_path    TEXT,               -- master used (masterdark rows only)
    dark_exptime REAL,
    n_extracted  INTEGER,
    n_sat_cols   INTEGER,            -- columns with any saturated pixel
    peak_flux    REAL,               -- max optimal flux (ADU)
    median_flux  REAL,
    -- wavelength ----------------------------------------------------------
    anchor_status TEXT,              -- 'halpha+o2' | 'halpha_only' | 'none'
    x_halpha REAL, halpha_snr REAL, halpha_width_px INTEGER,
    x_o2b REAL, o2b_snr REAL,
    x_o2a REAL, o2a_snr REAL,
    disp_a_per_px REAL,              -- per-frame measurement (signed)
    disp_source  TEXT,
    -- products ------------------------------------------------------------
    snippet_json TEXT,               -- [x, flux] pairs around Halpha
    spectrum_fits TEXT,              -- products-relative output path
    contamination_flag TEXT,         -- 'C4_Be_Halpha' on tet CrB rows
    debt_median_rel_diff REAL,       -- flanking vs masterdark (both rows
                                     -- of a frame carry the same number)
    status       TEXT,               -- 'ok' | error text
    UNIQUE (obs_rowid, method)
);
CREATE TABLE IF NOT EXISTS g_build_meta (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""


def ensure_schema(con: sqlite3.Connection) -> None:
    """Create the g_* tables when absent (idempotent — IF NOT EXISTS)."""
    con.executescript(G_SCHEMA)


def existing_keys(con: sqlite3.Connection) -> set:
    """(obs_rowid, method) pairs already recorded — the resume key set."""
    try:
        return set(con.execute(
            "SELECT obs_rowid, method FROM g_extractions"))
    except sqlite3.OperationalError:
        return set()


def insert_extraction(con: sqlite3.Connection, row: dict) -> None:
    """Insert one g_extractions row from a plain dict (missing keys become
    NULL).  ``INSERT OR REPLACE`` on the (obs_rowid, method) key makes a
    deliberate re-run of a frame overwrite its old record instead of
    stacking duplicates."""
    cols = [r[1] for r in con.execute(
        "PRAGMA table_info(g_extractions)")]
    vals = [row.get(c) for c in cols]
    con.execute(
        f"INSERT OR REPLACE INTO g_extractions ({','.join(cols)}) "
        f"VALUES ({','.join('?' * len(cols))})", vals)


def set_meta(con: sqlite3.Connection, key: str, value: str) -> None:
    con.execute("INSERT OR REPLACE INTO g_build_meta (key, value) "
                "VALUES (?, ?)", (key, str(value)))
