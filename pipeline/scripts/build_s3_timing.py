#!/usr/bin/env python
"""Build the S3 time axis: audits, mid-exposure BJD_TDB, clock validation.

WHAT THIS SCRIPT DOES (stage S3 of the shared pipeline)
-------------------------------------------------------
Establishes — with evidence, not assumption — what instant every header
time stamp refers to, then computes mid-exposure BJD_TDB for EVERY
canonical science frame and validates the observatory clock against an
astrophysical standard.  It AUGMENTS the manifest with NEW tables only
(S0/S0b/S1/S2 tables are never touched):

* ``s3_dateobs_audit``   — full-archive JD vs DATE-OBS consistency, one
                           row per readout family + an outlier list
                           (``s3_dateobs_outliers``).
* ``s3_header_audit``    — sampled FITS headers across era families x
                           years: JD-HELIO semantics (start vs mid),
                           TELUT behavior, per-frame residuals against
                           our own heliocentric computation.
* ``s3_cadence``         — back-to-back series cadence per readout mode x
                           exposure time; the StackPro dead-time bound.
* ``frame_times``        — the product: one row per canonical science
                           frame (keyed by path) with jd_utc_mid,
                           bjd_tdb, the correction terms, the coordinates
                           used, and method identifiers.
* ``s3_clock_points`` / ``s3_clock_eclipses`` — AG LMi eclipse photometry
                           and the fitted O-C against the VSX ephemeris:
                           the observatory clock bound.
* ``s3_build_meta``      — timestamp, code version, ephemeris used, git
                           commit, the VSX ephemeris + its source.

Unless ``--skip-report`` is given (or a partial ``--stage`` ran), it then
renders ``docs/pipeline/s3_timing.html`` + figures from the database.

READ-ONLY DISCIPLINE
--------------------
The archive under DATA/ASTRO/rlmt-archive is immutable: this script opens
FITS files read-only and writes nothing outside the manifest DB and docs/.

IDEMPOTENCE / SAFETY
--------------------
Every s3_* table is built under a temporary name and swapped into place in
one transaction, so a crashed build never leaves a half-written table.
The two slow stages (``audit-headers``, ``clock``) accumulate per-frame
rows with INSERT OR REPLACE and skip already-measured paths, so re-running
after an interruption resumes instead of restarting.

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s3_timing.py            # everything
    ... build_s3_timing.py --stage frame-times         # one stage only

Stages: audit-scan, cadence, audit-headers, frame-times, clock, report.
"""

from __future__ import annotations

import argparse
import collections
import sqlite3
import subprocess
import sys
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Make the pipeline package importable no matter where the script is
# invoked from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import timing as tm                          # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: Which manifest rows count as "canonical science frames" for frame_times:
#: canonical, and not classified as a calibration kind by IMAGETYP.  The
#: NULL-imagetyp rows (3,901 — header-error files and blank-header science)
#: are INCLUDED so that every non-calibration canonical frame has a
#: frame_times row, even if only to record why no BJD exists for it.
SCIENCE_WHERE = ("is_canonical = 1 AND "
                 "(imagetyp IS NULL OR imagetyp LIKE 'Light%')")

#: Frames per (readout family, calendar year) stratum read in the header
#: audit: the longest exposure (JD-HELIO start-vs-mid discrimination needs
#: EXPTIME >> the residual noise), the shortest, and the lexicographically
#: first (an arbitrary-but-deterministic "typical" pick).
HEADER_SAMPLE_PER_STRATUM = 3

#: In-series gap ceiling for the cadence stage: gaps beyond this multiple
#: of a series' median gap are pauses (clouds, refocus), not cadence.
CADENCE_GAP_CEILING = 3.0

#: Minimum frames for a (target, night, exptime) run to count as a series.
CADENCE_MIN_RUN = 5

# ---------------------------------------------------------------------------
# The AG LMi clock standard (stage ``clock``).
# ---------------------------------------------------------------------------
#: VSX ephemeris of AG LMi, fetched live via astroquery when possible; the
#: values below were recorded from Vizier B/vsx/vsx (OID 167169) on
#: 2026-08-18 and serve as the offline fallback.  VSX minimum epochs are
#: HELIOCENTRIC Julian Dates (HJD, UTC scale) — the clock comparison is
#: therefore done in HJD_UTC, like against like.
AGLMI = {
    "name": "AG LMi", "vsx_oid": 167169,
    "ra_deg": 161.2035, "dec_deg": 33.35331,
    "epoch_hjd": 2458211.194, "period_d": 1.3590176,
    "type": "EA", "mag_range": "V 10.55 - 11.16",
    "source": "VSX via Vizier B/vsx/vsx, recorded 2026-08-18",
}

#: Ephemeris precision floor, stated with the O-C: VSX quotes the epoch to
#: 0.001 d (+-43 s quantization) and the period to 1e-7 d (+-0.5e-7 d per
#: cycle drift).  Both are computed against the actual cycle count below.
EPOCH_QUANT_S = 0.0005 * 86400.0            # +-43.2 s
PERIOD_QUANT_D = 0.5e-7                      # +- half the last quoted digit

#: Independent period for the O-C interpretation: the Gaia DR3 eclipsing-
#: binary solution for the same star (recorded from Vizier I/358/veb on
#: 2026-08-18).  P_gaia = 1/frequency = 1/0.7358287960636182 d; its
#: uncertainty e_Freq/Freq^2 = 3.0e-5 d dwarfs VSX's quoted last digit —
#: the honest scale of how well this star's period is actually known.
#: The VSX-vs-Gaia period difference alone predicts an O-C drift of
#: (P_gaia - P_vsx) * E at cycle E — the yardstick that separates
#: "stale survey ephemeris" from "observatory clock error".
GAIA_EB = {
    "source_id": 737974651032299776,
    "period_d": 1.0 / 0.7358287960636182,          # 1.3590123 d
    "period_err_d": 1.6305177e-5 / 0.7358287960636182 ** 2,   # 3.0e-5 d
    "source": ("Gaia DR3 vari_eclipsing_binary via Vizier I/358/veb, "
               "recorded 2026-08-18"),
}

#: Nights measured by the clock stage.  Eclipse nights = the four nights
#: whose sampled phase coverage reaches |phase| < 0.02 (computed from the
#: manifest against the VSX ephemeris); baseline nights supply the
#: out-of-eclipse reference level per (readout config, filter) and are
#: subsampled to spare archive IO.
CLOCK_ECLIPSE_NIGHTS = ("2023-03-18", "2024-02-22", "2024-02-29",
                        "2024-03-01")
CLOCK_BASELINE_NIGHTS = ("2023-03-09", "2024-02-23", "2024-03-09")
CLOCK_BASELINE_EVERY = 4          # keep every 4th frame of baseline nights

#: Photometry parameters for the clock stage (rough differential
#: photometry of a V~10.5 star with a 0.6 mag eclipse — aperture sizes
#: are deliberately generous; sub-percent precision is not the goal).
CLOCK_APER_PX = 15.0              # aperture radius (GSENSE bin1 ~0.54"/px)
CLOCK_ANNULUS_PX = (25.0, 35.0)   # local sky annulus radii
CLOCK_N_COMP = 10                 # ensemble size (fixed sky positions)
CLOCK_MATCH_MAX_PX = 25.0         # WCS position sanity for forced apertures
CLOCK_OOE_PHASE = 0.08            # |phase| beyond this = out of eclipse
CLOCK_FIT_PHASE = 0.12            # points inside this enter the dip fit
#: S2 detector facts: High Gain clips at 3,496 ADU (12-bit); a StackPro
#: frame is a 16x sum.  Comparison-star apertures whose peak approaches
#: the clip are excluded from the ensemble.
CLOCK_PEAK_VETO = {"High Gain": 3200.0, "High Gain StackPro": 16 * 3200.0}


# ---------------------------------------------------------------------------
# Small shared helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[S3] {msg}", flush=True)


def git_commit() -> str:
    """Short git commit of the repo, or '' when git is unavailable."""
    try:
        return subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        return ""


def swap_table(con: sqlite3.Connection, name: str, create_sql: str,
               rows: list[tuple], insert_sql: str) -> None:
    """Build table ``name`` under a temp name and swap it in atomically."""
    tmp = f"{name}__new"
    con.execute(f"DROP TABLE IF EXISTS {tmp}")
    con.execute(create_sql.format(table=tmp))
    con.executemany(insert_sql.format(table=tmp), rows)
    # Commit the temp-table load first (sqlite3 auto-opened a transaction
    # for the INSERTs), then swap in one atomic transaction: a crash can
    # leave a stale __new table behind, never a half-replaced real one.
    con.commit()
    con.execute("BEGIN")
    con.execute(f"DROP TABLE IF EXISTS {name}")
    con.execute(f"ALTER TABLE {tmp} RENAME TO {name}")
    con.commit()


def write_meta(con: sqlite3.Connection, extra: dict) -> None:
    """Merge key/value pairs into s3_build_meta (created on first use)."""
    con.execute("""CREATE TABLE IF NOT EXISTS s3_build_meta
                   (key TEXT PRIMARY KEY, value TEXT)""")
    base = {
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "code_version": tm.S3_CODE_VERSION,
        "git_commit": git_commit(),
    }
    for k, v in {**base, **extra}.items():
        con.execute("INSERT OR REPLACE INTO s3_build_meta VALUES (?, ?)",
                    (k, str(v)))
    con.commit()


def open_header(archive: Path, rel_path: str):
    """Read the science header of one archive file (read-only).

    fpack files put the pixels (and the interesting header) in extension
    1; plain files in the primary HDU.  Returns the astropy Header, or
    None when the file is unreadable.
    """
    from astropy.io import fits
    try:
        with fits.open(archive / rel_path) as hdul:
            if len(hdul) > 1 and hdul[1].data is not None:
                return hdul[1].header
            return hdul[0].header
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Stage: audit-scan — JD vs DATE-OBS over the whole archive (DB only)
# ---------------------------------------------------------------------------
def stage_audit_scan(con: sqlite3.Connection) -> None:
    """Compare header JD against DATE-OBS for every canonical frame.

    Both cards were captured into the catalog at ingest, so this scan
    costs no archive IO.  If JD is the same instant DATE-OBS names, the
    difference is bounded by the cards' write precision (~2 ms); anything
    larger is an outlier worth naming individually.
    """
    rows = con.execute("""
        SELECT path, night, readoutm, date_obs, jd FROM frames
        WHERE is_canonical = 1 AND jd IS NOT NULL
          AND date_obs IS NOT NULL AND date_obs != ''""").fetchall()
    per_family: dict[str, list[float]] = collections.defaultdict(list)
    outliers: list[tuple] = []
    n_unparseable = 0
    for path, night, readoutm, date_obs, jd in rows:
        jd_from_date = tm.parse_date_obs(date_obs)
        if jd_from_date is None:
            n_unparseable += 1
            continue
        diff_s = (jd - jd_from_date) * 86400.0
        per_family[(readoutm or "").strip() or "(blank)"].append(diff_s)
        if abs(diff_s) > 0.5:
            outliers.append((path, night, readoutm, round(diff_s, 3)))
    audit_rows = []
    for fam, diffs in sorted(per_family.items(), key=lambda kv: -len(kv[1])):
        a = np.array(diffs)
        audit_rows.append((
            fam, len(a), float(np.median(a)), float(np.percentile(a, 1)),
            float(np.percentile(a, 99)), float(np.max(np.abs(a))),
            int(np.sum(np.abs(a) > 0.1))))
    swap_table(con, "s3_dateobs_audit", """
        CREATE TABLE {table} (
            readoutm TEXT PRIMARY KEY, n_frames INTEGER, median_s REAL,
            p1_s REAL, p99_s REAL, max_abs_s REAL, n_gt_100ms INTEGER)""",
               audit_rows, "INSERT INTO {table} VALUES (?,?,?,?,?,?,?)")
    swap_table(con, "s3_dateobs_outliers", """
        CREATE TABLE {table} (
            path TEXT PRIMARY KEY, night TEXT, readoutm TEXT, diff_s REAL)""",
               outliers, "INSERT INTO {table} VALUES (?,?,?,?)")
    write_meta(con, {"dateobs_scan_frames": len(rows),
                     "dateobs_scan_unparseable": n_unparseable,
                     "dateobs_scan_outliers": len(outliers)})
    log(f"audit-scan: {len(rows):,} frames, {len(audit_rows)} families, "
        f"{len(outliers)} outliers > 0.5 s, "
        f"{n_unparseable} unparseable DATE-OBS")


# ---------------------------------------------------------------------------
# Stage: cadence — back-to-back series gaps; the StackPro dead-time bound
# ---------------------------------------------------------------------------
def stage_cadence(con: sqlite3.Connection) -> None:
    """Measure inter-frame gaps in continuous same-config series.

    For every (target, night, exptime) run of >= CADENCE_MIN_RUN light
    frames in one readout mode, the gaps between consecutive exposure
    STARTS bound the frame's true wall-clock span: span <= gap.  The
    minimum StackPro (gap - EXPTIME) across the archive is therefore an
    upper bound on the total internal dead time of a StackPro frame —
    the number that caps the mid-time policy's worst-case error.
    """
    rows = con.execute(f"""
        SELECT readoutm, canonical_target, night, exptime, jd FROM frames
        WHERE {SCIENCE_WHERE} AND jd IS NOT NULL AND exptime > 0
          AND readoutm IS NOT NULL AND readoutm != ''
        ORDER BY readoutm, canonical_target, night, jd""").fetchall()
    runs: dict[tuple, list[float]] = collections.defaultdict(list)
    for readoutm, target, night, exptime, jd in rows:
        runs[(readoutm, target, night, round(exptime, 3))].append(jd)
    per_mode_exp: dict[tuple, list[float]] = collections.defaultdict(list)
    for (readoutm, _t, _n, exptime), jds in runs.items():
        if len(jds) < CADENCE_MIN_RUN:
            continue
        gaps = np.diff(np.array(sorted(jds))) * 86400.0
        gaps = gaps[gaps > 0]
        if not len(gaps):
            continue
        med = np.median(gaps)
        per_mode_exp[(readoutm, exptime)].extend(
            gaps[gaps < CADENCE_GAP_CEILING * med])
    out = []
    for (readoutm, exptime), gaps in sorted(per_mode_exp.items()):
        a = np.array(gaps)
        if len(a) < 10:
            continue
        out.append((readoutm, exptime, len(a),
                    float(np.min(a) - exptime),
                    float(np.percentile(a, 5) - exptime),
                    float(np.median(a) - exptime)))
    swap_table(con, "s3_cadence", """
        CREATE TABLE {table} (
            readoutm TEXT, exptime_s REAL, n_gaps INTEGER,
            min_overhead_s REAL, p5_overhead_s REAL, median_overhead_s REAL,
            PRIMARY KEY (readoutm, exptime_s))""",
               out, "INSERT INTO {table} VALUES (?,?,?,?,?,?)")
    sp = [r for r in out if tm.is_stackpro(r[0])]
    bound = min((r[3] for r in sp), default=None)
    write_meta(con, {"stackpro_deadtime_bound_measured_s":
                     f"{bound:.3f}" if bound is not None else "n/a"})
    log(f"cadence: {len(out)} (mode, exptime) cells; StackPro minimum "
        f"overhead {bound:.3f} s (module constant "
        f"{tm.STACKPRO_DEADTIME_BOUND_S})" if bound is not None
        else f"cadence: {len(out)} cells; no StackPro series found")


# ---------------------------------------------------------------------------
# Stage: audit-headers — sampled FITS headers: JD-HELIO and TELUT semantics
# ---------------------------------------------------------------------------
def pick_header_sample(con: sqlite3.Connection) -> list[tuple]:
    """Deterministic stratified sample: per (readout family, year), the
    longest exposure, the shortest, and the first path — canonical frames
    with JD, DATE-OBS, coordinates, and positive EXPTIME."""
    base_where = (f"{SCIENCE_WHERE} AND jd IS NOT NULL AND exptime > 0 "
                  "AND date_obs IS NOT NULL AND date_obs != '' "
                  "AND ra_deg IS NOT NULL")
    strata = con.execute(f"""
        SELECT DISTINCT COALESCE(NULLIF(TRIM(readoutm), ''), '(blank)'),
               substr(night, 1, 4)
        FROM frames WHERE {base_where} AND night IS NOT NULL""").fetchall()
    picked: dict[str, tuple] = {}
    for family, year in strata:
        fam_sql = ("TRIM(COALESCE(readoutm, '')) = ''" if family == "(blank)"
                   else "TRIM(readoutm) = :fam")
        for order in ("exptime DESC, path", "exptime ASC, path", "path"):
            row = con.execute(f"""
                SELECT path, era_id, readoutm, night, exptime, jd, date_obs,
                       ra_deg, dec_deg
                FROM frames WHERE {base_where} AND {fam_sql}
                  AND substr(night, 1, 4) = :year
                ORDER BY {order} LIMIT 1""",
                {"fam": family, "year": year}).fetchone()
            if row and row[0] not in picked:
                picked[row[0]] = (family, year) + row
            if len([1 for v in picked.values()
                    if v[0] == family and v[1] == year]
                   ) >= HEADER_SAMPLE_PER_STRATUM:
                break
    return list(picked.values())


def stage_audit_headers(con: sqlite3.Connection, archive: Path,
                        ephemeris: str) -> None:
    """Read the sampled headers and test JD-HELIO / TELUT semantics.

    For each sampled frame we compute our own heliocentric correction and
    form two residuals: header JD-HELIO minus our HJD at exposure START,
    and minus our HJD at MID-exposure.  Whichever residual is ~0 tells us
    which instant the header's heliocentric stamp refers to — and, since
    the correction itself is tiny over an exposure, a mid-based JD-HELIO
    on a start-based JD proves the JD is the start.  TELUT (the telescope
    UTC clock at header write) minus DATE-OBS should grow ~1:1 with
    EXPTIME if DATE-OBS is the start (the header is written after
    readout) — a second, ephemeris-free semantics check.
    """
    con.execute("""CREATE TABLE IF NOT EXISTS s3_header_audit (
        path TEXT PRIMARY KEY, family TEXT, year TEXT, era_id INTEGER,
        night TEXT, exptime_s REAL, jd_header REAL, date_obs TEXT,
        jd_minus_dateobs_s REAL, jd_helio_header REAL, telut TEXT,
        telut_minus_dateobs_s REAL, ra_deg REAL, dec_deg REAL,
        helio_resid_start_s REAL, helio_resid_mid_s REAL)""")
    sample = pick_header_sample(con)
    done = {r[0] for r in con.execute(
        "SELECT path FROM s3_header_audit").fetchall()}
    todo = [s for s in sample if s[2] not in done]
    log(f"audit-headers: {len(sample)} sampled, {len(todo)} to read")
    n_read = 0
    for family, year, path, era_id, _ro, night, exptime, jd, date_obs, \
            ra_deg, dec_deg in todo:
        hdr = open_header(archive, path)
        if hdr is None:
            continue
        jd_h = hdr.get("JD")
        jd_helio = hdr.get("JD-HELIO")
        telut = str(hdr.get("TELUT") or "")
        date_obs_h = str(hdr.get("DATE-OBS") or date_obs)
        jd_from_date = tm.parse_date_obs(date_obs_h)
        jd_minus_dateobs = ((jd_h - jd_from_date) * 86400.0
                            if jd_h is not None and jd_from_date is not None
                            else None)
        telut_jd = tm.parse_date_obs(telut)
        telut_minus_dateobs = ((telut_jd - jd_from_date) * 86400.0
                               if telut_jd is not None
                               and jd_from_date is not None else None)
        # TELUT '1899-12-30...' means "telescope clock not set" — discard.
        if telut_minus_dateobs is not None and \
                abs(telut_minus_dateobs) > 7200.0:
            telut_minus_dateobs = None
        resid_start = resid_mid = None
        if jd_helio is not None and jd_h is not None:
            hjd_start, _ = tm.hjd_utc_from_utc(jd_h, ra_deg, dec_deg,
                                               ephemeris=ephemeris)
            mid, _m = tm.jd_utc_mid(jd_h, exptime, family)
            hjd_mid, _ = tm.hjd_utc_from_utc(mid, ra_deg, dec_deg,
                                             ephemeris=ephemeris)
            resid_start = (float(jd_helio) - float(hjd_start)) * 86400.0
            resid_mid = (float(jd_helio) - float(hjd_mid)) * 86400.0
        con.execute("INSERT OR REPLACE INTO s3_header_audit VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (path, family, year, era_id, night, exptime, jd_h,
                     date_obs_h, jd_minus_dateobs, jd_helio, telut,
                     telut_minus_dateobs, ra_deg, dec_deg,
                     resid_start, resid_mid))
        n_read += 1
        if n_read % 20 == 0:
            con.commit()
            log(f"  ... {n_read}/{len(todo)} headers")
    con.commit()
    n_total = con.execute(
        "SELECT count(*) FROM s3_header_audit").fetchone()[0]
    write_meta(con, {"header_audit_rows": n_total})
    log(f"audit-headers: {n_read} read this run, {n_total} total")


# ---------------------------------------------------------------------------
# Stage: frame-times — mid-exposure BJD_TDB for every canonical science frame
# ---------------------------------------------------------------------------
def stage_frame_times(con: sqlite3.Connection, ephemeris: str,
                      chunk: int = 25000) -> None:
    """Compute and write the ``frame_times`` table (atomic swap).

    One row per canonical science frame, keyed by path.  Frames without a
    JD get a row with method 'no_jd'; frames without coordinates get their
    UTC mid-time but a NULL BJD with method 'no_coords' — the row count of
    this table always equals the canonical-science row count, so nothing
    silently falls off the time axis.
    """
    rows = con.execute(f"""
        SELECT path, obs_rowid, era_id, readoutm, jd, exptime,
               ra_deg, dec_deg
        FROM frames WHERE {SCIENCE_WHERE} ORDER BY path""").fetchall()
    log(f"frame-times: {len(rows):,} canonical science frames")
    out: list[tuple] = []
    # First pass: the pure mid-time policy per frame.
    mids, ras, decs, computable = [], [], [], []
    for i, (path, obs_rowid, era_id, readoutm, jd, exptime,
            ra_deg, dec_deg) in enumerate(rows):
        mid, method = tm.jd_utc_mid(jd, exptime, readoutm)
        out.append([path, obs_rowid, era_id, jd, exptime, mid,
                    None, None, None, ra_deg, dec_deg, method,
                    None, tm.S3_CODE_VERSION])
        if mid is not None and ra_deg is not None and dec_deg is not None:
            computable.append(i)
            mids.append(mid)
            ras.append(ra_deg)
            decs.append(dec_deg)
        elif mid is not None:
            out[i][12] = "no_coords"
        else:
            out[i][12] = "no_jd"
    # Second pass: vectorized BJD_TDB in chunks (astropy handles ~20k
    # frames/second; chunking just keeps peak memory flat).
    bjd_method = f"bary_{ephemeris}_winer_framecenter"
    for start in range(0, len(computable), chunk):
        idx = computable[start:start + chunk]
        bjd, ltt_s, tdb_utc_s = tm.bjd_tdb_from_utc(
            np.array(mids[start:start + chunk]),
            np.array(ras[start:start + chunk]),
            np.array(decs[start:start + chunk]), ephemeris=ephemeris)
        bjd = np.atleast_1d(bjd)
        ltt_s = np.atleast_1d(ltt_s)
        tdb_utc_s = np.atleast_1d(tdb_utc_s)
        for j, i in enumerate(idx):
            out[i][6] = float(bjd[j])
            out[i][7] = float(ltt_s[j])
            out[i][8] = float(tdb_utc_s[j])
            out[i][12] = bjd_method
        log(f"  ... BJD for {min(start + chunk, len(computable)):,}"
            f"/{len(computable):,}")
    swap_table(con, "frame_times", """
        CREATE TABLE {table} (
            path TEXT PRIMARY KEY,
            obs_rowid INTEGER,
            era_id INTEGER,
            jd_utc_start REAL,     -- header JD (== DATE-OBS, audited)
            exptime_s REAL,
            jd_utc_mid REAL,       -- start + EXPTIME/2 (mid_method)
            bjd_tdb REAL,          -- the S3 product
            bary_ltt_s REAL,       -- light-travel term applied, seconds
            tdb_minus_utc_s REAL,  -- scale offset applied, seconds
            ra_deg REAL,           -- position the correction pointed at
            dec_deg REAL,          --   (frame center; recompute for
                                   --    off-center targets)
            mid_method TEXT,
            bjd_method TEXT,
            code_version TEXT)""",
               [tuple(r) for r in out],
               "INSERT INTO {table} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)")
    n_bjd = sum(1 for r in out if r[6] is not None)
    write_meta(con, {"frame_times_rows": len(out),
                     "frame_times_with_bjd": n_bjd,
                     "ephemeris": ephemeris})
    log(f"frame-times: wrote {len(out):,} rows ({n_bjd:,} with BJD_TDB)")


# ---------------------------------------------------------------------------
# Stage: clock — AG LMi eclipse photometry against the VSX ephemeris
# ---------------------------------------------------------------------------
def fetch_vsx_ephemeris() -> dict:
    """Try to fetch the AG LMi ephemeris live from VSX (Vizier B/vsx);
    fall back to the recorded constants.  Returns the dict actually used,
    with its provenance in ['source']."""
    try:
        from astroquery.vizier import Vizier
        res = Vizier(columns=["**"]).query_constraints(
            catalog="B/vsx/vsx", Name=AGLMI["name"])
        row = res[0][0]
        eph = dict(AGLMI)
        eph.update({
            "ra_deg": float(row["RAJ2000"]),
            "dec_deg": float(row["DEJ2000"]),
            "epoch_hjd": float(row["Epoch"]),
            "period_d": float(row["Period"]),
            "source": ("VSX via Vizier B/vsx/vsx, fetched live "
                       f"{datetime.now(timezone.utc).date().isoformat()}"),
        })
        return eph
    except Exception as e:
        eph = dict(AGLMI)
        eph["source"] += f" (live fetch failed: {type(e).__name__})"
        return eph


def clock_frame_list(con: sqlite3.Connection) -> list[tuple]:
    """The AG LMi frames the clock stage measures: all frames of the
    eclipse nights + every Nth frame of the baseline nights."""
    nights = CLOCK_ECLIPSE_NIGHTS + CLOCK_BASELINE_NIGHTS
    marks = ",".join("?" for _ in nights)
    rows = con.execute(f"""
        SELECT path, night, readoutm, filter, exptime, jd
        FROM frames
        WHERE {SCIENCE_WHERE} AND target_key = 'aglmi'
          AND night IN ({marks}) AND jd IS NOT NULL AND exptime > 0
        ORDER BY night, jd""", nights).fetchall()
    keep = []
    per_night_count: dict[str, int] = collections.defaultdict(int)
    for row in rows:
        night = row[1]
        if night in CLOCK_BASELINE_NIGHTS:
            if per_night_count[night] % CLOCK_BASELINE_EVERY:
                per_night_count[night] += 1
                continue
            per_night_count[night] += 1
        keep.append(row)
    return keep


def pick_comparison_stars(archive: Path, ref_path: str, target_radec,
                          n_comp: int) -> list[tuple[float, float]]:
    """Pick the fixed comparison-star sky positions from one reference
    frame: the brightest unsaturated detections away from the target and
    the frame edges, converted to ICRS through the frame's own WCS."""
    import sep
    import warnings
    from astropy.io import fits
    from astropy.wcs import WCS, FITSFixedWarning
    with fits.open(archive / ref_path) as hdul:
        hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None \
            else hdul[0]
        data = np.ascontiguousarray(hdu.data, dtype=np.float32)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(hdu.header)
    sep.set_extract_pixstack(1_000_000)
    bkg = sep.Background(data)
    sub = data - bkg.back()
    objs = sep.extract(sub, 10.0, err=bkg.globalrms, minarea=9)
    tx, ty = wcs.world_to_pixel_values(*target_radec)
    h, w = data.shape
    order = np.argsort(-objs["flux"])
    picks: list[tuple[float, float]] = []
    for i in order:
        x, y = float(objs["x"][i]), float(objs["y"][i])
        if min(x, y, w - x, h - y) < 80:            # stay off the edges
            continue
        if (x - tx) ** 2 + (y - ty) ** 2 < 60.0 ** 2:   # not the target
            continue
        if objs["peak"][i] + bkg.globalback > 2800.0:   # High Gain veto
            continue
        ra, dec = wcs.pixel_to_world_values(x, y)
        picks.append((float(ra), float(dec)))
        if len(picks) >= n_comp:
            break
    return picks


def measure_clock_frame(archive: Path, path: str, target_radec,
                        comps: list[tuple[float, float]],
                        readoutm: str) -> dict | None:
    """Forced aperture photometry at fixed sky positions on one frame.

    Returns target flux, summed ensemble flux (unsaturated comps only),
    and their errors — or None when the frame has no usable WCS or the
    target lands off-frame.
    """
    import sep
    import warnings
    from astropy.io import fits
    from astropy.wcs import WCS, FITSFixedWarning
    try:
        with fits.open(archive / path) as hdul:
            hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None \
                else hdul[0]
            data = np.ascontiguousarray(hdu.data, dtype=np.float32)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", FITSFixedWarning)
                wcs = WCS(hdu.header)
            if not wcs.has_celestial:
                return None
            egain = hdu.header.get("EGAIN")
    except Exception:
        return None
    h, w = data.shape
    positions = [target_radec] + list(comps)
    xy = np.array([wcs.world_to_pixel_values(ra, dec)
                   for ra, dec in positions], dtype=float)
    if not (CLOCK_MATCH_MAX_PX <= xy[0, 0] < w - CLOCK_MATCH_MAX_PX
            and CLOCK_MATCH_MAX_PX <= xy[0, 1] < h - CLOCK_MATCH_MAX_PX):
        return None                       # target off-frame / bad solve
    sep.set_extract_pixstack(1_000_000)
    bkg = sep.Background(data)
    sub = data - bkg.back()
    gain = float(egain) if egain else None
    flux, fluxerr, _flag = sep.sum_circle(
        sub, xy[:, 0], xy[:, 1], CLOCK_APER_PX, err=bkg.globalrms,
        gain=gain, bkgann=CLOCK_ANNULUS_PX)
    # Peak-based saturation veto per aperture (crude local max).
    veto = CLOCK_PEAK_VETO.get((readoutm or "").strip(), 3200.0)
    keep = []
    for k in range(1, len(positions)):
        x, y = int(round(xy[k, 0])), int(round(xy[k, 1]))
        if not (2 <= x < w - 2 and 2 <= y < h - 2):
            continue
        if float(data[y - 2:y + 3, x - 2:x + 3].max()) >= veto:
            continue
        if flux[k] <= 0:
            continue
        keep.append(k)
    if len(keep) < 3 or flux[0] <= 0:
        return None
    ens = float(np.sum(flux[keep]))
    ens_err = float(np.sqrt(np.sum(fluxerr[keep] ** 2)))
    return {"flux_t": float(flux[0]), "fluxerr_t": float(fluxerr[0]),
            "flux_ens": ens, "fluxerr_ens": ens_err, "n_comp": len(keep)}


def stage_clock(con: sqlite3.Connection, archive: Path,
                ephemeris: str) -> None:
    """Measure the AG LMi eclipse nights and fit the clock offset."""
    eph = fetch_vsx_ephemeris()
    log(f"clock: ephemeris {eph['epoch_hjd']} + E x {eph['period_d']} "
        f"({eph['source']})")
    con.execute("""CREATE TABLE IF NOT EXISTS s3_clock_points (
        path TEXT PRIMARY KEY, night TEXT, readoutm TEXT, filter TEXT,
        exptime_s REAL, jd_utc_mid REAL, hjd_utc_mid REAL, phase REAL,
        flux_t REAL, fluxerr_t REAL, flux_ens REAL, fluxerr_ens REAL,
        n_comp INTEGER, dmag REAL, dmag_err REAL)""")
    frames = clock_frame_list(con)
    done = {r[0] for r in con.execute(
        "SELECT path FROM s3_clock_points").fetchall()}
    todo = [f for f in frames if f[0] not in done]
    log(f"clock: {len(frames)} frames selected, {len(todo)} to measure")
    target_radec = (eph["ra_deg"], eph["dec_deg"])
    # Reference frame for comparison-star selection: the first frame of
    # the best-covered eclipse night (deterministic).
    comps_meta = con.execute("""SELECT value FROM s3_build_meta
        WHERE key = 'clock_comps'""").fetchone()
    if comps_meta:
        comps = [tuple(map(float, pair.split(",")))
                 for pair in comps_meta[0].split(";")]
    else:
        ref = next(f for f in frames if f[1] == "2024-02-22")
        comps = pick_comparison_stars(archive, ref[0], target_radec,
                                      CLOCK_N_COMP)
        write_meta(con, {"clock_comps":
                         ";".join(f"{ra:.6f},{dec:.6f}"
                                  for ra, dec in comps),
                         "clock_comps_ref": ref[0]})
    log(f"clock: {len(comps)} fixed comparison positions")
    n_done = 0
    for path, night, readoutm, filt, exptime, jd in todo:
        mid, _method = tm.jd_utc_mid(jd, exptime, readoutm)
        meas = measure_clock_frame(archive, path, target_radec, comps,
                                   readoutm)
        if meas is None:
            # Record the failure so a re-run does not retry forever.
            con.execute("INSERT OR REPLACE INTO s3_clock_points "
                        "(path, night, readoutm, filter, exptime_s, "
                        " jd_utc_mid) VALUES (?,?,?,?,?,?)",
                        (path, night, readoutm, filt, exptime, mid))
            continue
        hjd, _ = tm.hjd_utc_from_utc(mid, *target_radec,
                                     ephemeris=ephemeris)
        phase = float(tm.fold_phase(float(hjd), eph["epoch_hjd"],
                                    eph["period_d"]))
        dmag = -2.5 * np.log10(meas["flux_t"] / meas["flux_ens"])
        dmag_err = (2.5 / np.log(10)) * np.sqrt(
            (meas["fluxerr_t"] / meas["flux_t"]) ** 2
            + (meas["fluxerr_ens"] / meas["flux_ens"]) ** 2)
        con.execute("INSERT OR REPLACE INTO s3_clock_points VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (path, night, readoutm, filt, exptime, mid, float(hjd),
                     phase, meas["flux_t"], meas["fluxerr_t"],
                     meas["flux_ens"], meas["fluxerr_ens"],
                     meas["n_comp"], float(dmag), float(dmag_err)))
        n_done += 1
        if n_done % 20 == 0:
            con.commit()
            log(f"  ... {n_done}/{len(todo)} frames measured")
    con.commit()
    fit_clock(con, eph)


def fit_clock(con: sqlite3.Connection, eph: dict) -> None:
    """Baseline-correct the photometry, fit the dip, write the verdict."""
    rows = con.execute("""
        SELECT night, readoutm, filter, hjd_utc_mid, phase, dmag, dmag_err
        FROM s3_clock_points WHERE dmag IS NOT NULL""").fetchall()
    if not rows:
        log("clock: no usable photometry; nothing to fit")
        return
    # Out-of-eclipse baseline per (readout config, filter): raw
    # uncalibrated ensembles differ between the StackPro-R 2023 config
    # and the High-Gain-gri 2024 config, so each gets its own zero.
    groups = collections.defaultdict(list)
    for night, readoutm, filt, hjd, ph, dm, de in rows:
        groups[((readoutm or "").strip(), (filt or "").strip())].append(
            (night, hjd, ph, dm, de))
    pts = []
    for (config, filt), members in groups.items():
        ooe = [dm for _n, _h, ph, dm, _e in members
               if abs(ph) > CLOCK_OOE_PHASE]
        if len(ooe) < 3:
            continue                      # no baseline -> group unusable
        base = float(np.median(ooe))
        for night, hjd, ph, dm, de in members:
            pts.append((night, config, filt, hjd, ph, dm - base, de))
    fit_pts = [p for p in pts if abs(p[4]) <= CLOCK_FIT_PHASE]
    phases = np.array([p[4] for p in fit_pts])
    dmags = np.array([p[5] for p in fit_pts])
    errs = np.array([max(p[6], 0.01) for p in fit_pts])
    fit = tm.fit_eclipse_offset(phases, dmags, errs)
    results: list[tuple] = []
    # Global fit row + per-night rows where a night alone constrains it.
    def summarize(tag: str, sel_nights) -> None:
        sel = [p for p in fit_pts if p[0] in sel_nights] \
            if sel_nights else fit_pts
        if len(sel) < 8:
            return
        ph_a = np.array([p[4] for p in sel])
        f = tm.fit_eclipse_offset(ph_a,
                                  np.array([p[5] for p in sel]),
                                  np.array([max(p[6], 0.01) for p in sel]))
        if f["ph0"] is None:
            results.append((tag, len(sel), float(ph_a.min()),
                            float(ph_a.max()), None, None, None, None,
                            None, "no_dip_found"))
            return
        oc_s = f["ph0"] * eph["period_d"] * 86400.0
        oc_err_s = f["ph0_err"] * eph["period_d"] * 86400.0
        # Cycle count at the epoch of these points, for the drift term.
        mean_hjd = float(np.mean([p[3] for p in sel]))
        cycles = abs(mean_hjd - eph["epoch_hjd"]) / eph["period_d"]
        eph_sys_s = EPOCH_QUANT_S + PERIOD_QUANT_D * cycles * 86400.0
        bound_s = abs(oc_s) + oc_err_s + eph_sys_s
        results.append((tag, len(sel), float(ph_a.min()),
                        float(ph_a.max()), f["depth"], f["width"],
                        oc_s, oc_err_s, bound_s, "ok"))
    summarize("global", None)
    for night in CLOCK_ECLIPSE_NIGHTS:
        summarize(night, {night})
    swap_table(con, "s3_clock_eclipses", """
        CREATE TABLE {table} (
            tag TEXT PRIMARY KEY, n_points INTEGER,
            phase_min REAL, phase_max REAL, depth_mag REAL, width_phase REAL,
            o_minus_c_s REAL, o_minus_c_err_s REAL, clock_bound_s REAL,
            status TEXT)""",
               results, "INSERT INTO {table} VALUES (?,?,?,?,?,?,?,?,?,?)")
    # The O-C interpretation yardstick: how much offset the Gaia-vs-VSX
    # period disagreement alone predicts at our mean observation epoch.
    mean_hjd_all = float(np.mean([p[3] for p in fit_pts])) if fit_pts \
        else eph["epoch_hjd"]
    cycles_all = abs(mean_hjd_all - eph["epoch_hjd"]) / eph["period_d"]
    gaia_drift_s = (GAIA_EB["period_d"] - eph["period_d"]) \
        * cycles_all * 86400.0
    gaia_envelope_s = GAIA_EB["period_err_d"] * cycles_all * 86400.0
    write_meta(con, {
        "vsx_name": eph["name"], "vsx_epoch_hjd": eph["epoch_hjd"],
        "vsx_period_d": eph["period_d"], "vsx_source": eph["source"],
        "vsx_type": eph["type"], "vsx_mag_range": eph["mag_range"],
        "clock_fit_points": len(fit_pts),
        "clock_mean_cycle": f"{cycles_all:.1f}",
        "gaia_period_d": f"{GAIA_EB['period_d']:.7f}",
        "gaia_period_err_d": f"{GAIA_EB['period_err_d']:.2e}",
        "gaia_source": GAIA_EB["source"],
        "gaia_predicted_oc_s": f"{gaia_drift_s:.0f}",
        "gaia_oc_envelope_s": f"{gaia_envelope_s:.0f}",
    })
    g = next((r for r in results if r[0] == "global"), None)
    if g and g[6] is not None:
        log(f"clock: global O-C = {g[6]:+.0f} +- {g[7]:.0f} s; "
            f"|clock error| bound {g[8]:.0f} s")
    else:
        log("clock: global fit did not converge — see s3_clock_eclipses")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
STAGES = ("audit-scan", "cadence", "audit-headers", "frame-times", "clock",
          "report")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=("Build the S3 time axis: header audits, mid-exposure "
                     "BJD_TDB for every canonical science frame, and the "
                     "AG LMi clock validation. Augments the manifest with "
                     "s3_* tables and frame_times; never modifies earlier "
                     "stages' tables; never writes into the archive."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="manifest database to augment")
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE,
                   help="immutable archive root (read-only)")
    p.add_argument("--stage", choices=STAGES, action="append",
                   help="run only the given stage(s); default: all, "
                        "in dependency order")
    p.add_argument("--ephemeris", default=None,
                   help="solar-system ephemeris (default: first of "
                        f"{tm.EPHEMERIS_PREFERENCE} that loads)")
    p.add_argument("--skip-report", action="store_true",
                   help="do not render the HTML report at the end")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.manifest.exists():
        print(f"ERROR: manifest not found: {args.manifest}",
              file=sys.stderr)
        return 2
    stages = args.stage or [s for s in STAGES if s != "report"]
    ephemeris = args.ephemeris or tm.resolve_ephemeris()
    log(f"ephemeris: {ephemeris}")
    with closing(sqlite3.connect(args.manifest)) as con:
        if "audit-scan" in stages:
            stage_audit_scan(con)
        if "cadence" in stages:
            stage_cadence(con)
        if "audit-headers" in stages:
            stage_audit_headers(con, args.archive, ephemeris)
        if "frame-times" in stages:
            stage_frame_times(con, ephemeris)
        if "clock" in stages:
            stage_clock(con, args.archive, ephemeris)
    if ("report" in stages) or (not args.stage and not args.skip_report):
        log("rendering evidence report ...")
        from macro_core import report_s3
        path = report_s3.render_report(args.manifest)
        log(f"report -> {path}")
    log("done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
