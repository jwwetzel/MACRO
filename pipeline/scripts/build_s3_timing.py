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
                           exposure time; the StackPro dead-time bound
                           (plus ``s3_cadence_outliers``, the pairs whose
                           stamps are too close together to be real).
* ``s3_clock_drift``     — TELUT (telescope clock) vs DATE-OBS
                           (acquisition clock) sampled across every era:
                           the archive's only RELATIVE clock check.
* ``frame_times``        — the product: one row per canonical science
                           frame (keyed by path) with jd_utc_mid,
                           bjd_tdb, the correction terms, the coordinates
                           used, method identifiers, whether the
                           start-of-exposure semantics are proven for
                           that era, and any raw-vs-reduced stamp
                           disagreement (plus ``s3_time_outliers``, the
                           exposures whose BJD was withdrawn).
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

Stages: audit-scan, cadence, audit-headers, drift, frame-times, clock,
report.
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
from macro_core import fitsgeom                              # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: Which manifest rows count as "canonical science frames" for frame_times:
#: canonical, not classified as a calibration kind by IMAGETYP, and not a
#: path S0b already catalogued as a calibration.  The NULL-imagetyp rows
#: (3,901 — header-error files and blank-header science) are INCLUDED so
#: that every non-calibration canonical frame has a frame_times row, even
#: if only to record why no BJD exists for it.
#:
#: The calib_frames clause is not redundant: 853 frames whose path says
#: dark/bias/flat carry IMAGETYP = 'Light Frame' in their headers, and 232
#: of those are the very paths S0b catalogued in calib_frames (164 master
#: flats, 60 raw flats, 8 master darks).  A MASTER is a stack — its header
#: JD is not an exposure instant at all — so a BJD for it would be
#: meaningless and indistinguishable from a real science time.  S0b's
#: classification is the authority; IMAGETYP alone is not.
SCIENCE_WHERE = ("is_canonical = 1 AND "
                 "(imagetyp IS NULL OR imagetyp LIKE 'Light%') AND "
                 "path NOT IN (SELECT path FROM calib_frames)")

#: Frames per (readout family, calendar year) stratum read in the header
#: audit: the longest exposure (JD-HELIO start-vs-mid discrimination needs
#: EXPTIME >> the residual noise), the shortest, and the lexicographically
#: first (an arbitrary-but-deterministic "typical" pick).
HEADER_SAMPLE_PER_STRATUM = 3

#: Frames per ERA read by the ``drift`` stage.  The relative-drift check
#: needs the SAME comparison repeated across the whole 2023-2026 baseline,
#: so it samples by era (the pinned registry) rather than by family/year.
DRIFT_SAMPLE_PER_ERA = 4

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
#: The phase cuts (out-of-eclipse baseline, fit window) and the coverage
#: gate live in macro_core.timing, NOT here: the report re-applies the
#: same baseline convention when it draws the folded light curve, and two
#: copies of one policy in two modules is exactly how a figure silently
#: stops matching the fit it illustrates.  Referenced below as
#: tm.CLOCK_OOE_PHASE / tm.CLOCK_FIT_PHASE.
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
    # IMMEDIATE, not a bare BEGIN.  A deferred transaction takes its READ
    # snapshot on the first statement and only asks for the write lock
    # later; in WAL mode SQLite refuses that upgrade instantly with
    # SQLITE_BUSY — no busy handler, so ``busy_timeout`` above buys
    # nothing — if any other connection committed in between.  With the
    # S1 plate-solve batch committing into this same database in bursts,
    # that is not a rare race: it is the normal case.  Asking for the
    # write lock up front puts the wait where the busy handler can serve
    # it.  (The bare BEGIN survived until now only because the first
    # statement below is a DROP of an EXISTING table, i.e. already a
    # write; on the very first build, when there is nothing to drop, the
    # old form would have raced.)
    con.execute("BEGIN IMMEDIATE")
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
    #: Per family, how many DATE-OBS strings carry a fractional-seconds
    #: field at all.  A family that writes none stamps WHOLE SECONDS, so
    #: its time axis has a 1 s granularity (up to 0.5 s of systematic if
    #: the writer truncates rather than rounds) — three orders of
    #: magnitude coarser than the ~10 ms the MaxIm families deliver, and
    #: invisible in the JD-vs-DATE-OBS comparison because BOTH cards are
    #: written from the same coarse value.
    per_family_frac: dict[str, int] = collections.defaultdict(int)
    outliers: list[tuple] = []
    n_unparseable = 0
    for path, night, readoutm, date_obs, jd in rows:
        jd_from_date = tm.parse_date_obs(date_obs)
        if jd_from_date is None:
            n_unparseable += 1
            continue
        fam = (readoutm or "").strip() or "(blank)"
        diff_s = (jd - jd_from_date) * 86400.0
        per_family[fam].append(diff_s)
        if "." in str(date_obs):
            per_family_frac[fam] += 1
        if abs(diff_s) > 0.5:
            outliers.append((path, night, readoutm, round(diff_s, 3)))
    audit_rows = []
    for fam, diffs in sorted(per_family.items(), key=lambda kv: -len(kv[1])):
        a = np.array(diffs)
        n_frac = per_family_frac[fam]
        # Stated stamping resolution: 1 s when NO frame in the family
        # writes a fraction, else the millisecond class the cards show.
        stamp_res_s = 1.0 if n_frac == 0 else 0.001
        audit_rows.append((
            fam, len(a), float(np.median(a)), float(np.percentile(a, 1)),
            float(np.percentile(a, 99)), float(np.max(np.abs(a))),
            int(np.sum(np.abs(a) > 0.1)), int(n_frac), float(stamp_res_s)))
    swap_table(con, "s3_dateobs_audit", """
        CREATE TABLE {table} (
            readoutm TEXT PRIMARY KEY, n_frames INTEGER, median_s REAL,
            p1_s REAL, p99_s REAL, max_abs_s REAL, n_gt_100ms INTEGER,
            n_with_fractional_s INTEGER, stamp_resolution_s REAL)""",
               audit_rows, "INSERT INTO {table} VALUES (?,?,?,?,?,?,?,?,?)")
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

    For every (target, night, exptime) run of >= tm.CADENCE_MIN_RUN light
    frames in one readout mode, the gaps between consecutive exposure
    STARTS bound the frame's true wall-clock span: span <= gap.

    THE ESTIMATOR MATTERS MORE THAN THE DATA HERE.  A raw minimum over
    ~15,000 StackPro gaps is an extreme order statistic that one bad time
    stamp owns outright — and does: the three smallest overheads in the
    archive are single frames landing < 1.2 s from a neighbour inside
    11-13 s cadences, and the same estimator returns NEGATIVE overheads
    (gap < EXPTIME) in several plain High Gain cells.  So this stage
    computes, per cell, BOTH statistics:

    * ``min_overhead_s`` — minimum over gaps that survive the
      impossibly-short cut (:data:`tm.CADENCE_MIN_GAP_FRACTION` of the
      series' own median), with ``n_short_discarded`` counting what the
      cut removed and ``s3_cadence_outliers`` naming every discarded
      pair, so the exclusion is auditable rather than invisible;
    * ``raw_min_overhead_s`` — the unfiltered minimum, kept precisely so
      the report can show what the old number was and why it is wrong;
    * ``regular_overhead_s`` — the SMALLEST (median gap - EXPTIME) over
      the cell's REGULAR series (see :func:`tm.series_cadence`).  This is
      the statistic the dead-time bound is taken from: a median cannot be
      moved by one bad stamp, and "the camera repeatedly delivers a frame
      every median-gap seconds" is an exact physical constraint on
      everything it does per frame, sub-read dead time included.
    """
    rows = con.execute(f"""
        SELECT readoutm, canonical_target, night, exptime, jd, path FROM frames
        WHERE {SCIENCE_WHERE} AND jd IS NOT NULL AND exptime > 0
          AND readoutm IS NOT NULL AND readoutm != ''
        ORDER BY readoutm, canonical_target, night, jd""").fetchall()
    runs: dict[tuple, list[tuple]] = collections.defaultdict(list)
    for readoutm, target, night, exptime, jd, path in rows:
        runs[(readoutm, target, night, round(exptime, 3))].append((jd, path))
    # Per (mode, exptime) cell: kept gaps, discard count, and the regular
    # series' overheads.  Outliers are collected globally.
    kept_by_cell: dict[tuple, list[float]] = collections.defaultdict(list)
    n_short_by_cell: dict[tuple, int] = collections.defaultdict(int)
    raw_by_cell: dict[tuple, list[float]] = collections.defaultdict(list)
    regular_by_cell: dict[tuple, list[float]] = collections.defaultdict(list)
    outliers: list[tuple] = []
    for (readoutm, target, night, exptime), members in runs.items():
        if len(members) < tm.CADENCE_MIN_RUN:
            continue
        members = sorted(members)                 # ascending JD, path rides
        jds = [m[0] for m in members]
        stats = tm.series_cadence(jds, exptime)
        if stats["median_gap_s"] is None:
            continue
        cell = (readoutm, exptime)
        kept_by_cell[cell].extend(float(g) for g in stats["kept_s"])
        n_short_by_cell[cell] += len(stats["short_idx"])
        # The raw (unfiltered) in-band gaps, for the "what the naive
        # estimator would have said" column.
        raw = np.diff(np.array(jds)) * 86400.0
        med = stats["median_gap_s"]
        raw_by_cell[cell].extend(
            float(g) for g in raw
            if 0 < g < tm.CADENCE_GAP_CEILING * med)
        if stats["overhead_s"] is not None:
            regular_by_cell[cell].append(float(stats["overhead_s"]))
        for i in stats["short_idx"]:
            # Name both frames of every impossibly-short pair: an
            # out-of-sequence frame arriving 0.74 s after its neighbour
            # is itself a time-stamp defect worth flagging on the axis.
            # ``series_regular`` says how strong the evidence is: in a
            # REGULAR series the median IS the machine cycle time, so a
            # gap at half of it cannot be real; in a bursty series (a few
            # snapshots an hour apart) the median is not a cadence at all
            # and a short gap may simply be a burst.
            gap_s = (jds[i + 1] - jds[i]) * 86400.0
            outliers.append((members[i][1], members[i + 1][1], readoutm,
                             exptime, night, target, float(gap_s),
                             float(med), float(gap_s / med),
                             int(stats["regular"])))
    out = []
    for cell in sorted(set(kept_by_cell) | set(regular_by_cell)):
        readoutm, exptime = cell
        a = np.array(kept_by_cell[cell])
        if len(a) < 10:
            continue
        raw_a = np.array(raw_by_cell[cell])
        reg = regular_by_cell[cell]
        out.append((readoutm, exptime, len(a), int(n_short_by_cell[cell]),
                    float(np.min(a) - exptime),
                    float(np.min(raw_a) - exptime) if len(raw_a) else None,
                    float(np.percentile(a, 5) - exptime),
                    float(np.median(a) - exptime),
                    len(reg),
                    float(min(reg)) if reg else None))
    swap_table(con, "s3_cadence", """
        CREATE TABLE {table} (
            readoutm TEXT, exptime_s REAL,
            n_gaps INTEGER,              -- gaps kept as genuine cadence
            n_short_discarded INTEGER,   -- impossibly-short, see outliers
            min_overhead_s REAL,         -- min over KEPT gaps
            raw_min_overhead_s REAL,     -- min with no short cut (naive)
            p5_overhead_s REAL, median_overhead_s REAL,
            n_regular_series INTEGER,    -- series with a machine cadence
            regular_overhead_s REAL,     -- min (median gap - EXPTIME)
            PRIMARY KEY (readoutm, exptime_s))""",
               out, "INSERT INTO {table} VALUES (?,?,?,?,?,?,?,?,?,?)")
    swap_table(con, "s3_cadence_outliers", """
        CREATE TABLE {table} (
            path_a TEXT, path_b TEXT, readoutm TEXT, exptime_s REAL,
            night TEXT, canonical_target TEXT, gap_s REAL,
            series_median_gap_s REAL, gap_over_median REAL,
            series_regular INTEGER,
            PRIMARY KEY (path_a, path_b))""",
               outliers, "INSERT OR REPLACE INTO {table} "
                         "VALUES (?,?,?,?,?,?,?,?,?,?)")
    # SENSITIVITY, measured rather than asserted.  A bound whose value
    # depends on the analyst's choice of cut is not a bound; the page
    # says so about the filtered-minimum estimator, so it owes the same
    # test of its own.  Re-derive the StackPro ceiling with each cut
    # moved well away from its default and record the spread.
    sweep: list[float] = []
    for kwargs in ({"regular_spread": 0.05}, {"regular_spread": 0.30},
                   {"regular_min_gaps": 3}, {"regular_min_gaps": 20},
                   {"min_gap_fraction": 0.3}, {"min_gap_fraction": 0.7},
                   {"gap_ceiling": 2.0}, {"gap_ceiling": 5.0}):
        alt = [s["overhead_s"]
               for (ro, _t, _n, exp), members in runs.items()
               if tm.is_stackpro(ro) and len(members) >= tm.CADENCE_MIN_RUN
               for s in [tm.series_cadence(
                   [m[0] for m in sorted(members)], exp, **kwargs)]
               if s["overhead_s"] is not None]
        if alt:
            sweep.append(min(alt))
    sp = [r for r in out if tm.is_stackpro(r[0])]
    bound = min((r[9] for r in sp if r[9] is not None), default=None)
    naive = min((r[5] for r in sp if r[5] is not None), default=None)
    filt = min((r[4] for r in sp), default=None)
    write_meta(con, {
        "stackpro_deadtime_bound_measured_s":
            f"{bound:.3f}" if bound is not None else "n/a",
        "stackpro_naive_min_overhead_s":
            f"{naive:.3f}" if naive is not None else "n/a",
        "stackpro_filtered_min_overhead_s":
            f"{filt:.3f}" if filt is not None else "n/a",
        "cadence_short_gaps_discarded": len(outliers),
        "cadence_short_gaps_in_regular_series":
            sum(1 for o in outliers if o[9]),
        "stackpro_bound_sweep_n": len(sweep),
        "stackpro_bound_sweep_min":
            f"{min(sweep):.3f}" if sweep else "n/a",
        "stackpro_bound_sweep_max":
            f"{max(sweep):.3f}" if sweep else "n/a",
    })
    if bound is None:
        log(f"cadence: {len(out)} cells; no regular StackPro series found")
        return
    log(f"cadence: {len(out)} (mode, exptime) cells; "
        f"{len(outliers)} impossibly-short gaps discarded; StackPro "
        f"dead-time bound {bound:.3f} s (naive min would have said "
        f"{naive:.3f} s)")
    # A silent drift between the measured bound and the constant the
    # policy quotes is exactly the defect this stage exists to prevent.
    if abs(bound - tm.STACKPRO_DEADTIME_BOUND_S) > 0.01:
        log(f"  WARNING: macro_core.timing.STACKPRO_DEADTIME_BOUND_S = "
            f"{tm.STACKPRO_DEADTIME_BOUND_S} does not match the measured "
            f"{bound:.3f} s — update the constant and re-run.")


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
    # Columns added after the first build; ALTER is the cheap, safe way to
    # extend an accumulating (non-swapped) table.  Rows written before the
    # extension carry NULL and are re-read below.
    have = {r[1] for r in con.execute("PRAGMA table_info(s3_header_audit)")}
    for col, decl in (("swcreate", "TEXT"), ("naxis1", "INTEGER"),
                      ("naxis2", "INTEGER"), ("pixscale_arcsec", "REAL"),
                      ("pixscale_source", "TEXT"), ("corner_ltt_s", "REAL")):
        if col not in have:
            con.execute(f"ALTER TABLE s3_header_audit ADD COLUMN "
                        f"{col} {decl}")
    con.commit()
    sample = pick_header_sample(con)
    # "Already measured" means measured WITH the current column set: a row
    # missing the geometry columns is re-read rather than skipped.
    done = {r[0] for r in con.execute(
        "SELECT path FROM s3_header_audit "
        "WHERE corner_ltt_s IS NOT NULL OR naxis1 IS NOT NULL").fetchall()}
    # Frames that fell out of the sample (e.g. calibration paths now
    # excluded from SCIENCE_WHERE) must not linger as evidence.
    keep = {s[2] for s in sample}
    if keep:                       # an empty IN () is a SQL syntax error,
        con.execute(               # and would also wipe a good table
            "DELETE FROM s3_header_audit WHERE path NOT IN (%s)"
            % ",".join("?" * len(keep)), list(keep))
        con.commit()
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
        # Field geometry, for the frame-center caveat: the stored BJD
        # points at the frame CENTER, and an object at a corner differs
        # by up to corner_ltt_s.  Measured from each sampled header's own
        # WCS/optics rather than hand-typed once for the whole archive.
        # TRUE geometry, not the header's NAXIS cards: for a tile-compressed
        # file whose header astropy could not fully translate, NAXIS1 is the
        # BINTABLE row length in bytes (8), and a pixel scale or a corner
        # light-time computed from an 8-pixel field is nonsense.  See
        # macro_core.fitsgeom and docs/pipeline/s0e_geometry_fix.html.
        g_n1, g_n2 = fitsgeom.resolve_geometry_or_none(hdr)
        scale, scale_source = tm.pixel_scale_arcsec(
            g_n1, g_n2,
            cdelt1=hdr.get("CDELT1"), cdelt2=hdr.get("CDELT2"),
            cd1_1=hdr.get("CD1_1"), cd1_2=hdr.get("CD1_2"),
            xpixsz=hdr.get("XPIXSZ"), focallen=hdr.get("FOCALLEN"),
            secpix1=hdr.get("SECPIX1"))
        corner_ltt = tm.field_corner_light_time_s(g_n1, g_n2, scale)
        con.execute("INSERT OR REPLACE INTO s3_header_audit "
                    "(path, family, year, era_id, night, exptime_s, "
                    " jd_header, date_obs, jd_minus_dateobs_s, "
                    " jd_helio_header, telut, telut_minus_dateobs_s, "
                    " ra_deg, dec_deg, helio_resid_start_s, "
                    " helio_resid_mid_s, swcreate, naxis1, naxis2, "
                    " pixscale_arcsec, pixscale_source, corner_ltt_s) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (path, family, year, era_id, night, exptime, jd_h,
                     date_obs_h, jd_minus_dateobs, jd_helio, telut,
                     telut_minus_dateobs, ra_deg, dec_deg,
                     resid_start, resid_mid,
                     str(hdr.get("SWCREATE") or ""),
                     g_n1, g_n2,
                     scale, scale_source, corner_ltt))
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
# Stage: drift — TELUT vs DATE-OBS across the whole baseline (relative clock)
# ---------------------------------------------------------------------------
def pick_drift_sample(con: sqlite3.Connection) -> list[tuple]:
    """Deterministic sample for the relative-drift check: up to
    :data:`DRIFT_SAMPLE_PER_ERA` frames from EVERY era.

    Sampling by era (the pinned registry) rather than by family/year is
    the point: a drift bound needs the same comparison repeated along the
    whole 2023-2026 baseline, not a stratified spread over configurations.
    Within an era the frames are taken at evenly spaced ranks of the
    path-ordered list, so the sample is reproducible and not clustered in
    one night.
    """
    picks: list[tuple] = []
    eras = con.execute(f"""
        SELECT era_id, count(*) FROM frames
        WHERE {SCIENCE_WHERE} AND jd IS NOT NULL AND exptime > 0
          AND date_obs IS NOT NULL AND date_obs != ''
        GROUP BY era_id ORDER BY era_id""").fetchall()
    for era_id, n in eras:
        # Evenly spaced ranks: for n frames and k picks, ranks
        # n*(2i+1)/(2k) put the samples at the centres of k equal blocks.
        k = min(DRIFT_SAMPLE_PER_ERA, n)
        for i in range(k):
            offset = (n * (2 * i + 1)) // (2 * k)
            row = con.execute(f"""
                SELECT path, era_id, night, exptime, jd, date_obs
                FROM frames WHERE {SCIENCE_WHERE} AND era_id = ?
                  AND jd IS NOT NULL AND exptime > 0
                  AND date_obs IS NOT NULL AND date_obs != ''
                ORDER BY path LIMIT 1 OFFSET ?""",
                (era_id, offset)).fetchone()
            if row is not None:
                picks.append(row)
    return picks


def stage_drift(con: sqlite3.Connection, archive: Path) -> None:
    """Sample TELUT vs DATE-OBS per era: the only two-independent-clock
    comparison the archive offers.

    DATE-OBS/JD are written by the acquisition software from the
    acquisition PC's clock.  TELUT is the TELESCOPE control system's UTC,
    read when the header is written — a genuinely different machine.  If
    the two clocks drifted apart over the 2023-2026 baseline, the
    residual

        TELUT - (DATE-OBS + EXPTIME)

    (i.e. what is left after the expected "header written one exposure
    later" offset) would trend with epoch.  Its epoch-to-epoch spread is
    therefore a REAL bound on relative clock behaviour — which the
    internal JD-vs-DATE-OBS agreement is not, because those two cards are
    written from the SAME clock and can only ever agree with each other.

    A row counts as INFORMATIVE only when its TELUT could physically be a
    header-write clock read: present, not the '1899-12-30' sentinel, not a
    verbatim copy of DATE-OBS, and no EARLIER than the end of the exposure
    it belongs to.  That last test is a validity check on the card, not a
    fit to the answer — the header cannot be written before the exposure
    finishes, so a TELUT that precedes the exposure end is a copy of
    something else.  It is what disqualifies the 2026 pyscope eras, whose
    TELUT sits 0 or exactly 1 s after DATE-OBS whether the exposure was
    0.25 s or 90 s: it does not track EXPTIME, so it is not a clock read.
    """
    con.execute("""CREATE TABLE IF NOT EXISTS s3_clock_drift (
        path TEXT PRIMARY KEY, era_id INTEGER, night TEXT,
        jd_utc_start REAL, exptime_s REAL, date_obs TEXT, telut TEXT,
        telut_minus_dateobs_s REAL, resid_s REAL, informative INTEGER,
        swcreate TEXT)""")
    sample = pick_drift_sample(con)
    done = {r[0] for r in con.execute(
        "SELECT path FROM s3_clock_drift").fetchall()}
    todo = [s for s in sample if s[0] not in done]
    log(f"drift: {len(sample)} sampled across eras, {len(todo)} to read")
    n_read = 0
    for path, era_id, night, exptime, jd, date_obs in todo:
        hdr = open_header(archive, path)
        if hdr is None:
            continue
        telut = str(hdr.get("TELUT") or "")
        date_obs_h = str(hdr.get("DATE-OBS") or date_obs)
        jd_from_date = tm.parse_date_obs(date_obs_h)
        telut_jd = tm.parse_date_obs(telut)
        delta = resid = None
        informative = 0
        if telut_jd is not None and jd_from_date is not None:
            delta = (telut_jd - jd_from_date) * 86400.0
            # The '1899-12-30' sentinel means "telescope clock not set".
            if abs(delta) <= 7200.0:
                exp_s = float(exptime or 0.0)
                resid = delta - exp_s
                # A TELUT byte-identical to DATE-OBS is a COPY, not a
                # second clock read.  So is a TELUT that lands before the
                # exposure could possibly have ended (allowing one whole
                # second for the coarsest stamping in the archive).
                informative = int(telut.strip() != date_obs_h.strip()
                                  and delta >= exp_s - 1.0)
            else:
                delta = None
        con.execute("INSERT OR REPLACE INTO s3_clock_drift VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?)",
                    (path, era_id, night, jd, exptime, date_obs_h, telut,
                     delta, resid, informative,
                     str(hdr.get("SWCREATE") or "")))
        n_read += 1
        if n_read % 25 == 0:
            con.commit()
            log(f"  ... {n_read}/{len(todo)} headers")
    con.commit()
    stats = con.execute("""
        SELECT count(*), count(DISTINCT era_id), min(resid_s), max(resid_s)
        FROM s3_clock_drift WHERE informative = 1
          AND resid_s IS NOT NULL""").fetchone()
    n_total = con.execute("SELECT count(*) FROM s3_clock_drift").fetchone()[0]
    write_meta(con, {"drift_rows": n_total,
                     "drift_informative_rows": stats[0],
                     "drift_eras": stats[1],
                     "drift_resid_min_s": stats[2],
                     "drift_resid_max_s": stats[3]})
    if stats[0]:
        log(f"drift: {n_total} rows ({stats[0]} informative across "
            f"{stats[1]} eras); TELUT-(DATE-OBS+EXPTIME) spans "
            f"{stats[2]:+.2f} to {stats[3]:+.2f} s")
    else:
        log(f"drift: {n_total} rows, none informative")


# ---------------------------------------------------------------------------
# Stage: frame-times — mid-exposure BJD_TDB for every canonical science frame
# ---------------------------------------------------------------------------
def families_with_start_evidence(con: sqlite3.Connection) -> set[str]:
    """Readout families for which convention 1 (header JD = exposure
    START) is actually PROVEN by the header audit.

    Two independent probes count as proof, and a family needs at least
    one of them:

    * a header JD-HELIO that matches our own heliocentric JD evaluated at
      start + EXPTIME/2 (so the base stamp must be the start);
    * a TELUT that is a genuinely different value from DATE-OBS and sits
      one exposure later (the ephemeris-free version of the same test).

    The 2026 pyscope eras satisfy NEITHER: they write no JD-HELIO card at
    all and copy DATE-OBS verbatim into TELUT.  Returning a set (rather
    than assuming universality) is what lets every frame_times row say on
    its face which case it is.

    Both probes are required to be DISCRIMINATING, not merely present:
    they separate "start" from "mid" by EXPTIME/2, so a 1-millisecond
    exposure proves nothing at all and must not be allowed to certify a
    family.  Hence the EXPTIME >= 2 s condition (a >= 1 s separation
    against ~0.15 s of residual noise).  This is what keeps the HDR
    family — whose only long-exposure JD-HELIO residuals are 836 s wide —
    out of the proven set on the strength of one 1 ms frame.
    """
    rows = con.execute("""
        SELECT family,
               sum(CASE WHEN helio_resid_mid_s IS NOT NULL
                         AND exptime_s >= 2.0
                         AND abs(helio_resid_mid_s) <= 1.0
                        THEN 1 ELSE 0 END),
               sum(CASE WHEN telut_minus_dateobs_s IS NOT NULL
                         AND exptime_s >= 2.0
                         AND telut_minus_dateobs_s >= exptime_s - 1.0
                        THEN 1 ELSE 0 END)
        FROM s3_header_audit GROUP BY family""").fetchall()
    return {fam for fam, n_helio, n_telut in rows if (n_helio or n_telut)}


def sibling_jd_drift(con: sqlite3.Connection, on_axis: set[str]
                     ) -> tuple[dict[str, float], int]:
    """Path -> |JD disagreement| (seconds) with the SAME exposure's other
    copy, plus the number of such pairs.

    Section 1's JD-vs-DATE-OBS audit can only catch a re-stamp that moved
    ONE card; when the reduction pipeline re-stamped BOTH, the copy stays
    internally consistent and looks perfect.  S0b already measured those
    cases by matching raw to reduced on the file stem, so the drift is
    known — it just never reached the time axis.  Both sides of a pair get
    the value, so a consumer can see it on whichever row they hold.

    ONLY pairs whose BOTH copies are on the axis are flagged.  Most links
    point at a reduced copy that S0 collapsed as a duplicate (not
    canonical, so not in ``frame_times``): there the raw row is the sole,
    authoritative stamp and a flag on it would be pure noise.  The pairs
    that matter are the ones where the same exposure genuinely appears
    TWICE on the shared axis.
    """
    out: dict[str, float] = {}
    n_pairs = 0
    for raw_path, reduced_path, drift_s in con.execute(
            "SELECT raw_path, reduced_path, jd_drift_s FROM "
            "raw_reduced_links WHERE jd_drift_s IS NOT NULL"):
        if raw_path not in on_axis or reduced_path not in on_axis:
            continue
        n_pairs += 1
        drift = abs(float(drift_s))
        # Keep the WORST disagreement if a path appears in several links.
        for p in (raw_path, reduced_path):
            if drift > out.get(p, -1.0):
                out[p] = drift
    return out, n_pairs


def stage_frame_times(con: sqlite3.Connection, ephemeris: str,
                      chunk: int = 25000) -> None:
    """Compute and write the ``frame_times`` table (atomic swap).

    One row per canonical science frame, keyed by path.  Frames without a
    JD get a row with method 'no_jd'; frames without coordinates get their
    UTC mid-time but a NULL BJD with method 'no_coords' — the row count of
    this table always equals the canonical-science row count, so nothing
    silently falls off the time axis.

    Two provenance columns ride along, because a time stamp that LOOKS
    fine is the dangerous kind:

    * ``start_evidence`` — whether this frame's readout family actually
      has header evidence that JD is the exposure start;
    * ``sibling_jd_drift_s`` — how far this exposure's OTHER copy
      (raw vs reduced) disagrees about when it happened.  Beyond
      :data:`tm.JD_SIBLING_DISAGREE_S` the reduced copy's BJD is
      WITHDRAWN (NULL, method
      :data:`tm.BJD_JD_DISAGREES`): two rows on the shared axis claiming
      the same photons at times 70 minutes apart is not a caveat, it is a
      defect, and the S0 rule says timing comes from the raw copy.
    """
    rows = con.execute(f"""
        SELECT path, obs_rowid, era_id, readoutm, jd, exptime,
               ra_deg, dec_deg
        FROM frames WHERE {SCIENCE_WHERE} ORDER BY path""").fetchall()
    n_calib_excluded = con.execute("""
        SELECT count(*) FROM frames
        WHERE is_canonical = 1
          AND (imagetyp IS NULL OR imagetyp LIKE 'Light%')
          AND path IN (SELECT path FROM calib_frames)""").fetchone()[0]
    log(f"frame-times: {len(rows):,} canonical science frames "
        f"({n_calib_excluded:,} header-mislabelled calibrations excluded)")
    proven = families_with_start_evidence(con)
    drifts, n_sibling_pairs = sibling_jd_drift(con, {r[0] for r in rows})
    # Which side of a re-stamped pair loses its BJD: the REDUCED copy.
    reduced_paths = {r[0] for r in con.execute(
        "SELECT reduced_path FROM raw_reduced_links")}
    out: list[list] = []
    # First pass: the pure mid-time policy per frame.
    mids, ras, decs, computable = [], [], [], []
    n_unverified = 0
    for i, (path, obs_rowid, era_id, readoutm, jd, exptime,
            ra_deg, dec_deg) in enumerate(rows):
        mid, method = tm.jd_utc_mid(jd, exptime, readoutm)
        family = (readoutm or "").strip() or "(blank)"
        evidence = (tm.START_VERIFIED if family in proven
                    else tm.START_UNVERIFIED)
        n_unverified += evidence == tm.START_UNVERIFIED
        out.append([path, obs_rowid, era_id, jd, exptime, mid,
                    None, None, None, ra_deg, dec_deg, method,
                    None, evidence, drifts.get(path), tm.S3_CODE_VERSION])
        if mid is not None and ra_deg is not None and dec_deg is not None:
            computable.append(i)
            mids.append(mid)
            ras.append(ra_deg)
            decs.append(dec_deg)
        elif mid is not None:
            out[i][12] = tm.BJD_NO_COORDS
        else:
            out[i][12] = tm.BJD_NO_JD
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
    # Third pass: withdraw the BJD of any REDUCED copy whose raw parent
    # disagrees about the epoch by more than the stated threshold.
    withdrawn: list[tuple] = []
    for r in out:
        drift = r[14]
        if drift is not None and drift > tm.JD_SIBLING_DISAGREE_S \
                and r[0] in reduced_paths and r[6] is not None:
            withdrawn.append((r[0], r[2], r[3], float(drift), r[6]))
            r[6] = r[7] = r[8] = None
            r[12] = tm.BJD_JD_DISAGREES
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
            start_evidence TEXT,   -- is 'JD = start' proven for this era?
            sibling_jd_drift_s REAL, -- raw-vs-reduced stamp disagreement
            code_version TEXT)""",
               [tuple(r) for r in out],
               "INSERT INTO {table} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)")
    swap_table(con, "s3_time_outliers", """
        CREATE TABLE {table} (
            path TEXT PRIMARY KEY, era_id INTEGER, jd_utc_start REAL,
            sibling_jd_drift_s REAL, withdrawn_bjd_tdb REAL)""",
               withdrawn, "INSERT INTO {table} VALUES (?,?,?,?,?)")
    n_bjd = sum(1 for r in out if r[6] is not None)
    n_flagged = sum(1 for r in out
                    if r[14] is not None
                    and r[14] > tm.JD_SIBLING_DISAGREE_S)
    write_meta(con, {"frame_times_rows": len(out),
                     "frame_times_with_bjd": n_bjd,
                     "frame_times_calib_excluded": n_calib_excluded,
                     "frame_times_unverified_start": n_unverified,
                     "frame_times_sibling_pairs": n_sibling_pairs,
                     "frame_times_sibling_flagged": n_flagged,
                     "frame_times_bjd_withdrawn": len(withdrawn),
                     "ephemeris": ephemeris})
    log(f"frame-times: wrote {len(out):,} rows ({n_bjd:,} with BJD_TDB); "
        f"{n_unverified:,} rows with unverified start semantics; "
        f"{len(withdrawn)} BJDs withdrawn for raw/reduced JD disagreement")


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
               if abs(ph) > tm.CLOCK_OOE_PHASE]
        if len(ooe) < 3:
            continue                      # no baseline -> group unusable
        base = float(np.median(ooe))
        for night, hjd, ph, dm, de in members:
            pts.append((night, config, filt, hjd, ph, dm - base, de))
    fit_pts = [p for p in pts if abs(p[4]) <= tm.CLOCK_FIT_PHASE]
    results: list[tuple] = []
    # Nights whose coverage actually brackets the eclipse: only these may
    # enter the global fit (see the gate inside summarize).
    good_nights: set[str] = set()

    def summarize(tag: str, sel_nights) -> None:
        """Fit one selection of points and append its row.

        A row is ALWAYS appended — even when the fit is refused — so that
        a configured night can never disappear from the table with no
        trace, and so the report's text can be derived from the table
        instead of from the hand-written night tuple.
        """
        # ``sel_nights is None`` means "every point"; an EMPTY set means
        # "no night qualified" and must select nothing, not everything.
        sel = fit_pts if sel_nights is None else \
            [p for p in fit_pts if p[0] in sel_nights]
        if len(sel) < tm.CLOCK_MIN_FIT_POINTS:
            results.append((tag, len(sel),
                            float(min(p[4] for p in sel)) if sel else None,
                            float(max(p[4] for p in sel)) if sel else None,
                            None, None, None, None, None,
                            tm.CLOCK_STATUS_TOO_FEW))
            return
        ph_a = np.array([p[4] for p in sel])
        # COVERAGE GATE.  A symmetric template fitted to a one-sided arc
        # converges happily and returns a confident centre that is really
        # the slope of the flank it was given — nothing about the
        # eclipse's midpoint.  Demand real sampling on BOTH sides of
        # phase zero before calling the result a measurement.
        n_before, n_after = tm.phase_coverage(ph_a)
        if min(n_before, n_after) < tm.CLOCK_MIN_SIDE_POINTS:
            results.append((tag, len(sel), float(ph_a.min()),
                            float(ph_a.max()), None, None, None, None,
                            None, tm.CLOCK_STATUS_ONE_SIDED))
            return
        f = tm.fit_eclipse_offset(ph_a,
                                  np.array([p[5] for p in sel]),
                                  np.array([max(p[6], 0.01) for p in sel]))
        if f["ph0"] is None:
            results.append((tag, len(sel), float(ph_a.min()),
                            float(ph_a.max()), None, None, None, None,
                            None, tm.CLOCK_STATUS_NO_DIP))
            return
        oc_s = f["ph0"] * eph["period_d"] * 86400.0
        oc_err_s = f["ph0_err"] * eph["period_d"] * 86400.0
        # Cycle count at the epoch of these points, for the drift term.
        mean_hjd = float(np.mean([p[3] for p in sel]))
        cycles = abs(mean_hjd - eph["epoch_hjd"]) / eph["period_d"]
        # THE EPHEMERIS TERM.  clock_error = (O-C) - (ephemeris error), so
        # the bound can never be tighter than how well the star's own
        # ephemeris is known.  VSX's quoted last digit (0.5e-7 d/cycle) is
        # a TYPOGRAPHIC precision, not a measurement uncertainty; the Gaia
        # DR3 solution for the same star puts the real period uncertainty
        # at 3.0e-5 d, ~600x larger.  Taking the max of the two keeps the
        # printed bound honest instead of advertising a ceiling the
        # ephemeris cannot support.
        vsx_term_s = EPOCH_QUANT_S + PERIOD_QUANT_D * cycles * 86400.0
        gaia_term_s = EPOCH_QUANT_S + \
            GAIA_EB["period_err_d"] * cycles * 86400.0
        eph_sys_s = max(vsx_term_s, gaia_term_s)
        bound_s = abs(oc_s) + oc_err_s + eph_sys_s
        results.append((tag, len(sel), float(ph_a.min()),
                        float(ph_a.max()), f["depth"], f["width"],
                        oc_s, oc_err_s, bound_s, tm.CLOCK_STATUS_OK))
        if sel_nights:
            good_nights.update(sel_nights)

    # Per-night rows FIRST: the global fit may only use nights that
    # passed the coverage gate.  Folding a one-sided night in with a good
    # one does not average two measurements, it contaminates one.
    for night in CLOCK_ECLIPSE_NIGHTS:
        summarize(night, {night})
    summarize("global", good_nights)
    fit_pts = [p for p in fit_pts if p[0] in good_nights]
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
        # What the tight-but-wrong version of the bound would have been:
        # kept so the report can show the difference between propagating
        # VSX's quoted last digit and propagating a credible period error.
        "clock_vsx_quant_term_s":
            f"{EPOCH_QUANT_S + PERIOD_QUANT_D * cycles_all * 86400.0:.0f}",
        "clock_eph_term_s":
            f"{EPOCH_QUANT_S + GAIA_EB['period_err_d'] * cycles_all * 86400.0:.0f}",
        "clock_nights_gated": ";".join(sorted(good_nights)) or "none",
        "clock_nights_configured": ";".join(CLOCK_ECLIPSE_NIGHTS),
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
STAGES = ("audit-scan", "cadence", "audit-headers", "drift", "frame-times",
          "clock", "report")

#: Other MACRO stages (the S1 plate-solve batch, above all) write to this
#: same database concurrently.  Five minutes of patience on a locked
#: writer is cheap; a crashed build in the middle of a long night is not.
BUSY_TIMEOUT_MS = 300_000


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
        # Wait out a concurrent writer instead of failing on it.
        con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
        if "audit-scan" in stages:
            stage_audit_scan(con)
        if "cadence" in stages:
            stage_cadence(con)
        if "audit-headers" in stages:
            stage_audit_headers(con, args.archive, ephemeris)
        if "drift" in stages:
            stage_drift(con, args.archive)
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
