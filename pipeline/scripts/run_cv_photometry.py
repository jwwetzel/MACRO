#!/usr/bin/env python
"""Production ensemble photometry for the five staged CV targets.

WHAT THIS SCRIPT DOES
---------------------
Reads the CV staging manifest (``stage_cv_timeseries`` in the S0 manifest),
measures every staged science frame of ST LMi, VV Pup, AN UMa, EU UMa and
YZ Cnc, registers the frames to a per-(target, era) reference — through the
S1 plate solutions where they exist, through astroalign triangles where
they do not — solves a Honeycutt (1992) inhomogeneous ensemble SEPARATELY
for every (target, era, filter) series, and records the diagnostics that
say how far each result may be trusted.  Everything lands in one database:

    products/phot/cv_timeseries.sqlite

This is a PRODUCTION run of the S4 machinery in ``pipeline/macro_phot``,
not a new photometry code.  The aperture arithmetic, the reference-frame
double-image QC, the Honeycutt solver, the comparison-star stability
iteration and the error model are all imported unchanged.  What this script
adds is what production scale demands and the prototype never faced:

* **Series separation.**  Every row in every table carries a ``series_key``
  of the form ``target|eNN|filter``.  Era mixing inside one light curve is
  impossible by construction, not by convention.
* **Provenance per era.**  The staging manifest stages RAW frames.  Some
  eras have a complete server-reduced counterpart, some have none.  One
  provenance is chosen per (target, era) and applied to all of it
  (``macro_phot.series.choose_provenance``); eras with no reduced tree are
  calibrated locally from the staged era-matched masters.
* **The S2 saturation vetoes, per readout mode.**  High Gain clips at
  3,496 ADU and Mode0 at 65,535 — a single global threshold cannot serve
  both, and the prototype's 55,000 vetoed literally nothing on the entire
  2024 High Gain season.  Saturated measurements are FLAGGED and withheld
  from the ensemble; they are still written to the light curve so the next
  analyst can see them rather than wonder where they went.
* **Resumability and parallelism.**  Extraction and matching are chunked,
  committed per batch, and safe to re-run; a killed run resumes where it
  stopped.  Worker processes are capped (``--workers``, hard maximum 6)
  because this machine is also running an S1 batch solve.

WHAT IT DELIBERATELY DOES NOT DO
--------------------------------
It does not interpret a single light curve.  No period, no eclipse, no
state change, no colour.  The deliverable is a measurement set with its
diagnostics attached; the science verdicts belong to the next stage under
a different discipline.

TABLES WRITTEN
--------------
``cv_series``       one row per (target, era, filter): provenance, frame
                    census, match rate, ensemble convergence, comparison
                    and check-star counts, check RMS, chi2 inflation.
``cv_frames``       one row per staged science frame: BJD_TDB, provenance,
                    calibration recipe, extraction stats, registration
                    method and QC, the frame's ensemble zero point.
``cv_detections``   every source detection, its aperture flux, its
                    saturation flag, and its star identity once matched.
``cv_ref``          the reference frame per (target, era) + its QC.
``cv_ref_stars``    the per-(target, era) star catalog, with sky positions
                    and Gaia identities where available.
``cv_field_tie``    how the target was identified in each (target, era).
``cv_stars``        per (series, star): role, mean magnitude, RMS, chi2.
``cv_zeropoints``   the ensemble zero-point timeline, per series per frame.
``cv_lightcurve``   the deliverable: every measurement of the target, the
                    comparison stars and the check stars, in BJD_TDB.
``cv_error_model``  per series: check-star RMS floor, median chi2, the
                    error inflation factor.
``cv_allan``        Allan deviation ladder of one long night of one check
                    star, per target.
``cv_selection``    the frame ledger: staged vs usable vs excluded, with
                    the reason for every exclusion.
``cv_build_meta``   timestamp, code version, git commit, every constant.

STAGES (run in this order; each is resumable and safe to repeat)
----------------------------------------------------------------
    init                     build the series registry and frame worklist
    extract [--limit N]      measure N pending frames (repeat until 0)
    match   [--limit N]      choose references; register + match N frames
    field                    identify the target in each (target, era)
    ensemble                 solve every (target, era, filter) ensemble
    errors                   error model + Allan ladders
    status                   progress summary (safe at any time)

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_photometry.py init
    $P pipeline/scripts/run_cv_photometry.py extract --limit 1200 --workers 6
    $P pipeline/scripts/run_cv_photometry.py match   --limit 1200 --workers 6
    $P pipeline/scripts/run_cv_photometry.py field
    $P pipeline/scripts/run_cv_photometry.py ensemble
    $P pipeline/scripts/run_cv_photometry.py errors
    $P pipeline/scripts/run_cv_photometry.py status
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Make the pipeline packages importable however the script is invoked.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import S4_CODE_VERSION                        # noqa: E402
from macro_phot import calib as cal                           # noqa: E402
from macro_phot import ensemble as ens                        # noqa: E402
from macro_phot import errors as err                          # noqa: E402
from macro_phot import extract as ext                         # noqa: E402
from macro_phot import gaia as gg                             # noqa: E402
from macro_phot import photometry as ph                       # noqa: E402
from macro_phot import register as reg                        # noqa: E402
from macro_phot import series as sr                           # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare commands Just Work)
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
DEFAULT_WCS_ROOT = REPO_ROOT / "products" / "astrom" / "wcs"
DEFAULT_RECON_DIR = REPO_ROOT / "products" / "detector" / "recon"

#: The campaign's code version, recorded into cv_build_meta.
CV_CODE_VERSION = "CV-S4 v1.0 (2026-08-18)"

#: Concurrency ceiling.  The machine is simultaneously running a 10-worker
#: S1 batch solve and a large rclone transfer; this build stays a polite
#: guest.  --workers above this is silently clamped, and the clamp is
#: recorded in the build metadata so a fast-looking run cannot later be
#: mistaken for an unthrottled one.
MAX_WORKERS = 6

#: The manifest is shared with a running S1 batch.  Five minutes of
#: patience beats a spurious 'database is locked' in the middle of a
#: two-hour campaign.
BUSY_TIMEOUT_MS = 300_000

#: How many frames one worker batch carries before the parent commits.
#: Small enough that a kill loses seconds of work, large enough that the
#: process pool is not re-primed constantly.
WRITE_BATCH = 100


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def connect(db_path: Path) -> sqlite3.Connection:
    """Open the PRODUCTS database (this build's own; we are its only writer)."""
    con = sqlite3.connect(db_path, timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    # Durability is traded for speed deliberately: every table here is
    # REGENERABLE from the archive by re-running the stage, so the worst
    # a power cut can cost is one batch of work, not a scientific result.
    con.execute("PRAGMA synchronous=NORMAL")
    return con


def connect_manifest(path: Path) -> sqlite3.Connection:
    """Open the S0 manifest READ-ONLY, patiently (a batch solve is writing)."""
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True,
                          timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    return con


def git_commit() -> str:
    """HEAD's short hash with an honest '-dirty' marker (S4 discipline)."""
    try:
        head = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              cwd=REPO_ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT,
                               capture_output=True, text=True,
                               check=True).stdout.strip()
        return f"{head}-dirty" if dirty else head
    except Exception:
        return ""


def meta_write(con: sqlite3.Connection, extra: dict = ()) -> None:
    """Upsert build metadata; every stage that writes keeps it current."""
    con.execute("""CREATE TABLE IF NOT EXISTS cv_build_meta
                   (key TEXT PRIMARY KEY, value TEXT)""")
    rows = {
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "cv_code_version": CV_CODE_VERSION,
        "s4_code_version": S4_CODE_VERSION,
        "git_commit": git_commit(),
        "aperture_radius_arcsec": str(ph.APERTURE_RADIUS_ARCSEC),
        "sky_annulus_arcsec": str(ph.SKY_ANNULUS_ARCSEC),
        "detect_sigma": str(ph.DETECT_SIGMA),
        "clip_sigma": str(ens.CLIP_SIGMA),
        "weight_floor_mag": str(ens.WEIGHT_FLOOR_MAG),
        "min_comps": str(ens.MIN_COMPS),
        "n_check_stars": str(ens.N_CHECK_STARS),
        "min_series_frames": str(sr.MIN_SERIES_FRAMES),
        "gain_nominal_e_per_adu": str(sr.NOMINAL_GAIN_E_PER_ADU),
        "gain_bracket_e_per_adu": str(sr.GAIN_BRACKET_E_PER_ADU),
        "s2_mode_veto_adu": json.dumps(sr.S2_MODE_VETO_ADU),
        "s2_mode_ceiling_adu": json.dumps(sr.S2_MODE_CEILING_ADU),
        "timing_source": "frame_times.bjd_tdb (S3); header JD-HELIO refused",
        "provenance_rule": ("one pixel provenance per (target, era); "
                            "never mixed inside a series"),
        "max_workers": str(MAX_WORKERS),
        **dict(extra),
    }
    con.executemany("INSERT OR REPLACE INTO cv_build_meta VALUES (?,?)",
                    list(rows.items()))
    con.commit()


def swap_in(con: sqlite3.Connection, name: str, create_sql: str,
            rows: list[tuple]) -> None:
    """Atomic-swap rebuild of one derived table (the S0b discipline)."""
    tmp = f"{name}__new"
    con.execute(f"DROP TABLE IF EXISTS {tmp}")
    con.execute(create_sql.format(t=tmp))
    if rows:
        width = len(rows[0])
        con.executemany(
            f"INSERT INTO {tmp} VALUES ({','.join('?' * width)})", rows)
    con.commit()
    con.execute("BEGIN")
    con.execute(f"DROP TABLE IF EXISTS {name}")
    con.execute(f"ALTER TABLE {tmp} RENAME TO {name}")
    con.commit()


def fnum(x) -> float | None:
    """float(x) when finite, else None — the SQL-safe cast used everywhere.

    NaN survives a round trip through SQLite as NULL on the way in but as
    a float on the way out in some drivers; forcing it to None at the
    boundary keeps 'no value' a single, unambiguous thing.
    """
    if x is None:
        return None
    v = float(x)
    return v if np.isfinite(v) else None


def clamp_workers(n: int) -> int:
    """Requested workers, clamped to the concurrency ceiling (>= 1)."""
    return max(1, min(int(n), MAX_WORKERS))


# ---------------------------------------------------------------------------
# Stage: init
# ---------------------------------------------------------------------------
def _recon_transform(recon_dir: Path, era_id: int) -> tuple | None:
    """The S2-measured (flat median, dark median, pedestal) of one era.

    ``products/detector/recon/eraNN.npz`` holds the S2 reconstruction
    experiment's fitted flat image F, dark image D and per-pair pedestal
    for the server reduction of that era.  Their medians turn a RAW
    saturation veto into the reduced frame's own ADU scale
    (:func:`macro_phot.series.veto_in_reduced_adu`).  Returns None when the
    era was never reconstructed, and the caller then records that the raw
    threshold was used unmapped.
    """
    path = Path(recon_dir) / f"era{int(era_id)}.npz"
    if not path.exists():
        return None
    try:
        z = np.load(path, allow_pickle=True)
        f_med = float(np.nanmedian(z["F"]))
        d_med = float(np.nanmedian(z["D"]))
        ped = float(np.nanmedian(z["pedestal"]))
        if not (np.isfinite(f_med) and f_med > 0):
            return None
        return f_med, d_med, ped
    except Exception:
        return None


def _staged_masters(mcon: sqlite3.Connection) -> dict:
    """Every staged master calibration, grouped by (era, kind).

    Returned as ``{(era_id, kind): [(jd, abs_path, filter, exptime), ...]}``
    where kind is 'dark' or 'flat'.  ``master_bias`` frames are not used:
    the archived master DARKS of these eras are bias-inclusive (they are
    the 'master_dark' of a bias+dark stack, which is why they exist at one
    file per exposure time), so subtracting a bias as well would remove it
    twice.
    """
    out: dict = {}
    rows = mcon.execute("""
        SELECT role, CAST(era_id AS INT), jd, abs_path, filter, exptime
        FROM stage_cv_timeseries
        WHERE role IN ('master_dark', 'master_flat')""").fetchall()
    for role, era, jd, abs_path, filt, expt in rows:
        kind = "dark" if role == "master_dark" else "flat"
        out.setdefault((era, kind), []).append((jd, abs_path, filt, expt))
    return out


def _pick_frame_masters(masters: dict, era: int, filt: str | None,
                        jd: float | None, exptime: float | None
                        ) -> tuple[str | None, str | None]:
    """The master dark and flat that serve ONE science frame.

    The dark must match the science exposure time (no scaling — these
    masters include the bias, and scaling a bias-inclusive dark corrupts
    the bias); the flat must match the filter exactly, case-sensitively,
    because 'g' and 'G' are different filters in different eras and
    quietly folding them together would flat-field with the wrong band.
    Ties in either list are broken by nearest observation date, then by
    path, so the choice is reproducible.
    """
    darks = [(c[0], c[1]) for c in masters.get((era, "dark"), [])
             if sr.dark_exptime_matches(c[3], exptime)]
    flats = [(c[0], c[1]) for c in masters.get((era, "flat"), [])
             if filt is not None and c[2] == filt]
    d = sr.pick_master(darks, jd)
    f = sr.pick_master(flats, jd)
    return (d[1] if d else None), (f[1] if f else None)


def cmd_init(args) -> None:
    """Build the series registry and the frame worklist. Idempotent.

    Every staged science row becomes a ``cv_frames`` row whose columns
    already encode the decisions this campaign is accountable for: which
    series it belongs to, which pixels it may read, which saturation veto
    applies to those pixels, whether its geometry admits aperture
    photometry at all, and its BJD_TDB from S3 (never the header's
    heliocentric UTC JD).  Re-running preserves work already done:
    measurements are keyed on frame_id and inserted with OR IGNORE.
    """
    args.db.parent.mkdir(parents=True, exist_ok=True)
    mcon = connect_manifest(args.manifest)
    con = connect(args.db)

    con.execute("""CREATE TABLE IF NOT EXISTS cv_frames (
        frame_id INTEGER PRIMARY KEY,        -- obs_rowid in the S0 manifest
        series_key TEXT, target_key TEXT, era_id INTEGER, filter TEXT,
        night TEXT, bjd_tdb REAL, jd_header REAL, exptime REAL, airmass REAL,
        readoutm TEXT,
        provenance TEXT,                     -- server_reduced/local_master/raw
        pixel_path TEXT,                     -- archive-relative, as measured
        raw_path TEXT,                       -- always the staged raw frame
        master_dark TEXT, master_flat TEXT, calib_recipe TEXT,
        veto_adu REAL,                       -- saturation veto in THESE units
        veto_basis TEXT,                     -- how that number was obtained
        qc_flags TEXT, pointing_offset_deg REAL,
        naxis1 INTEGER, naxis2 INTEGER,
        status TEXT DEFAULT 'pending',
        exclude_reason TEXT,
        plate_scale REAL, scale_basis TEXT, aper_px REAL,
        n_detected INTEGER, n_saturated INTEGER, n_nonfinite INTEGER,
        bkg_adu REAL, bkg_rms REAL, fwhm_px REAL,
        has_wcs INTEGER,
        reg_method TEXT, ali_nmatch INTEGER, ali_rms_px REAL,
        ali_scale REAL, ali_rot_deg REAL, match_rate REAL,
        zp REAL, zp_err REAL, n_star_used INTEGER)""")
    con.execute("""CREATE INDEX IF NOT EXISTS idx_frames_series
                   ON cv_frames (series_key, status)""")
    con.execute("""CREATE INDEX IF NOT EXISTS idx_frames_era
                   ON cv_frames (target_key, era_id, status)""")
    con.execute("""CREATE TABLE IF NOT EXISTS cv_detections (
        frame_id INTEGER, det_id INTEGER,
        x REAL, y REAL, flux REAL, fluxerr REAL, fwhm REAL, peak REAL,
        sep_flag INTEGER, saturated INTEGER,
        star_id INTEGER,
        PRIMARY KEY (frame_id, det_id))""")
    con.execute("""CREATE INDEX IF NOT EXISTS idx_det_star
                   ON cv_detections (star_id)""")
    # Idempotent schema top-up so init stays the ONE place the schema is
    # declared: a database created by an earlier run of this script has no
    # n_nonfinite column, and re-creating the table would throw away
    # measurements that are still perfectly good.
    for ddl in ("ALTER TABLE cv_frames ADD COLUMN n_nonfinite INTEGER",
                "ALTER TABLE cv_frames ADD COLUMN scale_basis TEXT"):
        try:
            con.execute(ddl)
        except sqlite3.OperationalError:
            pass          # column already present — the normal case

    # ---- era facts: readout mode (for the veto) and the S2 reduction map
    eras = {int(r[0]): r[1] for r in
            mcon.execute("SELECT era_id, readoutm FROM eras").fetchall()}
    masters = _staged_masters(mcon)

    # ---- the staged science rows, with BJD_TDB and the reduced link
    rows = mcon.execute("""
        SELECT s.obs_rowid, s.target_key, CAST(s.era_id AS INT), s.filter,
               s.night, s.jd, s.exptime, s.path, s.qc_flags,
               s.pointing_offset_deg,
               l.reduced_path, t.bjd_tdb, f.airmass, f.naxis1, f.naxis2
        FROM stage_cv_timeseries s
        LEFT JOIN raw_reduced_links l ON l.raw_rowid = s.obs_rowid
        LEFT JOIN frame_times      t ON t.obs_rowid = s.obs_rowid
        LEFT JOIN frames           f ON f.obs_rowid = s.obs_rowid
        WHERE s.role = 'science'
        ORDER BY s.obs_rowid""").fetchall()
    print(f"init: {len(rows)} staged science rows", flush=True)

    # ---- provenance is decided per (target, era), never per frame
    per_era: dict = {}
    for r in rows:
        key = (r[1], r[2])
        n, nl = per_era.get(key, (0, 0))
        per_era[key] = (n + 1, nl + (1 if r[10] else 0))
    provenance: dict = {}
    for key, (n, nl) in per_era.items():
        era = key[1]
        has_masters = bool(masters.get((era, "dark")) or
                           masters.get((era, "flat")))
        provenance[key] = sr.choose_provenance(n, nl, has_masters)

    # ---- the saturation veto, expressed in the units actually measured
    veto_for: dict = {}
    for key, (prov, _why) in provenance.items():
        era = key[1]
        mode = eras.get(era)
        raw_veto = sr.veto_adu(mode)
        if raw_veto is None:
            veto_for[key] = (None, f"no S2 ceiling measured for mode {mode!r}")
        elif prov == "server_reduced":
            tr = _recon_transform(args.recon_dir, era)
            if tr is None:
                veto_for[key] = (float(raw_veto), (
                    f"S2 raw veto {raw_veto} ADU for {mode!r}; era not "
                    f"reconstructed, so the reduced-frame offset is UNKNOWN "
                    f"and the raw threshold is used unmapped"))
            else:
                f_med, d_med, ped = tr
                mapped = sr.veto_in_reduced_adu(raw_veto, f_med, d_med, ped)
                veto_for[key] = (mapped, (
                    f"S2 raw veto {raw_veto} ADU for {mode!r} mapped through "
                    f"the measured era-{era} reduction (F={f_med:.4f}, "
                    f"D={d_med:.1f} ADU, pedestal={ped:.0f})"))
        else:
            veto_for[key] = (float(raw_veto), (
                f"S2 raw veto {raw_veto} ADU for {mode!r}, applied to raw "
                f"pixels (the local dark median is subtracted from it at "
                f"extraction time when a master dark is applied)"))

    n_new = 0
    ledger: dict = {}
    con.execute("BEGIN")
    for (fid, tk, era, filt, night, jd, expt, raw_path, qcf, poff,
         red_path, bjd, airmass, nax1, nax2) in rows:
        key = (tk, era)
        prov, _why = provenance[key]
        skey = sr.series_key(tk, era, filt)
        veto, veto_basis = veto_for[key]

        # Which file's pixels this frame contributes, and the local recipe.
        mdark = mflat = None
        recipe = None
        if prov == "server_reduced":
            pixel_path = red_path
        else:
            pixel_path = raw_path
            if prov == "local_master":
                mdark, mflat = _pick_frame_masters(masters, era, filt, jd,
                                                   expt)
                recipe = ("dark+flat" if (mdark and mflat) else
                          "dark_only" if mdark else
                          "flat_only" if mflat else "none")

        # Admission checks that need no pixels: geometry and provenance.
        scale = None
        aper = None
        status, reason = "pending", None
        if pixel_path is None:
            status, reason = "excluded", (
                "no pixel file for this provenance: the frame has no "
                "server-reduced counterpart while its era uses the reduced "
                "tree, and mixing raw pixels into the series is forbidden")
        # NOTE (2026-08-18): there is deliberately NO geometry test here.
        # An earlier version of this stage rejected frames whose manifest
        # NAXIS said the image was too small to hold an aperture, and it
        # threw out all 207 EU UMa era-80 frames as '8-pixel readout
        # strips'.  They are nothing of the kind: opening the files shows
        # 4,800 x 3,211 raw and 4,787 x 3,193 reduced images.  The 8 is the
        # BINTABLE ROW LENGTH of a tile-compressed header that the S0 scan
        # read without translating — the same artifact macro_core.fitsgeom
        # exists to resolve.  Geometry is now judged ONLY at extraction,
        # against the pixels themselves (macro_phot.series.geometry_verdict
        # on the resolved header), because a decision this destructive must
        # never rest on a second-hand number.
        if bjd is None and status == "pending":
            status, reason = "excluded", (
                "no BJD_TDB in frame_times: the timing rule forbids "
                "substituting the header's heliocentric UTC JD")

        cur = con.execute("""INSERT OR IGNORE INTO cv_frames
            (frame_id, series_key, target_key, era_id, filter, night,
             bjd_tdb, jd_header, exptime, airmass, readoutm, provenance,
             pixel_path, raw_path, master_dark, master_flat, calib_recipe,
             veto_adu, veto_basis, qc_flags, pointing_offset_deg,
             naxis1, naxis2, status, exclude_reason)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (fid, skey, tk, era, filt, night, bjd, jd, expt, airmass,
             eras.get(era), prov, pixel_path, raw_path, mdark, mflat, recipe,
             veto, veto_basis, qcf, poff,
             int(nax1) if nax1 else None, int(nax2) if nax2 else None,
             status, reason))
        n_new += cur.rowcount
        if cur.rowcount == 0:
            # Resumed build: refresh the DECISION columns (they may have
            # improved upstream) without touching measurement columns or
            # the status of already-processed frames.
            con.execute("""UPDATE cv_frames SET bjd_tdb=?, provenance=?,
                pixel_path=?, master_dark=?, master_flat=?, calib_recipe=?,
                veto_adu=?, veto_basis=?, qc_flags=?, series_key=?
                WHERE frame_id=?""",
                (bjd, prov, pixel_path, mdark, mflat, recipe, veto,
                 veto_basis, qcf, skey, fid))
        led = ledger.setdefault(skey, {"staged": 0, "excluded": 0})
        led["staged"] += 1
        if status == "excluded":
            led["excluded"] += 1
    con.commit()

    # ---- the series registry and the selection ledger
    series_rows = []
    for skey, led in sorted(ledger.items()):
        tk, era, filt = sr.parse_series_key(skey)
        prov, why = provenance[(tk, era)]
        series_rows.append((skey, tk, era, filt, eras.get(era), prov, why,
                            led["staged"], led["excluded"]))
    swap_in(con, "cv_selection",
            """CREATE TABLE {t} (series_key TEXT PRIMARY KEY,
               target_key TEXT, era_id INTEGER, filter TEXT, readoutm TEXT,
               provenance TEXT, provenance_reason TEXT,
               n_staged INTEGER, n_excluded_at_init INTEGER)""", series_rows)
    meta_write(con, {"manifest_path": str(args.manifest),
                     "archive_root": str(args.archive)})
    total = con.execute("SELECT count(*) FROM cv_frames").fetchone()[0]
    pend = con.execute("SELECT count(*) FROM cv_frames "
                       "WHERE status='pending'").fetchone()[0]
    print(f"init: worklist holds {total} frames ({n_new} new); "
          f"{pend} pending, {total - pend} excluded up front", flush=True)
    for skey, tk, era, filt, mode, prov, why, ns, nx in series_rows:
        if nx:
            print(f"  {skey:22s} {prov:15s} {nx}/{ns} excluded at init",
                  flush=True)
    con.close()
    mcon.close()


# ---------------------------------------------------------------------------
# Stage: extract (parallel, chunked, resumable)
# ---------------------------------------------------------------------------
def _extract_one(job: dict) -> dict:
    """Measure ONE frame.  Runs in a worker process; touches no database.

    Reads the frame's pixels (server-reduced, or raw plus this era's local
    master recipe), converts the fixed 4-arcsec sky aperture to pixels
    using THIS frame's own plate-scale cards, detects sources and measures
    aperture fluxes, and flags every detection at or above the readout
    mode's S2 saturation veto.  Returns a plain dict — the parent writes.
    """
    out = {"frame_id": job["frame_id"]}
    try:
        data, meta = ext.read_reduced(Path(job["pixel_path"]))
        recipe = "server_reduced"
        veto = job["veto_adu"]
        if job["provenance"] != "server_reduced":
            dark = (cal.read_master(Path(job["master_dark"]))
                    if job.get("master_dark") else None)
            flat = (cal.read_master(Path(job["master_flat"]))
                    if job.get("master_flat") else None)
            data, recipe = cal.apply_masters(data, dark, flat)
            if dark is not None and veto is not None:
                # The veto is a level on the RAW pixels; the pixels being
                # measured have had the dark removed, so the threshold
                # moves with them.  (The flat is median-normalized, so it
                # does not shift levels.)
                veto = float(veto) - float(np.median(dark))
        scale, scale_basis = sr.resolve_plate_scale(meta)
        aper = ph.aperture_radius_px(scale)
        if aper is None:
            raise ValueError("no plate scale: neither XPIXSZ/FOCALLEN nor a "
                             "CD matrix nor CDELT is usable in this header")
        ok, why = sr.geometry_verdict(meta["naxis1"], meta["naxis2"], aper)
        if not ok:
            raise ValueError(why)
        # The gain handed to sep is the S2 NOMINAL value, never the header
        # card: EGAIN reads 0.247 on Mode0 and 0 on the iKon, both of which
        # would silently corrupt the photon-shot term of every error bar.
        stats, dets = ext.measure_frame(data, sr.NOMINAL_GAIN_E_PER_ADU, aper)
        sat = sr.saturated_mask(dets["peak"], stats["bkg_adu"], veto)
        # Drop detections with no measurable position or flux.  Local
        # flat-fielding sets DEAD flat pixels to NaN on purpose (dividing
        # by them would manufacture a spectacular fake star), and any
        # aperture overlapping one returns a NaN flux.  SQLite stores a NaN
        # float as NULL, so such a row later reads back as "a star with no
        # brightness" and, sorted into the astroalign bright pool, crashed
        # the matcher outright — the honest fix is to never record a
        # measurement that does not exist.  The count of what was dropped
        # is kept (n_nonfinite) rather than silently absorbed.
        keep = (np.isfinite(dets["x"]) & np.isfinite(dets["y"])
                & np.isfinite(dets["flux"]))
        n_bad = int((~keep).sum())
        if n_bad:
            dets = {k: np.asarray(v)[keep] for k, v in dets.items()}
            sat = sat[keep]
            stats = dict(stats, n_detected=int(keep.sum()))
        out.update(ok=True, recipe=recipe, plate_scale=scale,
                   scale_basis=scale_basis, aper_px=aper,
                   veto_used=veto, stats=stats, dets=dets,
                   saturated=sat.astype(int),
                   n_saturated=int(sat.sum()), n_nonfinite=n_bad)
    except Exception as e:
        out.update(ok=False, error=f"{type(e).__name__}: {e}"[:200])
    return out


def cmd_extract(args) -> None:
    """Measure up to --limit pending frames with a pool of worker processes."""
    con = connect(args.db)
    todo = con.execute("""SELECT frame_id, pixel_path, provenance,
                                 master_dark, master_flat, veto_adu
                          FROM cv_frames WHERE status='pending'
                          ORDER BY frame_id LIMIT ?""",
                       (args.limit,)).fetchall()
    if not todo:
        print("extract: nothing pending", flush=True)
        con.close()
        return
    jobs = [{"frame_id": r[0],
             "pixel_path": str(args.archive / r[1]),
             "provenance": r[2],
             "master_dark": r[3], "master_flat": r[4],
             "veto_adu": r[5]} for r in todo]
    nw = clamp_workers(args.workers)
    t0 = time.time()
    n_ok = n_fail = 0
    pending_rows: list = []

    def flush(rows):
        """Write one batch inside ONE short transaction (polite to the S1 run)."""
        if not rows:
            return
        con.execute("BEGIN")
        for res in rows:
            fid = res["frame_id"]
            if not res["ok"]:
                con.execute("""UPDATE cv_frames SET status='failed_extract',
                               exclude_reason=? WHERE frame_id=?""",
                            (res["error"], fid))
                continue
            s, d = res["stats"], res["dets"]
            con.execute("""UPDATE cv_frames SET status='extracted',
                calib_recipe=?, plate_scale=?, scale_basis=?, aper_px=?,
                veto_adu=?, n_detected=?, n_saturated=?, n_nonfinite=?,
                bkg_adu=?, bkg_rms=?, fwhm_px=?
                WHERE frame_id=?""",
                (res["recipe"], fnum(res["plate_scale"]), res["scale_basis"],
                 fnum(res["aper_px"]),
                 fnum(res["veto_used"]), s["n_detected"], res["n_saturated"],
                 res["n_nonfinite"],
                 fnum(s["bkg_adu"]), fnum(s["bkg_rms"]), fnum(s["fwhm_px"]),
                 fid))
            con.execute("DELETE FROM cv_detections WHERE frame_id=?", (fid,))
            con.executemany("""INSERT INTO cv_detections
                (frame_id, det_id, x, y, flux, fluxerr, fwhm, peak,
                 sep_flag, saturated, star_id)
                VALUES (?,?,?,?,?,?,?,?,?,?,NULL)""",
                [(fid, i, float(d["x"][i]), float(d["y"][i]),
                  float(d["flux"][i]), float(d["fluxerr"][i]),
                  float(d["fwhm"][i]), float(d["peak"][i]),
                  int(d["flag"][i]), int(res["saturated"][i]))
                 for i in range(len(d["x"]))])
        con.commit()

    with cf.ProcessPoolExecutor(max_workers=nw) as pool:
        for res in pool.map(_extract_one, jobs, chunksize=4):
            pending_rows.append(res)
            n_ok += 1 if res["ok"] else 0
            n_fail += 0 if res["ok"] else 1
            if len(pending_rows) >= WRITE_BATCH:
                flush(pending_rows)
                pending_rows = []
                done = n_ok + n_fail
                print(f"  extract {done}/{len(jobs)} "
                      f"({done / max(time.time() - t0, 1e-9):.1f} frames/s)",
                      flush=True)
    flush(pending_rows)
    left = con.execute("SELECT count(*) FROM cv_frames "
                       "WHERE status='pending'").fetchone()[0]
    print(f"extract: {n_ok} ok, {n_fail} failed, {time.time() - t0:.0f}s "
          f"({nw} workers); pending {left}", flush=True)
    con.close()


# ---------------------------------------------------------------------------
# Stage: match (reference selection + registration, parallel, chunked)
# ---------------------------------------------------------------------------
def _has_wcs(wcs_root: Path, raw_path: str) -> bool:
    return reg.sidecar_path(wcs_root, raw_path).exists()


def _ensure_reference(con, args, target_key: str, era_id: int):
    """Choose (once), quality-control, and return the (target, era) reference.

    Two rules stack on the prototype's:

    1.  Candidates are walked in ``rank_references`` order and each must
        pass the DOUBLE-IMAGE test (:func:`macro_phot.photometry.
        paired_fraction`) — a guiding jump doubles every star into an
        equal-brightness pair, INFLATES the detection count, and would
        otherwise be rewarded by 'most detections among the sharpest'.
    2.  A PLATE-SOLVED candidate is preferred over an unsolved one of
        equal standing, because a solved reference unlocks the sky-chained
        registration route for every solved frame in the era and lets the
        target be identified from its catalogue coordinates without a
        network query.  Solved candidates are therefore tried first, and
        the whole ranking is retried only if none of them survives QC.
    """
    con.execute("""CREATE TABLE IF NOT EXISTS cv_ref (
        target_key TEXT, era_id INTEGER, ref_frame_id INTEGER,
        ref_raw_path TEXT, ref_has_wcs INTEGER,
        n_stars INTEGER, fwhm_px REAL, tol_px REAL, plate_scale REAL,
        doubled_frac REAL, n_cand_rejected INTEGER,
        PRIMARY KEY (target_key, era_id))""")
    con.execute("""CREATE TABLE IF NOT EXISTS cv_ref_stars (
        target_key TEXT, era_id INTEGER, star_id INTEGER,
        x REAL, y REAL, flux REAL,
        ra_deg REAL, dec_deg REAL,
        gaia_source_id INTEGER, gaia_gmag REAL,
        PRIMARY KEY (target_key, era_id, star_id))""")
    row = con.execute("""SELECT ref_frame_id, tol_px FROM cv_ref
                         WHERE target_key=? AND era_id=?""",
                      (target_key, era_id)).fetchone()
    if not row:
        stats = con.execute("""SELECT frame_id, n_detected, fwhm_px, raw_path
                               FROM cv_frames
                               WHERE target_key=? AND era_id=?
                                 AND status IN ('extracted','matched')""",
                            (target_key, era_id)).fetchall()
        if not stats:
            return None, None, None
        solved = {s[0]: (1 if _has_wcs(args.wcs_root, s[3]) else 0)
                  for s in stats}
        fwhm_of = {s[0]: s[2] for s in stats}
        ranking = ph.rank_references([(s[0], s[1], s[2], solved[s[0]])
                                      for s in stats])
        # Solved candidates go first, but only if they are deep enough to
        # be worth it — the reference's star list is the ceiling on what
        # any frame in the era can match.
        ordered = sr.order_reference_candidates(
            ranking, {s[0]: s[1] for s in stats}, solved)
        ref_id = None
        doubled = None
        n_rejected = 0
        for cand in ordered:
            cdets = con.execute("""SELECT x, y, flux FROM cv_detections
                                   WHERE frame_id=? ORDER BY det_id""",
                                (cand,)).fetchall()
            if len(cdets) < ph.MIN_STARS_FOR_ALIGN:
                n_rejected += 1
                continue
            cxy = np.array([[d[0], d[1]] for d in cdets], dtype=float)
            cfl = np.array([d[2] for d in cdets], dtype=float)
            frac = ph.paired_fraction(
                cxy, cfl, ph.REF_PAIR_RADIUS_FWHM * float(fwhm_of[cand]))
            if frac <= ph.REF_DOUBLED_MAX_FRAC:
                ref_id, doubled = cand, frac
                break
            n_rejected += 1
            print(f"  reference QC: {target_key}/era{era_id} frame {cand} "
                  f"REJECTED — paired fraction {frac:.2f} > "
                  f"{ph.REF_DOUBLED_MAX_FRAC:g} (double-imaged)", flush=True)
        if ref_id is None:
            return None, None, None
        fwhm = fwhm_of[ref_id]
        tol = ph.match_tolerance_px(fwhm)
        meta = con.execute("""SELECT raw_path, plate_scale FROM cv_frames
                              WHERE frame_id=?""", (ref_id,)).fetchone()
        dets = con.execute("""SELECT det_id, x, y, flux FROM cv_detections
                              WHERE frame_id=? ORDER BY det_id""",
                           (ref_id,)).fetchall()
        con.execute("BEGIN")
        con.executemany("""INSERT OR REPLACE INTO cv_ref_stars
            (target_key, era_id, star_id, x, y, flux)
            VALUES (?,?,?,?,?,?)""",
            [(target_key, era_id, d, x, y, fl) for d, x, y, fl in dets])
        con.execute("""INSERT OR REPLACE INTO cv_ref
                       VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                    (target_key, era_id, ref_id, meta[0],
                     solved.get(ref_id, 0), len(dets), fnum(fwhm), fnum(tol),
                     fnum(meta[1]), fnum(doubled), n_rejected))
        # The reference matches itself by identity, so its own frame joins
        # the ensemble on the same footing as every other frame.
        con.execute("""UPDATE cv_detections SET star_id = det_id
                       WHERE frame_id=?""", (ref_id,))
        con.execute("""UPDATE cv_frames SET status='matched',
                       reg_method='reference', ali_nmatch=?, ali_rms_px=0.0,
                       ali_scale=1.0, ali_rot_deg=0.0, match_rate=1.0,
                       has_wcs=? WHERE frame_id=?""",
                    (len(dets), solved.get(ref_id, 0), ref_id))
        con.commit()
        print(f"  reference {target_key}/era{era_id}: frame {ref_id}, "
              f"{len(dets)} stars, FWHM {fwhm:.2f} px, tol {tol:.2f} px, "
              f"wcs={bool(solved.get(ref_id))}, doubled_frac={doubled:.3f}",
              flush=True)
    else:
        ref_id, tol = row
    stars = con.execute("""SELECT star_id, x, y FROM cv_ref_stars
                           WHERE target_key=? AND era_id=? ORDER BY star_id""",
                        (target_key, era_id)).fetchall()
    return ref_id, tol, stars


def _match_one(job: dict) -> dict:
    """Register ONE frame against its series reference. Worker process.

    Tries the sky chain first when both ends are plate-solved, because it
    is cheap and cannot starve; falls back to the astroalign triangle
    ladder when the chain is missing or its match fails
    :func:`macro_phot.series.wcs_match_ok`.  The method actually used is
    returned and stored per frame — a run in which the chain silently
    stopped working is then visible as a column, not as a mystery.
    """
    out = {"frame_id": job["frame_id"]}
    xy = np.asarray(job["xy"], dtype=float)
    ref_xy = np.asarray(job["ref_xy"], dtype=float)
    tol = float(job["tol"])
    try:
        idx = None
        method = None
        scale = rot = None
        if job["method"] == "wcs":
            fw = reg.load_wcs(Path(job["frame_wcs"]))
            rw = reg.load_wcs(Path(job["ref_wcs"]))
            if fw is not None and rw is not None:
                moved = reg.chain_to_reference(fw, rw, xy)
                good = np.isfinite(moved).all(axis=1)
                cand = np.full(len(xy), -1, dtype=int)
                if good.any():
                    cand[good] = ph.match_one_to_one(ref_xy, moved[good], tol)
                n = int((cand >= 0).sum())
                if sr.wcs_match_ok(n, len(xy), len(ref_xy)):
                    idx, method = cand, "wcs"
                    out["moved"] = moved
        if idx is None:
            if len(xy) < ph.MIN_STARS_FOR_ALIGN:
                raise ValueError("too few detections to align")
            bright = np.asarray(job["bright"], dtype=float)
            ref_bright = np.asarray(job["ref_bright"], dtype=float)
            # seed = frame_id pins astroalign's otherwise-unseedable RANSAC,
            # so star identities reproduce bit-for-bit on a re-run.
            try:
                tf = ext.find_series_transform(bright, ref_bright,
                                               seed=job["frame_id"])
            except Exception:
                # DEPTH-MATCHED RETRY.  The ladder already builds both
                # bright pools the same way; what it cannot fix is the two
                # catalogs reaching DIFFERENT DEPTHS.  A reference is
                # chosen as the deepest frame of its era, and ST LMi's 2024
                # season runs from 16 to 815 detections per frame — so a
                # cloudy 50-star frame's "brightest 100" is its whole
                # catalog while the reference's "brightest 100" reaches two
                # magnitudes fainter, the two control sets barely overlap,
                # and astroalign exhausts its triangles on stars the frame
                # never saw.  Truncating BOTH catalogs to the shallower
                # one's length puts the two control sets back on the same
                # brightness range.  Tried second, so frames that already
                # matched keep exactly the transform they had.
                n = min(len(bright), len(ref_bright))
                if n < ph.MIN_STARS_FOR_ALIGN:
                    raise
                tf = ext.find_series_transform(bright[:n], ref_bright[:n],
                                               seed=job["frame_id"])
                method_suffix = "_depthmatched"
                job = dict(job, _suffix=method_suffix)
            moved = tf(xy)
            idx = ph.match_one_to_one(ref_xy, moved, tol)
            method = ("astroalign" if job["method"] == "astroalign"
                      else "astroalign_fallback") + job.get("_suffix", "")
            scale, rot = float(tf.scale), float(np.degrees(tf.rotation))
            out["moved"] = moved
        matched = idx >= 0
        moved = out.pop("moved")
        if matched.any():
            d = np.hypot(moved[matched, 0] - ref_xy[idx[matched], 0],
                         moved[matched, 1] - ref_xy[idx[matched], 1])
            rms = float(np.sqrt(np.mean(d ** 2)))
        else:
            rms = None
        out.update(ok=True, method=method, idx=idx.tolist(),
                   n_match=int(matched.sum()), rms=rms,
                   scale=scale, rot=rot,
                   rate=sr.match_rate(int(matched.sum()), len(xy),
                                      len(ref_xy)))
    except Exception as e:
        out.update(ok=False, error=f"{type(e).__name__}: {e}"[:200])
    return out


def cmd_match(args) -> None:
    """Register up to --limit extracted frames onto their series reference."""
    con = connect(args.db)
    eras = con.execute("""SELECT DISTINCT target_key, era_id FROM cv_frames
                          WHERE status IN ('extracted','matched')
                          ORDER BY 1, 2""").fetchall()
    nw = clamp_workers(args.workers)
    budget = args.limit
    t0 = time.time()
    n_ok = n_fail = 0
    # ONE pool for the whole invocation.  macOS starts worker processes with
    # 'spawn', which re-imports this module and numpy/sep/astropy in every
    # child — about a second each.  Creating a pool per 100-frame block paid
    # that toll dozens of times over; creating it once amortizes it to
    # nothing.  The pool is closed by the surrounding 'with' on every exit
    # path, including an exception.
    pool = cf.ProcessPoolExecutor(max_workers=nw)
    for tk, era in eras:
        if budget <= 0:
            break
        ref_id, tol, stars = _ensure_reference(con, args, tk, era)
        if ref_id is None:
            print(f"  {tk}/era{era}: no usable reference frame", flush=True)
            continue
        ref_ids = np.array([s[0] for s in stars])
        ref_xy = np.array([[s[1], s[2]] for s in stars], dtype=float)
        max_pool = max(p for p, _ in ext.ALIGN_ATTEMPTS)
        # Bright pools are built the SAME way on both sides — unclipped,
        # flux-descending — because astroalign's control points are the
        # first N entries and an asymmetric pool poisons them from the top.
        ref_bright = np.array([[r[0], r[1]] for r in con.execute(
            """SELECT x, y FROM cv_detections
               WHERE frame_id=? AND saturated=0 ORDER BY flux DESC LIMIT ?""",
            (ref_id, max_pool))], dtype=float)
        ref_raw = con.execute("SELECT raw_path FROM cv_frames WHERE frame_id=?",
                              (ref_id,)).fetchone()[0]
        ref_wcs = reg.sidecar_path(args.wcs_root, ref_raw)
        ref_solved = ref_wcs.exists()

        # Frames S0 already proved were pointed elsewhere share no stars
        # with the reference; astroalign burns ~15 s per frame discovering
        # that.  The skip is a recorded status, never a silent drop.
        con.execute("""UPDATE cv_frames SET status='skipped_pointing',
            exclude_reason=?
            WHERE target_key=? AND era_id=? AND status='extracted'
              AND qc_flags LIKE '%pointing_gt1deg%'""",
            ("S0 flagged pointing_gt1deg: the frame does not contain the "
             "reference field", tk, era))
        con.commit()

        todo = con.execute("""SELECT frame_id, raw_path FROM cv_frames
                              WHERE target_key=? AND era_id=?
                                AND status='extracted'
                              ORDER BY frame_id LIMIT ?""",
                           (tk, era, budget)).fetchall()
        if not todo:
            continue
        print(f"  {tk}/era{era}: matching {len(todo)} frames against "
              f"reference {ref_id} ({len(stars)} stars, "
              f"ref_wcs={ref_solved})", flush=True)
        # Build jobs in blocks so detection arrays for at most WRITE_BATCH
        # frames are ever resident at once.
        for i0 in range(0, len(todo), WRITE_BATCH):
            block = todo[i0:i0 + WRITE_BATCH]
            jobs = []
            det_ids: dict = {}
            for fid, raw_path in block:
                # 'flux IS NOT NULL' is defensive, not decorative: SQLite
                # stores a NaN float as NULL, so any measurement that could
                # not be made reads back as None and would sort as None
                # inside the bright pool below.
                dets = con.execute("""SELECT det_id, x, y, flux, saturated
                                      FROM cv_detections WHERE frame_id=?
                                        AND flux IS NOT NULL
                                      ORDER BY det_id""", (fid,)).fetchall()
                det_ids[fid] = [d[0] for d in dets]
                xy = [[d[1], d[2]] for d in dets]
                bright = [[d[1], d[2]] for d in
                          sorted((d for d in dets if not d[4]),
                                 key=lambda d: -d[3])[:max_pool]]
                fw = reg.sidecar_path(args.wcs_root, raw_path)
                has = fw.exists()
                jobs.append({
                    "frame_id": fid, "xy": xy, "bright": bright,
                    "ref_xy": ref_xy.tolist(),
                    "ref_bright": ref_bright.tolist(), "tol": tol,
                    "method": sr.registration_method(has, ref_solved),
                    "frame_wcs": str(fw), "ref_wcs": str(ref_wcs),
                    "has_wcs": 1 if has else 0})
            results = list(pool.map(_match_one, jobs, chunksize=2))
            con.execute("BEGIN")
            for res, job in zip(results, jobs):
                fid = res["frame_id"]
                if not res["ok"]:
                    con.execute("""UPDATE cv_frames SET status='failed_match',
                        exclude_reason=?, has_wcs=? WHERE frame_id=?""",
                        (res["error"], job["has_wcs"], fid))
                    n_fail += 1
                    continue
                ids = det_ids[fid]
                idx = res["idx"]
                con.executemany(
                    """UPDATE cv_detections SET star_id=?
                       WHERE frame_id=? AND det_id=?""",
                    [(int(ref_ids[idx[k]]), fid, int(ids[k]))
                     for k in range(len(ids)) if idx[k] >= 0])
                con.execute("""UPDATE cv_frames SET status='matched',
                    reg_method=?, ali_nmatch=?, ali_rms_px=?, ali_scale=?,
                    ali_rot_deg=?, match_rate=?, has_wcs=?
                    WHERE frame_id=?""",
                    (res["method"], res["n_match"], fnum(res["rms"]),
                     fnum(res["scale"]), fnum(res["rot"]), fnum(res["rate"]),
                     job["has_wcs"], fid))
                n_ok += 1
            con.commit()
            budget -= len(block)
            print(f"    {n_ok + n_fail} matched so far "
                  f"({(n_ok + n_fail) / max(time.time() - t0, 1e-9):.1f}/s)",
                  flush=True)
            if budget <= 0:
                break
    pool.shutdown()
    left = con.execute("SELECT count(*) FROM cv_frames "
                       "WHERE status='extracted'").fetchone()[0]
    print(f"match: {n_ok} ok, {n_fail} failed, {time.time() - t0:.0f}s "
          f"({nw} workers); extracted-but-unmatched {left}", flush=True)
    con.close()


# ---------------------------------------------------------------------------
# Stage: field (identify the target in each (target, era))
# ---------------------------------------------------------------------------

def _tie_upsert(con, row: tuple) -> None:
    """Record ONE (target, era) field tie immediately, not at the end.

    Written the moment it is known so a later network failure cannot undo
    an identification that already cost a slow archive query.
    """
    con.execute("INSERT OR REPLACE INTO cv_field_tie VALUES "
                "(?,?,?,?,?,?,?,?,?,?,?,?,?)", row)
    con.commit()


def cmd_field(args) -> None:
    """Find THE TARGET among each (target, era) reference catalog.

    Two routes, and the cheap one is tried first:

    ``wcs``
        The reference frame carries an S1 plate solution.  The target's
        catalogue coordinates convert straight to reference pixels, every
        reference star gets a sky position from the same solution, and the
        nearest star within the tolerance is the target.  No network, no
        triangle fit, and the identification inherits S1's own astrometric
        residual.
    ``gaia``
        The reference is unsolved (whole eras of these polars are).  The
        prototype's route runs: a Gaia DR3 cone about the target, a
        parity-tolerant similarity fit of reference pixels to the Gaia
        tangent plane, then the same nearest-star rule.

    Identifying the target matters for two independent reasons, and both
    are load-bearing: the target's light curve is the deliverable, and the
    target must be EXCLUDED from the comparison-star pool — a polar's
    orbital modulation must never help set the zero point it is measured
    against.
    """
    con = connect(args.db)
    # RESUMABLE by design.  The Gaia archive is a shared public service
    # whose transient 500s cost minutes each (measured during this build),
    # and rebuilding the whole table on every run would throw away every
    # tie already paid for whenever the last one failed.  So the table
    # persists and each (target, era) is upserted; --force recomputes.
    con.execute("""CREATE TABLE IF NOT EXISTS cv_field_tie (
        target_key TEXT, era_id INTEGER, method TEXT, parity TEXT,
        scale_fit REAL, rot_deg REAL, target_star_id INTEGER,
        n_gaia INTEGER, n_gaia_matched INTEGER,
        target_ra REAL, target_dec REAL, coord_source TEXT, status TEXT,
        PRIMARY KEY (target_key, era_id))""")
    done = set() if args.force else {
        (r[0], r[1]) for r in con.execute(
            "SELECT target_key, era_id FROM cv_field_tie WHERE status='ok'")}
    refs = con.execute("""SELECT r.target_key, r.era_id, r.ref_frame_id,
                                 r.ref_raw_path, r.ref_has_wcs, r.plate_scale
                          FROM cv_ref r ORDER BY 1, 2""").fetchall()
    for tk, era, ref_id, ref_raw, has_wcs, scale in refs:
        if (tk, era) in done:
            print(f"  field {tk}/era{era}: already tied — skipping "
                  f"(use --force to redo)", flush=True)
            continue
        ra0, dec0, coord_src = gg.resolve_target(tk)
        stars = con.execute("""SELECT star_id, x, y FROM cv_ref_stars
                               WHERE target_key=? AND era_id=?
                               ORDER BY star_id""", (tk, era)).fetchall()
        if not stars:
            continue
        sids = np.array([s[0] for s in stars])
        ref_xy = np.array([[s[1], s[2]] for s in stars], dtype=float)
        method = None
        radec = None
        n_gaia = n_matched = 0
        gaia_idx = None
        gaia_tab = None
        parity = scale_fit = rot = None

        if has_wcs:
            w = reg.load_wcs(reg.sidecar_path(args.wcs_root, ref_raw))
            if w is not None:
                try:
                    sky = w.all_pix2world(ref_xy[:, 0], ref_xy[:, 1], 0)
                    radec = np.column_stack([np.asarray(sky[0], dtype=float),
                                             np.asarray(sky[1], dtype=float)])
                    method = "wcs"
                except Exception:
                    radec = None
            if radec is not None and args.gaia_identities:
                # A solved reference already knows where its stars are, so
                # the Gaia tie needs no triangle fit at all — the two star
                # lists are cross-matched directly on the sky.  This is not
                # required for the photometry (which is instrumental
                # throughout), but it hands the next analyst a catalogue
                # identity and a G magnitude for every comparison star,
                # which is what makes a comparison set checkable by someone
                # who was not here.  A failed cone query is not fatal: the
                # field tie stands on the WCS, and the Gaia columns simply
                # stay empty.  It is OPT-IN (--gaia-identities) because the
                # archive is a shared public service whose latency this
                # campaign measured at 6 to 68 seconds per identical query,
                # with intermittent 500s on top: a convenience column must
                # not be allowed to hold the whole campaign hostage.  An
                # UNSOLVED reference has no such choice — Gaia is the only
                # route to a target identification there — so that branch
                # always queries.
                try:
                    sc = scale or 0.45
                    n1, n2 = con.execute(
                        """SELECT max(x), max(y) FROM cv_detections
                           WHERE frame_id=?""", (ref_id,)).fetchone()
                    radius = (1.2 * sc * float(np.hypot(n1 or 0, n2 or 0))
                              / 2.0 / 3600.0)
                    gaia_tab = gg.cone_query(ra0, dec0, radius)
                    n_gaia = len(gaia_tab["ra"])
                    # Project both catalogues onto the same tangent plane
                    # (arcsec) and reuse the pure one-to-one matcher, so
                    # the WCS route and the astroalign route share exactly
                    # one matching rule.
                    gx, gy = gg.tangent_project(gaia_tab["ra"],
                                                gaia_tab["dec"], ra0, dec0)
                    rx, ry = gg.tangent_project(radec[:, 0], radec[:, 1],
                                                ra0, dec0)
                    gaia_idx = ph.match_one_to_one(
                        np.column_stack([gx, gy]),
                        np.column_stack([rx, ry]),
                        gg.GAIA_MATCH_TOL_ARCSEC)
                    n_matched = int((gaia_idx >= 0).sum())
                except Exception as e:
                    print(f"  gaia (identities only) failed for "
                          f"{tk}/era{era}: {type(e).__name__}", flush=True)
                    gaia_tab, gaia_idx, n_gaia, n_matched = None, None, 0, 0
        if method is None:
            # ---- Gaia route (unsolved reference)
            max_pool = max(p for p, _ in ext.ALIGN_ATTEMPTS)
            ref_bright = np.array([[r[0], r[1]] for r in con.execute(
                """SELECT x, y FROM cv_detections WHERE frame_id=?
                   AND saturated=0 ORDER BY flux DESC LIMIT ?""",
                (ref_id, max_pool))], dtype=float)
            n1, n2 = con.execute("""SELECT max(x), max(y) FROM cv_detections
                                    WHERE frame_id=?""", (ref_id,)).fetchone()
            sc = scale or 0.45
            radius = 1.2 * sc * float(np.hypot(n1 or 0, n2 or 0)) / 2 / 3600.0
            fit_radius = sc * float(min(n1 or 0, n2 or 0)) / 2.0
            try:
                gaia_tab = gg.cone_query(ra0, dec0, radius)
                n_gaia = len(gaia_tab["ra"])
            except Exception as e:
                print(f"  gaia: cone failed for {tk}/era{era}: {e}",
                      flush=True)
                _tie_upsert(con, (tk, era, "none", None, None, None, None,
                                  0, 0, ra0, dec0, coord_src,
                                  f"query_failed: {type(e).__name__}"))
                continue
            fit = gg.identify_reference(ref_xy, gaia_tab, ra0, dec0,
                                        ref_bright_xy=ref_bright,
                                        fit_radius_arcsec=fit_radius,
                                        seed=ref_id)
            if fit is None:
                _tie_upsert(con, (tk, era, "none", None, None, None, None,
                                  n_gaia, 0, ra0, dec0, coord_src,
                                  "fit_failed"))
                continue
            method = "gaia"
            radec = fit["ref_radec"]
            gaia_idx = fit["gaia_idx"]
            n_matched = fit["n_matched"]
            parity = fit["parity"]
            scale_fit = fit["scale_arcsec_per_px"]
            rot = fit["rot_deg"]

        # ---- write sky positions (and Gaia identities where we have them)
        con.execute("BEGIN")
        for i, sid in enumerate(sids):
            g_id = g_mag = None
            if gaia_idx is not None and gaia_idx[i] >= 0:
                g_id = int(gaia_tab["source_id"][gaia_idx[i]])
                g_mag = float(gaia_tab["gmag"][gaia_idx[i]])
            con.execute("""UPDATE cv_ref_stars SET ra_deg=?, dec_deg=?,
                gaia_source_id=?, gaia_gmag=?
                WHERE target_key=? AND era_id=? AND star_id=?""",
                (fnum(radec[i, 0]), fnum(radec[i, 1]), g_id, g_mag,
                 tk, era, int(sid)))
        con.commit()

        # ---- the target: the reference star nearest the catalogue position
        sep_deg = np.hypot((radec[:, 0] - ra0) * np.cos(np.radians(dec0)),
                           radec[:, 1] - dec0)
        j = int(np.argmin(sep_deg))
        sep_arcsec = float(sep_deg[j] * 3600.0)
        target_sid = (int(sids[j]) if sep_arcsec <= gg.TARGET_ID_TOL_ARCSEC
                      else None)
        status = "ok" if target_sid is not None else (
            f"target not in reference catalog (nearest star "
            f"{sep_arcsec:.1f}\" away, tolerance "
            f"{gg.TARGET_ID_TOL_ARCSEC:g}\")")
        _tie_upsert(con, (tk, era, method, parity, fnum(scale_fit),
                          fnum(rot), target_sid, n_gaia, n_matched, ra0,
                          dec0, coord_src, status))
        print(f"  field {tk}/era{era}: via {method}; target star "
              f"{target_sid} at {sep_arcsec:.2f}\"; "
              f"gaia matched {n_matched}/{len(sids)}", flush=True)

    meta_write(con)
    con.close()


# ---------------------------------------------------------------------------
# Stage: ensemble
# ---------------------------------------------------------------------------
def _series_matrix(con, skey: str):
    """Assemble the (stars x frames) magnitude/error matrices of ONE series.

    Only frames of THIS series key contribute — one target, one era, one
    filter — so era mixing is structurally impossible rather than merely
    discouraged.  SATURATED detections are dropped at load: a measurement
    at the digitization ceiling is not a faint version of a bright star,
    it is a censored number, and fitting it would drag the zero point.
    They remain in ``cv_detections`` (flagged) and are re-emitted into the
    light curve so the next analyst sees the censoring instead of a gap.
    """
    frames = con.execute("""SELECT frame_id, exptime, bjd_tdb, night
                            FROM cv_frames
                            WHERE series_key=? AND status='matched'
                            ORDER BY bjd_tdb""", (skey,)).fetchall()
    if not frames:
        return None
    tk, era, _filt = sr.parse_series_key(skey)
    stars = con.execute("""SELECT star_id FROM cv_ref_stars
                           WHERE target_key=? AND era_id=?
                           ORDER BY star_id""", (tk, era)).fetchall()
    if not stars:
        return None
    sid_row = {s[0]: i for i, s in enumerate(stars)}
    fid_col = {f[0]: j for j, f in enumerate(frames)}
    mag = np.full((len(stars), len(frames)), np.nan)
    sig = np.full((len(stars), len(frames)), np.nan)
    expt = {f[0]: f[1] for f in frames}
    # Chunked IN-lists: a single 2,000-frame IN clause exceeds SQLite's
    # default variable limit, and string-building one giant literal makes
    # the query planner miserable.
    fids = [f[0] for f in frames]
    for i0 in range(0, len(fids), 400):
        block = fids[i0:i0 + 400]
        qs = ",".join("?" * len(block))
        rows = con.execute(f"""SELECT frame_id, star_id, flux, fluxerr
                               FROM cv_detections
                               WHERE frame_id IN ({qs})
                                 AND star_id IS NOT NULL
                                 AND saturated = 0""", block).fetchall()
        for fid, sid, flux, fluxerr in rows:
            i, j = sid_row[sid], fid_col[fid]
            mag[i, j] = ph.instrumental_mag(np.array([flux]), expt[fid])[0]
            sig[i, j] = ph.mag_error(np.array([flux]),
                                     np.array([fluxerr]))[0]
    return (np.array([s[0] for s in stars]), np.array(fids),
            np.array([f[2] for f in frames], dtype=float), mag, sig)


def cmd_ensemble(args) -> None:
    """Solve the Honeycutt ensemble of every (target, era, filter) series."""
    con = connect(args.db)
    tie = {(r[0], r[1]): r[2] for r in con.execute(
        "SELECT target_key, era_id, target_star_id FROM cv_field_tie")}
    keys = [r[0] for r in con.execute(
        """SELECT DISTINCT series_key FROM cv_frames WHERE status='matched'
           ORDER BY 1""")]
    star_rows, series_rows, zp_rows, lc_rows = [], [], [], []
    for skey in keys:
        tk, era, filt = sr.parse_series_key(skey)
        packed = _series_matrix(con, skey)
        if packed is None:
            continue
        sids, fids, bjd, mag, sig = packed
        admitted, why = sr.series_admission(len(fids))
        rate = con.execute("""SELECT avg(match_rate) FROM cv_frames
                              WHERE series_key=? AND status='matched'""",
                           (skey,)).fetchone()[0]
        n_staged, n_sat = con.execute(
            """SELECT count(*), sum(coalesce(n_saturated,0)) FROM cv_frames
               WHERE series_key=?""", (skey,)).fetchone()
        prov = con.execute("SELECT provenance FROM cv_selection "
                           "WHERE series_key=?", (skey,)).fetchone()
        prov = prov[0] if prov else None
        if not admitted:
            series_rows.append((skey, tk, era, filt, prov, n_staged,
                                len(fids), fnum(rate), n_sat or 0, 0, 0,
                                0, 0, 0, None, 0, 0, None, None, None, None,
                                None, None, "not_solved", why))
            print(f"  {skey:24s} NOT SOLVED — {why}", flush=True)
            continue

        target_sid = tie.get((tk, era))
        target_row = None
        if target_sid is not None and (sids == target_sid).any():
            target_row = int(np.flatnonzero(sids == target_sid)[0])
        sel = ens.select_comps(mag, sig, target_row=target_row)
        sol = sel.solution
        mean, rms, nobs, chi2nu = ens.star_stats(mag, sig, sol.zp)
        roles = sel.role
        n_comp = int((roles == "comp").sum())
        n_check = int((roles == "check").sum())
        n_drop = int((roles == "dropped_unstable").sum())

        gaia_cols = {r[0]: (r[1], r[2]) for r in con.execute(
            """SELECT star_id, gaia_source_id, gaia_gmag FROM cv_ref_stars
               WHERE target_key=? AND era_id=?""", (tk, era))}
        for i, sid in enumerate(sids):
            if nobs[i] == 0:
                continue
            g_id, g_mag = gaia_cols.get(int(sid), (None, None))
            star_rows.append((skey, tk, era, filt, int(sid), roles[i],
                              fnum(mean[i]), fnum(rms[i]), int(nobs[i]),
                              fnum(chi2nu[i]), g_id, g_mag))

        # ---- zero points onto the frames and into the timeline table.
        # A frame whose comparison stars were all absent or saturated gets
        # NaN from the solver and is stored as NULL, never 0.0: a
        # fabricated zero point looks exactly like a good one downstream.
        con.execute("BEGIN")
        con.executemany("""UPDATE cv_frames SET zp=?, zp_err=?, n_star_used=?
                           WHERE frame_id=?""",
                        [(fnum(sol.zp[j]), fnum(sol.zp_err[j]),
                          int(sol.n_star_used[j]), int(fids[j]))
                         for j in range(len(fids))])
        con.commit()
        zp_rows.extend([(skey, int(fids[j]), float(bjd[j]),
                         fnum(sol.zp[j]), fnum(sol.zp_err[j]),
                         int(sol.n_star_used[j]))
                        for j in range(len(fids))])

        # ---- the deliverable light curves: target, comps and checks.
        # Saturated points are re-attached here (they were withheld from
        # the FIT) so the series is complete and censoring is visible.
        of_interest = {int(sids[i]): roles[i] for i in range(len(sids))
                       if roles[i] in ("target", "comp", "check")}
        if of_interest:
            col = {int(f): j for j, f in enumerate(fids)}
            # Exposure times are fetched ONCE per series, not once per
            # measurement: the naive per-row lookup turned a two-second
            # write into a quarter of a million single-row queries.
            expt = {r[0]: r[1] for r in con.execute(
                """SELECT frame_id, exptime FROM cv_frames
                   WHERE series_key=? AND status='matched'""", (skey,))}
            for i0 in range(0, len(fids), 400):
                block = [int(f) for f in fids[i0:i0 + 400]]
                qs = ",".join("?" * len(block))
                for fid, sid, flux, ferr, sat in con.execute(
                        f"""SELECT frame_id, star_id, flux, fluxerr, saturated
                            FROM cv_detections WHERE frame_id IN ({qs})
                              AND star_id IS NOT NULL""", block):
                    role = of_interest.get(int(sid))
                    if role is None:
                        continue
                    j = col[fid]
                    m = ph.instrumental_mag(np.array([flux]), expt[fid])[0]
                    me = ph.mag_error(np.array([flux]), np.array([ferr]))[0]
                    zp = sol.zp[j]
                    lc_rows.append((skey, int(sid), role, fid,
                                    float(bjd[j]), fnum(m), fnum(me),
                                    fnum(zp),
                                    fnum(m - zp) if np.isfinite(zp) else None,
                                    int(sat)))

        # ---- how often the TARGET itself hit the digitization ceiling.
        # The campaign's contract is that a saturated target measurement is
        # flagged, never silently fitted; this is the number that lets the
        # next analyst see at a glance whether that mattered here.
        n_tgt = n_tgt_sat = 0
        for row in lc_rows:
            if row[0] == skey and row[2] == "target":
                n_tgt += 1
                n_tgt_sat += int(row[9])

        # ---- per-series diagnostics
        chk = np.flatnonzero(roles == "check")
        chk_rms = fnum(np.nanmedian(rms[chk])) if chk.size else None
        infl = err.inflation_factor(chi2nu[chk]) if chk.size else float("nan")
        zp_ok = sol.zp[np.isfinite(sol.zp)]
        series_rows.append((
            skey, tk, era, filt, prov, n_staged, len(fids), fnum(rate),
            n_sat or 0, n_tgt, n_tgt_sat,
            n_comp, n_check, n_drop, sel.n_passes, sol.n_iter,
            int(sol.converged), fnum(sel.comp_rms_median), chk_rms,
            fnum(infl), fnum(np.std(zp_ok)) if zp_ok.size else None,
            target_sid, sr.check_star_verdict(n_check), "solved", why))
        print(f"  {skey:24s} n={len(fids):5d} rate={rate or float('nan'):.2f} "
              f"comps={n_comp:3d} checks={n_check} dropped={n_drop:3d} "
              f"passes={sel.n_passes} "
              f"iters={sol.n_iter:3d} conv={int(sol.converged)} "
              f"checkRMS={chk_rms if chk_rms is not None else float('nan'):.4f} "
              f"infl={infl:.2f}", flush=True)

    swap_in(con, "cv_stars",
            """CREATE TABLE {t} (series_key TEXT, target_key TEXT,
               era_id INTEGER, filter TEXT, star_id INTEGER, role TEXT,
               mean_mag REAL, rms REAL, nobs INTEGER, chi2nu REAL,
               gaia_source_id INTEGER, gaia_gmag REAL)""", star_rows)
    swap_in(con, "cv_zeropoints",
            """CREATE TABLE {t} (series_key TEXT, frame_id INTEGER,
               bjd_tdb REAL, zp REAL, zp_err REAL, n_star_used INTEGER)""",
            zp_rows)
    swap_in(con, "cv_lightcurve",
            """CREATE TABLE {t} (series_key TEXT, star_id INTEGER, role TEXT,
               frame_id INTEGER, bjd_tdb REAL, inst_mag REAL,
               inst_mag_err REAL, zp REAL, mag REAL, saturated INTEGER)""",
            lc_rows)
    swap_in(con, "cv_series",
            """CREATE TABLE {t} (series_key TEXT PRIMARY KEY,
               target_key TEXT, era_id INTEGER, filter TEXT, provenance TEXT,
               n_staged INTEGER, n_frames_used INTEGER, match_rate REAL,
               n_saturated_dets INTEGER,
               n_target_points INTEGER, n_target_saturated INTEGER,
               n_comp INTEGER, n_check INTEGER,
               n_dropped_unstable INTEGER, comp_passes INTEGER,
               ens_niter INTEGER,
               ens_converged INTEGER, comp_rms_median REAL,
               check_rms_median REAL, chi2_inflation REAL, zp_std REAL,
               target_star_id INTEGER, check_verdict TEXT,
               status TEXT, note TEXT)""", series_rows)
    con.execute("""CREATE INDEX IF NOT EXISTS idx_lc_series
                   ON cv_lightcurve (series_key, star_id, bjd_tdb)""")
    con.commit()
    meta_write(con)
    con.close()


# ---------------------------------------------------------------------------
# Stage: errors
# ---------------------------------------------------------------------------
def cmd_errors(args) -> None:
    """Empirical error model per series, plus one Allan ladder per target.

    Every number here comes from CHECK stars — stars held OUT of the
    zero-point solve.  A comparison star cannot validate the solution it
    helped define, and a polar cannot validate its own error bars because
    its variability is the signal.
    """
    con = connect(args.db)
    model_rows = []
    for skey, in con.execute("SELECT series_key FROM cv_series "
                             "WHERE status='solved' ORDER BY 1"):
        checks = con.execute("""SELECT mean_mag, rms, chi2nu FROM cv_stars
                                WHERE series_key=? AND role='check'
                                  AND rms IS NOT NULL""", (skey,)).fetchall()
        rms = np.array([c[1] for c in checks], dtype=float)
        chi = np.array([c[2] for c in checks if c[2] is not None], dtype=float)
        infl = err.inflation_factor(chi)
        # The photon-noise floor the formal errors PREDICT for those same
        # stars, so the report can put the claim and the outcome on one
        # plot without a hand-typed number in between.
        pred = con.execute("""SELECT avg(inst_mag_err) FROM cv_lightcurve
                              WHERE series_key=? AND role='check'""",
                           (skey,)).fetchone()[0]
        model_rows.append((skey, len(checks),
                           fnum(np.min(rms)) if rms.size else None,
                           fnum(np.median(rms)) if rms.size else None,
                           fnum(np.median(chi)) if chi.size else None,
                           fnum(infl), fnum(pred)))
    swap_in(con, "cv_error_model",
            """CREATE TABLE {t} (series_key TEXT PRIMARY KEY,
               n_check INTEGER, check_rms_min REAL, check_rms_med REAL,
               chi2nu_med REAL, inflation REAL, predicted_err_mean REAL)""",
            model_rows)

    # ---- Allan deviation: the best-observed check star of each TARGET,
    # on its longest uninterrupted run.  One exhibit per target keeps the
    # correlated-noise question answerable per camera era without
    # producing hundreds of near-identical ladders.
    allan_rows = []
    for (tk,) in con.execute("SELECT DISTINCT target_key FROM cv_series "
                             "WHERE status='solved' ORDER BY 1"):
        best = con.execute("""
            SELECT l.series_key, l.star_id, f.night, count(*) n
            FROM cv_lightcurve l
            JOIN cv_frames f ON f.frame_id = l.frame_id
            WHERE l.role='check' AND l.saturated=0 AND l.mag IS NOT NULL
              AND f.target_key=?
            GROUP BY l.series_key, l.star_id, f.night
            ORDER BY n DESC LIMIT 1""", (tk,)).fetchone()
        if not best:
            continue
        skey, sid, night, _n = best
        pts = con.execute("""SELECT l.bjd_tdb, l.mag FROM cv_lightcurve l
            JOIN cv_frames f ON f.frame_id = l.frame_id
            WHERE l.series_key=? AND l.star_id=? AND f.night=?
              AND l.saturated=0 AND l.mag IS NOT NULL
            ORDER BY l.bjd_tdb""", (skey, sid, night)).fetchall()
        t = np.array([p[0] for p in pts], dtype=float)
        m = np.array([p[1] for p in pts], dtype=float)
        a, b = err.longest_run(t)
        t, m = t[a:b], m[a:b]
        if len(m) < 8:
            continue
        dt = float(np.median(np.diff(t)) * 86400.0)
        taus, adevs, npairs = err.allan_deviation(m, dt)
        allan_rows.extend([(tk, skey, int(sid), night, float(x), float(y),
                            int(n)) for x, y, n in zip(taus, adevs, npairs)])
        print(f"  allan {tk}: {skey} star {sid} night {night}: {len(m)} pts, "
              f"cadence {dt:.0f}s, {len(taus)} rungs", flush=True)
    swap_in(con, "cv_allan",
            """CREATE TABLE {t} (target_key TEXT, series_key TEXT,
               star_id INTEGER, night TEXT, tau_s REAL, adev_mag REAL,
               n_pairs INTEGER)""", allan_rows)
    meta_write(con)
    for r in model_rows:
        print("  error_model", r, flush=True)
    con.close()


# ---------------------------------------------------------------------------
# Stage: status
# ---------------------------------------------------------------------------
def cmd_status(args) -> None:
    """Progress at a glance; safe to run at any time, including mid-batch."""
    con = connect(args.db)
    print("frames by status:")
    for st, n in con.execute("""SELECT status, count(*) FROM cv_frames
                                GROUP BY 1 ORDER BY 2 DESC"""):
        print(f"  {st:>22s}: {n}")
    print("\nper (target, era):")
    for row in con.execute("""SELECT target_key, era_id, provenance,
            count(*) n,
            sum(status='matched') matched,
            sum(status='extracted') extracted,
            sum(status='pending') pending,
            sum(status LIKE 'failed%') failed,
            sum(status='excluded' OR status='skipped_pointing') excluded
        FROM cv_frames GROUP BY 1,2,3 ORDER BY 1,2"""):
        print(f"  {row[0]:6s} era{row[1]:<3d} {row[2] or '-':15s} "
              f"n={row[3]:5d} matched={row[4]:5d} extracted={row[5]:5d} "
              f"pending={row[6]:5d} failed={row[7]:4d} excluded={row[8]:4d}")
    try:
        print("\nregistration methods:")
        for m, n in con.execute("""SELECT coalesce(reg_method,'-'), count(*)
                                   FROM cv_frames WHERE status='matched'
                                   GROUP BY 1 ORDER BY 2 DESC"""):
            print(f"  {m:>22s}: {n}")
    except sqlite3.OperationalError:
        pass
    for tbl in ("cv_series", "cv_stars", "cv_lightcurve", "cv_zeropoints",
                "cv_error_model", "cv_allan"):
        try:
            n = con.execute(f"SELECT count(*) FROM {tbl}").fetchone()[0]
            print(f"{tbl}: {n} rows")
        except sqlite3.OperationalError:
            print(f"{tbl}: not built yet")
    con.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--db", type=Path, default=DEFAULT_DB)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    p.add_argument("--wcs-root", type=Path, default=DEFAULT_WCS_ROOT)
    p.add_argument("--recon-dir", type=Path, default=DEFAULT_RECON_DIR)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init")
    for name in ("extract", "match"):
        sp = sub.add_parser(name)
        sp.add_argument("--limit", type=int, default=1000,
                        help="max frames this invocation (chunked runs)")
        sp.add_argument("--workers", type=int, default=4,
                        help=f"worker processes (clamped to {MAX_WORKERS})")
    fp = sub.add_parser("field")
    fp.add_argument("--force", action="store_true",
                    help="recompute ties that already succeeded")
    fp.add_argument("--gaia-identities", action="store_true",
                    help="also query Gaia for star identities on eras whose "
                         "reference is already plate-solved (optional; the "
                         "tie itself does not need it)")
    for name in ("ensemble", "errors", "status"):
        sub.add_parser(name)
    args = p.parse_args()
    {"init": cmd_init, "extract": cmd_extract, "match": cmd_match,
     "field": cmd_field, "ensemble": cmd_ensemble, "errors": cmd_errors,
     "status": cmd_status}[args.cmd](args)


if __name__ == "__main__":
    main()
