#!/usr/bin/env python
"""Build the S4 ensemble-photometry prototype (AN UMa + VV Pup).

WHAT THIS SCRIPT DOES (stage S4 of the shared pipeline)
-------------------------------------------------------
Reads the S0 manifest, measures every server-REDUCED frame of the two
reduction-complete polars, matches stars across frames WITHOUT assuming a
WCS, solves a Honeycutt (1992) inhomogeneous ensemble per (target, era,
filter), ties the field to Gaia DR3, and derives the empirical error model
(the S5 seed).  Everything lands in ONE photometry database:

    products/phot/anuma_vvpup_prototype.sqlite

Tables written (the manifest DB is READ-ONLY to this script):

* ``phot_frames``      — one row per selected science frame: provenance
                         (raw + reduced path), extraction stats, alignment
                         QC, and the frame's ensemble zero point.
* ``phot_detections``  — every source detection on every frame, with its
                         aperture flux and (after matching) its star_id.
* ``phot_ref``         — the reference frame per (target, era) + match tol.
* ``phot_ref_stars``   — the per-(target, era) star catalog: reference
                         positions, and Gaia identity once tied.
* ``phot_gaia_tie``    — the Gaia fit per (target, era): parity, fitted
                         plate scale, match count, the target's star_id.
* ``phot_stars``       — ensemble results per (target, era, filter, star):
                         role (comp/check/target/...), mean mag, RMS, chi2.
* ``phot_series``      — one row per (target, era, filter) ensemble:
                         convergence, comp census, Gaia offset.
* ``phot_selection``   — the frame-selection ledger (raw vs linked vs
                         selected) so excluded frames are counted, not lost.
* ``phot_error_model`` — per (era, filter): check-star RMS floor, median
                         chi2, the error inflation factor.
* ``phot_allan``       — Allan deviation ladder of one long night of one
                         check star (the correlated-noise exhibit).
* ``s4_build_meta``    — timestamp, code version, git commit, constants.

RESUMABLE STAGES (each its own subcommand; run them in this order)
------------------------------------------------------------------
    init                 create/refresh the DB and the frame worklist
    extract [--limit N]  measure N pending frames (repeat until 0 left)
    match   [--limit N]  pick references; align + star-match N frames
    gaia                 tie each (target, era) field to Gaia DR3
    ensemble             solve every (target, era, filter) ensemble
    errors               build the empirical error model + Allan ladder
    report               render docs/pipeline/s4_photometry.html + figures
    status               one-line progress summary (safe anytime)

``extract`` and ``match`` are chunked and idempotent: each invocation
processes at most ``--limit`` unprocessed frames and commits per frame, so
a killed run resumes exactly where it stopped.  The derived stages
(``gaia`` onward) rebuild their tables from scratch each run via the
temp-table + atomic-swap pattern (the S0b discipline).

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s4_photometry.py init
    # repeat until 'pending 0':
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s4_photometry.py extract --limit 400
    ... match --limit 800 ... gaia ... ensemble ... errors ... report
"""

from __future__ import annotations

import argparse
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

# Make the pipeline package importable no matter where the script is invoked
# from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import S4_CODE_VERSION                       # noqa: E402
from macro_phot import photometry as ph                      # noqa: E402
from macro_phot import ensemble as ens                       # noqa: E402
from macro_phot import errors as err                         # noqa: E402
from macro_phot import extract as ext                        # noqa: E402
from macro_phot import gaia as gg                            # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare commands Just Work).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "anuma_vvpup_prototype.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: The prototype targets (decision 2026-08-18: reduction-complete polars).
TARGETS = ("anuma", "vvpup")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def connect(db_path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(db_path)
    con.execute("PRAGMA journal_mode=WAL")   # safe interleaved read/write
    return con


def git_commit() -> str:
    """HEAD's short hash, with an honest '-dirty' marker.

    A bare hash claims 'this tree reproduces the DB'.  When uncommitted or
    untracked files exist that claim is false — the first S4 build recorded
    a commit that contained none of the S4 code — so any non-empty
    ``git status --porcelain`` appends '-dirty' and the report renders the
    weakness instead of hiding it.
    """
    try:
        head = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT,
            capture_output=True, text=True, check=True).stdout.strip()
        porcelain = subprocess.run(
            ["git", "status", "--porcelain"], cwd=REPO_ROOT,
            capture_output=True, text=True, check=True).stdout.strip()
        return f"{head}-dirty" if porcelain else head
    except Exception:
        return ""


def meta_write(con: sqlite3.Connection, extra: dict = ()) -> None:
    """Upsert build metadata (kept current by every stage that writes)."""
    con.execute("""CREATE TABLE IF NOT EXISTS s4_build_meta
                   (key TEXT PRIMARY KEY, value TEXT)""")
    rows = {
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "code_version": S4_CODE_VERSION,
        "git_commit": git_commit(),
        "aperture_radius_arcsec": str(ph.APERTURE_RADIUS_ARCSEC),
        "sky_annulus_arcsec": str(ph.SKY_ANNULUS_ARCSEC),
        "detect_sigma": str(ph.DETECT_SIGMA),
        "clip_sigma": str(ens.CLIP_SIGMA),
        "weight_floor_mag": str(ens.WEIGHT_FLOOR_MAG),
        "ref_doubled_max_frac": str(ph.REF_DOUBLED_MAX_FRAC),
        # RANSAC seed policy: astroalign's shuffle is routed through
        # macro_phot.extract.seeded_ransac with these seeds, so the whole
        # match/gaia geometry reproduces bit-for-bit on re-run.
        "align_seed_policy": ("frame_id for match; ref_frame_id for gaia; "
                              "via macro_phot.extract.seeded_ransac"),
        **dict(extra),
    }
    con.executemany("INSERT OR REPLACE INTO s4_build_meta VALUES (?,?)",
                    rows.items())
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
    con.commit()          # close the implicit insert transaction first
    con.execute("BEGIN")
    con.execute(f"DROP TABLE IF EXISTS {name}")
    con.execute(f"ALTER TABLE {tmp} RENAME TO {name}")
    con.commit()


# ---------------------------------------------------------------------------
# Stage: init
# ---------------------------------------------------------------------------
def cmd_init(args) -> None:
    """Create the DB and the frame worklist from the S0 manifest.

    Selection rule (recorded per (target, era) in ``phot_selection``):
    canonical, error-free Light frames of the prototype targets, outside
    the reduced tree, that HAVE a reduced counterpart through
    ``raw_reduced_links`` — the S4 prototype reads server-reduced pixels
    only and never mixes reduction provenances within a (target, era)
    series.  Existing extraction work is preserved on re-run (INSERT OR
    IGNORE), so init is safe to repeat after a manifest refresh.
    """
    args.db.parent.mkdir(parents=True, exist_ok=True)
    mcon = sqlite3.connect(f"file:{args.manifest}?mode=ro", uri=True)
    con = connect(args.db)
    con.execute("""CREATE TABLE IF NOT EXISTS phot_frames (
        frame_id INTEGER PRIMARY KEY,      -- raw obs_rowid from the manifest
        target_key TEXT, era_id INTEGER, filter TEXT, night TEXT,
        jd REAL, exptime REAL, airmass REAL,
        raw_path TEXT, reduced_path TEXT,
        pointing_offset_deg REAL, qc_flags TEXT,   -- carried from S0
        status TEXT DEFAULT 'pending',     -- pending/extracted/matched/failed_*
        plate_scale REAL, aper_px REAL, egain REAL, pltsolvd INTEGER,
        n_detected INTEGER, bkg_adu REAL, bkg_rms REAL, fwhm_px REAL,
        ali_nmatch INTEGER, ali_rms_px REAL, ali_scale REAL, ali_rot_deg REAL,
        zp REAL, zp_err REAL, n_star_used INTEGER)""")
    con.execute("""CREATE TABLE IF NOT EXISTS phot_detections (
        frame_id INTEGER, det_id INTEGER,
        x REAL, y REAL, flux REAL, fluxerr REAL, fwhm REAL, peak REAL,
        flag INTEGER, clipped INTEGER,
        star_id INTEGER,                   -- ref-star id once matched
        PRIMARY KEY (frame_id, det_id))""")
    con.execute("""CREATE INDEX IF NOT EXISTS idx_det_star
                   ON phot_detections (star_id)""")

    # The selection query mirrors the S0b canonical-science definition.
    sel = mcon.execute("""
        SELECT f.obs_rowid, f.target_key, f.era_id, f.filter, f.night,
               f.jd, f.exptime, f.airmass, f.path, l.reduced_path,
               f.pltsolvd, f.pointing_offset_deg, f.qc_flags
        FROM frames f
        JOIN raw_reduced_links l ON l.raw_rowid = f.obs_rowid
        WHERE f.target_key IN (%s)
          AND f.is_canonical = 1 AND f.error IS NULL AND f.tree != 'reduced'
          AND (f.imagetyp LIKE 'Light%%' OR f.imagetyp IS NULL
               OR f.imagetyp = '')
        ORDER BY f.obs_rowid""" % ",".join("?" * len(TARGETS)),
        TARGETS).fetchall()
    n_new = 0
    for row in sel:
        (fid, tk, era, filt, night, jd, expt, am, rpath, redpath, solved,
         poff, qcf) = row
        cur = con.execute("""INSERT OR IGNORE INTO phot_frames
            (frame_id, target_key, era_id, filter, night, jd, exptime,
             airmass, raw_path, reduced_path, pltsolvd,
             pointing_offset_deg, qc_flags)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (fid, tk, era, filt, night, jd, expt, am, rpath, redpath,
             1 if solved else 0, poff, qcf))
        n_new += cur.rowcount
        if cur.rowcount == 0:
            # Existing row (resumed build): refresh the carried S0 QC
            # columns, which cost nothing and may have improved upstream.
            con.execute("""UPDATE phot_frames SET pointing_offset_deg=?,
                           qc_flags=? WHERE frame_id=?""", (poff, qcf, fid))
    con.commit()

    # Selection ledger: raw canonical vs reduced-linked, per (target, era) —
    # the report shows exactly what the provenance rule excluded.
    ledger = mcon.execute("""
        SELECT f.target_key, f.era_id,
               count(*) AS n_raw,
               count(l.reduced_rowid) AS n_linked
        FROM frames f
        LEFT JOIN raw_reduced_links l ON l.raw_rowid = f.obs_rowid
        WHERE f.target_key IN (%s)
          AND f.is_canonical = 1 AND f.error IS NULL AND f.tree != 'reduced'
          AND (f.imagetyp LIKE 'Light%%' OR f.imagetyp IS NULL
               OR f.imagetyp = '')
        GROUP BY f.target_key, f.era_id""" % ",".join("?" * len(TARGETS)),
        TARGETS).fetchall()
    swap_in(con, "phot_selection",
            """CREATE TABLE {t} (target_key TEXT, era_id INTEGER,
               n_raw INTEGER, n_linked INTEGER)""",
            [tuple(r) for r in ledger])
    meta_write(con, {"manifest_path": str(args.manifest)})
    total = con.execute("SELECT count(*) FROM phot_frames").fetchone()[0]
    print(f"init: worklist holds {total} frames ({n_new} new)")
    con.close(); mcon.close()


# ---------------------------------------------------------------------------
# Stage: extract (chunked, resumable)
# ---------------------------------------------------------------------------
def cmd_extract(args) -> None:
    """Measure up to --limit pending frames; commit per frame (resumable)."""
    con = connect(args.db)
    todo = con.execute("""SELECT frame_id, reduced_path, exptime
                          FROM phot_frames WHERE status = 'pending'
                          ORDER BY frame_id LIMIT ?""",
                       (args.limit,)).fetchall()
    t0 = time.time()
    n_done = n_fail = 0
    for fid, redpath, exptime in todo:
        path = args.archive / redpath
        try:
            data, meta = ext.read_reduced(path)
            scale = ph.plate_scale_arcsec_per_px(meta["xpixsz"],
                                                 meta["focallen"])
            aper = ph.aperture_radius_px(scale)
            if aper is None:
                raise ValueError("no plate scale in header")
            stats, dets = ext.measure_frame(data, meta["egain"], aper)
            con.execute("BEGIN")
            con.execute("""UPDATE phot_frames SET status='extracted',
                plate_scale=?, aper_px=?, egain=?, n_detected=?,
                bkg_adu=?, bkg_rms=?, fwhm_px=? WHERE frame_id=?""",
                (scale, aper, meta["egain"], stats["n_detected"],
                 stats["bkg_adu"], stats["bkg_rms"], stats["fwhm_px"], fid))
            con.execute("DELETE FROM phot_detections WHERE frame_id=?",
                        (fid,))
            con.executemany("""INSERT INTO phot_detections
                (frame_id, det_id, x, y, flux, fluxerr, fwhm, peak, flag,
                 clipped, star_id)
                VALUES (?,?,?,?,?,?,?,?,?,?,NULL)""",
                [(fid, i, float(dets["x"][i]), float(dets["y"][i]),
                  float(dets["flux"][i]), float(dets["fluxerr"][i]),
                  float(dets["fwhm"][i]), float(dets["peak"][i]),
                  int(dets["flag"][i]), int(dets["clipped"][i]))
                 for i in range(len(dets["x"]))])
            con.commit()
            n_done += 1
        except Exception as e:            # a bad file must not stall the run
            con.rollback()
            con.execute("""UPDATE phot_frames SET status=?
                           WHERE frame_id=?""",
                        (f"failed_extract:{type(e).__name__}", fid))
            con.commit()
            n_fail += 1
    left = con.execute("SELECT count(*) FROM phot_frames "
                       "WHERE status='pending'").fetchone()[0]
    print(f"extract: {n_done} ok, {n_fail} failed, "
          f"{time.time() - t0:.0f}s; pending {left}")
    con.close()


# ---------------------------------------------------------------------------
# Stage: match (reference choice + astroalign, chunked)
# ---------------------------------------------------------------------------
def _ensure_reference(con, target_key: str, era_id: int) -> tuple:
    """Choose (once), QC, and return the reference for one (target, era).

    Candidates are walked in :func:`macro_phot.photometry.rank_references`
    order and each must pass DOUBLE-IMAGE quality control before adoption:
    a guiding jump doubles every star into an equal-brightness pair, which
    INFLATES the detection count — so 'most detections among the sharpest'
    actively rewarded doubled frames, and the first build adopted doubled
    references for both VV Pup series (83-89% of catalog stars paired,
    vs 0.7% on a clean frame).  The pure detector
    (:func:`macro_phot.photometry.paired_fraction`) rejects any candidate
    with more than ``REF_DOUBLED_MAX_FRAC`` of its stars carrying a
    similar-flux companion within ``REF_PAIR_RADIUS_FWHM`` x FWHM; the
    accepted fraction and the number of rejected candidates are recorded
    in ``phot_ref`` as the audit trail.
    """
    con.execute("""CREATE TABLE IF NOT EXISTS phot_ref (
        target_key TEXT, era_id INTEGER, ref_frame_id INTEGER,
        n_stars INTEGER, fwhm_px REAL, tol_px REAL,
        doubled_frac REAL, n_cand_rejected INTEGER,
        PRIMARY KEY (target_key, era_id))""")
    con.execute("""CREATE TABLE IF NOT EXISTS phot_ref_stars (
        target_key TEXT, era_id INTEGER, star_id INTEGER,
        x REAL, y REAL, flux REAL,
        gaia_source_id INTEGER, gaia_gmag REAL, ra_deg REAL, dec_deg REAL,
        PRIMARY KEY (target_key, era_id, star_id))""")
    row = con.execute("""SELECT ref_frame_id, tol_px FROM phot_ref
                         WHERE target_key=? AND era_id=?""",
                      (target_key, era_id)).fetchone()
    if row:
        ref_id, tol = row
    else:
        stats = con.execute("""SELECT frame_id, n_detected, fwhm_px, pltsolvd
                               FROM phot_frames
                               WHERE target_key=? AND era_id=?
                                 AND status IN ('extracted','matched')""",
                            (target_key, era_id)).fetchall()
        ranking = ph.rank_references([tuple(s) for s in stats])
        fwhm_of = {s[0]: s[2] for s in stats}
        ref_id = None
        doubled = None
        n_rejected = 0
        for cand in ranking:
            cdets = con.execute("""SELECT x, y, flux FROM phot_detections
                                   WHERE frame_id=? ORDER BY det_id""",
                                (cand,)).fetchall()
            cxy = np.array([[d[0], d[1]] for d in cdets], dtype=float)
            cfl = np.array([d[2] for d in cdets], dtype=float)
            frac = ph.paired_fraction(
                cxy, cfl, ph.REF_PAIR_RADIUS_FWHM * float(fwhm_of[cand]))
            if frac <= ph.REF_DOUBLED_MAX_FRAC:
                ref_id, doubled = cand, frac
                break
            n_rejected += 1
            print(f"reference QC: {target_key}/era{era_id} frame {cand} "
                  f"REJECTED — paired fraction {frac:.2f} "
                  f"(> {ph.REF_DOUBLED_MAX_FRAC:g}: double-imaged)")
        if ref_id is None:
            return None, None, None
        fwhm = fwhm_of[ref_id]
        tol = ph.match_tolerance_px(fwhm)
        # The reference star catalog = the reference frame's own detections
        # (star_id := its det_id — a stable, auditable identity).
        dets = con.execute("""SELECT det_id, x, y, flux FROM phot_detections
                              WHERE frame_id=? ORDER BY det_id""",
                           (ref_id,)).fetchall()
        con.executemany("""INSERT OR REPLACE INTO phot_ref_stars
            (target_key, era_id, star_id, x, y, flux)
            VALUES (?,?,?,?,?,?)""",
            [(target_key, era_id, d, x, y, fl) for d, x, y, fl in dets])
        con.execute("""INSERT OR REPLACE INTO phot_ref
                       VALUES (?,?,?,?,?,?,?,?)""",
                    (target_key, era_id, ref_id, len(dets), fwhm, tol,
                     doubled, n_rejected))
        # The reference matches itself by identity — mark it 'matched' so
        # its own measurements join the ensemble like everyone else's.
        con.execute("""UPDATE phot_detections SET star_id = det_id
                       WHERE frame_id=?""", (ref_id,))
        con.execute("""UPDATE phot_frames SET status='matched', ali_nmatch=?,
                       ali_rms_px=0.0, ali_scale=1.0, ali_rot_deg=0.0
                       WHERE frame_id=?""", (len(dets), ref_id))
        con.commit()
    stars = con.execute("""SELECT star_id, x, y FROM phot_ref_stars
                           WHERE target_key=? AND era_id=?
                           ORDER BY star_id""",
                        (target_key, era_id)).fetchall()
    return ref_id, tol, stars


def cmd_match(args) -> None:
    """Align up to --limit extracted frames to their series reference."""
    import astroalign as aa
    con = connect(args.db)
    series = con.execute("""SELECT DISTINCT target_key, era_id
                            FROM phot_frames ORDER BY 1, 2""").fetchall()
    t0 = time.time()
    n_done = n_fail = 0
    budget = args.limit
    for tk, era in series:
        if budget <= 0:
            break
        ref_id, tol, stars = _ensure_reference(con, tk, era)
        if ref_id is None:
            continue
        ref_ids = np.array([s[0] for s in stars])
        ref_xy = np.array([[s[1], s[2]] for s in stars], dtype=float)
        # Bright reference pool for the triangle fit: the reference
        # frame's own UNCLIPPED detections, flux-descending — exactly the
        # rule applied to the frame side.  Symmetry is load-bearing:
        # astroalign's control points are the FIRST N entries, and a
        # saturated star kept on one side but excluded on the other
        # poisons the control sets from the very top (the VV Pup dense
        # field failed wholesale on exactly this asymmetry).
        max_pool = max(p for p, _ in ext.ALIGN_ATTEMPTS)
        ref_id_row = con.execute("""SELECT ref_frame_id FROM phot_ref
            WHERE target_key=? AND era_id=?""", (tk, era)).fetchone()
        bright_ref = np.array([[r[0], r[1]] for r in con.execute(
            """SELECT x, y FROM phot_detections
               WHERE frame_id=? AND clipped=0
               ORDER BY flux DESC LIMIT ?""",
            (ref_id_row[0], max_pool))], dtype=float)
        # Fast-skip the frames S0 already proved are pointed elsewhere: a
        # 10-degree offset field shares no stars with the reference, and
        # astroalign burns ~15 s exhausting triangles to say so.  The skip
        # is recorded as its own status — visible, never silent.
        con.execute("""UPDATE phot_frames
                       SET status='skipped_pointing'
                       WHERE target_key=? AND era_id=? AND status='extracted'
                         AND qc_flags LIKE '%pointing_gt1deg%'""", (tk, era))
        con.commit()
        todo = con.execute("""SELECT frame_id FROM phot_frames
                              WHERE target_key=? AND era_id=?
                                AND status='extracted'
                              ORDER BY frame_id LIMIT ?""",
                           (tk, era, budget)).fetchall()
        for (fid,) in todo:
            dets = con.execute("""SELECT det_id, x, y, flux, clipped
                                  FROM phot_detections WHERE frame_id=?
                                  ORDER BY det_id""", (fid,)).fetchall()
            xy = np.array([[d[1], d[2]] for d in dets], dtype=float)
            try:
                if len(dets) < ph.MIN_STARS_FOR_ALIGN:
                    raise ValueError("too few detections to align")
                bright = np.array(
                    [[d[1], d[2]] for d in sorted(
                        (d for d in dets if not d[4]),
                        key=lambda d: -d[3])[:max_pool]])
                # seed = frame_id: pins astroalign's RANSAC so identical
                # inputs give identical transforms (and star identities)
                # on every re-run — the regenerable-products discipline.
                tf = ext.find_series_transform(bright, bright_ref, seed=fid)
                moved = tf(xy)            # ALL detections into ref pixels
                idx = ph.match_one_to_one(ref_xy, moved, tol)
                matched = idx >= 0
                # Alignment residual: matched pairs only.
                if matched.any():
                    d2 = np.hypot(
                        moved[matched, 0] - ref_xy[idx[matched], 0],
                        moved[matched, 1] - ref_xy[idx[matched], 1])
                    rms = float(np.sqrt(np.mean(d2 ** 2)))
                else:
                    rms = None
                con.execute("BEGIN")
                con.executemany(
                    """UPDATE phot_detections SET star_id=?
                       WHERE frame_id=? AND det_id=?""",
                    [(int(ref_ids[idx[i]]), fid, int(dets[i][0]))
                     for i in range(len(dets)) if matched[i]])
                con.execute("""UPDATE phot_frames SET status='matched',
                    ali_nmatch=?, ali_rms_px=?, ali_scale=?, ali_rot_deg=?
                    WHERE frame_id=?""",
                    (int(matched.sum()), rms, float(tf.scale),
                     float(np.degrees(tf.rotation)), fid))
                con.commit()
                n_done += 1
            except Exception as e:
                con.rollback()
                con.execute("UPDATE phot_frames SET status=? "
                            "WHERE frame_id=?",
                            (f"failed_match:{type(e).__name__}", fid))
                con.commit()
                n_fail += 1
            budget -= 1
            if budget <= 0:
                break
    left = con.execute("SELECT count(*) FROM phot_frames "
                       "WHERE status='extracted'").fetchone()[0]
    print(f"match: {n_done} ok, {n_fail} failed, "
          f"{time.time() - t0:.0f}s; unmatched-extracted {left}")
    con.close()


# ---------------------------------------------------------------------------
# Stage: gaia
# ---------------------------------------------------------------------------
def cmd_gaia(args) -> None:
    """Tie each (target, era) reference catalog to Gaia DR3."""
    con = connect(args.db)
    refs = con.execute("""SELECT r.target_key, r.era_id, r.ref_frame_id,
                                 f.plate_scale, f.fwhm_px
                          FROM phot_ref r
                          JOIN phot_frames f ON f.frame_id = r.ref_frame_id
                          ORDER BY 1, 2""").fetchall()
    tie_rows = []
    for tk, era, ref_id, scale, fwhm in refs:
        ra0, dec0, coord_src = gg.resolve_target(tk)
        stars = con.execute("""SELECT star_id, x, y, flux FROM phot_ref_stars
                               WHERE target_key=? AND era_id=?
                               ORDER BY star_id""", (tk, era)).fetchall()
        star_ids = np.array([s[0] for s in stars])
        ref_xy = np.array([[s[1], s[2]] for s in stars], dtype=float)
        # Unclipped flux-descending bright pool — the SAME construction the
        # frame matcher uses (symmetric pools; saturated stars excluded).
        max_pool = max(p for p, _ in ext.ALIGN_ATTEMPTS)
        ref_bright = np.array([[r[0], r[1]] for r in con.execute(
            """SELECT x, y FROM phot_detections
               WHERE frame_id=? AND clipped=0
               ORDER BY flux DESC LIMIT ?""", (ref_id, max_pool))],
            dtype=float)
        # Cone radius: half the frame diagonal at the header plate scale,
        # padded 20% for pointing scatter about the target.
        n1, n2 = con.execute("""SELECT max(x), max(y) FROM phot_detections
                                WHERE frame_id=?""", (ref_id,)).fetchone()
        radius = 1.2 * scale * float(np.hypot(n1, n2)) / 2.0 / 3600.0
        # The FIT pool is confined to the frame's inscribed circle — the
        # cone is wider than the footprint by construction, and bright
        # off-frame stars would crowd the triangle fit with unmatchables.
        fit_radius = scale * float(min(n1, n2)) / 2.0
        try:
            gaia = gg.cone_query(ra0, dec0, radius)
        except Exception as e:
            print(f"gaia: cone query failed for {tk}/era{era}: {e}")
            tie_rows.append((tk, era, None, None, None, 0, 0, None,
                             None, ra0, dec0, coord_src, "query_failed"))
            continue
        # seed = the reference frame id (same policy as the match stage):
        # the Gaia tie must reproduce bit-for-bit on re-run.
        fit = gg.identify_reference(ref_xy, gaia, ra0, dec0,
                                    ref_bright_xy=ref_bright,
                                    fit_radius_arcsec=fit_radius,
                                    seed=ref_id)
        if fit is None:
            tie_rows.append((tk, era, None, None, None, len(gaia["ra"]), 0,
                             None, None, ra0, dec0, coord_src, "fit_failed"))
            continue
        # Write the Gaia identity onto the reference star catalog.
        idx = fit["gaia_idx"]
        con.execute("BEGIN")
        for i, sid in enumerate(star_ids):
            if idx[i] >= 0:
                con.execute("""UPDATE phot_ref_stars
                    SET gaia_source_id=?, gaia_gmag=?, ra_deg=?, dec_deg=?
                    WHERE target_key=? AND era_id=? AND star_id=?""",
                    (int(gaia["source_id"][idx[i]]),
                     float(gaia["gmag"][idx[i]]),
                     float(fit["ref_radec"][i, 0]),
                     float(fit["ref_radec"][i, 1]), tk, era, int(sid)))
            else:
                con.execute("""UPDATE phot_ref_stars
                    SET ra_deg=?, dec_deg=?
                    WHERE target_key=? AND era_id=? AND star_id=?""",
                    (float(fit["ref_radec"][i, 0]),
                     float(fit["ref_radec"][i, 1]), tk, era, int(sid)))
        con.commit()
        # The TARGET star: reference star nearest the resolved coordinates
        # (tolerance = the Gaia match tolerance; a polar in a deep low
        # state may genuinely be absent from the reference frame).
        sep_deg = np.hypot(
            (fit["ref_radec"][:, 0] - ra0)
            * np.cos(np.radians(dec0)), fit["ref_radec"][:, 1] - dec0)
        j = int(np.argmin(sep_deg))
        target_sep = float(sep_deg[j] * 3600.0)
        target_sid = (int(star_ids[j])
                      if target_sep <= gg.TARGET_ID_TOL_ARCSEC else None)
        tie_rows.append((tk, era, fit["parity"],
                         fit["scale_arcsec_per_px"], fit["rot_deg"],
                         len(gaia["ra"]), fit["n_matched"], target_sid,
                         target_sep, ra0, dec0, coord_src, "ok"))
        print(f"gaia: {tk}/era{era}: {fit['n_matched']} matched of "
              f"{len(star_ids)} ref stars; scale {fit['scale_arcsec_per_px']:.4f}"
              f" \"/px (header {scale:.4f}); target star "
              f"{target_sid} ({sep_deg[j] * 3600:.2f}\" away)")
    swap_in(con, "phot_gaia_tie",
            """CREATE TABLE {t} (target_key TEXT, era_id INTEGER,
               parity TEXT, scale_fit REAL, rot_deg REAL,
               n_gaia INTEGER, n_matched INTEGER, target_star_id INTEGER,
               target_sep_arcsec REAL,
               target_ra REAL, target_dec REAL, coord_source TEXT,
               status TEXT)""", tie_rows)
    meta_write(con)
    con.close()


# ---------------------------------------------------------------------------
# Stage: ensemble
# ---------------------------------------------------------------------------
def _series_matrix(con, tk: str, era: int, filt: str):
    """Assemble the (stars x frames) magnitude/error matrices of one series.

    Only 'matched' frames of the series' filter contribute; clipped
    (non-linear) detections are dropped at load.  Returns None for an
    empty series.
    """
    frames = con.execute("""SELECT frame_id, exptime, jd, night
                            FROM phot_frames
                            WHERE target_key=? AND era_id=? AND filter=?
                              AND status='matched' ORDER BY jd""",
                         (tk, era, filt)).fetchall()
    if not frames:
        return None
    stars = con.execute("""SELECT star_id FROM phot_ref_stars
                           WHERE target_key=? AND era_id=?
                           ORDER BY star_id""", (tk, era)).fetchall()
    sid_row = {s[0]: i for i, s in enumerate(stars)}
    fid_col = {f[0]: j for j, f in enumerate(frames)}
    S, F = len(stars), len(frames)
    mag = np.full((S, F), np.nan)
    sig = np.full((S, F), np.nan)
    ph_ids = ",".join(str(f[0]) for f in frames)
    rows = con.execute(f"""SELECT frame_id, star_id, flux, fluxerr
                           FROM phot_detections
                           WHERE star_id IS NOT NULL AND clipped = 0
                             AND frame_id IN ({ph_ids})""").fetchall()
    expt = {f[0]: f[1] for f in frames}
    for fid, sid, flux, fluxerr in rows:
        i, j = sid_row[sid], fid_col[fid]
        m = ph.instrumental_mag(np.array([flux]), expt[fid])[0]
        e = ph.mag_error(np.array([flux]), np.array([fluxerr]))[0]
        mag[i, j], sig[i, j] = m, e
    star_ids = np.array([s[0] for s in stars])
    frame_ids = np.array([f[0] for f in frames])
    return star_ids, frame_ids, mag, sig


def cmd_ensemble(args) -> None:
    """Solve every (target, era, filter) Honeycutt ensemble."""
    con = connect(args.db)
    tie = {(r[0], r[1]): r for r in con.execute(
        """SELECT target_key, era_id, target_star_id, status
           FROM phot_gaia_tie""")}
    series = con.execute("""SELECT DISTINCT target_key, era_id, filter
                            FROM phot_frames WHERE status='matched'
                            ORDER BY 1, 2, 3""").fetchall()
    star_rows, series_rows = [], []
    for tk, era, filt in series:
        packed = _series_matrix(con, tk, era, filt)
        if packed is None:
            continue
        star_ids, frame_ids, mag, sig = packed
        target_sid = tie.get((tk, era), (None,) * 4)[2]
        target_row = (int(np.flatnonzero(star_ids == target_sid)[0])
                      if target_sid is not None
                      and (star_ids == target_sid).any() else None)
        sel = ens.select_comps(mag, sig, target_row=target_row)
        sol = sel.solution
        # Per-star statistics against the FINAL fixed zero points — every
        # star, not only comps, so the target and checks get honest stats.
        mean, rms, nobs, chi2nu = ens.star_stats(mag, sig, sol.zp)
        gaia_cols = dict((r[0], (r[1], r[2])) for r in con.execute(
            """SELECT star_id, gaia_source_id, gaia_gmag FROM phot_ref_stars
               WHERE target_key=? AND era_id=?""", (tk, era)))
        for i, sid in enumerate(star_ids):
            if nobs[i] == 0:
                continue                  # never seen in this filter
            gsid, gmag = gaia_cols.get(int(sid), (None, None))
            star_rows.append((tk, era, filt, int(sid), sel.role[i],
                              float(mean[i]) if np.isfinite(mean[i]) else None,
                              float(rms[i]) if np.isfinite(rms[i]) else None,
                              int(nobs[i]),
                              float(chi2nu[i]) if np.isfinite(chi2nu[i])
                              else None, gsid, gmag))
        # Frame zero points land on phot_frames (a frame has one filter,
        # hence exactly one ensemble).  A frame whose comps were all
        # absent/clipped has NaN ZP from the solver and is stored as NULL —
        # every downstream consumer filters on `zp IS NOT NULL`, so its
        # measurements are excluded rather than 'corrected' by a
        # fabricated zero point (the 0.0 bug of the first build).
        con.execute("BEGIN")
        con.executemany("""UPDATE phot_frames SET zp=?, zp_err=?,
                           n_star_used=? WHERE frame_id=?""",
                        [(float(sol.zp[j])
                          if np.isfinite(sol.zp[j]) else None,
                          float(sol.zp_err[j])
                          if np.isfinite(sol.zp_err[j]) else None,
                          int(sol.n_star_used[j]), int(frame_ids[j]))
                         for j in range(len(frame_ids))])
        con.commit()
        # Gaia zero-point offset over comp stars with Gaia G.
        comp_mask = sel.role == "comp"
        g = np.array([gaia_cols.get(int(s), (None, np.nan))[1] or np.nan
                      for s in star_ids])
        off, mad, n_off = gg.median_offset(g[comp_mask], mean[comp_mask])
        roles = sel.role
        series_rows.append((
            tk, era, filt, len(frame_ids), int((roles == "comp").sum()),
            int((roles == "check").sum()),
            int((roles == "dropped_unstable").sum()),
            sol.n_iter, int(sol.converged), sel.comp_rms_median,
            float(np.nanstd(sol.zp)),
            off if np.isfinite(off) else None,
            mad if np.isfinite(mad) else None, n_off,
            "gaia_median_offset" if n_off else "instrumental"))
        print(f"ensemble {tk}/era{era}/{filt}: {len(frame_ids)} frames, "
              f"comps {(roles == 'comp').sum()}, checks "
              f"{(roles == 'check').sum()}, dropped "
              f"{(roles == 'dropped_unstable').sum()}, iters {sol.n_iter}, "
              f"converged {sol.converged}, comp rms {sel.comp_rms_median:.4f}")
    swap_in(con, "phot_stars",
            """CREATE TABLE {t} (target_key TEXT, era_id INTEGER,
               filter TEXT, star_id INTEGER, role TEXT, mean_mag REAL,
               rms REAL, nobs INTEGER, chi2nu REAL,
               gaia_source_id INTEGER, gaia_gmag REAL)""", star_rows)
    swap_in(con, "phot_series",
            """CREATE TABLE {t} (target_key TEXT, era_id INTEGER,
               filter TEXT, n_frames INTEGER, n_comp INTEGER,
               n_check INTEGER, n_dropped INTEGER, ens_niter INTEGER,
               ens_converged INTEGER, comp_rms_median REAL, zp_std REAL,
               gaia_offset REAL, gaia_offset_mad REAL, gaia_offset_n INTEGER,
               zp_source TEXT)""", series_rows)
    meta_write(con)
    con.close()


# ---------------------------------------------------------------------------
# Stage: errors (the S5 seed)
# ---------------------------------------------------------------------------
def cmd_errors(args) -> None:
    """Empirical error model: RMS floors, inflation factors, Allan ladder."""
    con = connect(args.db)
    # ---- inflation + floors per (era, filter), CHECK stars pooled over
    # targets (both polars share era 76; the check stars are field stars,
    # so pooling is legitimate and doubles the statistics).
    combos = con.execute("""SELECT DISTINCT era_id, filter FROM phot_stars
                            ORDER BY 1, 2""").fetchall()
    model_rows = []
    for era, filt in combos:
        checks = con.execute("""SELECT rms, chi2nu FROM phot_stars
            WHERE era_id=? AND filter=? AND role='check'
              AND rms IS NOT NULL""", (era, filt)).fetchall()
        rms = np.array([c[0] for c in checks], dtype=float)
        chi = np.array([c[1] for c in checks if c[1] is not None],
                       dtype=float)
        infl = err.inflation_factor(chi)
        model_rows.append((
            int(era), filt, len(checks),
            float(np.min(rms)) if rms.size else None,
            float(np.median(rms)) if rms.size else None,
            float(np.median(chi)) if chi.size else None,
            infl if np.isfinite(infl) else None))
    swap_in(con, "phot_error_model",
            """CREATE TABLE {t} (era_id INTEGER, filter TEXT,
               n_check INTEGER, check_rms_min REAL, check_rms_med REAL,
               chi2nu_med REAL, inflation REAL)""", model_rows)

    # ---- Allan deviation: the longest uninterrupted run of the
    # best-observed check star, chosen automatically (bright + many points).
    best = con.execute("""
        SELECT s.target_key, s.era_id, s.filter, s.star_id, f.night,
               count(*) AS n
        FROM phot_stars s
        JOIN phot_ref_stars r ON r.target_key = s.target_key
             AND r.era_id = s.era_id AND r.star_id = s.star_id
        JOIN phot_detections d ON d.star_id = s.star_id
        JOIN phot_frames f ON f.frame_id = d.frame_id
             AND f.target_key = s.target_key AND f.era_id = s.era_id
             AND f.filter = s.filter
        WHERE s.role = 'check' AND d.clipped = 0
        GROUP BY s.target_key, s.era_id, s.filter, s.star_id, f.night
        ORDER BY n DESC LIMIT 1""").fetchone()
    allan_rows = []
    if best:
        tk, era, filt, sid, night, _n = best
        pts = con.execute("""
            SELECT f.jd, d.flux, f.exptime, f.zp FROM phot_detections d
            JOIN phot_frames f ON f.frame_id = d.frame_id
            WHERE d.star_id=? AND d.clipped=0 AND f.target_key=?
              AND f.era_id=? AND f.filter=? AND f.night=? AND f.zp IS NOT NULL
            ORDER BY f.jd""", (sid, tk, era, filt, night)).fetchall()
        jd = np.array([p[0] for p in pts])
        m = np.array([ph.instrumental_mag(np.array([p[1]]), p[2])[0] - p[3]
                      for p in pts])
        a, b = err.longest_run(jd)
        jd, m = jd[a:b], m[a:b]
        ok = np.isfinite(m)
        jd, m = jd[ok], m[ok]
        if len(m) >= 8:
            dt = float(np.median(np.diff(jd)) * 86400.0)
            taus, adevs, npairs = err.allan_deviation(m, dt)
            allan_rows = [(tk, int(era), filt, int(sid), night,
                           float(t), float(ad), int(np_))
                          for t, ad, np_ in zip(taus, adevs, npairs)]
            print(f"allan: {tk}/era{era}/{filt} star {sid} night {night}: "
                  f"{len(m)} pts, cadence {dt:.0f}s, "
                  f"{len(taus)} tau rungs")
    swap_in(con, "phot_allan",
            """CREATE TABLE {t} (target_key TEXT, era_id INTEGER,
               filter TEXT, star_id INTEGER, night TEXT, tau_s REAL,
               adev_mag REAL, n_pairs INTEGER)""", allan_rows)
    meta_write(con)
    for r in model_rows:
        print("error_model", r)
    con.close()


# ---------------------------------------------------------------------------
# Stage: report / status
# ---------------------------------------------------------------------------
def cmd_report(args) -> None:
    from macro_phot.report_s4 import render_report
    out = render_report(args.db)
    print(f"report: {out}")


def cmd_status(args) -> None:
    con = connect(args.db)
    for row in con.execute("""SELECT status, count(*) FROM phot_frames
                              GROUP BY status ORDER BY 2 DESC"""):
        print(f"{row[0]:>28}: {row[1]}")
    con.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--db", type=Path, default=DEFAULT_DB,
                   help="photometry database (default: the real one)")
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="S0 manifest database (read-only)")
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE,
                   help="archive root holding the reduced/ tree")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("init")
    for name in ("extract", "match"):
        sp = sub.add_parser(name)
        sp.add_argument("--limit", type=int, default=400,
                        help="max frames this invocation (chunked runs)")
    sub.add_parser("gaia")
    sub.add_parser("ensemble")
    sub.add_parser("errors")
    sub.add_parser("report")
    sub.add_parser("status")
    args = p.parse_args()
    {"init": cmd_init, "extract": cmd_extract, "match": cmd_match,
     "gaia": cmd_gaia, "ensemble": cmd_ensemble, "errors": cmd_errors,
     "report": cmd_report, "status": cmd_status}[args.cmd](args)


if __name__ == "__main__":
    main()
