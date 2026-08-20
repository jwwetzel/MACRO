#!/usr/bin/env python
"""Run the S2 detector-truth campaign against the RLMT archive.

WHAT THIS SCRIPT DOES (stage S2 of the shared pipeline)
-------------------------------------------------------
Reads pixels from the archive for four detector probes, AUGMENTS the S0/S0b
manifest database with new ``s2_*`` tables plus the ``detector_params``
table, writes regenerable pixel products under ``products/detector/``, and
renders the evidence report ``docs/pipeline/s2_detector.html``.  Existing
tables are never modified.

* ``ceiling``      — per-mode science-frame pixel histograms; clip/pileup
                     location; the ceiling + saturation-veto memo numbers.
* ``ptc``          — 2023-06-07 repeated darks + repeated star fields:
                     difference-pair photon transfer (gain, read noise,
                     StackPro variance suppression), amp-glow check.
* ``reconstruct``  — per-era per-pixel fits of raw = F*reduced + D across
                     raw<->reduced pairs: the effective master dark/flat the
                     unaudited reduction actually applied; era 47 graded
                     against its archived master bias/dark.
* ``noise``        — the EMPIRICAL noise model: same-scene consecutive
                     science-frame pairs from every readout mode (not just
                     the one PTC night) turned into a measured
                     counts-vs-variance table per mode.  The CV
                     time-series error model reads this table; it assumes
                     no gain, no Poisson law, no formula.
* ``linearity``    — the 2024-05-20 Vega exposure ladder + every other
                     archival ladder the manifest surfaces: counts-vs-
                     exptime residuals per mode.
* ``params``       — distills all of the above into ``detector_params``
                     (one row per (era_group, quantity) with value,
                     uncertainty, method, provenance).
* ``report``       — renders the S2 evidence page from the database.

RESUMABILITY (the 10-minute-batch discipline)
---------------------------------------------
Every pixel-reading subcommand records finished work in its s2_* table and
skips it on re-invocation, so the campaign is driven as repeated short
calls (``--batch`` caps frames per call).  Interrupting anything mid-batch
loses at most one uncommitted batch.

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/run_s2_campaign.py ceiling --batch 200
    ... (repeat until it reports nothing left to do, same for the others)
    ... ptc / reconstruct / linearity / params / report
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import subprocess
import sys
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits

# Make the pipeline package importable no matter where the script is invoked
# from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core.inventory import exptime_bin                  # noqa: E402
from rlmt_diagnostics import S2_CODE_VERSION                  # noqa: E402
from rlmt_diagnostics import ceiling as ceil                  # noqa: E402
from rlmt_diagnostics import linearity as lin                 # noqa: E402
from rlmt_diagnostics import noise as noisemod                # noqa: E402
from rlmt_diagnostics import ptc                              # noqa: E402
from rlmt_diagnostics import reconstruct as rec               # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
PRODUCTS_DIR = REPO_ROOT / "products" / "detector"

#: Readout modes the ceiling stage samples, with per-mode frame targets.
#: Majors get 150 frames; the sparse modes take what exists.  The blank-
#: READOUTM 2026 frames (current camera, headers unwritten — S0's finding)
#: are sampled as their own explicit group.
CEILING_MODES: dict[str, int] = {
    "High Gain": 150,
    "High Gain StackPro": 150,
    "Low Gain": 120,
    "Mode0": 150,
    "Fast": 150,
    "1MHz High Sensitivity 16-bit": 150,
    "5MHz High Sensitivity 16-bit": 120,
    "(blank 2026)": 120,
}

#: PTC pairing caps: consecutive same-(scene, mode, exptime) frame pairs.
PTC_NIGHT = "2023-06-07"
PTC_MAX_PAIRS_PER_GROUP = 15

#: Reconstruction: eras with at least this many usable raw<->reduced links
#: enter the experiment; at most this many pairs are read per era.
RECON_MIN_LINKS = 20
RECON_MAX_PAIRS = 36

#: How far the archived master dark may be scaled in exposure time before
#: the ground-truth comparison drops the dark term instead of extrapolating
#: it.  4x is generous for a linear dark-current model and still refuses the
#: 3,100x extrapolation era 76's 0.01 s master would have demanded.
RECON_TRUTH_MAX_SCALE = 4.0

#: Empirical noise model: same-scene consecutive science pairs per mode.
#: 16 pairs is enough to give every level bin several independent pairs
#: (the curve builder refuses thin bins) while keeping one mode's pass
#: inside the ten-minute batch discipline.
NOISE_MAX_PAIRS_PER_MODE = 16

#: ... and at most this many from any ONE (night, target, filter, exptime)
#: scene, so a single well-observed night cannot own a mode's curve.
NOISE_MAX_PAIRS_PER_SCENE = 3

#: Two frames are "the same scene seconds apart" only if their JD gap is
#: below this many days (6 minutes).  Longer gaps let the sky rotate, the
#: transparency drift and the target move — all of which inflate the
#: difference variance with things that are not detector noise.
NOISE_MAX_GAP_DAYS = 6.0 / (24.0 * 60.0)

#: A scene must hold at least this many frames before it is paired: two
#: frames are one pair and no cross-check.
NOISE_MIN_SCENE_FRAMES = 3

#: Linearity: cap on auto-discovered ladders (the Vega ladder is always
#: included when present).
MAX_LADDERS = 12

#: A ladder with at least this many rungs may hold a single frame per rung
#: (the archive's dedicated 2023-10 exposure sequences shot each of their
#: 11-14 exposure times exactly once; the rung count itself provides the
#: redundancy a 3-rung ladder gets from repeated frames).
LONG_LADDER_RUNGS = 8


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:                                          # pragma: no cover
        return ""


def read_image(archive: Path, rel_path: str) -> tuple[np.ndarray, fits.Header]:
    """Read one archive frame's pixels + header (fpack .fz or plain FITS).

    fpack files carry data in HDU 1 (CompImageHDU); plain files in HDU 0.
    Returns the data as stored (uint16 for the int16+BZERO raw frames).
    """
    with fits.open(archive / rel_path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                return np.asarray(hdu.data), hdu.header
    raise ValueError(f"no 2-D image HDU in {rel_path}")


#: SQL predicate for "this is a science exposure".
#:
#: The obvious ``imagetyp LIKE 'Light%'`` is WRONG for this archive and the
#: ceiling stage already knew it: the 2026 camera writes no IMAGETYP card at
#: all (the same unwritten-header condition that leaves READOUTM blank), so a
#: strict test silently drops every blank-2026 frame.  Measured on the
#: current manifest, canonical rawimage frames with a null/empty IMAGETYP
#: number 2,212 and ALL of them are blank-READOUTM 2026 frames — so allowing
#: the blank widens the net for exactly that camera and for nothing else.
#: The ceiling stage's sample was built with this rule; noise and linearity
#: now use the same one, so all three describe the same frame population.
SCIENCE_IMAGETYP = ("(f.imagetyp LIKE 'Light%' OR f.imagetyp IS NULL "
                    "OR trim(f.imagetyp) = '')")


def mode_where(mode: str) -> tuple[str, tuple]:
    """SQL fragment selecting one ceiling mode group's science frames."""
    if mode == "(blank 2026)":
        return ("(f.readoutm IS NULL OR trim(f.readoutm) = '')", ())
    return ("f.readoutm = ?", (mode,))


def open_db(path: Path) -> sqlite3.Connection:
    # The manifest is shared with sibling pipeline stages that may hold
    # their own connections: never change the journal mode (that needs an
    # exclusive lock), just wait politely when a writer is ahead of us.
    # 300 s: sibling stages (the S1 batch solver runs ten workers) can hold
    # the write lock for minutes at a time; waiting is always cheaper than
    # failing a batch that has already read pixels off the archive.
    con = sqlite3.connect(path, timeout=300.0)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


def ensure_tables(con: sqlite3.Connection) -> None:
    """Create every S2 table (new tables only — S0/S0b are never touched)."""
    con.executescript("""
    CREATE TABLE IF NOT EXISTS s2_ceiling_frames (
        obs_rowid INTEGER PRIMARY KEY, mode TEXT, night TEXT, exptime REAL,
        max_adu INTEGER, n_at_max INTEGER, p999_adu REAL);
    CREATE TABLE IF NOT EXISTS s2_ceiling_hist (
        mode TEXT, adu INTEGER, count INTEGER, PRIMARY KEY (mode, adu));
    CREATE TABLE IF NOT EXISTS s2_ceiling_modes (
        mode TEXT PRIMARY KEY, n_frames INTEGER, n_pixels INTEGER,
        hard_max_adu INTEGER, clip_adu INTEGER, spike_count INTEGER,
        tail_level REAL, ratio REAL, veto_adu INTEGER, bits INTEGER,
        adc_full_scale INTEGER, unused_codes INTEGER);
    CREATE TABLE IF NOT EXISTS s2_ptc_points (
        pair_id TEXT, mode TEXT, kind TEXT, exptime REAL,
        level REAL, var REAL, n_pix INTEGER);
    CREATE TABLE IF NOT EXISTS s2_ptc_pairs (
        pair_id TEXT PRIMARY KEY, mode TEXT, kind TEXT, exptime REAL,
        scene TEXT, path_a TEXT, path_b TEXT, n_points INTEGER);
    CREATE TABLE IF NOT EXISTS s2_ptc_fits (
        mode TEXT, kind TEXT, gain_e_per_adu REAL, gain_err REAL,
        read_noise_adu REAL, read_noise_adu_err REAL, read_noise_e REAL,
        slope REAL, intercept REAL, n_points INTEGER,
        PRIMARY KEY (mode, kind));
    CREATE TABLE IF NOT EXISTS s2_ampglow (
        obs_rowid INTEGER PRIMARY KEY, mode TEXT, exptime REAL,
        center_med REAL, edge_med REAL, edge_excess REAL,
        hottest_corner_med REAL, hottest_corner_excess REAL);
    CREATE TABLE IF NOT EXISTS s2_recon_eras (
        era_id INTEGER PRIMARY KEY, mode TEXT, n_links INTEGER,
        n_pairs_used INTEGER, exptime_med REAL, pedestal REAL,
        flat_median REAL, flat_mad_sigma REAL, dark_median REAL,
        dark_mad_sigma REAL, fit_fraction REAL, rms_median REAL,
        truth_master TEXT, truth_offset REAL, truth_resid_rms REAL,
        truth_resid_mad REAL, truth_n_pix INTEGER, npz_path TEXT);
    CREATE TABLE IF NOT EXISTS s2_linearity_ladders (
        ladder_id TEXT PRIMARY KEY, mode TEXT, night TEXT, target_key TEXT,
        n_rungs INTEGER, n_frames INTEGER, rate_adu_per_s REAL,
        max_abs_resid_pct REAL);
    CREATE TABLE IF NOT EXISTS s2_linearity_rungs (
        ladder_id TEXT, exptime REAL, n_frames INTEGER, flux_med REAL,
        peak_med REAL, resid_pct REAL, PRIMARY KEY (ladder_id, exptime));
    CREATE TABLE IF NOT EXISTS detector_params (
        era_group TEXT, quantity TEXT, value REAL, uncertainty REAL,
        method TEXT, provenance TEXT, PRIMARY KEY (era_group, quantity));
    CREATE TABLE IF NOT EXISTS s2_build_meta (key TEXT PRIMARY KEY, value TEXT);
    -- Per-(mode, egain) near-ceiling frame-max statistics.  The adversarial
    -- review showed the High Gain "mound" pools two egain epochs whose clip
    -- levels are cleanly separated (1.054: 3,526-3,584; 1.057: 3,427-3,546),
    -- so the pooled ceiling must be readable per epoch.
    CREATE TABLE IF NOT EXISTS s2_ceiling_egain (
        mode TEXT, egain REAL, n_frames INTEGER, min_max_adu INTEGER,
        median_max_adu REAL, max_max_adu INTEGER, PRIMARY KEY (mode, egain));
    -- The empirical noise model (CV-P15-noise-model).  s2_noise_pairs and
    -- s2_noise_points are the MEASUREMENTS (one row per level bin of one
    -- same-scene frame pair); s2_noise_curve is the distilled per-mode
    -- lookup table the photometry error model actually reads.
    CREATE TABLE IF NOT EXISTS s2_noise_pairs (
        pair_id TEXT PRIMARY KEY, mode TEXT, night TEXT, target_key TEXT,
        filter TEXT, egain REAL, exptime REAL, era_id INTEGER,
        path_a TEXT, path_b TEXT, gap_s REAL, n_points INTEGER);
    CREATE TABLE IF NOT EXISTS s2_noise_points (
        pair_id TEXT, mode TEXT, level REAL, var REAL, n_pix INTEGER);
    CREATE TABLE IF NOT EXISTS s2_noise_curve (
        mode TEXT, bin_index INTEGER, level_adu REAL, var_adu2 REAL,
        var_mad_adu2 REAL, sigma_adu REAL, n_points INTEGER,
        n_pairs INTEGER, n_pix REAL, PRIMARY KEY (mode, bin_index));
    """)
    # Columns added after the first field campaign (schema migration for an
    # existing s2_ceiling_frames): the argmax position, which arbitrates
    # hot-pixel clusters vs true ceilings.  ALTER fails harmlessly when the
    # column already exists.
    for col in ("max_y", "max_x"):
        try:
            con.execute(f"ALTER TABLE s2_ceiling_frames ADD COLUMN {col} INTEGER")
        except sqlite3.OperationalError:
            pass
    # Frame-max-cluster evidence per mode (kept even when the cluster is
    # REJECTED as a hot pixel — the report cites the rejection from here).
    for col, typ in (("cluster_adu", "REAL"), ("cluster_diversity", "REAL"),
                     ("cluster_n_pos", "INTEGER")):
        try:
            con.execute(f"ALTER TABLE s2_ceiling_modes ADD COLUMN {col} {typ}")
        except sqlite3.OperationalError:
            pass
    # The F-D degeneracy diagnostic per era (adversarial-review addition):
    # Pearson correlation of the per-pixel F and D estimates, computed from
    # the stored npz products.  Strongly negative = poor level diversity =
    # per-pixel F values are noise and only the median carries meaning.
    try:
        con.execute("ALTER TABLE s2_recon_eras ADD COLUMN fd_corr REAL")
    except sqlite3.OperationalError:
        pass
    # The reduced-vs-raw CROP, measured per era by patch alignment.  The
    # report used to state "the 2026 eras crop by 13-18 pixels" as prose;
    # prose is not evidence, and the geometry repair moved which eras those
    # are.  Recording the measured offset makes the sentence a query.
    for col in ("crop_dy", "crop_dx", "n_pairs_cropped"):
        try:
            con.execute(f"ALTER TABLE s2_recon_eras ADD COLUMN {col} INTEGER")
        except sqlite3.OperationalError:
            pass
    con.commit()


def write_meta(con: sqlite3.Connection, manifest: Path) -> None:
    for k, v in [("built_utc", utcnow()), ("code_version", S2_CODE_VERSION),
                 ("git_commit", git_commit()),
                 ("manifest_path", str(manifest))]:
        con.execute("INSERT OR REPLACE INTO s2_build_meta VALUES (?, ?)", (k, v))
    con.commit()


# ---------------------------------------------------------------------------
# Subcommand: ceiling
# ---------------------------------------------------------------------------
def cmd_ceiling(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Sample science frames per mode, accumulate pixel histograms."""
    todo: list[tuple] = []
    for mode, target in CEILING_MODES.items():
        cond, params = mode_where(mode)
        rows = con.execute(f"""
            SELECT f.obs_rowid, f.path, f.night, f.exptime FROM frames f
            WHERE {cond} AND f.is_canonical = 1 AND f.tree = 'rawimage'
              AND {SCIENCE_IMAGETYP}
            ORDER BY f.obs_rowid""", params).fetchall()
        if not rows:
            continue
        # Deterministic spread across the era: every k-th frame.
        step = max(1, len(rows) // target)
        sample = rows[::step][:target]
        done = {r[0] for r in con.execute(
            "SELECT obs_rowid FROM s2_ceiling_frames WHERE mode = ?", (mode,))}
        todo += [(mode,) + r for r in sample if r[0] not in done]

    if not todo:
        print("[S2:ceiling] nothing left to do.")
        return 0
    todo = todo[:batch]
    print(f"[S2:ceiling] processing {len(todo)} frames ...")
    hist_acc: dict[str, np.ndarray] = {}
    frame_rows, skipped = [], 0
    for mode, rowid, path, night, exptime in todo:
        try:
            data, _ = read_image(archive, path)
        except Exception as e:
            print(f"[S2:ceiling]   SKIP {path}: {e}")
            skipped += 1
            # Record a sentinel so the frame is not retried forever.
            frame_rows.append((rowid, mode, night, exptime, -1, 0, -1.0,
                               None, None))
            continue
        flat = np.asarray(data).ravel()
        if flat.dtype.kind == "f":                    # a rare float frame
            flat = np.clip(flat, 0, 65535).astype(np.uint16)
        h = np.bincount(flat, minlength=65536)
        hist_acc[mode] = ceil.merge_hist(hist_acc.get(mode,
                                                      np.zeros(1, np.int64)), h)
        st = ceil.frame_top_stats(data, None)
        frame_rows.append((rowid, mode, night, exptime, st["max_adu"],
                           st["n_at_max"], st["p999_adu"],
                           st["max_y"], st["max_x"]))
    # One transaction per batch: frames + histogram increments together.
    con.executemany("INSERT OR REPLACE INTO s2_ceiling_frames "
                    "(obs_rowid, mode, night, exptime, max_adu, n_at_max, "
                    "p999_adu, max_y, max_x) VALUES (?,?,?,?,?,?,?,?,?)",
                    frame_rows)
    for mode, h in hist_acc.items():
        nz = np.flatnonzero(h)
        con.executemany("""
            INSERT INTO s2_ceiling_hist (mode, adu, count) VALUES (?,?,?)
            ON CONFLICT(mode, adu) DO UPDATE SET count = count + excluded.count
            """, [(mode, int(a), int(h[a])) for a in nz])
    con.commit()
    print(f"[S2:ceiling] batch done ({len(frame_rows)} frames, "
          f"{skipped} unreadable). Re-run until 'nothing left to do'.")
    return 0


def cmd_ceilpos(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Backfill argmax positions for frame-max-cluster candidate frames.

    Only frames whose maximum falls inside a mode's frame-max cluster
    window need a position (they are the cluster's diversity evidence);
    frames sampled before the max_y/max_x columns existed get re-read here.
    """
    n_done = 0
    for mode in CEILING_MODES:
        maxes = [r[0] for r in con.execute(
            "SELECT max_adu FROM s2_ceiling_frames WHERE mode=? AND max_adu>0",
            (mode,))]
        cl = ceil.frame_max_cluster(maxes)
        if cl is None:
            continue
        lo = cl["clip_adu"] * (1 - ceil.CLUSTER_REL_WINDOW)
        hi = cl["clip_adu"] * (1 + ceil.CLUSTER_REL_WINDOW)
        rows = con.execute("""
            SELECT c.obs_rowid, f.path FROM s2_ceiling_frames c
            JOIN frames f ON f.obs_rowid = c.obs_rowid
            WHERE c.mode = ? AND c.max_adu BETWEEN ? AND ?
              AND c.max_y IS NULL""", (mode, lo, hi)).fetchall()
        for rowid, path in rows:
            if n_done >= batch:
                print("[S2:ceilpos] batch cap reached; re-run to continue.")
                return 0
            try:
                data, _ = read_image(archive, path)
            except Exception as e:
                print(f"[S2:ceilpos]   SKIP {path}: {e}")
                continue
            st = ceil.frame_top_stats(data, None)
            con.execute("UPDATE s2_ceiling_frames SET max_y=?, max_x=? "
                        "WHERE obs_rowid=?", (st["max_y"], st["max_x"], rowid))
            con.commit()
            n_done += 1
    print(f"[S2:ceilpos] done ({n_done} positions this run; 0 = complete).")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: ptc
# ---------------------------------------------------------------------------
def _ptc_groups(con: sqlite3.Connection) -> list[dict]:
    """Same-(scene, mode, exptime) frame series on the PTC night.

    Scene identity: the filename family for the cmos_tests/latency darks
    and lights (everything before the trailing counter), the target_key for
    rawimage science.  Consecutive-in-JD frames within a series are paired.
    """
    rows = con.execute("""
        SELECT f.obs_rowid, f.path, f.basename, f.tree, f.readoutm,
               f.exptime, f.imagetyp, f.jd, f.target_key
        FROM frames f
        WHERE f.night = ? AND f.is_canonical = 1
          AND f.naxis1 = 4096 AND f.readoutm LIKE 'High Gain%'
        ORDER BY f.jd""", (PTC_NIGHT,)).fetchall()
    from collections import defaultdict
    groups: dict[tuple, list] = defaultdict(list)
    for (rowid, path, base, tree, mode, expt, ityp, jd, tkey) in rows:
        kind = "dark" if (ityp or "").startswith("Dark") else "light"
        if tree in ("cmos_tests", "latency"):
            # Albireo_dark.64.01.fts.fz -> scene 'Albireo_dark.64' etc.
            scene = ".".join(base.split(".")[:2])
        elif tree == "rawimage" and kind == "light" and tkey:
            scene = tkey
        else:
            continue
        groups[(mode, ceil.mode_group(mode), kind,
                exptime_bin(expt), scene)].append((jd, rowid, path))
    out = []
    for (mode, _label, kind, ebin, scene), members in sorted(groups.items()):
        if len(members) < 2 or ebin is None:
            continue
        members.sort()
        pairs = [(members[i], members[i + 1])
                 for i in range(0, len(members) - 1, 2)]
        out.append({"mode": mode, "kind": kind, "exptime": ebin,
                    "scene": scene, "pairs": pairs[:PTC_MAX_PAIRS_PER_GROUP]})
    return out


def cmd_ptc(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Difference-pair photon transfer on the 2023-06-07 series."""
    done = {r[0] for r in con.execute("SELECT pair_id FROM s2_ptc_pairs")}
    ceilings = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes WHERE clip_adu IS NOT NULL"))
    n_done = 0
    for g in _ptc_groups(con):
        clip = ceilings.get(g["mode"])
        level_max = (ptc.PTC_LEVEL_CEILING_FRACTION * clip) if clip else None
        for (jd_a, id_a, path_a), (jd_b, id_b, path_b) in g["pairs"]:
            pair_id = f"{id_a}-{id_b}"
            if pair_id in done:
                continue
            if n_done >= batch:
                print("[S2:ptc] batch cap reached; re-run to continue.")
                return 0
            try:
                a, _ = read_image(archive, path_a)
                b, _ = read_image(archive, path_b)
            except Exception as e:
                print(f"[S2:ptc]   SKIP pair {pair_id}: {e}")
                con.execute("INSERT OR REPLACE INTO s2_ptc_pairs VALUES "
                            "(?,?,?,?,?,?,?,0)",
                            (pair_id, g["mode"], g["kind"], g["exptime"],
                             g["scene"], path_a, path_b))
                con.commit()
                continue
            points = ptc.pair_ptc_points(a.astype(np.float64),
                                         b.astype(np.float64),
                                         level_max=level_max)
            # REPLACE this pair's points, never append to them.  s2_ptc_pairs
            # is keyed on pair_id and so is naturally idempotent, but the
            # points table has no key — so a pair measured twice (two
            # campaign processes overlapping, or a batch re-run after a kill
            # that landed between the two writes) silently doubled its
            # points, and a doubled point is a doubled WEIGHT in the fit.
            # Deleting first makes re-measuring a pair a no-op, which is
            # what "resumable" has to mean.
            con.execute("DELETE FROM s2_ptc_points WHERE pair_id = ?",
                        (pair_id,))
            con.executemany(
                "INSERT INTO s2_ptc_points VALUES (?,?,?,?,?,?,?)",
                [(pair_id, g["mode"], g["kind"], g["exptime"],
                  p["level"], p["var"], p["n_pix"]) for p in points])
            con.execute("INSERT OR REPLACE INTO s2_ptc_pairs VALUES "
                        "(?,?,?,?,?,?,?,?)",
                        (pair_id, g["mode"], g["kind"], g["exptime"],
                         g["scene"], path_a, path_b, len(points)))
            con.commit()
            n_done += 1
            # Amp-glow check rides along on the longest darks (>= 100 s).
            if g["kind"] == "dark" and g["exptime"] >= 100:
                for rid, pth, img in ((id_a, path_a, a), (id_b, path_b, b)):
                    if con.execute("SELECT 1 FROM s2_ampglow WHERE obs_rowid=?",
                                   (rid,)).fetchone():
                        continue
                    m = ptc.amp_glow_metric(img)
                    con.execute("INSERT OR REPLACE INTO s2_ampglow VALUES "
                                "(?,?,?,?,?,?,?,?)",
                                (rid, g["mode"], g["exptime"], m["center_med"],
                                 m["edge_med"], m["edge_excess"],
                                 m["hottest_corner_med"],
                                 m["hottest_corner_excess"]))
                con.commit()
    print(f"[S2:ptc] done ({n_done} new pairs this run; "
          "0 new pairs = nothing left to do).")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: noise  (the empirical per-mode counts-vs-variance model)
# ---------------------------------------------------------------------------
def _noise_pairs_for_mode(con: sqlite3.Connection, mode: str) -> list[dict]:
    """Same-scene consecutive science pairs for one readout mode.

    A "scene" is one (night, target, filter, camera gain, exposure bin,
    geometry) group — everything that must match before two frames can be
    differenced meaningfully.  Within a scene the frames are ordered in
    time and adjacent ones are paired, provided their gap stays under
    :data:`NOISE_MAX_GAP_DAYS`.

    Selection is FAIR, not greedy: scenes are visited richest-first and
    each contributes at most :data:`NOISE_MAX_PAIRS_PER_SCENE` pairs, in
    round-robin passes, until the mode has
    :data:`NOISE_MAX_PAIRS_PER_MODE`.  A curve pooled from many scenes is a
    statement about the DETECTOR; a curve pooled from one night's 16 pairs
    is a statement about that night.
    """
    cond, params = mode_where(mode)
    rows = con.execute(f"""
        SELECT f.obs_rowid, f.path, f.night, f.target_key,
               coalesce(f.filter, ''), coalesce(f.egain, -1), f.exptime,
               f.jd, f.naxis1, f.naxis2, f.era_id
        FROM frames f
        WHERE {cond} AND f.is_canonical = 1 AND f.tree = 'rawimage'
          AND {SCIENCE_IMAGETYP} AND f.target_key IS NOT NULL
          AND f.exptime > 0 AND f.jd IS NOT NULL
        ORDER BY f.jd""", params).fetchall()
    from collections import defaultdict
    scenes: dict[tuple, list] = defaultdict(list)
    for (rowid, path, night, tkey, filt, eg, expt, jd, n1, n2, era) in rows:
        ebin = exptime_bin(expt)
        if ebin is None or ebin <= 0:
            continue
        scenes[(night, tkey, filt, eg, ebin, n1, n2, era)].append(
            (jd, rowid, path))
    # Candidate pairs per scene, in time order.
    per_scene: list[tuple[tuple, list]] = []
    for key, members in scenes.items():
        if len(members) < NOISE_MIN_SCENE_FRAMES:
            continue
        members.sort()
        pairs = [(members[i], members[i + 1])
                 for i in range(0, len(members) - 1, 2)
                 if members[i + 1][0] - members[i][0] <= NOISE_MAX_GAP_DAYS]
        if pairs:
            per_scene.append((key, pairs))
    # Richest scenes first; ties broken by the scene key so the choice is
    # reproducible run to run.
    per_scene.sort(key=lambda kv: (-len(kv[1]), str(kv[0])))
    out: list[dict] = []
    for slot in range(NOISE_MAX_PAIRS_PER_SCENE):
        for key, pairs in per_scene:
            if len(out) >= NOISE_MAX_PAIRS_PER_MODE:
                return out
            if slot >= len(pairs):
                continue
            (jd_a, id_a, path_a), (jd_b, id_b, path_b) = pairs[slot]
            night, tkey, filt, eg, ebin, _n1, _n2, era = key
            out.append({
                "pair_id": f"n{id_a}-{id_b}", "mode": mode, "night": night,
                "target_key": tkey, "filter": filt, "egain": eg,
                "exptime": ebin, "era_id": era, "path_a": path_a,
                "path_b": path_b, "gap_s": (jd_b - jd_a) * 86400.0})
    return out


def cmd_noise(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Measure counts-vs-variance for every readout mode.

    This is the same difference-pair arithmetic the PTC uses, applied to
    ordinary science frames across the WHOLE archive instead of one
    dedicated night — which is what makes it available for every mode.  It
    buys breadth at the cost of purity: science pairs contain stars that
    move a fraction of a pixel between exposures, so the bright end of each
    curve is an UPPER bound on detector noise.  That limitation is recorded
    with the curve (``curve_shape_index`` reads >1 when it bites) rather
    than hidden.
    """
    done = {r[0] for r in con.execute("SELECT pair_id FROM s2_noise_pairs")}
    # Per-mode ceilings gate the level axis: pixels near the clip have an
    # artificially collapsed variance and would drag every fit down.
    ceilings = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes WHERE clip_adu IS NOT NULL"))
    n_done = 0
    for mode in CEILING_MODES:
        clip = ceilings.get(mode)
        level_max = (ptc.PTC_LEVEL_CEILING_FRACTION * clip) if clip else None
        for pr in _noise_pairs_for_mode(con, mode):
            if pr["pair_id"] in done:
                continue
            if n_done >= batch:
                print("[S2:noise] batch cap reached; re-run to continue.")
                return 0
            try:
                a, _ = read_image(archive, pr["path_a"])
                b, _ = read_image(archive, pr["path_b"])
            except Exception as e:
                print(f"[S2:noise]   SKIP {pr['pair_id']}: {e}")
                # Sentinel row: an unreadable pair is recorded once, never
                # retried, and never counted as a measurement (n_points 0).
                con.execute("INSERT OR REPLACE INTO s2_noise_pairs VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?,0)",
                            (pr["pair_id"], mode, pr["night"],
                             pr["target_key"], pr["filter"], pr["egain"],
                             pr["exptime"], pr["era_id"], pr["path_a"],
                             pr["path_b"], pr["gap_s"]))
                con.commit()
                continue
            if a.shape != b.shape:
                print(f"[S2:noise]   SKIP {pr['pair_id']}: shape mismatch")
                continue
            points = ptc.pair_ptc_points(a.astype(np.float64),
                                         b.astype(np.float64),
                                         level_max=level_max)
            # Same idempotence rule as the PTC points (see cmd_ptc): a pair
            # re-measured must REPLACE its rows, or a duplicated point
            # becomes a duplicated weight in the mode's curve.
            con.execute("DELETE FROM s2_noise_points WHERE pair_id = ?",
                        (pr["pair_id"],))
            con.executemany(
                "INSERT INTO s2_noise_points VALUES (?,?,?,?,?)",
                [(pr["pair_id"], mode, p["level"], p["var"], p["n_pix"])
                 for p in points])
            con.execute("INSERT OR REPLACE INTO s2_noise_pairs VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?)",
                        (pr["pair_id"], mode, pr["night"], pr["target_key"],
                         pr["filter"], pr["egain"], pr["exptime"],
                         pr["era_id"], pr["path_a"], pr["path_b"],
                         pr["gap_s"], len(points)))
            con.commit()
            n_done += 1
    # Falling out of the loop means every mode's pair list is exhausted —
    # the batch cap returns early above, so this line is the real "done".
    print(f"[S2:noise] nothing left to do. ({n_done} new pairs this run)")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: reconstruct
# ---------------------------------------------------------------------------
def _era_links(con: sqlite3.Connection, era_id: int) -> list[tuple]:
    """Usable raw<->reduced pairs for one era (unambiguous methods only)."""
    return con.execute("""
        SELECT l.raw_path, l.reduced_path, f.exptime, f.naxis1, f.naxis2
        FROM raw_reduced_links l JOIN frames f ON f.obs_rowid = l.raw_rowid
        WHERE f.era_id = ? AND l.raw_rowid IS NOT NULL
          AND l.match_method IN ('same_basename_jd', 'stem_jd',
                                 'stem_jd_drift', 'target_jd')
        ORDER BY l.raw_rowid""", (era_id,)).fetchall()


def recon_eras(con: sqlite3.Connection) -> list[tuple[int, str, int]]:
    """Eras that qualify for the reconstruction experiment."""
    return [tuple(r) for r in con.execute("""
        SELECT f.era_id, max(e.readoutm), count(*) AS n
        FROM raw_reduced_links l
        JOIN frames f ON f.obs_rowid = l.raw_rowid
        JOIN eras e ON e.era_id = f.era_id
        WHERE l.match_method IN ('same_basename_jd', 'stem_jd',
                                 'stem_jd_drift', 'target_jd')
        GROUP BY f.era_id HAVING n >= ? ORDER BY n DESC""",
        (RECON_MIN_LINKS,))]


def cmd_reconstruct(con: sqlite3.Connection, archive: Path,
                    only_era: int | None) -> int:
    """Fit D and F per pixel for every qualifying era (one era per call
    unless --era is given; the s2_recon_eras table records completion)."""
    PRODUCTS_DIR.joinpath("recon").mkdir(parents=True, exist_ok=True)
    ceilings = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes WHERE clip_adu IS NOT NULL"))
    done = {r[0] for r in con.execute("SELECT era_id FROM s2_recon_eras")}
    for era_id, mode, n_links in recon_eras(con):
        if only_era is not None and era_id != only_era:
            continue
        if era_id in done:
            continue
        print(f"[S2:recon] era {era_id} ({mode!r}, {n_links} links) ...")
        links = _era_links(con, era_id)
        # Dominant geometry: the (naxis1, naxis2) most links share.
        from collections import Counter
        geom, _ = Counter((l[3], l[4]) for l in links).most_common(1)[0]
        links = [l for l in links if (l[3], l[4]) == geom]
        step = max(1, len(links) // RECON_MAX_PAIRS)
        chosen = links[::step][:RECON_MAX_PAIRS]

        regions = None
        red_stack, raw_stack, exptimes, pedestals = [], [], [], []
        crops: list[tuple[int, int]] = []       # measured reduced-vs-raw crop
        for raw_path, red_path, expt, _n1, _n2 in chosen:
            try:
                raw_img, _ = read_image(archive, raw_path)
                red_img, red_hdr = read_image(archive, red_path)
            except Exception as e:
                print(f"[S2:recon]   SKIP pair {raw_path}: {e}")
                continue
            if raw_img.shape != red_img.shape:
                # The 2026 pipeline crops its reduced output by a few
                # rows/columns; measure the crop offset and align the raw
                # frame onto the reduced grid.  Anything that is not a
                # small crop (resampled/stacked product) stays skipped.
                off = rec.find_crop_offset(raw_img, red_img)
                if off is None:
                    continue
                # How many rows/columns the reduction dropped, in total —
                # the number the report used to assert in prose.
                crops.append((raw_img.shape[0] - red_img.shape[0],
                              raw_img.shape[1] - red_img.shape[1]))
                raw_img = raw_img[off["dy"]:off["dy"] + red_img.shape[0],
                                  off["dx"]:off["dx"] + red_img.shape[1]]
            if regions is None:
                regions = rec.sample_regions(*red_img.shape)
            ped = float(red_hdr.get("PEDESTAL", 0) or 0)
            pick = lambda im: np.concatenate(
                [np.asarray(im[ys, xs], dtype=np.float64).ravel()
                 for _nm, ys, xs in regions])
            red_stack.append(pick(red_img) - ped)
            raw_stack.append(pick(raw_img))
            exptimes.append(expt)
            pedestals.append(ped)
        if len(red_stack) < rec.RECON_MIN_PAIRS:
            print(f"[S2:recon]   era {era_id}: only {len(red_stack)} readable "
                  "pairs — recorded as unfittable.")
            con.execute("INSERT OR REPLACE INTO s2_recon_eras (era_id, mode, "
                        "n_links, n_pairs_used) VALUES (?,?,?,?)",
                        (era_id, mode, n_links, len(red_stack)))
            con.commit()
            continue

        fit = rec.fit_pixel_lines(np.stack(red_stack), np.stack(raw_stack),
                                  sat_adu=ceilings.get(mode))
        summ = rec.summarize_reconstruction(fit["F"], fit["D"], fit["rms"])

        # Ground truth for era 47: the archived master bias (+ scaled dark).
        truth_name, tr = None, {"offset": None, "resid_rms": None,
                                "resid_mad_sigma": None, "n_pix": None}
        exp_med = float(np.median([e for e in exptimes if e is not None]
                                  or [0.0]))
        truth_stamp = None
        masters = con.execute("""
            SELECT path, kind, exptime FROM calib_frames
            WHERE era_id = ? AND is_master = 1 AND kind IN ('bias','dark')
            ORDER BY kind""", (era_id,)).fetchall()
        bias_row = next((m for m in masters if m[1] == "bias"), None)
        # The master dark must be the one CLOSEST IN EXPOSURE to the pairs,
        # not simply the first row the query returns.
        #
        # The truth model is  truth = bias + (t_pairs / t_dark) * (dark - bias),
        # so the scale factor is a ratio of exposure times and a badly chosen
        # dark multiplies its own noise by that ratio.  Era 76 caught this:
        # its master set opens with a 0.01 s dark, which against a 31 s
        # median exposure is a 3,100x extrapolation — every pixel where the
        # dark differs from the bias by one ADU acquired 3,100 ADU of
        # invented signal, and the era graded at 2,638 ADU RMS against a
        # MAD-sigma of 5.7 (the giveaway: the bulk of pixels agreed fine and
        # a handful of extrapolated outliers owned the RMS).
        #
        # Closest in LOG exposure, because the damage is multiplicative; and
        # beyond RECON_TRUTH_MAX_SCALE either way the dark term is dropped
        # rather than extrapolated, leaving an honest bias-only truth.
        darks = [m for m in masters if m[1] == "dark" and m[2] and m[2] > 0]
        dark_row = None
        if darks and exp_med > 0:
            dark_row = min(darks, key=lambda m: abs(np.log(m[2] / exp_med)))
            scale = exp_med / dark_row[2]
            if not (1.0 / RECON_TRUTH_MAX_SCALE <= scale
                    <= RECON_TRUTH_MAX_SCALE):
                print(f"[S2:recon]   era {era_id}: nearest master dark is "
                      f"{dark_row[2]:g}s against {exp_med:g}s pairs "
                      f"({scale:.3g}x) — beyond the "
                      f"{RECON_TRUTH_MAX_SCALE:g}x extrapolation limit, so "
                      "the truth is bias-only.")
                dark_row = None
        if bias_row:
            try:
                bias_img, _ = read_image(archive, bias_row[0])
                truth = np.asarray(bias_img, dtype=np.float64)
                truth_name = bias_row[0]
                if dark_row and dark_row[2]:
                    dark_img, _ = read_image(archive, dark_row[0])
                    # Master dark includes bias; add the dark-current term
                    # scaled to the pairs' median exposure time.
                    truth = truth + (exp_med / float(dark_row[2])) * (
                        np.asarray(dark_img, dtype=np.float64) - truth)
                    truth_name += f" + {dark_row[0]} x {exp_med:g}s"
                # Geometry gate: the master must share the era's (ny, nx)
                # frame shape, or the region slices would not correspond.
                truth_pix = np.concatenate(
                    [truth[ys, xs].ravel() for _nm, ys, xs in regions]) \
                    if truth.shape == (geom[1], geom[0]) else None
                if truth_pix is not None:
                    tr = rec.residual_vs_truth(fit["D"], truth_pix)
                    truth_stamp = truth_pix
            except Exception as e:
                print(f"[S2:recon]   era {era_id}: truth unreadable: {e}")

        npz_path = PRODUCTS_DIR / "recon" / f"era{era_id}.npz"
        # Atomic write: temp name, then replace.  The temp name must end in
        # '.npz' or numpy silently appends the extension and the rename
        # target never exists.
        # The temp name carries THIS process's pid.  A shared "eraNN.tmp.npz"
        # is not actually atomic when two campaign processes reconstruct the
        # same era: the first one's replace() consumes the shared temp file
        # and the second dies with FileNotFoundError on a path it just
        # wrote.  (Observed, not hypothesised.)  The name must still end in
        # '.npz' or numpy appends the extension and the rename target never
        # exists.
        tmp = npz_path.with_suffix(f".tmp{os.getpid()}.npz")
        np.savez_compressed(
            tmp, F=fit["F"], D=fit["D"], n_used=fit["n_used"], rms=fit["rms"],
            regions=np.array([(nm, ys.start, ys.stop, xs.start, xs.stop)
                              for nm, ys, xs in regions], dtype=object),
            pedestal=np.array(pedestals), exptime=np.array(
                [e if e is not None else np.nan for e in exptimes]),
            truth=(truth_stamp if truth_stamp is not None else np.array([])))
        tmp.replace(npz_path)

        # Columns are named, never positional: this table has been ALTERed
        # twice (fd_corr, then the crop columns), and a bare
        # `INSERT ... VALUES (?,...)` breaks the moment the column count
        # moves — silently in review, loudly at 3 a.m. mid-campaign.
        con.execute("""INSERT OR REPLACE INTO s2_recon_eras
            (era_id, mode, n_links, n_pairs_used, exptime_med, pedestal,
             flat_median, flat_mad_sigma, dark_median, dark_mad_sigma,
             fit_fraction, rms_median, truth_master, truth_offset,
             truth_resid_rms, truth_resid_mad, truth_n_pix, npz_path,
             crop_dy, crop_dx, n_pairs_cropped)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (era_id, mode, n_links, len(red_stack), exp_med,
             float(np.median(pedestals)), summ["flat_median"],
             summ["flat_mad_sigma"], summ["dark_median"],
             summ["dark_mad_sigma"], summ["fit_fraction"], summ["rms_median"],
             truth_name, tr["offset"], tr["resid_rms"], tr["resid_mad_sigma"],
             tr["n_pix"], str(npz_path.relative_to(REPO_ROOT)),
             int(np.median([c[0] for c in crops])) if crops else None,
             int(np.median([c[1] for c in crops])) if crops else None,
             len(crops) or None))
        con.commit()
        print(f"[S2:recon]   era {era_id}: F~{summ['flat_median']:.4f} "
              f"D~{summ['dark_median']:.1f} ADU, rms {summ['rms_median']:.2f}")
        if only_era is None:
            # One era per invocation keeps every call under the time cap.
            print("[S2:recon] one era per call — re-run for the next.")
            return 0
    print("[S2:recon] nothing left to do.")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: linearity
# ---------------------------------------------------------------------------
def cmd_linearity(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Fit every archival exposure ladder the manifest can surface."""
    # A ladder rung is only comparable to its neighbours under the SAME
    # filter and the SAME camera gain: the archive mixes filters mid-visit
    # (the Vega 0.1 s rung is HaGrism among OGGrism rungs) and the 2026
    # Fast camera toggles EGAIN mid-sequence (HDR pairs, 16x apart) — both
    # would fake enormous "non-linearity" if pooled.
    cands = con.execute(f"""
        SELECT f.night, f.target_key, f.readoutm,
               coalesce(f.filter, '') AS filt,
               coalesce(f.egain, -1) AS eg,
               count(DISTINCT round(coalesce(f.exptime, -1), 4)) AS nex,
               count(*) AS n
        FROM frames f
        WHERE f.is_canonical = 1 AND f.tree = 'rawimage'
          AND {SCIENCE_IMAGETYP} AND f.target_key IS NOT NULL
        GROUP BY f.night, f.target_key, f.readoutm, filt, eg
        HAVING nex >= 3 AND n >= 6
        ORDER BY nex DESC, n DESC""").fetchall()
    # Ordering, in two layers.
    #
    # QUALITY (within a mode): the Vega BeStar ladder is the roadmap's
    # named case and goes first; after it, prefer ladders with the most
    # FRAMES PER RUNG (the rung median beats transients/clouds; the field
    # campaign showed dedicated 14-rung single-frame sequences drowning in
    # sky variation), then the most rungs.
    #
    # FAIRNESS (across modes): take every mode's best candidate before any
    # mode's second.  The earlier pass ranked globally, which handed all
    # twelve slots to the two richest modes and left Low Gain, 5 MHz and
    # blank-2026 with no ladder at all — exactly the modes CV-P15 needs
    # covered.  Round-robin guarantees each mode is TRIED; whether its
    # frames support a fit is then a finding, not a scheduling accident.
    ordered = lin.fair_ladder_order(
        cands,
        mode_of=lambda r: ceil.mode_group(r[2]),
        quality=lambda r: (0 if "vega" in (r[1] or "") else 1,
                           -r[6] / r[5], -r[5]))
    done = {r[0] for r in con.execute(
        "SELECT ladder_id FROM s2_linearity_ladders")}
    n_run = 0
    for night, tkey, mode, filt, eg, _nex, _n in ordered[:MAX_LADDERS * 6]:
        if n_run >= min(batch, MAX_LADDERS):
            break
        # The mode LABEL is the canonical group name ("(blank 2026)" for the
        # current camera's unwritten READOUTM cards), because every other S2
        # table keys on that label — a ladder stored under a raw NULL/'' mode
        # could never join s2_ceiling_modes for its saturation veto.  The raw
        # value is kept only to re-select the frames.
        label = ceil.mode_group(mode)
        ladder_id = f"{night}|{tkey}|{label}|{filt}|{eg:g}"
        if ladder_id in done:
            continue
        # ... and the re-selection must use the same NULL-tolerant predicate
        # the mode grouping used: `readoutm = NULL` matches nothing in SQL.
        mode_cond, mode_params = mode_where(label)
        frames = con.execute(f"""
            SELECT f.path, f.basename, f.exptime FROM frames f
            WHERE f.night = ? AND f.target_key = ? AND {mode_cond}
              AND coalesce(f.filter, '') = ? AND coalesce(f.egain, -1) = ?
              AND f.is_canonical = 1 AND f.tree = 'rawimage'
              AND {SCIENCE_IMAGETYP} ORDER BY f.jd""",
            (night, tkey) + mode_params + (filt, eg)).fetchall()
        # Rung assignment: filename token beats a rounded-to-zero header.
        rungs: dict[float, list[str]] = {}
        for path, base, expt in frames:
            t = lin.effective_exptime(expt, base)
            tb = exptime_bin(t)
            if tb is None or tb <= 0:
                continue
            rungs.setdefault(tb, []).append(path)
        min_per = (1 if len(rungs) >= LONG_LADDER_RUNGS
                   else lin.MIN_FRAMES_PER_RUNG)
        rungs = {t: ps for t, ps in rungs.items() if len(ps) >= min_per}
        if len(rungs) < lin.MIN_RUNGS:
            continue
        print(f"[S2:linearity] {ladder_id}: rungs {sorted(rungs)}")
        rung_rows = []
        for t in sorted(rungs):
            fluxes, peaks = [], []
            for path in rungs[t][:8]:          # 8 frames per rung suffice
                try:
                    img, _ = read_image(archive, path)
                except Exception:
                    continue
                ph = lin.brightest_box_flux(np.asarray(img, dtype=np.float64))
                fluxes.append(ph["flux"])
                peaks.append(ph["peak_adu"])
            if fluxes:
                rung_rows.append((t, len(fluxes), float(np.median(fluxes)),
                                  float(np.median(peaks))))
        fit = lin.fit_ladder([r[0] for r in rung_rows],
                             [r[2] for r in rung_rows])
        if fit is None:
            con.execute("INSERT OR REPLACE INTO s2_linearity_ladders "
                        "(ladder_id, mode, night, target_key, n_rungs, "
                        "n_frames) VALUES (?,?,?,?,?,?)",
                        (ladder_id, label, night, tkey, len(rung_rows),
                         sum(r[1] for r in rung_rows)))
            con.commit()
            n_run += 1
            continue
        resid = dict(zip(fit["exptimes"], fit["resid_pct"]))
        con.execute("INSERT OR REPLACE INTO s2_linearity_ladders VALUES "
                    "(?,?,?,?,?,?,?,?)",
                    (ladder_id, label, night, tkey, fit["n_rungs"],
                     sum(r[1] for r in rung_rows), fit["rate_adu_per_s"],
                     fit["max_abs_resid_pct"]))
        con.executemany("INSERT OR REPLACE INTO s2_linearity_rungs VALUES "
                        "(?,?,?,?,?,?)",
                        [(ladder_id, t, n, fx, pk, resid.get(t))
                         for t, n, fx, pk in rung_rows])
        con.commit()
        n_run += 1
    print(f"[S2:linearity] done ({n_run} ladders this run; 0 = complete).")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: params  (distill everything into detector_params)
# ---------------------------------------------------------------------------
def cmd_params(con: sqlite3.Connection) -> int:
    """Derive s2_ceiling_modes, s2_ptc_fits and the detector_params table."""
    prov_meta = f"{S2_CODE_VERSION}; run_s2_campaign.py"
    # These three tables are pure distillations of the measurement tables —
    # rebuild them from scratch so a re-run never leaves stale rows behind.
    con.execute("DELETE FROM detector_params")
    con.execute("DELETE FROM s2_ceiling_modes")
    con.execute("DELETE FROM s2_ceiling_egain")

    def put(group, qty, val, unc, method, prov):
        con.execute("INSERT OR REPLACE INTO detector_params VALUES "
                    "(?,?,?,?,?,?)", (group, qty, val, unc, method,
                                      f"{prov}; {prov_meta}"))

    def egain_split(mode: str, veto: int) -> int:
        """Store per-egain near-ceiling frame-max stats; return group count.

        "Near-ceiling" = frames whose maximum reaches the mode's veto
        threshold (a saturated-star frame); grouping their maxima by the
        era's EGAIN exposes epoch drift the pooled histogram hides (the
        review's High Gain finding: two cleanly separated egain
        populations were being read as one wide "mound").
        """
        rows = con.execute("""
            SELECT round(e.egain, 3), c.max_adu
            FROM s2_ceiling_frames c
            JOIN frames f ON f.obs_rowid = c.obs_rowid
            JOIN eras e ON e.era_id = f.era_id
            WHERE c.mode = ? AND c.max_adu >= ? AND e.egain > 0""",
            (mode, veto)).fetchall()
        groups: dict[float, list[int]] = {}
        for eg, mx in rows:
            groups.setdefault(float(eg), []).append(int(mx))
        for eg, maxes in sorted(groups.items()):
            con.execute("INSERT OR REPLACE INTO s2_ceiling_egain VALUES "
                        "(?,?,?,?,?,?)",
                        (mode, eg, len(maxes), min(maxes),
                         float(np.median(maxes)), max(maxes)))
        # Only groups with a few frames count as evidence of a split.
        return sum(1 for m in groups.values() if len(m) >= 5)

    # --- ceilings from the accumulated histograms -------------------------
    modes = [r[0] for r in con.execute(
        "SELECT DISTINCT mode FROM s2_ceiling_hist")]
    for mode in modes:
        rows = con.execute("SELECT adu, count FROM s2_ceiling_hist "
                           "WHERE mode = ?", (mode,)).fetchall()
        n_frames = con.execute("SELECT count(*) FROM s2_ceiling_frames "
                               "WHERE mode = ? AND max_adu >= 0",
                               (mode,)).fetchone()[0]
        hist = np.zeros(65536, dtype=np.int64)
        for adu, count in rows:
            hist[adu] = count
        clip = ceil.find_clip(hist)
        hard_max = int(np.flatnonzero(hist)[-1]) if hist.sum() else None
        n_pixels = int(hist.sum())
        prov = f"{n_frames} science frames, {n_pixels} px histogram"
        if clip:
            veto = ceil.veto_threshold(clip["clip_adu"])
            bits = ceil.bit_depth_reading(clip["clip_adu"])
            con.execute("""INSERT OR REPLACE INTO s2_ceiling_modes
                (mode, n_frames, n_pixels, hard_max_adu, clip_adu,
                 spike_count, tail_level, ratio, veto_adu, bits,
                 adc_full_scale, unused_codes) VALUES
                (?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (mode, n_frames, n_pixels, hard_max, clip["clip_adu"],
                         clip["spike_count"], clip["tail_level"],
                         clip["ratio"], veto, bits["bits"],
                         bits["adc_full_scale"], bits["unused_codes"]))
            # Epoch honesty: the pooled histogram can hide egain drift
            # (the review's High Gain finding) — record per-egain stats and
            # say so in the method string when more than one epoch exists.
            n_epochs = egain_split(mode, veto)
            epoch_note = (f"; pooled over {n_epochs} egain epochs — "
                          "per-epoch stats in s2_ceiling_egain"
                          if n_epochs > 1 else "")
            put(mode, "ceiling_adu", clip["clip_adu"], 1.0,
                "pileup spike in science-frame pixel histogram"
                + epoch_note, prov)
            put(mode, "saturation_veto_adu", veto, None,
                f"floor({ceil.VETO_FRACTION} x ceiling, "
                f"{ceil.VETO_GRANULARITY_ADU} ADU) "
                "(exact derivation — no uncertainty)", prov)
            put(mode, "adc_bits", bits["bits"], None,
                "smallest ADC range containing the clip (consistency "
                "reading, exact by derivation; ADC-vs-full-well awaits the "
                "October hardware readback)", prov)
        else:
            # Histogram mound below the density threshold (sparse
            # saturation): fall back to per-frame-maximum clustering.
            maxes = [r[0] for r in con.execute(
                "SELECT max_adu FROM s2_ceiling_frames "
                "WHERE mode = ? AND max_adu > 0", (mode,))]
            cl = ceil.frame_max_cluster(maxes)
            diversity, n_pos = None, None
            if cl:
                # Diversity gate: the cluster is only a ceiling if its
                # members' maxima land at DIFFERENT places on the sensor
                # (Low Gain's fake cluster is one stable hot feature).
                lo = cl["clip_adu"] * (1 - ceil.CLUSTER_REL_WINDOW)
                hi = cl["clip_adu"] * (1 + ceil.CLUSTER_REL_WINDOW)
                pos = con.execute(
                    "SELECT max_y, max_x FROM s2_ceiling_frames "
                    "WHERE mode=? AND max_adu BETWEEN ? AND ? "
                    "AND max_y IS NOT NULL", (mode, lo, hi)).fetchall()
                diversity, n_pos = ceil.position_diversity(pos), len(pos)
                if n_pos < 10 or diversity < ceil.DIVERSITY_MIN_FRAC:
                    print(f"[S2:params] {mode}: frame-max cluster at "
                          f"{cl['clip_adu']} REJECTED (diversity "
                          f"{diversity:.2f} over {n_pos} positions) — "
                          "hot-pixel signature, not a ceiling.")
                    cl = {**cl, "rejected": True}
            if cl and not cl.get("rejected"):
                veto = ceil.veto_threshold(cl["clip_adu"])
                bits = ceil.bit_depth_reading(cl["clip_adu"])
                con.execute("""INSERT OR REPLACE INTO s2_ceiling_modes
                    (mode, n_frames, n_pixels, hard_max_adu, clip_adu,
                     veto_adu, bits, adc_full_scale, unused_codes,
                     cluster_adu, cluster_diversity, cluster_n_pos) VALUES
                    (?,?,?,?,?,?,?,?,?,?,?,?)""",
                            (mode, n_frames, n_pixels, hard_max,
                             cl["clip_adu"], veto, bits["bits"],
                             bits["adc_full_scale"], bits["unused_codes"],
                             cl["clip_adu"], diversity, n_pos))
                n_epochs = egain_split(mode, veto)
                epoch_note = (f"; pooled over {n_epochs} egain epochs — "
                              "per-epoch stats in s2_ceiling_egain"
                              if n_epochs > 1 else "")
                put(mode, "ceiling_adu", cl["clip_adu"], cl["mad_adu"],
                    f"frame-maximum cluster ({100 * cl['cluster_frac']:.0f}%"
                    f" of frames; argmax diversity {diversity:.2f}; "
                    "unc = MAD-sigma of cluster members" + epoch_note, prov)
                put(mode, "saturation_veto_adu", veto, None,
                    f"floor({ceil.VETO_FRACTION} x ceiling, "
                    f"{ceil.VETO_GRANULARITY_ADU} ADU) "
                    "(exact derivation — no uncertainty)", prov)
                put(mode, "adc_bits", bits["bits"], None,
                    "smallest ADC range containing the clip (consistency "
                    "reading, exact by derivation; ADC-vs-full-well awaits "
                    "the October hardware readback)", prov)
            else:
                con.execute("""INSERT OR REPLACE INTO s2_ceiling_modes
                    (mode, n_frames, n_pixels, hard_max_adu, cluster_adu,
                     cluster_diversity, cluster_n_pos) VALUES
                    (?,?,?,?,?,?,?)""",
                            (mode, n_frames, n_pixels, hard_max,
                             cl["clip_adu"] if cl else None, diversity,
                             n_pos))
                put(mode, "observed_max_adu", hard_max, None,
                    "no pileup detected; observed maximum only (exact "
                    "observed value — no uncertainty)"
                    + ("; frame-max cluster rejected as hot pixel"
                       if cl else ""), prov)

    # --- PTC --------------------------------------------------------------
    # Straight-line fits per (mode, kind) go into s2_ptc_fits as recorded
    # facts.  Neither slope is adopted as THE gain: the dark slope is
    # biased shallow (part of the hot-pixel population barely fluctuates
    # between consecutive frames, so its variance grows sub-Poisson) and
    # the light slope is biased steep (sub-pixel scene motion inflates the
    # difference variance wherever the image has gradients).  The two
    # therefore BRACKET the true gain; the header EGAIN sits inside the
    # bracket and is recorded as the nominal value.  Read noise, by
    # contrast, is measured cleanly: the variance floor of the shortest
    # darks has no scene and no dark-current shot noise.
    con.execute("DELETE FROM s2_ptc_fits")
    fits_by_mode: dict[tuple, dict] = {}
    for mode, kind in con.execute(
            "SELECT DISTINCT mode, kind FROM s2_ptc_points"):
        pts = con.execute("SELECT level, var, n_pix FROM s2_ptc_points "
                          "WHERE mode = ? AND kind = ?",
                          (mode, kind)).fetchall()
        fit = ptc.fit_ptc([p[0] for p in pts], [p[1] for p in pts],
                          [p[2] for p in pts])
        if fit is None:
            continue
        fits_by_mode[(mode, kind)] = fit
        con.execute("INSERT OR REPLACE INTO s2_ptc_fits VALUES "
                    "(?,?,?,?,?,?,?,?,?,?)",
                    (mode, kind, fit["gain_e_per_adu"], fit["gain_err"],
                     fit["read_noise_adu"], fit["read_noise_adu_err"],
                     fit["read_noise_e"], fit["slope"], fit["intercept"],
                     fit["n_points"]))
    # Nominal header gain per mode: the EGAIN of the mode's biggest era.
    egain_of = {}
    for mode, eg in con.execute("""
            SELECT readoutm, egain FROM eras WHERE egain > 0
            GROUP BY readoutm HAVING n_frames = max(n_frames)"""):
        egain_of[mode] = eg
    dark_pts_of: dict[str, list] = {}
    for mode, in con.execute(
            "SELECT DISTINCT mode FROM s2_ptc_points WHERE kind = 'dark'"):
        dark_pts_of[mode] = con.execute(
            "SELECT exptime, level, var FROM s2_ptc_points "
            "WHERE mode = ? AND kind = 'dark'", (mode,)).fetchall()
    brackets: dict[str, tuple[float, float]] = {}
    # sorted(): "High Gain" precedes "High Gain StackPro", so the StackPro
    # electron conversion can fall back on the already-computed HG bracket.
    for mode, dpts in sorted(dark_pts_of.items()):
        rn = ptc.read_noise_from_dark_points(dpts)
        if rn is None:
            continue
        prov = (f"variance floor of {rn['exptime']:g}s dark pairs, "
                f"night {PTC_NIGHT}, {rn['n_points']} level bins")
        # Gain bracket FIRST (read_noise_e's uncertainty depends on it):
        # dark slope = upper bound, sky-level light slope = lower bound,
        # header EGAIN recorded as the nominal in between.
        fd = fits_by_mode.get((mode, "dark"))
        if fd:
            put(mode, "gain_upper_bound_e_per_adu", fd["gain_e_per_adu"],
                fd["gain_err"], "dark-pair PTC slope (sub-Poisson-biased: "
                "quiet hot pixels)", f"{fd['n_points']} dark points")
        lpts = con.execute(
            "SELECT level, var, n_pix FROM s2_ptc_points WHERE mode = ? "
            "AND kind = 'light' AND level < ?",
            (mode, ptc.GAIN_LOWER_BOUND_LEVEL_ADU)).fetchall()
        fl = ptc.fit_ptc([p[0] for p in lpts], [p[1] for p in lpts],
                         [p[2] for p in lpts])
        if fl:
            put(mode, "gain_lower_bound_e_per_adu", fl["gain_e_per_adu"],
                fl["gain_err"], "sky-level light-pair PTC slope (motion-"
                "inflated variance)", f"{fl['n_points']} light points "
                f"below {ptc.GAIN_LOWER_BOUND_LEVEL_ADU:g} ADU")
        if fd and fl:
            brackets[mode] = (fl["gain_e_per_adu"], fd["gain_e_per_adu"])
        # Read noise in ADU: the floor's statistical error alone understates
        # it (review finding) — fold in the measured dark-current shot term
        # still hiding inside the shortest floor (8 s darks are not 0 s).
        shot = ptc.dark_shot_fraction(dpts)
        rn_stat = rn["read_noise_adu_err"] or 0.0
        rn_shot = shot["rn_bias_adu"] if shot else 0.0
        rn_unc = float(np.hypot(rn_stat, rn_shot))
        shot_note = (f" + dark-shot floor systematic ({shot['t_short']:g}s "
                     f"floor carries {100 * shot['frac_of_floor']:.1f}% of "
                     "the RN variance in dark current, measured from the "
                     f"{shot['t_short']:g}s-vs-{shot['t_long']:g}s floors)"
                     if shot else "")
        put(mode, "read_noise_adu", rn["read_noise_adu"], rn_unc,
            "shortest-dark pair variance floor; unc = statistical"
            + shot_note, prov)
        put(mode, "bias_offset_adu", rn["offset_adu"], rn["offset_adu_err"],
            "minimum dark-pair level (the bias pedestal); unc = half-spread "
            "of the floor bins' levels", prov)
        if mode in egain_of:
            # Electron conversion: the value uses the NOMINAL header EGAIN,
            # but the campaign only BRACKETS the true gain — so the honest
            # uncertainty is the bracket's half-width propagated through
            # RN(ADU), not the (negligible) ADU statistical error (review
            # finding: +/-0.001 e- was false precision).  The High Gain
            # bracket covers the StackPro family too (same sensor and
            # sub-read gain; SP's own sky never reaches the fit window).
            br = brackets.get(mode) or (
                brackets.get("High Gain")
                if mode.startswith("High Gain") else None)
            if br:
                e_lo, e_hi = (rn["read_noise_adu"] * br[0],
                              rn["read_noise_adu"] * br[1])
                put(mode, "read_noise_e",
                    rn["read_noise_adu"] * egain_of[mode],
                    (e_hi - e_lo) / 2.0,
                    "RN(ADU) x nominal header EGAIN; unc = half-width of "
                    f"RN(ADU) x measured gain bracket [{br[0]:.2f}, "
                    f"{br[1]:.2f}] e-/ADU (dominant systematic until the "
                    "October flat-field PTC)", prov)
            else:
                put(mode, "read_noise_e",
                    rn["read_noise_adu"] * egain_of[mode],
                    rn_unc * egain_of[mode],
                    "RN(ADU) x nominal header EGAIN (statistical only — "
                    "no gain bracket measured for this mode)", prov)
        if mode in egain_of:
            put(mode, "gain_e_per_adu_nominal", egain_of[mode], None,
                "header EGAIN (inside the archival PTC bracket; hardware "
                "PTC = October confirmation item; header constant — no "
                "measurement uncertainty)", "eras table")
    # StackPro N_sub: three independent ratios against plain High Gain
    # (bias offset, read-noise variance, saturation ceiling — all x N_sub
    # if StackPro frames are sums of N_sub sub-exposures).
    ceilings_now = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes "
        "WHERE clip_adu IS NOT NULL"))
    if "High Gain StackPro" in dark_pts_of and "High Gain" in dark_pts_of:
        sig = ptc.stackpro_signature(
            dark_pts_of["High Gain StackPro"], dark_pts_of["High Gain"],
            ceilings_now.get("High Gain StackPro"),
            ceilings_now.get("High Gain"))
        if sig:
            put("High Gain StackPro", "nsub", sig["nsub"], sig["max_misfit"],
                "consensus of offset/read-noise-variance/ceiling ratios "
                "vs High Gain",
                "; ".join(f"{k}={v:.2f}" for k, v in sig.items()
                          if k.endswith("_ratio")))
    # Amp glow: the archive's longest darks, hottest-corner excess.
    for mode, expt, corner, edge in con.execute("""
            SELECT mode, exptime, hottest_corner_excess, edge_excess
            FROM (SELECT mode, exptime,
                         hottest_corner_excess, edge_excess FROM s2_ampglow)
            GROUP BY mode, exptime
            HAVING hottest_corner_excess = max(hottest_corner_excess)"""):
        n = con.execute("SELECT count(*) FROM s2_ampglow WHERE mode=?",
                        (mode,)).fetchone()[0]
        # Spread of the same statistic across this mode's darks at the same
        # exposure = the honest uncertainty of a max-over-darks value.
        vals = [r[0] for r in con.execute(
            "SELECT hottest_corner_excess FROM s2_ampglow "
            "WHERE mode = ? AND exptime = ?", (mode, expt))]
        unc = ((max(vals) - min(vals)) / 2.0) if len(vals) > 1 else None
        put(mode, f"amp_glow_corner_adu_{expt:g}s", corner, unc,
            "hottest-corner median minus center median, longest darks; "
            "unc = half-range across darks",
            f"max over {n} darks; edge-band excess {edge:.1f} ADU")

    # --- the empirical noise model ----------------------------------------
    # Pure distillation of s2_noise_points: rebuild from scratch so a re-run
    # never leaves a stale curve behind.
    con.execute("DELETE FROM s2_noise_curve")
    for mode, in con.execute("SELECT DISTINCT mode FROM s2_noise_points"):
        pts = con.execute("SELECT level, var, n_pix, pair_id "
                          "FROM s2_noise_points WHERE mode = ?",
                          (mode,)).fetchall()
        curve = noisemod.empirical_noise_curve(pts)
        if not curve:
            continue
        con.executemany("INSERT OR REPLACE INTO s2_noise_curve VALUES "
                        "(?,?,?,?,?,?,?,?,?)",
                        [(mode, i, c["level"], c["var"], c["var_mad"],
                          c["sigma"], c["n_points"], c["n_pairs"], c["n_pix"])
                         for i, c in enumerate(curve)])
        n_pairs = con.execute("SELECT count(*) FROM s2_noise_pairs "
                              "WHERE mode = ? AND n_points > 0",
                              (mode,)).fetchone()[0]
        n_scenes = con.execute(
            "SELECT count(DISTINCT night || '|' || target_key) "
            "FROM s2_noise_pairs WHERE mode = ? AND n_points > 0",
            (mode,)).fetchone()[0]
        prov = (f"{n_pairs} same-scene science pairs across {n_scenes} "
                f"(night, target) scenes; {len(curve)} measured level bins "
                f"spanning {curve[0]['level']:,.0f}-{curve[-1]['level']:,.0f}"
                " ADU (table: s2_noise_curve)")
        fl = noisemod.noise_floor(curve)
        if fl:
            put(mode, "noise_floor_adu", fl["floor_sigma_adu"],
                fl["floor_sigma_err_adu"],
                "MEASURED sigma at the bottom of the counts-vs-variance "
                f"curve (level {fl['level_adu']:,.0f} ADU, "
                f"{fl['n_bins']} bins); this is the whole floor a science "
                "frame carries — read noise PLUS dark PLUS bias structure — "
                "not a zero-second read noise; unc = half-spread across the "
                "pooled floor bins", prov)
            # The crossover contrasts the floor against the rest of the
            # curve, so it means nothing when the floor window IS the whole
            # curve (a mode measured at a single sky level — the 5 MHz
            # iKon's science pairs all sit on a flat ~6,000 ADU sky).
            x = (None if fl["is_whole_curve"]
                 else noisemod.crossover_level(curve, fl["floor_var_adu2"]))
            if x is not None:
                put(mode, "noise_crossover_adu", x, None,
                    f"level where the MEASURED variance reaches "
                    f"{noisemod.CROSSOVER_VAR_MULTIPLE:g}x the floor "
                    "(interpolated between measured bins — no gain and no "
                    "Poisson law assumed; exact by interpolation)", prov)
        slope = noisemod.curve_shape_index(curve)
        if slope is not None:
            put(mode, "noise_curve_logslope", slope, None,
                "log-log slope of measured variance vs level: 1.0 = "
                "shot-noise-dominated, ~0 = floor-dominated, >1 = the "
                "pair difference is also measuring scene motion, so the "
                "bright end is an UPPER bound on detector noise "
                "(descriptive statistic — no uncertainty attached)", prov)

    # --- reconstruction ---------------------------------------------------
    # Backfill the F-D degeneracy diagnostic from the stored npz products
    # (cheap: eight small files; no archive pixels re-read).
    for era_id, npz_path in con.execute(
            "SELECT era_id, npz_path FROM s2_recon_eras "
            "WHERE npz_path IS NOT NULL"):
        try:
            npz = np.load(REPO_ROOT / npz_path, allow_pickle=True)
            corr = rec.flat_dark_correlation(npz["F"], npz["D"])
        except Exception as e:                         # pragma: no cover
            print(f"[S2:params] era {era_id}: fd_corr unreadable: {e}")
            continue
        con.execute("UPDATE s2_recon_eras SET fd_corr = ? WHERE era_id = ?",
                    (None if np.isnan(corr) else float(corr), era_id))
    for (era_id, flat_med, f_mad, dark_med, d_mad, rms_med, t_rms,
         t_mad) in con.execute(
            """SELECT era_id, flat_median, flat_mad_sigma, dark_median,
                      dark_mad_sigma, rms_median, truth_resid_rms,
                      truth_resid_mad FROM s2_recon_eras
               WHERE n_pairs_used >= ?""", (rec.RECON_MIN_PAIRS,)):
        grp = f"era {era_id}"
        # The MAD-sigmas are the per-pixel SPREAD of the recovered values,
        # not a standard error of the median — stated in the method string
        # so a programmatic reader knows which kind of number it holds.
        put(grp, "recon_flat_median", flat_med, f_mad,
            "per-pixel raw-vs-reduced slope; unc = MAD-sigma of per-pixel "
            "values (spread incl. F-D degeneracy noise, not a standard "
            "error — see s2_recon_eras.fd_corr)", "s2_recon_eras")
        put(grp, "recon_dark_median_adu", dark_med, d_mad,
            "per-pixel raw-vs-reduced intercept; unc = MAD-sigma of "
            "per-pixel values (spread, not a standard error)",
            "s2_recon_eras")
        put(grp, "recon_residual_rms_adu", rms_med, None,
            "median per-pixel line-fit RMS (summary statistic — no "
            "uncertainty attached)", "s2_recon_eras")
        if t_rms is not None:
            put(grp, "recon_vs_master_rms_adu", t_rms, t_mad,
                "reconstructed D vs archived master (offset removed); "
                "unc = MAD-sigma of the same residuals",
                "s2_recon_eras.truth_master")

    # --- linearity --------------------------------------------------------
    # Per mode: the cleanest ladder's maximum |residual| over UNSATURATED
    # rungs only.  A rung whose peak pixel sits above the mode's veto
    # threshold measures the ceiling (flux loss to clipping), not detector
    # linearity, so it is excluded from the linearity statistic — its
    # residual still lives in s2_linearity_rungs as the ceiling cross-check.
    # min over ladders: the CLEANEST ladder bounds the mode's real
    # non-linearity (dirtier ladders add clouds/tracking, not detector).
    best_by_mode: dict[str, tuple] = {}
    for ladder_id, mode in con.execute(
            "SELECT ladder_id, mode FROM s2_linearity_ladders "
            "WHERE rate_adu_per_s IS NOT NULL"):
        veto = con.execute("SELECT veto_adu FROM s2_ceiling_modes "
                           "WHERE mode = ?", (mode,)).fetchone()
        veto_adu = veto[0] if veto and veto[0] is not None else float("inf")
        rungs = con.execute(
            "SELECT resid_pct, peak_med FROM s2_linearity_rungs "
            "WHERE ladder_id = ? AND resid_pct IS NOT NULL",
            (ladder_id,)).fetchall()
        clean_signed = [r for r, pk in rungs if pk is None or pk < veto_adu]
        clean = [abs(r) for r in clean_signed]
        if len(clean) < 3:
            continue
        worst = max(clean)
        # Scatter of the clean rungs' signed residuals: the honest scale of
        # a single-ladder bound (sky-transparency drift folds in here).
        scatter = float(np.std(clean_signed, ddof=1))
        if mode not in best_by_mode or worst < best_by_mode[mode][0]:
            best_by_mode[mode] = (worst, ladder_id, len(clean), scatter)
    for mode, (worst, ladder_id, n_clean, scatter) in best_by_mode.items():
        put(mode, "linearity_max_dev_pct", worst, scatter,
            "best archival exposure ladder, max |residual| over "
            "unsaturated rungs; SINGLE-LADDER consistency bound (the "
            "median-rate fit zeroes one rung by construction and residuals "
            "include sky-transparency drift — not a measured detector "
            "non-linearity); unc = std of clean-rung residuals",
            f"{ladder_id} ({n_clean} rungs below the saturation veto)")
    con.commit()
    n = con.execute("SELECT count(*) FROM detector_params").fetchone()[0]
    print(f"[S2:params] detector_params now holds {n} rows.")
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=("Run the S2 detector campaign (ceiling memo, PTC, "
                     "master reconstruction, linearity) against the RLMT "
                     "archive. Resumable: every subcommand records finished "
                     "work in the manifest DB and skips it when re-invoked."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("stage", choices=["ceiling", "ceilpos", "ptc", "noise",
                                     "reconstruct", "linearity", "params",
                                     "report"])
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="S0/S0b manifest database to augment")
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE,
                   help="archive root holding the pixel trees")
    p.add_argument("--batch", type=int, default=250,
                   help="max frames/pairs/ladders processed this invocation")
    p.add_argument("--era", type=int, default=None,
                   help="reconstruct: process only this era")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.manifest.exists():
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 2
    if args.stage == "report":
        from rlmt_diagnostics import report_s2
        path = report_s2.render_report(args.manifest)
        print(f"[S2] report -> {path}")
        return 0
    with closing(open_db(args.manifest)) as con:
        ensure_tables(con)
        write_meta(con, args.manifest)
        if args.stage == "ceiling":
            return cmd_ceiling(con, args.archive, args.batch)
        if args.stage == "ceilpos":
            return cmd_ceilpos(con, args.archive, args.batch)
        if args.stage == "ptc":
            return cmd_ptc(con, args.archive, args.batch)
        if args.stage == "noise":
            return cmd_noise(con, args.archive, args.batch)
        if args.stage == "reconstruct":
            return cmd_reconstruct(con, args.archive, args.era)
        if args.stage == "linearity":
            return cmd_linearity(con, args.archive, args.batch)
        if args.stage == "params":
            return cmd_params(con)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
