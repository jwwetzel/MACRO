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
from rlmt_diagnostics import badpix                           # noqa: E402
from rlmt_diagnostics import camera as cam                    # noqa: E402
from rlmt_diagnostics import ceiling as ceil                  # noqa: E402
from rlmt_diagnostics import flatptc                          # noqa: E402
from rlmt_diagnostics import linearity as lin                 # noqa: E402
from rlmt_diagnostics import noise as noisemod                # noqa: E402
from rlmt_diagnostics import ptc                              # noqa: E402
from rlmt_diagnostics import reconstruct as rec               # noqa: E402
from rlmt_diagnostics import saturation as satmod             # noqa: E402
from rlmt_diagnostics import starphot                         # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: The CV time-series photometry database (read-only here).  Its
#: ``cv_detections`` table holds, for every matched star in every frame,
#: the aperture flux AND the peak pixel — which makes ~7,800 frames of
#: comparison-star photometry a ready-made linearity experiment
#: (stage ``peaklin``; committee finding DE.F3).
DEFAULT_CV_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
PRODUCTS_DIR = REPO_ROOT / "products" / "detector"

#: Where S2 COMPUTES.  Since the 2026-10-03 review S2 no longer writes the
#: shared manifest while it works: every s2_* table and ``detector_params``
#: live in this stage-owned database, with the manifest attached READ-ONLY
#: for the S0/S0b tables the stages query (``frames``, ``eras``, ...).  The
#: ``promote`` subcommand copies the finished tables into the manifest —
#: one deliberate step, run at the integrated rebuild, instead of eight
#: stages each holding the shared write lock.
DEFAULT_DETECTOR_DB = PRODUCTS_DIR / "detector.sqlite"

#: Schema name the manifest is attached under in a work connection.
MANIFEST_SCHEMA = "mf"

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

# ---------------------------------------------------------------------------
# Flat-pair photon transfer (stage ``flats``) and sky-pair photon transfer
# (stage ``skypairs``).
# ---------------------------------------------------------------------------
#: Archive trees searched for calibration frames.  The census is a RULE,
#: not a list of directories, so the flats the October 2026 re-opening
#: delivers for the QHY600 are picked up by the same command.
FLAT_TREES = ("calib", "Calibrations", "iKon", "rawimage", "cmos_tests",
              "latency")

#: Trees whose DARK frames may serve as zero-signal pairs (read noise) for
#: a configuration that has no bias frames (the AC4040 has none at all).
ZERO_DARK_TREES = ("calib", "Calibrations", "cmos_tests", "iKon")

#: Two flats are a pair only if taken within this many days (10 minutes):
#: a twilight sky or a warming panel is the same scene only briefly.
FLAT_MAX_GAP_DAYS = 10.0 / (24.0 * 60.0)

#: Pairs read per flat group (one configuration x directory x filter x
#: exposure x temperature).  Spread evenly through the group's sequence.
#: Four per group is plenty: one pair of 4096^2 frames already fixes its
#: variance to 0.03%, and what limits the gain is pair-to-pair and
#: group-to-group systematics, which more GROUPS sample and more pairs per
#: group do not.  (The first AC4040 groups were read at 12 per group
#: before this was lowered; those pairs are kept and used.)
FLAT_MAX_PAIRS_PER_GROUP = 4

#: Zero-signal (bias or shortest-dark) pairs read per configuration and
#: temperature group.
ZERO_MAX_PAIRS = 12

#: Set-point grouping: CCD-TEMP rounded to the nearest this many degrees
#: (the ASI was run at -10 C and at 0 C; the two are fitted apart as well
#: as together — committee finding DE.F5).
TEMP_GROUP_STEP_C = 5.0

#: A "flat" whose level stands less than this many zero-level sigmas above
#: the bias is not a flat (25 frames typed 'Flat Field' in
#: Calibrations/2025-01 are 2 s darks).
FLAT_MIN_SIGNAL_SIGMAS = 10.0

#: A pair whose tile-to-tile level ratio scatters by more than this (MAD-
#: sigma, relative) did not see the same illumination pattern twice.
FLAT_MAX_RATIO_SCATTER = 0.03

#: Sky pairs read per configuration.  Scenes are stratified on exposure
#: time (the only sky-level proxy the manifest holds) so the level axis is
#: spanned, not sampled at one sky brightness.
SKY_MAX_PAIRS_PER_CONFIG = 12

#: A configuration needs at least this many canonical science frames
#: before a sky-pair PTC is attempted for it ...
SKY_MIN_CONFIG_FRAMES = 1000

#: ... except these, attempted whatever their size: the ASI's 61 UNBINNED
#: frames are the only direct measurement of the native-pixel gain, which
#: is what the 2x2 average-vs-sum question (D2) is decided against.
SKY_ALWAYS_CONFIGS = ("ASI Mode0 1x1",)

#: Filters never used for sky pairs: through a grism the "sky" is a
#: dispersed lozenge with sharp edges, not a flat field.
SKY_EXCLUDED_FILTERS = ("hrg", "lrg", "hagrism", "oggrism", "hag", "ogg")

# ---------------------------------------------------------------------------
# Star-based differential linearity (stage ``starlin``).
# ---------------------------------------------------------------------------
#: The dedicated star sets the committee named (DE.F3, SN-S2-linearity).
#: Each family is a RULE selecting frames; frames are then grouped into
#: sets by (family, night, filter) — one field, one filter, several
#: exposure times and/or readout modes.
#:
#: * ``albireo``          — 2023-06-07 cmos_tests, the Albireo pointing:
#:   High Gain I at 8/32/64/128 s plus StackPro I at 32/64 s (an exposure
#:   ladder across a 16x range AND a mode comparison).
#: * ``albireo_flanking`` — the flanking field: ~11 High Gain and 10
#:   StackPro frames, all I 32 s (StackPro, whose sub-reads never approach
#:   the clip, is the linear reference for the single read).
#: * ``sn2023ixf``        — the 0.5 s and 2 s frames of 2023-05-23/24, per
#:   filter: a 4x exposure ratio on one field minutes apart.
STARLIN_FAMILIES: dict[str, str] = {
    "albireo": ("f.tree = 'cmos_tests' AND f.basename LIKE 'Albi%' "
                "AND f.basename NOT LIKE '%flanking%' "
                "AND f.imagetyp LIKE 'Light%'"),
    "albireo_flanking": ("f.tree = 'cmos_tests' "
                         "AND f.basename LIKE '%flanking%' "
                         "AND f.imagetyp LIKE 'Light%'"),
    "sn2023ixf": ("f.tree = 'rawimage' AND f.is_canonical = 1 "
                  "AND f.target_key = '2023ixf' AND f.imagetyp LIKE 'Light%' "
                  "AND f.readoutm LIKE 'High Gain%' AND f.exptime <= 2.1 "
                  "AND f.night IN ('2023-05-23', '2023-05-24') "
                  "AND f.filter IN ('G', 'R', 'I')"),
    # The QHY600 (camera on the telescope since 2026-03): no flats, no
    # ladders — but V426 Oph is a rich field observed all night at one
    # exposure, so seeing alone sweeps each star's peak across the scale.
    # Fixed-exposure time series, analysed like the CV comparison stars.
    "qhy_v426oph": ("f.tree = 'rawimage' AND f.is_canonical = 1 "
                    "AND f.target_key = 'v426oph' AND f.readoutm = 'Fast' "
                    "AND f.night = '2026-05-31' AND f.filter = 'g' "
                    "AND f.exptime BETWEEN 9.9 AND 10.1"),
    # V426 Oph at 10 s reaches only 0.3 of scale; V2400 Oph at 30 s (also
    # a galactic field, observed all night) carries stars to the clip.
    "qhy_v2400oph": ("f.tree = 'rawimage' AND f.is_canonical = 1 "
                     "AND f.target_key = 'v2400oph' AND f.readoutm = 'Fast' "
                     "AND f.night = '2026-06-02' AND f.filter = 'g' "
                     "AND f.exptime BETWEEN 29.9 AND 30.1"),
}

#: A FIXED-exposure set (one mode, one exposure) is sampled down to this
#: many frames, evenly spread through the night: what it needs is seeing
#: and airmass variety, which an even spread keeps and 1,500 frames do not
#: improve on.
STARLIN_MAX_FRAMES_FIXED = 30

#: Stars tracked per set (the brightest in the set's reference frame).
#: 400 keeps the analysis on stars measured at high signal-to-noise.  A
#: deeper list (1,500) was TRIED so that the 64 s and 128 s Albireo frames
#: — in which every one of the 400 brightest stars is above the reference
#: regime — would have faint stars to anchor their frame factors.  It
#: failed its own control: with faint stars in the glare of a saturated
#: third-magnitude double carrying the normalisation across a 16x exposure
#: range, the reference-regime bins themselves read +9% where they must
#: read zero.  Those two frames therefore contribute nothing (they have
#: fewer than the minimum number of reference stars), and
#: ``_params_caps`` refuses any source whose reference bins are not flat.
STARLIN_MAX_STARS = 400

# ---------------------------------------------------------------------------
# Peak-at-target saturation census (stage ``peakcensus``).
# ---------------------------------------------------------------------------
#: Targets censused: target_key -> (RA deg, Dec deg, ICRS).  T CrB is the
#: committee's named case (OA.E3: "all 224 imaging frames").
CENSUS_TARGETS: dict[str, tuple[float, float]] = {
    "tcrb": (239.875675, 25.920172),
}

#: Filters that are grisms by NAME (never direct images).  Slots whose
#: identity is only known from pixels ('6', 'W') are handled through the
#: S2c ``frame_dispersion`` verdict, not through this list.
CENSUS_GRISM_FILTERS = ("hrg", "lrg", "hagrism", "oggrism", "hag", "ogg")

# ---------------------------------------------------------------------------
# Comparison-star residual vs own peak (stage ``peaklin``).
# ---------------------------------------------------------------------------
#: A CV series enters the peak-linearity analysis with at least this many
#: matched frames and this many matched stars.
PEAKLIN_MIN_FRAMES = 30
PEAKLIN_MIN_STARS = 10

#: Brightest stars kept per series (the linearity question lives at the
#: bright end; the faint ones only anchor the frame factors, and a few
#: hundred of those are plenty).
PEAKLIN_MAX_STARS = 600

#: A star must be measured in at least this fraction of a series' frames.
PEAKLIN_MIN_COVERAGE = 0.5

#: Null-injection trials per mode (on the mode's largest series).
PEAKLIN_INJECTION_TRIALS = 40

# ---------------------------------------------------------------------------
# Bad-pixel masks (stage ``badpix``).
# ---------------------------------------------------------------------------
#: Science frames accumulated per mask (distinct night x target scenes).
BADPIX_SCIENCE_FRAMES = 36

#: Dark frames accumulated for the independent dark-based check, and for
#: the temporal-noise (RTS) mask.
BADPIX_DARK_FRAMES = 10

#: Science frames shorter than this are not used for a mask: a warm pixel
#: needs integration time to stand out.
BADPIX_MIN_EXPTIME_S = 20.0

#: A camera/geometry/temperature group needs this many candidate science
#: frames before a mask is attempted.
BADPIX_MIN_GROUP_FRAMES = 1000

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


def open_work_db(target: Path, manifest: Path) -> sqlite3.Connection:
    """Open the S2 work database with the manifest attached read-only.

    ``main`` is the stage-owned database (:data:`DEFAULT_DETECTOR_DB`); the
    manifest is attached as schema :data:`MANIFEST_SCHEMA` with
    ``mode=ro``, so SQLite itself refuses any write to it.  Unqualified
    table names resolve to ``main`` first and to the attached manifest
    second — which is exactly the split wanted: ``frames`` / ``eras`` /
    ``calib_frames`` are found in the manifest, every ``s2_*`` table and
    ``detector_params`` in the work database, and no stage needed a single
    SQL change.

    When ``target`` IS the manifest (the pre-review behaviour, still
    available for a manifest-only rebuild) nothing is attached.
    """
    target = Path(target).resolve()
    manifest = Path(manifest).resolve()
    if target == manifest:
        return open_db(manifest)
    target.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(target.as_uri(), uri=True, timeout=300.0)
    con.execute("PRAGMA busy_timeout = 300000")
    con.execute(f"ATTACH DATABASE ? AS {MANIFEST_SCHEMA}",
                (manifest.as_uri() + "?mode=ro",))
    return con


def has_manifest_schema(con: sqlite3.Connection) -> bool:
    """True when the connection is a work DB with the manifest attached."""
    return any(r[1] == MANIFEST_SCHEMA
               for r in con.execute("PRAGMA database_list"))


def s2_table_names(con: sqlite3.Connection, schema: str) -> list[str]:
    """The S2-owned tables present in one schema of a connection.

    S2 owns every ``s2_*`` table (NOT ``s2c_*`` — that is the dispersion
    stage) and ``detector_params``.
    """
    rows = con.execute(
        f"SELECT name FROM {schema}.sqlite_master WHERE type = 'table'")
    return sorted(r[0] for r in rows
                  if r[0] == "detector_params" or r[0].startswith("s2_"))


def cmd_seed(con: sqlite3.Connection) -> int:
    """Copy the manifest's existing S2 tables into the work database.

    Run once, before the first stage, so the resumable stages find the
    work that S2 v1.2 already recorded (ceiling histograms, PTC pairs,
    reconstruction eras ...) and skip it instead of re-reading thousands
    of frames.  Rows already present in the work database are never
    overwritten (``INSERT OR IGNORE`` on keyed tables; an un-keyed table
    is seeded only while empty), so re-running ``seed`` is harmless.
    """
    if not has_manifest_schema(con):
        print("[S2:seed] target is the manifest itself — nothing to seed.")
        return 0
    main_tables = set(s2_table_names(con, "main"))
    for name in s2_table_names(con, MANIFEST_SCHEMA):
        if name not in main_tables:
            # A table this code version does not declare: carry it over
            # verbatim rather than silently dropping measurements.
            con.execute(f"CREATE TABLE main.{name} AS "
                        f"SELECT * FROM {MANIFEST_SCHEMA}.{name} WHERE 0")
        cols_main = [r[1] for r in con.execute(
            f"PRAGMA main.table_info({name})")]
        cols_mf = [r[1] for r in con.execute(
            f"PRAGMA {MANIFEST_SCHEMA}.table_info({name})")]
        common = [c for c in cols_mf if c in cols_main]
        has_pk = any(r[5] for r in con.execute(
            f"PRAGMA main.table_info({name})"))
        n_before = con.execute(
            f"SELECT count(*) FROM main.{name}").fetchone()[0]
        if not has_pk and n_before:
            print(f"[S2:seed]   {name}: {n_before} rows already present "
                  "(un-keyed table) — left alone.")
            continue
        collist = ", ".join(common)
        con.execute(f"INSERT OR IGNORE INTO main.{name} ({collist}) "
                    f"SELECT {collist} FROM {MANIFEST_SCHEMA}.{name}")
        n_after = con.execute(
            f"SELECT count(*) FROM main.{name}").fetchone()[0]
        print(f"[S2:seed]   {name}: {n_after - n_before} rows seeded "
              f"({n_after} total).")
    con.commit()
    return 0


def cmd_promote(target: Path, manifest: Path) -> int:
    """Publish the work database's S2 tables into the manifest.

    The ONE step in which S2 writes the shared manifest.  Each S2 table is
    replaced wholesale (drop, re-create with the work database's exact
    schema, copy rows) inside a single transaction, so the manifest either
    holds the complete new S2 state or the complete old one.  Run by the
    chair at the integrated rebuild; never by a work package.
    """
    target, manifest = Path(target).resolve(), Path(manifest).resolve()
    if target == manifest:
        print("[S2:promote] target is the manifest — nothing to promote.")
        return 0
    if not target.exists():
        print(f"ERROR: work database not found: {target}", file=sys.stderr)
        return 2
    con = sqlite3.connect(manifest, timeout=300.0)
    try:
        con.execute("PRAGMA busy_timeout = 300000")
        con.execute("ATTACH DATABASE ? AS det", (str(target),))
        con.execute("BEGIN IMMEDIATE")
        for name in s2_table_names(con, "det"):
            ddl = con.execute(
                "SELECT sql FROM det.sqlite_master WHERE type = 'table' "
                "AND name = ?", (name,)).fetchone()[0]
            con.execute(f"DROP TABLE IF EXISTS main.{name}")
            con.execute(ddl)                 # unqualified -> created in main
            con.execute(f"INSERT INTO main.{name} SELECT * FROM det.{name}")
            n = con.execute(f"SELECT count(*) FROM main.{name}").fetchone()[0]
            print(f"[S2:promote]   {name}: {n} rows")
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()
    print(f"[S2:promote] manifest updated from {target}.")
    return 0


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
    -- ------------------------------------------------------------------
    -- Flat-pair photon transfer per camera/configuration (review F-4).
    -- s2_flat_frames: every calibration frame considered, with the pixel
    --   verdict that admitted or refused it (a file NAMED flat is not
    --   thereby a flat — 25 "Flat Field" frames of 2025-01 are 2 s darks).
    -- s2_flat_pairs:  one row per differenced pair; kind = flat | bias |
    --   dark | sky (sky = star-masked science pairs, for configurations
    --   with no flats on disk).
    -- s2_flat_points: the (signal, variance) points each pair yields.
    -- s2_flat_fits:   every estimator's fit per configuration, chi2 and
    --   dof included (standing rule 1).
    -- s2_camera_configs: the adopted gain / read noise / full scale per
    --   configuration — the table downstream code reads.
    -- ------------------------------------------------------------------
    CREATE TABLE IF NOT EXISTS s2_flat_frames (
        path TEXT PRIMARY KEY, obs_rowid INTEGER, config TEXT, role TEXT,
        night TEXT, filter TEXT, exptime REAL, date_obs TEXT, jd REAL,
        ccd_temp REAL, temp_group INTEGER, group_key TEXT);
    CREATE TABLE IF NOT EXISTS s2_flat_pairs (
        pair_id TEXT PRIMARY KEY, config TEXT, kind TEXT, group_key TEXT,
        filter TEXT, exptime REAL, temp_group INTEGER, night TEXT,
        path_a TEXT, path_b TEXT, gap_s REAL,
        level_a REAL, level_b REAL, ratio REAL, ratio_scatter REAL,
        bias_adu REAL, bias_basis TEXT, n_tiles INTEGER, n_used INTEGER,
        n_points INTEGER, rho_x REAL, rho_x_err REAL, rho_y REAL,
        rho_y_err REAL, status TEXT);
    CREATE TABLE IF NOT EXISTS s2_flat_points (
        pair_id TEXT, config TEXT, kind TEXT, signal_adu REAL,
        var_adu2 REAL, var_err_adu2 REAL, n_tiles INTEGER,
        n_dropped INTEGER, ratio REAL);
    CREATE TABLE IF NOT EXISTS s2_flat_fits (
        config TEXT, estimator TEXT, gain REAL, gain_err REAL,
        gain_err_formal REAL, gain_err_boot REAL, gain_boot_bias REAL,
        intercept REAL, intercept_err REAL, read_noise_adu REAL,
        read_noise_adu_err REAL, chi2 REAL, dof INTEGER, chi2nu REAL,
        quad_coeff REAL, quad_coeff_err REAL, quad_z REAL, gain_quad REAL,
        n_points INTEGER, n_pairs INTEGER, signal_lo REAL, signal_hi REAL,
        PRIMARY KEY (config, estimator));
    -- Per-pair sky gains: K_i = S / (V - V0) for each signal-dominated
    -- sky pair, with the zero reference (bias, V0) it was computed
    -- against and whether it survived the robust clip.
    CREATE TABLE IF NOT EXISTS s2_sky_gains (
        pair_id TEXT PRIMARY KEY, config TEXT, nsub INTEGER,
        signal_adu REAL, var_adu2 REAL, bias_adu REAL, v0_adu2 REAL,
        k_i REAL, used INTEGER, reject_reason TEXT);
    CREATE TABLE IF NOT EXISTS s2_camera_configs (
        config TEXT PRIMARY KEY, camera TEXT, mode TEXT, n_native INTEGER,
        header_egain REAL, n_frames INTEGER, first_night TEXT,
        last_night TEXT, eras TEXT,
        gain_e_per_adu REAL, gain_err REAL, gain_stat_err REAL,
        gain_sys_err REAL, gain_basis TEXT, gain_rel_err REAL,
        meets_target INTEGER,
        read_noise_adu REAL, read_noise_adu_err REAL, read_noise_e REAL,
        read_noise_e_err REAL, rn_basis TEXT, bias_adu REAL,
        clip_adu REAL, full_scale_e REAL, full_scale_e_err REAL,
        binning_verdict TEXT, binning_ratio REAL, binning_ratio_err REAL,
        rho_nn REAL, rho_nn_err REAL, status TEXT, note TEXT);
    -- ------------------------------------------------------------------
    -- Star-based differential linearity (review F-5, DE.F3).
    -- s2_starlin_frames/meas: what was read and measured (one row per
    --   star x frame x aperture) — the raw material, kept so the analysis
    --   in `params` never needs the archive again.
    -- s2_linearity_points: per star x frame fractional deviation from
    --   linear response against that star's own raw peak.
    -- s2_linearity_curve:  the binned deviation curve per source and mode.
    -- s2_linearity_injection: signed recovered-minus-injected bias of the
    --   estimator on each source's real geometry (standing rule 3).
    -- s2_linearity_caps:   ONE recommended cap per mode.
    -- ------------------------------------------------------------------
    CREATE TABLE IF NOT EXISTS s2_starlin_frames (
        set_id TEXT, path TEXT, family TEXT, mode TEXT, config TEXT,
        egain REAL, filter TEXT, exptime REAL, jd REAL, night TEXT,
        is_coord_ref INTEGER, fwhm_med REAL, n_detected INTEGER,
        n_matched INTEGER, shift_x REAL, shift_y REAL, status TEXT,
        PRIMARY KEY (set_id, path));
    CREATE TABLE IF NOT EXISTS s2_starlin_meas (
        set_id TEXT, path TEXT, star_id INTEGER, aperture_px REAL,
        x REAL, y REAL, flux REAL, flux_err REAL, peak_raw REAL, sky REAL,
        fwhm REAL, PRIMARY KEY (set_id, path, star_id, aperture_px));
    CREATE TABLE IF NOT EXISTS s2_linearity_points (
        source TEXT, set_id TEXT, mode TEXT, path TEXT, star_id INTEGER,
        peak_raw REAL, peak_frac REAL, dev REAL, dev_err REAL,
        is_ref INTEGER);
    CREATE TABLE IF NOT EXISTS s2_linearity_curve (
        source TEXT, mode TEXT, lo REAL, hi REAL, peak_frac REAL,
        dev_pct REAL, dev_err_pct REAL, n_points INTEGER,
        n_groups INTEGER, measured INTEGER,
        PRIMARY KEY (source, mode, lo));
    CREATE TABLE IF NOT EXISTS s2_linearity_injection (
        source TEXT, mode TEXT, case_name TEXT, lo REAL, hi REAL,
        injected_pct REAL, recovered_pct REAL, bias_pct REAL,
        bias_err_pct REAL, n_trials INTEGER,
        PRIMARY KEY (source, mode, case_name, lo));
    CREATE TABLE IF NOT EXISTS s2_peaklin_series (
        series_key TEXT PRIMARY KEY, mode TEXT, provenance TEXT,
        era_id INTEGER, filter TEXT, n_frames INTEGER, n_stars INTEGER,
        n_anchored INTEGER, n_points INTEGER, n_points_above_ref INTEGER,
        max_peak_frac REAL, z_err_median REAL, zero_adu REAL,
        scale_to_raw REAL, status TEXT);
    CREATE TABLE IF NOT EXISTS s2_linearity_source_caps (
        mode TEXT, source TEXT, cap_fraction REAL, limited_by TEXT,
        worst_dev_pct REAL, precision_pct REAL, first_bad_lo REAL,
        first_bad_dev_pct REAL, first_bad_err_pct REAL,
        measured_to_fraction REAL, PRIMARY KEY (mode, source));
    CREATE TABLE IF NOT EXISTS s2_linearity_caps (
        mode TEXT PRIMARY KEY, bias_adu REAL, ceiling_adu REAL,
        cap_fraction REAL, cap_adu REAL, limited_by TEXT,
        worst_dev_pct REAL, precision_pct REAL, slope_pct_per_scale REAL,
        slope_err REAL, first_bad_lo REAL, first_bad_dev_pct REAL,
        measured_to_fraction REAL, sources TEXT, status TEXT, note TEXT);
    -- ------------------------------------------------------------------
    -- Peak-at-target saturation census (review OA.E3; standing rule 4).
    -- ------------------------------------------------------------------
    CREATE TABLE IF NOT EXISTS s2_target_peaks (
        obs_rowid INTEGER PRIMARY KEY, path TEXT, target_key TEXT,
        night TEXT, mode TEXT, config TEXT, filter TEXT, exptime REAL,
        is_canonical INTEGER, dup_group INTEGER, n_bin INTEGER,
        locate_method TEXT, x REAL, y REAL, locate_offset_px REAL,
        peak_raw REAL, sky REAL, n_at_peak INTEGER, fwhm_px REAL,
        dispersion_verdict TEXT, status TEXT);
    CREATE TABLE IF NOT EXISTS s2_target_verdicts (
        obs_rowid INTEGER PRIMARY KEY, target_key TEXT, mode TEXT,
        filter TEXT, bias_adu REAL, clip_adu REAL, native_factor REAL,
        peak_native REAL, peak_fraction REAL, cap_fraction REAL,
        verdict TEXT);
    -- ------------------------------------------------------------------
    -- Bad-pixel masks per camera (review DE cross-cutting 5, OA.E8).
    -- ------------------------------------------------------------------
    CREATE TABLE IF NOT EXISTS s2_badpix_masks (
        mask_key TEXT PRIMARY KEY, camera TEXT, naxis1 INTEGER,
        naxis2 INTEGER, temp_group INTEGER, n_native INTEGER,
        n_science INTEGER, n_hot INTEGER, n_rail INTEGER, n_rail1 INTEGER,
        n_rail2 INTEGER, n_rail3 INTEGER, n_noisy INTEGER, n_bad INTEGER,
        bad_fraction REAL, n_dark INTEGER, dark_exptime REAL,
        n_dark_hot INTEGER, n_dark_and_science INTEGER,
        rail1_median_adu REAL, rail1_expected_adu REAL,
        pedestal_adu REAL, npz_path TEXT);
    CREATE TABLE IF NOT EXISTS s2_badpix_frames (
        mask_key TEXT, path TEXT, role TEXT, exptime REAL, sigma REAL,
        threshold REAL, n_fired INTEGER, n_rail INTEGER,
        PRIMARY KEY (mask_key, path));
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
    #
    # ALIAS-SAFE since the 2026-10-03 review (DE.F4): a scene's rows are
    # first collapsed to one row per DATE-OBS, and a pair needs a real time
    # gap of at least ~one exposure.  The v1.2 code paired list neighbours
    # as they came, and for the 2026 camera 7 of 16 "pairs" were a frame
    # against its own ``_wcs`` copy (gap 0 s, variance exactly 0).  The
    # scene-size test also counts DISTINCT exposures now: three rows that
    # are one exposure written three times are not a three-frame scene.
    per_scene: list[tuple[tuple, list]] = []
    for key, members in scenes.items():
        if len(noisemod.distinct_exposures(members)) < NOISE_MIN_SCENE_FRAMES:
            continue
        pairs = noisemod.consecutive_pairs(members, NOISE_MAX_GAP_DAYS,
                                           exptime_s=key[4])
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
    # Purge pairs a previous code version recorded that the alias rule now
    # refuses: any measured pair with no time gap between its frames is a
    # frame differenced against its own copy.  Its points go too (they are
    # exact zeros and would poison the curve's lower envelope).
    stale = [r[0] for r in con.execute(
        "SELECT pair_id FROM s2_noise_pairs WHERE gap_s <= 0")]
    if stale:
        con.executemany("DELETE FROM s2_noise_points WHERE pair_id = ?",
                        [(pid,) for pid in stale])
        con.executemany("DELETE FROM s2_noise_pairs WHERE pair_id = ?",
                        [(pid,) for pid in stale])
        con.commit()
        print(f"[S2:noise] purged {len(stale)} alias pairs (gap 0 s) "
              "recorded by S2 v1.2.")
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
            if np.array_equal(a, b):
                # Belt to the DATE-OBS braces: two files holding the same
                # pixels are one exposure whatever their headers claim.
                # Recorded as a sentinel (n_points 0) so it is never
                # retried and never counted as a measurement.
                print(f"[S2:noise]   SKIP {pr['pair_id']}: identical pixels "
                      "(alias of one exposure)")
                con.execute("DELETE FROM s2_noise_points WHERE pair_id = ?",
                            (pr["pair_id"],))
                con.execute("INSERT OR REPLACE INTO s2_noise_pairs VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?,0)",
                            (pr["pair_id"], mode, pr["night"],
                             pr["target_key"], pr["filter"], pr["egain"],
                             pr["exptime"], pr["era_id"], pr["path_a"],
                             pr["path_b"], pr["gap_s"]))
                con.commit()
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
# Subcommand: flats  (flat-pair photon transfer per camera/configuration)
# ---------------------------------------------------------------------------
def calibration_role(imagetyp: str | None, basename: str) -> str | None:
    """What a calibration-tree frame claims to be: flat, bias, dark, None.

    Masters (already averaged — their variance is not a single frame's)
    and ``_cal`` files (pipeline-reduced copies) are refused outright.
    IMAGETYP decides bias and dark; "flat" is accepted from IMAGETYP or
    from the filename, because the 2025 panel flats were saved as
    ``Light Frame`` under names like ``Flat_g_07.fts``.  The claim is only
    a claim: :func:`cmd_flats` verifies every flat in pixels.
    """
    b = basename.lower()
    if b.startswith("master") or "_cal" in b:
        return None
    it = (imagetyp or "").strip().lower()
    if it.startswith("bias"):
        return "bias"
    if it.startswith("dark"):
        return "dark"
    if it.startswith("flat") or "flat" in b:
        return "flat"
    return None


def temp_group_of(ccd_temp) -> int | None:
    """CCD-TEMP rounded to the nearest :data:`TEMP_GROUP_STEP_C` degrees."""
    if ccd_temp is None:
        return None
    return int(round(float(ccd_temp) / TEMP_GROUP_STEP_C) * TEMP_GROUP_STEP_C)


def _spread(items: list, n: int) -> list:
    """At most ``n`` items, evenly spread through the sequence."""
    if len(items) <= n:
        return list(items)
    idx = np.unique(np.linspace(0, len(items) - 1, n).round().astype(int))
    return [items[i] for i in idx]


def _flat_census(con: sqlite3.Connection) -> list[dict]:
    """Every calibration frame the flat-PTC stage may use, grouped.

    One query against the manifest, then pure classification: role from
    :func:`calibration_role`, configuration from
    :func:`rlmt_diagnostics.camera.config_key`, and a group key
    (configuration | directory | filter | exposure | temperature group)
    inside which consecutive frames are pairable.  Recorded in
    ``s2_flat_frames`` so the report can say exactly which frames each
    configuration's gain rests on.
    """
    trees = ", ".join(f"'{t}'" for t in FLAT_TREES)
    rows = con.execute(f"""
        SELECT f.obs_rowid, f.path, f.basename, f.tree, f.imagetyp,
               f.instrume, f.readoutm, f.xbinning, f.naxis1, f.naxis2,
               f.egain, f.filter, f.exptime, f.date_obs, f.jd, f.night,
               f.ccd_temp
        FROM frames f
        WHERE f.is_canonical = 1 AND f.tree IN ({trees})
          AND f.jd IS NOT NULL AND f.naxis1 >= 1000 AND f.naxis2 >= 1000
          AND (f.imagetyp LIKE 'Flat%' OR f.imagetyp LIKE 'Bias%'
               OR f.imagetyp LIKE 'Dark%' OR lower(f.basename) LIKE '%flat%')
        ORDER BY f.jd, f.path""").fetchall()
    out = []
    for (rowid, path, base, tree, ityp, instr, mode, xb, n1, n2, eg, filt,
         expt, dobs, jd, night, temp) in rows:
        role = calibration_role(ityp, base)
        if role is None:
            continue
        if role == "dark" and tree not in ZERO_DARK_TREES:
            continue
        config = cam.config_key(instr, mode, xb, n1, eg)
        tg = temp_group_of(temp)
        ebin = exptime_bin(expt)
        directory = path.rsplit("/", 1)[0]
        # Zero-signal frames group without filter/directory: a bias is a
        # bias whatever wheel slot was in the beam.
        gkey = (f"{config}|{role}|{ebin}|T{tg}|{int(n1)}x{int(n2)}"
                if role != "flat" else
                f"{config}|flat|{directory}|{filt or ''}|{ebin}|T{tg}"
                f"|{int(n1)}x{int(n2)}")
        out.append({"obs_rowid": rowid, "path": path, "config": config,
                    "role": role, "night": night, "filter": filt or "",
                    "exptime": ebin, "date_obs": dobs, "jd": jd,
                    "ccd_temp": temp, "temp_group": tg, "group_key": gkey})
    return out


def _zero_groups(census: list[dict]) -> dict[tuple, list[dict]]:
    """Per (config, temp group): the frames that serve as zero-signal pairs.

    Bias frames when the configuration has any; otherwise its darks at the
    SHORTEST exposure on disk (the dark-current shot term in those is
    measured and carried by the existing ``ptc`` stage; here they are the
    read-noise anchor and the bias level for the flats).
    """
    by: dict[tuple, dict[str, list[dict]]] = {}
    for fr in census:
        if fr["role"] in ("bias", "dark"):
            by.setdefault((fr["config"], fr["temp_group"]), {}) \
              .setdefault(fr["role"], []).append(fr)
    out: dict[tuple, list[dict]] = {}
    for key, roles in by.items():
        if roles.get("bias") and len(roles["bias"]) >= 4:
            out[key] = roles["bias"]
        elif roles.get("dark"):
            t_min = min(d["exptime"] for d in roles["dark"]
                        if d["exptime"] is not None)
            out[key] = [d for d in roles["dark"] if d["exptime"] == t_min]
    return out


def _central_level(img: np.ndarray) -> float:
    """Median of the central quarter of a frame (the level a pair is at)."""
    ny, nx = img.shape
    return float(np.median(img[ny // 4: 3 * ny // 4, nx // 4: 3 * nx // 4]))


def _store_pair(con, pair: dict, points: list[dict]) -> None:
    """Write one measured pair and its points (replace, never append)."""
    con.execute("DELETE FROM s2_flat_points WHERE pair_id = ?",
                (pair["pair_id"],))
    con.executemany(
        "INSERT INTO s2_flat_points VALUES (?,?,?,?,?,?,?,?,?)",
        [(pair["pair_id"], pair["config"], pair["kind"], pt["signal"],
          pt["var"], pt["var_err"], pt["n_tiles"], pt["n_dropped"],
          pt["ratio"]) for pt in points])
    cols = ("pair_id", "config", "kind", "group_key", "filter", "exptime",
            "temp_group", "night", "path_a", "path_b", "gap_s", "level_a",
            "level_b", "ratio", "ratio_scatter", "bias_adu", "bias_basis",
            "n_tiles", "n_used", "n_points", "rho_x", "rho_x_err", "rho_y",
            "rho_y_err", "status")
    con.execute(
        f"INSERT OR REPLACE INTO s2_flat_pairs ({', '.join(cols)}) "
        f"VALUES ({', '.join('?' * len(cols))})",
        tuple(pair.get(c) for c in cols))
    con.commit()


def _config_with_header(config: str, header) -> str:
    """Refine a configuration key with what only the header knows.

    The iKon's pre-amplifier setting (GAIN card) is not in the manifest
    but changes the gain by its factor; pairs are filed under the key
    with the preamp appended (``iKon 1MHz 4x``).
    """
    if config.startswith("iKon") and len(config.split()) == 2:
        g = str(header.get("GAIN", "") or "").strip()
        if g:
            return f"{config} {g}"
    return config


def _measure_signal_pair(con, archive: Path, pair: dict, bias: float | None,
                         clip: float | None, sky: bool) -> None:
    """Read one flat (or sky) pair, verify it in pixels, store its points.

    The pixel checks, in order — each failure is a recorded status, never
    a silent skip: unreadable; shape mismatch; identical pixels (an alias
    the DATE-OBS rule missed); no signal above the bias (a "flat" that is
    a dark); level within 15% of the mode's ceiling (variance collapses
    at the clip); and an illumination pattern that changed between the two
    frames (tile ratios scatter).
    """
    try:
        a, hdr_a = read_image(archive, pair["path_a"])
        b, _ = read_image(archive, pair["path_b"])
    except Exception as e:
        print(f"[S2:{pair['kind']}]   SKIP {pair['pair_id']}: {e}")
        _store_pair(con, {**pair, "status": "unreadable"}, [])
        return
    pair["config"] = _config_with_header(pair["config"], hdr_a)
    if a.shape != b.shape:
        _store_pair(con, {**pair, "status": "shape_mismatch"}, [])
        return
    af, bf = a.astype(np.float64), b.astype(np.float64)
    pair["level_a"], pair["level_b"] = _central_level(af), _central_level(bf)
    if np.array_equal(a, b):
        _store_pair(con, {**pair, "status": "identical_pixels"}, [])
        return
    b0 = bias if bias is not None else 0.0
    if clip and max(pair["level_a"], pair["level_b"]) \
            > ptc.PTC_LEVEL_CEILING_FRACTION * clip:
        _store_pair(con, {**pair, "status": "near_ceiling"}, [])
        return
    valid = None
    if sky:
        valid = flatptc.star_free_mask((af + bf) / 2.0)
    st = flatptc.tile_pair_stats(af, bf, bias=b0, valid=valid)
    pair["n_tiles"], pair["n_used"] = st["n_tiles"], st["n_used"]
    if st["n_used"]:
        r = st["ratio"]
        pair["ratio"] = float(np.median(r))
        pair["ratio_scatter"] = float(
            1.4826 * np.median(np.abs(r - np.median(r))) / np.median(r))
    points = flatptc.pair_points(st)
    status = "ok"
    if not points:
        status = "no_usable_tiles"
    elif (not sky and pair.get("zero_sigma")
          and np.median(st["signal"])
          < FLAT_MIN_SIGNAL_SIGMAS * pair["zero_sigma"]):
        status = "no_signal"
    elif (not sky and pair.get("ratio_scatter") is not None
          and pair["ratio_scatter"] > FLAT_MAX_RATIO_SCATTER):
        status = "scene_changed"
    if status == "ok" and not sky:
        # Pixel-correlation diagnostic on the central megapixel (the
        # statistic is a per-tile mean; the centre is where a flat is
        # flattest and the cost stays bounded).
        ny, nx = af.shape
        cy, cx, h = ny // 2, nx // 2, 512
        ac = flatptc.diff_autocorrelation(
            af[cy - h:cy + h, cx - h:cx + h], bf[cy - h:cy + h, cx - h:cx + h],
            bias=b0, edge_fraction=0.0)
        pair.update({"rho_x": ac["rho_x"], "rho_x_err": ac["rho_x_err"],
                     "rho_y": ac["rho_y"], "rho_y_err": ac["rho_y_err"]})
    pair["n_points"] = len(points) if status == "ok" else 0
    _store_pair(con, {**pair, "status": status},
                points if status == "ok" else [])


def cmd_flats(con: sqlite3.Connection, archive: Path, batch: int,
              only: str | None = None) -> int:
    """Flat-pair photon transfer for every camera/configuration with flats.

    Three passes, all resumable (a pair recorded in ``s2_flat_pairs`` is
    never re-read):

    1. the census (``s2_flat_frames``), rebuilt from the manifest each
       call — it is one query and no pixels;
    2. zero-signal pairs (bias, or shortest darks) per configuration and
       temperature group: the read noise and the bias level;
    3. flat pairs, differenced with the level ratio and the bias from
       pass 2 (:func:`rlmt_diagnostics.flatptc.tile_pair_stats`).

    Fitting is not done here: ``params`` distils ``s2_flat_points`` into
    ``s2_flat_fits`` and ``s2_camera_configs``.

    ``only`` restricts the pixel passes to configurations whose key
    contains that text (the census is always complete) — a scheduling
    aid for a slow disk, never a change to what is eventually measured.
    """
    census = _flat_census(con)
    con.execute("DELETE FROM s2_flat_frames")
    con.executemany(
        "INSERT OR REPLACE INTO s2_flat_frames VALUES "
        "(?,?,?,?,?,?,?,?,?,?,?,?)",
        [(f["path"], f["obs_rowid"], f["config"], f["role"], f["night"],
          f["filter"], f["exptime"], f["date_obs"], f["jd"], f["ccd_temp"],
          f["temp_group"], f["group_key"]) for f in census])
    con.commit()
    done = {r[0] for r in con.execute("SELECT pair_id FROM s2_flat_pairs")}
    ceilings = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes "
        "WHERE clip_adu IS NOT NULL"))
    n_done = 0

    # ---- pass 2: zero-signal pairs --------------------------------------
    for (config, tg), frames_ in sorted(
            _zero_groups(census).items(), key=lambda kv: str(kv[0])):
        if only and only not in config:
            continue
        by_group: dict[str, list] = {}
        for fr in frames_:
            by_group.setdefault(fr["group_key"], []).append(
                (fr["jd"], fr["obs_rowid"], fr["path"]))
        cands = []
        for gkey, members in sorted(by_group.items()):
            kind = gkey.split("|")[1]
            for (ja, ia, pa), (jb, ib, pb) in noisemod.consecutive_pairs(
                    members, FLAT_MAX_GAP_DAYS):
                cands.append((gkey, kind, ia, ib, pa, pb, (jb - ja) * 86400.0))
        for gkey, kind, ia, ib, pa, pb, gap in _spread(cands, ZERO_MAX_PAIRS):
            pair_id = f"z{ia}-{ib}"
            if pair_id in done:
                continue
            if n_done >= batch:
                print("[S2:flats] batch cap reached; re-run to continue.")
                return 0
            fr0 = next(f for f in frames_ if f["path"] == pa)
            pair = {"pair_id": pair_id, "config": config, "kind": kind,
                    "group_key": gkey, "filter": fr0["filter"],
                    "exptime": fr0["exptime"], "temp_group": tg,
                    "night": fr0["night"], "path_a": pa, "path_b": pb,
                    "gap_s": gap}
            try:
                a, hdr_a = read_image(archive, pa)
                b, _ = read_image(archive, pb)
            except Exception as e:
                print(f"[S2:flats]   SKIP {pair_id}: {e}")
                _store_pair(con, {**pair, "status": "unreadable"}, [])
                continue
            pair["config"] = _config_with_header(config, hdr_a)
            n_done += 1
            if a.shape != b.shape or np.array_equal(a, b):
                _store_pair(con, {**pair, "status": "shape_mismatch"
                                  if a.shape != b.shape
                                  else "identical_pixels"}, [])
                continue
            rn = flatptc.pair_read_noise(a.astype(np.float64),
                                         b.astype(np.float64))
            if rn is None:
                _store_pair(con, {**pair, "status": "no_usable_tiles"}, [])
                continue
            pair.update({"level_a": _central_level(a),
                         "level_b": _central_level(b), "ratio": 1.0,
                         "n_tiles": rn["n_tiles"], "n_used": rn["n_tiles"],
                         "n_points": 1, "status": "ok"})
            _store_pair(con, pair, [{
                "signal": 0.0, "var": rn["var"], "var_err": rn["var_err"],
                "n_tiles": rn["n_tiles"], "n_dropped": 0, "ratio": 1.0}])
            print(f"[S2:flats]   {config} T{tg} {kind}: RN "
                  f"{rn['read_noise_adu']:.3f} ADU at level "
                  f"{rn['level']:.1f}")

    # ---- the bias level and zero-level sigma each flat pair will use ----
    zero: dict[tuple, tuple[float, float, str]] = {}
    for config, tg, lvl, var, n in con.execute("""
            SELECT p.config, p.temp_group,
                   avg((p.level_a + p.level_b) / 2.0), avg(q.var_adu2),
                   count(*)
            FROM s2_flat_pairs p JOIN s2_flat_points q USING (pair_id)
            WHERE p.kind IN ('bias', 'dark') AND p.status = 'ok'
            GROUP BY p.config, p.temp_group"""):
        zero[(config, tg)] = (lvl, float(np.sqrt(max(var, 0.0))),
                              f"{n} zero-signal pairs, T group {tg}")

    def zero_for(config: str, tg) -> tuple[float | None, float | None, str]:
        if (config, tg) in zero:
            return zero[(config, tg)]
        # Header-refined keys ('iKon 1MHz 4x') serve their base key.
        other = [v for (c, t), v in zero.items()
                 if c.startswith(config) and t == tg]
        other += [v for (c, _t), v in zero.items() if c.startswith(config)]
        if other:
            return other[0][0], other[0][1], other[0][2] + " (other T group)"
        return None, None, "none (free-intercept only)"

    # ---- pass 3: flat pairs ----------------------------------------------
    groups: dict[str, list[dict]] = {}
    for fr in census:
        if fr["role"] == "flat":
            groups.setdefault(fr["group_key"], []).append(fr)
    for gkey, frames_ in sorted(groups.items()):
        if only and only not in frames_[0]["config"]:
            continue
        members = [(f["jd"], f["obs_rowid"], f["path"]) for f in frames_]
        pairs = noisemod.consecutive_pairs(members, FLAT_MAX_GAP_DAYS,
                                           exptime_s=frames_[0]["exptime"])
        fr0 = frames_[0]
        bias, zsig, basis = zero_for(fr0["config"], fr0["temp_group"])
        clip = ceilings.get(cam.mode_label_of(fr0["config"]))
        for (ja, ia, pa), (jb, ib, pb) in _spread(pairs,
                                                  FLAT_MAX_PAIRS_PER_GROUP):
            pair_id = f"f{ia}-{ib}"
            if pair_id in done:
                continue
            if n_done >= batch:
                print("[S2:flats] batch cap reached; re-run to continue.")
                return 0
            pair = {"pair_id": pair_id, "config": fr0["config"],
                    "kind": "flat", "group_key": gkey,
                    "filter": fr0["filter"], "exptime": fr0["exptime"],
                    "temp_group": fr0["temp_group"], "night": fr0["night"],
                    "path_a": pa, "path_b": pb, "gap_s": (jb - ja) * 86400.0,
                    "bias_adu": bias, "bias_basis": basis,
                    "zero_sigma": zsig}
            _measure_signal_pair(con, archive, pair, bias, clip, sky=False)
            n_done += 1
            row = con.execute(
                "SELECT status, level_a, ratio FROM s2_flat_pairs "
                "WHERE pair_id = ?", (pair_id,)).fetchone()
            print(f"[S2:flats]   {fr0['config']} {fr0['filter']}: {row[0]} "
                  f"(level {row[1] or 0:,.0f}, r {row[2] or 0:.3f})")
    print(f"[S2:flats] nothing left to do. ({n_done} new pairs this run)")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: skypairs  (sky-pair photon transfer, every configuration)
# ---------------------------------------------------------------------------
def _sky_candidates(con: sqlite3.Connection) -> dict[str, list[dict]]:
    """Per configuration: sky pairs stratified on exposure time.

    A sky pair is two DISTINCT consecutive science exposures of one scene
    (alias-safe: :func:`rlmt_diagnostics.noise.consecutive_pairs`).  The
    photon-transfer signal is the SKY — stars are masked — so the level
    axis is spanned by choosing scenes across the configuration's whole
    range of exposure times: the scenes are ordered by exposure and
    :data:`SKY_MAX_PAIRS_PER_CONFIG` of them taken evenly through that
    order, one pair per scene.
    """
    rows = con.execute(f"""
        SELECT f.obs_rowid, f.path, f.basename, f.instrume, f.readoutm,
               f.xbinning, f.naxis1, f.naxis2, f.egain, f.filter, f.exptime,
               f.jd, f.night, f.target_key, f.ccd_temp
        FROM frames f
        WHERE f.is_canonical = 1 AND f.tree = 'rawimage'
          AND {SCIENCE_IMAGETYP} AND f.target_key IS NOT NULL
          AND f.exptime > 0 AND f.jd IS NOT NULL
          AND f.naxis1 >= 1000 AND f.naxis2 >= 1000
        ORDER BY f.jd, f.path""").fetchall()
    from collections import defaultdict
    scenes: dict[tuple, list] = defaultdict(list)
    n_config: dict[str, int] = defaultdict(int)
    for (rowid, path, base, instr, mode, xb, n1, n2, eg, filt, expt, jd,
         night, tkey, temp) in rows:
        if (filt or "").strip().lower() in SKY_EXCLUDED_FILTERS:
            continue
        ebin = exptime_bin(expt)
        if ebin is None or ebin <= 0:
            continue
        config = cam.config_key(instr, mode, xb, n1, eg)
        n_config[config] += 1
        scenes[(config, night, tkey, filt or "", ebin, int(n1), int(n2),
                temp_group_of(temp))].append((jd, rowid, path))
    per_config: dict[str, list[dict]] = defaultdict(list)
    for key, members in scenes.items():
        config, night, tkey, filt, ebin, n1, n2, tg = key
        if (n_config[config] < SKY_MIN_CONFIG_FRAMES
                and config not in SKY_ALWAYS_CONFIGS):
            continue
        pairs = noisemod.consecutive_pairs(members, NOISE_MAX_GAP_DAYS,
                                           exptime_s=ebin)
        if not pairs:
            continue
        (ja, ia, pa), (jb, ib, pb) = pairs[0]
        per_config[config].append({
            "pair_id": f"s{ia}-{ib}", "config": config, "kind": "sky",
            "group_key": f"{config}|sky|{night}|{tkey}|{filt}|{ebin}",
            "filter": filt, "exptime": ebin, "temp_group": tg,
            "night": night, "path_a": pa, "path_b": pb,
            "gap_s": (jb - ja) * 86400.0})
    out = {}
    for config, cands in per_config.items():
        cands.sort(key=lambda c: (c["exptime"], c["night"], c["pair_id"]))
        out[config] = _spread(cands, SKY_MAX_PAIRS_PER_CONFIG)
    return out


def cmd_skypairs(con: sqlite3.Connection, archive: Path, batch: int,
                 only: str | None = None) -> int:
    """Sky-pair photon transfer: the gain for configurations without flats.

    Same tile arithmetic as ``flats``, applied to star-masked consecutive
    science frames.  It is run for EVERY configuration, including the ones
    that do have flats — those are the validation set: ``params`` compares
    the sky-pair gain with the flat-pair gain camera by camera, and only
    that measured agreement licenses the sky-pair number for a
    configuration (StackPro, the 2023 High Gain epoch, the QHY600) whose
    flats do not exist.
    """
    done = {r[0] for r in con.execute("SELECT pair_id FROM s2_flat_pairs")}
    ceilings = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes "
        "WHERE clip_adu IS NOT NULL"))
    zero = {c: lvl for c, lvl in con.execute("""
        SELECT config, avg((level_a + level_b) / 2.0) FROM s2_flat_pairs
        WHERE kind IN ('bias', 'dark') AND status = 'ok' GROUP BY config""")}
    # A header-refined zero key ('iKon 1MHz 4x') serves its base key too.
    for c in list(zero):
        zero.setdefault(" ".join(c.split()[:2]), zero[c])
    n_done = 0
    for config, cands in sorted(_sky_candidates(con).items()):
        if only and only not in config:
            continue
        clip = ceilings.get(cam.mode_label_of(config))
        bias = zero.get(config)
        for pair in cands:
            if pair["pair_id"] in done:
                continue
            if n_done >= batch:
                print("[S2:skypairs] batch cap reached; re-run to continue.")
                return 0
            pair.update({"bias_adu": bias, "bias_basis": (
                "zero-signal pairs" if bias is not None
                else "none (raw level; slope is bias-independent)")})
            _measure_signal_pair(con, archive, pair, bias, clip, sky=True)
            n_done += 1
            row = con.execute(
                "SELECT status, level_a, n_points FROM s2_flat_pairs "
                "WHERE pair_id = ?", (pair["pair_id"],)).fetchone()
            print(f"[S2:skypairs]   {config} {pair['exptime']:g}s "
                  f"{pair['filter']}: {row[0]} (level {row[1] or 0:,.0f}, "
                  f"{row[2] or 0} points)")
    print(f"[S2:skypairs] nothing left to do. ({n_done} new pairs this run)")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: starlin  (star photometry of the dedicated linearity sets)
# ---------------------------------------------------------------------------
def _starlin_sets(con: sqlite3.Connection) -> dict[str, list[dict]]:
    """The dedicated linearity sets: {set_id: [frame dicts, time-ordered]}.

    A set is one family x night x filter (one field through one filter).
    Same-DATE-OBS aliases are collapsed, and a set must offer a contrast
    to be worth reading: at least two distinct (mode, exposure) settings,
    or at least three frames of one setting (seeing then supplies the
    peak contrast).
    """
    from collections import defaultdict
    sets: dict[str, list[dict]] = defaultdict(list)
    for family, where in STARLIN_FAMILIES.items():
        rows = con.execute(f"""
            SELECT f.obs_rowid, f.path, f.instrume, f.readoutm, f.xbinning,
                   f.naxis1, f.egain, f.filter, f.exptime, f.jd, f.night
            FROM frames f WHERE {where} AND f.jd IS NOT NULL
            ORDER BY f.jd, f.path""").fetchall()
        for (rowid, path, instr, mode, xb, n1, eg, filt, expt, jd,
             night) in rows:
            sets[f"{family}|{night}|{filt or ''}"].append({
                "obs_rowid": rowid, "path": path, "family": family,
                "mode": ceil.mode_group(mode),
                "config": cam.config_key(instr, mode, xb, n1, eg),
                "egain": eg, "filter": filt or "",
                "exptime": exptime_bin(expt),
                "jd": jd, "night": night})
    out = {}
    for set_id, frames_ in sets.items():
        keep = {m[2] for m in noisemod.distinct_exposures(
            [(f["jd"], f["obs_rowid"], f["path"]) for f in frames_])}
        frames_ = [f for f in frames_ if f["path"] in keep]
        settings = {(f["mode"], f["exptime"]) for f in frames_}
        if len(settings) == 1:
            frames_ = _spread(frames_, STARLIN_MAX_FRAMES_FIXED)
        if len(settings) >= 2 or len(frames_) >= 3:
            out[set_id] = frames_
    return out


def cmd_starlin(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Measure every star in every frame of the dedicated linearity sets.

    Pixels only — no linearity arithmetic happens here.  Per set: the
    coordinate-reference frame (the longest single-read exposure) defines
    the star list; every frame is then measured at all of
    :data:`rlmt_diagnostics.starphot.APERTURES_PX` and matched to that
    list by a voted translation.  ``params`` turns the stored
    measurements into deviation curves, so a change to the estimator
    never costs a re-read of the archive.
    """
    done = {(r[0], r[1]) for r in con.execute(
        "SELECT set_id, path FROM s2_starlin_frames")}
    n_done = 0
    for set_id, frames_ in sorted(_starlin_sets(con).items()):
        single = [f for f in frames_ if "StackPro" not in f["mode"]] or frames_
        ref = max(single, key=lambda f: (f["exptime"], -f["jd"]))
        todo = [f for f in frames_ if (set_id, f["path"]) not in done]
        if not todo:
            continue
        # The star list: from the stored reference measurements when the
        # reference frame is already done, else by measuring it now.
        ref_rows = con.execute(
            "SELECT star_id, x, y FROM s2_starlin_meas WHERE set_id = ? "
            "AND path = ? AND aperture_px = ? ORDER BY star_id",
            (set_id, ref["path"], starphot.APERTURES_PX[0])).fetchall()
        if (set_id, ref["path"]) not in done or not ref_rows:
            todo = [ref] + [f for f in todo if f["path"] != ref["path"]]
            ref_xy = None
        else:
            ref_xy = np.array([(r[1], r[2]) for r in ref_rows])
        print(f"[S2:starlin] {set_id}: {len(todo)} frames to measure "
              f"(reference {ref['path'].rsplit('/', 1)[-1]})")
        for fr in todo:
            if n_done >= batch:
                print("[S2:starlin] batch cap reached; re-run to continue.")
                return 0
            is_ref = fr["path"] == ref["path"]
            base_row = (set_id, fr["path"], fr["family"], fr["mode"],
                        fr["config"], fr["egain"], fr["filter"],
                        fr["exptime"], fr["jd"], fr["night"], int(is_ref))
            cols = ("set_id, path, family, mode, config, egain, filter, "
                    "exptime, jd, night, is_coord_ref, fwhm_med, n_detected, "
                    "n_matched, shift_x, shift_y, status")
            try:
                img, _ = read_image(archive, fr["path"])
            except Exception as e:
                print(f"[S2:starlin]   SKIP {fr['path']}: {e}")
                con.execute(f"INSERT OR REPLACE INTO s2_starlin_frames "
                            f"({cols}) VALUES ({','.join('?' * 17)})",
                            base_row + (None, 0, 0, None, None, "unreadable"))
                con.commit()
                continue
            n_done += 1
            m = starphot.measure_stars(img)
            xy = np.c_[m["x"], m["y"]]
            if is_ref or ref_xy is None:
                ref_xy = xy[:STARLIN_MAX_STARS]
                idx = np.arange(len(ref_xy))
                shift = (0.0, 0.0)
            else:
                shift = starphot.match_shift(ref_xy, xy)
                if shift is None:
                    con.execute(
                        f"INSERT OR REPLACE INTO s2_starlin_frames ({cols}) "
                        f"VALUES ({','.join('?' * 17)})",
                        base_row + (m["fwhm_med"], len(xy), 0, None, None,
                                    "no_match"))
                    con.commit()
                    continue
                idx = starphot.match_stars(ref_xy, xy, shift)
            rows = []
            for star_id, j in enumerate(idx):
                if j < 0:
                    continue
                for k, r_ap in enumerate(m["apertures"]):
                    rows.append((set_id, fr["path"], star_id, r_ap,
                                 float(m["x"][j]), float(m["y"][j]),
                                 float(m["flux"][k][j]),
                                 float(m["flux_err"][k][j]),
                                 float(m["peak_raw"][j]), float(m["sky"][j]),
                                 float(m["fwhm"][j])))
            con.execute("DELETE FROM s2_starlin_meas WHERE set_id = ? "
                        "AND path = ?", (set_id, fr["path"]))
            con.executemany("INSERT INTO s2_starlin_meas VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?)", rows)
            n_matched = int((np.asarray(idx) >= 0).sum())
            con.execute(
                f"INSERT OR REPLACE INTO s2_starlin_frames ({cols}) "
                f"VALUES ({','.join('?' * 17)})",
                base_row + (m["fwhm_med"], len(xy), n_matched,
                            shift[0], shift[1], "ok"))
            con.commit()
            print(f"[S2:starlin]   {fr['mode']} {fr['exptime']:g}s: "
                  f"{n_matched} stars matched (FWHM "
                  f"{m['fwhm_med'] or float('nan'):.1f} px)")
    print(f"[S2:starlin] nothing left to do. ({n_done} frames this run)")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: peakcensus  (peak-at-target saturation census)
# ---------------------------------------------------------------------------
def _target_xy_from_wcs(header, ra_deg: float, dec_deg: float,
                        shape: tuple[int, int]) -> tuple[float, float] | None:
    """Pixel position of a sky coordinate from a celestial WCS, or None.

    None when the header carries no celestial solution, when astropy
    cannot build one, or when the position falls outside the frame (a
    stale or wrong solution must not be trusted to point at empty sky).
    """
    import warnings
    from astropy.wcs import WCS
    if header.get("CRVAL1") is None or header.get("CRPIX1") is None:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            w = WCS(header, naxis=2)
            if not w.has_celestial:
                return None
            x, y = w.all_world2pix([[ra_deg, dec_deg]], 0)[0]
    except Exception:
        return None
    ny, nx = shape
    if not (np.isfinite(x) and np.isfinite(y) and 0 <= x < nx and 0 <= y < ny):
        return None
    return float(x), float(y)


def cmd_peakcensus(con: sqlite3.Connection, archive: Path, batch: int) -> int:
    """Measure the target's own peak pixel in every imaging frame.

    For each target in :data:`CENSUS_TARGETS`: every ``rawimage`` science
    row with a non-grism filter name — canonical or not, because the
    committee's count (224 for T CrB) is of rows, and a duplicate row that
    silently differed from its twin would be a finding in itself.  The
    target is located by, in order of trust: the frame's own header WCS;
    the S1 plate solution on disk; the brightest source near the frame
    centre.  The method used is recorded per frame.  Verdicts are issued
    by ``params`` (they depend on the linearity cap).
    """
    import warnings
    done = {r[0] for r in con.execute("SELECT obs_rowid FROM s2_target_peaks")}
    has_disp = con.execute(
        "SELECT count(*) FROM sqlite_master WHERE name = 'frame_dispersion'"
    ).fetchone()[0] or (has_manifest_schema(con) and con.execute(
        f"SELECT count(*) FROM {MANIFEST_SCHEMA}.sqlite_master "
        "WHERE name = 'frame_dispersion'").fetchone()[0])
    grism = ", ".join(f"'{g}'" for g in CENSUS_GRISM_FILTERS)
    wcs_root = REPO_ROOT / "products" / "astrom" / "wcs"
    n_done = 0
    for target_key, (ra, dec) in CENSUS_TARGETS.items():
        rows = con.execute(f"""
            SELECT f.obs_rowid, f.path, f.night, f.instrume, f.readoutm,
                   f.xbinning, f.naxis1, f.egain, f.filter, f.exptime,
                   f.is_canonical, f.dup_group
            FROM frames f
            WHERE f.tree = 'rawimage' AND f.target_key = ?
              AND {SCIENCE_IMAGETYP}
              AND lower(coalesce(f.filter, '')) NOT IN ({grism})
            ORDER BY f.is_canonical DESC, f.jd, f.path""",
            (target_key,)).fetchall()
        print(f"[S2:peakcensus] {target_key}: {len(rows)} imaging rows, "
              f"{sum(1 for r in rows if r[0] not in done)} to measure.")
        for (rowid, path, night, instr, mode, xb, n1, eg, filt, expt,
             is_canon, dup) in rows:
            if rowid in done:
                continue
            if n_done >= batch:
                print("[S2:peakcensus] batch cap reached; re-run to continue.")
                return 0
            verdict = None
            if has_disp:
                r = con.execute("SELECT verdict FROM frame_dispersion "
                                "WHERE obs_rowid = ?", (rowid,)).fetchone()
                verdict = r[0] if r else None
            rec_ = {"obs_rowid": rowid, "path": path,
                    "target_key": target_key, "night": night,
                    "mode": ceil.mode_group(mode),
                    "config": cam.config_key(instr, mode, xb, n1, eg),
                    "filter": filt or "", "exptime": exptime_bin(expt),
                    "is_canonical": is_canon, "dup_group": dup,
                    "n_bin": int(xb or 1), "dispersion_verdict": verdict}
            cols = ("obs_rowid", "path", "target_key", "night", "mode",
                    "config", "filter", "exptime", "is_canonical",
                    "dup_group", "n_bin", "locate_method", "x", "y",
                    "locate_offset_px", "peak_raw", "sky", "n_at_peak",
                    "fwhm_px", "dispersion_verdict", "status")

            def store(**kw):
                rec_.update(kw)
                con.execute(
                    f"INSERT OR REPLACE INTO s2_target_peaks "
                    f"({', '.join(cols)}) VALUES "
                    f"({', '.join('?' * len(cols))})",
                    tuple(rec_.get(c) for c in cols))
                con.commit()

            try:
                img, hdr = read_image(archive, path)
            except Exception as e:
                print(f"[S2:peakcensus]   SKIP {path}: {e}")
                store(status="unreadable")
                continue
            n_done += 1
            data = np.asarray(img, dtype=np.float64)
            ny, nx = data.shape
            xy, method = _target_xy_from_wcs(hdr, ra, dec, data.shape), \
                "wcs_header"
            if xy is None:
                wcs_file = wcs_root / (path[:-3] if path.endswith(".fz")
                                       else path)
                wcs_file = wcs_file.with_suffix(".wcs")
                if wcs_file.exists():
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            whdr = fits.getheader(wcs_file)
                        xy = _target_xy_from_wcs(whdr, ra, dec, data.shape)
                        method = "wcs_s1"
                    except Exception:
                        xy = None
            if xy is None:
                m = starphot.measure_stars(data)
                j = satmod.brightest_central_source(
                    m["x"], m["y"], m["flux"][0] if m["flux"].size else [],
                    data.shape)
                if j is None:
                    store(status="target_not_found", locate_method="none")
                    continue
                xy, method = (float(m["x"][j]), float(m["y"][j])), \
                    "brightest_central"
            # Re-centre on the star: the brightest 3x3 BOX-SUM within the
            # search radius (a single hot pixel cannot out-sum a star).
            R = satmod.RECENTRE_RADIUS_PX
            x0, y0 = int(round(xy[0])), int(round(xy[1]))
            ya, yb = max(y0 - R - 1, 0), min(y0 + R + 2, ny)
            xa, xb_ = max(x0 - R - 1, 0), min(x0 + R + 2, nx)
            box = data[ya:yb, xa:xb_]
            if box.shape[0] < 3 or box.shape[1] < 3:
                store(status="target_off_frame", locate_method=method)
                continue
            bs = sum(box[dy:box.shape[0] - 2 + dy, dx:box.shape[1] - 2 + dx]
                     for dy in range(3) for dx in range(3))
            by, bx = np.unravel_index(int(np.argmax(bs)), bs.shape)
            cy, cx = ya + by + 1, xa + bx + 1
            H = satmod.STAMP_HALF
            stamp = data[max(cy - H, 0):cy + H + 1, max(cx - H, 0):cx + H + 1]
            st = satmod.stamp_peak_stats(stamp)
            store(status="ok", locate_method=method, x=float(cx), y=float(cy),
                  locate_offset_px=float(np.hypot(cx - xy[0], cy - xy[1])),
                  peak_raw=st["peak_raw"], sky=st["sky"],
                  n_at_peak=st["n_at_peak"],
                  fwhm_px=(st["fwhm_px"] if np.isfinite(st["fwhm_px"])
                           else None))
            # Keep the stamp for the report's contact sheet (pixels, not
            # claims): one small npz per target, appended frame by frame.
            stamp_dir = PRODUCTS_DIR / "census"
            stamp_dir.mkdir(parents=True, exist_ok=True)
            pad = np.full((2 * H + 1, 2 * H + 1), np.nan, dtype=np.float32)
            pad[:stamp.shape[0], :stamp.shape[1]] = stamp
            np.save(stamp_dir / f"{target_key}_{rowid}.npy", pad)
    print(f"[S2:peakcensus] nothing left to do. ({n_done} frames this run)")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: peaklin  (comparison-star residual vs own peak, per mode)
# ---------------------------------------------------------------------------
def _store_curve(con, source: str, mode: str, curve: list[dict]) -> None:
    """Replace one (source, mode) deviation curve."""
    con.execute("DELETE FROM s2_linearity_curve WHERE source = ? AND mode = ?",
                (source, mode))
    con.executemany(
        "INSERT INTO s2_linearity_curve VALUES (?,?,?,?,?,?,?,?,?,?)",
        [(source, mode, c["lo"], c["hi"], c["peak_frac"], c["dev_pct"],
          None if not np.isfinite(c["dev_err_pct"]) else c["dev_err_pct"],
          c["n_points"], c["n_groups"], int(c["measured"])) for c in curve])


def _store_injection(con, source: str, mode: str, case: str,
                     rows: list[dict]) -> None:
    """Replace one (source, mode, case) injection-bias table."""
    con.execute("DELETE FROM s2_linearity_injection WHERE source = ? "
                "AND mode = ? AND case_name = ?", (source, mode, case))
    con.executemany(
        "INSERT INTO s2_linearity_injection VALUES (?,?,?,?,?,?,?,?,?,?)",
        [(source, mode, case, r["lo"], r["hi"], r["injected_pct"],
          r["recovered_pct"], r["bias_pct"],
          None if not np.isfinite(r["bias_err_pct"]) else r["bias_err_pct"],
          r["n_trials"]) for r in rows])


def rolloff_5pct(pf: np.ndarray) -> np.ndarray:
    """The injected test non-linearity: linear to 70% of scale, then a
    straight roll-off reaching -5% at the ceiling.  A shape a real sensor
    could plausibly have, and one that crosses the 1% criterion at 76% —
    inside the range where the downstream caps disagree (70/80/92%)."""
    return np.where(pf > 0.7, -0.05 * (pf - 0.7) / 0.3, 0.0)


def cmd_peaklin(con: sqlite3.Connection, cv_db: Path) -> int:
    """Residual-vs-own-peak for the CV comparison stars, per readout mode.

    Every solved CV series is a fixed-field time series: seeing moves each
    star's peak pixel up and down from frame to frame while its true flux
    stays put.  :func:`rlmt_diagnostics.linearity.ensemble_deviation`
    turns that into fractional deviation from linear response against the
    star's own peak, with the frame-to-frame terms (transparency, aperture
    loss) taken out by the stars whose peaks stay in the linear reference
    regime.

    Peaks are converted to a fraction of the RAW usable scale: locally
    reduced frames are raw minus a master dark (zero point 0); server-
    reduced frames are mapped back through the measured reduction of
    their era (``s2_recon_eras``: raw = F x (reduced - pedestal) + D).

    Writes one curve per series (source ``cv:<series>``), one pooled curve
    per mode (source ``cv``), and for each mode the null and roll-off
    injection tables computed on its largest series' real geometry.
    """
    if not Path(cv_db).exists():
        print(f"[S2:peaklin] CV database not found: {cv_db} — skipped.")
        return 0
    cv = sqlite3.connect(f"{Path(cv_db).resolve().as_uri()}?mode=ro", uri=True)
    modes = {m: (clip, veto) for m, clip, veto in con.execute(
        "SELECT mode, clip_adu, veto_adu FROM s2_ceiling_modes "
        "WHERE clip_adu IS NOT NULL")}
    recon = {e: (ped or 0.0, f, d) for e, ped, f, d in con.execute(
        "SELECT era_id, pedestal, flat_median, dark_median "
        "FROM s2_recon_eras WHERE flat_median IS NOT NULL")}
    con.execute("DELETE FROM s2_peaklin_series")
    con.execute("DELETE FROM s2_linearity_curve WHERE source LIKE 'cv%'")
    con.execute("DELETE FROM s2_linearity_injection WHERE source LIKE 'cv%'")
    pooled: dict[str, dict[str, list]] = {}
    largest: dict[str, tuple[int, dict]] = {}
    series = cv.execute(
        "SELECT series_key, era_id, filter, provenance FROM cv_series "
        "WHERE status = 'solved' ORDER BY series_key").fetchall()
    for s_idx, (skey, era, filt, prov) in enumerate(series):
        frames_ = cv.execute("""
            SELECT frame_id, readoutm, exptime, night, bkg_adu,
                   dark_median_adu
            FROM cv_frames WHERE series_key = ? AND status = 'matched'
            ORDER BY bjd_tdb""", (skey,)).fetchall()
        mode = ceil.mode_group(frames_[0][1]) if frames_ else None
        row = {"series_key": skey, "mode": mode, "provenance": prov,
               "era_id": era, "filter": filt, "n_frames": len(frames_)}

        def record(status: str, **kw):
            row.update(kw)
            cols = ("series_key", "mode", "provenance", "era_id", "filter",
                    "n_frames", "n_stars", "n_anchored", "n_points",
                    "n_points_above_ref", "max_peak_frac", "z_err_median",
                    "zero_adu", "scale_to_raw", "status")
            row["status"] = status
            con.execute(
                f"INSERT OR REPLACE INTO s2_peaklin_series "
                f"({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})",
                tuple(row.get(c) for c in cols))

        if len(frames_) < PEAKLIN_MIN_FRAMES:
            record("too_few_frames")
            continue
        if mode not in modes:
            record("no_ceiling_for_mode")
            continue
        clip, _veto = modes[mode]
        if prov == "server_reduced":
            if era not in recon:
                record("no_reduction_map")
                continue
            zero, scale, bias = recon[era]
        else:
            # raw minus master dark: zero point 0, unit scale; the bias is
            # the dark level that was subtracted.
            dk = [f[5] for f in frames_ if f[5] is not None]
            zero, scale = 0.0, 1.0
            bias = float(np.median(dk)) if dk else 0.0
        fid = {f[0]: j for j, f in enumerate(frames_)}
        dets = cv.execute(f"""
            SELECT d.frame_id, d.star_id, d.flux, d.fluxerr, d.peak
            FROM cv_detections d JOIN cv_frames f USING (frame_id)
            WHERE f.series_key = ? AND f.status = 'matched'
              AND d.star_id IS NOT NULL AND d.flux > 0""",
            (skey,)).fetchall()
        if not dets:
            record("no_detections")
            continue
        star_ids = sorted({d[1] for d in dets})
        sidx = {sid: i for i, sid in enumerate(star_ids)}
        shape = (len(star_ids), len(frames_))
        flux = np.full(shape, np.nan)
        ferr = np.full(shape, np.nan)
        peak = np.full(shape, np.nan)
        for frame_id, star_id, fl, fe, pk in dets:
            i, j = sidx[star_id], fid[frame_id]
            flux[i, j], ferr[i, j], peak[i, j] = fl, fe, pk
        # Coverage and brightness cuts (by median flux — a property of the
        # star, not of any one noisy measurement).
        cover = np.isfinite(flux).mean(axis=1)
        keep = np.flatnonzero(cover >= PEAKLIN_MIN_COVERAGE)
        if keep.size < PEAKLIN_MIN_STARS:
            record("too_few_stars", n_stars=int(keep.size))
            continue
        with lin._quiet_nan_warnings():
            medflux = np.nanmedian(flux[keep], axis=1)
        keep = keep[np.argsort(medflux)[::-1][:PEAKLIN_MAX_STARS]]
        flux, ferr, peak = flux[keep], ferr[keep], peak[keep]
        bkg = np.array([f[4] if f[4] is not None else np.nan
                        for f in frames_])
        expt = np.array([f[2] if f[2] else np.nan for f in frames_])
        # Peak pixel in frame units (sep's peak is above background), then
        # as raw signal above bias, then as a fraction of the usable scale.
        raw_signal = scale * (peak + bkg[None, :] - zero)
        pfrac = raw_signal / (clip - bias)
        res = lin.ensemble_deviation(flux, ferr, expt, pfrac)
        nights = sorted({f[3] for f in frames_})
        nidx = np.array([nights.index(f[3]) for f in frames_])
        # Independent unit = star x night, made unique across series.
        group = (s_idx * 10_000_000
                 + np.arange(flux.shape[0])[:, None] * 1000 + nidx[None, :])
        ok = np.isfinite(res["dev"]) & np.isfinite(pfrac)
        curve = lin.deviation_curve(pfrac[ok], res["dev"][ok], group[ok])
        _store_curve(con, f"cv:{skey}", mode, curve)
        with lin._quiet_nan_warnings():
            record("ok", n_stars=int(flux.shape[0]),
                   n_anchored=int(np.isfinite(res["base"]).sum()),
                   n_points=int(ok.sum()),
                   n_points_above_ref=int(
                       (ok & (pfrac >= lin.REF_PEAK_FRACTION)).sum()),
                   max_peak_frac=float(np.nanmax(pfrac[ok])) if ok.any()
                   else None,
                   z_err_median=float(np.nanmedian(res["z_err"])),
                   zero_adu=zero, scale_to_raw=scale)
        acc = pooled.setdefault(mode, {"pf": [], "dev": [], "grp": []})
        acc["pf"].append(pfrac[ok])
        acc["dev"].append(res["dev"][ok])
        acc["grp"].append(group[ok])
        n_above = int((ok & (pfrac >= lin.REF_PEAK_FRACTION)).sum())
        if mode not in largest or n_above > largest[mode][0]:
            with np.errstate(invalid="ignore", divide="ignore"):
                model = res["base"][:, None] * res["z"][None, :] * expt[None, :]
                ppf = pfrac / model
            largest[mode] = (n_above, {
                "series": skey, "flux_err": ferr, "exptime": expt,
                "base": res["base"], "z": res["z"], "ppf": ppf})
        con.commit()
        print(f"[S2:peaklin]   {skey} ({mode}): {int(ok.sum()):,} points, "
              f"{n_above:,} above the reference regime")
    for mode, acc in sorted(pooled.items()):
        curve = lin.deviation_curve(np.concatenate(acc["pf"]),
                                    np.concatenate(acc["dev"]),
                                    np.concatenate(acc["grp"]))
        _store_curve(con, "cv", mode, curve)
    # Injection on each mode's richest series: real stars, real frame
    # factors, real error floors, real seeing — only the fluxes replaced.
    for mode, (_n, geo) in sorted(largest.items()):
        anchored = np.isfinite(geo["base"])
        # Each star's own typical peak-to-flux ratio stands in where a
        # measurement is missing.
        with lin._quiet_nan_warnings():
            ppf_star = np.nanmedian(geo["ppf"], axis=1)
        ppf = np.where(np.isfinite(geo["ppf"]), geo["ppf"],
                       ppf_star[:, None])
        good = anchored & np.isfinite(ppf_star)
        if int(good.sum()) < PEAKLIN_MIN_STARS:
            continue
        zz = np.where(np.isfinite(geo["z"]), geo["z"], 1.0)
        ee = np.where(np.isfinite(geo["flux_err"][good]),
                      geo["flux_err"][good],
                      np.nanmedian(geo["flux_err"][good]))
        tt = np.where(np.isfinite(geo["exptime"]), geo["exptime"], 1.0)
        for case, fn in (("null", lambda pf: 0.0 * pf),
                         ("rolloff_5pct", rolloff_5pct)):
            rows = lin.injection_bias(
                ee, tt, geo["base"][good], zz, ppf[good], fn,
                n_trials=PEAKLIN_INJECTION_TRIALS)
            _store_injection(con, f"cv:{geo['series']}", mode, case, rows)
    con.commit()
    cv.close()
    print(f"[S2:peaklin] done ({len(series)} series, "
          f"{len(pooled)} modes pooled).")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: badpix  (bad-pixel mask per camera / geometry / temperature)
# ---------------------------------------------------------------------------
def _badpix_groups(con: sqlite3.Connection) -> dict[tuple, dict]:
    """Per (camera, naxis1, naxis2, temp group): science and dark samples.

    Science sample: one frame per (night, target) scene, exposure at least
    :data:`BADPIX_MIN_EXPTIME_S`, non-grism, spread evenly over the
    group's date range — different pointings, so a star never repeats at
    one pixel.  Dark sample: the group's most populous dark exposure of at
    least that length, for the independent check and the RTS mask.
    """
    from collections import defaultdict
    rows = con.execute(f"""
        SELECT f.path, f.instrume, f.readoutm, f.xbinning, f.naxis1,
               f.naxis2, f.filter, f.exptime, f.jd, f.night, f.target_key,
               f.ccd_temp
        FROM frames f
        WHERE f.is_canonical = 1 AND f.tree = 'rawimage'
          AND {SCIENCE_IMAGETYP} AND f.target_key IS NOT NULL
          AND f.jd IS NOT NULL AND f.naxis1 >= 1000 AND f.naxis2 >= 1000
          AND f.exptime >= ?
        ORDER BY f.jd, f.path""", (BADPIX_MIN_EXPTIME_S,)).fetchall()
    groups: dict[tuple, dict] = defaultdict(
        lambda: {"science": {}, "n": 0, "xbin": 1})
    for (path, instr, mode, xb, n1, n2, filt, expt, jd, night, tkey,
         temp) in rows:
        if (filt or "").strip().lower() in SKY_EXCLUDED_FILTERS:
            continue
        camera = cam.camera_of(instr, n1, mode)
        if camera == "unknown":
            continue
        key = (camera, int(n1), int(n2), temp_group_of(temp))
        g = groups[key]
        g["n"] += 1
        g["xbin"] = int(xb or 1)
        g["science"].setdefault((night, tkey), (path, exptime_bin(expt)))
    out = {}
    for key, g in groups.items():
        if g["n"] < BADPIX_MIN_GROUP_FRAMES:
            continue
        camera, n1, n2, tg = key
        darks = con.execute("""
            SELECT f.path, f.exptime, f.instrume, f.readoutm, f.ccd_temp
            FROM frames f
            WHERE f.is_canonical = 1 AND f.imagetyp LIKE 'Dark%'
              AND f.naxis1 = ? AND f.naxis2 = ? AND f.exptime >= ?
              AND lower(f.basename) NOT LIKE 'master%'
            ORDER BY f.jd, f.path""",
            (n1, n2, BADPIX_MIN_EXPTIME_S)).fetchall()
        by_exp: dict[float, list[str]] = defaultdict(list)
        for path, expt, instr, mode, temp in darks:
            if (cam.camera_of(instr, n1, mode) == camera
                    and temp_group_of(temp) == tg
                    and "StackPro" not in (mode or "")):
                by_exp[exptime_bin(expt)].append(path)
        dark_exp, dark_paths = None, []
        if by_exp:
            dark_exp = max(by_exp, key=lambda e: (len(by_exp[e]) >= 8,
                                                  e, len(by_exp[e])))
            dark_paths = badpix.sample_evenly(by_exp[dark_exp],
                                              BADPIX_DARK_FRAMES)
        out[key] = {
            "science": badpix.sample_evenly(
                list(g["science"].values()), BADPIX_SCIENCE_FRAMES),
            "n_candidates": g["n"], "xbin": g["xbin"],
            "dark_exptime": dark_exp, "darks": dark_paths}
    return out


def orientation_tag(header) -> str:
    """Image-orientation epoch of a frame (``FM``/``N``) from FLIPSTAT.

    Part of a mask's identity: the first ASI mask built without it came
    out EMPTY — every hot pixel fired in exactly half of the sample —
    because the acquisition software's flip setting changed between May
    and November 2025, so one sensor pixel is stored at two file positions.
    """
    return badpix.orientation_of(header.get("FLIPSTAT", ""))


def cmd_badpix(con: sqlite3.Connection, archive: Path, batch: int,
               only: str | None = None) -> int:
    """Build the bad-pixel mask products.

    A group is camera x geometry x temperature; within it the frames are
    split by image orientation (:func:`orientation_tag`) and each
    orientation with at least 8 readable science frames gets its own mask
    ``products/detector/badpix/<key>.npz`` and ``s2_badpix_masks`` row.
    A group is all-or-nothing and is skipped once any of its masks exists.
    """
    done_groups = {r[0] for r in con.execute(
        "SELECT DISTINCT camera || '_' || naxis1 || 'x' || naxis2 || '_T' "
        "|| temp_group FROM s2_badpix_masks")}
    out_dir = REPO_ROOT / badpix.BADPIX_SUBDIR
    n_groups = 0
    for key, g in sorted(_badpix_groups(con).items(), key=lambda kv: str(kv[0])):
        camera, n1, n2, tg = key
        gkey = badpix.mask_key(camera, n1, n2, tg)
        if gkey in done_groups or (only and only not in gkey):
            continue
        if n_groups >= batch:
            print("[S2:badpix] batch cap reached; re-run to continue.")
            return 0
        n_native = 4 if (camera in ("ASI", "QHY600") and g["xbin"] == 2) else 1
        print(f"[S2:badpix] {gkey}: {len(g['science'])} science frames, "
              f"{len(g['darks'])} darks ({g['dark_exptime']}s), "
              f"n_native {n_native}")
        accs: dict[str, badpix.PersistenceAccumulator] = {}
        reps: dict[str, tuple[float, np.ndarray]] = {}
        paths: dict[str, list[str]] = {}
        frame_rows: dict[str, list] = {}
        for path, expt in g["science"]:
            try:
                img, hdr = read_image(archive, path)
                tag = orientation_tag(hdr)
                acc = accs.setdefault(
                    tag, badpix.PersistenceAccumulator((n2, n1), n_native))
                st = acc.add(img)
            except Exception as e:
                print(f"[S2:badpix]   SKIP {path}: {e}")
                continue
            frame_rows.setdefault(tag, []).append(
                ("science", path, expt, st["sigma"], st["threshold"],
                 st["n_fired"], st["n_rail"]))
            paths.setdefault(tag, []).append(path)
            # Representative frame for the rail evidence: the QUIETEST one
            # (lowest excess sigma = darkest, emptiest sky).
            if tag not in reps or st["sigma"] < reps[tag][0]:
                reps[tag] = (st["sigma"], np.asarray(img, dtype=np.float32))
        # Dark stack, split by orientation the same way.
        daccs: dict[str, badpix.PersistenceAccumulator] = {}
        stacks: dict[str, list] = {}
        for path in g["darks"]:
            try:
                img, hdr = read_image(archive, path)
                tag = orientation_tag(hdr)
                st = daccs.setdefault(
                    tag, badpix.PersistenceAccumulator((n2, n1), n_native)
                ).add(img)
            except Exception as e:
                print(f"[S2:badpix]   SKIP dark {path}: {e}")
                continue
            stacks.setdefault(tag, []).append(np.asarray(img, dtype=np.uint16))
            frame_rows.setdefault(tag, []).append(
                ("dark", path, g["dark_exptime"], st["sigma"],
                 st["threshold"], st["n_fired"], st["n_rail"]))
        n_groups += 1
        for tag, acc in sorted(accs.items()):
            mkey = f"{gkey}_{tag}"
            if acc.n_frames < 8:
                print(f"[S2:badpix]   {mkey}: only {acc.n_frames} readable "
                      "science frames in this orientation — no mask.")
                continue
            mask = acc.mask()
            rail_k = acc.rail_k() if n_native > 1 else None
            n_dark, n_dark_hot, n_both, n_noisy = 0, None, None, 0
            dacc = daccs.get(tag)
            if dacc is not None and dacc.n_frames >= 4:
                n_dark = dacc.n_frames
                agree = badpix.mask_agreement(
                    (mask & badpix.BAD_HOT) > 0,
                    (dacc.mask() & badpix.BAD_HOT) > 0)
                n_dark_hot, n_both = agree["n_b"], agree["n_both"]
                if n_dark >= badpix.NOISY_MIN_FRAMES:
                    cube = np.stack(stacks[tag])
                    sigs = []
                    for y0 in range(0, n2, 256):   # bounded working copy
                        blk = cube[:, y0:y0 + 256, :].astype(np.float32)
                        med = np.median(blk, axis=0)
                        sigs.append(1.4826 * np.median(
                            np.abs(blk - med[None]), axis=0))
                    sig = np.concatenate(sigs, axis=0)
                    noisy = sig > badpix.NOISY_FACTOR * max(
                        float(np.median(sig)), 1e-6)
                    mask[noisy] |= badpix.BAD_NOISY
                    n_noisy = int(noisy.sum())
                    del cube
            # Rail evidence for D2: where rail-1 pixels actually sit vs
            # where (65535 + 3 x neighbour level)/4 says they should.
            rail1_med = rail1_exp = ped = None
            if n_native > 1:
                rep = reps[tag][1]
                _exc, med_img = badpix.neighbour_excess(rep)
                sel = ((mask & badpix.BAD_RAIL) > 0) & (rail_k == 1)
                if sel.any():
                    rail1_med = float(np.median(rep[sel]))
                    rail1_exp = float(np.median(
                        (badpix.NATIVE_FULL_SCALE + 3.0 * med_img[sel]) / 4.0))
                    ped = float(np.median(med_img[sel]))
            n_hot = int(((mask & badpix.BAD_HOT) > 0).sum())
            n_rail = int(((mask & badpix.BAD_RAIL) > 0).sum())
            rk = [int((((mask & badpix.BAD_RAIL) > 0) & (rail_k == k)).sum())
                  if rail_k is not None else 0 for k in (1, 2, 3)]
            n_bad = int((mask > 0).sum())
            meta = {"mask_key": mkey, "camera": camera, "naxis1": n1,
                    "naxis2": n2, "temp_group": tg, "orientation": tag,
                    "n_native": n_native, "n_science": acc.n_frames,
                    "n_dark": n_dark, "dark_exptime": g["dark_exptime"],
                    "hot_sigma": badpix.HOT_SIGMA,
                    "hot_min_excess_adu": badpix.HOT_MIN_EXCESS_ADU,
                    "persist_fraction": badpix.PERSIST_FRACTION,
                    "rail_tol_adu": badpix.RAIL_TOL_ADU,
                    "noisy_factor": badpix.NOISY_FACTOR,
                    "code_version": S2_CODE_VERSION,
                    "science_paths": paths[tag]}
            npz = badpix.save_mask(out_dir, mkey, mask, acc.hot, rail_k, meta)
            con.executemany(
                "INSERT OR REPLACE INTO s2_badpix_frames "
                "(mask_key, role, path, exptime, sigma, threshold, n_fired, "
                "n_rail) VALUES (?,?,?,?,?,?,?,?)",
                [(mkey,) + r for r in frame_rows.get(tag, [])])
            con.execute("""INSERT OR REPLACE INTO s2_badpix_masks
                (mask_key, camera, naxis1, naxis2, temp_group, n_native,
                 n_science, n_hot, n_rail, n_rail1, n_rail2, n_rail3,
                 n_noisy, n_bad, bad_fraction, n_dark, dark_exptime,
                 n_dark_hot, n_dark_and_science, rail1_median_adu,
                 rail1_expected_adu, pedestal_adu, npz_path)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (mkey, camera, n1, n2, tg, n_native, acc.n_frames, n_hot,
                 n_rail, rk[0], rk[1], rk[2], n_noisy, n_bad,
                 n_bad / float(n1 * n2), n_dark, g["dark_exptime"],
                 n_dark_hot, n_both, rail1_med, rail1_exp, ped,
                 str(npz.relative_to(REPO_ROOT))))
            con.commit()
            print(f"[S2:badpix]   {mkey}: {acc.n_frames} frames, {n_bad:,} "
                  f"bad ({100.0 * n_bad / (n1 * n2):.4f}%): hot {n_hot:,}, "
                  f"rail {n_rail:,} (k=1/2/3: {rk[0]}/{rk[1]}/{rk[2]}), "
                  f"noisy {n_noisy:,}")
    print(f"[S2:badpix] nothing left to do. ({n_groups} groups this run)")
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
# params helpers: flat-pair PTC -> s2_flat_fits, s2_camera_configs
# ---------------------------------------------------------------------------
#: A flats-only (free-intercept) fit is attempted only when the flat
#: points span at least this ratio of signal; below it the slope is a
#: property of the noise in the fit (the same rule the noise curve uses).
FLAT_ONLY_MIN_SPAN = 1.5

#: Sky-pair gain: a pair qualifies when its sky signal is at least this
#: multiple of the zero-signal variance (shot variance comparable to or
#: above the read variance — below that K_i = S/(V - V0) is a ratio of
#: two noises), and its two frames' levels agree to this fraction.
#: DISCLOSURE: the multiple was first set to 2.0 and lowered to 1.5 after
#: the per-pair gains of the two configurations whose truth is known from
#: flats (AC4040 High Gain e1.057, ASI Mode0) had been looked at: their
#: K_i are flat down to S ~ 1.5 V0 and fall away only below ~0.7 V0, and
#: at 2.0 the AC4040 validation set had three pairs, one short of
#: SKY_MIN_PAIRS.  It was tuned on validation configurations only, never
#: on a configuration whose gain the sky method decides.
SKY_MIN_SIGNAL_OVER_V0 = 1.5
SKY_MAX_LEVEL_CHANGE = 0.10

#: Robust clip of the per-pair gains (MAD-sigmas) and the fewest pairs
#: that make a sky gain.
SKY_CLIP_SIGMA = 5.0
SKY_MIN_PAIRS = 4

#: Longest exposure a sky pair may have to serve as the PSEUDO-BIAS of a
#: camera that has no bias or dark frames at all (the QHY600).
PSEUDO_BIAS_MAX_EXPTIME_S = 1.0

#: Fewest flat pairs for a configuration's flat-pair gain to be adopted.
FLAT_MIN_PAIRS = 3

#: Flats taken THROUGH A GRISM are excluded from every gain fit.  A
#: dispersed flat is not flat: it carries the sharp edges of the spectrum
#: and of the grism's vignetting lozenge, and a sub-pixel flexure shift
#: between the two exposures turns those edges into variance that is not
#: photon noise.  Measured on the iKon (the only camera with grism flats):
#: per-pair K = 0.88-0.96 through OGGrism/HaGrism against 0.975-0.985
#: through y, g, i and r on the same nights.  DISCLOSURE: this exclusion was
#: written AFTER those per-pair values were seen; it is the same physical
#: rule the sky-pair stage applied from the start (SKY_EXCLUDED_FILTERS).
#: The pairs stay in s2_flat_pairs; only the fits ignore them.

#: A flat pair whose DIFFERENCE image has a nearest-neighbour correlation
#: above this is not a pair of flats: something structured changed between
#: the two frames (the 2024-05-15 iKon grism "flats" read 0.04-0.19, with
#: per-pair gains scattered over 10%, against 0.002-0.02 for every clean
#: pair of the same camera).  Such pairs stay in s2_flat_pairs with their
#: numbers but are excluded from every fit.
FLAT_MAX_RHO = 0.03

#: The flat groups' additive variance excess over read noise counts as
#: detected above this many sigma (see the adoption rule in _params_flat).
EXCESS_SIGNIFICANCE = 2.0

#: A within-group fit weaker than this relative error is not adopted on
#: its own (it then has too little within-group signal spread).
WITHIN_MAX_REL_ERR = 0.05

#: Above this significance the additive excess is beyond argument and the
#: pinned fit is simply biased; between EXCESS_SIGNIFICANCE and this, the
#: within-group fit is still adopted but the difference between the two
#: estimators is carried as a systematic (a 2-sigma excess is a reason to
#: prefer the safer estimator, not a licence to ignore the other).
EXCESS_CLEAR = 5.0


def _points(con, config: str, kinds: tuple[str, ...],
            extra: str = "", params: tuple = ()) -> list[dict]:
    """PTC points of one configuration for the given pair kinds."""
    marks = ", ".join("?" * len(kinds))
    grism = ", ".join(f"'{g}'" for g in SKY_EXCLUDED_FILTERS)
    return [{"signal": r[0], "var": r[1], "var_err": r[2], "pair_id": r[3],
             "group": r[4]}
            for r in con.execute(f"""
                SELECT q.signal_adu, q.var_adu2, q.var_err_adu2, q.pair_id,
                       p.group_key
                FROM s2_flat_points q JOIN s2_flat_pairs p USING (pair_id)
                WHERE p.config = ? AND p.status = 'ok'
                  AND (p.rho_x IS NULL OR max(p.rho_x, p.rho_y) <= ?)
                  AND NOT (p.kind = 'flat' AND lower(coalesce(p.filter, ''))
                           IN ({grism}))
                  AND p.kind IN ({marks}) {extra}""",
                (config, FLAT_MAX_RHO) + kinds + params)]


def _store_fit(con, config: str, estimator: str, fit: dict | None) -> None:
    """Record one estimator's fit (chi2 and dof included, always)."""
    if fit is None:
        return
    fit = {**{"gain_err_boot": None, "gain_boot_bias": None}, **fit}
    con.execute(
        "INSERT OR REPLACE INTO s2_flat_fits VALUES "
        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (config, estimator, fit["gain"], flatptc.adopted_gain_error(fit),
         fit["gain_err_formal"], fit["gain_err_boot"], fit["gain_boot_bias"],
         fit["intercept"], fit["intercept_err"], fit["read_noise_adu"],
         fit["read_noise_adu_err"], fit["chi2"], fit["dof"], fit["chi2nu"],
         fit["quad_coeff"], fit["quad_coeff_err"], fit["quad_z"],
         fit["gain_quad"], fit["n_points"], fit["n_pairs"],
         fit["signal_lo"], fit["signal_hi"]))


def _config_census(con) -> dict[str, dict]:
    """Every camera/configuration in the archive, with its frame census.

    Built from the manifest's ``frames`` table by grouping on the header
    fields :func:`rlmt_diagnostics.camera.config_key` reads — so a
    configuration with NO calibration frames still gets a row (and a
    verdict saying so) instead of being absent from the table.
    """
    out: dict[str, dict] = {}
    for (instr, mode, xb, n1, eg, n, first, last, eras) in con.execute("""
            SELECT f.instrume, f.readoutm, f.xbinning, f.naxis1,
                   round(f.egain, 4), count(*), min(f.night), max(f.night),
                   group_concat(DISTINCT f.era_id)
            FROM frames f
            WHERE f.is_canonical = 1 AND f.tree = 'rawimage'
              AND f.naxis1 >= 1000 AND f.naxis2 >= 1000
            GROUP BY f.instrume, f.readoutm, f.xbinning, f.naxis1,
                     round(f.egain, 4)"""):
        config = cam.config_key(instr, mode, xb, n1, eg)
        if config.startswith("unknown"):
            continue
        c = out.setdefault(config, {"n_frames": 0, "first": first,
                                    "last": last, "eras": set(),
                                    "egain": eg})
        c["n_frames"] += n
        c["first"] = min(filter(None, (c["first"], first)), default=None)
        c["last"] = max(filter(None, (c["last"], last)), default=None)
        c["eras"].update(e for e in (eras or "").split(",") if e)
    return out


def _zero_references(con) -> dict[str, dict]:
    """Per configuration: single-read bias level and zero-signal variance.

    In order of preference:

    1. the configuration's own zero-signal pairs (bias / shortest dark),
       divided by the number of StackPro sub-reads when it is a sum;
    2. a header-refined sibling (``iKon 1MHz 4x`` serves ``iKon 1MHz``);
    3. for StackPro, the plain High Gain configuration of the same EGAIN
       epoch (same sensor, same gain, one read);
    4. for a camera with no zero-signal frames at all (the QHY600): the
       shortest-exposure sky pair of that camera, if at most
       :data:`PSEUDO_BIAS_MAX_EXPTIME_S` long — a PSEUDO-bias whose level
       and variance are upper bounds on the true ones.  Flagged as such.
    """
    out: dict[str, dict] = {}
    for config, expt, lvl, var, n in con.execute("""
            SELECT p.config, min(p.exptime),
                   avg((p.level_a + p.level_b) / 2.0), avg(q.var_adu2),
                   count(*)
            FROM s2_flat_pairs p JOIN s2_flat_points q USING (pair_id)
            WHERE p.kind IN ('bias', 'dark') AND p.status = 'ok'
            GROUP BY p.config"""):
        nsub = (ptc.stackpro_nsub_for_exptime(expt)
                if "StackPro" in config else 1)
        out[config] = {"bias": lvl / nsub, "v0": var / nsub,
                       "basis": f"{n} zero-signal pairs"
                       + (f" (/{nsub} sub-reads)" if nsub > 1 else "")}
    for config in list(out):
        base = " ".join(config.split()[:2])
        if config.startswith("iKon") and base not in out:
            out[base] = {**out[config], "basis": out[config]["basis"]
                         + f" of {config}"}
    for epoch in ("e1.054", "e1.057"):
        sp, hg = f"AC4040 StackPro {epoch}", f"AC4040 High Gain {epoch}"
        if sp not in out and hg in out:
            out[sp] = {**out[hg], "basis": out[hg]["basis"] + f" of {hg}"}
        if hg not in out and sp in out:
            out[hg] = {**out[sp], "basis": out[sp]["basis"] + f" of {sp}"}
    # Pseudo-bias for cameras with no zero-signal frames.
    cams_done = {cam.camera_name(c) for c in out}
    for config, expt, lvl, var in con.execute("""
            SELECT p.config, p.exptime, (p.level_a + p.level_b) / 2.0,
                   avg(q.var_adu2)
            FROM s2_flat_pairs p JOIN s2_flat_points q USING (pair_id)
            WHERE p.kind = 'sky' AND p.status = 'ok' AND p.exptime <= ?
            GROUP BY p.pair_id ORDER BY p.exptime, 3""",
            (PSEUDO_BIAS_MAX_EXPTIME_S,)).fetchall():
        camera = cam.camera_name(config)
        if camera in cams_done:
            continue
        cams_done.add(camera)
        ref = {"bias": lvl, "v0": var, "pseudo": True,
               "basis": f"PSEUDO-bias: the {expt:g} s sky pair of {config} "
                        "(no bias or dark frames exist for this camera; "
                        "level and variance are upper bounds)"}
        for other, in con.execute(
                "SELECT DISTINCT config FROM s2_flat_pairs WHERE kind = 'sky'"
                ).fetchall():
            if cam.camera_name(other) == camera and other not in out:
                out[other] = ref
    return out


def _sky_gain(con, config: str, zref: dict | None) -> dict | None:
    """Robust sky-pair gain of one configuration; fills ``s2_sky_gains``.

    WHY NOT A LINE FIT.  The first version fitted one free-intercept line
    through every sky pair of a configuration.  Its chi2/dof ran to 10^4
    and its validation against flats failed (ASI: 1.81 against 1.05),
    because sky pairs do not share an intercept: StackPro frames carry
    1, 4, 8 or 16 pedestals depending on exposure, the QHY changed OFFSET
    after its first two nights, and a few frames in ``rawimage`` are
    already bias-subtracted.  Looked at pair by pair, though, every
    signal-dominated pair gave a sensible gain.  So that is the estimator:

        K_i = S_i / (V_i - N_i V0),   S_i = level_i - N_i bias

    per pair (N_i = StackPro sub-reads, else 1), for pairs with
    ``S_i >= SKY_MIN_SIGNAL_OVER_V0 x N_i V0`` and a level change under
    :data:`SKY_MAX_LEVEL_CHANGE`; then the MEDIAN of the K_i after a
    :data:`SKY_CLIP_SIGMA` MAD clip, with the standard error of a median
    from their scatter.  None with fewer than :data:`SKY_MIN_PAIRS`.
    """
    con.execute("DELETE FROM s2_sky_gains WHERE config = ?", (config,))
    if zref is None:
        return None
    rows = con.execute("""
        SELECT p.pair_id, p.exptime, p.ratio, coalesce(p.bias_adu, 0.0),
               sum(q.signal_adu * q.n_tiles) / sum(q.n_tiles),
               sum(q.var_adu2 * q.n_tiles) / sum(q.n_tiles)
        FROM s2_flat_pairs p JOIN s2_flat_points q USING (pair_id)
        WHERE p.config = ? AND p.kind = 'sky' AND p.status = 'ok'
        GROUP BY p.pair_id""", (config,)).fetchall()
    recs = []
    for pid, expt, ratio, bias_used, sig, var in rows:
        nsub = (ptc.stackpro_nsub_for_exptime(expt)
                if "StackPro" in config else 1)
        bias, v0 = nsub * zref["bias"], nsub * zref["v0"]
        s_true = sig + bias_used - bias          # undo the bias used then
        reason = None
        if ratio is None or abs(ratio - 1.0) > SKY_MAX_LEVEL_CHANGE:
            reason = "level_changed"
        elif s_true < SKY_MIN_SIGNAL_OVER_V0 * v0:
            reason = "read_noise_dominated"
        elif var <= v0:
            reason = "variance_below_floor"
        k_i = (s_true / (var - v0)) if var > v0 else None
        recs.append([pid, config, nsub, s_true, var, bias, v0, k_i,
                     0 if reason else 1, reason])
    ks = np.array([r[7] for r in recs if r[8]], dtype=np.float64)
    result = None
    if ks.size >= SKY_MIN_PAIRS:
        keep = np.ones(ks.size, dtype=bool)
        for _ in range(5):
            med = np.median(ks[keep])
            mad = 1.4826 * np.median(np.abs(ks[keep] - med))
            new = np.abs(ks - med) <= SKY_CLIP_SIGMA * max(mad, 1e-6)
            if np.array_equal(new, keep):
                break
            keep = new
        j = 0
        for r in recs:
            if r[8]:
                if not keep[j]:
                    r[8], r[9] = 0, "outlier"
                j += 1
        good = ks[keep]
        if good.size >= SKY_MIN_PAIRS:
            med = float(np.median(good))
            mad = 1.4826 * float(np.median(np.abs(good - med)))
            spread = mad if mad > 0 else float(good.std(ddof=1))
            result = {"gain": med,
                      "gain_err": 1.2533 * spread / np.sqrt(good.size),
                      "n_pairs": int(good.size),
                      "n_rejected": int(len(recs) - good.size),
                      "basis": zref["basis"],
                      "pseudo": bool(zref.get("pseudo"))}
    con.executemany("INSERT OR REPLACE INTO s2_sky_gains VALUES "
                    "(?,?,?,?,?,?,?,?,?,?)", recs)
    return result


def _params_flat(con, put) -> None:
    """Distil flat/sky pairs into fits and the per-configuration table.

    THE ADOPTION RULE.  Rules 2-3 were written before any fit was looked
    at.  Rule 1 was REVISED once, after the first AC4040 fit: the single
    line through flat and bias points gave chi2/dof = 500/82 with a
    strongly negative quadratic term, the signature of an additive
    variance in the twilight flats (faint stars, present twice in a
    difference of two sky frames) that a bias pair does not have.  The
    revision adds an estimator immune to that (free additive term per
    flat group) and a test for when it must be used; it was fixed before
    any other camera's flats had been fitted.

    1. A configuration with at least :data:`FLAT_MIN_PAIRS` usable flat
       pairs gets two fits:

       * ``flat-within`` — free additive term per flat group
         (:func:`rlmt_diagnostics.flatptc.fit_ptc_within_groups`);
       * ``flat+zero``   — one line through flat points and zero-signal
         (bias / shortest-dark) points.

       The additive excess of the flat groups over the read-noise variance
       is then measured (:func:`rlmt_diagnostics.flatptc.additive_excess`).
       If it is significant (> :data:`EXCESS_SIGNIFICANCE` sigma), or no
       pinned fit exists, ``flat-within`` is adopted.  If it is not
       significant, the more precise of the two is adopted and their
       difference enters the systematic error.  If ``flat-within`` does
       not exist or is weaker than :data:`WITHIN_MAX_REL_ERR` (flats at a
       single level have no within-group leverage) ``flat+zero`` is
       adopted and the note says the additive term could not be tested.

       Statistical error: the larger of the formal and the pair-bootstrap
       error.  Systematic error, in quadrature: the estimator difference
       above (when applicable); the EXCESS filter-to-filter scatter of the
       within-group gain (a gain cannot depend on wavelength — what leaks
       past the clip can); and 4 x |rho| x K for the measured
       nearest-neighbour pixel correlation rho of the difference images
       (4 x sigma_rho x K when rho is consistent with zero).
    2. A configuration with no flats adopts the SKY-pair free-intercept
       fit, multiplied by the mean flat/sky gain ratio measured on the
       configurations that have both, with the scatter of those ratios
       (or, with a single validation configuration, its own |ratio - 1|
       and error) added as the systematic.  With no validation
       configuration at all the sky value is recorded but NOT adopted.
    3. Otherwise the configuration has no gain measurement and says so.
    """
    con.execute("DELETE FROM s2_flat_fits")
    con.execute("DELETE FROM s2_camera_configs")
    census = _config_census(con)
    measured = {r[0] for r in con.execute(
        "SELECT DISTINCT config FROM s2_flat_pairs WHERE status = 'ok'")}
    ceilings = dict(con.execute(
        "SELECT mode, clip_adu FROM s2_ceiling_modes "
        "WHERE clip_adu IS NOT NULL"))
    egain_clip = {(m, round(e, 3)): v for m, e, v in con.execute(
        "SELECT mode, egain, median_max_adu FROM s2_ceiling_egain")}
    fits: dict[tuple, dict] = {}
    for config in sorted(measured | set(census)):
        flat = _points(con, config, ("flat",))
        zero = _points(con, config, ("bias", "dark"))
        n_flat_pairs = len({p_["pair_id"] for p_ in flat})
        if n_flat_pairs >= FLAT_MIN_PAIRS:
            fits[(config, "flat+zero")] = flatptc.fit_ptc_line(flat + zero) \
                if zero else None
            sig = [p_["signal"] for p_ in flat]
            if min(sig) > 0 and max(sig) / min(sig) >= FLAT_ONLY_MIN_SPAN:
                fits[(config, "flat")] = flatptc.fit_ptc_line(flat)
            fits[(config, "flat-within")] = \
                flatptc.fit_ptc_within_groups(flat)
            # Sub-fits: per temperature group and per filter, as
            # consistency checks (never adopted on their own).
            for tg, in con.execute(
                    "SELECT DISTINCT temp_group FROM s2_flat_pairs WHERE "
                    "config = ? AND kind = 'flat' AND status = 'ok' "
                    "AND temp_group IS NOT NULL", (config,)).fetchall():
                sub = _points(con, config, ("flat", "bias", "dark"),
                              "AND p.temp_group = ?", (tg,))
                fits[(config, f"flat+zero@T{tg}")] = flatptc.fit_ptc_line(sub)
            for filt, in con.execute(
                    "SELECT DISTINCT filter FROM s2_flat_pairs WHERE "
                    "config = ? AND kind = 'flat' AND status = 'ok'",
                    (config,)).fetchall():
                only = _points(con, config, ("flat",), "AND p.filter = ?",
                               (filt,))
                fits[(config, f"flat+zero:{filt}")] = \
                    flatptc.fit_ptc_line(only + zero)
                fits[(config, f"flat-within:{filt}")] = \
                    flatptc.fit_ptc_within_groups(only, n_boot=200)
    for (config, est), fit in fits.items():
        _store_fit(con, config, est, fit)

    # ---- sky-method validation: configurations with both -----------------
    ratios = []
    zrefs = _zero_references(con)
    adopted: dict[str, dict] = {}
    sky: dict[str, dict] = {}
    for config in sorted(measured):
        adopted[config] = _adopt_flat_gain(con, config, fits)
        sg = _sky_gain(con, config, zrefs.get(config))
        if sg:
            sky[config] = sg
            con.execute(
                "INSERT OR REPLACE INTO s2_flat_fits (config, estimator, "
                "gain, gain_err, n_pairs) VALUES (?,?,?,?,?)",
                (config, "sky-pairs", sg["gain"], sg["gain_err"],
                 sg["n_pairs"]))
    # ---- sky-method validation: configurations with flats AND sky --------
    for config in sorted(measured):
        ad, sg = adopted[config], sky.get(config)
        if ad["gain"] is not None and sg:
            r = ad["gain"] / sg["gain"]
            r_err = r * float(np.hypot(ad["stat"] / ad["gain"],
                                       sg["gain_err"] / sg["gain"]))
            ratios.append((config, r, r_err))
            put(config, "sky_method_ratio", r, r_err,
                "flat-pair gain / sky-pair gain on this configuration "
                "(validation of the sky-pair method; 1.0 = the two agree)",
                "s2_flat_fits, s2_sky_gains")
    if ratios:
        rr = np.array([r for _c, r, _e in ratios])
        corr = float(rr.mean())
        corr_sys = (float(rr.std(ddof=1)) if rr.size > 1
                    else float(np.hypot(abs(rr[0] - 1.0), ratios[0][2])))
    else:
        corr, corr_sys = None, None

    # A census key the header-refined measurements supersede ('iKon 1MHz'
    # when every measured pair is 'iKon 1MHz 4x' AND the sky pairs — read
    # from science frames — carry that same refined key): its frame census
    # moves to the refined key.
    for base in list(census):
        if not base.startswith("iKon"):
            continue
        refined = [c for c in measured if c.startswith(base + " ")
                   and c in sky]
        if base not in measured and len(refined) == 1:
            census[refined[0]] = census.pop(base)

    # ---- per-configuration adoption --------------------------------------
    for config in sorted(measured | set(census)):
        meta = census.get(config, {})
        mode = cam.mode_label_of(config)
        n_native = cam.n_native_per_pixel(config)
        sg = sky.get(config)
        zero_rows = con.execute("""
            SELECT q.var_adu2, q.var_err_adu2,
                   (p.level_a + p.level_b) / 2.0
            FROM s2_flat_points q JOIN s2_flat_pairs p USING (pair_id)
            WHERE p.config = ? AND p.kind IN ('bias', 'dark')
              AND p.status = 'ok'""", (config,)).fetchall()
        # Nearest-neighbour pixel correlation of the flat-pair difference
        # images: mean over pairs, standard error from the pair scatter.
        rhos = np.array([r[0] for r in con.execute("""
            SELECT (rho_x + rho_y) / 2.0 FROM s2_flat_pairs
            WHERE config = ? AND kind = 'flat' AND status = 'ok'
              AND rho_x IS NOT NULL AND rho_y IS NOT NULL
              AND max(rho_x, rho_y) <= ?
              AND lower(coalesce(filter, '')) NOT IN
                  ('hrg', 'lrg', 'hagrism', 'oggrism', 'hag', 'ogg')""",
            (config, FLAT_MAX_RHO))])
        rho_nn = float(rhos.mean()) if rhos.size else None
        rho_err = (float(rhos.std(ddof=1) / np.sqrt(rhos.size))
                   if rhos.size > 1 else None)
        gain = gain_stat = gain_sys = None
        basis, status, note = None, "no_measurement", ""
        ad = adopted.get(config)
        if ad and ad["gain"] is not None:
            gain, gain_stat = ad["gain"], ad["stat"]
            note += ad["note"]
            sys2 = ad["sys_estimator"] ** 2
            if rho_nn is not None:
                # Pixel correlation.  Variance that a correlated readout
                # (or charge sharing) moves from a pixel into covariance
                # with its neighbours is missing from the per-pixel
                # variance, so the PTC slope overstates the conversion
                # gain by (1 + 2 rho_x + 2 rho_y) = (1 + 4 rho_nn), rho_nn
                # being the mean of the two lag-1 coefficients.  A
                # SIGNIFICANT rho (> 2 sigma) is corrected for and its
                # uncertainty carried; an insignificant one is carried as
                # a systematic of its own size.  (The first version only
                # carried it — which on the iKon, rho_nn ~ 0.01, meant a
                # 4% error bar around a value known to be 4% high.)
                if rho_err is not None and abs(rho_nn) > 2.0 * rho_err:
                    factor = 1.0 + 4.0 * rho_nn
                    gain, gain_stat = gain / factor, gain_stat / factor
                    sys2 = (sys2 / factor ** 2
                            + (4.0 * rho_err * gain) ** 2)
                    note += (f"pixel correlation rho_nn = {rho_nn:.4f} +/- "
                             f"{rho_err:.4f} corrected for (PTC slope / "
                             f"{factor:.4f}); ")
                else:
                    rho_bound = abs(rho_nn) if rho_err is None else \
                        max(abs(rho_nn), rho_err)
                    sys2 += (4.0 * rho_bound * gain) ** 2
            gain_sys = float(np.sqrt(sys2))
            basis = ad["basis"]
            status = "measured"
        elif sg and corr is not None:
            # Validation ratio: the SAME camera's, when one of its other
            # configurations has flats (what biases a sky pair — pixel
            # scale, star density per pixel, read-noise floor — is a
            # property of the camera); its whole departure from 1 is then
            # carried as the systematic, on top of its own error.  With
            # no same-camera validation: the all-camera mean ratio and the
            # camera-to-camera scatter.
            same = [(r, e) for c, r, e in ratios
                    if cam.camera_name(c) == cam.camera_name(config)]
            if same:
                c_use = float(np.mean([r for r, _e in same]))
                c_sys = float(np.hypot(abs(c_use - 1.0),
                                       np.mean([e for _r, e in same])))
                c_txt = (f"flat/sky ratio {c_use:.4f} measured on the same "
                         f"camera ({len(same)} configuration with flats)")
            else:
                c_use, c_sys = corr, corr_sys
                c_txt = (f"flat/sky ratio {corr:.3f} +/- {corr_sys:.3f} "
                         f"(mean and scatter over {len(ratios)} OTHER "
                         "cameras that have flats)")
            gain = sg["gain"] * c_use
            gain_stat = sg["gain_err"] * c_use
            gain_sys = float(gain * c_sys)
            basis = (f"sky-pair gain (median of {sg['n_pairs']} signal-"
                     f"dominated pairs; zero reference: {sg['basis']}) x "
                     + c_txt)
            status = "measured_sky"
            note += ("no flats on disk for this configuration — a flat-"
                     "pair PTC at the October re-opening is the "
                     "confirmation; ")
            if sg["pseudo"]:
                status = "provisional_sky"
                note += ("no bias or dark frames either: the zero "
                         "reference is a short sky pair; ")
        elif sg:
            note += (f"sky-pair gain {sg['gain']:.3f} recorded but NOT "
                     "adopted (no validation configuration); ")
        gain_err = (float(np.hypot(gain_stat, gain_sys))
                    if gain is not None else None)
        # Read noise from the zero-signal pairs (inverse-variance mean of
        # the pair variances; error = the larger of formal and scatter).
        rn_adu = rn_adu_err = rn_e = rn_e_err = bias = None
        rn_basis = None
        if zero_rows:
            v = np.array([r[0] for r in zero_rows])
            e = np.array([r[1] for r in zero_rows])
            w = 1.0 / e ** 2
            vm = float((w * v).sum() / w.sum())
            formal = float(1.0 / np.sqrt(w.sum()))
            scatter = (float(v.std(ddof=1) / np.sqrt(v.size))
                       if v.size > 1 else formal)
            verr = max(formal, scatter)
            rn_adu = float(np.sqrt(max(vm, 0.0)))
            rn_adu_err = verr / (2.0 * rn_adu) if rn_adu > 0 else None
            bias = float(np.mean([r[2] for r in zero_rows]))
            rn_basis = f"{v.size} zero-signal pairs (bias / shortest dark)"
            if gain is not None:
                rn_e = rn_adu * gain
                rn_e_err = float(np.hypot(rn_adu_err * gain,
                                          rn_adu * gain_err))
        if bias is None and config in zrefs:
            # A camera with no zero-signal frames (QHY600): the level of
            # the zero reference (a pseudo-bias for that camera) stands in,
            # so peak fractions can be formed; flagged in the note.
            bias = zrefs[config]["bias"]
            if zrefs[config].get("pseudo"):
                note += (f"bias level {bias:.1f} ADU from the pseudo-bias "
                         "(short sky pair); ")
        clip = egain_clip.get((mode, round(meta.get("egain") or 0, 3))) \
            if config.startswith("AC4040") else None
        clip = clip or ceilings.get(mode)
        fs_e = fs_e_err = None
        if gain is not None and clip is not None and bias is not None:
            fs_e, fs_e_err = flatptc.full_scale_electrons(
                clip, bias, gain, gain_err, 1.0)
        verdict = b_ratio = b_err = None
        hdr_gain = meta.get("egain")
        if (n_native > 1 and gain is not None and hdr_gain
                and config.startswith("ASI") and 0.1 < hdr_gain < 0.6):
            bv = flatptc.binning_verdict(gain, gain_err, hdr_gain, n_native)
            verdict, b_ratio, b_err = (bv["verdict"], bv["ratio"],
                                       bv["ratio_err"])
        rel = (gain_err / gain) if gain else None
        con.execute("""INSERT OR REPLACE INTO s2_camera_configs
            (config, camera, mode, n_native, header_egain, n_frames,
             first_night, last_night, eras, gain_e_per_adu, gain_err,
             gain_stat_err, gain_sys_err, gain_basis, gain_rel_err,
             meets_target, read_noise_adu, read_noise_adu_err, read_noise_e,
             read_noise_e_err, rn_basis, bias_adu, clip_adu, full_scale_e,
             full_scale_e_err, binning_verdict, binning_ratio,
             binning_ratio_err, rho_nn, rho_nn_err, status, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
                    ?,?,?,?)""",
            (config, cam.camera_name(config), mode, n_native, hdr_gain,
             meta.get("n_frames"), meta.get("first"), meta.get("last"),
             ",".join(sorted(meta.get("eras", []), key=lambda x: float(x))),
             gain, gain_err, gain_stat, gain_sys, basis, rel,
             None if rel is None
             else int(rel <= flatptc.GAIN_TARGET_REL_ERR),
             rn_adu, rn_adu_err, rn_e, rn_e_err, rn_basis, bias, clip,
             fs_e, fs_e_err, verdict, b_ratio, b_err, rho_nn, rho_err,
             status, note.rstrip("; ")))
        # detector_params rows for this configuration.
        if gain is not None:
            put(config, "gain_e_per_adu", gain, gain_err,
                f"{basis}; unc = stat {gain_stat:.4f} (+) sys "
                f"{gain_sys:.4f} in quadrature", "s2_flat_fits")
        if rn_adu is not None:
            put(config, "read_noise_adu", rn_adu, rn_adu_err,
                "tile variance of zero-signal pair differences / 2; unc = "
                "larger of formal and pair-to-pair scatter", rn_basis)
            put(config, "bias_level_adu", bias, None,
                "mean central level of the zero-signal pairs", rn_basis)
        if rn_e is not None:
            put(config, "read_noise_e", rn_e, rn_e_err,
                "RN(ADU) x MEASURED gain; unc propagates both", rn_basis)
        if fs_e is not None:
            put(config, "full_scale_e", fs_e, fs_e_err,
                "(clip - bias) x measured gain: the largest signal the "
                "configuration can record", f"clip {clip:g} ADU")
        if verdict is not None:
            put(config, "binning_gain_ratio", b_ratio, b_err,
                f"measured K / header native EGAIN; 1 = on-camera SUM, "
                f"{n_native} = on-camera AVERAGE -> verdict: {verdict}",
                "s2_camera_configs")

    # ---- StackPro inherits the single-read gain --------------------------
    # A StackPro frame is a SUM of 1-16 High Gain reads of the same pixels
    # through the same amplifier: e-/ADU is unchanged by construction.
    # Its own sky pairs cannot measure it to 3% (sixteen reads of noise sit
    # under every sky level the archive holds), so the configuration takes
    # the same-epoch High Gain value, and the few StackPro sky pairs with
    # signal above their read variance are recorded as the consistency
    # check that the inheritance is not contradicted.
    for epoch in ("e1.054", "e1.057"):
        sp, hg = f"AC4040 StackPro {epoch}", f"AC4040 High Gain {epoch}"
        src = con.execute(
            "SELECT gain_e_per_adu, gain_err, gain_stat_err, gain_sys_err, "
            "status FROM s2_camera_configs WHERE config = ? "
            "AND gain_e_per_adu IS NOT NULL", (hg,)).fetchone()
        tgt = con.execute("SELECT gain_e_per_adu, read_noise_adu, bias_adu, "
                          "clip_adu FROM s2_camera_configs WHERE config = ?",
                          (sp,)).fetchone()
        if not src or not tgt or tgt[0] is not None:
            continue
        ks = [r[0] for r in con.execute(
            "SELECT k_i FROM s2_sky_gains WHERE config = ? AND k_i IS NOT "
            "NULL AND signal_adu >= v0_adu2 AND (reject_reason IS NULL OR "
            "reject_reason = 'read_noise_dominated')", (sp,))]
        chk = (f"{len(ks)} StackPro sky pairs with S >= V0 give K_i "
               f"{min(ks):.2f}-{max(ks):.2f} (median {np.median(ks):.3f})"
               if ks else "no StackPro sky pair has S >= V0")
        g, ge, gst, gsy, st = src
        rn_e = rn_e_err = fs_e = fs_e_err = None
        if tgt[1] is not None:
            rn_e = tgt[1] * g
            rn_e_err = tgt[1] * ge
        con.execute("""UPDATE s2_camera_configs SET gain_e_per_adu = ?,
            gain_err = ?, gain_stat_err = ?, gain_sys_err = ?,
            gain_rel_err = ?, meets_target = ?, gain_basis = ?, status = ?,
            read_noise_e = ?, read_noise_e_err = ?, note = ?
            WHERE config = ?""",
            (g, ge, gst, gsy, ge / g,
             int(ge / g <= flatptc.GAIN_TARGET_REL_ERR),
             f"inherited from {hg} ({st}): a StackPro frame is a sum of "
             "single reads, which leaves e-/ADU unchanged",
             "inherited", rn_e, rn_e_err,
             chk + "; bias, read-noise variance and ceiling scale with the "
             "number of 2 s sub-reads (1-16) — the values in this row are "
             "for 16 (exposures >= 32 s)", sp))
        put(sp, "gain_e_per_adu", g, ge,
            f"inherited from {hg} (sum of single reads); consistency: {chk}",
            "s2_camera_configs, s2_sky_gains")
        if rn_e is not None:
            put(sp, "read_noise_e", rn_e, rn_e_err,
                "RN(ADU) of the 16-read sum x inherited gain", sp)

    # ---- per-MODE rows (the keys existing consumers read) ----------------
    # A mode whose configurations all carry a measured gain gets ONE
    # mode-level value: the single configuration's, or the inverse-
    # variance mean when there are several and they agree; when they
    # disagree by more than 2 sigma the mode-level row is withheld and the
    # disagreement is the finding.
    by_mode: dict[str, list] = {}
    for config, mode, g, ge, rne, rnee, st in con.execute("""
            SELECT config, mode, gain_e_per_adu, gain_err, read_noise_e,
                   read_noise_e_err, status FROM s2_camera_configs
            WHERE gain_e_per_adu IS NOT NULL
              -- only configurations science frames were taken in vote:
              -- the iKon's 1x-preamp flats (2024-04-09/10 engineering)
              -- measure a real but unused configuration
              AND coalesce(n_frames, 0) > 0"""):
        # Sub-populations that are not the mode's standard configuration
        # (the ASI's unbinned and alternate-gain frames) do not vote.
        if config.startswith("ASI") and config != "ASI Mode0 2x2":
            continue
        by_mode.setdefault(mode, []).append((config, g, ge, rne, rnee, st))
    for mode, rows in by_mode.items():
        g = np.array([r[1] for r in rows])
        ge = np.array([r[2] for r in rows])
        w = 1.0 / ge ** 2
        gm = float((w * g).sum() / w.sum())
        gme = float(1.0 / np.sqrt(w.sum()))
        names = ", ".join(r[0] for r in rows)
        if len(rows) > 1:
            z = float(abs(g[0] - g[1]) / np.hypot(ge[0], ge[1])) \
                if len(rows) == 2 else float(
                    np.sqrt((w * (g - gm) ** 2).sum() / (len(rows) - 1)))
            if z > 2.0:
                put(mode, "gain_epoch_disagreement_sigma", z, None,
                    "configurations of this mode disagree on the gain — "
                    "no mode-level gain issued; use the per-configuration "
                    "rows", names)
                continue
            # The epochs agree: the spread between them is carried.
            gme = float(max(gme, g.std(ddof=1) / np.sqrt(len(rows))))
        put(mode, "gain_e_per_adu", gm, gme,
            "measured gain for the mode: inverse-variance mean over its "
            "configuration(s) (flat-pair PTC where flats exist, validated "
            "sky-pair PTC otherwise) — replaces the v1.2 [lower, upper] "
            "bracket", names)
        rne = [(r[3], r[4]) for r in rows if r[3] is not None and r[4]]
        if rne:
            v = np.array([x[0] for x in rne])
            e = np.array([x[1] for x in rne])
            w2 = 1.0 / e ** 2
            put(mode, "read_noise_e", float((w2 * v).sum() / w2.sum()),
                float(max(1.0 / np.sqrt(w2.sum()),
                          v.std(ddof=1) / np.sqrt(v.size)
                          if v.size > 1 else 0.0)),
                "RN(ADU) from zero-signal pairs x MEASURED gain, inverse-"
                "variance mean over the mode's configuration(s); replaces "
                "the v1.2 value whose uncertainty was the half-width of a "
                "gain bracket", names)
def _adopt_flat_gain(con, config: str, fits: dict) -> dict:
    """Apply rule 1 of the adoption rule to one configuration's flat fits.

    Returns ``{"gain", "stat", "sys_estimator", "basis", "note"}`` (gain
    None when the configuration has no usable flat fit) and records the
    measured additive excess in ``detector_params``-ready form in the
    note.
    """
    fw = fits.get((config, "flat-within"))
    fj = fits.get((config, "flat+zero"))
    none = {"gain": None, "stat": None, "sys_estimator": 0.0,
            "basis": None, "note": ""}
    fw_ok = (fw is not None and flatptc.adopted_gain_error(fw) / fw["gain"]
             <= WITHIN_MAX_REL_ERR)
    if not fw_ok and fj is None:
        return none
    zero = con.execute("""
        SELECT q.var_adu2, q.var_err_adu2
        FROM s2_flat_points q JOIN s2_flat_pairs p USING (pair_id)
        WHERE p.config = ? AND p.kind IN ('bias', 'dark')
          AND p.status = 'ok'""", (config,)).fetchall()
    excess = None
    if fw_ok and zero:
        v = np.array([z[0] for z in zero])
        zv = float(v.mean())
        zv_err = (float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1
                  else float(zero[0][1]))
        excess = flatptc.additive_excess(fw["groups"], zv, zv_err)

    def desc(name: str, f: dict) -> str:
        return (f"{name} ({f['n_pairs']} pairs, chi2/dof "
                f"{f['chi2']:.1f}/{f['dof']})")

    # Between-filter consistency of the within-group gain.  A gain does
    # not depend on wavelength, but what leaks past the clip does (stars
    # are brighter against a B sky than an I sky), so filter-to-filter
    # scatter BEYOND the per-filter errors is an empirical bound on what
    # the contamination still does to the slope.  It is added to the
    # systematic as the excess scatter, not divided by sqrt(n).
    per_filter = [(f["gain"], flatptc.adopted_gain_error(f))
                  for (c, est), f in fits.items()
                  if c == config and est.startswith("flat-within:")
                  and f is not None and f["n_pairs"] >= FLAT_MIN_PAIRS
                  and flatptc.adopted_gain_error(f) / f["gain"]
                  <= WITHIN_MAX_REL_ERR]
    filter_sys = 0.0
    note = ""
    if len(per_filter) >= 3:
        g = np.array([x[0] for x in per_filter])
        e = np.array([x[1] for x in per_filter])
        filter_sys = float(np.sqrt(max(g.var(ddof=1) - np.mean(e ** 2),
                                       0.0)))
        note += (f"per-filter within-group gains ({len(per_filter)} "
                 f"filters): {g.min():.4f}-{g.max():.4f}, excess scatter "
                 f"{filter_sys:.4f}; ")
    if excess is not None:
        note += (f"additive variance of the flat groups over read noise: "
                 f"{excess['excess']:+.2f} +/- {excess['excess_err']:.2f} "
                 f"ADU^2 ({excess['n_groups']} groups); ")
    if fw_ok and (fj is None or (excess is not None and excess["z"]
                                 is not None
                                 and excess["z"] > EXCESS_SIGNIFICANCE)):
        if fj is not None:
            note += (f"pinned fit gives {fj['gain']:.4f} "
                     f"({100 * (fj['gain'] / fw['gain'] - 1):+.2f}%), "
                     "biased by that additive term; ")
        clear = (fj is None or (excess is not None
                                and excess["z"] > EXCESS_CLEAR))
        est_sys = 0.0 if clear else abs(fw["gain"] - fj["gain"])
        return {"gain": fw["gain"], "stat": flatptc.adopted_gain_error(fw),
                "sys_estimator": float(np.hypot(filter_sys, est_sys)),
                "basis": "flat-pair PTC, " + desc(
                    "free additive term per flat group", fw),
                "note": note}
    if fw_ok and fj is not None:
        ew, ej = (flatptc.adopted_gain_error(fw),
                  flatptc.adopted_gain_error(fj))
        best, name = ((fw, "free additive term per flat group")
                      if ew <= ej else
                      (fj, "joint fit with zero-signal pairs"))
        return {"gain": best["gain"], "stat": min(ew, ej),
                "sys_estimator": float(np.hypot(
                    abs(fw["gain"] - fj["gain"]), filter_sys)),
                "basis": "flat-pair PTC, " + desc(name, best),
                "note": note + (f"within-group {fw['gain']:.4f} vs pinned "
                                f"{fj['gain']:.4f}; ")}
    note += ("flats carry no within-group signal spread — an additive "
             "term in them could not be tested; ")
    return {"gain": fj["gain"], "stat": flatptc.adopted_gain_error(fj),
            "sys_estimator": 0.0,
            "basis": "flat-pair PTC, " + desc(
                "joint fit with zero-signal pairs", fj),
            "note": note}


# ---------------------------------------------------------------------------
# params helpers: linearity curves -> ONE cap per mode
# ---------------------------------------------------------------------------
#: Injection trials on each dedicated star set.
STARLIN_INJECTION_TRIALS = 200

#: A recommended cap is quoted rounded DOWN to this many ADU (a number a
#: human can put in an observing plan); small-scale modes round finer.
CAP_GRANULARITY_ADU = 100
CAP_GRANULARITY_SMALL_ADU = 50


def _mode_scale(con) -> dict[str, dict]:
    """Per mode: bias level and ceiling(s) the peak fractions are built on.

    ``bias`` is the mean zero-signal level of the mode's configurations
    (``s2_camera_configs``), falling back to the v1.2 dark-floor pedestal
    (``detector_params.bias_offset_adu``).  ``clip`` is the mode's
    measured ceiling; ``clip_by_egain`` holds the per-EGAIN-epoch median
    frame maximum for the AC4040, whose two epochs clip ~46 ADU apart.
    """
    out: dict[str, dict] = {}
    for mode, clip in con.execute(
            "SELECT mode, clip_adu FROM s2_ceiling_modes "
            "WHERE clip_adu IS NOT NULL"):
        out[mode] = {"clip": float(clip), "bias": None, "clip_by_egain": {}}
    for mode, bias in con.execute(
            "SELECT mode, avg(bias_adu) FROM s2_camera_configs "
            "WHERE bias_adu IS NOT NULL GROUP BY mode"):
        if mode in out:
            out[mode]["bias"] = float(bias)
    for mode, bias in con.execute(
            "SELECT era_group, value FROM detector_params "
            "WHERE quantity = 'bias_offset_adu'"):
        if mode in out and out[mode]["bias"] is None:
            out[mode]["bias"] = float(bias)
    for mode, eg, med in con.execute(
            "SELECT mode, egain, median_max_adu FROM s2_ceiling_egain"):
        if mode in out and mode.startswith("High Gain"):
            out[mode]["clip_by_egain"][round(eg, 3)] = float(med)
    return out


def _frame_scale(scale: dict[str, dict], mode: str, egain,
                 exptime) -> tuple[float, float] | None:
    """(bias, ceiling) in raw ADU for one frame, or None if unknown.

    Single-read modes: the mode's bias and its ceiling (per EGAIN epoch on
    the AC4040).  StackPro: N x the single-read High Gain bias and N x the
    single-read High Gain clip of the same epoch, with N the number of
    2-second sub-reads at this exposure
    (:func:`rlmt_diagnostics.ptc.stackpro_nsub_for_exptime`) — a 16 s
    StackPro frame saturates at 8 x 3.5 kADU, not at the 56 kADU of the
    32 s-and-longer frames the mode's pooled ceiling describes.
    """
    if mode == "High Gain StackPro":
        hg = scale.get("High Gain")
        if not hg or hg["bias"] is None:
            return None
        n = ptc.stackpro_nsub_for_exptime(exptime)
        clip1 = hg["clip_by_egain"].get(round(egain or 0, 3), hg["clip"])
        return n * hg["bias"], n * clip1
    sc = scale.get(mode)
    if not sc or sc["bias"] is None:
        return None
    return sc["bias"], sc["clip_by_egain"].get(round(egain or 0, 3),
                                              sc["clip"])


def _params_starlin(con, scale: dict[str, dict]) -> None:
    """Turn the stored star measurements into deviation points and curves."""
    con.execute("DELETE FROM s2_linearity_points")
    con.execute("DELETE FROM s2_linearity_curve WHERE source NOT LIKE 'cv%' "
                "AND source != 'combined'")
    con.execute("DELETE FROM s2_linearity_injection "
                "WHERE source NOT LIKE 'cv%'")
    sets = [r[0] for r in con.execute(
        "SELECT DISTINCT set_id FROM s2_starlin_frames WHERE status = 'ok'")]
    pooled: dict[tuple, dict] = {}
    biggest: dict[str, tuple[int, dict]] = {}
    for set_no, set_id in enumerate(sorted(sets)):
        family = set_id.split("|")[0]
        frames_ = con.execute("""
            SELECT path, mode, egain, exptime, fwhm_med
            FROM s2_starlin_frames WHERE set_id = ? AND status = 'ok'
            ORDER BY jd""", (set_id,)).fetchall()
        fw = [f[4] for f in frames_ if f[4]]
        if len(frames_) < 2 or not fw:
            continue
        aperture = starphot.choose_aperture(max(fw))
        fidx = {f[0]: j for j, f in enumerate(frames_)}
        meas = con.execute("""
            SELECT path, star_id, flux, flux_err, peak_raw
            FROM s2_starlin_meas WHERE set_id = ? AND aperture_px = ?""",
            (set_id, aperture)).fetchall()
        if not meas:
            continue
        n_stars = max(m[1] for m in meas) + 1
        shape = (n_stars, len(frames_))
        flux = np.full(shape, np.nan)
        ferr = np.full(shape, np.nan)
        peak = np.full(shape, np.nan)
        for path, sid, fl, fe, pk in meas:
            if path in fidx:
                flux[sid, fidx[path]] = fl
                ferr[sid, fidx[path]] = fe
                peak[sid, fidx[path]] = pk
        pfrac = np.full(shape, np.nan)
        for j, (_p, mode, eg, t_exp, _f) in enumerate(frames_):
            fs_ = _frame_scale(scale, mode, eg, t_exp)
            if fs_ is None:
                continue
            pfrac[:, j] = lin.peak_fraction(peak[:, j], fs_[0], fs_[1])
        expt = np.array([f[3] for f in frames_], dtype=np.float64)
        flux = np.where(flux > 0, flux, np.nan)
        res = lin.ensemble_deviation(flux, ferr, expt, pfrac)
        rows = []
        for j, (path, mode, _eg, _t, _f) in enumerate(frames_):
            for i in range(n_stars):
                if np.isfinite(res["dev"][i, j]) and np.isfinite(pfrac[i, j]):
                    rows.append((family, set_id, mode, path, i,
                                 float(peak[i, j]), float(pfrac[i, j]),
                                 float(res["dev"][i, j]),
                                 float(res["dev_err"][i, j]),
                                 int(res["is_ref"][i, j])))
                    acc = pooled.setdefault((family, mode),
                                            {"pf": [], "dev": [], "grp": []})
                    acc["pf"].append(pfrac[i, j])
                    acc["dev"].append(res["dev"][i, j])
                    acc["grp"].append(set_no * 1_000_000 + i)
        con.executemany("INSERT INTO s2_linearity_points VALUES "
                        "(?,?,?,?,?,?,?,?,?,?)", rows)
        n_above = sum(1 for r in rows if r[6] >= lin.REF_PEAK_FRACTION)
        if family not in biggest or n_above > biggest[family][0]:
            with np.errstate(invalid="ignore", divide="ignore"):
                model = res["base"][:, None] * res["z"][None, :] * expt[None, :]
                ppf = pfrac / model
            biggest[family] = (n_above, {
                "flux_err": ferr, "exptime": expt, "base": res["base"],
                "z": res["z"], "ppf": ppf,
                "modes": [f[1] for f in frames_]})
    for (family, mode), acc in sorted(pooled.items()):
        curve = lin.deviation_curve(np.array(acc["pf"]),
                                    np.array(acc["dev"]),
                                    np.array(acc["grp"]))
        _store_curve(con, family, mode, curve)
    # Signed injection bias on each family's richest set (rule 3).
    for family, (_n, geo) in sorted(biggest.items()):
        with lin._quiet_nan_warnings():
            ppf_star = np.nanmedian(geo["ppf"], axis=1)
            err_star = np.nanmedian(geo["flux_err"], axis=1)
        good = np.isfinite(geo["base"]) & np.isfinite(ppf_star)
        if int(good.sum()) < 5:
            continue
        ppf = np.where(np.isfinite(geo["ppf"]), geo["ppf"],
                       ppf_star[:, None])[good]
        ee = np.where(np.isfinite(geo["flux_err"]), geo["flux_err"],
                      err_star[:, None])[good]
        zz = np.where(np.isfinite(geo["z"]), geo["z"], 1.0)
        for case, fn in (("null", lambda pf: 0.0 * pf),
                         ("rolloff_5pct", rolloff_5pct)):
            rows = lin.injection_bias(ee, geo["exptime"], geo["base"][good],
                                      zz, ppf, fn,
                                      n_trials=STARLIN_INJECTION_TRIALS)
            # The simulation pools every frame of the set; it is filed
            # under the family's tested (single-read) mode.
            tested = next((m for m in geo["modes"] if "StackPro" not in m),
                          geo["modes"][0])
            _store_injection(con, family, tested, case, rows)


#: No cap is issued above this fraction of the scale, whatever a curve
#: says: within 5% of the ceiling "peak fraction" is itself uncertain (the
#: AC4040 clip is a per-pixel mound tens of ADU wide; a reduced frame's
#: mapping back to raw carries the flat-field's percent-level structure).
CAP_MAX_FRACTION = 0.95


def _params_caps(con, put, scale: dict[str, dict]) -> None:
    """Judge every source's curve on its own and issue ONE cap per mode.

    THE RULE.  Each independent source (a dedicated star set; the pooled
    CV series) gets its own cap from
    :func:`rlmt_diagnostics.linearity.recommend_cap`.  Then:

    * if ANY source's cap is set by measured non-linearity, the mode's cap
      is the LOWEST such cap — a cap has to hold in every data set that
      tested it; averaging a source that sees a roll-off with one that
      does not would hide the roll-off inside a wider error bar (the first
      version of this function did exactly that, by cap-walking an
      inverse-variance COMBINED curve);
    * otherwise the cap is the HIGHEST level any source measured linear to
      ("linear as far as measured");
    * never above :data:`CAP_MAX_FRACTION`.

    The combined curve is still stored (source ``combined``) because it is
    the right thing to PLOT; it decides nothing.
    """
    con.execute("DELETE FROM s2_linearity_curve WHERE source = 'combined'")
    con.execute("DELETE FROM s2_linearity_caps")
    con.execute("DELETE FROM s2_linearity_source_caps")
    for mode, sc in sorted(scale.items()):
        rows = con.execute("""
            SELECT source, lo, hi, peak_frac, dev_pct, dev_err_pct,
                   n_points, n_groups, measured
            FROM s2_linearity_curve
            WHERE mode = ? AND source NOT LIKE 'cv:%' AND source != 'combined'
            ORDER BY source, lo""", (mode,)).fetchall()
        bias, clip = sc["bias"], sc["clip"]
        if mode == "High Gain StackPro":
            # Peak fractions of StackPro frames were built per frame on
            # N x the single-read scale; quote the cap on the 16-read one.
            hg = scale.get("High Gain")
            if hg and hg["bias"] is not None:
                n = ptc.STACKPRO_MAX_NSUB
                bias, clip = n * hg["bias"], n * hg["clip"]
        by_src: dict[str, list[dict]] = {}
        for src, lo, hi, pf, dv, de, npnt, ngrp, meas in rows:
            by_src.setdefault(src, []).append({
                "lo": lo, "hi": hi, "peak_frac": pf, "dev_pct": dv,
                "dev_err_pct": de if de is not None else float("nan"),
                "n_points": npnt, "n_groups": ngrp, "measured": bool(meas)})
        if not by_src or bias is None:
            con.execute(
                "INSERT OR REPLACE INTO s2_linearity_caps "
                "(mode, bias_adu, ceiling_adu, sources, status, note) "
                "VALUES (?,?,?,?,?,?)",
                (mode, bias, clip, "", "no_data",
                 "no peak-resolved linearity measurement exists for this "
                 "mode in the archive"))
            continue
        recs = {}
        for src, curve in sorted(by_src.items()):
            # CONTROL: the reference-regime bins must read zero.  They are
            # zero "by construction" only if the normalisation worked; a
            # source whose measured bins between 0.1 and the reference
            # limit depart from zero by more than the criterion has a
            # broken normalisation and may not issue or lower a cap.
            ctrl = [c for c in curve if c["measured"] and c["lo"] >= 0.1
                    and c["hi"] <= lin.REF_PEAK_FRACTION + 1e-9]
            if any(abs(c["dev_pct"]) > lin.CAP_CRITERION_PCT for c in ctrl):
                con.execute(
                    "INSERT OR REPLACE INTO s2_linearity_source_caps "
                    "(mode, source, limited_by, worst_dev_pct) "
                    "VALUES (?,?,?,?)",
                    (mode, src, "failed_reference_control",
                     max((c["dev_pct"] for c in ctrl), key=abs)))
                continue
            r = lin.recommend_cap(curve)
            recs[src] = r
            fb = r["first_bad"]
            con.execute(
                "INSERT OR REPLACE INTO s2_linearity_source_caps VALUES "
                "(?,?,?,?,?,?,?,?,?,?)",
                (mode, src, r["cap_fraction"], r["limited_by"],
                 r["worst_dev_pct"], r["precision_pct"],
                 fb["lo"] if fb else None, fb["dev_pct"] if fb else None,
                 fb["dev_err_pct"] if fb and np.isfinite(fb["dev_err_pct"])
                 else None,
                 max((c["hi"] for c in curve if c["measured"]),
                     default=None)))
        # Combined curve, for the figure only.
        by_bin: dict[tuple, list] = {}
        for src, curve in by_src.items():
            for c in curve:
                if c["measured"] and np.isfinite(c["dev_err_pct"]) \
                        and c["dev_err_pct"] > 0:
                    by_bin.setdefault((c["lo"], c["hi"]), []).append(c)
        combined = []
        for (lo, hi), items in sorted(by_bin.items()):
            dv = np.array([i["dev_pct"] for i in items])
            w = 1.0 / np.array([i["dev_err_pct"] for i in items]) ** 2
            err = float(1.0 / np.sqrt(w.sum()))
            if len(items) > 1:
                err = float(max(err, dv.std(ddof=1) / np.sqrt(len(items))))
            combined.append({
                "lo": lo, "hi": hi,
                "peak_frac": float(np.average(
                    [i["peak_frac"] for i in items], weights=w)),
                "dev_pct": float((w * dv).sum() / w.sum()),
                "dev_err_pct": err,
                "n_points": int(sum(i["n_points"] for i in items)),
                "n_groups": int(sum(i["n_groups"] for i in items)),
                "measured": True})
        _store_curve(con, "combined", mode, combined)
        if not recs:
            con.execute(
                "INSERT OR REPLACE INTO s2_linearity_caps "
                "(mode, bias_adu, ceiling_adu, sources, status, note) "
                "VALUES (?,?,?,?,?,?)",
                (mode, bias, clip, ", ".join(sorted(by_src)), "no_data",
                 "every source failed its reference-regime control"))
            continue
        nonlin = {s_: r for s_, r in recs.items()
                  if r["limited_by"] == "nonlinearity"}
        if nonlin:
            deciding = min(nonlin, key=lambda s_: (
                nonlin[s_]["cap_fraction"]
                if nonlin[s_]["cap_fraction"] is not None else 0.0))
            status = "cap_set_by_measured_nonlinearity"
        else:
            with_cap = {s_: r for s_, r in recs.items()
                        if r["cap_fraction"] is not None}
            if not with_cap:
                con.execute(
                    "INSERT OR REPLACE INTO s2_linearity_caps "
                    "(mode, bias_adu, ceiling_adu, sources, status, note) "
                    "VALUES (?,?,?,?,?,?)",
                    (mode, bias, clip, ", ".join(sorted(recs)), "no_data",
                     "no bin above the reference regime was measured"))
                continue
            deciding = max(with_cap,
                           key=lambda s_: with_cap[s_]["cap_fraction"])
            status = "linear_as_far_as_measured"
        rec_ = recs[deciding]
        cap_f = rec_["cap_fraction"]
        if cap_f is None:
            status = "nonlinear_from_first_tested_bin"
        elif cap_f > CAP_MAX_FRACTION:
            cap_f = CAP_MAX_FRACTION
            status = "linear_to_cap_maximum"
        cap_adu = None
        if cap_f is not None:
            gran = (CAP_GRANULARITY_SMALL_ADU if clip < 10000
                    else CAP_GRANULARITY_ADU)
            cap_adu = float(np.floor((bias + cap_f * (clip - bias)) / gran)
                            * gran)
        fb = rec_["first_bad"]
        src_txt = ", ".join(sorted(recs))
        con.execute("""INSERT OR REPLACE INTO s2_linearity_caps
            (mode, bias_adu, ceiling_adu, cap_fraction, cap_adu, limited_by,
             worst_dev_pct, precision_pct, slope_pct_per_scale, slope_err,
             first_bad_lo, first_bad_dev_pct, measured_to_fraction, sources,
             status, note) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (mode, bias, clip, cap_f, cap_adu, rec_["limited_by"],
             rec_["worst_dev_pct"], rec_["precision_pct"],
             rec_["slope_pct_per_scale"], rec_["slope_err"],
             fb["lo"] if fb else None, fb["dev_pct"] if fb else None,
             max((c["hi"] for c in by_src[deciding] if c["measured"]),
                 default=None),
             src_txt, status,
             f"deciding source: {deciding}; criterion |deviation| <= "
             f"{lin.CAP_CRITERION_PCT:g}%; reference regime peak < "
             f"{lin.REF_PEAK_FRACTION:g} of scale"))
        if cap_adu is not None:
            put(mode, "linearity_cap_adu", cap_adu, None,
                f"highest peak level (raw ADU) below which the measured "
                f"response deviation stays within "
                f"{lin.CAP_CRITERION_PCT:g}% in every source that tested "
                f"it ({status.replace('_', ' ')}; deciding source "
                f"{deciding}); = bias + {cap_f:.2f} x (ceiling - bias), "
                "rounded down", src_txt)
            put(mode, "linearity_cap_fraction", cap_f, None,
                "the same cap as a fraction of the usable scale "
                "(ceiling - bias); exact bin edge — no uncertainty",
                src_txt)
        if rec_["slope_pct_per_scale"] is not None:
            put(mode, "linearity_slope_pct_per_scale",
                rec_["slope_pct_per_scale"], rec_["slope_err"],
                "weighted straight-line slope of response deviation (%) "
                "against peak fraction over the deciding source's bins up "
                "to the cap", deciding)


def _source_epoch(con, source: str, mode: str) -> float | None:
    """The EGAIN epoch (rounded to 3 decimals) a linearity source belongs to.

    Dedicated star sets: the majority EGAIN of their frames of that mode.
    The CV source: the majority EGAIN of the eras of its series of that
    mode.  None when the frames carry no EGAIN.
    """
    if source == "cv":
        rows = con.execute("""
            SELECT round(e.egain, 3), sum(s.n_points)
            FROM s2_peaklin_series s JOIN eras e ON e.era_id = s.era_id
            WHERE s.mode = ? AND s.status = 'ok' AND e.egain > 0
            GROUP BY 1 ORDER BY 2 DESC""", (mode,)).fetchall()
    else:
        rows = con.execute("""
            SELECT round(egain, 3), count(*) FROM s2_starlin_frames
            WHERE family = ? AND mode = ? AND status = 'ok' AND egain > 0
            GROUP BY 1 ORDER BY 2 DESC""", (source, mode)).fetchall()
    return float(rows[0][0]) if rows else None


def _params_epoch_caps(con, put, scale: dict[str, dict]) -> None:
    """High Gain caps per EGAIN epoch (1.054: 2023-02 to 07; 1.057 after).

    The two epochs clip 46 ADU apart (s2_ceiling_egain), so they are two
    configurations (DE.F9) and a cap measured in one does not license the
    other.  Each source is assigned to the epoch of its frames
    (:func:`_source_epoch`) and the per-mode rule of :func:`_params_caps`
    is applied within each epoch.  The mode-level High Gain cap remains
    the minimum over both (any consumer that does not know the epoch gets
    the safe value).
    """
    sc = scale.get("High Gain")
    if not sc or sc["bias"] is None:
        return
    by_epoch: dict[float, list] = {}
    for src, cap_f, limited, worst, prec in con.execute("""
            SELECT source, cap_fraction, limited_by, worst_dev_pct,
                   precision_pct FROM s2_linearity_source_caps
            WHERE mode = 'High Gain'
              AND limited_by != 'failed_reference_control'""").fetchall():
        eg = _source_epoch(con, src, "High Gain")
        if eg is not None:
            by_epoch.setdefault(eg, []).append((src, cap_f, limited, worst,
                                                prec))
    for eg, rows in sorted(by_epoch.items()):
        nonlin = [r for r in rows if r[2] == "nonlinearity"]
        if nonlin:
            dec = min(nonlin, key=lambda r: r[1] if r[1] is not None else 0)
            status = "cap_set_by_measured_nonlinearity"
        else:
            capped = [r for r in rows if r[1] is not None]
            if not capped:
                continue
            dec = max(capped, key=lambda r: r[1])
            status = "linear_as_far_as_measured"
        cap_f = min(dec[1], CAP_MAX_FRACTION) if dec[1] is not None else None
        clip = sc["clip_by_egain"].get(round(eg, 3), sc["clip"])
        key = f"AC4040 High Gain e{eg:.3f}"
        cap_adu = (float(np.floor((sc["bias"] + cap_f * (clip - sc["bias"]))
                                  / CAP_GRANULARITY_SMALL_ADU)
                         * CAP_GRANULARITY_SMALL_ADU)
                   if cap_f is not None else None)
        srcs = ", ".join(sorted(r[0] for r in rows))
        con.execute("""INSERT OR REPLACE INTO s2_linearity_caps
            (mode, bias_adu, ceiling_adu, cap_fraction, cap_adu, limited_by,
             worst_dev_pct, precision_pct, sources, status, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
            (key, sc["bias"], clip, cap_f, cap_adu, dec[2], dec[3], dec[4],
             srcs, status, f"deciding source: {dec[0]}; EGAIN epoch "
             f"{eg:.3f}; ceiling = that epoch's median saturated frame "
             "maximum (s2_ceiling_egain)"))
        if cap_adu is not None:
            put(key, "linearity_cap_adu", cap_adu, None,
                f"per-EGAIN-epoch High Gain cap ({status.replace('_', ' ')};"
                f" deciding source {dec[0]}) = bias + {cap_f:.2f} x "
                f"(epoch ceiling {clip:.0f} - bias), rounded down", srcs)
            put(key, "linearity_cap_fraction", cap_f, None,
                "the same cap as a fraction of the epoch's usable scale "
                "(exact bin edge — no uncertainty)", srcs)


def _params_targets(con, scale: dict[str, dict]) -> None:
    """Issue a saturation verdict for every censused target frame."""
    con.execute("DELETE FROM s2_target_verdicts")
    caps = dict(con.execute(
        "SELECT mode, cap_fraction FROM s2_linearity_caps"))
    egain_of = dict(con.execute(
        "SELECT config, header_egain FROM s2_camera_configs"))
    for (rowid, tkey, mode, filt, n_bin, peak, n_at, fwhm, disp,
         status, config, t_exp) in con.execute("""
            SELECT obs_rowid, target_key, mode, filter, n_bin, peak_raw,
                   n_at_peak, fwhm_px, dispersion_verdict, status, config,
                   exptime
            FROM s2_target_peaks""").fetchall():
        fs_ = _frame_scale(scale, mode, egain_of.get(config), t_exp)
        row = [rowid, tkey, mode, filt, None, None, None, None, None,
               caps.get(mode), None]
        if status != "ok":
            row[-1] = status
        elif disp == "dispersed":
            # A grism frame has no direct image of the target: the "peak"
            # found is a piece of spectrum, and no photometric verdict is
            # meaningful.  Recorded, not judged.
            row[-1] = "dispersed"
        elif fs_ is None:
            # No measured ceiling for this mode (Low Gain): only a flat
            # top can be recognised.
            row[-1] = ("flat_topped"
                       if (n_at or 0) >= satmod.FLAT_TOP_MIN_PIXELS
                       else "no_ceiling_for_mode")
        else:
            bias, clip = fs_
            native, factor = satmod.native_equivalent_peak(
                peak, bias, fwhm, n_bin or 1)
            frac = (native - bias) / (clip - bias)
            row[4:9] = [bias, clip, factor, native, frac]
            row[-1] = satmod.saturation_verdict(
                native, bias, clip, caps.get(mode), n_at or 0)
        con.execute("INSERT OR REPLACE INTO s2_target_verdicts VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?)", row)


def _params_badpix(con, put) -> None:
    """One detector_params row per bad-pixel mask product."""
    for (mkey, n_bad, frac, n_sci, n_hot, n_rail, n_noisy, npz) in \
            con.execute("""SELECT mask_key, n_bad, bad_fraction, n_science,
                                  n_hot, n_rail, n_noisy, npz_path
                           FROM s2_badpix_masks""").fetchall():
        put(f"mask {mkey}", "bad_pixel_count", n_bad, None,
            f"pixels firing the neighbour-excess (> {badpix.HOT_SIGMA:g} "
            f"sigma) or rail test in >= "
            f"{100 * badpix.PERSIST_FRACTION:.0f}% of {n_sci} science "
            f"frames from different fields, plus temporally noisy pixels "
            f"from a dark stack (hot {n_hot}, rail {n_rail}, noisy "
            f"{n_noisy}); exact count — no uncertainty", npz)
        put(f"mask {mkey}", "bad_pixel_fraction", frac, None,
            "bad pixels / all pixels (exact ratio — no uncertainty)", npz)


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
                "quiet hot pixels) — S2 v1.2 bracket, SUPERSEDED by the "
                "flat-pair gain_e_per_adu", f"{fd['n_points']} dark points")
        lpts = con.execute(
            "SELECT level, var, n_pix FROM s2_ptc_points WHERE mode = ? "
            "AND kind = 'light' AND level < ?",
            (mode, ptc.GAIN_LOWER_BOUND_LEVEL_ADU)).fetchall()
        fl = ptc.fit_ptc([p[0] for p in lpts], [p[1] for p in lpts],
                         [p[2] for p in lpts])
        if fl:
            put(mode, "gain_lower_bound_e_per_adu", fl["gain_e_per_adu"],
                fl["gain_err"], "sky-level light-pair PTC slope (motion-"
                "inflated variance) — S2 v1.2 bracket, SUPERSEDED by the "
                "flat-pair gain_e_per_adu", f"{fl['n_points']} light points "
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

    # --- flat-pair / sky-pair photon transfer (review F-4, D2) -------------
    # Runs AFTER the v1.2 PTC block on purpose: the mode-level
    # ``read_noise_e`` it writes replaces the bracket-based value above
    # (same key), while the v1.2 bracket rows stay as the record of what
    # the dark/star-field method could and could not do.
    _params_flat(con, put)

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
    # --- peak-resolved linearity, ONE cap per mode (review F-5) ------------
    scale = _mode_scale(con)
    _params_starlin(con, scale)
    _params_caps(con, put, scale)
    _params_epoch_caps(con, put, scale)
    # The pyscope-era QHY frames (blank READOUTM) are the same sensor at
    # the same GAIN 56 as the MaxIm-era 'Fast' frames — only OFFSET
    # differs, which moves the bias and not the response — and carry no
    # linearity data of their own: they take the Fast cap, labelled so.
    fast = con.execute("SELECT cap_fraction, worst_dev_pct, precision_pct, "
                       "sources FROM s2_linearity_caps WHERE mode = 'Fast' "
                       "AND cap_fraction IS NOT NULL").fetchone()
    blank = cam.mode_label_of("QHY600 pyscope 2x2")
    sc_b = scale.get(blank)
    if fast and sc_b and sc_b["bias"] is not None:
        cap_adu = float(np.floor((sc_b["bias"] + fast[0] * (
            sc_b["clip"] - sc_b["bias"])) / CAP_GRANULARITY_ADU)
            * CAP_GRANULARITY_ADU)
        con.execute("""INSERT OR REPLACE INTO s2_linearity_caps
            (mode, bias_adu, ceiling_adu, cap_fraction, cap_adu, limited_by,
             worst_dev_pct, precision_pct, sources, status, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
            (blank, sc_b["bias"], sc_b["clip"], fast[0], cap_adu, "data",
             fast[1], fast[2], fast[3], "inherited_from_Fast",
             "same QHY600 sensor and GAIN 56; no linearity data of its own"))
        put(blank, "linearity_cap_adu", cap_adu, None,
            "inherited from mode Fast (same sensor and gain setting)",
            fast[3])
        put(blank, "linearity_cap_fraction", fast[0], None,
            "inherited from mode Fast (exact bin edge)", fast[3])
    # --- per-frame target saturation verdicts (review OA.E3, rule 4) ------
    _params_targets(con, scale)
    # --- bad-pixel mask products ------------------------------------------
    _params_badpix(con, put)
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
    p.add_argument("stage", choices=["seed", "ceiling", "ceilpos", "ptc",
                                     "noise", "flats", "skypairs",
                                     "starlin", "peaklin",
                                     "peakcensus", "badpix",
                                     "reconstruct", "linearity", "params",
                                     "report", "promote"])
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="S0/S0b manifest database (read-only unless it is "
                        "also --db, or the stage is 'promote')")
    p.add_argument("--db", type=Path, default=DEFAULT_DETECTOR_DB,
                   help="S2 work database: every s2_* table and "
                        "detector_params are written here")
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE,
                   help="archive root holding the pixel trees")
    p.add_argument("--batch", type=int, default=250,
                   help="max frames/pairs/ladders processed this invocation")
    p.add_argument("--era", type=int, default=None,
                   help="reconstruct: process only this era")
    p.add_argument("--only", default=None,
                   help="flats/skypairs/badpix: restrict this invocation to "
                        "configurations (or mask keys) containing this text")
    p.add_argument("--cv-db", type=Path, default=DEFAULT_CV_DB,
                   help="peaklin: CV photometry database (read-only)")
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
    if args.stage == "promote":
        return cmd_promote(args.db, args.manifest)
    with closing(open_work_db(args.db, args.manifest)) as con:
        ensure_tables(con)
        if args.stage == "seed":
            return cmd_seed(con)
        write_meta(con, args.manifest)
        if args.stage == "flats":
            return cmd_flats(con, args.archive, args.batch, args.only)
        if args.stage == "skypairs":
            return cmd_skypairs(con, args.archive, args.batch, args.only)
        if args.stage == "starlin":
            return cmd_starlin(con, args.archive, args.batch)
        if args.stage == "peaklin":
            return cmd_peaklin(con, args.cv_db)
        if args.stage == "peakcensus":
            return cmd_peakcensus(con, args.archive, args.batch)
        if args.stage == "badpix":
            return cmd_badpix(con, args.archive, args.batch, args.only)
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
