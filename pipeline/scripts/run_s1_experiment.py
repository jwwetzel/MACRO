#!/usr/bin/env python
"""Run the S1 stratified astrometry go/no-go experiment.

WHAT THIS SCRIPT DOES (stage S1 of the shared pipeline)
-------------------------------------------------------
Decides, with evidence, whether batch re-solving the unsolved backlog is
viable — per stratum — before anyone spends days of compute on ~59k frames.
The CV paper's polars (ST LMi / VV Pup / EU UMa / AN UMa Sloan series) are
73–95% unsolved; this experiment samples every homogeneous (camera family ×
exposure band × project) cell, solves the samples with a local
astrometry.net, and records the outcomes in the manifest so the report and
the Week-2 review argue from a table, not an anecdote.

It AUGMENTS the manifest with NEW tables only (S0/S0b tables untouched):

* ``s1_strata``           — the design: definition, population, sample, seed
* ``s1_populations``      — the candidate-universe census (what was excluded
                            and why: measured spectra, sub-frame windows)
* ``s1_gate_comparison``  — the FILTER-label universe cross-tabbed against
                            the measured-dispersion universe, per frame
                            class: what the gate correction moves, and which
                            way
* ``s1_solve_experiment`` — one row per sampled frame with the solve outcome
* ``s1_failure_autopsy``  — image-statistics post-mortem of failures
* ``s1_build_meta``       — timestamp, versions, tooling inventory

and, once ``snapshot`` has been run, a frozen copy of a previous experiment
under ``s1_baseline_*`` so the report can render a before/after delta from
the database rather than from anybody's memory.

SUBCOMMANDS (run in this order)
-------------------------------
    snapshot  freeze the CURRENT experiment as the comparison baseline
              (run BEFORE a design change; refuses to overwrite an
              existing baseline without --force)
    design    build strata, draw the reproducible samples (idempotent:
              re-running replaces the s1 tables wholesale)
    run       solve pending frames until --batch-seconds elapses; SAFE TO
              RE-RUN — it resumes where it stopped (each invocation is one
              resumable batch, sized for the 10-minute shell cap)
    autopsy   image-statistics post-mortem + thumbnails for failures
    report    render docs/pipeline/s1_astrometry.html from the database
    status    one-line-per-stratum progress table

USAGE (a student's quick start)
-------------------------------
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/run_s1_experiment.py snapshot      # only when the
                                                            # design changes
    $PY pipeline/scripts/run_s1_experiment.py design
    $PY pipeline/scripts/run_s1_experiment.py run --batch-seconds 420
    ...repeat run until status shows no pending...
    $PY pipeline/scripts/run_s1_experiment.py autopsy
    $PY pipeline/scripts/run_s1_experiment.py report

The experiment NEVER batch-solves beyond its samples: the batch decision
belongs to the Week-2 review, reading the S1 report.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import time
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

# Make the pipeline package importable no matter where the script is invoked
# from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import astrom                                 # noqa: E402
from macro_core.astrom import S1_CODE_VERSION                 # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
DEFAULT_INDEX_DIR = Path.home() / "astrometry-indices"
DEFAULT_CONFIG = DEFAULT_INDEX_DIR / "astrometry.cfg"
PRODUCTS_DIR = REPO_ROOT / "products" / "s1"
THUMB_DIR = PRODUCTS_DIR / "thumbs"

#: Saturation rail per readout family for the autopsy, in ADU.  The raw
#: trees are 16-bit; the GSENSE High Gain mode clips near 3.5 kADU (the SN
#: team's measurement, S2 will confirm) — everything else rails at 65535.
GSENSE_HIGHGAIN_SAT_ADU = 3500.0
DEFAULT_SAT_ADU = 65000.0

#: SQLite lock patience, in milliseconds — the project-wide standard (see
#: rescan_geometry.py, build_s3_timing.py, run_s1_batch.py).  This stage
#: does NOT own the database: the S1b production batch commits into
#: ``s1_batch`` in bursts while this experiment runs, and ``design``'s
#: table swap is a single exclusive transaction that would otherwise die
#: on the driver's 5-second default the moment the batch happened to be
#: committing.  Five minutes of patience costs nothing; losing a design
#: swap (or a solved batch of frames) costs a re-run.
BUSY_TIMEOUT_MS = 300_000


def connect(path, read_only: bool = False) -> sqlite3.Connection:
    """Open the manifest with this project's standard lock patience.

    ``timeout`` (a driver-level sleep-and-retry) and ``busy_timeout`` (the
    SQLite-level equivalent) are set to the same value: some code paths
    inside the driver consult one, some the other, and a mismatch is how a
    connection ends up patient in theory and impatient in practice.
    """
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=BUSY_TIMEOUT_MS / 1000)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    return con


def saturation_adu_for(readoutm) -> float:
    """ADU level treated as 'saturated' in the autopsy statistics."""
    r = (readoutm or "").strip().lower()
    return GSENSE_HIGHGAIN_SAT_ADU if r.startswith("high gain") \
        else DEFAULT_SAT_ADU


# ---------------------------------------------------------------------------
# The candidate-universe base query lives in macro_core.astrom (BASE_SQL /
# fetch_candidates) — ONE definition shared with the report renderer, which
# re-runs it to compute per-stratum night coverage.
# ---------------------------------------------------------------------------
fetch_candidates = astrom.fetch_candidates


# ---------------------------------------------------------------------------
# snapshot — freeze the current experiment as the before/after baseline
# ---------------------------------------------------------------------------
#: The tables a snapshot copies, and the prefix it copies them under.  The
#: baseline is what makes the report's delta section possible: without a
#: frozen copy, ``design`` overwrites the previous experiment and the only
#: record of the old rates is whatever a human happened to write down.
SNAPSHOT_TABLES = ("s1_strata", "s1_populations", "s1_solve_experiment",
                   "s1_failure_autopsy", "s1_build_meta")
BASELINE_PREFIX = "s1_baseline_"


def table_exists(con, name: str) -> bool:
    """True when ``name`` is a table in this database."""
    return con.execute(
        "SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?",
        (name,)).fetchone()[0] > 0


def cmd_snapshot(args) -> int:
    """Copy the live s1_* tables to s1_baseline_* — the comparison baseline.

    Deliberately a SEPARATE command rather than something ``design`` does
    automatically.  If ``design`` snapshotted on every run, the second run
    would overwrite the baseline with the very design it is supposed to be
    compared against, and the delta would silently become zero.  Taking
    the snapshot is therefore an explicit act, and it refuses to clobber an
    existing baseline without ``--force``.
    """
    con = connect(args.manifest)
    with closing(con):
        cur = con.cursor()
        existing = [t for t in SNAPSHOT_TABLES
                    if table_exists(con, BASELINE_PREFIX + t[len("s1_"):])]
        if existing and not args.force:
            print("snapshot: a baseline already exists "
                  f"({len(existing)} tables) — refusing to overwrite it.  "
                  "Re-run with --force only if you really mean to discard "
                  "the comparison baseline.")
            return 1
        missing = [t for t in SNAPSHOT_TABLES if not table_exists(con, t)]
        if missing:
            print(f"snapshot: nothing to freeze — missing {missing}")
            return 1
        cur.execute("BEGIN IMMEDIATE")          # see cmd_design's note
        for t in SNAPSHOT_TABLES:
            dest = BASELINE_PREFIX + t[len("s1_"):]
            cur.execute(f"DROP TABLE IF EXISTS {dest}")
            # CREATE ... AS SELECT copies data and column names, not
            # constraints — which is exactly right for a frozen record:
            # the baseline is evidence, never something to write into.
            cur.execute(f"CREATE TABLE {dest} AS SELECT * FROM {t}")
        cur.execute(f"DROP TABLE IF EXISTS {BASELINE_PREFIX}meta")
        cur.execute(f"""CREATE TABLE {BASELINE_PREFIX}meta
                        (key TEXT PRIMARY KEY, value TEXT)""")
        cur.executemany(
            f"INSERT INTO {BASELINE_PREFIX}meta VALUES (?,?)",
            [("frozen_utc", datetime.now(timezone.utc).isoformat()),
             ("label", args.label),
             ("n_sampled", str(con.execute(
                 "SELECT count(*) FROM s1_solve_experiment"
             ).fetchone()[0])),
             ("n_solved", str(con.execute(
                 "SELECT count(*) FROM s1_solve_experiment "
                 "WHERE status='solved'").fetchone()[0])),
             ("code_version", str(con.execute(
                 "SELECT value FROM s1_build_meta WHERE key='code_version'"
             ).fetchone()[0]))])
        con.commit()
        print(f"snapshot: froze {len(SNAPSHOT_TABLES)} tables as "
              f"{BASELINE_PREFIX}* — label {args.label!r}")
    return 0


def carry_over_map(con) -> dict:
    """Prior per-frame solve results, keyed by obs_rowid.

    Read from the BASELINE when one exists (the frozen record of the
    previous experiment), else from the live table.  Only finished rows
    are carried — a 'pending' row carries no information.
    """
    src = (BASELINE_PREFIX + "solve_experiment"
           if table_exists(con, BASELINE_PREFIX + "solve_experiment")
           else "s1_solve_experiment")
    if not table_exists(con, src):
        return {}
    cols = ["obs_rowid", "status", "used_hint", "solve_time_s", "solved_ra",
            "solved_dec", "pixscale_arcsec", "rotation_deg", "n_matched",
            "rms_arcsec", "log_tail"]
    return {r[0]: dict(zip(cols, r)) for r in con.execute(
        f"SELECT {', '.join(cols)} FROM {src} WHERE status != 'pending'")}


# ---------------------------------------------------------------------------
# design — build strata, sample, write the s1 design tables
# ---------------------------------------------------------------------------
def cmd_design(args) -> int:
    con = connect(args.manifest)
    with closing(con):
        rows = fetch_candidates(con)
        # -- Census of the exclusions ------------------------------------
        # The report must show what astrometry CANNOT apply to, with
        # counts, before showing what it can.  Every count below is
        # computed by calling the SAME pure gate the design uses, in the
        # SAME order the gate applies them, so a class here can never
        # describe a frame the gate treated differently.
        n_total = len(rows)
        n_spectrum = sum(astrom.is_measured_spectrum(r["dispersion"])
                         for r in rows)
        n_vocab = sum((not astrom.is_measured_spectrum(r["dispersion"]))
                      and astrom.is_calib_vocab_filter(r["filter"])
                      for r in rows)
        n_window = sum((not astrom.is_measured_spectrum(r["dispersion"]))
                       and (not astrom.is_calib_vocab_filter(r["filter"]))
                       and astrom.is_window_geometry(r["naxis1"],
                                                     r["naxis2"])
                       for r in rows)
        candidates = [r for r in rows if astrom.is_solvable_candidate(r)]
        # -- The label-vs-measurement cross-tab --------------------------
        # One row per (label class x dispersion class) cell, carrying both
        # gates' verdicts and the movement between them.  This is the
        # evidence for the correction: the report reads it instead of
        # asserting numbers, and a cell that moves frames cannot hide.
        cells: dict[tuple, int] = {}
        for r in rows:
            key = ("grism_label" if astrom.is_grism_filter(r["filter"])
                   else "plain_label",
                   astrom.dispersion_class(r["dispersion"]),
                   astrom.is_solvable_candidate_by_label(r),
                   astrom.is_solvable_candidate(r),
                   astrom.gate_movement(r))
            cells[key] = cells.get(key, 0) + 1
        gate_rows = [(lab, disp, "included" if old else "excluded",
                      "included" if new else "excluded", move, n)
                     for (lab, disp, old, new, move), n
                     in sorted(cells.items())]
        # Headline totals the report interpolates directly.
        n_label_solvable = sum(astrom.is_solvable_candidate_by_label(r)
                               for r in rows)
        n_moved_in = sum(1 for r in rows
                         if astrom.gate_movement(r) == "moved_in")
        n_moved_out = sum(1 for r in rows
                          if astrom.gate_movement(r) == "moved_out")
        n_in_direct = sum(1 for r in rows
                          if astrom.gate_movement(r) == "moved_in"
                          and astrom.dispersion_class(r["dispersion"])
                          == astrom.DIRECT_VERDICT)
        n_in_indet = sum(1 for r in rows
                         if astrom.gate_movement(r) == "moved_in"
                         and astrom.dispersion_class(r["dispersion"])
                         == astrom.INDETERMINATE_VERDICT)
        n_unmeasured = sum(
            1 for r in rows
            if astrom.dispersion_class(r["dispersion"])
            == astrom.UNMEASURED_CLASS)
        # Classify every candidate; count strata populations + the residue.
        by_stratum: dict[str, list[dict]] = {s.stratum_id: []
                                             for s in astrom.STRATA}
        n_residue = 0
        for r in candidates:
            sid = astrom.classify_stratum(r)
            if sid is None:
                n_residue += 1
            else:
                by_stratum[sid].append(r)
        # The SAME classification under the retired label gate — needed for
        # the per-stratum before/after population column.
        label_pop: dict[str, int] = {s.stratum_id: 0 for s in astrom.STRATA}
        for r in rows:
            sid = astrom.classify_stratum_by_label(r)
            if sid is not None:
                label_pop[sid] += 1
        # Prior results, kept for frames that survive into the new sample.
        prior = carry_over_map(con) if args.carry_over else {}
        n_carried = 0
        # Draw the reproducible sample per stratum.
        sample_rows = []
        strata_rows = []
        for s in astrom.STRATA:
            pop = by_stratum[s.stratum_id]
            ids = {r["obs_rowid"]: r for r in pop}
            picked = astrom.sample_frames(list(ids), args.n_per_stratum,
                                          astrom.SAMPLE_SEED, s.stratum_id)
            strata_rows.append((s.stratum_id, s.population, s.description,
                                len(pop), len(picked), astrom.SAMPLE_SEED,
                                label_pop[s.stratum_id]))
            for order, rowid in enumerate(picked):
                r = ids[rowid]
                # Carry-over: a solve outcome is a property of the FRAME
                # and the solver configuration, neither of which this
                # change touched.  Re-solving a frame that was already
                # solved under the identical configuration would spend
                # CPU to reproduce a known answer AND would let solver
                # nondeterminism leak into a delta that is supposed to
                # isolate the gate change.  Frames new to the sample are
                # left 'pending' and solved normally.
                p = prior.get(rowid) if args.carry_over else None
                if p is not None:
                    n_carried += 1
                sample_rows.append((
                    r["obs_rowid"], s.stratum_id, r["path"],
                    r["target_key"], r["canonical_target"], r["readoutm"],
                    r["xbinning"], r["filter"], r["exptime"],
                    r["naxis1"], r["naxis2"], r["ra_deg"], r["dec_deg"],
                    r["night"], order,
                    p["status"] if p else "pending",
                    p["used_hint"] if p else None,
                    p["solve_time_s"] if p else None,
                    p["solved_ra"] if p else None,
                    p["solved_dec"] if p else None,
                    p["pixscale_arcsec"] if p else None,
                    p["rotation_deg"] if p else None,
                    p["n_matched"] if p else None,
                    p["rms_arcsec"] if p else None,
                    p["log_tail"] if p else None))
        # Tooling inventory recorded as build facts (the report quotes it).
        idx_files = sorted(Path(args.index_dir).glob("index-*.fits"))
        idx_bytes = sum(p.stat().st_size for p in idx_files)
        try:
            sf_version = subprocess.run(
                ["solve-field", "--version"], capture_output=True,
                text=True, timeout=10).stdout.strip()
        except Exception:
            sf_version = "unavailable"
        git = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True,
                             cwd=REPO_ROOT).stdout.strip()
        # Atomic swap: build every s1 table under a _tmp name, then swap
        # inside one transaction — the same idempotence contract as S0b.
        cur = con.cursor()
        # IMMEDIATE, not a plain BEGIN, and the distinction is not cosmetic
        # when another stage is writing this database.  A plain (deferred)
        # BEGIN takes a READ snapshot on its first statement — here the
        # DROP TABLE IF EXISTS sweep, which reads the schema — and only
        # asks for the write lock later, at the first CREATE.  In WAL mode
        # SQLite refuses that upgrade *immediately* with SQLITE_BUSY if any
        # other connection has committed since the snapshot was taken:
        # waiting could deadlock two upgraders, so the busy handler is
        # never consulted and ``busy_timeout`` buys nothing.  That is
        # exactly what happened on the first re-run attempt against the
        # live S1b batch — "database is locked" in under a second, despite
        # five minutes of configured patience.  BEGIN IMMEDIATE asks for
        # the write lock up front, where the busy handler DOES apply, so
        # the swap simply waits its turn behind the batch's commits.
        cur.execute("BEGIN IMMEDIATE")
        for t in ("s1_strata", "s1_populations", "s1_gate_comparison",
                  "s1_solve_experiment", "s1_failure_autopsy",
                  "s1_build_meta"):
            cur.execute(f"DROP TABLE IF EXISTS {t}_s1_tmp")
        cur.execute("""CREATE TABLE s1_strata_s1_tmp (
            stratum_id TEXT PRIMARY KEY, population TEXT, description TEXT,
            n_population INTEGER, n_sample INTEGER, seed INTEGER,
            n_population_label_gate INTEGER)""")
        cur.executemany(
            "INSERT INTO s1_strata_s1_tmp VALUES (?,?,?,?,?,?,?)",
            strata_rows)
        # The key is the FULL cell, not just (label, dispersion): the two
        # OTHER gates (the 'dark' filter-wheel glitch and the sub-frame
        # geometry cut) are shared by both universes and can split a
        # (label x dispersion) cell into an included and an excluded half.
        # Keying on the pair alone collapsed those halves into a UNIQUE
        # violation — which is the schema catching a real distinction, so
        # the schema widened rather than the data being deduplicated away.
        cur.execute("""CREATE TABLE s1_gate_comparison_s1_tmp (
            label_class TEXT, dispersion_class TEXT,
            label_gate TEXT, measured_gate TEXT, movement TEXT,
            n_frames INTEGER,
            PRIMARY KEY (label_class, dispersion_class,
                         label_gate, measured_gate))""")
        cur.executemany(
            "INSERT INTO s1_gate_comparison_s1_tmp VALUES (?,?,?,?,?,?)",
            gate_rows)
        cur.execute("""CREATE TABLE s1_populations_s1_tmp (
            class TEXT PRIMARY KEY, n_frames INTEGER, note TEXT)""")
        cur.executemany(
            "INSERT INTO s1_populations_s1_tmp VALUES (?,?,?)",
            [("unsolved_total", n_total,
              "canonical rawimage Light frames with pltsolvd != 1"),
             ("excluded_measured_spectrum", n_spectrum,
              "S2c MEASURED dispersion traces on this frame — it is a "
              "spectrum, whatever its FILTER string says"),
             ("excluded_calib_vocab_filter", n_vocab,
              "FILTER = 'dark'/'bias'/'flat' header glitch (S0b)"),
             # The note describes the GATE, not a shape anyone has seen.
             # It used to end "(8x3211 strips)", naming the phantom
             # geometry as though it were a real readout mode; the S0e
             # repair showed that signature was a tile-compressed
             # BINTABLE's row length misread as an image width, and this
             # count went to zero.  The gate stays (a genuinely narrow
             # sub-frame would still be unsolvable) — only the claim about
             # what has been seen through it goes.
             ("excluded_window_geometry", n_window,
              f"either axis < {astrom.MIN_SOLVABLE_NAXIS} px — too little "
              "sky for quad matching (see docs/pipeline/s0e_geometry_fix"
              ".html: the 8x3211 signature was a metadata artifact)"),
             ("solvable_candidates", len(candidates),
              "frames the solver can be pointed at"),
             ("candidates_unstratified", n_residue,
              "solvable but in no stratum (small heterogeneous residue)"),
             # --- the retired label gate, for comparison only -----------
             ("label_gate_solvable_candidates", n_label_solvable,
              "what the RETIRED FILTER-label gate would have called the "
              "candidate universe (comparison only — nothing gates on it)"),
             ("gate_moved_in_total", n_moved_in,
              "frames the label gate excluded that the measurement gate "
              "keeps"),
             ("gate_moved_in_direct", n_in_direct,
              "...of those, MEASURED DIRECT IMAGES carrying a grism-looking "
              "FILTER — images the label rule deleted unseen"),
             ("gate_moved_in_indeterminate", n_in_indet,
              "...of those, S2c-indeterminate frames: kept because "
              "exclusion must be earned by a positive measurement"),
             ("gate_moved_out_total", n_moved_out,
              "frames the label gate INCLUDED that S2c measures as "
              "dispersed — the contamination this correction removes"),
             ("included_unmeasured", n_unmeasured,
              "candidates S2c never measured: kept (absence of a "
              "measurement is not evidence of a spectrum); this rule moves "
              "ZERO frames — every unmeasured frame was already in")])
        cur.execute("""CREATE TABLE s1_solve_experiment_s1_tmp (
            obs_rowid INTEGER PRIMARY KEY, stratum_id TEXT, path TEXT,
            target_key TEXT, canonical_target TEXT, readoutm TEXT,
            xbinning REAL, filter TEXT, exptime REAL,
            naxis1 REAL, naxis2 REAL, ra_hint_deg REAL, dec_hint_deg REAL,
            night TEXT, sample_order INTEGER, status TEXT,
            used_hint INTEGER, solve_time_s REAL,
            solved_ra REAL, solved_dec REAL, pixscale_arcsec REAL,
            rotation_deg REAL, n_matched INTEGER, rms_arcsec REAL,
            log_tail TEXT)""")
        cur.executemany("""INSERT INTO s1_solve_experiment_s1_tmp
            (obs_rowid, stratum_id, path, target_key, canonical_target,
             readoutm, xbinning, filter, exptime, naxis1, naxis2,
             ra_hint_deg, dec_hint_deg, night, sample_order, status,
             used_hint, solve_time_s, solved_ra, solved_dec,
             pixscale_arcsec, rotation_deg, n_matched, rms_arcsec,
             log_tail)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            sample_rows)
        # ``dispersion_class`` is recorded on every autopsy row so the
        # report can PROVE, from the table itself, that no frame in the
        # failure taxonomy is a measured spectrum — the exact defect this
        # version repairs.  A claim that can only be checked by re-running
        # a join is a claim a reader has to take on trust.
        cur.execute("""CREATE TABLE s1_failure_autopsy_s1_tmp (
            obs_rowid INTEGER PRIMARY KEY, stratum_id TEXT,
            n_sources INTEGER, n_psf_sources INTEGER,
            median_elongation REAL, median_a_px REAL,
            bright_median_a_px REAL, saturated_fraction REAL,
            bkg_rms REAL, saturation_adu REAL,
            diagnosis TEXT, thumb_path TEXT, dispersion_class TEXT)""")
        cur.execute("""CREATE TABLE s1_build_meta_s1_tmp
                       (key TEXT PRIMARY KEY, value TEXT)""")
        cur.executemany("INSERT INTO s1_build_meta_s1_tmp VALUES (?,?)", [
            ("built_utc", datetime.now(timezone.utc).isoformat()),
            ("code_version", S1_CODE_VERSION),
            ("git_commit", git),
            ("solve_field_version", sf_version),
            ("index_dir", str(args.index_dir)),
            ("index_n_files", str(len(idx_files))),
            ("index_bytes", str(idx_bytes)),
            ("index_series", ",".join(sorted({p.name.split("-")[1]
                                              for p in idx_files}))),
            ("sample_seed", str(astrom.SAMPLE_SEED)),
            ("n_per_stratum", str(args.n_per_stratum)),
            ("solve_timeout_s", str(astrom.SOLVE_TIMEOUT_S)),
            ("solve_cpulimit_s", str(astrom.SOLVE_CPU_LIMIT_S)),
            ("hint_radius_deg", str(astrom.HINT_RADIUS_DEG)),
            ("downsample", str(astrom.SOLVE_DOWNSAMPLE)),
            ("workers", str(args.workers)),
            ("candidate_gate", "measured dispersion (frame_dispersion."
                               "verdict = 'dispersed' excludes)"),
            ("carry_over", "on" if args.carry_over else "off"),
            ("n_carried_over", str(n_carried)),
        ])
        for t in ("s1_strata", "s1_populations", "s1_gate_comparison",
                  "s1_solve_experiment", "s1_failure_autopsy",
                  "s1_build_meta"):
            cur.execute(f"DROP TABLE IF EXISTS {t}")
            cur.execute(f"ALTER TABLE {t}_s1_tmp RENAME TO {t}")
        con.commit()
        n_samp = len(sample_rows)
        print(f"design: {n_total:,} unsolved -> {len(candidates):,} "
              f"solvable candidates -> {n_samp} sampled frames across "
              f"{len(astrom.STRATA)} strata (seed {astrom.SAMPLE_SEED})")
        print(f"design: gate = MEASURED dispersion; the retired label gate "
              f"would give {n_label_solvable:,} candidates "
              f"({n_moved_in:,} frames move IN, {n_moved_out:,} move OUT)")
        print(f"design: carried over {n_carried} prior solve results; "
              f"{n_samp - n_carried} frames pending")
        for row in strata_rows:
            delta = row[3] - row[6]
            print(f"  {row[0]:<24} pop {row[3]:>6,}  sample {row[4]:>3}"
                  f"  (label gate {row[6]:>6,}, {delta:+d})")
    return 0


# ---------------------------------------------------------------------------
# run — solve pending frames until the batch budget elapses (resumable)
# ---------------------------------------------------------------------------
def _solve_task(row: dict, archive_root: Path, config: str,
                scratch_root: Path) -> tuple[int, dict]:
    """Worker: solve one frame in a private temp dir, always cleaned up."""
    work = Path(tempfile.mkdtemp(prefix=f"s1_{row['obs_rowid']}_",
                                 dir=scratch_root))
    try:
        res = astrom.solve_one_frame(
            str(archive_root / row["path"]), str(work), config,
            row["readoutm"], row["xbinning"],
            row["ra_hint_deg"], row["dec_hint_deg"])
        return row["obs_rowid"], res
    finally:
        shutil.rmtree(work, ignore_errors=True)


def cmd_run(args) -> int:
    t_start = time.monotonic()
    scratch_root = Path(tempfile.mkdtemp(prefix="s1_batch_"))
    con = connect(args.manifest)
    n_done = 0
    try:
        while time.monotonic() - t_start < args.batch_seconds:
            # Fetch the next chunk of pending frames (small chunks so a
            # killed batch loses at most one chunk of already-solved work).
            cols = ["obs_rowid", "path", "readoutm", "xbinning",
                    "ra_hint_deg", "dec_hint_deg"]
            rows = [dict(zip(cols, r)) for r in con.execute(
                f"""SELECT {', '.join(cols)} FROM s1_solve_experiment
                    WHERE status = 'pending'
                    ORDER BY stratum_id, sample_order LIMIT ?""",
                (args.workers * 3,))]
            if not rows:
                print("run: nothing pending — experiment complete")
                break
            with concurrent.futures.ThreadPoolExecutor(args.workers) as ex:
                futures = [ex.submit(_solve_task, r, args.archive,
                                     args.config, scratch_root)
                           for r in rows]
                for fut in concurrent.futures.as_completed(futures):
                    rowid, res = fut.result()
                    # Single-writer discipline: only this thread touches
                    # the DB; workers return plain dicts.
                    con.execute("""UPDATE s1_solve_experiment SET
                        status=?, used_hint=?, solve_time_s=?, solved_ra=?,
                        solved_dec=?, pixscale_arcsec=?, rotation_deg=?,
                        n_matched=?, rms_arcsec=?, log_tail=?
                        WHERE obs_rowid=?""", (
                        res["status"], res["used_hint"],
                        res["solve_time_s"], res["solved_ra"],
                        res["solved_dec"], res["pixscale_arcsec"],
                        res["rotation_deg"], res["n_matched"],
                        res["rms_arcsec"], res["log_tail"], rowid))
                    n_done += 1
                con.commit()
            elapsed = time.monotonic() - t_start
            print(f"run: {n_done} frames this batch, "
                  f"{elapsed:.0f}s elapsed", flush=True)
        n_left = con.execute("SELECT count(*) FROM s1_solve_experiment "
                             "WHERE status='pending'").fetchone()[0]
        print(f"run: batch done ({n_done} frames); {n_left} still pending")
    finally:
        con.close()
        shutil.rmtree(scratch_root, ignore_errors=True)
    return 0


# ---------------------------------------------------------------------------
# autopsy — image statistics + thumbnails for a sample of failures
# ---------------------------------------------------------------------------
def cmd_autopsy(args) -> int:
    import numpy as np
    from astropy.io import fits
    THUMB_DIR.mkdir(parents=True, exist_ok=True)
    con = connect(args.manifest)
    with closing(con):
        # EVERY failure is autopsied — no per-stratum cap.  (v1.0 capped
        # at 8 per stratum, which left 25 of 58 failures undiagnosed in
        # exactly the NO-GO strata where the taxonomy matters most.)
        # ``bad_solve`` rows — false-positive WCS caught by the
        # acceptance gate — are failures too and get the same treatment.
        # The S2c verdict travels with each failure so the taxonomy can be
        # audited for spectra without a second query.
        has_disp = astrom.has_dispersion_table(con)
        disp_sel = ("d.verdict" if has_disp else "NULL")
        disp_join = ("LEFT JOIN frame_dispersion d USING (obs_rowid)"
                     if has_disp else "")
        rows = con.execute(f"""
            SELECT s.obs_rowid, s.stratum_id, s.path, s.readoutm, {disp_sel}
            FROM s1_solve_experiment s {disp_join}
            WHERE s.status IN ('unsolved', 'timeout', 'bad_solve')
            ORDER BY s.stratum_id, s.sample_order""").fetchall()
        con.execute("DELETE FROM s1_failure_autopsy")
        for rowid, sid, relpath, readoutm, verdict in rows:
            disp_class = astrom.dispersion_class(verdict)
            src = args.archive / relpath
            try:
                with fits.open(src) as hdul:
                    # fpack files put pixels in HDU 1; plain FITS in HDU 0.
                    hdu = hdul[1] if len(hdul) > 1 and \
                        hdul[0].data is None else hdul[0]
                    data = np.asarray(hdu.data, dtype=np.float32)
                sat = saturation_adu_for(readoutm)
                m = astrom.image_metrics(data, sat)
                # Thumbnail: asinh stretch between the 5th and 99.5th
                # percentiles — faint stars visible, saturation obvious.
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt
                lo, hi = np.percentile(data, [5.0, 99.5])
                stretched = np.arcsinh((data - lo)
                                       / max(hi - lo, 1e-3) * 10.0)
                fig, ax = plt.subplots(figsize=(3.2, 3.2), dpi=80)
                ax.imshow(stretched, cmap="gray", origin="lower")
                ax.axis("off")
                thumb = THUMB_DIR / f"{rowid}.png"
                fig.savefig(thumb, bbox_inches="tight", pad_inches=0.02)
                plt.close(fig)
                con.execute("""INSERT OR REPLACE INTO s1_failure_autopsy
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    rowid, sid, m["n_sources"], m["n_psf_sources"],
                    m["median_elongation"], m["median_a_px"],
                    m["bright_median_a_px"], m["saturated_fraction"],
                    m["bkg_rms"], sat,
                    astrom.diagnose_failure(m["n_sources"],
                                            m["n_psf_sources"],
                                            m["median_elongation"],
                                            m["saturated_fraction"],
                                            m["bright_median_a_px"]),
                    str(thumb), disp_class))
            except Exception as exc:    # noqa: BLE001 — recorded, not hidden
                con.execute("""INSERT OR REPLACE INTO s1_failure_autopsy
                    (obs_rowid, stratum_id, diagnosis, thumb_path,
                     dispersion_class)
                    VALUES (?,?,?,NULL,?)""",
                            (rowid, sid,
                             f"unreadable: {type(exc).__name__}",
                             disp_class))
            con.commit()
        n = con.execute(
            "SELECT count(*) FROM s1_failure_autopsy").fetchone()[0]
        # The invariant this stage now guarantees, checked out loud: a
        # measured spectrum can no longer reach the failure taxonomy,
        # because the candidate gate excluded it from the universe.
        n_spec = con.execute(
            "SELECT count(*) FROM s1_failure_autopsy "
            "WHERE dispersion_class = ?",
            (astrom.DISPERSED_VERDICT,)).fetchone()[0]
        print(f"autopsy: {n} failures examined (all of them — no cap); "
              f"{n_spec} of them are measured spectra "
              f"(must be 0 under the measured-dispersion gate)")
    return 0


# ---------------------------------------------------------------------------
# status / report
# ---------------------------------------------------------------------------
def cmd_status(args) -> int:
    con = connect(args.manifest, read_only=True)
    with closing(con):
        for sid, tot, pend, ok, bad, badsol, tmo, err, med in con.execute("""
            SELECT stratum_id, count(*),
                   sum(status = 'pending'), sum(status = 'solved'),
                   sum(status = 'unsolved'), sum(status = 'bad_solve'),
                   sum(status = 'timeout'), sum(status = 'error'),
                   round(avg(solve_time_s), 1)
            FROM s1_solve_experiment GROUP BY stratum_id
            ORDER BY stratum_id"""):
            print(f"{sid:<24} n={tot:>3} pending={pend:>3} solved={ok:>3} "
                  f"unsolved={bad:>3} bad_solve={badsol:>3} "
                  f"timeout={tmo:>3} error={err:>3} avg_t={med}s")
    return 0


def cmd_report(args) -> int:
    from macro_core.report_s1 import render_report
    path = render_report(Path(args.manifest))
    print(f"report: {path}")
    return 0


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="S1 stratified astrometry go/no-go experiment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    p.add_argument("--config", type=str, default=str(DEFAULT_CONFIG))
    p.add_argument("--index-dir", type=Path, default=DEFAULT_INDEX_DIR)
    p.add_argument("--workers", type=int, default=astrom.DEFAULT_WORKERS)
    sub = p.add_subparsers(dest="command", required=True)
    sn = sub.add_parser("snapshot",
                        help="freeze the current experiment as the "
                             "before/after baseline")
    sn.add_argument("--label", type=str, default="previous design",
                    help="what this baseline IS, in a few words")
    sn.add_argument("--force", action="store_true",
                    help="overwrite an existing baseline (discards it)")
    d = sub.add_parser("design", help="build strata + draw samples")
    d.add_argument("--n-per-stratum", type=int,
                   default=astrom.N_PER_STRATUM)
    # Carry-over defaults ON: a re-design that only changes WHICH frames
    # are sampled must not re-spend CPU reproducing solve outcomes it
    # already has for the frames it keeps.
    d.add_argument("--no-carry-over", dest="carry_over",
                   action="store_false",
                   help="re-solve every sampled frame from scratch, even "
                        "ones with a recorded prior outcome")
    d.set_defaults(carry_over=True)
    r = sub.add_parser("run", help="solve pending frames (resumable batch)")
    r.add_argument("--batch-seconds", type=int, default=420)
    sub.add_parser("autopsy", help="post-mortem failures")
    sub.add_parser("status", help="progress table")
    sub.add_parser("report", help="render the evidence report")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    return {"snapshot": cmd_snapshot, "design": cmd_design, "run": cmd_run,
            "autopsy": cmd_autopsy, "status": cmd_status,
            "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
