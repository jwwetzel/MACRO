#!/usr/bin/env python
"""Measure filter identity FROM THE PIXELS, frame by frame (stage S2c).

WHAT THIS SCRIPT DOES
---------------------
The FILTER header card cannot be trusted to say whether a frame is a
spectrum.  The wheel's grism slots were named three different ways across
three years, and the earliest name — the bare slot number ``6`` — turns out
to be dispersed on some targets and direct on others.  Rather than argue
about the label, this campaign opens every candidate frame, extracts its
sources, and measures whether the light was dispersed.

The physics and the decision rules live in ``rlmt_diagnostics.dispersion``
(pure, unit-tested).  This script is only the plumbing: build a queue, run
it resumably across a worker pool, and write the numbers into ONE new
manifest table, ``frame_dispersion``.  No existing table is modified, and
the archive is opened strictly read-only.

THE QUEUE
---------
Two populations, both needed:

* ``candidate`` — every frame whose FILTER is a grism name or is disputed:
  ``6``, ``W``, ``hrg``, ``lrg``, ``HaGrism``, ``OGGrism``, ``HaG``, and the
  stragglers ``w`` / ``lrgblue``.  These are the frames the projects need
  verdicts on.
* ``control``   — a random sample of frames whose FILTER is an ordinary
  photometric band nobody disputes (``g r i V R I B L``).  These are the
  ground truth for the DIRECT side of the calibration, and the only way to
  measure the classifier's false-positive rate honestly.  Without them a
  "100% of grism frames read as dispersed" claim would be unfalsifiable.

Note that the labelled grism frames (``hrg`` etc.) serve double duty: they
are candidates AND the ground truth for the DISPERSED side.

MEASURE FIRST, JUDGE LATER
--------------------------
Every row stores the raw measured numbers alongside the verdict.  The
measurement costs a frame decompression and a source extraction; the verdict
costs nothing.  Keeping them separate means the thresholds can be
recalibrated over the whole archive — and the whole archive reclassified —
without touching a single pixel again.  That is what ``reclassify`` does.

SUBCOMMANDS
-----------
    build       construct the queue (refuses to clobber progress; --rebuild)
    run         measure pending frames; SAFE TO RE-RUN — a killed run loses
                only the frames in flight, which stay pending
    status      progress + per-label verdict tallies (read-only)
    reclassify  recompute verdicts from the STORED numbers, no pixel reads
                (v1.2: trace verdict, then lozenge templates, morphology
                class and the two-witness verdict)
    calibrate   print the known-label separation table (read-only)

    -- S2c v1.2, finding F-6 (the second witness) --
    truth-load       load the tracked truth set into s2c_truth
    queue-sn         queue every never-measured SN 2023ixf frame (SN-G0d)
    morph-queue      mark the rows that get a background map (scope in
                     the function's docstring)
    morph-run        measure pending background maps (resumable)
    morph-calibrate  template_r of the undisputed-label populations
    truth-score      both confusion matrices with Wilson intervals

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/run_s2c_dispersion.py build
    $PY pipeline/scripts/run_s2c_dispersion.py run --workers 6
    $PY pipeline/scripts/run_s2c_dispersion.py status
    $PY pipeline/scripts/run_s2c_dispersion.py calibrate

``run --max-seconds N`` returns cleanly after N seconds so the campaign can
be driven from an environment with a command timeout; just call it again.

CONCURRENCY NOTE
----------------
The manifest is a WAL database that another stage (the S1 astrometry batch)
may be writing at the same time.  Readers therefore never block, and this
script keeps its own write transactions short — results are flushed in small
batches — so the two writers interleave instead of colliding.  Every
connection sets a five-minute busy timeout as a backstop.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import random
import sqlite3
import sys
import time
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

# Make the pipeline package importable regardless of the working directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from rlmt_diagnostics import dispersion as dsp                  # noqa: E402
from rlmt_diagnostics.dispersion import DISPERSION_CODE_VERSION  # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: Five minutes.  The other writer's transactions are short; this is a
#: backstop against a long checkpoint, not an expected wait.
BUSY_TIMEOUT_MS = 300_000

#: FILTER values whose identity is in question or is a known grism name.
CANDIDATE_FILTERS = ("6", "W", "w", "hrg", "lrg", "HaGrism", "OGGrism",
                     "HaG", "lrgblue")

#: FILTER values nobody disputes are direct imaging — the control ground
#: truth.  Deliberately excludes narrowband (ha/oiii/sii) and the luminance
#: family, whose own identity is not the subject of this study.
CONTROL_FILTERS = ("g", "r", "i", "V", "R", "I", "B", "L")

#: How many control frames to draw per direct label.  400 x 8 = 3,200 is
#: enough to put a useful bound on the false-positive rate while adding only
#: ~13% to the campaign's runtime.
CONTROL_PER_FILTER = 400

#: Seed for the control draw, so the campaign is reproducible.
CONTROL_SEED = 20260818

#: A SECOND, disjoint draw from the same undisputed labels — the holdout.
#:
#: The control sample above cannot honestly be used to quote an error rate,
#: because the thresholds were moved in response to frames inside it: the
#: PA-scatter gate went from 20 deg to 5 after a control ``r`` frame, the
#: sparsity gate was invented after a control ``L`` frame of M57, and the
#: calibration-frame exclusion was added after eleven control ``B`` frames
#: came back "dispersed".  Every one of those frames is still scored in the
#: control totals.  A number fitted on the same data it is quoted over is a
#: lower bound on the error, not an estimate of it, and the first version of
#: this campaign published it as though it were the latter.
#:
#: The holdout closes that hole the only way it can be closed: a fresh draw
#: under a different seed, EXPLICITLY EXCLUDING every obs_rowid already in
#: the table, measured once with the thresholds frozen, and never consulted
#: while tuning.  If the holdout rate matches the control rate, the fitting
#: cost nothing measurable; if it is worse, the control number was optimistic
#: and the report must say by how much.  Either answer is worth having.
HOLDOUT_SEED = 20260819
HOLDOUT_PER_FILTER = 300

#: Frames that are not observations of the sky, and must not be measured.
#:
#: The first smoke run learned this the hard way: eleven of eighteen
#: known-direct ``B`` frames came back "dispersed", and every one of them
#: was a twilight flat.  A flat field has no stars — the only things the
#: extractor finds are dust shadows and detector column defects, which are
#: perfectly straight, perfectly parallel (they ARE columns), and therefore
#: an ideal forgery of a grism's shared dispersion axis.  Calibration frames
#: cannot answer the question this campaign asks, so they never enter it.
NON_SCIENCE_IMAGETYP = ("Flat Field", "FLAT", "Dark Frame", "Bias Frame",
                        "DARK", "BIAS")

#: The calibration subtree, excluded for the same reason (its master frames
#: often carry no IMAGETYP card at all, so the card test alone misses them).
NON_SCIENCE_TREES = ("calib",)

#: How many frames a worker chunk fetches, and how often results are
#: flushed.  Small enough that a kill loses little and that each write
#: transaction is brief.
CHUNK = 120

#: Queue priority — frames are measured in this order, lowest number first.
#:
#: The campaign takes hours against a shared spinning disk, so the order is
#: chosen so that a run interrupted at ANY point has already answered the
#: questions that were asked.  The disputed labels come first because they
#: are the entire question; the control follows immediately because without
#: a measured false-positive rate no verdict on the disputed labels can be
#: believed; the small transitional vocabulary comes next because it bridges
#: the naming epochs; and the 16k-frame hrg/lrg bulk — whose labels are not
#: in dispute and which only refines the census — comes last.
PRIORITY = {
    "6": 1, "W": 1, "w": 1,                       # the disputed slots
    # control frames get priority 2, assigned in build()
    "HaGrism": 3, "OGGrism": 3, "HaG": 3, "lrgblue": 3,
    "hrg": 4, "lrg": 4,
}
CONTROL_PRIORITY = 2
#: The holdout is measured last: it must not exist as a temptation while the
#: thresholds are still moving.
HOLDOUT_PRIORITY = 5
DEFAULT_PRIORITY = 4


def utcnow() -> str:
    """ISO-8601 UTC timestamp for log lines and DB facts."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def connect(path, read_only: bool = False) -> sqlite3.Connection:
    """Open the manifest with the shared concurrency settings."""
    if read_only:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
    else:
        con = sqlite3.connect(str(path), timeout=300)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    return con


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
#: The stored schema.  ``status`` drives resumability; the measurement
#: columns are what makes reclassification free.
SCHEMA = """
CREATE TABLE IF NOT EXISTS frame_dispersion (
    obs_rowid        INTEGER PRIMARY KEY,
    path             TEXT NOT NULL,
    tree             TEXT,
    filter           TEXT,
    night            TEXT,
    canonical_target TEXT,
    exptime          REAL,
    xbinning         REAL,
    era_id           INTEGER,
    population       TEXT NOT NULL,   -- candidate | control
    priority         INTEGER NOT NULL DEFAULT 4,   -- lower is measured first
    status           TEXT NOT NULL,   -- pending | measured | unreadable
    -- measured numbers (NULL until measured)
    n_detected       INTEGER,
    n_sources        INTEGER,
    n_bright         INTEGER,
    median_ab        REAL,
    max_ab           REAL,
    pa_median        REAL,
    pa_scatter       REAL,
    n_trace          INTEGER,
    trace_frac       REAL,
    trace_ab         REAL,
    trace_a_px       REAL,
    trace_pa         REAL,
    trace_pa_scatter REAL,
    detect_sigma     REAL,
    height           INTEGER,
    width            INTEGER,
    -- judgement (recomputable from the columns above)
    verdict          TEXT,
    strength_class   TEXT,
    reason           TEXT,
    -- bookkeeping
    measure_s        REAL,
    error            TEXT,
    code_version     TEXT,
    measured_at      TEXT
);
CREATE INDEX IF NOT EXISTS ix_fdisp_status ON frame_dispersion(status);
CREATE INDEX IF NOT EXISTS ix_fdisp_filter ON frame_dispersion(filter);
CREATE INDEX IF NOT EXISTS ix_fdisp_verdict ON frame_dispersion(verdict);
CREATE TABLE IF NOT EXISTS s2c_build_meta (
    key TEXT PRIMARY KEY, value TEXT
);
"""

#: Frame columns copied into the queue so later analysis never needs a join.
_FRAME_COLS = ("obs_rowid", "path", "tree", "filter", "night",
               "canonical_target", "exptime", "xbinning", "era_id")


def cmd_build(args) -> int:
    con = connect(args.manifest)
    with closing(con):
        con.executescript(SCHEMA)
        done = con.execute(
            "SELECT count(*) FROM frame_dispersion "
            "WHERE status != 'pending'").fetchone()[0]
        if done and not args.rebuild:
            print(f"build: frame_dispersion already holds {done:,} measured "
                  "frames — refusing to clobber progress. "
                  "Pass --rebuild to requeue everything.")
            return 1
        if args.rebuild:
            con.execute("DELETE FROM frame_dispersion")

        cols = ", ".join(_FRAME_COLS)
        # Science-only clause, shared by both populations: a real canonical
        # frame, readable, pointed at the sky rather than at a flat screen.
        it_marks = ",".join("?" * len(NON_SCIENCE_IMAGETYP))
        tr_marks = ",".join("?" * len(NON_SCIENCE_TREES))
        science = (f"is_canonical = 1 AND error IS NULL "
                   f"AND (imagetyp IS NULL OR imagetyp NOT IN ({it_marks})) "
                   f"AND (tree IS NULL OR tree NOT IN ({tr_marks}))")
        science_params = (*NON_SCIENCE_IMAGETYP, *NON_SCIENCE_TREES)

        # -- candidates: every disputed or grism-named frame ---------------
        marks = ",".join("?" * len(CANDIDATE_FILTERS))
        cand = con.execute(
            f"""SELECT {cols} FROM frames
                WHERE {science} AND filter IN ({marks})""",
            (*science_params, *CANDIDATE_FILTERS)).fetchall()

        # -- controls: a reproducible random draw per undisputed label -----
        rng = random.Random(CONTROL_SEED)
        ctrl = []
        for filt in CONTROL_FILTERS:
            pool = con.execute(
                f"""SELECT {cols} FROM frames
                    WHERE {science} AND filter = ?
                    ORDER BY obs_rowid""",
                (*science_params, filt)).fetchall()
            # Sample without replacement from the FULL pool, so the control
            # spans every night and era the label ever appeared in rather
            # than clustering on whichever rows happen to sort first.
            take = min(CONTROL_PER_FILTER, len(pool))
            ctrl.extend(rng.sample(pool, take))

        # ``filter`` is the 4th column of _FRAME_COLS; priority keys off it
        # for candidates and is fixed for the control population.
        fi = _FRAME_COLS.index("filter")
        rows = ([(*r, "candidate", PRIORITY.get(r[fi], DEFAULT_PRIORITY),
                  "pending") for r in cand]
                + [(*r, "control", CONTROL_PRIORITY, "pending")
                   for r in ctrl])
        con.executemany(
            f"INSERT OR REPLACE INTO frame_dispersion "
            f"({cols}, population, priority, status) "
            f"VALUES ({','.join('?' * (len(_FRAME_COLS) + 3))})", rows)
        for k, v in (("built_at", utcnow()),
                     ("code_version", DISPERSION_CODE_VERSION),
                     ("archive_root", str(args.archive)),
                     ("n_candidate", str(len(cand))),
                     ("n_control", str(len(ctrl))),
                     ("control_per_filter", str(CONTROL_PER_FILTER)),
                     ("control_seed", str(CONTROL_SEED))):
            con.execute("INSERT OR REPLACE INTO s2c_build_meta VALUES (?,?)",
                        (k, v))
        con.commit()
        print(f"build: queued {len(cand):,} candidate + {len(ctrl):,} control "
              f"= {len(rows):,} frames")
    return 0


def cmd_holdout(args) -> int:
    """Queue the out-of-sample holdout: a fresh, disjoint control draw.

    Idempotent by construction — it inserts only obs_rowids that are not
    already in ``frame_dispersion`` at all, so running it twice adds nothing
    and can never disturb a measured row.
    """
    con = connect(args.manifest)
    with closing(con):
        con.executescript(SCHEMA)
        have = con.execute(
            "SELECT count(*) FROM frame_dispersion "
            "WHERE population = 'holdout'").fetchone()[0]
        if have and not args.rebuild:
            print(f"holdout: already queued ({have:,} frames) — nothing to "
                  "do.  Pass --rebuild to draw a fresh one.")
            return 0
        if args.rebuild:
            con.execute("DELETE FROM frame_dispersion "
                        "WHERE population = 'holdout'")

        cols = ", ".join(_FRAME_COLS)
        it_marks = ",".join("?" * len(NON_SCIENCE_IMAGETYP))
        tr_marks = ",".join("?" * len(NON_SCIENCE_TREES))
        # Identical science clause to cmd_build — the holdout must be drawn
        # from exactly the population the control was, or it is not a
        # holdout, it is a different experiment.
        science = (f"is_canonical = 1 AND error IS NULL "
                   f"AND (imagetyp IS NULL OR imagetyp NOT IN ({it_marks})) "
                   f"AND (tree IS NULL OR tree NOT IN ({tr_marks}))")
        science_params = (*NON_SCIENCE_IMAGETYP, *NON_SCIENCE_TREES)

        rng = random.Random(HOLDOUT_SEED)
        rows, per_filter = [], []
        for filt in CONTROL_FILTERS:
            pool = con.execute(
                f"""SELECT {cols} FROM frames
                    WHERE {science} AND filter = ?
                      AND obs_rowid NOT IN (SELECT obs_rowid
                                            FROM frame_dispersion)
                    ORDER BY obs_rowid""",
                (*science_params, filt)).fetchall()
            take = min(HOLDOUT_PER_FILTER, len(pool))
            drawn = rng.sample(pool, take)
            per_filter.append((filt, len(pool), take))
            rows.extend((*r, "holdout", HOLDOUT_PRIORITY, "pending")
                        for r in drawn)

        con.executemany(
            f"INSERT OR IGNORE INTO frame_dispersion "
            f"({cols}, population, priority, status) "
            f"VALUES ({','.join('?' * (len(_FRAME_COLS) + 3))})", rows)
        for k, v in (("holdout_built_at", utcnow()),
                     ("holdout_seed", str(HOLDOUT_SEED)),
                     ("holdout_per_filter", str(HOLDOUT_PER_FILTER)),
                     ("n_holdout", str(len(rows)))):
            con.execute("INSERT OR REPLACE INTO s2c_build_meta VALUES (?,?)",
                        (k, v))
        con.commit()
        for filt, avail, take in per_filter:
            note = "  (pool exhausted)" if take < HOLDOUT_PER_FILTER else ""
            print(f"  {filt:<3} {take:>5} drawn of {avail:>7} unused{note}")
        print(f"holdout: queued {len(rows):,} frames, disjoint from the "
              f"control sample, seed {HOLDOUT_SEED}")
    return 0


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
def _frame_task(task: dict) -> tuple[int, dict]:
    """Measure ONE frame.  Runs in a worker process; must never raise.

    A frame that cannot be read is a fact about the archive, not a crash:
    it is recorded as ``unreadable`` with the exception text, and the batch
    moves on.
    """
    t0 = time.time()
    out = {"status": "measured", "error": None}
    try:
        shape, bmap = dsp.measure_file_with_map(task["abs_path"])
        verdict = dsp.classify_frame(shape)
        out.update(shape.as_dict())
        out.update(_morph_fields(bmap))
        # The two-witness verdict is set by reclassify (templates need
        # the whole population); until then the trace verdict stands.
        out["verdict"] = verdict.verdict
        out["verdict_traces"] = verdict.verdict
        out["strength_class"] = verdict.strength_class
        out["reason"] = verdict.reason
    except Exception as exc:                      # noqa: BLE001 — recorded
        out["status"] = "unreadable"
        out["error"] = f"{type(exc).__name__}: {exc}"[:300]
    out["measure_s"] = round(time.time() - t0, 3)
    out["code_version"] = DISPERSION_CODE_VERSION
    out["measured_at"] = utcnow()
    return task["obs_rowid"], out


#: Every column a worker may fill, in UPDATE order.
_RESULT_COLS = ["status", "n_detected", "n_sources", "n_bright",
                "median_ab", "max_ab", "pa_median", "pa_scatter",
                "n_trace", "trace_frac", "trace_ab", "trace_a_px",
                "trace_pa", "trace_pa_scatter", "detect_sigma",
                "height", "width", "verdict", "strength_class", "reason",
                "measure_s", "error", "code_version", "measured_at",
                "verdict_traces", "morph_status", "bg_corner", "bg_edge_tb",
                "bg_peak", "bg_sigma", "bg_contrast", "bg_edge_ratio",
                "bg_area_frac", "bg_map"]


def _flush(con, results: list[tuple[int, dict]]) -> None:
    """Write a batch of results in ONE short transaction.

    Short is the point: another stage may be writing this WAL database, and
    a long-held write lock is how two cooperating jobs turn into a deadlock.
    """
    sets = ", ".join(f"{c} = ?" for c in _RESULT_COLS)
    con.executemany(
        f"UPDATE frame_dispersion SET {sets} WHERE obs_rowid = ?",
        [([r.get(c) for c in _RESULT_COLS] + [rid]) for rid, r in results])
    con.commit()


def cmd_run(args) -> int:
    started = time.time()
    con = connect(args.manifest)
    n_done = 0
    with closing(con):
        ensure_morph_schema(con)
        total = con.execute(
            "SELECT count(*) FROM frame_dispersion").fetchone()[0]
        if not total:
            print("run: queue is empty — run 'build' first.")
            return 1
        print(f"run: {DISPERSION_CODE_VERSION}  workers={args.workers}  "
              f"started {utcnow()}", flush=True)
        # ONE pool for the whole invocation.  Building a fresh pool per chunk
        # re-pays the interpreter start plus the numpy/astropy/sep import cost
        # for every worker, every chunk — measured at roughly half the total
        # runtime on the first production attempt.  The workers are stateless,
        # so a long-lived pool is safe.
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.workers) as pool:
            while True:
                if args.max_seconds and (time.time() - started
                                         > args.max_seconds):
                    print(f"run: time budget reached; {n_done:,} measured "
                          "this invocation. Re-run to continue.", flush=True)
                    break
                if args.limit and n_done >= args.limit:
                    break
                rows = con.execute(
                    """SELECT obs_rowid, path FROM frame_dispersion
                       WHERE status = 'pending'
                       ORDER BY priority, obs_rowid LIMIT ?""",
                    (CHUNK,)).fetchall()
                if not rows:
                    print("run: no pending frames left — campaign complete.",
                          flush=True)
                    break
                tasks = [{"obs_rowid": r[0], "path": r[1],
                          "abs_path": str(args.archive / r[1])} for r in rows]
                results = list(pool.map(_frame_task, tasks, chunksize=4))
                _flush(con, results)
                n_done += len(results)
                elapsed = time.time() - started
                rate = n_done / elapsed if elapsed else 0.0
                remaining = con.execute(
                    "SELECT count(*) FROM frame_dispersion "
                    "WHERE status = 'pending'").fetchone()[0]
                eta_min = (remaining / rate / 60.0) if rate else float("nan")
                print(f"  {utcnow()}  measured {n_done:,} this run  "
                      f"({rate:.1f} frame/s)  pending {remaining:,}  "
                      f"ETA {eta_min:.0f} min", flush=True)
    return 0


# ---------------------------------------------------------------------------
# reclassify — recompute verdicts from stored numbers, no pixel reads
# ---------------------------------------------------------------------------
def cmd_reclassify(args) -> int:
    """Re-run the judgement over every measured row.

    This is what makes threshold calibration honest: the thresholds can be
    set AFTER looking at the measured distributions of the known-label
    populations, and applied to the whole archive in seconds, without the
    temptation to tune them by re-measuring a convenient subset.
    """
    con = connect(args.manifest)
    with closing(con):
        cols = ["obs_rowid", "n_detected", "n_sources", "n_bright",
                "median_ab", "max_ab", "pa_median", "pa_scatter", "n_trace",
                "trace_frac", "trace_ab", "trace_a_px", "trace_pa",
                "trace_pa_scatter", "detect_sigma", "height", "width"]
        rows = con.execute(
            f"SELECT {','.join(cols)} FROM frame_dispersion "
            "WHERE status = 'measured'").fetchall()
        updates = []
        for r in rows:
            d = dict(zip(cols, r))
            rid = d.pop("obs_rowid")
            # Integer columns can come back NULL on an empty frame; the
            # dataclass wants real ints there.
            for k in ("n_detected", "n_sources", "n_bright", "n_trace"):
                d[k] = int(d[k] or 0)
            for k in ("height", "width"):
                d[k] = int(d[k] or 0)
            d["detect_sigma"] = float(d["detect_sigma"] or dsp.DETECT_SIGMA)
            shape = dsp.FrameShape(**d)
            v = dsp.classify_frame(shape)
            updates.append((v.verdict, v.strength_class, v.reason,
                            DISPERSION_CODE_VERSION, rid))
        ensure_morph_schema(con)
        con.executemany(
            "UPDATE frame_dispersion SET verdict_traces = ?, "
            "strength_class = ?, reason = ?, code_version = ? "
            "WHERE obs_rowid = ?", updates)
        # v1.2: the second witness, then the two-witness verdict.
        judge_morphology(con)
        # Every verdict now comes from this code, so the stage's recorded
        # code version must say so too — otherwise provenance keeps reading
        # the build-time version and the stage can never go fresh.
        for k, v in (("code_version", DISPERSION_CODE_VERSION),
                     ("reclassified_at", utcnow())):
            con.execute("INSERT OR REPLACE INTO s2c_build_meta (key, value) "
                        "VALUES (?, ?)", (k, v))
        con.commit()
        print(f"reclassify: {len(updates):,} rows re-judged under "
              f"{DISPERSION_CODE_VERSION}")
    return 0


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------
def cmd_status(args) -> int:
    con = connect(args.manifest, read_only=True)
    with closing(con):
        print(f"S2c dispersion campaign — {utcnow()}")
        print(f"{'population':<12}{'total':>9}{'measured':>10}"
              f"{'unreadable':>12}{'pending':>10}")
        for pop, tot, meas, bad, pend in con.execute(
                """SELECT population, count(*),
                          sum(status = 'measured'),
                          sum(status = 'unreadable'),
                          sum(status = 'pending')
                   FROM frame_dispersion GROUP BY population"""):
            print(f"{pop:<12}{tot:>9,}{meas:>10,}{bad:>12,}{pend:>10,}")
        print()
        print(f"{'filter':<10}{'n':>7}{'dispersed':>11}{'direct':>9}"
              f"{'indet':>8}{'unread':>8}")
        for filt, n, disp, direct, indet, unread in con.execute(
                """SELECT filter, count(*),
                          sum(verdict = 'dispersed'),
                          sum(verdict = 'direct'),
                          sum(verdict = 'indeterminate'),
                          sum(status = 'unreadable')
                   FROM frame_dispersion GROUP BY filter
                   ORDER BY count(*) DESC"""):
            print(f"{str(filt):<10}{n:>7,}{(disp or 0):>11,}"
                  f"{(direct or 0):>9,}{(indet or 0):>8,}{(unread or 0):>8,}")
    return 0


# ---------------------------------------------------------------------------
# calibrate — how well do the KNOWN labels separate?
# ---------------------------------------------------------------------------
def cmd_calibrate(args) -> int:
    """Score the classifier against the labels nobody disputes.

    Prints, for each known label, the measured distribution and the verdict
    tally.  The disputed labels ('6', 'W') are printed too but scored
    against nothing — they are the question.
    """
    import numpy as np
    con = connect(args.manifest, read_only=True)
    with closing(con):
        print(f"S2c calibration — {DISPERSION_CODE_VERSION} — {utcnow()}\n")
        print(f"{'label':<10}{'truth':<11}{'n':>7}{'agree%':>8}"
              f"{'trace_ab p50':>14}{'pa_scat p50':>13}{'medab p50':>11}")
        for filt in (dsp.KNOWN_DISPERSED_FILTERS + dsp.KNOWN_DIRECT_FILTERS
                     + ("6", "W")):
            rows = con.execute(
                """SELECT verdict, trace_ab, trace_pa_scatter, median_ab
                   FROM frame_dispersion
                   WHERE filter = ? AND status = 'measured'""",
                (filt,)).fetchall()
            if not rows:
                continue
            truth = dsp.expected_verdict(filt) or "—"
            agree = (100.0 * sum(r[0] == truth for r in rows) / len(rows)
                     if truth != "—" else float("nan"))
            def _p50(i):
                v = [r[i] for r in rows if r[i] is not None]
                return np.median(v) if v else float("nan")
            print(f"{filt:<10}{truth:<11}{len(rows):>7,}{agree:>8.1f}"
                  f"{_p50(1):>14.1f}{_p50(2):>13.1f}{_p50(3):>11.2f}")
    return 0


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def cmd_report(args) -> int:
    """Render the evidence report from the manifest (read-only)."""
    from rlmt_diagnostics.report_s2c import render_report
    path = render_report(args.manifest)
    print(f"report: wrote {path}")
    return 0


# ---------------------------------------------------------------------------
# The second witness (F-6): background morphology, truth set, SN coverage
# ---------------------------------------------------------------------------
#: Columns added to frame_dispersion by S2c v1.2.  ``verdict`` stays THE
#: verdict every consumer reads; it is now the two-witness verdict, and
#: the trace-only verdict it was before is kept beside it.
MORPH_COLUMNS = (
    ("verdict_traces", "TEXT"),   # trace-only verdict (the v1.1 rule)
    ("verdict_basis", "TEXT"),    # which witness(es) decided
    ("morph_status", "TEXT"),     # NULL (not queued) | pending | measured
                                  # | unreadable
    ("bg_corner", "REAL"), ("bg_edge_tb", "REAL"), ("bg_peak", "REAL"),
    ("bg_sigma", "REAL"), ("bg_contrast", "REAL"),
    ("bg_edge_ratio", "REAL"), ("bg_area_frac", "REAL"),
    ("bg_template_r", "REAL"),    # best lozenge-template correlation
    ("bg_template", "TEXT"),      # which template (class, shape, rotation)
    ("bg_map", "BLOB"),           # normalised 48x48 map, float16
    ("morph_class", "TEXT"), ("morph_reason", "TEXT"),
)

TRUTH_SCHEMA = """
CREATE TABLE IF NOT EXISTS s2c_truth (
    truth_idx INTEGER PRIMARY KEY,
    path      TEXT NOT NULL,
    stratum   TEXT,                -- FILTER class / v1.1 verdict at draw
    label     TEXT NOT NULL,       -- dispersed | direct | unusable
    note      TEXT
);
"""

#: The truth set: tracked source file (paths + eyeball labels).
TRUTH_CSV = PIPELINE_ROOT / "rlmt_diagnostics" / "s2c_truth_labels.csv"

#: How the truth set was DRAWN (2026-10-03, seed TRUTH_SEED) from the
#: measured rows: (stratum, predicate, size).  Kept as the record of the
#: design; the draw itself is frozen in TRUTH_CSV.  The weights follow
#: DS.F8 — populations whose verdicts nobody had checked get the most.
TRUTH_STRATA = (
    ("named_grism/dispersed", 20), ("named_grism/indeterminate", 30),
    ("named_grism/direct", 30), ("slot6/dispersed", 15),
    ("slot6/indeterminate", 20), ("slot6/direct", 20),
    ("slotW/dispersed", 8), ("slotW/indeterminate", 8), ("slotW/direct", 9),
    ("ordinary/dispersed", 20), ("ordinary/indeterminate", 8),
    ("ordinary/direct", 12))
TRUTH_SEED = 20261003

#: Template members: at most this many maps per template group, drawn by
#: a fixed seed, so a group with thousands of frames cannot make the
#: leave-one-target-out medians slow.
TEMPLATE_MAX_MEMBERS = 200
TEMPLATE_SEED = 20261004

#: Named-grism frames the trace test calls dispersed are NOT all given a
#: background map (the sky can only demote them, and the truth set
#: measures how often that is right); this many per (card class, frame
#: shape) are, to build the named-grism templates.
NAMED_TEMPLATE_SAMPLE = 60


def filter_class(filt) -> str:
    """Template group of a FILTER card: the two named units' vocabularies
    pooled, the two disputed slots, or 'ordinary'."""
    f = (filt or "").strip().lower()
    if f in ("hrg", "hagrism", "hag"):
        return "hrg"
    if f in ("lrg", "oggrism", "lrgblue"):
        return "lrg"
    return {"6": "slot6", "w": "slotW"}.get(f, "ordinary")


def ensure_morph_schema(con) -> None:
    """Add the v1.2 columns and the truth table when absent (idempotent;
    ALTER TABLE ADD COLUMN does not rewrite the table)."""
    have = {r[1] for r in con.execute("PRAGMA table_info(frame_dispersion)")}
    for col, typ in MORPH_COLUMNS:
        if col not in have:
            con.execute(f"ALTER TABLE frame_dispersion ADD COLUMN {col} {typ}")
    con.executescript(TRUTH_SCHEMA)
    con.commit()


def _morph_fields(bmap) -> dict:
    """Morphology numbers + normalised-map BLOB for one background map
    (template correlation and class are set by ``reclassify``, which
    sees the whole population)."""
    m = dsp.morphology_metrics(bmap)
    return {"morph_status": "measured", "bg_corner": m["corner"],
            "bg_edge_tb": m["edge_tb"], "bg_peak": m["peak"],
            "bg_sigma": m["sigma"], "bg_contrast": m["contrast"],
            "bg_edge_ratio": m["edge_ratio"], "bg_area_frac": m["area_frac"],
            "bg_map": dsp.map_to_blob(dsp.normalized_map(bmap))}


def cmd_truth_load(args) -> int:
    """Load the tracked truth set into ``s2c_truth`` (replaces it)."""
    import csv
    con = connect(args.manifest)
    with closing(con):
        ensure_morph_schema(con)
        with open(TRUTH_CSV, newline="") as fh:
            rows = [r for r in csv.DictReader(
                line for line in fh if not line.startswith("#"))]
        con.execute("DELETE FROM s2c_truth")
        con.executemany(
            "INSERT INTO s2c_truth VALUES (?,?,?,?,?)",
            [(int(r["truth_idx"]), r["path"], r["stratum"], r["label"],
              r["note"]) for r in rows])
        missing = con.execute(
            "SELECT count(*) FROM s2c_truth t LEFT JOIN frame_dispersion d "
            "USING (path) WHERE d.obs_rowid IS NULL").fetchone()[0]
        con.commit()
    print(f"truth-load: {len(rows)} labels; {missing} paths not in "
          "frame_dispersion")
    return 0 if missing == 0 else 2


def cmd_queue_sn(args) -> int:
    """Queue every SN 2023ixf Gate-0 frame that S2c has never measured
    (finding SN-G0d: nine in ten broadband frames had only been ASSUMED
    direct from the filter name).  Population ``sn_campaign``."""
    con = connect(args.manifest)
    with closing(con):
        ensure_morph_schema(con)
        cols = ", ".join(f"f.{c}" for c in _FRAME_COLS)
        rows = con.execute(f"""
            SELECT {cols} FROM sn_g0_frames s JOIN frames f USING (obs_rowid)
            LEFT JOIN frame_dispersion d USING (obs_rowid)
            WHERE d.obs_rowid IS NULL""").fetchall()
        con.executemany(
            f"INSERT INTO frame_dispersion ({', '.join(_FRAME_COLS)}, "
            "population, priority, status, morph_status) "
            f"VALUES ({', '.join('?' * len(_FRAME_COLS))}, 'sn_campaign', "
            "1, 'pending', 'pending')", rows)
        con.commit()
    print(f"queue-sn: {len(rows)} SN 2023ixf frames queued")
    return 0


def cmd_morph_queue(args) -> int:
    """Mark the rows whose background morphology is to be measured.

    Scope, and why (the sky can PROMOTE a frame to dispersed and DEMOTE a
    trace-dispersed frame to indeterminate):
      * every measured row whose trace verdict is NOT dispersed — the sky
        is the only witness that can promote it;
      * every slot '6' / 'W' row — the disputed slots, re-issued in full;
      * every control / holdout / SN row — the undisputed-label
        populations whose false-positive rates the report publishes;
      * every truth-set frame;
      * a fixed-seed sample of named-grism trace-dispersed frames per
        (card class, frame shape), to build their templates.
    Not queued: the remaining named-grism frames the trace test already
    calls dispersed; for them v1.2 equals v1.1 ("traces, no sky
    evidence").  The report states the count.
    """
    con = connect(args.manifest)
    with closing(con):
        ensure_morph_schema(con)
        con.execute("""
            UPDATE frame_dispersion SET morph_status = 'pending'
            WHERE status = 'measured' AND morph_status IS NULL AND (
                  coalesce(verdict_traces, verdict) != 'dispersed'
               OR lower(filter) IN ('6', 'w')
               OR population IN ('control', 'holdout', 'sn_campaign')
               OR path IN (SELECT path FROM s2c_truth))""")
        rng = random.Random(TEMPLATE_SEED)
        groups: dict = {}
        for rid, filt, h, w in con.execute("""
                SELECT obs_rowid, filter, height, width FROM frame_dispersion
                WHERE status = 'measured' AND morph_status IS NULL
                  AND coalesce(verdict_traces, verdict) = 'dispersed'
                ORDER BY obs_rowid"""):
            if filter_class(filt) in ("hrg", "lrg"):
                groups.setdefault((filter_class(filt), h, w), []).append(rid)
        pick = []
        for key in sorted(groups, key=str):
            ids = groups[key]
            pick += ids if len(ids) <= NAMED_TEMPLATE_SAMPLE else \
                rng.sample(ids, NAMED_TEMPLATE_SAMPLE)
        con.executemany("UPDATE frame_dispersion SET morph_status = "
                        "'pending' WHERE obs_rowid = ?", [(r,) for r in pick])
        con.commit()
        n = con.execute("SELECT count(*) FROM frame_dispersion WHERE "
                        "morph_status = 'pending'").fetchone()[0]
    print(f"morph-queue: {n:,} rows pending a background map "
          f"({len(pick)} named-grism template members)")
    return 0


def _morph_task(task: dict) -> tuple[int, dict]:
    """Worker: one pixel read, the background map.  Never raises."""
    try:
        from macro_grism.fits_io import load_frame
        data, _h, _l = load_frame(task["abs_path"])
        return task["obs_rowid"], _morph_fields(dsp.background_map(data))
    except Exception as exc:                      # noqa: BLE001 — recorded
        return task["obs_rowid"], {
            "morph_status": "unreadable",
            "morph_reason": f"{type(exc).__name__}: {exc}"[:200]}


def cmd_morph_run(args) -> int:
    """Measure pending background maps (resumable; one short transaction
    per chunk)."""
    started = time.time()
    con = connect(args.manifest)
    n_done = 0
    with closing(con):
        with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.workers) as pool:
            while True:
                if args.max_seconds and (time.time() - started
                                         > args.max_seconds):
                    break
                rows = con.execute(
                    "SELECT obs_rowid, path FROM frame_dispersion WHERE "
                    "morph_status = 'pending' ORDER BY priority, obs_rowid "
                    "LIMIT ?", (CHUNK,)).fetchall()
                if not rows:
                    print("morph-run: nothing pending.", flush=True)
                    break
                tasks = [{"obs_rowid": r[0],
                          "abs_path": str(args.archive / r[1])} for r in rows]
                res = list(pool.map(_morph_task, tasks, chunksize=4))
                for rid, out in res:
                    keys = list(out)
                    con.execute(
                        f"UPDATE frame_dispersion SET "
                        f"{', '.join(k + ' = ?' for k in keys)} "
                        "WHERE obs_rowid = ?", [out[k] for k in keys] + [rid])
                con.commit()
                n_done += len(res)
                print(f"  {utcnow()}  maps {n_done:,} "
                      f"({n_done / (time.time() - started):.2f}/s)",
                      flush=True)
    return 0


def judge_morphology(con) -> int:
    """Second half of ``reclassify``: lozenge templates (leave-one-target-
    out), morphology class, and the two-witness verdict, for every
    measured row.  Rows without a background map get
    combine_verdicts(trace verdict, None) — the trace verdict, with the
    basis saying there was no sky evidence."""
    import numpy as np
    rows = con.execute("""
        SELECT obs_rowid, verdict_traces, filter, height, width,
               coalesce(canonical_target, ''), bg_contrast, bg_edge_ratio,
               bg_map
        FROM frame_dispersion WHERE status = 'measured'""").fetchall()
    maps, members = {}, {}
    for rid, vt, filt, h, w, tgt, s_, e_, blob in rows:
        if blob is None:
            continue
        maps[rid] = dsp.blob_to_map(blob)
        cls = filter_class(filt)
        if (vt == dsp.VERDICT_DISPERSED and cls != "ordinary"
                and (s_ or 0) >= dsp.MORPH_TEMPLATE_MIN_CONTRAST):
            members.setdefault((cls, h, w), []).append((rid, tgt))
    rng = random.Random(TEMPLATE_SEED)
    for key in sorted(members, key=str):
        if len(members[key]) > TEMPLATE_MAX_MEMBERS:
            members[key] = rng.sample(members[key], TEMPLATE_MAX_MEMBERS)
    stacks = {k: (np.array([maps[r] for r, _ in v]),
                  np.array([t for _, t in v]), np.array([r for r, _ in v]))
              for k, v in members.items()}
    cache: dict = {}

    def templates_without(tgt: str, rid: int) -> dict:
        # A frame never votes on its own template; neither do frames of
        # its target.  (Excluding the frame itself only matters when the
        # target is blank, which is rare; it is handled exactly.)
        key = tgt if tgt else ("#", rid)
        if key not in cache:
            out = {}
            for g, (stk, tg, ids) in stacks.items():
                keep = (tg != tgt) if tgt else (ids != rid)
                out[g] = dsp.build_template(list(stk[keep]))
            if tgt:
                cache[key] = out
            else:
                return out
        return cache[key]

    updates = []
    for rid, vt, filt, h, w, tgt, s_, e_, blob in rows:
        if rid in maps:
            r, key, rot = dsp.template_correlation(
                maps[rid], templates_without(tgt, rid))
            cls, why = dsp.classify_morphology(
                {"contrast": s_, "edge_ratio": e_, "template_r": r})
            tname = None if key is None else \
                f"{key[0]} {key[2]}x{key[1]} rot{rot}"
        else:
            r, tname, cls, why = None, None, None, None
        verdict, basis = dsp.combine_verdicts(vt, cls)
        updates.append((r, tname, cls, why, verdict, basis, rid))
    con.executemany(
        "UPDATE frame_dispersion SET bg_template_r = ?, bg_template = ?, "
        "morph_class = ?, morph_reason = ?, verdict = ?, verdict_basis = ? "
        "WHERE obs_rowid = ?", updates)
    sizes = ", ".join(f"{k[0]} {k[2]}x{k[1]}: {len(v)}"
                      for k, v in sorted(members.items(), key=str))
    print(f"reclassify: template members — {sizes}")
    return len(updates)


def confusion(con, column: str) -> dict:
    """{truth label: {verdict: count}} on the truth set, for the verdict
    held in ``column`` ('verdict_traces' = v1.1, 'verdict' = v1.2)."""
    out: dict = {}
    for label, v in con.execute(f"""
            SELECT t.label, d.{column} FROM s2c_truth t
            JOIN frame_dispersion d USING (path)"""):
        out.setdefault(label, {}).setdefault(v, 0)
        out[label][v] += 1
    return out


def cmd_truth_score(args) -> int:
    """Print both confusion matrices with Wilson 95% intervals (the
    report renders the same numbers from the same function)."""
    con = connect(args.manifest, read_only=True)
    with closing(con):
        for name, col in (("v1.1 traces only", "verdict_traces"),
                          ("v1.2 traces + sky", "verdict")):
            cm = confusion(con, col)
            print(f"\n{name}")
            for label in ("dispersed", "direct", "unusable"):
                row = cm.get(label, {})
                n = sum(row.values())
                cells = []
                for v in ("dispersed", "indeterminate", "direct"):
                    k = row.get(v, 0)
                    lo, hi = dsp.wilson_interval(k, n)
                    cells.append(f"{v} {k}/{n} [{lo:.3f},{hi:.3f}]")
                print(f"  truth {label:9s}: " + "; ".join(cells))
    return 0


def cmd_morph_calibrate(args) -> int:
    """Template correlation and contrast of the undisputed-label
    populations (control/holdout = ordinary filters; named grisms called
    dispersed by the trace test) — the evidence the thresholds rest on."""
    import numpy as np
    con = connect(args.manifest, read_only=True)
    with closing(con):
        for name, pred in (
                ("ordinary filters (control+holdout)",
                 "population IN ('control','holdout')"),
                ("named grisms, trace-dispersed",
                 "lower(filter) IN ('hrg','lrg','hagrism','oggrism','hag') "
                 "AND verdict_traces = 'dispersed'")):
            v = np.array([r[0] for r in con.execute(
                f"SELECT bg_template_r FROM frame_dispersion WHERE {pred} "
                f"AND bg_template_r IS NOT NULL AND bg_contrast >= "
                f"{dsp.MORPH_MIN_CONTRAST}")])
            if len(v):
                q = np.percentile(v, [0, 1, 50, 99, 100])
                print(f"{name:40s} n={len(v):5d} template_r min/p1/p50/p99/"
                      f"max = " + "/".join(f"{x:.2f}" for x in q))
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="construct the measurement queue")
    b.add_argument("--rebuild", action="store_true",
                   help="drop all progress and requeue every frame")
    b.set_defaults(func=cmd_build)

    h = sub.add_parser("holdout",
                       help="queue a fresh, disjoint out-of-sample control")
    h.add_argument("--rebuild", action="store_true",
                   help="discard the existing holdout and draw a new one")
    h.set_defaults(func=cmd_holdout)

    r = sub.add_parser("run", help="measure pending frames (resumable)")
    r.add_argument("--workers", type=int, default=6,
                   help="worker processes (the house cap is 6)")
    r.add_argument("--limit", type=int, default=0,
                   help="stop after roughly N frames (0 = no limit)")
    r.add_argument("--max-seconds", type=float, default=0.0,
                   help="return cleanly after N seconds (0 = no limit)")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("status", help="progress and per-label tallies")
    s.set_defaults(func=cmd_status)

    rc = sub.add_parser("reclassify", help="re-judge from stored numbers")
    rc.set_defaults(func=cmd_reclassify)

    c = sub.add_parser("calibrate", help="known-label separation table")
    c.set_defaults(func=cmd_calibrate)

    rp = sub.add_parser("report", help="render the S2c evidence report")
    rp.set_defaults(func=cmd_report)

    for name, func, helptext in (
            ("truth-load", cmd_truth_load, "load the tracked truth set"),
            ("queue-sn", cmd_queue_sn, "queue unmeasured SN 2023ixf frames"),
            ("morph-queue", cmd_morph_queue, "queue background maps"),
            ("morph-calibrate", cmd_morph_calibrate,
             "template_r of the labelled populations"),
            ("truth-score", cmd_truth_score, "confusion matrices")):
        sp = sub.add_parser(name, help=helptext)
        sp.set_defaults(func=func)
    mr = sub.add_parser("morph-run", help="measure pending background maps")
    mr.add_argument("--workers", type=int, default=4)
    mr.add_argument("--max-seconds", type=float, default=0.0)
    mr.set_defaults(func=cmd_morph_run)

    args = p.parse_args(argv)
    # Guard the house rule: never more than six workers against this disk.
    if getattr(args, "workers", 0) > 6:
        print("run: refusing more than 6 workers (shared disk).")
        return 2
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
