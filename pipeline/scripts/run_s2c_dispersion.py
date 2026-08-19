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
    calibrate   print the known-label separation table (read-only)

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
        shape = dsp.measure_file(task["abs_path"])
        verdict = dsp.classify_frame(shape)
        out.update(shape.as_dict())
        out["verdict"] = verdict.verdict
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
                "measure_s", "error", "code_version", "measured_at"]


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
        con.executemany(
            "UPDATE frame_dispersion SET verdict = ?, strength_class = ?, "
            "reason = ?, code_version = ? WHERE obs_rowid = ?", updates)
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

    args = p.parse_args(argv)
    # Guard the house rule: never more than six workers against this disk.
    if getattr(args, "workers", 0) > 6:
        print("run: refusing more than 6 workers (shared disk).")
        return 2
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
