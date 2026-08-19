#!/usr/bin/env python
"""Run the S1b production astrometry batch over the stratified backlog.

WHAT THIS SCRIPT DOES (stage S1b of the shared pipeline)
--------------------------------------------------------
The S1 go/no-go experiment (``run_s1_experiment.py``, report
``docs/pipeline/s1_astrometry.html``) measured per-stratum solve rates
and the review accepted the verdicts.  This script executes the batch
those verdicts authorize: solve the ~38k stratified backlog frames with
the SAME wrapper, the SAME per-frame caps and the SAME acceptance gate
as the experiment, in priority order (CV polars → dwarf → QC-gated SN →
facility backlog), resumably, writing:

* one NEW manifest table, ``s1_batch`` — one row per queued frame with
  its status (pending | solved | bad_solve | failed | skipped_qc), the
  accepted WCS numbers, and (for QC-gated strata) the pre-gate
  statistics.  Plus ``s1_batch_meta`` build facts.  Existing tables are
  never modified.
* sidecar WCS files under ``products/astrom/wcs/`` mirroring the
  archive tree — the archive itself is NEVER written to.

Below-GO strata (SN 2023ixf, dwarf, cv_gsense_misc, fast_fullframe) run
behind the QC pre-gate: the experiment's failure-autopsy statistics
classify each frame first, and starved / defocused / trailed / flooded
frames are recorded as ``skipped_qc`` (with their diagnosis) instead of
burning the 75 s solve budget.  Policy + thresholds live in
``macro_core.batch`` / ``macro_core.astrom`` — nothing is decided here.

SUBCOMMANDS
-----------
    build     construct the queue from the manifest (refuses to clobber
              an existing queue unless --rebuild)
    run       solve pending frames in priority order; SAFE TO RE-RUN —
              a killed run loses nothing but in-flight frames, which
              stay pending and are picked up on resume
    status    per-population progress table + ETA (read-only)
    verify    cross-check s1_batch rows against the sidecar WCS files

USAGE (a student's quick start)
-------------------------------
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/run_s1_batch.py build
    $PY pipeline/scripts/run_s1_batch.py run            # until done
    $PY pipeline/scripts/run_s1_batch.py status         # any time
    $PY pipeline/scripts/run_s1_batch.py verify         # afterwards

The full run takes hours: launch it detached (nohup/caffeinate) and
watch the log with tail -f.  ``run --limit N`` and ``run --smoke N``
(N frames spread across every stratum) exist for testing.
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

# Make the pipeline package importable no matter where the script is
# invoked from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import astrom, batch                       # noqa: E402
from macro_core.batch import S1B_CODE_VERSION              # noqa: E402


# How long to keep retrying a write that loses the race for the manifest's
# write lock, and how long to wait between attempts.  The connection's own
# ``timeout=`` already waits for the lock, but that wait is per-attempt and
# raises OperationalError once it expires — which on 2026-08-18 killed a
# multi-hour batch outright after some other stage held the lock past the
# 120 s ceiling.  A long job must never die because a sibling stage was
# briefly slow, so the write is retried with backoff and only gives up when
# the lock has genuinely been held for many minutes.
WRITE_RETRIES = 8
WRITE_BACKOFF_S = 15.0


def execute_resilient(con: sqlite3.Connection, sql: str, params) -> None:
    """Run one write, tolerating a manifest lock held by another stage.

    Retries only on SQLite's lock errors ("database is locked" / "database is
    busy"); every other OperationalError is a real fault and propagates
    immediately, so genuine bugs still fail loudly.
    """
    for attempt in range(1, WRITE_RETRIES + 1):
        try:
            con.execute(sql, params)
            return
        except sqlite3.OperationalError as exc:
            msg = str(exc).lower()
            if "locked" not in msg and "busy" not in msg:
                raise                       # not contention — a real error
            if attempt == WRITE_RETRIES:
                raise
            # Linear backoff: the holder is usually a bulk insert from a
            # sibling stage, which finishes in tens of seconds, not hours.
            time.sleep(WRITE_BACKOFF_S)

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
DEFAULT_INDEX_DIR = Path.home() / "astrometry-indices"
DEFAULT_CONFIG = DEFAULT_INDEX_DIR / "astrometry.cfg"
WCS_ROOT = REPO_ROOT / "products" / "astrom" / "wcs"


def utcnow() -> str:
    """ISO-8601 UTC timestamp for log lines and DB facts."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# ---------------------------------------------------------------------------
# build — construct the queue table from the manifest
# ---------------------------------------------------------------------------
def cmd_build(args) -> int:
    con = sqlite3.connect(args.manifest, timeout=60)
    with closing(con):
        exists = con.execute(
            "SELECT count(*) FROM sqlite_master "
            "WHERE type='table' AND name='s1_batch'").fetchone()[0]
        if exists and not args.rebuild:
            n_done = con.execute(
                "SELECT count(*) FROM s1_batch "
                "WHERE status != 'pending'").fetchone()[0]
            print(f"build: s1_batch already exists ({n_done:,} frames "
                  "finished) — refusing to clobber progress. "
                  "Pass --rebuild to drop and requeue everything.")
            return 1
        # The queue = the experiment's candidate universe, classified by
        # the SAME pure gates and stratum rules, with the batch policy
        # (population, priority, QC gating) attached.  Residue frames in
        # no stratum are excluded — no measured rate, no batch.
        rows = astrom.fetch_candidates(con)
        queue = batch.build_queue_rows(rows)
        # Tooling facts recorded for the eventual report.
        try:
            sf_version = subprocess.run(
                ["solve-field", "--version"], capture_output=True,
                text=True, timeout=10).stdout.strip()
        except Exception:
            sf_version = "unavailable"
        git = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True,
                             cwd=REPO_ROOT).stdout.strip()
        cur = con.cursor()
        cur.execute("BEGIN")
        cur.execute("DROP TABLE IF EXISTS s1_batch")
        cur.execute("DROP TABLE IF EXISTS s1_batch_meta")
        cur.execute("""CREATE TABLE s1_batch (
            obs_rowid INTEGER PRIMARY KEY,
            stratum_id TEXT NOT NULL,
            population TEXT NOT NULL,
            priority INTEGER NOT NULL,
            qc_gated INTEGER NOT NULL,
            path TEXT NOT NULL,
            night TEXT, target_key TEXT, canonical_target TEXT,
            readoutm TEXT, xbinning REAL, filter TEXT, exptime REAL,
            ra_hint_deg REAL, dec_hint_deg REAL,
            status TEXT NOT NULL DEFAULT 'pending',
            fail_kind TEXT,
            qc_n_sources INTEGER, qc_n_psf_sources INTEGER,
            qc_median_elongation REAL, qc_bright_median_a_px REAL,
            qc_saturated_fraction REAL, qc_diagnosis TEXT,
            used_hint INTEGER, solve_time_s REAL,
            solved_ra REAL, solved_dec REAL,
            pixscale_arcsec REAL, rotation_deg REAL,
            n_matched INTEGER, rms_arcsec REAL,
            wcs_path TEXT,
            log_tail TEXT,
            finished_utc TEXT)""")
        cur.execute("CREATE INDEX ix_s1b_status ON s1_batch(status)")
        cur.execute("CREATE INDEX ix_s1b_pop ON s1_batch(population)")
        cur.execute(
            "CREATE INDEX ix_s1b_queue ON s1_batch(status, priority)")
        cur.executemany("""INSERT INTO s1_batch
            (obs_rowid, stratum_id, population, priority, qc_gated,
             path, night, target_key, canonical_target, readoutm,
             xbinning, filter, exptime, ra_hint_deg, dec_hint_deg,
             status)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,'pending')""",
            [(q["obs_rowid"], q["stratum_id"], q["population"],
              q["priority"], q["qc_gated"], q["path"], q["night"],
              q["target_key"], q["canonical_target"], q["readoutm"],
              q["xbinning"], q["filter"], q["exptime"],
              q["ra_deg"], q["dec_deg"]) for q in queue])
        cur.execute("""CREATE TABLE s1_batch_meta
                       (key TEXT PRIMARY KEY, value TEXT)""")
        cur.executemany("INSERT INTO s1_batch_meta VALUES (?,?)", [
            ("built_utc", utcnow()),
            ("code_version", S1B_CODE_VERSION),
            ("git_commit", git),
            ("solve_field_version", sf_version),
            ("config", str(args.config)),
            ("solve_timeout_s", str(astrom.SOLVE_TIMEOUT_S)),
            ("solve_cpulimit_s", str(astrom.SOLVE_CPU_LIMIT_S)),
            ("hint_radius_deg", str(astrom.HINT_RADIUS_DEG)),
            ("downsample", str(astrom.SOLVE_DOWNSAMPLE)),
            ("wcs_root", str(WCS_ROOT)),
            ("n_queued", str(len(queue))),
        ])
        con.commit()
        print(f"build: {len(queue):,} frames queued "
              f"({len(rows):,} unsolved candidates scanned)")
        for pop in batch.POPULATIONS:
            sids = [p.stratum_id for p in batch.STRATUM_POLICY
                    if p.population == pop]
            n = sum(1 for q in queue if q["population"] == pop)
            n_gated = sum(1 for q in queue
                          if q["population"] == pop and q["qc_gated"])
            print(f"  {pop:<12} {n:>7,} frames "
                  f"({n_gated:,} QC-gated)  strata: {', '.join(sids)}")
    return 0


# ---------------------------------------------------------------------------
# run — the worker task (QC gate, solve, sidecar) + the resumable loop
# ---------------------------------------------------------------------------
def _read_pixels(src: Path):
    """Pixel array of one archive frame (fpack puts data in HDU 1)."""
    import numpy as np
    from astropy.io import fits
    with fits.open(src) as hdul:
        hdu = hdul[1] if len(hdul) > 1 and hdul[0].data is None \
            else hdul[0]
        return np.asarray(hdu.data, dtype=np.float32)


def _write_sidecar(scratch_wcs: Path, rel_path: str) -> str:
    """Atomically install one accepted .wcs into the sidecar tree.

    Returns the sidecar path RELATIVE to the wcs root (what s1_batch
    stores).  Atomic: copy to a dot-tmp neighbor, fsync, rename — a
    killed run can never leave a half-written sidecar under its final
    name (resume overwrites whole files the same way).
    """
    rel = batch.sidecar_rel_path(rel_path)
    dest = WCS_ROOT / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.parent / f".{dest.name}.tmp"
    shutil.copyfile(scratch_wcs, tmp)
    with open(tmp, "rb") as fh:            # flush the copy to disk
        os.fsync(fh.fileno())
    os.replace(tmp, dest)                  # atomic on one filesystem
    return rel


def _frame_task(row: dict, archive_root: Path, config: str,
                scratch_root: Path) -> tuple[int, dict]:
    """Worker: QC-gate (if gated) then solve ONE frame; returns a plain
    result dict — the main thread is the only DB writer.

    The result dict always carries ``status`` (a terminal status from
    ``batch.TERMINAL_STATUSES``) plus whichever evidence fields the path
    taken produced (QC statistics, WCS numbers, sidecar path, log tail).
    """
    src = archive_root / row["path"]
    res: dict = {"status": "failed", "fail_kind": "error",
                 "log_tail": None}
    # -- 1. QC pre-gate (below-GO strata only) --------------------------
    if row["qc_gated"]:
        try:
            data = _read_pixels(src)
            m = astrom.image_metrics(
                data, batch.saturation_adu_for(row["readoutm"]))
            passed, diagnosis = batch.qc_pregate(m)
            res.update({
                "qc_n_sources": m["n_sources"],
                "qc_n_psf_sources": m["n_psf_sources"],
                "qc_median_elongation": m["median_elongation"],
                "qc_bright_median_a_px": m["bright_median_a_px"],
                "qc_saturated_fraction": m["saturated_fraction"],
                "qc_diagnosis": diagnosis})
            if not passed:
                res["status"] = "skipped_qc"
                res["fail_kind"] = None
                return row["obs_rowid"], res
        except Exception as exc:           # noqa: BLE001 — recorded
            res["status"] = "skipped_qc"
            res["fail_kind"] = None
            res["qc_diagnosis"] = \
                f"unreadable: {type(exc).__name__}: {exc}"[:200]
            return row["obs_rowid"], res
    # -- 2. solve in a private scratch dir, same caps as the experiment -
    work = Path(tempfile.mkdtemp(prefix=f"s1b_{row['obs_rowid']}_",
                                 dir=scratch_root))
    try:
        s = astrom.solve_one_frame(
            str(src), str(work), config,
            row["readoutm"], row["xbinning"],
            row["ra_hint_deg"], row["dec_hint_deg"])
        # Map the experiment's status vocabulary onto the batch's:
        # solved / bad_solve pass through; unsolved / timeout / error
        # collapse to 'failed' with the detail kept in fail_kind.
        if s["status"] in ("solved", "bad_solve"):
            res["status"] = s["status"]
            res["fail_kind"] = None
        else:
            res["status"] = "failed"
            res["fail_kind"] = s["status"]
        res.update({k: s[k] for k in (
            "used_hint", "solve_time_s", "solved_ra", "solved_dec",
            "pixscale_arcsec", "rotation_deg", "n_matched",
            "rms_arcsec", "log_tail")})
        # -- 3. sidecar: ONLY gate-accepted solutions leave scratch ----
        if res["status"] == "solved":
            wcs_file = work / batch.scratch_wcs_name(row["path"])
            res["wcs_path"] = _write_sidecar(wcs_file, row["path"])
        return row["obs_rowid"], res
    except Exception as exc:               # noqa: BLE001 — recorded
        res["status"] = "failed"
        res["fail_kind"] = "error"
        res["log_tail"] = f"{type(exc).__name__}: {exc}"[:400]
        return row["obs_rowid"], res
    finally:
        shutil.rmtree(work, ignore_errors=True)


#: Columns a worker needs — fetched per chunk, passed as plain dicts.
_TASK_COLS = ["obs_rowid", "path", "stratum_id", "population",
              "qc_gated", "readoutm", "xbinning",
              "ra_hint_deg", "dec_hint_deg"]

#: Every result field _frame_task may produce, in UPDATE order.
_RESULT_COLS = ["status", "fail_kind",
                "qc_n_sources", "qc_n_psf_sources",
                "qc_median_elongation", "qc_bright_median_a_px",
                "qc_saturated_fraction", "qc_diagnosis",
                "used_hint", "solve_time_s", "solved_ra", "solved_dec",
                "pixscale_arcsec", "rotation_deg", "n_matched",
                "rms_arcsec", "wcs_path", "log_tail"]


def _fetch_chunk(con, n: int, populations, smoke: bool) -> list[dict]:
    """Next ``n`` pending frames.

    Normal mode: strict priority order — the paper-critical strata
    drain first.  Smoke mode: round-robin across strata (lowest rowids
    of each) so a small test batch exercises EVERY population, QC gate
    included, instead of chewing the head of the CV queue only.
    """
    pop_clause, params = "", []
    if populations:
        pop_clause = ("AND population IN (%s)"
                      % ",".join("?" * len(populations)))
        params = list(populations)
    if not smoke:
        rows = con.execute(
            f"""SELECT {', '.join(_TASK_COLS)} FROM s1_batch
                WHERE status = 'pending' {pop_clause}
                ORDER BY priority, obs_rowid LIMIT ?""",
            (*params, n)).fetchall()
        return [dict(zip(_TASK_COLS, r)) for r in rows]
    # Smoke: rank within each stratum, interleave ranks across strata.
    rows = con.execute(
        f"""SELECT {', '.join(_TASK_COLS)} FROM (
                SELECT *, ROW_NUMBER() OVER (
                    PARTITION BY stratum_id ORDER BY obs_rowid) AS rk
                FROM s1_batch WHERE status = 'pending' {pop_clause})
            ORDER BY rk, priority LIMIT ?""",
        (*params, n)).fetchall()
    return [dict(zip(_TASK_COLS, r)) for r in rows]


def cmd_run(args) -> int:
    t_start = time.monotonic()
    scratch_root = Path(tempfile.mkdtemp(prefix="s1b_run_"))
    con = sqlite3.connect(args.manifest, timeout=120)
    # WAL lets `status` read while the run commits — and a killed run
    # recovers to its last committed chunk, nothing half-written.
    con.execute("PRAGMA journal_mode=WAL")
    budget = args.batch_seconds if args.batch_seconds > 0 else float("inf")
    limit = args.limit if args.limit > 0 else float("inf")
    n_done = 0
    counts: dict[str, int] = {}
    print(f"[{utcnow()}] run: workers={args.workers} "
          f"limit={args.limit or 'none'} "
          f"populations={','.join(args.populations) or 'all'} "
          f"smoke={bool(args.smoke)}", flush=True)
    try:
        while (time.monotonic() - t_start) < budget and n_done < limit:
            # Small chunks: a killed run loses at most one chunk of
            # uncommitted results — those frames stay pending, resume
            # re-solves them (sidecar overwrite is atomic + idempotent).
            want = int(min(args.workers * 3, limit - n_done))
            rows = _fetch_chunk(con, want, args.populations,
                                bool(args.smoke))
            if not rows:
                print(f"[{utcnow()}] run: nothing pending — done",
                      flush=True)
                break
            by_id = {r["obs_rowid"]: r for r in rows}
            with concurrent.futures.ThreadPoolExecutor(args.workers) as ex:
                futs = [ex.submit(_frame_task, r, args.archive,
                                  str(args.config), scratch_root)
                        for r in rows]
                for fut in concurrent.futures.as_completed(futs):
                    rowid, res = fut.result()
                    # Single-writer discipline: only this thread
                    # touches the DB.  The WHERE guards the transition
                    # contract — only pending rows may be finished.
                    sets = ", ".join(f"{c}=?" for c in _RESULT_COLS)
                    execute_resilient(
                        con,
                        f"""UPDATE s1_batch SET {sets}, finished_utc=?
                            WHERE obs_rowid=? AND status='pending'""",
                        [res.get(c) for c in _RESULT_COLS]
                        + [utcnow(), rowid])
                    n_done += 1
                    st = res["status"]
                    counts[st] = counts.get(st, 0) + 1
                    # One line per frame — tail -f–friendly progress.
                    t = res.get("solve_time_s")
                    t_str = f" t={t:.1f}s" if t is not None else ""
                    extra = res.get("qc_diagnosis") or \
                        res.get("fail_kind") or ""
                    print(f"[{utcnow()}] {st:<10} "
                          f"{by_id[rowid]['stratum_id']:<22} "
                          f"rowid={rowid}{t_str}"
                          + (f"  [{extra}]" if st != "solved" and extra
                             else ""), flush=True)
            con.commit()
            n_pend = con.execute("SELECT count(*) FROM s1_batch "
                                 "WHERE status='pending'").fetchone()[0]
            summary = " ".join(f"{k}={v}" for k, v in sorted(counts.items()))
            print(f"[{utcnow()}] chunk done: {n_done} this run "
                  f"({summary}); {n_pend:,} pending", flush=True)
        print(f"[{utcnow()}] run: finished ({n_done} frames, "
              f"{time.monotonic() - t_start:.0f}s)", flush=True)
    finally:
        con.commit()
        con.close()
        shutil.rmtree(scratch_root, ignore_errors=True)
    return 0


# ---------------------------------------------------------------------------
# status — per-population progress + ETA (read-only)
# ---------------------------------------------------------------------------
def _median(vals: list) -> float | None:
    """Plain median (interpolating) without importing numpy for a CLI."""
    v = sorted(x for x in vals if x is not None)
    if not v:
        return None
    mid = len(v) // 2
    return v[mid] if len(v) % 2 else (v[mid - 1] + v[mid]) / 2.0


def cmd_status(args) -> int:
    con = sqlite3.connect(f"file:{args.manifest}?mode=ro", uri=True,
                          timeout=30)
    with closing(con):
        try:
            con.execute("SELECT 1 FROM s1_batch LIMIT 1")
        except sqlite3.OperationalError:
            print("status: no s1_batch table — run `build` first")
            return 1
        print(f"{'population':<12} {'stratum':<22} {'total':>7} "
              f"{'pend':>7} {'solved':>7} {'skipQC':>6} {'badWCS':>6} "
              f"{'failed':>6} {'rate':>6} {'med_t':>6}")
        pend_by_stratum: dict[str, int] = {}
        med_by_stratum: dict[str, float] = {}
        for (pop, sid, tot, pend, ok, skip, bad, fail) in con.execute("""
                SELECT population, stratum_id, count(*),
                       sum(status='pending'), sum(status='solved'),
                       sum(status='skipped_qc'), sum(status='bad_solve'),
                       sum(status='failed')
                FROM s1_batch GROUP BY population, stratum_id
                ORDER BY min(priority)"""):
            times = [r[0] for r in con.execute(
                "SELECT solve_time_s FROM s1_batch WHERE stratum_id=? "
                "AND solve_time_s IS NOT NULL", (sid,))]
            live_med = _median(times)
            med = batch.stratum_median_s(sid, live_med, len(times))
            pend_by_stratum[sid] = pend
            med_by_stratum[sid] = med
            attempted = ok + bad + fail       # QC skips are not attempts
            rate = f"{100.0 * ok / attempted:.0f}%" if attempted else "-"
            print(f"{pop:<12} {sid:<22} {tot:>7,} {pend:>7,} {ok:>7,} "
                  f"{skip:>6,} {bad:>6,} {fail:>6,} {rate:>6} "
                  f"{med:>5.1f}s")
        n_pend = sum(pend_by_stratum.values())
        eta_s = batch.eta_seconds(pend_by_stratum, med_by_stratum,
                                  args.workers)
        n_wcs = con.execute("SELECT count(*) FROM s1_batch "
                            "WHERE status='solved'").fetchone()[0]
        print(f"\ntotal: {n_pend:,} pending, {n_wcs:,} accepted WCS; "
              f"ETA ~{eta_s / 3600.0:.1f} h at {args.workers} workers "
              f"(strata medians: live when >= "
              f"{batch.MIN_FRAMES_FOR_LIVE_MEDIAN} frames, else the S1 "
              "experiment's priors)")
    return 0


# ---------------------------------------------------------------------------
# verify — s1_batch rows and sidecar files must agree
# ---------------------------------------------------------------------------
def cmd_verify(args) -> int:
    con = sqlite3.connect(f"file:{args.manifest}?mode=ro", uri=True,
                          timeout=30)
    with closing(con):
        # Direction 1: every 'solved' row names a sidecar that exists.
        missing = []
        n_solved = 0
        for (rowid, wcs_rel) in con.execute(
                "SELECT obs_rowid, wcs_path FROM s1_batch "
                "WHERE status='solved'"):
            n_solved += 1
            if wcs_rel is None or not (WCS_ROOT / wcs_rel).is_file():
                missing.append(rowid)
        # Direction 2: every sidecar on disk belongs to a 'solved' row
        # (an orphan from a killed run is benign — its row is pending
        # and will be re-solved — but it must be REPORTED, not hidden).
        known = {r[0] for r in con.execute(
            "SELECT wcs_path FROM s1_batch WHERE wcs_path IS NOT NULL")}
        orphans = []
        if WCS_ROOT.is_dir():
            for p in WCS_ROOT.rglob("*.wcs"):
                if str(p.relative_to(WCS_ROOT)) not in known:
                    orphans.append(str(p))
        print(f"verify: {n_solved:,} solved rows; "
              f"{len(missing)} missing sidecars; "
              f"{len(orphans)} orphan sidecars")
        for rowid in missing[:20]:
            print(f"  MISSING sidecar for solved row {rowid}")
        for p in orphans[:20]:
            print(f"  ORPHAN {p}")
        return 0 if not missing else 1


# ---------------------------------------------------------------------------
# CLI plumbing
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="S1b production astrometry batch driver",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    p.add_argument("--workers", type=int, default=astrom.DEFAULT_WORKERS)
    sub = p.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="construct the queue table")
    b.add_argument("--rebuild", action="store_true",
                   help="drop an existing s1_batch and requeue")
    r = sub.add_parser("run", help="solve pending frames (resumable)")
    r.add_argument("--batch-seconds", type=int, default=0,
                   help="stop after this many seconds (0 = run to empty)")
    r.add_argument("--limit", type=int, default=0,
                   help="stop after this many frames (0 = no cap)")
    r.add_argument("--smoke", type=int, default=0, metavar="N",
                   help="test mode: N frames spread across ALL strata "
                        "instead of strict priority order")
    r.add_argument("--populations", type=lambda s: s.split(","),
                   default=[], metavar="P1,P2",
                   help="restrict to these populations "
                        f"({','.join(batch.POPULATIONS)})")
    e = sub.add_parser("enqueue",
                       help="ADD newly-solvable frames without dropping the "
                            "queue (non-destructive; use after an S0 rebuild)")
    e.add_argument("--dry-run", action="store_true",
                   help="report what would be added and change nothing")
    sub.add_parser("status", help="progress + ETA (read-only)")
    sub.add_parser("verify", help="cross-check rows vs sidecar files")
    args = p.parse_args(argv)
    if args.command == "run" and args.smoke:
        # --smoke N is sugar for a mixed --limit N run.
        args.limit = args.smoke
    if args.command == "run":
        bad = set(args.populations) - set(batch.POPULATIONS)
        if bad:
            p.error(f"unknown population(s): {', '.join(sorted(bad))}")
    return args


def cmd_enqueue(args) -> int:
    """Add frames that have BECOME solvable, without touching finished work.

    ``build --rebuild`` DROPs ``s1_batch``, which throws away every solved,
    failed and skipped_qc verdict the batch has earned.  That is the right
    behaviour when the queue design changes, and the wrong behaviour when
    the frame universe merely GREW — which is exactly what the S0e geometry
    repair did: 18k frames that the solvability gate rejected on phantom
    8-pixel geometry are full fields and were never unsolvable.

    So this subcommand is purely additive.  It recomputes the candidate
    universe with the same pure gates, and INSERTs only those queue rows
    whose ``obs_rowid`` is not already present, as ``pending``.  Rows
    already in the table — finished or pending — are left exactly as they
    are (``INSERT OR IGNORE`` on the primary key guarantees it).

    Frames in no stratum are reported and NOT queued, same as ``build``:
    nothing outside a measured stratum gets solved.  That report is the
    point — after the S0e repair a large block of EU UMa frames lands in
    that gap, and silently dropping them is how the artifact stayed
    invisible the first time.
    """
    con = sqlite3.connect(args.manifest, timeout=60)
    with closing(con):
        con.execute("PRAGMA busy_timeout = 300000")
        exists = con.execute(
            "SELECT count(*) FROM sqlite_master "
            "WHERE type='table' AND name='s1_batch'").fetchone()[0]
        if not exists:
            print("enqueue: no s1_batch table — run `build` first.")
            return 1
        rows = astrom.fetch_candidates(con)
        queue = batch.build_queue_rows(rows)
        # The residue: solvable frames that classify_stratum leaves unlabelled.
        n_solvable = sum(1 for r in rows if astrom.is_solvable_candidate(r))
        n_residue = n_solvable - len(queue)
        have = {r[0] for r in con.execute("SELECT obs_rowid FROM s1_batch")}
        fresh = [q for q in queue if q["obs_rowid"] not in have]
        print(f"enqueue: {len(rows):,} candidates, {n_solvable:,} solvable, "
              f"{len(queue):,} in a stratum, {len(have):,} already queued")
        print(f"enqueue: {len(fresh):,} NEW frames to add")
        if n_residue:
            print(f"enqueue: WARNING — {n_residue:,} solvable frames are in "
                  "NO stratum and will NOT be queued")
        by_stratum = {}
        for q in fresh:
            by_stratum[q["stratum_id"]] = by_stratum.get(q["stratum_id"], 0) + 1
        for sid, n in sorted(by_stratum.items(), key=lambda kv: -kv[1]):
            print(f"    {sid:24s} +{n:,}")
        # Ordering metadata that has DRIFTED from the policy table.  Adding a
        # stratum renumbers the ranks below it, and s1_batch stores the rank
        # on every row — so without this, two frames of one stratum could
        # carry different ranks purely because of when each was enqueued, and
        # the queue order would stop being a function of the policy at all.
        # Only priority / population / qc_gated are re-derived: these are
        # ordering metadata, never results.  Status, WCS and timings are
        # untouched.
        stored = [{"obs_rowid": r[0], "population": r[1], "priority": r[2],
                   "qc_gated": r[3]}
                  for r in con.execute("SELECT obs_rowid, population, "
                                       "priority, qc_gated FROM s1_batch")]
        drift = [(d["priority"], d["population"], d["qc_gated"],
                  d["obs_rowid"]) for d in batch.policy_drift(stored, queue)]
        if drift:
            print(f"enqueue: {len(drift):,} existing rows carry stale "
                  "priority/population/gating — re-syncing from the policy")
        if args.dry_run:
            print("enqueue: --dry-run, nothing written")
            return 0
        if drift:
            con.executemany(
                "UPDATE s1_batch SET priority = ?, population = ?, "
                "qc_gated = ? WHERE obs_rowid = ?", drift)
            con.commit()
        if not fresh:
            return 0
        cols = ["obs_rowid", "stratum_id", "population", "priority",
                "qc_gated", "path", "night", "target_key", "canonical_target",
                "readoutm", "xbinning", "filter", "exptime",
                "ra_hint_deg", "dec_hint_deg"]
        con.executemany(
            "INSERT OR IGNORE INTO s1_batch (%s) VALUES (%s)"
            % (", ".join(cols), ",".join("?" * len(cols))),
            [[q.get("ra_deg") if c == "ra_hint_deg" else
              q.get("dec_deg") if c == "dec_hint_deg" else q.get(c)
              for c in cols] for q in fresh])
        con.commit()
        total = con.execute("SELECT count(*) FROM s1_batch").fetchone()[0]
        pend = con.execute("SELECT count(*) FROM s1_batch "
                           "WHERE status='pending'").fetchone()[0]
        print(f"enqueue: done — s1_batch now {total:,} rows, {pend:,} pending")
    return 0


def main(argv=None) -> int:
    args = parse_args(argv)
    return {"build": cmd_build, "run": cmd_run, "enqueue": cmd_enqueue,
            "status": cmd_status, "verify": cmd_verify}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
