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
                            and why: grism spectra, 8-px photometry windows)
* ``s1_solve_experiment`` — one row per sampled frame with the solve outcome
* ``s1_failure_autopsy``  — image-statistics post-mortem of failures
* ``s1_build_meta``       — timestamp, versions, tooling inventory

SUBCOMMANDS (run in this order)
-------------------------------
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
# design — build strata, sample, write the s1 design tables
# ---------------------------------------------------------------------------
def cmd_design(args) -> int:
    con = sqlite3.connect(args.manifest)
    with closing(con):
        rows = fetch_candidates(con)
        # Census of the exclusions: the report must show what astrometry
        # CANNOT apply to, with counts, before showing what it can.
        n_total = len(rows)
        n_grism = sum(astrom.is_grism_filter(r["filter"]) for r in rows)
        n_vocab = sum((not astrom.is_grism_filter(r["filter"]))
                      and astrom.is_calib_vocab_filter(r["filter"])
                      for r in rows)
        n_window = sum((not astrom.is_grism_filter(r["filter"]))
                       and (not astrom.is_calib_vocab_filter(r["filter"]))
                       and astrom.is_window_geometry(r["naxis1"],
                                                     r["naxis2"])
                       for r in rows)
        candidates = [r for r in rows if astrom.is_solvable_candidate(r)]
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
        # Draw the reproducible sample per stratum.
        sample_rows = []
        strata_rows = []
        for s in astrom.STRATA:
            pop = by_stratum[s.stratum_id]
            ids = {r["obs_rowid"]: r for r in pop}
            picked = astrom.sample_frames(list(ids), args.n_per_stratum,
                                          astrom.SAMPLE_SEED, s.stratum_id)
            strata_rows.append((s.stratum_id, s.population, s.description,
                                len(pop), len(picked), astrom.SAMPLE_SEED))
            for order, rowid in enumerate(picked):
                r = ids[rowid]
                sample_rows.append((
                    r["obs_rowid"], s.stratum_id, r["path"],
                    r["target_key"], r["canonical_target"], r["readoutm"],
                    r["xbinning"], r["filter"], r["exptime"],
                    r["naxis1"], r["naxis2"], r["ra_deg"], r["dec_deg"],
                    r["night"], order, "pending"))
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
        cur.execute("BEGIN")
        for t in ("s1_strata", "s1_populations", "s1_solve_experiment",
                  "s1_failure_autopsy", "s1_build_meta"):
            cur.execute(f"DROP TABLE IF EXISTS {t}_s1_tmp")
        cur.execute("""CREATE TABLE s1_strata_s1_tmp (
            stratum_id TEXT PRIMARY KEY, population TEXT, description TEXT,
            n_population INTEGER, n_sample INTEGER, seed INTEGER)""")
        cur.executemany("INSERT INTO s1_strata_s1_tmp VALUES (?,?,?,?,?,?)",
                        strata_rows)
        cur.execute("""CREATE TABLE s1_populations_s1_tmp (
            class TEXT PRIMARY KEY, n_frames INTEGER, note TEXT)""")
        cur.executemany(
            "INSERT INTO s1_populations_s1_tmp VALUES (?,?,?)",
            [("unsolved_total", n_total,
              "canonical rawimage Light frames with pltsolvd != 1"),
             ("excluded_grism", n_grism,
              "FILTER names a grism — slitless spectra, never solvable"),
             ("excluded_calib_vocab_filter", n_vocab,
              "FILTER = 'dark'/'bias'/'flat' header glitch (S0b)"),
             ("excluded_window_geometry", n_window,
              f"either axis < {astrom.MIN_SOLVABLE_NAXIS} px — high-speed "
              "photometry windows (8×3211 strips), no sky for quads"),
             ("solvable_candidates", len(candidates),
              "frames the solver can be pointed at"),
             ("candidates_unstratified", n_residue,
              "solvable but in no stratum (small heterogeneous residue)")])
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
             ra_hint_deg, dec_hint_deg, night, sample_order, status)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", sample_rows)
        cur.execute("""CREATE TABLE s1_failure_autopsy_s1_tmp (
            obs_rowid INTEGER PRIMARY KEY, stratum_id TEXT,
            n_sources INTEGER, n_psf_sources INTEGER,
            median_elongation REAL, median_a_px REAL,
            bright_median_a_px REAL, saturated_fraction REAL,
            bkg_rms REAL, saturation_adu REAL,
            diagnosis TEXT, thumb_path TEXT)""")
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
        ])
        for t in ("s1_strata", "s1_populations", "s1_solve_experiment",
                  "s1_failure_autopsy", "s1_build_meta"):
            cur.execute(f"DROP TABLE IF EXISTS {t}")
            cur.execute(f"ALTER TABLE {t}_s1_tmp RENAME TO {t}")
        con.commit()
        n_samp = len(sample_rows)
        print(f"design: {n_total:,} unsolved -> {len(candidates):,} "
              f"solvable candidates -> {n_samp} sampled frames across "
              f"{len(astrom.STRATA)} strata (seed {astrom.SAMPLE_SEED})")
        for row in strata_rows:
            print(f"  {row[0]:<24} pop {row[3]:>6,}  sample {row[4]:>3}")
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
    con = sqlite3.connect(args.manifest, timeout=60)
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
    con = sqlite3.connect(args.manifest, timeout=60)
    with closing(con):
        # EVERY failure is autopsied — no per-stratum cap.  (v1.0 capped
        # at 8 per stratum, which left 25 of 58 failures undiagnosed in
        # exactly the NO-GO strata where the taxonomy matters most.)
        # ``bad_solve`` rows — false-positive WCS caught by the
        # acceptance gate — are failures too and get the same treatment.
        rows = con.execute("""
            SELECT obs_rowid, stratum_id, path, readoutm
            FROM s1_solve_experiment
            WHERE status IN ('unsolved', 'timeout', 'bad_solve')
            ORDER BY stratum_id, sample_order""").fetchall()
        con.execute("DELETE FROM s1_failure_autopsy")
        for rowid, sid, relpath, readoutm in rows:
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
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    rowid, sid, m["n_sources"], m["n_psf_sources"],
                    m["median_elongation"], m["median_a_px"],
                    m["bright_median_a_px"], m["saturated_fraction"],
                    m["bkg_rms"], sat,
                    astrom.diagnose_failure(m["n_sources"],
                                            m["n_psf_sources"],
                                            m["median_elongation"],
                                            m["saturated_fraction"],
                                            m["bright_median_a_px"]),
                    str(thumb)))
            except Exception as exc:    # noqa: BLE001 — recorded, not hidden
                con.execute("""INSERT OR REPLACE INTO s1_failure_autopsy
                    (obs_rowid, stratum_id, diagnosis, thumb_path)
                    VALUES (?,?,?,NULL)""",
                            (rowid, sid,
                             f"unreadable: {type(exc).__name__}"))
            con.commit()
        n = con.execute(
            "SELECT count(*) FROM s1_failure_autopsy").fetchone()[0]
        print(f"autopsy: {n} failures examined (all of them — no cap)")
    return 0


# ---------------------------------------------------------------------------
# status / report
# ---------------------------------------------------------------------------
def cmd_status(args) -> int:
    con = sqlite3.connect(f"file:{args.manifest}?mode=ro", uri=True)
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
    d = sub.add_parser("design", help="build strata + draw samples")
    d.add_argument("--n-per-stratum", type=int,
                   default=astrom.N_PER_STRATUM)
    r = sub.add_parser("run", help="solve pending frames (resumable batch)")
    r.add_argument("--batch-seconds", type=int, default=420)
    sub.add_parser("autopsy", help="post-mortem failures")
    sub.add_parser("status", help="progress table")
    sub.add_parser("report", help="render the evidence report")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    return {"design": cmd_design, "run": cmd_run, "autopsy": cmd_autopsy,
            "status": cmd_status, "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
