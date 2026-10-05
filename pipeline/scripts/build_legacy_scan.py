#!/usr/bin/env python
"""L0 — header scan of the pre-MACRO legacy archive into ``legacy.sqlite``.

WHAT THIS SCRIPT DOES (stage L0 of the legacy census)
-----------------------------------------------------
1. **Walks** ``legacy-archive/`` once and records EVERY file it finds in the
   ``disk_files`` table (path, kind), so that "what is on disk" is a
   table, not a memory of a ``find`` command.
2. **Scans the header** of every ``*.fz`` file in parallel (headers only —
   no pixel is decompressed) into the ``scan`` table: one row per file, with
   the cards listed in ``macro_legacy.scan`` (camera identity, readout
   configuration, temperatures, optics, every time card, filter, exposure,
   target, image type).  Unreadable files get a row too, with ``error`` set:
   files-on-disk = rows, always.
3. Records each distinct header keyword set once in ``keysets`` (so "does
   any header carry card X?" is a query) and the run's provenance in
   ``scan_meta``.

It derives nothing: reconciliation, camera timeline, target census and the
go/no-go gates are ``build_legacy_census.py``, which reads these tables.

WHY ONLY ``*.fz``
-----------------
A dedupe/compress job may still be running on the archive, deleting each
uncompressed ``.fts`` whose ``.fts.fz`` twin has been verified.  The ``.fz``
files are the stable population; uncompressed files are recorded in
``disk_files`` (they are named exclusions of the census if they have no
``.fz`` twin) but never opened.  Files that vanish between the walk and the
scan are tolerated and recorded.

WHY THE SCAN CAN START BEFORE THE WALK FINISHES
-----------------------------------------------
On this spinning disk a full directory walk takes the better part of an
hour while the dedupe job is reading the same tree.  The transfer manifest
(``legacy_manifest.csv``, one row per file the Drive crawl found) already
names almost every file, so by default the scan is SEEDED from it: each
manifest path + ``.fz`` is tried at once, while the walk runs in a
background thread.  A seeded path that does not exist produces no row (it is
not on disk).  When the walk completes, ``disk_files`` is written and any
``.fz`` on disk that the seed did not cover is scanned.  The end state is
identical to walk-then-scan; ``--no-seed`` forces that order.

IDEMPOTENCE / SAFETY
--------------------
The archive is opened read-only and never written.  The scan is resumable:
a path already in ``scan`` with no error is skipped, so
an interrupted run continues where it stopped.  ``--rewalk`` refreshes
``disk_files`` (do this once the dedupe job has finished).

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_legacy_scan.py            # walk (if needed) + scan
    … build_legacy_scan.py --rewalk                      # refresh the file list
    … build_legacy_scan.py --limit 500                   # smoke test
"""

from __future__ import annotations

import argparse
import csv
import os
import sqlite3
import subprocess
import sys
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

# Make the pipeline package importable no matter where the script is run from.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_legacy import LEGACY_SCAN_VERSION                 # noqa: E402
from macro_legacy import scan as lscan                       # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/legacy-archive")
DEFAULT_DB = REPO_ROOT / "products" / "legacy" / "legacy.sqlite"
#: The transfer manifest written by the Drive crawl (file_id, rel_path, name).
DEFAULT_MANIFEST = DEFAULT_ARCHIVE.parent / "legacy_manifest.csv"

#: Uncompressed FITS suffixes (recorded, never opened).
FITS_SUFFIXES = (".fts", ".fit", ".fits")
#: Rows are committed in batches of this size (a crash loses at most one).
BATCH = 2000


# ---------------------------------------------------------------------------
# Step 1 — walk the tree
# ---------------------------------------------------------------------------
def file_kind(name: str) -> str:
    """Classify a file name: ``fz`` (scanned), ``fits`` (uncompressed), ``other``."""
    low = name.lower()
    if low.endswith(".fz"):
        return "fz"
    if low.endswith(FITS_SUFFIXES):
        return "fits"
    return "other"


def walk_archive(root: Path) -> list[tuple[str, str]]:
    """Return ``(relative path, kind)`` for every file under ``root``.

    Directory listings only — deliberately NO ``stat`` per file.  The
    archive sits on a spinning disk where a per-file stat of ~400k entries
    costs the better part of an hour (and far more while the dedupe job is
    walking the same tree); names alone are enough to classify a file, and
    the size of each scanned file is recorded by the scan worker, which has
    to touch the file anyway.  A directory that vanishes mid-walk is skipped
    by ``os.walk`` itself.
    """
    out = []
    for dirpath, _dirs, files in os.walk(root):
        rel_dir = os.path.relpath(dirpath, root)
        for name in files:
            rel = name if rel_dir == "." else os.path.join(rel_dir, name)
            out.append((rel, file_kind(name)))
    return out


def write_disk_files(con: sqlite3.Connection, rows) -> None:
    """Replace ``disk_files`` atomically with a fresh walk."""
    with con:
        con.execute("DROP TABLE IF EXISTS disk_files")
        con.execute("CREATE TABLE disk_files (path TEXT PRIMARY KEY, "
                    "kind TEXT)")
        con.executemany("INSERT INTO disk_files VALUES (?,?)", rows)


# ---------------------------------------------------------------------------
# Step 2 — scan headers
# ---------------------------------------------------------------------------
#: ``keyset`` (the full keyword list) is stored once per signature in the
#: ``keysets`` table, not on each of ~200k rows.
SCAN_DB_COLUMNS = [c for c in lscan.SCAN_COLUMNS if c != "keyset"]


def ensure_scan_tables(con: sqlite3.Connection) -> None:
    """Create ``scan``, ``keysets`` and ``scan_meta`` if absent."""
    cols = ", ".join(f'"{c}"' for c in SCAN_DB_COLUMNS)
    con.execute(f"CREATE TABLE IF NOT EXISTS scan ({cols}, PRIMARY KEY(path))")
    con.execute("CREATE TABLE IF NOT EXISTS keysets "
                "(keyset_sig TEXT PRIMARY KEY, keywords TEXT)")
    con.execute("CREATE TABLE IF NOT EXISTS scan_meta "
                "(key TEXT PRIMARY KEY, value TEXT)")


def todo_paths(con: sqlite3.Connection) -> list[str]:
    """``.fz`` paths on disk that have no clean (error-free) row yet."""
    return [r[0] for r in con.execute("""
        SELECT d.path FROM disk_files d
        LEFT JOIN scan s ON s.path = d.path
        WHERE d.kind = 'fz'
          AND (s.path IS NULL OR s.error IS NOT NULL)
        ORDER BY d.path""")]


def seed_paths(con: sqlite3.Connection, manifest_csv: Path) -> list[str]:
    """Manifest paths + ``.fz`` that have no clean row yet (see module doc).

    The manifest lists the files as the Drive holds them (uncompressed
    names); the archive holds each as ``<name>.fz`` after fpack.  A manifest
    entry that is already ``.fz`` is used as is.
    """
    have = {r[0] for r in con.execute(
        "SELECT path FROM scan WHERE error IS NULL")}
    out = []
    with open(manifest_csv, newline="") as fh:
        for rec in csv.DictReader(fh):
            rel = rec["rel_path"]
            rel = rel if rel.lower().endswith(".fz") else rel + ".fz"
            if rel not in have:
                out.append(rel)
    return sorted(set(out))


def run_scan(con: sqlite3.Connection, root: Path, paths: list[str],
             workers: int, skip_missing: bool = False) -> int:
    """Scan ``paths`` in parallel, committing in batches.  Returns row count.

    ``skip_missing`` is for the manifest-seeded pass, where a path may
    simply not exist on disk: such a path gets NO row (it is not a file on
    disk; the reconciliation names it).  In the walk-driven pass a file that
    vanished after the walk keeps its error row.
    """
    ins = "INSERT OR REPLACE INTO scan (%s) VALUES (%s)" % (
        ", ".join(f'"{c}"' for c in SCAN_DB_COLUMNS),
        ",".join("?" * len(SCAN_DB_COLUMNS)))
    batch, keysets, n, n_seen, t0 = [], {}, 0, 0, time.time()

    def flush():
        with con:
            con.executemany(ins, batch)
            con.executemany("INSERT OR IGNORE INTO keysets VALUES (?,?)",
                            list(keysets.items()))
        batch.clear()
        keysets.clear()

    jobs = [(str(root), p) for p in paths]
    with ProcessPoolExecutor(workers) as pool:
        for row in pool.map(lscan.scan_one, jobs, chunksize=32):
            n_seen += 1
            if n_seen % 10000 == 0:
                rate = n_seen / (time.time() - t0)
                print(f"  scanned {n_seen:,}/{len(paths):,}  "
                      f"({rate:,.0f} files/s)", flush=True)
            if skip_missing and (row.get("error") or "").startswith(
                    "FileNotFoundError"):
                continue
            if row.get("keyset_sig"):
                keysets[row["keyset_sig"]] = row["keyset"]
            batch.append([row.get(c) for c in SCAN_DB_COLUMNS])
            n += 1
            if len(batch) >= BATCH:
                flush()
    if batch:
        flush()
    return n


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
def git_commit() -> str:
    """Short commit of the repo, with '+dirty' when the tree has changes."""
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO_ROOT, capture_output=True, text=True,
                             check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT,
                               capture_output=True, text=True).stdout.strip()
        return sha + ("+dirty" if dirty else "")
    except Exception:                      # noqa: BLE001 — provenance only
        return ""


def write_meta(con: sqlite3.Connection, **items) -> None:
    with con:
        con.executemany("INSERT OR REPLACE INTO scan_meta VALUES (?,?)",
                        [(k, str(v)) for k, v in items.items()])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--rewalk", action="store_true",
                    help="refresh disk_files even if it already exists")
    ap.add_argument("--walk-only", action="store_true",
                    help="refresh disk_files and stop (no header is opened)")
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                    help="transfer manifest used to seed the scan")
    ap.add_argument("--no-seed", action="store_true",
                    help="walk first, then scan (no manifest seeding)")
    ap.add_argument("--limit", type=int, default=None,
                    help="scan at most this many files (smoke test)")
    args = ap.parse_args(argv)

    if not args.archive.is_dir():
        print(f"archive not found: {args.archive}", file=sys.stderr)
        return 1
    args.db.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(args.db, timeout=120)
    # Write-ahead logging: appends sequentially instead of journalling every
    # replaced page — on the spinning disk this database lives on, replacing
    # the derived tables is otherwise tens of minutes.
    con.execute("PRAGMA journal_mode=WAL")
    try:
        ensure_scan_tables(con)
        have_walk = con.execute(
            "SELECT count(*) FROM sqlite_master WHERE name='disk_files'"
        ).fetchone()[0]
        need_walk = args.rewalk or args.walk_only or not have_walk

        # The walk runs in a background thread so the seeded scan can use
        # the same wall-clock time (directory listing is I/O-bound; the GIL
        # is released while it waits on the disk).
        walked: dict = {}
        walker = None
        if need_walk:
            print(f"walking {args.archive} …", flush=True)
            walker = threading.Thread(
                target=lambda: walked.update(rows=walk_archive(args.archive)))
            walker.start()

        n = 0
        seed = (not args.no_seed and not args.walk_only
                and args.manifest.exists())
        if seed:
            paths = seed_paths(con, args.manifest)
            if args.limit:
                paths = paths[:args.limit]
            print(f"seeded from {args.manifest.name}: {len(paths):,} headers "
                  "to try", flush=True)
            n += run_scan(con, args.archive, paths, args.workers,
                          skip_missing=True) if paths else 0

        if walker is not None:
            walker.join()
            write_disk_files(con, walked["rows"])
            write_meta(con, walk_utc=datetime.now(timezone.utc).isoformat(
                timespec="seconds"), walk_n_files=len(walked["rows"]))
            print(f"  walk done: {len(walked['rows']):,} files on disk",
                  flush=True)
        if args.walk_only:
            return 0

        paths = todo_paths(con)
        if args.limit:
            paths = paths[:args.limit]
        print(f"headers still to scan after the walk: {len(paths):,}",
              flush=True)
        n += run_scan(con, args.archive, paths, args.workers) if paths else 0
        # A row whose file is no longer on disk (it can only have come from
        # an earlier walk) is not a file on disk: drop it so rows = files.
        with con:
            gone = con.execute(
                "DELETE FROM scan WHERE path NOT IN "
                "(SELECT path FROM disk_files WHERE kind='fz')").rowcount
        if gone:
            print(f"  dropped {gone:,} rows whose file is no longer on disk")
        total, errs = con.execute(
            "SELECT count(*), sum(error IS NOT NULL) FROM scan").fetchone()
        write_meta(con,
                   scan_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   scan_code_version=LEGACY_SCAN_VERSION,
                   scan_git_commit=git_commit(),
                   archive_root=str(args.archive))
        print(f"DONE: scanned {n:,} this run; scan table holds {total:,} rows, "
              f"{errs or 0:,} with errors")
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
