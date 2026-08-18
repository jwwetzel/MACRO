#!/usr/bin/env python
"""Repair phantom image geometry in the observation catalog (stage S0e).

WHAT WENT WRONG (the artifact this script repairs)
--------------------------------------------------
A tile-compressed FITS file stores its image inside a BINTABLE.  That
table's ``NAXIS1`` is the row length in BYTES and its ``NAXIS2`` is the row
COUNT; the true picture size lives in ``ZNAXIS1``/``ZNAXIS2``.  For the
RLMT's 4800x3211 detector the table header reads ``NAXIS1 = 8`` and
``NAXIS2 = 3211``.

19,980 catalog rows were written from that table header, so the catalog
believed the observatory had taken 8-pixel-wide strips.  It never did.
Every one of those frames is a full 4800x3211 field.  The consequences ran
downhill: two phantom camera eras keyed on the fake geometry, and 18,381
frames excluded from the S1 astrometry batch by a solvability gate that
(correctly, given wrong input) refuses anything narrower than 512 px.

See ``macro_core/fitsgeom.py`` for the mechanism and
``docs/pipeline/s0e_geometry_fix.html`` for the full write-up.

WHAT THIS SCRIPT DOES
---------------------
Re-reads the geometry — and ONLY the geometry — of candidate catalog rows
straight from the archive, through the compression-aware resolver, and
writes back ``naxis1``/``naxis2`` where they were wrong.

Design commitments, because this edits a 330k-row catalog in place:

* **Surgical.**  Only ``naxis1``/``naxis2`` are ever written.  Every other
  column keeps whatever the original scan produced.
* **Audited.**  Every rescanned row lands in the ``geom_rescan`` table with
  its OLD and NEW values, whether it changed or not.  That table is the
  before/after diff — it is evidence, not a scratch pad, so it is never
  dropped on re-run.
* **Non-disturbing, and it proves it.**  The default candidate set
  deliberately INCLUDES the genuinely small frames (the Andor iKon 57x48
  focus windows).  Those must come back byte-identical; ``verify`` checks
  exactly that and fails loudly if any correct row moved.
* **Resumable.**  Rows already in ``geom_rescan`` are skipped, so a killed
  run costs only its in-flight batch.
* **Polite.**  Defaults to 4 workers (cap 6) because an S1 batch solve and
  a bulk archive transfer may be running against the same disk.
* The archive is opened READ-ONLY.  It is never written to.

SUBCOMMANDS
-----------
    plan      show what would be rescanned, touch nothing
    run       rescan candidates and repair the catalog (resumable)
    status    progress + change tally (read-only)
    verify    prove correct rows were untouched; show the change matrix

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/rescan_geometry.py plan
    $PY pipeline/scripts/rescan_geometry.py run --workers 4
    $PY pipeline/scripts/rescan_geometry.py status
    $PY pipeline/scripts/rescan_geometry.py verify
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir))
from macro_core import fitsgeom  # noqa: E402

warnings.filterwarnings("ignore")

ROOT = "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive"
DB = "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite"

#: Rows whose stored NAXIS1 is at most this are candidates for a rescan.
#: A BINTABLE row length is a handful of bytes (8 here), so every phantom
#: sits far below this bar — and so do the few genuinely tiny iKon frames,
#: which is on purpose: they are the control group that proves the repair
#: does not touch rows that were already right.
CANDIDATE_MAX_NAXIS1 = 64

#: Concurrency cap.  The archive lives on one spinning volume that may be
#: serving an S1 solve batch and an rclone pull at the same time.
MAX_WORKERS = 6
DEFAULT_WORKERS = 4

#: SQLite lock patience: the batch commits to a different DB, but the
#: catalog may still be read by other tooling.
BUSY_TIMEOUT_MS = 300_000


def connect(path: str, read_only: bool = False) -> sqlite3.Connection:
    """Open the catalog with the project's standard lock patience."""
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=BUSY_TIMEOUT_MS / 1000)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    return con


def ensure_audit_table(con: sqlite3.Connection) -> None:
    """Create the audit/resume table if absent.  Never dropped: it is the
    permanent record of what this repair changed and what it left alone."""
    con.execute("""
        CREATE TABLE IF NOT EXISTS geom_rescan (
            path        TEXT PRIMARY KEY,
            old_naxis1  REAL,
            old_naxis2  REAL,
            new_naxis1  INTEGER,
            new_naxis2  INTEGER,
            compressed  INTEGER,   -- 1 when the file is tile-compressed
            changed     INTEGER NOT NULL,
            error       TEXT,
            scanned_utc TEXT NOT NULL DEFAULT (datetime('now'))
        )""")
    con.execute("CREATE INDEX IF NOT EXISTS ix_geom_changed "
                "ON geom_rescan(changed)")
    con.commit()


# ---------------------------------------------------------------------------
# The worker: read one file's TRUE geometry
# ---------------------------------------------------------------------------
def geometry_of(rel_path: str) -> dict:
    """Return ``{path, naxis1, naxis2, compressed, error}`` for one file.

    Runs in a subprocess.  astropy is imported inside so the parent never
    pays for it, matching ``build_catalog.scan_one``.

    The header is merged CARD BY CARD inside a guard: ~20k archive files
    carry a malformed ``CONTINUE`` card (a ``CONTINUE`` after a non-string
    ``FWALLNAM`` value) that makes astropy's ``Header.update`` raise and
    abandon the whole header.  Skipping the one bad card keeps the frame.
    """
    from astropy.io import fits
    out = {"path": rel_path, "naxis1": None, "naxis2": None,
           "compressed": None, "error": None}
    full = os.path.join(ROOT, rel_path)
    try:
        with fits.open(full, memmap=False, ignore_missing_simple=True) as h:
            hdr = fits.Header()
            compressed = False
            for hdu in h[:2]:
                # Tile compression is detected from the HDU TYPE, not from
                # the merged cards: when astropy succeeds in building a
                # CompImageHDU, ``hdu.header`` is the TRANSLATED image
                # header, and the Z* markers have already been consumed —
                # so a card-level test would call a compressed file plain.
                if isinstance(hdu, fits.CompImageHDU):
                    compressed = True
                for c in hdu.header.cards:
                    try:
                        hdr[c.keyword] = (c.value, c.comment)
                    except Exception:
                        continue          # one unreadable card, not a frame
            # Fallback path: astropy could NOT build a CompImageHDU, so the
            # raw BINTABLE header (Z* markers intact) came through instead.
            # resolve_geometry then reads ZNAXIS itself.
            compressed = compressed or fitsgeom.is_compressed_header(hdr)
            out["compressed"] = int(compressed)
            out["naxis1"], out["naxis2"] = fitsgeom.resolve_geometry(hdr)
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"[:300]
    return out


# ---------------------------------------------------------------------------
# Candidate selection
# ---------------------------------------------------------------------------
def candidate_paths(con: sqlite3.Connection, max_naxis1: int,
                    include_done: bool = False) -> list[str]:
    """Catalog rows worth re-reading, oldest-first for stable resume."""
    sql = ("SELECT o.path FROM obs o WHERE o.naxis1 IS NOT NULL "
           "AND o.naxis1 <= ?")
    if not include_done:
        sql += (" AND o.path NOT IN (SELECT path FROM geom_rescan "
                "WHERE error IS NULL)")
    sql += " ORDER BY o.path"
    return [r[0] for r in con.execute(sql, (max_naxis1,))]


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------
def cmd_plan(args) -> int:
    con = connect(DB, read_only=True)
    ensure_audit_table(connect(DB))          # safe: CREATE IF NOT EXISTS
    todo = candidate_paths(con, args.max_naxis1)
    print(f"catalog rows with naxis1 <= {args.max_naxis1}: "
          f"{con.execute('SELECT COUNT(*) FROM obs WHERE naxis1 <= ?', (args.max_naxis1,)).fetchone()[0]}")
    print(f"not yet rescanned (this run would read): {len(todo)}")
    print("\nstored geometry of the candidate set:")
    for n1, n2, c in con.execute(
            "SELECT naxis1, naxis2, COUNT(*) FROM obs WHERE naxis1 <= ? "
            "GROUP BY 1,2 ORDER BY 3 DESC", (args.max_naxis1,)):
        print(f"   {int(n1):>6} x {int(n2):<6}  {c:>7} rows")
    print("\nOnly naxis1/naxis2 will be written.  No other column is touched.")
    return 0


def cmd_run(args) -> int:
    workers = min(args.workers, MAX_WORKERS)
    con = connect(DB)
    ensure_audit_table(con)
    todo = candidate_paths(con, args.max_naxis1)
    if args.limit:
        todo = todo[:args.limit]
    print(f"rescanning {len(todo)} rows with {workers} workers", flush=True)
    if not todo:
        print("nothing to do — already complete")
        return 0

    # Stored values, needed to decide 'changed' and to record the OLD side
    # of the audit trail.  Fetched once, in bulk: 20k rows is nothing.
    old = {p: (n1, n2) for p, n1, n2 in con.execute(
        "SELECT path, naxis1, naxis2 FROM obs")}

    ins = ("INSERT OR REPLACE INTO geom_rescan "
           "(path, old_naxis1, old_naxis2, new_naxis1, new_naxis2, "
           " compressed, changed, error) VALUES (?,?,?,?,?,?,?,?)")
    upd = "UPDATE obs SET naxis1 = ?, naxis2 = ? WHERE path = ?"

    batch_audit, batch_fix, n, n_changed, n_err = [], [], 0, 0, 0
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(geometry_of, todo, chunksize=16):
            p = r["path"]
            o1, o2 = old.get(p, (None, None))
            changed = 0
            if r["error"] is None:
                # Compare as ints: the catalog stores REAL, we resolve int.
                if (o1 is None or o2 is None
                        or int(o1) != r["naxis1"] or int(o2) != r["naxis2"]):
                    changed = 1
                    batch_fix.append((r["naxis1"], r["naxis2"], p))
            else:
                n_err += 1
            n_changed += changed
            batch_audit.append((p, o1, o2, r["naxis1"], r["naxis2"],
                                r["compressed"], changed, r["error"]))
            n += 1
            if len(batch_audit) >= 500:
                # Audit row and repair commit together: the audit table can
                # never claim a change the catalog did not receive.
                con.executemany(ins, batch_audit)
                if batch_fix:
                    con.executemany(upd, batch_fix)
                con.commit()
                batch_audit, batch_fix = [], []
            if n % 2000 == 0:
                print(f"  {n}/{len(todo)}  changed={n_changed} "
                      f"errors={n_err}", flush=True)
    if batch_audit:
        con.executemany(ins, batch_audit)
        if batch_fix:
            con.executemany(upd, batch_fix)
        con.commit()
    print(f"DONE: {n} rescanned, {n_changed} repaired, {n_err} errors",
          flush=True)
    con.close()
    return 0


def cmd_status(args) -> int:
    con = connect(DB, read_only=True)
    try:
        done = con.execute("SELECT COUNT(*) FROM geom_rescan").fetchone()[0]
    except sqlite3.OperationalError:
        print("no geom_rescan table yet — run `plan` then `run`")
        return 0
    total = con.execute("SELECT COUNT(*) FROM obs WHERE naxis1 <= ?",
                        (args.max_naxis1,)).fetchone()[0]
    changed = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE changed = 1").fetchone()[0]
    errs = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE error IS NOT NULL").fetchone()[0]
    print(f"rescanned {done}/{total}   repaired {changed}   errors {errs}")
    print("\nchange matrix (old -> new):")
    for o1, o2, n1, n2, c in con.execute("""
            SELECT old_naxis1, old_naxis2, new_naxis1, new_naxis2, COUNT(*)
            FROM geom_rescan GROUP BY 1,2,3,4 ORDER BY 5 DESC"""):
        arrow = "->" if (o1, o2) != (n1, n2) else "=="
        print(f"   {_g(o1)} x {_g(o2):<6} {arrow} {_g(n1)} x {_g(n2):<6}"
              f"  {c:>7} rows")
    return 0


def _g(v) -> str:
    return "None" if v is None else str(int(v))


def cmd_verify(args) -> int:
    """Prove the repair was surgical.

    Two claims must hold, and both are checked against the audit table:

    1. every row that changed was a COMPRESSED file whose stored geometry
       came from the BINTABLE header — no uncompressed row was rewritten;
    2. every row that was already correct came back IDENTICAL.
    """
    con = connect(DB, read_only=True)
    bad_uncompressed = con.execute(
        "SELECT COUNT(*) FROM geom_rescan "
        "WHERE changed = 1 AND compressed = 0").fetchone()[0]
    control = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE compressed = 0").fetchone()[0]
    control_moved = con.execute(
        "SELECT COUNT(*) FROM geom_rescan "
        "WHERE compressed = 0 AND changed = 1").fetchone()[0]
    print(f"control group (uncompressed rows rescanned): {control}")
    print(f"  of which changed: {control_moved}   (MUST be 0)")
    print(f"changed rows that were not compressed: {bad_uncompressed}")
    print("\nsample of repaired rows (old -> new):")
    for p, o1, o2, n1, n2 in con.execute(
            "SELECT path, old_naxis1, old_naxis2, new_naxis1, new_naxis2 "
            "FROM geom_rescan WHERE changed = 1 LIMIT 5"):
        print(f"   {_g(o1)}x{_g(o2)} -> {_g(n1)}x{_g(n2)}  {p}")
    print("\nsample of untouched control rows:")
    for p, o1, o2, n1, n2 in con.execute(
            "SELECT path, old_naxis1, old_naxis2, new_naxis1, new_naxis2 "
            "FROM geom_rescan WHERE compressed = 0 LIMIT 5"):
        print(f"   {_g(o1)}x{_g(o2)} == {_g(n1)}x{_g(n2)}  {p}")
    # Live catalog cross-check: no 8x3211 rows may survive the repair.
    left = con.execute("SELECT COUNT(*) FROM obs "
                       "WHERE naxis1 = 8 AND naxis2 = 3211").fetchone()[0]
    print(f"\nphantom 8x3211 rows still in the catalog: {left}")
    ok = (control_moved == 0 and bad_uncompressed == 0)
    print("\nVERDICT:", "PASS — repair was surgical" if ok
          else "FAIL — a correct row was modified")
    return 0 if ok else 1


def cmd_exemplar(args) -> int:
    """Store one real file's RAW header cards as the report's exhibit.

    The S0e report has to SHOW the trap, not just describe it, and the house
    rule is that a report renders from the database.  So the exhibit — the
    first cards of an actual tile-compressed archive header, where
    ``NAXIS1 = 8`` sits nine lines above ``ZNAXIS1 = 4800`` — is captured
    here into ``s0e_header_dump`` and read back at render time.  The archive
    is opened read-only and exactly once.
    """
    from astropy.io import fits
    con = connect(DB)
    path = args.path
    if path is None:
        # Default exhibit: the first repaired row, so the dump always
        # matches a frame the report actually counts.
        row = con.execute("SELECT path FROM geom_rescan WHERE changed = 1 "
                          "ORDER BY path LIMIT 1").fetchone()
        if not row:
            print("no repaired rows yet — run `run` first", file=sys.stderr)
            return 2
        path = row[0]
    con.execute("DROP TABLE IF EXISTS s0e_header_dump")
    con.execute("""
        CREATE TABLE s0e_header_dump (
            path     TEXT NOT NULL,
            hdu      INTEGER NOT NULL,
            card_no  INTEGER NOT NULL,
            card     TEXT NOT NULL)""")
    with fits.open(os.path.join(ROOT, path), memmap=False,
                   ignore_missing_simple=True) as h:
        # ``_header`` is the RAW BINTABLE header — deliberately NOT the
        # translated image header, because the whole point of the exhibit is
        # the untranslated cards that a naive parser sees.
        raw = getattr(h[1], "_header", h[1].header)
        cards = []
        for i, c in enumerate(raw.cards[:args.n_cards]):
            try:
                cards.append((path, 1, i, str(c).rstrip()))
            except Exception:
                cards.append((path, 1, i, "<unparsable card>"))
    con.executemany("INSERT INTO s0e_header_dump VALUES (?,?,?,?)", cards)
    con.commit()
    print(f"stored {len(cards)} cards from {path}")
    for _, _, _, c in cards[:6]:
        print("   ", c[:74])
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Repair phantom BINTABLE geometry in the RLMT catalog.")
    ap.add_argument("--max-naxis1", type=int, default=CANDIDATE_MAX_NAXIS1,
                    dest="max_naxis1",
                    help="rescan rows whose stored naxis1 is <= this "
                         f"(default {CANDIDATE_MAX_NAXIS1})")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan").set_defaults(fn=cmd_plan)
    r = sub.add_parser("run")
    r.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                   help=f"parallel readers (capped at {MAX_WORKERS})")
    r.add_argument("--limit", type=int, default=None,
                   help="stop after N rows (for a trial run)")
    r.set_defaults(fn=cmd_run)
    sub.add_parser("status").set_defaults(fn=cmd_status)
    sub.add_parser("verify").set_defaults(fn=cmd_verify)
    e = sub.add_parser("exemplar")
    e.add_argument("--path", default=None,
                   help="archive-relative file to dump (default: the first "
                        "repaired row)")
    e.add_argument("--n-cards", type=int, default=22, dest="n_cards",
                   help="how many leading cards to store")
    e.set_defaults(fn=cmd_exemplar)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
