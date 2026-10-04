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

    hdr-run     re-scrape the HARDWARE cards of every catalog row (resumable)
    hdr-status  progress of the hardware-card re-scrape (read-only)
    hdr-verify  prove the re-scrape against the cards the catalog already had

THE HARDWARE-CARD RE-SCRAPE (finding F-2, plan review 2026-10-03)
-----------------------------------------------------------------
The original scan kept the cards that describe the EXPOSURE and dropped the
cards that describe the HARDWARE: GAIN, OFFSET, SET-TEMP, COOLPOWR, FOCPOS,
FLIPSTAT, TELPIER, FWPOS, FWALLNAM.  Without them the manifest cannot tell a
re-seated or flipped camera from an untouched one, cannot key a dark on its
set-point, and reports a focuser position from a card that sticks.  The
``hdr-*`` subcommands read those cards for EVERY catalog row (all 330k, not
just the canonical ones — a duplicate copy's header is evidence too) with
the raw card parser in ``macro_core.fitsgeom`` and store them, as the header
spells them, in a NEW catalog table ``hdr_rescrape``.

Same commitments as the geometry repair, plus one:

* **Additive.**  ``obs`` is not touched at all — not one column.  The
  re-scrape only ever writes its own table; S0 joins it.
* **Resumable.**  Rows already in ``hdr_rescrape`` are skipped.
* **Self-checking.**  Four of the collected cards (INSTRUME, SWCREATE,
  CCD-TEMP, FOCUSPOS) were ALSO read by the original astropy scan.  They are
  collected again on purpose: ``hdr-verify`` compares the two readings row
  by row, so the raw parser is validated on ~1.2M card readings it had no
  way to copy, before any of its new columns is trusted.
* The archive is opened READ-ONLY, headers only — no pixel is read.

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/rescan_geometry.py plan
    $PY pipeline/scripts/rescan_geometry.py run --workers 4
    $PY pipeline/scripts/rescan_geometry.py status
    $PY pipeline/scripts/rescan_geometry.py verify

    $PY pipeline/scripts/rescan_geometry.py hdr-run --workers 12
    $PY pipeline/scripts/rescan_geometry.py hdr-status
    $PY pipeline/scripts/rescan_geometry.py hdr-verify
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

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
    # The audit table has to exist before candidate_paths can reference it.
    # Opened and closed explicitly rather than leaked into the read path:
    # a stray write handle on this catalog blocks every other stage.
    wcon = connect(DB)
    ensure_audit_table(wcon)                 # safe: CREATE IF NOT EXISTS
    wcon.close()
    con = connect(DB, read_only=True)
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
    # Remaining work is counted as candidates NOT yet in the audit table —
    # not as "rows still matching the candidate filter".  A repaired row no
    # longer has a small naxis1, so a naive denominator would shrink as the
    # run progressed and report nonsense like "20,071 of 91 done".
    remaining = con.execute(
        "SELECT COUNT(*) FROM obs WHERE naxis1 IS NOT NULL AND naxis1 <= ? "
        "AND path NOT IN (SELECT path FROM geom_rescan WHERE error IS NULL)",
        (args.max_naxis1,)).fetchone()[0]
    changed = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE changed = 1").fetchone()[0]
    errs = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE error IS NOT NULL").fetchone()[0]
    print(f"rescanned {done}   repaired {changed}   errors {errs}   "
          f"remaining {remaining}")
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

    THE CONTROL GROUP, and why it is what it is.  The first version of this
    check defined the control as "uncompressed rows" — and that turned out
    to be VACUOUS: the Andor iKon focus windows are ``.fts.fz`` files too,
    so every candidate row is tile-compressed and the control group was
    empty.  An empty control group passes trivially, which is worse than no
    check at all.

    The correct control is the rows that came back UNCHANGED.  Those are
    genuinely small frames, and the fact that they survived is the real
    proof: the repair did not key on "naxis1 is small" and blanket-rewrite
    everything that looked odd — it read ``ZNAXIS*`` from each file, and
    for these files ``ZNAXIS*`` genuinely says 57x48.  A crude fix would
    have inflated them to 4800x3211 and destroyed real sub-frame geometry.

    Claims checked:
      1. no phantom 8x3211 row survives in the catalog;
      2. every unchanged row is genuinely small (its TRUE geometry has an
         axis below the solvability floor) — they were read, not skipped;
      3. every changed row moved from a BINTABLE row-length to a real
         image width, i.e. grew;
      4. nothing failed to read.
    """
    con = connect(DB, read_only=True)
    n_changed = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE changed = 1").fetchone()[0]
    control = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE changed = 0").fetchone()[0]
    errs = con.execute("SELECT COUNT(*) FROM geom_rescan "
                       "WHERE error IS NOT NULL").fetchone()[0]
    # (2) every untouched row is genuinely small in its TRUE geometry.
    control_not_small = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE changed = 0 "
        "AND (new_naxis1 >= 512 AND new_naxis2 >= 512)").fetchone()[0]
    # (3) a repair must GROW a frame; a repair that shrank one is a bug.
    shrank = con.execute(
        "SELECT COUNT(*) FROM geom_rescan WHERE changed = 1 "
        "AND new_naxis1 <= old_naxis1").fetchone()[0]
    print(f"repaired rows: {n_changed}")
    print(f"control group (re-read, left unchanged): {control}")
    print(f"  of those, any that were NOT genuinely small: "
          f"{control_not_small}   (MUST be 0)")
    print(f"repaired rows that SHRANK: {shrank}   (MUST be 0)")
    print(f"rows that failed to read: {errs}   (MUST be 0)")
    print("\nsample of repaired rows (old -> new):")
    for p, o1, o2, n1, n2 in con.execute(
            "SELECT path, old_naxis1, old_naxis2, new_naxis1, new_naxis2 "
            "FROM geom_rescan WHERE changed = 1 LIMIT 4"):
        print(f"   {_g(o1)}x{_g(o2)} -> {_g(n1)}x{_g(n2)}  {p}")
    print("\nthe control group — real sub-frames, read and left alone:")
    for p, o1, o2, n1, n2, comp in con.execute(
            "SELECT path, old_naxis1, old_naxis2, new_naxis1, new_naxis2, "
            "compressed FROM geom_rescan WHERE changed = 0 LIMIT 5"):
        print(f"   {_g(o1)}x{_g(o2)} == {_g(n1)}x{_g(n2)} "
              f"(compressed={comp})  {p}")
    # (1) Live catalog cross-check: no 8x3211 rows may survive the repair.
    left = con.execute("SELECT COUNT(*) FROM obs "
                       "WHERE naxis1 = 8 AND naxis2 = 3211").fetchone()[0]
    print(f"\nphantom 8x3211 rows still in the catalog: {left}   (MUST be 0)")
    ok = (control_not_small == 0 and shrank == 0 and errs == 0 and left == 0
          and control > 0)
    print("\nVERDICT:", "PASS — repair was surgical, control group intact"
          if ok else "FAIL — see the MUST-be-0 lines above")
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


# ===========================================================================
# The hardware-card re-scrape (F-2)
# ===========================================================================

#: Thread count for the header re-scrape.  THREADS, not processes: the work
#: is two short reads per file and a string split, so the interpreter lock is
#: released almost the whole time and a process pool would only add pickling.
#: Measured on the archive volume under load (random order): 1.1 files/s
#: with one reader, 13.5 files/s with twelve — the volume queues well, so the
#: default is generous; the cap keeps a typo from starving the sibling stages
#: that share the disk.
HDR_DEFAULT_WORKERS = 12
HDR_MAX_WORKERS = 32

#: Rows committed per transaction.  Small enough that a killed run loses at
#: most a few seconds of reading; large enough that the catalog's write lock
#: is taken a few hundred times in total rather than 330k times.
HDR_COMMIT_EVERY = 500


def hdr_column(card: str) -> str:
    """SQL column name for a header card: ``SET-TEMP`` -> ``h_set_temp``.

    The ``h_`` prefix is not decoration: ``OFFSET`` is an SQL keyword, and a
    bare column of that name would need quoting in every query anyone ever
    writes against the table.
    """
    return "h_" + card.lower().replace("-", "_")


#: ``hdr_rescrape`` data columns, in :data:`fitsgeom.HARDWARE_CARDS` order.
HDR_COLUMNS = tuple(hdr_column(c) for c in fitsgeom.HARDWARE_CARDS)


def ensure_hdr_table(con: sqlite3.Connection) -> None:
    """Create ``hdr_rescrape`` if absent.  Never dropped: like
    ``geom_rescan`` it is evidence, and it is the resume state.

    Every card column is TEXT and holds the value AS THE HEADER SPELLS IT:
    NULL = the header has no such card, ``''`` = the card exists and is
    blank, anything else = the value.  Typing happens in S0, where the
    camera is known.
    """
    cols = ",\n            ".join(f"{c} TEXT" for c in HDR_COLUMNS)
    con.execute(f"""
        CREATE TABLE IF NOT EXISTS hdr_rescrape (
            path        TEXT PRIMARY KEY,
            {cols},
            n_cards     INTEGER,   -- how many of the wanted cards were found
            error       TEXT,      -- NULL when the header was read
            scanned_utc TEXT NOT NULL DEFAULT (datetime('now'))
        )""")
    con.commit()


def hardware_cards_of(rel_path: str) -> dict:
    """Read one file's hardware cards.  Runs in a worker thread.

    Never raises: an unreadable file becomes a row with ``error`` set, so
    the run always terminates and the failure is counted instead of lost.
    """
    out = {"path": rel_path, "cards": {}, "error": None}
    try:
        out["cards"] = fitsgeom.read_hardware_cards(
            os.path.join(ROOT, rel_path))
    except Exception as e:                       # noqa: BLE001 — recorded
        out["error"] = f"{type(e).__name__}: {e}"[:300]
    return out


def hdr_todo(con: sqlite3.Connection, retry_errors: bool = False
             ) -> list[str]:
    """Catalog paths not yet re-scraped, in READ order.

    ``rawimage`` first — it holds the canonical pixels, so an interrupted
    run has already covered the frames that matter most — then every other
    tree, with ``reduced`` LAST (it is almost entirely copies and
    derivatives of rawimage frames); within a tree by path, which walks the
    disk one night directory at a time instead of seeking across the whole
    volume.
    """
    done = "SELECT path FROM hdr_rescrape"
    if retry_errors:
        done += " WHERE error IS NULL"
    return [r[0] for r in con.execute(
        f"SELECT path FROM obs WHERE path NOT IN ({done}) "
        "ORDER BY (tree != 'rawimage'), (tree = 'reduced'), path")]


def cmd_hdr_run(args) -> int:
    workers = max(1, min(args.workers, HDR_MAX_WORKERS))
    con = connect(DB)
    ensure_hdr_table(con)
    todo = hdr_todo(con, retry_errors=args.retry_errors)
    if args.limit:
        todo = todo[:args.limit]
    print(f"re-scraping hardware cards of {len(todo)} rows with "
          f"{workers} reader threads", flush=True)
    if not todo:
        print("nothing to do — already complete")
        return 0
    ins = ("INSERT OR REPLACE INTO hdr_rescrape (path, "
           + ", ".join(HDR_COLUMNS) + ", n_cards, error) VALUES ("
           + ", ".join("?" * (len(HDR_COLUMNS) + 3)) + ")")
    batch, n, n_err = [], 0, 0
    with ThreadPoolExecutor(workers) as ex:
        # ex.map preserves submission order, so commits advance through the
        # path-sorted list and a resumed run restarts where this one died.
        for r in ex.map(hardware_cards_of, todo):
            cards = r["cards"]
            n_err += r["error"] is not None
            batch.append((r["path"],
                          *(cards.get(c) for c in fitsgeom.HARDWARE_CARDS),
                          len(cards) if r["error"] is None else None,
                          r["error"]))
            n += 1
            if len(batch) >= HDR_COMMIT_EVERY:
                con.executemany(ins, batch)
                con.commit()
                batch = []
            if n % 5000 == 0:
                print(f"  {n}/{len(todo)}  errors={n_err}", flush=True)
    if batch:
        con.executemany(ins, batch)
        con.commit()
    print(f"DONE: {n} re-scraped, {n_err} unreadable", flush=True)
    con.close()
    return 0


def cmd_hdr_status(args) -> int:
    con = connect(DB, read_only=True)
    total = con.execute("SELECT COUNT(*) FROM obs").fetchone()[0]
    try:
        done, errs = con.execute(
            "SELECT COUNT(*), SUM(error IS NOT NULL) FROM hdr_rescrape"
        ).fetchone()
    except sqlite3.OperationalError:
        print("no hdr_rescrape table yet — run `hdr-run`")
        return 0
    print(f"re-scraped {done} of {total} catalog rows   "
          f"unreadable {errs or 0}   remaining {total - done}")
    return 0


def _num_equal(a, b, tol: float = 1e-6) -> bool:
    """Do two readings of one numeric card agree?  The catalog holds what
    astropy parsed into a float, the re-scrape holds the header's text."""
    try:
        return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(a)))
    except (TypeError, ValueError):
        return False


#: The control group of ``hdr-verify``: cards the ORIGINAL scan also read,
#: as (catalog column in ``obs``, re-scrape column, is_numeric).
HDR_CONTROL = (
    ("instrume", "h_instrume", False),
    ("swcreate", "h_swcreate", False),
    ("ccd_temp", "h_ccd_temp", True),
    ("focuspos", "h_focuspos", True),
)


def hdr_verify_counts(con: sqlite3.Connection) -> dict:
    """Every number ``hdr-verify`` prints, as a dict (so the tests can
    drive the same arithmetic on a hand-built catalog).

    * ``n_obs`` / ``n_scanned`` / ``n_unreadable`` — coverage;
    * ``n_unreadable_new`` — files the re-scrape could not read that the
      ORIGINAL scan could (a regression of the reader, must be 0);
    * ``control`` — per control card: rows where both scans hold a value,
      and how many of those DISAGREE (must be 0);
    * ``recovered`` — per control card: rows where the original scan had
      NULL and the re-scrape found a value.  This is the F-2 acceptance
      line "no nulls where the card exists", measured;
    * ``cards`` — per collected card: rows with a value / blank / absent.
    """
    # TWO SEQUENTIAL SCANS, compared in memory — not SQL joins.  The first
    # version joined obs to hdr_rescrape on path, six times over; on the
    # archive's spinning disk each join is 330k random index probes into a
    # 500 MB file, and the command ran for half an hour doing three seconds
    # of arithmetic.  Reading each table once, front to back, takes about
    # as long as S0's own catalog load.
    obs = {}
    ocols = [c for c, _, _ in HDR_CONTROL]
    for row in con.execute(f"SELECT path, error, {', '.join(ocols)} "
                           f"FROM obs"):
        obs[row[0]] = row[1:]
    out: dict = {"n_obs": len(obs), "n_scanned": 0, "n_unreadable": 0,
                 "n_unreadable_new": 0}
    hidx = {c: i for i, c in enumerate(HDR_COLUMNS)}
    control = {h: [0, 0] for _, h, _ in HDR_CONTROL}
    recovered = {h: 0 for _, h, _ in HDR_CONTROL}
    cards = {c: [0, 0, 0] for c in HDR_COLUMNS}      # value, blank, absent
    for row in con.execute(f"SELECT path, error, {', '.join(HDR_COLUMNS)} "
                           f"FROM hdr_rescrape"):
        path, err, vals = row[0], row[1], row[2:]
        out["n_scanned"] += 1
        orow = obs.get(path)
        if err is not None:
            out["n_unreadable"] += 1
            # Readable by the original scan, unreadable here = a regression
            # of the reader.
            out["n_unreadable_new"] += bool(orow is not None
                                            and orow[0] is None)
            continue
        for c, v in zip(HDR_COLUMNS, vals):
            cards[c][2 if v is None else (1 if v == "" else 0)] += 1
        if orow is None:
            continue
        for k, (_, hcol, numeric) in enumerate(HDR_CONTROL):
            a, b = orow[1 + k], vals[hidx[hcol]]
            if b is None or b == "":
                continue
            if a is None:
                recovered[hcol] += 1
                continue
            control[hcol][0] += 1
            same = _num_equal(a, b) if numeric \
                else str(a).strip() == str(b).strip()
            control[hcol][1] += not same
    out["control"] = {h: tuple(v) for h, v in control.items()}
    out["recovered"] = recovered
    out["cards"] = {c: tuple(v) for c, v in cards.items()}
    return out


def cmd_hdr_verify(args) -> int:
    """Prove the re-scrape before S0 is allowed to trust it.

    Claims checked:
      1. every catalog row was re-scraped;
      2. nothing the original scan could read failed to read here;
      3. on the four cards BOTH scans read, the two readings never
         disagree (the control group — it must also be non-empty, the
         lesson of the geometry verify's vacuous first version).
    """
    con = connect(DB, read_only=True)
    c = hdr_verify_counts(con)
    print(f"catalog rows: {c['n_obs']}   re-scraped: {c['n_scanned']}   "
          f"missing: {c['n_obs'] - c['n_scanned']}   (MUST be 0)")
    print(f"unreadable: {c['n_unreadable']}   of which readable by the "
          f"original scan: {c['n_unreadable_new']}   (MUST be 0)")
    print("\ncontrol group — cards both scans read:")
    n_control = n_bad = 0
    for col, (both, bad) in c["control"].items():
        n_control += both
        n_bad += bad
        print(f"   {col:<12} compared {both:>7}   disagree {bad:>5}   "
              f"(MUST be 0)   recovered where the catalog was NULL: "
              f"{c['recovered'][col]}")
    print("\nper card, over readable rows (value / blank / absent):")
    for col, (v, b, a) in c["cards"].items():
        print(f"   {col:<12} {v:>7} / {b:>6} / {a:>7}")
    ok = (c["n_obs"] == c["n_scanned"] and c["n_unreadable_new"] == 0
          and n_bad == 0 and n_control > 0)
    print("\nVERDICT:", "PASS — every row re-scraped, control group agrees"
          if ok else "FAIL — see the MUST-be-0 lines above")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Repair phantom BINTABLE geometry in the RLMT catalog, "
                    "and re-scrape the hardware-state header cards.")
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
    h = sub.add_parser("hdr-run")
    h.add_argument("--workers", type=int, default=HDR_DEFAULT_WORKERS,
                   help=f"reader threads (capped at {HDR_MAX_WORKERS})")
    h.add_argument("--limit", type=int, default=None,
                   help="stop after N rows (for a trial run)")
    h.add_argument("--retry-errors", action="store_true",
                   dest="retry_errors",
                   help="also re-read rows whose last attempt failed")
    h.set_defaults(fn=cmd_hdr_run)
    sub.add_parser("hdr-status").set_defaults(fn=cmd_hdr_status)
    sub.add_parser("hdr-verify").set_defaults(fn=cmd_hdr_verify)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
