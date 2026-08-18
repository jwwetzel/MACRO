#!/usr/bin/env python
"""Forecast the era consequences of the S0e geometry repair — on a COPY.

WHY A COPY
----------
``build_s0_manifest.py`` rebuilds the manifest into a temp file and
``os.replace``s it over the live path.  That swap replaces the whole FILE.
Running it while an S1 batch solve is writing ``s1_batch`` would pull the
database out from under the batch, so this script NEVER touches the live
manifest: it copies it, rebuilds against the copy, and diffs the two era
tables.  The live manifest is opened read-only, once, to be copied.

WHAT IT ANSWERS
---------------
1. Which era ids lose all their frames (the phantoms) once the geometry is
   right?
2. Which real eras absorb those frames?
3. Does any PUBLISHED era id change its meaning underneath a citation?

Question 3 is the one that matters.  Era ids are a pinned registry
(``build_s0_manifest.load_prior_era_ids``): reports, five strategy
documents and the ops request all cite "era 76", "era 47".  A published id
that silently came to mean a DIFFERENT camera configuration would corrupt
every one of those citations.  The registry is keyed on the configuration
tuple itself — ``(READOUTM, NAXIS1, NAXIS2, XBINNING, EGAIN)`` — so a
frame whose geometry is corrected moves to the id that already owns its
TRUE configuration, and the phantom id it left keeps its (now empty)
definition.  This script PROVES that rather than assuming it.

The forecast is written to the catalog DB as ``s0e_era_forecast`` so the
S0e report can render it from a table instead of from prose.

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/s0e_era_forecast.py --workdir <scratch dir>
"""

from __future__ import annotations

import argparse
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
LIVE_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
CATALOG = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite")
PYTHON = sys.executable
BUSY_TIMEOUT_MS = 300_000

#: The era-table columns that FORM the registry key.  Two eras are "the
#: same configuration" exactly when these five agree.
KEY_COLS = ("readoutm", "naxis1", "naxis2", "xbinning", "egain")


def era_map(db: Path) -> dict:
    """``{era_id: (key_tuple, n_frames, first_night, last_night)}``."""
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    try:
        return {int(r[0]): (tuple(r[1:6]), r[6], r[7], r[8])
                for r in con.execute(
                    "SELECT era_id, readoutm, naxis1, naxis2, xbinning, "
                    "egain, n_frames, first_night, last_night FROM eras")}
    finally:
        con.close()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", type=Path, required=True,
                    help="scratch directory for the manifest copy")
    ap.add_argument("--skip-rebuild", action="store_true",
                    help="reuse an existing rebuilt copy (for re-analysis)")
    args = ap.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)

    before_db = args.workdir / "manifest-before.sqlite"
    after_db = args.workdir / "manifest-after.sqlite"

    if not args.skip_rebuild:
        # 'before' = the live manifest as it stands (phantom geometry).
        # 'after'  = a second copy, rebuilt against the REPAIRED catalog.
        # Both are copies; the live file is only ever read.
        print(f"[S0e] copying live manifest -> {before_db}")
        shutil.copy2(LIVE_MANIFEST, before_db)
        print(f"[S0e] copying live manifest -> {after_db}")
        shutil.copy2(LIVE_MANIFEST, after_db)
        # Rebuilding with --out pointing at the COPY makes the copy both the
        # era-registry source (load_prior_era_ids) and the swap target, so
        # every published id is pinned exactly as a live rebuild would pin
        # it — and the live file is untouched.
        cmd = [PYTHON, str(REPO_ROOT / "pipeline" / "scripts"
                           / "build_s0_manifest.py"),
               "--catalog", str(CATALOG),
               "--out", str(after_db),
               "--eras-csv", str(args.workdir / "eras-after.csv"),
               "--skip-report"]
        print("[S0e] rebuilding against the copy:\n      " + " ".join(cmd))
        r = subprocess.run(cmd, cwd=str(REPO_ROOT))
        if r.returncode != 0:
            print("[S0e] rebuild FAILED", file=sys.stderr)
            return r.returncode

    before, after = era_map(before_db), era_map(after_db)
    print(f"\n[S0e] eras before: {len(before)}   after: {len(after)}")

    rows = []
    for era_id in sorted(set(before) | set(after)):
        b = before.get(era_id)
        a = after.get(era_id)
        bkey, bn = (b[0], b[1]) if b else (None, None)
        akey, an = (a[0], a[1]) if a else (None, None)
        if bkey is not None and akey is not None and bkey != akey:
            verdict = "REDEFINED"          # the serious case
        elif b is None:
            verdict = "NEW"
        elif a is None:
            verdict = "DROPPED"
        elif bn and not an:
            verdict = "RETIRED (emptied)"
        elif an and bn and an > bn:
            verdict = "ABSORBED (+frames)"
        elif an == bn:
            verdict = "unchanged"
        else:
            verdict = "shrank"
        rows.append((era_id, str(bkey), str(akey), bn, an, verdict))

    # --- the headline checks -------------------------------------------
    redefined = [r for r in rows if r[5] == "REDEFINED"]
    retired = [r for r in rows if r[5] == "RETIRED (emptied)"]
    absorbed = [r for r in rows if r[5] == "ABSORBED (+frames)"]

    print("\n=== ERA CHANGES ===")
    for r in rows:
        if r[5] != "unchanged":
            print(f"  era {r[0]:>3}  {r[5]:<20} "
                  f"n_frames {r[3]} -> {r[4]}\n"
                  f"        key before: {r[1]}\n"
                  f"        key after : {r[2]}")
    print(f"\nretired (emptied): {[r[0] for r in retired]}")
    print(f"absorbed (+frames): {[r[0] for r in absorbed]}")
    verdict_line = (
        "FAIL — these published ids changed meaning: "
        + str([r[0] for r in redefined]) if redefined
        else "PASS — no published era id was redefined")
    print(f"PUBLISHED-ID VERDICT: {verdict_line}")

    # --- persist the forecast so the report renders from a table --------
    con = sqlite3.connect(CATALOG)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    con.execute("DROP TABLE IF EXISTS s0e_era_forecast")
    con.execute("""
        CREATE TABLE s0e_era_forecast (
            era_id        INTEGER PRIMARY KEY,
            key_before    TEXT,
            key_after     TEXT,
            n_before      INTEGER,
            n_after       INTEGER,
            verdict       TEXT NOT NULL)""")
    con.executemany("INSERT INTO s0e_era_forecast VALUES (?,?,?,?,?,?)", rows)
    con.commit()
    print(f"\n[S0e] wrote s0e_era_forecast ({len(rows)} rows) to the catalog")

    # ------------------------------------------------------------------
    # The S1 re-queue population, measured on the REBUILT manifest.
    # ------------------------------------------------------------------
    # These are the frames the astrometry batch excluded as "sub-frame
    # photometry windows" and must now reconsider.  Computed by running the
    # project's OWN pure gates (macro_core.astrom) over the before- and
    # after-manifests, so the number on the report is the number the batch
    # driver would actually queue — not an independent re-derivation that
    # could drift from it.
    sys.path.insert(0, str(REPO_ROOT / "pipeline"))
    from macro_core import astrom

    def candidates(db: Path) -> dict:
        c = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
        c.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
        try:
            return {r["obs_rowid"]: r for r in astrom.fetch_candidates(c)}
        finally:
            c.close()

    b_rows, a_rows = candidates(before_db), candidates(after_db)
    newly = [a for rid, a in a_rows.items()
             if astrom.is_solvable_candidate(a)
             and rid in b_rows and not astrom.is_solvable_candidate(b_rows[rid])]
    print(f"[S0e] newly solvable frames: {len(newly)}")

    con.execute("DROP TABLE IF EXISTS s0e_requeue")
    con.execute("""
        CREATE TABLE s0e_requeue (
            obs_rowid        INTEGER PRIMARY KEY,
            path             TEXT,
            canonical_target TEXT,
            target_key       TEXT,
            readoutm         TEXT,
            filter           TEXT,
            night            TEXT,
            stratum_id       TEXT)""")
    con.executemany(
        "INSERT INTO s0e_requeue VALUES (?,?,?,?,?,?,?,?)",
        [(r["obs_rowid"], r["path"], r["canonical_target"], r["target_key"],
          r["readoutm"], r["filter"], r["night"],
          astrom.classify_stratum(r)) for r in newly])
    con.commit()
    n_unstrat = con.execute("SELECT COUNT(*) FROM s0e_requeue "
                            "WHERE stratum_id IS NULL").fetchone()[0]
    print(f"[S0e] wrote s0e_requeue ({len(newly)} rows); "
          f"{n_unstrat} carry NO stratum and would NOT be queued")
    con.close()
    return 1 if redefined else 0


if __name__ == "__main__":
    raise SystemExit(main())
