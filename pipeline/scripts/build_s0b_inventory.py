#!/usr/bin/env python
"""Build the S0b inventory: raw<->reduced links + calibration coverage.

WHAT THIS SCRIPT DOES (stage S0b of the shared pipeline)
--------------------------------------------------------
Reads the S0 manifest database, applies the pure S0b logic from
``macro_core.inventory`` to every frame, and AUGMENTS the manifest with four
NEW tables (the S0 tables — frames, aliases, eras, project_counts,
build_meta — are never touched):

* ``raw_reduced_links``  — one row per (raw canonical frame, reduced
                           counterpart) pair, with the match method that
                           proved the pair; reduced frames with no raw
                           parent appear with NULL raw columns (orphans).
* ``calib_frames``       — every calibration frame (bias/dark/flat,
                           normalized kind), with its S0 era_id, exposure
                           time, filter, night, temperatures, and a
                           master-product flag.
* ``calib_coverage``     — the matrix: for each era with canonical science
                           frames, one row per requirement (bias; dark per
                           science exposure time; flat per science filter)
                           with have-counts and a status.
* ``calib_gaps``         — the October shopping list: every requirement not
                           met, with an acquisition spec, the number of
                           science frames blocked, and the projects hit.
* ``s0b_build_meta``     — timestamp, code version, git commit, constants.

Unless ``--skip-report`` is given, it then renders the S0b evidence report
(``docs/pipeline/s0b_calibration_inventory.html`` + figures) from the
database it just wrote — the report reads ONLY the database, so every number
on the page is reproducible from the file alone.

IDEMPOTENCE / SAFETY
--------------------
Each S0b table is built under a temporary name and swapped into place inside
a single transaction (DROP old + RENAME new), so a crashed build never
leaves a half-written table and re-running is always safe.  This is also the
designed ingest path for the October re-opening: after each archive sync,
re-run S0 then S0b and the inventory (and its shopping list) updates in
place.

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s0b_inventory.py

That is the whole thing: the default points at the real manifest.  Add
``--help`` for every option.
"""

from __future__ import annotations

import argparse
import sqlite3
import subprocess
import sys
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# Make the pipeline package importable no matter where the script is invoked
# from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import S0B_CODE_VERSION                      # noqa: E402
from macro_core import inventory as inv                      # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

#: The frames-table columns S0b actually needs (the full table is wide;
#: loading only these keeps the build fast and the memory footprint small).
FRAME_COLUMNS = [
    "obs_rowid", "path", "tree", "basename", "jd", "night", "target_key",
    "canonical_target", "imagetyp", "filter", "exptime", "era_id",
    "is_canonical", "dup_group", "error", "readoutm", "camtemp", "ccd_temp",
]


# ---------------------------------------------------------------------------
# Step 1 — load the manifest frames (plus eras and project_counts)
# ---------------------------------------------------------------------------
def load_frames(con: sqlite3.Connection) -> pd.DataFrame:
    """Read the S0 frames table (selected columns) into a DataFrame."""
    cols = ", ".join(FRAME_COLUMNS)
    df = pd.read_sql_query(f"SELECT {cols} FROM frames", con)
    # Normalized calibration kind for EVERY row, once — the linkage, census,
    # and coverage steps all key off this single classification.
    df["calib_kind"] = [inv.calib_kind(it, bn)
                        for it, bn in zip(df["imagetyp"], df["basename"])]
    return df


def load_project_targets(con: sqlite3.Connection
                         ) -> tuple[dict, frozenset[str]]:
    """Project membership per target key, from the S0 project_counts table.

    Returns ``(project_of_key, dw_projects)`` for
    :func:`macro_core.inventory.projects_of_target`: explicit target keys
    map to their project sets; the ``__dw_survey__`` sentinel becomes the
    Dwarf-survey project set, applied by prefix rule.
    """
    rows = con.execute(
        "SELECT DISTINCT project, target_key FROM project_counts").fetchall()
    project_of_key: dict[str, set[str]] = {}
    dw_projects: set[str] = set()
    for project, tkey in rows:
        if tkey == "__dw_survey__":
            dw_projects.add(project)
        else:
            project_of_key.setdefault(tkey, set()).add(project)
    frozen = {k: frozenset(v) for k, v in project_of_key.items()}
    return frozen, frozenset(dw_projects)


# ---------------------------------------------------------------------------
# Step 2 — raw<->reduced links
# ---------------------------------------------------------------------------
def build_links(df: pd.DataFrame) -> pd.DataFrame:
    """Link every reduced-tree row to its raw canonical parent.

    Two mechanisms, mirroring the module docstring of
    ``macro_core.inventory``:

    * Rows S0 already deduplicated (same basename AND same JD as a canonical
      frame in another tree) link through their dup_group —
      method ``same_basename_jd``.
    * Everything else walks the pure match ladder
      (:func:`macro_core.inventory.link_reduced`): stem_jd, stem_jd_drift,
      target_jd, target_jd_ambiguous, or orphan.
    """
    red = df[df["tree"] == "reduced"]
    # Raw side of every link: canonical frames OUTSIDE the reduced tree.
    raw = df[(df["tree"] != "reduced") & (df["is_canonical"] == 1)]

    # ---- 2a. dup_group heads for the exact copies -------------------------
    # Map dup_group -> (raw obs_rowid, path) for groups whose canonical row
    # is a non-reduced frame.  (Exploration: all 79,719 non-canonical
    # reduced rows resolve this way; the code still handles the general
    # case — a reduced row whose group head is itself reduced simply falls
    # through to the ladder below.)
    head_of_group = dict(zip(raw["dup_group"], raw["obs_rowid"]))

    # ---- 2b. lookup tables for the pure match ladder ----------------------
    raw_by_stem: dict[str, list[tuple]] = {}
    raw_by_tjd: dict[tuple, list[tuple]] = {}
    for rid, bn, jd, night, tkey in zip(raw["obs_rowid"], raw["basename"],
                                        raw["jd"], raw["night"],
                                        raw["target_key"]):
        jd = None if pd.isna(jd) else float(jd)
        entry = (int(rid), jd, night)
        raw_by_stem.setdefault(inv.frame_stem(bn), []).append(entry)
        if tkey is not None and jd is not None:
            raw_by_tjd.setdefault((tkey, round(jd, 7)), []).append(entry)

    # ---- 2c. walk every reduced row ---------------------------------------
    rows = []
    for (rrid, rpath, rbn, rjd, rnight, rtkey, rgroup) in zip(
            red["obs_rowid"], red["path"], red["basename"], red["jd"],
            red["night"], red["target_key"], red["dup_group"]):
        rjd = None if pd.isna(rjd) else float(rjd)
        head = head_of_group.get(rgroup)
        if head is not None:
            # Exact copy: S0's own (basename, JD) dedup already proved it.
            matches = [(int(head), "same_basename_jd", 0.0)]
        else:
            matches = inv.link_reduced(
                inv.reduced_stem(rbn), rjd, rnight, rtkey,
                raw_by_stem, raw_by_tjd)
        if not matches:
            # Orphan: recorded with NULL raw columns, characterized in the
            # report — never silently dropped.
            rows.append({
                "reduced_rowid": int(rrid), "reduced_path": rpath,
                "raw_rowid": None, "raw_path": None,
                "match_method": "orphan", "jd": rjd, "night": rnight,
                "target_key": rtkey, "jd_drift_s": None,
            })
            continue
        for raw_id, method, drift in matches:
            rows.append({
                "reduced_rowid": int(rrid), "reduced_path": rpath,
                "raw_rowid": raw_id, "raw_path": None,   # filled below
                "match_method": method, "jd": rjd, "night": rnight,
                "target_key": rtkey, "jd_drift_s": drift,
            })
    links = pd.DataFrame(rows)
    # Resolve raw paths in one vectorized pass (dict.get, per the S0 lesson:
    # .map(dict) round-trips through a Series and corrupts odd keys).
    path_of_rowid = dict(zip(raw["obs_rowid"], raw["path"]))
    links["raw_path"] = links["raw_rowid"].map(path_of_rowid.get)
    return links


# ---------------------------------------------------------------------------
# Step 3 — calibration census
# ---------------------------------------------------------------------------
def build_calib_frames(df: pd.DataFrame) -> pd.DataFrame:
    """One row per canonical calibration frame, kind normalized.

    Duplicate copies are excluded the same way science dedup works: the
    canonical row represents the exposure.  era_id is S0's own assignment,
    joined straight from the frames table (never recomputed here).
    """
    cal = df[(df["calib_kind"].notna()) & (df["is_canonical"] == 1)].copy()
    cal["is_master"] = [1 if inv.is_master(bn) else 0
                        for bn in cal["basename"]]
    cal["exptime_bin"] = [inv.exptime_bin(x) for x in cal["exptime"]]
    out = cal[["obs_rowid", "path", "tree", "basename", "night", "jd",
               "era_id", "readoutm", "exptime", "exptime_bin", "filter",
               "camtemp", "ccd_temp"]].copy()
    out["kind"] = cal["calib_kind"]
    out["is_master"] = cal["is_master"]
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Step 4 — the coverage matrix and the shopping list
# ---------------------------------------------------------------------------
def select_science(df: pd.DataFrame) -> pd.DataFrame:
    """The science universe of the coverage matrix.

    Canonical, error-free, non-reduced-tree frames that
    :func:`macro_core.inventory.is_science` accepts (Light frames plus the
    blank-IMAGETYP 2026 nights).  The reduced tree is excluded because its
    canonical rows are overwhelmingly renamed copies of rawimage frames
    (section 1 of the S0b report) — counting them would double-count
    science and smear it across the packaging-artifact eras the
    decompressed copies occupy.
    """
    science_ok = pd.Series(
        [inv.is_science(it, k)
         for it, k in zip(df["imagetyp"], df["calib_kind"])],
        index=df.index)
    mask = ((df["is_canonical"] == 1)
            & (df["tree"] != "reduced")
            & df["error"].isna()
            & science_ok)
    return df[mask]


def build_coverage(df: pd.DataFrame,
                   project_of_key: dict[str, frozenset[str]],
                   dw_projects: frozenset[str],
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute calib_coverage and calib_gaps for every science era.

    For each era holding canonical science frames:

    * one ``bias`` row (biases are exposure- and filter-independent);
    * one ``dark`` row per science exposure-time bin, darks matched by the
      documented tolerance (:func:`macro_core.inventory.dark_matches`)
      against the bin value;
    * one ``flat`` row per science filter (exact filter-name match).

    Masters count separately from raw frames; the spec (and therefore the
    shopping list) is satisfied by raw frames only.  ``scaled_dark_ok`` is
    set on dark rows ONLY when the era has bias frames (a bias-subtracted
    dark can be exposure-scaled; without a bias the note would be wrong).
    """
    sci = select_science(df)
    cal = df[(df["calib_kind"].notna()) & (df["is_canonical"] == 1)].copy()
    cal["is_master"] = [1 if inv.is_master(bn) else 0
                        for bn in cal["basename"]]

    def _era_meta(era_frames: pd.DataFrame) -> str:
        # Camera description string from the era's own header values.
        r = era_frames.iloc[0]
        readout = r["readoutm"] if pd.notna(r["readoutm"]) and str(
            r["readoutm"]).strip() else "(blank READOUTM)"
        return readout

    cov_rows, gap_rows = [], []
    for era_id, sgrp in sci.groupby("era_id"):
        era_id = int(era_id)
        cgrp = cal[cal["era_id"] == era_id]
        n_sci_era = len(sgrp)
        nights = sgrp["night"].dropna()
        era_desc = _era_meta(sgrp)

        def _projects(frames: pd.DataFrame) -> str:
            hit: set[str] = set()
            for tk in frames["target_key"].dropna().unique():
                hit |= inv.projects_of_target(tk, project_of_key, dw_projects)
            return ",".join(sorted(hit))

        def _emit(req_kind: str, req_key, n_sci: int, n_raw: int,
                  n_master: int, spec_n: int, scaled_dark_ok, frames,
                  gap_eligible: bool = True):
            # gap_eligible=False keeps the coverage cell (the matrix stays
            # complete — nothing hidden) but suppresses the shopping-list
            # row: used for header-glitch filter strings that collide with
            # the calibration vocabulary ('flat dark' is not acquirable).
            status = inv.coverage_status(n_raw, n_master, spec_n)
            cov_rows.append({
                "era_id": era_id, "req_kind": req_kind, "req_key": req_key,
                "n_science": n_sci, "n_calib_raw": n_raw,
                "n_calib_master": n_master, "spec_n": spec_n,
                "status": status, "scaled_dark_ok": scaled_dark_ok,
            })
            if status != "ok" and gap_eligible:
                gap_rows.append({
                    "era_id": era_id, "camera": era_desc,
                    "first_night": nights.min() if len(nights) else None,
                    "last_night": nights.max() if len(nights) else None,
                    "need_kind": req_kind,
                    "spec": inv.gap_spec(req_kind, req_key, n_raw, spec_n),
                    "have_raw": n_raw, "have_master": n_master,
                    "status": status,
                    "n_science_frames_blocked": n_sci,
                    "projects_affected": _projects(frames),
                })

        # ---- bias: one requirement per era --------------------------------
        biases = cgrp[cgrp["calib_kind"] == "bias"]
        n_bias_raw = int((biases["is_master"] == 0).sum())
        n_bias_master = int((biases["is_master"] == 1).sum())
        _emit("bias", None, n_sci_era, n_bias_raw, n_bias_master,
              inv.SPEC_N_BIAS, None, sgrp)
        era_has_bias = (n_bias_raw + n_bias_master) > 0

        # ---- darks: one requirement per science exposure-time bin ---------
        darks = cgrp[cgrp["calib_kind"] == "dark"]
        sbins = pd.Series([inv.exptime_bin(x) for x in sgrp["exptime"]],
                          index=sgrp.index)
        for b in sorted({x for x in sbins if x is not None}):
            sub = sgrp[sbins == b]
            match = pd.Series([inv.dark_matches(x, b)
                               for x in darks["exptime"]], index=darks.index)
            n_raw = int(((darks["is_master"] == 0) & match).sum())
            n_master = int(((darks["is_master"] == 1) & match).sum())
            # Scaled-dark suitability note ONLY where the era has biases.
            scaled = 1 if (era_has_bias and len(darks) > 0) else None
            _emit("dark", inv.fmt_exptime(b), len(sub), n_raw, n_master,
                  inv.SPEC_N_DARK, scaled, sub)

        # ---- flats: one requirement per science filter --------------------
        flats = cgrp[cgrp["calib_kind"] == "flat"]
        filt = sgrp["filter"].fillna("(blank)")
        for f in sorted(filt.unique()):
            sub = sgrp[filt == f]
            fmatch = flats["filter"] == f
            n_raw = int(((flats["is_master"] == 0) & fmatch).sum())
            n_master = int(((flats["is_master"] == 1) & fmatch).sum())
            # A FILTER string that collides with the calibration vocabulary
            # ('dark'/'bias'/'flat' — a filter-wheel/header glitch) stays in
            # the matrix but never becomes a shopping-list row.
            _emit("flat", f, len(sub), n_raw, n_master,
                  inv.SPEC_N_FLAT, None, sub,
                  gap_eligible=not inv.is_calib_vocab_filter(f))

    coverage = pd.DataFrame(cov_rows)
    gaps = pd.DataFrame(gap_rows).sort_values(
        ["n_science_frames_blocked", "era_id"],
        ascending=[False, True]).reset_index(drop=True)
    return coverage, gaps


# ---------------------------------------------------------------------------
# Step 5 — atomic augmentation of the manifest
# ---------------------------------------------------------------------------
#: The tables this build owns.  ONLY these are ever dropped/replaced; the S0
#: tables are protected by the assertion in write_inventory().
S0B_TABLES = ("raw_reduced_links", "calib_frames", "calib_coverage",
              "calib_gaps", "s0b_build_meta")
S0_PROTECTED = frozenset({"frames", "aliases", "eras", "project_counts",
                          "build_meta"})


def write_inventory(manifest_path: Path, links: pd.DataFrame,
                    calib: pd.DataFrame, coverage: pd.DataFrame,
                    gaps: pd.DataFrame) -> None:
    """Swap the four S0b tables into the manifest, atomically per table set.

    Each table is written under a ``_s0b_tmp`` name first; one transaction
    then drops the old tables and renames the new ones into place — a reader
    (or a crash) can never observe a half-written S0b table, and the S0
    tables are never touched (asserted, not assumed).
    """
    assert not (set(S0B_TABLES) & S0_PROTECTED), \
        "an S0b table name collides with a protected S0 table"
    meta = pd.DataFrame([
        {"key": "built_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"key": "code_version", "value": S0B_CODE_VERSION},
        {"key": "git_commit", "value": _git_commit()},
        {"key": "dark_match_rel_tol", "value": str(inv.DARK_MATCH_REL_TOL)},
        {"key": "dark_match_abs_tol_s", "value": str(inv.DARK_MATCH_ABS_TOL)},
        {"key": "spec_n_bias", "value": str(inv.SPEC_N_BIAS)},
        {"key": "spec_n_dark", "value": str(inv.SPEC_N_DARK)},
        {"key": "spec_n_flat", "value": str(inv.SPEC_N_FLAT)},
    ])
    tables = {
        "raw_reduced_links": links,
        "calib_frames": calib,
        "calib_coverage": coverage,
        "calib_gaps": gaps,
        "s0b_build_meta": meta,
    }
    # closing() is required: sqlite3's own context manager only manages the
    # TRANSACTION — it does NOT close the file handle (the S0 lesson).
    with closing(sqlite3.connect(manifest_path)) as con:
        for name, frame in tables.items():
            con.execute(f"DROP TABLE IF EXISTS {name}_s0b_tmp")
            frame.to_sql(f"{name}_s0b_tmp", con, index=False)
        # The swap: one transaction, so every reader sees either the old
        # complete set or the new complete set, never a mixture.
        con.execute("BEGIN")
        for name in tables:
            con.execute(f"DROP TABLE IF EXISTS {name}")
            con.execute(f"ALTER TABLE {name}_s0b_tmp RENAME TO {name}")
        con.execute("COMMIT")
        # Indexes for the report's (and downstream stages') queries.
        con.execute("CREATE INDEX IF NOT EXISTS ix_links_raw "
                    "ON raw_reduced_links(raw_rowid)")
        con.execute("CREATE INDEX IF NOT EXISTS ix_links_red "
                    "ON raw_reduced_links(reduced_rowid)")
        con.execute("CREATE INDEX IF NOT EXISTS ix_calib_era "
                    "ON calib_frames(era_id)")
        con.execute("CREATE INDEX IF NOT EXISTS ix_cov_era "
                    "ON calib_coverage(era_id)")
        con.commit()


def _git_commit() -> str:
    """Best-effort short git hash of the repo (empty string off-repo)."""
    try:
        return subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build the S0b inventory (raw<->reduced links, calibration "
            "census, era coverage matrix, October shopping list) by "
            "augmenting the S0 manifest database with new tables, then "
            "render the evidence report. Safe to re-run: each S0b table is "
            "rebuilt and swapped in atomically; S0 tables are never "
            "modified. Re-run after every archive sync — this is the "
            "ingest path for the October observations."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="S0 manifest database to augment")
    p.add_argument("--skip-report", action="store_true",
                   help="build the tables only; do not render the HTML "
                        "report/figures afterwards")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.manifest.exists():
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 2

    print(f"[S0b] reading manifest {args.manifest} ...")
    with closing(sqlite3.connect(
            f"file:{args.manifest}?mode=ro", uri=True)) as con:
        df = load_frames(con)
        project_of_key, dw_projects = load_project_targets(con)
    print(f"[S0b]   {len(df):,} manifest rows; "
          f"{int(df['calib_kind'].notna().sum()):,} rows classify as "
          "calibration frames")

    print("[S0b] linking reduced tree to raw canonical frames ...")
    links = build_links(df)
    by_method = links["match_method"].value_counts()
    print("[S0b]   " + "; ".join(f"{k}: {v:,}" for k, v in by_method.items()))

    print("[S0b] building calibration census ...")
    calib = build_calib_frames(df)
    print(f"[S0b]   {len(calib):,} canonical calibration frames "
          f"({int((calib['is_master'] == 1).sum())} masters)")

    print("[S0b] computing era coverage matrix and gaps ...")
    coverage, gaps = build_coverage(df, project_of_key, dw_projects)
    print(f"[S0b]   {len(coverage):,} coverage cells across "
          f"{coverage['era_id'].nunique()} science eras; "
          f"{len(gaps):,} gaps")

    print(f"[S0b] writing inventory tables -> {args.manifest}")
    write_inventory(args.manifest, links, calib, coverage, gaps)

    if not args.skip_report:
        print("[S0b] rendering evidence report ...")
        from macro_core import report_s0b
        report_path = report_s0b.render_report(args.manifest)
        print(f"[S0b] report -> {report_path}")

    print("[S0b] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
