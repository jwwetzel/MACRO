#!/usr/bin/env python
"""Build the S0c per-project staging manifests — the working sets, no copies.

WHAT THIS SCRIPT DOES (stage S0c of the shared pipeline)
--------------------------------------------------------
Reads the S0/S0b manifest database and, for each of the five projects
(the selection rules are DATA in ``macro_core.staging.PROJECT_SELECTIONS``):

* selects the project's **science frames** (canonical, error-free Light
  frames of the project's published target list, per its filter rule);
* attaches the **calibration frames of every camera era that science
  touches** (bias/dark/flat, raw and recovered masters alike, from the S0b
  census), each row carrying an explicit ``match_basis``;
* writes the combined list to a ``stage_<project>`` table in the manifest
  (atomic swap, S0/S0b tables never touched) **and** to
  ``<repo>/<Project>/data/stage_manifest.csv`` (atomic os.replace) with a
  ``data/README.md`` explaining the no-copy law and the columns.

THE NO-COPY LAW (why there is no 'copy frames' option)
------------------------------------------------------
S0 exists because copies proliferated — 132k duplicate rows in the catalog.
No frame is ever copied into a project directory: the staging manifest IS
the working set, and stages read the immutable archive directly through it.
The only materialization on offer is ``--symlink-farm`` (default OFF), a
disposable browsable view for humans:

* ``<Project>/data/frames/<role>/<night>_<basename> -> archive file``
* symlinks only, pointing INTO the archive; delete the farm any time and
  regenerate it from the CSV;
* Dropbox does NOT sync symlink targets — the farm is a local convenience,
  never a transport mechanism.

Unless ``--skip-report`` is given, the S0c evidence report
(``docs/pipeline/s0c_staging.html`` + figures) is rendered from the database
just written — every number on the page is a SQL query result.

IDEMPOTENCE / SAFETY
--------------------
Stage tables are built under temporary names and swapped in inside one
transaction; CSVs are written to a temp file and os.replace()d.  Re-running
after an archive sync (the October ingest path: S0 → S0b → S0c) refreshes
everything in place.  The archive itself is opened by nobody here: this
script reads ONLY the manifest database and writes ONLY under the project
``data/`` directories, the manifest, and ``docs/``.

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s0c_staging.py

Add ``--symlink-farm`` for the browsable view, ``--help`` for everything.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
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

from macro_core import S0C_CODE_VERSION                      # noqa: E402
from macro_core import manifest as mf                        # noqa: E402
from macro_core import staging as stg                        # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

#: frames-table columns the science selection needs (the table is wide;
#: loading only these keeps the build fast).
FRAME_COLUMNS = [
    "obs_rowid", "path", "tree", "basename", "size", "jd", "night",
    "target_key", "canonical_target", "imagetyp", "filter", "exptime",
    "era_id", "is_canonical", "dup_group", "error", "qc_flags",
    "pointing_offset_deg",
    # Coordinates: needed only by the cone-candidate clause (a selection with
    # cone_radius_deg set), but loaded always so the frame loader has one
    # shape regardless of which projects this run stages.
    "ra_deg", "dec_deg",
]

#: Tables this build owns (plus the five stage_<project> tables, derived
#: from PROJECT_SELECTIONS below).  S0/S0b tables are protected.
S0C_FIXED_TABLES = ("s0c_stage_files", "s0c_build_meta")
PROTECTED_TABLES = frozenset({
    "frames", "aliases", "eras", "project_counts", "build_meta",
    "raw_reduced_links", "calib_frames", "calib_coverage", "calib_gaps",
    "s0b_build_meta",
})


# ---------------------------------------------------------------------------
# Step 1 — load the manifest
# ---------------------------------------------------------------------------
def load_frames(con: sqlite3.Connection) -> pd.DataFrame:
    """Read the S0 frames table (selected columns) into a DataFrame."""
    cols = ", ".join(f'"{c}"' for c in FRAME_COLUMNS)
    return pd.read_sql_query(f"SELECT {cols} FROM frames", con)


def load_calib(con: sqlite3.Connection) -> pd.DataFrame:
    """Read the S0b calibration census, with each frame's catalog size.

    ``calib_frames`` carries no size column; the join back to ``frames``
    on obs_rowid recovers it (every catalog row has a size — asserted).
    """
    df = pd.read_sql_query(
        """SELECT c.obs_rowid, c.path, c.tree, c.night, c.jd, c.era_id,
                  c.exptime, c."filter", c.kind, c.is_master, f.size
           FROM calib_frames c JOIN frames f ON f.obs_rowid = c.obs_rowid""",
        con)
    assert df["size"].notna().all(), \
        "calib frame with no catalog size — S0 scan incomplete?"
    return df


# ---------------------------------------------------------------------------
# Step 2 — build one project's staging manifest
# ---------------------------------------------------------------------------
def reference_positions(sci: pd.DataFrame,
                        min_frames: int = 3) -> dict[str, tuple[float, float]]:
    """Median (RA, Dec) per staged target, from the staged science rows.

    The cone clause needs a position per target and must NOT invent one: it
    uses the positions of the frames the project already staged, so a target
    with too few usable coordinates simply has no cone and matches nothing.
    ``min_frames`` mirrors ``manifest.MIN_SOLVED_FOR_REFERENCE`` in spirit —
    one or two frames are not a position.  RA-wrap-aware via S0's own
    :func:`macro_core.manifest.median_radec`.
    """
    refs: dict[str, tuple[float, float]] = {}
    have = sci[sci["ra_deg"].notna() & sci["dec_deg"].notna()
               & sci["target_key"].notna()]
    for key, grp in have.groupby("target_key"):
        if len(grp) >= min_frames:
            refs[str(key)] = mf.median_radec(grp["ra_deg"], grp["dec_deg"])
    return refs


def build_project_stage(sel: stg.ProjectSelection, frames: pd.DataFrame,
                        calib: pd.DataFrame, archive_root: str,
                        build_id: str) -> pd.DataFrame:
    """Science + cone-candidate + era-matched calibration rows for a project.

    Row order is deterministic (science by night/jd/path, then calibration
    by role/era/night/path) so re-running the build on unchanged inputs
    reproduces the CSV byte-for-byte after the build-id line.
    """
    # ---- science: apply the pure selection predicate row by row ----------
    mask = [
        stg.is_staged_science(
            sel, tk, it, err, canon, tree, flt, bn or "")
        for tk, it, err, canon, tree, flt, bn in zip(
            frames["target_key"], frames["imagetyp"], frames["error"],
            frames["is_canonical"], frames["tree"], frames["filter"],
            frames["basename"])
    ]
    sci = frames[pd.Series(mask, index=frames.index)]
    sci_rows = [stg.science_row(sel, rec, archive_root, build_id)
                for rec in sci.to_dict("records")]

    # ---- cone candidates: name-less frames pointed at a staged target ----
    # Only runs for a selection that asked for it (BeStar_Grism).  The
    # reference positions come from the science rows just selected, so the
    # cone can never point somewhere the project does not already observe.
    cone_rows: list[dict] = []
    if sel.cone_radius_deg is not None:
        refs = reference_positions(sci)
        eligible = [
            stg.is_cone_candidate(
                sel, tk, it, err, canon, tree, flt, ra, dec, bn or "")
            for tk, it, err, canon, tree, flt, ra, dec, bn in zip(
                frames["target_key"], frames["imagetyp"], frames["error"],
                frames["is_canonical"], frames["tree"], frames["filter"],
                frames["ra_deg"], frames["dec_deg"], frames["basename"])
        ]
        for rec in frames[pd.Series(eligible, index=frames.index)] \
                .to_dict("records"):
            hit = stg.cone_match(rec["ra_deg"], rec["dec_deg"], refs,
                                 sel.cone_radius_deg)
            if hit is not None:
                cone_rows.append(stg.cone_candidate_row(
                    rec, hit[0], hit[1], archive_root, build_id))

    # ---- calibration: every frame of every era the science touches -------
    eras = sorted({int(e) for e in sci["era_id"].dropna().unique()})
    cal = calib[calib["era_id"].isin(eras)]
    cal_rows = [stg.calib_row(rec, archive_root, build_id)
                for rec in cal.to_dict("records")]

    out = pd.DataFrame(sci_rows + cone_rows + cal_rows,
                       columns=stg.STAGE_CSV_COLUMNS)
    # Deterministic order: science first, then cone candidates, then
    # calibration grouped by role/era.  A stable sort on a synthetic rank
    # (science 0, unresolved 1, calibration 2) does all three at once.
    rank = {stg.ROLE_SCIENCE: 0, stg.ROLE_SCIENCE_UNRESOLVED: 1}
    out["_sci"] = [rank.get(r, 2) for r in out["role"]]
    out = out.sort_values(
        ["_sci", "role", "era_id", "night", "jd", "path"],
        na_position="last", kind="mergesort").drop(columns="_sci")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Step 3 — atomic writes: manifest tables, CSVs, READMEs
# ---------------------------------------------------------------------------
def write_stage_tables(manifest_path: Path,
                       stages: dict[str, pd.DataFrame],
                       stage_files: pd.DataFrame,
                       build_id: str) -> None:
    """Swap the S0c tables into the manifest atomically.

    Same discipline as S0b: build under ``_s0c_tmp`` names, then one
    transaction drops the old set and renames the new set into place — a
    reader never observes a half-written staging table, and the S0/S0b
    tables are never touched (asserted, not assumed).
    """
    tables: dict[str, pd.DataFrame] = {
        stg.stage_table_name(p): df for p, df in stages.items()}
    meta = pd.DataFrame([
        {"key": "built_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"key": "code_version", "value": S0C_CODE_VERSION},
        {"key": "build_id", "value": build_id},
        {"key": "git_commit", "value": _git_commit()},
        {"key": "staging_sha256_12", "value": _staging_source_hash()},
        {"key": "archive_root", "value": stg.DEFAULT_ARCHIVE_ROOT},
        {"key": "checksum_note", "value": stg.CHECKSUM_NOTE},
        {"key": "match_basis_science", "value": stg.MATCH_BASIS_SCIENCE},
        {"key": "match_basis_calib", "value": stg.MATCH_BASIS_CALIB},
    ])
    tables["s0c_stage_files"] = stage_files
    tables["s0c_build_meta"] = meta
    assert not (set(tables) & PROTECTED_TABLES), \
        "an S0c table name collides with a protected S0/S0b table"

    # closing() is required: sqlite3's own context manager only manages the
    # TRANSACTION — it does NOT close the file handle (the S0 lesson).
    # The generous busy timeout lets this build coexist with a concurrent
    # S1/S2 campaign writing its own tables: sqlite serializes the writers
    # instead of failing fast with 'database is locked'.
    with closing(sqlite3.connect(manifest_path, timeout=300.0)) as con:
        for name, frame in tables.items():
            con.execute(f"DROP TABLE IF EXISTS {name}_s0c_tmp")
            frame.to_sql(f"{name}_s0c_tmp", con, index=False)
        con.execute("BEGIN")
        for name in tables:
            con.execute(f"DROP TABLE IF EXISTS {name}")
            con.execute(f"ALTER TABLE {name}_s0c_tmp RENAME TO {name}")
        con.execute("COMMIT")
        # Indexes for the report's and downstream stages' common queries.
        for name in (stg.stage_table_name(p) for p in stages):
            con.execute(f"CREATE INDEX IF NOT EXISTS ix_{name}_role "
                        f"ON {name}(role)")
            con.execute(f"CREATE INDEX IF NOT EXISTS ix_{name}_era "
                        f"ON {name}(era_id)")
        con.commit()


def write_csv_atomic(csv_path: Path, df: pd.DataFrame) -> None:
    """Write one stage_manifest.csv via a temp file + os.replace.

    ``os.replace`` is atomic on the same filesystem, so a crash mid-write
    can never leave a truncated CSV where a stage expects a complete one.
    """
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = csv_path.with_suffix(".csv.tmp")
    df.to_csv(tmp, index=False, quoting=csv.QUOTE_MINIMAL)
    os.replace(tmp, csv_path)


#: The regeneration command printed in every data/README.md.
#:
#: ABSOLUTE, and quoted.  The first version printed the repo-relative script
#: path with no working-directory instruction, so a reader who ran it from
#: the directory the README lives in got a file-not-found — the one
#: invocation a student auditor is guaranteed to copy.  The script itself is
#: location-independent, so only the documented command was ever broken.
#: The quotes are load-bearing: the repo path contains spaces.
REGEN_COMMAND = (f'/opt/miniconda3/envs/rlmt-checks/bin/python \\\n'
                 f'        "{REPO_ROOT}/pipeline/scripts/'
                 f'build_s0c_staging.py"')


def readme_text(sel: stg.ProjectSelection, n_science: int, n_calib: int,
                build_id: str, n_cone: int = 0) -> str:
    """The data/README.md for one project — law, columns, regeneration."""
    cols = "\n".join(f"| `{c}` | " + {
        "path": "archive-relative POSIX path — the frame's identity",
        "abs_path": "absolute archive path (QUOTE IT: the root has spaces)",
        "role": (f"`{stg.ROLE_SCIENCE}`, "
                 f"`{stg.ROLE_SCIENCE_UNRESOLVED}` (cone candidate — NOT "
                 "science until a project adjudicates it), "
                 "`bias`/`dark`/`flat`, or `master_*` products"),
        "match_basis": (f"`{stg.MATCH_BASIS_SCIENCE}` (science: the rule "
                        f"below), `{stg.MATCH_BASIS_CONE}` (no target name; "
                        "matched by coordinates) or "
                        f"`{stg.MATCH_BASIS_CALIB}` (calibration: "
                        "same S0 era as this project's science)"),
        "tree": "top-level archive tree holding the canonical copy",
        "era_id": "S0 pinned camera-era registry id",
        "night": "local-noon-to-noon night label",
        "jd": "header JD = **UTC exposure START** (BJD_TDB is stage S3's "
              "job — never use this for timing)",
        "filter": "cataloged filter string",
        "exptime": "header EXPTIME (s)",
        "canonical_target": "S0 alias-merged display name (science rows)",
        "target_key": "S0 normalized target key (science rows)",
        "dup_group": "S0 global duplicate-group id",
        "qc_flags": "S0 QC flags — flags mark, they never delete",
        "pointing_offset_deg": "offset from the target's reference position",
        "size_bytes": "integrity surrogate (see note below)",
        "obs_rowid": "catalog/manifest join key",
        "stage_build_id": "S0c build that emitted the row",
    }[c] + " |" for c in stg.STAGE_CSV_COLUMNS)
    # The cone paragraph appears ONLY for a selection that enables the cone
    # clause, so the other four READMEs never document a rule they do not run.
    cone_note = "" if sel.cone_radius_deg is None else f"""

**Cone candidates (`role = '{stg.ROLE_SCIENCE_UNRESOLVED}'`).** Frames that
pass every other gate but carry NO target name enter the working set when
their coordinates fall within **{sel.cone_radius_deg:g}°** of a staged
target's reference position (the median position of that target's own staged
science frames). They are `match_basis = '{stg.MATCH_BASIS_CONE}'`, their
`target_key` is the CANDIDATE match and their `pointing_offset_deg` is the
measured separation. **They are not science.** A stage that wants them must
ask for the role by name and adjudicate them first — that adjudication is
this project's Step 0, and it now happens inside the manifest instead of by
querying `frames` behind S0c's back."""
    return f"""# {sel.project} — staging manifest (`stage_manifest.csv`)

**THE NO-COPY LAW.** No frame is ever copied into this directory. S0 exists
because copies proliferated (132k duplicate rows in the archive catalog);
this manifest **is** the working set. Every pipeline stage reads the
immutable archive directly through the paths below — the archive is
read-only, always.

**This file is regenerable, not precious** (`*/data/` is gitignored). Run
this from anywhere — the path is absolute and quoted because the repo path
contains spaces:

    {REGEN_COMMAND}

**Selection rule (science rows).** {sel.rule}
Source: {sel.source}{cone_note}

**Calibration rows.** For every camera era the science frames touch, ALL of
that era's calibration frames from the S0b census are included (raw frames
and recovered `Calibrations/` masters alike), `match_basis =
'{stg.MATCH_BASIS_CALIB}'`. Staging deliberately over-includes; each stage
narrows by kind/exposure/filter with the S0b coverage matrix as its guide.

**This build ({build_id}):** {n_science:,} science rows +
{n_cone:,} cone-candidate rows + {n_calib:,} calibration rows.

## Columns

| column | meaning |
|---|---|
{cols}

**Integrity note.** {stg.CHECKSUM_NOTE}

## Optional symlink farm (`frames/`)

`build_s0c_staging.py --symlink-farm` materializes
`frames/<role>/<night>_<basename>` symlinks into the archive for humans who
want a browsable view. The farm is **disposable** (delete and regenerate at
will) and **Dropbox does not sync symlink targets** — it is a local
convenience, never a transport mechanism.
"""


def write_readme(data_dir: Path, text: str) -> None:
    """Atomic README write (same temp + replace discipline as the CSV)."""
    tmp = data_dir / "README.md.tmp"
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, data_dir / "README.md")


# ---------------------------------------------------------------------------
# Step 4 — the optional, disposable symlink farm
# ---------------------------------------------------------------------------
def build_symlink_farm(data_dir: Path, df: pd.DataFrame) -> int:
    """Materialize data/frames/<role>/ symlinks into the archive.

    The farm is rebuilt from scratch each time (old links removed first) so
    it can never drift from the CSV.  Only symlinks are ever created —
    nothing is copied, nothing in the archive is touched.  Returns the
    number of links created.
    """
    farm = data_dir / "frames"
    # Remove ONLY symlinks (and then-empty dirs) — belt and braces against
    # someone having parked a real file inside the farm.
    if farm.exists():
        for role_dir in sorted(farm.iterdir()):
            if not role_dir.is_dir():
                continue
            for link in sorted(role_dir.iterdir()):
                if link.is_symlink():
                    link.unlink()
            try:
                role_dir.rmdir()
            except OSError:
                pass  # non-link content left behind on purpose
    n = 0
    for rec in df.to_dict("records"):
        role_dir = farm / rec["role"]
        role_dir.mkdir(parents=True, exist_ok=True)
        link = role_dir / stg.farm_link_name(rec["night"],
                                             Path(rec["path"]).name)
        if not link.exists() and not link.is_symlink():
            link.symlink_to(rec["abs_path"])
            n += 1
    return n


def _git_commit() -> str:
    """Short git hash of the repo, marked ``-dirty`` when the tree is not.

    THE 17ef904 LESSON (2026-08-18 review).  The first S0c build stamped
    ``git_commit='17ef904'`` — the commit BEFORE the one that introduced
    ``staging.py``.  The build had run from a dirty working tree, and a bare
    ``rev-parse HEAD`` records the last COMMIT, not the code that ran.  For a
    product whose entire claim is that no number is hand-typed, a provenance
    field that points at a tree without the selection rules in it is worse
    than no field.  ``git status --porcelain`` is cheap; when it prints
    anything, the hash is suffixed ``-dirty`` and
    :func:`_staging_source_hash` records what actually ran.
    """
    try:
        head = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""
    try:
        dirty = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "status", "--porcelain"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return head            # HEAD known, cleanliness unknown: say no more
    return f"{head}-dirty" if dirty else head


def _staging_source_hash() -> str:
    """SHA-256 (first 12 hex) of the file that encodes the selection rules.

    The commit hash identifies the repo; this identifies the ONE file whose
    contents decide which frames every project gets.  Recorded next to the
    commit so a rebuild's rows can be traced to the exact rule text that
    produced them even when the build ran from an uncommitted tree.
    """
    try:
        source = (PIPELINE_ROOT / "macro_core" / "staging.py").read_bytes()
    except OSError:
        return ""
    return hashlib.sha256(source).hexdigest()[:12]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Build the S0c per-project staging manifests: for each of the "
            "five projects, a provenance-complete frame list (science + "
            "era-matched calibration) written as a stage_<project> table "
            "in the manifest AND <Project>/data/stage_manifest.csv. No "
            "frames are ever copied — the manifest IS the working set. "
            "Safe to re-run; part of the S0 -> S0b -> S0c ingest chain."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                   help="S0/S0b manifest database to read and augment")
    p.add_argument("--archive-root", default=stg.DEFAULT_ARCHIVE_ROOT,
                   help="immutable archive root the abs_path column points "
                        "into (read-only; never written)")
    p.add_argument("--project", action="append", default=None,
                   metavar="NAME",
                   help="stage only this project (repeatable; default all "
                        "five)")
    p.add_argument("--symlink-farm", action="store_true",
                   help="also materialize data/frames/<role>/ symlinks "
                        "into the archive (disposable browsable view; "
                        "Dropbox does not sync symlink targets)")
    p.add_argument("--skip-report", action="store_true",
                   help="build tables/CSVs only; skip the HTML report")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.manifest.exists():
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 2
    selections = list(stg.PROJECT_SELECTIONS)
    if args.project:
        try:
            selections = [stg.selection_for(p) for p in args.project]
        except KeyError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 2

    build_id = (f"{S0C_CODE_VERSION} @ "
                f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%MZ')}")

    print(f"[S0c] reading manifest {args.manifest} ...")
    with closing(sqlite3.connect(
            f"file:{args.manifest}?mode=ro", uri=True)) as con:
        frames = load_frames(con)
        calib = load_calib(con)
    print(f"[S0c]   {len(frames):,} manifest rows; "
          f"{len(calib):,} calibration census rows")

    stages: dict[str, pd.DataFrame] = {}
    file_rows = []
    for sel in selections:
        df = build_project_stage(sel, frames, calib,
                                 args.archive_root, build_id)
        # Count by role explicitly.  The old `len(df) - n_sci` idiom silently
        # counted cone candidates as calibration the moment a third role
        # existed — the reason SCIENCE_ROLES is a shared constant.
        n_sci = int((df["role"] == stg.ROLE_SCIENCE).sum())
        n_cone = int((df["role"] == stg.ROLE_SCIENCE_UNRESOLVED).sum())
        n_cal = int((~df["role"].isin(stg.SCIENCE_ROLES)).sum())
        assert n_sci + n_cone + n_cal == len(df), "a row escaped its role"
        stages[sel.project] = df

        data_dir = REPO_ROOT / sel.project / "data"
        csv_path = data_dir / "stage_manifest.csv"
        write_csv_atomic(csv_path, df)
        write_readme(data_dir,
                     readme_text(sel, n_sci, n_cal, build_id, n_cone))
        n_links = 0
        if args.symlink_farm:
            n_links = build_symlink_farm(data_dir, df)

        file_rows.append({
            "project": sel.project,
            "stage_table": stg.stage_table_name(sel.project),
            "csv_path": str(csv_path.relative_to(REPO_ROOT)),
            "n_rows": len(df), "n_science": n_sci, "n_calib": n_cal,
            "n_cone": n_cone,
            "n_eras": int(df["era_id"].dropna().nunique()),
            "n_symlinks": n_links,
            "selection_rule": sel.rule, "selection_source": sel.source,
        })
        print(f"[S0c]   {sel.project}: {n_sci:,} science + {n_cone:,} cone "
              f"+ {n_cal:,} calib rows -> {csv_path.relative_to(REPO_ROOT)}"
              + (f" (+{n_links:,} farm links)" if args.symlink_farm else ""))

    stage_files = pd.DataFrame(file_rows)
    # A --project run rebuilds SOME stage tables; the registry must keep the
    # other projects' rows (their tables are untouched), or the report would
    # silently forget them.  Merge: this run's rows win, existing rows for
    # projects outside this run survive.
    built = {r["project"] for r in file_rows}
    with closing(sqlite3.connect(
            f"file:{args.manifest}?mode=ro", uri=True)) as con:
        has_registry = con.execute(
            "SELECT count(*) FROM sqlite_master "
            "WHERE type='table' AND name='s0c_stage_files'").fetchone()[0]
        if has_registry:
            existing = pd.read_sql_query(
                "SELECT * FROM s0c_stage_files", con)
            keep = existing[~existing["project"].isin(built)]
            if len(keep):
                stage_files = pd.concat([stage_files, keep],
                                        ignore_index=True)
    stage_files = stage_files.sort_values("project").reset_index(drop=True)
    # A registry row carried from a build that pre-dates the cone clause has
    # no n_cone value; it had no cone rows either, so zero is the truth.
    if "n_cone" in stage_files:
        stage_files["n_cone"] = stage_files["n_cone"].fillna(0).astype(int)
    print(f"[S0c] writing stage tables -> {args.manifest}")
    write_stage_tables(args.manifest, stages, stage_files, build_id)

    if not args.skip_report:
        print("[S0c] rendering evidence report ...")
        from macro_core import report_s0c
        report_path = report_s0c.render_report(args.manifest)
        print(f"[S0c] report -> {report_path}")

    print("[S0c] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
