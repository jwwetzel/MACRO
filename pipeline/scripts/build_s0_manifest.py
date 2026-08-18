#!/usr/bin/env python
"""Build the S0 manifest database for the MACRO/RLMT archive.

WHAT THIS SCRIPT DOES (stage S0 of the shared pipeline, ROADMAP.md sec. 1.1)
---------------------------------------------------------------------------
Reads the observation catalog (READ-ONLY — this script never writes to it),
applies the pure S0 logic from ``macro_core.manifest`` to every row, and
writes a fresh manifest database with five tables:

* ``frames``          — one row per catalog row, plus: basename, night label,
                        era_id, duplicate group id, canonical flag, resolved
                        target name, pointing offset, QC flags.
* ``aliases``         — every raw target name → its canonical name, with the
                        exact normalization rules that fired and the result
                        of the coordinate-cone audit.
* ``eras``            — the camera-era table keyed on (READOUTM, geometry,
                        binning, EGAIN), also exported as CSV.
* ``project_counts``  — per-project canonical-frame counts next to the
                        numbers each strategy document claims (section 7 of
                        the report renders this reconciliation).
* ``build_meta``      — timestamp, catalog path, code version, git commit.

Unless ``--skip-report`` is given, it then renders the chain-of-evidence
report (``docs/pipeline/s0_manifest.html`` + figures) from the database it
just wrote — the report reads ONLY the manifest, never the catalog, so every
number on the page is reproducible from the manifest alone.

IDEMPOTENCE / SAFETY
--------------------
The manifest is rebuilt from scratch on every run and swapped into place
atomically (write to a temp file in the same directory, then ``os.replace``),
so a crashed or interrupted build can never leave a half-written manifest,
and re-running is always safe.

USAGE (a student's quick start)
-------------------------------
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s0_manifest.py

That is the whole thing: the defaults point at the real catalog and the real
repo.  Add ``--help`` for every option.
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import subprocess
import sys
import tempfile
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# Make the pipeline package importable no matter where the script is invoked
# from: the package root is the parent of this script's directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import S0_CODE_VERSION                       # noqa: E402
from macro_core import manifest as m                         # noqa: E402

# ---------------------------------------------------------------------------
# Default locations (real paths, so the bare command Just Works).
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_CATALOG = Path(
    "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite")
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ERAS_CSV = REPO_ROOT / "products" / "manifest" / "eras.csv"


# ---------------------------------------------------------------------------
# Step 1 — load the catalog (read-only)
# ---------------------------------------------------------------------------
def load_catalog(catalog_path: Path) -> pd.DataFrame:
    """Read every ``obs`` row into a DataFrame, catalog opened read-only.

    The ``mode=ro`` URI guarantees the ground-truth catalog cannot be
    modified even by accident (ROADMAP convention 1).
    """
    uri = f"file:{catalog_path}?mode=ro"
    with sqlite3.connect(uri, uri=True) as con:
        # rowid gives every catalog row a stable integer identity that we
        # carry into the manifest (useful for tracing a frame back).
        df = pd.read_sql_query("SELECT rowid AS obs_rowid, * FROM obs", con)
    return df


# ---------------------------------------------------------------------------
# Step 2 — alias resolution (normalization rules + cone-gated synonyms)
# ---------------------------------------------------------------------------
def resolve_aliases(df: pd.DataFrame) -> tuple[pd.DataFrame, dict, dict]:
    """Resolve every raw target name to a canonical alias group.

    Returns
    -------
    aliases : DataFrame
        One row per distinct raw ``target_best`` value: the canonical name,
        frame count, rules applied, and the coordinate-cone audit result.
    key_of_raw : dict
        raw name → final normalized key (None for blank names).
    display_of_key : dict
        final key → canonical display name (the most frequent raw variant).
    """
    # ---- 2a. run the pure normalizer over every distinct raw name --------
    raw_counts = df["target_best"].fillna("").value_counts()
    norm = {raw: m.normalize_target(raw if raw else None)
            for raw in raw_counts.index}

    # ---- 2b. median plate-solved coordinates per raw name ----------------
    # Used both to gate synonym merges and to audit every merge afterwards.
    solved = df[(df["pltsolvd"] == 1)
                & df["ra_deg"].notna() & df["dec_deg"].notna()]
    coords_of_raw: dict[str, tuple[float, float]] = {}
    for raw, grp in solved.groupby(solved["target_best"].fillna("")):
        coords_of_raw[raw] = m.median_radec(grp["ra_deg"], grp["dec_deg"])

    # ---- 2c. cone-gate every synonym-table merge -------------------------
    # A synonym merge joins names no string rule can relate, so it must be
    # *verified* on the sky: the merged-in side and the destination side
    # must agree within CONE_RADIUS_DEG where both have solved coordinates.
    # A failed gate refuses the merge (the names stay separate) and the
    # refusal is recorded in the aliases table.
    refused_synonyms: set[str] = set()
    for raw, info in norm.items():
        if "synonym" not in info.rules:
            continue
        src_key = info.pre_synonym_key       # key before the synonym fired
        dst_key = info.key                   # key after
        # Coordinates of everything that formed the destination group
        # WITHOUT the synonym rule (the natives of dst_key):
        native_coords = [coords_of_raw[r] for r, i in norm.items()
                         if i.key == dst_key and "synonym" not in i.rules
                         and r in coords_of_raw]
        here = coords_of_raw.get(raw)
        if here is None or not native_coords:
            # No solved coordinates on one side → the gate cannot run;
            # the merge stands (it is an explicit, documented table entry)
            # and the audit column will show NULL for this name.
            continue
        ra0, dec0 = m.median_radec([c[0] for c in native_coords],
                                   [c[1] for c in native_coords])
        sep = m.angular_separation_deg(here[0], here[1], ra0, dec0)
        if sep > m.CONE_RADIUS_DEG:
            # Gate FAILED: undo the merge for every raw name that mapped
            # through this same synonym source key.
            refused_synonyms.add(src_key)

    key_of_raw: dict[str, str | None] = {}
    for raw, info in norm.items():
        if "synonym" in info.rules and info.pre_synonym_key in refused_synonyms:
            key_of_raw[raw] = info.pre_synonym_key   # merge refused
        else:
            key_of_raw[raw] = info.key

    # ---- 2d. canonical display name per key ------------------------------
    # Vote with the *cleaned* form of each raw name (junk and leaked tokens
    # already stripped), weighted by row count: 'PHECDA lrg 0-25s' votes
    # for 'PHECDA', so the display name is always a real name even when the
    # leaked variants outnumber the clean one.  Count ties break
    # lexicographically for run-to-run determinism.
    #
    # Synonym destinations are special: the SYNONYM_TABLE documents a merge
    # DIRECTION (e.g. alphalyr -> vega), so the displayed name must be the
    # destination's own — a NATIVE name, one whose pre-synonym key already
    # equals the group key.  Without this, a merged-in name with more rows
    # would out-vote the destination ('Alpha Lyr' at 794 rows would defeat
    # 'Vega'), inverting the documented arrow.  We therefore keep a second,
    # natives-only ballot and let it override for synonym destinations.
    display_of_key: dict[str, str] = {}
    cleaned_votes: dict[str, dict[str, int]] = {}
    native_votes: dict[str, dict[str, int]] = {}
    for raw, n in raw_counts.items():
        key = key_of_raw.get(raw)
        if key is None:
            continue
        cleaned = norm[raw].cleaned or raw
        bucket = cleaned_votes.setdefault(key, {})
        bucket[cleaned] = bucket.get(cleaned, 0) + int(n)
        # Native = the synonym table did not move this name here: its key
        # before the synonym rule already equals the final group key.
        if norm[raw].pre_synonym_key == key:
            nb = native_votes.setdefault(key, {})
            nb[cleaned] = nb.get(cleaned, 0) + int(n)
    synonym_destinations = set(m.SYNONYM_TABLE.values())
    for key, bucket in cleaned_votes.items():
        if key in synonym_destinations and native_votes.get(key):
            # Documented merge direction wins: vote among natives only.
            bucket = native_votes[key]
        display_of_key[key] = max(bucket.items(),
                                  key=lambda kv: (kv[1], _tie(kv[0])))[0]

    # ---- 2e. per-alias cone audit ----------------------------------------
    # For every raw name with solved coordinates, measure its distance to
    # the *final group's* pooled median position.  1 = inside the cone,
    # 0 = outside (worth an eyebrow), NULL = no solved coordinates.
    group_coords: dict[str, tuple[float, float]] = {}
    for key in set(k for k in key_of_raw.values() if k is not None):
        members = [coords_of_raw[r] for r, k in key_of_raw.items()
                   if k == key and r in coords_of_raw]
        if members:
            group_coords[key] = m.median_radec([c[0] for c in members],
                                               [c[1] for c in members])

    rows = []
    for raw, n in raw_counts.items():
        info = norm[raw]
        key = key_of_raw.get(raw)
        method = ",".join(info.rules) if info.rules else "identity"
        if "synonym" in info.rules and info.pre_synonym_key in refused_synonyms:
            method += ",synonym_refused_by_cone"
        cone: float | None = None
        if key is not None and raw in coords_of_raw and key in group_coords:
            here, there = coords_of_raw[raw], group_coords[key]
            sep = m.angular_separation_deg(here[0], here[1],
                                           there[0], there[1])
            cone = 1 if sep <= m.CONE_RADIUS_DEG else 0
        rows.append({
            "raw_name": raw if raw else None,
            "canonical_target": display_of_key.get(key) if key else None,
            "target_key": key,
            "n_frames": int(n),
            "method": method,
            "cone_check_passed": cone,
        })
    aliases = pd.DataFrame(rows)
    return aliases, key_of_raw, display_of_key


def _tie(s: str) -> tuple:
    """Deterministic tie-break helper: shorter, then lexicographic reverse
    so that the comparison inside ``resolve_aliases`` prefers the higher
    count first and stays stable for equal counts."""
    return (-len(s), s)


# ---------------------------------------------------------------------------
# Step 3 — frames table: dedup, canonical choice, eras, nights, pointing, QC
# ---------------------------------------------------------------------------
def load_prior_era_ids(db_path: Path) -> dict:
    """Read the era registry a previous manifest build published, if any.

    Returns ``{era_key_tuple: era_id}`` from the ``eras`` table of the
    manifest at ``db_path``, or ``{}`` when no prior build exists (first
    run, or the table is absent).  Keys are re-derived through
    :func:`macro_core.manifest.era_key` from the stored components, so the
    normalization (whitespace strip, int casts, EGAIN rounding) is byte-for-
    byte the same one the assignment step uses — a stored key always maps
    onto its own registry entry.
    """
    if not Path(db_path).exists():
        return {}
    try:
        with closing(sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)) as con:
            rows = con.execute(
                "SELECT era_id, readoutm, naxis1, naxis2, xbinning, egain "
                "FROM eras").fetchall()
    except sqlite3.Error:
        # No eras table (or unreadable DB) — behave as a first build.
        return {}
    return {m.era_key(r, n1, n2, xb, eg): int(eid)
            for eid, r, n1, n2, xb, eg in rows}


def build_frames(df: pd.DataFrame, key_of_raw: dict,
                 display_of_key: dict,
                 prior_era_ids: dict | None = None,
                 ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute every derived column of the ``frames`` table.

    ``prior_era_ids`` maps era keys (as built by :func:`macro_core.manifest
    .era_key`) to the era_id a PREVIOUS manifest build published for them.
    Passing it pins those ids forever — see the registry comment at the era
    assignment step for why renumbering is forbidden.

    Returns the enriched frames DataFrame and the eras table.
    """
    # ---- 3a. basename and night label ------------------------------------
    df = df.copy()
    df["basename"] = [m.basename_of(p) for p in df["path"]]
    df["night"] = [m.night_label(j) for j in df["jd"]]

    # ---- 3b. alias key + canonical target for every row ------------------
    raw_series = df["target_best"].fillna("")
    df["target_key"] = raw_series.map(key_of_raw)
    df["canonical_target"] = df["target_key"].map(display_of_key)

    # ---- 3c. global duplicate groups on (basename, jd) -------------------
    # Frames without a JD (unreadable headers) become singleton groups keyed
    # by their catalog rowid — dup_key() implements that rule.
    keys = [m.dup_key(bn, jd, rid) for bn, jd, rid
            in zip(df["basename"], df["jd"], df["obs_rowid"])]
    df["dup_group"] = pd.factorize(pd.Series(keys, dtype=object))[0]

    # ---- 3d. canonical member per duplicate group (tree policy) ----------
    # Rank every row by its effective tree priority: the archive-wide
    # default, or a documented per-target exception (NGC 5548 → macalester).
    # The rank of every tree comes from m.tree_rank() — the SAME function
    # choose_canonical() uses and the unit tests exercise — evaluated once
    # per distinct tree and broadcast with a plain dict .map, so the
    # vectorized selection here cannot drift from the tested pure logic
    # (test_manifest.py asserts the two agree member-by-member).
    trees_seen = df["tree"].unique()
    default_rank = {t: m.tree_rank(t) for t in trees_seen}
    # NOTE: .map(dict.get) — a per-element callable lookup — is used instead
    # of .map(dict) throughout this script: pandas converts a plain dict to
    # a Series first, and that round-trip corrupts lookups for keys pandas
    # cannot index cleanly (the era-key tuples containing None were silently
    # dropped that way — the shipped era_id bug).  dict.get never converts.
    rank = df["tree"].map(default_rank.get).astype(int)
    for key, prio in m.TREE_PRIORITY_EXCEPTIONS.items():
        mask = df["target_key"] == key
        if mask.any():
            exc_rank = {t: m.tree_rank(t, prio) for t in trees_seen}
            rank.loc[mask] = df.loc[mask, "tree"].map(exc_rank.get).astype(int)
    df["_rank"] = rank
    # Within each group: best rank wins; ties broken by lexicographically
    # smallest path (earliest night directory — the SN July copies lose).
    order = df.sort_values(["dup_group", "_rank", "path"],
                           kind="mergesort")  # stable sort → deterministic
    winners = order.groupby("dup_group", sort=False).head(1).index
    df["is_canonical"] = 0
    df.loc[winners, "is_canonical"] = 1
    df.drop(columns=["_rank"], inplace=True)

    # ---- 3e. era assignment (READOUTM, geometry, binning, EGAIN) ---------
    # Unreadable rows (header error) carry no camera keys → era NULL.
    ok = df["error"].isna()
    ekeys = [m.era_key(r, n1, n2, xb, eg) if o else None
             for o, r, n1, n2, xb, eg
             in zip(ok, df["readoutm"], df["naxis1"], df["naxis2"],
                    df["xbinning"], df["egain"])]
    df["_era_key"] = pd.Series(ekeys, index=df.index, dtype=object)
    # Era ids are a REGISTRY, not a ranking (2026-08-18).  The first build
    # ordered ids by each configuration's first appearance on sky, and those
    # numbers are now published: reports, all five strategy documents, the
    # ops request, and the s1_*/detector_params/phot_* tables all cite
    # "era 76", "era 47", ....  Renumbering would silently re-point every
    # one of those references at a different camera.  So: ids already issued
    # by a previous build are PINNED via ``prior_era_ids``; only keys new to
    # this build receive fresh ids, appended AFTER the existing maximum, in
    # first-appearance (min JD) order among themselves.  The original
    # oldest-camera-is-era-1 property therefore holds only within the first
    # build's keys — a documented, deliberate trade for citation stability.
    # (Trigger: the Calibrations/ recovery added mid-timeline configurations
    # that would have spliced into a pure JD ordering and renumbered every
    # era after them.)
    firsts = (df[ok].groupby("_era_key")["jd"].min().sort_values())
    era_id_of_key = dict(prior_era_ids or {})
    next_id = max(era_id_of_key.values(), default=0) + 1
    for k in firsts.index:
        if k not in era_id_of_key:
            era_id_of_key[k] = next_id
            next_id += 1
    # REGRESSION-CRITICAL: look keys up with dict.get, never .map(dict).
    # pandas turns a tuple-keyed dict into a MultiIndex-backed Series, and
    # any None inside a key tuple becomes a NaN level whose lookup MISSES —
    # in the first shipped build every missing-EGAIN era (29 of 83, 5,779
    # frames) silently got era_id NULL this way.  dict.get is a per-element
    # hash lookup and cannot suffer index conversion.
    df["era_id"] = df["_era_key"].map(era_id_of_key.get).astype("Int64")
    # Contract check, enforced at build time: every error-free frame MUST
    # carry an era_id (only the handful of unreadable-header rows may not).
    n_unassigned = int((df["era_id"].isna() & ok).sum())
    if n_unassigned:
        raise AssertionError(
            f"era assignment incomplete: {n_unassigned} error-free frames "
            "received no era_id — the key lookup regressed")

    # Era summary table: canonical error-free frames per era + night span.
    canon = df[(df["is_canonical"] == 1) & ok]
    era_rows = []
    for ekey, era_id in era_id_of_key.items():
        sub = canon[canon["_era_key"] == ekey]
        if len(sub) == 0:
            # Every copy of this configuration's frames lost dedup (rare);
            # fall back to all rows so the era is still documented.
            sub = df[df["_era_key"] == ekey]
        nights = sub["night"].dropna()
        era_rows.append({
            "era_id": era_id,
            "readoutm": ekey[0], "naxis1": ekey[1], "naxis2": ekey[2],
            "xbinning": ekey[3], "egain": ekey[4],
            "n_frames": int(len(sub)),
            "first_night": nights.min() if len(nights) else None,
            "last_night": nights.max() if len(nights) else None,
        })
    eras = pd.DataFrame(era_rows).sort_values("era_id").reset_index(drop=True)
    df.drop(columns=["_era_key"], inplace=True)

    # ---- 3f. pointing validation -----------------------------------------
    # Reference position per target: median (ra_deg, dec_deg) over its
    # plate-solved canonical frames, only when at least
    # MIN_SOLVED_FOR_REFERENCE of them exist (spec).  RA-wrap-aware median.
    refs: dict[str, tuple[float, float]] = {}
    solved_canon = df[(df["is_canonical"] == 1) & (df["pltsolvd"] == 1)
                      & df["ra_deg"].notna() & df["dec_deg"].notna()
                      & df["target_key"].notna()]
    for key, grp in solved_canon.groupby("target_key"):
        if len(grp) >= m.MIN_SOLVED_FOR_REFERENCE:
            refs[key] = m.median_radec(grp["ra_deg"], grp["dec_deg"])

    # Offset for EVERY frame that has coordinates and a target reference —
    # canonical and duplicate alike (a duplicate's pointing is still real).
    offsets = np.full(len(df), np.nan)
    has = df["ra_deg"].notna() & df["dec_deg"].notna() \
        & df["target_key"].isin(refs.keys())
    idx = np.flatnonzero(has.to_numpy())
    ra = df["ra_deg"].to_numpy()
    dec = df["dec_deg"].to_numpy()
    tkey = df["target_key"].to_numpy(dtype=object)
    for i in idx:
        ra0, dec0 = refs[tkey[i]]
        offsets[i] = m.angular_separation_deg(ra[i], dec[i], ra0, dec0)
    df["pointing_offset_deg"] = offsets

    # ---- 3g. QC flags ----------------------------------------------------
    df["qc_flags"] = [
        m.qc_flags(err, ext, am, jd, radeg, tk, off)
        for err, ext, am, jd, radeg, tk, off
        in zip(df["error"], df["exptime"], df["airmass"], df["jd"],
               df["ra_deg"], df["target_key"], df["pointing_offset_deg"])
    ]
    return df, eras


# ---------------------------------------------------------------------------
# Step 4 — per-project reconciliation counts (Output C)
# ---------------------------------------------------------------------------
def build_project_counts(df: pd.DataFrame,
                         display_of_key: dict) -> pd.DataFrame:
    """Count canonical frames per strategy claim, like-for-like.

    Each claim in ``manifest.STRATEGY_CLAIMS`` names a counting metric; this
    function implements those metrics against the manifest and records the
    manifest value next to the claimed value.  Section 7 of the report
    renders the diff — the whole point of S0's existence.
    """
    is_light = df["imagetyp"].fillna("").str.startswith("Light")
    err_free = df["error"].isna()
    canonical = df["is_canonical"] == 1

    rows = []
    for (project, tkey, metric,
         claimed_frames, claimed_nights, source) in m.STRATEGY_CLAIMS:
        # Target selector: one alias key, or the whole Dw survey family.
        if tkey == "__dw_survey__":
            tsel = df["target_key"].fillna("").str.startswith("dw1")
            display = "Dw survey (19 fields)"
            in_primary = df["tree"] == "rawimage"
        else:
            tsel = df["target_key"] == tkey
            display = display_of_key.get(tkey, tkey)
            # The tree this target's own strategy counted (rawimage, or the
            # documented exception — NGC 5548's macalester superset).
            in_primary = df["tree"] == m.primary_tree(tkey)

        # Metric selector — every rule documented in manifest.py.
        if metric == "rows_all_trees":
            sel = tsel & is_light                       # raw rows, no dedup
        elif metric == "unique_light":
            sel = tsel & is_light & err_free & canonical
        elif metric == "grism_light":
            sel = tsel & is_light & err_free & canonical \
                & df["filter"].isin(m.GRISM_FILTERS)
        elif metric == "grism4_light":
            sel = tsel & is_light & err_free & canonical \
                & df["filter"].isin(m.GRISM4_FILTERS)
        else:                                           # defensive: typo in the table
            raise ValueError(f"unknown metric {metric!r}")

        # Two views of the same selection: like-for-like with the claim
        # (primary tree only — what the strategy's rule counted), and the
        # fully global canonical count (nothing hidden).  rows_all_trees is
        # by definition global, so both views coincide there.
        if metric == "rows_all_trees":
            like = sel
        else:
            like = sel & in_primary
        n_frames = int(like.sum())
        n_nights = int(df.loc[like, "night"].dropna().nunique())
        n_frames_global = int(sel.sum())
        n_nights_global = int(df.loc[sel, "night"].dropna().nunique())
        rows.append({
            "project": project,
            "target": display,
            "target_key": tkey,
            "metric": metric,
            "claimed_frames": claimed_frames,
            "claimed_nights": claimed_nights,
            "manifest_frames": n_frames,
            "manifest_nights": n_nights,
            "manifest_frames_global": n_frames_global,
            "manifest_nights_global": n_nights_global,
            "diff_frames": (n_frames - claimed_frames
                            if claimed_frames is not None else None),
            "diff_nights": (n_nights - claimed_nights
                            if claimed_nights is not None else None),
            "source": source,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Step 5 — atomic write
# ---------------------------------------------------------------------------
def write_manifest(out_path: Path, frames: pd.DataFrame, aliases: pd.DataFrame,
                   eras: pd.DataFrame, project_counts: pd.DataFrame,
                   catalog_path: Path) -> None:
    """Write all five tables to a temp file, then atomically swap it in.

    ``os.replace`` on the same filesystem is atomic, so a reader never sees
    a half-built manifest and an interrupted build changes nothing.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # Temp file must live in the SAME directory for os.replace to be atomic.
    fd, tmp_name = tempfile.mkstemp(prefix="rlmt-manifest.", suffix=".tmp",
                                    dir=out_path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        # closing() is required: sqlite3's own context manager only manages
        # the TRANSACTION (commit/rollback) — it does NOT close the file
        # handle, and os.replace over a still-open handle is a portability
        # hazard (fails on Windows).  closing() guarantees the connection is
        # closed BEFORE the swap below.
        with closing(sqlite3.connect(tmp)) as con, con:
            frames.to_sql("frames", con, index=False)
            aliases.to_sql("aliases", con, index=False)
            eras.to_sql("eras", con, index=False)
            project_counts.to_sql("project_counts", con, index=False)
            # Build provenance: enough to reproduce or audit this file.
            meta = pd.DataFrame([
                {"key": "built_utc",
                 "value": datetime.now(timezone.utc).isoformat()},
                {"key": "catalog_path", "value": str(catalog_path)},
                {"key": "code_version", "value": S0_CODE_VERSION},
                {"key": "git_commit", "value": _git_commit()},
                {"key": "night_shift_days", "value": str(m.NIGHT_SHIFT_DAYS)},
                {"key": "cone_radius_deg", "value": str(m.CONE_RADIUS_DEG)},
                {"key": "pointing_outlier_deg",
                 "value": str(m.POINTING_OUTLIER_DEG)},
            ])
            meta.to_sql("build_meta", con, index=False)
            # Indexes for the report's (and downstream stages') queries.
            cur = con.cursor()
            cur.execute("CREATE INDEX ix_frames_dup ON frames(dup_group)")
            cur.execute("CREATE INDEX ix_frames_tgt ON frames(canonical_target)")
            cur.execute("CREATE INDEX ix_frames_key ON frames(target_key)")
            cur.execute("CREATE INDEX ix_frames_night ON frames(night)")
            cur.execute("CREATE INDEX ix_frames_era ON frames(era_id)")
            con.commit()
        # mkstemp creates 0600 files; open the permissions up BEFORE the
        # swap so the live path never exists in a locked-down state.
        os.chmod(tmp, 0o644)
        os.replace(tmp, out_path)          # the atomic swap
    finally:
        if tmp.exists():                   # only on failure paths
            tmp.unlink()


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
            "Build the S0 manifest database (global dedup, alias table, "
            "camera eras, night labels, pointing validation, QC flags) from "
            "the RLMT observation catalog, then render the chain-of-evidence "
            "report. Safe to re-run: the manifest is rebuilt from scratch "
            "and swapped in atomically."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG,
                   help="observation catalog (opened READ-ONLY)")
    p.add_argument("--out", type=Path, default=DEFAULT_MANIFEST,
                   help="manifest database to (re)build")
    p.add_argument("--eras-csv", type=Path, default=DEFAULT_ERAS_CSV,
                   help="where to write the era table as CSV")
    p.add_argument("--skip-report", action="store_true",
                   help="build the database only; do not render the HTML "
                        "report/figures afterwards")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.catalog.exists():
        print(f"ERROR: catalog not found: {args.catalog}", file=sys.stderr)
        return 2

    print(f"[S0] reading catalog {args.catalog} ...")
    df = load_catalog(args.catalog)
    print(f"[S0]   {len(df):,} catalog rows")

    print("[S0] resolving target aliases ...")
    aliases, key_of_raw, display_of_key = resolve_aliases(df)
    n_raw = aliases["raw_name"].notna().sum()
    n_canon = aliases["target_key"].nunique()
    print(f"[S0]   {n_raw:,} raw names -> {n_canon:,} canonical targets")

    print("[S0] building frames table (dedup, eras, nights, pointing, QC) ...")
    # Pin era ids already published by the previous build (registry rule —
    # see the era-assignment comment in build_frames).  Read BEFORE the
    # atomic swap replaces the file.
    prior_era_ids = load_prior_era_ids(args.out)
    if prior_era_ids:
        print(f"[S0] era registry: pinning {len(prior_era_ids)} ids "
              "from the previous build")
    frames, eras = build_frames(df, key_of_raw, display_of_key,
                                prior_era_ids=prior_era_ids)
    n_groups = frames["dup_group"].nunique()
    n_can = int((frames["is_canonical"] == 1).sum())
    print(f"[S0]   {n_groups:,} duplicate groups; {n_can:,} canonical frames; "
          f"{len(eras)} camera eras")

    print("[S0] computing project reconciliation counts ...")
    counts = build_project_counts(frames, display_of_key)

    print(f"[S0] writing manifest -> {args.out}")
    write_manifest(args.out, frames, aliases, eras, counts, args.catalog)
    args.eras_csv.parent.mkdir(parents=True, exist_ok=True)
    eras.to_csv(args.eras_csv, index=False)
    print(f"[S0] wrote era table CSV -> {args.eras_csv}")

    if not args.skip_report:
        print("[S0] rendering chain-of-evidence report ...")
        from macro_core import report_s0
        report_path = report_s0.render_report(args.out)
        print(f"[S0] report -> {report_path}")

    print("[S0] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
