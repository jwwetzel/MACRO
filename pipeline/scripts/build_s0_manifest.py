#!/usr/bin/env python
"""Build the S0 manifest database for the MACRO/RLMT archive.

WHAT THIS SCRIPT DOES (stage S0 of the shared pipeline, ROADMAP.md sec. 1.1)
---------------------------------------------------------------------------
Reads the observation catalog (READ-ONLY — this script never writes to it),
applies the pure S0 logic from ``macro_core.manifest`` to every row, and
writes a fresh manifest database with five tables:

* ``frames``          — one row per catalog row, plus: basename, night label,
                        era_id, duplicate group id, canonical flag and the
                        evidence for it (``dup_basis``), resolved target
                        name, pointing offset, QC flags, and the typed
                        hardware-state columns from the header re-scrape
                        (camera, gain/offset setting, set-point, cooler
                        power, true focuser position, flip state, pier side,
                        wheel slot and wheel map).
* ``aliases``         — every raw target name → its canonical name, with the
                        exact normalization rules that fired and the result
                        of the coordinate-cone audit.
* ``eras``            — the camera-era table keyed on (READOUTM, geometry,
                        binning, EGAIN), also exported as CSV.  An era whose
                        frames are all reduced copies of another era's
                        frames is marked as that era's ALIAS.
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
The five S0 tables are rebuilt from scratch on every run and swapped into
the live manifest inside ONE SQLite transaction (new tables are written
under temporary names, then ``DROP`` old + ``RENAME`` new + ``COMMIT``).  A
reader therefore sees either the complete old manifest or the complete new
one, an interrupted build changes nothing, every other stage's tables are
left exactly where they are, and a sibling stage writing its own table at
the same moment is serialized by SQLite instead of losing its write.

That last point is why the swap is a transaction and no longer a file
replace (changed 2026-10-03).  The old mechanism built a new FILE and
``os.replace``d it over the live one.  It needed a carry step to keep the
sibling tables, it left the previous file's ``-wal`` beside the new file
(the 2026-08-19 "malformed database schema" incident), and — the defect
that had not bitten yet — any process holding the old file open kept
writing to an inode that no longer had a name.  The file-replace path
survives as ``--fresh-file`` for building a manifest where none exists or
for a deliberate clean-file rebuild with nobody else attached.

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
    with closing(sqlite3.connect(uri, uri=True, timeout=300.0)) as con:
        # rowid gives every catalog row a stable integer identity that we
        # carry into the manifest (useful for tracing a frame back).
        df = pd.read_sql_query("SELECT rowid AS obs_rowid, * FROM obs", con)
        # The hardware-card re-scrape (F-2) lives in its own catalog table,
        # written by ``rescan_geometry.py hdr-run``.  It is joined here, by
        # path, as ``h_*`` text columns; ``h_scanned`` marks the rows whose
        # header was actually read.  A catalog with no such table (a first
        # build, a test fixture) simply yields frames with NULL hardware
        # columns and ``hdr_scanned = 0`` — the absence is visible, never
        # papered over.
        has_hdr = con.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' "
            "AND name='hdr_rescrape'").fetchone() is not None
        if has_hdr:
            hdr = pd.read_sql_query(
                "SELECT * FROM hdr_rescrape WHERE error IS NULL", con)
            hdr = hdr.drop(columns=["n_cards", "error", "scanned_utc"])
            hdr["h_scanned"] = 1
            df = df.merge(hdr, on="path", how="left", validate="one_to_one")
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


#: An era is another era's ALIAS when at least this fraction of its rows
#: are non-canonical copies of frames in that one other era.  Not 1.0: a
#: reduced tree always holds a few products with no raw parent (stacks,
#: frames whose JD was rewritten during reduction), and those must not stop
#: a 25,000-row alias from being named.  Not much lower either — below
#: this the era holds real exposures of its own.
ERA_ALIAS_MIN_FRACTION = 0.95

#: ``frames.hdr_scanned`` values.
HDR_NOT_SCANNED = 0      # no header re-scrape for this row: hardware NULL
HDR_SCANNED = 1          # this file's own header was read
HDR_INHERITED = 2        # values copied from an EXACT copy of the file


def attach_hardware_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add the typed hardware-state columns to the frames DataFrame (F-2).

    Input: the catalog rows, carrying the re-scrape's ``h_*`` text columns
    when ``load_catalog`` found the ``hdr_rescrape`` table (and none of them
    otherwise).  Output: the same rows with the columns listed in
    ``manifest.HARDWARE_FRAME_COLUMNS`` and WITHOUT the raw ``h_*`` columns.

    Three things happen, in this order:

    1.  **Backfill.**  ``instrume``, ``swcreate`` and ``ccd_temp`` exist in
        the original scan but are NULL wherever astropy gave up on the
        header.  Where the raw-card re-scrape read them, they are filled —
        the acceptance line "no nulls where the card exists".  A value the
        original scan DID hold is never overwritten.
    2.  **Typing.**  ``manifest.hardware_columns`` turns the header text
        into the typed columns.
    3.  **Inheritance for exact copies.**  A row that was not re-scraped
        but is an exact copy (same basename, same JD — ``dup_basis =
        'same_basename_jd'``) of a canonical row that WAS takes that row's
        values and is marked ``hdr_scanned = 2``.  A byte copy of a file
        has the file's header; a reduced DERIVATIVE does not, and never
        inherits.  With a completed re-scrape this step changes nothing —
        it exists so a partial scan degrades to labelled inheritance
        instead of to silent NULLs.
    """
    hcols = [c for c in df.columns if c.startswith("h_")]
    work = df[hcols].astype(object).where(df[hcols].notna(), None) \
        if hcols else pd.DataFrame(index=df.index)

    # ---- 1. backfill the three cards the original scan also read ---------
    if hcols:
        for col, hcol in (("instrume", "h_instrume"),
                          ("swcreate", "h_swcreate")):
            if col in df and hcol in work:
                fill = work[hcol].map(lambda v: v if v else None)
                df[col] = df[col].where(df[col].notna(), fill)
        if "ccd_temp" in df and "h_ccd_temp" in work:
            fill = work["h_ccd_temp"].map(m.card_float)
            df["ccd_temp"] = df["ccd_temp"].where(df["ccd_temp"].notna(),
                                                  fill)

    # ---- 2. typing --------------------------------------------------------
    records = work.to_dict("records") if hcols else [{}] * len(df)
    instr = df["instrume"] if "instrume" in df else [None] * len(df)
    typed = []
    for rec, ins in zip(records, instr):
        rec = dict(rec)
        # Camera identity uses the BACKFILLED instrume, so a row the
        # original scan read but the re-scrape has not reached yet still
        # gets its camera.
        rec["h_instrume"] = ins if isinstance(ins, str) else None
        typed.append(m.hardware_columns(rec))
    hw = pd.DataFrame(typed, index=df.index,
                      columns=[c for c, _ in m.HARDWARE_FRAME_COLUMNS])
    df = pd.concat([df.drop(columns=hcols), hw], axis=1)

    # ---- 3. exact copies inherit from their scanned canonical row --------
    need = (df["hdr_scanned"] == HDR_NOT_SCANNED) \
        & (df["dup_basis"] == m.DUP_SAME_BASENAME)
    if need.any():
        heads = df[(df["is_canonical"] == 1)
                   & (df["hdr_scanned"] == HDR_SCANNED)].set_index("dup_group")
        cols = [c for c, _ in m.HARDWARE_FRAME_COLUMNS if c != "hdr_scanned"]
        got = df.loc[need, "dup_group"].isin(heads.index)
        idx = got[got].index
        if len(idx):
            src = heads.loc[df.loc[idx, "dup_group"], cols]
            df.loc[idx, cols] = src.to_numpy()
            df.loc[idx, "hdr_scanned"] = HDR_INHERITED
    df["fwpos"] = df["fwpos"].astype("Int64")
    df["hdr_scanned"] = df["hdr_scanned"].astype(int)
    return df


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

    # ---- 3c. global duplicate groups on (exposure root, jd) --------------
    # F-1: a reduced or re-packaged copy of a frame is the same exposure.
    # Every basename is taken apart by the pure m.exposure_name(); a file
    # carrying a processing suffix is then resolved to the raw frame it
    # derives from by m.exposure_root(), which needs to know whether an
    # UNPROCESSED file of a given stem exists at the SAME header JD.  That
    # question is answered from one set built here — the same lookup the
    # unit tests drive with a lambda.
    names = [m.exposure_name(bn) for bn in df["basename"]]
    jds = [None if (j is None or (isinstance(j, float) and np.isnan(j)))
           else float(j) for j in df["jd"]]
    raw_at = {(n.stem, j) for n, j in zip(names, jds)
              if j is not None and n.n_processing == 0}
    roots, hows = [], []
    for n, j in zip(names, jds):
        if j is None:
            # No JD: nothing can be proven about this file (dup_key makes it
            # a singleton whatever its name), so no parent is looked for.
            roots.append(n.stem)
            hows.append(m.ROOT_IS_SELF)
            continue
        root, how = m.exposure_root(n, lambda stem, j=j: (stem, j) in raw_at)
        roots.append(root)
        hows.append(how)
    df["_how"] = hows
    df["_nproc"] = [n.n_processing for n in names]
    # Frames without a JD (unreadable headers) become singleton groups keyed
    # by their catalog rowid — dup_key() implements that rule.
    keys = [m.dup_key(root, jd, rid) for root, jd, rid
            in zip(roots, df["jd"], df["obs_rowid"])]
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
    # Within each group: the least-processed copy wins (a raw frame beats
    # its own ``_calibrated`` twin whatever tree either sits in), then best
    # tree rank, then lexicographically smallest path (earliest night
    # directory — the SN July copies lose).  Same key, same order, as the
    # pure m.choose_canonical().
    order = df.sort_values(["dup_group", "_nproc", "_rank", "path"],
                           kind="mergesort")  # stable sort → deterministic
    winners = order.groupby("dup_group", sort=False).head(1).index
    df["is_canonical"] = 0
    df.loc[winners, "is_canonical"] = 1
    # Why each non-canonical row is a duplicate — the audit trail of F-1.
    head_name = dict(zip(df.loc[winners, "dup_group"],
                         df.loc[winners, "basename"]))
    head_how = dict(zip(df.loc[winners, "dup_group"],
                        df.loc[winners, "_how"]))
    df["dup_basis"] = [
        m.dup_basis(bn, head_name[g], canon == 1, how, head_how[g])
        for bn, g, canon, how in zip(df["basename"], df["dup_group"],
                                     df["is_canonical"], df["_how"])]
    df.drop(columns=["_rank", "_nproc", "_how"], inplace=True)

    # ---- 3d'. typed hardware-state columns (F-2) --------------------------
    df = attach_hardware_columns(df)

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

    # Era ALIASES (F-1).  The QHY600's reduced files are trimmed of their
    # overscan border (4787x3193 against the raw 4800x3211), and geometry is
    # part of the era key, so the reduced copies of one camera's frames
    # formed "eras" of their own: 79 beside 78, 82 beside 81.  With the
    # copies now deduplicated against their raw parents, such an era holds
    # no exposure of its own.  For every era we therefore record how many of
    # its rows are non-canonical copies whose canonical frame lives in ONE
    # other era; at or above ERA_ALIAS_MIN_FRACTION the era is that era's
    # alias.  The id stays in the registry (pinned, never reused) — only its
    # meaning is now stated instead of implied.
    # Plain float arrays with NaN for "no era" (header-error rows): the
    # nullable Int64 column raises on NA comparisons, and NaN != NaN would
    # otherwise count two era-less rows as "foreign" to each other.
    own_era = df["era_id"].astype("float64").to_numpy()
    is_head = (df["is_canonical"] == 1).to_numpy()
    head_era = dict(zip(df["dup_group"].to_numpy()[is_head],
                        own_era[is_head]))
    parent_era = np.array([head_era[g] for g in df["dup_group"]],
                          dtype="float64")
    n_rows_of = df[ok].groupby("era_id").size().to_dict()
    n_canon_of = df[ok & (df["is_canonical"] == 1)].groupby(
        "era_id").size().to_dict()
    era_alias: dict[int, tuple] = {}
    foreign = (ok.to_numpy() & ~is_head & ~np.isnan(parent_era)
               & ~np.isnan(own_era) & (parent_era != own_era))
    for era_id in np.unique(own_era[foreign]):
        parents, counts = np.unique(parent_era[foreign & (own_era == era_id)],
                                    return_counts=True)
        # Modal foreign parent era; ties break to the smaller id (argmax
        # returns the first maximum of the id-sorted unique values).
        best = int(np.argmax(counts))
        frac = float(counts[best]) / float(n_rows_of[int(era_id)])
        if frac >= ERA_ALIAS_MIN_FRACTION:
            era_alias[int(era_id)] = (int(parents[best]), round(frac, 6))

    # Era summary table: canonical error-free frames per era + night span.
    canon = df[(df["is_canonical"] == 1) & ok]
    era_rows = []
    for ekey, era_id in era_id_of_key.items():
        sub = canon[canon["_era_key"] == ekey]
        # n_frames is the CANONICAL count, always — including zero.  Until
        # 2026-10-03 an era whose every frame lost dedup fell back to its
        # raw row count "so the era is still documented"; with F-1 that
        # made the two reduced-alias eras claim 25,715 and 1,682 frames
        # they do not own (and tripped the S0 report's own guard, which
        # requires every era above the figure threshold to hold canonical
        # frames).  The row count now has its own column, n_rows, and only
        # the NIGHT SPAN falls back to all rows, so the era stays dated.
        n_frames = int(len(sub))
        if len(sub) == 0:
            sub = df[df["_era_key"] == ekey]
        nights = sub["night"].dropna()
        alias_of, alias_frac = era_alias.get(era_id, (None, None))
        era_rows.append({
            "era_id": era_id,
            "readoutm": ekey[0], "naxis1": ekey[1], "naxis2": ekey[2],
            "xbinning": ekey[3], "egain": ekey[4],
            "n_frames": n_frames,
            "first_night": nights.min() if len(nights) else None,
            "last_night": nights.max() if len(nights) else None,
            "n_rows": int(n_rows_of.get(era_id, 0)),
            "n_canonical": int(n_canon_of.get(era_id, 0)),
            "alias_of_era": alias_of,
            "alias_fraction": alias_frac,
        })
    eras = pd.DataFrame(era_rows).sort_values("era_id").reset_index(drop=True)
    # Nullable integer: an era that is nobody's alias holds NULL, not NaN-
    # coerced-to-float (which would write 78.0 into the registry).
    eras["alias_of_era"] = eras["alias_of_era"].astype("Int64")
    df.drop(columns=["_era_key"], inplace=True)

    # ---- 3f. pointing validation -----------------------------------------
    # Reference position per target: plate-solved canonical frames when at
    # least MIN_SOLVED_FOR_REFERENCE of them exist, otherwise the median of
    # that target's HEADER coordinates.  The fallback exists because NGC 5548
    # has zero plate-solved frames out of 279 and therefore had no reference,
    # no offsets, and no pointing flags at all — see manifest.pointing_
    # reference for the full argument.  The basis is recorded per frame so
    # the weaker evidence is visible instead of implicit.
    refs: dict[str, tuple[float, float]] = {}
    basis_of_key: dict[str, str] = {}
    canon = df[(df["is_canonical"] == 1) & df["target_key"].notna()]
    for key, grp in canon.groupby("target_key"):
        ra0, dec0, basis = m.pointing_reference(
            grp["ra_deg"], grp["dec_deg"], grp["pltsolvd"])
        basis_of_key[key] = basis
        if ra0 is not None:
            refs[key] = (ra0, dec0)

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
    # One column, one meaning: how THIS frame's target reference was derived
    # ('plate_solved', 'header_median', or 'none' for an unreferenced target).
    df["pointing_ref_basis"] = [
        basis_of_key.get(k, m.POINTING_REF_NONE) for k in tkey]

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
#: Suffix of the temporary tables the in-place swap writes before renaming.
_TMP_SUFFIX = "_s0_tmp"

#: Indexes on ``frames`` for the report's (and downstream stages') queries.
#: Created INSIDE the swap transaction, so no reader ever meets an
#: un-indexed frames table.
_FRAME_INDEXES = (
    ("ix_frames_dup", "dup_group"),
    ("ix_frames_tgt", "canonical_target"),
    ("ix_frames_key", "target_key"),
    ("ix_frames_night", "night"),
    ("ix_frames_era", "era_id"),
    ("ix_frames_rowid", "obs_rowid"),
)


def _build_meta(catalog_path: Path) -> pd.DataFrame:
    """Build provenance: enough to reproduce or audit this build."""
    return pd.DataFrame([
        {"key": "built_utc",
         "value": datetime.now(timezone.utc).isoformat()},
        {"key": "catalog_path", "value": str(catalog_path)},
        {"key": "code_version", "value": S0_CODE_VERSION},
        {"key": "git_commit", "value": _git_commit()},
        {"key": "night_shift_days", "value": str(m.NIGHT_SHIFT_DAYS)},
        {"key": "cone_radius_deg", "value": str(m.CONE_RADIUS_DEG)},
        {"key": "pointing_outlier_deg",
         "value": str(m.POINTING_OUTLIER_DEG)},
        {"key": "era_alias_min_fraction",
         "value": str(ERA_ALIAS_MIN_FRACTION)},
    ])


def write_manifest(out_path: Path, frames: pd.DataFrame, aliases: pd.DataFrame,
                   eras: pd.DataFrame, project_counts: pd.DataFrame,
                   catalog_path: Path, fresh_file: bool = False) -> None:
    """Write the five S0 tables into the manifest, atomically.

    Two mechanisms, chosen by whether a manifest already exists:

    * **in place** (the default whenever ``out_path`` exists) — see
      :func:`_swap_in_place`: one transaction inside the live database.
    * **fresh file** (``out_path`` absent, or ``fresh_file=True``) — see
      :func:`_write_fresh_file`: build a new file, carry the sibling tables
      across, ``os.replace``.  Only safe when no other process has the
      manifest open; the CLI flag says so.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tables = {"frames": frames, "aliases": aliases, "eras": eras,
              "project_counts": project_counts,
              "build_meta": _build_meta(catalog_path)}
    if out_path.exists() and not fresh_file:
        _swap_in_place(out_path, tables)
    else:
        _write_fresh_file(out_path, tables)


def _swap_in_place(out_path: Path, tables: dict) -> None:
    """Replace the S0 tables inside the live manifest in ONE transaction.

    WHY NOT A FILE REPLACE.  The manifest is shared: S1, S2, S3, the grism
    track and the plan ledger all keep tables in it and several of them are
    open at any moment (the file runs in WAL mode).  Replacing the FILE
    under them has three failure modes, two of which this project has
    already paid for:

    1.  every sibling table has to be copied across or it is destroyed
        (2026-08-18: S1's and the detector memo's evidence tables);
    2.  the old file's ``-wal`` survives beside the new file and SQLite
        replays it against a database it does not describe (2026-08-19);
    3.  a process that still holds the OLD file open goes on writing to an
        inode with no name — its work vanishes without an error.

    A transaction has none of them.  The new tables are written first under
    temporary names (slow, but invisible to every consumer); then a single
    ``BEGIN IMMEDIATE … COMMIT`` drops the old five, renames the new five
    and rebuilds the indexes.  Readers see the old manifest until the
    commit and the new one after it — never a mixture, never a missing
    table — and a concurrent writer simply waits its turn.
    """
    # ---- phase 0: who else lives here ------------------------------------
    with closing(sqlite3.connect(out_path, timeout=600.0)) as con:
        con.execute("PRAGMA busy_timeout = 600000")
        siblings = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table' "
            "AND name NOT LIKE 'sqlite_%'")
            if r[0] not in S0_OWNED_TABLES
            and not r[0].endswith(_TMP_SUFFIX)]
        meta = tables["build_meta"]
        tables = dict(tables)
        # Same key as the file-replace path wrote, so every reader of
        # build_meta keeps working: the sibling tables that outlived this
        # rebuild.  Their rows are keyed to the PREVIOUS frame universe —
        # the staleness is each stage's own provenance record's to report.
        tables["build_meta"] = pd.concat([meta, pd.DataFrame([
            {"key": "carried_tables", "value": ",".join(siblings)},
            {"key": "swap", "value": "in_place_transaction"}])],
            ignore_index=True)
        # ---- phase 1: write the new tables under temporary names ----------
        # pandas commits each table itself; a crash here leaves only
        # ``*_s0_tmp`` tables behind, which the next build drops first.
        for name, frame in tables.items():
            con.execute(f'DROP TABLE IF EXISTS "{name}{_TMP_SUFFIX}"')
            con.commit()
            frame.to_sql(f"{name}{_TMP_SUFFIX}", con, index=False,
                         chunksize=20000)
    # isolation_level=None: explicit transactions only, so the BEGIN
    # IMMEDIATE below is the one and only transaction on this connection.
    with closing(sqlite3.connect(out_path, timeout=600.0,
                                 isolation_level=None)) as con:
        con.execute("PRAGMA busy_timeout = 600000")
        # ---- phase 2: the swap --------------------------------------------
        con.execute("BEGIN IMMEDIATE")
        try:
            for name in tables:
                con.execute(f'DROP TABLE IF EXISTS "{name}"')
                con.execute(f'ALTER TABLE "{name}{_TMP_SUFFIX}" '
                            f'RENAME TO "{name}"')
            for ix, col in _FRAME_INDEXES:
                con.execute(f"CREATE INDEX {ix} ON frames({col})")
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise
    if siblings:
        print(f"[S0] swapped in place; {len(siblings)} sibling table(s) "
              "untouched")
        print("[S0]   NOTE: sibling rows are keyed to the PREVIOUS frame "
              "universe — re-run those stages if frames changed.")


def _write_fresh_file(out_path: Path, tables: dict) -> None:
    """Build the manifest as a NEW file and ``os.replace`` it into place.

    The original S0 mechanism, kept for the two cases it is right for: the
    first build (there is nothing to swap into) and an explicit
    ``--fresh-file`` rebuild.  It carries every sibling table across (see
    :func:`carry_sibling_tables`) and removes the superseded file's WAL
    sidecars.  DO NOT use it while another process has the manifest open.
    """
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
            for name, frame in tables.items():
                frame.to_sql(name, con, index=False)
            cur = con.cursor()
            for ix, col in _FRAME_INDEXES:
                cur.execute(f"CREATE INDEX {ix} ON frames({col})")
            con.commit()
        # Carry downstream stages' tables forward in a SEPARATE connection:
        # SQLite forbids ATTACH inside an open transaction, and the block
        # above ran under one (sqlite3's connection context manager).
        with closing(sqlite3.connect(tmp)) as con:
            carried = carry_sibling_tables(con, out_path)
            con.execute("INSERT INTO build_meta VALUES (?, ?)",
                        ("carried_tables", ",".join(carried)))
            con.execute("INSERT INTO build_meta VALUES (?, ?)",
                        ("swap", "fresh_file_replace"))
            con.commit()
        if carried:
            print(f"[S0] carried {len(carried)} downstream table(s) forward: "
                  f"{', '.join(carried)}")
            print("[S0]   NOTE: carried rows are keyed to the PREVIOUS frame "
                  "universe — re-run those stages if frames changed.")
        # mkstemp creates 0600 files; open the permissions up BEFORE the
        # swap so the live path never exists in a locked-down state.
        os.chmod(tmp, 0o644)
        os.replace(tmp, out_path)          # the atomic swap
        # The swap replaces the FILE, but SQLite's write-ahead log lives in
        # sidecar files named after it.  If any process had the OLD manifest
        # open in WAL mode, its -wal/-shm survive the rename and SQLite then
        # replays that log against the NEW database — which on 2026-08-19
        # produced "malformed database schema (reduced/2026-06-13/...)" and
        # made a freshly built, internally perfect manifest unreadable.
        # The log describes a database that no longer exists, so it must go
        # with it.  (Its contents are not lost work: carry_sibling_tables
        # read THROUGH the log when it copied the sibling tables.)
        for sidecar in (Path(str(out_path) + "-wal"), Path(str(out_path) + "-shm")):
            if sidecar.exists():
                print(f"[S0] removing stale {sidecar.name} left by a previous "
                      "reader/writer of the replaced manifest")
                sidecar.unlink()
    finally:
        if tmp.exists():                   # only on failure paths
            tmp.unlink()


# Tables S0 owns and rewrites from the catalog on every build.  Everything
# else in the manifest belongs to a downstream stage (S0b, S0c, S1, S2, S3,
# G, ...) and must survive an S0 rebuild — see carry_sibling_tables.
S0_OWNED_TABLES = ("frames", "aliases", "eras", "project_counts",
                   "build_meta")


def carry_sibling_tables(con: sqlite3.Connection, live_path: Path) -> list:
    """Copy tables S0 does not own from the live manifest into the new one.

    Used by the ``--fresh-file`` path only (the default in-place swap never
    moves a sibling table, so it has nothing to carry).  That path builds a
    fresh database in a temp file and atomically swaps it over the live
    manifest.  The swap replaces the whole FILE, so any table a later stage
    had added was silently destroyed: the 2026-08-18 ingest wiped
    ``s1_strata``, ``s1_solve_experiment``, ``s1_failure_autopsy`` and
    ``detector_params`` — the accepted evidence behind the astrometry
    go/no-go verdict and the detector memo — and the batch driver had to pin
    those numbers as code constants to survive.  Copying the sibling tables
    forward (exact original schema, then rows, then their indexes) makes a
    rebuild non-destructive.

    Returns the carried table names, which the caller records in
    ``build_meta`` and warns about: carried rows are keyed to the PREVIOUS
    frame universe, so a rebuild that changed ``frames`` leaves them STALE
    until their own stage re-runs.  Preserving stale evidence beats deleting
    it — the staleness is visible and fixable, deletion was neither.
    """
    if not Path(live_path).exists():
        return []                      # first build: nothing to carry
    # A concurrent stage may hold the write lock; wait rather than fail.
    con.execute("PRAGMA busy_timeout = 300000")
    con.execute("ATTACH DATABASE ? AS prev", (str(live_path),))
    try:
        tables = con.execute(
            "SELECT name, sql FROM prev.sqlite_master "
            "WHERE type = 'table' AND name NOT LIKE 'sqlite_%'").fetchall()
        carried = []
        for name, ddl in tables:
            if name in S0_OWNED_TABLES or not ddl:
                continue
            con.execute(ddl)           # replay the exact original schema
            con.execute(f'INSERT INTO main."{name}" '
                        f'SELECT * FROM prev."{name}"')
            carried.append(name)
        # Indexes belonging to carried tables (skip auto-indexes, which have
        # no DDL, and any index whose table we did not carry).
        for idx_sql, tbl in con.execute(
                "SELECT sql, tbl_name FROM prev.sqlite_master "
                "WHERE type = 'index' AND sql IS NOT NULL").fetchall():
            if tbl in carried:
                con.execute(idx_sql)
        # The DDL/INSERTs above opened an implicit transaction, and SQLite
        # refuses to DETACH while one is open ("database prev is locked").
        # Close it explicitly: commit the carry on success, roll it back on
        # failure — either way DETACH then runs on a quiet connection.
        con.commit()
    except Exception:
        con.rollback()
        raise
    finally:
        con.execute("DETACH DATABASE prev")
    return carried


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
    p.add_argument("--fresh-file", action="store_true",
                   help="rebuild as a NEW file and replace the live one "
                        "(sibling tables are carried across) instead of "
                        "swapping the S0 tables in place.  Only safe when "
                        "no other process has the manifest open")
    p.add_argument("--dry-run", action="store_true",
                   help="compute everything and print the dedup / era "
                        "summary, but write nothing")
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

    by_basis = frames["dup_basis"].value_counts()
    print("[S0]   duplicate evidence: " + "; ".join(
        f"{k}: {v:,}" for k, v in by_basis.items()))
    scanned = frames["hdr_scanned"].value_counts().to_dict()
    print(f"[S0]   header re-scrape: {scanned.get(HDR_SCANNED, 0):,} rows "
          f"scanned, {scanned.get(HDR_INHERITED, 0):,} inherited from an "
          f"exact copy, {scanned.get(HDR_NOT_SCANNED, 0):,} not scanned")
    for _, e in eras[eras["alias_of_era"].notna()].iterrows():
        print(f"[S0]   era {int(e['era_id'])} is a reduced ALIAS of era "
              f"{int(e['alias_of_era'])} ({e['alias_fraction']:.1%} of its "
              f"{int(e['n_rows']):,} rows; {int(e['n_canonical']):,} "
              "canonical left)")

    print("[S0] computing project reconciliation counts ...")
    counts = build_project_counts(frames, display_of_key)

    if args.dry_run:
        print("[S0] --dry-run: nothing written.")
        return 0

    print(f"[S0] writing manifest -> {args.out}")
    write_manifest(args.out, frames, aliases, eras, counts, args.catalog,
                   fresh_file=args.fresh_file)
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
