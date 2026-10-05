#!/usr/bin/env python
"""L0/L1/L2 — the legacy-archive census, its gates and its evidence page.

WHAT THIS SCRIPT DOES
---------------------
Reads the header scan (``scan``, ``disk_files``, ``keysets`` in
``products/legacy/legacy.sqlite``, written by ``build_legacy_scan.py``),
applies the pure rules of ``macro_legacy.census`` and writes the derived
tables into the same database:

L0 — what is on disk
* ``reconciliation``       one row per logical path (manifest ∪ disk) with
                           its status; ``reconciliation_summary`` the counts.
                           Acceptance: files-on-disk = rows + named exclusions.
* ``collision_audit``      what the first (BAD) crawl manifest would have
                           overwritten, and proof from the headers that the
                           files on disk are distinct exposures.
* ``frames``               one row per scanned file with everything derived:
                           camera, kind, filter band, target key, JD, night,
                           duplicate group, canonical flag, exclusion reason.
* ``cameras``              camera-identity census (the detector engineer's
                           list: INSTRUME, READOUTM, gain cards, set-points,
                           binning, geometry, software).
* ``camera_stints``        the camera timeline: first/last night per stint.
* ``mech_epochs``          rotation epochs (> 0.3° steps) inside each stint.

L1 — what it could support
* ``targets``, ``series``, ``runs``   nights × filters × longest same-filter
                           run, per target × filter × camera.
* ``overlap``              the RLMT-era project targets found in the legacy
                           archive, by name and by position.
* ``calib_census``, ``flat_pairs``, ``calib_by_camera``   calibration frames
                           per camera, and whether flat PAIRS exist (so gain
                           can be measured rather than read from a header).
* ``eclipsers``, ``eclipser_nights``  eclipsing/contact systems (VSX type and
                           ephemeris from ``build_legacy_external.py``) with
                           their minimum-bearing nights per season.
* ``time_audit``           the header-time convention audit per camera and
                           acquisition software.

L2 — the decision
* ``gates``                G0–G3 exactly as pre-registered, with the value
                           and threshold of each clause.
* ``census_meta``          provenance, including the sha256 of the
                           pre-registration note the gates implement.

Unless ``--skip-report`` is given it then renders
``docs/Legacy_Rigel/legacy_census.html`` from the database it just wrote.

ORDER OF OPERATIONS
-------------------
    build_legacy_scan.py        (slow; headers only)
    build_legacy_census.py      (fast; writes targets/series, gates G0/G2/G3)
    build_legacy_external.py    (VSX + SIMBAD for the targets found; cached)
    build_legacy_census.py      (again: now G1 and the positional overlap)

The census runs without the external tables — it then reports G1 and the
positional overlap as NOT EVALUATED rather than guessing — so the external
query can be driven by the census's own target list.

IDEMPOTENCE / SAFETY
--------------------
Every derived table is rebuilt from ``scan`` on each run and replaced inside
one transaction; ``scan``, ``disk_files``, ``keysets`` and the external
tables are never modified here.  The archive itself is not touched at all.
The RLMT manifest is opened read-only.

USAGE
-----
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_legacy_census.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sqlite3
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_legacy import LEGACY_CENSUS_VERSION               # noqa: E402
from macro_legacy import census as lc                        # noqa: E402
from macro_legacy import scan as lscan                       # noqa: E402

# ---------------------------------------------------------------------------
# Default locations
# ---------------------------------------------------------------------------
REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "legacy" / "legacy.sqlite"
ARCHIVE_PARENT = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO")
DEFAULT_MANIFEST = ARCHIVE_PARENT / "legacy_manifest.csv"
DEFAULT_BAD_MANIFEST = ARCHIVE_PARENT / "legacy_manifest_BAD_collisions.csv"
DEFAULT_RLMT = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
PREREG_NOTE = REPO_ROOT / "Legacy_Rigel" / "notes" / "GO_NOGO_PREREGISTERED.md"

#: Position match radius for the overlap query — the S0 cone (0.2°), so the
#: two archives agree on what "the same field" means.
OVERLAP_CONE_DEG = 0.2

#: RLMT projects whose targets are TIMING targets under gate G2 (the polars
#: and YZ Cnc): they additionally need one run ≥ 1 h.
TIMING_PROJECTS = frozenset({"CV_TimeSeries"})

#: Human-supplied (question, reader, venue) for a target that passes G3.
#: EMPTY BY DESIGN: the pre-registration lets G3 produce a go only if the
#: question and reader "can be named", and naming them is a decision for
#: the committee and James, recorded here by an explicit edit — never
#: inferred by the script.
G3_NAMED_QUESTIONS: dict[str, tuple[str, str, str]] = {}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""


def git_commit() -> str:
    """Short commit of the repo, '+dirty' if the tree has changes (read-only)."""
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO_ROOT, capture_output=True, text=True,
                             check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT,
                               capture_output=True, text=True).stdout.strip()
        return sha + ("+dirty" if dirty else "")
    except Exception:                      # noqa: BLE001 — provenance only
        return ""


def table_exists(con: sqlite3.Connection, name: str) -> bool:
    return bool(con.execute("SELECT count(*) FROM sqlite_master WHERE "
                            "type='table' AND name=?", (name,)).fetchone()[0])


def join_counts(values, limit: int = 8) -> str:
    """'value (n), value (n), …' — most common first, for census cells."""
    c = Counter(v for v in values if v is not None and not
                (isinstance(v, float) and np.isnan(v)))
    parts = [f"{_short(v)} ({n:,})" for v, n in c.most_common(limit)]
    if len(c) > limit:
        parts.append(f"… +{len(c) - limit} more")
    return ", ".join(parts)


def _short(v) -> str:
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def software_family(swcreate, telescop) -> str:
    """Acquisition-software label: 'MaxIm DL 6.08', or 'Talon (no SWCREATE)'.

    MaxIm writes ``SWCREATE = 'MaxIm DL Version 6.08 150819 2UT5A'``; the
    trailing build/serial tokens are dropped so one version is one family.
    The 2015 Rigel frames carry no SWCREATE at all; their card set (JD,
    RAWHENC, FOCUSPOS, HCOMSTAT …) is the Talon observatory software.
    """
    if isinstance(swcreate, str) and swcreate.strip():
        tokens = swcreate.split()
        if "Version" in tokens:
            i = tokens.index("Version")
            return " ".join(tokens[:i] + tokens[i + 1:i + 2])
        return " ".join(tokens[:3])
    if isinstance(telescop, str) and "rigel" in telescop.lower():
        return "Talon (no SWCREATE)"
    return "(no SWCREATE)"


# ---------------------------------------------------------------------------
# L0a — reconciliation against the transfer manifests
# ---------------------------------------------------------------------------
def load_manifest_paths(path: Path) -> list[tuple[str, str]]:
    """``[(file_id, rel_path), …]`` from a crawl manifest CSV."""
    with open(path, newline="") as fh:
        return [(r["file_id"], r["rel_path"]) for r in csv.DictReader(fh)]


def build_reconciliation(con, manifest: list[tuple[str, str]]) -> pd.DataFrame:
    """One row per logical path found in the manifest or on disk."""
    disk = pd.read_sql_query("SELECT path, kind FROM disk_files", con)
    errors = {r[0] for r in con.execute(
        "SELECT path FROM scan WHERE error IS NOT NULL")}
    in_manifest = {rel for _fid, rel in manifest}

    facts: dict[str, dict] = defaultdict(
        lambda: dict(has_fz=False, has_fits=False, fz_error=False,
                     is_fits_name=False))
    for p, kind in zip(disk["path"], disk["kind"]):
        f = facts[lc.logical_path(p)]
        if kind == "fz":
            f["has_fz"], f["is_fits_name"] = True, True
            f["fz_error"] = p in errors
        elif kind == "fits":
            f["has_fits"], f["is_fits_name"] = True, True
    for rel in in_manifest:
        facts[rel]["is_fits_name"] = facts[rel]["is_fits_name"] or \
            rel.lower().endswith((".fts", ".fit", ".fits"))

    rows = []
    for logical, f in facts.items():
        inm = logical in in_manifest
        rows.append(dict(
            logical_path=logical, in_manifest=int(inm),
            has_fz=int(f["has_fz"]), has_fits=int(f["has_fits"]),
            status=lc.reconcile_one(inm, f["has_fz"], f["has_fits"],
                                    f["fz_error"], f["is_fits_name"])))
    return pd.DataFrame(rows)


def reconciliation_summary(rec: pd.DataFrame, con) -> pd.DataFrame:
    """Counts that must close: files, rows, named exclusions.

    Two identities are checked and stored as rows of the summary:

    * FILE identity:  files on disk = scan rows + uncompressed twins of a
      scanned file + uncompressed-only files + non-FITS files.
    * MANIFEST identity:  manifest paths = scanned + unreadable +
      uncompressed-only + not-on-disk (each restricted to manifest paths).
    """
    n_disk = con.execute("SELECT count(*) FROM disk_files").fetchone()[0]
    n_fz = con.execute(
        "SELECT count(*) FROM disk_files WHERE kind='fz'").fetchone()[0]
    n_fits = con.execute(
        "SELECT count(*) FROM disk_files WHERE kind='fits'").fetchone()[0]
    n_other = con.execute(
        "SELECT count(*) FROM disk_files WHERE kind='other'").fetchone()[0]
    n_scan = con.execute("SELECT count(*) FROM scan").fetchone()[0]
    n_scan_err = con.execute(
        "SELECT count(*) FROM scan WHERE error IS NOT NULL").fetchone()[0]
    n_twin = int(((rec["has_fz"] == 1) & (rec["has_fits"] == 1)).sum())
    n_fits_only = int((rec["status"] == lc.REC_UNCOMPRESSED_ONLY).sum())
    m = rec[rec["in_manifest"] == 1]
    out = [
        ("files on disk (walk)", n_disk),
        ("  .fz files = scan rows", n_fz),
        ("    of which header unreadable (named exclusion)", n_scan_err),
        ("  uncompressed twin of a scanned .fz (named exclusion)", n_twin),
        ("  uncompressed with no .fz twin (named exclusion)", n_fits_only),
        ("  non-FITS files (named exclusion)", n_other),
        ("FILE IDENTITY residual (must be 0)",
         n_disk - (n_fz + n_twin + n_fits_only + n_other)),
        ("scan rows minus .fz files on disk (must be 0)", n_scan - n_fz),
        ("uncompressed files on disk (twins + twinless; check)",
         n_fits - (n_twin + n_fits_only)),
        ("manifest paths (crawl)", len(m)),
        ("  scanned", int((m["status"] == lc.REC_SCANNED).sum())),
        ("  .fz unreadable", int((m["status"] == lc.REC_UNREADABLE).sum())),
        ("  uncompressed only",
         int((m["status"] == lc.REC_UNCOMPRESSED_ONLY).sum())),
        ("  not on disk", int((m["status"] == lc.REC_NOT_ON_DISK).sum())),
        ("MANIFEST IDENTITY residual (must be 0)",
         len(m) - int(m["status"].isin(
             [lc.REC_SCANNED, lc.REC_UNREADABLE, lc.REC_UNCOMPRESSED_ONLY,
              lc.REC_NOT_ON_DISK]).sum())),
        ("on disk but not in the manifest (FITS)",
         int(((rec["in_manifest"] == 0)
              & (rec["status"] != lc.REC_NOT_FITS)).sum())),
    ]
    return pd.DataFrame(out, columns=["quantity", "n"]).assign(
        ord=range(len(out)))


def build_collision_audit(manifest, bad_manifest, frames: pd.DataFrame
                          ) -> tuple[pd.DataFrame, dict]:
    """What the BAD manifest would have overwritten — and whether it did.

    The first crawl stripped the year from every path, so e.g.
    ``2015/day004/frm00400.fts`` and ``2017/day004/frm00400.fts`` both
    became ``day004/frm00400.fts``.  Had the download used it, later files
    would have overwritten earlier ones.  Three checks:

    1. the BAD and good manifests list the same Drive file ids (the crawl
       itself lost nothing);
    2. re-stripping the good manifest reproduces the collision groups;
    3. for every collision group, the files now on disk carry DISTINCT
       header identities (DATE-OBS) — i.e. they are different exposures,
       each in its own year directory, not one file written N times.
    """
    good_ids = {fid for fid, _ in manifest}
    bad_ids = {fid for fid, _ in bad_manifest}
    bad_paths = Counter(rel for _, rel in bad_manifest)
    groups: dict[str, list[str]] = defaultdict(list)
    for _fid, rel in manifest:
        groups[lc.strip_first_component(rel)].append(rel)
    collide = {k: v for k, v in groups.items() if len(v) > 1}

    by_logical = frames.set_index("logical_path")[["date_obs", "error"]]
    rows = []
    for stripped, members in collide.items():
        present = [m for m in members if m in by_logical.index]
        stamps = [by_logical.at[m, "date_obs"] for m in present]
        stamps = [s for s in stamps if isinstance(s, str) and s]
        rows.append(dict(
            stripped_path=stripped, n_manifest=len(members),
            n_scanned=len(present), n_with_date=len(stamps),
            n_distinct_date_obs=len(set(stamps)),
            bad_manifest_rows=bad_paths.get(stripped, 0)))
    audit = pd.DataFrame(rows, columns=[
        "stripped_path", "n_manifest", "n_scanned", "n_with_date",
        "n_distinct_date_obs", "bad_manifest_rows"])
    meta = dict(
        collision_good_ids=len(good_ids), collision_bad_ids=len(bad_ids),
        collision_ids_only_in_good=len(good_ids - bad_ids),
        collision_ids_only_in_bad=len(bad_ids - good_ids),
        collision_bad_distinct_paths=len(bad_paths),
        collision_bad_rows=len(bad_manifest),
        collision_groups=len(collide),
        collision_paths_involved=sum(len(v) for v in collide.values()),
        collision_paths_at_risk=sum(len(v) - 1 for v in collide.values()),
        collision_groups_all_distinct=int(
            (audit["n_distinct_date_obs"] == audit["n_with_date"]).sum())
        if len(audit) else 0,
        collision_groups_identical_stamps=int(
            (audit["n_distinct_date_obs"] < audit["n_with_date"]).sum())
        if len(audit) else 0,
    )
    return audit, meta


# ---------------------------------------------------------------------------
# L0b — frames: everything derived from one header
# ---------------------------------------------------------------------------
def derive_frames(scan: pd.DataFrame) -> pd.DataFrame:
    """Add the derived columns to the scan rows (no row is dropped)."""
    f = scan.copy()
    # A column that is NULL on every row arrives from SQLite as dtype
    # object; force the numeric cards numeric so .round()/.median() work
    # whether or not any frame carries the card.
    for c in NUMERIC_COLUMNS:
        if c in f:
            f[c] = pd.to_numeric(f[c], errors="coerce")
    f["logical_path"] = [lc.logical_path(p) for p in f["path"]]
    f["basename"] = [p.rsplit("/", 1)[-1] for p in f["path"]]
    # An unreadable file has no INSTRUME; it must not be "identified".
    readable = f["error"].isna()
    f["tree"] = [p.split("/", 1)[0] if p.split("/", 1)[0] in
                 ("archivar", "old") else "dayNNN" for p in f["path"]]
    f["prefix"] = [b[:3].lower() for b in f["basename"]]

    # --- camera identity ---------------------------------------------------
    # Binning and native format first: the sensor signature (native pixel
    # pitch, native width) is the check on the INSTRUME label.
    f["binning"] = [lc.binning_of(a, b)
                    for a, b in zip(f["xbinning"], f["xfactor"])]
    f["native_x"] = f["naxis1"] * f["binning"]
    f["native_y"] = f["naxis2"] * f["binning"]
    f["pixel_um"] = f["xpixsz"] / f["binning"]
    f["camera_label"] = [lc.camera_name(v) for v in f["instrume"]]
    ident = [lc.camera_identity(i, p, x) for i, p, x in
             zip(f["instrume"], f["pixel_um"], f["native_x"])]
    f["camera"] = [c for c, _ in ident]
    f["camera_basis"] = [b for _, b in ident]
    f.loc[~readable, ["camera", "camera_basis", "camera_label"]] = None
    unknown = sorted(set(f.loc[f["camera"].isna() & f["error"].isna(),
                               "instrume"]))
    if unknown:
        raise SystemExit(
            "INSTRUME values no CAMERA_RULES entry recognises — add a rule, "
            f"do not guess: {unknown}")
    f["software"] = [software_family(s, t)
                     for s, t in zip(f["swcreate"], f["telescop"])]

    # --- time --------------------------------------------------------------
    f["jd_start"] = [lc.parse_date_obs(s) for s in f["date_obs"]]
    f["date_decimals"] = [lc.date_obs_decimals(s) for s in f["date_obs"]]
    f["night"] = [lc.night_label(j) for j in f["jd_start"]]
    f["jd_start"] = f["jd_start"].astype(float)       # None -> NaN, one dtype
    f["jd_card_diff_s"] = (f["jd"] - f["jd_start"]) * 86400.0
    f["lst_resid_s"] = [lc.lst_residual_seconds(js, lst)
                        for js, lst in zip(f["jd_start"], f["lst"])]
    f["path_day_offset"] = [lc.path_date_offset_days(p, js)
                            for p, js in zip(f["path"], f["jd_start"])]

    # --- what the frame is ---------------------------------------------------
    f["exptime_s"] = f["exptime"].where(f["exptime"].notna(), f["exposure"])
    f["kind"] = [lc.frame_kind(it, b, o, e) for it, b, o, e in
                 zip(f["imagetyp"], f["basename"], f["object"], f["exptime_s"])]
    # Only Talon's single letters name photometric bands (see filter_info).
    info = [lc.filter_info(v, talon=sw.startswith("Talon"))
            for v, sw in zip(f["filter"], f["software"])]
    f["band"] = [i.band for i in info]
    f["filter_system"] = [i.system for i in info]
    f["tieable"] = [int(i.tieable) for i in info]
    keys = [lc.target_key(o) for o in f["object"]]
    f["target_key"] = [k for k, _ in keys]
    f["target_name"] = [c for _, c in keys]

    # --- pointing (requested coordinates first: they name the TARGET) -------
    def ra_deg(row):
        for col in ("objra", "ra", "objctra"):
            v = lc.parse_sexagesimal(row[col])
            if v is not None:
                return v * 15.0
        return None

    def dec_deg(row):
        for col in ("objdec", "dec", "objctdec"):
            v = lc.parse_sexagesimal(row[col])
            if v is not None:
                return v
        return None
    coords = f[["objra", "ra", "objctra", "objdec", "dec", "objctdec"]]
    f["ra_deg"] = [ra_deg(r) for r in coords.to_dict("records")]
    f["dec_deg"] = [dec_deg(r) for r in coords.to_dict("records")]

    # --- the dedup rule ------------------------------------------------------
    key = [lc.dup_key(c, d, e, b, n1, n2, i) for c, d, e, b, n1, n2, i in
           zip(f["camera"], f["date_obs"], f["exptime_s"], f["band"],
               f["naxis1"], f["naxis2"], f.index)]
    f["dup_group"] = pd.factorize(pd.Series(key, dtype=object))[0]
    rank = [lc.canonical_rank(c, p) for c, p in zip(f["calstat"], f["path"])]
    f["_rank"] = rank
    best = f.sort_values("_rank").groupby("dup_group").head(1).index
    f["is_canonical"] = 0
    f.loc[best, "is_canonical"] = 1
    f = f.drop(columns="_rank")
    f["n_in_group"] = f.groupby("dup_group")["path"].transform("size")

    # --- why a frame is not in the science counts ---------------------------
    # One reason per frame, first match wins; NULL = a canonical science
    # frame.  These are the frame-level named exclusions.
    reason = np.full(len(f), None, dtype=object)
    reason[f["error"].notna().to_numpy()] = "header_unreadable"
    m = pd.isna(pd.Series(reason)) & f["jd_start"].isna().to_numpy()
    reason[m.to_numpy()] = "no_usable_date_obs"
    m = pd.isna(pd.Series(reason)) & (f["is_canonical"] == 0).to_numpy()
    reason[m.to_numpy()] = "duplicate_copy"
    m = pd.isna(pd.Series(reason)) & (f["kind"] != lc.KIND_LIGHT).to_numpy()
    reason[m.to_numpy()] = "not_science:" + f.loc[m.to_numpy(), "kind"]
    m = pd.isna(pd.Series(reason)) & f["target_key"].isna().to_numpy()
    reason[m.to_numpy()] = "no_target_name"
    f["exclusion"] = reason
    f["is_science"] = pd.isna(f["exclusion"]).astype(int)
    return f


# ---------------------------------------------------------------------------
# L0c — cameras, stints, mechanical epochs
# ---------------------------------------------------------------------------
def build_cameras(f: pd.DataFrame) -> pd.DataFrame:
    """Camera-identity census: one row per camera, every identity card."""
    rows = []
    ok = f[f["error"].isna()]
    for cam, g in ok.groupby("camera", dropna=False):
        dated = g[g["night"].notna()]
        canon = g[g["is_canonical"] == 1]
        rows.append(dict(
            camera=cam,
            n_files=len(g), n_canonical=len(canon),
            n_science=int(g["is_science"].sum()),
            n_no_date=int(g["jd_start"].isna().sum()),
            first_night=dated["night"].min() if len(dated) else None,
            last_night=dated["night"].max() if len(dated) else None,
            n_nights=dated["night"].nunique(),
            instrume=join_counts(g["instrume"].fillna("(none)")),
            identity_basis=join_counts(g["camera_basis"]),
            n_label_stale=int((g["camera_basis"]
                               == "sensor (label stale)").sum()),
            telescop=join_counts(g["telescop"]),
            focallen_mm=join_counts(g["focallen"]),
            aptdia_mm=join_counts(g["aptdia"]),
            native_pixel_um=join_counts(g["pixel_um"].round(2)),
            native_format=join_counts(
                [f"{int(a)}x{int(b)}" for a, b in
                 zip(g["native_x"], g["native_y"])
                 if pd.notna(a) and pd.notna(b)]),
            geometry=join_counts(
                [f"{int(a)}x{int(b)}" for a, b in zip(g["naxis1"], g["naxis2"])
                 if pd.notna(a) and pd.notna(b)]),
            binning=join_counts(g["binning"]),
            readoutm=join_counts(g["readoutm"]),
            egain=join_counts(g["egain"].round(3)),
            gain=join_counts(g["gain"]),
            offset=join_counts(g["offset"]),
            set_temp=join_counts(g["set_temp"].round(0)),
            ccd_temp_median=float(g["ccd_temp"].median())
            if g["ccd_temp"].notna().any() else None,
            camtemp_median=float(g["camtemp"].median())
            if g["camtemp"].notna().any() else None,
            software=join_counts(g["software"]),
            flipstat=join_counts(g["flipstat"]),
            calstat=join_counts(g["calstat"].fillna("(none)")),
            bitpix=join_counts(g["bitpix"]),
        ))
    return pd.DataFrame(rows).sort_values("first_night", na_position="last")


def build_stints(f: pd.DataFrame) -> pd.DataFrame:
    """The camera timeline: maximal blocks of consecutive nights per camera.

    Nights are ordered; a stint ends when the NEXT night with any frame is
    dominated by a different camera.  A night on which two cameras wrote
    frames is assigned to the camera with more frames and counted in
    ``n_mixed_nights`` — a swap during the night, or a mislabelled header.
    """
    d = f[f["night"].notna() & f["error"].isna() & (f["is_canonical"] == 1)]
    per = d.groupby(["night", "camera"]).size().reset_index(name="n")
    top = per.sort_values(["night", "n"], ascending=[True, False]) \
             .groupby("night").head(1).set_index("night")
    n_cams = per.groupby("night")["camera"].nunique()
    total = per.groupby("night")["n"].sum()
    stints: list[dict] = []
    for night in sorted(top.index):
        cam = top.at[night, "camera"]
        if stints and stints[-1]["camera"] == cam:
            s = stints[-1]
            s["last_night"] = night
            s["n_nights"] += 1
            s["n_frames"] += int(total[night])
            s["n_mixed_nights"] += int(n_cams[night] > 1)
        else:
            stints.append(dict(camera=cam, first_night=night, last_night=night,
                               n_nights=1, n_frames=int(total[night]),
                               n_mixed_nights=int(n_cams[night] > 1)))
    out = pd.DataFrame(stints)
    out.insert(0, "stint", range(1, len(out) + 1))
    return out


def build_mech_epochs(f: pd.DataFrame, stints: pd.DataFrame) -> pd.DataFrame:
    """Rotation epochs inside each camera stint (telescope engineer F1).

    Evidence: the plate-solution rotation the scheduler wrote into CROTA2,
    folded modulo 180° (a pier flip is not a re-mount), reduced to a nightly
    median over science frames, then segmented at steps > 0.3°.  Nights
    with fewer than ``ROTATION_MIN_FRAMES`` solved frames are not used.
    A stint with no usable rotation gets one epoch with NULL rotation: the
    camera identity still bounds it.
    """
    d = f[(f["is_science"] == 1) & f["crota2"].notna()].copy()
    d["rot"] = [lc.fold_rotation(v) for v in d["crota2"]]
    rows = []
    for s in stints.itertuples():
        g = d[(d["camera"] == s.camera) & (d["night"] >= s.first_night)
              & (d["night"] <= s.last_night)]
        nightly = g.groupby("night")["rot"].agg(["median", "size"])
        nightly = nightly[nightly["size"] >= lc.ROTATION_MIN_FRAMES]
        if nightly.empty:
            rows.append(dict(stint=s.stint, camera=s.camera,
                             first_night=s.first_night,
                             last_night=s.last_night, rotation_deg=None,
                             n_nights_measured=0, n_excursion_nights=0))
            continue
        for a, b, med, n, exc in lc.rotation_epochs(
                list(zip(nightly.index, nightly["median"]))):
            rows.append(dict(stint=s.stint, camera=s.camera, first_night=a,
                             last_night=b, rotation_deg=round(med, 3),
                             n_nights_measured=n, n_excursion_nights=exc))
    out = pd.DataFrame(rows)
    out.insert(0, "mech_epoch", range(1, len(out) + 1))
    return out


# ---------------------------------------------------------------------------
# L1a — calibration census
# ---------------------------------------------------------------------------
CALIB_KINDS = (lc.KIND_BIAS, lc.KIND_DARK, lc.KIND_FLAT)


def build_calib_census(f: pd.DataFrame) -> pd.DataFrame:
    """Calibration frames per camera × kind × binning × band."""
    c = f[(f["is_canonical"] == 1) & f["kind"].isin(CALIB_KINDS)
          & f["night"].notna()]
    rows = []
    for (cam, kind, binning, band), g in c.groupby(
            ["camera", "kind", "binning", "band"], dropna=False):
        rows.append(dict(
            camera=cam, kind=kind, binning=binning,
            band=band if kind == lc.KIND_FLAT else "—",
            n_frames=len(g), n_nights=g["night"].nunique(),
            first_night=g["night"].min(), last_night=g["night"].max(),
            exptimes=join_counts(g["exptime_s"].round(2), limit=6),
            set_temp=join_counts(g["set_temp"].round(0), limit=4)))
    out = pd.DataFrame(rows, columns=[
        "camera", "kind", "binning", "band", "n_frames", "n_nights",
        "first_night", "last_night", "exptimes", "set_temp"])
    if out.empty:
        return out
    # Bias and dark have no meaningful band: collapse their band rows.
    agg = out.groupby(["camera", "kind", "binning", "band"], dropna=False,
                      as_index=False).agg(
        n_frames=("n_frames", "sum"), n_nights=("n_nights", "max"),
        first_night=("first_night", "min"), last_night=("last_night", "max"),
        exptimes=("exptimes", "first"), set_temp=("set_temp", "first"))
    return agg


def build_flat_pairs(f: pd.DataFrame) -> pd.DataFrame:
    """Flat pairs per camera × binning × band × night (detector engineer).

    A photon-transfer gain needs PAIRS of flats at the same illumination;
    the header proxy is "same night, same filter, same exposure".  Only
    nights with at least one pair are stored.
    """
    c = f[(f["is_canonical"] == 1) & (f["kind"] == lc.KIND_FLAT)
          & f["night"].notna()]
    rows = []
    for (cam, binning, band, night), g in c.groupby(
            ["camera", "binning", "band", "night"], dropna=False):
        pairs, levels = lc.count_flat_pairs(list(g["exptime_s"]))
        if pairs:
            rows.append(dict(camera=cam, binning=binning, band=band,
                             night=night, n_flats=len(g), n_pairs=pairs,
                             n_levels=levels))
    return pd.DataFrame(rows, columns=["camera", "binning", "band", "night",
                                       "n_flats", "n_pairs", "n_levels"])


def calibration_index(f: pd.DataFrame) -> tuple[dict, dict]:
    """Sorted night ordinals of flats and of bias/dark, for the window test.

    Returns ``(flats, zeros)``: ``flats[(camera, binning, band)]`` and
    ``zeros[(camera, binning)]`` are sorted lists of night ordinals.
    """
    c = f[(f["is_canonical"] == 1) & f["night"].notna()]
    flats: dict = defaultdict(set)
    zeros: dict = defaultdict(set)
    for cam, binning, band, kind, night in zip(
            c["camera"], c["binning"], c["band"], c["kind"], c["night"]):
        if kind == lc.KIND_FLAT:
            flats[(cam, binning, band)].add(lc.night_ordinal(night))
        elif kind in (lc.KIND_BIAS, lc.KIND_DARK):
            zeros[(cam, binning)].add(lc.night_ordinal(night))
    return ({k: sorted(v) for k, v in flats.items()},
            {k: sorted(v) for k, v in zeros.items()})


# ---------------------------------------------------------------------------
# L1b — targets, series, runs
# ---------------------------------------------------------------------------
def build_runs(sci: pd.DataFrame) -> pd.DataFrame:
    """Every run (pre-registration §2) of every series on every night."""
    rows = []
    for (tkey, band, cam, night), g in sci.groupby(
            ["target_key", "band", "camera", "night"]):
        # One acquisition-software version per camera-night in practice;
        # the most common is recorded so a run can be tied to its time-audit
        # epoch (camera × software).
        sw = g["software"].mode().iloc[0]
        for a, b, n in lc.split_runs(list(g["jd_start"])):
            rows.append((tkey, band, cam, night, sw, a, b, n, (b - a) * 24.0))
    return pd.DataFrame(rows, columns=["target_key", "band", "camera", "night",
                                       "software", "jd_first", "jd_last",
                                       "n_frames", "run_hours"])


def build_series(sci: pd.DataFrame, runs: pd.DataFrame, flats: dict,
                 zeros: dict, seasons: dict) -> pd.DataFrame:
    """One row per series = target × band × camera (never pooled)."""
    longest = runs.groupby(["target_key", "band", "camera"])["run_hours"].max()
    rows = []
    for (tkey, band, cam), g in sci.groupby(["target_key", "band", "camera"]):
        nights = sorted(set(g["night"]))
        # A night is calibrated if ANY binning used on it has both a flat
        # (same camera, binning, band) and a bias/dark within the window.
        calibrated = set()
        for night, binning in set(zip(g["night"], g["binning"])):
            if lc.is_calibrated_night(
                    lc.night_ordinal(night),
                    flats.get((cam, binning, band), []),
                    zeros.get((cam, binning), [])):
                calibrated.add(night)
        per_night = runs[(runs["target_key"] == tkey) & (runs["band"] == band)
                         & (runs["camera"] == cam)].groupby("night")[
            "run_hours"].max()
        rows.append(dict(
            target_key=tkey, band=band, camera=cam,
            tieable=int(g["tieable"].iloc[0]),
            filter_system=g["filter_system"].iloc[0],
            n_frames=len(g), n_nights=len(nights),
            n_calibrated_nights=len(calibrated),
            n_calstat_nights=g.loc[g["calstat"].fillna("").str.strip() != "",
                                   "night"].nunique(),
            n_seasons=len({seasons[tkey][n] for n in nights}),
            first_night=nights[0], last_night=nights[-1],
            longest_run_h=float(longest[(tkey, band, cam)]),
            n_nights_run_1h=int((per_night >= 1.0).sum()),
            n_nights_run_3h=int((per_night >= 3.0).sum()),
            median_exptime_s=float(g["exptime_s"].median())))
    return pd.DataFrame(rows)


def build_targets(sci: pd.DataFrame, series: pd.DataFrame, seasons: dict
                  ) -> pd.DataFrame:
    """One row per alias-merged target."""
    rows = []
    best = series.sort_values("n_nights", ascending=False) \
                 .groupby("target_key").head(1).set_index("target_key")
    for tkey, g in sci.groupby("target_key"):
        nights = set(g["night"])
        rows.append(dict(
            target_key=tkey,
            target_name=g["target_name"].mode().iloc[0],
            n_raw_names=g["object"].nunique(),
            n_frames=len(g), n_nights=len(nights),
            n_seasons=len(set(seasons[tkey].values())),
            first_night=min(nights), last_night=max(nights),
            cameras=join_counts(g["camera"], limit=7),
            bands=join_counts(g["band"], limit=10),
            observers=g["observer"].nunique(),
            ra_deg=float(g["ra_deg"].median())
            if g["ra_deg"].notna().any() else None,
            dec_deg=float(g["dec_deg"].median())
            if g["dec_deg"].notna().any() else None,
            best_series_band=best.at[tkey, "band"],
            best_series_camera=best.at[tkey, "camera"],
            best_series_nights=int(best.at[tkey, "n_nights"]),
            longest_run_h=float(series.loc[series["target_key"] == tkey,
                                           "longest_run_h"].max())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# L1c — overlap with the RLMT-era targets
# ---------------------------------------------------------------------------
def angular_sep_deg(ra1, dec1, ra2, dec2):
    """Great-circle separation in degrees (vectorised haversine)."""
    r1, d1, r2, d2 = (np.radians(np.asarray(x, dtype=float))
                      for x in (ra1, dec1, ra2, dec2))
    a = (np.sin((d2 - d1) / 2) ** 2
         + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2)
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))


def build_overlap(con, sci: pd.DataFrame, series: pd.DataFrame,
                  seasons: dict
                  ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """RLMT-era project targets in the legacy archive, by name and position.

    ``rlmt_targets`` (written by ``build_legacy_external.py``) lists each
    project target with its key, project and catalogue position.  A legacy
    science frame belongs to the target if its alias key equals the
    target's key (``name``) or its requested coordinates lie within the S0
    cone of the catalogue position (``position``) — a frame can match on
    both.  Position matching is what finds a target observed under another
    name (``NGC 5457`` for M101) and what refuses a same-named field
    elsewhere.

    Returns ``(overlap, overlap_series, overlap_nights)``: one row per
    project target; one row per (project target, legacy series) with the G2
    clause values; and one row per (project target, night, filter, camera)
    — the matched frames on the calendar, which is what the page plots.
    """
    tg = pd.read_sql_query("SELECT * FROM rlmt_targets", con)
    o_rows, s_rows, n_rows = [], [], []
    for t in tg.itertuples():
        by_name = sci["target_key"] == t.target_key
        if pd.notna(t.ra_deg):
            sep = angular_sep_deg(sci["ra_deg"], sci["dec_deg"],
                                  t.ra_deg, t.dec_deg)
            by_pos = pd.Series(sep <= OVERLAP_CONE_DEG,
                               index=sci.index).fillna(False)
        else:
            by_pos = pd.Series(False, index=sci.index)
        g = sci[by_name | by_pos]
        o_rows.append(dict(
            project=t.project, target=t.target, target_key=t.target_key,
            n_frames=len(g), n_by_name=int(by_name.sum()),
            n_by_position=int(by_pos.sum()),
            n_name_not_position=int((by_name & ~by_pos).sum()),
            n_nights=g["night"].nunique(),
            first_night=g["night"].min() if len(g) else None,
            last_night=g["night"].max() if len(g) else None,
            legacy_names=join_counts(g["target_name"], limit=6),
            cameras=join_counts(g["camera"], limit=7),
            bands=join_counts(g["band"], limit=10)))
        if g.empty:
            continue
        for (night, band, cam), h in g.groupby(["night", "band", "camera"]):
            n_rows.append(dict(project=t.project, target=t.target,
                               target_key=t.target_key, night=night,
                               band=band, camera=cam, n_frames=len(h)))
        # Seasons are computed on the POOLED nights of the matched frames
        # (one target, possibly several legacy names).
        season_of = lc.assign_seasons(g["night"])
        timing = t.project in TIMING_PROJECTS
        for (band, cam), h in g.groupby(["band", "camera"]):
            nights = sorted(set(h["night"]))
            longest = 0.0
            for _night, hh in h.groupby("night"):
                for a, b, _n in lc.split_runs(list(hh["jd_start"])):
                    longest = max(longest, b - a)
            n_seasons = len({season_of[n] for n in nights})
            tieable = bool(h["tieable"].iloc[0])
            s_rows.append(dict(
                project=t.project, target=t.target, target_key=t.target_key,
                band=band, camera=cam, tieable=int(tieable),
                filter_system=h["filter_system"].iloc[0],
                n_frames=len(h), n_nights=len(nights), n_seasons=n_seasons,
                first_night=nights[0], last_night=nights[-1],
                longest_run_h=longest * 24.0, is_timing_target=int(timing),
                median_exptime_s=float(h["exptime_s"].median()),
                g2_pass=int(lc.g2_series_passes(
                    len(nights), n_seasons, tieable, timing, longest))))
    cols = ["project", "target", "target_key", "band", "camera", "tieable",
            "filter_system", "n_frames", "n_nights", "n_seasons",
            "first_night", "last_night", "longest_run_h", "is_timing_target",
            "median_exptime_s", "g2_pass"]
    night_cols = ["project", "target", "target_key", "night", "band",
                  "camera", "n_frames"]
    return (pd.DataFrame(o_rows), pd.DataFrame(s_rows, columns=cols),
            pd.DataFrame(n_rows, columns=night_cols))


def name_overlap_all(rlmt_path: Path, targets: pd.DataFrame) -> pd.DataFrame:
    """Every RLMT science target key that also names a legacy target.

    Wider than the project list: any of the ~2,000 RLMT-era target keys.
    Name-only (no positions), so it is a candidate list, not a claim.
    """
    con = sqlite3.connect(f"file:{rlmt_path}?mode=ro", uri=True)
    try:
        r = pd.read_sql_query("""
            SELECT target_key, min(canonical_target) AS rlmt_name,
                   count(*) AS rlmt_frames, count(DISTINCT night) AS rlmt_nights
            FROM frames WHERE is_canonical = 1 AND imagetyp LIKE 'Light%'
              AND target_key IS NOT NULL GROUP BY target_key""", con)
    finally:
        con.close()
    m = r.merge(targets[["target_key", "target_name", "n_frames", "n_nights",
                         "n_seasons", "first_night", "last_night",
                         "longest_run_h"]], on="target_key")
    return m.rename(columns={"n_frames": "legacy_frames",
                             "n_nights": "legacy_nights",
                             "n_seasons": "legacy_seasons"})


# ---------------------------------------------------------------------------
# L1d — eclipsing systems
# ---------------------------------------------------------------------------
def heliocentric_offset_days(jd: np.ndarray, ra_deg: float, dec_deg: float
                             ) -> np.ndarray:
    """HJD − JD(UTC) for a direction, in days (first-order, < 1 s error).

    Catalogue epochs of minimum are heliocentric; header times are
    geocentric UTC.  The Römer delay (≤ 8.3 min) matters against the 30-min
    margin of the predicted-minimum clause, so run limits are converted.
    The UTC→TT offset (~69 s) is deliberately NOT applied: VSX epochs are
    HJD on the UTC scale by convention.
    """
    from astropy.coordinates import get_sun
    from astropy.time import Time
    t = Time(np.asarray(jd, dtype=float), format="jd", scale="utc")
    sun = get_sun(t)
    # Earth→Sun unit vector scaled by distance, dotted with the target
    # direction: light reaches the Sun later than the Earth by r·n/c when
    # the target is on the far side.
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    n = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                  np.sin(dec)])
    s = sun.cartesian.xyz.to_value("au")           # geocentric Sun, AU
    earth_from_sun = -s
    return (earth_from_sun.T @ n) * 499.004784 / 86400.0


def build_eclipsers(con, targets: pd.DataFrame, runs: pd.DataFrame,
                    seasons: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Minimum-bearing nights per eclipsing system (pre-registration §2).

    Uses ``vsx_matches`` (type, period, epoch per legacy target).  For each
    system every run of every series is tested under clause (a) run ≥ P/2
    and clause (b) a predicted minimum inside the run with 30 min margins;
    a night is minimum-bearing if any of its runs satisfies either.
    """
    vsx = pd.read_sql_query("SELECT * FROM vsx_matches", con)
    vsx = vsx[vsx["vsx_name"].notna()]
    vsx["is_eclipsing"] = [lc.is_eclipsing_type(t) for t in vsx["vsx_type"]]
    vsx["is_transit"] = ["EP" in (t or "").upper().replace("/", "+").split("+")
                         for t in vsx["vsx_type"]]
    tmeta = targets.set_index("target_key")
    sys_rows, night_rows = [], []
    for v in vsx[vsx["is_eclipsing"] | vsx["is_transit"]].itertuples():
        if v.target_key not in tmeta.index:
            continue
        r = runs[runs["target_key"] == v.target_key]
        period = v.period if pd.notna(v.period) and v.period > 0 else None
        epoch = v.epoch if pd.notna(v.epoch) and v.epoch > 0 else None
        ra, dec = v.vsx_ra_deg, v.vsx_dec_deg
        if len(r) and epoch and period:
            off0 = heliocentric_offset_days(r["jd_first"].to_numpy(), ra, dec)
            off1 = heliocentric_offset_days(r["jd_last"].to_numpy(), ra, dec)
        else:
            off0 = off1 = np.zeros(len(r))
        per_night: dict[str, dict] = {}
        for run, o0, o1 in zip(r.itertuples(), off0, off1):
            days = run.run_hours / 24.0
            a = lc.run_guarantees_minimum(days, period)
            minima = lc.predicted_minima_in_run(
                run.jd_first + o0, run.jd_last + o1, epoch, period)
            d = per_night.setdefault(run.night, dict(
                guaranteed=False, predicted=False, n_predicted=0,
                longest_h=0.0))
            d["guaranteed"] |= a
            d["predicted"] |= bool(minima)
            d["n_predicted"] = max(d["n_predicted"], len(minima))
            d["longest_h"] = max(d["longest_h"], run.run_hours)
        season_of = seasons[v.target_key]
        for night, d in sorted(per_night.items()):
            night_rows.append(dict(
                target_key=v.target_key, night=night,
                season=season_of[night], longest_run_h=d["longest_h"],
                guaranteed=int(d["guaranteed"]), predicted=int(d["predicted"]),
                n_predicted_minima=d["n_predicted"],
                bearing=int(d["guaranteed"] or d["predicted"])))
        nights = pd.DataFrame([n for n in night_rows
                               if n["target_key"] == v.target_key])

        def n_seasons(col):
            return int(nights.loc[nights[col] == 1, "season"].nunique()) \
                if len(nights) else 0

        def n_nights(col):
            return int((nights[col] == 1).sum()) if len(nights) else 0
        bearing = nights[nights["bearing"] == 1] if len(nights) else nights
        sys_rows.append(dict(
            target_key=v.target_key,
            target_name=tmeta.at[v.target_key, "target_name"],
            vsx_name=v.vsx_name, vsx_type=v.vsx_type,
            is_eclipsing=int(v.is_eclipsing), is_transit=int(v.is_transit),
            period_d=period, epoch_hjd=epoch, match=v.match_method,
            sep_arcsec=v.sep_arcsec,
            n_nights=int(tmeta.at[v.target_key, "n_nights"]),
            n_seasons_observed=int(tmeta.at[v.target_key, "n_seasons"]),
            n_nights_guaranteed=n_nights("guaranteed"),
            n_nights_predicted=n_nights("predicted"),
            n_nights_bearing=n_nights("bearing"),
            n_seasons_guaranteed=n_seasons("guaranteed"),
            n_seasons_predicted=n_seasons("predicted"),
            n_seasons_bearing=n_seasons("bearing"),
            first_bearing_night=bearing["night"].min() if len(bearing) else None,
            last_bearing_night=bearing["night"].max() if len(bearing) else None,
            bands=tmeta.at[v.target_key, "bands"],
            cameras=tmeta.at[v.target_key, "cameras"]))
    sys_cols = ["target_key", "target_name", "vsx_name", "vsx_type",
                "is_eclipsing", "is_transit", "period_d", "epoch_hjd", "match",
                "sep_arcsec", "n_nights", "n_seasons_observed",
                "n_nights_guaranteed", "n_nights_predicted",
                "n_nights_bearing", "n_seasons_guaranteed",
                "n_seasons_predicted", "n_seasons_bearing",
                "first_bearing_night", "last_bearing_night", "bands",
                "cameras"]
    night_cols = ["target_key", "night", "season", "longest_run_h",
                  "guaranteed", "predicted", "n_predicted_minima", "bearing"]
    return (pd.DataFrame(sys_rows, columns=sys_cols),
            pd.DataFrame(night_rows, columns=night_cols))


# ---------------------------------------------------------------------------
# L1e — the header-time convention audit
# ---------------------------------------------------------------------------
def build_time_audit(f: pd.DataFrame) -> pd.DataFrame:
    """One row per camera × acquisition software.

    Columns answer, in order: does every frame have a usable DATE-OBS; how
    fine is the stamp; does the JD card agree with it; does the mount's LST
    agree with it (UTC, not local time; no gross offset); is the directory
    named for the UT date; and does the stamp mark the START of the
    exposure (the overlap test of ``stamp_hypothesis_violations``).
    """
    ok = f[f["error"].isna()]
    rows = []
    for (cam, sw), g in ok.groupby(["camera", "software"], dropna=False):
        dated = g[g["jd_start"].notna()]
        both = dated[dated["jd_card_diff_s"].notna()]
        n_both, n_agree, frac = lc.card_agreement(list(both["jd_card_diff_s"]))
        lst = dated["lst_resid_s"].dropna()
        # Overlap test on canonical frames, consecutive on the camera.
        canon = dated[(dated["is_canonical"] == 1)
                      & dated["exptime_s"].notna()].sort_values("jd_start")
        tot = Counter()
        for _night, h in canon.groupby("night"):
            tot.update(lc.stamp_hypothesis_violations(
                list(zip(h["jd_start"], h["exptime_s"]))))
        off = dated["path_day_offset"].dropna()
        dec = dated["date_decimals"].dropna()
        slope, icpt, n_fit = lst_exposure_slope(dated["exptime_s"],
                                                dated["lst_resid_s"])
        rows.append(dict(
            camera=cam, software=sw, n_files=len(g),
            n_usable_date=len(dated), n_no_usable_date=len(g) - len(dated),
            first_night=dated["night"].min() if len(dated) else None,
            last_night=dated["night"].max() if len(dated) else None,
            date_decimals=int(dec.mode().iloc[0]) if len(dec) else None,
            n_with_jd_card=n_both, n_jd_agree=n_agree,
            frac_jd_agree=frac if n_both else None,
            jd_diff_median_s=float(both["jd_card_diff_s"].median())
            if n_both else None,
            jd_diff_maxabs_s=float(both["jd_card_diff_s"].abs().max())
            if n_both else None,
            n_with_lst=len(lst),
            lst_resid_median_s=float(lst.median()) if len(lst) else None,
            lst_resid_p05_s=float(lst.quantile(0.05)) if len(lst) else None,
            lst_resid_p95_s=float(lst.quantile(0.95)) if len(lst) else None,
            frac_lst_within_60s=float((lst.abs() < 60).mean())
            if len(lst) else None,
            lst_slope_vs_exptime=None if np.isnan(slope) else slope,
            lst_intercept_s=None if np.isnan(icpt) else icpt,
            n_lst_fit=n_fit,
            n_path_dated=len(off),
            frac_path_ut_date=float((off == 0).mean()) if len(off) else None,
            n_path_offset_other=int((off.abs() > 1).sum()) if len(off) else 0,
            n_pairs=tot["pairs"], n_informative_pairs=tot["informative"],
            viol_start=tot["viol_start"], viol_mid=tot["viol_mid"],
            viol_end=tot["viol_end"], n_overlapping_pairs=tot["viol_start_all"],
        ))
    out = pd.DataFrame(rows).sort_values("first_night", na_position="last")
    out["convention"] = [convention_verdict(r) for r in out.itertuples()]
    out["card_pass"] = [int(bool(r.n_with_jd_card)
                            and r.frac_jd_agree >= lc.TIME_AGREE_FRACTION)
                        for r in out.itertuples()]
    return out


#: Verdict strings of the overlap test.
CONV_START = "start of exposure"
CONV_MID = "MIDDLE of exposure"
CONV_END = "END of exposure"


def convention_verdict(r) -> str:
    """What the overlap test says DATE-OBS marks, in words.

    Each hypothesis (start / middle / end) predicts a minimum spacing of
    consecutive stamps; the true convention is the one that is (almost)
    never violated while the others are violated often.  A stamp at
    mid-exposure violates BOTH the start and the end hypothesis and never
    the middle one — that signature is what identifies it.  When the
    exposures never differ enough, or overheads hide the difference, the
    verdict is 'not identifiable from headers' — never a default.
    """
    if not r.n_informative_pairs:
        return "not identifiable from headers (no unequal consecutive exposures)"
    tol = max(1, int(0.001 * r.n_informative_pairs))
    start_ok, mid_ok, end_ok = (r.viol_start <= tol, r.viol_mid <= tol,
                                r.viol_end <= tol)
    if start_ok and r.viol_end > 10 * max(r.viol_start, 1):
        return CONV_START
    if end_ok and r.viol_start > 10 * max(r.viol_end, 1):
        return CONV_END
    if mid_ok and not start_ok and not end_ok:
        return CONV_MID
    if start_ok and mid_ok and end_ok:
        return "not identifiable from headers (no hypothesis violated)"
    return "inconsistent (every hypothesis violated)"


def lst_exposure_slope(exptime, resid) -> tuple[float, float, int]:
    """Least-squares ``resid = intercept + slope · exptime`` (robustly clipped).

    The scheduler samples the mount's LST as the exposure begins.  If
    DATE-OBS is the start, the LST residual does not depend on exposure
    time (slope 0); if DATE-OBS is the middle, the residual is −exptime/2
    (slope −0.5); if the end, −exptime (slope −1).  This is independent of
    the overlap test: it uses a second card, not the spacing of frames.

    Frames with |residual| > 1000 s (stale LST cards) and zero exposures
    are dropped, then one 5-σ(MAD) clip is applied about a first fit.
    Returns ``(slope, intercept, n used)``; NaNs when fewer than 20 frames
    or a single exposure time.
    """
    x = np.asarray(exptime, dtype=float)
    y = np.asarray(resid, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (np.abs(y) < 1000)
    x, y = x[ok], y[ok]
    if len(x) < 20 or np.ptp(x) == 0:
        return float("nan"), float("nan"), int(len(x))
    slope, icpt = np.polyfit(x, y, 1)
    res = y - (icpt + slope * x)
    mad = np.median(np.abs(res - np.median(res))) * 1.4826
    keep = np.abs(res - np.median(res)) <= 5 * max(mad, 0.5)
    if keep.sum() >= 20 and np.ptp(x[keep]) > 0:
        slope, icpt = np.polyfit(x[keep], y[keep], 1)
    return float(slope), float(icpt), int(keep.sum())


def series_time_defects(sci: pd.DataFrame) -> pd.DataFrame:
    """Duplicate and overlapping stamps within the runs of each series.

    The pre-registered pass forbids "duplicate or non-monotonic timestamps
    within a run".  Per series: ``n_dup_stamps`` = canonical frames sharing
    a DATE-OBS with another frame of the same series (distinct exposures
    cannot share a start time); ``n_overlaps`` = consecutive frames of the
    series whose stamps are closer than the earlier exposure (they would
    have been exposing at once).
    """
    rows = []
    for (tkey, band, cam), g in sci.groupby(["target_key", "band", "camera"]):
        g = g.sort_values("jd_start")
        dup = int(g["date_obs"].duplicated(keep=False).sum())
        ov = 0
        for _night, h in g.groupby("night"):
            ov += lc.stamp_hypothesis_violations(
                list(zip(h["jd_start"], h["exptime_s"])))["viol_start_all"]
        rows.append((tkey, band, cam, dup, ov))
    return pd.DataFrame(rows, columns=["target_key", "band", "camera",
                                       "n_dup_stamps", "n_overlaps"])


# ---------------------------------------------------------------------------
# L2 — the gates
# ---------------------------------------------------------------------------
def evaluate_gates(summary: pd.DataFrame, stints: pd.DataFrame,
                   series: pd.DataFrame, overlap_series, eclipsers,
                   time_audit: pd.DataFrame, sci: pd.DataFrame,
                   recompute_g1=None) -> tuple[pd.DataFrame, dict]:
    """Evaluate G0–G3 as pre-registered.  Returns ``(gates table, meta)``.

    Each gate row carries the clause, its measured value, the threshold and
    pass (1/0) or NULL when it could not be evaluated (external catalogue
    not yet fetched) — an unevaluated gate is never treated as a pass.

    ``recompute_g1`` is a callable ``failing_epochs -> eclipsers table``
    used only for the amended reading (see the end of this function).
    """
    rows = []

    def add(gate, clause, value, threshold, passed, detail=""):
        rows.append(dict(gate=gate, clause=clause, value=str(value),
                         threshold=threshold,
                         passed=None if passed is None else int(passed),
                         detail=detail))

    # ---- G1: eclipsing systems with minimum-bearing nights in >= 3 seasons
    if eclipsers is None:
        g1 = None
        add("G1", "eclipsing systems with minimum-bearing nights in ≥ 3 "
            "seasons", "not evaluated", f"≥ {lc.G1_MIN_SYSTEMS}", None,
            "vsx_matches absent — run build_legacy_external.py")
        g1_series = series.iloc[0:0]
    else:
        ecl = eclipsers[eclipsers["is_eclipsing"] == 1]
        n_union = int((ecl["n_seasons_bearing"] >= lc.G1_MIN_SEASONS).sum())
        n_guar = int((ecl["n_seasons_guaranteed"] >= lc.G1_MIN_SEASONS).sum())
        n_pred = int((ecl["n_seasons_predicted"] >= lc.G1_MIN_SEASONS).sum())
        wide = eclipsers[(eclipsers["is_eclipsing"] == 1)
                         | (eclipsers["is_transit"] == 1)]
        n_wide = int((wide["n_seasons_bearing"] >= lc.G1_MIN_SEASONS).sum())
        g1 = lc.g1_passes(list(ecl["n_seasons_bearing"]))
        add("G1", "eclipsing systems with minimum-bearing nights in ≥ "
            f"{lc.G1_MIN_SEASONS} seasons (union of clauses a, b)", n_union,
            f"≥ {lc.G1_MIN_SYSTEMS}", g1,
            f"clause (a) run ≥ P/2 alone: {n_guar}; clause (b) predicted "
            f"minimum alone: {n_pred}; wider reading incl. planet transits: "
            f"{n_wide}")
        keys = set(ecl.loc[ecl["n_seasons_bearing"] >= lc.G1_MIN_SEASONS,
                           "target_key"]) if g1 else set()
        g1_series = series[series["target_key"].isin(keys)]

    # ---- G2: overlap
    if overlap_series is None:
        g2 = None
        add("G2", "RLMT-era target series passing", "not evaluated", "≥ 1",
            None, "rlmt_targets absent — run build_legacy_external.py")
        g2_pass = None
    else:
        g2_pass = overlap_series[overlap_series["g2_pass"] == 1]
        g2 = len(g2_pass) > 0
        names = sorted(set(g2_pass["target"]))
        add("G2", f"RLMT-era target series with ≥ {lc.G2_MIN_NIGHTS} nights "
            f"over ≥ {lc.G2_MIN_SEASONS} seasons in a tieable filter "
            "(timing targets: + a run ≥ 1 h)", len(g2_pass), "≥ 1", g2,
            "targets: " + (", ".join(names) if names else "none")
            + " — conditional on a pixel saturation test (standing rule 4)")

    # ---- G3: a long calibrated series
    g3_pass = series[series["n_calibrated_nights"]
                     >= lc.G3_MIN_CALIBRATED_NIGHTS]
    g3 = len(g3_pass) > 0
    best = series.sort_values("n_calibrated_nights", ascending=False).head(1)
    add("G3", f"series with ≥ {lc.G3_MIN_CALIBRATED_NIGHTS} calibrated nights",
        len(g3_pass), "≥ 1", g3,
        "largest: " + (f"{best['target_key'].iloc[0]} / "
                       f"{best['band'].iloc[0]} / {best['camera'].iloc[0]} = "
                       f"{int(best['n_calibrated_nights'].iloc[0])} calibrated "
                       f"of {int(best['n_nights'].iloc[0])} nights"
                       if len(best) else "no series"))
    named = [k for k in set(g3_pass["target_key"]) if k in G3_NAMED_QUESTIONS]
    add("G3", "passing targets with a named question and reader",
        len(named), "≥ 1 for a go under G3", bool(named) if g3 else None,
        "G3_NAMED_QUESTIONS is empty by design; naming is a human decision")

    # ---- G0: validity
    sm = dict(zip(summary["quantity"], summary["n"]))
    files_ok = (sm["FILE IDENTITY residual (must be 0)"] == 0
                and sm["scan rows minus .fz files on disk (must be 0)"] == 0
                and sm["MANIFEST IDENTITY residual (must be 0)"] == 0)
    add("G0", "files-on-disk = rows + named exclusions (both identities)",
        f"file residual {sm['FILE IDENTITY residual (must be 0)']}, "
        f"manifest residual {sm['MANIFEST IDENTITY residual (must be 0)']}",
        "0 and 0", files_ok)
    timeline_ok = len(stints) > 0 and stints["first_night"].notna().all()
    add("G0", "camera/mechanical timeline with first/last night per camera",
        f"{len(stints)} stints, {stints['camera'].nunique()} cameras",
        "exists", timeline_ok)

    # Time-convention pass for every camera/software epoch that contributes
    # frames to a PASSING gate.  An epoch is (camera, acquisition software);
    # it contributes if a science frame of a contributing series was written
    # in it.
    contrib = []
    if g1:
        contrib.append(g1_series[["target_key", "band", "camera"]])
    if g2:
        contrib.append(g2_pass[["target_key", "band", "camera"]])
    if g3:
        contrib.append(g3_pass[["target_key", "band", "camera"]])
    failing: set = set()
    if contrib:
        c = pd.concat(contrib).drop_duplicates()
        used_frames = sci.merge(c, on=["target_key", "band", "camera"])
        used = set(zip(used_frames["camera"], used_frames["software"]))
        epochs = time_audit[[(a, b) in used for a, b in
                             zip(time_audit["camera"], time_audit["software"])]]
        defects = series_time_defects(used_frames)
        n_dup = int(defects["n_dup_stamps"].sum())
        n_ov = int(defects["n_overlaps"].sum())
        is_start = epochs["convention"] == CONV_START
        is_utc = epochs["lst_resid_median_s"].abs() < 60
        is_card = epochs["card_pass"] == 1
        card_ok, start_ok, utc_ok = (bool(is_card.all()), bool(is_start.all()),
                                     bool(is_utc.all()))
        failing = set(zip(epochs.loc[~(is_start & is_utc & is_card), "camera"],
                          epochs.loc[~(is_start & is_utc & is_card),
                                     "software"]))
        time_ok = card_ok and start_ok and utc_ok and n_dup == 0 and n_ov == 0
        add("G0", "time convention: JD card within 2 s of DATE-OBS on ≥ 99% "
            "(contributing camera epochs)",
            f"{int(is_card.sum())} of {len(epochs)} epochs", "all", card_ok)
        add("G0", "time convention: DATE-OBS identified as exposure START by "
            "the overlap test (contributing camera epochs)",
            f"{int(is_start.sum())} of {len(epochs)} epochs", "all", start_ok,
            "; ".join(f"{a} / {b}: {v}" for a, b, v in
                      zip(epochs["camera"], epochs["software"],
                          epochs["convention"])))
        add("G0", "time convention: DATE-OBS is UTC (|median LST residual| "
            "< 60 s; contributing camera epochs)",
            f"{int(is_utc.sum())} of {len(epochs)} epochs", "all", utc_ok)
        add("G0", "no duplicate timestamps within contributing series",
            n_dup, "0", n_dup == 0)
        add("G0", "no overlapping (non-monotonic) timestamps within "
            "contributing series", n_ov, "0", n_ov == 0,
            "an overlap is the signature of a stamp that is not the start")
    else:
        time_ok = True
        add("G0", "time-convention pass for contributing camera epochs",
            "no gate passes — no contributing epoch", "all", True,
            "vacuous: the clause binds only epochs that feed a passing gate")
    g0 = files_ok and timeline_ok and time_ok

    undecided = g1 is None or g2 is None
    if undecided:
        outcomes = ["NOT YET EVALUATED — external catalogue tables absent"]
    else:
        outcomes = lc.decide(g0, bool(g1), bool(g2), bool(g3), bool(named))
    meta = dict(g0=int(g0), g1="" if g1 is None else int(g1),
                g2="" if g2 is None else int(g2), g3=int(g3),
                g3_question_named=int(bool(named)),
                outcome=" + ".join(outcomes),
                failing_time_epochs="; ".join(f"{a} / {b}"
                                              for a, b in sorted(failing)))

    # ---- The remedy the rule itself asks for ("fix the census"), reported
    # as an AMENDED reading beside the written one (pre-registration §5):
    # drop every run taken in a camera epoch that fails the time-convention
    # clause and re-evaluate G1 on what is left.  G0's time clause is then
    # satisfied by construction for the surviving epochs.
    meta["outcome_amended"] = ""
    if g1 and failing and files_ok and timeline_ok and recompute_g1 is not None:
        ecl2 = recompute_g1(failing)
        ecl2 = ecl2[ecl2["is_eclipsing"] == 1]
        n2 = int((ecl2["n_seasons_bearing"] >= lc.G1_MIN_SEASONS).sum())
        g1b = lc.g1_passes(list(ecl2["n_seasons_bearing"]))
        names = ", ".join(sorted(ecl2.loc[
            ecl2["n_seasons_bearing"] >= lc.G1_MIN_SEASONS, "target_name"]))
        add("G1*", "AMENDED: G1 re-evaluated without runs from camera epochs "
            "that fail the time-convention clause", n2,
            f"≥ {lc.G1_MIN_SYSTEMS}", g1b,
            f"excluded epochs: {meta['failing_time_epochs']}; systems: "
            f"{names or 'none'}")
        meta["g1_amended"] = int(g1b)
        meta["g1_amended_systems"] = names
        meta["outcome_amended"] = " + ".join(
            lc.decide(True, g1b, bool(g2), bool(g3), bool(named)))
    return pd.DataFrame(rows), meta


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def write_tables(con: sqlite3.Connection, tables: dict[str, pd.DataFrame],
                 meta: dict) -> None:
    """Replace every derived table and the meta in ONE transaction."""
    con.execute("BEGIN")
    try:
        for name, df in tables.items():
            con.execute(f'DROP TABLE IF EXISTS "{name}"')
            df.to_sql(name, con, index=False)
        con.execute("DROP TABLE IF EXISTS census_meta")
        con.execute("CREATE TABLE census_meta (key TEXT PRIMARY KEY, value TEXT)")
        con.executemany("INSERT INTO census_meta VALUES (?,?)",
                        [(k, str(v)) for k, v in meta.items()])
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    for name, cols in INDEXES.items():
        if name in tables:
            con.execute(f'CREATE INDEX IF NOT EXISTS "ix_{name}" '
                        f'ON "{name}" ({cols})')
    con.commit()


#: Indexes that make the report's queries (and a reader's) fast.
INDEXES = {"frames": "target_key, band, camera", "runs": "target_key",
           "series": "target_key", "reconciliation": "status"}

#: Scan columns carried into ``frames`` (the rest stay in ``scan``).
FRAME_SCAN_COLUMNS = [
    "path", "size", "error", "telescop", "instrume", "readoutm", "swcreate",
    "imagetyp", "filter", "object", "observer", "calstat", "calstart",
    "date_obs", "lst", "flipstat", "objra", "ra", "objctra", "objdec", "dec",
    "objctdec", "xbinning", "xfactor", "xpixsz", "gain", "egain", "offset",
    "set_temp", "ccd_temp", "camtemp", "focallen", "aptdia", "focuspos",
    "exptime", "exposure", "jd", "jd_helio", "airmass", "crota2", "fwhmh",
    "fwhmv", "naxis1", "naxis2", "bitpix", "keyset_sig",
]
#: Scan columns that hold numbers (see ``derive_frames``).
NUMERIC_COLUMNS = ([lscan.col(k) for k in lscan.NUM_KEYS]
                   + ["size", "naxis1", "naxis2", "bitpix"])
#: Helper columns of ``frames`` not worth storing (recoverable from ``scan``).
FRAME_DROP = ["objra", "ra", "objctra", "objdec", "dec", "objctdec",
              "xbinning", "xfactor", "exposure"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--bad-manifest", type=Path, default=DEFAULT_BAD_MANIFEST)
    ap.add_argument("--rlmt-manifest", type=Path, default=DEFAULT_RLMT)
    ap.add_argument("--skip-report", action="store_true")
    args = ap.parse_args(argv)

    con = sqlite3.connect(args.db, timeout=120)
    # Write-ahead logging: appends sequentially instead of journalling every
    # replaced page — on the spinning disk this database lives on, replacing
    # the derived tables is otherwise tens of minutes.
    con.execute("PRAGMA journal_mode=WAL")
    try:
        if not table_exists(con, "scan") or not table_exists(con, "disk_files"):
            print("scan/disk_files absent — run build_legacy_scan.py first",
                  file=sys.stderr)
            return 1
        cols = ", ".join(f'"{c}"' for c in FRAME_SCAN_COLUMNS)
        scan = pd.read_sql_query(f"SELECT {cols} FROM scan", con)
        print(f"scan rows: {len(scan):,}", flush=True)

        # ---- L0 -----------------------------------------------------------
        frames = derive_frames(scan)
        manifest = load_manifest_paths(args.manifest)
        bad = load_manifest_paths(args.bad_manifest)
        rec = build_reconciliation(con, manifest)
        summary = reconciliation_summary(rec, con)
        collisions, coll_meta = build_collision_audit(manifest, bad, frames)
        cameras = build_cameras(frames)
        stints = build_stints(frames)
        mech = build_mech_epochs(frames, stints)
        print(f"canonical frames: {int(frames['is_canonical'].sum()):,}; "
              f"science: {int(frames['is_science'].sum()):,}; "
              f"cameras: {len(cameras)}; stints: {len(stints)}", flush=True)

        # ---- L1 -----------------------------------------------------------
        sci = frames[frames["is_science"] == 1]
        seasons = {k: lc.assign_seasons(g["night"])
                   for k, g in sci.groupby("target_key")}
        flats, zeros = calibration_index(frames)
        runs = build_runs(sci)
        series = build_series(sci, runs, flats, zeros, seasons)
        targets = build_targets(sci, series, seasons)
        calib = build_calib_census(frames)
        pairs = build_flat_pairs(frames)
        audit = build_time_audit(frames)
        defects = series_time_defects(sci)
        series = series.merge(defects, on=["target_key", "band", "camera"])
        overlap_all = name_overlap_all(args.rlmt_manifest, targets)
        print(f"targets: {len(targets):,}; series: {len(series):,}; "
              f"runs: {len(runs):,}", flush=True)

        tables = dict(
            reconciliation=rec, reconciliation_summary=summary,
            collision_audit=collisions, cameras=cameras,
            camera_stints=stints, mech_epochs=mech, calib_census=calib,
            flat_pairs=pairs, runs=runs, series=series, targets=targets,
            time_audit=audit, overlap_all=overlap_all)

        overlap_series = eclipsers = None
        if table_exists(con, "rlmt_targets"):
            overlap, overlap_series, overlap_nights = build_overlap(
                con, sci, series, seasons)
            tables.update(overlap=overlap, overlap_series=overlap_series,
                          overlap_nights=overlap_nights)
        if table_exists(con, "vsx_matches"):
            eclipsers, ecl_nights = build_eclipsers(con, targets, runs, seasons)
            tables.update(eclipsers=eclipsers, eclipser_nights=ecl_nights)

        # ---- L2 -----------------------------------------------------------
        def recompute_g1(failing_epochs):
            """Eclipser table from runs outside the failing camera epochs."""
            keep = [(c, w) not in failing_epochs
                    for c, w in zip(runs["camera"], runs["software"])]
            return build_eclipsers(con, targets, runs[keep], seasons)[0]

        gates, gate_meta = evaluate_gates(
            summary, stints, series, overlap_series, eclipsers, audit, sci,
            recompute_g1 if eclipsers is not None else None)
        tables["gates"] = gates
        tables["frames"] = frames.drop(columns=FRAME_DROP)

        meta = dict(
            built_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            code_version=LEGACY_CENSUS_VERSION, git_commit=git_commit(),
            prereg_path=str(PREREG_NOTE.relative_to(REPO_ROOT)),
            prereg_sha256=sha256_of(PREREG_NOTE),
            manifest_sha256=sha256_of(args.manifest),
            bad_manifest_sha256=sha256_of(args.bad_manifest),
            run_gap_min=lc.RUN_GAP_DAYS * 1440,
            season_gap_days=lc.SEASON_GAP_DAYS,
            calib_window_nights=lc.CALIB_WINDOW_NIGHTS,
            rotation_step_deg=lc.ROTATION_STEP_DEG,
            overlap_cone_deg=OVERLAP_CONE_DEG,
            **coll_meta, **gate_meta)
        write_tables(con, tables, meta)
        print("gates:")
        for g in gates.itertuples():
            print(f"  {g.gate}: {g.clause} = {g.value} "
                  f"[{ {1: 'PASS', 0: 'FAIL'}.get(g.passed, 'n/a') }]")
        print(f"OUTCOME (as written): {gate_meta['outcome']}")
        if gate_meta.get("outcome_amended"):
            print(f"OUTCOME (amended, failing time epochs excluded): "
                  f"{gate_meta['outcome_amended']}")
    finally:
        con.close()

    if not args.skip_report:
        from macro_legacy import report
        out = report.render_report(args.db)
        print(f"report: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
