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
* ``mech_epoch``         — the HARDWARE history beneath the eras: one row per
                           mechanical state of the telescope (camera, flip
                           state, wheel map, measured sky rotation), with
                           its first and last night and why it began.
* ``night_mech_epoch``   — which epoch each (camera, night) belongs to, and
                           whether that placement is certain.
* ``frame_mech_epoch``   — the same, joined down to every frame (the table
                           a stage actually joins on ``obs_rowid``).
* ``mech_rotation_null`` — every same-target rotation difference measured
                           INSIDE an epoch: the false-alarm test of the
                           0.3-degree step rule.
* ``mech_injection``     — injected-step recovery per epoch and step size:
                           what the segmentation would have found, and how
                           late, had the camera been turned.
* ``calib_coverage_mech`` / ``calib_gaps_mech`` — the coverage matrix and
                           shopping list re-counted per (era, mechanical
                           epoch) under the rule that no calibration
                           crosses a boundary that invalidates it, keyed
                           additionally on gain, offset and set-point.
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

import numpy as np
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

#: Link method recorded for a reduced row S0 ALREADY grouped with its raw
#: parent, by the evidence S0 used (``frames.dup_basis``).  The vocabulary
#: is the one S0b has always published — S2's reconstruction experiment
#: selects links by these names — so a ``_calibrated`` twin is ``stem_jd``
#: whether S0's dedup found it (v1.1) or this stage's ladder did (v1.0).
#: Only the counter form is new, and it is named apart because it is the
#: weakest evidence in the table.
LINK_METHOD_OF_DUP_BASIS = {
    "same_basename_jd": "same_basename_jd",
    "packaging": "stem_jd",
    "processing_suffix": "stem_jd",
    "suffix_counter": "stem_jd_counter",
}

#: The hardware-state columns S0 v1.1 adds to ``frames`` (F-2).  Loaded
#: when present; a manifest built before the header re-scrape simply has
#: none of them, every frame's camera is then unknown, and the mechanical
#: tables come out empty rather than wrong.
HARDWARE_COLUMNS = [
    "camera", "flipstat", "fwallnam", "swcreate", "telpier",
    "hdr_gain", "hdr_offset", "set_temp",
]

#: Injected rotation steps for the recovery test, degrees.  Brackets the
#: 0.3-degree rule from below (0.2: must NOT be found) to far above it.
INJECTION_STEPS_DEG = (0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0)
INJECTION_TRIALS = 200
INJECTION_SEED = 20261003
#: An epoch needs at least this many nights to host an injection (the step
#: is placed strictly inside it, with nights on both sides).
INJECTION_MIN_NIGHTS = 10


# ---------------------------------------------------------------------------
# Step 1 — load the manifest frames (plus eras and project_counts)
# ---------------------------------------------------------------------------
def load_frames(con: sqlite3.Connection) -> pd.DataFrame:
    """Read the S0 frames table (selected columns) into a DataFrame."""
    have = {r[1] for r in con.execute("PRAGMA table_info(frames)")}
    cols = ", ".join(FRAME_COLUMNS
                     + [c for c in HARDWARE_COLUMNS if c in have])
    df = pd.read_sql_query(f"SELECT {cols} FROM frames", con)
    for c in HARDWARE_COLUMNS:
        if c not in df:
            df[c] = None
    # frames.dup_basis exists from S0 v1.1; an older manifest only ever
    # grouped exact copies, which is what the default says.
    if "dup_basis" in have:
        df["dup_basis"] = pd.read_sql_query(
            "SELECT dup_basis FROM frames", con)["dup_basis"]
    else:
        df["dup_basis"] = "same_basename_jd"
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
    # frames.dup_basis exists from S0 v1.1; without it every grouped row
    # is an exact copy, the only kind S0 v1.0 could group.
    basis = (red["dup_basis"] if "dup_basis" in red
             else ["same_basename_jd"] * len(red))
    for (rrid, rpath, rbn, rjd, rnight, rtkey, rgroup, rbasis) in zip(
            red["obs_rowid"], red["path"], red["basename"], red["jd"],
            red["night"], red["target_key"], red["dup_group"], basis):
        rjd = None if pd.isna(rjd) else float(rjd)
        head = head_of_group.get(rgroup)
        if head is not None:
            # S0's own dedup already grouped this row with a raw frame; the
            # method is the evidence S0 recorded for it (an exact copy, or
            # since S0 v1.1 a processed/re-packaged twin).
            matches = [(int(head), LINK_METHOD_OF_DUP_BASIS.get(
                rbasis, "same_basename_jd"), 0.0)]
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
    # Appended AFTER the original columns (F-3 / DE.F5), so every existing
    # reader of this table keeps its column positions: where the frame sits
    # in the hardware history, and the camera settings it was taken at.
    for col in ("camera", "mech_epoch", "detector_epoch", "epoch_certain",
                "hdr_gain", "hdr_offset", "set_temp"):
        out[col] = cal[col] if col in cal else None
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Step 3b — mechanical epochs (F-3)
# ---------------------------------------------------------------------------
def load_rotations(con: sqlite3.Connection) -> pd.DataFrame:
    """Position angles S1b measured on plate-solved frames.

    One row per solved frame: its catalog id, the alias key it was solved
    under, the rotation, and where it actually pointed.  Returns an empty
    frame when the S1b table does not exist yet (a manifest that has never
    been plate-solved has no rotation evidence, and the segmentation then
    rests on the headers alone — which it says, per epoch).
    """
    cols = ["obs_rowid", "target_key", "rotation_deg", "solved_ra",
            "solved_dec"]
    has = con.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                      "AND name='s1_batch'").fetchone()
    if not has:
        return pd.DataFrame(columns=cols)
    return pd.read_sql_query(
        f"SELECT {', '.join(cols)} FROM s1_batch WHERE status = 'solved' "
        "AND rotation_deg IS NOT NULL AND solved_ra IS NOT NULL "
        "AND solved_dec IS NOT NULL", con)


def _modal(values: pd.Series):
    """Most frequent non-null value (ties: smallest), or None."""
    vals = values.dropna()
    if not len(vals):
        return None
    counts = vals.value_counts()
    top = counts[counts == counts.iloc[0]].index
    return sorted(top)[0]


def build_night_states(df: pd.DataFrame, rot: pd.DataFrame
                       ) -> list[inv.NightState]:
    """One :class:`inventory.NightState` per (camera, night) with frames.

    Built from the CANONICAL, error-free frames — science and calibration
    alike, so that a day of bench darks is a node on the timeline too and
    every calibration frame can be placed.  Header values are the night's
    MODE over the frames that carry the card; rotation samples are joined
    from the S1b solves by catalog id.
    """
    use = df[(df["is_canonical"] == 1) & df["error"].isna()
             & df["camera"].notna() & df["night"].notna()]
    rot_of: dict[tuple, list] = {}
    if len(rot):
        key = use.set_index("obs_rowid")[["camera", "night"]]
        joined = rot.join(key, on="obs_rowid", how="inner")
        for cam, night, tk, r, ra, dec in zip(
                joined["camera"], joined["night"], joined["target_key"],
                joined["rotation_deg"], joined["solved_ra"],
                joined["solved_dec"]):
            rot_of.setdefault((cam, night), []).append(
                inv.RotSample(tk if isinstance(tk, str) else None,
                              float(r), float(ra), float(dec)))
    states = []
    for (cam, night), grp in use.groupby(["camera", "night"], sort=True):
        # One modal wheel map per SLOT COUNT: maps of different lengths are
        # different descriptions of the wheels (see NightState) and are
        # tracked side by side, never compared with each other.
        maps = grp["fwallnam"].dropna()
        by_len: dict[int, pd.Series] = {}
        for n_slots, sub in maps.groupby(maps.str.count(r"\|") + 1):
            by_len[int(n_slots)] = sub
        states.append(inv.NightState(
            camera=cam, night=night, n_frames=int(len(grp)),
            flipstat=_modal(grp["flipstat"]),
            wheel_map=_modal(grp["fwallnam"]),
            swcreate=_modal(grp["swcreate"]),
            telpier=_modal(grp["telpier"]),
            rot=tuple(rot_of.get((cam, night), ())),
            wheel_maps=tuple(_modal(by_len[n]) for n in sorted(by_len))))
    return states


def injection_test(states: list[inv.NightState], epochs: pd.DataFrame
                   ) -> pd.DataFrame:
    """Would the segmentation have found a step, had there been one?

    The null distribution says how often the rule fires when nothing
    changed.  This is the other half: how often it FAILS to fire when
    something did — the dangerous direction, because an undetected
    re-seat lets a flat be applied across it.

    For each real epoch of at least :data:`INJECTION_MIN_NIGHTS` nights, and
    each step size in :data:`INJECTION_STEPS_DEG`, a step is added to every
    rotation sample from a randomly chosen interior night onward, and the
    epoch's own nights are re-segmented.  Recorded per (epoch, step):

    * ``frac_detected`` — a rotation boundary appeared at or after the
      injected night;
    * ``frac_exact`` — it appeared ON the injected night;
    * ``median_delay_nights`` / ``p90_delay_nights`` — how many of the
      camera's nights passed, mis-filed under the old epoch, before the
      step was seen (over detected trials).

    The rotation samples, target lists and night spacing are the REAL ones,
    so the answer reflects each epoch's actual coverage — an epoch observed
    as a different target every night has no same-target pairs and can
    only ever see the coarse test.
    """
    rng = np.random.default_rng(INJECTION_SEED)
    by_epoch: dict[str, list[inv.NightState]] = {}
    epoch_of = {}
    for _, e in epochs.iterrows():
        epoch_of[(e["camera"], e["first_night"], e["last_night"])] = \
            e["mech_epoch"]
    for st in states:
        for (cam, lo, hi), eid in epoch_of.items():
            if st.camera == cam and lo <= st.night <= hi:
                by_epoch.setdefault(eid, []).append(st)
                break
    rows = []
    for eid, sts in by_epoch.items():
        sts = sorted(sts, key=lambda s: s.night)
        if len(sts) < INJECTION_MIN_NIGHTS:
            continue
        for step in INJECTION_STEPS_DEG:
            delays = []
            for _ in range(INJECTION_TRIALS):
                k = int(rng.integers(2, len(sts) - 1))
                shifted = sts[:k] + [
                    inv.NightState(
                        camera=s.camera, night=s.night, n_frames=s.n_frames,
                        flipstat=s.flipstat, wheel_map=s.wheel_map,
                        rot=tuple(inv.RotSample(r.target,
                                                (r.rot_deg + step) % 360.0,
                                                r.ra_deg, r.dec_deg)
                                  for r in s.rot))
                    for s in sts[k:]]
                found, _, _ = inv.segment_mech_epochs(shifted)
                # Boundaries the injection created: rotation-caused epochs
                # beginning at or after the injected night.
                hits = [e["first_night"] for e in found
                        if "rotation" in e["boundary_cause"]
                        and e["first_night"] >= sts[k].night]
                if hits:
                    first = min(hits)
                    delays.append(sum(1 for s in sts[k:]
                                      if s.night < first))
            n_det = len(delays)
            rows.append({
                "mech_epoch": eid, "n_nights": len(sts),
                "step_deg": step, "n_trials": INJECTION_TRIALS,
                "n_detected": n_det,
                "frac_detected": n_det / INJECTION_TRIALS,
                "frac_exact": sum(1 for d in delays if d == 0)
                / INJECTION_TRIALS,
                "median_delay_nights": float(np.median(delays))
                if delays else None,
                "p90_delay_nights": float(np.percentile(delays, 90))
                if delays else None,
            })
    return pd.DataFrame(rows, columns=[
        "mech_epoch", "n_nights", "step_deg", "n_trials", "n_detected",
        "frac_detected", "frac_exact", "median_delay_nights",
        "p90_delay_nights"])


def build_mech_tables(df: pd.DataFrame, rot: pd.DataFrame) -> dict:
    """Segment the hardware history and join it down to every frame.

    Returns the five mechanical tables as DataFrames keyed by table name.
    ``frame_mech_epoch`` covers EVERY manifest row whose camera and night
    are known — duplicates included, since a copy was taken on the same
    night by the same camera as its original.
    """
    states = build_night_states(df, rot)
    epochs, nights, null = inv.segment_mech_epochs(states)
    epoch_cols = ["seq", "mech_epoch", "camera", "first_night", "last_night",
                  "n_nights", "n_frames", "boundary_cause", "detector_epoch",
                  "rotation_deg", "rotation_mad_deg", "n_rot_nights",
                  "n_rot_frames", "step_deg", "step_err_deg",
                  "step_n_targets", "step_basis", "step_marginal",
                  "flipstat", "wheel_map", "swcreate", "telpier",
                  "n_transitions", "n_tested_fine", "n_tested_coarse",
                  "n_gap_nights"]
    for e in epochs:
        # A rotation step is MARGINAL when it does not clear the threshold
        # by two standard errors: such a boundary is kept (an extra boundary
        # only narrows which flats are shared) but flagged for a human.
        fine = inv.CAUSE_ROT_FINE in e["boundary_cause"]
        e["step_marginal"] = int(
            fine and e["step_deg"] is not None
            and abs(e["step_deg"]) - 2.0 * (e["step_err_deg"] or 0.0)
            <= inv.MECH_ROT_STEP_DEG)
    ep = pd.DataFrame(epochs, columns=epoch_cols)
    ni = pd.DataFrame(nights, columns=[
        "camera", "night", "mech_epoch", "detector_epoch", "certain",
        "basis", "n_frames", "rot_median_deg", "rot_mad_deg", "n_rot",
        "flipstat", "wheel_map"])
    nu = pd.DataFrame(null, columns=[
        "camera", "night", "mech_epoch", "target", "delta_deg",
        "delta_se_deg", "dec_deg"])
    # Join down to frames on (camera, night).
    key = ni[["camera", "night", "mech_epoch", "detector_epoch", "certain"]]
    fm = df[["obs_rowid", "camera", "night"]].merge(
        key, on=["camera", "night"], how="inner")
    fm = fm.rename(columns={"certain": "epoch_certain"})
    inj = injection_test(states, ep)
    return {"mech_epoch": ep, "night_mech_epoch": ni,
            "frame_mech_epoch": fm, "mech_rotation_null": nu,
            "mech_injection": inj}


def attach_epochs(df: pd.DataFrame, frame_epoch: pd.DataFrame
                  ) -> pd.DataFrame:
    """Add ``mech_epoch`` / ``detector_epoch`` / ``epoch_certain`` columns
    to the frames DataFrame (NULL / 0 where the frame has no epoch)."""
    out = df.merge(
        frame_epoch[["obs_rowid", "mech_epoch", "detector_epoch",
                     "epoch_certain"]], on="obs_rowid", how="left")
    out["epoch_certain"] = out["epoch_certain"].fillna(0).astype(int)
    return out


# ---------------------------------------------------------------------------
# Step 3c — coverage per mechanical epoch (TE.F1 close, DE.F5 key)
# ---------------------------------------------------------------------------
def _setting_class(gain, offset, temp) -> tuple:
    """(gain, offset, set-point) of a science frame, as a grouping key.

    The set-point is rounded to the nearest degree: real set-points are
    integers (-10, 0, -20), and the rounding only absorbs the driver's
    float fuzz.  Unknown values stay ``None`` — a class of their own.
    """
    def _n(v, nd=0):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return None
        return round(float(v), nd)
    g = None if (gain is None or (isinstance(gain, float)
                                   and np.isnan(gain))) else str(gain)
    return (g, _n(offset), _n(temp))


def _setting_label(cls: tuple) -> str:
    """Human spelling of a setting class for a requirement key."""
    g, o, t = cls
    parts = []
    if t is not None:
        parts.append(f"{t:+.0f}C")
    if g is not None:
        parts.append(f"gain {g}")
    if o is not None:
        parts.append(f"offset {o:g}")
    return ", ".join(parts) if parts else "settings unrecorded"


def build_coverage_mech(df: pd.DataFrame,
                        project_of_key: dict[str, frozenset[str]],
                        dw_projects: frozenset[str],
                        ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Coverage and gaps per (era, mechanical epoch) — the hardware truth.

    The era-level matrix (:func:`build_coverage`) asks "does a calibration
    frame with the same header signature exist anywhere in the archive?".
    This one asks the question a reduction actually poses: "does one exist
    that may be APPLIED to these frames?" — which adds three conditions:

    * the boundary rule, :func:`inventory.calib_valid_for`: a flat only
      inside the science frames' own mechanical epoch (and never to or
      from an uncertain-gap night), a dark or bias only inside the same
      detector epoch;
    * the settings key, :func:`inventory.settings_match` (DE.F5): darks
      and biases must share the science frames' gain setting, offset
      setting and cooler set-point (+/- 2 C).  A calibration frame whose
      header does not record a setting the science frame records is
      counted as ``n_unverified`` — visible, and never as a match;
    * science frames are grouped by (era, mech_epoch, certainty) and, for
      darks and biases, by setting class, so a -10 C season and a 0 C
      season of one era are two requirements.

    ``status`` is :func:`inventory.coverage_status` over the VERIFIED
    frames, with one extra value: ``unverified`` when the spec would be met
    only by counting frames whose settings cannot be checked.
    """
    sci = select_science(df)
    cal = df[(df["calib_kind"].notna()) & (df["is_canonical"] == 1)].copy()
    cal["is_master"] = [1 if inv.is_master(bn) else 0
                        for bn in cal["basename"]]
    cal_by_era = {int(e): g for e, g in cal.groupby("era_id")}

    def _projects(frames: pd.DataFrame) -> str:
        hit: set[str] = set()
        for tk in frames["target_key"].dropna().unique():
            hit |= inv.projects_of_target(tk, project_of_key, dw_projects)
        return ",".join(sorted(hit))

    def _tag(mech, det, certain) -> inv.EpochTag:
        none = lambda v: None if (v is None or v != v) else v
        return inv.EpochTag(none(mech), none(det), bool(certain))

    cov_rows, gap_rows = [], []
    # dropna=False: science with NO epoch (unknown camera) is a group too —
    # the group for which nothing can be shown valid.
    for (era_id, mech, certain), sgrp in sci.groupby(
            ["era_id", "mech_epoch", "epoch_certain"], dropna=False):
        era_id = int(era_id)
        mech = None if (mech is None or mech != mech) else mech
        det = sgrp["detector_epoch"].iloc[0]
        sci_tag = _tag(mech, det, certain)
        cgrp = cal_by_era.get(era_id, cal.iloc[0:0])
        nights = sgrp["night"].dropna()
        readout = sgrp["readoutm"].iloc[0]
        camera = (readout if isinstance(readout, str) and readout.strip()
                  else "(blank READOUTM)")
        # Which calibration frames pass the BOUNDARY rule for this group.
        cal_tags = [_tag(m_, d_, c_) for m_, d_, c_ in zip(
            cgrp["mech_epoch"], cgrp["detector_epoch"],
            cgrp["epoch_certain"])]
        ok_flat = np.array([inv.calib_valid_for("flat", sci_tag, t)
                            for t in cal_tags], dtype=bool)
        ok_det = np.array([inv.calib_valid_for("dark", sci_tag, t)
                           for t in cal_tags], dtype=bool)
        kind = cgrp["calib_kind"].to_numpy()
        master = cgrp["is_master"].to_numpy() == 1

        def _emit(req_kind, req_key, frames, n_raw, n_master, n_unver,
                  spec_n, gap_eligible=True):
            status = inv.coverage_status(n_raw, n_master, spec_n)
            if status != "ok" and n_raw + n_unver >= spec_n:
                status = "unverified"
            cov_rows.append({
                "era_id": era_id, "mech_epoch": mech,
                "epoch_certain": int(certain),
                "req_kind": req_kind, "req_key": req_key,
                "n_science": int(len(frames)), "n_calib_raw": int(n_raw),
                "n_calib_master": int(n_master),
                "n_unverified": int(n_unver), "spec_n": spec_n,
                "status": status})
            if status != "ok" and gap_eligible:
                gap_rows.append({
                    "era_id": era_id, "mech_epoch": mech,
                    "epoch_certain": int(certain), "camera": camera,
                    "first_night": nights.min() if len(nights) else None,
                    "last_night": nights.max() if len(nights) else None,
                    "need_kind": req_kind,
                    "spec": inv.gap_spec(req_kind, req_key, n_raw, spec_n),
                    "have_raw": int(n_raw), "have_master": int(n_master),
                    "have_unverified": int(n_unver), "status": status,
                    "n_science_frames_blocked": int(len(frames)),
                    "projects_affected": _projects(frames)})

        # ---- bias and dark: per setting class ----------------------------
        classes = [_setting_class(g, o, t) for g, o, t in zip(
            sgrp["hdr_gain"], sgrp["hdr_offset"], sgrp["set_temp"])]
        cls_series = pd.Series(classes, index=sgrp.index, dtype=object)
        cal_set = list(zip(cgrp["hdr_gain"], cgrp["hdr_offset"],
                           cgrp["set_temp"]))
        for cls in sorted(set(classes), key=repr):
            sub = sgrp[cls_series == cls]
            verdict = np.array(
                [inv.settings_match(cls[0], cls[1], cls[2], g, o, t)
                 for g, o, t in cal_set], dtype=object)
            match = verdict == inv.SETTINGS_MATCH
            unver = verdict == inv.SETTINGS_UNVERIFIED
            label = _setting_label(cls)
            is_bias = (kind == "bias") & ok_det
            _emit("bias", label, sub,
                  (is_bias & match & ~master).sum(),
                  (is_bias & match & master).sum(),
                  (is_bias & unver).sum(), inv.SPEC_N_BIAS)
            sbins = pd.Series([inv.exptime_bin(x) for x in sub["exptime"]],
                              index=sub.index)
            cal_exp = cgrp["exptime"].to_numpy()
            for b in sorted({x for x in sbins if x is not None}):
                ssub = sub[sbins == b]
                emat = np.array([inv.dark_matches(x, b) for x in cal_exp],
                                dtype=bool)
                is_dark = (kind == "dark") & ok_det & emat
                _emit("dark", f"{inv.fmt_exptime(b)} @ {label}", ssub,
                      (is_dark & match & ~master).sum(),
                      (is_dark & match & master).sum(),
                      (is_dark & unver).sum(), inv.SPEC_N_DARK)

        # ---- flat: per filter, inside the mechanical epoch only ----------
        filt = sgrp["filter"].fillna("(blank)")
        cal_filt = cgrp["filter"].to_numpy()
        for f in sorted(filt.unique()):
            sub = sgrp[filt == f]
            is_flat = (kind == "flat") & ok_flat & (cal_filt == f)
            _emit("flat", f, sub, (is_flat & ~master).sum(),
                  (is_flat & master).sum(), 0, inv.SPEC_N_FLAT,
                  gap_eligible=not inv.is_calib_vocab_filter(f))

    cov_cols = ["era_id", "mech_epoch", "epoch_certain", "req_kind",
                "req_key", "n_science", "n_calib_raw", "n_calib_master",
                "n_unverified", "spec_n", "status"]
    gap_cols = ["era_id", "mech_epoch", "epoch_certain", "camera",
                "first_night", "last_night", "need_kind", "spec", "have_raw",
                "have_master", "have_unverified", "status",
                "n_science_frames_blocked", "projects_affected"]
    coverage = pd.DataFrame(cov_rows, columns=cov_cols)
    gaps = pd.DataFrame(gap_rows, columns=gap_cols)
    if len(gaps):
        gaps = gaps.sort_values(
            ["n_science_frames_blocked", "era_id"],
            ascending=[False, True]).reset_index(drop=True)
    return coverage, gaps


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
              "calib_gaps", "mech_epoch", "night_mech_epoch",
              "frame_mech_epoch", "mech_rotation_null", "mech_injection",
              "calib_coverage_mech", "calib_gaps_mech", "s0b_build_meta")
S0_PROTECTED = frozenset({"frames", "aliases", "eras", "project_counts",
                          "build_meta"})


def null_summary(null: pd.DataFrame) -> dict:
    """The false-alarm numbers of the same-target rotation test.

    Computed from ``mech_rotation_null`` and stored in ``s0b_build_meta`` so
    the report can quote them without re-deriving them: how many
    within-epoch same-target differences were measured, their robust
    sigma, their 95th/99th percentile and maximum in absolute value, and
    the fraction exceeding the step threshold.

    CAVEAT, stated where the number is made: this distribution is
    CENSORED.  A night whose MEDIAN same-target difference exceeded the
    threshold became a boundary and left the null; what remains
    under-represents the tail.  The fraction above threshold is therefore a
    LOWER bound on the per-target false-alarm rate, and the injection test,
    not this number, is the evidence for the rule's power.
    """
    if not len(null):
        return {"mech_null_n": "0"}
    d = np.abs(null["delta_deg"].to_numpy(dtype=float))
    signed = null["delta_deg"].to_numpy(dtype=float)
    mad = float(np.median(np.abs(signed - np.median(signed))))
    return {
        "mech_null_n": str(len(d)),
        "mech_null_sigma_deg": f"{1.4826 * mad:.4f}",
        "mech_null_p95_abs_deg": f"{np.percentile(d, 95):.4f}",
        "mech_null_p99_abs_deg": f"{np.percentile(d, 99):.4f}",
        "mech_null_max_abs_deg": f"{d.max():.4f}",
        "mech_null_frac_above_step":
            f"{float((d > inv.MECH_ROT_STEP_DEG).mean()):.5f}",
    }


def write_inventory(manifest_path: Path, links: pd.DataFrame,
                    calib: pd.DataFrame, coverage: pd.DataFrame,
                    gaps: pd.DataFrame, mech: dict | None = None,
                    coverage_mech: pd.DataFrame | None = None,
                    gaps_mech: pd.DataFrame | None = None) -> None:
    """Swap the S0b tables into the manifest, atomically per table set.

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
        {"key": "mech_rot_step_deg", "value": str(inv.MECH_ROT_STEP_DEG)},
        {"key": "mech_rot_coarse_deg",
         "value": str(inv.MECH_ROT_COARSE_DEG)},
        {"key": "mech_rot_max_abs_dec",
         "value": str(inv.MECH_ROT_MAX_ABS_DEC)},
        {"key": "mech_rot_min_frames", "value": str(inv.MECH_ROT_MIN_FRAMES)},
        {"key": "set_temp_tol_c", "value": str(inv.SET_TEMP_TOL_C)},
        {"key": "injection_trials", "value": str(INJECTION_TRIALS)},
        {"key": "injection_seed", "value": str(INJECTION_SEED)},
    ])
    tables = {
        "raw_reduced_links": links,
        "calib_frames": calib,
        "calib_coverage": coverage,
        "calib_gaps": gaps,
    }
    if mech is not None:
        tables.update(mech)
        tables["calib_coverage_mech"] = coverage_mech
        tables["calib_gaps_mech"] = gaps_mech
        meta = pd.concat([meta, pd.DataFrame(
            [{"key": k, "value": v} for k, v in
             null_summary(mech["mech_rotation_null"]).items()])],
            ignore_index=True)
    tables["s0b_build_meta"] = meta
    # closing() is required: sqlite3's own context manager only manages the
    # TRANSACTION — it does NOT close the file handle (the S0 lesson).
    # The generous timeout lets this build coexist with sibling stages
    # writing their own tables: SQLite serializes the writers.
    with closing(sqlite3.connect(manifest_path, timeout=600.0)) as con:
        con.execute("PRAGMA busy_timeout = 600000")
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
        if mech is not None:
            con.execute("CREATE INDEX IF NOT EXISTS ix_fme_rowid "
                        "ON frame_mech_epoch(obs_rowid)")
            con.execute("CREATE INDEX IF NOT EXISTS ix_nme_cam_night "
                        "ON night_mech_epoch(camera, night)")
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
# Step 6 — the mechanical-epoch evidence figure
# ---------------------------------------------------------------------------
def fig_mech_epochs(manifest_path: Path) -> Path:
    """Draw the evidence behind ``mech_epoch`` from the tables just written.

    Three panels, one question each:

    (a) WHAT DID THE SKY DO ON THE DETECTOR?  Every night's median position
        angle against date, one marker style per camera, with each epoch
        boundary drawn and labelled by its cause.  The angle is folded to
        (-90, 270] so that the two orientations the cameras have been
        mounted in (near 0 and near 180) both sit away from the wrap.
    (b) HOW OFTEN DOES THE STEP RULE FIRE ON NOTHING?  The distribution of
        same-target rotation differences INSIDE epochs, with the +/-0.3
        threshold.
    (c) WOULD IT HAVE SEEN A REAL STEP?  Injected-step recovery fraction
        against step size, one line per epoch that could host the test.

    Reads ONLY the manifest tables (the house rule: a figure is a query
    result) and writes ``docs/pipeline/figures/s0b/s0b_mech_epochs.png``.
    """
    import datetime as _dt

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    from macro_core import plotstyle as ps

    with closing(sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True,
                                 timeout=600.0)) as con:
        nights = pd.read_sql_query(
            "SELECT camera, night, rot_median_deg, n_rot, certain "
            "FROM night_mech_epoch WHERE rot_median_deg IS NOT NULL", con)
        epochs = pd.read_sql_query(
            "SELECT * FROM mech_epoch ORDER BY seq", con)
        null = pd.read_sql_query(
            "SELECT delta_deg FROM mech_rotation_null", con)
        inj = pd.read_sql_query("SELECT * FROM mech_injection", con)

    def _date(s):
        return _dt.date.fromisoformat(s)

    def _fold(a):
        """Fold an angle into (-90, 270]."""
        return ((np.asarray(a, dtype=float) + 90.0) % 360.0) - 90.0

    out_dir = REPO_ROOT / "docs" / "pipeline" / "figures" / "s0b"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "s0b_mech_epochs.png"
    cams = sorted(nights["camera"].unique(),
                  key=lambda c: nights[nights["camera"] == c]["night"].min())
    with plt.rc_context(ps.STYLE):
        fig = plt.figure(figsize=(ps.COL_DOUBLE * 1.5, 7.6))
        gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0])
        ax = fig.add_subplot(gs[0, :])
        for i, cam in enumerate(cams):
            sub = nights[nights["camera"] == cam]
            kw = ps.series(i)
            ax.plot([_date(n) for n in sub["night"]],
                    _fold(sub["rot_median_deg"]), linestyle="none",
                    markersize=3.2, label=f"{cam} ({len(sub)} nights)", **kw)
        for _, e in epochs.iterrows():
            ax.axvline(_date(e["first_night"]), color=ps.RULE, linewidth=0.7,
                       linestyle=":" if e["step_marginal"] else "-",
                       alpha=0.7, zorder=0)
        ax.set_ylabel("nightly median position angle (deg, folded)")
        ax.set_title(f"(a) {len(epochs)} mechanical epochs: every vertical "
                     "rule is a boundary (dotted = marginal step)")
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        ax.legend(loc="center left", fontsize=7, ncol=1)

        axn = fig.add_subplot(gs[1, 0])
        if len(null):
            d = null["delta_deg"].to_numpy(dtype=float)
            lim = max(0.6, float(np.abs(d).max()) * 1.05)
            axn.hist(d, bins=np.linspace(-lim, lim, 61), color=ps.ACCENT,
                     edgecolor=ps.PAPER, linewidth=0.3)
            for sgn in (-1, 1):
                axn.axvline(sgn * inv.MECH_ROT_STEP_DEG,
                            **ps.reference_kw(color=ps.BAD))
            frac = float((np.abs(d) > inv.MECH_ROT_STEP_DEG).mean())
            axn.set_yscale("log")
            axn.set_title(f"(b) same-target differences inside epochs: "
                          f"N = {len(d)}, {100 * frac:.1f}% beyond "
                          f"+/-{inv.MECH_ROT_STEP_DEG} deg")
        axn.set_xlabel("rotation difference, night minus epoch (deg)")
        axn.set_ylabel("target-nights")

        axi = fig.add_subplot(gs[1, 1])
        for i, (eid, grp) in enumerate(inj.groupby("mech_epoch", sort=True)):
            grp = grp.sort_values("step_deg")
            kw = ps.line_series(i)
            axi.plot(grp["step_deg"], grp["frac_detected"], markersize=3.5,
                     label=f"{eid} ({int(grp['n_nights'].iloc[0])} n)", **kw)
        axi.axvline(inv.MECH_ROT_STEP_DEG, **ps.reference_kw(color=ps.BAD))
        axi.set_xscale("log")
        axi.set_ylim(-0.03, 1.03)
        axi.set_xlabel("injected rotation step (deg)")
        axi.set_ylabel("fraction recovered")
        axi.set_title("(c) injected-step recovery, per epoch")
        if len(inj):
            axi.legend(fontsize=5.5, ncol=2, loc="center right")
        fig.tight_layout()
        fig.savefig(out, dpi=ps.WEB_DPI)
        plt.close(fig)
    return out


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
    p.add_argument("--skip-figures", action="store_true",
                   help="do not draw the mechanical-epoch evidence figure")
    p.add_argument("--dry-run", action="store_true",
                   help="compute everything and print the summary, but "
                        "write nothing")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.manifest.exists():
        print(f"ERROR: manifest not found: {args.manifest}", file=sys.stderr)
        return 2

    print(f"[S0b] reading manifest {args.manifest} ...")
    with closing(sqlite3.connect(
            f"file:{args.manifest}?mode=ro", uri=True, timeout=600.0)) as con:
        df = load_frames(con)
        project_of_key, dw_projects = load_project_targets(con)
        rot = load_rotations(con)
    print(f"[S0b]   {len(df):,} manifest rows; "
          f"{int(df['calib_kind'].notna().sum()):,} rows classify as "
          "calibration frames")

    print("[S0b] linking reduced tree to raw canonical frames ...")
    links = build_links(df)
    by_method = links["match_method"].value_counts()
    print("[S0b]   " + "; ".join(f"{k}: {v:,}" for k, v in by_method.items()))

    print("[S0b] segmenting mechanical epochs ...")
    mech = build_mech_tables(df, rot)
    ep = mech["mech_epoch"]
    print(f"[S0b]   {len(rot):,} plate-solved rotations; "
          f"{len(mech['night_mech_epoch']):,} camera-nights -> "
          f"{len(ep)} mechanical epochs "
          f"({ep['detector_epoch'].nunique()} detector epochs); "
          f"{int((mech['night_mech_epoch']['certain'] == 0).sum())} "
          "uncertain-gap nights")
    for _, e in ep.iterrows():
        rot_txt = ("no rotation evidence" if pd.isna(e["rotation_deg"])
                   else f"rot {e['rotation_deg']:.2f} deg")
        step_txt = ("" if pd.isna(e["step_deg"])
                    else f", step {e['step_deg']:+.2f}"
                         f"+/-{e['step_err_deg']:.2f}"
                         + (" MARGINAL" if e["step_marginal"] else ""))
        print(f"[S0b]     {e['mech_epoch']:<20} -> {e['last_night']}  "
              f"{e['boundary_cause']:<28} {rot_txt}{step_txt}")
    df = attach_epochs(df, mech["frame_mech_epoch"])

    print("[S0b] building calibration census ...")
    calib = build_calib_frames(df)
    print(f"[S0b]   {len(calib):,} canonical calibration frames "
          f"({int((calib['is_master'] == 1).sum())} masters)")

    print("[S0b] computing era coverage matrix and gaps ...")
    coverage, gaps = build_coverage(df, project_of_key, dw_projects)
    print(f"[S0b]   {len(coverage):,} coverage cells across "
          f"{coverage['era_id'].nunique()} science eras; "
          f"{len(gaps):,} gaps")

    print("[S0b] re-counting coverage per mechanical epoch ...")
    coverage_mech, gaps_mech = build_coverage_mech(
        df, project_of_key, dw_projects)
    by_status = coverage_mech["status"].value_counts()
    print(f"[S0b]   {len(coverage_mech):,} cells; "
          + "; ".join(f"{k}: {v:,}" for k, v in by_status.items())
          + f"; {len(gaps_mech):,} gaps")

    if args.dry_run:
        print("[S0b] --dry-run: nothing written.")
        return 0

    print(f"[S0b] writing inventory tables -> {args.manifest}")
    write_inventory(args.manifest, links, calib, coverage, gaps,
                    mech=mech, coverage_mech=coverage_mech,
                    gaps_mech=gaps_mech)

    if not args.skip_figures:
        print("[S0b] drawing the mechanical-epoch figure ...")
        fig_path = fig_mech_epochs(args.manifest)
        print(f"[S0b] figure -> {fig_path}")

    if not args.skip_report:
        print("[S0b] rendering evidence report ...")
        from macro_core import report_s0b
        report_path = report_s0b.render_report(args.manifest)
        print(f"[S0b] report -> {report_path}")

    print("[S0b] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
