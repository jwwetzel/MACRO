"""Batch plumbing shared by the G runners: reduce frames in worker
processes, cache the 1-D spectra, record ``g_frames`` rows.

The archive lives on a spinning disk and a frame costs 1-6 s to read and
decompress, so every runner works from a CACHE: a frame is reduced once
(``products/grism/spec1d/<stem>.npz``) and every later stage — line
measurement, zero points, variance checks, equivalent widths — reads the
cached arrays.  The cache key is the archive-relative path plus the sky
method; ``g_frames.code_version`` records which code produced it, and a
version bump re-reduces.

Workers never touch SQLite (one writer: the parent), and never write
anywhere but the cache directory.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

from . import config as gconfig
from . import db as gdb
from .reduce import reduce_frame

#: Bump when a change would alter cached spectra.
REDUCE_VERSION = "G v2.0 (2026-10-03)"

DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
SPEC_DIR = gconfig.PRODUCTS / "spec1d"

#: Arrays cached per spectrum.
SPEC_KEYS = ("flux", "var", "box", "box_var", "bg", "peak", "n_sat",
             "n_rej", "inside", "fwhm_px", "fwhm_x")

#: Manifest columns a runner selects for its worklist.
MANIFEST_COLS = ("obs_rowid, path, filter, night, jd, exptime, "
                 "target_best, pointing_offset_deg")


def spec_path(rel_path: str, sky: str = "poly",
              spec_dir: Path = SPEC_DIR) -> Path:
    """Cache file for one archive frame and sky method."""
    stem = Path(rel_path).name
    for suffix in (".fz", ".fts", ".fits", ".fit"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    night_dir = Path(rel_path).parent.name
    return spec_dir / f"{night_dir}__{stem}__{sky}.npz"


def load_spec(rel_path: str, sky: str = "poly",
              spec_dir: Path = SPEC_DIR) -> Optional[dict]:
    """Cached spectrum arrays as a dict, or None when not cached."""
    p = spec_path(rel_path, sky, spec_dir)
    if not p.exists():
        return None
    with np.load(p) as z:
        return {k: z[k] for k in z.files}


def _work(task: dict) -> dict:
    """Worker: reduce one frame, write its cache file, return the
    ``g_frames`` row (plain dict).  Every failure becomes a row with the
    error text in ``status`` — a bad frame must leave a record, not a
    crashed batch."""
    rel = task["path"]
    row = {"path": rel, "obs_rowid": task.get("obs_rowid"),
           "sample": task.get("sample"), "target": task.get("target_best"),
           "star": task.get("star"), "filter": task.get("filter"),
           "grism": gconfig.grism_unit(task.get("filter")),
           "night": task.get("night"),
           "mech_epoch": gconfig.mech_epoch_id(task["night"])
           if task.get("night") else None,
           "jd": task.get("jd"), "exptime": task.get("exptime"),
           "pointing_offset_deg": task.get("pointing_offset_deg"),
           "code_version": REDUCE_VERSION}
    sky = task.get("sky", "poly")
    try:
        r = reduce_frame(str(Path(task["archive"]) / rel),
                         night=task.get("night"), sky=sky,
                         diagnostics=bool(task.get("diagnostics")))
    except Exception as exc:                         # noqa: BLE001
        row["status"] = f"error: {type(exc).__name__}: {exc}"[:300]
        return row
    h, tr = r["header"], r["trace"]
    row.update(
        status=r["status"], layout=r["layout"],
        ny=int(r["shape"][0]), nx=int(r["shape"][1]),
        detector_key=r["detector"].key,
        instrume=h.get("INSTRUME"), readoutm=h.get("READOUTM"),
        flipstat=h.get("FLIPSTAT"), telpier=h.get("TELPIER"),
        ccd_temp=h.get("CCD-TEMP"), set_temp=h.get("SET-TEMP"),
        focpos=h.get("FOCPOS"), focuspos=h.get("FOCUSPOS"),
        airmass=h.get("AIRMASS"),
        hdr_ra=h.get("OBJCTRA") or h.get("RA"),
        hdr_dec=h.get("OBJCTDEC") or h.get("DEC"),
        trace_height=tr["height"], trace_slope=tr["slope"],
        trace_u=tr["u"], trace_rms_px=tr["rms_px"],
        trace_n=tr["n_centroids"],
        trace_coeffs=json.dumps([float(c) for c in tr["coeffs"]]))
    if row.get("jd") is None and h.get("JD") is not None:
        row["jd"] = h.get("JD")
    if row.get("exptime") is None:
        row["exptime"] = h.get("EXPTIME")
    if tr["extent"] is not None:
        row["trace_x0"], row["trace_x1"] = tr["extent"]
    if r["status"] == "ok":
        out = spec_path(rel, sky, Path(task["spec_dir"]))
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_suffix(".tmp.npz")
        arrays = {k: r["spec"][k] for k in SPEC_KEYS}
        arrays.update(r.get("diag", {}))       # 2-D-derived diagnostics
        np.savez_compressed(tmp, **arrays)
        os.replace(tmp, out)
        row.update(peak_adu=r["peak_adu"], n_sat_cols=r["n_sat_cols"],
                   snr_median=r["snr_median"],
                   fwhm_px=r["fwhm_px_median"],
                   sky_adu=r["sky_median_adu"],
                   spec_file=str(out.relative_to(gconfig.PRODUCTS))
                   if str(out).startswith(str(gconfig.PRODUCTS))
                   else str(out))
    return row


def reduce_batch(con, tasks: Iterable[dict], archive: Path = DEFAULT_ARCHIVE,
                 spec_dir: Path = SPEC_DIR, workers: int = 4,
                 redo: bool = False, log=print) -> int:
    """Reduce every task not already in ``g_frames`` (at the current
    code version), writing rows as results arrive.  Returns the number
    of frames reduced this call.

    A task is a dict with at least ``path`` and ``night``; ``sample`` and
    ``star`` label the row.  An existing row is KEPT (resume) unless
    ``redo`` — but its ``sample``/``star`` labels are refreshed, because
    one frame can serve two samples (a T CrB frame is both science and
    an identity-gate member).
    """
    tasks = list(tasks)
    done = {} if redo else {
        r[0]: (r[1], r[2]) for r in con.execute(
            "SELECT path, code_version, spec_file FROM g_frames")}
    todo = []
    for t in tasks:
        prev = done.get(t["path"])
        if prev is not None and prev[0] == REDUCE_VERSION:
            # Already reduced at this code version.  Redo only when this
            # task needs the 2-D diagnostics and the cache lacks them.
            need = False
            if t.get("diagnostics") and prev[1]:
                cached = load_spec(t["path"], t.get("sky", "poly"),
                                   spec_dir)
                need = cached is not None and "skyvar" not in cached
            if not need:
                # Refresh the labels (one frame can serve two samples).
                if t.get("sample"):
                    con.execute(
                        "UPDATE g_frames SET sample = ?, star = "
                        "coalesce(?, star) WHERE path = ?",
                        (t["sample"], t.get("star"), t["path"]))
                continue
        t = dict(t)
        t["archive"], t["spec_dir"] = str(archive), str(spec_dir)
        todo.append(t)
    con.commit()                 # release the label-refresh transaction
    if not todo:
        return 0
    n = 0
    if workers > 1:
        from multiprocessing import get_context
        with get_context("spawn").Pool(workers) as pool:
            for row in pool.imap_unordered(_work, todo, chunksize=1):
                # One row = one transaction.  Several G runners share
                # this database; a write transaction held open across
                # ten slow frame reads starves the others into
                # "database is locked".
                gdb.upsert(con, "g_frames", row)
                con.commit()
                n += 1
                if n % 10 == 0:
                    log(f"  reduced {n}/{len(todo)}")
    else:
        for t in todo:
            gdb.upsert(con, "g_frames", _work(t))
            con.commit()
            n += 1
            if n % 10 == 0:
                log(f"  reduced {n}/{len(todo)}")
    con.commit()
    return n


def default_manifest() -> Path:
    """The manifest a G runner should read: the grism track's snapshot
    (``run_g_snapshot.py``) when it exists, else the shared manifest.
    Both are opened read-only by the callers."""
    snap = gconfig.PRODUCTS / "manifest_snapshot.sqlite"
    if snap.exists():
        return snap
    return gconfig.REPO_ROOT / "products" / "manifest" / \
        "rlmt-manifest.sqlite"
