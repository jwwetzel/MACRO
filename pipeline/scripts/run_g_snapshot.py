#!/usr/bin/env python
"""Snapshot the manifest tables the grism track reads.

WHY
---
This wave runs a dozen work packages in one working tree.  The shared
manifest (600 MB, on a spinning disk) is being read by all of them and
REBUILT by one (``foundation-s0`` swaps it atomically when S0 re-runs).
Two consequences for the grism track:

* a full scan of ``frames`` costs a minute under that contention, and the
  G runners would otherwise issue dozens of them;
* a rebuild in the middle of a run could change which frames are
  canonical between two stages of the same analysis.

So the G runners read a SNAPSHOT: the rows of the manifest that concern
grism work, copied once into ``products/grism/manifest_snapshot.sqlite``
together with the source file's size, mtime and build stamps.  Every
number this package reports is then traceable to one named state of the
manifest, and the integrated rebuild re-runs this script first.

The manifest is opened READ-ONLY (``mode=ro``); nothing is written to it.

WHAT IS COPIED
--------------
* ``frames``            — every row whose FILTER card can be a grism
                          (hrg/lrg/HaGrism/OGGrism/HaG, slot '6', 'W'),
                          plus every row of the named science targets in
                          any filter (T CrB, NGC 5548) — canonical or not.
* ``frame_dispersion``  — whole table (S2c verdicts; F-6 re-issues them).
* ``g_extractions``, ``g_gate_calib``, ``g_build_meta`` — the v1 grism
                          library's results (the "before" of step 7).
* ``calib_frames``, ``eras``, ``detector_params``, ``build_meta``,
  ``s2c_build_meta``, ``mech_epoch`` — small reference tables.
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_grism import config as gconfig                # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
SNAPSHOT = gconfig.PRODUCTS / "manifest_snapshot.sqlite"

#: FILTER cards (lower-case) that can be a grism.
GRISM_CARDS = ("hrg", "lrg", "hagrism", "oggrism", "hag", "6", "w")

#: Science targets copied in every filter (lower(target_best) LIKE).
TARGET_LIKES = ("t crb%", "ngc 5548%", "ngc5548%")

WHOLE_TABLES = ("frame_dispersion", "g_extractions", "g_gate_calib",
                "g_build_meta", "calib_frames", "eras", "detector_params",
                "build_meta", "s2c_build_meta", "mech_epoch")


def snapshot(manifest: Path, out: Path) -> dict:
    """Write the snapshot atomically (temp file + rename) and return the
    row counts per table."""
    tmp = out.with_suffix(".tmp.sqlite")
    if tmp.exists():
        tmp.unlink()
    out.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(str(tmp))
    con.execute(f"ATTACH DATABASE 'file:{manifest}?mode=ro' AS m")
    counts = {}
    cards = ",".join(f"'{c}'" for c in GRISM_CARDS)
    likes = " OR ".join(f"lower(target_best) LIKE '{t}'"
                        for t in TARGET_LIKES)
    con.execute(f"""
        CREATE TABLE frames AS SELECT * FROM m.frames
        WHERE lower(filter) IN ({cards}) OR {likes}""")
    con.execute("CREATE INDEX ix_snap_frames_path ON frames(path)")
    con.execute("CREATE INDEX ix_snap_frames_rowid ON frames(obs_rowid)")
    con.execute("CREATE INDEX ix_snap_frames_tgt ON frames(target_best)")
    counts["frames"] = con.execute(
        "SELECT count(*) FROM frames").fetchone()[0]
    have = {r[0] for r in con.execute(
        "SELECT name FROM m.sqlite_master WHERE type='table'")}
    for t in WHOLE_TABLES:
        if t not in have:
            counts[t] = None                 # absent in this manifest
            continue
        # frame_dispersion's background-map BLOBs are not needed here.
        cols = "*" if t != "frame_dispersion" else ", ".join(
            r[1] for r in con.execute(f"PRAGMA m.table_info({t})")
            if r[1] != "bg_map")
        con.execute(f"CREATE TABLE {t} AS SELECT {cols} FROM m.{t}")
        counts[t] = con.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
    if counts.get("frame_dispersion"):
        con.execute("CREATE INDEX ix_snap_fd ON frame_dispersion(obs_rowid)")
    st = os.stat(manifest)
    con.execute("CREATE TABLE snapshot_meta (key TEXT PRIMARY KEY, "
                "value TEXT)")
    meta = {"source": str(manifest), "source_size": st.st_size,
            "source_mtime_utc": time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)),
            "snapshot_utc": time.strftime(
                "%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    meta.update({f"rows_{k}": v for k, v in counts.items()})
    con.executemany("INSERT INTO snapshot_meta VALUES (?, ?)",
                    [(k, str(v)) for k, v in meta.items()])
    con.commit()
    con.execute("DETACH DATABASE m")
    con.close()
    os.replace(tmp, out)
    return counts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    ap.add_argument("--out", default=str(SNAPSHOT))
    args = ap.parse_args(argv)
    counts = snapshot(Path(args.manifest), Path(args.out))
    for k, v in counts.items():
        print(f"  {k:18s} {v}")
    print(f"snapshot: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
