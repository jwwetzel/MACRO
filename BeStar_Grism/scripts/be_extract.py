#!/usr/bin/env python
"""be_extract — BE-S3: per-frame trace and optimal extraction of every Be-star
frame whose Step-0 disposition admits it, on the SHARED grism library.

WHY IT IS THIN
--------------
The committee ruled (U2, SYNTHESIS §4 "Inherit G-1…G-5") that the Be paper
inherits the grism calibration library rather than keeping a second copy of
it.  This script therefore only (1) builds the work list from ``be_frames``
and (2) hands it to ``macro_grism.store.reduce_batch`` — the same function the
T CrB validation uses: HDU-resolving reader, per-frame trace fit (never reused
between frames), Horne extraction with the library's background model, the
measured gain/saturation of the frame's detector (G-2), the detector
package's bad-pixel mask.  Results go to THIS project's own grism-schema
database ``BeStar_Grism/products/be_grism.sqlite`` (table ``g_frames``) and
spectrum cache ``BeStar_Grism/products/spec1d/`` — never into the grism
worker's files.  A frame that fails leaves a row with the error, not a gap.

Frames dispositioned ``extract_if_identity`` (S2c indeterminate) are reduced
too; whether they enter the series is decided by the identity gate in
``be_measure.py`` (BE-S0-dispositions).

USAGE
    /opt/miniconda3/envs/rlmt-checks/bin/python BeStar_Grism/scripts/be_extract.py [--workers 4]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

sys.path.insert(0, str(C.REPO / "pipeline"))
from macro_grism import db as gdb  # noqa: E402
from macro_grism import store  # noqa: E402

BE_GRISM_DB = C.PROJECT / "products" / "be_grism.sqlite"
SPEC_DIR = C.PROJECT / "products" / "spec1d"


def tasks() -> list[dict]:
    con = C.be_db()
    fr = pd.read_sql("""SELECT obs_rowid, path, filter, night, jd, exptime, label, role,
                               canonical_target, pointing_offset_deg
                        FROM be_frames WHERE disposition != 'exclude'
                        """, con)
    return [dict(path=r.path, obs_rowid=int(r.obs_rowid), filter=r.filter, night=r.night,
                 jd=r.jd, exptime=r.exptime, target_best=r.canonical_target, star=r.label,
                 sample=f"be_{r.role}", pointing_offset_deg=r.pointing_offset_deg)
            for r in fr.itertuples()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    con = gdb.connect_grism(BE_GRISM_DB)
    t = tasks()
    print(f"{len(t)} grism frames (HaGrism = hrg unit, OGGrism = lrg unit) to reduce at {store.REDUCE_VERSION}")
    n = store.reduce_batch(con, t, spec_dir=SPEC_DIR, workers=a.workers)
    gdb.set_g_meta(con, "be_reduce_version", store.REDUCE_VERSION)
    con.commit()
    print(f"reduced {n} this call;",
          dict(con.execute("SELECT status, COUNT(*) FROM g_frames GROUP BY status").fetchall()))


if __name__ == "__main__":
    main()
