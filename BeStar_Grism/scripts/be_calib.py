#!/usr/bin/env python
"""be_calib — BE-S2: the calibration each Be-star frame actually receives, stated
from the archive's own inventory, per mechanical state.

WHY
---
The strategy's Step 2 asked the observatory for biases, darks and flat pairs at
the October 2026 re-opening.  No frame can be taken before then, and the
camera now mounted (QHY600) has no calibration frame of any kind in the
archive; no grism (hrg/lrg) flat exists for ANY camera.  The chair's ruling:
state the method, apply what the archive allows, and move the acquisition to
the 2027 backlog (BE-X3-calib-acquisition).  This script emits the evidence:
the calibration inventory per (mechanical state, kind) from the manifest's
``calib_frames`` (read-only), the measured detector parameters each state uses
(``s2_camera_configs``), and the label each state's error floor carries.

THE METHOD (what the paper says, and why each choice is forced or chosen)
------------------------------------------------------------------------
* Bias/dark: where the archive holds darks of the same mechanical state (ASI,
  both states) the grism library's master-dark arm is available; everywhere
  the per-column FLANKING-BAND background is the primary subtraction, because
  it removes bias, dark current and sky together at the trace and is the only
  option for the QHY.  The two arms are compared where both exist (G-4).
* Flat field: none.  No grism flat exists, and a slitless EW is a ratio of the
  line window to its own adjacent continuum along one trace: pixel-response
  structure enters only through its variation between the line and continuum
  windows.  That residual is not measured on archival data; the dither test
  that would measure it is BE-X1 (2027 backlog), and the paper says so.
* Gain, read noise, saturation: measured, per configuration (F-4/F-5,
  ``s2_camera_configs`` / ``s2_linearity_caps``); never the header EGAIN.
* Error-floor label per state: a state with standards nights (>= 2025-12-05)
  gets its per-night floor from the standards' night-to-night scatter (BE-S10);
  a state without standards is LOWER-BOUND-ONLY (intra-night scatter).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_calib.py"


def main() -> None:
    import argparse
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    con = C.be_db()
    era = pd.read_sql("SELECT DISTINCT mech_state, mech_epoch, camera, gain_config, gain_e_per_adu, "
                      "gain_err, gain_status, read_noise_e, clip_adu FROM be_era_table", con)
    std = pd.read_sql("""SELECT mech_state, COUNT(DISTINCT night) AS std_nights FROM be_frames
                         WHERE role='standard' AND disposition!='exclude' AND night >= ?
                         GROUP BY mech_state""", con, params=(C.STANDARDS_EPOCH,))
    man = C.manifest_ro()
    cal = pd.read_sql("""SELECT mech_epoch, kind, COUNT(*) AS n,
                         SUM(filter IN ('hrg','lrg','HaGrism','OGGrism')) AS n_grism,
                         SUM(is_master) AS n_master, MIN(night) AS first, MAX(night) AS last
                         FROM calib_frames GROUP BY mech_epoch, kind""", man)
    man.close()
    inv = era[["mech_state", "mech_epoch"]].merge(cal, on="mech_epoch", how="left")
    piv = (inv.pivot_table(index=["mech_state", "mech_epoch"], columns="kind", values="n",
                           aggfunc="sum", fill_value=0).reset_index())
    piv = (era[["mech_state", "mech_epoch"]].drop_duplicates()
           .merge(piv, on=["mech_state", "mech_epoch"], how="left").fillna(0))
    for k in ("bias", "dark", "flat"):
        if k not in piv:
            piv[k] = 0
    gflat = inv[inv.kind == "flat"].groupby("mech_epoch").n_grism.sum()
    piv["grism_flats"] = piv.mech_epoch.map(gflat).fillna(0).astype(int)
    for k in ("bias", "dark", "flat"):
        piv[k] = piv[k].astype(int)
    piv = piv.merge(era[["mech_epoch", "gain_config", "gain_e_per_adu", "gain_err", "gain_status"]]
                    .drop_duplicates("mech_epoch"), on="mech_epoch", how="left")
    st = piv.groupby("mech_state").agg(bias=("bias", "sum"), dark=("dark", "sum"), flat=("flat", "sum"),
                                       grism_flats=("grism_flats", "sum")).reset_index()
    st = st.merge(std, on="mech_state", how="left").fillna({"std_nights": 0})
    st["dark_subtraction"] = ["master-dark arm available + flanking band" if d > 0 else "flanking band only"
                              for d in st.dark]
    st["flat_field"] = "none (no hrg/lrg flat; BE-X1 dither test deferred)"
    st["error_floor_label"] = ["standards night-to-night (BE-S10)" if n > 0 else "LOWER-BOUND-ONLY (intra-night)"
                               for n in st.std_nights]
    st.to_sql("be_calib_states", con, if_exists="replace", index=False)
    piv.to_sql("be_calib_inventory", con, if_exists="replace", index=False)
    con.commit()
    body = (C.md_table(st) + "\n#### Per F-3 epoch (archive `calib_frames`, any filter)\n\n" + C.md_table(piv, 3) +
            "\nMethod (stated in the paper): flanking-band background everywhere, master-dark arm where the state "
            "has darks; no flat-field division (no hrg/lrg flat exists for any camera; the Andor's HaGrism/OGGrism "
            "flats are not applied, for one method across states — the dither test that would bound the pixel-response residual is BE-X1, 2027 backlog); gain/read noise/saturation from the "
            "measured `s2_camera_configs`/`s2_linearity_caps`. The QHY600 now mounted has no bias, dark or flat "
            "in the archive; its acquisition is BE-X3-calib-acquisition (2027 backlog). Era B (ASI) before the "
            "standards epoch carries LOWER-BOUND-ONLY errors; ASI after the monsoon and the QHY take their floor "
            "from the standards (BE-S10).\n")
    C.write_md(C.NOTES / "calibration" / "calibration.md",
               "Calibration from the archive, per mechanical state (BE-S2)", body, SCRIPT)
    print(st.to_string())


if __name__ == "__main__":
    main()
