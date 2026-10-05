#!/usr/bin/env python
"""be_response — BE-S6: the relative spectral response per (grism, mechanical
state), anchored on Vega/CALSPEC and transferred across the ASI→QHY boundary by
the two stability standards; seasons without a standard get none.

THE CHAIN (strategy §4 Step 6, as amended in §10)
------------------------------------------------
* QHY re-seated (S5) and the Andor lrg ladder night (S1): directly on Vega,
  R(λ) = median over Vega frames of [count rate(λ)] / F_CALSPEC(λ), CALSPEC
  ``alpha_lyr_stis_012`` (Bohlin; MAST CALSPEC, fetched and cached here).
* ASI post-monsoon (S3): R_S3(λ) = R_S5(λ) x [C_S3(λ)/C_S5(λ)] of the SAME
  star (η Hya, θ Vir separately), each count-rate spectrum normalised at
  6563 Å.  The two transfers' disagreement (rms over λ) is published as the
  transfer systematic.
* ASI pre-monsoon (S2, season 1): no standard was observed — no response;
  season 1 is EW only (BE-S8).
Every R is RELATIVE (normalised to 1 at 6563 Å): grey transparency is not
removed, and nothing here is quoted as absolute flux.  Masked before the
ratio: Hα (6520–6610), the O2 B and A bands and the 7200 Å water band.

The EW does not use R (it is normalised by its own local continuum); R is the
product of BE-S6 and the input any pseudo-band flux would need.

Writes ``be_response`` and ``be_response_summary``; the BE-S6 acceptance test
(θ Vir EW vs the per-frame water band depth) is ``be_h2o_test`` in be_series.
"""
from __future__ import annotations

import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

sys.path.insert(0, str(C.REPO / "pipeline"))
from macro_grism import store  # noqa: E402
from macro_grism import wavelength as gw  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_response.py"
CALSPEC_URL = ("https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/"
               "alpha_lyr_stis_012.fits")
CALSPEC = C.PROJECT / "products" / "external" / "calspec" / "alpha_lyr_stis_012.fits"
SPEC_DIR = C.PROJECT / "products" / "spec1d"
LIB_DB = C.REPO / "products" / "grism" / "grism.sqlite"
GRID = {"hrg": np.arange(5900.0, 7400.0, 10.0), "lrg": np.arange(4600.0, 7900.0, 25.0)}
MASKS = ((6520.0, 6610.0), (6860.0, 6940.0), (7160.0, 7330.0), (7580.0, 7700.0))
REF = 6563.0


def calspec():
    from astropy.io import fits
    if not CALSPEC.exists():
        CALSPEC.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(CALSPEC_URL, timeout=120) as r:
            CALSPEC.write_bytes(r.read())
    with fits.open(CALSPEC) as hl:
        d = hl[1].data
        return np.asarray(d["WAVELENGTH"], float), np.asarray(d["FLUX"], float)


def binned(w, f, grid):
    """Median of f in bins centred on grid; NaN where empty."""
    half = np.diff(grid).mean() / 2
    out = np.full(len(grid), np.nan)
    for i, g in enumerate(grid):
        s = (w >= g - half) & (w < g + half) & np.isfinite(f)
        if s.sum() >= 2:
            out[i] = np.median(f[s])
    return out


def masked(grid, y):
    y = y.copy()
    for lo, hi in MASKS:
        y[(grid >= lo) & (grid <= hi)] = np.nan
    ok = np.isfinite(y)
    if ok.sum() > 3:
        y[~ok] = np.interp(grid[~ok], grid[ok], y[ok])
    return y


def star_spectra(con, sols, label, grism, state):
    m = pd.read_sql("""SELECT path, night FROM be_frame_meas
                       WHERE qc='pass' AND label=? AND grism=? AND mech_state=?""",
                    con, params=(label, grism, state))
    meta = pd.read_sql("SELECT path, x_halpha, night, grism FROM be_frame_meas", con).set_index("path")
    ex = pd.read_sql("SELECT path, exptime FROM be_frames", con).set_index("path").exptime
    out = []
    for p in m.path:
        sp = store.load_spec(p, "poly", SPEC_DIR)
        r = meta.loc[p]
        sol = gw.solution_for(sols, grism, r.night)
        if sp is None or sol is None or not np.isfinite(r.x_halpha):
            continue
        w = gw.wavelength_axis(len(sp["flux"]), r.x_halpha, sol)
        o = np.argsort(w)
        c = binned(w[o], np.asarray(sp["flux"], float)[o] / ex.loc[p], GRID[grism])
        k = np.argmin(np.abs(GRID[grism] - REF))
        if np.isfinite(c[k]) and c[k] > 0:
            out.append(c / c[k])
    return np.nanmedian(np.array(out), axis=0) if out else None, len(out)


def main() -> None:
    import argparse
    import sqlite3
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    con = C.be_db()
    sols = gw.load_solutions(sqlite3.connect(f"file:{LIB_DB}?mode=ro", uri=True))
    cw, cf = calspec()
    rows, summ = [], []
    resp = {}
    for grism, state in (("hrg", "S5 QHY re-seated"), ("lrg", "S5 QHY re-seated"), ("lrg", "S1 Andor")):
        vega, n = star_spectra(con, sols, "alf Lyr", grism, state)
        if vega is None:
            summ.append(dict(grism=grism, mech_state=state, method="Vega/CALSPEC", n_frames=0,
                             note="no QC-passing Vega frame"))
            continue
        g = GRID[grism]
        ref = binned(cw, cf, g)
        R = masked(g, vega / (ref / ref[np.argmin(np.abs(g - REF))]))
        resp[(grism, state)] = R
        rows += [dict(grism=grism, mech_state=state, wave=float(x), R=float(y), method="Vega/CALSPEC")
                 for x, y in zip(g, R)]
        summ.append(dict(grism=grism, mech_state=state, method="Vega/CALSPEC", n_frames=n,
                         R_6000=float(np.interp(6000, g, R)), R_7000=float(np.interp(7000, g, R)),
                         transfer_rms=np.nan, note="direct"))
    for grism in ("hrg", "lrg"):
        if (grism, "S5 QHY re-seated") not in resp:
            continue
        g = GRID[grism]
        trans = []
        for star in ("eta Hya", "tet Vir"):
            a, na = star_spectra(con, sols, star, grism, "S3 ASI post-monsoon (flipped)")
            b, nb = star_spectra(con, sols, star, grism, "S5 QHY re-seated")
            if a is not None and b is not None:
                trans.append(masked(g, resp[(grism, "S5 QHY re-seated")] * a / b))
        if not trans:
            continue
        R = np.nanmean(trans, axis=0)
        rms = float(np.sqrt(np.nanmean((trans[0] - trans[1]) ** 2)) / 2) if len(trans) == 2 else np.nan
        rows += [dict(grism=grism, mech_state="S3 ASI post-monsoon (flipped)", wave=float(x), R=float(y),
                      method="transfer via eta Hya + tet Vir") for x, y in zip(g, R)]
        summ.append(dict(grism=grism, mech_state="S3 ASI post-monsoon (flipped)",
                         method="transfer via eta Hya + tet Vir", n_frames=len(trans),
                         R_6000=float(np.interp(6000, g, R)), R_7000=float(np.interp(7000, g, R)),
                         transfer_rms=rms, note="two transfer stars; their half-difference is the rms"))
    for grism in ("hrg", "lrg"):
        summ.append(dict(grism=grism, mech_state="S2 ASI pre-monsoon", method="none", n_frames=0,
                         note="season 1: no standard observed — EW only (BE-S8)"))
    pd.DataFrame(rows).to_sql("be_response", con, if_exists="replace", index=False)
    sm = pd.DataFrame(summ)
    sm.to_sql("be_response_summary", con, if_exists="replace", index=False)
    con.commit()
    C.write_md(C.NOTES / "results" / "response.md",
               "Relative response per (grism, mechanical state) (BE-S6)",
               C.md_table(sm, 3) + "\nR normalised to 1 at 6563 Å; CALSPEC alpha_lyr_stis_012 (MAST). "
               "Relative only: no absolute flux is quoted anywhere in the paper.\n", SCRIPT)
    print(sm.to_string())


if __name__ == "__main__":
    main()
