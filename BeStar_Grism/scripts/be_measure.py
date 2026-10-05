#!/usr/bin/env python
"""be_measure — per-frame measurement of every reduced Be-star spectrum:
wavelength zero point on the FIXED library dispersion (BE-S4), delivered
resolution (BE-S5), the pixel identity test (BE-S0-dispositions), the QC gates
(BE-S1), the telluric H2O band depth (BE-S6 regressor) and the Hα equivalent
width on the paper-wide windows (BE-S7).

INPUTS
  be_grism.sqlite g_frames + spec1d/*.npz   (be_extract.py, shared library)
  products/grism/grism.sqlite g_dispersion  (G-1 fixed solutions, READ-ONLY)
  bestar.sqlite be_frames                   (role, airmass, header cards)

THE MEASUREMENTS, AND WHY EACH IS DONE THIS WAY
----------------------------------------------
Wavelength (U2): the dispersion is the library's one solution per (grism,
  mechanical epoch); per frame only the zero point is fitted — the pixel of Hα.
  Emission stars: the library's emission-line centroid.  Absorption stars
  (standards, null stars, and any science frame where no emission peak is
  found): the deepest absorption feature whose predicted O2-B edge (from the
  fixed solution) is found within 4% of where it should be.  The O2-B
  separation check is stored for EVERY frame; its scatter per (grism, epoch) is
  the acceptance statistic of BE-S4.
Resolution (G-5 method): for a slitless spectrum the line-spread function is
  the image of the star along the dispersion, so LSF FWHM = cross-dispersion
  FWHM (px) x local dispersion at Hα (Å/px), per frame.
Identity (pixel test, never the header alone): ACCEPT when an Hα anchor is
  found and, where the O2-B edge is measurable, it lies within 4% of the fixed
  solution's prediction (a wrong trace or a wrong line fails that); REJECT
  otherwise — and REJECT a strong BeSS emitter (season-median BeSS EW < -5 Å)
  whose Hα is found only in absorption: the brightest trace is another star.  Indeterminate-S2c frames enter the series only on ACCEPT.
EW (strategy §4 Step 7): line 6520–6610 Å, continuum = straight line fitted
  (inverse-variance) through 6480–6520 and 6620–6660 Å; EW = Σ(1 − F/Fc)Δλ,
  positive = absorption.  Error = propagated pixel variance plus the
  continuum-fit covariance.  Columns with a pixel at or above the measured
  saturation cap inside the line window are counted (saturation is judged at
  the target, per frame, in native pixels — standing rule 4).
H2O (PH memo, BE-S6): depth of the 7160–7320 Å water band against a straight
  continuum through 7100–7150 and 7330–7360 Å, where the spectrum covers it.
Continuum (standing rule 6): continuum count rate at 6563 Å (ADU/s) — the
  input of the relative continuum light curve in be_series.py.

QC GATES (BE-S1, strategy §4 Step 1; each logged with its reason)
  no_trace | no_solution | no_anchor | identity_reject | saturated_in_window |
  low_snr (integrated S/N in the line window < 100) | airmass (> 2.5) |
  wide_psf (> 1.5 x star-night median) | focus_off (|FOCPOS - night median| > 100)

OUTPUT: bestar.sqlite ``be_frame_meas`` (one row per reduced frame).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

sys.path.insert(0, str(C.REPO / "pipeline"))
from macro_grism import linecal as lc  # noqa: E402
from macro_grism import store  # noqa: E402
from macro_grism import wavelength as gw  # noqa: E402

BE_GRISM_DB = C.PROJECT / "products" / "be_grism.sqlite"
SPEC_DIR = C.PROJECT / "products" / "spec1d"
LIB_DB = C.REPO / "products" / "grism" / "grism.sqlite"

WIN_LINE = (6520.0, 6610.0)
WIN_CONT = ((6480.0, 6520.0), (6620.0, 6660.0))
H2O_BAND = (7160.0, 7320.0)
H2O_CONT = ((7100.0, 7150.0), (7330.0, 7360.0))
H2O_VALID = (-0.05, 0.5)  # physical range of the band depth (declared, not tuned)
O2B_TOL = 0.04
SNR_MIN = 100.0
AIRMASS_MAX = 2.5
PSF_FACTOR = 1.5
FOCUS_TOL = 100.0
STRONG_EW = -5.0          # Å, BeSS season-median EW of a strong emitter


def line_fit(w, f, v, wins):
    """Inverse-variance straight line through the continuum windows.
    Returns (coef, cov) or None."""
    m = np.zeros_like(w, bool)
    for lo, hi in wins:
        m |= (w >= lo) & (w <= hi)
    m &= np.isfinite(f) & np.isfinite(v) & (v > 0)
    if m.sum() < 8:
        return None
    for lo, hi in wins:                          # both flanks must be present
        if not np.any(m & (w >= lo) & (w <= hi)):
            return None
    x = w[m] - 6563.0
    A = np.column_stack([np.ones_like(x), x])
    W = 1.0 / v[m]
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    coef = cov @ (A.T @ (W * f[m]))
    # scale covariance by the reduced chi2 of the flanks (the flanks carry
    # stellar and telluric structure beyond photon noise)
    r = f[m] - A @ coef
    chi = float(np.sum(r ** 2 * W) / max(m.sum() - 2, 1))
    return coef, cov * max(chi, 1e-9)


def ew_measure(w, f, v):
    fit = line_fit(w, f, v, WIN_CONT)
    lm = (w >= WIN_LINE[0]) & (w <= WIN_LINE[1]) & np.isfinite(f) & np.isfinite(v)
    if fit is None or lm.sum() < 20:
        return None
    coef, cov = fit
    x = w[lm] - 6563.0
    fc = coef[0] + coef[1] * x
    if np.any(fc <= 0):
        return None
    dl = np.abs(np.gradient(w))[lm]
    ew = float(np.sum((1 - f[lm] / fc) * dl))
    # d EW / d coef: sum f/fc^2 * dl * [1, x]
    g = np.array([np.sum(f[lm] / fc ** 2 * dl), np.sum(f[lm] / fc ** 2 * dl * x)])
    var = float(np.sum(v[lm] / fc ** 2 * dl ** 2) + g @ cov @ g)
    snr = float(np.sum(f[lm]) / np.sqrt(np.sum(v[lm])))
    return dict(ew=ew, ew_err=float(np.sqrt(var)), cont_6563=float(coef[0]), snr_win=snr)


def vr_measure(w, f, cont, disp, lsf_a):
    """V/R of a double-peaked emission profile (BE-VR-hold; only used if the
    delivered LSF admits it): the profile, normalised by the continuum level
    at 6563 Å and smoothed by 1 px, is searched for its two highest local
    maxima within 6550–6576 Å.  They count as a resolved pair only when they
    are separated by > 2 px and the minimum between them dips >= 1% of the
    continuum below the fainter peak, AND they are farther apart than this
    frame's delivered LSF FWHM (closer peaks are not resolved); V/R = F_V/F_R,
    the ratio of the normalised peak intensities.  A single
    peak returns vr = NaN with ``vr_status='single'`` — never a value."""
    from scipy.signal import argrelextrema
    LSF_REF_A = lsf_a if np.isfinite(lsf_a) else np.inf
    sel = (w >= 6550.0) & (w <= 6576.0)
    if sel.sum() < 10 or not cont > 0:
        return {"vr_status": "no_data"}
    x, y = w[sel], np.convolve(f[sel] / cont, np.ones(3) / 3, mode="same")
    mx = argrelextrema(y, np.greater)[0]
    mx = mx[(mx > 0) & (mx < len(y) - 1)]
    if len(mx) < 2:
        return {"vr_status": "single"}
    top = mx[np.argsort(y[mx])[-2:]]
    iv, ir = sorted(top)
    dip = y[iv:ir + 1].min()
    if ir - iv <= 2 or min(y[iv], y[ir]) - dip < 0.01 or min(y[iv], y[ir]) <= 1.0:
        return {"vr_status": "single"}
    sep = float(x[ir] - x[iv])
    if not sep > LSF_REF_A:
        # peaks closer than the delivered LSF are not resolved: no value
        return {"vr_status": "unresolved", "peak_sep_a": sep}
    return {"vr_status": "double", "vr": float(y[iv] / y[ir]), "peak_sep_a": sep}


def o2b_edge_fwhm(w, f, edge_wave):
    """BE-S5 acceptance quantity: the O2-B blue edge's 20–80% width in Å,
    converted to a Gaussian-equivalent FWHM (a step convolved with a Gaussian
    rises 20→80% over 0.7153 FWHM).  The band head is not a perfect step, so
    this is an upper-bound LSF proxy; it is reported beside the G-5 library
    method (cross-dispersion FWHM x dispersion), never instead of it."""
    blue = (w >= edge_wave - 40) & (w <= edge_wave - 8)
    red = (w >= edge_wave) & (w <= edge_wave + 15)
    if blue.sum() < 5 or red.sum() < 5:
        return np.nan
    top = np.median(f[blue])
    floor = np.min(f[red])
    if not top > floor:
        return np.nan
    seg = (w >= edge_wave - 10) & (w <= edge_wave + 6)
    x, y = w[seg], (f[seg] - floor) / (top - floor)
    try:
        x80 = x[np.where(y >= 0.8)[0][-1]]          # last point still above 80%
        x20 = x[np.where(y <= 0.2)[0][0]]           # first point below 20%
    except IndexError:
        return np.nan
    width = x20 - x80
    return float(width / 0.7153) if width > 0 else np.nan


def h2o_depth(w, f, v):
    fit = line_fit(w, f, v, H2O_CONT)
    bm = (w >= H2O_BAND[0]) & (w <= H2O_BAND[1]) & np.isfinite(f)
    if fit is None or bm.sum() < 10:
        return np.nan
    coef = fit[0]
    fc = coef[0] + coef[1] * (w[bm] - 6563.0)
    d = float(1.0 - np.sum(f[bm]) / np.sum(fc))
    # A band depth outside [H2O_VALID] is not a water measurement (the band
    # or a flank runs off the trace end, or the continuum fit failed): NaN.
    return d if H2O_VALID[0] <= d <= H2O_VALID[1] else np.nan


def absorption_anchor(flux, sol):
    """Deepest absorption candidate whose O2-B edge is where the fixed
    solution predicts it."""
    norm, _, _ = lc.normalize(flux)
    best = None
    for c in lc.absorption_candidates(norm, n=6):
        chk = gw.check_o2b(flux, c["x"], sol)
        if chk is None or abs(chk["frac_dev"]) > O2B_TOL:
            continue
        if best is None or c["depth"] > best[0]["depth"]:
            best = (c, chk)
    return best


def measure_one(row, sols, meta):
    out = dict(path=row.path, obs_rowid=row.obs_rowid)
    if row.status != "ok":
        out["qc"] = "no_trace"
        return out
    sp = store.load_spec(row.path, "poly", SPEC_DIR)
    sol = gw.solution_for(sols, row.grism, row.night)
    if sp is None or sol is None:
        out["qc"] = "no_solution" if sp is not None else "no_spectrum"
        return out
    flux, var = np.asarray(sp["flux"], float), np.asarray(sp["var"], float)
    nx = len(flux)
    anchor, x_ha, chk = None, None, None
    if meta.role == "science":
        zp = gw.halpha_zero_point(flux, var)
        if zp is not None and zp["snr"] >= gw.PEAK_MIN_SNR:
            anchor, x_ha = "emission", zp["x"]
            chk = gw.check_o2b(flux, x_ha, sol)
    if anchor is None:
        ab = absorption_anchor(flux, sol)
        if ab is not None:
            anchor, x_ha, chk = "absorption", ab[0]["x"], ab[1]
    out.update(anchor=anchor, x_halpha=x_ha, sol_epoch=sol.mech_epoch,
               disp_halpha=sol.local_disp(x_ha) if x_ha is not None else np.nan)
    if chk is not None:
        out.update(o2b_frac_dev=chk["frac_dev"], o2b_dwave=chk["dwave"], o2b_depth=chk["depth"])
    if x_ha is None:
        out["qc"] = "no_anchor"
        out["identity"] = "REJECT"
        return out
    out["identity"] = "REJECT" if (chk is not None and abs(chk["frac_dev"]) > O2B_TOL) else "ACCEPT"
    # A strong BeSS-verified emitter whose spectrum shows Hα in ABSORPTION is
    # the wrong star in the trace (e.g. a brighter Pleiad beside Pleione): the
    # pixels contradict the identity.  Weak emitters (BeSS EW > STRONG_EW) can
    # legitimately show an absorption-dominated profile and are not tested.
    if (meta.role == "science" and anchor == "absorption"
            and np.isfinite(meta.bess_ew_med) and meta.bess_ew_med < STRONG_EW):
        out["identity"] = "REJECT"
        out["identity_reason"] = "absorption_where_BeSS_shows_strong_emission"
    w = gw.wavelength_axis(nx, x_ha, sol)
    order = np.argsort(w)
    w, f, v = w[order], flux[order], var[order]
    nsat = np.asarray(sp["n_sat"], float)[order] if "n_sat" in sp else np.zeros(nx)
    m = ew_measure(w, f, v)
    if m:
        out.update(m)
        out["cont_rate"] = m["cont_6563"] / row.exptime if row.exptime else np.nan
    if m and anchor == "emission":
        out.update(vr_measure(w, f, m["cont_6563"], out["disp_halpha"],
                              (row.fwhm_px or np.nan) * out["disp_halpha"]))
    lm = (w >= WIN_LINE[0]) & (w <= WIN_LINE[1])
    out["n_sat_win"] = int(np.sum(nsat[lm] > 0))
    out["h2o_depth"] = h2o_depth(w, f, v)
    if chk is not None:
        out["o2b_fwhm_a"] = o2b_edge_fwhm(w, f, float(chk["wave"]))
    out["lsf_fwhm_a"] = (row.fwhm_px or np.nan) * out["disp_halpha"]
    return out


def main() -> None:
    import argparse
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    import sqlite3
    gcon = sqlite3.connect(f"file:{BE_GRISM_DB}?mode=ro", uri=True, timeout=600)
    g = pd.read_sql("""SELECT path, obs_rowid, grism, night, exptime, status, fwhm_px, peak_adu,
                              sat_cap_adu, n_sat_cols, snr_median, trace_height, hot_mask
                       FROM g_frames""", gcon)
    lib = sqlite3.connect(f"file:{LIB_DB}?mode=ro", uri=True, timeout=600)
    sols = gw.load_solutions(lib)
    lib_meta = dict(lib.execute("SELECT key, substr(value,1,80) FROM g_meta").fetchall())
    con = C.be_db()
    fr = pd.read_sql("""SELECT obs_rowid, main_id, label, role, filter, night, season, mech_state,
                               mech_epoch, airmass_calc, focpos, ccd_temp, bjd_tdb, disposition,
                               standards_epoch FROM be_frames""", con)
    fr = fr.merge(pd.read_sql("SELECT main_id, bess_ew_med FROM be_sample", con), on="main_id",
                  how="left").set_index("obs_rowid")
    rows = []
    for r in g.itertuples():
        meta = fr.loc[r.obs_rowid]
        d = measure_one(r, sols, meta)
        d.update(peak_adu=r.peak_adu, sat_cap_adu=r.sat_cap_adu, fwhm_px=r.fwhm_px,
                 trace_height=r.trace_height, grism=r.grism)
        rows.append(d)
    m = pd.DataFrame(rows).merge(fr.reset_index(), on="obs_rowid", how="left")
    for c in ("focpos",):
        m[c] = pd.to_numeric(m[c], errors="coerce")

    # --- QC gates (first failing gate is the logged reason) ---------------
    key = ["main_id", "grism", "night"]
    m["psf_night_med"] = m.groupby(key).fwhm_px.transform("median")
    m["focus_night_med"] = m.groupby(key).focpos.transform("median")
    reasons = []
    for r in m.itertuples():
        q = getattr(r, "qc", None)
        if isinstance(q, str):
            reasons.append(q)
        elif r.identity != "ACCEPT":
            reasons.append("identity_reject")
        elif not np.isfinite(getattr(r, "ew", np.nan)):
            reasons.append("ew_window_off_spectrum")
        elif r.n_sat_win > 0:
            reasons.append("saturated_in_window")
        elif r.snr_win < SNR_MIN:
            reasons.append("low_snr")
        elif not (r.airmass_calc < AIRMASS_MAX):
            reasons.append("airmass")
        elif r.fwhm_px > PSF_FACTOR * r.psf_night_med:
            reasons.append("wide_psf")
        elif np.isfinite(r.focpos) and abs(r.focpos - r.focus_night_med) > FOCUS_TOL:
            reasons.append("focus_off")
        else:
            reasons.append("pass")
    m["qc"] = reasons
    m["peak_frac_cap"] = m.peak_adu / m.sat_cap_adu
    m = m.drop(columns=["psf_night_med", "focus_night_med"])
    m.to_sql("be_frame_meas", con, if_exists="replace", index=False)
    for k, v in (("library_g_meta_keys", ",".join(lib_meta)),
                 ("library_solutions", ";".join(f"{a}|{b}" for a, b in sorted(sols))),
                 ("reduce_version", store.REDUCE_VERSION)):
        con.execute("INSERT OR REPLACE INTO be_meta VALUES (?,?)", (k, v))
    con.commit()
    print(m.qc.value_counts().to_string())
    print(m.groupby(["role", "anchor"], dropna=False).size().to_string())


if __name__ == "__main__":
    main()
