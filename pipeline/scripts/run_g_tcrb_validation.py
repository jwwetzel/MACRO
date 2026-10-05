#!/usr/bin/env python
"""G-track: the corrected grism library applied to the T CrB series.

This runner replaces the 2026-08-18 validation pass.  That pass solved a
dispersion per frame, used a header gain four times too small, masked
every pixel above 16.3 kADU and rejected frames on a stale header
pointing; the committee review of 2026-10-03 withdrew all four (OA.E1,
OA.E2, DE.F1, TE.F3, PH.P7).  What runs here is the corrected library on
EVERY T CrB hrg/lrg frame, from the 1-D spectrum cache written by
``run_g_reduce.py --sample tcrb`` (and the calibrator cache of
``run_g_dispersion.py``) — no stage below reads the archive.

STAGES (each independent; ``--all`` chains them in order)
---------------------------------------------------------
``--gate``      (runs first) G-3 / TCRB-A0: the pixel identity gate and
                the Halpha identification it implies   -> g_identity
``--zero``      TCRB-A3: per-frame zero point (Halpha emission) under the
                FIXED detector-coordinate dispersion of the frame's
                (grism, mechanical epoch); the O2-B edge as the check (G-1:
                its wavelength constant across the series) -> g_zero_point
``--variance``  G-2: measured vs predicted sky variance in the trace
                flanks; pooled photon-transfer fit          -> g_variance_check
``--null``      G-4: empty-aperture residual of each sky model, negative-
                continuum census, sky-method and boxcar-vs-optimal
                differences                                 -> g_null_sky
``--lsf``       G-5: delivered line-spread function per frame (T CrB and
                calibrators), focus offset and CCD-TEMP     -> g_lsf
``--sat``       TCRB-A6: saturation triage against the measured cap with
                the bad-pixel mask applied                  -> g_saturation
``--ew``        Halpha equivalent width, corrected library ('v2') and the
                August library's spectra and dispersions ('v1'), same
                estimator                                    -> g_ew
``--summary``   print every acceptance number (the report renders the
                same numbers from the same tables)
``--report``    render docs/pipeline/g_grism.html

Reads ``products/grism/grism.sqlite`` (own tables), the spectrum cache,
and — for the v1 comparison only — the August ``g_extractions`` rows in
the manifest snapshot and their FITS products under
``products/grism/spectra`` (read-only).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from typing import Optional
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_grism import config as gconfig                # noqa: E402
from macro_grism import db as gdb                        # noqa: E402
from macro_grism import ew as gew                        # noqa: E402
from macro_grism import extract as gext                  # noqa: E402
from macro_grism import gate as ggate                    # noqa: E402
from macro_grism import store as gstore                  # noqa: E402
from macro_grism.summary_g import ptc_fit                # noqa: E402
from macro_grism import wavelength as gwave              # noqa: E402

#: Acceptance thresholds, quoted from SYNTHESIS section 3 and the ledger.
G1_SEP_TOL = 0.02          # O2-B edge constant to 1-2% of the lever arm
G2_VAR_TOL = 0.20          # predicted vs measured variance within 20%
G4_METHOD_TOL = 0.03       # sky-method / extraction-method difference < 3%

#: Nominal hrg-minus-lrg focuser offset (counts) under MaxIm's filter-
#: offset table (TE.F4), and how far from nominal a frame may sit before
#: it is flagged off-nominal.  (Under pyscope the offsets are zero —
#: ops package — but the T CrB series is entirely MaxIm-era.)
NOMINAL_FOCUS_OFFSET = {"hrg": 650.0, "lrg": 0.0}
FOCUS_TOL = 150.0

#: Sensor temperature (C) above which a frame is "warm" relative to the
#: -10 C master darks and bad-pixel masks (DE.F5).
WARM_CCD_TEMP = -5.0

#: A column is "negative continuum" when its optimal flux is below zero by
#: more than this many sigma (inside the trace extent, Halpha excluded).
NEG_SIGMA = 3.0

#: Saturation triage: a saturated column within this many A of Halpha
#: (the EW window plus its continuum bands) DISCARDS the frame for EW
#: work; anywhere else on the trace it only FLAGS it.
SAT_DISCARD_HALFWIDTH_A = 65.0

C_KMS = 299792.458


def tcrb_rows(con, extra: str = "") -> list[dict]:
    cur = con.execute(f"""
        SELECT path, obs_rowid, sample, target, grism, night, mech_epoch,
               jd, exptime, detector_key, focpos, ccd_temp, airmass,
               pointing_offset_deg, snr_median, peak_adu, n_sat_cols,
               fwhm_px, sky_adu, trace_x0, trace_x1, sat_cap_adu, hot_mask,
               status
        FROM g_frames WHERE sample = 'tcrb' {extra}
        ORDER BY night, path""")
    cols = [c[0] for c in cur.description]
    return [dict(zip(cols, r)) for r in cur]


def science_flux(spec: dict, key: str = "flux") -> np.ndarray:
    """Flux inside the trace extent (NaN elsewhere)."""
    return np.where(spec["inside"], spec[key], np.nan)


def mad_sigma(v) -> float:
    v = np.asarray(v, dtype=float)
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


# ---------------------------------------------------------------------------
# --zero (TCRB-A3)
# ---------------------------------------------------------------------------
def run_zero(con) -> None:
    sols = gwave.load_solutions(con)
    chosen = dict(con.execute("SELECT path, x_halpha FROM g_identity "
                              "WHERE x_halpha IS NOT NULL"))
    out = []
    for r in tcrb_rows(con, "AND status = 'ok'"):
        sol = gwave.solution_for(sols, r["grism"], r["night"])
        spec = gstore.load_spec(r["path"])
        row = {"path": r["path"], "sample": r["sample"],
               "target": r["target"], "grism": r["grism"],
               "mech_epoch": r["mech_epoch"], "night": r["night"],
               "anchor": "none", "status": "no_solution"}
        if sol is None or spec is None:
            out.append(row)
            continue
        flux = science_flux(spec)
        peak = None
        if chosen.get(r["path"]) is not None:
            # The Halpha identification of the gate (the candidate that
            # matches T CrB's spectrum), re-centroided here.
            for c in gwave.emission_candidates(flux):
                if abs(c["x"] - chosen[r["path"]]) < 2.0:
                    peak = c
        z = gwave.halpha_zero_point(flux, spec["var"], peak=peak)
        if z is None:
            row["status"] = "no_halpha"
            out.append(row)
            continue
        row.update(anchor="halpha_em", x_halpha=z["x"],
                   x_halpha_err=z["x_err"], halpha_snr=z["snr"],
                   halpha_fwhm_px=z["fwhm_px"],
                   halpha_amp_frac=z["amp_frac"],
                   disp_at_halpha=sol.local_disp(z["x"]), status="ok")
        o = gwave.check_o2b(flux, z["x"], sol)
        if o is None:
            row["status"] = "ok_no_o2b"
        else:
            row.update(x_o2b=o["x"], x_o2b_err=o["x_err"],
                       o2b_depth=o["depth"], o2b_wave=o["wave"],
                       o2b_wave_err=o["wave_err"], o2b_dwave=o["dwave"],
                       sep_px=o["sep_px"], sep_pred_px=o["sep_pred_px"],
                       sep_frac_dev=o["frac_dev"])
        out.append(row)
    con.execute("DELETE FROM g_zero_point WHERE sample = 'tcrb'")
    for row in out:
        gdb.upsert(con, "g_zero_point", row)
    con.commit()
    print(f"zero points: {len(out)} frames")


# ---------------------------------------------------------------------------
# --gate (G-3, TCRB-A0)
# ---------------------------------------------------------------------------
def _candidate_spectra(r: dict, sols: dict):
    """Every Halpha hypothesis of one frame: for each of its strongest
    sharp emission features (gwave.emission_candidates), the zero point
    and the frame's spectrum on the fingerprint grid under the FIXED
    dispersion of its (grism, epoch) — raw-normalised and high-passed.
    Returns (grid, has_trace, [(zero_point, raw, highpassed), ...]).
    Nothing here reads a header value."""
    grid = np.arange(ggate.FP_WAVE[0], ggate.FP_WAVE[1],
                     ggate.FP_STEP_A[r["grism"]])
    if r["status"] != "ok":
        return grid, False, []
    sol = gwave.solution_for(sols, r["grism"], r["night"])
    spec = gstore.load_spec(r["path"])
    if sol is None or spec is None:
        return grid, True, []
    flux = np.where(spec["n_sat"] > 0, np.nan, science_flux(spec))
    out = []
    for cand in gwave.emission_candidates(flux):
        z = gwave.halpha_zero_point(flux, spec["var"], peak=cand)
        if z is None:
            continue
        wave = gwave.wavelength_axis(len(flux), z["x"], sol)
        raw = ggate.normalise(ggate.resample(wave, flux, grid))
        out.append((z, raw, ggate.highpass(grid, raw)))
    return grid, True, out


def _template(members, grism, exclude_night=None):
    use = [f for n, f in members[grism] if n != exclude_night]
    return np.nanmedian(np.array(use), axis=0) if use else None


def run_gate(con) -> None:
    """The pixel gate on every T CrB frame and every control frame, and
    the Halpha identification it implies.

    1. Halpha hypotheses: each frame's strongest sharp emission features
       (on the lrg the zero-order image is often the strongest);
    2. a first T CrB template from the strongest feature of every T CrB
       frame (a median: a minority of wrong picks cannot shape it);
    3. per frame — target or control alike — the hypothesis whose
       spectrum best matches that template (leave-own-night-out for
       T CrB), i.e. each control frame gets its best chance to mimic
       T CrB;
    4. thresholds from T CrB alone, per grism: two-fold cross-validated
       by night parity (out-of-sample completeness), then from all T CrB
       frames for the adopted verdicts;
    5. the pre-registered thresholds applied as well, for the record.
    The chosen Halpha pixel is stored (``x_halpha``) and ``--zero`` uses
    it.
    """
    sols = gwave.load_solutions(con)
    cur = con.execute("""
        SELECT path, sample, target, grism, night, mech_epoch, status,
               pointing_offset_deg FROM g_frames
        WHERE sample IN ('tcrb', 'identity_control') AND grism IS NOT NULL
        ORDER BY night, path""")
    cols = [c[0] for c in cur.description]
    rows = [dict(zip(cols, r)) for r in cur]
    cands = {r["path"]: _candidate_spectra(r, sols) for r in rows}
    v1 = {}
    mcon = gdb.connect_manifest_ro(gstore.default_manifest())
    for path, verdict, reason in mcon.execute(
            "SELECT path, gate_verdict, gate_reason FROM g_extractions "
            "WHERE method = 'flanking'"):
        v1[path] = (verdict, reason)
    # 2. first template: strongest feature of each T CrB frame.
    first = defaultdict(list)
    for r in rows:
        grid, _has, cs = cands[r["path"]]
        if r["sample"] == "tcrb" and cs:
            first[r["grism"]].append((r["night"], cs[0][2]))
    # 3. choose the hypothesis per frame.
    chosen = {}
    for r in rows:
        grid, has, cs = cands[r["path"]]
        if not cs:
            chosen[r["path"]] = None
            continue
        excl = r["night"] if r["sample"] == "tcrb" else None
        tpl = _template(first, r["grism"], excl)
        scores = [ggate.fingerprint_r(grid, c[2], tpl) if tpl is not None
                  else None for c in cs]
        best = max(range(len(cs)), key=lambda i: -2.0 if scores[i] is None
                   else scores[i])
        chosen[r["path"]] = cs[best]
    hp, raw = defaultdict(list), defaultdict(list)
    for r in rows:
        c = chosen[r["path"]]
        if r["sample"] == "tcrb" and c is not None:
            hp[r["grism"]].append((r["night"], c[2]))
            raw[r["grism"]].append((r["night"], c[1]))
    # Per-frame statistics against the final (leave-own-night-out)
    # templates.
    stats = {}
    for r in rows:
        grid, has, _cs = cands[r["path"]]
        c = chosen[r["path"]]
        excl = r["night"] if r["sample"] == "tcrb" else None
        if c is None:
            stats[r["path"]] = (has, None, None, None, None, None)
            continue
        z, fr, f = c
        stats[r["path"]] = (
            has, z["snr"],
            ggate.fingerprint_r(grid, f, _template(hp, r["grism"], excl)),
            ggate.fingerprint_r(grid, fr, _template(raw, r["grism"], excl)),
            ggate.tio_step(grid, fr), z["x"])
    # 4. thresholds from T CrB only.
    thr, cv = {}, {}
    for grism in ("hrg", "lrg"):
        t = [(r, stats[r["path"]]) for r in rows if r["sample"] == "tcrb"
             and r["grism"] == grism and stats[r["path"]][2] is not None]
        thr[grism] = ggate.derive_thresholds([s[2] for _, s in t],
                                             [s[4] for _, s in t])
        k = n = 0
        for parity in (0, 1):
            train = [s for r, s in t if int(r["night"][-1]) % 2 == parity]
            test = [s for r, s in t if int(r["night"][-1]) % 2 != parity]
            rmin, smax = ggate.derive_thresholds([s[2] for s in train],
                                                 [s[4] for s in train])
            for s_ in test:
                n += 1
                k += ggate.fingerprint_verdict(True, s_[1], s_[2], s_[4],
                                               rmin, smax)[0] == "ACCEPT"
        cv[grism] = (k, n)
        gdb.set_g_meta(con, f"gate_thresholds_{grism}", json.dumps(
            {"r_min": thr[grism][0], "tio_step_max": thr[grism][1],
             "cv_accept": k, "cv_n": n}))
    out = []
    for r in rows:
        has, snr, rr, rr_raw, step, x_ha = stats[r["path"]]
        rmin, smax = thr[r["grism"]]
        verdict, reason = ggate.fingerprint_verdict(has, snr, rr, step,
                                                    rmin, smax)
        pre = ggate.fingerprint_verdict(has, snr, rr_raw, step)[0]
        old = v1.get(r["path"], (None, None))
        out.append({
            "path": r["path"], "sample": r["sample"],
            "claimed": r["target"], "grism": r["grism"],
            "mech_epoch": r["mech_epoch"], "night": r["night"],
            "pointing_offset_deg": r["pointing_offset_deg"],
            "old_verdict": old[0], "old_reason": old[1],
            "cc_peak": rr, "cc_raw": rr_raw, "halpha_em_snr": snr,
            "x_halpha": x_ha, "n_candidates": len(cands[r["path"]][2]),
            "tio_step": step, "verdict": verdict, "reason": reason,
            "verdict_prereg": pre,
            "truth": "target" if r["sample"] == "tcrb" else "non_target"})
    con.execute("DELETE FROM g_identity")
    for row in out:
        gdb.upsert(con, "g_identity", row)
    con.commit()
    print(f"identity gate: {len(out)} frames; thresholds "
          f"{ {g: tuple(round(v, 3) for v in t) for g, t in thr.items()} }; "
          f"cross-validated T CrB acceptance {cv}")


# ---------------------------------------------------------------------------
# --variance (G-2)
# ---------------------------------------------------------------------------
def run_variance(con) -> None:
    out = []
    for r in tcrb_rows(con, "AND status = 'ok'"):
        spec = gstore.load_spec(r["path"])
        if spec is None or "skyvar" not in spec or len(spec["skyvar"]) < 5:
            continue
        sv = spec["skyvar"]
        ratio = sv[:, 1] / sv[:, 2]
        row = {"path": r["path"], "sample": r["sample"],
               "grism": r["grism"], "mech_epoch": r["mech_epoch"],
               "detector_key": r["detector_key"], "n_bins": int(len(sv)),
               "level_lo": float(sv[:, 0].min()),
               "level_hi": float(sv[:, 0].max()),
               "var_ratio_median": float(np.median(ratio)),
               "var_ratio_lo": float(np.percentile(ratio, 16)),
               "var_ratio_hi": float(np.percentile(ratio, 84)),
               "bins_json": json.dumps(np.round(sv, 2).tolist())}
        if sv[:, 0].max() - sv[:, 0].min() > 30.0:
            k, ke, rn, rne = ptc_fit(sv[:, 0], sv[:, 1])
            row.update(ptc_gain=k, ptc_gain_err=ke, ptc_rn_adu=rn,
                       ptc_rn_err=rne)
        out.append(row)
    con.execute("DELETE FROM g_variance_check WHERE sample = 'tcrb'")
    for row in out:
        gdb.upsert(con, "g_variance_check", row)
    con.commit()
    print(f"variance checks: {len(out)} frames")


# ---------------------------------------------------------------------------
# --null (G-4, TCRB-A2)
# ---------------------------------------------------------------------------
def run_null(con) -> None:
    out, method_rows = [], []
    for r in tcrb_rows(con, "AND status = 'ok'"):
        spec = gstore.load_spec(r["path"])
        if spec is None or "null_poly_p200" not in spec:
            continue
        inside = spec["inside"]
        cont = float(np.nanmedian(spec["box"][inside]))
        for method in ("poly", "flanking"):
            for off in gext.NULL_OFFSETS:
                tag = f"null_{method}_{'m' if off < 0 else 'p'}{abs(off)}"
                b = spec[tag][inside]
                b = b[np.isfinite(b)]
                if len(b) < 100:
                    continue
                out.append({
                    "path": r["path"], "method": method,
                    "row_offset": int(off), "sample": r["sample"],
                    "grism": r["grism"], "mech_epoch": r["mech_epoch"],
                    "null_median_adu": float(np.median(b)),
                    "null_p16": float(np.percentile(b, 16)),
                    "null_p84": float(np.percentile(b, 84)),
                    "cont_median_adu": cont,
                    "null_frac": float(np.median(b) / cont) if cont else None,
                    "n_cols": int(len(b))})
        # Negative continuum and the two method differences.
        f, v = spec["flux"], spec["var"]
        sel = inside & np.isfinite(f) & np.isfinite(v) & (v > 0)
        neg = {}
        for name, flux in (("poly", f), ("flanking", spec["flux_flanking"])):
            s = sel & np.isfinite(flux)
            neg[name] = float(np.mean(flux[s] < -NEG_SIGMA * np.sqrt(v[s])))
        method_rows.append({
            "path": r["path"], "grism": r["grism"],
            "neg_frac_poly": neg["poly"], "neg_frac_flanking": neg["flanking"],
            "sky_method_diff": gext.median_relative_difference(
                np.where(inside, f, np.nan),
                np.where(inside, spec["flux_flanking"], np.nan)),
            "box_opt_diff": gext.median_relative_difference(
                np.where(inside, f, np.nan),
                np.where(inside, spec["box"], np.nan))})
    con.execute("DELETE FROM g_null_sky WHERE sample = 'tcrb'")
    for row in out:
        gdb.upsert(con, "g_null_sky", row)
    con.execute("DELETE FROM g_method_diff")
    for row in method_rows:
        gdb.upsert(con, "g_method_diff", row)
    con.commit()
    print(f"null-aperture rows: {len(out)}; method rows {len(method_rows)}")


# ---------------------------------------------------------------------------
# --lsf (G-5)
# ---------------------------------------------------------------------------
def _lsf_row(r: dict, spec: dict, x_ha: float, disp: float,
             line_fwhm, offset) -> Optional[dict]:
    """The delivered LSF of one frame at its Halpha column.

    For a point source a slitless spectrograph's line-spread function is
    the stellar image along the dispersion; the image ACROSS the
    dispersion is measured directly as the trace's cross-dispersion FWHM
    (focus is chromatic: the trace is narrowest near the middle of the
    chip and 2-3x broader at its ends, so the FWHM is interpolated to the
    Halpha column).  It is converted with the LOCAL dispersion."""
    fx, fw = spec["fwhm_x"], spec["fwhm_px"]
    good = np.isfinite(fw)
    if good.sum() < 2:
        return None
    xd = float(np.interp(x_ha, fx[good], fw[good]))
    lsf_a = xd * disp
    nominal = NOMINAL_FOCUS_OFFSET.get(r["grism"], 0.0)
    return {"path": r["path"], "sample": r["sample"], "target": r["target"],
            "grism": r["grism"], "mech_epoch": r["mech_epoch"],
            "night": r["night"], "focpos": r["focpos"],
            "focus_offset": offset, "ccd_temp": r["ccd_temp"],
            "fwhm_xd_px": xd, "fwhm_xd_min_px": float(np.nanmin(fw)),
            "line_fwhm_px": line_fwhm, "lsf_fwhm_px": xd,
            "lsf_fwhm_a": lsf_a, "lsf_fwhm_kms": lsf_a / 6562.8 * C_KMS,
            "resolving_power": 6562.8 / lsf_a,
            "off_nominal": int(offset is not None
                               and abs(offset - nominal) > FOCUS_TOL)}


def nightly_lrg_focus(con) -> dict:
    """Median lrg focuser position per night from every reduced frame —
    the reference the hrg offset is measured from (TE.F4)."""
    by = defaultdict(list)
    for night, fp in con.execute(
            "SELECT night, focpos FROM g_frames WHERE grism = 'lrg' "
            "AND focpos IS NOT NULL"):
        by[night].append(fp)
    return {n: float(np.median(v)) for n, v in by.items()}


def run_lsf(con) -> None:
    sols = gwave.load_solutions(con)
    ref = nightly_lrg_focus(con)
    out = []
    # T CrB: Halpha from the zero point.
    zero = {r[0]: r[1:] for r in con.execute(
        "SELECT path, x_halpha, halpha_fwhm_px FROM g_zero_point "
        "WHERE x_halpha IS NOT NULL")}
    for r in tcrb_rows(con, "AND status = 'ok'"):
        if r["path"] not in zero:
            continue
        sol = gwave.solution_for(sols, r["grism"], r["night"])
        spec = gstore.load_spec(r["path"])
        x_ha, ha_fwhm = zero[r["path"]]
        off = (r["focpos"] - ref[r["night"]]
               if r["focpos"] is not None and r["night"] in ref else None)
        row = _lsf_row(r, spec, x_ha, sol.local_disp(x_ha), ha_fwhm, off)
        if row:
            out.append(row)
    # Calibrators: Halpha centre from the line measurements.
    cur = con.execute("""
        SELECT f.path, f.sample, f.star AS target, f.grism, f.night,
               f.mech_epoch, f.focpos, f.ccd_temp, m.x
        FROM g_frames f JOIN g_line_meas m USING (path)
        WHERE f.sample = 'calibrator' AND m.line = 'Halpha'""")
    cols = [c[0] for c in cur.description]
    for vals in cur.fetchall():
        r = dict(zip(cols, vals))
        sol = gwave.solution_for(sols, r["grism"], r["night"])
        spec = gstore.load_spec(r["path"])
        if sol is None or spec is None:
            continue
        off = (r["focpos"] - ref[r["night"]]
               if r["focpos"] is not None and r["night"] in ref else None)
        row = _lsf_row(r, spec, r["x"], sol.local_disp(r["x"]), None, off)
        if row:
            out.append(row)
    con.execute("DELETE FROM g_lsf")
    for row in out:
        gdb.upsert(con, "g_lsf", row)
    con.commit()
    print(f"LSF rows: {len(out)}")


# ---------------------------------------------------------------------------
# --sat (TCRB-A6)
# ---------------------------------------------------------------------------
def run_sat(con) -> None:
    sols = gwave.load_solutions(con)
    zero = {r[0]: r[1] for r in con.execute(
        "SELECT path, x_halpha FROM g_zero_point WHERE x_halpha IS NOT NULL")}
    out = []
    for r in tcrb_rows(con, "AND status = 'ok'"):
        spec = gstore.load_spec(r["path"])
        inside = spec["inside"]
        sat_cols = inside & (spec["n_sat"] > 0)
        n_ap = int(inside.sum()) * (2 * gext.APERTURE_HALFWIN + 1)
        hot_frac = float(spec["n_hot"][inside].sum() / n_ap) if n_ap else None
        near = None
        sol = gwave.solution_for(sols, r["grism"], r["night"])
        if r["path"] in zero and sol is not None:
            wave = gwave.wavelength_axis(len(inside), zero[r["path"]], sol)
            near = int((sat_cols & (np.abs(wave - 6562.8)
                                    <= SAT_DISCARD_HALFWIDTH_A)).sum())
        verdict = ("discard" if near else
                   "flag" if sat_cols.any() else "clean")
        out.append({"path": r["path"], "grism": r["grism"],
                    "night": r["night"], "peak_adu": r["peak_adu"],
                    "sat_cap_adu": r["sat_cap_adu"],
                    "n_sat_cols": int(sat_cols.sum()),
                    "n_sat_cols_halpha": near, "hot_mask": r["hot_mask"],
                    "masked_frac": hot_frac, "verdict": verdict})
    con.execute("DELETE FROM g_saturation")
    for row in out:
        gdb.upsert(con, "g_saturation", row)
    con.commit()
    print(f"saturation triage: {len(out)} frames")


# ---------------------------------------------------------------------------
# --ew (old library vs new, one estimator)
# ---------------------------------------------------------------------------
def v1_rows(mcon) -> dict:
    """The August library's flanking-method rows with an Halpha anchor."""
    cur = mcon.execute("""
        SELECT path, filter, x_halpha, disp_a_per_px, gate_verdict,
               spectrum_fits, n_sat_cols
        FROM g_extractions
        WHERE method = 'flanking' AND target = 'T CrB'
          AND spectrum_fits IS NOT NULL AND x_halpha IS NOT NULL""")
    cols = [c[0] for c in cur.description]
    return {r[0]: dict(zip(cols, r)) for r in cur}


def run_ew(con, mcon) -> None:
    from astropy.io import fits
    sols = gwave.load_solutions(con)
    zero = {r[0]: r[1] for r in con.execute(
        "SELECT path, x_halpha FROM g_zero_point WHERE x_halpha IS NOT NULL")}
    gate = {r[0]: r[1] for r in con.execute(
        "SELECT path, verdict FROM g_identity")}
    v1 = v1_rows(mcon)
    # v1's fallback for frames without an O2 pair: the per-grism median of
    # its own per-frame dispersions (as its parquet stage did).
    v1_med = {}
    for g in ("hrg", "lrg"):
        d = [r["disp_a_per_px"] for r in v1.values()
             if r["filter"] == g and r["disp_a_per_px"]]
        if d:
            v1_med[g] = float(np.median(d))
    out = []
    for r in tcrb_rows(con, "AND status = 'ok'"):
        base = {"path": r["path"], "grism": r["grism"], "night": r["night"],
                "jd": r["jd"], "mech_epoch": r["mech_epoch"]}
        sol = gwave.solution_for(sols, r["grism"], r["night"])
        if sol is not None and r["path"] in zero:
            spec = gstore.load_spec(r["path"])
            x_ha = zero[r["path"]]
            wave = gwave.wavelength_axis(len(spec["flux"]), x_ha, sol)
            flux = np.where(spec["n_sat"] > 0, np.nan, science_flux(spec))
            m = gew.equivalent_width(wave, flux, spec["var"])
            mb = gew.equivalent_width(wave, science_flux(spec, "box"),
                                      spec["box_var"])
            row = dict(base, library="v2", gate=gate.get(r["path"]),
                       x_halpha=x_ha, disp_a_per_px=sol.local_disp(x_ha),
                       disp_source="fixed (g_dispersion)",
                       sky_method="poly", n_sat_cols=r["n_sat_cols"],
                       status="ok" if m else "no_ew")
            if m:
                row.update(ew_a=m["ew_a"], ew_err_a=m["ew_err_a"],
                           ew_px=m["ew_px"], cont_adu=m["cont"],
                           cont_snr=m["cont_snr"], line_peak_adu=m["peak"],
                           n_masked_line=m["n_masked"],
                           ew_box_a=mb["ew_a"] if mb else None)
            out.append(row)
        if r["path"] in v1:
            o = v1[r["path"]]
            fpath = gconfig.PRODUCTS / o["spectrum_fits"]
            disp = o["disp_a_per_px"] or v1_med.get(o["filter"])
            if not fpath.exists() or not disp:
                continue
            with fits.open(fpath) as h:
                t = h["FLANKING"].data
                f1 = np.array(t["FLUX_ADU"], dtype=float)
                var1 = np.array(t["VAR_ADU2"], dtype=float)
                nsat = np.array(t["N_SAT"])
            wave1 = 6562.8 + disp * (np.arange(len(f1)) - o["x_halpha"])
            m = gew.equivalent_width(wave1, f1, var1)
            row = dict(base, library="v1", gate=o["gate_verdict"],
                       x_halpha=o["x_halpha"], disp_a_per_px=abs(disp),
                       disp_source=("per-frame O2 pair"
                                    if o["disp_a_per_px"] else
                                    "v1 grism median"),
                       sky_method="flanking", n_sat_cols=o["n_sat_cols"],
                       status="ok" if m else "no_ew")
            if m:
                lo = int(max(0, o["x_halpha"] - 40))
                row.update(ew_a=m["ew_a"], ew_err_a=m["ew_err_a"],
                           ew_px=m["ew_px"], cont_adu=m["cont"],
                           cont_snr=m["cont_snr"], line_peak_adu=m["peak"],
                           n_masked_line=int((nsat[lo:lo + 80] > 0).sum()))
            out.append(row)
    con.execute("DELETE FROM g_ew")
    for row in out:
        gdb.upsert(con, "g_ew", row)
    con.commit()
    print(f"EW rows: {len(out)}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--db", default=str(gconfig.GRISM_DB))
    ap.add_argument("--manifest", default=str(gstore.default_manifest()))
    stages = ("zero", "gate", "variance", "null", "lsf", "sat", "ew",
              "summary", "report")
    for stage in stages + ("all",):
        ap.add_argument(f"--{stage}", action="store_true")
    args = ap.parse_args(argv)
    con = gdb.connect_grism(args.db)
    if args.gate or args.all:
        run_gate(con)
    if args.zero or args.all:
        run_zero(con)
    if args.variance or args.all:
        run_variance(con)
    if args.null or args.all:
        run_null(con)
    if args.lsf or args.all:
        run_lsf(con)
    if args.sat or args.all:
        run_sat(con)
    if args.ew or args.all:
        run_ew(con, gdb.connect_manifest_ro(args.manifest))
    if args.summary or args.all:
        from macro_grism import summary_g
        print(summary_g.text_summary(con))
    if args.report or args.all:
        from macro_grism.report_g import render_report
        print(f"report: {render_report()}")
    gdb.set_g_meta(con, "validation_code", ggate.G_CODE_VERSION)
    con.commit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
