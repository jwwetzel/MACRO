"""Every acceptance number of the grism library, computed from the grism
database — one function, used by the runner's ``--summary`` and by the
report page, so the two can never disagree.

``numbers(con)`` returns a nested dict; ``text_summary(con)`` formats it.
Nothing here is typed: each value is a query or a statistic of one.
"""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from rlmt_diagnostics.dispersion import wilson_interval

#: The two disputed high-resolution dispersions and the lrg values on
#: file before 2026-10-03 — the hypotheses the D1 test compares.
HRG_DISPUTED = (0.47, 1.59)
LEVER_O2B_A = 6866.2 - 6562.8     # Halpha to the O2-B edge (A)


def mad_sigma(v) -> float:
    v = np.asarray(v, dtype=float)
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def _q(con, sql, args=()):
    return con.execute(sql, args).fetchall()


def d1_numbers(con) -> dict:
    """D1: the identification under the adopted seed vs the rivals, and
    the depth where 1.59 A/px puts Hbeta."""
    out = {}
    for grism in ("hrg", "lrg"):
        rows = _q(con, """SELECT n_match, n_stellar, rms_px, rivals_json,
                                 hbeta_depth_159, hbeta_depth_159_err,
                                 refined, spt
                          FROM g_line_id WHERE grism = ?""", (grism,))
        if not rows:
            continue
        best_rival, wins, ties = [], 0, 0
        for nm, ns, rms, rj, *_ in rows:
            rv = json.loads(rj)
            br = max(v["n_match"] for v in rv.values())
            best_rival.append(br)
            wins += nm > br
            ties += nm == br
        nm = np.array([r[0] for r in rows])
        d = {"n_frames": len(rows),
             "adopted_median_matches": float(np.median(nm)),
             "rival_median_matches": float(np.median(best_rival)),
             "adopted_wins": int(wins), "ties": int(ties),
             "refined_frames": int(sum(r[6] or 0 for r in rows))}
        # The B stars carry the narrow He I / Si II / Ne I lines whose
        # spacings no wrong scale reproduces: the decisive subset.
        b = [(r[0], max(v["n_match"] for v in json.loads(r[3]).values()))
             for r in rows if r[7] == "B"]
        if b:
            d["b_frames"] = len(b)
            d["b_adopted_median"] = float(np.median([x[0] for x in b]))
            d["b_rival_median"] = float(np.median([x[1] for x in b]))
            d["b_adopted_wins"] = int(sum(x[0] > x[1] for x in b))
            d["b_rival_wins"] = int(sum(x[0] < x[1] for x in b))
        if grism == "hrg":
            hb = np.array([r[4] for r in rows if r[4] is not None])
            err = np.array([r[5] for r in rows if r[4] is not None])
            if len(hb):
                d["hbeta159_depth_median"] = float(np.median(hb))
                d["hbeta159_noise_median"] = float(np.median(err))
                d["hbeta159_n_detected_5sigma"] = int((hb > 5 * err).sum())
                d["hbeta159_n"] = int(len(hb))
        out[grism] = d
    return out


def dispersion_rows(con) -> list[dict]:
    cur = con.execute("""SELECT grism, mech_epoch, disp_a_per_px, disp_err,
            disp_minus1000, disp_plus1000, x_ref, degree, n_frames,
            n_stars, n_lines, n_meas, n_clipped, rms_px, rms_a, chi2, dof,
            rchi2, o2b_wave_eff, o2b_wave_err, stars, lines, status, note
            FROM g_dispersion ORDER BY grism, mech_epoch""")
    cols = [c[0] for c in cur.description]
    return [dict(zip(cols, r)) for r in cur]


def lsf_table(con) -> list[dict]:
    """Median delivered LSF per (grism, epoch, sample)."""
    rows = _q(con, """SELECT grism, mech_epoch, sample, lsf_fwhm_px,
                             lsf_fwhm_a, lsf_fwhm_kms, resolving_power,
                             off_nominal FROM g_lsf""")
    groups = defaultdict(list)
    for r in rows:
        groups[(r[0], r[1], r[2])].append(r[3:])
    out = []
    for (g, e, s), v in sorted(groups.items()):
        a = np.array(v, dtype=float)
        out.append({"grism": g, "mech_epoch": e, "sample": s, "n": len(a),
                    "fwhm_px": float(np.median(a[:, 0])),
                    "fwhm_a": float(np.median(a[:, 1])),
                    "fwhm_a_p16": float(np.percentile(a[:, 1], 16)),
                    "fwhm_a_p84": float(np.percentile(a[:, 1], 84)),
                    "fwhm_kms": float(np.median(a[:, 2])),
                    "R": float(np.median(a[:, 3])),
                    "off_nominal": int(np.nansum(a[:, 4]))})
    return out


def lsf_regressions(con) -> list[dict]:
    """LSF (px) against focus offset and CCD-TEMP, T CrB frames, per
    grism: least-squares slope, Pearson r, n."""
    out = []
    for grism in ("hrg", "lrg"):
        rows = np.array(_q(con, """SELECT lsf_fwhm_px, focus_offset, ccd_temp
                FROM g_lsf WHERE sample = 'tcrb' AND grism = ?""", (grism,)),
                        dtype=float)
        if not len(rows):
            continue
        for name, col in (("focus_offset", 1), ("ccd_temp", 2)):
            ok = np.isfinite(rows[:, 0]) & np.isfinite(rows[:, col])
            if ok.sum() < 10 or np.ptp(rows[ok, col]) == 0:
                continue
            c = np.polyfit(rows[ok, col], rows[ok, 0], 1)
            r = float(np.corrcoef(rows[ok, col], rows[ok, 0])[0, 1])
            out.append({"grism": grism, "regressor": name,
                        "slope_px_per_unit": float(c[0]), "r": r,
                        "n": int(ok.sum()),
                        "range": (float(rows[ok, col].min()),
                                  float(rows[ok, col].max()))})
    return out


def numbers(con) -> dict:
    n = {"d1": d1_numbers(con), "dispersion": dispersion_rows(con),
         "lsf": lsf_table(con), "lsf_regressions": lsf_regressions(con)}
    # ---- G-1 on T CrB: O2-B edge constancy --------------------------------
    g1 = {}
    for grism in ("hrg", "lrg"):
        tot = _q(con, """SELECT count(*), sum(x_halpha IS NOT NULL),
                                sum(o2b_wave IS NOT NULL)
                         FROM g_zero_point WHERE sample = 'tcrb'
                           AND grism = ?""", (grism,))[0]
        rows = np.array(_q(con, """SELECT o2b_dwave, sep_px, sep_frac_dev,
                                          o2b_wave_err
                FROM g_zero_point WHERE sample = 'tcrb' AND grism = ?
                  AND o2b_wave IS NOT NULL""", (grism,)), dtype=float)
        d = {"n_frames": tot[0], "n_halpha": tot[1], "n_o2b": tot[2]}
        if len(rows):
            d.update(dwave_median=float(np.median(rows[:, 0])),
                     dwave_scatter=mad_sigma(rows[:, 0]),
                     dwave_scatter_frac=mad_sigma(rows[:, 0]) / LEVER_O2B_A,
                     sep_px_median=float(np.median(rows[:, 1])),
                     sep_px_scatter_frac=mad_sigma(rows[:, 1])
                     / float(np.median(rows[:, 1])),
                     n_dev_gt_2pct=int((np.abs(rows[:, 2]) > 0.02).sum()),
                     edge_err_median=float(np.median(rows[:, 3])))
        g1[grism] = d
    n["g1_tcrb"] = g1
    # ---- G-2 ----------------------------------------------------------------
    g2 = {}
    for grism in ("hrg", "lrg"):
        rows = _q(con, """SELECT var_ratio_median, bins_json
                          FROM g_variance_check WHERE sample = 'tcrb'
                            AND grism = ?""", (grism,))
        if not rows:
            continue
        ratio = np.array([r[0] for r in rows])
        bins = np.vstack([np.array(json.loads(r[1])) for r in rows])
        k, ke, rn, rne = ptc_fit(bins[:, 0], bins[:, 1])
        g2[grism] = {"n_frames": len(rows),
                     "ratio_median": float(np.median(ratio)),
                     "ratio_p16": float(np.percentile(ratio, 16)),
                     "ratio_p84": float(np.percentile(ratio, 84)),
                     "n_within_20pct": int((np.abs(ratio - 1) <= 0.2).sum()),
                     "ptc_gain": k, "ptc_gain_err": ke, "ptc_floor_adu": rn,
                     "ptc_floor_err": rne, "n_bins": int(len(bins)),
                     "level_range": (float(bins[:, 0].min()),
                                     float(bins[:, 0].max()))}
    n["g2"] = g2
    # ---- G-3 / A0 -----------------------------------------------------------
    g3 = {}
    for truth in ("target", "non_target"):
        rows = _q(con, """SELECT verdict, reason, pointing_offset_deg,
                                 old_verdict, old_reason, grism
                          FROM g_identity WHERE truth = ?""", (truth,))
        acc = sum(1 for r in rows if r[0] == "ACCEPT")
        lo, hi = wilson_interval(acc, len(rows))
        reasons = defaultdict(int)
        for r in rows:
            reasons[r[1]] += 1
        d = {"n": len(rows), "accept": acc, "wilson95": (lo, hi),
             "reasons": dict(reasons),
             "by_grism": {g: (sum(1 for r in rows if r[5] == g
                                  and r[0] == "ACCEPT"),
                              sum(1 for r in rows if r[5] == g))
                          for g in ("hrg", "lrg")}}
        if truth == "target":
            hdr_bad = [r for r in rows if (r[2] or 0) > 1.0]
            d["header_off_n"] = len(hdr_bad)
            d["header_off_accept"] = sum(1 for r in hdr_bad
                                         if r[0] == "ACCEPT")
            d["v1_header_rejects"] = sum(1 for r in rows
                                         if r[4] == "header_off_target")
            d["v1_header_rejects_now_accepted"] = sum(
                1 for r in rows if r[4] == "header_off_target"
                and r[0] == "ACCEPT")
            # Counts by S2c verdict are attached by the report (manifest).
        g3[truth] = d
    n["g3"] = g3
    # ---- G-4 / A2 -----------------------------------------------------------
    g4 = {}
    for grism in ("hrg", "lrg"):
        d = {}
        for method in ("poly", "flanking"):
            v = np.array([r[0] for r in _q(con, """SELECT null_frac FROM
                g_null_sky WHERE sample = 'tcrb' AND grism = ? AND method = ?
                AND null_frac IS NOT NULL""", (grism, method))])
            if len(v):
                d[method] = {"n": len(v), "median": float(np.median(v)),
                             "scatter": mad_sigma(v),
                             "n_gt_3pct": int((np.abs(v) > 0.03).sum())}
        rows = np.array(_q(con, """SELECT neg_frac_poly, neg_frac_flanking,
                sky_method_diff, box_opt_diff FROM g_method_diff
                WHERE grism = ?""", (grism,)), dtype=float)
        if len(rows):
            d["neg_frac_poly_median"] = float(np.nanmedian(rows[:, 0]))
            d["neg_frac_poly_max"] = float(np.nanmax(rows[:, 0]))
            d["neg_frac_flanking_median"] = float(np.nanmedian(rows[:, 1]))
            d["neg_frac_flanking_max"] = float(np.nanmax(rows[:, 1]))
            d["sky_method_diff_median"] = float(np.nanmedian(rows[:, 2]))
            d["sky_method_diff_p84"] = float(np.nanpercentile(rows[:, 2], 84))
            d["box_opt_diff_median"] = float(np.nanmedian(rows[:, 3]))
            d["box_opt_diff_p84"] = float(np.nanpercentile(rows[:, 3], 84))
            d["n_frames"] = int(len(rows))
        g4[grism] = d
    n["g4"] = g4
    # ---- A6 -----------------------------------------------------------------
    a6 = {}
    for grism in ("hrg", "lrg"):
        rows = _q(con, """SELECT verdict, masked_frac, peak_adu, sat_cap_adu
                          FROM g_saturation WHERE grism = ?""", (grism,))
        if not rows:
            continue
        cnt = defaultdict(int)
        for r in rows:
            cnt[r[0]] += 1
        mf = np.array([r[1] for r in rows if r[1] is not None])
        a6[grism] = {"n": len(rows), "verdicts": dict(cnt),
                     "masked_frac_median": float(np.median(mf))
                     if len(mf) else None,
                     "peak_max": float(max(r[2] or 0 for r in rows)),
                     "cap": float(rows[0][3])}
    n["a6"] = a6
    # ---- EW, old library vs new ---------------------------------------------
    ew = {}
    for grism in ("hrg", "lrg"):
        rows = np.array(_q(con, """SELECT a.ew_a, a.ew_err_a, b.ew_a,
                b.ew_err_a, a.disp_a_per_px, b.disp_a_per_px
                FROM g_ew a JOIN g_ew b USING (path)
                WHERE a.library = 'v1' AND b.library = 'v2' AND a.grism = ?
                  AND a.ew_a IS NOT NULL AND b.ew_a IS NOT NULL""",
                           (grism,)), dtype=float)
        allv2 = np.array(_q(con, """SELECT ew_a, ew_err_a FROM g_ew
                WHERE library = 'v2' AND grism = ? AND ew_a IS NOT NULL
                  AND gate = 'ACCEPT'""", (grism,)), dtype=float)
        d = {}
        if len(allv2):
            d["v2_n"] = int(len(allv2))
            d["v2_median"] = float(np.median(allv2[:, 0]))
            d["v2_p16"] = float(np.percentile(allv2[:, 0], 16))
            d["v2_p84"] = float(np.percentile(allv2[:, 0], 84))
            d["v2_err_median"] = float(np.median(allv2[:, 1]))
        if len(rows):
            ratio = rows[:, 2] / rows[:, 0]
            d.update(paired_n=int(len(rows)),
                     v1_median=float(np.median(rows[:, 0])),
                     v2_median_paired=float(np.median(rows[:, 2])),
                     ratio_median=float(np.median(ratio)),
                     ratio_p16=float(np.percentile(ratio, 16)),
                     ratio_p84=float(np.percentile(ratio, 84)),
                     v1_err_median=float(np.median(rows[:, 1])),
                     v2_err_median_paired=float(np.median(rows[:, 3])),
                     v1_scatter_frac=mad_sigma(rows[:, 0])
                     / abs(float(np.median(rows[:, 0]))),
                     v2_scatter_frac=mad_sigma(rows[:, 2])
                     / abs(float(np.median(rows[:, 2]))),
                     v1_disp_range=(float(rows[:, 4].min()),
                                    float(rows[:, 4].max())),
                     v2_disp_range=(float(rows[:, 5].min()),
                                    float(rows[:, 5].max())))
        ew[grism] = d
    n["ew"] = ew
    return n


def ptc_fit(level: np.ndarray, var: np.ndarray):
    """Photon-transfer line var = RN^2 + level / K by least squares over
    sky-level bins; (K, K_err, RN, RN_err) with errors from the fit
    covariance scaled by the residual scatter."""
    a = np.vstack([np.ones(len(level)), level]).T
    coef, *_ = np.linalg.lstsq(a, var, rcond=None)
    dof = max(1, len(level) - 2)
    s2 = float(np.sum((var - a @ coef) ** 2) / dof)
    cov = s2 * np.linalg.inv(a.T @ a)
    k = 1.0 / coef[1]
    k_err = np.sqrt(cov[1, 1]) / coef[1] ** 2
    rn = np.sqrt(max(coef[0], 0.0))
    rn_err = np.sqrt(cov[0, 0]) / (2.0 * max(rn, 1e-9))
    return float(k), float(k_err), float(rn), float(rn_err)


def text_summary(con) -> str:
    """The numbers as indented text (for the runner's --summary)."""
    def fmt(v, ind=0):
        pad = "  " * ind
        if isinstance(v, dict):
            return "\n".join(f"{pad}{k}:" + ("\n" + fmt(x, ind + 1)
                                            if isinstance(x, (dict, list))
                                            else f" {_f(x)}")
                             for k, x in v.items())
        if isinstance(v, list):
            return "\n".join(f"{pad}- " + ", ".join(
                f"{k}={_f(x)}" for k, x in row.items()) if isinstance(row,
                dict) else f"{pad}- {_f(row)}" for row in v)
        return pad + _f(v)
    return fmt(numbers(con))


def _f(x) -> str:
    if isinstance(x, float):
        return f"{x:.4g}"
    if isinstance(x, tuple):
        return "(" + ", ".join(_f(v) for v in x) + ")"
    return str(x)
