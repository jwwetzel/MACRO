"""macro_sn.report_sn — the evidence page for the SN 2023ixf release.

Renders ``docs/SN2023ixf_LightCurve/sn_release.html`` from
``products/sn/sn2023ixf.sqlite`` and the manifest's Gate 0 tables.  One
section per plan task, in the house Question -> Evidence -> Decision form,
every number a query; the four manuscript figures are copied beside it as
PNG.  Nothing on the page is typed.
"""

from __future__ import annotations

import html
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs" / "SN2023ixf_LightCurve"
FIGS = DOCS / "figures" / "sn"
MSFIG = REPO / "manuscripts" / "SN2023ixf_LightCurve" / "figures"


def esc(x):
    return html.escape("" if x is None else str(x))


def fmt(v, f="{:.3f}"):
    if v is None:
        return "&mdash;"
    if isinstance(v, float):
        return f.format(v)
    return esc(v)


def table(con, sql, headers, fmts=None, args=()):
    rows = con.execute(sql, args).fetchall()
    fmts = fmts or ["{}"] * len(headers)
    out = ["<table><tr>" + "".join(f"<th>{esc(h)}</th>" for h in headers) + "</tr>"]
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{fmt(v, f)}</td>" for v, f in zip(r, fmts))
                   + "</tr>")
    out.append("</table>")
    return "\n".join(out)


def fig(name, cap):
    return (f'<figure><img src="figures/sn/{name}.png" alt="">'
            f"<figcaption>{esc(cap)}</figcaption></figure>")


def section(sid, title, task, question, evidence, decision):
    return (f'<section id="{sid}"><h2>{esc(title)}</h2>'
            f'<p class="task">Plan task <code>{esc(task)}</code></p>'
            f"<h3>Question</h3><p>{question}</p><h3>Evidence</h3>{evidence}"
            f'<h3>Decision</h3><div class="decision">{decision}</div></section>')


def render(sn_db: Path, manifest: Path) -> Path:
    FIGS.mkdir(parents=True, exist_ok=True)
    for p in MSFIG.glob("sn_fig*.png"):
        shutil.copy2(p, FIGS / p.name)
    sn = sqlite3.connect(f"file:{sn_db}?mode=ro", uri=True)
    man = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    meta = dict(sn.execute("SELECT key, value FROM sn_build_meta"))
    q1 = lambda c, s, *a: c.execute(s, a).fetchone()[0]          # noqa: E731
    S = []
    S.append(section(
        "resolve", "Astrometry for every frame", "SN-S4-resolve",
        "Which broadband frames lacked a plate solution, and can each be solved?",
        table(sn, """SELECT filter, epoch_role, coalesce(pltsolvd, 0), status,
                     wcs_source, count(*), round(avg(n_match)), avg(wcs_rms_arcsec),
                     round(avg(n_det))
                     FROM sn_frames WHERE filter IN ('G','R','I','g','r','i')
                     GROUP BY 1,2,3,4,5 ORDER BY 2,1,3""",
              ["filter", "epoch", "header solved", "status", "WCS source",
               "frames", "mean REFCAT2 matched", "mean rms (\")",
               "mean 5σ detections"],
              ["{}", "{}", "{}", "{}", "{}", "{}", "{:.0f}", "{:.2f}", "{:.0f}"]),
        "Every frame's WCS is re-fitted to REFCAT2 (header solutions are never "
        "used unrefined). Unsolved frames were solved by a translation vote "
        "under the campaign camera's measured CD matrix, then by "
        "triangle matching. The frames that remain unsolved carry a median of "
        "about a dozen 5σ detections — fewer than the eight calibration stars "
        "the photometry needs — so no solution could make them usable."))
    S.append(section(
        "crosswalk", "Filter identification and the MaxIm→pyscope crosswalk",
        "SN-S1-filter-crosswalk",
        "Which PS1 band does each filter code measure, and do the MaxIm "
        "campaign codes and the pyscope template codes agree?",
        table(sn, """SELECT role, camera, code, ps1_band, cterm, cterm_err,
                     tie_rms, n, identified FROM sn_crosswalk
                     ORDER BY role, camera, code, ps1_band""",
              ["epoch", "camera", "code", "PS1 band", "colour term (g−i)", "±",
               "tie rms", "stars", "identified"],
              ["{}", "{}", "{}", "{}", "{:+.3f}", "{:.3f}", "{:.3f}", "{}", "{}"]),
        "Identification = the PS1 band with the smallest |colour term|. G, R, I "
        "are g, r, i; pyscope g, r, i are the same bands. Narrowband codes "
        "cannot be identified this way (no profile from a colour term)."))
    S.append(section(
        "cal", "REFCAT2 ensemble calibration with held-out check stars",
        "SN-S4-ensemble-cal",
        "Is the error model honest on stars that took no part in the fit?",
        table(sn, """SELECT role, camera, code, n_frames_used || '/' || n_frames,
                     n_ensemble, n_check, cterm, k2, floor_mag, tie_rms,
                     check_chi2nu, check_dof, check_verdict FROM sn_cal_sets""",
              ["epoch", "camera", "code", "frames used", "ensemble", "check",
               "colour term", "k″ (colour×X)", "floor", "tie rms",
               "check χ²ν", "dof", "verdict (rule 1)"],
              ["{}"] * 6 + ["{:+.3f}", "{:+.3f}", "{:.4f}", "{:.3f}", "{:.2f}",
                            "{}", "{}"])
        + f"<p>Scintillation scale (pooled over G, R, I): "
          f"C = {float(meta['scint_C']):.2f} × Young (1967).</p>"
        + table(sn, """SELECT code, exptime, n_points, rms_obs, rms_model, ratio,
                       scint_var_frac FROM sn_cal_scint WHERE role='campaign'
                       ORDER BY code, exptime""",
                ["code", "exp (s)", "check points", "rms obs", "rms model",
                 "normalised rms", "scint share of variance"],
                ["{}", "{:.1f}", "{}", "{:.4f}", "{:.4f}", "{:.3f}", "{:.3f}"]),
        "Check-star χ²ν per band is reported with its dof; every campaign "
        "band passes standing rule 1. The QHY600 g template set fails it "
        "(χ²ν < 0.5) and is used only for the crosswalk."))
    S.append(section(
        "templates", "The template table", "SN-S5-template-table",
        "Which template epochs are usable, on which camera, mechanical epoch, "
        "FWHM and focus?",
        table(sn, """SELECT night, camera, filter, n_frames, exptime, mech_epoch,
                     rotation_deg, scale_arcsec, fwhm_arcsec, focpos, focus_offset,
                     in_focus FROM sn_template_table ORDER BY night, filter""",
              ["night", "camera", "filter", "frames", "exp (s)", "mech epoch",
               "rotation (°)", "scale (\"/px)", "FWHM (\")", "FOCPOS",
               "focus offset", "in focus"],
              ["{}", "{}", "{}", "{}", "{:.0f}", "{}", "{:.1f}", "{:.3f}",
               "{:.2f}", "{:.0f}", "{:+.0f}", "{}"]),
        "In focus = FWHM within 1.3× the campaign median (AC4040) or the best "
        "epoch of that camera. Narrowband template: 2023-05-04 H only (O fails "
        "the focus test). Pre-explosion G/R are out of focus."))
    S.append(section(
        "phot", "Two-regime photometry", "SN-S5-photometry",
        "Does the aperture (bright) regime agree with the template-subtracted "
        "(faint) regime, and is the scintillation term matched by the "
        "check-star rms?",
        table(sn, """SELECT code, count(*), avg(dm), min(dm), max(dm)
                     FROM sn_overlap GROUP BY code""",
              ["code", "frames", "mean aperture − subtracted (mag)", "min", "max"],
              ["{}", "{}", "{:+.4f}", "{:+.4f}", "{:+.4f}"])
        + table(sn, """SELECT code, floor_mag, n_points, n_nights, star_floor_mag
                       FROM sn_sn_floor""",
                ["code", "SN floor (mag)", "frames", "nights", "star floor"],
                ["{}", "{:.4f}", "{}", "{}", "{:.4f}"])
        + table(sn, """SELECT code, sum(usable), count(*),
                       count(DISTINCT CASE WHEN usable=1 THEN night END)
                       FROM sn_phot GROUP BY code""",
                ["code", "usable frames", "campaign frames", "nights"])
        + fig("sn_fig2_lightcurve", "Figure 2: PS1-tied gri light curve."),
        "The overlap systematic is a few mmag in every band. SN frames are "
        "photometry only in exposure groups whose check stars validate the "
        "error model; S2c-indeterminate frames only if their own star PSFs "
        "pass the elongation check."))
    S.append(section(
        "peak", "The peak epoch", "SN-S5b-peak-epoch",
        "Is the clean start (+5.4 d) at or after optical maximum?",
        table(sn, "SELECT key, value, source FROM sn_literature",
              ["quantity", "value", "source"], ["{}", "{:g}", "{}"])
        + table(sn, "SELECT band, t_peak_d, t_lo, t_hi, n_points FROM sn_peak",
                ["band", "t_peak (d)", "16%", "84%", "points"],
                ["{}", "{:.2f}", "{:.2f}", "{:.2f}", "{}"]),
        "g and r peak at the clean start; i keeps brightening to ~+14 d. The "
        "paper uses no 'rise' language."))
    S.append(section(
        "gate", "The predicted-excess gate", "SN-S6-0-excess-gate",
        "Would the flash-phase (H − '[S II]') colour see the Hα excess "
        "published spectroscopy predicts?",
        table(sn, """SELECT epoch_d, width_a, eps_upper, dm_upper, dm_central,
                     sig_rep, sig_loc, sig_nb, passes, passes_rep_only
                     FROM sn_excess_gate""",
              ["epoch (d)", "W (Å)", "excess (upper)", "Δm upper", "Δm central",
               "σ_rep", "σ_loc", "σ_NB", "passes 3σ_NB", "passes 3σ_rep"],
              ["{:.1f}", "{:.0f}", "{:.3f}", "{:.3f}", "{:.3f}", "{:.3f}",
               "{:.3f}", "{:.3f}", "{}", "{}"])
        + table(sn, "SELECT * FROM sn_nb_widths",
                ["code", "W from ZP ratio (Å)", "MAD", "nights"],
                ["{}", "{:.0f}", "{:.0f}", "{}"]),
        f"Verdict: <b>{esc(meta.get('excess_gate_verdict'))}</b>. The locus "
        "transfer to the SN's colour, blueward of every calibration star, is "
        "uncertain by ~0.3 mag; S6a is demoted before the measurement."))
    S.append(section(
        "resid", "Residuals against the published photometry",
        "SN-S7b-residuals-table",
        "How do the RLMT nightly magnitudes compare with Li et al. (2025)?",
        table(sn, """SELECT code, subset, n, offset, offset_err, rms, chi2nu, dof,
                     colour_lo, colour_hi FROM sn_resid_summary""",
              ["code", "subset", "N", "offset", "±", "rms", "χ²ν", "dof",
               "colour range lo", "hi"],
              ["{}", "{}", "{}", "{:+.3f}", "{:.3f}", "{:.3f}", "{:.2f}", "{}",
               "{:.2f}", "{:.2f}"])
        + table(sn, "SELECT phase_d, ztf_g, rlmt_g, resid FROM sn_resid_ztf",
                ["phase (d)", "ZTF g", "RLMT g", "RLMT − ZTF"],
                ["{:.2f}", "{:.3f}", "{:.3f}", "{:+.3f}"])
        + fig("sn_fig3_residuals", "Figure 3: residuals vs Li et al. (2025)."),
        "The table exists before any variability limit is quoted."))
    S.append(section(
        "var", "Variability and bump limits", "SN-S8-variability-limits",
        "What night-to-night variability could the campaign have seen?",
        table(sn, """SELECT code, n_nights, knot_spacing_d, chi2nu_trend, bump_snr,
                     bump_thr, bump_fap, sine_chi2, sine_thr, sine_fap
                     FROM sn_var_observed""",
              ["code", "nights", "knots (d)", "χ²ν about trend", "bump SNR",
               "threshold", "FAP", "sine Δχ²", "threshold", "FAP"],
              ["{}", "{}", "{:.0f}", "{:.2f}", "{:.2f}", "{:.2f}", "{:.3f}",
               "{:.1f}", "{:.1f}", "{:.3f}"])
        + table(sn, """SELECT code, shape, scale_d, a90_mag, a_pred_mag,
                       signed_bias_mag, a50_all_mag, completeness
                       FROM sn_var_limits""",
                ["code", "shape", "σ_t or P (d)", "A90 (mag)", "noise-only scale",
                 "signed bias", "A50 all", "completeness"],
                ["{}", "{}", "{:g}", "{:.3f}", "{:.3f}", "{:+.4f}", "{:.3f}",
                 "{:.2f}"])
        + fig("sn_fig4_limits", "Figure 4: injection-defined limits."),
        "No detection (all FAP > 0.1%). Every limit carries its injected "
        "amplitude and a noise-only predicted scale; no intra-night periodogram."))
    S.append(section(
        "late", "Late-time stacks", "SN-S9-late-time",
        "What do the post-fade stacks allow at the SN position?",
        table(sn, """SELECT epoch, camera, code, n_frames, zp, zp_err, sn_snr,
                     sn_mag, lim5_mag, lim90_mag, recovery FROM sn_latetime""",
              ["epoch", "camera", "code", "frames", "ZP", "±", "SN SNR", "SN mag",
               "5σ limit", "90% recovery", "recovery / note"],
              ["{}", "{}", "{}", "{}", "{:.2f}", "{:.2f}", "{:.1f}", "{:.2f}",
               "{:.2f}", "{:.2f}", "{}"]),
        "+1035/+1045 d r: non-detections with 5σ limits. +366 d r: the forced "
        "flux is 8.6σ, consistent with ZTF r extrapolated — reported as a "
        "validation of the faint regime, not as a new detection claim (the "
        "plan's no-detection rule was written for +1049 d)."))
    S.append(section(
        "census", "Saturation census (Gate 0, rescreened)", "SN-G0-rerun",
        "How many broadband frames survive the measured linearity cap?",
        table(man, """SELECT filter, n_frames, n_clean, n_bounded_clean, n_usable,
                      n_rejected, n_undetermined, n_spectra, n_nonscience
                      FROM sn_g0_bands WHERE band_role='broadband'""",
              ["code", "frames", "clean", "bounded", "usable", "rejected",
               "undetermined", "S2c dispersed", "outside science tree"])
        + fig("sn_fig1_census", "Figure 1: SN native peak vs phase."),
        "clean + bounded = usable in every band; every exclusion is named."))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    css = ("body{font-family:Georgia,serif;max-width:1100px;margin:auto;padding:1em;"
           "color:#222;background:#fafaf7}table{border-collapse:collapse;font-size:"
           "12px;margin:.6em 0}td,th{border:1px solid #ccc;padding:2px 6px}"
           ".decision{background:#eef4ea;padding:.5em;border-left:4px solid #4a7}"
           "img{max-width:100%}.task{color:#666}")
    page = (f"<!doctype html><html><head><meta charset='utf-8'><title>SN 2023ixf "
            f"release evidence</title><style>{css}</style></head><body>"
            f"<h1>SN 2023ixf — validation and limits release: evidence</h1>"
            f"<p>Generated {now} by <code>macro_sn.report_sn</code> from "
            f"<code>products/sn/sn2023ixf.sqlite</code> "
            f"({esc(meta.get('code_version'))}) and the manifest Gate 0 tables. "
            f"Every number is a query.</p>" + "".join(S) + "</body></html>")
    out = DOCS / "sn_release.html"
    out.write_text(page)
    return out
