"""macro_sn.paper_sn — the SN 2023ixf manuscript's figures and numbers.

Everything the manuscript states as a number is emitted here, from the two
databases (``products/sn/sn2023ixf.sqlite`` and the manifest's Gate 0
tables), as a LaTeX macro in ``manuscripts/SN2023ixf_LightCurve/numbers.tex``;
every figure is drawn here from the same tables in the house print style
(``macro_core.plotstyle``).  Nothing in the prose is typed.

Figures (the ruled five-figure cap; four are drawn because the fifth —
the flash-phase differential colour — was demoted by its predicted-excess
gate before the measurement was made):

1. ``sn_fig1_census.pdf``   — the saturation census: the SN's own native
   peak in every campaign broadband frame against phase, with the measured
   linearity cap and the clip.
2. ``sn_fig2_lightcurve.pdf`` — the PS1-tied gri light curve (per frame and
   nightly) over the published one.
3. ``sn_fig3_residuals.pdf`` — RLMT minus Li et al. (2025) per band, with
   the colour-term validity range marked.
4. ``sn_fig4_limits.pdf``   — injection-defined variability limits beside
   the noise-only scale.
"""

from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

from macro_core import plotstyle as ps                           # noqa: E402

REPO = Path(__file__).resolve().parents[2]
MS = REPO / "manuscripts" / "SN2023ixf_LightCurve"
FIG = MS / "figures"
COL1, COL2 = 3.4, 7.0
BAND_COLOR = {"G": ps.CYCLE[2], "R": ps.CYCLE[1], "I": ps.CYCLE[0]}
BAND_MARK = {"G": "o", "R": "s", "I": "D"}
PS1 = {"G": "g", "R": "r", "I": "i"}


def _q(con, sql, *a):
    return con.execute(sql, a).fetchall()


def _q1(con, sql, *a):
    r = con.execute(sql, a).fetchone()
    return None if r is None else r[0]


# ===========================================================================
# Figures
# ===========================================================================
def fig_census(man, sn) -> Path:
    cap = _q1(man, "SELECT value FROM detector_params WHERE era_group = "
                   "'AC4040 High Gain e1.054' AND quantity = 'linearity_cap_adu'")
    clip = _q1(man, "SELECT clip_adu FROM s2_ceiling_modes WHERE mode='High Gain'")
    fig, ax = plt.subplots(figsize=(COL1, 2.5))
    for code in ("G", "R", "I"):
        rows = np.array(_q(man, """SELECT phase_d, peak_adu, saturation_class
            FROM sn_g0_census WHERE epoch_role='campaign' AND filter=? AND
            tree='rawimage' AND quality='wcs' AND peak_adu IS NOT NULL""", code),
            dtype=object)
        ph, pk = rows[:, 0].astype(float), rows[:, 1].astype(float)
        ok = rows[:, 2] == "clean"
        ax.scatter(ph[ok], pk[ok], s=7, marker=BAND_MARK[code],
                   color=BAND_COLOR[code], lw=0.3, edgecolor="k",
                   label=f"{code} below cap", zorder=3)
        ax.scatter(ph[~ok], pk[~ok], s=7, marker=BAND_MARK[code],
                   facecolor="none", edgecolor=BAND_COLOR[code], lw=0.6,
                   label=f"{code} at/above cap", zorder=2)
    ax.axhline(cap, **ps.reference_kw(style="--"))
    ax.axhline(clip, **ps.reference_kw(style="-"))
    ax.text(1.0, cap * 1.04, f"linearity cap {cap:,.0f} ADU", ha="left",
            va="bottom", fontsize=6)
    ax.text(1.0, clip * 1.04, f"clip {clip:,.0f} ADU", ha="left", va="bottom",
            fontsize=6)
    ax.set_yscale("log")
    ax.set_ylim(100, 6000)
    ax.set_xlim(0.5, 51)
    ax.set_xlabel("days after first light")
    ax.set_ylabel("SN native peak (raw ADU)")
    ax.legend(ncol=3, loc="lower center", fontsize=5.2, columnspacing=0.8)
    out = FIG / "sn_fig1_census.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


def _li25(sn_db_dir: Path):
    rows = []
    for ln in (sn_db_dir / "external" / "li2025_table1.dat").read_text().splitlines():
        if ln.strip():
            rows.append((float(ln[10:17]), ln[18].strip(), float(ln[20:26]),
                         float(ln[28:33])))
    return rows


def fig_lightcurve(sn, li) -> Path:
    off = {"G": 0.0, "R": 0.6, "I": 1.2}
    fig, ax = plt.subplots(figsize=(COL1, 3.4))
    for code in ("G", "R", "I"):
        b = PS1[code]
        L = np.array([(p, m, e) for p, bb, m, e in li if bb == b and p < 55])
        ax.errorbar(L[:, 0], L[:, 1] + off[code], yerr=L[:, 2], fmt="none",
                    ecolor=ps.MUTED, elinewidth=0.4, alpha=0.6, zorder=1)
        ax.scatter(L[:, 0], L[:, 1] + off[code], s=4, color=ps.MUTED,
                   marker="x", lw=0.5, zorder=1,
                   label="Li et al. (2025)" if code == "G" else None)
        f = np.array(_q(sn, "SELECT phase_d, m_ps1 FROM sn_phot WHERE usable=1 "
                            "AND code=? AND m_ps1 IS NOT NULL", code))
        ax.scatter(f[:, 0], f[:, 1] + off[code], s=3, color=BAND_COLOR[code],
                   alpha=0.35, lw=0, zorder=2)
        n = np.array(_q(sn, "SELECT phase_d, m_ps1, err FROM sn_nightly WHERE "
                            "code=? AND m_ps1 IS NOT NULL", code))
        ax.errorbar(n[:, 0], n[:, 1] + off[code], yerr=n[:, 2], fmt=BAND_MARK[code],
                    ms=3, color=BAND_COLOR[code], mec="k", mew=0.3, elinewidth=0.6,
                    zorder=3, label=f"RLMT {b}" + (f" + {off[code]:.1f}"
                                                   if off[code] else ""))
    ax.invert_yaxis()
    ax.set_xlim(3, 52)
    ax.set_xlabel("days after first light")
    ax.set_ylabel("PS1 AB magnitude (+ offset)")
    ax.legend(loc="lower left", fontsize=6)
    out = FIG / "sn_fig2_lightcurve.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


def fig_residuals(sn) -> Path:
    fig, axs = plt.subplots(3, 1, figsize=(COL1, 4.2), sharex=True)
    for ax, code in zip(axs, ("G", "R", "I")):
        r = np.array(_q(sn, "SELECT phase_d, resid, err, colour_in_range FROM "
                            "sn_resid WHERE code=? ORDER BY phase_d", code))
        s = _q(sn, "SELECT offset, offset_err, rms FROM sn_resid_summary WHERE "
                   "code=? AND subset='all'", code)[0]
        out_r = r[:, 3] == 0
        if out_r.any():
            ax.axvspan(r[out_r, 0].min() - 0.5, r[out_r, 0].max() + 0.5,
                       color=ps.tint(ps.MUTED, 0.8), lw=0, zorder=0)
        ax.errorbar(r[:, 0], r[:, 1] * 1e3, yerr=r[:, 2] * 1e3,
                    fmt=BAND_MARK[code], ms=3, color=BAND_COLOR[code], mec="k",
                    mew=0.3, elinewidth=0.6, zorder=3)
        ax.axhline(0, **ps.reference_kw(style=":"))
        ax.axhline(s[0] * 1e3, **ps.reference_kw(color=BAND_COLOR[code],
                                                 style="--"))
        if code == "G":
            z = np.array(_q(sn, "SELECT phase_d, resid, ztf_err FROM sn_resid_ztf"))
            if len(z):
                ax.scatter(z[:, 0], z[:, 1] * 1e3, marker="^", s=14,
                           facecolor="none", edgecolor="k", lw=0.6,
                           label="vs ZTF g (alerts)", zorder=4)
                ax.legend(loc="lower right", fontsize=5.2)
        ax.text(0.02, 0.92, f"{PS1[code]}: offset {s[0]*1e3:+.0f}$\\pm${s[1]*1e3:.0f}"
                f" mmag, rms {s[2]*1e3:.0f} mmag", transform=ax.transAxes,
                va="top", fontsize=6)
    axs[1].set_ylabel("RLMT $-$ Li et al. (2025), PS1 (mmag)")
    axs[2].text(0.98, 0.04, "shaded: SN colour outside the colour-term training range",
                transform=axs[2].transAxes, ha="right", va="bottom", fontsize=5.2)
    axs[-1].set_xlabel("days after first light")
    out = FIG / "sn_fig3_residuals.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


def fig_limits(sn) -> Path:
    sp_ticks = {"bump": [1, 2, 4], "sine": [2.5, 5, 10, 20]}
    fig, axs = plt.subplots(1, 2, figsize=(COL2 * 0.62, 2.4), sharey=True)
    for ax, shape, lab in ((axs[0], "bump", "Gaussian bump $\\sigma_t$ (d)"),
                           (axs[1], "sine", "sinusoid period (d)")):
        for k, code in enumerate(("G", "R", "I")):
            r = np.array(_q(sn, "SELECT scale_d, a90_mag, a_pred_mag FROM "
                                "sn_var_limits WHERE code=? AND shape=? "
                                "ORDER BY scale_d", code, shape), dtype=float)
            x = r[:, 0] * (1 + 0.06 * (k - 1))
            ok = np.isfinite(r[:, 1])
            ax.plot(x[ok], r[ok, 1] * 1e3, BAND_MARK[code], ms=4,
                    color=BAND_COLOR[code], mec="k", mew=0.3,
                    label=f"{PS1[code]}: 90% recovery")
            if (~ok).any():
                ax.plot(x[~ok], np.full((~ok).sum(), 450.0), "v", ms=4,
                        mfc="none", mec=BAND_COLOR[code], mew=0.7)
            ax.plot(r[:, 0], r[:, 2] * 1e3, "-", lw=0.7, color=BAND_COLOR[code],
                    alpha=0.6)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ticks = sp_ticks[shape]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.minorticks_off()
        ax.set_xlabel(lab)
    axs[0].set_ylabel("amplitude (mmag)")
    axs[0].legend(fontsize=5.5, loc="lower left")
    axs[1].text(0.03, 0.97, "lines: noise-only scale\n"
                "open: 90% not reached\nby 400 mmag",
                transform=axs[1].transAxes, ha="left", va="top", fontsize=5.2)
    axs[0].set_ylim(10, 1500)
    out = FIG / "sn_fig4_limits.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


# ===========================================================================
# Numbers
# ===========================================================================
class Num:
    """Accumulates \\newcommand lines with their provenance."""

    def __init__(self):
        self.lines = []

    def add(self, name, value, fmt="{}", src=""):
        if value is None or (isinstance(value, float) and not np.isfinite(value)):
            txt = r"\NumMissing"
        else:
            txt = fmt.format(value)
        self.lines.append(f"\\newcommand{{\\{name}}}{{{txt}}}  % [{src}]")


def _int(v):
    return f"{int(v):,}".replace(",", r"\,")


def emit_numbers(man, sn) -> Path:
    N = Num()
    S, G = "sn2023ixf.sqlite", "manifest Gate 0"
    a = N.add
    v = dict(_q(sn, "SELECT key, value FROM sn_build_meta"))
    a("SNNumCampFrames", _int(_q1(man, "SELECT count(*) FROM sn_g0_frames WHERE epoch_role='campaign'")), "{}", G)
    a("SNNumCampNights", _q1(man, "SELECT count(DISTINCT night) FROM sn_g0_frames WHERE epoch_role='campaign'"), "{}", G)
    b = _q(man, """SELECT sum(n_frames), sum(n_usable), sum(n_clean),
                  sum(n_bounded_clean), sum(n_rejected), sum(n_undetermined),
                  sum(n_nonscience) FROM sn_g0_bands WHERE band_role='broadband'""")[0]
    for nm, x in zip(("Broad", "Usable", "Clean", "Bounded", "Rejected", "Undet",
                      "NonSci"), b):
        a(f"SNNumG{nm}", _int(x), "{}", G)
    a("SNLinCap", _int(_q1(man, "SELECT value FROM detector_params WHERE era_group='AC4040 High Gain e1.054' AND quantity='linearity_cap_adu'")), "{}", "detector_params")
    a("SNClip", _int(_q1(man, "SELECT clip_adu FROM s2_ceiling_modes WHERE mode='High Gain'")), "{}", "s2_ceiling_modes")
    a("SNGain", _q1(man, "SELECT value FROM detector_params WHERE era_group='AC4040 High Gain e1.054' AND quantity='gain_e_per_adu'"), "{:.3f}", "detector_params")
    a("SNGainErr", _q1(man, "SELECT uncertainty FROM detector_params WHERE era_group='AC4040 High Gain e1.054' AND quantity='gain_e_per_adu'"), "{:.3f}", "detector_params")
    a("SNFirstClean", _q1(man, "SELECT min(first_clean_phase_d) FROM sn_g0_bands WHERE band_role='broadband'"), "{:.1f}", G)
    a("SNNumSlotSix", _q1(man, "SELECT count(*) FROM sn_g0_frames WHERE filter='6'"), "{}", G)
    # astrometry
    a("SNNumSolvedNew", _q1(sn, "SELECT count(*) FROM sn_frames WHERE status='measured' AND wcs_source!='header' AND filter IN ('G','R','I') AND epoch_role='campaign'"), "{}", S)
    a("SNNumUnsolvable", _q1(sn, "SELECT count(*) FROM sn_frames WHERE status!='measured' AND filter IN ('G','R','I') AND epoch_role='campaign' AND tree='rawimage'"), "{}", S)
    a("SNWcsRms", _q1(sn, "SELECT avg(wcs_rms_arcsec) FROM sn_frames WHERE status='measured' AND epoch_role='campaign'"), "{:.2f}", S)
    # calibration
    a("SNScintC", float(v["scint_C"]), "{:.2f}", S)
    for code in ("G", "R", "I"):
        r = _q(sn, "SELECT * FROM sn_cal_sets WHERE role='campaign' AND code=?", code)[0]
        d = dict(zip([c[1] for c in sn.execute("PRAGMA table_info(sn_cal_sets)")], r))
        a(f"SNCterm{code}", d["cterm"], "{:+.3f}", S)
        a(f"SNCtermErr{code}", d["cterm_err"], "{:.3f}", S)
        a(f"SNKtwo{code}", d["k2"], "{:+.3f}", S)
        a(f"SNFloor{code}", d["floor_mag"] * 1e3, "{:.0f}", S)
        a(f"SNNens{code}", d["n_ensemble"], "{}", S)
        a(f"SNNchk{code}", d["n_check"], "{}", S)
        a(f"SNChk{code}", d["check_chi2nu"], "{:.2f}", S)
        a(f"SNChkDof{code}", _int(d["check_dof"]), "{}", S)
        a(f"SNTieRms{code}", d["tie_rms"] * 1e3, "{:.0f}", S)
        a(f"SNColLo{code}", d["colour_lo"], "{:.2f}", S)
        a(f"SNColHi{code}", d["colour_hi"], "{:.2f}", S)
        a(f"SNFramesUsed{code}", d["n_frames_used"], "{}", S)
        a(f"SNFramesCal{code}", d["n_frames"], "{}", S)
        a(f"SNSNFloor{code}", _q1(sn, "SELECT floor_mag FROM sn_sn_floor WHERE code=?", code) * 1e3, "{:.0f}", S)
        a(f"SNUse{code}", _q1(sn, "SELECT count(*) FROM sn_phot WHERE usable=1 AND code=?", code), "{}", S)
        a(f"SNUseN{code}", _q1(sn, "SELECT count(DISTINCT night) FROM sn_phot WHERE usable=1 AND code=?", code), "{}", S)
        a(f"SNOverlap{code}", _q1(sn, "SELECT avg(dm) FROM (SELECT dm FROM sn_overlap WHERE code=? ORDER BY dm LIMIT 2 - (SELECT count(*) FROM sn_overlap WHERE code=?) % 2 OFFSET (SELECT (count(*) - 1) / 2 FROM sn_overlap WHERE code=?))", code, code, code) * 1e3, "{:.1f}", S)
        a(f"SNOverlapN{code}", _q1(sn, "SELECT count(*) FROM sn_overlap WHERE code=?", code), "{}", S)
        ratios = [x[0] for x in _q(sn, "SELECT ratio FROM sn_cal_scint WHERE role='campaign' AND code=? AND n_points>=30", code)]
        a(f"SNScintLo{code}", min(ratios), "{:.2f}", S)
        a(f"SNScintHi{code}", max(ratios), "{:.2f}", S)
        # crosswalk
        a(f"SNXband{code}", _q1(sn, "SELECT ps1_band FROM sn_crosswalk WHERE role='campaign' AND code=? AND identified=1", code), "{}", S)
        a(f"SNXqhy{code}", _q1(sn, "SELECT cterm FROM sn_crosswalk WHERE role='template_post' AND camera='QHY600' AND code=? AND identified=1", PS1[code]), "{:+.3f}", S)
        # residuals
        rs = _q(sn, "SELECT n, offset, offset_err, rms, chi2nu, dof FROM sn_resid_summary WHERE code=? AND subset='all'", code)[0]
        for nm, x, f in zip(("N", "Off", "OffErr", "Rms", "Chi", "Dof"), rs,
                            ("{}", "{:+.0f}", "{:.0f}", "{:.0f}", "{:.2f}", "{}")):
            a(f"SNRes{nm}{code}", x * 1e3 if nm in ("Off", "OffErr", "Rms") else x, f, S)
        ri = _q(sn, "SELECT offset, offset_err, rms, n FROM sn_resid_summary WHERE code=? AND subset='colour inside training range'", code)
        if ri:
            a(f"SNResInOff{code}", ri[0][0] * 1e3, "{:+.0f}", S)
            a(f"SNResInRms{code}", ri[0][2] * 1e3, "{:.0f}", S)
            a(f"SNResInN{code}", ri[0][3], "{}", S)
        ru = _q(sn, "SELECT offset FROM sn_resid_summary WHERE code=? AND subset LIKE 'all, published%'", code)
        a(f"SNResRawOff{code}", ru[0][0] * 1e3, "{:+.0f}", S)
        # variability
        o = _q(sn, "SELECT n_nights, chi2nu_trend, bump_snr, bump_fap, sine_chi2, sine_fap, knot_spacing_d, bump_thr, sine_thr FROM sn_var_observed WHERE code=?", code)[0]
        for nm, x, f in zip(("Nights", "TrendChi", "BumpSnr", "BumpFap", "SineChi",
                             "SineFap", "Knots", "BumpThr", "SineThr"), o,
                            ("{}", "{:.2f}", "{:.1f}", "{:.2f}", "{:.1f}", "{:.2f}",
                             "{:.0f}", "{:.1f}", "{:.1f}")):
            a(f"SNVar{nm}{code}", x, f, S)
        sl = [x[0] for x in _q(sn, "SELECT a90_mag FROM sn_var_limits WHERE code=? AND shape='sine' AND scale_d<=5 AND a90_mag IS NOT NULL", code)]
        a(f"SNSineTen{code}", _q1(sn, "SELECT a90_mag FROM sn_var_limits WHERE code=? AND shape='sine' AND scale_d=10", code) * 1e3, "{:.0f}", S)
        a(f"SNSineLo{code}", min(sl) * 1e3, "{:.0f}", S)
        a(f"SNSineHi{code}", max(sl) * 1e3, "{:.0f}", S)
        bl = [x[0] for x in _q(sn, "SELECT a90_mag FROM sn_var_limits WHERE code=? AND shape='bump' AND a90_mag IS NOT NULL", code)]
        a(f"SNBumpBest{code}", (min(bl) * 1e3) if bl else None, "{:.0f}", S)
        a(f"SNPred{code}", _q1(sn, "SELECT a_pred_mag FROM sn_var_limits WHERE code=? AND shape='sine' AND scale_d=5", code) * 1e3, "{:.0f}", S)
    z = _q(sn, "SELECT resid FROM sn_resid_ztf ORDER BY phase_d")
    a("SNNumZtf", len(z), "{}", S)
    a("SNZtfLo", min(x[0] for x in z) * 1e3 if z else None, "{:+.0f}", S)
    a("SNZtfHi", max(x[0] for x in z) * 1e3 if z else None, "{:+.0f}", S)
    a("SNNumLi", _q1(sn, "SELECT count(*) FROM sn_resid"), "{}", S)
    a("SNPsfIndPass", _q1(sn, "SELECT sum(passes) FROM sn_psfcheck WHERE s2c_class='indeterminate'"), "{}", S)
    a("SNPsfIndN", _q1(sn, "SELECT count(*) FROM sn_psfcheck WHERE s2c_class='indeterminate'"), "{}", S)
    a("SNPsfDirPass", _q1(sn, "SELECT sum(passes) FROM sn_psfcheck WHERE s2c_class='direct'"), "{}", S)
    a("SNPsfDirN", _q1(sn, "SELECT count(*) FROM sn_psfcheck WHERE s2c_class='direct'"), "{}", S)
    a("SNPsfMaxElong", float(v["psf_max_elong"]), "{:.1f}", S)
    a("SNNumIndetBroad", _q1(man, "SELECT count(*) FROM sn_g0_census WHERE epoch_role='campaign' AND band_role='broadband' AND dispersion_class='indeterminate'"), "{}", G)
    a("SNNumDispBroad", _q1(man, "SELECT count(*) FROM sn_g0_census WHERE epoch_role='campaign' AND band_role='broadband' AND dispersion_class='dispersed'"), "{}", G)
    a("SNNumUsePhot", _q1(sn, "SELECT count(*) FROM sn_phot WHERE usable=1 AND code IN ('G','R','I')"), "{}", S)
    a("SNNumCloud", _q1(sn, "SELECT count(*) FROM sn_phot WHERE code IN ('G','R','I') AND census_class IN ('clean','bounded_clean') AND exclusion LIKE '%cloud%'"), "{}", S)
    a("SNNumNoWcs", _q1(sn, "SELECT count(*) FROM sn_phot WHERE code IN ('G','R','I') AND census_class IN ('clean','bounded_clean') AND exclusion LIKE '%no WCS%'"), "{}", S)
    a("SNNumNotValid", _q1(sn, "SELECT count(*) FROM sn_phot WHERE code IN ('G','R','I') AND census_class IN ('clean','bounded_clean') AND exclusion LIKE '%not validated%' AND exclusion NOT LIKE '%cloud%' AND exclusion NOT LIKE '%no WCS%'"), "{}", S)
    # peak
    for b in ("g", "r", "i"):
        r = _q(sn, "SELECT t_peak_d, t_lo, t_hi FROM sn_peak WHERE band=?", b)[0]
        a(f"SNPeak{b}", r[0], "{:.1f}", "sn_peak (Li+25 data)")
        a(f"SNPeakLo{b}", r[1], "{:.1f}", "sn_peak")
        a(f"SNPeakHi{b}", r[2], "{:.1f}", "sn_peak")
    # gate
    w = dict((r[0], r[1:]) for r in _q(sn, "SELECT code, width_a, width_mad, n_nights FROM sn_nb_widths"))
    a("SNWidthH", w["H"][0], "{:.0f}", S); a("SNWidthHErr", w["H"][1], "{:.0f}", S)
    a("SNWidthS", w["1"][0], "{:.0f}", S); a("SNWidthSErr", w["1"][1], "{:.0f}", S)
    g = _q(sn, "SELECT epoch_d, max(dm_upper), min(dm_upper), max(sig_rep), max(sig_loc), min(sig_loc), max(sig_nb) FROM sn_excess_gate GROUP BY epoch_d ORDER BY epoch_d")
    for (ep, dmx, dmn, srep, sl_hi, sl_lo, snb), tag in zip(g, ("A", "B", "C")):
        a(f"SNGateDm{tag}", dmx * 1e3, "{:.0f}", S)
        a(f"SNGateDmMin{tag}", dmn * 1e3, "{:.0f}", S)
        a(f"SNGateLoc{tag}", sl_hi * 1e3, "{:.0f}", S)
    a("SNGateRep", g[0][3] * 1e3, "{:.0f}", S)
    a("SNGateVerdict", v.get("excess_gate_verdict"), "{}", S)
    # late time
    lt = _q(sn, "SELECT epoch, code, sn_snr, sn_mag, lim5_mag, lim90_mag, zp_err, n_frames FROM sn_latetime WHERE lim5_mag IS NOT NULL ORDER BY epoch")
    for ep, code, snr, mg, l5, l90, zpe, nfr in lt:
        tag = {"+366 d": "A", "+1035 d": "B", "+1045 d": "C"}[ep] + code
        a(f"SNLateMagErr{tag}", float(np.hypot(1.0857 / snr, zpe)) if snr and snr > 0 else None, "{:.2f}", S)
        a(f"SNLateNfr{tag}", nfr, "{}", S)
        a(f"SNLateSnr{tag}", snr, "{:.1f}", S)
        a(f"SNLateMag{tag}", mg, "{:.2f}", S)
        a(f"SNLateLim{tag}", l5, "{:.2f}", S)
        a(f"SNLateNinety{tag}", l90, "{:.2f}", S)
    # templates
    a("SNTplHFwhm", _q1(sn, "SELECT fwhm_arcsec FROM sn_template_table WHERE night='2023-05-04' AND filter='H'"), "{:.1f}", S)
    a("SNTplRFwhm", _q1(sn, "SELECT fwhm_arcsec FROM sn_template_table WHERE night='2023-05-04' AND filter='R'"), "{:.1f}", S)
    a("SNTplRFoc", _q1(sn, "SELECT focus_offset FROM sn_template_table WHERE night='2023-05-04' AND filter='R'"), "{:+.0f}", S)
    a("SNTplGFoc", _q1(sn, "SELECT focus_offset FROM sn_template_table WHERE night='2023-05-04' AND filter='G'"), "{:+.0f}", S)
    a("SNNumTplEpochs", _q1(sn, "SELECT count(DISTINCT night) FROM sn_template_table"), "{}", S)
    a("SNNumRefcat", _q1(sn, "SELECT count(DISTINCT star_id) FROM sn_star_phot"), "{}", S)
    used = _q(sn, """SELECT c.code, c.exptime, c.ratio FROM sn_cal_scint c
        WHERE c.role='campaign' AND c.code IN ('G','R','I') AND c.n_points >= 30
        AND EXISTS (SELECT 1 FROM sn_phot p WHERE p.usable=1 AND p.code=c.code
                    AND round(p.exptime,1)=round(c.exptime,1))""")
    a("SNScintUsedLo", min(x[2] for x in used), "{:.2f}", S)
    a("SNScintUsedHi", max(x[2] for x in used), "{:.2f}", S)
    a("SNScintUsedN", len(used), "{}", S)
    bad = _q(sn, """SELECT code, exptime, ratio FROM sn_cal_scint WHERE
        role='campaign' AND code IN ('G','R','I') AND n_points >= 30 AND
        (ratio < 0.8 OR ratio > 1.25)""")
    a("SNScintBadN", len(bad), "{}", S)
    a("SNScintBadList", ", ".join(f"{c} {t:g}\\,s ({r:.2f})" for c, t, r in bad), "{}", S)
    allsl = [x[0] for x in _q(sn, "SELECT a90_mag FROM sn_var_limits WHERE shape='sine' AND scale_d<=5 AND a90_mag IS NOT NULL")]
    a("SNSineAllLo", min(allsl) * 1e3, "{:.0f}", S)
    a("SNSineAllHi", max(allsl) * 1e3, "{:.0f}", S)
    a("SNPeakVTeja", _q1(sn, "SELECT value FROM sn_literature WHERE key LIKE 'peak_V_Teja%'"), "{:g}", "Teja et al. 2023")
    a("SNPeakVLi", _q1(sn, "SELECT value FROM sn_literature WHERE key LIKE 'peak_V_Li%'"), "{:g}", "Li et al. 2025")
    a("SNHaMaxPhase", _q1(sn, "SELECT value FROM sn_literature WHERE key='f_halpha_max_phase'"), "{:g}", "Bostroem et al. 2023")
    a("SNHaFade", _q1(sn, "SELECT value FROM sn_literature WHERE key='f_halpha_fade'"), "{:g}", "Smith et al. 2023")
    ep = [r[0] for r in _q(sn, "SELECT DISTINCT epoch_d FROM sn_excess_gate ORDER BY epoch_d")]
    for e, tag in zip(ep, ("A", "B", "C")):
        a(f"SNGateEpoch{tag}", e, "{:.1f}", S)
    a("SNTzeroMJD", _q1(sn, "SELECT value FROM sn_literature WHERE key='t0_mjd_li25'"), "{:.3f}", "Li et al. 2025")
    # method constants, read from the code that applies them
    import importlib.util
    from macro_sn import snio, snphot as sp
    spec = importlib.util.spec_from_file_location(
        "run_sn_photometry", REPO / "pipeline" / "scripts" / "run_sn_photometry.py")
    R = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(R)
    M = "code constant"
    a("SNAperFac", snio.APER_FWHM, "{:g}", M)
    a("SNAnnIn", snio.ANNULUS_FWHM[0], "{:g}", M)
    a("SNAnnOut", snio.ANNULUS_FWHM[1], "{:g}", M)
    a("SNFapPct", 100 * sp.FAP_LEVEL, "{:g}", M)
    a("SNRecovPct", 100 * sp.RECOVERY_LEVEL, "{:g}", M)
    a("SNInFocus", R.IN_FOCUS_FACTOR, "{:g}", M)
    a("SNEnsSnr", R.ENS_MIN_SNR, "{:g}", M)
    a("SNMinEns", R.MIN_ENSEMBLE, "{}", M)
    a("SNCloud", R.CLOUD_MAG, "{:g}", M)
    a("SNMaxFwhm", R.MAX_FWHM_PX, "{:g}", M)
    a("SNNboot", _int(R.N_BOOT), "{}", M)
    a("SNNinj", R.N_INJ, "{}", M)
    a("SNGateFactor", R.GATE_FACTOR, "{:g}", M)
    a("SNCheckMod", sp.CHECK_MODULUS, "{}", M)
    a("SNScintWinLo", sp.SCINT_MATCH_WINDOW[0], "{:g}", M)
    a("SNScintWinHi", sp.SCINT_MATCH_WINDOW[1], "{:g}", M)
    a("SNCalRlo", sp.CAL_R_RANGE[0], "{:g}", M)
    a("SNCalRhi", sp.CAL_R_RANGE[1], "{:g}", M)
    a("SNMinMatch", sp.WCS_MIN_MATCH, "{}", M)
    a("SNMinMatchShift", sp.WCS_MIN_MATCH_SHIFT, "{}", M)
    a("SNNoiseAp", R.N_NOISE_AP, "{}", M)
    a("SNSubHalfArcmin", (2 * R.SUB_HALF + 1) * snio.GRID_SCALE_ARCSEC / 60, "{:.1f}", M)
    a("SNNrepStars", float(v.get("nb_rep_n", "nan")), "{:.0f}", S)
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    head = ["%% numbers.tex -- GENERATED FILE.  DO NOT EDIT.",
            "%% Emitted by pipeline/macro_sn/paper_sn.py (run_sn_photometry.py paper)",
            f"%% built {now} from products/sn/sn2023ixf.sqlite and the manifest's Gate 0 tables.",
            r"\providecommand{\NumMissing}{\textbf{[NUMBER MISSING]}}", ""]
    out = MS / "numbers.tex"
    out.write_text("\n".join(head + N.lines) + "\n")
    return out


def build(sn_db: Path, manifest: Path) -> list:
    FIG.mkdir(parents=True, exist_ok=True)
    sn = sqlite3.connect(f"file:{sn_db}?mode=ro", uri=True)
    man = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    li = _li25(sn_db.parent)
    with ps.context("print"):
        outs = [fig_census(man, sn), fig_lightcurve(sn, li), fig_residuals(sn),
                fig_limits(sn)]
    outs.append(emit_numbers(man, sn))
    return outs
