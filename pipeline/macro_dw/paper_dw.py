"""macro_dw.paper_dw — the Dwarf-Galaxy Hα manuscript's figures and numbers.

Everything the manuscript states as a number is emitted here, from
``products/dwarf/dwarf.sqlite`` (and the read-only manifest), as a LaTeX
macro in ``manuscripts/DwarfGalaxy_AGN_Survey/numbers.tex``; every figure is
drawn here from the same tables in the house print style
(``macro_core.plotstyle``).  Nothing in the prose is typed.

Figure set (Seat 6, committee/reviews/2026-10-05-seat6-dwarf-abstract.md):
Fig 4 (the Hα result) is built in both venues; Figs 2 (depth and
detectability), 3 (atlas) and 5 (NGC 5238 Hα map) only when the
pre-declared venue trigger (committee/work/dwarf/PREDECLARED.md) returns AJ.
Numbering in the files follows the approved list, not the order of
appearance in an RNAAS note.
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
from macro_dw import dwcore as core                              # noqa: E402

REPO = Path(__file__).resolve().parents[2]
MS = REPO / "manuscripts" / "DwarfGalaxy_AGN_Survey"
FIG = MS / "figures"
STACKS = REPO / "products" / "dwarf" / "stacks"
COL1, COL2 = 3.4, 7.0
CAT_STYLE = {"in": (ps.CYCLE[0], "o", "H$\\alpha$ in band"),
             "unknown": (ps.CYCLE[1], "s", "no published velocity"),
             "out": (ps.MUTED, "D", "H$\\alpha$ out of band (assumed filter)")}


def _q(con, sql, *a):
    return con.execute(sql, a).fetchall()


def _q1(con, sql, *a):
    r = con.execute(sql, a).fetchone()
    return None if r is None else r[0]


def _meta(con) -> dict:
    return dict(_q(con, "SELECT key, value FROM dw_build_meta"))


# ===========================================================================
# Figure 4 — the Hα result (both venues)
# ===========================================================================
def fig_halpha(con) -> Path:
    rows = [dict(r) for r in _q(con, "SELECT * FROM dw_halpha")]
    order = {"in": 0, "unknown": 1, "out": 2}
    rows.sort(key=lambda r: (order[r["band"]], r["field"]))
    fig, ax = plt.subplots(figsize=(COL1, 4.0))
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        col, mk, _ = CAT_STYLE[r["band"]]
        if r["verdict"] in ("no_continuum", "no_halpha_frames"):
            ax.text(0.02, yi, "no " + ("R continuum" if r["verdict"] == "no_continuum"
                                       else "H$\\alpha$ frame") + " passed QC",
                    transform=ax.get_yaxis_transform(), fontsize=5, va="center",
                    color=ps.MUTED)
            continue
        if r["verdict"] == "detected":
            ax.errorbar(r["flux_5238scale"], yi, xerr=r["flux_err"], fmt=mk,
                        color=col, mec="k", mew=0.4, ms=4.5, lw=0.8)
        else:
            ax.plot(r["flux_lim_3sig"], yi, "<", ms=5.5, mfc="none", mec=col,
                    mew=0.9)
        if r["F_pred"]:
            ax.plot(r["F_pred"], yi, "x", ms=5, color="k", mew=0.9)
            if r["verdict"] != "detected":
                ax.plot([r["F_pred"], r["flux_lim_3sig"]], [yi, yi], "-",
                        color=col, lw=0.5, alpha=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([r["candidate"] for r in rows], fontsize=6)
    for k, (cat, (col, mk, lab)) in enumerate(CAT_STYLE.items()):
        idx = [yi for yi, r in zip(y, rows) if r["band"] == cat]
        if idx:
            ax.axhspan(min(idx) - 0.5, max(idx) + 0.5, color=col, alpha=0.07, lw=0)
    ax.set_xscale("log")
    ax.set_xlabel(r"H$\alpha$ flux (erg s$^{-1}$ cm$^{-2}$, NGC 5238 scale)")
    h = [plt.Line2D([], [], ls="", marker="<", mfc="none", mec="k",
                    label="3$\\sigma$ upper limit (10$''$ aperture)"),
         plt.Line2D([], [], ls="", marker="o", color=ps.CYCLE[0], mec="k",
                    label="detection ($\\pm1\\sigma$)"),
         plt.Line2D([], [], ls="", marker="x", color="k", label="FUV/SFR-predicted flux")]
    h += [plt.Line2D([], [], ls="", marker="s", color=c, alpha=0.3, label=l)
          for c, _, l in CAT_STYLE.values()]
    ax.legend(handles=h, fontsize=5.2, loc="lower center", frameon=False,
              bbox_to_anchor=(0.42, 1.0), ncol=2, columnspacing=0.8)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    out = FIG / "dw_fig4_halpha.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


# ===========================================================================
# Figures 2, 3, 5 — AJ only
# ===========================================================================
def fig_depth(con) -> Path:
    fig, ax = plt.subplots(figsize=(COL1, 2.8))
    fields = [r[0] for r in _q(con, "SELECT DISTINCT field FROM dw_injection_contours")]
    for lev, ls in ((0.5, "-"), (0.9, ":")):
        allc = []
        for f in fields:
            c = np.array(_q(con, "SELECT re_arcsec, mu0_lim FROM dw_injection_contours "
                                 "WHERE field=? AND level=? ORDER BY re_arcsec", f, lev),
                         dtype=float)
            allc.append(c[:, 1])
            ax.plot(c[:, 0], c[:, 1], ls, color=ps.MUTED, lw=0.4, alpha=0.5)
        med = np.nanmedian(np.array(allc), axis=0)
        ax.plot(c[:, 0], med, ls, color="k", lw=1.2,
                label=f"{int(lev * 100)}% recovery (median field)")
    s = [dict(r) for r in _q(con, "SELECT * FROM dw_sersic WHERE sersic_ok=1")]
    hal = {r["field"]: r["band"] for r in _q(con, "SELECT field, band FROM dw_halpha")}
    for r in s:
        col, mk, _ = CAT_STYLE.get(hal.get(r["field"], "out"), CAT_STYLE["out"])
        ax.plot(r["re_arcsec"], r["mu0_L"], mk, color=col, mec="k", mew=0.4, ms=4)
    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel(r"$r_e$ (arcsec)")
    ax.set_ylabel(r"$\mu_{0,L}$ (mag arcsec$^{-2}$, $r$-tied)")
    ax.legend(fontsize=5.5, loc="lower left")
    out = FIG / "dw_fig2_depth.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


def _cut(name, ra, dec, half_px):
    from astropy.io import fits
    from astropy.wcs import WCS
    p = STACKS / f"{name}.fits"
    if not p.exists():
        return None
    with fits.open(p) as h:
        d, w = np.asarray(h[0].data, float), WCS(h[0].header)
    x, y = (int(round(float(v))) for v in w.all_world2pix(ra, dec, 0))
    return d[y - half_px:y + half_px + 1, x - half_px:x + half_px + 1]


def fig_atlas(con, radec: dict) -> Path:
    rows = [dict(r) for r in _q(con, "SELECT * FROM dw_halpha ORDER BY field")]
    n = len(rows)
    fig, axs = plt.subplots(n, 3, figsize=(COL1, 1.05 * n))
    half = int(45 / 0.54)
    for i, r in enumerate(rows):
        ra, de = radec[r["field"]]
        for j, (suf, lab) in enumerate((("L", "L"), ("R", "R"),
                                        ("HaSub", "H$\\alpha$$-$cont"))):
            ax = axs[i, j]
            c = _cut(f"{r['field']}_{suf}", ra, de, half)
            ax.set_xticks([])
            ax.set_yticks([])
            if c is None:
                continue
            lo, hi = np.nanpercentile(c, [5, 99.5 if j < 2 else 99.0])
            ax.imshow(c, origin="lower", cmap=ps.IMAGE_GREY + "_r", vmin=lo, vmax=hi)
            circ = plt.Circle((half, half), 10 / 0.54, fill=False, color=ps.CYCLE[0],
                              lw=0.5)
            ax.add_patch(circ)
            if i == 0:
                ax.set_title(lab, fontsize=6)
        axs[i, 0].set_ylabel(r["candidate"], fontsize=5)
        axs[i, 2].text(0.97, 0.05, r["verdict"].replace("_", " "), fontsize=4.5,
                       ha="right", transform=axs[i, 2].transAxes, color=ps.CYCLE[1])
    fig.subplots_adjust(wspace=0.03, hspace=0.05, left=0.2, right=0.99,
                        top=0.98, bottom=0.01)
    out = FIG / "dw_fig3_atlas.pdf"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return out


def fig_n5238(con, radec) -> Path:
    cal = dict(_q(con, "SELECT * FROM dw_ha_cal")[0])
    ra, de = radec["NGC5238"]
    half = int(150 / 0.54)
    c = _cut("NGC5238_HaSub", ra, de, half)
    cl = _cut("NGC5238_L", ra, de, half)
    fig, axs = plt.subplots(1, 2, figsize=(COL1, 1.8))
    for ax, im, lab in ((axs[0], cl, "L"), (axs[1], c, "H$\\alpha$ $-$ continuum")):
        lo, hi = np.nanpercentile(im, [5, 99.7])
        ext = np.array([-half, half, -half, half]) * 0.54
        ax.imshow(im, origin="lower", cmap=ps.IMAGE_GREY + "_r", vmin=lo, vmax=hi, extent=ext)
        ax.add_patch(plt.Circle((0, 0), cal["r_cal_arcsec"], fill=False,
                                color=ps.CYCLE[0], lw=0.7))
        ax.set_title(lab, fontsize=6.5)
        ax.set_xlabel(r"$\Delta\alpha$ ($''$)", fontsize=6)
    axs[0].set_ylabel(r"$\Delta\delta$ ($''$)", fontsize=6)
    axs[1].set_yticklabels([])
    out = FIG / "dw_fig5_ngc5238.pdf"
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


def _sci(v, digits=1):
    """2.3e-15 -> 2.3\\times10^{-15} (math mode content)."""
    if v is None or not np.isfinite(v):
        return None
    e = int(np.floor(np.log10(abs(v))))
    m = v / 10 ** e
    return f"{m:.{digits}f}\\times10^{{{e}}}"


def emit_numbers(con, man) -> Path:
    N = Num()
    a = N.add
    D, M = "dwarf.sqlite", "manifest"
    meta = _meta(con)
    # ---- data -----------------------------------------------------------------
    a("DWNframesStaged", _int(_q1(con, "SELECT count(*) FROM dw_frames")), "{}", D)
    a("DWNdirect", _int(_q1(con, "SELECT count(*) FROM dw_solve")), "{}", D)
    a("DWNsolved", _int(_q1(con, "SELECT count(*) FROM dw_solve WHERE status='solved'")), "{}", D)
    a("DWNscience", _int(_q1(con, "SELECT count(*) FROM dw_disposition WHERE disposition='science'")), "{}", D)
    a("DWNqcRej", _int(_q1(con, "SELECT count(*) FROM dw_disposition WHERE disposition='qc_rejected'")), "{}", D)
    a("DWNnonlinSky", _int(_q1(con, "SELECT count(*) FROM dw_disposition WHERE qc_reasons LIKE '%nonlinear_sky%'")), "{}", D)
    a("DWNunsolved", _int(_q1(con, "SELECT count(*) FROM dw_disposition WHERE disposition='unsolved'")), "{}", D)
    a("DWNmispointed", _int(_q1(con, "SELECT count(*) FROM dw_disposition WHERE disposition='mispointed'")), "{}", D)
    a("DWwcsRms", _q1(con, "SELECT avg(rms_arcsec) FROM dw_solve WHERE status='solved'"), "{:.2f}", D)
    a("DWpixScale", _q1(con, "SELECT avg(scale_arcsec) FROM dw_solve WHERE status='solved'"), "{:.3f}", D)
    a("DWfwhmHmed", _q1(con, "SELECT avg(fwhm_frames_med) FROM dw_stacks WHERE filter='H' AND night IS NULL") * 0.54, "{:.1f}", D)
    a("DWlinCap", _int(_q(man, "SELECT value FROM detector_params WHERE era_group='AC4040 High Gain e1.054' AND quantity='linearity_cap_adu'")[0][0]), "{}", M)
    a("DWnHaFields", _q1(con, "SELECT count(*) FROM dw_halpha"), "{}", D)
    # ---- flats ----------------------------------------------------------------
    v = _q(con, "SELECT max(rms) FROM dw_flat_validation WHERE filter='L' AND test LIKE 'half%'")[0][0]
    a("DWflatLresid", 100 * v, "{:.2f}", D)
    a("DWflatRresid", 100 * _q1(con, "SELECT max(rms) FROM dw_flat_validation WHERE filter='R' AND test LIKE 'half%'"), "{:.2f}", D)
    a("DWflatHresid", 100 * _q1(con, "SELECT max(rms) FROM dw_flat_validation WHERE filter='H' AND test LIKE 'half%'"), "{:.1f}", D)
    a("DWflatLepoch", 100 * _q1(con, "SELECT rms FROM dw_flat_validation WHERE filter='L' AND test='NGC5238_FebMay_on_June'"), "{:.2f}", D)
    a("DWnFlatL", _q1(con, "SELECT n_build+n_test FROM dw_flat_validation WHERE filter='L' AND test='half_B_on_A'"), "{}", D)
    g = dict(_q(con, "SELECT 'x', rho_moonangle FROM dw_gradient_moon WHERE filter='H'"))
    a("DWgradMoonRhoH", g.get("x"), "{:+.2f}", D)
    a("DWgradMoonPH", _q1(con, "SELECT p_moonangle FROM dw_gradient_moon WHERE filter='H'"), "{:.2f}", D)
    # ---- calibration ------------------------------------------------------------
    for b in ("L", "R", "H"):
        a(f"DWcterm{b}", _q1(con, "SELECT sum(cterm_gr*n_meas)/sum(n_meas) FROM dw_zp_colour WHERE filter=?", b), "{:+.2f}", D)
    a("DWzpRmsR", _q1(con, "SELECT rms FROM dw_zp_colour WHERE filter='R' AND readoutm='High Gain'"), "{:.3f}", D)
    a("DWpsOneMax", _q1(con, "SELECT max(abs(median_ps1_minus_refcat2_r)) FROM dw_ps1_check"), "{:.3f}", D)
    cal = dict(_q(con, "SELECT * FROM dw_ha_cal")[0])
    a("DWcalRadius", cal["r_cal_arcsec"], "{:.0f}", D)
    a("DWcalNights", cal["n_nights"], "{}", D)
    a("DWcalScatter", 100 * cal["S25_night_rel_scatter"], "{:.0f}", D)
    a("DWcalPubErr", 100 * cal["S25_rel_err_pub"], "{:.0f}", D)
    a("DWWeff", cal["W_eff_A"], "{:.0f}", D)
    a("DWdlamCal", cal["dlam_ngc5238_A"], "{:.0f}", D)
    a("DWfluxPub", _sci(5.50e-13, 2), "{}", "LVGDB")
    # ---- Hα results -------------------------------------------------------------
    rows = [dict(r) for r in _q(con, "SELECT * FROM dw_halpha")]
    inb = [r for r in rows if r["band"] == "in"]
    a("DWNinBand", len(inb), "{}", D)
    a("DWNoutBand", sum(r["band"] == "out" for r in rows), "{}", D)
    a("DWNnoVel", sum(r["band"] == "unknown" for r in rows), "{}", D)
    a("DWNdetInBand", sum(r["verdict"] == "detected" for r in inb), "{}", D)
    a("DWNdetAll", sum(r["verdict"] == "detected" for r in rows), "{}", D)
    a("DWNnondetInBand", sum(r.get("flux_lim_3sig") is not None for r in inb), "{}", D)
    a("DWNinformative", sum(r["informative_limit"] or 0 for r in inb), "{}", D)
    a("DWNnoMeas", sum(r["verdict"] in ("no_continuum", "no_halpha_frames") for r in rows), "{}", D)
    det = [r for r in inb if r["verdict"] == "detected"]
    a("DWdetField", ", ".join(r["candidate"] for r in det) or None, "{}", D)
    a("DWdetFlux", _sci(det[0]["flux_5238scale"]) if det else None, "{}", D)
    a("DWdetFluxErr", _sci(det[0]["flux_err"]) if det else None, "{}", D)
    a("DWdetSnr", det[0]["snr"] if det else None, "{:.0f}", D)
    a("DWdetPred", _sci(det[0]["F_pred"]) if det and det[0]["F_pred"] else None, "{}", D)
    a("DWdetOverPred", (det[0]["flux_5238scale"] / det[0]["F_pred"])
      if det and det[0]["F_pred"] else None, "{:.2f}", D)
    sysf = float(np.hypot(cal["S25_night_rel_scatter"], cal["S25_rel_err_pub"]))
    a("DWscaleSys", 100 * sysf, "{:.0f}", D)
    a("DWdetFluxSys", _sci(det[0]["flux_5238scale"] * sysf) if det else None, "{}", D)
    a("DWNinBandSpec", sum(1 for r in inb if r.get("spec_emission")), "{}", D)
    # literature counts from the DW-N1 novelty table (the same CSV the
    # pipeline reads), never typed
    import csv
    with open(REPO / "DwarfGalaxy_AGN_Survey" / "notes" / "novelty" / "out"
              / "novelty_table.csv") as fh:
        nov = list(csv.DictReader(fh))
    a("DWNfieldsAll", len(nov), "{}", "novelty_table.csv")
    a("DWNclassified", sum(1 for r in nov if r["v_hel_kms"]), "{}", "novelty_table.csv")
    a("DWWassumed", core.W_ASSUMED_A, "{:.0f}", "dwcore")
    a("DWdetSigma", core.HA_DETECT_SIGMA, "{:.0f}", "dwcore")
    lims = [r["flux_lim_3sig"] for r in rows if r["flux_lim_3sig"]]
    a("DWlimMin", _sci(min(lims)), "{}", D)
    a("DWlimMax", _sci(max(lims)), "{}", D)
    a("DWlimMedian", _sci(float(np.median(lims))), "{}", D)
    lin = [r["flux_lim_3sig"] for r in inb if r["flux_lim_3sig"]]
    a("DWlimInMin", _sci(min(lin)) if lin else None, "{}", D)
    a("DWlimInMax", _sci(max(lin)) if lin else None, "{}", D)
    a("DWNlimBelowPred", sum(1 for r in inb if r.get("F_pred") and r.get("flux_lim_3sig")
                             and r["flux_lim_3sig"] < r["F_pred"]), "{}", D)
    for r in rows:
        key = r["field"].replace("+", "p").replace("Dw", "")
        key = "".join(ch for ch in key if ch.isalpha() or ch.isdigit())
        key = "".join("ABCDEFGHIJ"[int(ch)] if ch.isdigit() else ch for ch in key)
        a(f"DWsnr{key}", r.get("snr"), "{:+.1f}", D)
        a(f"DWlim{key}", _sci(r["flux_lim_3sig"]) if r.get("flux_lim_3sig") else None, "{}", D)
        a(f"DWratio{key}", (r["F_pred"] / r["flux_lim_3sig"])
          if r.get("F_pred") and r.get("flux_lim_3sig") else None, "{:.2f}", D)
        a(f"DWpred{key}", _sci(r["F_pred"]) if r.get("F_pred") else None, "{}", D)
    a("DWvenue", meta.get("venue"), "{}", D)
    a("DWvenueN", meta.get("venue_trigger_n"), "{}", D)
    # ---- depth / structure --------------------------------------------------------
    a("DWmuLimLmed", _q1(con, "SELECT avg(mu_lim_3sig_10as) FROM dw_depth WHERE band='L'"), "{:.1f}", D)
    a("DWmuLimHmed", _q1(con, "SELECT avg(mu_lim_3sig_10as) FROM dw_depth WHERE band='H'"), "{:.1f}", D)
    m50 = _q1(con, "SELECT avg(mu0_lim) FROM dw_injection_contours WHERE level=0.5 AND re_arcsec=10")
    a("DWmuFifty", m50, "{:.1f}", D)
    a("DWNsersic", _q1(con, "SELECT count(*) FROM dw_sersic WHERE sersic_ok=1"), "{}", D)
    a("DWNdetL", _q1(con, "SELECT count(*) FROM dw_sersic WHERE detected=1"), "{}", D)
    # ---- NGC 5548 ----------------------------------------------------------------
    a("DWzoMedRms", 100 * float(meta["zeroorder_comp_rms_median"])
      if meta.get("zeroorder_comp_rms_median") not in (None, "None") else None, "{:.0f}", D)
    a("DWzoMinRms", 100 * _q1(con, "SELECT min(comp_rms_min) FROM dw_zeroorder_nights"), "{:.1f}", D)
    a("DWNslotSix", _q1(con, "SELECT count(*) FROM dw_frames WHERE filter='6'"), "{}", D)
    a("DWNslotSixMis", _q1(con, "SELECT count(*) FROM dw_disposition WHERE filter='6' AND disposition='mispointed'"), "{}", D)
    # ---- NGC 5238 completeness (data release) -----------------------------------------
    a("DWlcNstars", meta.get("lc_n_stars"), "{}", D)
    out = MS / "numbers.tex"
    head = ["%% numbers.tex -- GENERATED FILE.  DO NOT EDIT.",
            "%% Emitted by pipeline/macro_dw/paper_dw.py (run_dw_paper.py paper)",
            f"%% built {datetime.now(timezone.utc):%Y-%m-%dT%H:%M:%SZ} from products/dwarf/dwarf.sqlite and the manifest.",
            r"\providecommand{\NumMissing}{\textbf{[NUMBER MISSING]}}", ""]
    out.write_text("\n".join(head + N.lines) + "\n")
    return out


def build(db: Path, manifest: Path, radec: dict) -> list:
    FIG.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    man = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    with ps.context("print"):
        outs = [fig_halpha(con)]
        if _meta(con).get("venue") == "AJ":
            outs += [fig_depth(con), fig_atlas(con, radec), fig_n5238(con, radec)]
    outs.append(emit_numbers(con, man))
    return outs
