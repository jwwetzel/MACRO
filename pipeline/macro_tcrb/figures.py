"""macro_tcrb.figures — the five-figure set of the T CrB manuscript.

Seat 6 (committee/reviews/2026-10-05-seat6-tcrb-abstract.md) approved the
abstract with the EW-against-phase figure merged into the centrepiece, so
the set is five, each tied to an abstract sentence:

1. ``tcrb_fig1_series.pdf``  — EW (both grisms) with same-date ARAS EWs,
   line flux, and the AAVSO B curve; orbital phase on the top axis with the
   ellipsoidal extrema marked.
2. ``tcrb_fig2_montage.pdf`` — representative spectra with the EW windows.
3. ``tcrb_fig3_lsf.pdf``     — dispersion residuals and the measured LSF.
4. ``tcrb_fig4_floor.pdf``   — the theta CrB floor (line-free window) and
   its Hα (Be) behaviour.
5. ``tcrb_fig5_aras.pdf``    — RLMT vs degraded-ARAS EW on the same dates.

Everything is read from the project database and the grism library; house
style from ``macro_core.plotstyle``.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                  # noqa: E402

from macro_core import plotstyle as ps                           # noqa: E402
from macro_tcrb import db, spec                                  # noqa: E402

MS = db.REPO / "manuscripts" / "TCrB_Monitoring"
FIG = MS / "figures"
GR_STYLE = {"hrg": (ps.ACCENT, "o"), "lrg": (ps.BAD, "s")}
MJD0 = 2400000.5


def _ephem():
    import csv
    for r in csv.DictReader(open(db.NOVELTY / "external" /
                                 "ephemerides.csv")):
        if r["ephemeris"].startswith("Munari"):
            return float(r["t0_hjd"]), float(r["period_d"])
    raise RuntimeError("Munari ephemeris missing")


def _save(fig, name: str) -> Path:
    FIG.mkdir(parents=True, exist_ok=True)
    p = FIG / name
    fig.savefig(p)
    fig.savefig(p.with_suffix(".png"), dpi=ps.PNG_DPI)
    plt.close(fig)
    return p


def fig_series(con) -> Path:
    t0, P = _ephem()
    fig, ax = plt.subplots(3, 1, figsize=(ps.COL_DOUBLE, 5.6), sharex=True,
                           gridspec_kw={"height_ratios": [2.2, 1.4, 1.2]})
    for gr, (col, mk) in GR_STYLE.items():
        a = np.array(db.q(con, """SELECT jd, ew, ew_err, line_flux,
            line_flux_err FROM tcrb_flux WHERE grism = ? ORDER BY jd""", gr),
                     float)
        if a.size == 0:
            continue
        ax[0].errorbar(a[:, 0] - MJD0, a[:, 1], a[:, 2], ms=3.2, lw=0.6,
                       capsize=0, label=f"RLMT {gr}",
                       **ps.measurement_kw(col, mk))
        ok = np.isfinite(a[:, 3])
        ax[1].errorbar(a[ok, 0] - MJD0, a[ok, 3], a[ok, 4], ms=3.2, lw=0.6,
                       capsize=0, **ps.measurement_kw(col, mk))
    ar = np.array(db.q(con, """SELECT jd_mid, ew_deg_lrg FROM tcrb_aras_ew
        WHERE date_obs BETWEEN '2025-02-15' AND '2025-06-30' AND
        ew_deg_lrg IS NOT NULL ORDER BY jd_mid"""), float)
    if ar.size:
        ax[0].plot(ar[:, 0] - MJD0, ar[:, 1], ".", color=ps.FAINT, ms=3,
                   label="ARAS (degraded, same windows)", zorder=1)
    b = np.array(db.q(con, """SELECT CAST(jd + 0.5 AS INT) AS d,
        AVG(mag), COUNT(*) FROM tcrb_aavso WHERE band = 'B' AND own = 0
        AND jd BETWEEN 2460720 AND 2460860 GROUP BY d ORDER BY d"""), float)
    if b.size:
        ax[2].plot(b[:, 0] - 0.5 - MJD0, b[:, 1], ".", color=ps.INK, ms=3,
                   label="AAVSO B (nightly median)")
        ax[2].invert_yaxis()
    ax[0].set_ylabel(r"H$\alpha$ EW (Å)")
    ax[1].set_ylabel("line flux\n(erg cm$^{-2}$ s$^{-1}$)")
    ax[2].set_ylabel("B (mag)")
    ax[2].set_xlabel("MJD")
    ax[0].legend(fontsize=6, loc="best")
    # top axis: orbital phase; ellipsoidal minima at 0.25 / 0.75 and
    # maximum at 0.5 (T0 = giant at maximum velocity, Munari et al. 2025)
    lo, hi = ax[2].get_xlim()
    sec = ax[0].secondary_xaxis(
        "top", functions=(lambda m: ((m + MJD0 - t0) / P),
                          lambda ph: ph * P + t0 - MJD0))
    sec.set_xlabel("orbital cycle (Munari et al. 2025 ephemeris)")
    k0 = np.floor((lo + MJD0 - t0) / P)
    for k in (k0, k0 + 1):
        for ph, sty in ((0.25, ":"), (0.5, "--"), (0.75, ":")):
            m = (k + ph) * P + t0 - MJD0
            if lo < m < hi:
                for a_ in ax:
                    a_.axvline(m, **ps.reference_kw(style=sty))
    fig.tight_layout()
    return _save(fig, "tcrb_fig1_series.pdf")


def fig_montage(con, g) -> Path:
    fig, ax = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 3.4))
    for k, gr in enumerate(("hrg", "lrg")):
        rows = db.q(con, """SELECT e.path, e.night, e.x_halpha FROM
            tcrb_ew_frames e WHERE e.grism = ? AND e.status = 'ok'
            ORDER BY e.jd""", gr)
        if not rows:
            continue
        pick = [rows[int(i)] for i in np.linspace(0, len(rows) - 1, 5)]
        coeffs, xref = json.loads(db.q1(g, """SELECT coeffs_json FROM
            g_dispersion WHERE grism = ? AND mech_epoch = 'ASI-pre'""", gr)),\
            db.q1(g, """SELECT x_ref FROM g_dispersion WHERE grism = ? AND
            mech_epoch = 'ASI-pre'""", gr)
        for j, (path, night, xha) in enumerate(pick):
            sf = db.q1(g, "SELECT spec_file FROM g_frames WHERE path = ?",
                       path)
            d = np.load(db.REPO / "products" / "grism" / sf)
            lam = spec.wavelength_scale(np.arange(d["flux"].size), coeffs,
                                        xref, xha)
            sel = (lam > 6300) & (lam < 6900)
            f = d["flux"][sel]
            f = f / np.nanmedian(f[(lam[sel] > 6470) & (lam[sel] < 6520)])
            ax[k].plot(lam[sel], f + 1.2 * j, lw=0.6, color=ps.INK)
            ax[k].text(6310, 1.2 * j + 1.25, night, fontsize=6)
        for lo, hi in (spec.BLUE, spec.RED):
            ax[k].axvspan(lo, hi, color=ps.GRID, zorder=0)
        h = db.q1(con, "SELECT half_width_a FROM tcrb_lsf_adopted WHERE "
                       "grism = ?", gr)
        ax[k].axvspan(spec.HALPHA - h, spec.HALPHA + h, color=ps.tint(
            ps.ACCENT, 0.8), zorder=0)
        ax[k].set_xlabel(r"wavelength (Å)")
        ax[k].set_title(gr, fontsize=8)
    ax[0].set_ylabel("normalised flux + offset")
    fig.tight_layout()
    return _save(fig, "tcrb_fig2_montage.pdf")


def fig_lsf(con, g) -> Path:
    fig, ax = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 2.8))
    for gr, (col, mk) in GR_STYLE.items():
        r = np.array(db.q(g, """SELECT wave_ref, resid_a FROM
            g_dispersion_resid WHERE grism = ? AND mech_epoch = 'ASI-pre'
            AND resid_a IS NOT NULL""", gr), float)
        if r.size:
            ax[0].plot(r[:, 0], r[:, 1], **ps.measurement_kw(col, mk, 3),
                       label=gr)
        l_ = np.array(db.q(g, """SELECT focus_offset, lsf_fwhm_a FROM g_lsf
            WHERE sample = 'tcrb' AND grism = ? AND lsf_fwhm_a > 0""", gr),
                      float)
        if l_.size:
            ax[1].plot(l_[:, 0], l_[:, 1], **ps.measurement_kw(col, mk, 3),
                       label=gr)
    ax[0].axhline(0, **ps.reference_kw())
    ax[0].set_xlabel(r"line wavelength (Å)")
    ax[0].set_ylabel(r"dispersion-fit residual (Å)")
    ax[1].set_xlabel("focus offset (focuser counts)")
    ax[1].set_ylabel(r"LSF FWHM at H$\alpha$ (Å)")
    ax[1].set_yscale("log")
    ax[0].legend(fontsize=6)
    fig.tight_layout()
    return _save(fig, "tcrb_fig3_lsf.pdf")


def fig_floor(con) -> Path:
    """theta CrB, 2025 Mode0, good frames only: (a) nightly line-free-window
    means per (grism, exposure), offset to the exposure median, with the
    adopted nightly floor; (b) nightly Halpha (hrg) at both exposures."""
    fig, ax = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 2.8))
    fl = dict(db.q(con, "SELECT grism, floor_ew FROM tcrb_floor"))
    k = 0
    for gr, (col, mk) in GR_STYLE.items():
        for ex in sorted({r[0] for r in db.q(con, """SELECT exptime FROM
                tcrb_tet_ew WHERE grism = ? AND grism_epoch = 'ASI-pre'""",
                gr)}):
            a = np.array(db.q(con, """SELECT AVG(jd), AVG(ew_linefree),
                AVG(ew_ha) FROM tcrb_tet_ew WHERE grism = ? AND exptime = ?
                AND grism_epoch = 'ASI-pre' AND continuum_fault = 0
                GROUP BY night""", gr, ex), float)
            if not a.size:
                continue
            y = a[:, 1] - np.median(a[:, 1])
            mk2 = mk if k % 2 == 0 else ("^" if gr == "hrg" else "D")
            ax[0].plot(a[:, 0] - MJD0, y, label=f"{gr} {ex:g} s",
                       **ps.measurement_kw(col if k % 2 == 0 else
                                           ps.tint(col), mk2, 3.2))
            if gr == "hrg":
                ax[1].plot(a[:, 0] - MJD0, a[:, 2], label=f"hrg {ex:g} s",
                           **ps.measurement_kw(col if k % 2 == 0 else
                                               ps.tint(col), mk2, 3.2))
            k += 1
    for gr in ("lrg", "hrg"):        # the narrower band drawn on top
        col = GR_STYLE[gr][0]
        if gr in fl:
            ax[0].axhspan(-fl[gr], fl[gr], color=ps.tint(col, 0.8),
                          zorder=0)
    ax[0].axhline(0, **ps.reference_kw())
    ax[0].set_ylabel("line-free window EW − median (Å)")
    ax[1].set_ylabel(r"$\theta$ CrB H$\alpha$ EW (Å)")
    for a_ in ax:
        a_.set_xlabel("MJD")
        a_.legend(fontsize=6)
    fig.tight_layout()
    return _save(fig, "tcrb_fig4_floor.pdf")


def fig_aras(con) -> Path:
    fig, ax = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 2.9))
    for gr, (col, mk) in GR_STYLE.items():
        a = np.array(db.q(con, """SELECT rlmt_ew, rlmt_err, aras_ew_deg,
            aras_ew_native FROM tcrb_xval WHERE grism = ?""", gr), float)
        if not a.size:
            continue
        ax[0].errorbar(a[:, 2], a[:, 0], a[:, 1], lw=0.5, capsize=0,
                       label=gr, **ps.measurement_kw(col, mk, 3))
        ax[1].plot(a[:, 3], a[:, 2] / a[:, 3] - 1,
                   **ps.measurement_kw(col, mk, 3), label=gr)
    lim = ax[0].get_xlim()
    ax[0].plot(lim, lim, **ps.reference_kw())
    ax[0].set_xlabel(r"ARAS EW, degraded to the RLMT LSF (Å)")
    ax[0].set_ylabel(r"RLMT EW, same date (Å)")
    ax[1].axhline(0, **ps.reference_kw())
    ax[1].set_xlabel(r"ARAS EW at native resolution (Å)")
    ax[1].set_ylabel("degraded / native $-$ 1")
    ax[0].legend(fontsize=6)
    fig.tight_layout()
    return _save(fig, "tcrb_fig5_aras.pdf")


def build_all() -> list[Path]:
    ps.apply("print")
    con = db.connect()
    g = db.grism_ro()
    out = []
    for fn, a in ((fig_series, (con,)), (fig_montage, (con, g)),
                  (fig_lsf, (con, g)), (fig_floor, (con,)),
                  (fig_aras, (con,))):
        try:
            out.append(fn(*a))
        except Exception as exc:                          # reported
            print(f"  figure {fn.__name__} failed: {type(exc).__name__}: "
                  f"{exc}")
    return out
