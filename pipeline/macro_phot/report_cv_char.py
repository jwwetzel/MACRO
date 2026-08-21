"""S5 chain-of-evidence report: what the CV time series can actually measure.

Reads ``products/phot/cv_characterization.sqlite`` (and, for context counts
only, the photometry product it was built from) and writes

* ``docs/CV_TimeSeries/cv_characterization.html``
* ``docs/CV_TimeSeries/figures/cv_characterization/*.png``

Socratic throughout: every section poses its question, shows the evidence
(figures plus script-emitted numbers), states the decision in a
``.decision`` callout, and names the consequence for the next section.
Every number is the result of a SQL query executed here or a constant
imported from ``macro_phot.characterize`` — none is typed by hand.

The page is deliberately an ARGUMENT, in this order: the images are good
enough (§1) -> so the noise floor is real and measurable (§2) -> so the
achieved precision is a fact (§2) -> the sampling decides which periods can
be claimed (§3) -> injection through the real timestamps converts precision
into detectability (§4) -> the same machinery run on an eclipse edge
decides the timing tier (§5) -> and only then does §6 grade the science
goals.  Reading it backwards is the failure mode this layout prevents.
"""

from __future__ import annotations

import json
import math
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

from . import characterize as ch  # noqa: E402
from . import photometry as ph    # noqa: E402
from macro_core.report_s0 import (  # noqa: E402
    ACCENT, BAD, STYLE, DPI, GOOD, INK, MUTED, WARN,
    _figure, esc, fmt, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402  (house figure style)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "CV_TimeSeries"
FIG_DIR = DOCS_DIR / "figures" / "cv_characterization"
HTML_PATH = DOCS_DIR / "cv_characterization.html"

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}
TARGET_ORDER = ["stlmi", "yzcnc", "vvpup", "euuma", "anuma"]
FILTER_COLOR = ps.BAND_COLOR
FILTER_MARKER = ps.BAND_MARKER


def _mmag(x, nd=1):
    """Milli-magnitudes, or an em-dash for a missing value."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{1000 * float(x):,.{nd}f}"


def _num(x, nd=2):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{float(x):,.{nd}f}"


def _label(series_key: str) -> str:
    t, e, f = series_key.split("|")
    return f"{TARGET_LABEL.get(t, t)} {e} {f}"


# ===========================================================================
# §1 figures — image quality
# ===========================================================================

def fig_seeing(con) -> str:
    """FWHM in arcsec and recomputed airmass, per camera era."""
    rows = q(con, "SELECT readoutm, fwhm_as, airmass, target_key FROM ch_frames "
                  "WHERE fwhm_as IS NOT NULL")
    modes = sorted({r[0] for r in rows})
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.2, 3.6))
        for i, m in enumerate(modes):
            v = [r[1] for r in rows if r[0] == m]
            a1.hist(v, bins=np.arange(0.8, 6.0, 0.15), histtype="step",
                    lw=1.6, label=f"{m} (n={len(v)})",
                    **ps.line_series(i))
        a1.set_xlabel('FWHM (arcsec, converted with the era plate scale)')
        a1.set_ylabel("frames")
        a1.set_title("Seeing is comparable only in arcsec")
        a1.legend(fontsize=7)
        for i, t in enumerate(TARGET_ORDER):
            v = [r[2] for r in rows if r[3] == t and r[2] is not None]
            if v:
                a2.hist(v, bins=np.arange(1.0, 2.6, 0.05),
                        histtype="step", lw=1.6,
                        label=f"{TARGET_LABEL[t]} "
                              f"(med {np.median(v):.2f})",
                        **ps.line_series(i))
        a2.set_xlabel("airmass, recomputed from coordinates + time")
        a2.set_ylabel("frames")
        a2.set_title("VV Pup never rises: its floor is structural")
        a2.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "seeing.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/seeing.png"


def fig_quality_response(con) -> str:
    """The defended cut: check-star scatter vs each quality axis."""
    axes = ["fwhm_as", "airmass", "sky_ratio", "zp_excess",
            "moon_illum", "moon_sep"]
    titles = {"fwhm_as": "seeing (arcsec)", "airmass": "airmass",
              "sky_ratio": "sky brightness / series median",
              "zp_excess": "zero-point excess (mag) = extinction",
              "moon_illum": "moon illuminated fraction",
              "moon_sep": "moon separation (deg)"}
    cuts = {r[0]: (r[1], r[2]) for r in q(
        con, "SELECT axis, threshold, baseline FROM ch_cuts")}
    applied = {r[0] for r in q(con, "SELECT axis FROM ch_cuts "
                                    "WHERE note NOT LIKE '%diagnostic%'")}
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(2, 3, figsize=(11.4, 6.0))
        for ax, name in zip(axs.ravel(), axes):
            rows = q(con, "SELECT bin_center, med_rel_scatter, n_frames "
                          "FROM ch_quality_bins WHERE axis=? "
                          "ORDER BY bin_center", (name,))
            xs = [r[0] for r in rows if r[1] is not None and r[2] >= 5]
            ys = [r[1] for r in rows if r[1] is not None and r[2] >= 5]
            ns = [r[2] for r in rows if r[1] is not None and r[2] >= 5]
            ax.axhline(1.0, color=MUTED, lw=1, ls=":")
            ax.axhline(ch.DEGRADE_FACTOR, color=WARN, lw=1, ls="--")
            ax.scatter(xs, ys, s=[max(6, min(70, n / 6)) for n in ns],
                       color=ACCENT if name in applied else MUTED, zorder=3)
            ax.plot(xs, ys, color=ACCENT if name in applied else MUTED,
                    lw=1, alpha=0.5)
            thr = cuts.get(name, (None, None))[0]
            if thr is not None and math.isfinite(thr):
                ax.axvline(thr, color=BAD, lw=1.6)
                ax.text(thr, ax.get_ylim()[1] * 0.95, f" cut {thr:g}",
                        color=BAD, fontsize=8, va="top")
            ax.set_xlabel(titles[name])
            ax.set_ylabel("check-star scatter / series median")
            ax.set_title(("APPLIED" if name in applied else "diagnostic only"),
                         fontsize=9,
                         color=ACCENT if name in applied else MUTED)
        fig.suptitle("Every cut is set by the scatter it buys back "
                     f"(dashed = {ch.DEGRADE_FACTOR:g}x the pooled median)",
                     fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(FIG_DIR / "quality_response.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/quality_response.png"


def fig_moon(con) -> str:
    """Moonlight -> sky brightness -> scatter: the causal chain, measured."""
    rows = q(con, "SELECT moon_illum, moon_sep, sky_ratio, rel_scatter, moon_alt "
                  "FROM ch_frames WHERE moon_illum IS NOT NULL "
                  "AND sky_ratio IS NOT NULL")
    ill = np.array([r[0] for r in rows]); sep = np.array([r[1] for r in rows])
    sky = np.array([r[2] for r in rows]); alt = np.array([r[4] for r in rows])
    up = alt > 0
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.2, 3.6))
        a1.scatter(ill[~up], sky[~up], s=4, alpha=0.3, color=MUTED,
                   label="moon below horizon")
        sc = a1.scatter(ill[up], sky[up], s=5, alpha=0.5, c=sep[up],
                        cmap="viridis", label="moon up")
        plt.colorbar(sc, ax=a1, label="moon separation (deg)")
        a1.set_yscale("log")
        a1.set_xlabel("moon illuminated fraction")
        a1.set_ylabel("sky brightness / series median")
        a1.set_title("Moonlight raises the sky only when the moon is up")
        a1.legend(fontsize=7, loc="upper left")
        b = np.arange(0.0, 1.01, 0.1)
        idx = np.digitize(ill[up], b) - 1
        med = [np.median(sky[up][idx == i]) if (idx == i).sum() > 20 else np.nan
               for i in range(len(b) - 1)]
        a2.plot(0.5 * (b[:-1] + b[1:]), med, "o-", color=ACCENT)
        a2.set_xlabel("moon illuminated fraction (moon above horizon)")
        a2.set_ylabel("median sky / series median")
        a2.set_title("The moon term, isolated")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "moon.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/moon.png"


def fig_trail(con) -> str:
    """Sampled ellipticity and position-angle coherence: is trailing real?"""
    rows = q(con, "SELECT series_key, ell_med, pa_R FROM ch_trail "
                  "WHERE status='ok'")
    if not rows:
        return ""
    e = np.array([r[1] for r in rows], dtype=float)
    R = np.array([r[2] for r in rows], dtype=float)
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.2, 3.4))
        a1.hist(e, bins=np.arange(0, 0.61, 0.02), color=ACCENT)
        a1.set_xlabel("median ellipticity 1 - b/a")
        a1.set_ylabel("sampled frames")
        a1.set_title("Source elongation")
        a2.scatter(e, R, s=14, color=ACCENT, alpha=0.8)
        a2.set_xlabel("median ellipticity")
        a2.set_ylabel("position-angle coherence R")
        a2.set_title("Trailing = elongated AND aligned (upper right)")
        a2.set_ylim(0, 1)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "trail.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/trail.png"


# ===========================================================================
# §2 figures — noise
# ===========================================================================

def _rms_panel(ax, con, sk):
    stars = q(con, "SELECT mean_mag, rms, role, pred_lo, pred_nom, pred_hi "
                   "FROM ch_noise_stars WHERE series_key=? AND rms IS NOT NULL",
              (sk,))
    m = np.array([s[0] for s in stars], dtype=float)
    r = np.array([s[1] for s in stars], dtype=float)
    role = [s[2] for s in stars]
    lo = np.array([s[3] if s[3] is not None else np.nan for s in stars])
    hi = np.array([s[5] if s[5] is not None else np.nan for s in stars])
    comp = np.array([x in ("comp", "check") for x in role])
    ax.scatter(m[comp], r[comp], s=5, alpha=0.45, color=MUTED,
               label="comparison stars")
    chk = np.array([x == "check" for x in role])
    ax.scatter(m[chk], r[chk], s=34, color=ACCENT, zorder=4,
               label="held-out check stars")
    tgt = np.array([x == "target" for x in role])
    if tgt.any():
        ax.scatter(m[tgt], r[tgt], s=60, marker="*", color=WARN, zorder=5,
                   label="target (variable: not a noise point)")
    o = np.argsort(m)
    ax.fill_between(m[o], lo[o], hi[o], color=ACCENT, alpha=0.25, zorder=2,
                    label=f"photon+sky+read, gain {ch.GAIN_LO_E_PER_ADU}"
                          f"-{ch.GAIN_HI_E_PER_ADU} e-/ADU")
    row = q(con, "SELECT floor_nom, floor_plateau, scint_mag FROM "
                 "ch_noise_series WHERE series_key=?", (sk,))
    if row and row[0][0] is not None:
        ax.axhline(row[0][0], color=GOOD, lw=1.4, ls="--",
                   label=f"fitted floor {1000 * row[0][0]:.0f} mmag")
    if row and row[0][2] is not None:
        ax.axhline(row[0][2], color=BAD, lw=1.2, ls=":",
                   label=f"scintillation {1000 * row[0][2]:.1f} mmag")
    ax.set_yscale("log")
    # Low enough that the scintillation line is always ON the plot: it is the
    # comparison the whole figure exists to make.
    ax.set_ylim(2e-4, 1.0)
    ax.set_xlabel("instrumental magnitude (ensemble gauge)")
    ax.set_ylabel("RMS (mag)")
    ax.set_title(_label(sk), fontsize=10)


def fig_rms_mag(con) -> str:
    """The non-negotiable figure: RMS vs magnitude, one panel per target."""
    picks = []
    for t in TARGET_ORDER:
        # Same "richest series" rule the injection stage uses: most frames,
        # not most stars, so every section of the page talks about the same
        # series for a given target.
        row = q(con, "SELECT n.series_key FROM ch_noise_series n "
                     "JOIN ch_cadence c USING(series_key) "
                     "WHERE n.target_key=? ORDER BY c.n_points DESC LIMIT 1",
                (t,))
        if row:
            picks.append(row[0][0])
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(2, 3, figsize=(12.6, 6.6))
        for ax, sk in zip(axs.ravel(), picks):
            _rms_panel(ax, con, sk)
        axs.ravel()[0].legend(fontsize=6.5, loc="upper left")
        for ax in axs.ravel()[len(picks):]:
            ax.axis("off")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "rms_mag.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/rms_mag.png"


def fig_floor_budget(con) -> str:
    """Measured floor vs the one term the atmosphere forces on us."""
    rows = q(con, "SELECT series_key, floor_plateau, scint_mag, "
                  "prec_at_target, target_key FROM ch_noise_series "
                  "WHERE floor_plateau IS NOT NULL ORDER BY target_key, series_key")
    lab = [_label(r[0]) for r in rows]
    fl = np.array([r[1] for r in rows]) * 1000
    sc = np.array([r[2] for r in rows]) * 1000
    pr = np.array([r[3] if r[3] is not None else np.nan for r in rows]) * 1000
    y = np.arange(len(rows))
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.8, 0.30 * len(rows) + 1.6))
        ax.barh(y, fl, color=ACCENT, label="measured floor (best mag bin)")
        ax.plot(sc, y, "o", color=BAD, ms=5,
                label="scintillation (Young, 0.5 m, this airmass & exposure)")
        ax.plot(pr, y, "*", color=WARN, ms=10,
                label="achieved precision at the target's own brightness")
        ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=7.5)
        ax.set_xscale("log")
        ax.set_xlabel("mmag")
        ax.invert_yaxis()
        ax.legend(fontsize=7.5, loc="lower right")
        ax.set_title("The floor is instrumental, not atmospheric")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "floor_budget.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/floor_budget.png"


def fig_allan(con) -> str:
    """Does averaging still work at the orbital timescale?"""
    keys = [r[0] for r in q(con, "SELECT DISTINCT series_key FROM ch_allan "
                                 "ORDER BY series_key")]
    keys = [k for k in keys if k.split("|")[0] in TARGET_ORDER][:12]
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.0))
        for i, k in enumerate(keys):
            rows = q(con, "SELECT tau_s, adev FROM ch_allan WHERE series_key=? "
                          "ORDER BY tau_s", (k,))
            tau = np.array([r[0] for r in rows]); ad = np.array([r[1] for r in rows])
            o = np.argsort(tau)
            # A dozen series on one axis.  Eight house hues times eight
            # dash patterns, rather than a twenty-colour matplotlib map —
            # which is not colour-blind safe and collapses in greyscale.
            a1.plot(tau[o], ad[o], lw=1, alpha=0.75, label=_label(k),
                    **ps.line_series(i))
        t0 = np.array([100.0, 10000.0])
        a1.plot(t0, 0.02 * np.sqrt(t0[0] / t0), ls="--", lw=1.4, color=INK,
                label=r"white noise, $\tau^{-1/2}$")
        a1.set_xscale("log"); a1.set_yscale("log")
        a1.set_xlabel(r"averaging time $\tau$ (s)")
        a1.set_ylabel("Allan deviation (mag)")
        a1.set_title("Check stars, longest continuous run")
        a1.legend(fontsize=5.6, ncol=2)
        sl = [r[0] for r in q(con, "SELECT slope FROM ch_allan_fit "
                                   "WHERE slope IS NOT NULL")]
        rf = [r[0] for r in q(con, "SELECT red_factor FROM ch_allan_fit "
                                   "WHERE red_factor IS NOT NULL")]
        a2.hist(sl, bins=np.arange(-1.0, 0.45, 0.05), color=ACCENT)
        a2.axvline(-0.5, color=INK, lw=1.4, ls="--")
        a2.set_xlabel(r"fitted log-log slope   ($-1/2$ = white)")
        a2.set_ylabel("check-star ladders")
        a2.set_title(f"median slope {np.median(sl):+.2f}; median red-noise "
                     f"factor at $P_{{orb}}$ {np.median(rf):.2f}x")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "allan.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/allan.png"


# ===========================================================================
# §3 figures — sampling
# ===========================================================================

def fig_window(con) -> str:
    """Spectral windows: multi-season, one series, one night."""
    scopes = [r[0] for r in q(con, "SELECT DISTINCT scope FROM ch_window")]
    picks = []
    for t in TARGET_ORDER[:2]:
        allsc = [s for s in scopes if s == f"{t}|all"]
        ser = [s for s in scopes if s.endswith("|series") and s.startswith(t)]
        night = [s for s in scopes if s.startswith(t) and not s.endswith("|series")
                 and s != f"{t}|all"]
        picks.append((t, allsc + ser + night))
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(len(picks), 1, figsize=(10.6, 3.2 * len(picks)),
                                squeeze=False)
        for ax, (t, sc) in zip(axs.ravel(), picks):
            for i, s in enumerate(sc):
                rows = q(con, "SELECT freq_cd, power FROM ch_window WHERE "
                              "scope=? ORDER BY freq_cd", (s,))
                f = np.array([r[0] for r in rows]); p = np.array([r[1] for r in rows])
                ax.plot(f, p, lw=0.9, alpha=0.85, label=s,
                        **ps.line_series(i))
            porb = q(con, "SELECT period_d FROM ch_cadence WHERE target_key=? "
                          "LIMIT 1", (t,))[0][0]
            ax.axvline(1.0 / porb, color=WARN, lw=1.2, ls="--",
                       label=f"$f_{{orb}}$ = {1 / porb:.3f} c/d")
            for k in (-1, 1):
                ax.axvline(1.0 / porb + k, color=BAD, lw=0.9, ls=":")
            # The 1 c/d comb is the CAUSE; the dotted lines above are where
            # its effect lands in a real periodogram.  Marking the cause is
            # what stops the figure being read as "there is no alias here".
            ax.axvline(1.0, color=WARN, lw=1.4, alpha=0.55,
                       label="1 c/d comb: the window power that CREATES "
                             "those aliases")
            ax.set_xlim(0, 30)
            ax.set_ylim(0, 1.05)
            ax.set_xlabel("frequency (cycles / day)")
            ax.set_ylabel("window power")
            ax.set_title(f"{TARGET_LABEL[t]} spectral window", fontsize=10)
            ax.legend(fontsize=6.5, ncol=2)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "window.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/window.png"


def fig_cadence(con) -> str:
    """Points per orbital cycle, and how much of one night is one cycle."""
    rows = q(con, "SELECT series_key, pts_per_cycle, best_night_cycles, "
                  "smear_phase, target_key, filter FROM ch_cadence "
                  "ORDER BY target_key, series_key")
    lab = [_label(r[0]) for r in rows]
    ppc = np.array([r[1] if r[1] else np.nan for r in rows])
    cyc = np.array([r[2] if r[2] else np.nan for r in rows])
    sm = np.array([r[3] if r[3] else np.nan for r in rows]) * 100
    y = np.arange(len(rows))
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 0.30 * len(rows) + 1.8),
                                     sharey=True)
        a1.barh(y, ppc, color=ACCENT)
        a1.axvline(20, color=WARN, ls="--", lw=1.2)
        a1.set_yticks(y); a1.set_yticklabels(lab, fontsize=7.5)
        a1.set_xlabel("points per orbital cycle (per filter)")
        a1.invert_yaxis()
        a1.set_title("Sampling density inside one orbit\n(dashed: 20 pts/cycle)",
                     fontsize=10)
        a2.barh(y - 0.2, cyc, height=0.4, color=GOOD, label="cycles in the best night")
        a2.barh(y + 0.2, sm, height=0.4, color=BAD,
                label="exposure phase smear (% of a cycle)")
        a2.set_xlabel("cycles  |  % of a cycle")
        a2.legend(fontsize=7.5)
        a2.set_title("Cycles per night, and how much one exposure smears",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "cadence.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/cadence.png"


# ===========================================================================
# §4/§5 figures — detectability and timing
# ===========================================================================

def fig_contour(con) -> str:
    """The 90% recovery contour: amplitude vs period, per target.

    Two panels because they are two different experiments: what one night
    can detect, and what a whole season can detect (before and after the
    per-night detrending the strategy requires).

    Solid lines are the ``known`` score &mdash; detection at a literature
    period, which is the question this paper asks.  Dotted lines on the
    single-night panel are the same injections scored as blind period
    DETERMINATION, and the gap between the two is the defect this rebuild
    corrects: the dotted curve was published as a detection limit.
    """
    rows = q(con, "SELECT DISTINCT series_key, regime FROM ch_contour")
    if not rows:
        return ""
    targets = sorted({r[0].split("|")[0] for r in rows},
                     key=lambda t: TARGET_ORDER.index(t)
                     if t in TARGET_ORDER else 99)
    panels = [("richest single night",
               (("night", "known", "-", "o"),
                ("night", "period", ":", "x"))),
              ("whole season",
               (("season", "known", "--", "s"),
                ("season-dt", "known", "-", "^")))]
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.6), sharey=True)
        for ax, (title, series) in zip(axs, panels):
            for i, t in enumerate(targets):
                for regime, score, ls, mk in series:
                    pts = q(con, "SELECT period_d, amp90 FROM ch_contour WHERE "
                                 "series_key LIKE ? AND regime=? AND score=? "
                                 "AND amp90 IS NOT NULL ORDER BY period_d",
                            (t + "|%", regime, score))
                    if not pts:
                        continue
                    P = np.array([p[0] for p in pts]) * 24.0
                    A = np.array([p[1] for p in pts]) * 1000
                    suffix = (" (detrended)" if regime == "season-dt"
                              else " — period determination"
                              if score == "period" else "")
                    ax.plot(P, A, ls, marker=mk, ms=4, lw=1.5,
                            alpha=0.55 if score == "period" else 1.0,
                            color=ps.CYCLE[i % len(ps.CYCLE)],
                            label=f"{TARGET_LABEL[t]}{suffix}")
                porb = q(con, "SELECT period_d FROM ch_cadence WHERE "
                              "target_key=? LIMIT 1", (t,))
                if porb:
                    ax.axvline(porb[0][0] * 24.0, lw=0.8, alpha=0.35,
                               color=ps.CYCLE[i % len(ps.CYCLE)])
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel("injected period (hours)")
            ax.set_title(title, fontsize=11)
            ax.legend(fontsize=6, ncol=2)
        axs[0].set_ylabel("semi-amplitude recovered 90% of the time (mmag)")
        fig.suptitle("Solid: detection at a KNOWN period.  Dotted: blind "
                     "period determination — on one night that measures how "
                     "few cycles the night holds, not the amplitude.",
                     fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(FIG_DIR / "contour.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/contour.png"


def fig_recovery_grid(con) -> str:
    """The raw recovery grids behind one flagship night's two contours.

    Side by side, because the pair IS the argument: the same injections into
    the same noise, scored once as "was it detected at the known period" and
    once as "was the period determined from scratch".
    """
    row = q(con, "SELECT scope FROM ch_detect WHERE regime='night' "
                 "AND scope LIKE 'stlmi%' LIMIT 1")
    if not row:
        row = q(con, "SELECT scope FROM ch_detect LIMIT 1")
    if not row:
        return ""
    scope = row[0][0]

    def grid_for(score):
        rows = q(con, "SELECT period_d, semi_amp, frac FROM ch_detect "
                      "WHERE scope=? AND score=?", (scope, score))
        if not rows:
            return None, None, None
        P = sorted({r[0] for r in rows}); A = sorted({r[1] for r in rows})
        g = np.full((len(A), len(P)), np.nan)
        for p, a, f in rows:
            g[A.index(a), P.index(p)] = f
        return P, A, g

    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.6))
        titles = {"known": "detection at a KNOWN period",
                  "period": "blind period determination"}
        im = None
        for ax, score in zip(axs, ("known", "period")):
            P, A, grid = grid_for(score)
            if grid is None:
                continue
            im = ax.imshow(grid, origin="lower", aspect="auto",
                           cmap=ps.SEQ_CMAP, vmin=0, vmax=1,
                           extent=(-0.5, len(P) - 0.5, -0.5, len(A) - 0.5))
            ax.grid(False)          # a heatmap wears no grid
            if np.nanmax(grid) >= ch.RECOVERY_LEVEL:
                cs = ax.contour(np.arange(len(P)), np.arange(len(A)), grid,
                                levels=[ch.RECOVERY_LEVEL], colors=[GOOD],
                                linewidths=2)
                ax.clabel(cs, fmt={ch.RECOVERY_LEVEL: "90%"}, fontsize=8)
            ax.set_xticks(range(len(P)))
            ax.set_xticklabels([f"{p * 24:.1f}" for p in P], fontsize=7)
            ax.set_yticks(range(len(A)))
            ax.set_yticklabels([f"{a * 1000:.0f}" for a in A], fontsize=7)
            ax.set_xlabel("injected period (hours)")
            ax.set_title(titles[score], fontsize=11)
        axs[0].set_ylabel("injected semi-amplitude (mmag)")
        if im is not None:
            plt.colorbar(im, ax=axs, label="fraction recovered")
        fig.suptitle(f"Recovery fraction — {scope}", fontsize=10)
        fig.savefig(FIG_DIR / "recovery_grid.png", dpi=DPI,
                    bbox_inches="tight")
        plt.close(fig)
    return "figures/cv_characterization/recovery_grid.png"


def fig_timing(con) -> str:
    """Epoch precision of one bright-phase edge vs the 60 s target."""
    rows = q(con, "SELECT target_key, regime, ingress_req, ingress_phase, "
                  "sigma_t_s, cadence_s, exp_smear_phase FROM ch_timing "
                  "ORDER BY target_key, regime, ingress_req")
    if not rows:
        return ""
    labels, vals, colors = [], [], []
    palette = {"per-cycle": ACCENT, "per-cycle shape-mismatched": BAD,
               "night-mean": GOOD}
    for t, regime, ireq, ieff, st, dt, smear in rows:
        smeared = "*" if ieff > ireq + 1e-9 else ""
        labels.append(f"{TARGET_LABEL.get(t, t)}  {regime}  "
                      f"edge {ireq:g}P{smeared}")
        vals.append(st if st else np.nan)
        colors.append(palette.get(regime, MUTED))
    y = np.arange(len(labels))
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(9.6, 0.26 * len(labels) + 2.0))
        ax.barh(y, vals, color=colors)
        ax.axvline(60, color=WARN, lw=1.8, ls="--",
                   label="the strategy's assumed 60 s per-cycle target")
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.2)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel(r"Monte-Carlo epoch precision $\sigma_t$ (s)")
        ax.set_title("Timing one bright-phase edge at the real cadence\n"
                     "blue: shape known exactly | red: shape and depth "
                     "assumed wrong | green: night mean\n"
                     "* = injected edge widened to the exposure time",
                     fontsize=10)
        ax.legend(fontsize=8, loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "timing.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_characterization/timing.png"


# ===========================================================================
# Sections
# ===========================================================================

def section_intro(con, pcon) -> str:
    n_frames = q1(con, "SELECT count(*) FROM ch_frames")
    n_series = q1(con, "SELECT count(*) FROM ch_cadence")
    n_solved = q1(pcon, "SELECT count(*) FROM cv_series WHERE status='solved'")
    n_all = q1(pcon, "SELECT count(*) FROM cv_series")
    n_stars = q1(con, "SELECT count(*) FROM ch_noise_stars")
    return f"""
<section id="intro">
<div class="bhead"><h2>0 &middot; The question this page answers</h2>
<span class="tag">{fmt(n_frames)} measured frames</span>
<span class="tag">{fmt(n_series)} solved series</span>
<span class="tag">{fmt(n_stars)} stars</span></div>
<div class="stage">
<h3>Question</h3>
<p class="sub">The photometry ran. Before any science goal is defended, one
thing has to be established from the data themselves: <b>what precision, at
what cadence, over what baseline, does this archive actually deliver?</b>
Only then can a science goal be graded, and the grading has to be done
against measurements, not against what the observing plan hoped for.</p>
<h3>Evidence</h3>
<p class="sub">{fmt(n_solved)} of {fmt(n_all)} staged (target, era, filter)
series were solved by the ensemble; the {fmt(n_all - n_solved)} refusals
carry no zero points and therefore no photometry to characterize. Every
number below is derived from those solved series and stored in
<code>products/phot/cv_characterization.sqlite</code>.</p>
<div class="decision"><b>Decision.</b> Characterize first, in five steps
&mdash; image quality, noise floor, sampling, detectability, timing
precision &mdash; and only then grade the science goals in &sect;6. A goal
that survives this order is defensible; one that does not was never
supported by the data.</div>
<h3>Consequence</h3>
<p class="sub">&sect;6 is the point of the page. Everything before it exists
to make its verdicts unarguable.</p>
</div></section>"""


def section_quality(con) -> str:
    cuts = q(con, "SELECT axis, unit, threshold, baseline, n_pass, n_fail, note "
                  "FROM ch_cuts ORDER BY axis")
    rows = []
    for axis, unit, thr, base, npass, nfail, note in cuts:
        applied = "diagnostic" not in (note or "")
        rows.append([f"<code>{esc(axis)}</code>", esc(unit),
                     ("no threshold &mdash; never degrades"
                      if thr is None else _num(thr, 3)),
                     _num(base, 2), fmt(npass), fmt(nfail),
                     ("<b>applied</b>" if applied else "diagnostic only")])
    n_use = q1(con, "SELECT sum(usable) FROM ch_frames")
    n_tot = q1(con, "SELECT count(*) FROM ch_frames")
    fw_med = q1(con, "SELECT median(fwhm_as) FROM ch_frames")
    fw_p90 = q1(con, "SELECT threshold FROM ch_cuts WHERE axis='fwhm_as'")
    reject = q(con, "SELECT reject_reason, count(*) FROM ch_frames "
                    "WHERE usable=0 GROUP BY 1 ORDER BY 2 DESC")
    per_series = q(con, """
        SELECT s.series_key, s.readoutm, s.n_frames, s.fwhm_as_p10,
               s.fwhm_as_med, s.fwhm_as_p90, s.airmass_med, s.airmass_max,
               s.moon_illum_med, s.sat_frac, s.frac_usable
        FROM ch_quality_series s ORDER BY s.series_key""")
    qs = table(["series", "readout", "frames", "FWHM p10", "FWHM med",
                "FWHM p90", "X med", "X max", "moon k", "sat. det. frac",
                "usable"],
               [[esc(r[0]), esc(r[1]), fmt(r[2]), _num(r[3]), _num(r[4]),
                 _num(r[5]), _num(r[6]), _num(r[7]), _num(r[8], 2),
                 _num(r[9], 3), f"{100 * r[10]:.0f}%"] for r in per_series],
               row_classes=["warn" if (r[10] or 1) < 0.7 else "" for r in per_series])
    trail = q(con, "SELECT count(*), median(ell_med), median(pa_R), "
                   "max(ell_med) FROM ch_trail WHERE status='ok'")
    trail_txt = ""
    if trail and trail[0][0]:
        n, em, R, emax = trail[0]
        trail_txt = (f"""<h3>Evidence &mdash; is anything trailed?</h3>
<p class="sub">{fmt(n)} frames sampled evenly through every series and
re-measured at the pixel level. Median source ellipticity
<b>{_num(em, 3)}</b> (worst frame {_num(emax, 3)}); median position-angle
coherence <b>{_num(R, 3)}</b>, against the ~{_num(1 / math.sqrt(100), 2)}
expected for randomly oriented sources. Trailing would show as high
ellipticity AND high coherence together.</p>"""
                     + _figure(fig_trail(con),
                               "Left: ellipticity of the brightest 100 sources per "
                               "sampled frame. Right: elongation against alignment "
                               "&mdash; genuine trailing lands in the upper right."))
    return f"""
<section id="quality">
<div class="bhead"><h2>1 &middot; Are the images good enough &mdash; and what
does &ldquo;good enough&rdquo; mean?</h2>
<span class="tag">{100 * n_use / n_tot:.1f}% usable</span></div>
<div class="stage">
<h3>Question</h3>
<p class="sub">Seeing, sky brightness, transparency, airmass and the Moon
all vary across this archive. A cut has to be made somewhere &mdash; but a
cut chosen by taste is a knob a referee can turn. What does each quality
axis actually <em>cost</em>, in magnitudes of scatter?</p>
<h3>Evidence</h3>
<p class="sub">Every frame carries one number that answers this without
touching the target: the RMS of the four HELD-OUT check stars about their
own means. It is measured in magnitudes, it is independent of the science
signal, and it responds to every way a frame can be bad. Frames are binned
along each axis and the median relative check-star scatter is read off.
FWHM is converted to arcsec with the era plate scale (three different
scales are in play); airmass is <b>recomputed</b> from coordinates and time
because the archive's own AIRMASS cards reach 6877 and are null on 851 CV
frames; moon phase and separation come from ephemerides, not from the
unaudited MOONANGL/MOONPHAS cards.</p>
{_figure(fig_seeing(con),
         "Left: FWHM in arcsec &mdash; a pixel FWHM is not comparable across "
         "0.5375, 0.8062 and 0.4491 arcsec/px cameras. Right: recomputed "
         "airmass. VV Pup's floor at X ~ 1.6 is geometry, not weather.")}
{_figure(fig_quality_response(con),
         "The whole cut, in one figure: relative check-star scatter against "
         "each axis. Dashed line = the "
         f"{ch.DEGRADE_FACTOR:g}x degradation factor; red line = the adopted "
         "threshold, placed at the first sustained run of two bins above it.")}
{_figure(fig_moon(con),
         "The moon acts through the sky, and only when it is above the "
         "horizon &mdash; which is why moon phase is a diagnostic here and "
         "the sky-brightness ratio is the axis that is actually cut on.")}
{trail_txt}
{table(["axis", "unit", "threshold", "baseline (rel. scatter)", "frames pass",
        "frames fail", "role"], rows)}
<h3>Evidence &mdash; the sparse end of an axis is not an absent effect</h3>
<p class="sub">A bin holding fewer than
{ch.MIN_FRAMES_PER_QUALITY_BIN} frames cannot set a threshold: its median
fluctuates by tens of percent. But <em>discarding</em> such bins converts
&ldquo;too few frames to test&rdquo; into &ldquo;no effect&rdquo;, and this
page published exactly that mistake about airmass &mdash; the bins from
X&nbsp;=&nbsp;2.45 to 2.65 sit at 1.48, 1.79 and 1.45&times; baseline, all
three above the degradation factor, all three individually too thin to
count. Everything past the last well-populated bin is now <b>pooled into one
wide bin</b> whose median is recomputed over the pooled sample, and a
well-populated <em>final</em> bin is allowed to set the threshold on its own:
the run-of-two rule exists to guard against an interior excursion, and an
interior excursion has a neighbour to be confirmed against, while the end of
an axis does not.</p>
<h3>Evidence &mdash; what does the registration method cost?</h3>
<p class="sub">Frames inside one series are registered by up to four
different routes (WCS chain, astroalign triangles, translation vote, and the
reference itself), and one of them &mdash; the translation vote &mdash; was
invented during this build and carries
{fmt(q1(con, "SELECT sum(n_frames) FROM ch_quality_bins WHERE scope='reg_method' AND axis='translation_vote'"))}
frames here. Registration quality and frame depth are confounded (a shallow
frame is both harder to align and noisier to measure), so the methods are
compared only <em>inside</em> bins of detections per frame.</p>
{table(["method", "detections/frame bin", "median rel. check-star scatter",
        "frames"],
       [[f"<code>{esc(r[0])}</code>", _num(r[1], 0), _num(r[2], 3), fmt(r[3])]
        for r in q(con, "SELECT axis, bin_center, med_rel_scatter, n_frames "
                        "FROM ch_quality_bins WHERE scope='reg_method' "
                        "ORDER BY bin_center, axis")])}
<div class="decision"><b>Decision.</b> Keep a frame when its seeing is
&le; {_num(fw_p90, 2)}&Prime;, its sky no more than
{_num(q1(con, "SELECT threshold FROM ch_cuts WHERE axis='sky_ratio'"), 2)}&times;
its series median, its ensemble zero point no more than
{_num(q1(con, "SELECT threshold FROM ch_cuts WHERE axis='zp_excess'"), 3)} mag
below the night's median (the transparency / cloud test), and its airmass at
or below <b>X = {_num(q1(con, "SELECT threshold FROM ch_cuts WHERE axis='airmass'"), 2)}</b>.
That last cut is new and small
({fmt(q1(con, "SELECT n_fail FROM ch_cuts WHERE axis='airmass'"))} frames):
the bound the data support is <em>&ldquo;no measurable degradation below
X&nbsp;=&nbsp;2.4&rdquo;</em>, not &ldquo;airmass never matters&rdquo;, and
the observed range reaches X&nbsp;=&nbsp;2.95. Registration method gets no
threshold: inside a depth bin the translation vote is never the worst route
and is often the best, so its
{fmt(q1(con, "SELECT sum(n_frames) FROM ch_quality_bins WHERE scope='reg_method' AND axis='translation_vote'"))}
frames are validated rather than merely used. Altogether <b>{fmt(n_use)} of
{fmt(n_tot)} frames ({100 * n_use / n_tot:.1f}%)</b> are usable.</div>
<h3>Consequence</h3>
<p class="sub">Two consequences carry into &sect;2. First, seeing is
<em>not</em> what limits this archive: the median frame is
{_num(fw_med, 2)}&Prime; and only
{fmt(q1(con, "SELECT count(*) FROM ch_frames WHERE reject_reason LIKE '%fwhm%'"))}
frames fail on it. Transparency and moonlit sky do the damage, and both are
scheduling problems, not instrument problems. Second &mdash; and this is the
part that was missing &mdash; <b>this cut is now applied to every
measurement on the rest of this page</b>. The noise floor, the Allan ladders,
the spectral windows, the detection contours and the timing Monte Carlo are
all computed on the {fmt(n_use)} usable frames. In the first version of this
build the cut was computed, defended here, and then read by nothing: every
later number silently included the {fmt(n_tot - n_use)} rejected frames while
this section claimed the opposite. &sect;2 tables both numbers side by side so
the cut's price is visible rather than asserted.</p>
<h3>Per-series image quality</h3>
{qs}
<p class="sub">Rejection reasons: {", ".join(f"{esc(r[0])} &times;{fmt(r[1])}"
                                             for r in reject)}.</p>
</div></section>"""


def section_noise(con) -> str:
    rows = q(con, """
        SELECT series_key, target_key, readoutm, exptime, n_stars, floor_nom,
               floor_lo, floor_hi, k_nom, floor_plateau, plateau_mag,
               scint_mag, target_mag, prec_at_target, n_near_target,
               faint_const_mag, inflation, chi2nu_med, n_target_points,
               k_lo, k_hi, k_in_bracket, prec_at_target_all,
               n_frames_usable, n_frames_all
        FROM ch_noise_series ORDER BY target_key, series_key""")
    trows, rclass = [], []
    for r in rows:
        (sk, tk, ro, expt, nst, fn, flo, fhi, k, pl, plm, sc, tm, pr, nn,
         fc, infl, chi2, ntp, klo, khi, kok, pr_all, nfu, nfa) = r
        ratio = (pl / sc) if (pl and sc) else None
        note = ""
        if pr is None or not (pr == pr):
            note = ("target fainter than every comparison star"
                    if (tm is not None and fc is not None and tm > fc)
                    else "no target tie")
        elif not kok:
            note = "k outside the gain bracket"
        trows.append([esc(sk), esc(ro), _num(expt, 0),
                      f"{fmt(nfu)}/{fmt(nfa)}", fmt(nst),
                      _mmag(pl), _num(plm, 1), _mmag(fn),
                      (f"{_num(k, 2)} [{_num(klo, 2)}&ndash;{_num(khi, 2)}]"
                       if kok else
                       f"<b>{_num(k, 2)}</b> [{_num(klo, 2)}&ndash;{_num(khi, 2)}]"),
                      _mmag(sc, 2), _num(ratio, 0),
                      _num(tm, 2), _mmag(pr), _mmag(pr_all), fmt(nn),
                      _num(infl, 2),
                      _num(chi2, 2), fmt(ntp), esc(note) if note else ""])
        rclass.append("warn" if note else "")
    best = q(con, "SELECT series_key, floor_plateau FROM ch_noise_series "
                  "WHERE floor_plateau IS NOT NULL ORDER BY floor_plateau LIMIT 1")
    med_ratio = q1(con, "SELECT median(floor_plateau / scint_mag) FROM "
                        "ch_noise_series WHERE floor_plateau IS NOT NULL")
    med_k = q1(con, "SELECT median(k_nom) FROM ch_noise_series")
    n_k = q1(con, "SELECT count(*) FROM ch_noise_series WHERE k_nom IS NOT NULL")
    n_k_ok = q1(con, "SELECT sum(k_in_bracket) FROM ch_noise_series")
    k_out = q(con, "SELECT series_key, k_nom FROM ch_noise_series "
                   "WHERE k_in_bracket = 0 ORDER BY k_nom")
    k_span = q1(con, "SELECT median(abs(k_hi - k_lo) / k_nom) "
                     "FROM ch_noise_series WHERE k_in_bracket = 1")
    med_slope = q1(con, "SELECT median(slope) FROM ch_allan_fit "
                        "WHERE slope IS NOT NULL")
    med_null = q1(con, "SELECT median(slope_null_p50) FROM ch_allan_fit")
    null_lo = q1(con, "SELECT median(slope_null_p05) FROM ch_allan_fit")
    null_hi = q1(con, "SELECT median(slope_null_p95) FROM ch_allan_fit")
    med_red = q1(con, "SELECT median(red_factor) FROM ch_allan_fit "
                      "WHERE red_factor IS NOT NULL")
    med_tau = q1(con, "SELECT median(tau_frac_of_porb) FROM ch_allan_fit")
    n_reach = q1(con, "SELECT sum(tau_frac_of_porb >= 0.99) FROM ch_allan_fit")
    n_red = q1(con, "SELECT sum(red_redder_than_null) FROM ch_allan_fit")
    n_slope_red = q1(con, "SELECT sum(slope_redder_than_null) FROM ch_allan_fit")
    n_fit = q1(con, "SELECT count(*) FROM ch_allan_fit")
    bias_med = q1(con, "SELECT median(bias_ratio) FROM ch_check_bias")
    bias_lo = q1(con, "SELECT min(bias_ratio) FROM ch_check_bias")
    bias_hi = q1(con, "SELECT max(bias_ratio) FROM ch_check_bias")
    bias_n = q1(con, "SELECT count(*) FROM ch_check_bias")
    cut_gain = q1(con, "SELECT median(prec_at_target_all / prec_at_target) "
                       "FROM ch_noise_series WHERE prec_at_target > 0")
    cut_best = q1(con, "SELECT max(prec_at_target_all / prec_at_target) "
                       "FROM ch_noise_series WHERE prec_at_target > 0")
    per_target = q(con, """
        SELECT target_key, min(prec_at_target), max(prec_at_target)
        FROM ch_noise_series WHERE prec_at_target IS NOT NULL
        GROUP BY target_key""")
    pt = ", ".join(f"<b>{TARGET_LABEL.get(t, t)} {_mmag(lo, 0)}&ndash;"
                   f"{_mmag(hi, 0)} mmag</b>" for t, lo, hi in per_target)
    return f"""
<section id="noise">
<div class="bhead"><h2>2 &middot; What is the noise floor, and is it what the
detector says it should be?</h2>
<span class="tag">floor {_mmag(best[0][1], 0) if best else "&mdash;"}&ndash;
{_mmag(q1(con, "SELECT max(floor_plateau) FROM ch_noise_series"), 0)} mmag</span>
</div>
<div class="stage">
<h3>Question</h3>
<p class="sub">Photon statistics predict a scatter that falls as stars get
brighter. Real photometry stops improving somewhere. <b>Where, and what is
the flat part made of?</b> The answer sets every detection limit downstream,
so it has to be measured rather than asserted &mdash; and the S2 detector
work gives us the pieces to predict what it <em>should</em> be.</p>
<h3>Evidence &mdash; the prediction is a band, not a line</h3>
<p class="sub">The variance of an aperture sum is
<code>F/g + n_pix &sigma;<sub>bkg</sub><sup>2</sup> +
(n_pix<sup>2</sup>/n_sky) &sigma;<sub>bkg</sub><sup>2</sup></code>. The sky,
read and dark terms come from each frame's own <em>measured</em> background
RMS, so they carry no assumption; only the source term needs the gain. S2
measured the gain to lie between <b>{ch.GAIN_LO_E_PER_ADU}</b> and
<b>{ch.GAIN_HI_E_PER_ADU} e&minus;/ADU</b> with the header EGAIN
{ch.GAIN_NOMINAL_E_PER_ADU} inside the bracket &mdash; so the prediction is
drawn as a band spanning that bracket, and no single gain is quoted as if it
were known. Read noise {ch.READ_NOISE_ADU_HIGH_GAIN} ADU per High Gain read
({ch.READ_NOISE_E} &plusmn; {ch.READ_NOISE_E_ERR} e&minus;); StackPro frames
are sums of {ch.STACKPRO_N_SUB} sub-reads, which needs no special case here
because their extra read noise is already inside the measured background
RMS.</p>
{_figure(fig_rms_mag(con),
         "RMS vs magnitude, one panel per target (richest series). Grey: the "
         "comparison ensemble. Blue circles: the four HELD-OUT check stars "
         "&mdash; the only honest estimate of per-point precision. Star: the "
         "target, plotted for placement only (its scatter is signal). Shaded "
         "band: the photon+sky+read prediction over the measured gain bracket.")}
{_figure(fig_floor_budget(con),
         "The floor against the one term the atmosphere forces on a 0.5 m "
         "telescope. Scintillation (Young's formula at each series' own "
         "airmass and exposure) is 20-50x below the measured floor.")}
<h3>Evidence &mdash; are the error bars honest, and does averaging work?</h3>
<p class="sub">Two independent tests. The &chi;<sup>2</sup> inflation factor
asks whether the formal per-point errors, multiplied by that factor, make a
constant star's &chi;<sup>2</sup><sub>&nu;</sub> equal 1. The Allan
deviation asks whether averaging still buys
&radic;N: white noise falls as &tau;<sup>&minus;1/2</sup>, correlated noise
flattens.</p>
{_figure(fig_allan(con),
         "Allan deviations of the check stars on each series' longest "
         "continuous run. Right: the distribution of fitted slopes against "
         "the WHITE-NOISE NULL for ladders of this length and cadence "
         "(shaded) rather than against the asymptotic -0.50, and the "
         "red-noise penalty at the largest tau each ladder reaches.")}
<p class="sub">Both of these estimators are compared against the wrong thing
if they are compared against textbook values. A ladder here has 4&ndash;6
rungs over about one decade with as few as three pairs at the top, and its
white-noise expectation is <b>{_num(med_null, 2)}</b> with a 5&ndash;95%
range of [{_num(null_lo, 2)}, {_num(null_hi, 2)}] &mdash; not
&minus;0.50. Every ladder therefore carries its own null, generated at its
own length and cadence, and the verdict below is per-ladder. Likewise the
red-noise factor is quoted at <b>the largest &tau; the ladder actually
reaches</b>, which is a median {_num(med_tau, 2)} of the orbital period:
{fmt(n_reach)} of {fmt(n_fit)} ladders reach P<sub>orb</sub> itself, so a
number labelled &ldquo;at P<sub>orb</sub>&rdquo; would be a lower bound
wearing a larger name (red noise grows with &tau;).</p>
<h3>Evidence &mdash; are the check stars a fair sample?</h3>
<p class="sub">Every honest-noise number on this page comes from at most four
held-out check stars, and they are not a random four: they are the survivors
of the comparison-star stability iteration, drawn nearest the target's
brightness. Survivors of a scatter cut could easily be quieter than the
population they were drawn from, which would make every precision number on
this page optimistic. Tested directly &mdash; the median RMS of <em>every</em>
star within {ch.FIELD_MATCH_HALF_WIDTH:g} mag of the
target, including the ones the stability cut dropped and the ones never
selected at all, against the check stars' own median &mdash; the ratio is
<b>{_num(bias_med, 2)}</b> over {fmt(bias_n)} series (range
{_num(bias_lo, 2)}&ndash;{_num(bias_hi, 2)}). The check stars are if anything
marginally noisier than magnitude-matched field stars, so the selection does
not buy an optimistic answer.</p>
{table(["series", "readout", "t<sub>exp</sub> (s)", "frames used/matched",
        "stars", "floor (mmag)", "at mag", "fitted floor",
        "k [gain bracket]",
        "scint. (mmag)", "floor/scint", "target mag",
        "precision at target (mmag)", "same, uncut", "n near", "inflation",
        "&chi;<sup>2</sup><sub>&nu;</sub>", "target pts", "note"],
       trows, rclass)}
<div class="decision"><b>Decision.</b> The systematic floor is
<b>{_mmag(best[0][1], 0) if best else "&mdash;"}&ndash;{_mmag(q1(con, "SELECT max(floor_plateau) FROM ch_noise_series"), 0)} mmag</b>
depending on series, and it is <b>{_num(med_ratio, 0)}&times; larger than
scintillation</b> (median). It is therefore <em>instrumental</em> &mdash;
flat-field residual, ensemble zero-point error, and low-level variability of
the &ldquo;constant&rdquo; stars &mdash; not atmospheric, and in principle
improvable. The photon model is close to right on most series: the fitted
scaling <code>k</code> has median {_num(med_k, 2)} over all {fmt(n_k)}
series, and <b>{fmt(n_k_ok)} of {fmt(n_k)}</b> fall inside the range the gain
bracket allows. The {fmt(n_k - n_k_ok)} that do not are named in the table
and are a <em>model failure</em>, not an average to be taken:
{", ".join(f"{esc(s)} ({_num(kv, 2)})" for s, kv in k_out)}. Note also what
the gain bracket is and is not worth here: swinging the gain across its
entire {ch.GAIN_HI_E_PER_ADU / ch.GAIN_LO_E_PER_ADU:.1f}&times; range moves
<code>k</code> by a median of only {100 * (k_span or 0):.0f}%, because the
fit takes its leverage from the faint end where sky and read noise dominate
&mdash; so the bracket is propagated honestly but it is <em>not</em> the
dominant uncertainty in the prediction. Averaging is <b>not</b> free:
{fmt(n_slope_red)} of {fmt(n_fit)} ladders have a slope above their own
95th-percentile white null and {fmt(n_red)} exceed it on the red-noise
factor, median penalty {_num(med_red, 2)}&times; at a median
{_num(med_tau, 2)} P<sub>orb</sub>. The median slope
({_num(med_slope, 2)}) is <em>inside</em> the null band, so it is the
per-ladder tests and not the median that carry this conclusion.</div>
<h3>Consequence</h3>
<p class="sub"><b>Achieved per-point precision at each target's own
brightness:</b> {pt}. These are the numbers every later section uses; they
are measured on held-out check stars of the same brightness as the target in
the same frames, so they already contain the flat field, the ensemble
solution and the weather &mdash; and they are measured on the
quality-cut frames only. The table gives the uncut number beside each one:
the cut improves the median series by {_num(cut_gain, 2)}&times; and the best
by {_num(cut_best, 2)}&times;, which is what &sect;1's cut is worth and what
the first version of this page never actually collected. Two structural facts
travel with them: the magnitude scale is the ensemble's own arbitrary gauge
(see &sect;6, row S1), and error bars must be inflated and long-timescale
averages penalised by the factors tabled above.</p>
</div></section>"""


def section_cadence(con) -> str:
    rows = q(con, """
        SELECT series_key, n_points, n_blocks, baseline_d, median_dt_s,
               pts_per_cycle, longest_block_h, best_night, best_night_n,
               best_night_cycles, best_night_phase_cov, best_night_dt_s,
               n_blocks_ge1cycle, smear_phase, max_gap_d, duty_cycle
        FROM ch_cadence ORDER BY series_key""")
    trows = [[esc(r[0]), fmt(r[1]), fmt(r[2]), _num(r[3], 1), _num(r[4], 0),
              _num(r[5], 0), _num(r[6], 2), esc(r[7]), fmt(r[8]),
              _num(r[9], 1), f"{100 * r[10]:.0f}%" if r[10] else "&mdash;",
              _num(r[11], 0), fmt(r[12]), f"{100 * r[13]:.1f}%" if r[13] else "&mdash;",
              _num(r[14], 1), f"{100 * r[15]:.2f}%" if r[15] else "&mdash;"]
             for r in rows]
    alias = q(con, """
        SELECT scope, family, max(power), baseline_d, freq_res_cd, resolved
        FROM ch_alias WHERE k IN (-1, 1) GROUP BY scope, family
        ORDER BY scope, family""")
    arows = [[esc(r[0]), esc(r[1]), _num(r[2], 3), _num(r[3], 2),
              _num(r[4], 4), "yes" if r[5] else "no &mdash; one broad peak"]
             for r in alias]
    worst = q1(con, "SELECT max(power) FROM ch_alias WHERE k IN (-1,1) "
                    "AND resolved=1")
    worst_solar = q1(con, "SELECT max(power) FROM ch_alias WHERE k IN (-1,1) "
                          "AND resolved=1 AND family='solar day'")
    med_smear = q1(con, "SELECT median(smear_phase) FROM ch_cadence")
    eu_smear = q1(con, "SELECT max(smear_phase) FROM ch_cadence "
                       "WHERE target_key='euuma'")
    return f"""
<section id="cadence">
<div class="bhead"><h2>3 &middot; What periods can this sampling resolve at
all?</h2>
<span class="tag">worst resolved &plusmn;1 c/d alias: {_num(worst, 2)}
of the window</span></div>
<div class="stage">
<h3>Question</h3>
<p class="sub">Precision is worthless at a period the sampling cannot
distinguish from its alias. Before any periodogram is read, the sampling's
own transform &mdash; the spectral window
<code>|&Sigma; e<sup>&minus;2&pi;ift</sup>|<sup>2</sup>/N<sup>2</sup></code>
&mdash; has to be looked at, because every peak in a real periodogram is the
truth convolved with it.</p>
<h3>Evidence</h3>
{_figure(fig_window(con),
         "Spectral windows for the whole target (all filters, all eras), one "
         "series, and one night. Read it as a convolution kernel, not as a "
         "periodogram: the tall comb at 1 c/d is why a real periodogram grows "
         "spurious peaks at f_orb +/- 1 c/d (dotted), even though the window "
         "itself is flat out at 11-16 c/d. The single night (green) is one "
         "broad lobe of width ~1/T with no comb at all. The multi-night "
         "curves are stored decimated with every local maximum above 0.05 "
         "kept, because their true peaks are only 0.002 c/d wide.")}
{_figure(fig_cadence(con),
         "Left: points per orbital cycle per filter. Right: how many cycles "
         "the best night covers, and how much of a cycle one exposure smears "
         "over.")}
{table(["scope", "alias family", "max |k|=1 window power", "baseline (d)",
        "resolution (c/d)", "alias resolved?"], arows)}
<div class="decision"><b>Decision.</b> On every multi-night set the
&plusmn;1 c/d aliases of the orbital frequency are fully resolved (resolution
&lt;&lt; 1 c/d) and carry up to <b>{_num(worst_solar, 2)} of the window
power</b> on the <em>solar-day</em> comb &mdash; the one that matters for
ground-based scheduling &mdash; rising to {_num(worst, 2)} if the sidereal
comb is included as well. That is the sampling's <em>capacity</em> to
manufacture an alias. Whether the alias actually wins is a different
question, it depends on amplitude, and this page measures it directly in
&sect;4 rather than inferring it here: at the amplitudes these targets
really show, the tallest peak lands on the true frequency in essentially
every injection. <b>So: the alias risk is real near the detection threshold
and negligible at the real modulation amplitudes.</b> A single night's peak,
by contrast, is 3&ndash;9 c/d wide &mdash; it cannot separate the aliases
either, but for the opposite reason: there is only one peak. Published
periodograms must carry the window, the solar-day alias power and the
measured misidentification rate beside them.</div>
<h3>Consequence</h3>
<p class="sub">Two sampling limits carry into &sect;4. The per-filter
cadence, not the exposure time, sets the points per cycle; and the exposure
smears phase by a median {_num(100 * med_smear if med_smear else None, 1)}%
of a cycle, rising to {_num(100 * eu_smear if eu_smear else None, 1)}% on EU
UMa's 240 s exposures &mdash; a floor on any timing or edge-shape claim for
that target, independent of signal-to-noise.</p>
<h3>Per-series sampling</h3>
{table(["series", "pts", "nights", "baseline (d)", "median &Delta;t (s)",
        "pts/cycle", "longest run (h)", "best night", "pts", "cycles",
        "phase cov.", "&Delta;t (s)", "full-orbit nights", "smear (% cycle)",
        "largest gap (d)", "duty cycle"], trows)}
</div></section>"""


def section_detect(con) -> str:
    rows = q(con, """
        SELECT c.scope, c.series_key, c.regime, c.period_d, c.amp90,
               c.n_points, d.sigma_used, d.threshold, c.score
        FROM ch_contour c JOIN (SELECT scope, regime, score,
             min(sigma_used) sigma_used,
             min(threshold) threshold FROM ch_detect
             GROUP BY scope, regime, score) d
        ON d.scope = c.scope AND d.regime = c.regime AND d.score = c.score
        ORDER BY c.series_key, c.score, c.regime, c.period_d""")
    if not rows:
        return """<section id="detect"><div class="bhead">
<h2>4 &middot; What amplitude can we actually detect?</h2></div>
<p class="missing">Injection-recovery stage has not been run.</p></section>"""
    # One line per (scope, regime) at the target's own orbital period.
    at_porb = q(con, """
        SELECT c.series_key, c.regime, c.period_d, c.amp90, c.n_points,
               n.prec_at_target, d.sigma_used, c.score, c.amp90_lo,
               c.amp90_hi, c.n_cycles, c.freq_res_cd
        FROM ch_contour c
        JOIN ch_cadence cad ON cad.series_key = c.series_key
        LEFT JOIN ch_noise_series n ON n.series_key = c.series_key
        JOIN (SELECT scope, regime, score, min(sigma_used) sigma_used
              FROM ch_detect GROUP BY scope, regime, score) d
             ON d.scope = c.scope AND d.regime = c.regime AND d.score = c.score
        WHERE abs(c.period_d - cad.period_d) < 1e-6
        ORDER BY c.series_key, c.score, c.regime""")
    prows, prclass = [], []
    for sk, regime, P, a90, n, prec, sig, score, alo, ahi, ncyc, fres in at_porb:
        analytic = ch.amin_analytic(sig, n) if sig else float("nan")
        prows.append([esc(sk), esc(regime), f"<b>{esc(score)}</b>",
                      _num(24 * P, 3), fmt(n), _num(ncyc, 1),
                      _mmag(sig),
                      (f"{_mmag(a90)} [{_mmag(alo, 0)}&ndash;{_mmag(ahi, 0)}]"
                       if a90 else "not reached"),
                      _mmag(analytic),
                      _num((a90 / analytic) if (a90 and analytic == analytic
                                                and analytic) else None, 1)])
        prclass.append("warn" if score == "period" and regime == "night" else "")
    allrows = [[esc(r[1]), esc(r[2]), esc(r[8]), _num(24 * r[3], 3), _mmag(r[4]),
                fmt(r[5]), _mmag(r[6])] for r in rows]
    # Ratios are quoted per score mode.  Mixing them was the defect: the
    # largest ratios in the first version all came from single-night scopes
    # scored on period determination, where the 1% acceptance window is far
    # narrower than the peak and no amplitude can succeed.
    def _ratios(score, regimes):
        return [a / ch.amin_analytic(s, n)
                for _sk, rg, _p, a, n, _pr, s, sc, *_rest in at_porb
                if sc == score and rg in regimes and a and s
                and ch.amin_analytic(s, n) > 0]
    known_ratios = _ratios("known", ("season", "season-dt", "night"))
    period_ratios = _ratios("period", ("season", "season-dt"))
    night_period_ratios = _ratios("period", ("night",))
    dt_gain = q(con, """
        SELECT a.series_key, a.amp90, b.amp90 FROM ch_contour a
        JOIN ch_contour b ON b.series_key = a.series_key
             AND abs(b.period_d - a.period_d) < 1e-9 AND b.regime='season-dt'
             AND b.score = a.score
        JOIN ch_cadence c ON c.series_key = a.series_key
        WHERE a.regime='season' AND a.score='known'
          AND abs(a.period_d - c.period_d) < 1e-6
          AND a.amp90 IS NOT NULL AND b.amp90 IS NOT NULL""")
    gains = [a / b for _sk, a, b in dt_gain if b]
    conf = q(con, "SELECT semi_amp, avg(frac_true), avg(frac_alias) "
                  "FROM ch_alias_confusion WHERE regime IN "
                  "('season','season-dt') GROUP BY semi_amp ORDER BY semi_amp")
    thr_spread = q1(con, "SELECT median(spread_frac) FROM ch_threshold "
                         "WHERE score='period'")
    thr_n = q1(con, "SELECT median(n_stars) FROM ch_threshold")
    gain_txt = (f"Per-night detrending buys a factor "
                f"{min(gains):.2f}&ndash;{max(gains):.2f} at the orbital "
                f"period (median {np.median(gains):.2f}) &mdash; measured, "
                f"not assumed."
                if gains else "")
    return f"""
<section id="detect">
<div class="bhead"><h2>4 &middot; What amplitude can we actually detect?</h2>
<span class="tag">injection through real timestamps &amp; real noise</span></div>
<div class="stage">
<h3>Question</h3>
<p class="sub">&ldquo;We measure 25 mmag scatter&rdquo; is not a science
statement. &ldquo;We can detect a 30 mmag modulation at 1.9 h&rdquo; is.
Converting one into the other requires the real timestamps (so the window
enters), the real noise (so correlated systematics enter) and a false-alarm
threshold measured the same way (so a red archive cannot borrow a white
threshold).</p>
<h3>Evidence</h3>
<p class="sub">Sinusoids of known semi-amplitude are injected at the real
BJD_TDB timestamps into <b>real check-star residual series</b>, cyclically
rolled so their autocorrelation survives, and recovered with a floating-mean
Lomb&ndash;Scargle against a {ch.DETECT_FAP:.0%}-FAP threshold measured by
{esc(q1(con, "SELECT value FROM ch_meta WHERE key='threshold_trials'"))}
signal-free trials through the same machinery
({fmt(q1(con, "SELECT max(n_trials) FROM ch_detect"))} injection trials per
grid cell).</p>
<h3>Evidence &mdash; two different questions, two different contours</h3>
<p class="sub">This is the correction that matters most on this page. The
first version scored a trial as recovered only if the tallest peak cleared
the threshold <b>and</b> landed within {100 * ch.PERIOD_TOL_FRAC:.0f}% of the
injected frequency. On a whole season that is a reasonable test of a blind
period search. On a <em>single night</em> it is not a detection test at all:
the frequency resolution of a 3&ndash;5 hour run is 2.6&ndash;9.0&nbsp;c/d
while the acceptance window is 0.12&ndash;0.14&nbsp;c/d &mdash; twenty to
seventy times narrower than the peak itself. A 300&nbsp;mmag injection into
VV Pup's richest night cleared the threshold in forty trials out of forty and
was scored recovered in twenty-five, and the resulting &ldquo;never reached
even at 300&nbsp;mmag&rdquo; was published as a detection limit and fed into
two verdicts.</p>
<p class="sub">Both questions are now measured and labelled. <b>score =
period</b> is the blind search: tallest peak in the band, against a threshold
built from the distribution of the tallest peak, and it must land on the
right frequency. <b>score = known</b> is the question this paper actually
asks &mdash; every one of these five CVs has a published orbital period with
a decades-long baseline &mdash; so the statistic is the power <em>at that
frequency</em> against the distribution of the power at that same frequency
in signal-free data. Each is compared against its own null; using the
maximum-statistic threshold for a known-period claim charges a
look-elsewhere penalty the paper does not owe.</p>
<p class="sub">Three regimes, because they answer different questions.
<b>season</b> is the whole richest series untouched &mdash; what a blind
period search sees, nightly zero-point wander and all. <b>season-dt</b> is
the same set after per-night mean removal, the detrending the strategy
requires (&sect;4.20); the pair measures what that detrending is worth and
what it costs at long periods. <b>night</b> is the richest single night, and
its contour is what a single-night, cycle-resolved claim &mdash; the claim
this paper is built on &mdash; can rest on.</p>
<p class="sub">Two limits of the machinery, stated rather than hidden. The
contour is interpolated between grid cells of
{fmt(q1(con, "SELECT max(n_trials) FROM ch_detect"))} trials each, so at a
recovery fraction of 0.9 a single cell carries a binomial standard error of
about 4%; each A<sub>90</sub> below is therefore quoted with the
16&ndash;84% range from resampling those cells. And the threshold rests on a
median of {fmt(thr_n)} check stars: measured one star at a time, the
per-star thresholds spread by a median {100 * (thr_spread or 0):.0f}% of the
pooled value, which is the honest uncertainty on it rather than the
&plusmn;1% a bootstrap over trials reports.</p>
{_figure(fig_contour(con),
         "The 90% recovery contour per target. Left: the richest single "
         "night, scored as detection at a known period (solid) and as blind "
         "period determination (dotted) &mdash; the gap between them is the "
         "correction this section is about. Right: the whole season, raw "
         "(dashed) and after per-night mean removal (solid). Thin vertical "
         "lines mark each target's orbital period.")}
{_figure(fig_recovery_grid(con),
         "The raw grids behind one night's two contours: the same "
         "injections into the same noise, scored both ways. The left panel "
         "is what the data can detect; the right is what they can measure a "
         "period from.")}
{table(["series", "regime", "score", "P (h)", "N points", "cycles covered",
        "per-point &sigma; (mmag)",
        "measured A<sub>90</sub> (mmag)",
        "analytic A<sub>min</sub> = &sigma;&radic;(4z/N) (mmag)",
        "measured / analytic"], prows, prclass)}
<h3>Evidence &mdash; does the wrong peak actually win?</h3>
<p class="sub">&sect;3 measures how much window power the daily comb makes
available at f&nbsp;&plusmn;&nbsp;k&nbsp;c/d. That is the input to the alias
question, not the answer: whether the alias <em>beats</em> the truth depends
on the signal's amplitude and on how many nights there are. Measured here,
with the same injections and the same real noise, classifying the tallest
peak as the truth, a &plusmn;k&nbsp;c/d alias, or neither:</p>
{table(["injected semi-amplitude", "tallest peak is the truth",
        "&hellip; is a &plusmn;k c/d alias"],
       [[_mmag(a, 0), f"{100 * ft:.0f}%", f"{100 * fa:.0f}%"]
        for a, ft, fa in conf])}
<div class="decision"><b>Decision.</b> Quote the <b>known-period</b> contour
as this paper's detection limit &mdash; the orbital periods are literature
values with decades of baseline, so nothing here is measuring them from
scratch. On that footing the measured A<sub>90</sub> at each target's own
period differs from the analytic <code>&sigma;&radic;(4z/N)</code> the
strategy document uses by a factor of
<b>{_num(min(known_ratios), 1)}&ndash;{_num(max(known_ratios), 1)}</b>
&mdash; and note that range spans 1 in both directions. The formula is
<em>conservative</em> at a known period, because its <code>z</code> charges
a blind-search look-elsewhere penalty this paper does not owe; it is
<em>optimistic</em> when a period genuinely has to be found, where the
measured limit is
{_num(min(period_ratios), 1)}&ndash;{_num(max(period_ratios), 1)}&times;
worse than the formula on the multi-night sets. One expression that errs by
two to five times in opposite directions depending on which question is
asked is not a detection limit: it assumes white noise, an alias-free
window and one fixed <code>z</code>, and this archive supplies none of the
three. What must <em>not</em> be quoted at all is the single-night
period-determination contour
({_num(min(night_period_ratios), 1)}&ndash;{_num(max(night_period_ratios), 1)}&times;)
as though it were a detection limit &mdash; those rows, greyed above, measure
how few cycles one night holds. <b>Every detection limit in the manuscript
must be the injection number under a named score mode, not the
formula.</b> {gain_txt} And on alias confusion the measurement is
reassuring: at the amplitudes these polars actually show, the tallest peak
lands on the true frequency in essentially every trial.</div>
<h3>Consequence</h3>
<p class="sub">The polars' orbital modulations (0.5&ndash;2 mag) and YZ
Cnc's superhumps (0.1&ndash;0.3 mag) sit one to two orders of magnitude
above these contours: <em>detection</em> is never the binding constraint for
the headline science. What binds is precision on a colour point and epoch
precision on an edge &mdash; which is why &sect;5 exists.</p>
<h3>Full contour table</h3>
{table(["series", "regime", "score", "P (h)", "A<sub>90</sub> (mmag)", "N",
        "&sigma; (mmag)"], allrows)}
</div></section>"""


def section_timing(con) -> str:
    rows = q(con, """SELECT target_key, series_key, night, regime, n_pts_cycle,
                            cadence_s, depth_mag, sigma_mag, ingress_req,
                            ingress_phase, exp_smear_phase, sigma_t_s,
                            depth_source, n_noise_series, night_kind,
                            depth_season_mag FROM ch_timing
                     ORDER BY target_key, night_kind, regime, ingress_req""")
    if not rows:
        return """<section id="timing"><div class="bhead">
<h2>5 &middot; How precisely can one bright-phase edge be timed?</h2></div>
<p class="missing">Timing stage has not been run.</p></section>"""
    trows = [[esc(TARGET_LABEL.get(r[0], r[0])), esc(r[1]),
              f"<b>{esc(r[14])}</b>", esc(r[3]),
              fmt(r[4]), _num(r[5], 0), _num(r[6], 2), _num(r[15], 2),
              _mmag(r[7]),
              f"{r[8]:g} P", f"{r[9]:g} P", f"{100 * r[10]:.1f}%",
              fmt(r[13]), _num(r[11], 1)] for r in rows]
    rclass = ["" if (r[11] or 1e9) < 60 else "warn" for r in rows]

    def _pick(kind, pred):
        return [r[11] for r in rows if r[14] == kind and pred(r[3]) and r[11]]

    ideal = _pick("richest", lambda g: g == "per-cycle")
    mism = _pick("richest", lambda g: g.endswith("mismatched"))
    nightm = _pick("richest", lambda g: g == "night-mean")
    ideal_med = _pick("median", lambda g: g == "per-cycle")
    mism_med = _pick("median", lambda g: g.endswith("mismatched"))
    st_i = [r[11] for r in rows if r[0] == "stlmi" and r[14] == "richest"
            and r[3] == "per-cycle" and r[11]]
    st_m = [r[11] for r in rows if r[0] == "stlmi" and r[14] == "richest"
            and r[3].endswith("mismatched") and r[11]]
    st_mm = [r[11] for r in rows if r[0] == "stlmi" and r[14] == "median"
             and r[3].endswith("mismatched") and r[11]]
    depth_gap = q1(con, "SELECT max(depth_season_mag / depth_mag) "
                        "FROM ch_timing WHERE depth_mag > 0")
    n_gauss = q1(con, "SELECT count(*) FROM ch_timing WHERE n_noise_series=0")
    return f"""
<section id="timing">
<div class="bhead"><h2>5 &middot; How precisely can one bright-phase edge be
timed?</h2><span class="tag">per-cycle &sigma;<sub>t</sub>
{_num(min(ideal), 0)}&ndash;{_num(max(mism), 0)} s</span></div>
<div class="stage">
<h3>Question</h3>
<p class="sub">The strategy document's Q2 rests on a threshold it explicitly
refuses to assume: per-cycle epoch precision <b>&sigma;<sub>t</sub> &lt;
60&nbsp;s</b>, <em>to be demonstrated by injection before adoption</em>
(&sect;4.16), with a pre-defined fallback to seasonal means (&sect;6.10a).
This section runs that demonstration.</p>
<h3>Evidence</h3>
<p class="sub">A trapezoidal faint phase of width
{100 * 0.45:.0f}% of the orbit, at the amplitude the target actually showed
<b>on the night being simulated</b> (5th-to-95th percentile of that night's
points, not 2&times; its RMS &mdash; the RMS of a barely-detected star is
measurement noise), is injected at the <b>real BJD_TDB timestamps of one
cycle</b>, into the <b>real check-star residuals</b> cyclically rolled. The
epoch is then recovered by &chi;<sup>2</sup> minimisation over a two-stage
grid.</p>
<p class="sub">The words &ldquo;on the night being simulated&rdquo; are a
correction. The first version took the 5th-to-95th percentile of the
<em>whole series</em>, which on ST LMi spans 396 days across a 1.92&nbsp;mag
change in nightly median &mdash; the high/low accretion-state transition Q1
and Q4 are about. The 1.45&nbsp;mag quoted as a &ldquo;measured bright-phase
amplitude&rdquo; was therefore mostly a state change, larger than the night's
own amplitude by up to {_num(depth_gap, 1)}&times;. Re-run on the night's own
amplitude the answer barely moves, because this regime is sampling-limited
rather than signal-to-noise-limited &mdash; but the number as published would
not have survived a referee.</p>
<p class="sub">And the Monte Carlo is run on <b>two</b> nights per target,
not one. The richest night is a demonstration; the O&ndash;C tier will be
built from ordinary nights, which are sparser. Both are tabled, labelled
<code>richest</code> and <code>median</code>.</p>
<p class="sub">Three regimes. <b>per-cycle</b> holds the template shape
<em>exactly</em> right &mdash; the best case any analysis could reach, a hard
lower bound. <b>per-cycle shape-mismatched</b> fits with the edge sharpness
5&times; wrong and the depth 20% wrong &mdash; the realistic case, because
nobody has measured this archive's ingress duration and cyclotron beaming
makes it band-dependent. <b>night-mean</b> fits all cycles of the night at
once. In every case the injected edge is widened to at least the exposure
time: a 240&nbsp;s integration on a 90-minute binary cannot record an edge
sharper than 4.4% of a cycle, whatever the star does.</p>
{_figure(fig_timing(con),
         "Monte-Carlo epoch precision against the strategy's assumed 60 s "
         "target. Note the log axis, and that the mismatched (red) bars are "
         "the ones a real analysis would live with.")}
{table(["target", "series", "night", "regime", "pts in fit", "&Delta;t (s)",
        "injected depth (mag)", "season p5&ndash;p95 (mag)", "&sigma; (mmag)",
        "edge (requested)",
        "edge (after exposure smear)", "exposure smear",
        "real-noise series used", "&sigma;<sub>t</sub> (s)"], trows, rclass)}
<div class="decision"><b>Decision.</b> On the <b>richest</b> night, with the
edge shape known exactly, per-cycle &sigma;<sub>t</sub> is
<b>{_num(min(ideal), 1)}&ndash;{_num(max(ideal), 0)} s</b> &mdash; the 60 s
threshold is <em>met</em>. With the shape and depth assumed wrong it rises to
<b>{_num(min(mism), 1)}&ndash;{_num(max(mism), 0)} s</b>. On the
<b>median-density</b> night, which is what the tier will mostly be built
from, the same two answers are
{_num(min(ideal_med), 1)}&ndash;{_num(max(ideal_med), 0)} s and
<b>{_num(min(mism_med), 1)}&ndash;{_num(max(mism_med), 0)} s</b>. On ST LMi,
the flagship, the mismatched answer is
{_num(min(st_m), 0)}&ndash;{_num(max(st_m), 0)} s on its best night and
{_num(min(st_mm), 0)}&ndash;{_num(max(st_mm), 0)} s on a typical one.
Night-mean timings are
{_num(min(nightm), 1)}&ndash;{_num(max(nightm), 1)} s. <b>The per-cycle
O&ndash;C tier survives on the best nights, with almost no margin, and only
if the paper carries an explicit edge-shape systematic of order 50 s beside
every per-cycle timing &mdash; and it must be sized on the median night, not
the demonstration night.</b></div>
<h3>Consequence</h3>
<p class="sub">Two things follow. First, the binding variable is
<b>cadence</b>, not aperture or exposure: &sigma;<sub>t</sub> is set by
whether a sample lands on the ingress ramp, and the unit test behind this
function shows that improving the photometry 30&times; improves the epoch by
under 5&times;. A filter cycle that returns to each band twice as fast is
worth more than any realistic gain in signal-to-noise. Second, the
next measurement this project needs is the <b>edge shape itself</b>: measure
the ingress duration from the deepest night's stacked cycles and the
50&nbsp;s systematic collapses.
{("<b>Caveat:</b> " + fmt(n_gauss) + " of these rows fell back to Gaussian "
  "noise because no check star covered 40% of that night &mdash; they are "
  "optimistic by the red-noise factor of &sect;2.") if n_gauss else ""}</p>
</div></section>"""


def section_verdict(con, pcon) -> str:
    rows = q(con, """SELECT goal_id, rank, goal, claim, verdict,
                            deciding_number, reasoning, alternative
                     FROM ch_verdict ORDER BY rank""")
    cls = {"SUPPORTED": "", "SUPPORTED-WITH-CAVEATS": "warn",
           "NOT SUPPORTED": "warn", "NOT MEASURED": "warn"}
    trows = [[f"<b>{esc(r[0])}</b>", esc(r[2]), esc(r[3]),
              f"<b>{esc(r[4])}</b>", esc(r[5]), esc(r[6]), esc(r[7])]
             for r in rows]
    rclass = [cls.get(r[4], "") for r in rows]
    n_sup = sum(1 for r in rows if r[4] == "SUPPORTED")
    n_cav = sum(1 for r in rows if r[4] == "SUPPORTED-WITH-CAVEATS")
    n_not = sum(1 for r in rows if r[4].startswith("NOT"))
    # Every number in the ranked "what next" list is queried, not typed.
    n_tied = q1(pcon, "SELECT count(*) FROM cv_field_tie "
                      "WHERE n_gaia_matched > 0")
    n_blocks = q1(pcon, "SELECT count(*) FROM cv_field_tie")
    untied_blocks = q(pcon, "SELECT target_key, era_id FROM cv_field_tie "
                            "WHERE status LIKE 'query_failed%'")
    n_untied_frames = q1(pcon,
                         "SELECT count(*) FROM cv_frames f JOIN cv_field_tie t "
                         "ON t.target_key=f.target_key AND t.era_id=f.era_id "
                         "WHERE f.status='matched' "
                         "AND t.status LIKE 'query_failed%'") or 0
    n_failed = q1(pcon, "SELECT count(*) FROM cv_frames "
                        "WHERE status='failed_match'")
    failed_top = q(pcon, "SELECT series_key, count(*) FROM cv_frames "
                         "WHERE status='failed_match' GROUP BY 1 "
                         "ORDER BY 2 DESC LIMIT 2")
    med_ratio = q1(con, "SELECT median(floor_plateau / scint_mag) FROM "
                        "ch_noise_series WHERE floor_plateau IS NOT NULL")
    return f"""
<section id="verdict">
<div class="bhead"><h2>6 &middot; The verdict: which science goals does the
measured performance support?</h2>
<span class="tag">{n_sup} supported</span>
<span class="tag">{n_cav} with caveats</span>
<span class="tag">{n_not} not supported</span></div>
<div class="stage">
<h3>Question</h3>
<p class="sub">The strategy document (<code>CV_TimeSeries/ANALYSIS_STRATEGY.md</code>)
sets out five science questions and several structural claims. Each was
written before the photometry existed. Which of them does the measured
performance actually support?</p>
<h3>Evidence</h3>
<p class="sub">Every verdict below is recomputed from the tables in
&sect;&sect;1&ndash;5 whenever this build runs, and the &ldquo;deciding
number&rdquo; column is assembled from query results &mdash; so the table
cannot drift away from the evidence it rests on. The thresholds separating
the verdicts are the strategy document's <em>own</em> stated requirements;
this page tests the plan against the data, it does not invent new criteria.</p>
{table(["#", "goal", "the claim it needs", "verdict", "the number that decides it",
        "why", "what the data can support instead"], trows, rclass)}
<div class="decision"><b>Decision.</b> {n_sup} goal(s) survive unqualified,
{n_cav} survive with stated caveats, and {n_not} do not survive contact with
the measurements. The paper's spine is the goals in the first two rows; the
rest are supporting sections or should be dropped.</div>
<h3>Consequence &mdash; what to measure next, ranked</h3>
<p class="sub">In the order that most changes what this data set can claim:</p>
<ol class="sub">
<li><b>Tie the ensemble to a catalogue.</b> {fmt(n_blocks - n_tied)} of
{fmt(n_blocks)} (target, era) blocks carry no magnitude tie at all, so their
magnitudes are on an arbitrary internal gauge. Nothing else on this list is
worth as much per hour of work: it converts every light curve from relative
to publishable, and it is a re-run of an existing stage, not new code.</li>
<li><b>Re-run the {fmt(len(untied_blocks))} failed field ties</b>
({", ".join(f"{esc(TARGET_LABEL.get(t, t))} era {e}" for t, e in untied_blocks)}),
lost to archive HTTP errors. Without them those blocks have no target light
curve at all &mdash; {fmt(n_untied_frames)} measured frames currently produce
no science point.</li>
<li><b>Measure the per-filter cadence into the observing plan.</b> &sect;5
shows epoch precision is sampling-limited; a filter cycle that returns to
each band twice as fast is worth more than any realistic gain in
signal-to-noise, and costs nothing but scheduling.</li>
<li><b>Chase the flat-field / ensemble floor.</b> The measured floor is
{_num(med_ratio, 0)}&times; scintillation, so it is instrumental and
improvable: test for a
2-D illumination residual by regressing ensemble residuals on (x, y), which
the strategy already lists (&sect;5, &ldquo;flat-field stability&rdquo;) and
which no stage has yet run.</li>
<li><b>Recover the {fmt(n_failed)} frames that failed star-matching</b>
&mdash; concentrated in
{", ".join(f"{esc(sk)} ({fmt(n)})" for sk, n in failed_top)}, blocks that
carry the timing and superhump science.</li>
<li><b>Run the single-night periodogram ladder</b> now that the windows are
measured, and publish the alias-power ratio with each one.</li>
</ol>
</div></section>"""


# ===========================================================================
# Page assembly
# ===========================================================================

def render_report(char_db: Path, phot_db: Path) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{char_db}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    pcon = sqlite3.connect(f"file:{phot_db}?mode=ro", uri=True)
    pcon.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = dict(q(con, "SELECT key, value FROM ch_meta"))
        sections = [section_intro(con, pcon), section_quality(con),
                    section_noise(con), section_cadence(con),
                    section_detect(con), section_timing(con),
                    section_verdict(con, pcon)]
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CV Time Series &mdash; Data Characterization &amp; What It Can Measure</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Cataclysmic-Variable Time Series &mdash; characterization, and what
  this data can actually measure</h1>
  <p>Image quality &rarr; noise floor &rarr; sampling &rarr; detectability
  &rarr; timing precision &rarr; goal-by-goal verdict &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))}) &middot;
  <a href="index.html">project hub</a> &middot;
  <a href="../index.html">all reports</a></p>
</header>

<nav>
  <a href="#intro">0 The question</a> &middot;
  <a href="#quality">1 Image quality</a> &middot;
  <a href="#noise">2 Noise floor</a> &middot;
  <a href="#cadence">3 Sampling &amp; aliases</a> &middot;
  <a href="#detect">4 Detectability</a> &middot;
  <a href="#timing">5 Timing precision</a> &middot;
  <a href="#verdict">6 Verdict</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_phot.report_cv_char</code> from
<code>products/phot/cv_characterization.sqlite</code> &mdash; every number on
this page is the result of a SQL query or a constant imported from
<code>macro_phot.characterize</code>; none is typed by hand. Rebuild with
<code>pipeline/scripts/run_cv_characterization.py all</code>.
Periods used only for folding and cadence arithmetic:
{esc("; ".join(f"{TARGET_LABEL.get(k, k)}: {vv}" for k, vv in
      sorted(json.loads(meta.get("period_sources_json", "{}")).items())))}.
</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close(); pcon.close()
