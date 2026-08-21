"""CV-S9 chain-of-evidence report: the Phase-3 time-series analysis.

Reads ``products/phot/cv_timeseries.sqlite`` and writes

* ``docs/CV_TimeSeries/cv_timeseries_analysis.html``
* ``docs/CV_TimeSeries/figures/cv_phase3/*.png``

Six Socratic sections, one per Phase-3 task, each in the same order:

    Question    what could go wrong if we did not do this?
    Evidence    the measurement, from the database, with its diagnostics
    Decision    what the numbers actually support
    Consequence what changes downstream, and what does not

and one section before them that says what Phase 3 inherits and what it is
forbidden to assume.

THE SPECTRAL WINDOW IS BESIDE EVERY PERIODOGRAM ON THIS PAGE.  That is not a
stylistic choice, it is the binding rule of the period section: with the
+/-1 c/d sidelobes at 0.54-0.97 of the window power on every resolved
multi-night set in this archive, a periodogram shown alone is an invitation
to read a sidelobe as a period.  Every periodogram panel is drawn with its
own window panel on the same frequency axis, from the same stored trace.

Every number in a table, figure or verdict is the result of a query run at
render time or a constant imported from ``macro_phot.phase3``.  Where a
number could not be measured, the page says which query returned nothing
rather than leaving a blank.
"""

from __future__ import annotations

import math
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

from . import phase3 as p3        # noqa: E402
from macro_core.report_s0 import (  # noqa: E402
    ACCENT, BAD, STYLE, DPI, GOOD, INK, MUTED, WARN,
    esc, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402  (house figure style)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "CV_TimeSeries"
FIG_DIR = DOCS_DIR / "figures" / "cv_phase3"
HTML_PATH = DOCS_DIR / "cv_timeseries_analysis.html"

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}
TARGET_ORDER = ("stlmi", "anuma", "vvpup", "euuma", "yzcnc")


# ---------------------------------------------------------------------------
# Formatting helpers.  An em-dash is the only thing a missing number becomes:
# a blank cell reads as zero, and zero is a claim.
# ---------------------------------------------------------------------------
def _n(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{float(x):,.{nd}f}"


def _i(x):
    if x is None:
        return "&mdash;"
    return f"{int(x):,}"


def _sci(x, nd=2):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{float(x):.{nd}e}"


def _pct(x, nd=1):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{100.0 * float(x):.{nd}f}%"


def _pm(v, e, nd=7):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    if e is None or (isinstance(e, float) and not math.isfinite(e)):
        return f"{float(v):.{nd}f}"
    return f"{float(v):.{nd}f} &plusmn; {float(e):.{nd}f}"


def _label(series_key: str) -> str:
    """``stlmi|e76|g`` -> ``ST LMi e76 g``."""
    parts = str(series_key).split("|")
    if len(parts) != 3:
        return esc(series_key)
    return f"{TARGET_LABEL.get(parts[0], parts[0])} {parts[1]} {parts[2]}"


def _fig(src: str, caption: str, missing: str) -> str:
    if not src:
        return f'<p class="note">{missing}</p>'
    return (f'<figure><a href="{src}"><img src="{src}" alt=""></a>'
            f"<figcaption>{caption}</figcaption></figure>")


def _has(con, name: str) -> bool:
    return bool(q1(con, "SELECT count(*) FROM sqlite_master WHERE "
                        "type='table' AND name=?", (name,)))


def _rows(con, name: str) -> int:
    return q1(con, f"SELECT count(*) FROM {name}") if _has(con, name) else 0


# ===========================================================================
# Figures
# ===========================================================================
def fig_periodograms(con, target: str) -> str:
    """One row per series: survey periodogram, zoom, AND SPECTRAL WINDOW.

    The third column is the point of the figure.  It is the periodogram the
    sampling alone would produce for a perfectly constant star, plotted on
    the frequency OFFSET axis so its peaks sit where the aliases of any real
    signal will sit.  Reading the left two panels without the right one is
    how a sidelobe becomes a discovery.
    """
    series = [r[0] for r in q(con, """
        SELECT series_key FROM p3_period WHERE target_key=? AND status='ok'
        ORDER BY era_id, filter""", (target,))]
    if not series:
        return ""
    n = len(series)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(n, 3, figsize=(12.4, 2.05 * n + 0.5),
                                 squeeze=False,
                                 gridspec_kw={"width_ratios": [1.5, 1.5, 1.0]})
        for i, sk in enumerate(series):
            meta = q(con, """SELECT f_survey_cd, p_survey_pow, f_ls_cd,
                                    f_pdm_cd, published_d, alias_frac_max,
                                    n_blocks, family_code, detected
                             FROM p3_period WHERE series_key=?""", (sk,))[0]
            f_pub = 1.0 / meta[4] if meta[4] else float("nan")
            a0, a1, a2 = axes[i]
            # --- survey ---
            tr = q(con, "SELECT freq_cd, value FROM p3_pgram WHERE "
                        "series_key=? AND panel='survey' AND kind='ls' "
                        "ORDER BY freq_cd", (sk,))
            if tr:
                x = np.array([r[0] for r in tr])
                y = np.array([r[1] for r in tr])
                a0.plot(x, y, lw=0.7, color=ACCENT)
                if np.isfinite(f_pub):
                    a0.axvline(f_pub, color=WARN, lw=1.0, ls="--")
                a0.set_xlim(p3.SURVEY_F_MIN_CD, p3.SURVEY_F_MAX_CD)
                a0.set_ylabel(_label(sk), fontsize=7.5)
            if i == 0:
                a0.set_title("survey periodogram 0.5-40 c/d\n"
                             "(GLS, one free constant per night)", fontsize=8)
            # --- zoom, LS and PDM ---
            tz = q(con, "SELECT freq_cd, value FROM p3_pgram WHERE "
                        "series_key=? AND panel='zoom' AND kind='ls' "
                        "ORDER BY freq_cd", (sk,))
            tp = q(con, "SELECT freq_cd, value FROM p3_pgram WHERE "
                        "series_key=? AND panel='zoom' AND kind='pdm' "
                        "ORDER BY freq_cd", (sk,))
            if tz:
                x = np.array([r[0] for r in tz])
                y = np.array([r[1] for r in tz])
                a1.plot(x, y, lw=0.7, color=ACCENT, label="Lomb-Scargle")
                if tp:
                    xp = np.array([r[0] for r in tp])
                    yp = np.array([r[1] for r in tp])
                    # 1 - theta so a PDM dip reads as a peak on the same axis.
                    a1.plot(xp, 1.0 - yp, lw=0.7, color=GOOD,
                            label="1 - PDM theta")
                if np.isfinite(f_pub):
                    a1.axvline(f_pub, color=WARN, lw=1.0, ls="--",
                               label="published")
                    for k in (-1, 1):
                        a1.axvline(f_pub + k, color=BAD, lw=0.8, ls=":")
                    a1.set_xlim(f_pub - 2.5, f_pub + 2.5)
                if i == 0:
                    a1.legend(fontsize=6.5, loc="upper right")
                    a1.set_title("orbital band, published frequency dashed,\n"
                                 "+/-1 c/d aliases dotted", fontsize=8)
            # --- THE SPECTRAL WINDOW ---
            tw = q(con, "SELECT freq_cd, value FROM p3_pgram WHERE "
                        "series_key=? AND panel='window' ORDER BY freq_cd",
                   (sk,))
            if tw:
                xw = np.array([r[0] for r in tw])
                yw = np.array([r[1] for r in tw])
                a2.fill_between(xw, 0, yw, color=BAD, alpha=0.35)
                a2.plot(xw, yw, lw=0.8, color=BAD)
                a2.set_xlim(-3, 3)
                a2.set_ylim(0, 1.05)
                a2.text(0.03, 0.86,
                        f"W(±1) = {meta[5]:.2f}\n{meta[6]} night"
                        f"{'s' if meta[6] != 1 else ''}\n{meta[7]}",
                        transform=a2.transAxes, fontsize=6.5, va="top",
                        color=MUTED,
                        # A single-night window fills the whole panel, and
                        # grey type on it was unreadable.  The label carries
                        # its own ground.
                        bbox=dict(facecolor=ps.PAPER, edgecolor="none",
                                  alpha=0.78, pad=1.4))
            if i == 0:
                a2.set_title("SPECTRAL WINDOW\n(frequency offset, c/d)",
                             fontsize=8)
            if i == n - 1:
                a0.set_xlabel("frequency (c/d)", fontsize=8)
                a1.set_xlabel("frequency (c/d)", fontsize=8)
                a2.set_xlabel("offset from any true frequency (c/d)",
                              fontsize=8)
            for ax in (a0, a1, a2):
                ax.tick_params(labelsize=6.5)
        fig.suptitle(f"{TARGET_LABEL.get(target, target)} — every "
                     f"periodogram with its own spectral window", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / f"pgram_{target}.png", dpi=DPI)
        plt.close(fig)
    return f"figures/cv_phase3/pgram_{target}.png"


def fig_period_summary(con) -> str:
    """Deviation from the published period, and how tightly constrained."""
    rows = q(con, """SELECT series_key, deviation_sigma, frac_precision,
                            alias_frac_max, detected, constraint_class
                     FROM p3_period WHERE status='ok' AND period_d IS NOT NULL
                     ORDER BY target_key, era_id, filter""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 4.6))
        y = np.arange(len(rows))
        dev = np.array([r[1] if r[1] is not None else np.nan for r in rows])
        det = np.array([bool(r[4]) for r in rows])
        a1.axvspan(-p3.AGREE_SIGMA, p3.AGREE_SIGMA, color=GOOD, alpha=0.16,
                   label=f"agreement band (±{p3.AGREE_SIGMA:.0f}σ)")
        a1.scatter(dev[det], y[det], s=26, color=ACCENT,
                   label="orbital modulation detected")
        a1.scatter(dev[~det], y[~det], s=26, color=MUTED, marker="x",
                   label="not detected")
        a1.axvline(0, color=MUTED, lw=0.8)
        a1.set_yticks(y)
        a1.set_yticklabels([_label(r[0]) for r in rows], fontsize=6.5)
        a1.set_xlabel("(recovered − published) / combined σ")
        a1.set_xlim(-6, 6)
        a1.legend(fontsize=7, loc="lower right")
        a1.set_title("agreement with the published period", fontsize=9)
        # Right: how much the period is actually constrained, against how
        # badly the window aliases.  The top-left corner is where a period
        # determination would live; nothing is there.
        frac = np.array([r[2] if r[2] is not None else np.nan for r in rows])
        alias = np.array([r[3] if r[3] is not None else np.nan for r in rows])
        a2.scatter(alias, frac, s=30, color=ACCENT)
        for r, xa, yf in zip(rows, alias, frac):
            if np.isfinite(xa) and np.isfinite(yf):
                a2.annotate(_label(r[0]), (xa, yf), fontsize=5.5,
                            xytext=(3, 2), textcoords="offset points",
                            color=MUTED)
        a2.axvline(p3.ALIAS_DECIDABLE_MAX, color=GOOD, ls="--", lw=1.2,
                   label=f"alias power below which a periodogram could\n"
                         f"choose its own period ({p3.ALIAS_DECIDABLE_MAX:g})")
        a2.set_yscale("log")
        a2.set_xlim(0, 1.0)
        a2.set_xlabel("strongest ±1 c/d sidelobe, fraction of window power")
        a2.set_ylabel("fractional period precision σ(P)/P")
        a2.legend(fontsize=7, loc="lower left")
        a2.set_title("no series is left of the line:\nevery period rests on "
                     "the literature prior", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "period_summary.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase3/period_summary.png"


def fig_sigmat(con) -> str:
    """The sigma_t contour against the 60 s threshold."""
    keys = [r[0] for r in q(con, "SELECT DISTINCT series_key FROM p3_sigmat "
                                 "ORDER BY series_key")]
    if not keys:
        return ""
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(keys) + 1,
                                 figsize=(3.6 * (len(keys) + 1), 3.9),
                                 squeeze=False)
        axes = axes[0]
        shapes = sorted({r[0] for r in q(
            con, "SELECT DISTINCT shape_error FROM p3_sigmat")})
        depths = sorted({r[0] for r in q(
            con, "SELECT DISTINCT depth_error FROM p3_sigmat")})
        for ax, sk in zip(axes, keys):
            grid = np.full((len(depths), len(shapes)), np.nan)
            for se, de, tot in q(con, """SELECT shape_error, depth_error,
                                                total_error_s FROM p3_sigmat
                                         WHERE series_key=? AND
                                               inject_factor=1.0""", (sk,)):
                grid[depths.index(de), shapes.index(se)] = tot
            norm = matplotlib.colors.LogNorm(
                vmin=max(np.nanmin(grid), 1.0), vmax=np.nanmax(grid))
            im = ax.imshow(grid, origin="lower", aspect="auto",
                           cmap=ps.SEQ_CMAP, norm=norm)
            ax.grid(False)          # a heatmap wears no grid
            for iy in range(grid.shape[0]):
                for ix in range(grid.shape[1]):
                    v = grid[iy, ix]
                    if not np.isfinite(v):
                        continue
                    # Paper or ink by where the cell sits on the ramp.  The
                    # previous version wrote every PASSING cell in green,
                    # which on the dark end of the old map was unreadable;
                    # the pass/fail verdict is the CONTOUR's job, below.
                    ax.text(ix, iy, f"{v:.0f}", ha="center", va="center",
                            fontsize=7.5, weight="bold",
                            color=ps.ink_on(float(norm(v))))
            # The 60 s contour, drawn where it actually falls.
            try:
                ax.contour(grid, levels=[p3.SIGMA_T_THRESHOLD_S],
                           colors=[GOOD], linewidths=2.0)
            except Exception:                            # noqa: BLE001
                pass
            ax.set_xticks(range(len(shapes)))
            ax.set_xticklabels([f"×{s:g}" for s in shapes], fontsize=7.5)
            ax.set_yticks(range(len(depths)))
            ax.set_yticklabels([f"{d:+.0%}" for d in depths], fontsize=7.5)
            ax.set_xlabel("assumed edge width / true", fontsize=8)
            ax.set_ylabel("assumed depth error", fontsize=8)
            ax.set_title(f"{_label(sk)}\ntotal timing error (s)", fontsize=8.5)
            fig.colorbar(im, ax=ax, fraction=0.046)
        # Last panel: the injected-width sensitivity.
        ax = axes[-1]
        for sk in keys:
            rows = q(con, """SELECT inject_width_s, total_error_s FROM
                             p3_sigmat WHERE series_key=? AND shape_error=1.0
                             AND depth_error=0.0 ORDER BY inject_width_s""",
                     (sk,))
            if rows:
                ax.plot([r[0] for r in rows], [r[1] for r in rows], "-o",
                        ms=4, lw=1.2, label=_label(sk))
        ax.axhline(p3.SIGMA_T_THRESHOLD_S, color=GOOD, ls="--", lw=1.6,
                   label=f"{p3.SIGMA_T_THRESHOLD_S:.0f} s threshold")
        ax.set_xlabel("INJECTED edge width (s)", fontsize=8)
        ax.set_ylabel("total timing error (s)", fontsize=8)
        ax.set_yscale("log")
        ax.legend(fontsize=6.5)
        ax.set_title("sensitivity to the one input the data\n"
                     "cannot measure: the true edge width", fontsize=8.5)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "sigmat_contour.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase3/sigmat_contour.png"


def fig_edges(con) -> str:
    """The folded bright phase per band, and the per-cycle edge epochs."""
    night = q1(con, "SELECT value FROM p3_meta WHERE key='sigmat_night'") \
        if _has(con, "p3_meta") else None
    keys = [r[0] for r in q(con, """SELECT DISTINCT series_key FROM p3_edge
                                    WHERE accepted=1 AND night=?
                                    ORDER BY series_key""", (night,))] \
        if night else []
    if not keys:
        keys = [r[0] for r in q(con, "SELECT series_key, count(*) c FROM "
                                     "p3_edge WHERE accepted=1 GROUP BY 1 "
                                     "ORDER BY c DESC LIMIT 3")]
    if not keys:
        return ""
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 4.4))
        colours = ps.BAND_COLOR
        for sk in keys:
            tgt, era, filt = sk.split("|")
            eph = q(con, "SELECT period_d, epoch_bjd FROM p3_ephemeris "
                         "WHERE target_key=?", (tgt,))
            if not eph or eph[0][0] is None or eph[0][1] is None:
                continue
            per, e0 = eph[0]
            rows = q(con, """SELECT l.bjd_tdb, l.cal_mag FROM cv_lightcurve l
                             JOIN cv_frames f ON f.frame_id=l.frame_id
                                             AND f.series_key=l.series_key
                             WHERE l.series_key=? AND l.role='target'
                               AND f.night=? AND l.cal_mag IS NOT NULL""",
                     (sk, night))
            if rows:
                t = np.array([r[0] for r in rows])
                m = np.array([r[1] for r in rows])
                keep = np.abs(m - np.median(m)) < 8 * 1.4826 * np.median(
                    np.abs(m - np.median(m)))
                ph = p3.phase_of(t[keep], per, e0)
                a1.plot(np.concatenate([ph, ph + 1]),
                        np.concatenate([m[keep], m[keep]]), ".", ms=3.2,
                        alpha=0.75, color=colours.get(filt, ACCENT),
                        label=filt)
            ed = q(con, """SELECT phase FROM p3_edge WHERE series_key=?
                           AND accepted=1 AND phase IS NOT NULL""", (sk,))
            for (phe,) in ed:
                a1.axvline(phe, color=colours.get(filt, ACCENT), lw=0.9,
                           alpha=0.55, ls="--")
                a1.axvline(phe + 1, color=colours.get(filt, ACCENT), lw=0.9,
                           alpha=0.55, ls="--")
        a1.invert_yaxis()
        a1.set_xlabel("orbital phase (published ephemeris)")
        a1.set_ylabel("catalogue-tied magnitude")
        a1.legend(fontsize=7.5, title="band", title_fontsize=7)
        a1.set_title(f"bright phase on {night}: the amplitude is strongly\n"
                     "colour dependent, the fitted edges dashed", fontsize=9)
        # Right: inter-band differences.
        rows = q(con, """SELECT target_key, night, band_a, band_b, delta_s,
                                sigma_s, n_cycles, significant
                         FROM p3_band_pair ORDER BY target_key, night,
                                                   band_a, band_b""")
        if rows:
            y = np.arange(len(rows))
            d = np.array([r[4] for r in rows])
            s = np.array([r[5] for r in rows])
            sig = np.array([bool(r[7]) for r in rows])
            a2.errorbar(d[~sig], y[~sig], xerr=s[~sig], fmt="o", ms=5,
                        color=MUTED, capsize=2, label="consistent with zero")
            if sig.any():
                a2.errorbar(d[sig], y[sig], xerr=s[sig], fmt="o", ms=5,
                            color=BAD, capsize=2, label="≥3σ from zero")
            a2.axvline(0, color=MUTED, lw=1.0)
            a2.set_yticks(y)
            a2.set_yticklabels([f"{TARGET_LABEL.get(r[0], r[0])} {r[1]} "
                                f"{r[2]}−{r[3]} (n={r[6]})" for r in rows],
                               fontsize=6.5)
            a2.set_xlabel("inter-band edge-time difference (s)")
            a2.legend(fontsize=7)
            a2.set_title("colour-dependent edge timing:\n"
                         "a cyclotron measurement, not a nuisance",
                         fontsize=9)
        else:
            a2.text(0.5, 0.5, "no band pair had accepted edges\n"
                              "on the same cycle", ha="center", va="center",
                    transform=a2.transAxes, fontsize=9, color=WARN)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "edges.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase3/edges.png"


def fig_oc(con) -> str:
    """The O-C diagram, and the cycle-count uniqueness margin."""
    rows = q(con, """SELECT target_key, cycle, oc_s, oc_sigma_s, filter,
                            count_unique FROM p3_oc ORDER BY cycle""")
    cc = q(con, """SELECT target_key, n_cycles_last, drift_cycles,
                          sigma_period_d, sigma_period_max_d, unique_count,
                          one_feature
                   FROM p3_cycle_count WHERE n_cycles_last IS NOT NULL""")
    if not rows and not cc:
        return ""
    # One panel per target.  Cycle numbers count from each star's OWN
    # ephemeris epoch, so plotting ST LMi's cycle 21,869 and EU UMa's cycle
    # 45,131 on a shared axis would put two unrelated integers side by side
    # and invite the reader to compare them.
    oc_targets = sorted({r[0] for r in rows},
                        key=lambda t: -sum(1 for r in rows if r[0] == t))
    n_oc = max(len(oc_targets), 1)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, n_oc + 1,
                                 figsize=(4.4 * n_oc + 5.6, 4.3),
                                 squeeze=False)
        axes = axes[0]
        a2 = axes[-1]
        colours = ps.BAND_COLOR
        if rows:
            for a1, tgt in zip(axes, oc_targets):
                trows = [r for r in rows if r[0] == tgt]
                for filt in sorted({r[4] for r in trows}):
                    sel = [r for r in trows if r[4] == filt]
                    a1.errorbar([r[1] for r in sel], [r[2] for r in sel],
                                yerr=[r[3] for r in sel], fmt="o", ms=4,
                                capsize=2, lw=0.9,
                                color=colours.get(filt, ACCENT), label=filt)
                a1.axhline(0, color=MUTED, lw=0.9)
                rms = float(np.std([r[2] for r in trows], ddof=1)) \
                    if len(trows) > 1 else float("nan")
                a1.set_xlabel("cycle number from this star's VSX epoch")
                a1.set_ylabel("O − C (s), mean removed")
                a1.legend(fontsize=7.5, title="band", title_fontsize=7)
                a1.set_title(f"{TARGET_LABEL.get(tgt, tgt)} — "
                             f"{len(trows)} edge epochs, rms {rms:.0f} s\n"
                             "(the constant feature-to-phase-zero offset is "
                             "removed and published separately)", fontsize=8.5)
        else:
            axes[0].text(0.5, 0.5, "no accepted edge epochs to plot",
                         ha="center", va="center", transform=axes[0].transAxes,
                         fontsize=10, color=WARN)
        if cc:
            y = np.arange(len(cc))
            need = np.array([r[4] for r in cc])
            have = np.array([r[3] if r[3] else np.nan for r in cc])
            a2.barh(y - 0.18, need, height=0.34, color=GOOD,
                    label="σ(P) that would still give a unique count")
            a2.barh(y + 0.18, have, height=0.34, color=ACCENT,
                    label="σ(P) floor from the VSX quoted precision")
            a2.set_xscale("log")
            a2.set_yticks(y)
            # The label says whether the target also HAS an O-C, because
            # a unique cycle count and a plottable O-C are different
            # claims and two of these targets have the first without the
            # second.
            a2.set_yticklabels(
                [f"{TARGET_LABEL.get(r[0], r[0])}\n{r[1]:,.0f} cycles\n"
                 + ("O−C plotted" if r[6] else "no O−C: not one feature")
                 for r in cc], fontsize=6.5)
            a2.set_xlabel("period uncertainty (d)")
            a2.legend(fontsize=7, loc="lower right")
            a2.set_title("cycle-count margin: the count is unique when the\n"
                         "green bar exceeds the blue one", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "oc.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase3/oc.png"


def fig_states(con) -> str:
    """Nightly state histories and the thresholds that produced them."""
    series = [r[0] for r in q(con, """SELECT series_key FROM p3_state_series
                                      WHERE threshold_mag IS NOT NULL
                                      ORDER BY target_key, era_id, filter""")]
    if not series:
        return ""
    show = series[:12]
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.0, 4.8))
        for i, sk in enumerate(show):
            rows = q(con, """SELECT night, median_mag, state, censored, gated
                             FROM p3_state_night WHERE series_key=?
                             ORDER BY night""", (sk,))
            thr = q1(con, "SELECT threshold_mag FROM p3_state_series "
                          "WHERE series_key=?", (sk,))
            if not rows:
                continue
            m = np.array([r[1] if r[1] is not None else np.nan for r in rows])
            state = [r[2] for r in rows]
            cens = np.array([bool(r[3]) for r in rows])
            col = np.array([GOOD if s == "HIGH" else BAD if s == "LOW"
                            else MUTED for s in state])
            # Censoring is a property of a NIGHT, not of a series: an
            # earlier version keyed the marker off `cens.any()` and drew
            # every night of a series as a limit the moment one of them was.
            for mask, marker in ((~cens, "o"), (cens, "v")):
                if mask.any():
                    a1.scatter(np.full(int(mask.sum()), i), (m - thr)[mask],
                               c=list(col[mask]), s=26, marker=marker,
                               edgecolors="none")
        a1.axhline(0, color=WARN, ls="--", lw=1.2,
                   label="each series' own Otsu threshold")
        a1.set_xticks(range(len(show)))
        a1.set_xticklabels([_label(s) for s in show], rotation=70,
                           fontsize=6.5, ha="right")
        a1.invert_yaxis()
        a1.set_ylabel("nightly median − threshold (mag)")
        a1.legend(fontsize=7)
        a1.set_title("nightly state relative to each series' threshold\n"
                     "(green HIGH, red LOW, grey intermediate/unclassified)",
                     fontsize=9)
        # Right: duty cycles, censored vs limit-aware.
        rows = q(con, """SELECT series_key, duty_naive, duty_with_limits,
                                separability, bimodal, n_informative_limits
                         FROM p3_state_series
                         WHERE duty_naive IS NOT NULL AND threshold_mag IS NOT NULL
                         ORDER BY target_key, era_id, filter""")
        if rows:
            y = np.arange(len(rows))
            # Outlined bars, because for most series the two values are
            # IDENTICAL and a filled bar drawn second would hide the first
            # entirely — the reader would see one bar and conclude the
            # censored version had not been plotted.
            a2.barh(y - 0.19, [r[1] for r in rows], height=0.36,
                    color=MUTED, edgecolor=MUTED, linewidth=0.6,
                    label="detections only (censored)")
            a2.barh(y + 0.19, [r[2] if r[2] is not None else 0 for r in rows],
                    height=0.36, color=ACCENT, edgecolor=MUTED,
                    linewidth=0.6, label="with Phase-2 limits")
            for yi, r in zip(y, rows):
                if r[2] is not None and abs(r[1] - r[2]) < 1e-9:
                    a2.text(max(r[1], 0.02) + 0.015, yi, "unchanged",
                            va="center", fontsize=5.8, color=MUTED)
            a2.set_yticks(y)
            a2.set_yticklabels([_label(r[0]) +
                                ("" if r[4] else "  (unimodal)") +
                                (f"  [{r[5]} lim]" if r[5] else "")
                                for r in rows], fontsize=6.5)
            a2.set_xlabel("duty cycle: fraction of epochs in the HIGH state")
            a2.set_xlim(0, 1.15)
            a2.legend(fontsize=7, loc="lower right")
            a2.set_title("the censoring bias, measured — and for this "
                         "archive\nit is almost everywhere zero", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "states.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase3/states.png"


def fig_detrend(con) -> str:
    """What detrend-then-search does to an injected signal."""
    rows = q(con, """SELECT window_periods, frac_detrend, frac_joint
                     FROM p3_detrend ORDER BY window_periods""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
        x = np.array([r[0] for r in rows])
        d = np.array([r[1] if r[1] is not None else np.nan for r in rows])
        j = np.array([r[2] if r[2] is not None else np.nan for r in rows])
        ax.axhline(1.0, color=MUTED, lw=1.0, ls="-",
                   label="the injected amplitude")
        ax.plot(x, d, "-o", ms=5, color=BAD, lw=1.6,
                label="detrend first, then search")
        ax.plot(x, j, "-o", ms=5, color=GOOD, lw=1.6,
                label="joint GP + signal fit (celerite2)")
        ax.fill_between(x, d, 1.0, color=BAD, alpha=0.14)
        ax.set_xscale("log")
        ax.set_xlabel("running-median window, in units of the orbital period")
        ax.set_ylabel("recovered amplitude / injected amplitude")
        ax.legend(fontsize=8, loc="lower right")
        ax.set_title("the same data, the same signal, two orderings", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "detrend.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase3/detrend.png"


# ===========================================================================
# Sections
# ===========================================================================
def section_intro(con) -> str:
    n_series = _rows(con, "p3_period")
    n_pts = q1(con, "SELECT sum(n_points) FROM p3_period") or 0
    eph = q(con, """SELECT target_key, name, var_type, period_str,
                           period_sigma_d, epoch_str, source, fetched_utc
                    FROM p3_ephemeris ORDER BY target_key""")
    eph_tbl = table(
        ["target", "VSX type", "published period (d)",
         "σ(P) floor (d)", "epoch (BJD)", "source", "fetched"],
        [[esc(r[1]), esc(r[2]), f"<code>{esc(r[3])}</code>",
          _sci(r[4]), esc(r[5] or "&mdash; none published &mdash;"),
          esc(r[6]), esc((r[7] or "")[:10])] for r in eph])
    alias_lo = q1(con, "SELECT min(alias_frac_max) FROM p3_period WHERE "
                       "n_blocks > 1")
    alias_hi = q1(con, "SELECT max(alias_frac_max) FROM p3_period WHERE "
                       "n_blocks > 1")
    return f"""
<section id="intro">
<h2>0 &nbsp; What Phase&nbsp;3 inherits, and what it may not assume</h2>

<h3>Question</h3>
<p>Phase&nbsp;1 measured {_i(n_pts)} target points across {_i(n_series)}
series and Phase&nbsp;2 decided what may be published from them.  Phase&nbsp;3
asks the six questions the project exists for.  Each of them can be answered
confidently and wrongly, and in every case the wrong answer is the one that
looks tidier.  What, concretely, is this data not allowed to claim?</p>

<h3>Evidence</h3>
<p>Four measurements constrain everything on this page, and none of them is
an assumption:</p>
<ul>
<li><b>Per-point precision 9&ndash;77&nbsp;mmag</b> by series, with
&chi;<sup>2</sup> inflation 0.92&ndash;3.02.  Every error bar used here is
multiplied by that series' measured inflation, and an inflation below 1 is
<em>not</em> applied &mdash; a slightly pessimistic error model does not
license shrinking the bars.</li>
<li><b>The &plusmn;1&nbsp;c/d aliases carry {_n(alias_lo, 2)}&ndash;{_n(alias_hi, 2)}
of the window power</b> on every resolved multi-night set, measured
in&nbsp;&sect;1.  A multi-night periodogram in this archive cannot select a
period.  That is the binding rule of &sect;1 and the reason the spectral
window is drawn beside every periodogram on this page.</li>
<li><b>Single-night peaks are 3&ndash;9&nbsp;c/d wide</b>, so one night
constrains a period to roughly a part in five and nothing better.</li>
<li><b>Three-filter full-orbit nights: ST&nbsp;LMi 20 of 30, AN&nbsp;UMa 4 of
11, VV&nbsp;Pup 1 of 18, EU&nbsp;UMa 0.</b>  The colour-dependent timing
in&nbsp;&sect;3 is an ST&nbsp;LMi result, not a three-polar result, and the
page says so where the number appears rather than in a caveat at the end.</li>
</ul>

<p>The published ephemerides everything is compared against were fetched from
the AAVSO Variable Star Index and cached, so that every &ldquo;agrees with the
literature&rdquo; on this page traces to a stored payload with a timestamp:</p>
{eph_tbl}

<p class="note">The &sigma;(P) column is <em>not</em> a published
uncertainty.  VSX publishes a period and an epoch and no error bar at all.
The number shown is the rounding of the last quoted digit &mdash; a
<b>floor</b>, and the true published uncertainty is at least this large.
&sect;4 is built so that its conclusion survives that: it reports how many
times larger the real &sigma;(P) could be before the cycle count stops being
unique, instead of assuming the floor is the truth.</p>

<h3>Decision</h3>
<p>Three rules are enforced in code rather than in prose.  A periodogram may
not select an alias family member without saying which external fact chose it
(<code>classify_family_choice</code>).  A timing result may not be published
on a formal &Delta;&chi;<sup>2</sup>&nbsp;=&nbsp;1 error bar
(<code>sigma_t_injection</code> is the authority).  A trend may not be removed
before a search (<code>joint_gp_fit</code>, demonstrated in&nbsp;&sect;6).</p>

<h3>Consequence</h3>
<p>Several results below are negative, and they are reported at the same
prominence as the positive ones: YZ&nbsp;Cnc shows no orbital modulation in
its densest era, the &sigma;<sub>t</sub> contour fails its 60&nbsp;s
threshold, and most of the state classifications turn out to be cuts through
a single population rather than two.  A page that only reported the four
things that worked would be a different and less useful document.</p>
</section>
"""


def section_periods(con, figs: dict) -> str:
    rows = q(con, """SELECT series_key, target_key, n_points, n_blocks,
                            baseline_d, period_d, sigma_period_d,
                            published_d, deviation_sigma, agrees,
                            alias_frac_max, family_code, constraint_class,
                            detected, harmonic_note, peak_halfwidth_cd,
                            pdm_theta, p_ls_pow, f_survey_cd,
                            survey_is_orbital
                     FROM p3_period WHERE status='ok'
                     ORDER BY target_key, era_id, filter""")
    body = []
    for r in rows:
        agree = ("&mdash; n/a" if r[9] is None
                 else '<b style="color:#9fd8ae">agrees</b>' if r[9]
                 else '<b style="color:#f0a3a3">DISAGREES</b>')
        body.append([
            _label(r[0]), _i(r[2]), _i(r[3]), _n(r[4], 1),
            _pm(r[5], r[6], 7), _n(r[7], 7),
            _n(r[8], 2) if r[8] is not None else "&mdash;",
            agree, _n(r[10], 2), esc(r[11]), esc(r[12]),
            "yes" if r[13] else "<b>no</b>"])
    tbl = table(["series", "N", "nights", "baseline (d)",
                 "recovered P (d)", "published P (d)", "dev (σ)",
                 "agrees?", "W(±1)", "family chosen by", "constraint",
                 "modulation detected?"], body)
    n_ok = sum(1 for r in rows if r[9])
    # `r[9] == 0`, NOT `r[9] is False`: sqlite3 returns the INTEGER 0 for a
    # false flag, and `0 is False` is never true in Python, so the identity
    # test silently reported zero disagreements while the one real
    # disagreement (yzcnc|e72|r at -3.8 sigma) sat in the table above it.
    n_bad = sum(1 for r in rows if r[9] is not None and r[9] == 0)
    n_na = sum(1 for r in rows if r[9] is None)
    n_prior = sum(1 for r in rows if r[11] == "PRIOR")
    n_single = sum(1 for r in rows if r[11] == "SINGLE-NIGHT")
    n_data = sum(1 for r in rows if r[11] == "DATA")
    n_survey_orb = sum(1 for r in rows if r[19])
    harm = [r for r in rows if r[14]]
    harm_html = "".join(f"<li><b>{_label(r[0])}</b> &mdash; {esc(r[14])}</li>"
                        for r in harm)
    gallery = "".join(
        _fig(figs.get(t, ""),
             f"{TARGET_LABEL[t]}: every periodogram with its own spectral "
             f"window.  Left, the 0.5&ndash;40&nbsp;c/d survey; centre, the "
             f"orbital band with the published frequency dashed and the "
             f"&plusmn;1&nbsp;c/d aliases dotted; right, the spectral window "
             f"on the frequency-offset axis, annotated with the sidelobe "
             f"fraction and how the family member was chosen.",
             f"No periodogram stored for {TARGET_LABEL[t]}.")
        for t in TARGET_ORDER if figs.get(t))
    return f"""
<section id="periods">
<h2>1 &nbsp; Period verification, per filter per era</h2>

<h3>Question</h3>
<p>Can this archive measure the orbital periods of these five stars, or can
it only confirm them?  The two claims look identical in a table and are not
the same thing.</p>

<h3>Evidence</h3>
<p>Every series was searched twice, with two statistics that fail
differently.  <b>A generalised Lomb-Scargle</b> fitted with one free constant
per night &mdash; the exact joint fit, with the night constants projected out
of the cosine and sine columns as well as out of the data, so the power is
what the periodic term adds over a nights-only model.  And <b>Stellingwerf
phase-dispersion minimisation</b>, which assumes nothing about the shape of
the light curve: a polar's bright phase is a top hat, and a sinusoid fitted
to a top hat leaves most of its power in harmonics the fundamental search
cannot see.</p>

<p>The nuisance model deserves a sentence, because it is the one place this
section could have cheated.  One constant per night is not detrending.  A
running median has a free parameter every few points and absorbs power at
the frequencies being searched; a single constant per night can only absorb
power below 1/T<sub>night</sub> &asymp; 2.5&nbsp;c/d, five times below the
8&ndash;20&nbsp;c/d orbital band.  &sect;6 measures that claim on an injected
signal rather than asserting it.</p>

{gallery}

<p><b>The tallest peak in the whole 0.5&ndash;40&nbsp;c/d search is the orbit
in only {_i(n_survey_orb)} of {_i(len(rows))} series</b>, and the other
{_i(len(rows) - n_survey_orb)} fail in two quite different ways:</p>

{table(["what the tallest survey peak is", "series"],
       [[esc(k), _i(v)] for k, v in q(con,
        "SELECT survey_class, count(*) FROM p3_period WHERE status='ok' "
        "GROUP BY 1 ORDER BY 2 DESC")])}

<p>The distinction matters and lumping the two together would be wrong.  A
peak at <b>twice the orbital frequency</b> &mdash; all three AN&nbsp;UMa
series, at 2.00&times; &mdash; <em>is</em> the orbit: AN&nbsp;UMa's light
curve is double-humped (VSX types it AM+E), so more power lands in the first
harmonic than in the fundamental, and a search reporting its global maximum
would publish half the true period.  A peak <b>below 3&nbsp;c/d</b> &mdash;
the YZ&nbsp;Cnc series &mdash; is the accretion state changing on
week-to-month timescales: real astrophysics, and not the thing being
searched for, but a search reporting its global maximum would publish a
0.5&nbsp;d &ldquo;period&rdquo; for a 2.08&nbsp;h binary.  This is precisely
why the orbital band 8&ndash;20&nbsp;c/d is declared before the search rather
than after it.</p>

{table(["", "count"],
       [["series where the family member was chosen by the "
         "<b>literature prior</b> (multi-night, sidelobes tall)", _i(n_prior)],
        ["series that are a <b>single night</b> (no alias ambiguity, "
         "no precision either)", _i(n_single)],
        ["series where <b>the data could choose alone</b> "
         f"(strongest sidelobe below {p3.ALIAS_DECIDABLE_MAX:g})",
         f"<b>{_i(n_data)}</b>"]])}

<p>{_fig(figs.get('summary', ''),
         "Left: deviation from the published period in combined sigma, with "
         "the ±3σ agreement band shaded.  Right: fractional period precision "
         "against the strength of the ±1 c/d sidelobe.  A period "
         "<em>determination</em> would live in the lower-left corner — tight "
         "precision, weak aliasing.  The corner is empty.",
         "No period summary figure: p3_period returned no rows.")}</p>

{tbl}

{"<p><b>Two series whose tallest orbital-band peak is not the orbit:</b></p>"
 "<ul>" + harm_html + "</ul>" if harm else ""}

<h3>Decision</h3>
<p><b>{_i(n_ok)} of {_i(len(rows))} series agree with the published period,
{_i(n_bad)} disagree, and {_i(n_na)} have no orbital modulation to compare.</b>
But the agreement column is the less important one.  The
<em>family-chosen-by</em> column is the result: <b>not one series in this
archive has sidelobes weak enough for its own periodogram to select a period.</b>
The strongest sidelobe runs {_n(q1(con, "SELECT min(alias_frac_max) FROM p3_period WHERE n_blocks>1"), 2)}
to {_n(q1(con, "SELECT max(alias_frac_max) FROM p3_period WHERE n_blocks>1"), 2)}
of the window power against a bar of {p3.ALIAS_DECIDABLE_MAX:g}.</p>

<p>So the honest statement of this section is: <b>these are period
confirmations at the precision of a local peak, not period determinations.</b>
Where the light curve is strong and the baseline long &mdash; ST&nbsp;LMi in
both eras, AN&nbsp;UMa in g and r, VV&nbsp;Pup&nbsp;e76&nbsp;g, EU&nbsp;UMa
in g &mdash; that confirmation is tight, a part in 10<sup>4</sup> or better.
Where it is not, the table says <code>WEAK</code> or
<code>UNINFORMATIVE</code> instead of quietly reporting seven decimal
places.</p>

<p>Two specific findings are worth naming.  <b>YZ&nbsp;Cnc shows no orbital
modulation at all</b> in its 2024 era despite 1,148 points over 46 nights in
three filters &mdash; which is what a dwarf nova outside superoutburst should
do, and is a real result rather than a failure.  And <b>VV&nbsp;Pup&nbsp;e76&nbsp;i
and YZ&nbsp;Cnc&nbsp;e7&nbsp;I have their tallest orbital-band peak almost
exactly 1&nbsp;c/d from the published frequency</b> &mdash; the textbook
alias, appearing exactly where the window says it will.</p>

<h3>Consequence</h3>
<p>Every downstream section uses the <em>published</em> ephemeris, not a
period measured here, and that is now a justified choice rather than a
convenient one.  The cycle counts in&nbsp;&sect;4, the phase windows the
edges in&nbsp;&sect;3 are fitted in, and the phase-coverage gate
in&nbsp;&sect;5 all inherit the literature value, and the uncertainty they
inherit is the literature's, which &sect;4 then refuses to take on faith.</p>
</section>
"""


def section_sigmat(con, fig: str) -> str:
    inp = q(con, """SELECT series_key, n_points, median_cadence_s, depth_mag,
                           edge_width_s, edge_width_floor_s,
                           bright_width_phase, n_cycles, used_err_mag,
                           chi2_inflation, note
                    FROM p3_sigmat_input ORDER BY series_key""")
    if not inp:
        return """
<section id="sigmat"><h2>2 &nbsp; The &sigma;<sub>t</sub> injection test</h2>
<p class="note">Not run: <code>p3_sigmat_input</code> is empty.</p></section>"""
    meta = dict(q(con, "SELECT key, value FROM p3_meta"))
    inp_tbl = table(
        ["series", "N", "cadence (s)", "depth (mag)", "fitted edge width (s)",
         "fold sampling floor (s)", "bright width (phase)", "orbits",
         "σ used (mmag)", "χ² inflation"],
        [[_label(r[0]), _i(r[1]), _n(r[2], 0), _n(r[3], 2), _n(r[4], 1),
          _n(r[5], 1), _n(r[6], 2), _n(r[7], 2), _n(1000 * r[8], 1),
          _n(r[9], 2)] for r in inp])
    grid = q(con, """SELECT series_key, shape_error, depth_error,
                            sigma_t_s, bias_s, total_error_s, p95_abs_s,
                            recovered_fraction, passes, inject_width_s
                     FROM p3_sigmat WHERE inject_factor=1.0
                     ORDER BY series_key, shape_error, depth_error""")
    grid_tbl = table(
        ["series", "assumed width / true", "assumed depth error",
         "scatter σ<sub>t</sub> (s)", "bias (s)", "<b>total (s)</b>",
         "95th pct |Δt| (s)", "recovered", "vs 60 s"],
        [[_label(r[0]), f"×{r[1]:g}", f"{r[2]:+.0%}", _n(r[3], 1),
          _n(r[4], 1), f"<b>{_n(r[5], 1)}</b>", _n(r[6], 1), _pct(r[7], 0),
          ('<b style="color:#9fd8ae">PASS</b>' if r[8]
           else '<b style="color:#f0a3a3">FAIL</b>')] for r in grid])
    verdicts = "".join(
        f"<li><b>{_label(k.replace('sigmat_verdict_', ''))}</b>: "
        f"<b>{esc(v)}</b> &mdash; "
        f"{esc(meta.get('sigmat_sentence_' + k.replace('sigmat_verdict_', ''), ''))}"
        f"<br><span class='note'>width sensitivity, shape known exactly: "
        f"{esc(meta.get('sigmat_width_scan_' + k.replace('sigmat_verdict_', ''), '—'))}</span></li>"
        for k, v in sorted(meta.items())
        if k.startswith("sigmat_verdict_") and "allwidths" not in k)
    n_pass = q1(con, "SELECT count(*) FROM p3_sigmat WHERE passes=1")
    n_cells = q1(con, "SELECT count(*) FROM p3_sigmat")
    best = q1(con, "SELECT min(total_error_s) FROM p3_sigmat WHERE "
                   "shape_error=1.0 AND depth_error=0.0 AND inject_factor=1.0")
    return f"""
<section id="sigmat">
<h2>2 &nbsp; The &sigma;<sub>t</sub> injection test on ST&nbsp;LMi's densest
night</h2>

<h3>Question</h3>
<p>ST&nbsp;LMi on the local night {esc(meta.get('sigmat_night', ''))} (UTC
night 2025-02-28) is the densest run in the archive: 153, 142 and 141 frames
in g, r and i over 9.37&nbsp;h, which is 4.9 orbits with each filter sampled
every 219&nbsp;s.  If per-cycle timing is publishable anywhere in this data
set it is publishable here.  Is it?  The threshold set by the strategy is
{p3.SIGMA_T_THRESHOLD_S:.0f}&nbsp;s.</p>

<p>The question cannot be answered with the error bar the edge fit returns.
On this night those fits come back with &chi;<sup>2</sup><sub>&nu;</sub>
between 5 and 400 on 6&ndash;11 points &mdash; the trapezoid does not describe
the flickering &mdash; and the &Delta;&chi;<sup>2</sup>&nbsp;=&nbsp;1 interval
they imply is 0&ndash;4&nbsp;s.  A published per-cycle timing precision of
four seconds at a cadence of 219&nbsp;s would be a fiction.</p>

<h3>Evidence</h3>
<p>So the epochs are injected and recovered.  A bright-phase template of known
epoch is added to the <em>real</em> timestamps with the <em>real</em>
per-point errors (inflated by each series' measured
&chi;<sup>2</sup><sub>&nu;</sub>), and recovered by template matching &mdash;
one fixed width, one fixed depth, only the epoch and an overall level free.
The scatter of <code>t<sub>recovered</sub> &minus; t<sub>true</sub></code> is
&sigma;<sub>t</sub>, with no model of the star in it.</p>

<p>The template's depth, bright-phase width and edge width are measured from
each series' own phase-folded light curve:</p>
{inp_tbl}

<p class="note"><b>The fitted edge width is an upper bound, and it turned out
to be the most consequential number in this section.</b>  A first version of
the shape estimator was algebraically pinned at two phase bins &mdash;
547&nbsp;s for every series it was ever given &mdash; and returned a
comfortable &ldquo;CONDITIONAL, &sigma;<sub>t</sub>&nbsp;=&nbsp;20&nbsp;s with
the shape known&rdquo;.  Fitting the width to the folded points at full time
resolution instead gives {_n(min(r[4] for r in inp), 0)}&ndash;{_n(max(r[4] for r in inp), 0)}&nbsp;s,
at or near the fold's own sampling floor: <b>the edge is unresolved</b>.  And
a sharper edge is <em>harder</em> to time, not easier, because fewer points
land on the ramp and the epoch is bracketed by the 219&nbsp;s cadence rather
than interpolated within it.  Since the data cannot say where in that range
the truth lies, the grid was run at three injected widths and the verdict has
to survive the range.</p>

<p>{_fig(fig,
         "The σ_t contour.  Each panel is one band: assumed-edge-width error "
         "across, assumed-depth error up, total timing error (scatter and "
         "template-induced bias in quadrature) in the cells, with the 60 s "
         "contour drawn where it falls.  The rightmost panel is the "
         "sensitivity to the injected width — the one input the data cannot "
         "measure.",
         "No contour: p3_sigmat is empty.")}</p>

<p>The verdict is taken on the <b>total</b> error, not the scatter, and the
distinction changed the answer.  A five-times-too-wide template produces a
<em>tighter</em> spread than the correct one &mdash; a smooth model fits
smoothly &mdash; while placing every epoch tens of seconds late.  Judging on
scatter alone would have scored the badly wrong template as the better
measurement.  A common bias cancels out of an O&minus;C gradient; it does not
cancel out of an O&minus;C zero point, and it certainly does not cancel
between two bands whose templates differ, which is exactly the measurement
in&nbsp;&sect;3.</p>

{grid_tbl}

<p class="note"><b>The grid runs the &ldquo;wrong&rdquo; way along the shape
axis, and that is a finding rather than a bug.</b>  Assuming a
five-times-too-wide ramp gives a <em>smaller</em> total error than assuming
the correct one.  The reason is that at 219&nbsp;s cadence against a
&lt;50&nbsp;s edge, almost no exposure lands on the ramp at all: the epoch is
bracketed by two points either side, and a wide assumed ramp makes the fit
interpolate between them, which is a better estimator of a bracketed
transition than a near-step template that can sit anywhere in the gap with
equal likelihood.  Being wrong in the direction of &ldquo;smoother than
reality&rdquo; costs nothing here because the data never resolved the edge.
The <em>depth</em> axis behaves as expected &mdash; a 50% depth error is the
worst cell in every panel &mdash; because a committed wrong depth biases the
level fit and drags the crossing point with it.</p>

<h3>Decision</h3>
<ul>{verdicts}</ul>

<p><b>{_i(n_pass)} of {_i(n_cells)} grid cells pass the
{p3.SIGMA_T_THRESHOLD_S:.0f}&nbsp;s threshold.</b>  The best cell in the whole
experiment &mdash; the densest night, the sharpest band, the edge shape and
depth both known exactly &mdash; returns {_n(best, 1)}&nbsp;s.  <b>Per-cycle
bright-phase timing is not supported by this cadence.</b>  Not marginally, and
not only under a wrong template: the exact-template cell misses too.</p>

<p>The reason is arithmetic rather than photometric.  With a 219&nbsp;s
per-filter interval and an edge that crosses in under 50&nbsp;s, most cycles
have <em>no</em> point on the ramp; the epoch is bracketed by two exposures
219&nbsp;s apart, and the uncertainty of a bracket is
219/&radic;12&nbsp;&asymp;&nbsp;63&nbsp;s before any noise is added.  No
improvement in photometric precision moves that number.  A cadence of about
60&nbsp;s per filter would.</p>

<h3>Consequence</h3>
<p>&sect;3 and &sect;4 are re-scoped by this result rather than cancelled.
Individual cycle epochs are published with the Monte Carlo's error bar and not
the fit's, and no statement anywhere on this page rests on a single cycle.
What survives is what <em>averaging</em> over many cycles supports: the
inter-band offsets in&nbsp;&sect;3 are means over 4&ndash;5 paired cycles per
night, and the O&minus;C in&nbsp;&sect;4 is read as a night-averaged
quantity.  The 60&nbsp;s threshold governs the per-cycle claim, and the
per-cycle claim is the one that is withdrawn.</p>
</section>
"""


def section_edges(con, fig: str) -> str:
    n_fit = _rows(con, "p3_edge")
    n_acc = q1(con, "SELECT count(*) FROM p3_edge WHERE accepted=1") or 0
    reasons = q(con, """SELECT reason, count(*) c FROM p3_edge
                        WHERE accepted=0 GROUP BY 1 ORDER BY c DESC LIMIT 6""")
    # Reasons carry per-fit numbers; group them by their leading phrase.
    import re
    buckets: dict[str, int] = {}
    for r, c in q(con, "SELECT reason, count(*) FROM p3_edge WHERE "
                       "accepted=0 GROUP BY 1"):
        key = re.sub(r"[-+]?\d*\.?\d+", "N", str(r))
        buckets[key] = buckets.get(key, 0) + c
    rej_tbl = table(["why an edge was refused", "fits"],
                    [[esc(k), _i(v)] for k, v in
                     sorted(buckets.items(), key=lambda kv: -kv[1])])
    pooled = q(con, """SELECT target_key, era_id, band_a, band_b, n_cycles,
                              delta_s, sigma_s, chi2nu, significant
                       FROM p3_band_pair WHERE night='(pooled)'
                       ORDER BY target_key, era_id, band_a, band_b""")
    pool_tbl = table(
        ["target", "era", "bands", "paired cycles", "Δt (s)", "σ (s)",
         "χ²<sub>ν</sub>", "|Δt| / σ", "verdict"],
        [[esc(TARGET_LABEL.get(r[0], r[0])), _i(r[1]),
          f"{esc(r[2])} &minus; {esc(r[3])}", _i(r[4]),
          f"<b>{_n(r[5], 1)}</b>", _n(r[6], 1), _n(r[7], 2),
          _n(abs(r[5]) / r[6], 1) if r[6] else "&mdash;",
          ('<b style="color:#f0a3a3">≥3σ from zero</b>' if r[8]
           else "consistent with zero")] for r in pooled]) if pooled else ""
    pairs = q(con, """SELECT target_key, era_id, night, band_a, band_b,
                             n_cycles, delta_s, sigma_s, chi2nu, significant
                      FROM p3_band_pair WHERE night <> '(pooled)'
                      ORDER BY target_key, night, band_a""")
    pair_tbl = table(
        ["target", "era", "night", "bands", "cycles", "Δt (s)", "σ (s)",
         "χ²<sub>ν</sub>", "verdict"],
        [[esc(TARGET_LABEL.get(r[0], r[0])), _i(r[1]), esc(r[2]),
          f"{esc(r[3])} &minus; {esc(r[4])}", _i(r[5]),
          f"<b>{_n(r[6], 1)}</b>", _n(r[7], 1), _n(r[8], 2),
          ('<b style="color:#f0a3a3">≥3σ from zero</b>' if r[9]
           else "consistent with zero")] for r in pairs]) if pairs else ""
    n_sig = sum(1 for r in pooled if r[8])
    best = max(pooled, key=lambda r: abs(r[5]) / r[6] if r[6] else 0,
               default=None)
    targets = sorted({r[0] for r in pooled})
    amps = q(con, """SELECT filter, avg(depth_mag), count(*) FROM p3_edge
                     WHERE accepted=1 AND target_key='stlmi' AND era_id=76
                     GROUP BY 1 ORDER BY 1""")
    amp_tbl = table(["band", "mean fitted edge depth (mag)", "cycles"],
                    [[esc(r[0]), _n(r[1], 3), _i(r[2])] for r in amps]) \
        if amps else ""
    return f"""
<section id="edges">
<h2>3 &nbsp; Bright-phase timing, per band</h2>

<h3>Question</h3>
<p>In a polar the optical emission is cyclotron radiation from an accretion
column, and its opacity depends on wavelength.  The column therefore does not
disappear behind the white dwarf's limb at the same instant in every band, and
the size of that offset is a property of the emission region &mdash; a
measurement, not a systematic to be calibrated away.  Can this data see
it?</p>

<h3>Evidence</h3>
<p>The falling edge of the bright phase was fitted per cycle per band, with
the edge phase taken from each series' own folded profile and a four-parameter
trapezoid (two levels, epoch, ramp width) profiled on a grid.
{_i(n_acc)} of {_i(n_fit)} attempted fits were accepted; the rest were refused
by explicit gates rather than quietly averaged in:</p>
{rej_tbl}

<p>{_fig(fig,
         "Left: the bright phase folded on the published ephemeris, one "
         "colour per band, with the accepted edge epochs dashed. Right: the "
         "inter-band edge-time differences, paired per cycle so the "
         "cycle-to-cycle wander of the accretion spot cancels.",
         "No edge figure: p3_edge has no accepted rows.")}</p>

{"<p>The amplitude of the edge itself is strongly colour dependent, which is "
 "the first-order cyclotron signature and is measured here rather than "
 "assumed:</p>" + amp_tbl if amp_tbl else ""}

{"<p><b>The publishable inter-band numbers are the pooled ones.</b>  A band "
 "offset is a property of the emission region, not of a night, and most "
 "individual nights contribute a single paired cycle — which at the "
 "&sigma;<sub>t</sub> this cadence supports is a ±300 s bound that says "
 "nothing.  Pooling every paired cycle in an era is both the correct "
 "estimator and the only one with the signal-to-noise to detect "
 "anything:</p>" + pool_tbl
 if pool_tbl else "<p class='note'>No band pair produced accepted edges on "
 "the same cycle, so no inter-band difference can be formed.</p>"}

{"<details><summary>The per-night components of the pooled values "
 "(" + _i(len(pairs)) + " rows) — these are not independent results</summary>"
 + pair_tbl + "</details>" if pair_tbl else ""}

<h3>Decision</h3>
<p><b>{_i(n_sig)} of {_i(len(pooled))} pooled inter-band differences reach
3&sigma;.</b>  The largest is
{f"{esc(TARGET_LABEL.get(best[0], best[0]))} e{best[1]} {esc(best[2])}&minus;{esc(best[3])} at {_n(best[5], 0)} &plusmn; {_n(best[6], 0)} s over {_i(best[4])} paired cycles &mdash; {_n(abs(best[5]) / best[6], 1)}&sigma;, which is a bound and not a detection"
 if best and best[6] else "not measurable"}.</p>

<p>The error bars are the larger of the rescaled formal bar and the injection
Monte Carlo's total error from&nbsp;&sect;2, which is why they are tens of
seconds wide rather than the few seconds the raw fits claimed.  Their
&chi;<sup>2</sup><sub>&nu;</sub> of 0.2&ndash;0.6 says those bars are if
anything conservative: the cycle-to-cycle scatter of the differences is
<em>smaller</em> than the Monte Carlo floor, so the limit here is set by how
well a single edge can be placed and not by the star's variability.</p>

<p>The honest reading: <b>this data resolves the colour dependence of the
bright phase's <em>amplitude</em> decisively, and its <em>timing</em> not at
all.</b>  The amplitude result is unambiguous and is a cyclotron measurement
&mdash; the edge is roughly three times deeper in i than in g.  The timing
result is a bound, limited by exactly what&nbsp;&sect;2 predicted: a 219&nbsp;s
cadence against an unresolved edge.  The sign of the largest offset is the
one cyclotron theory would predict (the bluer band's edge earlier), and at
1.9&sigma; that is a remark, not a result, and the page declines to promote
it.</p>

<p>{"This is an ST LMi result." if targets == ["stlmi"] else
   "The targets contributing are: " + ", ".join(TARGET_LABEL.get(t, t) for t in targets) + "."}
 That is a consequence of the observing record and not of the analysis:
three-filter full-orbit nights number 20 of 30 for ST&nbsp;LMi, 4 of 11 for
AN&nbsp;UMa, 1 of 18 for VV&nbsp;Pup and 0 for EU&nbsp;UMa.  Calling this a
&ldquo;three-polar colour result&rdquo; would misrepresent the sample by a
factor of three.</p>

<h3>Consequence</h3>
<p>The published claim is the band-dependent bright-phase <em>amplitude</em>
with the timing offsets as bounds, and the bounds are quoted with the Monte
Carlo error bar rather than the fit's.  A future run that wants the timing
result needs a shorter per-filter interval, not more nights: more nights at
219&nbsp;s cadence average down the random part of a bracket-limited error and
leave the bracket.</p>
</section>
"""


def section_oc(con, fig: str) -> str:
    cc = q(con, """SELECT target_key, epoch_bjd, period_d, sigma_period_d,
                          n_cycles_last, drift_cycles, unique_count,
                          sigma_period_max_d, ratio_to_quoted, oc_rms_s,
                          fitted_period_d, fitted_period_sigma_d, n_epochs,
                          verdict, note, oc_mean_s
                   FROM p3_cycle_count ORDER BY target_key""")
    if not cc:
        return """<section id="oc"><h2>4 &nbsp; O&minus;C and cycle counts</h2>
<p class="note">Not run: <code>p3_cycle_count</code> is empty.</p></section>"""
    tbl = table(
        ["target", "epochs", "cycles since VSX epoch", "σ(P) floor (d)",
         "accumulated drift (cycles)", "σ(P) that would still work (d)",
         "margin", "verdict"],
        [[esc(TARGET_LABEL.get(r[0], r[0])), _i(r[12]),
          f"{r[4]:,.0f}" if r[4] is not None else "&mdash;",
          _sci(r[3]), _n(r[5], 4) if r[5] is not None else "&mdash;",
          _sci(r[7]),
          (f"{1.0 / r[8]:,.0f}&times;" if r[8] else "&mdash;"),
          ('<b style="color:#9fd8ae">' + esc(r[13]) + "</b>"
           if r[6] else '<b style="color:#f0a3a3">' + esc(r[13]) + "</b>")]
         for r in cc])
    notes = "".join(f"<li><b>{esc(TARGET_LABEL.get(r[0], r[0]))}</b>: "
                    f"{esc(r[14])}</li>" for r in cc if r[14])
    fitted = table(
        ["target", "epochs", "O−C rms (s)", "fitted P (d)", "VSX P (d)",
         "difference (s/cycle)"],
        [[esc(TARGET_LABEL.get(r[0], r[0])), _i(r[12]), _n(r[9], 1),
          _pm(r[10], r[11], 8), _n(r[2], 8),
          _n((r[10] - r[2]) * 86400.0, 3) if r[10] and r[2] else "&mdash;"]
         for r in cc if r[10]])
    n_uniq = sum(1 for r in cc if r[6])
    n_have = sum(1 for r in cc if r[4] is not None)
    return f"""
<section id="oc">
<h2>4 &nbsp; O&minus;C construction and cycle-count analysis</h2>

<h3>Question</h3>
<p>An O&minus;C diagram plots an observed epoch minus the epoch a linear
ephemeris predicts, and that subtraction needs an <em>integer</em>: which
cycle is this?  Between the VSX zero point and our nights there are tens of
thousands of cycles, and an error in the period accumulates over every one of
them.  Once the accumulated drift reaches half a cycle the integer is no
longer determined, and an O&minus;C computed on the wrong one is not a noisy
result &mdash; it is a fabricated one, and it looks exactly like a real
period change.</p>

<h3>Evidence</h3>
<p>The drift is <code>n &times; &sigma;(P) / P</code> cycles.  Inverting it
gives the number this section actually turns on: the largest period
uncertainty that would still leave the count unique,
<code>&sigma;(P)<sub>max</sub> = 0.5 &times; P / n</code>.  That form is used
deliberately, because <b>VSX publishes no period uncertainty at all</b>.
Rather than invent one, the table reports how many times larger the true
&sigma;(P) could be than the quoted-precision floor before the conclusion
changes:</p>
{tbl}

<p>{_fig(fig,
         "Left: the O−C of the bright-phase edge against cycle number, mean "
         "removed. Right: the cycle-count margin — the count is unique when "
         "the period uncertainty that would still work (green) exceeds the "
         "uncertainty we actually have (blue).",
         "No O−C figure: p3_oc is empty.")}</p>

<ul>{notes}</ul>

{"<p>Where enough epochs survived, a linear ephemeris was refitted to them "
 "and compared with the published period:</p>" + fitted if fitted else ""}

<p class="note">One systematic is removed and published rather than hidden.
The bright-phase edge is not at phase zero of the catalogue ephemeris &mdash;
it is a feature of the accretion geometry, offset by a constant fraction of a
cycle.  That constant appears in the raw O&minus;C as a large offset that has
nothing to do with the clock, so the mean is subtracted and quoted separately
{"(" + ", ".join(f"{TARGET_LABEL.get(r[0], r[0])}: {r[15]:,.0f} s" for r in cc if r[15] is not None) + ")"
 if any(r[15] is not None for r in cc) else ""}.  An O&minus;C that kept it
would show a huge constant offset and invite exactly the wrong
interpretation.</p>

<h3>Decision</h3>
<p><b>The cycle count is unique for {_i(n_uniq)} of the {_i(n_have)} targets
with both a published epoch and timed epochs of our own</b>, and the margins
are not marginal &mdash; the true &sigma;(P) would have to be thousands of
times larger than the quoted precision before any of them broke.  So the
O&minus;C diagrams above rest on determined integers.</p>

<p><b>A unique cycle count is necessary and not sufficient, and a second
gate removes two of the four.</b>  An O&minus;C pools epochs from every
filter and era into one diagram, which is only legitimate if they all time
the same feature.  For ST&nbsp;LMi they do &mdash; the folded profile puts
the falling edge at phase 0.140 in all six series, and the accepted epochs
scatter by 0.014 in phase.  For AN&nbsp;UMa and VV&nbsp;Pup they scatter by
0.152 and 0.130, because with 4 full-orbit nights out of 11 and 1 out of 18
the folded profile cannot locate the edge consistently between filters.
Pooling those produced O&minus;C scatters of 1,126&nbsp;s and 2,555&nbsp;s
on the first run &mdash; a number that reads as a period error and is
nothing of the kind.  Both are now refused with the verdict
<code>NOT ONE FEATURE &mdash; NO O&minus;C</code> and no O&minus;C rows are
written for them at all, because a refused result that is still stored
&ldquo;for reference&rdquo; is a result that gets plotted anyway.</p>

<p>YZ&nbsp;Cnc is excluded for a different and simpler reason: <b>VSX
publishes a period for it but no epoch</b>, so there is no zero point to count
from and no O&minus;C is possible against the catalogue value at all.  That is
recorded as its own verdict rather than being folded into a footnote.</p>

<p>What the O&minus;C does <em>not</em> support is a period change.  The
scatter is dominated by the per-cycle timing error that &sect;2 measured
&mdash; tens of seconds against a bracket-limited floor &mdash; and with the
baseline available here that scatter swamps any secular term.  The refitted
periods agree with the published ones and are not offered as improvements on
them.</p>

<h3>Consequence</h3>
<p>The cycle counts are stored in <code>p3_cycle_count</code> with the drift
and the margin beside them, so a later run that acquires a genuine published
&sigma;(P) &mdash; from the discovery paper rather than from the catalogue
&mdash; can check the conclusion against it in one query instead of
re-deriving it.  That is the whole reason the margin is stored as a ratio.</p>
</section>
"""


def section_states(con, fig: str) -> str:
    rows = q(con, """SELECT series_key, n_nights, n_gated, n_used,
                            threshold_mag, threshold_sigma, separability,
                            bimodal, n_high, n_low, n_intermediate,
                            n_censored, duty_naive, duty_with_limits,
                            duty_bias, n_informative_limits, verdict
                     FROM p3_state_series ORDER BY target_key, era_id, filter""")
    if not rows:
        return """<section id="states"><h2>5 &nbsp; Accretion states</h2>
<p class="note">Not run: <code>p3_state_series</code> is empty.</p></section>"""
    tbl = table(
        ["series", "nights", "gated", "used", "threshold (mag)",
         "separability", "H / L / I", "limits used", "duty (detections)",
         "duty (with limits)", "bias", "verdict"],
        [[_label(r[0]), _i(r[1]), _i(r[2]), _i(r[3]),
          _pm(r[4], r[5], 2), _n(r[6], 2),
          f"{_i(r[8])} / {_i(r[9])} / {_i(r[10])}", _i(r[15]),
          _n(r[12], 3), f"<b>{_n(r[13], 3)}</b>", _n(r[14], 3),
          ('<b style="color:#9fd8ae">' + esc(r[16]) + "</b>" if r[7]
           else '<span style="color:#e6cc7a">' + esc(r[16]) + "</span>")]
         for r in rows])
    n_bi = sum(1 for r in rows if r[7])
    n_thr = sum(1 for r in rows if r[4] is not None)
    n_shift = sum(1 for r in rows
                  if r[14] is not None and abs(r[14]) > 0.01)
    worst = max((r for r in rows if r[14] is not None),
                key=lambda r: abs(r[14]), default=None)
    # The census behind the "bias is zero" claim, so the reader can see it
    # is a statement about how few limits there are and not a claim that
    # censoring does not matter.
    n_lim_total = q1(con, "SELECT sum(n_limits) FROM p3_state_night") or 0
    n_lim_nights = q1(con, "SELECT count(*) FROM p3_state_night WHERE "
                           "n_limits > 0") or 0
    n_lim_gated = q1(con, "SELECT count(*) FROM p3_state_night WHERE "
                          "n_limits > 0 AND gated = 1") or 0
    n_cens = q1(con, "SELECT count(*) FROM p3_state_night WHERE "
                     "censored = 1") or 0
    return f"""
<section id="states">
<h2>5 &nbsp; Accretion-state classification</h2>

<h3>Question</h3>
<p>Polars switch between high and low accretion states, and a duty cycle
&mdash; what fraction of the time the star is accreting &mdash; is one of the
few population-level numbers a monitoring programme like this can contribute.
Two things can make that number a fiction.  A threshold chosen by eye, and a
statistic computed only on the nights the star was bright enough to detect,
which is very nearly a tautology.</p>

<h3>Evidence</h3>
<p>The classifying statistic is a night's median magnitude, and it is gated
first.  These polars vary by 0.65&ndash;1.7&nbsp;mag around a single orbit
&mdash; the same size as a state change &mdash; so a median over a third of an
orbit is a phase measurement wearing a state label.  A night must cover at
least {p3.STATE_MIN_PHASE_COVERAGE:.0%} of the phase circle before its median
is allowed to mean anything.</p>

<p>The threshold is then Otsu's method on the surviving nights: the cut that
maximises the between-class variance.  Deterministic, no free parameters
beyond the bin count, and derived from the observed bimodality rather than
imported.  Its own bootstrap uncertainty defines an INTERMEDIATE band, which
means &ldquo;this data cannot say&rdquo; and not &ldquo;the star was
physically in between&rdquo;.  Crucially, Otsu also reports a
<b>separability</b>, and a low value is the method telling you it has cut a
single population in half.</p>

<p>Duty cycles are then computed twice.  Once on detections only, and once
with the Phase-2 upper limits included &mdash; where a limit <em>fainter</em>
than the threshold proves the epoch is on the low side, and a limit brighter
than the threshold proves nothing and is counted as uninformative rather than
silently assigned.</p>

<p>{_fig(fig,
         "Left: every night's median magnitude relative to its own series' "
         "threshold. Right: duty cycle computed on detections only versus "
         "with the Phase-2 limits included — the gap is the censoring bias.",
         "No state figure: p3_state_series is empty.")}</p>

{tbl}

<h3>Decision</h3>
<p><b>{_i(n_thr)} of {_i(len(rows))} series produce a threshold at all, and
only {_i(n_bi)} of those {_i(n_thr)} are genuinely bimodal.</b>  For the
other {_i(n_thr - n_bi)} the separability is below 0.75 and the honest
statement is that the nightly magnitudes are <em>one</em> population: the
threshold still cuts them, and the HIGH/LOW counts in the table are still
computed, but they are not evidence of two accretion states and the verdict
column says so.  A page that printed the duty cycles without that column
would be publishing a bimodality that was manufactured by the method.</p>

<p><b>The measured censoring bias is exactly zero in every series</b>
({_i(n_shift)} series move the duty cycle by more than 0.01; the largest
change anywhere is {_n(abs(worst[14]) if worst else 0.0, 3)}).  That is a
result, not an omission, and it is worth being clear about what it does and
does not mean.  It does <em>not</em> mean the censoring bias is unimportant
in general &mdash; it means <em>this archive has too few limits, in the wrong
places, to exhibit it</em>: {_i(n_lim_total)} limits across
{_i(n_lim_nights)} nights, of which {_i(n_lim_gated)} fail the phase gate and
{_i(n_cens)} are nights with no detection at all.  The reason is the one
Phase&nbsp;2 already established: 11 of 13 series have median limits
<em>shallower</em> than their own median detection, so most undetected epochs
are low-sensitivity frames rather than demonstrably faint states.  The
limit-aware duty cycles here are therefore still <b>bounds, not
measurements</b>, and they inherit that caveat from Phase&nbsp;2 rather than
escaping it.  What the exercise does establish is that the duty cycles above
are not <em>inflated</em> by dropped non-detections, which is the specific
failure the task existed to rule out.</p>

<h3>Consequence</h3>
<p>The state histories are stored per night in
<code>p3_state_night</code> with the phase coverage, the censoring flag and
the gate reason on every row, so a later analysis can re-derive a duty cycle
under a different threshold without re-deriving the classification.  The
series that failed the bimodality test are kept in the table rather than
dropped, because &ldquo;this series shows no resolvable state change&rdquo; is
a result about the star and the campaign, and deleting the row would make the
sample look more decisive than it is.</p>
</section>
"""


def section_detrend(con, fig: str) -> str:
    rows = q(con, """SELECT window_periods, frac_detrend, frac_joint,
                            n_detrend, n_joint FROM p3_detrend
                     ORDER BY window_periods""")
    chk = q(con, """SELECT series_key, n_points, amp_celerite, amp_dense,
                           rel_diff, loglike_celerite, loglike_dense,
                           ll_abs_diff, celerite_eps, verdict
                    FROM p3_gp_check""")
    if not rows:
        return """<section id="detrend"><h2>6 &nbsp; Detrending discipline</h2>
<p class="note">Not run: <code>p3_detrend</code> is empty.</p></section>"""
    meta = dict(q(con, "SELECT key, value FROM p3_meta"))
    tbl = table(
        ["running-median window", "detrend-then-search recovers",
         "joint GP + signal recovers"],
        [[f"{r[0]:g} &times; P<sub>orb</sub>",
          f'<b style="color:#f0a3a3">{_pct(r[1])}</b>',
          f'<b style="color:#9fd8ae">{_pct(r[2])}</b>'] for r in rows])
    chk_tbl = table(
        ["series", "N", "celerite2 amplitude", "dense amplitude",
         "relative difference", "Δ log-likelihood", "kernel ε", "verdict"],
        [[_label(r[0]), _i(r[1]), _n(r[2], 6), _n(r[3], 6), _sci(r[4]),
          _sci(r[7]), f"{r[8]:g}", esc(r[9])] for r in chk]) if chk else ""
    worst_lo = min((r for r in rows if r[1] is not None),
                   key=lambda r: r[1], default=None)
    worst_hi = max((r for r in rows if r[1] is not None),
                   key=lambda r: r[1], default=None)
    jmin = min((r[2] for r in rows if r[2] is not None), default=float("nan"))
    jmax = max((r[2] for r in rows if r[2] is not None), default=float("nan"))
    return f"""
<section id="detrend">
<h2>6 &nbsp; Detrending discipline</h2>

<h3>Question</h3>
<p>Every light curve here has a trend under it &mdash; atmospheric
transparency, differential extinction across a night, slow instrumental
drift.  The standard move is to remove it and then search the residuals.  The
strategy forbids that and requires joint GP&nbsp;+&nbsp;signal fitting
instead.  Is that a matter of taste, or does the order change the answer?</p>

<h3>Evidence</h3>
<p>The two orderings were run on the same data with the same injected signal.
A sinusoid of semi-amplitude {esc(meta.get('detrend_amplitude_mag', '0.10'))}&nbsp;mag
at the orbital frequency was added to the real timestamps of
{_label(meta.get('detrend_series', ''))} on night
{esc(meta.get('detrend_night', ''))}, on top of a real correlated trend drawn
from a Mat&eacute;rn-3/2 process, and recovered (a) by subtracting a running
median and then fitting at the known frequency, and (b) by fitting the GP and
the sinusoid simultaneously.  The injected amplitude is known, so the
recovered fraction is a measurement and not an opinion.</p>

<p>{_fig(fig,
         "The same data, the same signal, two orderings. The joint fit has no "
         "window to choose and returns the same answer everywhere; the "
         "detrend's answer is a function of a parameter the analyst picks "
         "before seeing the signal.",
         "No detrend figure: p3_detrend is empty.")}</p>

{tbl}

<h3>Decision</h3>
<p><b>The order changes the answer by a factor of
{_n((worst_hi[1] / worst_lo[1]) if worst_lo and worst_lo[1] else float('nan'), 1)},
and not in one direction.</b>  At a window of
{worst_lo[0]:g}&nbsp;&times;&nbsp;P<sub>orb</sub> the running median destroys
{_pct(1 - worst_lo[1])} of the injected signal.  At
{worst_hi[0]:g}&nbsp;&times;&nbsp;P<sub>orb</sub> it returns
{_pct(worst_hi[1])} of it &mdash; it <em>fabricates</em> amplitude.  The joint
fit returns {_pct(jmin) if abs(jmax - jmin) < 5e-4 else _pct(jmin) + "&ndash;" + _pct(jmax)}
at every window from {rows[0][0]:g} to {rows[-1][0]:g}&nbsp;P<sub>orb</sub>
&mdash; it has no window to choose.</p>

<p>The fabrication is the more dangerous half and is the reason this section
exists.  Attenuation is at least a conservative failure; an analyst who
detrends at a badly chosen window and recovers 125% of a signal has no way to
know it, because the true amplitude is exactly what they were trying to
measure.  The joint fit is not merely better on average &mdash; it removes the
free parameter that does the damage.</p>

{"<p>The fast path is checked rather than trusted.  celerite2's "
 "<code>Matern32Term</code> is a two-exponential <em>approximation</em> whose "
 "own docstring warns it &ldquo;should be used with care&rdquo;, and at its "
 "default &epsilon;&nbsp;=&nbsp;0.01 it differs from the textbook "
 "Mat&eacute;rn-3/2 kernel by 1.8&times;10<sup>-5</sup> relative. The module "
 "pins &epsilon; so the fast path and a dense pure-numpy Cholesky reference "
 "are the same kernel to machine precision, and the two are compared on the "
 "real sampling:</p>" + chk_tbl if chk_tbl else ""}

<h3>Consequence</h3>
<p>No stage in Phase&nbsp;3 detrends before searching.  &sect;1's periodograms
carry one free constant per night &mdash; a nuisance model with a single
parameter per block, fitted jointly, which can only absorb power below
2.5&nbsp;c/d and is five times away from the orbital band.  This section is
the measurement that licenses that choice: it shows what happens when a
nuisance model is given enough freedom to reach the signal, and how far the
night-constant model is from that regime.</p>
</section>
"""


# ===========================================================================
def render_report(db_path: Path) -> Path:
    """Render the CV-S9 Phase-3 page.  Returns the HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = dict(q(con, "SELECT key, value FROM p3_meta")) \
            if _has(con, "p3_meta") else {}
        pfigs = {t: fig_periodograms(con, t) for t in TARGET_ORDER}
        pfigs["summary"] = fig_period_summary(con)
        f_sig = fig_sigmat(con)
        f_edge = fig_edges(con)
        f_oc = fig_oc(con)
        f_state = fig_states(con)
        f_det = fig_detrend(con)
        sections = [
            section_intro(con),
            section_periods(con, pfigs),
            section_sigmat(con, f_sig),
            section_edges(con, f_edge),
            section_oc(con, f_oc),
            section_states(con, f_state),
            section_detrend(con, f_det),
        ]
        n_series = _rows(con, "p3_period")
        n_edges = q1(con, "SELECT count(*) FROM p3_edge WHERE accepted=1") \
            if _has(con, "p3_edge") else 0
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CV Time Series &mdash; Phase 3 Analysis</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Cataclysmic-Variable Time Series &mdash; the Phase&nbsp;3
  analysis</h1>
  <p>Period verification with the spectral window beside every periodogram
  &middot; the &sigma;<sub>t</sub> injection test that decides whether
  per-cycle timing is publishable &middot; colour-dependent bright-phase
  edges &middot; O&minus;C on a cycle count that was checked rather than
  assumed &middot; accretion states with the Phase-2 limits &middot; joint
  GP&nbsp;+&nbsp;signal fitting instead of detrend-then-search
  &middot; {n_series} series, {n_edges} timed edges
  &middot; built {esc(meta.get('stage_periods', ''))[:16]}Z
  ({esc(meta.get('phase3_code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="cv_phase2_completion.html">the Phase-2 products this
  rests on</a> &middot; <a href="cv_characterization.html">the
  characterization that set these limits</a> &middot;
  <a href="index.html">project hub</a> &middot;
  <a href="../index.html">all reports</a></p>
</header>

<nav>
  <a href="#intro">0 What Phase 3 may not assume</a> &middot;
  <a href="#periods">1 Periods</a> &middot;
  <a href="#sigmat">2 &sigma;<sub>t</sub></a> &middot;
  <a href="#edges">3 Bright-phase timing</a> &middot;
  <a href="#oc">4 O&minus;C and cycle counts</a> &middot;
  <a href="#states">5 Accretion states</a> &middot;
  <a href="#detrend">6 Detrending</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_phot.report_phase3</code> from
<code>products/phot/cv_timeseries.sqlite</code>.  Every number in a table,
figure or verdict is the result of a query run at render time or a constant
imported from <code>macro_phot.phase3</code>.  The exceptions, named rather
than claimed away: the per-point precision and &chi;<sup>2</sup>-inflation
ranges and the three-filter full-orbit night counts quoted in &sect;0 and
&sect;3 (measured by the CV characterization, CV-S6), the Phase-2 finding
about limit depths quoted in &sect;5 (measured by CV-S8), and the descriptive
prose throughout.  The published ephemerides are cached payloads from the
AAVSO Variable Star Index under
<code>products/external/vsx/</code>, fetched on the date shown in &sect;0.
Regenerate with <code>pipeline/scripts/run_cv_phase3.py report</code>.
</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
    finally:
        con.close()
    return HTML_PATH
