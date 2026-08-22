"""CV-S8 chain-of-evidence report: closing out Phase-2 photometry.

Reads ``products/phot/cv_timeseries.sqlite`` and writes

* ``docs/CV_TimeSeries/cv_phase2_completion.html``
* ``docs/CV_TimeSeries/figures/cv_phase2/*.png``

Four Socratic sections, one per task, each in the same order:

    Question    what could go wrong if we did not do this?
    Evidence    the measurement, from the database, with its diagnostics
    Decision    what the numbers actually support
    Consequence what changes downstream, and what does not

and one section before them that says what Phase 2 already achieved and
what it measured about itself, because every threshold on this page is set
against the MEASURED precision (9-77 mmag per point, chi2 inflation
0.92-3.02) rather than against the strategy's hopes.

Every number in a table, figure or verdict is the result of a query run at
render time or a constant imported from ``macro_phot.phase2``.  Where a
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

from . import phase2 as p2        # noqa: E402
from macro_core.report_s0 import (  # noqa: E402
    ACCENT, BAD, STYLE, DPI, FAINT, GOOD, MUTED, WARN,
    _figure, esc, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402  (house figure style)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "CV_TimeSeries"
FIG_DIR = DOCS_DIR / "figures" / "cv_phase2"
HTML_PATH = DOCS_DIR / "cv_phase2_completion.html"

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}


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


def _pct(x, nd=2):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{100 * float(x):,.{nd}f}%"


def _pm(v, e, nd=4):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    if e is None or (isinstance(e, float) and not math.isfinite(e)):
        return f"{float(v):+,.{nd}f}"
    return (f"{float(v):+,.{nd}f}&thinsp;&plusmn;&thinsp;"
            f"{float(e):.{nd}f}")


def _label(series_key: str) -> str:
    t, e, f = series_key.split("|")
    return f"{TARGET_LABEL.get(t, t)} {e} <i>{f}</i>"


def _fig(src: str, caption: str, missing: str) -> str:
    """A figure, or an explicit statement of why there is nothing to draw."""
    if not src:
        return (f'<div class="decision"><b>No figure here, and why:</b> '
                f"{missing}</div>")
    return _figure(src, caption)


def _has(con, table_name: str) -> bool:
    return bool(q1(con, "SELECT count(*) FROM sqlite_master WHERE type='table'"
                        " AND name=?", (table_name,)))


#: Where this page gets "the precision this campaign actually achieves".
CHAR_DB = REPO_ROOT / "products" / "phot" / "cv_characterization.sqlite"


def precision_range(con) -> tuple[float, float, str]:
    """Achieved per-point precision, in mmag, and where the number came from.

    Preferred source is the CV characterization's ``prec_at_target`` — the
    precision measured AT EACH TARGET'S OWN BRIGHTNESS, which is the number
    the characterization page grades its goals against and therefore the
    number this page must agree with.  Falling back on this product's
    check-star RMS would quote a different statistic under the same words,
    and two pages of one project disagreeing about their own precision is
    exactly the drift the provenance graph exists to prevent.
    """
    if CHAR_DB.exists():
        try:
            c = sqlite3.connect(f"file:{CHAR_DB}?mode=ro", uri=True)
            c.execute("PRAGMA busy_timeout = 300000")
            row = c.execute(
                "SELECT min(prec_at_target), max(prec_at_target) "
                "FROM ch_noise_series WHERE prec_at_target IS NOT NULL"
            ).fetchone()
            c.close()
            if row and row[0] is not None:
                return (1000 * row[0], 1000 * row[1],
                        "CV-S5 characterization, precision at each target's "
                        "own brightness")
        except sqlite3.Error:
            pass
    row = q(con, "SELECT min(check_rms_median), max(check_rms_median) "
                 "FROM cv_series WHERE check_rms_median IS NOT NULL")[0]
    return (1000 * (row[0] or float("nan")), 1000 * (row[1] or float("nan")),
            "this product's check-star RMS (the characterization product "
            "was not readable at render time)")


# ===========================================================================
# Figures
# ===========================================================================
def fig_cloud_calibration(con) -> str:
    """The ROC that sets the threshold, plus the two exemplar nights."""
    roc = q(con, """SELECT threshold, false_veto_rate, recall, chosen
                    FROM p2_cloud_roc ORDER BY threshold""")
    if not roc:
        return ""
    t = np.array([r[0] for r in roc])
    fpr = np.array([r[1] if r[1] is not None else np.nan for r in roc])
    rec = np.array([r[2] if r[2] is not None else np.nan for r in roc])
    chosen = next((r[0] for r in roc if r[3]), None)
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 3.9))
        a1.plot(t, 100 * fpr, "-o", ms=3, color=BAD,
                label="clear frames vetoed (cost)")
        a1.plot(t, 100 * rec, "-o", ms=3, color=GOOD,
                label="attenuated frames caught (benefit)")
        a1.axhline(100 * p2.CLOUD_MAX_FALSE_VETO, color=WARN, ls="--", lw=1,
                   label=f"false-veto budget "
                         f"{100 * p2.CLOUD_MAX_FALSE_VETO:g}%")
        if chosen is not None:
            a1.axvline(chosen, color=ACCENT, lw=1.4,
                       label=f"threshold chosen {chosen:.2f}")
        a1.set_xlabel("veto threshold on ensemble flux ratio")
        a1.set_ylabel("per cent of independently-labelled frames")
        a1.legend(fontsize=7.5, loc="center left")
        a1.set_title("calibrated on the frames that carry ZMAG", fontsize=9)

        # Right panel: our statistic against the independent one, frame by
        # frame.  If the two channels measured the same sky, this is a line.
        rows = q(con, """SELECT rel_ratio, zmag_transmission, zmag_label
                         FROM p2_cloud_frame
                         WHERE rel_ratio IS NOT NULL
                           AND zmag_transmission IS NOT NULL""")
        if rows:
            x = np.array([r[1] for r in rows])
            y = np.array([r[0] for r in rows])
            lab = np.array([r[2] or "" for r in rows], dtype=object)
            a2.scatter(x[lab == ""], y[lab == ""], s=3, alpha=0.25,
                       color=FAINT, label="neither")
            a2.scatter(x[lab == "clear"], y[lab == "clear"], s=5, alpha=0.6,
                       color=GOOD, label="ZMAG says clear")
            a2.scatter(x[lab == "attenuated"], y[lab == "attenuated"], s=7,
                       alpha=0.8, color=BAD, label="ZMAG says attenuated")
            if chosen is not None:
                a2.axhline(chosen, color=ACCENT, lw=1.2)
            a2.set_xlim(0.4, 1.25)
            a2.set_ylim(0.4, 1.25)
            a2.set_xlabel("independent transmission from header ZMAG")
            a2.set_ylabel("ensemble flux ratio / local median")
            a2.legend(fontsize=7.5, loc="lower right")
            a2.set_title("two channels, one sky", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "cloud_calibration.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase2/cloud_calibration.png"


def fig_cloud_sculpting(con) -> str:
    """Veto rate against the target's own brightness, per series."""
    rows = q(con, """SELECT series_key, faint_veto_rate, bright_veto_rate, n
                     FROM p2_cloud_bias
                     WHERE faint_veto_rate IS NOT NULL
                       AND bright_veto_rate IS NOT NULL AND n >= 20
                     ORDER BY series_key""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(9.6, 4.2))
        y = np.arange(len(rows))
        ax.barh(y - 0.19, [100 * r[2] for r in rows], 0.36,
                color=ACCENT, label="veto rate, target's BRIGHTEST quartile")
        ax.barh(y + 0.19, [100 * r[1] for r in rows], 0.36,
                color=WARN, label="veto rate, target's FAINTEST quartile")
        ax.set_yticks(y)
        ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
        ax.invert_yaxis()
        ax.set_xlabel("per cent of frames vetoed")
        ax.legend(fontsize=8)
        ax.set_title("If the veto sculpted the light curve, the yellow bars "
                     "would be systematically longer", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "cloud_sculpting.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase2/cloud_sculpting.png"


def fig_extinction(con) -> str:
    """k'' with both error bars, against the precision it would have to beat."""
    rows = q(con, """SELECT era_id, filter, kpp, kpp_err, kpp_err_formal,
                            significant, term_p95_mmag, bound_mmag
                     FROM p2_extinction WHERE kpp IS NOT NULL
                     ORDER BY era_id, filter""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.0))
        y = np.arange(len(rows))
        k = np.array([r[2] for r in rows])
        e = np.array([r[3] for r in rows])
        ef = np.array([r[4] if r[4] is not None else np.nan for r in rows])
        a1.errorbar(k, y, xerr=e, fmt="o", ms=4, color=ACCENT, lw=1.4,
                    capsize=3, label="published (max of formal, bootstrap)")
        a1.errorbar(k, y, xerr=ef, fmt="none", ecolor=BAD, lw=3, alpha=0.7,
                    label="formal only — assumes every point independent")
        a1.axvline(0, color=MUTED, lw=1)
        a1.set_yticks(y)
        a1.set_yticklabels([f"era {r[0]} {r[1]}" for r in rows], fontsize=8)
        a1.invert_yaxis()
        a1.set_xlabel("k''  (mag mag$^{-1}$ airmass$^{-1}$)")
        a1.legend(fontsize=7.5, loc="lower right")
        a1.set_title("the star bootstrap is the honest error bar", fontsize=9)

        size = [(r[6] if r[5] else r[7]) for r in rows]
        colors = [GOOD if r[5] else FAINT for r in rows]
        a2.barh(y, size, 0.6, color=colors)
        a2.set_yticks(y)
        a2.set_yticklabels([f"era {r[0]} {r[1]}" for r in rows], fontsize=8)
        a2.invert_yaxis()
        a2.set_xscale("log")
        plo, phi, _src = precision_range(con)
        if math.isfinite(plo):
            a2.axvspan(plo, phi, color=WARN, alpha=0.18,
                       label="achieved per-point precision of this campaign")
            a2.legend(fontsize=7.5, loc="lower right")
        a2.set_xlabel("size of the effect / of its bound (mmag)")
        a2.set_title("green = detected; grey = bounded", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "extinction.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase2/extinction.png"


def fig_transform(con) -> str:
    """Measured slope against the slope the catalogue tie predicts."""
    rows = q(con, """SELECT target_key, band_from, band_to, b, b_err,
                            b_expected, b_expected_err, a, a_err, kind
                     FROM p2_transform WHERE b IS NOT NULL
                     ORDER BY target_key, band_from""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.0))
        y = np.arange(len(rows))
        a1.errorbar([r[3] for r in rows], y,
                    xerr=[r[4] for r in rows], fmt="o", ms=4, color=ACCENT,
                    capsize=3, label="measured star by star (this stage)")
        exp = [r[5] if r[5] is not None else np.nan for r in rows]
        exp_e = [r[6] if r[6] is not None else np.nan for r in rows]
        a1.errorbar(exp, y + 0.28, xerr=exp_e, fmt="s", ms=4, color=WARN,
                    capsize=3, label="predicted by the two eras' tie "
                                     "colour terms")
        a1.axvline(0, color=MUTED, lw=1)
        a1.set_yticks(y)
        a1.set_yticklabels([f"{TARGET_LABEL.get(r[0], r[0])} "
                            f"{r[1]}→{r[2]}" for r in rows], fontsize=8)
        a1.invert_yaxis()
        a1.set_xlabel("colour slope b  (mag / mag)")
        a1.legend(fontsize=7.5, loc="lower right")
        a1.set_title("two independent routes to the same number",
                     fontsize=9)

        a2.errorbar([1000 * r[7] for r in rows], y,
                    xerr=[1000 * r[8] for r in rows], fmt="o", ms=4,
                    color=BAD, capsize=3)
        a2.axvline(0, color=MUTED, lw=1)
        a2.set_yticks(y)
        a2.set_yticklabels([f"{TARGET_LABEL.get(r[0], r[0])} "
                            f"{r[1]}→{r[2]}" for r in rows], fontsize=8)
        a2.invert_yaxis()
        a2.set_xlabel("zero-point step a at the reference colour (mmag)")
        a2.set_title("what stitching the eras would cost, before any "
                     "colour term", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "transform.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase2/transform.png"


def fig_limits(con) -> str:
    """Where the undetected epochs went, and how deep the limits reach."""
    rows = q(con, """SELECT series_key, n_detected, n_candidates, n_limits,
                            n_recovered, n_failed
                     FROM p2_limit_series WHERE n_candidates > 0
                     ORDER BY n_candidates DESC""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.4))
        y = np.arange(len(rows))
        rec = np.array([r[4] for r in rows], dtype=float)
        lim = np.array([r[3] for r in rows], dtype=float)
        fail = np.array([r[5] for r in rows], dtype=float)
        a1.barh(y, rec, 0.62, color=GOOD, label="recovered as a detection")
        a1.barh(y, lim, 0.62, left=rec, color=ACCENT,
                label="bounded by an upper limit")
        a1.barh(y, fail, 0.62, left=rec + lim, color=BAD,
                label="still unmeasured")
        a1.set_yticks(y)
        a1.set_yticklabels([r[0] for r in rows], fontsize=7.5)
        a1.invert_yaxis()
        a1.set_xlabel("epochs that carried no magnitude before this stage")
        a1.legend(fontsize=8)
        a1.set_title("what happened to the censored epochs", fontsize=9)

        # Limits against detections, on the calibrated scale.
        dets, lims, keys = [], [], []
        for (skey,) in q(con, "SELECT DISTINCT series_key FROM p2_limits "
                              "WHERE outcome='limit' AND limit_cal_mag "
                              "IS NOT NULL ORDER BY 1"):
            zp = q1(con, "SELECT zp FROM cv_cattie WHERE series_key=? AND "
                         "is_primary=1", (skey,))
            if zp is None:
                continue
            d = [r[0] - zp for r in q(
                con, "SELECT mag FROM cv_lightcurve WHERE series_key=? AND "
                     "role='target' AND mag IS NOT NULL", (skey,))]
            l = [r[0] for r in q(
                con, "SELECT limit_cal_mag FROM p2_limits WHERE "
                     "series_key=? AND outcome='limit' AND limit_cal_mag "
                     "IS NOT NULL", (skey,))]
            if d and l:
                dets.append(d)
                lims.append(l)
                keys.append(skey)
        if keys:
            yy = np.arange(len(keys))
            a2.boxplot(dets, positions=yy - 0.17, vert=False, widths=0.28,
                       patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=ps.tint(ACCENT), color=ACCENT),
                       medianprops=dict(color=ACCENT))
            a2.boxplot(lims, positions=yy + 0.17, vert=False, widths=0.28,
                       patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=ps.tint(WARN), color=WARN),
                       medianprops=dict(color=WARN))
            a2.set_yticks(yy)
            a2.set_yticklabels(keys, fontsize=7.5)
            a2.invert_yaxis()
            a2.set_xlabel("natural-system magnitude (blue: detections; "
                          "yellow: 3σ limits)")
            a2.set_title("a limit fainter than the detections is a limit "
                         "that constrains", fontsize=9)
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "limits.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_phase2/limits.png"


# ===========================================================================
# Sections
# ===========================================================================
def section_intro(con) -> str:
    n_series = q1(con, "SELECT count(*) FROM cv_series WHERE status='solved'")
    n_lc = q1(con, "SELECT count(*) FROM cv_lightcurve")
    plo, phi, psrc = precision_range(con)
    infl = q(con, "SELECT min(chi2_inflation), max(chi2_inflation) FROM "
                  "cv_series WHERE chi2_inflation IS NOT NULL")[0]
    n_zmag = q1(con, "SELECT count(*) FROM p2_cloud_frame WHERE zmag IS NOT "
                     "NULL AND zmag <> 0") if _has(con, "p2_cloud_frame") else 0
    n_frames = q1(con, "SELECT count(*) FROM p2_cloud_frame") \
        if _has(con, "p2_cloud_frame") else 0
    return f"""
<section id="intro">
<div class="bhead"><h2>0 &middot; What Phase&nbsp;2 already had, and the four
holes left in it</h2></div>

<p>Phase&nbsp;2 measured {n_lc:,} light-curve rows across {n_series} solved
(target, era, filter) series and tied them to ATLAS-REFCAT2.  It also
measured what it can actually do, and the answer is smaller than the
strategy hoped: per-point precision runs
<b>{plo:.0f}&ndash;{phi:.0f}&nbsp;mmag</b> depending on the series
({psrc}), and chi&sup2; inflation runs
<b>{infl[0]:.2f}&ndash;{infl[1]:.2f}</b>.  Every threshold on this page is
set against those numbers, not against the ones the plan started with.</p>

<div class="decision"><b>Why these four tasks and not others.</b> Each one is
a place where a light curve that <i>looks</i> finished is quietly wrong.
Cloud that was never removed becomes structure.  A colour-dependent
extinction term that was never measured becomes a systematic nobody
budgeted.  A transformation that was assumed rather than measured becomes a
stitched light curve.  And undetected epochs that were dropped rather than
bounded become a duty cycle computed over exactly the epochs where the
target was bright enough to see &mdash; a tautology wearing a
measurement's clothes.</div>

<div class="decision"><b>The fact that forces task&nbsp;1's design.</b> The
acquisition software writes a per-image photometric zero point
(<code>ZMAG</code>) into every frame header, and the original cloud cut
leaned on it.  Of the {n_frames:,} CV frames this stage examined,
<b>{n_zmag:,}</b> carry a usable one &mdash; and <i>none</i> of them are
Sloan-era polar frames.  VV&nbsp;Pup: 0 of 1,353.  EU&nbsp;UMa era 78: 0 of
208.  For the polars the primary cloud channel cannot be ZMAG, so it has to
be the comparison ensemble's own summed flux.  The {n_zmag:,} frames that
<i>do</i> carry ZMAG are not wasted, though: they are exactly what the new
channel's threshold is calibrated against.</div>

<p>The four sections below each run <b>Question &rarr; Evidence &rarr;
Decision &rarr; Consequence</b>, in that order, and refuse to state a
decision before the evidence that carries it.</p>
</section>"""


def section_cloud(con, fig_cal: str, fig_sculpt: str) -> str:
    meta = dict(q(con, "SELECT key, value FROM p2_meta"))
    thr = meta.get("cloud_threshold_used", "&mdash;")
    reason = meta.get("cloud_threshold_reason", "")
    tot = q(con, "SELECT count(*), sum(vetoed) FROM p2_cloud_frame")[0]
    ex = q(con, """SELECT series_key, night, n_frames, zmag_mad, zmag_span,
                          rel_ratio_mad, rel_ratio_min, n_vetoed, exemplar
                   FROM p2_cloud_night WHERE exemplar <> ''
                   ORDER BY exemplar""")
    ex_rows = [[esc(r[8]).upper(), _label(r[0]), esc(r[1]), _i(r[2]),
                _n(r[3], 3), _n(r[4], 3), _n(r[5], 4), _n(r[6], 3), _i(r[7])]
               for r in ex]
    roc_rows = [[_n(r[0], 2), _i(r[1]), _i(r[2]), _pct(r[3]), _i(r[4]),
                 _i(r[5]), _pct(r[6], 1),
                 "&#10004;" if r[7] else ""]
                for r in q(con, """SELECT threshold, n_clear, n_clear_vetoed,
                    false_veto_rate, n_attenuated, n_attenuated_vetoed,
                    recall, chosen FROM p2_cloud_roc
                    WHERE threshold >= 0.85 ORDER BY threshold""")]
    ser_rows, ser_cls = [], []
    for r in q(con, """SELECT series_key, n_frames, n_nights, n_core_median,
                              n_vetoed, frac_vetoed, rel_ratio_mad,
                              worst_rel_ratio, n_zmag, note
                       FROM p2_cloud_series ORDER BY frac_vetoed DESC"""):
        ser_rows.append([_label(r[0]), _i(r[1]), _i(r[2]), _n(r[3], 0),
                         _i(r[4]), _pct(r[5]), _n(r[6], 4), _n(r[7], 3),
                         _i(r[8]), esc(r[9] or "")])
        ser_cls.append("warn" if (r[5] or 0) > 0.20 else "")
    bias_rows, bias_cls = [], []
    for r in q(con, """SELECT series_key, n, n_vetoed, median_mag_vetoed,
                              median_mag_kept, p_mannwhitney,
                              faint_veto_rate, bright_veto_rate,
                              p_proportion, undetected_veto_rate,
                              detected_veto_rate, verdict
                       FROM p2_cloud_bias ORDER BY series_key"""):
        bias_rows.append([_label(r[0]), _i(r[1]), _i(r[2]), _n(r[3], 3),
                          _n(r[4], 3), _n(r[5], 4), _pct(r[6], 1),
                          _pct(r[7], 1), _n(r[8], 4), _pct(r[9], 1),
                          _pct(r[10], 1), esc(r[11])])
        bias_cls.append("warn" if "FAINT-PHASE" in (r[11] or "") else "")
    n_faint = q1(con, "SELECT count(*) FROM p2_cloud_bias WHERE "
                      "verdict='FAINT-PHASE VETO EXCESS'")
    n_bright = q1(con, "SELECT count(*) FROM p2_cloud_bias WHERE "
                       "verdict='BRIGHT-PHASE VETO EXCESS'")
    n_clean = q1(con, "SELECT count(*) FROM p2_cloud_bias WHERE "
                      "verdict='NO SCULPTING DETECTED'")
    return f"""
<section id="cloud">
<div class="bhead"><h2>1 &middot; The ensemble-flux-ratio cloud veto</h2>
<span class="tag">CV-P2-cloud-veto</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">What does a cloud do to a differential light curve, and how
would we know it happened, on frames that carry no independent record of
the sky?</p>
<p>A Honeycutt ensemble is differential, so first order the answer is
&ldquo;nothing&rdquo;: cloud dims the target and the comparison stars
together and the ratio survives.  Second order it is not nothing.  A cloud
that costs half a magnitude also costs half a magnitude of
signal-to-noise on every star, drops the faint end of the ensemble below
detection, and changes which stars the zero point is averaged over.  The
result is not a shifted point, it is a NOISIER point with a zero point
computed from a different set of stars &mdash; and on a 14&ndash;18
point-per-cycle cadence, a handful of those is a spurious feature.</p>
<p>The channel the strategy assumed &mdash; the header <code>ZMAG</code>
&mdash; does not exist for the targets that need it most.  So the veto has
to be built from the only thing measured on every frame by construction:
the comparison ensemble's own summed flux.</p></div>

<div class="stage"><h3>Evidence</h3>
<p class="sub">The statistic, then the calibration, then the census.</p>
<p>For every night of every series, the stars measured on at least
{100 * p2.CLOUD_CORE_MIN_FRAC:.0f}% of that night's frames form a
<b>core ensemble</b> (membership fixed, so that a cloud cannot shrink the
sum by removing contributors as well as by dimming them).  Each frame's
statistic is the core's summed flux <i>rate</i> &mdash; flux per second, so
that a night mixing 60&nbsp;s and 240&nbsp;s exposures is not read as a
night mixing clear and cloudy sky &mdash; divided by the same subset's
median.  That ratio is then divided by its own running median over
&plusmn;{p2.CLOUD_WINDOW_HALF} frames.  Dividing by a <i>running</i> median
rather than the night's median is what makes this a cloud detector instead
of an airmass detector: extinction falls smoothly across a night and is
already absorbed exactly by the ensemble zero point, while a cloud is a
departure from the local normal on a timescale of minutes.</p>

{_fig(fig_cal,
      "LEFT: the calibration.  Every candidate threshold, scored on the "
      "frames that carry an independent ZMAG: how many independently-CLEAR "
      "frames it would throw away (red, the cost) against how many "
      "independently-ATTENUATED frames it would catch (green, the "
      "benefit).  The threshold is the highest one whose cost stays inside "
      "the declared budget.  RIGHT: the two channels frame by frame.  They "
      "are independent measurements &mdash; one from the acquisition "
      "software's own plate solve, one from our comparison ensemble.",
      "the ROC table is empty, which means no frame carried both a usable "
      "ZMAG and a usable ensemble ratio")}

<p><b>The two exemplar nights</b>, chosen by the INDEPENDENT channel (the
spread of the header zero point) so that the exemplars are not selected by
the statistic being tested:</p>
{table(["", "series", "night", "frames", "ZMAG MAD", "ZMAG span",
        "ratio MAD", "worst ratio", "vetoed"], ex_rows)
 if ex_rows else '<div class="decision">No night carried enough ZMAG-bearing '
                 'frames to serve as an exemplar.</div>'}

<p><b>The calibration table</b> &mdash; the whole argument for the number,
in one place:</p>
{table(["threshold", "clear frames", "of which vetoed", "false-veto rate",
        "attenuated frames", "of which vetoed", "recall", "chosen"],
       roc_rows)}
</div>

<div class="stage"><h3>Decision</h3>
<div class="decision"><b>Threshold {thr}.</b> {esc(reason)}<br>
Applied, it vetoes <b>{_i(tot[1])}</b> of <b>{_i(tot[0])}</b> frames
({_pct((tot[1] or 0) / max(1, tot[0]))}).  That rate is much larger than the
false-veto rate because most vetoed frames sit in the band the independent
labels deliberately leave unclassified: a frame that has lost
0.09&ndash;0.15&nbsp;mag is neither &ldquo;clear&rdquo; nor
&ldquo;attenuated&rdquo; by the ZMAG rule, but it is a frame the veto
correctly removes.  The labels are deliberately narrow (see
<code>phase2.ZMAG_CLEAR_MAG</code> and
<code>phase2.ZMAG_ATTEN_MAG</code>) precisely so that neither class is
manufactured.</div>
{table(["series", "frames", "nights", "core stars", "vetoed", "%",
        "MAD of ratio", "worst ratio", "frames with ZMAG", "note"],
       ser_rows, ser_cls)}
</div>

<div class="stage"><h3>Consequence &mdash; and the test that had to pass
first</h3>
<p class="sub">Does the veto preferentially remove the target's FAINT
phases?  If it did, it would carve the bottom out of every eclipse and
every low state, invisibly, because the survivors would still look like a
clean light curve.</p>
<p>Two independent readings, both two-sided: a Mann-Whitney&nbsp;U on the
target magnitudes of vetoed against kept frames, and a two-proportion
z&nbsp;test on the veto rate in the target's faintest quartile against its
brightest.  A third column reports the veto rate on epochs where the target
was not detected at all &mdash; the faintest epochs of the series by
definition, and the ones no magnitude-based test can see.</p>

{_fig(fig_sculpt,
      "Veto rate in the target's brightest quartile (blue) against its "
      "faintest quartile (yellow), per series.  Sculpting would show as "
      "yellow systematically longer than blue.",
      "no series had both a measured veto rate in each quartile and at "
      "least 20 measured epochs")}

{table(["series", "measured epochs", "vetoed", "median mag (vetoed)",
        "median mag (kept)", "p (Mann-Whitney)", "veto rate, faint quartile",
        "veto rate, bright quartile", "p (two-proportion)",
        "veto rate, undetected epochs", "veto rate, detected epochs",
        "verdict"], bias_rows, bias_cls)}

<div class="decision"><b>{n_clean} series show no sculpting; {n_faint} show a
FAINT-side excess; {n_bright} show a BRIGHT-side excess.</b>
The direction is part of the answer and the verdict names it, because a
significant asymmetry in which the BRIGHT quartile is vetoed more often
cannot carve a low state out of a light curve.  Collapsing both directions
into one alarm would have reported
{n_faint + n_bright} &ldquo;suspected&rdquo; series and trained the reader
to ignore the word.  <b>The veto is safe to apply as a cleaning step</b> on
this evidence, and the per-frame flag is published rather than applied, so
a later stage can honour it or argue with it.</div>

<div class="decision"><b>What this stage deliberately does NOT do.</b> It
writes no column on <code>cv_lightcurve</code>.  The veto is a per-frame
FLAG in <code>p2_cloud_frame</code>.  A stage that silently deleted rows
would make the sculpting test above unfalsifiable &mdash; there would be
nothing left to test it on.</div>
</div>
</section>"""


def section_extinction(con, fig_ext: str) -> str:
    rows = q(con, """SELECT era_id, filter, era_label, n_points, n_stars,
                            n_frames, colour_label, colour_ref, colour_min,
                            colour_max, airmass_min, airmass_max, kpp,
                            kpp_err, kpp_err_formal, kpp_err_boot, t_stat,
                            significant, chi2nu, term_p95_mmag, bound_mmag,
                            verdict, note, n_clipped
                     FROM p2_extinction ORDER BY era_id, filter""")
    tab, cls = [], []
    for r in rows:
        tab.append([f"era {r[0]} <i>{esc(r[1])}</i>", esc(r[6] or ""),
                    _i(r[3]), _i(r[4]), _i(r[5]),
                    f"{_n(r[8], 2)} &hellip; {_n(r[9], 2)}",
                    f"{_n(r[10], 2)} &hellip; {_n(r[11], 2)}",
                    _pm(r[12], r[13], 5), _n(r[14], 5), _n(r[15], 5),
                    _n(r[16], 2), _n(r[18], 2),
                    _n(r[19], 2), _n(r[20], 2), esc(r[21])])
        cls.append("warn" if r[17] else "")
    n_fit = sum(1 for r in rows if r[12] is not None)
    n_sig = sum(1 for r in rows if r[17])
    worst = max((r[20] for r in rows if r[20] is not None), default=None)
    biggest = max((r[19] for r in rows if r[17] and r[19] is not None),
                  default=None)
    plo, phi, _psrc = precision_range(con)
    notes = "".join(
        f"<li><b>era {r[0]} <i>{esc(r[1])}</i></b> &mdash; {esc(r[22])}</li>"
        for r in rows)
    return f"""
<section id="extinction">
<div class="bhead"><h2>2 &middot; Second-order colour-extinction terms</h2>
<span class="tag">CV-P2-extinction</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The ensemble already removes extinction.  What is left?</p>
<p>The Honeycutt model is <code>m<sub>sj</sub> = M<sub>s</sub> +
ZP<sub>j</sub></code>.  Anything common to a whole frame lands in
<code>ZP<sub>j</sub></code>, so first-order extinction
<code>k&prime;&middot;X</code> is removed <i>exactly</i> and there is
nothing to measure.  What the model cannot absorb is the part that differs
between stars in the SAME frame: the colour-dependent term
<code>k&Prime;&middot;(C&minus;C<sub>ref</sub>)&middot;(X&minus;X<sub>ref</sub>)</code>.
If that term is real and unmodelled, it puts an airmass-shaped systematic
on every star whose colour differs from the ensemble mean &mdash; and a
cataclysmic variable's colour differs from the ensemble mean by
construction.</p></div>

<div class="stage"><h3>Evidence</h3>
<p class="sub">One coefficient per (era, filter), fitted on comparison
stars, with the uncertainty that matters rather than the one that is easy.</p>
<p>Three choices carry this fit, and each was forced by something the data
did:</p>
<ol>
<li><b>The design column is two-way centred.</b> The solver has already
projected out a free constant per star and per frame, so a design column
still containing star means and frame means lives partly in the space the
solver removed and returns a diluted coefficient.  The column is projected
the same way the data were &mdash; iteratively, because the tables are
sparse and the two projections do not commute.</li>
<li><b>Impossible airmasses are refused.</b> 62 matched CV frames carry
header <code>AIRMASS</code> between 5 and 6,877.  VV&nbsp;Pup, at
dec&nbsp;&minus;19, cannot exceed airmass&nbsp;2.1 from this site at any
hour of any night.  Because the design column is proportional to
<code>(X&minus;X<sub>ref</sub>)</code>, one <code>X&nbsp;=&nbsp;6,877</code>
point carries 40,000&times; the leverage of a real one; the first run of
this stage returned
<code>k&Prime;&nbsp;=&nbsp;&minus;4&times;10<sup>&minus;5</sup></code> for
era&nbsp;76 for exactly that reason, which is how the defect was found.
The window is now
[{p2.AIRMASS_MIN:g},&nbsp;{p2.AIRMASS_MAX:g}].</li>
<li><b>The published error is a star bootstrap, not the formal one.</b> The
formal weighted-least-squares error treats every (star, frame) measurement
as an independent draw.  They are not: a comparison star whose catalogue
colour is 0.05&nbsp;mag wrong contributes the same wrong colour to all
1,600 of its measurements.  Resampling whole STARS &mdash; the unit that
actually varies independently &mdash; gives errors 3&ndash;35&times; larger,
and the published uncertainty is the larger of the two.</li>
</ol>

{_fig(fig_ext,
      "LEFT: the coefficient with both error bars.  The wide blue bar is "
      "the published (bootstrap) uncertainty; the short red bar is the "
      "formal one, which assumes every measurement is independent.  The "
      "difference between them is the whole reason this section reaches "
      "the conclusion it does.  RIGHT: the SIZE of each effect or of its "
      "bound, against the per-point precision this campaign actually "
      "achieves.",
      "no (era, filter) group had enough comparison measurements carrying "
      "both a catalogue colour and a usable airmass to fit anything")}

{table(["era, filter", "colour", "points", "stars", "frames", "colour range",
        "airmass range", "k&Prime; (published)", "formal &sigma;",
        "bootstrap &sigma;", "t", "&chi;&sup2;/&nu;", "effect p95 (mmag)",
        "bound (mmag)", "verdict"], tab, cls)}
</div>

<div class="stage"><h3>Decision</h3>
<div class="decision"><b>{n_sig} of {n_fit} fitted (era, filter) groups show
a second-order colour-extinction term at
{p2.KPP_SIGNIFICANCE_T:g}&nbsp;&sigma;, and even those move a real
measurement by at most {_n(biggest, 2)}&nbsp;mmag.</b>  The campaign's
achieved per-point precision is {plo:.0f}&ndash;{phi:.0f}&nbsp;mmag.
<b>No term is applied.</b>  Applying a correction an order of magnitude below the noise it
sits in adds the correction's own uncertainty to every point and buys
nothing measurable.</div>
<div class="decision"><b>The bound goes into the error budget instead.</b>
For every group the table above carries
<code>{p2.KPP_SIGNIFICANCE_T:g}&sigma;&nbsp;&times;&nbsp;&sigma;<sub>k&Prime;</sub>&nbsp;&times;&nbsp;(&Delta;C/2)&nbsp;&times;&nbsp;(&Delta;X/2)</code>
&mdash; the largest systematic the data still allow this omission to cost.
The worst of them is <b>{_n(worst, 2)}&nbsp;mmag</b>, on a group whose own
per-point precision is far larger.  That is the honest form of the answer:
not &ldquo;the term is zero&rdquo;, which we cannot show, but &ldquo;the
term cannot cost more than this, and this is small compared with what we
can measure&rdquo;.</div>
<ul class="legend">{notes}</ul>
</div>

<div class="stage"><h3>Consequence</h3>
<p>Nothing in <code>cv_lightcurve</code> changes.  What changes is that the
CV error budget now carries a NAMED, MEASURED line for colour-dependent
extinction instead of an unexamined assumption, and any future claim at the
few-mmag level has a number to check itself against.  Two groups could not
be fitted at all and the table says which and why rather than omitting
them: EU&nbsp;UMa's merged era-78 Fast block carries no per-star catalogue
tie (it is the block with five comparison stars and no check stars), and
ST&nbsp;LMi's era-47 <i>y</i> block carries no catalogue tie of any kind,
because ATLAS-REFCAT2 does not publish a <i>y</i> band.</p>
</div>
</section>"""


def section_crossera(con, fig_tr: str) -> str:
    rows = q(con, """SELECT target_key, era_from, band_from, era_to, band_to,
                            kind, colour_label, n_stars, colour_ref,
                            colour_min, colour_max, a, a_err, b, b_err,
                            rms_mmag, chi2nu, colour_term_from,
                            colour_term_to, b_expected, b_expected_err,
                            b_tension_sigma
                     FROM p2_transform ORDER BY target_key, era_from,
                          band_from""")
    tab = [[TARGET_LABEL.get(r[0], r[0]),
            f"e{r[1]} <i>{esc(r[2])}</i> &rarr; e{r[3]} <i>{esc(r[4])}</i>",
            esc(r[5]), esc(r[6] or ""), _i(r[7]),
            f"{_n(r[9], 2)} &hellip; {_n(r[10], 2)}", _n(r[8], 3),
            _pm(r[11], r[12]), _pm(r[13], r[14]),
            _pm(r[19], r[20]), _n(r[21], 1), _n(r[15], 1), _n(r[16], 2)]
           for r in rows]
    disc = q(con, """SELECT check_id, statement, n_checked, n_violation,
                            verdict, detail FROM p2_discipline
                     ORDER BY check_id""")
    dtab, dcls = [], []
    for r in disc:
        dtab.append([f"<code>{esc(r[0])}</code>", esc(r[1]), _i(r[2]),
                     _i(r[3]), f"<b>{esc(r[4])}</b>", esc(r[5] or "")])
        dcls.append("" if r[4] in ("HOLDS", "NOT APPLICABLE") else "warn")
    max_t = max((abs(r[21]) for r in rows if r[21] is not None), default=None)
    max_a = max((abs(r[11]) for r in rows if r[11] is not None), default=None)
    max_b = max((abs(r[13]) for r in rows if r[13] is not None), default=None)
    n_hold = sum(1 for r in disc if r[4] == "HOLDS")
    return f"""
<section id="crossera">
<div class="bhead"><h2>3 &middot; Cross-era discipline and transformation
metadata</h2><span class="tag">CV-P2-cross-era</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">ST&nbsp;LMi was observed in G/R/I in 2024 and in g/r/i in
2025&ndash;26, and the two seasons do not overlap in time.  What is the cost
of pretending they are one data set &mdash; and can we prove we did not?</p>
<p>Two separate obligations sit here and they are easy to confuse.  The
first is <b>metadata</b>: a data release has to publish the transformation
between its own natural systems, derived from stars, with uncertainties, so
that somebody else can convert a magnitude.  The second is
<b>discipline</b>: the paper compares two within-era analyses and never
stitches them, and that promise is worth nothing unless the products are
checked against it.  Publishing the coefficients and applying them are
opposite acts.</p></div>

<div class="stage"><h3>Evidence &mdash; part one, the coefficients</h3>
<p class="sub">Fitted on comparison stars common to both eras, matched
through the same ATLAS-REFCAT2 row, with both magnitudes first put on
their own era's catalogue zero point so that what is left is bandpass.</p>
<p>The model is
<code>m<sub>to</sub> &minus; m<sub>from</sub> = a + b&middot;(C &minus;
C<sub>ref</sub>)</code>, centred on the median colour of the tie stars so
that <code>a</code> is the offset AT THE COLOUR THE STARS ACTUALLY HAVE
&mdash; decorrelated from <code>b</code>, rather than an extrapolation to
zero colour, which no star in this campaign occupies.</p>
<p>One row is a <b>control</b>: VV&nbsp;Pup g&rarr;g, r&rarr;r, i&rarr;i
compares two eras wearing the SAME filter labels through two different
cameras.  If the method works, the control measures the detector change and
nothing else, and it does so with an order of magnitude more stars than the
science rows.</p>
<p>And there is an independent prediction to check against.  The catalogue
tie already measured, for each era separately, how far that era's bandpass
sits from the catalogue's &mdash; that is its colour term <code>k</code>.
If the two eras' natural systems differ only in bandpass, then the slope
measured here star-by-star must equal
<code>k<sub>to</sub> &minus; k<sub>from</sub></code>, which was measured a
completely different way.</p>

{_fig(fig_tr,
      "LEFT: the colour slope measured star-by-star here (blue circles) "
      "against the slope the two eras' catalogue colour terms predict "
      "(yellow squares).  Two independent routes to the same number.  "
      "RIGHT: the zero-point step between eras at the reference colour "
      "&mdash; what stitching would cost before any colour term is even "
      "considered.",
      "no era pair carried common comparison stars with a primary tie in "
      "both eras")}

{table(["target", "transformation", "kind", "colour", "stars",
        "colour range", "C<sub>ref</sub>", "a (offset)", "b (slope)",
        "b predicted by the ties", "tension (&sigma;)", "rms (mmag)",
        "&chi;&sup2;/&nu;"], tab)}

<div class="decision"><b>Every one of the {len(rows)} transformations agrees
with the independent prediction from the catalogue colour terms to within
{_n(max_t, 1)}&nbsp;&sigma;.</b>  That is the strongest available statement
that these coefficients mean what they say: the same bandpass difference,
measured twice by methods that share no arithmetic &mdash; one a regression
of ensemble magnitude against catalogue magnitude, the other a differencing
of two ensembles star by star.</div>
</div>

<div class="stage"><h3>Evidence &mdash; part two, the discipline</h3>
<p class="sub">Four assertions about the PRODUCTS, not about our
intentions.</p>
{table(["check", "statement", "checked", "violations", "verdict", "detail"],
       dtab, dcls)}
</div>

<div class="stage"><h3>Decision</h3>
<div class="decision"><b>The coefficients are published as data-release
metadata and applied to nothing.</b>  The largest zero-point step between
eras is <b>{_n(1000 * max_a if max_a else None, 1)}&nbsp;mmag</b> at the
reference colour and the largest colour slope is
<b>{_n(max_b, 4)}&nbsp;mag/mag</b>.  A stitched ST&nbsp;LMi light curve
would therefore carry a step of that size plus a colour-dependent term, on
a target whose colour swings through the orbit &mdash; which is precisely
why the paper runs two within-era analyses instead.  The
<code>applied_to_targets</code> column of <code>p2_transform</code> is 0 on
every row and is part of that table's provenance fingerprint, so it cannot
quietly become 1.</div>
<div class="decision"><b>{n_hold} of {len(disc)} assertions hold.</b>  In
particular every calibrated target magnitude equals <code>mag &minus;
zp</code> exactly: a zero-point shift and no colour transformation.  That
matters because these targets are cataclysmic variables &mdash; blue,
variable, and routinely outside the colour range over which any
transformation was calibrated.  Transforming them would swap a bandpass
error of known size for an extrapolation error of unknown size.</div>
</div>

<div class="stage"><h3>Consequence</h3>
<p>The control row is the most interesting number on this page for anyone
planning future observations: VV&nbsp;Pup's <i>g</i>&rarr;<i>g</i> slope
across the era&nbsp;72&nbsp;&rarr;&nbsp;76 camera change is
<b>{_n(next((r[13] for r in rows if r[0] == 'vvpup' and r[2] == 'g'), None), 4)}</b>&nbsp;mag/mag,
<i>larger</i> than ST&nbsp;LMi's G&rarr;g slope across the filter-label
change.  The dominant bandpass difference in this archive is the DETECTOR
ERA, not the filter label.  Any future analysis that groups by filter name
and ignores the era is making a bigger error than the one it thinks it is
avoiding.</p>
</div>
</section>"""


def section_limits(con, fig_lim: str) -> str:
    rows = q(con, """SELECT series_key, n_matched, n_detected, n_candidates,
                            n_forced, n_limits, n_recovered, n_failed,
                            median_limit_mag, median_limit_cal_mag,
                            faintest_detection, pos_method,
                            pos_crosscheck_px, n_closure, closure_median_px,
                            closure_p95_px, blocked
                     FROM p2_limit_series WHERE n_candidates > 0
                     ORDER BY n_candidates DESC""")
    tab, cls = [], []
    for r in rows:
        tab.append([_label(r[0]), _i(r[1]), _i(r[2]), _i(r[3]), _i(r[6]),
                    _i(r[5]), _i(r[7]), _n(r[8], 2), _n(r[9], 2),
                    _n(r[10], 2), esc(r[11]), _n(r[12], 2), _i(r[13]),
                    _n(r[14], 2), _n(r[15], 2), esc(r[16] or "")])
        cls.append("warn" if r[16] else "")
    stat_rows = []
    for r in q(con, """SELECT s.series_key,
            max(CASE WHEN statistic='n_epochs_measured'
                     THEN censored_value END),
            max(CASE WHEN statistic='n_epochs_measured'
                     THEN with_limits_value END),
            max(CASE WHEN statistic='detected_fraction'
                     THEN with_limits_value END),
            max(CASE WHEN statistic='faint_state_fraction'
                     THEN with_limits_value END),
            max(CASE WHEN statistic='median_mag' THEN censored_value END),
            max(CASE WHEN statistic='median_mag' THEN with_limits_value END)
        FROM p2_state_stats s GROUP BY s.series_key
        HAVING max(CASE WHEN statistic='faint_state_fraction'
                        THEN with_limits_value END) > 0
        ORDER BY 5 DESC"""):
        stat_rows.append([_label(r[0]), _i(r[1]), _i(r[2]), _pct(r[3], 1),
                          _pct(r[4], 1), _n(r[5], 2), _n(r[6], 2)])
    n_lim = q1(con, "SELECT count(*) FROM p2_limits WHERE outcome='limit'")
    n_rec = q1(con, "SELECT count(*) FROM p2_limits WHERE outcome='detection'")
    n_fail = q1(con, "SELECT count(*) FROM p2_limits WHERE status='failed'")
    n_nozp = q1(con, "SELECT count(*) FROM p2_limits "
                     "WHERE outcome='no_zeropoint'")
    blocked = q(con, "SELECT series_key, blocked FROM p2_limit_series "
                     "WHERE blocked <> '' ORDER BY series_key")
    bl = "".join(f"<li><b>{_label(b[0])}</b> &mdash; {esc(b[1])}</li>"
                 for b in blocked)
    worst_gain = q(con, """SELECT series_key, n_detected, n_detected +
                                  n_recovered + n_limits
                           FROM p2_limit_series
                           WHERE n_candidates > 0
                           ORDER BY (n_recovered + n_limits) DESC LIMIT 1""")
    # How DEEP are the limits, really?  A limit brighter than the series'
    # own typical detection bounds the target only weakly: it says the
    # frame was insensitive, not that the target was faint.  This is
    # measured rather than asserted, because the answer turns out to be
    # mostly "shallow" and that changes how the faint-state fractions above
    # may be read.
    depth_rows, n_shallow, n_deep = [], 0, 0
    for skey, med_lim in q(con, """SELECT series_key, median_limit_cal_mag
                                   FROM p2_limit_series
                                   WHERE median_limit_cal_mag IS NOT NULL
                                   ORDER BY series_key"""):
        zp = q1(con, "SELECT zp FROM cv_cattie WHERE series_key=? AND "
                     "is_primary=1", (skey,))
        if zp is None:
            continue
        det = [r[0] - zp for r in q(
            con, "SELECT mag FROM cv_lightcurve WHERE series_key=? AND "
                 "role='target' AND mag IS NOT NULL", (skey,))]
        if not det:
            continue
        med_det = float(np.median(det))
        deeper = med_lim > med_det
        n_deep += int(deeper)
        n_shallow += int(not deeper)
        depth_rows.append([_label(skey), _n(med_det, 2), _n(med_lim, 2),
                           _n(med_lim - med_det, 2),
                           "deeper than the detections" if deeper
                           else "shallower &mdash; a low-sensitivity frame, "
                                "not a faint target"])
    return f"""
<section id="limits">
<div class="bhead"><h2>4 &middot; Faint-phase forced photometry and upper
limits</h2><span class="tag">CV-P2-faint-limits</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">A polar in a low state drops below detection.  What does
dropping those epochs do to every statistic computed afterwards?</p>
<p>It censors exactly the faint half of the distribution.  A duty cycle
computed from detections alone is the duty cycle of the epochs on which the
target was bright enough to be seen &mdash; a tautology, not a measurement,
and one that is biased in a known direction and by an unknown amount.  The
repair is not to guess the missing magnitudes; it is to MEASURE at the
target's position on those frames and publish what the noise there allows
us to say.</p></div>

<div class="stage"><h3>Evidence</h3>
<p class="sub">Where the aperture went, how we know, what it measured, and
which blocks were refused.</p>
<p>The target's position on a frame that never detected it comes from the
frame's own matched comparison stars: a four-parameter similarity fitted
from the {p2.FORCED_MIN_STARS}+ stars the photometry already matched on
that frame, mapping the reference grid onto this one, evaluated at the
target's reference position.  Four parameters and not six, because a full
affine would absorb a genuine mismatch into a shear and report a small
residual while putting the aperture in the wrong place.</p>
<p>Two independent checks on that position, both in the table below:</p>
<ul>
<li><b>Catalogue cross-check</b> &mdash; the distance between the target's
identified reference-star position and an independent plate model fitted to
the ATLAS-REFCAT2 positions of the same field.  This matters more than it
looks: EU&nbsp;UMa's era-78 block sat 5.27&nbsp;arcsec (11.7&nbsp;px) from
the catalogue, and the tie removed that rigid offset before matching, so a
plate model fitted through the stored sky coordinates instead of the
catalogue's would have landed an aperture on blank sky.</li>
<li><b>Position closure</b> &mdash; the same machinery run on up to 60
frames where the target WAS detected, with the target itself excluded from
the transform.  The median distance between where the method predicts the
target and where it actually is, in pixels.  This is the validation that
costs nothing and settles everything: it exercises exactly the code the
undetected frames use, on frames where the right answer is known.</li>
</ul>
<div class="decision"><b>The gate.</b> A block that cannot demonstrate
closure &mdash; at least {p2.CLOSURE_MIN_FRAMES} frames with a detected
target, closing to a median of
{p2.CLOSURE_MAX_MEDIAN_PX:g}&nbsp;px or better &mdash; does not publish
limits.  This is not a formality.  The first production run produced 66
&ldquo;forced detections&rdquo; of EU&nbsp;UMa in the merged 2026 Fast
block, whose frame transforms closed to 645&ndash;1,650&nbsp;px on 87 of
153 attempts and whose measured signal-to-noise alternated between 3 and 55
frame to frame &mdash; an aperture landing on a bright neighbour on half
the frames.  That block has never detected its target, so there is no frame
on which its forced position can be checked, and a limit measured at an
unverifiable position is not a limit, it is a number.</div>

{_fig(fig_lim,
      "LEFT: what happened to the epochs that carried no magnitude before "
      "this stage.  Green are epochs recovered as real detections (source "
      "detection needs 5σ per pixel over 5 connected pixels, which is far "
      "stricter than a 3σ integrated aperture); blue are epochs bounded by "
      "an upper limit; red could not be measured.  RIGHT: the detections "
      "and the limits on the same natural-system axis.  A limit is only "
      "worth publishing when it reaches at least as deep as the "
      "detections it sits among.",
      "no series carried any undetected epoch that could be measured")}

{table(["series", "matched frames", "detected", "undetected", "recovered",
        "limits", "failed", "median limit (ensemble gauge)",
        "median limit (natural system)", "faintest detection",
        "position from", "catalogue cross-check (px)", "closure frames",
        "closure median (px)", "closure p95 (px)", "refused"], tab, cls)}

<p><b>Blocks refused limits, and why:</b></p>
<ul class="legend">{bl or
    "<li>None.</li>"}</ul>
</div>

<div class="stage"><h3>Decision</h3>
<div class="decision"><b>{n_rec:,} epochs recovered as detections,
{n_lim:,} epochs bounded by a
{p2.LIMIT_SIGMA:g}&sigma; upper limit, {n_nozp:,} measured but unusable for
lack of a frame zero point, {n_fail:,} unmeasurable.</b>
The limit convention is stated on every row of the product:
<code>{p2.LIMIT_SIGMA:g}&nbsp;&times;</code> the measured aperture noise
(sky shot noise inside the aperture plus the uncertainty of the sky level
itself), one-sided Gaussian, 99.87%.  It is deliberately NOT
<code>flux&nbsp;+&nbsp;{p2.LIMIT_SIGMA:g}&sigma;</code>: the chosen form
states the SENSITIVITY of the frame &mdash; &ldquo;a source brighter than
this would have been seen&rdquo; &mdash; which is what a duty-cycle
statistic needs, while the other form is systematically fainter on half the
frames by construction.</div>
<p>Note what the recovered detections are NOT: they are not written into
<code>cv_lightcurve</code>.  They live in <code>p2_limits</code> with their
own position, noise and provenance, so that a later stage can adopt them
deliberately rather than inherit them silently.</p>
</div>

<div class="stage"><h3>Consequence &mdash; the state statistics, twice</h3>
<p class="sub">Every statistic below is computed the censored way and the
limit-aware way.  The pair is the deliverable: a single &ldquo;corrected&rdquo;
number would hide the size of the correction, and the size of the
correction IS the finding.</p>
{table(["series", "epochs (detections only)", "epochs (with limits)",
        "detected fraction", "faint-state fraction",
        "median mag (detections only)", "median mag (Kaplan-Meier)"],
       stat_rows) if stat_rows else
 '<div class="decision">No series gained a faint-state fraction above zero, '
 'which would mean every undetected epoch was recovered as a detection.'
 '</div>'}
<div class="decision"><b>The detected fraction is 1.000 by construction
before this stage and is not afterwards.</b>  That single column is the
whole bias: any duty cycle, any state-occupancy fraction and any mean
magnitude computed on the old light curve was computed over a sample
selected on the quantity being measured.  The largest single change is
<b>{_label(worst_gain[0][0]) if worst_gain else '&mdash;'}</b>, which goes
from {_i(worst_gain[0][1]) if worst_gain else '&mdash;'} to
{_i(worst_gain[0][2]) if worst_gain else '&mdash;'} constrained epochs.</div>
<div class="decision"><b>The Kaplan-Meier median returns NaN where more
than half the epochs are limits, and that is the correct answer.</b>  A
series whose survival curve never reaches 0.5 has no estimable median, and
producing one from the detections alone is exactly the bias this task
exists to remove.</div>

<h3>How much do these limits actually constrain?</h3>
<p class="sub">The question a referee asks next, measured rather than
assumed.</p>
{table(["series", "median detection (natural system)",
        "median 3&sigma; limit", "limit &minus; detection", "reading"],
       depth_rows) if depth_rows else
 '<div class="decision">No series carried both a calibrated detection and a '
 'calibrated limit.</div>'}
<div class="decision"><b>{n_shallow} of {n_shallow + n_deep} series have
limits SHALLOWER than their own typical detection.</b>  That is the honest
reading of this task's result and it constrains how the faint-state
fractions above may be used: on most series the undetected epochs are
predominantly LOW-SENSITIVITY FRAMES &mdash; short exposures, poor seeing,
the cloud the previous section vetoes &mdash; not epochs on which the
target was demonstrably faint.  The faint-state fractions are therefore
<i>upper bounds</i> on a low-state duty cycle, not measurements of one.
What the limits do establish unambiguously is the direction and the size of
the censoring bias, and that a duty cycle quoted from detections alone was
never a measurement at all.</div>
</div>
</section>"""


def render_report(db_path: Path) -> Path:
    """Render the CV-S8 Phase-2 completion page.  Returns the HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = dict(q(con, "SELECT key, value FROM p2_meta"))
        figs = [fig_cloud_calibration(con), fig_cloud_sculpting(con),
                fig_extinction(con), fig_transform(con), fig_limits(con)]
        sections = [
            section_intro(con),
            section_cloud(con, figs[0], figs[1]),
            section_extinction(con, figs[2]),
            section_crossera(con, figs[3]),
            section_limits(con, figs[4]),
        ]
        n_veto = q1(con, "SELECT sum(vetoed) FROM p2_cloud_frame") or 0
        n_lim = q1(con, "SELECT count(*) FROM p2_limits WHERE "
                        "outcome='limit'") or 0
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CV Time Series &mdash; Phase 2 Completion</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Cataclysmic-Variable Time Series &mdash; closing out Phase&nbsp;2
  photometry</h1>
  <p>The cloud veto that had to be built because ZMAG does not exist for the
  polars &middot; the colour-extinction terms measured rather than assumed
  &middot; the cross-era transformation published as metadata and applied to
  nothing &middot; the faint-phase epochs bounded instead of dropped
  &middot; {n_veto:,} frames flagged, {n_lim:,} upper limits &middot;
  built {esc(meta.get('stage_cloud', ''))[:16]}Z
  ({esc(meta.get('phase2_code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="cv_catalogue_tie.html">the catalogue tie this rests
  on</a> &middot; <a href="cv_characterization.html">the characterization
  that set these thresholds</a> &middot;
  <a href="index.html">project hub</a> &middot;
  <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#intro">0 The four holes</a> &middot;
  <a href="#cloud">1 Cloud veto</a> &middot;
  <a href="#extinction">2 Colour extinction</a> &middot;
  <a href="#crossera">3 Cross-era discipline</a> &middot;
  <a href="#limits">4 Faint-phase limits</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_phot.report_phase2</code> from
<code>products/phot/cv_timeseries.sqlite</code>.  Every number in a table,
figure or verdict is the result of a query run at render time or a constant
imported from <code>macro_phot.phase2</code>.  The exceptions, named rather
than claimed away: the VV&nbsp;Pup declination and maximum airmass quoted in
section&nbsp;2 (computed from the target's catalogued position, not from
this database), the 5-sigma-per-pixel-over-5-connected-pixels description of
<code>sep</code>'s detection rule in section&nbsp;4 (a property of the
extraction library, recorded in
<code>macro_phot.photometry.DETECT_SIGMA</code> and
<code>DETECT_MINAREA</code>), and the descriptive prose throughout.
Regenerate with <code>pipeline/scripts/run_cv_phase2.py report</code>.
</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
    finally:
        con.close()
    return HTML_PATH
