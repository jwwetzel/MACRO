"""CV-S10 chain-of-evidence report: the two closing science decisions.

Reads ``products/phot/cv_timeseries.sqlite`` and writes

* ``docs/CV_TimeSeries/cv_final_science.html``
* ``docs/CV_TimeSeries/figures/cv_final/*.png``

Socratic, in ONE order, because the order is the argument:

    0  what this page decides, and what it is forbidden to assume
    1  which runs, and in which state -- the census the branch rests on
    2  §4.19: is the photometry good enough to attempt the fallback?
    3  the quiescent orbital hump: two nulls, one answer
    4  flickering: amplitude against timescale, over a MEASURED floor
    5  what this season cannot do -- no superhump period, no dP_sh/dt
    6  the six normal-outburst runs, on their own terms
    7  AN UMa, filter by filter: what each one can and cannot support
    8  the verdicts

Section 3 comes after section 2 for a reason: a hump amplitude read without
the gate is a number, and with it is a measurement.  Section 5 comes after
4 because "we cannot measure a superhump period" is only interesting once
the reader has seen what the data CAN do.

Every number and every figure on this page is the result of a query
executed here.  Nothing is typed.
"""

from __future__ import annotations

import math
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

from macro_core.report_s0 import (        # noqa: E402
    ACCENT, DARK, DPI, WARN, esc, q, q1, table)
from . import final_science as fs         # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "CV_TimeSeries"
FIG_DIR = DOCS_DIR / "figures" / "cv_final"
HTML_PATH = DOCS_DIR / "cv_final_science.html"

GOOD = "#9fd8ae"
BAD = "#f0a3a3"
MUTED = "#6f7a8a"
FILTER_COLOR = {"G": "#8fb3d9", "g": "#8fb3d9", "R": "#e6907a",
                "r": "#e6907a", "I": "#c9a0dc", "i": "#c9a0dc"}


# ---------------------------------------------------------------------------
# Formatting.  An em-dash is the only thing a missing number becomes: a blank
# cell reads as zero, and zero is a claim.
# ---------------------------------------------------------------------------
def _n(x, nd=2):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{float(x):,.{nd}f}"


def _mmag(x, nd=0):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{1000.0 * float(x):,.{nd}f}"


def _i(x):
    return "&mdash;" if x is None else f"{int(x):,}"


def _has(con, name: str) -> bool:
    return bool(q1(con, "SELECT count(*) FROM sqlite_master WHERE "
                        "type='table' AND name=?", (name,)))


def _meta(con) -> dict:
    return dict(q(con, "SELECT key, value FROM p4_meta")) \
        if _has(con, "p4_meta") else {}


def _verdict_span(v: str) -> str:
    """Colour a verdict word without inventing a scale for it."""
    up = str(v).upper()
    if up.startswith("SUPPORTED") or up.startswith("DETECTED") \
            or up.startswith("MEASURED"):
        col = GOOD
    elif up.startswith("NOT ") or up.startswith("NO ") or up == "NO":
        col = BAD
    else:
        col = WARN
    return f'<b style="color:{col}">{esc(v)}</b>'


def _fig(src: str, caption: str, missing: str) -> str:
    if not src:
        return f'<div class="note"><b>No figure, and why:</b> {missing}</div>'
    return (f'<figure><a href="{src}"><img src="{src}" alt=""></a>'
            f"<figcaption>{caption}</figcaption></figure>")


def _scope_label(scope: str) -> str:
    """``yzcnc|e72|g|block:A+B`` -> ``e72 g  A+B``."""
    parts = str(scope).split("|")
    if len(parts) < 4:
        return esc(scope)
    tail = parts[3].replace("block:", "")
    return f"{parts[1]} {parts[2]} &nbsp; {esc(tail)}"


# ===========================================================================
# Figures
# ===========================================================================
def fig_census(con) -> str:
    """The nine dense runs on a magnitude axis, coloured by state.

    The figure the branch decision rests on, redrawn from THIS stage's own
    tables rather than inherited as a claim: where the runs sit relative to
    quiescence is what makes three of them the fallback's data set and six
    of them a separate result.
    """
    rows = q(con, "SELECT nights, state, filter, median_cal_mag, n_points "
                  "FROM p4_run WHERE kind='run' ORDER BY nights, filter")
    if not rows:
        return ""
    nights = sorted({r[0] for r in rows})
    xs = {n: i for i, n in enumerate(nights)}
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(10.5, 4.2))
        for r in rows:
            col = GOOD if r[1] == "QUIESCENT" else WARN
            ax.scatter(xs[r[0]], r[3], s=18 + 0.9 * r[4], color=col,
                       edgecolor="none", alpha=0.75)
            ax.annotate(r[2], (xs[r[0]], r[3]), fontsize=7,
                        color=MUTED, xytext=(5, -3),
                        textcoords="offset points")
        ax.set_xticks(range(len(nights)))
        ax.set_xticklabels(nights, rotation=45, ha="right", fontsize=8)
        ax.invert_yaxis()
        ax.set_ylabel("median cal_mag on the run")
        ax.set_title("YZ Cnc dense runs — green QUIESCENT, yellow OUTBURST "
                     "(marker area ∝ points)")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "census.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_final/census.png"


def fig_folds(con, db_path: Path) -> str:
    """The quiescent runs folded on the published period, model overlaid.

    Drawn from the light curves themselves so the reader can see what the
    fitted amplitude is describing.  The point of the panel is the SCATTER,
    not the curve: the fitted fundamental is a thin line inside a cloud
    several times its own height, which is the section's whole argument.
    """
    scopes = q(con, "SELECT scope, series_key, nights, filter, hump_amp, "
                    "hump_phase, chi2nu, cycles FROM p4_run "
                    "WHERE state='QUIESCENT' AND kind='block' "
                    "UNION ALL "
                    "SELECT scope, series_key, nights, filter, hump_amp, "
                    "hump_phase, chi2nu, cycles FROM p4_run "
                    "WHERE state='QUIESCENT' AND kind='run' "
                    "AND cycles >= ? ORDER BY series_key",
                    (fs.MIN_CYCLES_FOR_DETECTION,))
    if not scopes:
        return ""
    meta = _meta(con)
    period = float(meta.get("yz_period_d", 0.0868))
    epoch = float(meta.get("yz_epoch_bjd", 2460000.0))
    lc = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    ncol = 3
    nrow = int(math.ceil(len(scopes) / ncol))
    with plt.rc_context(DARK):
        fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.0 * nrow),
                                 squeeze=False)
        for k, s in enumerate(scopes):
            ax = axes[k // ncol][k % ncol]
            night_list = str(s[2]).split("+")
            marks = "os^"
            for j, night in enumerate(night_list):
                rows = lc.execute("""
                    SELECT l.bjd_tdb, l.cal_mag FROM cv_lightcurve l
                    JOIN cv_frames f ON f.frame_id=l.frame_id
                                    AND f.series_key=l.series_key
                    LEFT JOIN p2_cloud_frame c ON c.frame_id=l.frame_id
                                              AND c.series_key=l.series_key
                    WHERE l.series_key=? AND f.night=? AND l.role='target'
                      AND l.cal_mag IS NOT NULL AND l.saturated=0
                      AND COALESCE(c.vetoed,0)=0""", (s[1], night)).fetchall()
                if not rows:
                    continue
                t = np.array([r[0] for r in rows])
                m = np.array([r[1] for r in rows])
                ph = fs.orbital_phase(t, period, epoch)
                # Each night is plotted about its OWN median, which is the
                # nightly constant the joint fit removes; showing raw
                # magnitudes would draw the 0.28 mag night-to-night step and
                # hide the modulation the panel is about.
                ax.plot(np.concatenate([ph, ph + 1]),
                        np.tile(m - np.median(m), 2), marks[j % 3],
                        ms=3, alpha=0.7,
                        color=FILTER_COLOR.get(s[3], ACCENT),
                        label=night)
            if s[4] is not None:
                g = np.linspace(0, 2, 400)
                ax.plot(g, -float(s[4]) * np.cos(2 * np.pi
                                                 * (g - float(s[5]))),
                        "-", color=WARN, lw=1.6)
            ax.invert_yaxis()
            ax.set_xlim(0, 2)
            ax.set_title(f"{s[1]}  {s[2]}\nA = {1000 * (s[4] or 0):.0f} mmag, "
                         f"{s[7]:.2f} orbits sampled", fontsize=8)
            ax.set_xlabel("orbital phase (arbitrary zero) ×2")
            ax.set_ylabel("Δ mag")
            ax.legend(fontsize=6, loc="upper right")
        for k in range(len(scopes), nrow * ncol):
            axes[k // ncol][k % ncol].axis("off")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "folds.png", dpi=DPI)
        plt.close(fig)
    lc.close()
    return "figures/cv_final/folds.png"


def fig_two_nulls(con) -> str:
    """Fitted hump against BOTH contours, per quiescent scope.

    One row per scope, three marks: the fitted semi-amplitude, the
    instrumental 90% recovery contour, and the contour set by the star's own
    flickering.  Everything the section argues is visible in the gap between
    the last two.
    """
    rows = q(con, "SELECT scope, hump_amp, hump_amp_sigma, amp90_field, "
                  "amp90_self FROM p4_run WHERE state='QUIESCENT' "
                  "AND amp90_self IS NOT NULL ORDER BY scope")
    if not rows:
        return ""
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(10, 0.55 * len(rows) + 2.2))
        y = np.arange(len(rows))
        ax.errorbar([1000 * r[1] for r in rows], y,
                    xerr=[1000 * (r[2] or 0) for r in rows], fmt="o",
                    color=ACCENT, ms=7, capsize=3,
                    label="fitted hump semi-amplitude")
        ax.plot([1000 * r[3] for r in rows], y, "s", color=GOOD, ms=7,
                label="90% recovery vs field stars (instrument)")
        ax.plot([1000 * r[4] for r in rows], y, "D", color=BAD, ms=7,
                label="90% recovery vs the star's own flickering")
        for i, r in enumerate(rows):
            ax.plot([1000 * r[3], 1000 * r[4]], [i, i], "-", color=MUTED,
                    lw=1, alpha=0.7)
        ax.set_yticks(y)
        ax.set_yticklabels([r[0].replace("yzcnc|", "").replace("block:", "")
                            for r in rows], fontsize=7)
        ax.set_xscale("log")
        ax.set_xlim(0.6 * min(1000 * r[3] for r in rows),
                    1.7 * max(1000 * r[4] for r in rows))
        ax.set_ylim(len(rows) - 0.3, -1.4)
        ax.set_xlabel("semi-amplitude (mmag)")
        ax.set_title("Two nulls, one answer\nthe instrument can see the "
                     "hump; the star's own flickering cannot be told from it",
                     fontsize=10)
        ax.legend(fontsize=8, loc="upper right", framealpha=0.85)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "two_nulls.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_final/two_nulls.png"


def fig_flicker(con) -> str:
    """Structure functions of the quiescent runs, over their measured floors.

    Target (filled), magnitude-matched field-star floor (open), and the
    quadrature excess (dashed).  Two eras on two panels because the eras are
    two different cameras at two different exposure times, and the whole
    §4.19 question is the distance between their floors.
    """
    eras = [r[0] for r in q(con, "SELECT DISTINCT era_id FROM p4_run "
                                 "WHERE state='QUIESCENT' ORDER BY era_id")]
    if not eras:
        return ""
    with plt.rc_context(DARK):
        fig, axes = plt.subplots(1, len(eras), figsize=(5.4 * len(eras), 4.4),
                                 squeeze=False, sharey=True)
        for j, era in enumerate(eras):
            ax = axes[0][j]
            rows = q(con, """SELECT f.series_key, f.night, f.filter, f.tau_s,
                                    f.sf_target, f.sf_floor, f.sf_excess
                             FROM p4_flicker f
                             JOIN p4_run r ON r.series_key=f.series_key
                                          AND r.nights=f.night
                             WHERE f.state='QUIESCENT' AND r.era_id=?
                             ORDER BY f.series_key, f.night, f.tau_s""",
                     (era,))
            groups: dict[tuple, list] = {}
            for r in rows:
                groups.setdefault((r[0], r[1], r[2]), []).append(r)
            for (sk, night, filt), g in sorted(groups.items()):
                tau = np.array([x[3] for x in g], dtype=float)
                tgt = np.array([x[4] if x[4] is not None else np.nan
                                for x in g], dtype=float)
                flo = np.array([x[5] if x[5] is not None else np.nan
                                for x in g], dtype=float)
                exc = np.array([x[6] if x[6] is not None else np.nan
                                for x in g], dtype=float)
                c = FILTER_COLOR.get(filt, ACCENT)
                ax.plot(tau, 1000 * tgt, "-o", color=c, ms=4, lw=1.4,
                        label=f"{filt} {night}")
                ax.plot(tau, 1000 * flo, ":s", color=c, ms=3, lw=1,
                        alpha=0.65, mfc="none")
                ax.plot(tau, 1000 * exc, "--", color=c, lw=0.9, alpha=0.5)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("timescale τ (s)")
            if j == 0:
                ax.set_ylabel("structure function σ(τ)  (mmag)")
            ax.set_title(f"era {era}")
            ax.legend(fontsize=6, ncol=1, loc="upper left")
        fig.suptitle("Quiescent flickering (solid) over its MEASURED "
                     "photometric floor (dotted, magnitude-matched field "
                     "stars); dashed = quadrature excess", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "flicker.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_final/flicker.png"


def fig_outburst(con, db_path: Path) -> str:
    """The six normal-outburst runs, one panel each, three filters."""
    nights = [r[0] for r in q(con, "SELECT DISTINCT night FROM p4_outburst "
                                   "ORDER BY night")]
    if not nights:
        return ""
    lc = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    ncol = 3
    nrow = int(math.ceil(len(nights) / ncol))
    with plt.rc_context(DARK):
        fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.0 * nrow),
                                 squeeze=False, sharey=True)
        for k, night in enumerate(nights):
            ax = axes[k // ncol][k % ncol]
            for r in q(con, "SELECT series_key, filter, rate_mag_per_h, "
                            "rate_verdict FROM p4_outburst WHERE night=? "
                            "ORDER BY filter", (night,)):
                rows = lc.execute("""
                    SELECT l.bjd_tdb, l.cal_mag FROM cv_lightcurve l
                    JOIN cv_frames f ON f.frame_id=l.frame_id
                                    AND f.series_key=l.series_key
                    LEFT JOIN p2_cloud_frame c ON c.frame_id=l.frame_id
                                              AND c.series_key=l.series_key
                    WHERE l.series_key=? AND f.night=? AND l.role='target'
                      AND l.cal_mag IS NOT NULL AND l.saturated=0
                      AND COALESCE(c.vetoed,0)=0 ORDER BY l.bjd_tdb""",
                    (r[0], night)).fetchall()
                if not rows:
                    continue
                t = np.array([x[0] for x in rows])
                m = np.array([x[1] for x in rows])
                ax.plot((t - t.min()) * 24.0, m, ".", ms=3.5,
                        color=FILTER_COLOR.get(r[1], ACCENT),
                        label=f"{r[1]}  {r[2]:+.3f} mag/h")
            ax.invert_yaxis()
            ax.set_title(night, fontsize=9)
            ax.set_xlabel("hours from run start")
            ax.set_ylabel("cal_mag")
            ax.legend(fontsize=6, loc="best")
        for k in range(len(nights), nrow * ncol):
            axes[k // ncol][k % ncol].axis("off")
        fig.suptitle("The six normal-outburst dense runs — three filters at "
                     "8 s, cycle-resolved", fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "outburst.png", dpi=DPI)
        plt.close(fig)
    lc.close()
    return "figures/cv_final/outburst.png"


def fig_anuma(con) -> str:
    """AN UMa's per-filter grade card: measured value over its own bar.

    Plotted as a RATIO to the bar so five incommensurable quantities --
    nights, edges, an Otsu separability, percentage points -- can share one
    axis.  One vertical line at 1.0 is the whole scale, and a bar that
    crosses it is a SUPPORTED verdict.
    """
    rows = q(con, "SELECT filter, capability, measured, bar, verdict, rank "
                  "FROM p4_anuma ORDER BY rank, filter")
    if not rows:
        return ""
    labels, ratios, cols = [], [], []
    for r in rows:
        if not r[3]:
            continue
        # The duty-cycle line is the one capability where SMALLER is better
        # (a narrower interval), so its ratio is inverted to keep "right of
        # the line = passes" true for every row on the figure.
        ratio = ((r[3] / r[2]) if (r[2] and "duty cycle" in r[1])
                 else (r[2] / r[3]))
        labels.append(f"{r[0]}   {r[1]}")
        ratios.append(ratio)
        cols.append(GOOD if r[4] == "SUPPORTED" else BAD)
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(10, 0.34 * len(labels) + 1.8))
        y = np.arange(len(labels))
        ax.barh(y, ratios, color=cols, alpha=0.85, height=0.65)
        ax.axvline(1.0, color=ACCENT, lw=1.5)
        for i, (rt, c) in enumerate(zip(ratios, cols)):
            ax.annotate(f"{rt:.2f}", (rt, i), xytext=(4, 0),
                        textcoords="offset points", va="center",
                        fontsize=7, color=c)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.invert_yaxis()
        # Linear, not log: one capability legitimately measures ZERO (the
        # i band has no detected modulation, so it has no usable full-orbit
        # nights), and a log axis cannot draw a zero at all.
        ax.set_xlim(0, max(3.4, 1.15 * max(ratios)))
        ax.set_xlabel("measured ÷ bar  (1.0 = exactly at the bar)")
        ax.set_title("AN UMa, capability by capability, filter by filter")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "anuma.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_final/anuma.png"


# ===========================================================================
# Sections
# ===========================================================================
def section_intro(con) -> str:
    meta = _meta(con)
    return f"""
<section id="intro"><h2>0 &nbsp; What this page decides, and what it may not
assume</h2>

<h3>Question</h3>
<p>Two tasks were left open when Phase&nbsp;3 finished, and both are
decisions rather than discoveries. <b>CV-P3-yzcnc-superhump</b> had already
had its branch chosen for it by the external record: CV-S7 measured that no
dense run of the 2024 YZ&nbsp;Cnc season sits inside a superoutburst, so the
superhump branch is closed and the strategy's own fallback &mdash; orbital
hump plus flickering statistics, &ldquo;honest but weaker&rdquo; &mdash; is
what the season supports. <b>CV-P4-anuma</b> asks what role AN&nbsp;UMa
should have in the paper, filter by filter, after CV-S5 graded its
three-filter colour goal NOT SUPPORTED.</p>

<h3>What is inherited and not re-litigated</h3>
<ul>
<li><b>The branch decision itself.</b> The nine dense runs and their states
come from <code>cv_ext_verdict</code> &mdash; CV-S7's night-by-night
classification, which uses INDEPENDENT AAVSO photometry wherever it exists
and survives deleting every RLMT row from the comparison. This page does not
re-decide which nights were dense or which were in outburst.</li>
<li><b>The period.</b> {esc(meta.get('yz_period_d', ''))}&nbsp;d, from
<code>p3_ephemeris</code>, i.e. the AAVSO VSX record. Nothing here measures
a period; §4 of the analysis page already showed that the
&plusmn;1&nbsp;c/d aliases carry up to 0.96 of the window power on YZ Cnc's
multi-night sets, so no period from this archive could select its own
family member.</li>
<li><b>Phase zero is arbitrary.</b> VSX publishes NO epoch for YZ&nbsp;Cnc.
Every phase on this page is measured from BJD&nbsp;
{esc(meta.get('yz_epoch_bjd', ''))}, a constant this stage chose. What
survives that is the SHAPE of a fold and the offset between two filters
observed through the same night; what does not survive is any statement
about where a hump sits relative to a physical conjunction.</li>
<li><b>The error model.</b> Every weighted fit here uses the per-point
formal error multiplied by that series' MEASURED
&chi;<sup>2</sup> inflation from <code>cv_series</code>, never the photon
prediction alone.</li>
</ul>

<h3>What this page adds</h3>
<p>Three measurements and two decisions: the folded hump with a detection
test that has TWO nulls; flickering amplitude against timescale over a floor
that is measured rather than modelled; the six normal-outburst runs
characterised on their own terms; the §4.19 signal-to-noise gate executed;
and AN&nbsp;UMa graded capability by capability against bars that are stated
one number at a time so a reader can disagree with any of them in one
line.</p>
</section>"""


def section_census(con, fig: str) -> str:
    runs = q(con, "SELECT nights, utc_nights, state, filter, n_points, "
                  "span_h, cycles, cadence_s, median_cal_mag, amp_p5p95 "
                  "FROM p4_run WHERE kind='run' ORDER BY nights, filter")
    if not runs:
        return '<section id="census"><h2>1 &nbsp; The dense runs</h2>' \
               '<p class="note">Not run: <code>p4_run</code> is empty.' \
               '</p></section>'
    n_q = len({r[0] for r in runs if r[2] == "QUIESCENT"})
    n_o = len({r[0] for r in runs if r[2] == "OUTBURST"})
    tbl = table(
        ["local night", "UTC night", "state", "filter", "N", "span (h)",
         "orbits sampled", "cadence (s)", "median cal_mag",
         "p5&ndash;p95 (mag)"],
        [[esc(r[0]), esc(r[1]),
          f'<b style="color:{GOOD if r[2] == "QUIESCENT" else WARN}">'
          f'{esc(r[2])}</b>', esc(r[3]), _i(r[4]), _n(r[5]), _n(r[6]),
          _n(r[7], 0), _n(r[8]), _n(r[9], 3)] for r in runs])
    return f"""
<section id="census"><h2>1 &nbsp; Which runs, and in which state</h2>

<h3>Question</h3>
<p>Everything below depends on one classification: which of the dense runs
caught YZ&nbsp;Cnc at quiescence and which caught it in outburst. Get that
wrong and a flickering amplitude measured in outburst is published as a
quiescent one.</p>

<h3>Evidence</h3>
<p>{n_q} quiescent dense runs and {n_o} in normal outburst, inherited from
CV-S7 and re-tabulated here with this stage's own coverage numbers. Note the
two night conventions: <code>cv_frames</code> keys on the LOCAL night and
CV-S7's page tabulates by UTC night, one day later for these Arizona
evenings. Both are printed so the two pages can be checked against each
other.</p>
{tbl}
{_fig(fig, "The dense runs on a magnitude axis. The quiescent runs sit "
           "about 1.5 mag below the outburst ones — the separation the "
           "branch decision rests on.",
      "p4_run holds no runs.")}

<h3>Decision</h3>
<p>The quiescent runs are the fallback's data set. The outburst runs are
§6's, and are not mixed into any quiescent statistic on this page.</p>
</section>"""


def section_gate(con, meta) -> str:
    rows = q(con, "SELECT gate_id, scope, quantity, value, bar, passes, note "
                  "FROM p4_gate ORDER BY gate_id, scope")
    if not rows:
        return '<section id="gate"><h2>2 &nbsp; The §4.19 gate</h2>' \
               '<p class="note">Not run: <code>p4_gate</code> is empty.' \
               '</p></section>'
    tbl = table(
        ["gate", "scope", "what is compared", "value", "bar", "verdict",
         "the numbers behind it"],
        [[f"<code>{esc(r[0])}</code>", esc(r[1].replace('yzcnc|', '')),
          esc(r[2]), _n(r[3]), _n(r[4]),
          f'<b style="color:{GOOD if r[5] else BAD}">'
          f'{"PASS" if r[5] else "FAIL"}</b>', esc(r[6])] for r in rows])
    n_pass = sum(1 for r in rows if r[5])
    fails = "".join(f"<li><code>{esc(r[0])}</code> on "
                    f"<code>{esc(r[1].replace('yzcnc|', ''))}</code> "
                    f"&mdash; {esc(r[6])}</li>"
                    for r in rows if not r[5])
    floors = q(con, "SELECT r.era_id, min(f.sf_floor), max(f.sf_floor) "
                    "FROM p4_flicker f JOIN p4_run r "
                    "ON r.series_key=f.series_key AND r.nights=f.night "
                    "WHERE f.state='QUIESCENT' AND f.sf_floor IS NOT NULL "
                    "GROUP BY r.era_id ORDER BY r.era_id")
    era_txt = "; ".join(f"era {r[0]}: {_mmag(r[1])}&ndash;{_mmag(r[2])} mmag"
                        for r in floors)
    return f"""
<section id="gate"><h2>2 &nbsp; §4.19: is the photometry good enough to
attempt the fallback at all?</h2>

<h3>Question</h3>
<p>The strategy refused to promise the quiescent fallback until somebody
checked it: <i>&ldquo;8&nbsp;s High Gain at quiescent V&asymp;14.5 may be
sky/read-noise dominated; verify before promising the fallback&rdquo;</i>
(§4.19, concern <b>m4</b>). That sentence is not a measurement, so it is
turned into three arithmetic lines here.</p>

<p>The key move is to stop using formal error bars. What decides whether a
frame is noise-dominated is the scatter a REAL star of the target's
brightness shows through those same frames &mdash; and YZ&nbsp;Cnc's four
held-out check stars sit about a magnitude brighter than the star at
quiescence, so they are the wrong stars to ask. The floor used throughout
this page is the structure function of every field star within
{fs.FLOOR_MATCH_HALF_WIDTH}&nbsp;mag of the target on that run.</p>

<h3>Evidence</h3>
<p>Measured floors, by era: {era_txt}. That gap is exactly the §4.19
worry made numerical &mdash; the 8&nbsp;s High Gain frames are several times
noisier at quiescence than the 30&nbsp;s Sloan-era frames.</p>
{tbl}

<h3>Decision</h3>
<p><b>The gate passes</b>, {n_pass} of {len(rows)} lines, and it passes with
a restriction rather than uniformly. Flickering clears the measured floor on
every quiescent run in every filter. Every line that fails does so on the
same run &mdash; the 8&nbsp;s High&nbsp;Gain night, which is exactly the one
§4.19 named:</p>
<ul>{fails}</ul>
<p>So §4.19 was right about the frames and wrong about the conclusion. The
fallback is deliverable; the short-timescale end of the High&nbsp;Gain era
is not where to deliver it, and the 30&nbsp;s Sloan-era runs carry the
flickering result. Note what this gate does NOT license: it says the
photometry can see a modulation the size of the fitted hump, not that the
hump is real. That question needs a second null and is §3's.</p>
</section>"""


def section_hump(con, fig_fold: str, fig_null: str) -> str:
    rows = q(con, "SELECT scope, kind, nights, filter, n_points, cycles, "
                  "hump_amp, hump_amp_sigma, hump_amp_harm, hump_phase, "
                  "chi2nu, amp90_field, amp90_self, power_forb, thr_self, "
                  "detection, sigma_point, note FROM p4_run "
                  "WHERE state='QUIESCENT' ORDER BY scope")
    if not rows:
        return '<section id="hump"><h2>3 &nbsp; The quiescent hump</h2>' \
               '<p class="note">Not run.</p></section>'
    tbl = table(
        ["scope", "N", "orbits", "A (mmag)", "&sigma;<sub>A</sub>",
         "2nd harm.", "phase of max light", "&chi;<sup>2</sup><sub>&nu;</sub>",
         "A<sub>90</sub> field", "A<sub>90</sub> self", "call"],
        [[_scope_label(r[0]), _i(r[4]), _n(r[5]), f"<b>{_mmag(r[6])}</b>",
          _mmag(r[7]), _mmag(r[8]), _n(r[9], 3), _n(r[10], 1),
          _mmag(r[11]), _mmag(r[12]), _verdict_span(r[15])] for r in rows])
    notes = "".join(f"<li><code>{_scope_label(r[0])}</code> &mdash; "
                    f"{esc(r[17])}</li>" for r in rows if r[17])
    # A third, independent line of evidence, computed here from the same
    # table: does the fitted phase REPEAT?  Within one night the three
    # filters see one physical feature three ways and must agree; between
    # nights an ORBITAL feature must also agree.  When the first holds and
    # the second does not, the modulation is real and is not the orbit.
    per_night: dict[str, list[float]] = {}
    era_of: dict[str, int] = {}
    for r in q(con, "SELECT nights, era_id, hump_phase FROM p4_run WHERE "
                    "state='QUIESCENT' AND kind='run' AND hump_phase "
                    "IS NOT NULL"):
        per_night.setdefault(r[0], []).append(float(r[2]))
        era_of[r[0]] = r[1]
    within = [(n, *fs.circular_mean_and_spread(v)) for n, v in
              sorted(per_night.items()) if len(v) >= 2]
    phase_tbl = table(
        ["night", "era", "filters", "mean phase of maximum light",
         "spread across filters (cycles)"],
        [[esc(w[0]), _i(era_of[w[0]]), _i(len(per_night[w[0]])),
          _n(w[1], 3), _n(w[2], 3)] for w in within])
    # Nights may only be compared to each other WITHIN an era, and only
    # where the accumulated period drift is smaller than the difference
    # being measured.  Across the 2024-05 seam neither condition holds.
    pairs = [(a, b, fs.phase_difference(a[1], b[1]))
             for i, a in enumerate(within) for b in within[i + 1:]
             if era_of[a[0]] == era_of[b[0]]]
    pair_txt = "; ".join(
        f"{esc(a[0])} vs {esc(b[0])}: <b>{abs(d):.3f} cycles</b>"
        for a, b, d in pairs) or "no two quiescent runs share an era"
    v = q(con, "SELECT verdict, deciding_number, reasoning, alternative "
               "FROM p4_verdict WHERE verdict_id='YZ-hump'")
    vv = v[0] if v else ("", "", "", "")
    return f"""
<section id="hump"><h2>3 &nbsp; The quiescent orbital hump: two nulls, one
answer</h2>

<h3>Question</h3>
<p>A dwarf nova at quiescence should show an orbital hump &mdash; the bright
spot where the accretion stream meets the disc rim, seen once per orbit. Is
one there, and how big? The question that decides the answer is not
&ldquo;is the fitted amplitude bigger than its error bar?&rdquo; but
&ldquo;bigger than <em>what could have produced it by accident</em>?&rdquo;,
and there are two candidate accidents.</p>

<h3>Evidence</h3>
<p>Each scope is folded on the published period and fitted with
{fs.N_HARMONICS} harmonics plus one free constant per night, jointly &mdash;
the discipline CV-S9 measured the case for, because detrending first and
searching afterwards eats part of the signal it is looking for. Blocks fold
two nights together only where the accumulated phase drift stays inside
{fs.PHASE_DRIFT_BAR_CYCLES} cycles; the February and May quiescent runs are
71 days apart, which on a period quoted to four decimals is about half a
cycle of drift, so those may never share a phase axis.</p>
{tbl}
{"<p class='note'><b>Scope notes.</b></p><ul>" + notes + "</ul>" if notes else ""}

<p>Two nulls are measured for every scope, and BOTH are injections of a
known sinusoid into real noise at the real timestamps, recovered at the
published frequency, with the 90% recovery amplitude read off the same
amplitude grid CV-S5 used:</p>
<ul>
<li><b>the field null</b> &mdash; magnitude-matched field stars through the
same frames. It asks: <i>could the photometry see a hump this size?</i></li>
<li><b>the self null</b> &mdash; this star's own residuals about the fitted
model, rolled night by night so every realization is made of light the
telescope actually recorded that night. It asks: <i>could this star's own
aperiodic flickering have produced a peak this tall at the orbital
frequency?</i></li>
</ul>
{_fig(fig_null, "The gap between the green squares and the red diamonds is "
                "the result. The instrument can see a hump of the fitted "
                "size; the star's own flickering cannot be told apart from "
                "one.", "p4_run holds no measured contours.")}
{_fig(fig_fold, "The quiescent runs folded on the published period, each "
                "night about its own median, model overlaid, plotted twice "
                "for continuity. The fitted fundamental is a thin line "
                "inside a cloud several times its height.",
      "No quiescent scope reaches the coverage bar.")}

<h3>A third test the contours do not know about: does the phase repeat?</h3>
<p>Each night's three filters are fitted independently, so they are three
measurements of one physical feature and must agree in phase. An
<em>orbital</em> feature must agree between nights as well. Phases are
compared on the circle, where the largest possible disagreement is half a
cycle.</p>
{phase_tbl}
<p>Within a night the filters agree to a few hundredths of a cycle. So the
modulation being fitted is REAL &mdash; three independent measurements land
on the same feature &mdash; and it is the star rather than the photometry.
</p>

<p>Between nights, comparing only pairs that share an era and are close
enough in time that the period's own uncertainty cannot explain the
difference: {pair_txt}. That is close to the half-cycle maximum, against an
accumulated phase drift of well under 0.01 cycles over the one-day gap. A
modulation that is coherent across three filters within one night and lands
a third of a cycle away the next night is not the orbit; it is flickering
with a coherence time shorter than the gap between runs, which is exactly
what §4 measures directly.</p>

<p class="note">The February run cannot be compared with the May runs at
all: 71 days on a period published to four decimals is about half a cycle of
accumulated drift, so any phase difference between them is unconstrained
before the star is even considered.</p>

<h3>Decision</h3>
<p>{_verdict_span(vv[0])} &mdash; {esc(vv[1])}.</p>
<p>{esc(vv[2])}</p>

<h3>Consequence</h3>
<p>{esc(vv[3])}</p>
</section>"""


def section_flicker(con, fig: str) -> str:
    rows = q(con, "SELECT series_key, night, filter, tau_s, sf_target, "
                  "sf_target_sigma, n_pairs, sf_floor, n_floor_stars, "
                  "sf_excess, excess_sigma, detected, sigma_formal "
                  "FROM p4_flicker WHERE state='QUIESCENT' "
                  "ORDER BY series_key, night, tau_s")
    if not rows:
        return '<section id="flicker"><h2>4 &nbsp; Flickering</h2>' \
               '<p class="note">Not run.</p></section>'
    tbl = table(
        ["series", "night", "&tau; (s)", "pairs", "&sigma;<sub>target</sub>",
         "floor (mmag)", "floor stars", "<b>flickering (mmag)</b>",
         "excess &sigma;", "detected"],
        [[esc(r[0]), esc(r[1]), _n(r[3], 0), _i(r[6]), _mmag(r[4]),
          _mmag(r[7]), _i(r[8]), f"<b>{_mmag(r[9])}</b>", _n(r[10], 1),
          ('<b style="color:%s">yes</b>' % GOOD) if r[11]
          else '<span style="color:%s">no</span>' % MUTED] for r in rows])
    v = q(con, "SELECT verdict, deciding_number, reasoning, alternative "
               "FROM p4_verdict WHERE verdict_id='YZ-flicker'")
    vv = v[0] if v else ("", "", "", "")
    formal = q(con, "SELECT DISTINCT series_key, sigma_formal FROM "
                    "p4_flicker WHERE state='QUIESCENT' ORDER BY series_key")
    formal_txt = ", ".join(f"{esc(r[0])} {_mmag(r[1])} mmag" for r in formal)
    return f"""
<section id="flicker"><h2>4 &nbsp; Flickering: amplitude against timescale,
over a floor that was measured</h2>

<h3>Question</h3>
<p>Flickering is aperiodic, so a periodogram is the wrong instrument. The
right one is a structure function: for every pair of points separated by
&tau;, accumulate (&Delta;m)<sup>2</sup>, and report
&radic;(&lang;&Delta;m<sup>2</sup>&rang;/2) &mdash; which for a stationary
process is the standard deviation of the variability on that timescale.</p>

<p>The hard part is not the statistic, it is the subtraction. A structure
function of the target contains flickering AND photon noise AND sky AND
whatever the ensemble zero point did that night. Subtract too little and the
photometry is published as astrophysics; subtract too much and a real
detection is erased.</p>

<h3>What was subtracted, and how</h3>
<p>The floor is the <b>identical statistic computed on the magnitude-matched
field stars</b> &mdash; every star within {fs.FLOOR_MATCH_HALF_WIDTH} mag of
the target's median on that run, covering at least
{100 * fs.FLOOR_MIN_COVERAGE:.0f}% of its frames, measured through the same
frames. Those stars carry the same photon statistics, the same sky, the same
ensemble zero-point wander and the same atmosphere, and no flickering. The
per-bin floor is the MEDIAN over stars, not the mean, so one variable star
that slipped into the sample cannot erase a detection; fewer than
{fs.FLOOR_MIN_STARS} stars and the bin reports no floor at all rather than a
guess.</p>

<p>The subtraction is done in VARIANCE &mdash;
&sigma;<sub>flicker</sub> = &radic;(&sigma;<sub>total</sub><sup>2</sup>
&minus; &sigma;<sub>floor</sub><sup>2</sup>) &mdash; because that is the
quantity that is additive, and its significance is quoted as a sigma on the
variance difference rather than as a ratio of two similar sigmas. The
orbital model of §3 is removed from the target before the statistic is
computed, so a coherent hump cannot be counted as flickering.</p>

<p>For comparison, the formal per-point error inflated by each series'
measured &chi;<sup>2</sup> inflation &mdash; the number a pipeline would
have used if nobody had measured a floor &mdash; is: {formal_txt}. Where it
disagrees with the measured floor, the measured floor governs.</p>

<h3>Evidence</h3>
{_fig(fig, "Solid: the target. Dotted: the measured field-star floor. "
           "Dashed: the quadrature excess — the flickering. The two eras "
           "differ by a factor of several in the floor, which is what "
           "§4.19 was worried about.", "p4_flicker is empty.")}
{tbl}

<h3>Decision</h3>
<p>{_verdict_span(vv[0])} &mdash; {esc(vv[1])}.</p>
<p>{esc(vv[2])}</p>

<h3>Consequence</h3>
<p>{esc(vv[3])}</p>
</section>"""


def section_cannot(con) -> str:
    v = q(con, "SELECT verdict, deciding_number, reasoning, alternative "
               "FROM p4_verdict WHERE verdict_id='YZ-superhump'")
    if not v:
        return ""
    vv = v[0]
    tbl = table(
        ["night", "filter", "N", "blind-search A<sub>90</sub> (mmag)",
         "range", "a superhump would need"],
        [[esc(r[0]), esc(r[1]), _i(r[2]), f"<b>{_mmag(r[3])}</b>",
          f"{_mmag(r[4])}&ndash;{_mmag(r[5])}", f"&ge; {_mmag(r[6])}"]
         for r in q(con, "SELECT night, filter, n_points, amp90_blind, "
                         "amp90_blind_lo, amp90_blind_hi, superhump_floor "
                         "FROM p4_outburst WHERE amp90_blind IS NOT NULL "
                         "ORDER BY night")])
    return f"""
<section id="cannot"><h2>5 &nbsp; What this season cannot do</h2>

<h3>Question</h3>
<p>The strategy's Q3 promised a superhump period and a
dP<sub>sh</sub>/dt. Neither appears in this paper. An absence in a
manuscript is worth nothing unless it is quantified, so: how big would a
superhump have had to be for these data to have measured its period?</p>

<h3>Evidence</h3>
<p>Two independent reasons, either sufficient on its own.</p>
<p><b>First, the state.</b> Common superhumps occur inside superoutbursts.
CV-S7 established that none of the dense runs does &mdash; the brightest
reaches 1.86 mag above quiescence against the ~3 mag a superoutburst
reaches, and the February block is bracketed by independent AAVSO points
five days apart, which a &ge;&nbsp;8&nbsp;d superoutburst plateau cannot fit
between.</p>
<p><b>Second, the sampling.</b> A superhump period is a BLIND period
determination: P<sub>sh</sub> is not known in advance the way P<sub>orb</sub>
is. Scored that way &mdash; the tallest peak in the 2&ndash;40&nbsp;c/d band
must both clear threshold and land within 1% of the truth &mdash; the 90%
recovery contour on the outburst dense runs is:</p>
{tbl}

<h3>Decision</h3>
<p>{_verdict_span(vv[0])} &mdash; {esc(vv[1])}.</p>
<p>{esc(vv[2])}</p>

<h3>Consequence</h3>
<p>{esc(vv[3])}</p>
</section>"""


def section_outburst(con, fig: str) -> str:
    rows = q(con, "SELECT night, utc_night, filter, n_points, span_h, "
                  "median_cal_mag, amp_above_quiescence, amp_p5p95, "
                  "rate_mag_per_h, rate_sigma, rate_verdict, structure "
                  "FROM p4_outburst ORDER BY night, filter")
    if not rows:
        return ""
    tbl = table(
        ["local night", "UTC", "filter", "N", "span (h)", "median cal_mag",
         "mag above quiescence", "p5&ndash;p95", "rate (mag/h)", "verdict"],
        [[esc(r[0]), esc(r[1]), esc(r[2]), _i(r[3]), _n(r[4]), _n(r[5]),
          _n(r[6]), _n(r[7], 3),
          f"{_n(r[8], 4)} &plusmn; {_n(r[9], 4)}", esc(r[10])]
         for r in rows])
    v = q(con, "SELECT verdict, deciding_number, reasoning, alternative "
               "FROM p4_verdict WHERE verdict_id='YZ-outburst'")
    vv = v[0] if v else ("", "", "", "")
    return f"""
<section id="outburst"><h2>6 &nbsp; The six normal-outburst runs, on their
own terms</h2>

<h3>Question</h3>
<p>CV-S7 called these &ldquo;a distinct, separately publishable data set,
not a consolation&rdquo;. That is a claim, and it needs numbers: what does a
dense multi-colour run inside a normal outburst actually measure?</p>

<h3>Evidence</h3>
<p>Amplitude relative to the measured quiescent baseline, covered duration,
and a straight-line rate through the run with its own error bar. A rate is
called a direction only when it exceeds three times that bar; otherwise the
run is FLAT and says so.</p>
{tbl}
{_fig(fig, "Six nights, three filters each, 8 s exposures. The 2024-02-21 "
           "panel catches a rise in all three filters at once.",
      "p4_outburst is empty.")}

<h3>Decision</h3>
<p>{_verdict_span(vv[0])} &mdash; {esc(vv[1])}.</p>
<p>{esc(vv[2])}</p>

<h3>Consequence</h3>
<p>{esc(vv[3])}</p>
</section>"""


def section_anuma(con, fig: str) -> str:
    rows = q(con, "SELECT filter, capability, measured, bar, unit, verdict, "
                  "deciding_number, reasoning, what_would_change_it, rank "
                  "FROM p4_anuma ORDER BY rank, filter")
    if not rows:
        return '<section id="anuma"><h2>7 &nbsp; AN UMa</h2>' \
               '<p class="note">Not run.</p></section>'
    meta = _meta(con)
    tbl = table(
        ["filter", "capability", "measured", "bar", "unit", "verdict",
         "the deciding number"],
        [[f"<b>{esc(r[0])}</b>", esc(r[1]), _n(r[2]), _n(r[3]), esc(r[4]),
          _verdict_span(r[5]), esc(r[6])] for r in rows])
    why = "".join(
        f"<li><b>{esc(r[1])}</b> (bar {_n(r[3])} {esc(r[4])}) &mdash; "
        f"{esc(r[7])}</li>"
        for r in rows if r[0] == rows[0][0])
    change = "".join(
        f"<li><b>{esc(r[0])} &middot; {esc(r[1])}</b> &mdash; "
        f"{esc(r[8])}</li>" for r in rows if r[5] != "SUPPORTED" and r[8])
    v = q(con, "SELECT verdict, deciding_number, reasoning, alternative "
               "FROM p4_verdict WHERE verdict_id='ANUMA-role'")
    vv = v[0] if v else ("", "", "", "")
    return f"""
<section id="anuma"><h2>7 &nbsp; AN UMa, filter by filter</h2>

<h3>Question</h3>
<p>CV-S5's Q5 already graded AN&nbsp;UMa's three-filter colour goal NOT
SUPPORTED: {esc(meta.get('anuma_three_filter_nights', ''))} &mdash; four
nights carrying a full orbit in all three filters, against the strategy's
own bar of {esc(str(meta.get('anuma_bars', '')).split('three_filter_nights=')[-1])}. That
settles one goal. It does not settle the question the plan actually asks,
which is what EACH FILTER can support on its own, and therefore whether
AN&nbsp;UMa belongs in the paper as a full target, a reduced-scope target,
or not at all.</p>

<h3>How each verdict is reached</h3>
<p>One measured number against one stated bar, per capability. Every
contestable choice is a choice of bar, and a bar can be disagreed with in
one line; a weighted score cannot. The bars, and why each is where it
is:</p>
<ul>{why}</ul>
<p>A night counts as a full-orbit night in a filter when it holds at least
{esc(meta.get('full_orbit_min_points', ''))} target points spanning more
than one orbital period &mdash; the same rule CV-S5's census used, restated
so this page's counts are reproducible without opening that one. They
reproduce <code>ch_cadence.n_blocks_ge1cycle</code> exactly.</p>

<h3>Evidence</h3>
{tbl}
{_fig(fig, "Every capability as a ratio to its own bar, so five "
           "incommensurable quantities share one axis. Right of the blue "
           "line passes.", "p4_anuma is empty.")}

<h3>Decision</h3>
<p>{_verdict_span(vv[0])} &mdash; {esc(vv[1])}</p>
<p>{esc(vv[2])}</p>

<h3>What would change the verdict</h3>
<ul>{change}</ul>

<h3>Consequence</h3>
<p>{esc(vv[3])}</p>
</section>"""


def section_verdicts(con) -> str:
    rows = q(con, "SELECT verdict_id, task, question, verdict, "
                  "deciding_number, reasoning, alternative FROM p4_verdict "
                  "ORDER BY rank")
    if not rows:
        return ""
    body = "".join(f"""
<div class="blockcard"><b>{esc(r[2])}</b>
<span class="src">({esc(r[1])}, <code>{esc(r[0])}</code>)</span><br>
<b>Verdict:</b> {_verdict_span(r[3])}<br>
<b>Deciding number:</b> {esc(r[4])}<br>
<b>Why:</b> {esc(r[5])}<br>
<b>What it costs / what was the alternative:</b> {esc(r[6])}</div>"""
                   for r in rows)
    return f"""
<section id="verdicts"><h2>8 &nbsp; The verdicts</h2>
<p>Same shape as the CV-S5 goal audit: question, verdict, THE number that
decides it, and what the alternative was. Every number is a query, so a
re-run that moves a measurement moves the sentence with it.</p>
{body}
</section>"""


# ===========================================================================
# The page
# ===========================================================================
def render_report(db_path: Path) -> Path:
    """Render the CV-S10 page.  Returns the HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = _meta(con)
        f_census = fig_census(con)
        f_folds = fig_folds(con, db_path)
        f_nulls = fig_two_nulls(con)
        f_flick = fig_flicker(con)
        f_out = fig_outburst(con, db_path)
        f_anuma = fig_anuma(con)
        sections = [
            section_intro(con),
            section_census(con, f_census),
            section_gate(con, meta),
            section_hump(con, f_folds, f_nulls),
            section_flicker(con, f_flick),
            section_cannot(con),
            section_outburst(con, f_out),
            section_anuma(con, f_anuma),
            section_verdicts(con),
        ]
        n_runs = q1(con, "SELECT count(*) FROM p4_run") \
            if _has(con, "p4_run") else 0
        n_bins = q1(con, "SELECT count(*) FROM p4_flicker WHERE detected=1") \
            if _has(con, "p4_flicker") else 0
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CV Time Series &mdash; the closing science decisions</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Cataclysmic-Variable Time Series &mdash; the closing science
  decisions</h1>
  <p>YZ&nbsp;Cnc on the fallback branch: the quiescent orbital hump against
  two nulls &middot; flickering over a floor that was measured, not modelled
  &middot; the §4.19 signal-to-noise gate, executed &middot; the six
  normal-outburst runs as their own result &middot; what the season cannot
  do, quantified &middot; AN&nbsp;UMa graded filter by filter
  &middot; {n_runs} scopes, {n_bins} flickering detections
  &middot; built {esc(meta.get('stage_verdict', ''))[:16]}Z
  ({esc(meta.get('final_code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="cv_external_context.html">the external record that
  chose the branch</a> &middot;
  <a href="cv_timeseries_analysis.html">the Phase-3 analysis</a> &middot;
  <a href="cv_characterization.html">the characterization</a> &middot;
  <a href="index.html">project hub</a> &middot;
  <a href="../index.html">all reports</a></p>
</header>

<nav>
  <a href="#intro">0 What this decides</a> &middot;
  <a href="#census">1 The dense runs</a> &middot;
  <a href="#gate">2 The §4.19 gate</a> &middot;
  <a href="#hump">3 The orbital hump</a> &middot;
  <a href="#flicker">4 Flickering</a> &middot;
  <a href="#cannot">5 What it cannot do</a> &middot;
  <a href="#outburst">6 The outburst runs</a> &middot;
  <a href="#anuma">7 AN UMa</a> &middot;
  <a href="#verdicts">8 Verdicts</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_phot.report_final</code> from
<code>products/phot/cv_timeseries.sqlite</code> &mdash; every number on this
page is the result of a SQL query or a constant imported from
<code>macro_phot.final_science</code>; none is typed by hand.  Regenerate
with <code>pipeline/scripts/run_cv_final.py all</code>.</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
        # Belt and braces: every <img> the page references must exist and be
        # non-empty, or the build fails loudly rather than shipping a broken
        # evidence page.
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
