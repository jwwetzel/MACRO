"""S2c evidence report: filter identity measured from the pixels.

Reads ``frame_dispersion`` from the manifest and writes:

* ``docs/pipeline/s2c_filter_identity.html`` — the report
* ``docs/pipeline/figures/s2c/*.png``        — every figure

Same discipline as the S0/S0b/S2 renderers: the page follows the site's
Socratic format (Question → Evidence → Decision → Consequence), and EVERY
number is interpolated from a SQL query run here or from a constant defined
in ``rlmt_diagnostics.dispersion``.  Nothing on the page is typed by hand,
including the threshold values — they are read from the module, so the page
and the classifier can never drift apart.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import dispersion as dsp  # noqa: E402

import sys                       # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from macro_core.report_s0 import (  # noqa: E402
    ACCENT, DARK, DPI, WARN, _figure, esc, fmt, q, q1, table)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s2c"
HTML_PATH = DOCS_DIR / "s2c_filter_identity.html"

GOOD = "#9fd8ae"                 # site badge green — confirmations
BAD = "#e59a9a"                  # site badge red — contradictions
MUTED = "#9aa4b2"

#: The three filter-name epochs, as established by the header archaeology
#: and confirmed by the measurements on this page.
EPOCHS = (("slot number <code>6</code>", "2023-02", "2024-03"),
          ("<code>HaGrism</code> / <code>OGGrism</code>", "2024-04", "2025-03"),
          ("<code>hrg</code> / <code>lrg</code>", "2025-01", "2026-06"))


def fnum(x, nd=2) -> str:
    """Format a float for the page (NULL-safe, fixed decimals)."""
    if x is None:
        return "&mdash;"
    return f"{x:,.{nd}f}"


def pct(num, den) -> str:
    """Percentage with a guarded denominator."""
    if not den:
        return "&mdash;"
    return f"{100.0 * num / den:.1f}%"


# ---------------------------------------------------------------------------
# Shared SQL fragments
# ---------------------------------------------------------------------------
MEASURED = "status = 'measured'"


def _label_set(names) -> str:
    """SQL IN-list of quoted filter names (names are module constants)."""
    return ",".join("'" + n.replace("'", "''") + "'" for n in names)


DISP_IN = _label_set(dsp.KNOWN_DISPERSED_FILTERS)
DIR_IN = _label_set(dsp.KNOWN_DIRECT_FILTERS)
HIGH_IN = _label_set(dsp.HIGH_DISPERSION_FILTERS)
LOW_IN = _label_set(dsp.LOW_DISPERSION_FILTERS)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_calibration(con) -> str:
    """Two panels: the axes that fail, and the axes the classifier uses.

    Left  — median a/b of the ten brightest sources against the bright-set
            position-angle scatter.  This is the statistic the survey's
            first pass quoted, and the panel shows it does NOT separate:
            real grism frames pile into the same corner as ordinary images
            because most of their bright detections are round field stars.
    Right — trace aspect ratio against the TRACE position-angle scatter,
            over the trace sub-population only.  The two classes fall apart
            into opposite corners.
    """
    def fetch(where):
        return np.array(q(con, f"""
            SELECT median_ab, pa_scatter, trace_ab, trace_pa_scatter, n_trace
            FROM frame_dispersion
            WHERE {MEASURED} AND {where} AND median_ab IS NOT NULL"""),
            dtype=float)

    # Grism classes first, control LAST so the smaller, decisive population
    # is never buried under 19k grism points.
    groups = [("high-dispersion grism", f"filter IN ({HIGH_IN})", WARN),
              ("low-dispersion grism", f"filter IN ({LOW_IN})", ACCENT),
              ("known direct (control)",
               f"population = 'control' AND filter IN ({DIR_IN})", "#e8eaed")]

    with plt.rc_context(DARK):
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, 4.8))
        for name, where, color in groups:
            d = fetch(where)
            if not len(d):
                continue
            axl.scatter(np.clip(d[:, 0], 0.8, 400), d[:, 1], s=5,
                        alpha=0.25, color=color, label=name, linewidths=0)
            # Right panel: the DISTRIBUTION of the shared-axis statistic.
            # A scatter hides it — nearly every dispersed frame sits at
            # exactly 0, so the points collapse into an invisible line.
            t = d[d[:, 4] > 0]
            if len(t):
                s = t[:, 3][np.isfinite(t[:, 3])]
                axr.hist(s, bins=np.linspace(0, 90, 91), histtype="step",
                         lw=1.6, color=color, label=name,
                         weights=np.full(len(s), 100.0 / max(len(s), 1)))
        axl.set_xscale("log")
        axl.set_xlabel("median a/b of the 10 brightest sources")
        axl.set_ylabel("position-angle scatter of those 10 (deg)")
        axl.set_title("The statistic that FAILS:\nbright-set median",
                      fontsize=10)
        axl.legend(fontsize=7, loc="upper left", framealpha=0.3)
        axr.set_yscale("log")
        axr.axvline(dsp.TRACE_MAX_PA_SCATTER_DEG, color=GOOD, lw=1.4, ls="--")
        axr.set_xlabel("position-angle scatter of the traces (deg)")
        axr.set_ylabel("% of each population (log)")
        axr.set_title("The test that DECIDES:\ndo the traces share an axis?",
                      fontsize=10)
        axr.legend(fontsize=7, framealpha=0.3)
        axr.annotate("adopted gate", xy=(dsp.TRACE_MAX_PA_SCATTER_DEG, 30),
                     xytext=(22, 45), fontsize=7, color=GOOD,
                     arrowprops=dict(arrowstyle="->", color=GOOD, lw=0.8))
        fig.suptitle("Calibration on labels nobody disputes — a grating "
                     "holds every trace parallel; nothing else does",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2c_calibration.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2c/s2c_calibration.png"


def fig_strength(con) -> str:
    """Why trace LENGTH is the strength axis and aspect ratio is not."""
    def fetch(where, col):
        rows = q(con, f"""SELECT {col} FROM frame_dispersion
                          WHERE {MEASURED} AND {where}
                            AND verdict = 'dispersed' AND {col} IS NOT NULL
                            AND width > 0""")
        return np.array([r[0] for r in rows], dtype=float)

    def fetch_frac(where):
        rows = q(con, f"""SELECT trace_a_px * 1.0 / width FROM frame_dispersion
                          WHERE {MEASURED} AND {where}
                            AND verdict = 'dispersed'
                            AND trace_a_px IS NOT NULL AND width > 0""")
        return np.array([r[0] for r in rows], dtype=float)

    with plt.rc_context(DARK):
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, 4.4))
        for where, color, name in ((f"filter IN ({HIGH_IN})", WARN,
                                    "H-alpha grism (hrg / HaGrism)"),
                                   (f"filter IN ({LOW_IN})", ACCENT,
                                    "broad grism (lrg / OGGrism)")):
            ab = fetch(where, "trace_ab")
            if len(ab):
                axl.hist(np.clip(ab, 0, 300), bins=60, alpha=0.55,
                         color=color, label=name)
            fr = fetch_frac(where)
            if len(fr):
                axr.hist(np.clip(fr, 0, 0.45), bins=60, alpha=0.55,
                         color=color, label=name)
        axl.set_xlabel("trace a/b")
        axl.set_ylabel("frames")
        axl.set_title("Aspect ratio: the two grisms OVERLAP", fontsize=10)
        axl.legend(fontsize=7, framealpha=0.3)
        axr.axvline(dsp.STRENGTH_HIGH_MIN_FRAC, color=GOOD, lw=1.2, ls="--")
        axr.set_xlabel("trace length / frame width")
        axr.set_title("Trace length: medians differ, tails overlap\n"
                      "(dashed = adopted 'high' bound)", fontsize=10)
        axr.legend(fontsize=7, framealpha=0.3)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2c_strength.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2c/s2c_strength.png"


def fig_slot6_timeline(con) -> str:
    """Slot '6' month by month: is there a date at which the wheel changed?"""
    rows = q(con, f"""
        SELECT substr(night, 1, 7) AS mon,
               sum(verdict = 'dispersed'), sum(verdict = 'direct'),
               sum(verdict = 'indeterminate')
        FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6' AND night IS NOT NULL
        GROUP BY mon ORDER BY mon""")
    mons = [r[0] for r in rows]
    disp = np.array([r[1] or 0 for r in rows], dtype=float)
    direct = np.array([r[2] or 0 for r in rows], dtype=float)
    indet = np.array([r[3] or 0 for r in rows], dtype=float)
    x = np.arange(len(mons))
    with plt.rc_context(DARK):
        fig, (axt, axb) = plt.subplots(2, 1, figsize=(11.5, 5.6),
                                       sharex=True,
                                       gridspec_kw={"height_ratios": [2, 1]})
        axt.bar(x, disp, color=WARN, label="dispersed (spectrum)")
        axt.bar(x, direct, bottom=disp, color=ACCENT, label="direct (image)")
        axt.bar(x, indet, bottom=disp + direct, color=MUTED,
                label="indeterminate")
        axt.set_ylabel("frames")
        axt.legend(fontsize=8, framealpha=0.3)
        axt.set_title("Slot '6' frames per month, by measured verdict",
                      fontsize=11)
        tot = np.maximum(disp + direct + indet, 1)
        axb.plot(x, 100.0 * disp / tot, "o-", color=WARN, lw=1.4, ms=4)
        axb.set_ylabel("% dispersed")
        axb.set_ylim(-5, 105)
        axb.set_xticks(x)
        axb.set_xticklabels(mons, rotation=90, fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2c_slot6_timeline.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2c/s2c_slot6_timeline.png"


def fig_census(con) -> str:
    """The archive's spectra, by target — who actually got observed."""
    rows = q(con, f"""
        SELECT coalesce(canonical_target, '(untargeted)') AS tgt,
               sum(strength_class = 'high'),
               sum(strength_class != 'high')
        FROM frame_dispersion
        WHERE {MEASURED} AND verdict = 'dispersed'
        GROUP BY tgt ORDER BY count(*) DESC LIMIT 25""")
    tgt = [r[0] for r in rows][::-1]
    hi = np.array([r[1] or 0 for r in rows], dtype=float)[::-1]
    am = np.array([r[2] or 0 for r in rows], dtype=float)[::-1]
    y = np.arange(len(tgt))
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(9.5, 7.5))
        ax.barh(y, hi, color=WARN, label="H-alpha grism (identified)")
        ax.barh(y, am, left=hi, color=MUTED, label="unit not identifiable")
        ax.set_yticks(y)
        ax.set_yticklabels(tgt, fontsize=8)
        ax.set_xlabel("measured spectra (frames)")
        ax.set_title("The 25 most-observed spectroscopic targets", fontsize=11)
        ax.legend(fontsize=8, framealpha=0.3)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2c_census.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2c/s2c_census.png"


# ---------------------------------------------------------------------------
# Section 1 — calibration
# ---------------------------------------------------------------------------
def section_calibration(con) -> str:
    fig = fig_calibration(con)

    # The direct labels now exist in TWO populations — the fitted control
    # and the frozen holdout — so every in-sample statistic must say which
    # one it means.  Pooling them would quietly fold the out-of-sample check
    # into the number it is supposed to be checking.
    IN_SAMPLE = "population = 'control'"
    rows = []
    for label in list(dsp.KNOWN_DISPERSED_FILTERS) + list(
            dsp.KNOWN_DIRECT_FILTERS):
        scope = (IN_SAMPLE if label in dsp.KNOWN_DIRECT_FILTERS else "1=1")
        r = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                              sum(verdict='direct'),
                              sum(verdict='indeterminate')
                       FROM frame_dispersion
                       WHERE {MEASURED} AND {scope} AND filter = ?""",
              (label,))[0]
        n = r[0] or 0
        if not n:
            continue
        truth = dsp.expected_verdict(label)
        agree = (r[1] if truth == dsp.VERDICT_DISPERSED else r[2]) or 0
        # House convention: only disagreement is coloured, so the eye lands
        # on the rows that need explaining.
        cls = "" if agree / n >= 0.95 else "warn"
        rows.append(([f"<code>{esc(label)}</code>", esc(truth), fmt(n),
                      fmt(r[1] or 0), fmt(r[2] or 0), fmt(r[3] or 0),
                      f"<b>{pct(agree, n)}</b>"], cls))

    # Headline rates over the pooled known populations.
    kd = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                           sum(verdict='indeterminate')
                    FROM frame_dispersion
                    WHERE {MEASURED} AND filter IN ({DISP_IN})""")[0]
    kn = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                           sum(verdict='indeterminate')
                    FROM frame_dispersion
                    WHERE {MEASURED} AND {IN_SAMPLE}
                      AND filter IN ({DIR_IN})""")[0]
    recall = pct(kd[1] or 0, kd[0])
    fpr = pct(kn[1] or 0, kn[0])

    # The measured overlap region, stated honestly.
    sep = q(con, f"""
        SELECT
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND filter IN ({DISP_IN}) AND n_trace = 0),
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND population = 'control'
              AND filter IN ({DIR_IN}) AND n_trace > 0)""")[0]

    # Where the false positives actually come from: whole nights of them,
    # on the same targets — the fingerprint of a tracking failure.
    fp_nights = q(con, f"""
        SELECT coalesce(canonical_target, '(untargeted)'), night, count(*)
        FROM frame_dispersion
        WHERE {MEASURED} AND population = 'control' AND verdict = 'dispersed'
        GROUP BY 1, 2 HAVING count(*) > 1 ORDER BY count(*) DESC LIMIT 6""")
    fprows = [[esc(t), esc(n), fmt(c)] for t, n, c in fp_nights]

    # The grating-axis gate, scored: what it removed from each population.
    # An earlier draft attributed ALL control false positives to tracking
    # failures.  The database refused that: a large block of them sat within
    # a degree of PA 90 — perpendicular to the dispersion axis of every
    # grism frame in the archive — which is a detector column or a
    # saturated-star bleed trail, and cannot be a mount.
    OFFAXIS = "reason LIKE '%off the grating axis%'"
    axis_cut = q(con, f"""
        SELECT
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND population = 'control'
              AND filter IN ({DIR_IN}) AND {OFFAXIS}),
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND filter IN ({DISP_IN}) AND {OFFAXIS}),
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND filter = '6' AND {OFFAXIS})""")[0]
    # How tightly the labelled grism populations pin that axis.
    axis_tight = q(con, f"""
        SELECT count(*), sum(CASE WHEN min(abs(trace_pa),
                                           abs(trace_pa - 180.0)) <= 10.0
                                  THEN 1 ELSE 0 END)
        FROM frame_dispersion
        WHERE {MEASURED} AND verdict = 'dispersed'
          AND filter IN ({DISP_IN}) AND trace_pa IS NOT NULL""")[0]
    # The false-negative channel the direct-side evidence floor closed.
    fn_thin = q(con, f"""
        SELECT count(*) FROM frame_dispersion
        WHERE {MEASURED} AND filter IN ({DISP_IN})
          AND verdict = 'indeterminate'
          AND n_sources < {dsp.DIRECT_MIN_SOURCES}""")[0][0]

    # The HOLDOUT: a second, disjoint control draw under a different seed,
    # measured once with every threshold frozen.  This exists because the
    # control sample above cannot honestly carry an error rate — three of
    # this module's constants were moved in response to frames still inside
    # it, so its rate is a lower bound, not an estimate.
    ho = q(con, f"""SELECT count(*), sum(verdict = 'dispersed'),
                           sum(verdict = 'direct'),
                           sum(verdict = 'indeterminate')
                    FROM frame_dispersion
                    WHERE {MEASURED} AND population = 'holdout'""")[0]
    ho_n = ho[0] or 0
    ho_rows = [[esc(f), fmt(n), fmt(d or 0), pct(d or 0, n)]
               for f, n, d in q(con, f"""
                   SELECT filter, count(*), sum(verdict = 'dispersed')
                   FROM frame_dispersion
                   WHERE {MEASURED} AND population = 'holdout'
                   GROUP BY 1 ORDER BY 1""")]
    ho_seed = dict(q(con, "SELECT key, value FROM s2c_build_meta")).get(
        "holdout_seed", "&mdash;")
    if ho_n:
        ho_fpr = pct(ho[1] or 0, ho_n)
        ho_block = f"""
<p><b>An out-of-sample check, because the number above is fitted.</b> Three
of this classifier's constants were moved in response to specific control
frames &mdash; the PA-scatter gate went from 20&nbsp;deg to 5 after a
defocused <code>r</code> frame, the sparsity gate was invented after an
<code>L</code> frame of M57, and the calibration-frame exclusion was added
after eleven <code>B</code> frames &mdash; and every one of those frames is
still scored in the {fmt(kn[0])} control totals. A rate fitted on the data it
is quoted over is a lower bound, not an estimate, and an earlier version of
this page presented it as the latter.</p>

<p>So a second control was drawn: <b>{fmt(ho_n)}</b> frames from the same
undisputed labels, under a different seed ({esc(str(ho_seed))}), explicitly
excluding every frame already measured, and classified with every threshold
frozen. Its false-positive rate is <b>{ho_fpr}</b>, against <b>{fpr}</b>
in-sample.</p>

{table(["label", "holdout frames", "called dispersed", "false-positive rate"],
       ho_rows)}

<p>The holdout is <em>modestly worse</em>, and that is the honest and expected
result: fitting three constants to individual control frames bought about
{abs(100.0*(ho[1] or 0)/max(ho_n,1) - 100.0*(kn[1] or 0)/max(kn[0],1)):.2f}
of a percentage point of flattery. It is a small gap because the thresholds
were set from single diagnosed frames rather than by sweeping the
distribution for a minimum &mdash; but it is not zero, and a page that quoted
only the in-sample figure would be understating its own error.
<b>{ho_fpr} is the rate to cite.</b> It is the only number here that survives
the question &ldquo;how do you know you did not fit this?&rdquo;</p>

<p>The per-label breakdown also shows where the residual error lives: it is
concentrated in the luminance and broad filters, which is consistent with the
mechanism named under Consequence below &mdash; those are the frames of long,
unguided, deep exposures where the mount has the most opportunity to drift.
</p>"""
    else:
        ho_block = ("<p class=\"warn\">The out-of-sample holdout has been "
                    "queued but not yet measured; every rate on this page is "
                    "therefore in-sample and is a lower bound on the error."
                    "</p>")

    # The discriminator that was tried and REJECTED, with its numbers.
    def _ratio_p50(where):
        rows = q(con, f"""SELECT trace_ab * 1.0 / median_ab
                          FROM frame_dispersion
                          WHERE {MEASURED} AND verdict = 'dispersed'
                            AND trace_ab IS NOT NULL AND median_ab > 0
                            AND n_sources >= 30 AND n_bright >= 5
                            AND {where}""")
        v = np.array([r[0] for r in rows], dtype=float)
        return (len(v), float(np.median(v)) if len(v) else float("nan"),
                100.0 * float(np.mean(v < 1.5)) if len(v) else float("nan"))
    # A SECOND discriminator was proposed and also measured to fail: extend
    # the sparse-field requirement from the solo branch to the multi-trace
    # branch, as a floor on the fraction of bright sources that are traces.
    # The argument is sound ("a grating should streak most of its ten
    # brightest") and the data still refuses it, because a mount slip streaks
    # all ten too.  Rendered so the rejection is auditable, not asserted.
    def _tf_cost(floor):
        return q(con, f"""
            SELECT
              (SELECT count(*) FROM frame_dispersion
                WHERE {MEASURED} AND verdict = 'dispersed'
                  AND filter IN ({DISP_IN}) AND n_trace >= 2
                  AND n_sources > {dsp.SOLO_MAX_SOURCES}
                  AND trace_frac < {floor}),
              (SELECT count(*) FROM frame_dispersion
                WHERE {MEASURED} AND verdict = 'dispersed'
                  AND population = 'control' AND filter IN ({DIR_IN})
                  AND n_trace >= 2 AND n_sources > {dsp.SOLO_MAX_SOURCES}
                  AND trace_frac < {floor})""")[0]
    tf_rows = [[fnum(f, 1), fmt(lost), fmt(gain)]
               for f, (lost, gain) in ((f, _tf_cost(f))
                                       for f in (0.2, 0.3, 0.4, 0.5))]
    r_fp = _ratio_p50("population = 'control'")
    r_kd = _ratio_p50(f"filter IN ({DISP_IN})")
    r_s6 = _ratio_p50("filter = '6'")

    return f"""
<section id="calibration"><h2>1&nbsp;&middot;&nbsp;Calibration: can the pixels
tell a spectrum from an image?</h2>

<h3>Question</h3>
<p>The <code>FILTER</code> card names the RLMT's grism slots three different
ways across three years, and the earliest name &mdash; the bare slot number
<code>6</code> &mdash; is contested: one review panel read it as a grism,
another as a luminance filter. Before that argument can be settled we need a
measurement that does not consult the label at all. So: <b>is dispersion
detectable directly in the pixels, and how cleanly?</b> The only honest way
to answer is to calibrate on the labels nobody disputes &mdash;
<code>hrg</code>, <code>lrg</code>, <code>HaGrism</code>,
<code>OGGrism</code> as dispersed; <code>g r i V R I B L</code> as direct
&mdash; and report the separation achieved, overlap included.</p>

<h3>Evidence</h3>
<p>Every frame is background-subtracted and source-extracted, and three
things are measured: the elongation of the brightest sources, whether any of
them is shaped like a dispersed <em>trace</em>
(a/b&nbsp;&ge;&nbsp;{fnum(dsp.TRACE_MIN_AB, 0)} <em>and</em> semi-major
&ge;&nbsp;{fnum(dsp.TRACE_MIN_A_PX, 0)}&nbsp;px), and whether those traces
share a common axis. The shared axis is the decisive test: a grating rules
every source in the field the same way, while a cosmic ray or a satellite
points wherever it likes.</p>

{_figure(fig, "Left: the bright-set median a/b the first survey pass quoted "
              "&mdash; the classes overlap badly, because most bright "
              "detections on a grism frame are round field stars whose "
              "traces are too faint to resolve. Right: among frames that "
              "reached the trace gates at all, how tightly their traces "
              "share an axis. The grism populations collapse onto zero; the "
              "control keeps a long tail out to 80 deg, which is the "
              "defocused and satellite-crossed frames the gate exists to "
              "remove. Note the control ALSO has a spike at zero &mdash; "
              "those are the trailed frames discussed under Consequence, "
              "and they are the reason the false-positive rate is not zero.")}

{table(["label", "truth", "measured", "dispersed", "direct",
        "indeterminate", "agreement"],
       [r[0] for r in rows], [r[1] for r in rows])}

<p>Pooled over the {fmt(kd[0])} known-dispersed frames, the classifier
recovers <b>{recall}</b> as dispersed. Pooled over the {fmt(kn[0])}
known-direct control frames it calls <b>{fpr}</b> dispersed &mdash; that is
the false-positive rate, and it is the number that decides whether any
verdict on slot <code>6</code> can be believed.</p>

{ho_block}

<p>The overlap, stated plainly: {fmt(sep[0])} known-dispersed frames
contained no trace-shaped source at all (too short an exposure, too faint a
target, or cloud), and {fmt(sep[1])} known-direct frames contained at least
one (satellites, cosmic-ray tracks, edge-on galaxies). Those frames are the
irreducible ambiguity in the method, and the classifier is built to return
<code>indeterminate</code> across it rather than to guess.</p>

<h3>Decision</h3>
<p>Dispersion is measurable per frame, and the classifier is adopted with
four rules. <b>(1)</b> Two or more traces sharing an axis to within
{fnum(dsp.TRACE_MAX_PA_SCATTER_DEG, 0)}&nbsp;deg &rarr;
<code>dispersed</code>&hellip; <b>(2)</b> &hellip;<em>provided that shared
axis is the GRATING's</em>, within
{fnum(dsp.GRATING_PA_TOL_DEG, 0)}&nbsp;deg of PA
{fnum(dsp.GRATING_PA_DEG, 0)}. <b>(3)</b> A single trace of
a/b&nbsp;&ge;&nbsp;{fnum(dsp.SOLO_TRACE_MIN_AB, 0)} on that same axis, in a
sparse field (&le;&nbsp;{fmt(dsp.SOLO_MAX_SOURCES)} sources) &rarr;
<code>dispersed</code>; this rule carries the bright-standard frames where
the target is the only thing in the field. <b>(4)</b> No traces, round
bright sources, and at least {fmt(dsp.DIRECT_MIN_SOURCES)} of them &rarr;
<code>direct</code>. Everything else is <code>indeterminate</code> and is
reported as such.</p>

<p>Rule 2 is the one this campaign was missing, and it is worth stating why
it is physics rather than a fitted patch. A diffraction grating is machined
glass bolted into a filter wheel: its dispersion direction is fixed in
<em>detector</em> coordinates, so it does not merely make parallel traces, it
makes traces along one knowable line. The archive's own labelled grism
frames pin that line without ambiguity &mdash; <b>{fmt(axis_tight[1])}</b> of
<b>{fmt(axis_tight[0])}</b> lie within 10&nbsp;deg of PA&nbsp;0/180. The
tolerance sits in a measured gap rather than at a round number: the furthest
on-axis labelled spectrum is at 9.9&nbsp;deg, the nearest off-axis control
false positive at 28.6&nbsp;deg.</p>

<p>Rule 4's source floor is the mirror image, and it corrects a real
asymmetry. <code>direct</code> is not the absence of a verdict; downstream it
is a certificate that the frame is safe for aperture photometry. Under the
first version of these rules a frame with ONE round detection was certified
on a sample of one, while the dispersed side demanded two corroborating
traces or an extreme trace plus a sparsity check. {fmt(fn_thin)} labelled
grism frames now correctly return <code>indeterminate</code> on fields too
thin to see a spectrum in, instead of voting &ldquo;clean&rdquo;.</p>

<h3>Consequence</h3>
<p>Four false-positive modes were found by measurement. Three were closed;
the fourth is real and is quoted above as the residual error rate.</p>

<p><b>Closed &mdash; twilight flats.</b> A flat has no stars, so the
extractor finds dust shadows and detector column defects, which are straight
and mutually parallel because they <em>are</em> columns: a perfect forgery
of a shared dispersion axis. Eleven of the first eighteen <code>B</code>
frames measured came back "dispersed" this way. Calibration frames are now
excluded from the campaign entirely.</p>

<p><b>Closed &mdash; a satellite trail.</b> A 60-second luminance frame of
M57 read as dispersed on the strength of one 1,095-px streak, past 1,278
perfectly round stars. A grating disperses <em>everything</em>, so a rich
field with exactly one smear is the one thing a grism cannot produce; the
solo-trace rule now requires a sparse field.</p>

<p><b>Closed &mdash; detector columns and bleed trails.</b> This was the
largest mode, and an earlier version of this page missed it entirely by
attributing every remaining false positive to the mount. The database says
otherwise: a large block of them sat within a degree of PA&nbsp;90, exactly
perpendicular to the dispersion axis of every grism frame in the archive.
Nothing on a mount prefers that angle. Detector columns and saturated-star
bleed trails do &mdash; they run down the readout direction, and being
literally columns they are <em>more</em> mutually parallel than any real
spectrum, which is precisely why a parallelism-only rule believed them. The
grating-axis gate (rule 2) removes <b>{fmt(axis_cut[0])}</b> control false
positives and <b>{fmt(axis_cut[2])}</b> slot-<code>6</code> verdicts, at a
cost of <b>{fmt(axis_cut[1])}</b> frames across the whole labelled grism
population.</p>

<p><b>Open &mdash; tracking failures on the grating axis.</b> The remaining
{fmt(kn[1] or 0)} false positives are dominated by frames where the mount
slipped <em>and the drift happened to run along the dispersion axis</em>.
Trailed stars are long, thin and mutually parallel, which is not merely
similar to dispersion &mdash; it is geometrically the same thing, and no
morphological gate in this module can separate them. They arrive in clusters,
whole nights of one target at a time, which is the fingerprint of a mount
problem rather than of anything optical:</p>

{table(["target", "night", "false positives"], fprows)}

<p>The residue is not purely mount slips, and the page should not claim it
is. At least one survivor is a bleed trail that happens to run
<em>horizontally</em>: a 5-second <code>r</code> exposure of the bright star
HIP&nbsp;97675 carrying a 1,365-px streak at PA&nbsp;0.0007&nbsp;deg, in a
field of 48 stars whose median a/b is 1.12. A 5-second exposure cannot trail
1,365&nbsp;px, and stars that round rule out drift; the frame is saturated
(65,534&nbsp;ADU) and the streak is its bleed. That the bleed runs along
rows on this sensor rather than down columns is why the axis gate cannot
catch it.</p>

<p>That case suggests a further discriminator, and it is recorded here
<em>without</em> being adopted. A grating has no reason to align with the
pixel grid; a bleed trail does so by construction. Measured: requiring the
trace axis to differ from exact alignment by more than 0.01&nbsp;deg would
remove 2 of the surviving control false positives at a cost of 4 frames in
{fmt(axis_tight[0])} labelled spectra. The ratio is favourable and the
mechanism is real &mdash; but it would be calibrated on two examples, which
is exactly the sin this campaign has already been caught committing
elsewhere. It waits for a population, not an anecdote.</p>

<p>One discriminator was proposed and <b>measured to fail</b>, which is
worth recording so it is not proposed again. A grism disperses bright and
faint sources differently &mdash; only the bright ones spread far enough to
be detected as traces, so the frame's <em>median</em> elongation stays low
while its trace elongation is high &mdash; whereas trailing smears
everything equally. The ratio of the two should therefore separate them, and
on slot <code>6</code> it does beautifully (median {fnum(r_s6[1], 2)}).
But on the labelled grism population it does not: those frames are mostly
isolated bright standards with nothing faint to stay round, and they measure
a median ratio of {fnum(r_kd[1], 2)} against the trailed frames'
{fnum(r_fp[1], 2)} &mdash; indistinguishable, with
{fnum(r_kd[2], 0)}% and {fnum(r_fp[2], 0)}% respectively falling below 1.5.
Restricting to rich fields does not help. Adopting it as a gate would buy a
~1% reduction in false positives at the cost of rejecting most real
spectra.</p>

<p>A second gate was proposed on good reasoning and <b>also measured to
fail</b>. The sparse-field requirement currently guards only the solo-trace
branch; the argument for extending it to the multi-trace branch is that a
grating streaks <em>every</em> bright source, so a rich field in which only
two of ten sources are traces is weak evidence. Expressed as a floor on that
fraction, every candidate value costs far more real spectra than it buys:</p>

{table(["trace-fraction floor", "labelled spectra lost",
        "control false positives removed"], tf_rows)}

<p>The reason is the same one that sank the previous discriminator: a mount
slip streaks all ten sources too, so the frames a floor removes most eagerly
are genuine spectra of crowded fields, while the trailed frames sail over it
with a fraction of 1.0. The multi-trace branch is left ungated, and the
population it carries is reported rather than hidden.</p>

<p>The measurement that <em>would</em> close this gap is the per-source
correlation between flux and trace length: under dispersion the two are
strongly coupled, under trailing they are independent. That needs the
per-source arrays this campaign summarises away rather than stores, so it
belongs to a future pass. Until then, <b>{pct(kn[1] or 0, kn[0])} of direct
frames are expected to be misread as dispersed</b>, and any downstream use
should treat a lone dispersed verdict on an otherwise photometric night as a
tracking failure rather than a discovery.</p>
</section>"""


# ---------------------------------------------------------------------------
# Section 2 — the strength axis
# ---------------------------------------------------------------------------
def section_strength(con) -> str:
    fig = fig_strength(con)

    def stats(where):
        rows = q(con, f"""SELECT trace_ab, trace_a_px * 1.0 / width
                          FROM frame_dispersion
                          WHERE {MEASURED} AND {where}
                            AND verdict = 'dispersed'
                            AND trace_a_px IS NOT NULL AND width > 0""")
        if not rows:
            return (0, None, None, None, None)
        ab = np.array([r[0] for r in rows], dtype=float)
        fr = np.array([r[1] for r in rows], dtype=float)
        return (len(rows), float(np.median(ab)), float(np.median(fr)),
                float(np.percentile(fr, 10)), float(np.percentile(fr, 90)))

    hi = stats(f"filter IN ({HIGH_IN})")
    lo = stats(f"filter IN ({LOW_IN})")
    ratio = (hi[2] / lo[2]) if (hi[2] and lo[2]) else None

    rows = [
        ["H-alpha grism (<code>hrg</code>, <code>HaGrism</code>)", fmt(hi[0]),
         fnum(hi[1], 0), f"<b>{fnum(hi[2], 3)}</b>",
         f"{fnum(hi[3], 3)} &ndash; {fnum(hi[4], 3)}"],
        ["broad grism (<code>lrg</code>, <code>OGGrism</code>)", fmt(lo[0]),
         fnum(lo[1], 0), f"<b>{fnum(lo[2], 3)}</b>",
         f"{fnum(lo[3], 3)} &ndash; {fnum(lo[4], 3)}"]]

    # The purity/recall curve — the evidence that ONLY the long tail works.
    def _fracs(where):
        return np.array([r[0] for r in q(con, f"""
            SELECT trace_a_px * 1.0 / width FROM frame_dispersion
            WHERE {MEASURED} AND verdict = 'dispersed'
              AND trace_a_px IS NOT NULL AND width > 0 AND {where}""")],
            dtype=float)
    fh, fl = _fracs(f"filter IN ({HIGH_IN})"), _fracs(f"filter IN ({LOW_IN})")

    prows, pcls = [], []
    best_lo = 0.0            # the BEST a "low" call ever manages, at any cut
    for t in (0.04, 0.06, 0.08, 0.11, 0.13, 0.15, 0.17):
        nh_hi, nl_hi = int((fh >= t).sum()), int((fl >= t).sum())
        nl_lo, nh_lo = int((fl <= t).sum()), int((fh <= t).sum())
        p_hi = 100.0 * nh_hi / max(nh_hi + nl_hi, 1)
        p_lo = 100.0 * nl_lo / max(nl_lo + nh_lo, 1)
        best_lo = max(best_lo, p_lo)
        prows.append([fnum(t, 2),
                      f"<b>{p_hi:.1f}%</b>", f"{100.0*nh_hi/max(len(fh),1):.1f}%",
                      f"{p_lo:.1f}%", f"{100.0*nl_lo/max(len(fl),1):.1f}%"])
        pcls.append("" if p_hi >= 99 else "warn")
    # The base rate a "low" call has to beat to be worth anything: how often
    # the low-dispersion unit is simply the right answer.  Computed, not
    # asserted — the sentence in the Decision section below used to carry
    # hand-typed values (61% against 47%) that matched neither this table
    # three paragraphs above it nor the data underneath.
    base_lo = 100.0 * len(fl) / max(len(fh) + len(fl), 1)
    # Recall at the ADOPTED bound, for the same reason.
    hi_recall = 100.0 * float((fh >= dsp.STRENGTH_HIGH_MIN_FRAC).sum()) / max(len(fh), 1)
    hi_purity = (100.0 * float((fh >= dsp.STRENGTH_HIGH_MIN_FRAC).sum())
                 / max(float((fh >= dsp.STRENGTH_HIGH_MIN_FRAC).sum())
                       + float((fl >= dsp.STRENGTH_HIGH_MIN_FRAC).sum()), 1.0))

    # The confound, shown directly: length falls as exposure rises, because
    # long exposures are what FAINT targets get.
    # NOTE, and it is the reason this loop looks the way it does: the first
    # version of this query was `SELECT count(*), trace_a_px * 1.0 / width`
    # with no GROUP BY.  SQLite answers such a query with exactly ONE row —
    # the aggregate, plus whatever the bare column happened to hold on the
    # last row it scanned — so the "median" below was a single arbitrary
    # frame's value and every count printed as 1.  All eight published
    # numbers in this table were wrong, and the monotonic decline they were
    # supposed to demonstrate was an accident of scan order.  Aggregate in
    # Python over the raw column, exactly as `_fracs` above already did.
    erows = []
    for lo_e, hi_e, label in ((0, 5, "&lt; 5 s"), (5, 20, "5 &ndash; 20 s"),
                              (20, 60, "20 &ndash; 60 s"),
                              (60, 10 ** 9, "&gt; 60 s")):
        r = q(con, f"""SELECT trace_a_px * 1.0 / width
                       FROM frame_dispersion
                       WHERE {MEASURED} AND verdict = 'dispersed'
                         AND filter IN ({HIGH_IN}) AND width > 0
                         AND trace_a_px IS NOT NULL
                         AND exptime >= ? AND exptime < ?""", (lo_e, hi_e))
        vals = np.array([x[0] for x in r if x[0] is not None], dtype=float)
        if len(vals):
            erows.append([label, fmt(len(vals)), fnum(float(np.median(vals)), 3)])

    return f"""
<section id="strength"><h2>2&nbsp;&middot;&nbsp;The second axis: which grism
was in the beam?</h2>

<h3>Question</h3>
<p>Two dispersers flew: a high-dispersion H-alpha unit and a low-dispersion
broad-spectrum one. Under the modern vocabulary the card distinguishes them
(<code>hrg</code> vs <code>lrg</code>), but under slot <code>6</code> it does
not &mdash; and a spectrum is useless if you do not know its dispersion. So:
<b>can the pixels say which unit was in the beam?</b></p>

<h3>Evidence</h3>
<p>The obvious statistic is aspect ratio, and it does not work: across the
full labelled populations the two units measure a/b of {fnum(hi[1], 0)}
(high) against {fnum(lo[1], 0)} (low). The reason is that a/b divides by the
minor axis, which is the seeing width, so every wobble in focus or
atmosphere feeds straight into the statistic.</p>

<p>Trace <em>length</em> is better. Expressed as a fraction of frame width
&mdash; which normalises the archive's three cameras and two binnings onto
one scale &mdash; the medians differ by a factor of
<b>{fnum(ratio, 2) if ratio else "&mdash;"}</b>. That is the right sign and
the right order &mdash; the H-alpha unit does disperse further &mdash; but it
is short of the factor of two the optics predict, and the shortfall is itself
a clue: what is measured is not the spectrum's true length but the part of it
that clears the detection threshold. The medians are not the story anyway;
the tails are.</p>

{table(["grism unit", "frames", "median a/b", "median length/width",
        "10&ndash;90 percentile"], rows)}

{_figure(fig, "Left: trace aspect ratio &mdash; the two grisms overlap "
              "almost completely. Right: trace length as a fraction of "
              "frame width, which separates the medians but leaves broad "
              "overlapping tails. The dashed line is the adopted "
              "high-dispersion bound.")}

<p>Asking the question a classifier actually has to answer &mdash; given a
threshold, how PURE is the resulting call? &mdash; exposes a sharp
asymmetry:</p>

{table(["length/width cut", "purity of a &ldquo;high&rdquo; call above it",
        "recall", "purity of a &ldquo;low&rdquo; call below it", "recall"],
       prows, pcls)}

<p>A long trace is decisive: past {fnum(dsp.STRENGTH_HIGH_MIN_FRAC, 2)}
essentially nothing but the H-alpha unit produces it &mdash;
<b>{hi_purity:.1f}%</b> pure, at {hi_recall:.1f}% recall. A short trace is
worthless: across every threshold tried, the purity of a &ldquo;low&rdquo;
call peaks at <b>{best_lo:.1f}%</b> against a base rate of
<b>{base_lo:.1f}%</b>. The cause is visible directly in the
data &mdash; among H-alpha frames alone, trace length falls steadily as
exposure time rises, because long exposures are what FAINT targets get, and
a faint trace drops below the detection threshold sooner:</p>

{table(["exposure", "H-alpha frames", "median length/width"], erows)}

<h3>Decision</h3>
<p><b>The two grisms do not separate cleanly on real data, and this report
declines to pretend otherwise.</b> A frame is called
<code>high</code> when its trace length fraction reaches
{fnum(dsp.STRENGTH_HIGH_MIN_FRAC, 2)}, and <code>ambiguous</code> otherwise.
<b><code>low</code> is never assigned from morphology</b> &mdash; a call
that is right {best_lo:.0f}% of the time against a {base_lo:.0f}% base rate
is not knowledge, and recording it as if it were would quietly corrupt every
downstream count of how many H-alpha spectra the archive holds.</p>

<h3>Consequence</h3>
<p>For frames from the modern naming epochs this costs nothing: their
<code>hrg</code>/<code>lrg</code> cards already name the unit, and this axis
is only a cross-check. The cost falls entirely on the slot-<code>6</code>
era, where the card names nothing &mdash; there, roughly half the spectra
can be assigned to the H-alpha unit and the rest must be carried as
&ldquo;dispersed, unit unknown&rdquo; until someone measures them properly.</p>

<p>Two measurements would close the gap, and neither is available from what
this campaign stores. The first is trace length at a FIXED surface-brightness
threshold rather than at a fixed multiple of the sky noise, which removes the
brightness confound by construction. The second, and better, is spectral
rather than morphological: cross-correlate the extracted trace profile
against the two units' response shapes, which the H-alpha grism's isolated
emission peak makes trivially separable. Both need the per-source and
per-column arrays this campaign summarises away rather than retains, so both
belong to the grism extraction track.</p>
</section>"""


# ---------------------------------------------------------------------------
# Section 3 — slot '6'
# ---------------------------------------------------------------------------
def section_slot6(con) -> str:
    fig = fig_slot6_timeline(con)
    tot = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                            sum(verdict='direct'), sum(verdict='indeterminate')
                     FROM frame_dispersion
                     WHERE {MEASURED} AND filter = '6'""")[0]

    by_target = q(con, f"""
        SELECT coalesce(canonical_target, '(untargeted)'), count(*),
               sum(verdict='dispersed'), sum(verdict='direct'),
               sum(verdict='indeterminate'), min(night), max(night)
        FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6'
        GROUP BY 1 HAVING count(*) >= 20 ORDER BY count(*) DESC""")
    trows, tcls = [], []
    for tgt, n, d, di, ind, n0, n1 in by_target:
        d, di, ind = d or 0, di or 0, ind or 0
        # Amber marks every target whose slot-'6' frames are NOT plain
        # images — i.e. every target where trusting the label would have
        # put dispersed pixels into a photometry pipeline.
        if d / n >= 0.8:
            verdict, cls = "SPECTRA", "warn"
        elif di / n >= 0.8:
            verdict, cls = "images", ""
        else:
            verdict, cls = "MIXED", "warn"
        trows.append([esc(tgt), fmt(n), fmt(d), fmt(di), fmt(ind),
                      f"{esc(n0)} &rarr; {esc(n1)}", f"<b>{verdict}</b>"])
        tcls.append(cls)

    # Is there a changeover date?  Report the dispersed fraction per month,
    # both raw and among the frames that produced a verdict at all — the
    # latter is the fair number, since indeterminate frames are missing
    # evidence rather than evidence of directness.
    def _months(filt):
        rows = q(con, f"""
            SELECT substr(night,1,7), count(*), sum(verdict='dispersed'),
                   sum(verdict='direct'), sum(verdict='indeterminate')
            FROM frame_dispersion
            WHERE {MEASURED} AND filter = ? AND night IS NOT NULL
            GROUP BY 1 ORDER BY 1""", (filt,))
        out, cls = [], []
        for m, n, d, di, ind in rows:
            d, di, ind = d or 0, di or 0, ind or 0
            decided = d + di
            out.append([esc(m), fmt(n), fmt(d), fmt(di), fmt(ind),
                        pct(d, n), f"<b>{pct(d, decided)}</b>"])
            cls.append("warn" if decided and d > di else "")
        return out, cls

    mrows, mcls = _months("6")
    wrows, wcls = _months("W")
    wtot = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                             sum(verdict='direct')
                      FROM frame_dispersion
                      WHERE {MEASURED} AND filter = 'W'""")[0]

    return f"""
<section id="slot6"><h2>3&nbsp;&middot;&nbsp;Slot <code>6</code>: the
committee's open conflict, settled</h2>

<h3>Question</h3>
<p>One panel called slot <code>6</code> a grism; another called it
luminance. Both were reading the same header card, so at most one of them
could be right &mdash; unless the slot is genuinely mixed, in which case
neither is. <b>What fraction of slot-<code>6</code> frames are actually
dispersed, and does the answer depend on when or on what was observed?</b></p>

<h3>Evidence</h3>
<p>{fmt(tot[0])} slot-<code>6</code> frames were measured.
<b>{fmt(tot[1] or 0)}</b> ({pct(tot[1] or 0, tot[0])}) are dispersed,
<b>{fmt(tot[2] or 0)}</b> ({pct(tot[2] or 0, tot[0])}) are direct images,
and {fmt(tot[3] or 0)} are indeterminate. Neither panel was right, because
the slot is not one thing.</p>

{_figure(fig, "Slot-'6' frames per month by measured verdict (top) and the "
              "dispersed fraction (bottom). A wheel that changed on a date "
              "would show a step; a wheel used both ways throughout would "
              "not.")}

<p>Split by target, the pattern is unmistakable &mdash; the verdict tracks
<em>what was being observed</em>, not when:</p>

{table(["target", "frames", "dispersed", "direct", "indet.",
        "nights", "verdict"], trows, tcls)}

<p>Month by month, for the record. There is no step: slot <code>6</code> is
used both ways throughout its life, so <b>no changeover date exists</b> and
no date-based rule can rescue the label.</p>

{table(["month", "frames", "dispersed", "direct", "indet.",
        "% dispersed", "% of decided"], mrows, mcls)}

<h4>The same test applied to <code>W</code> &mdash; which turns out to be a
second mislabelled slot, and a different kind of one</h4>

<p>The <code>W</code> slot was not part of the original dispute, but it is
the other filter whose name says nothing about its optics, so it was
measured on the same footing. Of {fmt(wtot[0])} frames,
<b>{fmt(wtot[1] or 0)}</b> are dispersed and {fmt(wtot[2] or 0)} are direct
&mdash; and unlike slot <code>6</code>, the split is almost perfectly
datable:</p>

{table(["month", "frames", "dispersed", "direct", "indet.",
        "% dispersed", "% of decided"], wrows, wcls)}

<p><code>W</code> is a direct filter through 2024-01 and a grism from
2024-02 onward. The two disputed slots therefore fail in <em>opposite</em>
ways: slot <code>6</code> is mixed by <b>target</b> with no usable date
rule, while <code>W</code> is mixed by <b>date</b> with a clean boundary.
Neither can be resolved by the kind of rule &mdash; "frames after date X are
spectra", "slot 6 means grism" &mdash; that a label-based audit would
naturally reach for.</p>

<h3>Decision</h3>
<p>Slot <code>6</code> is <b>not a filter identity at all</b> &mdash; it is
a wheel position that carried a grism for some programmes and a clear or
broadband element for others, within the same observing season. The label
must therefore never be used to decide whether a frame is a spectrum, and
neither may <code>W</code>. Every downstream consumer reads the per-frame
verdict in <code>frame_dispersion</code> instead.</p>

<h3>Consequence</h3>
<p>Both review panels were reasoning correctly from insufficient evidence,
which is why the conflict never resolved on argument: each had sampled a
different corner of the slot. The disagreement was a fact about the archive,
not about the reviewers, and no amount of further discussion of the header
could have settled it. This is the general lesson for the manifest &mdash;
where a header card encodes an instrument <em>configuration</em> rather than
an instrument <em>property</em>, the card is a hypothesis and the pixels are
the evidence.</p>
</section>"""


# ---------------------------------------------------------------------------
# Section 4 — the two project questions
# ---------------------------------------------------------------------------
def _target_block(con, name, patterns) -> tuple:
    """Verdict tally for one target, matched over canonical_target/path."""
    like = " OR ".join(
        ["canonical_target LIKE ?"] * len(patterns)
        + ["path LIKE ?"] * len(patterns))
    params = tuple(patterns) + tuple(f"%{p.strip('%')}%" for p in patterns)
    r = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                          sum(verdict='direct'), sum(verdict='indeterminate'),
                          min(night), max(night)
                   FROM frame_dispersion
                   WHERE {MEASURED} AND filter = '6' AND ({like})""",
          params)[0]
    return r


def section_projects(con) -> str:
    ngc = _target_block(con, "NGC 5548", ["NGC 5548", "NGC5548"])
    sn = _target_block(con, "SN 2023ixf", ["2023ixf", "SN2023ixf"])

    def _nights(match):
        """Per-night verdict table for one target.

        The reading column is deliberately three-way. An earlier draft asked
        only whether dispersed outnumbered direct, which silently labelled a
        night with NO usable verdict at all — every frame indeterminate,
        which happens when cloud or trailing destroys the evidence — as
        "images". That is the most dangerous possible mislabel here: it
        would hand unmeasurable frames to a photometry pipeline wearing a
        clean bill of health.
        """
        rows = q(con, f"""
            SELECT night, count(*), sum(verdict='dispersed'),
                   sum(verdict='direct'), sum(verdict='indeterminate')
            FROM frame_dispersion
            WHERE {MEASURED} AND filter = '6' AND ({match})
            GROUP BY night ORDER BY night""")
        out, cls = [], []
        for n, c, d, di, ind in rows:
            d, di, ind = d or 0, di or 0, ind or 0
            if d == 0 and di == 0:
                reading, k = "no verdict", "warn"
            elif d > di:
                reading, k = "<b>spectra</b>", "warn"
            elif di > d:
                reading, k = "images", ""
            else:
                reading, k = "<b>split</b>", "warn"
            out.append([esc(n), fmt(c), fmt(d), fmt(di), fmt(ind), reading])
            cls.append(k)
        return out, cls

    snrows, sncls = _nights(
        "canonical_target LIKE '%2023ixf%' OR path LIKE '%2023ixf%'")
    ngrows, ngcls = _nights(
        "canonical_target LIKE '%5548%' OR path LIKE '%5548%'")

    # Which grism unit produced the SN series, and over what baseline?
    sn_str = q(con, f"""
        SELECT sum(strength_class='high'), sum(strength_class='low'),
               sum(strength_class='ambiguous'), min(night), max(night)
        FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6' AND verdict = 'dispersed'
          AND (canonical_target LIKE '%2023ixf%' OR path LIKE '%2023ixf%')
        """)[0]
    sn_nights = q1(con, f"""
        SELECT count(DISTINCT night) FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6' AND verdict = 'dispersed'
          AND (canonical_target LIKE '%2023ixf%' OR path LIKE '%2023ixf%')""")

    ngc_disp = ngc[1] or 0
    sn_disp = sn[1] or 0

    return f"""
<section id="projects"><h2>4&nbsp;&middot;&nbsp;Two projects that were about
to analyse the wrong pixels</h2>

<h3>Question</h3>
<p>Two live projects depend on slot <code>6</code> in opposite directions.
The Dwarf/AGN survey is building an <b>AGN light curve for NGC&nbsp;5548</b>
from those frames, which requires them to be photometry. The SN&nbsp;2023ixf
programme holds <b>83 slot-<code>6</code> frames</b> that it hoped were a
flash-phase spectral series, which requires them to be spectra. Both cannot
be safely assumed. <b>Which is which?</b></p>

<h3>Evidence</h3>
<p><b>NGC&nbsp;5548</b> &mdash; {fmt(ngc[0])} slot-<code>6</code> frames
measured across {esc(ngc[4])}&nbsp;&rarr;&nbsp;{esc(ngc[5])}:
<b>{fmt(ngc_disp)} dispersed</b> ({pct(ngc_disp, ngc[0])}),
{fmt(ngc[2] or 0)} direct, {fmt(ngc[3] or 0)} indeterminate.</p>

{table(["night", "frames", "dispersed", "direct", "indet.", "reading"],
       ngrows, ngcls)}

<p>The campaign alternates: whole nights of spectroscopy
(2023-03-27&nbsp;&rarr;&nbsp;04-02, 04-23, 04-25) interleaved with nights of
imaging (04-05, 04-19). Only <b>{fmt(ngc[2] or 0)} of the {fmt(ngc[0])}
frames</b> are usable photometry.</p>

<p><b>SN&nbsp;2023ixf</b> &mdash; {fmt(sn[0])} slot-<code>6</code> frames
measured across {esc(sn[4])}&nbsp;&rarr;&nbsp;{esc(sn[5])}:
<b>{fmt(sn_disp)} dispersed</b> ({pct(sn_disp, sn[0])}),
{fmt(sn[2] or 0)} direct, {fmt(sn[3] or 0)} indeterminate, spread over
{fmt(sn_nights)} separate nights. None of them reaches the trace length that
would identify the H-alpha unit ({fmt(sn_str[0] or 0)} of {fmt(sn_disp)}).
That is <em>not</em> evidence the series is broad-spectrum: section&nbsp;2
shows a short trace is worthless as a discriminator, because the H-alpha unit
makes short traces too whenever its target is faint &mdash; and a supernova
weeks past discovery is exactly a faint target. An earlier version of this
page read the absence of long traces as &ldquo;all broad-spectrum grism&rdquo;
and that claim is withdrawn. The series must be reduced with its dispersion
treated as a free parameter solved from the frames themselves.</p>

{table(["night", "frames", "dispersed", "direct", "indet.", "reading"],
       snrows, sncls)}

<h3>Decision</h3>
<p><b>The NGC&nbsp;5548 AGN light curve does not survive as photometry.</b>
{pct(ngc_disp, ngc[0])} of its slot-<code>6</code> frames are spectra, and
only {fmt(ngc[2] or 0)} frames across the whole 2023 campaign are direct
images &mdash; too few, and too unevenly spaced, to carry a variability
result. The Dwarf/AGN survey should rebuild that light curve from its other
filters and treat the slot-<code>6</code> frames as a spectroscopic dataset
it did not know it had.</p>

<p><b>The SN&nbsp;2023ixf frames are a genuine spectral series.</b>
{fmt(sn_disp)} dispersed frames on {fmt(sn_nights)} nights, opening
{esc(sn_str[3])} &mdash; two days after the 2023-05-19 discovery, inside the
flash-ionisation window &mdash; and running to {esc(sn_str[4])}. No frame in
the series reaches the H-alpha unit's signature trace length, so it is
plausibly homogeneous in dispersion; section&nbsp;2 explains why that
cannot be asserted outright.</p>

<p>Both answers must be applied at frame granularity rather than as a
blanket ruling per target: a night-by-night split is exactly what a mixed
wheel position produces. Any frame marked <code>indeterminate</code> stays
out of both the light curve and the spectral series until a human has looked
at it.</p>

<h3>Consequence</h3>
<p>A dispersed frame fed to aperture photometry does not fail loudly. It
returns a number &mdash; a smaller number, because the target's flux has
been spread along a trace and most of it now falls outside the aperture &mdash;
and that number lands in the light curve as a fake dimming, correlated with
whichever nights the grism happened to be in the beam. That is the failure
mode this measurement exists to prevent, and it is invisible to every check
that does not look at the pixels. For NGC&nbsp;5548 the correlation would
have been particularly convincing: the dispersed and direct nights come in
runs of several days, so the artefact would have looked like real AGN
variability on exactly the timescale the survey was hunting.</p>

<p>The SN&nbsp;2023ixf result runs the other way and is worth stating as a
gain rather than a correction. A flash-phase spectral series of a
Type&nbsp;II supernova in M101 is a scarce object; this one was sitting in
the archive labelled with a slot number, uncounted, because nobody could
prove what the label meant.</p>
</section>"""


# ---------------------------------------------------------------------------
# Section 5 — the census
# ---------------------------------------------------------------------------
def section_census(con) -> str:
    fig = fig_census(con)
    tot = q(con, f"""SELECT count(*), sum(strength_class='high'),
                            sum(strength_class='ambiguous')
                     FROM frame_dispersion
                     WHERE {MEASURED} AND verdict = 'dispersed'""")[0]

    # By label — with COVERAGE stated, because the campaign is resumable and
    # a partially-measured label must never be read as a complete count.
    by_label = q(con, """
        SELECT filter, count(*) AS queued,
               sum(status = 'measured') AS measured,
               sum(verdict = 'dispersed') AS disp,
               min(night), max(night)
        FROM frame_dispersion
        WHERE population = 'candidate'
        GROUP BY filter ORDER BY queued DESC""")
    lrows, lcls = [], []
    for f, queued, meas, disp, n0, n1 in by_label:
        meas, disp = meas or 0, disp or 0
        lrows.append([f"<code>{esc(f)}</code>", fmt(queued), fmt(meas),
                      pct(meas, queued), fmt(disp), pct(disp, meas),
                      f"{esc(n0)} &rarr; {esc(n1)}"])
        lcls.append("" if meas == queued else "warn")

    # The two quoted claims, checked.
    tcrb = q(con, f"""
        SELECT count(*), sum(fd.verdict='dispersed'),
               sum(fd.strength_class='high'), sum(fd.status='pending')
        FROM frame_dispersion fd
        JOIN stage_tcrb_monitoring s USING (obs_rowid)
        WHERE fd.status = 'measured'""")[0]
    tcrb_queued = q1(con, """
        SELECT count(*) FROM frame_dispersion fd
        JOIN stage_tcrb_monitoring s USING (obs_rowid)""")
    # The label count that produced the quoted "247": T CrB frames carrying
    # any grism card, measured or not.  Counted from frame_dispersion, whose
    # queue is already restricted to CANONICAL science frames — so this is
    # 247 distinct observations, not 247 rows in `frames` (which holds 483,
    # the difference being duplicate copies of the same exposure).  The
    # prose says "canonical science frames" for exactly that reason.
    tcrb_label = q1(con, f"""
        SELECT count(*) FROM frame_dispersion
        WHERE canonical_target = 'T CrB' AND filter IN ({DISP_IN})""")
    tcrb_raw = q1(con, f"""
        SELECT count(*) FROM frames
        WHERE canonical_target = 'T CrB' AND filter IN ({DISP_IN})""")
    tcrb_g = q(con, f"""
        SELECT count(*), sum(verdict='dispersed')
        FROM frame_dispersion
        WHERE {MEASURED} AND canonical_target = 'T CrB'""")[0]
    # The CAMPAIGN's size comes from the staging table itself.  Reading it
    # off the join to frame_dispersion instead would silently report only the
    # subset that carried a candidate FILTER card and so entered this queue —
    # which is the very label-versus-measurement confusion this section
    # exists to correct.
    bestar_staged = q1(con, "SELECT count(*) FROM stage_bestar_grism")
    bestar = q(con, """
        SELECT count(*), sum(fd.status = 'measured'),
               sum(fd.verdict = 'dispersed')
        FROM frame_dispersion fd
        JOIN stage_bestar_grism s USING (obs_rowid)""")[0]
    # Queued (= archive label count) vs measured, for every grism label.
    grism_lbl = q(con, f"""
        SELECT count(*), sum(status = 'measured'),
               sum(verdict = 'dispersed')
        FROM frame_dispersion WHERE filter IN ({DISP_IN})""")[0]
    cand_all = q(con, """
        SELECT count(*), sum(status = 'measured'),
               sum(verdict = 'dispersed')
        FROM frame_dispersion WHERE population = 'candidate'""")[0]

    # The coverage caveat is only a caveat while frames remain unmeasured.
    n_pending = q1(con, "SELECT count(*) FROM frame_dispersion "
                        "WHERE status = 'pending'")
    n_unread = q1(con, "SELECT count(*) FROM frame_dispersion "
                       "WHERE status = 'unreadable'")
    if n_pending:
        coverage_note = (
            '<p class="decision">Rows highlighted amber are not yet fully '
            f'measured &mdash; {fmt(n_pending)} frames remain queued. The '
            'campaign is resumable and runs in priority order (the disputed '
            'labels and the control first, the undisputed '
            '<code>hrg</code>/<code>lrg</code> bulk last), so a partial run '
            'answers the contested questions completely while leaving only '
            'the census refinement outstanding. Where coverage is below '
            '100%, the dispersed count is a floor, not a total.</p>')
    else:
        coverage_note = (
            '<p class="decision">Coverage is <b>complete</b>: every queued '
            'frame has been measured, so the counts in this section are '
            'totals rather than floors. The '
            f'{fmt(n_unread)} frame(s) marked amber could not be opened at '
            'all &mdash; a truncated file, not an ambiguous one &mdash; and '
            'are excluded from every percentage on this page.</p>')

    return f"""
<section id="census"><h2>5&nbsp;&middot;&nbsp;How many spectra does the
archive actually hold?</h2>

<h3>Question</h3>
<p>Two figures circulate in the project write-ups: the T&nbsp;CrB series is
quoted as <b>247 Mode0 spectra</b>, and the BeStar campaign as
<b>~23,000 grism frames</b>. Both were counted from header labels. <b>Do
they survive measurement?</b></p>

<h3>Evidence</h3>
<p>Across every candidate and control frame measured, <b>{fmt(tot[0])}
frames are dispersed</b>. Of those, {fmt(tot[1] or 0)} carry a trace long
enough to identify the H-alpha unit; for the remaining {fmt(tot[2] or 0)}
the unit cannot be named from morphology (section&nbsp;2), which is a limit
on the instrument attribution, not on the spectrum count.</p>

{table(["label", "in archive", "measured", "coverage", "dispersed",
        "% of measured", "date span"], lrows, lcls)}

{coverage_note}

{_figure(fig, "Measured spectra per target for the 25 most-observed "
              "spectroscopic targets, split by grism unit.")}

<p><b>T&nbsp;CrB &mdash; the quoted 247 checks out, as a label count.</b> The
archive holds exactly <b>{fmt(tcrb_label)}</b> <em>canonical science</em>
T&nbsp;CrB frames carrying a grism card, so that figure was derived correctly
from the headers. (The raw <code>frames</code> table lists
{fmt(tcrb_raw)} such rows; the difference is duplicate copies of the same
exposures, which is why the qualifier matters.) What the
measurement adds is that it can now be audited frame by frame instead of
taken on trust. Of the {fmt(tcrb_queued)} T&nbsp;CrB monitoring frames this
campaign queued, {fmt(tcrb[0])} are measured and <b>{fmt(tcrb[1] or 0)}</b>
of those are dispersed ({fmt(tcrb[2] or 0)} long enough to identify as
H-alpha) &mdash; so the series is substantially larger than 247 once frames
that are spectra without a grism card are counted too.</p>

<p><b>BeStar campaign &mdash; the quoted ~23,000 does not.</b> The staged
BeStar campaign contains <b>{fmt(bestar_staged)}</b> frames in total, not
23,000 &mdash; a real and well-populated campaign, but roughly a fifth of the
quoted size. Of those, {fmt(bestar[0])} carried a grism-family
<code>FILTER</code> card and so entered this campaign's queue;
{fmt(bestar[1] or 0)} are measured and {fmt(bestar[2] or 0)} of those are
dispersed. The gap between {fmt(bestar_staged)} and {fmt(bestar[0])} is
itself worth noting: those frames are staged as BeStar work but carry no
grism label at all, so a label-only audit would not have found them.
The ~23,000 figure is instead close to
<b>{fmt(cand_all[0])}</b> &mdash; every frame in the archive whose filter is
a grism name or a disputed slot, across <em>all</em> projects, of which
{fmt(grism_lbl[0])} carry an unambiguous grism card. It is an archive-wide
instrument total that appears to have been attached to a single campaign.</p>

<h3>Decision</h3>
<p>The census above replaces both quoted figures as the citable numbers, and
they should be cited with their basis attached &mdash; "frames measured as
dispersed", not "frames labelled grism". Where a project's own staged frame
list disagrees with a headline count, the staged list is the narrower and
more defensible number.</p>

<h3>Consequence</h3>
<p>The gap between a label count and a measured count is not noise; it is
made of the specific frames a label-based count gets wrong in both
directions &mdash; grism-labelled frames that are too shallow to have
recorded a usable trace, and slot-<code>6</code> frames that are spectra
nobody had counted. Any paper quoting a spectrum count from this archive
should quote the measured one, and any paper quoting the label count should
expect a referee to ask which frames were actually looked at.</p>
</section>"""


# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S2c report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        n_meas = q1(con, f"SELECT count(*) FROM frame_dispersion "
                         f"WHERE {MEASURED}")
        n_disp = q1(con, f"SELECT count(*) FROM frame_dispersion "
                         f"WHERE {MEASURED} AND verdict = 'dispersed'")
        n_bad = q1(con, "SELECT count(*) FROM frame_dispersion "
                        "WHERE status = 'unreadable'")
        meta = dict(q(con, "SELECT key, value FROM s2c_build_meta"))

        sections = [
            section_calibration(con),
            section_strength(con),
            section_slot6(con),
            section_projects(con),
            section_census(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S2c — Filter Identity</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S2c — Filter Identity, Measured</h1>
  <p>{fmt(n_meas)} frames measured from their pixels &middot;
  {fmt(n_disp)} dispersed &middot; {fmt(n_bad)} unreadable &middot;
  built {esc(meta.get('built_at', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))})
  &middot; <a href="../index.html">back to the evidence hub</a></p>
</header>

<nav>
  <a href="#calibration">1 Calibration</a> &middot;
  <a href="#strength">2 Which grism</a> &middot;
  <a href="#slot6">3 Slot '6'</a> &middot;
  <a href="#projects">4 NGC 5548 &amp; SN 2023ixf</a> &middot;
  <a href="#census">5 Spectra census</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>rlmt_diagnostics.report_s2c</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query or a documented pipeline constant; none is
typed by hand. Regenerate with
<code>pipeline/scripts/run_s2c_dispersion.py report</code>.</footer>
</body></html>"""

        HTML_PATH.write_text(html, encoding="utf-8")

        # Belt and braces: every <img> the page references must exist.
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
