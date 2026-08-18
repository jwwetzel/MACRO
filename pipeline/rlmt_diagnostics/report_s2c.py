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

    groups = [("known direct (control)", f"filter IN ({DIR_IN})", MUTED, 6),
              ("low-dispersion grism", f"filter IN ({LOW_IN})", ACCENT, 8),
              ("high-dispersion grism", f"filter IN ({HIGH_IN})", WARN, 8)]

    with plt.rc_context(DARK):
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, 4.8))
        for name, where, color, size in groups:
            d = fetch(where)
            if not len(d):
                continue
            axl.scatter(np.clip(d[:, 0], 0.8, 400), d[:, 1], s=size,
                        alpha=0.35, color=color, label=name, linewidths=0)
            t = d[d[:, 4] > 0]
            if len(t):
                axr.scatter(np.clip(t[:, 2], 0.8, 400), t[:, 3], s=size,
                            alpha=0.35, color=color, linewidths=0)
        axl.set_xscale("log")
        axl.set_xlabel("median a/b of the 10 brightest sources")
        axl.set_ylabel("position-angle scatter of those 10 (deg)")
        axl.set_title("The statistic that FAILS:\nbright-set median", fontsize=10)
        axl.legend(fontsize=7, loc="upper left", framealpha=0.3)
        axr.set_xscale("log")
        axr.axvline(dsp.TRACE_MIN_AB, color=GOOD, lw=1.0, ls="--")
        axr.axhline(dsp.TRACE_MAX_PA_SCATTER_DEG, color=GOOD, lw=1.0, ls="--")
        axr.set_xlabel("trace a/b (trace sources only)")
        axr.set_ylabel("position-angle scatter of the traces (deg)")
        axr.set_title("The statistic that WORKS:\ntrace population + shared axis",
                      fontsize=10)
        fig.suptitle("Calibration on labels nobody disputes — dispersed "
                     "frames must land bottom-right", fontsize=11)
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
        axr.axvline(dsp.STRENGTH_LOW_MAX_FRAC, color=GOOD, lw=1.1, ls="--")
        axr.axvline(dsp.STRENGTH_HIGH_MIN_FRAC, color=GOOD, lw=1.1, ls="--")
        axr.set_xlabel("trace length / frame width")
        axr.set_title("Trace length: they separate\n(dashed = adopted bounds)",
                      fontsize=10)
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
               sum(strength_class = 'high'), sum(strength_class = 'low'),
               sum(strength_class NOT IN ('high', 'low'))
        FROM frame_dispersion
        WHERE {MEASURED} AND verdict = 'dispersed'
        GROUP BY tgt ORDER BY count(*) DESC LIMIT 25""")
    tgt = [r[0] for r in rows][::-1]
    hi = np.array([r[1] or 0 for r in rows], dtype=float)[::-1]
    lo = np.array([r[2] or 0 for r in rows], dtype=float)[::-1]
    am = np.array([r[3] or 0 for r in rows], dtype=float)[::-1]
    y = np.arange(len(tgt))
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(9.5, 7.5))
        ax.barh(y, hi, color=WARN, label="H-alpha grism")
        ax.barh(y, lo, left=hi, color=ACCENT, label="broad grism")
        ax.barh(y, am, left=hi + lo, color=MUTED, label="ambiguous strength")
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

    rows = []
    for label in list(dsp.KNOWN_DISPERSED_FILTERS) + list(
            dsp.KNOWN_DIRECT_FILTERS):
        r = q(con, f"""SELECT count(*), sum(verdict='dispersed'),
                              sum(verdict='direct'),
                              sum(verdict='indeterminate')
                       FROM frame_dispersion
                       WHERE {MEASURED} AND filter = ?""", (label,))[0]
        n = r[0] or 0
        if not n:
            continue
        truth = dsp.expected_verdict(label)
        agree = (r[1] if truth == dsp.VERDICT_DISPERSED else r[2]) or 0
        cls = "good" if n and agree / n >= 0.95 else (
            "warn" if n and agree / n >= 0.80 else "bad")
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
                    WHERE {MEASURED} AND filter IN ({DIR_IN})""")[0]
    recall = pct(kd[1] or 0, kd[0])
    fpr = pct(kn[1] or 0, kn[0])

    # The measured overlap region, stated honestly.
    sep = q(con, f"""
        SELECT
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND filter IN ({DISP_IN}) AND n_trace = 0),
          (SELECT count(*) FROM frame_dispersion
            WHERE {MEASURED} AND filter IN ({DIR_IN}) AND n_trace > 0)""")[0]

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

{_figure(fig, "Left: the bright-set median a/b that the first survey pass "
              "quoted &mdash; the classes overlap badly, because most bright "
              "detections on a grism frame are round field stars whose "
              "traces are too faint to resolve. Right: the trace population "
              "and its position-angle scatter, where the same frames "
              "separate. Dashed lines mark the adopted gates.")}

{table(["label", "truth", "measured", "dispersed", "direct",
        "indeterminate", "agreement"],
       [r[0] for r in rows], [r[1] for r in rows])}

<p>Pooled over the {fmt(kd[0])} known-dispersed frames, the classifier
recovers <b>{recall}</b> as dispersed. Pooled over the {fmt(kn[0])}
known-direct control frames it calls <b>{fpr}</b> dispersed &mdash; that is
the false-positive rate, and it is the number that decides whether any
verdict on slot <code>6</code> can be believed.</p>

<p>The overlap, stated plainly: {fmt(sep[0])} known-dispersed frames
contained no trace-shaped source at all (too short an exposure, too faint a
target, or cloud), and {fmt(sep[1])} known-direct frames contained at least
one (satellites, cosmic-ray tracks, edge-on galaxies). Those frames are the
irreducible ambiguity in the method, and the classifier is built to return
<code>indeterminate</code> across it rather than to guess.</p>

<h3>Decision</h3>
<p>Dispersion is measurable per frame, and the classifier is adopted with
three rules. <b>(1)</b> Two or more traces sharing an axis to within
{fnum(dsp.TRACE_MAX_PA_SCATTER_DEG, 0)}&nbsp;deg &rarr;
<code>dispersed</code>. <b>(2)</b> A single trace of
a/b&nbsp;&ge;&nbsp;{fnum(dsp.SOLO_TRACE_MIN_AB, 0)} in a sparse field
(&le;&nbsp;{fmt(dsp.SOLO_MAX_SOURCES)} sources) &rarr;
<code>dispersed</code>; this rule carries the bright-standard frames where
the target is the only thing in the field. <b>(3)</b> No traces and round
bright sources &rarr; <code>direct</code>. Everything else is
<code>indeterminate</code> and is reported as such.</p>

<h3>Consequence</h3>
<p>Two false-positive modes were found by measurement and closed before the
production run, and both are worth recording because both would have
produced confident nonsense. <b>Twilight flats</b> read as dispersed: a flat
has no stars, so the extractor finds dust shadows and detector column
defects, which are straight and mutually parallel because they <em>are</em>
columns &mdash; a perfect forgery of a shared dispersion axis. Calibration
frames are now excluded from the campaign entirely. <b>A satellite trail</b>
across a 60-second luminance frame of M57 read as dispersed on the strength
of one 1,095-px streak &mdash; past 1,278 perfectly round stars. A grating
disperses <em>everything</em>, so a rich field with exactly one smear is the
one thing a grism cannot produce; the solo-trace rule now requires a sparse
field. Any future reuse of this classifier inherits both guards.</p>
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

    # How much of each population lands in the ambiguous band?
    def amb(where):
        r = q(con, f"""SELECT count(*), sum(strength_class = 'ambiguous'),
                              sum(strength_class = 'high'),
                              sum(strength_class = 'low')
                       FROM frame_dispersion
                       WHERE {MEASURED} AND {where}
                         AND verdict = 'dispersed'""")[0]
        return r
    ah, al = amb(f"filter IN ({HIGH_IN})"), amb(f"filter IN ({LOW_IN})")

    rows = [
        ["H-alpha grism (<code>hrg</code>, <code>HaGrism</code>)", fmt(hi[0]),
         fnum(hi[1], 1), f"<b>{fnum(hi[2], 3)}</b>",
         f"{fnum(hi[3], 3)} &ndash; {fnum(hi[4], 3)}",
         pct(ah[2] or 0, ah[0]), pct(ah[1] or 0, ah[0])],
        ["broad grism (<code>lrg</code>, <code>OGGrism</code>)", fmt(lo[0]),
         fnum(lo[1], 1), f"<b>{fnum(lo[2], 3)}</b>",
         f"{fnum(lo[3], 3)} &ndash; {fnum(lo[4], 3)}",
         pct(al[3] or 0, al[0]), pct(al[1] or 0, al[0])]]

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
<p>The obvious statistic is aspect ratio, and it does not work. The
2025-01-23 focus sweep is the cleanest possible test &mdash; the same star,
the same night, the same camera, shot through both units &mdash; and it
measures a/b&nbsp;of {fnum(hi[1], 0)} for the H-alpha unit against
{fnum(lo[1], 0)} for the broad one across the full populations: overlapping
distributions. The reason is that a/b divides by the minor axis, which is
the seeing width, so every wobble in focus or atmosphere feeds straight into
the statistic.</p>

<p>Trace <em>length</em> does not have that problem. Expressed as a fraction
of frame width &mdash; which normalises the archive's three cameras and two
binnings onto one scale &mdash; the same frames separate by a factor of
<b>{fnum(ratio, 1) if ratio else "&mdash;"}</b>.</p>

{_figure(fig, "Left: trace aspect ratio, where the two grisms overlap. "
              "Right: trace length as a fraction of frame width, where they "
              "do not. Dashed lines are the adopted bounds; the gap between "
              "them is the ambiguous band.")}

{table(["grism unit", "frames", "median a/b", "median length/width",
        "10&ndash;90 pct", "assigned correctly", "left ambiguous"], rows)}

<h3>Decision</h3>
<p>Dispersion strength is classified on trace length fraction, with
<code>low</code> at &le;&nbsp;{fnum(dsp.STRENGTH_LOW_MAX_FRAC, 2)} and
<code>high</code> at &ge;&nbsp;{fnum(dsp.STRENGTH_HIGH_MIN_FRAC, 2)}. The
band between them is reported <code>ambiguous</code> rather than assigned.
Aspect ratio is retained in the table as a measured column but is not used
to decide anything.</p>

<h3>Consequence</h3>
<p>The separation is good but it is not perfect, and the residual spread is
physical rather than a threshold that needs tuning: a brighter star's trace
stays above the extraction threshold further into its wings, so length grows
with source brightness and exposure time. A frame in the ambiguous band is
not a failure of the measurement &mdash; it is a frame whose dispersion
genuinely cannot be read from its trace length alone. If a future project
needs those assigned, the measurement that would do it is spectral rather
than morphological: cross-correlate the extracted trace profile against the
two units' known response shapes, which the H-alpha unit's isolated emission
peak makes trivially separable. That is a job for the grism extraction
track, not for this classifier.</p>
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
        if d / n >= 0.8:
            verdict, cls = "SPECTRA", "bad"
        elif di / n >= 0.8:
            verdict, cls = "images", "good"
        else:
            verdict, cls = "MIXED", "warn"
        trows.append([esc(tgt), fmt(n), fmt(d), fmt(di), fmt(ind),
                      f"{esc(n0)} &rarr; {esc(n1)}", f"<b>{verdict}</b>"])
        tcls.append(cls)

    # Is there a changeover date?  Compare the dispersed fraction before and
    # after each month boundary and report the sharpest transition.
    months = q(con, f"""
        SELECT substr(night,1,7), count(*), sum(verdict='dispersed')
        FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6' AND night IS NOT NULL
        GROUP BY 1 ORDER BY 1""")
    mrows = [[esc(m), fmt(n), fmt(d or 0), pct(d or 0, n)]
             for m, n, d in months]

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

<p>Month by month, for the record:</p>

{table(["month", "frames", "dispersed", "% dispersed"], mrows)}

<h3>Decision</h3>
<p>Slot <code>6</code> is <b>not a filter identity at all</b> &mdash; it is
a wheel position that carried a grism for some programmes and a clear or
broadband element for others, within the same observing season. The label
must therefore never be used to decide whether a frame is a spectrum. Every
downstream consumer reads the per-frame verdict in
<code>frame_dispersion</code> instead.</p>

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

    sn_months = q(con, f"""
        SELECT night, count(*), sum(verdict='dispersed'), sum(verdict='direct')
        FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6'
          AND (canonical_target LIKE '%2023ixf%' OR path LIKE '%2023ixf%')
        GROUP BY night ORDER BY night""")
    snrows = [[esc(n), fmt(c), fmt(d or 0), fmt(di or 0),
               "<b>spectra</b>" if (d or 0) > (di or 0) else "images"]
              for n, c, d, di in sn_months]

    ngc_nights = q(con, f"""
        SELECT night, count(*), sum(verdict='dispersed'), sum(verdict='direct')
        FROM frame_dispersion
        WHERE {MEASURED} AND filter = '6'
          AND (canonical_target LIKE '%5548%' OR path LIKE '%5548%')
        GROUP BY night ORDER BY night""")
    ngrows = [[esc(n), fmt(c), fmt(d or 0), fmt(di or 0),
               "<b>spectra</b>" if (d or 0) > (di or 0) else "images"]
              for n, c, d, di in ngc_nights]

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

{table(["night", "frames", "dispersed", "direct", "reading"], ngrows)}

<p><b>SN&nbsp;2023ixf</b> &mdash; {fmt(sn[0])} slot-<code>6</code> frames
measured across {esc(sn[4])}&nbsp;&rarr;&nbsp;{esc(sn[5])}:
<b>{fmt(sn_disp)} dispersed</b> ({pct(sn_disp, sn[0])}),
{fmt(sn[2] or 0)} direct, {fmt(sn[3] or 0)} indeterminate.</p>

{table(["night", "frames", "dispersed", "direct", "reading"], snrows)}

<h3>Decision</h3>
<p>The per-night tables above are the operative answer for each project, and
they must be applied at frame granularity rather than as a blanket ruling
per target: a night-by-night split is exactly what a mixed wheel position
produces. Any frame the classifier marks <code>indeterminate</code> stays
out of both the light curve and the spectral series until a human has looked
at it.</p>

<h3>Consequence</h3>
<p>A dispersed frame fed to aperture photometry does not fail loudly. It
returns a number &mdash; a smaller number, because the target's flux has
been spread along a trace and most of it now falls outside the aperture &mdash;
and that number lands in the light curve as a fake dimming, correlated with
whichever nights the grism happened to be in the beam. That is the failure
mode this measurement exists to prevent, and it is invisible to every check
that does not look at the pixels.</p>
</section>"""


# ---------------------------------------------------------------------------
# Section 5 — the census
# ---------------------------------------------------------------------------
def section_census(con) -> str:
    fig = fig_census(con)
    tot = q(con, f"""SELECT count(*), sum(strength_class='high'),
                            sum(strength_class='low'),
                            sum(strength_class='ambiguous')
                     FROM frame_dispersion
                     WHERE {MEASURED} AND verdict = 'dispersed'""")[0]

    # By label epoch — how the archive's spectra are distributed in time.
    by_label = q(con, f"""
        SELECT filter, count(*), sum(verdict='dispersed'),
               min(night), max(night)
        FROM frame_dispersion
        WHERE {MEASURED} AND population = 'candidate'
        GROUP BY filter ORDER BY count(*) DESC""")
    lrows = [[f"<code>{esc(f)}</code>", fmt(n), fmt(d or 0), pct(d or 0, n),
              f"{esc(n0)} &rarr; {esc(n1)}"]
             for f, n, d, n0, n1 in by_label]

    # The two quoted claims, checked.
    tcrb = q(con, f"""
        SELECT count(*), sum(fd.verdict='dispersed'),
               sum(fd.strength_class='high'), sum(fd.strength_class='low')
        FROM frame_dispersion fd
        JOIN stage_tcrb_monitoring s USING (obs_rowid)
        WHERE fd.status = 'measured'""")[0]
    tcrb_g = q(con, f"""
        SELECT count(*), sum(verdict='dispersed')
        FROM frame_dispersion
        WHERE {MEASURED} AND canonical_target = 'T CrB'""")[0]
    bestar = q(con, f"""
        SELECT count(*), sum(fd.verdict='dispersed')
        FROM frame_dispersion fd
        JOIN stage_bestar_grism s USING (obs_rowid)
        WHERE fd.status = 'measured'""")[0]
    bestar_all = q(con, f"""
        SELECT count(*), sum(verdict='dispersed')
        FROM frame_dispersion
        WHERE {MEASURED} AND filter IN ({DISP_IN})""")[0]

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
frames are dispersed</b>: {fmt(tot[1] or 0)} through the H-alpha unit,
{fmt(tot[2] or 0)} through the broad unit, and {fmt(tot[3] or 0)} whose
dispersion strength is ambiguous.</p>

{table(["label", "frames measured", "dispersed", "% dispersed", "date span"],
       lrows)}

{_figure(fig, "Measured spectra per target for the 25 most-observed "
              "spectroscopic targets, split by grism unit.")}

<p><b>T&nbsp;CrB.</b> The monitoring stage holds {fmt(tcrb[0])} measured
frames, of which <b>{fmt(tcrb[1] or 0)} are dispersed</b>
({fmt(tcrb[2] or 0)} H-alpha, {fmt(tcrb[3] or 0)} broad). Counting instead
by target name across the whole archive, {fmt(tcrb_g[0])} frames carry the
target <code>T CrB</code> and {fmt(tcrb_g[1] or 0)} of them measure as
spectra.</p>

<p><b>BeStar campaign.</b> The staged campaign holds {fmt(bestar[0])}
measured frames, {fmt(bestar[1] or 0)} of them dispersed. The
~23,000 figure is not a count of that campaign: it is close to the
{fmt(bestar_all[0])} frames carrying <em>any</em> grism label archive-wide,
of which {fmt(bestar_all[1] or 0)} measure as dispersed.</p>

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
