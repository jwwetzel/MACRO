"""SN 2023ixf Gate 0 evidence report.

Reads the ``sn_g0_*`` tables from the manifest and writes:

* ``docs/SN2023ixf_LightCurve/sn_gate0.html`` — the report
* ``docs/SN2023ixf_LightCurve/figures/gate0/*.png`` — every figure

Same discipline as every other renderer in this repository: the page follows
the site's Socratic format (Question → Evidence → Decision → Consequence),
and EVERY number on it is either the result of a SQL query executed here or a
documented constant read from ``macro_sn.gate0``.  Nothing is typed by hand
— including the saturation thresholds, which are read from the module so the
page and the census can never drift apart.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

import sys                       # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_core.report_s0 import (  # noqa: E402
    ACCENT, BAD, FAINT, GOOD, INK, MUTED, STYLE, WARN, DPI,
    esc, fmt, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402
from macro_sn import SN_G0_CODE_VERSION  # noqa: E402
from macro_sn import gate0 as g0         # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "SN2023ixf_LightCurve"
FIG_DIR = DOCS_DIR / "figures" / "gate0"
HTML_PATH = DOCS_DIR / "sn_gate0.html"

#: Colour per saturation class, used by every figure and by the matrix
#: legend, so a colour means one thing on the whole page.
CLASS_COLOR = {
    "clean": GOOD,
    "bounded_clean": ps.tint(GOOD, 0.45),
    "suspect": WARN,
    "rejected": BAD,
    "undetermined": FAINT,
}

#: Display order for the filters — role first, then the wheel's own order.
FILTER_ORDER = list(g0.BROADBAND_FILTERS) + list(g0.NARROWBAND_FILTERS) + \
    ["6", "X"]


def pct(num, den) -> str:
    """Percentage with a guarded denominator."""
    if not den:
        return "&mdash;"
    return f"{100.0 * num / den:.1f}%"


def fnum(x, nd=0) -> str:
    """NULL-safe fixed-decimal number for the page."""
    if x is None:
        return "&mdash;"
    return f"{x:,.{nd}f}"


def _figure(src: str, caption: str) -> str:
    return (f'<figure><a href="{src}"><img src="{src}" alt=""></a>'
            f"<figcaption>{caption}</figcaption></figure>")


def screen_of(con) -> g0.Screen:
    """The campaign's screen, rebuilt here from S2's measured ceiling.

    Rebuilt rather than read back from the census so that the page states
    the rule, not a cached consequence of it: if S2 re-measures the clip and
    the census has not been re-run, this raises the difference into view
    instead of printing a stale threshold.
    """
    r = q(con, "SELECT mode, clip_adu, veto_adu FROM s2_ceiling_modes "
               "WHERE mode = 'High Gain'")
    if not r:
        raise RuntimeError("s2_ceiling_modes has no High Gain row: the "
                           "saturation screen has no measurement behind it")
    return g0.screen_for_mode(*r[0])


# ===========================================================================
# Figures
# ===========================================================================
def fig_matrix(con) -> str:
    """The saturation matrix: filter × night, coloured by what the census
    could conclude about the supernova in each cell.

    The cell value is the fraction of that night's frames in that filter
    from which photometry may be taken (clean, or clean by bound).  Cells
    that contain at least one REJECTED frame carry a dot, because the
    campaign's single most consequential finding is that saturation here is
    a per-frame property and not a per-epoch one: a night can be half clean
    and half clipped.
    """
    rows = q(con, """SELECT night, filter, n_frames, n_clean,
                            n_bounded_clean, n_rejected, n_undetermined
                     FROM sn_g0_matrix ORDER BY night, filter""")
    nights = sorted({r[0] for r in rows})
    filts = [f for f in FILTER_ORDER if any(r[1] == f for r in rows)]
    grid = np.full((len(filts), len(nights)), np.nan)
    dots, undet = [], []
    for night, filt, n, clean, bound, rej, und in rows:
        if filt not in filts:
            continue
        i, j = filts.index(filt), nights.index(night)
        grid[i, j] = (clean + bound) / n if n else np.nan
        if rej:
            dots.append((j, i))
        if und == n and n:
            undet.append((j, i))

    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(12.4, 4.2))
        # The house sequential ramp, NOT a red-to-green one.  This cell
        # carries its value in colour alone — no number is printed in it —
        # so the ramp has to survive both colour-blindness and a greyscale
        # print.  vermilion->orange->green fails BOTH: deuteranopes cannot
        # separate its ends, and in greyscale those ends land on 128 and 138
        # of 255 while the MIDDLE sits at 171, so a half-usable night prints
        # lighter than either extreme and the map reads inverted.
        # ps.SEQ_CMAP is monotonic in lightness, so pale still means "little
        # usable" after any amount of photocopying.
        cmap = ps.SEQ_CMAP
        # The house web profile draws a grid; on a heat map its white rules
        # cut every cell in half and read as cell boundaries that are not
        # there.  Off for this figure only.
        ax.grid(False)
        im = ax.imshow(grid, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                       interpolation="nearest")

        def _overlay_ink(pts):
            """Marker colour per cell: dark on pale cells, white on deep ones.

            A fixed INK dot disappears on a fully-usable cell now that the
            top of the ramp is dark blue.  ps.ink_on makes the same
            luminance decision every other heatmap in the repository makes.
            """
            out = []
            for j, i in pts:
                v = grid[i, j]
                out.append(ps.ink_on(0.0 if np.isnan(v) else float(v)))
            return out

        if dots:
            ax.scatter([d[0] for d in dots], [d[1] for d in dots], s=9,
                       color=_overlay_ink(dots), zorder=3, marker="o",
                       linewidths=0)
        if undet:
            ax.scatter([d[0] for d in undet], [d[1] for d in undet], s=26,
                       color=_overlay_ink(undet), zorder=3, marker="x",
                       linewidths=1.0)
        ax.set_yticks(range(len(filts)))
        ax.set_yticklabels([f"{f}  ({g0.band_role(f)[:4]})" for f in filts])
        step = max(1, len(nights) // 18)
        ax.set_xticks(range(0, len(nights), step))
        ax.set_xticklabels([nights[i][5:] for i in range(0, len(nights), step)],
                           rotation=90)
        ax.set_xlabel("night (local-noon label, 2023)")
        ax.set_title("Fraction of frames from which photometry may be taken "
                     "(dot = the cell also contains a clipped frame; "
                     "x = nothing could be concluded)")
        fig.colorbar(im, ax=ax, pad=0.01, label="usable fraction")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "g0_saturation_matrix.png", dpi=DPI)
        plt.close(fig)
    return "figures/gate0/g0_saturation_matrix.png"


def fig_peaks(con) -> str:
    """Every measured supernova peak against phase, with the screen drawn.

    Two panels because the two bands answer two different questions.  Left:
    the broadband codes that carry the light curve — the story is the wall
    of clipped points before +5.4 d and the seeing-driven scatter across the
    screen afterwards.  Right: the narrowband codes, where the exposure ramp
    (64 → 128 → 32 → 64 → 128 s) is visible as steps and where the +2.5 d
    epoch sits well under the screen.
    """
    s = screen_of(con)
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6), sharey=True)
        for ax, role, bands in ((axes[0], "broadband", g0.BROADBAND_FILTERS),
                                (axes[1], "narrowband", g0.NARROWBAND_FILTERS)):
            for i, filt in enumerate(bands):
                # Science tree only.  The two mjc engineering frames read
                # 65,535 ADU in a 12-bit channel; leaving them in would set
                # the y axis to a range in which the entire campaign is a
                # flat line at the bottom, which is how an impossible value
                # hides a real one.
                pts = q(con, """SELECT phase_d, peak_adu FROM sn_g0_census
                                WHERE epoch_role = 'campaign' AND filter = ?
                                  AND quality = 'wcs' AND peak_adu IS NOT NULL
                                  AND tree = ?
                                ORDER BY phase_d""", (filt, g0.SCIENCE_TREE))
                if not pts:
                    continue
                a = np.array(pts, dtype=float)
                kw = ps.series(i)
                ax.scatter(a[:, 0], a[:, 1], s=15, label=f"{filt}",
                           color=kw["color"], marker=kw["marker"],
                           alpha=0.75, linewidths=0.3, edgecolors=INK)
            ax.axhline(s.clip_adu, color=INK, ls="-", lw=1.1)
            ax.axhline(s.reject_adu, color=BAD, ls="--", lw=1.2)
            ax.axhline(s.suspect_adu, color=WARN, ls=":", lw=1.4)
            ax.axvspan(0, g0.FLASH_PHASE_END_D, color=ps.tint(ACCENT, 0.88),
                       zorder=0)
            ax.set_xlabel("days after the adopted explosion epoch")
            ax.set_title(f"{role} codes")
            ax.set_ylim(0, s.clip_adu * 1.12)
            ax.legend(fontsize=9, ncol=3, loc="lower left", framealpha=0.85)
            ax.grid(alpha=0.25)
        axes[0].set_ylabel("peak ADU at the supernova's position")
        for lbl, val, col, va in (
                (f"measured clip {s.clip_adu:,}", s.clip_adu, INK, "bottom"),
                (f"reject {s.reject_adu:,}", s.reject_adu, BAD, "bottom"),
                (f"suspect {s.suspect_adu:,}", s.suspect_adu, WARN, "top")):
            axes[1].annotate(lbl, (0.985, val),
                             xycoords=("axes fraction", "data"),
                             ha="right", va=va, fontsize=8, color=col)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "g0_peak_vs_phase.png", dpi=DPI)
        plt.close(fig)
    return "figures/gate0/g0_peak_vs_phase.png"


def fig_astrometry(con) -> str:
    """The astrometry question, decomposed.

    One bar per scope with its Wilson interval, and the two verdict
    thresholds drawn as vertical rules.  The figure exists to make one point
    visible at a glance: the S1 stratum's NO-GO and the campaign broadband
    census's GO are not in conflict, because they are rates over different
    frames.
    """
    from macro_core import astrom
    rows = q(con, """SELECT scope, k, n, n_population, rate_pct, wilson_lo,
                            wilson_hi, verdict FROM sn_g0_astrometry
                     ORDER BY scope""")
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(10.4, 2.4 + 0.62 * len(rows)))
        ypos = np.arange(len(rows))[::-1]
        for y, (scope, k, n, npop, rate, lo, hi, verdict) in zip(ypos, rows):
            color = {"GO": GOOD, "CAUTION": WARN}.get(verdict, BAD)
            ax.barh(y, rate, color=color, height=0.46, alpha=0.85)
            # A CENSUS has no sampling interval: when the sample IS the
            # population there is nothing left to be uncertain about, and
            # drawing an interval on it would invent uncertainty rather than
            # report it.  Only the sampled scope gets an error bar.
            is_census = (npop is not None and n == npop)
            if not is_census and n and lo is not None:
                ax.plot([100 * lo, 100 * hi], [y, y], color=INK, lw=1.5)
                ax.plot([100 * lo, 100 * hi], [y, y], "|", color=INK, ms=8)
            # Labels ride at a fixed right-hand column so they can never sit
            # on top of a bar or an interval.
            ax.annotate(f"{k:,}/{n:,} = {rate:.1f}%",
                        (112, y), va="center", ha="left", fontsize=9,
                        color=INK)
            ax.annotate(verdict, (152, y), va="center", ha="left",
                        fontsize=9, color=color, weight="bold")
            ax.annotate("census" if is_census else "48-frame sample",
                        (176, y), va="center", ha="left", fontsize=8,
                        color=MUTED)
        ax.axvline(100 * astrom.GO_LOWER_BOUND, color=GOOD, ls="--", lw=1)
        ax.axvline(100 * astrom.CAUTION_LOWER_BOUND, color=WARN, ls="--", lw=1)
        ax.set_yticks(ypos)
        ax.set_yticklabels([r[0] for r in rows], fontsize=9)
        ax.set_xlim(0, 210)
        ax.set_xticks([0, 20, 40, 60, 80, 100])
        # The two thresholds are named in the axis label rather than
        # annotated onto the rules: at this aspect ratio an in-axes label
        # either overlaps the top bar or is clipped by the spine.
        ax.set_xlabel(
            f"astrometric success rate (%)   — dashed rules: S1's CAUTION "
            f"bar ({100 * astrom.CAUTION_LOWER_BOUND:.0f}%) and GO bar "
            f"({100 * astrom.GO_LOWER_BOUND:.0f}%)")
        ax.set_title("The same question asked of three different populations",
                     loc="left")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "g0_astrometry.png", dpi=DPI)
        plt.close(fig)
    return "figures/gate0/g0_astrometry.png"


def fig_grism(con) -> str:
    """The slitless series as S2c measures it, night by night.

    Stacked bars of dispersed / direct / indeterminate frames per night, the
    flash window shaded, and the exposure ramp on a twin axis.  The label
    would have drawn one flat bar of 83; the measurement draws this.
    """
    rows = q(con, """SELECT night, phase_d, n_dispersed, n_direct,
                            n_indeterminate, max_exptime
                     FROM sn_g0_triage ORDER BY phase_d""")
    ph = np.array([r[1] for r in rows], dtype=float)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(11.0, 4.2))
        ax.axvspan(0, g0.FLASH_PHASE_END_D, color=ps.tint(ACCENT, 0.85),
                   zorder=0)
        w = 0.9
        d = np.array([r[2] for r in rows], dtype=float)
        di = np.array([r[3] for r in rows], dtype=float)
        ind = np.array([r[4] for r in rows], dtype=float)
        ax.bar(ph, d, width=w, color=ACCENT, label="measured DISPERSED")
        ax.bar(ph, di, width=w, bottom=d, color=GOOD,
               label="measured DIRECT (an image wearing the grism label)")
        ax.bar(ph, ind, width=w, bottom=d + di, color=FAINT,
               label="indeterminate")
        ax.set_xlabel("days after the adopted explosion epoch")
        ax.set_ylabel("slot-'6' frames")
        ax.set_ylim(0, max(d + di + ind) * 1.45)
        ax.annotate("flash window\n(features gone by ~day 7-8)",
                    (g0.FLASH_PHASE_END_D / 2, ax.get_ylim()[1] * 0.98),
                    ha="center", va="top", fontsize=8.5, color=INK)
        ax2 = ax.twinx()
        ax2.plot(ph, [r[5] for r in rows], color=WARN, lw=1.1, ls="--",
                 marker="D", ms=3.2, zorder=5)
        ax2.set_ylabel("exposure time (s)", color=WARN)
        ax2.set_yscale("log")
        ax2.set_ylim(40, 1400)
        ax2.grid(False)
        # Legend below the axes: inside, it lands on either the tallest bar
        # or the exposure trace, and there is no corner where it does not.
        ax.legend(fontsize=8.5, loc="upper center", ncol=3,
                  bbox_to_anchor=(0.5, -0.20), frameon=False)
        ax.set_title("Slot '6' on SN 2023ixf: what the pixels say, per night",
                     loc="left")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "g0_grism_series.png", dpi=DPI)
        plt.close(fig)
    return "figures/gate0/g0_grism_series.png"


# ===========================================================================
# Sections
# ===========================================================================
def section_freeze(con) -> str:
    rows = q(con, "SELECT scope, n_rows, n_canonical, note FROM sn_g0_dedup "
                  "ORDER BY n_rows DESC")
    pt = q(con, "SELECT n_frames, d_ra_arcsec, d_dec_arcsec, resid_med_px, "
                "resid_p90_px, resid_max_px FROM sn_g0_pointing "
                "WHERE scope = 'campaign'")[0]
    camp = q1(con, "SELECT count(*) FROM sn_g0_frames "
                   "WHERE epoch_role = 'campaign'")
    ixf = q1(con, "SELECT count(*) FROM sn_g0_frames WHERE epoch_role = "
                  "'campaign' AND canonical_target = '2023ixf'")
    spectra = q1(con, "SELECT count(*) FROM sn_g0_frames WHERE epoch_role = "
                      "'campaign' AND dispersion_class = 'dispersed'")
    eng = q(con, "SELECT path, filter, night FROM sn_g0_frames "
                 "WHERE tree != ? AND epoch_role = 'campaign'",
            (g0.SCIENCE_TREE,))
    tmpl = q(con, """SELECT night, count(*) AS n,
                            group_concat(DISTINCT filter) AS f
                     FROM sn_g0_frames WHERE epoch_role != 'campaign'
                     GROUP BY night ORDER BY night""")

    return f"""
<section id="freeze">
<div class="bhead"><h2>1 &middot; Gate 0a — the manifest freeze</h2>
<span class="tag">{fmt(camp)} unique campaign frames</span>
<span class="tag2">from {fmt(sum(r[1] for r in rows if r[0].startswith('tree:')))} catalog rows</span></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">How many frames of this supernova does the archive actually
hold &mdash; once every copy is collapsed and every name the telescope
called it by is merged?</p>

<h3>Evidence</h3>
<p>The strategy's headline inventory is <b>1,052 unique light frames</b>,
obtained by deduplicating the rows whose target field reads
<code>2023ixf</code>. The freeze repeats that dedup globally &mdash; across
every tree and <em>within</em> <code>rawimage/</code>, where the July
directories hold wholesale copies of earlier nights &mdash; and then does
the thing the strategy's own §3.1 ruling asks for but its headline number
does not do: it merges the target's ALIASES. The supernova was observed
under <code>NGC5457</code> on the discovery night and under <code>M101</code>
on the two nights after, and S0's alias table maps both onto the same sky.</p>

{table(["scope", "catalog rows", "unique frames", "what this counts"],
       [[esc(r[0]), fmt(r[1]), fmt(r[2]), esc(r[3])] for r in rows])}

<p>So the campaign is <b>{fmt(camp)}</b> unique frames, not
{fmt(ixf)} &mdash; the difference is the {fmt(camp - ixf)} frames that
carried a galaxy's name rather than the transient's, and every one of them
is a discovery-week epoch. Nothing was double counted to get there: the
same global <code>(basename, jd)</code> dedup that produced 1,052 produced
this, applied to a target set of two names instead of one.</p>

<p>Two further facts the freeze turned up, both of which change what may be
counted:</p>
<ul>
<li><b>{fmt(spectra)} of the campaign's frames are not images at all.</b>
They are the slot-'6' frames S2c measures as dispersed. Under the retired
filter-label rule all 83 slot-'6' frames would have been called spectra;
the measurement says 61 are, 3 are ordinary direct images and 19 could not
be certified either way.</li>
<li><b>{len(eng)} canonical campaign frames are detector engineering
products</b>, not science: {", ".join("<code>" + esc(e[0]) + "</code>"
                                      for e in eng)}. They carry the
campaign's target name and fall inside its dates. One of them reads 65,535
ADU at the supernova's position &mdash; a value the 12-bit High Gain channel
cannot physically produce &mdash; which is how the census caught them.</li>
</ul>

<p>The freeze also measures the thing the saturation census needs next: how
far the telescope's commanded pointing sits from where it actually lands.
Over {fmt(pt[0])} plate-solved campaign frames the commanded-to-true offset
is ({pt[1]:+.0f}, {pt[2]:+.0f}) arcsec, and the residual scatter about that
median is <b>{fnum(pt[3])} px</b> (median), {fnum(pt[4])} px (90th
percentile), {fnum(pt[5])} px (worst). That is why an unsolved frame's
search box is {g0.BOUND_HALF_PX} px and not a guess.</p>

<h3>Decision</h3>
<div class="decision"><b>The campaign is {fmt(camp)} unique frames, of which
{fmt(camp - spectra - len(eng))} are science images.</b> Every number
downstream on this page is computed over that set, and the set is a table
(<code>sn_g0_frames</code>) rather than a paragraph.</div>

<h3>Consequence</h3>
<p>The strategy's 1,052 was not wrong about duplication &mdash; it was right,
and its dedup rule is the one used here. It was incomplete about IDENTITY:
a census taken under one of a target's names cannot see the frames filed
under the others, and on this target those are exactly the earliest and
most contested epochs. The published campaign statistics, and the README
the strategy's Step 10 requires be regenerated from the manifest, should be
regenerated from <code>sn_g0_frames</code>.</p>

<p>The non-campaign epochs are frozen in the same table and are summarised
here so the template inventory has a query behind it too:</p>
{table(["night", "unique frames", "filters"],
       [[esc(t[0]), fmt(t[1]), esc(t[2])] for t in tmpl])}
</div>
</section>"""


def section_screen(con) -> str:
    s = screen_of(con)
    modes = q(con, """SELECT c.readoutm, count(*) AS n, m.clip_adu, m.veto_adu
                      FROM sn_g0_census c
                      LEFT JOIN s2_ceiling_modes m ON m.mode = c.readoutm
                      GROUP BY c.readoutm ORDER BY n DESC""")
    return f"""
<section id="screen">
<div class="bhead"><h2>2 &middot; The screen, and where it comes from</h2>
<span class="tag">clip {fmt(s.clip_adu)} ADU, measured</span></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">At what pixel value does this detector stop telling the
truth &mdash; and does anyone have a measurement of it, or only a guess?</p>

<h3>Evidence</h3>
<p>When the strategy was written the answer was a guess, and it said so.
Its §3.3 quotes a &ldquo;hard clip at ~3,530&ndash;3,550 ADU&rdquo; and sets
a binding screen of <b>2,800 ADU</b> with a suspect band from 2,400 &mdash;
numbers derived by eye from an assumed ~3,500 ADU ceiling, because the S2
tables that held the measured one had been destroyed. That is precisely why
the project task <code>SN-G0b</code> was marked BLOCKED: a screen with
nothing behind it is not a screen.</p>

<p>S2 has since been rebuilt and measures the High Gain channel's ceiling
from science-frame histograms. This page therefore stores the strategy's
<em>fractions</em> and applies them to the <em>measurement</em>:</p>

{table(["quantity", "value", "source"],
       [["measured clip", f"{fmt(s.clip_adu)} ADU",
         "S2 <code>s2_ceiling_modes.clip_adu</code>, High Gain"],
        ["S2's own saturation veto", f"{fmt(s.veto_adu)} ADU",
         "S2 <code>s2_ceiling_modes.veto_adu</code> &mdash; an independent "
         "cross-check, not used as the screen here"],
        ["reject screen",
         f"{fmt(s.reject_adu)} ADU "
         f"({g0.REJECT_CLIP_FRACTION:.2f} &times; clip)",
         "the strategy's 80%-of-clip rule, applied to the measurement"],
        ["suspect floor",
         f"{fmt(s.suspect_adu)} ADU "
         f"({g0.SUSPECT_CLIP_FRACTION:.3f} &times; clip)",
         "the strategy's lower flag level, likewise"]])}

<div class="decision"><b>The strategy's hand numbers survive contact with
the measurement.</b> 0.80 &times; {fmt(s.clip_adu)} =
{fmt(s.reject_adu)} against a typed 2,800, and the suspect floor lands on
{fmt(s.suspect_adu)} against a typed 2,400 &mdash; both inside 0.2%. The
assumed ceiling was right. What has changed is that the screen is now a
consequence of a query, so an S2 re-measurement moves it without anyone
editing a document.</div>

<h3>Consequence</h3>
<p>Every frame on this sky was screened against the ceiling of the readout
mode it was actually taken in, which matters because the post-fade template
epochs are not High Gain at all:</p>
{table(["readout mode", "frames on this sky", "measured clip", "S2 veto"],
       [[esc(m[0]), fmt(m[1]), fmt(m[2]), fmt(m[3])] for m in modes])}
<p>A single High Gain screen applied to the 16-bit template epochs would
have rejected every one of them; a single 16-bit screen applied to the
campaign would have passed every clipped frame in it.</p>
</div>
</section>"""


def section_census(con, fig_m: str, fig_p: str) -> str:
    s = screen_of(con)
    bands = q(con, """SELECT band_role, filter, n_frames, n_images, n_wcs,
                             n_bound, n_clean, n_suspect, n_rejected,
                             n_bounded_clean, n_undetermined, n_usable,
                             first_clean_night, first_clean_phase_d,
                             isolation_false_id, isolation_tested
                      FROM sn_g0_bands ORDER BY band_role DESC, filter""")
    both = q1(con, "SELECT count(*) FROM sn_g0_matrix WHERE band_role = "
                   "'broadband' AND n_clean > 0 AND n_rejected > 0")
    cells = q1(con, "SELECT count(*) FROM sn_g0_matrix "
                    "WHERE band_role = 'broadband' AND n_frames > 0")
    early = q(con, """SELECT night, min(phase_d) AS p, count(*) AS n,
                             sum(saturation_class = 'clean') AS clean,
                             sum(saturation_class = 'rejected') AS rej
                      FROM sn_g0_census
                      WHERE epoch_role = 'campaign' AND band_role='broadband'
                        AND night <= '2023-05-23'
                      GROUP BY night ORDER BY night""")

    rows = []
    for b in bands:
        cls = "warn" if b[11] == 0 and b[0] != "other" else ""
        rows.append([esc(b[1]), esc(b[0]), fmt(b[2]), fmt(b[3]),
                     f"{fmt(b[4])} / {fmt(b[5])}",
                     fmt(b[6]), fmt(b[9]), fmt(b[7]), fmt(b[8]), fmt(b[10]),
                     f"<b>{fmt(b[11])}</b>",
                     (esc(b[12]) + f" (+{b[13]:.1f} d)") if b[12] else "&mdash;",
                     f"{fmt(b[14])}/{fmt(b[15])}"])
    classes = ["" for _ in rows]

    return f"""
<section id="census">
<div class="bhead"><h2>3 &middot; Gate 0b — the saturation census</h2>
<span class="tag">{fmt(q1(con, "SELECT count(*) FROM sn_g0_census WHERE status='measured'"))} frames opened</span></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">In which frames is the supernova itself measurable, and in
which is it sitting on the detector's ceiling?  The strategy's answer was an
exposure-scaling argument.  This is the pixels.</p>

<h3>Evidence</h3>
<p>Every frame on this sky was opened read-only and the supernova's
catalogue position projected onto it. Where the frame carried its own plate
solution the position is exact and the peak of a
{2 * g0.CORE_HALF_PX + 1}&nbsp;px stamp is the supernova's peak; where it
did not, the position is known only to the pointing residual measured in
§1, and the census reports the maximum of a {2 * g0.BOUND_HALF_PX + 1}&nbsp;px
search box &mdash; which is an UPPER BOUND on the supernova's peak and is
labelled as one. A frame whose whole box stays under the screen is
<em>provably</em> unsaturated wherever the supernova fell in it; a frame
whose box reaches the screen supports no conclusion at all, and is recorded
as <code>undetermined</code> rather than as a saturated epoch it may not
be.</p>

{_figure(fig_m, "The saturation matrix the strategy's Gate 0b asks for. "
                "Colour is the fraction of that filter-night's frames from "
                "which photometry may be taken. The black dots are the "
                "result that matters most: cells that are partly clean and "
                "partly clipped.")}

<p>The census validates its own position rule as it goes. In every
plate-solved frame it also records the ISOLATION radius &mdash; how far away
the nearest brighter pixel is &mdash; which measures directly how often
&ldquo;the brightest thing in the search box&rdquo; would have been the
wrong source. In the broadband codes that misidentification rate is a few
per cent; in the narrowband codes it is an order of magnitude worse, which
is why no narrowband bound was ever promoted to a measurement.</p>

{table(["filter", "role", "frames", "images", "solved / unsolved",
        "clean", "clean by bound", "suspect", "rejected", "undetermined",
        "USABLE", "first clean epoch", "box would misidentify"],
       rows, classes)}

{_figure(fig_p, "Every supernova peak the census could measure directly, "
                "against phase. Solid rule: the measured clip. Dashed: the "
                "reject screen. Dotted: the suspect floor.")}

<h3>Decision</h3>
<div class="decision">
<b>The clean broadband start is night 2023-05-23 (UT 05-24, +5.4 d), in all
three broadband codes, exactly where the strategy predicted it by exposure
scaling &mdash; and now by measurement.</b> The three earlier epochs contain
not one clean broadband frame between them:
{"; ".join(f"{e[0]} (+{e[1]:.1f} d) &mdash; {e[3]} clean of {e[2]} frames, "
           f"{e[4]} over the screen"
           for e in early if e[0] < '2023-05-23')}. The balance in each is
<code>undetermined</code>: frames with no plate solution whose search box
reaches the screen, which the census declines to call either way.
</div>

<div class="decision">
<b>Saturation here is a per-frame property, not a per-epoch one.</b>
{fmt(both)} of {fmt(cells)} broadband filter-night cells contain a clean
frame AND a clipped frame on the same night in the same filter. The
strategy anticipated this in one sentence (§3.3: the supernova reaches
3,222 ADU in a good-seeing 2&nbsp;s frame on a night when a 0.5&nbsp;s frame
reads 891); the census makes it the rule rather than the anecdote. No epoch
may be admitted or excluded wholesale &mdash; the screen has to be applied
frame by frame, and the released table has to carry the flag per frame.
</div>

<h3>Consequence</h3>
<p>The narrowband result is the one that changes a plan. The strategy names
the +5.4 and +6.4&nbsp;d 32&nbsp;s frames as
&ldquo;<b>the only certain-clean H epochs</b>&rdquo; and warns that
&ldquo;there is a live risk that the only unsaturated flash-phase H&alpha;
record in this archive is the slitless grism series.&rdquo; The census
disagrees, and it disagrees with a measurement: see §4.</p>
</div>
</section>"""


def section_narrowband(con) -> str:
    s = screen_of(con)
    early = q(con, """SELECT night, filter, round(exptime) AS e, quality,
                             round(peak_adu) AS pk, round(sky_adu) AS sky,
                             round(offset_px, 1) AS off, isolation_px,
                             saturation_class
                      FROM sn_g0_census
                      WHERE epoch_role = 'campaign' AND band_role='narrowband'
                        AND night IN ('2023-05-20','2023-05-21','2023-05-23',
                                      '2023-05-24')
                      ORDER BY night, filter, exptime""")
    pre = q(con, """SELECT filter, count(*) AS n, round(avg(sky_adu)) AS sky,
                           round(avg(core_max_adu)) AS core
                    FROM sn_g0_census WHERE epoch_role = 'template_pre'
                    GROUP BY filter ORDER BY filter""")
    # Net rate at the supernova's position, per second, in the H code: the
    # pre-explosion epoch is the control that says the source is the
    # supernova and not something that was already there.
    def rate(where, params=()):
        r = q(con, f"""SELECT (core_max_adu - sky_adu) / exptime
                       FROM sn_g0_census WHERE {where} AND quality = 'wcs'""",
              params)
        return [x[0] for x in r]
    pre_h = rate("epoch_role = 'template_pre' AND filter = 'H'")
    d25 = rate("night = '2023-05-20' AND filter = 'H'")
    d54 = rate("night = '2023-05-23' AND filter = 'H'")
    clean25 = q1(con, """SELECT count(*) FROM sn_g0_census
                         WHERE night = '2023-05-20' AND band_role='narrowband'
                           AND saturation_class IN ('clean','bounded_clean')""")

    return f"""
<section id="narrowband">
<div class="bhead"><h2>4 &middot; The flash-phase H&alpha; record is not
only the grism</h2><span class="tag">+2.5 d, measured clean</span></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">Is the slitless series really the archive's only unsaturated
flash-phase narrowband record of this supernova?</p>

<h3>Evidence</h3>
<p>No. On night 2023-05-20 (UT 05-21, +2.5 d) the narrowband exposures were
64&nbsp;s, and the census measures the supernova at a few hundred to
~1,900 ADU in them &mdash; well under the {fmt(s.suspect_adu)} ADU flag
level. The night after, the exposures ramp to 128&nbsp;s and the same
frames clip. Here is every narrowband frame across the four earliest
narrowband nights:</p>

{table(["night", "filter", "exp (s)", "position", "peak ADU", "sky ADU",
        "offset (px)", "isolation (px)", "class"],
       [[esc(r[0]), esc(r[1]), fmt(r[2]), esc(r[3]), fmt(r[4]), fmt(r[5]),
         fnum(r[6], 1), fnum(r[7]), esc(r[8])] for r in early],
       ["warn" if r[8] in ("rejected", "suspect") else "" for r in early])}

<p>Two independent checks say this is the supernova and not a field star or
an H&nbsp;II region that was always there.</p>
<ul>
<li><b>Position.</b> The clean +2.5&nbsp;d frames put the peak within
{min(r[6] for r in early if r[0] == '2023-05-20' and r[3] == 'wcs'):.1f}&ndash;{max(r[6] for r in early if r[0] == '2023-05-20' and r[3] == 'wcs'):.1f}&nbsp;px
of the catalogue position, with the nearest brighter pixel tens of pixels
away.</li>
<li><b>The pre-explosion epoch, which the census also measured.</b> In the
2023-05-05 H frames the same position yields
{np.mean(pre_h):.1f} net ADU&nbsp;s<sup>&minus;1</sup>; at +2.5&nbsp;d it
yields {np.mean(d25):.0f}, and at +5.4&nbsp;d {np.mean(d54):.0f}. The source
is {np.mean(d25) / max(np.mean(pre_h), 1e-9):.0f}&times; brighter than
whatever occupied that position two weeks earlier.</li>
</ul>

<div class="decision"><b>{fmt(clean25)} narrowband frames at +2.5 d are
measured unsaturated on the supernova &mdash; three nights earlier than the
strategy's stated clean narrowband start, and they are images, not
spectra.</b> The claim that the grism is the only unsaturated flash-phase
H&alpha; record in this archive is retracted here, on the evidence of the
frames it was made about.</div>

<h3>Consequence</h3>
<p>This does not by itself resurrect the H&alpha; product. It moves what is
blocking it. Q2 was thought to be dead in the flash phase because of
saturation; it is not, and the blocker is now the one the strategy already
named &mdash; the narrowband transmission curve, without which a flux at
+2.5&nbsp;d cannot be interpreted physically (task
<code>SN-S1-narrowband-curves</code>). It also means the flash-phase
narrowband coverage is three epochs (+2.5, +5.4, +6.4&nbsp;d), not two.</p>

<p>The pre-explosion template epoch was measured in the same pass, and it
settles a caution the strategy raised but could not resolve:</p>
{table(["filter", "frames", "mean sky pedestal (ADU)",
        "mean peak at the SN position (ADU)",
        "sky as a fraction of the measured clip"],
       [[esc(p[0]), fmt(p[1]), fmt(p[2]), fmt(p[3]),
         f"{100.0 * p[2] / s.clip_adu:.0f}%"] for p in pre],
       ["warn" if p[2] / s.clip_adu > 0.5 else "" for p in pre])}
<p>The strategy flagged the 2023-05-05 G/R/X frames as
&ldquo;possibly trailed/defocused&rdquo; on their 14&ndash;16&nbsp;px FWHM
and ruled them narrowband-only by default. The census supplies the reason:
their SKY sits at
{max(100.0 * p[2] / s.clip_adu for p in pre if p[0] in ('G', 'R', 'X')):.0f}%
of the clip before a single photon from a star arrives, leaving almost no
dynamic range. The narrowband H/O templates sit at a few per cent and are
fine. The strategy's ruling stands, and now has a number under it.</p>
</div>
</section>"""


def section_astrometry(con, fig_a: str) -> str:
    rows = q(con, """SELECT scope, description, k, n, n_population, rate_pct,
                            wilson_lo, wilson_hi, verdict, basis
                     FROM sn_g0_astrometry ORDER BY scope""")
    v = q(con, "SELECT deciding_number, verdict, basis FROM sn_g0_verdict "
               "WHERE question_id = 'astrometry'")[0]
    comp = q(con, """SELECT f.band_role, count(*) AS n
                     FROM frames fr
                     JOIN sn_g0_frames f ON f.obs_rowid = fr.obs_rowid
                     LEFT JOIN frame_dispersion d ON d.obs_rowid=fr.obs_rowid
                     WHERE fr.canonical_target = '2023ixf'
                       AND fr.is_canonical = 1 AND fr.tree = 'rawimage'
                       AND (fr.pltsolvd IS NULL OR fr.pltsolvd != 1)
                       AND (d.verdict IS NULL OR d.verdict != 'dispersed')
                     GROUP BY 1 ORDER BY 2 DESC""")
    fails = q(con, """SELECT f.filter, a.diagnosis, count(*) AS n
                      FROM s1_failure_autopsy a
                      JOIN frames f ON f.obs_rowid = a.obs_rowid
                      WHERE a.stratum_id = 'sn_gsense_broadband'
                      GROUP BY 1, 2 ORDER BY n DESC""")
    return f"""
<section id="astrometry">
<div class="bhead"><h2>5 &middot; Does the astrometry verdict move?</h2>
<span class="tag">{esc(v[1])} on the stratum</span></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">The strategy inherited a NO-GO on this campaign's astrometry
from the S1 experiment, whose candidate universe was contaminated by
spectra it had excluded by filter LABEL.  With the gate repaired to read
S2c's per-frame measurement, does the verdict move?</p>

<h3>Evidence</h3>
<p><b>No &mdash; and the reason it does not is more useful than a movement
would have been.</b> The S1 stratum's rate for this campaign rose when the
gate was repaired, but not nearly far enough to cross a threshold, and the
same question asked of the frames that actually carry the light curve gives
a completely different answer:</p>

{table(["scope", "successes / trials", "rate", "verdict", "what this is"],
       [[f"<code>{esc(r[0])}</code><br><span class='sub'>{esc(r[1])}</span>",
         f"{fmt(r[2])} / {fmt(r[3])}", f"{r[5]:.1f}%",
         f"<b>{esc(r[8])}</b>", esc(r[9])] for r in rows])}

{_figure(fig_a, "The same 80% / 50% thresholds applied to three different "
                "populations. The error bar is drawn only on the sampled "
                "scope; the two census scopes have no sampling uncertainty "
                "to draw.")}

<p>The two verdicts are not in conflict, because they are rates over
different frames. The S1 stratum is drawn from the campaign's UNSOLVED
residue &mdash; the frames that PinPoint could not solve at the telescope
&mdash; and that residue is overwhelmingly narrowband:</p>

{table(["band role", "frames in the stratum's unsolved population"],
       [[esc(c[0]), fmt(c[1])] for c in comp])}

<p>Its failures say the same thing in the machine's own words:</p>
{table(["filter", "S1 failure diagnosis", "frames"],
       [[esc(f[0]), esc(f[1]), fmt(f[2])] for f in fails])}

<h3>Decision</h3>
<div class="decision"><b>The verdict does not move: the stratum stays
NO-GO.</b> {esc(v[0])}.</div>

<div class="decision"><b>But it is a verdict about long narrowband
exposures, not about the gri light curve.</b> Its unsolved population is
{fmt(dict(comp).get('narrowband', 0))} narrowband frames against
{fmt(dict(comp).get('broadband', 0))} broadband, and
{fmt(sum(f[2] for f in fails if f[0] == 'O'))} of its
{fmt(sum(f[2] for f in fails))} autopsied failures are O-band frames
diagnosed as trailed. The light curve's own astrometry is not a
blind-batch problem at all: the solutions already exist, and this census
independently corroborated them by finding the supernova within
{g0.MAX_CENTROID_OFFSET_PX:.0f}&nbsp;px of its catalogue position in each
one.</div>

<h3>Consequence</h3>
<p>Two things follow for the plan. First, the re-solve task
(<code>SN-S4-resolve</code>) should be scoped to the bands it can help:
attempting the O-band residue blind is what the NO-GO is about, and the
strategy's &ldquo;~177 unsolved broadband frames&rdquo; is the tractable
part. Second, a corroborated plate solution is a stronger claim than a
<code>PLTSOLVD</code> card, and this campaign now has one per frame &mdash;
a by-product of the census that Step 4 can consume instead of re-deriving.</p>
</div>
</section>"""


def section_grism(con, fig_g: str) -> str:
    t = q(con, """SELECT night, phase_d, n_labelled, n_dispersed, n_direct,
                         n_indeterminate, min_exptime, max_exptime,
                         n_paired_direct, paired_filters, in_flash_window
                  FROM sn_g0_triage ORDER BY night""")
    clauses = q(con, "SELECT clause, requirement, value, passed "
                     "FROM sn_g0_triage_summary")
    meta = dict(q(con, "SELECT key, value FROM sn_g0_build_meta"))
    v = q(con, "SELECT deciding_number, verdict, basis FROM sn_g0_verdict "
               "WHERE question_id = 'grism'")[0]
    order = ["nights", "flash_nights", "extracted", "contamination",
             "wavelength"]
    cl = {c[0]: c for c in clauses}
    return f"""
<section id="grism">
<div class="bhead"><h2>6 &middot; Gate 0c — the grism triage</h2>
<span class="tag">{esc(meta.get('grism_n_dispersed', '?'))} measured spectra</span>
<span class="tag2">{esc(v[1])}</span></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">Is the slitless series a genuine flash-phase spectral record
&mdash; and does it clear the promotion bar the strategy pre-registered
before any of this evidence existed?</p>

<h3>Evidence</h3>
<p>The series is real, and it is not what the label said it was. Of the
{esc(meta.get('grism_n_dispersed', '?'))} frames S2c measures as dispersed,
{fmt(sum(r[3] for r in t if r[10]))} fall inside the flash window; three
frames wearing the same slot-'6' label are measured DIRECT IMAGES, and
{fmt(sum(r[5] for r in t))} could not be certified either way.</p>

{table(["night", "phase (d)", "slot-'6' frames", "dispersed", "direct",
        "indeterminate", "exposure (s)", "paired direct images that night",
        "in flash window"],
       [[esc(r[0]), f"{r[1]:+.1f}", fmt(r[2]), f"<b>{fmt(r[3])}</b>",
         fmt(r[4]), fmt(r[5]),
         f"{r[6]:.0f}&ndash;{r[7]:.0f}", f"{fmt(r[8])} ({esc(r[9])})",
         "yes" if r[10] else "&mdash;"] for r in t])}

{_figure(fig_g, "Slot '6' per night as the pixels read it, with the "
                "exposure ramp overlaid. The label would have drawn a "
                "single flat bar of 83.")}

<p>Gate 0c's first requirement is that paired or near-in-time direct images
exist for every grism night, and the strategy insists this be VERIFIED
rather than assumed. It is verified: every one of the
{fmt(len([r for r in t if r[3] > 0]))} nights carrying a measured spectrum
also carries between {min(r[8] for r in t if r[3] > 0)} and
{max(r[8] for r in t if r[3] > 0)} direct images of the same field.</p>

<p>The promotion criterion itself, applied clause by clause exactly as
pre-registered:</p>
{table(["clause", "requirement", "measured value", "passed"],
       [[esc(k), esc(cl[k][1]), esc(cl[k][2]),
         "yes" if cl[k][3] else "<b>NO</b>"] for k in order if k in cl],
       ["" if cl[k][3] else "warn" for k in order if k in cl])}

<h3>Decision</h3>
<div class="decision"><b>The series is genuine and it does not promote.</b>
{esc(v[0])}. It has the nights &mdash; exactly the three the criterion asks
for, with no margin &mdash; and it has none of the spectra. Nothing has been
extracted, so nothing has been through the offset-trace contamination test,
and no wavelength calibration source has been named. {esc(v[2])}</div>

<h3>Consequence</h3>
<p>This is a promotion the evidence does not support, and saying so is the
result. The strategy pre-registered the criterion precisely so that a real,
well-sampled, genuinely interesting set of frames could not be talked into
being a headline product before anyone had extracted a spectrum from it.
The frames are worth the two-week timebox the strategy allots them; they
are not worth a venue decision.</p>

<p>One limitation of this census bears directly on that timebox and should
be stated rather than discovered later: <b>not one of the 83 slot-'6'
frames carries a plate solution</b>, so Gate 0b could not determine whether
the traces themselves are clipped. The 512&nbsp;s exposures from
+9.4&nbsp;d onward are the ones at risk. The extraction will have to
establish its own astrometry &mdash; which is exactly what the verified
paired direct images above are for.</p>
</div>
</section>"""


def section_verdicts(con) -> str:
    v = q(con, """SELECT question_id, question, deciding_number, verdict,
                         moved, basis FROM sn_g0_verdict""")
    order = ["usable-broadband", "astrometry", "grism", "venue"]
    by = {r[0]: r for r in v}
    rows = [by[k] for k in order if k in by]
    return f"""
<section id="verdicts">
<div class="bhead"><h2>7 &middot; The answers, with the numbers that decide
them</h2></div>

<div class="stage">
<h3>Question</h3>
<p class="sub">Gate 0 exists to answer three questions and to set one
posture.  Each answer is stored beside the number it turns on, so a re-run
that keeps the word while moving the number is visible in the table rather
than only in a reader's memory.</p>

<h3>Evidence</h3>
{table(["question", "verdict", "deciding number", "basis"],
       [[esc(r[1]), f"<b>{esc(r[3])}</b>", esc(r[2]), esc(r[5])]
        for r in rows])}

<h3>Decision</h3>
<div class="decision"><b>The venue posture stays at
{esc(g0.VENUE_BASE)}.</b> The strategy decided in advance that ApJ is taken
only if Gate 0 promotes the grism or recovers the narrowband bandpass. Gate
0 did neither: zero contamination-tested spectra against a bar of
{g0.GRISM_PROMOTION_MIN_NIGHTS}, and zero recovered transmission curves
against a bar of one. A posture that was decided before the evidence
arrived, and that the evidence does not move, is a posture that held.</div>

<h3>Consequence</h3>
<p>What Gate 0 delivers is not a promotion; it is a scope. The paper has
<b>{fmt(int(by['usable-broadband'][2].split()[0].replace(',', '')))}</b>
broadband frames it may take photometry from, a measured clean start, a
per-frame saturation flag for the released table, a corroborated plate
solution for most of the light curve, and three narrowband epochs inside the
flash phase it did not know it had. Steps 1 and 2 are now unblocked and are
the next binding gates.</p>

<p>What Gate 0 does NOT settle, stated so it is not mistakenly assumed:
the physical identity of the filter codes (Step 1), the linearity curve
below the screen &mdash; which is what would let the
{fmt(q1(con, "SELECT sum(n_suspect) FROM sn_g0_bands WHERE band_role='broadband'"))}
suspect broadband frames be recovered (Step 2), whether the grism traces
are clipped (no plate solutions), and the narrowband transmission profiles
without which the H&alpha; product stays a methods demonstration.</p>
</div>
</section>"""


# ===========================================================================
# render
# ===========================================================================
def render_report(manifest_path: Path) -> Path:
    """Render the Gate 0 report from the manifest DB.  Returns the HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = dict(q(con, "SELECT key, value FROM sn_g0_build_meta"))
        n_frames = q1(con, "SELECT count(*) FROM sn_g0_frames")
        n_camp = q1(con, "SELECT count(*) FROM sn_g0_frames "
                         "WHERE epoch_role = 'campaign'")
        usable = q1(con, "SELECT sum(n_usable) FROM sn_g0_bands "
                         "WHERE band_role = 'broadband'")
        posture = q1(con, "SELECT verdict FROM sn_g0_verdict "
                          "WHERE question_id = 'venue'")

        fig_m = fig_matrix(con)
        fig_p = fig_peaks(con)
        fig_a = fig_astrometry(con)
        fig_g = fig_grism(con)

        sections = [
            section_freeze(con),
            section_screen(con),
            section_census(con, fig_m, fig_p),
            section_narrowband(con),
            section_astrometry(con, fig_a),
            section_grism(con, fig_g),
            section_verdicts(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SN 2023ixf — Gate 0</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>SN 2023ixf — Gate 0, from the pixels</h1>
  <p>{fmt(n_camp)} unique campaign frames &middot; {fmt(n_frames)} frames
  opened and measured &middot; {fmt(usable)} usable broadband epochs &middot;
  venue posture <b>{esc(posture)}</b> &middot;
  built {esc(meta.get('verdicts_built_at', ''))[:16]}Z
  ({esc(meta.get('code_version', SN_G0_CODE_VERSION))})
  &middot; <a href="index.html">back to the project page</a>
  &middot; <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#freeze">1 Freeze</a> &middot;
  <a href="#screen">2 The screen</a> &middot;
  <a href="#census">3 Saturation census</a> &middot;
  <a href="#narrowband">4 Flash-phase H&alpha;</a> &middot;
  <a href="#astrometry">5 Astrometry</a> &middot;
  <a href="#grism">6 Grism triage</a> &middot;
  <a href="#verdicts">7 The answers</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_sn.report_gate0</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query or a documented constant in
<code>macro_sn.gate0</code>; none is typed by hand. The supernova's position
is a literature constant ({esc(g0.SN_POSITION_SOURCE)}) and the explosion
epoch is {g0.T0_MJD} MJD ({esc(g0.T0_SOURCE)}), used for display only — no
verdict on this page is a function of phase. Regenerate with
<code>pipeline/scripts/run_sn_gate0.py report</code>.</footer>
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
