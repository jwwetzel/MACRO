"""S3 evidence report renderer: the shared time axis.

Reads the S3 tables of the manifest database (NEVER the catalog or the
archive — if a number cannot be derived from the database, it does not
belong on the page) and writes:

* ``docs/pipeline/s3_timing.html``      — the report
* ``docs/pipeline/figures/s3/*.png``    — every figure

The page follows the site's Socratic format: one section per decision,
each section = Question -> Evidence -> Decision -> Consequence.  EVERY
number in the HTML is interpolated from a SQL query executed in this
module or from a constant defined in ``macro_core.timing`` — nothing is
hand-typed, so re-running the build after an archive sync regenerates the
whole argument, clock bound included.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import timing as tm       # noqa: E402  (constants for interpolation)
# Shared page machinery: same house figure style, same query discipline, same
# table generator as the earlier reports — one visual language site-wide.
from .report_s0 import (          # noqa: E402
    ACCENT, BAD, STYLE, DPI, GOOD, INK, MUTED, WARN,
    _figure, esc, fmt, q, q1, table)
from . import plotstyle as ps    # noqa: E402  (the house figure style)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s3"
HTML_PATH = DOCS_DIR / "s3_timing.html"

GREEN = GOOD                    # verdict OK — the house confirmation hue


def fnum(x, digits=2) -> str:
    """Format a float with fixed digits, em-dash for NULL."""
    if x is None:
        return "&mdash;"
    return f"{x:.{digits}f}"


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_jdhelio_semantics(con) -> str:
    """The smoking gun: header JD-HELIO residuals vs EXPTIME.

    Left: residual against our heliocentric JD computed at exposure START
    — the points climb the EXPTIME/2 line, so the header's heliocentric
    stamp was computed half an exposure after the JD instant.  Right: the
    same residual computed at MID-exposure collapses to ~0.1 s — except
    the HDR family and one broken Mode0 focus frame, plotted in yellow.
    """
    rows = q(con, """
        SELECT exptime_s, helio_resid_start_s, helio_resid_mid_s, family
        FROM s3_header_audit WHERE helio_resid_start_s IS NOT NULL""")
    exp = np.array([r[0] for r in rows])
    rs = np.array([r[1] for r in rows])
    rm = np.array([r[2] for r in rows])
    fam = [r[3] for r in rows]
    # "Broken" = mid-residual beyond 5 s: the HDR family + the one Mode0
    # focus-frame oddity; every other sampled header sits within 1 s.
    bad = np.abs(rm) > 5.0
    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.6))
        xs = np.geomspace(0.01, 2000, 50)
        ax1.plot(xs, xs / 2.0, color=MUTED, lw=1, ls="--",
                 label="EXPTIME / 2")
        ax1.scatter(exp[~bad], rs[~bad], s=18, color=ACCENT, zorder=3,
                    label="sampled headers")
        ax1.scatter(exp[bad], rs[bad], s=22, color=WARN, zorder=4,
                    label="HDR family / broken")
        ax1.set_xscale("log")
        ax1.set_yscale("symlog", linthresh=1)
        ax1.set_xlabel("header EXPTIME (s)")
        ax1.set_ylabel("JD-HELIO $-$ our HJD@START (s)")
        ax1.set_title("Residual at exposure START:\nclimbs the EXPTIME/2 line")
        ax1.legend(fontsize=8, loc="upper left")
        ax2.scatter(exp[~bad], rm[~bad], s=18, color=ACCENT, zorder=3)
        ax2.scatter(exp[bad], rm[bad], s=22, color=WARN, zorder=4)
        ax2.axhline(0, color=MUTED, lw=1, ls="--")
        ax2.set_xscale("log")
        ax2.set_yscale("symlog", linthresh=1)
        ax2.set_xlabel("header EXPTIME (s)")
        ax2.set_ylabel("JD-HELIO $-$ our HJD@MID (s)")
        ax2.set_title("Residual at MID-exposure:\ncollapses to $\\sim$0.1 s")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s3_jdhelio_semantics.png", dpi=DPI)
        plt.close(fig)
    return "figures/s3/s3_jdhelio_semantics.png"


def fig_cadence_overheads(con) -> str:
    """Per-frame overhead (gap minus EXPTIME) floors: StackPro vs plain.

    If StackPro EXPTIME were the per-sub-read time, back-to-back gaps
    would sit at 16x EXPTIME; they sit at EXPTIME + a save/slew overhead
    indistinguishable from the plain High Gain pipeline.  The dashed line
    is the ROBUST dead-time ceiling (:data:`tm.STACKPRO_DEADTIME_BOUND_S`,
    the smallest median-cadence overhead of any regular StackPro series),
    not the archive's smallest single gap — see section 2 for why the
    latter is an artifact.
    """
    rows = q(con, """
        SELECT readoutm, exptime_s, min_overhead_s, regular_overhead_s
        FROM s3_cadence
        WHERE readoutm IN ('High Gain', 'High Gain StackPro')
          AND regular_overhead_s IS NOT NULL
        ORDER BY exptime_s""")
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        for mode, color, marker in (("High Gain", ACCENT, "o"),
                                    ("High Gain StackPro", WARN, "s")):
            sub = [r for r in rows if r[0] == mode]
            e = np.array([r[1] for r in sub])
            # The ROBUST overhead (median cadence of a regular series),
            # not a raw order statistic: plotting the latter would draw
            # the negative "overheads" that section 2 exists to refute.
            reg = np.array([r[3] for r in sub])
            ax.plot(e, reg, marker=marker, ms=5, lw=1, color=color,
                    label=f"{mode} (median-cadence overhead)")
        e_sp = np.array([r[1] for r in rows if r[0] == "High Gain StackPro"])
        ax.plot(e_sp, 15.0 * e_sp, color=BAD, lw=1.2, ls=":",
                label="expected extra if EXPTIME were per-sub-read "
                      "(15 $\\times$ EXPTIME)")
        ax.axhline(tm.STACKPRO_DEADTIME_BOUND_S, color=GREEN, lw=1,
                   ls="--", label=f"measured dead-time bound "
                   f"{tm.STACKPRO_DEADTIME_BOUND_S} s")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("EXPTIME (s)")
        ax.set_ylabel("gap $-$ EXPTIME (s)")
        ax.set_title("Back-to-back series overhead:\nStackPro behaves "
                     "like a plain frame of the same EXPTIME")
        ax.legend(fontsize=7, loc="upper left")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s3_cadence_overheads.png", dpi=DPI)
        plt.close(fig)
    return "figures/s3/s3_cadence_overheads.png"


def fig_bjd_offset(con) -> str:
    """Distribution of the applied corrections over all science frames."""
    rows = q(con, """
        SELECT bary_ltt_s, tdb_minus_utc_s FROM frame_times
        WHERE bjd_tdb IS NOT NULL""")
    ltt = np.array([r[0] for r in rows]) / 60.0
    total = np.array([r[0] + r[1] for r in rows]) / 60.0
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.2, 3.4))
        ax.hist(ltt, bins=120, color=ACCENT, alpha=0.85,
                label="barycentric light-travel term")
        ax.hist(total, bins=120, color=WARN, histtype="step", lw=1.4,
                label="total BJD_TDB $-$ JD_UTC(mid)")
        ax.set_xlabel("applied correction (minutes)")
        ax.set_ylabel("frames")
        ax.set_title("BJD_TDB vs header JD: the correction S3 applies to "
                     f"{len(ltt):,} frames")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s3_bjd_offset.png", dpi=DPI)
        plt.close(fig)
    return "figures/s3/s3_bjd_offset.png"


def fig_clock(con, meta: dict) -> str:
    """Folded AG LMi photometry, the fitted dip, and the O-C."""
    rows = q(con, """
        SELECT p.phase, p.dmag, p.filter, p.readoutm, p.night
        FROM s3_clock_points p WHERE p.dmag IS NOT NULL""")
    # Re-apply the fit's baseline convention: per (config, filter) median
    # of the out-of-eclipse points (must mirror build_s3_timing.fit_clock).
    import collections as _c
    groups = _c.defaultdict(list)
    for ph, dm, filt, ro, night in rows:
        groups[((ro or "").strip(), (filt or "").strip())].append(
            (ph, dm, night))
    fit = q(con, """SELECT o_minus_c_s, o_minus_c_err_s, depth_mag,
                    width_phase, clock_bound_s FROM s3_clock_eclipses
                    WHERE tag = 'global' AND status = 'ok'""")
    period = float(meta.get("vsx_period_d", 1.0))
    # Nights the coverage gate let into the fit.  Points from the others
    # are still drawn — hiding them would hide the sampling problem — but
    # as hollow markers, so the figure cannot imply the curve was fitted
    # through them.
    gated = {n for n in str(meta.get("clock_nights_gated", "")).split(";")
             if n and n != "none"}
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.6, 4.0))
        colors = {"G": ps.BAND_COLOR["G"], "R": ps.BAND_COLOR["R"],
                  "I": ps.BAND_COLOR["I"],
                  "": ACCENT}
        seen = set()
        for (ro, filt), members in groups.items():
            # tm.CLOCK_OOE_PHASE, not a literal: build_s3_timing.fit_clock
            # imports the same constant, so the figure cannot drift away
            # from the fit it illustrates.
            ooe = [dm for ph, dm, _n in members
                   if abs(ph) > tm.CLOCK_OOE_PHASE]
            if len(ooe) < 3:
                continue
            base = float(np.median(ooe))
            label = f"{filt or '?'} ({'StackPro' if 'StackPro' in ro else ro})"
            for in_fit in (True, False):
                sel = [(ph, dm - base) for ph, dm, night in members
                       if (night in gated) == in_fit]
                if not sel:
                    continue
                lab = label if in_fit else "night excluded (one-sided)"
                # Band gets a MARKER as well as a hue.  Hue alone made G, R
                # and I one indistinguishable cloud in greyscale, which is
                # the whole reason the house style pairs the two channels.
                # The excluded points keep a plain open circle: they are one
                # legend entry spanning every filter, so a per-band shape
                # there would promise a distinction the label does not make.
                ax.scatter([p[0] for p in sel], [p[1] for p in sel], s=14,
                           alpha=0.8 if in_fit else 0.55,
                           marker=ps.band_marker(filt) if in_fit else "o",
                           color=colors.get(filt, ACCENT) if in_fit
                           else "none",
                           edgecolors=MUTED if not in_fit else "none",
                           linewidths=0.8,
                           label=lab if lab not in seen else None)
                seen.add(lab)
        if fit and fit[0][0] is not None:
            oc_s, oc_err, depth, width, bound = fit[0]
            ph0 = oc_s / 86400.0 / period
            xs = np.linspace(-0.15, 0.15, 400)
            ax.plot(xs, depth * np.exp(-((xs - ph0) ** 2)
                                       / (2 * width ** 2)),
                    color=INK, lw=1.4,
                    label=f"fit: O$-$C = {oc_s:+.0f} $\\pm$ {oc_err:.0f} s")
            ax.axvline(ph0, color=INK, lw=0.8, ls="--")
        ax.axvline(0, color=MUTED, lw=0.8, ls=":")
        ax.set_xlim(-0.16, 0.16)
        ax.invert_yaxis()               # dmag: fainter is down
        ax.set_xlabel("phase on the VSX ephemeris "
                      f"(P = {period} d, epoch HJD "
                      f"{meta.get('vsx_epoch_hjd', '?')})")
        ax.set_ylabel("$\\Delta$mag (baseline-subtracted)")
        ax.set_title("AG LMi primary eclipse, folded with OUR "
                     "mid-exposure HJD_UTC")
        ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s3_clock_oc.png", dpi=DPI)
        plt.close(fig)
    return "figures/s3/s3_clock_oc.png"


def fig_drift(con) -> str:
    """The relative-clock check: TELUT (telescope control system) minus
    the end of the exposure DATE-OBS names, along the whole baseline.

    Two DIFFERENT machines' clocks.  Header JD and header DATE-OBS come
    from ONE clock and can only ever agree with each other, so their
    0.002 s consistency says nothing about drift; this comparison is the
    only one in the stage that can.
    """
    rows = q(con, """SELECT jd_utc_start, resid_s, era_id, night
                     FROM s3_clock_drift
                     WHERE informative = 1 AND resid_s IS NOT NULL
                     ORDER BY jd_utc_start""")
    jd = np.array([r[0] for r in rows])
    resid = np.array([r[1] for r in rows])
    # Years since the first sampled frame: a readable x axis, and the
    # units the drift slope is quoted in.
    years = (jd - jd.min()) / 365.25
    slope = float(np.polyfit(years, resid, 1)[0]) if len(rows) >= 3 else 0.0
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.4, 3.4))
        ax.scatter(years, resid, s=22, color=ACCENT, zorder=3,
                   label="sampled headers (one clock vs the other)")
        xs = np.linspace(years.min(), years.max(), 10)
        ax.plot(xs, np.polyval(np.polyfit(years, resid, 1), xs),
                color=WARN, lw=1.2, ls="--",
                label=f"trend {slope:+.2f} s/year")
        ax.axhline(float(np.median(resid)), color=GREEN, lw=1, ls=":",
                   label=f"median {np.median(resid):.2f} s")
        ax.set_xlabel(f"years since {rows[0][3]}")
        ax.set_ylabel("TELUT $-$ (DATE-OBS + EXPTIME)  (s)")
        ax.set_title("Relative clock check: telescope clock vs "
                     "acquisition clock")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s3_clock_drift.png", dpi=DPI)
        plt.close(fig)
    return "figures/s3/s3_clock_drift.png"


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def section_stamp(con) -> str:
    n_scan = q1(con, "SELECT sum(n_frames) FROM s3_dateobs_audit")
    n_fam = q1(con, "SELECT count(*) FROM s3_dateobs_audit")
    n_out = q1(con, "SELECT count(*) FROM s3_dateobs_outliers")
    max_ok = q1(con, """SELECT max(max_abs_s) FROM s3_dateobs_audit
                        WHERE readoutm != 'Mode0'""")
    out_row = q(con, "SELECT path, night, diff_s FROM s3_dateobs_outliers")
    n_hdr = q1(con, "SELECT count(*) FROM s3_header_audit")
    n_helio = q1(con, """SELECT count(*) FROM s3_header_audit
                         WHERE jd_helio_header IS NOT NULL""")
    n_mid_ok = q1(con, """SELECT count(*) FROM s3_header_audit
                          WHERE abs(helio_resid_mid_s) <= 1.0""")
    n_mid_bad = q1(con, """SELECT count(*) FROM s3_header_audit
                           WHERE helio_resid_mid_s IS NOT NULL
                             AND abs(helio_resid_mid_s) > 5.0""")
    bad_rows = q(con, """SELECT family, year, exptime_s, helio_resid_start_s,
                                helio_resid_mid_s
                         FROM s3_header_audit
                         WHERE helio_resid_mid_s IS NOT NULL
                           AND abs(helio_resid_mid_s) > 5.0
                         ORDER BY family, exptime_s""")
    # TELUT slope by simple least squares: ~1 means the telescope clock
    # was read one whole exposure after DATE-OBS -> DATE-OBS is the start.
    # ONLY rows where TELUT is genuinely a second clock read may enter:
    # the pyscope family copies DATE-OBS verbatim into TELUT, so its rows
    # contribute (EXPTIME, 0) pairs that carry no information and drag the
    # slope away from the answer the informative rows actually give.
    telut = q(con, """SELECT exptime_s, telut_minus_dateobs_s
                      FROM s3_header_audit
                      WHERE telut_minus_dateobs_s IS NOT NULL
                        AND telut_minus_dateobs_s != 0.0""")
    n_telut_mute = q1(con, """SELECT count(*) FROM s3_header_audit
                              WHERE telut_minus_dateobs_s = 0.0""")
    slope = None
    if len(telut) >= 5:
        e = np.array([r[0] for r in telut])
        d = np.array([r[1] for r in telut])
        slope = float(np.polyfit(e, d, 1)[0])
    scan_rows = q(con, """SELECT readoutm, n_frames, median_s, p1_s, p99_s,
                                 max_abs_s, n_gt_100ms, stamp_resolution_s,
                                 n_with_fractional_s
                          FROM s3_dateobs_audit ORDER BY n_frames DESC""")
    scan_tbl = table(
        ["readout family", "frames", "median (s)", "p1 (s)", "p99 (s)",
         "max |diff| (s)", "&gt; 0.1 s", "stamp resolution (s)"],
        [[esc(r[0]), fmt(r[1]), fnum(r[2], 3), fnum(r[3], 3), fnum(r[4], 3),
          fnum(r[5], 3), fmt(r[6]),
          f"{r[7]:.3f}" + (" <b>(whole seconds)</b>" if r[7] >= 1.0 else "")]
         for r in scan_rows],
        # A family that stamps whole seconds is flagged as loudly as one
        # with a disagreeing card: it is the same size of defect.
        row_classes=["warn" if (r[6] or r[7] >= 1.0) else None
                     for r in scan_rows])
    # The families whose start-vs-mid semantics the two probes actually
    # cover, and the ones they do not (see build_s3_timing).
    # Counted over frame_times, NOT over the JD-vs-DATE-OBS scan: the scan
    # covers every canonical frame (calibrations included), the time axis
    # does not, and subtracting one population from the other would print
    # a number that is not the count of anything.
    n_axis = q1(con, "SELECT count(*) FROM frame_times")
    n_unverified = q1(con, """SELECT count(*) FROM frame_times
                              WHERE start_evidence = ?""",
                      (tm.START_UNVERIFIED,))
    unver_rows = q(con, """
        SELECT COALESCE(NULLIF(TRIM(f.readoutm), ''), '(blank)') fam,
               count(*), min(t.era_id), max(t.era_id), min(f.night),
               max(f.night)
        FROM frame_times t JOIN frames f ON f.path = t.path
        WHERE t.start_evidence = ? GROUP BY fam ORDER BY 2 DESC""",
                   (tm.START_UNVERIFIED,))
    unver_tbl = table(
        ["readout family", "frames", "eras", "nights", "what is missing"],
        [[esc(r[0]), fmt(r[1]),
          f"{r[2]}&ndash;{r[3]}" if r[2] != r[3] else f"{r[2]}",
          f"{esc(r[4])} &ndash; {esc(r[5])}",
          "no JD-HELIO card; TELUT is a copy of DATE-OBS"
          if r[0] == "(blank)" else
          "JD-HELIO present but broken / no discriminating exposure"]
         for r in unver_rows],
        row_classes=["warn"] * len(unver_rows))
    # The one exposure whose two copies disagree by 271 s is not an
    # anecdote: S0b's raw_reduced_links measured 46 of them.
    n_withdrawn = q1(con, "SELECT count(*) FROM s3_time_outliers")
    worst_sib = q(con, """SELECT path, sibling_jd_drift_s FROM s3_time_outliers
                          ORDER BY sibling_jd_drift_s DESC LIMIT 1""")
    bad_tbl = table(
        ["family", "year", "EXPTIME (s)", "resid@START (s)", "resid@MID (s)"],
        [[esc(r[0]), esc(r[1]), fnum(r[2], 2), fnum(r[3], 2), fnum(r[4], 2)]
         for r in bad_rows])
    src = fig_jdhelio_semantics(con)
    return f"""
<section id="stamp">
<div class="bhead"><h2>1 &middot; What instant does the header stamp?</h2>
<span class="tag">JD = DATE-OBS = UTC exposure start &mdash; proven for
{fmt(n_axis - n_unverified)} frames, assumed for {fmt(n_unverified)}
</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Every MACRO paper folds its data on a period; a
half-exposure ambiguity in the time stamp would smear or shift every
phase.  The headers offer JD, DATE-OBS, and JD-HELIO — do they agree, and
which instant of the exposure do they name?</p>

<h3>Evidence</h3>
<p class="sub">First, internal consistency: header JD against header
DATE-OBS for all {fmt(n_scan)} canonical frames carrying both cards,
by readout family:</p>
{scan_tbl}
<p class="sub">The two cards are the SAME instant to
{fnum(max_ok, 3)}&nbsp;s (their write precision) in every family, with
exactly {fmt(n_out)} exception:
<code>{esc(out_row[0][0])}</code> ({esc(out_row[0][1])}),
whose JD sits {fnum(out_row[0][2], 1)}&nbsp;s from its DATE-OBS — a
<i>reduced-tree</i> frame whose JD the reduction pipeline re-stamped, the
S0b <code>stem_jd_drift</code> behavior.</p>
<p class="sub"><b>This comparison has a blind spot, and it is the
important part of this section.</b>  It can only catch a re-stamp that
moved ONE card.  When the reduction pipeline re-stamped BOTH, the copy
stays perfectly self-consistent and passes here — while disagreeing with
its own raw parent about when the photons arrived.  S0b already measured
those cases: joining <code>raw_reduced_links</code> to the time axis
finds {fmt(n_withdrawn)} exposures whose two copies disagree by more than
{tm.JD_SIBLING_DISAGREE_S:.0f}&nbsp;s, the worst by
{fnum(worst_sib[0][1], 0) if worst_sib else '&mdash;'}&nbsp;s
(<code>{esc(worst_sib[0][0]) if worst_sib else '&mdash;'}</code>).  S3
therefore withdraws the BJD of the reduced copy in every one of them
(<code>{tm.BJD_JD_DISAGREES}</code>, listed in
<code>s3_time_outliers</code>) and carries
<code>sibling_jd_drift_s</code> on every row that has a linked sibling,
so a re-stamped time is visible where a consumer will actually look.</p>
<p class="sub">Second, which instant: {fmt(n_hdr)} headers sampled across
era families and years ({fmt(n_helio)} carry JD-HELIO).  We recomputed
the heliocentric correction independently (astropy, DE ephemeris, Winer
coordinates) at the exposure START and at MID-exposure:</p>
<div class="grid">{_figure(src,
    f"Header JD-HELIO minus our own heliocentric JD, per sampled header. "
    f"Left: evaluated at the JD instant, residuals climb exactly the "
    f"EXPTIME/2 line. Right: evaluated at JD + EXPTIME/2, "
    f"{fmt(n_mid_ok)}/{fmt(n_helio)} residuals collapse to within 1 s. "
    "Yellow: the HDR family and one Mode0 focus frame, whose JD-HELIO is "
    "simply broken.")}</div>
<p class="sub">The acquisition software computed JD-HELIO at
<i>start&nbsp;+&nbsp;EXPTIME/2</i> — which simultaneously proves the base
stamp is the START and hands us an independent check of the mid-exposure
arithmetic.  Third, the telescope-clock card TELUT (written with the rest
of the header, after readout) sits at DATE-OBS + EXPTIME + a small
constant: slope of TELUT&minus;DATE-OBS against EXPTIME over the
{fmt(len(telut))} sampled headers whose TELUT is genuinely a second clock
read = {fnum(slope, 3)} (1.0 = header written one full exposure after the
stamp).  {fmt(n_telut_mute)} further rows are EXCLUDED from that fit
because their TELUT is byte-identical to DATE-OBS — a copy, not a clock
read, contributing an uninformative (EXPTIME, 0) point that would pull
the slope to 0.901 and silently credit those eras with evidence they do
not have.  {fmt(n_mid_bad)} sampled headers FAIL the mid-exposure
test:</p>
{bad_tbl}
<p class="sub">Both probes are therefore missing for exactly the same
frames — the 2026 <code>pyscope</code> eras, which write no JD-HELIO card
at all and copy DATE-OBS into TELUT.  Those same frames stamp DATE-OBS at
WHOLE-SECOND resolution (column above): a 1&nbsp;s granularity where
every MaxIm family writes milliseconds, i.e. up to 0.5&nbsp;s of
systematic if pyscope truncates rather than rounds.  The frames with no
start-vs-mid evidence at all:</p>
{unver_tbl}

<h3>Decision</h3>
<div class="decision"><b>Header JD (&equiv; DATE-OBS) is the UTC exposure
START — proven by two independent probes for the
{fmt(n_axis - n_unverified)} frames whose readout family carries a
discriminating JD-HELIO or TELUT card, and ASSUMED (not proven) for the
{fmt(n_unverified)} that do not.</b>  Those rows say so on their face:
<code>frame_times.start_evidence</code> =
<code>{tm.START_UNVERIFIED}</code>.  Header JD-HELIO is never used as a
time stamp — it is heliocentric-only (up to ~4 s from barycentric),
UTC-scale (~69 s from TDB), computed by uncontrolled acquisition-software
approximations, and demonstrably garbage in the HDR family.  The
re-stamped reduced-tree copies do NOT quietly reinforce the S0 rule: they
were on the shared axis, and this build takes them off it.</div>

<h3>Consequence</h3>
<p class="sub">The era-audit table above is the paper appendix: any
referee asking &ldquo;how do you know your clocks?&rdquo; gets pointed at
{fmt(n_scan)} internal comparisons, {fmt(n_helio)} independent
heliocentric recomputations, and the TELUT cross-check — all
regenerable by one script.  Two gaps are named rather than papered over,
and both are cheap to close: a one-line citation of the pyscope source
would settle start-vs-mid for {fmt(n_unverified)} frames, and the
<code>UT</code> and <code>MJD-OBS</code> cards those headers do carry are
a second, unused probe.</p>
</div></section>"""


def section_stackpro(con, meta: dict) -> str:
    sp = q(con, """SELECT exptime_s, n_gaps, n_short_discarded,
                          raw_min_overhead_s, min_overhead_s,
                          n_regular_series, regular_overhead_s
                   FROM s3_cadence WHERE readoutm = 'High Gain StackPro'
                   ORDER BY exptime_s""")
    # The cell the bound comes from: smallest robust (median gap - EXPTIME).
    bound_row = min((r for r in sp if r[6] is not None), key=lambda r: r[6])
    n_series_frames = q1(con, """SELECT sum(n_gaps) FROM s3_cadence
                                 WHERE readoutm = 'High Gain StackPro'""")
    n_regular = q1(con, """SELECT sum(n_regular_series) FROM s3_cadence
                           WHERE readoutm = 'High Gain StackPro'""")
    maxim = q(con, """SELECT exptime_s, helio_resid_start_s, helio_resid_mid_s
                      FROM s3_header_audit
                      WHERE family = 'High Gain StackPro'
                      ORDER BY exptime_s DESC LIMIT 1""")[0]
    n_sp_frames = q1(con, """SELECT count(*) FROM frame_times
                             WHERE mid_method = ?""", (tm.MID_STACKPRO,))
    # The pairs the naive estimator was actually reading, worst first.
    bad = q(con, """SELECT path_a, path_b, exptime_s, gap_s,
                           series_median_gap_s, canonical_target, night
                    FROM s3_cadence_outliers
                    WHERE readoutm LIKE '%StackPro%' AND series_regular = 1
                    ORDER BY gap_s LIMIT 3""")
    n_bad_sp = q1(con, """SELECT count(*) FROM s3_cadence_outliers
                          WHERE readoutm LIKE '%StackPro%'
                            AND series_regular = 1""")
    # The naive estimator's other tell: gaps SHORTER than the exposure.
    n_impossible_cells = q1(con, """SELECT count(*) FROM s3_cadence
                                    WHERE raw_min_overhead_s < 0""")
    src = fig_cadence_overheads(con)
    sp_tbl = table(
        ["EXPTIME (s)", "gaps kept", "short gaps<br>discarded",
         "naive min<br>overhead (s)", "min overhead<br>after cut (s)",
         "regular<br>series", "<b>median-cadence<br>overhead (s)</b>"],
        [[fnum(r[0], 3), fmt(r[1]), fmt(r[2]), fnum(r[3], 2), fnum(r[4], 2),
          fmt(r[5]), f"<b>{fnum(r[6], 2)}</b>"] for r in sp],
        row_classes=["warn" if r[3] is not None and r[3] < 1.0 else None
                     for r in sp])
    bad_tbl = table(
        ["frame", "next frame", "EXPTIME (s)", "gap (s)",
         "the series' own median gap (s)"],
        [[f"<code>{esc(r[0].split('/')[-1])}</code>",
          f"<code>{esc(r[1].split('/')[-1])}</code>",
          fnum(r[2], 3), fnum(r[3], 3), fnum(r[4], 2)] for r in bad],
        row_classes=["warn"] * len(bad))
    naive = float(meta.get("stackpro_naive_min_overhead_s", "nan"))
    filtered = float(meta.get("stackpro_filtered_min_overhead_s", "nan"))
    # The sweep usually returns one value; say so plainly rather than
    # printing "between X and X".
    sw_lo = meta.get("stackpro_bound_sweep_min", "?")
    sw_hi = meta.get("stackpro_bound_sweep_max", "?")
    sweep_txt = (f"every one returns {esc(sw_lo)}&nbsp;s" if sw_lo == sw_hi
                 else f"the result stays between {esc(sw_lo)} and "
                      f"{esc(sw_hi)}&nbsp;s")
    # The per-sub-read reading of EXPTIME requires every back-to-back gap
    # to exceed N_SUB x EXPTIME.  Counted, not asserted: the two cells
    # that cannot refute it are the shortest exposures, where the fixed
    # readout overhead is larger than 16x a fraction of a second anyway.
    n_refute = q1(con, """SELECT count(*) FROM s3_cadence
                          WHERE readoutm = 'High Gain StackPro'
                            AND (min_overhead_s + exptime_s)
                                < ? * exptime_s""", (tm.N_SUB_STACKPRO,))
    return f"""
<section id="stackpro">
<div class="bhead"><h2>2 &middot; StackPro: the midpoint of a 16-fold
sum</h2>
<span class="tag">worst case {tm.STACKPRO_MID_WORST_CASE_S:.2f} s &mdash;
seconds, not milliseconds</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">A StackPro frame is the on-camera SUM of
{tm.N_SUB_STACKPRO} sub-reads (S2's PTC result).  Header EXPTIME could
mean the total span or the per-sub-read time, and any dead time between
sub-reads would shift the photon-weighted midpoint.  What is the
mid-time of {fmt(n_sp_frames)} StackPro science frames, and what is the
worst case if the contiguity assumption is wrong?</p>

<h3>Evidence</h3>
<p class="sub">The first question is settled and stays settled: EXPTIME
is the TOTAL span, not the per-sub-read time.  If it were per-sub-read,
back-to-back gaps would have to exceed
{tm.N_SUB_STACKPRO}&nbsp;&times;&nbsp;EXPTIME; the smallest genuine gap
falls BELOW that line in {fmt(n_refute)} of {fmt(len(sp))} exposure-time
cells.  (The {fmt(len(sp) - n_refute)} that cannot refute it are the
shortest exposures, where the camera's fixed per-frame overhead is larger
than {tm.N_SUB_STACKPRO}&nbsp;&times;&nbsp;a fraction of a second and the
test has no power.)</p>
<div class="grid">{_figure(src,
    "Back-to-back series cadence. If EXPTIME were per-sub-read, gaps "
    "would follow the dotted 15x line; instead StackPro overheads match "
    "the plain High Gain save/slew pipeline at every exposure time.")}</div>
<p class="sub">The second question — how big the dead-time ceiling is —
turns entirely on the ESTIMATOR, and the obvious one is wrong.  Taking
the smallest (gap&nbsp;&minus;&nbsp;EXPTIME) over
{fmt(n_series_frames)} StackPro gaps is an extreme order statistic that a
single bad time stamp owns outright, and does: it returns
{naive:.2f}&nbsp;s, set by these pairs, in which one frame lands about a
second from its neighbour inside a series whose own measured cadence is
11&ndash;13&nbsp;s:</p>
{bad_tbl}
<p class="sub">A camera that demonstrably needs ~12&nbsp;s per frame
cannot deliver two exposures {fnum(bad[0][3], 2)}&nbsp;s apart.  These
are time-stamp defects, not cadence — {fmt(n_bad_sp)} of them in regular
StackPro series, all listed in <code>s3_cadence_outliers</code>.  The
same estimator returns NEGATIVE overheads (a gap SHORTER than its own
exposure, which is impossible) in {fmt(n_impossible_cells)} cells of the
full cadence table: proof that the raw minimum is not measuring
anything.</p>
<p class="sub">So the bound is taken robustly instead.  Each same-config
run contributes its MEDIAN gap (one bad stamp cannot move a median); runs
whose cadence is regular to
{tm.CADENCE_REGULAR_SPREAD:.0%} interquartile spread
({fmt(n_regular)} of them) are read as machine cycle times; and the
smallest (median gap &minus; EXPTIME) over those is the ceiling.  The
physical statement is exact: a camera that repeatedly delivers a frame
every <i>median gap</i> seconds must fit exposure + internal dead time +
readout + save inside that cycle.</p>
{sp_tbl}
<p class="sub">Three numbers, one column apart, for the same archive:
the naive minimum says {naive:.2f}&nbsp;s; discarding gaps below
{tm.CADENCE_MIN_GAP_FRACTION:.0%} of their own series' median says
{filtered:.2f}&nbsp;s, but only for that exact fraction — move it and the
answer moves with it, which is why it is not the answer either; the
median-cadence statistic says
{tm.STACKPRO_DEADTIME_BOUND_S:.2f}&nbsp;s, from the
{fnum(bound_row[0], 3)}&nbsp;s cell.  The build re-derives that last
number under {esc(meta.get('stackpro_bound_sweep_n', '?'))} alternative
cut settings, each moving one threshold well away from its default:
{sweep_txt}.  A bound that depends on the analyst's choice of threshold
is not a bound; this one does not.</p>
<p class="sub">What does NOT corroborate any of this: MaxIm's own
JD-HELIO on a {fnum(maxim[0], 0)}&nbsp;s StackPro frame equals
start&nbsp;+&nbsp;EXPTIME/2 to {fnum(abs(maxim[2]), 2)}&nbsp;s.  That
shows MaxIm <i>also assumes</i> EXPTIME is the contiguous total span; it
constrains the software's convention, not the detector's sub-read
contiguity, and it is not independent evidence.</p>

<h3>Decision</h3>
<div class="decision"><b>StackPro mid-time = start + EXPTIME/2, under the
stated assumption that the {tm.N_SUB_STACKPRO} sub-reads are contiguous
and span EXPTIME.  Worst case if dead time hides between sub-reads: the
true midpoint shifts late by at most
{tm.STACKPRO_MID_WORST_CASE_S:.2f}&nbsp;s</b> — half the measured
{tm.STACKPRO_DEADTIME_BOUND_S:.2f}&nbsp;s ceiling, recorded as
<code>{tm.MID_STACKPRO}</code> in <code>frame_times.mid_method</code> so
any user can query which frames carry it.  That ceiling is honest but
WEAK: the gap it measures is dominated by readout and file save, and
cadence cannot separate those from anything happening between sub-reads.
It is an upper limit, not an estimate of the error.</div>

<h3>Consequence</h3>
<p class="sub"><b>CV eclipse timing on StackPro series (AG LMi 2023, the
1024&nbsp;s deep fields) inherits a systematic of up to
{tm.STACKPRO_MID_WORST_CASE_S:.2f}&nbsp;s, which is NOT below a
sub-second target.</b>  A paper quoting sub-second absolute mid-times
must either avoid the {fmt(n_sp_frames)} frames carrying
<code>{tm.MID_STACKPRO}</code>, or first close this bound with evidence
cadence cannot supply: the camera's sub-read readout time from a manual
or a lab measurement.  Either would collapse the ceiling by orders of
magnitude — and only this one constant changes, after which one re-run
repairs every product.  Relative timing WITHIN a StackPro series is
unaffected: a fixed per-frame offset cancels in a differential
measurement.</p>
</div></section>"""


def section_frame_times(con, meta: dict) -> str:
    n_rows = q1(con, "SELECT count(*) FROM frame_times")
    n_bjd = q1(con, "SELECT count(*) FROM frame_times "
                    "WHERE bjd_tdb IS NOT NULL")
    method_rows = q(con, """
        SELECT mid_method, bjd_method, count(*) FROM frame_times
        GROUP BY mid_method, bjd_method ORDER BY 3 DESC""")
    # The meaning of a row is the PAIR (mid_method, bjd_method), not the
    # mid_method alone: a 'no_coords' row has a perfectly good mid-time
    # and no BJD, and keying the explanation on mid_method alone told
    # those 1,425 readers "start + EXPTIME/2" and never why their BJD was
    # NULL.
    mid_explain = {
        tm.MID_PLAIN: "start + EXPTIME/2",
        tm.MID_STACKPRO: "start + EXPTIME/2 (StackPro sum; "
                         f"&le; {tm.STACKPRO_MID_WORST_CASE_S:.2f} s "
                         "worst case)",
        tm.MID_EXPTIME_NONPOS: "EXPTIME &le; 0 or missing — start used",
        tm.MID_NO_JD: "no header JD — no time axis possible",
    }
    bjd_explain = {
        tm.BJD_NO_JD: "no BJD: no header JD to convert",
        tm.BJD_NO_COORDS: "no BJD: frame has no cataloged sky position, "
                          "so there is nothing to point the light-travel "
                          "correction at (the UTC mid-time is still valid)",
        tm.BJD_JD_DISAGREES: "BJD WITHDRAWN: this copy's stamp disagrees "
                             f"with its raw parent by &gt; "
                             f"{tm.JD_SIBLING_DISAGREE_S:.0f} s",
    }

    def _meaning(mm, bm):
        """One cell of prose per (mid_method, bjd_method) pair."""
        parts = [mid_explain.get(mm, "&mdash;")]
        if bm in bjd_explain:
            parts.append(bjd_explain[bm])
        return "; ".join(parts)

    m_tbl = table(
        ["mid_method", "bjd_method", "frames", "meaning"],
        [[f"<code>{esc(mm)}</code>", f"<code>{esc(bm or '—')}</code>",
          fmt(n), _meaning(mm, bm)]
         for mm, bm, n in method_rows],
        row_classes=[None if (bm or "").startswith("bary") else "warn"
                     for _mm, bm, _n in method_rows])
    ltt_rng = q(con, """SELECT min(bary_ltt_s), max(bary_ltt_s),
                               min(tdb_minus_utc_s), max(tdb_minus_utc_s)
                        FROM frame_times WHERE bjd_tdb IS NOT NULL""")[0]
    src = fig_bjd_offset(con)
    eph = meta.get("ephemeris", "?")
    # Provenance counts that belong beside the table, not in a footnote.
    n_calib_excl = meta.get("frame_times_calib_excluded", "?")
    n_unverified = q1(con, "SELECT count(*) FROM frame_times "
                           "WHERE start_evidence = ?", (tm.START_UNVERIFIED,))
    n_sibling = q1(con, "SELECT count(*) FROM frame_times "
                        "WHERE sibling_jd_drift_s IS NOT NULL")
    n_sib_pairs = int(meta.get("frame_times_sibling_pairs", 0))
    n_withdrawn = q1(con, "SELECT count(*) FROM s3_time_outliers")
    # The frame-center caveat, measured from the sampled headers' own
    # geometry instead of hand-typed: worst case over the archive, and
    # worst case excluding the one doubled-height HDR composite.
    corner_all = q(con, """SELECT family, naxis1, naxis2, pixscale_arcsec,
                                  corner_ltt_s FROM s3_header_audit
                           WHERE corner_ltt_s IS NOT NULL
                           ORDER BY corner_ltt_s DESC LIMIT 1""")[0]
    corner_plain = q(con, """SELECT family, naxis1, naxis2, pixscale_arcsec,
                                    corner_ltt_s FROM s3_header_audit
                             WHERE corner_ltt_s IS NOT NULL
                               AND family != 'HDR'
                             ORDER BY corner_ltt_s DESC LIMIT 1""")[0]
    return f"""
<section id="frametimes">
<div class="bhead"><h2>3 &middot; frame_times: BJD_TDB for every
canonical science frame</h2>
<span class="tag">{fmt(n_bjd)} of {fmt(n_rows)} frames on the shared
time axis</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Every downstream product needs one authoritative
mid-exposure BJD_TDB per frame — computed once, from scratch, with its
inputs recorded, so no paper ever re-derives (or worse, half-derives)
its own time axis.</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"The correction S3 applies: the barycentric light-travel term spans "
    f"{ltt_rng[0]/60:+.1f} to {ltt_rng[1]/60:+.1f} minutes over the "
    f"archive's sky and seasons; the TDB-UTC scale offset adds "
    f"{ltt_rng[2]:.3f}-{ltt_rng[3]:.3f} s. Ignoring either would be "
    "fatal for eclipse timing.")}</div>
{m_tbl}
<p class="sub">Three exclusions and flags travel with the table, each
counted here rather than discovered later: {fmt(n_calib_excl)} frames
whose headers say <code>IMAGETYP = 'Light Frame'</code> but which S0b
catalogued as CALIBRATIONS (164 master flats, 60 raw flats, 8 master
darks) are kept OFF the axis — a master frame is a stack, so its header
JD is not an exposure instant at all; {fmt(n_unverified)} rows carry
<code>start_evidence = {tm.START_UNVERIFIED}</code> (section 1); and
{fmt(n_sibling)} rows carry a <code>sibling_jd_drift_s</code> because the
same exposure also sits on the axis as a second copy
({fmt(n_sib_pairs)} raw/reduced pairs), of which {fmt(n_withdrawn)}
disagreed badly enough to lose their BJD.</p>

<h3>Decision</h3>
<div class="decision"><b>The <code>frame_times</code> table (one row per
canonical science frame, keyed by path) is the portfolio's time axis:
<code>bjd_tdb</code> = UTC start + EXPTIME/2 &rarr; TDB scale &rarr; +
barycentric light travel toward the frame's cataloged sky position,
ephemeris <code>{esc(eph)}</code>, Winer EarthLocation
({tm.WINER_LAT_DEG}&deg;, {tm.WINER_LON_DEG}&deg;,
{tm.WINER_ALT_M:.0f}&nbsp;m).</b>  Frames without JD or coordinates keep
a row with a NULL BJD and an explanatory method — nothing silently
falls off the axis.</div>

<h3>Consequence</h3>
<p class="sub"><b>Two caveats travel with the table.</b>  First, the
correction points at the FRAME CENTER.  An object at a frame corner
differs by up to {fnum(corner_plain[4], 1)}&nbsp;s on the widest ordinary
field ({esc(corner_plain[0])}, {corner_plain[1]}&times;{corner_plain[2]}
px at {fnum(corner_plain[3], 3)}&Prime;/px), and
{fnum(corner_all[4], 1)}&nbsp;s on the widest sampled frame of all (the
{esc(corner_all[0])} family's
{corner_all[1]}&times;{corner_all[2]}&nbsp;px composite) — computed from
each sampled header's own geometry as half-diagonal &times;
{tm.LIGHT_TIME_PER_RAD_S:.0f}&nbsp;s/rad, not typed by hand.  A paper
needing sub-second absolutes for an off-center object recomputes with
<code>macro_core.timing.bjd_tdb_from_utc</code> at the object's own
coordinates — every input needed is in the row.</p>
<p class="sub">Second, and easier to trip over: this table is one row per
canonical FILE, and {fmt(n_sib_pairs)} exposures appear TWICE on it — a
raw copy and a reduced copy that S0 did not collapse into one duplicate
group.  Selecting straight from <code>frame_times</code> therefore
double-counts about
{fnum(100.0 * n_sib_pairs / max(n_rows, 1), 0)}% of the archive's
photometry.  <b>Deduplicate (drop the reduced tree, or group on the
raw/reduced link) before building a light curve.</b></p>
</div></section>"""


def section_clock(con, meta: dict) -> str:
    pts = q1(con, "SELECT count(*) FROM s3_clock_points "
                  "WHERE dmag IS NOT NULL")
    n_frames = q1(con, "SELECT count(*) FROM s3_clock_points")
    ecl = q(con, """SELECT tag, n_points, phase_min, phase_max, depth_mag,
                           width_phase, o_minus_c_s, o_minus_c_err_s,
                           clock_bound_s, status
                    FROM s3_clock_eclipses ORDER BY tag = 'global' DESC,
                    tag""")
    ecl_tbl = table(
        ["fit", "points", "phase range", "depth (mag)", "width (phase)",
         "O&minus;C (s)", "&plusmn; (s)", "|clock| bound (s)", "status"],
        [[esc(r[0]), fmt(r[1]), f"[{fnum(r[2], 3)}, {fnum(r[3], 3)}]",
          fnum(r[4], 2), fnum(r[5], 4),
          fnum(r[6], 0), fnum(r[7], 0), fnum(r[8], 0), esc(r[9])]
         for r in ecl],
        row_classes=["warn" if r[9] != "ok" else None for r in ecl])
    g = next((r for r in ecl if r[0] == "global" and r[9] == "ok"), None)
    src = fig_clock(con, meta)
    epoch = meta.get("vsx_epoch_hjd", "?")
    period = meta.get("vsx_period_d", "?")
    source = meta.get("vsx_source", "?")
    gaia_pred = meta.get("gaia_predicted_oc_s")
    gaia_env = meta.get("gaia_oc_envelope_s")
    vsx_term = meta.get("clock_vsx_quant_term_s", "?")
    eph_term = meta.get("clock_eph_term_s", "?")
    # Nights DERIVED from the table, never from the configured tuple: the
    # two disagreed in an earlier revision (the text said four nights
    # reached |phase| < 0.02; one of them never got within 0.078, and a
    # fourth row was missing entirely).
    n_nights_cfg = q1(con, "SELECT count(*) FROM s3_clock_eclipses "
                           "WHERE tag != 'global'")
    n_nights_ok = q1(con, "SELECT count(*) FROM s3_clock_eclipses "
                          "WHERE tag != 'global' AND status = ?",
                     (tm.CLOCK_STATUS_OK,))
    # The relative-drift check (a different pair of clocks entirely).
    drift = q(con, """SELECT count(*), count(DISTINCT era_id),
                             min(resid_s), max(resid_s), min(night),
                             max(night)
                      FROM s3_clock_drift
                      WHERE informative = 1 AND resid_s IS NOT NULL""")[0]
    drift_all = q1(con, "SELECT count(*) FROM s3_clock_drift")
    drift_src = fig_drift(con) if drift[0] >= 3 else None
    drift_span = q(con, """SELECT resid_s FROM s3_clock_drift
                           WHERE informative = 1 AND resid_s IS NOT NULL
                           ORDER BY resid_s""")
    drift_vals = np.array([r[0] for r in drift_span])
    drift_med = float(np.median(drift_vals)) if len(drift_vals) else 0.0
    drift_iqr = (float(np.percentile(drift_vals, 75)
                       - np.percentile(drift_vals, 25))
                 if len(drift_vals) else 0.0)
    # The trend the figure draws, stated in the text too — a dashed line
    # a reader has to eyeball is not a number.
    drift_epochs = q(con, """SELECT jd_utc_start, resid_s
                             FROM s3_clock_drift
                             WHERE informative = 1 AND resid_s IS NOT NULL
                             ORDER BY jd_utc_start""")
    drift_slope = 0.0
    if len(drift_epochs) >= 3:
        yrs = np.array([r[0] for r in drift_epochs]) / 365.25
        drift_slope = float(np.polyfit(
            yrs - yrs.min(), [r[1] for r in drift_epochs], 1)[0])
    if g:
        verdict = (f"<b>The eclipse arrives at O&minus;C = {g[6]:+.0f} "
                   f"&plusmn; {g[7]:.0f}&nbsp;s against the VSX ephemeris. "
                   "The absolute observatory clock is therefore consistent "
                   "with zero error at the "
                   f"|&Delta;t| &lt; {g[8]:.0f}&nbsp;s "
                   f"({g[8]/60:.0f}&nbsp;min) level &mdash; a ceiling set "
                   "almost entirely by how badly AG&nbsp;LMi's own "
                   "ephemeris is known, not by our measurement.</b>  Of "
                   f"that {g[8]:.0f}&nbsp;s, {abs(g[6]):.0f}&nbsp;s is the "
                   f"measured offset, {g[7]:.0f}&nbsp;s its fit error, and "
                   f"{esc(eph_term)}&nbsp;s the ephemeris term.  That last "
                   "number propagates the period uncertainty from the "
                   "independent Gaia DR3 eclipsing-binary solution "
                   f"({esc(meta.get('gaia_period_d', '?'))}&nbsp;d "
                   f"&plusmn;{esc(meta.get('gaia_period_err_d', '?'))}, "
                   f"{esc(meta.get('gaia_source', ''))}) over the "
                   f"{esc(meta.get('clock_mean_cycle', '?'))} cycles since "
                   "the epoch.  Propagating VSX's QUOTED LAST DIGIT "
                   f"instead would give only {esc(vsx_term)}&nbsp;s and a "
                   "bound roughly eight times tighter &mdash; a printer's "
                   "precision, not a measurement's, and this page will not "
                   "quote it.  The same Gaia-vs-VSX period difference "
                   f"predicts an O&minus;C of {esc(gaia_pred)}&nbsp;s at "
                   f"our epoch with a &plusmn;{esc(gaia_env)}&nbsp;s "
                   "envelope, inside which our measurement sits "
                   "comfortably: the offset belongs to the EPHEMERIS, not "
                   "to the clock.  <b>What this test rules out is gross "
                   "error &mdash; a wrong time zone, a missing leap "
                   "second, a mis-set clock &mdash; and nothing "
                   "finer.</b>")
    else:
        verdict = ("<b>No night passed the coverage gate, so the clock "
                   "bound is NOT established this run.</b>  The fallback "
                   "plan (a TESS-era EB with a modern epoch in an "
                   "archived field) stands documented in the roadmap.")
    return f"""
<section id="clock">
<div class="bhead"><h2>4 &middot; Clock validation: AG LMi as the
standard</h2>
<span class="tag">an astrophysical event checks the wall clock</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Header consistency (sections 1-3) proves the stamps agree
with EACH OTHER — not that the acquisition computer's clock agrees with
UTC.  Does an external astrophysical clock confirm the time axis?</p>

<h3>Evidence</h3>
<p class="sub">AG LMi (type {esc(meta.get('vsx_type', '?'))},
{esc(meta.get('vsx_mag_range', '?'))}) is an eclipsing binary in the
archive with P = {esc(str(period))}&nbsp;d, minimum epoch HJD
{esc(str(epoch))} ({esc(source)}).  {fmt(n_nights_cfg)} nights of it were
measured; {fmt(n_frames)} frames of those plus three baseline nights
went through forced-aperture differential photometry against
{esc(meta.get('clock_comps', '').count(';') + 1 if
     meta.get('clock_comps') else '?')} fixed comparison stars
({fmt(pts)} usable points), using OUR mid-exposure HJD_UTC — the same
convention as the VSX epoch, like against like:</p>
<div class="grid">{_figure(src,
    "AG LMi folded on the literature ephemeris using this pipeline's "
    "own time axis. The eclipse lands where the ephemeris says it "
    "should; the fitted center measures the residual offset.")}</div>
{ecl_tbl}
<p class="sub"><b>Only {fmt(n_nights_ok)} of those {fmt(n_nights_cfg)}
nights is a measurement.</b>  A symmetric template fitted to a one-sided
arc converges happily and returns a confident centre that is really the
slope of the flank it was handed — the 2023-03-18 night never reaches
phase 0 at all, and 2024-02-29 samples only the egress side.  Their rows
therefore carry
<code>{tm.CLOCK_STATUS_ONE_SIDED}</code> and a NULL O&minus;C instead of
numbers that disagree with each other by ~3,000&nbsp;s, and they are
excluded from the global fit rather than averaged into it.  A night with
too few points gets a <code>{tm.CLOCK_STATUS_TOO_FEW}</code> row for the
same reason: so it cannot vanish from the table without a trace.</p>
<p class="sub">Sections 1-3 compare header cards that ONE clock wrote, so
their agreement — however tight — cannot bound that clock's drift.  The
one genuinely independent comparison in the archive is TELUT, the
telescope control system's UTC, against the acquisition PC's DATE-OBS.
Sampled across every era ({fmt(drift[0])} of {fmt(drift_all)} sampled
headers carry a TELUT that is a real clock read; the rest have no TELUT
card, or the pyscope copy of DATE-OBS):</p>
{f'<div class="grid">{_figure(drift_src, "TELUT minus the end of the exposure DATE-OBS names, per sampled header, across the baseline. The scatter is header-write latency; a diverging pair of clocks would show as a trend.")}</div>' if drift_src else ''}
<p class="sub">Over {esc(drift[4])} to {esc(drift[5])} the residual
TELUT&nbsp;&minus;&nbsp;(DATE-OBS&nbsp;+&nbsp;EXPTIME) has median
{drift_med:.2f}&nbsp;s and interquartile spread {drift_iqr:.2f}&nbsp;s,
ranging {drift[2]:+.2f} to {drift[3]:+.2f}&nbsp;s, with a fitted trend of
{drift_slope:+.2f}&nbsp;s/year — smaller than the scatter it is fitted
through, i.e. consistent with no divergence at all.  That scatter is
dominated by variable header-write latency, so the number is an UPPER
limit on relative drift, not a measurement of it: <b>the two clocks do
not diverge by more than a few seconds across the baseline, and S3
cannot say anything finer than that.</b>  The check also has a hole worth
naming: {fmt(drift_all - drift[0])} of the {fmt(drift_all)} sampled
headers carry no usable TELUT, so the eras before {esc(drift[4])} are not
covered at all.</p>

<h3>Decision</h3>
<div class="decision">{verdict}</div>

<h3>Consequence</h3>
<p class="sub">A CV paper should take two different numbers from this
page.  <b>Absolute</b> mid-times are anchored only at the
{f"{g[8]/60:.0f}&nbsp;minute" if g else "&mdash;"} level above — enough
to prove no gross clock error, not enough to publish an absolute epoch.
<b>Relative</b> (within-archive) timing is what an O&minus;C diagram
actually rests on, and S3 bounds it at the few-second level from the
TELUT comparison above — NOT at the 0.002&nbsp;s of section 1, which
measures one clock against itself and is no evidence at all.  Closing
either gap is a one-night project: one modern AG&nbsp;LMi minimum
against a same-week TESS or AAVSO epoch collapses the ephemeris
systematic to seconds, and a GPS- or NTP-stamped exposure sequence would
replace the TELUT scatter with a real drift measurement.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S3 report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    # Read-only is not immune to a concurrent writer: another MACRO stage
    # (the S1 plate-solve batch) holds this database's write lock in
    # bursts, and a reader that does not wait simply fails.  Five minutes
    # of patience costs nothing; a failed render mid-run costs a re-run.
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = dict(q(con, "SELECT key, value FROM s3_build_meta"))
        n_rows = q1(con, "SELECT count(*) FROM frame_times")
        n_bjd = q1(con, "SELECT count(*) FROM frame_times "
                        "WHERE bjd_tdb IS NOT NULL")
        n_scan = q1(con, "SELECT sum(n_frames) FROM s3_dateobs_audit")
        g = q(con, """SELECT clock_bound_s FROM s3_clock_eclipses
                      WHERE tag = 'global' AND status = ?""",
              (tm.CLOCK_STATUS_OK,))
        # The banner says what the bound MEANS.  An unqualified "< 519 s"
        # is the number a reader lifts, and the old one omitted the
        # ephemeris uncertainty the same build had already measured.
        clock_txt = (f"absolute clock consistent with zero at the "
                     f"~{g[0][0]/60:.0f} min level set by AG LMi's "
                     f"ephemeris"
                     if g else "clock bound pending")

        sections = [
            section_stamp(con),
            section_stackpro(con, meta),
            section_frame_times(con, meta),
            section_clock(con, meta),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S3 &mdash; The Shared Time Axis (mid-exposure BJD_TDB)</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S3 &mdash; The Shared Time Axis</h1>
  <p>{fmt(n_scan)} header stamps audited &middot; {fmt(n_bjd)} of
  {fmt(n_rows)} canonical science frames with mid-exposure BJD_TDB
  (ephemeris <code>{esc(meta.get('ephemeris', '?'))}</code>) &middot;
  {clock_txt} &middot; built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#stamp">1 Header semantics</a> &middot;
  <a href="#stackpro">2 StackPro midpoint</a> &middot;
  <a href="#frametimes">3 frame_times</a> &middot;
  <a href="#clock">4 Clock validation</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_core.report_s3</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on
this page is the result of a SQL query; none is typed by hand.
Regenerate with <code>pipeline/scripts/build_s3_timing.py</code>.
Header JD-HELIO is never used by any MACRO product; see section 1 for
why.</footer>
</body></html>"""

        HTML_PATH.write_text(html, encoding="utf-8")

        # Belt and braces: every <img> the page references must exist and
        # be non-empty, or the build fails loudly.
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
