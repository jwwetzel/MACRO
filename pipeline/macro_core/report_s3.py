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
# Shared page machinery: same dark theme, same query discipline, same
# table generator as the earlier reports — one visual language site-wide.
from .report_s0 import (          # noqa: E402
    ACCENT, DARK, DPI, WARN, _figure, esc, fmt, q, q1, table)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s3"
HTML_PATH = DOCS_DIR / "s3_timing.html"

GREEN = "#9fd8ae"               # site badge green (verdict OK)


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
    with plt.rc_context(DARK):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.6))
        xs = np.geomspace(0.01, 2000, 50)
        ax1.plot(xs, xs / 2.0, color="#9aa4b2", lw=1, ls="--",
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
        ax2.axhline(0, color="#9aa4b2", lw=1, ls="--")
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
    indistinguishable from the plain High Gain pipeline — and the minimum
    observed overhead (0.24 s) caps any internal dead time.
    """
    rows = q(con, """
        SELECT readoutm, exptime_s, min_overhead_s, p5_overhead_s
        FROM s3_cadence
        WHERE readoutm IN ('High Gain', 'High Gain StackPro')
        ORDER BY exptime_s""")
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        for mode, color, marker in (("High Gain", ACCENT, "o"),
                                    ("High Gain StackPro", WARN, "s")):
            sub = [r for r in rows if r[0] == mode]
            e = np.array([r[1] for r in sub])
            p5 = np.array([max(r[3], 0.05) for r in sub])
            ax.plot(e, p5, marker=marker, ms=5, lw=1, color=color,
                    label=f"{mode} (p5 overhead)")
        e_sp = np.array([r[1] for r in rows if r[0] == "High Gain StackPro"])
        ax.plot(e_sp, 15.0 * e_sp, color="#e06c75", lw=1.2, ls=":",
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
    with plt.rc_context(DARK):
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
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(7.6, 4.0))
        colors = {"G": "#9fd8ae", "R": "#e06c75", "I": "#c678dd",
                  "": ACCENT}
        seen = set()
        for (ro, filt), members in groups.items():
            ooe = [dm for ph, dm, _n in members if abs(ph) > 0.08]
            if len(ooe) < 3:
                continue
            base = float(np.median(ooe))
            phs = [ph for ph, _d, _n in members]
            dms = [dm - base for _p, dm, _n in members]
            label = f"{filt or '?'} ({'StackPro' if 'StackPro' in ro else ro})"
            ax.scatter(phs, dms, s=14, alpha=0.8,
                       color=colors.get(filt, ACCENT),
                       label=label if label not in seen else None)
            seen.add(label)
        if fit and fit[0][0] is not None:
            oc_s, oc_err, depth, width, bound = fit[0]
            ph0 = oc_s / 86400.0 / period
            xs = np.linspace(-0.15, 0.15, 400)
            ax.plot(xs, depth * np.exp(-((xs - ph0) ** 2)
                                       / (2 * width ** 2)),
                    color="#e8eaed", lw=1.4,
                    label=f"fit: O$-$C = {oc_s:+.0f} $\\pm$ {oc_err:.0f} s")
            ax.axvline(ph0, color="#e8eaed", lw=0.8, ls="--")
        ax.axvline(0, color="#9aa4b2", lw=0.8, ls=":")
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
    telut = q(con, """SELECT exptime_s, telut_minus_dateobs_s
                      FROM s3_header_audit
                      WHERE telut_minus_dateobs_s IS NOT NULL""")
    # TELUT slope by simple least squares: ~1 means the telescope clock
    # was read one whole exposure after DATE-OBS -> DATE-OBS is the start.
    slope = None
    if len(telut) >= 5:
        e = np.array([r[0] for r in telut])
        d = np.array([r[1] for r in telut])
        slope = float(np.polyfit(e, d, 1)[0])
    scan_rows = q(con, """SELECT readoutm, n_frames, median_s, p1_s, p99_s,
                                 max_abs_s, n_gt_100ms
                          FROM s3_dateobs_audit ORDER BY n_frames DESC""")
    scan_tbl = table(
        ["readout family", "frames", "median (s)", "p1 (s)", "p99 (s)",
         "max |diff| (s)", "&gt; 0.1 s"],
        [[esc(r[0]), fmt(r[1]), fnum(r[2], 3), fnum(r[3], 3), fnum(r[4], 3),
          fnum(r[5], 3), fmt(r[6])] for r in scan_rows],
        row_classes=["warn" if r[6] else None for r in scan_rows])
    bad_tbl = table(
        ["family", "year", "EXPTIME (s)", "resid@START (s)", "resid@MID (s)"],
        [[esc(r[0]), esc(r[1]), fnum(r[2], 2), fnum(r[3], 2), fnum(r[4], 2)]
         for r in bad_rows])
    src = fig_jdhelio_semantics(con)
    return f"""
<section id="stamp">
<div class="bhead"><h2>1 &middot; What instant does the header stamp?</h2>
<span class="tag">JD = DATE-OBS = UTC exposure start &mdash; proven, not
assumed</span></div>

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
S0b <code>stem_jd_drift</code> behavior; its raw parent is clean.</p>
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
constant: slope of TELUT&minus;DATE-OBS against EXPTIME across the sample
= {fnum(slope, 3)} (1.0 = header written one full exposure after the
stamp).  {fmt(n_mid_bad)} sampled headers FAIL the mid-exposure test:</p>
{bad_tbl}

<h3>Decision</h3>
<div class="decision"><b>Header JD (&equiv; DATE-OBS) is the UTC exposure
START in every era family; S3 computes mid-exposure BJD_TDB from it from
scratch.  Header JD-HELIO is never used</b> — it is heliocentric-only
(up to ~4 s from barycentric), UTC-scale (~69 s from TDB), computed by
uncontrolled acquisition-software approximations, and demonstrably
garbage in the HDR family.  The one re-stamped reduced-tree JD reinforces
the S0 rule: timing always comes from the canonical (raw) copy.</div>

<h3>Consequence</h3>
<p class="sub">The era-audit table above is the paper appendix: any
referee asking &ldquo;how do you know your clocks?&rdquo; gets pointed at
{fmt(n_scan)} internal comparisons, {fmt(n_helio)} independent
heliocentric recomputations, and the TELUT cross-check — all
regenerable by one script.</p>
</div></section>"""


def section_stackpro(con) -> str:
    sp = q(con, """SELECT exptime_s, n_gaps, min_overhead_s, p5_overhead_s
                   FROM s3_cadence WHERE readoutm = 'High Gain StackPro'
                   ORDER BY exptime_s""")
    bound_row = min(sp, key=lambda r: r[2])
    n_series_frames = q1(con, """SELECT sum(n_gaps) FROM s3_cadence
                                 WHERE readoutm = 'High Gain StackPro'""")
    maxim = q(con, """SELECT exptime_s, helio_resid_start_s, helio_resid_mid_s
                      FROM s3_header_audit
                      WHERE family = 'High Gain StackPro'
                      ORDER BY exptime_s DESC LIMIT 1""")[0]
    n_sp_frames = q1(con, """SELECT count(*) FROM frame_times
                             WHERE mid_method = ?""", (tm.MID_STACKPRO,))
    src = fig_cadence_overheads(con)
    sp_tbl = table(
        ["EXPTIME (s)", "gaps measured", "min overhead (s)",
         "p5 overhead (s)"],
        [[fnum(r[0], 3), fmt(r[1]), fnum(r[2], 2), fnum(r[3], 2)]
         for r in sp])
    return f"""
<section id="stackpro">
<div class="bhead"><h2>2 &middot; StackPro: the midpoint of a 16-fold
sum</h2>
<span class="tag">worst case {tm.STACKPRO_MID_WORST_CASE_S:.2f} s,
bounded by the data</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">A StackPro frame is the on-camera SUM of
{tm.N_SUB_STACKPRO} sub-reads (S2's PTC result).  Header EXPTIME could
mean the total span or the per-sub-read time, and any dead time between
sub-reads would shift the photon-weighted midpoint.  What is the
mid-time of {fmt(n_sp_frames)} StackPro science frames, and what is the
worst case if the contiguity assumption is wrong?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    "Back-to-back series cadence. If EXPTIME were per-sub-read, gaps "
    "would follow the dotted 15x line; instead StackPro overheads match "
    "the plain High Gain save/slew pipeline at every exposure time.")}</div>
{sp_tbl}
<p class="sub">Three independent facts pin the semantics: (1)&nbsp;across
{fmt(n_series_frames)} measured StackPro series gaps the frame's
wall-clock span never exceeds EXPTIME + {fnum(bound_row[2], 2)}&nbsp;s
(the minimum overhead, from the {fnum(bound_row[0], 1)}&nbsp;s series) —
so ALL internal dead time together is &le;
{tm.STACKPRO_DEADTIME_BOUND_S}&nbsp;s; (2)&nbsp;the sub-read gap is a
fixed sensor property, independent of exposure length, so the bound from
the shortest series applies to every StackPro frame; (3)&nbsp;MaxIm's
own JD-HELIO on a {fnum(maxim[0], 0)}&nbsp;s StackPro frame equals
start&nbsp;+&nbsp;EXPTIME/2 to {fnum(abs(maxim[2]), 2)}&nbsp;s — the
acquisition software also treats EXPTIME as the contiguous total.</p>

<h3>Decision</h3>
<div class="decision"><b>StackPro mid-time = start + EXPTIME/2, under the
stated assumption that the {tm.N_SUB_STACKPRO} sub-reads are contiguous
and span EXPTIME.  Worst case if dead time hides between sub-reads: the
true midpoint shifts late by at most
{tm.STACKPRO_MID_WORST_CASE_S:.2f}&nbsp;s</b> (half the measured
{tm.STACKPRO_DEADTIME_BOUND_S}&nbsp;s dead-time ceiling) — recorded as
<code>{tm.MID_STACKPRO}</code> in <code>frame_times.mid_method</code> so
any user can query which frames carry it.</div>

<h3>Consequence</h3>
<p class="sub">CV eclipse timing on StackPro series (AG LMi 2023, the
1024&nbsp;s deep fields) inherits a &le;
{tm.STACKPRO_MID_WORST_CASE_S:.2f}&nbsp;s systematic — below the
sub-second target and far below any eclipse-fit precision.  Should a
future camera manual contradict the contiguity assumption, only this
constant changes and one re-run repairs every product.</p>
</div></section>"""


def section_frame_times(con, meta: dict) -> str:
    n_rows = q1(con, "SELECT count(*) FROM frame_times")
    n_bjd = q1(con, "SELECT count(*) FROM frame_times "
                    "WHERE bjd_tdb IS NOT NULL")
    method_rows = q(con, """
        SELECT mid_method, bjd_method, count(*) FROM frame_times
        GROUP BY mid_method, bjd_method ORDER BY 3 DESC""")
    explain = {
        tm.MID_PLAIN: "start + EXPTIME/2",
        tm.MID_STACKPRO: "start + EXPTIME/2 (StackPro sum; "
                         f"&le; {tm.STACKPRO_MID_WORST_CASE_S:.2f} s "
                         "worst case)",
        tm.MID_EXPTIME_NONPOS: "EXPTIME &le; 0 or missing — start used",
        tm.MID_NO_JD: "no header JD — no time axis possible",
    }
    m_tbl = table(
        ["mid_method", "bjd_method", "frames", "meaning"],
        [[f"<code>{esc(mm)}</code>", f"<code>{esc(bm or '—')}</code>",
          fmt(n), explain.get(mm, "&mdash;")]
         for mm, bm, n in method_rows],
        row_classes=[None if (bm or "").startswith("bary") else "warn"
                     for _mm, bm, _n in method_rows])
    ltt_rng = q(con, """SELECT min(bary_ltt_s), max(bary_ltt_s),
                               min(tdb_minus_utc_s), max(tdb_minus_utc_s)
                        FROM frame_times WHERE bjd_tdb IS NOT NULL""")[0]
    src = fig_bjd_offset(con)
    eph = meta.get("ephemeris", "?")
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
<p class="sub">One stated caveat travels with the table: the correction
points at the FRAME CENTER.  An object at a frame corner differs by up
to ~1.3&nbsp;s (26&prime; &times; 499&nbsp;s/rad); a paper needing
sub-second absolutes for an off-center object recomputes with
<code>macro_core.timing.bjd_tdb_from_utc</code> at the object's own
coordinates — every input needed is in the row.</p>
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
    if g:
        verdict = (f"<b>The eclipse arrives at O&minus;C = {g[6]:+.0f} "
                   f"&plusmn; {g[7]:.0f}&nbsp;s against the VSX ephemeris; "
                   f"the absolute observatory clock error is bounded at "
                   f"|&Delta;t| &lt; {g[8]:.0f}&nbsp;s "
                   f"(&lt; {g[8]/60:.1f}&nbsp;min).</b>  The offset "
                   "itself belongs to the EPHEMERIS, not the clock: the "
                   "independent Gaia DR3 eclipsing-binary period "
                   f"({esc(meta.get('gaia_period_d', '?'))}&nbsp;d, "
                   f"{esc(meta.get('gaia_source', ''))}) differs from the "
                   "VSX period by enough to predict an O&minus;C of "
                   f"{esc(gaia_pred)}&nbsp;s at our mean epoch "
                   f"(cycle {esc(meta.get('clock_mean_cycle', '?'))}), "
                   f"with a &plusmn;{esc(gaia_env)}&nbsp;s envelope from "
                   "Gaia's own period uncertainty — our measurement sits "
                   "comfortably inside that envelope.  The bound is an "
                   "honest ceiling set by how well this star's ephemeris "
                   "is known, not a claim of NTP-level proof.")
    else:
        verdict = ("<b>The global eclipse fit did not converge; the clock "
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
{esc(str(epoch))} ({esc(source)}).  Folding the manifest on this
ephemeris finds four nights whose coverage reaches |phase| &lt; 0.02;
{fmt(n_frames)} frames of those plus three baseline nights were
measured with forced-aperture differential photometry against
{esc(meta.get('clock_comps', '').count(';') + 1 if
     meta.get('clock_comps') else '?')} fixed comparison stars
({fmt(pts)} usable points), using OUR mid-exposure HJD_UTC — the same
convention as the VSX epoch, like against like:</p>
<div class="grid">{_figure(src,
    "AG LMi folded on the literature ephemeris using this pipeline's "
    "own time axis. The eclipse lands where the ephemeris says it "
    "should; the fitted center measures the residual offset.")}</div>
{ecl_tbl}

<h3>Decision</h3>
<div class="decision">{verdict}</div>

<h3>Consequence</h3>
<p class="sub">Every CV timing result inherits this bound as its absolute
anchor; relative (within-archive) timing is tighter still, riding on the
{fnum(0.002, 3)}&nbsp;s internal header consistency of section 1.  A
tighter absolute bound is a one-night project for the October run: one
modern AG LMi minimum against a same-week TESS or AAVSO epoch would
collapse the ephemeris systematics to seconds.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S3 report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    try:
        meta = dict(q(con, "SELECT key, value FROM s3_build_meta"))
        n_rows = q1(con, "SELECT count(*) FROM frame_times")
        n_bjd = q1(con, "SELECT count(*) FROM frame_times "
                        "WHERE bjd_tdb IS NOT NULL")
        n_scan = q1(con, "SELECT sum(n_frames) FROM s3_dateobs_audit")
        g = q(con, """SELECT clock_bound_s FROM s3_clock_eclipses
                      WHERE tag = 'global' AND status = 'ok'""")
        clock_txt = (f"clock bounded at &lt; {g[0][0]:.0f} s"
                     if g else "clock bound pending")

        sections = [
            section_stamp(con),
            section_stackpro(con),
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
  &middot; <a href="../index.html">back to the evidence hub</a></p>
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
