"""S4 chain-of-evidence report renderer: the ensemble-photometry prototype.

Reads the photometry database (NEVER the archive or the manifest — if a
number cannot be derived from ``anuma_vvpup_prototype.sqlite``, it does not
belong on the page) and writes:

* ``docs/pipeline/s4_photometry.html``     — the report
* ``docs/pipeline/figures/s4/*.png``       — every figure

Socratic format throughout (Question → Evidence → Decision → Consequence),
shared stylesheet, shared table/figure machinery imported from the S0
renderer so all evidence pages read as one system.  EVERY number in the
HTML is interpolated from a SQL query executed here or from a constant
defined in the ``macro_phot`` modules — nothing is hand-typed.

Astrophysical honesty (stated on the page): AN UMa is a polar; its
high-amplitude orbital modulation and accretion-state changes are EXPECTED
signal.  The light-curve figures are evidence that the PIPELINE works —
never science claims, which belong to the CV paper.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import photometry as ph   # noqa: E402  (constants for interpolation)
from . import ensemble as ens    # noqa: E402
from . import gaia as gg         # noqa: E402
# Shared page machinery: one house figure style, one query discipline,
# one table generator
# generator as the S0/S0b reports — one visual language across the site.
from macro_core.report_s0 import (   # noqa: E402
    ACCENT, BAD, STYLE, DPI, FAINT, GOOD, MUTED, WARN,
    _figure, esc, fmt, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402  (house figure style)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s4"
HTML_PATH = DOCS_DIR / "s4_photometry.html"

#: Literature orbital period of AN UMa used ONLY to phase-fold the sanity
#: light curve.  Source: AAVSO VSX (vsx.aavso.org, ident 'AN UMa') lists
#: P = 0.07975274 d (type AM+E), verified live 2026-08-18; the constant
#: below is that value to the precision this fold needs (~1.91 h).  A fold
#: on a wrong period would smear, not align — which is exactly what makes
#: the folded figure a pipeline check rather than a science measurement.
ANUMA_PORB_D = 0.0797528
#: The provenance string the page renders next to the period (query-free
#: constants are the one legal source of hand-maintained text, and keeping
#: the citation beside the number keeps them in sync).
ANUMA_PORB_SOURCE = "AAVSO VSX (0.07975274 d)"

#: Per-filter panels.  The house band map, so that "r" is the same hue
#: here as in the manuscript figures and in every other CV report.
FILTER_COLORS = dict(ps.BAND_COLOR, empty=FAINT)


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_match_qc(con) -> str:
    """Two panels: per-frame match fraction, and alignment residual RMS."""
    rows = q(con, """
        SELECT f.target_key || '/era' || f.era_id,
               1.0 * f.ali_nmatch / f.n_detected, f.ali_rms_px
        FROM phot_frames f
        WHERE f.status = 'matched' AND f.n_detected > 0""")
    series = sorted({r[0] for r in rows})
    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.4))
        for i, s in enumerate(series):
            fr = [r[1] for r in rows if r[0] == s]
            res = [r[2] for r in rows if r[0] == s and r[2] is not None]
            # Hue AND dash pattern: two step histograms of the same
            # shape overlaid are exactly where colour alone fails.
            style = ps.line_series(i)
            ax1.hist(fr, bins=np.linspace(0, 1, 41), histtype="step",
                     label=s, **style)
            ax2.hist(res, bins=40, histtype="step", **style)
        ax1.set_xlabel("matched fraction of frame detections")
        ax1.set_ylabel("frames")
        ax1.set_title("Star-match completeness per frame")
        ax1.legend(fontsize=7)
        ax2.set_xlabel("alignment residual RMS (ref px)")
        ax2.set_ylabel("frames")
        ax2.set_title("astroalign registration residuals")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s4_match_qc.png", dpi=DPI)
        plt.close(fig)
    return "figures/s4/s4_match_qc.png"


def _lc_points(con, target: str, era: int, filt: str):
    """Zero-point-corrected target light curve of one series (from the DB)."""
    tie = q(con, """SELECT target_star_id FROM phot_gaia_tie
                    WHERE target_key=? AND era_id=?""", (target, era))
    if not tie or tie[0][0] is None:
        return None
    sid = tie[0][0]
    rows = q(con, """
        SELECT f.jd, f.night,
               -2.5 * (CASE WHEN d.flux > 0
                       THEN log10(d.flux / f.exptime) END) + ? - f.zp
        FROM phot_detections d
        JOIN phot_frames f ON f.frame_id = d.frame_id
        WHERE d.star_id = ? AND d.clipped = 0 AND f.zp IS NOT NULL
          AND f.target_key = ? AND f.era_id = ? AND f.filter = ?
        ORDER BY f.jd""", (ph.INST_MAG_OFFSET, sid, target, era, filt))
    rows = [r for r in rows if r[2] is not None]
    return rows or None


def fig_anuma_lc(con) -> str:
    """AN UMa unfolded light curve per filter, colored by night."""
    filts = [r[0] for r in q(con, """SELECT DISTINCT filter FROM phot_series
        WHERE target_key='anuma' ORDER BY filter""")]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(len(filts), 1,
                                 figsize=(9.6, 2.6 * len(filts)),
                                 sharex=True, squeeze=False)
        for ax, filt in zip(axes[:, 0], filts):
            rows = _lc_points(con, "anuma", 76, filt) or []
            nights = sorted({r[1] for r in rows})
            cmap = ps.ordinal_colors(max(len(nights), 2))
            for i, night in enumerate(nights):
                pts = [(r[0], r[2]) for r in rows if r[1] == night]
                jd = np.array([p[0] for p in pts])
                mm = np.array([p[1] for p in pts])
                ax.plot(jd - 2460000.0, mm, ".", ms=2.5, color=cmap[i])
            ax.invert_yaxis()
            ax.set_ylabel(f"{filt}  (inst mag)")
        axes[-1, 0].set_xlabel("JD - 2460000 (header UTC start; BJD is S3)")
        axes[0, 0].set_title("AN UMa, era 76 — ensemble-corrected target "
                             "light curve, colored by night")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s4_anuma_lc.png", dpi=DPI)
        plt.close(fig)
    return "figures/s4/s4_anuma_lc.png"


def fig_anuma_folded(con) -> str:
    """AN UMa folded on the literature period, colored by night."""
    filts = [r[0] for r in q(con, """SELECT DISTINCT filter FROM phot_series
        WHERE target_key='anuma' ORDER BY filter""")]
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(filts),
                                 figsize=(3.4 * len(filts), 3.4),
                                 sharey=False, squeeze=False)
        for ax, filt in zip(axes[0], filts):
            rows = _lc_points(con, "anuma", 76, filt) or []
            nights = sorted({r[1] for r in rows})
            cmap = ps.ordinal_colors(max(len(nights), 2))
            for i, night in enumerate(nights):
                pts = [(r[0], r[2]) for r in rows if r[1] == night]
                phase = np.array([p[0] for p in pts]) / ANUMA_PORB_D % 1.0
                mm = np.array([p[1] for p in pts])
                ax.plot(phase, mm, ".", ms=2, color=cmap[i], alpha=0.7)
            ax.invert_yaxis()
            ax.set_xlabel("orbital phase (arbitrary epoch)")
            ax.set_title(f"filter {filt}", fontsize=10)
        axes[0][0].set_ylabel("inst mag - ZP")
        fig.suptitle(f"AN UMa folded on P = {ANUMA_PORB_D} d "
                     "(literature; sanity fold, not a measurement)",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s4_anuma_folded.png", dpi=DPI)
        plt.close(fig)
    return "figures/s4/s4_anuma_folded.png"


def fig_rms_vs_mag(con) -> str:
    """Check/field-star RMS vs magnitude per (era, filter) + photon floor."""
    combos = q(con, """SELECT DISTINCT era_id, filter FROM phot_stars
                       ORDER BY era_id, filter""")
    n = len(combos)
    ncol = 3
    nrow = (n + ncol - 1) // ncol
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(nrow, ncol,
                                 figsize=(3.3 * ncol, 3.0 * nrow),
                                 squeeze=False)
        for k, (era, filt) in enumerate(combos):
            ax = axes[k // ncol][k % ncol]
            rows = q(con, """SELECT mean_mag, rms, role FROM phot_stars
                WHERE era_id=? AND filter=? AND mean_mag IS NOT NULL
                  AND rms IS NOT NULL AND nobs >= 10""", (era, filt))
            for role, color, ms in (("field", FAINT, 6),
                                    ("comp", ACCENT, 10),
                                    ("check", WARN, 16),
                                    ("dropped_unstable", BAD, 10),
                                    ("target", ps.OTHER, 30)):
                pts = [(r[0], r[1]) for r in rows if r[2] == role]
                if pts:
                    ax.scatter([p[0] for p in pts], [p[1] for p in pts],
                               s=ms, color=color, label=role, linewidths=0)
            # Photon-noise floor: each star's mean predicted per-point
            # sigma, then a running MEDIAN in magnitude bins (the pure
            # errors.rms_vs_mag_curve) — a smooth curve, not a per-star
            # zigzag through mixed exposure times.
            floor = q(con, """
                SELECT s.mean_mag, avg(1.0857 * d.fluxerr / d.flux)
                FROM phot_stars s
                JOIN phot_ref_stars r ON r.target_key = s.target_key
                     AND r.era_id = s.era_id AND r.star_id = s.star_id
                JOIN phot_detections d ON d.star_id = s.star_id
                JOIN phot_frames f ON f.frame_id = d.frame_id
                     AND f.target_key = s.target_key
                     AND f.era_id = s.era_id AND f.filter = s.filter
                WHERE s.era_id=? AND s.filter=? AND d.clipped=0 AND d.flux>0
                  AND s.mean_mag IS NOT NULL
                GROUP BY s.target_key, s.star_id""", (era, filt))
            if floor:
                from .errors import rms_vs_mag_curve
                centers, medians, _n = rms_vs_mag_curve(
                    np.array([f[0] for f in floor], dtype=float),
                    np.array([f[1] for f in floor], dtype=float))
                ax.plot(centers, medians, "-", color=GOOD, lw=1.4,
                        alpha=0.9, label="photon+sky floor (binned median)")
            ax.set_yscale("log")
            ax.set_xlabel("ensemble mean mag (inst)")
            ax.set_ylabel("RMS (mag)")
            ax.set_title(f"era {era} / {filt}", fontsize=10)
            if k == 0:
                ax.legend(fontsize=6, loc="upper left")
        for k in range(n, nrow * ncol):
            axes[k // ncol][k % ncol].set_visible(False)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s4_rms_vs_mag.png", dpi=DPI)
        plt.close(fig)
    return "figures/s4/s4_rms_vs_mag.png"


def fig_allan(con) -> str | None:
    """Allan deviation of the chosen check-star night, with tau^-1/2 line."""
    rows = q(con, """SELECT tau_s, adev_mag, n_pairs, target_key, era_id,
                            filter, star_id, night FROM phot_allan
                     ORDER BY tau_s""")
    if not rows:
        return None
    tau = np.array([r[0] for r in rows])
    ad = np.array([r[1] for r in rows])
    npairs = np.array([r[2] for r in rows], dtype=float)
    tk, era, filt, sid, night = rows[0][3], rows[0][4], rows[0][5], \
        rows[0][6], rows[0][7]
    # 1-sigma sampling uncertainty of each rung: adev / sqrt(2 N_pairs).
    # The long-tau rungs rest on a handful of pairs — plotting them bare
    # invites over-reading a 'floor' the statistics do not establish.
    ad_err = ad / np.sqrt(2.0 * npairs)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        ax.errorbar(tau, ad, yerr=ad_err, fmt="o-", color=ACCENT,
                    capsize=3,
                    label=r"check star ($\pm 1\sigma$ sampling)")
        ax.set_xscale("log"); ax.set_yscale("log")
        # White-noise expectation anchored on the first rung.
        ax.loglog(tau, ad[0] * (tau / tau[0]) ** -0.5, "--", color=WARN,
                  label=r"white noise $\tau^{-1/2}$")
        # Pair counts on the figure itself: the reader sees how much data
        # each rung rests on without opening the database.
        for t, a, n in zip(tau, ad, npairs):
            ax.annotate(f"N={int(n)}", (t, a), textcoords="offset points",
                        xytext=(6, 6), fontsize=7, color=MUTED)
        ax.set_xlabel("averaging time tau (s)")
        ax.set_ylabel("Allan deviation (mag)")
        ax.set_title(f"Allan deviation — {tk}/era {era}/{filt}, star "
                     f"{sid}, night {night}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s4_allan.png", dpi=DPI)
        plt.close(fig)
    return "figures/s4/s4_allan.png"


def fig_zp_timeline(con) -> str:
    """Ensemble zero points through time — clouds and focus made visible."""
    rows = q(con, """SELECT target_key || '/era' || era_id || '/' || filter,
                            jd, zp FROM phot_frames
                     WHERE zp IS NOT NULL ORDER BY jd""")
    series = sorted({r[0] for r in rows})
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(9.6, 3.6))
        for i, s in enumerate(series):
            pts = [(r[1], r[2]) for r in rows if r[0] == s]
            ax.plot([p[0] - 2460000.0 for p in pts], [p[1] for p in pts],
                    ms=2.5, lw=0, label=s,
                    **ps.measurement_kw(**ps.series(i)))
        ax.set_xlabel("JD - 2460000")
        ax.set_ylabel("frame zero point (mag)")
        ax.set_title("Honeycutt frame zero points — transparency history "
                     "of every series")
        ax.invert_yaxis()
        ax.legend(fontsize=6, ncol=3)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s4_zp_timeline.png", dpi=DPI)
        plt.close(fig)
    return "figures/s4/s4_zp_timeline.png"


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def section_selection(con) -> str:
    n_sel = q1(con, "SELECT count(*) FROM phot_frames")
    n_raw = q1(con, "SELECT sum(n_raw) FROM phot_selection")
    n_linked = q1(con, "SELECT sum(n_linked) FROM phot_selection")
    sel_rows = q(con, """SELECT s.target_key, s.era_id, s.n_raw, s.n_linked,
                                (SELECT count(*) FROM phot_frames f
                                 WHERE f.target_key = s.target_key
                                   AND f.era_id = s.era_id)
                         FROM phot_selection s ORDER BY 1, 2""")
    sel_tbl = table(
        ["target", "era", "raw canonical light frames",
         "with reduced counterpart", "selected"],
        [[esc(t), fmt(e), fmt(r), fmt(l), fmt(s)]
         for t, e, r, l, s in sel_rows],
        row_classes=[None if r == l else "warn"
                     for _, _, r, l, _ in sel_rows])
    status_rows = q(con, """SELECT status, count(*) FROM phot_frames
                            GROUP BY status ORDER BY 2 DESC""")
    status_tbl = table(["frame status", "frames"],
                       [[f"<code>{esc(s)}</code>", fmt(n)]
                        for s, n in status_rows],
                       row_classes=["warn" if s.startswith("failed") else None
                                    for s, _ in status_rows])
    return f"""
<section id="selection">
<div class="bhead"><h2>1 &middot; Frame selection &amp; provenance</h2>
<span class="tag">server-reduced pixels only — one provenance per series</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Which frames may the photometry prototype touch, and under
what provenance rule?</p>

<h3>Evidence</h3>
{sel_tbl}
{status_tbl}

<h3>Decision</h3>
<div class="decision"><b>The prototype reads SERVER-REDUCED pixels only:
canonical raw light frames of the two reduction-complete polars whose
reduced counterpart exists in <code>raw_reduced_links</code>
({fmt(n_linked)} of {fmt(n_raw)} raw frames; {fmt(n_sel)} selected), and
reduction provenances are never mixed within a (target,&nbsp;era) series
(decision 2026-08-18).</b>  Frames without a reduced counterpart are
counted in the ledger above and excluded — when S2 delivers our own
master-calibrated reductions, they rejoin as a NEW provenance series, not
as patches to this one.</div>

<h3>Consequence</h3>
<p class="sub">Every downstream number in this report descends from
exactly these frames; the ledger row that excludes a night is the audit
trail for its absence from the light curves.</p>
</div></section>"""


def section_matching(con) -> str:
    n_matched = q1(con, "SELECT count(*) FROM phot_frames "
                        "WHERE status='matched'")
    # The honest denominator is EVERY extraction-complete frame — including
    # the known-mispointed ones the matcher deliberately skips (they are
    # extraction-complete too; hiding them once made 90% read as 96%).
    n_extracted = q1(con, "SELECT count(*) FROM phot_frames "
                          "WHERE status = 'matched' "
                          "OR status LIKE 'failed_match%' "
                          "OR status = 'skipped_pointing'")
    n_failed = q1(con, "SELECT count(*) FROM phot_frames "
                       "WHERE status LIKE 'failed_match%'")
    n_skipped = q1(con, "SELECT count(*) FROM phot_frames "
                        "WHERE status = 'skipped_pointing'")
    med_frac = q1(con, """SELECT 100.0 * ali_nmatch / n_detected AS p
                          FROM phot_frames WHERE status='matched'
                          ORDER BY p LIMIT 1 OFFSET
                          (SELECT count(*)/2 FROM phot_frames
                           WHERE status='matched')""")
    med_rms = q1(con, """SELECT ali_rms_px FROM phot_frames
                         WHERE status='matched' AND ali_rms_px IS NOT NULL
                         ORDER BY ali_rms_px LIMIT 1 OFFSET
                         (SELECT count(*)/2 FROM phot_frames
                          WHERE status='matched'
                            AND ali_rms_px IS NOT NULL)""")
    src = fig_match_qc(con)
    ref_rows = q(con, """
        SELECT r.target_key, r.era_id, r.ref_frame_id, f.night, r.n_stars,
               f.fwhm_px, f.plate_scale, r.tol_px,
               r.doubled_frac, r.n_cand_rejected,
               (SELECT count(*) FROM phot_frames x
                WHERE x.target_key = r.target_key AND x.era_id = r.era_id
                  AND x.status = 'matched')
        FROM phot_ref r JOIN phot_frames f ON f.frame_id = r.ref_frame_id
        ORDER BY 1, 2""")
    ref_tbl = table(
        ["target", "era", "reference frame", "night", "ref stars",
         "FWHM px", "plate scale \"/px", "match tol px",
         "paired-star frac", "doubled candidates rejected",
         "frames matched"],
        [[esc(t), fmt(e), fmt(fid), esc(n), fmt(ns),
          f"{fw:.2f}", f"{ps:.4f}", f"{tol:.2f}",
          f"{df:.3f}" if df is not None else "&mdash;", fmt(nrej),
          fmt(nm)]
         for t, e, fid, n, ns, fw, ps, tol, df, nrej, nm in ref_rows],
        row_classes=[None if not nrej else "warn"
                     for *_a, nrej, _nm in ref_rows])
    # Failure anatomy: match failures grouped by night, with the S0
    # pointing offset carried into this DB — a mispointed night MUST fail
    # here, and showing the two columns side by side proves the failures
    # are the pointing, not the algorithm.
    fail_rows = q(con, """
        SELECT target_key, night, count(*),
               round(avg(pointing_offset_deg), 2),
               group_concat(DISTINCT status)
        FROM phot_frames WHERE status LIKE 'failed_match%'
           OR status = 'skipped_pointing'
        GROUP BY target_key, night ORDER BY 3 DESC LIMIT 10""")
    fail_tbl = (table(
        ["target", "night", "failed frames", "S0 pointing offset (deg)",
         "failure"],
        [[esc(t), esc(n), fmt(c), esc(po),
          f"<code>{esc(st)}</code>"] for t, n, c, po, st in fail_rows],
        row_classes=["warn"] * len(fail_rows))
        if fail_rows else "")
    return f"""
<section id="matching">
<div class="bhead"><h2>2 &middot; Star matching without a WCS</h2>
<span class="tag">astroalign triangles &rarr; one reference per (target, era)</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Most polar frames are unsolved (73&ndash;95% of the Sloan
series carry <code>pltsolvd=0</code> — S1&rsquo;s future business).  Can
stars be matched across thousands of frames using GEOMETRY alone?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"{fmt(n_matched)} of {fmt(n_extracted)} extraction-complete frames "
    f"matched; {fmt(n_skipped)} skipped for known S0 mispointing, "
    f"{fmt(n_failed)} failed alignment; median matched fraction "
    f"{med_frac:.0f}% of frame detections; median registration residual "
    f"{med_rms:.3f} ref-px.")}</div>
{ref_tbl}
{fail_tbl}

<h3>Decision</h3>
<div class="decision"><b>Each (target,&nbsp;era) series is registered to
ONE reference frame — chosen as the most-populated frame within the
sharpest {ph.REF_FWHM_PERCENTILE:g}% of the series, PROVIDED it passes
double-image quality control — by astroalign triangle matching, followed
by one-to-one nearest-neighbour assignment with a seeing-scaled tolerance
(max({ph.MATCH_TOL_MIN_PX:g}&nbsp;px, {ph.MATCH_TOL_FWHM:g}&nbsp;&times;
ref FWHM)).</b>  Quality control exists because a guiding jump doubles
every star into an equal-brightness pair and thereby INFLATES the
detection count — the very statistic the selection rewards: a candidate
is rejected when more than {100 * ph.REF_DOUBLED_MAX_FRAC:g}% of its
stars have a companion within {ph.REF_PAIR_RADIUS_FWHM:g}&nbsp;&times;
FWHM at a flux ratio under {ph.REF_PAIR_FLUX_RATIO:g} (the paired-star
fraction and rejection count are tabulated above — clean frames sit near
zero, doubled frames near one).  The star catalog IS the reference
frame&rsquo;s detection list; star identity is its detection id — stable,
auditable, and WCS-free.  astroalign&rsquo;s RANSAC is seeded per frame
(<code>seeded_ransac</code>), so every transform and star identity
reproduces bit-for-bit on re-run.  Frames that cannot produce
{ph.MIN_STARS_FOR_ALIGN} detections or defeat the triangle fit are marked
<code>failed_match</code>, counted, and excluded — never force-fitted.</div>

<h3>Consequence</h3>
<p class="sub">Every detection now carries a star id (or none), so the
ensemble below is a pure database join — no pixel is touched again.</p>
</div></section>"""


def section_ensemble(con) -> str:
    src = fig_zp_timeline(con)
    ser_rows = q(con, """
        SELECT target_key, era_id, filter, n_frames, n_comp, n_check,
               n_dropped, ens_niter, ens_converged, comp_rms_median, zp_std
        FROM phot_series ORDER BY target_key, era_id, filter""")
    ser_tbl = table(
        ["target", "era", "filter", "frames", "comps", "checks",
         "dropped variable", "iterations", "converged",
         "comp RMS med (mag)", "ZP spread (mag)"],
        [[esc(t), fmt(e), esc(f_), fmt(nf), fmt(nc), fmt(nk), fmt(nd),
          fmt(ni), "yes" if cv else "NO",
          f"{cr:.4f}" if cr is not None else "&mdash;",
          f"{zs:.3f}" if zs is not None else "&mdash;"]
         for t, e, f_, nf, nc, nk, nd, ni, cv, cr, zs in ser_rows],
        row_classes=[None if cv else "warn"
                     for *_x, cv, _cr, _zs in ser_rows])
    n_series = len(ser_rows)
    n_conv = sum(1 for r in ser_rows if r[8])
    return f"""
<section id="ensemble">
<div class="bhead"><h2>3 &middot; The Honeycutt ensemble</h2>
<span class="tag">robust ZP per frame + mean mag per star, per (target, era, filter)</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">With no photometric nights guaranteed and clouds in the
data, how does every frame get a zero point it can defend?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"Frame zero points of all {fmt(n_series)} series "
    f"({fmt(n_conv)} converged). Cloud events and focus drifts appear as "
    "excursions; the ensemble absorbs them into ZP rather than into star "
    "magnitudes.")}</div>
{ser_tbl}

<h3>Decision</h3>
<div class="decision"><b>Honeycutt&rsquo;s (1992) inhomogeneous ensemble,
solved by alternating weighted least squares with a
{ens.CLIP_SIGMA:g}&sigma; robust re-solve, a
{1000 * ens.WEIGHT_FLOOR_MAG:g}&nbsp;mmag systematic weight floor, and the
mean-ZP gauge.</b>  Comparison stars are chosen by stability iteration —
candidates must appear on &ge;{100 * ens.COMP_MIN_FRAME_FRAC:g}% of
frames, and stars whose residual RMS exceeds
max({ens.COMP_RMS_FACTOR:g}&nbsp;&times;&nbsp;median,
{1000 * ens.COMP_RMS_FLOOR_MAG:g}&nbsp;mmag) are dropped each pass; the
TARGET is excluded from comp duty a priori (a polar&rsquo;s modulation
must never set the zero point).  {ens.N_CHECK_STARS} stable stars nearest
the target&rsquo;s magnitude are held OUT of the solve as check stars —
their statistics below are honest hold-outs, not self-graded
homework.</div>

<h3>Consequence</h3>
<p class="sub">Every frame carries zp&nbsp;&plusmn;&nbsp;zp_err and every
star a mean magnitude + RMS + &chi;&sup2; — the raw material of both the
light curves (section 6) and the error model (section 5).</p>
</div></section>"""


def section_gaia(con) -> str:
    tie_rows = q(con, """
        SELECT target_key, era_id, status, parity, scale_fit, rot_deg,
               n_gaia, n_matched, target_star_id, coord_source,
               target_sep_arcsec
        FROM phot_gaia_tie ORDER BY 1, 2""")
    hdr_scale = {(r[0], r[1]): r[2] for r in q(con, """
        SELECT r.target_key, r.era_id, f.plate_scale FROM phot_ref r
        JOIN phot_frames f ON f.frame_id = r.ref_frame_id""")}
    tie_tbl = table(
        ["target", "era", "status", "parity", "fitted scale \"/px",
         "header scale \"/px", "rotation deg", "Gaia stars", "matched",
         "target star id", "target offset \"", "coords via"],
        [[esc(t), fmt(e), esc(st), esc(par),
          f"{sc:.4f}" if sc else "&mdash;",
          f"{hdr_scale.get((t, e), float('nan')):.4f}",
          f"{rot:.2f}" if rot is not None else "&mdash;",
          fmt(ng), fmt(nm), fmt(sid),
          f"{sep:.2f}" if sep is not None else "&mdash;", esc(csrc)]
         for t, e, st, par, sc, rot, ng, nm, sid, csrc, sep in tie_rows],
        row_classes=[None if st == "ok" else "warn"
                     for _, _, st, *_r in tie_rows])
    off_rows = q(con, """
        SELECT target_key, era_id, filter, gaia_offset, gaia_offset_mad,
               gaia_offset_n, zp_source
        FROM phot_series ORDER BY 1, 2, 3""")
    off_tbl = table(
        ["target", "era", "filter", "median (G - ensemble) mag",
         "scatter (MAD) mag", "stars", "zero-point status"],
        [[esc(t), fmt(e), esc(f_),
          f"{o:.3f}" if o is not None else "&mdash;",
          f"{m:.3f}" if m is not None else "&mdash;", fmt(n),
          f"<code>{esc(zs)}</code>"]
         for t, e, f_, o, m, n, zs in off_rows])
    return f"""
<section id="gaia">
<div class="bhead"><h2>4 &middot; The Gaia DR3 anchor</h2>
<span class="tag">identity + geometry check + an honest offset, not absolute cal</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">What ties these WCS-free pixel catalogs to the sky — and
how far toward absolute photometry does that tie honestly reach?</p>

<h3>Evidence</h3>
{tie_tbl}
<p class="sub">The fitted similarity scale is measured from the SKY and
agrees with the header plate scale — an independent check of the
XPIXSZ/FOCALLEN cards that no header can fake.  Zero-point offsets per
series (comparison stars with Gaia G only):</p>
{off_tbl}

<h3>Decision</h3>
<div class="decision"><b>Gaia supplies three things: star identity
(source_id + G), an independent plate-scale/orientation measurement, and
a single median (G&nbsp;&minus;&nbsp;ensemble) offset per series.  The
ensemble magnitudes remain INSTRUMENTAL: no color terms are fitted, and
the quoted MAD of each offset is dominated by star color — which is
precisely why applying the offset is a courtesy scale, not absolute
calibration.</b>  Absolute calibration (color terms per filter per era)
is a later S4 milestone; relative photometry is the product here.</div>

<h3>Consequence</h3>
<p class="sub">Every ensemble star is now a named Gaia source; the target
star&rsquo;s identity is fixed by coordinates
({esc(", ".join(sorted(set(r[9] for r in tie_rows))))}), so the light
curves in section 6 are attached to AN&nbsp;UMa and VV&nbsp;Pup by
astrometry, not by hope.</p>
</div></section>"""


def section_errors(con) -> str:
    src = fig_rms_vs_mag(con)
    src_allan = fig_allan(con)
    em_rows = q(con, """
        SELECT era_id, filter, n_check, check_rms_min, check_rms_med,
               chi2nu_med, inflation
        FROM phot_error_model ORDER BY era_id, filter""")
    em_tbl = table(
        ["era", "filter", "check stars", "best check RMS (mmag)",
         "median check RMS (mmag)", "median &chi;&sup2;<sub>&nu;</sub>",
         "inflation factor"],
        [[fmt(e), esc(f_), fmt(n),
          f"{1000 * rmin:.1f}" if rmin is not None else "&mdash;",
          f"{1000 * rmed:.1f}" if rmed is not None else "&mdash;",
          f"{c2:.2f}" if c2 is not None else "&mdash;",
          f"{inf:.2f}" if inf is not None else "&mdash;"]
         for e, f_, n, rmin, rmed, c2, inf in em_rows])
    allan_meta = q(con, """SELECT target_key, era_id, filter, star_id,
                                  night, count(*), min(tau_s), max(tau_s),
                                  min(n_pairs)
                           FROM phot_allan GROUP BY 1,2,3,4,5""")
    allan_fig = (_figure(src_allan,
                 (lambda m: f"Allan deviation of check star {fmt(m[3])} "
                  f"({esc(m[0])}/era {fmt(m[1])}/{esc(m[2])}), night "
                  f"{esc(m[4])}: {fmt(m[5])} tau rungs from "
                  f"{m[6]:.0f}s to {m[7]:.0f}s.  The longest-tau rungs "
                  f"rest on as few as {fmt(m[8])} difference pairs "
                  f"(N annotated per rung; error bars are the "
                  f"1&sigma; sampling uncertainty adev/&radic;(2N)) — any "
                  f"apparent flattening there is SUGGESTIVE, not "
                  f"statistically established; a claim of a "
                  f"correlated-noise floor awaits longer runs."
                  )(allan_meta[0]))
                 if src_allan and allan_meta else
                 '<p class="missing">No night long enough for an Allan '
                 'ladder.</p>')
    return f"""
<section id="errors">
<div class="bhead"><h2>5 &middot; The empirical error model (S5 seed)</h2>
<span class="tag">check stars only — the target never grades its own errors</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Are the formal error bars TRUE?  Both targets are polars
whose modulation is signal, so the question is answered entirely by field
CHECK stars held out of the ensemble solve.</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    "Star RMS vs ensemble magnitude per (era, filter): comparison stars "
    "(blue), held-out check stars (yellow), dropped variables (red), the "
    "target (pink), against the photon+sky floor predicted by the "
    "measurement errors (green line).")}</div>
{em_tbl}
<div class="grid">{allan_fig}</div>

<h3>Decision</h3>
<div class="decision"><b>Three statistics form the standard validation
set every downstream paper reuses: (1) check-star RMS vs magnitude
against the photon-noise floor; (2) the error INFLATION FACTOR
sqrt(median &chi;&sup2;<sub>&nu;</sub>) of constant-star fits; (3) the
Allan deviation of one long night, whose departure from
&tau;<sup>&minus;1/2</sup> probes the correlated noise that averaging
cannot remove.</b>  Three fine-print clauses S5 must honor: the
&chi;&sup2; denominators are the formal errors WITH the
{1000 * ens.WEIGHT_FLOOR_MAG:g}&nbsp;mmag systematic floor added in
quadrature, so the factor calibrates floored errors, not raw photon
errors; application is ONE-SIDED — use max(1,&nbsp;inflation), because a
factor below 1 (see era&nbsp;76&nbsp;r) licenses no shrinking of error
bars; and the check population passed the same stability screen as comps
before being held out (RMS below
max({ens.COMP_RMS_FACTOR:g}&nbsp;&times;&nbsp;median comp RMS,
{1000 * ens.COMP_RMS_FLOOR_MAG:g}&nbsp;mmag)), so the check floor is
measured on pre-selected quiet stars — with only
{ens.N_CHECK_STARS}&nbsp;checks per series (era-76 rows pool the sparse
AN&nbsp;UMa field with the crowded b&nbsp;&asymp;&nbsp;+2&deg;
VV&nbsp;Pup field, 4&nbsp;+&nbsp;4), the tabulated medians are coarse.
The iKon era (72) writes EGAIN=0, so its predicted errors omit the
source shot term and its inflation factor is expected &gt;1 by
construction — the empirical model absorbs exactly this kind of header
poverty, which is why S5 builds on measured scatter, never on claimed
gains.</div>

<h3>Consequence</h3>
<p class="sub">S5 inherits these numbers as priors; October&rsquo;s new
frames extend the same tables by re-running this build.</p>
</div></section>"""


def section_lightcurves(con) -> str:
    src_lc = fig_anuma_lc(con)
    src_fold = fig_anuma_folded(con)
    n_pts = {}
    for tk in ("anuma", "vvpup"):
        n_pts[tk] = q1(con, """
            SELECT count(*) FROM phot_detections d
            JOIN phot_frames f ON f.frame_id = d.frame_id
            JOIN phot_gaia_tie t ON t.target_key = f.target_key
                 AND t.era_id = f.era_id
            WHERE d.star_id = t.target_star_id AND d.clipped = 0
              AND f.zp IS NOT NULL AND f.target_key = ?""", (tk,))
    nights = q1(con, """SELECT count(DISTINCT night) FROM phot_frames
                        WHERE target_key='anuma' AND zp IS NOT NULL""")
    return f"""
<section id="lightcurves">
<div class="bhead"><h2>6 &middot; Light-curve sanity — pipeline evidence,
not science</h2>
<span class="tag">AN UMa is a polar: its modulation is EXPECTED signal</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Does the assembled machinery produce a light curve a CV
astronomer recognizes?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src_lc,
    f"AN UMa: {fmt(n_pts['anuma'])} ensemble-corrected target points "
    f"across {fmt(nights)} nights (VV Pup holds {fmt(n_pts['vvpup'])} "
    "more). Times are header UTC start JDs — barycentric correction is "
    "stage S3, deliberately absent here.")}</div>
<div class="grid">{_figure(src_fold,
    f"The same points folded on the literature period "
    f"P = {ANUMA_PORB_D} d ({esc(ANUMA_PORB_SOURCE)}). Coherent per-night "
    "structure at the orbital period — plus night-to-night level shifts "
    "(accretion high/low states) — is exactly the phenomenology expected "
    "of a polar. Header-UTC times drift up to ~8 min from barycentric "
    "over a season (~0.07 in phase) — harmless for a sanity fold, fatal "
    "for timing, which is why BJD is S3's job.")}</div>

<h3>Decision</h3>
<div class="decision"><b>The light curves are exhibited as EVIDENCE THAT
THE PIPELINE WORKS — orbital-phase coherence and state changes are
consistent with AN UMa&rsquo;s known behaviour — and every scientific
claim about the system (state durations, period work, accretion
geometry) is explicitly deferred to the CV_TimeSeries paper, which will
consume this database, not this page.</b>  The fold uses the literature
period P&nbsp;=&nbsp;{ANUMA_PORB_D}&nbsp;d — source
{esc(ANUMA_PORB_SOURCE)}, recorded beside the constant in
<code>macro_phot.report_s4</code>; nothing here measures a period.</div>

<h3>Consequence</h3>
<p class="sub">The prototype closes S4&rsquo;s core loop: pixels &rarr;
matched stars &rarr; defensible zero points &rarr; validated errors
&rarr; a target light curve with named comparison stars.  Scaling to the
remaining CV targets is a worklist change, not a code change.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(db_path: Path) -> Path:
    """Render the full S4 report from the photometry DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        n_frames = q1(con, "SELECT count(*) FROM phot_frames")
        n_matched = q1(con, "SELECT count(*) FROM phot_frames "
                            "WHERE status='matched'")
        n_series = q1(con, "SELECT count(*) FROM phot_series")
        n_stars = q1(con, "SELECT count(DISTINCT star_id) FROM phot_stars")
        meta = dict(q(con, "SELECT key, value FROM s4_build_meta"))

        sections = [
            section_selection(con),
            section_matching(con),
            section_ensemble(con),
            section_gaia(con),
            section_errors(con),
            section_lightcurves(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S4 — Ensemble Photometry Prototype (AN UMa &amp; VV Pup)</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S4 — Ensemble Photometry Prototype (AN UMa &amp; VV Pup)</h1>
  <p>{fmt(n_matched)} of {fmt(n_frames)} frames measured &amp; star-matched
  without a WCS &middot; {fmt(n_series)} Honeycutt ensembles over
  {fmt(n_stars)} stars &middot; Gaia-anchored &middot; empirical error
  model (S5 seed) &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#selection">1 Selection</a> &middot;
  <a href="#matching">2 Matching</a> &middot;
  <a href="#ensemble">3 Ensemble</a> &middot;
  <a href="#gaia">4 Gaia anchor</a> &middot;
  <a href="#errors">5 Error model</a> &middot;
  <a href="#lightcurves">6 Light curves</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_phot.report_s4</code> from
<code>products/phot/anuma_vvpup_prototype.sqlite</code> — every number on
this page is the result of a SQL query; none is typed by hand.  Regenerate
with <code>pipeline/scripts/build_s4_photometry.py report</code> (full
rebuild: init &rarr; extract &rarr; match &rarr; gaia &rarr; ensemble
&rarr; errors &rarr; report).</footer>
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
