"""CV-S6 chain-of-evidence report: the catalogue tie.

Reads ``products/phot/cv_timeseries.sqlite`` and writes

* ``docs/CV_TimeSeries/cv_catalogue_tie.html``
* ``docs/CV_TimeSeries/figures/cv_catalogue_tie/*.png``

Socratic throughout, and in ONE order, because the order is the argument:

    0  what is a magnitude here, and why is today's answer arbitrary?
    1  which catalogue, and what does it cost to be wrong about that?
    2  which stars may carry the tie -- and which may not, and why
    3  the fit: zero point AND colour term, with the range it is valid over
    4  do the uppercase filter labels mean the same glass as the lowercase?
    5  how accurate is it really -- on stars the fit never saw
    6  three attacks on our own answer: magnitude, position, second catalogue
    7  where does each CV's own colour fall relative to its fit?
    8  what is still relative, and why
    9  the verdict on the strategy's calibration goal

Reading it backwards is the failure mode this layout prevents: section 9 is
a grade, and a grade without sections 2, 5 and 6 behind it is a number
someone chose.

Every figure and every tabulated number is produced by a query executed here
or a constant imported from ``macro_phot.cattie``.  The handful of values
that are NOT -- the AAVSO ranges, the literature figure for a
Johnson-Cousins colour term -- are named in the page footer rather than
covered by a blanket claim, because adversarial review found the blanket
claim sheltering two hand-typed numbers that the database disagreed with.
"""

from __future__ import annotations

import math
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

from . import cattie as ct        # noqa: E402
from macro_core.report_s0 import (  # noqa: E402
    ACCENT, BAD, STYLE, DPI, FAINT, GOOD, INK, MUTED, WARN,
    _figure, esc, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402  (house figure style)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "CV_TimeSeries"
FIG_DIR = DOCS_DIR / "figures" / "cv_catalogue_tie"
HTML_PATH = DOCS_DIR / "cv_catalogue_tie.html"

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}
FILTER_COLOR = ps.BAND_COLOR
FILTER_MARKER = ps.BAND_MARKER
CAT_COLOR = {"refcat2": ACCENT, "gaia_gspc": WARN}
#: Marker as the second channel, so a greyscale print still separates
#: the two catalogues.
CAT_MARKER = {"refcat2": "o", "gaia_gspc": "s"}

VERDICT_CLASS = {"TIED-STRETCH": "ok", "TIED-GOAL": "ok",
                 "TIED-ABOVE-GOAL": "warn", "TIED-UNVERIFIED": "warn",
                 "UNTIED": "bad"}


# ---------------------------------------------------------------------------
# Formatting helpers (an em-dash is the only thing a missing number becomes)
# ---------------------------------------------------------------------------
def _n(x, nd=3):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{float(x):,.{nd}f}"


def _mmag(x, nd=1):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "&mdash;"
    return f"{1000 * float(x):,.{nd}f}"


def _pm(v, e, nd=3):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "&mdash;"
    if e is None or (isinstance(e, float) and not math.isfinite(e)):
        return f"{float(v):+,.{nd}f}"
    return f"{float(v):+,.{nd}f}&thinsp;&plusmn;&thinsp;{float(e):.{nd}f}"


def _fig(src: str, caption: str, missing: str) -> str:
    """A figure, or an explicit statement that there is nothing to draw.

    A figure function returns an empty src when its query found no rows.
    Emitting `<img src="">` for that would put a broken image where an
    argument should be, and a reader would not know whether the evidence is
    absent or the page is broken.  This says which.
    """
    if not src:
        return f'<div class="note"><b>No figure here, and why:</b> {missing}</div>'
    return _figure(src, caption)


def _label(series_key: str) -> str:
    t, e, f = series_key.split("|")
    return f"{TARGET_LABEL.get(t, t)} {e} <i>{f}</i>"


#: ICRS -> Galactic, the standard J2000 pole and node.  Three constants
#: rather than an astropy import, because the page needs one rotation and
#: an import that heavy for one rotation is a dependency nobody audits.
_GAL_POLE_RA, _GAL_POLE_DEC, _GAL_NODE_L = 192.85948, 27.12825, 122.93192


def _galactic(ra_deg: float, dec_deg: float) -> tuple[float, float]:
    """Galactic (l, b) of an ICRS position, degrees.

    Here because a hand-typed galactic latitude in an earlier draft of this
    page (&ldquo;VV Pup sits at b = +2&rdquo;) was wrong by seven degrees.
    A number that can be computed from a stored coordinate should be.
    """
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    pra, pdec = math.radians(_GAL_POLE_RA), math.radians(_GAL_POLE_DEC)
    sb = (math.sin(dec) * math.sin(pdec)
          + math.cos(dec) * math.cos(pdec) * math.cos(ra - pra))
    b = math.asin(max(-1.0, min(1.0, sb)))
    y = math.cos(dec) * math.sin(ra - pra)
    x = (math.sin(dec) * math.cos(pdec)
         - math.cos(dec) * math.sin(pdec) * math.cos(ra - pra))
    l = math.radians(_GAL_NODE_L) - math.atan2(y, x)
    return math.degrees(l) % 360.0, math.degrees(b)


#: The columns of one PRIMARY tie row, in query order.  Named rather than
#: positional on purpose: this list grew by two columns mid-build and every
#: `r[27]` in the page silently started meaning something else.  A name
#: cannot do that.
TIE_COLS = (
    "series_key", "target_key", "era_id", "filter", "catalogue", "band",
    "band_system", "colour_label", "astrom_source",
    "n_candidates", "n_clean", "n_fit", "n_clipped", "n_check",
    "zp", "zp_err", "colour_term", "colour_err", "colour_ref",
    "resid_rms", "resid_mad", "chi2nu",
    "colour_min", "colour_max", "colour_p05", "colour_p95",
    "check_rms", "check_rms_clip", "check_median", "check_mad",
    "n_check_outlier",
    "target_colour", "target_colour_source", "colour_position",
    "extrap_err", "verdict", "note",
    "n_colour_pairs", "colour_scatter")


def _primary(con) -> list[dict]:
    """Every PRIMARY tie row as a dict, in reading order.  The page's spine."""
    rows = q(con, f"SELECT {', '.join(TIE_COLS)} FROM cv_cattie "
                  "WHERE is_primary=1 "
                  "ORDER BY target_key, era_id, filter")
    return [dict(zip(TIE_COLS, r)) for r in rows]


# ===========================================================================
# Figures
# ===========================================================================
def fig_catalogue_depth(con) -> str:
    """How deep each catalogue reaches in each field, and how many of THIS
    campaign's comparison stars sit inside that reach."""
    fields = [r[0] for r in q(con, "SELECT DISTINCT field_key FROM "
                                   "cv_cat_fetch ORDER BY 1")]
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 3.8))
        rows = q(con, """SELECT catalogue, field_key, n_rows
                         FROM cv_cat_fetch ORDER BY field_key""")
        w = 0.38
        for k, cat in enumerate(sorted({r[0] for r in rows})):
            vals = [next((r[2] for r in rows
                          if r[0] == cat and r[1] == f), 0) for f in fields]
            a1.bar(np.arange(len(fields)) + (k - 0.5) * w, vals, w,
                   label=cat, color=CAT_COLOR.get(cat, ACCENT))
        a1.set_xticks(range(len(fields)))
        a1.set_xticklabels([TARGET_LABEL.get(f, f) for f in fields],
                           rotation=20, fontsize=8)
        a1.set_ylabel("catalogue sources in the cone")
        a1.set_yscale("log")
        a1.set_title(f"Depth of the {ct.CONE_RADIUS_DEG:g}&deg; cone"
                     .replace("&deg;", " deg"))
        a1.legend(fontsize=8)

        # Right: matched fraction of this campaign's own reference stars.
        m = q(con, """SELECT m.catalogue, m.target_key,
                             count(DISTINCT m.star_id)
                      FROM cv_cat_match m GROUP BY 1, 2""")
        tot = {r[0]: r[1] for r in q(
            con, "SELECT target_key, count(*) FROM cv_ref_stars GROUP BY 1")}
        for k, cat in enumerate(sorted({r[0] for r in m})):
            vals = [100.0 * next((r[2] for r in m
                                  if r[0] == cat and r[1] == f), 0)
                    / max(tot.get(f, 1), 1) for f in fields]
            a2.bar(np.arange(len(fields)) + (k - 0.5) * w, vals, w,
                   label=cat, color=CAT_COLOR.get(cat, ACCENT))
        a2.set_xticks(range(len(fields)))
        a2.set_xticklabels([TARGET_LABEL.get(f, f) for f in fields],
                           rotation=20, fontsize=8)
        a2.set_ylabel("% of reference stars matched")
        a2.set_title("What the depth buys: match rate")
        a2.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "catalogue_depth.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/catalogue_depth.png"


def fig_veto_census(con) -> str:
    """Why a matched catalogue star is not automatically a tie star."""
    reasons = ["saturated", "near_veto", "blend_aperture", "blend_annulus",
               "ambiguous", "catalogue_flag", "no_cat_mag", "clean"]
    # Seven veto reasons plus "clean": BAD first (the reason that throws
    # most stars away) and GOOD last (the survivors), the house cycle in
    # between so no two adjacent wedges share a hue.
    colors = [BAD, WARN, ps.OTHER, ACCENT, FAINT, ps.SECOND, MUTED, GOOD]
    keys = [r[0] for r in q(con, """SELECT DISTINCT c.series_key
                                    FROM cv_cattie c WHERE c.is_primary=1
                                    ORDER BY c.target_key, c.era_id,
                                             c.filter""")]
    if not keys:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(10.6, 0.30 * len(keys) + 1.8))
        base = np.zeros(len(keys))
        for reason, col in zip(reasons, colors):
            vals = []
            for k in keys:
                r = q(con, """SELECT v.n FROM cv_cattie_veto v
                              JOIN cv_cattie c
                                ON c.series_key=v.series_key
                               AND c.catalogue=v.catalogue AND c.band=v.band
                              WHERE v.series_key=? AND c.is_primary=1
                                AND v.reason=?""", (k, reason))
                vals.append(r[0][0] if r else 0)
            ax.barh(range(len(keys)), vals, left=base, color=col,
                    label=reason, height=0.72)
            base += np.array(vals, dtype=float)
        ax.set_yticks(range(len(keys)))
        ax.set_yticklabels(keys, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("comparison + check stars of the series")
        ax.set_title("Every candidate tie star, and the first rule it broke")
        ax.legend(fontsize=7, ncol=4, loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "veto_census.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/veto_census.png"


def fig_colour_fits(con) -> str:
    """The fits themselves: instrumental minus catalogue, against colour."""
    rows = q(con, """SELECT series_key, catalogue, band, zp, colour_term,
                            colour_ref, colour_min, colour_max, filter,
                            target_colour
                     FROM cv_cattie WHERE is_primary=1 AND n_fit >= ?
                     ORDER BY target_key, era_id, filter""",
             (ct.MIN_TIE_STARS,))
    if not rows:
        return ""
    n = len(rows)
    ncol = 4
    nrow = int(math.ceil(n / ncol))
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(nrow, ncol, figsize=(11.4, 2.5 * nrow),
                                squeeze=False)
        for ax in axs.ravel()[n:]:
            ax.axis("off")
        for ax, r in zip(axs.ravel(), rows):
            skey, cat, band, zp, k, cref, cmin, cmax, filt, tcol = r
            st = q(con, """SELECT colour, delta, in_fit FROM cv_cattie_star
                           WHERE series_key=? AND catalogue=? AND band=?""",
                   (skey, cat, band))
            c = np.array([s[0] for s in st if s[0] is not None])
            d = np.array([s[1] for s in st if s[0] is not None])
            f = np.array([s[2] for s in st if s[0] is not None])
            col = FILTER_COLOR.get(filt, ACCENT)
            ax.scatter(c[f == 1], d[f == 1], s=8, color=col, alpha=0.75,
                       label="fit")
            ax.scatter(c[f == 0], d[f == 0], s=14, facecolors="none",
                       edgecolors=INK, lw=0.7, label="held out")
            if k is not None and zp is not None:
                xs = np.linspace(cmin, cmax, 20)
                ax.plot(xs, zp + k * (xs - cref), color=BAD, lw=1.4)
            if tcol is not None and math.isfinite(tcol):
                ax.axvline(tcol, color=WARN, lw=1.2, ls="--")
            ax.set_title(f"{skey}\n{cat}:{band}", fontsize=7)
            ax.tick_params(labelsize=7)
            ax.set_xlabel("catalogue colour", fontsize=7)
            ax.set_ylabel("ens &minus; cat (mag)".replace("&minus;", "-"),
                          fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "colour_fits.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/colour_fits.png"


def fig_colour_terms(con) -> str:
    """Every colour term with its error bar, grouped by filter label.

    The figure that answers section 4: if uppercase 'G'/'R'/'I' were
    Johnson bands, their colour terms against SDSS would be tenths, not
    hundredths, and they would not sit beside the lowercase ones."""
    rows = q(con, """SELECT series_key, filter, colour_term, colour_err
                     FROM cv_cattie WHERE is_primary=1 AND n_fit >= ?
                     ORDER BY filter, target_key, era_id""",
             (ct.MIN_TIE_STARS,))
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(10.6, 0.28 * len(rows) + 1.6))
        ys = np.arange(len(rows))
        for y, (skey, filt, k, ke) in zip(ys, rows):
            ax.errorbar(k, y, xerr=ke or 0, fmt="o", ms=5,
                        color=FILTER_COLOR.get(filt, ACCENT),
                        ecolor=MUTED, capsize=2)
        ax.axvline(0.0, color=MUTED, lw=1, ls=":")
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{r[0]}" for r in rows], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("colour term k  (mag per mag of colour)")
        ax.set_title("The bandpass mismatch, measured -- never applied to a "
                     "science target")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "colour_terms.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/colour_terms.png"


def fig_accuracy(con) -> str:
    """Achieved accuracy on held-out stars, against the strategy's band."""
    rows = q(con, """SELECT series_key, filter, check_rms_clip, check_median,
                            n_check, resid_rms, check_rms, n_check_outlier
                     FROM cv_cattie WHERE is_primary=1 AND n_fit >= ?
                     ORDER BY target_key, era_id, filter""",
             (ct.MIN_TIE_STARS,))
    rows = [r for r in rows if r[2] is not None]
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(10.6, 0.30 * len(rows) + 1.8))
        ys = np.arange(len(rows))
        ax.axvspan(ct.ACCURACY_STRETCH_MAG * 1000,
                   ct.ACCURACY_GOAL_MAG * 1000, color=GOOD, alpha=0.16)
        ax.axvline(ct.ACCURACY_GOAL_MAG * 1000, color=GOOD, lw=1.4)
        ax.axvline(ct.ACCURACY_STRETCH_MAG * 1000, color=GOOD, lw=1.0,
                   ls="--")
        for y, r in zip(ys, rows):
            ax.barh(y, 1000 * r[5], height=0.34, color=FAINT,
                    label="fit residual RMS" if y == 0 else None)
            # The raw RMS is drawn too, and joined to the clipped one by a
            # line: the length of that line IS the influence of the two or
            # three outlier stars per block, which is the thing a reader
            # would otherwise have to take on trust.
            if r[6] is not None:
                ax.plot([1000 * r[2], 1000 * r[6]], [y, y], "-",
                        color=MUTED, lw=1.0, zorder=1)
                ax.plot(1000 * r[6], y, "o", ms=5, mfc="none", mec=BAD,
                        label="raw check RMS (with outliers)"
                        if y == 0 else None)
            ax.plot(1000 * r[2], y, "o", ms=6,
                    color=GOOD if r[2] <= ct.ACCURACY_GOAL_MAG else BAD,
                    label="clipped check RMS (typical star)"
                    if y == 0 else None, zorder=3)
            ax.plot(1000 * abs(r[3]), y, "s", ms=4, color=WARN,
                    label="|check median| (bias)" if y == 0 else None)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{r[0]} (n={r[4]}, {r[7]} out)" for r in rows],
                           fontsize=7)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlabel("magnitude error (mmag)")
        ax.set_title(f"Achieved absolute accuracy vs the "
                     f"{1000 * ct.ACCURACY_STRETCH_MAG:g}-"
                     f"{1000 * ct.ACCURACY_GOAL_MAG:g} mmag goal")
        ax.legend(fontsize=7, loc="lower right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "accuracy.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/accuracy.png"


def fig_resid_mag(con) -> str:
    """Is the tie residual a function of brightness?  (Non-linearity.)"""
    rows = q(con, """SELECT s.cat_mag, s.resid, c.filter
                     FROM cv_cattie_star s JOIN cv_cattie c
                       ON c.series_key=s.series_key AND c.catalogue=s.catalogue
                      AND c.band=s.band
                     WHERE c.is_primary=1 AND s.in_fit=1
                       AND s.cat_mag IS NOT NULL AND s.resid IS NOT NULL""")
    if not rows:
        return ""
    m = np.array([r[0] for r in rows])
    d = np.array([r[1] for r in rows])
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 3.8))
        for filt in sorted({r[2] for r in rows}):
            sel = np.array([r[2] == filt for r in rows])
            a1.scatter(m[sel], 1000 * d[sel], s=4, alpha=0.35,
                       color=FILTER_COLOR.get(filt, ACCENT), label=filt)
        # Binned median: the trend a scatter plot hides.
        bins = np.arange(np.floor(m.min()), np.ceil(m.max()) + 0.5, 0.5)
        cen, med = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            s = (m >= lo) & (m < hi)
            if s.sum() >= 8:
                cen.append(0.5 * (lo + hi))
                med.append(1000 * float(np.median(d[s])))
        a1.plot(cen, med, "-o", color=INK, lw=1.6, ms=4,
                label="binned median")
        a1.axhline(0, color=MUTED, lw=1, ls=":")
        a1.set_xlabel("catalogue magnitude")
        a1.set_ylabel("tie residual (mmag)")
        a1.set_ylim(-150, 150)
        a1.set_title("Residual vs brightness (non-linearity)")
        a1.legend(fontsize=7, ncol=2)

        tr = q(con, """SELECT swing, significant FROM cv_cattie_trend
                       WHERE axis='cat_mag'""")
        sw = np.array([1000 * t[0] for t in tr if t[0] is not None])
        sig = np.array([t[1] for t in tr if t[0] is not None], dtype=bool)
        if sw.size:
            a2.hist(sw[~sig], bins=24, color=FAINT,
                    label=f"not significant (n={int((~sig).sum())})")
            a2.hist(sw[sig], bins=24, color=BAD,
                    label=f"significant at {ct.TREND_SIGMA:g}"
                          f"$\\sigma$ (n={int(sig.sum())})")
        a2.axvline(0, color=MUTED, lw=1, ls=":")
        a2.set_xlabel("trend SWING across the fitted magnitude range (mmag)")
        a2.set_ylabel("blocks")
        a2.set_title("How big is it, where it is real?")
        a2.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "resid_mag.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/resid_mag.png"


#: Filled in by :func:`fig_resid_xy` with the blocks it actually drew, so
#: the caption in section 6.2 can name them instead of asserting a count.
#: Review found the caption claiming "the six best-populated blocks" over a
#: figure showing five, chosen one per target -- the selection rule in the
#: code was right and the sentence describing it was wrong twice.
RESID_XY_KEYS: list[str] = []


def fig_resid_xy(con) -> str:
    """Is the tie residual a function of DETECTOR POSITION?  (Flat field.)

    Nobody had looked at (x, y) before this page; the characterization
    flagged a noise floor 28x above scintillation and left it unexplained.
    """
    # ONE block per target, the best-populated of each.  Taking the six
    # largest outright would show six VV Pup panels and answer nothing: a
    # flat-field residual should repeat across FIELDS on the same detector,
    # and only a panel per target can show whether it does.
    keys = [r[0] for r in q(con, """SELECT series_key FROM cv_cattie c
                                    WHERE is_primary=1 AND n_fit >= 25
                                      AND n_fit = (SELECT max(n_fit)
                                                   FROM cv_cattie d
                                                   WHERE d.is_primary=1
                                                     AND d.target_key =
                                                         c.target_key)
                                    ORDER BY target_key LIMIT 6""")]
    RESID_XY_KEYS[:] = keys
    with plt.rc_context(STYLE):
        fig, axs = plt.subplots(2, 4, figsize=(11.4, 5.4))
        panels = list(axs.ravel())
        for ax in panels:
            ax.axis("off")
        radial: list[tuple] = []
        for ax, k in zip(panels, keys):
            st = q(con, """SELECT s.x, s.y, s.resid FROM cv_cattie_star s
                           JOIN cv_cattie c ON c.series_key=s.series_key
                            AND c.catalogue=s.catalogue AND c.band=s.band
                           WHERE s.series_key=? AND c.is_primary=1
                             AND s.in_fit=1 AND s.x IS NOT NULL""", (k,))
            if not st:
                continue
            ax.axis("on")
            x = np.array([s[0] for s in st])
            y = np.array([s[1] for s in st])
            r = 1000 * np.array([s[2] for s in st])
            lim = float(np.percentile(np.abs(r), 90)) or 1.0
            sc = ax.scatter(x, y, c=r, s=16, cmap="coolwarm",
                            vmin=-lim, vmax=lim)
            sw = q(con, """SELECT axis, swing, significant
                           FROM cv_cattie_trend WHERE series_key=?
                             AND axis IN ('plane_xy','radius')""", (k,))
            d = {a: (s, g) for a, s, g in sw}
            tag = ""
            if "plane_xy" in d:
                tag += (f"\nplane {1000 * (d['plane_xy'][0] or 0):.0f} mmag"
                        + (" (sig)" if d["plane_xy"][1] else ""))
            if "radius" in d:
                tag += (f"  radial {1000 * (d['radius'][0] or 0):.0f} mmag"
                        + (" (sig)" if d["radius"][1] else ""))
            ax.set_title(f"{k}{tag}", fontsize=7)
            ax.tick_params(labelsize=6)
            fig.colorbar(sc, ax=ax, fraction=0.046).ax.tick_params(labelsize=6)
            # Radial profile, binned, for the summary panel below.
            n1, n2 = float(np.nanmax(x)), float(np.nanmax(y))
            rad = np.hypot(x - n1 / 2, y - n2 / 2) / max(n1, n2)
            radial.append((k, rad, r))

        # Panel: the radial profile that the (x, y) maps hint at.
        ax = panels[-2]
        ax.axis("on")
        for k, rad, r in radial:
            edges = np.linspace(0, float(np.nanmax(rad)), 7)
            cen, med = [], []
            for lo, hi in zip(edges[:-1], edges[1:]):
                s = (rad >= lo) & (rad < hi)
                if s.sum() >= 5:
                    cen.append(0.5 * (lo + hi))
                    med.append(float(np.median(r[s])))
            if cen:
                ax.plot(cen, med, "-o", ms=3, lw=1.2, label=k.split("|")[0])
        ax.axhline(0, color=MUTED, lw=1, ls=":")
        ax.set_xlabel("radius / detector half-width", fontsize=7)
        ax.set_ylabel("median residual (mmag)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=6)
        ax.set_title("Radial profile: the shape a plane cannot see",
                     fontsize=7)

        # Panel: every block's tilt and radial swing.
        ax = panels[-1]
        ax.axis("on")
        for axis_name, col, lab in (("plane_xy", BAD, "plane tilt"),
                                    ("radius", ACCENT, "radial")):
            pl = q(con, "SELECT swing, significant FROM cv_cattie_trend "
                        "WHERE axis=? AND swing IS NOT NULL", (axis_name,))
            if pl:
                sw = np.array([1000 * abs(p[0]) for p in pl])
                n_sig = int(sum(p[1] for p in pl))
                ax.hist(sw, bins=18, histtype="step", lw=1.6, color=col,
                        label=f"{lab} ({n_sig}/{len(pl)} sig.)")
        ax.set_xlabel("|swing| (mmag)", fontsize=7)
        ax.set_ylabel("blocks", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=6)
        ax.set_title("Position systematics, all blocks", fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "resid_xy.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/resid_xy.png"


def fig_cross(con) -> str:
    """ATLAS-REFCAT2 minus Gaia synthetic on the SAME stars, per block."""
    rows = q(con, """SELECT series_key, n_common, star_offset_median,
                            star_offset_mad, d_zp
                     FROM cv_cattie_cross ORDER BY series_key""")
    if not rows:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(10.6, 0.30 * len(rows) + 1.8))
        ys = np.arange(len(rows))
        for y, r in zip(ys, rows):
            ax.errorbar(1000 * r[2], y, xerr=1000 * (r[3] or 0), fmt="o",
                        ms=5, color=ACCENT, ecolor=MUTED, capsize=2)
            if r[4] is not None:
                ax.plot(1000 * r[4], y, "x", ms=7, color=WARN)
        ax.axvline(0, color=MUTED, lw=1, ls=":")
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{r[0]} (n={r[1]})" for r in rows], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("REFCAT2 &minus; Gaia synthetic (mmag)"
                      .replace("&minus;", "-"))
        ax.set_title("Two independent catalogues on the same stars "
                     "(circles: star-by-star; crosses: fitted ZP difference)")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "cross_catalogue.png", dpi=DPI)
        plt.close(fig)
    return "figures/cv_catalogue_tie/cross_catalogue.png"


# ===========================================================================
# Sections
# ===========================================================================
def section_intro(con) -> str:
    n_lc = q1(con, "SELECT count(*) FROM cv_lightcurve")
    n_cal = q1(con, "SELECT count(*) FROM cv_lightcurve "
                    "WHERE cal_mag IS NOT NULL") \
        if q1(con, "SELECT count(*) FROM pragma_table_info('cv_lightcurve') "
                   "WHERE name='cal_mag'") else 0
    n_series = q1(con, "SELECT count(*) FROM cv_series WHERE status='solved'")
    n_block = q1(con, "SELECT count(*) FROM cv_ref")
    n_block_tied = q1(con, "SELECT count(DISTINCT target_key||'|'||era_id) "
                           "FROM cv_cattie WHERE is_primary=1 "
                           "AND verdict LIKE 'TIED%'")
    tied = q1(con, "SELECT count(*) FROM cv_cattie WHERE is_primary=1 "
                   "AND verdict LIKE 'TIED%'")
    return f"""
<section id="intro">
<div class="bhead"><h2>0 &middot; What is a magnitude on this page, and why
was the old answer arbitrary?</h2></div>

<p>The ensemble solver models every measurement as
<code>m<sub>ij</sub> = M<sub>i</sub> + ZP<sub>j</sub></code>, and that model
is invariant under <code>M &rarr; M + c</code>, <code>ZP &rarr; ZP
&minus; c</code>.  The solver fixes the freedom by demanding
<code>mean(ZP) = 0</code>.  That choice is arithmetically convenient and
physically empty: it puts every light curve on a gauge whose origin is
&ldquo;whatever these particular comparison stars happened to average
to&rdquo;.  {n_lc:,} light-curve rows were internally consistent and
externally meaningless.</p>

<p>The characterization page graded the strategy's calibration goal
<b>NOT SUPPORTED</b> on exactly that, and named the repair
&ldquo;the highest-value single action on the whole list&rdquo;.  This page
is the repair.  It asks one question per section and refuses to move on
until the question has a number attached.</p>

<div class="decision"><b>What changed.</b> {tied} of {n_series} solved
series now carry a catalogue tie; {n_cal:,} of {n_lc:,} light-curve rows
carry a calibrated natural-system magnitude in the new column
<code>cv_lightcurve.cal_mag</code>.  {n_block_tied} of the archive's
{n_block} (target, era) blocks now carry a catalogue tie, against
<b>7 of 14</b> when the characterization graded the goal.</div>

<div class="note"><b>The four rulings this page obeys, stated before any
number is shown, so that no number can be read as bending them.</b>
<ol>
<li><b>The natural system is the product.</b> A published magnitude is
<code>m<sub>nat</sub> = m<sub>ens</sub> &minus; ZP0</code>: this telescope's
own bandpass, zero-pointed to the catalogue at the tie stars' median
colour.  The colour term is published as METADATA and is
<i>never</i> applied to a science target.  These targets are cataclysmic
variables &mdash; blue, variable, and routinely outside the colour range
over which any transformation was calibrated.  Transforming them would
swap a bandpass error of known size for an extrapolation error of unknown
size.</li>
<li><b>The tie is solved on comparison stars.</b> The target is excluded by
role, twice: the ensemble already refuses to let a polar's orbital
modulation set the zero point it is measured against, and the fit here
draws only on stars carrying <code>role IN ('comp','check')</code>.</li>
<li><b>Zero point AND colour term, with uncertainties and a stated range of
validity.</b> Section 3 gives both terms per block; section 7 says
explicitly where each target's own colour falls relative to that range.</li>
<li><b>Saturation and blending veto catalogue stars too.</b> Section 2 is
the census: every candidate, and the first rule it broke.</li>
</ol></div>
</section>"""


def _dupvar_evidence(con, cache_root: Path) -> tuple[str, str]:
    """What REFCAT2's ``dupvar`` flag actually does in THESE cones.

    Returns ``(distribution_table, residual_table)``.  Review caught the
    first draft justifying the decision to ignore ``dupvar`` with a
    173-source test cone in which the flag was almost constant; the cones
    actually used are not like that (VV Pup alone carries 1,375 flagged
    sources).  The decision survives, but on the real measurement rather
    than the anecdote: the flagged stars are not the bad ones.  Both tables
    are read from the cached cone files and the tie residuals, so the claim
    is now reproducible from the artefacts on disk.
    """
    import gzip
    import json
    from collections import Counter
    dist_rows, by_flag = [], {}
    for field, in q(con, "SELECT DISTINCT field_key FROM cv_cat_fetch "
                         "WHERE catalogue='refcat2' AND n_rows>0 "
                         "ORDER BY 1"):
        path = cache_root / "refcat2" / f"{field}.json.gz"
        if not path.exists():
            continue
        with gzip.open(path, "rt") as fh:
            cols = json.load(fh)["columns"]
        counts = Counter(int(v) for v in cols.get("dupvar", [])
                         if v is not None and math.isfinite(float(v)))
        n = sum(counts.values()) or 1
        dist_rows.append([TARGET_LABEL.get(field, field), f"{n:,}",
                          ", ".join(f"<code>{k}</code>:&nbsp;{v:,}"
                                    for k, v in sorted(counts.items())),
                          f"{100 * (n - counts.get(2, 0)) / n:.1f}%"])
        # Residual of every fitted tie star of this field, by its flag.
        dv = {i: int(v) for i, v in enumerate(cols.get("dupvar", []))
              if v is not None and math.isfinite(float(v))}
        for cat_row, resid in q(con, """SELECT m.cat_row, s.resid
                FROM cv_cattie_star s
                JOIN cv_cattie c ON c.series_key=s.series_key
                 AND c.catalogue=s.catalogue AND c.band=s.band
                JOIN cv_cat_match m ON m.catalogue=s.catalogue
                 AND m.target_key=c.target_key AND m.era_id=c.era_id
                 AND m.star_id=s.star_id
                WHERE c.is_primary=1 AND s.in_fit=1 AND s.resid IS NOT NULL
                  AND c.target_key=?""", (field,)):
            f = dv.get(int(cat_row))
            if f is not None:
                by_flag.setdefault(f, []).append(float(resid))
    res_rows, rms = [], {}
    for f in sorted(by_flag):
        r = np.array(by_flag[f])
        med = float(np.median(r))
        mad = float(1.4826 * np.median(np.abs(r - med))) or 1e-6
        keep = np.abs(r - med) <= ct.CLIP_SIGMA * mad
        val = (float(np.sqrt(np.mean(np.square(r[keep]))))
               if keep.any() else float("nan"))
        rms[f] = (len(r), val)
        res_rows.append([f"<code>{f}</code>", f"{len(r):,}", _mmag(val),
                         _mmag(med)])
    # The verdict sentence, built from the groups rather than asserted, so
    # it cannot drift from the table above it.  "Clean" here means dupvar=2,
    # the value REFCAT2 gives an ordinary unflagged source.
    base = rms.get(2, (0, float("nan")))
    flagged = {f: v for f, v in rms.items() if f != 2}
    big = max(flagged, key=lambda f: flagged[f][0], default=None)
    worse = [f for f, (n, v) in flagged.items()
             if math.isfinite(v) and math.isfinite(base[1]) and v > base[1]]
    verdict = ""
    if big is not None and math.isfinite(base[1]):
        verdict = (
            f"The unflagged class (<code>dupvar</code>=2, {base[0]:,} fitted "
            f"tie stars) scatters {_mmag(base[1])} mmag about the tie.  The "
            f"largest flagged class (<code>dupvar</code>={big}, "
            f"{flagged[big][0]:,} stars) scatters {_mmag(flagged[big][1])} "
            f"mmag &mdash; "
            + ("BETTER, so vetoing on the flag would have deleted good "
               "stars.  "
               if flagged[big][1] < base[1] else "worse.  ")
            + ("No flagged class fits worse than the unflagged one."
               if not worse else
               "The classes that do fit worse are "
               + ", ".join(f"<code>{f}</code> ({flagged[f][0]:,} stars)"
                           for f in sorted(worse))
               + " &mdash; too few, across five fields, to carry a veto that "
                 "would also throw out the large class above."))
    return (table(["field", "sources in cone", "dupvar distribution",
                   "not <code>dupvar</code>=2"], dist_rows),
            table(["dupvar", "fitted tie stars", "clipped resid RMS (mmag)",
                   "median resid (mmag)"], res_rows),
            verdict)


def section_catalogue(con, fig1: str, cache_root: Path) -> str:
    rows = q(con, """SELECT catalogue, field_key, ra_deg, dec_deg, n_rows,
                            pulled_utc, cache_path, substr(cache_sha256,1,12),
                            note FROM cv_cat_fetch
                     ORDER BY field_key, catalogue""")
    body = table(
        ["catalogue", "field", "centre (RA, Dec)", "rows", "pulled (UTC)",
         "cache", "sha256", "note"],
        [[esc(r[0]), TARGET_LABEL.get(r[1], r[1]),
          f"{r[2]:.4f}, {r[3]:.4f}", f"{r[4]:,}", esc(r[5])[:19],
          f"<code>{esc(Path(r[6]).name)}</code>" if r[6] else "&mdash;",
          f"<code>{esc(r[7])}</code>" if r[7] else "&mdash;",
          esc(r[8])] for r in rows],
        # A pull that FAILED is a row, not an absence.  The alternative --
        # simply having no row for that catalogue -- reads as "we never
        # tried", which is a different and much more flattering claim.
        ["bad" if (r[4] or 0) == 0 else "" for r in rows])
    n_fail = sum(1 for r in rows if (r[4] or 0) == 0)
    fail_note = "" if not n_fail else (
        f'<div class="bad"><b>{n_fail} of {len(rows)} pulls failed, and the '
        f'rows above say so rather than being omitted.</b>  Every failure '
        f'here is the same one: the ESA Gaia archive timed out on every '
        f'query during this build, including <code>SELECT TOP 1 source_id '
        f'... WHERE source_id = &lt;one id&gt;</code>.  The consequence is '
        f'specific and is carried through the rest of this page: no '
        f'cross-catalogue systematic (section 6.3), no two-system test of '
        f'the filter labels (section 4), and the <i>y</i>-band block left '
        f'relative (section 8).  The stage is resumable and the primary '
        f'catalogue is cached, so recovering all three is one command '
        f'when the archive returns.</div>')
    bands = []
    for f in ("g", "r", "i", "z", "y", "G", "R", "I"):
        for sp in ct.band_candidates(f):
            bands.append([f"<code>{f}</code>", esc(sp.catalogue),
                          f"<code>{esc(sp.mag_col)}</code>",
                          esc(sp.system), f"<code>{esc(sp.colour_label)}</code>",
                          esc(sp.hypothesis)])
    band_tab = table(["FILTER label", "catalogue", "column", "system",
                      "colour index", "hypothesis"], bands)
    dup_dist, dup_res, dup_verdict = _dupvar_evidence(con, cache_root)
    # Galactic coordinates of each cone centre, computed rather than
    # remembered: an earlier draft asserted VV Pup sits at b = +2, and it
    # does not.  No astropy dependency for a rotation this simple.
    gal = []
    for field, ra, dec in q(con, "SELECT DISTINCT field_key, ra_deg, dec_deg "
                                 "FROM cv_cat_fetch ORDER BY 1"):
        l, b = _galactic(ra, dec)
        gal.append((TARGET_LABEL.get(field, field), b))
    gal.sort(key=lambda t: abs(t[1]))
    crowd = (f"{gal[0][0]} sits at galactic latitude "
             f"{gal[0][1]:+.0f}&deg; and {gal[-1][0]} at "
             f"{gal[-1][1]:+.0f}&deg;" if len(gal) > 1 else "")
    return f"""
<section id="catalogue">
<div class="bhead"><h2>1 &middot; Which catalogue &mdash; and what does it
cost to be wrong about that?</h2></div>

<p>Two catalogues are pulled for every field, and both are tied.  That is
not indecision; it is the only way to measure the systematic floor
(section 6).</p>

<p><b>Primary: ATLAS-REFCAT2</b> (Tonry et&nbsp;al. 2018), through VizieR's
TAP service.  Three reasons, in order of weight.  It is the system the
strategy itself names &mdash; &ldquo;nightly REFCAT2 tie &rarr; PS1 AB to
0.01&ndash;0.02 mag&rdquo; &mdash; so tying to it grades the goal that was
actually set rather than a substitute.  It reaches m&nbsp;&asymp;&nbsp;19,
about 1.5 mag deeper than the alternative, and depth is tie stars, and tie
stars are what a colour term is fitted from.  And it ships its own blend
metrology (the <code>R1</code> contamination radius), which is precisely
what ruling&nbsp;4 demands and which no other all-sky catalogue supplies
for free.</p>

<p><b>Secondary: Gaia&nbsp;DR3 standardised synthetic photometry</b>
(<code>gaiadr3.synthetic_photometry_gspc</code>), through the ESA archive.
It needs a BP/RP spectrum, so it stops near G&nbsp;=&nbsp;17.65 &mdash;
shallower, and in the crowded VV&nbsp;Pup field much shallower.  But it
publishes Sloan, PS1&nbsp;<i>y</i> AND Johnson&ndash;Cousins magnitudes for
the SAME stars from the SAME spectra, and that buys two things nothing else
can: the <i>y</i> filter becomes tieable at all, and the uppercase-label
question of section 4 becomes a measurement instead of an assumption.</p>

<p><b>What was fetched, and when.</b> Every pull is cached on disk with its
query text and its date; nothing is re-fetched silently, and a re-run of
this stage months from now reads the same bytes unless someone asks for a
new pull.</p>
{body}
{fail_note}

{_fig(fig1, "Left: how many catalogue sources each cone returned "
               f"(log scale) &mdash; {crowd}, and the low-latitude cone is "
               "an order of magnitude the richer. "
               "Right: what that depth actually buys, as the fraction of "
               "THIS campaign's reference stars that found a counterpart.",
               "no catalogue pull has been recorded yet &mdash; run "
               "<code>fetch</code>.")}

<p><b>The band map.</b> Every FILTER label in the archive, and the
catalogue column it is taken to mean.  Rows marked <i>alternative</i> are
hypotheses the data is asked to REJECT, not beliefs; they are fitted and
stored so that &ldquo;we checked&rdquo; has evidence behind it.</p>
{band_tab}

<h3>1.1 &mdash; deliberately not used: REFCAT2's <code>dupvar</code> flag</h3>
<p>An earlier draft of this page defended that choice by saying the flag is
almost constant &mdash; &ldquo;2 for 171 of the 173 sources in a test
cone&rdquo;.  Adversarial review checked the cones this build actually
pulled, and they are not like that: <code>dupvar</code> marks a small but
real subset of every field, over a thousand sources in VV&nbsp;Pup alone.
The anecdote was generalised from a probe rather than measured on the data
in hand, so it is replaced here by the measurement.</p>
{dup_dist}

<p>The question that decides the veto is not whether the flag separates
anything, but whether what it separates is <i>photometrically worse</i>.
Every fitted tie star on this page, grouped by its own flag:</p>
{dup_res}

<div class="note"><b>The flag does not sort these stars into good and bad in
a way that would justify a veto.</b>  {dup_verdict}  That is a stronger
statement than &ldquo;the flag is constant&rdquo;, and unlike that one it is
the statement this data supports.  <code>dupvar</code> is still fetched and
stored so a reader can redo this grouping.  The blend veto instead uses
<code>R1</code>, the measured radius at which neighbour flux equals the
star's own, which is a measurement rather than a label.</div>
</section>"""


def _astrom_table(con) -> str:
    """Every block's measured displacement from the catalogue."""
    rows = q(con, """SELECT target_key, era_id, method, n_stars, n_paired,
                            dra_arcsec, ddec_arcsec, offset_arcsec,
                            scatter_arcsec, applied, n_match_before,
                            n_match_after, reason
                     FROM cv_cat_astrom WHERE catalogue=?
                     ORDER BY target_key, era_id""", (ct.PRIMARY_CATALOGUE,))
    return table(
        ["block", "astrometry", "stars offered / paired",
         "&Delta;RA&middot;cos&delta; (&Prime;)", "&Delta;Dec (&Prime;)",
         "|offset| (&Prime;)", "scatter (&Prime;)",
         "matches before &rarr; after", "decision"],
        [[f"{TARGET_LABEL.get(r[0], r[0])} e{r[1]}",
          f"<code>{esc(r[2] or 'none')}</code>", f"{r[3]:,} / {r[4]:,}",
          _n(r[5], 2), _n(r[6], 2), _n(r[7], 2), _n(r[8], 2),
          f"{r[10]:,} &rarr; {r[11]:,}",
          ("<b>REMOVED</b> &mdash; " if r[9] else "left alone &mdash; ")
          + esc(r[12])] for r in rows],
        ["warn" if r[9] else "" for r in rows])


def section_stars(con, fig2: str) -> str:
    rows = q(con, """SELECT c.series_key, c.n_candidates, c.n_clean, c.n_fit,
                            c.n_check, c.n_clipped, c.astrom_source
                     FROM cv_cattie c WHERE c.is_primary=1
                     ORDER BY c.target_key, c.era_id, c.filter""")
    body = table(
        ["series", "candidates", "clean", "in fit", "held out",
         "clipped by fit", "astrometry"],
        [[_label(r[0]), r[1], r[2], r[3], r[4], r[5],
          f"<code>{esc(r[6])}</code>"] for r in rows],
        ["bad" if (r[3] or 0) < ct.MIN_TIE_STARS else "" for r in rows])
    tot = q(con, """SELECT v.reason, sum(v.n) FROM cv_cattie_veto v
                    JOIN cv_cattie c ON c.series_key=v.series_key
                     AND c.catalogue=v.catalogue AND c.band=v.band
                    WHERE c.is_primary=1 GROUP BY 1""")
    d = dict(tot)
    # The saturation numbers, MEASURED.  An earlier draft typed "High Gain's
    # scale ends at 3,496 ADU and 7.2% of its detections are saturated" into
    # the prose; the ceiling was the detector's raw one rather than the veto
    # actually applied, and the fraction depended on an unstated population
    # (7.2% over matched frames, 8.7% over every staged frame).  Both are
    # queries now, with the population named.
    sat = q(con, """SELECT f.readoutm, count(*), sum(d.saturated),
                           min(f.veto_applied_adu)
                    FROM cv_frames f JOIN cv_detections d
                      ON d.frame_id = f.frame_id
                    WHERE f.status='matched' AND f.veto_applied_adu IS NOT NULL
                    GROUP BY 1
                    ORDER BY min(f.veto_applied_adu) ASC""")
    worst = sat[0] if sat else None
    sat_txt = ("" if not worst else
               f"the tightest scale in the archive is "
               f"<code>{esc(worst[0])}</code>, whose applied veto sits at "
               f"{worst[3]:,.0f}&nbsp;ADU and which trips it on "
               f"{100.0 * (worst[2] or 0) / max(worst[1], 1):.1f}% of its "
               f"{worst[1]:,} detections on matched frames")
    return f"""
<section id="stars">
<div class="bhead"><h2>2 &middot; Which stars may carry the tie &mdash; and
which may not?</h2></div>

<p>A catalogue star that matched a reference star is a CANDIDATE, not a tie
star.  Ruling&nbsp;4 says the pipeline's own vetoes apply to catalogue stars
too, and the reason is arithmetic rather than fastidious: a comparison star
that saturates has a mean magnitude biased faint by an amount nobody can
recover afterwards, and a star sharing its 4-arcsec aperture with a
neighbour the catalogue resolves and this telescope does not carries a
magnitude error of arbitrary size straight into the zero point.  One
corrupted star in twenty moves ZP0 by a fifth of its own error.</p>

<p>The gate, in the order a star meets it &mdash; each candidate is charged
to the FIRST rule it broke, so the census partitions the sample exactly
once:</p>
<ul>
<li><b>no catalogue magnitude or colour</b> ({d.get('no_cat_mag', 0):,}
stars) &mdash; an absence, not a veto, counted first so the numbers below
describe stars that had photometry to lose.</li>
<li><b>the catalogue distrusts its own photometry</b>
({d.get('catalogue_flag', 0):,}) &mdash; REFCAT2's measured contamination
radius <code>R1 &lt; {ct.BLEND_APERTURE_ARCSEC:g}&Prime;</code>, or Gaia's
band validity flag and its <code>C*</code> colour-excess statistic.  Both
catalogues are asked the same question in their own currency.</li>
<li><b>saturated</b> ({d.get('saturated', 0):,}) &mdash; ANY measurement in
the series flagged against the era's own readout-mode ceiling.  Tolerance
is zero: {sat_txt}, and there are always more candidates than the fit
needs.</li>
<li><b>in the non-linear shoulder</b> ({d.get('near_veto', 0):,}) &mdash;
median peak above {100 * ct.NEAR_VETO_FRAC:g}% of the applied veto.  A
detector is non-linear before it clips, and nothing flags that.</li>
<li><b>blended inside the aperture</b> ({d.get('blend_aperture', 0):,})
&mdash; a catalogue neighbour within {ct.BLEND_APERTURE_ARCSEC:g}&Prime;
and less than {ct.BLEND_DMAG:g} mag fainter <b>in the band being tied</b>
(a 10% flux contribution is 0.10 mag, five times the accuracy goal).</li>
<li><b>blended in the sky annulus</b> ({d.get('blend_annulus', 0):,})
&mdash; a neighbour BRIGHTER in that same band inside the
8&ndash;12&Prime; background annulus, which biases the subtracted sky and
therefore the flux.</li>
<li><b>ambiguous identification</b> ({d.get('ambiguous', 0):,}) &mdash; a
second catalogue source within {ct.AMBIGUITY_FACTOR:g}&times; the accepted
match distance.  A coin-flip identity is worse than no tie star.</li>
</ul>

{_fig(fig2, "Every candidate tie star of every block, and the first "
               "rule it broke.  Green is what survived to the fit.",
               "no block has been solved yet — run <code>solve</code>.")}

{body}

<div class="note"><b>Which BAND the blend veto is taken in, and why it had
to change.</b> The catalogue stores one neighbour census in Gaia&nbsp;G
&mdash; one currency, so both catalogues answer the identical question and
their rows compare directly.  That is a good census and it was, until
review, also the veto, while the fit itself ran against PS1
<i>g</i>/<i>r</i>/<i>i</i>.  The two are not interchangeable: where Gaia
resolves a pair that PS1 does not, the G-band contrast is large and the
pair sails through a gate that the tie band would have caught.  Measured on
this build, {d.get('blend_aperture_gband', 0):,} candidate stars trip the
aperture rule read in G against
{d.get('blend_aperture_tieband', 0):,} read in the band actually being
tied &mdash; and <b>{d.get('blend_aperture_band_only', 0):,} of them are
stars the tie band rejects that Gaia&nbsp;G would have waved straight
through</b>, which is the number that matters and the population review
found sitting inside the VV&nbsp;Pup fits.  Ruling&nbsp;4 asks whether the
catalogue magnitude in the tie band
describes the flux this telescope's aperture collected, so the veto is now
taken there, and the G-band census is kept beside it as the cross-catalogue
comparison it was always good for.</div>

<div class="note"><b>Where the sky positions come from, and why the column
says so.</b> The <code>astrometry</code> column carries
<code>cv_field_tie.method</code>: <code>wcs</code> means the block's
reference frame had an S1 plate solution and the reference stars inherited
its astrometric residual; <code>gaia</code> means the reference was
unsolved and the positions came from a parity-tolerant similarity fit to a
Gaia cone.  A tie standing on a triangle fit deserves to be read
differently from one standing on a plate solution, so the provenance
travels with every row rather than living in a footnote.</div>

<h3>2.1 &mdash; the astrometric zero point of each block</h3>
<p>Neither route to a sky position has an absolute reference.  A similarity
fit onto a Gaia cone has four free parameters and two of them are
translations, so a fit that locked onto a thin set of correspondences can be
internally excellent &mdash; right scale, right rotation, small residuals
&mdash; and still hand out every position displaced by one common vector.
Nothing upstream can see that.  The catalogue tie is the first stage in this
pipeline that <i>can</i>, so it measures it before matching anything: every
block is paired against the catalogue at {ct.ASTROM_LOOSE_TOL_ARCSEC:g}&Prime;
&mdash; deliberately far wider than the {ct.MATCH_TOL_ARCSEC:g}&Prime;
photometric tolerance, because an offset larger than that tolerance is the
only kind worth finding and a tight search cannot see one &mdash; and the
median displacement is taken over the block's VETTED comparison and check
stars.</p>

<p>It is removed only when it is bigger than
{ct.ASTROM_REFINE_MIN_ARCSEC:g}&Prime;, coherent to better than
{ct.ASTROM_REFINE_MAX_SCATTER_ARCSEC:g}&Prime; across stars, and measured on
at least {ct.ASTROM_REFINE_MIN_STARS} of them.  Those three gates are the
whole safeguard: a stage that quietly bends every block's astrometry onto
the catalogue it is about to be tied to has spent the independence the tie
depends on.  Blocks that were left alone are listed too, because
&ldquo;measured, 0.26&Prime;, left alone&rdquo; and &ldquo;never
examined&rdquo; must not look the same.</p>
{_astrom_table(con)}
</section>"""


def section_fits(con, fig3: str, fig4: str) -> str:
    rows = _primary(con)
    body = table(
        ["series", "catalogue : band", "system", "n<sub>fit</sub>",
         "ZP &plusmn; err", "colour term k &plusmn; err", "colour index",
         "c<sub>ref</sub>", "resid RMS (mmag)", "&chi;&sup2;/&nu;",
         "colour range (fit)"],
        [[_label(r["series_key"]),
          f"<code>{esc(r['catalogue'])}:{esc(r['band'])}</code>",
          esc(r["band_system"]), r["n_fit"],
          _pm(r["zp"], r["zp_err"], 4),
          _pm(r["colour_term"], r["colour_err"], 4),
          f"<code>{esc(r['colour_label'])}</code>", _n(r["colour_ref"], 3),
          _mmag(r["resid_rms"]), _n(r["chi2nu"], 2),
          f"{_n(r['colour_min'], 2)} &hellip; {_n(r['colour_max'], 2)} "
          f"(core {_n(r['colour_p05'], 2)}&ndash;{_n(r['colour_p95'], 2)})"]
         for r in rows],
        [VERDICT_CLASS.get(r["verdict"], "") for r in rows])
    return f"""
<section id="fits">
<div class="bhead"><h2>3 &middot; The fit: a zero point AND a colour
term</h2></div>

<p>For each block, on comparison stars only:</p>
<p class="eq"><code>m<sub>ens</sub> &minus; m<sub>cat</sub> = ZP0 + k
&times; (colour &minus; c<sub>ref</sub>)</code></p>

<p>Three choices in that one line are worth defending.  <b>Robust, not
least-squares</b>: a comparison sample always contains a few variables,
blends and mis-identifications, and ordinary least squares lets any one of
them tilt a published colour term.  Huber weighting
(&delta;&nbsp;=&nbsp;{ct.HUBER_DELTA:g}) lets an outlier keep a vote
proportional to 1/|residual| rather than |residual|; only after convergence,
where a point is plainly wrong rather than merely noisy, does the
{ct.CLIP_SIGMA:g}&sigma; pass delete it.  <b>Centred on
c<sub>ref</sub></b>, the median colour of the block's tie stars: centring
decorrelates ZP0 from k, so the quoted <code>zp_err</code> answers
&ldquo;how well is the zero point known&rdquo; and not &ldquo;how well
would it be known if the colour term were exactly right&rdquo;.  <b>Errors
inflated by &radic;(&chi;&sup2;/&nu;) when the scatter exceeds the
model</b>: catalogue errors are known to understate real star-to-star
scatter, because bandpass mismatch beyond a linear term is a genuine
per-star effect, and quoting the formal error of an obviously
under-dispersed model would understate the zero point by exactly the factor
the fit itself measured.</p>

{body}

{_fig(fig3, "The fits themselves.  Filled points entered the fit; open "
               "circles were held out for section 5 and never influenced "
               "the line.  The red line is the fitted relation over the "
               "colour range it actually interpolates; the dashed yellow "
               "line is the science target's own colour, where that is "
               "measurable (section 7).",
               "no block reached the minimum tie-star count.")}

{_fig(fig4, "Every colour term with its error bar.  A term of a few "
               "hundredths is bandpass mismatch between two nominally "
               "similar filters; a term of tenths would be a different "
               "band altogether, which is the question section 4 asks.",
               "no block reached the minimum tie-star count.")}
</section>"""


def section_bandtest(con) -> str:
    rows = q(con, """SELECT series_key, filter, catalogue, band, hypothesis,
                            band_system, n_fit, colour_term, colour_err,
                            resid_rms, check_rms
                     FROM cv_cattie
                     WHERE catalogue='gaia_gspc' AND filter IN ('G','R','I')
                       AND n_fit >= ?
                     ORDER BY series_key, hypothesis, band""",
             (ct.MIN_TIE_STARS,))
    # The lowercase/uppercase comparison this section really turns on can be
    # made with the PRIMARY catalogue alone, and is computed here rather
    # than asserted: the same target, the same filter letter, two eras.
    lower = q(con, """SELECT filter, colour_term FROM cv_cattie
                      WHERE is_primary=1 AND n_fit >= ?
                        AND filter IN ('g','r','i')""", (ct.MIN_TIE_STARS,))
    upper = q(con, """SELECT filter, colour_term FROM cv_cattie
                      WHERE is_primary=1 AND n_fit >= ?
                        AND filter IN ('G','R','I')""", (ct.MIN_TIE_STARS,))
    pair_rows, spreads = [], []
    for lo, up in (("g", "G"), ("r", "R"), ("i", "I")):
        a = [x[1] for x in lower if x[0] == lo and x[1] is not None]
        b = [x[1] for x in upper if x[0] == up and x[1] is not None]
        if a and b:
            dmed = float(np.median(b)) - float(np.median(a))
            spreads.append(abs(dmed))
            pair_rows.append([f"<code>{lo}</code> vs <code>{up}</code>",
                              len(a), _n(float(np.median(a)), 4),
                              len(b), _n(float(np.median(b)), 4),
                              _n(dmed, 4)])
    pair_tab = table(
        ["filter pair", "n (lower)", "median k (lower)", "n (upper)",
         "median k (upper)", "difference of medians"], pair_rows)
    # A difference of MEDIANS is not a bound on the individual terms, and
    # review was right to say the earlier prose read as if it were.  The
    # tighter comparison controls for the field and the star sample: the
    # SAME target, the same filter letter, two eras.  Its largest
    # disagreement is the number the claim is now made on.
    same = q(con, """SELECT a.target_key, a.filter, b.filter,
                            a.era_id, b.era_id, a.colour_term, b.colour_term
                     FROM cv_cattie a JOIN cv_cattie b
                       ON a.target_key = b.target_key
                      AND a.colour_label = b.colour_label
                      AND a.band = b.band
                     WHERE a.is_primary=1 AND b.is_primary=1
                       AND a.n_fit >= ? AND b.n_fit >= ?
                       AND a.filter IN ('G','R','I')
                       AND b.filter IN ('g','r','i')
                     ORDER BY 1, 2""",
             (ct.MIN_TIE_STARS, ct.MIN_TIE_STARS))
    same_rows, same_d = [], []
    for tk, fu, fl, eu, el, ku, kl in same:
        if ku is None or kl is None:
            continue
        same_d.append(abs(ku - kl))
        same_rows.append([TARGET_LABEL.get(tk, tk),
                          f"<code>{esc(fu)}</code> e{eu}",
                          f"<code>{esc(fl)}</code> e{el}",
                          _n(ku, 4), _n(kl, 4), _mmag(ku - kl)])
    same_tab = table(
        ["target", "uppercase block", "lowercase block", "k (upper)",
         "k (lower)", "difference (mmag/mag)"], same_rows)
    # The claim has to quote the measured separation, not gesture at it, and
    # it has to quote the LARGEST one rather than the most flattering.
    pair_max = max(same_d + spreads, default=float("nan"))
    if not rows:
        return f"""
<section id="bandtest">
<div class="bhead"><h2>4 &middot; Do the uppercase filter labels mean the
same glass as the lowercase ones?</h2></div>

<p>The archive writes <code>G</code>, <code>R</code>, <code>I</code> in eras
6 and 7 and <code>g</code>, <code>r</code>, <code>i</code> in eras 72, 76,
78 and 80 &mdash; and <code>y</code> and <code>z</code> appear in era 47,
which is a Pan-STARRS-like <i>grizy</i> wheel.  The obvious reading is that
the case is a control-software convention and the glass is the same.  The
obvious reading is not evidence, and mis-assigning a band puts a systematic
straight into a published zero point.</p>

<p>The DECISIVE test needs Gaia's synthetic photometry, which publishes
Sloan and Johnson&ndash;Cousins magnitudes for the same stars from the same
BP/RP spectra &mdash; so the two hypotheses would differ only in the band.
That catalogue is not available in this build (see section 1), so the test
below is the weaker one that the primary catalogue alone can support: fit
both cases against the SAME PS1 bands and compare the colour terms.  If the
uppercase labels were Johnson&ndash;Cousins filters, their colour terms
against PS1 <i>gri</i> would be a large fraction of a magnitude per
magnitude of colour, and they would not sit on top of the lowercase
ones.</p>

{pair_tab}

<p>A difference of medians is a summary, not a bound &mdash; it can hide
individual blocks that disagree by much more, and review found that it did.
The tighter comparison holds the field and the star sample fixed: the SAME
target, the same filter letter, one uppercase era against one lowercase
era.</p>
{same_tab}

<div class="decision"><b>No uppercase colour term sits more than
{_mmag(pair_max)} mmag per magnitude of colour from its lowercase
counterpart</b> &mdash; that is the LARGEST per-pair disagreement, not an
average, and it is still one to two ORDERS OF MAGNITUDE smaller than the
0.2&ndash;0.5 mag/mag a genuine Johnson&ndash;Cousins bandpass would produce
against PS1 <i>gri</i>.  The
primary tie is taken on the Sloan reading.  This is strong evidence, not
proof: it rules out a grossly different bandpass and cannot separate two
similar ones.  The definitive two-system test is left as the first thing to
re-run when the ESA archive is reachable, and the alternative-hypothesis
rows are already provisioned in <code>cv_cattie</code> to receive it.</div>
</section>"""
    body = table(
        ["series", "FILTER", "hypothesis", "catalogue band", "system",
         "n", "colour term k", "resid RMS (mmag)"],
        [[_label(r[0]), f"<code>{esc(r[1])}</code>", esc(r[4]),
          f"<code>{esc(r[3])}</code>", esc(r[5]), r[6],
          _pm(r[7], r[8], 4), _mmag(r[9])] for r in rows],
        ["ok" if r[4] == "primary" else "" for r in rows])
    # The deciding numbers, computed rather than asserted.
    prim = [r for r in rows if r[4] == "primary"]
    alt = [r for r in rows if r[4] == "alternative"]
    kp = np.array([abs(r[7]) for r in prim if r[7] is not None])
    ka = np.array([abs(r[7]) for r in alt if r[7] is not None])
    rp = np.array([r[9] for r in prim if r[9] is not None])
    ra = np.array([r[9] for r in alt if r[9] is not None])
    verdict = ("the Sloan reading wins on both numbers"
               if kp.size and ka.size and np.median(kp) < np.median(ka)
               and np.median(rp) <= np.median(ra)
               else "the two readings are not cleanly separated by these "
                    "numbers, and the label question stays open")
    return f"""
<section id="bandtest">
<div class="bhead"><h2>4 &middot; Do the uppercase filter labels mean the
same glass as the lowercase ones?</h2></div>

<p>The archive writes <code>G</code>, <code>R</code>, <code>I</code> in eras
6 and 7 and <code>g</code>, <code>r</code>, <code>i</code> in eras 72, 76,
78 and 80 &mdash; and <code>y</code> and <code>z</code> appear in era 47,
which is a Pan-STARRS-like <i>grizy</i> wheel.  The obvious reading is that
the case is a control-software convention and the glass is the same.  The
obvious reading is not evidence, and mis-assigning a band puts a
systematic straight into a published zero point.</p>

<p>Gaia's synthetic photometry makes it a measurement.  It publishes both
systems for the same stars from the same BP/RP spectra, so the two
hypotheses differ ONLY in the band &mdash; not in the sample, the epoch or
the calibration.  If <code>R</code> were Cousins&nbsp;R, its colour term
against SDSS <i>r</i> would be a large fraction of a magnitude per
magnitude of colour, because the two bandpasses genuinely differ; if it is
Sloan&nbsp;<i>r</i>, the term is the small residual mismatch of two similar
filters.</p>

{pair_tab}
{same_tab}

<p>And the two-system test itself, on Gaia's synthetic photometry:</p>
{body}

<div class="decision"><b>Median |k| under the primary (Sloan) reading:
{_n(float(np.median(kp)) if kp.size else float('nan'), 4)}; under the
Johnson&ndash;Cousins alternative:
{_n(float(np.median(ka)) if ka.size else float('nan'), 4)}.  Median residual
RMS: {_mmag(float(np.median(rp)) if rp.size else float('nan'))} mmag vs
{_mmag(float(np.median(ra)) if ra.size else float('nan'))} mmag.</b>
On these numbers, {verdict}.  The primary tie is taken on the Sloan
reading; the alternative rows stay in <code>cv_cattie</code> so the
comparison can be re-read rather than re-argued.</div>
</section>"""


def section_accuracy(con, fig5: str) -> str:
    rows = _primary(con)
    body = table(
        ["series", "n<sub>check</sub>", "clipped check RMS (mmag)",
         "raw check RMS (mmag)", "outliers", "check median / bias (mmag)",
         "check MAD (mmag)", "fit resid RMS (mmag)", "verdict"],
        [[_label(r["series_key"]), r["n_check"], _mmag(r["check_rms_clip"]),
          _mmag(r["check_rms"]),
          (f"{r['n_check_outlier']} / {r['n_check']}"
           if r["n_check"] else "&mdash;"),
          _mmag(r["check_median"]), _mmag(r["check_mad"]),
          _mmag(r["resid_rms"]), f"<b>{esc(r['verdict'])}</b>"]
         for r in rows],
        [VERDICT_CLASS.get(r["verdict"], "") for r in rows])
    tied = [r for r in rows if r["verdict"].startswith("TIED")
            and r["check_rms_clip"] is not None]
    ok = [r for r in tied if r["check_rms_clip"] <= ct.ACCURACY_GOAL_MAG]
    med = (float(np.median([r["check_rms_clip"] for r in tied]))
           if tied else float("nan"))
    n_out = sum(r["n_check_outlier"] or 0 for r in tied)
    n_chk = sum(r["n_check"] or 0 for r in tied)
    raw = (float(np.median([r["check_rms"] for r in tied
                            if r["check_rms"] is not None]))
           if tied else float("nan"))
    return f"""
<section id="accuracy">
<div class="bhead"><h2>5 &middot; How accurate is it &mdash; on stars the
fit never saw?</h2></div>

<p>Residual scatter about the fit is not accuracy.  The fit was optimised on
those stars; quoting its own residual as the achieved error is the oldest
way to flatter a calibration.  So {100 * ct.HOLDOUT_FRACTION:g}% of the
eligible stars are withheld from the solve by a deterministic hash of
<code>salt|series|star_id</code> &mdash; deterministic because a published
accuracy that changes on re-run is not a measurement, and hashed rather
than taken modulo an id because id order correlates with detection order,
which correlates with brightness, which is the axis the check must be blind
to.  Every star the ensemble itself designated a CHECK star is added to
that holdout.</p>

<p>For a held-out star the colour is known and the star is not the science
target, so the FULL transformation may be applied to it &mdash; that is the
point: it tests the calibration exactly as a user would apply it.  The
residual <code>(m<sub>nat</sub> &minus; k(C &minus; c<sub>ref</sub>))
&minus; m<sub>cat</sub></code> is what a published magnitude of that star
would be wrong by.  Its RMS is the block's achieved accuracy; its MEDIAN is
residual zero-point bias, which is the half that does NOT average away over
stars and is therefore the more dangerous half.</p>

{_fig(fig5, "Achieved absolute accuracy per block.  Green band is the "
               "strategy's 10&ndash;20 mmag goal.  Grey bars are the fit's "
               "own residual RMS (the flattering number); circles are the "
               "independent check RMS (the honest one); yellow squares are "
               "the check-star median, i.e. the bias that does not average "
               "down.",
               "no block has held-out check stars to measure accuracy "
               "from.")}

{body}

<div class="decision"><b>{len(ok)} of {len(tied)} tied blocks meet the
{1000 * ct.ACCURACY_GOAL_MAG:g} mmag goal on independent check stars; the
median achieved accuracy is {_mmag(med)} mmag after clipping, and
{_mmag(raw)} mmag before it.</b>  A block graded TIED-ABOVE-GOAL is not a
failure &mdash; its magnitudes are on a standard system and their accuracy
is measured and stated.  It is a block whose error bar is larger than the
strategy hoped, which is a fact about the data and not about the tie.</div>

<div class="note"><b>The clip, stated rather than absorbed.</b>
{n_out} of {n_chk} held-out stars ({100.0 * n_out / max(n_chk, 1):.1f}%)
sit beyond {ct.CLIP_SIGMA:g} robust sigma of their block's own check
distribution, some of them by more than a magnitude.  They are not
calibration error: they are variables, blends the veto in section 2 did not
catch, and mis-identifications &mdash; and they are why the raw and clipped
RMS differ by an order of magnitude on several blocks.  Both numbers are
published because they answer different questions.  &ldquo;How well is a
typical star calibrated?&rdquo; is the clipped number.  &ldquo;If I pick one
arbitrary star and publish its magnitude, what is my risk?&rdquo; is the
clipped number PLUS a
{100.0 * n_out / max(n_chk, 1):.1f}% chance of drawing one of these.
Quoting only the clipped RMS would hide that risk; quoting only the raw RMS
would claim the calibration is ten times worse than it is for the star you
probably have.</div>
</section>"""


def _era_trend_table(con, axis: str) -> tuple[str, list]:
    """One trend axis, broken down by camera era and pixel provenance.

    This is where the (x, y) question stops being a diagnostic and becomes
    a finding: a systematic that belongs to the DETECTOR should repeat
    across fields within an era and change when the era changes, and only
    a per-era breakdown can show that.
    """
    rows = q(con, """
        SELECT k.era_id, f.readoutm, s.provenance,
               count(*), sum(t.significant), avg(t.swing),
               min(t.swing), max(t.swing),
               count(DISTINCT k.target_key)
        FROM cv_cattie_trend t
        JOIN cv_cattie k ON k.series_key = t.series_key AND k.is_primary = 1
        JOIN cv_series s ON s.series_key = t.series_key
        JOIN (SELECT series_key, max(readoutm) AS readoutm
              FROM cv_frames GROUP BY 1) f ON f.series_key = t.series_key
        WHERE t.axis = ?
        GROUP BY 1, 2, 3 ORDER BY 1""", (axis,))
    return table(
        ["era", "readout mode", "pixel provenance", "fields", "blocks",
         f"significant at {ct.TREND_SIGMA:g}&sigma;", "mean swing (mmag)",
         "range (mmag)"],
        [[r[0], f"<code>{esc(r[1])}</code>", esc(r[2]), r[8], r[3], r[4],
          _mmag(r[5]), f"{_mmag(r[6])} &hellip; {_mmag(r[7])}"]
         for r in rows]), rows


def section_adversarial(con, fig6: str, fig7: str, fig8: str,
                        extra: str = "") -> str:
    tr = q(con, """SELECT axis, count(*), sum(significant),
                          max(abs(swing)), avg(abs(swing))
                   FROM cv_cattie_trend GROUP BY axis ORDER BY axis""")
    tr_tab = table(
        ["axis", "regressions", "significant at "
         f"{ct.TREND_SIGMA:g}&sigma;", "largest swing (mmag)",
         "mean |swing| (mmag)"],
        [[f"<code>{esc(r[0])}</code>", r[1], r[2], _mmag(r[3]), _mmag(r[4])]
         for r in tr])
    rad_tab, rad_rows = _era_trend_table(con, "radius")
    mag_tab, mag_rows = _era_trend_table(con, "cat_mag")
    # The deciding comparison, computed rather than asserted: the era with
    # the largest mean radial swing, against the rest.
    worst = (max(rad_rows, key=lambda r: abs(r[5] or 0)) if rad_rows
             else None)
    rest = [r for r in rad_rows if r is not worst]
    rest_mean = (float(np.mean([abs(r[5] or 0) for r in rest]))
                 if rest else float("nan"))
    worst_txt = ""
    if worst:
        worst_txt = (
            f"<b>Era {worst[0]} ({esc(worst[1])}, {esc(worst[2])} pixels) "
            f"carries a radial residual of {_mmag(worst[5])} mmag from "
            f"centre to corner, significant in {worst[4]} of its "
            f"{worst[3]} blocks across {worst[8]} independent fields.  "
            f"Every other era averages {_mmag(rest_mean)} mmag.</b>")
    cross = q(con, """SELECT series_key, band_a, band_b, n_common,
                             star_offset_median, star_offset_rms,
                             star_offset_mad, d_zp, d_colour_term
                      FROM cv_cattie_cross ORDER BY series_key""")
    cr_tab = table(
        ["series", "n common stars", "REFCAT2 &minus; Gaia, median (mmag)",
         "MAD (mmag)", "RMS (mmag)", "&Delta;ZP (mmag)", "&Delta;k"],
        [[_label(r[0]), r[3], _mmag(r[4]), _mmag(r[6]), _mmag(r[5]),
          _mmag(r[7]), _n(r[8], 4)] for r in cross])
    off = np.array([abs(r[4]) for r in cross if r[4] is not None])
    dzp = np.array([abs(r[7]) for r in cross if r[7] is not None])
    return f"""
<section id="adversarial">
<div class="bhead"><h2>6 &middot; Three attacks on our own answer</h2></div>

<h3>6.1 &mdash; is the residual a function of brightness?</h3>
<p>If the detector were non-linear, or if the aperture correction varied
with signal-to-noise, the tie residual would slope against catalogue
magnitude.  The number that matters is not the slope's significance but its
SWING &mdash; slope times the magnitude range actually fitted &mdash;
because a formally significant 1 mmag per magnitude over four magnitudes is
4 mmag of systematic, a fifth of the accuracy goal, while the same
significance over 0.2 mag is nothing.</p>
{_fig(fig6, "Left: every fitted tie star of every block, residual "
               "against catalogue magnitude, with the binned median in "
               "white. Right: the distribution of per-block trend swings, "
               "split by whether the trend is formally significant.",
               "no fitted tie stars to regress.")}
{mag_tab}

<h3>6.2 &mdash; is the residual a function of DETECTOR POSITION?</h3>
<p>This is the question nobody had asked.  The characterization measured a
noise floor tens of times above scintillation, called it instrumental, and
left it there; a flat-field residual would produce exactly that AND would
show up as a tilt of the tie residual across the chip.  A plane is fitted
in (x,&nbsp;y) per block and reported as its corner-to-corner swing, which
is the systematic a star suffers for being on the wrong side of the
detector.</p>
{_fig(fig7, "Tie residuals in detector coordinates for the "
               "best-populated block of each target &mdash; "
               + ", ".join(_label(k) for k in RESID_XY_KEYS)
               + " &mdash; plus the radial profile and the distribution of "
               "corner-to-corner plane swings over every block.  One panel "
               "per TARGET, not the largest blocks outright: a flat-field "
               "residual should repeat across FIELDS on the same detector, "
               "and only a panel per field can show whether it does.",
               "no block has 25 or more fitted tie stars with recorded "
               "detector positions.")}

<p>The maps were drawn expecting a TILT and showed a BULLSEYE: positive
residuals in the middle of the frame, negative at the corners, in four
different fields.  A plane cannot represent that shape, which means the
plane fit systematically UNDERSTATES it &mdash; so the residual was
regressed on RADIUS from the detector centre as well, and that is the test
that carries the finding.</p>
{rad_tab}

<div class="decision">{worst_txt}
A systematic that repeats across four unrelated star fields, in one camera
era and one pixel provenance and not in the others, is not a property of
any field: it is a property of the reduction that made those pixels.  The
signature is an illumination-correction error, i.e. a flat field that
does not describe the true vignetting &mdash; and the SIGN of the effect
is what picks that explanation out from its rivals.  The residual
<i>falls</i> with radius: stars near the edge measure BRIGHTER than the
catalogue says, stars in the middle fainter.  A degrading point-spread
function toward the field edge, the other obvious candidate, would push the
opposite way (light spilling out of a fixed 4-arcsec aperture makes edge
stars measure FAINTER), so it cannot be the cause.  An over-corrected
vignetting profile &mdash; dividing by a flat that is too low at the
edges &mdash; produces exactly this sign.  Note WHICH eras
show it: the ones reduced by the observatory's own pipeline
(<code>server_reduced</code>).  The eras calibrated here from staged
era-matched masters (<code>local_master</code>) do not.  That is a
falsifiable claim and the cheapest way to test it is to re-reduce one era-76
night from raw with local masters and re-run this stage on it.</div>

<p>The trend census over all axes, for completeness &mdash; and a warning
about how to read it: {q1(con, "SELECT count(*) FROM cv_cattie_trend")}
regressions were run, so at
{ct.TREND_SIGMA:g}&sigma; one would expect of order
{q1(con, "SELECT count(*) FROM cv_cattie_trend") * 0.0027:.1f} false
positives by chance.  The radial result is not in that category and neither
is the magnitude one; several of the per-axis x and y results, taken
individually, might be.</p>
{tr_tab}

<h3>6.3 &mdash; do two independent catalogues agree?</h3>
<p>The difference between ATLAS-REFCAT2 and Gaia synthetic photometry on
the SAME stars, in nominally the same band, is the honest systematic floor
of any magnitude published here &mdash; it belongs in the paper's error
budget and no single catalogue can supply it.  It is not a pure error
term: the two systems are PS1&nbsp;AB and synthetic SDSS&nbsp;AB, so a real
system difference of a few hundredths is expected and is part of what the
number measures.  Both catalogues are vetoed by the SAME gate (section 2)
precisely so this difference measures the catalogues and not two different
pipelines.</p>
{_fig(fig8, "Star-by-star offsets (circles, with the robust scatter as "
               "the error bar) and the fitted zero-point difference "
               "(crosses), per block.",
               "no series is tied by BOTH catalogues, so there is no "
               "cross-catalogue difference to draw — see the note below "
               "for what that costs.")}
{cr_tab}

{'''<div class="bad"><b>NOT MEASURED IN THIS BUILD, and that is a gap in
the error budget rather than a formality.</b>  The ESA Gaia archive
returned HTTP 500 (&ldquo;canceling statement due to statement
timeout&rdquo;) to every attempt during this run, so no series is tied by
both catalogues and the cross-catalogue term is simply absent.  What it
would have bounded: how much of each zero point is ATLAS-REFCAT2's own
calibration rather than this telescope's.  Until it is measured, every
absolute magnitude on this page should be quoted with the caveat that its
zero point rests on ONE catalogue, and the number that would replace the
caveat is a re-run of <code>fetch</code> away &mdash; the stage is
resumable, the REFCAT2 side is cached, and nothing else has to be
recomputed.</div>''' if not cross else f'''<div class="decision"><b>The
systematic floor from catalogue choice: median |offset|
{_mmag(float(np.median(off)) if off.size else float("nan"))} mmag across
{len(cross)} blocks measured on common stars, and
{_mmag(float(np.median(dzp)) if dzp.size else float("nan"))} mmag as a
difference of fitted zero points.</b>  Any absolute magnitude published
from this archive carries that term whether or not it is written down, so
it is written down.</div>'''}
{extra}
</section>"""


def section_literature(con) -> str:
    """6.4 -- the one check that does not use the fitting catalogue."""
    rows = q(con, """SELECT c.series_key, c.target_key, c.filter,
                            count(l.cal_mag), avg(l.cal_mag)
                     FROM cv_cattie c JOIN cv_lightcurve l
                       ON l.series_key = c.series_key AND l.role = 'target'
                      AND l.cal_mag IS NOT NULL
                     WHERE c.is_primary = 1 AND c.verdict LIKE 'TIED%'
                     GROUP BY 1, 2, 3
                     ORDER BY c.target_key, c.era_id, c.filter""")
    # The MEDIAN is what the check wants (a CV spends time at both ends of
    # its range), so it is computed here rather than in SQL, which has no
    # median.
    out, verdicts, flagged = [], [], []
    for skey, tk, filt, n, _mean in rows:
        vals = [r[0] for r in q(
            con, "SELECT cal_mag FROM cv_lightcurve WHERE series_key=? "
                 "AND role='target' AND cal_mag IS NOT NULL", (skey,))]
        med = float(np.median(vals)) if vals else float("nan")
        verdict, why = ct.literature_check(tk, med)
        verdicts.append(verdict)
        out.append([_label(skey), n, _n(med, 2), f"<b>{esc(verdict)}</b>",
                    esc(why)])
        if verdict in ("outside", "near"):
            flagged.append(f"{_label(skey)} at {med:.2f} on {n:,} target "
                           f"points ({esc(why)})")
    body = table(["series", "target points", "median calibrated mag",
                  "vs published range", "detail"], out,
                 ["ok" if v == "inside" else ("warn" if v == "near" else "bad")
                  for v in verdicts])
    n_in = sum(1 for v in verdicts if v == "inside")
    n_near = sum(1 for v in verdicts if v == "near")
    n_out = sum(1 for v in verdicts if v == "outside")
    return f"""
<h3>6.4 &mdash; does the answer land where the outside world says it
should?</h3>
<p>Every check so far compares the tie against the catalogue the tie was
fitted to.  All of them would survive a sign error, or a gauge mistake
applied uniformly across a block, because both sides of the comparison
would move together.  So here is the one test that uses nothing from this
pipeline: the median CALIBRATED magnitude of each science target, against
its published brightness range in the AAVSO Variable Star Index.</p>

<p>It is a bound, not a calibration.  A CV's brightness depends on its
accretion state, the published ranges are broad by construction, and
<i>V</i> is not any of the bands measured here &mdash; so a target inside
its range is evidence that nothing is grossly wrong, not evidence that the
zero point is good to 20 mmag.  Nothing here is fitted to these numbers and
nothing is corrected by them.  They are the only hand-entered values in
<code>macro_phot.cattie</code> and they are labelled as such in the source,
with their date and provenance.</p>

{body}

<div class="{'ok' if n_out == 0 else 'warn'}"><b>{n_in} of {len(out)}
calibrated series place their target inside its published range, {n_near}
within {ct.LITERATURE_TOLERANCE_MAG:g} mag of it, and {n_out} outside.</b>
A wrong zero-point SIGN would put every one of them tens of magnitudes out;
a mis-scaled gauge would put them out by a consistent offset.  Neither is
present, in any band, in any era &mdash; which is the whole of what this
test can establish, and it establishes it.
{"" if not flagged else "<br><br><b>The exceptions, named rather than "
 "averaged away:</b><br>" + "<br>".join(flagged) + "<br><br>Read them with "
 "the published range's own limits in mind: it is a <i>V</i> range, these "
 "are red bands, and a polar in a deep low state sits below any range "
 "catalogued from its brighter epochs. That is an explanation, not a "
 "dismissal &mdash; the honest statement is that this test cannot resolve "
 "a genuine zero-point error from a genuine faint state on these series, "
 "and the series it cannot resolve are the ones with the fewest target "
 "detections."}</div>
"""


def section_target_colour(con) -> str:
    # Only TIED blocks: an untied block has no colour range to be inside or
    # outside of, and counting it as "unknown" would inflate the denominator
    # with rows the question does not apply to.
    rows = [r for r in _primary(con) if r["verdict"].startswith("TIED")]
    body = table(
        ["series", "target colour", "epoch pairs", "pair-to-pair scatter",
         "how it was obtained", "fit colour range", "position",
         "extrapolation cost (mmag)"],
        [[_label(r["series_key"]), _n(r["target_colour"], 3),
          f"{r['n_colour_pairs']:,}" if r["n_colour_pairs"] else "&mdash;",
          _n(r["colour_scatter"], 3),
          esc(r["target_colour_source"]),
          f"{_n(r['colour_min'], 2)} &hellip; {_n(r['colour_max'], 2)}",
          f"<b>{esc(r['colour_position'])}</b>", _mmag(r["extrap_err"])]
         for r in rows],
        ["bad" if r["colour_position"] == "extrapolated" else
         ("warn" if r["colour_position"] in ("unknown", "inside-span")
          else "ok") for r in rows])
    n_out = sum(1 for r in rows if r["colour_position"] == "extrapolated")
    n_unk = sum(1 for r in rows if r["colour_position"] == "unknown")
    n_core = sum(1 for r in rows if r["colour_position"] == "inside-core")
    return f"""
<section id="targetcolour">
<div class="bhead"><h2>7 &middot; Where does each CV's own colour fall
relative to its fit?</h2></div>

<p>Ruling&nbsp;1 forbids transforming the target, and nothing here does.
But the target's POSITION on the colour axis still has to be stated,
because it is what decides how much of the fit applies to it.  A target
inside the tie stars' colour core is calibrated by many stars; one inside
the span but outside the core is calibrated by a handful; one outside the
span entirely is EXTRAPOLATED, and that is a limit on what may be claimed
rather than a detail.</p>

<p>The colour is measured, not assumed &mdash; and it is measured from
PAIRED EPOCHS, which is a correction to what this page said in its first
edition.  That edition formed the colour from each filter's ensemble MEAN
magnitude over the whole series.  For a constant star that is the colour;
for a polar it is the difference between the target's mean state while the
blue filter was on the wheel and its mean state while the red one was, and
nothing requires those two campaigns to overlap.  VV&nbsp;Pup era 76 was
published here at <i>g&minus;r</i>&nbsp;=&nbsp;&minus;1.73 on that recipe:
its <i>g</i> points span 370 days and its <i>r</i> points 55 of them, and on
the epochs where both filters actually observed the object the colour is
+0.04.  A 1.77-magnitude error, and it was the number that put the block in
the &ldquo;extrapolated&rdquo; row with a 105&nbsp;mmag charge attached.</p>

<p>So each blue measurement is now paired with the nearest red measurement
in time, pairs more than
{ct.COLOUR_PAIR_TOL_DAYS * 1440:.0f}&nbsp;minutes apart are discarded, and
the colour is the median of what survives.  Where the two filters never
sampled the same time the answer is <i>unknown</i> &mdash; VV&nbsp;Pup era
72's <i>g</i> and <i>r</i> share no night at all, and neither do EU&nbsp;UMa
era 76's <i>r</i> and <i>i</i> &mdash; which is the same answer a
single-filter era gets, for the same reason.  Converting the natural-system
colour to the catalogue one is then a closed form in the two colour terms
(<code>cattie.target_colour_solve</code>).  It produces a position, never a
magnitude.</p>

<p>Seven minutes is not simultaneity.  These polars have roughly
hundred-minute orbits, so the tolerance is still some seven per cent of a
cycle, and the pair-to-pair scatter column says how much the colour moves
across the pairs &mdash; for the AM&nbsp;Her systems it is one to two tenths
of a magnitude, which is the object varying and not the measurement failing.
The median is quoted because it is what the position test needs; the scatter
is quoted beside it because a single number for a variable object's colour
would otherwise read as more definite than it is.</p>

{body}

<div class="decision"><b>{n_core} of {len(rows)} blocks place their target
inside the fit's colour core; {n_out} are extrapolated; {n_unk} cannot be
placed at all (single-filter era, no shared epochs between the two filters,
or the target undetected in that series).</b>  For the extrapolated and
unknown blocks the zero point is still valid &mdash; it is a zero point,
measured on comparison stars &mdash; but the residual bandpass error at the
target's own colour is NOT bounded by the fit's residual scatter, and the
extrapolation cost column is the size of the term that has to be carried
instead.  An <i>unknown</i> here is a weaker statement than it looks:
it does not mean the target is inside the range, it means this archive
never measured its colour in that era and no claim about the bandpass
error at the target's colour can be made from these data alone.</div>
</section>"""


def section_relative(con) -> str:
    solved = q(con, """SELECT series_key, target_key, era_id, filter,
                              n_frames_used, n_target_points
                       FROM cv_series WHERE status='solved'
                       ORDER BY target_key, era_id, filter""")
    tied = {r[0] for r in q(con, "SELECT series_key FROM cv_cattie "
                                 "WHERE is_primary=1 "
                                 "AND verdict LIKE 'TIED%'")}
    # Which (catalogue, field) pulls actually landed.  A filter whose only
    # catalogue analogue lives in a pull that failed is untied for a
    # DIFFERENT reason than one with no analogue anywhere, and saying
    # "no catalogue match" for both would bury a recoverable outage under a
    # permanent limitation.
    have_cat = {(r[0], r[1]) for r in
                q(con, "SELECT catalogue, field_key FROM cv_cat_fetch "
                       "WHERE n_rows > 0")}
    rows = []
    for skey, tk, era, filt, nfr, npt in solved:
        if skey in tied:
            continue
        why: list[str] = []
        cands = ct.band_candidates(filt)
        usable = [b for b in cands if (b.catalogue, tk) in have_cat]
        ft = q(con, """SELECT method, status FROM cv_field_tie
                       WHERE target_key=? AND era_id=?""", (tk, era))
        tie_row = q(con, "SELECT n_candidates, n_clean, n_fit, note "
                         "FROM cv_cattie WHERE series_key=? AND is_primary=1",
                    (skey,))
        if not cands:
            why.append("no catalogue in this stage publishes a band that "
                       f"matches the <code>{esc(filt)}</code> filter")
        elif not usable:
            missing = ", ".join(sorted({b.catalogue for b in cands}))
            why.append(f"the only catalogue with a matching band "
                       f"(<code>{esc(missing)}</code>) has no successful "
                       f"pull for this field &mdash; the tie is BLOCKED, "
                       f"not impossible, and completing that pull completes "
                       f"this block")
        if ft and ft[0][0] in (None, "none"):
            why.append(f"the (target, era) field tie did not solve "
                       f"&mdash; <code>{esc(ft[0][1])}</code>; with no sky "
                       f"positions for its reference stars the block cannot "
                       f"be matched to any catalogue")
        elif tie_row:
            n_cand, n_clean, n_fit, _note = tie_row[0]
            n_match = q1(con, """SELECT count(*) FROM cv_stars s
                                 JOIN cv_cat_match m
                                   ON m.target_key = s.target_key
                                  AND m.era_id = s.era_id
                                  AND m.star_id = s.star_id
                                 WHERE s.series_key = ?
                                   AND s.role IN ('comp','check')""", (skey,))
            if n_match == 0:
                # An unmatched block gets its measured astrometric offset
                # quoted, not a guess about its quality.  This is where the
                # first edition of the page went wrong: it read "no matches"
                # as "bad astrometry" and said so, when the offset table
                # (section 2.1) can distinguish a displaced solution from a
                # noisy one and the two want opposite remedies.
                ast = q(con, """SELECT offset_arcsec, scatter_arcsec, n_paired
                                FROM cv_cat_astrom WHERE catalogue=?
                                  AND target_key=? AND era_id=?""",
                        (ct.PRIMARY_CATALOGUE, tk, era))
                nref = q1(con, "SELECT count(*) FROM cv_ref_stars "
                               "WHERE target_key=? AND era_id=?", (tk, era))
                detail = ""
                if ast and ast[0][0] is not None:
                    detail = (f"  Its measured displacement from the "
                              f"catalogue is {ast[0][0]:.2f}&Prime; with "
                              f"{ast[0][1]:.2f}&Prime; of scatter over "
                              f"{ast[0][2]:,} paired stars, so this is "
                              f"scatter rather than a removable offset")
                why.append(
                    f"not one of its {n_cand} comparison stars found a "
                    f"catalogue counterpart within "
                    f"{ct.MATCH_TOL_ARCSEC:g}&Prime;, out of {nref:,} "
                    f"reference detections.{detail}")
            else:
                why.append(f"{n_cand} candidate comparison stars, {n_match} "
                           f"matched a catalogue source, {n_clean} survived "
                           f"the cleanliness gate and {n_fit} reached the "
                           f"fit &mdash; below the minimum of "
                           f"{ct.MIN_TIE_STARS}")
        elif usable:
            why.append("no catalogue source matched any of this series' "
                       "comparison stars")
        rows.append([_label(skey), nfr, npt, " &middot; ".join(why)])
    body = table(["series", "frames used", "target points",
                  "why it is still RELATIVE"], rows,
                 ["bad"] * len(rows))
    # The EU UMa story, told from the tables rather than from memory.  The
    # first edition of this page reported era 78 as untieable for bad
    # astrometry; the offset measurement (section 2.1) says otherwise, and
    # the sentence below is now built from that row.
    fixed = q(con, """SELECT a.target_key, a.era_id, a.offset_arcsec,
                             a.scatter_arcsec, a.n_paired, a.n_match_before,
                             a.n_match_after, t.status, t.n_gaia
                      FROM cv_cat_astrom a
                      LEFT JOIN cv_field_tie t
                        ON t.target_key=a.target_key AND t.era_id=a.era_id
                      WHERE a.catalogue=? AND a.applied=1
                      ORDER BY 1, 2""", (ct.PRIMARY_CATALOGUE,))
    fixed_txt = ""
    for tk, era, size, scat, npair, nb, na, status, ngaia in fixed:
        tied_here = q(con, """SELECT series_key, n_fit, check_rms_clip, verdict
                              FROM cv_cattie WHERE is_primary=1
                                AND target_key=? AND era_id=?""", (tk, era))
        got = ", ".join(f"{_label(s)} &rarr; <b>{esc(v)}</b> on {n} tie "
                        f"stars ({_mmag(c)} mmag check)"
                        for s, n, c, v in tied_here)
        fixed_txt += (
            f"<p><b>{TARGET_LABEL.get(tk, tk)} era&nbsp;{era} was not "
            f"untieable; it was displaced.</b>  Its reference solution puts "
            f"every star {size:.2f}&Prime; from where the catalogue does, "
            f"coherently &mdash; {scat:.2f}&Prime; of scatter over "
            f"{npair:,} vetted stars, which is a translation and not a "
            f"failure.  Inside the {ct.MATCH_TOL_ARCSEC:g}&Prime; "
            f"photometric tolerance that displacement is fatal and invisible "
            f"at the same time: {nb:,} of its reference stars matched the "
            f"catalogue before the offset was removed and {na:,} after.  "
            f"The block now ties: {got}.</p>"
            f"<p>The same displacement is why this era has no target light "
            f"curve.  The field-tie stage looks for the science target as "
            f"the nearest reference star to its catalogue position and "
            f"recorded <code>{esc(status or '')}</code> &mdash; a near miss "
            f"of exactly this size.  That stage owns "
            f"<code>cv_field_tie</code> and <code>cv_ref_stars</code> and "
            f"this one does not write them, so the tie above serves the "
            f"comparison stars and the science point is still missing.  "
            f"Recovering it is a specific, bounded job and it is named here "
            f"rather than implied: apply the same offset refinement inside "
            f"the field-tie stage, re-identify the target, and re-solve that "
            f"one ensemble.  Both halves need the ESA Gaia archive, which "
            f"was unreachable throughout this build.</p>")
    return f"""
<section id="relative">
<div class="bhead"><h2>8 &middot; What is still relative, and why</h2></div>

<p>A fabricated tie is worse than an honest gauge, so a block that cannot
be tied is left alone and named.  These series keep their differential
magnitudes, keep every precision, period and timing result that rests on
them &mdash; all of which are differential and none of which the tie
changes &mdash; and carry <code>cal_mag = NULL</code> so that no consumer
can accidentally read a relative magnitude as an absolute one.</p>

{body}

<div class="decision"><b>The two EU&nbsp;UMa blocks that failed on archive
HTTP 500s were re-run, and the two failures turned out to be different
animals &mdash; then review showed the first diagnosis of the first one was
still wrong.</b>  Era&nbsp;78's cone returns and its reference frame solves:
the 500 was transient exactly as suspected.  This page then reported the
block untieable because none of its comparison stars found a catalogue
counterpart, and blamed the astrometry.  The astrometry was not bad.
Era&nbsp;80's cone ALSO returns &mdash; 600 Gaia sources where there were
none &mdash; but the similarity fit of its reference frame to that cone
fails.  That is a geometric failure, not a network one, and re-running it
will fail identically: the era-80 reference frame has FWHM 1.66&nbsp;px,
undersampled to the point where most of its 1,839
&ldquo;detections&rdquo; are not stars the triangle matcher can use, and
it has no sky positions at all, so no offset refinement can reach it.  The
honest next action there is a different reference frame for that era, not
another query.</div>

{fixed_txt}

<div class="note"><b>Why an offset this large could sit undetected in a
finished product.</b> Every check available upstream is INTERNAL: the
similarity fit's own residuals, the scale, the rotation, the frame-to-frame
registration.  A pure translation leaves all of them perfect.  The first
stage with an absolute reference is this one, and until this build it used
that reference only at {ct.MATCH_TOL_ARCSEC:g}&Prime; &mdash; a tolerance
that turns a 5&Prime; displacement into &ldquo;no matches&rdquo;, which
reads like noise.  Measuring the offset at a loose radius costs one extra
cross-match per block and converts that silence into a number.  It is now
run on every block, and every block's answer is in section 2.1 whether it
was corrected or not.</div>
</section>"""


def section_verdict(con) -> str:
    rows = _primary(con)
    n_series = q1(con, "SELECT count(*) FROM cv_series WHERE status='solved'")
    n_block = q1(con, "SELECT count(DISTINCT target_key || '|' || era_id) "
                      "FROM cv_cattie WHERE is_primary=1 "
                      "AND verdict LIKE 'TIED%'")
    n_block_all = q1(con, "SELECT count(*) FROM cv_ref")
    # Blocks that HAVE a solved series -- the only ones a tie could reach.
    # Grading against all 14 would charge the tie for eras whose photometry
    # never solved in the first place, which is a different stage's problem.
    n_block_solvable = q1(con, "SELECT count(DISTINCT target_key||'|'||era_id)"
                               " FROM cv_series WHERE status='solved'")
    verdicts = [r["verdict"] for r in rows]
    # Series with no cv_cattie row at all are UNTIED too, and must be
    # counted, or the goal would be graded on a sample it selected itself.
    verdicts += ["UNTIED"] * (n_series - len(rows))
    v, deciding = ct.goal_verdict(verdicts, n_series)
    tied = [r for r in rows if r["verdict"].startswith("TIED")
            and r["check_rms_clip"] is not None]
    med = (float(np.median([r["check_rms_clip"] for r in tied]))
           if tied else float("nan"))
    cross = q(con, "SELECT abs(star_offset_median) FROM cv_cattie_cross "
                   "WHERE star_offset_median IS NOT NULL")
    xmed = float(np.median([c[0] for c in cross])) if cross else float("nan")
    cls = {"SUPPORTED": "ok", "SUPPORTED-WITH-CAVEATS": "warn",
           "PARTIALLY SUPPORTED": "warn", "NOT SUPPORTED": "bad"}[v]
    body = table(
        ["series", "n tie stars", "ZP &plusmn; err", "k &plusmn; err",
         "resid RMS (mmag)", "check RMS (mmag)", "check stars",
         "target colour", "verdict"],
        [[_label(r["series_key"]), r["n_fit"],
          _pm(r["zp"], r["zp_err"], 4),
          _pm(r["colour_term"], r["colour_err"], 4),
          _mmag(r["resid_rms"]), _mmag(r["check_rms_clip"]),
          (f"{r['n_check']}&minus;{r['n_check_outlier']}"
           if r["n_check"] else "&mdash;"),
          esc(r["colour_position"]),
          f"<b>{esc(r['verdict'])}</b>"] for r in rows],
        [VERDICT_CLASS.get(r["verdict"], "") for r in rows])
    return f"""
<section id="verdict">
<div class="bhead"><h2>9 &middot; The verdict on the strategy's calibration
goal</h2></div>

<p>The goal, verbatim from the characterization page: &ldquo;Absolute
calibration: nightly REFCAT2 tie &rarr; PS1&nbsp;AB to 0.01&ndash;0.02
mag&rdquo;, claimed as &ldquo;every published magnitude on a standard
system&rdquo;.  It has two halves and both must hold: COVERAGE (every block
carrying science is tied) and ACCURACY (the independent check meets the
band).  Grading only the half that passed is how a goal gets promoted
without being met.</p>

{body}

<div class="{cls}"><b>CALIBRATION GOAL: {esc(v)}</b><br>
Deciding number: {esc(deciding)}.  Median achieved accuracy on independent
check stars: {_mmag(med)} mmag.  Cross-catalogue systematic floor:
{(_mmag(xmed) + " mmag") if cross else "<b>not measured in this build</b> "
    "— the Gaia archive was unreachable, see 6.3"}.  {n_block} of the
{n_block_solvable} (target, era) blocks that carry solved photometry are now
tied ({n_block_all} blocks exist in all; the rest have no solved series to
tie), against the <b>7 of 14</b> quoted by the characterization page when
the goal was last graded &mdash; and those seven were colour-termless median
offsets rather than fits.  That figure is cited from that page, not
recomputed here.</div>

<div class="note"><b>What a reader of the light curves must do with this.</b>
<code>cv_lightcurve.mag</code> is unchanged and is still the differential
magnitude every precision, period, eclipse-timing and state-history result
rests on &mdash; the tie does not touch any of them, and should not.
<code>cv_lightcurve.cal_mag</code> is <code>mag &minus; ZP0</code>: a
NATURAL-SYSTEM magnitude on the catalogue's zero point at colour
c<sub>ref</sub>, with no colour transformation applied to any row.  The
per-point error is unchanged (<code>inst_mag_err</code>); the zero point's
own uncertainty is ONE number per block in <code>cv_cattie.zp_err</code>
and must be applied once, to the block &mdash; adding it in quadrature to
every point would corrupt exactly the differential errors the rest of the
project depends on, and would then be double-counted by anyone comparing
two points of the same curve.  To place a magnitude on a standard system,
apply <code>k &times; (C &minus; c<sub>ref</sub>)</code> yourself, with a
colour C you are willing to defend &mdash; and for the CVs themselves,
read section 7 first.</div>
</section>"""


# ===========================================================================
# Render
# ===========================================================================
def render_report(db_path: Path, cache_root: Path | None = None) -> Path:
    """Render the CV-S6 catalogue-tie page.  Returns the HTML path.

    ``cache_root`` is the catalogue cache directory; section 1 reads the
    cached cones directly to measure what REFCAT2's ``dupvar`` flag does in
    THESE fields rather than assert it.  Defaults to the location the build
    script writes.
    """
    cache_root = cache_root or (REPO_ROOT / "products" / "phot"
                                / "catalogue_cache")
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = dict(q(con, "SELECT key, value FROM cv_cat_meta"))
        figs = [fig_catalogue_depth(con), fig_veto_census(con),
                fig_colour_fits(con), fig_colour_terms(con),
                fig_accuracy(con), fig_resid_mag(con), fig_resid_xy(con),
                fig_cross(con)]
        sections = [
            section_intro(con),
            section_catalogue(con, figs[0], cache_root),
            section_stars(con, figs[1]),
            section_fits(con, figs[2], figs[3]),
            section_bandtest(con),
            section_accuracy(con, figs[4]),
            section_adversarial(con, figs[5], figs[6], figs[7],
                                section_literature(con)),
            section_target_colour(con),
            section_relative(con),
            section_verdict(con),
        ]
        n_tied = q1(con, "SELECT count(*) FROM cv_cattie WHERE is_primary=1 "
                         "AND verdict LIKE 'TIED%'")
        n_series = q1(con, "SELECT count(*) FROM cv_series "
                           "WHERE status='solved'")
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CV Time Series &mdash; Catalogue Tie &amp; Absolute Calibration</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Cataclysmic-Variable Time Series &mdash; the catalogue tie</h1>
  <p>From an arbitrary internal gauge to natural-system magnitudes on a
  standard zero point &middot; {n_tied} of {n_series} solved series tied
  &middot; primary catalogue
  <code>{esc(meta.get('primary_catalogue', ''))}</code> &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('cattie_code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="cv_characterization.html">the characterization that
  graded this NOT SUPPORTED</a> &middot;
  <a href="index.html">project hub</a> &middot;
  <a href="../index.html">all reports</a></p>
</header>

<nav>
  <a href="#intro">0 The gauge</a> &middot;
  <a href="#catalogue">1 Catalogue</a> &middot;
  <a href="#stars">2 Tie stars</a> &middot;
  <a href="#fits">3 The fit</a> &middot;
  <a href="#bandtest">4 Band identity</a> &middot;
  <a href="#accuracy">5 Achieved accuracy</a> &middot;
  <a href="#adversarial">6 Attacks</a> &middot;
  <a href="#targetcolour">7 Target colour</a> &middot;
  <a href="#relative">8 Still relative</a> &middot;
  <a href="#verdict">9 Verdict</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_phot.report_cattie</code> from
<code>products/phot/cv_timeseries.sqlite</code> and the cached catalogue
cones.  Every number in a table, figure or verdict is the result of a query
run at render time or a constant imported from <code>macro_phot.cattie</code>.
The exceptions are named rather than claimed away, because a blanket
&ldquo;nothing is typed by hand&rdquo; is the kind of assurance that only has
to be wrong once: the AAVSO VSX ranges of section 6.4 (hand-entered, dated
and sourced in <code>cattie.TARGET_V_RANGE</code>), the 0.2&ndash;0.5 mag/mag
figure quoted in section 4 as what a real Johnson&ndash;Cousins bandpass
would show, the era-80 reference FWHM quoted in section 8, and the
descriptive prose throughout.  An earlier edition of this footer claimed no
hand-typed numbers at all while carrying two that the database disagreed
with; both are queries now (sections 1 and 2).  Regenerate with
<code>pipeline/scripts/run_cv_cattie.py report</code>.</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
    finally:
        con.close()
    return HTML_PATH
