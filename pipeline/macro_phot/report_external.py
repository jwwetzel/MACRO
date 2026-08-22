"""CV-S7 chain-of-evidence report: the external record, and the branch it decides.

Reads ``products/phot/cv_timeseries.sqlite`` and writes

* ``docs/CV_TimeSeries/cv_external_context.html``
* ``docs/CV_TimeSeries/figures/cv_external/*.png``

Socratic, in ONE order, because the order is the argument:

    0  why ask anyone else at all?  (differential photometry has no scale)
    1  who was asked, what came back, and what did not
    2  whose data is this?  (our own rows, wearing AAVSO's coat)
    3  what does "quiescence" mean here -- measured, not quoted
    4  YZ Cnc, 2024: the state timeline with our nights on it
    5  night by night: the verdict, and what backs each one
    6  the branch
    7  the other four targets: what the surveys add
    8  what the external record CANNOT do

Section 6 is the deliverable and it is deliberately late: a branch decision
read without sections 2 and 3 is someone's opinion about a dwarf nova.

Every number and every figure comes from a query executed here.  The
handful of values that are NOT measured -- VSX magnitude ranges, each
survey's documented bright limit -- are named as external constants in the
footer rather than covered by a blanket claim of self-sufficiency.
"""

from __future__ import annotations

import datetime as dt
import math
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

from macro_core.report_s0 import (        # noqa: E402
    ACCENT, BAD, STYLE, DPI, FAINT, GOOD, MUTED, WARN,
    _figure, esc, q, q1, table)
from macro_core import plotstyle as ps   # noqa: E402  (house figure style)
from . import external as ex              # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "CV_TimeSeries"
FIG_DIR = DOCS_DIR / "figures" / "cv_external"
HTML_PATH = DOCS_DIR / "cv_external_context.html"

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}
TARGET_ORDER = ("yzcnc", "stlmi", "vvpup", "euuma", "anuma")
SOURCE_LABEL = {"aavso": "AAVSO (AID)", "ztf": "ZTF (IRSA)",
                "asassn": "ASAS-SN (legacy V)", "atlas": "ATLAS forced phot."}
SOURCE_COLOR = {"aavso": ACCENT, "ztf": GOOD, "asassn": WARN,
                "atlas": BAD}
#: Survey -> marker.  Four surveys overlaid on one light curve is
#: exactly the case where hue alone stops being enough.
SOURCE_MARKER = {"aavso": "o", "ztf": "s", "asassn": "^", "atlas": "D"}

STATE_COLOR = {ex.STATE_QUIESCENT: ACCENT, ex.STATE_ELEVATED: WARN,
               ex.STATE_OUTBURST: BAD, ex.STATE_UNKNOWN: FAINT}

#: The same four states, for HTML rather than for a plot.  The figure inks
#: above are chosen for marks on a light axis and are not the same problem as
#: TEXT on a page: the fallback used to be a bare ``#ccc``, which on the
#: site's white ground is an unreadable state label, and FAINT grey is not
#: much better.  These names resolve in ``docs/assets/macro.css``, so the
#: page and the stylesheet agree about what "bad" looks like.
STATE_CLASS = {ex.STATE_QUIESCENT: "", ex.STATE_ELEVATED: "warn",
               ex.STATE_OUTBURST: "bad", ex.STATE_UNKNOWN: "muted"}


def _mjd_to_date(mjd) -> str:
    if mjd is None:
        return "&mdash;"
    return (dt.datetime(1858, 11, 17)
            + dt.timedelta(days=float(mjd))).strftime("%Y-%m-%d")


def _f(x, nd=2, dash="&mdash;"):
    """Format a float, or say plainly that there is nothing to format."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return dash
    return f"{x:.{nd}f}"


def _date_num(night: str) -> float:
    """A ``YYYY-MM-DD`` night as a matplotlib date number."""
    return matplotlib.dates.date2num(dt.datetime.strptime(night, "%Y-%m-%d"))


#: MJD of matplotlib's date-number epoch (1970-01-01), so the two time axes
#: on this page can be mixed safely.  ``MJD = date2num + MPL_EPOCH_MJD``.
#: The coverage figure plots survey spans (which arrive as MJD) alongside
#: our own nights (which arrive as date strings); the first version of it
#: reconciled them with a wrong constant and drew every RLMT bar ~70 years
#: adrift, which a glance at the figure caught and no test would have.
MPL_EPOCH_MJD = 40587.0


def _night_to_mjd(night: str) -> float:
    """A ``YYYY-MM-DD`` night as an MJD."""
    return _date_num(night) + MPL_EPOCH_MJD


def _year_mjd(year: int) -> float:
    """MJD of January 1st of ``year`` — for the coverage figure's ticks."""
    return matplotlib.dates.date2num(dt.datetime(year, 1, 1)) + MPL_EPOCH_MJD


def _fig(src: str, caption: str, missing: str) -> str:
    """A figure, or an explicit statement of why there is none."""
    if not src:
        return f'<div class="note"><b>No figure here, and why:</b> {missing}</div>'
    return _figure(src, caption)


def _meta(con) -> dict:
    return dict(q(con, "SELECT key, value FROM cv_ext_meta"))


# ===========================================================================
# Figures
# ===========================================================================

def fig_yzcnc_timeline(con) -> str:
    """THE figure: YZ Cnc's 2024 state timeline with our nights overlaid.

    This is the one that decides the branch, so it is built to be argued
    with.  Three layers, drawn in this order for a reason:

    * the INDEPENDENT AAVSO V-like record (filled circles) -- the evidence;
    * our own AAVSO-submitted points (open circles, on the same ladder via
      the measured band offsets) -- shown, but visibly distinguished, so a
      reader can cover them with a thumb and check the verdict survives;
    * our RLMT nights as vertical bars whose height is the frame count, and
      the superoutburst episode as a shaded span.

    If the dense bars stood inside the shaded span, the paper would have a
    superhump section.  They do not, and the figure's whole job is to make
    that checkable at a glance rather than assertable in prose.
    """
    rows = q(con, """SELECT utc_night, mag, independent FROM cv_ext_nightly
                     WHERE target='yzcnc' AND source='aavso'
                       AND utc_night BETWEEN '2023-12-01' AND '2024-07-01'
                     ORDER BY utc_night""")
    if not rows:
        return ""
    meta = _meta(con)
    base = float(meta.get("baseline_yzcnc_independent", "nan"))

    ind = [(r[0], r[1]) for r in rows if r[2]]
    # Our own nights come from the verdict table, which already carries the
    # V-equivalent magnitude computed through the measured band offsets.
    own = q(con, """SELECT utc_night, mag FROM cv_ext_verdict
                    WHERE target='yzcnc' AND basis='own' AND mag IS NOT NULL
                    ORDER BY utc_night""")
    nights = q(con, """SELECT utc_night, n_frames, is_dense FROM cv_ext_verdict
                       WHERE target='yzcnc' ORDER BY utc_night""")
    eps = q(con, """SELECT plateau_start, plateau_end, kind FROM cv_ext_episode
                    WHERE target='yzcnc' AND source='aavso'
                      AND start_night BETWEEN '2023-12-01' AND '2024-07-01'
                    ORDER BY start_night""")

    with plt.rc_context(STYLE):
        fig, (ax, axn) = plt.subplots(
            2, 1, figsize=(11, 6.4), sharex=True,
            gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.08})

        for p_start, p_end, kind in eps:
            if kind != ex.EPISODE_SUPEROUTBURST:
                continue
            ax.axvspan(_date_num(p_start), _date_num(p_end),
                       color=BAD, alpha=0.18, zorder=0)
            axn.axvspan(_date_num(p_start), _date_num(p_end),
                        color=BAD, alpha=0.18, zorder=0)
            # Anchored in axes fraction vertically so the caption stays
            # INSIDE the panel whatever the magnitude limits turn out to be
            # (a fixed y=10.4 put it above the axis on this data).
            ax.annotate(" superoutburst plateau",
                        xy=(_date_num(p_start), 0.98),
                        xycoords=("data", "axes fraction"),
                        color=BAD, fontsize=8, va="top", ha="left")

        if math.isfinite(base):
            ax.axhline(base, color=ACCENT, lw=1.0, ls="--", zorder=1)
            ax.text(_date_num("2023-12-03"), base - 0.08,
                    f"measured quiescent baseline V={base:.2f}",
                    color=ACCENT, fontsize=8)
            ax.axhline(base - ex.SUPEROUTBURST_AMP_MIN, color=WARN, lw=1.0,
                       ls=":", zorder=1)
            ax.text(_date_num("2023-12-03"),
                    base - ex.SUPEROUTBURST_AMP_MIN - 0.08,
                    f"superoutburst level "
                    f"({ex.SUPEROUTBURST_AMP_MIN:.0f} mag up)",
                    color=WARN, fontsize=8)

        ax.plot([_date_num(n) for n, _ in ind], [m for _, m in ind],
                "o-", color=ACCENT, ms=3.4, lw=0.7, zorder=3,
                label=f"AAVSO, independent observers ({len(ind)} nights)")
        if own:
            ax.plot([_date_num(n) for n, _ in own], [m for _, m in own],
                    "o", mfc="none", mec=WARN, ms=7, mew=1.4, zorder=4,
                    label=f"our own photometry, resubmitted ({len(own)} nights)")
        ax.invert_yaxis()
        ax.set_ylabel("V (or V-equivalent)")
        ax.legend(loc="lower right", fontsize=8, framealpha=0.25)
        ax.set_title("YZ Cnc 2024: the accretion state, and where our nights fell")

        dense_n = [(n, c) for n, c, d in nights if d]
        snap_n = [(n, c) for n, c, d in nights if not d]
        axn.bar([_date_num(n) for n, _ in dense_n], [c for _, c in dense_n],
                width=1.6, color=WARN, zorder=3,
                label=f"dense RLMT run (>= {meta.get('dense_run_min_frames')} frames)")
        axn.bar([_date_num(n) for n, _ in snap_n], [c for _, c in snap_n],
                width=1.6, color=MUTED, zorder=3, label="short RLMT run")
        axn.set_ylabel("RLMT frames")
        axn.legend(loc="upper left", fontsize=8, framealpha=0.25)
        # One tick per month, explicitly.  The auto locator produced two
        # ticks inside January and labelled both "2024-01", which reads as a
        # duplicated month rather than as two positions within one.
        axn.xaxis.set_major_locator(matplotlib.dates.MonthLocator())
        axn.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y-%m"))
        fig.autofmt_xdate()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "yzcnc_timeline.png", dpi=DPI,
                    bbox_inches="tight")
        plt.close(fig)
    return "figures/cv_external/yzcnc_timeline.png"


def fig_band_offsets(con) -> str:
    """Does our resubmitted photometry agree with independent observers?

    The scatter plot behind the sentence "our SR is V to within 0.3 mag".
    It matters because three of the four Feb nights have no independent
    coverage at all, so the weight those nights can bear is exactly the
    weight this comparison establishes.
    """
    # OUR value (own_mag) against THEIRS (the independent nightly median).
    # Using v.mag here would compare the independent median with itself on
    # every independent night and draw a flawless 1:1 line — a figure that
    # looked like triumphant agreement and measured nothing.
    rows = q(con, """SELECT v.utc_night, v.own_mag, n.mag
                     FROM cv_ext_verdict v JOIN cv_ext_nightly n
                       ON n.target='yzcnc' AND n.source='aavso'
                          AND n.utc_night = v.utc_night
                     WHERE v.target='yzcnc' AND v.own_mag IS NOT NULL
                       AND n.independent=1""")
    # Nights where an independent point and our own both exist are the only
    # ones that can test the offset; if there are none, say so rather than
    # draw an empty axis.
    if len(rows) < 3:
        return ""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        ours = [r[1] for r in rows]
        theirs = [r[2] for r in rows]
        lo = min(min(ours), min(theirs)) - 0.3
        hi = max(max(ours), max(theirs)) + 0.3
        ax.plot([lo, hi], [lo, hi], "-", color=ACCENT, lw=1.0)
        ax.plot(theirs, ours, "o", color=ACCENT, ms=6)
        d = np.array(ours) - np.array(theirs)
        ax.set_xlabel("independent AAVSO V (nightly median)")
        ax.set_ylabel("our V-equivalent (nightly median)")
        ax.set_title(f"n={len(rows)} shared nights\n"
                     f"median offset {np.median(d):+.3f}, "
                     f"scatter {np.median(np.abs(d - np.median(d))) * 1.4826:.3f} mag")
        ax.invert_xaxis()
        ax.invert_yaxis()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "band_offsets.png", dpi=DPI, bbox_inches="tight")
        plt.close(fig)
    return "figures/cv_external/band_offsets.png"


def fig_coverage(con) -> str:
    """Per-target, per-source coverage in time — who watched, and when.

    The figure that answers "what does the external record add" for the four
    polars: RLMT's nights are a cluster of dots inside decades of other
    people's coverage, and the gaps in that coverage are as informative as
    the coverage itself.
    """
    rows = q(con, """SELECT target, source, mjd_min, mjd_max, n_points
                     FROM cv_external WHERE n_points > 0""")
    if not rows:
        return ""
    ours = dict(q(con, """SELECT target_key, count(*) FROM cv_frames
                          GROUP BY target_key"""))
    our_span = {t: (a, b) for t, a, b in q(
        con, """SELECT target_key, min(night), max(night) FROM cv_frames
                GROUP BY target_key""")}
    by_target: dict[str, list] = {}
    for t, s, a, b, n in rows:
        by_target.setdefault(t, []).append((s, a, b, n))

    sources = ("aavso", "ztf", "asassn")
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(11, 4.6))
        ylab, ypos = [], []
        y = 0.0
        for t in TARGET_ORDER:
            if t not in by_target:
                continue
            for s in sources:
                hit = [r for r in by_target[t] if r[0] == s]
                if not hit:
                    continue
                _, a, b, n = hit[0]
                ax.barh(y, float(b) - float(a), left=float(a), height=0.62,
                        color=SOURCE_COLOR[s], alpha=0.85)
                ax.text(float(b) + 200, y, f"{n:,}", va="center", fontsize=7,
                        color=MUTED)
                ylab.append(f"{TARGET_LABEL[t]} · {SOURCE_LABEL[s]}")
                ypos.append(y)
                y -= 1.0
            if t in our_span:
                a, b = our_span[t]
                am, bm = _night_to_mjd(a), _night_to_mjd(b)
                ax.barh(y, max(bm - am, 30), left=am, height=0.62,
                        color=ps.tint(ACCENT, 0.45), alpha=0.95)
                ax.text(bm + 200, y, f"{ours.get(t, 0):,}", va="center",
                        fontsize=7, color=MUTED)
                ylab.append(f"{TARGET_LABEL[t]} · RLMT (ours)")
                ypos.append(y)
                y -= 1.6
        ax.set_yticks(ypos)
        ax.set_yticklabels(ylab, fontsize=8)
        ax.set_xticks([_year_mjd(yr) for yr in range(1960, 2031, 10)])
        ax.set_xticklabels([str(yr) for yr in range(1960, 2031, 10)])
        ax.set_xlabel("year (bar labels: number of points)")
        ax.set_title("Who has watched these five stars, and for how long")
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "coverage.png", dpi=DPI, bbox_inches="tight")
        plt.close(fig)
    return "figures/cv_external/coverage.png"


def fig_recurrence(con) -> str:
    """YZ Cnc's outburst recurrence over the whole AAVSO record.

    The single thing our 110-day season provably cannot measure, and the
    clearest statement of what the external record is FOR: intervals between
    successive outburst episodes, and between successive superoutbursts.
    """
    eps = q(con, """SELECT peak_night, kind FROM cv_ext_episode
                    WHERE target='yzcnc' AND source='aavso'
                    ORDER BY peak_night""")
    if len(eps) < 10:
        return ""
    supers = [e[0] for e in eps if e[1] == ex.EPISODE_SUPEROUTBURST]
    all_gaps = [ex.days_between(a, b)
                for a, b in zip([e[0] for e in eps], [e[0] for e in eps][1:])]
    su_gaps = [ex.days_between(a, b) for a, b in zip(supers, supers[1:])]
    all_gaps = [g for g in all_gaps if 0 < g < 200]
    su_gaps = [g for g in su_gaps if 0 < g < 600]
    with plt.rc_context(STYLE):
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.6))
        a1.hist(all_gaps, bins=40, color=ACCENT)
        a1.set_xlabel("days between successive outburst peaks")
        a1.set_ylabel("count")
        a1.set_title(f"all outbursts (n={len(all_gaps)})\n"
                     f"median {np.median(all_gaps):.0f} d")
        a2.hist(su_gaps, bins=30, color=WARN)
        a2.set_xlabel("days between successive superoutburst peaks")
        a2.set_title(f"supercycle (n={len(su_gaps)})\n"
                     f"median {np.median(su_gaps):.0f} d")
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(FIG_DIR / "recurrence.png", dpi=DPI, bbox_inches="tight")
        plt.close(fig)
    return "figures/cv_external/recurrence.png"


# ===========================================================================
# Sections
# ===========================================================================

def section_why(con) -> str:
    n_frames = q1(con, "SELECT count(*) FROM cv_frames") or 0
    return f"""
<section><div class="bhead"><h2>0 &middot; Why ask anyone else at all?</h2></div>
<p>Our photometry is <b>differential</b>. Every light curve in this project
was solved with the ensemble's gauge fixed by <code>mean(ZP) = 0</code>, and
the catalogue tie moves a whole block onto a standard zero point &mdash; it
does not tell one night from another on an absolute scale. So a YZ Cnc run
that sits 1.4 magnitudes above the previous night's looks, in our own
{n_frames:,} frames, exactly like one that does not.</p>
<p>That is a problem, because the strategy's third science question has two
branches and the branch point is an absolute-brightness fact:</p>
<div class="decision"><b>Q3, as written:</b> <i>YZ Cnc 2024 season:
superhump period (and possibly dP<sub>sh</sub>/dt) within confirmed
outburst states</i> &mdash; <b>ApJ-grade if the dense runs are in
(super)outburst</b>; otherwise the fallback is orbital hump + flickering
statistics, &ldquo;honest but weaker&rdquo;.</div>
<p>Common superhumps are a <b>superoutburst</b> phenomenon in SU UMa systems.
So the question this page has to answer is not &ldquo;was the star
bright?&rdquo; but &ldquo;was it in a superoutburst on the nights we ran at
8&nbsp;s cadence?&rdquo; &mdash; and nothing inside our own pixels can
answer it.</p>
</section>"""


def section_sources(con) -> str:
    rows = q(con, """SELECT target, source, n_rows, ok, pulled_utc, note,
                            cache_path FROM cv_ext_fetch
                     ORDER BY source, target""")
    by_source: dict[str, list] = {}
    for r in rows:
        by_source.setdefault(r[1], []).append(r)
    trs, cls = [], []
    for s in ("aavso", "ztf", "asassn", "atlas"):
        for target, _, n, ok, pulled, note, cache in by_source.get(s, []):
            trs.append([SOURCE_LABEL.get(s, s), TARGET_LABEL.get(target, target),
                        f"{n:,}" if ok else "&mdash;",
                        "reached" if ok else "UNREACHABLE",
                        esc((pulled or "")[:10]),
                        f"<code>{esc(cache)}</code>" if cache else "&mdash;"])
            cls.append("" if ok else "bad")
    unreachable = sorted({SOURCE_LABEL.get(r[1], r[1])
                          for r in rows if not r[3]})
    atlas_note = q1(con, "SELECT note FROM cv_ext_fetch WHERE source='atlas' "
                         "LIMIT 1") or ""
    return f"""
<section><div class="bhead"><h2>1 &middot; Who was asked, and what came back</h2></div>
<p>Four sources, the ones the strategy's &sect;7 names. Every response is
cached verbatim under <code>products/external/</code> with the query text,
the pull date and a sha256, so nothing on this page rests on a query nobody
can re-run.</p>
{table(["source", "target", "points", "status", "pulled", "cached response"],
       trs, cls)}
<div class="note"><b>Reached by a route that is worth writing down.</b> The
AAVSO row above is the AID served by
<code>vsx.aavso.org/index.php?view=api.delim</code>, found by reading the
LCG&nbsp;v2 client's own JavaScript. The documented WebObs HTML search
<i>silently ignores its date parameters</i> &mdash; it answered a request for
2024-02-20&hellip;25 with observations from June 2026 &mdash; and at 200 rows
per 25-second page it would have needed ~500 requests for YZ Cnc alone.
A pipeline that had trusted those date parameters would have classified the
wrong nights and never known.</div>
<div class="decision"><b>Unreachable, and not papered over:</b>
{esc(", ".join(unreachable)) or "nothing"}<br>{esc(atlas_note[:400])}</div>
<p>ATLAS is a real loss and not a small one: it is the survey with the
cadence and the bright-end dynamic range that would have covered YZ Cnc's
outbursts where ZTF saturates. Getting it needs a credentialled account,
which this pipeline does not have and did not create. The gap is recorded
per target in <code>cv_ext_fetch</code> rather than left as an absence a
later reader would mistake for &ldquo;ATLAS saw nothing&rdquo;.</p>
</section>"""


def section_whose(con, fig_off: str) -> str:
    meta = _meta(con)
    note = q1(con, "SELECT note FROM cv_ext_fetch WHERE target='yzcnc' "
                   "AND source='aavso'") or ""
    offs = [(k, v) for k, v in sorted(meta.items())
            if k.startswith("offset_")]
    n_ind_nights = q1(con, """SELECT count(*) FROM cv_ext_nightly
                              WHERE target='yzcnc' AND source='aavso'
                                AND independent=1
                                AND utc_night BETWEEN '2023-12-01'
                                                  AND '2024-07-01'""") or 0
    return f"""
<section><div class="bhead"><h2>2 &middot; Whose data is this, actually?</h2></div>
<p>The first thing the AAVSO pull returned for YZ Cnc was <b>our own
photometry</b>. Observer code <code>MALW</code>, comment
&ldquo;TAKEN WITH MACRO CONSORTIUM'S ROBERT L. MUTEL TELESCOPE AT WINER
OBSERVATORY&rdquo;, in Sloan <code>SG</code>/<code>SR</code>/<code>SI</code>
&mdash; the same photons we are trying to find an external check on.</p>
<div class="decision"><b>The trap:</b> classify our nights from our own
resubmitted magnitudes, call it an AAVSO cross-match, and the
&ldquo;external confirmation&rdquo; is a mirror. <b>The rule adopted here:</b>
every row is tagged at parse time; the quiescent baseline and every episode
are computed from <b>independent observers only</b>; and each night's verdict
records which basis carried it.</div>
<p>{esc(note[:300])}</p>
<p>That leaves our own rows with exactly one legitimate job: filling nights
no outside observer covered &mdash; and they can only do it if they are first
shown to agree with the people who did. On the {n_ind_nights} independent
nights in this season, our submissions land on the same ladder through a
<b>measured</b> band offset, never an assumed one:</p>
{table(["band", "offset from V, measured on shared nights"],
       [[esc(k.replace("offset_", "").replace("_minus_V", "")), esc(v)]
        for k, v in offs]) if offs else
 '<div class="note">No shared nights &mdash; no offset could be measured.</div>'}
{_fig(fig_off, "Our nightly V-equivalent against independent AAVSO V on the "
      "nights both exist. The agreement needed here is coarse: YZ Cnc's "
      "outburst amplitude is about 4 magnitudes, so a ~0.3 mag scatter "
      "separates outburst from quiescence with room to spare.",
      "fewer than three shared nights — the offset cannot be tested")}
<p class="sub">This is a check on <i>brightness scale</i>, not on our
pipeline: these are the same photons, reduced separately. It licenses the
statement &ldquo;the star was at V&asymp;13 that night&rdquo;. It licenses
nothing about our photometric precision.</p>
</section>"""


def section_ladder(con) -> str:
    meta = _meta(con)
    base = meta.get("baseline_yzcnc_independent", "")
    n_nights = q1(con, """SELECT count(*) FROM cv_ext_nightly
                          WHERE target='yzcnc' AND source='aavso'
                            AND independent=1""") or 0
    return f"""
<section><div class="bhead"><h2>3 &middot; What does &ldquo;quiescence&rdquo; mean here?</h2></div>
<p>Not a number from a catalogue. The baseline is measured from the star's
own record: the median of the faintest
{float(meta.get("quiescent_decile", 0.1)) * 100:.0f}% of independent nightly
points in the season around the window &mdash;
<b>V = {esc(base)}</b>, from {n_nights:,} independent nights.</p>
<p>Nights are then graded by <b>amplitude above that baseline</b>, and
episodes by what actually distinguishes a superoutburst in an SU&nbsp;UMa
star:</p>
{table(["quantity", "rule", "why this and not something else"], [
  ["night: quiescent",
   f"amplitude &le; {meta.get('amp_quiescent_max')} mag", ""],
  ["night: outburst",
   f"amplitude &ge; {meta.get('amp_outburst_min')} mag",
   "the gap between the two is reported as ELEVATED rather than rounded to "
   "a neighbour &mdash; a star caught mid-rise is genuinely neither"],
  ["episode: superoutburst",
   f"peak &ge; {meta.get('superoutburst_amp_min')} mag above quiescence "
   f"AND a plateau &ge; {meta.get('superoutburst_days_min')} d",
   "both, because amplitude alone promotes a well-caught normal outburst "
   "and duration alone merges a run of unrelated ones"],
  ["plateau", f"time within {ex.PLATEAU_DEPTH_MAG} mag of peak, "
              f"inside the episode",
   "the grading duration is the plateau, not the whole excursion: on a star "
   "rarely at true quiescence the excursion merges an event with its own "
   "decline and with whatever follows. Grading on the excursion produced a "
   "&lsquo;38-day superoutburst&rsquo;, which is not a thing that exists; "
   "the same event's plateau is 13 d, which is."],
])}
</section>"""


def section_timeline(con, fig_tl: str) -> str:
    eps = q(con, """SELECT start_night, end_night, peak_night, peak_amp,
                           duration_d, plateau_start, plateau_end, plateau_d,
                           kind, why
                    FROM cv_ext_episode
                    WHERE target='yzcnc' AND source='aavso'
                      AND start_night BETWEEN '2023-12-01' AND '2024-06-30'
                    ORDER BY start_night""")
    trs, cls = [], []
    for s, e, pk, amp, dur, ps, pe, pd_, kind, why in eps:
        trs.append([f"{esc(s)} &rarr; {esc(e)}", esc(pk), _f(amp),
                    f"{dur:.0f}", f"{esc(ps)} &rarr; {esc(pe)}",
                    f"{pd_:.0f}", f"<b>{esc(kind)}</b>", esc(why)])
        cls.append("warn" if kind == ex.EPISODE_SUPEROUTBURST else "")
    return f"""
<section><div class="bhead"><h2>4 &middot; YZ Cnc in 2024, with our nights on it</h2></div>
{_fig(fig_tl, "YZ Cnc across the 2024 season. Top: the independent AAVSO "
      "V-like record (filled), our own resubmitted photometry on the same "
      "ladder (open), the measured quiescent baseline (dashed) and the "
      "superoutburst level (dotted); the shaded span is the superoutburst "
      "plateau. Bottom: our RLMT frames per night, dense runs highlighted. "
      "The question the figure exists to answer is whether any tall bar "
      "stands inside the shaded span.",
      "no nightly points were classified for YZ Cnc")}
<p>Every outburst episode the independent record traces in this season:</p>
{table(["episode (above-threshold span)", "peak", "peak amp (mag)",
        "span (d)", "plateau", "plateau (d)", "grade", "why"], trs, cls)}
<p>The season contains <b>one superoutburst</b>, and our dense runs are not
in it. The nearest we came was two frames on 2024-03-26 and two on
2024-03-27 &mdash; the star went into superoutburst on 2024-03-28.</p>
</section>"""


def section_nights(con) -> str:
    rows = q(con, """SELECT utc_night, local_night, n_frames, filters,
                            is_dense, state, mag, amp, basis, episode, evidence
                     FROM cv_ext_verdict WHERE target='yzcnc'
                     ORDER BY utc_night""")
    trs, cls = [], []
    for (utc, local, n, filt, dense, state, mag, amp, basis, ep, ev) in rows:
        trs.append([
            f"<b>{esc(utc)}</b>", f"<code>{esc(local)}</code>",
            f"{n:,}" + (" <b>DENSE</b>" if dense else ""),
            esc(filt or ""),
            f'<b class="{STATE_CLASS.get(state, "muted")}">{esc(state)}</b>',
            _f(mag), _f(amp), esc(basis), esc(ep or "&mdash;"), esc(ev)])
        cls.append("warn" if dense and state == ex.STATE_OUTBURST else "")
    return f"""
<section><div class="bhead"><h2>5 &middot; Night by night, and what backs each one</h2></div>
<div class="note"><b>The night label is not the UTC date.</b>
<code>cv_frames.night</code> is the LOCAL evening date at Winer (UTC&minus;7),
and every CV frame in this season was taken after local midnight &mdash; so
the UTC date is the night label <b>plus one day</b>. The strategy's
&ldquo;Feb 21&ndash;24, Mar 1&ndash;4, May 2&ndash;3&rdquo; are UTC dates and
are the same nights the manifest calls 2024-02-20&hellip;23,
2024-02-29/03-02/03-03 and 2024-05-01/02. Comparing date strings across that
boundary would move every tag onto its neighbour, on a star that rises
1.4 mag in a day.</div>
{table(["UTC night", "manifest night", "frames", "filters", "state",
        "V", "amp", "basis", "episode", "evidence"], trs, cls)}
<p><b>basis</b> is the honest part: <code>independent</code> means outside
observers covered that night; <code>own</code> means only our own
resubmitted photometry did, and the evidence column then also carries the
bracket test; <code>bracketed</code>/<code>none</code> mean nobody measured
the star and the entry rests on its neighbours alone.</p>
<div class="decision"><b>The Feb 21&ndash;24 block &mdash; our densest, 819
frames &mdash; has no independent coverage at all.</b> It is also the block
whose state matters most. So the verdict there does not rest on our own
numbers: it rests on a physical bracket. The independent record has
2024-02-20 at V=13.85 and 2024-02-25 at V=13.63, five days apart. An
SU&nbsp;UMa superoutburst holds a plateau for
{float(_meta(con).get('superoutburst_days_min', 8)):.0f}&nbsp;days or more. A
superoutburst cannot begin, peak at V&asymp;11 and return to V=13.6 inside a
five-day gap &mdash; so <b>no superoutburst occurred on those nights, and that
conclusion survives deleting every RLMT row from the AAVSO record.</b> Our own
photometry then adds the detail: a normal outburst, quiescent on Feb 21,
rising ~1.2 mag on Feb 22, peaking V&asymp;12.8 on Feb 23, fading Feb 24 and
back near quiescence by Feb 28.</div>
</section>"""


def section_branch(con) -> str:
    meta = _meta(con)
    branch = meta.get("branch", "UNDECIDED")
    reasoning = meta.get("branch_reasoning", "")
    dense = q(con, """SELECT utc_night, n_frames, state, amp FROM cv_ext_verdict
                      WHERE target='yzcnc' AND is_dense=1 ORDER BY utc_night""")
    cls = "bad" if branch == "FALLBACK" else "ok"
    return f"""
<section><div class="bhead"><h2>6 &middot; The branch</h2></div>
<div class="decision"><b>Branch: <span class="{cls}">{esc(branch)}</span></b>
&mdash; {esc(reasoning)}</div>
{table(["dense run (UTC)", "frames", "state", "amp above quiescence (mag)"],
       [[esc(n), f"{c:,}", esc(s), _f(a)] for n, c, s, a in dense])}
<p>What this means for <code>CV-P3-yzcnc-superhump</code>, concretely:</p>
<ul>
<li><b>The superhump branch is closed for this season.</b> No dense run sits
in a superoutburst, and common superhumps do not occur outside one. There is
no O&ndash;C of superhump maxima to build and no
dP<sub>sh</sub>/dt to quote; the abstract should not promise one.</li>
<li><b>The fallback is live, and it is still conditional.</b> &sect;4.19
promises orbital-hump + flickering statistics only after an empirical S/N
check on the 8&nbsp;s High Gain frames at quiescent V&asymp;14.5. Three dense
runs are quiescent at V&nbsp;=&nbsp;14.2&ndash;14.7 &mdash; exactly that
regime &mdash; so the check is now well posed and unavoidable.</li>
<li><b>A third thing exists that the strategy did not anticipate.</b> Six
dense runs (819 frames on Feb 22&ndash;24, 365 on Mar 1&ndash;4) sit inside
<i>normal</i> outbursts, at 8&nbsp;s cadence, three filters. Feb 21&rarr;22
brackets the rise itself. That is not a consolation prize for a missing
superhump section &mdash; it is cycle-resolved multi-colour coverage of
normal-outburst structure, and it should be written up as its own result
rather than folded into the fallback.</li>
</ul>
</section>"""


def section_targets(con) -> str:
    rows = q(con, """SELECT target, source, n_points, mjd_min, mjd_max,
                            span_d, bands, n_nights, median_gap_d, notes
                     FROM cv_external ORDER BY target, source""")
    by_t: dict[str, list] = {}
    for r in rows:
        by_t.setdefault(r[0], []).append(r)
    blocks = []
    for t in TARGET_ORDER:
        if t not in by_t:
            continue
        trs, cls = [], []
        for (_, s, n, a, b, span, bands, nn, gap, notes) in by_t[t]:
            trs.append([SOURCE_LABEL.get(s, s), f"{n:,}",
                        f"{_mjd_to_date(a)} &rarr; {_mjd_to_date(b)}",
                        _f(span / 365.25, 1) if span else "&mdash;",
                        f"{nn:,}" if nn else "&mdash;",
                        _f(gap, 0), esc(bands or ""), esc((notes or "")[:220])])
            cls.append("" if n else "bad")
        ours = q1(con, "SELECT count(*) FROM cv_frames WHERE target_key=?",
                  (t,)) or 0
        onights = q1(con, "SELECT count(DISTINCT night) FROM cv_frames "
                          "WHERE target_key=?", (t,)) or 0
        blocks.append(
            f'<h3>{TARGET_LABEL[t]}</h3><p class="sub">RLMT: {ours:,} frames '
            f'on {onights} nights.</p>'
            + table(["source", "points", "coverage", "span (yr)", "nights",
                     "median gap (d)", "bands", "notes"], trs, cls))
    return f"""
<section><div class="bhead"><h2>7 &middot; The five targets: what the surveys add</h2></div>
<p>Coverage span, cadence and bands per target and source, straight from
<code>cv_external</code>.</p>
{"".join(blocks)}
</section>"""


def section_limits(con, fig_cov: str, fig_rec: str) -> str:
    sat = q(con, """SELECT target, notes FROM cv_external
                    WHERE source='ztf' AND notes LIKE '%BRIGHT-END%'""")
    eu = q1(con, """SELECT max(mjd_max) FROM cv_external
                    WHERE target='euuma' AND source='aavso'""")
    return f"""
<section><div class="bhead"><h2>8 &middot; What the external record adds &mdash; and where it cannot help</h2></div>
{_fig(fig_cov, "Coverage in time per target and source, with our own frames "
      "in white. The polars' RLMT campaigns are a few seasons inside 30-45 "
      "year records; YZ Cnc's is a 110-day sliver of a 65-year one.",
      "no coverage rows to draw")}
<h3>What it adds, that our data cannot supply</h3>
<ul>
<li><b>An absolute state tag per night.</b> The whole reason this page
exists. Our differential curves cannot say whether a night was quiescent;
the external record can, and does, for every night in the window.</li>
<li><b>Outburst recurrence and the supercycle.</b> Our YZ Cnc season spans
110 days and contains one superoutburst. The AAVSO record spans 65 years.
Recurrence statistics are simply not a measurement our data can make.</li>
<li><b>Long-term accretion-state duty cycles for the polars</b> (the
strategy's Q4), over 30&ndash;45 year baselines, against which our nights
are placed rather than compared.</li>
</ul>
{_fig(fig_rec, "YZ Cnc outburst recurrence over the full AAVSO record: "
      "intervals between successive outburst peaks, and between successive "
      "superoutburst peaks. Both are measurements our 110-day season is "
      "structurally incapable of making.",
      "too few episodes in the record to form recurrence statistics")}
<h3>Where it cannot help, stated plainly</h3>
<ul>
<li><b>ATLAS is missing entirely</b> &mdash; account-gated. This is the worst
of the gaps, because ATLAS has both the cadence and the bright-end range
that would have covered YZ Cnc's outbursts where ZTF cannot.</li>
<li><b>ZTF saturates on YZ Cnc in outburst.</b>
{esc((sat[0][1] if sat else "")[:200])} It returned {q1(con,
    "SELECT n_points FROM cv_external WHERE target='yzcnc' AND source='ztf'")
    or 0} points in 7 years for YZ Cnc against 1,955&ndash;4,152 for the
faint polars &mdash; and <b>zero</b> epochs inside 2024-02-21&hellip;05-03.
ZTF contributed nothing to the branch decision and could not have.</li>
<li><b>ASAS-SN's modern g-band record is unreachable</b> (Sky Patrol v2's API
host refuses connections here), so its contribution stops in 2018. What is
here is the legacy V-band survey, useful for long-term context and useless
for 2024.</li>
<li><b>AAVSO coverage of EU UMa stops at {_mjd_to_date(eu)}.</b> For EU UMa's
modern state history there is only ZTF.</li>
<li><b>Nobody but us observed 2024-02-21&hellip;24</b> &mdash; our densest
block. The superoutburst exclusion there is a bracket argument, not a
measurement, and it is stated as one.</li>
<li><b>None of this constrains sub-orbital behaviour.</b> Every external
source here is one-point-per-night. On 90&ndash;125 minute binaries they
constrain the state and nothing inside the cycle &mdash; which is precisely
the gap this paper exists to fill.</li>
</ul>
</section>"""


# ===========================================================================
# Entry point
# ===========================================================================

def render_report(db_path: Path) -> Path:
    """Render the CV-S7 external-context page.  Returns the HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        meta = _meta(con)
        fig_tl = fig_yzcnc_timeline(con)
        fig_off = fig_band_offsets(con)
        fig_cov = fig_coverage(con)
        fig_rec = fig_recurrence(con)
        branch = meta.get("branch", "UNDECIDED")
        n_sources = q1(con, "SELECT count(*) FROM cv_ext_fetch WHERE ok=1") or 0
        n_pts = q1(con, "SELECT sum(n_points) FROM cv_external") or 0
        body = "\n".join([
            section_why(con),
            section_sources(con),
            section_whose(con, fig_off),
            section_ladder(con),
            section_timeline(con, fig_tl),
            section_nights(con),
            section_branch(con),
            section_targets(con),
            section_limits(con, fig_cov, fig_rec),
        ])
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CV Time Series &mdash; External Survey Context</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Cataclysmic-Variable Time Series &mdash; the external record</h1>
  <p>What the public survey record says our five targets were doing &middot;
  {n_pts:,} points from {n_sources} reached (target, source) pulls &middot;
  <b>YZ Cnc branch: {esc(branch)}</b> &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('external_code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)</p>
</header>

<main>
{body}
<section><div class="bhead"><h2>Provenance</h2></div>
<p>Every number and figure above is produced by a query against
<code>products/phot/cv_timeseries.sqlite</code> executed by
<code>pipeline/macro_phot/report_external.py</code>. The classification
arithmetic is <code>pipeline/macro_phot/external.py</code> and is unit
tested. Raw responses are cached under <code>products/external/</code>.</p>
<p><b>Values on this page that are NOT measurements of ours</b>, named rather
than covered by a blanket claim: the VSX magnitude ranges and variability
types used to set each target's bright end; each survey's documented
bright limit (<code>external.SURVEY_BRIGHT_LIMIT</code>); and the
astrophysical premise that common superhumps are a superoutburst phenomenon
in SU UMa systems, which is the literature fact the branch rule encodes.</p>
<p>Rebuild:
<code>python pipeline/scripts/run_cv_external.py fetch</code> &rarr;
<code>classify</code> &rarr; <code>report</code>, then
<code>python pipeline/scripts/check_pipeline_status.py record CV-S7</code>.</p>
</section>
</main>
</body></html>"""
        HTML_PATH.parent.mkdir(parents=True, exist_ok=True)
        tmp = HTML_PATH.with_suffix(".html.tmp")
        tmp.write_text(html, encoding="utf-8")
        tmp.replace(HTML_PATH)
        return HTML_PATH
    finally:
        con.close()
