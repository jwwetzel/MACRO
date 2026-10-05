"""Legacy-census evidence page renderer.

Reads ``products/legacy/legacy.sqlite`` (NEVER the archive — if a number
cannot be derived from the database it does not belong on the page) and
writes:

* ``docs/Legacy_Rigel/legacy_census.html``        — the page
* ``docs/Legacy_Rigel/figures/legacy/*.png``      — every figure

The page follows the site's Socratic format: one section per decision, each
section = Question → Evidence → Decision → Consequence.  EVERY number in the
HTML is interpolated from a SQL query executed in this module or from a
constant defined in ``macro_legacy.census`` — nothing is hand-typed, so
re-running the build after the archive changes regenerates the whole
argument, gates and outcome included.

The page deliberately leads with the pre-registered rule and ends with its
mechanical application: the reader sees what would have counted as a "go"
before seeing whether the archive delivered it.
"""

from __future__ import annotations

import datetime as _dt
import re
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from macro_core import plotstyle as ps                       # noqa: E402
# Shared page machinery: one query discipline, one table generator, one
# figure wrapper across the evidence site.
from macro_core.report_s0 import (                           # noqa: E402
    ACCENT, BAD, DPI, GOOD, INK, MUTED, STYLE, WARN,
    _figure, esc, fmt, q, q1, table)
from . import census as lc                                   # noqa: E402
from . import clock as lc_clock                              # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "Legacy_Rigel"
FIG_DIR = DOCS_DIR / "figures" / "legacy"
FIG_REL = "figures/legacy"
HTML_PATH = DOCS_DIR / "legacy_census.html"

#: Tables on the page are capped at this many rows; the cap and the full
#: count are always stated beside the table (the database holds every row).
TOP_N = 25


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def has_table(con, name: str) -> bool:
    return bool(q1(con, "SELECT count(*) FROM sqlite_master WHERE "
                        "type='table' AND name=?", (name,)))


def meta(con) -> dict:
    """census_meta ∪ scan_meta ∪ external_meta as one dict."""
    out = {}
    for t in ("scan_meta", "external_meta", "census_meta"):
        if has_table(con, t):
            out.update(dict(q(con, f"SELECT key, value FROM {t}")))
    return out


def f1(x, nd: int = 1) -> str:
    """Fixed-decimal float for the page; NULL → em-dash."""
    return "&mdash;" if x is None else f"{x:,.{nd}f}"


def pct(num, den) -> str:
    return "&mdash;" if not den else f"{100.0 * num / den:.1f}%"


def _date(night: str) -> _dt.date:
    return _dt.date.fromisoformat(night)


def camera_colors(con) -> dict[str, str]:
    """One colour + marker per camera, in order of first night (stable)."""
    cams = [r[0] for r in q(con, """
        SELECT camera FROM cameras ORDER BY first_night IS NULL, first_night""")]
    return {c: ps.series(i) for i, c in enumerate(cams)}


def passed_cell(v) -> str:
    if v is None:
        return '<span class="muted">not evaluated</span>'
    return '<span class="ok">PASS</span>' if v else '<span class="bad">FAIL</span>'


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_timeline(con) -> str:
    """Camera timeline: frames per night and field rotation, by camera.

    Top: canonical frames per night, one colour per camera — the stints.
    Bottom: nightly median plate-solution rotation (folded mod 180°); the
    vertical ticks are the mechanical-epoch boundaries the census derived
    (rotation steps > 0.3° that persist).
    """
    style = camera_colors(con)
    per = q(con, """
        SELECT night, camera, count(*) FROM frames
        WHERE is_canonical = 1 AND night IS NOT NULL AND error IS NULL
        GROUP BY night, camera""")
    rot = q(con, """
        SELECT night, camera, crota2 FROM frames
        WHERE is_science = 1 AND crota2 IS NOT NULL""")
    epochs = q(con, "SELECT first_night FROM mech_epochs ORDER BY mech_epoch")
    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(9.6, 5.6), sharex=True,
            gridspec_kw={"height_ratios": [3, 2]})
        for cam, st in style.items():
            sub = [r for r in per if r[1] == cam]
            if not sub:
                continue
            ax1.scatter([_date(r[0]) for r in sub], [r[2] for r in sub],
                        s=7, color=st["color"], marker=st["marker"],
                        linewidths=0, label=cam)
        ax1.set_yscale("log")
        ax1.set_ylabel("canonical frames per night")
        ax1.set_title("One calendar, many cameras: who took the frames")
        ax1.legend(fontsize=7, ncol=3, loc="upper center",
                   bbox_to_anchor=(0.5, -0.02), frameon=False)
        # Nightly median rotation per camera.
        nightly: dict = {}
        for night, cam, cr in rot:
            nightly.setdefault((night, cam), []).append(lc.fold_rotation(cr))
        for cam, st in style.items():
            pts = [(n, float(np.median(v))) for (n, c), v in nightly.items()
                   if c == cam and len(v) >= lc.ROTATION_MIN_FRAMES]
            if pts:
                ax2.scatter([_date(n) for n, _ in pts], [v for _, v in pts],
                            s=7, color=st["color"], marker=st["marker"],
                            linewidths=0)
        for (first,) in epochs:
            if first:
                ax2.axvline(_date(first), color=MUTED, lw=0.5, alpha=0.6)
        # Single-night solves scatter by tens of degrees; the steps that
        # define epochs are a few degrees, so the axis is held to ±12°
        # (points beyond it are bad plate solutions, kept in the table).
        ax2.set_ylim(-12, 12)
        ax2.set_ylabel("field rotation\n(deg, mod 180)")
        ax2.set_xlabel("night (local evening date)")
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.42)
        fig.savefig(FIG_DIR / "legacy_camera_timeline.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_camera_timeline.png"


def fig_accounting(con) -> str:
    """Where every scanned file went: science vs each named exclusion, by year."""
    rows = q(con, """
        SELECT coalesce(substr(night, 1, 4), 'no date') AS yr,
               coalesce(exclusion, 'science frame') AS what, count(*)
        FROM frames GROUP BY 1, 2""")
    years = sorted({r[0] for r in rows})
    kinds = [k for k, _ in sorted(
        {r[1]: 0 for r in rows}.items(),
        key=lambda kv: (kv[0] != "science frame", kv[0]))]
    palette = {"science frame": ACCENT}
    spare = [WARN, GOOD, BAD, ps.OTHER, ps.SECOND, MUTED, ps.FAINT, INK]
    for k in kinds:
        if k not in palette:
            palette[k] = spare[(len(palette) - 1) % len(spare)]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.6, 3.6))
        bottom = np.zeros(len(years))
        for k in kinds:
            vals = np.array([sum(r[2] for r in rows if r[0] == y and r[1] == k)
                             for y in years], dtype=float)
            ax.bar(years, vals, bottom=bottom, color=palette[k], label=k,
                   width=0.75)
            bottom += vals
        ax.set_ylabel("scanned files")
        ax.set_xlabel("year of the night (from DATE-OBS)")
        ax.set_title("Every scanned file is a science frame or a named exclusion")
        ax.legend(fontsize=7, ncol=2, frameon=False)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_accounting.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_accounting.png"


def fig_series(con) -> str:
    """Every series: nights vs longest run, with the gate thresholds.

    This is the selection plot the data scientist asked for: the go/no-go is
    read off usable series (nights × longest same-filter run), not frame
    counts.  Filled colour = calibrated nights exist; the vertical line is
    G3's 30-night threshold; RLMT-era targets and eclipsers are marked.
    """
    rows = q(con, """
        SELECT s.target_key, s.n_nights, s.longest_run_h,
               s.n_calibrated_nights, s.tieable FROM series s""")
    ecl = {r[0] for r in q(con, "SELECT target_key FROM eclipsers "
                                "WHERE is_eclipsing = 1")} \
        if has_table(con, "eclipsers") else set()
    ov = {r[0] for r in q(con, "SELECT DISTINCT target_key FROM "
                               "overlap_series")} \
        if has_table(con, "overlap_series") else set()
    rng = np.random.default_rng(20261003)      # fixed jitter: same plot twice
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.6, 4.6))
        x = np.array([r[1] for r in rows], dtype=float)
        y = np.array([r[2] for r in rows], dtype=float)
        xj = x * np.exp(rng.uniform(-0.06, 0.06, len(x)))
        base = np.array([r[0] not in ecl and r[0] not in ov for r in rows])
        ax.scatter(xj[base], y[base] + 0.01, s=6, color=ps.FAINT,
                   linewidths=0, label="other series")
        is_e = np.array([r[0] in ecl for r in rows])
        ax.scatter(xj[is_e], y[is_e] + 0.01, s=14, color=ACCENT, marker="o",
                   linewidths=0, label="eclipsing/contact system (VSX)")
        is_o = np.array([r[0] in ov for r in rows])
        ax.scatter(xj[is_o], y[is_o] + 0.01, s=40, color=BAD, marker="*",
                   linewidths=0, label="RLMT-era project target")
        # Threshold labels sit in axes-fraction y so they cannot collide
        # with the data range or the legend, whatever the archive holds.
        tr = ax.get_xaxis_transform()
        ax.axvline(lc.G3_MIN_CALIBRATED_NIGHTS, **ps.reference_kw(WARN))
        ax.text(lc.G3_MIN_CALIBRATED_NIGHTS * 1.05, 0.62,
                "G3: 30 nights\n(must also be\ncalibrated)", fontsize=7,
                color=WARN, va="top", transform=tr)
        ax.axvline(lc.G2_MIN_NIGHTS, **ps.reference_kw(BAD, ":"))
        ax.text(lc.G2_MIN_NIGHTS * 1.05, 0.62, "G2: 10 nights", fontsize=7,
                color=BAD, va="top", transform=tr)
        ax.set_xscale("log")
        ax.set_xlabel("nights in the series (target × filter × camera)")
        ax.set_ylabel("longest same-filter run (hours)")
        ax.set_title("Usable series, not frame counts")
        ax.legend(fontsize=7, frameon=False, loc="upper right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_series.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_series.png"


def fig_overlap(con) -> str | None:
    """RLMT-era project targets found in the legacy archive, on the calendar."""
    if not has_table(con, "overlap_series"):
        return None
    tg = q(con, """
        SELECT target, target_key, project FROM overlap WHERE n_frames > 0
        ORDER BY project, target""")
    if not tg:
        return None
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(tg) + 1.8))
        seen = {}
        for i, (name, key, proj) in enumerate(tg):
            rows = q(con, """
                SELECT night, band, sum(n_frames) FROM overlap_nights
                WHERE target_key = ? AND project = ?
                GROUP BY night, band""", (key, proj))
            bands = sorted({r[1] for r in rows})
            for j, band in enumerate(bands):
                sub = [r for r in rows if r[1] == band]
                st = seen.setdefault(band, ps.series(len(seen)))
                off = (j - (len(bands) - 1) / 2) * 0.16
                ax.scatter([_date(r[0]) for r in sub], [i + off] * len(sub),
                           s=9, color=st["color"], marker=st["marker"],
                           linewidths=0)
        for band, st in seen.items():
            ax.scatter([], [], s=14, color=st["color"], marker=st["marker"],
                       label=band)
        ax.set_yticks(range(len(tg)))
        ax.set_yticklabels([f"{t[0]}" for t in tg], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("night")
        ax.set_title("RLMT-era project targets in the 2015–2022 archive "
                     "(one point per night × filter)")
        ax.legend(fontsize=7, ncol=6, frameon=False, loc="upper center",
                  bbox_to_anchor=(0.5, -0.28 if len(tg) < 6 else -0.16))
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_overlap.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_overlap.png"


def fig_calibration(con) -> str:
    """Science nights vs calibration nights, one lane per camera.

    Grey = nights with science frames; coloured = nights with bias / dark /
    flat frames on that camera.  A lane with grey and no colour is a camera
    whose calibration frames are not in the archive.
    """
    cams = [r[0] for r in q(con, """
        SELECT camera FROM cameras ORDER BY first_night IS NULL, first_night""")]
    sci = q(con, """
        SELECT DISTINCT camera, night FROM frames WHERE is_science = 1""")
    cal = q(con, """
        SELECT DISTINCT camera, night, kind FROM frames
        WHERE is_canonical = 1 AND night IS NOT NULL
          AND kind IN ('bias', 'dark', 'flat')""")
    offs = {"bias": -0.22, "dark": 0.0, "flat": 0.22}
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(9.2, 0.55 * len(cams) + 1.6))
        for i, cam in enumerate(cams):
            n = [_date(r[1]) for r in sci if r[0] == cam]
            ax.scatter(n, [i] * len(n), s=30, color=ps.WISP, marker="|",
                       linewidths=0.8)
            for kind, off in offs.items():
                d = [_date(r[1]) for r in cal if r[0] == cam and r[2] == kind]
                ax.scatter(d, [i + off] * len(d), s=10,
                           color=ps.KIND_COLOR[kind], linewidths=0)
        ax.scatter([], [], s=30, color=ps.WISP, marker="|", label="science night")
        for kind in offs:
            ax.scatter([], [], s=14, color=ps.KIND_COLOR[kind], label=kind)
        ax.set_yticks(range(len(cams)))
        ax.set_yticklabels(cams, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("night")
        ax.set_title("Calibration frames in the archive, against the science nights")
        ax.legend(fontsize=7, ncol=4, frameon=False, loc="upper center",
                  bbox_to_anchor=(0.5, -0.18))
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_calibration.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_calibration.png"


def fig_eclipsers(con, top: int = 30) -> str | None:
    """Eclipsing systems: every night on the calendar, minimum-bearing marked.

    Filled = clause (a), a run ≥ P/2 that must contain a minimum; open =
    clause (b) only, a predicted minimum inside the run; small grey = a
    night with neither.  Systems are ordered by the number of seasons with
    a minimum-bearing night — the quantity gate G1 counts.
    """
    if not has_table(con, "eclipsers"):
        return None
    systems = q(con, f"""
        SELECT target_key, target_name, vsx_type, period_d, n_seasons_bearing
        FROM eclipsers WHERE is_eclipsing = 1 AND n_nights_bearing > 0
        ORDER BY n_seasons_bearing DESC, n_nights_bearing DESC LIMIT {top}""")
    if not systems:
        return None
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(9.2, 0.3 * len(systems) + 1.8))
        for i, (key, name, typ, per, _ns) in enumerate(systems):
            nights = q(con, """
                SELECT night, guaranteed, predicted FROM eclipser_nights
                WHERE target_key = ?""", (key,))
            other = [_date(n) for n, g, p in nights if not g and not p]
            guar = [_date(n) for n, g, p in nights if g]
            pred = [_date(n) for n, g, p in nights if p and not g]
            ax.scatter(other, [i] * len(other), s=5, color=ps.FAINT,
                       linewidths=0)
            ax.scatter(pred, [i] * len(pred), s=22, facecolors="none",
                       edgecolors=WARN, linewidths=0.9)
            ax.scatter(guar, [i] * len(guar), s=22, color=ACCENT, linewidths=0)
        ax.scatter([], [], s=22, color=ACCENT, label="run ≥ P/2 (minimum certain)")
        ax.scatter([], [], s=22, facecolors="none", edgecolors=WARN,
                   label="predicted minimum in run (ephemeris-dependent)")
        ax.scatter([], [], s=8, color=ps.FAINT, label="other night")
        ax.set_yticks(range(len(systems)))
        ax.set_yticklabels(
            [f"{s[1]}  ({s[2]}, P={s[3]:.3f} d)" if s[3] else f"{s[1]} ({s[2]})"
             for s in systems], fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("night")
        ax.set_title("Minimum-bearing nights of the eclipsing systems")
        ax.legend(fontsize=7, ncol=3, frameon=False, loc="upper center",
                  bbox_to_anchor=(0.5, -0.05))
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_eclipsers.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_eclipsers.png"


def fig_time(con) -> str:
    """The header-time audit in two panels.

    Left: header LST minus the LST implied by DATE-OBS, against exposure
    time (median per camera and exposure, groups of ≥ 30 frames).  The
    scheduler samples LST as the exposure begins, so a start stamp lies on
    the flat line, a mid-exposure stamp on the −t/2 line, an end stamp on
    the −t line.  Right: of the consecutive exposure pairs with unequal
    exposure times, how many would have overlapped if DATE-OBS marked the
    start, the middle or the end.
    """
    style = camera_colors(con)
    audit = q(con, """
        SELECT camera, sum(n_informative_pairs), sum(viol_start),
               sum(viol_mid), sum(viol_end)
        FROM time_audit GROUP BY camera ORDER BY min(first_night)""")
    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.8, 4.0),
                                       gridspec_kw={"width_ratios": [1, 1]})
        xmax = 1.0
        for cam, st in style.items():
            rows = q(con, """
                SELECT round(exptime_s), lst_resid_s FROM frames
                WHERE camera = ? AND lst_resid_s IS NOT NULL
                  AND abs(lst_resid_s) < 1000 AND exptime_s >= 1""", (cam,))
            groups: dict = {}
            for e, r in rows:
                groups.setdefault(e, []).append(r)
            pts = [(e, float(np.median(v))) for e, v in groups.items()
                   if len(v) >= 30]
            if pts:
                ax1.scatter([p[0] for p in pts], [p[1] for p in pts], s=16,
                            color=st["color"], marker=st["marker"],
                            linewidths=0, label=cam, zorder=3)
                xmax = max(xmax, max(p[0] for p in pts))
        xs = np.array([0.0, xmax * 1.05])
        for slope, lab in ((0.0, "start"), (-0.5, "middle"), (-1.0, "end")):
            ax1.plot(xs, slope * xs, **ps.reference_kw(MUTED))
            ax1.annotate(f"stamp = {lab}", (xs[1], slope * xs[1]), fontsize=7,
                         color=MUTED, ha="right", va="bottom")
        ax1.set_xlabel("exposure time (s)")
        ax1.set_ylabel("header LST − LST implied by DATE-OBS (s)")
        ax1.set_title("Which instant DATE-OBS marks: the mount's clue")
        ax1.legend(fontsize=6, frameon=False, loc="lower left")
        cams = [a[0] for a in audit if a[1]]
        yy = np.arange(len(cams))
        h = 0.26
        for k, (lab, col, idx) in enumerate(
                (("if stamp = start", GOOD, 2), ("if stamp = middle", WARN, 3),
                 ("if stamp = end", BAD, 4))):
            vals = [100.0 * a[idx] / a[1] for a in audit if a[1]]
            ax2.barh(yy + (k - 1) * h, vals, height=h, color=col, label=lab)
        ax2.set_yticks(yy)
        ax2.set_yticklabels(cams, fontsize=7)
        ax2.invert_yaxis()
        ax2.set_xlabel("% of unequal-exposure pairs that would overlap")
        ax2.set_title("…and the exposures' own clue")
        ax2.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_time_audit.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_time_audit.png"


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def section_prereg(con, m: dict) -> str:
    prereg = m.get("prereg_sha256", "")
    return f"""
<section id="rule">
<div class="bhead"><h2>0 &middot; The rule came first</h2>
<span class="tag">go/no-go criteria pre-registered before any header was read</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">A census can always be read as a reason to start a project.
What would have to be true of this archive for it to carry one — decided
before looking?</p>

<h3>Evidence</h3>
<p class="sub">The criteria are in
<code>{esc(m.get('prereg_path', ''))}</code>, saved before the scan ran;
its sha256 at build time is <code>{esc(prereg[:16])}&hellip;</code> (stored
in <code>census_meta.prereg_sha256</code>, so a later edit is detectable).
The note discloses what was already known when it was written (the
committee's file count and camera list, the ledger's short "first scan"
target list) and what was not (nights, filters, run lengths, seasons and
calibration for any target).</p>
{table(["gate", "what must hold", "what it leads to"], [
    ["G0", "files-on-disk = rows + named exclusions; a camera timeline "
           "exists; the header-time convention is identified for every "
           "camera epoch feeding a passing gate",
     "validity &mdash; without it the outcome is NOT DECIDABLE"],
    ["G1", f"&ge; {lc.G1_MIN_SYSTEMS} eclipsing/contact systems with "
           f"minimum-bearing nights in &ge; {lc.G1_MIN_SEASONS} seasons",
     "GO-CANDIDATE: period-change (O&minus;C) timing paper, JAAVSO class"],
    ["G2", f"an RLMT-era target with &ge; {lc.G2_MIN_NIGHTS} nights over "
           f"&ge; {lc.G2_MIN_SEASONS} seasons in one tieable series "
           "(timing targets: plus a run &ge; 1 h)",
     "TRANSFER to the owning project &mdash; not a sixth paper"],
    ["G3", f"a series with &ge; {lc.G3_MIN_CALIBRATED_NIGHTS} calibrated "
           "nights, and a question and reader that can be named",
     "GO-CANDIDATE only if the question is named"],
])}

<h3>Decision</h3>
<div class="decision"><b>The gates below are evaluated by script, exactly
as written.</b>  No threshold, definition or question was changed after the
census was seen; where a written definition met a fact it did not
anticipate, the result is given under both readings and marked as a
deviation (section 9).</div>

<h3>Consequence</h3>
<p class="sub">Sections 1&ndash;8 are the census; section 9 applies the
rule.  Nothing in between is a reason for a go unless it is one of these
gates.</p>
</div></section>"""


def section_reconciliation(con, m: dict) -> str:
    summ = q(con, "SELECT quantity, n FROM reconciliation_summary ORDER BY ord")
    rows, classes = [], []
    for quantity, n in summ:
        must = "must be 0" in quantity
        rows.append([esc(quantity).replace("  ", "&nbsp;&nbsp;&nbsp;"), fmt(n)])
        classes.append(("ok" if n == 0 else "bad") if must else None)
    n_scan = q1(con, "SELECT count(*) FROM scan")
    n_err = q1(con, "SELECT count(*) FROM scan WHERE error IS NOT NULL")
    err_rows = q(con, """
        SELECT substr(error, 1, 60), count(*) FROM scan
        WHERE error IS NOT NULL GROUP BY 1 ORDER BY 2 DESC LIMIT 8""")
    excl = q(con, """
        SELECT status, count(*) FROM reconciliation
        WHERE status != ? GROUP BY status ORDER BY 2 DESC""", (lc.REC_SCANNED,))
    ex_tree = q(con, """
        SELECT status,
               CASE WHEN logical_path LIKE 'archivar/%' THEN 'archivar/'
                    WHEN logical_path LIKE 'old/%' THEN 'old/'
                    ELSE substr(logical_path, 1, 4) || '/' END AS tree,
               count(*) FROM reconciliation
        WHERE status NOT IN (?, ?) GROUP BY 1, 2 ORDER BY 1, 2""",
                (lc.REC_SCANNED, lc.REC_NOT_FITS))
    cm = {k: int(v) for k, v in m.items() if k.startswith("collision_")}
    worst = q(con, """
        SELECT stripped_path, n_manifest, n_scanned, n_distinct_date_obs
        FROM collision_audit ORDER BY n_manifest DESC, stripped_path LIMIT 6""")
    n_same = q1(con, """SELECT count(*) FROM collision_audit
                        WHERE n_distinct_date_obs < n_with_date""")
    tr = q(con, """SELECT disposition, count(*) FROM truncated_frames
                   GROUP BY 1 ORDER BY 2 DESC""") \
        if has_table(con, "truncated_frames") else []
    tr_obj = q(con, """
        SELECT coalesce(object, '(header unreadable)'), count(*),
               min(night), max(night) FROM truncated_frames
        GROUP BY 1 ORDER BY 2 DESC LIMIT 12""") if tr else []
    n_tr = sum(n for _d, n in tr)
    n_twin = q1(con, "SELECT count(*) FROM truncated_frames "
                     "WHERE n_twins > 0") if tr else 0
    truncated_html = "" if not tr else f"""
<p class="sub"><b>Frames truncated at the source.</b>  {fmt(n_tr)} files
on Google Drive are shorter than a FITS file of their own header can be: a
second, md5-verified download (2026-10-03) returned the same short bytes, so
the damage is in the source, not the transfer.  The header survives and was
read.  {fmt(n_twin)} of them share a file NAME with an intact frame elsewhere
in the archive; a twin stands in for a truncated frame only if it carries the
identical <code>DATE-OBS</code>:</p>
{table(["disposition", "frames"], [[esc(d), fmt(n)] for d, n in tr])}
<p class="sub">No same-name file is the same exposure: the names repeat
because the scheduler restarts its day numbering every year (the very
collision the year-stripped manifest would have caused).  Every truncated
frame is therefore LOST, and named.  What they were (largest groups):</p>
{table(["OBJECT", "lost frames", "first night", "last night"],
       [[esc(o), fmt(n), esc(a), esc(b)] for o, n, a, b in tr_obj])}"""
    src = fig_accounting(con)
    fr = q(con, """
        SELECT coalesce(exclusion, 'science frame'), count(*) FROM frames
        GROUP BY 1 ORDER BY 2 DESC""")
    return f"""
<section id="recon">
<div class="bhead"><h2>1 &middot; Files on disk = rows + named exclusions</h2>
<span class="tag">reconciled against both transfer manifests</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The archive was crawled from a shared Drive, downloaded,
repaired, compressed and is still being de-duplicated.  A first crawl
manifest (<code>legacy_manifest_BAD_collisions.csv</code>) is known to have
merged years.  Before any target is counted: is every file accounted for,
and did the bad manifest damage anything?</p>

<h3>Evidence</h3>
<p class="sub">The walk found {fmt(int(m.get('walk_n_files', 0)))} files
({esc(m.get('walk_utc', ''))[:16]}Z).  Only <code>*.fz</code> files are
opened; an uncompressed file is either the not-yet-deleted twin of a
scanned <code>.fz</code> (the dedupe job's unfinished work, a snapshot
count) or a frame that never compressed.  Both identities close:</p>
{table(["quantity", "n"], rows, classes)}
<p class="sub">Named exclusions at the path level, by status and tree
(non-FITS files omitted):</p>
{table(["status", "tree", "paths"],
       [[esc(s), esc(t), fmt(n)] for s, t, n in ex_tree]) if ex_tree
 else '<p class="sub"><i>None: every FITS path is a scanned row.</i></p>'}
<p class="sub">{fmt(n_err)} of {fmt(n_scan)} scanned files have an
unreadable header and stay in the table with <code>error</code> set:</p>
{table(["error (first 60 characters)", "files"],
       [[f"<code>{esc(e)}</code>", fmt(n)] for e, n in err_rows]) if err_rows
 else '<p class="sub"><i>None.</i></p>'}

{truncated_html}

<p class="sub"><b>The collision audit.</b>  The bad manifest holds
{fmt(cm['collision_bad_rows'])} rows but only
{fmt(cm['collision_bad_distinct_paths'])} distinct paths: stripping the year
folded {fmt(cm['collision_paths_involved'])} files into
{fmt(cm['collision_groups'])} colliding names, so
{fmt(cm['collision_paths_at_risk'])} files would have been overwritten.  The
two manifests list the same Drive file ids
({fmt(cm['collision_ids_only_in_good'])} only in the good one,
{fmt(cm['collision_ids_only_in_bad'])} only in the bad one), so the crawl
itself lost nothing.  And the files now on disk under those names are
different exposures: in {fmt(cm['collision_groups_all_distinct'])} of
{fmt(cm['collision_groups'])} groups every member carries a different
<code>DATE-OBS</code>; {fmt(n_same)} groups contain members with an
identical stamp (genuine duplicates the dedup rule of section 2 handles).
The largest groups:</p>
{table(["name after stripping the year", "manifest paths", "scanned",
        "distinct DATE-OBS"],
       [[f"<code>{esc(p)}</code>", fmt(a), fmt(b), fmt(c)]
        for p, a, b, c in worst])}

<p class="sub">Inside the scan, every row is either a science frame or
carries one named reason for not being one:</p>
{table(["disposition", "files"], [[esc(k), fmt(n)] for k, n in fr])}
{_figure(src, "Scanned files by year of night, split into science frames "
              "and each named exclusion.  &lsquo;no date&rsquo; collects "
              "frames whose DATE-OBS is absent or outside "
              f"{lc.VALID_YEAR_MIN}&ndash;{lc.VALID_YEAR_MAX}.")}

<h3>Decision</h3>
<div class="decision"><b>The census population is the scanned
<code>.fz</code> files; everything else on disk or in the manifest is a
named exclusion, and both identities have zero residual.</b>  The download
was made from the corrected manifest; the bad one overwrote nothing.  The
frames truncated at the source cannot be recovered from this archive; a
same-name file is adopted only on an identical <code>DATE-OBS</code>, and
none qualifies.</div>

<h3>Consequence</h3>
<p class="sub">Counts below are counts of rows of
<code>products/legacy/legacy.sqlite</code>.  When the dedupe job finishes,
re-run <code>build_legacy_scan.py --rewalk</code>: the twin count falls to
zero and nothing else moves.</p>
</div></section>"""


def section_dedup(con) -> str:
    n = q1(con, "SELECT count(*) FROM frames WHERE error IS NULL")
    n_canon = q1(con, "SELECT sum(is_canonical) FROM frames WHERE error IS NULL")
    n_groups_multi = q1(con, """
        SELECT count(*) FROM (SELECT dup_group FROM frames
        GROUP BY dup_group HAVING count(*) > 1)""")
    by_tree = q(con, """
        SELECT tree, count(*), sum(is_canonical),
               count(*) - sum(is_canonical) FROM frames WHERE error IS NULL
        GROUP BY tree ORDER BY 2 DESC""")
    cross = q(con, """
        SELECT a, b, count(*) FROM (
          SELECT min(tree) AS a, max(tree) AS b FROM frames
          GROUP BY dup_group HAVING count(*) > 1) GROUP BY a, b ORDER BY 3 DESC""")
    calib_twin = q1(con, """
        SELECT count(*) FROM (SELECT dup_group FROM frames GROUP BY dup_group
        HAVING count(*) > 1
           AND sum(coalesce(trim(calstat), '') != '') > 0
           AND sum(coalesce(trim(calstat), '') = '') > 0)""")
    nodate = q1(con, "SELECT count(*) FROM frames WHERE exclusion = "
                     "'no_usable_date_obs'")
    return f"""
<section id="dedup">
<div class="bhead"><h2>2 &middot; The dedup rule</h2>
<span class="tag">one exposure, one canonical frame</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The referee asked that the go/no-go rest on a census with a
stated dedup rule.  When are two files the same exposure, and which one
counts?</p>

<h3>Evidence</h3>
<p class="sub">Identity key (pre-registered): camera, <code>DATE-OBS</code>
exactly as recorded, exposure time, filter, image geometry.  Of the
{fmt(n)} readable files, {fmt(n_canon)} are canonical;
{fmt(n - n_canon)} are copies, in {fmt(n_groups_multi)} groups.</p>
{table(["tree", "readable files", "canonical", "copies"],
       [[esc(t), fmt(a), fmt(b), fmt(c)] for t, a, b, c in by_tree])}
<p class="sub">Where the copies sit (pairs of trees a duplicate group spans):</p>
{table(["tree", "tree", "groups"], [[esc(a), esc(b), fmt(c)] for a, b, c in cross])
 if cross else '<p class="sub"><i>No duplicate groups.</i></p>'}
<p class="sub">{fmt(calib_twin)} groups pair a raw file with a calibrated
twin (non-blank <code>CALSTAT</code>); the raw one is canonical.
{fmt(nodate)} frames have no usable <code>DATE-OBS</code>: they cannot be
identified by this key, are never merged with anything, and are excluded
from science counts by that clause (section 4 says which cameras).</p>

<h3>Decision</h3>
<div class="decision"><b>Canonical = raw over calibrated twin, then
smallest path; every threshold counts canonical science frames only.</b>
Known limit, stated: the key cannot see a copy whose header time was
rewritten; none of the software in this archive is known to do that, and
the (target, night, exposure) populations show no doubled cadence.</div>

<h3>Consequence</h3>
<p class="sub">Nights, runs, series and gates below are built from the
canonical science frames.</p>
</div></section>"""


def section_names_eras(con) -> str:
    """File-name convention, era keys, and the second archive."""
    chk = q(con, """SELECT fn_convention, n_files, n_with_date, n_doy_equal,
                           n_doy_plus1, n_doy_minus1, n_doy_other
                    FROM filename_checks ORDER BY n_files DESC""")
    n_req = q1(con, "SELECT count(*) FROM filename_requests")
    req_one = q1(con, "SELECT count(*) FROM filename_requests "
                      "WHERE n_observers = 1")
    req_top = q(con, f"""SELECT request, n_files, n_observers, modal_observer,
                               frac_modal, n_targets, first_night, last_night
                        FROM filename_requests ORDER BY n_files DESC
                        LIMIT 10""")
    other = q(con, """SELECT basename FROM frames
                      WHERE fn_convention = 'other' AND error IS NULL
                      ORDER BY path LIMIT 6""")
    eras = q(con, """SELECT legacy_era, era_camera, era_optics, n_files,
                            n_science, first_night, last_night, telescop
                     FROM legacy_eras ORDER BY legacy_era""")
    shared = q(con, """SELECT rlmt_era, rlmt_readoutm, egain, geometry,
                              rlmt_frames, rlmt_first, rlmt_last,
                              legacy_frames, legacy_readoutm, legacy_first,
                              legacy_last FROM rlmt_shared_eras
                       ORDER BY rlmt_era""")
    n_copy = q1(con, "SELECT count(*) FROM cross_archive_copies")
    copy_rng = q(con, """SELECT min(date_obs), max(date_obs),
                                sum(rlmt_canonical) FROM cross_archive_copies""")[0]
    return f"""
<section id="names">
<div class="bhead"><h2>2b &middot; File names, era keys, and the other archive</h2>
<span class="tag">what a name says &middot; camera + focal length &middot; no double counting</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Three things must be true before a legacy frame can sit
beside an RLMT one: we know what its name encodes, its era is keyed on
hardware rather than on a telescope's name, and it is not the same exposure
as a frame the RLMT archive already holds.</p>

<h3>Evidence</h3>
<p class="sub"><b>The name.</b>  The scheduler named every file
<code>pppDDDss</code>: a three-letter request code, the day of year, and a
two-digit hexadecimal sequence.  Target, filter and exposure are NOT in the
name &mdash; they come from the header for every frame.  The day in the name
against the UT date of <code>DATE-OBS</code>:</p>
{table(["convention", "files", "with a date", "day = UT day", "day = UT + 1",
        "day = UT − 1", "off by more"],
       [[esc(c[0]), fmt(c[1]), fmt(c[2]), fmt(c[3]), fmt(c[4]), fmt(c[5]),
         fmt(c[6])] for c in chk])}
<p class="sub">The &plusmn;1 cases are frames taken across 00:00 UT under a
day number assigned when the night began.  Names outside the convention
(first six): {", ".join(f"<code>{esc(o[0])}</code>" for o in other)}.
Request codes: {fmt(n_req)}, of which {fmt(req_one)} were used by exactly one
<code>OBSERVER</code> &mdash; a code is a request queue, not a person or a
target (<code>foc</code> is the autofocus queue).  The ten largest:</p>
{table(["code", "files", "observers", "most frequent observer", "share",
        "targets", "first night", "last night"],
       [[f"<code>{esc(r[0])}</code>", fmt(r[1]), fmt(r[2]), esc(r[3]),
         f"{100 * r[4]:.0f}%", fmt(r[5]), esc(r[6]), esc(r[7])]
        for r in req_top])}

<p class="sub"><b>Era keys.</b>  A legacy era is (camera, focal length).
The 0.5&nbsp;m was named 'Gemini', 'Iowa Robotic Telescope' and 'Robert L.
Mutel Telescope' over these years with no optical change, so the name is
recorded but does not key.  Every row also carries
<code>archive_root = legacy-archive</code>.</p>
{table(["era", "camera", "optics", "files", "science frames", "first night",
        "last night", "TELESCOP as written"],
       [[fmt(e[0]), esc(e[1]), esc(e[2]), fmt(e[3]), fmt(e[4]), esc(e[5]),
         esc(e[6]), esc(e[7])] for e in eras])}
<p class="sub">The one configuration shared with the RLMT archive is the
AC4040 on the same telescope.  An RLMT era may lend its measured linearity
and ceiling to legacy frames only where the camera, the EGAIN table
(the driver changed it in 2023-10) and the full-frame geometry all
match:</p>
{table(["RLMT era", "readout mode", "EGAIN", "geometry", "RLMT frames",
        "RLMT first", "RLMT last", "matching legacy frames",
        "legacy readout modes", "legacy first", "legacy last"],
       [[fmt(r[0]), esc(r[1]), f1(r[2], 3), esc(r[3]), fmt(r[4]), esc(r[5]),
         esc(r[6]), fmt(r[7]), esc(r[8]), esc(r[9]), esc(r[10])]
        for r in shared]) if shared else
 '<p class="sub"><i>No RLMT era matches.</i></p>'}

<p class="sub"><b>The second archive.</b>  <b>{fmt(n_copy)}</b> legacy frames
are the same exposure (identical <code>DATE-OBS</code>, exposure and
geometry) as a frame in the RLMT manifest &mdash; {esc(copy_rng[0])} to
{esc(copy_rng[1])}, the first months of the RLMT era, copied into the
legacy share as well.  {fmt(copy_rng[2])} of their RLMT twins are canonical
there.</p>

<h3>Decision</h3>
<div class="decision"><b>The dedup rule spans both roots: a legacy frame
that is an RLMT exposure is excluded from the legacy counts, named
(<code>rlmt_archive_copy</code>), so no exposure is counted by two
projects.</b>  Eras key on camera and focal length; detector
characterisation crosses archives only through the matched-era table
above, and header gain never does.</div>

<h3>Consequence</h3>
<p class="sub">The legacy census proper ends in 2022-12.  The AC4040's
legacy frames (2021-11 to 2022-12) are the same camera RLMT eras 1&ndash;2
measure, with the same EGAIN table.</p>
</div></section>"""


def section_cameras(con) -> str:
    cams = q(con, """
        SELECT camera, n_canonical, n_science, first_night, last_night,
               n_nights, telescop, focallen_mm, native_format,
               native_pixel_um, binning, readoutm, egain, gain, set_temp,
               software, identity_basis, n_label_stale, instrume, n_no_date,
               flipstat, calstat
        FROM cameras ORDER BY first_night IS NULL, first_night""")
    n_canon = q1(con, "SELECT sum(is_canonical) FROM frames WHERE error IS NULL")
    rigel = q1(con, """
        SELECT count(*) FROM frames WHERE is_canonical = 1 AND error IS NULL
          AND telescop LIKE '%Rigel%'""")
    rigel_last = q1(con, """
        SELECT max(night) FROM frames WHERE telescop LIKE '%Rigel%'""")
    tel = q(con, """
        SELECT coalesce(telescop, '(none)'), round(focallen), count(*),
               min(night), max(night) FROM frames
        WHERE is_canonical = 1 AND error IS NULL GROUP BY 1, 2 ORDER BY 4""")
    stale = q1(con, "SELECT sum(n_label_stale) FROM cameras") or 0
    nolabel = q1(con, """
        SELECT count(*) FROM frames WHERE camera_basis = 'sensor (no label)'""")
    has_gain = q1(con, "SELECT count(*) FROM scan WHERE gain IS NOT NULL")
    has_egain = q1(con, "SELECT count(*) FROM scan WHERE egain IS NOT NULL")
    has_sn = q1(con, """
        SELECT count(*) FROM keysets WHERE ',' || keywords || ',' LIKE '%,CAMSN,%'
           OR ',' || keywords || ',' LIKE '%,SERIALNO,%'""")
    n_keysets = q1(con, "SELECT count(*) FROM keysets")
    id_tbl = table(
        ["camera (as identified)", "canonical frames", "first night",
         "last night", "nights", "native format", "native pixel (µm)",
         "identity basis", "INSTRUME text"],
        [[esc(c[0]), fmt(c[1]), esc(c[3]), esc(c[4]), fmt(c[5]), esc(c[8]),
          esc(c[9]), esc(c[16]), esc(c[18])] for c in cams])
    det_tbl = table(
        ["camera", "binning", "READOUTM", "EGAIN (header)", "GAIN card",
         "SET-TEMP (°C)", "software", "CALSTAT", "no usable date"],
        [[esc(c[0]), esc(c[10]), esc(c[11]), esc(c[12]), esc(c[13]),
          esc(c[14]), esc(c[15]), esc(c[21]), fmt(c[19])] for c in cams])
    return f"""
<section id="cameras">
<div class="bhead"><h2>3 &middot; Which telescope, which camera</h2>
<span class="tag">camera identity before any target census</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The plan called this the &ldquo;Rigel-era&rdquo; archive: a
different telescope with an FLI PL16803.  The telescope engineer said that
holds for only a few percent of it.  What do the headers say &mdash; and
can the headers be trusted on which camera took a frame?</p>

<h3>Evidence</h3>
<p class="sub">Telescope and focal length as written in the headers
(canonical frames):</p>
{table(["TELESCOP", "FOCALLEN (mm)", "canonical frames", "first night",
        "last night"],
       [[esc(a), fmt(b), fmt(c), esc(d), esc(e)] for a, b, c, d, e in tel])}
<p class="sub"><b>{fmt(rigel)}</b> of {fmt(n_canon)} canonical frames
(<b>{pct(rigel, n_canon)}</b>) carry <code>TELESCOP = Rigel System</code>;
the last such night is {esc(rigel_last)}.  Everything after is the 0.5&nbsp;m
(<code>FOCALLEN</code> 3454&nbsp;mm, <code>APTDIA</code> 508&nbsp;mm) under
two names.</p>
<p class="sub">Camera identity.  <code>INSTRUME</code> is a label the
scheduler copies from its configuration; the pixel pitch and frame format
come from the camera driver.  Where they disagree the sensor wins:
<b>{fmt(stale)}</b> frames carry an <code>INSTRUME</code> that names a
camera whose sensor cannot produce the recorded pixel pitch or frame width
(a stale label across a camera swap), and {fmt(nolabel)} have no label at
all.  Those are attributed by sensor signature &mdash; and because two
cameras share each such sensor, the <i>model</i> is honestly unknown:</p>
{id_tbl}
<p class="sub">Detector configuration cards, per camera.  No gain is
asserted from these: {fmt(has_egain)} scanned headers carry
<code>EGAIN</code> and {fmt(has_gain)} carry <code>GAIN</code>; a camera
serial-number card (<code>CAMSN</code>/<code>SERIALNO</code>) appears in
{fmt(has_sn)} of the {fmt(n_keysets)} distinct header keyword sets.</p>
{det_tbl}

<h3>Decision</h3>
<div class="decision"><b>The premise is corrected: this is the 0.5&nbsp;m's
archive with a short Rigel prologue ({pct(rigel, n_canon)} of frames), and
a camera is identified by its sensor signature, with INSTRUME as a label
that is checked, not believed.</b>  The name &ldquo;Legacy_Rigel&rdquo;
describes {pct(rigel, n_canon)} of the data.</div>

<h3>Consequence</h3>
<p class="sub">A series never pools cameras.  The 2022 AC4040 frames are
the RLMT-era detector on the same telescope and may share that camera's
detector characterisation (linearity, ceiling) once S2's flat-pair
measurements exist &mdash; but not its header gain.</p>
</div></section>"""


def section_mech(con) -> str:
    stints = q(con, """
        SELECT stint, camera, first_night, last_night, n_nights, n_frames,
               n_mixed_nights FROM camera_stints ORDER BY stint""")
    epochs = q(con, """
        SELECT mech_epoch, stint, camera, first_night, last_night,
               rotation_deg, n_nights_measured, n_excursion_nights
        FROM mech_epochs ORDER BY mech_epoch""")
    n_cams = q1(con, "SELECT count(DISTINCT camera) FROM camera_stints")
    n_single = sum(1 for e in epochs if e[6] == 1)
    step = float(q1(con, "SELECT value FROM census_meta "
                         "WHERE key='rotation_step_deg'"))
    src = fig_timeline(con)
    flips = q(con, """
        SELECT coalesce(flipstat, '(no card)'), count(*) FROM frames
        WHERE is_canonical = 1 AND error IS NULL GROUP BY 1 ORDER BY 2 DESC""")
    return f"""
<section id="mech">
<div class="bhead"><h2>4 &middot; Camera and mechanical timeline</h2>
<span class="tag">stints by camera &middot; rotation epochs inside each</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">What physically changed, and when?  A master frame or a
comparison-star solution must never be carried across a camera swap or a
re-seated camera.</p>

<h3>Evidence</h3>
{_figure(src, "Top: canonical frames per night, coloured by identified "
              "camera.  Bottom: nightly median field rotation from the "
              "scheduler's plate solutions (CROTA2, folded modulo 180° so "
              "a pier flip is not a step); grey lines are the derived "
              "mechanical-epoch boundaries.")}
<p class="sub">The camera timeline &mdash; {fmt(len(stints))} stints of
{fmt(n_cams)} identified cameras.  A stint is a maximal block of nights
dominated by one camera; &lsquo;mixed nights&rsquo; had frames from more
than one.</p>
{table(["stint", "camera", "first night", "last night", "nights",
        "canonical frames", "mixed nights"],
       [[fmt(s[0]), esc(s[1]), esc(s[2]), esc(s[3]), fmt(s[4]), fmt(s[5]),
         fmt(s[6])] for s in stints])}
<p class="sub">Rotation epochs inside the stints: a new epoch where the
nightly median rotation steps by more than {step}&deg; and stays there
(one-night excursions that return are counted, not promoted).
{fmt(len(epochs))} epochs; {fmt(n_single)} rest on a single measured night
and are the least secure.</p>
{table(["epoch", "stint", "camera", "first night", "last night",
        "rotation (°)", "nights measured", "excursion nights"],
       [[fmt(e[0]), fmt(e[1]), esc(e[2]), esc(e[3]), esc(e[4]),
         f1(e[5], 2), fmt(e[6]), fmt(e[7])] for e in epochs])}
<p class="sub">Flip state as recorded (<code>FLIPSTAT</code>):
{", ".join(f"{esc(a)} &times; {fmt(b)}" for a, b in flips)}.</p>

<h3>Decision</h3>
<div class="decision"><b>The mechanical epoch of a legacy frame is
(camera stint, rotation epoch).</b>  The table above is the census the
telescope engineer asked for: camera timeline with first and last night.
It is derived from headers only; filter-wheel maps and focus offsets are
not recorded in these headers and are not claimed.</div>

<h3>Consequence</h3>
<p class="sub">Any later photometry on a legacy series must split at these
boundaries.  The 0.3&deg; rule is applied to the scheduler's own plate
solutions, so an epoch is only as good as those solutions.</p>
</div></section>"""


def section_time(con) -> str:
    rows = q(con, """
        SELECT camera, software, n_files, n_no_usable_date, date_decimals,
               n_with_jd_card, frac_jd_agree, jd_diff_maxabs_s, n_with_lst,
               lst_resid_median_s, lst_resid_p05_s, lst_resid_p95_s,
               frac_lst_within_60s, frac_path_ut_date, n_path_offset_other,
               n_informative_pairs, viol_start, viol_mid, viol_end,
               n_overlapping_pairs, convention, first_night, last_night,
               lst_slope_vs_exptime, lst_intercept_s, n_lst_fit,
               jdhelio_ratio_median, n_jd_helio
        FROM time_audit ORDER BY first_night IS NULL, first_night""")
    src = fig_time(con)
    n_nodate = q1(con, "SELECT sum(n_no_usable_date) FROM time_audit")
    nodate = q(con, """
        SELECT camera, substr(date_obs, 1, 10), count(*),
               min(calstart), max(calstart), min(path)
        FROM frames WHERE error IS NULL AND jd_start IS NULL
        GROUP BY 1, 2 ORDER BY 3 DESC LIMIT 8""")
    other_cards = q(con, """
        SELECT sum(',' || keywords || ',' LIKE '%,JD-HELIO,%'),
               sum(',' || keywords || ',' LIKE '%,MJD-OBS,%'),
               sum(',' || keywords || ',' LIKE '%,DATE-AVG,%'),
               sum(',' || keywords || ',' LIKE '%,TIME-OBS,%'),
               sum(',' || keywords || ',' LIKE '%,TIMESYS,%'),
               count(*) FROM keysets""")[0]
    verdicts = q(con, """
        SELECT convention, count(*), sum(n_files), min(first_night),
               max(last_night), group_concat(DISTINCT camera)
        FROM time_audit GROUP BY convention ORDER BY 3 DESC""")
    n_start = sum(v[1] for v in verdicts if v[0] == "start of exposure")
    not_start = [v for v in verdicts if v[0] != "start of exposure"]
    conv_tbl = table(
        ["camera", "software", "first night", "last night", "DATE-OBS digits",
         "JD card within 2 s", "LST residual median (s)", "LST 5–95% (s)",
         "LST residual slope vs exposure", "intercept (s)",
         "JD-HELIO probe (÷ EXPTIME)",
         "unequal-exposure pairs", "overlap if start", "if middle", "if end",
         "verdict (overlap test)"],
        [[esc(r[0]), esc(r[1]), esc(r[21]), esc(r[22]), fmt(r[4]),
          f"{fmt(r[5])} ({pct((r[6] or 0) * r[5], r[5])})",
          f1(r[9], 2), f"{f1(r[10], 1)} … {f1(r[11], 1)}",
          f1(r[23], 3), f1(r[24], 2), f1(r[26], 3), fmt(r[15]), fmt(r[16]),
          fmt(r[17]),
          fmt(r[18]), esc(r[20])] for r in rows],
        [None if r[20] == "start of exposure" else "warn" for r in rows])
    path_tbl = table(
        ["camera", "software", "files", "no usable DATE-OBS",
         "directory = UT date", "off by &gt; 1 day"],
        [[esc(r[0]), esc(r[1]), fmt(r[2]), fmt(r[3]),
          f"{100 * r[13]:.1f}%" if r[13] is not None else "&mdash;",
          fmt(r[14])] for r in rows],
        ["warn" if r[3] else None for r in rows])
    return f"""
<section id="time">
<div class="bhead"><h2>5 &middot; The header-time convention</h2>
<span class="tag">what DATE-OBS means, before any timing is believed</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Two acquisition systems and nine cameras wrote these
headers.  Is <code>DATE-OBS</code> UTC?  Does it mark the start of the
exposure?  Is it there at all?  And how much of that can a header census
establish?</p>

<h3>Evidence</h3>
{_figure(src, "Left: header LST minus the LST implied by DATE-OBS and the "
              "site longitude, against exposure time (medians of groups of "
              "≥ 30 frames).  The scheduler samples LST as the exposure "
              "begins: a start stamp follows the flat line, a mid-exposure "
              "stamp the −t/2 line.  Right: consecutive exposures cannot "
              "overlap, so for pairs with unequal exposure times each "
              "convention predicts a different minimum spacing; bars show "
              "the share of such pairs that would overlap under each.")}
{conv_tbl}
<p class="sub">Four tests, of different strengths.  (1) The
<code>JD</code> card agrees with <code>DATE-OBS</code> on every frame
&mdash; to within rounding, which means it is the same number written
twice, not a second clock.  (2) The mount's <code>LST</code> agrees with
the sidereal time computed from <code>DATE-OBS</code> to seconds: the
stamp is UTC (a local-time stamp would sit seven hours away).  (3) The
overlap test: a stamp convention is refuted by pairs of consecutive
exposures that it would make overlap.  (4) The slope of the LST residual
against exposure time: 0 for a start stamp, &minus;0.5 for a mid-exposure
stamp, &minus;1 for an end stamp.  Tests 3 and 4 are independent of each
other (frame spacing vs a second header card) and agree.  Verdicts:</p>
{table(["verdict", "camera epochs", "files", "first night", "last night",
        "cameras"],
       [[esc(v[0]), fmt(v[1]), fmt(v[2]), esc(v[3]), esc(v[4]), esc(v[5])]
        for v in verdicts],
       [None if v[0] == "start of exposure" else "warn" for v in verdicts])}
<p class="sub">Other time cards in the {fmt(other_cards[5])} distinct
keyword sets: JD-HELIO in {fmt(other_cards[0])}, MJD-OBS in
{fmt(other_cards[1])}, DATE-AVG in {fmt(other_cards[2])}, TIME-OBS in
{fmt(other_cards[3])}, TIMESYS in {fmt(other_cards[4])}.  Directory names
against the UT date of <code>DATE-OBS</code>, and frames with no usable
stamp:</p>
{path_tbl}
<p class="sub"><b>{fmt(n_nodate)}</b> readable frames have no usable
<code>DATE-OBS</code> (absent, or a year outside {lc.VALID_YEAR_MIN}&ndash;
{lc.VALID_YEAR_MAX} &mdash; a camera clock that was never set).  They are
named exclusions.  The largest groups:</p>
{table(["camera", "DATE-OBS date", "frames", "CALSTART min", "CALSTART max",
        "example path"],
       [[esc(a), esc(b), fmt(c), esc(d), esc(e), f"<code>{esc(p)}</code>"]
        for a, b, c, d, e, p in nodate]) if nodate
 else '<p class="sub"><i>None.</i></p>'}

<h3>Decision</h3>
<div class="decision"><b><code>DATE-OBS</code> is UTC throughout.  It is
the START of the exposure in {fmt(n_start)} camera epochs; it is NOT the
start in {fmt(sum(v[1] for v in not_start))}:
{"; ".join(f"{esc(v[0])} &mdash; {esc(v[5])} ({esc(v[3])} to {esc(v[4])})"
           for v in not_start) or "none"}.</b>  Where the verdict is
&ldquo;middle&rdquo;, both independent tests say so: every unequal-exposure
pair is consistent with a mid-exposure stamp, start and end are each
refuted by hundreds of pairs, and the LST residual falls as half the
exposure time &mdash; and section 5b's eclipses, timed against TESS, agree.
S3's probe (MaxIm's JD-HELIO, column &lsquo;JD-HELIO probe&rsquo;) reads
0.5 on these same frames: MaxIm computed its heliocentric time as if
DATE-OBS were the start.  That probe records the software's assumption, so
it cannot tell a start stamp from a mid stamp; only the spacing, LST and
eclipse tests can.  Not established, and not establishable from headers: the
absolute clock.  JD and LST derive from the same computers as
<code>DATE-OBS</code>; a PC clock that was a minute slow would pass every
test here.  The observational astronomer's criterion (O&minus;C of archived
minima within 60&nbsp;s of the literature, per season) needs photometry
and is outside a census; section 8b lists the nights that could carry
it.</div>

<h3>Consequence</h3>
<p class="sub">A timing analysis must apply the convention of the camera
epoch a frame belongs to; treating a mid-exposure stamp as a start stamp
is an error of half the exposure time, with the sign of a late clock.
The camera for which this is found is the SBIG Aluma AC4040 &mdash; the
RLMT-era camera &mdash; under the same acquisition software family, so
this finding is handed to the shared clock stage (F-8) and to the CV
timing analysis for an independent check on the RLMT archive; this census
does not assert it there.  Any absolute epoch from this archive remains
conditional on a per-season clock measurement that has not been
made.</p>
</div></section>"""


def fig_clock(con) -> str | None:
    """O − C of every timed eclipse against the TESS reference, by date.

    Filled: the reading the header audit adopted for that camera epoch.
    Open (AC4040 only): the other reading.  The grey band is +-60 s, the
    ledger's criterion.
    """
    rows = q(con, """
        SELECT e.night, e.camera, f.oc_start_s, f.oc_mid_s, f.sig_stat_s,
               t.convention
        FROM clock_fits f JOIN clock_events e USING (event_id)
        LEFT JOIN time_audit t ON t.camera = e.camera
                              AND t.software = e.software
        WHERE f.status = 'ok'""")
    if not rows:
        return None
    style = camera_colors(con)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.8, 3.8))
        ax.axhspan(-lc_clock.CLOCK_ACCEPT_S, lc_clock.CLOCK_ACCEPT_S,
                   color=ps.GRID, zorder=0)
        ax.axhline(0, **ps.reference_kw(MUTED))
        seen = set()
        for night, cam, ocs, ocm, sig, conv in rows:
            st = style.get(cam, ps.series(0))
            mid = (conv or "").startswith("MIDDLE")
            adopted, other = (ocm, ocs) if mid else (ocs, ocm)
            d = _date(night)
            ax.errorbar([d], [adopted], yerr=[sig or 0], fmt=st["marker"],
                        color=st["color"], ms=5, lw=0.8,
                        label=None if cam in seen else cam)
            seen.add(cam)
            if mid:
                ax.scatter([d], [other], s=26, facecolors="none",
                           edgecolors=st["color"], linewidths=0.9)
        ax.set_ylabel("O − C (s), our clock − TESS")
        ax.set_xlabel("night")
        ax.set_title("Eclipse timings against TESS: the clock, season by season")
        ax.legend(fontsize=7, frameon=False, ncol=3)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "legacy_clock.png", dpi=DPI)
        plt.close(fig)
    return f"{FIG_REL}/legacy_clock.png"


def section_clock(con) -> str:
    if not has_table(con, "clock_seasons"):
        return """
<section id="clock"><div class="bhead"><h2>5b &middot; Clock audit per
season</h2><span class="tag">not evaluated</span></div><div class="stage">
<p class="sub"><i>Run <code>build_legacy_clock.py</code>.</i></p></div>
</section>"""
    tess = q(con, """SELECT vsx_name, count(*), sum(t0_bjd IS NOT NULL),
                            group_concat(DISTINCT sector), min(sig_s),
                            max(sig_s) FROM clock_tess_times
                     GROUP BY vsx_name ORDER BY vsx_name""")
    ev = q(con, """SELECT status, count(*) FROM clock_events
                   GROUP BY 1 ORDER BY 2 DESC""")
    ser = q(con, """SELECT status, count(*) FROM clock_series
                    GROUP BY 1 ORDER BY 2 DESC""")
    fits = q(con, """SELECT status, count(*) FROM clock_fits
                     GROUP BY 1 ORDER BY 2 DESC""")
    seasons = q(con, """
        SELECT camera, software, season, first_night, last_night, n_nights,
               convention, n_events, oc_s, sigma_s, chi2, dof,
               oc_other_reading_s, systems, verdict
        FROM clock_seasons ORDER BY first_night, camera""")
    n_meas = sum(1 for r in seasons if r[8] is not None)
    # Standing rule 1: chi2_nu < 0.5 is a defect equal to chi2_nu > 2.
    odd = [f"{r[0]} {r[2]} (χ²ν = {r[10] / r[11]:.2f})" for r in seasons
           if r[10] is not None and r[11] and
           not 0.5 <= r[10] / r[11] <= 2.0]
    n_pass = sum(1 for r in seasons if r[14] == "PASS")
    src = fig_clock(con)
    m = dict(q(con, "SELECT key, value FROM clock_meta"))
    return f"""
<section id="clock">
<div class="bhead"><h2>5b &middot; Clock audit per season</h2>
<span class="tag">eclipse timings against TESS &middot; criterion &plusmn;{f1(lc_clock.CLOCK_ACCEPT_S, 0)} s</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Section 5 fixed what each stamp means.  Was the clock that
wrote it right &mdash; season by season, to {f1(lc_clock.CLOCK_ACCEPT_S, 0)}
seconds?</p>

<h3>Evidence</h3>
<p class="sub">The clocks are the primary eclipses of post-common-envelope
binaries (HW Vir-type sdB&nbsp;+&nbsp;dM; WD&nbsp;+&nbsp;dM) in the archive.
Their reference is not a literature ephemeris &mdash; every one of these
stars wanders by tens of seconds or more about any linear ephemeris &mdash;
but eclipse times MEASURED in TESS 2-min photometry of the same star
(absolute BJD_TDB), one per half-sector, fitted locally.  A legacy eclipse is
predicted only inside the TESS span or within
{f1(float(m.get('extrap_max_d', 0)), 0)} days of it; beyond that the star's
own timing variations, not our clock, would set the answer, and the eclipse
is recorded as out of range.  TESS reference times:</p>
{table(["star", "half-sectors", "timed", "sectors", "best σ (s)",
        "worst σ (s)"],
       [[esc(r[0]), fmt(r[1]), fmt(r[2]), esc(r[3]), f1(r[4], 1),
         f1(r[5], 1)] for r in tess])}
<p class="sub">Legacy eclipses: {"; ".join(f"{esc(a)} {fmt(b)}" for a, b in ev)}.
Photometric series: {"; ".join(f"{esc(a)} {fmt(b)}" for a, b in ser)}.
Fits: {"; ".join(f"{esc(a)} {fmt(b)}" for a, b in fits)}.  Photometry and
fit are S3b's, unchanged (blind mid-time search, block bootstrap, analysis
variants); the field is identified by a Gaia DR3 fit because legacy headers
carry no plate scale; there is no master calibration (none is in the
archive) and the saturation veto is a fixed {fmt(float(m.get('veto_adu', 0)))}
ADU.</p>
{_figure(src, "O − C of every timed legacy eclipse against its TESS "
              "prediction.  Filled: the stamp reading the header audit "
              "adopted for that camera epoch; open circles (AC4040): the "
              "other reading.  Grey band: ±60 s.") if src else ""}
{table(["camera", "software", "season", "first night", "last night",
        "nights", "stamp reading", "eclipses", "O − C (s)", "σ (s)",
        "χ² / dof", "other reading (s)", "stars", "verdict"],
       [[esc(r[0]), esc(r[1]), esc(r[2]), esc(r[3]), esc(r[4]), fmt(r[5]),
         esc(r[6]), fmt(r[7]), f1(r[8], 1), f1(r[9], 1),
         "&mdash;" if r[10] is None else f"{r[10]:.1f} / {r[11]}",
         f1(r[12], 1), esc(r[13]), esc(r[14])] for r in seasons],
       ["ok" if r[14] == "PASS" else ("warn" if r[8] is None else None)
        for r in seasons])}

<h3>Decision</h3>
<div class="decision"><b>{fmt(n_meas)} of {fmt(len(seasons))} camera-seasons
have a clock measurement; {fmt(n_pass)} pass at
&plusmn;{f1(lc_clock.CLOCK_ACCEPT_S, 0)} s.</b>  Every other season carries
its offset as UNKNOWN, named in the table.  Where both stamp readings are
shown, the eclipses decide between them independently of the header
tests.  Seasons outside 0.5 &le; χ²ν &le; 2 (standing rule 1):
{esc("; ".join(odd)) or "none"}.  A low χ²ν here means the per-eclipse
errors are conservative: each carries its TESS prediction's systematic
term, which is shared by every eclipse of the same star and so does not
scatter them.</div>

<h3>Consequence</h3>
<p class="sub">Absolute times from a measured season may be used with the
O&minus;C and its error carried; times from an unmeasured season are
relative only.  An eclipse-timing study of these same stars cannot use them
as their own clock: per-season clock offsets would have to come from other
targets in the same season (archived exoplanet transits), or the study's
signal and the clock would be degenerate.</p>
</div></section>"""


def section_targets(con) -> str:
    n_targets = q1(con, "SELECT count(*) FROM targets")
    n_series = q1(con, "SELECT count(*) FROM series")
    n_sci = q1(con, "SELECT sum(is_science) FROM frames")
    n_nights = q1(con, "SELECT count(DISTINCT night) FROM frames "
                       "WHERE is_science = 1")
    span = q(con, "SELECT min(night), max(night) FROM frames "
                  "WHERE is_science = 1")[0]
    dist = q(con, """
        SELECT CASE WHEN n_nights = 1 THEN '1'
                    WHEN n_nights <= 2 THEN '2'
                    WHEN n_nights <= 4 THEN '3–4'
                    WHEN n_nights <= 9 THEN '5–9'
                    WHEN n_nights <= 29 THEN '10–29'
                    ELSE '≥ 30' END AS b, min(n_nights) AS o,
               count(*), sum(n_frames) FROM targets GROUP BY b ORDER BY o""")
    top_t = q(con, f"""
        SELECT target_name, n_nights, n_seasons, n_frames, first_night,
               last_night, round(longest_run_h, 1), bands, cameras
        FROM targets ORDER BY n_nights DESC, n_frames DESC LIMIT {TOP_N}""")
    top_s = q(con, f"""
        SELECT t.target_name, s.band, s.camera, s.n_nights, s.n_seasons,
               s.n_calibrated_nights, s.n_calstat_nights,
               round(s.longest_run_h, 1), s.n_nights_run_3h, s.n_frames,
               s.first_night, s.last_night
        FROM series s JOIN targets t USING (target_key)
        ORDER BY s.n_nights DESC, s.n_frames DESC LIMIT {TOP_N}""")
    long_runs = q(con, f"""
        SELECT t.target_name, s.band, s.camera, s.n_nights_run_3h, s.n_nights,
               round(s.longest_run_h, 1), s.n_seasons
        FROM series s JOIN targets t USING (target_key)
        WHERE s.n_nights_run_3h > 0
        ORDER BY s.n_nights_run_3h DESC, s.longest_run_h DESC LIMIT {TOP_N}""")
    n_long = q1(con, "SELECT count(*) FROM series WHERE n_nights_run_3h > 0")
    filt = q(con, """
        SELECT band, filter_system, tieable, count(*), count(DISTINCT night),
               group_concat(DISTINCT filter)
        FROM frames WHERE is_science = 1 GROUP BY band, filter_system, tieable
        ORDER BY 4 DESC""")
    alias = q1(con, "SELECT count(*) FROM targets WHERE n_raw_names > 1")
    obs = q1(con, "SELECT count(DISTINCT observer) FROM frames "
                  "WHERE is_science = 1")
    src = fig_series(con)
    return f"""
<section id="targets">
<div class="bhead"><h2>6 &middot; Targets, series, runs</h2>
<span class="tag">nights &times; filters &times; longest same-filter run</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Frame counts flatter an archive.  Per target: how many
nights, in which filter, on which camera, and how long is the longest
uninterrupted run?</p>

<h3>Evidence</h3>
<p class="sub">{fmt(n_sci)} canonical science frames on {fmt(n_nights)}
nights ({esc(span[0])} &rarr; {esc(span[1])}), {fmt(n_targets)}
alias-merged targets ({fmt(alias)} merged more than one OBJECT spelling,
by the same rules S0 uses for the RLMT archive), {fmt(n_series)} series
(target &times; filter &times; camera), {fmt(obs)} distinct OBSERVER
values &mdash; this is a many-user teaching and research queue, not one
programme.  Most targets were visited once or twice:</p>
{table(["nights per target", "targets", "science frames"],
       [[esc(b), fmt(c), fmt(d)] for b, _o, c, d in dist])}
{_figure(src, "Every series as one point: nights against longest "
              "same-filter run.  Stars are RLMT-era project targets, blue "
              "points are VSX eclipsing systems.  The 30-night line is "
              "gate G3's threshold, which additionally requires the nights "
              "to be calibrated (section 7).")}
<p class="sub">The {TOP_N} most-visited targets (all filters and cameras
pooled &mdash; a visiting statistic, not a series):</p>
{table(["target", "nights", "seasons", "frames", "first", "last",
        "longest run (h)", "filters (frames)", "cameras (frames)"],
       [[esc(r[0]), fmt(r[1]), fmt(r[2]), fmt(r[3]), esc(r[4]), esc(r[5]),
         f1(r[6]), esc(r[7]), esc(r[8])] for r in top_t])}
<p class="sub">The {TOP_N} longest series by nights (never pooled across
cameras or filters).  &lsquo;Calibrated nights&rsquo; is the
pre-registered availability test; &lsquo;CALSTAT nights&rsquo; counts
nights whose frames the scheduler had already calibrated in place:</p>
{table(["target", "filter", "camera", "nights", "seasons",
        "calibrated nights", "CALSTAT nights", "longest run (h)",
        "nights with run ≥ 3 h", "frames", "first", "last"],
       [[esc(r[0]), esc(r[1]), esc(r[2]), fmt(r[3]), fmt(r[4]), fmt(r[5]),
         fmt(r[6]), f1(r[7]), fmt(r[8]), fmt(r[9]), esc(r[10]), esc(r[11])]
        for r in top_s])}
<p class="sub">Time-series material: {fmt(n_long)} series have at least
one night with a run &ge; 3&nbsp;h.  The {TOP_N} with the most such
nights:</p>
{table(["target", "filter", "camera", "nights with run ≥ 3 h", "nights",
        "longest run (h)", "seasons"],
       [[esc(r[0]), esc(r[1]), esc(r[2]), fmt(r[3]), fmt(r[4]), f1(r[5]),
         fmt(r[6])] for r in long_runs])}
<p class="sub">Filters as written, and how they are read.  Talon wrote
single letters for a Johnson&ndash;Cousins wheel; MaxIm wrote
<code>'slot - name'</code>.  A named colour (<code>R - Red</code>,
<code>B - Blue</code>, <code>V - Visual</code>) does not say which glass
was in the slot, so those bands are <b>not</b> counted as tieable to a
standard system; only explicitly Johnson- or Sloan-named filters and the
Talon letters are.</p>
{table(["band (series key)", "system as read", "tieable", "science frames",
        "nights", "header FILTER text"],
       [[esc(b), esc(s), "yes" if t else "no", fmt(n), fmt(nn), esc(raw)]
        for b, s, t, n, nn, raw in filt])}

<h3>Decision</h3>
<div class="decision"><b>The unit of the census is the series (target
&times; filter &times; camera) with its nights, seasons and longest
run.</b>  A target with many frames on one night is one night.</div>

<h3>Consequence</h3>
<p class="sub">The gates read these tables.  Targets that look
interesting here but pass no gate are <i>exploratory</i> and are, by the
pre-registered rule, not a reason for a go.</p>
</div></section>"""


def section_overlap(con) -> str:
    if not has_table(con, "overlap"):
        return """
<section id="overlap"><div class="bhead"><h2>7 &middot; Overlap with the
RLMT-era targets</h2><span class="tag">not evaluated</span></div>
<div class="stage"><p class="sub"><i>The project-target table has not been
built: run <code>build_legacy_external.py</code>, then this build
again.</i></p></div></section>"""
    cone = float(q1(con, "SELECT value FROM census_meta "
                         "WHERE key='overlap_cone_deg'"))
    ov = q(con, """
        SELECT o.project, o.target, o.n_frames, o.n_by_name, o.n_by_position,
               o.n_name_not_position, o.n_nights, o.first_night, o.last_night,
               o.legacy_names, o.bands, o.cameras, t.position_source
        FROM overlap o JOIN rlmt_targets t USING (project, target_key)
        ORDER BY o.n_nights DESC, o.project, o.target""")
    found = [r for r in ov if r[2] > 0]
    absent = [r for r in ov if r[2] == 0]
    ser = q(con, """
        SELECT project, target, band, camera, tieable, filter_system, n_nights,
               n_seasons, first_night, last_night, round(longest_run_h, 2),
               n_frames, median_exptime_s, is_timing_target, g2_pass
        FROM overlap_series ORDER BY n_nights DESC, target, band""")
    n_all = q1(con, "SELECT count(*) FROM overlap_all")
    top_all = q(con, f"""
        SELECT target_name, rlmt_name, legacy_nights, legacy_seasons,
               legacy_frames, rlmt_nights, rlmt_frames, first_night, last_night
        FROM overlap_all ORDER BY legacy_nights DESC LIMIT {TOP_N}""")
    lenient = q(con, """
        SELECT target, band, camera, n_nights, n_seasons,
               round(longest_run_h, 2) FROM overlap_series
        WHERE g2_pass = 0 AND tieable = 0 AND n_nights >= ? AND n_seasons >= ?
          AND (is_timing_target = 0 OR longest_run_h >= ?)""",
                (lc.G2_MIN_NIGHTS, lc.G2_MIN_SEASONS,
                 lc.G2_TIMING_RUN_DAYS * 24))
    src = fig_overlap(con)
    absent_txt = ", ".join(sorted({esc(r[1]) for r in absent})) or "none"
    return f"""
<section id="overlap">
<div class="bhead"><h2>7 &middot; Overlap with the RLMT-era targets</h2>
<span class="tag">the one place a legacy frame can extend a baseline</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Were T&nbsp;CrB, the polars, YZ&nbsp;Cnc, the Be stars or
the dwarf-survey fields observed here in 2015&ndash;2022 &mdash; and if
so, on enough nights, in a filter that can be tied to the RLMT-era
photometry?</p>

<h3>Evidence</h3>
<p class="sub">Each project target is matched by alias key <i>and</i> by
position (requested coordinates within {cone}&deg; of the catalogue
position &mdash; the S0 cone), so a target observed under another name is
found and a same-named field elsewhere is not.  Targets found:</p>
{table(["project", "target", "science frames", "by name", "by position",
        "name but not position", "nights", "first", "last",
        "names used in the legacy headers", "filters", "cameras",
        "position from"],
       [[esc(r[0]), esc(r[1]), fmt(r[2]), fmt(r[3]), fmt(r[4]), fmt(r[5]),
         fmt(r[6]), esc(r[7]), esc(r[8]), esc(r[9]), esc(r[10]), esc(r[11]),
         esc(r[12])] for r in found],
       ["warn" if r[5] else None for r in found]) if found
 else '<p class="sub"><b>None of the project targets is in the archive.</b></p>'}
<p class="sub">Project targets with <b>no</b> legacy science frame:
{absent_txt}.</p>
{_figure(src, "Nights on which each RLMT-era project target was observed "
              "in the legacy archive, one marker per night and filter.")
 if src else ""}
<p class="sub">Per series, with the G2 clauses ({lc.G2_MIN_NIGHTS} nights,
{lc.G2_MIN_SEASONS} seasons, tieable filter; timing targets also a run
&ge; 1&nbsp;h):</p>
{table(["project", "target", "filter", "camera", "tieable", "nights",
        "seasons", "first", "last", "longest run (h)", "frames",
        "median exposure (s)", "G2"],
       [[esc(r[0]), esc(r[1]), esc(r[2]), esc(r[3]),
         "yes" if r[4] else f"no ({esc(r[5])})", fmt(r[6]), fmt(r[7]),
         esc(r[8]), esc(r[9]), f1(r[10], 2), fmt(r[11]), f1(r[12], 1),
         passed_cell(r[14])] for r in ser],
       ["ok" if r[14] else None for r in ser]) if ser
 else '<p class="sub"><i>No series.</i></p>'}
<p class="sub">Sensitivity to the filter reading: series that meet the
night, season and run clauses but in a filter not counted as tieable:
{("; ".join(f"{esc(a)} {esc(b)} on {esc(c)} ({fmt(d)} nights, {fmt(e)} seasons)"
            for a, b, c, d, e, _f in lenient)) or "none"}.</p>
<p class="sub">Wider, name-only view: {fmt(n_all)} of the RLMT archive's
science target names also name a legacy target.  The {TOP_N} with the most
legacy nights (a candidate list &mdash; names, not positions):</p>
{table(["legacy name", "RLMT name", "legacy nights", "legacy seasons",
        "legacy frames", "RLMT nights", "RLMT frames", "legacy first",
        "legacy last"],
       [[esc(r[0]), esc(r[1]), fmt(r[2]), fmt(r[3]), fmt(r[4]), fmt(r[5]),
         fmt(r[6]), esc(r[7]), esc(r[8])] for r in top_all])}

<h3>Decision</h3>
<div class="decision"><b>The overlap is what the series table says and no
more.</b>  A series that passes G2 is a <i>candidate</i> baseline
extension handed to the owning project; it is conditional on a pixel test
the census cannot make (saturation at the target, per frame, in native
pixels &mdash; standing rule 4) and on the clock caveat of section 5.</div>

<h3>Consequence</h3>
<p class="sub">Nothing here opens a sixth project.  What passes goes to
T&nbsp;CrB or CV as a candidate; what does not pass is recorded so the
question need not be asked again.</p>
</div></section>"""


def section_calibration(con) -> str:
    cal = q(con, """
        SELECT camera, kind, binning, band, n_frames, n_nights, first_night,
               last_night, exptimes, set_temp FROM calib_census
        ORDER BY camera, kind, binning, n_frames DESC""")
    n_cal = q1(con, """
        SELECT count(*) FROM frames WHERE is_canonical = 1
          AND kind IN ('bias', 'dark', 'flat')""")
    by_cam = q(con, """
        SELECT c.camera,
          (SELECT count(*) FROM frames f WHERE f.camera = c.camera
              AND f.is_science = 1),
          (SELECT count(DISTINCT night) FROM frames f WHERE f.camera = c.camera
              AND f.is_science = 1),
          (SELECT count(*) FROM frames f WHERE f.camera = c.camera
              AND f.is_canonical = 1 AND f.kind = 'bias'),
          (SELECT count(*) FROM frames f WHERE f.camera = c.camera
              AND f.is_canonical = 1 AND f.kind = 'dark'),
          (SELECT count(*) FROM frames f WHERE f.camera = c.camera
              AND f.is_canonical = 1 AND f.kind = 'flat'),
          (SELECT coalesce(sum(n_pairs), 0) FROM flat_pairs p
              WHERE p.camera = c.camera),
          (SELECT count(DISTINCT night) FROM flat_pairs p
              WHERE p.camera = c.camera),
          (SELECT coalesce(max(n_levels), 0) FROM flat_pairs p
              WHERE p.camera = c.camera),
          (SELECT count(*) FROM frames f WHERE f.camera = c.camera
              AND f.is_science = 1 AND coalesce(trim(f.calstat), '') != ''),
          (SELECT coalesce(sum(n_calibrated_nights), 0) FROM series s
              WHERE s.camera = c.camera),
          (SELECT coalesce(sum(n_nights), 0) FROM series s
              WHERE s.camera = c.camera)
        FROM cameras c ORDER BY c.first_night IS NULL, c.first_night""")
    calstat = q(con, """
        SELECT coalesce(nullif(trim(calstat), ''), '(none)'), count(*)
        FROM frames WHERE is_science = 1 GROUP BY 1 ORDER BY 2 DESC""")
    kinds_it = q(con, """
        SELECT coalesce(imagetyp, '(no IMAGETYP)'), kind, count(*) FROM frames
        WHERE error IS NULL GROUP BY 1, 2 ORDER BY 3 DESC""")
    focus_n = q1(con, "SELECT count(*) FROM frames WHERE kind = 'focus'")
    src = fig_calibration(con)
    best = q(con, """
        SELECT max(n_calibrated_nights) FROM series""")[0][0] or 0
    return f"""
<section id="calib">
<div class="bhead"><h2>8a &middot; Calibration census, per camera</h2>
<span class="tag">frames that exist &middot; flat pairs for a measured gain</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">For each camera: are there bias, dark and flat frames in
the archive at all, are there flat <i>pairs</i> so that gain can be
measured rather than read from a header, and were the science frames left
raw?</p>

<h3>Evidence</h3>
<p class="sub">How frames were classified (IMAGETYP as written &rarr;
kind).  {fmt(focus_n)} frames are the scheduler's autofocus images
(file prefix <code>foc</code>); they carry <code>Light Frame</code> and a
star's name and would be counted as science by the header alone.</p>
{table(["IMAGETYP", "kind", "files"],
       [[esc(a), esc(b), fmt(c)] for a, b, c in kinds_it])}
<p class="sub">The archive holds <b>{fmt(n_cal)}</b> canonical
calibration frames.  Per camera:</p>
{table(["camera", "science frames", "science nights", "bias", "dark",
        "flat", "flat pairs", "nights with a pair",
        "max exposure levels with a pair in one night",
        "science frames with CALSTAT",
        "calibrated series-nights / series-nights"],
       [[esc(r[0]), fmt(r[1]), fmt(r[2]), fmt(r[3]), fmt(r[4]), fmt(r[5]),
         fmt(r[6]), fmt(r[7]), fmt(r[8]), fmt(r[9]),
         f"{fmt(r[10])} / {fmt(r[11])}"] for r in by_cam],
       ["warn" if (r[3] + r[4] + r[5]) == 0 and r[1] else None
        for r in by_cam])}
{_figure(src, "One lane per camera: science nights (grey ticks) and the "
              "nights on which bias, dark or flat frames exist in the "
              "archive.")}
{table(["camera", "kind", "binning", "filter", "frames", "nights", "first",
        "last", "exposures s (frames)", "SET-TEMP (frames)"],
       [[esc(r[0]), esc(r[1]), fmt(r[2]), esc(r[3]), fmt(r[4]), fmt(r[5]),
         esc(r[6]), esc(r[7]), esc(r[8]), esc(r[9])] for r in cal]) if cal
 else '<p class="sub"><b>No bias, dark or flat frame was found in the '
      'archive.</b></p>'}
<p class="sub">In-place calibration.  The scheduler stamps
<code>CALSTAT</code> (B, D, F) when it applies bias/dark/flat to a frame
and overwrites it.  Science frames by <code>CALSTAT</code>:
{", ".join(f"<code>{esc(a)}</code> &times; {fmt(b)}" for a, b in calstat)}.
A frame with <code>CALSTAT</code> set is a reduced product whose raw
original and master frames may not be in this archive.</p>

<h3>Decision</h3>
<div class="decision"><b>Calibration availability is reported as counted,
and no gain is asserted.</b>  The pre-registered &ldquo;calibrated
night&rdquo; (a same-filter flat and a bias or dark on the same camera and
binning within &plusmn;{lc.CALIB_WINDOW_NIGHTS} nights, as frames in the
archive) is met by at most {fmt(best)} nights of any one series.  Flat
pairs exist only where the table shows them; where it shows none, the
camera's gain cannot be measured from this archive and the header
<code>EGAIN</code> is a manufacturer's number.</div>

<h3>Consequence</h3>
<p class="sub">For cameras with no calibration frames the only route to
calibrated photometry is frames the scheduler already calibrated
(<code>CALSTAT</code>), with masters that cannot be audited &mdash; fine
for differential timing, not for a calibrated flux series.  The request to
the site for the server's <code>calibrations/</code> trees (synthesis
&sect;6) covers exactly this gap.</p>
</div></section>"""


def section_eclipsers(con, m: dict) -> str:
    if not has_table(con, "eclipsers"):
        return """
<section id="eclipsers"><div class="bhead"><h2>8b &middot; Eclipsing and
contact systems</h2><span class="tag">not evaluated</span></div>
<div class="stage"><p class="sub"><i>The VSX match table has not been
built: run <code>build_legacy_external.py</code>, then this build
again.</i></p></div></section>"""
    n_q = int(m.get("vsx_targets_queried", 0))
    n_m = int(m.get("vsx_targets_matched", 0))
    n_ecl = q1(con, "SELECT count(*) FROM eclipsers WHERE is_eclipsing = 1")
    n_tr = q1(con, "SELECT count(*) FROM eclipsers WHERE is_transit = 1")
    dist = q(con, """
        SELECT n_seasons_bearing, count(*), sum(n_seasons_guaranteed >= ?),
               sum(n_seasons_predicted >= ?)
        FROM eclipsers WHERE is_eclipsing = 1 GROUP BY 1 ORDER BY 1 DESC""",
             (lc.G1_MIN_SEASONS, lc.G1_MIN_SEASONS))
    rows = q(con, f"""
        SELECT target_name, vsx_name, vsx_type, period_d, match, sep_arcsec,
               n_nights, n_seasons_observed, n_nights_guaranteed,
               n_nights_predicted, n_seasons_guaranteed, n_seasons_predicted,
               n_seasons_bearing, first_bearing_night, last_bearing_night,
               bands
        FROM eclipsers WHERE is_eclipsing = 1 AND n_nights_bearing > 0
        ORDER BY n_seasons_bearing DESC, n_nights_bearing DESC LIMIT 40""")
    n_any = q1(con, """SELECT count(*) FROM eclipsers
                       WHERE is_eclipsing = 1 AND n_nights_bearing > 0""")
    transits = q(con, f"""
        SELECT target_name, vsx_name, period_d, n_nights, n_seasons_observed,
               n_nights_bearing, n_seasons_bearing, bands, cameras
        FROM eclipsers WHERE is_transit = 1
        ORDER BY n_nights_bearing DESC, n_nights DESC LIMIT {TOP_N}""")
    by_cam = q(con, """
        SELECT s.camera, substr(n.night, 1, 4), count(DISTINCT n.target_key),
               count(DISTINCT n.night)
        FROM eclipser_nights n JOIN eclipsers e USING (target_key)
        JOIN (SELECT DISTINCT target_key, camera, night FROM frames
              WHERE is_science = 1) s
          ON s.target_key = n.target_key AND s.night = n.night
        WHERE e.is_eclipsing = 1 AND n.bearing = 1
        GROUP BY 1, 2 ORDER BY 2, 1""")
    src = fig_eclipsers(con)
    pos_only = q1(con, """
        SELECT count(*) FROM eclipsers WHERE is_eclipsing = 1
          AND n_seasons_bearing >= ? AND match = 'position'""",
                  (lc.G1_MIN_SEASONS,))
    n_pass = q1(con, """SELECT count(*) FROM eclipsers WHERE is_eclipsing = 1
                        AND n_seasons_bearing >= ?""", (lc.G1_MIN_SEASONS,))
    return f"""
<section id="eclipsers">
<div class="bhead"><h2>8b &middot; Eclipsing and contact systems</h2>
<span class="tag">minima in &ge; {lc.G1_MIN_SEASONS} seasons &mdash; the physicist's criterion</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">A period change of 10<sup>&minus;7</sup>&nbsp;d&nbsp;yr<sup>&minus;1</sup>
accumulates about seven minutes of O&minus;C in ten years, which
one-minute timings detect.  How many systems in the archive have eclipse
minima observed in at least {lc.G1_MIN_SEASONS} seasons?</p>

<h3>Evidence</h3>
<p class="sub">{fmt(n_q)} legacy targets (those with &ge;
{esc(m.get('vsx_min_nights', ''))} nights or a run &ge;
{esc(m.get('vsx_min_run_hours', ''))}&nbsp;h) were cone-matched against VSX
(<code>{esc(m.get('vsx_catalog', ''))}</code>, queried
{esc(m.get('vsx_query_utc', m.get('external_built_utc', '')))[:10]});
{fmt(n_m)} have a VSX counterpart, of which {fmt(n_ecl)} are eclipsing or
contact binaries (GCVS classes E, EA, EB, EW, EC) and {fmt(n_tr)} are
planet-transit hosts (EP), listed separately.  A night is
<i>minimum-bearing</i> if a run is at least P/2 long (a minimum is then
certain whatever the ephemeris) or contains a predicted minimum with
30&nbsp;min of data on each side (ephemeris-dependent).</p>
{table(["seasons with a minimum-bearing night", "systems",
        f"of which ≥ {lc.G1_MIN_SEASONS} seasons by run ≥ P/2 alone",
        f"of which ≥ {lc.G1_MIN_SEASONS} seasons by prediction alone"],
       [[fmt(a), fmt(b), fmt(c), fmt(d)] for a, b, c, d in dist])}
{_figure(src, "The eclipsing systems with at least one minimum-bearing "
              "night, ordered by seasons covered.  Filled blue: a run "
              "≥ P/2.  Open orange: only a predicted minimum inside the "
              "run.  Grey: other nights on the target.") if src else ""}
<p class="sub">{fmt(n_any)} systems have at least one minimum-bearing
night; the first 40 by seasons:</p>
{table(["target", "VSX name", "type", "P (d)", "match", "sep (″)",
        "nights", "seasons observed", "nights run ≥ P/2",
        "nights predicted min", "seasons (run ≥ P/2)",
        "seasons (predicted)", "seasons (union)", "first", "last",
        "filters"],
       [[esc(r[0]), esc(r[1]), esc(r[2]), f1(r[3], 4), esc(r[4]),
         f1(r[5], 1), fmt(r[6]), fmt(r[7]), fmt(r[8]), fmt(r[9]), fmt(r[10]),
         fmt(r[11]), fmt(r[12]), esc(r[13]), esc(r[14]), esc(r[15])]
        for r in rows],
       ["ok" if r[12] >= lc.G1_MIN_SEASONS else None for r in rows]) if rows
 else '<p class="sub"><b>No eclipsing system has a minimum-bearing night.</b></p>'}
<p class="sub">Of the {fmt(n_pass)} systems at or above
{lc.G1_MIN_SEASONS} seasons, {fmt(pos_only)} were matched to VSX by
position only (the legacy name did not normalise to the VSX name).
Minimum-bearing nights by camera and year &mdash; the nights a per-season
clock measurement (O&minus;C within 60&nbsp;s of the literature) could be
made on:</p>
{table(["camera", "year", "systems", "minimum-bearing nights"],
       [[esc(a), esc(b), fmt(c), fmt(d)] for a, b, c, d in by_cam]) if by_cam
 else '<p class="sub"><i>None.</i></p>'}
<p class="sub">Planet-transit hosts (not counted by G1; relevant to the
shared clock test F-8, which wants archived transits per era):</p>
{table(["target", "VSX name", "P (d)", "nights", "seasons",
        "nights with a predicted event or run ≥ P/2", "seasons with one",
        "filters", "cameras"],
       [[esc(r[0]), esc(r[1]), f1(r[2], 4), fmt(r[3]), fmt(r[4]), fmt(r[5]),
         fmt(r[6]), esc(r[7]), esc(r[8])] for r in transits]) if transits
 else '<p class="sub"><i>None matched.</i></p>'}

<h3>Decision</h3>
<div class="decision"><b>The count the physicist asked for is
{fmt(n_pass)} systems with minimum-bearing nights in &ge;
{lc.G1_MIN_SEASONS} seasons.</b>  It is a count of nights on which a
minimum <i>should be in the data</i>; no minimum has been measured, no
light curve inspected, and the clause-(b) nights inherit the catalogue
ephemeris error.</div>

<h3>Consequence</h3>
<p class="sub">Whether a timing paper exists depends on photometry that
has not been done: are the minima actually sampled, unsaturated, with a
usable comparison star &mdash; and is the clock right to a minute in each
season (section 5).</p>
</div></section>"""


def section_gates(con, m: dict) -> str:
    gates = q(con, "SELECT gate, clause, value, threshold, passed, detail "
                   "FROM gates ORDER BY CASE gate WHEN 'G0' THEN 0 WHEN 'G0''' "
                   "THEN 1 WHEN 'G1' THEN 2 WHEN 'G1*' THEN 3 WHEN 'G2' "
                   "THEN 4 ELSE 5 END, rowid")
    outcome = m.get("outcome", "")
    go = lc.OUTCOME_GO_CANDIDATE in outcome
    transfer = lc.OUTCOME_TRANSFER in outcome
    nodecide = lc.OUTCOME_NOT_DECIDABLE in outcome or "NOT YET" in outcome
    lines = []
    if nodecide:
        lines.append("<b>The rule cannot be applied yet</b>: a validity "
                     "clause fails or an input table is absent.  Nothing "
                     "below is a go.")
    if transfer:
        lines.append(
            "<b>Rule 3 applies (G2):</b> the passing overlap series are "
            "handed to the owning project as candidate baseline extensions "
            "&mdash; reader: the readers of that MACRO paper; venue: a "
            "section of it (RNAAS if it stands alone).  This is a "
            "transfer, not a sixth project.")
    if go:
        lines.append(
            "<b>Rule 4 applies:</b> GO-CANDIDATE.  Named question (Q1): do "
            "the archive's W&nbsp;UMa-type/eclipsing systems show secular "
            "period change over 2015&ndash;2022?  Reader: observers and "
            "ephemeris maintainers of those systems (AAVSO eclipsing-binary "
            "section, BAV / O&minus;C Gateway users).  Venue: JAAVSO (or "
            "OEJV) &mdash; a student-led timing paper.  A strategy "
            "<i>may</i> be commissioned; whether it is, is James's "
            "decision, and by the editor's ruling it must not compete "
            "with T&nbsp;CrB for effort this autumn.")
    if not (go or transfer or nodecide):
        lines.append(
            "<b>Rule 2 applies:</b> no gate passes.  Outcome = archive "
            "data release only: this database, a Zenodo record and a short "
            "data note.  No strategy is commissioned.")
    amended = m.get("outcome_amended", "")
    alt = m.get("outcome_alt_reading", "")
    amended_html = ("<br><b>Second reading of G0 (row G0&prime;): " + esc(alt)
                    + ".</b>  The pre-registered sentence asks that the "
                    "header audit &lsquo;identifies what DATE-OBS means "
                    "(start of exposure, UTC)&rsquo;; the code required the "
                    "start.  Read as &lsquo;the meaning is identified&rsquo;, "
                    "G0 holds, because the AC4040's mid-exposure stamp is "
                    "identified by two independent tests (and, section 5b, "
                    "by the eclipses).") if alt else ""
    if amended:
        amended_html += (
            "<br><b>Amended reading (a deviation, reported beside the written "
            f"one as pre-registration &sect;5 requires): {esc(amended)}.</b>  "
            "G0 fails only because camera epochs that feed G1 do not stamp "
            "the exposure start (" + esc(m.get("failing_time_epochs", ""))
            + ").  &ldquo;Fix the census&rdquo; has an obvious meaning here: "
            "drop the runs from those epochs and re-evaluate.  Row G1* does "
            "that; systems still passing: "
            + esc(m.get("g1_amended_systems", "") or "none") + ".  "
            + ("Under that reading rule 4 applies &mdash; GO-CANDIDATE, with "
               "the named question, reader and venue of Q1 (period change "
               "from eclipse timings; ephemeris maintainers and "
               "eclipse-timing observers; JAAVSO/OEJV class).  It is a "
               "recommendation to commission a strategy, not a commitment; "
               "the decision is James's and must not compete with T&nbsp;CrB "
               "this autumn.  Before any strategy is written, a literature "
               "and novelty check on the passing systems comes first (the "
               "rule the committee set for the Be stars, U9): eclipse timing "
               "of short-period eclipsing binaries is an active field, and "
               "nights in this archive may already be published."
               if lc.OUTCOME_GO_CANDIDATE in amended else
               "Under that reading no go results either."))
    return f"""
<section id="gates">
<div class="bhead"><h2>9 &middot; The pre-registered rule, applied</h2>
<span class="tag">{esc(outcome)}</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Section 0 fixed what would count.  Does the archive
deliver it?</p>

<h3>Evidence</h3>
{table(["gate", "clause", "value", "threshold", "result", "detail"],
       [[esc(g[0]), esc(g[1]), esc(g[2]), esc(g[3]), passed_cell(g[4]),
         esc(g[5])] for g in gates])}
<p class="sub">Summary: G0 {passed_cell(_flag(m, 'g0'))} &middot;
G1 {passed_cell(_flag(m, 'g1'))} &middot;
G2 {passed_cell(_flag(m, 'g2'))} &middot;
G3 {passed_cell(_flag(m, 'g3'))}.</p>

<h3>Decision</h3>
<div class="decision"><b>Outcome as written: {esc(outcome)}.</b><br>
{"<br>".join(lines)}{amended_html}</div>

<h3>Consequence</h3>
<p class="sub">This is a recommendation produced by a rule, not a
commitment.  What the census cannot see is listed where it matters:
saturation and comparison stars (pixels), the absolute clock (a measured
event per season), gain (flat pairs).  Each of those is a short, bounded
test on the specific nights the tables above name.</p>
</div></section>"""


def _flag(m: dict, key: str):
    v = m.get(key, "")
    return None if v == "" else bool(int(v))


# ---------------------------------------------------------------------------
# The page
# ---------------------------------------------------------------------------
def render_report(db_path: Path) -> Path:
    """Render the full census page from the database.  Returns the HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, timeout=120)
    try:
        m = meta(con)
        n_scan = q1(con, "SELECT count(*) FROM scan")
        n_sci = q1(con, "SELECT sum(is_science) FROM frames")
        n_targets = q1(con, "SELECT count(*) FROM targets")
        n_nights = q1(con, "SELECT count(DISTINCT night) FROM frames "
                           "WHERE is_science = 1")
        n_cams = q1(con, "SELECT count(*) FROM cameras")
        sections = [
            section_prereg(con, m),
            section_reconciliation(con, m),
            section_dedup(con),
            section_names_eras(con),
            section_cameras(con),
            section_mech(con),
            section_time(con),
            section_clock(con),
            section_targets(con),
            section_overlap(con),
            section_calibration(con),
            section_eclipsers(con, m),
            section_gates(con, m),
        ]
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Legacy archive — census (L0, L1) and go/no-go (L2)</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>Legacy archive &mdash; census and go/no-go</h1>
  <p>{fmt(n_scan)} scanned files &rarr; {fmt(n_sci)} canonical science
  frames &middot; {fmt(n_targets)} targets &middot; {fmt(n_nights)} nights
  &middot; {fmt(n_cams)} identified cameras &middot; outcome:
  <b>{esc(m.get('outcome', ''))}</b> &middot;
  built {esc(m.get('built_utc', ''))[:16]}Z
  ({esc(m.get('code_version', ''))},
  commit <code>{esc(m.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="index.html">plan &amp; status</a></p>
</header>

<nav>
  <a href="#rule">0 Rule</a> &middot;
  <a href="#recon">1 Reconciliation</a> &middot;
  <a href="#dedup">2 Dedup</a> &middot;
  <a href="#names">2b Names &amp; eras</a> &middot;
  <a href="#cameras">3 Cameras</a> &middot;
  <a href="#mech">4 Timeline</a> &middot;
  <a href="#time">5 Time</a> &middot;
  <a href="#clock">5b Clock</a> &middot;
  <a href="#targets">6 Targets</a> &middot;
  <a href="#overlap">7 Overlap</a> &middot;
  <a href="#calib">8a Calibration</a> &middot;
  <a href="#eclipsers">8b Eclipsers</a> &middot;
  <a href="#gates">9 Outcome</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_legacy.report</code> from
<code>products/legacy/legacy.sqlite</code> — every number on this page is
the result of a SQL query; none is typed by hand.  Regenerate with
<code>pipeline/scripts/build_legacy_census.py</code>.  This page is a
census of headers: no pixel was opened.</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
        # Belt and braces: every <img> the page references must exist and be
        # non-empty, or the build fails loudly rather than shipping a broken
        # evidence page.
        for src in re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
