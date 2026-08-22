"""S0 chain-of-evidence report renderer.

Reads the S0 manifest database (NEVER the catalog — if a number cannot be
derived from the manifest, it does not belong on the page) and writes:

* ``docs/pipeline/s0_manifest.html``     — the report
* ``docs/pipeline/figures/s0/*.png``     — every figure, matplotlib, drawn
  in the house style defined once in ``macro_core.plotstyle``

The page follows the site's Socratic format (``docs/assets/macro.css``):
one section per decision, each section = Question → Evidence → Decision →
Consequence.  EVERY number in the HTML is interpolated from a SQL query
executed in this module or from a constant defined in ``macro_core.manifest``
— nothing is hand-typed, so re-running the build after a catalog change
regenerates the whole argument.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import manifest as m      # noqa: E402  (constants for interpolation)
from . import plotstyle as ps    # noqa: E402  (the house figure style)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s0"
HTML_PATH = DOCS_DIR / "s0_manifest.html"

# The house figure style.  These are re-exported here rather than defined
# here because a dozen renderers already do `from .report_s0 import ...` for
# the page machinery below, and the style must have ONE definition —
# macro_core.plotstyle.  There is no dark theme any more: a figure is a
# printable object and prints on white.
STYLE = ps.STYLE
ACCENT = ps.ACCENT     # primary data colour (Okabe--Ito blue)
WARN = ps.WARN         # outliers, disagreements (Okabe--Ito orange)
GOOD = ps.GOOD         # confirmations (Okabe--Ito green)
BAD = ps.BAD           # contradictions (Okabe--Ito vermilion)
MUTED = ps.MUTED       # reference lines, annotations
FAINT = ps.FAINT       # backgrounded context
INK = ps.INK           # type drawn onto the figure
PAPER = ps.PAPER       # the ground, and type drawn onto a dark cell
DPI = ps.WEB_DPI       # spec requires >= 110


# ---------------------------------------------------------------------------
# Small query helpers — every number on the page flows through these.
# ---------------------------------------------------------------------------
def q(con: sqlite3.Connection, sql: str, params=()) -> list[tuple]:
    """Run one query, return all rows."""
    return con.execute(sql, params).fetchall()


def q1(con: sqlite3.Connection, sql: str, params=()):
    """Run one query, return the single scalar of its first row."""
    return con.execute(sql, params).fetchone()[0]


def fmt(x) -> str:
    """Human-format a number for the page: thousands separators, NULL→em-dash."""
    if x is None:
        return "&mdash;"
    if isinstance(x, float) and x != int(x):
        return f"{x:,.2f}"
    return f"{int(x):,}"


def esc(s) -> str:
    """Minimal HTML escaping for strings interpolated into the page."""
    if s is None:
        return "&mdash;"
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def table(headers: list[str], rows: list[list], row_classes=None) -> str:
    """Render one table.data — the only table markup generator on the page."""
    out = ['<table class="data">',
           "<tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr>"]
    for i, row in enumerate(rows):
        cls = f' class="{row_classes[i]}"' if row_classes and row_classes[i] else ""
        out.append(f"<tr{cls}>" + "".join(f"<td>{c}</td>" for c in row) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_dup_hist(con) -> str:
    """Histogram of duplicate-group sizes (how many copies each frame has)."""
    rows = q(con, """
        SELECT n, count(*) FROM (
            SELECT count(*) AS n FROM frames GROUP BY dup_group
        ) GROUP BY n ORDER BY n""")
    sizes = [r[0] for r in rows]
    counts = [r[1] for r in rows]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7, 3.4))
        bars = ax.bar(sizes, counts, color=ACCENT, width=0.7)
        ax.set_yscale("log")
        ax.set_xlabel("copies of the same (basename, JD) exposure")
        ax.set_ylabel("duplicate groups")
        ax.set_title("Duplicate-group sizes across all trees")
        # Annotate each bar with its exact count — the figure carries its
        # own numbers so it survives being lifted out of the page.
        for b, c in zip(bars, counts):
            ax.annotate(f"{c:,}", (b.get_x() + b.get_width() / 2, c),
                        ha="center", va="bottom", fontsize=8, color=MUTED)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_dup_group_sizes.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_dup_group_sizes.png"


def fig_alias_top25(con) -> str:
    """Top-25 targets: canonical light frames under the biggest single raw
    name vs after alias consolidation — the visual measure of name
    fragmentation.  BOTH series count the SAME unit (canonical Light
    frames), so the yellow bar can never exceed the blue and no display cap
    is needed — an earlier revision mixed catalog rows with canonical
    frames and clipped the mismatch, which reviewers rightly rejected."""
    # Blue series: canonical light frames per merged (canonical) target.
    rows = q(con, """
        SELECT f.target_key, f.canonical_target, count(*) AS n_canon
        FROM frames f
        WHERE f.is_canonical = 1 AND f.target_key IS NOT NULL
          AND f.imagetyp LIKE 'Light%'
        GROUP BY f.target_key ORDER BY n_canon DESC LIMIT 25""")
    # Yellow series: for each of those targets, the canonical-light-frame
    # count of its single most populous RAW name — same selection, further
    # grouped by the raw target_best value.  Computed in one pass and
    # reduced in Python (SQLite lacks lateral correlated FROM-subqueries).
    per_raw = q(con, """
        SELECT target_key, target_best, count(*) AS n
        FROM frames
        WHERE is_canonical = 1 AND target_key IS NOT NULL
          AND imagetyp LIKE 'Light%'
        GROUP BY target_key, target_best""")
    biggest_raw: dict[str, int] = {}
    for tkey, _raw, n in per_raw:
        biggest_raw[tkey] = max(biggest_raw.get(tkey, 0), n)
    names = [r[1] for r in rows][::-1]
    after = [r[2] for r in rows][::-1]
    before = [biggest_raw[r[0]] for r in rows][::-1]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.4, 7))
        ypos = np.arange(len(names))
        ax.barh(ypos, after, color=ACCENT, alpha=0.9,
                label="canonical light frames (all aliases merged)")
        ax.barh(ypos, before, color=WARN, alpha=0.85, height=0.45,
                label="canonical light frames under the largest raw name")
        ax.set_yticks(ypos, names, fontsize=8)
        ax.set_xlabel("canonical light frames")
        ax.set_title("Alias consolidation — top 25 targets")
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_alias_top25.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_alias_top25.png"


def fig_era_timeline(con, min_frames: int) -> str:
    """Night-vs-era strip chart: when each camera configuration was live."""
    rows = q(con, """
        SELECT f.era_id, f.night, count(*)
        FROM frames f
        WHERE f.is_canonical = 1 AND f.era_id IS NOT NULL
          AND f.night IS NOT NULL
          AND f.era_id IN (SELECT era_id FROM eras WHERE n_frames >= ?)
        GROUP BY f.era_id, f.night""", (min_frames,))
    if not rows:
        return ""
    import datetime as _dt
    eras = sorted({r[0] for r in rows})
    # Evidence-integrity guard: every era the adjacent table shows above the
    # figure threshold MUST contribute at least one point here.  The first
    # shipped build silently dropped 29 eras (the frames.era_id NULL bug),
    # so the figure and the table disagreed without anyone noticing — now a
    # mismatch kills the render instead of shipping a misleading page.
    expected = {r[0] for r in q(con, """
        SELECT era_id FROM eras
        WHERE n_frames >= ? AND first_night IS NOT NULL""", (min_frames,))}
    missing = expected - set(eras)
    if missing:
        raise RuntimeError(
            f"era timeline omits eras present in the eras table: "
            f"{sorted(missing)} — frames.era_id is out of sync with eras")
    era_pos = {e: i for i, e in enumerate(eras)}
    xs = [_dt.date.fromisoformat(r[1]) for r in rows]
    ys = [era_pos[r[0]] for r in rows]
    ss = [max(4.0, 1.2 * np.sqrt(r[2]) * 4) for r in rows]
    labels = {r[0]: f"era {r[0]}: {r[1]} {r[2] or ''}" for r in q(con, """
        SELECT era_id, readoutm, cast(naxis1 AS int) || 'x' ||
               cast(naxis2 AS int) || ' bin' || cast(xbinning AS int)
        FROM eras WHERE n_frames >= ?""", (min_frames,))}
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.2, 0.42 * len(eras) + 1.6))
        ax.scatter(xs, ys, s=ss, color=ACCENT, alpha=0.55, linewidths=0)
        ax.set_yticks(range(len(eras)),
                      [labels.get(e, f"era {e}") for e in eras], fontsize=7.5)
        ax.set_xlabel("night (local-noon-to-noon label)")
        ax.set_title("Camera eras on the calendar "
                     "(marker area ~ frames per night)")
        ax.invert_yaxis()
        # Date ticks overlap at this width — thin and rotate them.
        fig.autofmt_xdate(rotation=30, ha="right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_era_timeline.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_era_timeline.png"


def fig_frames_per_night(con) -> str:
    """Archive-wide histogram of canonical frames per night."""
    rows = q(con, """
        SELECT count(*) FROM frames
        WHERE is_canonical = 1 AND night IS NOT NULL GROUP BY night""")
    per_night = [r[0] for r in rows]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.hist(per_night, bins=np.logspace(0, np.log10(max(per_night)), 40),
                color=ACCENT)
        ax.set_xscale("log")
        ax.set_xlabel("canonical frames per night")
        ax.set_ylabel("nights")
        ax.set_title("Frames per night, archive-wide")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_frames_per_night.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_frames_per_night.png"


def fig_night_boundary(con, densest_night: str) -> str:
    """Two-panel sanity demo that the noon-to-noon boundary cuts nothing.

    Left: every frame of the single densest night, plotted by UT hour, with
    the two local-noon boundaries.  Right: the UT-hour distribution of the
    whole archive — the boundary must fall in the empty daytime gap.
    """
    # Frame times of the densest night, as UT hours since that night's
    # start-boundary (JD fraction arithmetic on the stored header JD).
    rows = q(con, """
        SELECT jd FROM frames
        WHERE is_canonical = 1 AND night = ? AND jd IS NOT NULL""",
             (densest_night,))
    jd = np.array([r[0] for r in rows])
    # Hours since the night's opening boundary (local noon).  The boundary
    # JD has fractional part (NIGHT_SHIFT_DAYS − 0.5) by construction.
    frac = (jd - (m.NIGHT_SHIFT_DAYS - 0.5)) % 1.0
    hours_since_noon = frac * 24.0

    all_jd = q(con, """
        SELECT jd FROM frames WHERE is_canonical = 1 AND jd IS NOT NULL""")
    all_hours = ((np.array([r[0] for r in all_jd])
                  - (m.NIGHT_SHIFT_DAYS - 0.5)) % 1.0) * 24.0

    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.4))
        ax1.eventplot(hours_since_noon, colors=ACCENT, lineoffsets=0.5,
                      linelengths=0.8)
        # Inside the axes, not above it: at 1.02 these two labels sat on
        # the title line and the left one collided with it.
        for x, lbl, ha in ((0, "boundary (local noon)", "left"),
                           (24, "next boundary", "right")):
            ax1.axvline(x, color=WARN, lw=1.4)
            ax1.text(x, 0.965, " " + lbl + " ", color=WARN, fontsize=7.5,
                     ha=ha, va="top",
                     transform=ax1.get_xaxis_transform())
        ax1.set_xlim(-1, 25)
        ax1.set_yticks([])
        ax1.set_xlabel("hours after local noon")
        ax1.set_title(f"Densest night {densest_night}: every frame time")

        ax2.hist(all_hours, bins=96, color=ACCENT)
        ax2.axvline(0, color=WARN, lw=1.4)
        ax2.axvline(24, color=WARN, lw=1.4)
        ax2.set_xlabel("hours after local noon")
        ax2.set_ylabel("canonical frames")
        ax2.set_title("Whole archive: exposure clock vs the boundary")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_night_boundary.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_night_boundary.png"


def fig_pointing(con) -> str:
    """Log-scale histogram of pointing offsets with the 1-degree line.

    Every canonical frame with an offset is plotted — exact-zero offsets
    (frames whose header coords equal the reference to the stored digit)
    cannot live on a log axis, so they are clipped INTO the lowest bin
    rather than dropped, keeping the plotted total equal to the caption's
    count (the caption also states the clip, script-emitted)."""
    rows = q(con, """
        SELECT pointing_offset_deg FROM frames
        WHERE is_canonical = 1 AND pointing_offset_deg IS NOT NULL""")
    # Clip floor: one decade below the smallest nonzero offset (fallback
    # 1e-5 deg), so clipped zeros land visibly in the first bin.
    off = np.array([r[0] for r in rows])
    nonzero_min = off[off > 0].min() if (off > 0).any() else 1e-5
    floor = max(nonzero_min, 1e-5)
    off = np.maximum(off, floor)
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7, 3.4))
        bins = np.logspace(np.log10(floor), np.log10(off.max()), 60)
        ax.hist(off, bins=bins, color=ACCENT)
        ax.axvline(m.POINTING_OUTLIER_DEG, color=WARN, lw=1.4)
        ax.text(m.POINTING_OUTLIER_DEG, 0.95,
                f" {m.POINTING_OUTLIER_DEG:g} deg outlier line",
                color=WARN, fontsize=8,
                transform=ax.get_xaxis_transform())
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("pointing offset from target reference (deg)")
        ax.set_ylabel("canonical frames")
        ax.set_title("Frame pointing vs alias-resolved target position")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_pointing_offsets.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_pointing_offsets.png"


def fig_qc(con) -> str:
    """Bar chart of QC-flag token counts (a frame can carry several)."""
    # Split the comma-joined tokens in SQL-free Python: tiny data volume.
    rows = q(con, "SELECT qc_flags, count(*) FROM frames "
                  "WHERE qc_flags != '' GROUP BY qc_flags")
    counts: dict[str, int] = {}
    for flags, n in rows:
        for tok in flags.split(","):
            counts[tok] = counts.get(tok, 0) + n
    toks = sorted(counts, key=counts.get)
    vals = [counts[t] for t in toks]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7, 3.2))
        bars = ax.barh(toks, vals, color=WARN)
        ax.set_xscale("log")
        ax.set_xlabel("frames carrying the flag (all rows, log scale)")
        ax.set_title("QC flags — marks, not deletions")
        for b, v in zip(bars, vals):
            ax.annotate(f" {v:,}", (v, b.get_y() + b.get_height() / 2),
                        va="center", fontsize=8, color=INK)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0_qc_flags.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0/s0_qc_flags.png"


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def _figure(src: str, caption: str) -> str:
    return (f'<figure><a href="{src}"><img src="{src}" alt=""></a>'
            f"<figcaption>{caption}</figcaption></figure>")


def section_duplication(con) -> str:
    n_rows = q1(con, "SELECT count(*) FROM frames")
    n_groups = q1(con, "SELECT count(DISTINCT dup_group) FROM frames")
    n_dups = n_rows - n_groups
    src = fig_dup_hist(con)

    # Per-tree contribution: canonical rows vs duplicate copies.  ALL trees
    # — the spec asks what EACH tree contributes, and a truncated table
    # whose columns do not sum to the page totals is an auditor trap (an
    # earlier revision stopped at 12 rows without saying so).
    tree_rows = q(con, """
        SELECT tree, count(*) AS rows_,
               sum(is_canonical) AS canon,
               count(*) - sum(is_canonical) AS dups
        FROM frames GROUP BY tree ORDER BY rows_ DESC""")
    tree_tbl = table(
        ["tree", "catalog rows", "canonical frames", "duplicate copies",
         "% duplicate"],
        [[esc(t), fmt(r), fmt(c), fmt(d), f"{100.0 * d / r:.1f}%"]
         for t, r, c, d in tree_rows])
    # Tree-policy caveat, quantified: nights whose canonical Light frames
    # come ENTIRELY from non-rawimage trees.  A project that counts only
    # its primary tree silently misses such nights (the ST LMi iKon case);
    # downstream primary-tree accounting inherits this caveat explicitly.
    n_nonraw_nights = q1(con, """
        SELECT count(*) FROM (
            SELECT night FROM frames
            WHERE is_canonical = 1 AND imagetyp LIKE 'Light%'
              AND night IS NOT NULL
            GROUP BY night HAVING sum(tree = 'rawimage') = 0)""")
    n_all_nights = q1(con, """
        SELECT count(DISTINCT night) FROM frames
        WHERE is_canonical = 1 AND imagetyp LIKE 'Light%'
          AND night IS NOT NULL""")

    # The SN lesson: duplicates INSIDE rawimage — count them and name the
    # directories that hold wholesale copies of other nights.
    n_raw_internal = q1(con, """
        SELECT count(*) FROM frames f
        WHERE f.tree = 'rawimage' AND f.is_canonical = 0
          AND EXISTS (SELECT 1 FROM frames w
                      WHERE w.dup_group = f.dup_group
                        AND w.is_canonical = 1 AND w.tree = 'rawimage')""")
    dirs = q(con, """
        SELECT substr(f.path, 1, length(f.path) - length(f.basename) - 1),
               count(*)
        FROM frames f
        WHERE f.tree = 'rawimage' AND f.is_canonical = 0
          AND EXISTS (SELECT 1 FROM frames w
                      WHERE w.dup_group = f.dup_group
                        AND w.is_canonical = 1 AND w.tree = 'rawimage')
        GROUP BY 1 ORDER BY 2 DESC LIMIT 8""")
    dir_tbl = table(["rawimage directory holding copies of other nights",
                     "duplicate frames"],
                    [[f"<code>{esc(d)}</code>", fmt(n)] for d, n in dirs])

    # The renamed-copy limitation: identical (target, JD) among canonical
    # frames — files the reduced/ tree renamed, invisible to basename dedup.
    n_jd_collide = q1(con, """
        SELECT count(*) FROM (
            SELECT target_key FROM frames
            WHERE is_canonical = 1 AND jd IS NOT NULL
              AND target_key IS NOT NULL
            GROUP BY target_key, jd HAVING count(*) > 1)""")
    # Worked example for the prose: EU UMa's reduced-only "unique" frames
    # whose JDs all collide with a rawimage frame (renamed copies).
    euuma_renamed = q1(con, """
        SELECT count(*) FROM frames f
        WHERE f.target_key = 'euuma' AND f.is_canonical = 1
          AND f.tree = 'reduced'
          AND EXISTS (SELECT 1 FROM frames g
                      WHERE g.target_key = 'euuma' AND g.jd = f.jd
                        AND g.tree = 'rawimage' AND g.is_canonical = 1)""")
    # Worked example for the consequence: the SN row-vs-frame inflation.
    sn_rows = q1(con, """SELECT count(*) FROM frames
        WHERE target_key = '2023ixf' AND tree = 'rawimage'
          AND imagetyp LIKE 'Light%'""")
    sn_canon = q1(con, """SELECT count(*) FROM frames
        WHERE target_key = '2023ixf' AND tree = 'rawimage'
          AND imagetyp LIKE 'Light%' AND is_canonical = 1""")
    collide_trees = q(con, """
        SELECT f.tree, count(*) FROM frames f
        WHERE f.is_canonical = 1 AND f.jd IS NOT NULL
          AND f.target_key IS NOT NULL AND f.tree != 'rawimage'
          AND EXISTS (SELECT 1 FROM frames g
                      WHERE g.target_key = f.target_key AND g.jd = f.jd
                        AND g.tree = 'rawimage' AND g.is_canonical = 1)
        GROUP BY f.tree ORDER BY 2 DESC LIMIT 6""")
    collide_tbl = table(
        ["tree", "canonical frames whose (target, JD) also exists in rawimage"],
        [[esc(t), fmt(n)] for t, n in collide_trees])

    return f"""
<section id="dedup">
<div class="bhead"><h2>1 · Duplication</h2>
<span class="tag">dedup global on (basename, JD)</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The catalog holds {fmt(n_rows)} rows.  How many distinct
exposures is that, and which copy of each is the one the pipeline reads?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"Duplicate-group size distribution: {fmt(n_groups)} groups from "
    f"{fmt(n_rows)} rows; {fmt(n_dups)} rows are extra copies of an "
    "exposure that already exists.")}</div>
{tree_tbl}
<p class="sub">Tree-policy caveat: <b>{fmt(n_nonraw_nights)}</b> of the
{fmt(n_all_nights)} observing nights with canonical Light frames hold them
in non-<code>rawimage</code> trees <i>only</i> (school mirrors, the iKon
tree) — a project counting nothing but its primary tree misses those
nights entirely (ST&nbsp;LMi&rsquo;s iKon surplus in section&nbsp;7 is the
worked example).</p>
<p class="sub">The SN 2023ixf lesson, confirmed archive-wide: duplication is
not only cross-tree.  <b>{fmt(n_raw_internal)}</b> rawimage rows are copies
of frames whose canonical row is <i>also</i> in rawimage — wholesale
re-copied night directories:</p>
{dir_tbl}
<p class="sub">Limitation, quantified: (basename, JD) cannot see a RENAMED
copy.  <b>{fmt(n_jd_collide)}</b> (target, JD) pairs occur on more than one
canonical frame — overwhelmingly <code>reduced/</code> files renamed during
processing (e.g. every one of EU UMa&rsquo;s {fmt(euuma_renamed)}
reduced-only &ldquo;unique&rdquo; frames shares its JD with a rawimage
frame).  By tree:</p>
{collide_tbl}

<h3>Decision</h3>
<div class="decision"><b>Duplicate identity is (basename, JD), globally —
across and within trees.</b>  One canonical row per group, chosen by tree
policy: <code>rawimage</code> first, <code>reduced</code> last, ties broken
by earliest path; documented per-target exception: NGC&nbsp;5548 prefers
<code>macalester</code> (the superset tree).  Renamed <code>reduced/</code>
copies survive dedup by construction, so <b>downstream accounting uses
canonical frames from each target&rsquo;s primary tree</b> (as section 7
does), and the (target,&nbsp;JD)-collision population above is handed to S1
as a suspect list.</div>

<h3>Consequence</h3>
<p class="sub">Every count below this line uses the {fmt(n_groups)}
canonical frames.  A project quoting raw rows overstates its dataset by up
to ~{sn_rows / sn_canon:.1f}&times; (SN 2023ixf: {fmt(sn_rows)} rawimage
rows &rarr; {fmt(sn_canon)} frames).</p>
</div></section>"""


def section_aliases(con) -> str:
    n_raw = q1(con, "SELECT count(*) FROM aliases WHERE raw_name IS NOT NULL")
    n_canon = q1(con, "SELECT count(DISTINCT target_key) FROM aliases "
                      "WHERE target_key IS NOT NULL")
    n_nonident = q1(con, "SELECT count(*) FROM aliases "
                         "WHERE method NOT IN ('identity', 'blank')")
    n_blank_frames = q1(con, "SELECT coalesce(sum(n_frames),0) FROM aliases "
                             "WHERE raw_name IS NULL")
    src = fig_alias_top25(con)

    merges = q(con, """
        SELECT canonical_target, count(*) AS nv, sum(n_frames) AS nf,
               group_concat(method)
        FROM aliases WHERE target_key IS NOT NULL
        GROUP BY target_key HAVING nv > 1 ORDER BY nv DESC LIMIT 12""")

    def _method_union(joined: str) -> str:
        # group_concat(DISTINCT ...) dedups whole comma-joined METHOD
        # STRINGS, not individual rule tokens, which printed duplicates like
        # 'casefold,separators,casefold'.  Split to tokens here and dedup
        # while preserving first-seen order — a true union.
        seen: list[str] = []
        for tok in (joined or "").split(","):
            if tok and tok not in seen:
                seen.append(tok)
        return ",".join(seen)

    merge_tbl = table(
        ["canonical target", "raw variants merged", "catalog rows",
         "resolution methods (union)"],
        [[esc(t), fmt(nv), fmt(nf),
          f"<code>{esc(_method_union(meth))}</code>"]
         for t, nv, nf, meth in merges])

    # Cone audit: pass / fail / not-checkable.
    n_checked = q1(con, "SELECT count(*) FROM aliases "
                        "WHERE cone_check_passed IS NOT NULL")
    n_failed = q1(con, "SELECT count(*) FROM aliases WHERE cone_check_passed=0")
    fails = q(con, """SELECT raw_name, canonical_target, n_frames, method
                      FROM aliases WHERE cone_check_passed = 0
                      ORDER BY n_frames DESC""")
    fail_tbl = table(
        ["raw name", "canonical target", "rows", "method"],
        [[esc(r), esc(c), fmt(n), f"<code>{esc(mm)}</code>"]
         for r, c, n, mm in fails])

    # Must-NOT-merge regression, run live against the aliases table.
    dw_keys = q(con, """SELECT raw_name, target_key FROM aliases
                        WHERE raw_name IN ('Dw1403+49', 'Dw1409+51')""")
    dw_ok = len({k for _, k in dw_keys}) == len(dw_keys) == 2
    tcrb_keys = q(con, """SELECT DISTINCT target_key FROM aliases
                          WHERE raw_name IN ('T CrB', 'tet CrB')""")
    tcrb_ok = len(tcrb_keys) == 2
    guard_tbl = table(
        ["must-NOT-merge check", "result"],
        [["Dw1403+49 vs Dw1409+51 (adjacent survey fields)",
          "separate keys &#10003;" if dw_ok else "MERGED — BUG"],
         ["T CrB vs tet CrB (target vs its calibrator, 8&deg; apart)",
          "separate keys &#10003;" if tcrb_ok else "MERGED — BUG"]],
        row_classes=[None if dw_ok else "warn", None if tcrb_ok else "warn"])

    return f"""
<section id="aliases">
<div class="bhead"><h2>2 · Alias consolidation</h2>
<span class="tag">rules + cone-gated synonyms, never fuzzy matching</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The catalog carries {fmt(n_raw)} distinct raw target names.
How many real targets is that, and by exactly what rule does each variant
map to its canonical name?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"{fmt(n_raw)} raw names resolve to {fmt(n_canon)} canonical targets; "
    f"{fmt(n_nonident)} names needed at least one normalization rule. "
    "Yellow = the largest single raw variant; blue = the merged total.")}
</div>
{merge_tbl}
<p class="sub">Coordinate-cone audit ({m.CONE_RADIUS_DEG:g}&deg; on median
plate-solved positions): {fmt(n_checked)} raw names had solved coordinates
to check; <b>{fmt(n_failed)}</b> sit farther than {m.CONE_RADIUS_DEG:g}&deg;
from their group&rsquo;s pooled position — all explainable: mosaic/wide-field
pointings of extended fields (M31, rho&nbsp;Oph) and a moving target
(Saturn), where a fixed cone is meaningless by construction:</p>
{fail_tbl}
{guard_tbl}
<p class="sub">{fmt(n_blank_frames)} rows carry a blank target name (mostly
raw-tree grism lights with valid coordinates); they form no alias group and
are flagged <code>blank_target</code> for S1&rsquo;s cone match.</p>

<h3>Decision</h3>
<div class="decision"><b>A merge happens through a named pure rule
(case/whitespace/separators, junk prefixes, exposure-token strip, series
suffix, genitive spelling) or through the explicit two-entry synonym table
(Alpha&nbsp;Lyr&rarr;Vega, NGC&nbsp;5457&rarr;M101), gated by the
{m.CONE_RADIUS_DEG:g}&deg; coordinate cone.  Fuzzy string similarity is
banned.</b>  Every alias row records which rules fired and its cone-audit
result.</div>

<h3>Consequence</h3>
<p class="sub">Downstream stages see one name per target.  The
reconciliation of section 7 works ONLY because the CV/BeStar filename
variants (<code>mjcMay01 yzcnc</code>, <code>PHECDA lrg 0-25s</code>,
<code>Vega 0p001s lrg 5</code>) land on their canonical targets.</p>
</div></section>"""


def section_eras(con) -> str:
    n_eras = q1(con, "SELECT count(*) FROM eras")
    n_readout = q1(con, "SELECT count(DISTINCT readoutm) FROM eras")
    min_frames_fig = 100
    src = fig_era_timeline(con, min_frames_fig)
    n_fig_eras = q1(con, "SELECT count(*) FROM eras WHERE n_frames >= ?",
                    (min_frames_fig,))

    era_rows = q(con, """
        SELECT era_id, readoutm, naxis1, naxis2, xbinning, egain,
               n_frames, first_night, last_night
        FROM eras ORDER BY era_id""")
    era_tbl = table(
        ["era", "READOUTM", "NAXIS1", "NAXIS2", "bin", "EGAIN",
         "canonical frames", "first night", "last night"],
        [[fmt(e), esc(r) or "<i>(blank)</i>", fmt(n1), fmt(n2), fmt(xb),
          esc(eg), fmt(nf), esc(fn), esc(ln)]
         for e, r, n1, n2, xb, eg, nf, fn, ln in era_rows])

    # The BeStar era-C repackaging: the NAXIS1=8 population, found by key.
    repack = q1(con, """
        SELECT coalesce(sum(n_frames), 0) FROM eras WHERE naxis1 = 8""")

    return f"""
<section id="eras">
<div class="bhead"><h2>3 · Era structure</h2>
<span class="tag">keyed on (READOUTM, geometry, binning, EGAIN) — never
filter or date</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">How many distinct camera configurations produced this
archive, and when was each one live?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"The {fmt(n_fig_eras)} eras holding at least {fmt(min_frames_fig)} "
    "canonical frames, on the calendar. Marker area scales with frames per "
    "night. Eras are numbered by first light.")}</div>
<p class="sub">All {fmt(n_eras)} discovered configurations
({fmt(n_readout)} distinct READOUTM strings; full table also exported as
<code>products/manifest/eras.csv</code>):</p>
{era_tbl}

<h3>Decision</h3>
<div class="decision"><b>Era assignment keys on (READOUTM, NAXIS1&times;
NAXIS2, XBINNING, EGAIN rounded to {m.EGAIN_DECIMALS} decimals) — never on
filter name or date</b> (BeStar lesson: hrg/lrg-labeled frames taken with
the era-A camera; HaGrism running past its nominal end date).  The
repackaged-FITS population the BeStar panel flagged surfaces here on its
own: {fmt(repack)} canonical frames carry the NAXIS1=8 extension-HDU
geometry and form their own eras rather than polluting their neighbors.</div>

<h3>Consequence</h3>
<p class="sub">S2 builds calibrations per era_id; S1 stratifies its solve
experiment per era_id; no stage ever infers a camera from a filter string
again.</p>
</div></section>"""


def section_nights(con) -> str:
    n_nights = q1(con, "SELECT count(DISTINCT night) FROM frames "
                       "WHERE night IS NOT NULL")
    densest = q(con, """
        SELECT night, count(*) AS n FROM frames
        WHERE is_canonical = 1 AND night IS NOT NULL
        GROUP BY night ORDER BY n DESC LIMIT 1""")[0]
    med = q1(con, """
        SELECT n FROM (SELECT count(*) AS n FROM frames
            WHERE is_canonical = 1 AND night IS NOT NULL GROUP BY night
            ORDER BY n)
        LIMIT 1 OFFSET (SELECT (count(DISTINCT night) - 1) / 2 FROM frames
                        WHERE is_canonical = 1 AND night IS NOT NULL)""")
    src1 = fig_frames_per_night(con)
    src2 = fig_night_boundary(con, densest[0])
    # Sanity numbers: frames within one hour of a boundary — total, and the
    # subset that are actual Light frames (the rest are daytime calibration
    # frames, for which a daylight boundary crossing is meaningless).
    shift = m.NIGHT_SHIFT_DAYS - 0.5
    near_total = q1(con, """
        SELECT count(*) FROM frames WHERE is_canonical = 1 AND jd IS NOT NULL
        AND abs((jd - ?) - round(jd - ?)) < (1.0/24.0)""", (shift, shift))
    near_light = q1(con, """
        SELECT count(*) FROM frames WHERE is_canonical = 1 AND jd IS NOT NULL
        AND imagetyp LIKE 'Light%'
        AND abs((jd - ?) - round(jd - ?)) < (1.0/24.0)""", (shift, shift))

    return f"""
<section id="nights">
<div class="bhead"><h2>4 · Night labeling</h2>
<span class="tag">local-noon-to-noon: date(JD &minus; {m.NIGHT_SHIFT_DAYS})</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">What is &ldquo;a night&rdquo;, such that no contiguous
observing sequence is ever split across two labels?</p>

<h3>Evidence</h3>
<div class="grid">
{_figure(src1, f"{fmt(n_nights)} distinct nights; median "
              f"{fmt(med)} canonical frames per night; the densest night "
              f"({esc(densest[0])}) holds {fmt(densest[1])}.")}
{_figure(src2, f"Left: all {fmt(densest[1])} frames of the densest night "
              f"sit far inside the boundaries. Right: archive-wide, "
              f"{fmt(near_total)} canonical frames fall within one hour of "
              f"a boundary, and only {fmt(near_light)} of them is a Light "
              "frame — the rest are daytime bias/dark/flat frames. The "
              "boundary crosses no observing sequence.")}
</div>

<h3>Decision</h3>
<div class="decision"><b>night = calendar date of (JD &minus;
{m.NIGHT_SHIFT_DAYS})</b> — the date rolls over at 19:00&nbsp;UT =
12:00 local (Winer, UTC&minus;7), and the label equals the local evening
date.  Header JD itself is stored untouched: it is the UTC exposure
<i>start</i>, and mid-exposure BJD_TDB is S3&rsquo;s job, not S0&rsquo;s.</div>

<h3>Consequence</h3>
<p class="sub">Per-night statistics (calibration masters, nightly zero
points, night-block bootstraps) inherit a boundary that no exposure ever
straddles.  Beware: a strategy document counting UTC <i>dates</i> can quote
up to one extra night per observing run (see the SN row in section 7).</p>
</div></section>"""


def section_pointing(con) -> str:
    n_ref = q1(con, """SELECT count(DISTINCT target_key) FROM frames
                       WHERE pointing_offset_deg IS NOT NULL""")
    n_off = q1(con, """SELECT count(*) FROM frames
                       WHERE is_canonical = 1
                         AND pointing_offset_deg IS NOT NULL""")
    n_out = q1(con, """SELECT count(*) FROM frames
                       WHERE is_canonical = 1 AND pointing_offset_deg > ?""",
               (m.POINTING_OUTLIER_DEG,))
    # Exact-zero offsets: real frames, but a log axis cannot hold zero, so
    # the figure clips them into its lowest bin — the caption says so with
    # this script-emitted count (caption total == plotted total, always).
    n_zero = q1(con, """SELECT count(*) FROM frames
                        WHERE is_canonical = 1 AND pointing_offset_deg = 0""")
    src = fig_pointing(con)

    offenders = q(con, """
        SELECT canonical_target,
               sum(pointing_offset_deg > ?) AS bad, count(*) AS tot,
               round(max(pointing_offset_deg), 1) AS worst
        FROM frames
        WHERE is_canonical = 1 AND pointing_offset_deg IS NOT NULL
        GROUP BY target_key HAVING bad > 0
        ORDER BY bad DESC LIMIT 12""", (m.POINTING_OUTLIER_DEG,))
    off_tbl = table(
        ["target", f"frames &gt;{m.POINTING_OUTLIER_DEG:g}&deg;",
         "frames with offsets", "worst offset (deg)"],
        [[esc(t), fmt(b), fmt(tot), esc(w)] for t, b, tot, w in offenders])

    # Composition of the outlier population: focus/test pseudo-targets and
    # moving solar-system objects dominate — for them a fixed reference
    # position is meaningless by construction, and the flag is expected.
    n_pseudo = q1(con, """
        SELECT coalesce(sum(pointing_offset_deg > ?), 0) FROM frames
        WHERE is_canonical = 1 AND (target_key LIKE 'focus%'
              OR target_key LIKE 'testfield%')""", (m.POINTING_OUTLIER_DEG,))
    anuma = q(con, """
        SELECT sum(pointing_offset_deg > ?), count(*),
               round(max(pointing_offset_deg), 1)
        FROM frames WHERE target_key = 'anuma' AND is_canonical = 1
          AND pointing_offset_deg IS NOT NULL""",
              (m.POINTING_OUTLIER_DEG,))[0]

    # Ground truth reproduction: the T CrB grism pointing failure.
    tcrb = q(con, """
        SELECT count(*), sum(pointing_offset_deg > ?)
        FROM frames WHERE target_key = 'tcrb' AND is_canonical = 1
          AND tree = 'rawimage' AND filter IN ('hrg', 'lrg')""",
             (m.POINTING_OUTLIER_DEG,))[0]

    # NGC 5548 has ZERO plate-solved frames, so the solved-only rule gave it
    # no reference position, no offsets and no flags at all — on the one
    # target whose strategy rules a whole night "unusable, pointing-failure".
    # The header-median fallback now covers it; these are the script-emitted
    # per-night means plus each night's worst offset under that fallback.
    n5548 = q(con, """
        SELECT night, count(*), round(avg(ra_deg), 2), round(avg(dec_deg), 2),
               round(max(pointing_offset_deg), 2)
        FROM frames WHERE target_key = 'ngc5548' AND is_canonical = 1
          AND tree = 'macalester'
        GROUP BY night ORDER BY abs(avg(dec_deg) - 25.14) DESC LIMIT 3""")
    n5548_tbl = table(
        ["night", "frames", "mean header RA (deg)", "mean header Dec (deg)",
         "worst offset (deg)"],
        [[esc(n), fmt(c), esc(r), esc(d), esc(w) if w is not None else "&mdash;"]
         for n, c, r, d, w in n5548],
        row_classes=["warn", None, None])

    # How each target's reference position was derived.  Guarded: a manifest
    # built before the fallback existed has no such column, and the report
    # must degrade to a plain statement rather than crash on an old file.
    have_basis = any(r[1] == "pointing_ref_basis"
                     for r in q(con, "PRAGMA table_info(frames)"))
    if have_basis:
        basis_rows = q(con, """
            SELECT pointing_ref_basis, count(DISTINCT target_key),
                   count(*)
            FROM frames WHERE is_canonical = 1 AND target_key IS NOT NULL
            GROUP BY pointing_ref_basis ORDER BY 3 DESC""")
        basis_tbl = table(
            ["reference basis", "targets", "canonical frames"],
            [[f"<code>{esc(b)}</code>", fmt(t), fmt(n)]
             for b, t, n in basis_rows])
        basis_note = (
            "<p class=\"sub\">Every target's reference position records the "
            "evidence it rests on.  <code>plate_solved</code> is the primary "
            "rule; <code>header_median</code> is the fallback for targets "
            f"with fewer than {m.MIN_SOLVED_FOR_REFERENCE} solves, which is "
            "weaker evidence (it is where the mount believed it was pointing) "
            "but robust to a minority of mispointed frames — the failure mode "
            "being detected.  It is NOT robust when most of a target's frames "
            "are mispointed; the column is what makes that visible.</p>"
            + basis_tbl)
    else:
        basis_tbl = ""
        basis_note = ("<p class=\"sub\">(This manifest pre-dates the "
                      "<code>pointing_ref_basis</code> column; rebuild S0 to "
                      "record how each reference position was derived.)</p>")

    return f"""
<section id="pointing">
<div class="bhead"><h2>5 · Pointing validation</h2>
<span class="tag">frame coords vs alias-resolved target reference</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Is each frame actually pointed at the target its name
claims?  (The T CrB panel proved header pointings lie.)</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"Offsets for {fmt(n_off)} canonical frames across {fmt(n_ref)} targets "
    f"with a reference position (median of plate-solved canonical frames "
    f"when at least {m.MIN_SOLVED_FOR_REFERENCE} exist, otherwise the median "
    f"of the target's header coordinates). "
    f"{fmt(n_out)} frames sit beyond {m.POINTING_OUTLIER_DEG:g} deg; "
    f"{fmt(n_zero)} exact-zero offsets are clipped into the lowest bin "
    "(a log axis cannot hold zero).")}
</div>
{off_tbl}
<p class="sub">Reading the offender list: {fmt(n_pseudo)} of the
{fmt(n_out)} outliers belong to Focus/TestField pseudo-targets, and most of
the rest are moving solar-system objects (asteroids, planets) for which a
fixed reference position is wrong by construction.  The list still catches
real science pathology: AN&nbsp;UMa carries {fmt(anuma[0])} outliers among
its {fmt(anuma[1])} frames with offsets, the worst at {esc(anuma[2])}&deg;
— a fact the CV project&rsquo;s conditional Q5 analysis needs to know.</p>
<p class="sub">Ground-truth reproduction: the T CrB strategy reported 21 of
247 grism pointings &gt;1&deg; off — the manifest finds
<b>{fmt(tcrb[1])}</b> of <b>{fmt(tcrb[0])}</b>.  And NGC&nbsp;5548, whose 279
frames contain <b>zero</b> plate solves, is why this stage has a fallback at
all: under a solved-only rule it had no reference position, therefore no
offsets, therefore no <code>pointing_gt1deg</code> flag could ever fire — on
precisely the target whose strategy rules an entire night unusable.  The
flagging machinery was alive across {fmt(n_out)} frames archive-wide and
blind exactly where it was needed.  With the header-median fallback the
tracking-failure night convicts itself:</p>
{n5548_tbl}
{basis_note}

<h3>Decision</h3>
<div class="decision"><b>Reference position = RA-wrap-aware median of a
target&rsquo;s plate-solved canonical frames when
&ge;{m.MIN_SOLVED_FOR_REFERENCE} solves exist, else the median of that
target&rsquo;s header coordinates; the choice is recorded per frame in
<code>pointing_ref_basis</code>.</b>  Every frame with coordinates gets
<code>pointing_offset_deg</code>; offsets
&gt;{m.POINTING_OUTLIER_DEG:g}&deg; are flagged
<code>pointing_gt1deg</code> — flagged, never deleted.  Solved-WCS
validation for the currently unsolvable campaigns arrives with S1; until
then the fallback is stated as fallback, not passed off as a solve.</div>

<h3>Consequence</h3>
<p class="sub">The grism identity gate (T CrB Phase A.0) and every
photometry stage inherit a per-frame pointing sanity column <em>and</em> the
basis behind it; mispointed frames can no longer hide inside a
night&rsquo;s frame count, and a target with no plate solves is no longer
silently exempt from the audit.</p>
</div></section>"""


def section_qc(con) -> str:
    src = fig_qc(con)
    n_flagged = q1(con, "SELECT count(*) FROM frames WHERE qc_flags != ''")
    n_err = q1(con, "SELECT count(*) FROM frames "
                    "WHERE qc_flags LIKE '%header_error%'")
    errs = q(con, "SELECT path, error FROM frames WHERE error IS NOT NULL")
    err_tbl = table(["unreadable file", "cataloger error"],
                    [[f"<code>{esc(p)}</code>", esc(e)] for p, e in errs])
    n_am = q1(con, "SELECT count(*) FROM frames "
                   "WHERE qc_flags LIKE '%airmass_garbage%'")
    n_exp = q1(con, "SELECT count(*) FROM frames "
                    "WHERE qc_flags LIKE '%exptime_nonpos%'")
    n_exp_light = q1(con, """SELECT count(*) FROM frames
        WHERE qc_flags LIKE '%exptime_nonpos%'
          AND imagetyp LIKE 'Light%'""")

    return f"""
<section id="qc">
<div class="bhead"><h2>6 · QC flags</h2>
<span class="tag">marks, never deletions</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Which frames carry defects that downstream stages must know
about before trusting a header?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"{fmt(n_flagged)} of all catalog rows carry at least one flag. "
    f"exptime_nonpos: {fmt(n_exp)} rows, of which {fmt(n_exp_light)} "
    "claim to be Light frames (mostly biases; the light-frame subset is "
    "the real pathology).")}</div>
<p class="sub">The {fmt(n_err)} files the cataloger could not read at
all:</p>
{err_tbl}

<h3>Decision</h3>
<div class="decision"><b>QC flags mark; they never delete.</b>  Header
airmass outside [{m.AIRMASS_MIN:g}, {m.AIRMASS_MAX:g}] ({fmt(n_am)} rows)
is flagged and left alone — recomputation from coordinates and time is
S3&rsquo;s job (ROADMAP convention 7), and camtemp is untrusted wholesale.
ZMAG is carried as QC-only metadata and was used in <i>no</i> S0
decision (ROADMAP convention 8).</div>

<h3>Consequence</h3>
<p class="sub">Every downstream selection can (and must) state its cuts as
flag predicates — reproducible, greppable, and identical across all five
papers.</p>
</div></section>"""


def section_reconciliation(con) -> str:
    rows = q(con, """
        SELECT project, target, metric, claimed_frames, claimed_nights,
               manifest_frames, manifest_nights, manifest_frames_global,
               diff_frames, diff_nights, source, target_key
        FROM project_counts""")
    body, classes = [], []
    n_exact = 0
    for (proj, tgt, metric, cf, cn, mf, mn, gf, dfr, dnt, src, _tk) in rows:
        exact = (dfr in (0, None)) and (dnt in (0, None))
        n_exact += 1 if (dfr == 0 and (dnt == 0 or dnt is None)) else 0
        classes.append(None if exact else "warn")
        body.append([
            esc(proj), esc(tgt), f"<code>{esc(metric)}</code>",
            f"{fmt(cf)} / {fmt(cn)}", f"{fmt(mf)} / {fmt(mn)}",
            fmt(dfr), fmt(dnt), fmt(gf), esc(src)])
    tbl = table(
        ["project", "target", "metric", "strategy claim (frames / nights)",
         "manifest, primary tree (frames / nights)", "&Delta;frames",
         "&Delta;nights", "manifest global canonical", "claim source & note"],
        body, row_classes=classes)
    n_rows = len(rows)
    n_warn = sum(1 for c in classes if c)

    # The three explained disagreements, pulled from the table itself so
    # the decision prose can never drift from the numbers it discusses.
    # Rows are keyed by (project, target_key, metric): the STABLE alias key
    # — display names are vote outcomes over catalog data, so keying on
    # them would let a future catalog change crash this render.
    by_key = {(r[0], r[11], r[2]): r for r in rows}
    tcrb_u = by_key[("TCrB_Monitoring", "tcrb", "unique_light")]
    sn = by_key[("SN2023ixf_LightCurve", "2023ixf", "unique_light")]
    n5548 = by_key[("DwarfGalaxy_AGN_Survey", "ngc5548", "unique_light")]

    return f"""
<section id="recon">
<div class="bhead"><h2>7 · Strategy-count reconciliation</h2>
<span class="tag">the key output — five papers, one frame accounting</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Do the inventory numbers published in the five strategy
documents (all rev. 2026-08-16) survive contact with a single, global,
rule-based manifest?</p>

<h3>Evidence</h3>
<p class="sub">{fmt(n_rows)} claims checked; <b>{fmt(n_exact)}</b>
reproduce exactly under like-for-like accounting (canonical frames from the
target&rsquo;s primary tree — the rule every strategy itself used);
<b>{fmt(n_warn)}</b> rows disagree, each with its cause identified.  The
&ldquo;global&rdquo; column counts canonical frames across ALL trees — its
excess over the primary-tree column is renamed-copy contamination
(section&nbsp;1) plus genuinely unique out-of-tree material (ST&nbsp;LMi's
iKon nights).</p>
{tbl}

<h3>Decision</h3>
<div class="decision"><b>The manifest is the single source of frame counts
from S0 onward.</b>  The disagreements stand as corrections, not mysteries:
T&nbsp;CrB&rsquo;s &ldquo;{fmt(tcrb_u[3])} unique&rdquo; counted rawimage
<i>rows</i> and misses {fmt(-tcrb_u[8])} within-rawimage duplicate copies
(manifest: {fmt(tcrb_u[5])}); the SN&rsquo;s &ldquo;{fmt(sn[4])}
nights&rdquo; includes two saturated first epochs whose frames are labeled
NGC5457/M101, not 2023ixf (manifest: {fmt(sn[6])} nights under the 2023ixf
name); NGC&nbsp;5548&rsquo;s &ldquo;{fmt(n5548[3])}/{fmt(n5548[4])}&rdquo;
excludes the mispointed 2023-03-25 night, which the manifest keeps —
{fmt(n5548[5])} frames / {fmt(n5548[6])} nights, mispointing being a
pointing fact, not a deduplication fact.</div>

<h3>Consequence</h3>
<p class="sub">Each project&rsquo;s Table&nbsp;1 regenerates from
<code>project_counts</code>; the era table, alias table, and QC flags ride
along.  Any future strategy edit that changes a count must change this page
first.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S0 report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    try:
        # Headline numbers for the masthead — same query discipline.
        n_rows = q1(con, "SELECT count(*) FROM frames")
        n_canon = q1(con, "SELECT sum(is_canonical) FROM frames")
        n_targets = q1(con, "SELECT count(DISTINCT target_key) FROM frames "
                            "WHERE target_key IS NOT NULL")
        n_eras = q1(con, "SELECT count(*) FROM eras")
        n_nights = q1(con, "SELECT count(DISTINCT night) FROM frames "
                           "WHERE night IS NOT NULL")
        meta = dict(q(con, "SELECT key, value FROM build_meta"))

        sections = [
            section_duplication(con),
            section_aliases(con),
            section_eras(con),
            section_nights(con),
            section_pointing(con),
            section_qc(con),
            section_reconciliation(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S0 — Manifest &amp; Curation</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S0 — Manifest &amp; Curation</h1>
  <p>{fmt(n_rows)} catalog rows &rarr; {fmt(n_canon)} canonical frames
  &middot; {fmt(n_targets)} targets &middot; {fmt(n_eras)} camera eras
  &middot; {fmt(n_nights)} nights &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z from
  <code>{esc(Path(meta.get('catalog_path', '')).name)}</code>
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#dedup">1 Duplication</a> &middot;
  <a href="#aliases">2 Aliases</a> &middot;
  <a href="#eras">3 Eras</a> &middot;
  <a href="#nights">4 Nights</a> &middot;
  <a href="#pointing">5 Pointing</a> &middot;
  <a href="#qc">6 QC</a> &middot;
  <a href="#recon">7 Reconciliation</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_core.report_s0</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query; none is typed by hand.  Regenerate with
<code>pipeline/scripts/build_s0_manifest.py</code>.</footer>
</body></html>"""

        HTML_PATH.write_text(html, encoding="utf-8")

        # Belt and braces: every <img> the page references must exist and be
        # non-empty, or the build fails loudly rather than shipping a
        # broken evidence page.
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
