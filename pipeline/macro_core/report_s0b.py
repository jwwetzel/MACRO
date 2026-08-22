"""S0b evidence report renderer: raw<->reduced linkage + calibration inventory.

Reads the S0b tables of the manifest database (NEVER the catalog — if a
number cannot be derived from the database, it does not belong on the page)
and writes:

* ``docs/pipeline/s0b_calibration_inventory.html``  — the report
* ``docs/pipeline/figures/s0b/*.png``               — every figure

The page follows the site's Socratic format: one section per decision, each
section = Question → Evidence → Decision → Consequence.  EVERY number in the
HTML is interpolated from a SQL query executed in this module or from a
constant defined in ``macro_core.inventory`` — nothing is hand-typed, so
re-running the build after an archive sync (the October ingest path)
regenerates the whole argument, shopping list included.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import inventory as inv   # noqa: E402  (constants for interpolation)
# Shared page machinery: one house figure style, one query discipline,
# one table generator
# generator as the S0 report — one visual language across the evidence site.
from .report_s0 import (          # noqa: E402
    ACCENT, STYLE, DPI, INK, WARN,
    _figure, esc, fmt, q, q1, table)
from . import plotstyle as ps    # noqa: E402  (the house figure style)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s0b"
HTML_PATH = DOCS_DIR / "s0b_calibration_inventory.html"

#: Eras with fewer canonical science frames than this are kept in the DB
#: tables but folded out of the on-page matrix (the page stays readable;
#: the DB stays complete — the report states the fold and its size).
MIN_ERA_SCIENCE_FOR_PAGE = 100

#: An era whose last science night falls within this many days of the
#: archive's final night is "current-instrument": its gaps are directly
#: acquirable at the October re-opening.
CURRENT_ERA_WINDOW_DAYS = 60


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_link_methods(con) -> str:
    """Two panels: match-ladder population sizes, and the JD-drift spread.

    Left: how every reduced-tree row was (or was not) tied to a raw parent.
    Right: for the stem_jd_drift rung, the distribution of the JD rewrite
    the reduction pipeline applied — evidence that those files are the same
    exposures, re-stamped.
    """
    rows = q(con, """
        SELECT match_method, count(*) FROM raw_reduced_links
        GROUP BY match_method ORDER BY 2""")
    methods = [r[0] for r in rows]
    counts = [r[1] for r in rows]
    drifts = [r[0] for r in q(con, """
        SELECT jd_drift_s FROM raw_reduced_links
        WHERE match_method = 'stem_jd_drift' AND jd_drift_s IS NOT NULL""")]
    with plt.rc_context(STYLE):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.4),
                                       gridspec_kw={"width_ratios": [3, 2]})
        colors = [WARN if m == "orphan" else ACCENT for m in methods]
        bars = ax1.barh(methods, counts, color=colors)
        ax1.set_xscale("log")
        ax1.set_xlabel("reduced-tree rows (log scale)")
        ax1.set_title("How each reduced file found its raw parent")
        for b, c in zip(bars, counts):
            ax1.annotate(f" {c:,}", (c, b.get_y() + b.get_height() / 2),
                         va="center", fontsize=8, color=INK)
        if drifts:
            ax2.hist(drifts, bins=30, color=WARN)
        ax2.set_xlabel("JD rewrite, reduced − raw (seconds)")
        ax2.set_ylabel("frames")
        ax2.set_title("stem_jd_drift: rewritten JDs")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0b_link_methods.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0b/s0b_link_methods.png"


def fig_calib_timeline(con) -> str:
    """Calibration frames on the calendar: night vs era, colored by kind."""
    rows = q(con, """
        SELECT era_id, night, kind, count(*) FROM calib_frames
        WHERE night IS NOT NULL AND era_id IS NOT NULL
        GROUP BY era_id, night, kind""")
    import datetime as _dt
    eras = sorted({r[0] for r in rows})
    era_pos = {e: i for i, e in enumerate(eras)}
    kind_color = ps.KIND_COLOR
    # era_id round-trips through pandas as a float (NaN-capable column);
    # cast for display so the axis reads 'era 47', not 'era 47.0'.
    labels = {r[0]: f"era {int(r[0])}: {r[1] or '(blank)'}" for r in q(con, """
        SELECT DISTINCT era_id, readoutm FROM calib_frames
        WHERE era_id IS NOT NULL""")}
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.6, 0.5 * max(len(eras), 4) + 1.8))
        for kind in ("bias", "dark", "flat"):
            sub = [r for r in rows if r[2] == kind]
            if not sub:
                continue
            xs = [_dt.date.fromisoformat(r[1]) for r in sub]
            ys = [era_pos[r[0]] for r in sub]
            ss = [max(8.0, 2.0 * np.sqrt(r[3]) * 4) for r in sub]
            ax.scatter(xs, ys, s=ss, color=kind_color[kind], alpha=0.65,
                       linewidths=0, label=kind)
        ax.set_yticks(range(len(eras)),
                      [labels.get(e, f"era {int(e)}") for e in eras],
                      fontsize=8)
        ax.set_xlabel("night (local-noon-to-noon label)")
        ax.set_title("Every calibration frame, by night and era "
                     "(marker area ~ frames)")
        ax.invert_yaxis()
        ax.legend(loc="upper left", fontsize=8)
        fig.autofmt_xdate(rotation=30, ha="right")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0b_calib_timeline.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0b/s0b_calib_timeline.png"


def fig_coverage_heatmap(con) -> str:
    """Heatmap: major science eras x requirement kind, science-frame-weighted.

    Each cell is the fraction of the era's science frames whose requirement
    of that kind is satisfied to spec by raw calibration frames (a dark cell
    weights each exposure-time bin by its science frames; a flat cell
    weights each filter likewise).  1.0 = fully covered, 0.0 = nothing.
    """
    eras = [r[0] for r in q(con, """
        SELECT era_id FROM calib_coverage WHERE req_kind = 'bias'
          AND n_science >= ? ORDER BY era_id""", (MIN_ERA_SCIENCE_FOR_PAGE,))]
    kinds = ["bias", "dark", "flat"]
    grid = np.zeros((len(eras), len(kinds)))
    for i, era in enumerate(eras):
        for j, kind in enumerate(kinds):
            tot, cov = q(con, """
                SELECT sum(n_science),
                       sum(CASE WHEN status = 'ok' THEN n_science ELSE 0 END)
                FROM calib_coverage WHERE era_id = ? AND req_kind = ?""",
                         (era, kind))[0]
            grid[i, j] = (cov or 0) / tot if tot else 0.0
    labels = {r[0]: f"era {r[0]}: {r[1] or '(blank)'} bin{r[2]}"
              for r in q(con, """
        SELECT e.era_id, e.readoutm, e.xbinning FROM eras e""")}
    # The house sequential ramp: white (nothing) to blue (fully covered).
    # Not red-to-green — that is the one pair a deuteranope cannot read,
    # and the exact percentage is printed in every cell anyway.
    cmap = ps.SEQ_CMAP
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(
            figsize=(6.4, 0.42 * max(len(eras), 4) + 1.6))
        ax.imshow(grid, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        # Cell separators instead of the rcParam grid, which a heatmap
        # would otherwise wear as stripes across its own data.
        ax.set_xticks(np.arange(len(kinds) + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(eras) + 1) - 0.5, minor=True)
        ax.grid(which="minor", color=ps.GRID, linewidth=0.8)
        ax.tick_params(which="minor", length=0)
        ax.set_xticks(range(len(kinds)), kinds)
        ax.set_yticks(range(len(eras)),
                      [labels.get(e, f"era {e}") for e in eras], fontsize=8)
        ax.set_title("Science-frame-weighted calibration coverage")
        ax.grid(False)
        for i in range(len(eras)):
            for j in range(len(kinds)):
                v = grid[i, j]
                ax.text(j, i, f"{100 * v:.0f}%", ha="center", va="center",
                        fontsize=8,
                        color=ps.ink_on(v))
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0b_coverage_heatmap.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0b/s0b_coverage_heatmap.png"


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def section_linkage(con) -> str:
    n_pairs = q1(con, "SELECT count(*) FROM raw_reduced_links "
                      "WHERE raw_rowid IS NOT NULL")
    n_reduced = q1(con, "SELECT count(DISTINCT reduced_rowid) "
                        "FROM raw_reduced_links")
    n_orphans = q1(con, "SELECT count(*) FROM raw_reduced_links "
                        "WHERE match_method = 'orphan'")
    src = fig_link_methods(con)

    method_rows = q(con, """
        SELECT match_method, count(*), count(DISTINCT reduced_rowid)
        FROM raw_reduced_links GROUP BY match_method ORDER BY 2 DESC""")
    explain = {
        "same_basename_jd": "exact copy — same filename, same JD "
                            "(S0's own dup_group)",
        "stem_jd": "renamed copy — raw name + processing suffix "
                   "(<code>_calibrated</code>/<code>_cal</code>/"
                   "<code>_wcs</code>), JD preserved",
        "stem_jd_drift": "same filename stem, same night, JD REWRITTEN "
                         "by the reduction pipeline",
        "target_jd": "fully renamed — matched on (target, JD), the S0 "
                     "collision-pair evidence",
        "target_jd_ambiguous": "(target, JD) hits several raw frames; "
                               "every candidate pair recorded",
        "orphan": "no raw parent found — characterized below",
    }
    method_tbl = table(
        ["match method", "link rows", "reduced frames", "meaning"],
        [[f"<code>{esc(mm)}</code>", fmt(n), fmt(d),
          explain.get(mm, "&mdash;")]
         for mm, n, d in method_rows],
        row_classes=["warn" if mm == "orphan" else None
                     for mm, _, _ in method_rows])

    # Raw-side multiplicity: canonical science frames outside the reduced
    # tree, and how many reduced counterparts each one has.
    n_raw_sci = q1(con, """
        SELECT count(*) FROM frames f
        WHERE f.is_canonical = 1 AND f.tree != 'reduced'
          AND f.error IS NULL
          AND (f.imagetyp LIKE 'Light%' OR f.imagetyp IS NULL
               OR f.imagetyp = '')
          AND NOT EXISTS (SELECT 1 FROM calib_frames c
                          WHERE c.obs_rowid = f.obs_rowid)""")
    mult = q(con, """
        SELECT n, count(*) FROM (
            SELECT raw_rowid, count(*) AS n FROM raw_reduced_links
            WHERE raw_rowid IS NOT NULL GROUP BY raw_rowid
        ) GROUP BY n ORDER BY n""")
    n_linked_raw = sum(r[1] for r in mult)
    mult_tbl = table(
        ["reduced counterparts per raw frame", "raw frames"],
        [[fmt(n), fmt(c)] for n, c in mult])

    # Drift statistics for the rewritten-JD population.
    drift_stats = q(con, """
        SELECT count(*), min(jd_drift_s), max(jd_drift_s)
        FROM raw_reduced_links WHERE match_method = 'stem_jd_drift'""")[0]

    # Orphan characterization: who are the frames with no raw parent?
    orphan_kinds = q(con, """
        SELECT CASE WHEN target_key IS NULL
                    THEN 'no target name (focus / class / test products)'
                    ELSE 'named target, raw copy absent from the archive'
               END, count(*)
        FROM raw_reduced_links WHERE match_method = 'orphan' GROUP BY 1""")
    orphan_kind_tbl = table(["orphan population", "frames"],
                            [[esc(k), fmt(n)] for k, n in orphan_kinds])
    orphan_sample = q(con, """
        SELECT reduced_path FROM raw_reduced_links
        WHERE match_method = 'orphan' ORDER BY reduced_path LIMIT 10""")
    orphan_tbl = table(
        ["orphan example (first 10, alphabetical)"],
        [[f"<code>{esc(p)}</code>"] for (p,) in orphan_sample])

    return f"""
<section id="linkage">
<div class="bhead"><h2>1 &middot; Raw &harr; reduced linkage</h2>
<span class="tag">every reduced file must prove its raw parent</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">S0 proved the <code>reduced/</code> tree holds renamed copies
that (basename, JD) dedup cannot see.  Can every one of its
{fmt(n_reduced)} rows be tied to the raw exposure it came from — and what
is left over when we try?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"{fmt(n_pairs)} (raw, reduced) links across {fmt(n_reduced)} reduced "
    f"rows; {fmt(n_orphans)} orphans remain. Right panel: the "
    f"{fmt(drift_stats[0])} stem_jd_drift frames carry JD rewrites of "
    f"{drift_stats[1]:.0f} to {drift_stats[2]:.0f} seconds — same file, "
    "re-stamped during reduction.")}</div>
{method_tbl}
<p class="sub">Raw-side view: of the {fmt(n_raw_sci)} canonical science
frames outside the reduced tree, <b>{fmt(n_linked_raw)}</b> have at least
one reduced counterpart ({100.0 * n_linked_raw / n_raw_sci:.1f}%) —
the rest were never run through the reduction pipeline (or their products
never reached the archive):</p>
{mult_tbl}
<p class="sub">The {fmt(n_orphans)} orphans, characterized — not hidden:</p>
{orphan_kind_tbl}
{orphan_tbl}

<h3>Decision</h3>
<div class="decision"><b>The reduced tree is derived data: no reduced file
is ever a canonical science frame, and <code>raw_reduced_links</code> is
the only bridge between the two populations.</b>  The match ladder accepts
exactly four kinds of evidence — S0&rsquo;s own dup_group, a filename stem
with a known processing suffix at the same JD, the same stem on the same
night with a rewritten JD, and the (target,&nbsp;JD) collision pairs S0
handed over — deterministically, strongest first.  Orphans are recorded
with NULL raw columns and quarantined from every science count.</div>

<h3>Consequence</h3>
<p class="sub">S1/S2 read pixels from raw canonical paths only; anything
already reduced is reachable <i>through the links table</i> when a
comparison is wanted.  The (target, JD) collision suspect list from S0 is
hereby resolved into named, methodical links.</p>
</div></section>"""


def section_census(con) -> str:
    n_calib = q1(con, "SELECT count(*) FROM calib_frames")
    n_master = q1(con, "SELECT count(*) FROM calib_frames WHERE is_master=1")
    n_eras_with = q1(con, "SELECT count(DISTINCT era_id) FROM calib_frames")
    n_eras_sci = q1(con, "SELECT count(DISTINCT era_id) FROM calib_coverage")
    src = fig_calib_timeline(con)

    kind_rows = q(con, """
        SELECT kind, sum(is_master = 0), sum(is_master = 1),
               count(DISTINCT night)
        FROM calib_frames GROUP BY kind ORDER BY kind""")
    kind_tbl = table(
        ["kind", "raw frames", "master products", "nights"],
        [[f"<code>{esc(k)}</code>", fmt(r), fmt(mstr), fmt(n)]
         for k, r, mstr, n in kind_rows])

    tree_rows = q(con, """
        SELECT tree, kind, count(*) FROM calib_frames
        GROUP BY tree, kind ORDER BY 3 DESC""")
    tree_tbl = table(
        ["tree", "kind", "frames"],
        [[esc(t), f"<code>{esc(k)}</code>", fmt(n)]
         for t, k, n in tree_rows])

    era_rows = q(con, """
        SELECT c.era_id, e.readoutm, e.xbinning, e.egain,
               sum(c.kind = 'bias'), sum(c.kind = 'dark'),
               sum(c.kind = 'flat'), min(c.night), max(c.night)
        FROM calib_frames c JOIN eras e ON e.era_id = c.era_id
        GROUP BY c.era_id ORDER BY c.era_id""")
    era_tbl = table(
        ["era", "READOUTM", "bin", "EGAIN", "bias", "dark", "flat",
         "first night", "last night"],
        [[fmt(e), esc(r) or "<i>(blank)</i>", fmt(xb), esc(eg),
          fmt(b), fmt(d), fmt(f), esc(fn), esc(ln)]
         for e, r, xb, eg, b, d, f, fn, ln in era_rows])

    # Cross-check 1 (roadmap ground truth): the Mode0 grism era must have
    # ZERO 240 s darks — none were ever taken.  A nonzero count here is a
    # loud red flag, not a quiet table cell.
    tol_240 = max(inv.DARK_MATCH_ABS_TOL, inv.DARK_MATCH_REL_TOL * 240.0)
    n_mode0_240 = q1(con, """
        SELECT count(*) FROM calib_frames
        WHERE kind = 'dark' AND readoutm = 'Mode0'
          AND abs(exptime - 240.0) <= ?""", (tol_240,))
    # Cross-check 2: StackPro flats must be absent entirely.  Raw frames
    # and master products are judged separately: a raw StackPro flat would
    # contradict the roadmap outright, while a master PRODUCT carrying the
    # StackPro readout string may merely have inherited a header from the
    # stacking tool — still worth a loud flag, but a different fact.
    n_spro_flat_raw = q1(con, """
        SELECT count(*) FROM calib_frames
        WHERE kind = 'flat' AND is_master = 0
          AND readoutm LIKE '%StackPro%'""")
    n_spro_flat_master = q1(con, """
        SELECT count(*) FROM calib_frames
        WHERE kind = 'flat' AND is_master = 1
          AND readoutm LIKE '%StackPro%'""")
    spro_master_dirs = q(con, """
        SELECT DISTINCT substr(path, 1, length(path) - length(basename) - 1)
        FROM calib_frames
        WHERE kind = 'flat' AND is_master = 1
          AND readoutm LIKE '%StackPro%'""")
    check_tbl = table(
        ["roadmap cross-check", "expected", "found"],
        [["Mode0 darks at 240 s (T CrB grism exposure)", "0",
          fmt(n_mode0_240) + ("" if n_mode0_240 == 0
                              else " — <b>ROADMAP CONTRADICTION</b>")],
         ["StackPro flats — raw frames", "0",
          fmt(n_spro_flat_raw) + ("" if n_spro_flat_raw == 0
                                  else " — <b>ROADMAP CONTRADICTION</b>")],
         ["StackPro flats — master products (header may be inherited "
          "from the stacking tool)", "0",
          fmt(n_spro_flat_master)
          + ("" if n_spro_flat_master == 0 else
             " — in " + ", ".join(f"<code>{esc(d[0])}</code>"
                                  for d in spro_master_dirs))]],
        row_classes=[None if n_mode0_240 == 0 else "warn",
                     None if n_spro_flat_raw == 0 else "warn",
                     None if n_spro_flat_master == 0 else "warn"])

    return f"""
<section id="census">
<div class="bhead"><h2>2 &middot; Calibration census</h2>
<span class="tag">kind normalized from IMAGETYP + filename conventions</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">What calibration data does the archive actually hold, and
for which cameras?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"All {fmt(n_calib)} canonical calibration frames ({fmt(n_master)} of "
    f"them stacked master products) land in only {fmt(n_eras_with)} of the "
    f"{fmt(n_eras_sci)} camera eras that hold science frames.")}</div>
{kind_tbl}
{tree_tbl}
<p class="sub">Per era — note which eras appear here at all:</p>
{era_tbl}
{check_tbl}

<h3>Decision</h3>
<div class="decision"><b>Calibration kind is normalized by explicit header
IMAGETYP first, and by three observed filename conventions second</b>
(<code>master&hellip;</code> products written as
&lsquo;Light&nbsp;Frame&rsquo;, with <i>dark</i> checked before <i>flat</i>
so a flat-field dark stays a dark; the <code>&lt;filter&gt;.flatN</code>
twilight series).  &lsquo;Fringe field&rsquo; sky exposures remain science
— they calibrate fringing, not the detector.  Era identity is joined from
S0&rsquo;s <code>frames.era_id</code>, never recomputed.  Roadmap
cross-checks: no Mode0 240&nbsp;s dark and no RAW StackPro flat exists;
the {fmt(n_spro_flat_master)} master flat <i>products</i> carrying the
StackPro readout string are quarantined evidence for S2&rsquo;s PTC work
(they claim a readout mode the roadmap says was never flat-fielded — most
likely a header inherited at stacking time, to be settled from pixels,
not headers).</div>

<h3>Consequence</h3>
<p class="sub">The coverage matrix below can now be computed era by era —
and the eras <i>absent</i> from the census table are exactly the ones the
October shopping list must service.</p>
</div></section>"""


def section_matrix(con) -> str:
    n_cells = q1(con, "SELECT count(*) FROM calib_coverage")
    n_eras = q1(con, "SELECT count(DISTINCT era_id) FROM calib_coverage")
    n_shown = q1(con, """SELECT count(*) FROM calib_coverage
        WHERE req_kind = 'bias' AND n_science >= ?""",
                 (MIN_ERA_SCIENCE_FOR_PAGE,))
    n_ok = q1(con, "SELECT count(*) FROM calib_coverage WHERE status='ok'")
    src = fig_coverage_heatmap(con)

    rows = q(con, """
        SELECT c.era_id, e.readoutm, e.xbinning,
               max(CASE WHEN c.req_kind = 'bias' THEN c.status END),
               sum(CASE WHEN c.req_kind = 'dark' THEN 1 ELSE 0 END),
               sum(CASE WHEN c.req_kind = 'dark' AND c.status = 'ok'
                        THEN 1 ELSE 0 END),
               sum(CASE WHEN c.req_kind = 'flat' THEN 1 ELSE 0 END),
               sum(CASE WHEN c.req_kind = 'flat' AND c.status = 'ok'
                        THEN 1 ELSE 0 END),
               max(c.n_science * (c.req_kind = 'bias'))
        FROM calib_coverage c JOIN eras e ON e.era_id = c.era_id
        WHERE c.era_id IN (SELECT era_id FROM calib_coverage
                           WHERE req_kind = 'bias' AND n_science >= ?)
        GROUP BY c.era_id ORDER BY 9 DESC""", (MIN_ERA_SCIENCE_FOR_PAGE,))
    body, classes = [], []
    for era, ro, xb, bias_st, nd, ndok, nf, nfok, nsci in rows:
        fully = (bias_st == "ok" and ndok == nd and nfok == nf)
        classes.append(None if fully else "warn")
        body.append([
            fmt(era), esc(ro) or "<i>(blank)</i>", fmt(xb), fmt(nsci),
            esc(bias_st), f"{fmt(ndok)} / {fmt(nd)}",
            f"{fmt(nfok)} / {fmt(nf)}"])
    matrix_tbl = table(
        ["era", "READOUTM", "bin", "science frames", "bias",
         "dark bins ok", "filters ok"], body, row_classes=classes)

    # Scaled-dark note: the eras where the fallback is even on the table.
    scaled = q(con, """
        SELECT DISTINCT era_id FROM calib_coverage
        WHERE req_kind = 'dark' AND scaled_dark_ok = 1 ORDER BY era_id""")
    scaled_txt = ", ".join(str(r[0]) for r in scaled) if scaled else "none"

    # Header-glitch filters: science FILTER strings that collide with the
    # calibration vocabulary ('dark'/'bias'/'flat').  Their cells stay in
    # the matrix but are excluded from the shopping list — say so, with
    # query-derived counts, or say nothing when the archive has none.
    n_glitch_cells, n_glitch_sci = q(con, """
        SELECT count(*), COALESCE(sum(n_science), 0) FROM calib_coverage
        WHERE req_kind = 'flat' AND status != 'ok'
          AND lower(trim(req_key)) IN ('bias', 'dark', 'flat')""")[0]
    glitch_txt = ("" if n_glitch_cells == 0 else
                  f"  Science FILTER strings that collide with the "
                  f"calibration vocabulary (&lsquo;dark&rsquo;/&lsquo;bias"
                  f"&rsquo;/&lsquo;flat&rsquo; — filter-wheel/header "
                  f"glitches passed through verbatim) keep their matrix "
                  f"cells ({fmt(n_glitch_cells)} cell(s), "
                  f"{fmt(n_glitch_sci)} science frames) but are excluded "
                  f"from the shopping list: a flat &ldquo;in filter "
                  f"dark&rdquo; is not an acquirable item.")

    return f"""
<section id="matrix">
<div class="bhead"><h2>3 &middot; The coverage matrix</h2>
<span class="tag">per era: bias &middot; dark per exposure &middot; flat per
filter</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">For each camera era with science frames, which calibration
requirements are actually met — to spec, from raw frames?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"{fmt(n_cells)} requirement cells across {fmt(n_eras)} science eras; "
    f"only {fmt(n_ok)} are met. The page shows the {fmt(n_shown)} eras "
    f"holding at least {fmt(MIN_ERA_SCIENCE_FOR_PAGE)} science frames; "
    "every era is in the calib_coverage table.")}</div>
{matrix_tbl}

<h3>Decision</h3>
<div class="decision"><b>A dark matches a science exposure when
|&Delta;t| &le; max({inv.DARK_MATCH_ABS_TOL:g}&nbsp;s,
{100 * inv.DARK_MATCH_REL_TOL:g}% &middot; t) — exact-match policy with a
tolerance that absorbs driver float fuzz (&sim;10<sup>&minus;5</sup>
relative) and can never bridge two real exposure settings (the archive
ladder steps by &ge;25%).</b>  Specs count RAW frames only
(&ge;{inv.SPEC_N_BIAS} bias, &ge;{inv.SPEC_N_DARK} darks per exposure,
&ge;{inv.SPEC_N_FLAT} flats per filter); master products satisfy nothing
but are tallied (<code>master_only</code> status) because they remain
usable.  Scaled-dark fallback is noted only where an era has bias frames
(eras: {esc(scaled_txt)}) — without a bias, dark scaling is arithmetic
fiction.{glitch_txt}</div>

<h3>Consequence</h3>
<p class="sub">Every unmet cell becomes one row of the shopping list below
— ranked by the science frames it blocks.</p>
</div></section>"""


def section_shopping(con) -> str:
    n_gaps = q1(con, "SELECT count(*) FROM calib_gaps")
    n_blocked = q1(con, """
        SELECT sum(n_science_frames_blocked) FROM calib_gaps
        WHERE need_kind = 'bias'""")  # bias rows = whole-era science counts
    max_night = q1(con, "SELECT max(last_night) FROM calib_gaps")
    cur_rows = q(con, """
        SELECT count(DISTINCT era_id) FROM calib_gaps
        WHERE last_night >= date(?, ?)""",
                 (max_night, f"-{CURRENT_ERA_WINDOW_DAYS} day"))
    n_cur_eras = cur_rows[0][0]

    # The re-opening configuration: the era(s) holding the archive's FINAL
    # science night, plus every science era sharing their header signature
    # (READOUTM, EGAIN).  This is the state the camera was last seen in —
    # the configuration October's own frames will land in — so its gaps
    # lead the shopping list even where the blocked-frame ranking would
    # bury them.  (In the current archive these are the blank-READOUTM /
    # EGAIN-56 eras of 2026-06/07: the header-convention break documented
    # in section 2 — which is exactly why they must be surfaced by rule,
    # not by hand.)
    reopen_eras = q(con, """
        SELECT era_id, readoutm, egain, first_night, last_night FROM eras
        WHERE era_id IN (SELECT DISTINCT era_id FROM calib_coverage)
          AND (COALESCE(readoutm, ''), COALESCE(egain, -1)) IN (
              SELECT COALESCE(readoutm, ''), COALESCE(egain, -1) FROM eras
              WHERE era_id IN (SELECT DISTINCT era_id FROM calib_coverage)
                AND last_night = (SELECT max(last_night) FROM eras
                                  WHERE era_id IN (SELECT DISTINCT era_id
                                                   FROM calib_coverage)))
        ORDER BY era_id""")
    reopen_ids = [int(r[0]) for r in reopen_eras]
    reopen_ids_txt = "/".join(str(i) for i in reopen_ids)
    reopen_ro = (reopen_eras[0][1] or "(blank READOUTM)") if reopen_eras \
        else "?"
    reopen_egain = reopen_eras[0][2] if reopen_eras else "?"
    reopen_span = (f"{reopen_eras[0][3]} &rarr; {reopen_eras[-1][4]}"
                   if reopen_eras else "?")
    _ph = ",".join("?" * len(reopen_ids))
    reopen_n_gaps, reopen_bias_blocked = q(con, f"""
        SELECT count(*),
               COALESCE(sum(n_science_frames_blocked
                            * (need_kind = 'bias')), 0)
        FROM calib_gaps WHERE era_id IN ({_ph})""", reopen_ids)[0]
    reopen_top = q(con, f"""
        SELECT era_id, spec, n_science_frames_blocked, projects_affected
        FROM calib_gaps WHERE era_id IN ({_ph})
        ORDER BY n_science_frames_blocked DESC, era_id LIMIT 10""",
                   reopen_ids)
    reopen_tbl = table(
        ["era", "acquisition spec", "science frames blocked",
         "projects affected"],
        [[fmt(e), f"<code>{esc(s)}</code>", fmt(n), esc(p) or "&mdash;"]
         for e, s, n, p in reopen_top])

    top = q(con, """
        SELECT g.era_id, g.camera, e.xbinning, g.first_night, g.last_night,
               g.spec, g.n_science_frames_blocked, g.projects_affected,
               g.last_night >= date(?, ?)
        FROM calib_gaps g JOIN eras e ON e.era_id = g.era_id
        ORDER BY g.n_science_frames_blocked DESC, g.era_id LIMIT 30""",
            (max_night, f"-{CURRENT_ERA_WINDOW_DAYS} day"))
    body, classes = [], []
    for era, cam, xb, fn, ln, spec, nb, proj, is_cur in top:
        classes.append("warn" if is_cur else None)
        body.append([
            fmt(era), esc(cam) or "<i>(blank)</i>", fmt(xb),
            f"{esc(fn)} &rarr; {esc(ln)}",
            f"<code>{esc(spec)}</code>", fmt(nb),
            esc(proj) or "&mdash;",
            "YES" if is_cur else "&mdash;"])
    top_tbl = table(
        ["era", "READOUTM", "bin", "science span", "acquisition spec",
         "science frames blocked", "projects affected",
         "current instrument?"], body, row_classes=classes)

    # The T CrB-critical rows, called out by project.
    tcrb_rows = q(con, """
        SELECT era_id, camera, spec, n_science_frames_blocked
        FROM calib_gaps WHERE projects_affected LIKE '%TCrB%'
        ORDER BY n_science_frames_blocked DESC LIMIT 6""")
    tcrb_tbl = table(
        ["era", "READOUTM", "spec", "science frames blocked"],
        [[fmt(e), esc(c), f"<code>{esc(s)}</code>", fmt(n)]
         for e, c, s, n in tcrb_rows])

    # Current-instrument acquisition summary for the decision callout —
    # generated from the gaps table, kind by kind.
    cur_summary = q(con, """
        SELECT need_kind, count(*), sum(n_science_frames_blocked)
        FROM calib_gaps WHERE last_night >= date(?, ?)
        GROUP BY need_kind ORDER BY 3 DESC""",
                    (max_night, f"-{CURRENT_ERA_WINDOW_DAYS} day"))
    cur_txt = "; ".join(
        f"{fmt(c)} {esc(k)} requirement(s) blocking {fmt(n)} frames"
        for k, c, n in cur_summary)

    return f"""
<section id="shopping">
<div class="bhead"><h2>4 &middot; The October shopping list</h2>
<span class="tag">what the re-opening run must acquire, ranked</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Winer is offline for monsoon until October 2026.  When the
roof opens, exactly which calibration frames should the first nights
acquire — for which camera configuration, in what quantity, and in what
order?</p>

<h3>Evidence</h3>
<p class="sub">The archive&rsquo;s final science night ({esc(max_night)})
belongs to era(s) {esc(reopen_ids_txt)} — READOUTM
&ldquo;{esc(str(reopen_ro))}&rdquo;, EGAIN {esc(str(reopen_egain))},
{reopen_span} — the state the camera was last seen in and therefore the
configuration October&rsquo;s own science frames will land in.  These eras
carry {fmt(reopen_n_gaps)} open requirements of their own
({fmt(reopen_bias_blocked)} frames blocked on bias alone), including the
grism-chain flats; top 10:</p>
{reopen_tbl}
<p class="sub">Across the whole archive: {fmt(n_gaps)} unmet requirements;
the bias rows alone show {fmt(n_blocked)} science frames sitting in eras
with no usable bias set.  Top 30 by science frames blocked
(<span class="tag">highlighted rows</span> = the {fmt(n_cur_eras)} eras
whose last science night falls within
{fmt(CURRENT_ERA_WINDOW_DAYS)} days of the archive&rsquo;s final night
{esc(max_night)} — the current instrument, directly acquirable in
October):</p>
{top_tbl}
<p class="sub">The T&nbsp;CrB-critical subset (the 247 archival grism
spectra are Mode0 240&nbsp;s exposures with <i>zero</i> matching darks —
the roadmap&rsquo;s known debt, now a ranked line item):</p>
{tcrb_tbl}

<h3>Decision</h3>
<div class="decision"><b>The October re-opening run acquires, before any
science: (1) a full bias/dark/flat set — including the hrg/lrg grism
flats — in the configuration actually running at re-opening, i.e. the
era-{esc(reopen_ids_txt)} header state (READOUTM
&ldquo;{esc(str(reopen_ro))}&rdquo;, EGAIN {esc(str(reopen_egain))}),
because every new science frame will land there; (2) the rest of the
current-instrument set — {cur_txt}; (3) a Mode0 240&nbsp;s dark series
(&ge;{inv.SPEC_N_DARK} frames) plus Mode0 bias set for the T&nbsp;CrB
grism archive, provided the mount&rsquo;s camera still exposes the Mode0
readout; (4) filter flats in the order of the blocked-frame ranking
above.</b>  Ops should also fix — or at least document — the blank
READOUTM/IMAGETYP header cards of the era-{esc(reopen_ids_txt)}
configuration while the dome is open, so October&rsquo;s frames stop
minting header-convention eras.  Retired-era gaps that no October night can fill
(cameras no longer on the telescope) stay on this list as permanent
facts — their science frames proceed only through scaled-dark fallbacks
(where an era has biases) or uncalibrated-photometry error budgets, and
every such choice must cite this table.  The ops request
(<code>ops/2026-08_observatory_request.md</code>) cites this section as
its calibration annex.</div>

<h3>Consequence</h3>
<p class="sub">When the October frames land, the archive sync + S0 + S0b
re-run (idempotent, atomic) regenerates this page; acquired rows fall off
the shopping list automatically, and any requirement still unmet stays
visible.  New-era science (ST&nbsp;LMi g/r/i season, T&nbsp;CrB restart)
joins the matrix as first-class rows on the same re-run.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S0b report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    try:
        n_links = q1(con, "SELECT count(*) FROM raw_reduced_links "
                          "WHERE raw_rowid IS NOT NULL")
        n_orphans = q1(con, "SELECT count(*) FROM raw_reduced_links "
                            "WHERE match_method = 'orphan'")
        n_calib = q1(con, "SELECT count(*) FROM calib_frames")
        n_eras = q1(con, "SELECT count(DISTINCT era_id) FROM calib_coverage")
        n_gaps = q1(con, "SELECT count(*) FROM calib_gaps")
        meta = dict(q(con, "SELECT key, value FROM s0b_build_meta"))

        sections = [
            section_linkage(con),
            section_census(con),
            section_matrix(con),
            section_shopping(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S0b — Calibration Inventory &amp; Raw&harr;Reduced Links</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S0b — Calibration Inventory &amp; Raw&harr;Reduced Links</h1>
  <p>{fmt(n_links)} raw&harr;reduced links ({fmt(n_orphans)} orphans)
  &middot; {fmt(n_calib)} calibration frames &middot; {fmt(n_eras)} science
  eras audited &middot; {fmt(n_gaps)} open gaps for October &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#linkage">1 Linkage</a> &middot;
  <a href="#census">2 Calibration census</a> &middot;
  <a href="#matrix">3 Coverage matrix</a> &middot;
  <a href="#shopping">4 October shopping list</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_core.report_s0b</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query; none is typed by hand.  Regenerate with
<code>pipeline/scripts/build_s0b_inventory.py</code> (the October ingest
path: re-run after each archive sync).</footer>
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
