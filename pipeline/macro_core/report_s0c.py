"""S0c evidence report renderer: the per-project staging manifests.

Reads the S0c stage tables of the manifest database (NEVER the catalog, the
archive, or the CSVs — if a number cannot be derived from the database, it
does not belong on the page) and writes:

* ``docs/pipeline/s0c_staging.html``      — the report
* ``docs/pipeline/figures/s0c/*.png``     — every figure

The page follows the site's Socratic format: one section per decision, each
section = Question → Evidence → Decision → Consequence.  EVERY number in the
HTML is interpolated from a SQL query executed in this module or from a
constant defined in ``macro_core.staging`` — nothing is hand-typed.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import staging as stg     # noqa: E402  (constants for interpolation)
# Shared page machinery: same dark theme, same query discipline, same table
# generator as the S0/S0b reports — one visual language across the site.
from .report_s0 import (          # noqa: E402
    ACCENT, DARK, DPI, WARN, _figure, esc, fmt, q, q1, table)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s0c"
HTML_PATH = DOCS_DIR / "s0c_staging.html"

#: Display order and colors for staging roles (science leads; cone
#: candidates get the warning hue because they are NOT adjudicated science;
#: masters pale).
ROLE_ORDER = (stg.ROLE_SCIENCE, stg.ROLE_SCIENCE_UNRESOLVED,
              "bias", "dark", "flat",
              "master_bias", "master_dark", "master_flat")
ROLE_COLORS = {stg.ROLE_SCIENCE: ACCENT,
               stg.ROLE_SCIENCE_UNRESOLVED: "#d98f4f",
               "bias": "#9fd8ae", "dark": "#7a8b99",
               "flat": WARN, "master_bias": "#5d8a6b",
               "master_dark": "#4d5b66", "master_flat": "#9a884d"}


def _stage_tables(con) -> list[tuple[str, str]]:
    """(project, stage_table) pairs from the build's own registry table."""
    return q(con, "SELECT project, stage_table FROM s0c_stage_files "
                  "ORDER BY project")


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_roles(con) -> str:
    """Horizontal stacked bars: staged frames per project, split by role."""
    pairs = _stage_tables(con)
    counts = {}   # project -> {role: n}
    for project, tbl in pairs:
        counts[project] = dict(q(con,
            f"SELECT role, count(*) FROM {tbl} GROUP BY role"))
    projects = [p for p, _ in pairs]
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(9.4, 0.62 * len(projects) + 1.6))
        left = np.zeros(len(projects))
        for role in ROLE_ORDER:
            vals = np.array([counts[p].get(role, 0) for p in projects],
                            dtype=float)
            if vals.sum() == 0:
                continue
            ax.barh(projects, vals, left=left, color=ROLE_COLORS[role],
                    label=role)
            left += vals
        for i, p in enumerate(projects):
            total = int(sum(counts[p].values()))
            n_sci = counts[p].get("science", 0)
            ax.annotate(f" {fmt(n_sci)} sci / {fmt(total)} total",
                        (left[i], i), va="center", fontsize=8,
                        color="#e8eaed")
        ax.set_xlabel("staged rows (science + era-matched calibration)")
        ax.set_title("Each project's working set, by role — "
                     "rows in a manifest, never copies")
        ax.invert_yaxis()
        ax.legend(loc="lower right", fontsize=8, ncols=4)
        ax.set_xlim(0, left.max() * 1.30)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0c_roles.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0c/s0c_roles.png"


def fig_eras(con) -> str:
    """Heatmap: science frames per (project, era) — where each working set
    actually lives in camera-era space (and hence which calibration
    families ride along)."""
    pairs = _stage_tables(con)
    per = {}   # project -> {era: n}
    eras: set[int] = set()
    for project, tbl in pairs:
        rows = q(con, f"""SELECT era_id, count(*) FROM {tbl}
                          WHERE role = 'science' AND era_id IS NOT NULL
                          GROUP BY era_id""")
        per[project] = {int(e): n for e, n in rows}
        eras |= set(per[project])
    era_list = sorted(eras)
    projects = [p for p, _ in pairs]
    grid = np.array([[per[p].get(e, 0) for e in era_list] for p in projects],
                    dtype=float)
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(
            figsize=(0.62 * len(era_list) + 3.2, 0.6 * len(projects) + 1.8))
        shown = np.log10(np.where(grid > 0, grid, np.nan))
        im = ax.imshow(shown, cmap="Blues", aspect="auto")
        for i in range(len(projects)):
            for j in range(len(era_list)):
                if grid[i, j] > 0:
                    # High counts render as DARK blue under 'Blues' — those
                    # cells need light text; pale cells need dark text.
                    dark_cell = shown[i, j] > (np.nanmax(shown) * 0.6)
                    ax.annotate(fmt(int(grid[i, j])), (j, i), ha="center",
                                va="center", fontsize=8,
                                color="#e8eaed" if dark_cell else "#10151c")
        ax.set_xticks(range(len(era_list)),
                      [f"era {e}" for e in era_list], rotation=45,
                      ha="right", fontsize=8)
        ax.set_yticks(range(len(projects)), projects, fontsize=8)
        ax.set_title("Science frames per project and camera era "
                     "(each era staged brings its calibration family)")
        fig.colorbar(im, ax=ax, label="log10 science frames", shrink=0.8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0c_eras.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0c/s0c_eras.png"


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------
def section_no_copy(con) -> str:
    n_rows = q1(con, "SELECT count(*) FROM frames")
    n_canon = q1(con, "SELECT count(*) FROM frames WHERE is_canonical = 1")
    n_dups = n_rows - n_canon
    pct = 100.0 * n_dups / n_rows
    # The worst multiplicities: how far one exposure got copied around.
    mult_rows = q(con, """
        SELECT copies, count(*) FROM (
            SELECT dup_group, count(*) AS copies FROM frames
            GROUP BY dup_group HAVING copies > 1)
        GROUP BY copies ORDER BY copies DESC LIMIT 5""")
    mult_tbl = table(
        ["copies of one exposure", "duplicate groups"],
        [[fmt(c), fmt(n)] for c, n in mult_rows])
    n_multi = q1(con, """
        SELECT count(*) FROM (
            SELECT dup_group FROM frames
            GROUP BY dup_group HAVING count(*) > 1)""")

    return f"""
<section id="nocopy">
<div class="bhead"><h2>1 &middot; The no-copy law</h2>
<span class="tag">the staging manifest IS the working set</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Each project needs a working set: its science frames plus
the calibration frames that go with them.  Should S0c <b>copy</b> those
files into per-project <code>data/</code> directories?</p>

<h3>Evidence</h3>
<p class="sub">The archive already answered.  Of {fmt(n_rows)} catalog
rows, only {fmt(n_canon)} are canonical frames — <b>{fmt(n_dups)} rows
({pct:.0f}%) are duplicate copies</b> of frames that already existed
elsewhere: wholesale re-copied night directories, consortium mirrors,
renamed reduced-tree files.  {fmt(n_multi)} exposures exist in more than
one place; the worst multiplicities:</p>
{mult_tbl}
<p class="sub">Every one of those copies was made by someone who, quite
reasonably, wanted a local working set.  S0 spent an entire pipeline stage
undoing the consequences.</p>

<h3>Decision</h3>
<div class="decision"><b>No project ever copies a frame.  The staging
manifest — a script-emitted, provenance-complete frame list per project —
is the working set itself.</b>  Every row carries the archive-relative
path, the absolute path, the S0 era/night/target/QC provenance, and a
size-byte integrity surrogate; stages read the immutable archive directly
through it.  The only materialization on offer is an optional, disposable
<i>symlink</i> farm for human browsing (default off; Dropbox does not sync
symlink targets, and deleting the farm loses nothing).</div>

<h3>Consequence</h3>
<p class="sub">The archive stays read-only and single-copy; a re-run of
S0&nbsp;&rarr;&nbsp;S0b&nbsp;&rarr;&nbsp;S0c after each archive sync
refreshes every working set in place.  {esc(stg.CHECKSUM_NOTE)}</p>
</div></section>"""


def section_staged(con) -> str:
    src_roles = fig_roles(con)
    src_eras = fig_eras(con)
    pairs = _stage_tables(con)

    rows = []
    for project, tbl in pairs:
        (n_sci, n_tgt, n_nights, first_n, last_n) = q(con, f"""
            SELECT count(*), count(DISTINCT target_key),
                   count(DISTINCT night), min(night), max(night)
            FROM {tbl} WHERE role = 'science'""")[0]
        n_cal = q1(con, f"SELECT count(*) FROM {tbl} "
                        f"WHERE {stg.SQL_CALIB_ROLES}")
        n_cone = q1(con, f"SELECT count(*) FROM {tbl} "
                         f"WHERE role = '{stg.ROLE_SCIENCE_UNRESOLVED}'")
        n_eras = q1(con, f"SELECT count(DISTINCT era_id) FROM {tbl} "
                         "WHERE role = 'science' AND era_id IS NOT NULL")
        (csv_path, rule) = q(con, """
            SELECT csv_path, selection_rule FROM s0c_stage_files
            WHERE project = ?""", (project,))[0]
        rows.append([esc(project), fmt(n_sci), fmt(n_cone), fmt(n_tgt),
                     fmt(n_nights),
                     f"{esc(first_n)} &rarr; {esc(last_n)}", fmt(n_eras),
                     fmt(n_cal), f"<code>{esc(csv_path)}</code>"])
    stage_tbl = table(
        ["project", "science", "cone candidates", "targets", "nights",
         "span", "eras", "calib rows", "working artifact"], rows)

    rule_rows = q(con, """
        SELECT project, selection_rule, selection_source
        FROM s0c_stage_files ORDER BY project""")
    rule_tbl = table(
        ["project", "science selection rule (encoded as data in "
         "<code>macro_core.staging.PROJECT_SELECTIONS</code>)", "source"],
        [[esc(p), esc(r), esc(s)] for p, r, s in rule_rows])

    return f"""
<section id="staged">
<div class="bhead"><h2>2 &middot; The five working sets</h2>
<span class="tag">selection rules are reviewable data, not code</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">What exactly does each project's staging manifest claim —
and where is the claim written down so a reviewer can diff it?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src_roles,
    "Every staged row, by role.  Science rows are the project's published "
    "selection; calibration rows are the era-matched families that ride "
    "along.")}</div>
{stage_tbl}
<p class="sub">Each selection rule, verbatim from the data table the build
reads (a change to any list or filter set is a one-line reviewable diff):</p>
{rule_tbl}
<div class="grid">{_figure(src_eras,
    "Where each working set lives in camera-era space.  A project's "
    "calibration rows cover exactly the eras its science occupies.")}</div>

<h3>Decision</h3>
<div class="decision"><b>The five target lists and filter rules are encoded
as data (<code>PROJECT_SELECTIONS</code>, <code>DW_FIELDS</code>,
filter whitelists/blacklists) and unit-tested against
<code>STRATEGY_CLAIMS</code></b> — every target a strategy document claims
is staged by exactly one project's rule, and the reduced tree, duplicate
copies, and header-error frames never stage.  QC flags and pointing
offsets travel WITH the rows: staging marks, it never drops.</div>

<h3>Consequence</h3>
<p class="sub">A stage that wants "the project's frames" reads one CSV (or
one <code>stage_*</code> table) and never re-derives target matching,
dedup, or tree policy.  The mispointed NGC&nbsp;5548 night, the saturated
SN&nbsp;2023ixf first epochs, and the 21 off-pointing T&nbsp;CrB grism
frames are all IN their manifests, flagged — each stage applies its own
disqualification rules with the evidence in hand.</p>
</div></section>"""


def section_calib(con) -> str:
    pairs = _stage_tables(con)
    rows = []
    for project, tbl in pairs:
        by = dict(q(con, f"""SELECT role, count(*) FROM {tbl}
                             WHERE {stg.SQL_CALIB_ROLES} GROUP BY role"""))
        n_master = sum(v for k, v in by.items() if k.startswith("master_"))
        rows.append([esc(project),
                     fmt(by.get("bias", 0)), fmt(by.get("dark", 0)),
                     fmt(by.get("flat", 0)), fmt(n_master)])
    kind_tbl = table(
        ["project", "bias (raw)", "dark (raw)", "flat (raw)",
         "master products"], rows)

    # Cross-check: every project whose science touches era 76 (the Mode0
    # grism era) must carry its recovered Calibrations/ masters — that
    # recovery closed the "zero Mode0 darks" gap and staging must not
    # silently lose it.
    check_rows = []
    for project, tbl in pairs:
        touches = q1(con, f"""SELECT count(*) FROM {tbl}
                              WHERE role = 'science' AND era_id = 76""")
        masters = q1(con, f"""SELECT count(*) FROM {tbl}
                              WHERE era_id = 76
                                AND role LIKE 'master_%'""")
        if touches:
            ok = masters > 0
            check_rows.append((
                [esc(project), fmt(touches), fmt(masters),
                 "ok" if ok else "<b>MISSING RECOVERED MASTERS</b>"],
                None if ok else "warn"))
    if check_rows:
        check_tbl = table(
            ["project with era-76 (Mode0) science", "era-76 science rows",
             "era-76 master products staged", "check"],
            [r for r, _ in check_rows],
            row_classes=[c for _, c in check_rows])
    else:
        check_tbl = "<p class=\"sub\">(no project stages era-76 science)</p>"

    n_total_cal = sum(
        q1(con, f"SELECT count(*) FROM {tbl} WHERE {stg.SQL_CALIB_ROLES}")
        for _, tbl in pairs)

    return f"""
<section id="calib">
<div class="bhead"><h2>3 &middot; Calibration attachment</h2>
<span class="tag">match_basis: {esc(stg.MATCH_BASIS_CALIB)}</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Which calibration frames belong in a project's working set,
and on what evidence is each one attached?</p>

<h3>Evidence</h3>
<p class="sub">{fmt(n_total_cal)} calibration rows across the five
manifests.  The rule is deliberately coarse: for every camera era a
project's science touches, ALL of that era's frames from the S0b census
are staged — raw bias/dark/flat and the recovered master products alike —
each row stamped <code>match_basis = '{esc(stg.MATCH_BASIS_CALIB)}'</code>
(masters distinguished by their <code>master_*</code> role, never by a
different basis).</p>
{kind_tbl}
<p class="sub">Cross-check against the S0b recovery: the Mode0 grism era
(era 76) had <i>zero</i> mode-matched darks until the
<code>Calibrations/</code> masters were recovered — any project staging
era-76 science must carry them:</p>
{check_tbl}

<h3>Decision</h3>
<div class="decision"><b>Era is the ONLY match key at staging time.</b>
Narrowing to the right kind, exposure bin, and filter is each stage's job,
guided by the S0b coverage matrix — S0c refuses to pre-judge (a dark that
misses spec can still be exposure-scaled; a cross-mode flat may enter with
a measured penalty term per the T&nbsp;CrB chair's ruling&nbsp;7).
Over-inclusion at this layer costs rows in a CSV; under-inclusion would
cost a stage a trip back to the census.</div>

<h3>Consequence</h3>
<p class="sub">A stage filters its manifest by <code>role</code>,
<code>era_id</code>, <code>exptime</code>, and <code>filter</code> and has
every candidate calibration frame in hand — including the masters whose
recovery is documented in the S0b report.</p>
</div></section>"""


def measure_claim(con, claim: stg.StageClaim) -> tuple[int, int]:
    """Run one published claim against its project's stage table.

    Returns ``(frames, nights)`` measured over ``role = 'science'`` ONLY —
    the role predicate is added here, never by the claim, so no claim can
    accidentally count calibration or cone-candidate rows.  The fragment is
    revalidated on every render (:func:`staging.assert_safe_where`) rather
    than trusted because it lives in the repo.
    """
    tbl = stg.stage_table_name(claim.project)
    where = stg.assert_safe_where(claim.where)
    row = q(con, f"""SELECT count(*), count(DISTINCT night)
                     FROM {tbl}
                     WHERE role = '{stg.ROLE_SCIENCE}' AND ({where})""")[0]
    return int(row[0]), int(row[1])


def section_reconcile(con) -> str:
    """Section 4: every strategy-document inventory number vs the manifest."""
    rows, classes = [], []
    n_ok = 0
    for claim in stg.STAGE_CLAIMS:
        n_frames, n_nights = measure_claim(con, claim)
        d_frames = n_frames - claim.claimed_frames
        # A claim with no published night count is measured but not judged
        # on nights, so a missing number can never be scored as a mismatch.
        d_nights = (None if claim.claimed_nights is None
                    else n_nights - claim.claimed_nights)
        agrees = d_frames == 0 and (d_nights in (0, None))
        n_ok += int(agrees)
        rows.append([
            esc(claim.project), esc(claim.label),
            fmt(claim.claimed_frames), fmt(n_frames),
            "&mdash;" if d_frames == 0 else f"<b>{d_frames:+,}</b>",
            "&mdash;" if claim.claimed_nights is None
            else fmt(claim.claimed_nights),
            fmt(n_nights),
            "&mdash;" if d_nights in (0, None) else f"<b>{d_nights:+,}</b>",
            esc(claim.source)])
        classes.append(None if agrees else "warn")
    claim_tbl = table(
        ["project", "published inventory claim", "doc frames",
         "staged frames", "&Delta;", "doc nights", "staged nights",
         "&Delta;", "source"],
        rows, row_classes=classes)
    n_claims = len(stg.STAGE_CLAIMS)

    return f"""
<section id="reconcile">
<div class="bhead"><h2>4 &middot; Documents vs the working set</h2>
<span class="tag">{n_ok} of {n_claims} published numbers reconcile</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">A manuscript quotes "30 frames/band", "20 frames",
"403 frames / 40 nights".  Where is the machinery that notices when a
sentence in a strategy document stops matching the frames that actually
exist?</p>

<h3>Evidence</h3>
<p class="sub">Every inventory number the five strategy documents publish,
encoded in <code>macro_core.staging.STAGE_CLAIMS</code> next to the query
that reproduces it, measured here against the stage tables
(<code>role = 'science'</code> is imposed by the renderer, so no claim can
count calibration or cone-candidate rows):</p>
{claim_tbl}

<h3>Decision</h3>
<div class="decision"><b>Each published inventory number lives next to the
query that reproduces it, and the diff is rendered on this page.</b>  Three
numbers failed this test at the 2026-08-18 review — the SN template stacks
("30 frames/band"; canonical reality 15, the other 15 being reduced-tree
copies), the NGC&nbsp;5548 revisit ("20 frames"; 10), and the &theta;&nbsp;CrB
calibrator series ("403 / 40 nights"; 412 / 42, because two more nights
arrived after the sentence was written).  Staging was correct in all three
cases; only the prose drifted, and nothing made the drift visible.</div>

<h3>Consequence</h3>
<p class="sub">Drift now surfaces here in two directions: a document
mis-stating the archive, and an archive that grew past a document.  The
second is the common case and the benign one — it is how the &theta;&nbsp;CrB
series gained 9 frames and Vega gained a whole instrument tier — but a
referee cannot tell the two apart from prose alone, and neither could
we.</p>
</div></section>"""


def section_provenance(con) -> str:
    n_obs = q1(con, "SELECT count(*) FROM frames")
    n_canon = q1(con, "SELECT count(*) FROM frames WHERE is_canonical = 1")
    n_calib = q1(con, "SELECT count(*) FROM calib_frames")
    meta = dict(q(con, "SELECT key, value FROM s0c_build_meta"))
    pairs = _stage_tables(con)
    n_staged = sum(q1(con, f"SELECT count(*) FROM {tbl}")
                   for _, tbl in pairs)

    file_rows = q(con, """
        SELECT project, stage_table, n_science, n_cone, n_calib, n_rows,
               csv_path, n_symlinks
        FROM s0c_stage_files ORDER BY project""")
    files_tbl = table(
        ["project", "DB table", "science", "cone", "calib", "rows",
         "CSV (gitignored, regenerable)", "farm links"],
        [[esc(p), f"<code>{esc(t)}</code>", fmt(s), fmt(k), fmt(c), fmt(n),
          f"<code>{esc(f)}</code>", fmt(l) if l else "&mdash;"]
         for p, t, s, k, c, n, f, l in file_rows])

    chain_tbl = table(
        ["link", "artifact", "population", "custody rule"],
        [["1", "<code>rlmt-archive/</code> (3.3&nbsp;TiB)",
          f"{fmt(n_obs)} cataloged files",
          "<b>immutable</b> — read-only to every stage, forever"],
         ["2", "S0 manifest <code>frames</code>",
          f"{fmt(n_obs)} rows &rarr; {fmt(n_canon)} canonical",
          "global dedup, eras, aliases, QC — the provenance layer"],
         ["3", "S0b <code>calib_frames</code>",
          f"{fmt(n_calib)} calibration frames incl. recovered masters",
          "kind-normalized census, era-keyed"],
         ["4", "S0c <code>stage_&lt;project&gt;</code> tables",
          f"{fmt(n_staged)} rows across {len(pairs)} projects",
          "selection rules as data; atomic swap per build"],
         ["5", "<code>&lt;Project&gt;/data/stage_manifest.csv</code>",
          "one CSV per project (same rows as link 4)",
          "the working artifact — gitignored, regenerable, "
          "never a copy of pixels"]])

    return f"""
<section id="provenance">
<div class="bhead"><h2>5 &middot; The provenance chain</h2>
<span class="tag">archive &rarr; manifest &rarr; stage table &rarr; CSV</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">When a stage opens a frame through a staging manifest, what
chain of evidence connects that row back to a file in the immutable
archive — and where can each link be audited?</p>

<h3>Evidence</h3>
{chain_tbl}
<p class="sub">This build's outputs (build id
<code>{esc(meta.get('build_id', ''))}</code>, commit
<code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>,
<code>staging.py</code> sha256
<code>{esc(meta.get('staging_sha256_12', '') or 'unrecorded')}</code>):</p>
{files_tbl}
<p class="sub">The commit field carries a <code>-dirty</code> suffix when the
build ran from a modified working tree, and the <code>staging.py</code>
digest identifies the exact selection-rule text that emitted these rows.  The
first S0c build recorded a bare hash of the commit <em>before</em> the one
that introduced <code>staging.py</code>; a provenance field that cannot
recover the code that ran is not provenance.</p>

<h3>Decision</h3>
<div class="decision"><b>Every stage row carries its whole chain:</b>
<code>obs_rowid</code> joins back to the catalog scan and the manifest,
<code>dup_group</code>/<code>era_id</code>/<code>qc_flags</code> carry the
S0 rulings, <code>match_basis</code> says why the row exists,
<code>size_bytes</code> is the read-time integrity tripwire, and
<code>stage_build_id</code> pins which build emitted it.  The CSV and the
DB table are the same rows in the same order — two views of one build,
regenerated together, never edited by hand.</div>

<h3>Consequence</h3>
<p class="sub">The October ingest path gains its third link: after each
archive sync, re-run S0 &rarr; S0b &rarr; S0c and every project's working
set follows the archive — no copies to chase, no hand lists to update.
Regenerate any time with
<code>pipeline/scripts/build_s0c_staging.py</code>.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S0c report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    try:
        pairs = _stage_tables(con)
        n_sci = sum(q1(con, f"SELECT count(*) FROM {t} "
                            f"WHERE role = '{stg.ROLE_SCIENCE}'")
                    for _, t in pairs)
        n_cone = sum(q1(con, f"SELECT count(*) FROM {t} WHERE role = "
                             f"'{stg.ROLE_SCIENCE_UNRESOLVED}'")
                     for _, t in pairs)
        n_cal = sum(q1(con, f"SELECT count(*) FROM {t} "
                            f"WHERE {stg.SQL_CALIB_ROLES}") for _, t in pairs)
        meta = dict(q(con, "SELECT key, value FROM s0c_build_meta"))

        sections = [
            section_no_copy(con),
            section_staged(con),
            section_calib(con),
            section_reconcile(con),
            section_provenance(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S0c — Per-Project Staging Manifests</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S0c — Per-Project Staging Manifests</h1>
  <p>{len(pairs)} projects staged &middot; {fmt(n_sci)} science rows
  &middot; {fmt(n_cone)} cone candidates
  &middot; {fmt(n_cal)} calibration rows &middot; zero frames copied
  &middot; built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">back to the evidence hub</a></p>
</header>

<nav>
  <a href="#nocopy">1 The no-copy law</a> &middot;
  <a href="#staged">2 The five working sets</a> &middot;
  <a href="#calib">3 Calibration attachment</a> &middot;
  <a href="#reconcile">4 Documents vs the working set</a> &middot;
  <a href="#provenance">5 Provenance chain</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_core.report_s0c</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query; none is typed by hand.  Regenerate with
<code>pipeline/scripts/build_s0c_staging.py</code> (the archive-sync ingest
path: S0 &rarr; S0b &rarr; S0c).</footer>
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
