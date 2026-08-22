"""Project page renderer — the plan, the progress, and what's next.

Writes ``docs/<Project>/index.html`` for every project in
:data:`macro_core.project_plan.PROJECTS`, plus the full evidence index
``docs/evidence.html``.

These pages are the **Plan & Status** view of each project.  The three
navigation layers around them — the project row, the view tabs, the sticky
question rail — are put there afterwards by :mod:`macro_core.site`, which
re-reads what this module wrote and wraps it.  Nothing here needs to know
that; keep emitting the same markup and the chrome will fit.

WHY THIS EXISTS
---------------
The project pages used to be hand-written readiness paragraphs.  They went
stale the way hand-written things do — quoting frame counts that had since
moved, detector constants whose tables had since been destroyed, and an
astrometry claim ("EU UMa's frames can never be plate-solved") that a later
geometry re-characterisation overturned outright.  Worse, none of them
answered the question a reader actually arrives with: **what is the whole
plan, and where in it are we?**

So every page here is generated, and the two halves come from two places
that cannot drift:

* the PLAN from :mod:`macro_core.project_plan` — phases, tasks, the strategy
  section each was derived from;
* every NUMBER from a SQL query against the manifest, and every stage
  verdict from :mod:`macro_core.provenance`.

Nothing on these pages is typed by hand except prose that is explicitly
about the plan.  A detector constant whose backing table was destroyed is
shown as UNBACKED rather than quoted as fact — the pages say what has a
query behind it and what does not, which is the only way a reader can
calibrate how much to trust the rest.

PAGE STRUCTURE (fixed, in this order)
-------------------------------------
a. what the paper will claim, and where it will go
b. PROGRESS AT A GLANCE — counts, a bar, and live headline numbers
c. THE PLAN — every phase, every task, status, product, evidence
d. WHAT'S BLOCKING — each blocker and what would clear it
e. NEXT UP — the next actionable tasks
f. EVIDENCE — the pipeline reports this project rests on, each with its
   CURRENT provenance verdict
g. a footer stating when the page was generated and from what
"""

from __future__ import annotations

import re
import sqlite3
from pathlib import Path
from typing import Mapping, Optional, Sequence

from . import project_plan as pp
from . import provenance as pv
from .report_s0 import esc, fmt, q, q1, table

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"

#: The full evidence index — every project and every pipeline stage, each
#: with its live verdict.
#:
#: This used to be ``docs/index.html``, i.e. the front door.  It is the wrong
#: front door: it opens on two long tables and sixteen stage rows, which is
#: exactly the "quite a load of stuff to crawl through" a colleague met when
#: the site was circulated.  ``docs/index.html`` is now the LANDING page
#: written by :mod:`macro_core.site`, which orients a reader in thirty
#: seconds and links here for the ones who want the audit.  Same content,
#: one door further in.  Both files sit in ``docs/`` so every relative link
#: computed against this path is unchanged.
HUB_PATH = DOCS_DIR / "evidence.html"

#: Stage table per project, for the live headline numbers.  Keyed off the
#: manifest's own registry (s0c_stage_files) at render time; this map is
#: only the fallback ordering for projects that have no staged rows yet.
NO_STAGE_TABLE = ""

#: Colour class per provenance verdict.  DESTROYED is its own class: a
#: reader must be able to tell "the numbers moved" from "there is no table".
VERDICT_CLASS = {
    pv.FRESH: "v-fresh",
    pv.STALE: "v-stale",
    pv.STALE_UPSTREAM: "v-wait",
    pv.NEVER_RUN: "v-never",
    pv.OUTPUT_MISSING: "v-gone",
    "DESTROYED": "v-gone",
    "UNKNOWN": "v-wait",
}

#: The project pages once carried their whole presentation here, inlined,
#: so this renderer "owns its own presentation and cannot break the pipeline
#: reports that share the stylesheet".  That reasoning inverted the moment
#: the site became one system: two definitions of a status chip is two
#: chances for `done` to be a different green on two pages a reader compares
#: side by side, and the inline copy was the dark-ground one, so it survived
#: the move to a light theme by turning pale green text onto white.
#:
#: Every class this block defined — `.chip`, `.bar`, `.stat`, `.phase`,
#: `.blockcard`, `.nextcard`, `.unbacked`, `.src` — now lives once, in
#: `docs/assets/macro.css`, with the colour-blind-safe status palette.  This
#: constant stays, empty, because it is a documented seam: a page-specific
#: rule that genuinely belongs to one renderer goes here, and anything that
#: describes a SHARED component belongs in the stylesheet instead.
PAGE_CSS = ""


# ---------------------------------------------------------------------------
# Render context — computed ONCE and shared by all seven pages
# ---------------------------------------------------------------------------
class Context:
    """Everything the pages need, gathered once.

    The freshness computation fingerprints every declared resource in the
    DAG; doing it per page would hash the manifest seven times over a
    spinning archive drive for seven identical answers.
    """

    def __init__(self, con: sqlite3.Connection, repo_root: Path):
        self.con = con
        self.repo_root = repo_root
        self.generated = pp.utcnow()
        self.freshness, self.fingerprints = pp.stage_freshness(con, repo_root)
        self.recorded = pp.read_statuses(con)
        self.recorded_evidence = pp.read_evidence(con)
        # Frame-level lookups, as SETS rather than SQL joins: the manifest
        # has no index on stage_*.obs_rowid, and the join form of these two
        # questions takes minutes on the archive drive while the set form
        # takes under a second.
        self.solved: set[int] = {
            r[0] for r in q(con, "SELECT obs_rowid FROM s1_batch "
                                 "WHERE status = 'solved'")}
        self.solved |= {r[0] for r in q(con, "SELECT obs_rowid FROM frames "
                                             "WHERE pltsolvd = 1")}
        self.timed: set[int] = {
            r[0] for r in q(con, "SELECT obs_rowid FROM frame_times")}
        # S2c's per-frame verdicts, so a page can say what a FILTER actually
        # does rather than what its name suggests.  Same set-lookup reason.
        # strength_class comes along because a verdict's CONFIDENCE changes
        # what it licenses: portfolio-wide, most 'dispersed' verdicts are
        # 'low', and a MIXED flag driven entirely by low-strength frames is
        # a different statement from one driven by confident ones.
        self.dispersion: dict[int, tuple[str, str, str]] = {
            r[0]: (r[1] or "(blank)", r[2], r[3] or "n/a")
            for r in q(con, "SELECT obs_rowid, filter, verdict, "
                            "strength_class FROM frame_dispersion "
                            "WHERE verdict IS NOT NULL")}
        self.stage_tables: dict[str, str] = dict(
            q(con, "SELECT project, stage_table FROM s0c_stage_files"))
        # Which stages the manifest can prove once produced something, so
        # DESTROYED is only ever said about a stage that actually ran.
        self.ever_ran = pp.stages_ever_run(con)
        # The strategy documents themselves, read once.  The reconciliation
        # table checks its own transcriptions against them at render time.
        self.documents = pp.read_cited_documents(repo_root)

    def verdict(self, stage_key: str) -> tuple[str, str]:
        """``(verdict, why)`` for one stage — always with ``ever_ran``, so no
        caller can accidentally label a never-built stage DESTROYED."""
        return pp.evidence_verdict(stage_key, self.freshness,
                                   self.fingerprints, self.ever_ran)

    def statuses(self, tasks: Sequence[pp.Task]) -> dict[str, str]:
        return pp.overlay_statuses(tasks, self.recorded)

    def evidence(self, task: pp.Task) -> str:
        """Recorded evidence wins over the ledger's — the CLI's --evidence
        is how a person points at what they actually produced."""
        return self.recorded_evidence.get(task.id, task.evidence)


# ---------------------------------------------------------------------------
# Small formatting helpers
# ---------------------------------------------------------------------------
def _chip(cls: str, text: str) -> str:
    return f'<span class="chip {cls}">{esc(text)}</span>'


def _status_chip(status: str) -> str:
    return _chip(status, pp.STATUS_LABEL.get(status, status))


def _verdict_chip(verdict: str) -> str:
    return _chip(VERDICT_CLASS.get(verdict, "v-wait"), verdict)


def _link(target_page: Path, repo_rel: str, text: str = "") -> str:
    """A link from a project page to a repo-relative artifact.

    Only paths under ``docs/`` are linkable from a published page; anything
    else (a product database, an ops memo) is shown as a code path, because
    a link a reader cannot follow is worse than an honest filename.
    """
    label = text or Path(repo_rel).name
    if repo_rel.startswith("docs/"):
        rel = Path(repo_rel).relative_to("docs")
        depth = len(target_page.relative_to(DOCS_DIR).parts) - 1
        prefix = "../" * depth
        return f'<a href="{prefix}{rel.as_posix()}">{esc(label)}</a>'
    return f"<code>{esc(repo_rel)}</code>"


def _report_link(page: Path, stage_key: str, ctx: "Context") -> str:
    """A link to a stage's published report, WITH that stage's verdict.

    The verdict never used to travel with the link.  It lived in a different
    column of a table further down the page, so a reader clicking through
    from the Plan table arrived at ``s2_detector.html`` — which states the
    High Gain ceiling, the veto threshold and the read noise as measured
    values, five times over, with no banner — carrying no warning at all
    that every table behind those numbers has been wiped.  The chip is
    cheap; a reader quoting an orphaned constant is not.
    """
    report = pp.stage_report_path(stage_key)
    if not report:
        return "&mdash;"
    verdict, _ = ctx.verdict(stage_key)
    link = _link(page, report)
    if verdict == pv.FRESH:
        return link
    return f'{link} {_verdict_chip(verdict)}'


def _source_link(source: pp.Source) -> str:
    """Cite the strategy section a task came from, linking to GitHub for
    documents that are not published on this site."""
    if source.document.startswith("docs/"):
        return f'<span class="src">{esc(source.document)} {esc(source.section)}</span>'
    url = ("https://github.com/jwwetzel/MACRO/blob/main/"
           + source.document)
    return (f'<span class="src"><a href="{url}">{esc(source.document)}</a> '
            f'{esc(source.section)}</span>')


# ---------------------------------------------------------------------------
# (b) live headline numbers
# ---------------------------------------------------------------------------
def _stage_table_of(project_key: str, ctx: Context) -> Optional[str]:
    return ctx.stage_tables.get(project_key)


def _target_rows(ctx: Context, stage_table: str) -> list[tuple]:
    """``(target, frames, nights, wcs, bjd)`` per staged science target.

    Every number here is counted from the staging table — the project's
    actual working set — not transcribed from a strategy document.
    """
    rows = q(ctx.con,
             f"SELECT obs_rowid, canonical_target, night FROM {stage_table} "
             f"WHERE role = 'science'")
    per: dict[str, dict] = {}
    for rowid, target, night in rows:
        name = target or "(unnamed)"
        rec = per.setdefault(name, {"n": 0, "nights": set(),
                                    "wcs": 0, "bjd": 0})
        rec["n"] += 1
        if night:
            rec["nights"].add(night)
        if rowid in ctx.solved:
            rec["wcs"] += 1
        if rowid in ctx.timed:
            rec["bjd"] += 1
    out = [(name, r["n"], len(r["nights"]), r["wcs"], r["bjd"])
           for name, r in per.items()]
    out.sort(key=lambda t: -t[1])
    return out


#: A staged target name carrying no more than this many frames is almost
#: certainly an unmerged alias of a neighbouring row, not a separate object.
#: Two of the SN project's five "targets" were single frames named
#: '2023ixf1' and '2023ixf2'.
_SPLINTER_MAX_FRAMES = 2


def _alias_note(rows: Sequence[tuple]) -> str:
    """Name the target rows that look like unmerged aliases.

    Derived, not asserted: it flags rows by frame count, says why they are
    suspicious, and leaves the astronomy to the alias task rather than
    silently merging anything.
    """
    splinters = [(name, n) for name, n, *_ in rows
                 if n <= _SPLINTER_MAX_FRAMES]
    if not splinters:
        return ""
    listed = ", ".join(f"<code>{esc(name)}</code> ({fmt(n)} frame"
                       f"{'' if n == 1 else 's'})" for name, n in splinters)
    return (f'<div class="blockcard"><b>Unmerged alias candidates.</b> '
            f'{listed} — each holds at most {_SPLINTER_MAX_FRAMES} frames, '
            f'which is what an alias fragment looks like, not a target. '
            f'They are counted in "catalog target names" above because that '
            f'stat counts names; they should not be read as separate '
            f'objects. Merging them is the alias task\'s job, and until it '
            f'runs this page will keep saying so.</div>')


def _stat(n: str, k: str) -> str:
    return f'<div class="stat"><span class="n">{n}</span>' \
           f'<span class="k">{esc(k)}</span></div>'


def section_progress(project: pp.Project, ctx: Context) -> str:
    tasks = project.tasks
    statuses = ctx.statuses(tasks)
    counts = pp.status_counts(tasks, statuses)
    done, total = pp.progress_fraction(counts)

    order = (pp.DONE, pp.IN_PROGRESS, pp.REDO_NEEDED, pp.BLOCKED, pp.PENDING)
    segments = "".join(
        f'<span class="{s}" style="width:{100.0 * counts[s] / total:.4f}%"></span>'
        for s in order if counts[s])
    legend = " &middot; ".join(
        f"<b>{counts[s]}</b> {pp.STATUS_LABEL[s]}" for s in order if counts[s])

    stage_table = _stage_table_of(project.key, ctx)
    stats = [_stat(f"{done}/{total}", "plan tasks complete")]
    body = ""
    if stage_table:
        rows = _target_rows(ctx, stage_table)
        n_sci = sum(r[1] for r in rows)
        n_wcs = sum(r[3] for r in rows)
        n_bjd = sum(r[4] for r in rows)
        n_cal = q1(ctx.con, f"SELECT count(*) FROM {stage_table} "
                            f"WHERE role NOT IN ('science', "
                            f"'science_unresolved')")
        n_nights = q1(ctx.con, f"SELECT count(DISTINCT night) FROM "
                               f"{stage_table} WHERE role = 'science'")
        n_eras = q1(ctx.con, f"SELECT count(DISTINCT era_id) FROM "
                             f"{stage_table} WHERE role = 'science' "
                             f"AND era_id IS NOT NULL")
        stats += [
            _stat(fmt(n_sci), "staged science frames"),
            _stat(fmt(n_nights), "nights"),
            # NOT "targets".  This counts distinct canonical_target STRINGS,
            # and on the SN page two of the five were one-frame alias
            # splinters of the supernova itself ('2023ixf1', '2023ixf2')
            # while two more ('M101', 'pinwheel galaxy') were two catalog
            # names for one host.  Name fragmentation is a named failure
            # mode in these strategies; rendering it as "5 targets" turned a
            # bookkeeping artifact into a fact about the sky.
            _stat(fmt(len(rows)), "catalog target names"),
            _stat(f"{100.0 * n_wcs / n_sci:.0f}%" if n_sci else "&mdash;",
                  f"carry a WCS ({fmt(n_wcs)}/{fmt(n_sci)})"),
            _stat(f"{100.0 * n_bjd / n_sci:.0f}%" if n_sci else "&mdash;",
                  f"carry BJD_TDB ({fmt(n_bjd)}/{fmt(n_sci)})"),
            _stat(fmt(n_cal), "era-matched calibration rows"),
            _stat(fmt(n_eras), "camera eras spanned"),
        ]
        body = table(
            ["Target", "Staged science frames", "Nights", "With WCS",
             "With BJD_TDB"],
            [[esc(name), fmt(n), fmt(nights), f"{fmt(w)} ({100.0 * w / n:.0f}%)",
              f"{fmt(b)} ({100.0 * b / n:.0f}%)"]
             for name, n, nights, w, b in rows])
        body = (f'<p class="sub">Every row is counted from '
                f'<code>{esc(stage_table)}</code> at render time — the '
                f'project\'s actual working set, not a number transcribed '
                f'from a strategy document. Each row is one '
                f'<code>canonical_target</code> STRING, which is not the '
                f'same as one object.</p>{body}{_alias_note(rows)}')
    else:
        body = ('<p class="sub">This project has no staging table in the '
                'manifest yet, so there are no pipeline-emitted numbers to '
                'show. Any frame count you have seen for it comes from the '
                'archive page, not from a query.</p>')

    claims = _claims_table(project, ctx)
    filters = _filter_identity(project, ctx, stage_table)

    return f"""
<section id="progress">
<div class="bhead"><h2>Progress at a glance</h2>
<span class="tag">{done} of {total} plan tasks complete</span></div>
<div class="stage">
<div class="bar">{segments}</div>
<p class="legend">{legend}</p>
<div class="stats">{"".join(stats)}</div>
{body}
{filters}
{claims}
</div>
</section>"""


#: S2c's PUBLISHED classification rule, adopted here verbatim.
#:
#: This renderer used to flag a filter MIXED whenever it held at least one
#: dispersed frame and at least one direct one — no threshold, no minimum
#: count.  S2c's own report page, on this same site, has always used ≥80%
#: dispersed → SPECTRA, ≥80% direct → images, otherwise MIXED, over filters
#: with at least 20 measured frames.  So the same word meant two different
#: things two pages apart, and the looser one was the one shown to project
#: readers: Johnson I was rendered bold MIXED on the strength of ONE
#: dispersed frame out of 68 (96% direct — "images" by the published rule),
#: Sloan i on four of 58, and NGC 5548's W slot on one of 47.  Every one of
#: those dispersed verdicts carried strength_class 'low'.  A photometry
#: project was being told its Sloan i and Johnson I slots were contaminated
#: with spectra on the evidence of single low-confidence frames.
#:
#: One word, one meaning, site-wide.  If S2c's rule changes, this constant
#: and the test that pins it to `report_s2c` change together.
S2C_MAJORITY = 0.8

#: S2c reports per-filter verdicts only where at least this many frames were
#: actually measured; below it the fractions are noise.
S2C_MIN_MEASURED = 20


def classify_filter(counts: Mapping[str, int]) -> tuple[str, str]:
    """``(verdict, css class)`` for one filter, by S2c's published rule.

    Pure, so the test can drive it with the exact counts that produced the
    old false alarms and pin the answer.
    """
    total = sum(counts.values())
    if total < S2C_MIN_MEASURED:
        return ("too few measured", "")
    dispersed = counts.get("dispersed", 0)
    direct = counts.get("direct", 0)
    if dispersed / total >= S2C_MAJORITY:
        return ("SPECTRA", "warn")
    if direct / total >= S2C_MAJORITY:
        return ("images", "")
    return ("MIXED", "warn")


def _filter_identity(project: pp.Project, ctx: Context,
                     stage_table: Optional[str]) -> str:
    """What this project's FILTERS actually do, measured frame by frame.

    Every strategy in this portfolio contains at least one rule of the form
    "filter x is a spectrum, so exclude it" — asserted from the filter's
    NAME.  S2c measures it instead, and a slot that turns out to be mixed
    (dispersed on one target, direct on another) is a fact no name can carry.
    Rendered from the measurements, so the page cannot assert an identity the
    database has not confirmed — and by :data:`S2C_MAJORITY`, so "MIXED"
    here means what "MIXED" means on S2c's own page.
    """
    if not stage_table:
        return ""
    per = _dispersion_by_filter(ctx, stage_table)
    if not per:
        return ""

    kinds = ("dispersed", "direct", "indeterminate")
    body, classes = [], []
    for filt in sorted(per, key=lambda f: -sum(per[f]["verdicts"].values())):
        counts = per[filt]["verdicts"]
        total = sum(counts.values())
        verdict, cls = classify_filter(counts)
        low = per[filt]["low"]
        # A MIXED (or SPECTRA) verdict resting on low-strength frames is a
        # weaker statement than the same verdict resting on confident ones,
        # and the reader cannot tell without this column.
        strength = (f"{fmt(low)} of {fmt(counts.get('dispersed', 0))} low"
                    if counts.get("dispersed") else "&mdash;")
        body.append([f"<code>{esc(filt)}</code>", fmt(total)]
                    + [fmt(counts.get(k, 0)) for k in kinds]
                    + [strength,
                       f"<b>{esc(verdict)}</b>" if cls else esc(verdict)])
        classes.append(cls)

    return (f'<h3>What this project&rsquo;s filters actually do</h3>'
            f'<p class="sub">Per-frame verdicts from S2c, restricted to this '
            f'project&rsquo;s staged science frames. The measurement is '
            f'source elongation, not the filter&rsquo;s name — which is the '
            f'point, since every strategy here carries an exclusion rule '
            f'keyed on the name. The verdict column uses <b>S2c&rsquo;s own '
            f'published rule</b> — ≥{S2C_MAJORITY:.0%} dispersed is SPECTRA, '
            f'≥{S2C_MAJORITY:.0%} direct is images, anything else is MIXED, '
            f'and only filters with at least {S2C_MIN_MEASURED} measured '
            f'frames get a verdict at all — so the word means the same thing '
            f'here as it does on <a href="../pipeline/s2c_filter_identity.html">'
            f'the S2c page</a>. The strength column matters: a verdict built '
            f'from low-confidence frames is a weaker claim, and most '
            f'dispersed verdicts in this archive are low. Frames S2c has not '
            f'reached yet simply do not appear.</p>'
            + table(["Filter", "Measured", "Dispersed", "Direct",
                     "Indeterminate", "Dispersed strength", "Verdict"],
                    body, row_classes=classes))


def _dispersion_by_filter(ctx: Context, stage_table: str) -> dict[str, dict]:
    """``{filter: {"verdicts": {verdict: n}, "low": n}}`` over staged science
    frames.  Shared by the filter panel and the claim panel so the two can
    never disagree about the same frames."""
    rows = q(ctx.con, f"SELECT obs_rowid FROM {stage_table} "
                      f"WHERE role = 'science'")
    per: dict[str, dict] = {}
    for (rowid,) in rows:
        hit = ctx.dispersion.get(rowid)
        if hit is None:
            continue
        filt, verdict, strength = hit
        rec = per.setdefault(filt, {"verdicts": {}, "low": 0})
        rec["verdicts"][verdict] = rec["verdicts"].get(verdict, 0) + 1
        if verdict == "dispersed" and strength == "low":
            rec["low"] += 1
    return per


def _claim_measurement(project: pp.Project, ctx: Context,
                       stage_table: Optional[str]) -> str:
    """The measurement the CLAIM rests on, printed under the claim.

    A claim paragraph is what a reader quotes.  When it rests on a filter
    slot being dispersed, the honest place for S2c's verdict is right there,
    not a hundred lines down in a table the quoter never reached.
    """
    if not project.claim_filters or not stage_table:
        return ""
    per = _dispersion_by_filter(ctx, stage_table)
    cards = []
    for filt in project.claim_filters:
        rec = per.get(filt)
        if not rec:
            cards.append(
                f'<div class="blockcard">Slot <code>{esc(filt)}</code>: '
                f'S2c has not measured any of this project&rsquo;s staged '
                f'frames in this slot yet, so the claim above rests on the '
                f'slot&rsquo;s NAME.</div>')
            continue
        counts = rec["verdicts"]
        total = sum(counts.values())
        verdict, cls = classify_filter(counts)
        d = counts.get("dispersed", 0)
        cards.append(
            f'<div class="{"blockcard" if cls else "nextcard"}">'
            f'<b>What the claim rests on, measured.</b> Slot '
            f'<code>{esc(filt)}</code>: of {fmt(total)} staged science '
            f'frames S2c has measured, <b>{fmt(d)}</b> are dispersed, '
            f'{fmt(counts.get("direct", 0))} measure as <b>direct '
            f'imaging</b>, and {fmt(counts.get("indeterminate", 0))} are '
            f'undecided — S2c&rsquo;s verdict for the slot is '
            f'<b>{esc(verdict)}</b>'
            + (f', and {fmt(rec["low"])} of the {fmt(d)} dispersed verdicts '
               f'are low-strength' if d else '')
            + '. The claim above is worth what these numbers say it is '
              'worth, and no more.</div>')
    return "".join(cards)


#: The quoted fragment inside a project_counts `source` string — S0's own
#: transcription of the sentence it read out of the strategy.
_QUOTED_RE = re.compile(r"'([^']{4,})'")


def _transcription_is_current(source: str, document: str) -> Optional[bool]:
    """Does S0's quoted fragment still appear in the strategy document?

    ``None`` when the source string carries no quote to check.

    THIS IS THE GUARD THE RECONCILIATION TABLE WAS MISSING.  Every
    ``project_counts`` row is a transcription made when S0 last ran, and S0
    is STALE — so the table whose entire job is catching drift between the
    strategy and the pipeline was itself aging silently.  It reported
    θ CrB / grism_light as claims 403 / holds 403 / diff —, i.e. perfect
    agreement, on a number the strategy had already retracted: §3 now reads
    "412 unique rawimage frames, 42 nights" and labels "403 / 40 nights" as
    the previous revision, while the staging table holds 412 across 42
    nights.  Both ends had moved; only the frozen middle said they agreed.

    Checking the QUOTE rather than the number is what makes this work.  The
    bare digits "403" still appear in the document — inside the sentence
    retracting them — so a numeric search would report everything fine.  The
    quoted phrase "403 unique rawimage frames, 40 nights" does not.
    """
    m = _QUOTED_RE.search(source or "")
    if not m:
        return None
    return m.group(1) in (document or "")


def _claims_table(project: pp.Project, ctx: Context) -> str:
    """The strategy's own claimed counts against the manifest's.

    S0 records both; showing the diff is how a reader sees, without taking
    anyone's word for it, which strategy numbers the pipeline reproduces and
    which it does not.  Every row now carries the verdict of the stage that
    WROTE it, and is re-checked against the live strategy text — a
    reconciliation table that cannot go stale is the only kind worth
    publishing.
    """
    rows = q(ctx.con,
             "SELECT target, metric, claimed_frames, manifest_frames, "
             "diff_frames, source FROM project_counts WHERE project = ? "
             "ORDER BY target", (project.key,))
    if not rows:
        return ""

    s0_verdict, _ = ctx.verdict("S0")
    document = ctx.documents.get(project.strategy, "")

    body, classes, n_superseded = [], [], 0
    for target, metric, claimed, manifest, diff, source in rows:
        current = _transcription_is_current(source, document)
        if current is False:
            n_superseded += 1
            status = ('<b>SUPERSEDED</b><br><span class="src">the sentence '
                      'S0 transcribed is no longer in the strategy</span>')
            cls = "warn"
        elif diff:
            status = f"<b>{fmt(diff)}</b>"
            cls = "warn"
        elif current is None:
            status = 'agrees <span class="src">(quote unverifiable)</span>'
            cls = ""
        else:
            status = "agrees"
            cls = ""
        body.append([esc(target), f"<code>{esc(metric)}</code>",
                     fmt(claimed), fmt(manifest), status])
        classes.append(cls)

    note = ""
    if n_superseded:
        note = (f'<div class="unbacked"><b>{n_superseded} row(s) transcribe '
                f'a sentence the strategy no longer contains.</b> S0 copied '
                f'these numbers out of the strategy the last time it ran, '
                f'and the strategy has been revised since. A row marked '
                f'SUPERSEDED is not a disagreement between the plan and the '
                f'pipeline — it is a stale copy sitting between them, and '
                f'the number in the "strategy claims" column is not what '
                f'the strategy currently claims. Re-run S0 to refresh '
                f'them.</div>')
    elif s0_verdict != pv.FRESH:
        note = (f'<p class="sub">S0, which wrote these rows, is currently '
                f'<b>{esc(s0_verdict)}</b>. Every quoted sentence still '
                f'matches the strategy text, so the transcriptions are '
                f'current even though the stage is not.</p>')

    grid = table(["Target", "Metric", "Strategy claims", "Manifest holds",
                  "Reconciliation"], body, row_classes=classes)
    return (f'<h3>What the strategy claims vs what the manifest holds '
            f'{_verdict_chip(s0_verdict)}</h3>'
            f'<p class="sub">Reconciliation rows written by S0 — and '
            f're-checked at render time by looking for the sentence S0 says '
            f'it transcribed in the strategy document as it stands now. A '
            f'non-zero diff is not necessarily an error — it is a question '
            f'that has an answer in the <code>source</code> column of '
            f'<code>project_counts</code>. A SUPERSEDED row is a different '
            f'thing: the question itself has expired.</p>{grid}{note}')


# ---------------------------------------------------------------------------
# (c) the plan
# ---------------------------------------------------------------------------
def section_plan(project: pp.Project, ctx: Context, page: Path) -> str:
    statuses = ctx.statuses(project.tasks)
    blocks = []
    for phase in project.phases:
        rows = []
        classes = []
        for t in phase.tasks:
            status = statuses[t.id]
            ev = ctx.evidence(t)
            ev_cell = _link(page, ev) if ev else "&mdash;"
            unmet = pp.unmet_dependencies(t, statuses)
            gate = ""
            if unmet and status in pp.OPEN_STATUSES:
                names = ", ".join(f"<code>{esc(d)}</code>" for d in unmet)
                gate = (f'<br><span class="src"><b>Gated:</b> waiting on '
                        f'{names}</span>')
            rows.append([
                _status_chip(status),
                # The id is shown because it is the handle: it is what a
                # person types into `update_project_plan.py set`, and a plan
                # you cannot address is a plan you cannot update.
                f"<b>{esc(t.title)}</b><br><code>{esc(t.id)}</code><br>"
                f"{_source_link(t.source)}{gate}",
                esc(t.produces),
                f"<code>{esc(t.stage)}</code> "
                f"{_report_link(page, t.stage, ctx)}",
                ev_cell,
            ])
            classes.append("warn" if status in (pp.BLOCKED, pp.REDO_NEEDED)
                           or unmet else "")
        blocks.append(f"""
<div class="phase">
<h3>{esc(phase.name)}</h3>
<p class="sub">{esc(phase.intent)}</p>
{table(["Status", "Task", "What it produces", "Stage &amp; its report",
        "Evidence"], rows, row_classes=classes)}
</div>""")

    strategy = (f'derived from <a href="https://github.com/jwwetzel/MACRO/'
                f'blob/main/{project.strategy}">{esc(project.strategy)}</a>'
                if project.strategy else
                "derived from this page's own recorded decisions — there is "
                "no committee strategy for this archive yet")

    return f"""
<section id="plan">
<div class="bhead"><h2>The plan</h2>
<span class="tag">{len(project.tasks)} tasks in {len(project.phases)} phases</span></div>
<div class="stage">
<p class="sub">Every task below was read out of the project's execution
order — {strategy} — and names the pipeline stage its result rests on. That
last column is not decoration: when a stage stops being fresh, every
<span class="chip done">done</span> task resting on it becomes
<span class="chip redo_needed">redo needed</span>, automatically, with the
reason attached.</p>
{"".join(blocks)}
</div>
</section>"""


# ---------------------------------------------------------------------------
# (d) blockers
# ---------------------------------------------------------------------------
def section_blocking(project: pp.Project, ctx: Context) -> str:
    statuses = ctx.statuses(project.tasks)
    blocked = pp.open_blockers(project.tasks, statuses)
    if not blocked:
        body = ('<p class="sub">Nothing in this plan is blocked. Every '
                'remaining task is startable.</p>')
    else:
        body = "".join(
            f'<div class="blockcard"><b>{esc(t.title)}</b> '
            f'<span class="src">({esc(t.id)}, stage <code>{esc(t.stage)}'
            f'</code>)</span><br>{esc(t.blocker)}</div>'
            for t in blocked)
    return f"""
<section id="blocking">
<div class="bhead"><h2>What's blocking</h2>
<span class="tag">{len(blocked)} blocked task{"" if len(blocked) == 1 else "s"}</span></div>
<div class="stage">
<p class="sub">Each entry states the reason and what would clear it. A
blocker with no clearing condition is a decision nobody has made yet, and
it is written as one.</p>
{body}
</div>
</section>"""


# ---------------------------------------------------------------------------
# (e) next up
# ---------------------------------------------------------------------------
def _next_cards(tasks: Sequence[pp.Task],
                statuses: Mapping[str, str]) -> str:
    return "".join(
        f'<div class="nextcard">{_status_chip(statuses[t.id])} '
        f'<b>{esc(t.title)}</b> '
        f'<span class="src">({esc(t.phase)} &middot; {esc(t.id)} &middot; '
        f'stage <code>{esc(t.stage)}</code>)</span>'
        f'<br>{esc(t.produces)}<br>{_source_link(t.source)}</div>'
        for t in tasks)


def section_next(project: pp.Project, ctx: Context, page: Path) -> str:
    statuses = ctx.statuses(project.tasks)
    # Two questions, kept apart because they have different answers.  "What
    # must be restored" is a pipeline re-run; "what moves the paper forward"
    # is science. Collapsing them into one list buried the second under the
    # first every time an upstream stage went stale.
    restore = pp.next_up(project.tasks, statuses, limit=4,
                         include=(pp.REDO_NEEDED,))
    forward = pp.next_up(project.tasks, statuses, limit=3,
                         include=(pp.IN_PROGRESS, pp.PENDING))

    gated = pp.gated_tasks(project.tasks, statuses)

    parts = []
    if restore:
        parts.append(
            '<h3>First, restore the evidence</h3>'
            '<p class="sub">These tasks were done and are not backed any '
            'more — the stage their result rested on is no longer FRESH. '
            'One row per stage: re-running the stage clears every task that '
            'rests on it, so this is a list of RE-RUNS, not of re-decisions. '
            'The ordered commands are in <code>check_pipeline_status.py '
            'plan</code>.</p>'
            + _next_cards(restore, statuses))
    if forward:
        parts.append(
            '<h3>Then, the next steps forward</h3>'
            + _next_cards(forward, statuses))
    if not parts:
        parts.append(
            '<p class="sub">Nothing is actionable: every open task is '
            'blocked or gated. Read the sections above — the next move is '
            'to clear a blocker, not to start a task.</p>')

    if gated:
        cards = []
        for task, unmet in gated:
            names = ", ".join(
                f"<code>{esc(d)}</code> "
                f"({esc(pp.STATUS_LABEL.get(statuses.get(d, ''), '?'))})"
                for d in unmet)
            cards.append(
                f'<div class="blockcard">{_status_chip(statuses[task.id])} '
                f'<b>{esc(task.title)}</b> '
                f'<span class="src">({esc(task.phase)} &middot; '
                f'{esc(task.id)})</span><br>Waiting on {names}.'
                + (f'<br><span class="src">{esc(task.forbids)}</span>'
                   if task.forbids else "")
                + '</div>')
        parts.append(
            '<h3>Started or startable, but out of order</h3>'
            '<p class="sub">These are deliberately NOT recommended above. '
            'Their status says they are available; the execution order this '
            'plan was derived from says something else has to land first. '
            'That distinction used to be invisible — dependencies lived '
            'only inside blocker prose that nothing read, so the page '
            'offered production photometry as the actionable front while '
            'the detector work its own strategy puts in front of it sat '
            'blocked on destroyed tables.</p>' + "".join(cards))

    return f"""
<section id="next">
<div class="bhead"><h2>Next up</h2>
<span class="tag">the actionable front of the plan</span></div>
<div class="stage">
<p class="sub">In plan order, excluding blocked work.</p>
{"".join(parts)}
<p class="sub">When one completes:
<code>python pipeline/scripts/update_project_plan.py set &lt;task-id&gt; done
--evidence &lt;page-or-product&gt;</code>, then re-render. The page is never
edited by hand.</p>
</div>
</section>"""


# ---------------------------------------------------------------------------
# (f) evidence, with current provenance verdicts
# ---------------------------------------------------------------------------
def section_evidence(project: pp.Project, ctx: Context, page: Path) -> str:
    seen: list[str] = []
    for t in project.tasks:
        if t.stage not in seen:
            seen.append(t.stage)
    order = {k: i for i, k in enumerate(pv.topological_order(pv.STAGES))}
    seen.sort(key=lambda k: order.get(k, 999))

    rows, classes = [], []
    n_bad = 0
    for key in seen:
        stage = pv.STAGE_BY_KEY[key]
        verdict, why = ctx.verdict(key)
        if verdict != pv.FRESH:
            n_bad += 1
        rows.append([
            f"<code>{esc(key)}</code>",
            esc(stage.title),
            _verdict_chip(verdict),
            _report_link(page, key, ctx),
            esc(why) if why else "&mdash;",
        ])
        classes.append("" if verdict == pv.FRESH else "warn")

    return f"""
<section id="evidence">
<div class="bhead"><h2>Evidence</h2>
<span class="tag">{n_bad} of {len(seen)} stages are not FRESH</span></div>
<div class="stage">
<p class="sub">The pipeline stages this project's plan rests on, each with
the verdict <code>macro_core.provenance</code> returns right now.
<b>DESTROYED</b> is separated from <b>STALE</b> deliberately: stale means the
numbers moved, destroyed means the table the number came from is not in the
database at all. A reader is never invited to trust an orphaned report.</p>
{table(["Stage", "What it produces", "Verdict now", "Report", "Why"],
       rows, row_classes=classes)}
{_unbacked_note(seen, ctx)}
</div>
</section>"""


def _unbacked_note(stage_keys: Sequence[str], ctx: Context) -> str:
    """Name, in words, every constant this project would want to quote that
    currently has no table behind it.

    This is the correctness sweep made permanent: rather than a person
    remembering not to quote the High Gain ceiling, the page derives the
    warning from the same MISSING fingerprints that produce the verdicts.
    """
    dead: list[str] = []
    for key in stage_keys:
        stage = pv.STAGE_BY_KEY[key]
        for w in stage.writes:
            if w.startswith("table:") and ctx.fingerprints.get(w) == "MISSING":
                dead.append(w[len("table:"):])
    if not dead:
        return ""
    return (f'<div class="unbacked"><b>Constants with no query behind '
            f'them.</b> These tables are absent from the manifest right now: '
            f'<code>{esc(", ".join(sorted(set(dead))))}</code>. Every '
            f'detector number this project would otherwise quote — a '
            f'saturation ceiling, a veto threshold, a StackPro sub-read '
            f'count, a read noise — lived in them. Until the producing stage '
            f're-runs, those values are <b>unbacked prose</b>, not '
            f'measurements.<br><br><b>No page in this plan quotes them</b> — '
            f'and that is the limit of what this panel can promise. The '
            f'stage reports linked above are a different matter: they were '
            f'rendered while the tables still existed and still print those '
            f'constants as measured values, because nothing has re-rendered '
            f'them since the wipe. That is why every link to a report whose '
            f'stage is not FRESH now carries its verdict beside it. Treat an '
            f'orphaned report as a historical document, not a reference.'
            f'</div>')


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def _decisions(project: pp.Project) -> str:
    """Standing decisions this project's page must keep carrying.

    They come from the ledger, not from the previous version of this file:
    a generated page overwrites whatever prose was there, so any ruling that
    must survive has to live somewhere a render cannot destroy.
    """
    if not project.decisions:
        return ""
    blocks = "".join(
        f'<div class="decision"><b>{esc(head)}.</b> {esc(text)}</div>'
        for head, text in project.decisions)
    return (f'<h3>Standing decisions</h3><p class="sub">Rulings this project '
            f'is already operating under. They are held in the plan ledger '
            f'so a page regeneration cannot lose them.</p>{blocks}')


def _masthead_extra(counts: Mapping[str, int]) -> str:
    """The two counts a bare 'done/total' fraction would hide.

    A project whose evidence was destroyed reads as 0/32 — which is true but
    invites the wrong conclusion, that nothing was ever done. Saying "7 need
    redoing" in the same breath is the difference between "not started" and
    "not backed any more".
    """
    bits = []
    if counts.get(pp.REDO_NEEDED):
        bits.append(f"{counts[pp.REDO_NEEDED]} need redoing")
    if counts.get(pp.BLOCKED):
        bits.append(f"{counts[pp.BLOCKED]} blocked")
    return (" (" + ", ".join(bits) + ")") if bits else ""


def _page_path(project: pp.Project) -> Path:
    return DOCS_DIR / project.key / "index.html"


def render_project(project: pp.Project, ctx: Context) -> Path:
    page = _page_path(project)
    page.parent.mkdir(parents=True, exist_ok=True)
    statuses = ctx.statuses(project.tasks)
    counts = pp.status_counts(project.tasks, statuses)
    done, total = pp.progress_fraction(counts)

    strategy_link = (
        f'Strategy: <a href="https://github.com/jwwetzel/MACRO/blob/main/'
        f'{project.strategy}">ANALYSIS_STRATEGY.md</a> &middot; '
        if project.strategy else
        '<b>No committee strategy exists for this archive yet.</b> &middot; ')

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MACRO — {esc(project.title)}</title>
<link rel="stylesheet" href="../assets/macro.css">
<style>{PAGE_CSS}</style>
</head><body>

<header>
  <h1>{esc(project.title)} — Plan &amp; Progress</h1>
  <p>{done} of {total} plan tasks complete{_masthead_extra(counts)}
  &middot; {esc(project.venue)}<br>
  {strategy_link}<a href="../index.html">&larr; the front page</a></p>
</header>

<nav>
  <a href="#claim">Claim</a> &middot;
  <a href="#progress">Progress</a> &middot;
  <a href="#plan">The plan</a> &middot;
  <a href="#blocking">Blocking</a> &middot;
  <a href="#next">Next up</a> &middot;
  <a href="#evidence">Evidence</a>
</nav>

<section id="claim">
<div class="bhead"><h2>What this paper will claim</h2>
<span class="tag">{esc(project.venue)}</span></div>
<div class="stage">
<div class="decision">{esc(project.claim)}</div>
{_claim_measurement(project, ctx, _stage_table_of(project.key, ctx))}
{_decisions(project)}
</div>
</section>

{section_progress(project, ctx)}
{section_plan(project, ctx, page)}
{section_blocking(project, ctx)}
{section_next(project, ctx, page)}
{section_evidence(project, ctx, page)}

<footer>Generated {esc(ctx.generated)} by
<code>macro_core.report_projects</code> from the plan ledger in
<code>macro_core.project_plan</code> ({esc(pp.PLAN_CODE_VERSION)}), the
progress recorded in <code>{esc(pp.STATUS_TABLE)}</code>, and live queries
against <code>products/manifest/rlmt-manifest.sqlite</code>. Stage verdicts
come from <code>macro_core.provenance</code>
({esc(pv.PROVENANCE_CODE_VERSION)}). Regenerate with
<code>python pipeline/scripts/update_project_plan.py render</code>; record a
completed step with <code>… update_project_plan.py set &lt;task-id&gt;
done</code>.<br><br><b>Every number in the tables above is queried at render
time.</b> That is the claim this page can make, and the earlier and stronger
one — "no number on this page is typed by hand" — was false: task
descriptions quote figures from the strategy documents, and a quoted figure
can disagree with the table beside it when the two count under different
rules. Where a quoted number gates a decision, the accounting rule is stated
with it. Prose is prose; the tables are the measurement.</footer>
</body></html>"""
    page.write_text(html, encoding="utf-8")
    return page


# ---------------------------------------------------------------------------
# The hub
# ---------------------------------------------------------------------------
#: Stages listed on the hub, in build order, with the one-line scope the hub
#: has always shown.  The STATUS column is no longer typed: it is the live
#: provenance verdict, which is how the hub stopped claiming that eight
#: stages were "done" while two of them had no tables at all.
#: The scope strings COMPLEMENT each stage's own title (which is printed
#: beside them from the DAG) rather than restating it — they say what the
#: stage decides, which is what a reader scanning for "can I trust this?"
#: actually needs.
HUB_STAGES: tuple[tuple[str, str], ...] = (
    ("S0e", "rewrites the NAXIS values a raw-header parser got wrong; the "
            "event that invalidated S0"),
    ("S0", "global dedup, alias table, era tagging, pointing validation, "
           "night labels"),
    ("S0b", "raw↔reduced links, calibration census per era, coverage matrix, "
            "the October shopping list"),
    ("S0c", "per-project working sets, by reference — no project ever copies "
            "a frame"),
    ("S1", "the stratified go/no-go that set the batch-solve acceptance "
           "threshold"),
    ("S1b", "the production solve queue: which frames actually have a WCS"),
    ("S2", "ceiling, StackPro PTC, master reconstruction, linearity — every "
           "detector constant any paper would quote"),
    ("S2c", "measures source elongation per frame: is this FILTER a "
            "spectrum, or was that only its name?"),
    ("S3", "mid-exposure BJD_TDB, DATE-OBS audit, the clock bound, cadence"),
    ("S4", "the Honeycutt ensemble core and the empirical error model"),
    ("CV-S4", "the five staged CV targets, taken through the production "
              "photometry chain"),
    ("G", "grism identity gate, extraction, wavelength, response, EW"),
)


def render_hub(ctx: Context) -> Path:
    stage_rows, stage_classes = [], []
    for key, scope in HUB_STAGES:
        verdict, _why = ctx.verdict(key)
        stage_rows.append([
            f"<code>{esc(key)}</code> {esc(pv.STAGE_BY_KEY[key].title)}",
            esc(scope),
            _verdict_chip(verdict),
            _report_link(HUB_PATH, key, ctx),
        ])
        stage_classes.append("" if verdict == pv.FRESH else "warn")

    proj_rows = []
    for project in pp.PROJECTS:
        statuses = ctx.statuses(project.tasks)
        counts = pp.status_counts(project.tasks, statuses)
        done, total = pp.progress_fraction(counts)
        blocked = counts[pp.BLOCKED]
        nxt = pp.next_up(project.tasks, statuses, limit=1)
        # The hub bar carries the SAME five segments as a project page.  A
        # two-colour done/not-done bar would hide the distinction that
        # matters most here: work that was finished and then stopped being
        # backed is not the same as work never started.
        order = (pp.DONE, pp.IN_PROGRESS, pp.REDO_NEEDED, pp.BLOCKED,
                 pp.PENDING)
        bar = ('<div class="bar" style="width:150px;margin:4px 0">'
               + "".join(
                   f'<span class="{s}" style="width:'
                   f'{100.0 * counts[s] / total:.4f}%"></span>'
                   for s in order if counts[s])
               + "</div>")
        redo = counts[pp.REDO_NEEDED]
        proj_rows.append([
            f'<a href="{project.key}/index.html">{esc(project.title)}</a>',
            esc(project.venue),
            f"<b>{done}/{total}</b>"
            + (f' <span class="chip redo_needed">{redo} redo</span>'
               if redo else "") + bar,
            f"{blocked}" if blocked else "&mdash;",
            esc(nxt[0].title) if nxt else "&mdash;",
        ])

    n_fresh = sum(1 for k, _ in HUB_STAGES if ctx.verdict(k)[0] == pv.FRESH)
    all_tasks = pp.all_tasks()
    all_status = ctx.statuses(all_tasks)
    all_counts = pp.status_counts(all_tasks, all_status)
    a_done, a_total = pp.progress_fraction(all_counts)

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MACRO — the full evidence index</title>
<link rel="stylesheet" href="assets/macro.css">
<style>{PAGE_CSS}</style>
</head><body>

<header>
  <h1>The full evidence index</h1>
  <p>Every project and every pipeline stage, each with the verdict it carries
     right now. This is the audit view — the long way round. For the short
     one, start at <a href="index.html">the front page</a>.<br>
     Robert L. Mutel Telescope (0.5 m, Winer Observatory) &middot; 2023–2026
     archive: 3.34 TiB, ~329k FITS, 628 nights &middot; Every analytical
     decision on this site is backed by a plot and a script-emitted number.
     Code: <a href="https://github.com/jwwetzel/MACRO">github.com/jwwetzel/MACRO</a></p>
</header>

<nav>
  <a href="#projects">Projects</a> &middot;
  <a href="#pipeline">Pipeline</a>
</nav>

<section id="projects">
  <div class="bhead"><h2>Projects</h2>
    <span class="tag">{a_done} of {a_total} plan tasks complete across six projects</span></div>
  <div class="stage">
    <p class="sub">Each project page carries its whole plan — every phase,
    every task, what it produces, and where it stands — derived from that
    project's committee strategy and re-checked against the database on every
    render. The progress column is a count of plan tasks, not an opinion.</p>
{table(["Project", "Venue posture", "Plan progress", "Blocked", "Next up"],
       proj_rows)}
  </div>
</section>

<section id="pipeline">
  <div class="bhead"><h2>Shared pipeline</h2>
    <span class="tag">{n_fresh} of {len(HUB_STAGES)} stages are FRESH right now</span></div>
  <div class="stage">
    <p class="sub">Stages in build order. The verdict column is computed by
    <code>macro_core.provenance</code> at render time — it is not a
    hand-maintained "done" list, which is what let this table claim eight
    finished stages while two of them had no tables in the database at all.
    <b>DESTROYED</b> means the stage's declared outputs are absent, not merely
    stale.</p>
{table(["Stage", "Scope", "Verdict now", "Report"],
       stage_rows, row_classes=stage_classes)}
    <p class="sub">Full DAG, reasons and the ordered re-run plan:
    <code>python pipeline/scripts/check_pipeline_status.py status</code>.</p>
  </div>
</section>

<footer>Generated {esc(ctx.generated)} by
<code>macro_core.report_projects</code> &middot; project progress from
<code>macro_core.project_plan</code> ({esc(pp.PLAN_CODE_VERSION)}) and
<code>{esc(pp.STATUS_TABLE)}</code>; stage verdicts from
<code>macro_core.provenance</code> ({esc(pv.PROVENANCE_CODE_VERSION)}).
Reports regenerate with the pipeline; numbers are never edited by hand.</footer>
</body></html>"""
    HUB_PATH.write_text(html, encoding="utf-8")
    return HUB_PATH


def render_all(manifest_path: Path,
               projects: Optional[Sequence[str]] = None) -> list[Path]:
    """Render every project page and the hub.  Returns the paths written."""
    pp.validate()
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True,
                          timeout=300.0)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        ctx = Context(con, REPO_ROOT)
        written = []
        for project in pp.PROJECTS:
            if projects and project.key not in projects:
                continue
            written.append(render_project(project, ctx))
        written.append(render_hub(ctx))
        return written
    finally:
        con.close()
