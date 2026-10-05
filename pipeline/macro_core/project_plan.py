"""The PLAN LEDGER — every project's whole arc, as reviewable data.

WHY THIS MODULE EXISTS
----------------------
``provenance.py`` answers "is what we already did still true?".  It cannot
answer the question a person actually opens a project page to ask: **what is
the plan, how much of it is done, and what is the next thing to do?**

Before this module, that answer lived in three places that could not be
reconciled: each project's ``ANALYSIS_STRATEGY.md`` (the plan, in prose,
hundreds of lines long, with the execution order buried in §9), the pipeline
status check (the machinery, with no idea which project cares), and a
hand-written readiness paragraph on each ``docs/<Project>/index.html`` (a
snapshot of somebody's memory on the day they typed it).  A reader could
learn that S1 is stale and that a project "is ready for production analysis"
in the same sitting and never find out those two sentences contradict.

So this module makes the PLAN a DATA STRUCTURE, in the same spirit and the
same house pattern as :data:`macro_core.staging.PROJECT_SELECTIONS`:

* :data:`PROJECTS` declares, per project, its ordered :class:`Phase` list;
  each phase holds ordered :class:`Task` rows.
* Every task names the pipeline :class:`~macro_core.provenance.Stage` it
  depends on, so **staleness propagates into the plan**: a task marked
  ``done`` whose stage is no longer FRESH is not done any more, and
  :func:`derive_sync` says so with the reason attached.
* Every task cites the ANALYSIS_STRATEGY.md section it came from, so the
  ledger cannot quietly invent work the committee never planned — and a
  reviewer can diff it against the strategy in one pass.
* The code holds the PLAN.  The manifest's :data:`STATUS_TABLE` holds the
  PROGRESS, append-only with a UTC stamp, so a status history is replayable
  and "when did this become true?" has an answer.

WHAT IS PURE HERE
-----------------
Everything that decides anything: :func:`overlay_statuses`,
:func:`status_counts`, :func:`next_up`, :func:`open_blockers`,
:func:`derive_sync` and :func:`validate` see nothing but plain values, so
``pipeline/tests/test_project_plan.py`` drives them with hand-built
fixtures.  Only :func:`ensure_status_table`, :func:`record_status`,
:func:`read_statuses`, :func:`read_history` and :func:`stage_freshness`
touch a database or a disk.

THE STATUS RULE (the part that is easy to get wrong)
----------------------------------------------------
A ``done`` task is a claim about a database.  If the stage that produced its
evidence has since gone STALE, NEVER_RUN or lost its outputs, the claim is
no longer backed and the plan must say ``redo_needed`` — not ``done`` with a
small footnote.  :func:`derive_sync` therefore only ever moves a task
between ``done`` and ``redo_needed``; ``pending`` / ``in_progress`` /
``blocked`` are human judgements and no automatic rule may overwrite them.

HOW A TASK CLOSES  (plan review of 2026-10-03, SYNTHESIS §0)
------------------------------------------------------------
Before the review a task could leave the plan only by being done, so a task
the data could never support sat ``pending`` or ``blocked`` for ever — on
NGC 5548, behind a blocker ("the S1b batch reaches these frames") that could
not clear because the frames are spectra.  The committee ruled that a task
closes in exactly one of THREE ways, each recorded with its evidence:

* ``done``      — acceptance criterion met, evidence linked.
* ``dropped``   — ruled impossible or pointless on the data that exist.
* ``deferred``  — needs frames that do not exist yet; it moves to the
                  project's *2027 backlog*, outside the paper's critical path.

``dropped`` and ``deferred`` are not opinions a person may type: each must
carry a :class:`Ruling` — the committee document, the section and the
finding ids — and :func:`validate` refuses one without it.  Neither counts
toward completion in either direction: the headline fraction is
``done / in-scope`` (:func:`progress_fraction`), with the dropped count and
the backlog reported beside it (:func:`scope_summary`), so dropping work can
never make a project look more finished than the work that remains.

BLOCKERS ARE COMPUTED WHERE THEY CAN BE  (ruling U4)
----------------------------------------------------
Five tasks were blocked for six weeks on "S2's tables were destroyed" while
the tables sat in the database with rows in them.  A prose blocker goes
stale silently; a computed one cannot.  So wherever a gate is a fact about
the database it is DATA, evaluated when it is read:

* ``depends_on``   another task must be ``done`` — read from the status
                   table, and it may name a task in another group (the
                   shared foundation gates every project).
* ``needs_fresh``  a pipeline stage must be FRESH — read from the DAG.
* ``probes``       a read-only scalar query must come out a stated way.

:func:`computed_blockers` turns all three into sentences.  ``blocker`` prose
survives only for what no query can see: a closed dome, an unanswered
e-mail, a decision that is James's to make.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import re
import sqlite3
from dataclasses import dataclass, field, replace
from typing import Iterable, Mapping, Optional, Sequence

from . import provenance as pv

__all__ = [
    "PLAN_CODE_VERSION",
    "STATUS_TABLE",
    "DONE", "IN_PROGRESS", "BLOCKED", "PENDING", "REDO_NEEDED",
    "DROPPED", "DEFERRED", "CLOSED_BY_RULING", "BACKLOG_NAME",
    "LEDGER_STATUSES", "ALL_STATUSES", "OPEN_STATUSES", "STATUS_LABEL",
    "PlanError",
    "Source", "Ruling", "Probe", "ProbeResult", "Task", "Phase", "Project",
    "Change",
    "SYNTHESIS_2026_10_03", "COMMITTEE_SEATS", "RULING_ACTIONS",
    "AMENDMENTS_SECTION", "ruling_action",
    "PROJECTS", "PROJECT_BY_KEY", "FOUNDATION", "groups", "group_of",
    "all_tasks", "ledger_tasks", "tasks_of", "task_by_id", "project_of",
    "overlay_statuses", "status_counts", "progress_fraction",
    "scope_summary", "backlog", "dropped_tasks",
    "next_up", "open_blockers", "derive_sync", "validate",
    "stage_gates", "computed_blockers", "acceptance_gaps", "probe_is_met",
    "run_probes", "verify_rulings",
    "ruling_problems", "read_committee_documents", "read_status_stamps",
    "last_sync_utc", "manuscript_path", "manuscript_title", "with_title",
    "stage_report_path", "evidence_verdict", "stages_ever_run",
    "citation_problems", "verify_citations", "read_cited_documents",
    "split_numbered_sections", "unmet_dependencies", "gated_tasks",
    "ensure_status_table", "record_status", "read_statuses", "read_history",
    "stage_freshness", "utcnow",
]

#: Bumped when a change here would alter what a stored status row MEANS
#: (a new status value, a changed sync rule).  Recorded in the note of every
#: row this module writes, so a later reader can tell v1 progress from v2.
#:
#: v2.0 (2026-10-03, committee plan review): two closing statuses
#: (``dropped``, ``deferred``) that require a :class:`Ruling`; completion is
#: counted over IN-SCOPE tasks; a third sync rule lets a dated ruling
#: supersede a status recorded before it; computed gates (``needs_fresh``,
#: ``probes``, cross-group ``depends_on``); the shared-foundation group.
PLAN_CODE_VERSION = "PLAN v2.0 (2026-10-03)"

#: The one table this module owns inside the manifest database.
STATUS_TABLE = "project_plan_status"


class PlanError(RuntimeError):
    """Raised when the ledger cannot be read honestly — a duplicate task id,
    an unknown status, a stage that is not in the provenance DAG."""


# ===========================================================================
# 1.  STATUS VOCABULARY
# ===========================================================================

DONE = "done"
IN_PROGRESS = "in_progress"
BLOCKED = "blocked"
PENDING = "pending"

#: NOT a ledger value.  Only :func:`derive_sync` produces it, and only for a
#: task whose ``done`` claim rests on a stage that is no longer fresh.  It is
#: deliberately distinct from ``pending``: pending work was never done;
#: redo_needed work WAS done and stopped being true, which is a different
#: sentence to put in front of a reader.
REDO_NEEDED = "redo_needed"

#: CLOSED BY RULING — the committee decided the data that exist cannot
#: support this task, or that nobody needs its result (SYNTHESIS §0).  It
#: leaves the completion count altogether: not done, and no longer owed.
DROPPED = "dropped"

#: DEFERRED TO FUTURE DATA — the task needs frames that do not exist yet
#: (post-re-opening).  It moves to the project's 2027 backlog, which is
#: outside the paper's critical path and outside the completion count.
DEFERRED = "deferred"

#: The two statuses only a committee ruling can assign.  Both REQUIRE a
#: :class:`Ruling`; :func:`validate` and the ``set`` command enforce it.
CLOSED_BY_RULING: tuple[str, ...] = (DROPPED, DEFERRED)

#: What the deferred list is called wherever it is shown.
BACKLOG_NAME = "2027 backlog"

#: What a task may declare in code.
LEDGER_STATUSES: tuple[str, ...] = (DONE, IN_PROGRESS, BLOCKED, PENDING,
                                    DROPPED, DEFERRED)

#: Every value the status table may hold.
ALL_STATUSES: tuple[str, ...] = LEDGER_STATUSES + (REDO_NEEDED,)

#: Statuses that mean "there is work here that nobody has finished".
OPEN_STATUSES: tuple[str, ...] = (IN_PROGRESS, REDO_NEEDED, PENDING)

STATUS_LABEL: dict[str, str] = {
    DONE: "done",
    IN_PROGRESS: "in progress",
    BLOCKED: "blocked",
    PENDING: "pending",
    REDO_NEEDED: "redo needed",
    DROPPED: "dropped",
    DEFERRED: "2027 backlog",
}

#: The binding document of the 2026-10-03 plan review.  Every ``dropped`` /
#: ``deferred`` status and every task the review added or changed cites it.
SYNTHESIS_2026_10_03 = "committee/reviews/2026-10-03/SYNTHESIS.md"

#: Seat prefix -> the memo that seat filed, beside the synthesis.  A finding
#: id is ``<seat>.<id>`` (``DS.F1``, ``RF.B2``); a bare seat (``ED``) cites a
#: plan-hardening item its memo did not number.
COMMITTEE_SEATS: dict[str, str] = {
    "DS": "data-scientist.md",
    "OA": "observational-astronomer.md",
    "PH": "physicist.md",
    "DE": "detector-engineer.md",
    "TE": "telescope-engineer.md",
    "ED": "journal-editor.md",
    "RF": "hostile-referee.md",
}


# ===========================================================================
# 2.  THE LEDGER TYPES
# ===========================================================================

@dataclass(frozen=True)
class Source:
    """Where a task came from.

    ``document`` is a repo-relative path that must EXIST, and ``section``
    names the place inside it.  A task with no source is a task somebody
    invented, and the whole point of deriving the ledger from the five
    committee strategies is that nobody gets to do that silently.

    CHECKING THE DOCUMENT WAS NEVER ENOUGH.  The original validator asserted
    only ``task.source.document == project.strategy``, so any ``section``
    string whatsoever passed — and one did: T CrB's mid-exposure BJD_TDB
    task cited "§4 Phase A/B (BJD_TDB rule)" on a public page, hyperlinked
    to GitHub, against a document containing no BJD rule and, in its 213
    lines, not one occurrence of "BJD", "TDB", "barycentric" or
    "mid-exposure".  The work was real; the provenance was fiction, on the
    one page whose whole premise is that nobody invents work silently.
    :func:`citation_problems` closes that hole — see it for the rule.
    """

    document: str
    section: str

    def __str__(self) -> str:                      # pragma: no cover - trivial
        return f"{self.document} {self.section}"


@dataclass(frozen=True)
class Ruling:
    """A committee ruling that binds a task: where it is written, and why.

    ``section``   the place in the synthesis — ``"§4 TCrB_Monitoring"``,
                  ``"§3 F-4"``, ``"§1 U6"``.  Checked exactly like a
                  :class:`Source`: the ``§N`` must be a numbered section of
                  the document and every contentful word must appear in it.
    ``findings``  the ids behind the ruling — ``"U6"`` / ``"D1"`` for a
                  synthesis ruling, ``"DS.F1"`` for a seat's finding,
                  ``"ED"`` for an un-numbered plan-hardening item in that
                  seat's memo.  :func:`ruling_problems` resolves each one
                  against the memo it names, so a ruling cannot cite a
                  finding nobody filed.
    ``date``      the day it was made.  It is what lets ``sync`` tell a
                  status recorded BEFORE the ruling (superseded by it) from
                  one recorded after (a later human judgement, left alone).
    ``action``    the committee's verb, ONLY where it cannot be read off the
                  task (see :func:`ruling_action`): a HOLD looks exactly
                  like any other blocked task, and a task added to a project
                  with no strategy cites the same decision a changed one
                  does.  Left empty everywhere else, so the verb cannot
                  disagree with the status it describes.
    """

    section: str
    findings: tuple[str, ...] = ()
    document: str = SYNTHESIS_2026_10_03
    date: str = "2026-10-03"
    action: str = ""

    @property
    def source(self) -> Source:
        return Source(self.document, self.section)

    def __str__(self) -> str:
        cited = f" ({', '.join(self.findings)})" if self.findings else ""
        return f"{self.document} {self.section}{cited}, {self.date}"


@dataclass(frozen=True)
class Probe:
    """A gate that is a fact about the database, asked of the database.

    ``sql`` is ONE read-only query returning ONE integer.  ``met_when`` says
    which answer satisfies it: ``"zero"`` (nothing left of the defect) or
    ``"positive"`` (the thing now exists).  A probe whose table is absent is
    unmet and says so — an absent table is not evidence of anything.

    ``kind`` says which side of the task the fact sits on:

    * ``"gate"``    it must hold BEFORE the task may start.  Unmet, it is a
                    computed blocker (:func:`computed_blockers`).
    * ``"accept"``  it holds once the task is DONE — the task's acceptance
                    criterion, where that criterion is a number the manifest
                    can produce.  Unmet, it blocks nothing (a task is not
                    blocked by its own incompleteness); but a task recorded
                    ``done`` with an unmet acceptance probe is a
                    contradiction, and :func:`acceptance_gaps` reports it.

    This is the cure for the blocker that read "S2's tables were destroyed"
    for six weeks after they were rebuilt: nobody has to remember to update
    a sentence when the sentence is a query.
    """

    label: str
    sql: str
    met_when: str = "zero"
    kind: str = "gate"


@dataclass(frozen=True)
class ProbeResult:
    """One evaluated :class:`Probe`: the number, and whether it opens."""

    label: str
    value: Optional[int]
    met: bool
    error: str = ""
    kind: str = "gate"

    def __str__(self) -> str:
        if self.error:
            return f"{self.label}: {self.error}"
        return f"{self.label}: {self.value:,}"


@dataclass(frozen=True)
class Task:
    """One unit of the plan.

    ``id``        stable, human-typable, never reused (it is the key in the
                  status table, so renaming one orphans its history).
    ``title``     what a person would call the work.
    ``produces``  ONE line: the artifact that exists afterwards.  A task
                  whose product cannot be named in one line is two tasks.
    ``stage``     the ``provenance.STAGES`` key this task's result rests on.
                  This is the edge that makes staleness propagate into the
                  plan; a task with the wrong stage is a task that will lie.
    ``source``    the strategy section it was derived from.
    ``status``    the LEDGER status (see :data:`LEDGER_STATUSES`).  The
                  manifest's status table overrides it once work happens.
    ``evidence``  repo-relative page or product path, for ``done`` tasks.
    ``blocker``   why it cannot start, for ``blocked`` tasks — and what
                  would clear it.  Required when status is ``blocked``.
    ``depends_on`` task ids that must be ``done`` before this one may be
                  RECOMMENDED.  Before this field existed, dependencies
                  lived only as prose inside ``blocker`` strings, which
                  nothing read — so :func:`next_up` ranked on status alone
                  and cheerfully offered "production photometry on ST LMi"
                  as the actionable front while the two detector tasks its
                  own strategy puts in front of it sat blocked.  A plan that
                  cannot represent "after" will recommend work its execution
                  order forbids.
    ``forbids``   the sentence in the strategy that PROHIBITS running this
                  task before its dependencies land, quoted.  Shown next to
                  the task wherever it is offered, so a reader who ignores
                  the gate at least does so knowingly.
    ``ruling``    the committee ruling that added, changed or closed this
                  task.  MANDATORY for ``dropped`` and ``deferred`` — those
                  two statuses are rulings, not opinions.
    ``accept``    ONE line: the test that closes the task.  Mandatory for
                  every open task that carries a ruling, because the review
                  found tasks marked done against a description rather than
                  a criterion ("measured linearity curves" for ladders that
                  never reached a veto).  :func:`_build` appends it to
                  ``produces`` so it shows wherever the task does.
    ``needs_fresh`` pipeline stages that must be FRESH before this task may
                  be recommended — a gate read from the DAG, not typed.
    ``probes``    database facts that must hold first (see :class:`Probe`).
    ``project`` / ``phase`` are filled in by :func:`_build` from the
    containing structures, so they can never disagree with it.
    """

    id: str
    title: str
    produces: str
    stage: str
    source: Source
    status: str = PENDING
    evidence: str = ""
    blocker: str = ""
    depends_on: tuple[str, ...] = ()
    forbids: str = ""
    ruling: Optional[Ruling] = None
    accept: str = ""
    needs_fresh: tuple[str, ...] = ()
    probes: tuple[Probe, ...] = ()
    project: str = ""
    phase: str = ""


@dataclass(frozen=True)
class Phase:
    """An ordered group of tasks with one intent."""

    name: str
    intent: str
    tasks: tuple[Task, ...]


@dataclass(frozen=True)
class Project:
    """One paper (or one candidate for one).

    ``claim``   the paragraph the page opens with: what the paper will claim
                and where it will go.  It is here, in the ledger, rather than
                in the renderer, because it is a decision — and decisions
                belong with the plan they justify.
    ``strategy`` repo-relative path of the governing document, or ``""`` for
                a project that has none yet (Legacy_Rigel).
    ``ruling``  the committee ruling that re-scoped the project, if one did.
                Carried on the project because a ruling that changes a title
                or a venue is not a task and must still be citable.
    ``paper_title`` the manuscript's ``\\title{...}``, in LaTeX, or ``""`` for
                a project whose manuscript this ledger does not own (the CV
                draft is edited by hand) or that has none.  It lives here
                because a title is a CLAIM: three skeletons went on
                advertising papers their own strategies had killed — a
                "Three-Year Quiescent Baseline", "Early" SN photometry,
                "Photometric Monitoring of NGC 5548" — because the title was
                typed once into a .tex file nothing re-read.
                ``update_project_plan.py titles`` compares and re-emits it.
    ``decisions`` standing decisions this project's page must keep carrying
                — ``(heading, text)`` pairs.  They live HERE, in reviewed and
                diffable code, precisely because the page is now generated:
                a decision that existed only as prose on a hand-written page
                would be destroyed by the first render.  Legacy_Rigel's
                separate-archive-root ruling is the case in point.
    ``claim_filters`` filter slots the CLAIM depends on.  A claim that rests
                on "the slitless grism series in slot '6'" is a claim about
                a measurement, and S2c has measured it; naming the slot here
                makes the renderer print the live verdict directly beneath
                the claim.  The alternative — a frame count typed into the
                claim prose — is what let the SN page assert an 83-frame
                grism series in its opening paragraph while its own filter
                table, a hundred lines below, measured 3 of those frames as
                direct imaging and 19 as undecided, and never connected the
                two sentences.
    """

    key: str
    title: str
    claim: str
    venue: str
    strategy: str
    phases: tuple[Phase, ...]
    decisions: tuple[tuple[str, str], ...] = ()
    claim_filters: tuple[str, ...] = ()
    ruling: Optional[Ruling] = None
    paper_title: str = ""

    @property
    def tasks(self) -> tuple[Task, ...]:
        return tuple(t for ph in self.phases for t in ph.tasks)


@dataclass(frozen=True)
class Change:
    """One status transition proposed by :func:`derive_sync`."""

    task_id: str
    old: str
    new: str
    reason: str


def _build(project: Project) -> Project:
    """Stamp every task with its project and phase.

    Done once at import so the two fields cannot drift from the structure
    that contains them — the alternative (typing ``project=`` on 131 tasks)
    is a typo waiting to mis-file a task onto the wrong page.

    The acceptance criterion is appended to ``produces`` here, once, so
    every view that prints what a task produces also prints the test that
    closes it.  It is composed rather than typed twice: ``accept`` stays the
    single source, and a renderer that later grows its own column for it
    only has to stop reading the composed line.
    """
    def stamp(t: Task, phase: str) -> Task:
        produces = t.produces
        if t.accept and _ACCEPT_MARK not in produces:
            produces = f"{produces.rstrip()} {_ACCEPT_MARK} {t.accept}"
        return replace(t, project=project.key, phase=phase, produces=produces)

    phases = tuple(
        replace(ph, tasks=tuple(stamp(t, ph.name) for t in ph.tasks))
        for ph in project.phases)
    return replace(project, phases=phases)


#: How an acceptance criterion is introduced inside ``produces``.
_ACCEPT_MARK = "Accept:"

#: The two phases every project files its closed-by-ruling tasks under, so a
#: page shows them apart from live work without the renderer having to know
#: the statuses exist.  :func:`validate` holds the ledger to it.
DROPPED_PHASE = "Dropped by committee ruling"
BACKLOG_PHASE = f"{BACKLOG_NAME} — deferred to future data"


def _dropped_phase(*tasks: Task) -> Phase:
    """The phase a project's ``dropped`` tasks are filed under."""
    return Phase(
        DROPPED_PHASE,
        "Ruled impossible or pointless on the data that exist. Each row "
        "names its ruling. These are outside the completion count: not "
        "done, and no longer owed.",
        tuple(tasks))


def _backlog_phase(*tasks: Task) -> Phase:
    """The phase a project's ``deferred`` tasks are filed under."""
    return Phase(
        BACKLOG_PHASE,
        "Needs frames that do not exist yet. Off this paper's critical "
        "path and outside its completion count; picked up when the "
        "post-re-opening data exist.",
        tuple(tasks))


def _ruled(section: str, *findings: str, action: str = "") -> Ruling:
    """A ruling of the 2026-10-03 plan review."""
    return Ruling(section, tuple(findings), action=action)


#: The committee's verbs (SYNTHESIS §4 uses exactly these).
RULING_ACTIONS: tuple[str, ...] = (
    "ADD", "CHANGE", "DROP", "DEFER", "HOLD", "CLOSE", "UNBLOCK")

#: The numbered section of each strategy that records the rulings binding
#: it.  A task whose SOURCE is this section did not exist before the review.
AMENDMENTS_SECTION = "§10"


def ruling_action(task: Task) -> str:
    """What the committee did to this task, in the committee's own verb.

    Derived wherever the task already says it, so the verb cannot drift
    from the facts it summarises:

    * ``dropped`` / ``deferred``            -> DROP / DEFER
    * source is the amendments section      -> ADD   (it is new)
    * ``done``                              -> CLOSE (closed by the ruling)
    * the ruling's section says UNBLOCK     -> UNBLOCK
    * otherwise                             -> CHANGE

    An explicit :attr:`Ruling.action` wins over all of that; it exists for
    the two cases the derivation cannot see (HOLD; ADD in a project that
    has no strategy to cite).  A task with no ruling has no verb.
    """
    if task.ruling is None:
        return ""
    if task.ruling.action:
        return task.ruling.action
    if task.status == DROPPED:
        return "DROP"
    if task.status == DEFERRED:
        return "DEFER"
    if task.source.section.startswith(AMENDMENTS_SECTION):
        return "ADD"
    if task.status == DONE:
        return "CLOSE"
    if "UNBLOCK" in task.ruling.section:
        return "UNBLOCK"
    return "CHANGE"


# ===========================================================================
# 3.  THE PLAN — six projects and the shared foundation
# ===========================================================================
# Every task below was read out of the cited document's execution order, or
# added by the plan review of 2026-10-03 (committee/reviews/2026-10-03/).
# Where a phase has already been executed it is marked done and cites the
# evidence page.  Where the review changed, closed or added a task it
# carries the `Ruling` that did so; a task the review dropped or deferred is
# filed under its project's "Dropped by committee ruling" or "2027 backlog"
# phase and is outside the completion count.
#
# NOTHING BELOW IS BLOCKED ON "S2's tables were destroyed" ANY MORE (ruling
# U4).  That sentence was true for one day in August and sat in five
# blockers for six weeks after S2 was rebuilt.  The tasks it held are
# unblocked; what really gates each of them is stated as a dependency on a
# shared-foundation task, which is read from the status table every time it
# is shown, and the tables themselves are probed.

#: The section of each strategy that records the 2026-10-03 rulings.  Tasks
#: the review ADDED cite it (it lists every one by id, with the finding
#: behind it), so they obey the same rule as every other task: the source is
#: the project's own strategy, and the citation checker can resolve it.
_AMENDMENTS = AMENDMENTS_SECTION


#: Ruling U4: these S2 tables were named as "destroyed" by blockers that
#: outlived the rebuild by six weeks.  The probe is the blocker's
#: replacement — it asks the database instead of remembering.
def _s2_table_probe(table: str) -> Probe:
    return Probe(f"{table} rows (U4: the table a stale blocker called "
                 f"destroyed)", f"SELECT count(*) FROM {table}", "positive")



_CV = "CV_TimeSeries/ANALYSIS_STRATEGY.md"

#: Phase 2 production photometry may not run before these three land.  The
#: edge is real and written down; before it was DATA rather than prose, the
#: page recommended ST LMi photometry as the actionable front while all
#: three sat blocked.
_CV_DETECTOR_GATE = ("CV-P15-linearity-ladders", "CV-P15-noise-model",
                     "CV-P2-vetoes")

_CV_GATE_RULE = (
    "\u00a75 row 1 is explicit: until the linearity ladder exists, N_sub and the "
    "StackPro variance model are hypotheses and NO mixed-mode fit is legal. "
    "\u00a79's execution order puts the ladders and the noise model between Phase 1 "
    "and Phase 2 photometry, and step 10's per-mode saturation vetoes gate "
    "every frame this task would use.")


def _cv_src(section: str) -> Source:
    return Source(_CV, section)


CV_TIMESERIES = Project(
    key="CV_TimeSeries",
    title="Cataclysmic-Variable Time Series (reopened)",
    claim=(
        "REOPENED 2026-10-03 — major revision; not complete. What the data "
        "deliver is one polar's colour–phase curve, a timing null and a "
        "release: ST LMi's g−r swings strongly through the bright phase "
        "with the same morphology in two instrument eras a year apart, "
        "state-tagged against the survey record; bright-phase edge timing "
        "bounds the change of the spin/spot-longitude period (not the "
        "orbital one); YZ Cnc contributes an excluded-amplitude superhump "
        "statement; VV Pup, EU UMa and AN UMa are a coverage audit. The "
        "draft's sentence that the bright-phase edge 'shows no offset' "
        "between bands is contradicted by its own released edges: g leads "
        "i on nearly every night both bands were timed. Whether that is "
        "astrophysical or an estimator bias is NOT decided (D3) — it is "
        "tested by per-band injection, and the paper reports whichever the "
        "test supports. Fifteen tasks (CV-R1 … CV-R15) block release."),
    venue="AJ — major revision open (reopened 2026-10-03; the venue was "
          "ApJ)",
    strategy=_CV,
    ruling=_ruled("§4 CV_TimeSeries", "U3", "D3"),
    decisions=(
        ("Plan review of 2026-10-03 — the paper is not finished",
         "Major revision (U3): the plan read 34 of 34 and the committee "
         "found four blockers in the abstract and conclusions. (1) The "
         "band-offset 'null' is a product of the error model — "
         "max(χ²ν, 1) on a transported budget doubled the error bars — "
         "and is re-tested with paired, scatter-based statistics. (2) The "
         "'bias floor' is replaced by the signed matched-cell bias per "
         "band and era, and D3 (astrophysics or estimator bias?) is "
         "settled by injection, not by vote. (3) The O−C is refitted with "
         "per-band constants and night-level epochs. (4) The superhump "
         "'measurement of absence' is restated as the amplitude range "
         "actually excluded. (5) The colour–phase result, which the "
         "abstract never quantified, becomes the lead. (6) The paper meets "
         "the literature, is cut to ≤ 18 pages with a ≤ 250-word abstract, "
         "and goes to AJ. (7) Nothing is citable until the chain is "
         "rebuilt FRESH at a clean commit and the macro diff is filed. "
         "(8) Authorship, ORCIDs, a DOI and an outside reader need James."),
    ),
    phases=(
        Phase("Phase 0 — Curation",
              "Blocks everything else: one canonical, alias-merged, "
              "tree-pinned view per target, and the AAVSO call that decides "
              "whether Q3 exists.",
              (
                  Task("CV-P0-curation-sql",
                       "Canonical per-target frame view",
                       "One alias-merged, rawimage-pinned, photometric-filter "
                       "frame list per target, staged as stage_cv_timeseries.",
                       "S0c", _cv_src("§4 Phase 0 step 1"), DONE,
                       evidence="docs/pipeline/s0c_staging.html"),
                  Task("CV-P0-aavso-yzcnc",
                       "AAVSO cross-match for YZ Cnc, 2024-02-21 → 2024-05-03",
                       "An outburst-state tag per dense run — the branch point "
                       "that decides whether Q3 is superhumps or flickering.",
                       "CV-S7", _cv_src("§4 Phase 0 step 2"), DONE,
                       evidence="docs/CV_TimeSeries/cv_external_context.html"),
                  Task("CV-P0-survey-context",
                       "Pull the survey record for all five targets",
                       "Cached ZTF/ATLAS/ASAS-SN/AAVSO/TESS/Gaia/eROSITA "
                       "series that the state histories are plotted over.",
                       "CV-S7", _cv_src("§4 Phase 0 step 3"), DONE,
                       evidence="docs/CV_TimeSeries/cv_external_context.html"),
              )),
        Phase("Phase 0.5 — Astrometry go/no-go",
              "The pipeline's first bottleneck: ~4,600 of ~5,500 polar "
              "Sloan-era frames had no WCS, so no ensemble photometry.",
              (
                  Task("CV-P05-solve-experiment",
                       "Stratified 200-frame re-solve experiment",
                       "A measured per-target success rate, and the "
                       "acceptance threshold set FROM it rather than "
                       "asserted.",
                       "S1", _cv_src("§4 Phase 0.5 step 3a"), DONE,
                       evidence="docs/pipeline/s1_astrometry.html"),
                  Task("CV-P05-batch-solve",
                       "Production batch solve of the unsolved pool",
                       "A WCS for every solvable CV frame; the per-target "
                       "solved fraction that gates Phase 2.",
                       "S1b", _cv_src("§4 Phase 0.5 step 3b"), IN_PROGRESS,
                       evidence="docs/pipeline/s1_astrometry.html"),
                  Task("CV-P05-geometry-requeue",
                       "Re-queue the frames excluded by the NAXIS artifact",
                       "The 18,381 frames whose recorded geometry was a "
                       "tile-compressed BINTABLE row length — EU UMa's 207 "
                       "among them — returned to the solve queue.",
                       "S0e", _cv_src("§5 row 'Astrometry'"), IN_PROGRESS,
                       evidence="docs/pipeline/s0e_geometry_fix.html"),
              )),
        Phase("Phase 1 — Timing foundation",
              "Mid-exposure BJD_TDB from scratch; header JD is UTC "
              "exposure start and is never used.",
              (
                  Task("CV-P1-bjd",
                       "Mid-exposure BJD_TDB for every frame",
                       "frame_times: barycentric mid-exposure timestamps "
                       "with a JPL ephemeris at Winer's EarthLocation.",
                       "S3", _cv_src("§4 Phase 1 step 4"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("CV-P1-era-audit",
                       "DATE-OBS start-vs-mid convention audit, per era",
                       "s3_dateobs_audit: the convention verified "
                       "independently in MaxIm and pyscope frames.",
                       "S3", _cv_src("§4 Phase 1 step 5"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("CV-P1-clock-validation",
                       "Observatory clock validation on an archived EB",
                       "A published clock bound: |Δt| < 4,517 s (~75 min) "
                       "once the ephemeris period uncertainty is carried.",
                       "S3", _cv_src("§4 Phase 1 step 6"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("CV-P1-stackpro-midtime",
                       "StackPro mid-time semantics",
                       "The worst-case StackPro mid-time error (5.67 s) and "
                       "the frame list sub-second timing must avoid.",
                       "S3", _cv_src("§6 failure mode 7"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
              )),
        Phase("Phase 1.5 — Detector truth",
              "§5 rows 1–2: until the ladders exist, N_sub and the StackPro "
              "variance model are hypotheses and no mixed-mode fit is legal.",
              (
                  # These three carried blockers reading "S2's evidence
                  # tables were destroyed" (ruling U4).  S2 re-ran on
                  # 2026-08-20 and all three were recorded done against its
                  # page; the ledger now says what the status table has
                  # said since.  What the review found is not that the
                  # tables are missing but that they do not measure enough
                  # — DE.F2, DE.F3, DS.F7 — and that is CV-R10 and the
                  # shared foundation's F-4/F-5, not a blocker here.
                  Task("CV-P15-linearity-ladders",
                       "Linearity ladder per readout mode",
                       "The best archival linearity ladder per readout "
                       "mode and the veto thresholds in force. NOT a "
                       "linearity measurement near any veto: the Mode0 "
                       "ladder peaks at ~1% of full scale. The "
                       "residual-vs-peak test that would be one is F-5, "
                       "reported in the paper by CV-R10.",
                       "S2", _cv_src("§5 row 'Linearity, per readout mode'"),
                       DONE, evidence="docs/pipeline/s2_detector.html",
                       probes=(_s2_table_probe("s2_linearity_ladders"),)),
                  Task("CV-P15-noise-model",
                       "Empirical noise model per readout mode",
                       "Check-star RMS vs magnitude per mode per night, with "
                       "the model overplotted — the non-negotiable figure. "
                       "It ran on a nominal gain with a bracket; the "
                       "measured gains replace the bracket in CV-R10.",
                       "S2", _cv_src("§5 row 'Noise model, per mode'"),
                       DONE, evidence="docs/pipeline/s2_detector.html",
                       probes=(_s2_table_probe("s2_ptc_fits"),)),
              )),
        Phase("Phase 2 — Photometry",
              "Ensemble differential photometry tied to ATLAS-REFCAT2, per "
              "target, in the order §9 sets: ST LMi first (richest), then "
              "VV Pup (hardest), EU UMa, YZ Cnc.",
              (
                  # Stage re-pointed S4 -> CV-S4 (2026-08-20).  S4 was the
                  # two-target prototype that first proved the solver; the
                  # production run re-solved the same core over all 8,716
                  # staged frames and carries the synthetic-recovery test in
                  # its own suite.  Leaving the task pinned to the prototype
                  # made it perpetually redo_needed for a stage nothing
                  # downstream reads any more — the evidence moved, so the
                  # citation follows it.
                  Task("CV-P2-ensemble-core",
                       "Honeycutt ensemble core, proven on a prototype",
                       "A validated inhomogeneous-ensemble solver: synthetic "
                       "zero-point pattern recovered to <6 mmag with check "
                       "stars held out.",
                       "CV-S4", _cv_src("§4 Phase 2 steps 7–9"), DONE,
                       evidence="docs/pipeline/s4_photometry.html"),
                  Task("CV-P2-vetoes",
                       "Per-mode saturation vetoes",
                       "A peak-ADU veto per readout mode applied to every "
                       "frame, with veto counts reported per target.",
                       "S2", _cv_src("§4 Phase 2 step 10"), DONE,
                       evidence="docs/pipeline/s2_detector.html",
                       probes=(_s2_table_probe("s2_ceiling_modes"),)),
                  Task("CV-P2-stlmi", "Production photometry: ST LMi",
                       "Calibrated per-frame light curves for the flagship "
                       "target, both filter eras kept separate.",
                       "CV-S4", _cv_src("§9 execution order"), IN_PROGRESS,
                       depends_on=_CV_DETECTOR_GATE,
                       forbids=_CV_GATE_RULE),
                  Task("CV-P2-vvpup", "Production photometry: VV Pup",
                       "Per-camera calibrated light curves (iKon and Mode0 "
                       "never jointly solved — camera and epoch are "
                       "confounded).",
                       "CV-S4", _cv_src("§4 step 13a"), PENDING,
                       depends_on=_CV_DETECTOR_GATE,
                       forbids=_CV_GATE_RULE),
                  Task("CV-P2-euuma", "Production photometry: EU UMa",
                       "Calibrated light curves with the 4.4% phase smear of "
                       "the dominant 240 s mode propagated.",
                       "CV-S4", _cv_src("§9 execution order"), PENDING,
                       depends_on=_CV_DETECTOR_GATE,
                       forbids=_CV_GATE_RULE),
                  Task("CV-P2-yzcnc", "Production photometry: YZ Cnc",
                       "Calibrated light curves for the dense 2024 blocks, "
                       "every outburst-night frame saturation-audited.",
                       "CV-S4", _cv_src("§9 execution order"), PENDING,
                       depends_on=_CV_DETECTOR_GATE,
                       forbids=_CV_GATE_RULE),
                  Task("CV-P2-cloud-veto",
                       "Ensemble-flux-ratio cloud veto",
                       "A primary frame-quality veto that works without "
                       "zmag — which does not exist for the polar Sloan era.",
                       "CV-S4", _cv_src("§4 Phase 2 step 11"), PENDING),
                  Task("CV-P2-extinction",
                       "Second-order colour-extinction terms",
                       "k″·(g−r)·X solved inside the ensemble, per camera "
                       "for VV Pup (X ≥ 1.57 always).",
                       "CV-S4", _cv_src("§4 Phase 2 step 12"), PENDING),
                  # Added 2026-08-19.  The ledger was derived before this
                  # work existed, but the strategy has always carried the
                  # goal: §4 step 9 ends "Tie the ensemble zero point to
                  # REFCAT2 g,r,i per night -> PS1 AB system to ~0.01-0.02
                  # mag absolute", and §5's table row "Absolute calibration"
                  # repeats it.  Without a task, completed work had nowhere
                  # to be recorded and the page under-reported the project.
                  Task("CV-P2-cattie",
                       "Tie the ensembles to ATLAS-REFCAT2",
                       "Zero point and colour term per (target, era, filter) "
                       "solved on comparison stars, validated on held-out "
                       "stars; natural-system magnitudes published, the CV "
                       "itself never transformed.",
                       "CV-S6", _cv_src("§4 Phase 2 step 9"), PENDING),
                  Task("CV-P2-cross-era",
                       "Cross-era discipline and transformation metadata",
                       "G/R/I → g/r/i coefficients ±σ published as data-"
                       "release metadata; the CV itself never transformed.",
                       "CV-S4", _cv_src("§4 Phase 2 step 13"), PENDING),
                  Task("CV-P2-faint-limits",
                       "Faint-phase forced photometry and upper limits",
                       "Uncensored state statistics: forced photometry at "
                       "the solved position, limits for non-detections.",
                       "CV-S4", _cv_src("§4 Phase 2 step 14"), PENDING),
              )),
        Phase("Phase 3 — Time-series analysis",
              "Confirmation and alias hygiene, not discovery: three "
              "independent period methods must agree, and every threshold is "
              "demonstrated by injection before it is adopted.",
              (
                  Task("CV-P3-periods", "Period verification, per filter per era",
                       "LS + PDM + conditional entropy agreeing on one alias "
                       "family, with the published spectral window per season.",
                       "CV-S4", _cv_src("§4 Phase 3 step 15"), PENDING),
                  Task("CV-P3-sigma-t",
                       # The night is 2025-02-27 by the roadmap's
                       # local-noon convention, which is what the database
                       # and the paper carry; the strategy's "2025-02-28"
                       # is the same night by its UTC date (RF minor 5).
                       "σ_t injection test on ST LMi (night 2025-02-27)",
                       "The demonstrated per-cycle timing precision that "
                       "decides whether the per-cycle O–C tier exists at all.",
                       "CV-S4", _cv_src("§4 Phase 3 step 16"), PENDING),
                  Task("CV-P3-bright-phase",
                       "Bright-phase timing, per band",
                       "Per-cycle ingress/egress/centroid epochs with emcee "
                       "uncertainties; band-dependent egress as a result.",
                       "CV-S4", _cv_src("§4 Phase 3 step 16"), PENDING),
                  Task("CV-P3-oc", "O–C construction and cycle-count analysis",
                       "Seasonal O–C against the 40-yr baseline plus an "
                       "explicit ambiguity analysis across the 289-d gap.",
                       "CV-S4", _cv_src("§4 Phase 3 step 17"), PENDING),
                  Task("CV-P3-states", "Accretion-state classification",
                       "A two-component mixture model on nightly means — "
                       "state boundaries from the model, never by eye.",
                       "CV-S4", _cv_src("§4 Phase 3 step 18"), PENDING),
                  Task("CV-P3-yzcnc-superhump",
                       "YZ Cnc superhump analysis (or the fallback)",
                       "P_sh and a Kato-style O–C per filter if the dense "
                       "runs were in outburst; orbital hump + flickering "
                       "statistics if not.",
                       "CV-S10", _cv_src("§4 Phase 3 step 19"), DONE,
                       evidence="docs/CV_TimeSeries/cv_final_science.html"),
                  Task("CV-P3-detrending", "Detrending discipline",
                       "Per-night systematics fit JOINTLY with the periodic "
                       "model — low-order airmass polynomial, or a celerite2 "
                       "Matérn-3/2 GP whose timescale prior is bounded below "
                       "at 3× the candidate period — plus the with/without "
                       "comparison figures. Never pre-whiten a night "
                       "spanning under three cycles with a free smooth "
                       "trend: EU UMa's nights average ~1.5 cycles and a "
                       "free spline eats the orbit.",
                       "CV-S4", _cv_src("§4 Phase 3 step 20"), PENDING),
                  Task("CV-P3-injection",
                       "Detection limits and injection–recovery",
                       "90%-recovery amplitude–period contours computed at "
                       "the real timestamps, per target per filter.",
                       "CV-S4", _cv_src("§4 Phase 3 step 21"), PENDING),
              )),
        Phase("Phase 4 — Decisions, figures, draft",
              "The conditional target call, then the figure set, then the "
              "manuscript.",
              (
                  Task("CV-P4-anuma", "AN UMa go/no-go, per filter",
                       "A decision: colour analysis only if ≥8 full-orbit "
                       "three-filter nights survive curation (currently ~7).",
                       "CV-S10", _cv_src("§2 Q5"), DONE,
                       evidence="docs/CV_TimeSeries/cv_final_science.html"),
                  Task("CV-P4-figures", "The 13-figure set",
                       "Every figure of §7, each regenerable from the "
                       "photometry product.",
                       "CV-S4", _cv_src("§7 Figure list"), PENDING),
                  Task("CV-P4-draft", "Manuscript draft",
                       "manuscripts/CV_TimeSeries/main.tex filled out "
                       "against the §8 outline.",
                       "CV-S4", _cv_src("§8 Manuscript outline"), PENDING),
              )),
        Phase("Phase 5 — Major revision (plan review of 2026-10-03)",
              "All fifteen block release. The draft's numbers reproduce; "
              "four of the inferences drawn from them do not survive, and "
              "the paper has not met the literature.",
              (
                  Task("CV-R1-band-offset",
                       "Re-test the band offset with paired statistics",
                       "A paired, distribution-free test of the g−i (and "
                       "g−r, r−i) bright-phase edge offset on same-cycle "
                       "and same-night edges, a combined-era estimate, and "
                       "scatter-based errors beside the budget-based ones.",
                       "CV-S9", _cv_src(f"{_AMENDMENTS} CV-R1-band-offset"),
                       PENDING,
                       ruling=_ruled("§4 CV-R1-band-offset",
                                     "U3", "DS.F1", "RF.B1"),
                       accept="the abstract's band-offset sentence matches "
                              "a test with a stated p-value and trials "
                              "factor; 'k of N significant' appears "
                              "nowhere (standing rule 2)."),
                  Task("CV-R2-bias-injection",
                       "Per-band, per-era injection with the signed bias "
                       "carried",
                       "Injection grids per band with per-band ramp "
                       "widths, and on a 2024 High Gain night; the SIGNED "
                       "matched-cell bias carried per band and era; no "
                       "max(χ²ν, 1); no transported budget. This is the "
                       "test that decides D3.",
                       "CV-S9", _cv_src(f"{_AMENDMENTS} CV-R2-bias-injection"),
                       PENDING,
                       ruling=_ruled("§4 CV-R2-bias-injection",
                                     "D3", "RF.B2", "RF.B3", "DS.F3"),
                       accept="differential g−i estimator bias ± error "
                              "tabulated; no epoch carries a transported "
                              "budget; the paper states whichever reading "
                              "of D3 the test supports and no more."),
                  Task("CV-R3-oc-refit",
                       "Refit the O−C: per-band constants, night-level "
                       "epochs",
                       "The O−C refitted with per-band constants, "
                       "night-level epochs (one per night, not one per "
                       "band per night) and an era-offset nuisance term, "
                       "with the period-change bound quoted each way and "
                       "with 2024 dropped.",
                       "CV-S9", _cv_src(f"{_AMENDMENTS} CV-R3-oc-refit"),
                       PENDING, depends_on=("CV-R2-bias-injection",),
                       ruling=_ruled("§4 CV-R3-oc-refit", "DS.F2", "RF.M3"),
                       accept="the bound is quoted for every variant; "
                              "χ²ν reported per band with its dof "
                              "(standing rule 1); if the bound moves by "
                              "more than 20% the abstract quotes the "
                              "range."),
                  Task("CV-R4-superhump-wording",
                       "State the superhump result as what was excluded",
                       "The YZ Cnc result restated as an excluded "
                       "semi-amplitude range, in how many of the "
                       "run-filters, with the published floor cited; "
                       "'measurement of absence' and 'is real — and is "
                       "not orbital' struck. A recovery contour above the "
                       "smallest expected signal means that signal would "
                       "NOT have been found.",
                       "CV-S10",
                       _cv_src(f"{_AMENDMENTS} CV-R4-superhump-wording"),
                       PENDING,
                       ruling=_ruled("§4 CV-R4-superhump-wording",
                                     "PH.P1", "RF.B4", "RF.M7"),
                       accept="no sentence claims superhumps were "
                              "excluded below the recovery contour."),
                  Task("CV-R5-pdot-physics",
                       "Say what the period-change bound constrains",
                       "The bound relabelled as a spin/spot-longitude "
                       "bound, set beside three physical scales (secular "
                       "orbital evolution, an asynchronous polar relaxing, "
                       "spot-longitude drift), with the longitude "
                       "stability quoted in degrees and split by accretion "
                       "state.",
                       "CV-S9", _cv_src(f"{_AMENDMENTS} CV-R5-pdot-physics"),
                       PENDING, depends_on=("CV-R3-oc-refit",),
                       ruling=_ruled("§4 CV-R5-pdot-physics",
                                     "PH.P2", "RF.M3"),
                       accept="the bound is labelled spin/spot with its "
                              "three scales beside it; longitude "
                              "stability is in degrees, per state."),
                  Task("CV-R6-colour-result",
                       "Quantify the colour–phase result",
                       "ST LMi's g−r swing quantified — amplitude, phase "
                       "of the extrema, two-era repeatability with errors "
                       "— shown to survive a ≤ 120 s pairing window, with "
                       "its cyclotron context. The abstract promised "
                       "'colour results' and the paper stated none.",
                       "CV-S9", _cv_src(f"{_AMENDMENTS} CV-R6-colour-result"),
                       PENDING,
                       ruling=_ruled("§4 CV-R6-colour-result",
                                     "ED.E2", "RF.M1", "PH.P6"),
                       accept="amplitude and phase quoted with errors for "
                              "both eras, unchanged within errors at a "
                              "≤ 120 s pairing window."),
                  Task("CV-R7-edge-fit-disclosure",
                       "Disclose the edge fit",
                       "The edge model written down, the search grid, the "
                       "per-edge χ²ν distribution, and a best, a "
                       "representative and a worst fit shown.",
                       "CV-S9",
                       _cv_src(f"{_AMENDMENTS} CV-R7-edge-fit-disclosure"),
                       PENDING,
                       ruling=_ruled("§4 CV-R7-edge-fit-disclosure",
                                     "RF.M2"),
                       accept="model equation, grid and χ²ν distribution "
                              "are in the paper, per band and era with "
                              "dof."),
                  Task("CV-R8-counts",
                       "Make the counts mean what they say",
                       "Staged and measured frame counts stated "
                       "separately; a signal-to-noise floor on what is "
                       "called a measurement; EU UMa i and r out of the "
                       "count; the frames with FWHM beyond the aperture "
                       "cut or shown harmless.",
                       "CV-S5", _cv_src(f"{_AMENDMENTS} CV-R8-counts"),
                       PENDING,
                       ruling=_ruled("§4 CV-R8-counts",
                                     "RF.M5", "RF.M6", "OA"),
                       accept="every counted series has a quoted "
                              "precision; staged and measured counts are "
                              "separate macros."),
                  Task("CV-R9-reduction",
                       "Describe and test the reduction",
                       "A reduction paragraph (flats, darks, the reduced "
                       "tree, the aperture); one night per era "
                       "photometered from raw+masters and from the "
                       "server-reduced frames; the flat-field ramp tested "
                       "as catalogue-tie residual against detector x; a "
                       "mechanical epoch per series.",
                       "CV-S4", _cv_src(f"{_AMENDMENTS} CV-R9-reduction"),
                       PENDING, depends_on=("F-3",),
                       ruling=_ruled("§4 CV-R9-reduction",
                                     "OA.E5", "TE.F5"),
                       accept="raw and reduced differential light curves "
                              "agree within the check-star scatter; tie "
                              "residuals are flat in x to < 0.5%, or the "
                              "ramp is quantified in the photometry "
                              "section; no fold combines mechanical "
                              "states without an offset."),
                  Task("CV-R10-instrument-section",
                       "Rewrite the instrument section on measured numbers",
                       "Measured gains and read noise in place of the "
                       "bracket; the StackPro and clock-card sentences "
                       "corrected; a rolling-shutter bound; a "
                       "residual-vs-peak linearity panel.",
                       "S2",
                       _cv_src(f"{_AMENDMENTS} CV-R10-instrument-section"),
                       PENDING, depends_on=("F-4", "F-5"),
                       ruling=_ruled("§4 CV-R10-instrument-section",
                                     "DE.F2", "DE.F3", "DE.F7"),
                       accept="the gain paragraph cites detector_params; "
                              "error-inflation factors are re-quoted per "
                              "mode with the measured gains; linearity "
                              "slope < 1% to the adopted cap is shown."),
                  Task("CV-R11-clock",
                       "Print the clock residual; test the clock on "
                       "transits",
                       "The one eclipse residual the paper called 'weak' "
                       "and did not print, printed and explained; the "
                       "absolute clock tested per camera era on archived "
                       "transits; what that implies for the common "
                       "~1,000 s edge offset of ST LMi and EU UMa.",
                       "S3", _cv_src(f"{_AMENDMENTS} CV-R11-clock"),
                       PENDING, depends_on=("F-8",),
                       ruling=_ruled("§4 CV-R11-clock", "RF.M4", "OA.E6"),
                       accept="|O−C| < 120 s per era, or the offset is "
                              "measured and carried into every absolute "
                              "epoch."),
                  Task("CV-R12-literature",
                       "Meet the literature",
                       "Every external dataset, catalogue and ephemeris "
                       "cited, and a comparison-with-previous-work "
                       "subsection for ST LMi and YZ Cnc.",
                       "R-CV-S11", _cv_src(f"{_AMENDMENTS} CV-R12-literature"),
                       PENDING,
                       ruling=_ruled("§4 CV-R12-literature",
                                     "ED.E1", "RF.M8"),
                       accept="about forty verified references or more; "
                              "no external resource used without a "
                              "citation; the build log still clean."),
                  Task("CV-R13-restructure",
                       "Restructure: the colour result leads",
                       "An abstract of ≤ 250 words that leads with the "
                       "colour result; ≤ 18 pages; process prose moved to "
                       "an appendix; internal language removed; venue AJ.",
                       "R-CV-S11",
                       _cv_src(f"{_AMENDMENTS} CV-R13-restructure"),
                       PENDING,
                       depends_on=("CV-R1-band-offset",
                                   "CV-R4-superhump-wording",
                                   "CV-R5-pdot-physics",
                                   "CV-R6-colour-result"),
                       ruling=_ruled("§4 CV-R13-restructure",
                                     "ED.E2", "ED.E3", "RF"),
                       accept="the abstract is approved by seat 6 "
                              "(standing rule 5); ≤ 18 pages; every "
                              "main-text figure supports a sentence of "
                              "the abstract."),
                  Task("CV-R14-rebuild",
                       "Rebuild the chain FRESH and file the macro diff",
                       "The CV chain rebuilt FRESH at a clean, tagged "
                       "commit, with a macro-by-macro diff of numbers.tex "
                       "filed. 'No macro moved' is an acceptable outcome; "
                       "it has to be shown.",
                       "R-CV-S11", _cv_src(f"{_AMENDMENTS} CV-R14-rebuild"),
                       PENDING, depends_on=("F-10",),
                       ruling=_ruled("§4 CV-R14-rebuild",
                                     "U5", "ED.E5", "RF.M9"),
                       accept="zero STALE CV stages at a non-dirty "
                              "commit; every moved macro explained."),
                  Task("CV-R15-release-readiness",
                       "Release readiness: authors, ORCIDs, DOI, outside "
                       "reader",
                       "An authorship policy agreed by the consortium, "
                       "real ORCIDs, a reserved Zenodo DOI, AAS data- and "
                       "software-availability statements, and one outside "
                       "reader who works on polars.",
                       "R-CV-S11",
                       _cv_src(f"{_AMENDMENTS} CV-R15-release-readiness"),
                       BLOCKED,
                       blocker="Needs James (SYNTHESIS §6): the authorship "
                               "policy is the consortium's to agree, the "
                               "ORCIDs and the DOI are his to obtain, and "
                               "the outside reader is his to ask. A draft "
                               "policy and draft e-mails are prepared; "
                               "nothing clears until he takes them out.",
                       ruling=_ruled("§4 CV-R15-release-readiness",
                                     "ED.E4"),
                       accept="AUTHORSHIP.md agreed; no placeholder "
                              "ORCID; DOI reserved; the outside reader's "
                              "one-sentence summary of the result matches "
                              "ours."),
              )),
    ),
)


_TCRB = "TCrB_Monitoring/ANALYSIS_STRATEGY.md"

def _tcrb_src(section: str) -> Source:
    return Source(_TCRB, section)


TCRB_MONITORING = Project(
    key="TCrB_Monitoring",
    title="T CrB Pre-Eruption Monitoring",
    claim=(
        "An Hα equivalent-width series of T CrB through the pre-eruption "
        "dip and recovery, Feb–Jun 2025, from the slitless grism at ~2-day "
        "cadence, cross-validated against ARAS. Equivalent width is the "
        "observable; line flux is EW × a contemporaneous continuum light "
        "curve, never a θ CrB response (θ CrB is a variable Be star and the "
        "absolute-flux step was dropped), and both are shown against "
        "orbital phase because the giant's ellipsoidal continuum modulates "
        "any raw EW. How many frames are spectra is reported by S2c "
        "verdict, not by filter name. Photometric anchors are B-only — R "
        "and I are clipped at the target. Profile morphology and wing "
        "velocities are HELD until the hrg dispersion disagreement (D1) is "
        "settled by test. Flickering and the period search have left this "
        "paper: no ≥2 h run is possible before mid-January 2027. Whether "
        "the series is new is decided FIRST, by the novelty table "
        "(TCRB-N1). The paper makes zero eruption-date predictions."),
    venue="ApJ only if TCRB-N1 shows a density or homogeneity advantage "
          "over ARAS; otherwise AJ or a PASP validation note",
    paper_title=(r"H$\alpha$ Equivalent-Width Monitoring of the Recurrent "
                 r"Nova T~CrB through the 2025 Pre-Eruption Dip and Recovery "
                 r"with a Slitless Grism on the Robert L. Mutel Telescope"),
    strategy=_TCRB,
    ruling=_ruled("§4 TCrB_Monitoring", "U2", "U6", "D1"),
    decisions=(
        ("The strategy states no timing rule, and that is a gap, not a "
         "silence to fill",
         "Mid-exposure BJD_TDB is real, built, and published (S3) — but "
         "this project's ANALYSIS_STRATEGY.md contains no BJD, TDB, "
         "barycentric or mid-exposure rule anywhere in it. The "
         "timing task therefore cites §9's facility-level products, which "
         "is where shared machinery is authorised, and NOT a Phase A/B "
         "timing rule, which does not exist. It was cited that way once, "
         "on this page, hyperlinked to a document that did not contain it; "
         "the citation checker now refuses that. When the committee next "
         "revises the strategy, a per-project timing rule belongs in it."),
        ("Plan review of 2026-10-03 — what now binds this paper",
         "Execute with amendments; highest priority in the portfolio. "
         "(1) The novelty table comes first and can re-scope the paper to a "
         "validation note. (2) Grism dispersion is hardware: one solution "
         "per (grism, mechanical epoch) from hot stars, per frame only a "
         "zero point (U2). (3) Whether hrg is ≈0.47 or ≈1.59 Å/px is "
         "decided by test (D1); profile and velocity work is held until it "
         "returns, neither dropped nor kept. (4) Absolute flux from a θ CrB "
         "response is dropped; line flux is EW × continuum. (5) Anchors "
         "are B-only. (6) The period search is dropped and the archival "
         "flickering limits shrink to one sentence and a table; the "
         "2026–27 flickering runs and the splice onto the QHY (a second "
         "instrument) are the 2027 backlog (U6). (7) Six figures, not ten."),
    ),
    phases=(
        Phase("Novelty gate — first, and alone",
              "Days of reading that decide the venue and whether this is a "
              "science paper or a validation note. Nothing that builds the "
              "EW series is recommended until it reports.",
              (
                  Task("TCRB-N1-novelty-table",
                       "Is the RLMT Hα series new? Spectra per month vs "
                       "ARAS and Asiago",
                       "Spectra per month, RLMT vs ARAS vs Asiago/Munari, "
                       "Feb–Jun 2025, with resolution, and the published Hα "
                       "EW series covering the same window.",
                       "S0c", _tcrb_src(f"{_AMENDMENTS} TCRB-N1-novelty-table"),
                       PENDING,
                       ruling=_ruled("§4 TCRB-N1-novelty-table", "ED"),
                       accept="table on disk with sources; if ARAS alone is "
                              "denser and higher-resolution, the paper "
                              "re-scopes to a validation/methods note and "
                              "says so."),
              )),
        Phase("Phase 0 — Gates",
              "What must be true of the instrument before a spectrum is "
              "measured. The S2-gated tasks are unblocked: the tables a "
              "stale blocker called destroyed are in the database.",
              (
                  Task("TCRB-P0-staging", "Working set staged by reference",
                       "stage_tcrb_monitoring: T CrB and θ CrB science "
                       "frames with their era-matched calibration.",
                       "S0c", _tcrb_src("§9 catalog hygiene"), DONE,
                       evidence="docs/pipeline/s0c_staging.html"),
                  Task("TCRB-P0-timing", "Mid-exposure BJD_TDB",
                       "frame_times entries for the imaging and grism "
                       "series — a FACILITY-level product, built once for "
                       "every project; T CrB's own strategy states no "
                       "timing rule (see the standing decision).",
                       "S3", _tcrb_src("§9 facility-level products"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("TCRB-P0-mech-epoch",
                       "Assign every staged frame a mechanical epoch",
                       "A mech_epoch per staged frame (camera, rotation, "
                       "flip, wheel map), and the 2025 grism series "
                       "confirmed to sit inside one.",
                       "S0", _tcrb_src(f"{_AMENDMENTS} TCRB-P0-mech-epoch"),
                       PENDING, depends_on=("F-3",),
                       ruling=_ruled("§4 TCRB-P0-mech-epoch", "TE.F1"),
                       accept="table emitted; the 2025 series confirmed "
                              "single-epoch; no master, trace prior or "
                              "zero-order prediction crosses a boundary."),
                  Task("TCRB-P0-gain-ptc",
                       "Measured gain and read noise for the grism camera",
                       "Flat-pair gain and read noise for Mode0 (and for "
                       "the QHY once October flats exist) in "
                       "detector_params, read by the grism variance — the "
                       "header EGAIN is the unbinned gain and is wrong for "
                       "average-binned frames.",
                       "S2", _tcrb_src(f"{_AMENDMENTS} TCRB-P0-gain-ptc"),
                       PENDING, depends_on=("F-4",),
                       ruling=_ruled("§4 TCrB_Monitoring gain-ptc", "D2", "DE.F1"),
                       accept="K ± 3% and read noise in detector_params; "
                              "the grism variance and saturation threshold "
                              "read them from there."),
                  Task("TCRB-P0-temp-split",
                       "Split the grism frames by detector temperature",
                       "The grism series split by CCD-TEMP, with the dark "
                       "treatment stated per temperature group — a third of "
                       "the 240 s frames were taken warm against −10 °C "
                       "master darks.",
                       "S0", _tcrb_src(f"{_AMENDMENTS} TCRB-P0-temp-split"),
                       PENDING, depends_on=("F-2",),
                       ruling=_ruled("§4 TCrB_Monitoring temp-split", "DE.F5"),
                       accept="dark treatment stated per temperature group; "
                              "flanking-band adequacy shown separately for "
                              "the warm frames."),
                  Task("TCRB-P0-eruption-block",
                       "The eruption-response block, as an artefact",
                       "The pyscope eruption-response schedule committed "
                       "under ops/ with a dry-run log. The strategy says it "
                       "was to be loaded during the closure; nothing on "
                       "disk shows that it was.",
                       "OPS", _tcrb_src(f"{_AMENDMENTS} TCRB-P0-eruption-block"),
                       PENDING,
                       ruling=_ruled("§4 TCRB-P0-eruption-block", "TE"),
                       accept="the schedule file is committed under ops/ "
                              "with a dry-run log beside it."),
                  Task("TCRB-P0-filter-forensics",
                       "Map the single-character filter codes",
                       "A filter-mapping table for codes 6,1,G,H,L,O,W — 62 "
                       "frames are unusable for calibrated work until it "
                       "exists.",
                       "S0", _tcrb_src("§4 Phase 0 P0-2"), PENDING),
                  Task("TCRB-P0-bitdepth",
                       "High Gain bit depth — closed from the archive",
                       "The High Gain clip, 12-bit-consistent and measured "
                       "per EGAIN epoch from the archive itself "
                       "(s2_ceiling_modes). The one-afternoon hardware test "
                       "is retired with the camera, which left the beam in "
                       "March 2024; it is re-opened only if the AC4040 "
                       "turns out to be on site.",
                       "S2", _tcrb_src("§4 Phase 0 P0-3(i)"), DONE,
                       evidence="docs/pipeline/s2_detector.html",
                       ruling=_ruled("§4 TCrB_Monitoring P0-bitdepth",
                                     "U4", "DE.F6", "OA.E3"),
                       probes=(_s2_table_probe("s2_ceiling_modes"),)),
                  Task("TCRB-P0-ladders",
                       "One linearity cap per readout mode, from "
                       "residual-vs-peak",
                       "One adopted photometry cap per mode, set where the "
                       "comparison-star residual-vs-own-peak slope leaves "
                       "the error floor — replacing the planned 70% cap, "
                       "which nobody measured, and the three other caps "
                       "(70%, 80%, 92%) in force across the projects.",
                       "S2", _tcrb_src("§4 Phase 0 P0-3(ii)"), PENDING,
                       depends_on=("F-5",),
                       ruling=_ruled("§4 TCrB_Monitoring UNBLOCK ladders",
                                     "U4", "DE.F3", "DE.F6"),
                       accept="residual-vs-peak slope < 1% up to the "
                              "adopted cap in each mode, judged at the "
                              "target in native pixels.",
                       probes=(_s2_table_probe("s2_linearity_ladders"),)),
                  Task("TCRB-P0-calib-acquisition",
                       "State the dark method; put the bench-dark question "
                       "to the site",
                       "The method for the 240 s Mode0 frames, stated as "
                       "the method: the existing −10 °C masters plus "
                       "flanking-band subtraction, with a penalty term for "
                       "the warm frames. Mode0 darks cannot be acquired on "
                       "the QHY now mounted, so the acquisition this task "
                       "used to wait for is a question for the site — does "
                       "the ASI still exist, can it be cooled on a bench?",
                       "G", _tcrb_src("§4 Phase 0 P0-3(iii)"), PENDING,
                       depends_on=("TCRB-P0-temp-split",),
                       ruling=_ruled("§4 TCrB_Monitoring "
                                     "P0-calib-acquisition", "TE", "DE"),
                       accept="method paragraph and its adequacy figure "
                              "exist; the bench-dark question is in the "
                              "rev. 3 observatory request."),
                  Task("TCRB-P0-shutter-timing",
                       "Rolling-shutter skew and exposure accuracy at short "
                       "exposures",
                       "Row-sequential read skew and driver exposure "
                       "accuracy, bounded on the archival 0.085 s θ CrB "
                       "frames. Those are electronic-shutter CMOS frames: "
                       "what needs testing is skew, not shutter travel.",
                       "S2", _tcrb_src("§4 Phase 0 P0-3(iv)"), PENDING,
                       ruling=_ruled("§4 TCrB_Monitoring shutter-timing",
                                     "U4", "DE"),
                       accept="skew and exposure-time error stated as "
                              "bounds in seconds before the eruption plan "
                              "relies on 0.1 s."),
                  Task("TCRB-P0-zmag-provenance", "ZMAG provenance per filter",
                       "Which catalog bandpass PinPoint solved each filter's "
                       "ZMAG against — without it the QC cut of step B4 is "
                       "applied to a quantity nobody has defined.",
                       "S0", _tcrb_src("§4 Phase 0 P0-3(v)"), PENDING),
                  Task("TCRB-P0-resolve", "Re-solve the unsolved imaging",
                       "WCS for the ~40% of imaging frames that lack it; "
                       "grism frames get the identity gate instead.",
                       "S1b", _tcrb_src("§4 Phase 0 P0-4"), IN_PROGRESS),
              )),
        Phase("Phase A — Grism (the paper's spine)",
              "From per-frame identity to the cross-validated EW series, on "
              "the shared grism library (G-1…G-5); no frame that fails the "
              "gate enters anything.",
              (
                  Task("TCRB-A0-identity-gate",
                       "Per-frame identity gate, from the pixels",
                       "A pass/fail per frame from the spectrum itself — "
                       "the Hα emission plus the TiO 7054 step at a fixed "
                       "pixel offset — with the discard fraction published "
                       "and the counts reported by S2c verdict. The header "
                       "is never the sole reason: frames it called >1° off "
                       "are on target with a stale RA/Dec card.",
                       "G", _tcrb_src("§4 Phase A step 0"), IN_PROGRESS,
                       depends_on=("G-3", "F-6"),
                       ruling=_ruled("§4 TCrB_Monitoring A0",
                                     "OA.E2", "RF", "DS"),
                       accept="zero rejections whose sole reason is a "
                              "header; false-accept rate measured on "
                              "non-T CrB grism frames; N reported by "
                              "verdict, on unique frames."),
                  Task("TCRB-A0b-early-spectra",
                       "Extract or reject the slot '6' and 'W' spectra, "
                       "2023-05 → 2024-03",
                       "An extract-or-reject decision per frame for the "
                       "measured-dispersed T CrB frames under filter codes "
                       "'6' and 'W' — possibly the only RLMT Hα points "
                       "before and inside the dip.",
                       "G", _tcrb_src(f"{_AMENDMENTS} TCRB-A0b"), PENDING,
                       depends_on=("G-1", "F-6"),
                       ruling=_ruled("§4 TCRB-A0b", "DS", "OA"),
                       accept="EW ± error per frame, or a documented "
                              "saturation/contamination failure."),
                  Task("TCRB-A1-calibrator-characterization",
                       "Characterise θ CrB before using θ CrB",
                       "Its own Hα/Hβ profile and stability across all 412 "
                       "frames — it is a Be/shell star, not a clean standard.",
                       "G", _tcrb_src("§4 Phase A step 1"), PENDING),
                  Task("TCRB-A2-extraction",
                       "Trace and optimal extraction",
                       "Boxcar and Horne-extracted 1-D spectra per frame, "
                       "with the variance built from the MEASURED gain and "
                       "the background from a sky-lozenge template before "
                       "the flanking bands. The extraction in progress used "
                       "the header gain (≈4× wrong for average-binned "
                       "frames) and masked valid pixels above a 16.3 kADU "
                       "'rail' that does not exist.",
                       "G", _tcrb_src("§4 Phase A step 2"), IN_PROGRESS,
                       depends_on=("G-2", "G-4"),
                       ruling=_ruled("§3 G-2", "D2", "DE.F1", "OA.E7"),
                       accept="predicted variance within 20% of the "
                              "flanking-band variance; no negative "
                              "continuum; boxcar-vs-optimal difference "
                              "< 3%."),
                  Task("TCRB-A3-wavelength",
                       "Fixed dispersion per (grism, mechanical epoch); "
                       "per-frame zero point only",
                       "One dispersion and sign per grism per mechanical "
                       "epoch, taken from the shared library, and a zero "
                       "point per frame (Hα, with O₂-B as the consistency "
                       "check), with the km/s uncertainty quoted beside "
                       "every velocity. The per-frame dispersions on disk "
                       "(−1.8 to +2.0 Å/px) are not physical and are "
                       "retired.",
                       "G", _tcrb_src("§4 Phase A step 3"), PENDING,
                       depends_on=("G-1",),
                       ruling=_ruled("§4 TCrB_Monitoring A3",
                                     "U2", "D1", "OA.E1", "PH.P7", "TE.F3"),
                       accept="the Hα–O₂B pixel separation is constant to "
                              "< 1–2% across every T CrB frame."),
                  Task("TCRB-A5a-detection-rule",
                       "Pre-register the EW-change detection rule",
                       "The rule for calling an EW change real — threshold "
                       "and consecutive-night requirement — written into "
                       "the strategy before the series is looked at.",
                       "STRAT", _tcrb_src(f"{_AMENDMENTS} TCRB-A5a"),
                       PENDING,
                       ruling=_ruled("§4 TCrB_Monitoring A5", "DS"),
                       accept="the rule is in the strategy, dated, before "
                              "TCRB-A5-ew runs."),
                  Task("TCRB-A5-ew",
                       "Hα equivalent width on Munari's convention",
                       "THE primary observable: a season-long EW series that "
                       "splices onto the published curves, on windows shown "
                       "to survive this instrument's resolution.",
                       "G", _tcrb_src("§4 Phase A step 5"), PENDING,
                       depends_on=("TCRB-N1-novelty-table",
                                   "TCRB-A3-wavelength",
                                   "TCRB-A5a-detection-rule", "G-5"),
                       forbids="SYNTHESIS §4: the novelty table is first, "
                               "and the detection rule is pre-registered "
                               "\"before A5 runs\". An EW series measured "
                               "before either is a result looking for its "
                               "own threshold.",
                       ruling=_ruled("§4 TCrB_Monitoring A5", "PH", "DS"),
                       accept="ARAS spectra degraded to the measured LSF "
                              "reproduce their native EW within the stated "
                              "tolerance; windows documented and tested "
                              "against ±10 Å shifts."),
                  Task("TCRB-A5b-line-flux",
                       "Line flux as EW × continuum, against orbital phase",
                       "F(Hα) = EW × the continuum flux from "
                       "contemporaneous AAVSO or zero-order photometry, "
                       "with the photometric source named per epoch, and "
                       "EW and flux both plotted against orbital phase — "
                       "the giant's ellipsoidal continuum puts a guaranteed "
                       "geometric oscillation into any raw EW series. "
                       "Replaces the dropped θ CrB response.",
                       "G", _tcrb_src(f"{_AMENDMENTS} TCRB-A5b"), PENDING,
                       depends_on=("TCRB-A5-ew", "TCRB-D1-external-pulls"),
                       ruling=_ruled("§4 TCRB-A5b", "OA", "PH"),
                       accept="a continuum light curve accompanies the EW "
                              "series (standing rule 6); photometric source "
                              "named per epoch; both series shown vs "
                              "orbital phase."),
                  Task("TCRB-A6-saturation-triage",
                       "Scripted saturation triage of every grism frame",
                       "Flagged and discarded frames with the fraction "
                       "reported, judged against the measured clip with a "
                       "hot-pixel mask — not against 16.3 kADU, which is "
                       "the signature of one saturated native pixel "
                       "averaged with three normal ones.",
                       "G", _tcrb_src("§4 Phase A step 6"), PENDING,
                       depends_on=("G-2",),
                       ruling=_ruled("§4 TCrB_Monitoring A6",
                                     "D2", "DE.F1", "OA.E8"),
                       accept="n_sat_cols recomputed from detector_params; "
                              "masked fraction reported; saturation judged "
                              "at the target in native pixels."),
                  Task("TCRB-A7-error-floor",
                       "Rebuild the empirical error floor",
                       "θ CrB night-to-night EW scatter from the matched "
                       "2025 Mode0 subset, split-half and cross-checked on "
                       "a line-free window, with a focus-offset regressor "
                       "and the extraction-method difference carried as an "
                       "EW floor, plus a separate 240 s smear term.",
                       "G", _tcrb_src("§4 Phase A step 7"), PENDING,
                       depends_on=("G-5",),
                       ruling=_ruled("§4 TCrB_Monitoring A7",
                                     "DS", "TE.F4", "RF"),
                       accept="χ²ν of the θ CrB continuum-window series in "
                              "[0.7, 1.4] with its dof; EW floor ≥ the "
                              "measured extraction-method difference; "
                              "off-nominal-focus nights flagged."),
                  Task("TCRB-A8-cross-validation",
                       "ARAS/published EW cross-validation",
                       "Agreement within 10–15% or an explanation. The "
                       "spectra-per-month table this step used to deliver "
                       "is now TCRB-N1 and comes first.",
                       "G", _tcrb_src("§4 Phase A step 8"), PENDING),
                  Task("TCRB-A9-profiles",
                       "Profile morphology and wing velocities — held "
                       "pending D1",
                       "Line-profile and differential wing-velocity "
                       "evolution from hrg, if and only if the delivered "
                       "resolution supports it.",
                       "G", _tcrb_src(f"{_AMENDMENTS} TCRB-A9-profiles"),
                       BLOCKED, depends_on=("G-1",),
                       blocker="HELD by ruling D1 — neither dropped nor "
                               "kept. The seats disagree on the hrg "
                               "dispersion by a factor 3.4 (≈0.47 vs "
                               "≈1.59 Å/px), which is the difference "
                               "between resolving T CrB's Hα and not. "
                               "Clears when G-1 returns its test on Vega "
                               "and θ CrB (≥3 lines, residual < 1 px): the "
                               "task then opens on the measured LSF, or is "
                               "dropped by it.",
                       ruling=_ruled("§2 D1", "D1", "OA.E1", "PH.P8", "RF",
                                     action="HOLD"),
                       accept="delivered FWHM in km/s stated before any "
                              "velocity is plotted; nights at off-nominal "
                              "focus excluded."),
              )),
        Phase("Phase B — Photometry (supporting, B-only)",
              "The 2023–2024 anchors that survive a peak-at-target census: "
              "B, plus singletons. R and I are flat-topped at the clip.",
              (
                  Task("TCRB-B0-peak-census",
                       "Peak-at-target census of every imaging frame",
                       "Peak ADU at the target per frame, in "
                       "native-pixel-equivalent units, with a clipped/clean "
                       "verdict — the census that decides which anchors "
                       "exist.",
                       "S2", _tcrb_src(f"{_AMENDMENTS} TCRB-B0-peak-census"),
                       PENDING,
                       ruling=_ruled("§4 TCrB_Monitoring Phase B", "OA.E3"),
                       accept="every imaging frame carries a "
                              "peak-at-target and a verdict; the anchor "
                              "set is restated from it."),
                  Task("TCRB-B1-calibration",
                       "Calibrate per readout mode with measured penalties",
                       "Mode-matched masters where they exist; a MEASURED "
                       "cross-mode penalty term where they never will.",
                       "S2", _tcrb_src("§4 Phase B step 1 (ruling 7)"),
                       PENDING, depends_on=("F-4",),
                       ruling=_ruled("§4 TCrB_Monitoring UNBLOCK B1",
                                     "U4", "TE.F11"),
                       accept="the cross-mode penalty is measured as "
                              "ensemble residual vs detector position and "
                              "carried in the systematic budget."),
                  Task("TCRB-B2-aperture", "Aperture photometry with growth curves",
                       "Per-frame photometry at r = 1.5×FWHM, cross-checked "
                       "against fixed and wide apertures.",
                       "S4", _tcrb_src("§4 Phase B step 2"), PENDING,
                       depends_on=("TCRB-B0-peak-census",)),
                  Task("TCRB-B3-ensemble",
                       "REFCAT2 ensemble with propagated colour extrapolation",
                       "Natural-system and transformed magnitudes, with the "
                       "M4III colour-extrapolation uncertainty separated.",
                       "S4", _tcrb_src("§4 Phase B step 3"), PENDING),
                  Task("TCRB-B4-zmag-qc", "ZMAG demoted to QC only",
                       "A per-frame QC flag cutting frames >0.5 mag below "
                       "the per-filter season mode, with M1's 0.3 mag flag "
                       "recorded alongside so the paper can tabulate how "
                       "sensitive the results are to the threshold. ZMAG is "
                       "never used as calibration.",
                       "S4", _tcrb_src("§4 Phase B step 4 (ruling 4)"),
                       PENDING, depends_on=("TCRB-P0-zmag-provenance",),
                       forbids="Ruling 4 settled a four-way disagreement "
                               "about the cut threshold (0.3 vs 0.5 vs 1.0 "
                               "mag); applying it before P0-3(v) establishes "
                               "what ZMAG is measured against would cut "
                               "frames on an undefined quantity."),
                  Task("TCRB-B5-errors", "Empirical per-frame errors",
                       "Check-star RMS per night per mode — no Poisson "
                       "arithmetic on StackPro, ever.",
                       "S4", _tcrb_src("§4 Phase B step 5"), PENDING),
                  Task("TCRB-B6-precision-budget",
                       "State the precision achieved",
                       "The achieved per-frame precision and season-level "
                       "systematic floor, stated from the check-star rms. "
                       "The 3–5 mmag nightly-mean promise is deleted: three "
                       "B nights cannot carry it.",
                       "S4", _tcrb_src("§4 Phase B step 6"), PENDING,
                       depends_on=("TCRB-B5-errors",),
                       ruling=_ruled("§4 TCrB_Monitoring Phase B", "OA"),
                       accept="precision quoted per mode from check stars; "
                              "no nightly-mean promise anywhere in the "
                              "paper."),
              )),
        Phase("Phase C — Time series",
              "What the archive can bound, in one sentence and one table. "
              "The period search is dropped and the flickering runs are the "
              "2027 backlog.",
              (
                  Task("TCRB-C1-flickering-limits",
                       "Archival flickering limits — a sentence and a table",
                       "One sentence and a one-table statement of the "
                       "detrended rms and 95% σ_flick upper limit from the "
                       "surviving B snippets, with the window function. No "
                       "figure, no state comparison, no periodogram.",
                       "S4", _tcrb_src("§4 Phase C step 1"), PENDING,
                       ruling=_ruled("§4 TCrB_Monitoring C1",
                                     "U6", "DS", "OA"),
                       accept="each limit carries the injected amplitude "
                              "it would have recovered and the predicted "
                              "scale beside it (standing rule 2)."),
                  Task("TCRB-C4-uncertainties",
                       "Uncertainties on every quoted number",
                       "emcee posteriors, a red-noise inflation factor from "
                       "binned comparison-star rms, and every headline value "
                       "quoted as value ± stat ± sys with the linearity and "
                       "cross-mode penalty terms listed separately.",
                       "S4", _tcrb_src("§4 Phase C step 4"), PENDING),
              )),
        Phase("Phase D — Co-analysis and release",
              "External data, the centrepiece figure, the archive deposit.",
              (
                  Task("TCRB-D1-external-pulls",
                       "Pull AAVSO / ASAS-SN / ARAS / Swift / TESS",
                       "A cached external context set with pull dates, "
                       "fetched now rather than in October.",
                       "S0c", _tcrb_src("§4 Phase D"), PENDING),
                  Task("TCRB-D2-centerpiece",
                       "The centrepiece figure",
                       "The RLMT Hα EW series over the AAVSO B curve through "
                       "high state → dip → recovery.",
                       "G", _tcrb_src("§7 Figure 8"), PENDING),
                  Task("TCRB-N2-eruption-contingency",
                       "Eruption contingency skeleton",
                       "A pre-written RNAAS/ATel skeleton with the 2025 EW "
                       "series as its baseline, every number "
                       "script-emitted, ready for the night the nova goes.",
                       "G", _tcrb_src(f"{_AMENDMENTS} "
                                      f"TCRB-N2-eruption-contingency"),
                       PENDING, depends_on=("TCRB-A5-ew",),
                       ruling=_ruled("§4 TCRB-N2-eruption-contingency", "ED"),
                       accept="the skeleton compiles with script-emitted "
                              "numbers."),
                  Task("TCRB-D3-release",
                       "Release the reduced 1-D spectra",
                       "A Zenodo deposit plus a submission to the ARAS T CrB "
                       "database.",
                       "G", _tcrb_src("§4 Phase D"), PENDING),
                  Task("TCRB-D4-figures", "The six-figure set",
                       "Six figures, each mapped to a sentence of the "
                       "abstract: the centrepiece (EW over the AAVSO B "
                       "curve), the spectra montage, the calibrator floor, "
                       "the ARAS cross-validation, EW and line flux against "
                       "orbital phase, and the measured LSF. Flickering "
                       "limits are a table. Cut from ten.",
                       "G", _tcrb_src("§7 Figure list"), PENDING,
                       ruling=_ruled("§4 TCrB_Monitoring D4", "ED"),
                       accept="six figures or fewer, each cited in support "
                              "of an abstract sentence; the ≤ 250-word "
                              "abstract is approved by seat 6 before they "
                              "are built (standing rule 5)."),
                  Task("TCRB-D5-draft", "Write the manuscript",
                       "manuscripts/TCrB_Monitoring/main.tex through the "
                       "§8 outline. Its first action — changing the title, "
                       "because the skeleton's 'Three-Year Quiescent "
                       "Baseline' claimed exactly what the panel rejected — "
                       "was carried out on 2026-10-03.",
                       "G", _tcrb_src("§8 Manuscript outline"), PENDING,
                       depends_on=("TCRB-D4-figures",)),
              )),
        _dropped_phase(
            Task("TCRB-A4-response",
                 "Instrument response from θ CrB → absolute flux",
                 "Nothing: dropped. θ CrB is variable, Hα-contaminated, "
                 "exposed 0.6–5 s against the target's 240 s, and its "
                 "second-order blue light lands on the first-order red "
                 "where T CrB has none. Line flux now comes from "
                 "TCRB-A5b (EW × continuum).",
                 "G", _tcrb_src("§4 Phase A step 4"), DROPPED,
                 ruling=_ruled("§4 TCrB_Monitoring DROP A4", "OA", "PH")),
            Task("TCRB-C3-period-search",
                 "Per-season period search, P ≤ span/3",
                 "Nothing: dropped. B exists on three nights; the seasons "
                 "are 23 d and 42 d clumps in which nothing physical is "
                 "expected at 1.5–14 d and most night gaps are exactly one "
                 "day. A periodogram would only advertise emptiness.",
                 "S4", _tcrb_src("§4 Phase C step 3"), DROPPED,
                 ruling=_ruled("§4 TCrB_Monitoring DROP C3",
                               "U6", "DS", "OA", "ED")),
        ),
        _backlog_phase(
            Task("TCRB-P0-restart",
                 "Restart T CrB observations — on a second instrument",
                 "The 2026–27 series on the QHY600: a new instrument that "
                 "needs its own trace, wavelength and response solution, "
                 "spliced to 2025 only through same-night θ CrB + Vega in "
                 "hrg/lrg on ≥3 nights, with the splice offset reported "
                 "with its error. October gives at most two weeks at "
                 "airmass 1.7–2.9; the real restart is the morning "
                 "apparition from mid-January 2027.",
                 "OPS", _tcrb_src("§4 Phase 0 P0-1"), DEFERRED,
                 evidence="ops/2026-08_observatory_request.md",
                 ruling=_ruled("§4 TCrB_Monitoring DEFER C2",
                               "U6", "U10", "TE.F7", "OA")),
            Task("TCRB-C2-2026-runs",
                 "Weekly ≥2 hr B flickering runs",
                 "The first monitoring-grade flickering data RLMT would "
                 "have: ≥6 weekly ≥2 h runs. No such run is possible "
                 "before mid-January 2027, and six complete in mid-March "
                 "at the earliest — so flickering is paper 2.",
                 "OPS", _tcrb_src("§4 Phase C step 2"), DEFERRED,
                 evidence="ops/2026-08_observatory_request.md",
                 ruling=_ruled("§4 TCrB_Monitoring DEFER C2",
                               "U6", "DS", "OA", "TE", "ED", "RF")),
        ),
    ),
)


_BE = "BeStar_Grism/ANALYSIS_STRATEGY.md"


def _be_src(section: str) -> Source:
    return Source(_BE, section)


#: Ruling U9: the BeSS novelty check runs first and alone, and the gate it
#: feeds decides whether a science paper exists.  Every pipeline task waits
#: on the gate — as data, so the page cannot recommend extraction while the
#: question "are these stars even in emission?" is open.
_BE_NOVELTY_GATE = ("BE-N1-gate",)

_BE_GATE_RULE = (
    "SYNTHESIS §4 (U9): BE-S-1a-bess first and alone; if fewer than two "
    "targets are verified-active emitters the science paper stops and the "
    "standards/precision material moves to the instrument section. A "
    "database check costing days decides the venue and whether there is a "
    "science paper at all — it sat pending for seven weeks while pipeline "
    "tasks progressed.")


BESTAR_GRISM = Project(
    key="BeStar_Grism",
    title="Be-Star Grism Campaign",
    claim=(
        "RE-DRAWN by the chair's ruling on BE-N1-gate (2026-10-04): the "
        "planned ten stars held no BeSS-verified active emitter, so the "
        "paper is the 19 BeSS-verified active Be stars with ≥10 RLMT grism "
        "nights in a season (QQ Gem among them). Multi-season Hα "
        "equivalent-width monitoring of bright Be stars with a 0.5 m "
        "slitless grism — a survey blind spot, since every target "
        "saturates ZTF/ASAS-SN/ATLAS/Gaia — validated one-to-one against "
        "BeSS and anchored to stability standards. EW changes are claimed "
        "only from the standards epoch (2025-12-05) onward; earlier seasons "
        "are descriptive. Season 1 is EW only. The short-period search runs "
        "only for stars with ≥3 nights of >2 h span, with a global "
        "false-alarm probability. V/R is held until the hrg dispersion "
        "disagreement (D1) is settled by test."),
    venue="Decided by the evidence after the re-draw: AJ/PASP "
          "(methods-plus-monitoring) unless a standards-epoch event with a "
          "TESS-coincident onset survives the detection rule",
    # Unchanged in wording by the review (it is conditional on BE-N1-gate
    # either way); re-emitted because the skeleton typed "\\alpha", which
    # LaTeX reads as a line break followed by the word "alpha".
    # Retitled by seat 6 (2026-10-05, edit 7): "Multi-Semester" and
    # "Variability" claimed more than one standards season and an unknown
    # event count can support.
    paper_title=(r"H$\alpha$ Equivalent-Width Monitoring of 19 Bright Be "
                 r"Stars with a 0.5 m Slitless-Grism Telescope, Validated "
                 r"against BeSS"),
    strategy=_BE,
    ruling=_ruled("§4 BeStar_Grism", "U9", "U2", "D1"),
    decisions=(
        ("Plan review of 2026-10-03 — what now binds this paper",
         "Execute with amendments. (1) The BeSS check is first and alone, "
         "and BE-N1-gate stops the science paper below two verified-active "
         "emitters (U9). (2) Injection–recovery runs for every star, "
         "through detrending, before any periodogram is opened. (3) The "
         "grism library G-1…G-5 is inherited: fixed dispersion per (grism, "
         "mechanical epoch), measured gain, pixel-based identity, lozenge "
         "background, focus and temperature regressors. (4) The instrument "
         "has five mechanical states, not three eras, and its gains are "
         "measured. (5) V/R is held pending D1. (6) Detection claims exist "
         "only from the standards epoch; pre-standards seasons are "
         "descriptive. (7) The short tier is three qualifying stars with a "
         "global FAP; the λ Eri short search is dropped. (8) Season-2 "
         "observing and the dither test are the 2027 backlog. (9) Six "
         "figures, not twelve."),
        ("Sample re-draw of 2026-10-04 — the chair's ruling on BE-N1-gate",
         "The gate failed for the planned ten (0 BeSS-verified active "
         "emitters) and passed for the campaign as observed. The paper is "
         "re-drawn to the BeSS-verified active Be stars with ≥10 RLMT grism "
         "nights in a season, taken from the novelty report's tables; QQ "
         "Gem (HD 46264) is one of them. η Hya and θ Vir stay the "
         "standards; the stars BeSS shows constant become null tests inside "
         "the standards epoch. The camera now mounted has no grism flats: "
         "BE-S2 states the method from the archive and the October "
         "acquisition is the 2027 backlog (BE-X3). Strategy §10, Sample "
         "re-draw."),
    ),
    phases=(
        Phase("Step −1 — Feasibility gates",
              "Run before any per-target pipeline effort: they set the "
              "sample, the venue, and whether there is a science paper.",
              (
                  Task("BE-S-1a-bess", "BeSS emission-state check, all ten targets",
                       "A per-target active/inactive table with the BeSS "
                       "spectra dates that bracket our seasons, and the "
                       "verified-active count that fixes the venue.",
                       "S0c", _be_src("§4 Step −1(a)"), PENDING,
                       ruling=_ruled("§4 BE-S-1a-bess", "U9", "ED", "RF"),
                       accept="per-target table with sources on disk; the "
                              "verified-active count stated."),
                  Task("BE-N1-gate",
                       "The novelty gate: is there a science paper?",
                       "A recorded decision. Fewer than two verified-active "
                       "emitters ⇒ the science paper stops and the "
                       "standards/precision material moves to the "
                       "instrument section.",
                       "S0c", _be_src(f"{_AMENDMENTS} BE-N1-gate"), PENDING,
                       depends_on=("BE-S-1a-bess",),
                       ruling=_ruled("§4 BE-N1-gate", "U9", "ED"),
                       accept="the decision and the count behind it are "
                              "recorded, with the venue that follows."),
                  Task("BE-S-1b-lameri-injection",
                       "Injection–recovery for every sample star, through "
                       "detrending",
                       "A 90%-recovery contour per re-drawn sample star on "
                       "its real nightly timestamps, through the detrending "
                       "(free offset per mechanical epoch and the airmass, "
                       "temperature and focus regressors), in units of the "
                       "night-to-night scatter — slow tier for every star, "
                       "short tier for the stars the admission rule admits.",
                       "S0c", _be_src("§4 Step −1(b)"), PENDING,
                       ruling=_ruled("§4 BE-S-1b", "DS"),
                       accept="a contour per star is published BEFORE any "
                              "periodogram is opened, with the signed "
                              "matched-cell bias (standing rule 3)."),
                  Task("BE-S-1c-erac-geometry",
                       "Characterise the era-C FITS repackaging",
                       "The HDU layout of the post-2026-04-25 files, and "
                       "repaired geometry in the catalog — the 'NAXIS1=8' "
                       "frames are full-frame images, not 8-pixel strips.",
                       "S0e", _be_src("§3.5"), DONE,
                       evidence="docs/pipeline/s0e_geometry_fix.html"),
                  Task("BE-S-1c-hrg-bandpass",
                       "Settle the hrg O₂ B-band coverage question",
                       "Whether a telluric anchor at 6867 Å lies on the "
                       "chip for hrg. Two seats say it does — and place it "
                       "at two different pixels, which is the dispersion "
                       "disagreement D1; the answer is whichever G-1 "
                       "measures.",
                       "G", _be_src("§4 Step −1(c)"), PENDING,
                       depends_on=("G-1",),
                       ruling=_ruled("§2 D1", "D1", "OA.E1", "PH"),
                       accept="O₂-B identified on one Be-star hrg frame at "
                              "the pixel offset G-1's solution predicts."),
              )),
        Phase("Steps 0–2 — Master table, QC and calibration",
              "One master frame table keyed on the measured instrument "
              "state, then the gates every frame must pass.",
              (
                  Task("BE-S0-master-table", "Master frame table",
                       "be_frames: the grism-whitelisted frames of the "
                       "re-drawn sample and its standards, with mechanical "
                       "epoch, night, BJD_TDB and S2c verdict attached "
                       "(re-issued 2026-10-05; stage_bestar_grism held the "
                       "planned ten).",
                       "S0c", _be_src("§4 Step 0"), DONE,
                       evidence="docs/pipeline/s0c_staging.html"),
                  Task("BE-S0-cone-match",
                       "Adjudicate the blank-target rows",
                       "The measured answer: the blank-target pool is not "
                       "hidden core-ten inventory, so the §3.2 totals do not "
                       "move.",
                       "S0c", _be_src("§3.2"), DONE,
                       evidence="docs/pipeline/s0c_staging.html"),
                  Task("BE-S0-filter-identity",
                       "Measure, don't assume, which filters disperse",
                       "A per-frame direct/dispersed/indeterminate verdict "
                       "from source elongation, replacing the assumption "
                       "that a filter NAME implies a spectrum.",
                       "S2c", _be_src("§3.4"), IN_PROGRESS,
                       evidence="docs/pipeline/s2c_filter_identity.html"),
                  Task("BE-S0-dispositions",
                       "Disposition the indeterminate and direct frames, "
                       "and QQ Gem",
                       "A disposition column on the staged set for every "
                       "frame S2c did not measure dispersed, and an "
                       "explicit disposition of the orphaned QQ Gem series "
                       "(comparison pool, sample extension, or dropped "
                       "with the reason).",
                       "S2c", _be_src(f"{_AMENDMENTS} BE-S0-dispositions"),
                       PENDING, depends_on=("F-6",),
                       ruling=_ruled("§4 BeStar_Grism Disposition",
                                     "DS", "ED"),
                       accept="no staged frame enters extraction without a "
                              "disposition; QQ Gem's is recorded."),
                  Task("BE-S0-era-table",
                       "The instrument table: five mechanical states, "
                       "measured gains",
                       "One canonical table of the campaign's mechanical "
                       "states (Andor; ASI before and after the 2025 "
                       "monsoon flip; QHY night 1; QHY, with the grisms "
                       "re-seated at the swap) with a measured gain for "
                       "each — replacing the three-era table and its "
                       "header gains.",
                       "S0", _be_src(f"{_AMENDMENTS} BE-S0-era-table"),
                       PENDING, depends_on=("F-3", "F-4"),
                       ruling=_ruled("§4 BeStar_Grism Era table",
                                     "TE.F1", "DE.F1"),
                       accept="one row per mechanical state with first and "
                              "last night and a measured (or explicitly "
                              "unmeasured) gain; no calibration crosses a "
                              "boundary."),
                  Task("BE-S0-header-rescrape",
                       "Re-scrape the headers: temperature, focus, flip, "
                       "wheel, gain",
                       "Per-frame CCD-TEMP, SET-TEMP, COOLPOWR, GAIN, "
                       "OFFSET, FOCPOS, FLIPSTAT, TELPIER, FWPOS and "
                       "FWALLNAM — the catalog's camtemp is NULL for every "
                       "grism frame and its focus column is a stuck MaxIm "
                       "card.",
                       "S0", _be_src("§4 Step 0"), PENDING,
                       depends_on=("F-2",),
                       ruling=_ruled("§3 F-2", "DE.F5", "TE.F8"),
                       accept="no nulls in the grism eras where the card "
                              "exists in the header."),
                  Task("BE-S0-times", "BJD_TDB for every grism frame",
                       "Barycentric timestamps before any period analysis.",
                       "S3", _be_src("§3.3"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("BE-S1-qc-gates", "Frame QC gates",
                       "A logged pass/reject per frame, whose trace-peak "
                       "distribution IS the paper's saturation statement — "
                       "with a saturated-hot-pixel mask applied before the "
                       "trace is fitted.",
                       "G", _be_src("§4 Step 1"), PENDING,
                       depends_on=_BE_NOVELTY_GATE + ("G-2", "G-3"),
                       forbids=_BE_GATE_RULE),
                  Task("BE-S2-calibration", "Calibration method from the "
                                            "archive — no grism flats exist",
                       "The calibration each frame actually receives, stated "
                       "and applied: gain and saturation per mechanical "
                       "state from detector_params, background from the "
                       "flanking bands in place of darks, no flat-field "
                       "division (the dither test is BE-X1), and the era-B "
                       "label. The camera now mounted has no grism flat, "
                       "bias or dark in the archive; the October "
                       "acquisition is BE-X3.",
                       "OPS", _be_src("§4 Step 2"), PENDING,
                       evidence="ops/2026-08_observatory_request.md",
                       ruling=_ruled("§4 BeStar_Grism", "U10", "DE", "RF"),
                       accept="the method is stated with the archive's "
                              "calibration inventory per mechanical state; "
                              "era B is labelled lower-bound-only unless its "
                              "EW floor is measured from standards."),
              )),
        Phase("Steps 3–6 — Extraction and the calibration chain",
              "Trace, wavelength, delivered resolution, and the response "
              "chain rebuilt on verified overlaps — on the shared grism "
              "library, not a second copy of it.",
              (
                  Task("BE-S3-extraction", "Per-frame trace and optimal extraction",
                       "Horne-extracted spectra with a per-frame trace fit — "
                       "never reused between frames.",
                       "G", _be_src("§4 Step 3"), PENDING,
                       depends_on=_BE_NOVELTY_GATE + ("G-2", "G-4",
                                                      "BE-S0-dispositions"),
                       forbids=_BE_GATE_RULE),
                  Task("BE-S4-wavelength", "Wavelength calibration per "
                                           "(grism, mechanical epoch)",
                       "The fixed dispersion per grism per mechanical "
                       "epoch from the shared library plus a per-frame "
                       "zero point, with the achieved telluric-anchor rms "
                       "in Å reported.",
                       "G", _be_src("§4 Step 4"), PENDING,
                       depends_on=_BE_NOVELTY_GATE + ("G-1",),
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism Inherit",
                                     "U2", "PH.P7", "RF"),
                       accept="standard-star EW scatter reported per "
                              "mechanical epoch; dispersion stable to "
                              "< 2% within one."),
                  Task("BE-S5-resolution", "Measure the delivered resolution",
                       "The measured line-spread function per grism, epoch "
                       "and focus offset — quoted as measured, never as "
                       "expected. It decides whether the held V/R task "
                       "opens.",
                       "G", _be_src("§4 Step 5"), PENDING,
                       depends_on=_BE_NOVELTY_GATE + ("G-1", "G-5"),
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§3 G-5", "TE.F4", "PH.P8"),
                       accept="FWHM of O₂-B per night published as an LSF "
                              "table; off-nominal-focus nights flagged."),
                  Task("BE-VR-hold",
                       "V/R decomposition — held pending D1",
                       "A two-component V/R series on the strongest sample "
                       "emitters, if and only if the delivered hrg "
                       "resolution supports it; "
                       "otherwise the O₂-referenced asymmetry moment "
                       "stands in.",
                       "G", _be_src(f"{_AMENDMENTS} BE-VR-hold"), BLOCKED,
                       depends_on=("G-1",),
                       blocker="HELD by ruling D1 — neither dropped nor "
                               "kept. At ≈0.47 Å/px V/R is in scope; at "
                               "≈1.59 Å/px the strategy's own FWHM < 4 Å "
                               "rule excludes it. Clears when G-1 returns "
                               "its ≥3-line test; the task then opens on "
                               "the measured LSF or is dropped by it.",
                       ruling=_ruled("§4 BeStar_Grism HOLD",
                                     "D1", "OA.E1", "PH.P8", action="HOLD"),
                       accept="measured FWHM at Hα stated before any V/R "
                              "value is plotted."),
                  Task("BE-S6-response-chain", "Response and telluric correction",
                       "A sensitivity curve per (grism, mechanical epoch, "
                       "season) — era C on Vega, era B transferred via "
                       "η Hya/θ Vir, era A on the 2024-05-20 Vega ladder, "
                       "season 1 relative only — with the 7200 Å H₂O band "
                       "depth carried per frame as a telluric regressor: "
                       "water inside the EW window scales with PWV, not "
                       "airmass, and is as large as the claimed floor.",
                       "G", _be_src("§4 Step 6"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism BE-S6", "PH"),
                       accept="EW of θ Vir shows no residual correlation "
                              "with the per-frame 7200 Å band depth."),
              )),
        Phase("Steps 7–10 — Measurement and error calibration",
              "Equivalent widths on fixed windows, flux tiers only where a "
              "standard exists, and error floors that say which seasons "
              "can carry a claim.",
              (
                  Task("BE-S7-ew", "Equivalent widths on the paper-wide windows",
                       "Per-frame and nightly EW with the photospheric "
                       "template subtracted before any disk-radius "
                       "conversion.",
                       "G", _be_src("§4 Step 7"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE),
                  Task("BE-S8-flux-tiers", "Continuum/flux calibration — "
                                           "season 1 is EW only",
                       "A pseudo-r flux series that breaks the "
                       "EW/continuum degeneracy where a standard exists. "
                       "Season 1 has no standard, so season 1 publishes EW "
                       "only.",
                       "G", _be_src("§4 Step 8"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism BE-S8", "ED"),
                       accept="no flux, absolute or pseudo-r, is quoted for "
                              "any season-1 epoch; wherever the continuum "
                              "varies a continuum light curve is carried "
                              "(standing rule 6)."),
                  Task("BE-S9-era-crosscal", "Cross-calibration across "
                                             "mechanical epochs",
                       "One free offset per (star, band) per mechanical "
                       "epoch from overlaps, with its degrees of freedom "
                       "stated — two transfer stars constrain one number.",
                       "G", _be_src("§4 Step 9"), PENDING,
                       depends_on=_BE_NOVELTY_GATE + ("BE-S0-era-table",),
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism Era table",
                                     "TE.F1", "DS"),
                       accept="standards bridge each boundary, or the "
                              "boundary is declared unbridged."),
                  Task("BE-S10-error-floors", "Error calibration — claims "
                                              "only from the standards epoch",
                       "Per-night errors from standards from 2025-12-05 "
                       "onward, where a detection rule can exist; earlier "
                       "seasons carry intra-night lower bounds and are "
                       "descriptive, with no EW-change claims.",
                       "G", _be_src("§4 Step 10"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism BE-S10",
                                     "U9", "DS", "RF"),
                       accept="every χ²ν reported per star and epoch with "
                              "its dof (standing rule 1); no ΔEW claim "
                              "before the standards epoch."),
              )),
        Phase("Steps 11–13 — Time series, events, external",
              "Two tiers, both with published windows and closed "
              "completeness contours.",
              (
                  Task("BE-S11-timeseries", "Slow- and short-tier searches",
                       "GLS on nightly medians for every sample star, and "
                       "the 0.3–2 d search only on sample stars with ≥3 "
                       "nights of >2 h span, with a GLOBAL false-alarm "
                       "probability over the whole band and the focus and "
                       "H₂O regressors among the systematics. The λ Eri "
                       "short search is dropped.",
                       "G", _be_src("§4 Step 11"), PENDING,
                       depends_on=_BE_NOVELTY_GATE
                       + ("BE-S-1b-lameri-injection",),
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism BE-S11",
                                     "DS", "ED", "RF"),
                       accept="FAP is global (night-block bootstrap, "
                              "max-statistic); every null carries the "
                              "amplitude it would have recovered and a "
                              "predicted scale (standing rule 2)."),
                  Task("BE-S12-event-timing", "Outburst/disk-event epochs",
                       "Onset epochs to ~1–3 d from template fits — the "
                       "money figure if one is TESS-coincident.",
                       "G", _be_src("§4 Step 12"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE),
                  Task("BE-S13-external", "BeSS and TESS co-analysis",
                       "The resolution-matched one-to-one BeSS validation "
                       "plot on the re-drawn sample — the figure the paper "
                       "stands on — and the TESS sectors (MAST) of every "
                       "sample star.",
                       "G", _be_src("§4 Step 13"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE),
                  Task("BE-figures", "The six-figure set",
                       "Six figures, each mapped to a sentence of the "
                       "abstract: campaign and instrument states, the "
                       "stability standards and error floor, the BeSS "
                       "validation (the credibility figure), the main EW "
                       "curves, the event analysis, and the periodograms "
                       "with their injection–recovery contours. Cut from "
                       "twelve.",
                       "G", _be_src("§7 Figure list"), PENDING,
                       depends_on=_BE_NOVELTY_GATE,
                       forbids=_BE_GATE_RULE,
                       ruling=_ruled("§4 BeStar_Grism Figures", "ED"),
                       accept="six figures or fewer; the ≤ 250-word "
                              "abstract is approved by seat 6 before they "
                              "are built (standing rule 5)."),
                  Task("BE-draft", "Write the manuscript",
                       "manuscripts/BeStar_Grism/main.tex through the "
                       "manuscript outline.",
                       "G", _be_src("§8 Manuscript outline"), PENDING,
                       depends_on=("BE-figures",)),
              )),
        _backlog_phase(
            Task("BE-X1-dither-test",
                 "Dither test in lieu of grism flats",
                 "The same star at three detector positions through each "
                 "grism — the only flat-field test a slitless grism "
                 "admits. Needs new frames.",
                 "OPS", _be_src(f"{_AMENDMENTS} BE-X1-dither-test"),
                 DEFERRED,
                 ruling=_ruled("§4 BeStar_Grism Dither", "U9", "OA")),
            Task("BE-X2-season2-observing",
                 "Season-2 observing: restart, nightly standards",
                 "λ Eri, 69 Ori and 5 Cnc from the autumn, with a standard "
                 "on every science night (HR 1544 in the autumn, η Hya "
                 "from December). Carried in the rev. 3 observatory "
                 "request; none of it is on this paper's critical path.",
                 "OPS", _be_src(f"{_AMENDMENTS} BE-X2-season2-observing"),
                 DEFERRED,
                 ruling=_ruled("§4 BeStar_Grism DEFER", "U10", "OA")),
            Task("BE-X3-calib-acquisition",
                 "Calibration frames for the mounted camera (October)",
                 "Biases, darks and flat pairs at six or more levels to "
                 "~60 kADU in the as-found QHY configuration — the same "
                 "frames give gain and linearity — and whether the ASI "
                 "(era-B) camera survives. Era-B frames cannot be acquired "
                 "on the QHY. Carried in the rev. 3 observatory request.",
                 "OPS", _be_src(f"{_AMENDMENTS} BE-X3-calib-acquisition"),
                 DEFERRED,
                 ruling=_ruled("§4 BeStar_Grism DEFER", "U10", "OA")),
        ),
    ),
)


_SN = "SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md"


def _sn_src(section: str) -> Source:
    return Source(_SN, section)


SN2023IXF = Project(
    key="SN2023ixf_LightCurve",
    title="SN 2023ixf: +5.4 to +50 d Validation and Limits Release",
    claim=(
        "A homogeneous, single-instrument nightly gri record of SN 2023ixf "
        "from +5.4 to +50 d, tied to PS1/REFCAT2, published as residuals "
        "against the world dataset with night-to-night variability limits "
        "and the saturation matrix that defines where it is trustworthy: "
        "the photometric pipeline's external validation, not an early-time "
        "physics paper. The unconditional clean start is +5.4 d and may be "
        "AT peak, so no 'rise' is claimed and 'early' has left the title; "
        "earlier epochs appear only in the saturated-frame inventory. The "
        "slot-'6' slitless grism series was triaged and NOT PROMOTED — "
        "zero extracted spectra, no wavelength source, no contamination "
        "test — and its frames go in the release; how much of that series "
        "is a grism series at all is the measurement in the panel "
        "immediately below. The one physics upside left is a flash-phase "
        "Hα excess from the (H − '[S II]') differential colour at three "
        "epochs, behind a predicted-excess gate. The model fit is dropped."),
    venue="AJ/PASP — decided 2026-10-03 on Gate 0's pre-registered rule; "
          "the ApJ upside is closed",
    paper_title=(r"Multi-Band Photometry of SN~2023ixf in M101 from $+5.4$ "
                 r"to $+50$ Days: A Validation and Limits Release from the "
                 r"Robert L. Mutel Telescope"),
    strategy=_SN,
    ruling=_ruled("§4 SN2023ixf_LightCurve", "U7"),
    #: The claim rests on what slot '6' is, so the renderer prints S2c's
    #: live verdict for it directly beneath the claim.  That is what stopped
    #: a paragraph asserting "the 83-frame grism series" from sitting a
    #: hundred lines above a table measuring 3 of those 83 frames as direct
    #: imaging, 19 as undecided, and every dispersed verdict at 'low'
    #: strength — with nothing connecting the two statements.  The series
    #: is NOT PROMOTED now; the panel stays because the release still has
    #: to say what those frames are.
    claim_filters=("6",),
    decisions=(
        ("Plan review of 2026-10-03 — what now binds this paper",
         "Re-scoped (U7). (1) A +5.4 → +50 d validation and limits release "
         "for AJ/PASP; 'early' leaves the title. (2) The grism series is "
         "NOT PROMOTED; its frames go in the release. (3) The model "
         "consistency fit (S7) is dropped. (4) Gate 0 is re-run on a fresh "
         "manifest before anything cites it, with its 438/439 arithmetic "
         "fixed, its unnamed rule named, and every broadband frame "
         "S2c-measured. (5) The linearity audit is the key detector task. "
         "(6) Variability is judged on nightly means with a whole-night "
         "bootstrap; the intra-night periodogram is dropped. (7) The Hα "
         "curve is split: a flash-phase differential colour behind a "
         "predicted-excess gate, and an ejecta-phase curve that is dropped "
         "unless a transmission curve arrives by 2026-11-15. (8) A "
         "residuals-vs-published table precedes any variability limit. "
         "(9) Five figures, not twelve."),
    ),
    phases=(
        Phase("Gate 0 — everything is downstream of this",
              "The frozen manifest, the saturation census and the grism "
              "triage. Its verdicts stand — venue AJ/PASP, grism NOT "
              "PROMOTED — and the gate is re-run on a fresh manifest "
              "before anything cites its numbers.",
              (
                  Task("SN-G0a-manifest", "Manifest freeze with global dedup",
                       "stage_sn2023ixf_lightcurve: one row per unique frame "
                       "after (basename, jd) dedup across AND within trees.",
                       "S0c", _sn_src("§4 Step 0a"), DONE,
                       evidence="docs/pipeline/s0c_staging.html"),
                  Task("SN-S0-alias-recovery",
                       "Recover the alias-hidden template epochs",
                       "The 140-frame 2026-03-21/22 'pinwheel galaxy' epoch "
                       "and the two 2023ixf1/2 post-fade frames folded into "
                       "the working set — the deepest post-fade material.",
                       "S0", _sn_src("§3.4 templates"), DONE,
                       evidence="docs/pipeline/s0_manifest.html"),
                  Task("SN-P-timing", "Mid-exposure BJD_TDB",
                       "Barycentric timestamps for every campaign frame.",
                       "S3", _sn_src("§4 Step 8"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("SN-G0b-saturation-census",
                       "Full-campaign SN saturation census",
                       "THE saturation matrix (filter × night): pixel-level "
                       "peak ADU of the SN in every frame, screened at 80% "
                       "of the measured clip. It decides the true clean "
                       "start per band and whether Q2 survives.",
                       # Bound to SN-G0, the stage that actually produces
                       # the census, not to S2, which only supplies the
                       # ceiling it screens against.  The old binding was
                       # written when the blocker WAS S2; keeping it would
                       # leave this 'done' backed by a stage that cannot go
                       # stale when the census's own inputs move.
                       "SN-G0", _sn_src("§4 Step 0b"), DONE,
                       evidence="docs/SN2023ixf_LightCurve/sn_gate0.html"),
                  Task("SN-G0c-grism-triage",
                       "Grism triage — closed: NOT PROMOTED",
                       "The verdict on the filter-'6' frames: NOT PROMOTED. "
                       "The timebox lapsed with zero extracted spectra, no "
                       "wavelength source and no contamination test; the "
                       "frames are published in the release as they are.",
                       # Re-bound G -> SN-G0 (2026-10-03): the verdict is a
                       # row of sn_g0_verdict, which SN-G0 writes.  Bound to
                       # G, this 'done' would have gone stale whenever the
                       # T CrB extraction moved and never when Gate 0 did.
                       "SN-G0", _sn_src("§4 Step 0c"), DONE,
                       evidence="docs/SN2023ixf_LightCurve/sn_gate0.html",
                       ruling=_ruled("§4 SN2023ixf_LightCurve SN-G0c",
                                     "U7", "DS", "ED", "RF", "TE")),
                  Task("SN-venue-decision", "The venue decision — AJ/PASP",
                       "AJ/PASP, decided by the rule the strategy "
                       "pre-registered: ApJ only if Gate 0 promoted the "
                       "grism or recovered a narrowband bandpass, and it "
                       "did neither.",
                       # Re-bound S4 -> SN-G0 (2026-10-03): the decision is
                       # sn_g0_verdict's 'venue' row, not a photometry
                       # product.
                       "SN-G0", _sn_src("§2 venue call"), DONE,
                       evidence="docs/SN2023ixf_LightCurve/sn_gate0.html",
                       ruling=_ruled("§4 SN2023ixf_LightCurve "
                                     "SN-venue-decision", "U7", "ED")),
                  Task("SN-G0-rerun",
                       "Re-run Gate 0 on the fresh manifest",
                       "Gate 0 rebuilt on the corrected S0, with the "
                       "usable-broadband arithmetic fixed (clean + bounded "
                       "must equal usable, per band) and the one unnamed "
                       "exclusion rule named.",
                       "SN-G0", _sn_src(f"{_AMENDMENTS} SN-G0-rerun"),
                       PENDING, depends_on=("F-1",),
                       needs_fresh=("S0", "S0c"),
                       ruling=_ruled("§4 SN2023ixf_LightCurve Re-run Gate",
                                     "U7", "DS.F9", "ED", "RF"),
                       accept="clean + bounded = usable in every band; "
                              "every exclusion rule has a name; SN-G0 "
                              "reads FRESH."),
                  Task("SN-G0d-s2c-broadband",
                       "S2c-measure every broadband frame",
                       "A measured direct/dispersed verdict for all the "
                       "campaign's broadband frames — nine in ten were "
                       "never measured, only assumed direct from the "
                       "filter name.",
                       "S2c", _sn_src(f"{_AMENDMENTS} SN-G0d-s2c-broadband"),
                       PENDING, depends_on=("F-6",),
                       ruling=_ruled("§4 SN2023ixf_LightCurve S2c-measure",
                                     "DS"),
                       accept="zero unmeasured frames in the photometry "
                              "set."),
              )),
        Phase("Steps 1–2 — The blocking audits",
              "Nothing photometric is publishable before the filters are "
              "identified and the linearity curve exists.",
              (
                  Task("SN-S1-slot6-dispersion",
                       "Measure whether filter '6' disperses, per frame",
                       "A per-frame direct/dispersed verdict for the '6' "
                       "frames — the measurement behind the grism claim.",
                       "S2c", _sn_src("§3.2 filter table / §4 Step 0c"),
                       DONE,
                       evidence="docs/pipeline/s2c_filter_identity.html"),
                  Task("SN-S1-filter-crosswalk",
                       "Empirical broadband filter identification",
                       "A colour-term regression per filter code plus the "
                       "MaxIm→pyscope crosswalk — which also gates template "
                       "adoption, not just the appendix.",
                       "S4", _sn_src("§4 Step 1"), PENDING),
                  Task("SN-S1-narrowband-curves",
                       "Recover the narrowband transmission profiles",
                       "The actual H/O/'1' filter curves, from a "
                       "manufacturer sheet or a monochromator scan. Only "
                       "the ejecta-phase curve (SN-S6-halpha-curve) still "
                       "needs them: the flash-phase estimator needs the "
                       "filters' equivalent widths, which the zero-point "
                       "ratios already give.",
                       "OPS", _sn_src("§4 Step 1(d) / Step 6"), BLOCKED,
                       blocker="External, and dated. A colour-term "
                               "regression cannot measure a narrowband "
                               "profile; the curve has to come from MACRO "
                               "records or a member campus, and nobody "
                               "holds it yet. Clears if it arrives by "
                               "2026-11-15; if it does not, this task and "
                               "the ejecta-phase curve are dropped on "
                               "that date by the ruling already made.",
                       ruling=_ruled("§4 SN2023ixf_LightCurve S6b",
                                     "U7", "ED", "RF"),
                       accept="a transmission curve on file by "
                              "2026-11-15, or the dated drop is recorded."),
                  Task("SN-S2-linearity", "Linearity audit",
                       "The key detector task of the portfolio: deviation "
                       "vs peak ADU up to the clip from the 0.5 s/2 s "
                       "pairs, cross-checked on the twilight-flat ramp and "
                       "the 2023-06-07 Albireo set, with a gain per EGAIN "
                       "epoch and StackPro on its own curve (sub-read "
                       "clipping is invisible in the sum).",
                       # UNBLOCKED 2026-08-21: S2 re-ran and measures the
                       # High Gain clip, so the value this audit refines
                       # has a query behind it.  Pending rather than in
                       # progress: nothing has been fitted.
                       "S2", _sn_src("§4 Step 2"), PENDING,
                       ruling=_ruled("§4 SN2023ixf_LightCurve "
                                     "SN-S2-linearity", "DE.F3", "DE.F9"),
                       accept="deviation vs peak ADU published to the "
                              "clip, with the screening cap justified by "
                              "a slope < 1% or moved; K per EGAIN epoch; "
                              "StackPro-vs-High-Gain flux ratio vs peak."),
              )),
        Phase("Steps 4–6 — Calibration, photometry, the Hα excess",
              "Two photometric regimes with a published overlap test; the "
              "narrowband product differential, behind a gate.",
              (
                  Task("SN-S4-resolve", "Re-solve the unsolved broadband frames",
                       "WCS for the ~177 unsolved broadband frames.",
                       "S1b", _sn_src("§4 Step 4"), IN_PROGRESS),
                  Task("SN-S4-ensemble-cal", "REFCAT2 ensemble calibration",
                       "Per-frame zero points, one campaign colour term per "
                       "filter, nightly extinction plus a second-order "
                       "(colour × airmass) term — validated on held-out "
                       "check stars that are not ensemble members.",
                       "S4", _sn_src("§4 Step 4"), PENDING,
                       depends_on=("SN-G0-rerun",),
                       ruling=_ruled("§4 SN2023ixf_LightCurve S4", "DS"),
                       accept="χ²ν of the held-out check stars reported "
                              "per band with its dof; < 0.5 is a defect "
                              "equal to > 2 (standing rule 1)."),
                  Task("SN-S5-template-table",
                       "The template table",
                       "Every template epoch listed with camera, "
                       "mechanical epoch, FWHM and focus offset. The "
                       "2023-05-04 G/R/X templates are ~300 focus counts "
                       "out; the late ones are two other cameras at two "
                       "other pixel scales.",
                       "S0", _sn_src(f"{_AMENDMENTS} SN-S5-template-table"),
                       PENDING, depends_on=("F-3",),
                       ruling=_ruled("§4 SN2023ixf_LightCurve Template",
                                     "TE.F8"),
                       accept="the table exists; narrowband templates are "
                              "taken only from the in-focus H/O frames."),
                  Task("SN-S5-photometry", "Two-regime photometry",
                       "Aperture photometry bright, template-subtracted "
                       "forced PSF photometry faint, with the overlap "
                       "agreement published as a systematic and an "
                       "explicit scintillation term in the error model "
                       "for the 0.5–2 s exposures.",
                       "S4", _sn_src("§4 Step 5"), PENDING,
                       depends_on=("SN-S2-linearity",
                                   "SN-S5-template-table"),
                       ruling=_ruled("§4 SN2023ixf_LightCurve S5", "OA"),
                       accept="the scintillation term is matched by the "
                              "check-star rms."),
                  Task("SN-S5b-peak-epoch",
                       "Read the peak epoch from the literature",
                       "The epoch of optical maximum taken from the "
                       "published light curves. If the clean start at "
                       "+5.4 d is at peak, the paper contains no rise and "
                       "says so.",
                       "S4", _sn_src(f"{_AMENDMENTS} SN-S5b-peak-epoch"),
                       PENDING,
                       ruling=_ruled("§4 SN2023ixf_LightCurve Peak epoch",
                                     "PH"),
                       accept="peak epoch cited; no 'rise' language "
                              "survives if the clean start is at peak."),
                  Task("SN-S6-0-excess-gate",
                       "The predicted-excess gate for the flash-phase Hα",
                       "The Hα excess the flash-phase estimator SHOULD "
                       "see at +2.5, +5.4 and +6.4 d, from published line "
                       "luminosities and filter widths of ~65 Å taken from "
                       "the zero-point ratios, against the narrowband "
                       "systematic.",
                       "S4", _sn_src(f"{_AMENDMENTS} SN-S6-0-excess-gate"),
                       PENDING,
                       ruling=_ruled("§4 SN2023ixf_LightCurve "
                                     "predicted-excess gate", "PH"),
                       accept="if the predicted excess is < 3× the "
                              "narrowband systematic at every epoch for "
                              "every plausible width, S6a is demoted "
                              "now, before the work."),
                  Task("SN-S6a-flash-colour",
                       "Flash-phase Hα excess from the (H − '[S II]') "
                       "colour",
                       "The SN's (H − '1') colour minus the field-star "
                       "locus at the flash epochs, with errors, compared "
                       "with synthetic photometry of published spectra "
                       "through top-hats. A Type II at < +10 d has no "
                       "[S II] emission, so that filter is pure continuum "
                       "at the same exposures on the same stars.",
                       "S4", _sn_src(f"{_AMENDMENTS} SN-S6a-flash-colour"),
                       PENDING, depends_on=("SN-S6-0-excess-gate",),
                       forbids="SYNTHESIS §4: the flash-phase estimator "
                               "sits \"behind a predicted-excess gate\" — "
                               "the expected signal is computed before the "
                               "measurement is attempted.",
                       ruling=_ruled("§4 SN2023ixf_LightCurve S6a", "PH"),
                       accept="differential excess ± error at the three "
                              "epochs, each beside its predicted scale "
                              "(standing rule 2)."),
                  Task("SN-S6-halpha-curve",
                       "Ejecta-phase Hα-band curve (S6b)",
                       "Nightly H-band flux over the census-approved "
                       "epochs after the flash phase, forward-modelled "
                       "through the bandpass — never presented as a raw "
                       "'Hα light curve'.",
                       "S4", _sn_src("§4 Step 6"), BLOCKED,
                       depends_on=("SN-S1-narrowband-curves",),
                       blocker="Waits on the transmission curve, and on a "
                               "date: by ruling it is dropped unless the "
                               "curve arrives by 2026-11-15. Once the "
                               "broad Hα wing fills the filter, a width "
                               "from zero-point ratios is not enough — "
                               "without the measured profile this is a "
                               "methods demonstration, not a physical "
                               "curve.",
                       ruling=_ruled("§4 SN2023ixf_LightCurve S6b",
                                     "U7", "PH", "ED", "RF"),
                       accept="a curve forward-modelled through a "
                              "MEASURED profile, or the dated drop "
                              "recorded on 2026-11-15."),
              )),
        Phase("Steps 8–10 — Limits, late time, release",
              "Residuals against the world first; then limits rather than "
              "detections; the manifest as the single source of truth.",
              (
                  Task("SN-S7b-residuals-table",
                       "Residuals against the published photometry",
                       "A per-band offset and RMS against the published "
                       "world dataset, with the colour-term validity range "
                       "— the paper's validation result.",
                       "S4", _sn_src(f"{_AMENDMENTS} SN-S7b-residuals-table"),
                       PENDING, depends_on=("SN-S5-photometry",),
                       ruling=_ruled("§4 SN2023ixf_LightCurve "
                                     "residuals-vs-published", "RF"),
                       accept="the table exists before any variability "
                              "limit is quoted."),
                  Task("SN-S8-variability-limits", "Variability and bump limits",
                       "Night-to-night variability limits on nightly "
                       "means, with a whole-night bootstrap and "
                       "injection-defined limits. The intra-night "
                       "periodogram is dropped: ~6 points a night cannot "
                       "carry one.",
                       "S4", _sn_src("§4 Step 8"), PENDING,
                       depends_on=("SN-S7b-residuals-table",),
                       forbids="SYNTHESIS §4: the residuals-vs-published "
                               "table is added \"before any variability "
                               "limit\".",
                       ruling=_ruled("§4 SN2023ixf_LightCurve S8", "DS"),
                       accept="every limit carries the injected amplitude "
                              "it would have recovered and a predicted "
                              "scale (standing rule 2); no intra-night "
                              "periodogram."),
                  Task("SN-S9-late-time", "Late-time stacks and limits",
                       "Per-stack 5σ upper limits — a late-time detection is "
                       "forbidden.",
                       "S4", _sn_src("§4 Step 9"), PENDING),
                  Task("SN-S10-release", "Machine-readable release",
                       "Per-frame and nightly tables with census flags, the "
                       "saturation matrix, the un-promoted grism frames, "
                       "pipeline on GitHub, Zenodo deposit, README "
                       "regenerated from the manifest.",
                       "S4", _sn_src("§4 Step 10"), PENDING),
                  Task("SN-figures", "The five-figure set",
                       "Five figures, each mapped to a sentence of the "
                       "abstract: the saturation census, the gri light "
                       "curve, the residuals against the world dataset, "
                       "the variability limits with their injection "
                       "contours, and the flash-phase differential colour "
                       "if it passes its gate. Cut from twelve.",
                       "S4", _sn_src("§7 Figure list"), PENDING,
                       ruling=_ruled("§4 SN2023ixf_LightCurve Figures",
                                     "ED"),
                       accept="five figures or fewer; the ≤ 250-word "
                              "abstract is approved by seat 6 before they "
                              "are built (standing rule 5)."),
                  Task("SN-draft", "Write the manuscript",
                       "manuscripts/SN2023ixf_LightCurve/main.tex through "
                       "the outline — including replacing the skeleton's "
                       "sec:obs pointer at the T CrB strategy with a "
                       "shared-facility paragraph both papers cite. The "
                       "title was regenerated on 2026-10-03: 'Early' is "
                       "gone.",
                       "S4", _sn_src("§8 Manuscript outline"), PENDING,
                       depends_on=("SN-figures",)),
              )),
        _dropped_phase(
            Task("SN-S7-model-consistency", "Model consistency fits",
                 "Nothing: dropped. There are no clean data before +5.4 d "
                 "and no UV, so a consistency fit adds a referee target "
                 "and no information. One overlay on the published curves "
                 "stands in for it.",
                 "S4", _sn_src("§4 Step 7"), DROPPED,
                 ruling=_ruled("§4 SN2023ixf_LightCurve DROP S7",
                               "U7", "ED", "RF")),
        ),
    ),
)


_DW = "DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md"


def _dw_src(section: str) -> Source:
    return Source(_DW, section)


#: Ruling U1 in one probe: how many NGC 5548 slot-'6' frames S2c has NOT
#: measured as dispersed.  The seat that opened the frames found the same
#: sky lozenge and zero-order stop in "direct", "dispersed" and
#: "indeterminate" alike; F-6's background-morphology test re-issues the
#: verdicts, and this number going to zero is what closes DW-P02-slot6.
_DW_SLOT6_PROBE = Probe(
    "NGC 5548 slot-'6' frames not yet measured dispersed (F-6 re-issues "
    "them; the ruling expects zero)",
    "SELECT count(*) FROM frame_dispersion WHERE filter = '6' "
    "AND canonical_target LIKE '%5548%' "
    "AND coalesce(verdict, '') != 'dispersed'",
    "zero", kind="accept")


DWARF_AGN = Project(
    key="DwarfGalaxy_AGN_Survey",
    title="Dwarf-Galaxy Candidates: Hα Vetting",
    claim=(
        "Hα detections and non-detections for the 13 of 19 published "
        "DESI-Legacy dwarf-candidate fields that have Hα imaging, plus "
        "NGC 5238 — vetting, not discovery, restricted to candidates the "
        "per-field depth can actually reach. A detection through a ~65 Å "
        "filter is a VELOCITY statement (the candidate is nearby); a "
        "non-detection cannot tell a quenched dwarf from a background "
        "galaxy out of band and is never tabulated as a star-formation "
        "limit without that caveat. Calibrated fluxes are conditional on "
        "the filter curve John Cannon holds. NGC 5548 has left the title: "
        "filter slot '6' is a grism on every night, so the broadband AGN "
        "light curve does not exist; what remains of it is a two-day "
        "triage that either earns a paragraph or is dropped. The "
        "field-star period search and eclipse timing are dropped. "
        "Explicitly not a lag, not a period, not a distance, and not a "
        "confirmation — and nothing is stacked until the literature "
        "cross-match says the candidates are not already classified."),
    venue="AJ at best: detection/non-detection with broadband depth; "
          "calibrated Hα fluxes only if Cannon's filter curve arrives",
    # Provisional, as the strategy requires: "the title is fixed only after
    # 3.1" (the literature cross-match).  What is settled is what it may no
    # longer say.
    # Retitled by Seat 6 on 2026-10-05 (the candidates are faint
    # HI-selected galaxies, 12 of 19 background; RNAAS note).
    paper_title=(r"H$\alpha$ Imaging of Faint \ion{H}{1}-Selected Galaxies "
                 r"and NGC~5238 with a 0.5\,m Telescope: Detections and "
                 r"Limits"),
    strategy=_DW,
    ruling=_ruled("§4 DwarfGalaxy_AGN_Survey", "U1"),
    decisions=(
        ("Plan review of 2026-10-03 — what now binds this paper",
         "Re-scoped (U1). (1) NGC 5548 broadband photometry is dead: slot "
         "'6' is a grism on every night, the aperture, ensemble, "
         "variability-statistics and value-gate tasks are dropped, and "
         "those frames are not plate-solved. (2) Two short triage tasks "
         "(one-night broad-Hα extraction; zero-order photometry) each earn "
         "a paragraph or are dropped. (3) The novelty cross-match and a "
         "depth-versus-expected-L(Hα) table come before any stacking. "
         "(4) The noise model is unblocked, at the StackPro read-noise "
         "floor; sky-limited depth is no longer promised. (5) No AC4040 "
         "flat will ever exist: a night-sky superflat from the June 2023 L "
         "frames is the flat. (6) The period search and eclipse timing "
         "are dropped, and the transient search with them unless its "
         "threshold and trials are pre-declared with injection; the "
         "completeness map is kept for NGC 5238 only; W leaves NGC 5238's "
         "surface photometry. (7) Six figures, not fourteen. (8) Cannon "
         "must be contacted — by James."),
    ),
    phases=(
        Phase("Phase 0 — Gate tasks",
              "Nothing else starts until these finish, and no frame count "
              "enters the manuscript until 0.3 closes.",
              (
                  Task("DW-P0-staging", "Working set staged by reference",
                       "stage_dwarfgalaxy_agn_survey: the Dw fields, "
                       "NGC 5238 and NGC 5548 with era-matched calibration.",
                       "S0c", _dw_src("§4 Phase 0.1"), DONE,
                       evidence="docs/pipeline/s0c_staging.html"),
                  Task("DW-P01-dedup-pointing",
                       "Deduplicate and validate pointing",
                       "A per-frame disposition with pointing columns — "
                       "including the 2023-03-25 NGC 5548 tracking failure, "
                       "which unique-JD counting alone missed.",
                       "S0", _dw_src("§4 Phase 0.1"), IN_PROGRESS,
                       evidence="docs/pipeline/s0_manifest.html"),
                  Task("DW-P02-slot6-dispersion",
                       "Close slot '6' with a per-night verdict table",
                       "A per-night verdict table for the NGC 5548 "
                       "campaign, re-issued by the background-morphology "
                       "test: the measurement that decided the AGN section "
                       "is not photometry.",
                       "S2c", _dw_src("§4 Phase 0.2(b)"), IN_PROGRESS,
                       evidence="docs/pipeline/s2c_filter_identity.html",
                       depends_on=("F-6",),
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey "
                                     "DW-P02-slot6", "U1", "OA.E4", "TE.F9"),
                       accept="one verdict per night, from pixels, with "
                              "Wilson intervals from the labelled truth "
                              "set; no frame left 'indeterminate'.",
                       probes=(_DW_SLOT6_PROBE,)),
                  Task("DW-P02-filter-dossier",
                       "The filter-slot dossier, science and calibration",
                       "An appendix table mapping every slot letter to a "
                       "physical bandpass per epoch, from config logs AND a "
                       "colour-locus regression — both required, and they "
                       "must agree. The focus-offset evidence is in it: "
                       "slot '6' carried one thick element for the whole "
                       "AC4040 era; slot 'W' changed element in 2024-02.",
                       "S2c", _dw_src("§4 Phase 0.2"), BLOCKED,
                       blocker="External, and it needs James. Half the "
                               "dossier is John Cannon's: the spring-2023 "
                               "wheel configs, the Hα transmission curve, "
                               "the W flats and Dw1643+07's provenance — "
                               "and it is an authorship question as well. "
                               "He has not been contacted; an e-mail is "
                               "drafted and unsent. Clears on a dated "
                               "reply. No closure of the observatory ever "
                               "blocked this.",
                       ruling=_ruled("§6 Cannon", "ED", "RF"),
                       accept="a dated reply from Cannon is on file."),
                  Task("DW-P03-plate-solve", "Plate-solve the Dw fields and "
                                             "NGC 5238",
                       "WCS, FWHM and a pointing-jitter distribution for "
                       "the unsolved Dw-field and NGC 5238 rows. The "
                       "NGC 5548 slot-'6' frames are out of this task: "
                       "they are spectra, and spectra do not plate-solve.",
                       "S1b", _dw_src("§4 Phase 0.3"), PENDING,
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey P03",
                                     "U1", "U4", "DS", "RF"),
                       accept="every staged direct-imaging frame carries "
                              "a solve verdict; no spectrum is in the "
                              "queue."),
                  Task("DW-P04-noise-model", "StackPro/High Gain noise model",
                       "Effective gain and read noise per mode from the "
                       "flat-pair photon transfer, the sum-versus-average "
                       "contradiction about StackPro resolved, and every "
                       "depth promise re-stated at the StackPro read-noise "
                       "floor — 512 s StackPro Hα is likely read-noise "
                       "limited, not sky limited.",
                       "S2", _dw_src("§4 Phase 0.4"), PENDING,
                       depends_on=("F-4",),
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DW-P04",
                                     "U4", "DE.F6", "DE.F9"),
                       accept="depth and χ² numbers use the measured "
                              "StackPro read noise per pixel per frame; "
                              "the strategy's 'average' is corrected or "
                              "confirmed against pixels.",
                       probes=(_s2_table_probe("s2_ptc_fits"),)),
              )),
        Phase("Phase 1 — Calibration frames",
              "An executable path per band, including for the bands that "
              "will never have a flat.",
              (
                  Task("DW-P11-flats", "Night-sky superflat from the June "
                                       "2023 L frames",
                       "A night-sky superflat built from the June 2023 L "
                       "frames (~400 frames, 19 fields, two weeks) — the "
                       "only flat those data will ever have, because the "
                       "AC4040 left the telescope in March 2024 and no "
                       "flat can be taken for it now.",
                       "S0b", _dw_src("§4 Phase 1.1"), PENDING,
                       evidence="docs/pipeline/s0b_calibration_inventory.html",
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DW-P11",
                                     "OA", "TE", "DE"),
                       accept="residual large-scale structure < 0.3% of "
                              "sky."),
                  Task("DW-P12-fringe-moon", "Fringing and moonlight gradients",
                       "A measured answer per filter rather than an "
                       "assumption, plus per-moon-regime delta flats if the "
                       "Hα residuals demand them.",
                       "S4", _dw_src("§4 Phase 1.2"), PENDING),
              )),
        Phase("Phase 2 — Photometric calibration",
              "Shared across the campaigns; ZMAG is QC only.",
              (
                  Task("DW-P21-zmag-qc",
                       "Per-image catalog ZMAG is QC only",
                       "A written rule and a per-frame QC flag: PinPoint's "
                       "ZMAG carries unmodelled colour terms and is absent "
                       "for every NGC 5548 frame, so it never becomes "
                       "science calibration.",
                       "S0", _dw_src("§4 Phase 2.1"), PENDING),
                  Task("DW-P2-ensemble-zp", "REFCAT2 ensemble zero points",
                       "ZP plus a linear colour term per filter per "
                       "readout-mode epoch, with PS1 as the independent "
                       "check.",
                       "S4", _dw_src("§4 Phase 2.2"), PENDING),
                  Task("DW-P2-band-transformations",
                       "Publish the band transformations",
                       "Every depth and surface-brightness limit quoted in "
                       "an IDENTIFIED band — never as 'L'.",
                       "S4", _dw_src("§4 Phase 2.3"), BLOCKED,
                       depends_on=("DW-P02-filter-dossier",),
                       blocker="A transformation needs a band, and the "
                               "band comes from the filter dossier, which "
                               "waits on Cannon. Clears when "
                               "DW-P02-filter-dossier does."),
              )),
        Phase("Phase 3 — Dwarf fields (the lead science)",
              "Cross-match first — it fixes the counts, the framing and the "
              "title — and then the novelty gate, before anything is "
              "stacked.",
              (
                  Task("DW-P31-crossmatch", "Literature cross-match per candidate",
                       "The table that defines the remaining novelty: "
                       "ELVES/HSC status, existing Kaisin & Karachentsev Hα, "
                       "and DESI spectroscopy per candidate.",
                       "S0c", _dw_src("§4 Phase 3.1"), PENDING),
                  Task("DW-N1-novelty",
                       "The novelty gate: are the candidates already "
                       "classified?",
                       "A per-candidate verdict — informative, already "
                       "classified, or uninformative at our depth — and "
                       "the decision it forces: if most candidates already "
                       "have an Hα or spectroscopic classification, the "
                       "paper stops.",
                       "S0c", _dw_src(f"{_AMENDMENTS} DW-N1-novelty"),
                       PENDING, depends_on=("DW-P31-crossmatch",),
                       ruling=_ruled("§4 DW-N1-novelty", "ED"),
                       accept="per-candidate table with sources; the "
                              "go/stop decision recorded before any "
                              "stack exists."),
                  Task("DW-P32-qc-gates", "Uniform QC gates",
                       "A pipeline-emitted disposition table — the only "
                       "rejection authority; no hand-written night lists "
                       "anywhere.",
                       "S4", _dw_src("§4 Phase 3.2"), PENDING),
                  Task("DW-P33-stacks", "Weighted coadds per field per filter",
                       "Deep stacks with a constant per-frame sky — never a "
                       "mesh background, which eats LSB flux.",
                       "S4", _dw_src("§4 Phase 3.3"), PENDING,
                       depends_on=("DW-N1-novelty", "DW-P36-0-depth-table"),
                       forbids="SYNTHESIS §4: the literature cross-match "
                               "gate is added \"before any stacking\". "
                               "Stacking fields whose candidates are "
                               "already classified, or that the depth "
                               "cannot reach, is effort with a known "
                               "answer."),
                  Task("DW-P34-depth-detectability",
                       "Measured depth and the detectability gate",
                       "Román+2020 limits, synthetic-dwarf recovery "
                       "contours, and the per-candidate table that decides "
                       "which candidates we may say anything about.",
                       "S4", _dw_src("§4 Phase 3.4"), PENDING),
                  Task("DW-P35-sersic", "Sérsic structure and colours",
                       "imfit μ_0, r_e and ellipticity per detected "
                       "candidate, colours where R exists.",
                       "S4", _dw_src("§4 Phase 3.5"), PENDING),
                  Task("DW-P36-0-depth-table",
                       "Hα depth against expected L(Hα), before stacking",
                       "The predicted Hα limit per field (from the H-band "
                       "zero point, a ~65 Å width and the frames in hand) "
                       "against the L(Hα) each candidate would have if "
                       "star-forming — so fields where the limit cannot "
                       "separate a dIrr from a dSph are labelled "
                       "uninformative before any effort is spent on them.",
                       "S4", _dw_src(f"{_AMENDMENTS} DW-P36-0"), PENDING,
                       depends_on=("DW-P04-noise-model",),
                       ruling=_ruled("§4 DW-P36-0", "PH"),
                       accept="predicted limit vs expected L(Hα) "
                              "tabulated per candidate; uninformative "
                              "fields labelled."),
                  Task("DW-P36-halpha-limits", "Hα detections and "
                                               "non-detections",
                       "THE lead result: a continuum-subtracted Hα "
                       "detection or non-detection for each of the fields "
                       "with Hα imaging. A detection is a velocity "
                       "statement; a non-detection carries the out-of-band "
                       "caveat and is not an SFR limit. Calibrated fluxes, "
                       "and SFR at an explicit fiducial distance, are "
                       "quoted only if Cannon's transmission curve arrives.",
                       "S4", _dw_src("§4 Phase 3.6"), PENDING,
                       depends_on=("DW-P36-0-depth-table",),
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey Hα "
                                     "products", "PH", "ED", "RF"),
                       accept="every non-detection carries the "
                              "out-of-band caveat and the predicted scale "
                              "beside it (standing rule 2); no flux is "
                              "quoted without the filter curve."),
                  Task("DW-P37-photometry-mode",
                       "Aperture photometry with growth-curve corrections",
                       "Aperture + curve-of-growth photometry (photutils/"
                       "sep) as the declared mode — no PSF fitting, which at "
                       "5.4″ FWHM in uncrowded fields buys nothing.",
                       "S4", _dw_src("§4 Phase 3.7"), PENDING),
              )),
        Phase("Phase 4 — NGC 5548: a triage, not a light curve",
              "The broadband campaign does not exist — the frames are "
              "slitless spectra. Two short tasks decide whether anything "
              "of it earns a paragraph.",
              (
                  Task("DW-P4x-broad-halpha-triage",
                       "One-night broad-Hα extraction triage",
                       "One night of the slot-'6' frames extracted: the "
                       "broad Hα equivalent width and its nightly scatter, "
                       "reported before any more effort. Two days. No "
                       "lags, ever.",
                       "G", _dw_src(f"{_AMENDMENTS} DW-P4x"), PENDING,
                       depends_on=("G-1",),
                       ruling=_ruled("§4 DW-P4x", "U1", "PH", "OA"),
                       accept="EW and its nightly scatter for one night; "
                              "it earns a paragraph or the task is "
                              "dropped."),
                  Task("DW-P4y-zero-order-photometry",
                       "Zero-order differential photometry feasibility",
                       "Whether the AGN's zero-order image can be "
                       "photometered against the field stars' zero orders "
                       "in the same grism frames.",
                       "G", _dw_src(f"{_AMENDMENTS} DW-P4y"), PENDING,
                       ruling=_ruled("§4 DW-P4y", "U1", "PH"),
                       accept="comparison-star zero-order rms < 2%, or "
                              "the idea is dropped."),
                  Task("DW-P46-bjd",
                       "BJD_TDB for everything in the paper",
                       "Barycentric mid-exposure times for every frame the "
                       "paper quotes a time for — the headers carry UTC "
                       "JD, which is not what a paper may quote.",
                       "S3", _dw_src("§4 Phase 4.6"), PENDING),
              )),
        Phase("Phase 5 — NGC 5238",
              "The archive's densest time-series field: stacks, an Hα map, "
              "and the completeness map that says what its light curves "
              "could have shown.",
              (
                  Task("DW-P51-ngc5238", "NGC 5238 stacks and Hα map",
                       "Per-band deep stacks and a continuum-subtracted Hα "
                       "star-formation map of an actively interacting "
                       "dwarf. W is out of the surface photometry: it has "
                       "one 0.4 s flat.",
                       "S4", _dw_src("§4 Phase 5.1"), PENDING,
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey W dropped",
                                     "RF"),
                       accept="no W-band surface-brightness number is "
                              "quoted."),
                  Task("DW-P54-completeness", "Injection–recovery completeness map",
                       "The map that converts non-detections into a "
                       "publishable statement — for NGC 5238 only, the one "
                       "field with enough epochs to carry it.",
                       "S4", _dw_src("§4 Phase 5.4"), PENDING,
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey P54",
                                     "DS"),
                       accept="90% contour published with the signed "
                              "matched-cell bias (standing rule 3)."),
                  Task("DW-P55-detrending",
                       "Covariate detrending discipline",
                       "Frame-level regression against airmass, FWHM, x/y "
                       "drift and sky with at most two Sys-Rem-like "
                       "components, and multiband Lomb–Scargle with floating "
                       "per-band offsets when the L and R+Hα blocks merge. "
                       "NEVER GP-detrend in time before a period search — "
                       "trend and signal are fitted simultaneously.",
                       "S4", _dw_src("§4 Phase 5.5"), PENDING),
                  Task("DW-figures", "The six-figure set",
                       "Six figures, each mapped to a sentence of the "
                       "abstract: the survey footprint, the "
                       "depth/detectability gate, the candidate atlas, the "
                       "Hα results, the NGC 5238 Hα map, and its "
                       "completeness map. Cut from fourteen; the NGC 5548 "
                       "panels are gone.",
                       "S4", _dw_src("§7 Figure list"), PENDING,
                       ruling=_ruled("§4 DwarfGalaxy_AGN_Survey Figures",
                                     "ED"),
                       accept="six figures or fewer; the ≤ 250-word "
                              "abstract is approved by seat 6 before they "
                              "are built (standing rule 5)."),
                  Task("DW-draft", "Write the manuscript",
                       "manuscripts/DwarfGalaxy_AGN_Survey/main.tex through "
                       "the manuscript outline. Every title and abstract "
                       "count is drawn from the Phase 3.1 and Phase 0 "
                       "tables, never typed independently. NGC 5548 left "
                       "the title on 2026-10-03.",
                       "S4", _dw_src("§8 Manuscript outline"), PENDING,
                       depends_on=("DW-figures",)),
              )),
        _dropped_phase(
            Task("DW-P41-aperture", "Host-aware aperture and sky",
                 "Nothing: dropped. An aperture on NGC 5548 presumes "
                 "direct images, and slot '6' is a grism on every night.",
                 "S1b", _dw_src("§4 Phase 4.1"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P41",
                               "U1", "OA.E4", "TE.F9", "ED.E7")),
            Task("DW-P42-ensemble", "Differential ensemble light curve",
                 "Nothing: dropped. There is no broadband series to form "
                 "an ensemble on; the other filters total four frames on "
                 "one night.",
                 "S1b", _dw_src("§4 Phase 4.2–4.3"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P41",
                               "U1", "OA.E4", "TE.F9")),
            Task("DW-P44-statistics", "Variability statistics",
                 "Nothing: dropped. F_var and a structure function from "
                 "two or three imaging epochs are not statistics.",
                 "S4", _dw_src("§4 Phase 4.4"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P41",
                               "U1", "PH", "DS")),
            Task("DW-P45-gate", "The Phase 4.5 value gate",
                 "Nothing: dropped. The gate asked whether an RLMT light "
                 "curve adds to ZTF/ASAS-SN/ATLAS; there is no RLMT light "
                 "curve to compare.",
                 "S4", _dw_src("§4 Phase 4.5"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P41",
                               "U1", "DS", "RF")),
            Task("DW-P52-subtraction", "Image-subtraction transient search",
                 "Nothing: dropped — unless the detection threshold and "
                 "the trials count (fields × epochs × pixels) are "
                 "pre-declared and injected point sources are recovered. "
                 "A variable-star census in dwarf fields is a different "
                 "paper for a different reader.",
                 "S4", _dw_src("§4 Phase 5.2"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P53",
                               "DS", "ED")),
            Task("DW-P53-period-search", "Field-star period search",
                 "Nothing: dropped. Five to sixteen nightly epochs over "
                 "≤ 32 d cannot support a period search.",
                 "S4", _dw_src("§4 Phase 5.3"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P53",
                               "DS", "OA", "ED")),
            Task("DW-P56-eclipse-timing",
                 "Eclipse timing, only if eclipsing binaries emerge",
                 "Nothing: dropped, with the period search that would "
                 "have had to find the binaries first.",
                 "S4", _dw_src("§4 Phase 5.6"), DROPPED,
                 ruling=_ruled("§4 DwarfGalaxy_AGN_Survey DROP P53",
                               "DS", "OA", "ED")),
        ),
    ),
)


#: The legacy archive has no committee strategy, so its tasks derive from
#: the standing decisions in THIS module.  They used to cite
#: "docs/Legacy_Rigel/index.html §Status" — a heading on a page this
#: renderer generates, which the first render deleted, leaving eleven
#: citations pointing at nothing and two tasks naming the page as the
#: evidence for its own claims.  A generated page cannot be a source.
_RIG = "pipeline/macro_core/project_plan.py"


#: The decision headings the legacy tasks derive from.
_RIG_S0D = "§Multi-archive keying (stage S0d)"
_RIG_ROOT = "§Why a separate archive root"
_RIG_SCAN = "§What the first scan shows"
_RIG_U8 = "§Census only (committee ruling U8, 2026-10-03)"


def _rig_src(section: str) -> Source:
    return Source(_RIG, section)


#: THE KEY STAYS ``Legacy_Rigel``; only the displayed title changes.  The
#: key is a directory (``docs/Legacy_Rigel/``, ``Legacy_Rigel/notes/``), a
#: provenance resource, a navigation label and a link target on the public
#: site; renaming it would break every one of those to correct a word a
#: reader never sees.  The task ids keep their ``RIG-`` prefix for the same
#: reason ids are never renamed anywhere: the status history is keyed on
#: them.
LEGACY_RIGEL = Project(
    key="Legacy_Rigel",
    title="Legacy Archive 2015–2022 (census only)",
    claim=(
        "A CENSUS, not a paper — and no science commitment before it. This "
        "is the pre-MACRO Winer/Iowa archive: about 210,000 files across "
        "~1,100 nights (a file count from the transfer itself — NOT a "
        "pipeline-emitted frame count, because nothing here has reached "
        "the manifest yet). It was filed as the 'Rigel' archive on the "
        "premise that it is a different telescope with one camera; the "
        "headers say otherwise. Only the first few per cent, through "
        "spring 2015, is the Rigel system with the FLI ProLine PL16803. "
        "From late 2015 it is the 0.5 m telescope itself, with six cameras "
        "in turn — the last of which, the SBIG Aluma AC4040, is the 2023 "
        "RLMT camera. The census (cameras first, then targets) decides "
        "whether there is a paper; the go/no-go criteria are "
        "pre-registered before the census is seen and must name a "
        "question, a reader and a venue, or the outcome is a data-release "
        "note."),
    venue="none — census only (L0–L1); a pre-registered go/no-go, else a "
          "data-release note",
    strategy="",
    ruling=_ruled("§4 Legacy archive", "U8", "TE.F10"),
    decisions=(
        ("Why a separate archive root",
         "The legacy data lives in legacy-archive/, never inside "
         "rlmt-archive/. Two reasons, both about not corrupting what already "
         "works. That tree is a byte-verified mirror of the Linode "
         "'testimages' bucket (rclone check: 0 differences), and foreign "
         "files would break that equivalence permanently. And the archive "
         "spans seven cameras and two filename conventions, so merging "
         "them would let legacy frames leak into RLMT-era analyses. "
         "Consequence: the pipeline reads two roots and tags every frame "
         "with its origin."),
        ("What the first scan shows",
         "A programme dominated by contact binaries and eclipsing systems "
         "(W UMa, XY Leo, TU Boo, RW Com, CC Com, RZ Com, AW Vir, TX Cnc), "
         "plus HAT-P-12, a comet and an asteroid — only a small fraction "
         "overlaps the RLMT-era target lists. That is why this is a "
         "candidate SIXTH project rather than a supplement to the five."),
        ("Multi-archive keying (stage S0d)",
         "The catalog gains archive_root and telescope columns, and era "
         "keying on CAMERA and FOCAL LENGTH — not on a telescope name — so "
         "that a configuration is never shared by accident and is shared "
         "on purpose where it is physically the same detector: the 2022 "
         "AC4040 frames may share linearity and ceiling characterisation "
         "with RLMT eras 1–2. A legacy filename parser completes it. The "
         "full target census lands when the transfer completes."),
        ("Census only (committee ruling U8, 2026-10-03)",
         "The plan review corrected the premise and bounded the work. "
         "Premise: 'Legacy_Rigel' is a misnomer — the Rigel system and its "
         "PL16803 account for roughly 4% of the files, through spring "
         "2015; the rest is the 0.5 m with six cameras (Apogee F47, SBIG "
         "6303e, Andor Aspen CG42, Andor iKon-L 936, SBIG STXL-6303, SBIG "
         "Aluma AC4040). The displayed name is corrected; the directory "
         "key is kept so that links survive. Scope: L0–L1 only. L0 "
         "reconciles the transfer's collision list and identifies cameras "
         "and mechanical epochs before any target is counted. L1 runs the "
         "overlap query FIRST (T CrB, ST LMi, YZ Cnc, 2015–2022), then "
         "nights × filters × longest run per target, a calibration census "
         "with flat pairs, a per-season clock audit, and contact binaries "
         "with minima in three or more seasons. L2: the go/no-go criteria "
         "are pre-registered before the census is seen and must name a "
         "question, a reader and a venue; otherwise the outcome is a "
         "data-release note. It does not compete with T CrB for effort."),
    ),
    phases=(
        Phase("Phase L0 — Ingest",
              "Make the pipeline multi-archive without letting legacy "
              "frames leak into RLMT-era analyses — and find out which "
              "camera took each frame before counting anything.",
              (
                  Task("RIG-L0-transfer", "Complete the archive transfer",
                       "A byte-verified legacy-archive/ mirror, separate "
                       "from the rlmt-archive/ tree whose rclone-check "
                       "equivalence must not be broken.",
                       "S0", _rig_src(_RIG_ROOT), IN_PROGRESS),
                  Task("RIG-L0-dedup-reconcile",
                       "Reconcile the collision list",
                       "Files on disk = manifest rows + named exclusions, "
                       "with the transfer's own collision list "
                       "(legacy_manifest_BAD_collisions.csv) resolved and "
                       "the dedup rule stated.",
                       "S0", _rig_src(_RIG_U8), PENDING,
                       ruling=_ruled("§4 Legacy archive L0", "U8", "DS", "RF",
                                     action="ADD"),
                       accept="files-on-disk = rows + named exclusions; "
                              "the census states its dedup rule."),
                  Task("RIG-L0-multi-archive", "Multi-archive catalog keying",
                       "archive_root and telescope columns, and era keying "
                       "on camera + focal length, so the 2022 AC4040 "
                       "frames can share detector characterisation with "
                       "RLMT eras 1–2 and nothing else is shared by "
                       "accident.",
                       "S0", _rig_src(_RIG_S0D), IN_PROGRESS),
                  Task("RIG-L0-filename-parser", "Legacy filename parser",
                       "Target, filter and exposure recovered from the "
                       "legacy filename conventions.",
                       "S0", _rig_src(_RIG_S0D), IN_PROGRESS),
                  Task("RIG-L0-header-scan",
                       "Camera-identity header scan of the legacy tree",
                       "A catalog row per legacy frame carrying INSTRUME, "
                       "camera serial, telescope and focal length, "
                       "READOUTM, GAIN/EGAIN, SET-TEMP, CCD-TEMP, binning, "
                       "geometry and SWCREATE — a camera-identity census "
                       "before any target census.",
                       "S0", _rig_src(_RIG_U8), PENDING,
                       ruling=_ruled("§4 Legacy archive L0", "U8", "DE"),
                       accept="every frame is assigned to a camera from "
                              "its header, not from its date."),
                  Task("RIG-L0-mech-epoch",
                       "Camera and mechanical-epoch timeline",
                       "The seven cameras' timeline — first and last night "
                       "of each, with the mechanical epochs inside them.",
                       "S0", _rig_src(_RIG_U8), PENDING,
                       depends_on=("RIG-L0-header-scan",),
                       ruling=_ruled("§4 Legacy archive L0", "U8", "TE.F10",
                                     action="ADD"),
                       accept="camera timeline table with first/last "
                              "night per camera."),
              )),
        Phase("Phase L1 — Census",
              "The overlap query first — it is the only task with science "
              "leverage — then what each target could actually carry.",
              (
                  Task("RIG-L1-overlap", "Overlap with the RLMT-era targets "
                                         "— run first",
                       "T CrB, ST LMi and YZ Cnc in 2015–2022 before "
                       "anything else: the only place a legacy frame could "
                       "extend an RLMT baseline. With it, the count of "
                       "contact binaries holding minima in three or more "
                       "seasons — the one physically motivated product "
                       "visible from here.",
                       "S0", _rig_src(_RIG_U8), PENDING,
                       depends_on=("RIG-L0-header-scan",),
                       ruling=_ruled("§4 Legacy archive L1",
                                     "U8", "OA", "DS", "PH"),
                       accept="per-target overlap table, and the count of "
                              "contact binaries with minima in ≥ 3 "
                              "seasons."),
                  Task("RIG-L1-target-census", "Target census",
                       "Frames and nights per target across ~1,100 nights, "
                       "alias-merged the way S0 merges RLMT targets.",
                       "S0", _rig_src(_RIG_S0D), PENDING,
                       depends_on=("RIG-L0-dedup-reconcile",
                                   "RIG-L0-header-scan")),
                  Task("RIG-L1-night-census",
                       "Usable-series census: nights × filters × longest run",
                       "Per target, nights × filters × the longest "
                       "same-filter run — so the go/no-go is made on "
                       "usable series, not on frame counts.",
                       "S0", _rig_src(_RIG_U8), PENDING,
                       depends_on=("RIG-L0-dedup-reconcile",),
                       ruling=_ruled("§4 Legacy archive L1", "U8", "DS"),
                       accept="one row per target with nights, filters "
                              "and longest same-filter run."),
                  Task("RIG-L1-calibration-census", "Legacy calibration census",
                       "Flats, darks and biases per filter per run per "
                       "camera — and flat PAIRS per camera, so gain is "
                       "measured rather than read from a header.",
                       "S0b", _rig_src(_RIG_U8), PENDING,
                       depends_on=("RIG-L0-header-scan",),
                       ruling=_ruled("§4 Legacy archive L1",
                                     "U8", "DE", "OA"),
                       accept="per-camera table of calibration frames, "
                              "with the flat pairs available for a gain "
                              "measurement."),
                  Task("RIG-L1-clock-audit", "Per-season clock audit",
                       "The header-time convention established per "
                       "camera, and each season's clock checked against "
                       "archived contact-binary minima.",
                       "S3", _rig_src(_RIG_U8), PENDING,
                       depends_on=("RIG-L0-header-scan",),
                       ruling=_ruled("§4 Legacy archive L1",
                                     "U8", "OA", "RF", action="ADD"),
                       accept="O−C within 60 s of the literature per "
                              "season, or the offset is carried."),
              )),
        Phase("Phase L2 — Decision",
              "Whether this becomes a sixth project, on criteria written "
              "down before anyone has seen the census.",
              (
                  Task("RIG-L2-prereg",
                       "Pre-register the go/no-go criteria",
                       "The criteria written and saved BEFORE the census "
                       "is seen: at least one target with ≥ 30 nights in "
                       "one filter with calibration, or an overlap that "
                       "extends an RLMT-era baseline — and a named "
                       "question, reader and venue.",
                       "STRAT", _rig_src(_RIG_U8), PENDING,
                       ruling=_ruled("§4 Legacy archive L2",
                                     "U8", "DS", "ED", action="ADD"),
                       accept="the criteria file predates every census "
                              "table."),
                  Task("RIG-L2-gonogo", "Sixth-project go/no-go",
                       "A decision on the pre-registered criteria, with "
                       "the census as its evidence. If it cannot name a "
                       "question, a reader and a venue, the outcome is a "
                       "data-release note — and contact-binary O−C "
                       "extensions are student papers, said in advance.",
                       "STRAT", _rig_src(_RIG_U8), PENDING,
                       depends_on=("RIG-L2-prereg", "RIG-L1-overlap",
                                   "RIG-L1-night-census",
                                   "RIG-L1-calibration-census"),
                       forbids="SYNTHESIS §4 (U8): the go/no-go criteria "
                               "are \"pre-registered before the census is "
                               "seen\". A criterion chosen after the "
                               "census is a description of the census.",
                       ruling=_ruled("§4 Legacy archive L2",
                                     "U8", "DS", "ED"),
                       accept="the decision names a question, a reader "
                              "and a venue, or records 'data-release "
                              "note'."),
                  Task("RIG-L2-strategy", "Commission a committee strategy",
                       "An ANALYSIS_STRATEGY.md for this archive — only if "
                       "the go/no-go says go; otherwise a recorded decision "
                       "not to write one.",
                       "STRAT", _rig_src(_RIG_SCAN), PENDING,
                       depends_on=("RIG-L2-gonogo",)),
                  Task("RIG-L2-page", "Publish the census on this page",
                       "This page's science section, replacing the ingest "
                       "status it currently carries.",
                       "WEB", _rig_src(_RIG_S0D), PENDING),
              )),
    ),
)


#: The shared foundation — Wave 0 of the 2026-10-03 plan review (SYNTHESIS
#: §3).  Every project stands on these fifteen tasks, and none of the six
#: owns them, so they are a GROUP of their own rather than six copies.  The
#: ledger had no pipeline area to put them in: ``PROJECTS`` is the six
#: papers, and the site's navigation, its provenance resources and its
#: "across six projects" counts are all keyed on that tuple.  So the
#: foundation is held beside ``PROJECTS`` (see :func:`groups`), reachable by
#: id from every command and by ``depends_on`` from every project — which
#: is the part that matters: a project task gated on F-4 reads "waiting on
#: F-4", computed from the status table, where it used to carry a sentence
#: about destroyed tables that nobody was updating.
_FND = SYNTHESIS_2026_10_03


def _fnd(task_id: str, title: str, produces: str, stage: str, accept: str,
         findings: tuple[str, ...], **kwargs) -> Task:
    """One shared-foundation task.  Its source and its ruling are the same
    row of the synthesis table, so the id is the citation."""
    return Task(task_id, title, produces, stage,
                Source(_FND, f"§3 {task_id}"), PENDING,
                ruling=Ruling(f"§3 {task_id}", findings, action="ADD"),
                accept=accept, **kwargs)


FOUNDATION = _build(Project(
    key="Shared_Foundation",
    title="Shared foundation (Wave 0)",
    claim=(
        "Not a paper: the fifteen facility tasks every paper stands on, "
        "set by the plan review of 2026-10-03. The manifest tells the "
        "hardware truth (four cameras and at least eight mechanical "
        "states, not one camera in six modes); gain, read noise, linearity "
        "and saturation are measured per camera rather than read from a "
        "header; staleness means something again; the absolute clock is "
        "checked; and one grism calibration library replaces per-frame "
        "wavelength solutions that were not physical. Two disagreements "
        "between seats are settled here by test, not by vote: the hrg "
        "dispersion (D1, by G-1) and the Mode0 gain (D2, by F-4 and G-2)."),
    venue="none — shared pipeline evidence (docs/pipeline); the instrument "
          "paper it could become is James's decision (D4)",
    strategy=_FND,
    ruling=_ruled("§3 Shared foundation", "U2", "U4", "U5", "D1", "D2"),
    phases=(
        Phase("Wave 0 — manifest, detector, provenance, clock",
              "F-1 … F-10: what must be true of the archive and the "
              "detectors before any project's number is citable.",
              (
                  _fnd("F-1", "QHY-era dedup: reduced twins linked to raw "
                              "parents",
                       "`_calibrated` reduced twins linked to their raw "
                       "parents, so eras 79 and 82 become reduced-only "
                       "aliases and stop being counted as canonical frames.",
                       "S0",
                       "canonical 2026 frames drop by ≈27k; eras 79/82 "
                       "are reduced aliases.",
                       ("TE.F6",),
                       probes=(Probe(
                           "canonical frames that are `_calibrated` "
                           "reduced twins",
                           "SELECT count(*) FROM frames WHERE "
                           "is_canonical = 1 AND path LIKE "
                           "'%\\_calibrated.%' ESCAPE '\\'", "zero",
                           kind="accept"),)),
                  _fnd("F-2", "Header re-scrape into the manifest",
                       "INSTRUME, GAIN, OFFSET, SET-TEMP, CCD-TEMP, "
                       "COOLPOWR, FOCPOS, FLIPSTAT, TELPIER, FWPOS, "
                       "FWALLNAM and SWCREATE as columns of `frames`.",
                       "S0",
                       "the columns exist, with no nulls where the card "
                       "exists in the header.",
                       ("DE.F5", "TE.F8"),
                       probes=(Probe(
                           "`frames` columns named FOCPOS (the true focus "
                           "card; `focuspos` is a stuck MaxIm value)",
                           "SELECT count(*) FROM pragma_table_info"
                           "('frames') WHERE lower(name) = 'focpos'",
                           "positive", kind="accept"),)),
                  _fnd("F-3", "The mechanical-epoch table beneath the eras",
                       "A `mech_epoch` table — camera, rotation step "
                       "> 0.3°, flip, wheel map — beneath the header-keyed "
                       "era registry, with calibration validity bounded "
                       "by it.",
                       "S0",
                       "table emitted; a test asserts no master is "
                       "applied across a boundary.",
                       ("TE.F1",), depends_on=("F-2",),
                       probes=(Probe(
                           "`mech_epoch` tables in the manifest",
                           "SELECT count(*) FROM sqlite_master WHERE "
                           "type = 'table' AND name = 'mech_epoch'",
                           "positive", kind="accept"),)),
                  _fnd("F-4", "Flat-pair photon transfer per camera and "
                              "configuration",
                       "Measured gain and read noise per camera and "
                       "configuration in `detector_params` — one table "
                       "every consumer reads — with the alias pairs "
                       "removed from `s2_noise_pairs`. Settles D2: is "
                       "Mode0 ≈1.0 e⁻/ADU, as the flat pairs say, or the "
                       "header's 0.2467?",
                       "S2",
                       "K ± 3% and read noise for the AC4040 (both EGAIN "
                       "epochs), the iKon and the ASI; the QHY flagged "
                       "'October'.",
                       ("D2", "DE.F1", "DE.F2", "DE.F4", "DS.F7")),
                  _fnd("F-5", "One saturation and linearity policy",
                       "One cap per mode from a measured residual-vs-peak "
                       "curve, applied to native-pixel peaks for "
                       "average-binned data, with a bad-pixel mask per "
                       "camera.",
                       "S2",
                       "slope < 1% up to the adopted cap in each mode.",
                       ("DE.F3", "OA.E8"), depends_on=("F-4",)),
                  _fnd("F-6", "S2c: background morphology and a labelled "
                              "truth set",
                       "The sky-lozenge background-morphology test added "
                       "to the dispersion classifier, and a labelled "
                       "truth set of 200 frames with its confusion matrix.",
                       "S2c",
                       "NGC 5548 slot '6' re-issued; Wilson intervals "
                       "published.",
                       ("OA.E4", "DS.F8")),
                  _fnd("F-7", "Provenance: data digests apart from page "
                              "digests",
                       "A freshness signal that means something: report "
                       "pages no longer stale the data stages they "
                       "describe, the status check exits non-zero only on "
                       "data staleness, ledger blockers are computed, and "
                       "a status listing states its sync age.",
                       "S0e",
                       "at least one stage reads FRESH for a verifiable "
                       "reason.",
                       ("U5", "DS.F5", "DS.F6")),
                  _fnd("F-8", "Absolute clock from archived transits, per "
                              "era",
                       "An absolute clock check per camera era from "
                       "archived exoplanet transits, replacing a bound of "
                       "over an hour — and the one eclipse residual the "
                       "paper withheld, printed and explained.",
                       "S3",
                       "|O−C| < 120 s per era, or the offset is measured "
                       "and carried.",
                       ("OA.E6", "DS.F10", "RF.M4")),
                  _fnd("F-9", "Evidence snapshot",
                       "sha256 and row counts of every product, per "
                       "release, committed; a test mode in which an "
                       "absent product is a failure, not a skip. Whether "
                       "manuscripts are tracked is James's decision.",
                       "S0",
                       "the snapshot manifest is in the repo; the CI "
                       "mode exists.",
                       ("DS.F4",)),
                  _fnd("F-10", "Clean-tree rebuild to all-FRESH",
                       "The whole DAG rebuilt end to end from a clean "
                       "tree at a tagged commit.",
                       "WEB",
                       "zero STALE stages; the tag is recorded.",
                       ("U5", "RF.M9"),
                       depends_on=("F-1", "F-2", "F-3", "F-4", "F-5", "F-6",
                                   "F-7", "F-8", "F-9", "G-1", "G-2", "G-3",
                                   "G-4", "G-5")),
              )),
        Phase("Wave 0 — the grism library",
              "G-1 … G-5: one calibration library that T CrB, the Be "
              "stars, SN 2023ixf and NGC 5548 all inherit.",
              (
                  _fnd("G-1", "Fixed dispersion per (grism, mechanical "
                              "epoch) — resolves D1",
                       "One dispersion and sign per grism per mechanical "
                       "epoch, from hot stars on disk (Vega, θ CrB): Hα "
                       "plus at least two of telluric 6277, 6867, 7186, "
                       "7594 and Hβ. It decides D1 — hrg ≈0.47 or "
                       "≈1.59 Å/px — and with it whether profile and V/R "
                       "work is in scope. The measured LSF is then "
                       "published per grism and epoch.",
                       "G",
                       "≥ 3 lines, residual < 1 px; Hα–O₂B separation "
                       "constant to < 1–2%.",
                       ("U2", "D1", "OA.E1", "PH.P7", "TE.F3", "RF")),
                  _fnd("G-2", "Grism variance and saturation from "
                              "detector_params",
                       "The Horne variance and the saturation threshold "
                       "read from the measured detector table; the "
                       "16.3 kADU 'rail' replaced by a hot-pixel mask and "
                       "the true 65,535 clip.",
                       "G",
                       "predicted variance within 20% of the "
                       "flanking-band variance.",
                       ("D2", "DE.F1"), depends_on=("F-4",)),
                  _fnd("G-3", "Pixel-based identity gate",
                       "A spectral-fingerprint identity gate: the header "
                       "is never the sole reason a frame is rejected.",
                       "G",
                       "zero header-only rejections; the false-accept "
                       "rate is measured.",
                       ("OA.E2", "TE.F2", "RF")),
                  _fnd("G-4", "Sky-lozenge background template for lrg",
                       "A lozenge template fitted before the flanking "
                       "bands, so the sharp-edged dispersed sky stops "
                       "being over-subtracted.",
                       "G",
                       "no negative continuum; extraction-method "
                       "difference < 3%.",
                       ("OA.E7",)),
                  _fnd("G-5", "Focus and temperature regressors; the "
                              "measured LSF",
                       "Line width regressed on focus offset and "
                       "CCD-TEMP, and a line-spread function measured per "
                       "grism, epoch and focus.",
                       "G",
                       "LSF table published; off-nominal nights flagged.",
                       ("TE.F4", "DE.F5", "PH.P8"), depends_on=("F-2",)),
              )),
    ),
))


#: The six project plans, in the order the hub lists them.
PROJECTS: tuple[Project, ...] = tuple(_build(p) for p in (
    TCRB_MONITORING, CV_TIMESERIES, SN2023IXF, BESTAR_GRISM, DWARF_AGN,
    LEGACY_RIGEL))

PROJECT_BY_KEY: dict[str, Project] = {p.key: p for p in PROJECTS}


# ===========================================================================
# 4.  PURE HELPERS
# ===========================================================================

def groups() -> tuple[Project, ...]:
    """Every group in the ledger: the shared foundation, then the projects.

    A function, not a constant, so a test that swaps ``PROJECTS`` for a
    fixture is validated against the fixture and not against a tuple frozen
    at import.
    """
    return (FOUNDATION,) + tuple(PROJECTS)


def all_tasks() -> tuple[Task, ...]:
    """Every PROJECT task, in project then phase then ledger order.

    The six papers' tasks and nothing else — which is what every "N of M
    plan tasks across six projects" sentence on the site counts.  The shared
    foundation is deliberately not in here: its fifteen tasks belong to no
    paper, and folding them into a project total would move six fractions
    for work none of the six owns.  :func:`ledger_tasks` is everything.
    """
    return tuple(t for p in PROJECTS for t in p.tasks)


def ledger_tasks() -> tuple[Task, ...]:
    """Every task the ledger holds: the shared foundation plus the projects.

    This is the set a status may be recorded against and the set ``sync``
    and :func:`validate` sweep.
    """
    return tuple(t for g in groups() for t in g.tasks)


def tasks_of(project_key: str) -> tuple[Task, ...]:
    """Every task of one group (a project, or the shared foundation)."""
    return group_of(project_key).tasks


def group_of(key: str) -> Project:
    """Look a group up by key, or raise."""
    for g in groups():
        if g.key == key:
            return g
    raise PlanError(f"unknown project or group: {key!r}")


def task_by_id(task_id: str) -> Task:
    """Look one task up, or raise — a typo'd id must never write a status
    row for a task that does not exist."""
    for t in ledger_tasks():
        if t.id == task_id:
            return t
    raise PlanError(f"unknown task id: {task_id!r}")


def project_of(task_id: str) -> Project:
    """The group a task belongs to (a project, or the shared foundation)."""
    return group_of(task_by_id(task_id).project)


def overlay_statuses(tasks: Sequence[Task],
                     recorded: Mapping[str, str]) -> dict[str, str]:
    """Effective status per task id: the recorded value when one exists,
    the ledger's declared value otherwise.

    The ledger is the plan's opening position; the table is what happened.
    A recorded value always wins, INCLUDING when it moves a task backwards —
    ``sync`` writing ``redo_needed`` over a ledger ``done`` is the single
    most important thing this function has to let through.

    The map also carries the status of every task these tasks DEPEND ON,
    when that task lives elsewhere in the ledger.  A dependency may cross
    groups — T CrB's wavelength solution waits on the shared foundation's
    G-1 — and a caller that overlays one project's tasks would otherwise
    see that dependency as "no status", which :func:`unmet_dependencies`
    must read as unmet: the gate would stay shut for ever, however done
    G-1 became.  Extra keys are harmless to every counter, because counts
    iterate the TASKS, never the map.
    """
    out: dict[str, str] = {}
    for t in tasks:
        status = recorded.get(t.id, t.status)
        if status not in ALL_STATUSES:
            raise PlanError(f"task {t.id}: unknown status {status!r}")
        out[t.id] = status
    wanted = {d for t in tasks for d in t.depends_on} - set(out)
    if wanted:
        for t in ledger_tasks():
            if t.id in wanted:
                status = recorded.get(t.id, t.status)
                if status not in ALL_STATUSES:
                    raise PlanError(f"task {t.id}: unknown status {status!r}")
                out[t.id] = status
    return out


def status_counts(tasks: Sequence[Task],
                  statuses: Mapping[str, str]) -> dict[str, int]:
    """Count per status, with every status present (zeros included) so a
    caller cannot render a bar that silently omits a category."""
    counts = {s: 0 for s in ALL_STATUSES}
    for t in tasks:
        counts[statuses[t.id]] += 1
    return counts


def progress_fraction(counts: Mapping[str, int]) -> tuple[int, int]:
    """``(done, in-scope total)``.

    ``redo_needed`` counts as NOT done — that is the whole reason the state
    exists.  ``dropped`` and ``deferred`` count as NEITHER: they are not in
    the denominator at all (SYNTHESIS §0).  Leaving them in would make every
    ruling that removes impossible work look like a project falling behind;
    counting them as done would let a project finish by deleting its plan.
    The two counts this fraction leaves out are in :func:`scope_summary`,
    and a view that shows the fraction owes its reader those as well.
    """
    in_scope = sum(n for s, n in counts.items() if s not in CLOSED_BY_RULING)
    return counts.get(DONE, 0), in_scope


def scope_summary(counts: Mapping[str, int]) -> dict[str, int]:
    """The whole accounting, so no view can show a fraction and hide the
    rest: ``done`` of ``in_scope``, plus ``dropped`` (shown separately) and
    ``deferred`` (the 2027 backlog), summing to ``total``."""
    done, in_scope = progress_fraction(counts)
    dropped = counts.get(DROPPED, 0)
    deferred = counts.get(DEFERRED, 0)
    return {"done": done, "in_scope": in_scope, "dropped": dropped,
            "deferred": deferred, "total": in_scope + dropped + deferred}


def dropped_tasks(tasks: Sequence[Task],
                  statuses: Mapping[str, str]) -> tuple[Task, ...]:
    """Every task closed as ``dropped``, in plan order."""
    return tuple(t for t in tasks if statuses[t.id] == DROPPED)


def backlog(tasks: Sequence[Task],
            statuses: Mapping[str, str]) -> tuple[Task, ...]:
    """The 2027 backlog: every ``deferred`` task, in plan order."""
    return tuple(t for t in tasks if statuses[t.id] == DEFERRED)


def next_up(tasks: Sequence[Task], statuses: Mapping[str, str],
            limit: int = 3,
            include: Sequence[str] = OPEN_STATUSES,
            stage_states: Optional[Mapping[str, str]] = None,
            probe_results: Optional[Mapping[str, Sequence[ProbeResult]]]
            = None) -> tuple[Task, ...]:
    """The next actionable tasks, in plan order.

    Actionable means: in progress, needing redo, or pending.  Blocked tasks
    are deliberately excluded — a blocker is not a next action, it is a
    reason there is no next action here, and it gets its own section.
    Redo-needed sorts first: work that silently stopped being true outranks
    work that was never started.

    Redo-needed tasks are COLLAPSED to one per stage.  When an S0 rebuild
    invalidates six tasks that all rest on S3, the next action is "re-run
    S3" once, not six identical rows; listing all six would push every real
    science action off a three-item list and make this section useless
    exactly when the pipeline needs attention most.

    Tasks with UNMET DEPENDENCIES are skipped too, for the same reason
    blocked ones are.  Status alone cannot see an execution order: ST LMi
    production photometry sat here as the recommended next action while
    §5 row 1 of its own strategy said "no mixed-mode fits before" the two
    detector tasks that were blocked on destroyed tables.  A page that
    recommends what its cited strategy forbids is worse than a page with an
    empty Next-up list.  Such tasks are not lost — :func:`gated_tasks`
    returns them together with the dependencies holding them.

    When ``stage_states`` and/or ``probe_results`` are supplied, the gates
    read from the DAG and from the database (``needs_fresh``, ``probes``)
    are honoured the same way.  They are optional because evaluating them
    costs a fingerprint pass over the whole DAG; a caller that has not paid
    for that still gets the dependency gates, which cost nothing.
    """
    rank = {REDO_NEEDED: 0, IN_PROGRESS: 1, PENDING: 2}
    actionable, seen_stages = [], set()
    for t in tasks:
        status = statuses[t.id]
        if status not in rank or status not in include:
            continue
        if computed_blockers(t, statuses, stage_states, probe_results):
            continue
        if status == REDO_NEEDED:
            if t.stage in seen_stages:
                continue
            seen_stages.add(t.stage)
        actionable.append(t)
    actionable.sort(key=lambda t: rank[statuses[t.id]])
    return tuple(actionable[:limit])


def open_blockers(tasks: Sequence[Task],
                  statuses: Mapping[str, str]) -> tuple[Task, ...]:
    """Every currently blocked task, in plan order."""
    return tuple(t for t in tasks if statuses[t.id] == BLOCKED)


def derive_sync(tasks: Sequence[Task], statuses: Mapping[str, str],
                stage_states: Mapping[str, str],
                recorded_at: Optional[Mapping[str, str]] = None,
                ) -> tuple[Change, ...]:
    """Statuses that the DATABASE, or a dated RULING, has already decided.

    Exactly three rules, and no others may be added without a reason written
    here — an automatic rule that can overwrite a human judgement is how a
    tracking system starts lying in the other direction:

    1. ``done`` + stage not FRESH  -> ``redo_needed``.  The evidence this
       task's claim rests on is stale, missing or was never recorded, so the
       claim is not backed any more.
    2. ``redo_needed`` + stage FRESH -> ``done``.  The stage was re-run; the
       claim is backed again, and the plan should not keep nagging.

    3. A status RECORDED BEFORE a ruling's date, on a task whose ledger
       status that ruling set, gives way to the ledger status.  The status
       table outranks the ledger because it is what happened after the plan
       was written; a ruling is what happened after the status was recorded,
       so it outranks the table by the same argument.  Without this rule a
       committee decision to close SN-G0c as NOT PROMOTED would be silently
       overridden by the ``in_progress`` somebody recorded six weeks before
       the committee met.  It needs ``recorded_at`` (the latest stamp per
       task) and is OFF without it; a status recorded on or after the
       ruling's date is a later human judgement and is never touched.

    Rule 3 is applied first and rules 1–2 act on its result, so a task a
    ruling closes as ``done`` on a stale stage becomes ``redo_needed`` in
    one step with both reasons, not ``done`` now and ``redo_needed`` on the
    next sync.

    Apart from rule 3, ``pending`` / ``in_progress`` / ``blocked`` /
    ``dropped`` / ``deferred`` are never touched: no fingerprint can tell
    whether a person has started something, and none can reverse a ruling.
    """
    changes: list[Change] = []
    for t in tasks:
        old = statuses[t.id]
        now, reasons = old, []
        stamp = (recorded_at or {}).get(t.id)
        if (t.ruling is not None and stamp is not None
                and stamp[:10] < t.ruling.date and now != t.status):
            now = t.status
            reasons.append(
                f"ruling {t.ruling} supersedes the {old!r} recorded "
                f"{stamp}, before it")
        state = stage_states.get(t.stage)
        if state is not None:
            if now == DONE and state != pv.FRESH:
                now = REDO_NEEDED
                reasons.append(
                    f"stage {t.stage} is {state} — the evidence this task's "
                    f"'done' rests on is no longer backed")
            elif now == REDO_NEEDED and state == pv.FRESH:
                now = DONE
                reasons.append(
                    f"stage {t.stage} is FRESH again — the claim is backed")
        if now != old:
            changes.append(Change(t.id, old, now, "; ".join(reasons)))
    return tuple(changes)


#: Words that describe a citation's SHAPE rather than its content.  They
#: carry no evidence that the cited passage exists, so requiring them to
#: appear in the document would only produce false alarms.
_STRUCTURAL_WORDS = frozenset({
    "phase", "phases", "step", "steps", "sec", "section", "row", "rows",
    "rule", "and", "the", "for", "its", "per", "via",
})

#: A token worth checking: letters/digits/underscore/hyphen, three or more
#: characters.  Hyphens are INSIDE the token on purpose so "P0-3" survives
#: as one checkable handle rather than dissolving into "P0" and "3".
_TOKEN_RE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_-]{2,}")

#: The "§4" in a citation, naming a top-level numbered section.
_ANCHOR_RE = re.compile(r"§\s*(\d+)")

#: A top-level markdown heading that opens a numbered section: "## 4. ...".
_HEADING_RE = re.compile(r"^##\s+(\d+)\s*\.", re.MULTILINE)


def split_numbered_sections(text: str) -> dict[str, str]:
    """Map ``"4"`` to the text of ``## 4. …`` up to the next ``## `` heading.

    Deliberately coarse: a citation of "§3.2" is checked against the whole
    of section 3, because sub-heading conventions differ across the five
    strategies and a checker that is wrong about the structure is worse
    than one that is merely generous about the scope.
    """
    out: dict[str, str] = {}
    marks = [(m.group(1), m.start()) for m in _HEADING_RE.finditer(text)]
    for i, (num, start) in enumerate(marks):
        end = marks[i + 1][1] if i + 1 < len(marks) else len(text)
        out[num] = text[start:end]
    return out


def citation_problems(items: Sequence[tuple[str, "Source"]],
                      documents: Mapping[str, str]) -> tuple[str, ...]:
    """Every citation that does not resolve, as readable one-line problems.

    PURE: it is handed the document TEXT, so the tests drive it with strings
    and it never touches a disk.

    THE RULE, in two parts:

    1. Every ``§N`` in the citation must match a ``## N.`` heading in the
       cited document.  A citation into a section that does not exist is
       not a citation.
    2. Every CONTENTFUL word of the citation — three or more characters,
       not a structural word like "Phase" or "step", not a bare number —
       must appear somewhere inside the sections it cites.

    Part 2 is what catches an invented rule.  "§4 Phase A/B (BJD_TDB rule)"
    passes part 1 — §4 exists — and fails part 2 on ``BJD_TDB``, which is
    the whole substance of the citation and appears nowhere in the document.
    That is exactly the shape of the failure: the scaffolding of a real
    citation wrapped around a claim nobody wrote down.

    Documents that are not numbered markdown (the ledger itself, for a
    project with no committee strategy) are checked differently: the
    citation's contentful words must appear in the document verbatim.  A
    project without a strategy still may not cite a section that is gone —
    which is how eleven Legacy_Rigel tasks came to cite a "§Status" heading
    that a page regeneration had already deleted.
    """
    problems: list[str] = []
    for task_id, source in items:
        text = documents.get(source.document)
        if text is None:
            problems.append(
                f"{task_id}: cites {source.document!r}, which does not exist")
            continue

        sections = split_numbered_sections(text)
        anchors = _ANCHOR_RE.findall(source.section)
        if anchors:
            missing = [a for a in anchors if a not in sections]
            if missing:
                problems.append(
                    f"{task_id}: cites {source.document} "
                    f"§{', §'.join(missing)}, but that document has "
                    f"no such numbered section "
                    f"(it has §{', §'.join(sorted(sections, key=int))})")
                continue
            haystack = "\n".join(sections[a] for a in anchors).lower()
        else:
            haystack = text.lower()

        for token in _TOKEN_RE.findall(source.section):
            if token.lower() in _STRUCTURAL_WORDS or token.isdigit():
                continue
            if token.lower() not in haystack:
                problems.append(
                    f"{task_id}: cites {source.document} "
                    f"{source.section!r}, but {token!r} appears nowhere in "
                    f"the cited section — a citation whose substance is not "
                    f"in the document is an invented source")
    return tuple(problems)


def read_cited_documents(repo_root) -> dict[str, str]:
    """Read every document the ledger cites.  The I/O half of the check."""
    docs: dict[str, str] = {}
    names = [t.source.document for t in ledger_tasks()]
    names += [t.ruling.document for t in ledger_tasks() if t.ruling]
    names += [g.ruling.document for g in groups() if g.ruling]
    for name in names:
        if name in docs:
            continue
        path = os.path.join(str(repo_root), name)
        if os.path.exists(path):
            with open(path, encoding="utf-8", errors="replace") as fh:
                docs[name] = fh.read()
    return docs


def verify_citations(repo_root) -> tuple[str, ...]:
    """Resolve every citation in the ledger against the real documents:
    each task's source, and each ruling's place in the synthesis."""
    items = [(t.id, t.source) for t in ledger_tasks()]
    items += [(f"{t.id} (ruling)", t.ruling.source)
              for t in ledger_tasks() if t.ruling]
    items += [(f"{g.key} (ruling)", g.ruling.source)
              for g in groups() if g.ruling]
    return citation_problems(items, read_cited_documents(repo_root))


#: A finding id: a synthesis ruling (U1…U10, D1…D4), a seat's numbered
#: finding (DS.F1, RF.B2, RF.M9), or a bare seat for an un-numbered item.
_FINDING_RE = re.compile(r"^(?:(U|D)(\d+)|([A-Z]{2})(?:\.([A-Z]\d+))?)$")


def ruling_problems(items: Sequence[tuple[str, Ruling]],
                    documents: Mapping[str, str]) -> tuple[str, ...]:
    """Every finding id that does not resolve, as readable problems.

    PURE: ``documents`` maps a document name to its text — the synthesis
    under its repo-relative path, each memo under its file name.

    * ``U6`` / ``D1`` must appear in the synthesis (as ``U6`` in the rulings
      table, ``D1`` in the disagreements section).
    * ``DS.F1`` must name a known seat, and ``F1`` must appear in that
      seat's memo.  A finding nobody filed cannot justify a ruling.
    * ``ED`` (a bare seat) must name a known seat whose memo exists.

    This is the same hole :func:`citation_problems` closes for strategies,
    one level up: a ``dropped`` status is only as good as the ruling behind
    it, and a ruling that cites an invented finding is an opinion with a
    footnote.
    """
    problems: list[str] = []
    for owner, ruling in items:
        synthesis = documents.get(ruling.document)
        if synthesis is None:
            problems.append(f"{owner}: ruling cites {ruling.document!r}, "
                            f"which does not exist")
            continue
        for finding in ruling.findings:
            m = _FINDING_RE.match(finding)
            if not m:
                problems.append(f"{owner}: {finding!r} is not a finding id "
                                f"(expected U6, D1, DS.F1 or a bare seat)")
                continue
            if m.group(1):
                if not re.search(rf"\b{re.escape(finding)}\b", synthesis):
                    problems.append(f"{owner}: ruling {finding} appears "
                                    f"nowhere in {ruling.document}")
                continue
            seat, local = m.group(3), m.group(4)
            memo_name = COMMITTEE_SEATS.get(seat)
            if memo_name is None:
                problems.append(f"{owner}: {finding!r} names no committee "
                                f"seat (known: {', '.join(COMMITTEE_SEATS)})")
                continue
            memo = documents.get(memo_name)
            if memo is None:
                problems.append(f"{owner}: {finding} cites the memo "
                                f"{memo_name}, which does not exist")
            elif local and not re.search(rf"\b{re.escape(local)}\b", memo):
                problems.append(f"{owner}: {finding} — {local} appears "
                                f"nowhere in {memo_name}")
    return tuple(problems)


def read_committee_documents(repo_root) -> dict[str, str]:
    """The synthesis and the seven memos beside it.  I/O half of
    :func:`ruling_problems`."""
    docs: dict[str, str] = {}
    rulings = [t.ruling for t in ledger_tasks() if t.ruling]
    rulings += [g.ruling for g in groups() if g.ruling]
    for ruling in rulings:
        if ruling.document in docs:
            continue
        path = os.path.join(str(repo_root), ruling.document)
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8", errors="replace") as fh:
            docs[ruling.document] = fh.read()
        for memo in COMMITTEE_SEATS.values():
            memo_path = os.path.join(os.path.dirname(path), memo)
            if os.path.exists(memo_path):
                with open(memo_path, encoding="utf-8",
                          errors="replace") as fh:
                    docs[memo] = fh.read()
    return docs


def verify_rulings(repo_root) -> tuple[str, ...]:
    """Resolve every finding id in the ledger against the real memos."""
    items = [(t.id, t.ruling) for t in ledger_tasks() if t.ruling]
    items += [(g.key, g.ruling) for g in groups() if g.ruling]
    return ruling_problems(items, read_committee_documents(repo_root))


def unmet_dependencies(task: Task,
                       statuses: Mapping[str, str]) -> tuple[str, ...]:
    """Dependency ids of ``task`` that are not ``done``."""
    return tuple(d for d in task.depends_on if statuses.get(d) != DONE)


def stage_gates(task: Task, stage_states: Optional[Mapping[str, str]],
                ) -> tuple[tuple[str, str], ...]:
    """``(stage, verdict)`` for each ``needs_fresh`` stage that is not FRESH.

    A stage with no verdict in the map is reported as ``UNKNOWN`` rather
    than skipped: a gate nobody evaluated is not a gate that opened.  With
    no map at all there is nothing to report — the caller has not asked.
    """
    if stage_states is None:
        return ()
    out = []
    for key in task.needs_fresh:
        state = stage_states.get(key, "UNKNOWN")
        if state != pv.FRESH:
            out.append((key, state))
    return tuple(out)


def probe_is_met(probe: Probe, value: Optional[int]) -> bool:
    """Does this number open the gate?  ``None`` (no answer) never does."""
    if value is None:
        return False
    if probe.met_when == "zero":
        return value == 0
    if probe.met_when == "positive":
        return value > 0
    raise PlanError(f"probe {probe.label!r}: unknown met_when "
                    f"{probe.met_when!r} (expected 'zero' or 'positive')")


def computed_blockers(
        task: Task, statuses: Mapping[str, str],
        stage_states: Optional[Mapping[str, str]] = None,
        probe_results: Optional[Mapping[str, Sequence[ProbeResult]]] = None,
) -> tuple[str, ...]:
    """Why this task may not start, as sentences nobody typed.

    Three kinds of gate, each read from where the truth lives:

    * a dependency that is not ``done``        (the status table);
    * a ``needs_fresh`` stage that is not FRESH (the provenance DAG);
    * a probe whose number does not open it     (the manifest itself).

    PURE: the stage verdicts and probe values are handed in.  A gate whose
    evidence was not supplied is not reported — except dependencies, which
    need nothing but the statuses and are therefore always evaluated.
    """
    def _state(dep: str) -> str:
        return STATUS_LABEL.get(statuses.get(dep, ""), "not in the ledger")

    out = [f"waiting on {d} ({_state(d)})"
           for d in unmet_dependencies(task, statuses)]
    out += [f"needs stage {key} FRESH; it is {state}"
            for key, state in stage_gates(task, stage_states)]
    for result in (probe_results or {}).get(task.id, ()):
        if result.kind == "gate" and not result.met:
            out.append(f"database gate not met — {result}")
    return tuple(out)


def acceptance_gaps(
        task: Task,
        probe_results: Optional[Mapping[str, Sequence[ProbeResult]]],
) -> tuple[ProbeResult, ...]:
    """This task's ACCEPTANCE probes that the manifest does not yet satisfy.

    For an open task this is progress information — "27,261 reduced twins
    are still counted canonical" is how far F-1 has to go.  For a task
    recorded ``done`` it is a contradiction between the status table and
    the database, which is the exact failure the plan review found in the
    other direction (blocked, while the database said open).
    """
    return tuple(r for r in (probe_results or {}).get(task.id, ())
                 if r.kind == "accept" and not r.met)


def gated_tasks(
        tasks: Sequence[Task],
        statuses: Mapping[str, str],
) -> tuple[tuple[Task, tuple[str, ...]], ...]:
    """Open, unblocked tasks whose declared dependencies are not done yet.

    These are the ones a status-only reading calls actionable and the
    execution order calls premature.  They are surfaced rather than hidden:
    "CV-P2-stlmi is running, and the two detector tasks §5 puts in front of
    it are blocked" is a sentence a reader needs.
    """
    out = []
    for t in tasks:
        if statuses.get(t.id) not in OPEN_STATUSES:
            continue
        unmet = unmet_dependencies(t, statuses)
        if unmet:
            out.append((t, unmet))
    return tuple(out)


def validate() -> None:
    """Raise :class:`PlanError` on any structural defect in the ledger.

    Called by the tests and by the CLI before it writes anything: a ledger
    that cannot be trusted must fail loudly at the top of a command rather
    than render a plausible-looking page.
    """
    seen: set[str] = set()
    for project in groups():
        if not project.phases:
            raise PlanError(f"{project.key}: no phases")
        for task in project.tasks:
            if task.id in seen:
                raise PlanError(f"duplicate task id: {task.id}")
            seen.add(task.id)
            if task.status not in LEDGER_STATUSES:
                raise PlanError(
                    f"{task.id}: ledger status {task.status!r} is not one of "
                    f"{LEDGER_STATUSES}")
            if task.stage not in pv.STAGE_BY_KEY:
                raise PlanError(
                    f"{task.id}: stage {task.stage!r} is not in the "
                    f"provenance DAG")
            if not task.produces.strip():
                raise PlanError(f"{task.id}: no product named")
            if not task.source.document or not task.source.section:
                raise PlanError(f"{task.id}: incomplete source citation")
            if task.status == BLOCKED and not task.blocker.strip():
                raise PlanError(
                    f"{task.id}: blocked with no blocker text — a blocker "
                    f"nobody can read is not a blocker")
            if task.status == DONE and not task.evidence.strip():
                raise PlanError(
                    f"{task.id}: done with no evidence link")
            # The two closing statuses are RULINGS.  One typed without the
            # ruling behind it is how a plan gets shorter by wishing.
            if task.status in CLOSED_BY_RULING and task.ruling is None:
                raise PlanError(
                    f"{task.id}: {task.status} with no ruling — a task "
                    f"leaves the plan only by being done or by a committee "
                    f"ruling that is cited (SYNTHESIS §0)")
            # ...and they are FILED where a reader will see them as closed.
            # A dropped task left among live phases reads as live work.
            want = {DROPPED: DROPPED_PHASE, DEFERRED: BACKLOG_PHASE}
            if task.status in want and task.phase != want[task.status]:
                raise PlanError(
                    f"{task.id}: {task.status} tasks are filed under the "
                    f"phase {want[task.status]!r}, not {task.phase!r}")
            if (task.phase in want.values()
                    and task.status not in CLOSED_BY_RULING):
                raise PlanError(
                    f"{task.id}: filed under {task.phase!r} but its ledger "
                    f"status is {task.status!r}")
            if task.ruling is not None:
                verb = task.ruling.action
                if verb and verb not in RULING_ACTIONS:
                    raise PlanError(
                        f"{task.id}: ruling action {verb!r} is not one of "
                        f"{RULING_ACTIONS}")
                forced = {DROPPED: "DROP", DEFERRED: "DEFER"}
                if verb and forced.get(task.status, verb) != verb:
                    raise PlanError(
                        f"{task.id}: ruling action {verb!r} contradicts "
                        f"its status {task.status!r}")
                if verb in forced.values() and task.status not in forced:
                    raise PlanError(
                        f"{task.id}: ruling action {verb!r} on a task whose "
                        f"status is {task.status!r}")
            # A task the review added or changed, and that is still to be
            # done, states the test that closes it.
            if (task.ruling is not None
                    and task.status in (PENDING, IN_PROGRESS, BLOCKED)
                    and not task.accept.strip()):
                raise PlanError(
                    f"{task.id}: carries a ruling and is still open, but "
                    f"states no acceptance criterion")
            for key in task.needs_fresh:
                if key not in pv.STAGE_BY_KEY:
                    raise PlanError(
                        f"{task.id}: needs_fresh names {key!r}, which is "
                        f"not in the provenance DAG")
            for probe in task.probes:
                if probe.met_when not in ("zero", "positive"):
                    raise PlanError(
                        f"{task.id}: probe {probe.label!r} has met_when "
                        f"{probe.met_when!r}")
                if probe.kind not in ("gate", "accept"):
                    raise PlanError(
                        f"{task.id}: probe {probe.label!r} has kind "
                        f"{probe.kind!r} (expected 'gate' or 'accept')")
                if not probe.sql.lstrip().upper().startswith("SELECT"):
                    raise PlanError(
                        f"{task.id}: probe {probe.label!r} is not a SELECT "
                        f"— a gate may read the database, never write it")
            if project.strategy and task.source.document != project.strategy:
                raise PlanError(
                    f"{task.id}: cites {task.source.document!r} but this "
                    f"project is governed by {project.strategy!r}")
            if not project.strategy:
                # A project with no committee strategy derives its tasks
                # from its standing decisions, which live in this module.
                # It may NOT cite a page this renderer generates: the first
                # render deletes whatever heading was cited, which is how
                # eleven Legacy_Rigel tasks came to cite a "§Status" section
                # that no longer existed anywhere.
                headings = {f"§{h}" for h, _ in project.decisions}
                if task.source.section not in headings:
                    raise PlanError(
                        f"{task.id}: {project.key} has no strategy, so its "
                        f"tasks must cite one of its standing decisions "
                        f"({', '.join(sorted(headings)) or 'none declared'})"
                        f" — {task.source.section!r} is not one")
            # The circularity this forbids is a task citing the page THIS
            # RENDERER GENERATES from the task's own status —
            # docs/<Project>/index.html and nothing else (see
            # report_projects.page_path).  It was written as a ban on the
            # whole docs/<Project>/ directory, which also caught the
            # stage-built evidence pages that live there by convention
            # (cv_characterization.html, cv_catalogue_tie.html,
            # cv_external_context.html).  Those are the opposite of
            # circular: each is rendered from the product database by the
            # stage that did the work, and is exactly what a DONE task
            # should point at.  Narrowed to the generated page, so the rule
            # forbids what it always meant to forbid and no more.
            if task.evidence == f"docs/{project.key}/index.html":
                raise PlanError(
                    f"{task.id}: names its own generated project page as "
                    f"its evidence — a page cannot be the evidence for the "
                    f"claims it makes")

    # Dependencies are checked in a second pass: a task may legitimately
    # depend on one declared later in the same phase, so the ids only all
    # exist once the first pass has seen them.
    for task in ledger_tasks():
        for dep in task.depends_on:
            if dep not in seen:
                raise PlanError(
                    f"{task.id}: depends on {dep!r}, which is not a task")
            if dep == task.id:
                raise PlanError(f"{task.id}: depends on itself")
    _reject_dependency_cycles()


def _reject_dependency_cycles() -> None:
    """A cycle would make :func:`gated_tasks` report a permanent standstill
    that no amount of work could clear — better to fail at import."""
    edges = {t.id: t.depends_on for t in ledger_tasks()}
    state: dict[str, int] = {}

    def walk(node: str, trail: tuple[str, ...]) -> None:
        if state.get(node) == 2:
            return
        if state.get(node) == 1:
            cycle = " -> ".join(trail[trail.index(node):] + (node,))
            raise PlanError(f"dependency cycle: {cycle}")
        state[node] = 1
        for nxt in edges.get(node, ()):
            walk(nxt, trail + (node,))
        state[node] = 2

    for task_id in edges:
        walk(task_id, ())


#: The one line of a manuscript this module emits.
_TITLE_LINE_RE = re.compile(r"^\\title\{.*\}[ \t]*$", re.MULTILINE)


def manuscript_path(project: Project) -> str:
    """Repo-relative path of a project's manuscript source."""
    return f"manuscripts/{project.key}/main.tex"


def manuscript_title(tex: str) -> Optional[str]:
    """The argument of the first one-line ``\\title{...}`` in ``tex``."""
    match = _TITLE_LINE_RE.search(tex)
    if match is None:
        return None
    line = match.group(0).rstrip()
    return line[len("\\title{"):-1]


def with_title(tex: str, title: str) -> str:
    """``tex`` with its title line replaced, and NOTHING else touched.

    Raises if there is no one-line ``\\title{...}`` to replace: guessing
    where a title should go in somebody's manuscript is not this module's
    business.  A function replacement is used so that the backslashes in a
    LaTeX title are written as they are, not read as regex escapes.
    """
    if _TITLE_LINE_RE.search(tex) is None:
        raise PlanError("manuscript has no one-line \\title{...} to replace")
    return _TITLE_LINE_RE.sub(lambda _m: "\\title{" + title + "}", tex,
                              count=1)


# ===========================================================================
# 5.  STAGE / EVIDENCE HELPERS
# ===========================================================================

def stage_report_path(stage_key: str) -> Optional[str]:
    """The published page for a stage, if one exists.

    Derived from declared ``writes`` rather than a second hand-maintained
    map, so a renamed report cannot leave a dead link here.  Two lookups, in
    order: the stage's own writes (S0e publishes its page itself), then the
    matching report stage ``R-<key>`` (most builders write tables and a
    sibling renderer writes the page).  A stage with neither returns None
    and the caller shows an em-dash — an honest "no page yet" beats a link
    to somebody else's page that does not read this stage's tables.
    """
    def _own(key: str) -> Optional[str]:
        stage = pv.STAGE_BY_KEY.get(key)
        if stage is None:
            return None
        for w in stage.writes:
            if w.startswith("file:docs/") and w.endswith(".html"):
                return w[len("file:"):]
        return None

    return _own(stage_key) or _own(f"R-{stage_key}")


def evidence_verdict(stage_key: str, freshness: Mapping[str, pv.Freshness],
                     fingerprints: Mapping[str, str],
                     ever_ran: Iterable[str] = ()) -> tuple[str, str]:
    """``(verdict, why)`` for one stage, as an evidence page should show it.

    DESTROYED is separated out from STALE deliberately.  "Stale" invites a
    reader to think the old numbers are roughly right; "destroyed" tells
    them there is no table behind the number at all — which is the true
    state of every S1 and S2 constant since the S0 table swap.

    DESTROYED IS A CLAIM ABOUT THE PAST, so it needs evidence about the
    past.  Missing outputs alone cannot tell "this was built and then wiped"
    from "this has never been built once" — provenance itself calls both
    NEVER_RUN.  The first genuinely new stage added to the DAG would
    otherwise render on every project page in red, announcing that its
    outputs had been destroyed, which would be alarming and false on
    precisely the pages this machinery exists to keep honest.  So the label
    requires ``ever_ran`` — stages the manifest can show once produced
    something (see :func:`stages_ever_run`).  Without that evidence the
    honest answer is NEVER_RUN, with the absent outputs given as the reason.
    """
    fresh = freshness.get(stage_key)
    stage = pv.STAGE_BY_KEY[stage_key]
    gone = [w for w in stage.writes if fingerprints.get(w) == "MISSING"]
    if gone:
        names = ", ".join(sorted(g.split(":", 1)[-1] for g in gone))
        if stage_key in set(ever_ran):
            return ("DESTROYED",
                    f"it ran and its declared output(s) are now absent from "
                    f"the database: {names}")
        return (pv.NEVER_RUN,
                f"it has never been recorded as run, and its declared "
                f"output(s) are absent: {names}")
    if fresh is None:
        return ("UNKNOWN", "no freshness verdict was computed")
    why = fresh.reasons[0] if fresh.reasons else ""
    return (fresh.state, why)


def stages_ever_run(con: sqlite3.Connection) -> set[str]:
    """Stage keys the manifest can prove produced something at some point.

    Three kinds of evidence, all read out of ``stage_provenance``:

    1. the stage recorded a run of its own;
    2. its report stage ``R-<key>`` recorded a run — a report page cannot
       be rendered from tables that never existed;
    3. some other recorded stage listed one of this stage's declared
       outputs among its inputs.

    S1 and S2 have no rows of their own (the S0 table swap took those too),
    but ``R-S1`` and ``R-S2`` do, and both name the wiped tables as inputs.
    That is how the pages can still say DESTROYED about them, honestly,
    while a brand-new stage says NEVER_RUN.
    """
    if not _table_exists(con, "stage_provenance"):
        return set()
    rows = con.execute(
        "SELECT stage, inputs_json FROM stage_provenance").fetchall()
    ran = {r[0] for r in rows}
    consumed: set[str] = set()
    for _, inputs_json in rows:
        try:
            consumed |= set(json.loads(inputs_json or "{}"))
        except (ValueError, TypeError):       # pragma: no cover - defensive
            continue
    out = set()
    for stage in pv.STAGES:
        if stage.key in ran or f"R-{stage.key}" in ran:
            out.add(stage.key)
        elif any(w in consumed for w in stage.writes):
            out.add(stage.key)
    return out


def _table_exists(con: sqlite3.Connection, name: str) -> bool:
    return bool(con.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name = ?",
        (name,)).fetchone())


# ===========================================================================
# 6.  THE STATUS TABLE  (I/O)
# ===========================================================================

def utcnow() -> str:
    """ISO-8601 UTC stamp, seconds resolution — the same format the
    provenance records use, so the two can be read side by side."""
    return dt.datetime.now(dt.timezone.utc).replace(
        microsecond=0).isoformat().replace("+00:00", "Z")


#: The table is APPEND-ONLY.  A status change is an event with a time, and
#: overwriting one destroys the answer to "when did this become true?" —
#: which is exactly the question the working rhythm this table exists to
#: support keeps asking.  Current status = the row with the latest stamp.
_CREATE_SQL = f"""
CREATE TABLE IF NOT EXISTS {STATUS_TABLE} (
    task_id     TEXT NOT NULL,
    status      TEXT NOT NULL,
    evidence    TEXT NOT NULL DEFAULT '',
    note        TEXT NOT NULL DEFAULT '',
    updated_utc TEXT NOT NULL
)"""

_INDEX_SQL = (f"CREATE INDEX IF NOT EXISTS ix_plan_status_task "
              f"ON {STATUS_TABLE}(task_id, updated_utc)")


def ensure_status_table(con: sqlite3.Connection) -> None:
    """Create the status table and its index if they do not exist.

    Only the WRITE paths call this.  The readers below check for existence
    instead, because they are given READ-ONLY connections — a renderer must
    never be able to modify the database it is only inspecting, and DDL on a
    read-only connection would turn a missing table into a crash rather than
    into the honest answer "no progress has been recorded yet".
    """
    con.execute(_CREATE_SQL)
    con.execute(_INDEX_SQL)


def _status_table_exists(con: sqlite3.Connection) -> bool:
    row = con.execute("SELECT 1 FROM sqlite_master WHERE type = 'table' "
                      "AND name = ?", (STATUS_TABLE,)).fetchone()
    return row is not None


def record_status(con: sqlite3.Connection, task_id: str, status: str,
                  evidence: str = "", note: str = "",
                  when: Optional[str] = None) -> str:
    """Append one status row.  Returns the stamp written.

    Validates the task id against the ledger and the status against
    :data:`ALL_STATUSES` BEFORE writing: a status row for a task that does
    not exist is unreadable later, and an unknown status value would break
    every reader that iterates the vocabulary.
    """
    task_by_id(task_id)                       # raises on an unknown id
    if status not in ALL_STATUSES:
        raise PlanError(f"unknown status {status!r}; expected one of "
                        f"{ALL_STATUSES}")
    stamp = when or utcnow()
    ensure_status_table(con)
    con.execute(
        f"INSERT INTO {STATUS_TABLE} "
        f"(task_id, status, evidence, note, updated_utc) VALUES (?,?,?,?,?)",
        (task_id, status, evidence or "", note or "", stamp))
    return stamp


def read_statuses(con: sqlite3.Connection) -> dict[str, str]:
    """Current status per task id: the latest row per task.

    Ties on the stamp (two rows written inside the same second, which
    ``sync`` can do) break on rowid — insertion order — so replay is
    deterministic.
    """
    if not _status_table_exists(con):
        return {}
    rows = con.execute(
        f"SELECT task_id, status FROM {STATUS_TABLE} "
        f"ORDER BY updated_utc, rowid").fetchall()
    out: dict[str, str] = {}
    for task_id, status in rows:
        out[task_id] = status
    return out


def read_status_stamps(con: sqlite3.Connection) -> dict[str, str]:
    """Latest ``updated_utc`` per task id — WHEN the current status was
    recorded.  What sync rule 3 compares against a ruling's date."""
    if not _status_table_exists(con):
        return {}
    rows = con.execute(
        f"SELECT task_id, updated_utc FROM {STATUS_TABLE} "
        f"ORDER BY updated_utc, rowid").fetchall()
    return {task_id: stamp for task_id, stamp in rows}


def last_sync_utc(con: sqlite3.Connection) -> Optional[str]:
    """When ``sync`` last wrote anything, or None if it never has.

    ``sync`` stamps every row it writes with a note beginning ``sync:``.
    This is the age a status listing owes its reader: a count printed from
    statuses last reconciled six weeks ago is a count about six weeks ago —
    which is how the CV plan went on reading 34/34 after every stage behind
    it had gone stale.
    """
    if not _status_table_exists(con):
        return None
    row = con.execute(
        f"SELECT max(updated_utc) FROM {STATUS_TABLE} "
        f"WHERE note LIKE 'sync:%'").fetchone()
    return row[0] if row and row[0] else None


def run_probes(con: sqlite3.Connection, tasks: Sequence[Task],
               ) -> dict[str, tuple[ProbeResult, ...]]:
    """Evaluate every probe on ``tasks`` against the manifest.  READ-ONLY.

    A query that fails — a table that does not exist yet, a column a later
    build will add — is an UNMET gate with the error as its text, never an
    exception: "mech_epoch does not exist" is precisely the answer the F-3
    probe is there to give until F-3 lands.
    """
    out: dict[str, tuple[ProbeResult, ...]] = {}
    for t in tasks:
        results = []
        for probe in t.probes:
            try:
                row = con.execute(probe.sql).fetchone()
                value = int(row[0]) if row and row[0] is not None else None
                results.append(ProbeResult(
                    probe.label, value, probe_is_met(probe, value),
                    "" if value is not None else "the query returned nothing",
                    probe.kind))
            except sqlite3.Error as exc:
                results.append(ProbeResult(probe.label, None, False,
                                           f"not answerable yet ({exc})",
                                           probe.kind))
        if results:
            out[t.id] = tuple(results)
    return out


def read_evidence(con: sqlite3.Connection) -> dict[str, str]:
    """Latest non-empty recorded evidence per task id."""
    if not _status_table_exists(con):
        return {}
    rows = con.execute(
        f"SELECT task_id, evidence FROM {STATUS_TABLE} "
        f"ORDER BY updated_utc, rowid").fetchall()
    out: dict[str, str] = {}
    for task_id, evidence in rows:
        if evidence:
            out[task_id] = evidence
    return out


def read_history(con: sqlite3.Connection,
                 task_id: Optional[str] = None) -> list[tuple]:
    """The full status history, oldest first — the replay."""
    if not _status_table_exists(con):
        return []
    if task_id:
        return con.execute(
            f"SELECT task_id, status, evidence, note, updated_utc FROM "
            f"{STATUS_TABLE} WHERE task_id = ? ORDER BY updated_utc, rowid",
            (task_id,)).fetchall()
    return con.execute(
        f"SELECT task_id, status, evidence, note, updated_utc FROM "
        f"{STATUS_TABLE} ORDER BY updated_utc, rowid").fetchall()


# ===========================================================================
# 7.  FRESHNESS  (I/O — the bridge to provenance.py)
# ===========================================================================

def stage_freshness(con: sqlite3.Connection,
                    repo_root) -> tuple[dict[str, pv.Freshness],
                                        dict[str, str]]:
    """``(freshness by stage key, fingerprint token by resource key)``.

    This mirrors ``check_pipeline_status.evaluate``.  It is duplicated
    rather than imported because that function lives in a SCRIPT, not an
    importable package, and importing a script would execute its argparse
    module-level setup — the alternative (a `sys.path` insertion from inside
    a library) is worse.  The judgements themselves are not duplicated:
    every one comes from ``provenance``.
    """
    keys = sorted({k for s in pv.STAGES for k in s.reads + s.writes})
    fingerprints = pv.fingerprint_all(keys, con, repo_root)
    records = pv.read_records(con)
    raw: dict[str, pv.Freshness] = {}
    for stage in pv.STAGES:
        inputs = {k: fingerprints[k] for k in stage.reads}
        outputs = {k: fingerprints[k] for k in stage.writes}
        raw[stage.key] = pv.is_stale(stage, records.get(stage.key), inputs,
                                     outputs, _code_version(stage, repo_root))
    return pv.propagate_staleness(raw, pv.STAGES), fingerprints


def _code_version(stage: pv.Stage, repo_root) -> str:
    """Current code-version constant for one stage.

    Same rule as the status script: hand-authored stages carry their own
    literal; script-hosted constants are PARSED from source (never
    imported); a report stage is versioned by the stage it renders.
    """
    if stage.hand_authored:
        return stage.code_version
    if stage.version_file:
        path = os.path.join(str(repo_root), stage.version_file)
        if not os.path.exists(path):
            return f"({stage.version_symbol}: source file missing)"
        with open(path, encoding="utf-8", errors="replace") as fh:
            got = pv.read_version_constant(fh.read(), stage.version_symbol)
        return got or f"({stage.version_symbol}: unreadable)"
    versions = pv._code_versions()
    lookup = stage.key[2:] if stage.key.startswith("R-") else stage.key
    return versions.get(lookup, stage.code_version)
