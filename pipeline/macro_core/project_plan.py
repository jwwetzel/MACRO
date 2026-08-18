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
    "LEDGER_STATUSES", "ALL_STATUSES", "OPEN_STATUSES", "STATUS_LABEL",
    "PlanError",
    "Source", "Task", "Phase", "Project", "Change",
    "PROJECTS", "PROJECT_BY_KEY",
    "all_tasks", "tasks_of", "task_by_id", "project_of",
    "overlay_statuses", "status_counts", "progress_fraction",
    "next_up", "open_blockers", "derive_sync", "validate",
    "stage_report_path", "evidence_verdict", "stages_ever_run",
    "citation_problems", "verify_citations", "read_cited_documents",
    "split_numbered_sections", "unmet_dependencies", "gated_tasks",
    "ensure_status_table", "record_status", "read_statuses", "read_history",
    "stage_freshness", "utcnow",
]

#: Bumped when a change here would alter what a stored status row MEANS
#: (a new status value, a changed sync rule).  Recorded in the note of every
#: row this module writes, so a later reader can tell v1 progress from v2.
PLAN_CODE_VERSION = "PLAN v1.0 (2026-08-18)"

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

#: What a task may declare in code.
LEDGER_STATUSES: tuple[str, ...] = (DONE, IN_PROGRESS, BLOCKED, PENDING)

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
    """
    phases = tuple(
        replace(ph, tasks=tuple(
            replace(t, project=project.key, phase=ph.name) for t in ph.tasks))
        for ph in project.phases)
    return replace(project, phases=phases)


# ===========================================================================
# 3.  THE PLAN — six projects
# ===========================================================================
# Every task below was read out of the cited document's execution order.
# Where a phase has already been executed it is marked done and cites the
# evidence page; where the 2026-08 audit showed the evidence was destroyed,
# the task is marked blocked with the destruction named, or left `done` for
# `sync` to flip once the DAG confirms the stage is no longer fresh.

_CV = "CV_TimeSeries/ANALYSIS_STRATEGY.md"

#: Phase 2 production photometry may not run before these three land.  The
#: edge is real and written down; before it was DATA rather than prose, the
#: page recommended ST LMi photometry as the actionable front while all
#: three sat blocked on the destroyed S2 tables.
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
    title="Cataclysmic-Variable Time Series",
    claim=(
        "Single-night, cycle-resolved, accretion-state-tagged multi-colour "
        "light curves of three polars (ST LMi, VV Pup, EU UMa) and the SU UMa "
        "dwarf nova YZ Cnc — the per-cycle quasi-simultaneous colour coverage "
        "that phase-averaged survey folds cannot produce for 90–125 min "
        "binaries — plus bright-phase timing read as an accretion-spot "
        "longitude tracker against 40-year literature ephemerides, and "
        "YZ Cnc's dense 2024 season analysed for superhumps if AAVSO "
        "confirms the runs were in outburst. Colour science is two "
        "independent within-era analyses, never a stitched series across the "
        "2024-05 instrument seam. The full machine-readable light curves ship "
        "as a data product."),
    venue="ApJ (AASTeX 7 skeleton in manuscripts/CV_TimeSeries/main.tex)",
    strategy=_CV,
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
                       "S0c", _cv_src("§4 Phase 0 step 2"), PENDING),
                  Task("CV-P0-survey-context",
                       "Pull the survey record for all five targets",
                       "Cached ZTF/ATLAS/ASAS-SN/AAVSO/TESS/Gaia/eROSITA "
                       "series that the state histories are plotted over.",
                       "S0c", _cv_src("§4 Phase 0 step 3"), PENDING),
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
                  Task("CV-P15-linearity-ladders",
                       "Linearity ladder per readout mode",
                       "Measured linearity curves and veto thresholds for "
                       "High Gain, StackPro, Mode0 and the iKon.",
                       "S2", _cv_src("§5 row 'Linearity, per readout mode'"),
                       BLOCKED,
                       blocker="S2's evidence tables (s2_ceiling_modes, "
                               "s2_ptc_fits, s2_linearity_ladders, "
                               "detector_params) were destroyed by the S0 "
                               "table swap; the campaign must be re-run "
                               "before any ceiling or veto number has a "
                               "query behind it. Clears when S2 re-runs."),
                  Task("CV-P15-noise-model",
                       "Empirical noise model per readout mode",
                       "Check-star RMS vs magnitude per mode per night, with "
                       "the model overplotted — the non-negotiable figure.",
                       "S2", _cv_src("§5 row 'Noise model, per mode'"),
                       BLOCKED,
                       blocker="Same destroyed S2 tables; and §5 row 1 "
                               "forbids mixed-mode fits until the ladder "
                               "exists. Clears when S2 re-runs."),
              )),
        Phase("Phase 2 — Photometry",
              "Ensemble differential photometry tied to ATLAS-REFCAT2, per "
              "target, in the order §9 sets: ST LMi first (richest), then "
              "VV Pup (hardest), EU UMa, YZ Cnc.",
              (
                  Task("CV-P2-ensemble-core",
                       "Honeycutt ensemble core, proven on a prototype",
                       "A validated inhomogeneous-ensemble solver: synthetic "
                       "zero-point pattern recovered to <6 mmag with check "
                       "stars held out.",
                       "S4", _cv_src("§4 Phase 2 steps 7–9"), DONE,
                       evidence="docs/pipeline/s4_photometry.html"),
                  Task("CV-P2-vetoes",
                       "Per-mode saturation vetoes",
                       "A peak-ADU veto per readout mode applied to every "
                       "frame, with veto counts reported per target.",
                       "S2", _cv_src("§4 Phase 2 step 10"), BLOCKED,
                       blocker="The thresholds this task applies (High Gain "
                               "ceiling, StackPro clipping, Mode0 55k) have "
                               "no backing table since the S2 wipe. Clears "
                               "when S2 re-runs."),
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
                       "σ_t injection test on ST LMi 2025-02-28",
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
                       "CV-S4", _cv_src("§4 Phase 3 step 19"), BLOCKED,
                       blocker="Which branch this is depends entirely on "
                               "CV-P0-aavso-yzcnc; and the quiescent "
                               "fallback is only promised after the 8 s "
                               "High Gain S/N check passes. Clears when the "
                               "AAVSO cross-match reports."),
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
                       "CV-S4", _cv_src("§2 Q5"), PENDING),
                  Task("CV-P4-figures", "The 13-figure set",
                       "Every figure of §7, each regenerable from the "
                       "photometry product.",
                       "CV-S4", _cv_src("§7 Figure list"), PENDING),
                  Task("CV-P4-draft", "Manuscript draft",
                       "manuscripts/CV_TimeSeries/main.tex filled out "
                       "against the §8 outline.",
                       "CV-S4", _cv_src("§8 Manuscript outline"), PENDING),
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
        "Accretion diagnostics of T CrB across the pre-eruption dip and "
        "recovery, led by the 2025 slitless-grism Hα series — 247 spectra on "
        "60 nights at ~2-day cadence — with equivalent width as the "
        "season-long observable and absolute Hα flux only for the Mar–Apr "
        "2025 window where a contemporaneous θ CrB calibrator exists. "
        "Calibrated 2023–2024 photometric anchors and an archival flickering "
        "UPPER-LIMIT table support it; real flickering science needs the "
        "post-restart data, not the archive. The paper makes zero "
        "eruption-date predictions."),
    venue="ApJ, with a slitless-spectrophotometry methods appendix",
    strategy=_TCRB,
    decisions=(
        ("The strategy states no timing rule, and that is a gap, not a "
         "silence to fill",
         "Mid-exposure BJD_TDB is real, built, and published (S3) — but "
         "this project's ANALYSIS_STRATEGY.md contains no BJD, TDB, "
         "barycentric or mid-exposure rule anywhere in its 213 lines. The "
         "timing task therefore cites §9's facility-level products, which "
         "is where shared machinery is authorised, and NOT a Phase A/B "
         "timing rule, which does not exist. It was cited that way once, "
         "on this page, hyperlinked to a document that did not contain it; "
         "the citation checker now refuses that. When the committee next "
         "revises the strategy, a per-project timing rule belongs in it."),
    ),
    phases=(
        Phase("Phase 0 — Gates",
              "Nothing downstream starts until these finish; two of them "
              "are observatory actions, not analysis.",
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
                  Task("TCRB-P0-restart", "Restart T CrB observations",
                       "A nightly block (3×1 s r + short B + one lrg + one "
                       "hrg) resuming the 2025 series homogeneously.",
                       "OPS", _tcrb_src("§4 Phase 0 P0-1"), BLOCKED,
                       evidence="ops/2026-08_observatory_request.md",
                       blocker="Winer is closed for monsoon until October "
                               "2026. The request is submitted and first in "
                               "queue; clears at the re-opening."),
                  Task("TCRB-P0-filter-forensics",
                       "Map the single-character filter codes",
                       "A filter-mapping table for codes 6,1,G,H,L,O,W — 62 "
                       "frames are unusable for calibrated work until it "
                       "exists.",
                       "S0", _tcrb_src("§4 Phase 0 P0-2"), PENDING),
                  Task("TCRB-P0-bitdepth",
                       "High Gain bit-depth hardware check",
                       "The one-afternoon test that arbitrates the 12-bit vs "
                       "16-bit dispute and decides whether the archival R "
                       "series survives at all (§3 counts 121 R rawimage "
                       "rows; the staging table's working set is smaller "
                       "after the global dedup — see the table above for "
                       "the number this project actually analyses).",
                       "S2", _tcrb_src("§4 Phase 0 P0-3(i)"), BLOCKED,
                       blocker="S2's tables were destroyed by the S0 table "
                               "swap, so there is no measured ceiling to "
                               "compare a check against. Clears when S2 "
                               "re-runs (and the afternoon test happens)."),
                  Task("TCRB-P0-ladders",
                       "Exposure ladders and the 70% cap",
                       "A measured ceiling per readout mode and a hard "
                       "photometry cap at 70% of it.",
                       "S2", _tcrb_src("§4 Phase 0 P0-3(ii)"), BLOCKED,
                       blocker="Same destroyed S2 evidence tables. Clears "
                               "when S2 re-runs."),
                  Task("TCRB-P0-calib-acquisition",
                       "Acquire Mode0 240 s darks/biases and current flats",
                       "The facility calibration library the amended "
                       "cross-mode rule (ruling 7) depends on.",
                       "OPS", _tcrb_src("§4 Phase 0 P0-3(iii)"), BLOCKED,
                       evidence="ops/2026-08_observatory_request.md",
                       blocker="Monsoon closure: acquisition happens at the "
                               "October re-opening, Item B of the submitted "
                               "request."),
                  Task("TCRB-P0-shutter-timing",
                       "Shutter/exposure-timing uniformity at short exposures",
                       "Verified 0.1 s timing before the eruption-response "
                       "plan relies on it — tested on the archival 0.085 s "
                       "θ CrB frames.",
                       "S2", _tcrb_src("§4 Phase 0 P0-3(iv)"), BLOCKED,
                       blocker="Rides the same S2 campaign. Clears when S2 "
                               "re-runs."),
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
              "Nine steps from per-frame identity to the cross-validated "
              "EW series; no frame that fails the gate enters anything.",
              (
                  Task("TCRB-A0-identity-gate",
                       "Per-frame Gaia identity gate on the zero order",
                       "A pass/fail per frame with the discard fraction "
                       "published — 21 of 247 header pointings are >1° off.",
                       "G", _tcrb_src("§4 Phase A step 0"), IN_PROGRESS),
                  Task("TCRB-A1-calibrator-characterization",
                       "Characterise θ CrB before using θ CrB",
                       "Its own Hα/Hβ profile and stability across all 412 "
                       "frames — it is a Be/shell star, not a clean standard.",
                       "G", _tcrb_src("§4 Phase A step 1"), PENDING),
                  Task("TCRB-A2-extraction",
                       "Trace and optimal extraction",
                       "Boxcar and Horne-extracted 1-D spectra per frame, "
                       "with flanking-band background doubling as the Mode0 "
                       "dark-structure remover.",
                       "G", _tcrb_src("§4 Phase A step 2"), IN_PROGRESS),
                  Task("TCRB-A3-wavelength",
                       "Per-frame self-anchored wavelength solution",
                       "A zero point from zero-order geometry per frame, "
                       "with the km/s uncertainty quoted beside every "
                       "velocity. Never a global solution.",
                       "G", _tcrb_src("§4 Phase A step 3"), PENDING),
                  Task("TCRB-A4-response",
                       "Instrument response from θ CrB, Hα interpolated over",
                       "An absolute-flux calibration for Mar–Apr 2025 only; "
                       "the other 32 of 60 nights are EW-only, stated.",
                       "G", _tcrb_src("§4 Phase A step 4"), PENDING),
                  Task("TCRB-A5-ew",
                       "Hα equivalent width on Munari's convention",
                       "THE primary observable: a season-long EW series that "
                       "splices onto the published curves.",
                       "G", _tcrb_src("§4 Phase A step 5"), PENDING),
                  Task("TCRB-A6-saturation-triage",
                       "Scripted saturation triage of every grism frame",
                       "Flagged and discarded frames with the fraction "
                       "reported.",
                       "G", _tcrb_src("§4 Phase A step 6"), PENDING),
                  Task("TCRB-A7-error-floor",
                       "Rebuild the empirical error floor",
                       "θ CrB night-to-night EW scatter from the matched "
                       "2025 Mode0 subset, split by (exposure, mode), plus a "
                       "separate 240 s smear term.",
                       "G", _tcrb_src("§4 Phase A step 7"), PENDING),
                  Task("TCRB-A8-cross-validation",
                       "ARAS/published EW cross-validation",
                       "Agreement within 10–15% or an explanation — and the "
                       "spectra-per-month table that has to exist before the "
                       "'densest record' claim can appear.",
                       "G", _tcrb_src("§4 Phase A step 8"), PENDING),
              )),
        Phase("Phase B — Photometry (supporting)",
              "The 2023–2024 anchors, calibrated honestly under the amended "
              "cross-mode rule.",
              (
                  Task("TCRB-B1-calibration",
                       "Calibrate per readout mode with measured penalties",
                       "Mode-matched masters where they exist; a MEASURED "
                       "cross-mode penalty term where they never will.",
                       "S2", _tcrb_src("§4 Phase B step 1 (ruling 7)"),
                       BLOCKED,
                       blocker="The penalty term is measured against S2's "
                               "calibration library, whose tables are gone. "
                               "Clears when S2 re-runs."),
                  Task("TCRB-B2-aperture", "Aperture photometry with growth curves",
                       "Per-frame photometry at r = 1.5×FWHM, cross-checked "
                       "against fixed and wide apertures.",
                       "S4", _tcrb_src("§4 Phase B step 2"), PENDING),
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
                       "Hold the precision budget",
                       "The adopted budget stated and checked against what "
                       "was achieved: 5–8 mmag per frame, 3–5 mmag nightly "
                       "means, 5–10 mmag season-level systematic floor — and "
                       "no promise of mmag-class work anywhere in the paper.",
                       "S4", _tcrb_src("§4 Phase B step 6"), PENDING,
                       depends_on=("TCRB-B5-errors",)),
              )),
        Phase("Phase C — Time series",
              "What the archive can and cannot bound, said out loud.",
              (
                  Task("TCRB-C1-flickering-limits",
                       "Archival per-snippet flickering upper limits",
                       "A per-snippet detrended rms and 95% σ_flick upper "
                       "limit table — the honest ceiling of an 11-minute "
                       "StackPro snippet.",
                       "S4", _tcrb_src("§4 Phase C step 1"), PENDING),
                  Task("TCRB-C2-2026-runs",
                       "Weekly ≥2 hr B flickering runs",
                       "The first monitoring-grade flickering data RLMT has; "
                       "the results subsection exists only if ≥6 good runs "
                       "land.",
                       "OPS", _tcrb_src("§4 Phase C step 2"), BLOCKED,
                       evidence="ops/2026-08_observatory_request.md",
                       blocker="No dark-sky ≥2 hr window until the field "
                               "re-emerges after the October re-opening; the "
                               "runs accrue mainly in the 2027 season."),
                  Task("TCRB-C3-period-search",
                       "Per-season period search, P ≤ span/3",
                       "Periodograms with published windows, bootstrap FAPs "
                       "and injection–recovery contours — or an honest "
                       "emptiness.",
                       "S4", _tcrb_src("§4 Phase C step 3"), PENDING),
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
                  Task("TCRB-D3-release",
                       "Release the reduced 1-D spectra",
                       "A Zenodo deposit plus a submission to the ARAS T CrB "
                       "database.",
                       "G", _tcrb_src("§4 Phase D"), PENDING),
                  Task("TCRB-D4-figures", "The full ten-figure set",
                       "Figures 1–10: the three-band coverage map, the "
                       "filter/detector forensics, photometric performance, "
                       "the 2023–2024 light curve, the flickering-limits "
                       "panel, injection–recovery contours, the grism atlas, "
                       "the centrepiece, profile evolution, and the θ CrB "
                       "repeatability appendix figure.",
                       "G", _tcrb_src("§7 Figure list"), PENDING),
                  Task("TCRB-D5-draft", "Write the manuscript",
                       "manuscripts/TCrB_Monitoring/main.tex through the "
                       "§8 outline — and its FIRST action, which the outline "
                       "states outright: change the title, because the "
                       "skeleton's 'Three-Year Quiescent Baseline' claims "
                       "exactly what the panel rejected.",
                       "G", _tcrb_src("§8 Manuscript outline"), PENDING,
                       depends_on=("TCRB-D4-figures",)),
              )),
    ),
)


_BE = "BeStar_Grism/ANALYSIS_STRATEGY.md"


def _be_src(section: str) -> Source:
    return Source(_BE, section)


BESTAR_GRISM = Project(
    key="BeStar_Grism",
    title="Be-Star Grism Campaign",
    claim=(
        "Multi-season, ~3-day-cadence Hα equivalent-width monitoring of "
        "bright emission-line stars with a 0.5 m slitless grism — a survey "
        "blind spot, since every target saturates ZTF/ASAS-SN/ATLAS/Gaia — "
        "validated one-to-one against BeSS and anchored to stability "
        "standards. Deliverables are nightly EW curves for ~10 stars over "
        "2024-12 → 2026-07, outburst epochs to 1–3 d, a short-period search "
        "on the three targets with genuine multi-hour intra-night spans, and "
        "photospheric-corrected disk radii. The venue scales with the "
        "verified-active count and the abstract is not drafted until Step −1 "
        "reports."),
    venue="ApJ if ≥4 BeSS-verified active emitters; AJ/PASP otherwise",
    strategy=_BE,
    phases=(
        Phase("Step −1 — Feasibility gates",
              "Run before any per-target pipeline effort: they set the "
              "sample, the venue, and whether the reader can be trusted.",
              (
                  Task("BE-S-1a-bess", "BeSS emission-state check, all ten targets",
                       "A verified-active count that fixes the venue and the "
                       "science ranking before the abstract exists.",
                       "S0c", _be_src("§4 Step −1(a)"), PENDING),
                  Task("BE-S-1b-lameri-injection",
                       "λ Eri injection–recovery on the real timestamps",
                       "The completeness map that justifies λ Eri's "
                       "demotion to the slow tier.",
                       "S0c", _be_src("§4 Step −1(b)"), PENDING),
                  Task("BE-S-1c-erac-geometry",
                       "Characterise the era-C FITS repackaging",
                       "The HDU layout of the post-2026-04-25 files, and "
                       "repaired geometry in the catalog — the 'NAXIS1=8' "
                       "frames are full-frame images, not 8-pixel strips.",
                       "S0e", _be_src("§3.5"), DONE,
                       evidence="docs/pipeline/s0e_geometry_fix.html"),
                  Task("BE-S-1c-hrg-bandpass",
                       "Settle the hrg O₂ B-band coverage question",
                       "Either a telluric anchor at 6867 Å for hrg, or a "
                       "stated fallback to continuum cross-correlation.",
                       "G", _be_src("§4 Step −1(c)"), PENDING),
              )),
        Phase("Steps 0–2 — Master table, QC and calibration",
              "One master frame table keyed on measured instrument era, then "
              "the gates every frame must pass.",
              (
                  Task("BE-S0-master-table", "Master frame table",
                       "stage_bestar_grism: the grism-whitelisted science "
                       "set with era, night and calibration family attached.",
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
                  Task("BE-S0-header-rescrape",
                       "Re-scrape FITS headers for temperature",
                       "A per-frame detector temperature series — the "
                       "catalog's camtemp is NULL for every grism frame.",
                       "S0", _be_src("§4 Step 0"), PENDING),
                  Task("BE-S0-times", "BJD_TDB for every grism frame",
                       "Barycentric timestamps before any period analysis.",
                       "S3", _be_src("§3.3"), DONE,
                       evidence="docs/pipeline/s3_timing.html"),
                  Task("BE-S1-qc-gates", "Frame QC gates",
                       "A logged pass/reject per frame, whose trace-peak "
                       "distribution IS the paper's saturation statement.",
                       "G", _be_src("§4 Step 1"), PENDING),
                  Task("BE-S2-calibration", "Era-C calibration set",
                       "Real biases, darks and flats for the QHY/Fast era, "
                       "plus the answer on whether the era-B camera survives.",
                       "OPS", _be_src("§4 Step 2"), BLOCKED,
                       evidence="ops/2026-08_observatory_request.md",
                       blocker="The only external blocker of this project. "
                               "Winer is closed until October 2026; the "
                               "request is submitted (Item B) and first in "
                               "queue. The interim surrogate-dark scheme "
                               "carries era B meanwhile."),
              )),
        Phase("Steps 3–6 — Extraction and the calibration chain",
              "Trace, wavelength, delivered resolution, and the response "
              "chain rebuilt on verified overlaps.",
              (
                  Task("BE-S3-extraction", "Per-frame trace and optimal extraction",
                       "Horne-extracted spectra with a per-frame trace fit — "
                       "never reused between frames.",
                       "G", _be_src("§4 Step 3"), PENDING),
                  Task("BE-S4-wavelength", "Wavelength calibration per (filter, era)",
                       "A dispersion solution per filter/era plus a per-frame "
                       "zero point, with the achieved telluric-anchor rms "
                       "in Å reported.",
                       "G", _be_src("§4 Step 4"), PENDING),
                  Task("BE-S5-resolution", "Measure the delivered resolution",
                       "The number that decides scope: V/R decomposition if "
                       "FWHM at Hα < ~4 Å, otherwise the asymmetry moment.",
                       "G", _be_src("§4 Step 5"), PENDING),
                  Task("BE-S6-response-chain", "Response and telluric correction",
                       "A sensitivity curve per (filter, era, season): era C "
                       "on Vega, era B transferred via η Hya/θ Vir, era A on "
                       "the 2024-05-20 Vega ladder, season 1 relative only.",
                       "G", _be_src("§4 Step 6"), PENDING),
              )),
        Phase("Steps 7–10 — Measurement and error calibration",
              "Equivalent widths on fixed windows, three flux tiers, and two "
              "honestly separated error floors.",
              (
                  Task("BE-S7-ew", "Equivalent widths on the paper-wide windows",
                       "Per-frame and nightly EW with the photospheric "
                       "template subtracted before any disk-radius "
                       "conversion.",
                       "G", _be_src("§4 Step 7"), PENDING),
                  Task("BE-S8-flux-tiers", "Three-tier continuum/flux calibration",
                       "A pseudo-r flux series that breaks the EW/continuum "
                       "degeneracy — where Tier 2 exists.",
                       "G", _be_src("§4 Step 8"), PENDING),
                  Task("BE-S9-era-crosscal", "Era cross-calibration",
                       "One offset per (star, band) from overlaps, verified "
                       "star-independent to <1%.",
                       "G", _be_src("§4 Step 9"), PENDING),
                  Task("BE-S10-error-floors", "Two-floor error calibration",
                       "Per-night errors from standards where they exist "
                       "(2025-12 onward), labelled intra-night lower bounds "
                       "elsewhere; precision quoted per star in Å.",
                       "G", _be_src("§4 Step 10"), PENDING),
              )),
        Phase("Steps 11–13 — Time series, events, external",
              "Two tiers, both with published windows and closed "
              "completeness contours.",
              (
                  Task("BE-S11-timeseries", "Slow- and short-tier searches",
                       "GLS + PDM on nightly medians for every target, and "
                       "the 0.3–2 d search on Phecda/Spica/φ Leo only.",
                       "G", _be_src("§4 Step 11"), PENDING),
                  Task("BE-S12-event-timing", "Outburst/disk-event epochs",
                       "Onset epochs to ~1–3 d from template fits — the "
                       "money figure if one is TESS-coincident.",
                       "G", _be_src("§4 Step 12"), PENDING),
                  Task("BE-S13-external", "BeSS and TESS co-analysis",
                       "The resolution-matched one-to-one BeSS validation "
                       "plot — the figure the paper stands on.",
                       "G", _be_src("§4 Step 13"), PENDING),
                  Task("BE-figures", "The twelve-figure set",
                       "Figures 1–12: campaign overview, instrument eras, "
                       "calibration performance, the flat and noise panel, "
                       "the stability-standards panel, the BeSS validation "
                       "(the credibility figure), the main EW curves, the "
                       "event analysis, periodograms, injection–recovery "
                       "maps, disk physics, and the conditional V/R series.",
                       "G", _be_src("§7 Figure list"), PENDING),
                  Task("BE-draft", "Write the manuscript",
                       "manuscripts/BeStar_Grism/main.tex through the "
                       "manuscript outline.",
                       "G", _be_src("§8 Manuscript outline"), PENDING,
                       depends_on=("BE-figures",)),
              )),
    ),
)


_SN = "SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md"


def _sn_src(section: str) -> Source:
    return Source(_SN, section)


SN2023IXF = Project(
    key="SN2023ixf_LightCurve",
    title="SN 2023ixf Early Light Curve",
    claim=(
        "A homogeneous, single-instrument nightly gri and narrowband record "
        "of SN 2023ixf from +5.4 to +50 d, released with the saturation "
        "matrix that defines where it is trustworthy, validated against the "
        "published world dataset. The unconditional clean start is +5.4 d; "
        "'+1.6 d' appears only in the saturated-frame inventory. The upside — "
        "and the reason ApJ is on the table at all — is the slot-'6' "
        "slitless grism series from +3.5 d, which may be the archive's only "
        "unsaturated flash-phase Hα record, and an Hα-band evolution curve "
        "forward-modelled through a recovered transmission profile. How much "
        "of that series is in fact a grism series is MEASURED, not assumed: "
        "the panel immediately below is the measurement, and this upside is "
        "worth exactly what that panel says it is worth."),
    venue="AJ/PASP base case; ApJ upside decided at week 3, post-Gate 0",
    strategy=_SN,
    #: The claim rests entirely on slot '6' being dispersed, so the renderer
    #: prints S2c's live verdict for it directly beneath the claim.  That is
    #: what stopped a paragraph asserting "the 83-frame grism series" from
    #: sitting a hundred lines above a table measuring 3 of those 83 frames
    #: as direct imaging, 19 as undecided, and every dispersed verdict at
    #: 'low' strength — with nothing connecting the two statements.
    claim_filters=("6",),
    phases=(
        Phase("Gate 0 — everything is downstream of this",
              "Three parallel blocking activities: the frozen manifest, the "
              "saturation census, and the grism triage.",
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
                       "S2", _sn_src("§4 Step 0b"), BLOCKED,
                       blocker="The screen needs a MEASURED clip, and the "
                               "S2 tables that held it were destroyed by "
                               "the S0 table swap — the '~3.5 kADU' figure "
                               "currently has no query behind it. Clears "
                               "when S2 re-runs."),
                  Task("SN-G0c-grism-triage",
                       "Contamination-hardened grism triage (2-week timebox)",
                       "A promote-or-appendix verdict on the 83 filter-'6' "
                       "frames, gated on an offset-trace contamination test.",
                       "G", _sn_src("§4 Step 0c"), PENDING),
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
                       IN_PROGRESS,
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
                       "manufacturer sheet or a monochromator scan.",
                       "OPS", _sn_src("§4 Step 1(d) / Step 6"), BLOCKED,
                       blocker="A colour-term regression cannot measure a "
                               "narrowband profile. Without a sheet or a "
                               "scan, Q2 is pre-agreed to demote to a "
                               "methods demonstration. Clears when MACRO "
                               "records or a member campus supply the curve."),
                  Task("SN-S2-linearity", "Linearity audit",
                       "A published linearity curve from the 0.5 s/2 s "
                       "ladder we already own, and the star-side screen.",
                       "S2", _sn_src("§4 Step 2"), BLOCKED,
                       blocker="Same destroyed S2 tables; the clip value "
                               "this fixes is exactly what is unbacked. "
                               "Clears when S2 re-runs."),
              )),
        Phase("Steps 4–6 — Calibration, photometry, the Hα curve",
              "Two photometric regimes with a published overlap test; the "
              "narrowband product forward-modelled, never raw.",
              (
                  Task("SN-S4-resolve", "Re-solve the unsolved broadband frames",
                       "WCS for the ~177 unsolved broadband frames.",
                       "S1b", _sn_src("§4 Step 4"), IN_PROGRESS),
                  Task("SN-S4-ensemble-cal", "REFCAT2 ensemble calibration",
                       "Per-frame zero points, one campaign colour term per "
                       "filter, nightly extinction plus a second-order "
                       "(colour × airmass) term.",
                       "S4", _sn_src("§4 Step 4"), PENDING),
                  Task("SN-S5-photometry", "Two-regime photometry",
                       "Aperture photometry bright, template-subtracted "
                       "forced PSF photometry faint, with the overlap "
                       "agreement published as a systematic.",
                       "S4", _sn_src("§4 Step 5"), PENDING),
                  Task("SN-S6-halpha-curve", "Hα-band evolution curve",
                       "Nightly H-band flux over the census-approved epochs, "
                       "forward-modelled through the bandpass — never "
                       "presented as a raw 'Hα light curve'.",
                       "S4", _sn_src("§4 Step 6"), BLOCKED,
                       blocker="Two upstream gates: the saturation census "
                               "(which epochs are clean) and the "
                               "transmission curve (whether the product is "
                               "physical at all). Clears when SN-G0b and "
                               "SN-S1-narrowband-curves clear."),
              )),
        Phase("Steps 7–10 — Limits, late time, release",
              "Consistency rather than constraint; limits rather than "
              "detections; the manifest as the single source of truth.",
              (
                  Task("SN-S7-model-consistency", "Model consistency fits",
                       "Morag+23 / Martinez+24 tracks at published "
                       "parameters, labelled consistency — no independent t0.",
                       "S4", _sn_src("§4 Step 7"), PENDING),
                  Task("SN-S8-variability-limits", "Variability and bump limits",
                       "A computed spectral window and a 90%-recovery "
                       "amplitude-vs-period curve on an honest grid.",
                       "S4", _sn_src("§4 Step 8"), PENDING),
                  Task("SN-S9-late-time", "Late-time stacks and limits",
                       "Per-stack 5σ upper limits — a late-time detection is "
                       "forbidden.",
                       "S4", _sn_src("§4 Step 9"), PENDING),
                  Task("SN-venue-decision", "The venue decision",
                       "AJ/PASP or ApJ, decided at week 3 immediately after "
                       "Gate 0 lands — not at referee stage.",
                       "S4", _sn_src("§2 venue call"), PENDING),
                  Task("SN-S10-release", "Machine-readable release",
                       "Per-frame and nightly tables with census flags, the "
                       "saturation matrix, pipeline on GitHub, Zenodo "
                       "deposit, README regenerated from the manifest.",
                       "S4", _sn_src("§4 Step 10"), PENDING),
                  Task("SN-figures", "The twelve-figure set",
                       "Figures 1–12: campaign overview, filter "
                       "identification, the saturation census and linearity, "
                       "calibration honesty, the gri light curve, the "
                       "RLMT-vs-the-world centrepiece, the conditional "
                       "Hα-band and grism panels, nightly sampling, the "
                       "window function and injection–recovery, model "
                       "consistency, and the template-subtraction appendix.",
                       "S4", _sn_src("§7 Figure list"), PENDING),
                  Task("SN-draft", "Write the manuscript",
                       "manuscripts/SN2023ixf_LightCurve/main.tex through "
                       "the outline — including replacing the skeleton's "
                       "sec:obs pointer at the T CrB strategy with a "
                       "shared-facility paragraph both papers cite.",
                       "S4", _sn_src("§8 Manuscript outline"), PENDING,
                       depends_on=("SN-figures",)),
              )),
    ),
)


_DW = "DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md"


def _dw_src(section: str) -> Source:
    return Source(_DW, section)


DWARF_AGN = Project(
    key="DwarfGalaxy_AGN_Survey",
    title="Dwarf Galaxies + NGC 5548",
    claim=(
        "Hα detections and flux upper limits for the 13 of 19 published "
        "DESI-Legacy dwarf-candidate fields that have Hα imaging, with "
        "deep-imaging VETTING restricted to candidates the per-field "
        "injection–recovery limits can actually reach; NGC 5238 — the "
        "archive's strongest time-series dataset — as the second pillar; and "
        "a band-integrated NGC 5548 nightly light curve intercalibrated with "
        "ZTF/ASAS-SN/ATLAS, conditional on a coverage/precision table showing "
        "a real gain over the public surveys. Explicitly not a lag, not a "
        "period, not a distance, and not a confirmation."),
    venue="ApJ (NGC 5548 section conditional; pre-authorised RNAAS fallback)",
    strategy=_DW,
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
                       "Measure whether slot '6' disperses",
                       "A per-frame direct/dispersed verdict for the "
                       "NGC 5548 campaign — the measurement that decides "
                       "whether the AGN section is photometry at all.",
                       "S2c", _dw_src("§4 Phase 0.2(b)"), IN_PROGRESS,
                       evidence="docs/pipeline/s2c_filter_identity.html"),
                  Task("DW-P02-filter-dossier",
                       "The filter-slot dossier, science and calibration",
                       "An appendix table mapping every slot letter to a "
                       "physical bandpass per epoch, from config logs AND a "
                       "colour-locus regression — both required, and they "
                       "must agree.",
                       "S2c", _dw_src("§4 Phase 0.2"), BLOCKED,
                       blocker="Half the dossier is John Cannon's: the "
                               "spring-2023 wheel configs, the Hα "
                               "transmission curve, the W flats and "
                               "Dw1643+07's provenance. He has not been "
                               "contacted. Clears on that conversation — "
                               "which the monsoon closure does not block."),
                  Task("DW-P03-plate-solve", "Plate-solve everything",
                       "WCS, FWHM and a pointing-jitter distribution for "
                       "NGC 5548 for the first time; ~250 unsolved Dw rows "
                       "recovered.",
                       "S1b", _dw_src("§4 Phase 0.3"), BLOCKED,
                       blocker="NGC 5548 has zero WCS solutions, ever — 153 "
                               "frames, none solved — and the aperture, sky "
                               "and error-budget decisions of Phase 4 are "
                               "all deferred until this closes. Clears when "
                               "the S1b batch reaches these frames."),
                  Task("DW-P04-noise-model", "StackPro/High Gain noise model",
                       "Photon-transfer curves from the named 2023-06-07/08 "
                       "dark repeats, an effective gain and read noise, and "
                       "the ADU level where variance departs from linear.",
                       "S2", _dw_src("§4 Phase 0.4"), BLOCKED,
                       blocker="S2's PTC and detector tables were destroyed "
                               "by the S0 table swap; χ², F_var and every "
                               "depth number depend on this. Clears when S2 "
                               "re-runs."),
              )),
        Phase("Phase 1 — Calibration frames",
              "An executable path per band, including for the band that has "
              "no flat at any epoch.",
              (
                  Task("DW-P11-flats", "Flat strategy per campaign per filter",
                       "Delta sky flats where dithering allows; for slot '6' "
                       "either a dither flat or a COMPUTED jitter × gradient "
                       "bound carried in the error budget.",
                       "S0b", _dw_src("§4 Phase 1.1"), BLOCKED,
                       evidence="docs/pipeline/s0b_calibration_inventory.html",
                       blocker="No flat for slot '6' exists at any epoch, "
                               "and 'W' has a single 0.4 s frame. The "
                               "acquisition rides the October request and is "
                               "contingent on the retired cameras still "
                               "existing; if they are gone, the Phase 1.1 "
                               "decision trees stand as written."),
                  Task("DW-P12-fringe-moon", "Fringing and moonlight gradients",
                       "A measured answer per filter rather than an "
                       "assumption, plus per-moon-regime delta flats if the "
                       "Hα residuals demand them.",
                       "S4", _dw_src("§4 Phase 1.2"), PENDING),
              )),
        Phase("Phase 2 — Photometric calibration",
              "Shared across all three campaigns; ZMAG is QC only.",
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
                       blocker="A transformation needs a band, and the band "
                               "comes from the filter dossier. Clears with "
                               "DW-P02-filter-dossier."),
              )),
        Phase("Phase 3 — Dwarf fields (the lead science)",
              "Cross-match first — it fixes the counts, the framing and the "
              "title.",
              (
                  Task("DW-P31-crossmatch", "Literature cross-match per candidate",
                       "The table that defines the remaining novelty: "
                       "ELVES/HSC status, existing Kaisin & Karachentsev Hα, "
                       "and DESI spectroscopy per candidate.",
                       "S0c", _dw_src("§4 Phase 3.1"), PENDING),
                  Task("DW-P32-qc-gates", "Uniform QC gates",
                       "A pipeline-emitted disposition table — the only "
                       "rejection authority; no hand-written night lists "
                       "anywhere.",
                       "S4", _dw_src("§4 Phase 3.2"), PENDING),
                  Task("DW-P33-stacks", "Weighted coadds per field per filter",
                       "Deep stacks with a constant per-frame sky — never a "
                       "mesh background, which eats LSB flux.",
                       "S4", _dw_src("§4 Phase 3.3"), PENDING),
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
                  Task("DW-P36-halpha-limits", "Hα fluxes and upper limits",
                       "THE lead result: continuum-subtracted Hα flux "
                       "detections and limits for the 13 fields with Hα, "
                       "with SFR only at an explicit fiducial distance.",
                       "S4", _dw_src("§4 Phase 3.6"), BLOCKED,
                       blocker="No flux can be quoted without the actual Hα "
                               "filter transmission curve (centre, width) — "
                               "Cannon holds it. Clears with "
                               "DW-P02-filter-dossier."),
                  Task("DW-P37-photometry-mode",
                       "Aperture photometry with growth-curve corrections",
                       "Aperture + curve-of-growth photometry (photutils/"
                       "sep) as the declared mode — no PSF fitting, which at "
                       "5.4″ FWHM in uncrowded fields buys nothing.",
                       "S4", _dw_src("§4 Phase 3.7"), PENDING),
              )),
        Phase("Phase 4 — NGC 5548",
              "Statistics only, no models with free timescales, and the "
              "whole section conditional on demonstrating it adds anything.",
              (
                  Task("DW-P41-aperture", "Host-aware aperture and sky",
                       "An aperture chosen from the measured FWHM history, "
                       "sky off the host, and a host fraction quoted "
                       "alongside F_var.",
                       "S1b", _dw_src("§4 Phase 4.1"), BLOCKED,
                       blocker="Explicitly deferred until Phase 0.3: with "
                               "zero solved frames the campaign's seeing "
                               "history is unknown, so no aperture may be "
                               "committed. Clears with DW-P03-plate-solve."),
                  Task("DW-P42-ensemble", "Differential ensemble light curve",
                       "Nightly means with comparison-star null light "
                       "curves and χ²-validated empirical errors.",
                       "S1b", _dw_src("§4 Phase 4.2–4.3"), BLOCKED,
                       blocker="Needs WCS and the aperture decision above. "
                               "Clears with DW-P03-plate-solve."),
                  Task("DW-P44-statistics", "Variability statistics",
                       "χ² vs constancy, F_var with host dilution stated, "
                       "and a structure function with per-bin pair counts "
                       "published.",
                       "S4", _dw_src("§4 Phase 4.4"), PENDING),
                  Task("DW-P45-gate", "The Phase 4.5 value gate",
                       "A coverage/precision table against ZTF g / ASAS-SN g "
                       "/ ATLAS o. If it shows no gain, the RNAAS fallback "
                       "executes without regret.",
                       "S4", _dw_src("§4 Phase 4.5"), PENDING),
                  Task("DW-P46-bjd",
                       "BJD_TDB for everything in the paper",
                       "Barycentric mid-exposure times for every NGC 5548 "
                       "frame — the headers carry UTC JD, which is not what "
                       "a variability paper may quote.",
                       "S3", _dw_src("§4 Phase 4.6"), PENDING),
              )),
        Phase("Phase 5 — NGC 5238 and the variability census",
              "The same machinery, applied to the archive's densest "
              "time-series field.",
              (
                  Task("DW-P51-ngc5238", "NGC 5238 stacks and Hα map",
                       "Per-band deep stacks and a continuum-subtracted Hα "
                       "star-formation map of an actively interacting dwarf.",
                       "S4", _dw_src("§4 Phase 5.1"), PENDING),
                  Task("DW-P52-subtraction", "Image-subtraction transient search",
                       "A difference-imaging pass across the 21 epochs.",
                       "S4", _dw_src("§4 Phase 5.2"), PENDING),
                  Task("DW-P53-period-search", "Field-star period search",
                       "LS + PDM with bootstrap FAPs on real timestamps and "
                       "explicit alias adjudication.",
                       "S4", _dw_src("§4 Phase 5.3"), PENDING),
                  Task("DW-P54-completeness", "Injection–recovery completeness map",
                       "The map that converts non-detections into a "
                       "publishable statement.",
                       "S4", _dw_src("§4 Phase 5.4"), PENDING),
                  Task("DW-P55-detrending",
                       "Covariate detrending discipline",
                       "Frame-level regression against airmass, FWHM, x/y "
                       "drift and sky with at most two Sys-Rem-like "
                       "components, and multiband Lomb–Scargle with floating "
                       "per-band offsets when the L and R+Hα blocks merge. "
                       "NEVER GP-detrend in time before a period search — "
                       "trend and signal are fitted simultaneously.",
                       "S4", _dw_src("§4 Phase 5.5"), PENDING),
                  Task("DW-P56-eclipse-timing",
                       "Eclipse timing, only if eclipsing binaries emerge",
                       "Kwee–van Woerden minima for fully covered events and "
                       "template fits otherwise, with period refinement to "
                       "~10⁻⁴ d — and NO period-CHANGE claims, which this "
                       "baseline cannot support.",
                       "S4", _dw_src("§4 Phase 5.6"), PENDING),
                  Task("DW-figures", "The fourteen-figure set",
                       "Figures 1–14: survey footprint, StackPro "
                       "characterization, filter identification, flat-field "
                       "validation, the depth/detectability gate, the "
                       "candidate atlas, Sérsic parameters, the Hα results, "
                       "the NGC 5238 variability search, completeness maps, "
                       "and the NGC 5548 panels (12 conditional on the "
                       "Phase 4.5 gate).",
                       "S4", _dw_src("§7 Figure list"), PENDING),
                  Task("DW-draft", "Write the manuscript",
                       "manuscripts/DwarfGalaxy_AGN_Survey/main.tex through "
                       "the manuscript outline. Every title and abstract "
                       "count is drawn from the Phase 3.1 and Phase 0 "
                       "tables, never typed independently — the outline says "
                       "so outright, and the title is fixed only after 3.1.",
                       "S4", _dw_src("§8 Manuscript outline"), PENDING,
                       depends_on=("DW-figures",)),
              )),
    ),
)


#: Legacy_Rigel has no committee strategy, so its tasks derive from the
#: standing decisions in THIS module.  They used to cite
#: "docs/Legacy_Rigel/index.html \u00a7Status" — a heading on a page this
#: renderer generates, which the first render deleted, leaving eleven
#: citations pointing at nothing and two tasks naming the page as the
#: evidence for its own claims.  A generated page cannot be a source.
_RIG = "pipeline/macro_core/project_plan.py"


#: The decision heading most Legacy tasks derive from.
_RIG_S0D = "§Multi-archive keying (stage S0d)"


def _rig_src(section: str) -> Source:
    return Source(_RIG, section)


LEGACY_RIGEL = Project(
    key="Legacy_Rigel",
    title="Rigel-Era Legacy Archive (candidate)",
    claim=(
        "UNDECIDED — and that is the honest state. This is the pre-MACRO "
        "Winer/Iowa archive: about 210,120 files across ~1,100 nights, "
        "2015–2023 (a file count from the transfer itself — NOT a "
        "pipeline-emitted frame count, because nothing here has reached the "
        "manifest yet), on a different telescope (Rigel System, FLI ProLine "
        "PL16803, Johnson BVRI). The first frames scanned show a programme "
        "dominated by "
        "contact binaries and eclipsing systems, only a small fraction "
        "overlapping the RLMT-era target lists. NO committee strategy exists "
        "for it, so this page carries an INGEST plan, not a science plan: the "
        "census is what decides whether there is a sixth paper here and what "
        "it would claim. Every task below cites the archive page's recorded "
        "decisions rather than a strategy document, because there is none to "
        "cite."),
    venue="undecided — a go/no-go, not a venue posture",
    strategy="",
    decisions=(
        ("Why a separate archive root",
         "The legacy data lives in legacy-archive/, never inside "
         "rlmt-archive/. Two reasons, both about not corrupting what already "
         "works. That tree is a byte-verified mirror of the Linode "
         "'testimages' bucket (rclone check: 0 differences), and foreign "
         "files would break that equivalence permanently. And this is a "
         "different telescope, detector and filename convention, so merging "
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
         "The catalog gains archive_root and telescope columns, plus "
         "telescope-aware era keying so a Rigel configuration can never "
         "share a camera era with an RLMT one, and a Rigel filename parser. "
         "The full target census — and therefore the real science case — "
         "lands when the transfer completes."),
    ),
    phases=(
        Phase("Phase L0 — Ingest",
              "Make the pipeline multi-archive without letting legacy frames "
              "leak into RLMT-era analyses.",
              (
                  Task("RIG-L0-transfer", "Complete the archive transfer",
                       "A byte-verified legacy-archive/ mirror, separate "
                       "from the rlmt-archive/ tree whose rclone-check "
                       "equivalence must not be broken.",
                       "S0", _rig_src("§Why a separate archive root"),
                       IN_PROGRESS),
                  Task("RIG-L0-multi-archive", "Multi-archive catalog keying",
                       "archive_root and telescope columns, and "
                       "telescope-aware era keying so a Rigel configuration "
                       "can never share a camera era with an RLMT one.",
                       "S0", _rig_src(_RIG_S0D), IN_PROGRESS),
                  Task("RIG-L0-filename-parser", "Rigel filename parser",
                       "Target, filter and exposure recovered from the "
                       "legacy filename convention.",
                       "S0", _rig_src(_RIG_S0D), IN_PROGRESS),
                  Task("RIG-L0-header-scan", "Header scan of the legacy tree",
                       "A catalog row per legacy frame — the input S0 needs "
                       "before any census is possible.",
                       "S0", _rig_src("§What the first scan shows"), PENDING),
              )),
        Phase("Phase L1 — Census",
              "The full target census — and therefore the real science case "
              "— lands only when the transfer completes.",
              (
                  Task("RIG-L1-target-census", "Target census",
                       "Frames and nights per target across ~1,100 nights, "
                       "alias-merged the way S0 merges RLMT targets.",
                       "S0", _rig_src(_RIG_S0D), PENDING),
                  Task("RIG-L1-night-census", "Night and coverage census",
                       "Per-target baselines, cadences and full-cycle night "
                       "counts — what a contact-binary programme can carry.",
                       "S0", _rig_src(_RIG_S0D), PENDING),
                  Task("RIG-L1-overlap", "Overlap with the RLMT-era targets",
                       "The list of targets both archives observed — the "
                       "only place a legacy frame could legitimately extend "
                       "an RLMT baseline.",
                       "S0", _rig_src("§What the first scan shows"), PENDING),
                  Task("RIG-L1-calibration-census", "Legacy calibration census",
                       "Whether mode-matched calibration exists for the "
                       "Rigel-era configurations at all.",
                       "S0b", _rig_src(_RIG_S0D), PENDING),
              )),
        Phase("Phase L2 — Decision",
              "Whether this becomes a sixth project, and on what terms.",
              (
                  Task("RIG-L2-strategy", "Commission a committee strategy",
                       "An ANALYSIS_STRATEGY.md for this archive — or a "
                       "recorded decision not to write one.",
                       "STRAT", _rig_src("§What the first scan shows"), PENDING),
                  Task("RIG-L2-gonogo", "Sixth-project go/no-go",
                       "A decision, with the census as its evidence, on "
                       "whether the legacy archive carries a paper.",
                       "STRAT", _rig_src("§What the first scan shows"), PENDING),
                  Task("RIG-L2-page", "Publish the census on this page",
                       "This page's science section, replacing the ingest "
                       "status it currently carries.",
                       "WEB", _rig_src(_RIG_S0D), PENDING),
              )),
    ),
)


#: The six project plans, in the order the hub lists them.
PROJECTS: tuple[Project, ...] = tuple(_build(p) for p in (
    TCRB_MONITORING, CV_TIMESERIES, SN2023IXF, BESTAR_GRISM, DWARF_AGN,
    LEGACY_RIGEL))

PROJECT_BY_KEY: dict[str, Project] = {p.key: p for p in PROJECTS}


# ===========================================================================
# 4.  PURE HELPERS
# ===========================================================================

def all_tasks() -> tuple[Task, ...]:
    """Every task in the ledger, in project then phase then ledger order."""
    return tuple(t for p in PROJECTS for t in p.tasks)


def tasks_of(project_key: str) -> tuple[Task, ...]:
    """Every task of one project."""
    return PROJECT_BY_KEY[project_key].tasks


def task_by_id(task_id: str) -> Task:
    """Look one task up, or raise — a typo'd id must never write a status
    row for a task that does not exist."""
    for t in all_tasks():
        if t.id == task_id:
            return t
    raise PlanError(f"unknown task id: {task_id!r}")


def project_of(task_id: str) -> Project:
    """The project a task belongs to."""
    return PROJECT_BY_KEY[task_by_id(task_id).project]


def overlay_statuses(tasks: Sequence[Task],
                     recorded: Mapping[str, str]) -> dict[str, str]:
    """Effective status per task id: the recorded value when one exists,
    the ledger's declared value otherwise.

    The ledger is the plan's opening position; the table is what happened.
    A recorded value always wins, INCLUDING when it moves a task backwards —
    ``sync`` writing ``redo_needed`` over a ledger ``done`` is the single
    most important thing this function has to let through.
    """
    out: dict[str, str] = {}
    for t in tasks:
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
    """``(done, total)``.  ``redo_needed`` counts as NOT done — that is the
    whole reason the state exists."""
    total = sum(counts.values())
    return counts.get(DONE, 0), total


def next_up(tasks: Sequence[Task], statuses: Mapping[str, str],
            limit: int = 3,
            include: Sequence[str] = OPEN_STATUSES) -> tuple[Task, ...]:
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
    """
    rank = {REDO_NEEDED: 0, IN_PROGRESS: 1, PENDING: 2}
    actionable, seen_stages = [], set()
    for t in tasks:
        status = statuses[t.id]
        if status not in rank or status not in include:
            continue
        if unmet_dependencies(t, statuses):
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
                stage_states: Mapping[str, str]) -> tuple[Change, ...]:
    """Statuses that the DATABASE, not a person, has already decided.

    Exactly two rules, and no others may be added without a reason written
    here — an automatic rule that can overwrite a human judgement is how a
    tracking system starts lying in the other direction:

    1. ``done`` + stage not FRESH  -> ``redo_needed``.  The evidence this
       task's claim rests on is stale, missing or was never recorded, so the
       claim is not backed any more.
    2. ``redo_needed`` + stage FRESH -> ``done``.  The stage was re-run; the
       claim is backed again, and the plan should not keep nagging.

    ``pending`` / ``in_progress`` / ``blocked`` are never touched: no
    fingerprint can tell whether a person has started something.
    """
    changes: list[Change] = []
    for t in tasks:
        state = stage_states.get(t.stage)
        if state is None:
            continue
        now = statuses[t.id]
        if now == DONE and state != pv.FRESH:
            changes.append(Change(
                t.id, now, REDO_NEEDED,
                f"stage {t.stage} is {state} — the evidence this task's "
                f"'done' rests on is no longer backed"))
        elif now == REDO_NEEDED and state == pv.FRESH:
            changes.append(Change(
                t.id, now, DONE,
                f"stage {t.stage} is FRESH again — the claim is backed"))
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
    for task in all_tasks():
        name = task.source.document
        if name in docs:
            continue
        path = os.path.join(str(repo_root), name)
        if os.path.exists(path):
            with open(path, encoding="utf-8", errors="replace") as fh:
                docs[name] = fh.read()
    return docs


def verify_citations(repo_root) -> tuple[str, ...]:
    """Resolve every citation in the ledger against the real documents."""
    return citation_problems(
        [(t.id, t.source) for t in all_tasks()],
        read_cited_documents(repo_root))


def unmet_dependencies(task: Task,
                       statuses: Mapping[str, str]) -> tuple[str, ...]:
    """Dependency ids of ``task`` that are not ``done``."""
    return tuple(d for d in task.depends_on if statuses.get(d) != DONE)


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
    for project in PROJECTS:
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
            if task.evidence.startswith(f"docs/{project.key}/"):
                raise PlanError(
                    f"{task.id}: names its own project page as its "
                    f"evidence — a page cannot be the evidence for the "
                    f"claims it makes")

    # Dependencies are checked in a second pass: a task may legitimately
    # depend on one declared later in the same phase, so the ids only all
    # exist once the first pass has seen them.
    for task in all_tasks():
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
    edges = {t.id: t.depends_on for t in all_tasks()}
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
