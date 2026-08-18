"""S1b production astrometry batch: queue policy + pure batch logic.

ROADMAP §1.1 S1 → S1b: the go/no-go experiment (``macro_core.astrom``,
report ``docs/pipeline/s1_astrometry.html``) measured per-stratum solve
rates on 472 sampled frames and the Week-2 review accepted the verdicts.
This module turns those verdicts into the PRODUCTION batch policy that
``pipeline/scripts/run_s1_batch.py`` drives over the full 38k-frame
stratified backlog.

Everything here that decides something is a *pure function* (no I/O, no
globals mutated), unit-tested in ``pipeline/tests/test_batch.py``:

* the stratum → (population, priority, QC-gating) policy table,
* queue construction from candidate rows,
* the QC pre-gate that maps autopsy image statistics to skip/attempt,
* sidecar WCS path naming,
* status-transition legality,
* ETA arithmetic with the experiment's measured per-stratum medians.

The batch writes ONE new manifest table, ``s1_batch`` (plus its
``s1_batch_meta`` build-facts sibling) — existing tables are never
modified, and the archive is never written to: accepted WCS solutions
land as SIDECAR files under ``products/astrom/wcs/``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

from . import astrom

# --------------------------------------------------------------------------
# Version + policy constants (single source of truth, quoted by reports)
# --------------------------------------------------------------------------

#: Version string recorded into ``s1_batch_meta``.
S1B_CODE_VERSION = "S1b v1.0 (2026-08-18)"

#: The verdicts the Week-2 review accepted, per stratum, copied from the
#: git-tracked evidence report ``docs/pipeline/s1_astrometry.html``
#: (S1 v1.1, built 2026-08-18T07:28Z, commit c930015).  They are PINNED
#: here as code because the s1_* evidence tables do not survive manifest
#: rebuilds — the report is the durable record, this dict is its
#: machine-readable quotation.  A stratum's verdict decides its gating:
#: GO strata are solved directly; every below-GO stratum runs behind the
#: QC pre-gate ("Populations below GO stay viable as *filtered* batches"
#: — the report's own mandate), which skips starless/defocused/trailed
#: frames before the solver burns its 75 s budget on them.
S1_VERDICTS: dict[str, str] = {
    "cv_mode0_sloan_short": "GO",        # 48/48, 100% [93-100]
    "cv_mode0_sloan_long":  "GO",        # 48/48, 100% [93-100]
    "cv_ikon_sloan":        "GO",        # 48/48, 100% [93-100]
    "cv_gsense_misc":       "CAUTION",   # census 25/40 = 62.5%
    "sn_gsense_broadband":  "NO-GO",     # 22/48, 46% [33-60]
    "dwarf_gsense_deep":    "CAUTION",   # 42/48, 88% [75-94]
    "mode0_backlog_short":  "GO",        # 48/48, 100% [93-100]
    "mode0_backlog_long":   "GO",        # 46/48, 96% [86-99]
    "fast_fullframe":       "CAUTION",   # 41/48, 85% [73-93]
    "ikon_backlog":         "GO",        # 45/48, 94% [83-98]
}

#: Median per-frame wall seconds measured by the experiment, per stratum
#: (section 3 of the report).  Used as the ETA PRIOR until the batch has
#: enough of its own finished frames in a stratum to measure a live
#: median.  Failure-heavy strata carry failure-heavy medians — the ETA is
#: honest about the mix, not just the happy path.
PRIOR_MEDIAN_S: dict[str, float] = {
    "cv_mode0_sloan_short": 4.3,
    "cv_mode0_sloan_long":  4.6,
    "cv_ikon_sloan":        1.7,
    "cv_gsense_misc":       5.0,
    "sn_gsense_broadband":  5.0,
    "dwarf_gsense_deep":    6.6,
    "mode0_backlog_short":  3.9,
    "mode0_backlog_long":   4.0,
    "fast_fullframe":       4.1,
    "ikon_backlog":         2.1,
}

#: A live per-stratum median replaces the prior only after this many
#: finished frames — a median of three frames is noise, not evidence.
MIN_FRAMES_FOR_LIVE_MEDIAN = 20

#: Batch population labels (the CLI filter vocabulary), in priority
#: order: the paper-critical WCS lands first if the run is interrupted.
POP_CV = "cv_polars"
POP_DWARF = "dwarf"
POP_SN = "sn2023ixf"
POP_FACILITY = "facility"
POPULATIONS: tuple[str, ...] = (POP_CV, POP_DWARF, POP_SN, POP_FACILITY)


@dataclass(frozen=True)
class StratumPolicy:
    """The batch policy for one experiment stratum: which population it
    belongs to, where it sits in the queue, and whether the QC pre-gate
    stands between its frames and the solver."""
    stratum_id: str
    population: str
    priority: int        # global queue rank, 1 = solved first
    qc_gated: bool       # True = image-statistics pre-gate applies


#: The queue policy, in priority order.  Population order implements the
#: task's mandate — CV polars first, then dwarf, then the QC-gated SN
#: frames, then the facility backlog — and within a population the
#: paper-critical / highest-rate strata lead.  ``qc_gated`` follows the
#: verdict mechanically: every below-GO stratum is gated (see
#: ``S1_VERDICTS``), which covers the task's "SN frames ONLY behind the
#: QC pre-gate" and extends the same protection to the CAUTION strata
#: whose failures the autopsy attributed to starless/defocused frames.
STRATUM_POLICY: tuple[StratumPolicy, ...] = (
    # -- CV polars: the paper-gating population --------------------------
    StratumPolicy("cv_mode0_sloan_short", POP_CV, 1, False),
    StratumPolicy("cv_mode0_sloan_long",  POP_CV, 2, False),
    StratumPolicy("cv_ikon_sloan",        POP_CV, 3, False),
    StratumPolicy("cv_gsense_misc",       POP_CV, 4, True),
    # -- Dwarf/AGN deep fields ------------------------------------------
    StratumPolicy("dwarf_gsense_deep",    POP_DWARF, 5, True),
    # -- SN 2023ixf: NO-GO as a raw batch, viable as a filtered batch ---
    StratumPolicy("sn_gsense_broadband",  POP_SN, 6, True),
    # -- Facility backlog: the bulk ------------------------------------
    StratumPolicy("mode0_backlog_short",  POP_FACILITY, 7, False),
    StratumPolicy("mode0_backlog_long",   POP_FACILITY, 8, False),
    StratumPolicy("ikon_backlog",         POP_FACILITY, 9, False),
    StratumPolicy("fast_fullframe",       POP_FACILITY, 10, True),
)

#: Fast lookup: stratum_id → its policy row.
POLICY_BY_STRATUM: dict[str, StratumPolicy] = {
    p.stratum_id: p for p in STRATUM_POLICY}

#: Terminal statuses a queue row can reach.  ``failed`` carries a
#: ``fail_kind`` (unsolved | timeout | error) preserving the experiment's
#: finer vocabulary inside the batch's five-status contract.
STATUS_PENDING = "pending"
TERMINAL_STATUSES = frozenset({"solved", "bad_solve", "failed",
                               "skipped_qc"})
ALL_STATUSES = frozenset({STATUS_PENDING}) | TERMINAL_STATUSES
FAIL_KINDS = frozenset({"unsolved", "timeout", "error"})


def allowed_transition(old: str, new: str) -> bool:
    """Legality of a status transition in ``s1_batch``.

    Only ``pending → terminal`` is legal during a run; terminal states
    are immutable (a re-run must never silently overwrite evidence —
    re-attempting a frame requires an explicit requeue that first sets
    it back to pending, which is also legal here so the requeue path is
    covered by the same rule).
    """
    if old == STATUS_PENDING:
        return new in TERMINAL_STATUSES
    # Explicit requeue: any terminal state may return to pending.
    return old in TERMINAL_STATUSES and new == STATUS_PENDING


# --------------------------------------------------------------------------
# Queue construction (pure: candidate row dicts in, ordered queue out)
# --------------------------------------------------------------------------

def build_queue_rows(candidates: Sequence[dict]) -> list[dict]:
    """Classify candidate frames and return the batch queue, in order.

    Input rows are ``astrom.BASE_COLS`` dicts (the same candidate
    universe the experiment designed from).  Frames that fail the
    solvable-candidate gate or fall in no stratum (the 5,183-frame
    heterogeneous residue the report explicitly excludes from the batch)
    are dropped — nothing outside a measured stratum is queued.

    Output: one dict per queued frame with the policy fields attached,
    sorted by (priority, obs_rowid) — a deterministic queue order that
    puts the paper-critical strata first.
    """
    out: list[dict] = []
    for r in candidates:
        sid = astrom.classify_stratum(r)   # None for residue/unsolvable
        if sid is None:
            continue
        pol = POLICY_BY_STRATUM[sid]
        q = dict(r)                        # never mutate the caller's row
        q["stratum_id"] = sid
        q["population"] = pol.population
        q["priority"] = pol.priority
        q["qc_gated"] = int(pol.qc_gated)
        out.append(q)
    # Deterministic order: priority rank first, then rowid within it.
    out.sort(key=lambda q: (q["priority"], q["obs_rowid"]))
    return out


# --------------------------------------------------------------------------
# QC pre-gate: autopsy statistics → attempt or skip, BEFORE solving
# --------------------------------------------------------------------------

#: Saturation rail per readout family for the QC statistics, in ADU —
#: the same values the failure autopsy calibrated its thresholds
#: against (S2 measured the GSENSE High Gain clip at 3,496 ADU; the
#: autopsy's 3,500 rail is that fact rounded, kept identical here so
#: the batch gate and the autopsy evidence share one definition).
GSENSE_HIGHGAIN_SAT_ADU = 3500.0
DEFAULT_SAT_ADU = 65000.0


def saturation_adu_for(readoutm) -> float:
    """ADU level treated as 'saturated' in the QC statistics."""
    r = (readoutm or "").strip().lower()
    return GSENSE_HIGHGAIN_SAT_ADU if r.startswith("high gain") \
        else DEFAULT_SAT_ADU


def qc_pregate(metrics: dict) -> tuple[bool, str]:
    """The batch QC pre-gate: (attempt_solve?, recorded diagnosis).

    ``metrics`` is the dict ``astrom.image_metrics`` returns.  The gate
    IS the experiment's failure-autopsy diagnosis run in reverse: the
    autopsy showed that 55 of 59 sampled failures were starved /
    defocused / trailed / saturated frames — detectable from cheap
    source statistics — so the batch spends ~2 s classifying instead of
    75 s failing.  Thresholds live in ``macro_core.astrom`` (the
    AUTOPSY_* constants, calibrated on the autopsy table and quoted by
    the S1 report) — ONE definition for experiment and batch.

    A frame passes only when the diagnosis finds a healthy star field
    (the autopsy's 'unexplained' branch: stars present, nothing wrong) —
    every pathological diagnosis skips.  Unreadable frames skip too:
    a frame whose pixels cannot be read cannot be solved.
    """
    diagnosis = astrom.diagnose_failure(
        metrics.get("n_sources"), metrics.get("n_psf_sources"),
        metrics.get("median_elongation"), metrics.get("saturated_fraction"),
        metrics.get("bright_median_a_px"))
    if diagnosis.startswith("unexplained"):
        # The autopsy phrase describes a solver failure; at pre-gate
        # time the same statistics mean "healthy field, worth solving".
        return True, "stars present (QC pass)"
    return False, diagnosis


# --------------------------------------------------------------------------
# Sidecar WCS naming (pure path arithmetic — no filesystem access)
# --------------------------------------------------------------------------

def scratch_wcs_name(archive_rel_path: str) -> str:
    """Basename of the ``.wcs`` file solve-field leaves in the scratch
    dir for one frame, mirroring ``astrom.solve_one_frame``'s naming:
    a ``.fz`` member is funpacked to its stem first, then solve-field
    replaces the (remaining) extension with ``.wcs``.

    'a/b/x.fts.fz' → 'x.wcs';  'a/b/x.fts' → 'x.wcs';  'a/b/x' → 'x.wcs'
    """
    name = archive_rel_path.rsplit("/", 1)[-1]
    if name.endswith(".fz"):
        name = name[:-len(".fz")]          # funpack strips the .fz layer
    stem, dot, _ext = name.rpartition(".")
    return (stem if dot else name) + ".wcs"


def sidecar_rel_path(archive_rel_path: str) -> str:
    """Sidecar path (relative to ``products/astrom/wcs/``) for one
    frame's accepted WCS solution.

    The sidecar tree MIRRORS the archive tree — same directories, the
    frame's own basename with ``.wcs`` appended after stripping the
    compression layer — so a solution is findable from its frame path
    by pure string arithmetic, no database in hand:

    'rawimage/2024/x.fts.fz' → 'rawimage/2024/x.fts.wcs'
    """
    rel = archive_rel_path
    if rel.endswith(".fz"):
        rel = rel[:-len(".fz")]            # the sidecar names the FITS,
    return rel + ".wcs"                    # not its compression wrapper


# --------------------------------------------------------------------------
# ETA arithmetic
# --------------------------------------------------------------------------

def stratum_median_s(stratum_id: str, live_median: Optional[float],
                     n_live: int) -> float:
    """Per-frame seconds estimate for one stratum: the batch's own
    measured median once it rests on enough frames, the experiment's
    prior until then, a conservative default for an unknown stratum."""
    if live_median is not None and n_live >= MIN_FRAMES_FOR_LIVE_MEDIAN:
        return live_median
    return PRIOR_MEDIAN_S.get(stratum_id, 6.0)


def eta_seconds(pending_by_stratum: dict[str, int],
                median_by_stratum: dict[str, float],
                workers: int) -> float:
    """Wall-clock ETA: Σ (pending × per-frame median) ÷ workers.

    The same projection formula the S1 report used for its cost table
    (``astrom.projected_hours`` is the hours flavor); medians are per
    stratum so failure-heavy strata charge their own slower rate.
    """
    if workers <= 0:
        return 0.0
    total = 0.0
    for sid, n in pending_by_stratum.items():
        total += n * median_by_stratum.get(
            sid, PRIOR_MEDIAN_S.get(sid, 6.0))
    return total / workers
