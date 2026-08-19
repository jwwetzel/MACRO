"""Stage provenance — the machinery that answers "is what we already did
still true?" without a human having to remember why.

WHY THIS MODULE EXISTS
----------------------
Every stage in this pipeline consumes another stage's output.  When an
upstream stage is re-run — or when a re-characterization changes what an
upstream column MEANS — everything downstream of it silently becomes a
claim about a database that no longer exists.  Nothing in the repo
noticed.  The only defence was a person remembering the dependency graph,
and people publish stale results for exactly that reason.

So this module makes the dependency graph a DATA STRUCTURE and the
"has my input changed?" question an ARITHMETIC one:

* :data:`STAGES` declares, per stage, what it READS, what it WRITES, and
  the code-version constant its build script stamps into the database.
* :func:`fingerprint_resource` reduces one stage output to a cheap,
  deterministic ``(row count, digest)`` pair over the columns that
  downstream stages actually consume.
* ``stage_provenance`` records, at each run, the input fingerprints the
  stage SAW and the output fingerprint it PRODUCED.
* :func:`is_stale` compares the recorded inputs against today's inputs and
  returns FRESH / STALE(reason) / NEVER_RUN / OUTPUT_MISSING.

THE FINGERPRINT DESIGN RULE (the part that is easy to get wrong)
----------------------------------------------------------------
A fingerprint over *every* column is worse than no fingerprint at all.
Several tables here carry columns that change on every rebuild by
construction — ``stage_build_id`` (a timestamp), ``finished_utc``,
``solve_time_s``, ``log_tail``, ``size`` (a filesystem scan artifact).
Hashing those would report every stage permanently STALE, which trains
the reader to ignore the alarm.  Hashing nothing reports everything
permanently FRESH, which is how this project got here.

The rule adopted, and enforced by :data:`RESOURCES` carrying a written
``why`` for every single resource: **hash exactly the columns a
downstream stage reads to make a decision.**  Each spec below names them
and says which consumer reads them.  A resource with no spec raises
:class:`ProvenanceError` rather than guessing — a wrong answer must never
be cheaper to produce than an honest failure.

WHAT IS PURE HERE
-----------------
Everything that decides anything: :func:`digest_rows`, :func:`digest_bytes`,
:func:`compare_fingerprints`, :func:`is_stale`, :func:`topological_order`,
:func:`propagate_staleness` and :func:`rerun_plan` see nothing but plain
values, so ``pipeline/tests/test_provenance.py`` drives them with
hand-built fixtures.  Only :func:`fingerprint_resource`,
:func:`read_records` and :func:`record_run` touch a database or a disk.
"""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Optional, Sequence

__all__ = [
    "PROVENANCE_CODE_VERSION",
    "UNRECORDED",
    "ProvenanceError",
    "ResourceSpec",
    "RESOURCES",
    "Stage",
    "STAGES",
    "STAGE_BY_KEY",
    "producer_of",
    "Fingerprint",
    "digest_rows",
    "digest_bytes",
    "fingerprint_resource",
    "fingerprint_all",
    "read_version_constant",
    "FRESH",
    "STALE",
    "STALE_UPSTREAM",
    "NEVER_RUN",
    "OUTPUT_MISSING",
    "ALL_STATES",
    "Freshness",
    "Record",
    "compare_fingerprints",
    "is_stale",
    "topological_order",
    "propagate_staleness",
    "rerun_plan",
    "ensure_table",
    "record_run",
    "read_records",
    "recorded_run_times",
    "PROVENANCE_TABLE",
]

# Version stamp written into every stage_provenance row.  Bump it whenever
# a change here would alter a stored digest (a changed column list, a
# changed canonicalization rule) — otherwise a later reader would compare
# a v1 digest against a v2 digest and read a rule change as a data change.
PROVENANCE_CODE_VERSION = "P v1.0 (2026-08-18)"

#: Name of the table this module owns inside the manifest database.
PROVENANCE_TABLE = "stage_provenance"

#: Sentinel prefix for a BACKFILLED input whose state at run time is not
#: recoverable.  It exists because the alternative — writing today's
#: fingerprint and calling the stage fresh — would launder an unknown into
#: a reassurance.  Anything carrying this prefix compares as CHANGED.
UNRECORDED = "UNRECORDED"


class ProvenanceError(RuntimeError):
    """Raised when the graph or a fingerprint cannot be computed honestly
    (an undeclared resource, a cyclic DAG, a malformed record)."""


# ===========================================================================
# 1.  RESOURCES — what a stage can read or write, and how it is fingerprinted
# ===========================================================================

@dataclass(frozen=True)
class ResourceSpec:
    """One fingerprintable thing: a manifest table, a table in a sibling
    database, or a file on disk.

    ``key``       stable identifier used in the DAG and in stored records.
    ``kind``      ``'table'`` (manifest DB) | ``'db'`` (sibling sqlite) |
                  ``'file'`` (bytes on disk) | ``'stat'`` (size + mtime of a
                  file too large, or too actively written, to hash).
    ``name``      table name (tables) or repo-relative path (db/file).
    ``database``  for ``kind='db'``: repo-relative path of the sqlite file.
    ``order_by``  SQL ORDER BY that makes the row sequence deterministic —
                  without it two identical databases can hash differently
                  purely because of page layout.
    ``columns``   the exact columns hashed.  SQL expressions are allowed so
                  a volatile representation can be normalized (rounding a
                  REAL, collapsing an error string to a boolean).
    ``where``     OPTIONAL row-scope predicate.  See ROW SCOPE below.
    ``scope_of``  when this spec is a SLICE of another resource, the key of
                  the whole-table resource it slices.  Reporting only; it
                  lets the page say "G reads 723 of the 330,865 frames".
    ``why``       the auditor-facing justification: which downstream stage
                  reads these columns, and which columns were deliberately
                  LEFT OUT because they churn without meaning.

    ROW SCOPE — why ``where`` exists
    --------------------------------
    Choosing the right COLUMNS is only half of the fingerprint design; the
    other half is choosing the right ROWS.  A whole-table fingerprint over
    ``frames`` answers "did anything anywhere in the archive change?", and
    the honest answer to that is almost always yes.  Every consumer then
    goes stale together, so the cheapest true correction — repairing one
    frame's geometry — mandates re-running the grism extraction and the
    photometry prototype, neither of which can possibly be affected.  A
    gate that expensive is a gate people route around.

    So a stage that reads a NARROW, DECLARABLE slice reads a scoped
    resource instead: ``table:frames@grism`` hashes the T CrB / tet CrB
    grism frames, ``table:frames@s4proto`` the two prototype targets.  The
    predicate MUST be a superset of what the stage actually consumes — a
    slice that is too narrow would hide a real change, which is the one
    error this module may never make.  Each scoped spec's ``why`` states
    the consuming query it was widened from.
    """

    key: str
    kind: str
    name: str
    why: str
    order_by: str = ""
    columns: tuple[str, ...] = ()
    database: Optional[str] = None
    where: str = ""
    scope_of: str = ""


def _t(key: str, name: str, order_by: str, columns: Sequence[str],
       why: str, where: str = "", scope_of: str = "") -> ResourceSpec:
    """Shorthand constructor for a manifest-table resource."""
    return ResourceSpec(key=key, kind="table", name=name, order_by=order_by,
                        columns=tuple(columns), why=why, where=where,
                        scope_of=scope_of)


# --- the resource registry -------------------------------------------------
# Read this list as the answer to "what would have to change for stage X to
# be worth re-running?".  Every entry states it in words.
RESOURCES: dict[str, ResourceSpec] = {}


def _reg(spec: ResourceSpec) -> ResourceSpec:
    """Register a spec, refusing duplicates (two specs for one key would
    make the digest depend on import order)."""
    if spec.key in RESOURCES:
        raise ProvenanceError(f"duplicate resource spec: {spec.key}")
    RESOURCES[spec.key] = spec
    return spec


# ---- S0's external input --------------------------------------------------
# The one input that comes from outside this repo: the header-scan catalog
# on the ASTRO archive drive.  Declaring it is not bookkeeping — it is the
# only way the machinery can see the 2026-08 geometry re-characterization,
# which happened IN THIS FILE (the rescue pass rewrote its NAXIS values) and
# left every stage keyed on frames.era_id downstream of a changed input.
_reg(ResourceSpec(
    key="stat:rlmt-catalog",
    kind="stat",
    name="/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite",
    why=("S0's sole input: the 330,865-row FITS header scan.  Fingerprinted "
         "by (size, mtime) rather than by content because it is 400 MB on a "
         "spinning archive drive and sibling processes rewrite it while a "
         "status check runs — a content hash would be both slow and "
         "irreproducible.  This is an integrity SURROGATE, in the same "
         "spirit as S0c's size_bytes note: it detects 'the catalog was "
         "rewritten', which is the only question S0's freshness turns on, "
         "and it cannot detect an in-place edit preserving size and mtime."))

)

# ---- S0 outputs -----------------------------------------------------------
#: The frames columns hashed, shared by the whole-table resource and by every
#: scoped slice of it.  Defined ONCE so a slice can never drift from the whole
#: — two column lists for one table would let a change be visible in one
#: fingerprint and invisible in the other.
_FRAME_COLUMNS = (
    "obs_rowid", "tree", "is_canonical", "era_id",
    "CAST(naxis1 AS INTEGER)", "CAST(naxis2 AS INTEGER)",
    "lower(trim(coalesce(filter,'')))", "coalesce(imagetyp,'')",
    "round(coalesce(exptime,-1),4)", "coalesce(canonical_target,'')",
    "coalesce(target_key,'')", "coalesce(night,'')",
    "coalesce(pltsolvd,-1)", "CASE WHEN error IS NULL THEN 1 ELSE 0 END")

_reg(_t(
    "table:frames", "frames", "obs_rowid",
    # obs_rowid: row identity — a re-numbered catalog IS a different input.
    # tree/is_canonical: every downstream selection starts with these two.
    # era_id/naxis1/naxis2: the geometry artifact lives exactly here (S1's
    #   window gate, S0b's era-matched calibration spec, S2's mode groups).
    # filter: S1's grism gate, S0c's whitelists, the whole G track.
    # imagetyp: the Light-frame gate shared by S1/S2/S0c.
    # exptime: S0b's dark-exposure matching and S2's linearity ladders.
    # canonical_target/target_key: every per-project selection.
    # night: era spans, cadence, the calibration shopping list's date rows.
    # pltsolvd: S1's "the header already claims a solve" gate.
    # error IS NULL is hashed as a BOOLEAN, not as its text: the rescue pass
    #   rewrote the message strings without changing which frames are usable.
    _FRAME_COLUMNS,
    "The frame-level facts every downstream stage branches on.  LEFT OUT: "
    "path/basename/size (filesystem-scan detail, not a decision), "
    "pointing_offset_deg and qc_flags (recomputed floats that churn at the "
    "1e-12 level on every rebuild), and the raw error TEXT (normalized to a "
    "boolean above).  Consumers that read the WHOLE table: S0b, S0c, S1, "
    "S1b, S2, S2c, S3.  Consumers that read a narrow slice use the scoped "
    "resources below instead."))

# ---- scoped slices of frames ----------------------------------------------
# Each predicate is deliberately WIDER than the consuming query (see the ROW
# SCOPE note on ResourceSpec): wider can only cost an unnecessary re-run,
# narrower could hide a real change.
_reg(_t(
    "table:frames@grism", "frames", "obs_rowid", _FRAME_COLUMNS,
    "G's slice: the T CrB / tet CrB grism frames its worklist is drawn from "
    "(723 of 330,865 rows today).  Widened from run_g_tcrb_validation's "
    "build_worklist, which additionally filters on pointing_offset_deg — a "
    "recomputed float this module refuses to hash, so the slice ignores it "
    "and keeps every candidate frame in scope.  NOTE this slice DOES contain "
    "the 54 tet CrB calibrator frames sitting in phantom eras, so the "
    "geometry repair correctly marks G stale; the scoping is not a "
    "whitewash, it is a narrower true statement.",
    where=("is_canonical = 1 "
           "AND lower(trim(coalesce(filter,''))) IN ('hrg','lrg') "
           "AND coalesce(target_best,'') LIKE '%CrB%'"),
    scope_of="table:frames"))

_reg(_t(
    "table:frames@s4proto", "frames", "obs_rowid", _FRAME_COLUMNS,
    "S4's slice: the two prototype targets (TARGETS = ('anuma','vvpup') in "
    "build_s4_photometry.py), 2,632 of 330,865 rows.  Zero of them are in a "
    "phantom era, so an S0 rebuild that only repairs geometry leaves this "
    "digest identical — which is the correct answer, and the whole-table "
    "fingerprint could not give it.",
    where="is_canonical = 1 AND target_key IN ('anuma','vvpup')",
    scope_of="table:frames"))

_reg(_t(
    "table:frames@cv", "frames", "obs_rowid", _FRAME_COLUMNS,
    "CV-S4's slice: the five staged CV targets, 9,083 rows, of which 207 "
    "(EU UMa) ARE in phantom era 80.  Contrast with the S4 slice above: same "
    "columns, same rule, opposite answer — which is exactly the "
    "discrimination a whole-table fingerprint cannot make.",
    where=("is_canonical = 1 AND target_key IN "
           "('anuma','vvpup','stlmi','yzcnc','euuma')"),
    scope_of="table:frames"))

_reg(_t(
    "table:eras", "eras", "era_id",
    ["era_id", "coalesce(readoutm,'')", "naxis1", "naxis2", "xbinning",
     "coalesce(egain,-1)", "n_frames", "coalesce(first_night,'')",
     "coalesce(last_night,'')"],
    "The era registry is small and every column is consumed: S0b keys the "
    "calibration spec on (era, readoutm, geometry), S2 groups readout modes, "
    "the ops shopping list quotes era spans verbatim.  Nothing is left out."))

_reg(_t(
    "table:aliases", "aliases", "raw_name",
    ["raw_name", "coalesce(canonical_target,'')", "coalesce(target_key,'')"],
    "Only the name->canonical mapping is consumed downstream.  LEFT OUT: "
    "n_frames and cone_check_passed (derived counts that move whenever any "
    "unrelated frame is added, without changing what any name MEANS)."))

_reg(_t(
    "table:project_counts", "project_counts", "project, target, metric",
    ["project", "target", "metric", "claimed_frames", "manifest_frames",
     "diff_frames"],
    "The strategy-claim reconciliation the five ANALYSIS_STRATEGY.md files "
    "quote.  LEFT OUT: the free-text source column (prose, not a number)."))

# ---- S0b outputs ----------------------------------------------------------
_reg(_t(
    "table:raw_reduced_links", "raw_reduced_links", "reduced_rowid",
    ["reduced_rowid", "coalesce(raw_rowid,-1)", "match_method"],
    "S2's reconstruction experiment selects eras by link count and match "
    "method.  LEFT OUT: jd_drift_s (a float recomputed per build) and the "
    "path strings (redundant with the rowids)."))

_reg(_t(
    "table:calib_frames", "calib_frames", "obs_rowid",
    ["obs_rowid", "coalesce(era_id,-1)", "kind", "round(coalesce(exptime_bin,-1),4)",
     "lower(coalesce(filter,''))", "is_master"],
    "What calibration exists, in which era, of which kind — the input to "
    "coverage, to S0c's era_exact calibration staging, and to G's master-dark "
    "path.  LEFT OUT: path/night/camtemp (not read by any consumer's rule)."))

_reg(_t(
    "table:calib_coverage", "calib_coverage", "era_id, req_kind, coalesce(req_key,'')",
    ["era_id", "req_kind", "coalesce(req_key,'')", "n_science", "n_calib_raw",
     "n_calib_master", "spec_n", "status"],
    "The per-era pass/fail the S0b report and the observatory request read."))

_reg(_t(
    "table:calib_gaps", "calib_gaps",
    "era_id, need_kind, spec",
    ["era_id", "coalesce(camera,'')", "need_kind", "spec", "have_raw",
     "have_master", "status", "n_science_frames_blocked",
     "coalesce(projects_affected,'')"],
    "THE artifact-bearing table: ops/2026-08_observatory_request.md quotes "
    "these rows by number.  Every column in it is quoted, so every column is "
    "hashed."))

# ---- S0c outputs ----------------------------------------------------------
for _proj, _tbl in (("bestar", "stage_bestar_grism"),
                    ("cv", "stage_cv_timeseries"),
                    ("dwarf", "stage_dwarfgalaxy_agn_survey"),
                    ("sn2023ixf", "stage_sn2023ixf_lightcurve"),
                    ("tcrb", "stage_tcrb_monitoring")):
    _reg(_t(
        f"table:{_tbl}", _tbl, "path, role",
        ["path", "role", "coalesce(era_id,-1)",
         "lower(coalesce(filter,''))", "coalesce(match_basis,'')"],
        "Which file is staged in which role for this project.  LEFT OUT: "
        "stage_build_id — it embeds the build TIMESTAMP, so hashing it would "
        "mark this table changed on every rebuild even when the selection is "
        "byte-identical.  That column is the canonical example of the "
        "volatile-column trap this module exists to avoid."))

_reg(_t(
    "table:s0c_stage_files", "s0c_stage_files", "project",
    ["project", "stage_table", "n_rows", "n_science", "n_calib", "n_cone",
     "n_eras", "selection_rule"],
    "The staging summary the S0c report and each project README quote, "
    "including the selection_rule text — a changed rule is a changed input "
    "even when the row count happens to match."))

# ---- S1 outputs (experiment + batch) --------------------------------------
_reg(_t(
    "table:s1_solve_experiment", "s1_solve_experiment", "obs_rowid",
    ["obs_rowid", "stratum_id", "coalesce(status,'')", "coalesce(fail_kind,'')"],
    "Per-frame solve outcomes behind the stratum verdicts.  LEFT OUT: "
    "solve_time_s and log tails (wall-clock noise)."))
_reg(_t(
    "table:s1_strata", "s1_strata", "stratum_id",
    ["stratum_id"],
    "The stratum registry the S1 verdict table is built from.  Column list "
    "is deliberately minimal because the table is DESTROYED: the spec exists "
    "so the machinery can report it MISSING rather than skip it silently."))
_reg(_t(
    "table:s1_populations", "s1_populations", "rowid",
    ["rowid"],
    "Population roll-up behind the GO/CAUTION/NO-GO verdicts.  DESTROYED; "
    "spec kept so the loss is reported, not silently ignored."))
_reg(_t(
    "table:s1_failure_autopsy", "s1_failure_autopsy", "rowid",
    ["rowid"],
    "Failure-mode autopsy behind the S1 report's diagnosis section. "
    "DESTROYED; spec kept so the loss is reported."))

_reg(_t(
    "table:s1_batch", "s1_batch", "obs_rowid",
    ["obs_rowid", "stratum_id", "population", "qc_gated", "status",
     "coalesce(fail_kind,'')"],
    "Which frames were queued, in which stratum, and how each finished. "
    "LEFT OUT: solve_time_s, log_tail, finished_utc, and the solved WCS "
    "numbers — a re-solve of the same frame must not read as a new input."))

# ---- S2 outputs -----------------------------------------------------------
_reg(_t(
    "table:s2_ceiling_modes", "s2_ceiling_modes", "mode",
    ["mode", "n_frames", "hard_max_adu", "clip_adu", "veto_adu", "bits"],
    "The adopted per-mode ADU ceiling and saturation veto every photometric "
    "stage inherits.  DESTROYED by the S0 table swap."))
_reg(_t(
    "table:s2_ptc_fits", "s2_ptc_fits", "mode, kind",
    ["mode", "kind", "round(coalesce(gain_e_per_adu,-1),6)",
     "round(coalesce(read_noise_e,-1),6)", "n_points"],
    "Gain and read noise — the error model's foundation.  Rounded to 1e-6 so "
    "float formatting cannot masquerade as a measurement change.  DESTROYED."))
_reg(_t(
    "table:s2_recon_eras", "s2_recon_eras", "era_id",
    ["era_id", "mode", "n_links", "n_pairs_used",
     "round(coalesce(flat_median,-1),6)", "round(coalesce(dark_median,-1),6)"],
    "The reconstructed effective flat/dark per era.  DESTROYED."))
_reg(_t(
    "table:s2_linearity_ladders", "s2_linearity_ladders", "ladder_id",
    ["ladder_id", "mode", "n_rungs", "n_frames",
     "round(coalesce(max_abs_resid_pct,-1),4)"],
    "Linearity residuals per archival exposure ladder.  DESTROYED."))
_reg(_t(
    "table:detector_params", "detector_params", "era_group, quantity",
    ["era_group", "quantity", "round(coalesce(value,-1),6)",
     "round(coalesce(uncertainty,-1),6)", "method"],
    "The S2 deliverable every later error budget cites by name.  DESTROYED."))

# The adopted S2 constants as they exist TODAY.  They are not in the database
# at all: when the S0 rebuild replaced the manifest that carried the s2_*
# tables, the ceiling and veto numbers survived only as literals copied into
# macro_phot/series.py, and both photometry products read them from there.
# Declaring the FILE is the only way the graph can see a change to them — and
# the only way S4 can be prevented from reporting FRESH while the measurement
# underneath it is both destroyed and, for two modes, known to be suspect.
_reg(ResourceSpec(
    key="file:pipeline/macro_phot/series.py", kind="file",
    name="pipeline/macro_phot/series.py",
    why=("The adopted per-mode ADU ceiling and saturation veto "
         "(S2_MODE_CEILING_ADU / S2_MODE_VETO_ADU) that S4 and CV-S4 apply "
         "to every frame, plus the nominal gain.  Hashed WHOLE because the "
         "constants are module-level literals with no separate home; an "
         "unrelated edit to this file will therefore read as a changed "
         "input.  That is the conservative direction: a spurious re-run "
         "costs time, a missed one publishes a wrong error bar.")))

# ---- S2c outputs (per-frame dispersion measurement) ------------------------
# The evidence EVERY 'the filter label is not a spectrum test' verdict is
# gated on.  Undeclared, it fills in silently while the five strategy
# documents keep quoting counts the measurement is in the middle of
# overturning — the exact drift this module exists to end.
_reg(_t(
    "table:frame_dispersion", "frame_dispersion", "obs_rowid",
    ["obs_rowid", "lower(coalesce(filter,''))", "coalesce(canonical_target,'')",
     "coalesce(verdict,'')", "coalesce(strength_class,'')",
     "coalesce(status,'')", "round(coalesce(median_ab,-1),3)",
     "coalesce(code_version,'')"],
    "Per-frame direct/dispersed/indeterminate verdicts from source "
    "elongation.  The verdict and its strength class are what a strategy "
    "document acts on; median_ab is rounded to 1e-3 because it is a measured "
    "float.  LEFT OUT: measured_at and measure_s (wall-clock), and every "
    "intermediate detection count (derivable, and noisy per re-measurement). "
    "code_version IS hashed: a reclassified verdict from the same pixels is "
    "still a changed input to a human decision."))

# ---- S3 outputs -----------------------------------------------------------
_reg(_t(
    "table:frame_times", "frame_times", "path",
    # bjd_tdb rounded to 1e-7 d = 8.6 ms — far finer than any claim made on
    # it, coarse enough that a re-run with an identical ephemeris re-hashes
    # identically instead of drifting in the last binary digit.
    ["path", "coalesce(era_id,-1)", "round(coalesce(bjd_tdb,-1),7)",
     "coalesce(mid_method,'')", "coalesce(bjd_method,'')",
     "coalesce(start_evidence,'')"],
    "The per-frame BJD_TDB and the evidence class behind it — what every "
    "time-series paper stands on.  LEFT OUT: sibling_jd_drift_s and the "
    "intermediate correction terms (derivable, and float-noisy)."))
_reg(_t(
    "table:s3_header_audit", "s3_header_audit", "path",
    ["path", "coalesce(era_id,-1)", "naxis1", "naxis2",
     "round(coalesce(corner_ltt_s,-1),6)", "coalesce(pixscale_source,'')"],
    "The header-semantics + field-geometry audit.  Its naxis columns come "
    "from re-reading each FITS header, NOT from frames — which is why this "
    "table already carries TRUE geometry while frames does not."))
_reg(_t(
    "table:s3_clock_drift", "s3_clock_drift", "rowid",
    ["rowid"],
    "Clock-drift evidence rows; identity-only fingerprint (row count is the "
    "meaningful signal, the residuals are float-noisy by construction)."))
_reg(_t(
    "table:s3_cadence", "s3_cadence", "rowid",
    ["rowid"],
    "Cadence roll-up; identity-only for the same reason as s3_clock_drift."))

# ---- S4 outputs (a sibling database, not the manifest) ---------------------
_reg(ResourceSpec(
    key="db:phot:phot_selection", kind="db", name="phot_selection",
    database="products/phot/anuma_vvpup_prototype.sqlite",
    order_by="target_key, era_id",
    columns=("target_key", "era_id", "n_raw", "n_linked"),
    why=("Which (target, era) sets the photometry prototype actually "
         "measured, and how many frames each drew — the single fact that "
         "decides whether a geometry or filter re-characterization can touch "
         "S4 at all.  (Earlier this spec hashed `rowid` alone, which made it "
         "a row COUNT wearing a content hash's clothes: the selection could "
         "have changed from era 76 to era 80 without the digest moving.)")))
_reg(ResourceSpec(
    key="db:phot:phot_series", kind="db", name="phot_series",
    database="products/phot/anuma_vvpup_prototype.sqlite",
    order_by="rowid", columns=("count",),
    why=("Row count of the produced light curves.  Counted rather than "
         "hashed: the table is ~10^6 rows of floats on a spinning disk, and "
         "a full hash would cost minutes per status check for no extra "
         "discrimination — the selection fingerprint above already carries "
         "the identity of what was measured.")))

# ---- CV-S4 outputs (the OTHER photometry product) --------------------------
# products/phot/cv_timeseries.sqlite is the only photometry product that
# actually contains phantom-era frames — 207 EU UMa frames in era 80, plus 75
# in era 78, which the geometry repair proves are ONE configuration split by
# the artifact.  It was built at 19:08Z on 2026-08-18, after the S0 rebuild
# and during the geometry rescue.  Until it was declared here, `status` could
# return 0 — a green light to publish — with this product resting on frames
# whose era assignment is known wrong.
_reg(ResourceSpec(
    key="db:cvphot:cv_selection", kind="db", name="cv_selection",
    database="products/phot/cv_timeseries.sqlite",
    order_by="series_key",
    columns=("series_key", "target_key", "era_id", "filter",
             "coalesce(readoutm,'')", "coalesce(provenance,'')"),
    why=("Which (target, era, filter) series the CV product measured, and "
         "the pixel provenance chosen for each.  This is the fingerprint "
         "that carries the artifact: 'euuma|e80|g' and 'euuma|e78|g' are two "
         "rows here for one camera configuration, and cv_build_meta records "
         "the rule 'one pixel provenance per (target, era); never mixed "
         "inside a series' — a rule being enforced on a phantom boundary.")))
_reg(ResourceSpec(
    key="db:cvphot:cv_frames", kind="db", name="cv_frames",
    database="products/phot/cv_timeseries.sqlite",
    order_by="rowid", columns=("count",),
    why=("Row count of the per-frame CV measurements.  Counted rather than "
         "hashed, for the same reason as phot_series: the identity of what "
         "was measured lives in cv_selection above, and a full hash of a "
         "400 MB product on a spinning disk would cost minutes per check.")))

# ---- CV-S5 outputs (the characterization product) --------------------------
# products/phot/cv_characterization.sqlite answers "what can this data set
# actually measure?" — image quality, the measured noise floor, the spectral
# window, injection-recovery contours and the epoch-precision Monte Carlo.
# Every science verdict in docs/CV_TimeSeries/cv_characterization.html is
# recomputed from it, so a stale characterization is a stale set of verdicts:
# it has to sit inside the DAG or `status` could green-light a paper whose
# detection limits were measured against superseded photometry.
_reg(ResourceSpec(
    key="db:cvchar:ch_noise_series", kind="db", name="ch_noise_series",
    database="products/phot/cv_characterization.sqlite",
    order_by="series_key",
    columns=("series_key", "coalesce(readoutm,'')",
             "printf('%.6f', coalesce(floor_plateau,-1))",
             "printf('%.6f', coalesce(prec_at_target,-1))",
             "coalesce(n_target_points,-1)"),
    why=("The measured noise floor and the achieved precision at each "
         "target's own brightness, per series.  Hashed rather than counted "
         "because these ARE the numbers the science verdicts turn on: a "
         "changed floor changes what the paper may claim, and a row count "
         "would not notice.")))
_reg(ResourceSpec(
    key="db:cvchar:ch_frames", kind="db", name="ch_frames",
    database="products/phot/cv_characterization.sqlite",
    order_by="rowid", columns=("count",),
    why=("Row count of the per-frame quality characterization.  Counted "
         "rather than hashed for the same reason as cv_frames: the identity "
         "of what was characterized lives in ch_noise_series above.")))
_reg(ResourceSpec(
    key="db:cvchar:ch_cuts", kind="db", name="ch_cuts",
    database="products/phot/cv_characterization.sqlite",
    order_by="axis",
    columns=("axis", "printf('%.6f', coalesce(threshold,-1))",
             "n_pass", "n_fail"),
    why=("The quality thresholds and what they reject.  Hashed because "
         "these cuts are now APPLIED to every measurement downstream — "
         "noise, cadence, detection and timing all read `usable` — so a "
         "changed threshold silently changes every number on the page.  "
         "For one production run the cut was computed and read by nothing, "
         "which is exactly the failure this entry makes visible.")))
_reg(ResourceSpec(
    key="db:cvchar:ch_contour", kind="db", name="ch_contour",
    database="products/phot/cv_characterization.sqlite",
    order_by="scope, regime, score, period_d",
    columns=("scope", "regime", "score",
             "printf('%.6f', period_d)",
             "printf('%.6f', coalesce(amp90,-1))"),
    why=("The detection contours, keyed by SCORE MODE.  The score column is "
         "part of the identity: the same injections scored as blind period "
         "determination rather than as detection at a known period give "
         "answers 3-8x apart, and quoting one under the other's name is the "
         "defect this column exists to prevent.")))
_reg(ResourceSpec(
    key="db:cvchar:ch_verdict", kind="db", name="ch_verdict",
    database="products/phot/cv_characterization.sqlite",
    order_by="rank",
    columns=("goal_id", "verdict"),
    why=("The goal-by-goal verdicts themselves.  Cheapest possible hash — "
         "the id and the one word — because that word is the deliverable: "
         "if a rebuild moves a goal from SUPPORTED to NOT SUPPORTED, every "
         "page and plan that cites it is stale and `status` must say so.")))

# ---- CV-S6 outputs (the catalogue tie) ------------------------------------
# The stage that converts the CV light curves from an arbitrary internal
# gauge to natural-system magnitudes on a catalogue zero point.  Everything
# ABOVE it in this file measures; this one decides what a published
# magnitude MEANS, so a stale tie is a page full of magnitudes on a zero
# point that no longer exists.
_reg(ResourceSpec(
    key="db:cvphot:cv_cat_fetch", kind="db", name="cv_cat_fetch",
    database="products/phot/cv_timeseries.sqlite",
    order_by="catalogue, field_key",
    columns=("catalogue", "field_key", "n_rows",
             "coalesce(cache_sha256,'')"),
    why=("WHICH catalogue rows the tie was solved against, identified by "
         "the sha256 of the cached pull rather than by its date.  A "
         "re-pull that returns the same bytes is not a change; a re-pull "
         "that returns different bytes changes every zero point "
         "downstream, and the date alone could not tell the two apart.")))
_reg(ResourceSpec(
    key="db:cvphot:cv_cat_astrom", kind="db", name="cv_cat_astrom",
    database="products/phot/cv_timeseries.sqlite",
    order_by="catalogue, target_key, era_id",
    columns=("catalogue", "target_key", "era_id",
             "printf('%.3f', coalesce(dra_arcsec, -999))",
             "printf('%.3f', coalesce(ddec_arcsec, -999))",
             "coalesce(applied, -1)"),
    why=("The astrometric zero point measured for each block against the "
         "catalogue, and whether it was REMOVED.  Load-bearing because it "
         "changes which catalogue source each comparison star is matched "
         "to: EU UMa era 78 sat 5.2 arcsec away, matched nothing, and was "
         "published as untieable until this was measured.  A rebuild that "
         "starts or stops applying an offset changes the tie stars, hence "
         "the zero point, hence every absolute magnitude of that block, "
         "while leaving the row count identical.")))
_reg(ResourceSpec(
    key="db:cvphot:cv_cattie", kind="db", name="cv_cattie",
    database="products/phot/cv_timeseries.sqlite",
    order_by="series_key, catalogue, band",
    columns=("series_key", "catalogue", "band",
             "printf('%.6f', coalesce(zp, -999))",
             "printf('%.6f', coalesce(colour_term, -999))",
             "coalesce(verdict,'')"),
    why=("The zero point, the colour term and the verdict per block — the "
         "three numbers that decide what a published magnitude is and "
         "whether it may be published at all.  Hashed, not counted: a "
         "re-solve that moves a zero point by 30 mmag leaves the row "
         "count identical and invalidates every absolute magnitude, every "
         "cross-survey comparison and the duty-cycle goal that depends on "
         "them.  n_fit and the check statistics are deliberately LEFT OUT "
         "— they are diagnostics of the same fit, and a changed fit "
         "always moves zp or colour_term first.")))

# ---- G outputs ------------------------------------------------------------
_reg(_t(
    "table:g_extractions", "g_extractions", "obs_rowid, method",
    ["obs_rowid", "method", "coalesce(filter,'')", "coalesce(era_id,-1)",
     "coalesce(role,'')", "coalesce(gate_verdict,'')",
     "coalesce(anchor_status,'')", "coalesce(status,'')"],
    "Which frames were extracted, how, and whether the identity gate passed. "
    "LEFT OUT: snippet_json and every measured float (the spectra "
    "themselves) — a re-extraction that reaches the same verdicts on the "
    "same frames is not a changed input to anything downstream."))

# ---- published artifacts (files) ------------------------------------------
def _f(key: str, path: str, why: str) -> None:
    _reg(ResourceSpec(key=key, kind="file", name=path, why=why))


_f("file:docs/pipeline/s0_manifest.html", "docs/pipeline/s0_manifest.html",
   "S0 chain-of-evidence page: published dedup, alias, era and pointing "
   "numbers.  Whole-file hash — every byte of a published page is the claim.")
_f("file:docs/pipeline/s0b_calibration_inventory.html",
   "docs/pipeline/s0b_calibration_inventory.html",
   "S0b inventory page — the source the October observatory request quotes "
   "its era-resolved shopping list from.  Whole-file hash.")
_f("file:docs/pipeline/s0c_staging.html", "docs/pipeline/s0c_staging.html",
   "S0c staging page: the per-project selection rules and staged counts each "
   "ANALYSIS_STRATEGY.md is diffed against.  Whole-file hash.")
_f("file:docs/pipeline/s1_astrometry.html", "docs/pipeline/s1_astrometry.html",
   "S1 astrometry verdict page: the GO/CAUTION/NO-GO calls per stratum, "
   "including the frames written off as unsolvable.  Whole-file hash.")
_f("file:docs/pipeline/s2_detector.html", "docs/pipeline/s2_detector.html",
   "S2 detector-truth page: ceiling, gain, read noise and linearity — the "
   "numbers every later error budget inherits.  Whole-file hash.")
_f("file:docs/pipeline/s3_timing.html", "docs/pipeline/s3_timing.html",
   "S3 timing page: BJD_TDB method, clock drift and the frame-center "
   "caveat.  Whole-file hash.")
_f("file:docs/pipeline/s4_photometry.html", "docs/pipeline/s4_photometry.html",
   "S4 photometry page: the AN UMa / VV Pup prototype light curves and "
   "their empirical error model.  Whole-file hash.")
_f("file:docs/CV_TimeSeries/cv_characterization.html",
   "docs/CV_TimeSeries/cv_characterization.html",
   "CV-S5 characterization page: measured image quality, noise floor, "
   "spectral window, detectability contours, epoch precision, and the "
   "goal-by-goal verdict on ANALYSIS_STRATEGY.md.  Whole-file hash.")
_f("file:docs/CV_TimeSeries/cv_catalogue_tie.html",
   "docs/CV_TimeSeries/cv_catalogue_tie.html",
   "CV-S6 catalogue-tie page: which catalogue, which stars carried the tie, "
   "the zero point and colour term per block, the accuracy achieved on "
   "independent check stars, the cross-catalogue systematic, and the list "
   "of blocks left RELATIVE with the reason for each.  Whole-file hash — "
   "every byte of a published page is the claim.")
_f("file:pipeline/macro_phot/cattie.py", "pipeline/macro_phot/cattie.py",
   "The catalogue-tie arithmetic: the cleanliness gate, the robust "
   "ZP + colour-term fit, the colour-range rules and every verdict "
   "threshold.  A change here changes what a published magnitude is, so it "
   "is an INPUT to the tie stage and to its page.")
_f("file:pipeline/macro_phot/characterize.py",
   "pipeline/macro_phot/characterize.py",
   "The pure characterization arithmetic: the detector constants quoted from "
   "S2, the degradation factor that sets every quality cut, the false-alarm "
   "probability and the recovery level.  Changing any of them changes every "
   "verdict on the CV-S5 page, so the page is stale the moment this file "
   "moves.")
_f("file:docs/pipeline/s0e_geometry_fix.html",
   "docs/pipeline/s0e_geometry_fix.html",
   "S0e geometry-repair page: which frames carried BINTABLE geometry, which "
   "were repaired, and which small-NAXIS eras are genuine subframes.  Whole-"
   "file hash.")
_f("file:docs/pipeline/s2c_filter_identity.html",
   "docs/pipeline/s2c_filter_identity.html",
   "S2c filter-identity page: the per-frame dispersion verdicts that decide "
   "which FILTER labels denote a spectrum.  Whole-file hash.")
_f("file:ops/2026-08_observatory_request.md",
   "ops/2026-08_observatory_request.md",
   "The October shopping list, quoting calib_gaps rows by number.")
for _p in ("BeStar_Grism", "CV_TimeSeries", "DwarfGalaxy_AGN_Survey",
           "SN2023ixf_LightCurve", "TCrB_Monitoring"):
    _f(f"file:{_p}/ANALYSIS_STRATEGY.md", f"{_p}/ANALYSIS_STRATEGY.md",
       f"{_p} strategy document — quotes staged counts and filter rules.")
# ROADMAP.md and the public project pages quote stage status and per-target
# counts.  They were previously given a verdict ON the status page while
# being invisible TO its exit code — a page issuing judgements its own gate
# could not enforce.  Declaring them closes that gap.
_f("file:ROADMAP.md", "ROADMAP.md",
   "The roadmap quotes stage status and per-project frame counts drawn from "
   "project_counts and s0c_stage_files.  It contains no independent "
   "measurement, so every number in it moves when those tables move.")
for _p in ("BeStar_Grism", "CV_TimeSeries", "DwarfGalaxy_AGN_Survey",
           "Legacy_Rigel", "SN2023ixf_LightCurve", "TCrB_Monitoring"):
    _f(f"file:docs/{_p}/index.html", f"docs/{_p}/index.html",
       f"Public {_p} page — hand-authored, quoting staged counts and "
       f"per-target totals.  Whole-file hash.")
_f("file:docs/index.html", "docs/index.html",
   "The public index, linking and summarising every project page.  Declared "
   "so the gate covers it; this module never writes it.")


# ===========================================================================
# 2.  STAGES — the DAG
# ===========================================================================

@dataclass(frozen=True)
class Stage:
    """One pipeline stage.

    ``code_version`` is the CONSTANT the build script stamps into its
    ``*_build_meta`` table.  It is part of the freshness test: a stage whose
    code version moved since its recorded run is stale even if every input
    byte is identical, because the rules that read those bytes changed.

    ``build_cmd`` is the command (or commands, one per line) that ACTUALLY
    re-runs the stage.  It is not documentation: it is printed by ``plan``,
    by ``status`` and on the published page, and a person pastes it into a
    terminal.  Every line is therefore checked by
    ``test_provenance.py::test_every_build_command_is_runnable``, which runs
    each one with ``--help`` appended and requires exit 0 — so a command
    that argparse would reject cannot reach the plan.

    ``version_file`` / ``version_symbol`` give the code version for stages
    whose constant lives in a SCRIPT rather than an importable package
    (CV-S4, S2c).  The constant is read from the source text with
    :func:`read_version_constant` rather than by importing the script, which
    would execute a module another workflow may be editing.
    """

    key: str
    title: str
    code_version: str
    reads: tuple[str, ...]
    writes: tuple[str, ...]
    build_cmd: str
    meta_table: Optional[str] = None
    meta_version_key: str = "code_version"
    note: str = ""
    version_file: str = ""
    version_symbol: str = ""

    @property
    def build_lines(self) -> tuple[str, ...]:
        """``build_cmd`` split into individual commands, blanks dropped."""
        return tuple(line.strip() for line in self.build_cmd.splitlines()
                     if line.strip())

    @property
    def hand_authored(self) -> bool:
        """True for stages a person writes rather than a script produces.

        These are the ones whose ``record`` needs a human to say what was
        re-derived: nothing else can tell whether the work happened.
        """
        return self.code_version.startswith("(")


# The code-version constants live in the packages that own each stage.  They
# are imported lazily (inside the function) so that importing provenance.py
# never drags in numpy/astropy through a sibling package's __init__.
def _code_versions() -> dict[str, str]:
    """Current code-version constants, read from the modules that own them.

    Import failures are not swallowed: a stage whose version constant cannot
    be found is a broken declaration, and the caller must see it.
    """
    from macro_core import (S0_CODE_VERSION, S0B_CODE_VERSION,
                            S0C_CODE_VERSION)
    from macro_core.astrom import S1_CODE_VERSION
    from macro_core.batch import S1B_CODE_VERSION
    from macro_core.timing import S3_CODE_VERSION
    from macro_phot import S4_CODE_VERSION
    from rlmt_diagnostics import S2_CODE_VERSION
    from macro_grism.gate import G_CODE_VERSION
    return {"S0": S0_CODE_VERSION, "S0b": S0B_CODE_VERSION,
            "S0c": S0C_CODE_VERSION, "S1": S1_CODE_VERSION,
            "S1b": S1B_CODE_VERSION, "S2": S2_CODE_VERSION,
            "S3": S3_CODE_VERSION, "S4": S4_CODE_VERSION,
            "G": G_CODE_VERSION}


#: Commands that are typed by a person rather than run by a script.  They are
#: recognised by their leading '(' and skipped by the runnability test.
HAND_AUTHORED_CMD_PREFIX = "("

STAGES: tuple[Stage, ...] = (
    Stage(
        key="S0e", title="Geometry repair of the external catalog",
        code_version="(external tool)",
        reads=(),
        writes=("stat:rlmt-catalog",
                "file:docs/pipeline/s0e_geometry_fix.html"),
        build_cmd="python pipeline/scripts/rescan_geometry.py run",
        note="Rewrites NAXIS1/2 in the external header-scan catalog for "
             "frames whose rescued values were the tile-compressed "
             "BINTABLE's row length and row count.  Declared as a STAGE, "
             "not as an unexplained external mtime change: the catalog "
             "rewrite is the event that invalidated S0, and a graph that "
             "cannot name its cause sends the reader back to memory."),
    Stage(
        key="S0", title="Manifest (dedup, aliases, eras, pointing)",
        code_version="S0_CODE_VERSION",
        reads=("stat:rlmt-catalog",),
        writes=("table:frames", "table:frames@grism", "table:frames@s4proto",
                "table:frames@cv", "table:eras", "table:aliases",
                "table:project_counts"),
        build_cmd="python pipeline/scripts/build_s0_manifest.py",
        meta_table="build_meta",
        note="Reads the external rlmt-catalog.sqlite; every geometry value "
             "in frames.naxis1/2 arrives from that catalog.  It writes the "
             "scoped frame slices too — they are views of the same table, "
             "so the same stage produces them."),
    Stage(
        key="S0b", title="Calibration inventory (links, coverage, gaps)",
        code_version="S0B_CODE_VERSION",
        reads=("table:frames", "table:eras"),
        writes=("table:raw_reduced_links", "table:calib_frames",
                "table:calib_coverage", "table:calib_gaps"),
        build_cmd="python pipeline/scripts/build_s0b_inventory.py",
        meta_table="s0b_build_meta"),
    Stage(
        key="S0c", title="Per-project staging manifests",
        code_version="S0C_CODE_VERSION",
        reads=("table:frames", "table:eras", "table:calib_frames"),
        writes=("table:stage_bestar_grism", "table:stage_cv_timeseries",
                "table:stage_dwarfgalaxy_agn_survey",
                "table:stage_sn2023ixf_lightcurve",
                "table:stage_tcrb_monitoring", "table:s0c_stage_files"),
        build_cmd="python pipeline/scripts/build_s0c_staging.py",
        meta_table="s0c_build_meta"),
    Stage(
        key="S1", title="Astrometry solvability experiment",
        code_version="S1_CODE_VERSION",
        reads=("table:frames", "table:eras"),
        writes=("table:s1_strata", "table:s1_populations",
                "table:s1_solve_experiment", "table:s1_failure_autopsy"),
        # Three subcommands, in this order: `design` builds the strata and
        # draws the samples, `run` solves them (resumable — re-invoke until
        # it reports nothing pending), `autopsy` classifies the failures.
        # The bare script exits 2: its subparser is required=True.
        build_cmd=("python pipeline/scripts/run_s1_experiment.py design\n"
                   "python pipeline/scripts/run_s1_experiment.py run\n"
                   "python pipeline/scripts/run_s1_experiment.py autopsy"),
        meta_table="s1_build_meta"),
    Stage(
        key="S1b", title="Astrometry production batch",
        code_version="S1B_CODE_VERSION",
        reads=("table:frames", "table:eras"),
        writes=("table:s1_batch",),
        # `enqueue` ADDS newly-solvable frames without dropping the queue —
        # the right verb after an S0 rebuild, because existing solves are
        # kept.  `build --rebuild` would discard them.  The bare script
        # exits 2: its subparser is required=True.
        build_cmd=("python pipeline/scripts/run_s1_batch.py enqueue\n"
                   "python pipeline/scripts/run_s1_batch.py run"),
        meta_table="s1_batch_meta"),
    Stage(
        key="S2", title="Detector truth (ceiling, PTC, recon, linearity)",
        code_version="S2_CODE_VERSION",
        reads=("table:frames", "table:eras", "table:calib_frames",
               "table:raw_reduced_links"),
        writes=("table:s2_ceiling_modes", "table:s2_ptc_fits",
                "table:s2_recon_eras", "table:s2_linearity_ladders",
                "table:detector_params"),
        # One invocation runs ONE sub-stage; the five declared output tables
        # need five of them (plus `ceilpos`, which backfills the ceiling
        # cluster's position evidence).  Each is resumable and skips
        # finished work, so re-invoking until it reports nothing pending is
        # the intended usage.  The bare script exits 2: `stage` is a
        # required positional.
        build_cmd=("python pipeline/scripts/run_s2_campaign.py ceiling\n"
                   "python pipeline/scripts/run_s2_campaign.py ceilpos\n"
                   "python pipeline/scripts/run_s2_campaign.py ptc\n"
                   "python pipeline/scripts/run_s2_campaign.py reconstruct\n"
                   "python pipeline/scripts/run_s2_campaign.py linearity\n"
                   "python pipeline/scripts/run_s2_campaign.py params"),
        meta_table="s2_build_meta"),
    Stage(
        key="S2c", title="Per-frame dispersion (is this FILTER a spectrum?)",
        code_version="DISPERSION_CODE_VERSION",
        version_file="pipeline/rlmt_diagnostics/dispersion.py",
        version_symbol="DISPERSION_CODE_VERSION",
        reads=("table:frames", "table:eras"),
        writes=("table:frame_dispersion",),
        build_cmd=("python pipeline/scripts/run_s2c_dispersion.py build\n"
                   "python pipeline/scripts/run_s2c_dispersion.py run"),
        meta_table="s2c_build_meta",
        note="Measures source elongation per frame and returns "
             "direct/dispersed/indeterminate.  It is the evidence every "
             "REDERIVE verdict below is waiting on, so the strategy "
             "documents READ it: a verdict landing here must be able to "
             "make them stale."),
    Stage(
        key="S3", title="Timing (BJD_TDB, clock, cadence)",
        code_version="S3_CODE_VERSION",
        reads=("table:frames", "table:eras"),
        writes=("table:frame_times", "table:s3_header_audit",
                "table:s3_clock_drift", "table:s3_cadence"),
        build_cmd="python pipeline/scripts/build_s3_timing.py",
        meta_table="s3_build_meta"),
    Stage(
        key="S4", title="Ensemble photometry (AN UMa / VV Pup prototype)",
        code_version="S4_CODE_VERSION",
        # S4 reads a NARROW slice of frames (two targets), so it uses the
        # scoped resource — an S0 rebuild that repairs geometry elsewhere
        # leaves this digest identical, which is the true answer.
        #
        # It also reads the S2 detector constants, and until this edit it did
        # not say so.  The saturation veto it applies to every frame comes
        # from S2_MODE_VETO_ADU in macro_phot/series.py, copied there when the
        # S0 rebuild destroyed the s2_* tables.  With no declared edge, S4
        # reported FRESH in the same run in which S2 reported OUTPUT_MISSING:
        # a green verdict over a destroyed and admittedly-suspect input,
        # which is the dangerous direction of failure.  Declaring
        # table:detector_params and table:s2_ceiling_modes (both MISSING
        # today) makes S4 report the absent-input reason instead.
        reads=("table:frames@s4proto", "table:eras",
               "table:detector_params", "table:s2_ceiling_modes",
               "file:pipeline/macro_phot/series.py"),
        writes=("db:phot:phot_selection", "db:phot:phot_series"),
        build_cmd=("python pipeline/scripts/build_s4_photometry.py init\n"
                   "python pipeline/scripts/build_s4_photometry.py extract\n"
                   "python pipeline/scripts/build_s4_photometry.py match\n"
                   "python pipeline/scripts/build_s4_photometry.py gaia\n"
                   "python pipeline/scripts/build_s4_photometry.py ensemble\n"
                   "python pipeline/scripts/build_s4_photometry.py errors"),
        note="Opens the manifest READ-ONLY and writes its own database. "
             "`extract` and `match` are chunked (--limit, default 400) and "
             "resumable: re-invoke until `status` reports nothing pending."),
    Stage(
        key="CV-S4", title="Production CV photometry (five staged targets)",
        code_version="CV_CODE_VERSION",
        version_file="pipeline/scripts/run_cv_photometry.py",
        version_symbol="CV_CODE_VERSION",
        # The product that actually contains phantom-era frames: 207 EU UMa
        # frames in era 80 and 75 in era 78 — one camera configuration split
        # by the geometry artifact — measured as two separate series under a
        # recorded rule of "never mixed inside a series".  It reads the CV
        # staging table, the reduced links, the BJD stamps, its own frame
        # slice (for naxis1/naxis2 and the qc flags) and the same S2
        # constants as S4.
        reads=("table:stage_cv_timeseries", "table:raw_reduced_links",
               "table:frame_times", "table:frames@cv", "table:eras",
               "table:calib_frames", "table:detector_params",
               "table:s2_ceiling_modes",
               "file:pipeline/macro_phot/series.py"),
        writes=("db:cvphot:cv_selection", "db:cvphot:cv_frames"),
        build_cmd=("python pipeline/scripts/run_cv_photometry.py init\n"
                   "python pipeline/scripts/run_cv_photometry.py extract\n"
                   "python pipeline/scripts/run_cv_photometry.py match\n"
                   "python pipeline/scripts/run_cv_photometry.py field\n"
                   "python pipeline/scripts/run_cv_photometry.py ensemble\n"
                   "python pipeline/scripts/run_cv_photometry.py errors"),
        note="Built 2026-08-18T19:08Z, after the S0 rebuild and during the "
             "geometry rescue.  Undeclared until this audit, which meant "
             "`status` could exit 0 while this product rested on frames "
             "whose era assignment is known wrong."),
    Stage(
        key="CV-S5", title="CV data characterization (quality, noise, "
                           "sampling, detectability, timing)",
        code_version="CODE_VERSION",
        version_file="pipeline/scripts/run_cv_characterization.py",
        version_symbol="CODE_VERSION",
        # Reads the CV photometry product (not the archive, not the pixels
        # except for the sampled trailing audit) plus the pure-arithmetic
        # module whose constants set every cut and every threshold.
        reads=("db:cvphot:cv_selection", "db:cvphot:cv_frames",
               "file:pipeline/macro_phot/characterize.py"),
        writes=("db:cvchar:ch_noise_series", "db:cvchar:ch_frames",
                "db:cvchar:ch_cuts", "db:cvchar:ch_contour",
                "db:cvchar:ch_verdict"),
        meta_table="ch_meta",
        build_cmd=("python pipeline/scripts/run_cv_characterization.py quality\n"
                   "python pipeline/scripts/run_cv_characterization.py trail\n"
                   "python pipeline/scripts/run_cv_characterization.py noise\n"
                   "python pipeline/scripts/run_cv_characterization.py cadence\n"
                   "python pipeline/scripts/run_cv_characterization.py detect\n"
                   "python pipeline/scripts/run_cv_characterization.py timing\n"
                   "python pipeline/scripts/run_cv_characterization.py verdict"),
        note="The science verdicts on ANALYSIS_STRATEGY.md are recomputed "
             "from this product every time it is rebuilt, so a stale "
             "characterization is a stale set of verdicts."),
    Stage(
        key="CV-S6", title="CV catalogue tie (relative magnitudes -> "
                           "natural-system magnitudes on a standard zero "
                           "point)",
        code_version="CATTIE_CODE_VERSION",
        version_file="pipeline/scripts/run_cv_cattie.py",
        version_symbol="CATTIE_CODE_VERSION",
        # Reads the CV photometry product (the ensemble magnitudes it ties,
        # and the frame census that supplies the saturation veto), plus the
        # pure module whose constants set every veto, every tolerance and
        # every verdict threshold.  It also reads two external catalogues,
        # whose identity is declared as its OWN output (cv_cat_fetch) rather
        # than as an input: they are not files in this repo, and the cached
        # pull with its sha256 IS the reproducible record of what was read.
        reads=("db:cvphot:cv_selection", "db:cvphot:cv_frames",
               "file:pipeline/macro_phot/cattie.py"),
        writes=("db:cvphot:cv_cat_fetch", "db:cvphot:cv_cat_astrom",
                "db:cvphot:cv_cattie"),
        meta_table="cv_cat_meta",
        meta_version_key="cattie_code_version",
        build_cmd=("python pipeline/scripts/run_cv_cattie.py fieldfix\n"
                   "python pipeline/scripts/run_cv_cattie.py fetch\n"
                   "python pipeline/scripts/run_cv_cattie.py match\n"
                   "python pipeline/scripts/run_cv_cattie.py solve\n"
                   "python pipeline/scripts/run_cv_cattie.py validate\n"
                   "python pipeline/scripts/run_cv_cattie.py apply"),
        note="Writes ONE new column on cv_lightcurve (cal_mag) beside the "
             "untouched relative magnitude, and leaves it NULL for every "
             "block it could not tie honestly.  `fetch` caches each "
             "catalogue pull under products/phot/catalogue_cache/ with its "
             "query text and sha256 and never re-pulls silently, so a "
             "re-run is reproducible without the archive being reachable."),
    Stage(
        key="G", title="Grism extraction + identity gate (T CrB)",
        code_version="G_CODE_VERSION",
        reads=("table:frames@grism", "table:eras", "table:calib_frames"),
        writes=("table:g_extractions",),
        # EVERY action in this script is behind a flag; the bare command
        # creates two directories and exits without touching g_extractions.
        # Recording the stage after the bare command would stamp a run that
        # never happened — precisely the laundering this module exists to
        # prevent — so the declared command is --all.
        build_cmd="python pipeline/scripts/run_g_tcrb_validation.py --all",
        meta_table="g_build_meta"),
    # ---- reports: stages too.  A published page is an OUTPUT with inputs. --
    #
    # THE COMMANDS BELOW ARE THE ONES THAT RUN.  The first version of this
    # file emitted `python -c 'from macro_core.report_s0 import main; main()'`
    # for all seven, which fails twice over: macro_core is not importable from
    # the repo root (where every other command in the plan is run), and NOT
    # ONE of the seven report modules defines main() — the entry point is
    # render_report(path).  A plan is a deliverable, not a description, so
    # each line below is a command a person can paste, and
    # test_every_build_command_is_runnable executes each with --help to prove
    # it.  Three of the renderers have no CLI of their own; for those, the
    # `render` subcommand of check_pipeline_status.py calls render_report
    # directly, which keeps the fix inside this module's own files rather
    # than editing stage code another workflow is holding.
    Stage(
        key="R-S0", title="Report: S0 manifest page",
        code_version="S0_CODE_VERSION",
        reads=("table:frames", "table:eras", "table:aliases",
               "table:project_counts"),
        writes=("file:docs/pipeline/s0_manifest.html",),
        build_cmd="python pipeline/scripts/check_pipeline_status.py render R-S0"),
    Stage(
        key="R-S0b", title="Report: S0b calibration inventory page",
        code_version="S0B_CODE_VERSION",
        reads=("table:calib_coverage", "table:calib_gaps",
               "table:calib_frames", "table:eras"),
        writes=("file:docs/pipeline/s0b_calibration_inventory.html",),
        build_cmd="python pipeline/scripts/check_pipeline_status.py render R-S0b"),
    Stage(
        key="R-S0c", title="Report: S0c staging page",
        code_version="S0C_CODE_VERSION",
        reads=("table:s0c_stage_files", "table:stage_bestar_grism",
               "table:stage_cv_timeseries",
               "table:stage_dwarfgalaxy_agn_survey",
               "table:stage_sn2023ixf_lightcurve",
               "table:stage_tcrb_monitoring"),
        writes=("file:docs/pipeline/s0c_staging.html",),
        build_cmd="python pipeline/scripts/check_pipeline_status.py render R-S0c"),
    Stage(
        key="R-S1", title="Report: S1 astrometry page",
        code_version="S1_CODE_VERSION",
        reads=("table:s1_solve_experiment", "table:s1_strata",
               "table:s1_populations", "table:s1_failure_autopsy",
               "table:frames"),
        writes=("file:docs/pipeline/s1_astrometry.html",),
        build_cmd="python pipeline/scripts/run_s1_experiment.py report"),
    Stage(
        key="R-S2", title="Report: S2 detector page",
        code_version="S2_CODE_VERSION",
        reads=("table:s2_ceiling_modes", "table:s2_ptc_fits",
               "table:s2_recon_eras", "table:s2_linearity_ladders",
               "table:detector_params"),
        writes=("file:docs/pipeline/s2_detector.html",),
        build_cmd="python pipeline/scripts/run_s2_campaign.py report"),
    Stage(
        key="R-S2c", title="Report: S2c filter-identity page",
        code_version="DISPERSION_CODE_VERSION",
        version_file="pipeline/rlmt_diagnostics/dispersion.py",
        version_symbol="DISPERSION_CODE_VERSION",
        reads=("table:frame_dispersion",),
        writes=("file:docs/pipeline/s2c_filter_identity.html",),
        build_cmd="python pipeline/scripts/run_s2c_dispersion.py report"),
    Stage(
        key="R-S3", title="Report: S3 timing page",
        code_version="S3_CODE_VERSION",
        reads=("table:frame_times", "table:s3_header_audit",
               "table:s3_clock_drift", "table:s3_cadence"),
        writes=("file:docs/pipeline/s3_timing.html",),
        build_cmd="python pipeline/scripts/build_s3_timing.py --stage report"),
    Stage(
        key="R-S4", title="Report: S4 photometry page",
        code_version="S4_CODE_VERSION",
        reads=("db:phot:phot_selection", "db:phot:phot_series"),
        writes=("file:docs/pipeline/s4_photometry.html",),
        build_cmd="python pipeline/scripts/build_s4_photometry.py report"),
    Stage(
        key="R-CV-S5", title="Report: CV characterization + verdict page",
        code_version="CODE_VERSION",
        version_file="pipeline/scripts/run_cv_characterization.py",
        version_symbol="CODE_VERSION",
        reads=("db:cvchar:ch_noise_series", "db:cvchar:ch_frames",
               "db:cvchar:ch_cuts", "db:cvchar:ch_contour",
               "db:cvchar:ch_verdict",
               "file:pipeline/macro_phot/characterize.py"),
        writes=("file:docs/CV_TimeSeries/cv_characterization.html",),
        build_cmd="python pipeline/scripts/run_cv_characterization.py report"),
    Stage(
        key="R-CV-S6", title="Report: CV catalogue-tie + calibration page",
        code_version="CATTIE_CODE_VERSION",
        version_file="pipeline/scripts/run_cv_cattie.py",
        version_symbol="CATTIE_CODE_VERSION",
        reads=("db:cvphot:cv_cattie", "db:cvphot:cv_cat_fetch",
               "file:pipeline/macro_phot/cattie.py"),
        writes=("file:docs/CV_TimeSeries/cv_catalogue_tie.html",),
        build_cmd="python pipeline/scripts/run_cv_cattie.py report"),
    # ---- human-authored artifacts, declared so they cannot hide ------------
    Stage(
        key="OPS", title="Observatory request (October shopping list)",
        code_version="(hand-authored)",
        reads=("table:calib_gaps", "table:calib_coverage", "table:eras"),
        writes=("file:ops/2026-08_observatory_request.md",),
        build_cmd="(hand-authored — regenerate the numbers from calib_gaps)",
        note="Not machine-generated: its numbers were transcribed from "
             "calib_gaps.  Declared here precisely BECAUSE it is hand-typed; "
             "a hand-typed number is the kind that silently goes stale."),
    Stage(
        key="STRAT", title="Five ANALYSIS_STRATEGY.md documents + ROADMAP",
        code_version="(hand-authored)",
        # table:frame_dispersion is the edge that matters most here.  Each of
        # these documents contains at least one rule of the form "FILTER x is
        # a spectrum, so exclude it" — and S2c is measuring, frame by frame,
        # whether that is true.  Without this edge the classifier could
        # finish, overturn an exclusion list, and nothing in the DAG would go
        # stale: the same silent drift, one level up.
        reads=("table:project_counts", "table:s0c_stage_files",
               "table:frames", "table:frame_dispersion"),
        writes=tuple(f"file:{p}/ANALYSIS_STRATEGY.md" for p in
                     ("BeStar_Grism", "CV_TimeSeries",
                      "DwarfGalaxy_AGN_Survey", "SN2023ixf_LightCurve",
                      "TCrB_Monitoring")) + ("file:ROADMAP.md",),
        build_cmd=("(hand-authored — reconcile against project_counts and "
                   "the per-frame verdicts in frame_dispersion)"),
        note="Hand-authored; their per-target counts and filter rules are "
             "claims about frames and must be re-reconciled when it moves."),
    Stage(
        key="WEB", title="Public project pages (docs/*/index.html)",
        code_version="(hand-authored)",
        reads=("table:project_counts", "table:s0c_stage_files"),
        writes=tuple(f"file:docs/{p}/index.html" for p in
                     ("BeStar_Grism", "CV_TimeSeries",
                      "DwarfGalaxy_AGN_Survey", "Legacy_Rigel",
                      "SN2023ixf_LightCurve", "TCrB_Monitoring"))
               + ("file:docs/index.html",),
        build_cmd="(hand-authored — reconcile against project_counts)",
        note="The pages a reader outside this repo actually sees.  The "
             "status page used to print a verdict for them while the exit "
             "code ignored them entirely; declaring them makes the gate "
             "cover what the page judges."),
)

STAGE_BY_KEY: dict[str, Stage] = {s.key: s for s in STAGES}


def producer_of(resource: str) -> Optional[str]:
    """Which stage writes ``resource``?  None when nothing declares it.

    Used to turn the reads/writes lists into DAG edges without a second,
    hand-maintained edge list that could disagree with them.
    """
    for stage in STAGES:
        if resource in stage.writes:
            return stage.key
    return None


# ===========================================================================
# 3.  FINGERPRINTS
# ===========================================================================

@dataclass(frozen=True)
class Fingerprint:
    """``(n_rows, digest)`` for a present resource; ``present=False`` for a
    resource that does not exist (a dropped table, a deleted file).

    "Missing" is a first-class state, not ``n_rows=0``: an empty table and
    an absent table mean completely different things to an auditor, and the
    S0 table swap produced the second kind."""

    present: bool
    n_rows: Optional[int] = None
    digest: Optional[str] = None

    @property
    def token(self) -> str:
        """Compact string stored in ``stage_provenance`` (and compared)."""
        if not self.present:
            return "MISSING"
        return f"{self.n_rows}:{self.digest}"

    @classmethod
    def from_token(cls, token: str) -> "Fingerprint":
        """Inverse of :attr:`token`.  Raises on anything else — a record we
        cannot parse must not be read as agreement."""
        if token == "MISSING":
            return cls(present=False)
        if ":" not in token:
            raise ProvenanceError(f"malformed fingerprint token: {token!r}")
        n, _, d = token.partition(":")
        try:
            return cls(present=True, n_rows=int(n), digest=d)
        except ValueError as exc:                     # not an integer count
            raise ProvenanceError(
                f"malformed fingerprint token: {token!r}") from exc


def _canon(value) -> str:
    """One cell -> one canonical string.

    NULL is a distinct symbol (``\\x00``) rather than the empty string, so a
    NULL filter and an empty-string filter cannot collide.  Floats are
    rendered with ``repr`` after the SQL layer has already rounded them, so
    the text is stable across platforms.
    """
    if value is None:
        return "\x00"
    if isinstance(value, float):
        # -0.0 and 0.0 must hash identically; adding 0.0 normalizes the sign.
        return repr(value + 0.0)
    if isinstance(value, bytes):
        return value.hex()
    return str(value)


def digest_rows(rows: Iterable[Sequence]) -> tuple[int, str]:
    """Hash an ordered row stream -> ``(n_rows, sha256-hex[:16])``.

    Cells are joined with ``\\x1f`` (unit separator) and rows with ``\\x1e``
    (record separator) — two control characters that cannot appear in the
    values being hashed, so no ambiguity can be constructed by a value that
    happens to contain a delimiter.  16 hex characters = 64 bits: collision
    probability across the ~40 resources here is nil, and a short digest is
    readable in a terminal, which matters for a tool people must trust.
    """
    h = hashlib.sha256()
    n = 0
    for row in rows:
        h.update("\x1f".join(_canon(v) for v in row).encode("utf-8"))
        h.update(b"\x1e")
        n += 1
    return n, h.hexdigest()[:16]


def digest_bytes(data: bytes) -> tuple[int, str]:
    """Hash a file's bytes -> ``(n_bytes, sha256-hex[:16])``.

    For a published HTML page or a hand-authored Markdown file there are no
    "columns that matter" — every byte is the artifact — so the whole file
    is hashed and ``n_rows`` carries the byte count.
    """
    return len(data), hashlib.sha256(data).hexdigest()[:16]


def _where_clause(spec: ResourceSpec) -> str:
    """`` WHERE <predicate>`` for a scoped spec, or ``''``.

    Pure string assembly, kept in one place so that every kind of resource
    that supports scoping builds the clause the same way — a scoped table and
    a scoped sibling-database table must not diverge in how they filter, or
    two fingerprints of "the same" slice would disagree.
    """
    return f" WHERE {spec.where}" if spec.where else ""


def read_version_constant(source: str, symbol: str) -> Optional[str]:
    """Extract a module-level string constant from PYTHON SOURCE TEXT.

    Used for the two stages whose code-version constant lives in a *script*
    rather than an importable package (``CV_CODE_VERSION`` in
    run_cv_photometry.py, ``DISPERSION_CODE_VERSION`` in dispersion.py).

    It parses rather than imports, and that choice is deliberate on two
    counts.  Importing a script executes it — module-level ``sys.path``
    surgery, heavy astronomy imports, and, while sibling workflows are
    editing those very files, whatever half-saved state they are in.  A
    status check must be able to run at any moment without side effects.
    Parsing gives the same answer for the cost of reading a file.

    Returns ``None`` when the symbol is absent or is not a plain string
    literal — the caller then reports "version unknown" rather than
    inventing agreement.
    """
    import ast
    try:
        tree = ast.parse(source)
    except SyntaxError:
        # A file being written to right now can parse as garbage.  That is
        # "unknown", not "unchanged".
        return None
    for node in tree.body:                      # module level only: a
        if not isinstance(node, ast.Assign):    # constant defined inside a
            continue                            # function is not the API
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == symbol:
                value = node.value
                if isinstance(value, ast.Constant) and \
                        isinstance(value.value, str):
                    return value.value
                return None
    return None


def _table_exists(con: sqlite3.Connection, name: str) -> bool:
    """True when a table or view of this name exists in the connection."""
    row = con.execute(
        "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') "
        "AND name = ?", (name,)).fetchone()
    return row is not None


def fingerprint_resource(spec: ResourceSpec, con: sqlite3.Connection,
                         repo_root) -> Fingerprint:
    """Compute one resource's fingerprint.  Thin I/O around the pure hashers.

    ``con`` is the manifest connection; ``repo_root`` resolves file and
    sibling-database paths.  A resource that does not exist returns a
    ``present=False`` fingerprint — never an invented zero.
    """
    if spec.kind == "table":
        if not _table_exists(con, spec.name):
            return Fingerprint(present=False)
        sql = (f"SELECT {', '.join(spec.columns)} FROM {spec.name}"
               f"{_where_clause(spec)} ORDER BY {spec.order_by}")
        n, d = digest_rows(con.execute(sql))
        return Fingerprint(present=True, n_rows=n, digest=d)

    if spec.kind == "file":
        path = os.path.join(str(repo_root), spec.name)
        if not os.path.exists(path):
            return Fingerprint(present=False)
        with open(path, "rb") as fh:
            n, d = digest_bytes(fh.read())
        return Fingerprint(present=True, n_rows=n, digest=d)

    if spec.kind == "stat":
        # Size + modification time, NOT a content hash.  Used for the one
        # input that is neither ours nor stable: the 400 MB external scan
        # catalog on the ASTRO drive, which sibling processes rewrite while
        # this check runs.  Hashing it would be slow AND non-reproducible.
        # This is an integrity SURROGATE and is labelled as one wherever it
        # is printed: it catches "the catalog was rewritten", which is the
        # only question S0's freshness turns on.  It cannot catch an
        # in-place edit that preserves both size and mtime.
        path = spec.name if os.path.isabs(spec.name) \
            else os.path.join(str(repo_root), spec.name)
        if not os.path.exists(path):
            return Fingerprint(present=False)
        st = os.stat(path)
        n, d = digest_rows([(int(st.st_size), int(st.st_mtime))])
        return Fingerprint(present=True, n_rows=int(st.st_size), digest=d)

    if spec.kind == "db":
        path = os.path.join(str(repo_root), spec.database or "")
        if not os.path.exists(path):
            return Fingerprint(present=False)
        # Sibling databases are opened READ-ONLY: a status check must never
        # be able to modify a product it is only inspecting.
        side = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
        try:
            side.execute("PRAGMA busy_timeout = 300000")
            if not _table_exists(side, spec.name):
                return Fingerprint(present=False)
            if spec.columns == ("count",):
                # Count-only resources (see phot_series' why): the row count
                # IS the fingerprint, and the digest records that choice so
                # a later reader cannot mistake it for a content hash.
                n = side.execute(f"SELECT count(*) FROM {spec.name}"
                                 f"{_where_clause(spec)}").fetchone()[0]
                return Fingerprint(present=True, n_rows=int(n),
                                   digest="count-only")
            cols = [c for c in spec.columns if c != "*none*"]
            select = ", ".join(cols) if cols else "*"
            n, d = digest_rows(side.execute(
                f"SELECT {select} FROM {spec.name}{_where_clause(spec)} "
                f"ORDER BY {spec.order_by}"))
            return Fingerprint(present=True, n_rows=n, digest=d)
        finally:
            side.close()

    raise ProvenanceError(f"unknown resource kind: {spec.kind!r}")


def fingerprint_all(keys: Iterable[str], con: sqlite3.Connection,
                    repo_root) -> dict[str, str]:
    """Fingerprint tokens for a set of resource keys.

    An undeclared key raises: the DAG and the registry must agree, and a
    silent skip here would make a stage look fresh because one of its
    inputs was never checked.
    """
    out: dict[str, str] = {}
    for key in keys:
        spec = RESOURCES.get(key)
        if spec is None:
            raise ProvenanceError(f"resource {key!r} has no spec — declare "
                                  f"it in RESOURCES before using it")
        out[key] = fingerprint_resource(spec, con, repo_root).token
    return out


# ===========================================================================
# 4.  STALENESS  (pure)
# ===========================================================================

FRESH = "FRESH"
STALE = "STALE"
#: A stage whose OWN declared inputs still match what it recorded, but which
#: sits downstream of a stage that is not fresh.  It must still re-run — its
#: inputs are about to be rebuilt — but it is a different kind of "must", and
#: keeping the two apart is what makes the plan triageable.  A STALE stage
#: has evidence against it; a STALE_UPSTREAM stage is merely waiting.
STALE_UPSTREAM = "STALE_UPSTREAM"
NEVER_RUN = "NEVER_RUN"
OUTPUT_MISSING = "OUTPUT_MISSING"

#: Every verdict, in decreasing severity.  Iterating this rather than a
#: hand-written list in each caller means a new state cannot be silently
#: dropped from a summary line or a colour key.
ALL_STATES = (OUTPUT_MISSING, NEVER_RUN, STALE, STALE_UPSTREAM, FRESH)


@dataclass(frozen=True)
class Record:
    """One row of ``stage_provenance``: what a stage saw and produced."""

    stage: str
    run_utc: str
    code_version: str
    git_commit: str
    inputs: Mapping[str, str]
    outputs: Mapping[str, str]
    note: str = ""


@dataclass(frozen=True)
class Freshness:
    """The verdict for one stage, with the reasons spelled out.

    ``reasons`` is what makes this usable: "STALE" alone sends a person
    back to reasoning from memory, which is the failure mode being fixed.
    """

    stage: str
    state: str
    reasons: tuple[str, ...] = ()
    changed_inputs: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.state == FRESH


def compare_fingerprints(recorded: Mapping[str, str],
                         current: Mapping[str, str]) -> list[tuple[str, str]]:
    """Diff two fingerprint maps -> ``[(resource, human reason), ...]``.

    Four distinguishable outcomes, because they call for different actions:
    an input that VANISHED, an input that APPEARED, an input whose CONTENT
    changed, and an input that was never recorded at all (an older record
    written before that dependency was declared).
    """
    changes: list[tuple[str, str]] = []
    for key in sorted(set(recorded) | set(current)):
        was = recorded.get(key)
        now = current.get(key)
        if was is not None and was.startswith(UNRECORDED):
            # Backfill sentinel: we know this input moved after the stage ran
            # (its producer was rebuilt later), but not what it looked like at
            # the time.  Reporting it as "changed" is the honest reading —
            # the alternative, assuming it matched, is the lie that produced
            # this whole situation.
            why = was[len(UNRECORDED):].lstrip(":") or "state at last run unknown"
            changes.append((key, f"changed after this stage ran — {why}"))
        elif was is None:
            changes.append((key, "dependency was not recorded at last run"))
        elif now is None:
            changes.append((key, "dependency is no longer declared"))
        elif was == now:
            continue
        elif was == "MISSING":
            changes.append((key, f"appeared (now {now})"))
        elif now == "MISSING":
            changes.append((key, f"DESTROYED (was {was})"))
        else:
            was_n = was.split(":")[0]
            now_n = now.split(":")[0]
            if was_n != now_n:
                changes.append((key, f"rows {was_n} -> {now_n}"))
            else:
                changes.append((key, f"content changed ({was_n} rows, "
                                     f"digest {was.split(':')[1]} -> "
                                     f"{now.split(':')[1]})"))
    return changes


def is_stale(stage: Stage, record: Optional[Record],
             current_inputs: Mapping[str, str],
             current_outputs: Mapping[str, str],
             current_code_version: str) -> Freshness:
    """The core question, answered from data alone.

    Order of tests matters and is deliberate:

    1. no record at all  -> NEVER_RUN (nothing to compare; do not guess);
    2. an output missing -> OUTPUT_MISSING (the strongest verdict: whatever
       the inputs say, the result is gone — this is what the S0 table swap
       did to S1 and S2, and it must not be softened to "stale");
    3. code version moved -> STALE (same bytes, different rules);
    4. any input changed  -> STALE, listing which and how;
    5. otherwise          -> FRESH.
    """
    absent_out = sorted(k for k, v in current_outputs.items()
                        if v == "MISSING")
    if record is None:
        why = ["no stage_provenance row: this stage has never been recorded, "
               "so nothing about it is verified"]
        if absent_out:
            why.append("and its declared outputs are ABSENT from the "
                       "database: " + ", ".join(absent_out))
        return Freshness(stage.key, NEVER_RUN, tuple(why),
                         tuple(absent_out))

    reasons: list[str] = []

    # An input that no longer EXISTS is its own category of wrong: the
    # evidence this stage was built from is gone, so the stage cannot even
    # be re-run until its producer is rebuilt.  Checked before the diff
    # because "MISSING == MISSING" would otherwise read as agreement.
    for key in sorted(current_inputs):
        if current_inputs[key] == "MISSING":
            reasons.append(f"input {key} is ABSENT — the evidence this "
                           f"stage was built from no longer exists")

    gone = [k for k, v in current_outputs.items() if v == "MISSING"]
    if gone:
        reasons.append("output(s) absent: " + ", ".join(sorted(gone)))
        return Freshness(stage.key, OUTPUT_MISSING, tuple(reasons),
                         tuple(sorted(gone)))

    if record.code_version and current_code_version and \
            record.code_version != current_code_version:
        reasons.append(f"code version {record.code_version!r} -> "
                       f"{current_code_version!r}")

    changes = compare_fingerprints(record.inputs, current_inputs)
    for key, why in changes:
        reasons.append(f"input {key}: {why}")

    # An output that changed since the record without the stage being re-run
    # means somebody wrote to it out of band.  Worth saying out loud.
    for key, token in current_outputs.items():
        was = record.outputs.get(key)
        if was is not None and was != token:
            reasons.append(f"output {key} changed since it was recorded "
                           f"({was} -> {token}) — written out of band")

    if reasons:
        return Freshness(stage.key, STALE, tuple(reasons),
                         tuple(k for k, _ in changes))
    return Freshness(stage.key, FRESH)


def topological_order(stages: Sequence[Stage]) -> list[str]:
    """Stage keys in dependency order (producers before consumers).

    Kahn's algorithm over edges derived from reads/writes.  A cycle raises
    rather than returning a partial order: an order that silently omits
    stages would produce a re-run plan that misses work.
    """
    keys = [s.key for s in stages]
    writer = {}
    for s in stages:
        for w in s.writes:
            writer[w] = s.key
    # edge producer -> consumer
    incoming: dict[str, set[str]] = {k: set() for k in keys}
    outgoing: dict[str, set[str]] = {k: set() for k in keys}
    for s in stages:
        for r in s.reads:
            src = writer.get(r)
            if src is not None and src != s.key:
                incoming[s.key].add(src)
                outgoing[src].add(s.key)
    order: list[str] = []
    ready = sorted(k for k in keys if not incoming[k])
    pending = {k: set(v) for k, v in incoming.items()}
    while ready:
        node = ready.pop(0)
        order.append(node)
        for nxt in sorted(outgoing[node]):
            pending[nxt].discard(node)
            if not pending[nxt] and nxt not in order and nxt not in ready:
                ready.append(nxt)
        ready.sort()
    if len(order) != len(keys):
        missing = sorted(set(keys) - set(order))
        raise ProvenanceError(f"dependency cycle involving: {missing}")
    return order


def propagate_staleness(freshness: Mapping[str, Freshness],
                        stages: Sequence[Stage]) -> dict[str, Freshness]:
    """Mark a FRESH stage STALE_UPSTREAM when an ancestor is not fresh.

    A stage can be byte-for-byte consistent with its recorded inputs and
    still be wrong, because those inputs are themselves about to be
    rebuilt.  Without this pass the tool would tell a reader that S0c is
    fine while S0 — which S0c reads — is queued for a rebuild.

    The verdict is STALE_UPSTREAM rather than STALE, and the difference is
    the review finding that produced it.  Told only "STALE", a reader
    cannot tell the stage whose own evidence changed from the stage that is
    simply downstream of one, so the plan reads as "re-run everything" and
    the cheapest true correction looks as expensive as the worst.  The
    reason line here names the ancestors AND says the stage's own declared
    inputs currently match — which, combined with the scoped resources
    above, is often enough to conclude that after the rebuild it will have
    nothing to do.
    """
    order = topological_order(stages)
    by_key = {s.key: s for s in stages}
    writer = {w: s.key for s in stages for w in s.writes}
    out = dict(freshness)
    for key in order:
        stage = by_key[key]
        bad_parents = []
        for r in stage.reads:
            src = writer.get(r)
            if src and src != key and src in out and not out[src].ok:
                bad_parents.append(f"{src} ({out[src].state})")
        if bad_parents and out.get(key) and out[key].ok:
            out[key] = Freshness(
                key, STALE_UPSTREAM,
                (f"upstream not fresh: {', '.join(sorted(set(bad_parents)))}"
                 f" — this stage's own declared inputs still match what it "
                 f"recorded, so re-run it only after those ancestors, and "
                 f"only if its own fingerprints move",))
    return out


def rerun_plan(freshness: Mapping[str, Freshness],
               stages: Sequence[Stage]) -> list[str]:
    """The stage keys that must re-run, in dependency order.

    Everything that is not FRESH, ordered topologically so that a person
    executing the list top to bottom never rebuilds a stage before the
    stage it reads.
    """
    order = topological_order(stages)
    return [k for k in order
            if k in freshness and not freshness[k].ok]


# ===========================================================================
# 5.  PERSISTENCE  (thin I/O)
# ===========================================================================

DDL = f"""
CREATE TABLE IF NOT EXISTS {PROVENANCE_TABLE} (
    stage         TEXT NOT NULL,   -- Stage.key
    run_utc       TEXT NOT NULL,   -- when the stage produced these outputs
    code_version  TEXT,            -- the stage's own version constant
    git_commit    TEXT,            -- repo state at that run
    prov_version  TEXT,            -- PROVENANCE_CODE_VERSION (digest rules)
    inputs_json   TEXT NOT NULL,   -- {{resource: fingerprint token}} AS SEEN
    outputs_json  TEXT NOT NULL,   -- {{resource: fingerprint token}} PRODUCED
    note          TEXT,
    PRIMARY KEY (stage, run_utc)
)
"""


def ensure_table(con: sqlite3.Connection) -> None:
    """Create ``stage_provenance`` if absent.  Purely additive: this module
    never drops or rewrites another stage's table."""
    con.execute(DDL)
    con.commit()


def record_run(con: sqlite3.Connection, stage: str, run_utc: str,
               code_version: str, git_commit: str,
               inputs: Mapping[str, str], outputs: Mapping[str, str],
               note: str = "") -> None:
    """Write one provenance row (idempotent on (stage, run_utc))."""
    ensure_table(con)
    con.execute(
        f"INSERT OR REPLACE INTO {PROVENANCE_TABLE} "
        "(stage, run_utc, code_version, git_commit, prov_version, "
        " inputs_json, outputs_json, note) VALUES (?,?,?,?,?,?,?,?)",
        (stage, run_utc, code_version, git_commit, PROVENANCE_CODE_VERSION,
         json.dumps(dict(inputs), sort_keys=True),
         json.dumps(dict(outputs), sort_keys=True), note))
    con.commit()


def recorded_run_times(con: sqlite3.Connection, stage: str) -> list[str]:
    """Every ``run_utc`` already recorded for one stage, oldest first.

    ``record`` uses it to refuse two things that used to pass silently: a
    second row on a ``run_utc`` that already exists (``INSERT OR REPLACE``
    would overwrite the earlier run, losing exactly the history this table
    exists to keep) and a "new" run whose timestamp does not advance (which
    is what a re-run whose build_meta was not updated looks like — a stamp
    that precedes the work it claims to attest).
    """
    if not _table_exists(con, PROVENANCE_TABLE):
        return []
    return [r[0] for r in con.execute(
        f"SELECT run_utc FROM {PROVENANCE_TABLE} WHERE stage = ? "
        f"ORDER BY run_utc", (stage,))]


def read_records(con: sqlite3.Connection) -> dict[str, Record]:
    """Latest provenance record per stage (by run_utc), or {} if the table
    does not exist yet."""
    if not _table_exists(con, PROVENANCE_TABLE):
        return {}
    out: dict[str, Record] = {}
    for row in con.execute(
            f"SELECT stage, run_utc, code_version, git_commit, inputs_json, "
            f"outputs_json, coalesce(note,'') FROM {PROVENANCE_TABLE} "
            f"ORDER BY stage, run_utc"):
        stage, run_utc, cv, gc, ij, oj, note = row
        # ORDER BY run_utc means the last row seen per stage is the newest.
        out[stage] = Record(stage=stage, run_utc=run_utc,
                            code_version=cv or "", git_commit=gc or "",
                            inputs=json.loads(ij), outputs=json.loads(oj),
                            note=note)
    return out
