"""macro_core — shared MACRO Consortium pipeline core.

Subpackage layout (ROADMAP.md section 1):

* ``macro_core.manifest``  — S0 pure logic: target-name normalization, global
  deduplication, era keying, night labels, pointing math.  Every function here
  is side-effect free and unit-tested in ``pipeline/tests/test_manifest.py``.
* ``macro_core.report_s0`` — S0 chain-of-evidence report renderer: reads the
  manifest database and emits ``docs/pipeline/s0_manifest.html`` plus every
  figure under ``docs/pipeline/figures/s0/``.  No number in that report is
  hand-typed; each one is the result of a SQL query executed by this module.
* ``macro_core.inventory`` — S0b pure logic: raw<->reduced match ladder,
  calibration-kind normalization, dark exposure matching, coverage/gap
  arithmetic.  Unit-tested in ``pipeline/tests/test_inventory.py``.
* ``macro_core.report_s0b`` — S0b evidence report renderer: reads the S0b
  tables and emits ``docs/pipeline/s0b_calibration_inventory.html`` plus the
  figures under ``docs/pipeline/figures/s0b/``.  Same rule: every number is
  a query result.

* ``macro_core.staging``   — S0c pure logic: the five per-project staging
  selections (reviewable data), science/calibration row builders, the
  no-copy law's helpers.  Unit-tested in ``pipeline/tests/test_staging.py``.
* ``macro_core.report_s0c`` — S0c evidence report renderer: reads the stage
  tables and emits ``docs/pipeline/s0c_staging.html`` plus the figures under
  ``docs/pipeline/figures/s0c/``.  Same rule: every number is a query result.

The build entry points live in ``pipeline/scripts/build_s0_manifest.py``,
``pipeline/scripts/build_s0b_inventory.py``, and
``pipeline/scripts/build_s0c_staging.py``.
"""

# The version note recorded into the manifest's build_meta table.  Bump the
# string whenever the S0 logic changes in a way that alters manifest content,
# so downstream stages can tell which rules produced the file they read.
# v1.1 (plan review 2026-10-03): dedup on the EXPOSURE ROOT so reduced and
# re-packaged copies fold into their raw parent (F-1); typed hardware-state
# columns from the header re-scrape (F-2); era alias columns; the S0 tables
# are swapped into the live manifest in one transaction, not by file replace.
S0_CODE_VERSION = "S0 v1.1 (2026-10-03)"

# Same contract for the S0b inventory tables (recorded in s0b_build_meta).
# v1.1: header-glitch FILTER strings (calibration-vocabulary collisions)
# excluded from the shopping list; re-opening-configuration eras surfaced.
# v1.2 (plan review 2026-10-03): the mechanical-epoch layer beneath the eras
# (mech_epoch, night_mech_epoch, frame_mech_epoch) with its null and
# injection tests; coverage re-counted per mechanical epoch under the
# boundary rule and the gain/offset/set-point key (F-3, DE.F5).
S0B_CODE_VERSION = "S0b v1.2 (2026-10-03)"

# Same contract for the S0c staging manifests (recorded in s0c_build_meta and
# in every stage row's stage_build_id).
# v1.1 (plan review 2026-10-03): calibration rows are staged only where the
# mechanical-epoch rule lets them serve the project's science; every stage
# row carries mech_epoch / detector_epoch / epoch_certain (F-3).
S0C_CODE_VERSION = "S0c v1.1 (2026-10-03)"
