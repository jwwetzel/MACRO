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

The build entry points live in ``pipeline/scripts/build_s0_manifest.py`` and
``pipeline/scripts/build_s0b_inventory.py``.
"""

# The version note recorded into the manifest's build_meta table.  Bump the
# string whenever the S0 logic changes in a way that alters manifest content,
# so downstream stages can tell which rules produced the file they read.
S0_CODE_VERSION = "S0 v1.0 (2026-08-17)"

# Same contract for the S0b inventory tables (recorded in s0b_build_meta).
# v1.1: header-glitch FILTER strings (calibration-vocabulary collisions)
# excluded from the shopping list; re-opening-configuration eras surfaced.
S0B_CODE_VERSION = "S0b v1.1 (2026-08-17)"
