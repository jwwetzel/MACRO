"""macro_core — shared MACRO Consortium pipeline core.

Subpackage layout (ROADMAP.md section 1):

* ``macro_core.manifest``  — S0 pure logic: target-name normalization, global
  deduplication, era keying, night labels, pointing math.  Every function here
  is side-effect free and unit-tested in ``pipeline/tests/test_manifest.py``.
* ``macro_core.report_s0`` — S0 chain-of-evidence report renderer: reads the
  manifest database and emits ``docs/pipeline/s0_manifest.html`` plus every
  figure under ``docs/pipeline/figures/s0/``.  No number in that report is
  hand-typed; each one is the result of a SQL query executed by this module.

The S0 build entry point lives in ``pipeline/scripts/build_s0_manifest.py``.
"""

# The version note recorded into the manifest's build_meta table.  Bump the
# string whenever the S0 logic changes in a way that alters manifest content,
# so downstream stages can tell which rules produced the file they read.
S0_CODE_VERSION = "S0 v1.0 (2026-08-17)"
