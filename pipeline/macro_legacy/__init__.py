"""macro_legacy — census of the pre-MACRO (2015–2022) Winer/Iowa archive.

Sibling package to ``macro_core`` / ``macro_sn`` / ``macro_grism``, built to
the same standards: pure, unit-tested decision logic in logic modules; one
resumable CLI per stage under ``pipeline/scripts``; a report renderer that
derives EVERY published number from the census database.

SCOPE — A CENSUS AND NOTHING MORE (committee ruling U8, 2026-10-03)
--------------------------------------------------------------------
The legacy archive is a candidate, not a project.  This package answers
"what is on disk?" (L0), "what could it support?" (L1) and applies a
go/no-go rule that was written down BEFORE any census number existed (L2,
``Legacy_Rigel/notes/GO_NOGO_PREREGISTERED.md``).  It opens no pixel: every
statement here is a statement about headers and file names, and the report
says so wherever that limits what can be concluded (saturation, absolute
clock and gain are all pixel questions the census cannot settle).

* ``macro_legacy.scan``    — L0 header scan: which cards are read, how one
  header becomes one row, and the tolerant card merge.  The only module that
  touches FITS files.
* ``macro_legacy.census``  — L0/L1/L2 pure logic: reconciliation against the
  transfer manifests, camera identity, frame kinds, night labels, the dedup
  rule, target aliases, runs / seasons, calibration availability, the
  time-convention audit and the pre-registered gates.  No file or database
  access; unit tests in ``pipeline/tests/test_legacy_census.py``.
* ``macro_legacy.clock``   — the per-season clock audit's pure logic
  (TESS reference ephemerides, prediction ranges, O − C under both stamp
  readings, season verdicts).  Build: ``build_legacy_clock.py``.
* ``macro_legacy.report``  — the evidence page
  ``docs/Legacy_Rigel/legacy_census.html`` in the house Socratic format
  (Question → Evidence → Decision → Consequence).

Build entry points: ``pipeline/scripts/build_legacy_scan.py`` (the slow,
resumable header scan) and ``pipeline/scripts/build_legacy_census.py`` (the
fast derived tables, the gates and the page).

THE ARCHIVE IS READ-ONLY
------------------------
Nothing in this package writes under
``/Volumes/OWC StudioStack HDD/DATA/ASTRO``.  The census database lives in
``products/legacy/legacy.sqlite``; the RLMT manifest is opened ``mode=ro``.
"""

#: Recorded into ``scan_meta`` / ``census_meta`` so a reader of the database
#: can tell which rules produced it.  Bump when the logic changes content.
LEGACY_SCAN_VERSION = "L0-scan v1.0 (2026-10-03)"
# v1.1: truncated-at-source dispositions, filename parser, camera + focal
# length era keys with the RLMT sharing table, cross-archive copies
# excluded, per-season clock audit (macro_legacy.clock).
LEGACY_CENSUS_VERSION = "L-census v1.1 (2026-10-04)"
