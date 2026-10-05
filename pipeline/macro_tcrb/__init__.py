"""macro_tcrb — the T CrB Pre-Eruption Monitoring project's own stages.

Sibling package to ``macro_sn`` / ``macro_phot`` / ``macro_grism``, built to
the same standard: pure, unit-tested decision logic in logic modules; one
resumable CLI (``pipeline/scripts/run_tcrb.py``); every published number
emitted from a database (``products/tcrb/tcrb.sqlite`` plus read-only reads
of the manifest and the grism library) into
``manuscripts/TCrB_Monitoring/numbers.tex``.

Modules
-------
* ``macro_tcrb.db``      — the project database: connection, metadata, the
  one table writer every stage uses.
* ``macro_tcrb.detect``  — the PRE-REGISTERED EW-change detection rule
  (TCRB-A5a; ANALYSIS_STRATEGY.md §10, fixed 2026-10-05T02:40Z) as code.
* ``macro_tcrb.p0``      — Phase-0 gates: mechanical epochs, the
  temperature split and flanking-band adequacy, filter forensics, ZMAG
  provenance, short-exposure timing.
* ``macro_tcrb.phot``    — Phase B (B-only anchors) and the archival
  flickering limits (C1).
* ``macro_tcrb.spec``    — EW / line flux / cross-validation (A5, A5b, A7, A8)
  on top of the shared grism library (G-1…G-5).
* ``macro_tcrb.paper``   — figures and ``numbers.tex``.

Unit tests: ``pipeline/tests/test_tcrb.py``.
"""

#: Recorded in ``tcrb_meta`` and stamped on every table.  Bump when a rule
#: changes in a way that could alter a stored number.
TCRB_CODE_VERSION = "TCRB v1.0 (2026-10-05)"
