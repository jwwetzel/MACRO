"""macro_dw — the Dwarf-Galaxy Hα vetting project's own pipeline stages.

Sibling package to ``macro_core`` / ``macro_phot`` / ``macro_sn``, built to
the same standards: pure, unit-tested decision rules in one module, pixel
I/O in another, one resumable CLI under ``pipeline/scripts``, and a paper
module that emits every published number from the project database.

* ``macro_dw.dwcore``   — pure rules: frame disposition, QC gates, the
  superflat validation metric, Román et al. (2020) depth, Sérsic algebra,
  the NGC 5238 line-flux calibration arithmetic, the multiband
  periodogram with simultaneous covariate fit and its injection grid.
  Unit tests: ``pipeline/tests/test_dw.py``.
* ``macro_dw.dwio``     — frame I/O: dark/flat calibration, source masks,
  plate solutions refined on REFCAT2, multi-aperture forced photometry,
  resampling onto a sky grid, weighted stacks.
* ``macro_dw.paper_dw`` — the manuscript's figures and ``numbers.tex``.

Driver: ``pipeline/scripts/run_dw_paper.py``.  Database:
``products/dwarf/dwarf.sqlite`` (the manifest is opened read-only).
"""

#: Version of the DW stage (provenance graph key "DW").
DW_CODE_VERSION = "DW v1.0 (2026-10-05) — Hα paper pipeline"
