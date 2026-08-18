"""macro_phot — S4 ensemble-photometry core for the MACRO/RLMT archive.

Sibling package to ``macro_core`` (S0/S0b manifest + inventory).  Same house
rules: pure unit-testable functions in the logic modules, thin I/O wrappers
kept separate, chain-of-evidence report rendered FROM the database, atomic
output writes, no hand-typed number anywhere in the HTML.

Module layout:

* ``macro_phot.photometry`` — pure aperture/plate-scale/matching logic:
  plate-scale arithmetic, aperture scaling per era, FWHM estimation,
  reference-frame selection, one-to-one nearest-neighbour star matching,
  instrumental magnitudes and their photon-noise errors.  Unit-tested in
  ``pipeline/tests/test_phot.py``.
* ``macro_phot.ensemble``   — pure Honeycutt (1992, PASP 104, 435)
  inhomogeneous-ensemble solver: robust alternating least squares for
  per-frame zero points + per-star mean magnitudes, comparison-star
  selection by stability iteration, fixed-ZP evaluation of held-out stars.
  The REQUIRED synthetic-recovery test (injected ZP pattern) lives in
  ``test_phot.py``.
* ``macro_phot.errors``     — pure empirical-error-model arithmetic (the S5
  seed): RMS-vs-magnitude binning, chi-square of constant-star fits and the
  error inflation factor, Allan deviation of a time series.
* ``macro_phot.extract``    — the thin I/O layer: read one reduced frame,
  run sep background/detection/aperture photometry with parameters computed
  by the pure layer.
* ``macro_phot.gaia``       — Gaia DR3 cone query (I/O) + pure tangent-plane
  projection and parity-tolerant field identification used to tie the
  ensemble to the sky.
* ``macro_phot.report_s4``  — S4 chain-of-evidence report renderer: reads
  the photometry database and emits ``docs/pipeline/s4_photometry.html``
  plus every figure under ``docs/pipeline/figures/s4/``.

Build entry point: ``pipeline/scripts/build_s4_photometry.py`` (resumable
staged CLI — the extraction loop is chunked so no single invocation runs
long).
"""

# Version note recorded into the photometry DB's s4_build_meta table.  Bump
# whenever the S4 logic changes in a way that alters database content.
S4_CODE_VERSION = "S4 v1.0 (2026-08-17)"
