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

The PRODUCTION layer, added when the proven core was run over the whole
staged CV set (five targets, 14 (target, era) blocks, 8,716 frames) rather
than two prototype polars.  These modules hold the rules that only appear
at scale; the photometry itself is unchanged.

* ``macro_phot.series``     — pure production rules: the (target, era,
  filter) series key, the S2 per-readout-mode saturation vetoes and their
  mapping through a server reduction, the one-provenance-per-era decision,
  the geometry verdict that refuses 8-pixel readout strips, and the
  registration policy (when a plate solution may replace astroalign, and
  when its result must be disbelieved).  Unit-tested in
  ``pipeline/tests/test_series.py``.
* ``macro_phot.register``   — the sky-chained registration route: S1 WCS
  sidecars, frame pixels -> sky -> reference pixels, and the pixel-grid
  agreement test that licenses applying a raw-frame plate solution to the
  reduced version of the same exposure.
* ``macro_phot.calib``      — local master-calibration I/O for the eras
  that have no server-reduced tree at all (ST LMi and YZ Cnc 2024): read
  the era-matched master dark and flat, apply one uniform recipe to the
  whole era, and name the recipe on every frame.

Build entry points: ``pipeline/scripts/build_s4_photometry.py`` (the
two-polar prototype) and ``pipeline/scripts/run_cv_photometry.py`` (the
production CV campaign).  Both are resumable staged CLIs whose heavy loops
are chunked, so no single invocation runs long.
"""

# Version note recorded into the photometry DB's s4_build_meta table.  Bump
# whenever the S4 logic changes in a way that alters database content.
S4_CODE_VERSION = "S4 v1.0 (2026-08-17)"
