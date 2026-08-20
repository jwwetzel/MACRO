"""rlmt_diagnostics — S2 detector-truth extraction for the RLMT archive.

Sibling package to ``macro_core`` (same standards: pure, unit-tested logic
modules; one build/CLI entry point under ``pipeline/scripts``; a report
renderer that derives EVERY on-page number from the manifest database).

S2 answers the question every downstream paper's error model needs answered
first: *what does this detector actually do to a photon?*  Four probes, all
of them run against pixels the archive already holds (the dome is shut until
October — nothing here waits on new frames):

* ``rlmt_diagnostics.ceiling``     — where each readout mode's ADU scale
  actually ends (clip/pileup detection in science-frame histograms), the
  per-mode adopted ceiling + saturation-veto threshold, and the is-it-12-bit
  reading.  Roadmap C1/R1.
* ``rlmt_diagnostics.ptc``         — photon-transfer analysis of the
  2023-06-07 repeated darks and repeated Albireo exposures: gain (e-/ADU),
  read noise, and the StackPro variance-suppression signature.  Roadmap R2.
* ``rlmt_diagnostics.reconstruct`` — the master-reconstruction experiment:
  per-pixel linear fits of raw-vs-reduced pairs recover the effective dark
  D and flat F the reduction pipeline applied, turning unaudited reductions
  into audited ones (era 47's bias-complete iKon gives the ground truth).
* ``rlmt_diagnostics.linearity``   — counts-vs-exptime residuals from the
  2024-05-20 Vega exposure ladder and every other archival ladder the
  manifest can surface.
* ``rlmt_diagnostics.dispersion``  — S2c: per-frame FILTER IDENTITY measured
  from the pixels rather than read off the header.  Source elongation, the
  trace population, and whether those traces share a dispersion axis decide
  ``dispersed`` / ``direct`` / ``indeterminate`` per frame; trace length as a
  fraction of frame width separates the two grism units.  Settles the
  long-open slot-``6`` conflict, which turned out to be mixed rather than
  either of the two answers under debate.  Build entry point:
  ``pipeline/scripts/run_s2c_dispersion.py`` (table ``frame_dispersion``).
  Report renderer: ``rlmt_diagnostics.report_s2c`` →
  ``docs/pipeline/s2c_filter_identity.html``.  Unit tests:
  ``pipeline/tests/test_dispersion.py``.

Build entry point: ``pipeline/scripts/run_s2_campaign.py`` (resumable,
batched — every subcommand can be re-invoked and picks up where it left
off).  Report renderer: ``rlmt_diagnostics.report_s2`` →
``docs/pipeline/s2_detector.html``.  Unit tests:
``pipeline/tests/test_diagnostics.py``.
"""

# Version note recorded into the manifest's s2_build_meta table; bump when
# S2 logic changes in a way that alters stored results, so downstream stages
# can tell which rules produced the numbers they read.
# v1.1: adversarial-review round — era-79 verdict corrected (identity fit
# read as "same pixels", not "reduced = uncalibrated raw"), per-egain
# ceiling epochs recorded, detector_params uncertainties populated (incl.
# gain-bracket systematic on read_noise_e), 12-bit claim demoted to
# "consistent", dark-shot floor term quantified.
# v1.2: re-run after the S0e geometry repair, with the three products the CV
# project depends on made first-class: the empirical per-mode noise model
# (rlmt_diagnostics.noise + the `noise` subcommand — measured
# counts-vs-variance, no gain and no Poisson law assumed), round-robin
# ladder scheduling so EVERY readout mode is tried for linearity instead of
# the two richest taking all twelve slots, and blank-2026 ladders keyed on
# the canonical mode label so they can join their saturation veto.
S2_CODE_VERSION = "S2 v1.2 (2026-08-19)"
