"""macro_sn — the SN 2023ixf project's own pipeline stages.

Sibling package to ``macro_core`` / ``macro_phot`` / ``rlmt_diagnostics``,
built to the same standards: pure, unit-tested decision logic in a logic
module; one resumable CLI under ``pipeline/scripts``; a report renderer that
derives EVERY published number from the manifest database.

* ``macro_sn.gate0``        — Gate 0 logic (0a manifest freeze, 0b saturation
  census, 0c grism triage) plus the three verdict rules the gate exists to
  apply.  Pure functions only; no file or database access.  Unit tests:
  ``pipeline/tests/test_sn_gate0.py``.
* ``macro_sn.report_gate0`` — the Gate 0 evidence page,
  ``docs/SN2023ixf_LightCurve/sn_gate0.html``, in the house Socratic format
  (Question → Evidence → Decision → Consequence).

Build entry point: ``pipeline/scripts/run_sn_gate0.py`` (subcommands
``freeze`` / ``census`` / ``matrix`` / ``triage`` / ``verdicts`` / ``report``,
every one resumable and safe to re-invoke).

WHY GATE 0 IS ITS OWN STAGE
---------------------------
``SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md`` §4 Step 0 declares three
BLOCKING activities and forbids every downstream step until they land:
a globally deduplicated manifest freeze, a full-campaign saturation census
of the SN itself, and a contamination-hardened grism triage.  Each of the
three is a claim about frames, so each of the three has to be a query — and
a stage the provenance DAG can mark stale when its inputs move.

THE CORRECTION THIS STAGE CARRIES
---------------------------------
The strategy was written when filter slot ``6`` was believed to be a grism
because of its LABEL, and the S1 astrometry experiment excluded frames from
its candidate universe on the same label.  S2c has since measured dispersion
FRAME BY FRAME, and the slot turns out to be mixed: on this target 61 of its
83 frames are measured spectra, 3 are measured direct images and 19 are
indeterminate.  Every Gate 0 rule in this module therefore reads
``frame_dispersion`` and never a filter string — the same repair S1 v1.2
made one level up.
"""

# Recorded into sn_g0_build_meta and compared by the provenance DAG.  Bump
# it whenever a rule in gate0.py changes in a way that could alter a stored
# number, so a page built by the old rules can be told from one built by the
# new ones without reading either.
SN_G0_CODE_VERSION = "SN-G0 v1.0 (2026-08-20)"
