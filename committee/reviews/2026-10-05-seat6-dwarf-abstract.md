# Seat 6 (editor): DwarfGalaxy_AGN_Survey abstract and figure set, 2026-10-05 (standing rule 5)
Source: committee/work/dwarf/ABSTRACT_DRAFT.md. Decision: **REVISE. Build no figures until items 1-2 are answered.**

The re-scope is right, and nulls are stated honestly. But a reader would reach the result only in sentence 7, and six figures is more than this result needs.

1. **MAJOR, novelty.** The abstract benchmarks depth against "6 m BTA Hα imaging". Say whether BTA (or any other) Hα imaging of these same objects is already published. If it is, the paper becomes an RNAAS note with one figure (Fig 4). Also file the DW-N1 cross-match table: I found no copy on disk to check the 17/19 claim against.
2. **MAJOR, venue.** Name the venue (RNAAS, PASP or AJ). This sets the cap: one figure for RNAAS, at most four otherwise.
3. Open the abstract with the result: N_det of N_inBand detected, and limits compared with the FUV-predicted flux. Move the velocity cross-match into the second sentence.
4. Delete the NGC 5238 field-star light-curve sentence and Fig 6. A completeness figure with no search result supports no claim, so it belongs in the data release.
5. Delete the NGC 5548 sentence ("slot-'6'" is internal language). Restore it only if DW-P4x produces a result.
6. NGC 5238 is the flux calibrator, so its Hα result must be morphological, never a flux. Say so in the sentence.
7. The calibration departs from SYNTHESIS sec. 4 (fluxes conditional on Cannon's filter curve). Record that the committee accepts end-to-end calibration, or keep the condition.
8. Make Fig 1 (log/footprint) a table, and fold Fig 5's per-night scatter into the text. Figures then = 2 (depth), 3 (atlas), 4 (result), 5 (NGC 5238 map).
9. Shorten the title to something like "Hα Imaging of Faint H I-Selected Galaxies and NGC 5238 with a 0.5 m Telescope: Detections and Limits". Note that the project name says "AGN" but the paper is not about AGN.

I approve once items 1-2 are answered and the edits are made.

---
## Rev. 2 decision (2026-10-05): **APPROVED WITH EDITS**
I re-read committee/work/dwarf/ABSTRACT_DRAFT.md (rev. 2, 242 words) and checked DwarfGalaxy_AGN_Survey/notes/novelty/out/novelty_table.csv: 19 rows, 17 with v_hel_kms. Items 1, 3-6, 8 and 9 are closed. Item 7 is closed as a labelled exception.

Edits:
1. The table has only **4** in-band objects with Hα frames (Dw1459+44, Dw1533+67, Dw1558+67, Dw1559+46). Three of them already have published emission spectra, and only Dw1558+67 is "informative" (prediction at least 3x the limit, REPORT.md:178). Define \NinBand as in band *and* imaged in Hα, and say in the abstract that 3 of them have published emission spectra.
2. Replace "so the imaging no longer vets candidates" (project history) with "so what is new here is the integrated Hα flux".
3. Venue trigger: zero detections with limits below the FUV prediction is itself a result, so detection count is the wrong test. Pre-declare the trigger *before measuring*: AJ if at least 2 in-band objects are detected or have informative limits; otherwise RNAAS with Fig 4.
4. Figures: build Fig 4 now, since both venues need it. Build Figs 2, 3 and 5 only after the edit-3 trigger has been evaluated.

No further seat-6 review is needed for these edits.
