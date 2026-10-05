# Seat 6 (editor): BeStar_Grism abstract and figure set, 2026-10-05 (standing rule 5)
Source: BeStar_Grism/notes/abstract/abstract_v1.md. I checked notes/step0/sample.md (19 science rows, QQ Gem = HD 46264) and notes/injection/summary.md (12 admitted to the slow tier, 4 injected but not searched, no short tier).
Decision: **APPROVED WITH EDITS.** Re-review is needed only for edit 4. Build Figs 2-4 now; Figs 1, 5 and 6 after the edits.

1. **Length and lead.** The abstract is 252 words, over the limit. Cut sentence 3 (the "ten stars first planned" history) to the body; a reader cannot interpret it. Move the cadence comparison (now the last sentence) to sentence 2: it is the paper's answer to "why care".
2. **Cadence claim.** "On the stars where BeSS is sparse ... \NumCadenceRatio times" is a selected subset. Define "sparse" in the script before computing, state how many of the 19 it covers, and also give the median ratio over all 19. \NumCadenceRatio is defined nowhere yet.
3. **Fig 1** duplicates Table 1 as drawn. Make it RLMT against BeSS nights per star (mechanical states as ticks), so that it supports the cadence sentence.
4. **\EventSentence** is undefined. Write both branches (events > 0, events = 0) into the abstract file. The zero branch must read as a null with its threshold. Seat 6 approves that text.
5. **Pre-declared figure cuts.** If no standards-epoch event survives, Fig 5 is cut and the TESS sentence either states a result or goes. If \NumPeriods = 0, Fig 6 becomes a table of a90 limits (as ruled for the Dwarf paper). The abstract keeps the completeness sentence either way.
6. **Scope wording.** Say "\NumSearched of the \NumStars stars had enough standards-epoch nights to search". The a90 range must cover only those stars.
7. **Title.** "Multi-Semester" and "Variability" claim more than one standards season and an unknown event count support. Suggested: "Hα Equivalent-Width Monitoring of 19 Bright Be Stars with a 0.5 m Slitless-Grism Telescope, Validated against BeSS".
8. **Venue.** PASP, with ApJ/AJ only under the strategy's §10.6 trigger (a standards-epoch event with a TESS-coincident onset). Agreed.

---
## v2 decision (2026-10-05): **APPROVED.** Two word edits, no re-review needed. All six figures may be built under the pre-declared cuts.
I checked notes/abstract/abstract_v2.md and notes/cadence/cadence.md. The 19-star median ratio is 7.0, and gam Cas and 28 Tau fall below 1, so the "all 19" figure is honest. Edits 1-3 and 5-8 are closed.
- Branch B: "changes that large and that persistent are excluded" overclaims, because the sampling has gaps a two-night change can fall into. Write "none was seen on the \NumStdEpochNights{} standards-epoch nights of the \NumStdEpochStars{} stars".
- Branch A: add "after 2025 December 5" after "consecutive nights", so the detection is visibly tied to the standards epoch.
- NOTE: say "typically more than two weeks apart" (the definition is a median gap).
