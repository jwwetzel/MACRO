# Be-star paper — abstract v2, APPROVED by Seat 6 2026-10-05 (with its two word edits applied)

Responds to `committee/reviews/2026-10-05-seat6-bestar-abstract.md`. Venue: **PASP** (ApJ/AJ only under strategy §10.6:
a standards-epoch event with a TESS-coincident onset). Placeholders `\Num…` are emitted by script from `bestar.sqlite`.

**Title (edit 7):** Hα Equivalent-Width Monitoring of 19 Bright Be Stars with a 0.5 m Slitless-Grism Telescope,
Validated against BeSS

## Abstract

Be stars brighter than V ≈ 8 saturate the wide-field time-domain surveys, so their Hα disks are followed mostly by
amateur spectroscopists contributing to the BeSS database. Over the same seasons, a 0.5 m telescope with a slitless
grism sampled the Hα line of the 19 Be stars that BeSS shows in emission a median \NumCadenceAll{} times more often
than BeSS did, and \NumCadenceSparse{} times more often on the \NumSparse{} stars where BeSS spectra are typically
more than two weeks apart. We present \NumNights{} nights of Robert L. Mutel Telescope spectra of these stars from 2024–2026.
Equivalent widths are measured on fixed windows, with the dispersion fixed per grism and per mechanical state of the
instrument. The per-night error comes from two stability standards observed from 2025 December 5; earlier seasons
are reported but not used for claims of change. The night-to-night error floor is \NumFloor{} Å, and our equivalent
widths agree with resolution-matched BeSS measurements to \NumBessOffset{} ± \NumBessScatter{} Å over
\NumBessPairs{} contemporaneous pairs. [EVENT SENTENCE — one branch below.] [TESS SENTENCE — only under branch A,
and only if it states a result.] \NumSearched{} of the 19 stars had enough standards-epoch nights to search for
periods; with injection–recovery run through the detrending and a global false-alarm probability, the search finds
\NumPeriods{} significant periods, and at 90% completeness it would have recovered sinusoids of \NumAninetyRange{} Å
at periods of 5–20 d.

### Event sentence — branch A (\NumEventStars > 0)
\NumEventStars{} stars changed by more than three times the standards' night-to-night scatter on at least two
consecutive nights after 2025 December 5, with onsets timed to \NumOnsetPrecision{} d.
TESS sentence (A only, result required): \EventTessResult{} — e.g. "The onset in \EventStar{} falls inside TESS
Sector \NumSector{}, where the light curve \EventTessVerb{}." Deleted if no event lies inside a sector.

### Event sentence — branch B (\NumEventStars = 0), a null with its threshold
No star changed by more than three times the standards' night-to-night scatter (\NumDetThreshold{} Å) on two
consecutive nights; none was seen on the \NumStdEpochNights{} standards-epoch nights of the \NumStdEpochStars{} stars. (No TESS sentence; Figure 5 is cut.)

Word count: 235 (branch A without TESS sentence) / 246 (branch B), counted on the text above with placeholders as one word.
If branch A carries a TESS sentence (≈20 words), the clause "with injection–recovery run through the detrending and a global false-alarm probability" moves to the body to stay ≤ 250.

## Pre-declared figure cuts (edit 5)
- No standards-epoch event ⇒ Fig 5 cut; TESS sentence removed (branch B).
- \NumPeriods = 0 ⇒ Fig 6 becomes a table of a90 limits for the \NumSearched{} searched stars.
- Fig 1 (edit 3): RLMT nights against BeSS nights per star (`notes/cadence/cadence.md`), mechanical states as ticks.

Definitions (edit 2): `BeStar_Grism/scripts/be_cadence.py` header — ratio = RLMT nights / max(BeSS Hα nights, 1)
over our own season windows; BeSS-sparse = median gap between BeSS Hα dates > 14 d (or < 2 dates) in every season
window; fixed before the ratio was computed.
