# DwarfGalaxy_AGN_Survey — abstract rev. 2 (answers Seat 6 memo committee/reviews/2026-10-05-seat6-dwarf-abstract.md)

## Answers to the two majors
1. **Novelty.** No Hα *imaging* and no integrated Hα flux is published for any of the 19 objects:
   LVGDB (snapshot 2026-10-03) lists 9 of them with no F(Hα); the other 10 appear only in Nazarova+2025
   and the search papers it cites, none of which has Hα imaging. What IS published is spectroscopy:
   BTA/SCORPIO long-slit emission spectra for 5 objects (Dw1558+67, Dw1559+46, Dw1645+46, Dw1709+74,
   Dw1735+57; Karachentsev, Kaisina & Makarov 2025) and an SDSS fibre spectrum for Dw1459+44. So emission
   is already known in some in-band objects; what is new is the integrated (all-galaxy) Hα flux or its
   limit. "BTA depth" in rev. 1 meant the generic faintest BTA Hα limit of the Kaisin & Karachentsev
   programme (log F < −15.25), not imaging of these objects; rev. 2 drops the comparison from the abstract.
   The DW-N1 cross-match table is on disk: `DwarfGalaxy_AGN_Survey/notes/novelty/out/novelty_table.csv`
   (readable: `.../out/NOVELTY_TABLE.md`; per-object sources in `literature_status.csv` + `references.csv`;
   report `committee/work/novelty-dwarf/REPORT.md`). 17/19 = rows with a non-empty `v_hel_kms`.
2. **Venue.** AJ, short paper, ≤ 4 figures. Pre-declared fallback: if no in-band object is detected
   (N_det,in = 0) the result is a limits table and the paper becomes an RNAAS note with one figure (Fig 4).

## Calibration vs SYNTHESIS §4 (item 7)
The Cannon-curve condition is KEPT for everything that needs the band shape (band centre/edges, fluxes of
lines away from 6568 Å, band transformations — DW-P02-filter-dossier and DW-P2-band-transformations stay
blocked). The chair's brief of 2026-10-05 authorised one exception: an end-to-end line-flux scale set on
NGC 5238's published flux, stated as calibrating the throughput at NGC 5238's Hα wavelength only. Every
flux in the paper is labelled "NGC 5238 scale" and quoted only for objects whose Hα lies within
\DlamCal Å of NGC 5238's.

## Title
"Hα Imaging of Faint H I-Selected Galaxies and NGC 5238 with a 0.5 m Telescope: Detections and Limits"

## Abstract (≤250 words; \macros are script-emitted placeholders)

Of the \NinBand faint galaxies from GBT H I programme 23B-032 whose Hα falls inside our narrow-band filter,
we detect \NdetInBand in Hα and set 3σ limits of \LimRange erg s⁻¹ cm⁻² on \NnondetInBand; every limit is
compared with the flux predicted from the galaxy's GALEX FUV magnitude, and \NlimBelowPred limits lie below
their prediction. A literature cross-match shows that \Nclassified of the \NfieldsAll objects imaged already
have radial velocities (Nazarova et al. 2025), so the imaging no longer vets candidates and no integrated
Hα flux of any of them is published. The data are June–July 2023 Hα and broadband frames from the 0.5 m
Robert L. Mutel Telescope; \NfieldsHa fields have Hα. The filter transmission curve is unavailable, so we
calibrate the count rate end-to-end on the published Hα flux of NGC 5238: this sets the line-flux scale at
NGC 5238's Hα wavelength to \CalScatter, and implies an effective width of \Weff Å, but it does not fix the
band edges. With published velocities, Hα of \NoutBand objects falls outside the assumed band, and their
non-detections carry no star-formation information. For the two objects without velocities, a detection
would place them within cz ≲ \CzMax km s⁻¹; we \UnclassResult. Our median 3σ line-flux limit in 10″
apertures is \DepthMedian erg s⁻¹ cm⁻², and injected exponential discs are recovered at 50% to
\MuFifty mag arcsec⁻² (r_e = 10″). Because NGC 5238 is the flux calibrator, its Hα result is morphological
only: \NFiveTwoThreeEightResult.

## Figures (four)
2. Depth and detectability gate (Román+2020 limits; injection contours; each object's measured μ0, r_e).
3. Atlas of the 13 Hα fields (L, R, continuum-subtracted Hα), verdict per panel.
4. Hα result: flux or 3σ limit vs FUV-predicted flux, coded in band / out of band / no velocity.
5. NGC 5238 continuum-subtracted Hα map with the calibration aperture (per-night scatter in text).
Observing log/footprint = Table 1. NGC 5238 completeness (DW-P54) → data release only, no figure, no
abstract sentence. NGC 5548: no abstract sentence unless DW-P4x yields a result.
