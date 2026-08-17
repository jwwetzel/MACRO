# Consolidated RLMT Observing & Calibration Request — August 2026

**From:** James Wetzel (Coe College), for the MACRO archival-analysis program
**To:** RLMT operations / MACRO observing committee
**Date:** 2026-08-17 · **Priority:** time-critical items marked ⏰

This is one consolidated request covering every acquisition the five archival papers
(T CrB, CV time series, SN 2023ixf, Be-star grism, dwarf/AGN survey) are blocked on.
Rationale and full specifications: `ROADMAP.md` §5 (Week-1 sprint, item D1) in the
analysis repo, github.com/jwwetzel/MACRO.

---

## A. ⏰ Restart T CrB monitoring (dark since 2025-06-24)

T CrB has now gone unobserved by the RLMT for 14 months, deep inside the predicted
eruption window. The pre-eruption paper's flickering science exists only if frames
start flowing again, and the eruption-response plan is worthless if it isn't loaded
before the nova goes.

1. **Nightly** (every clear night, target airmass < 2):
   - 3 × 1 s in r (photometric anchor; short to stay below the High-Gain clip)
   - 1 short B exposure (per the T CrB strategy's exposure table)
   - 1 × 240 s lrg + 1 × 240 s hrg grism exposure (continues the 2025 Hα EW series
     homogeneously — same mode, same exposure)
2. **Weekly:** one ≥ 2 hr continuous B-band run at ≤ 60 s cadence (flickering block;
   the archival data cannot supply this — ≥ 6 such runs make the 2026 flickering
   section publishable)
3. **Load the eruption-response pyscope block now** (target of opportunity: the
   pre-agreed eruption cadence from `TCrB_Monitoring/ANALYSIS_STRATEGY.md` §9),
   verified loaded by **2026-09-01**.

## B. Calibration acquisitions (retire the archive's calibration debt)

These are one-time sets; each unblocks analyses across multiple papers.

1. **Mode0 240 s darks + biases** (the readout mode of all 247 archival T CrB grism
   spectra — zero mode-matched darks exist for them today): ≥ 15 darks at 240 s,
   ≥ 25 biases, camera at operating temperature.
2. **Flats in every currently active imaging mode** (High Gain, StackPro if still in
   use, and each grism slot): ≥ 9 twilight or panel flats per filter/mode
   combination. No StackPro or Low-Gain flat has ever been taken.
3. **Era-C grism calibration set:** bias/dark/flat matching the current (repackaged
   FITS) grism configuration, for the Be-star and T CrB reduction chains.
4. **One-afternoon hardware check** (can be run by any operator with the chair on
   the phone): confirm the High-Gain digitization ceiling (~3.5 kADU clip measured
   in the SN data — we need the definitive bit-depth answer on hardware) and answer
   the era-B camera question: is the ASI/era-B imaging train still installed, or
   was it retired?

## C. Queue request (season-long, low urgency)

- **ST LMi, one more season of g/r/i** time-series blocks (through spring 2027) —
  extends the CV paper's color baseline across its 2024-05 era seam.

---

*Requested by the archival-analysis program; contact James Wetzel with any
scheduling constraints. Item A supersedes lower-priority queue slots if conflicts
arise — the eruption window will not wait for us.*
