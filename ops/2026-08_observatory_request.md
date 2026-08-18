# Consolidated RLMT Observing & Calibration Request — August 2026

**From:** James Wetzel (Coe College), for the MACRO archival-analysis program
**To:** RLMT operations / MACRO observing committee
**Date:** 2026-08-17 (rev. 2, same day — monsoon reframing + era-resolved calibration list) · **Priority:** time-critical items marked ⏰

> **Operational reality (2026-08-17):** Winer Observatory is **offline for monsoon
> season until October 2026** — no new RLMT frames before re-opening. This request is
> therefore **submitted now to execute at the October re-opening**, and we ask that it
> be placed **first in the queue** for that re-opening. Nothing in it requires the dome
> open before October except the pyscope block-loading in Item A.3, which is
> software-side and must be verified loaded *before* re-opening.

This is one consolidated request covering every acquisition the five archival papers
(T CrB, CV time series, SN 2023ixf, Be-star grism, dwarf/AGN survey) are blocked on.
Rationale and full specifications: `ROADMAP.md` §5 (Week-1 sprint, item D1) in the
analysis repo, github.com/jwwetzel/MACRO. The calibration list in Item B is now
**era-resolved**, generated from the S0b calibration inventory
(`docs/pipeline/s0b_calibration_inventory.html`, §4).

---

## A. ⏰ Restart T CrB monitoring at the October re-opening (dark since 2025-06-24)

T CrB has now gone unobserved by the RLMT for 14 months, deep inside the predicted
eruption window, and the monsoon closure means no RLMT frame is possible before
October. The flickering science of the pre-eruption paper now exists only if frames
start flowing **the night the observatory re-opens**, and the eruption-response plan
is worthless if it isn't loaded before the nova goes — the nova will not wait for
the dome. (Coverage across the Aug–Oct gap itself is handled by external data — see
the Monsoon-gap plan below.)

1. **Nightly, from first light at re-opening** (every clear night, target airmass < 2):
   - 3 × 1 s in r (photometric anchor; short to stay below the High-Gain clip)
   - 1 short B exposure (per the T CrB strategy's exposure table)
   - 1 × 240 s lrg + 1 × 240 s hrg grism exposure (continues the 2025 Hα EW series
     homogeneously — same mode, same exposure)
2. **Weekly, resuming at the October re-opening:** one ≥ 2 hr continuous B-band run
   at ≤ 60 s cadence (flickering block; the archival data cannot supply this — ≥ 6
   such runs make the flickering section publishable). Visibility honesty, per the
   T CrB strategy's own P0-1: the field is up only in evening twilight through
   ~October, so runs resume at re-opening for that brief window, pause through solar
   conjunction, and accrue mainly when the field re-emerges in the 2027 observing
   season. The ≥ 6-run threshold is unchanged.
3. **Load the eruption-response pyscope block now, during the closure** (target of
   opportunity: the pre-agreed eruption cadence from
   `TCrB_Monitoring/ANALYSIS_STRATEGY.md` §9). This is software-side and does not
   need the dome open: **verified loaded before the October re-opening**, so the
   response triggers on night one if needed.

## B. Calibration acquisitions at re-opening (retire the archive's calibration debt — now era-resolved)

These are one-time sets; each unblocks analyses across multiple papers. The generic
list of the first revision is replaced by the **era-resolved shopping list** from the
S0b calibration inventory (`docs/pipeline/s0b_calibration_inventory.html`, §4). The
headline facts from that inventory: the **last calibration frame in the entire
archive is 2024-11-18** — every readout era active since 2025 has zero calibration
frames of any kind — and **82 of 83 science eras fail the ≥20-bias spec (164,769
science frames affected)**; only era 47 (the Andor iKon) is bias-complete.

0. **⏰ Recover the server-side master calibrations — executable NOW, no dome
   needed (2026-08-18 addition).** Reduced-frame FITS headers record full
   calibration provenance and *name the master files applied*, e.g.
   `master_bias_read0_g100_o30_2x2.fts`, `master_dark_240s_...`, living in server
   directories the headers cite verbatim: `/usr/local/telescope/rlmt/images/calibrations/`
   (2024-era) and `/mnt/ExtraImages/telescope/rlmt/images/calibrations/` (2025-era).
   None of these masters were ever synced to the archive bucket. **Request: copy
   both calibration directories (and any sibling calibration trees) into the
   `testimages` bucket now, during the closure.** If those directories survive,
   much of the archival-era shopping list below is satisfied without new sky time,
   and every applied master becomes auditable. The sky-time table below then
   shrinks to whatever the recovered directories don't cover.
   `rlmt-manifest.sqlite` `calib_gaps`, 2026-08-17 — do not hand-edit the numbers;
   camera and project labels are abbreviated here for width, verbatim strings in the
   inventory page)*:

   | Era | Camera / readout | Need | Science frames blocked | Papers affected |
   |---|---|---|---:|---|
   | 76 | Mode0 | bias × ≥20 (have 0) | 68,965 | BeStar, CV, TCrB |
   | 76 | Mode0 | flat g × ≥10 (have 0) | 40,031 | BeStar, CV |
   | 7 | High Gain | bias × ≥20 (have 0) | 22,079 | CV, TCrB |
   | 80 | Fast | bias × ≥20 (have 0) | 18,149 | BeStar, CV, TCrB |
   | 72 | 1MHz HiSens 16-bit | bias × ≥20 (have 0) | 17,414 | BeStar, CV, Dwarf |
   | 80 | Fast | flat g × ≥10 (have 0) | 13,821 | CV |
   | 1 | High Gain StackPro | bias × ≥20 (have 0) | 12,914 | Dwarf, TCrB |
   | 76 | Mode0 | dark 32 s × ≥15 (have 0) | 10,004 | — |
   | 78 | Fast | bias × ≥20 (have 0) | 8,788 | BeStar, CV, TCrB |
   | 76 | Mode0 | dark 60 s × ≥15 (have 0) | 7,964 | BeStar, CV |

   **T CrB-critical line item:** era 76 Mode0 **dark 240 s × ≥15 (have 0)** —
   blocking 2,699 frames (CV, TCrB), the readout mode of all 247 archival T CrB
   grism spectra. This is the known Mode0-dark debt of the first revision, now
   ranked in context.

   **First acquisition of the re-opening — the configuration actually running:**
   the archive's newest science frames (2026-06-28 → 2026-07-02) postdate the Fast
   eras (which end 2026-06-28) and sit in **eras 81–83**: blank READOUTM, EGAIN 56 —
   a header-convention break on the current instrument (S0b inventory §4). October's
   own science frames will land in this configuration, so take a **full
   bias/dark/flat set — including the hrg/lrg grism flats — in the as-found October
   configuration first** (the era-83 debt alone: bias × ≥20 blocking 1,831
   BeStar + TCrB frames; flat hrg 99; flat lrg 18). Please also fix — or at least
   document — the blank READOUTM/IMAGETYP header cards while the dome is open, so
   the October frames stop minting header-convention eras.

   Priority then follows the table: for the other still-installed modes (era 76
   Mode0; eras 78/80 Fast — the first revision's "era-C" grism set, which the
   eras-81–83 header state now postdates), take the frames directly at re-opening:
   ≥ 20 biases, ≥ 15 darks per listed exposure, ≥ 10 flats per listed filter,
   camera at operating temperature. These are the S0b pipeline **minima**; where a
   project strategy asks for more, **the larger number governs** — BeStar Step 2's
   standing era-C request is 50× bias + 50× dark per (mode, gain, exptime). For
   retired-hardware eras (7, 72, 1, …), the request stands only if the hardware
   answer from item 2 says the train still exists; otherwise those rows convert to
   measured cross-mode penalty terms on our side.
2. **One-afternoon hardware check** (can be run by any operator with the chair on
   the phone, any clear afternoon after re-opening): confirm the High-Gain
   digitization ceiling (~3.5 kADU clip measured in the SN data — we need the
   definitive bit-depth answer on hardware) and answer the era-B camera question:
   is the ASI/era-B imaging train still installed, or was it retired? The answer
   also dispositions the retired-era rows of the table above.

## C. Queue request (season-long, low urgency)

- **ST LMi, one more season of g/r/i** time-series blocks — season runs
  **October 2026 → spring 2027**, so it starts at the re-opening alongside Item A.
  Extends the CV paper's color baseline across its 2024-05 era seam.

## Monsoon-gap plan (Aug–Oct 2026) — what covers the closure

No RLMT action required; recorded here so the committee sees the gap is covered:

- **T CrB coverage through the monsoon window** comes from external data —
  AAVSO/ASAS-SN/ZTF/ATLAS photometry plus ARAS amateur spectra — already part of
  each strategy's external-data plan. **Those fetches start now, not in October.**
- **Eruption-during-monsoon contingency:** if T CrB erupts before re-opening, the
  eruption letter pivots to the external data (AAVSO alert stream, ARAS spectra,
  survey photometry), and RLMT joins the campaign post-restart. The pyscope
  eruption-response block (Item A.3) is loaded and verified during the closure
  regardless, so RLMT response begins the first open night.
- Everything archival in the analysis program (manifest, detector work from
  archival ladders, all five papers' reductions) proceeds through the closure.

---

*Requested by the archival-analysis program; contact James Wetzel with any
scheduling constraints. Item A supersedes lower-priority queue slots at the
October re-opening if conflicts arise — the eruption window will not wait for us.*
