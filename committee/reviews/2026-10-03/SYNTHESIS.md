# Chair's synthesis — plan review of 2026-10-03

Seven memos are on file beside this document. This synthesis is **binding**: it amends the
six strategies and the plan ledger. Finding ids below are `<seat>.<id>` (DS = data scientist,
OA = observational astronomer, PH = physicist, DE = detector engineer, TE = telescope
engineer, ED = editor, RF = referee).

## 0. What "complete" means from here

A plan task is closed in exactly one of three ways, each recorded in the ledger with evidence:

- **done** — acceptance criterion met, evidence linked.
- **dropped** — the committee ruled it impossible or pointless on the data that exist; the ruling is cited.
- **deferred (future data)** — needs frames that do not exist yet (post-re-opening). Moved to the
  project's *2027 backlog*, which is outside the paper's critical path and outside the completion count.

A project is complete when every in-scope task is done or dropped, its result has been through the
committee, and every seat records **satisfied** (seat 7: accept or minor-with-zero-majors).

## 1. Unanimous or near-unanimous rulings

| # | Ruling | Raised by |
|---|---|---|
| U1 | **NGC 5548 broadband photometry is dead.** Slot '6' is a grism on every night. Drop DW-P41/42/44/45; DW-P03 blocker can never clear. | all seven |
| U2 | **Grism dispersion is hardware.** One solution per (grism, mechanical epoch) from hot stars; per frame only a zero point. Per-frame dispersions on disk are unphysical (−1.8 to +2.0 Å/px). | OA.E1, PH.P7, TE.F3, RF |
| U3 | **The CV paper is not finished.** The band-offset "null" is contradicted by the paper's own edges (g−i ≈ −110 ± 28 s, 9/9 or 11/12 negative). Major revision. | DS.F1, RF.B1 |
| U4 | **Ledger blockers citing "S2 tables destroyed" are stale.** Blockers must be computed from the database. | DS, OA, DE, TE, ED, RF |
| U5 | **Whole DAG reads STALE for a cosmetic reason**, hiding real staleness. Separate page digests from data digests; clean-tree rebuild; macro diff. | DS.F5, ED.E5, RF.M9 |
| U6 | **T CrB flickering (C2) and period search (C3) leave paper 1.** No ≥2 h run is possible before mid-January 2027. | DS, OA, TE, ED, RF |
| U7 | **SN 2023ixf is a +5.4 → +50 d validation/limits release (AJ/PASP)**; grism NOT PROMOTED; S7 model fit dropped; "early" leaves the title. | all |
| U8 | **Legacy archive: census only (L0–L1)**, premise corrected (≈4% Rigel; rest is the 0.5 m with six cameras), go/no-go pre-registered. | all |
| U9 | **Be-star: BeSS novelty check first**; pre-standards seasons are descriptive only. | DS, ED, RF |
| U10 | **The October request must be rewritten (rev. 3)** for the camera actually mounted (QHY600) and the sky actually available. | OA, DE, TE |

## 2. Disagreements the chair resolves by test, not by vote

**D1 — hrg dispersion: 0.47 Å/px (OA) or 1.59 Å/px (code; PH and TE reasoned from it).**
These give R ≈ 2,500 vs R ≈ 200 and decide whether line profiles / V/R are in scope. OA's evidence
(no Hβ where 1.59 puts it on θ CrB; TiO 7054 step; O₂-B depth) is specific and checkable.
*Test:* on Vega and θ CrB hrg frames, locate Hα plus ≥2 of {telluric 6277, 6867, 7186, 7594; Hβ}.
≥3 lines, residual < 1 px ⇒ adopted. All profile/velocity/V-R tasks are **held pending D1**, neither
dropped nor kept until it returns. Measured LSF is then published per grism/epoch.

**D2 — Mode0 gain: header 0.2467 (code) or ≈1.0 e⁻/ADU (DE flat pairs).** DE's evidence wins unless a
Mode0 flat pair shows var/mean ≈ 4. *Test:* S2 flat-pair PTC per camera; grism variance and saturation
threshold read from `detector_params`; the 16.3 kADU "rail" is replaced by a hot-pixel mask + 65,535 clip.

**D3 — CV band offset: astrophysical or estimator bias?** Not decided here. *Test:* per-band injection
with per-band ramp widths, signed matched-cell bias carried (RF.B2). The paper reports whichever the
test supports; the sentence "no offset" goes either way.

**D4 — Instrument paper (ED cross-cutting 1).** The editor recommends a seventh, PASP instrument/pipeline
paper written first, with SN as its validation section. This changes the portfolio and is **James's
decision**. Until decided: instrument material is built once as `docs/pipeline` evidence, each science
paper carries a condensed instrument section, and SN proceeds as its own short release.

## 3. Shared foundation (Wave 0) — every project stands on this

| Id | Task | Accept |
|---|---|---|
| F-1 | QHY-era dedup: `_calibrated` reduced twins linked to raw parents (TE.F6) | canonical 2026 frames drop ≈27k; eras 79/82 become reduced aliases |
| F-2 | Header re-scrape: INSTRUME, GAIN, OFFSET, SET-TEMP, CCD-TEMP, COOLPOWR, FOCPOS, FLIPSTAT, TELPIER, FWPOS, FWALLNAM, SWCREATE (DE.F5, TE) | columns in `frames`, no nulls where the card exists |
| F-3 | `mech_epoch` table beneath eras (camera, rotation step > 0.3°, flip, wheel map) (TE.F1) | table emitted; test: no master applied across a boundary |
| F-4 | S2 flat-pair PTC per camera/config → `detector_params`; fix alias pairs in `s2_noise_pairs` (DE.F1/F2/F4, DS.F7) | K ± 3% and RN for AC4040 (both EGAIN epochs), iKon, ASI; QHY flagged "October" |
| F-5 | One saturation/linearity policy from residual-vs-peak per mode; native-pixel peaks for average-binned data; bad-pixel mask per camera (DE.F3, OA.E8) | slope < 1% up to adopted cap per mode |
| F-6 | S2c: sky-lozenge background-morphology test; labelled truth set of 200 frames with confusion matrix (OA.E4, DS.F8) | NGC 5548 slot '6' re-issued; Wilson intervals published |
| F-7 | Provenance: data digest vs page digest; status exits non-zero only on data staleness; ledger blockers computed; `show` syncs (DS.F5/F6) | at least one stage FRESH for a verifiable reason |
| F-8 | Absolute clock from archived exoplanet transits per era (OA.E6, DS.F10, RF.M4); print the −294 ± 59 s eclipse | \|O−C\| < 120 s per era or the offset is carried |
| F-9 | Evidence snapshot: sha256 + row counts of products per release committed; "product absent = failure" test mode; manuscripts tracked decision → James (DS.F4) | manifest file in repo; CI mode exists |
| F-10 | Clean-tree end-to-end rebuild to all-FRESH at a tagged commit | zero STALE; tag recorded |
| G-1 | Grism library: fixed dispersion per (grism, mech_epoch) — resolves D1 (U2) | ≥3 lines, residual < 1 px; Hα–O₂B separation constant < 1–2% |
| G-2 | Grism variance + saturation from `detector_params` (D2) | predicted vs flanking-band variance within 20% |
| G-3 | Pixel-based identity gate (spectral fingerprint), header never sole reason (OA.E2, RF) | zero header-only rejections; false-accept rate measured |
| G-4 | Sky-lozenge background template for lrg (OA.E7) | no negative continuum; method difference < 3% |
| G-5 | Focus-offset and CCD-TEMP regressors; measured LSF per grism/epoch/focus (TE.F4, DE.F5, PH.P8) | LSF table published; off-nominal nights flagged |

## 4. Per-project amendments

### CV_TimeSeries — major revision (reopened; was 34/34)
ADD, all blocking release:
- **CV-R1-band-offset** — paired/permutation test, combined-era estimate, scatter-based errors (DS.F1, RF.B1).
- **CV-R2-bias-injection** — per-band and per-era (incl. a 2024 High Gain night) injection; signed matched-cell bias; no `max(χ²ν,1)`; no transported budget (RF.B2/B3, DS.F3).
- **CV-R3-oc-refit** — per-band constants, night-level epochs (N = 17), era-offset nuisance; Ṗ quoted each way; sensitivity to dropping 2024 (DS.F2, RF.M3).
- **CV-R4-superhump-wording** — excluded-amplitude statement; cite floor; strike "measurement of absence" and "is real — and is not orbital" (PH.P1, RF.B4/M7).
- **CV-R5-pdot-physics** — relabel as spin/spot-longitude bound; three physical scales; longitude stability in degrees, by accretion state (PH.P2, RF.M3).
- **CV-R6-colour-result** — quantify the g−r swing (amplitude, phase, two-era repeatability); survive a ≤120 s pairing window; cyclotron context (ED.E2, RF.M1, PH.P6).
- **CV-R7-edge-fit-disclosure** — model written down, grid, χ²ν distribution, best and worst fit (RF.M2).
- **CV-R8-counts** — 8,716 staged / 7,802 measured; S/N floor on "measurement"; EU UMa i/r out of the count; FWHM > 5″ cut (RF.M5/M6, OA).
- **CV-R9-reduction** — reduction paragraph; raw-vs-reduced night per era; flat-field ramp test vs detector x; mech_epoch per series (OA.E5, TE.F5).
- **CV-R10-instrument-section** — measured gains replace the bracket; StackPro and clock-card sentences corrected; rolling-shutter bound; linearity panel (DE.F2/F3/F7).
- **CV-R11-clock** — print the eclipse residual; transit clock test (F-8) (RF.M4, OA.E6).
- **CV-R12-literature** — every dataset/catalogue/ephemeris cited; comparison-with-previous-work subsection; ≥ ~40 refs (ED.E1, RF.M8).
- **CV-R13-restructure** — abstract ≤ 250 words leading with the colour result; ≤ 18 pp; process prose to appendix; internal language removed; venue AJ (ED.E2/E3, RF minors).
- **CV-R14-rebuild** — chain FRESH at a clean commit; macro-by-macro diff of `numbers.tex` filed (ED.E5, RF.M9).
- **CV-R15-release-readiness** — authorship policy, ORCIDs, Zenodo DOI, availability statements, external reader (ED.E4) — *needs James*.

### TCrB_Monitoring — execute with amendments (highest priority)
- ADD **TCRB-N1-novelty-table** first: spectra/month RLMT vs ARAS vs Asiago, Feb–Jun 2025; re-scope rule stated (ED).
- ADD **TCRB-N2-eruption-contingency** skeleton (ED); **TCRB-P0-eruption-block** artefact (TE).
- ADD **TCRB-P0-mech-epoch**, **-gain-ptc**, **-temp-split** (TE, DE).
- CHANGE A0 (pixel gate, G-3; counts by S2c verdict; 402 unique frames not 471); A3 (fixed dispersion, G-1); A6 (saturation from measured clip); A7 (θ CrB floor split-half + line-free window; focus regressor; extraction-method difference as EW floor).
- ADD **TCRB-A0b** — extract-or-reject the slot '6' and 'W' spectra 2023-05 → 2024-03 (possible only pre-dip Hα points) (DS, OA).
- DROP A4 as written (θ CrB response → absolute flux). REPLACE with **TCRB-A5b**: F(Hα) = EW × continuum from contemporaneous AAVSO/zero-order photometry; EW and flux both plotted against orbital phase (ellipsoidal ±10%) (OA, PH).
- CHANGE A5: ARAS spectra degraded to measured LSF must reproduce native EW (PH). ADD pre-registered EW-change detection rule before A5 runs (DS).
- HOLD pending D1: profile morphology / wing velocities.
- CHANGE Phase B: anchors are **B-only** (23 frames) plus singletons — R and I are clipped; scripted peak-at-target census of all 224 imaging frames (OA.E3).
- CLOSE from archive: P0-bitdepth. UNBLOCK: -ladders, -shutter-timing (rolling-shutter skew), B1. P0-calib-acquisition → flanking-band method stated; bench-dark question to site.
- DROP C3; C1 → one sentence/table; DEFER C2 and the October/Jan 2027 splice (second instrument: QHY) to the 2027 backlog (U6).
- CHANGE D4: ten figures → six.

### SN2023ixf_LightCurve — re-scope (U7)
- Re-run Gate 0 on fresh S0; fix 438/439 arithmetic; name the unnamed rule; S2c-measure all 632 broadband frames.
- CLOSE SN-venue-decision (AJ/PASP) and SN-G0c (NOT PROMOTED; frames go in the release).
- SN-S2-linearity is the key detector task: add twilight-flat ramp, 2023-06-07 Albireo set, per-EGAIN-epoch gain, StackPro separately (DE).
- CHANGE S4: held-out check stars; S5: scintillation term; S8: nightly means, whole-night bootstrap; intra-night periodogram DROPPED.
- SPLIT S6: **S6a flash phase** via (H − "[S II]") differential colour with filter widths from zero-point ratios (~65 Å), behind a **predicted-excess gate** (PH); **S6b ejecta phase** dropped unless a transmission curve arrives by 2026-11-15.
- DROP S7. ADD residuals-vs-published table before any variability limit. Template table with camera/mech_epoch/FWHM/focus (TE). Peak epoch read from literature; "rise" language removed if clean start is at peak (PH). Figures twelve → five.

### BeStar_Grism — execute with amendments (U9)
- **BE-S-1a-bess first and alone**; ADD **BE-N1-gate** (< 2 verified-active emitters ⇒ material moves to instrument section; science paper stops).
- BE-S-1b injection–recovery per star through detrending, before periodograms are opened.
- Inherit G-1…G-5. Era table: five mechanical states, gains measured. HOLD V/R pending D1.
- BE-S6: telluric H₂O regressor (7200 Å band depth) (PH). BE-S8: season 1 EW only. BE-S10: detection rule only from the standards epoch (2025-12-05); earlier seasons descriptive.
- BE-S11 short tier: three qualifying stars, global FAP; λ Eri short search dropped. Dither test in lieu of grism flats → 2027 backlog. Disposition column for 308 indeterminate + 24 direct frames. QQ Gem dispositioned. Figures twelve → six.
- DEFER season-2 observing items to the 2027 backlog; include in ops rev. 3.

### DwarfGalaxy_AGN_Survey — re-scope (U1)
- NGC 5548 leaves the title. DROP P41/42/44/45, P03 for those frames. DW-P02-slot6 closes with a per-night verdict table (F-6). ADD **DW-P4x**: one-night broad-Hα extraction triage (2 days; EW and nightly scatter; no lags) and **DW-P4y** zero-order photometry feasibility — either earns a paragraph or is dropped.
- ADD **DW-N1-novelty** (literature cross-match per candidate) before any stacking; ADD **DW-P36-0** depth-vs-expected-L(Hα) table (PH).
- DW-P04 unblocked with StackPro RN ≈ 16.5 e⁻; resolve sum-vs-average; depth promises re-stated.
- DW-P11: no AC4040 flats will ever exist → night-sky superflat from the June 2023 L frames (OA).
- DROP P53, P56 (and P52 unless threshold/trials pre-declared with injection); keep P54 completeness for NGC 5238 only. W dropped from NGC 5238 surface photometry.
- Hα products: detection/non-detection is a velocity statement; non-detections carry the out-of-band caveat (PH). Fluxes conditional on Cannon's filter curve — *needs James to contact Cannon*. Figures fourteen → six.

### Legacy archive (was "Legacy_Rigel") — census only (U8)
- Premise and name corrected (TE.F10). L0: dedup/collision reconcile; camera-identity header scan (DE) ; mech-epoch census (seven cameras).
- L1: overlap query first (T CrB, ST LMi, YZ Cnc 2015–2022); per-target nights × filters × longest run; calibration census with flat pairs; per-season clock audit; contact binaries with minima in ≥3 seasons.
- L2: go/no-go criteria **pre-registered before the census is seen**; must name a question, a reader and a venue, else outcome = data release note.

## 5. Standing rules adopted (all projects)

1. Every χ²ν is reported per band/era/mode with its dof; χ²ν < 0.5 is a defect equal to χ²ν > 2. No `max(χ²ν, 1)`.
2. Every null carries the effect size it would have recovered and a "predicted scale" beside it. "k of N significant" is banned.
3. Injection grids report the signed matched-cell bias, never a median of absolute values.
4. Saturation is judged at the target, per frame, in native pixels.
5. A ≤ 250-word abstract with number placeholders is approved by seat 6 before figures are built.
6. EW is a ratio: wherever the continuum varies, carry a continuum light curve.

## 6. Items only James can close

- Decide D4 (instrument paper).
- Authorship policy with the consortium; real ORCIDs; Zenodo DOI; an outside reader for the CV draft.
- Contact Cannon (filter transmission curves; candidate dossier; authorship).
- Send ops request rev. 3 to Winer; ask whether the ASI and AC4040 cameras still exist; ask for the server `calibrations/` trees and a hardware change log.
- Whether `manuscripts/` should be tracked (private repo or branch) so drafts are versioned.
