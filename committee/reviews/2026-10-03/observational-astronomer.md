# Seat 2 — Observational Astronomer · plan review 2026-10-03

Read-only review. I opened 25 FITS frames, 6 extracted spectra and the manifest/CV databases; nothing was modified.
Frame paths are relative to `rlmt-archive/`. "Manifest" = `products/manifest/rlmt-manifest.sqlite`.

## 1. Existing work

**E1 — BLOCKER. The adopted hrg dispersion (1.59 Å/px) is wrong by ~3.4×; the O₂-"pair" proof is a TiO coincidence.**
`pipeline/macro_grism/wavelength.py:17-24` asserts the +650 px dip is O₂-A and the +193 px dip O₂-B. In an M4 III the +193 px dip is the TiO λ6651 head and the +650 px dip is O₂-B λ6867: (6867−6563)/(6651−6563) = 3.47, inside the 6% tolerance around 3.387 (`:83`). Evidence for ≈0.47 Å/px:
- θ CrB (B6, no TiO) `products/grism/spectra/mjc_tet_CrB_hrg_60s_2026-01-10T12-42-26_spec.fits`: Hα absorption at x≈2390, one telluric dip at +635 px, and **no Hβ/Hγ** at −1069/−1398 px where 1.59 Å/px puts them on-chip with 13–24 kADU of continuum.
- T CrB hrg (2025-03-20, 04-21): 40% flux step at +1040 px = TiO λ7054 at 0.47 Å/px; the "O₂" dip is 17–21% deep (B-band depth, not A-band); Hα is 13–20 px wide (≈300–450 km/s at 0.47; ≈1,000–1,500 km/s at 1.59).
- `g_extractions`: only 6/26 T CrB hrg and 2/37 lrg "halpha+o2" rows land near the adopted values; signs flip frame to frame and the same frames give different dispersions under `flanking` vs `masterdark`. On θ CrB the code reports an "Hα emission" anchor that does not exist.
*Close:* solve dispersion once per (grism, era) on hot stars already on disk (Vega 2026-04-15 hrg/lrg, θ CrB) using Hα + telluric 6277/6867/7186/7594; per-frame only the zero point. *Changes my mind:* Hβ visible at −1069 px in any A/B-star hrg frame.
Upside if I am right: hrg delivers R≈2,500–3,000 (21 km/s per px), far better than either strategy assumes.

**E2 — MAJOR. The "identity gate" is a header test and discards good T CrB spectra.**
`gate.py:231-242` rejects on header offset before looking at pixels; the content arm compares one coordinate of one star within ±400 px under two parities (`:73`). 55/193 validation rows are `header_off_target`. I opened two of them, `rawimage/2025-03-11/ext_T_CrB_{lrg,hrg}_240s_…`: bright late-type spectra with a strong emission line, header RA 05 29 Dec −17 — a position 39° below the horizon at that time. The header is stale; the telescope was on target (the docstring says as much, `:37-39`). *Close:* spectral-fingerprint gate (Hα emission + TiO 7054 step at fixed pixel offset; flux level for θ CrB), re-adjudicate all 21 header-bad T CrB frames, measure false-accepts on non-T CrB grism frames.

**E3 — MAJOR. T CrB photometric anchors: I is dead as well as R; only B survives.**
Strategy ruling 2 (`TCrB_Monitoring/ANALYSIS_STRATEGY.md:45`) plans anchors "on B and I". Pixels at the target:
- R 8 s High Gain `2023-05-13/mjc13301`: flat-topped 3,518–3,565 ADU over 8 px. StackPro 8 s `2023-05-20/mjc14007`: flat-topped 13.9–14.0 kADU.
- I 1 s `2024-03-12/xwg0720b`, `2024-03-19/xwg07911`: flat-topped 3,383–3,476. I 0.25 s `2024-03-22/xwg08214`: 2,950 (84% of clip, fails the 70% cap). The 4 s frames of 2024-03-05 look safe only because ZMAG is 18.6–19.1 (2.6 mag of cloud).
- Clean: B 16 s StackPro `2024-02-13/irm0440a` (6.8 k of a 56 k ceiling); R 0.25 s (1,064); G 1 s (1,648).
S2 already measured the ceiling (`s2_ceiling_modes`: High Gain clip 3,496, 12-bit, veto 3,200), yet the ledger still blocks five tasks on "S2's tables were destroyed" (`project_plan.py:675`). *Close:* scripted peak-at-target census of all 224 imaging frames; restate scope as B-only (23 frames, 2024-02-12→03-12) plus singletons.

**E4 — MAJOR. S2c's `direct` verdict is not a safe photometry certificate on slot '6'.**
NGC 5548 slot-'6' frames get three verdicts under one filter and exposure (direct 17 / dispersed 91 / indeterminate 35; mixed within 2023-04-05 and 04-22). I rendered `mrf09600` and `mrf1100b` (both "direct"), `mrf08706` ("dispersed"), `mrf08306` ("indeterminate"): all four show the same vignetted dispersed-sky lozenge plus zero-order stop image seen in T CrB lrg frames. Slot '6' is a grism on all 16 nights. *Close:* add a background-morphology test (T CrB lrg: corners 306–310 ADU, centre 1,170) and re-issue; I expect 143/143 dispersed.

**E5 — MAJOR (CV manuscript). The reduction is not described and its provenance contradicts ROADMAP convention 1.**
`cv_frames`: 5,407 of 7,802 matched frames are `server_reduced`; 2,395 are local raw+masters with flats 105–202 days from the science night. `manuscripts/CV_TimeSeries/main.tex` §2–3 never mentions flats, darks, the reduced tree, exposure times, the 4.0″ fixed aperture (`photometry.py:44`) or the annulus. *Close:* one reduction paragraph plus per-series recipe/flat-age columns; one night per era photometered from both raw+masters and reduced, differential light curves agreeing within σ_chk.

**E6 — MAJOR (CV timing). A common clock offset is not excluded at the size of the headline number.**
`p3_cycle_count`: ST LMi mean O−C 1,071 ± 97 s; EU UMa 1,060 ± 84 s — two independent ephemerides, same offset to 11 s — while the absolute clock is bounded only to ±75 min (`docs/pipeline/s3_timing.html` §4). I did rule out one mechanism: pyscope filenames carry the *scheduled* time and DATE-OBS lags them by 64→1,200 s within 2025-02-27, but O−C residuals stay flat across that night and match the 2024 MaxIm era, so DATE-OBS is the shutter time. **Satisfied on that point.** S3 does not audit this lag; record it. *Close:* time one archived transit per era (WASP-43 b 2024-02 and 2026-01, TrES-3 b 2024-03/05, WASP-52 2025-10); |O−C| < 120 s closes it.

**E7 — MINOR.** lrg flanking background over-subtracts (extracted flux to −1.7 kADU against a 6 kADU continuum; `debt_median_rel_diff` 7–16%) because the sky is a sharp-edged lozenge. Fit a lozenge template before flanking.
**E8 — MINOR.** Mode0 bin-2 is an *average*: hot pixels pile at ≈16.6 kADU = (65,535+3×304)/4 in every frame, so one saturated native pixel is invisible to the 60,200 veto. Also ~100 late-season T CrB grism frames were taken at 0 °C against a −10 °C master dark (2025-11-22, 20 frames); dark level is negligible (master median = pedestal 303), so only a per-temperature hot-pixel mask is needed.
**E9 — NOTE.** `ops/2026-08_observatory_request.md:90` still lists the 240 s Mode0 dark as "have 0"; `calib_gaps` now says `master_only`. No hrg/lrg flat exists in any era.

**Satisfied:** the High Gain 12-bit ceiling (confirmed in pixels); S3 mid-exposure BJD from DATE-OBS; CV cloud veto and the decision not to apply k″.

## 2. Plan hardening

What the sky allows (my astropy run, Winer, airmass < 2; hours with Sun < −12° / < −18°):

| Target | Oct 5 | Oct 15 | Nov 1 | Dec 15 | Jan 1 | Jan 15 | Mar 1 |
|---|---|---|---|---|---|---|---|
| T CrB | 1.1 / 0.6 | 0.6 / 0.2 | 0 | 0.2 / 0 | 1.4 / 0.9 | 2.4 / 1.9 | 4.8 / 4.3 |
| ST LMi | 0 | 0.3 / 0 | 1.6 / 1.2 | 5.1 / 4.6 | 6.3 / 5.8 | 7.2 / 6.7 | 9.2 |
| λ Eri | 4.0 | 4.8 | 6.1 | 6.2 | 6.2 | 6.2 | 2.8 |
| M101 | 0 | 0 | 0 | 3.0 / 2.5 | 4.2 | 5.1 | 7.6 |
| NGC 5548 | 0 | 0 | 0 | 1.9 / 1.4 | 3.1 | 4.1 | 6.5 |

**T CrB**
- CHANGE P0-restart: October gives ≤2 weeks of marginal data, none from ~Oct 25 to ~Dec 12; the real restart is the morning apparition. Grism only with Sun < −15° (twilight fills the sky lozenge). *Accept:* request states dates, Sun-altitude limit and readout mode.
- CHANGE the nightly block (`ops:34`): "3 × 1 s r to stay below the High-Gain clip" clips in High Gain (E3) and is scintillation-limited (~6 mmag at 1 s, ~12 mmag at 0.25 s). Use the as-found 16-bit mode, 5 × 10 s mildly defocused. *Accept:* target peak 15–40 kADU, check-star rms < 5 mmag.
- ADD a non-Be standard at T CrB's brightness: BD+33°2642 (CALSPEC, V 10.8, 7° away), one 240 s lrg+hrg per T CrB night. *Accept:* same-night standard on ≥80% of nights.
- ADD A3a: dispersion constants from hot stars (E1). *Accept:* ≥3 lines per grism, residual < 1 px.
- CHANGE A0 per E2. *Accept:* every reject justified by pixels.
- ADD: extract the T CrB slot-'6' and 'W' frames (S2c: 12/16 and 4/8 dispersed; 2023-05 → 2024-03). They may be the only RLMT Hα points before and inside the dip. *Accept:* EW ± error, or a documented saturation failure.
- DROP A4 (θ CrB response → absolute flux). θ CrB is variable, 0.6–5 s vs 240 s, and Hα-contaminated. Replace with F(Hα) = EW × continuum flux from contemporaneous AAVSO/own photometry. *Accept:* photometric source named per epoch.
- DROP C3 (period search): B exists on three nights. C1 reduces to two B snippets — a sentence, not a table. B6: delete the 3–5 mmag nightly promise.
- C2: ≥2 h dark runs are possible only from ~Jan 10; six weekly runs complete mid-March at the earliest, so the Dec–Jan submission target excludes flickering. Decide now: submit the spectroscopic paper without it.

**CV_TimeSeries**
- ADD reduction paragraph and raw-vs-reduced night (E5); ADD transit clock test (E6). *Accept:* as stated there.
- CHANGE ST LMi request (`ops:125`): it cannot start "at re-opening" (first useful ≥1.5 h night ≈ Nov 1). Ask for whole-orbit blocks (≥2 h), one filter per night alternating with gri nights, flats within 7 days.
- MINOR: 16 matched frames have FWHM > 5″ inside a 4″ aperture; cut or show target-minus-check is flat against FWHM.

**SN2023ixf**
- CHANGE scope/title: first clean broadband night is +5.4 d (`sn_g0_verdict`); G/R/I are clipped at 3.4–3.55 kADU on +1.6 to +3.5 d. Only the narrowband frames at +2.5 d reach the early phase, and their bandpasses are unknown. *Accept:* no "early/flash photometry" claim unless S1-narrowband-curves closes.
- CHANGE S5: 0.5–2 s exposures carry 5–11 mmag scintillation and under-exposed comparison stars. *Accept:* error model has an explicit scintillation term matched by check-star rms.
- G0c: slot-'6' frames reach the clip on every night (max peak 3,513–3,548) and M101's disc is dispersed across the trace. Keep the timebox; expect a methods footnote.
- ADD (if the 2023 camera survives, ops B.2): deep gri M101 templates from mid-December. *Accept:* field-star subtraction residual < 2%.

**BeStar_Grism**
- CHANGE BE-S4/S5: inherit E1; if hrg ≈ 0.47 Å/px, V/R is in scope (`ANALYSIS_STRATEGY.md:30` assumes it is not). *Accept:* measured FWHM of O₂-B per night.
- CLOSE BE-S-1c-hrg-bandpass: O₂-B sits ≈ +650 px from Hα, on-chip.
- ADD to the ops request (it contains no Be-star item): season-2 restart now — λ Eri, 69 Ori, 5 Cnc are up 4–6 h from October — with HR 1544 (π² Ori, beside λ Eri) as the autumn nightly standard and η Hya from December. *Accept:* a standard on every science night.
- ADD dither test in lieu of grism flats: same star at three detector positions. *Accept:* EW agrees within the floor.

**DwarfGalaxy_AGN_Survey**
- DROP P41, P42, P44, P45: NGC 5548 photometry is impossible (E4; other filters total four frames on one night). P03 plate-solving of those frames is pointless. Optional two-day triage: a broad-Hα EW series from the 143 grism frames.
- ADD: night-sky superflat from the June 2023 L frames (≈400 frames, 19 fields, two weeks) — the only flat those data will ever have if the camera is gone. *Accept:* residual large-scale structure < 0.3% of sky.
- DROP P53, P54, P56: 7–9 epochs over two weeks cannot support a period search or eclipse timing.

**Legacy_Rigel**
- ADD to L1: per-season clock audit from archived contact-binary minima (*accept:* O−C within 60 s of literature); flat/dark census per filter per run; short-exposure shutter check for the PL16803.
- CHANGE L1-overlap: query T CrB, ST LMi, YZ Cnc first — a 2015–2022 BVRI T CrB series would be worth more than the 2023–24 anchors.

## 3. Cross-cutting — top five

1. **Issue rev. 3 of the observatory request before the dome opens**: as-found mode named explicitly and a fixed set-point all season; the visibility table above; standards (BD+33°2642, HR 1544, Vega while it is an evening object); a compact PN (IC 418, NGC 7027) through both grisms; bias ×50 and darks at the exposures and temperature actually used; flats in every filter; a dome-lamp linearity ladder; one SDSS standard field for colour terms; fix READOUTM/IMAGETYP and the stale RA/Dec cards.
2. **One grism calibration library** — dispersion constants, sky-lozenge background template, pixel-based identity. T CrB, BeStar, SN and Dwarf all inherit E1, E2, E4, E7.
3. **Absolute clock from exoplanet transits per era**, replacing the 75-min AG LMi bound.
4. **Decide the reduced-tree question**: audit it per era and amend convention 1, or re-reduce the 5,407 CV frames from raw.
5. **Saturation is judged at the target, per frame, in native pixels** (E3, E8), and the ledger's stale S2 blockers are cleared.

## 4. Verdict

| Project | Verdict |
|---|---|
| TCrB_Monitoring | Execute with amendments — E1 and E2 first; photometry reduced to B-only; flickering split off |
| CV_TimeSeries | Execute with amendments — not ready to submit until E5 and E6 close |
| SN2023ixf_LightCurve | Execute with amendments — a +5 d onward gri light curve, AJ/PASP |
| BeStar_Grism | Execute with amendments — after the grism library is fixed; restart observing now |
| DwarfGalaxy_AGN_Survey | Re-scope — drop NGC 5548 photometry and the variability census |
| Legacy_Rigel | Execute L0–L1 only; no science commitment before the census |
