# Seat 4 — Detector Electrical Engineer · plan review 2026-10-03

Reviewed on disk: `products/manifest/rlmt-manifest.sqlite` (S2 tables, `detector_params`, `frames`, `eras`, `calib_gaps`), `pipeline/rlmt_diagnostics/`, `pipeline/macro_phot/{series,extract,photometry,characterize}.py`, `pipeline/macro_grism/extract.py`, `pipeline/macro_core/timing.py`, the plan ledger, `manuscripts/CV_TimeSeries/main.tex`, `ops/2026-08_observatory_request.md`, and 30-odd raw FITS frames (headers plus four quick flat-pair photon-transfer checks). Nothing was modified. The quick-look numbers below are mine, from a 1000² central region with 100-px boxes and no read-noise subtraction; they are good to a few percent and are meant to show what is measurable, not to be adopted.

**The archive holds four different cameras, not one camera in six modes.** READOUTM is a driver label; INSTRUME, geometry and the hot-pixel map identify the sensor:

| Era(s) | Camera (header evidence) | READOUTM | Dates |
|---|---|---|---|
| 1, 2, 6, 7, 8 | SBIG AC4040M (GSENSE4040 sCMOS), 4096², 9 µm | High Gain / StackPro / Low Gain | 2023-02 → 2024-03 |
| 47, 72, 74 | Andor iKon CCD, 2048², 13.5 µm, header `GAIN='4x'` | 1/3/5 MHz High Sensitivity | 2024-04 → 2024-12 |
| 76 | ZWO ASI (IMX455), 2×2, header GAIN 100 / OFFSET 30 | Mode0 | 2024-12 → 2026-03 |
| 78–83 | QHY600M (IMX455), 2×2, header GAIN 56 / OFFSET 76 | Fast, then blank | 2026-03 → |

Mode0 and Fast are different sensors: their saturated-hot-pixel positions in `s2_ceiling_frames` do not overlap.

---

## 1. Existing work

**F1 — BLOCKER (T CrB, BeStar grism). Mode0 header EGAIN is the unbinned gain; the frames are 2×2 average-binned, so the effective gain is about 1.0 e⁻/ADU, not 0.2467.**
Evidence:
- Flat-pair PTC on `Calibrations/2024-12-26/ZWO_Flat_g100_o30_day2-000{1..6}_{g,r,i}.fts` gives K = 1.01–1.06 e⁻/ADU across nine pairs (670–8,450 ADU), a factor 4.2 above header EGAIN. Bias-pair read noise is 2.1 ADU on a 303 ADU pedestal.
- Frame maxima in `s2_ceiling_frames` pile up at 16,615–16,636 ADU (Mode0) and 16,518–16,544 (Fast/blank) on repeating pixels. That is 65,535/4 plus three-quarters of the pedestal: one saturated native hot pixel averaged with three normal ones.

Consequences:
- `pipeline/macro_grism/extract.py:54-62` builds the Horne variance from EGAIN 0.2467 and 3.5 e⁻. The read term comes out 201 ADU² against roughly 4–8 measured, and the shot term is 4× too large. Weights are flattened and every grism error bar and `halpha_snr` in `g_extractions` is wrong by about 2×.
- `extract.py:64-67` sets `SATURATION_ADU = 16300` because "Mode0 raws rail at ~16.4k ADU (12-bit ADC scaled x4)". There is no such rail: the Mode0 histogram has 4–18 px per code around 16,383 and 29,175 px at 65,535. The 16.4k "rail" is the hot-pixel signature above. Valid trace pixels between 16.3k and 60k are being masked, and 40 of 193 extractions carry `n_sat_cols > 0`.
- `BeStar_Grism/ANALYSIS_STRATEGY.md:77-78` repeats 0.2467 and 1.0 as era gains.

Close: PTC per camera from on-disk flats written to `detector_params`; grism variance and saturation read from it; a test asserting predicted variance matches flanking-band variance within 20%.
Would change my mind: a Mode0 flat pair with var/mean ≈ 4.

**F2 — MAJOR (CV manuscript, S2). The gain bracket [0.60, 1.77] and the sentence that the closing measurement "requires calibration frames not yet acquired" (`main.tex:249-258`) are contradicted by the archive.**
S2 fits only dark and light pairs (`s2_ptc_fits`: all four fits return read noise 0.0 with negative intercepts, i.e. failed fits). It never uses flats. On disk: 230 High Gain twilight flats (`calib/twilightFlat/twflatall/`, 2024-01-15), 21 iKon 1 MHz flats plus 3,158 biases, and 18 Mode0 flats. Quick-look results:
- High Gain: K = 1.02–1.03 e⁻/ADU from four V pairs at 400–480 ADU, about 1.06 after removing the 3.93 ADU read noise. Header 1.057 is confirmed to a few percent.
- iKon 1 MHz 4×: K = 0.977 ± 0.005 e⁻/ADU, read noise 7.3 ADU.

So `NOMINAL_GAIN_E_PER_ADU = 1.057` everywhere (`series.py:127`) turns out roughly right for three cameras, but by accident, and the paper's read-noise uncertainty of ±2.3 e⁻ is a factor ~20 too pessimistic. For the iKon, `extract.py:129` drops the shot term entirely because EGAIN = 0.
Close: add a flat-pair PTC kind to S2, covering both EGAIN epochs of High Gain (era 2 has 17 `FLAT` frames from 2023); replace the bracket and rewrite the paragraph.

**F3 — MAJOR (S2, all projects). Linearity is not measured near any veto, and three different caps are in force.**
- `detector_params.linearity_max_dev_pct` is the *best* archival ladder per mode. The Mode0 winner peaks at 363–652 ADU, 1% of full scale. High Gain has 0 of 15 ladders under 5% residual and reports 4.2 ± 2.3%.
- The 0.92 fraction is an anecdote (`ceiling.py:63-69`).
- T CrB plans a 70% cap, SN an 80% reject (2,800 ADU), T CrB grism triage 80%, CV applied 92%.
- The ledger marks CV-P15-linearity-ladders done as "measured linearity curves".

Close: a differential test available now: for comparison stars in the CV series, regress ensemble residual against own peak ADU in bins up to the veto, per mode; adopt one cap where the slope exceeds the error floor. Add the SN 0.5 s vs 2 s pairs and the twilight-flat ramp for High Gain.
Satisfied if the residual-vs-peak slope is under 1% up to the adopted cap in each mode.

**F4 — MAJOR (S2). The "(blank 2026)" noise curve pairs frames with their own `_wcs` copies.**
7 of 16 rows in `s2_noise_pairs` have `gap_s = 0` (e.g. `…T04-26-21.fts.fz` vs `…T04-26-21_1_wcs.fts.fz`), and `s2_noise_curve` shows variance exactly 0.0 in six bins from 1,733 to 8,278 ADU. This is the current camera, the one all October data will come from. The other modes' curves are non-monotonic with log-slopes of 2.3–6.5, so above the floor they measure scene motion, not the detector.
Close: exclude alias and same-DATE-OBS pairs; replace the bright end with the flat PTC.

**F5 — MAJOR (S0/S0b era model). Eras ignore sensor temperature, gain/offset setting and camera identity.**
- Era 76 pools Mode0 at −10 °C (about 61,000 raw lights) with 0 °C (7,670, including all 6,686 from June 2025).
- 80 of the 247 T CrB 240 s grism frames were taken at 0 to −2 °C; the four 240 s master darks on disk are all at −10 °C.
- GAIN, OFFSET, SET-TEMP, COOLPOWR and the iKon preamp `GAIN` card are in the headers but not in the manifest.
- Era 76 also spans the 2025 monsoon gap and a MaxIm 6.40 → 6.30 change with no check that bias level or gain survived.

Close: re-scrape those cards into `frames`; key calibration matching on (camera, gain, offset, binning, set-point ±2 °C); show the bias level on both sides of the 2025 gap.

**F6 — MINOR (ledger truth).** `project_plan.py:675,684,1329` still block TCRB-P0-bitdepth, -ladders, -shutter-timing, -B1 and DW-P04 on "S2's tables were destroyed". The tables exist (S2 v1.2, 2026-08-20). The real blocker is that S2 does not measure gain or linearity (F2, F3).

**F7 — MINOR (CV manuscript).**
- `main.tex:260-262`: StackPro "correlates the noise between adjacent pixels". A per-pixel sum of 16 sub-reads does not, and nothing on disk measures a pixel correlation. What it does is raise read variance 16× (15.59 vs 3.93 ADU, confirmed) and hide sub-read clipping. Reword, or show the difference-image autocorrelation.
- `main.tex:327-329`: TELUT and DATE-OBS are not "independent clock cards"; both come from the same PC clock.
- `timing.py:56-70`: the error budget has no rolling-shutter term (both IMX455 cameras and the GSENSE read row-sequentially, order 0.1–0.5 s top to bottom) and no shutter term for the iKon. These are small against the paper's timing claims, but they should be stated as a bound.

**F8 — MINOR.** `photometry.py:61-65` keeps a global `PEAK_CLIP_ADU = 55000` with a comment that both cameras digitise 16 bits. It never fires on High Gain. Confirm it only selects alignment stars, or route it through `veto_adu()`.

**F9 — NOTE.**
- The High Gain clip differs by EGAIN epoch (median frame max 3,557 at 1.054 vs 3,511 at 1.057), so the two are separate configurations; the 3,200 veto is safe for both.
- With average binning, a binned peak understates the native peak by 7–16% for FWHM of 6–4 native pixels, so 0.92 leaves no margin in sub-arcsecond seeing.
- Raw frames are RICE-compressed 16-bit integers (lossless). Reduced frames are also stored as 16-bit integers, so they are quantised and re-clipped after flat division.
- 52 "Flat Field" frames in `Calibrations/2025-01` are 2 s darks (`ZWO_dark_2s_*`).
- The Dwarf strategy (line 43) calls StackPro an "average with offset stripping"; S2 measures a sum with a 1,484.5 ADU pedestal. One of them looked at a reduced frame.

**Satisfied:** the ceiling detection logic and its hot-pixel diversity gate; N_sub = 16 from three ratios; refusing to invent ceilings for Low Gain and 5 MHz; the DATE-OBS-is-start audit; the lossless-geometry resolver.

---

## 2. Plan hardening, per project

**T CrB (2/34)**
- ADD `TCRB-P0-gain-ptc`: flat-pair PTC for Mode0 and the QHY. *Accept:* K ± 3% and read noise in `detector_params`; grism variance uses it.
- ADD `TCRB-P0-temp-split`: split the 247 grism frames by CCD-TEMP. *Accept:* dark treatment stated per temperature group; flanking-band adequacy shown separately for the 80 warm frames.
- CHANGE `TCRB-P0-bitdepth`: the High Gain camera has been out of the beam since March 2024, so a hardware afternoon test may be impossible. The archive already answers it (12-bit-consistent, clip 3,511/3,557 by epoch, K ≈ 1.06, so about 3.6 ke⁻). *Accept:* closed from the archive; hardware check dropped unless the AC4040 is still on site.
- CHANGE `TCRB-P0-calib-acquisition`: Mode0 240 s darks cannot be taken unless the ASI is remounted, and one −10 °C set would not cover the 0 °C frames. Honest scope: the four −10 °C masters plus flanking-band subtraction, with a penalty term for the warm frames.
- CHANGE `TCRB-A6`: threshold from the measured 65,535 clip with a hot-pixel mask, not 16,300. *Accept:* `n_sat_cols` recomputed; masked fraction reported.
- CHANGE the ops request (`ops/…request.md:34`): "3 × 1 s in r … to stay below the High-Gain clip" is the retired camera's constraint. Re-derive exposures for the QHY at 65,535.
- `TCRB-P0-shutter-timing`: the 0.085 s θ CrB frames are electronic-shutter CMOS, so what needs testing is rolling-shutter skew and driver exposure accuracy, not shutter travel.

**CV (34/34)**
- Before release, fix F2 (gain paragraph and Table 1 read noise), F7 (three sentences) and F3 (add the residual-vs-peak panel or soften "measured linearity").
- ADD one regression: predicted vs empirical error inflation per mode with the measured gains. *Accept:* inflation factors re-quoted.
- None of this changes a timing result. It changes what the instrument section may claim.

**SN 2023ixf (5/20)**
- `SN-S2-linearity` is the single most important detector task in the portfolio and is correctly planned. ADD the twilight-flat ramp and the 2023-06-07 `cmos_tests/Albireo` set (High Gain and StackPro lights and darks, 8–128 s) as the cross-check. *Accept:* deviation vs peak ADU up to 3,400, with the 2,800 cap justified or moved.
- ADD: gain for the 1.054 epoch from the 2023 `FLAT` frames. *Accept:* K per epoch.
- StackPro frames get their own curve, since sub-read clipping under scintillation is invisible in the sum. *Accept:* StackPro vs High Gain flux ratio against peak for common stars.

**BeStar (4/24)**
- CHANGE `BE-S0-header-rescrape` to scrape GAIN, OFFSET, SET-TEMP, COOLPOWR and READOUT as well as CCD-TEMP. *Accept:* per-frame table with no nulls in eras B and C.
- CHANGE the era table: gains ≈ 1.0 (B, measured) and unknown (C, since 1.0 is a placeholder and 56 is the gain setting). "Three instruments" is literally true: iKon, ASI, QHY.
- `BE-S2-calibration`: acquire in the as-found QHY configuration first. *Accept:* flat pairs at six or more levels to 60k ADU, so the same frames give PTC and linearity.
- ADD saturated-hot-pixel mask from the 16.5k signature. *Accept:* mask per camera; applied before trace fitting.

**Dwarf/AGN (1/31)**
- `DW-P04-noise-model` is unblocked today: High Gain K ≈ 1.06, read noise 3.93 ADU, StackPro 15.59 ADU. *Accept:* depth and χ² numbers use StackPro read noise ≈ 16.5 e⁻ per pixel per frame.
- Hα at 512 s StackPro is likely read-noise-limited (16 reads per frame). State it and stop promising sky-limited depth. Resolve the sum-vs-average contradiction (F9) before P3.3 weights anything.
- `DW-P11-flats`: the AC4040 flats for slots '6' and 'W' depend on retired hardware. Assume they will never exist; delta-sky flats only.

**Legacy Rigel (0/11)**
- CHANGE `RIG-L0-header-scan` to capture INSTRUME, camera serial, READOUTM, GAIN/EGAIN, SET-TEMP, CCD-TEMP, binning, geometry and SWCREATE. *Accept:* a camera-identity census before any target census.
- ADD to `RIG-L1-calibration-census`: flat pairs per camera, so gain is measured, not read from headers.

---

## 3. Cross-cutting — top five

1. **A measured gain and read noise per camera/configuration from flat pairs, in one table that every consumer reads.** No more constants copied into `series.py`, `characterize.py` and `grism/extract.py`. Doable now for three of four cameras.
2. **Era identity = (camera, readout mode, gain, offset, binning, set-point).** Re-scrape the cards and stop treating READOUTM strings as modes of one instrument.
3. **One saturation and linearity policy**, set by a measured residual-vs-peak curve per mode and applied to native-pixel-equivalent peaks for average-binned data.
4. **October acquisition, rewritten for the camera that is actually mounted:** bias ×50, darks at the science exposures at the set-point, flat pairs at six levels per filter including hrg/lrg, with binning method, readout mode, gain and offset recorded in headers and pyscope fixed to write READOUTM, OFFSET, sub-second DATE-OBS and a real CCD-TEMP (one blank-era frame reads −0.0 for both CCD-TEMP and SET-TEMP).
5. **A bad-pixel mask per camera** (the 16.5k/33k/49k rail pixels in binned IMX455 frames, hot and RTS pixels for GSENSE), applied before detection, trace fitting and maxima statistics.

---

## 4. Verdict

| Project | Verdict |
|---|---|
| T CrB | **Execute with amendments.** F1 must be fixed before any more extraction (A2 is in progress on wrong variance and a wrong saturation mask). |
| CV | **Execute with amendments.** Three manuscript paragraphs and one figure panel (F2, F3, F7). Not satisfied on the instrument section; satisfied on ceilings, vetoes and N_sub. |
| SN 2023ixf | **Execute with amendments.** Linearity audit first, per EGAIN epoch, StackPro separately. |
| BeStar | **Execute with amendments.** Same F1 fix; era table corrected; QHY gain is unmeasured until October flats. |
| Dwarf/AGN | **Re-scope** the depth and Hα promises to the StackPro read-noise floor; noise-model task is unblocked. |
| Legacy Rigel | **Execute L0/L1 with the amended header scan.** No science decision until cameras are identified. |
