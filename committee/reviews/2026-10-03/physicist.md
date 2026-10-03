# Seat 3 — Physicist · plan review 2026-10-03

Read-only review of the tree at `27c962f-dirty`. Order-of-magnitude arithmetic is mine and shown. Literature values quoted from memory are marked *(verify)*.

## 1. Existing work

### CV_TimeSeries (`manuscripts/CV_TimeSeries/main.tex`)

**P1 — MAJOR. The YZ Cnc "measurement of absence" is logically inverted.**
`main.tex:1139-1144`, `1273-1277`, `1505-1511`: the blind-search 90% recovery contour is 82–237 mmag (`p4_outburst.amp90_blind`), the literature superhump semi-amplitude floor is 50 mmag (`final_science.py:100-105`), and the text concludes that because "the contour sits above that floor" a superhump "would have been visible". A recovery contour is the *smallest* amplitude the search recovers. A contour above the floor means superhumps of 50–82 mmag were invisible in every run, and anything below 237 mmag was invisible in the worst one. The search was sensitive only to large, freshly developed superhumps (full amplitude ≳0.16–0.47 mag). It is also moot physically: the runs are in normal outbursts, where no superhump is expected.
*Close:* state it as "sensitive to semi-amplitudes ≥ X in N of 5 run-filters; smaller superhumps not excluded", delete "measurement of absence" in all three places, and change the abstract's "null with a measured contour behind it" accordingly. *Would change my mind:* evidence the contour is an upper envelope, not a recovery threshold.

**P2 — MAJOR (framing of the headline number). The Ṗ bound is correct arithmetic with the wrong physical label and no scale.**
Arithmetic checks: 3.6×10⁻⁹ × 3.156×10⁷ s = 0.11 s/yr; P/|Ṗ| = 6833 s / 3.6×10⁻⁹ = 6.0×10⁴ yr; the 3σ value 2.5×10⁻⁹ agrees with my estimate (~2.2×10⁻⁹) for 36 epochs of 84 s over ±4,344 cycles. Good. But `main.tex:865-867` calls it "s of orbital period per year". The timed feature is the bright-phase falling edge, i.e. the accretion region rotating behind the white-dwarf limb. It clocks white-dwarf spin plus spot longitude, not the orbit, and the VSX period it is compared with is itself photometric. The paper gives no physical scale at all (no mention of gravitational radiation, synchronism or spot migration anywhere in the text). Scales: secular orbital evolution below the gap is |Ṗ| ~ 10⁻¹³–10⁻¹⁴ (τ ~ 10⁹ yr), 4–5 orders below the bound; an asynchronous polar relaxing like V1500 Cyg has Ṗ_spin ~ 4×10⁻⁸ *(verify)*, which the bound does exclude; an 84 s rms is 0.012 cycle = 4.4° of spot longitude.
*Close:* relabel as a bound on the bright-phase (spin/spot) period; add one paragraph giving the three scales; state the result the data actually support — spot longitude stable to ~4° over 1.9 yr, and whether that holds across accretion states (the state tags exist in `p3_state_night`). That is a more interesting sentence than the Ṗ limit.

**P3 — MINOR. The band-offset null has no expected size beside it.** `main.tex:903-906`, `1221-1223` say an offset "is expected on cyclotron grounds". Limb occultation of a compact region is achromatic to first order; chromatic shifts are bounded by the ingress duration (spot extent ~0.01–0.03 cycle, ~70–200 s) and are plausibly tens of seconds. Bounds of 134–313 s with 300 s exposures cannot test that. *Close:* quote the expected scale and say the bound is not yet constraining.

**P4 — MINOR.** `main.tex:688-694`: the folded-profile gradient at phase 0.16–0.19 is called independent corroboration of the 1,071 s offset. It is the same data folded on the same ephemeris; it checks the edge fitter, not the offset. Also, the offset has no physical meaning until the VSX fiducial (maximum? linear-polarisation spike?) is identified. *Close:* call it a consistency check; name the fiducial or say it is unknown.

**P5 — MINOR.** `main.tex:1069`, `NumSuperoutburstAmp = 3.0`: for YZ Cnc ~3 mag is roughly a *normal* outburst amplitude; superoutbursts are nearer 3.5–4 *(verify against AAVSO)*. The classification rests on the AAVSO bracketing, which is sound; the constant should be corrected or cited.

**P6 — NOTE.** The "headline physics figure" (colour vs phase) carries no physical inference. One paragraph tying the sign and size of the colour change to cyclotron harmonic number at ST LMi's ~12 MG field *(verify)* would answer "why would a physicist care".

**Satisfied:** timing chain, cycle-count uniqueness, the refusal to publish per-cycle errors, flickering floors, the 0.384-cycle hump phase shift.

### Grism stage G (`g_extractions`, `pipeline/macro_grism/`)

**P7 — MAJOR. Per-frame dispersion is not a physical quantity and the stored values prove it.** A grating's Å/px is fixed by hardware to ~1%. `g_extractions.disp_a_per_px` for accepted T CrB hrg frames is bimodal at 0.45–0.47 and 1.54–1.80; the ratio is 3.39, exactly the O₂ A/B lever ratio (`wavelength.py:17-19` describes the misidentification as solved; the stored rows still carry it). lrg scatters 0.23–2.0. Three mutually inconsistent lrg values are on disk: ~1.1 (stored mode), 1.9–2.2 (`wavelength.py:86`), and ~3.2 implied by S2c's "factor of two the optics predict" relative to hrg. EW in Å scales linearly with dispersion, so every lrg EW is currently uncertain by ×2–3.
*Close:* solve dispersion once per (grism, era) on hot stars (Vega, Spica, θ CrB) from Hβ, Hα, O₂-B, O₂-A; fit only a zero point per frame; require every T CrB frame's Hα–O₂ separation to agree within 1–2%.

**P8 — MAJOR (plan-level). Delivered resolution is R ≈ 160–300 on T CrB, not 1000–2000.** `halpha_width_px` is 13–26 px on hrg (`wavelength.py:57` says ~20); at 1.59 Å/px that is 21–41 Å, ~1000–1900 km/s. Trace FWHM of 5–9 px (`extract.py:44`) gives ≥8–14 Å even for short exposures. T CrB's Hα is a few hundred km/s wide: profiles, absorption reversals and wing velocities are unresolved. See §2.

### SN Gate 0 (`docs/SN2023ixf_LightCurve/sn_gate0.html` §4)

**Satisfied, with a sanity check passed.** r ≈ 12 at +2.5 d gives ~8 photons s⁻¹ Å⁻¹ into an effective ~700 cm²; a ~65 Å filter, 64 s and a 4-px FWHM predict a peak of ~1–2 kADU. Measured: 825–1,883 ADU. Right order.
**P9 — NOTE.** Strategy §3.3's "~3.7 ke⁻ full well" is a 12-bit ADC clip, not a well; linearity close to the clip is expected and the 46 "suspect" frames are probably recoverable.

## 2. Plan hardening

### TCrB_Monitoring — execute with amendments
- **DROP** hrg profile morphology and differential wing velocities (strategy line 80, Fig. 9; folds into TCRB-A5/D4). Impossible at R ~ 200 (P8). Honest scope: EW and line centroid relative to O₂-B. *Accept:* Fig. 9 removed or replaced by a measured-LSF figure.
- **ADD TCRB-A5b: EW → line flux via the red-giant continuum.** EW = F_line/F_cont, and F_cont at 6563 Å is the ellipsoidally modulated giant: ~±10% at 113.8 d *(verify R-band amplitude)*. The series spans 124 d = 1.09 ellipsoidal cycles, so raw EW contains a guaranteed geometric oscillation, and Rank 1 asks precisely whether recovery shows "oscillations". *Accept:* EW multiplied by a contemporaneous R/I continuum light curve (AAVSO, or zero-order differential photometry from the same grism frames) and both series plotted against orbital phase.
- **ADD TCRB-A3b:** fixed dispersion per grism (P7). *Accept:* per-frame Hα–O₂ separations constant to 2%.
- **CHANGE TCRB-A4:** θ CrB is a B6 star; second-order 3300–4500 Å light lands on first-order 6600–9000 Å, while T CrB has almost no blue flux. The response redward of Hα is therefore colour-dependent. *Accept:* second-order fraction measured (hot star with and without an order-sorting filter, or modelled) before any absolute flux is quoted. Note also that `g_extractions` finds an Hα *emission* peak at SNR 220–420 in θ CrB lrg frames: A1 will find the calibrator in emission.
- **CHANGE TCRB-A5:** Munari's windows were defined at echelle resolution; at 20–40 Å the 6651–6714 Å TiO heads bleed toward the red window. *Accept:* ARAS spectra degraded to the measured LSF reproduce their native EW within the stated tolerance.
- **NOTE:** October spectra are taken at airmass ≳2.5 in twilight and are not homogeneous with 2025; tag them.

### CV_TimeSeries — execute with amendments (P1, P2 before release)
- **ADD CV-P5-superhump-wording** (P1). *Accept:* no sentence claims superhumps were excluded below the contour.
- **ADD CV-P5-pdot-context** (P2). *Accept:* the bound is labelled spin/spot, with the three scales and the degrees-of-longitude statement.
- P3–P6 as wording changes.

### SN2023ixf_LightCurve — execute with amendments
- **CHANGE SN-S6: split it.** *S6a flash phase (+2.5/+5.4/+6.4 d): unblock.* The line is then narrower than the filter, so only the filter's equivalent width is needed, and that is already measurable: catalogue zero points give R − H = 3.2–3.4 mag (`frames.zmag`: R 21.5–21.8, H 18.3–18.4, "1" 18.2, O 18.0), i.e. W ≈ 1380 Å/20 ≈ 65 ± 15 Å for all three narrowbands. Expected signal: L_Hα ~ 10³⁹ erg s⁻¹ *(verify, Jacobson-Galán+23)* against L_λ ≈ 2×10³⁸ erg s⁻¹ Å⁻¹ at r = 12 gives EW ~ 5–10 Å, an excess of ~8–15% at +2.5 d falling to a few per cent by +6 d. That is detectable only differentially. *S6b ejecta phase stays blocked* on the transmission curve.
- **ADD SN-S6-0: use "1" ([S II]) and O as the continuum.** A Type II at <+10 d has no [S II] or [O III] emission; those filters are pure continuum, taken at the same exposures with the same star ensemble. (H − "1") against field stars is a near-differential Hα EW estimator, far better than g/i interpolation. From ~+15 d the broad Hα red wing (±8000 km/s = ±175 Å) enters the 6724 Å filter, so it is valid early only. *Accept:* (H − "1") colour of the SN minus the field-star locus at the three flash epochs, with errors, compared with synthetic photometry of WISeREP spectra through 50/65/80 Å top-hats.
- **ADD SN-S6-1: predicted-excess gate.** *Accept:* if the predicted excess is <3× the narrowband systematic at every epoch for every plausible width, S6 is demoted now rather than after the work.
- **CHANGE strategy line 12:** "g-band peak at ~+8–9 d" looks late; published optical maximum is nearer +5–7 d *(verify)*. If so the clean start at +5.4 d is at peak and the paper contains no rise. *Accept:* peak epoch read from the external light curves and all "rise" language removed.
- **SN-G0c grism:** worth the timebox, but the same order of magnitude applies: a 5–10 Å line at ≥20–30 Å resolution is a 15–30% bump on one resolution element, on an H II arm. The SNR > 10 criterion is right;

### BeStar_Grism — execute with amendments
- **DROP the V/R branch now** (strategy lines 30, 126; Fig. 12). Its own decision rule is FWHM < 4 Å; measured ≥8–14 Å. Keep the O₂-referenced centroid/asymmetry moment, which is position-independent in slitless data and plausibly good to ~10 km s⁻¹.
- **CLOSE BE-S-1c-hrg-bandpass from existing evidence:** O₂-B at +193 px and O₂-A at +650 px on era-B hrg (`wavelength.py:17-19`). *Accept:* one Be-star frame confirming, then done.
- **ADD to BE-S6: telluric H₂O.** Line 134 scales O₂ by airmass only. Water lines inside 6520–6610 Å total ~0.1–0.25 Å and scale with PWV, not airmass; that equals the claimed 0.1 Å floor and exceeds the <0.05 Å window-robustness criterion. hrg reaches the 7200 Å H₂O band. *Accept:* per-frame 7200 Å band depth carried as a regressor; EW of θ Vir shows no residual correlation with it.
- **CHANGE BE-S5/S4:** dispersion per (grism, era) as P7; quote resolution measured, not expected.
- Disk radii: Grundstrom & Gies only, inclination stated; no peak-separation radii.

### DwarfGalaxy_AGN_Survey — re-scope Phase 4; execute Phases 0–3 with amendments
- **Phase 4 as written is impossible.** S2c §4 and `frame_dispersion` show 91 of 143 NGC 5548 slot-6 frames are spectra, 17 direct images on two nights. DW-P41/P42/P44/P45 (`project_plan.py:1428-1443`) and strategy lines 40, 83–87 still assume a luminance light curve. F_var and a structure function from 2–3 epochs are not statistics. *Honest scope:* a slitless spectral series on ~11 nights inside the Xi et al. 2025 season.
- **ADD DW-P4x: broad-line series.** NGC 5548's broad Hα (FWHM ~5000–9000 km/s = 110–200 Å at 6676 Å) is resolved even at 30 Å. I estimate S/N ~5 per pixel on the continuum per 256 s frame and tens per night on the line. If the grism covers [O III] 5007, normalise to it (constant narrow line); if it is the Hα unit, there is no internal calibrator and the product is profile width/EW only. *Accept:* one night extracted; Hα EW and its nightly scatter reported before more effort. No lags.
- **ADD DW-P4y:** zero-order differential photometry feasibility (AGN versus field-star zero orders). *Accept:* comparison-star zero-order rms < 2% or the idea is dropped.
- **ADD DW-P36-0: Hα depth before stacking.** With H zero point 18.3 mag for 1 ADU s⁻¹ and W ≈ 65 Å, 1 ADU s⁻¹ ≈ 8×10⁻¹⁵ erg s⁻¹ cm⁻². Seven 512 s frames and a 10″-radius aperture give a 3σ limit of ~10⁻¹⁴ erg s⁻¹ cm⁻² before continuum-subtraction systematics: L ≈ 10³⁸ erg s⁻¹ at 10 Mpc, SFR ≈ 6×10⁻⁴ M☉ yr⁻¹, about ten O7 stars. That is 10–100× shallower than the BTA Hα programme *(verify)* and above the expected SFR of most M_B ≈ −10 candidates. *Accept:* predicted limit versus expected L_Hα per candidate tabulated; fields where the limit cannot separate dIrr from dSph are labelled uninformative.
- **CHANGE DW-P36 logic:** candidates have no velocities. A detection through a ~65 Å filter is a velocity statement (cz ≲ 1500 km/s) and is the valuable outcome. A non-detection cannot distinguish a quenched dwarf from a background galaxy out of band, and must not be tabulated as an "SFR limit" without that caveat. [N II] is negligible at these metallicities. DW-P36 no longer needs Cannon's curve for order-of-magnitude limits; it needs it for the centre wavelength.

### Legacy_Rigel — execute ingest only; no science plan yet
- **ADD to RIG-L1-overlap a physics criterion:** contact binaries with minima in ≥3 seasons. A dP/dt of 10⁻⁷ d yr⁻¹ accumulates ~7 min of O−C over ten years, which ~1 min timings detect. That, with RLMT-era extension, is the one physically motivated product visible from here. *Accept:* count of such systems in the census.

## 3. Cross-cutting — top five
1. **Dispersion is hardware.** One solution per (grism, era) from hot stars; per-frame zero point only (P7). Three projects depend on it.
2. **Publish the measured LSF** per grism/era/exposure regime and delete every plan item that needs R > 500 (T CrB profiles, Be V/R, SN narrow-line profiles).
3. **Narrowband equivalent widths from zero-point ratios** (~65 Å for H, O, "1"): a facility number that unblocks order-of-magnitude Hα work in SN and Dwarf without waiting for a data sheet.
4. **Every null gets its expected signal beside it** — a "predicted scale" column in each project's limits table (CV Ṗ and band offset, SN flash excess, Dwarf Hα). P1 would have been caught by this.
5. **EW is a ratio.** Wherever the continuum source varies (T CrB ellipsoidal, Be outburst continuum, AGN), carry a continuum light curve and report line flux as well.

## 4. Verdict
| Project | Verdict |
|---|---|
| TCrB_Monitoring | **Execute with amendments** — drop profile/velocity science; add continuum correction and fixed dispersion |
| CV_TimeSeries | **Execute with amendments** — P1 and P2 fixed before release; otherwise satisfied |
| SN2023ixf_LightCurve | **Execute with amendments** — split S6, add the [S II]-continuum estimator and predicted-excess gate |
| BeStar_Grism | **Execute with amendments** — drop V/R, add H₂O regressor |
| DwarfGalaxy_AGN_Survey | **Re-scope** Phase 4 (spectra, not photometry); Phases 0–3 execute with the depth gate |
| Legacy_Rigel | **Execute ingest only**; no go/no-go before the census |
