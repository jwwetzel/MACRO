# Seat 7 — The Difficult Referee · plan review 2026-10-03

Read cold: the CV manuscript sources, Figs 6 and 9, then the databases under `products/` (read-only). Line numbers are `main.tex` unless stated.

---

## 1. Existing work — referee report, CV manuscript

**Summary.** The paper presents RLMT photometry of five CVs and claims (i) an ST LMi O−C null with |Ṗ| < 3.6e-9, (ii) no band dependence of the bright-phase edge epoch, (iii) a "measured" superhump null and a hump upper limit for YZ Cnc, (iv) colour–phase "results" for ST LMi, and (v) a fully traceable release. Every number I recomputed reproduces; the inferences drawn from them do not all survive.

**Verdict: MAJOR REVISION.** Four blockers, all in the abstract or conclusions; none is fatal to the data set.

### Recomputation (from `products/phot/cv_timeseries.sqlite`)

| Published | Recomputed | Status |
|---|---|---|
| 36 epochs / 17 nights, RMS 84 s, χ²ν 0.92 (ν=35), median σ 84 s | 36 / 17, 83.96 s, 0.9225, 83.6 s (`p3_oc_night`) | match |
| Ṗ = 1.4σ; bound 3.6e-9; 0.11 s/yr; 6.0e4 yr | own weighted quadratic: 1.14e-9 ± 0.82e-9; 3.61e-9; 0.114; 5.99e4 | match |
| Refit P = 0.07908902 ± 6.4e-8 d | 0.07908902 ± 6.39e-8 | match |
| Tie median 25 / 64 mmag; 3 / 1 of 25 at goal | 24.9 / 63.8; 3 / 1 (`cv_cattie`) | match |
| 8,716 frames / 108 nights | 8,716 / 108 rows in `cv_frames` | count matches, **description does not** (M5) |
| g−i = −97 ± 51 s, 1.9σ, "no offset" | same 12 pairs, empirical scatter: **−108 ± 27 s, t = −4.0** | **number matches, conclusion does not** (B1) |
| All 307 macros vs `p5_number` | 0 differences | match |

### Blockers

**B1. The band-dependence null is contradicted by the paper's own released edges.** (l.118 "likewise shows no offset"; l.912 "none reaches even 2σ"; Conclusion 5.) From `p3_edge`, the 12 Mode0 cycles timed in both g and i give t_g − t_i with mean −108 s, standard error 27 s, 11 of 12 negative (t = −4.0, Wilcoxon p = 0.002; ≈0.007 after a three-pair trials factor). The published ±51 s exists only because `band_difference` uses the transported injection budget and `max(chi2nu, 1)` (`pipeline/macro_phot/phase3.py:1055`), while the pair's own χ²ν is 0.26 — the differences scatter half as much as their assigned errors. Corroboration: High Gain G−I is −137 ± 88 s, same sign; and the published per-night epochs give weighted means g = −42 ± 24 s, i = +44 ± 21 s (11 of 12 i epochs positive), i−g = +87 ± 32 s. *Closing test:* report the scatter-based paired statistic beside the budget-based one; fit per-band constants in the O−C (χ² falls 32.3 → 24.4 for 3 parameters); either claim the offset or show it is estimator bias (B2). *Changes my mind:* a per-band injection run showing a differential recovery bias of order 100 s.

**B2. The "bias floor of 5.5 s" (l.631) is not what the injection grid measured.** `p3_sigmat` at the matched cell (shape and depth known) has bias −76.7 s (g), −54.8 s (i), −52.1 s (r); the grid reaches −115 s. The 5.5 s is the median of |bias| over all 48 cells per band, including mis-specified cells where bias happens to cross zero (`pipeline/scripts/run_cv_phase3.py:1477–1491`). Consequences: (a) the systematic term in every published error bar is understated ~10×; (b) the bias is band-dependent (g−i ≈ −22 s) and depth-dependent, so it feeds B1 and could differ between the 2024 and 2025 eras, which is exactly a Ṗ-like signature; (c) "the residuals are the demonstrated timing error and nothing else" (l.1466) is unsupported. *Closing test:* carry the signed matched-cell bias per band and era; rerun the quadratic with an era-offset nuisance term and report how the Ṗ bound moves.

**B3. The "independent" check on the transported error budget is the same budget.** l.643 calls it "the fitter's Monte-Carlo sigma"; l.1216 "an independent per-night measurement" on 27 Mode0 epochs; Conclusion 4 rests on it. `sigma_t_mc_s` is one constant per series — 126.34 s for every g edge, 90.47 r, 96.80 i — read from the injection grid's most optimistic cell of the single night 2025-02-27 (`run_cv_phase3.py:1256–1260`). It is neither per-edge nor independent of the budget it is checking. *Closing test:* delete the claim, or run a real per-edge bootstrap and an injection grid on a 2024 night (the paper already says that is what would settle it).

**B4. The superhump "measurement of absence" is logically inverted.** l.1141–1144 and Conclusion 7: a 90% recovery contour of 82–237 mmag "sits above" the 50 mmag floor of published superhump semi-amplitudes, "which is what makes 'no superhump' a measurement". A detection threshold above the smallest expected signal means such signals would *not* have been found. "If one had [fallen in a superoutburst], a superhump would have been visible" (§5.3) holds only for amplitudes above 82–237 mmag, on 5 of 18 run-filters, and the 50 mmag figure is uncited. *Closing test:* state the excluded amplitude range, cite the floor, drop "twice over".

### Majors

**M1. The abstract promises "colour results" (l.97) and the paper states none.** The "headline physics figure" paragraph (l.824–832) describes the pairing rule and the tie bar only. Fig. 6 shows g−r swinging by ~0.75 mag; no amplitude, phase or interpretation is given. A 600 s pairing window is 0.09 cycle, wider than the edge itself, so colour excursions at ingress/egress may be non-simultaneity artefacts. *Close:* quantify the colour curve; show it survives a ≤120 s window.

**M2. Edge-fit quality is undisclosed.** The model is never written down (l.565 "a parametric model"). Accepted ST LMi edges have χ²ν median 38, maximum 1,890, 55 of 63 above 10, on 6–12 points at ~219 s cadence, with four trial widths. *Close:* state model, grid (±3 cadences about the catalogue prediction, 3.65 s step), χ²ν distribution, and show a representative and a worst fit.

**M3. The Ṗ bound is not interpreted and probably not orbital.** 3.6e-9 is orders of magnitude above any secular expectation and no literature value is compared. The timed feature is an accretion-spot edge that can migrate in longitude with accretion state; the 2024 lever arm is nine epochs in different bandpasses with transported errors and one fallback epoch; 2025-12-17 is a single epoch. *Close:* state what the bound constrains (spot-longitude stability) and its sensitivity to dropping 2024.

**M4. The clock validator's numbers are withheld.** l.338–343 says the eclipsing-binary bound "is weak". `s3_clock_eclipses`: one usable eclipse, O−C = −294 ± 59 s, clock bound 4,517 s. A 5σ, five-minute residual on the only external check of a timing paper must be printed and explained.

**M5. "8,716 usable light frames" (l.268; abstract l.73) is false as worded.** `cv_frames.status`: 7,802 matched on 104 nights; 528 failed registration, 287 "frame does not contain the reference field", 99 excluded. *Close:* say 8,716 staged, 7,802 measured.

**M6. The 6,676 "catalogue-tied measurements at 9–77 mmag" include noise.** EU UMa i: 24 points, catalogue magnitudes 19.0–26.0, all with σ > 0.2 mag, 21 below 3σ; EU UMa r: 25 of 45 below 3σ; g: 34 of 172 with σ > 0.2. `ch_noise_series` has no precision at all for EU UMa i and r, so the quoted range silently omits series whose points it counts. *Close:* an S/N floor on "measurement", or a precision range that covers every counted series.

**M7. YZ Cnc hump: "is real — and is not orbital" (l.1123)** for a modulation that fails the paper's own detection test in 0 of 6 scopes. Strike or demonstrate.

**M8. Bibliography.** Nine `\cite` calls, five science references. No citation for REFCAT2, ZTF, ASAS-SN, VSX, the ephemerides, any prior ST LMi/YZ Cnc/AN UMa timing or colour work, the ZTF polar folds asserted at l.172, or the superhump floor.

**M9. Every stage behind the paper reads STALE or STALE_UPSTREAM** (`check_pipeline_status.py`, today): `stage_cv_timeseries` changed after CV-S4; CV-S5, which supplies the abstract's precision, ran 2026-08-19, before the 2026-08-20 photometry rebuild and records `cv_selection` 34 → 33. Build commits are `2c897e7-dirty` and `f0bba55-dirty`. *Close:* clean-tree rebuild end to end, diff the 307 macros.

### Minors

1. Cycle-count uniqueness rests on an *assumed* σ_P (last digit of a VSX string, l.728); say so in the text, as the database note does.
2. "SUBSTITUTE FOR THE PLANNED FIGURE" (captions 6, 7), "CV-S9 graded…" (Fig. 9 footer) and "ANALYSIS_STRATEGY §4" are internal language a journal reader cannot resolve.
3. Fig. 6 lower-right ordinate label collides with tick labels.
4. Per-band error budget is not validated per band: i epochs scatter 45 s against 57–99 s assigned; g 109 s against 56–112 s. χ²ν = 0.92 is an average of two mis-calibrations.
5. Plan ledger names the injection night 2025-02-28; database and paper say 2025-02-27.

**Satisfied with:** traceability (307/307), the tie treatment, the saturation vetoes, the refusal to publish per-cycle error bars, the AN UMa grading.

---

## 2. Plan hardening, per project

**CV_TimeSeries (34/34 — not finished).**
- ADD `CV-R1-band-offset`: scatter-based paired test, per-band O−C constants. *Accept:* abstract states the offset or a demonstrated bias explanation.
- ADD `CV-R2-injection-2024`: injection grid on a High Gain night and per band, signed bias carried. *Accept:* no epoch carries a transported budget.
- ADD `CV-R3-clean-rebuild`: clean-commit rebuild of S0→CV-S11. *Accept:* zero STALE stages, macro diff published.
- CHANGE `CV-P3-yzcnc-superhump` verdict wording to an excluded-amplitude statement.
- DROP EU UMa i/r from the measurement count.

**TCrB_Monitoring.** As the referee who will get this paper:
- `TCRB-A3-wavelength` is not credible on disk: `g_extractions` per-frame dispersions on *accepted* T CrB frames span 0.45–1.81 Å/px (hrg) and 0.23–2.00 Å/px (lrg); the two grisms' means are indistinguishable (0.92 vs 0.83). The anchors are being misidentified. CHANGE: one dispersion per grism per era, fixed from the calibrator, per-frame zero point only. *Accept:* per-frame dispersion scatter < 2%.
- DROP wing-velocity evolution (Fig. 9 of the plan) unless measured resolution supports it; Hα is 15–22 px wide. *Accept:* delivered FWHM in km/s stated before any velocity is plotted.
- 8 of 35 accepted sample frames have saturated columns; flanking vs master-dark extraction differ by 3% (hrg) and 11% (lrg). ADD both as explicit EW systematics. *Accept:* EW error floor ≥ the measured extraction-method difference.
- 29 of 34 gate rejections are `header_off_target`, the very header the strategy calls untrustworthy. CHANGE the gate to reject on the Gaia pattern only. *Accept:* zero rejections whose sole reason is a header.
- θ CrB is a Be star and same-night on 17 of 60 nights; calibrator gate passes 10 of 18. Absolute flux is a Mar–Apr appendix at most.
- Honest scope: an EW series Feb–Jun 2025 with ARAS cross-validation. Flickering (C2) is impossible until 2027; DROP from this paper. The unique-frame count is 402, not 471 (`project_counts`: 69 within-tree duplicates).
- Ledger blockers citing "S2's tables were destroyed" are stale: `s2_*` tables exist (S2 run 2026-08-20). Re-evaluate, do not wait.

**SN2023ixf_LightCurve.**
- Gate 0 already answered: `sn_g0_verdict` grism NOT PROMOTED (0 extracted spectra, no wavelength source, no contamination test), venue AJ/PASP, no transmission curve. `SN-S6-halpha-curve` and the physical Q2 are pointless without the curve: DROP unless the curve arrives by a fixed date.
- 438 of 632 broadband frames clean, first at +5.4 d. "Early light curve" is not a supportable title. CHANGE title and scope to a post-+5.4 d validation and limits release.
- `SN-S7-model-consistency`: no data before +5.4 d, no UV. DROP, or one sentence.
- ADD: residuals against published photometry with colour-term validity range before any variability limit is quoted. *Accept:* per-band offset and RMS table exists before S8 starts.
- Gate 0 itself is STALE; re-run before anything cites it.

**BeStar_Grism.**
- Run `BE-S-1a-bess` before anything else; the ApJ/AJ decision and half the plan hang on a two-day task still pending after seven weeks.
- No grism-era calibration frames exist; era B camera may no longer exist. *Accept:* era-B EW floor measured from surrogate darks is published, or era B is labelled lower-bound-only.
- Inherit the T CrB dispersion finding: `BE-S4-wavelength` must demonstrate stable dispersion per (filter, era) before any EW. *Accept:* standard-star EW scatter reported per era.
- DROP `BE-S11` short-period tier for anything but the three qualifying stars; DROP season-to-season claims for 69 Ori and 5 Cnc (pre-standards). `BE-figures` and `BE-draft` should not start until Step −1 returns ≥2 verified emitters.

**DwarfGalaxy_AGN_Survey.**
- **Phase 4 (NGC 5548) is dead on disk.** `frame_dispersion`: of 143 slot-'6' frames, 91 dispersed, 35 indeterminate, 17 direct on 3 nights. Slot 6 is the grism. A 15-night broadband light curve does not exist; "zero WCS solutions ever" is the symptom, not a queue problem. DROP DW-P41/42/44/45; CHANGE `DW-P03` blocker text; mark `DW-P02-slot6-dispersion` done. *Accept:* ledger shows NGC 5548 removed from the paper's title and outline.
- Hα fluxes (`DW-P36`) need a transmission curve nobody has asked for. *Accept:* Cannon contacted this week, or Q1 becomes detection/non-detection only.
- NGC 5238 slot W: one 0.4 s flat. DROP W from surface photometry now.
- Honest scope: Hα detections/limits for 13 fields plus NGC 5238, vetting not discovery.

**Legacy_Rigel.** No strategy exists; correct. Keep to L0–L1 census only. ADD: a duplicate/collision audit (`legacy_manifest_BAD_collisions.csv` exists) and a header-time convention audit before any target census is believed. *Accept:* go/no-go is decided on a census with a stated dedup rule.

---

## 3. Cross-cutting — top five

1. **The whole DAG is stale and three build stamps are `-dirty`.** No number from any project is citable until a clean-tree rebuild passes. *Accept:* `check_pipeline_status.py` all FRESH at a tagged commit.
2. **Error bars that only inflate.** `max(chi2nu, 1)` hides over-estimated errors and manufactured the B1 null. Every stage reports χ²ν and a scatter-based alternative.
3. **Injection grids are summarised by medians over mis-specified cells.** Carry the signed matched-cell bias; never a median of absolute values.
4. **Grism wavelength solutions are per-frame and unstable** (factor 4–9 on disk). Fix dispersion per grism per era before T CrB or Be-star EWs are computed.
5. **Plan ledger drifts from disk** (stale S2 blockers, slot-6 task still "in progress", NGC 5548 still planned). The ledger should be regenerated from database verdict tables, not edited by hand.

---

## 4. Verdict

| Project | Verdict |
|---|---|
| CV_TimeSeries | **Major revision** — not complete; execute CV-R1–R3 before any external circulation |
| TCrB_Monitoring | **Execute with amendments** — EW series only; fix dispersion first |
| SN2023ixf_LightCurve | **Re-scope** — AJ/PASP validation and limits release; drop Q2 physical curve and model fits |
| BeStar_Grism | **Execute with amendments** — Step −1 gates before any pipeline work |
| DwarfGalaxy_AGN_Survey | **Re-scope** — remove NGC 5548; Hα limits conditional on the filter curve |
| Legacy_Rigel | **Execute census only**; no science plan |

I am not satisfied on the CV manuscript. I will be when B1–B4 are closed and M1–M9 are fixed or rebutted in writing. Nothing above requires new observations.
