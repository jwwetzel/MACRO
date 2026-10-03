# Seat 1 — Astronomical Data Scientist · plan review 2026-10-03

Reviewed on disk at `27c962f` (read-only). Queried `products/manifest/rlmt-manifest.sqlite` and `products/phot/cv_timeseries.sqlite`; ran `pytest pipeline/tests` (1784 passed, 1 skipped, 225 s); ran `check_pipeline_status.py` and `update_project_plan.py show`. Every number below is from those queries unless a file:line is given.

## 1. Existing work

**F1 — BLOCKER (CV paper). The headline inter-band null is contradicted by the paper's own epoch table.**
The abstract (`manuscripts/CV_TimeSeries/main.tex:117-125`), §timing (`:617`) and §limits (`:1221-1223`) state that the bright-phase edge "shows no offset" between bands: "0 of 32 band pairs significant", strongest pooled pair g−i = −97 ± 51 s (1.9σ), reported as an upper limit. The rule is a 3σ vote count on budget error bars (`pipeline/scripts/run_cv_phase3.py:1396-1419`).
From `p3_oc_night` (ST LMi, 36 epochs), pairing g and i epochs on the same night:

| test | result |
|---|---|
| nights with both g and i | 9 (2 High Gain 2024, 7 Mode0 2025) |
| sign of g−i | negative on **9 of 9** (sign test p = 0.004) |
| paired mean, empirical s.e. | **−123 ± 30 s** (t = −4.2, p = 0.003) |
| same, weighted by the paper's own budget errors | −106 ± 38 s (2.8σ) |
| per-band mean O−C | g −49 ± 30, r −14 ± 16, i +51 ± 13 s (monotone in wavelength; i positive in 11/12) |
| unpaired g vs i | Welch p = 0.008, Mann–Whitney p = 0.005 |

Both camera eras give the same sign and size (pipeline's own pooled rows: era 7 G−I −137 s, era 76 g−i −97 s) but are never combined. The pooled rows carry χ²ν = 0.26, 0.31, 0.18 on 12 cycles: the pair error bars are inflated ~2×, which is what turns a ~3–4σ offset into "1.9σ". The null is a product of the error model, not of the data. Bonferroni over three band pairs leaves p ≈ 0.01.
I do not claim it is astrophysical. The injection grid itself (`p3_sigmat`) shows edge-time bias from −115 s to +11 s under shape mismatch, while the budget carries a bias floor of 5.5 s (`p3_oc_night.sigma_floor_s`); fitted ramp widths are band-dependent and sit on grid nodes (g 97.8 s, r 273.9 s, i 223.9/447.9 s). So the offset is either a wavelength-dependent edge (the cyclotron effect the paper says it could not detect) or a ~100 s band-dependent estimator bias. Either reading falsifies a printed sentence.
*Close:* (a) add a paired, distribution-free test across nights and a combined-era estimate with scatter-based errors; (b) run the injection grid with per-band ramp widths to bound the estimator's band bias; (c) rewrite abstract/§5 accordingly. *Changes my mind:* a label-permutation test on same-cycle pairs giving p > 0.05, or (b) showing a ≥100 s differential bias (then it is a systematic, and the 5.5 s floor must go).

**F2 — MAJOR (CV paper). The O−C χ² and Ṗ rest on pooling bands that differ by ~100 s.**
χ²ν = 0.92 on 35 dof (`numbers.tex:161-162`) is a compensating mixture: about per-band means, i gives χ² = 3.8/11 and r 5.0/9 (errors over-stated ~1.7×), g 15.7/12. Per-band constants lower χ² by 7.9 for 3 parameters. Band mix changes with time (2024 epochs are I-heavy and positive; the last epoch is g-only), so a band offset aliases into the refitted period and the 1.4σ Ṗ. Also: 36 epochs on 17 nights are counted as 36 independent data.
*Close:* refit period and Ṗ with per-band constants and with night-level epochs (N = 17); report both. *Satisfied if* Ṗ limit moves < 20%.

**F3 — MAJOR (CV paper). "Monte-Carlo error on every edge" is one number per series.**
`p3_edge.sigma_t_mc_s` is constant within a band (g 126.34, i 96.80, r 90.47 s for every 2025 edge; NULL in 2024). The "independent" χ² of 0.52 (`main.tex:644-652`) therefore tests a per-series constant drawn from the same injection night, not per-edge errors. Per-edge fit χ²ν are 2–1890 on 6–12 points with ≥4 parameters; 6 of 63 O−C values are exact duplicates of another cycle (grid quantisation, up to ~32 s steps in era-7 I).
*Close:* state it as a per-series budget; add grid-step to the error budget; publish the per-edge χ²ν distribution.

**F4 — MAJOR (portfolio). Nothing that backs a claim is version-controlled or snapshot.**
`.gitignore:5,20` excludes `manuscripts/` and `products/`. `git ls-files manuscripts products` is empty. The 32-page draft, `numbers.tex` and both SQLite evidence bases live only in Dropbox, with live `-wal/-shm` files beside them. `numbers.tex` records "commit 1380616" but that commit contains none of its inputs. On a clean checkout the ~130 product/manuscript tests skip (`test_cv_products.py:32`, `test_manuscript_cv.py:32,65`) and the suite still reads green.
*Close:* content-hash manifest of products per release (sha256 + row counts) committed; manuscripts tracked; a CI mode where "product absent" is a failure, not a skip.

**F5 — MAJOR (pipeline). The freshness signal has saturated and carries no information.**
Every stage reads STALE/STALE_UPSTREAM. The root is S0e, stale only because `docs/pipeline/s0e_geometry_fix.html` was rewritten by the site rebuild (27c962f) — a cosmetic change that propagates to all 40 stages. Substantive staleness is hidden inside it: S2c ran 2026-08-19T00:55, S0 rebuilt `frames` at 01:02 (S2c, G, R-S0 report real digest changes). Meanwhile the plan ledger still prints CV 34/34 "done" because `sync` last ran 2026-08-20; under its own rule every one would flip to `redo_needed`.
*Close:* separate report-page digests from data digests; status must exit non-zero only on data staleness; `show` must run `sync` or print its age.

**F6 — MAJOR (plan ledger). Blockers contradict the database.**
Five tasks are blocked because "S2's tables were destroyed by the S0 table swap" (`project_plan.py:675` and siblings: TCRB-P0-bitdepth/-ladders/-shutter-timing, TCRB-B1, DW-P04). The tables exist: `s2_ptc_fits` 4 rows, `s2_linearity_ladders` 66, `detector_params` 82, built 2026-08-20; `carry_sibling_tables` (`build_s0_manifest.py:611`) now preserves them. Conversely the brief says SN Gate 0 was not re-run, yet `sn_g0_build_meta` is stamped 2026-08-21T03:33, after S1 (02:45).
*Close:* blockers derived from a query, not prose.

**F7 — MAJOR (S2). The PTC does not measure what downstream uses.**
`s2_ptc_fits` covers only High Gain and StackPro; read noise is 0.0 in all four fits (railed); "light" gains are 0.074 and 0.0028 e⁻/ADU. Mode0 (69,805 light frames), Fast (52,516) and 1 MHz (19,310) have no fit; CV photometry ran on a nominal 1.057 e⁻/ADU with a (0.6, 1.77) bracket (`cv_build_meta`).
*Close:* a PTC from the October flats per mode, or every Poisson-based error explicitly replaced by empirical scatter.

**F8 — MINOR (S2c). The dispersion classifier has a false-positive rate but no false-negative rate.**
Thresholds were tuned on the control sample (`rlmt_diagnostics/dispersion.py:162-206`), then a hold-out was drawn: 19/2215 "dispersed" (0.86%; L band 10/300 = 3.3%) vs control 21/3199. Nobody has established whether those are mislabelled spectra or errors, and there is no labelled truth set for the 448 "direct" + 1,584 "indeterminate" grism-labelled frames. 58,599 plain-label frames were never measured.
*Close:* eyeball-label 200 frames stratified by verdict; publish the confusion matrix with Wilson intervals.

**F9 — MINOR (SN Gate 0).** The verdict string "438 of 632 (371 measured clean + 68 clean by bound)" does not sum (439); R row: 106 + 21 = 127 vs `n_usable` 126. One exclusion rule is unnamed. `project_counts` still carries SN 1052 claimed vs 1054 vs 1129 campaign-unique (M101 alias).

**F10 — MAJOR (S3). Absolute time is unvalidated.**
The EB clock test has a ±4,121–4,165 s envelope (`s3_build_meta`). The "independent" header audit compares two cards written by the same PC clock. Relative timing is fine (drift residuals 1–18 s); any absolute epoch (the 1,071 s edge offset, SN phases, T CrB cross-matching to AAVSO) is not.
*Close:* one night on an EB with a TESS-era ephemeris (σ < 10 s) at re-opening, plus the NTP log.

Satisfied: manifest joins are sound (0 orphans/path mismatches across `frame_dispersion`, `s1_batch`, `frame_times`, `sn_g0_census`, `g_extractions`); global dedup reconciles 330,865 → 198,294; numbers are macro-generated; the paper's handling of aliases, the σ_t injection null and single-cycle epochs is exemplary.

## 2. Plan hardening, per project

**CV_TimeSeries** — not finished.
- ADD CV-R7-band-offset: paired/permutation test + combined-era estimate. *Accept:* abstract sentence matches a test with stated p.
- ADD CV-R7-band-bias-injection: per-band width injection. *Accept:* differential bias ± error tabulated.
- CHANGE CV-P3-oc: per-band constants; night-level N. *Accept:* Ṗ quoted both ways.
- ADD CV-freeze: rebuild from S0 on a fresh manifest, diff `numbers.tex`. *Accept:* zero macro changes or each change explained.

**TCrB_Monitoring**
- CHANGE headline count: of 247 grism-labelled frames, S2c measures 205 dispersed, 35 indeterminate, 7 direct. *Accept:* A0 gate reports N by verdict; "247 spectra" retired unless the 42 are adjudicated.
- ADD TCRB-A0b: slot '6' (12 dispersed/16, 14 nights 2023-05→2024-03) and 'W' (4/8) are measured spectra — a possible pre-2025 Hα baseline. *Accept:* extract-or-reject decision per frame.
- CHANGE TCRB-A7: θ CrB floor must be split-half and cross-checked on a line-free window before any ΔEW is called real. *Accept:* χ²ν of θ CrB continuum-window series in [0.7, 1.4].
- DROP-to-table TCRB-C1 (B: 23 frames on 3 nights; no run > 11 min) and TCRB-C3 (seasons of 23 and 42 d; nothing physical is expected at 1.5–14 d and 47/85 gaps are exactly 1 d). Honest scope: one upper-limit table with the window function; no periodogram figure.
- TCRB-C2 is impossible on disk (zero frames; accrues 2027). Mark as a future-data task outside the paper's critical path.
- Unblock the five S2-"destroyed" tasks (F6), then re-block the two that need October hardware.
- ADD pre-registration: write the EW-change detection rule (threshold, consecutive-night requirement) to the strategy before A5 runs.

**SN2023ixf_LightCurve**
- CHANGE SN-G0b: fix the 438/439 arithmetic; name the rule. *Accept:* clean + bounded = usable per band.
- ADD SN-G0d: S2c-measure all 632 broadband frames (≥90% are "unmeasured"). *Accept:* 0 unmeasured in the photometry set.
- CHANGE SN-S8: nightly means are the unit (5–7 frames/band/night); whole-night bootstrap. Intra-night periodogram to f = 10 d⁻¹ on ~6 points/night should be DROPPED — report injection-defined limits only.
- SN-S6/S1-narrowband: narrowband astrometric corroboration is 183/408 (45%), false-ID 62/80 in slot '1'. Until solved, narrowband epochs are not usable; keep demoted.
- SN-G0c grism: 0 extracted, no wavelength source, no slot-6 flat, timebox long expired (2026-08-21 + 2 wk). *Accept:* close as "not promoted" now; venue stays AJ/PASP.
- SN-S4-ensemble-cal: require held-out check stars (not ensemble members) for the χ²ν ≈ 1 claim.

**BeStar_Grism**
- BE-S-1b before anything else, and extend it: injection–recovery for every star, through detrending. *Accept:* 90% contour per star published before periodograms are opened.
- CHANGE BE-S11 short tier: three stars, ≥3 long nights each — trials factor over 0.3–2 d must be in the FAP (night-block bootstrap, max-statistic). *Accept:* global, not per-frequency, FAP.
- CHANGE BE-S10: the detection rule exists only from 2025-12-05 (standards epoch). Earlier seasons (Phecda, φ Leo, 53 Boo, HD 70340) have no contemporaneous floor. Honest scope: those series are descriptive, no ΔEW claims.
- ADD: 308 indeterminate + 24 direct of 3,883 staged frames need a disposition column before S3.
- BE-S9 era cross-cal: state the degrees of freedom — two transfer stars; an era offset is one number with ~2 constraints.

**DwarfGalaxy_AGN_Survey**
- RE-SCOPE Phase 4 (NGC 5548). Of 143 slot-'6' frames, **91 are measured dispersed (11 nights), 35 indeterminate, 17 direct (3 nights)**. DW-P03's blocker ("zero WCS… clears when the S1b batch reaches these frames", `project_plan.py:1319`) cannot clear: spectra do not plate-solve. DW-P41/42/44/45 (aperture, ensemble, F_var, structure function) are impossible on 3 imaging nights. Honest scope: either a slitless-spectrum flux series on ≤11 nights (no flat, no wavelength source — likely no) or one paragraph. *Accept:* DW-P02-slot6 closed with a per-night verdict table; Phase 4 tasks dropped or rewritten.
- DROP DW-P53/P56 for Dw fields: 5–16 nightly epochs over ≤ 32 d cannot support a period search or eclipse timing. Keep NGC 5238 only (21 nights, 544 frames). DW-P54 completeness map is still worth doing — it is the proof of the above.
- DW-P52 transient search: pre-declare the threshold and count trials (fields × epochs × pixels); require recovery of injected point sources.
- DW-P34 depth: 6 fields have L only, 11–13 frames — expect most to fail the detectability gate; say so in the plan.
- Six L-only fields and Dw1643+07 should be listed as non-contributing unless Cannon's dossier changes that.

**Legacy_Rigel**
- ADD RIG-L0-dedup-reconcile: `legacy_manifest_BAD_collisions.csv` exists beside the archive; the census is meaningless until collisions are resolved. *Accept:* files-on-disk = manifest rows + named exclusions.
- ADD RIG-L1-selection: per target, nights × filters × longest same-filter run, so the go/no-go is made on usable series, not frame counts.
- CHANGE RIG-L2-gonogo: pre-register the criteria (≥1 target with ≥30 nights in one filter with calibration, or overlap extending an RLMT-era baseline) before the census is seen.
- RIG-L1-overlap is the only task with science leverage (ST LMi/T CrB baselines); run it first.

## 3. Cross-cutting — top five

1. **Snapshot the evidence** (F4): hashed product releases, tracked manuscripts, "absent product = test failure" mode.
2. **Make staleness mean something** (F5/F6): data-digest vs page-digest; ledger blockers computed; `show` never prints an unsynced count.
3. **Error models validated per stratum, not globally.** A single χ²ν ≈ 1 hid a 4σ effect and 1.7× over-stated errors (F1/F2). Standing rule: every χ²ν is reported per band/era/mode with its dof, and χ²ν < 0.5 is treated as a defect equal to χ²ν > 2.
4. **Nulls need power, not vote counts.** Every "not detected" must carry the injected effect size it would have recovered, tested with the same estimator, and a combined test across strata. Ban "k of N significant".
5. **Truth sets and absolute references**: a labelled dispersion sample (F8), an absolute clock check (F10), per-mode PTC (F7) — all three are one night of October telescope time plus an afternoon of labelling.

## 4. Verdict

| Project | Verdict |
|---|---|
| CV_TimeSeries | **Execute with amendments** — not releasable; F1 is a BLOCKER, F2–F4 MAJOR. |
| TCrB_Monitoring | **Execute with amendments** — grism spine sound; drop the archival period/flickering analyses to a limits table; fix counts and blockers. |
| SN2023ixf_LightCurve | **Execute with amendments** — broadband +5.4→+50 d only; close the grism timebox. |
| BeStar_Grism | **Execute with amendments** — injection first; no ΔEW claims before the standards epoch. |
| DwarfGalaxy_AGN_Survey | **Re-scope** — NGC 5548 photometric section is unsupported by the frames; Dw-field time-domain tasks dropped. |
| Legacy_Rigel | **Execute Phase L0–L1 only**, with pre-registered go/no-go; no science commitment. |

What would move me to "satisfied" on the portfolio: F1 resolved by test, F4 snapshot in place, and a status page on which at least one stage is FRESH for a reason I can verify.
