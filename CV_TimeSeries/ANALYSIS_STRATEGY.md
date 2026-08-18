# ANALYSIS STRATEGY — Cataclysmic-Variable Time Series (RLMT Archive)

**MACRO Consortium analysis-strategy committee — Chair's synthesis**
**Date:** 2026-08-16 | **Target journal:** ApJ (AASTeX 7 skeleton: `manuscripts/CV_TimeSeries/main.tex`)
**Panel:** Precision Photometry / Detector Science; Time-Series & Statistical Methods; CV Specialist; Literature & Archival Scout.

---

## 1. Executive summary

This is a **four-target paper** — three polars (ST LMi, VV Pup, EU UMa) plus the SU UMa dwarf nova YZ Cnc — delivering **single-night, cycle-resolved, accretion-state-tagged multi-color light curves** that phase-averaged survey folds (ZTF/ASAS-SN/TESS) cannot produce for 90–125 min binaries. QQ Gem (a V≈7.6 Be star) and AG LMi (a detached eclipsing binary) are **cut** from this paper and spun off. The quantitative anchors are: bright-phase timing + accretion-state histories for the polars (ST LMi flagship: **3,157 raw light frames, 39 nights, two three-color seasons split by the 2024-05 instrument seam**), and superhump/orbital analysis of YZ Cnc's dense 2024 season (**1,920 frames**, 8 s cadence, consecutive-night blocks) conditional on AAVSO outburst-state confirmation. All counts here follow the canonical accounting rule of §3 (raw-tree, alias-merged — the referee round purged every tree-doubled number). Two structural facts shape the plan: (i) **~4,600 of ~5,500 polar Sloan-era frames are unsolved with no per-image zero points** — offline re-solving is the pipeline's first bottleneck and gets its own go/no-go (Phase 0.5); (ii) ST LMi's G/R/I and g/r/i seasons **do not overlap in time**, so the headline color analysis is two independent within-era analyses whose morphology is compared — never a stitched 744-d three-color series. All photometry is ensemble-differential tied to ATLAS-REFCAT2; all times are mid-exposure BJD_TDB computed from scratch (header JD is UTC exposure-start — never use it). The core bet: state-resolved cyclotron color light curves + a modern timing arc on 40-year baselines is an ApJ-grade unit; the paper also publishes the full machine-readable light curves as a data product.

---

## 2. Science questions, ranked by publishability

Ranked with honest venue calibration. The paper plans for ApJ by leading with Q1–Q3; Q4–Q5 are supporting material.

| # | Question | Verdict |
|---|---|---|
| **Q1** | **Accretion-state–resolved, orbit-phase-resolved cyclotron color curves of ST LMi, VV Pup, EU UMa — as *two independent within-era analyses per target* (G/R/I era; g/r/i era) whose color-phase *morphology* is compared across eras.** Color modulation (g−r, r−i vs. orbital phase) tracks cyclotron harmonics vs. field strength; states classified quantitatively (mixture model on nightly means) and anchored to the public survey record. No cross-era color zero-point comparison — the systems never overlap in time on ST LMi (r/i end 2025-05-09; the 2025-12→2026-02 tail is g-only). | **ApJ-grade.** The defensible novelty is **single-night, single-cycle, state-tagged quasi-simultaneous three-color coverage** — not phase resolution per se (ZTF sparse epochs fold fine on 90–125 min periods, and ST LMi/VV Pup carry 40 years of dedicated phase-resolved photometry/polarimetry: Cropper 1990; Ferrario et al. 1993 and successors). What no survey or historical campaign provides is per-cycle color curves inside a known accretion state on a known night. If schedulable before submission, **one more g/r/i season on ST LMi (up now through spring) properly repairs the era seam** — recommend to the observing queue. |
| **Q2** | **Bright-phase timing as an accretion-spot longitude tracker** (per-cycle ingress/egress/centroid timings, seasonal O–C vs. decades-baseline literature ephemerides). | **ApJ-grade as a *constraint/diagnostic*, not a discovery.** Committee ruling (see §"Resolved disagreements"): bright-phase boundaries move with spot longitude — they are accretion diagnostics, **not** clean orbital-period probes. We report spot-longitude drift and a period-stability constraint; no planet/Applegate claims. |
| **Q3** | **YZ Cnc 2024 season: superhump period (and possibly dP_sh/dt) within confirmed outburst states**; orbital hump + flickering statistics in quiescence. | **ApJ-grade if the dense runs are in (super)outburst** (AAVSO cross-match decides — highest-leverage first action). Fallback (quiescent runs): orbital modulation + flickering statistics — honest but weaker; would then be a supporting section, not an anchor. No multi-year period-change claims: the span is 110 d. |
| **Q4** | **Long-term accretion-state duty cycles** (high/low-state statistics) for the three polars, RLMT nights embedded in ZTF/ATLAS/ASAS-SN/AAVSO context. | Supporting section. Standing alone this is redundant with surveys (they own the long-term record); it earns its place only as the frame for Q1–Q2. |
| **Q5** | **AN UMa** (polar prototype sibling): 1,279 raw light frames, 14 nights, 2025-01→2026-01 (verified; memo 4's 2,558 double-counted reduced frames). 12 nights span ≥1 orbit but **only ~7 carry all three filters**; 2026 coverage is two nights (one g-only). | **Conditional add, criterion defined per-filter now:** include in the *color* analysis only if **≥8 full-orbit nights with all three filters** survive curation (currently ~7 — likely fails); otherwise include for timing/state history only, or hold for a follow-up (fresh X-ray context: MNRAS 541, 1913, 2025). Whatever passes: apply the pipeline, show the fold, stop — no fifth full analysis chain. Decide after curation, before writing. |
| — | QQ Gem Hα grism series (122 hrg + 41 lrg) | **Cut.** Be-star Hα monitoring — belongs in the grism program paper. |
| — | AG LMi eclipse timing (1,936 frames, 36 nights, 2023–2026) | **Cut.** Detached EB — spin off as an undergraduate EB timing/Wilson–Devinney note (AJ/PASP/JAAVSO scale). |

---

## 3. Data assessment (verified numbers — rebuilt from the canonical SQL, referee round 2026-08-16)

**Canonical accounting rule** (adopted for every number in the paper, and **actually run for every number below** — the previous draft's table quoted tree-doubled counts in violation of its own rule): `tree='rawimage'`, `imagetyp LIKE 'Light%'`, `error IS NULL`, alias-merged `target_best`, photometric filters only (exclude grisms — **all five archive spellings: `hrg`, `lrg`, `HaGrism`, `OGGrism`, `HaG`** — plus `empty`, `W`, `6`). The catalog double-counts across trees (`rawimage`/`reduced`/`macalester`/`iKon`: ST LMi returns 6,423 rows all-trees vs **3,157** rawimage) and fragments names (`mjcMay0x yzcnc` YZ Cnc aliases; `vv pup`, `EU Uma`, `ST Lmi`, `STLMi-z-series`, `ST-LMi-y-series`). Alias map first; then all counts by coordinate cone (5 arcmin on ra_deg/dec_deg) as a cross-check. **Every §3 number must be reproducible by replaying the published SQL** — a referee will run it.

### 3.1 Per-target coverage (alias-merged raw light frames; verified 2026-08-16, **counts corrected 2026-08-18**)

> **Correction (2026-08-18 referee round).** Three rows quoted the ALL-FILTER rawimage light count instead of the photometric count §3's own canonical rule demands: ST LMi 3,157 → **3,150** (3 frames in slot '6', 4 in 'W'), YZ Cnc 1,920 → **1,915** (3 in '6', 2 in 'W'), VV Pup 1,353 → **1,277** / 29 → **28** nights (the 76 `filter='empty'` frames its own note says are excluded). EU UMa and AN UMa were already correct. Every row is now reproduced by `stage_cv_timeseries` and diffed against this table in §4 of `docs/pipeline/s0c_staging.html`.
> A fifth grism spelling, `HaG` (the Andor/iKon label), also had to be added to the exclusion set: 10 dispersed iKon frames of ST LMi and YZ Cnc were staging as photometry, the `_hires` twins of `_lowres` frames the rule already excluded.

| Target | Class, P_orb | Frames | Nights | Full-orbit nights | Baseline | Notes |
|---|---|---:|---:|---:|---|---|
| ST LMi | polar, 113.9 min | **3,150** | **39** | ~34 (fewer per filter) | 2024-01 → 2026-02 (744 d span; **three-color coverage does NOT span it** — see filter-era note below) | Flagship. **7 nights ≥100 frames; 21 nights ≥50.** Anchor night **2025-02-28: 436 frames, 9.4 h (~5 orbits)** — the previously advertised 2024-03-03 "9.8 h night" contains 18 frames. All-filter median intra-night Δt 95 s; **per-filter median Δt 280 s (~24 points/orbit/filter)** — the number that matters for per-band timing and color curves. ~19 consecutive-night pairs; one **289-d gap** (2024-03-26 → 2025-01-09). |
| YZ Cnc | SU UMa DN, 2.086 h | **1,915** | 26 | 17 | 2024-01 → 2024-05 (**110 d — one season**) | Dense blocks: Feb 21–24 (**313/241/130/135** frames, 8 s exp), Mar 1–4, **May 2–3 (189/237 frames, 30 s, Andor iKon; May 1 has only 31; one stray frame May 20)**. Median Δt 66 s. Dense runs cycle R/G/I → **~100 frames/filter/night**. No multi-year baseline. |
| VV Pup | polar, 100.4 min | **1,277** | **28** | 18 | 2024-11-06 → 2026-01-14 (**434 d**) | δ=−19°: airmass ≥1.57 always (median 1.79; some airmass headers are garbage — recompute from coords+time). **Two cameras, cleanly separated in time: Andor iKon 558 frames (41%), 2024-11→2024-12 (8 nights); ASI Mode0 795 frames (59%), 2025-01→2026-01 (21 nights).** Camera and epoch are confounded → colors per-camera only (see §4.13a). **0 of 1,353 frames has zmag** — calibration entirely via field comps. 76 `filter='empty'` frames excluded. |
| EU UMa | polar, 90.1 min | **993** | **32** | 15 | 2025-01 → 2026-06 (510 d) | **96% of frames are 240 s (958/993)** — quote the 4.4% phase smear of the dominant mode, not a "120–300 s mix"; timing floor ≥100 s. Duty-cycle-limited sampling: honest comb at 86400/(240 s + overhead) ≈ 345 c/d. |
| AN UMa | polar, 114.8 min | 1,279 | 14 | 12 (≥1 orbit); **~7 with all three filters** | 2025-01 → 2026-01 | Conditional (Q5, per-filter criterion). 2026 coverage: two nights, one g-only. |

**ST LMi filter-era structure (drives the whole Q1 design):** Johnson-ish **G/R/I season: 2024 Jan–Mar (~20 nights)**; Sloan **g/r/i season: 2025 Jan–May (~13 nights)**; then a sparse **g-only tail: 2025-12-12 → 2026-02-09 (4 epochs, 23–30 frames each)**. r and i coverage **ends 2025-05-09**; 2025-02-22 is i-only and Feb 23–25 lack i. The two three-color systems have **zero overlap in time** on this target — colors cannot be compared across the seam (we refuse to color-transform the CV, §4.13), so Q1 is two within-era analyses compared by morphology.

### 3.2 Instruments — this is a three-camera, two-filter-system campaign

| Era | Camera / readout | Key properties | **Saturation ceiling** |
|---|---|---|---|
| MaxIm (2023-02 → ~2024-05) | DL Imaging GSENSE4040 CMOS, 4096², 9 µm, 0.54″/px; `High Gain` | ZP ≈ 21.7 | **~3,500 ADU** (tiny high-gain full well; veto peak > 3,200) |
| MaxIm | same, `StackPro` (~16 on-camera stacked reads — **N_sub unverified**) | non-Poisson correlated noise; candidate model variance = flux/g + N_sub·σ_read² — **to be measured, not assumed** (§5 rows 1–2); no mixed-mode fits before the ladder exists | ~55–56 k ADU, **but sub-read clipping at 3.5k×N_sub is non-obvious — verify** |
| pyscope (~2024-05 →) | ASI CMOS, `Mode0`, 2×2 on-camera bin, EGAIN 0.247 e⁻/ADU | targets in filename (`target_best` handles); carries **59% of VV Pup (795 frames, 2025-01→2026-01)** | 65,535 ADU; veto peak > 55 k (binned sums saturate on one hot sub-pixel) |
| — | Andor iKon CCD/EMCCD 2048², `1MHz High Sensitivity 16-bit`, −65 °C | carries **41% of VV Pup (558 frames, 2024-11→2024-12)** and **YZ Cnc's May 2024 block (~458 frames, 30 s)** + macalester/iKon tree | standard CCD; verify full well |

Filters: Johnson-ish G/R/I (MaxIm) vs Sloan g/r/i (pyscope). The 2024-05 seam sits **mid-baseline for ST LMi** — any multi-year series stitches two photometric systems and up to three detectors.

### 3.3 Quality metrics (canonical rule, science targets — the previous draft's global tree-mixed statistics are retired)

- FWHM median 4.93 px (~2.7″ on DL); p10–p90 3.0–7.4 px. Airmass median ~1.27 (except VV Pup).
- **Plate solves — the pipeline's first bottleneck, not a cleanup task.** Verified solved fractions for the science targets: ST LMi Sloan era **5%** (99/1,904), VV Pup **21%** (286/1,353), EU UMa **27%** (270/993), AN UMa **8%** (105/1,279); only YZ Cnc (MaxIm era) is healthy at **94%** (1,806/1,920). **~4,600 of ~5,500 polar Sloan-era frames are unsolved** → without offline re-solves there is no ensemble photometry for them. Cause of failure (2×2-binned ASI, short exposures?) is unknown — hence the Phase 0.5 go/no-go experiment (§4).
- **`zmag` is ZERO for the polars' Sloan-era rawimage frames** (VV Pup 0/1,353; EU UMa rawimage 0; ST LMi Sloan-era ~0). The zmag cloud veto (§4.11) therefore **exists only for MaxIm-era data (YZ Cnc, ST LMi 2024)**; for essentially all polar data the ensemble-flux-ratio veto is the *primary* method, not a fallback. Where zmag exists, within-night R-band zmag RMS median **0.147 mag** → per-image ZMAG is a quality veto, never a calibration. The paper states honestly: per-image ZPs are absent for the Sloan era; calibration is wholly differential + nightly REFCAT2 tie. (Any §3.3-style global statistics quoted in the paper must be re-derived under the canonical rule — the old "9,121/17,111" figure mixed trees and targets.)
- Timing headers (verified on `rawimage/2026-02-09/mem_ST_LMi_g_300s_...fts.fz`): header `JD` = **UTC exposure start** (catalog `jd` verified = UTC exposure start to <2 s); `JD-HELIO` ≈ mid-exposure heliocentric UTC (598 s offset on a 300 s frame — consistent but unverified across eras). **Use neither.**

---

## 4. Analysis method — step by step

### Phase 0 — Curation (blocks everything else)

1. Build `targets_cv.sql`: one canonical view per target — regex alias map on `target_best`, `tree='rawimage'`, `imagetyp LIKE 'Light%'`, `error IS NULL`, `jd IS NOT NULL`, photometric filters only. Tag each row with `era` (MaxIm/pyscope), `camera`, `readout`, `filter_system`. Cross-check counts with a 5-arcmin coordinate cone. **Publish this SQL in the paper's repo** (referee-proofing).
2. **AAVSO cross-match for YZ Cnc 2024-02-21 → 2024-05-03 immediately** — it decides whether the paper contains a superhump period derivative (Q3 branch point). Do this before pipeline work.
3. Pull survey context: ZTF forced photometry, ATLAS forced-photometry server, ASAS-SN Sky Patrol v2, AAVSO AID for all targets; TESScut for outburst epochs (treat 21″-pixel blending explicitly); Gaia DR3 parallaxes; eROSITA-DE DR1 X-ray states.

### Phase 0.5 — Astrometry go/no-go (gates the polar pipeline; do NOT skip to Phase 2)

3a. **Re-solve experiment first:** draw a random 200-frame sample per target from the unsolved pool (ST LMi Sloan, VV Pup, EU UMa, AN UMa), run offline astrometry.net with local indices, and **report the per-target success rate** before committing to any solved-fraction acceptance criterion. The failure cause on these frames (2×2-binned ASI, short exposures, field density?) is unknown — the §5 ">95% solved" criterion is an aspiration until this experiment supports it.
3b. **Go/no-go:** if the sample success rate is ≥95%, batch-solve everything and proceed. If it lands at 70–95%, quantify which nights/phases are lost and re-check per-filter full-orbit night counts before finalizing the target list. If <70%, escalate: the polar ensemble photometry plan needs redesign (e.g., solve one frame per pointing and propagate WCS by offset tracking) before any Phase 2 work.

### Phase 1 — Timing foundation

4. For every frame: mid-exposure = `DATE-OBS + EXPTIME/2` → **BJD_TDB** via `astropy.time.Time` + `light_travel_time` with Winer's `EarthLocation` and a JPL ephemeris (Eastman, Siverd & Gaudi 2010). The barycentric term is ±8 min/season ≈ 7% of an orbital cycle; UTC→TDB alone is ~69 s. Non-negotiable.
5. **Era audit:** verify the DATE-OBS convention (start vs mid) independently in MaxIm and pyscope frames — the start-vs-mid ambiguity is 150 s for 300 s exposures, larger than any statistical timing error. Document as a manuscript appendix.
6. **Clock validation:** time a well-determined eclipsing system observed in the same seasons and recover its published ephemeris — demonstrates second-level absolute clock authority rather than asserting it. Candidate: AG LMi (36 nights in-archive) — **but a 1.359-d detached EB with narrow eclipses may have few or zero in-eclipse nights among 36 sparse ones. First count in-eclipse nights against the corrected VSX period (1.3590176 d); if <3 usable minima, fall back to another archived TESS-era EB** before writing this into §5.

### Phase 2 — Photometry (shared pipeline with TCrB — see §9)

7. **Aperture photometry, `photutils`** (AstroImageJ acceptable for student verification runs). No PSF fitting: at 2.7″ seeing, ~0.5″/px, uncrowded fields, it adds only model-mismatch systematics on CMOS with per-pixel gain quirks. Aperture = 1.5×FWHM per frame (header `fwhm`), sky annulus 4–6×FWHM sigma-clipped; apertures fixed in **arcsec** per camera, re-derived per night — never one global pixel radius across three plate scales.
8. **Calibration frames:** nightly/nearest twilight flats per (camera, filter, binning) from `calib/`; fringe-frame subtraction for i (and any z) with fitted per-image fringe amplitude (scale, don't just divide), frames from `fringes/`.
9. **Ensemble differential photometry** (Honeycutt 1992 inhomogeneous ensemble): 8–15 comps per field from **ATLAS-REFCAT2** (12.5 < r < 16.5, 0.3 < g−r < 1.1, Gaia RUWE < 1.4, no neighbor within 3×FWHM contributing >1%, clean variability flags); inverse-variance combine; iteratively reject comps with residual RMS > 1.5× the noise model. Tie the ensemble zero point to REFCAT2 g,r,i per night → PS1 AB system to ~0.01–0.02 mag absolute.
10. **Detector vetoes per mode:** High Gain peak > 3,200 ADU; Mode0 peak > 55 k; StackPro > 45 k until the linearity ladder (see §5) says otherwise; propagate per-frame EGAIN. **Check every YZ Cnc outburst-night frame for saturation — the February 8 s High Gain frames (V≈10.5 clips even at 8 s and looks fine at 5% of the 16-bit range) AND the May 2–3 30 s frames (verified: Andor iKon `1MHz High Sensitivity`, 457 frames — outburst V≈10.5 at 30 s needs its own CCD full-well audit).**
11. **Frame-quality veto:** the **primary** cloud veto is the ensemble-flux-ratio test (drop frames whose comp-ensemble flux falls >0.5 mag below the nightly median) — it must work standalone because **zmag does not exist for essentially all polar Sloan-era frames** (§3.3). Where zmag exists (MaxIm era: YZ Cnc, ST LMi 2024), apply the zmag veto as a secondary cross-check.
12. **Airmass/color terms:** fit second-order extinction k″·(g−r)·X inside the ensemble solve. Mandatory for VV Pup (X ≥ 1.57 always; a blue CV vs red comps at X=1.8 in g is a 10–20 mmag airmass-locked systematic that masquerades as orbital-phase structure) and for all g-band work. **For VV Pup the term is solved per-camera on the ~558-frame iKon and ~795-frame Mode0 subsets separately** — the cameras never overlap in time, so a joint solution would alias camera into extinction.
13. **Cross-era discipline (supersedes "stitching"):** treat each (camera, filter) pair as its own photometric system; solve linear color transformations G/R/I → g/r/i from the ensemble stars (overdetermined) and tabulate coefficients ±σ **as metadata for the data release only**. **Never color-transform the CV itself** — CV colors are nonstellar and state-dependent. Consequence, stated plainly: **CV colors are not comparable across the 2024-05 seam at all.** All color-phase science is within-era; cross-era comparison is of *morphology* (phase of color extrema, amplitude ratios, state dependence), never zero points. Timing analyses run per contiguous instrument block first, then jointly, showing the answer doesn't move.
13a. **VV Pup camera handling (verified: cameras separate cleanly in time — iKon 2024-11→12, Mode0 2025-01→2026-01; zero interleaving):** camera and epoch are fully confounded, exactly like the ST LMi era seam. **Demote VV Pup to per-camera folded color curves; its cross-camera role is timing and state history only.** Do not build a joint VV Pup color-phase figure.
14. **Faint-phase handling:** forced photometry at the solved position; report upper limits for non-detections (polar low states at V ≳ 18) so state statistics aren't censored.

### Phase 3 — Time-series analysis

15. **Period verification** (all periods are known — this is confirmation + alias hygiene, not discovery): astropy `LombScargle` (floating mean, 10× oversampled grid, 2–200 c/d) **plus** PDM (Stellingwerf) **plus** conditional entropy, per filter per era; polars are strongly non-sinusoidal, so require all three methods to agree on the same alias family; cross-check multi-harmonic LS (`nterms=3`). **Single-night periodograms first** (ST LMi's **21 nights ≥50 frames — only 7 exceed 100** — give an alias-free period to 1–2%, anchored by 2025-02-28's 436-frame, 5-orbit night) to select the correct peak in the multi-season comb — not the other way around. Compute and publish the spectral window per season (window is era-dependent) and the chosen-peak/strongest-alias power ratio.
16. **Polars — bright-phase timing, per band** (ST LMi, VV Pup; EU UMa with smearing term): fold on published ephemerides; fit each cycle's bright phase with a sigmoid-edge model, uncertainties via `emcee`, cross-checked with Kwee–van Woerden on the clean egress branch. **Timing is done per filter and band-dependent egress epochs are reported as a cyclotron-geometry result** — combining filters for timing is unsafe because cyclotron beaming makes bright-phase boundaries color-dependent (the harmonics dominating g vs i beam differently). The per-filter cadence is **280 s median (~24 points/orbit/filter)** on ST LMi, so per-cycle σ_t < 60 s is **not assumed: demonstrate achievable σ_t by injection on the 2025-02-28 night (436 frames) BEFORE adopting any σ_t threshold**; if realistic per-cycle uncertainties are 100–200 s, drop the per-cycle O–C tier and carry seasonal-mean timings only, re-deriving the |dP/dt| envelope accordingly. (EU UMa: propagate the ≥100 s phase-smearing floor from 240 s exposures; include only if still informative.) One timing per night minimum via cross-correlation against a mean-cycle template; per-night bootstrap (resample **nights**, not points).
17. **O–C construction:** local epoch (T0 from our own data) + **literature period** for the cycle count; explicit cycle-count ambiguity analysis across gaps (ST LMi's **289-d gap** (2024-03-26 → 2025-01-09) needs δP < P²/2Δt ≈ **1.1×10⁻⁵ d** ≈ 1 s — only a decades-baseline literature ephemeris achieves this). Seasonal mean O–C vs. the 40-yr baseline → spot-longitude drift + |dP/dt| constraint. State detectable |dP/dt| up front (re-derived from the demonstrated σ_t, per step 16); report a constraint unless the data demand more.
18. **Accretion states:** nightly mean calibrated magnitude per filter; two-component mixture model on the nightly-mean distribution — state boundary from the mixture, not by eye. RLMT nights over-plotted on the ZTF/ATLAS/ASAS-SN/AAVSO record (Duffy et al. 2022 state framework).
19. **YZ Cnc:** per AAVSO state tags — if dense runs are in (super)outburst: per-run low-order polynomial detrend, LS + PDM for P_sh, Kato-style O–C of superhump maxima for P_sh and dP_sh/dt, resting the claim only on the consecutive-night blocks (Feb 21–24, Mar 1–4, **May 2–3**; the ~1.4 d P_sh–P_orb beat requires them). The dense runs cycle R/G/I at ~100 frames/filter/night — adequate for superhump-maximum timing, but **build the Kato-style O–C per filter, or include an explicit color-amplitude term** (superhump amplitude is color-dependent). If quiescent: orbital hump + flickering statistics, stated as the fallback it is — **but first run an empirical S/N check on the 8 s High Gain frames at quiescent V≈14.5** (they may be sky/read-noise dominated on a 0.5 m); the fallback is only promised if the check passes.
20. **Detrending discipline:** per-night systematics fit **jointly** with the periodic model (low-order airmass polynomial, or `celerite2` Matérn-3/2 GP with timescale prior bounded below at 3× the candidate period). Never pre-whiten a night spanning <3 cycles with a free smooth trend (EU UMa nights average ~1.5 cycles — a free spline eats the orbit). Show with/without comparisons.
21. **Detection limits + injection–recovery:** quote A_min = σ·sqrt(4z/N) with measured σ (check stars), z ≈ 18 (M≈1e5, FAP 0.1%), **computed from the true canonical-rule N per target and per filter** (the previous draft's limits were √2 optimistic from doubled counts). Indicative all-filter values with real N: **ST LMi (N=3,157) ~4–5 mmag; YZ Cnc (N=1,920) ~5–6; VV Pup (N=1,353) ~6–7; EU UMa (N=993) ~7–8; single 200-frame night ~12; per-filter subsets worse by ~√3.** Final numbers come from measured σ, not these placeholders. Verify with injection–recovery at the real timestamps (100 trials per amplitude–period cell; 90%-recovery contour figure). Use bootstrap FAPs, not analytic, wherever StackPro/red-noise data enter.

---

## 5. Calibration & validation plan

| Item | Method | Acceptance criterion |
|---|---|---|
| Linearity, per readout mode | Twilight-flat exposure ladder for High Gain, StackPro, Mode0 (+ iKon full-well check); verify N_sub ≈ 16 and the 3.5k/56k/65k ceilings; sub-read clipping test for StackPro. **Until the ladder exists, both N_sub and the StackPro variance model (flux/g + N_sub·σ_read²) are unverified hypotheses — no mixed-mode fits before then** | Measured linearity curves in the paper; veto thresholds and noise model set from data |
| Noise model, per mode | Empirical: check-star scatter binned by magnitude, per readout mode, per night; per-frame EGAIN propagated | RMS-vs-mag plots with noise model overplotted, per camera per band (**non-negotiable figure**) |
| Per-point errors | Pipeline errors validated against check-star RMS; fitted per-night error-scale inflation | Scale factors ~1.0–1.5; outlier nights flagged |
| Absolute calibration | Nightly REFCAT2 tie | ~0.01–0.02 mag on PS1 AB; nightly ZP scatter tabulated |
| Comp-ensemble stability | Comp residual light curves across the full campaign | Night-to-night RMS < 0.01 mag; long-term constancy across the 2024-05 era seam |
| Flat-field stability | Ensemble residuals vs (x, y) drift | If correlated at > 3 mmag → low-order 2-D illumination correction from dithered ensemble residuals |
| Cross-era consistency | Transformation coefficients ±σ; per-block vs joint timing solutions | Timing answers agree within errors across the seam |
| Clock authority | Count AG LMi in-eclipse nights first (corrected VSX period 1.3590176 d); if ≥3 usable minima, recover its published ephemeris; else substitute another archived TESS-era EB | Residual < 2 s systematic; audit appendix written; validator choice documented |
| Astrometry | **Phase 0.5 go/no-go first**: 200-frame random-sample re-solve per target, success rate reported. Then offline astrometry.net re-solves for the **~4,600 unsolved polar Sloan-era frames** (ST LMi Sloan 95% unsolved, VV Pup 79%, EU UMa 73%, AN UMa 92%) | Acceptance threshold set FROM the Phase 0.5 experiment (target ≥95% solved, but not asserted until demonstrated) |
| Expected precision (state in paper) | 8–10 mmag/frame at r ≤ 15; 20–25 mmag at r = 16.5; 0.05–0.10 mag at r = 18–19 | Demonstrated by the RMS-vs-mag figures |

---

## 6. Failure modes & mitigations

1. **Silent High Gain saturation at ~3,500 ADU** (looks fine at 5% of 16-bit range) → peak-pixel veto at 3,200 ADU; audit every YZ Cnc outburst frame **in both the February High Gain runs and the May Andor iKon runs**; report veto counts per target.
2. **Catalog double-counting / name fragmentation** → canonical curation SQL (tree-pinned, alias-mapped), coordinate-cone cross-check, SQL published. **This failure mode already bit this document once** — the first draft's §3 quoted tree-doubled counts; every number is now re-derived by actually running the rule, and no number enters the manuscript without a replayable query.
3. **Camera/era-seam steps masquerading as astrophysics** (2024-05 MaxIm→pyscope seam; VV Pup's iKon→Mode0 hand-off — both are epoch-confounded) → per-system light curves, ensemble-star long-term constancy test, per-block-then-joint timing, **within-era-only color science with morphology-level cross-era comparison** (§4.13, 4.13a).
4. **CV color-transformation abuse** → system-native magnitudes for the target, transformations as metadata only; no cross-seam color comparison at all.
5. **VV Pup airmass floor (X ≥ 1.57) + two epoch-confounded cameras with no zmag and 79% unsolved astrometry** → per-camera second-order color-extinction solves on half-sized subsets; constant-star validation figure; colors per-camera only, cross-camera role = timing/states; sanitize the garbage airmass headers (mean 13.7 — garbage) by recomputing airmass from coordinates + time.
5a. **Astrometry re-solve failure** (cause of the Sloan-era failures unknown) → Phase 0.5 sample experiment before any pipeline commitment; escalation path defined (WCS propagation by offset tracking) if astrometry.net does not converge.
6. **Wrong alias family** → single-night periodograms first, three independent methods, window functions per season, alias-ratio table.
7. **Timing sloppiness** (UTC start-of-exposure JD; unverified JD-HELIO; StackPro mid-time semantics for stacked reads) → own BJD_TDB recipe, era audit appendix, EB clock validation.
8. **Detrending eats the signal on short nights** → joint GP/polynomial + periodic fits only; GP timescale prior ≥ 3× period; with/without figures.
9. **Overclaiming** — the three hard lines: no YZ Cnc multi-year period change (110-d span); no orbital-period-change *discovery* from bright-phase timings (they track spot longitude); no planet/Applegate interpretation. Constraints, not discoveries.
10. **Superhump fallback** — if AAVSO says the dense runs were quiescent, Q3 degrades to orbital hump + flickering; the plan (and abstract) accommodate this now, not in the referee response — **contingent on the empirical S/N check of §4.19** (8 s High Gain at quiescent V≈14.5 may be sky/read-noise dominated; verify before promising the fallback).
10a. **Per-cycle timing precision shortfall** — if the 2025-02-28 injection test (§4.16) shows per-cycle σ_t of 100–200 s at 280 s per-filter cadence, the per-cycle O–C tier is dropped, seasonal means carry the timing section, and the |dP/dt| envelope is re-derived; band-dependent egress epochs become a result, not a nuisance.
11. **Censored state statistics** from faint-phase non-detections → forced photometry + upper limits.
12. **StackPro correlated noise breaks analytic FAPs/χ²** → empirical errors, bootstrap FAPs, mode-segregated fits (no mixed-mode timing fits without a cross-calibration night).

---

## 7. Figure list

1. **Coverage/cadence map** — all four targets: nights vs date, color = filter system, symbol = camera; survey epochs (ZTF/ASAS-SN) underlaid to show what RLMT adds: single-night, cycle-resolved sampling that sparse survey epochs cannot provide.
2. **RMS vs magnitude** — ensemble/check stars per camera and band, noise model overplotted; demonstrates the 8–25 mmag precision claims.
3. **Linearity & saturation** — twilight-ladder curves for High Gain / StackPro / Mode0 with adopted veto thresholds.
4. **Spectral windows + periodograms** — per target: LS/PDM/CE, window function per season, alias-ratio annotation; inset single-night periodogram for ST LMi.
5. **ST LMi phase-resolved light curves** — folded by accretion state (high/low), **as separate per-era panels (G/R/I 2024 | g/r/i 2025)**; bright-phase model fits overlaid; the g-only 2025-12→2026-02 tail shown in the state history (Fig 8), not here.
6. **Cyclotron color-phase diagram** — G−R/R−I (2024 era) and g−r/r−i (2025 era) vs orbital phase by state, **per-era panels with morphology comparison annotated — no cross-era color axis**; ST LMi + VV Pup **per-camera** (the headline physics figure).
7. **VV Pup + EU UMa folded light curves** — same layout as Fig. 5 (VV Pup split iKon | Mode0), EU UMa with smearing kernel noted.
8. **Long-term state histories** — nightly means over ZTF/ATLAS/ASAS-SN/AAVSO curves, mixture-model state boundaries marked; upper limits shown.
9. **O–C diagrams** — ST LMi (+VV Pup, EU UMa) seasonal **per-band** bright-phase timings vs literature ephemerides, historical timings included; detectable-|dP/dt| envelope drawn from the demonstrated (injection-tested) σ_t; band-dependent egress offsets shown as a panel, not averaged away.
10. **YZ Cnc season overview** — full 110-d curve with AAVSO/ASAS-SN outburst states shaded; dense-run insets.
11. **YZ Cnc superhump analysis** (conditional) — per-run periodograms + Kato-style O–C of superhump maxima; else orbital-hump/flickering statistics figure.
12. **Injection–recovery** — 90%-recovery amplitude–period contours per target at the real timestamps.
13. *(Appendix)* **Timing audit** — DATE-OBS/JD/JD-HELIO relationships per era; EB clock-validation residuals.

---

## 8. Manuscript outline (mapped to `manuscripts/CV_TimeSeries/main.tex`)

Current skeleton sections → planned content. Retitle to: *"Orbit-Resolved Multi-Color Photometry of Three Polars and the Dwarf Nova YZ Cnc from the Robert L. Mutel Telescope"*.

- **§1 Introduction** (`sec:intro`) — polars as cyclotron labs; **novelty framed precisely: single-night, cycle-resolved, state-tagged multi-color coverage vs. phase-averaged survey folds — NOT "surveys cannot phase-resolve these orbits" (they can; ZTF polar orbital light curves are published), and with an explicit comparison to the 40-year phase-resolved photometry/polarimetry literature on ST LMi and VV Pup (Cropper, Ferrario, and successors), stating what those campaigns lacked (state-tagging against the modern survey record; per-cycle quasi-simultaneous colors)**; target census incl. TESS CV period catalogues (arXiv:2603.03539, 2607.08727); scope statement (4 targets; AN UMa if Q5 passes per-filter).
- **§2 Observations and Data Reduction** (`sec:obs`) — RLMT + Winer; **instrument table** (era, camera, readout, filter, EGAIN, ceiling ADU, plate scale); curation rules + counts (Table: per-target coverage); calibration frames; pipeline summary (points to shared pipeline, §9); Figs 1–3.
  - **§2.x Timing** — BJD_TDB recipe + era audit (Fig 13 → appendix).
- **§3 Analysis** (`sec:analysis`) — ensemble photometry + REFCAT2 tie; cross-era transformations; period verification + windows (Fig 4); bright-phase timing model; state classification; detection limits + injection–recovery (Fig 12).
- **§4 Results** (`sec:results`) — per target: ST LMi (Figs 5–6, 8, 9), VV Pup + EU UMa (Figs 6–9), YZ Cnc (Figs 10–11); machine-readable light-curve tables.
- **§5 Discussion** (`sec:discussion`) — cyclotron color interpretation vs Ferrario et al. 1993 / Cropper 1990; spot-longitude drift + period-stability constraints vs the WD-binary timing literature (arXiv:2602.17800; OY Car ApJ study); YZ Cnc vs Kato superhump surveys; states vs Duffy et al. 2022; eROSITA context.
- **§6 Conclusions** (`sec:conclusions`) — deliverables incl. Zenodo data release.
- **Appendices** — A: timing audit; B: transformation coefficients; C: curation SQL + comp-star lists.
- `references.bib` — seed with: Cropper 1990; Ferrario et al. 1993; Howell et al. 1995; Honeycutt 1992; Eastman et al. 2010; Duffy et al. 2022; MNRAS 541, 1913 (2025, AN UMa); Kato superhump series; TESS CV catalogues (2026); Stellingwerf 1978; VanderPlas 2018.

---

## 9. Pipeline & dependency notes

- **Shared photometric pipeline with TCrB_Monitoring** (per project README): ensemble differential photometry core (photutils + REFCAT2 + Honeycutt solve) is common; CV-specific layers (per-mode saturation vetoes, StackPro noise model, bright-phase timing, state mixture model) live in this project. Keep the core in one repo imported by both — divergence is a maintenance failure mode.
- **Environment:** conda env `rlmt-checks` (`/opt/miniconda3/envs/rlmt-checks/bin/python`, astropy 6.1.7). Add: photutils, celerite2, emcee, astrometry.net + local indices (needed for the **~4,600 unsolved polar Sloan-era frames — 73–95% per target, the pipeline's first bottleneck**).
- **Data:** archive at `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/`; catalog `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite` (always quote the path — spaces). Curation SQL and derived light curves version-controlled with the manuscript.
- **External pulls:** ZTF forced photometry, ATLAS forced photometry, ASAS-SN Sky Patrol v2, AAVSO AID (YZ Cnc first), TESScut, Gaia DR3, eROSITA-DE DR1, ATLAS-REFCAT2 local copy for comp selection.
- **Execution order for a grad student:** Phase 0 (curation + AAVSO YZ Cnc check) → **Phase 0.5 (astrometry re-solve go/no-go — 200-frame samples, before anything else commits)** → Phase 1 (timing foundation; AG LMi in-eclipse count for the clock validator) → linearity ladders + noise model (§5 rows 1–2) → **σ_t injection test on ST LMi 2025-02-28 (sets the timing tier)** → Phase 2 photometry on ST LMi first (richest, validates pipeline) → VV Pup (hardest calibration, per-camera) → EU UMa → YZ Cnc → Phase 3 analysis → AN UMa go/no-go (per-filter criterion) → figures → draft. **In parallel: request one more g/r/i ST LMi season from the queue — it converts the era-seam workaround into a genuine repair if it lands before submission.**

### Planned observations (October 2026 restart) — 2026-08-17 update

Winer is offline for monsoon season until **October 2026**; the ST LMi request above is already submitted (`ops/2026-08_observatory_request.md`, Item C) to execute at the re-opening, first in queue. What this project expects from October:

- **ST LMi g/r/i time-series blocks, October 2026 → spring 2027** — the extra Sloan season that repairs the 2024-05 era seam (Q1). The "wait for it?" decision stays at the February consortium call, per §2 Q1.
- **Era-resolved calibration frames** for the polar-era readout modes (the S0b `calib_gaps` shopping list in the observatory request, Item B — Mode0 biases/flats rank at its top with CV among the blocked projects), acquired at re-opening.
- Nothing else: all archival analysis (Phases 0–3) proceeds through the monsoon closure on existing frames.
- **Pipeline ingest:** October nights are first-class planned inputs — after each post-October archive sync, S0 then S0b re-run (idempotent) and the manifest diff feeds the canonical counts; no new-frame number enters this paper except through that loop.

---

## Resolved disagreements (chair's calls)

1. **QQ Gem & AG LMi: cut (Memos 3–4) vs keep as monitoring (Memos 1–2).** *Cut both.* Memo 4's identifications are decisive (QQ Gem = HD 46264, V≈7.6 Be star; AG LMi = EA detached binary): neither is a CV, and keeping them hands a referee a credibility hit that taints the good targets. QQ Gem's grism series goes to the Be/grism program; AG LMi becomes a separate EB note — and its nights get reused here as the clock-validation source (§5), *if* the in-eclipse count passes (§4.6).
2. **AG LMi period discrepancy** (Memo 3: 0.6795 d; Memo 4: VSX 1.3590176 d — exactly 2×): moot for this paper, but flag to the EB spin-off: Memo 3 quoted the half-period alias; adopt VSX and re-verify.
3. **Frame counts.** Memos disagree (ST LMi 3,102 / 6,217 / 6,423) because of tree duplication and alias coverage. *Canonical rule adopted:* `tree='rawimage'` + alias-merged names + light frames + `error IS NULL`; paper publishes the SQL. **Referee-round postscript: the chair re-verified only AN UMa, and the first draft's §3.1 quoted the tree-doubled numbers for every other target anyway. The rule is now enforced by execution, not adoption — every §3 number was re-derived by running the SQL on 2026-08-16 (ST LMi 3,157/39; YZ Cnc 1,920/26; VV Pup 1,353/29; EU UMa 993/32; AN UMa 1,279/14).**
4. **O–C anchoring: literature ephemeris (Memo 2) vs local epochs + drift warning (Memo 4).** *Both, explicitly:* local T0 from our data, literature *period* for cycle counting (required by the 289-d gap arithmetic), plus an explicit cycle-count-ambiguity analysis; if extrapolated phases are off by >0.1 cycle, re-derive the local ephemeris before touching historical O–C.
5. **What bright-phase O–C means: period probe (Memos 2–3) vs accretion diagnostic (Memo 4).** *Memo 4 wins on physics:* bright-phase boundaries track spot longitude in self-eclipsing polars. Timing results are framed as spot-longitude drift + a period-*stability constraint*; no orbital-period-change discoveries, no planet/Applegate language.
6. **QQ Gem grism as bonus figure (Memo 2) vs full removal (Memos 3–4).** *Full removal* — a Be-star figure in a CV paper invites the scope question we just eliminated.
7. **AN UMa addition (Memo 4 only).** *Conditional yes* (Q5): same pipeline, include only if ≥8 full-orbit nights survive curation; decision before writing, not during review.
8. **QQ Gem "long-term state statistics" (Memos 1–2 fallback).** Superseded by the cut — no CV state statistics exist for a Be star.

**Bottom line for the team:** ST LMi's ~34 full-orbit nights (fewer per filter) and YZ Cnc's 8 s dense runs carry this paper; VV Pup and EU UMa make it an ensemble. The two highest-leverage next actions are the **AAVSO cross-match on YZ Cnc's February–May 2024 nights** (decides Q3) and the **Phase 0.5 astrometry re-solve experiment** (decides whether the polar pipeline runs at all as designed). At ~7,400 canonical polar frames plus a 1,920-frame dwarf-nova season, the real dataset — half the size the first draft believed — still amply supports the paper.

---

## Internal Referee Round (2026-08-16) — concerns and resolutions

All referee numbers independently re-verified against the catalog by the panel before adoption (counts, solve fractions, camera splits, filter end-dates, dense-run totals confirmed; one referee detail corrected: YZ Cnc's May frames are Andor iKon, not Mode0).

| Concern | Resolution |
|---|---|
| **M1** Table 3.1 = tree-doubled counts | §3 fully rebuilt from executed canonical SQL: 3,157/1,920/1,353/993/1,279; nights ≥100 = 7 (not 21); dense runs 313/241/130/135; ~19 pairs; detection limits recomputed √2 worse from true N (§4.21). |
| **M2** Sloan astrometry/ZP crisis understated | §3.3 restated: 73–95% unsolved per polar, zmag = 0 for polar Sloan frames; new Phase 0.5 go/no-go re-solve experiment gates the pipeline; §5 acceptance threshold now set by the experiment; ensemble-flux-ratio veto promoted to primary (§4.11). |
| **M3** Three-color claim breaks at the era seam | Q1 reframed as two independent within-era analyses, morphology-only cross-era comparison; Figs 5–6 redrawn per-era; "744-d three-color series" language removed; one more g/r/i season requested from the queue. |
| **M4** Per-cycle σ_t < 60 s unproven at 280 s/filter cadence | Timing now per band (band-dependent egress = cyclotron result); σ_t demonstrated by injection on 2025-02-28 (436 frames) before any threshold adopted; fallback tier (seasonal means only) pre-defined (§4.16, §6.10a); night bookkeeping fixed (2024-03-03 has 18 frames). |
| **M5** Survey-gap novelty overclaimed | Novelty rewritten everywhere as "single-night, cycle-resolved, state-tagged multi-color" vs. phase-averaged folds; explicit comparison to 40-yr phase-resolved literature added to §8 intro plan. |
| **M6** VV Pup camera attribution wrong; calibration confounded | Corrected: 59% Mode0 / 41% iKon, cleanly time-separated (verified — no interleaving) → colors per-camera only, cross-camera role = timing/states (§4.13a); per-camera extinction solves (§4.12); 0 zmag, 76 empty, 434 d, airmass headers garbage — all corrected. |
| **m1** Residual numerics | 289-d gap → δP ≲ 1.1×10⁻⁵ d (§4.17); 26 nights; ~34 full-orbit nights. |
| **m2** EU UMa exposure mix | Restated: 96% at 240 s, 4.4% smear, comb ≈ 345 c/d; zmag figure retired (rawimage = 0). |
| **m3** YZ Cnc May block + per-filter density | May 2–3 (May 1 = 31 frames); ~100 frames/filter/night stated; Kato O–C per filter or with color-amplitude term (§4.19). |
| **m4** Quiescent-fallback S/N unproven | Empirical S/N check on 8 s High Gain at V≈14.5 required before the fallback is promised (§4.19, §6.10). |
| **m5** Saturation audit scope | Extended to the May Andor iKon 30 s outburst frames (§4.10, §6.1) — referee's "Mode0" corrected to iKon per catalog. |
| **m6** AG LMi validator may have no eclipses | In-eclipse night count (corrected VSX period) required first; backup TESS-era EB validator named in §4.6/§5. |
| **m7** StackPro noise model asserted | §3.2 now marks N_sub and the variance model as to-be-measured; §5 row 1 forbids mixed-mode fits until the ladder exists. |
| **m8** AN UMa criterion filter-blind | Q5 criterion redefined per-filter now (≥8 full-orbit three-filter nights; currently ~7 — likely timing/states only). |
| **m9** §3.3 global stats mixed trees | Retired; §3.3 restated per-target under the canonical rule; paper quotes only replayable numbers. |
| **m10** Scope creep risk | Q4 held to one figure; YZ Cnc fallback held to one subsection; AN UMa capped at "apply pipeline, show fold, stop" (Q5). |
