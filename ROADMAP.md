# MACRO Consortium — Cross-Project Analysis Roadmap

**Committee Chair / Survey Scientist synthesis · 2026-08-16**
Inputs: the five finalized panel strategies (`TCrB_Monitoring`, `CV_TimeSeries`, `SN2023ixf_LightCurve`, `BeStar_Grism`, `DwarfGalaxy_AGN_Survey` — each `ANALYSIS_STRATEGY.md`, all rev. 2026-08-16, post-internal-referee). Cross-cutting claims re-verified against `"/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite"` on 2026-08-16 (T CrB last frame 2025-06-24 confirmed; filter-'6' 2023 usage re-queried; archive-wide unsolved light-frame backlog counted at ~59,100 rawimage rows).

This document governs shared infrastructure, sequencing, and resources. Per-project science decisions stay with the project strategies; where this roadmap and a project strategy conflict on shared machinery, this roadmap wins and the discrepancy is raised at the next consortium call.

---

## 0. Portfolio at a glance

| Project | Headline product | Venue posture | State |
|---|---|---|---|
| **TCrB_Monitoring** | Hα EW/profile series through the 2025 pre-eruption dip recovery (247 grism spectra / 60 nights) + photometric anchors | ApJ | Gated on P0 (restart obs, detector campaign, filter forensics) |
| **CV_TimeSeries** | State-tagged, cycle-resolved 3-color light curves of 3 polars + YZ Cnc superhumps (~7,400 + 1,920 frames) | ApJ | Gated on astrometry go/no-go + AAVSO cross-match |
| **SN2023ixf_LightCurve** | gri/narrowband campaign +5.4→+50 d; possible flash-phase grism series (1,052 unique frames / 35 nights) | AJ/PASP base, ApJ upside | Gated on Gate 0 (manifest, saturation census, grism triage) |
| **BeStar_Grism** | ~3-day-cadence Hα EW monitoring of bright Be stars, 2024-12→2026-07 (~3,400 core frames) | ApJ if ≥4 verified-active, else AJ/PASP | Gated on Step −1 (BeSS states, λ Eri injection, era-C FITS) |
| **DwarfGalaxy_AGN_Survey** | Hα fluxes/limits for 13 dwarf-candidate fields + NGC 5238 + NGC 5548 band-integrated LC | ApJ (NGC 5548 section conditional) | Gated on Phase 0 (dedup/pointing, filter dossier, Cannon) |

Every panel's internal referee round traced its worst errors to the same three root causes: **frame accounting** (duplicates, tree-doubling, alias fragmentation), **detector unknowns** (High Gain ceiling, StackPro noise, missing mode-matched calibrations), and **filter identity** (single-character MaxIm slot codes). These are facility problems, not project problems. They get built once, below.

---

## 1. The shared pipeline: `macro_pipeline`

One versioned monorepo (suggested location: `/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO/macro_pipeline/`), conda env `rlmt-checks` with pinned `environment.yml`. Subpackages: `macro_core` (catalog/manifest/astrometry/timing), `rlmt_diagnostics` (detector + calibration library and its documentation), `macro_phot` (ensemble photometry + validation), `macro_ts` (time-series toolkit), `macro_grism` (slitless spectroscopy). Project repos import these; **no project may fork shared code into its own tree** — divergence between papers is the referee trap the SN and T CrB panels both flagged (filter table, linearity curve, and flat solutions must be byte-identical across papers).

Binding conventions, adopted portfolio-wide (each is a lesson one panel paid for):

1. **Catalog is read-only ground truth for inventory; pixels come from `rawimage/`.** The `reduced/` tree is unaudited and demonstrably wrong in at least one case (BeStar Phecda pair).
2. **Every number in every paper is script-emitted** (BeStar's `make tables` policy, generalized). Published SQL appendix per paper.
3. **Dedup is global on (basename, jd) across and within trees** (SN lesson: wholesale duplicate nights live *inside* `rawimage/`).
4. **Era/camera assignment keys on (READOUTM, NAXIS geometry, EGAIN), never filter name or date** (BeStar lesson: mislabeled hrg/lrg frames, era-C repackaging).
5. **Night label = local-noon-to-noon, `date(JD − 0.7917)`** (Dwarf convention, adopted everywhere).
6. **Times: mid-exposure BJD_TDB computed from scratch; header JD (UTC exposure start) and JD-HELIO are never used.**
7. **Header airmass, camtemp, and era-C NAXIS are untrusted**; recompute airmass from coords+time, re-scrape temperatures from FITS headers, resolve HDU layout by inspection with hard-fail on unknown packaging.
8. **ZMAG is QC only, never calibration**; the primary cloud veto is the ensemble-flux-ratio test (zmag doesn't exist for most pyscope-era frames).

### 1.1 Stages

| Stage | Contents | Software | Built against (first consumer) |
|---|---|---|---|
| **S0 Manifest & curation** (`macro_core.manifest`) | Global dedup; alias/regex table + per-alias reconciliation report; tree policy (rawimage canonical, documented exceptions e.g. NGC 5548 macalester superset); coordinate-cone cross-checks; pointing-validation columns (header→target, then solved-WCS→target); era tagging; QC flags (exptime≤0, garbage airmass); night labels; header re-scrape hooks | sqlite3, astropy, pandas | All five simultaneously (week 1) |
| **S1 Astrometry** (`macro_core.astrom`) | Stratified go/no-go re-solve experiment (per era/binning/exptime, 200-frame samples) → batch `solve-field` (local 4100-series indices, SIP-2 with residual-curvature check) → escalation path: per-pointing WCS propagation by offset tracking | astrometry.net, astropy | CV (its Phase 0.5 is the experiment design; SN, Dwarf, NGC 5548 ride the same batch) |
| **S2 Detector & calibration library** (`rlmt_diagnostics`) | Ceiling/bit-depth reconciliation memo (see §4-R1); per-mode linearity from archival ladders (T CrB R/I ladders; SN 0.5 s/2 s pairs; BeStar half/full pairs + 2024-05-20 Vega ladder); **StackPro PTC** from the named 2023-06-07/08 repeated darks (N_sub, variance model, working ceiling); shutter-timing test at 0.085–0.1 s; master darks/biases per (mode, exptime); mode-matched flats where they exist + delta-sky-flat builder + measured cross-mode flat penalty terms; fringe maps; custody of the new calibration acquisitions | astropy, photutils | T CrB P0-3 + CV §5 rows 1–2 + Dwarf 0.4 are the same campaign — run once |
| **S3 Timing** (`macro_core.timing`) | Mid-exposure BJD_TDB (JPL ephemeris, Winer EarthLocation); DATE-OBS era audit (MaxIm vs pyscope; StackPro mid-time semantics); clock validation on an archived EB (AG LMi if ≥3 in-eclipse nights vs the corrected VSX period 1.3590176 d, else a TESS-era backup EB) | astropy.time | CV Phase 1; consumed by all |
| **S4 Ensemble photometry** (`macro_phot`) | Aperture photometry (apertures in arcsec per camera, per-night); ATLAS-REFCAT2 comp selection (mag/color/isolation/RUWE/variability cuts); Honeycutt-1992 inhomogeneous ensemble with robust ZP + linear color term + second-order k″·(color)·X extinction; per-mode saturation vetoes from S2; ensemble-flux-ratio cloud veto (ZMAG as secondary QC); forced photometry + upper limits; natural-system output + transformations-as-metadata (never transform the science target) | photutils, sep | SN Step 4–5 (cleanest testbed), then CV, T CrB Phase B, Dwarf Phase 2, BeStar Tier 3 |
| **S5 Error model & validation** (`macro_phot.validate`) | Empirical check-star RMS vs mag per (camera, mode, band); χ²-validated per-point errors with reported inflation factors; Allan deviation; the standard validation figure set (every paper's RMS-vs-mag and constant-star-χ² figures from one code path) | — | SN Sec. 5.1, then all |
| **S6 Time-series toolkit** (`macro_ts`) | LS + PDM + conditional entropy with published spectral windows; night-block bootstrap FAPs; alias tests (±1 c/d rule, Dawson–Fabrycky); injection–recovery engine (sinusoids/templates at real timestamps through the full pipeline, 90%-recovery contours that close); celerite2 joint GP+signal fits (never detrend-then-search); emcee/dynesty posteriors | astropy.timeseries, PyAstronomy, celerite2, emcee, dynesty | CV Phase 3 + BeStar Step 11 (nearly identical specs); T CrB C, SN Step 8, Dwarf 5.3–5.5 |
| **G Grism track** (`macro_grism`) | HDU-resolution FITS reader (hard-fail on unknown packaging); **per-frame identity gate** (Gaia zero-order field-pattern match, no WCS needed); trace + Horne optimal extraction with flanking-band background (doubles as Mode0 dark remover); wavelength (per-frame zero point self-anchored + telluric O₂ anchors; dispersion per filter/era); response vs CALSPEC with contaminated-region interpolation + airmass term; EW machinery (fixed windows, Vollmann–Eversberg + empirical floors); saturation triage | specreduce, specutils, synphot | T CrB Phase A builds it; BeStar extends (era chain, standards, tiers); SN Gate 0c consumes; the future **T CrB eruption letter reuses it verbatim** |

### 1.2 Where per-project forks happen

| Project | Consumes | Owns (never shared) |
|---|---|---|
| T CrB | S0–S6 + G | Archival-snippet flickering upper-limit table; Munari EW convention splice; θ CrB Be-characterization; AAVSO/ARAS co-analysis; eruption-response plan |
| CV | S0–S6 | Bright-phase sigmoid timing (per band); accretion-state mixture model; superhump Kato-style O–C; per-camera color discipline (VV Pup); cycle-count ambiguity analysis |
| SN 2023ixf | S0–S5 (+S6 for limits) + G (triage only) | Saturation matrix (filter × night); template subtraction (HOTPANTS vs sfft bake-off, crosswalk-gated templates); narrowband forward-modeling through measured bandpasses; S-corrections |
| BeStar | S0–S3 + G + S4 (Tier-3 imaging) + S6 | Era transfer chain (Vega/η Hya/θ Vir); two-floor error model from standards; BeSS validation; photospheric-corrected disk radii |
| Dwarf/AGN | S0–S5 (+S6 for census) | LSB stacking (SWarp/reproject, constant-sky policy); Román+2020 depth + synthetic-dwarf injection + detectability gate; imfit Sérsic; Hα continuum subtraction ([N II], R over-subtraction); host-aware AGN aperture + PyCALI intercalibration |

### 1.3 Cross-project synthesis catches (chair's findings, 2026-08-16)

These fell out of reading the five strategies side by side; each is assigned in §4/§5.

- **C1 — The High Gain ceiling is already effectively measured.** The T CrB panel left bit-depth as an open M1-vs-M2 dispute (P0-3(i) arbitration pending). But the SN panel *measured* a hard clip at ~3,530–3,550 ADU on the same GSENSE4040 High Gain channel, and the CV panel independently adopted ~3,500 ADU with a 3,200 ADU veto. Cross-project weight of evidence: **M2's low-ceiling reading is right; T CrB's 8 s R frames (peak 3,565 ADU) sit at the clip and the ruling-2 "no archival R survives" branch should be planned as the default** (anchors on B/I + P0-2 rescues). The one-afternoon hardware check remains scheduled — as confirmation, and to reconcile the three panels' numbers in one `rlmt_diagnostics` memo that all five papers cite.
- **C2 — Filter slot '6' is claimed as two different things in spring 2023.** SN panel: '6' = slitless grism (verified: 261 rawimage rows on 2023ixf, May 22–Jun 24, plain High Gain, 64–512 s, 0 solved, first frame labeled "M101 nova"). Dwarf panel: '6' = "likely luminance" for NGC 5548 (Mar 24–Apr 26, StackPro, 256 s, 0 solved) — and the chair's re-query shows the same code on NGC 4151, NGC 4725, NGC 5394 and standard stars through spring 2023, all StackPro. Either the wheel changed in early May 2023, or one panel is wrong — **and if NGC 5548's '6' is a grism, the AGN light curve (Dwarf Q4) is dead as photometry**. Resolution is one file-opening session (a dispersed trace is unmistakable): scheduled Week 1, feeds the filter dossier.
- **C3 — QQ Gem is orphaned.** CV cut it ("belongs in the grism program paper"); BeStar never adopted it (not in its target table; 122 hrg + 41 lrg frames idle). BeStar Step 0 must explicitly disposition QQ Gem (comparison pool / sample extension / drop with reason).
- **C4 — θ CrB's Be nature must travel with the shared response solutions.** T CrB's calibrator characterization (Hα-region interpolation flag) is metadata on `macro_grism` response products, so BeStar and the eruption letter can never accidentally use θ CrB Hα as an anchor.
- **C5 — AG LMi serves three roles**: cut from CV science, promoted to shared clock validator (S3), and queued as an undergraduate EB spin-off note. One dataset, one reduction, three consumers — reduce it once under `macro_phot`.
- **C6 — The astrometry batch should be scoped deliberately.** Projects need ~5,100–6,500 solves (CV polars ~4,600; SN ~177 broadband; Dwarf ~250 + NGC 5548 133 + NGC 5238 remainder). The archive-wide unsolved backlog is ~59,100 light frames — a facility-wide solve is ~2–3 days of Mac Studio time and worth doing *once the per-era success rates are known*, but it is not on any paper's critical path. Decision at the Week-2 review.

---

## 2. Sequencing

**One line: Facility gates first (Weeks 1–3, all projects' Phase-0s merged), then T CrB (grism spine, ApJ) in parallel with SN 2023ixf (photometric-core testbed), then CV, then BeStar, then Dwarf/AGN — with the T CrB eruption letter preempting everything if the nova goes.**

### Wave 0 — Facility gates + decision gates (now → ~2026-09-05)
All five strategies' Phase-0/Gate-0/Step-−1 tasks, merged and deduplicated (§5 sprint plan). Cheap (days each), and every science claim in the portfolio hangs on them. Includes the time-critical operational items: **T CrB observation restart** (dark since 2025-06-24, 14 months into Schaefer's window; the 2026 flickering science exists only if frames start flowing now) and the **consolidated calibration-frame acquisition request** (Mode0 240 s darks/biases; current-mode flats; era-C grism bias/dark/flat; the era-B-camera-survival question; ST LMi g/r/i season for the queue).

### Wave 1 — T CrB ∥ SN 2023ixf (Sep → Dec 2026)
- **T CrB leads on science urgency**: the eruption window is nearly exhausted, Munari-adjacent groups publish continuously, and the 2025 Hα recovery series is the portfolio's strongest unique asset. Building `macro_grism` here (identity gate, extraction, wavelength, response, EW) also creates the consortium's most reusable novel infrastructure — BeStar and the eruption letter inherit it whole. Target: submission Dec 2026–Jan 2027.
- **SN 2023ixf runs in parallel as the photometric core's validation vehicle**: cleanest detector story in the archive (plain High Gain, no StackPro), one 35-night season, healthiest plate-solve state, and a base-case AJ/PASP venue that keeps stakes low while S0/S1/S2/S4/S5 harden. It is also the most student-ready project (Gate 0 census and filter forensics are well-specified, bounded tasks). Venue decision at its week 3, per its strategy. Target: submission ~Jan 2027.
- Rationale for *not* starting CV first despite its size: its pipeline is blocked by the astrometry go/no-go and the StackPro PTC — both Wave-0 outputs — and its analysis layers (timing, states) are the deepest, benefiting most from a hardened core.

### Wave 2 — CV ∥ BeStar (Nov 2026 → Mar 2027)
- **CV** starts its photometry as soon as the Wave-0 batch re-solve lands (ST LMi first, per its own execution order), overlapping T CrB's writing phase. Strongest remaining ApJ unit; a grad-student-ready execution order already exists in its strategy. Target: submission Mar–Apr 2027 (later if the requested extra ST LMi g/r/i season is worth waiting for — decide at the February call).
- **BeStar** starts when `macro_grism` is mature and the era-C calibration frames + Step −1 verdicts are in. Its re-baselined Steps 0–4 (5–6 weeks) ride on T CrB-built extraction code. Venue set by the Step −1 verified-active count, per its strategy. Target: submission Feb–Mar 2027.

### Wave 3 — Dwarf/AGN (Feb → May 2027)
Gated on Cannon (configs, Hα transmission curve, W flats, Dw1643+07 provenance, co-authorship) — contact happens in Week 1, but the science schedule tolerates his latency. Its LSB fork shares the least pipeline surface, so it gains the most from a stable core. Hα-led framing insulates it from ELVES/HSC classification scooping. Target: submission Apr–May 2027.

### Standing preemption
If T CrB erupts: the pre-loaded response plan triggers, and the **eruption letter becomes the consortium's top priority** — its entire feasibility rests on `macro_grism` homogeneity with the 2025 series, which is another reason Track B leads. The pyscope eruption-response block must be verified loaded by **2026-09-01** (T CrB §9).

---

## 3. Resource estimates

### 3.1 Compute (Mac Studio, local)
Nothing in the portfolio needs external compute. The binding constraint is HDD I/O on the 3.34 TiB archive: stage per-project frame subsets to internal SSD scratch, read raw pixels once per stage, write intermediate products (extracted spectra, photometry tables) under each project's `products/`.

| Task | Size | Wall-clock estimate |
|---|---|---|
| S0 manifest (archive-wide, headers only) + re-scrape of ~25k grism headers | 329k rows / ~25k files | Hours; one overnight including re-scrape |
| S1 go/no-go samples (4 × 200 frames) | 800 solves | ~1–3 h (10–12 parallel workers) |
| S1 project-critical batch solve | ~5,100–6,500 frames | Overnight run |
| (Optional) facility-wide solve | ~59k frames | ~2–3 days; not on critical path |
| S2 linearity/PTC fits | Hundreds of frames | Hours |
| Aperture photometry per project | 1k–9k frames | Minutes–hours each |
| Grism extraction (T CrB 247 + θ CrB 403; BeStar core ~3.4k + standards) | ~4.5k frames | Hours per full pass |
| Injection–recovery (CV, BeStar, Dwarf census; full-pipeline trials) | ~10⁴–10⁵ trials each | The heavy item: 1–3 days multicore per project; design for resumable batches |
| MCMC/GP fits (emcee, celerite2) | Per target | Hours |
| Dwarf coadds + synthetic-dwarf injection | ~500 frames / 19 fields | Hours–1 day |

### 3.2 Human effort (FTE-weeks; grad-student grade with chair oversight)

| Project | Gates | Reduction/analysis | Writing | Total | Lead |
|---|---|---|---|---|---|
| Shared pipeline S0–S3 + S2 campaign | 2 | 3 (amortized) | — | ~5 | Detector-science seat + chair |
| T CrB | 2 | 8 (Phase A 5–6, B 2, C 1–2) | 4 | **14–18** | Chair + grad student |
| SN 2023ixf | 2 (incl. grism-triage timebox) | 6 | 3 | **10–13** | Grad student + undergrads |
| CV | 2 | 12 | 4 | **18–22** | Grad student |
| BeStar | 1 | 11 (Steps 0–4 ≈ 5–6 per its re-baseline) | 4 | **15–18** | Grad student |
| Dwarf/AGN | 2 | 9 | 3–4 | **13–16** | Grad student + Cannon liaison |

Portfolio ≈ 70–90 FTE-weeks. With two grad students, undergraduate task packets (AAVSO/BeSS cross-matches, saturation-census QC, filter-forensics record digging, EB spin-off note), and the chair on T CrB, the Wave schedule above is feasible inside ~9 months.

### 3.3 External data (fetch once, cache under `MACRO/external_data/<source>/` with pull dates)
- **Shared:** ATLAS-REFCAT2 local tiles for all fields (fetch first — every project's S4 needs them); Gaia DR3 (astrometry, RUWE, variability flags, XP synthetic photometry); PS1; CALSPEC (`alpha_lyr_stis`, B6 template); astrometry.net 4100-series indices (already local — verify coverage).
- **T CrB:** AAVSO AID (B/V/Vis 2023–2026), ASAS-SN Sky Patrol v2, ARAS T CrB spectra, Swift-XRT/UVOT, TESScut; ATel/AAVSO alert feed.
- **CV:** ZTF + ATLAS forced photometry, ASAS-SN, AAVSO (YZ Cnc first — decision gate), TESScut, eROSITA-DE DR1.
- **SN:** ZTF/ATLAS/ASAS-SN/AAVSO, WISeREP flux-calibrated spectra (feeds the narrowband forward model), MRT tables (Hosseinzadeh+23, Li+24, Chen+24), dustmaps.
- **BeStar:** BeSS (consortium `BeSS_Vis` repo — Step −1), TESScut sectors 2024–2026, AAVSO for the brightest.
- **Dwarf/AGN:** Karachentsev+2022/2025 via VizieR, Kaisin & Karachentsev Hα target lists, DESI target/redshift catalogs, Legacy Survey thumbnails + (μ₀, r_e), IRAS/Herschel/WISE cirrus cutouts, AGN STORM NGC 5548 host decomposition, ZTF/ASAS-SN/ATLAS NGC 5548 spring-2023 forced photometry (feeds the Phase 4.5 gate).

---

## 4. Cross-cutting risks and their resolvers

| # | Risk | Hits | Resolver / owner | When |
|---|---|---|---|---|
| R1 | **High Gain ceiling/bit depth** (~3.5 kADU clip; silent saturation at 5% of 16-bit range) | T CrB R-band anchors; CV YZ Cnc outbursts; SN early epochs | `rlmt_diagnostics` reconciliation memo (adopt C1 default: low ceiling; T CrB plans the no-archival-R branch now) + confirmatory afternoon hardware test + per-mode ladders. Owner: detector-science seat | Week 1–2 |
| R2 | **StackPro noise model unknown** (N_sub unverified; non-Poisson correlated noise; per-sub-read clipping) | T CrB archival snippets; CV MaxIm era; Dwarf NGC 5548/NGC 5238; every χ²/FAP on StackPro data | Single PTC campaign from the named 2023-06-07/08 dark repeats + star-based linearity; empirical errors + bootstrap FAPs everywhere until it lands. Owner: detector-science seat | Week 2 |
| R3 | **Missing mode-matched calibrations** (zero Mode0 darks/biases — under all 247 T CrB spectra; no StackPro/LowGain flats ever; no slot-'6' flat; no grism-era calib at all) | T CrB, BeStar, Dwarf | One consolidated acquisition request to Winer ops (chair signs, Week 1); archival frames get measured cross-mode penalty terms + the flanking-band adequacy test; delta-sky flats where dithered | Request Week 1; products this month |
| R4 | **Filter identity** (single-char MaxIm codes; slot maps epoch-dependent; **C2 slot-'6' conflict**) | T CrB P0-2; SN Step 1; Dwarf 0.2; NGC 5238 | One **filter-forensics dossier**, per-epoch, dual-track (config logs + Cannon + purchase records ∥ empirical color-term regressions on solved frames), covering science *and* calibration frames; the slot-'6' file-opening session settles C2. Owner: chair + archival seats | Weeks 1–3 |
| R5 | **Frame accounting** (every panel's first draft had wrong counts) | All | S0 manifest as sole source of counts; per-paper SQL appendices; READMEs regenerated from manifests (three are currently stale: T CrB, SN, Dwarf) | Week 1 |
| R6 | **Astrometry at scale** (~4,600 unsolved polar frames; failure cause on 2×2-binned ASI unknown) | CV primarily; SN, Dwarf | S1 stratified go/no-go experiment before commitment; escalation = WCS propagation by offset tracking; batch scoped per C6 | Weeks 1–2 |
| R7 | **Timing conventions** (header JD = UTC start; JD-HELIO unverified; StackPro mid-time semantics) | CV (fatal if wrong), all others | S3 shared BJD_TDB module + era audit appendix + EB clock validation (AG LMi count first) | Weeks 2–3 |
| R8 | **Scooping / time-criticality** (T CrB eruption; Munari-adjacent groups; ELVES/HSC classifications; Labadie-Bartz Be program) | T CrB, Dwarf, BeStar | Consolidated **weekly sweep** (one person, one hour: arXiv + ATels + AAVSO alerts for all five target lists); T CrB restart + response block by 2026-09-01; verify-at-submission flags on all target-status claims | Standing |
| R9 | **People/priority** (Cannon's program and data; single-machine, single-chair dependence) | Dwarf; all | Cannon contact Week 1 with co-authorship framing; pipeline + products in versioned repos (the archive is already a mirror; the catalog and derived products must be too); per-project leads named at the Week-2 review | Week 1–2 |
| R10 | **Calibrator astrophysics** (θ CrB is a Be/shell star; Spica is variable; standards coverage starts 2025-12) | T CrB, BeStar | C4 metadata flag on shared response products; BeStar's epoch-restricted detection rule; characterization-before-use ordering in both projects | With `macro_grism` build |

---

## 5. First two-week sprint (2026-08-17 → 2026-08-28)

Principle: every item below is a Phase-0/Gate-0 task some panel already declared blocking, merged here so no work is done twice. Owners in parentheses.

### Week 1 (Aug 17–21)
- **D1 — Ship the consolidated observatory request** (chair; half a day, highest leverage in the portfolio): (a) restart T CrB — nightly 3×1 s r + short B, one lrg + one hrg 240 s, weekly ≥2 hr B run at ≤60 s cadence; load the eruption-response pyscope block; (b) calibration acquisitions — Mode0 240 s darks + biases, flats in every active mode, era-C grism bias/dark/flat set, the era-B-camera-survival question; (c) queue request: one more g/r/i ST LMi season (up through spring).
- **D1 — Contact John Cannon** (chair): 2023 filter-wheel configs + Hα transmission curve (center/width) + W flats + Dw1643+07 provenance + co-authorship conversation.
- **D1–D3 — S0 manifest, archive-wide** (pipeline lead): global (basename, jd) dedup, alias table, era tagging, pointing-validation columns, night labels; re-emit each project's Table-1 counts and **diff against the five strategy documents** — discrepancies to the Week-2 review.
- **D2–D3 — Decision-gate lookups** (undergrad packets, parallel): AAVSO cross-match of YZ Cnc 2024-02-21→05-03 (decides CV Q3); BeSS emission-state check for the ten Be targets (fixes BeStar venue); start λ Eri injection–recovery setup.
- **D2 — File-opening session** (detector seat): era-C repackaged FITS sample (HDU layout) + hrg O₂ B-band coverage check + **the C2 slot-'6' session** — open NGC 5548 and 2023ixf filter-'6' frames side by side; a dispersed trace settles it in minutes. Outcome routes either to the Dwarf strategy (Q4 dies/survives as photometry) or closes the conflict.
- **D3–D4 — R1 ceiling work** (detector seat): the one-afternoon bit-depth hardware check; start per-mode exposure-ladder fits from archival ladders; draft the reconciliation memo adopting the C1 default.
- **D4–D5 — S1 go/no-go experiment** (grad student): 200-frame stratified samples for ST LMi-Sloan / VV Pup / EU UMa / AN UMa; report per-stratum success rates.
- **D4–D5 — SN Gate 0a+0b** (SN student): manifest freeze from S0 output; saturation-census script (peak-ADU at SN position, all 1,052 frames + the saturated early epochs); saturation matrix draft.

### Week 2 (Aug 24–28)
- **Batch re-solve** the project-critical ~5–6.5k frames if the experiment says go (overnight runs); on 70–95% success, quantify per-filter night losses; on <70%, invoke the WCS-propagation redesign before any CV Phase-2 work.
- **StackPro PTC** from the 2023-06-07/08 dark repeats, both modes; amp-glow check on 256 s darks; publishable gain/RN/ceiling numbers into `rlmt_diagnostics` (detector seat).
- **S3 timing module** + DATE-OBS era audit; count AG LMi in-eclipse nights against the corrected VSX period — pick the clock validator (time-series seat).
- **SN Gate 0c grism triage begins** (2-week timebox, hardened criteria: paired direct images verified, named wavelength source, offset-trace contamination test).
- **`macro_grism` identity-gate prototype** (chair + grad student): Gaia zero-order pattern match, run on the 21 known-bad T CrB pointings as the test set; design shared with BeStar Step-1 QC.
- **Filter dossier assembly**: color-term regressions on newly solved frames (broadband IDs); collate Cannon/config-log responses; slot-'6' verdict written up.
- **λ Eri injection–recovery run** completes (expected: confirms slow-tier demotion; its completeness map is a paper figure).
- **Fri Aug 28 — Week-2 review** (whole committee). Agenda: manifest-vs-strategy count diffs; astrometry verdict + C6 batch-scope decision; ceiling memo adopted (T CrB R-branch confirmed or reprieved); slot-'6' verdict; YZ Cnc outburst-state verdict; BeSS verified-active count (BeStar venue set); per-project leads named; Wave-1 kickoff authorized.

**Sprint exit criteria:** manifest is the single source of counts; astrometry path decided; detector memo v1 adopted; filter dossier v1 (broadband IDs + slot-'6') issued; all five projects' venue/scope branch points resolved or explicitly scheduled; T CrB frames flowing (or escalated to the consortium board).

---

## 6. Standing cadence

- **Weekly:** consolidated literature/alert sweep (R8); T CrB restart status until frames flow; sprint-board triage.
- **Biweekly:** pipeline-change review — any modification to `macro_pipeline` after a paper freezes its numbers triggers regeneration of that paper's tables (the script-emission policy makes this cheap).
- **Monthly consortium call:** wave-gate reviews (Wave-1 kickoff ~Sep; CV/BeStar starts ~Nov; Dwarf start ~Feb); eruption-response drill status; student assignments.
- **At every submission:** verify-at-submission flags (target statuses, competitive landscape tables) re-checked; README regenerated from the manifest; Zenodo/data-release checklist.
