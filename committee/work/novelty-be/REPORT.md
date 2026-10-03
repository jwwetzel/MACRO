# REPORT — package `novelty-be` (BE-S-1a-bess, BE-N1-gate, QQ Gem, θ CrB)

Date: 2026-10-03 · BeSS index pulled 2026-10-03 06:41–06:54 UTC · manifest S0 v1.0 (built 2026-08-19, commit d97c6d8), opened read-only.

Every number and table below is emitted by `BeStar_Grism/notes/novelty/bess_novelty_check.py` from
`BeStar_Grism/notes/novelty/novelty.sqlite` (tables `nv_*`); the Markdown tables are pasted verbatim from
`BeStar_Grism/notes/novelty/tables/*.md`. Literature statements carry their source and how far I verified it (§7).

---

## 1. Verdict in five lines

1. **Core ten: zero BeSS-verified active emitters.** Six of the ten are in BeSS; in 58 independent BeSS Hα spectra
   bracketing our seasons, none of 69 Ori, Phecda, 5 Cnc, φ Leo, HD 70340, 53 Boo shows emission above the continuum.
   λ Eri is the one partial exception: 1 of 7 spectra (2025-10-28, 15 d before our first night) shows a weak peak
   (F/Fc = 1.07), after which the line filled back toward pure absorption through our season.
2. **BE-N1 under the editor's rule: < 2 ⇒ the Be-star paper *as scoped in ANALYSIS_STRATEGY.md* stops as a science paper.**
3. **But the strategy drew the wrong sample.** The same grism archive holds **19 BeSS-verified active Be stars with
   ≥ 10 of our nights in a season** (plus 9 with one witness or mixed evidence) that the strategy left in an unadopted
   "comparison pool". The gate fails for the paper as planned and passes for the campaign as observed. Whether to
   re-draw the sample is the chair's / James's decision (§6).
4. **QQ Gem = HD 46264 is a verified-active Be star in both seasons** (BeSS EW −4.2 to −6.0 Å, peak F/Fc 1.6–1.7;
   37 of our nights vs 4 bracketing BeSS spectra). Disposition: **sample extension**, not drop.
5. **θ CrB is a Be star that was disk-less in Hα throughout our campaign**: 17 BeSS spectra by 12 observers,
   2025-01-17 → 2026-07-20, all pure absorption, EW 5.37 Å mean, sd 0.34 Å, no trend (+0.05 ± 0.04 Å/100 d). It is a
   usable *line* reference; its *continuum* constancy is not verified here.

---

## 2. What was done

| Step | Subcommand | Result |
|---|---|---|
| Our coverage | `targets` | 16,449 canonical raw-tree grism Light frames (matches ANALYSIS_STRATEGY §3.1 exactly), 294 nights, 422 distinct target strings after stripping the grism/exposure suffix |
| Identity | `resolve` | 260 names with ≥ 3 nights sent to CDS Sesame/SIMBAD; 62 SIMBAD objects surveyed = core ten + θ CrB + QQ Gem + every other grism target with ≥ 10 nights; 4 names unresolved (none Be-related) |
| BeSS index | `index` | BeSS SSA cone query (72″) per object, 1990 → 2026-10-03, bisecting on the service's 1000-row overflow: **92,635 spectrum records, 9,101 covering Hα, for 46 of 62 objects** |
| Evidence | `fetch`, `measure` | all 581 Hα spectra within ±60 d of any of our seasons downloaded; **578 measured**, 3 unusable (Pleione low-res files, kept with status) |
| Deliverables | `tables`, `figures` | 7 tables, 3 coverage figures, 45 per-star profile sheets |

**Rule fixed in the script header before any spectrum was measured.** *Verified active* = ≥ 10 of our grism nights in
a season AND ≥ 2 measured BeSS Hα spectra within ±60 d of that season, ALL with a 1 Å-smoothed peak F/Fc ≥ 1.05 at
≥ 5σ of the equally smoothed continuum. EW is ∫(1 − F/Fc) over 6543–6583 Å against a linear continuum from
6520–6540 and 6586–6606 Å; **positive = absorption, negative = emission**.

**Bias tests (`test_bess_novelty_check.py`, 10 tests, all pass; synthetic spectra, truth known).** Signed EW bias
< 0.05 Å for absorption and emission profiles at 0.05, 0.2 and 2 Å sampling; pull width within 0.7–1.4; emission
flag false-positive rate < 1% on pure absorption at S/N 30 (the peak is a max-statistic, hence the smoothing and the
5σ condition); a 15% peak is always detected at S/N 50, a 2% peak never. **Documented blind spot:** emission that only
fills the absorption core (peak < continuum) is not flagged; it shows in the EW, which is why EW and its scatter are
tabulated beside every verdict and the profiles are plotted.

**The honest error on a BeSS EW is the observer-to-observer scatter, not the formal error.** Formal errors are
0.01–0.05 Å; the scatter between observers on stars that did not change is 0.25–0.45 Å (θ CrB 0.34, 69 Ori 0.40,
Phecda 0.37). Every statement of change below is made against the scatter. This number also matters downstream: it is
the precision our own series must beat to add anything to BeSS.

---

## 3. Core ten — per-target active/inactive table (BE-S-1a deliverable)

### Core ten (staged science targets) — our seasons vs BeSS

| Star (SIMBAD) [manifest name] | SpT (SIMBAD) | SIMBAD type | Season | Nights (hrg/lrg/Ha/OG) | Our span | BeSS Hα all-time | in season | within ±60 d | gap before / after (d) | measured | in emission | EW median [min, max] (Å) | EW scatter (Å) | formal err (Å) | peak F/Fc median (max) | BeSS verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| f01 Ori [69 Ori] | B5Vn | Be* | 2024-25 | 2 (2/2/0/0) | 2025-03-24 → 2025-04-01 | 156 | 0 | 7 | 6 / 174 | 7 | 0 | 4.70 [4.44, 4.86] | 0.21 | 0.02 | 0.98 (0.99) | NO-EMISSION |
| f01 Ori [69 Ori] | B5Vn | Be* | 2025-26 | 49 (49/0/0/0) | 2025-10-29 → 2026-04-11 | 156 | 15 | 18 | 15 / 9 | 18 | 0 | 5.07 [4.35, 6.13] | 0.24 | 0.02 | 0.98 (0.99) | NO-EMISSION |
| gam UMa [PHECDA; gam UMa] | A0V | Em* | 2024-25 | 40 (30/32/4/4) | 2024-12-19 → 2025-04-17 | 123 | 2 | 3 | 242 / 36 | 3 | 0 | 8.55 [8.43, 8.70] | 0.18 | 0.03 | 0.95 (0.95) | NO-EMISSION |
| gam UMa [PHECDA; gam UMa] | A0V | Em* | 2025-26 | 7 (0/7/0/0) | 2026-02-20 → 2026-04-23 | 123 | 2 | 5 | 272 / 9 | 5 | 0 | 8.76 [8.29, 9.48] | 0.27 | 0.05 | 0.95 (0.95) | NO-EMISSION |
| lam Eri [Lam Eri] | B2III(e)p | Be* | 2025-26 | 47 (46/1/0/0) | 2025-11-12 → 2026-03-26 | 170 | 5 | 7 | 14 / — | 7 | 1 | 2.86 [1.28, 3.98] | 0.96 | 0.02 | 1.00 (1.07) | MIXED |
| 5 Cnc [5 Cnc] | B9.5Vn | * | 2024-25 | 21 (12/12/9/5) | 2024-12-19 → 2025-03-26 | 103 | 4 | 5 | 51 / 206 | 5 | 0 | 6.31 [4.76, 6.77] | 0.22 | 0.03 | 0.97 (1.01) | NO-EMISSION |
| 5 Cnc [5 Cnc] | B9.5Vn | * | 2025-26 | 24 (24/0/0/0) | 2025-11-06 → 2026-05-18 | 103 | 6 | 7 | 18 / — | 7 | 0 | 6.82 [6.52, 7.91] | 0.29 | 0.05 | 0.95 (0.96) | NO-EMISSION |
| phi Leo [phi Leo] | A5V | PM* | 2024-25 | 39 (30/26/7/5) | 2024-12-19 → 2025-04-17 | 146 | 6 | 7 | 36 / 247 | 7 | 0 | 6.81 [6.34, 7.06] | 0.26 | 0.02 | 0.98 (0.98) | NO-EMISSION |
| eta Hya [HR 3454] | B3V | bC* | 2025-26 | 34 (34/25/0/0) | 2025-12-04 → 2026-05-11 | 0 | 0 | 0 | — / — | 0 | 0 | — | — | — | — | UNWITNESSED |
| tet Vir [HR 4963] | A1IVs | ** | 2025-26 | 33 (33/32/0/0) | 2025-12-13 → 2026-06-12 | 0 | 0 | 0 | — / — | 0 | 0 | — | — | — | — | UNWITNESSED |
| HD 70340 [HD 70340] | A0V | ** | 2024-25 | 30 (17/17/13/11) | 2024-12-19 → 2025-04-09 | 36 | 2 | 2 | 239 / 358 | 2 | 0 | 8.61 [8.40, 8.81] | — | 0.03 | 0.95 (0.96) | NO-EMISSION |
| alf Vir [Spica] | B1V | bC* | 2025-26 | 29 (29/27/0/0) | 2026-01-12 → 2026-04-07 | 0 | 0 | 0 | — / — | 0 | 0 | — | — | — | — | UNWITNESSED |
| nu.02 Boo [53 Boo] | A2V | ** | 2024-25 | 23 (23/23/0/0) | 2025-02-26 → 2025-04-17 | 126 | 1 | 4 | 252 / 10 | 4 | 0 | 8.02 [7.74, 8.23] | 0.24 | 0.04 | 0.96 (0.96) | NO-EMISSION |

Columns "gap before / after" are days from our first / after our last night to the nearest BeSS Hα spectrum outside
the season. Figure: `BeStar_Grism/notes/novelty/figures/coverage_core_ten.png`; profiles:
`figures/profiles_<star>.png`.

**Reading, star by star** (literature sources in §7):

| Star | Classical Be? | In emission during our seasons? | Verdict |
|---|---|---|---|
| λ Eri | Yes — B2 Be, prototype of the class [L1, L2] | Weak and fading. Peak 1.07 on 2025-10-28; EW 1.28 → 3.98 Å by 2026-02-25, slope +1.75 ± 0.37 Å/100 d (N = 7, 4 observers; 4.8σ for one of six stars tested — treat as ≈ 3σ after trials). A small outburst decaying through our 47 nights. | **Weakly active; not verified by the pre-set rule** (1 of 7) |
| 69 Ori (f¹ Ori) | Yes historically — Be star that "began behaving as a normal star in November 1982" [L3, L4] | No. 25 spectra, 14 observers, peak ≤ 0.99, EW 4.99 ± 0.44 Å | **Inactive (disk-less)** |
| 5 Cnc | Reported Be with a weak disk in 1999 [L5]; SIMBAD type is plain star | No. 12 spectra, peak ≤ 1.01, EW 6.60 Å. One high-resolution spectrum (2024-12-23, R ≈ 40,000) has EW 4.76 Å, 1.8 Å below the rest — a single witness of possible core infill, 4 d after our first night; not a detection | **Inactive** (one unconfirmed infill) |
| Phecda (γ UMa) | Ae star in the Jaschek & Andrillat sense [L6]; SIMBAD "Em*" | No. 8 spectra, peak 0.95, EW 8.70 ± 0.36 Å, flat | **Inactive** — as the strategy assumed (null-test star) |
| φ Leo | A-type shell star, not a classical Be; variable circumstellar absorption, δ Sct [L7, L8] | No Hα emission. 7 spectra, EW 6.75 ± 0.26 Å | **Inactive in Hα** |
| 53 Boo (ν² Boo) | A-type shell star [L9], visual binary; not a classical Be | No. 4 spectra, EW 8.0 Å, peak 0.96 | **Inactive** |
| HD 70340 | In the BeSS catalogue; SIMBAD A0V double star; I found no literature on emission | No. 2 spectra, EW 8.6 Å, peak 0.95; nearest other BeSS spectra 239 d before / 358 d after | **Inactive** (thinly witnessed) |
| Spica | No — β Cep / ellipsoidal SB2 [L10] | not in BeSS | reference star |
| η Hya (HR 3454) | No | not in BeSS | standard |
| θ Vir (HR 4963) | No | not in BeSS | standard |

**Count of verified-active emitters in the core ten: 0** (strict rule); **1** if λ Eri's weak, fading episode is
admitted. Either way < 2.

Trend table (script-emitted; slope error from residual scatter, no χ² rescaling):

### BeSS Hα EW trend across the campaign (core ten + brief-named; ≥ 5 spectra)

| Star [manifest name] | N spectra | Observers | Span | EW mean (Å) | EW sd (Å) | slope (Å / 100 d) | ± (from residual scatter) | slope / σ | residual sd (Å) |
|---|---|---|---|---|---|---|---|---|---|
| f01 Ori [69 Ori] | 25 | 14 | 2025-02-01 → 2026-04-20 | 4.99 | 0.44 | +0.12 | 0.05 | +2.4 | 0.40 |
| gam UMa [PHECDA; gam UMa] | 8 | 6 | 2024-12-23 → 2026-05-29 | 8.70 | 0.36 | +0.05 | 0.06 | +0.8 | 0.37 |
| lam Eri [Lam Eri] | 7 | 4 | 2025-10-12 → 2026-02-25 | 2.83 | 0.95 | +1.75 | 0.37 | +4.8 | 0.44 |
| 5 Cnc [5 Cnc] | 12 | 8 | 2024-10-29 → 2026-05-03 | 6.60 | 0.73 | +0.24 | 0.08 | +2.9 | 0.56 |
| phi Leo [phi Leo] | 7 | 4 | 2024-11-13 → 2025-04-01 | 6.75 | 0.26 | +0.28 | 0.21 | +1.3 | 0.25 |
| tet CrB [tet CrB] | 17 | 12 | 2025-01-17 → 2026-07-20 | 5.37 | 0.34 | +0.05 | 0.04 | +1.1 | 0.34 |

---

## 4. θ CrB and QQ Gem (named in the brief)

### Targets named in the brief (θ CrB, QQ Gem)

| Star (SIMBAD) [manifest name] | SpT (SIMBAD) | SIMBAD type | Season | Nights (hrg/lrg/Ha/OG) | Our span | BeSS Hα all-time | in season | within ±60 d | gap before / after (d) | measured | in emission | EW median [min, max] (Å) | EW scatter (Å) | formal err (Å) | peak F/Fc median (max) | BeSS verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tet CrB [tet CrB] | B6Vnne | Be* | 2024-25 | 21 (21/21/0/0) | 2025-03-09 → 2025-04-16 | 317 | 0 | 6 | 51 / 24 | 6 | 0 | 5.24 [5.01, 5.53] | 0.25 | 0.02 | 0.98 (0.99) | NO-EMISSION |
| tet CrB [tet CrB] | B6Vnne | Be* | 2025-26 | 19 (19/0/0/0) | 2026-01-09 → 2026-06-01 | 317 | 6 | 11 | 146 / 19 | 11 | 0 | 5.41 [4.87, 6.29] | 0.35 | 0.03 | 0.98 (0.99) | NO-EMISSION |
| HD 46264 [QQ Gem] | Be | Be* | 2024-25 | 11 (8/7/3/2) | 2024-12-27 → 2025-03-31 | 28 | 0 | 2 | 0 / 256 | 2 | 2 | -5.27 [-6.00, -4.53] | — | 0.06 | 1.69 (1.73) | EMISSION |
| HD 46264 [QQ Gem] | Be | Be* | 2025-26 | 26 (26/0/0/0) | 2025-10-28 → 2026-04-24 | 28 | 2 | 2 | 305 / — | 2 | 2 | -4.44 [-4.64, -4.24] | — | 0.03 | 1.60 (1.61) | EMISSION |

**θ CrB (the T CrB calibrator).**
- *Be status:* confirmed — SIMBAD Be*, B6Vnne; in the BeSS catalogue (317 Hα spectra since 1991); a "bright Be-shell
  star" [L11]; a 2024 interferometric study describes Hα as "dominated by photospheric absorption" with only weak
  emission in the wings [L12]; an L-band study classes it among Be stars that have lost all or most of their disk [L13].
- *During our campaign:* 17 BeSS spectra, 12 observers, 2025-01-17 → 2026-07-20: **no emission above the continuum in
  any** (peak 0.965–0.986), EW mean 5.37 Å, sd 0.34 Å (= the inter-observer floor), slope +0.05 ± 0.04 Å/100 d.
  Figure `figures/profiles_tet_CrB.png`.
- *Consequences.* (i) The physicist's note that `g_extractions` finds an Hα **emission** peak at SNR 220–420 in θ CrB
  lrg frames (PH memo, TCRB-A4) is contradicted by 17 independent spectra: that "emission" is an artefact of the
  extraction/anchor code, as the observational astronomer argued (OA memo E1). This is evidence for the grism package,
  not a change I made. (ii) θ CrB is a sound Hα-*absorption* reference for 2025–26 (EW 5.4 Å, stable to the 0.3 Å floor).
  (iii) Its *continuum* is not certified: it is a historically variable Be star (0.7 mag fade in 1970 [L14]); I
  found no 2025–26 photometry in this package. The synthesis decision to drop A4 (θ CrB response → absolute flux)
  stands; nothing here rescues it.
- BeSS coverage gap: no BeSS Hα spectrum inside our 2024-25 θ CrB season (nearest 51 d before, 24 d after).

**QQ Gem (ROADMAP C3).**
- Identity: SIMBAD resolves QQ Gem to **HD 46264**, type Be*, spectral type "Be"; in the BeSS catalogue under the name
  "QQ Gem" (28 Hα spectra since 2009). Header pointings of our frames agree (98.3°, +17.0°). I found no dedicated
  paper; it appears only in survey lists (Hipparcos unsolved variable). No published period or EW series found.
- Our data: 141 canonical raw grism Light frames on 37 nights (hrg 109, lrg 22, HaGrism 4, OGGrism 6);
  2024-25: 11 nights; 2025-26: 26 nights, hrg only. (ROADMAP C3 quotes "122 hrg + 41 lrg"; that count is not
  the Step-0 selection — ledger note in §8.)
- BeSS during our seasons: 4 spectra, all in emission, EW −6.00, −4.53 (Nov–Dec 2024), −4.24, −4.64 Å (Dec 2025–Jan 2026),
  peak F/Fc 1.59–1.73. None inside our 2024-25 season; 2 inside 2025-26.
- **Disposition: SAMPLE EXTENSION.** It is verified active in both seasons, we have ≈ 9× BeSS's sampling (37 nights vs 4 bracketing spectra), and at
  V ≈ 7.65 it is the faintest and least-studied star of the set. Caveat for the pipeline: exposures are 15–300 s, so it
  sits in a different exposure regime from the bright core; and the two filename variants (`QQ Gem`,
  `QQ Gem hrg 1-2e+02s`) must be merged in the alias table.

---

## 5. The unadopted pool — where the active Be stars actually are

### BE-N1 count (rule fixed in the script header before any spectrum was measured)

Verified active = SIMBAD-resolved star with ≥ 10 of our grism nights in a season AND ≥ 2 measured BeSS Hα spectra within ±60 d of that season, ALL with smoothed peak F/Fc ≥ 1.05 at ≥ 5σ.

| Set | Objects | Verified active (BeSS) | n | One-witness or mixed | n |
|---|---|---|---|---|---|
| Core ten | 10 | — | 0 | lam Eri | 1 |
| Brief-named + pool (≥ 10 nights) | 52 | 12 Aur, 28 Tau, HD 22298, HD 34959, HD 45910, HD 46264, HD 60848, HD 65079, HD 65875, HD 698, bet CMi, eta Tau, gam Cas, phi Per, pi. Aqr, psi Per, psi09 Aur, ups Cyg, zet Tau | 19 | 56 Eri, HD 23800, HD 26906, HD 37352, HD 44637, HD 45314, HD 58050, HD 6226, HD 62367 | 9 |

Full per-season table for the 52 non-core objects: `BeStar_Grism/notes/novelty/tables/pool.md`
(machine-readable: `tables/target_season_bess.csv`). Figures: `figures/coverage_pool_a.png`, `coverage_pool_b.png`.

Three facts about the pool:

- **The RLMT grism target list is the BeSS catalogue.** Our manifest names reproduce BeSS's own idiosyncratic names
  ("PHECDA", "PLEIONE", "ALCYONE", "53 Boo", "69 Ori", "2 Ori", "CQ UMa", "HD 79066" …). 34 non-core objects with
  ≥ 10 nights are SIMBAD Be*. The observers were running a BeSS-style Be patrol; the strategy's "core ten" kept the
  stars with the most frames, which are disproportionately the quiescent ones and the standards.
- **Not every famous Be star was active.** κ Dra (42 nights), V442 And = HD 6226 (47 nights), φ And, 4 Her, OT Gem and
  θ CrB were in absorption or nearly so in BeSS during our seasons. κ Dra and HD 6226 are well-sampled
  quiescent-phase records, not emitters.
- **Novelty is where BeSS is thin, not where the star is famous.** Ranked by our nights per BeSS in-season spectrum:

### Active seasons ranked by RLMT nights per BeSS in-season Hα spectrum

| Star [manifest name] | SpT | Season | RLMT nights (hrg/lrg) | BeSS Hα in season | ratio nights/(BeSS+1) | BeSS EW median (Å) | BeSS verdict | core ten? |
|---|---|---|---|---|---|---|---|---|
| HD 23800 [HD 23800] | B1IV | 2025-26 | 16 (16/0) | 0 | 16.0 | -1.93 | EMISSION(1) | no |
| HD 62367 [HD 62367] | B8V | 2024-25 | 13 (5/5) | 0 | 13.0 | -7.48 | EMISSION(1) | no |
| HD 26906 [V586 Per] | B7 | 2025-26 | 12 (12/0) | 0 | 12.0 | -6.39 | EMISSION(1) | no |
| ups Cyg [ups Cyg] | B2Vne | 2025-26 | 12 (12/1) | 0 | 12.0 | -33.24 | EMISSION | no |
| HD 46264 [QQ Gem] | Be | 2024-25 | 11 (8/7) | 0 | 11.0 | -5.27 | EMISSION | no |
| HD 37352 [HD 37352] | A0 | 2024-25 | 11 (9/9) | 0 | 11.0 | -6.75 | EMISSION(1) | no |
| HD 34959 [V1369 Ori] | B7Ib/II | 2025-26 | 30 (30/0) | 2 | 10.0 | -6.18 | EMISSION | no |
| 56 Eri [56 Eri] | B2(V)nne | 2025-26 | 10 (10/0) | 0 | 10.0 | -21.67 | EMISSION(1) | no |
| HD 46264 [QQ Gem] | Be | 2025-26 | 26 (26/0) | 2 | 8.7 | -4.44 | EMISSION | no |
| lam Eri [Lam Eri] | B2III(e)p | 2025-26 | 47 (46/1) | 5 | 7.8 | 2.86 | MIXED | yes |
| psi Per [Psi Per] | B5Ve | 2025-26 | 23 (23/0) | 2 | 7.7 | -38.57 | EMISSION | no |
| HD 44637 [HD 44637] | B3III:[n]e | 2025-26 | 15 (15/0) | 1 | 7.5 | -37.66 | EMISSION(1) | no |
| HD 65079 [BT CMi] | B2V(n)(e?) | 2024-25 | 15 (10/11) | 1 | 7.5 | -9.78 | EMISSION | no |
| HD 45910 [AX Mon] | B2IIIe | 2025-26 | 31 (31/0) | 4 | 6.2 | -40.07 | EMISSION | no |
| psi09 Aur [HD 50658] | B8IIIe | 2024-25 | 12 (9/9) | 1 | 6.0 | 0.07 | EMISSION | no |
| 12 Aur [12 Aur] | B2Ve | 2024-25 | 11 (9/9) | 1 | 5.5 | -6.88 | EMISSION | no |
| bet CMi [bet CMi] | B8Ve | 2024-25 | 15 (11/10) | 2 | 5.0 | -1.94 | EMISSION | no |
| HD 60848 [BN Gem; HD 60848] | O8:V:pe | 2025-26 | 25 (25/0) | 4 | 5.0 | -7.25 | EMISSION | no |
| HD 698 [V742 Cas] | B7:Ib-II(e) | 2025-26 | 19 (19/0) | 3 | 4.8 | -10.13 | EMISSION | no |
| 12 Aur [12 Aur] | B2Ve | 2025-26 | 12 (12/0) | 2 | 4.0 | -6.21 | EMISSION | no |
| HD 65875 [V695 Mon] | B2.5Ve | 2024-25 | 12 (6/6) | 2 | 4.0 | -40.29 | EMISSION | no |
| bet CMi [bet CMi] | B8Ve | 2025-26 | 28 (28/0) | 7 | 3.5 | -2.27 | EMISSION | no |
| HD 22298 [CT Cam] | B2Vne | 2025-26 | 20 (20/0) | 6 | 2.9 | -17.89 | EMISSION | no |
| phi Per [HD10516; Phi Per] | B1.5V:e-shell | 2025-26 | 14 (14/0) | 4 | 2.8 | -37.40 | EMISSION | no |
| eta Tau [Alcyone] | B7III | 2025-26 | 14 (14/0) | 5 | 2.3 | -1.19 | EMISSION | no |
| pi. Aqr [pi Aqr] | B1III-IVe | 2025-26 | 21 (21/0) | 9 | 2.1 | -6.91 | EMISSION | no |
| HD 58050 [OT Gem] | B2Ve | 2024-25 | 11 (5/5) | 7 | 1.4 | 3.11 | MIXED | no |
| HD 6226 [V442 And] | — | 2025-26 | 42 (42/0) | 32 | 1.3 | 2.70 | MIXED | no |
| gam Cas [Gamma Cas; gam cas] | B0.5IVpe | 2025-26 | 15 (15/0) | 11 | 1.2 | -39.09 | EMISSION | no |
| zet Tau [Zeta Tau; zet Tau] | B1IVe_shell | 2025-26 | 30 (30/0) | 24 | 1.2 | -18.25 | EMISSION | no |
| 28 Tau [BU Tau; PLEIONE] | B8Vne | 2025-26 | 11 (0/11) | 15 | 0.7 | -30.72 | EMISSION | no |
| HD 45314 [HD 45314; PZ Gem] | O9:npe | 2025-26 | 18 (18/0) | 35 | 0.5 | -1.84 | MIXED | no |

  Top of the list (QQ Gem, V1369 Ori, ψ Per, AX Mon, BN Gem, V742 Cas, β CMi, CT Cam, 12 Aur, HD 44637) is where
  an RLMT series out-samples BeSS by 3–10×. The bottom (γ Cas, ζ Tau, Pleione, HD 45314, π Aqr) is where BeSS already
  has near-weekly coverage: those are the **validation** stars for the "one-to-one against BeSS" figure, not science.

Caveat that limits all of this: in 2025-26 nearly everything is **hrg only**, so the strategy's hrg/lrg
dual-confirmation is unavailable for the second season, and the standards epoch (2025-12-05 →) covers only
2025-26. By the synthesis rule BE-S10, 2024-25 series are descriptive only.

---

## 6. BE-N1 decision (recorded)

**Rule (journal-editor):** < 2 verified-active emitters ⇒ the Be paper stops as a science paper; standards/precision
material moves to the instrument section.

**Applied to the paper as planned (core ten, ANALYSIS_STRATEGY §1–§3): 0 verified (1 weak) < 2 ⇒ STOP.**
No per-target pipeline effort, figures or draft on the core-ten science case. Its durable content is (a) the
standards/precision material (η Hya, θ Vir, Vega, Spica, plus six stars now shown to be *constant in Hα to the BeSS
≈ 0.4 Å floor*: 69 Ori, Phecda, φ Leo, 53 Boo, HD 70340, θ CrB — these are free null-test stars), and (b) the λ Eri
decay, one star, one weak event.

**Not decided by this package — needs the chair / James:** whether to re-draw the sample from the full campaign
before the gate is treated as final. If re-drawn with the same rule, the count is 19 (+ QQ Gem already inside it),
and the paper's ≥ 4 threshold for the stronger venue is met several times over. My recommendation: do not stop the
project; **stop the core-ten plan and re-issue Step 0 with the sample defined by this table** — verified-active stars
with ≥ 10 nights in 2025-26, ranked by the novelty ratio — keeping the six constant stars as nulls. That is a
different paper from the one the strategy describes (many stars × ~20–30 nights, hrg only, one standards-covered
season), and the ≤ 250-word abstract (standing rule 5) should be written for it before anything else.

---

## 7. Literature and sources

Verification levels: **[A]** I read the abstract/page in this session; **[B]** bibliographic entry taken from a
reference list or search-result abstract, primary paper not opened (publisher returned 403/405); nothing below is
cited from memory alone.

| Id | Reference | Used for | Level |
|---|---|---|---|
| L0 | Neiner, de Batz, Cochard, Floquet, Mekkas & Desnoux 2011, *AJ* 142, 149, doi:10.1088/0004-6256/142/5/149 — "The Be Star Spectra (BeSS) Database"; data via `http://basebe.obspm.fr/cgi-bin/ssapBE.pl` | all BeSS records and spectra | A |
| L1 | Balona & James 2002, *MNRAS* 332, 714, doi:10.1046/j.1365-8711.2002.05336.x — "Short-period line profile and light variations in the Be star λ Eridani" | λ Eri Be status | B |
| L2 | Rivinius, Baade & Štefl 2003, *A&A* 411, 229, doi:10.1051/0004-6361:20031285 — "Non-radially pulsating Be stars" | λ Eri periods 0.702 d and 0.269 d | B |
| L3 | Bossi et al. 1981, *A&AS* 46, 173 (bibcode 1981A&AS...46..173B) — "Spectroscopic and photometric observations of the Be star 69 Orionis" | 69 Ori Be history | B |
| L4 | Goraya & Tur 1996, *Ap&SS* 236, 175, doi:10.1007/BF00645142 — "Spectrophotometric Study of Four Bright Be Stars" | 69 Ori "normal star since Nov 1982" (as quoted by the Wikipedia article citing L3/L4) | B |
| L5 | Ghosh et al. 1999, *A&AS* 134, 359, doi:10.1051/aas:1999144 — "Observations of Bn and An stars: New Be stars" | 5 Cnc reported as Be with weak disk | B |
| L6 | Jaschek & Andrillat 1998, *A&AS* 130, 507, doi:10.1051/aas:1998101 — "Ae and A type shell stars" | Phecda as Ae star | B (PDF refused, 403) |
| L7 | Eiroa et al. 2016, *A&A* 594, L1, arXiv:1609.04263 — "Exocomet signatures around the A-shell star φ Leo?" | φ Leo variable Ca II K absorption | A (abstract) |
| L8 | Eiroa et al. 2021, *A&A* 653, A115, arXiv:2106.16229 — "The A-shell star φ Leo revisited" | φ Leo: δ Sct, edge-on variable disk; Hα circumstellar part variable in one run only | A (abstract) |
| L9 | Hauck & Jaschek 2000, *A&A* 354, 157 (bibcode 2000A&A...354..157H) — "A-shell stars in the Geneva system" | ν² Boo as A-shell star | B |
| L10 | Tkachenko et al. 2013, *A&A* (arXiv:1307.1970) — "Spectral modelling of the α Virginis (Spica) binary system"; Harrington et al. 2009 (arXiv:0908.3336) — "Line-profile variability from tidal flows in Alpha Virginis (Spica)" | Spica: SB2, P ≈ 4.0145 d, β Cep primary 0.1738 d | A (search abstracts) |
| L11 | Rivinius, Štefl & Baade 2006, *A&A* 459, 137, doi:10.1051/0004-6361:20053008 — "Bright Be-shell stars" | θ CrB shell-star context | B |
| L12 | Klement et al. 2024, *ApJ* (arXiv:2312.08252), doi:10.3847/1538-4357/ad13ec — "The CHARA Array interferometric program on the multiplicity of classical Be stars" | θ CrB Hα "dominated by photospheric absorption" | A (search abstract) |
| L13 | Sabogal, Ubaque, García-Varela, Álvarez & Salas 2017, *PASP* 129, 014203 — "Evidence of Dissipation of Circumstellar Disks from L-band Spectra of Bright Galactic Be Stars" | θ CrB: no Hα emission "for a long time", transient emissions since 1999 | A (search abstract) |
| L14 | Roark 1971, *AJ* 76, 634, doi:10.1086/111176 — "Photometric variability of the Be star theta Corona Borealis" | 0.7 mag fade in 1970 | B |
| L15 | SIMBAD via CDS Sesame (queried 2026-10-03; responses cached) | identities, spectral types, object types | A |

**Asked for and not found** (stated rather than filled in): published Hα EW time series for any core-ten star in
2024–26; any literature on HD 70340's emission; any dedicated paper on QQ Gem; a published period for θ CrB, 69 Ori,
5 Cnc; TESS-sector coincidences (not attempted — that is BE strategy Step "external", not this package). The search
for a 2025–26 λ Eri outburst report returned nothing: the October 2025 episode in BeSS appears to be unpublished.

---

## 8. Finding ids — status

| Id | Status | Why |
|---|---|---|
| **BE-S-1a-bess** (ED, RF; SYNTHESIS U9) | **CLOSED** | Per-target active/inactive table with BeSS spectra counts in and bracketing each season, measured EW, sources (§3, `tables/core_ten.md`). |
| **BE-N1-gate** (ED) | **CLOSED — decision recorded: STOP for the core-ten paper** | 0 verified-active (1 weak) < 2 (§6). Re-scope question raised for the chair/James; not a decision this package can make. |
| **QQ Gem disposition** (ED "ADD"; ROADMAP C3) | **CLOSED** | Sample extension; verified active both seasons (§4). |
| **θ CrB Be status / variability** (RF TCrB bullet; OA "DROP A4"; PH TCRB-A4) | **PARTIAL** | Be status and Hα state verified from 17 spectra; Hα constant. Continuum constancy in 2025–26 NOT verified (no photometry in scope); historical variability cited. |
| RF "BE-figures/BE-draft should not start until Step −1 returns ≥ 2 verified emitters" | **Condition evaluated: not met for core ten** | §3. |
| DS BE-S10 (pre-standards seasons descriptive) | not mine; **supported** | 2024-25 has lrg+hrg but no standards; 2025-26 has standards but hrg only (§5). |

## 9. Files changed (all inside the two directories this package owns)

- `BeStar_Grism/notes/novelty/bess_novelty_check.py` — the script (new)
- `BeStar_Grism/notes/novelty/test_bess_novelty_check.py` — 10 tests (new; lives here because the package does not own `pipeline/tests/`)
- `BeStar_Grism/notes/novelty/novelty.sqlite` — tables `nv_nights, nv_alias, nv_targets, nv_bess, nv_bess_done, nv_halpha, nv_meta` (gitignored by `*.sqlite`)
- `BeStar_Grism/notes/novelty/tables/` — `core_ten.md, brief_named.md, pool.md, be_n1_count.md, novelty_ranking.md, ew_trends.md, aliases.md, target_season_bess.csv`
- `BeStar_Grism/notes/novelty/figures/` — `coverage_core_ten.png, coverage_pool_a.png, coverage_pool_b.png`, 45 × `profiles_<star>.png`
- `BeStar_Grism/notes/novelty/cache/` — 1,266 raw Sesame/BeSS responses and spectra, 210 MB (audit trail; ignored via the local `.gitignore`)
- `BeStar_Grism/notes/novelty/.gitignore`
- `committee/work/novelty-be/REPORT.md` — this file

Nothing else was touched: no git state, no manifest write, no archive access (the script reads only the manifest's
`frames`, `stage_bestar_grism` and `build_meta`).

Reproduce: `/opt/miniconda3/envs/rlmt-checks/bin/python BeStar_Grism/notes/novelty/bess_novelty_check.py all`
(offline once `cache/` exists). The first run was made with `NOVELTY_WORKDIR` on the local SSD and synced back,
because cache writes to the shared spinning disk took seconds each while other packages were running.

## 10. Ledger changes requested of the chair

1. `BE-S-1a-bess` → **done**; evidence = this report + `tables/core_ten.md`.
2. `BE-N1-gate` → **done (decision: core-ten science paper stops)**; add a new open item **BE-N2-resample**
   ("re-draw the Be sample from the verified-active pool, or confirm the stop") owned by the chair/James.
3. Hold `BE-S-1b` (λ Eri injection), `BE-S11`, `BE-figures`, `BE-draft` until BE-N2 is decided; under a stop they are dropped.
4. `ROADMAP C3` (QQ Gem) → **closed: sample extension**. Correct its frame count: the Step-0 selection gives 141 grism
   frames / 37 nights (109 hrg, 22 lrg, 4 HaGrism, 6 OGGrism), not "122 hrg + 41 lrg".
5. ANALYSIS_STRATEGY §3.2 inventory corrections (I did not edit the strategy):
   - **Phecda is 366 frames / 47 nights, not 333 / 40**: 33 lrg frames on 7 nights in 2026 are filed under
     `gam UMa` and were never alias-merged (S0 `aliases` table gap — foundation-s0's file, not mine).
   - 69 Ori's "two seasons" are 2 nights (2024-25) + 49 nights (2025-26).
   - Other unmerged alias pairs in the manifest: Zeta Tau/zet Tau, BN Gem/HD 60848, Ry Sct/Ry Scuti, PZ Gem/HD 45314,
     Phi Per/HD10516, X Per/hd24534, BD+59 553/hd17520, Gamma Cas/gam cas, phi And/HD6811, BU Tau/PLEIONE
     (`tables/aliases.md`).
6. For the grism package: θ CrB shows no Hα emission in any of 17 independent spectra; the `g_extractions`
   "Hα emission" anchor on θ CrB is spurious. θ CrB's Hα EW = 5.4 ± 0.3 Å is a usable external check value for the
   fixed-dispersion extraction (G-1), as are 69 Ori 5.0, Phecda 8.7, φ Leo 6.8 Å (same windows as this script).
7. TCRB plan: A4 stays dropped (θ CrB continuum unverified).
8. Suggest the BeSS cache move to `external_data/bess/` (ROADMAP §3.3 convention) when someone who owns that tree can do it.

## 11. Needs James

- **Decide BE-N2**: stop the Be paper, or re-scope it to the ~20 verified-active stars (a different, broader, hrg-only,
  one-standards-season paper).
- **Independence check:** are any BeSS uploads in 2024–26 from RLMT/MACRO observers? Observer names in the bracketing
  spectra read as ARAS/BeSS amateur observers (`nv_halpha.observer`), but only James can
  confirm no consortium member is among them. If one is, those spectra are not independent witnesses.
- Who chose the "core ten", and was the BeSS catalogue the observing list? (It reads that way; if so the original
  observing proposal may already state the science case the strategy was missing.)

## 12. Limitations (so nobody over-reads this)

- BeSS spectra are heterogeneous (R ≈ 500–40,000; tellurics uncorrected; continuum placed by each observer). The
  0.25–0.45 Å inter-observer scatter is the floor; nothing smaller is claimed.
- The emission flag is "above the continuum". Weak infill is visible only in EW; for late-B/A stars with 5–9 Å of
  photospheric absorption, an infill below ≈ 1 Å (3× the floor) would be missed. So "inactive" means "no emission
  ≳ 1 Å in EW and none above the continuum", not "no circumstellar gas" (φ Leo and 53 Boo are shell stars).
- "UNWITNESSED" ≠ inactive: X Per, RY Sct and HD 17520 have no BeSS Hα spectra near our seasons. BeSS holds only
  stars in its own catalogue, so zero records does not mean "not a Be star".
- Thresholds (≥ 10 nights, ±60 d, 5%, 5σ, ≥ 2 witnesses) were fixed before measurement and not tuned; with a ±30 d
  pad or a 10% peak the core-ten count is unchanged (it is zero under any of them — only λ Eri's single 7% peak is
  threshold-sensitive, and it fails the ≥ 2-witness condition regardless).
- SIMBAD spectral types are quoted as returned; I did not adjudicate them.
- Six trend fits were run; λ Eri's 4.8σ is quoted with that trials factor in mind.
