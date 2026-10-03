**`committee/work/novelty-tcrb/REPORT.md` was not written.** The harness refused the write (subagents must return findings as text), and I did not route around it. The full report is below for the chair to save; the generated table, script, tests and figures are on disk.

## Summary (≤ 250 words)

TCRB-N1 is closed: the RLMT series is not the densest Hα record, and the editor's re-scope rule fires.

- **Density:** in RLMT's own window (2025-02-21 → 06-24 UT, 124 d) ARAS holds 188 Hα-covering spectra on 97 dates. RLMT has 247 labelled frames on 60 dates, 205 on 59 dates once pixel-verified. ARAS is denser in all five months.
- **Resolution:** ARAS has 51 dates at R ≥ 7,000. RLMT's delivered R is about 1,074 or about 318, depending on which way D1 falls. The rule fires on both branches.
- **Homogeneity:** RLMT's 59 dates beat the best single ARAS observer-and-setup (36 dates) and the best single observer (49). This is an advantage in count only; whether it gives smaller EW scatter is untested.
- **Venue:** PASP validation/methods note, or the D4 instrument paper's validation section. AJ only via a pre-registered scatter test.
- **Published series:** no Hα EW or flux series covering Feb–Jun 2025 exists, only isolated ATel points and Munari's profile atlas to 2025-05-01.
- **Status:** T CrB has not erupted as of 2026-10-03 and has been fading again through 2026.
- **Ephemeris and dip dates:** verified. The window spans 1.09 ellipsoidal cycles. The strategy's "post-dip recovery" framing is wrong: recovery began May 2024, and RLMT follows a second dip (Sep 2024 – Feb 2025).
- **R/I ellipsoidal amplitude:** not verified. No published value was found; it is bracketed at ±8–18% in flux.

Not done: AAVSO photometry (the site returned a bot challenge) and the AVSpec/BAA spectral databases (form-only search). Both can only add to the non-RLMT side.

---

# REPORT — package `novelty-tcrb` (TCRB-N1)

Date: 2026-10-03. Owns `TCrB_Monitoring/notes/novelty/` and `committee/work/novelty-tcrb/`. Findings: ED TCRB-N1; PH/OA "(verify)" items for T CrB.

All numbers in §1–§2 are emitted by `build_novelty_table.py` into `NOVELTY_TABLE.md`.

## 1. What was done

- **ARAS:** pulled the T CrB listing from https://aras-database.github.io/database/tcrb.html — 2,214 spectra, 2012-04-01 → 2026-09-25, pulled 2026-10-03T06:42 UT, sha256 `8d8d603b…`. It is cached as served and parsed to one row per spectrum (observer, site, R, λmin, λmax).
- **Truncation guard:** the parser checks its row count against the count the page declares. A first download truncated silently at 933 rows, which is why the check exists.
- **RLMT:** read canonical T CrB grism frames from the manifest read-only, joined to `frame_dispersion` for the S2c pixel verdict.
- **Unit of comparison:** UT calendar dates, not frames. RLMT takes 1–4 back-to-back frames a night; ARAS uploads one co-add.
- **Literature:** curated into `external/*.csv`, one source URL per row.
- **Tests and figures:** 12 tests pass; two figures.

Biases considered:

| risk | handling |
|---|---|
| frames ≠ epochs | compare UT dates; frame counts shown, never compared |
| FILTER label ≠ spectrum | two tiers: labelled, and S2c `dispersed`; both are upper bounds until G-3 runs (21 of 247 frames have header pointings > 1° off) |
| pooling 21 ARAS observers | every single-observer, single-setup series ranked separately |
| "covers Hα" is a choice | strict variant (6430–6700 Å): ARAS 141 spectra on 90 dates; conclusion unchanged |
| RLMT resolving power unknown (D1) | rule evaluated on both branches at the narrowest measured line width (13 px) |
| local night vs UT date | both binned in UT; RLMT's 60 local nights are also 60 UT dates |
| observer-code case variants (`XDU`/`XDu`) | upper-cased before ranking, tested |

## 2. Evidence

Figures: `fig_novelty_timeline.png` (who observed on which night; the one to look at) and `fig_novelty_monthly.png` (2023–2026 context).

### 2.1 The editor's table

| UT month 2025 | RLMT frames | RLMT dates (labelled / dispersed) | ARAS Hα spectra | ARAS dates | ARAS observers | ARAS dates low / mid / high R |
|---|---|---|---|---|---|---|
| Feb | 11 | 3 / 3 | 33 | 20 | 9 | 16 / 1 / 10 |
| Mar | 36 | 10 / 10 | 51 | 25 | 14 | 24 / 1 / 14 |
| Apr | 86 | 18 / 18 | 50 | 24 | 13 | 21 / 1 / 11 |
| May | 65 | 15 / 15 | 47 | 24 | 14 | 18 / 4 / 15 |
| Jun | 49 | 14 / 13 | 28 | 19 | 13 | 16 / 0 / 6 |

Classes: low R < 2,000; mid 2,000–7,000; high ≥ 7,000. These are whole months, so ARAS Feb and Jun include days outside the RLMT window.

### 2.2 Head-to-head inside the RLMT window

| series | spectra / frames | dates | median gap (d) | max gap (d) | median R |
|---|---|---|---|---|---|
| RLMT hrg+lrg, labelled | 247 | 60 | 1 | 8 | D1 pending |
| RLMT hrg+lrg, S2c dispersed | 205 | 59 | 1 | 8 | D1 pending |
| RLMT hrg / lrg, S2c dispersed | 91 / 114 | 51 / 48 | 1 / 1 | 9 / 10 | — |
| ARAS, all Hα | 188 | 97 | 1 | 4 | 1,000 |
| ARAS, low-R only | 117 | 83 | 1 | 5 | 568 |
| ARAS, high-R only | 65 | 51 | 2 | 8 | 13,490 |
| ARAS best single series (JRF/MTP-US, low-R) | 36 | 36 | 1 | 34 | 484 |
| ARAS best single observer (JRF/MTP-US, two setups) | 62 | 49 | 1 | 20 | 527 |

- **Strategy count:** "247 spectra on 60 nights" holds as a label count. Pixel-verified it is 205 frames on 59 dates; 7 frames are S2c `direct` and 35 `indeterminate`. Both grisms are dispersed on the same date on 40 dates.
- **High-R comparison:** the ARAS high-R subset alone (51 dates, max gap 8 d) matches RLMT's hrg series (51 dates, max gap 9 d).
- **Cross-validation sample for A5:** 46 dates have both an RLMT dispersed frame and an ARAS Hα spectrum. All 59 RLMT dates have an ARAS spectrum within ±1 d; 13 have none on the same date.
- **A0b:** the 2023–24 slot `6`/`W` frames number 24 on 14 dates, 16 dispersed on 11 dates. They are not externally unique: ARAS has 18 dates in 2023-05 and 18 in 2024-03.

### 2.3 The rule, clause by clause

| clause | D1 → 0.47 Å/px (OA.E1) | D1 → 1.59 Å/px (PH.P8) |
|---|---|---|
| RLMT R at 13 px line width | ≈ 1,074 | ≈ 318 |
| ARAS denser, window and every month, vs RLMT labelled | yes | yes |
| ARAS pooled median R above RLMT | no (1,000 vs 1,074) | yes |
| ARAS dates individually above RLMT's R | 59 | 97 |
| editor's rule fires | **yes** | **yes** |

The pooled median is a weak statistic here: ARAS is bimodal and many entries are a round 1,000. I therefore judged the rule on the count of ARAS dates that individually out-resolve RLMT, with a threshold of more than half of RLMT's dates. That threshold is my judgement; the raw counts are printed so the chair can draw the line elsewhere. On the PH branch no placement changes the outcome.

### 2.4 Campaigns without a public log, and published Hα series

| campaign | R | coverage | source |
|---|---|---|---|
| Asiago/ANS (Munari et al. 2025, A&A 701, A176) | low-res B&C; echelles 12,000–27,800 | Hα profile at ~1-week cadence 2023-06-01 → 2025-05-01; data closed 2025-05-03; profile atlas, no Hα EW table in the text | https://arxiv.org/abs/2507.23323 |
| Tautenburg 2 m (ATel 17030) | 65,000 | "more or less regularly" for a year; EW doubled 2025-01-21 → 02-09 | https://www.astronomerstelegram.org/?read=17030 |
| Rozhen 2 m (ATel 17075) | echelle | EW 12.5 ± 0.6 Å (2025-01-18), 21.8 ± 1.0 Å (2025-03-09) | https://www.astronomerstelegram.org/?read=17075 |
| UCF 0.5 m (ATel 17110) | not stated | EW 10.2 / 5.5 / 13.1 Å on 2025-03-03 / 07 / 19 | https://www.astronomerstelegram.org/?read=17110 |
| Habtie et al. (ATel 17041) | 13,500 | Hα flux ×2.3 between 2025-01-15 and 02-07 | https://www.astronomerstelegram.org/?read=17041 |
| Mercator/HERMES (Planquart et al. 2025) | 86,000 | 100 spectra 2011-01 → 2023-06; ends before the window | https://arxiv.org/abs/2501.02984 |

No published Hα EW or flux series covering Feb–Jun 2025 was found. I checked arXiv:2512.19218, 2601.16190, 2504.20592, 2607.15245, 2605.20991 and 2604.02708; none contains an optical Hα series. The uniform measured series is therefore unpublished, but denser and better-resolved raw material for one is already public in ARAS.

Warning for A5: the telegram values differ by a factor of four within two days (5.5 Å on 03-07, 21.8 Å on 03-09). The "10–15% agreement" criterion cannot be tested against them. It needs same-date ARAS spectra re-measured with RLMT's own windows.

## 3. The "(verify)" items

### 3.1 Erupted? No, as of 2026-10-03

- ARAS spectra continue to 2026-09-25 under the banner "one spectrum/day until the next nova outburst".
- No T CrB telegram appears among ATel #18029–#18086 (30 days to 2026-10-02).
- The latest T CrB telegram found, ATel 17784 (2026-05-09), gives U = 11.85, B = 11.61, V = 10.25: about 1 mag fainter in U and 0.5 mag in B and V than 2024 (https://www.astronomerstelegram.org/?read=17784).
- arXiv:2609.28191 (2026-09-23) and 2609.06201 still treat the eruption as anticipated.
- Pei et al. 2026 condition on "no eruption by 2026 July 11" and report a renewed decline in 2026 (https://arxiv.org/abs/2607.15245).
- Wikipedia: "no nova has been observed as of October 2026".

The AAVSO light curve was not checked directly. The strategy's "recovered from the dip" is stale.

### 3.2 Ellipsoidal R/I amplitude — not verified; bracketed

No modern published R- or I-band amplitude was found. Published full amplitudes:

| band | full amplitude (mag) | source |
|---|---|---|
| V | ~0.3 | Schaefer 2023 (https://arxiv.org/abs/2303.04933); Merc et al. 2025 (https://arxiv.org/abs/2504.20592) |
| B, V | 0.35 (semi-amplitudes 0.172, 0.173) | Iłkiewicz et al. 2023 (https://arxiv.org/abs/2307.13838) |
| V | ~0.4 | ATel 16107; Zamanov et al. 2025 (https://arxiv.org/abs/2504.06029) |
| V | ~0.5 incl. hot spot | Munari et al. 2025 |
| I | ~0.35 observed, more than twice that expected for the giant alone | Shahbaz et al. 1997 (https://arxiv.org/abs/astro-ph/9703146) |
| J, K | 0.18, 0.15 | Maslennikova et al., via Hinkle et al. 2025 (https://arxiv.org/abs/2502.20664) |

Flux semi-amplitude is about 0.46 × the full amplitude in mag, so the continuum at 6563 Å swings between ±8% and ±18%. PH's ±10% sits at the low end of that bracket. TCRB-A5b must measure the amplitude from AAVSO Rc/Ic photometry.

### 3.3 Ephemeris — verified; three agree

| ephemeris | P (d) | T0 (HJD, giant at max velocity) | window phase |
|---|---|---|---|
| Munari et al. 2025 | 227.5528 ± 0.0002 | 2459978.37 ± 0.08 | 0.292 → 0.837 |
| Hinkle et al. 2025 | 227.5494 ± 0.0049 | 2455427.51 ± 0.10 | 0.292 → 0.837 |
| Fekel et al. 2000 (via Zamanov et al. 2025) | 227.5687 ± 0.0099 | 2447918.62 | 0.286 → 0.831 |

The window spans 1.09 ellipsoidal cycles, confirming PH. It runs from just after one minimum (phase 0.25) through maximum (0.5) to past the next minimum (0.75). I recommend adopting Munari et al. 2025.

### 3.4 Dip and super-active-phase dates — verified, framing corrected

| event | date | source |
|---|---|---|
| super-active phase ends | late Apr / early May 2023 (stepped drop JD 2460043–50) | Munari et al. 2025; ATel 16109; arXiv:2307.00255 |
| dip begins | Mar/Apr 2023 | ATel 16107 |
| first minimum | late Aug 2023 | Munari et al. 2025; Pei et al. 2025 |
| rebound; second minimum | mid-Jan 2024; late Mar 2024 | Munari et al. 2025 |
| recovery begins | May 2024 | Munari et al. 2025 |
| flare | He II peak 2024-11-09, B peak 2024-11-16 | Munari et al. 2025; ATel 16912 |
| second, shallower dip | Sep 2024 → Feb 2025 | Pei et al. 2025 (https://arxiv.org/abs/2512.19218) |
| Hα doubles | 2025-01-21 → 02-09 | ATel 17030, 17041 |
| renewed decline | 2026 | Pei et al. 2026; ATel 17784 |

RLMT does not cover "the post-dip recovery". It starts about two weeks after the February 2025 Hα jump that ended the second dip. RLMT's March 2024 `6`/`W` frames do fall on the second deep minimum.

### 3.5 Reference errors in `TCrB_Monitoring/ANALYSIS_STRATEGY.md` (not edited)

- arXiv:2405.11506 is Zamanov et al. 2024, not Munari.
- "Munari's convention" for Hα EW windows has no referent that I could find. Munari 2023 measures Hβ, He I and He II fluxes; Munari et al. 2025 gives an Hα profile atlas.
- The "A&A 2026 HST/Swift/NuSTAR/XMM" paper is Luna et al. (arXiv:2601.16190). I could not confirm the number `aa57435-25`.
- Confirmed as cited: arXiv:2303.04933, 2307.00255, 2507.23323, 2512.19218, 2504.20592, 2502.20664, 2607.05200, and the RNAAS DOI.

## 4. Findings

| id | status | why |
|---|---|---|
| ED TCRB-N1 | **CLOSED** | table on disk, generated and tested; acceptance clause triggered, re-scope stated |
| PH R-band amplitude (A5b) | **PARTIAL** | no published value; bracket ±8–18%; needs AAVSO Rc/Ic |
| PH 1.09 ellipsoidal cycles | **CLOSED** | reproduced on three ephemerides |
| target status | **CLOSED for 2026-10-03** | not erupted; re-verify at submission |
| dip / SAP dates | **CLOSED** | §3.4; "post-dip recovery" framing should change |
| OA A0b uniqueness | **REBUTTED** as a novelty claim | ARAS has 18 dates in each month; the task is still worthwhile |
| "densest homogeneous single-instrument" | **REBUTTED** on densest; **PARTIAL** on homogeneous | 59 vs 97 dates; 59 vs 36 for any single ARAS setup |

## 5. Venue

The T CrB paper is a PASP validation/methods note, or the D4 instrument paper's validation section. "ApJ-grade. The paper's spine" should be struck from strategy Rank 1.

AJ is reachable only via a test fixed before A5 runs:

1. D1 returns ≥ 3 lines with residual < 1 px, and G-3 leaves ≥ 50 dates.
2. EW is measured with identical windows on RLMT and on ARAS low-R spectra degraded to RLMT's measured LSF.
3. Both are detrended with the same smoother and residual scatter is compared by whole-night bootstrap.
4. AJ only if RLMT's scatter is smaller with a 95% interval excluding equality, and RLMT shows a feature the pooled ARAS series does not. The scatter ratio and the detectable effect size are reported either way.

## 6. Files changed

All are new, all under `TCrB_Monitoring/notes/novelty/`. No git state was touched and nothing outside the package was edited.

- `build_novelty_table.py`, `test_novelty_table.py`
- `NOVELTY_TABLE.md`, `novelty_monthly.csv`, `novelty_window.csv`, `aras_series_window.csv`, `aras_observers_window.csv`, `rlmt_grism_dates.csv`, `novelty_verdict.json`
- `fig_novelty_timeline.{png,pdf}`, `fig_novelty_monthly.{png,pdf}`
- `external/aras_tcrb.html.gz`, `aras_tcrb_pull.json`, `aras_tcrb_listing.csv`, `literature_campaigns.csv`, `published_halpha_points.csv`, `ephemerides.csv`

## 7. Limits

- **AVSpec and BAA not counted.** Both spectral databases are searchable only through web forms, which I did not submit. They can only add to the non-RLMT side.
- **AAVSO photometry not pulled.** The data endpoint returned a bot challenge and I did not attempt to get around it.
- **Asiago, Tautenburg and Rozhen** are represented by what their papers and telegrams state, not by counts.
- **ARAS R and λ ranges** are taken from the listing; no spectrum was opened.
- **Not read:** Stoyanov et al. 2025 (Bulg. Astron. J. 42, 29) and the Teyssier RNAAS flare note.

## 8. Ledger changes requested

1. `TCRB-N1-novelty-table` → done; evidence `TCrB_Monitoring/notes/novelty/NOVELTY_TABLE.md`.
2. Venue ApJ → PASP note or D4 validation section. Add `TCRB-N3-homogeneity-test` (§5) as the only AJ path.
3. Strategy text:
   - strike "densest";
   - replace "post-dip recovery" with the re-brightened state after the Sep 2024 – Feb 2025 second dip;
   - replace "247 spectra on 60 nights" with ≤ 205 pixel-verified frames on ≤ 59 dates;
   - fix the three references in §3.5.
4. `TCRB-A5`: cross-validate against same-date ARAS spectra (46 dates) re-measured with RLMT windows, and define the EW windows explicitly.
5. `TCRB-A5b`: measure the amplitude from AAVSO Rc/Ic; use the ±8–18% bracket until then; adopt the Munari et al. 2025 ephemeris.
6. `TCRB-A0b`: keep, but drop the uniqueness rationale.
7. New small task: count AVSpec and BAA T CrB spectra for Feb–Jun 2025.
8. `TCRB-N2`: note the 2026 fading, and that ARAS will supply the pre-eruption spectroscopic baseline regardless.

## 9. Needs James

- Download from AAVSO by browser: T CrB Rc/Ic/B/V photometry 2005–2026 and the AVSpec T CrB listing.
- D4: this result strengthens the case for T CrB as the instrument paper's validation section.
- Consider contacting F. Teyssier (ARAS). JRF/MTP-US's 49-date series is the natural comparison partner.
- Re-verify eruption status on the day of submission.