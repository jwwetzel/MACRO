**Summary (≤250 words)**

`REPORT.md` is not on disk: the harness refused the write (sub-agents may not write report files). Its full text is below for the chair to save at `committee/work/cv-literature/REPORT.md`. All other deliverables are written.

- **Bibliography:** `references.bib` grew from 10 to 113 entries. `verify_bib.py` checks each against Crossref, DataCite, arXiv, ADS scans or the publisher page: 112 verified, 1 manual (the AAVSO citation template), no mismatches. Unmodified `main.tex` compiles clean against it in a scratch copy.
- **Comparison subsection:** `comparison_section.tex` cites 62 keys, with no hand-typed numbers. Literature values come from `literature_scales.py`; 12 placeholders name the revision package that owns each.
- **Verify items:** PH.P2, P3, P5, P6 and RF minor 1 are closed in `verify_items.md`. CV-R12 is partial until `main.tex` takes `citation_map.md` (60 rows).

Findings that change the paper:
- The existing `ferrario1993` entry had a wrong title (fixed).
- `\software{}` credits photutils, which the pipeline never imports; photometry is `sep`.
- The VSX ST LMi ephemeris is not published: "from AAVSO data", November 2025, no error, fiducial undocumented.
- Cropper (1986) supplies the only published period error; cycle-count drift becomes 0.022 cycles, not 0.0014 (still unique).
- Bailey et al. (1985) measured a wavelength-dependent bright phase in ST LMi, predicting negative g−i, the sign of the committee's −110 s.
- The novelty claim is too strong: simultaneous multiband single-night curves already exist.
- YZ Cnc: superoutburst amplitude is 4.0 mag, not 3.0; a more precise period and a 0.5 mag prior hump detection exist.

Some values are second-hand or summary-grade; these are flagged in report §4.

---

**Full text for `committee/work/cv-literature/REPORT.md`**

# Package `cv-literature` — report

Date 2026-10-03. Findings owned: **CV-R12** (ED.E1, RF.M8) and the "(verify)" items **PH.P2, PH.P3, PH.P5, PH.P6, RF minor 1**. Files owned: `manuscripts/CV_TimeSeries/references.bib` and `committee/work/cv-literature/`. `main.tex` was not edited. No git state was touched.

## 1. What was done

1. **Bibliography: 10 → 113 entries, every one checked against its registry by a script.** `verify_bib.py` parses `references.bib` and compares each entry with its registry of record:
   - DOI → Crossref or DataCite: first author, year ±1, volume, first page, title.
   - e-print → arXiv API: first author and title.
   - bibcode without DOI → ADS scan service (a non-existent bibcode returns 403).
   - otherwise → the publisher's page.

   Result (`bib_verification.md`): **112 VERIFIED, 1 MANUAL, 0 mismatches, 0 duplicate keys**. Six offline tests (`test_verify_bib.py`) pin the parser and the comparison logic.
2. **`comparison_section.tex`** — a drafted "Comparison with previous work" subsection: five paragraphs for ST LMi, four for YZ Cnc, and an ephemeris-provenance table for all five targets. It cites **62** distinct keys. No number is typed: ours are `\Num…` macros, literature values are `\Lit…` macros, and values the revision packages have yet to produce are `\PH{…}` placeholders naming their owner (12 placeholders: CV-R1, R2, R3/R5, R4, R6, D3, cv-stats).
3. **`literature_scales.py`** — the record behind every literature number: 83 constants, each with BibTeX key, locator and a provenance grade (*direct / quoted / catalogue / summary*), and 43 derived scales. It emits `literature_macros.tex` (93 macros, in the style of `numbers.tex`) and `literature_scales.md`. It reads six of our macros read-only from `numbers.tex` so the derived scales track the pipeline.
4. **`verify_items.md`** — each "(verify)" answered with the cited value and how it was read (V1–V7), plus two items that surfaced on the way (V8, V9).
5. **`citation_map.md`** — 60 rows: for each uncited use (file:line) the key(s), and the correction where the sentence is wrong once the literature is read.

Compile evidence (tectonic, in the scratchpad, repo untouched):
- Unmodified `main.tex` against the new bib: 0 BibTeX warnings, 0 over/underfull boxes, 0 undefined citations.
- `\nocite{*}` of all entries: 113 `\bibitem`s, 0 warnings.
- `comparison_section.tex` with `numbers.tex` and `literature_macros.tex`: 0 undefined control sequences, 0 undefined citations, 62 `\bibitem`s.

## 2. Evidence

### 2.1 Verification table (summary; full table in `bib_verification.md`)

| Method | Entries | What it proves |
|---|---|---|
| Crossref (DOI) | 101 | metadata agree field by field |
| DataCite (Zenodo concept DOI) | 1 (`photutils`) | resolves, first author agrees; year/title not comparable for a concept DOI |
| arXiv API | 6 | first author and title agree |
| ADS scan served | 3 (`watson2006`, `bonnetbidaud1996`, `campbell1999`) | the bibcode exists; title/authors were then read off the scan by eye |
| publisher page | 1 (`walker1965`) | Konkoly archive page shows author and title |
| MANUAL | 1 (`aavso_aid`) | AAVSO's own citation template; aavso.org refuses scripted access |

Reproduce: `/opt/miniconda3/envs/rlmt-checks/bin/python committee/work/cv-literature/verify_bib.py` (cache: `bib_verification_cache.json`; `--refresh` to re-query).

The check caught three real errors:
- Two DOIs I recalled from memory for the Kato et al. survey papers resolved to unrelated PASJ papers.
- The **pre-existing** entry `ferrario1993` carried a title that is not the paper's ("Cyclotron emission from accreting magnetic white dwarfs"). The DOI's registered title is "Detection of cyclotron emission features in the infrared spectrum of ST LMi". Fixed. That paper is the ST LMi field-strength paper, not a general cyclotron reference, which changes how `main.tex:158` should use it.

### 2.2 What the literature says that the paper needs (numbers from `literature_scales.md`)

| # | Finding | Evidence | Bears on |
|---|---|---|---|
| L1 | **The ST LMi "catalogue ephemeris" is not a published ephemeris.** VSX's period and epoch were set on 2025-11-12 by the VSX moderator "from AAVSO data": no paper, no error, fiducial undocumented. AN UMa's epoch likewise. | VSX detail pages (OID 17253, 37171), VizieR `B/vsx` | `main.tex:557`, 684, 727; RF minor 1; PH.P4 |
| L2 | The only published σ_P is Cropper (1986): 0.07908908 ± 0.00000008 d. Drift over 21,869 cycles = **0.022 cycles** (paper: 0.0014). Count still unique, by 23×. Our refit agrees with Cropper at 0.6σ. | `cropper1986` §3, read on the ADS scan | RF minor 1 — closed |
| L3 | **Cropper's phase zero is the feature we time**: the linear-polarisation pulse "half way through the intensity's rapid decline". VSX zero falls at polarimetric phase 0.97 ± 0.18; our edge at 0.13 ± 0.18. The 1,071 s offset is consistent with the distance from bright-phase maximum to its end. | `cropper1986` §4; arithmetic on catalogue constants | PH.P4; CV-R3 |
| L4 | **A wavelength-dependent bright phase has been measured in ST LMi**: 0.29 P (white), 0.31 (J), 0.36 (H), 0.39 (K), simultaneous. Redder ends later: 68 s (white→J), 273 s (J→K); ≈ 37 s interpolated to g→i (order of magnitude). Predicted sign of g − i: **negative**. | `bailey1985` §7.3 | PH.P3 — closed; **D3 / CV-R1**: U3's −110 ± 28 s has the literature's sign |
| L5 | The bright phase changed length by 0.05 P twice in the record (1982→85; within 1982–83): one edge moves **171 s = 9°**. Our rms is 84 s = 4.4°. | `cropper1986` §4, `bailey1985` §8 | PH.P2; CV-R5 (state split now has a scale to beat) |
| L6 | Ṗ scales: secular 4×10⁻¹⁴–2×10⁻¹³; V1500 Cyg 3.86×10⁻⁸ (11× above our bound — excluded); DP Leo-like libration ≲ 4×10⁻¹¹ (100× below). | `schmidt1995` via `pavlenko2018`; `beuermann2014`; `knigge2011` | PH.P2 — closed |
| L7 | B = 12.1 ± 0.5 MG; g, r, i sit at cyclotron harmonics ≈ 18, 14, 12. Bright phase is red; amplitude 2.0 mag in R_c, ≈ 1.2 in V. | `campbell2008c`, `peacock1992`, `kafka2007` | PH.P6 — closed; CV-R6 |
| L8 | **The paper's novelty sentence is too strong.** Simultaneous multiband single-night light curves of ST LMi exist (Bailey et al. 1985; Peacock et al. 1992, three occasions, with polarimetry). What is new is nights × seasons × an independent state tag. | `peacock1992` abstract, `bailey1985` | `main.tex:162–176`, 1160–1163; CV-R13; ED.E2 |
| L9 | YZ Cnc superoutburst amplitude is **4.0 mag** (15.0→11.0), normal ≈ 3.0. `\NumSuperoutburstAmp` = 3.0 is a normal-outburst amplitude. | `kato2002`, `hakala2004`, `vanparadijs1994`, `patterson1979` | PH.P5 — closed |
| L10 | A YZ Cnc period 30× more precise, with an epoch, exists: 0.086924(7) d, HJD 2447518.7255(61). "Ephemeris drift below 0.01 cycles" (`main.tex:1121`) is 0.027 under the σ_P of the period used. Conclusion survives; number does not. | `vanparadijs1994` §6 | CV-R4 |
| L11 | **Prior detection of the YZ Cnc orbital hump: 0.5 mag full amplitude in quiescence**, maximum at phase 0.8, moving to phase 0.5 at late decline. Ours: 30–70 mmag semi-amplitude, 4–8× smaller. TESS: no orbital signal. | `vanparadijs1994` §6–7; `sun2026` | RF.M7; CV-R4 |
| L12 | Superhump peak semi-amplitude is 125–150 mmag (0.25–0.3 mag full); in YZ Cnc late superhumps survived into the next normal outburst (2011). The 50 mmag "floor" has no citation yet. | `smak2010`, `kato2012`, `dai2026`, `kato2014v` | PH.P1; CV-R4 |
| L13 | Two September-2026 preprints analyse YZ Cnc's TESS superoutbursts (Dai et al. arXiv:2609.06566; Sun et al. arXiv:2609.12461). Our superhump null must be framed against them. | arXiv API | CV-R4; ED |
| L14 | `\software{}` credits **photutils, which the pipeline does not import and the environment does not contain**; the photometry is `sep`. `astroalign`, `astroquery`, `pandas`, Gaia are used and uncredited. | `grep` over `pipeline/`; `import` test in `rlmt-checks` | `main.tex:1565–1571`; CV-R9 |

Circularity check (is the VSX period independent of our data?), query `SELECT target, source, n_points, n_independent, notes FROM cv_external` on `products/phot/cv_timeseries.sqlite` (opened `mode=ro`): ST LMi AAVSO has 9,014 rows, 0 flagged as RLMT resubmissions (YZ Cnc: 1,499 of 99,732).

### 2.3 Standing statistical rules (SYNTHESIS §5)

- Rule 2 (every null carries a predicted scale): supplied for all four of the paper's nulls — Ṗ (L6), band offset (L4), superhump (L12), orbital hump (L11).
- Rules 1 and 3 concern estimators this package does not run. The placeholders in the draft section ask CV-R1/R2/R3 for signed, scatter-based values and do not pre-judge D3.
- No fit, significance test or selection was performed here. The only arithmetic is in `literature_scales.py::derive()`, each line of which states its formula.

## 3. Finding status

| Finding | Status | Why |
|---|---|---|
| **CV-R12** (ED.E1, RF.M8) | **PARTIAL** | Everything this package can deliver is delivered and verified. ED.E1's closing condition also requires the citations and subsection to be *in* `main.tex` with `main.log` clean; that edit belongs to CV-R13. Compile-clean is demonstrated on a scratch copy. |
| **PH.P2** "(verify)" | **CLOSED** | V1500 Cyg rate confirmed (3.86×10⁻⁸); three scales tabulated; longitude scale from the star's own history added. The value was read in Pavlenko et al. (2018), not in Schmidt et al. (1995). |
| **PH.P3** | **CLOSED** | Expected scale and sign supplied from a measurement in this star (L4). |
| **PH.P5** "(verify)" | **CLOSED** | Confirmed; the constant is wrong and needs an emitter change (not mine). |
| **PH.P6** "(verify)" | **CLOSED** | 12.1 ± 0.5 MG; harmonic numbers computed; physical paragraph drafted. |
| **RF minor 1** | **CLOSED** | Published σ_P supplied; drift recomputed; wording drafted. |
| PH.P4 (not owned; brief item 1 covers its open half) | **PARTIAL** | Fiducials identified for VV Pup, EU UMa, AN UMa, YZ Cnc. ST LMi's VSX fiducial is undocumented and can only be closed by asking VSX. |

Nothing is REBUTTED. Nothing is BLOCKED.

## 4. Limits of this work

- **Second-hand values**, flagged *quoted* in `literature_scales.md`:
  - Schmidt et al. (1995) Ṗ, via Pavlenko et al. (2018).
  - Stockman et al. (1983) period, via Cropper (1986).
  - Ferrario et al. (1993) and Schmidt et al. (1983) field strengths, via Kafka et al. (2007), Campbell et al. (2008) and Cropper (1986).
  - Moffett & Barnes (1974) flickering amplitude.
  - Shafter & Hessman (1988) period.
  - Walker (1965) ephemeris, via VSX and Howell et al. (2006).
- **Summary-grade values** (from a fetched-page summary, not the text): the Robertson et al. (2008) ephemeris and hump phases, and Potter (2000)'s ±15–18° extent. Both are secondary comparisons; confirm at source before submission.
- Knigge et al. (2011) is cited for the gigayear timescale; I did not read a number off it. The secular Ṗ range is P/τ arithmetic.
- EU UMa: Mittaz et al. (1992) and Howell et al. (1995) are cited as the published studies; their period values were not re-read.
- The g→i chromatic estimate (37 s) is a log-linear interpolation of infrared durations into the optical, labelled order-of-magnitude wherever it appears.
- ADS-scan verification proves a bibcode exists, not that the metadata are right; for those three entries the metadata were read off the scan by eye.
- A late-plateau superhump amplitude was not read off Kato et al. (2012)'s figures.
- Not searched exhaustively: X-ray literature of ST LMi, eROSITA population papers, AN UMa polarimetry history. The strategy's seed list is covered except "OY Car ApJ study" and eROSITA.

## 5. Files changed

| File | Change |
|---|---|
| `manuscripts/CV_TimeSeries/references.bib` | 10 → 113 entries. The ten original entries are byte-identical except `ferrario1993`'s title (corrected). `manuscripts/` is git-ignored, so the pre-edit file is preserved as `committee/work/cv-literature/references.bib.before-cv-literature`. |
| `committee/work/cv-literature/verify_bib.py`, `test_verify_bib.py`, `literature_scales.py` | new |
| `bib_verification.md`, `bib_verification_cache.json`, `literature_macros.tex`, `literature_scales.md` | generated |
| `comparison_section.tex`, `verify_items.md`, `citation_map.md` | new |

Nothing else in the repository was written. Databases were opened read-only.

## 6. Changes wanted elsewhere (not made)

1. **`main.tex` (CV-R13 owner):**
   - Apply `citation_map.md`.
   - Add `\input{literature_macros}` beside `\input{numbers}`.
   - Place `comparison_section.tex` after "What these data add".
   - Reword ll. 162–176 and 1160–1163 (L8), 557–559 (L1), 727–729 (L2), 864–867 (L6), 1068–1070 (L9), 1097 and 1121 (L10), 1140–1144 (L12).
   - Fix `\software{}` and `\facilities{}` (L14).
2. **`numbers_cv.py` / `final_science.py` (emitter owners):**
   - `\NumSuperoutburstAmp` 3.0 → 4.0 with source `kato2002`.
   - Give `\NumSuperhumpFloorMmag` a named source and add the peak value.
   - Make the VSX provenance notes say "from AAVSO data, no published error" for ST LMi.
   - Use the `cropper1986` σ_P drift instead of the assumed one.
3. **CV-R4:**
   - Refold YZ Cnc on 0.086924 d and restate the between-night drift.
   - Count the run-filters whose contour lies below 125 mmag and the scopes whose red-noise contour lies below 250 mmag.
   - State where each outburst run sits relative to the preceding superoutburst.
4. **CV-R1/R2 (D3):** report g − i with its sign against the predicted "negative, tens of seconds to ~10² s".
5. **CV-R5:** report the high-/low-state edge-longitude difference against the 9° scale.
6. **CV-R3 (optional, pre-register before looking):** an O−C point at Cropper's 1982 epoch — same feature, 44-yr lever arm, |Ṗ| sensitivity ~10⁻¹¹. It is conditional on the cycle count (0.18 cycle at 1σ under Cropper's σ_P: unique at 1σ, not at 3σ). Candidate, not a result.
7. **Captions (`captions.tex:6, 9, 13, 14, 15`) and `tables.tex:81`:** survey and catalogue citations; verdict-row wording.
8. **Package `clock`:** supply the citation for the eclipsing binary's ephemeris (`main.tex:331–336`).

## 7. Ledger changes requested of the chair

- **CV-R12-literature:** evidence = the four deliverables in `committee/work/cv-literature/`; status *deliverables done; closes with CV-R13*.
- Add **CV-R12b-software-credits** (blocking, trivial): L14.
- Add **CV-R4** sub-item: van Paradijs et al. (1994) period/epoch and prior hump detection (L10, L11); frame against Dai et al. and Sun et al. (L13).
- Add to **CV-R5**: the 9° scale (L5) as the predicted scale for the state-split longitude.
- Add to **CV-R13**: novelty sentence (L8); "published orbital ephemerides" (L1).
- Record in **D3**: published prediction of the sign (L4).
- Amend **RF minor 1**: σ_P now sourced (L2).
- **CV-R15:** `verify_bib.py` exit status as a pre-submission gate.

## 8. Needs James

1. **One e-mail to VSX (Sebastián Otero):** what do the ST LMi and AN UMa epochs of the 2025-11-12 revisions mark, what data and method produced the ST LMi period, and is there an uncertainty? (L1)
2. **Decide the novelty wording** (L8): the abstract and introduction claim something Peacock et al. (1992) and Bailey et al. (1985) already did in part.
3. **Facility acknowledgments:** ZTF, ASAS-SN, Gaia/DPAC, VSX and SIMBAD/VizieR each prescribe wording, to be copied from their pages at submission. The year in `aavso_aid` must match the download year.
4. **Whether to pursue the 44-yr O−C point** (§6 item 6); it needs a pre-registered rule first.
5. **Instrument citation** (SYNTHESIS D4): nothing exists to cite for the telescope at `main.tex:223`.
6. `manuscripts/` is git-ignored: the bibliography now has 113 verified entries and no version history.