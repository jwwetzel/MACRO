**Summary (≤250 words)**

Package `ops` is delivered except for one file: the harness refused to let me write `committee/work/ops/REPORT.md`, so its full text is below for the chair to save. Nothing was sent and git was not touched.

- **Rev. 3** (`ops/2026-10_observatory_request_rev3.md`, marked DRAFT) is written for the QHY600 under pyscope, and its tables are script-generated.
- **Visibility** reproduces OA's table to 0.1 h. T CrB has a usable snapshot until 2026-10-23, then from 2026-12-16; a ≥2 h run first fits 2027-01-09. ST LMi reaches ≥1.5 h on 2026-10-30.
- **T CrB exposure:** 240 s is kept. The predicted brightest pixel on the QHY is 4–8% of the saturation veto, from a QHY/ASI bridge of about 0.8 measured on shared standards.
- **Imaging:** 5 × 10 s g, r, i, mildly defocused, with an in-focus fallback.
- **Calibration:** 29,148 QHY frames in the archive, zero raw calibration frames.
- **Eruption block** (`ops/eruption_block/`): three `.sch` files, CSV, README and telegram skeleton. The `.sch` syntax is not validated; pyscope is not installed here and its source could not be fetched.
- **`AUTHORSHIP.md`** and three unsent email drafts are in place.
- **Tests:** 59 pass.

New facts: the wheel has no B filter; QHY master calibrations exist on TCC21 but were never synced; filter focus offsets are zero under pyscope.

Deviations: rev. 2 is marked superseded in a new `ops/README.md`, not edited, because it is a provenance-digested output. I added `ops_exposure.py` beyond the brief.

Status: U10, OA, DE and TE.F7 items closed; eruption block, TCRB-N2 and ED.E4 partial.

---

**Full text for `committee/work/ops/REPORT.md`**

# Package `ops` — report

Date 2026-10-03. Nothing was sent to anyone. No git state was touched. The archive and the manifest were opened read-only.

## 1. What was delivered

All paths are under `/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO/`.

| # | Deliverable | Path | State |
|---|---|---|---|
| 1 | Observatory request rev. 3 | `ops/2026-10_observatory_request_rev3.md` | complete, marked DRAFT, not sent |
| 1a | Visibility generator and tests | `pipeline/scripts/ops_visibility.py`, `pipeline/tests/test_ops_visibility.py` | 34 tests passing |
| 1b | Exposure and calibration-census generator and tests | `pipeline/scripts/ops_exposure.py`, `pipeline/tests/test_ops_exposure.py` | 25 tests passing |
| 1c | Generated tables, CSVs, figures | `ops/generated/` | regenerable; commands in `ops/README.md` |
| 1d | Index marking rev. 2 superseded | `ops/README.md` | complete |
| 2 | Eruption block | `ops/eruption_block/` (`README.md`, three `.sch`, `tcrb_eruption_exposures.csv`, `telegram_skeleton.md`, `generated/`) | DRAFT; syntax not validated by pyscope |
| 3 | Authorship policy | `AUTHORSHIP.md` | one-page DRAFT for the consortium |
| 4 | Email drafts | `committee/work/ops/emails/01_winer_rev3_cover_note.md`, `02_cannon_filters_dossier_authorship.md`, `03_polar_reader_request.md` | unsent; reader name blank |
| 5 | This report | `committee/work/ops/REPORT.md` | **not written — the harness blocked the write** |

Only my own two test files were run.

## 2. Evidence

### 2.1 Visibility (OA §2, TE.F7)

`ops_visibility.py` scans every night from 2026-10-01 to 2027-04-30 at 2-minute sampling, for 27 targets and three Sun limits (−12°, −15°, −18°). Twenty-two targets take their position from the manifest median of their frames; five new standards use catalogue coordinates.

The altitude calculation is pinned to astropy `AltAz` at under 0.02° for a star and for the Sun. Eighteen cells of OA's independent table are reproduced to the 0.1 h quoted.

From `ops/generated/visibility_dates.md` (airmass < 2, Sun < −12°):

| Target | ≥ 0.25 h | ≥ 1.5 h | ≥ 2 h |
|---|---|---|---|
| T CrB | to 2026-10-23, then from 2026-12-16 | from 2027-01-02 | from 2027-01-09 |
| ST LMi | from 2026-10-14 | from 2026-10-30 | from 2026-11-06 |
| M101 | from 2026-11-11 | from 2026-11-26 | from 2026-12-03 |
| NGC 5548 | from 2026-11-25 | from 2026-12-10 | from 2026-12-17 |

- Grism-usable T CrB nights (Sun < −15°, ≥ 0.25 h) end 2026-10-18 and resume 2026-12-19.
- With the airmass limit relaxed to 3 for an eruption, T CrB is unreachable from 2026-11-09 to 2026-12-04 (`ops/eruption_block/generated/`).
- Figure: `ops/generated/fig_visibility.png`.

### 2.2 T CrB exposure for a 16-bit camera (DE §2, OA §2)

T CrB has no QHY frames, so the exposure was derived from pixels. The measurements are in `ops/generated/exposure_measurements.csv` (198 frames, one row each).

**Camera bridge.** Same star, same grism, same exposure on the ASI (era 76) and the QHY (era 78 under MaxIm; era 81 under pyscope):

| Grism | QHY (MaxIm) / ASI | Range over stars | Stars |
|---|---:|---|---:|
| hrg | 0.81 | 0.55–0.87 | 5 |
| lrg | 0.80 | 0.58–0.89 | 3 |

The pyscope-native days give 0.37 (θ CrB, 5 frames) and 0.80 (κ Dra, 8 frames), which is inconclusive.

**T CrB on the ASI at 240 s** (gate-accepted `g_extractions` rows, local sky subtracted), and the prediction for the QHY:

| Grism | Frames used | Continuum peak (ADU) | Brightest on-trace pixel, median / max (ADU) | Predicted brightest raw native pixel on QHY (ADU) | Fraction of 60,200 ADU veto |
|---|---:|---:|---:|---:|---:|
| hrg | 22 of 22 | 389 ± 106 | 1,352 / 2,315 | 2,453 | 4.1% |
| lrg | 29 of 32 | 1,489 ± 923 | 1,988 / 4,572 | 4,812 | 8.0% |

240 s is retained; saturation is not a constraint at quiescence.

**Imaging.** From the as-found zero points and image sizes (era 82: ZMAG g 23.01, r 22.22, i 21.13; FWHM 2.1–2.2″):
- In focus, the exposure that reaches 25 kADU is 4–9 s; at a 4″ defocus it is 14–23 s.
- Scintillation at airmass 1.8 is 11.7 mmag at 1 s and 3.7 mmag at 10 s.
- The request is 5 × 10 s in g, r, i defocused to about 4″, with an in-focus fallback (5 × 5 s g, 5 × 3 s r and i) and a first-night acceptance window of 15–40 kADU.

Figure: `ops/generated/fig_exposure.png`.

Inputs that are assumptions, stated in the script and the request: T CrB V = 10.0 and g, r, i = 10.7, 9.5, 8.3 (a ±0.5 mag bracket is printed); native-pixel factor 1.16 (DE.F9).

### 2.3 Calibration census (DE cross-cutting 4)

From `calibration_census` in `ops_exposure.py`:
- 29,148 QHY raw frames (2026-03-21 to 2026-07-02), with zero raw bias, dark or flat frames among them and zero raw flats in any of the nine filters.
- Twelve exposure times cover 96% of light frames and define the dark list: 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 240, 300 s.
- CCD-TEMP: 29,076 frames at −20 °C, 71 at 0 °C (all on 2026-06-28), 1 without a card.

### 2.4 Header facts used in the request

These were read directly from frames and are cited in §1 and §5 of rev. 3.

- **Detector settings:** `GAIN 56`; `OFFSET 76` on 2026-03-21 and `10` afterwards; `SET-TEMP −20`.
- **Under pyscope 0.3.1 (from 2026-06-28):** `READOUTM ''`, no `IMAGETYP`, `EGAIN 56`, `E-ADU 1.0`, `MAXADU 255`, `FULLWELL 65535`, whole-second `DATE-OBS`, `FOCOFFCG 0` for every filter, `FWALLOFF` all zero, `FOCTEMP` blank.
- **`FWALLNAM`** is written with `CONTINUE` cards that carry no string value; astropy raises `VerifyError` on it.
- **Wheel:** `FWPOS` 0 g, 1 lrg, 2 r, 3 i, 4 ha, 5 hrg, 6 lum, 7 sii, 8 oiii. There were seven named slots on the first pyscope frame of 2026-06-28 and nine afterwards. There is no B filter.
- **File-name time versus `DATE-OBS`** (era 81, 111 frames sampled): median 890 s, 90th percentile 1,207 s, maximum 1,461 s. The file name carries the scheduled start.
- **Reduced headers:** `BINCOMB = 'AVG'`, `OSCNMEAN` 172 ADU, `CALHOST = 'TCC21'`, `master_bias.fts` (40 frames), `master_dark_4s.fts` (40 frames), `master_flat_g.fts`. Masters exist on the control computer and were never synced.

## 3. Findings

| Finding | Status | Why |
|---|---|---|
| U10 — rev. 3 for the camera mounted and the sky available | CLOSED (document) | Written. Sending it is James's. |
| OA §2 visibility table | CLOSED | Script-generated, reproduces OA to 0.1 h, tested; adds the −15° grism limit, windows and first/last dates. |
| OA cross-cutting 1 | CLOSED | Every item is in rev. 3: as-found mode and fixed set-point (§1, §4); visibility (§2); BD+33°2642, HR 1544, Vega and a compact PN through both grisms (§6.B); bias × 50 and darks at the exposures and temperature used (§4); flats in every filter (§3a, §4); dome-lamp ladder (§4); SDSS field (§6.B); READOUTM/IMAGETYP and stale RA/Dec (§5). |
| OA §2 CHANGE items for ops (T CrB dates and Sun limit, nightly block, ST LMi, Be-star item, dither test, M101 templates) | CLOSED | Rev. 3 §6.A–E. |
| DE cross-cutting 4 | CLOSED | Rev. 3 §4 and §5. |
| DE §2, "re-derive exposures for the QHY at 65,535" | CLOSED | §2.2 above. |
| TE.F7 | CLOSED | Rev. 3 §1 ("new instrument"), §2, §9. |
| TE cross-cutting 2 (re-opening protocol, order a–f) | CLOSED | Rev. 3 §3, same order, with the new-epoch rule. |
| TE cross-cutting 3 (site log, calibration trees) | CLOSED (asked) | Rev. 3 §8 q. 2 and 4. The answers are the site's. |
| TE / DE camera-survival and bench-dark question | CLOSED (asked) | Rev. 3 §8 q. 1. |
| OA.E9 (rev. 2 says the 240 s Mode0 dark is "have 0") | CLOSED | Withdrawn in rev. 3 §9; `calib_gaps` says `master_only`. |
| TE TCRB-P0-eruption-block | PARTIAL | Schedule files and README exist. Acceptance also needs a dry-run log from the site, and the `.sch` syntax is unvalidated. |
| ED TCRB-N2-eruption-contingency | PARTIAL | A text skeleton with placeholders exists. The acceptance criterion ("compiles with script-emitted numbers") needs a LaTeX skeleton under `manuscripts/`, which this package does not own. |
| ED.E4 | PARTIAL | `AUTHORSHIP.md` exists as a draft. Consortium agreement, real ORCIDs, a Zenodo DOI and the AAS statements are James's or the CV package's. |

Nothing was rebutted. Nothing is blocked on my side.

## 4. Things found that the committee did not have

1. **There is no B filter in the wheel.** The only clean T CrB anchors are in B (OA.E3), so new photometry cannot continue them. Rev. 3 uses g, r, i and asks the site about it. The 2027 flickering runs are specified in g for the same reason.
2. **T CrB at 240 s is faint on the detector.** The hrg continuum peaks at about 390 ADU per binned pixel. There is no saturation argument for shortening the exposure.
3. **Masters exist on the control computer** for the QHY (bias, 4 s dark, g flat), although the archive holds no raw QHY calibration frames.
4. **Filter focus offsets are zero under pyscope.** The θ CrB hrg trace peak fell to 0.37 of its ASI value on 2026-06-29/30, against 0.82 under MaxIm. This could be defocus or cloud; it rests on one star and five frames, so it is flagged, not claimed.
5. **θ Crt bridges at 0.55–0.58**, well below the other standards (0.80–0.89), in both grisms. I do not have an explanation. The adopted ratio is a median over stars and the range is printed.
6. **The i-band ZMAG is not cleanly per-second.** Its slope against 2.5 log t is +0.22 (g −0.06, r +0.08). The imaging exposures are therefore explicitly a first-night check.
7. **`OFFSET` was 76 on the first QHY night and 10 afterwards.** Separately, the η Hya hrg frame of 2026-03-23 peaks at 28 kADU against about 6.5 kADU later. I did not establish which nights carry which offset, so the early QHY nights may not be the same configuration as the rest of era 78.
8. **File names lead `DATE-OBS` by a median of 15 minutes** under pyscope in June–July 2026. S3 does not audit this.

## 5. Deviations from the brief, and limits

- **`REPORT.md` was not written.** The harness refused the write for a report file. This text is the report; the chair needs to save it to `committee/work/ops/REPORT.md`.
- **`ops/2026-08_observatory_request.md` was not edited.** The brief says both "new files only" and "mark it superseded". The file is a declared output of provenance stage `OPS`, so editing it changes a digest. It is marked superseded in the new `ops/README.md` and in the header of rev. 3. If the chair wants a banner in the old file, it is a one-line edit.
- **A second script and two test files were added** beyond the one script the brief names. Without `ops_exposure.py` the exposure derivation and the calibration list would have been typed. All four files are new.
- **The `.sch` files are pyscope-style, not pyscope-validated.** pyscope is not installed locally and its source could not be fetched (the URLs I tried returned 404). Keyword spelling follows the `BLK*` fields pyscope writes into its headers and my recollection of the format. Each file says so in its first line, and the CSV is the authoritative content.
- **Prose dates in rev. 3 paraphrase the generated tables** ("late October", "mid-December", "second week of January"). If the tables are regenerated with different inputs, the prose must be re-read.
- **Header facts in rev. 3 §1 and §5 come from ad hoc header reads** of named frames, not from a manifest column. F-2 will make them queryable.
- **Counts are as read on 2026-10-03.** Other packages were rebuilding the manifest during this work.

## 6. Files changed

All are new; nothing existing was modified.

```
AUTHORSHIP.md
ops/2026-10_observatory_request_rev3.md
ops/README.md
ops/generated/  visibility.csv, visibility_table.md, visibility_windows.md, visibility_dates.md,
                fig_visibility.{png,pdf}, exposure_measurements.csv, exposure_bridge.md,
                exposure_grism.md, exposure_eruption.md, exposure_imaging.md,
                calibration_census.md, eruption_plan.md, fig_exposure.{png,pdf}
ops/eruption_block/  README.md, tcrb_eruption_1_peak.sch, tcrb_eruption_2_early_decline.sch,
                     tcrb_eruption_3_late_decline.sch, tcrb_eruption_exposures.csv,
                     telegram_skeleton.md, generated/*
pipeline/scripts/ops_visibility.py
pipeline/scripts/ops_exposure.py
pipeline/tests/test_ops_visibility.py
pipeline/tests/test_ops_exposure.py
committee/work/ops/emails/  01_winer_rev3_cover_note.md, 02_cannon_filters_dossier_authorship.md,
                            03_polar_reader_request.md
```

## 7. Ledger and provenance changes requested of the chair

1. **Provenance stage `OPS`** (`provenance.py:2104`): point `writes` at `ops/2026-10_observatory_request_rev3.md`. Set `build_cmd` to the commands in `ops/README.md`. Add to `reads`: `table:frames`, `table:g_extractions`, `table:detector_params`, `table:calib_gaps`, and the two `ops_*.py` scripts.
2. **TCRB-P0-restart:** reword to "evening block to about 2026-10-23; morning block from about 2026-12-16; second instrument (QHY), no splice to 2025 without at least three nights of same-night standards". Blocker: rev. 3 not yet sent.
3. **TCRB-P0-eruption-block (new):** in progress; evidence `ops/eruption_block/`; remaining is the site's dry-run log.
4. **TCRB-N2-eruption-contingency (new):** in progress; text skeleton exists; the LaTeX skeleton with emitted macros should be owned by the T CrB package.
5. **TCRB-C2-2026-runs:** move to the 2027 backlog with the computed date (first ≥ 2 h night 2027-01-09) and filter g, not B.
6. **T CrB Phase B:** record that B cannot be re-observed with the current wheel.
7. **CV ST LMi season request:** queue item, from 2026-10-30 (≥ 1.5 h) or 2026-11-06 (≥ 2 h); not on the paper's critical path.
8. **Be-star season 2 and the dither test:** 2027 backlog; included in rev. 3 §6.C as SYNTHESIS requires.
9. **CV-R15 / ED.E4:** add `AUTHORSHIP.md` as evidence; stays open pending James.
10. **For the F-2 / F-3 owners:** `OFFSET` 76 → 10 inside era 78 and the 7-slot → 9-slot wheel change on 2026-06-28 are configuration boundaries.
11. **For S3:** add the file-name-versus-`DATE-OBS` lag to the DATE-OBS audit.

## 8. Needs James

1. **Read and send rev. 3** to Winer with the cover note, and change the DRAFT status line. The dome may already be open, and the as-found flats (rev. 3 §3 step a) lose value every night anything is adjusted.
2. **Confirm the T CrB quiescent magnitudes** (V 10.0, g 10.7, r 9.5, i 8.3 in `ASSUMED` in `ops_exposure.py`), or supply better ones and re-run the report.
3. **Decide defocus versus in-focus** for the T CrB imaging if the site cannot offset focus per block.
4. **Eruption block:** supply the observer email and observing code, who is called on an alert, and the pre-agreed eruption author list (`AUTHORSHIP.md` §6.4).
5. **Send the Cannon email**; note the 2026-11-15 date for the filter curve.
6. **Choose the outside polar reader** and send that request only after CV-R1, R6 and R12 land. The draft quotes a colour amplitude that is not yet a script-emitted number.
7. **Take `AUTHORSHIP.md` to the consortium.** Its §10 lists seven questions it needs answered; ORCIDs and a Zenodo DOI follow.
8. **Ask for a B filter** if the T CrB photometric anchors matter (rev. 3 §8 q. 5).