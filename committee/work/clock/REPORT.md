The stamp-convention check, the full REPORT.md text and the summary follow. All numbers come from `products/clock/clock_transits.sqlite`, built by `pipeline/scripts/build_s3b_clock_transits.py`.

---

## 1 · Chair's correctness check: start or mid-exposure stamps (≤250 words)

The transits can't decide this: half-exposures on the AC4040 are 2–32 s, and the event errors are 46–1,200 s.

I used a header test instead, which needs no ephemeris. For consecutive frames with exposures eₐ and e_b, I compare the shortest stamp-to-stamp gap going a→b with the gap going b→a. All overheads cancel in that difference. It leaves eₐ−e_b if DATE-OBS marks the start of the exposure, and 0 if it marks the middle. The ratio is therefore 1 for a start stamp and 0 for a mid stamp; cells where the exposures differ by less than 8 s are dropped as timing jitter. Unit tests confirm it gets both cases right. The result is in table `s3b_stamp_convention`.

- **AC4040 plain High Gain, MaxIm 6.40 (era 7): START confirmed.** Median ratio 1.02 (IQR 0.94–1.33) over 50 cells.
- **AC4040 StackPro (eras 1 and 6): MID confirmed.** Era 1 gives 0.01 (−0.04 to 0.09, 22 cells); era 6 under 6.40 gives −0.06 (5 cells).
- **Era 2 (High Gain, MaxIm 6.30): undetermined** (0.59, IQR 0.26–0.90).
- **iKon, ASI and QHY: START**, ratios 0.98–1.03.

**Consequence for S3:** S3 adds EXPTIME/2 to every stamp. For StackPro frames that puts the time late by EXPTIME/2, up to 512 s. That affects 20,364 frames; S3 owns the fix.

**CV impact: none.** All ST LMi era-7 frames are plain High Gain under MaxIm 6.40, so they are start-stamped. Even if they were mid-stamped:
- Every era-7 O−C epoch night used 64 s in G, R and I, so the shift would be a uniform −32 s.
- The g−i band offset would change by exactly 0 s.
- The combined O−C would move by about −7 s (era 7 carries 21% of the weight), against ±97 s.
- The quadratic (Ṗ) term would move by about −6×10⁻⁷ s/cycle², against σ = 2.8×10⁻⁶ (0.2σ).

**Verdict:** start confirmed for era 7; mid confirmed for StackPro; undetermined for era 2.

---

## 2 · REPORT.md (for the chair to file at `committee/work/clock/REPORT.md`)

### Package `clock` — the absolute clock from archived transits (F-8)

**What was done**
- **Standards.** 27 clock standards were set up:
  - 25 ExoClock planets. The catalogue (776 planets) was fetched 2026-10-03 and cached with its sha256. Each planet also has an independent NASA Exoplanet Archive ephemeris, used only as a cross-check.
  - NSVS 07826147: Baştürk et al. 2026, arXiv:2602.09925, cross-checked against Pulley et al. 2025, MNRAS 544, 24.
  - GK Vir: Yates et al. 2026, arXiv:2602.17800, cross-checked against Parsons et al. 2010, MNRAS 407, 2362.
  - V2301 Oph, NN Ser, QS Vir and DE CVn were refused before any photometry, with reasons recorded in `s3b_refused`.
- **Census.** 74 series and 75 predicted events were found. 54 events were admitted. Admission allows a ±2,130 s window around the prediction so that the selection does not assume the clock it is testing.
- **Photometry.**
  - Frames were calibrated with era-matched master darks and flats, and stars matched by astroalign.
  - Each frame was measured in three apertures.
  - Times are the header JD (= DATE-OBS) + EXPTIME/2, recomputed as BJD_TDB at the target's coordinates.
  - 48 series and 9,187 of 9,308 frames were measured.
- **Fitting.**
  - Each event uses an exact limb-darkened occultation model (valid for any size ratio) and a blind search over the whole run.
  - Telescope re-pointing jumps get their own offsets, and the noise model is a 300-draw block bootstrap plus prayer-bead (the larger is adopted).
  - Each fit was repeated under 8–10 analysis variants, and injection-recovery was run at ±1,065, ±300, ±120 and 0 s.
- **Result.** 38 events were timed:
  - 16 timing grade (σ ≤ 120 s; only these enter the era means).
  - 17 low precision.
  - 5 one-sided eclipses timed with the shape fixed from the same star's complete eclipses; reported, never averaged.

**O−C per clock era** (`s3b_era`; O−C > 0 means our stamps are late)

| era | camera | O−C (s) | verdict |
|---|---|---|---|
| A | AC4040, 2023 (StackPro) | +40 ± 115 (WASP-43 b); **+24 if the mid-stamp correction is applied** | central value inside 120 s; 2σ reaches the criterion |
| B | AC4040, 2023–24 | +113 ± 106 (WASP-183 b, StackPro); **+81 mid-corrected** | same |
| C | iKon | −30 ± 47 (TrES-3 b) | same |
| D | ASI, Dec 2024–Jun 2025 | **no standard** | — |
| E | ASI, Oct 2025–Mar 2026 | **+7.3 ± 14.7** (NSVS 07826147, 9 eclipses, +11.4 ± 15 including its own timing-variation systematic; WASP-52 b in 3 bands, −18 ± 37) | PASS |
| F | QHY | **+1.7 ± 16.3** (GK Vir) | PASS |
| G | QHY under pyscope | **no standard** | — |

- **Injection test.** A ±1,065 s offset was recovered in 64/128 and 80/128 draws (the rest move the event out of the run). The signed bias is −8.9 and −5.3 s.
- **StackPro vs plain on the same transit.** These pairs (WASP-183 b: −70 ± 189 s; WASP-43 b: −172 ± 1,208 s) are too noisy to say anything.
- **One misfit.** WASP-43 b on 2024-03-04 fits at +11,941 s. Its error is 616 s, so it is graded low precision and stays out of the means.

**S3's AG LMi residual (RF.M4).** The −294 ± 59 s residual against VSX was eclipsed 1,584 cycles after the VSX epoch.
- Subtracting the measured era-B clock offset leaves −406 ± 121 s.
- That corresponds to a true period of 1.3590146 ± 0.0000009 d. This is 0.09σ from the Gaia DR3 period of 1.3590118 ± 0.0000301 d.
- So the residual comes from the stale VSX period, not from the clock.

**The shared ~1,065 s CV offset**
- EU UMa (era E): the clock is +7 ± 15 s while the offset is 1,060 s, 72σ apart.
- ST LMi, eras B, C and E: the clock is +113 ± 106, −30 ± 47 and +7 ± 15 s.
- **Conclusion: the offset is not a clock error.** It belongs to the catalogue ephemerides.
- ST LMi's era-D nights (11 nights, 2025-02 to 05) have no standard of their own.

**Time-stamp provenance** (`s3b_stamp_audit`)
- Every time used is within 1 ms of DATE-OBS.
- File names are never read for time. The time in the file name precedes DATE-OBS by a median of 866 s in era E and 174 s in era F, consistent with the ops package's 890 s.
- The start-versus-mid test is section 1 above; its consequence is a defect in S3 (see below).

**Stated bounds (DE.F7)** (`s3b_bounds`)
- Rolling-shutter skew: AC4040 ≤ 0.75 s, ASI ≤ 0.32 s, QHY ≤ 0.40 s.
- iKon shutter: ≤ 0.1 s, **assumed** (no vendor figure found).
- The report also states that TELUT and DATE-OBS are not independent clock cards.

**Finding status**
- **RF.M4 — CLOSED.** The residual is printed and explained.
- **DE.F7 (timing sentences) — CLOSED.** The iKon shutter bound is flagged as assumed.
- **OA.E6 — PARTIAL.** The 1,065 s clock hypothesis is excluded at ≥9σ in every era that has both CV data and a standard. The 120 s criterion is certified only for eras E and F; D and G have no standard.
- **DS.F10 — PARTIAL.** Absolute time is validated to ±15 s for Oct 2025 to Jun 2026 and ±47–115 s for 2023–24. Dec 2024 to Jun 2025 is unverified, and the site NTP log is still needed.
- **F-8 — PARTIAL.**

**New defect, for S3's owner**
- `macro_core/timing.py`'s StackPro policy (start + EXPTIME/2) is wrong. StackPro frames are mid-stamped: eras 1 and 6 confirmed, era 3 not testable, 20,364 frames in total. The mid-exposure time should equal DATE-OBS.
- The effect on SN 2023ixf phases is ≤ 512 s, negligible on a scale of days.
- CV and the YZ Cnc O−C are unaffected.
- Era 2 (High Gain, MaxIm 6.30) remains undetermined.

**Citation for cv-literature (main.tex:331–336)**
- The eclipsing binary is AG LMi. Its ephemeris is from VSX: Watson, Henden & Price 2006, already in the bibliography as `watson2006`. The entry is OID 167169, epoch HJD 2458211.194, P = 1.3590176 d, retrieved via VizieR B/vsx.
- The period check uses Gaia DR3 eclipsing binaries: Mowlavi et al. 2023, A&A 674, A16. This entry is not yet in `references.bib`.
- The manuscript text should now say the bound is limited by the stale period, and cite S3b.

**Files changed**
- New: `pipeline/macro_core/clock_transits.py`
- New: `pipeline/scripts/build_s3b_clock_transits.py`
- New: `pipeline/tests/test_clock_transits.py` (58 tests, all pass)
- `products/clock/` (`clock_transits.sqlite`, `ephemerides.json`, logs)
- `docs/pipeline/s3b_clock.html` and 9 figures in `docs/pipeline/figures/s3b/`

**For the chair**
- Register S3b in `provenance.py` and `site.py` (I don't own those files).
- Fix the StackPro policy in S3.
- Ledger: F-8 → partial; CV-R11 evidence → `docs/pipeline/s3b_clock.html`.
- Note for the incident record: my 02:08 photometry run was interrupted and has been fully re-run; the products are complete.

**For James**
- The site NTP or clock log.
- The iKon shutter open/close time.
- At re-opening, one NSVS 07826147 eclipse under pyscope (era G).

---

## 3 · Summary (≤250 words)

The archive's transits and eclipses rule out the observatory clock as the cause of the ~1,065 s CV offset.

- **Method.** 27 literature standards, 54 admitted events and 38 timed (16 timing grade). Each fit is blind to the prediction, with block-bootstrap errors, analysis variants and signed injection-recovery.
- **Clock by era.**
  - Eras E (ASI) and F (QHY) pass at |O−C| + 2σ < 120 s: +7.3 ± 14.7 s and +1.7 ± 16.3 s.
  - Eras A, B and C (AC4040, iKon) have central values of +40, +113 and −30 s, but their 2σ intervals reach the criterion.
  - Eras D and G have no standard.
- **The 1,065 s offset.** In every era that has both CV data and a standard, the clock sits ≥9σ away from it. The offset belongs to the catalogue ephemerides, not the clock.
- **AG LMi −294 s.** After removing the era-B clock offset, this is a stale VSX period: the implied period is within 0.09σ of Gaia DR3.
- **Stamp convention.** A new cadence test shows StackPro frames are mid-stamped, so S3 makes them late by EXPTIME/2 (20,364 frames). That is a fix for S3's owner.
- **CV impact.** None: ST LMi era 7 is plain High Gain and start-stamped. Even if it were mid-stamped, the g−i band offset would change by 0 s, the O−C by about −7 s against ±97 s, and Ṗ by 0.2σ.
- **Time stamps.** All times come from DATE-OBS (within 1 ms). File names lead DATE-OBS by a median of 866 s and are never used for time.
- **Status.** RF.M4 and DE.F7 are CLOSED. OA.E6, DS.F10 and F-8 are PARTIAL: eras D and G have no standard, and the site NTP log is still needed.