# Seat 5 — Telescope Systems Engineer · plan review 2026-10-03

Scope: read-only review of what is on disk. Evidence is from `products/manifest/rlmt-manifest.sqlite`
(tables `frames`, `eras`, `s1_batch`, `frame_dispersion`, `g_extractions`), raw FITS headers of ~90 frames,
one 56-frame sky-stack test, and the code/ledger lines cited. Queries are reproducible from the table and
column names given. My question throughout: *what physically changed on the telescope that night?*

## 1. Existing work — findings

**F1 · MAJOR · The era registry is a header history, not a hardware history.**
`manifest.py:219-245` keys eras on (READOUTM, NAXIS, XBINNING, EGAIN); `eras.csv` has 84 rows. The
telescope has had four cameras and at least eight mechanical/orientation states, and the registry sees
neither count. From the S1 solutions (`s1_batch.rotation_deg`, 49,380 solved) and headers:

| State | Nights | Evidence |
|---|---|---|
| AC4040 (GSENSE) | 2023-02 → 2024-03-26 | rot 179.0–179.4°, 0.539″/px |
| Andor iKon, state 1 | 2024-04 → 10-31 | rot 180.8° (−2.0° offset early season), 0.809″/px |
| Andor, re-seated | 11-01 → 11-07 | rot 180.0° |
| Andor, re-seated | 11-08 → 11-20 | rot 180.7°; wheel renamed (lower-case, hrg/lrg appear 11-08) |
| Andor, rotated | 11-30 → 12-10 | rot 204.3° (2,261 frames; includes 477 staged VV Pup frames) |
| ASI IMX455, pre-monsoon | 2024-12-14 → 2025-06-23 | rot +0.43°, `FLIPSTAT='Flip/Mirror'`, MaxIm 6.40 |
| ASI, post-monsoon | 2025-10-15 → 2026-03-16 | rot 179.75°, `FLIPSTAT` blank, MaxIm **6.30** |
| QHY600, night 1 | 2026-03-21 | rot 3.3° |
| QHY600 | 2026-03-22 → | rot 179.74°; pyscope-native headers from 06-28 |

Era 76 (70,603 frames) spans the 2025 monsoon boundary as one era. It is not one configuration: the raw
median-sky illumination pattern of the two halves correlates +0.78 directly and **+0.98 after a 180°
rotation** (14 unrelated g fields each side). The rotation differs from an exact software flip by ~0.7°
(+0.43° → −0.25°), so the camera was also re-seated. Every flat, dust-donut map, hot-pixel mask, amp-glow
model, grism trace prior and zero-order prediction is orientation-specific; the archive holds **no era-76
flat after 2025-01** (S0b: 70 flats, all 2024-12/2025-01).
*Close:* add a `mech_epoch` table (camera serial/driver, rotation step > 0.3°, FLIPSTAT, SWCREATE, wheel
map) beneath `eras`; S0b re-counts calibration coverage per mech_epoch; a test asserts no master is applied
across a boundary. *Would change my mind:* a hot-pixel map showing the same pixel coordinates before and
after 2025 monsoon would make the flip physical-only and leave darks valid (flats still not).

**F2 · MINOR · The grism identity gate's parity rationale is wrong.** `macro_grism/gate.py:15-19` says the
mount is a German equatorial that flips, and "the pier side is not in the headers". `TELPIER` is in every
pyscope-updated header I opened (2024-05 onward) and reads `pierEast` at hour angles from −7.8 h to +9.8 h;
of 398 solved nights, **zero** contain both parities. Orientation is a function of date (F1), not pier
side. The gate survives because it tries both parities, but that is a free binary parameter it does not
need. *Close:* fix parity per mech_epoch, re-run the 42 validation frames, verdicts unchanged.

**F3 · MAJOR · The per-frame grism dispersion is not a physical measurement.** `g_extractions` (ACCEPT rows,
flanking method): hrg `disp_a_per_px` = +0.46, −0.67, −1.78, −1.22, +1.55 … (range −1.80 to +1.59); lrg
−1.26 to +2.00; many rows have `x_o2a == x_o2b`. `wavelength.py:85-88` allows |disp| in (0.3, 3.0) of either
sign, while its own comment records hrg ≈ 1.59, lrg ≈ 1.9–2.2 Å/px. A grism at a fixed distance from the
detector has one dispersion and one sign per mech_epoch; a 650 µm focus move changes it by <1%. The O₂
pair finder is locking onto different features frame to frame. TCRB-A3 ("per-frame self-anchored
wavelength solution") inherits this design.
*Close:* solve dispersion and sign once per (grism, mech_epoch) from stacked high-S/N calibrator frames;
per-frame fit is zero-point only (Hα + O₂ B as a consistency check). Acceptance: Hα–O₂B pixel separation
constant to <1% across all 247 T CrB frames.

**F4 · MAJOR · Grism focus was not held constant on T CrB.** hrg normally runs +650 focuser counts from
lrg/g (the wheel's filter-offset table). On 10 of 51 T CrB nights it did not: offset 0 on 2025-03-09, -10,
-12, 04-14, 04-16; +1041 on 03-02; ~+570 or +428 on four others (14 hrg frames). Those nights have a
different line-spread function. Nothing in the T CrB or BeStar plan regresses line width on focus.
*Close:* regress `halpha_width_px` and cross-dispersion trace FWHM on (focus − nightly g focus); either EW
is shown invariant (then say so) or the off-nominal nights are flagged and "profile / wing velocity"
claims exclude them.

**F5 · MAJOR · CV manuscript: the server flat-field leaves a ~5% ramp in the densest epoch.** Median of 14
unrelated server-reduced g frames, pedestal removed, 32× binned: pre-monsoon-2025 epoch shows a monotonic
left→right sky ramp **+1.8% → −3.0%** (corr with its own L–R mirror −0.91, rms 1.6%); post-monsoon epoch
~2% (one edge). Reduced headers name the same file, `master_flat_g_read0_g100_o30_2x2.fts`, on both sides
of the flip; the masters were never synced, so which sky flats made them is unauditable. The recipe also
changed inside "era 76" (2025-03: "mode: ccd", dark `master_dark_-5s`; 2025-11: "mode: cmos", "No bias
subtraction applied"). The paper states the calibration goal was not met (`main.tex:1519`) but never
quantifies flat-field error. Differential light curves at fixed pointing are largely immune; the catalogue
tie, colour terms and check-star bias are not.
*Close:* ≥100-frame sky stack per (mech_epoch, filter); catalogue-tie residual vs detector x; one
paragraph and one number in §phot. *Caveat:* 14 frames per epoch; a real sky gradient is unlikely to
survive a 14-field median but is not excluded. *Would change my mind:* tie residuals flat in x to <0.5%.

**F6 · MAJOR · QHY-era reduced twins are counted as canonical.** `…_calibrated.fts.fz` in `reduced/` is not
matched to its raw parent: 25,442 era-79 and 1,681 era-82 reduced rows have `is_canonical=1` beside the
raw rows (1,503 JD-identical canonical pairs in ten May-2026 nights alone). Eras 78/79 and 81/82 are one
camera with and without the overscan border (4800×3211 raw, 4787×3193 trimmed). Convention 3 of the
roadmap is violated for every 2026 count. *Close:* strip `_calibrated` in the link rule; canonical frames
2026-03-21→07-02 drop by ~27k; eras 79/82 become reduced-only aliases.

**F7 · MAJOR · The ops request asks for something the telescope can no longer do, in a window that is
closing.** `ops/2026-08_observatory_request.md` Item A.1: "continues the 2025 Hα EW series homogeneously —
same mode, same exposure". The 2025 series is ASI/Mode0/MaxIm, unflipped-at-+0.4°. October's instrument is
QHY600/pyscope (GAIN 56, blank READOUTM), at 180°, with the wheel reloaded (hrg at `FWPOS` 5 on
2026-06-30; slot 4 in 2025's `FWALLNAM`) and lrg re-seated (trace PA −2° in era 76 → +7° in the QHY eras,
`frame_dispersion`). It is a new instrument that needs its own trace, wavelength and response solution.
Visibility (astropy, Winer): T CrB at airmass < 2 in nautical dark for **70 min on Oct 5, 42 on Oct 15,
14 on Oct 25, 0 from Nov 5 to mid-December**; a ≥2 h run is first possible in mid-January 2027 (morning).
The nightly block (≈12 min) fits for about three weeks, always at airmass 1.7–2.9, where slitless
differential refraction is at its worst. The OPS stage is also STALE (`check_pipeline_status.py`).

**F8 · MINOR · The manifest's focus and image-quality columns are not what they claim.** `focuspos` comes
from MaxIm's `FOCUSPOS`, which is stuck (19058, 9104) across 2025-10/11; the true value is pyscope's
`FOCPOS`. `fwhm` exists for ~530 rawimage frames after 2024-10; focus temperature for almost none. The S1 autopsy
label "defocused" (1,569 frames, 59 nights; `astrom.py:779-835`) is right in kind but the cause is
visible in the headers: 665 of 728 frames on 2026-06-25 sat at focuser 9263 all night, ~1,000 counts from
that month's norm, with median elongation 2.0. Elongated out-of-focus images mean astigmatism. On
2023-05-04 the SN pre-explosion G/R/X templates (FWHM 14–16 px) were ~300 counts from the focus that gave
5 px in H/O the same night at 10 °C.

**F9 · NOTE · Focus offsets identify wheel contents.** Slot '6' carries a +1400-count offset in every
month from 2023-02 to 2024-03: one thick element sat there the whole AC4040 era. With 91 of 143 NGC 5548
slot-6 frames measured dispersed, roadmap conflict C2 is settled mechanically: it was the grism in March
2023 too. Slot 'W' jumps from ~+200 to ~+2300 in 2024-02: its element changed.

**F10 · MAJOR · "Legacy_Rigel" is mis-premised.** `project_plan.py:1542` says the legacy archive is "a
different telescope (Rigel System, FLI ProLine PL16803)". Headers say that holds through 2015 day ≈120
(≈7.8k of 210k files). From late 2015 it is the 0.508 m, f = 3454 mm telescope ("Gemini", later "Iowa
Robotic Telescope") with six cameras: Apogee F47, SBIG 6303e, Andor Aspen CG42, Andor iKon-L 936, SBIG
STXL-6303, and in 2022 the SBIG Aluma AC4040 — the 2023 RLMT camera.

**F11 · NOTE · Stale blockers.** The ledger blocks TCRB-P0-bitdepth, -ladders, -shutter-timing, -B1 and
DW-P04 on "S2's tables were destroyed". `s2_*` tables are populated (built 2026-08-20, S2 v1.2).

## 2. Plan hardening, per project

**TCrB_Monitoring**
- ADD `TCRB-P0-mech-epoch`: assign every staged frame a mech_epoch. *Accept:* table emitted; 2025 series
  confirmed single-epoch (it is: 2025-02-20→06-23).
- CHANGE A3 to fixed dispersion per epoch (F3). CHANGE A7 to include the focus-offset regressor (F4).
- CHANGE P0-restart: the October frames are a *second instrument*. *Accept:* same-night θ CrB + Vega in
  hrg/lrg on ≥3 nights before any 2026 EW is spliced to 2025; splice offset reported with error.
- CHANGE C2-2026-runs: honest scope is zero runs before mid-January 2027; ≥6 runs is a spring-2027 goal.
  The flickering subsection should not be in a Dec 2026–Jan 2027 submission.
- CHANGE P0-calib-acquisition: Mode0 240 s darks need the ASI camera, which is no longer the installed
  camera. Ask whether it still exists and can be cooled on a bench (darks need no telescope). If not, the
  task is impossible and the flanking-band background is the method, stated as such.
- ADD `TCRB-P0-eruption-block`: no artefact on disk shows the pyscope block was loaded by 2026-09-01.
  *Accept:* the schedule file committed under `ops/` and a dry-run log.
- Unblock the four S2-gated tasks (F11).

**CV_TimeSeries** (result review)
- ADD the F5 closing test and paragraph. ADD mech_epoch to `cv_frames` and to the "camera" language:
  VV Pup has four mechanical states (Andor 180°, Andor 204°, ASI pre, ASI post), not two cameras; EU UMa
  three. *Accept:* per-state zero-point offsets tabulated; folds never combine states without one.
- Not satisfied until F5 is answered; everything else in my domain is acceptable.

**SN2023ixf_LightCurve**
- CHANGE the template plan: the 2023-05-04 G/R/X templates are defocused by ~300 counts (F8); late
  templates are Andor (0.81″/px, 2024) and QHY (0.45″/px, 2026, the 03-21/22 commissioning nights,
  one at 3.3° and one at 180°). *Accept:* template table lists camera, mech_epoch, FWHM, focus offset;
  narrowband templates only from H/O on 05-04.
- G0c grism triage: `sn_g0_triage_summary` shows 0 nights with Hα at S/N > 10 and no wavelength source.
  The slot-6 element is identified (F9) but uncalibrated; honest scope is an appendix. DROP promotion.
- Re-run Gate 0 after S0 goes fresh; single orientation (179.0–179.4°) all season is good news.

**BeStar_Grism**
- CHANGE §3.4 "three instruments": it is at least five states (Andor; ASI pre-monsoon; ASI post-monsoon
  flipped; QHY night 1; QHY) with the grisms re-seated at the QHY swap. *Accept:* free offsets per
  mech_epoch in S9; standards bridge each boundary or the boundary is declared unbridged.
- CHANGE S0-header-rescrape to pull `FOCPOS`, `FLIPSTAT`, `TELPIER`, `FWPOS`, `FWALLNAM`, `CCD-TEMP`.
- ADD the F4 focus regressor to S5 (delivered resolution) and S11 (systematics).
- S2-calibration: era-B (ASI) frames cannot be acquired on the QHY; same bench question as T CrB.

**DwarfGalaxy_AGN_Survey**
- DROP slot-6 broadband photometry of NGC 5548 (P41/P42 as written): the element is a grism (F9). Honest
  scope: R/L/G only (2023-02→06), or a clearly labelled zero-order series.
- CHANGE P11-flats: no AC4040 flat can be taken now; the camera left in March 2024. The task is a
  decision tree only. P02 dossier gains the focus-offset evidence (slot 6 constant; W changed 2024-02).
- Re-scope before execution.

**Legacy_Rigel**
- CHANGE the premise and the name (F10). CHANGE RIG-L0-multi-archive to key on camera + focal length, and
  allow the 2022 AC4040 frames to share detector characterisation with RLMT eras 1–2 (linearity, ceiling).
- ADD a legacy mech_epoch census (seven cameras). *Accept:* camera timeline table with first/last night.

## 3. Cross-cutting — top five

1. **Mechanical-epoch layer under the era registry** (F1), with calibration validity bounded by it.
2. **Re-opening protocol, in this order:** (a) before anyone touches the instrument, as-found flats in
   every filter including hrg/lrg, plus bias/darks; (b) through-focus sequence on a bright field to
   measure astigmatism and decide on collimation (F8); (c) focus-vs-temperature run and filter offsets;
   (d) Vega/θ CrB grism frames for trace PA and zero-order position; (e) 10-cycle wheel repeatability;
   (f) only then science. Any re-collimation or camera removal starts a new mech_epoch and repeats (a).
3. **Ask the site to log** every camera/wheel/focuser removal, collimation, MaxIm/pyscope version change,
   flip setting and wheel map, dated. Ask for the two server `calibrations/` trees (request Item B.0) —
   F5 cannot be audited without them.
4. **Fix QHY-era dedup** (F6) before any stage is re-run; all 2026 counts move.
5. **Header-derived image-quality series**: `FOCPOS`, focus residual from a per-epoch temperature model,
   S1's source size and elongation, per frame. Replaces morphological guessing with a prior and gives
   every project the focus regressor.

## 4. Verdict

| Project | Verdict |
|---|---|
| TCrB_Monitoring | Execute with amendments (F3, F4, F7 binding; flickering deferred to 2027) |
| CV_TimeSeries | Execute with amendments — **not yet satisfied**: F5 open, mech_epoch wording |
| SN2023ixf_LightCurve | Execute with amendments (template table; grism to appendix) |
| BeStar_Grism | Execute with amendments (five states, not three) |
| DwarfGalaxy_AGN_Survey | Re-scope (NGC 5548 slot 6; flats unobtainable) |
| Legacy_Rigel | Re-scope: continue ingest, correct the premise before any strategy is written |

Satisfied with: the S1 astrometry itself (it is what exposed the hardware history), the fpack geometry
repair, and the measured-dispersion classification of slot 6.
