# SN 2023ixf — RLMT +5.4 to +50 d validation and limits release

Generated 2026-10-04T23:19:36Z by `pipeline/scripts/run_sn_photometry.py release` (SN-PHOT v1.0 (2026-10-04), git ceea0fadd4) from `products/sn/sn2023ixf.sqlite` and the manifest's Gate 0 tables. Every number below is a query.

Campaign: 1129 unique frames (Gate 0 freeze, global (basename, jd) dedup, alias-merged) on 35 nights, night labels 2023-05-19 to 2023-07-06 (local-noon split).

Usable photometry per filter code (frames / nights):
- G: 112 / 32
- R: 54 / 24
- I: 55 / 19
- H: 8 / 4
- O: 21 / 17
- 1: 15 / 8

Files:

- `sn2023ixf_rlmt_frames.csv` — 1040 rows
- `sn2023ixf_rlmt_nightly.csv` — 104 rows
- `saturation_matrix.csv` — 207 rows
- `slot6_frames_not_promoted.csv` — 83 rows
- `calibration_sets.csv` — 14 rows
- `filter_crosswalk.csv` — 42 rows
- `template_table.csv` — 20 rows
- `residuals_vs_li2025.csv` — 58 rows
- `variability_limits.csv` — 21 rows
- `late_time.csv` — 6 rows

Magnitudes: `mag_natural` is the RLMT natural system zeroed to PS1 (REFCAT2) at (g-i)=0.8; `mag_ps1` applies the campaign colour term with the SN's own (g-i). Errors are photon, scintillation (Young 1967 law, scale fitted), floor and zero point in quadrature. Frames flagged `usable=0` carry the rule that excluded them in `exclusion`.

The slot-'6' series was triaged and NOT PROMOTED (Gate 0); its frames are listed as they are.

Pipeline: this repository (`pipeline/macro_sn/`, `pipeline/scripts/run_sn_gate0.py`, `pipeline/scripts/run_sn_photometry.py`). Archive DOI: pending (Zenodo deposit is an owner action).
