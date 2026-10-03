# Package `cv-stats` — the CV paper's statistics, redone to the committee's standard
Owns: `pipeline/macro_phot/{phase3,final_science,numbers_cv,figures_cv}.py` and sibling CV analysis modules, `pipeline/scripts/run_cv_{phase3,final,paper}.py`, their tests, `products/phot/cv_timeseries.sqlite`. You do NOT edit `manuscripts/CV_TimeSeries/main.tex` — you deliver numbers, figures, and proposed replacement text in your report; a later package integrates the manuscript. Findings: CV-R1…R8 (DS.F1–F3; RF.B1–B4, M1–M3, M5–M7, minors 1, 4, 5; PH.P1–P6; ED.E2 colour result).
Do:
1. R1: paired/permutation band-offset test on same-cycle and same-night edges, combined across eras, scatter-based errors, trials factor.
2. R2: injection grid per band with per-band ramp widths AND on a 2024 High Gain night; carry the signed matched-cell bias per band/era; remove `max(chi2nu,1)`; replace the constant "Monte-Carlo sigma" with a real per-edge bootstrap or delete the claim; add grid-quantisation to the budget. Decide D3: is the g−i offset astrophysical or estimator bias — say what the evidence supports and no more.
3. R3: O−C refit with per-band constants, night-level epochs (N=17), an era-offset nuisance term; Ṗ bound quoted each way and with 2024 dropped.
4. R4/R5: YZ Cnc excluded-amplitude statement (which amplitudes, in how many run-filters); spot-longitude stability in degrees, split by accretion state from `p3_state_night`.
5. R6: quantify ST LMi g−r (and other colours) vs phase — amplitude, phase of extrema, per-era repeatability with errors — and show it survives a ≤120 s pairing window.
6. R7: write the edge model down; publish the per-edge χ²ν distribution; best/representative/worst fit figure.
7. R8: frame-count language (staged vs measured), S/N floor on "measurement", EU UMa i/r, FWHM > 5″ cut effect.
Every standing rule in SYNTHESIS §5 applies. Update `numbers_cv.py` so every new number is a macro.
