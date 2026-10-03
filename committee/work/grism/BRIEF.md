# Package `grism` — one grism calibration library for four projects
Owns: `pipeline/macro_grism/*`, `pipeline/scripts/run_g_*.py`, `pipeline/rlmt_diagnostics/dispersion.py`, `pipeline/scripts/run_s2c_dispersion.py`, their tests, `products/grism/`. Findings: G-1…G-5, F-6, D1 (OA.E1, E2, E4, E7; PH.P7, P8; TE.F2, F3, F4; DE.F1 grism part; RF T CrB points; DS.F8).
Do NOT write the shared manifest in this wave: write to `products/grism/grism.sqlite` with the g_* / frame_dispersion schemas; the integrated rebuild will land them.
Do, in this order:
1. **D1** — settle hrg (and lrg) dispersion on hot stars on disk (Vega 2026-04-15, θ CrB, Spica, any A/B star): identify ≥3 lines per grism per mechanical epoch (Hα, Hβ if on-chip, telluric 6277/6867/7186/7594). Show the spectra with line ids. State plainly whether hrg is ≈0.47 or ≈1.59 Å/px and what the delivered resolution (LSF FWHM in Å and km/s) is. This single result decides whether profile/V-R science is in scope.
2. G-1 fixed dispersion + sign per (grism, mech epoch); per-frame zero point only; Hα–O₂B separation constancy across all T CrB frames. (Mechanical epochs: see TE.F1's table; coordinate through your report, the formal table is being built by `foundation-s0`.)
3. G-2 variance and saturation from measured gain (use DE's memo values ≈1.0 e⁻/ADU for Mode0 as provisional and read from a single config point so the `detector` package's final numbers drop in); replace the 16.3 kADU rail by hot-pixel mask + true clip.
4. G-3 pixel-based identity gate (spectral fingerprint); re-adjudicate every `header_off_target` frame; measure false-accept rate on non-target grism frames.
5. G-4 sky-lozenge background template; G-5 focus-offset and CCD-TEMP regressors and an LSF table.
6. F-6: add the background-morphology test to the S2c classifier; build a 200-frame stratified labelled truth set (render thumbnails, label them yourself by inspection, save the labels and thumbnails) and publish the confusion matrix with Wilson intervals; re-issue NGC 5548 slot '6', T CrB slot '6'/'W' verdicts.
7. Re-extract the T CrB validation sample with the corrected library and report how EW and its error moved.
