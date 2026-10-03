# Package `detector` — measured gain, read noise, linearity and saturation per camera
Owns: `pipeline/rlmt_diagnostics/{ptc,noise,linearity,ceiling}.py` (and new modules there except `dispersion.py`), `pipeline/scripts/run_s2_campaign.py`, their tests, `products/detector/` (scratch outputs), `docs/pipeline/s2_*.html` generators. Findings: F-4, F-5 (DE.F1–F4, F7–F9, DS.F7, OA.E3, OA.E8), SN-S2-linearity, TCRB-P0-gain-ptc.
Do NOT write the shared manifest in this wave: compute into `products/detector/detector.sqlite` with the same table schemas, and make the S2 stage able to write them into the manifest at the integrated rebuild.
Do:
1. Flat-pair photon-transfer per camera/configuration (AC4040 High Gain both EGAIN epochs, StackPro, iKon 1 MHz, ASI Mode0 2×2, QHY if any flats exist) with read noise from bias pairs; K ± error, RN, effective full scale. Decide D2 with evidence (is Mode0 ≈ 1.0 e⁻/ADU; is binning an average).
2. Fix `s2_noise_pairs` alias/same-DATE-OBS pairing (DE.F4).
3. Linearity: residual-vs-own-peak for comparison stars per mode up to each veto; High Gain twilight-flat ramp; the 2023-06-07 `cmos_tests/Albireo` set; SN 0.5 s vs 2 s pairs; StackPro vs High Gain flux ratio. Recommend ONE cap per mode with the slope criterion (< 1%).
4. Bad-pixel mask per camera (the 16.5k/33k/49k average-binned rail pixels; hot/RTS pixels), stored as a product, with an API other packages can call.
5. Scripted peak-at-target saturation census for all 224 T CrB imaging frames (OA.E3), native-pixel-equivalent.
6. A regenerated S2 evidence page section and figures (PTC curves, linearity residuals).
