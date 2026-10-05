<!-- emitted by BeStar_Grism/scripts/be_injection.py; do not edit -->
### Injection–recovery per sample star (BE-S-1b) — before any periodogram

| star | channel | n_nights | baseline_d | p_max_d | nuisance | threshold_power | a90_short | a90_mid | a90_long | frac_periods_closed | amp_bias_signed_all | freq_bias_signed_all | admitted | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| zet Tau | hrg | 20 | 121.790 | 40.597 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass;focpos;ccd_temp | 0.825 | 6.000 | 12.000 | 14.000 | 0.950 | 0.083 | -0.000 | 1 | admitted to slow-tier search |
| HD 45910 | hrg | 18 | 94.761 | 31.587 | offset[S3 ASI post-monsoon (flipped)];h2o_depth;airmass;focpos;ccd_temp | 0.856 | 8.000 | 12.000 | 12.000 | 0.975 | 0.081 | 0.000 | 1 | admitted to slow-tier search |
| HD 34959 | hrg | 17 | 124.792 | 41.597 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass;focpos | 0.889 | 8.000 | 12.000 | 20.000 | 0.950 | 0.054 | 0.001 | 1 | admitted to slow-tier search |
| psi Per | hrg | 15 | 114.854 | 38.285 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass;focpos | 0.922 | 8.000 | 16.000 | 24.000 | 0.950 | 0.048 | -0.001 | 1 | admitted to slow-tier search |
| bet CMi | hrg | 14 | 145.723 | 48.574 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass | 0.948 | 8.000 | 12.000 | 12.000 | 0.950 | 0.037 | -0.000 | 1 | admitted to slow-tier search |
| HD 22298 | hrg | 14 | 108.883 | 36.294 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass | 0.929 | 12.000 | 16.000 | 16.000 | 0.975 | 0.055 | 0.000 | 1 | admitted to slow-tier search |
| HD 46264 | hrg | 13 | 137.858 | 45.953 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass | 0.962 | 12.000 | 16.000 | 16.000 | 0.925 | 0.022 | 0.000 | 1 | admitted to slow-tier search |
| HD 60848 | hrg | 13 | 135.795 | 45.265 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth;airmass | 0.961 | 12.000 | 12.000 | 24.000 | 1.000 | 0.039 | -0.000 | 1 | admitted to slow-tier search |
| eta Tau | hrg | 11 | 62.880 | 20.960 | offset[S3 ASI post-monsoon (flipped)];h2o_depth;airmass | 0.972 | 16.000 | 24.000 | 16.000 | 0.650 | 0.028 | -0.000 | 1 | admitted to slow-tier search |
| HD 698 | hrg | 11 | 50.979 | 16.993 | offset[S3 ASI post-monsoon (flipped)];h2o_depth;airmass | 0.967 | 24.000 | 48.000 | — | 0.550 | 0.020 | 0.001 | 1 | admitted to slow-tier search |
| pi. Aqr | hrg | 10 | 181.362 | 60.454 | offset[S3 ASI post-monsoon (flipped)];offset[S5 QHY re-seated];h2o_depth | 0.993 | 24.000 | 24.000 | 28.000 | 0.450 | 0.019 | 0.000 | 1 | admitted to slow-tier search |
| phi Per | hrg | 7 | 41.982 | 13.994 | offset[S3 ASI post-monsoon (flipped)];h2o_depth | 0.999 | 56.000 | — | — | 0.050 | 0.033 | -0.000 | 0 | 5–9 nights: injected, not searched |
| 12 Aur | hrg | 2 | — | — | — | — | — | — | — | — | — | — | 0 | < 5 standards-epoch nights: no injection, no search |
| ups Cyg | hrg | 2 | — | — | — | — | — | — | — | — | — | — | 0 | < 5 standards-epoch nights: no injection, no search |
| gam Cas | hrg | 1 | — | — | — | — | — | — | — | — | — | — | 0 | < 5 standards-epoch nights: no injection, no search |

Amplitudes in units of the night-to-night scatter σ_night (white-noise injections; BE-S10 converts to Å). `a90_*` = median over the period grid of the smallest amplitude recovered in ≥ 90% of 20 phases, for P < 5 d / 5–20 d / ≥ 20 d (NaN = contour open in that band). Recovery: global peak within 1% of f_inj and power above the FAP = 1% max-statistic threshold (4000 null draws through the same nuisance model). Bias columns are SIGNED means over recovered trials (standing rule 3). Regressors: airmass,focpos,ccd_temp,h2o_depth. Search band 3/T – 1.2 d⁻¹. No science star has ≥ 3 nights of > 2 h span, so there is no short tier (`notes/step0/short_tier.md`). Periodograms of the real data opened at build time: no.

Built 2026-10-05T05:36:18+00:00, seed 20261005.
