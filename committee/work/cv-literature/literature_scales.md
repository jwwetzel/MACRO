# Literature constants and predicted scales (CV paper)

Emitted by `committee/work/cv-literature/literature_scales.py`. Do not edit by hand.

## 1. Constants as read

| name | value | source key | where | how read | note |
|---|---|---|---|---|---|
| `stlmi_vsx_epoch` | 2459298.4236 | `watson2006` | VSX OID 17253 (VizieR B/vsx; vsx.aavso.org), rev. 5 of 2025-11-12 by S. Otero: 'Period and epoch from AAVSO data' | catalogue record | HJD; fiducial NOT documented; no uncertainty published |
| `stlmi_vsx_period` | 0.07908912 | `watson2006` | same record | catalogue record | d; no uncertainty published |
| `stlmi_cropper_epoch` | 2445059.7024 | `cropper1986` | Sect. 3, ephemeris | read in the paper's own text | HJD of the linear-polarisation peak (epoch inherited from Stockman et al. 1983) |
| `stlmi_cropper_epoch_err` | 0.0003 | `cropper1986` | Sect. 3 | read in the paper's own text | d |
| `stlmi_cropper_period` | 0.07908908 | `cropper1986` | Sect. 3 | read in the paper's own text | d; 16 polarisation peaks, 1982-1985 |
| `stlmi_cropper_period_err` | 8e-08 | `cropper1986` | Sect. 3 | read in the paper's own text | d |
| `stlmi_sfst_period` | 0.079087 | `stockman1983` | quoted in Cropper (1986) Sect. 3 | as quoted by a later paper (not read at source) | d, +-0.000001 |
| `stlmi_bws_period` | 0.0790898 | `bailey1985` | Sect. 2.3 | read in the paper's own text | d; phase zero = end of bright phase |
| `stlmi_conj_epoch` | 2450596.9006 | `robertson2008` | ephemeris (period of Cropper 1986, epoch of Howell et al. 2000 via Kafka et al. 2007) | from a fetched page summary (lowest grade; confirm at source before submission) | HJD of inferior conjunction of the secondary |
| `stlmi_conj_minus_cropper_phase` | 0.17 | `kafka2007` | Sect. 2.1 | read in the paper's own text | conjunction phase zero = phase 0.17 on Cropper's polarimetric ephemeris |
| `stlmi_faint_fraction` | 0.64 | `cropper1986` | Sect. 4 | read in the paper's own text | fraction of the orbit the star is faint |
| `stlmi_onset_phase` | 0.69 | `cropper1986` | Sect. 4 | read in the paper's own text | bright-phase onset, polarimetric phase |
| `stlmi_fall_start_phase` | 0.98 | `cropper1986` | Sect. 4 | read in the paper's own text | rapid decline begins |
| `stlmi_fall_end_phase` | 0.05 | `cropper1986` | Sect. 4 | read in the paper's own text | intensity reaches faint level |
| `stlmi_linpol_peak_pct` | 11.5 | `cropper1986` | Sect. 4 | read in the paper's own text | per cent, at phase 0.0, 'half way through the intensity's rapid decline' |
| `stlmi_dur_1982` | 0.33 | `cropper1986` | Sect. 4 | read in the paper's own text | P; bright-phase duration in 1982 (SFST data) |
| `stlmi_dur_1985` | 0.38 | `cropper1986` | Sect. 4 | read in the paper's own text | P; bright-phase duration in 1985 |
| `stlmi_dur_white` | 0.29 | `bailey1985` | Sect. 7.3 | read in the paper's own text | P at half maximum, optical white light |
| `stlmi_dur_J` | 0.31 | `bailey1985` | Sect. 7.3 | read in the paper's own text | P at half maximum |
| `stlmi_dur_H` | 0.36 | `bailey1985` | Sect. 7.3 | read in the paper's own text | P at half maximum |
| `stlmi_dur_K` | 0.39 | `bailey1985` | Sect. 7.3 | read in the paper's own text | P at half maximum |
| `stlmi_dur_opt_1982jun` | 0.24 | `bailey1985` | Sect. 8 | read in the paper's own text | P at half height, optical, 1982 June |
| `stlmi_dur_opt_1983` | 0.29 | `bailey1985` | Sect. 8 | read in the paper's own text | P at half height, optical, 1983 Mar-May |
| `stlmi_fall_width_peacock` | 0.06 | `peacock1992` | results, R_c light curve of 1986 Jan 8 | read in the paper's own text | P_orb taken by the rapid fall (and by the rapid rise) |
| `stlmi_amp_Rc` | 2.0 | `peacock1992` | results, R_c light curve of 1986 Jan 8 | read in the paper's own text | mag, R_c orbital amplitude |
| `stlmi_amp_V` | 1.2 | `kafka2007` | Sect. 3, high-state discussion | read in the paper's own text | mag, V brightening in the bright phase |
| `stlmi_amp_B_quoted` | 0.6 | `peacock1992` | quoted in Campbell et al. (2008) Sect. 1 | as quoted by a later paper (not read at source) | mag, B |
| `stlmi_amp_R_quoted` | 1.4 | `peacock1992` | quoted in Campbell et al. (2008) Sect. 1 | as quoted by a later paper (not read at source) | mag, R |
| `stlmi_amp_I_quoted` | 2.0 | `peacock1992` | quoted in Campbell et al. (2008) Sect. 1 | as quoted by a later paper (not read at source) | mag, I |
| `stlmi_brightphase_width_kafka` | 0.25 | `kafka2007` | Sect. 3, high-state discussion | read in the paper's own text | phase units, V |
| `stlmi_hump_start_rob` | 0.55 | `robertson2008` | photometry, 2007 high state | from a fetched page summary (lowest grade; confirm at source before submission) | conjunction phase |
| `stlmi_peak_rob` | 0.78 | `robertson2008` | photometry, 2007 high state | from a fetched page summary (lowest grade; confirm at source before submission) | conjunction phase of peak |
| `stlmi_drop_rob` | 0.87 | `robertson2008` | photometry, 2007 high state | from a fetched page summary (lowest grade; confirm at source before submission) | conjunction phase of sharp drop |
| `stlmi_B_campbell` | 12.1 | `campbell2008c` | abstract; ST LMi modelling results | read in the paper's own text | MG, +-0.5; i=55 deg, beta=128 deg |
| `stlmi_B_ferrario` | 11.5 | `ferrario1993` | quoted in Kafka et al. (2007) Sect. 1 (+-0.5 MG); '~12.0 MG' in Campbell et al. (2008) | as quoted by a later paper (not read at source) | MG |
| `stlmi_B_second_pole` | 30.0 | `schmidt1983` | quoted in Cropper (1986) Sect. 1 (+-5 MG) and Kafka et al. (2007) | as quoted by a later paper (not read at source) | MG, second (normally non-accreting) pole |
| `stlmi_extent_dense` | 18.0 | `potter2000` | Stokes-imaging solution | from a fetched page summary (lowest grade; confirm at source before submission) | deg of magnetic longitude either side, high-density region; +-15 deg low-density |
| `stlmi_state_duration` | 24.1 | `duffy2022` | Table 1 | read in the paper's own text | d, +-7.6; mean duration of short-lived states |
| `stlmi_state_recurrence` | 273.0 | `duffy2022` | Table 1 | read in the paper's own text | d, +-90; 5 states (3 resolved) |
| `stlmi_extreme_low_V` | 18.5 | `kafka2007` | abstract | read in the paper's own text | mag; 'extreme low state' |
| `stlmi_low_V` | 17.5 | `kafka2007` | Sect. 1 (RoboScope, Kafka & Honeycutt 2005) | read in the paper's own text | mag; 1992-1997 low state; extreme low state V~18.5 |
| `v1500_pdot` | 3.86e-08 | `schmidt1995` | quoted in Pavlenko et al. (2018) Sect. 1 (polarimetry, 1987-1992) | as quoted by a later paper (not read at source) | dimensionless spin-period derivative |
| `v1500_tsync_yr` | 170.0 | `schmidt1995` | quoted in Pavlenko et al. (2018) Sect. 1 | as quoted by a later paper (not read at source) | yr |
| `dpleo_drift_deg` | 50.0 | `beuermann2014` | abstract | read in the paper's own text | deg of azimuth, 1979-2001 |
| `dpleo_lib_period_yr` | 60.0 | `beuermann2014` | abstract | read in the paper's own text | yr (approximate) |
| `dpleo_lib_amp_deg` | 25.0 | `beuermann2014` | abstract | read in the paper's own text | deg (approximate) |
| `yz_vsx_period` | 0.0868 | `shafter1988` | VSX OID 4857 (from Downes et al. 2001); value of Shafter & Hessman (1988) as quoted by Verbunt et al. (1999) and Dai et al. (2026) | as quoted by a later paper (not read at source) | d, +-0.0002; VSX gives no epoch |
| `yz_sh_period_err` | 0.0002 | `shafter1988` | as above | as quoted by a later paper (not read at source) | d |
| `yz_vp_period` | 0.086924 | `vanparadijs1994` | Sect. 6 | read in the paper's own text | d |
| `yz_vp_period_err` | 7e-06 | `vanparadijs1994` | Sect. 6 | read in the paper's own text | d |
| `yz_vp_epoch` | 2447518.7255 | `vanparadijs1994` | Sect. 6 | read in the paper's own text | HJD, +-0.0061; 'probably' superior conjunction of the white dwarf |
| `yz_vp_hump_full_mag` | 0.5 | `vanparadijs1994` | Sect. 6 | read in the paper's own text | mag, full orbital modulation in quiescence, maximum near phase 0.8 |
| `yz_vp_decline_full_mag` | 0.2 | `vanparadijs1994` | Sect. 7 | read in the paper's own text | mag, modulation at late decline, maximum near phase 0.5 |
| `yz_vp_hump_max_phase` | 0.8 | `vanparadijs1994` | Sect. 6 | read in the paper's own text | spectroscopic phase of hump maximum in quiescence |
| `yz_vp_decline_max_phase` | 0.5 | `vanparadijs1994` | Sect. 7 | read in the paper's own text | phase of maximum at late decline |
| `yz_vp_super_peak` | 10.6 | `vanparadijs1994` | Sect. 4 | read in the paper's own text | mag at superoutburst peak (1988-89) |
| `yz_vp_normal_peak` | 12.0 | `vanparadijs1994` | Sect. 4 | read in the paper's own text | mag at normal-outburst peak; ~3 d long; intervals 9-12 d |
| `yz_pat_tn` | 10.3 | `patterson1979` | Table I | read in the paper's own text | d, +-2.5; recurrence of short maxima |
| `yz_pat_ts` | 134.0 | `patterson1979` | Table I | read in the paper's own text | d, +-19; recurrence of long maxima |
| `yz_pat_super_minus_normal` | 1.0 | `patterson1979` | p. 804 | read in the paper's own text | mag; supermaxima ~1 mag brighter and ~5x longer |
| `yz_kato_max` | 11.0 | `kato2002` | Table 3 | read in the paper's own text | mag at (super)maximum |
| `yz_kato_min` | 15.0 | `kato2002` | Table 3 | read in the paper's own text | mag at minimum |
| `yz_hakala_quiescence` | 15.2 | `hakala2004` | observations section | read in the paper's own text | V mag in quiescence (AAVSO CCD V) |
| `yz_hakala_outburst` | 12.2 | `hakala2004` | observations section | read in the paper's own text | V mag reached in the normal outburst 4 d later |
| `yz_flicker_pp` | 0.75 | `moffett1974` | quoted in van Paradijs et al. (1994) Sect. 6 and Zhao et al. (2005) | as quoted by a later paper (not read at source) | mag peak to peak |
| `yz_psh_2007` | 0.09031 | `kato2009` | Sect. 6.34 | read in the paper's own text | d, +-0.00005 (timing); 0.09042(4) by PDM; 2007 Feb |
| `yz_eps_2007_pct` | 4.0 | `kato2009` | Sect. 6.34 | read in the paper's own text | per cent fractional superhump excess |
| `yz_psh_tess` | 0.09043 | `dai2026` | abstract | read in the paper's own text | d, +-0.00027; excess 4.03(31) per cent |
| `yz_sh_amp_tess_early` | 0.3 | `dai2026` | Sect. 4.2 | read in the paper's own text | mag (full), precursor and early plateau |
| `yz_sh_amp_tess_late` | 0.5 | `dai2026` | Sect. 4.2 | read in the paper's own text | mag (full), post-plateau |
| `sh_amp_max_lowi` | 0.25 | `smak2010` | abstract | read in the paper's own text | mag (full), maximum superhump amplitude at low inclination |
| `sh_amp_max_kato` | 0.25 | `kato2012` | Sect. 4.7 | read in the paper's own text | mag; sample mean of maximum superhump amplitude |
| `yz_tess_tn_lo` | 8.08 | `sun2026` | Sect. 3.16 | read in the paper's own text | d; mean normal-outburst recurrence, S71-S72 |
| `yz_tess_tn_hi` | 9.49 | `sun2026` | Sect. 3.16 | read in the paper's own text | d; mean normal-outburst recurrence, S44-S47 |
| `vvpup_epoch` | 2427889.6474 | `walker1965` | VSX OID 26642 remark 'The epoch of light maximum is given'; ephemeris reproduced in Howell et al. (2006) Sect. 3 | as quoted by a later paper (not read at source) | HJD of high-state maximum light |
| `vvpup_period` | 0.0697468256 | `walker1965` | same | as quoted by a later paper (not read at source) | d |
| `anuma_vsx_epoch` | 2456725.4053 | `watson2006` | VSX OID 37171, rev. 4 of 2025-11-12: 'Epoch from AAVSO data'; remark: 'The epoch of Min ... is given' | catalogue record | HJD of minimum light |
| `anuma_vsx_period` | 0.07975274 | `watson2006` | VSX OID 37171 (GCVS-era value) | catalogue record | d |
| `anuma_bb_period` | 0.07975282 | `bonnetbidaud1996` | abstract: (6890.6436 +- 0.0035) s; ephemeris quoted in Ok et al. (2025) Sect. 2: HJD 2443190.9921(2) | read in the paper's own text | d, +-0.00000004; T0 = linear-polarisation pulse of Liebert et al. (1982) |
| `anuma_ok_period` | 0.079752867 | `ok2025` | Eq. 1 | read in the paper's own text | d, +-0.000000012; T0(BJD) = 2443190.9926(2) |
| `euuma_vsx_epoch` | 2458231.7799 | `chen2020` | VSX OID 37268, rev. 5 of 2021-10-26: 'Updated from Chen 2020 catalog' | catalogue record | HJD; Chen et al. define T0 as the time of MINIMUM of the fitted ZTF light curve |
| `euuma_vsx_period` | 0.0626194 | `chen2020` | same | catalogue record | d |
| `cyc_lambda_100MG_um` | 1.071 | `wickramasinghe2000` | cyclotron fundamental, lambda_c = 10710 A (10^8 G / B) | read in the paper's own text | micron at B = 100 MG (= 2 pi m_e c^2 / e B) |

## 2. Derived (arithmetic in `derive()`)

| name | value | computed as |
|---|---|---|
| `vsx_minus_cropper_sigma` | 0.5 | (P_VSX - P_Cropper) / sigma_Cropper |
| `ours_minus_cropper_sigma` | -0.585652 | (P_refit - P_Cropper) / sqrt(sigma_refit^2 + sigma_Cropper^2) |
| `cropper_drift_cycles` | 0.0221209 | NumStLmiCycles x sigma_P(Cropper) / P: phase drift over the paper's cycle count under the only PUBLISHED period uncertainty |
| `cropper_drift_margin` | 22.6031 | 0.5 cycle / that drift |
| `cycles_cropper_to_vsx` | 180034 | (T_VSX - T_Cropper) / P_Cropper |
| `vsx_zero_on_cropper` | 0.97182 | fractional part of the above: polarimetric phase of VSX phase zero |
| `vsx_zero_on_cropper_err` | 0.182108 | cycles x sigma_P / P (1 sigma; epoch error and HJD-vs-BJD, <0.01 cycle, neglected) |
| `vsx_zero_on_conj` | 0.800734 | same on the inferior-conjunction ephemeris (Cropper period) |
| `vsx_zero_on_conj_err` | 0.111289 | cycles x sigma_P / P |
| `edge_on_cropper` | 0.12882 | VSX zero + NumStLmiOcOffsetCycles, polarimetric phase |
| `edge_on_conj` | 0.957734 | VSX zero + NumStLmiOcOffsetCycles, conjunction phase |
| `cycles_cropper_to_edge` | 201903 | cycles from Cropper's epoch to our last timed edge |
| `baseline_cropper_yr` | 43.7189 | that baseline in years |
| `pdot_gr_lo` | 4.33068e-14 | P / (5 Gyr): secular orbital evolution, slow end |
| `pdot_gr_hi` | 2.16534e-13 | P / (1 Gyr): secular orbital evolution, fast end |
| `v1500_over_bound` | 10.7222 | Pdot_spin(V1500 Cyg) / NumStLmiPdotLimit |
| `libration_pdot` | 3.57065e-11 | peak curvature of a DP Leo-like libration (A=25 deg, T=60 yr) expressed as an apparent Pdot at ST LMi's period: A (2 pi/T)^2 / 360 x P^2 |
| `bound_over_libration` | 100.822 | NumStLmiPdotLimit / that |
| `dpleo_rate_deg_yr` | 2.27273 | 50 deg / (2001 - 1979) |
| `rms_deg` | 4.42539 | NumStLmiOcRmsS / P x 360: O-C scatter as spot longitude |
| `rms_cycles` | 0.0122927 | NumStLmiOcRmsS / P |
| `fall_s_lo` | 409.998 | 0.06 P (Peacock et al. 1992) |
| `fall_s_hi` | 478.331 | (1.05 - 0.98) P (Cropper 1986) |
| `season_shift_s` | 170.832 | half the 1982->1985 change in bright-phase duration: how far one edge moves if the change is symmetric |
| `season_shift_deg` | 9 | the same in degrees |
| `chromatic_white_J_s` | 68.333 | half the white-light -> J difference in duration (Bailey et al. 1985) |
| `chromatic_J_K_s` | 273.332 | half the J -> K difference |
| `chromatic_g_i_s` | 37.1944 | the white->J edge shift scaled by ln(lambda_i/lambda_g)/ln(1.25/0.55): a log-linear interpolation, ORDER OF MAGNITUDE ONLY; sign: redder band ends later |
| `cyc_fundamental_um` | 8.85124 | 1.0710 um x (100 MG / 12.1 MG) |
| `harmonic_g` | 18.4017 | lambda_c / lambda_g |
| `harmonic_r` | 14.3456 | lambda_c / lambda_r |
| `harmonic_i` | 11.7703 | lambda_c / lambda_i |
| `yz_day_drift_sh` | 0.0265455 | phase drift over 1 d under Shafter & Hessman's published sigma_P (cycles) |
| `yz_day_drift_vp` | 0.000926443 | the same under van Paradijs et al.'s sigma_P |
| `yz_day_shift_between_periods` | 0.0164582 | phase slip per day from folding on 0.0868 d instead of 0.086924 d |
| `yz_vp_epoch_drift_2024` | 11.9338 | cycles of phase uncertainty in propagating the 1988 epoch to 2024: absolute phase is lost |
| `yz_super_amp` | 4 | 15.0 - 11.0 (Kato et al. 2002) |
| `yz_normal_amp` | 3 | 15.2 - 12.2 (Hakala et al. 2004) |
| `yz_super_minus_normal_vp` | 1.4 | 12.0 - 10.6 (van Paradijs et al. 1994) |
| `sh_semi_peak_mmag` | 125 | half of 0.25 mag, in mmag |
| `yz_sh_semi_tess_mmag` | 150 | half of 0.3 mag, in mmag |
| `yz_hump_semi_vp_mmag` | 250 | half of 0.5 mag, in mmag |
| `yz_decline_semi_vp_mmag` | 100 | half of 0.2 mag, in mmag |
