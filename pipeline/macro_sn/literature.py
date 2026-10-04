"""macro_sn.literature — every published number this paper uses, with its
source and the sentence it was read from.

Nothing here is a measurement made in this repository.  Each constant
names the paper (and, where it is a quotation, quotes it) so the methods
text, the code and the referee can all check the same line.  Verified
against the papers themselves (arXiv full texts) on 2026-10-04; no value
is from memory.
"""

#: First-light epoch adopted by Li et al. (2025, A&A 703, A168, VizieR
#: J/A+A/703/A168): "All phases are given relative to the time of the first
#: light estimated by Li et al. (2024) (MJD 60082.788)".  The Gate 0 phase
#: zero (macro_sn.gate0.T0_MJD = 60082.79) agrees to 0.002 d.
LI25_T0_MJD = 60082.788
LI25_REF = "Li et al. 2025, A&A 703, A168"

#: The optical maximum.  Two independent statements:
#:  Teja et al. (2023, ApJL 954, L12): the V light curve "reached a peak
#:  V-band magnitude of -18.06 +/- 0.07 mag around ~5 d after explosion".
#:  Li et al. (2025): "The V light curve of SN 2023ixf indicates a rise of
#:  ~6 days before a peak magnitude of M_V = -18.0 +/- 0.1 is reached."
PEAK_V_DAYS = {"Teja et al. 2023, ApJL 954, L12": 5.0,
               "Li et al. 2025, A&A 703, A168": 6.0}

#: The flash-phase H-alpha line flux.  Bostroem et al. (2023, ApJL 956, L5),
#: from their +1.36 d spectrum, "which has the maximum Halpha flux of our
#: spectral series": F_Halpha = 3.18e-13 erg cm^-2 s^-1 (integrated
#: 6300-6800 A after blackbody-continuum subtraction; L = 1.78e39 erg/s at
#: 6.85 Mpc).
F_HALPHA_MAX_CGS = 3.18e-13
F_HALPHA_MAX_PHASE_D = 1.36
F_HALPHA_REF = "Bostroem et al. 2023, ApJL 956, L5"

#: Its decline.  Smith et al. (2023, ApJ 956, 46): "the narrow/intermediate
#: Halpha line luminosity fades by a factor of ~3 during the first week of
#: observations" (and the EW by ~5).
F_HALPHA_FADE_FACTOR = 3.0
F_HALPHA_FADE_SPAN_D = 7.0
F_HALPHA_FADE_REF = "Smith et al. 2023, ApJ 956, 46"

#: Effective width of the PS1 r_P1 bandpass, SVO Filter Profile Service
#: entry PAN-STARRS/PS1.r ("WidthEff" = 1252.4 A; curve of Tonry et al.
#: 2012, ApJ 750, 99).  Used to turn zero-point RATIOS into narrowband
#: widths.  Cached as products/sn/external/svo_ps1_r.xml.
PS1_R_WIDTH_EFF_A = 1252.4068800721
PS1_R_REF = "SVO FPS PAN-STARRS/PS1.r (Tonry et al. 2012)"

#: Plausible narrowband equivalent widths (A), physicist's memo (PH, plan
#: review 2026-10-03): W ~ 65 +/- 15 A from the zero-point ratios.
PLAUSIBLE_WIDTHS_A = (50.0, 65.0, 80.0)

#: Flash-phase epochs the estimator is asked about (d after first light).
FLASH_EPOCHS_D = (2.5, 5.4, 6.4)

#: SDSS -> PS1 bandpass transformations for stellar SEDs, Tonry et al.
#: (2012, ApJ 750, 99) Table 6, LINEAR form y = B0 + B1 x with
#: x = (g - r)_SDSS:  (gP1 - gSDSS): B0 = -0.012, B1 = -0.139 (+/-0.007);
#: (rP1 - rSDSS): 0.000, -0.007 (+/-0.002); (iP1 - iSDSS): 0.004, -0.014
#: (+/-0.003).  "computed for stellar SEDs and use for other SEDs may be
#: less accurate" (Tonry et al.) — applied to the SN with that caveat.
TONRY12_SDSS_TO_PS1 = {"g": (-0.012, -0.139, 0.007),
                       "r": (0.000, -0.007, 0.002),
                       "i": (0.004, -0.014, 0.003)}
TONRY12_REF = "Tonry et al. 2012, ApJ 750, 99, Table 6"

#: Li et al. (2025) gri are SDSS AB magnitudes calibrated to APASS DR9:
#: "calibrated ... [to] the Sloan Digital Sky Survey (SDSS) photometric
#: system (Fukugita et al., 1996) in AB magnitudes ... using the AAVSO
#: Photometric All Sky Survey (APASS) DR9 Catalogue" (PSF photometry, no
#: template subtraction before day 625).
LI25_SYSTEM = "SDSS AB via APASS DR9, PSF photometry"
