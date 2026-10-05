#!/usr/bin/env python
"""Literature constants for the CV paper, and the arithmetic done on them -- emitted, never typed.

WHY THIS SCRIPT EXISTS
----------------------
The "Comparison with previous work" subsection (``comparison_section.tex``) and the answers to the
committee's "(verify)" items (``verify_items.md``) quote two kinds of number that are not measurements of
ours:

* **literature constants** -- an ephemeris, a field strength, a bright-phase duration -- each read from a
  named paper at a named place; and
* **predicted scales** -- arithmetic on those constants (how many cycles a published period error
  accumulates, what Pdot a DP Leo-like libration would mimic, which cyclotron harmonic lands in g).

The house rule is that no number reaches a manuscript by hand.  For measurements the rule is enforced by
the products database; for these there is no database, so this file IS the record: every constant is
declared once below with its BibTeX key, where in the paper it was read, and HOW it was read (the paper's
own text, or a quotation of it in a later paper -- the distinction is carried, because a second-hand
value is weaker evidence).  The arithmetic is done here and nowhere else.

OUTPUTS (both overwritten on every run; do not edit by hand)
------------------------------------------------------------
``literature_macros.tex``  ``\\newcommand{\\Lit...}`` macros used by ``comparison_section.tex``; each line
                           carries its source in a trailing comment, in the style of ``numbers.tex``.
``literature_scales.md``   the same values as a table with source, locator and provenance grade.

INPUT FROM THE PROJECT (read-only)
----------------------------------
Four of the derived scales need a value of OURS (the cycle count to the VSX epoch, the measured edge
offset, the refitted period and its error).  They are parsed read-only from the generated
``manuscripts/CV_TimeSeries/numbers.tex`` so that they track the pipeline; nothing is written there.
If CV-R3 changes those macros, re-run this script.

USAGE
-----
    cd <repo>
    /opt/miniconda3/envs/rlmt-checks/bin/python committee/work/cv-literature/literature_scales.py
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
NUMBERS_TEX = REPO / "manuscripts" / "CV_TimeSeries" / "numbers.tex"

SEC_PER_DAY = 86400.0
SEC_PER_YEAR = 365.25 * SEC_PER_DAY

#: Provenance grades, strongest first.
DIRECT = "read in the paper's own text"
QUOTED = "as quoted by a later paper (not read at source)"
CATALOGUE = "catalogue record"
SUMMARY = "from a fetched page summary (lowest grade; confirm at source before submission)"


@dataclass(frozen=True)
class Lit:
    """One literature constant: a value, where it comes from, and how it was read."""
    value: float | str
    key: str            # BibTeX key in references.bib
    where: str          # locator inside the source
    how: str            # one of the provenance grades above
    note: str = ""


# ==================================================================================================
# 1. LITERATURE CONSTANTS
# ==================================================================================================
C: dict[str, Lit] = {
    # ---- ST LMi: ephemerides -------------------------------------------------------------------
    "stlmi_vsx_epoch": Lit(2459298.4236, "watson2006", "VSX OID 17253 (VizieR B/vsx; vsx.aavso.org), "
                           "rev. 5 of 2025-11-12 by S. Otero: 'Period and epoch from AAVSO data'", CATALOGUE,
                           "HJD; fiducial NOT documented; no uncertainty published"),
    "stlmi_vsx_period": Lit(0.07908912, "watson2006", "same record", CATALOGUE, "d; no uncertainty published"),
    "stlmi_cropper_epoch": Lit(2445059.7024, "cropper1986", "Sect. 3, ephemeris", DIRECT,
                               "HJD of the linear-polarisation peak (epoch inherited from Stockman et al. 1983)"),
    "stlmi_cropper_epoch_err": Lit(0.0003, "cropper1986", "Sect. 3", DIRECT, "d"),
    "stlmi_cropper_period": Lit(0.07908908, "cropper1986", "Sect. 3", DIRECT, "d; 16 polarisation peaks, 1982-1985"),
    "stlmi_cropper_period_err": Lit(0.00000008, "cropper1986", "Sect. 3", DIRECT, "d"),
    "stlmi_sfst_period": Lit(0.079087, "stockman1983", "quoted in Cropper (1986) Sect. 3", QUOTED, "d, +-0.000001"),
    "stlmi_bws_period": Lit(0.0790898, "bailey1985", "Sect. 2.3", DIRECT, "d; phase zero = end of bright phase"),
    "stlmi_conj_epoch": Lit(2450596.9006, "robertson2008", "ephemeris (period of Cropper 1986, epoch of Howell "
                            "et al. 2000 via Kafka et al. 2007)", SUMMARY, "HJD of inferior conjunction of the secondary"),
    "stlmi_conj_minus_cropper_phase": Lit(0.17, "kafka2007", "Sect. 2.1", DIRECT,
                                          "conjunction phase zero = phase 0.17 on Cropper's polarimetric ephemeris"),
    # ---- ST LMi: the bright phase ---------------------------------------------------------------
    "stlmi_faint_fraction": Lit(0.64, "cropper1986", "Sect. 4", DIRECT, "fraction of the orbit the star is faint"),
    "stlmi_onset_phase": Lit(0.69, "cropper1986", "Sect. 4", DIRECT, "bright-phase onset, polarimetric phase"),
    "stlmi_fall_start_phase": Lit(0.98, "cropper1986", "Sect. 4", DIRECT, "rapid decline begins"),
    "stlmi_fall_end_phase": Lit(0.05, "cropper1986", "Sect. 4", DIRECT, "intensity reaches faint level"),
    "stlmi_linpol_peak_pct": Lit(11.5, "cropper1986", "Sect. 4", DIRECT,
                                 "per cent, at phase 0.0, 'half way through the intensity's rapid decline'"),
    "stlmi_dur_1982": Lit(0.33, "cropper1986", "Sect. 4", DIRECT, "P; bright-phase duration in 1982 (SFST data)"),
    "stlmi_dur_1985": Lit(0.38, "cropper1986", "Sect. 4", DIRECT, "P; bright-phase duration in 1985"),
    "stlmi_dur_white": Lit(0.29, "bailey1985", "Sect. 7.3", DIRECT, "P at half maximum, optical white light"),
    "stlmi_dur_J": Lit(0.31, "bailey1985", "Sect. 7.3", DIRECT, "P at half maximum"),
    "stlmi_dur_H": Lit(0.36, "bailey1985", "Sect. 7.3", DIRECT, "P at half maximum"),
    "stlmi_dur_K": Lit(0.39, "bailey1985", "Sect. 7.3", DIRECT, "P at half maximum"),
    "stlmi_dur_opt_1982jun": Lit(0.24, "bailey1985", "Sect. 8", DIRECT, "P at half height, optical, 1982 June"),
    "stlmi_dur_opt_1983": Lit(0.29, "bailey1985", "Sect. 8", DIRECT, "P at half height, optical, 1983 Mar-May"),
    "stlmi_fall_width_peacock": Lit(0.06, "peacock1992", "results, R_c light curve of 1986 Jan 8", DIRECT,
                                    "P_orb taken by the rapid fall (and by the rapid rise)"),
    "stlmi_amp_Rc": Lit(2.0, "peacock1992", "results, R_c light curve of 1986 Jan 8", DIRECT, "mag, R_c orbital amplitude"),
    "stlmi_amp_V": Lit(1.2, "kafka2007", "Sect. 3, high-state discussion", DIRECT, "mag, V brightening in the bright phase"),
    "stlmi_amp_B_quoted": Lit(0.6, "peacock1992", "quoted in Campbell et al. (2008) Sect. 1", QUOTED, "mag, B"),
    "stlmi_amp_R_quoted": Lit(1.4, "peacock1992", "quoted in Campbell et al. (2008) Sect. 1", QUOTED, "mag, R"),
    "stlmi_amp_I_quoted": Lit(2.0, "peacock1992", "quoted in Campbell et al. (2008) Sect. 1", QUOTED, "mag, I"),
    "stlmi_brightphase_width_kafka": Lit(0.25, "kafka2007", "Sect. 3, high-state discussion", DIRECT, "phase units, V"),
    "stlmi_hump_start_rob": Lit(0.55, "robertson2008", "photometry, 2007 high state", SUMMARY, "conjunction phase"),
    "stlmi_peak_rob": Lit(0.78, "robertson2008", "photometry, 2007 high state", SUMMARY, "conjunction phase of peak"),
    "stlmi_drop_rob": Lit(0.87, "robertson2008", "photometry, 2007 high state", SUMMARY, "conjunction phase of sharp drop"),
    # ---- ST LMi: field, geometry, states --------------------------------------------------------
    "stlmi_B_campbell": Lit(12.1, "campbell2008c", "abstract; ST LMi modelling results", DIRECT, "MG, +-0.5; i=55 deg, beta=128 deg"),
    "stlmi_B_ferrario": Lit(11.5, "ferrario1993", "quoted in Kafka et al. (2007) Sect. 1 (+-0.5 MG); "
                            "'~12.0 MG' in Campbell et al. (2008)", QUOTED, "MG"),
    "stlmi_B_second_pole": Lit(30.0, "schmidt1983", "quoted in Cropper (1986) Sect. 1 (+-5 MG) and Kafka et al. (2007)",
                               QUOTED, "MG, second (normally non-accreting) pole"),
    "stlmi_extent_dense": Lit(18.0, "potter2000", "Stokes-imaging solution", SUMMARY,
                              "deg of magnetic longitude either side, high-density region; +-15 deg low-density"),
    "stlmi_state_duration": Lit(24.1, "duffy2022", "Table 1", DIRECT, "d, +-7.6; mean duration of short-lived states"),
    "stlmi_state_recurrence": Lit(273.0, "duffy2022", "Table 1", DIRECT, "d, +-90; 5 states (3 resolved)"),
    "stlmi_extreme_low_V": Lit(18.5, "kafka2007", "abstract", DIRECT, "mag; 'extreme low state'"),
    "stlmi_low_V": Lit(17.5, "kafka2007", "Sect. 1 (RoboScope, Kafka & Honeycutt 2005)", DIRECT,
                       "mag; 1992-1997 low state; extreme low state V~18.5"),
    # ---- Spin-orbit synchronism and spot-longitude scales ---------------------------------------
    "v1500_pdot": Lit(3.86e-8, "schmidt1995", "quoted in Pavlenko et al. (2018) Sect. 1 (polarimetry, 1987-1992)",
                      QUOTED, "dimensionless spin-period derivative"),
    "v1500_tsync_yr": Lit(170.0, "schmidt1995", "quoted in Pavlenko et al. (2018) Sect. 1", QUOTED, "yr"),
    "dpleo_drift_deg": Lit(50.0, "beuermann2014", "abstract", DIRECT, "deg of azimuth, 1979-2001"),
    "dpleo_lib_period_yr": Lit(60.0, "beuermann2014", "abstract", DIRECT, "yr (approximate)"),
    "dpleo_lib_amp_deg": Lit(25.0, "beuermann2014", "abstract", DIRECT, "deg (approximate)"),
    # ---- YZ Cnc ---------------------------------------------------------------------------------
    "yz_vsx_period": Lit(0.0868, "shafter1988", "VSX OID 4857 (from Downes et al. 2001); value of Shafter & Hessman "
                         "(1988) as quoted by Verbunt et al. (1999) and Dai et al. (2026)", QUOTED, "d, +-0.0002; VSX gives no epoch"),
    "yz_sh_period_err": Lit(0.0002, "shafter1988", "as above", QUOTED, "d"),
    "yz_vp_period": Lit(0.086924, "vanparadijs1994", "Sect. 6", DIRECT, "d"),
    "yz_vp_period_err": Lit(0.000007, "vanparadijs1994", "Sect. 6", DIRECT, "d"),
    "yz_vp_epoch": Lit(2447518.7255, "vanparadijs1994", "Sect. 6", DIRECT,
                       "HJD, +-0.0061; 'probably' superior conjunction of the white dwarf"),
    "yz_vp_hump_full_mag": Lit(0.5, "vanparadijs1994", "Sect. 6", DIRECT, "mag, full orbital modulation in quiescence, maximum near phase 0.8"),
    "yz_vp_decline_full_mag": Lit(0.2, "vanparadijs1994", "Sect. 7", DIRECT, "mag, modulation at late decline, maximum near phase 0.5"),
    "yz_vp_hump_max_phase": Lit(0.8, "vanparadijs1994", "Sect. 6", DIRECT, "spectroscopic phase of hump maximum in quiescence"),
    "yz_vp_decline_max_phase": Lit(0.5, "vanparadijs1994", "Sect. 7", DIRECT, "phase of maximum at late decline"),
    "yz_vp_super_peak": Lit(10.6, "vanparadijs1994", "Sect. 4", DIRECT, "mag at superoutburst peak (1988-89)"),
    "yz_vp_normal_peak": Lit(12.0, "vanparadijs1994", "Sect. 4", DIRECT, "mag at normal-outburst peak; ~3 d long; intervals 9-12 d"),
    "yz_pat_tn": Lit(10.3, "patterson1979", "Table I", DIRECT, "d, +-2.5; recurrence of short maxima"),
    "yz_pat_ts": Lit(134.0, "patterson1979", "Table I", DIRECT, "d, +-19; recurrence of long maxima"),
    "yz_pat_super_minus_normal": Lit(1.0, "patterson1979", "p. 804", DIRECT, "mag; supermaxima ~1 mag brighter and ~5x longer"),
    "yz_kato_max": Lit(11.0, "kato2002", "Table 3", DIRECT, "mag at (super)maximum"),
    "yz_kato_min": Lit(15.0, "kato2002", "Table 3", DIRECT, "mag at minimum"),
    "yz_hakala_quiescence": Lit(15.2, "hakala2004", "observations section", DIRECT, "V mag in quiescence (AAVSO CCD V)"),
    "yz_hakala_outburst": Lit(12.2, "hakala2004", "observations section", DIRECT, "V mag reached in the normal outburst 4 d later"),
    "yz_flicker_pp": Lit(0.75, "moffett1974", "quoted in van Paradijs et al. (1994) Sect. 6 and Zhao et al. (2005)", QUOTED, "mag peak to peak"),
    "yz_psh_2007": Lit(0.09031, "kato2009", "Sect. 6.34", DIRECT, "d, +-0.00005 (timing); 0.09042(4) by PDM; 2007 Feb"),
    "yz_eps_2007_pct": Lit(4.0, "kato2009", "Sect. 6.34", DIRECT, "per cent fractional superhump excess"),
    "yz_psh_tess": Lit(0.09043, "dai2026", "abstract", DIRECT, "d, +-0.00027; excess 4.03(31) per cent"),
    "yz_sh_amp_tess_early": Lit(0.3, "dai2026", "Sect. 4.2", DIRECT, "mag (full), precursor and early plateau"),
    "yz_sh_amp_tess_late": Lit(0.5, "dai2026", "Sect. 4.2", DIRECT, "mag (full), post-plateau"),
    "sh_amp_max_lowi": Lit(0.25, "smak2010", "abstract", DIRECT, "mag (full), maximum superhump amplitude at low inclination"),
    "sh_amp_max_kato": Lit(0.25, "kato2012", "Sect. 4.7", DIRECT, "mag; sample mean of maximum superhump amplitude"),
    "yz_tess_tn_lo": Lit(8.08, "sun2026", "Sect. 3.16", DIRECT, "d; mean normal-outburst recurrence, S71-S72"),
    "yz_tess_tn_hi": Lit(9.49, "sun2026", "Sect. 3.16", DIRECT, "d; mean normal-outburst recurrence, S44-S47"),
    # ---- The other three polars: ephemerides of record ------------------------------------------
    "vvpup_epoch": Lit(2427889.6474, "walker1965", "VSX OID 26642 remark 'The epoch of light maximum is given'; "
                       "ephemeris reproduced in Howell et al. (2006) Sect. 3", QUOTED, "HJD of high-state maximum light"),
    "vvpup_period": Lit(0.0697468256, "walker1965", "same", QUOTED, "d"),
    "anuma_vsx_epoch": Lit(2456725.4053, "watson2006", "VSX OID 37171, rev. 4 of 2025-11-12: 'Epoch from AAVSO data'; "
                           "remark: 'The epoch of Min ... is given'", CATALOGUE, "HJD of minimum light"),
    "anuma_vsx_period": Lit(0.07975274, "watson2006", "VSX OID 37171 (GCVS-era value)", CATALOGUE, "d"),
    "anuma_bb_period": Lit(0.07975282, "bonnetbidaud1996", "abstract: (6890.6436 +- 0.0035) s; ephemeris quoted in "
                           "Ok et al. (2025) Sect. 2: HJD 2443190.9921(2)", DIRECT, "d, +-0.00000004; T0 = linear-polarisation pulse of Liebert et al. (1982)"),
    "anuma_ok_period": Lit(0.079752867, "ok2025", "Eq. 1", DIRECT, "d, +-0.000000012; T0(BJD) = 2443190.9926(2)"),
    "euuma_vsx_epoch": Lit(2458231.7799, "chen2020", "VSX OID 37268, rev. 5 of 2021-10-26: 'Updated from Chen 2020 catalog'",
                           CATALOGUE, "HJD; Chen et al. define T0 as the time of MINIMUM of the fitted ZTF light curve"),
    "euuma_vsx_period": Lit(0.0626194, "chen2020", "same", CATALOGUE, "d"),
    # ---- Physical constant ----------------------------------------------------------------------
    "cyc_lambda_100MG_um": Lit(1.0710, "wickramasinghe2000", "cyclotron fundamental, lambda_c = 10710 A (10^8 G / B)", DIRECT,
                               "micron at B = 100 MG (= 2 pi m_e c^2 / e B)"),
}

#: Nominal effective wavelengths of the Pan-STARRS1/REFCAT2 system the photometry is tied to (micron).
BAND_UM = {"g": 0.481, "r": 0.617, "i": 0.752}


# ==================================================================================================
# 2. OUR OWN VALUES, read-only from the generated macro file
# ==================================================================================================
def project_macro(name: str) -> float:
    """Numeric value of ``\\newcommand{\\<name>}{...}`` in numbers.tex (thin spaces and x10^n handled)."""
    m = re.search(r"\\newcommand\{\\" + name + r"\}\{(.*?)\}\s*%", NUMBERS_TEX.read_text(encoding="utf-8"))
    if not m:
        raise KeyError(f"{name} not found in {NUMBERS_TEX}")
    s = m.group(1).replace(r"\,", "").replace(" ", "")
    sci = re.fullmatch(r"([-\d.]+)\\times10\^\{(-?\d+)\}", s)
    return float(sci.group(1)) * 10 ** int(sci.group(2)) if sci else float(s)


# ==================================================================================================
# 3. THE ARITHMETIC
# ==================================================================================================
def derive() -> dict[str, tuple[float, str]]:
    """Return {name: (value, how it was computed)} for every predicted scale."""
    v = {k: c.value for k, c in C.items()}
    P_d = v["stlmi_vsx_period"]
    P_s = P_d * SEC_PER_DAY
    ours_cycles = project_macro("NumStLmiCycles")
    ours_offset = project_macro("NumStLmiOcOffsetCycles")
    ours_P = project_macro("NumStLmiFittedPeriodD")
    ours_sP = project_macro("NumStLmiFittedPeriodSigmaD")
    ours_pdot = project_macro("NumStLmiPdotLimit")
    ours_rms = project_macro("NumStLmiOcRmsS")
    d: dict[str, tuple[float, str]] = {}

    # --- period comparisons -----------------------------------------------------------------------
    d["vsx_minus_cropper_sigma"] = ((P_d - v["stlmi_cropper_period"]) / v["stlmi_cropper_period_err"],
                                    "(P_VSX - P_Cropper) / sigma_Cropper")
    d["ours_minus_cropper_sigma"] = ((ours_P - v["stlmi_cropper_period"])
                                     / math.hypot(ours_sP, v["stlmi_cropper_period_err"]),
                                     "(P_refit - P_Cropper) / sqrt(sigma_refit^2 + sigma_Cropper^2)")
    d["cropper_drift_cycles"] = (ours_cycles * v["stlmi_cropper_period_err"] / P_d,
                                 "NumStLmiCycles x sigma_P(Cropper) / P: phase drift over the paper's cycle count "
                                 "under the only PUBLISHED period uncertainty")
    d["cropper_drift_margin"] = (0.5 / d["cropper_drift_cycles"][0], "0.5 cycle / that drift")

    # --- where VSX phase zero, and our timed edge, fall on the two published ephemerides ----------
    n_c = (v["stlmi_vsx_epoch"] - v["stlmi_cropper_epoch"]) / v["stlmi_cropper_period"]
    n_h = (v["stlmi_vsx_epoch"] - v["stlmi_conj_epoch"]) / v["stlmi_cropper_period"]
    d["cycles_cropper_to_vsx"] = (n_c, "(T_VSX - T_Cropper) / P_Cropper")
    d["vsx_zero_on_cropper"] = (n_c % 1.0, "fractional part of the above: polarimetric phase of VSX phase zero")
    d["vsx_zero_on_cropper_err"] = (n_c * v["stlmi_cropper_period_err"] / v["stlmi_cropper_period"],
                                    "cycles x sigma_P / P (1 sigma; epoch error and HJD-vs-BJD, <0.01 cycle, neglected)")
    d["vsx_zero_on_conj"] = (n_h % 1.0, "same on the inferior-conjunction ephemeris (Cropper period)")
    d["vsx_zero_on_conj_err"] = (n_h * v["stlmi_cropper_period_err"] / v["stlmi_cropper_period"], "cycles x sigma_P / P")
    d["edge_on_cropper"] = ((n_c + ours_offset) % 1.0, "VSX zero + NumStLmiOcOffsetCycles, polarimetric phase")
    d["edge_on_conj"] = ((n_h + ours_offset) % 1.0, "VSX zero + NumStLmiOcOffsetCycles, conjunction phase")
    d["cycles_cropper_to_edge"] = (n_c + ours_cycles * P_d / v["stlmi_cropper_period"],
                                   "cycles from Cropper's epoch to our last timed edge")
    d["baseline_cropper_yr"] = (d["cycles_cropper_to_edge"][0] * P_d / 365.25, "that baseline in years")

    # --- the three physical scales for the Pdot bound (PH.P2) -------------------------------------
    d["pdot_gr_lo"] = (P_s / (5e9 * SEC_PER_YEAR), "P / (5 Gyr): secular orbital evolution, slow end")
    d["pdot_gr_hi"] = (P_s / (1e9 * SEC_PER_YEAR), "P / (1 Gyr): secular orbital evolution, fast end")
    d["v1500_over_bound"] = (v["v1500_pdot"] / ours_pdot, "Pdot_spin(V1500 Cyg) / NumStLmiPdotLimit")
    acc = v["dpleo_lib_amp_deg"] * (2 * math.pi / v["dpleo_lib_period_yr"]) ** 2       # deg / yr^2, peak
    d["libration_pdot"] = (acc / 360.0 * P_s ** 2 / SEC_PER_YEAR ** 2,
                           "peak curvature of a DP Leo-like libration (A=25 deg, T=60 yr) expressed as an "
                           "apparent Pdot at ST LMi's period: A (2 pi/T)^2 / 360 x P^2")
    d["bound_over_libration"] = (ours_pdot / d["libration_pdot"][0], "NumStLmiPdotLimit / that")
    d["dpleo_rate_deg_yr"] = (v["dpleo_drift_deg"] / 22.0, "50 deg / (2001 - 1979)")
    d["rms_deg"] = (ours_rms / P_s * 360.0, "NumStLmiOcRmsS / P x 360: O-C scatter as spot longitude")
    d["rms_cycles"] = (ours_rms / P_s, "NumStLmiOcRmsS / P")

    # --- scales for the timed edge itself (PH.P3) -------------------------------------------------
    d["fall_s_lo"] = (v["stlmi_fall_width_peacock"] * P_s, "0.06 P (Peacock et al. 1992)")
    d["fall_s_hi"] = (((v["stlmi_fall_end_phase"] + 1.0) - v["stlmi_fall_start_phase"]) * P_s,
                      "(1.05 - 0.98) P (Cropper 1986)")
    d["season_shift_s"] = (0.5 * (v["stlmi_dur_1985"] - v["stlmi_dur_1982"]) * P_s,
                           "half the 1982->1985 change in bright-phase duration: how far one edge moves "
                           "if the change is symmetric")
    d["season_shift_deg"] = (0.5 * (v["stlmi_dur_1985"] - v["stlmi_dur_1982"]) * 360.0, "the same in degrees")
    d["chromatic_white_J_s"] = (0.5 * (v["stlmi_dur_J"] - v["stlmi_dur_white"]) * P_s,
                                "half the white-light -> J difference in duration (Bailey et al. 1985)")
    d["chromatic_J_K_s"] = (0.5 * (v["stlmi_dur_K"] - v["stlmi_dur_J"]) * P_s, "half the J -> K difference")
    lam_um = math.log(1.25 / 0.55)          # ln(lambda_J / lambda_white), nominal 1.25 and 0.55 micron
    slope = 0.5 * (v["stlmi_dur_J"] - v["stlmi_dur_white"]) / lam_um
    d["chromatic_g_i_s"] = (slope * math.log(BAND_UM["i"] / BAND_UM["g"]) * P_s,
                            "the white->J edge shift scaled by ln(lambda_i/lambda_g)/ln(1.25/0.55): a log-linear "
                            "interpolation, ORDER OF MAGNITUDE ONLY; sign: redder band ends later")

    # --- cyclotron harmonics in our bands (PH.P6) -------------------------------------------------
    lam1 = v["cyc_lambda_100MG_um"] * 100.0 / v["stlmi_B_campbell"]
    d["cyc_fundamental_um"] = (lam1, "1.0710 um x (100 MG / 12.1 MG)")
    for b, lam in BAND_UM.items():
        d[f"harmonic_{b}"] = (lam1 / lam, f"lambda_c / lambda_{b}")

    # --- YZ Cnc -----------------------------------------------------------------------------------
    cyc_per_day = 1.0 / v["yz_vsx_period"]
    d["yz_day_drift_sh"] = (cyc_per_day * v["yz_sh_period_err"] / v["yz_vsx_period"],
                            "phase drift over 1 d under Shafter & Hessman's published sigma_P (cycles)")
    d["yz_day_drift_vp"] = ((1 / v["yz_vp_period"]) * v["yz_vp_period_err"] / v["yz_vp_period"],
                            "the same under van Paradijs et al.'s sigma_P")
    d["yz_day_shift_between_periods"] = (cyc_per_day * (v["yz_vp_period"] - v["yz_vsx_period"]) / v["yz_vsx_period"],
                                         "phase slip per day from folding on 0.0868 d instead of 0.086924 d")
    d["yz_vp_epoch_drift_2024"] = ((2460400.0 - v["yz_vp_epoch"]) / v["yz_vp_period"] * v["yz_vp_period_err"] / v["yz_vp_period"],
                                   "cycles of phase uncertainty in propagating the 1988 epoch to 2024: absolute phase is lost")
    d["yz_super_amp"] = (v["yz_kato_min"] - v["yz_kato_max"], "15.0 - 11.0 (Kato et al. 2002)")
    d["yz_normal_amp"] = (v["yz_hakala_quiescence"] - v["yz_hakala_outburst"], "15.2 - 12.2 (Hakala et al. 2004)")
    d["yz_super_minus_normal_vp"] = (v["yz_vp_normal_peak"] - v["yz_vp_super_peak"], "12.0 - 10.6 (van Paradijs et al. 1994)")
    d["sh_semi_peak_mmag"] = (500.0 * v["sh_amp_max_lowi"], "half of 0.25 mag, in mmag")
    d["yz_sh_semi_tess_mmag"] = (500.0 * v["yz_sh_amp_tess_early"], "half of 0.3 mag, in mmag")
    d["yz_hump_semi_vp_mmag"] = (500.0 * v["yz_vp_hump_full_mag"], "half of 0.5 mag, in mmag")
    d["yz_decline_semi_vp_mmag"] = (500.0 * v["yz_vp_decline_full_mag"], "half of 0.2 mag, in mmag")
    return d


# ==================================================================================================
# 4. EMISSION
# ==================================================================================================
def sci_tex(x: float, sig: int = 1) -> str:
    e = math.floor(math.log10(abs(x)))
    return f"{x / 10 ** e:.{sig}f} \\times 10^{{{e}}}"


#: macro name -> (formatted value builder, source comment).  Only what comparison_section.tex uses.
def macros(d: dict[str, tuple[float, str]]) -> list[tuple[str, str, str]]:
    v = {k: c.value for k, c in C.items()}
    g = lambda k: d[k][0]                                                             # noqa: E731
    lit = lambda k: f"literature: {C[k].key}, {C[k].where}"                            # noqa: E731
    der = lambda k: f"derived in literature_scales.py: {d[k][1]}"                      # noqa: E731
    return [
        ("LitStLmiCropperPeriodD", f"{v['stlmi_cropper_period']:.8f}", lit("stlmi_cropper_period")),
        ("LitStLmiCropperPeriodErrD", sci_tex(v["stlmi_cropper_period_err"], 0), lit("stlmi_cropper_period_err")),
        ("LitStLmiCropperEpoch", f"{v['stlmi_cropper_epoch']:.4f}", lit("stlmi_cropper_epoch")),
        ("LitStLmiBfield", f"{v['stlmi_B_campbell']:.1f}", lit("stlmi_B_campbell")),
        ("LitStLmiVsxMinusCropperSigma", f"{g('vsx_minus_cropper_sigma'):.1f}", der("vsx_minus_cropper_sigma")),
        ("LitStLmiRefitMinusCropperSigma", f"{abs(g('ours_minus_cropper_sigma')):.1f}", der("ours_minus_cropper_sigma")),
        ("LitStLmiCropperDriftCycles", f"{g('cropper_drift_cycles'):.3f}", der("cropper_drift_cycles")),
        ("LitStLmiCropperDriftMargin", f"{g('cropper_drift_margin'):.0f}", der("cropper_drift_margin")),
        ("LitStLmiVsxZeroOnCropper", f"{g('vsx_zero_on_cropper'):.2f}", der("vsx_zero_on_cropper")),
        ("LitStLmiVsxZeroOnCropperErr", f"{g('vsx_zero_on_cropper_err'):.2f}", der("vsx_zero_on_cropper_err")),
        ("LitStLmiEdgeOnCropper", f"{g('edge_on_cropper'):.2f}", der("edge_on_cropper")),
        ("LitStLmiEdgeOnConj", f"{g('edge_on_conj'):.2f}", der("edge_on_conj")),
        ("LitStLmiVsxZeroOnConjErr", f"{g('vsx_zero_on_conj_err'):.2f}", der("vsx_zero_on_conj_err")),
        ("LitStLmiCropperBaselineYr", f"{g('baseline_cropper_yr'):.0f}", der("baseline_cropper_yr")),
        ("LitStLmiFallRangeS", f"{g('fall_s_lo'):.0f}--{g('fall_s_hi'):.0f}", der("fall_s_lo") + " ... " + d["fall_s_hi"][1]),
        ("LitStLmiSeasonShiftS", f"{g('season_shift_s'):.0f}", der("season_shift_s")),
        ("LitStLmiSeasonShiftDeg", f"{g('season_shift_deg'):.0f}", der("season_shift_deg")),
        ("LitStLmiChromWhiteJS", f"{g('chromatic_white_J_s'):.0f}", der("chromatic_white_J_s")),
        ("LitStLmiChromJKS", f"{g('chromatic_J_K_s'):.0f}", der("chromatic_J_K_s")),
        ("LitStLmiChromGIS", f"{g('chromatic_g_i_s'):.0f}", der("chromatic_g_i_s")),
        ("LitStLmiHarmonicG", f"{g('harmonic_g'):.0f}", der("harmonic_g")),
        ("LitStLmiHarmonicR", f"{g('harmonic_r'):.0f}", der("harmonic_r")),
        ("LitStLmiHarmonicI", f"{g('harmonic_i'):.0f}", der("harmonic_i")),
        ("LitStLmiStateRecurrenceD", f"{v['stlmi_state_recurrence']:.0f} \\pm 90", lit("stlmi_state_recurrence")),
        ("LitStLmiStateDurationD", f"{v['stlmi_state_duration']:.1f} \\pm 7.6", lit("stlmi_state_duration")),
        ("LitVfifteenPdot", sci_tex(v["v1500_pdot"], 2), lit("v1500_pdot")),
        ("LitVfifteenOverBound", f"{g('v1500_over_bound'):.0f}", der("v1500_over_bound")),
        ("LitPdotSecularRange", f"{sci_tex(g('pdot_gr_lo'), 0)}$--${sci_tex(g('pdot_gr_hi'), 0)}", der("pdot_gr_lo")),
        ("LitLibrationPdot", sci_tex(g("libration_pdot"), 0), der("libration_pdot")),
        ("LitBoundOverLibration", f"{g('bound_over_libration'):.0f}", der("bound_over_libration")),
        ("LitDpLeoRateDegYr", f"{g('dpleo_rate_deg_yr'):.1f}", der("dpleo_rate_deg_yr")),
        ("LitStLmiRmsDeg", f"{g('rms_deg'):.1f}", der("rms_deg")),
        ("LitYzVpPeriodD", "0.086924", lit("yz_vp_period")),
        ("LitYzVpPeriodErrD", "0.000007", lit("yz_vp_period_err")),
        ("LitYzDayDriftSh", f"{g('yz_day_drift_sh'):.3f}", der("yz_day_drift_sh")),
        ("LitYzDayDriftVp", f"{g('yz_day_drift_vp'):.4f}", der("yz_day_drift_vp")),
        ("LitYzDaySlip", f"{g('yz_day_shift_between_periods'):.3f}", der("yz_day_shift_between_periods")),
        ("LitYzSuperAmp", f"{g('yz_super_amp'):.1f}", der("yz_super_amp")),
        ("LitYzNormalAmp", f"{g('yz_normal_amp'):.1f}", der("yz_normal_amp")),
        ("LitYzHumpSemiVpMmag", f"{g('yz_hump_semi_vp_mmag'):.0f}", der("yz_hump_semi_vp_mmag")),
        ("LitYzDeclineSemiVpMmag", f"{g('yz_decline_semi_vp_mmag'):.0f}", der("yz_decline_semi_vp_mmag")),
        ("LitShSemiPeakMmag", f"{g('sh_semi_peak_mmag'):.0f}", der("sh_semi_peak_mmag")),
        ("LitYzShSemiTessMmag", f"{g('yz_sh_semi_tess_mmag'):.0f}", der("yz_sh_semi_tess_mmag")),
        ("LitYzPshKato", "0.09031(5)", lit("yz_psh_2007")),
        ("LitYzPshTess", "0.09043(27)", lit("yz_psh_tess")),
        ("LitYzSupercycleD", "134 \\pm 19", lit("yz_pat_ts")),
        # -- literature values quoted verbatim in the prose -----------------------------------------
        ("LitStLmiFaintFraction", f"{v['stlmi_faint_fraction']:.2f}", lit("stlmi_faint_fraction")),
        ("LitStLmiFallPhases", f"{v['stlmi_fall_start_phase']:.2f}--{v['stlmi_fall_end_phase'] + 1:.2f}", lit("stlmi_fall_start_phase")),
        ("LitStLmiFallEndPhase", f"{v['stlmi_fall_end_phase']:.2f}", lit("stlmi_fall_end_phase")),
        ("LitStLmiFallWidthPeacock", f"{v['stlmi_fall_width_peacock']:.2f}", lit("stlmi_fall_width_peacock")),
        ("LitStLmiFallWidthRange", f"{v['stlmi_fall_width_peacock']:.2f}--{v['stlmi_fall_end_phase'] + 1 - v['stlmi_fall_start_phase']:.2f}",
         lit("stlmi_fall_width_peacock") + "; " + lit("stlmi_fall_start_phase")),
        ("LitStLmiDropPhaseRob", f"{v['stlmi_drop_rob']:.2f}", lit("stlmi_drop_rob")),
        ("LitStLmiBfieldFerrario", "12", lit("stlmi_B_ferrario")),
        ("LitStLmiCycFundamentalUm", f"{g('cyc_fundamental_um'):.0f}", der("cyc_fundamental_um")),
        ("LitStLmiSfstPeriod", "0.079087 \\pm 0.000001", lit("stlmi_sfst_period")),
        ("LitStLmiBwsPeriodD", f"{v['stlmi_bws_period']:.7f}", lit("stlmi_bws_period")),
        ("LitStLmiDurCropperEarly", f"{v['stlmi_dur_1982']:.2f}", lit("stlmi_dur_1982")),
        ("LitStLmiDurCropperLate", f"{v['stlmi_dur_1985']:.2f}", lit("stlmi_dur_1985")),
        ("LitStLmiDurBaileyEarly", f"{v['stlmi_dur_opt_1982jun']:.2f}", lit("stlmi_dur_opt_1982jun")),
        ("LitStLmiDurBaileyLate", f"{v['stlmi_dur_opt_1983']:.2f}", lit("stlmi_dur_opt_1983")),
        ("LitStLmiDurWhite", f"{v['stlmi_dur_white']:.2f}", lit("stlmi_dur_white")),
        ("LitStLmiDurJ", f"{v['stlmi_dur_J']:.2f}", lit("stlmi_dur_J")),
        ("LitStLmiDurH", f"{v['stlmi_dur_H']:.2f}", lit("stlmi_dur_H")),
        ("LitStLmiDurK", f"{v['stlmi_dur_K']:.2f}", lit("stlmi_dur_K")),
        ("LitStLmiExtentDeg", "15--18", lit("stlmi_extent_dense")),
        ("LitStLmiLowV", f"{v['stlmi_low_V']:.1f}", lit("stlmi_low_V")),
        ("LitStLmiExtremeLowV", f"{v['stlmi_extreme_low_V']:.1f}", lit("stlmi_extreme_low_V")),
        ("LitDpLeoAmpDeg", f"{v['dpleo_lib_amp_deg']:.0f}", lit("dpleo_lib_amp_deg")),
        ("LitDpLeoPeriodYr", f"{v['dpleo_lib_period_yr']:.0f}", lit("dpleo_lib_period_yr")),
        ("LitYzTnTess", f"{v['yz_tess_tn_lo']:.0f}--{v['yz_tess_tn_hi']:.0f}", lit("yz_tess_tn_lo")),
        ("LitYzTnPatterson", "10.3 \\pm 2.5", lit("yz_pat_tn")),
        ("LitYzKatoMax", f"{v['yz_kato_max']:.1f}", lit("yz_kato_max")),
        ("LitYzKatoMin", f"{v['yz_kato_min']:.1f}", lit("yz_kato_min")),
        ("LitShAmpMaxMag", f"{v['sh_amp_max_lowi']:.2f}", lit("sh_amp_max_lowi")),
        ("LitYzShAmpTessMag", f"{v['yz_sh_amp_tess_early']:.1f}", lit("yz_sh_amp_tess_early")),
        ("LitYzHumpFullVpMag", f"{v['yz_vp_hump_full_mag']:.1f}", lit("yz_vp_hump_full_mag")),
        ("LitYzDeclineFullVpMag", f"{v['yz_vp_decline_full_mag']:.1f}", lit("yz_vp_decline_full_mag")),
        ("LitYzHumpMaxPhase", f"{v['yz_vp_hump_max_phase']:.1f}", lit("yz_vp_hump_max_phase")),
        ("LitYzDeclineMaxPhase", f"{v['yz_vp_decline_max_phase']:.1f}", lit("yz_vp_decline_max_phase")),
        ("LitYzEpsPct", f"{v['yz_eps_2007_pct']:.1f}", lit("yz_eps_2007_pct")),
        ("LitYzShPeriod", "0.0868(2)", lit("yz_vsx_period")),
        ("LitYzVpPeriodParen", "0.086924(7)", lit("yz_vp_period")),
        ("LitYzFlickerPp", f"{v['yz_flicker_pp']:.2f}", lit("yz_flicker_pp")),
        ("LitStLmiBfieldErr", "0.5", lit("stlmi_B_campbell")),
        ("LitStLmiIrHarmonics", "4--7", "literature: campbell2008c, ST LMi modelling results (n = 4-7 at 2.25-1.30 um)"),
        # -- the ephemeris table of comparison_section.tex ------------------------------------------
        ("LitStLmiVsxEpoch", f"{v['stlmi_vsx_epoch']:.4f}", lit("stlmi_vsx_epoch")),
        ("LitVvPupEpoch", f"{v['vvpup_epoch']:.4f}", lit("vvpup_epoch")),
        ("LitVvPupPeriodD", f"{v['vvpup_period']:.10f}", lit("vvpup_period")),
        ("LitAnUmaVsxEpoch", f"{v['anuma_vsx_epoch']:.4f}", lit("anuma_vsx_epoch")),
        ("LitAnUmaOkPeriodD", "0.079752867(12)", lit("anuma_ok_period")),
        ("LitEuUmaVsxEpoch", f"{v['euuma_vsx_epoch']:.4f}", lit("euuma_vsx_epoch")),
        ("LitYzVpEpoch", f"{v['yz_vp_epoch']:.4f}", lit("yz_vp_epoch")),
        ("LitYzVpEpochDriftCycles", f"{g('yz_vp_epoch_drift_2024'):.0f}", der("yz_vp_epoch_drift_2024")),
    ]


def main() -> int:
    d = derive()
    rows = macros(d)
    names = [n for n, _, _ in rows]
    assert len(names) == len(set(names)), "duplicate macro name"
    assert all(re.fullmatch(r"[A-Za-z]+", n) for n in names), "TeX macro names are letters only"
    tex = ["%% literature_macros.tex -- GENERATED by committee/work/cv-literature/literature_scales.py.",
           "%% Do not edit.  External (literature) constants and arithmetic on them; none is a measurement",
           "%% of this programme.  Each line names its source, in the manner of numbers.tex."]
    tex += [f"\\newcommand{{\\{n}}}{{{val}}}  % [{src}]" for n, val, src in rows]
    (HERE / "literature_macros.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")
    # The manuscript inputs the same file from its own directory, so the
    # paper builds from manuscripts/CV_TimeSeries/ alone (CV-R12/R13).
    (NUMBERS_TEX.parent / "literature_macros.tex").write_text(
        "\n".join(tex) + "\n", encoding="utf-8")

    md = ["# Literature constants and predicted scales (CV paper)", "",
          "Emitted by `committee/work/cv-literature/literature_scales.py`. Do not edit by hand.", "",
          "## 1. Constants as read", "", "| name | value | source key | where | how read | note |", "|---|---|---|---|---|---|"]
    md += [f"| `{k}` | {c.value} | `{c.key}` | {c.where} | {c.how} | {c.note} |" for k, c in C.items()]
    md += ["", "## 2. Derived (arithmetic in `derive()`)", "", "| name | value | computed as |", "|---|---|---|"]
    md += [f"| `{k}` | {val:.6g} | {how} |" for k, (val, how) in d.items()]
    (HERE / "literature_scales.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"{len(C)} constants, {len(d)} derived scales, {len(rows)} macros -> literature_macros.tex, literature_scales.md")
    for k, (val, how) in d.items():
        print(f"  {k:32s} {val:14.6g}   {how[:90]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
