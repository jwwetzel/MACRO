"""macro_tcrb.paper — numbers.tex for the T CrB manuscript and the
eruption-contingency skeleton.

Everything the manuscripts state as a number is emitted here as a LaTeX
macro, from ``products/tcrb/tcrb.sqlite`` (and read-only reads of the
manifest and the grism library).  A number whose source table does not
exist yet is emitted as ``\\NumMissing`` — the PDF shows a bold flag, never
a typed stand-in.  Figures: ``macro_tcrb.figures``.
"""

from __future__ import annotations

import math
import sqlite3
from pathlib import Path
from typing import Callable, Optional

from macro_tcrb import db

MS = db.REPO / "manuscripts" / "TCrB_Monitoring"


def fmt(x, nd: int = 2, sign: bool = False) -> str:
    """LaTeX-safe number: thin-space thousands, fixed decimals."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return r"\NumMissing"
    if isinstance(x, int) or (isinstance(x, float) and nd == 0):
        s = f"{int(round(x)):,}".replace(",", r"\,")
        return ("+" + s) if sign and x > 0 else s
    s = f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"
    return s.replace("-", "$-$") if s.startswith("-") else s


class Macros:
    """Accumulates \\newcommand lines with their provenance."""

    def __init__(self, prefix: str = "TC"):
        self.prefix = prefix
        self.lines: list[str] = []
        self.missing: list[str] = []

    def add(self, name: str, fn: Callable[[], object], src: str,
            nd: int = 2, sign: bool = False) -> None:
        try:
            v = fn()
        except (sqlite3.Error, IndexError, KeyError, TypeError,
                ValueError, ZeroDivisionError, RuntimeError):
            v = None
        if isinstance(v, str):
            txt = v
        else:
            txt = fmt(v, nd, sign)
        if txt == r"\NumMissing":
            self.missing.append(name)
        self.lines.append(f"\\newcommand{{\\{self.prefix}{name}}}{{{txt}}}"
                          f"  % [{src}]")

    def text(self, header: str) -> str:
        return "\n".join([
            "%% numbers.tex -- GENERATED FILE.  DO NOT EDIT.",
            f"%% {header}",
            r"\providecommand{\NumMissing}{\textbf{[NUMBER MISSING]}}",
            ""] + self.lines) + "\n"


def q1(con, sql, *a):
    return db.q1(con, sql, *a)


# ===========================================================================
# numbers.tex
# ===========================================================================
def emit_numbers() -> Path:
    import json
    import numpy as np
    con = db.connect()
    M = Macros("TC")
    t = lambda sql, *a: (lambda: db.q1(con, sql, *a))       # noqa: E731
    nov = json.loads((db.NOVELTY / "novelty_verdict.json").read_text())

    # ---- novelty (TCRB-N1 product) --------------------------------------
    M.add("NovARASDates", lambda: nov["facts"]["aras_dates_ha"],
          "novelty_verdict.json", 0)
    M.add("NovARASBestObs", lambda: nov["verdict"]["best_single_observer_dates"],
          "novelty_verdict.json", 0)
    M.add("NovARASObservers", lambda: nov["facts"]["aras_observers_ha"],
          "novelty_verdict.json", 0)
    M.add("NovWindowDays", lambda: nov["facts"]["window_days"],
          "novelty_verdict.json", 0)
    M.add("NovLabelled", lambda: nov["facts"]["rlmt_frames_labelled"],
          "novelty_verdict.json", 0)
    M.add("NovDispersed", lambda: nov["facts"]["rlmt_frames_dispersed"],
          "novelty_verdict.json", 0)

    # ---- Phase 0 ----------------------------------------------------------
    M.add("NgrismFrames", t("SELECT COUNT(*) FROM tcrb_temp_split"),
          "tcrb_temp_split", 0)
    M.add("Ncold", t("SELECT n_frames FROM tcrb_temp_summary WHERE "
                     "temp_group='cold'"), "tcrb_temp_summary", 0)
    M.add("Nwarm", t("SELECT n_frames FROM tcrb_temp_summary WHERE "
                     "temp_group='warm'"), "tcrb_temp_summary", 0)
    M.add("LeakColdPct", lambda: 100 * db.q1(con, "SELECT leak_none_frac_med "
          "FROM tcrb_temp_summary WHERE temp_group='cold'"),
          "tcrb_temp_summary", 2)
    M.add("LeakWarmPct", lambda: 100 * db.q1(con, "SELECT leak_none_frac_med "
          "FROM tcrb_temp_summary WHERE temp_group='warm'"),
          "tcrb_temp_summary", 2)
    M.add("ShutterBound", t("SELECT value FROM tcrb_shutter_summary WHERE "
                            "quantity='delta_bound_s'"),
          "tcrb_shutter_summary", 2)
    M.add("ShutterMedian", t("SELECT value FROM tcrb_shutter_summary WHERE "
                             "quantity='delta_median_s'"),
          "tcrb_shutter_summary", 3, sign=True)
    M.add("ShutterNladders", t("SELECT value FROM tcrb_shutter_summary WHERE "
                               "quantity='n_ladders_constraining'"),
          "tcrb_shutter_summary", 0)
    M.add("SkewTrace", lambda: 1000 * db.q1(con, "SELECT value FROM "
          "tcrb_shutter_summary WHERE quantity='skew_across_trace_bound_s'"),
          "tcrb_shutter_summary", 0)
    M.add("Gain", lambda: db.q1(db.manifest_ro(), "SELECT value FROM "
          "detector_params WHERE era_group='ASI Mode0 2x2' AND "
          "quantity='gain_e_per_adu'"), "detector_params", 3)
    M.add("GainErr", lambda: db.q1(db.manifest_ro(), "SELECT uncertainty FROM "
          "detector_params WHERE era_group='ASI Mode0 2x2' AND "
          "quantity='gain_e_per_adu'"), "detector_params", 3)
    M.add("ReadNoise", lambda: db.q1(db.manifest_ro(), "SELECT value FROM "
          "detector_params WHERE era_group='ASI Mode0 2x2' AND "
          "quantity='read_noise_e'"), "detector_params", 2)

    # ---- Phase B / C1 -------------------------------------------------------
    M.add("NBframes", t("SELECT COUNT(*) FROM tcrb_bphot_frames"),
          "tcrb_bphot_frames", 0)
    M.add("NBnights", t("SELECT COUNT(DISTINCT night) FROM tcrb_bphot_frames"),
          "tcrb_bphot_frames", 0)
    M.add("BPrecision", t("SELECT per_frame_rms_mag FROM tcrb_bprecision"),
          "tcrb_bprecision", 3)
    M.add("BFloor", t("SELECT season_floor_mag FROM tcrb_bprecision"),
          "tcrb_bprecision", 3)
    M.add("BColourTerm", t("SELECT AVG(k_bv) FROM tcrb_bphot_frames"),
          "tcrb_bphot_frames", 2, sign=True)
    M.add("BCrossMode", t("SELECT AVG(crossmode_mag) FROM tcrb_bphot_frames"),
          "tcrb_bphot_frames", 3)
    M.add("FlickRecMin", t("SELECT MIN(quoted_limit_rms) FROM tcrb_flicker"),
          "tcrb_flicker", 2)
    M.add("FlickRecMax", t("SELECT MAX(quoted_limit_rms) FROM tcrb_flicker"),
          "tcrb_flicker", 2)
    M.add("FlickPubAmp", t("SELECT MAX(published_amp_B) FROM tcrb_flicker"),
          "tcrb_flicker", 2)
    _a_track(con, M)
    return _finish(con, M)


def _a_track(con, M) -> None:
    import numpy as np
    from astropy.time import Time
    t = lambda sql, *a: (lambda: db.q1(con, sql, *a))       # noqa: E731
    g = db.grism_ro()
    hl = lambda q: (lambda: db.q1(con, "SELECT value FROM tcrb_headline "
                                       "WHERE quantity = ?", q))  # noqa
    hs = lambda q, c: (lambda: db.q1(con, f"SELECT {c} FROM tcrb_headline "
                                          "WHERE quantity = ?", q))  # noqa
    # ---- gate, dates --------------------------------------------------
    M.add("Nframes", t("SELECT COUNT(*) FROM tcrb_ew_frames"),
          "tcrb_ew_frames", 0)
    M.add("Ngate", t("SELECT COUNT(*) FROM tcrb_ew_frames WHERE status="
                     "'ok'"), "tcrb_ew_frames", 0)

    def ut_dates():
        jds = [r[0] for r in db.q(con, "SELECT jd FROM tcrb_ew_nightly")]
        return sorted({Time(j, format="jd").iso[:10] for j in jds})
    M.add("Ndates", lambda: len(ut_dates()), "tcrb_ew_nightly", 0)
    M.add("FirstDate", lambda: Time(ut_dates()[0]).strftime("%Y %B %-d"),
          "tcrb_ew_nightly")
    M.add("LastDate", lambda: Time(ut_dates()[-1]).strftime("%Y %B %-d"),
          "tcrb_ew_nightly")
    for gr, G in (("hrg", "Hrg"), ("lrg", "Lrg")):
        M.add(f"Nights{G}", t("SELECT COUNT(*) FROM tcrb_ew_nightly WHERE "
                              "grism=?", gr), "tcrb_ew_nightly", 0)
        M.add(f"Disp{G}", lambda gr=gr: abs(db.q1(g, "SELECT disp_a_per_px "
              "FROM g_dispersion WHERE grism=? AND mech_epoch='ASI-pre'",
              gr)), "g_dispersion", 3)
        M.add(f"LSF{G}", t("SELECT fwhm_a FROM tcrb_lsf_adopted WHERE "
                           "grism=?", gr), "tcrb_lsf_adopted", 1)
        M.add(f"R{G}", lambda gr=gr: 6562.8 / db.q1(con, "SELECT fwhm_a FROM "
              "tcrb_lsf_adopted WHERE grism=?", gr), "tcrb_lsf_adopted", 0)
        M.add(f"HalfWidth{G}", t("SELECT half_width_a FROM tcrb_lsf_adopted"
                                 " WHERE grism=?", gr), "tcrb_lsf_adopted",
              0)
        M.add(f"WinMed{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "median_abs_frac FROM tcrb_window_survival WHERE grism=?", gr),
              "tcrb_window_survival", 1)
        M.add(f"WinP{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "p90_abs_frac FROM tcrb_window_survival WHERE grism=?", gr),
              "tcrb_window_survival", 1)
        M.add(f"ShiftM{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "AVG(ew_m10/ew-1) FROM tcrb_ew_frames WHERE grism=? AND "
              "status='ok'", gr), "tcrb_ew_frames", 1, sign=True)
        M.add(f"ShiftP{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "AVG(ew_p10/ew-1) FROM tcrb_ew_frames WHERE grism=? AND "
              "status='ok'", gr), "tcrb_ew_frames", 1, sign=True)
        M.add(f"NAras{G}", t("SELECT n_dates FROM tcrb_xval_summary WHERE "
                             "grism=?", gr), "tcrb_xval_summary", 0)
        M.add(f"ArasMean{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "mean_frac FROM tcrb_xval_summary WHERE grism=?", gr),
              "tcrb_xval_summary", 1, sign=True)
        M.add(f"ArasMeanStat{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "stat FROM tcrb_headline WHERE quantity=?",
              f"ARAS_mean_frac_{gr}"), "tcrb_headline", 1)
        M.add(f"ArasRms{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "rms_frac FROM tcrb_xval_summary WHERE grism=?", gr),
              "tcrb_xval_summary", 1)
        M.add(f"Floor{G}", t("SELECT floor_ew FROM tcrb_floor WHERE grism=?",
                             gr), "tcrb_floor", 2)
        M.add(f"FrameFloor{G}", t("SELECT linefree_frame_floor FROM "
                                  "tcrb_floor WHERE grism=?", gr),
              "tcrb_floor", 2)
        M.add(f"NightFloor{G}", t("SELECT night_floor FROM tcrb_floor WHERE "
                                  "grism=?", gr), "tcrb_floor", 2)
        M.add(f"FloorChi{G}", t("SELECT chi2_all FROM tcrb_floor WHERE "
                                "grism=?", gr), "tcrb_floor", 1)
        M.add(f"FloorDof{G}", t("SELECT dof_all FROM tcrb_floor WHERE "
                                "grism=?", gr), "tcrb_floor", 0)
        M.add(f"FloorChiB{G}", t("SELECT chi2_halfB FROM tcrb_floor WHERE "
                                 "grism=?", gr), "tcrb_floor", 1)
        M.add(f"FloorDofB{G}", t("SELECT dof_halfB FROM tcrb_floor WHERE "
                                 "grism=?", gr), "tcrb_floor", 0)
        M.add(f"MethodDiff{G}", t("SELECT method_diff FROM tcrb_floor WHERE "
                                  "grism=?", gr), "tcrb_floor", 2)
        M.add(f"Smear{G}", t("SELECT smear FROM tcrb_floor WHERE grism=?",
                             gr), "tcrb_floor", 2)
        M.add(f"OffFocus{G}", t("SELECT n_off_focus FROM tcrb_floor WHERE "
                                "grism=?", gr), "tcrb_floor", 0)
        M.add(f"EWMean{G}", hl(f"EW_season_mean_{gr}"), "tcrb_headline", 1)
        M.add(f"EWMeanStat{G}", hs(f"EW_season_mean_{gr}", "stat"),
              "tcrb_headline", 1)
        M.add(f"EWMeanSys{G}", hs(f"EW_season_mean_{gr}", "sys"),
              "tcrb_headline", 1)
        M.add(f"EWScatter{G}", hl(f"EW_intrinsic_scatter_{gr}"),
              "tcrb_headline", 1)
        M.add(f"EWScatterStat{G}", hs(f"EW_intrinsic_scatter_{gr}", "stat"),
              "tcrb_headline", 1)
        M.add(f"EWMin{G}", t("SELECT MIN(ew) FROM tcrb_ew_nightly WHERE "
                             "grism=?", gr), "tcrb_ew_nightly", 1)
        M.add(f"EWMax{G}", t("SELECT MAX(ew) FROM tcrb_ew_nightly WHERE "
                             "grism=?", gr), "tcrb_ew_nightly", 1)
        M.add(f"EWErr{G}", t("SELECT AVG(ew_err) FROM tcrb_ew_nightly WHERE "
                             "grism=?", gr), "tcrb_ew_nightly", 1)
        M.add(f"ChiConst{G}", t("SELECT chi2 FROM tcrb_detect WHERE grism=?",
                                gr), "tcrb_detect", 0)
        M.add(f"DofConst{G}", t("SELECT dof FROM tcrb_detect WHERE grism=?",
                                gr), "tcrb_detect", 0)
        M.add(f"PConst{G}", lambda gr=gr: "{:.1e}".format(db.q1(
            con, "SELECT p_const FROM tcrb_detect WHERE grism=?", gr)
        ).replace("e-0", r"\times10^{-").replace("e-", r"\times10^{-")
            + "}", "tcrb_detect")
        M.add(f"StepRec{G}", t("SELECT min_step_recovered FROM tcrb_detect "
                               "WHERE grism=?", gr), "tcrb_detect", 1)
        M.add(f"Eligible{G}", lambda gr=gr: 100 * db.q1(con, "SELECT "
              "eligible_fraction FROM tcrb_detect WHERE grism=?", gr),
              "tcrb_detect", 0)
        M.add(f"Nevents{G}", t("SELECT n_events FROM tcrb_detect WHERE "
                               "grism=?", gr), "tcrb_detect", 0)
        M.add(f"FluxChi{G}", lambda gr=gr: _chi_flux(con, gr)[0],
              "tcrb_flux", 0)
        M.add(f"FluxDof{G}", lambda gr=gr: _chi_flux(con, gr)[1],
              "tcrb_flux", 0)
        M.add(f"ContSame{G}", t("SELECT COUNT(*) FROM tcrb_flux WHERE "
                                "grism=? AND cont_source LIKE '%same%'", gr),
              "tcrb_flux", 0)
        M.add(f"ContInterp{G}", t("SELECT COUNT(*) FROM tcrb_flux WHERE "
                                  "grism=? AND cont_source LIKE '%interp%'",
                                  gr), "tcrb_flux", 0)
    M.add("StepThreshold", t("SELECT MAX(threshold_sigma) FROM tcrb_detect"),
          "tcrb_detect", 0)
    M.add("NevaluatedAras", t("SELECT COUNT(*) FROM tcrb_aras_ew WHERE "
                              "date_obs BETWEEN '2025-02-21' AND "
                              "'2025-06-25'"), "tcrb_aras_ew", 0)
    M.add("NArasDates", lambda: len({r[0] for r in db.q(
        con, "SELECT ut_date FROM tcrb_xval")}), "tcrb_xval", 0)
    M.add("PhaseMin", t("SELECT MIN(phase) FROM tcrb_flux"), "tcrb_flux", 2)
    M.add("PhaseMax", t("SELECT MAX(phase) FROM tcrb_flux"), "tcrb_flux", 2)
    M.add("RcBright", t("SELECT MIN(cont_mag) FROM tcrb_flux"), "tcrb_flux",
          2)
    M.add("RcFaint", t("SELECT MAX(cont_mag) FROM tcrb_flux"), "tcrb_flux",
          2)
    # ---- theta CrB -------------------------------------------------------
    M.add("Ntet", t("SELECT COUNT(*) FROM tcrb_tet_frames"),
          "tcrb_tet_frames", 0)
    M.add("NtetTraced", t("SELECT COUNT(*) FROM tcrb_tet_frames WHERE "
                          "status='ok'"), "tcrb_tet_frames", 0)
    M.add("TetHa", t("SELECT ha_ew_median FROM tcrb_a1 WHERE grism='hrg' "
                     "AND grism_epoch='ASI-pre' AND exptime=4.86"),
          "tcrb_a1", 2)
    M.add("TetHaSd", t("SELECT ha_nightly_robust_sd FROM tcrb_a1 WHERE "
                       "grism='hrg' AND grism_epoch='ASI-pre' AND "
                       "exptime=4.86"), "tcrb_a1", 2)
    M.add("TetHaChiShort", t("SELECT chi2 FROM tcrb_a1 WHERE grism='hrg' "
                             "AND grism_epoch='ASI-pre' AND exptime=2.43"),
          "tcrb_a1", 0)
    M.add("TetHaDofShort", t("SELECT dof FROM tcrb_a1 WHERE grism='hrg' "
                             "AND grism_epoch='ASI-pre' AND exptime=2.43"),
          "tcrb_a1", 0)
    # ---- profiles ----------------------------------------------------------
    ps_ = lambda q: (lambda: db.q1(con, "SELECT value FROM "  # noqa
                                        "tcrb_profiles_summary WHERE "
                                        "quantity=?", q))
    pu_ = lambda q: (lambda: db.q1(con, "SELECT uncertainty FROM "  # noqa
                                        "tcrb_profiles_summary WHERE "
                                        "quantity=?", q))
    M.add("LSFkms", ps_("lsf_fwhm_kms"), "tcrb_profiles_summary", 0)
    M.add("FWHMkms", hl("Halpha_FWHM_kms"), "tcrb_headline", 0)
    M.add("FWHMkmsStat", hs("Halpha_FWHM_kms", "stat"), "tcrb_headline", 0)
    M.add("FWHMkmsSys", hs("Halpha_FWHM_kms", "sys"), "tcrb_headline", 0)
    M.add("WingBlue", lambda: -db.q1(con, "SELECT value FROM "
          "tcrb_profiles_summary WHERE quantity='hwzi10_blue_kms_median'"),
          "tcrb_profiles_summary", 0)
    M.add("WingRed", ps_("hwzi10_red_kms_median"), "tcrb_profiles_summary",
          0)
    M.add("FWHMTrend", ps_("trend_fwhm_kms_per_100d"),
          "tcrb_profiles_summary", 0, sign=True)
    M.add("FWHMTrendErr", pu_("trend_fwhm_kms_per_100d"),
          "tcrb_profiles_summary", 0)
    M.add("NprofFrames", ps_("n_frames"), "tcrb_profiles_summary", 0)
    # ---- A0b ---------------------------------------------------------------
    M.add("AzeroBextract", t("SELECT COUNT(*) FROM tcrb_a0b WHERE decision"
                             "='extract'"), "tcrb_a0b", 0)
    M.add("AzeroBtotal", t("SELECT COUNT(*) FROM tcrb_a0b"), "tcrb_a0b", 0)
    M.add("AzeroBEWmay", t("SELECT AVG(ew) FROM tcrb_a0b WHERE decision="
                           "'extract' AND night LIKE '2023-05%'"),
          "tcrb_a0b", 1)
    M.add("AzeroBEWmar", t("SELECT AVG(ew) FROM tcrb_a0b WHERE decision="
                           "'extract' AND night LIKE '2024-03%'"),
          "tcrb_a0b", 1)
    M.add("AzeroBEWmarErr", t("SELECT AVG(ew_err) FROM tcrb_a0b WHERE "
                              "decision='extract' AND night LIKE '2024-03%'"),
          "tcrb_a0b", 1)
    # ---- B anchors (C4) ----------------------------------------------------
    for i, night in enumerate(r[0] for r in db.q(
            con, "SELECT night FROM tcrb_b_anchors ORDER BY night")):
        k = "ABC"[i]
        M.add(f"BAnchor{k}", hl(f"B_{night}"), "tcrb_headline", 3)
        M.add(f"BAnchorStat{k}", hs(f"B_{night}", "stat"), "tcrb_headline",
              3)
        M.add(f"BAnchorSys{k}", hs(f"B_{night}", "sys"), "tcrb_headline", 3)
    M.add("RedNoiseBeta", hl("red_noise_beta_B"), "tcrb_headline", 2)
    man = db.manifest_ro()
    M.add("BinGainRatio", lambda: db.q1(man, "SELECT value FROM "
          "detector_params WHERE era_group='ASI Mode0 2x2' AND "
          "quantity='binning_gain_ratio'"), "detector_params", 1)
    M.add("LinCapPct", lambda: 100 * db.q1(man, "SELECT value FROM "
          "detector_params WHERE era_group='Mode0' AND "
          "quantity='linearity_cap_fraction'"), "detector_params", 0)
    M.add("FalseCalls", t("SELECT MAX(expected_false) FROM tcrb_detect"),
          "tcrb_detect", 4)
    M.add("FlickNmin", t("SELECT MIN(n_frames) FROM tcrb_flicker"),
          "tcrb_flicker", 0)
    M.add("FlickNmax", t("SELECT MAX(n_frames) FROM tcrb_flicker"),
          "tcrb_flicker", 0)
    M.add("FlickSpanMax", t("SELECT MAX(span_min) FROM tcrb_flicker"),
          "tcrb_flicker", 0)
    # ---- orbit and the AAVSO B curve across the RLMT window ----------------
    import csv
    eph = [r for r in csv.DictReader(open(db.NOVELTY / "external" /
                                          "ephemerides.csv"))
           if r["ephemeris"].startswith("Munari")][0]
    M.add("Period", lambda: float(eph["period_d"]), "ephemerides.csv", 1)
    M.add("OrbitFrac", lambda: db.q1(con, "SELECT MAX(phase)-MIN(phase) "
          "FROM tcrb_flux"), "tcrb_flux", 2)
    M.add("EllCycles", lambda: 2 * db.q1(con, "SELECT MAX(phase)-MIN(phase)"
          " FROM tcrb_flux"), "tcrb_flux", 2)

    def bcurve():
        j0, j1 = db.q(con, "SELECT MIN(jd), MAX(jd) FROM tcrb_ew_nightly")[0]
        rows = db.q(con, """SELECT CAST(jd + 0.5 AS INT) d, mag FROM
            tcrb_aavso WHERE band = 'B' AND own = 0 AND jd BETWEEN ? AND ?
            """, j0 - 0.5, j1 + 0.5)
        by: dict = {}
        for d, m in rows:
            by.setdefault(d, []).append(m)
        med = {d: float(np.median(v)) for d, v in by.items() if len(v) >= 3}
        ks = sorted(med)
        first = float(np.median([med[k] for k in ks[:5]]))
        last = float(np.median([med[k] for k in ks[-5:]]))
        return first, last, min(med.values()), max(med.values()), len(ks)
    M.add("BcurveFirst", lambda: bcurve()[0], "tcrb_aavso", 2)
    M.add("BcurveLast", lambda: bcurve()[1], "tcrb_aavso", 2)
    M.add("BcurveBright", lambda: bcurve()[2], "tcrb_aavso", 2)
    M.add("BcurveFaint", lambda: bcurve()[3], "tcrb_aavso", 2)
    M.add("BcurveNights", lambda: bcurve()[4], "tcrb_aavso", 0)
    M.add("NaavsoRows", t("SELECT COUNT(*) FROM tcrb_aavso"), "tcrb_aavso", 0)
    M.add("NarasSpectra", t("SELECT COUNT(*) FROM tcrb_ext_fetch WHERE "
                            "source='aras' AND ok=1"), "tcrb_ext_fetch", 0)


def _chi_flux(con, gr):
    import numpy as np
    a = np.array(db.q(con, "SELECT line_flux, line_flux_err FROM tcrb_flux"
                           " WHERE grism=? AND line_flux IS NOT NULL", gr),
                 float)
    w = 1 / a[:, 1] ** 2
    m = np.sum(w * a[:, 0]) / np.sum(w)
    return float(np.sum(((a[:, 0] - m) / a[:, 1]) ** 2)), int(len(a) - 1)


def _finish(con, M) -> Path:
    MS.mkdir(parents=True, exist_ok=True)
    out = MS / "numbers.tex"
    out.write_text(M.text(f"built {db.utc_now()} (git {db.git_commit()}) "
                          "from products/tcrb/tcrb.sqlite, the manifest and "
                          "the grism library"))
    return out


# ===========================================================================
# Tables (emitted, never typed)
# ===========================================================================
def _num(x, nd=3, sign=False):
    return fmt(x, nd, sign)


def emit_tables() -> list[Path]:
    con = db.connect()
    out = []
    # ---- Table: filter codes ---------------------------------------------
    rows = db.q(con, """SELECT code, n_frames, n_dispersed, n_direct,
        dzmag_vs_R_same_night, lam_eff_A, lam_eff_sd_A, best_matching_band,
        mapping FROM tcrb_filter_map WHERE code IN
        ('6','W','L','G','B','R','I','O','H','1') """)
    order = ["6", "W", "B", "G", "L", "R", "I", "O", "H", "1"]
    adopted = {"6": "grism", "W": "clear, then grism (2024 Feb)",
               "B": "Johnson $B$", "G": "Sloan $g$", "L": "luminance",
               "R": "Sloan $r$", "I": "Sloan $i$", "O": "[\\ion{O}{3}]",
               "H": "H$\\alpha$", "1": "red narrow band ([\\ion{S}{2}])"}
    band_tex = {"b_jkc": "$B$", "g_sdss": "$g$", "v_jkc": "$V$",
                "r_sdss": "$r$", "r_jkc": "$R$", "i_sdss": "$i$",
                "refcat2_g": "$g$"}
    r_by = {r[0]: r for r in rows}
    L = [r"\begin{deluxetable*}{lrrrrlll}",
         r"\tablecaption{The single-character filter codes of the 2023--2024"
         r" camera, mapped from the pixels.\label{tab:filters}}",
         r"\tablehead{\colhead{Code} & \colhead{Frames} & "
         r"\colhead{Dispersed} & \colhead{$\Delta$ZMAG vs R} & "
         r"\colhead{$\lambda_{\rm eff}$ (\AA)} & \colhead{Nearest band} & "
         r"\colhead{Adopted}}",
         r"\startdata"]
    for c in order:
        r = r_by.get(c)
        if r is None:
            continue
        lam = ("" if r[5] is None else
               f"{_num(r[5], 0)}$\\pm${_num(r[6], 0)}")
        L.append(f"{c} & {_num(r[1], 0)} & {_num(r[2], 0)} & "
                 f"{'' if r[4] is None else _num(r[4], 2, True)} & {lam} & "
                 f"{band_tex.get(r[7], '')} & {adopted[c]} \\\\")
    L += [r"\enddata",
          r"\tablecomments{Frames and dispersed counts: all canonical frames"
          r" of the camera (2023 Feb--2024 Mar), dispersion judged per frame"
          r" from the pixels. $\Delta$ZMAG: same-night mean difference to R "
          r"on direct frames. $\lambda_{\rm eff}$: zero crossing of the "
          r"colour term against nine Gaia-synthetic bands, on T~CrB-field "
          r"frames (mean $\pm$ s.d. over frames).}",
          r"\end{deluxetable*}"]
    p = MS / "tab_filters.tex"
    p.write_text("\n".join(L) + "\n")
    out.append(p)
    # ---- Table: B anchors -------------------------------------------------
    rows = db.q(con, "SELECT * FROM tcrb_b_anchors ORDER BY night")
    L = [r"\begin{deluxetable}{lrrrrr}",
         r"\tablecaption{Archival $B$ anchors of T~CrB.\label{tab:anchors}}",
         r"\tablehead{\colhead{Night} & \colhead{$N$} & \colhead{$B$} & "
         r"\colhead{$\sigma_{\rm frame}$} & \colhead{$\sigma_{\rm col}$} & "
         r"\colhead{AAVSO $B$}}", r"\startdata"]
    for r in rows:
        L.append(f"{r[0]} & {r[1]} & {_num(r[2])} & {_num(r[4])} & "
                 f"{_num(r[5])} & {'' if r[7] is None else _num(r[7])} \\\\")
    L += [r"\enddata", r"\tablecomments{$B$: transformed to Johnson $B$ "
          r"through a Gaia-synthetic ensemble with the AAVSO $B-V$ of the "
          r"night. $\sigma_{\rm frame}$: check-star rms per frame; "
          r"$\sigma_{\rm col}$: colour-extrapolation error. AAVSO: median "
          r"within 0.1\,d of the RLMT frames.}", r"\end{deluxetable}"]
    p = MS / "tab_anchors.tex"
    p.write_text("\n".join(L) + "\n")
    out.append(p)
    # ---- Table: flickering limits -----------------------------------------
    rows = db.q(con, "SELECT * FROM tcrb_flicker ORDER BY night")
    L = [r"\begin{deluxetable}{lrrrrr}",
         r"\tablecaption{What the archival $B$ snippets could bound."
         r"\label{tab:flicker}}",
         r"\tablehead{\colhead{Night} & \colhead{$N$} & \colhead{Span} & "
         r"\colhead{rms$_{\rm obs}$} & \colhead{rms$_{\rm check}$} & "
         r"\colhead{rms$_{90}$}\\ \colhead{} & \colhead{} & \colhead{(min)} &"
         r" \colhead{(mag)} & \colhead{(mag)} & \colhead{(mag)}}",
         r"\startdata"]
    for r in rows:
        L.append(f"{r[0]} & {r[2]} & {_num(r[3], 1)} & {_num(r[6])} & "
                 f"{_num(r[8])} & {_num(r[14], 2)} \\\\")
    L += [r"\enddata", r"\tablecomments{rms$_{90}$: the injected red-noise "
          r"($\beta=2$) rms recovered above the check-star 95\% threshold "
          r"in 90\% of injections — the quoted limit. The published $B$ "
          r"amplitude is 0.07\,mag \citep{Maslennikova2023}.}",
          r"\end{deluxetable}"]
    p = MS / "tab_flicker.tex"
    p.write_text("\n".join(L) + "\n")
    out.append(p)
    return out
