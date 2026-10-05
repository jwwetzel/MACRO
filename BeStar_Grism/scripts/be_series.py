#!/usr/bin/env python
"""be_series — from per-frame measurements to the paper's results: nightly EW
series (BE-S7), error floors from the standards (BE-S10), cross-calibration
across mechanical states (BE-S9), the relative continuum light curve (BE-S8),
the telluric-H2O regression test (BE-S6), the event detection and onset timing
(BE-S12), the slow-tier period search (BE-S11), the BeSS one-to-one (BE-S13)
and the V/R series (BE-VR-hold, if the delivered LSF admits it).

PRE-DECLARED RULES (each cites its ruling; none was tuned on the results)
-----------------------------------------------------------------------
Channel.  The science series is the hrg unit (hrg + HaGrism).  lrg (R ~ 400,
  LSF ~ 16 Å against a 90 Å EW window) is measured and compared night by night
  with hrg (the strategy's hrg->lrg self-test, §4 Step 7); it carries no claim.
Nightly value (§4 Step 7).  Median of the QC-passing frames of a night; a night
  with < 3 frames is kept but flagged (``n_lt3``) and never enters detection.
Errors (BE-S10, SYNTHESIS §4).  Intra-night error of a nightly median =
  1.2533 * sd(frames)/sqrt(n) (with the frame errors as a floor).  From the
  standards epoch (>= 2025-12-05) the per-night error adds in quadrature the
  night-to-night jitter s of the STANDARDS (eta Hya, tet Vir) per mechanical
  state: s^2 = var(standard residuals) - mean(intra^2), pooled over standards.
  Before the standards epoch the error is the intra-night LOWER BOUND, labelled.
  Every chi2_nu is reported per star and state with its dof (standing rule 1).
Detection rule (strategy §5.3, standards epoch only).  A star has an EVENT when
  on >= 2 consecutive observed nights (gap <= 5 d, same mechanical state) its
  nightly EW departs from its own median in that state by more than 3x the
  per-night error, with the same sign.
Onset timing (§4 Step 12).  Piecewise template (constant, then linear ramp
  from t0) fitted to the state's nightly EWs; t0 on a fine grid, 68% interval
  where chi2 <= chi2_min + 1 (after scaling chi2 to chi2_nu = 1 if > 1).
Period search (BE-S11).  The joint-fit periodogram of be_tscore, nuisance = an
  offset per mechanical state + airmass, focus, CCD temperature and the H2O
  band depth; only the stars admitted by be_injection (>= 10 standards-epoch
  hrg nights); significance = night-label bootstrap of the max power over the
  band (10 000 draws), FAP < 1%.  Each null carries the injected 90% amplitude
  in Å and, as the "predicted scale" (standing rule 2), the star's own BeSS EW
  scatter over the campaign.
Continuum (BE-S8, standing rule 6).  Relative continuum magnitudes at 6563 Å,
  standards epoch only: m = a_star + z_night + k_state * airmass (ensemble of
  every star observed that night, Honeycutt 1992); a star's continuum "varies"
  if its residual rms exceeds 3x the standards' residual rms.  No flux is
  quoted for season 1.
BeSS one-to-one (BE-S13).  BeSS spectra within +-2 d of an RLMT hrg night,
  degraded to the median hrg LSF and measured on the same windows.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402
import be_tscore as T  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_series.py"
OUT = C.NOTES / "results"
STANDARD_IDS = ("* eta Hya", "* tet Vir")
SEED = 20261006
N_BOOT = 10000
EVENT_SIGMA = 3.0
EVENT_GAP_D = 5.0
PAIR_D = 2.0


def mad_sd(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


# ===========================================================================
# nightly series + floors
# ===========================================================================
def build_nightly(con) -> pd.DataFrame:
    m = pd.read_sql("""SELECT * FROM be_frame_meas WHERE qc='pass'""", con)
    m, slopes = h2o_correct(m)
    m.to_sql("be_frame_meas_corr", con, if_exists="replace", index=False)
    slopes.to_sql("be_h2o_slopes", con, if_exists="replace", index=False)
    m["unit"] = m.grism
    rows = []
    for (mid, unit, night), g in m.groupby(["main_id", "unit", "night"]):
        n = len(g)
        # robust intra-night scatter (one bad frame must not set a night's error)
        sd = mad_sd(g.ew) if n >= 3 else (g.ew.std(ddof=1) if n > 1 else np.nan)
        intra = max(1.2533 * sd / np.sqrt(n) if n > 1 else np.nan,
                    float(np.sqrt(np.mean(g.ew_err ** 2)) / np.sqrt(n)))
        rows.append(dict(main_id=mid, label=g.label.iloc[0], role=g.role.iloc[0], filter=unit,
                         night=night, season=g.season.iloc[0], mech_state=g.mech_state.iloc[0],
                         std_epoch=int(night >= C.STANDARDS_EPOCH), n=n, n_lt3=int(n < 3),
                         t=float(g.bjd_tdb.mean()), ew=float(g.ew.median()), ew_sd=sd, err_intra=intra,
                         h2o_depth=float(g.h2o_depth.median()), airmass=float(g.airmass_calc.median()),
                         focpos=float(pd.to_numeric(g.focpos, errors="coerce").median()),
                         ccd_temp=float(pd.to_numeric(g.ccd_temp, errors="coerce").median()),
                         cont_rate=float(g.cont_rate.median()), lsf_fwhm_a=float(g.lsf_fwhm_a.median())))
    nt = pd.DataFrame(rows)
    # Water is atmospheric: a night with no valid band depth for this star
    # takes the median depth of the other stars observed that night through
    # the same grism, else the median of its mechanical state (flagged).
    nt["h2o_imputed"] = nt.h2o_depth.isna().astype(int)
    same_night = nt.groupby(["filter", "night"]).h2o_depth.transform("median")
    same_state = nt.groupby(["filter", "mech_state"]).h2o_depth.transform("median")
    nt["h2o_depth"] = nt.h2o_depth.fillna(same_night).fillna(same_state)
    return nt


def floors(nt: pd.DataFrame) -> pd.DataFrame:
    """Night-to-night jitter of the standards per (unit, state), std epoch."""
    s = nt[(nt.main_id.isin(STANDARD_IDS)) & (nt.std_epoch == 1) & (nt.n_lt3 == 0)].copy()
    out = []
    for (unit, state), g in s.groupby(["filter", "mech_state"]):
        res, intra2 = [], []
        for _, gg in g.groupby("main_id"):
            if len(gg) < 3:
                continue
            res += list(gg.ew - gg.ew.mean())
            intra2 += list(gg.err_intra ** 2)
        if len(res) < 4:
            continue
        nstar = g.main_id.nunique()
        sd_r = mad_sd(res)                     # robust night-to-night scatter
        s2 = max(sd_r ** 2 - np.median(intra2), 0.0)
        out.append(dict(filter=unit, mech_state=state, n_nights=len(res), n_std=nstar,
                        sd_resid=float(sd_r), sd_resid_classical=float(np.std(res, ddof=nstar)),
                        median_intra=float(np.sqrt(np.median(intra2))), s_night=float(np.sqrt(s2))))
    return pd.DataFrame(out)


def apply_errors(nt, fl):
    nt = nt.merge(fl[["filter", "mech_state", "s_night"]], on=["filter", "mech_state"], how="left")
    std = nt.std_epoch == 1
    nt["err"] = np.where(std & nt.s_night.notna(), np.sqrt(nt.err_intra ** 2 + nt.s_night.fillna(0) ** 2),
                         nt.err_intra)
    nt["err_label"] = np.where(std & nt.s_night.notna(), "standards", "intra-night lower bound")
    return nt


def chi2_table(nt):
    out = []
    for (mid, label, role, unit, state, std), g in nt[nt.n_lt3 == 0].groupby(
            ["main_id", "label", "role", "filter", "mech_state", "std_epoch"]):
        if len(g) < 3:
            continue
        w = 1 / g.err ** 2
        mu = np.sum(w * g.ew) / np.sum(w)
        chi = float(np.sum(w * (g.ew - mu) ** 2))
        dof = len(g) - 1
        out.append(dict(star=label, role=role, filter=unit, mech_state=state, std_epoch=std,
                        n_nights=len(g), mean_ew=float(mu), rms=float(g.ew.std(ddof=1)),
                        med_err=float(g.err.median()), chi2=chi, dof=dof, chi2_nu=chi / dof,
                        err_label=g.err_label.iloc[0]))
    return pd.DataFrame(out)


# ===========================================================================
# cross-calibration (S9)
# ===========================================================================
def crosscal(nt):
    out = []
    states = sorted(nt.mech_state.dropna().unique())
    for a, b in zip(states[:-1], states[1:]):
        for (mid, label, role, unit), g in nt[nt.n_lt3 == 0].groupby(["main_id", "label", "role", "filter"]):
            ga, gb = g[g.mech_state == a], g[g.mech_state == b]
            if len(ga) < 2 or len(gb) < 2 or role != "standard":
                continue
            d = gb.ew.mean() - ga.ew.mean()
            e = np.hypot(mad_sd(ga.ew) / np.sqrt(len(ga)), mad_sd(gb.ew) / np.sqrt(len(gb)))
            d = gb.ew.median() - ga.ew.median()
            out.append(dict(boundary=f"{a} | {b}", star=label, filter=unit, n_a=len(ga), n_b=len(gb),
                            offset=float(d), err=float(e)))
    df = pd.DataFrame(out)
    summ = []
    for a, b in zip(states[:-1], states[1:]):
        bd = f"{a} | {b}"
        for unit in ("hrg", "lrg"):
            g = df[(df.boundary == bd) & (df["filter"] == unit)] if len(df) else df
            if len(g) == 0:
                summ.append(dict(boundary=bd, filter=unit, n_stars=0, offset=np.nan, err=np.nan,
                                 chi2=np.nan, dof=0, verdict="UNBRIDGED (no standard on both sides)"))
                continue
            w = 1 / g.err ** 2
            mu = float(np.sum(w * g.offset) / np.sum(w))
            chi = float(np.sum(w * (g.offset - mu) ** 2))
            stat = float(1 / np.sqrt(np.sum(w)))
            spread = float(g.offset.std(ddof=1)) if len(g) > 1 else np.nan
            summ.append(dict(boundary=bd, filter=unit, n_stars=len(g), offset=mu,
                             err=stat, star_spread=spread,
                             err_total=float(np.hypot(stat, spread)) if len(g) > 1 else stat,
                             chi2=chi, dof=len(g) - 1,
                             verdict="bridged" if len(g) >= 2 else "bridged by one star (dof 0)"))
    return df, pd.DataFrame(summ)


# ===========================================================================
# continuum (S8)
# ===========================================================================
def continuum(nt):
    d = nt[(nt.std_epoch == 1) & (nt.cont_rate > 0) & (nt["filter"] == "hrg")].copy()
    d["m"] = -2.5 * np.log10(d.cont_rate)
    stars = sorted(d.main_id.unique())
    nights = sorted(d.night.unique())
    st = sorted(d.mech_state.unique())
    si = {s: i for i, s in enumerate(stars)}
    ni = {n: i for i, n in enumerate(nights)}
    ki = {s: i for i, s in enumerate(st)}
    A = np.zeros((len(d), len(stars) + len(nights) + len(st)))
    for r, row in enumerate(d.itertuples()):
        A[r, si[row.main_id]] = 1
        A[r, len(stars) + ni[row.night]] = 1
        A[r, len(stars) + len(nights) + ki[row.mech_state]] = row.airmass
    y = d.m.values
    w = np.ones(len(y))
    for _ in range(5):                       # robust reweighting
        sol, *_ = np.linalg.lstsq(A * w[:, None], y * w, rcond=None)
        r = y - A @ sol
        s = mad_sd(r)
        w = np.where(np.abs(r) > 4 * s, 0.0, 1.0)
    d["cont_resid_mag"] = y - A @ sol
    per = (d.groupby(["main_id", "label", "role"]).cont_resid_mag
           .agg(mad_sd).rename("rms_mag").reset_index())
    std_rms = float(per[per.main_id.isin(STANDARD_IDS)].rms_mag.median())
    per["n_nights"] = d.groupby("main_id").size().reindex(per.main_id).values
    per["std_rms_mag"] = std_rms
    per["varies"] = (per.rms_mag > 3 * std_rms).astype(int)
    return d[["main_id", "label", "night", "t", "cont_resid_mag"]], per


# ===========================================================================
# events (S12)
# ===========================================================================
def events(nt):
    ev, onsets = [], []
    d = nt[(nt.role == "science") & (nt.std_epoch == 1) & (nt.n_lt3 == 0) & (nt["filter"] == "hrg")]
    for (mid, label, state), g in d.groupby(["main_id", "label", "mech_state"]):
        g = g.sort_values("t")
        if len(g) < 3:
            continue
        med = g.ew.median()
        z = (g.ew.values - med) / g.err.values
        t = g.t.values
        hits = []
        for i in range(len(g) - 1):
            if (t[i + 1] - t[i] <= EVENT_GAP_D and abs(z[i]) > EVENT_SIGMA and abs(z[i + 1]) > EVENT_SIGMA
                    and np.sign(z[i]) == np.sign(z[i + 1])):
                hits.append(i)
        if not hits:
            continue
        i0 = hits[0]
        ev.append(dict(main_id=mid, star=label, mech_state=state, n_nights=len(g),
                       first_night=g.night.iloc[i0], max_abs_z=float(np.max(np.abs(z))),
                       delta_ew=float(g.ew.values[i0] - med), n_hit_pairs=len(hits)))
        # onset: constant then ramp from t0
        tg = np.linspace(t[0] - 3, t[-1], 2000)
        chi = []
        for t0 in tg:
            X = np.column_stack([np.ones_like(t), np.clip(t - t0, 0, None)])
            W = 1 / g.err.values ** 2
            cov = np.linalg.pinv(X.T @ (X * W[:, None]))
            b = cov @ (X.T @ (W * g.ew.values))
            chi.append(float(np.sum(W * (g.ew.values - X @ b) ** 2)))
        chi = np.array(chi)
        k = int(np.argmin(chi))
        X = np.column_stack([np.ones_like(t), np.clip(t - tg[k], 0, None)])
        W = 1 / g.err.values ** 2
        bk = np.linalg.pinv(X.T @ (X * W[:, None])) @ (X.T @ (W * g.ew.values))
        # EW is negative for emission: a negative ramp is growing emission
        ev[-1]["ramp_A_per_d"] = float(bk[1])
        ev[-1]["sign"] = "emission growing" if bk[1] < 0 else "emission declining"
        ev[-1]["ew_change_over_ramp"] = float(bk[1] * (t[-1] - tg[k]))
        scale = max(chi[k] / max(len(g) - 3, 1), 1.0)
        ok = tg[chi <= chi[k] + scale]
        onsets.append(dict(main_id=mid, star=label, mech_state=state, t0_bjd=float(tg[k]),
                           t0_lo=float(ok.min()), t0_hi=float(ok.max()),
                           t0_halfwidth_d=float((ok.max() - ok.min()) / 2), chi2_min=float(chi[k]),
                           dof=len(g) - 3, chi2_scale=scale))
    return pd.DataFrame(ev), pd.DataFrame(onsets)


# ===========================================================================
# period search (S11)
# ===========================================================================
def search(con, nt):
    adm = pd.read_sql("SELECT main_id, star, channel, admitted FROM be_injection_summary WHERE admitted=1", con)
    inj = pd.read_sql("SELECT * FROM be_injection", con)
    bess = pd.read_sql("SELECT main_id, ew_native FROM be_bess_ew WHERE status='ok'", con)
    rng = np.random.default_rng(SEED)
    out, pgs = [], []
    for r in adm.itertuples():
        # exactly the nights the injection used (be_injection.nightly_design)
        g = nt[(nt.main_id == r.main_id) & (nt["filter"] == "hrg") & (nt.std_epoch == 1)
               & (nt.n_lt3 == 0)].sort_values("t")
        g = g[np.isfinite(g.ew)]
        if len(g) < 10:
            out.append(dict(main_id=r.main_id, star=r.star, n_nights=len(g), note="< 10 usable nights after QC"))
            continue
        t = g.t.values - 2460000.0
        # same order as be_injection.regressors_of: the ruled H2O regressor first
        reg = {"h2o_depth": g.h2o_depth.values, "airmass": g.airmass.values,
               "focpos": g.focpos.values, "ccd_temp": g.ccd_temp.values}
        N, names = T.nuisance_design(g.mech_state.values, reg)
        Q = T.projector(N)
        freqs = T.freq_grid(t, 1.2, p_max=np.ptp(t) / 3)
        y = g.ew.values / g.err.values.mean()
        P, a, b = T.power(t, y, Q, freqs)
        k = int(np.nanargmax(P))
        resid = Q @ y
        nullmax = T.max_power_null_bootstrap(t, resid, Q, freqs, N_BOOT, rng)
        fap = float((np.sum(nullmax >= P[k]) + 1) / (N_BOOT + 1))
        ii = inj[(inj.main_id == r.main_id) & (inj.channel == "hrg")]
        a90 = []
        for per, gg in ii.groupby("period_d"):
            okk = gg[gg.frac >= 0.9].sort_values("amp_sigma")
            if 5 <= per < 20 and len(okk):
                a90.append(okk.amp_sigma.iloc[0])
        # Alias forensics (strategy §4 Step 11, Dawson & Fabrycky 2010): a
        # peak within 3/T of an integer frequency has its alias partner inside
        # the excluded band P > T/3 — it may be the daily alias of slow
        # variability.  Test: add a quadratic trend per mechanical state to
        # the nuisance model and recompute the global FAP of the highest peak.
        T_ = float(np.ptp(t))
        near_daily = int(abs(freqs[k] - round(freqs[k])) < 3.0 / T_ and round(freqs[k]) >= 1)
        reg2 = dict(reg)
        for st_ in sorted(set(g.mech_state)):
            msk = (g.mech_state.values == st_).astype(float)
            tt = (t - t[msk > 0].mean()) * msk
            reg2[f"trend1[{st_}]"] = tt
            reg2[f"trend2[{st_}]"] = tt ** 2
        N2, _ = T.nuisance_design(g.mech_state.values, reg2, max_frac=0.6)
        Q2 = T.projector(N2)
        P2, _, _ = T.power(t, y, Q2, freqs)
        k2 = int(np.nanargmax(P2))
        null2 = T.max_power_null_bootstrap(t, Q2 @ y, Q2, freqs, 2000, rng)
        fap2 = float((np.sum(null2 >= P2[k2]) + 1) / 2001)
        sig = float(np.sqrt(np.mean(g.err ** 2)))
        bs = bess[bess.main_id == r.main_id].ew_native
        out.append(dict(main_id=r.main_id, star=r.star, n_nights=len(g), nuisance=";".join(names),
                        best_period_d=float(1 / freqs[k]), power=float(P[k]),
                        amp_A=float(np.hypot(a[k], b[k]) * g.err.values.mean()), fap_global=fap,
                        near_daily_alias=near_daily, fap_with_trend=fap2,
                        best_period_with_trend_d=float(1 / freqs[k2]),
                        significant=int(fap < 0.01),
                        period_claim=int(fap < 0.01 and fap2 < 0.01 and not near_daily),
                        sigma_night_A=sig,
                        a90_5_20d_A=float(np.median(a90) * sig) if a90 else np.nan,
                        predicted_scale_A=float(bs.std(ddof=1)) if len(bs) > 2 else np.nan,
                        n_bess=len(bs)))
        pgs.append(pd.DataFrame(dict(main_id=r.main_id, freq=freqs, power=P)))
    return pd.DataFrame(out), (pd.concat(pgs) if pgs else pd.DataFrame())


# ===========================================================================
# BeSS one-to-one (S13)
# ===========================================================================
def bess_pairs(con, nt, fwhm):
    b = pd.read_sql("SELECT * FROM be_bess_ew WHERE status='ok'", con)
    col = f"ew_fwhm{fwhm:g}"
    if col not in b:
        col = "ew_native"
    b["bjd"] = b.mjd + 2400000.5
    b["ew_match"] = b[col]
    rows = []
    for r in b.itertuples():
        g = nt[(nt.main_id == r.main_id) & (nt["filter"] == "hrg") & (nt.n_lt3 == 0)]
        if not len(g):
            continue
        dt_ = np.abs(g.t.values - r.bjd)
        k = int(np.argmin(dt_))
        if dt_[k] > PAIR_D:
            continue
        rr = g.iloc[k]
        rows.append(dict(main_id=r.main_id, star=rr.label, role=rr.role, night=rr.night,
                         mech_state=rr.mech_state, bess_id=r.bess_id, observer=r.observer,
                         dt_d=float(g.t.values[k] - r.bjd), ew_rlmt=rr.ew_cal, ew_rlmt_raw=rr.ew, offset_applied=rr.offset_applied, err_rlmt=rr.err,
                         ew_bess=r.ew_match, ew_bess_native=r.ew_native, bess_rp=r.rp))
    p = pd.DataFrame(rows)
    if len(p):
        p["diff"] = p.ew_rlmt - p.ew_bess
    return p, col


# ===========================================================================
# telluric regressor test (S6)
# ===========================================================================
H2O_SLOPE_STAR = "* eta Hya"     # derives the correction
H2O_TEST_STAR = "* tet Vir"      # tests it, independently


def h2o_regress(g):
    """Within-night, 4xMAD-clipped regression of EW on the water band depth.
    Returns (slope, err, n, r, depth_range) or None."""
    g = g[np.isfinite(g.h2o_depth) & np.isfinite(g.ew)].copy()
    if len(g) < 10:
        return None
    g["ew"] = g.ew - g.groupby("night").ew.transform("median")
    g["h2o_depth"] = g.h2o_depth - g.groupby("night").h2o_depth.transform("median")
    r0 = g.ew - g.ew.median()
    g = g[np.abs(r0) <= 4 * mad_sd(r0)]
    x = g.h2o_depth - g.h2o_depth.mean()
    y = g.ew - g.ew.mean()
    if np.sum(x * x) <= 0:
        return None
    slope = float(np.sum(x * y) / np.sum(x * x))
    res = y - slope * x
    se = float(np.sqrt(np.sum(res ** 2) / (len(g) - 2) / np.sum(x * x)))
    return slope, se, len(g), float(np.corrcoef(x, y)[0, 1]), float(np.ptp(g.h2o_depth))


def h2o_correct(m):
    """BE-S6 telluric H2O regressor applied per frame: in each (grism, state)
    the slope measured on eta Hya is removed from every frame,
    EW -> EW - slope * (depth - state median depth).  theta Vir, which did
    not set the slope, is the acceptance test (h2o_test)."""
    m = m.copy()
    m["ew_raw"] = m.ew
    slopes = []
    for (unit, state), g in m.groupby(["grism", "mech_state"]):
        fit = h2o_regress(g[g.main_id == H2O_SLOPE_STAR])
        med = g.h2o_depth.median()
        if fit is None:
            slopes.append(dict(filter=unit, mech_state=state, slope=np.nan, err=np.nan, n=0,
                               applied=0, note="no eta Hya frames: uncorrected"))
            continue
        sel = (m.grism == unit) & (m.mech_state == state) & m.h2o_depth.notna()
        m.loc[sel, "ew"] = m.loc[sel, "ew"] - fit[0] * (m.loc[sel, "h2o_depth"] - med)
        slopes.append(dict(filter=unit, mech_state=state, slope=fit[0], err=fit[1], n=fit[2],
                           applied=1, note="slope from eta Hya"))
    return m, pd.DataFrame(slopes)


def h2o_test(con):
    m = pd.read_sql("""SELECT label, grism, mech_state, night, ew_raw, ew, h2o_depth FROM be_frame_meas_corr
                       WHERE qc='pass' AND main_id=?""", con, params=(H2O_TEST_STAR,))
    out = []
    for (unit, state), g0 in m.groupby(["grism", "mech_state"]):
        for which in ("ew_raw", "ew"):
            g = g0.assign(ew=g0[which])
            fit = h2o_regress(g)
            if fit is None:
                continue
            out.append(dict(filter=unit, mech_state=state, ew_column="before correction" if which == "ew_raw"
                            else "after eta Hya correction", n_frames=fit[2], slope_A_per_depth=fit[0],
                            slope_err=fit[1], slope_sigma=fit[0] / fit[1], r=fit[3], depth_range=fit[4],
                            ew_change_over_range=fit[0] * fit[4]))
    return pd.DataFrame(out)


def vr_table(con, lsf):
    """BE-VR-hold: nightly V/R of the hrg emission profiles, standards epoch,
    for stars whose profile is resolved double-peaked on >= half their frames.
    Printed beside the measured hrg LSF (the hold's acceptance: FWHM at Hα
    stated before any V/R value)."""
    m = pd.read_sql("""SELECT main_id, label, night, mech_state, vr, vr_status, peak_sep_a
                       FROM be_frame_meas WHERE qc='pass' AND grism='hrg' AND role='science'
                         AND night >= ?""", con, params=(C.STANDARDS_EPOCH,))
    if "vr_status" not in m or not len(m):
        return pd.DataFrame(), pd.DataFrame()
    frac = m.groupby(["main_id", "label"]).vr_status.apply(lambda x: float((x == "double").mean()))
    keep = frac[frac >= 0.5].reset_index()
    d = m[m.main_id.isin(keep.main_id) & (m.vr_status == "double")]
    n = (d.groupby(["main_id", "label", "night", "mech_state"])
         .agg(vr=("vr", "median"), vr_sd=("vr", mad_sd), sep_a=("peak_sep_a", "median"), n=("vr", "size"))
         .reset_index())
    fw = float(lsf[lsf.grism == "hrg"].lsf_fwhm_a.median())
    summ = (n.groupby(["main_id", "label"]).agg(nights=("vr", "size"), vr_median=("vr", "median"),
                                                vr_min=("vr", "min"), vr_max=("vr", "max"),
                                                peak_sep_A=("sep_a", "median")).reset_index())
    summ = summ.merge(frac.rename("frac_double").reset_index(), on=["main_id", "label"])
    summ["lsf_fwhm_A"] = fw
    summ["peak_sep_over_lsf"] = summ.peak_sep_A / fw
    return n, summ


def hrg_lrg(nt):
    h = nt[(nt["filter"] == "hrg") & (nt.n_lt3 == 0)]
    l_ = nt[(nt["filter"] == "lrg") & (nt.n_lt3 == 0)]
    j = h.merge(l_, on=["main_id", "night"], suffixes=("_h", "_l"))
    j["diff"] = j.ew_l - j.ew_h
    return (j.groupby(["label_h", "role_h"]).agg(n_nights=("diff", "size"), median_diff=("diff", "median"),
                                                 mad_sd=("diff", mad_sd), ew_hrg=("ew_h", "median"))
            .reset_index().rename(columns={"label_h": "star", "role_h": "role"}))


def lsf_table(con):
    m = pd.read_sql("""SELECT grism, mech_state, night, lsf_fwhm_a, o2b_fwhm_a, focpos FROM be_frame_meas
                       WHERE qc IN ('pass','focus_off','wide_psf')""", con)
    m["focpos"] = pd.to_numeric(m.focpos, errors="coerce")
    n = (m.groupby(["grism", "mech_state", "night"])
         .agg(lsf_fwhm_a=("lsf_fwhm_a", "median"), o2b_fwhm_a=("o2b_fwhm_a", "median"),
              focpos=("focpos", "median"), n=("lsf_fwhm_a", "size"))
         .reset_index())
    # The focuser position is re-set nightly (temperature-tracking autofocus),
    # so a state-wide "nominal" focuser value does not exist; frames off their
    # OWN night's focus are already removed by the QC gate.  A night is
    # off-nominal when its delivered LSF exceeds 1.5x its state's median.
    n["lsf_state_med"] = n.groupby(["grism", "mech_state"]).lsf_fwhm_a.transform("median")
    n["off_nominal"] = (n.lsf_fwhm_a > 1.5 * n.lsf_state_med).astype(int)
    n["R"] = 6563.0 / n.lsf_fwhm_a
    return n


# ===========================================================================
def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["nightly", "results"])
    a = ap.parse_args(argv)
    con = C.be_db()
    if a.cmd == "nightly":
        nt = build_nightly(con)
        fl = floors(nt)
        nt = apply_errors(nt, fl)
        # BE-S9: offset of each state relative to the ASI post-monsoon state,
        # from the standards that straddle the boundary; the standards' spread
        # is carried as the offset's systematic error.  States no standard
        # bridges keep ew_cal = ew and are flagged unbridged.
        _, ccs = crosscal(nt)
        nt["ew_cal"], nt["offset_applied"], nt["bridged"] = nt.ew, 0.0, 0
        ref = "S3 ASI post-monsoon (flipped)"
        nt.loc[nt.mech_state == ref, "bridged"] = 1
        for r in ccs.itertuples():
            a, b = r.boundary.split(" | ")
            if a == ref and np.isfinite(r.offset):
                sel = (nt.mech_state == b) & (nt["filter"] == r.filter)
                nt.loc[sel, "ew_cal"] = nt.loc[sel, "ew"] - r.offset
                nt.loc[sel, "offset_applied"] = r.offset
                nt.loc[sel, "bridged"] = 1
        nt.to_sql("be_nightly_ew", con, if_exists="replace", index=False)
        fl.to_sql("be_floors", con, if_exists="replace", index=False)
        con.commit()
        print(fl.to_string())
        return
    nt = pd.read_sql("SELECT * FROM be_nightly_ew", con)
    chi = chi2_table(nt)
    cc, ccs = crosscal(nt)
    cont, contp = continuum(nt)
    ev, on = events(nt)
    # Robustness of the detection rule (reported beside it, not replacing
    # it): the BeSS-constant null stars scatter more than the standards.  Their
    # median per-star excess over the per-night error, per (filter, state), is
    # added in quadrature and the rule re-run.  The median is used so a null
    # star that turns out variable cannot set the floor.
    nf = []
    for (unit, state), g in chi[(chi.role == "null") & (chi.std_epoch == 1)].groupby(["filter", "mech_state"]):
        ex = np.sqrt(np.clip(g.rms ** 2 - g.med_err ** 2, 0, None))
        nf.append(dict(filter=unit, mech_state=state, n_null=len(g), null_excess_median=float(np.median(ex)),
                       null_rms_values=";".join(f"{a}:{b:.3f}" for a, b in zip(g.star, g.rms))))
    nf = pd.DataFrame(nf)
    nt_c = nt.merge(nf[["filter", "mech_state", "null_excess_median"]], on=["filter", "mech_state"], how="left")
    nt_c["err"] = np.sqrt(nt_c.err ** 2 + nt_c.null_excess_median.fillna(0) ** 2)
    ev_c, _ = events(nt_c)
    ev["survives_null_floor"] = ev.main_id.isin(set(ev_c.main_id)).astype(int) if len(ev) else []
    sr, pg = search(con, nt)
    lsf = lsf_table(con)
    fw = float(lsf[lsf.grism == "hrg"].lsf_fwhm_a.median())
    pairs, col = bess_pairs(con, nt, round(fw, 1))
    h2o = h2o_test(con)
    vr, vrs = vr_table(con, lsf)
    hl = hrg_lrg(nt)
    for name, df in (("be_chi2", chi), ("be_crosscal", cc), ("be_crosscal_summary", ccs),
                     ("be_continuum", cont), ("be_continuum_star", contp), ("be_events", ev),
                     ("be_onsets", on), ("be_search", sr), ("be_periodograms", pg), ("be_lsf", lsf),
                     ("be_bess_pairs", pairs), ("be_h2o_test", h2o), ("be_hrg_lrg", hl),
                     ("be_vr_nightly", vr), ("be_vr", vrs), ("be_null_floor", nf)):
        df.to_sql(name, con, if_exists="replace", index=False)
    con.execute("INSERT OR REPLACE INTO be_meta VALUES ('periodograms_opened', 'yes (be_series results)')")
    con.execute("INSERT OR REPLACE INTO be_meta VALUES ('bess_pair_column', ?)", (col,))
    con.commit()
    print("events:", len(ev), "| significant peaks:", int(sr.get("significant", pd.Series()).sum()),
          "| period claims:", int(sr.get("period_claim", pd.Series()).sum()),
          "| BeSS pairs:", len(pairs), col)


if __name__ == "__main__":
    main()
