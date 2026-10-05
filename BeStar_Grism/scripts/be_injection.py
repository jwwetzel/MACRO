#!/usr/bin/env python
"""be_injection — BE-S-1b: injection–recovery for every sample star, on its real
nightly timestamps and THROUGH the detrending, published before any periodogram
of the real data is opened.

WHAT IS PRE-DECLARED HERE (before any EW series exists)
------------------------------------------------------
* The slow-tier series of a star is its nightly EW per grism channel, on
  STANDARDS-EPOCH nights only (>= 2025-12-05): a period is a variability claim,
  and the detection rule exists only from that night (SYNTHESIS §4, BE-S10).
* A channel is ADMITTED to the slow-tier search when it has >= 10 such nights
  (the novelty gate's own season threshold).  Channels with 5–9 nights are
  injected so the reader sees why they are not searched; below 5, nothing.
* The short tier (0.3–2 d on frames) needs >= 3 nights of > 2 h span (strategy
  §3.3).  ``be_short_tier`` (be_step0.py) shows no science star qualifies, so
  no short-tier injection or search exists for this paper; that is stated, not
  hidden.
* Search band: f from 3/T (P < baseline/3, strategy §4 Step 11) to 1.2 d^-1.
* Detrending = the joint nuisance model of ``be_tscore``: a free EW offset per
  mechanical state, plus recomputed airmass, focus position and CCD temperature
  (and the H2O band depth once BE-S6 measures it — this run is REPEATED with it
  before the real periodograms are opened; ``be_meta.injection_regressors``
  records which set the published contour used).
* Recovery = the global periodogram peak lies within 1% of f_inj AND its power
  exceeds the FAP = 1% threshold of the max-statistic over the whole band.
* Noise: white Gaussian with unit night-to-night scatter, so amplitudes are in
  units of sigma_night; BE-S10 converts them to Angstrom with the measured
  per-star floor.  Real residuals may be redder; the real search uses the
  night-label bootstrap, not this Gaussian threshold.

OUTPUTS (``bestar.sqlite``)
---------------------------
be_injection          per (star, channel, period, amplitude): trials, recovered
                      fraction, SIGNED matched-cell bias of amplitude and
                      frequency (standing rule 3: never a median of |bias|).
be_injection_summary  per (star, channel): nights, baseline, nuisance columns,
                      threshold, the 90% amplitude per period band, whether the
                      90% contour closes on the grid, admission.
notes/injection/*.md  the emitted tables; the contour figure is panel (b) of
                      paper Figure 6 and is drawn by be_figures.py.

USAGE
    /opt/miniconda3/envs/rlmt-checks/bin/python BeStar_Grism/scripts/be_injection.py
"""
from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402
import be_tscore as T  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_injection.py"
OUT = C.NOTES / "injection"

SEED = 20261005
MIN_NIGHTS_INJECT = 5
MIN_NIGHTS_SEARCH = 10
F_MAX = 1.2                    # d^-1, strategy §4 Step 11
N_PERIODS = 40
N_PHASES = 20
AMPS = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0)
N_NULL = 4000
FAP = 0.01
F_TOL = 0.01


def nightly_design(con) -> pd.DataFrame:
    """Nightly rows (star, channel, night) on standards-epoch nights, with the
    timestamp and nightly-median regressors that the real series will carry.

    Once ``be_nightly_ew`` exists (be_series.py nightly — EWs measured, NO
    periodogram opened) the rows are exactly the nights the search will use:
    hrg unit, QC-passed, >= 3 frames; otherwise (the first, pre-extraction
    run) every non-excluded frame's night."""
    if con.execute("SELECT name FROM sqlite_master WHERE name='be_nightly_ew'").fetchone():
        g = pd.read_sql("""SELECT main_id, label, filter, night, t, mech_state AS state, airmass,
                                  focpos, ccd_temp FROM be_nightly_ew
                           WHERE role='science' AND std_epoch=1 AND n_lt3=0 AND filter='hrg'
                             AND ew IS NOT NULL""", con)
        return g
    fr = pd.read_sql("""SELECT main_id, label, filter, night, mech_state, bjd_tdb,
                               airmass_calc, focpos, ccd_temp
                        FROM be_frames WHERE role='science' AND disposition!='exclude'
                          AND night >= ?""", con, params=(C.STANDARDS_EPOCH,))
    for c in ("focpos", "ccd_temp"):
        fr[c] = pd.to_numeric(fr[c], errors="coerce")
    g = (fr.groupby(["main_id", "label", "filter", "night"])
         .agg(t=("bjd_tdb", "mean"), state=("mech_state", "first"),
              airmass=("airmass_calc", "median"), focpos=("focpos", "median"),
              ccd_temp=("ccd_temp", "median"))
         .reset_index())
    return g


def regressors_of(g: pd.DataFrame, h2o: pd.Series | None = None) -> dict:
    # The ruled telluric H2O regressor (SYNTHESIS §4 BE-S6) goes FIRST, so the
    # nuisance-column cap (be_tscore.nuisance_design) never drops it.
    reg = {}
    if h2o is not None:
        reg["h2o_depth"] = h2o.values
    reg.update({"airmass": g.airmass.values, "focpos": g.focpos.values, "ccd_temp": g.ccd_temp.values})
    return reg


def run_series(g: pd.DataFrame, rng, h2o=None) -> tuple[list[dict], dict]:
    t = (g.t.values - 2460000.0)
    N, names = T.nuisance_design(g.state.values, regressors_of(g, h2o))
    Q = T.projector(N)
    base = float(np.ptp(t))
    freqs = T.freq_grid(t, F_MAX, p_max=base / 3.0)
    nullmax = T.max_power_null_gaussian(t, Q, freqs, N_NULL, rng)
    thr = float(np.quantile(nullmax, 1 - FAP))
    periods = np.geomspace(1 / F_MAX, base / 3.0, N_PERIODS)
    rows = []
    for P in periods:
        f_inj = 1.0 / P
        phases = rng.uniform(0, 2 * np.pi, N_PHASES)
        sig = np.sin(2 * np.pi * f_inj * (t - t.mean())[:, None] + phases[None, :])   # n x m
        noise = rng.standard_normal(sig.shape)
        for A in AMPS:
            Y = A * sig + noise
            Pw, a, b = T.power(t, Y, Q, freqs)
            k = np.nanargmax(Pw, axis=0)
            f_rec = freqs[k]
            pmax = Pw[k, np.arange(Pw.shape[1])]
            amp = np.hypot(a[k, np.arange(a.shape[1])], b[k, np.arange(b.shape[1])])
            ok = (np.abs(f_rec - f_inj) / f_inj < F_TOL) & (pmax > thr)
            rows.append(dict(period_d=P, amp_sigma=A, n_trials=N_PHASES, n_rec=int(ok.sum()),
                             frac=float(ok.mean()),
                             amp_bias_signed=float(np.mean((amp[ok] - A) / A)) if ok.any() else np.nan,
                             freq_bias_signed=float(np.mean((f_rec[ok] - f_inj) / f_inj)) if ok.any() else np.nan))
            if ok.mean() >= 0.999 and A >= 2:     # well past the 90% contour at this period
                break
    summ = dict(n_nights=len(t), baseline_d=base, p_max_d=base / 3.0, n_freq=len(freqs),
                nuisance=";".join(names), n_nuisance=N.shape[1], threshold_power=thr)
    return rows, summ


def a90(df: pd.DataFrame) -> pd.Series:
    """Smallest grid amplitude with >= 90% recovery at each period (NaN = open)."""
    out = {}
    for P, g in df.groupby("period_d"):
        ok = g[g.frac >= 0.9].sort_values("amp_sigma")
        out[P] = ok.amp_sigma.iloc[0] if len(ok) else np.nan
    return pd.Series(out)


def main() -> None:
    import argparse
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    con = C.be_db()
    rng = np.random.default_rng(SEED)
    nightly = nightly_design(con)
    h2o_avail = con.execute("SELECT name FROM sqlite_master WHERE name='be_nightly_ew'").fetchone()
    allrows, summ = [], []
    for (mid, label, ch), g in nightly.groupby(["main_id", "label", "filter"]):
        g = g.sort_values("t").reset_index(drop=True)
        if len(g) < MIN_NIGHTS_INJECT:
            summ.append(dict(main_id=mid, star=label, channel=ch, n_nights=len(g),
                             admitted=0, note=f"< {MIN_NIGHTS_INJECT} standards-epoch nights: no injection, no search"))
            continue
        h2o = None
        if h2o_avail:
            e = pd.read_sql("SELECT night, h2o_depth FROM be_nightly_ew WHERE main_id=? AND filter=?",
                            con, params=(mid, ch)).set_index("night").h2o_depth
            h2o = g.night.map(e)
            if h2o.isna().any():
                h2o = None
        rows, s = run_series(g, rng, h2o)
        for r in rows:
            r.update(main_id=mid, star=label, channel=ch)
        allrows += rows
        df = pd.DataFrame(rows)
        A = a90(df)
        P = A.index.values
        s.update(main_id=mid, star=label, channel=ch,
                 admitted=int(len(g) >= MIN_NIGHTS_SEARCH),
                 a90_short=float(np.nanmedian(A[P < 5])) if (P < 5).any() else np.nan,
                 a90_mid=float(np.nanmedian(A[(P >= 5) & (P < 20)])) if ((P >= 5) & (P < 20)).any() else np.nan,
                 a90_long=float(np.nanmedian(A[P >= 20])) if (P >= 20).any() else np.nan,
                 frac_periods_closed=float(np.isfinite(A).mean()),
                 amp_bias_signed_all=float(np.nanmean(df.amp_bias_signed)),
                 freq_bias_signed_all=float(np.nanmean(df.freq_bias_signed)),
                 note="admitted to slow-tier search" if len(g) >= MIN_NIGHTS_SEARCH
                 else f"{MIN_NIGHTS_INJECT}–{MIN_NIGHTS_SEARCH - 1} nights: injected, not searched")
        summ.append(s)
        print(label, ch, len(g), f"thr={s['threshold_power']:.3f}", f"A90 mid={s['a90_mid']}",
              f"closed={s['frac_periods_closed']:.2f}")
    inj = pd.DataFrame(allrows)
    sm = pd.DataFrame(summ)
    inj.to_sql("be_injection", con, if_exists="replace", index=False)
    sm.to_sql("be_injection_summary", con, if_exists="replace", index=False)
    con.execute("CREATE TABLE IF NOT EXISTS be_meta (key TEXT PRIMARY KEY, value TEXT)")
    regs = "airmass,focpos,ccd_temp" + (",h2o_depth" if h2o_avail else "")
    for k, v in (("injection_built_utc", dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")),
                 ("injection_regressors", regs), ("injection_seed", SEED),
                 ("periodograms_opened", "no")):
        con.execute("INSERT OR REPLACE INTO be_meta VALUES (?,?)", (k, str(v)))
    con.commit()

    meta = dict(con.execute("SELECT key, value FROM be_meta").fetchall())
    cols = ["star", "channel", "n_nights", "baseline_d", "p_max_d", "nuisance", "threshold_power",
            "a90_short", "a90_mid", "a90_long", "frac_periods_closed", "amp_bias_signed_all",
            "freq_bias_signed_all", "admitted", "note"]
    body = C.md_table(sm.reindex(columns=cols).sort_values(["admitted", "n_nights"], ascending=False), 3)
    body += (f"\nAmplitudes in units of the night-to-night scatter σ_night (white-noise injections; BE-S10 "
             f"converts to Å). `a90_*` = median over the period grid of the smallest amplitude recovered in "
             f"≥ 90% of {N_PHASES} phases, for P < 5 d / 5–20 d / ≥ 20 d (NaN = contour open in that band). "
             f"Recovery: global peak within {F_TOL:.0%} of f_inj and power above the FAP = {FAP:.0%} "
             f"max-statistic threshold ({N_NULL} null draws through the same nuisance model). Bias columns "
             f"are SIGNED means over recovered trials (standing rule 3). Regressors: {regs}. Search band "
             f"3/T – {F_MAX} d⁻¹. No science star has ≥ 3 nights of > 2 h span, so there is no short tier "
             f"(`notes/step0/short_tier.md`). Periodograms of the real data opened at build time: "
             f"{meta.get('periodograms_opened')}.\n\nBuilt {meta['injection_built_utc']}, seed {SEED}.\n")
    C.write_md(OUT / "summary.md", "Injection–recovery per sample star (BE-S-1b) — before any periodogram",
               body, SCRIPT)
    cell = inj.pivot_table(index=["star", "channel", "period_d"], columns="amp_sigma", values="frac")
    cell.round(2).reset_index().to_csv(OUT / "recovery_grid.csv", index=False)
    print("written", OUT)


if __name__ == "__main__":
    main()
