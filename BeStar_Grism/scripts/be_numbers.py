#!/usr/bin/env python
"""be_numbers — every number the Be-star manuscript quotes, emitted as LaTeX
macros from ``bestar.sqlite`` (house rule: no number is typed), plus the
Markdown results tables the ledger links.

Writes ``manuscripts/BeStar_Grism/numbers.tex`` and
``BeStar_Grism/notes/results/*.md``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_numbers.py"
OUT = C.NOTES / "results"
TEX = C.REPO / "manuscripts" / "BeStar_Grism" / "numbers.tex"


def q(con, sql, *a):
    return pd.read_sql(sql, con, params=a)


def tex_escape(x):
    return str(x).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")


def tex_table(path, caption, label, df, fmt):
    cols = list(df.columns)
    out = [f"% emitted by {SCRIPT}; do not edit", r"\begin{deluxetable*}{" + "l" * len(cols) + "}",
           r"\tablecaption{" + caption + r"\label{" + label + "}}",
           r"\tablehead{" + " & ".join(r"\colhead{" + tex_escape(c) + "}" for c in cols) + "}",
           r"\startdata"]
    rows = []
    for r in df.itertuples(index=False):
        rows.append(" & ".join((tex_escape(v) if isinstance(v, str) else fmt.get(c, "{}").format(v))
                               if v == v else "\\nodata" for c, v in zip(cols, r)))
    out.append(" \\\\\n".join(rows))
    out += [r"\enddata", r"\end{deluxetable*}"]
    path.write_text("\n".join(out).replace("Å", r"\AA{}").replace("Hα", r"H$\alpha$") + "\n",
                    encoding="utf-8")


def write_tex_tables(con):
    d = TEX.parent
    s = q(con, """SELECT s.label AS Star, s.sptype AS SpT, s.verified_seasons AS `BeSS-active`,
                  ROUND(s.bess_ew_med,1) AS `BeSS EW`, c.rlmt_nights AS `RLMT nights`, c.bess_nights AS `BeSS nights`,
                  (SELECT COUNT(*) FROM be_nightly_ew n WHERE n.main_id=s.main_id AND n.filter='hrg' AND n.n_lt3=0)
                    AS `hrg nights`,
                  (SELECT COUNT(*) FROM be_nightly_ew n WHERE n.main_id=s.main_id AND n.filter='hrg' AND n.n_lt3=0
                    AND n.std_epoch=1) AS `std-epoch`
                  FROM be_sample s JOIN be_cadence_star c USING (main_id) WHERE s.role='science'
                  ORDER BY s.label""")
    s["Star"] = s.Star.map(tex_escape)
    tex_table(d / "tab_sample.tex", "The 19 BeSS-verified active Be stars. BeSS EW: median over the "
              "verifying spectra (novelty check windows). Nights: RLMT grism nights (all grisms) and BeSS Hα "
              "nights over our season windows; hrg: nights with $\\ge3$ QC-passing hrg frames; std-epoch: of "
              "those, on or after 2025 December 5.", "tab:sample", s, {"BeSS EW": "{:.1f}"})
    a = q(con, """SELECT star AS Star, n_nights AS Nights, ROUND(best_period_d,2) AS `Highest peak (d)`,
                  fap_global AS FAP, fap_with_trend AS `FAP with trend`, ROUND(a90_5_20d_A,2) AS `a90 (Å)`,
                  ROUND(predicted_scale_A,2) AS `BeSS EW sd (Å)` FROM be_search ORDER BY star""")
    a["Star"] = a.Star.map(tex_escape)
    tex_table(d / "tab_a90.tex", "Slow-tier period search (standards epoch, hrg). FAP: global, night-label "
              "bootstrap of the maximum power over 3/T--1.2 d$^{-1}$; FAP with trend: the same with a quadratic "
              "trend per mechanical state in the nuisance model. a90: median over 5--20 d of the smallest "
              "sinusoid amplitude recovered in $\\ge90\\%$ of injections through the same detrending. "
              "BeSS EW sd: the star's own BeSS scatter over the campaign (the predicted scale).",
              "tab:a90", a, {"FAP": "{:.4f}", "FAP with trend": "{:.3f}"})
    e = q(con, """SELECT e.star AS Star, e.sign AS Change, e.first_night AS `First flagged`,
                  ROUND(e.ew_change_over_ramp,2) AS `Ramp $\\Delta$EW (Å)`, ROUND(e.max_abs_z,1) AS `max |z|`, e.survives_null_floor AS `Null floor`,
                  ROUND(o.t0_bjd-2460000,1) AS `Onset BJD$-$2460000`, ROUND(o.t0_halfwidth_d,1) AS `$\\pm$ (d)`
                  FROM be_events e LEFT JOIN be_onsets o USING (main_id, mech_state) ORDER BY e.star""")
    e["Star"] = e.Star.map(tex_escape)
    tex_table(d / "tab_events.tex", "Stars whose nightly EW departs from their own median by more than three "
              "times the per-night error on two consecutive nights after 2025 December 5. Onset: constant-then-"
              "ramp template; interval where $\\chi^2\\le\\chi^2_{\\rm min}+s$, with $s=\\chi^2_\\nu$ when "
              "above one.", "tab:events", e, {})


def main() -> None:
    import argparse
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    con = C.be_db()
    meta = dict(con.execute("SELECT key, value FROM be_meta").fetchall())
    N = {}
    nt = q(con, "SELECT * FROM be_nightly_ew")
    sci = nt[(nt.role == "science") & (nt["filter"] == "hrg")]
    N["NumStars"] = int(q(con, "SELECT COUNT(*) n FROM be_sample WHERE role='science'").n[0])
    N["NumStarsWithData"] = int(sci.main_id.nunique())
    N["NumNights"] = int(len(sci[["main_id", "night"]].drop_duplicates()))
    N["NumFramesPass"] = int(q(con, "SELECT COUNT(*) n FROM be_frame_meas WHERE qc='pass' AND role='science'").n[0])
    N["NumFramesReduced"] = int(q(con, "SELECT COUNT(*) n FROM be_frame_meas").n[0])
    N["NumCadenceAll"] = f"{float(meta['cadence_median_ratio_all']):.0f}"
    N["NumCadenceSparse"] = f"{float(meta['cadence_median_ratio_sparse']):.0f}"
    N["NumSparse"] = meta["cadence_n_sparse"]
    fl = q(con, "SELECT * FROM be_floors WHERE filter='hrg'")
    N["NumFloor"] = (f"{fl.sd_resid.min():.2f}--{fl.sd_resid.max():.2f}" if len(fl) > 1
                     else f"{fl.sd_resid.iloc[0]:.2f}")
    std = sci[sci.std_epoch == 1]
    N["NumStdEpochNights"] = int(len(std[std.n_lt3 == 0]))
    N["NumStdEpochStars"] = int(std[std.n_lt3 == 0].main_id.nunique())
    N["NumDetThreshold"] = f"{3 * std[std.n_lt3 == 0].err.median():.2f}"
    p = q(con, "SELECT * FROM be_bess_pairs")
    if len(p):
        d = p["diff"]
        N["NumBessPairs"] = len(p)
        N["NumBessOffset"] = f"{np.median(d):+.2f}"
        N["NumBessScatter"] = f"{1.4826 * np.median(np.abs(d - np.median(d))):.2f}"
    ev = q(con, "SELECT * FROM be_events")
    on = q(con, "SELECT * FROM be_onsets")
    N["NumEventStars"] = int(ev.main_id.nunique()) if len(ev) else 0
    N["NumEventRobust"] = int(ev.survives_null_floor.sum()) if len(ev) else 0
    nf = q(con, "SELECT * FROM be_null_floor WHERE filter='hrg'")
    N["NumNullExcess"] = (f"{nf.null_excess_median.min():.2f}--{nf.null_excess_median.max():.2f}"
                          if len(nf) > 1 else f"{nf.null_excess_median.iloc[0]:.2f}")
    if len(on):
        N["NumOnsetPrecision"] = f"{on.t0_halfwidth_d.min():.0f}--{on.t0_halfwidth_d.max():.0f}"
    sr = q(con, "SELECT * FROM be_search")
    ok = sr[sr.best_period_d.notna()] if "best_period_d" in sr else sr.iloc[0:0]
    N["NumSearched"] = int(len(ok))
    N["NumPeriods"] = int(ok.period_claim.sum()) if len(ok) else 0
    N["NumPeaksAliased"] = int(ok.significant.sum() - ok.period_claim.sum()) if len(ok) else 0
    a = ok.a90_5_20d_A.dropna()
    if len(a):
        N["NumAninetyRange"] = f"{a.min():.1f}--{a.max():.1f}"
    tess = q(con, "SELECT * FROM be_tess WHERE status='ok' AND our_nights_in_sector>0")
    alltess = q(con, "SELECT * FROM be_tess WHERE status='ok'")
    n_cov = 0
    for o in on.itertuples():
        tt = alltess[alltess.main_id == o.main_id]
        lo, hi = o.t0_lo - 2400000.5, o.t0_hi - 2400000.5
        n_cov += int(((tt.t_max_mjd >= lo) & (tt.t_min_mjd <= hi)).any())
    N["NumEventTess"] = n_cov
    std_tess = tess[(tess.t_max_mjd + 2400000.5 >= pd.Timestamp(C.STANDARDS_EPOCH).to_julian_date())
                    & (tess.role == "science")]
    N["NumTessStdStars"] = int(std_tess.main_id.nunique())
    N["NumTessStars"] = int(tess[tess.role == "science"].main_id.nunique()) if len(tess) else 0
    lsf = q(con, "SELECT * FROM be_lsf")
    for g in ("hrg", "lrg"):
        x = lsf[lsf.grism == g]
        if len(x):
            N[f"NumLsf{g.capitalize()}"] = f"{x.lsf_fwhm_a.median():.1f}"
            N[f"NumR{g.capitalize()}"] = f"{(6563 / x.lsf_fwhm_a.median()):.0f}"
    ob = lsf[lsf.grism == "hrg"].groupby("mech_state").o2b_fwhm_a.median().dropna()
    if len(ob):
        N["NumOtwoBHrg"] = f"{ob.min():.1f}--{ob.max():.1f}"
    N["NumOffNominal"] = int(lsf[lsf.grism == "hrg"].off_nominal.sum())
    cc = q(con, "SELECT * FROM be_crosscal_summary WHERE filter='hrg' AND n_stars>0")
    if len(cc):
        N["NumQhyOffset"] = f"{cc.offset.iloc[0]:.2f}"
        N["NumQhyOffsetErr"] = f"{cc.err_total.iloc[0]:.2f}"
    h = q(con, "SELECT * FROM be_h2o_test WHERE filter='hrg'")
    after = h[h.ew_column.str.startswith("after")]
    before = h[h.ew_column.str.startswith("before")]
    N["NumHtwoOBeforeSig"] = f"{before.slope_sigma.abs().max():.1f}"
    N["NumHtwoOAfterSig"] = f"{after.slope_sigma.abs().max():.1f}"
    vr = q(con, "SELECT * FROM be_vr")
    N["NumVrStars"] = len(vr)
    cs = q(con, "SELECT * FROM be_continuum_star WHERE role='science'")
    N["NumContVary"] = int(cs.varies.sum())
    qc = q(con, "SELECT qc, COUNT(*) n FROM be_frame_meas WHERE role='science' GROUP BY qc").set_index("qc").n
    for k, mac in (("saturated_in_window", "NumSatWindow"), ("identity_reject", "NumIdReject"),
                   ("no_trace", "NumNoTrace"), ("focus_off", "NumFocusOff")):
        N[mac] = int(qc.get(k, 0))
    N["NumFramesScience"] = int(qc.sum())
    pr = q(con, "SELECT mech_state, COUNT(*) n, AVG(diff) m FROM be_bess_pairs GROUP BY mech_state")
    N["NumBessPairsPre"] = int(pr[pr.mech_state.str.startswith("S2")].n.sum())
    hl = q(con, "SELECT * FROM be_hrg_lrg")
    N["NumHrgLrgMedian"] = f"{hl.median_diff.median():+.2f}"
    if len(vr):
        N["NumVrMin"] = f"{vr.vr_min.min():.2f}"
        N["NumVrMax"] = f"{vr.vr_max.max():.2f}"
        N["NumVrSepMin"] = f"{vr.peak_sep_over_lsf.min():.1f}"
        N["NumVrSepMax"] = f"{vr.peak_sep_over_lsf.max():.1f}"
    # BeSS one-to-one split at |EW| = STRONG_SPLIT_A (declared here)
    STRONG_SPLIT_A = 20.0
    if len(p):
        w_ = p[np.abs(p.ew_bess) < STRONG_SPLIT_A]["diff"]
        s_ = p[np.abs(p.ew_bess) >= STRONG_SPLIT_A]["diff"]
        mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
        N["NumBessSplit"] = f"{STRONG_SPLIT_A:.0f}"
        N["NumBessScatterWeak"] = f"{mad(w_):.2f}"
        N["NumBessOffsetWeak"] = f"{np.median(w_):+.2f}"
        N["NumBessPairsWeak"] = len(w_)
        N["NumBessScatterStrong"] = f"{mad(s_):.2f}" if len(s_) else "--"
    # BeSS observer-to-observer scatter on the null stars (paper windows)
    bn = q(con, "SELECT e.main_id, e.ew_native FROM be_bess_ew e WHERE e.status='ok' AND e.role='null'")
    sds = bn.groupby("main_id").ew_native.agg(lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
                                              if len(x) >= 4 else np.nan).dropna()
    if len(sds):
        N["NumBessObsScatter"] = f"{sds.min():.2f}--{sds.max():.2f}"
    # null-star offsets across the ASI->QHY boundary (standards epoch medians)
    nn = nt[(nt.role == "null") & (nt["filter"] == "hrg") & (nt.std_epoch == 1) & (nt.n_lt3 == 0)]
    offs = []
    for mid, g in nn.groupby("main_id"):
        a_ = g[g.mech_state.str.startswith("S3")].ew
        b_ = g[g.mech_state.str.startswith("S5")].ew
        if len(a_) >= 2 and len(b_) >= 2:
            offs.append(b_.median() - a_.median())
    if offs:
        N["NumNullOffsetMin"] = f"{min(offs):.2f}"
        N["NumNullOffsetMax"] = f"{max(offs):.2f}"
        N["NumNullOffsetN"] = len(offs)
    write_tex_tables(con)
    lines = [f"% emitted by {SCRIPT} from BeStar_Grism/products/bestar.sqlite; do not edit"]
    for k, v in N.items():
        lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")
    TEX.parent.mkdir(parents=True, exist_ok=True)
    TEX.write_text("\n".join(lines) + "\n", encoding="utf-8")
    C.write_md(OUT / "numbers.md", "Manuscript numbers (macros in numbers.tex)",
               C.md_table(pd.DataFrame(list(N.items()), columns=["macro", "value"])), SCRIPT)

    tabs = {
        "floors.md": ("Error floors from the standards (BE-S10)", "SELECT * FROM be_floors"),
        "chi2.md": ("χ²ν per star, channel, mechanical state and epoch, with dof (standing rule 1)",
                    "SELECT * FROM be_chi2 ORDER BY role DESC, star, filter, mech_state"),
        "crosscal.md": ("Cross-calibration across mechanical states (BE-S9)",
                        "SELECT * FROM be_crosscal_summary"),
        "crosscal_stars.md": ("Per-standard offsets at each boundary", "SELECT * FROM be_crosscal"),
        "continuum.md": ("Relative continuum (standards epoch, hrg; BE-S8)", "SELECT * FROM be_continuum_star"),
        "null_floor.md": ("Null-star excess scatter per state (robustness of the detection rule)",
                          "SELECT * FROM be_null_floor"),
        "events.md": ("Events under the detection rule (BE-S12)", "SELECT * FROM be_events"),
        "onsets.md": ("Onset epochs (BE-S12)", "SELECT * FROM be_onsets"),
        "search.md": ("Slow-tier period search (BE-S11)", "SELECT * FROM be_search"),
        "lsf.md": ("Delivered LSF per night (BE-S5)", "SELECT * FROM be_lsf ORDER BY grism, night"),
        "h2o.md": ("θ Vir EW vs 7200 Å H₂O band depth (BE-S6 acceptance)", "SELECT * FROM be_h2o_test"),
        "hrg_lrg.md": ("hrg→lrg same-night self-test (strategy §4 Step 7)", "SELECT * FROM be_hrg_lrg"),
        "bess_pairs_summary.md": ("BeSS one-to-one per star (BE-S13)",
                                  "SELECT star, role, COUNT(*) n_pairs, AVG(diff) mean_diff, "
                                  "MIN(diff) min_diff, MAX(diff) max_diff FROM be_bess_pairs GROUP BY star"),
        "vr.md": ("V/R of resolved double-peaked hrg profiles (BE-VR-hold), with the measured LSF",
                  "SELECT * FROM be_vr"),
        "lsf_summary.md": ("Delivered resolution per (grism, state): G-5 method and O₂-B edge (BE-S5)",
                           "SELECT grism, mech_state, COUNT(*) nights, AVG(lsf_fwhm_a) lsf_xd_mean_A, "
                           "AVG(o2b_fwhm_a) o2b_edge_fwhm_mean_A, SUM(off_nominal) off_nominal_nights "
                           "FROM be_lsf GROUP BY 1,2"),
        "peaks.md": ("Trace-peak distribution as a fraction of the measured saturation cap — the saturation "
                     "statement (BE-S1)",
                     "SELECT role, grism, mech_state, COUNT(*) frames, "
                     "SUM(peak_frac_cap < 0.5) below_half_cap, SUM(peak_frac_cap >= 0.5 AND peak_frac_cap < 0.9) "
                     "half_to_0p9, SUM(peak_frac_cap >= 0.9) at_or_above_0p9, SUM(n_sat_win > 0) sat_in_window "
                     "FROM be_frame_meas WHERE peak_frac_cap IS NOT NULL GROUP BY 1,2,3"),
        "qc.md": ("Frame QC log (BE-S1)", "SELECT role, grism, qc, COUNT(*) frames FROM be_frame_meas "
                                         "GROUP BY 1,2,3 ORDER BY 1,2,4 DESC"),
        "o2b.md": ("O₂-B separation check per (grism, library epoch) (BE-S4)",
                   "SELECT grism, sol_epoch, mech_state, COUNT(o2b_frac_dev) n, AVG(o2b_frac_dev) mean_frac_dev, "
                   "AVG(o2b_frac_dev*o2b_frac_dev) ms FROM be_frame_meas WHERE qc='pass' GROUP BY 1,2,3"),
    }
    for fn, (title, sql) in tabs.items():
        try:
            df = q(con, sql)
        except Exception as exc:              # a table not built yet is stated, not hidden
            df = pd.DataFrame({"note": [f"not built: {exc}"]})
        if fn == "o2b.md" and len(df) and "ms" in df:
            df["rms_frac_dev"] = np.sqrt(df.ms - df.mean_frac_dev ** 2)
            df = df.drop(columns="ms")
        C.write_md(OUT / fn, title, C.md_table(df, 3), SCRIPT)
    print(N)


if __name__ == "__main__":
    main()
