#!/usr/bin/env python
"""be_figures — the six-figure set of the Be-star paper (BE-figures), drawn
only from ``bestar.sqlite`` and only under the cuts Seat 6 approved
(committee/reviews/2026-10-05-seat6-bestar-abstract.md):

  1  RLMT against BeSS Hα nights per star, mechanical states as ticks
  2  the standards' nightly EW and the per-state error floor
  3  RLMT vs resolution-matched BeSS EW, one-to-one
  4  nightly EW curves of the 19 stars (+ BeSS points, events marked)
  5  event zoom with onset fit and TESS where contemporaneous — CUT if no
     standards-epoch event survives
  6  periodograms with injection–recovery contours — becomes a TABLE of a90
     limits if no period is significant

House style: ``macro_core.plotstyle`` (white ground, Okabe–Ito, marker as a
second channel).  Output: ``manuscripts/BeStar_Grism/figures/`` (PDF + PNG)
and ``BeStar_Grism/figures/`` (PNG copies for the repo).
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

sys.path.insert(0, str(C.REPO / "pipeline"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from macro_core import plotstyle as ps  # noqa: E402

MS = C.REPO / "manuscripts" / "BeStar_Grism" / "figures"
REPO_FIG = C.PROJECT / "figures"
STATE_STYLE = {
    "S1 Andor": (ps.FAINT, "v"), "S2 ASI pre-monsoon": (ps.WARN, "s"),
    "S3 ASI post-monsoon (flipped)": (ps.ACCENT, "o"), "S4 QHY night 1": (ps.OTHER, "D"),
    "S5 QHY re-seated": (ps.GOOD, "^")}
JD0 = 2460000.0


def save(fig, name):
    MS.mkdir(parents=True, exist_ok=True)
    REPO_FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(MS / f"{name}.pdf")
    fig.savefig(MS / f"{name}.png", dpi=ps.PNG_DPI)
    shutil.copy(MS / f"{name}.png", REPO_FIG / f"{name}.png")
    plt.close(fig)


def std_epoch_bjd():
    return pd.Timestamp(C.STANDARDS_EPOCH).to_julian_date() - JD0


def fig1(con):
    st = pd.read_sql("SELECT * FROM be_cadence_star", con)
    fr = pd.read_sql("""SELECT DISTINCT main_id, label, night, mech_state FROM be_frames
                        WHERE role='science' AND disposition!='exclude'""", con)
    nv = C.novelty_ro()
    b = pd.read_sql("SELECT main_id, substr(date,1,10) d FROM nv_bess WHERE covers_ha=1 "
                    "AND date BETWEEN '2024-11-01' AND '2026-07-31'", nv)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 4.6), gridspec_kw={"width_ratios": [1, 1.6]})
    for sp, (lab, mk) in ((1, ("BeSS-sparse (median gap > 14 d)", "o")), (0, ("BeSS dense", "s"))):
        g = st[st["sparse"] == sp]
        a1.scatter(g.bess_nights.clip(lower=0.7), g.rlmt_nights, marker=mk, s=22,
                   color=ps.ACCENT if sp else ps.MUTED, label=lab, zorder=3)
    x = np.array([0.7, 60])
    a1.plot(x, x, color=ps.RULE, lw=0.8, ls="--", label="equal sampling")
    med = float(np.median(st.rlmt_nights / st.bess_nights.clip(lower=1)))
    a1.plot(x, med * x, color=ps.WARN, lw=0.8, ls=":", label=f"median ratio {med:.1f}")
    a1.set_xscale("log"), a1.set_yscale("log")
    a1.set_xlim(0.6, 60), a1.set_ylim(8, 70)
    a1.set_xlabel("BeSS Hα nights in our season windows (0 plotted at 0.7)")
    a1.set_ylabel("RLMT grism nights")
    a1.legend(fontsize=6.5, loc="lower right")
    order = st.sort_values("rlmt_nights").star.tolist()
    for i, lab in enumerate(order):
        g = fr[fr.label == lab]
        for state, gg in g.groupby("mech_state"):
            col, _ = STATE_STYLE.get(state, (ps.INK, "o"))
            t = pd.to_datetime(gg.night)
            a2.vlines(t, i + 0.05, i + 0.45, color=col, lw=0.9)
        mid = st.loc[st.star == lab, "main_id"].iloc[0]
        tb = pd.to_datetime(b[b.main_id == mid].d)
        a2.vlines(tb, i - 0.45, i - 0.05, color=ps.INK, lw=0.6)
    a2.set_yticks(range(len(order)))
    a2.set_yticklabels(order, fontsize=6)
    a2.axvline(pd.Timestamp(C.STANDARDS_EPOCH), color=ps.BAD, lw=0.8, ls="--")
    for state, (col, _) in STATE_STYLE.items():
        if state in set(fr.mech_state):
            a2.plot([], [], color=col, lw=1.5, label=f"RLMT, {state}")
    a2.plot([], [], color=ps.INK, lw=0.8, label="BeSS Hα spectrum")
    a2.plot([], [], color=ps.BAD, ls="--", lw=0.8, label="standards epoch begins")
    a2.set_ylim(-1, len(order) + 2.6)       # head-room for the legend
    a2.legend(fontsize=5.5, loc="upper left", ncol=3)
    a2.set_xlabel("Date (UT)")
    fig.autofmt_xdate()
    fig.tight_layout()
    save(fig, "fig1_cadence")


def fig2(con):
    nt = pd.read_sql("""SELECT * FROM be_nightly_ew WHERE filter='hrg' AND n_lt3=0
                        AND (role='standard' OR (role='null' AND std_epoch=1))""", con)
    fl = pd.read_sql("SELECT * FROM be_floors WHERE filter='hrg'", con)
    fig, axs = plt.subplots(2, 1, figsize=(ps.COL_DOUBLE, 4.8), sharex=True)
    for ax, roles in ((axs[0], ("standard",)), (axs[1], ("null",))):
        for i, (lab, g) in enumerate(nt[nt.role.isin(roles)].groupby("label")):
            g = g.copy()
            g["d"] = g.ew - g.groupby("mech_state").ew.transform("mean")
            sty = ps.series(i)
            ax.errorbar(g.t - JD0, g.d, yerr=g.err, fmt=sty["marker"], color=sty["color"], ms=3,
                        lw=0.6, label=lab)
        ax.axvspan(std_epoch_bjd(), nt.t.max() - JD0 + 5, color=ps.GRID, zorder=0)
        ax.axhline(0, color=ps.RULE, lw=0.6)
        ax.set_ylabel("EW − state mean (Å)")
        ax.legend(fontsize=6, ncol=4, loc="upper left")
    for r in fl.itertuples():
        g = nt[(nt.mech_state == r.mech_state) & (nt.std_epoch == 1)]
        if len(g):
            t0, t1 = g.t.min() - JD0, g.t.max() - JD0
            for ax in axs:
                ax.fill_between([t0, t1], -3 * r.sd_resid, 3 * r.sd_resid, color=ps.WARN, alpha=0.15, lw=0)
    axs[0].set_title("Standards (η Hya, θ Vir; Vega where observed): ±3× night-to-night scatter shaded",
                     fontsize=7)
    axs[1].set_title("BeSS-constant null-test stars, standards epoch", fontsize=7)
    axs[1].set_xlabel(f"BJD_TDB − {JD0:.0f}")
    fig.tight_layout()
    save(fig, "fig2_standards")


def fig3(con):
    p = pd.read_sql("SELECT * FROM be_bess_pairs", con)
    col = dict(con.execute("SELECT key, value FROM be_meta").fetchall()).get("bess_pair_column", "")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 3.4), gridspec_kw={"width_ratios": [1.3, 1]})
    for i, (role, g) in enumerate(p.groupby("role")):
        sty = ps.series(i)
        a1.errorbar(g.ew_bess, g.ew_rlmt, yerr=g.err_rlmt, fmt=sty["marker"], color=sty["color"], ms=3.5,
                    lw=0.6, label=f"{role} ({len(g)})")
    lim = [min(p.ew_bess.min(), p.ew_rlmt.min()) - 2, max(p.ew_bess.max(), p.ew_rlmt.max()) + 2]
    a1.plot(lim, lim, color=ps.RULE, lw=0.8, ls="--")
    a1.set_xlim(lim), a1.set_ylim(lim)
    fw = col.replace("ew_fwhm", "")
    a1.set_xlabel(f"BeSS EW, degraded to the hrg LSF ({fw} Å FWHM) (Å)")
    a1.set_ylabel("RLMT hrg nightly EW (Å)")
    a1.legend(fontsize=6.5)
    d = p["diff"]
    a2.hist(d, bins=30, color=ps.ACCENT, alpha=0.8)
    med, sd = float(np.median(d)), float(1.4826 * np.median(np.abs(d - np.median(d))))
    a2.axvline(med, color=ps.BAD, lw=0.9)
    a2.set_xlabel("RLMT − BeSS (Å)")
    a2.set_ylabel("pairs")
    a2.set_title(f"median {med:+.2f} Å, robust sd {sd:.2f} Å, N = {len(d)}", fontsize=7)
    fig.tight_layout()
    save(fig, "fig3_bess")


def fig4(con):
    nt = pd.read_sql("SELECT * FROM be_nightly_ew WHERE role='science' AND filter='hrg'", con)
    b = pd.read_sql("SELECT main_id, mjd, ew_native FROM be_bess_ew WHERE status='ok' AND role='science'", con)
    ev = pd.read_sql("SELECT * FROM be_events", con) if con.execute(
        "SELECT name FROM sqlite_master WHERE name='be_events'").fetchone() else pd.DataFrame()
    stars = pd.read_sql("SELECT main_id, label FROM be_sample WHERE role='science' ORDER BY label", con)
    cont = pd.read_sql("SELECT * FROM be_continuum", con)
    varying = set(pd.read_sql("SELECT main_id FROM be_continuum_star WHERE varies=1", con).main_id)
    n = len(stars)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(ps.COL_DOUBLE, 1.45 * nrow), squeeze=False)
    for ax, s in zip(axs.flat, stars.itertuples()):
        g = nt[nt.main_id == s.main_id]
        for std, mk, fc in ((0, "o", "none"), (1, "o", ps.ACCENT)):
            gg = g[g.std_epoch == std]
            ax.errorbar(gg.t - JD0, gg.ew_cal, yerr=gg.err, fmt=mk, ms=2.5, mfc=fc, color=ps.ACCENT, lw=0.5)
        bb = b[b.main_id == s.main_id]
        ax.plot(bb.mjd + 2400000.5 - JD0, bb.ew_native, "x", color=ps.INK, ms=3)
        if len(ev) and s.main_id in set(ev.main_id):
            e = ev[ev.main_id == s.main_id].iloc[0]
            t = g[g.night == e.first_night].t
            if len(t):
                ax.axvline(t.iloc[0] - JD0, color=ps.BAD, lw=0.8)
        ax.axvline(std_epoch_bjd(), color=ps.FAINT, lw=0.6, ls=":")
        if s.main_id in varying:            # standing rule 6: continuum carried
            cc = cont[cont.main_id == s.main_id]
            ax2 = ax.twinx()
            ax2.plot(cc.t - JD0, cc.cont_resid_mag, "s", ms=2, color=ps.WARN)
            ax2.invert_yaxis()
            ax2.tick_params(labelsize=4.5, colors=ps.WARN)
        ax.set_title(s.label, fontsize=6.5, pad=2)
        ax.tick_params(labelsize=5)
        ax.invert_yaxis()
    for ax in list(axs.flat)[n:]:
        ax.axis("off")
    fig.supxlabel(f"BJD_TDB − {JD0:.0f}", fontsize=7)
    fig.supylabel("Hα EW (Å; emission up; QHY offset applied)  —  orange squares: continuum (mag, right)", fontsize=6.5)
    fig.tight_layout()
    save(fig, "fig4_ew_curves")


def fig5(con):
    ev = pd.read_sql("SELECT * FROM be_events", con)
    if not len(ev):
        return "cut (no standards-epoch event)"
    on = pd.read_sql("SELECT * FROM be_onsets", con)
    nt = pd.read_sql("SELECT * FROM be_nightly_ew WHERE role='science' AND filter='hrg' AND std_epoch=1", con)
    tess = pd.read_sql("SELECT * FROM be_tess WHERE status='ok'", con)
    top = ev.sort_values("max_abs_z", ascending=False).head(4)
    fig, axs = plt.subplots(len(top), 1, figsize=(ps.COL_DOUBLE, 1.9 * len(top)), squeeze=False)
    for ax, e in zip(axs[:, 0], top.itertuples()):
        g = nt[(nt.main_id == e.main_id) & (nt.mech_state == e.mech_state)]
        ax.errorbar(g.t - JD0, g.ew, yerr=g.err, fmt="o", ms=3, color=ps.ACCENT, lw=0.6)
        o = on[(on.main_id == e.main_id) & (on.mech_state == e.mech_state)]
        if len(o):
            o = o.iloc[0]
            ax.axvspan(o.t0_lo - JD0, o.t0_hi - JD0, color=ps.WARN, alpha=0.3, lw=0)
        ax.invert_yaxis()
        ax.set_ylabel("EW (Å)")
        ax.set_title(f"{e.star}: {e.sign}, first flagged {e.first_night}", fontsize=7)
        # TESS only where a sector overlaps the plotted nights
        tt = tess[(tess.main_id == e.main_id) & (tess.t_max_mjd + 2400000.5 >= g.t.min())
                  & (tess.t_min_mjd + 2400000.5 <= g.t.max())]
        for r in tt.itertuples():
            from astropy.io import fits
            with fits.open(C.PROJECT / r.lc_file) as hl:
                d = hl[1].data
                q = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
                tb = d["TIME"][q] + 2457000.0 - JD0
                fl = d["PDCSAP_FLUX"][q] / np.nanmedian(d["PDCSAP_FLUX"][q])
            ax2 = ax.twinx()
            ax2.plot(tb, (fl - 1) * 1e3, ",", color=ps.FAINT, alpha=0.5)
            ax2.set_ylabel(f"TESS S{r.sector} (ppt)", fontsize=6)
    axs[-1, 0].set_xlabel(f"BJD_TDB − {JD0:.0f}")
    fig.tight_layout()
    save(fig, "fig5_events")
    return "built"


def fig6(con):
    sr = pd.read_sql("SELECT * FROM be_search", con)
    if not int(sr.get("period_claim", pd.Series([0])).sum()):
        return "table (no significant period)"
    pg = pd.read_sql("SELECT * FROM be_periodograms", con)
    inj = pd.read_sql("SELECT * FROM be_injection WHERE channel='hrg'", con)
    sig = sr[sr.period_claim == 1]
    fig, axs = plt.subplots(len(sig), 2, figsize=(ps.COL_DOUBLE, 1.8 * len(sig)), squeeze=False)
    for (a1, a2), r in zip(axs, sig.itertuples()):
        g = pg[pg.main_id == r.main_id]
        a1.plot(g.freq, g.power, color=ps.ACCENT, lw=0.7)
        a1.set_title(f"{r.star}: P = {r.best_period_d:.2f} d, FAP {r.fap_global:.3g}", fontsize=7)
        a1.set_xlabel("frequency (d⁻¹)"), a1.set_ylabel("power")
        ii = inj[inj.main_id == r.main_id].pivot_table(index="amp_sigma", columns="period_d", values="frac")
        a2.contourf(ii.columns, ii.index * r.sigma_night_A, ii.values, levels=[0, .5, .9, 1.01],
                    colors=[ps.PAPER, ps.WISP, ps.SECOND])
        a2.set_xscale("log"), a2.set_yscale("log")
        a2.set_xlabel("period (d)"), a2.set_ylabel("amplitude (Å)")
    fig.tight_layout()
    save(fig, "fig6_periods")
    return "built"


def main() -> None:
    import argparse
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    con = C.be_db()
    with ps.context("print"):
        fig1(con)
        fig2(con)
        fig3(con)
        fig4(con)
        s5 = fig5(con)
        s6 = fig6(con)
    for k, v in (("fig5", s5), ("fig6", s6)):
        con.execute("INSERT OR REPLACE INTO be_meta VALUES (?,?)", (f"figure_{k}", v))
    con.commit()
    print("fig5:", s5, "| fig6:", s6)


if __name__ == "__main__":
    main()
