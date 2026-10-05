#!/usr/bin/env python
"""be_cadence — the cadence comparison of the abstract's second sentence: how
often the RLMT grism sampled each sample star's Hα line, against how often
BeSS did, over the SAME nights (Seat 6 edit 2, 2026-10-05).

PRE-DECLARED (written before any ratio was computed)
---------------------------------------------------
* Window: for each (star, season), our first to our last grism night of that
  season.  RLMT nights = distinct nights with a frame not excluded by the
  Step-0 disposition.  BeSS nights = distinct UT dates of BeSS records covering
  Hα inside the same window (the full BeSS index, nv_bess, not only the spectra
  fetched — a record is a night of BeSS sampling whether or not we measured it).
* Ratio = RLMT nights / max(BeSS nights, 1).  A window with no BeSS spectrum
  is counted with one, so the ratio is a LOWER bound there, and it is flagged.
* "BeSS-sparse" (a property of BeSS alone, not of our sampling): a star whose
  season windows hold BeSS Hα spectra with a median gap > 14 d between
  consecutive dates, or fewer than two dates.  A star is sparse if all its
  season windows are.
* Reported: the median ratio over all sample stars, and over the sparse ones,
  with the count of sparse stars.

Writes ``be_cadence`` (per star-season), ``be_cadence_star`` (per star) and
``notes/cadence/cadence.md``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_cadence.py"
SPARSE_GAP_D = 14.0


def main() -> None:
    import argparse
    argparse.ArgumentParser(description=__doc__.split(chr(10))[0]).parse_args()
    con = C.be_db()
    fr = pd.read_sql("""SELECT main_id, label, season, night FROM be_frames
                        WHERE role='science' AND disposition!='exclude'""", con)
    nv = C.novelty_ro()
    bess = pd.read_sql("SELECT main_id, substr(date,1,10) AS d FROM nv_bess WHERE covers_ha=1", nv)
    rows = []
    for (mid, label, season), g in fr.groupby(["main_id", "label", "season"]):
        n0, n1 = g.night.min(), g.night.max()
        # our nights are local-noon labels; a BeSS UT date d belongs to night d-1..d
        dates = sorted(set(bess[(bess.main_id == mid) & (bess.d >= n0) & (bess.d <= n1)].d))
        gaps = np.diff(pd.to_datetime(dates).values).astype("timedelta64[D]").astype(float) if len(dates) > 1 else []
        med_gap = float(np.median(gaps)) if len(gaps) else np.nan
        rows.append(dict(main_id=mid, star=label, season=season, first=n0, last=n1,
                         rlmt_nights=g.night.nunique(), bess_nights=len(dates),
                         bess_median_gap_d=med_gap,
                         sparse=int(len(dates) < 2 or med_gap > SPARSE_GAP_D),
                         ratio=g.night.nunique() / max(len(dates), 1),
                         ratio_lower_bound=int(len(dates) == 0)))
    cs = pd.DataFrame(rows)
    st = (cs.groupby(["main_id", "star"])
          .agg(rlmt_nights=("rlmt_nights", "sum"), bess_nights=("bess_nights", "sum"),
               sparse=("sparse", "min"), seasons=("season", lambda x: ",".join(x)))
          .reset_index())
    st["ratio"] = st.rlmt_nights / st.bess_nights.clip(lower=1)
    st["ratio_lower_bound"] = (st.bess_nights == 0).astype(int)
    cs.to_sql("be_cadence", con, if_exists="replace", index=False)
    st.to_sql("be_cadence_star", con, if_exists="replace", index=False)
    summ = dict(n_stars=len(st), n_sparse=int(st["sparse"].sum()),
                median_ratio_all=float(st.ratio.median()),
                median_ratio_sparse=float(st[st["sparse"] == 1].ratio.median()) if st["sparse"].any() else np.nan,
                total_rlmt=int(st.rlmt_nights.sum()), total_bess=int(st.bess_nights.sum()))
    for k, v in summ.items():
        con.execute("INSERT OR REPLACE INTO be_meta VALUES (?,?)", (f"cadence_{k}", str(v)))
    con.commit()
    body = (C.md_table(st.sort_values("ratio", ascending=False)[
        ["star", "seasons", "rlmt_nights", "bess_nights", "ratio", "ratio_lower_bound", "sparse"]], 1)
        + "\n#### Per star-season\n\n" + C.md_table(cs[["star", "season", "first", "last", "rlmt_nights",
                                                        "bess_nights", "bess_median_gap_d", "ratio", "sparse"]], 1)
        + f"\nRule (pre-declared in the script header): ratio = RLMT nights / max(BeSS Hα nights, 1) over our own "
          f"season windows; BeSS-sparse = median gap between BeSS Hα dates > {SPARSE_GAP_D:g} d or < 2 dates in "
          f"every season window. **{summ['n_sparse']} of {summ['n_stars']} stars are BeSS-sparse; median ratio "
          f"{summ['median_ratio_sparse']:.1f} on those, {summ['median_ratio_all']:.1f} over all "
          f"{summ['n_stars']}** ({summ['total_rlmt']} RLMT vs {summ['total_bess']} BeSS nights in total).\n")
    C.write_md(C.NOTES / "cadence" / "cadence.md", "RLMT vs BeSS Hα sampling over the same seasons", body, SCRIPT)
    print(summ)


if __name__ == "__main__":
    main()
