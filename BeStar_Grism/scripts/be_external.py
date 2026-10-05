#!/usr/bin/env python
"""be_external — BE-S13 external data for the re-drawn sample: every BeSS Hα
spectrum of every sample star across the campaign, re-measured on the PAPER'S
windows, and every TESS light curve (MAST) whose sector overlaps our seasons.

WHY
---
The BeSS one-to-one is the figure the paper stands on (strategy §5.6), and it
is only a validation if both sides are measured the same way: identical EW
windows (strategy §4 Step 7: line 6520–6610 Å, continuum from a straight line
through 6480–6520 and 6620–6660 Å) and matched resolution (BeSS degraded to
our measured line-spread function, BE-S5).  The novelty package measured BeSS
on its own narrower windows for a different question (is there emission?);
those numbers are NOT reused for the validation.

TESS (strategy §4 Step 13): the event figure needs to know, per star, whether a
TESS sector was running while we observed.  SPOC 2-min / 20-s light curves are
used where MAST has them (bright Be stars saturate the FFIs; SPOC apertures
are built for that); nothing is fitted here, only fetched and summarised.

SOURCES (cited in the paper)
  BeSS: Neiner et al. 2011, AJ 142, 149 — SSA service basebe.obspm.fr; the
        index of records was pulled 2026-10-03 by the novelty package
        (novelty.sqlite nv_bess); spectra are fetched here by their acref.
  TESS: MAST (astroquery.mast.Observations), SPOC light curves (Jenkins et al.
        2016, SPIE 9913, 99133E).

SUBCOMMANDS
  bess-fetch     download the spectra (cached under products/external/bess)
  bess-measure   EW on the paper windows, native and degraded (--fwhm list, Å)
  tess           query + download overlapping SPOC light curves; summary table
  notes          emitted tables
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

SCRIPT = "BeStar_Grism/scripts/be_external.py"
EXT = C.PROJECT / "products" / "external"
OUT = C.NOTES / "external"
WIN_LINE = (6520.0, 6610.0)                       # strategy §4 Step 7
WIN_CONT = ((6480.0, 6520.0), (6620.0, 6660.0))
MIN_FLANK_A = 15.0                                # Å of each continuum window
CAMPAIGN = ("2024-11-01", "2026-07-31")
POLITE_S = 0.4


def get(url: str, path: Path, tries: int = 4) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MACRO-RLMT Be-star S13"})
            with urllib.request.urlopen(req, timeout=180) as r:
                path.write_bytes(r.read())
            time.sleep(POLITE_S)
            return path
        except Exception as exc:                  # network: retry, then fail loudly
            last = exc
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"GET failed: {url}: {last}")


def bess_rows() -> pd.DataFrame:
    con = C.be_db()
    sample = pd.read_sql("SELECT main_id, label, role FROM be_sample", con)
    nv = C.novelty_ro()
    b = pd.read_sql("SELECT main_id, bess_id, mjd, date, acref FROM nv_bess WHERE covers_ha=1 "
                    "AND date BETWEEN ? AND ?", nv, params=CAMPAIGN)
    return b.merge(sample, on="main_id")


def bess_path(r) -> Path:
    h = hashlib.sha1(r.acref.encode()).hexdigest()[:12]
    return EXT / "bess" / f"bess_{r.bess_id.replace('BeSS:', '')}_{h}.fits"


def cmd_bess_fetch(_a) -> None:
    rows = bess_rows()
    for r in rows.itertuples():
        get(r.acref, bess_path(r))
    print(f"bess-fetch: {len(rows)} spectra cached")


def read_bess(path: Path):
    """Same two BeSS layouts as the novelty package's reader."""
    from astropy.io import fits
    with fits.open(path) as hl:
        h0 = hl[0].header
        meta = {"observer": str(h0.get("OBSERVER", ""))[:60],
                "rp": h0.get("SPE_RPOW") or h0.get("BSS_ITRP") or h0.get("BSS_ESRP")}
        if hl[0].data is not None and np.ndim(hl[0].data) >= 1:
            f = np.asarray(hl[0].data, float).ravel()
            w = h0["CRVAL1"] + h0["CDELT1"] * (np.arange(f.size) + 1 - h0.get("CRPIX1", 1))
        elif len(hl) > 1 and hasattr(hl[1], "columns"):
            w, f = np.asarray(hl[1].data.field(0), float), np.asarray(hl[1].data.field(1), float)
        else:
            raise ValueError("unrecognised BeSS layout")
    ok = np.isfinite(w) & np.isfinite(f)
    o = np.argsort(w[ok])
    return w[ok][o], f[ok][o], meta


def degrade(w, f, fwhm_a: float):
    """Convolve onto a uniform grid with a Gaussian of FWHM ``fwhm_a`` (Å),
    in quadrature with the spectrum's own resolution assumed negligible."""
    dl = float(np.median(np.diff(w)))
    g = np.arange(w[0], w[-1], dl)
    fg = np.interp(g, w, f)
    s = fwhm_a / 2.3548 / dl
    k = np.arange(-int(4 * s) - 1, int(4 * s) + 2)
    ker = np.exp(-0.5 * (k / s) ** 2)
    ker /= ker.sum()
    return g, np.convolve(fg, ker, mode="same")


def ew_paper(w, f):
    """EW on the paper's windows: straight-line continuum through the two
    flanks (5–95% clipped), EW = ∫(1 − F/Fc) over the line window, Å,
    positive = absorption.  Returns (ew, err, cont_rms) or None."""
    cmask = ((w >= WIN_CONT[0][0]) & (w <= WIN_CONT[0][1])) | ((w >= WIN_CONT[1][0]) & (w <= WIN_CONT[1][1]))
    lmask = (w >= WIN_LINE[0]) & (w <= WIN_LINE[1])
    # Coverage rule (fixed before any comparison with our EWs): the whole
    # line window, and at least MIN_FLANK_A of EACH continuum window.  Most
    # BeSS uploads start at 6500 Å, so the blue flank is often 6500–6520 only;
    # the continuum is then fitted on the part that exists, and the one-to-one
    # carries that placement difference as part of what it measures.
    blue = np.ptp(w[(w >= WIN_CONT[0][0]) & (w <= WIN_CONT[0][1])]) if np.any(
        (w >= WIN_CONT[0][0]) & (w <= WIN_CONT[0][1])) else 0.0
    red = np.ptp(w[(w >= WIN_CONT[1][0]) & (w <= WIN_CONT[1][1])]) if np.any(
        (w >= WIN_CONT[1][0]) & (w <= WIN_CONT[1][1])) else 0.0
    if (w.min() > WIN_LINE[0] or w.max() < WIN_LINE[1] or blue < MIN_FLANK_A or red < MIN_FLANK_A
            or cmask.sum() < 8 or lmask.sum() < 5):
        return None
    wc, fc = w[cmask], f[cmask]
    p = np.polyfit(wc, fc, 1)
    r = fc - np.polyval(p, wc)
    lo, hi = np.percentile(r, [5, 95])
    k = (r >= lo) & (r <= hi)
    p = np.polyfit(wc[k], fc[k], 1)
    cont = np.polyval(p, w)
    if np.any(cont[lmask | cmask] <= 0):
        return None
    n = f / cont
    rms = float(np.std(n[cmask][k], ddof=2))
    dl = float(np.median(np.diff(w[lmask])))
    ew = float(np.trapz(1 - n[lmask], w[lmask]))
    err = float(np.hypot(rms * dl * np.sqrt(lmask.sum()), rms / np.sqrt(k.sum()) * (WIN_LINE[1] - WIN_LINE[0])))
    return ew, err, rms


def cmd_bess_measure(a) -> None:
    rows = bess_rows()
    out = []
    for r in rows.itertuples():
        p = bess_path(r)
        rec = dict(main_id=r.main_id, bess_id=r.bess_id, mjd=r.mjd, date=r.date, role=r.role)
        try:
            w, f, meta = read_bess(p)
            rec.update(observer=meta["observer"])
            try:
                rec["rp"] = float(meta["rp"]) if meta["rp"] not in (None, "") else np.nan
            except (TypeError, ValueError):
                rec["rp"] = np.nan
            m = ew_paper(w, f)
            rec.update(ew_native=m[0] if m else np.nan, ew_err_native=m[1] if m else np.nan,
                       status="ok" if m else "window_not_covered")
            for fw in a.fwhm:
                g, fg = degrade(w, f, fw)
                md = ew_paper(g, fg)
                rec[f"ew_fwhm{fw:g}"] = md[0] if md else np.nan
        except Exception as exc:
            rec.update(status=f"read_error: {type(exc).__name__}")
        out.append(rec)
    df = pd.DataFrame(out)
    con = C.be_db()
    df.to_sql("be_bess_ew", con, if_exists="replace", index=False)
    con.commit()
    print(df.status.value_counts().to_string())


def cmd_tess(_a) -> None:
    """SPOC light curves overlapping the campaign, per sample star."""
    from astropy.io import fits
    from astroquery.mast import Observations
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    con = C.be_db()
    sample = pd.read_sql("SELECT main_id, label, role, ra, dec FROM be_sample WHERE role!='standard'", con)
    t0, t1 = (pd.Timestamp(x).to_julian_date() - 2400000.5 for x in CAMPAIGN)
    rows = []
    for s in sample.itertuples():
        obs = Observations.query_criteria(coordinates=SkyCoord(s.ra * u.deg, s.dec * u.deg),
                                          radius=0.005 * u.deg, obs_collection="TESS",
                                          dataproduct_type="timeseries")
        for o in obs:
            if not (o["t_max"] >= t0 and o["t_min"] <= t1):
                continue
            rec = dict(main_id=s.main_id, star=s.label, role=s.role, obs_id=o["obs_id"],
                       sector=int(o["sequence_number"]), t_min_mjd=float(o["t_min"]),
                       t_max_mjd=float(o["t_max"]), provenance=str(o["provenance_name"]))
            prods = Observations.filter_products(Observations.get_product_list(o),
                                                 productSubGroupDescription="LC")
            if len(prods) == 0:
                rec["status"] = "no_lc_product"
                rows.append(rec)
                continue
            uri = prods[0]["dataURI"]
            path = EXT / "tess" / Path(uri).name
            get("https://mast.stsci.edu/api/v0.1/Download/file?uri=" + uri, path)
            with fits.open(path) as hl:
                d = hl[1].data
                q = (d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])
                fl = d["PDCSAP_FLUX"][q] / np.nanmedian(d["PDCSAP_FLUX"][q])
            rec.update(lc_file=str(path.relative_to(C.PROJECT)), n_cad=int(q.sum()),
                       rms_ppt=float(np.std(fl) * 1e3),
                       p2p_ppt=float(np.median(np.abs(np.diff(fl))) * 1e3), status="ok")
            rows.append(rec)
        print(s.label, sum(r["main_id"] == s.main_id for r in rows))
    df = pd.DataFrame(rows)
    # overlap with OUR nights of the same star
    fr = pd.read_sql("SELECT main_id, night, MIN(bjd_tdb) b FROM be_frames WHERE disposition!='exclude' "
                     "GROUP BY main_id, night", con)
    fr["mjd"] = fr.b - 2400000.5
    df["our_nights_in_sector"] = [int(((fr.main_id == r.main_id) & (fr.mjd >= r.t_min_mjd) &
                                       (fr.mjd <= r.t_max_mjd)).sum()) for r in df.itertuples()]
    df.to_sql("be_tess", con, if_exists="replace", index=False)
    con.execute("INSERT OR REPLACE INTO be_meta VALUES ('tess_pulled_utc', ?)",
                (dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),))
    con.commit()
    print(df.drop(columns=[c for c in ("lc_file", "obs_id") if c in df]).to_string())


def cmd_notes(_a) -> None:
    con = C.be_db()
    b = pd.read_sql("""SELECT s.role, s.label AS star, COUNT(e.bess_id) AS bess_spectra,
                       SUM(e.status='ok') AS measured, COUNT(DISTINCT e.observer) AS observers,
                       MIN(e.date) AS first, MAX(e.date) AS last,
                       AVG(e.ew_native) AS ew_mean_A
                       FROM be_sample s LEFT JOIN be_bess_ew e USING (main_id)
                       WHERE s.role!='standard' GROUP BY s.main_id ORDER BY s.role DESC, bess_spectra DESC""", con)
    body = C.md_table(b) + ("\nEW on the paper windows (line 6520–6610 Å; continuum 6480–6520 + 6620–6660 Å), "
                            "native BeSS resolution; positive = absorption. Campaign window "
                            f"{CAMPAIGN[0]} → {CAMPAIGN[1]}. Source: BeSS (Neiner et al. 2011), index pulled "
                            "2026-10-03.\n")
    if con.execute("SELECT name FROM sqlite_master WHERE name='be_tess'").fetchone():
        t = pd.read_sql("""SELECT star, role, sector, provenance, ROUND(t_min_mjd,1) AS t_min_mjd,
                           ROUND(t_max_mjd,1) AS t_max_mjd, n_cad, rms_ppt, our_nights_in_sector, status
                           FROM be_tess ORDER BY star, sector""", con)
        meta = dict(con.execute("SELECT key, value FROM be_meta").fetchall())
        body += ("\n### TESS sectors overlapping the campaign (MAST, SPOC light curves)\n\n" + C.md_table(t) +
                 f"\nPulled {meta.get('tess_pulled_utc')}. `our_nights_in_sector` = RLMT grism nights of the same "
                 "star inside the sector's time span — the coincidences the event analysis can use.\n")
    C.write_md(OUT / "external.md", "BeSS spectra of the re-drawn sample on the paper windows (BE-S13)", body, SCRIPT)
    print("written", OUT)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["bess-fetch", "bess-measure", "tess", "notes", "all"])
    ap.add_argument("--fwhm", type=float, nargs="*", default=[])
    a = ap.parse_args(argv)
    if a.cmd in ("bess-fetch", "all"):
        cmd_bess_fetch(a)
    if a.cmd in ("bess-measure", "all"):
        cmd_bess_measure(a)
    if a.cmd in ("tess", "all"):
        cmd_tess(a)
    if a.cmd in ("notes", "all"):
        cmd_notes(a)


if __name__ == "__main__":
    main()
