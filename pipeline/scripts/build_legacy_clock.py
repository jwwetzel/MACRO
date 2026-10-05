#!/usr/bin/env python
"""RIG-L1-clock-audit — the legacy archive's clock, per season, from eclipses.

WHAT THIS SCRIPT DOES
---------------------
For every observing season of the legacy archive it asks: was the
acquisition PC's clock right to 60 s?  The evidence is the primary eclipses
of post-common-envelope binaries (HW Vir-type sdB + dM; WD + dM) that the
census found in the archive, timed on OUR stamps and compared with
predictions from eclipse times MEASURED IN TESS photometry of the same stars
(absolute BJD_TDB; see ``macro_legacy.clock`` for why TESS and not a
literature ephemeris).  Stages, each resumable:

``tess``        find the 2-min SPOC light curves of every admitted clock
                target at MAST, download them once (cached under
                ``products/legacy/external/tess/``), and measure the
                primary-eclipse mid-time in each half-sector (orbit) by
                folding and fitting.  -> ``clock_tess_times``
``events``      every legacy run that covers a primary eclipse predicted
                from the TESS times with :data:`clock.RUN_COVER_P` on both
                sides; out-of-range predictions are recorded, not timed.
                -> ``clock_events``
``photometry``  per event: identify the field on one reference frame by a
                Gaia DR3 fit (the legacy headers carry no plate scale or
                parity), choose comparison stars, and measure every frame
                (S3b's worker, unchanged).  -> ``clock_flux``
``fit``         differential light curve, blind mid-time fit, block
                bootstrap, analysis variants, injection–recovery (S3b's
                ``analyse_event``, unchanged); O − C under both stamp
                readings.  -> ``clock_fits``
``summary``     O − C per season (camera × acquisition software × July–June
                season), with every season of the archive listed — a season
                without a usable eclipse says so.  -> ``clock_seasons``

INPUTS AND OUTPUTS
------------------
Reads ``products/legacy/legacy.sqlite`` (``eclipsers``, ``vsx_matches``,
``runs``, ``frames``, ``time_audit``) and writes the ``clock_*`` tables into
the same database.  The archive is opened READ-ONLY; no pixel is written.

The photometry and the fit are S3b's (``build_s3b_clock_transits.py``)
imported unchanged, so the two clock audits — RLMT and legacy — measure
with one method.  What differs is said where it differs: no master
calibration (the legacy archive holds none; most frames were calibrated in
place by the scheduler, ``CALSTAT``), a fixed saturation veto (no S2
measurement exists for these cameras), and Gaia astrometry in place of a
header WCS.

USAGE
-----
    /opt/miniconda3/envs/rlmt-checks/bin/python \\
        pipeline/scripts/build_legacy_clock.py                 # all stages
    … build_legacy_clock.py --stage fit --workers 8
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import sqlite3
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "scripts"))

from macro_core import clock_transits as ct                   # noqa: E402
from macro_core.timing import bjd_tdb_from_utc                 # noqa: E402
from macro_legacy import clock as lk                           # noqa: E402
import build_s3b_clock_transits as s3b                         # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "legacy" / "legacy.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/legacy-archive")
TESS_DIR = REPO_ROOT / "products" / "legacy" / "external" / "tess"
GAIA_DIR = REPO_ROOT / "products" / "legacy" / "external" / "gaia"

#: Saturation veto for legacy frames (frame ADU).  The cameras' 16-bit ADCs
#: clip at 65,535; most legacy frames were bias/dark-subtracted and
#: flat-divided in place by the scheduler, which moves a clipped pixel by up
#: to ~10 %.  50,000 ADU keeps every clipped pixel above the veto.  No S2
#: measurement exists for these cameras, so this is a policy, stated.
LEGACY_VETO_ADU = 50000.0

#: Binning passed to S3b's native-peak scaling: the legacy cameras bin ON
#: CHIP (charge summed before the ADC), so a binned pixel clips at the ADC
#: limit itself and no native-pixel rescaling applies.
ONCHIP_BINNING = 1

#: Gaia cone radius (deg) — covers the largest legacy field (AC4040,
#: 36.7' square) with margin.
GAIA_RADIUS_DEG = 0.40

#: Nominal HW Vir-type eclipse geometry (S3b's NSVS 07826147 template,
#: after For et al. 2010) used as the STARTING shape for every sdB + dM
#: system; duration and depth are refitted per event and the mid-time of a
#: symmetric eclipse does not depend on them (S3b's variants prove it).
SHAPE_HW = dict(a_over_r=5.90, k=0.916, b=0.350, u1=0.25, u2=0.20, f1=0.9,
                free=("a_over_r", "f1"), baseline_order=2)
#: WD + dM template (S3b's GK Vir shape, after Parsons et al. 2012).
SHAPE_WD = dict(a_over_r=107.0, k=9.1, b=0.0, u1=0.30, u2=0.20, f1=0.8,
                free=("a_over_r", "f1"), baseline_order=2)
MIN_DEPTH = 0.05


def log(msg: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def shape_for(vsx_type: str, period: float) -> ct.EventShape:
    base = SHAPE_WD if "WD" in vsx_type or "UV" in vsx_type else SHAPE_HW
    return ct.EventShape(period=period, **base)


def replace(con, name: str, df: pd.DataFrame) -> None:
    with con:
        con.execute(f'DROP TABLE IF EXISTS "{name}"')
        df.to_sql(name, con, index=False)


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
def clock_targets(con) -> pd.DataFrame:
    """Admitted clock targets: one row per VSX star (aliases merged).

    Every legacy target key matched to the same VSX object is the same star
    (e.g. 'HS0705+6700' and 'HS0705+67' are V0470 Cam), so targets are
    grouped by VSX name and carry the list of their legacy keys.
    """
    e = pd.read_sql_query("""
        SELECT e.target_key, e.vsx_name, e.vsx_type, e.period_d, e.epoch_hjd,
               v.vsx_ra_deg, v.vsx_dec_deg
        FROM eclipsers e JOIN vsx_matches v USING (target_key)
        WHERE e.is_eclipsing = 1 AND e.period_d IS NOT NULL""", con)
    e = e[[lk.is_clock_type(t) for t in e["vsx_type"]]]
    rows = []
    for name, g in e.groupby("vsx_name"):
        r = g.iloc[0]
        rows.append(dict(vsx_name=name, vsx_type=r.vsx_type,
                         period_vsx=float(r.period_d),
                         epoch_vsx=float(r.epoch_hjd) if pd.notna(
                             r.epoch_hjd) else None,
                         ra_deg=float(r.vsx_ra_deg),
                         dec_deg=float(r.vsx_dec_deg),
                         target_keys=",".join(sorted(g["target_key"]))))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Stage: tess
# ---------------------------------------------------------------------------
def tess_files(name: str, ra: float, dec: float) -> list[Path]:
    """SPOC 2-min light-curve files for one star (downloaded once)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.mast import Observations
    d = TESS_DIR / name.replace(" ", "_")
    d.mkdir(parents=True, exist_ok=True)
    have = sorted(d.rglob("*_lc.fits"))
    if have:
        return have
    obs = Observations.query_criteria(
        coordinates=SkyCoord(ra, dec, unit="deg"), radius=0.005 * u.deg,
        obs_collection="TESS", dataproduct_type="timeseries")
    obs = obs[[p == "SPOC" and float(e) == 120.0
               for p, e in zip(obs["provenance_name"], obs["t_exptime"])]]
    if len(obs) == 0:
        return []
    prods = Observations.get_product_list(obs)
    prods = prods[[str(s) == "LC" for s in prods["productSubGroupDescription"]]]
    Observations.download_products(prods, download_dir=str(d))
    return sorted(d.rglob("*_lc.fits"))


def read_tess(path: Path) -> dict:
    """Time (BJD_TDB), normalised PDCSAP flux, sector, for quality-0 points."""
    from astropy.io import fits
    with fits.open(path) as h:
        d = h[1].data
        sector = int(h[0].header["SECTOR"])
        t = np.asarray(d["TIME"], float) + 2457000.0
        f = np.asarray(d["PDCSAP_FLUX"], float)
        q = np.asarray(d["QUALITY"], int)
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0)
    t, f = t[ok], f[ok]
    return {"sector": sector, "t": t, "f": f / np.median(f)}


def time_tess_chunk(t, f, t_guess, period, shape, seed) -> dict | None:
    """Mid-time of the primary nearest the middle of one TESS chunk.

    Folds the chunk onto the cycle of its middle predicted primary, bins at
    :data:`clock.FOLD_BIN_S`, fits S3b's eclipse model blind, and takes the
    error from S3b's block bootstrap.  Returns None when the fit fails.
    """
    tmid = 0.5 * (t.min() + t.max())
    t_ref = t_guess + round((tmid - t_guess) / period) * period
    tf, keep = lk.fold_primary(t, t_ref, period)
    if keep.sum() < 200:
        return None
    tb, fb, eb, _n = lk.bin_series(tf[keep], f[keep],
                                   lk.FOLD_BIN_S / lk.SECONDS_PER_DAY)
    ex = np.full(len(tb), 120.0 / lk.SECONDS_PER_DAY)
    fit = ct.fit_event(tb, fb, eb, ex, shape)
    if fit["status"] != ct.STATUS_OK or fit["depth"] < MIN_DEPTH:
        return None
    boot = ct.bootstrap_t0(tb, fb, eb, ex, shape, fit, seed=seed)
    sig = boot["sigma"] * lk.SECONDS_PER_DAY
    if not np.isfinite(sig):
        return None
    return {"t0": fit["t0"], "sig_s": max(sig, 1e-3), "depth": fit["depth"],
            "n_points": int(keep.sum()), "chi2": fit["chi2"],
            "dof": fit["dof"]}


def stage_tess(con) -> None:
    targets = clock_targets(con)
    rows = []
    for tg in targets.itertuples():
        files = tess_files(tg.vsx_name, tg.ra_deg, tg.dec_deg)
        log(f"tess: {tg.vsx_name}: {len(files)} light-curve files")
        shape = shape_for(tg.vsx_type, tg.period_vsx)
        t_guess = tg.epoch_vsx
        for path in files:
            lc = read_tess(path)
            t, f = lc["t"], lc["f"]
            # Two orbits per sector, split at the mid-sector data gap.
            gaps = np.flatnonzero(np.diff(t) > 1.0)
            cut = gaps[np.argmax(np.diff(t)[gaps])] + 1 if len(gaps) else None
            chunks = [(t[:cut], f[:cut]), (t[cut:], f[cut:])] if cut \
                else [(t, f)]
            for k, (tc, fc) in enumerate(chunks):
                if len(tc) < 500:
                    continue
                seed = int(hashlib.sha1(
                    f"{tg.vsx_name}{lc['sector']}{k}".encode()
                ).hexdigest()[:6], 16)
                r = time_tess_chunk(tc, fc, t_guess, tg.period_vsx, shape,
                                    seed)
                rows.append(dict(
                    vsx_name=tg.vsx_name, sector=lc["sector"], orbit=k + 1,
                    t_start=float(tc.min()), t_end=float(tc.max()),
                    t0_bjd=None if r is None else r["t0"],
                    sig_s=None if r is None else r["sig_s"],
                    depth=None if r is None else r["depth"],
                    n_points=None if r is None else r["n_points"],
                    chi2=None if r is None else r["chi2"],
                    dof=None if r is None else r["dof"],
                    file=path.name))
    replace(con, "clock_targets", targets)
    replace(con, "clock_tess_times", pd.DataFrame(rows))
    log(f"tess: {sum(r['t0_bjd'] is not None for r in rows)} of {len(rows)} "
        "half-sector eclipse times measured")


# ---------------------------------------------------------------------------
# Stage: events
# ---------------------------------------------------------------------------
def tess_reference(con) -> dict:
    """Per star: TESS times, cycles and sector count for the predictions."""
    tt = pd.read_sql_query("SELECT * FROM clock_tess_times "
                           "WHERE t0_bjd IS NOT NULL", con)
    tg = pd.read_sql_query("SELECT * FROM clock_targets", con)
    out = {}
    for r in tg.itertuples():
        g = tt[tt["vsx_name"] == r.vsx_name]
        if g.empty:
            continue
        anchor = float(g["t0_bjd"].iloc[0])
        cyc = lk.cycle_numbers(g["t0_bjd"], anchor, r.period_vsx)
        out[r.vsx_name] = dict(anchor=anchor, period=r.period_vsx,
                               times=g["t0_bjd"].to_numpy(),
                               sigs=(g["sig_s"] / lk.SECONDS_PER_DAY
                                     ).to_numpy(),
                               cycles=cyc, n_sectors=g["sector"].nunique())
    return out


def stage_events(con) -> None:
    tg = pd.read_sql_query("SELECT * FROM clock_targets", con)
    ref = tess_reference(con)
    runs = pd.read_sql_query("SELECT * FROM runs", con)
    rows = []
    for t in tg.itertuples():
        keys = t.target_keys.split(",")
        r = runs[runs["target_key"].isin(keys)]
        for run in r.itertuples():
            b0, b1 = bjd_tdb_from_utc([run.jd_first, run.jd_last],
                                      t.ra_deg, t.dec_deg)[0]
            if t.vsx_name not in ref:
                rows.append(dict(event_id=None, vsx_name=t.vsx_name,
                                 target_key=run.target_key, night=run.night,
                                 camera=run.camera, software=run.software,
                                 band=run.band, cycle=None, t_pred=None,
                                 status="no TESS reference"))
                continue
            R = ref[t.vsx_name]
            pad = lk.RUN_COVER_P * R["period"]
            for e in ct.events_in_window(b0 + pad, b1 - pad, R["anchor"],
                                         R["period"]):
                t_appr = R["anchor"] + e * R["period"]
                p = lk.reference_prediction(e, t_appr, R["times"],
                                            R["sigs"], R["cycles"],
                                            R["n_sectors"])
                eid = (f"{t.vsx_name.replace(' ', '')}_{run.night}_"
                       f"{run.camera.split()[0]}_{run.band}_{e}")
                row = dict(event_id=eid, vsx_name=t.vsx_name,
                           target_key=run.target_key, night=run.night,
                           camera=run.camera, software=run.software,
                           band=run.band, cycle=int(e),
                           run_first=float(run.jd_first),
                           run_last=float(run.jd_last),
                           t_pred=p.get("t_pred"),
                           sig_pred_s=p.get("sig_stat_s"),
                           sig_sys_s=p.get("sig_sys_s"),
                           sys_basis=p.get("sys_basis"),
                           ref_range=p["range"], gap_d=p["gap_d"],
                           status=("admitted" if p["range"] != "out_of_range"
                                   else "prediction out of range"))
                rows.append(row)
    ev = pd.DataFrame(rows)
    replace(con, "clock_events", ev)
    log(f"events: {int((ev['status'] == 'admitted').sum())} admitted of "
        f"{len(ev)} candidate eclipses")


# ---------------------------------------------------------------------------
# Stage: photometry
# ---------------------------------------------------------------------------
def gaia_field(name: str, ra: float, dec: float) -> dict:
    """Gaia DR3 cone around a clock target (cached as JSON)."""
    from macro_phot import gaia
    GAIA_DIR.mkdir(parents=True, exist_ok=True)
    cache = GAIA_DIR / f"{name.replace(' ', '_')}.json"
    if cache.exists():
        d = json.loads(cache.read_text())
        return {k: np.asarray(v) for k, v in d.items()}
    try:
        g = gaia.cone_query(ra, dec, GAIA_RADIUS_DEG)
        source = "Gaia archive TAP (gaiadr3.gaia_source)"
    except Exception as exc:                         # noqa: BLE001
        # The ESA archive refuses large cones with 'statement timeout'
        # while it is being rebuilt for DR4; VizieR serves the same Gaia
        # DR3 table (I/355/gaiadr3, positions at epoch 2016.0 as in the
        # archive), so the fallback changes the server, not the catalogue.
        log(f"gaia: archive failed for {name} ({type(exc).__name__}); "
            "using VizieR I/355/gaiadr3")
        g = gaia_vizier(ra, dec, GAIA_RADIUS_DEG)
        source = "VizieR I/355/gaiadr3"
    cache.write_text(json.dumps({k: v.tolist() for k, v in g.items()}))
    (GAIA_DIR / f"{name.replace(' ', '_')}.source.txt").write_text(source)
    return g


def gaia_vizier(ra: float, dec: float, radius_deg: float) -> dict:
    """Gaia DR3 cone from VizieR: the brightest 2,000 stars with G < 19,
    in the same dict layout as ``macro_phot.gaia.cone_query``."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    from macro_phot import gaia
    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "+Gmag"],
               column_filters={"Gmag": f"<{gaia.GAIA_G_MAX}"}, row_limit=2000)
    t = v.query_region(SkyCoord(ra, dec, unit="deg"),
                       radius=radius_deg * u.deg, catalog="I/355/gaiadr3")[0]
    return {"source_id": np.asarray(t["Source"], dtype=np.int64),
            "ra": np.asarray(t["RA_ICRS"], dtype=float),
            "dec": np.asarray(t["DE_ICRS"], dtype=float),
            "gmag": np.asarray(t["Gmag"], dtype=float)}


def prepare(archive: Path, paths: list[str], ra: float, dec: float,
            gaia_cat: dict, seed: int) -> dict:
    """Reference frame, Gaia identification, target and comparison stars.

    Mirrors S3b's ``prepare_series`` except for the identification: legacy
    headers have no plate scale or parity, so the reference frame's
    detections are fitted to Gaia DR3 (``macro_phot.gaia.
    identify_reference``, both parities tried) and the target is the
    detection whose fitted sky position lies within
    ``gaia.TARGET_ID_TOL_ARCSEC`` of the catalogue position.
    """
    from macro_phot import gaia
    mid = len(paths) // 2
    order = sorted(range(len(paths)), key=lambda i: abs(i - mid))
    last = "no usable reference frame"
    for i in order[:s3b.MAX_REF_TRIES]:
        try:
            raw, hdr = s3b.read_frame(archive / paths[i])
            data, sub, bkg, objs = s3b.detect(raw)
        except Exception as e:                       # noqa: BLE001
            last = f"reference unreadable: {type(e).__name__}"
            continue
        if len(objs) < 8:
            last = "fewer than 8 sources on the reference"
            continue
        xy = np.column_stack([objs["x"], objs["y"]]).astype(float)[:800]
        h, w = data.shape
        scale_guess = None
        try:
            scale_guess = (206.265 * float(hdr["XPIXSZ"])
                           / float(hdr["FOCALLEN"]))
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            pass
        fit_r = 0.5 * min(h, w) * scale_guess if scale_guess else None
        idn = gaia.identify_reference(xy, gaia_cat, ra, dec,
                                      ref_bright_xy=xy,
                                      fit_radius_arcsec=fit_r, seed=seed)
        if idn is None:
            last = "Gaia fit did not converge"
            continue
        sky = idn["ref_radec"]
        d = np.hypot((sky[:, 0] - ra) * math.cos(math.radians(dec)),
                     sky[:, 1] - dec) * 3600.0
        it = int(np.argmin(d))
        if d[it] > gaia.TARGET_ID_TOL_ARCSEC:
            last = "target not identified on the reference"
            continue
        scale = idn["scale_arcsec_per_px"]
        fwhm = float(np.median(s3b.fwhm_of(objs[:50])))
        radii = [max(ct.APERTURE_MIN_PX, f * fwhm)
                 for f in ct.APERTURE_FWHM_FACTORS]
        ann = (max(ct.ANNULUS_FWHM_FACTORS[0] * fwhm, radii[-1] + 2.0),
               max(ct.ANNULUS_FWHM_FACTORS[1] * fwhm, radii[-1] + 6.0))
        margin = s3b.EDGE_MARGIN_PX + ann[1]
        comps = []
        for j in range(min(len(objs), 800)):
            if j == it or len(comps) >= ct.N_COMP_CANDIDATES:
                continue
            x, y = float(objs["x"][j]), float(objs["y"][j])
            if min(x, y, w - x, h - y) < margin:
                continue
            xi, yi = int(round(x)), int(round(y))
            peak = float(raw[max(yi - 2, 0):yi + 3,
                             max(xi - 2, 0):xi + 3].max())
            if peak >= s3b.COMP_PEAK_HEADROOM * LEGACY_VETO_ADU:
                continue
            dd = np.hypot(objs["x"] - x, objs["y"] - y)
            if (objs["flux"][(dd < ann[1]) & (dd > 0)]
                    > 0.1 * objs["flux"][j]).any():
                continue
            if np.hypot(x - objs["x"][it], y - objs["y"][it]) < ann[1]:
                continue
            comps.append(j)
        if len(comps) < ct.MIN_COMPS:
            last = "fewer than 2 clean comparison stars"
            continue
        idx = [it] + comps
        return {"ok": True, "ref_path": paths[i], "fwhm_px": fwhm,
                "scale": scale, "parity": idn["parity"],
                "n_gaia": idn["n_matched"], "target_sep": float(d[it]),
                "radii": radii, "annulus": ann,
                # Every index here is < 800 (the comparison loop stops at
                # 800 and the target comes from the Gaia fit of xy itself).
                "star_xy": xy[idx].tolist(),
                "ref_bright": np.column_stack(
                    [objs["x"][:800], objs["y"][:800]]).tolist()}
    return {"ok": False, "reason": last}


def _prepare_job(job: dict) -> dict:
    """Worker wrapper around :func:`prepare` (one series)."""
    warnings.filterwarnings("ignore")
    out = prepare(Path(job["archive"]), job["paths"], job["ra"], job["dec"],
                  {k: np.asarray(v) for k, v in job["gaia"].items()},
                  job["seed"])
    out["sid"] = job["sid"]
    return out


def _measure_job(task: dict) -> dict:
    """Worker wrapper around S3b's ``measure_one`` that keeps the series id."""
    warnings.filterwarnings("ignore")
    r = s3b.measure_one(task)
    r["sid"] = task["sid"]
    return r


def stage_photometry(con, archive: Path, workers: int) -> None:
    """Prepare every series, then measure every frame in ONE worker pool.

    Two passes so that the expensive part (decompressing and registering
    ~10^4 frames, 4096^2 for the AC4040) runs on one long-lived pool
    instead of a pool per series.  Resumable: a series already in
    ``clock_series`` is skipped.
    """
    ev = pd.read_sql_query("SELECT * FROM clock_events "
                           "WHERE status = 'admitted'", con)
    tg = pd.read_sql_query("SELECT * FROM clock_targets",
                           con).set_index("vsx_name")
    with con:
        con.execute("CREATE TABLE IF NOT EXISTS clock_series (series_id TEXT "
                    "PRIMARY KEY, vsx_name TEXT, night TEXT, camera TEXT, "
                    "software TEXT, band TEXT, n_frames INTEGER, status TEXT, "
                    "ref_path TEXT, parity TEXT, n_gaia INTEGER, "
                    "target_sep_arcsec REAL, scale_arcsec REAL, "
                    "fwhm_px REAL, n_comps INTEGER)")
        con.execute("CREATE TABLE IF NOT EXISTS clock_flux (series_id TEXT, "
                    "path TEXT, jd_start REAL, exptime_s REAL, "
                    "bjd_mid_startconv REAL, fwhm_px REAL, airmass REAL, "
                    "star INTEGER, x REAL, y REAL, matched INTEGER, "
                    "on_frame INTEGER, peak_adu REAL, flux0 REAL, flux1 REAL,"
                    " flux2 REAL, err0 REAL, err1 REAL, err2 REAL)")
    done = {r[0] for r in con.execute("SELECT series_id FROM clock_series")}
    # One photometric series per (star, night, camera, band, run): every
    # admitted eclipse inside it shares the series.
    ev["series_id"] = [f"{v.replace(' ', '')}_{n}_{c.split()[0]}_{b}_"
                       f"{rf:.4f}" for v, n, c, b, rf in zip(
                           ev["vsx_name"], ev["night"], ev["camera"],
                           ev["band"], ev["run_first"])]
    with con:
        con.execute("DROP TABLE IF EXISTS clock_event_series")
        ev[["event_id", "series_id"]].to_sql("clock_event_series", con,
                                              index=False)
    gaia_cats = {n: gaia_field(n, tg.loc[n, "ra_deg"], tg.loc[n, "dec_deg"])
                 for n in sorted(set(ev["vsx_name"]))}
    series, prep_jobs = {}, []
    for sid, g in ev.groupby("series_id"):
        if sid in done:
            continue
        r0 = g.iloc[0]
        t = tg.loc[r0.vsx_name]
        frames = pd.read_sql_query("""
            SELECT path, jd_start, exptime_s FROM frames
            WHERE is_science = 1 AND target_key = ? AND band = ?
              AND camera = ? AND night = ? AND jd_start BETWEEN ? AND ?
            ORDER BY jd_start""", con, params=(
            r0.target_key, r0.band, r0.camera, r0.night,
            r0.run_first - 1e-6, r0.run_last + 1e-6))
        seed = int(hashlib.sha1(sid.encode()).hexdigest()[:6], 16)
        series[sid] = dict(r0=r0, t=t, frames=frames, seed=seed)
        prep_jobs.append(dict(sid=sid, archive=str(archive),
                              paths=list(frames["path"]), ra=t.ra_deg,
                              dec=t.dec_deg, seed=seed,
                              gaia={k: v.tolist() for k, v in
                                    gaia_cats[r0.vsx_name].items()}))
    log(f"photometry: {len(prep_jobs)} series to prepare "
        f"({len(done)} already done)")
    if not prep_jobs:
        return
    with mp.get_context("spawn").Pool(workers) as pool:
        preps = {p["sid"]: p for p in pool.imap_unordered(_prepare_job,
                                                          prep_jobs)}
        tasks = []
        for sid, p in preps.items():
            S = series[sid]
            if not p["ok"]:
                continue
            tol = max(ct.APERTURE_MIN_PX, 0.8 * p["fwhm_px"])
            tasks += [{"sid": sid, "archive": str(archive), "path": path,
                       "dark": None, "flat": None,
                       "ref_bright": p["ref_bright"], "seed": S["seed"],
                       "star_xy": p["star_xy"], "tol_px": tol,
                       "radii": p["radii"], "annulus": p["annulus"]}
                      for path in S["frames"]["path"]]
        log(f"photometry: {sum(p['ok'] for p in preps.values())} series "
            f"prepared; measuring {len(tasks):,} frames")
        results: dict = {}
        for i, r in enumerate(pool.imap_unordered(_measure_job, tasks,
                                                  chunksize=2), 1):
            results.setdefault(r["sid"], {})[r["path"]] = r
            if i % 1000 == 0:
                log(f"  … {i:,}/{len(tasks):,} frames")
    for sid, p in preps.items():
        S, r0, t = series[sid], series[sid]["r0"], series[sid]["t"]
        frames = S["frames"]
        if not p["ok"]:
            with con:
                con.execute("INSERT INTO clock_series (series_id, vsx_name, "
                            "night, camera, software, band, n_frames, status)"
                            " VALUES (?,?,?,?,?,?,?,?)",
                            (sid, r0.vsx_name, r0.night, r0.camera,
                             r0.software, r0.band, len(frames), p["reason"]))
            continue
        mid = frames["jd_start"] + frames["exptime_s"] / 2 / 86400.0
        bjd = bjd_tdb_from_utc(mid.to_numpy(), t.ra_deg, t.dec_deg)[0]
        rows = []
        for fr, b in zip(frames.itertuples(), np.atleast_1d(bjd)):
            r = results.get(sid, {}).get(fr.path)
            if r is None or r["status"] != "ok":
                continue
            for st in r["stars"]:
                rows.append((sid, fr.path, fr.jd_start, fr.exptime_s,
                             float(b), r.get("fwhm_px"), r.get("airmass"),
                             st["star"], st["x"], st["y"], st["matched"],
                             st["on_frame"], st["peak"], *st["flux"],
                             *st["err"]))
        with con:
            con.executemany("INSERT INTO clock_flux VALUES "
                            f"({','.join('?' * 19)})", rows)
            con.execute("INSERT INTO clock_series VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        (sid, r0.vsx_name, r0.night, r0.camera, r0.software,
                         r0.band, len(frames), "ok", p["ref_path"],
                         p["parity"], p["n_gaia"], p["target_sep"],
                         p["scale"], p["fwhm_px"], len(p["star_xy"]) - 1))
    log(f"photometry: wrote {len(preps)} series")


# ---------------------------------------------------------------------------
# Stage: fit
# ---------------------------------------------------------------------------
def series_arrays(con, sid: str) -> dict | None:
    """S3b's ``load_series_arrays`` contract, from ``clock_flux``."""
    df = pd.read_sql_query("SELECT * FROM clock_flux WHERE series_id = ? "
                           "ORDER BY bjd_mid_startconv, star", con,
                           params=(sid,))
    if df.empty:
        return None
    fr = df.drop_duplicates("path")[["path", "bjd_mid_startconv",
                                      "exptime_s", "fwhm_px", "airmass"]]
    n_star = int(df["star"].max()) + 1
    if len(fr) < ct.MIN_FIT_POINTS or n_star < 1 + ct.MIN_COMPS:
        return None
    index = {p: i for i, p in enumerate(fr["path"])}
    n = len(fr)
    flux = np.full((n, n_star, 3), np.nan)
    err = np.full((n, n_star, 3), np.nan)
    peak = np.full((n, n_star), np.nan)
    good = np.zeros((n, n_star), dtype=bool)
    txy = np.full((n, 2), np.nan)
    for r in df.itertuples():
        i = index[r.path]
        if r.star == 0:
            txy[i] = (r.x, r.y)
        flux[i, r.star] = [r.flux0, r.flux1, r.flux2]
        err[i, r.star] = [r.err0, r.err1, r.err2]
        peak[i, r.star] = np.nan if r.peak_adu is None else r.peak_adu
        good[i, r.star] = bool(r.on_frame) and (bool(r.matched)
                                                or r.star == 0)
    flux = np.where(pd.isna(flux), np.nan, flux).astype(float)
    err = np.where(pd.isna(err), np.nan, err).astype(float)
    return {"paths": list(fr["path"]),
            "bjd": fr["bjd_mid_startconv"].to_numpy(float),
            "exptime_s": fr["exptime_s"].to_numpy(float),
            "fwhm": pd.to_numeric(fr["fwhm_px"]).to_numpy(float),
            "airmass": pd.to_numeric(fr["airmass"]).to_numpy(float),
            "flux": flux, "err": err, "peak": peak, "good": good,
            "txy": txy}


def stage_fit(con, workers: int) -> None:
    ev = pd.read_sql_query("""
        SELECT e.*, s.series_id FROM clock_events e
        JOIN clock_event_series s USING (event_id)
        WHERE e.status = 'admitted'""", con)
    tg = pd.read_sql_query("SELECT * FROM clock_targets",
                           con).set_index("vsx_name")
    ok_series = {r[0] for r in con.execute(
        "SELECT series_id FROM clock_series WHERE status = 'ok'")}
    jobs, dead, meta = [], [], {}
    for sid, g in ev.groupby("series_id"):
        if sid not in ok_series:
            dead += [(e, "photometry refused") for e in g["event_id"]]
            continue
        arr = series_arrays(con, sid)
        if arr is None:
            dead += [(e, ct.STATUS_TOO_FEW) for e in g["event_id"]]
            continue
        lcs = s3b.build_lightcurves(arr, LEGACY_VETO_ADU, ONCHIP_BINNING)
        ad = lcs["adopted"]
        if ad is None:
            dead += [(e, ct.STATUS_NO_COMPS) for e in g["event_id"]]
            continue
        for ev_row in g.itertuples():
            t = tg.loc[ev_row.vsx_name]
            shape = shape_for(t.vsx_type, t.period_vsx)
            aps = []
            for ap in lcs["apertures"]:
                if ap is None:
                    aps.append(None)
                    continue
                win = ap["use"] & s3b.event_window("eclipse", arr["bjd"],
                                                   ev_row.t_pred, shape)
                aps.append({"t": arr["bjd"][win], "flux": ap["flux"][win],
                            "err": ap["err"][win],
                            "ex": arr["exptime_s"][win] / 86400.0,
                            "fwhm": arr["fwhm"][win],
                            "airmass": arr["airmass"][win],
                            "steps": ct.jump_steps(arr["txy"][win, 0],
                                                   arr["txy"][win, 1])})
            if aps[ad] is None or len(aps[ad]["t"]) < ct.MIN_FIT_POINTS:
                dead.append((ev_row.event_id, ct.STATUS_TOO_FEW))
                continue
            meta[ev_row.event_id] = float(np.median(aps[ad]["ex"]) * 86400)
            jobs.append({
                "event_id": ev_row.event_id, "apertures": aps, "adopted": ad,
                "shape": dict(period=shape.period, a_over_r=shape.a_over_r,
                              k=shape.k, b=shape.b, u1=shape.u1, u2=shape.u2,
                              f1=shape.f1, free=shape.free,
                              baseline_order=shape.baseline_order),
                "depth_expected": shape.nominal_depth(),
                "min_depth": MIN_DEPTH,
                "seed": int(hashlib.sha1(ev_row.event_id.encode()
                                         ).hexdigest()[:6], 16)})
    log(f"fit: {len(jobs)} eclipses to fit ({len(dead)} refused before)")
    if workers > 1 and len(jobs) > 1:
        with mp.get_context("spawn").Pool(workers) as pool:
            results = list(pool.imap_unordered(s3b.analyse_event, jobs))
    else:
        results = [s3b.analyse_event(j) for j in jobs]
    evi = ev.set_index("event_id")
    rows = [dict(event_id=e, status=why) for e, why in dead]
    for r in results:
        g = r.get
        sig = [v for v in (g("sig_boot_s"), g("sig_bead_s"))
               if v is not None and np.isfinite(v)]
        sig_stat = max(sig) if sig else None
        status = r["status"]
        if status == ct.STATUS_OK and sig_stat is None:
            status = "bootstrap_failed"
        e = evi.loc[r["event_id"]]
        oc_s = oc_m = None
        if status == ct.STATUS_OK:
            oc_s, oc_m = lk.oc_both_conventions(g("t0"), e.t_pred,
                                                meta[r["event_id"]])
        rows.append(dict(
            event_id=r["event_id"], status=status, t0_bjd=g("t0"),
            sig_stat_s=sig_stat, sig_model_s=g("sig_model_s"),
            depth=g("depth"), n_points=g("n_points"),
            n_before=g("n_before"), n_after=g("n_after"),
            chi2=g("chi2"), dof=g("dof"),
            exptime_s=meta.get(r["event_id"]),
            oc_start_s=oc_s, oc_mid_s=oc_m))
    replace(con, "clock_fits", pd.DataFrame(rows))
    n_ok = sum(r["status"] == ct.STATUS_OK for r in rows)
    log(f"fit: {n_ok} of {len(rows)} eclipses timed")


# ---------------------------------------------------------------------------
# Stage: summary
# ---------------------------------------------------------------------------
def stage_summary(con) -> None:
    """O − C per season, every season of the archive listed.

    The ADOPTED reading for each camera epoch is the census verdict
    (``time_audit.convention``); the other reading is reported beside it,
    and the eclipses themselves say which one is right.  An event's error
    is its statistical (bootstrap/prayer-bead) and model (variant) terms
    plus the prediction's formal and systematic terms, in quadrature.
    """
    fr = pd.read_sql_query("""
        SELECT camera, software, night FROM frames WHERE is_science = 1""",
                           con)
    fr["season"] = [lk.season_label(n) for n in fr["night"]]
    seasons = fr.groupby(["camera", "software", "season"]).agg(
        n_frames=("night", "size"), n_nights=("night", "nunique"),
        first_night=("night", "min"), last_night=("night", "max")
    ).reset_index()
    conv = dict(((c, s), v) for c, s, v in con.execute(
        "SELECT camera, software, convention FROM time_audit"))
    ev = pd.read_sql_query("""
        SELECT e.event_id, e.vsx_name, e.night, e.camera, e.software,
               e.sig_pred_s, e.sig_sys_s, e.ref_range,
               f.status, f.sig_stat_s, f.sig_model_s, f.oc_start_s,
               f.oc_mid_s, f.exptime_s
        FROM clock_events e JOIN clock_fits f USING (event_id)""", con)
    ev["season"] = [lk.season_label(n) for n in ev["night"]]
    ok = ev[ev["status"] == ct.STATUS_OK].copy()
    ok["sigma_s"] = np.sqrt(ok["sig_stat_s"] ** 2
                            + ok["sig_model_s"].fillna(0) ** 2
                            + ok["sig_pred_s"].fillna(0) ** 2
                            + ok["sig_sys_s"].fillna(0) ** 2)
    rows = []
    for s in seasons.itertuples():
        g = ok[(ok["camera"] == s.camera) & (ok["software"] == s.software)
               & (ok["season"] == s.season)]
        c = conv.get((s.camera, s.software), "")
        mid = c.startswith("MIDDLE")
        base = dict(camera=s.camera, software=s.software, season=s.season,
                    first_night=s.first_night, last_night=s.last_night,
                    n_frames=s.n_frames, n_nights=s.n_nights,
                    convention=c, n_events=len(g),
                    n_events_tried=int(((ev["camera"] == s.camera)
                                        & (ev["software"] == s.software)
                                        & (ev["season"] == s.season)).sum()))
        if g.empty:
            rows.append({**base, "verdict": "NO CLOCK EVENT — offset not "
                         "measured; carried as unknown"})
            continue
        adopted = g["oc_mid_s"] if mid else g["oc_start_s"]
        other = g["oc_start_s"] if mid else g["oc_mid_s"]
        cm = ct.combine_oc(adopted, g["sigma_s"])
        co = ct.combine_oc(other, g["sigma_s"])
        rows.append({**base,
                     "oc_s": cm["wmean"], "sigma_s": cm["sigma_adopted"],
                     "chi2": cm["chi2"], "dof": cm["dof"],
                     "scatter_err_s": cm["scatter_err"],
                     "oc_other_reading_s": co["wmean"],
                     "systems": ", ".join(sorted(set(g["vsx_name"]))),
                     "verdict": lk.season_verdict(cm["wmean"],
                                                  cm["sigma_adopted"])})
    out = pd.DataFrame(rows).sort_values(["first_night", "camera"])
    replace(con, "clock_seasons", out)
    with con:
        con.execute("CREATE TABLE IF NOT EXISTS clock_meta "
                    "(key TEXT PRIMARY KEY, value TEXT)")
        con.executemany("INSERT OR REPLACE INTO clock_meta VALUES (?,?)", [
            ("built_utc", datetime.now(timezone.utc).isoformat(
                timespec="seconds")),
            ("accept_s", str(lk.CLOCK_ACCEPT_S)),
            ("extrap_max_d", str(lk.EXTRAP_MAX_D)),
            ("veto_adu", str(LEGACY_VETO_ADU))])
    n_meas = int(out["oc_s"].notna().sum()) if "oc_s" in out else 0
    log(f"summary: {n_meas} of {len(out)} seasons measured")


STAGES = ("tess", "events", "photometry", "fit", "summary")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--stage", choices=STAGES + ("all",), default="all")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args(argv)
    warnings.filterwarnings("ignore")
    con = sqlite3.connect(args.db, timeout=120)
    con.execute("PRAGMA journal_mode=WAL")
    try:
        todo = STAGES if args.stage == "all" else (args.stage,)
        for st in todo:
            if st == "tess":
                stage_tess(con)
            elif st == "events":
                stage_events(con)
            elif st == "photometry":
                stage_photometry(con, args.archive, args.workers)
            elif st == "fit":
                stage_fit(con, args.workers)
            elif st == "summary":
                stage_summary(con)
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
