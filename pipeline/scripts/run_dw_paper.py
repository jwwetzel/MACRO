#!/usr/bin/env python
"""run_dw_paper.py — the Dwarf-Galaxy Hα paper, frames to manuscript.

One resumable CLI for every remaining DwarfGalaxy_AGN_Survey task
(SYNTHESIS §4 DwarfGalaxy_AGN_Survey, U1; strategy §10).  Pure rules live in
``macro_dw.dwcore`` (unit-tested in ``pipeline/tests/test_dw.py``); pixels in
``macro_dw.dwio``; figures and numbers in ``macro_dw.paper_dw``.  Everything
lands in ONE database,

    products/dwarf/dwarf.sqlite

and the manifest is opened READ-ONLY.  The archive is never written.

SUBCOMMANDS, in dependency order (each rebuilds its own tables)
---------------------------------------------------------------
    fetch      REFCAT2 cones per field (VizieR J/ApJ/867/105), PS1 DR1 cone
               for the zero-point cross-check (VizieR II/349), cached with
               retrieval dates under products/dwarf/external/
    frames     dw_frames: every staged science frame with pointing columns,
               night label and BJD_TDB                      [DW-P01, DW-P46]
    flatprep   per-frame source masks and 64-px sky block maps (parallel)
    flats      superflats + held-out residual test, fringe and moonlight
               gradient tests                               [DW-P11, DW-P12]
    measure    calibrate, plate-solve on REFCAT2, multi-aperture forced
               photometry per frame (parallel, resumable)          [DW-P03]
    zp         ensemble zero points + colour terms per filter per readout
               mode, PS1 check, growth curves, ZMAG QC flag
                                              [DW-P2-ensemble-zp, P21, P37]
    qc         the disposition table (the only rejection authority)
                                                           [DW-P32, DW-P01]
    stacks     weighted coadds per field per filter, constant+plane sky
                                                           [DW-P33, DW-P51]
    sersic     structure of each object on its L stack            [DW-P35]
    depth      Román limits + synthetic-dwarf recovery           [DW-P34]
    halpha     NGC 5238 line-flux scale; Hα net flux / limits per field
                                                                  [DW-P36]
    lightcurves NGC 5238 field-star series, covariates, Sys-Rem   [DW-P55]
    completeness injection-recovery map, signed cell bias        [DW-P54]
    zeroorder  NGC 5548 zero-order differential photometry       [DW-P4y]
    paper      figures + numbers.tex                  [DW-figures, DW-draft]
    status     progress (read-only)

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/run_dw_paper.py fetch
    $PY pipeline/scripts/run_dw_paper.py frames
    $PY pipeline/scripts/run_dw_paper.py flatprep --workers 10
    ...
"""
from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import subprocess
import sys
import time
import urllib.request
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_dw import dwcore as core                  # noqa: E402

warnings.filterwarnings("ignore")

REPO = PIPELINE_ROOT.parent
MANIFEST = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
OUT = REPO / "products" / "dwarf"
DB = OUT / "dwarf.sqlite"
EXT = OUT / "external"
NOVELTY = REPO / "DwarfGalaxy_AGN_Survey" / "notes" / "novelty"

from macro_dw import DW_CODE_VERSION                 # noqa: E402

#: Galaxy positions not in the novelty table (NED, J2000).
NGC5238_RADEC = (203.677125, 51.613611)        # 13:34:42.51 +51:36:49.0
NGC5548_RADEC = (214.498042, 25.136778)        # 14:17:59.53 +25:08:12.4
#: Published integrated Hα flux of NGC 5238 (LVGDB, Kaisin & Karachentsev
#: BTA imaging; novelty package sources/lvgdb): log F = -12.26 -> 5.50e-13,
#: with the ±0.81e-13 error quoted there.
NGC5238_FHA = 5.50e-13
NGC5238_FHA_ERR = 0.81e-13
NGC5238_CZ = 229.0                              # km/s (LVGDB)

#: Bands that enter the work (W is out of NGC 5238 surface photometry by
#: ruling, and out of everything else for want of a flat it could trust).
DW_BANDS = ("L", "H", "R")
N5238_BANDS = ("G", "H", "L", "O", "R", "X")


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def connect(read_only: bool = False) -> sqlite3.Connection:
    OUT.mkdir(parents=True, exist_ok=True)
    uri = f"file:{DB}?mode=ro" if read_only else f"file:{DB}"
    con = sqlite3.connect(uri, uri=True, timeout=300)
    con.row_factory = sqlite3.Row
    return con


def manifest() -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{MANIFEST}?mode=ro", uri=True, timeout=300)
    con.row_factory = sqlite3.Row
    return con


def meta(con, **kw) -> None:
    con.execute("CREATE TABLE IF NOT EXISTS dw_build_meta "
                "(key TEXT PRIMARY KEY, value TEXT)")
    for k, v in kw.items():
        con.execute("INSERT OR REPLACE INTO dw_build_meta VALUES (?, ?)",
                    (k, str(v)))
    con.commit()


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def write_table(con, name: str, rows: list[dict]) -> None:
    """Replace table ``name`` with ``rows`` (columns from the first row)."""
    con.execute(f"DROP TABLE IF EXISTS {name}")
    if not rows:
        con.commit()
        return
    cols = list(rows[0].keys())
    con.execute(f"CREATE TABLE {name} ({', '.join(cols)})")
    con.executemany(f"INSERT INTO {name} VALUES ({','.join('?' * len(cols))})",
                    [tuple(_py(r.get(c)) for c in cols) for r in rows])
    con.commit()


def _py(v):
    if isinstance(v, (np.floating,)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, float) and not np.isfinite(v):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return int(v)
    return v


def literature() -> dict:
    """The novelty package's per-object rows, keyed by field name."""
    with open(NOVELTY / "out" / "novelty_table.csv") as fh:
        return {r["field"]: r for r in csv.DictReader(fh)}


def sexa(ra: str, dec: str) -> tuple[float, float]:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    return float(c.ra.deg), float(c.dec.deg)


def object_radec(target: str) -> tuple[float, float]:
    if target == "NGC5238":
        return NGC5238_RADEC
    if target == "NGC 5548":
        return NGC5548_RADEC
    r = literature()[target]
    return sexa(r["ra_hms"], r["dec_dms"])


# ===========================================================================
# fetch
# ===========================================================================
REFCAT_FMT = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"
              "-source=J/ApJ/867/105/refcat2&-c={ra:.5f}+{dec:+.5f}"
              "&-c.r=30&-c.u=arcmin&-out.max=100000"
              "&-out=RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,pmDE,Gmag,"
              "gmag,e_gmag,rmag,e_rmag,imag,e_imag,rp1,r1,dupvar&rmag=<19")
PS1_FMT = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"
           "-source=II/349/ps1&-c={ra:.5f}+{dec:+.5f}&-c.r=20&-c.u=arcmin"
           "&-out.max=100000&-out=RAJ2000,DEJ2000,gmag,e_gmag,rmag,e_rmag,"
           "imag,e_imag,Nr&rmag=<18")
#: Fields for the PS1 cross-check: the Dw field with most R frames and
#: NGC 5238 (REFCAT2 is itself PS1-based at these declinations, so the
#: check tests the REFCAT2 merge, not an independent system — stated).
PS1_FIELDS = ("Dw1403+49", "NGC5238")


def field_centres() -> dict:
    """Median header pointing per staged target (deg)."""
    m = manifest()
    rows = m.execute(
        "SELECT s.canonical_target t, f.ra_deg, f.dec_deg FROM "
        "stage_dwarfgalaxy_agn_survey s JOIN frames f USING(obs_rowid) "
        "WHERE s.role='science' AND f.ra_deg IS NOT NULL").fetchall()
    out = {}
    for t in sorted({r["t"] for r in rows}):
        ra = np.array([r["ra_deg"] for r in rows if r["t"] == t])
        de = np.array([r["dec_deg"] for r in rows if r["t"] == t])
        out[t] = (float(np.median(ra)), float(np.median(de)))
    return out


def _get(url: str, dest: Path) -> int:
    req = urllib.request.Request(url, headers={"User-Agent": "MACRO-RLMT/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r:
        data = r.read()
    dest.write_bytes(data)
    return len(data)


def cmd_fetch(args) -> int:
    EXT.mkdir(parents=True, exist_ok=True)
    log = EXT / "retrieval.json"
    got = json.loads(log.read_text()) if log.exists() else {}
    for t, (ra, de) in field_centres().items():
        key = t.replace(" ", "")
        jobs = [(f"refcat2_{key}", REFCAT_FMT.format(ra=ra, dec=de))]
        if key in PS1_FIELDS:
            jobs.append((f"ps1_{key}", PS1_FMT.format(ra=ra, dec=de)))
        for name, url in jobs:
            dest = EXT / f"{name}.tsv"
            if dest.exists() and not args.refresh:
                continue
            n = _get(url, dest)
            got[name] = {"url": url, "bytes": n, "retrieved": utcnow()}
            print(f"fetched {name} ({n} bytes)")
            log.write_text(json.dumps(got, indent=1))
    return 0


def load_tsv(path: Path) -> dict:
    """VizieR asu-tsv -> numpy columns (blank -> NaN)."""
    lines = [ln for ln in path.read_text().splitlines()
             if ln and not ln.startswith("#")]
    head = lines[0].split("\t")
    body = [ln.split("\t") for ln in lines[3:]]        # name, unit, dashes
    out = {}
    for j, h in enumerate(head):
        col = []
        for b in body:
            v = b[j].strip() if j < len(b) else ""
            try:
                col.append(float(v))
            except ValueError:
                col.append(np.nan)
        out[h.strip()] = np.array(col)
    return out


def refcat(target: str) -> dict:
    """REFCAT2 cone of a field, renamed to ra/dec/g/r/i (+ errors)."""
    c = load_tsv(EXT / f"refcat2_{target.replace(' ', '')}.tsv")
    ok = np.isfinite(c["rmag"]) & np.isfinite(c["gmag"])
    return {"ra": c["RA_ICRS"][ok], "dec": c["DE_ICRS"][ok],
            "g": c["gmag"][ok], "r": c["rmag"][ok], "i": c["imag"][ok],
            "er": c["e_rmag"][ok], "rp1": c["rp1"][ok], "dupvar": c["dupvar"][ok],
            "G": c["Gmag"][ok]}


# ===========================================================================
# frames  (DW-P01 pointing columns, DW-P46 BJD_TDB)
# ===========================================================================
#: Time-stamp convention (clock package, products/clock/clock_transits.sqlite
#: table s3b_stamp_convention, cadence test): AC4040 StackPro frames stamp
#: DATE-OBS (= header JD) at MID-exposure; plain High Gain frames of era 2
#: (MaxIm 6.30) are undetermined and are taken as start-stamped, so their
#: mid-times carry a +-EXPTIME/2 caveat (<= 256 s for 512 s Hα frames).
def mid_offset_s(readoutm: str, exptime: float) -> float:
    return 0.0 if "StackPro" in (readoutm or "") else 0.5 * float(exptime)


def bjd_tdb(jd_stamp: np.ndarray, offset_s: np.ndarray, ra: float,
            dec: float) -> np.ndarray:
    """Barycentric TDB Julian date at MID-exposure for a fixed target:
    the header stamp plus the convention's offset (``mid_offset_s``)."""
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation
    import astropy.units as u
    site = EarthLocation(lat=31.66556 * u.deg, lon=-110.60194 * u.deg,
                         height=1515 * u.m)          # Winer (header SITELAT/LONG)
    t = Time(jd_stamp + offset_s / 86400.0, format="jd", scale="utc",
             location=site)
    ltt = t.light_travel_time(SkyCoord(ra, dec, unit="deg"), "barycentric")
    return (t.tdb + ltt).jd


def cmd_frames(args) -> int:
    m = manifest()
    rows = [dict(r) for r in m.execute(
        "SELECT s.obs_rowid, s.path, s.canonical_target AS target, s.filter, "
        "s.night AS night_manifest, s.jd, s.exptime, s.era_id, s.mech_epoch, "
        "f.is_canonical, f.dup_group, f.dup_basis, f.readoutm, f.camera, "
        "f.ra_deg, f.dec_deg, f.airmass, f.fwhm AS fwhm_hdr, f.zmag, "
        "f.moonangl, f.moonphas, f.egain, f.pltsolvd, f.qc_flags, "
        "f.pointing_offset_deg AS offset_manifest_deg "
        "FROM stage_dwarfgalaxy_agn_survey s JOIN frames f USING(obs_rowid) "
        "WHERE s.role='science' ORDER BY s.canonical_target, s.jd")]
    centres = field_centres()
    for t in sorted({r["target"] for r in rows}):
        sel = [r for r in rows if r["target"] == t]
        ra0, de0 = object_radec(t)
        jd = np.array([r["jd"] for r in sel])
        off = np.array([mid_offset_s(r["readoutm"], r["exptime"]) for r in sel])
        bj = bjd_tdb(jd, off, ra0, de0)
        c_ra, c_de = centres[t]
        for r, b in zip(sel, bj):
            r["bjd_tdb_mid"] = float(b)
            r["time_convention"] = ("StackPro: stamp = mid-exposure (s3b)"
                                    if "StackPro" in (r["readoutm"] or "") else
                                    "High Gain era 2: start-stamp assumed (undetermined, +-EXPTIME/2)")
            r["night"] = core.night_label(r["jd"])
            if r["ra_deg"] is None:
                r["offset_field_deg"] = None
            else:
                d = np.hypot((r["ra_deg"] - c_ra) * np.cos(np.radians(c_de)),
                             r["dec_deg"] - c_de)
                r["offset_field_deg"] = float(d)
            r["obj_ra"], r["obj_dec"] = ra0, de0
            r["spectral"] = int(r["filter"] in core.SPECTRAL_FILTERS)
    con = connect()
    write_table(con, "dw_frames", rows)
    meta(con, frames_utc=utcnow(), code=DW_CODE_VERSION, git=git_commit())
    n_mis = sum(1 for r in rows if (r["offset_field_deg"] or 0) > core.MISPOINT_DEG)
    print(f"dw_frames: {len(rows)} frames, {n_mis} beyond "
          f"{core.MISPOINT_DEG} deg of their field centre")
    return 0


def worklist(con, imaging_only: bool = True, cameras=("AC4040",)) -> list[dict]:
    """Direct-imaging, canonical, on-field frames of the bands that enter."""
    rows = [dict(r) for r in con.execute(
        "SELECT * FROM dw_frames WHERE spectral=0 AND is_canonical=1")]
    out = []
    for r in rows:
        if (r["offset_field_deg"] or 0) > core.MISPOINT_DEG:
            continue
        if r["camera"] not in cameras and not (r["target"] == "NGC 5548"):
            continue
        if r["target"].startswith("Dw") and r["filter"] not in DW_BANDS:
            continue
        if r["target"] == "NGC5238" and r["filter"] not in N5238_BANDS + ("W",):
            continue
        out.append(r)
    return out


# ===========================================================================
# flatprep  (per-frame masks + block maps; DW-P11 input)
# ===========================================================================
def _flatprep_one(task: dict) -> dict:
    from macro_dw import dwio
    from macro_sn import snio
    out = {"obs_rowid": task["obs_rowid"], "ok": 0}
    try:
        raw, h = dwio.load(task["path"])
        img = dwio.dark_subtract(raw, task["readoutm"], task["exptime"])
        bp = snio.badpix_mask({"camera": "AC4040", "readoutm": task["readoutm"],
                               "flipstat": h.get("FLIPSTAT")}, raw.shape)
        t, dead = dwio.twilight(task["filter"] if task["filter"] in dwio.TWILIGHT
                                else "L")
        base = dead.copy()
        if bp is not None:
            base |= bp
        # Detect on the twilight-flattened frame so vignetting does not
        # bias the 2-sigma footprint toward the centre.
        mask = dwio.source_mask(img / t, base)
        blk = dwio.block_median(img, mask)
        np.savez_compressed(dwio.BLOCKDIR / f"{task['obs_rowid']}.npz",
                            blocks=blk.astype(np.float32))
        np.save(dwio.MASKDIR / f"{task['obs_rowid']}.npy", np.packbits(mask))
        out.update(ok=1, mask_frac=float(mask.mean()),
                   sky_med=float(np.nanmedian(blk)),
                   badpix=int(bp is not None))
    except Exception as e:                                   # pragma: no cover
        out["error"] = repr(e)[:300]
    return out


def cmd_flatprep(args) -> int:
    from multiprocessing import Pool
    from macro_dw import dwio
    dwio.BLOCKDIR.mkdir(parents=True, exist_ok=True)
    dwio.MASKDIR.mkdir(parents=True, exist_ok=True)
    con = connect()
    tasks = [r for r in worklist(con) if r["camera"] == "AC4040"]
    done = {p.stem for p in dwio.BLOCKDIR.glob("*.npz")}
    todo = [t for t in tasks if str(t["obs_rowid"]) not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"flatprep: {len(tasks)} frames, {len(todo)} to do")
    res = []
    with Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(_flatprep_one, todo, 2)):
            res.append(r)
            if k % 50 == 0:
                print(f"  {k}/{len(todo)}", flush=True)
    con.execute("CREATE TABLE IF NOT EXISTS dw_flatprep (obs_rowid INTEGER "
                "PRIMARY KEY, ok INTEGER, mask_frac REAL, sky_med REAL, "
                "badpix INTEGER, error TEXT)")
    con.executemany("INSERT OR REPLACE INTO dw_flatprep VALUES (?,?,?,?,?,?)",
                    [(r["obs_rowid"], r["ok"], r.get("mask_frac"),
                      r.get("sky_med"), r.get("badpix"), r.get("error"))
                     for r in res])
    con.commit()
    print(f"flatprep: {sum(r['ok'] for r in res)} ok, "
          f"{sum(1 - r['ok'] for r in res)} failed")
    return 0


# ===========================================================================
# flats  (DW-P11) + fringe and moonlight-gradient tests (DW-P12)
# ===========================================================================
#: The ruled L superflat is built from the June 2023 Dw L frames only;
#: NGC 5238's February-May L frames are the held-out epoch test.
L_FLAT_NIGHTS = ("2023-06-01", "2023-06-30")
#: Pixels covered by fewer frames than this fall back to twilight x sky.
L_MIN_COVER = 20
#: Pixel-level outlier rejection in the full-resolution pass.
CR_SIGMA = 6.0
#: Frames whose fitted sky is below this (ADU) carry no flat information
#: (the 32 s O frames); they are calibrated with the flat, never build it.
FLAT_MIN_SKY_ADU = 20.0
DETLIM: dict = {}
#: Every Dw field put its galaxy at the same chip position (~(1024, 1024))
#: and NGC 5238 mostly there too, so a superflat that only masks 2-sigma
#: footprints absorbs the galaxies' sub-threshold outskirts at that spot
#: and later carves a hole with a bright rim around every target (seen in
#: the first R stacks).  Each frame's own target is therefore excluded
#: inside this radius (arcsec) from every flat product.
FLAT_EXCLUDE_ARCSEC = {"NGC5238": 240.0, "dw": 120.0}
#: Bands bright enough in night sky for a FULL-RESOLUTION superflat (the
#: per-pixel noise of the frame mean, sigma_pix / sqrt(N), is below 1 %:
#: L 870 ADU sky x 323 frames, R 440 x 160, G 345 x 83, X 315 x 65,
#: W 580 x 46).  H (35 ADU) and O (9 ADU) are not: there the twilight
#: master supplies the pixel scales and the night sky only the large ones.
SUPERFLAT_BANDS = ("L", "R", "G", "X", "W")


def _smooth_blocks(m: np.ndarray) -> np.ndarray:
    from scipy.ndimage import median_filter, gaussian_filter
    from macro_dw import dwio
    f = dwio.fill_nan_blocks(m)
    f = median_filter(f, size=3, mode="nearest")
    return gaussian_filter(f, 1.0, mode="nearest")


def _fullres_chunk(job: dict) -> dict:
    """Masked, outlier-clipped per-half sums of (dark-subtracted, pedestal-
    corrected, sky-normalised) frames, divided by ``job['divide']`` flat
    (None for the L superflat itself)."""
    from macro_dw import dwio
    shape = (4096, 4096)
    acc = {k: np.zeros(shape, np.float32) for k in ("sA", "nA", "sB", "nB")}
    model = np.load(job["model"]) if job.get("model") else None
    div = np.load(job["divide"]) if job.get("divide") else None
    for fr in job["frames"]:
        raw, _ = dwio.load(fr["path"])
        img = dwio.dark_subtract(raw, fr["readoutm"], fr["exptime"])
        n = (img - fr["c"]) / fr["s0"]
        mask = np.unpackbits(np.load(dwio.MASKDIR / f"{fr['obs_rowid']}.npy")
                             )[:n.size].reshape(shape).astype(bool)
        yy, xx = np.ogrid[:shape[0], :shape[1]]
        mask |= (xx - fr["x_obj"]) ** 2 + (yy - fr["y_obj"]) ** 2 < fr["excl_px"] ** 2
        if model is not None:
            res = n - model
            s = core.robust_sigma(res[~mask][::101])
            mask |= np.abs(res) > CR_SIGMA * s
        if div is not None:
            n = n / div
        n[mask] = 0.0
        h = "A" if fr["half"] == 0 else "B"
        acc["s" + h] += n
        acc["n" + h] += (~mask)
    return acc


def _fringe_stat(R: np.ndarray, D: np.ndarray, good: np.ndarray) -> dict:
    """Pixel-scale (4-64 px) structure of a residual map R about 1, with its
    noise measured from the half-difference D: amplitude = sqrt(var(R_hp) -
    var(D_hp)/4) on 4x4-binned, 64-px-median-high-passed maps."""
    from scipy.ndimage import median_filter

    def hp(a):
        b = a.reshape(1024, 4, 1024, 4).mean(axis=(1, 3))
        return b - median_filter(b, size=16, mode="nearest")
    g = good.reshape(1024, 4, 1024, 4).all(axis=(1, 3))
    g[:64], g[-64:], g[:, :64], g[:, -64:] = False, False, False, False
    r, d = hp(np.where(good, R, 1.0)), hp(np.where(good, D, 0.0))
    sr, sd = core.robust_sigma(r[g]), core.robust_sigma(d[g])
    return {"hp_rms": sr, "hp_noise": sd / 2,
            "fringe_amp": float(np.sqrt(max(sr ** 2 - (sd / 2) ** 2, 0.0)))}


def cmd_flats(args) -> int:
    from multiprocessing import Pool
    from astropy.io import fits
    from scipy.stats import spearmanr
    from macro_dw import dwio
    dwio.FLATDIR.mkdir(parents=True, exist_ok=True)
    tmp = OUT / "tmp"
    tmp.mkdir(exist_ok=True)
    con = connect()
    frames = [dict(r) for r in con.execute(
        "SELECT d.*, s.x_obj, s.y_obj FROM dw_frames d JOIN dw_flatprep p "
        "USING(obs_rowid) LEFT JOIN dw_solve s USING(obs_rowid) WHERE p.ok=1")]
    # Unsolved frames take their target's median object position.
    for t in {f["target"] for f in frames}:
        xs = [f["x_obj"] for f in frames if f["target"] == t and f["x_obj"] is not None]
        ys = [f["y_obj"] for f in frames if f["target"] == t and f["y_obj"] is not None]
        for f in frames:
            if f["target"] == t and f["x_obj"] is None:
                f["x_obj"] = float(np.median(xs)) if xs else 1024.0
                f["y_obj"] = float(np.median(ys)) if ys else 1024.0
    for f in frames:
        f["excl_px"] = FLAT_EXCLUDE_ARCSEC["NGC5238" if f["target"] == "NGC5238"
                                           else "dw"] / core.PIX_ARCSEC_NOMINAL
    bc = (np.arange(4096 // core.FLAT_BLOCK) + 0.5) * core.FLAT_BLOCK
    rows_fit, rows_val, rows_grad = [], [], []
    global DETLIM
    DETLIM = detector_limits()
    S_L = None
    for filt in ("L", "H", "R", "G", "O", "X", "W"):
        fr = [f for f in frames if f["filter"] == filt]
        if not fr:
            continue
        T, dead = dwio.twilight(filt)
        Tb = dwio.block_median(T, dead)
        # Template the frame's sky is fitted against.  L: the twilight shape
        # with NO additive term — at ~2500 ADU of sky a pedestal mismatch of
        # ~10-20 ADU is < 1 %, and an additive term fitted on a template is
        # degenerate with the very large-scale correction being measured.
        # Other bands: twilight x the L night-sky correction (illumination is
        # common to the optics), with a per-frame additive term, which on a
        # FIXED template is identifiable and absorbs the dark masters'
        # pedestal drift (18 ADU between the 256 s and 512 s High Gain
        # masters) — a drift comparable to the 40 ADU Hα sky.
        Tt = Tb if filt == "L" else Tb * S_L
        R, keep = [], []
        for f in fr:
            b = np.load(dwio.BLOCKDIR / f"{f['obs_rowid']}.npz")["blocks"]
            b = np.where(np.hypot(bc[None, :] - f["x_obj"], bc[:, None] - f["y_obj"])
                         < f["excl_px"] + core.FLAT_BLOCK, np.nan, b)
            if filt == "L":
                s0, c = float(np.nanmedian(b / Tt)), 0.0
            else:
                s0, c = core.fit_sky_on_template(b, Tt)
            lim = DETLIM.get(f["readoutm"])
            ok = (np.isfinite(s0) and s0 >= FLAT_MIN_SKY_ADU
                  and (lim is None or np.nanmedian(b) + lim["bias"] <= lim["cap"]))
            f["s0"], f["c"] = s0, c
            R.append((b - c) / s0 / Tt if ok else np.full(b.shape, np.nan))
            keep.append(ok)
            rows_fit.append({"obs_rowid": f["obs_rowid"], "filter": filt,
                             "sky_s0": s0, "pedestal_c": c, "used": int(ok),
                             "c_over_sky": c / s0 if ok else None})
        keep = np.array(keep)
        fr = [f for f, k in zip(fr, keep) if k]
        R = np.array(R)[keep]
        nights = sorted({f["night"] for f in fr})
        half = np.array([nights.index(f["night"]) % 2 for f in fr])
        for f, h in zip(fr, half):
            f["half"] = int(h)
        if filt == "L":
            use = np.array([f["target"].startswith("Dw")
                            and L_FLAT_NIGHTS[0] <= f["night"] <= L_FLAT_NIGHTS[1]
                            for f in fr])
        else:
            use = np.ones(len(fr), bool)
        S = _smooth_blocks(np.nanmedian(R[use], axis=0))
        if filt == "L":
            S_L = S
        # -- held-out residual test (DW-P11 criterion) ----------------------
        for name, a, b in (("half_B_on_A", use & (half == 0), use & (half == 1)),
                           ("half_A_on_B", use & (half == 1), use & (half == 0))):
            if a.sum() < 3 or b.sum() < 3:
                continue
            Sa = _smooth_blocks(np.nanmedian(R[a], axis=0))
            resid = np.nanmedian(R[b] / Sa, axis=0)
            st = core.residual_structure(resid)
            rows_val.append({"filter": filt, "test": name, "n_build": int(a.sum()),
                             "n_test": int(b.sum()), **st})
        if filt == "L" and (~use).sum() >= 3:
            resid = np.nanmedian(R[~use] / S, axis=0)
            rows_val.append({"filter": filt, "test": "NGC5238_FebMay_on_June",
                             "n_build": int(use.sum()), "n_test": int((~use).sum()),
                             **core.residual_structure(resid)})
        dw = np.array([f["target"].startswith("Dw") for f in fr])
        if filt in ("H", "R") and dw.sum() >= 3 and (~dw).sum() >= 3:
            Sd = _smooth_blocks(np.nanmedian(R[dw], axis=0))  # noqa
            Sn = _smooth_blocks(np.nanmedian(R[~dw], axis=0))
            rows_val.append({"filter": filt, "test": "epoch_June_vs_FebMay",
                             "n_build": int(dw.sum()), "n_test": int((~dw).sum()),
                             **core.residual_structure(Sd / Sn)})
        # -- gradients vs moon (DW-P12) -------------------------------------
        for f, r in zip(fr, R):
            res, coef = core.remove_plane(r / S)
            g = float(np.hypot(coef[1], coef[2]) * 64)   # across 64 blocks
            rows_grad.append({"obs_rowid": f["obs_rowid"], "filter": filt,
                              "grad_frac_chip": g,
                              "grad_pa_deg": float(np.degrees(np.arctan2(coef[2], coef[1]))),
                              "resid_rms_after_plane": core.robust_sigma(res),
                              "moonangl": f["moonangl"], "moonphas": f["moonphas"]})
        # -- full-resolution pass -------------------------------------------
        Sup = dwio.upsample_blocks(S if filt == "L" else S * S_L)
        model = (T * Sup).astype(np.float32)
        np.save(tmp / f"model_{filt}.npy", model)
        sel = [f for f, u in zip(fr, use) if u]
        superflat = filt in SUPERFLAT_BANDS
        if not superflat:
            np.save(tmp / f"div_{filt}.npy", model)
        k = max(1, min(args.workers, len(sel) // 4))
        jobs = [{"frames": [dict(path=f["path"], readoutm=f["readoutm"],
                                 exptime=f["exptime"], c=f["c"], s0=f["s0"],
                                 obs_rowid=f["obs_rowid"], half=f["half"],
                                 x_obj=f["x_obj"], y_obj=f["y_obj"],
                                 excl_px=f["excl_px"])
                            for f in sel[i::k]],
                 "model": str(tmp / f"model_{filt}.npy"),
                 "divide": None if superflat else str(tmp / f"div_{filt}.npy")}
                for i in range(k)]
        with Pool(k) as pool:
            parts = pool.map(_fullres_chunk, jobs)
        acc = {key: sum(p[key] for p in parts) for key in parts[0]}
        nA, nB = acc["nA"], acc["nB"]
        mA = np.where(nA > 0, acc["sA"] / np.maximum(nA, 1), np.nan)
        mB = np.where(nB > 0, acc["sB"] / np.maximum(nB, 1), np.nan)
        ntot = nA + nB
        mean = np.where(ntot > 0, (acc["sA"] + acc["sB"]) / np.maximum(ntot, 1),
                        np.nan)
        if superflat:
            sf = mean / np.nanmedian(mean[1024:3072, 1024:3072])
            fallback = (ntot < (L_MIN_COVER if filt == "L" else 10)) \
                | ~np.isfinite(sf)
            flat = np.where(fallback, model, sf)
            Rm, Dm = sf / model, (mA - mB) / model
            good = ~fallback & ~dead & (nA > 5) & (nB > 5)
            extra = {"fallback_frac": float(fallback.mean()),
                     "median_cover": float(np.median(ntot))}
        else:
            flat = model
            Rm, Dm = mean, mA - mB
            good = ~dead & (nA > 5) & (nB > 5) & np.isfinite(Rm)
            extra = {}
        fs = _fringe_stat(np.nan_to_num(Rm, nan=1.0), np.nan_to_num(Dm), good)
        rows_val.append({"filter": filt, "test": "fringe_pixel_scale",
                         "n_build": len(sel), "n_test": len(sel),
                         "rms": fs["fringe_amp"], "rms_noplane": fs["hp_rms"],
                         "p95_p05": fs["hp_noise"], **extra})
        flat = np.where(dead, np.nan, flat).astype(np.float32)
        hdr = fits.Header()
        hdr["FILTER"] = filt
        hdr["RECIPE"] = ("night-sky superflat (June 2023 Dw L frames)"
                         if filt == "L" else
                         "night-sky superflat (all era-2 frames of the band)"
                         if superflat else
                         "twilight 2023-06-30 x night-sky large-scale")
        hdr["NFRAMES"] = len(sel)
        fits.writeto(dwio.FLATDIR / f"flat_{filt}.fits", flat, hdr,
                     overwrite=True)
        print(f"flat {filt}: {len(sel)} frames", flush=True)
    write_table(con, "dw_flat_fits", rows_fit)
    write_table(con, "dw_flat_validation", rows_val)
    write_table(con, "dw_gradients", rows_grad)
    # Moon correlation summary per filter (DW-P12).
    summ = []
    for filt in sorted({r["filter"] for r in rows_grad}):
        g = [r for r in rows_grad if r["filter"] == filt
             and r["moonangl"] is not None]
        if len(g) < 8:
            continue
        ga = np.array([r["grad_frac_chip"] for r in g])
        rho_a, p_a = spearmanr(ga, [r["moonangl"] for r in g])
        rho_p, p_p = spearmanr(ga, [r["moonphas"] for r in g])
        summ.append({"filter": filt, "n": len(g),
                     "grad_median": float(np.median(ga)),
                     "grad_p90": float(np.percentile(ga, 90)),
                     "rho_moonangle": float(rho_a), "p_moonangle": float(p_a),
                     "rho_moonphase": float(rho_p), "p_moonphase": float(p_p)})
    write_table(con, "dw_gradient_moon", summ)
    meta(con, flats_utc=utcnow())
    for r in rows_val:
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    return 0


# ===========================================================================
# measure  (DW-P03 plate solutions; per-star photometry for P2/P37/P54)
# ===========================================================================
#: The camera's CD matrix (PinPoint solution of a 2023-06-04 Dw1403+49 L
#: frame, "Flip/Mirror": 0.54"/px, roll 0.96 deg, both axes negative) —
#: the translation vote's fixed rotation and scale; the vote also tries the
#: 180-degree-rotated matrix (a meridian flip).  The fitted WCS then
#: measures both.
REF_CD = np.array([[-1.49989e-4, -2.5123e-6], [2.5121e-6, -1.49999e-4]])
#: Saturation veto in RAW native ADU per readout mode (detector_params).
SAT_VETO = {"High Gain": 3200.0, "High Gain StackPro": 51500.0,
            "1MHz High Sensitivity 16-bit": 59500.0}
#: Gain (e-/ADU) for photometric errors: the measured AC4040 High Gain
#: values (detector_params, 1.068-1.071); StackPro sums 16 sub-frames at
#: the same per-read gain.
GAIN = 1.07
PEAK_MEDIAN_FACTOR = 1.10
_CATS: dict = {}


def detector_limits() -> dict:
    """Linearity cap and bias pedestal (raw ADU) per readout mode, from the
    manifest's measured detector_params (AC4040; the High Gain cap is the
    EGAIN-1.054 epoch's, the epoch of every High Gain frame used here)."""
    m = manifest()
    q = lambda g, k: m.execute("SELECT value FROM detector_params WHERE "
                               "era_group=? AND quantity=?", (g, k)).fetchone()[0]
    return {"High Gain": {"cap": q("AC4040 High Gain e1.054", "linearity_cap_adu"),
                          "bias": q("High Gain", "bias_offset_adu")},
            "High Gain StackPro": {"cap": q("High Gain StackPro", "linearity_cap_adu"),
                                   "bias": q("High Gain StackPro", "bias_offset_adu")}}


def _cat(target: str) -> dict:
    if target not in _CATS:
        _CATS[target] = refcat(target)
    return _CATS[target]


def _measure_one(task: dict) -> dict:
    from macro_dw import dwio
    from macro_sn import snio
    t0 = time.time()
    out = {"obs_rowid": task["obs_rowid"], "status": "failed", "stars": None}
    try:
        raw, h = dwio.load(task["path"])
        if task["camera"] == "AC4040":
            img, bad = dwio.calibrate(raw, task["filter"], task["readoutm"],
                                      task["exptime"])
            bp = snio.badpix_mask({"camera": "AC4040",
                                   "readoutm": task["readoutm"],
                                   "flipstat": h.get("FLIPSTAT")}, raw.shape)
            if bp is not None:
                bad = bad | bp
                img[bp] = 0.0
            recipe = "dark+flat_" + task["filter"] + ("+badpix" if bp is not None else "")
        else:
            img, bad = raw - np.median(raw), np.zeros(raw.shape, bool)
            recipe = "none (no in-window masters for this camera)"
        cat = _cat(task["target"])
        sol = dwio.solve(img, bad, h, cat, REF_CD)
        obj = sol["obj"]
        out.update(n_det=int(sol["n_raw"]), recipe=recipe,
                   sky_adu=float(np.median(sol["bkg"].back()[::64, ::64])),
                   sky_rms=float(sol["bkg"].globalrms))
        if sol["wcs"] is None:
            out["status"] = "unsolved"
            return out
        w, sc = sol["wcs"], sol["scale"]
        # Bad pixels are REPAIRED (5x5 median) rather than masked for the
        # star photometry: a warm Hα frame's hot-pixel mask touches most
        # 6" apertures, and sep drops masked pixels from the sum without a
        # correction (flag 32), which silently loses flux.
        from scipy.ndimage import median_filter
        sub_fix = np.where(bad, median_filter(sol["sub"], size=5), sol["sub"])
        ph = dwio.photometry(sub_fix, np.zeros_like(bad), w, cat, GAIN, sc)
        veto = SAT_VETO.get(task["readoutm"], 51500.0)
        # Native peak on the 3x3-median-filtered RAW frame: isolated hot
        # pixels (a warm frame carries one per ~50 px) would otherwise put a
        # ~2000 ADU "peak" in nearly every 4 px search disc and veto every
        # star.  The filter lowers a 4-6 px-FWHM star's peak by <= 10 %,
        # which ``PEAK_MEDIAN_FACTOR`` restores before the cap is applied.
        from scipy.ndimage import median_filter as _mf
        peaks = snio.native_peaks(_mf(raw, size=3), ph["x"], ph["y"], 4.0) \
            * PEAK_MEDIAN_FACTOR
        # FWHM from Gaussian fits at the 40 brightest unsaturated catalogue
        # stars (sep's isophotal second moments truncate the wings and
        # under-read the FWHM by ~2x on these frames).
        f6 = ph["f" + dwio.rkey(dwio.MEAS_APER_ARCSEC)]
        ok = (peaks < veto) & np.isfinite(f6) & (f6 > 0)
        order = np.nonzero(ok)[0][np.argsort(-f6[ok])][:40]
        fw = snio.grid_fwhm(sol["sub"], ph["x"][order], ph["y"][order], half=15)
        xo, yo = w.all_world2pix(task["obj_ra"], task["obj_dec"], 0)
        cra, cdec = w.all_pix2world(2047.5, 2047.5, 0)
        out.update(status="solved", method=sol["method"], n_match=sol["n_match"],
                   rms_arcsec=sol["rms_arcsec"], scale_arcsec=sc, fwhm_px=fw,
                   x_obj=float(xo), y_obj=float(yo), ra_c=float(cra),
                   dec_c=float(cdec), wcs=dwio.wcs_to_text(w))
        st = {"star": ph["idx"], "x": ph["x"], "y": ph["y"], "peak": peaks}
        for r in dwio.APER_RADII_ARCSEC:
            for c in ("f", "e", "flag"):
                st[c + dwio.rkey(r)] = ph[c + dwio.rkey(r)]
        out["stars"] = st
    except Exception as e:                                   # pragma: no cover
        out["status"] = "error"
        out["error"] = repr(e)[:300]
    finally:
        out["t_s"] = time.time() - t0
    return out


SOLVE_COLS = ("obs_rowid", "status", "method", "n_det", "n_match", "rms_arcsec",
              "scale_arcsec", "fwhm_px", "sky_adu", "sky_rms", "x_obj", "y_obj",
              "ra_c", "dec_c", "recipe", "wcs", "error", "t_s")


def cmd_measure(args) -> int:
    from multiprocessing import Pool
    from macro_dw import dwio
    con = connect()
    con.execute(f"CREATE TABLE IF NOT EXISTS dw_solve ({', '.join(SOLVE_COLS)},"
                " PRIMARY KEY(obs_rowid))")
    rcols = [k + dwio.rkey(r) for r in dwio.APER_RADII_ARCSEC
             for k in ("f", "e", "flag")]
    con.execute("CREATE TABLE IF NOT EXISTS dw_phot (obs_rowid INTEGER, "
                "star INTEGER, x REAL, y REAL, peak REAL, "
                + ", ".join(f"{c} REAL" for c in rcols) + ")")
    con.execute("CREATE INDEX IF NOT EXISTS ix_dw_phot ON dw_phot(obs_rowid)")
    if args.refresh:
        con.execute("DELETE FROM dw_solve")
        con.execute("DELETE FROM dw_phot")
    done = {r[0] for r in con.execute("SELECT obs_rowid FROM dw_solve")}
    tasks = [t for t in worklist(con, cameras=("AC4040", "iKon"))
             if t["obs_rowid"] not in done]
    if args.limit:
        tasks = tasks[:args.limit]
    print(f"measure: {len(tasks)} to do", flush=True)
    with Pool(args.workers) as pool:
        for k, r in enumerate(pool.imap_unordered(_measure_one, tasks, 1)):
            con.execute(f"INSERT OR REPLACE INTO dw_solve VALUES "
                        f"({','.join('?' * len(SOLVE_COLS))})",
                        tuple(_py(r.get(c)) for c in SOLVE_COLS))
            st = r.get("stars")
            if st is not None:
                n = len(st["star"])
                con.executemany(
                    f"INSERT INTO dw_phot VALUES ({','.join('?' * (5 + len(rcols)))})",
                    [tuple(_py(v) for v in
                           [r["obs_rowid"], st["star"][i], st["x"][i], st["y"][i],
                            st["peak"][i]] + [st[c][i] for c in rcols])
                     for i in range(n)])
            if k % 25 == 0:
                con.commit()
                print(f"  {k}/{len(tasks)} {r['status']}", flush=True)
    con.commit()
    for row in con.execute("SELECT status, method, count(*) FROM dw_solve "
                           "GROUP BY 1, 2"):
        print(tuple(row))
    return 0


# ===========================================================================
# zp  (DW-P2-ensemble-zp, DW-P21 ZMAG QC, DW-P37 growth curves)
# ===========================================================================
#: Calibration stars: below the measured linearity cap in native pixels,
#: REFCAT2 12 < r < 19 (High Gain's 1800 ADU linearity cap leaves only
#: r >~ 17.5 stars linear on a sharp 256 s L frame, so the catalogue's
#: faint limit is the window's), 0.2 < g-r < 1.2, no REFCAT2 neighbour
#: within 5" (rp1), not flagged variable, S/N > 20 in the 6" aperture.
CAL_R = (12.0, 19.0)
CAL_GR = (0.2, 1.2)
CAL_RP1 = 5.0
CAL_SNR = 20.0
COLOUR_PIVOT = 0.6
#: Growth-curve stars: the same, S/N > 30 in 6"; correction 6" -> 13.5"
#: per frame from >= 5 stars, else the median correction of its
#: (filter, readout mode) group (flagged ``apcor_basis = 'group'``).
GC_SNR = 30.0


def load_phot(con, rowids) -> dict:
    """dw_phot rows for a set of frames, as numpy columns."""
    q = ("SELECT * FROM dw_phot WHERE obs_rowid IN (%s)"
         % ",".join(str(int(r)) for r in rowids))
    rows = con.execute(q).fetchall()
    keys = rows[0].keys() if rows else []
    return {k: np.array([r[k] for r in rows], float) for k in keys}


def fit_zp(dm: np.ndarray, col: np.ndarray, frame: np.ndarray,
           err: np.ndarray, clip: float = 3.0, iters: int = 5):
    """dm = ZP_frame + c (col - pivot): per-frame ZP, one colour term.
    Weighted LSQ, iterated sigma-clipping.  Returns (zp dict, c, sigma_c,
    rms, keep mask)."""
    fr = np.unique(frame)
    idx = np.searchsorted(fr, frame)
    keep = np.isfinite(dm) & np.isfinite(col)
    w = 1.0 / np.maximum(err, 0.005) ** 2
    c = 0.0
    for _ in range(iters):
        # alternate: ZP given c, then c given ZP
        for _ in range(10):
            zp = np.bincount(idx[keep], (w * (dm - c * (col - COLOUR_PIVOT)))[keep],
                             len(fr)) / np.maximum(np.bincount(idx[keep], w[keep],
                                                               len(fr)), 1e-30)
            x = (col - COLOUR_PIVOT)[keep]
            y = (dm - zp[idx])[keep]
            c = float(np.sum(w[keep] * x * y) / np.sum(w[keep] * x * x))
        res = dm - zp[idx] - c * (col - COLOUR_PIVOT)
        s = core.robust_sigma(res[keep])
        keep = np.isfinite(dm) & np.isfinite(col) & (np.abs(res) < clip * s)
    x = (col - COLOUR_PIVOT)[keep]
    resk = res[keep]
    sc = float(np.sqrt(np.sum(w[keep] * resk ** 2) / max(keep.sum() - len(fr) - 1, 1)
                       / np.sum(w[keep] * x * x)))
    nper = np.bincount(idx[keep], minlength=len(fr))
    return dict(zip(fr.tolist(), zp.tolist())), c, sc, float(core.robust_sigma(resk)), keep, nper


def cmd_zp(args) -> int:
    from macro_dw import dwio
    con = connect()
    frames = [dict(r) for r in con.execute(
        "SELECT d.*, s.fwhm_px, s.sky_adu FROM dw_frames d JOIN dw_solve s "
        "USING(obs_rowid) WHERE s.status='solved' AND d.camera='AC4040'")]
    km, k6, k135 = (dwio.rkey(dwio.MEAS_APER_ARCSEC), dwio.rkey(6.0),
                    dwio.rkey(13.5))
    rows_fr, rows_ct, rows_gc, rows_ps1 = [], [], [], []
    groups = sorted({(f["filter"], f["readoutm"]) for f in frames})
    for filt, mode in groups:
        fl = [f for f in frames if f["filter"] == filt and f["readoutm"] == mode]
        P = load_phot(con, [f["obs_rowid"] for f in fl])
        if not P:
            continue
        info = {f["obs_rowid"]: f for f in fl}
        veto = detector_limits()[mode]["cap"]       # linearity cap, raw ADU
        # Per-target catalogue lookups.
        tgt = np.array([info[int(r)]["target"] for r in P["obs_rowid"]])
        g = np.full(len(tgt), np.nan)
        r_ = np.full(len(tgt), np.nan)
        rp1 = np.full(len(tgt), np.nan)
        dv = np.full(len(tgt), np.nan)
        er = np.full(len(tgt), np.nan)
        for t in np.unique(tgt):
            m = tgt == t
            c = _cat(t)
            si = P["star"][m].astype(int)
            g[m], r_[m], rp1[m], dv[m], er[m] = (c["g"][si], c["r"][si],
                                                 c["rp1"][si], c["dupvar"][si],
                                                 c["er"][si])
        exp = np.array([info[int(r)]["exptime"] for r in P["obs_rowid"]])
        f6, e6 = P["f" + km], P["e" + km]
        snr = f6 / np.maximum(e6, 1e-9)
        good = ((P["peak"] < veto) & (P["flag" + km] == 0)
                & (r_ > CAL_R[0]) & (r_ < CAL_R[1])
                & (g - r_ > CAL_GR[0]) & (g - r_ < CAL_GR[1])
                & ((rp1 >= CAL_RP1) | ~np.isfinite(rp1)) & (dv != 1))
        # -- growth curve per frame (DW-P37) --------------------------------
        apcor = {}
        gcs = []
        for rid in np.unique(P["obs_rowid"]):
            m = (P["obs_rowid"] == rid) & good & (snr > GC_SNR) \
                & (P["flag" + k135] == 0) & (P["f" + k135] > 0)
            if m.sum() < 5:
                continue
            ratio = P["f" + k135][m] / f6[m]
            apcor[int(rid)] = float(-2.5 * np.log10(np.median(ratio)))
            gcs.append([np.median(P["f" + dwio.rkey(rr)][m] / P["f" + k135][m])
                        for rr in dwio.APER_RADII_ARCSEC])
        if gcs:
            gcs = np.array(gcs)
            rows_gc.append({"filter": filt, "readoutm": mode, "n_frames": len(gcs),
                            **{f"frac_{dwio.rkey(rr)}": float(np.median(gcs[:, j]))
                               for j, rr in enumerate(dwio.APER_RADII_ARCSEC)},
                            "apcor_median": float(np.median(list(apcor.values()))),
                            "apcor_rms": core.robust_sigma(np.array(list(apcor.values())))})
        ac_group = float(np.median(list(apcor.values()))) if apcor else 0.0
        ac = np.array([apcor.get(int(r), ac_group) for r in P["obs_rowid"]])
        minst = -2.5 * np.log10(np.where(f6 > 0, f6, np.nan) / exp) + ac
        merr = 1.0857 / snr
        sel = good & (snr > CAL_SNR) & np.isfinite(minst)
        dm = (minst - r_)[sel]
        zp, cterm, sc, rms, keep, nper = fit_zp(
            dm, (g - r_)[sel], P["obs_rowid"][sel],
            np.hypot(merr[sel], np.nan_to_num(er[sel])))
        rows_ct.append({"filter": filt, "readoutm": mode, "n_frames": len(zp),
                        "n_meas": int(keep.sum()), "cterm_gr": cterm,
                        "cterm_err": sc, "pivot_gr": COLOUR_PIVOT, "rms": rms})
        fr_ids = np.unique(P["obs_rowid"][sel])
        for k, rid in enumerate(fr_ids):
            f = info[int(rid)]
            if nper[k] < 3:            # no surviving calibration stars
                continue
            z = zp[rid] * -1.0          # m = r  =>  ZP = r + 2.5log(rate)
            rows_fr.append({"obs_rowid": int(rid), "filter": filt,
                            "readoutm": mode, "target": f["target"],
                            "night": f["night"], "zp": z, "n_star": int(nper[k]),
                            "apcor_6_to_13p5": apcor.get(int(rid), ac_group),
                            "apcor_basis": "frame" if int(rid) in apcor else "group",
                            "zmag_hdr": f["zmag"],
                            "zmag_minus_zp": (f["zmag"] - z) if f["zmag"] else None})
        # -- PS1 check (REFCAT2 merge vs PS1 DR1 direct) --------------------
        for t in PS1_FIELDS:
            path = EXT / f"ps1_{t}.tsv"
            m = sel & (tgt == t)
            if not path.exists() or m.sum() < 10:
                continue
            ps = load_tsv(path)
            c = _cat(t)
            si = P["star"][m].astype(int)
            from scipy.spatial import cKDTree
            tree = cKDTree(np.c_[ps["RAJ2000"] * np.cos(np.radians(ps["DEJ2000"])),
                                 ps["DEJ2000"]])
            d, j = tree.query(np.c_[c["ra"][si] * np.cos(np.radians(c["dec"][si])),
                                    c["dec"][si]], distance_upper_bound=1.0 / 3600)
            ok = np.isfinite(d)
            if ok.sum() < 10:
                continue
            jj = j[ok]
            d_r = ps["rmag"][jj] - r_[m][ok]
            rows_ps1.append({"filter": filt, "readoutm": mode, "field": t,
                             "n": int(ok.sum()),
                             "median_ps1_minus_refcat2_r": float(np.nanmedian(d_r)),
                             "rms": core.robust_sigma(d_r)})
        print(f"zp {filt:2s} {mode:20s} frames={len(zp):3d} c={cterm:+.3f}±{sc:.3f} "
              f"rms={rms:.3f}", flush=True)
    write_table(con, "dw_zp_frames", rows_fr)
    write_table(con, "dw_zp_colour", rows_ct)
    write_table(con, "dw_growth", rows_gc)
    write_table(con, "dw_ps1_check", rows_ps1)
    # ZMAG QC flag (DW-P21): ZMAG never calibrates; it is compared.
    zz = [r for r in rows_fr if r["zmag_minus_zp"] is not None]
    by = {}
    for r in zz:
        by.setdefault((r["filter"], r["readoutm"]), []).append(r["zmag_minus_zp"])
    flags = []
    for r in rows_fr:
        if r["zmag_minus_zp"] is None:
            fl = "absent"
        else:
            med = float(np.median(by[(r["filter"], r["readoutm"])]))
            fl = "consistent" if abs(r["zmag_minus_zp"] - med) <= 0.1 else "discrepant"
        flags.append({"obs_rowid": r["obs_rowid"], "zmag_flag": fl})
    write_table(con, "dw_zmag_qc", flags)
    meta(con, zp_utc=utcnow())
    return 0


# ===========================================================================
# qc  (DW-P32 gates; DW-P01 disposition)
# ===========================================================================
def unsolved_diagnosis(r: dict) -> str:
    """Machine diagnosis of a frame that did not solve, from its own
    numbers (thresholds stated in the strings)."""
    veto = SAT_VETO.get(r["readoutm"], 51500.0)
    sky, ndet = r.get("sky_adu") or 0.0, r.get("n_det") or 0
    if sky > 0.85 * (veto - 93.0):
        return "sky within 15% of the High Gain ceiling (twilight/moon)"
    if sky < 5.0:
        return "no sky signal (< 5 ADU: dome closed or overcast)"
    if ndet > 8000:
        return "detections > 8000: noise-dominated frame"
    return "too few REFCAT2 matches (clouds or trailing)"


def cmd_qc(args) -> int:
    con = connect()
    fr = [dict(r) for r in con.execute(
        "SELECT d.*, s.status AS solve_status, s.method, s.n_det, s.n_match, "
        "s.rms_arcsec, s.fwhm_px, s.sky_adu, s.ra_c, s.dec_c, z.zp, z.n_star, "
        "q.zmag_flag "
        "FROM dw_frames d LEFT JOIN dw_solve s USING(obs_rowid) "
        "LEFT JOIN dw_zp_frames z USING(obs_rowid) "
        "LEFT JOIN dw_zmag_qc q USING(obs_rowid)")]
    DL = detector_limits()
    # Night medians and clear-sky ZPs per (filter, readout mode).
    zps: dict = {}
    for r in fr:
        if r["zp"] is not None:
            zps.setdefault((r["filter"], r["readoutm"]), []).append(r)
    clear = {k: float(np.percentile([x["zp"] for x in v], core.QC_ZP_CLEAR_PCTL))
             for k, v in zps.items()}
    night_med = {}
    for k, v in zps.items():
        for n in {x["night"] for x in v}:
            night_med[k + (n,)] = float(np.median([x["zp"] for x in v
                                                   if x["night"] == n]))
    out = []
    for r in fr:
        k = (r["filter"], r["readoutm"])
        solved = None if r["spectral"] else (r["solve_status"] == "solved")
        in_scope = not (r["target"] == "NGC 5548" and not r["spectral"]) \
            and not (r["target"].startswith("Dw") and r["filter"] not in DW_BANDS)
        reasons = []
        lim = DL.get(r["readoutm"])
        if solved and in_scope:
            reasons = core.qc_reasons(
                r["fwhm_px"], r["zp"], night_med.get(k + (r["night"],)),
                clear.get(k),
                sky_raw_adu=(r["sky_adu"] + lim["bias"]) if lim and r["sky_adu"] is not None else None,
                lin_cap_adu=lim["cap"] if lim else None)
        disp = core.disposition(r["is_canonical"], r["filter"],
                                r["offset_field_deg"], solved, reasons)
        if disp == "science" and not in_scope:
            disp = "out_of_scope"
        if r["filter"] == "W" and disp == "science":
            disp = "detection_only (W: no surface photometry, ruling RF)"
        out.append({"obs_rowid": r["obs_rowid"], "target": r["target"],
                    "filter": r["filter"], "readoutm": r["readoutm"],
                    "night": r["night"], "exptime": r["exptime"],
                    "bjd_tdb_mid": r["bjd_tdb_mid"],
                    "disposition": disp, "qc_reasons": ",".join(reasons) or None,
                    "unsolved_diagnosis": (unsolved_diagnosis(r)
                                           if solved is False else None),
                    "solve_method": r["method"], "n_match": r["n_match"],
                    "rms_arcsec": r["rms_arcsec"], "fwhm_px": r["fwhm_px"],
                    "zp": r["zp"], "zp_clear": clear.get(k),
                    "zp_night_median": night_med.get(k + (r["night"],)),
                    "zmag_flag": r["zmag_flag"],
                    "offset_field_deg": r["offset_field_deg"],
                    "header_offset_manifest_deg": r["offset_manifest_deg"],
                    # Pointing column (b): solved chip centre minus the
                    # header's commanded centre, arcsec.
                    "wcs_minus_header_arcsec": (
                        float(3600 * np.hypot((r["ra_c"] - r["ra_deg"])
                                              * np.cos(np.radians(r["dec_c"])),
                                              r["dec_c"] - r["dec_deg"]))
                        if r["ra_c"] is not None and r["ra_deg"] is not None
                        else None)})
    write_table(con, "dw_disposition", out)
    meta(con, qc_utc=utcnow())
    for row in con.execute("SELECT disposition, count(*) FROM dw_disposition "
                           "GROUP BY 1 ORDER BY 2 DESC"):
        print(tuple(row))
    for row in con.execute("SELECT qc_reasons, count(*) FROM dw_disposition "
                           "WHERE qc_reasons IS NOT NULL GROUP BY 1"):
        print(tuple(row))
    return 0


# ===========================================================================
# stacks  (DW-P33 Dw fields; DW-P51 NGC 5238)
# ===========================================================================
#: Common photometric scale of every stack: m = 25 - 2.5 log10(sum of
#: pixel values), on the REFCAT2-r-tied instrumental system of each band.
STACK_ZP = 25.0
STACK_SCALE = 0.54                     #: arcsec/px of the stack grid
STACK_HALF = {"dw": 600, "NGC5238": 800}
#: Per-pixel rejection about the frame median: |v - med| > k sigma_frame
#: + f |med| (the second term keeps PSF-core differences between frames
#: from being clipped as cosmic rays).
CLIP_K, CLIP_F = 5.0, 0.2


def _frame_on_grid(fr: dict, grid, half: int):
    """One science frame: calibrated, constant+plane sky removed (fitted
    on the source-masked frame — never a mesh), scaled to STACK_ZP,
    resampled onto the grid.  Returns (image, sigma per grid px)."""
    from macro_dw import dwio
    raw, h = dwio.load(fr["path"])
    from scipy.ndimage import median_filter
    from macro_sn import snio
    img, bad = dwio.calibrate(raw, fr["filter"], fr["readoutm"], fr["exptime"])
    bp = snio.badpix_mask({"camera": "AC4040", "readoutm": fr["readoutm"],
                           "flipstat": h.get("FLIPSTAT")}, raw.shape)
    if bp is not None:
        bad = bad | bp
    mfile = dwio.MASKDIR / f"{fr['obs_rowid']}.npy"
    mask = np.unpackbits(np.load(mfile))[:img.size].reshape(img.shape).astype(bool) \
        if mfile.exists() else dwio.source_mask(img, bad)
    blk = dwio.block_median(img, mask | bad)
    ny, nx = blk.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    m = np.isfinite(blk)
    A = core.plane_design(xx[m] - nx / 2, yy[m] - ny / 2)
    coef, *_ = np.linalg.lstsq(A, blk[m], rcond=None)
    py, px = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    b = core.FLAT_BLOCK
    sky = (coef[0] + coef[1] * ((px + 0.5) / b - 0.5 - nx / 2)
           + coef[2] * ((py + 0.5) / b - 0.5 - ny / 2))
    img = img - sky
    s0 = core.robust_sigma(img[~(mask | bad)][::37])
    # Hot/RTS pixels not in the detector mask: sharp (more than half the
    # pixel's value above its 3x3 median) AND > 6 sigma.  A 5-px-FWHM star
    # core exceeds its 3x3 median by ~10 %, never by half.
    med3 = median_filter(img, size=3)
    d3 = img - med3
    bad = bad | ((d3 > 6 * s0) & (d3 > 0.5 * img))
    # Column pattern: the April darks do not match the June frames'
    # column offsets (visible as stripes in the faint-sky Hα stacks), so
    # each column's source-masked median is removed.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        col = np.nanmedian(np.where(mask | bad, np.nan, img), axis=0)
    img = img - np.nan_to_num(col)[None, :]
    sig_pix = core.robust_sigma(img[~(mask | bad)][::37])
    img[bad] = np.nan
    w = dwio.wcs_from_text(fr["wcs"])
    k = 10 ** (0.4 * (STACK_ZP - fr["zp"])) / fr["exptime"]
    out = dwio.resample(img, w, grid, half, STACK_SCALE) * k
    # Bilinear resampling smooths white noise; the sigma used for weights
    # and clipping is the native one, scaled (the stack's true noise is
    # measured afterwards, empirically, never from this number).
    a_img = abs(np.linalg.det(w.pixel_scale_matrix)) * 3600.0 ** 2
    return out, sig_pix * k * STACK_SCALE ** 2 / a_img


def _stack_job(job: dict) -> dict:
    from astropy.io import fits
    from macro_dw import dwio
    half = job["half"]
    grid = dwio.make_grid(job["ra"], job["dec"], half, STACK_SCALE)
    ims, sig = [], []
    for fr in job["frames"]:
        try:
            a, s_ = _frame_on_grid(fr, grid, half)
        except Exception:
            continue
        ims.append(a)
        sig.append(s_)
    if not ims:
        return {"name": job["name"], "n": 0}
    cube = np.array(ims, np.float32)
    sig = np.array(sig)
    w = 1.0 / sig ** 2
    med = np.nanmedian(cube, axis=0)
    bad = np.abs(cube - med) > (CLIP_K * sig[:, None, None] + CLIP_F * np.abs(med))
    good = np.isfinite(cube) & ~bad
    W = (good * w[:, None, None]).sum(0)
    sci = np.where(W > 0, np.nansum(np.where(good, cube, 0) * w[:, None, None], 0)
                   / np.maximum(W, 1e-30), np.nan)
    nfr = good.sum(0).astype(np.int16)
    hdr = grid.to_header()
    hdr["BUNIT"] = "flux (ZP=25)"
    hdr["MAGZP"] = STACK_ZP
    hdr["NFRAMES"] = len(ims)
    hdr["FILTER"] = job["filter"]
    hdr["SKYMODEL"] = "constant+plane per frame (source-masked blocks)"
    hdr["FWHMPX"] = float(np.median([f["fwhm_px"] for f in job["frames"]]))
    path = dwio.STACKDIR / f"{job['name']}.fits"
    fits.HDUList([fits.PrimaryHDU(sci.astype(np.float32), hdr),
                  fits.ImageHDU(W.astype(np.float32), name="WEIGHT"),
                  fits.ImageHDU(nfr, name="NFRAMES")]).writeto(path, overwrite=True)
    return {"name": job["name"], "target": job["target"], "filter": job["filter"],
            "night": job.get("night"), "n": len(ims),
            "exp_s": float(sum(f["exptime"] for f in job["frames"])),
            "fwhm_frames_med": hdr["FWHMPX"],
            "clip_frac": float(bad[np.isfinite(cube)].mean()),
            "path": str(path.relative_to(REPO))}


def stack_jobs(con) -> list[dict]:
    fr = [dict(r) for r in con.execute(
        "SELECT d.*, s.wcs, s.fwhm_px, z.zp FROM dw_frames d "
        "JOIN dw_disposition q USING(obs_rowid) JOIN dw_solve s USING(obs_rowid) "
        "JOIN dw_zp_frames z USING(obs_rowid) WHERE q.disposition='science'")]
    jobs = []
    keys = sorted({(f["target"], f["filter"]) for f in fr})
    for t, b in keys:
        sel = [f for f in fr if f["target"] == t and f["filter"] == b]
        ra, de = object_radec(t)
        half = STACK_HALF["NGC5238" if t == "NGC5238" else "dw"]
        name = f"{t.replace(' ', '')}_{b}"
        jobs.append({"name": name, "target": t, "filter": b, "frames": sel,
                     "ra": ra, "dec": de, "half": half})
        if t == "NGC5238" and b in ("H", "R"):
            for n in sorted({f["night"] for f in sel}):
                jobs.append({"name": f"{name}_{n}", "target": t, "filter": b,
                             "night": n, "ra": ra, "dec": de, "half": half,
                             "frames": [f for f in sel if f["night"] == n]})
    # Largest first, so the long NGC 5238 jobs do not finish last alone.
    return sorted(jobs, key=lambda j: -len(j["frames"]))


def cmd_stacks(args) -> int:
    from multiprocessing import Pool
    from macro_dw import dwio
    dwio.STACKDIR.mkdir(parents=True, exist_ok=True)
    # Every stack is rebuilt from the current disposition; a stack left over
    # from an earlier run (a band whose frames have since been rejected)
    # must not survive to be measured.
    if not args.limit:
        for old in dwio.STACKDIR.glob("*.fits"):
            old.unlink()
    con = connect()
    jobs = stack_jobs(con)
    if args.limit:
        jobs = [j for j in jobs if j["target"] != "NGC5238"][:args.limit]
    print(f"stacks: {len(jobs)} jobs", flush=True)
    with Pool(args.workers, maxtasksperchild=4) as pool:
        res = list(pool.imap_unordered(_stack_job, jobs))
    write_table(con, "dw_stacks", sorted(res, key=lambda r: r["name"]))
    meta(con, stacks_utc=utcnow())
    print(f"stacks: {sum(1 for r in res if r['n'])} written")
    return 0


# ===========================================================================
# halpha  (DW-P36; the NGC 5238 line-flux scale)
# ===========================================================================
HA_ANN_ARCSEC = (20.0, 40.0)
N_RANDOM_AP = 400
#: NGC 5238 curve of growth: radii (arcsec) and the plateau rule.
CAL_RADII = np.arange(10.0, 155.0, 5.0)
CAL_PLATEAU = 0.02
#: Colour-term difference sets how continuum colour moves k; the colour
#: uncertainty used when L-R is unmeasured (pivot assumed).
GR_FALLBACK_SIGMA = 0.3


def _colour_terms() -> dict:
    con = connect(read_only=True)
    out = {}
    for r in con.execute("SELECT filter, readoutm, cterm_gr, n_meas FROM dw_zp_colour"):
        out.setdefault(r["filter"], []).append((r["n_meas"], r["cterm_gr"]))
    return {f: float(np.average([c for _, c in v], weights=[n for n, _ in v]))
            for f, v in out.items()}


def _ha_prepare(target: str, seed: int) -> dict:
    """PSF-matched H, R (and L) stacks of one field, the continuum-
    subtracted image, its mask and the object's pixel position."""
    from macro_dw import dwio
    key = target.replace(" ", "")
    H, _, w, hh = dwio.read_stack(f"{key}_H")
    R, _, _, _ = dwio.read_stack(f"{key}_R")
    L = None
    if (dwio.STACKDIR / f"{key}_L.fits").exists():
        L, _, _, _ = dwio.read_stack(f"{key}_L")
    cat = _cat(target)
    fw = {b: dwio.stack_fwhm(im, w, cat) for b, im in (("H", H), ("R", R), ("L", L))
          if im is not None}
    target_fw = max(v for v in fw.values() if v)
    def match(im, f):
        return dwio._nan_gauss(im, np.sqrt(max(target_fw ** 2 - f ** 2, 0)) / 2.3548) \
            if f and target_fw - f > 0.05 else im
    H, R = match(H, fw["H"]), match(R, fw["R"])
    if L is not None:
        L = match(L, fw.get("L"))
    ra, de = object_radec(target)
    x0, y0 = (float(v) for v in w.all_world2pix(ra, de, 0))
    mask = dwio.stack_mask(R) | dwio.stack_mask(H) | ~np.isfinite(H) | ~np.isfinite(R)
    return {"H": H, "R": R, "L": L, "w": w, "hdr": hh, "x0": x0, "y0": y0,
            "mask": mask, "fwhm": fw, "fwhm_matched": target_fw,
            "rng": np.random.default_rng(seed)}


def _unmask_disc(mask, x0, y0, r):
    m = mask.copy()
    yy, xx = np.ogrid[:m.shape[0], :m.shape[1]]
    m[(xx - x0) ** 2 + (yy - y0) ** 2 <= r * r] = False
    return m


def _continuum_scale(P: dict, r_px: float, ann, ct: dict):
    """k and its uncertainty from the object's own L-R colour."""
    from macro_dw import dwio
    dc = ct["H"] - ct["R"]
    gr, sgr, basis = COLOUR_PIVOT, GR_FALLBACK_SIGMA, "pivot"
    if P["L"] is not None:
        mo = _unmask_disc(P["mask"], P["x0"], P["y0"], r_px)
        fL = dwio.aper_net(P["L"], P["x0"], P["y0"], r_px, ann, mo)
        fR = dwio.aper_net(P["R"], P["x0"], P["y0"], r_px, ann, mo)
        sL = core.robust_sigma(dwio.random_apertures(P["L"], P["mask"], r_px, ann,
                                                     150, P["rng"]))
        sR = core.robust_sigma(dwio.random_apertures(P["R"], P["mask"], r_px, ann,
                                                     150, P["rng"]))
        if fL > 10 * sL and fR > 10 * sR:
            col = -2.5 * np.log10(fL / fR)
            scol = 1.0857 * np.hypot(sL / fL, sR / fR)
            gr = COLOUR_PIVOT + col / (ct["L"] - ct["R"])
            sgr, basis = scol / abs(ct["L"] - ct["R"]), "L-R"
    k = 10 ** (-0.4 * dc * (gr - COLOUR_PIVOT))
    dk = k * 0.4 * np.log(10) * abs(dc) * sgr
    return k, dk, gr, sgr, basis


def _ha_measure(P: dict, r_arcsec: float, ct: dict) -> dict:
    from macro_dw import dwio
    r = r_arcsec / STACK_SCALE
    ann = tuple(a / STACK_SCALE for a in HA_ANN_ARCSEC)
    k, dk, gr, sgr, basis = _continuum_scale(P, r, ann, ct)
    D = P["H"] - k * P["R"]
    mo = _unmask_disc(P["mask"], P["x0"], P["y0"], r)
    net = dwio.aper_net(D, P["x0"], P["y0"], r, ann, mo)
    cont = dwio.aper_net(P["R"], P["x0"], P["y0"], r, ann, mo)
    rnd = dwio.random_apertures(D, P["mask"], r, ann, N_RANDOM_AP, P["rng"],
                                avoid=(P["x0"], P["y0"], 60 / STACK_SCALE))
    s_rand = core.robust_sigma(rnd)
    s_tot = float(np.hypot(s_rand, cont * dk))
    return {"k": k, "dk": dk, "gr_obj": gr, "gr_sigma": sgr, "k_basis": basis,
            "net": net, "cont_R": cont, "sigma_rand": s_rand,
            "rand_mean": float(np.mean(rnd)), "n_rand": len(rnd),
            "sigma_total": s_tot, "snr": net / s_tot if s_tot > 0 else np.nan,
            "D": D}


def cmd_halpha(args) -> int:
    from astropy.io import fits
    from macro_dw import dwio
    con = connect()
    ct = _colour_terms()
    # -- 1. NGC 5238: the end-to-end line-flux scale --------------------------
    P = _ha_prepare("NGC5238", 5238)
    ann_cog = []
    nets = []
    # Continuum colour of NGC 5238 in a 30" aperture (bright: S/N >> 10).
    k, dk, gr, sgr, basis = _continuum_scale(
        P, 30 / STACK_SCALE, tuple(a / STACK_SCALE for a in (60.0, 90.0)), ct)
    D = P["H"] - k * P["R"]
    for rr in CAL_RADII:
        r = rr / STACK_SCALE
        ann = ((rr + 20) / STACK_SCALE, (rr + 40) / STACK_SCALE)
        nets.append(dwio.aper_net(D, P["x0"], P["y0"], r, ann,
                                  _unmask_disc(P["mask"], P["x0"], P["y0"], r)))
        ann_cog.append(ann)
    nets = np.array(nets)
    rel = np.diff(nets) / np.maximum(np.abs(nets[1:]), 1e-30)
    j = int(np.argmax(np.abs(rel) < CAL_PLATEAU)) if (np.abs(rel) < CAL_PLATEAU).any() \
        else len(CAL_RADII) - 1
    r_cal = float(CAL_RADII[j])
    net_cal = float(nets[j])
    S25 = NGC5238_FHA / net_cal
    w_eff = core.effective_width(S25, STACK_ZP)
    # per-night scatter
    rows_n = []
    nights = [r["night"] for r in con.execute(
        "SELECT night FROM dw_stacks WHERE target='NGC5238' AND filter='H' "
        "AND night IS NOT NULL AND n > 0 INTERSECT SELECT night FROM dw_stacks "
        "WHERE target='NGC5238' AND filter='R' AND night IS NOT NULL AND n > 0")]
    for n in nights:
        Hn, _, _, _ = dwio.read_stack(f"NGC5238_H_{n}")
        Rn, _, _, _ = dwio.read_stack(f"NGC5238_R_{n}")
        Dn = Hn - k * Rn
        r = r_cal / STACK_SCALE
        ann = ((r_cal + 20) / STACK_SCALE, (r_cal + 40) / STACK_SCALE)
        v = dwio.aper_net(Dn, P["x0"], P["y0"], r, ann,
                          _unmask_disc(P["mask"], P["x0"], P["y0"], r))
        rows_n.append({"night": n, "net": v, "S25": NGC5238_FHA / v if v > 0 else None})
    sv = np.array([r["S25"] for r in rows_n if r["S25"]])
    cal = {"r_cal_arcsec": r_cal, "net_cal": net_cal, "k": k, "gr_obj": gr,
           "k_basis": basis, "S25": S25, "S25_rel_err_pub": NGC5238_FHA_ERR / NGC5238_FHA,
           "n_nights": len(sv), "S25_night_mean": float(np.mean(sv)),
           "S25_night_rel_scatter": float(np.std(sv, ddof=1) / np.mean(sv)),
           "W_eff_A": w_eff, "fwhm_matched_px": P["fwhm_matched"],
           "dlam_ngc5238_A": core.ha_offset_A(NGC5238_CZ)}
    write_table(con, "dw_ha_cal", [cal])
    write_table(con, "dw_ha_cal_cog", [{"r_arcsec": float(a), "net": float(b)}
                                       for a, b in zip(CAL_RADII, nets)])
    write_table(con, "dw_ha_cal_nights", rows_n)
    fits.writeto(dwio.STACKDIR / "NGC5238_HaSub.fits", (P["H"] - k * P["R"]).astype(np.float32),
                 P["w"].to_header(), overwrite=True)
    print({k_: (round(v, 4) if isinstance(v, float) else v) for k_, v in cal.items()})
    # -- 2. every Dw field with Hα ------------------------------------------
    lit = literature()
    rows = []
    fields = [r["target"] for r in con.execute(
        "SELECT DISTINCT target FROM dw_stacks WHERE filter='H' AND target LIKE 'Dw%' "
        "AND night IS NULL AND n > 0 ORDER BY target")]
    have_r = {r[0] for r in con.execute(
        "SELECT target FROM dw_stacks WHERE filter='R' AND night IS NULL AND n > 0")}
    for t in [t for t in lit if lit[t]["n_H"] not in ("", "0") and t not in fields + []]:
        fields.append(t)                    # Hα frames exist but none survived QC
    fields = sorted(set(fields))
    for i, t in enumerate(fields):
        L = lit[t]
        if t not in have_r or not (dwio.STACKDIR / f"{t.replace(' ', '')}_H.fits").exists():
            v = float(L["v_hel_kms"]) if L["v_hel_kms"] else None
            dl = core.ha_offset_A(v)
            rows.append({"field": t, "candidate": L["candidate"], "v_hel_kms": v,
                         "dlam_A": dl, "band": ("unknown" if dl is None else "in"
                                                if abs(dl) <= core.W_ASSUMED_A / 2 else "out"),
                         "spec_emission": L["spec_emission"] or None,
                         "verdict": ("no_continuum" if t not in have_r else "no_halpha_frames"),
                         "informative_limit": 0})
            print(f"{t}: no measurement ({rows[-1]['verdict']})")
            continue
        P = _ha_prepare(t, 1000 + i)
        m = _ha_measure(P, core.HA_APER_ARCSEC, ct)
        fits.writeto(dwio.STACKDIR / f"{t}_HaSub.fits", m.pop("D").astype(np.float32),
                     P["w"].to_header(), overwrite=True)
        v = float(L["v_hel_kms"]) if L["v_hel_kms"] else None
        dl = core.ha_offset_A(v)
        band = ("unknown" if dl is None else
                "in" if abs(dl) <= core.W_ASSUMED_A / 2 else "out")
        verdict = core.ha_verdict(m["snr"], band)
        flux = m["net"] * S25
        lim = core.upper_limit(m["net"], m["sigma_total"]) * S25
        fpred = (float(L["F_pred_FUV"]) if L["F_pred_FUV"] else
                 float(L["F_pred_P0"]) if L["F_pred_P0"] else None)
        rows.append({"field": t, "candidate": L["candidate"], "v_hel_kms": v,
                     "dlam_A": dl, "band": band, "spec_emission": L["spec_emission"] or None,
                     "n_H": con.execute("SELECT n FROM dw_stacks WHERE name=?",
                                        (f"{t}_H",)).fetchone()[0],
                     **{k_: m[k_] for k_ in ("k", "dk", "gr_obj", "gr_sigma", "k_basis",
                                             "net", "cont_R", "sigma_rand", "rand_mean",
                                             "n_rand", "sigma_total", "snr")},
                     "verdict": verdict,
                     "flux_5238scale": flux if verdict == "detected" else None,
                     "flux_err": m["sigma_total"] * S25,
                     "flux_lim_3sig": lim if verdict != "detected" else None,
                     "F_pred": fpred, "F_pred_basis": L["F_pred_basis"] or None,
                     "pred_over_lim": (fpred / lim) if (fpred and verdict != "detected") else None,
                     "F_lim_novelty": float(L["F_lim_3sig"]) if L["F_lim_3sig"] else None,
                     "informative_limit": int(bool(fpred and verdict != "detected"
                                                   and fpred >= 3 * lim)),
                     "fwhm_matched_px": P["fwhm_matched"]})
        print(f"{t}: band={band} snr={m['snr']:+.2f} k={m['k']:.3f}({m['k_basis']}) "
              f"lim={lim:.2e} pred={fpred}", flush=True)
    cols = []
    for r in rows:
        cols += [k_ for k_ in r if k_ not in cols]
    write_table(con, "dw_halpha", [{c_: r.get(c_) for c_ in cols} for r in rows])
    inb = [r for r in rows if r["band"] == "in"]
    n_ok = sum(1 for r in inb if r["verdict"] == "detected" or r["informative_limit"])
    venue = "AJ" if n_ok >= 2 else "RNAAS"
    meta(con, halpha_utc=utcnow(), venue_trigger_n=n_ok, venue=venue)
    print(f"venue trigger: {n_ok} in-band objects detected or informative -> {venue}")
    return 0


# ===========================================================================
# sersic  (DW-P35) and depth  (DW-P34)
# ===========================================================================
SERSIC_BIN = 3                         #: 3x3 binning -> 1.62"/px
SERSIC_HALF_ARCSEC = 90.0
DET_APER_ARCSEC = 10.0


def _sersic_fit(cut: np.ndarray, mask: np.ndarray, pix: float):
    """Sérsic + constant, least squares on unmasked pixels.  Returns
    (params dict, success)."""
    from scipy.optimize import least_squares
    ny, nx = cut.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    m = ~mask & np.isfinite(cut)
    c0 = (nx - 1) / 2

    def model(p):
        amp, re, n, x0, y0, e, th, b = p
        ct, st = np.cos(th), np.sin(th)
        dx, dy = xx - x0, yy - y0
        u = dx * ct + dy * st
        v = (-dx * st + dy * ct) / max(1 - e, 1e-3)
        r = np.hypot(u, v)
        bn = core.sersic_bn(n)
        return amp * np.exp(-bn * ((r / re) ** (1 / n) - 1)) + b

    def res(p):
        return (model(p) - cut)[m]
    peak = np.nanmax(cut[int(c0) - 3:int(c0) + 4, int(c0) - 3:int(c0) + 4])
    p0 = [max(peak / 3, 1e-6), 8.0 / pix, 1.0, c0, c0, 0.3, 0.0, 0.0]
    lo = [0, 1.0, 0.3, c0 - 6, c0 - 6, 0.0, -np.pi, -np.inf]
    hi = [np.inf, 60.0 / pix, 4.0, c0 + 6, c0 + 6, 0.9, np.pi, np.inf]
    try:
        r = least_squares(res, p0, bounds=(lo, hi), loss="soft_l1",
                          f_scale=core.robust_sigma(cut[m]) * 2)
    except Exception:
        return None, False
    J = r.jac
    dof = max(m.sum() - len(p0), 1)
    s2 = np.sum(r.fun ** 2) / dof
    try:
        cov = np.linalg.inv(J.T @ J) * s2
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        err = np.full(len(p0), np.nan)
    keys = ("amp", "re_px", "n", "x0", "y0", "ell", "theta", "bkg")
    out = {k: float(v) for k, v in zip(keys, r.x)}
    out.update({f"{k}_err": float(e) for k, e in zip(keys, err)})
    at_bound = any(np.isclose(r.x[i], b, rtol=1e-3) for i, b in
                   ((1, hi[1]), (2, lo[2]), (2, hi[2]), (5, hi[5])))
    return out, bool(r.success and not at_bound)


def cmd_sersic(args) -> int:
    from macro_dw import dwio
    con = connect()
    ct = _colour_terms()
    rows = []
    fields = [r["target"] for r in con.execute(
        "SELECT target FROM dw_stacks WHERE filter='L' AND target LIKE 'Dw%' "
        "AND night IS NULL AND n > 0 ORDER BY target")]
    lit = literature()
    for i, t in enumerate(fields):
        L, _, w, _ = dwio.read_stack(f"{t}_L")
        ra, de = object_radec(t)
        x0, y0 = (float(v) for v in w.all_world2pix(ra, de, 0))
        mask = dwio.stack_mask(L)
        rng = np.random.default_rng(2000 + i)
        r = DET_APER_ARCSEC / STACK_SCALE
        ann = tuple(a / STACK_SCALE for a in HA_ANN_ARCSEC)
        net = dwio.aper_net(L, x0, y0, r, ann, _unmask_disc(mask, x0, y0, r))
        sig = core.robust_sigma(dwio.random_apertures(L, mask, r, ann, 300, rng,
                                                      avoid=(x0, y0, 110)))
        snr = net / sig
        row = {"field": t, "candidate": lit[t]["candidate"], "L_net10": net,
               "L_sigma10": sig, "L_snr10": snr, "detected": int(snr >= 5),
               "m_L10": STACK_ZP - 2.5 * np.log10(net) if net > 0 else None}
        Rn = None
        if (dwio.STACKDIR / f"{t}_R.fits").exists():
            Rim, _, _, _ = dwio.read_stack(f"{t}_R")
            Rn = dwio.aper_net(Rim, x0, y0, r, ann, _unmask_disc(mask, x0, y0, r))
            sR = core.robust_sigma(dwio.random_apertures(Rim, mask, r, ann, 200, rng,
                                                         avoid=(x0, y0, 110)))
            if Rn > 0 and net > 0:
                row["L_minus_R10"] = -2.5 * np.log10(net / Rn)
                row["L_minus_R10_err"] = 1.0857 * np.hypot(sig / net, sR / Rn)
                row["gr_inferred"] = COLOUR_PIVOT + row["L_minus_R10"] / (ct["L"] - ct["R"])
        if snr >= 5:
            b = SERSIC_BIN
            h = int(SERSIC_HALF_ARCSEC / STACK_SCALE)
            cx, cy = int(round(x0)), int(round(y0))
            cut = L[cy - h:cy + h + 1, cx - h:cx + h + 1]
            mk = mask[cy - h:cy + h + 1, cx - h:cx + h + 1].copy()
            # unmask the object's own footprint: the segment touching centre
            from scipy.ndimage import label
            lab, _ = label(mk)
            lc = lab[h, h]
            if lc:
                mk[lab == lc] = False
            n = (cut.shape[0] // b) * b
            cb = cut[:n, :n].reshape(n // b, b, n // b, b).mean(axis=(1, 3))
            mb = mk[:n, :n].reshape(n // b, b, n // b, b).any(axis=(1, 3))
            pix = STACK_SCALE * b
            fit, ok = _sersic_fit(cb, mb, pix)
            if fit:
                re_as = fit["re_px"] * pix
                mu_e = STACK_ZP - 2.5 * np.log10(max(fit["amp"], 1e-30) / pix ** 2)
                row.update(sersic_ok=int(ok), re_arcsec=re_as,
                           re_err=fit["re_px_err"] * pix, n_sersic=fit["n"],
                           n_err=fit["n_err"], ell=fit["ell"], ell_err=fit["ell_err"],
                           mu_e_L=mu_e, mu0_L=core.sersic_mu0(mu_e, fit["n"]),
                           m_tot_L=core.sersic_total_mag(mu_e, re_as, fit["n"],
                                                         1 - fit["ell"]),
                           re_over_fwhm=re_as / (5.0 * STACK_SCALE))
        rows.append(row)
        print(t, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()
                  if k in ("L_snr10", "re_arcsec", "n_sersic", "mu0_L", "sersic_ok")},
              flush=True)
    cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in rows[0], k))
    write_table(con, "dw_sersic", [{c: r.get(c) for c in cols} for r in rows])
    meta(con, sersic_utc=utcnow())
    return 0


def _depth_one(job: dict) -> list[dict]:
    from macro_dw import dwio
    t, band = job["field"], job["band"]
    img, _, w, _ = dwio.read_stack(f"{t}_{band}")
    mask = dwio.stack_mask(img)
    ra, de = object_radec(t)
    x0, y0 = (float(v) for v in w.all_world2pix(ra, de, 0))
    rng = np.random.default_rng(job["seed"])
    ny, nx = img.shape
    box = int(round(core.ROMAN_BOX_ARCSEC / STACK_SCALE))
    means = []
    tries = 0
    while len(means) < 500 and tries < 20000:
        tries += 1
        x = int(rng.uniform(60, nx - 60 - box))
        y = int(rng.uniform(60, ny - 60 - box))
        if np.hypot(x - x0, y - y0) < 110 or mask[y:y + box, x:x + box].any():
            continue
        means.append(float(np.mean(img[y:y + box, x:x + box])))
    sbox = core.robust_sigma(np.array(means))
    out = [{"field": t, "band": band, "kind": "roman", "n_boxes": len(means),
            "sigma_box_mean": sbox,
            "mu_lim_3sig_10as": core.roman_mu_limit(sbox, STACK_ZP, STACK_SCALE)}]
    if band != "L":
        return out
    # Synthetic exponential discs, injected one at a time into the stack.
    ann_f = (1.5, 2.5)
    for re in core.INJ_RE_ARCSEC:
        rpx = max(re, 3.0) / STACK_SCALE
        ann = (max(ann_f[0] * rpx, rpx + 5), max(ann_f[1] * rpx, rpx + 15))
        rnd = dwio.random_apertures(img, mask, rpx, ann, 150, rng,
                                    avoid=(x0, y0, 110))
        srnd = core.robust_sigma(rnd)
        half = int(min(5 * re / STACK_SCALE, 300)) + int(ann[1]) + 2
        for mu0 in core.INJ_MU0:
            rec, bias = [], []
            for _ in range(core.INJ_PER_CELL):
                for _try in range(200):
                    x = rng.uniform(half + 5, nx - half - 5)
                    y = rng.uniform(half + 5, ny - half - 5)
                    if np.hypot(x - x0, y - y0) > 110 + half and \
                            not mask[int(y), int(x)]:
                        break
                cx, cy = int(x), int(y)
                st = core.exp_disc_image(2 * half + 1, x - cx + half, y - cy + half,
                                         mu0, re, STACK_ZP, STACK_SCALE,
                                         q=rng.uniform(0.5, 1.0),
                                         pa_deg=rng.uniform(0, 180))
                sub = img[cy - half:cy + half + 1, cx - half:cx + half + 1] + st
                ms = mask[cy - half:cy + half + 1, cx - half:cx + half + 1]
                xs, ys = x - cx + half, y - cy + half
                ms2 = _unmask_disc(ms, xs, ys, rpx)
                f_in = dwio.aper_net(sub, xs, ys, rpx, ann, ms2)
                f_bg = dwio.aper_net(sub - st, xs, ys, rpx, ann, ms2)
                f_true = dwio.aper_net(st, xs, ys, rpx, ann, None)
                rec.append(f_in / srnd >= core.DETECT_INJ_SIGMA)
                bias.append(((f_in - f_bg) - f_true) / max(f_true, 1e-30))
            out.append({"field": t, "band": band, "kind": "inject", "mu0": float(mu0),
                        "re_arcsec": float(re), "n": len(rec),
                        "frac": float(np.mean(rec)),
                        "signed_flux_bias": float(np.mean(bias)),
                        "sigma_rand_ap": srnd})
    return out


def cmd_depth(args) -> int:
    from multiprocessing import Pool
    con = connect()
    jobs = []
    for i, r in enumerate(con.execute(
            "SELECT target, filter FROM dw_stacks WHERE target LIKE 'Dw%' AND "
            "night IS NULL AND n > 0 ORDER BY target, filter")):
        jobs.append({"field": r["target"], "band": r["filter"], "seed": 3000 + i})
    with Pool(args.workers) as pool:
        res = [x for part in pool.map(_depth_one, jobs) for x in part]
    write_table(con, "dw_depth", [r for r in res if r["kind"] == "roman"])
    write_table(con, "dw_injection", [r for r in res if r["kind"] == "inject"])
    # 50 % / 90 % contours per field.
    cont = []
    for t in sorted({r["field"] for r in res if r["kind"] == "inject"}):
        g = [r for r in res if r["kind"] == "inject" and r["field"] == t]
        F = np.array([[next(x["frac"] for x in g if x["mu0"] == m and x["re_arcsec"] == re)
                       for re in core.INJ_RE_ARCSEC] for m in core.INJ_MU0])
        for lev in (0.5, 0.9):
            c = core.recovery_contour(F, core.INJ_MU0, lev)
            for re, mu in zip(core.INJ_RE_ARCSEC, c):
                cont.append({"field": t, "level": lev, "re_arcsec": float(re),
                             "mu0_lim": None if not np.isfinite(mu) else float(mu)})
    write_table(con, "dw_injection_contours", cont)
    meta(con, depth_utc=utcnow())
    for r in con.execute("SELECT band, count(*), round(min(mu_lim_3sig_10as),2), "
                         "round(avg(mu_lim_3sig_10as),2), round(max(mu_lim_3sig_10as),2) "
                         "FROM dw_depth GROUP BY band"):
        print(tuple(r))
    return 0


# ===========================================================================
# zeroorder  (DW-P4y: NGC 5548 zero-order differential photometry)
# ===========================================================================
#: Zero-order rule (measured on the frames, 2026-10-05): the AGN's zero
#: order is the brightest 8-sigma source with a/b < 5 (it is 14-31 px in
#: sep's a, elongation 3-4: host + undispersed AGN); first-order spectra
#: have a/b > 10 and hot-pixel blobs a < 2 px.  A zero order is therefore a
#: source with 2 < a < 40 px and a/b < 5.  Aperture 25 px, annulus 35-50 px.
ZO_AMIN, ZO_AMAX, ZO_ELONG = 2.0, 40.0, 5.0
ZO_R, ZO_ANN = 25.0, (35.0, 50.0)
ZO_MATCH_PX = 15.0
ZO_RMS_MAX = 0.02            #: the ledger's criterion


def _zo_one(task: dict) -> dict:
    import sep
    from macro_dw import dwio
    from macro_sn import snio
    out = {"obs_rowid": task["obs_rowid"], "night": task["night"], "src": None}
    try:
        raw, h = dwio.load(task["path"])
        img = np.ascontiguousarray(raw.astype(np.float32))
        bkg = sep.Background(img, bw=32, bh=32)
        sub = img - bkg.back()
        sep.set_extract_pixstack(5_000_000)
        obj = sep.extract(sub, 8.0, err=bkg.globalrms, minarea=20)
        el = obj["a"] / np.maximum(obj["b"], 1e-3)
        c = (obj["a"] > ZO_AMIN) & (obj["a"] < ZO_AMAX) & (el < ZO_ELONG)
        o = obj[c]
        f, e, fl = sep.sum_circle(sub, o["x"], o["y"], ZO_R, err=bkg.globalrms,
                                  gain=GAIN, bkgann=ZO_ANN)
        pk = snio.native_peaks(raw, o["x"], o["y"], 5.0)
        out["src"] = {"x": o["x"], "y": o["y"], "f": f, "e": e, "peak": pk}
    except Exception as ex:                                  # pragma: no cover
        out["error"] = repr(ex)[:200]
    return out


def cmd_zeroorder(args) -> int:
    from multiprocessing import Pool
    con = connect()
    tasks = [dict(r) for r in con.execute(
        "SELECT d.* FROM dw_frames d JOIN dw_disposition q USING(obs_rowid) "
        "WHERE d.target='NGC 5548' AND d.filter='6' AND q.disposition='spectrum' "
        "ORDER BY d.jd")]
    with Pool(args.workers) as pool:
        res = pool.map(_zo_one, tasks)
    veto = SAT_VETO["High Gain StackPro"]
    rows_f, rows_n = [], []
    for night in sorted({r["night"] for r in res}):
        fr = [r for r in res if r["night"] == night and r["src"] is not None
              and len(r["src"]["f"]) > 0]
        if len(fr) < 3:
            continue
        # The AGN zero order = the brightest compact source of the night's
        # first frame; every source is tracked by its offset from it.
        def rel(r):
            s_ = r["src"]
            k = int(np.argmax(s_["f"]))
            return s_["x"] - s_["x"][k], s_["y"] - s_["y"][k], k
        rx0, ry0, k0 = rel(fr[0])
        nsrc = len(rx0)
        F = np.full((nsrc, len(fr)), np.nan)
        E = np.full_like(F, np.nan)
        PK = np.full_like(F, np.nan)
        for j, r in enumerate(fr):
            rx, ry, _ = rel(r)
            for i in range(nsrc):
                d = np.hypot(rx - rx0[i], ry - ry0[i])
                if d.size and d.min() < ZO_MATCH_PX:
                    q = int(np.argmin(d))
                    F[i, j], E[i, j], PK[i, j] = (r["src"]["f"][q], r["src"]["e"][q],
                                                  r["src"]["peak"][q])
        ok = (np.isfinite(F).mean(axis=1) >= 0.8) & (np.nanmedian(F / E, axis=1) > 50) \
            & (np.nanmax(PK, axis=1) < veto)
        comps = [i for i in np.nonzero(ok)[0] if i != k0]
        rms_c = []
        for i in comps:
            others = [c for c in comps if c != i]
            if not others:
                continue
            ref = np.nansum(F[others], axis=0)
            dm = -2.5 * np.log10(F[i] / ref)
            dm = dm[np.isfinite(dm)]
            if dm.size >= 3:
                rms_c.append(float(np.std(dm, ddof=1)))
        agn_rms = None
        if comps:
            ref = np.nansum(F[comps], axis=0)
            dm = -2.5 * np.log10(F[k0] / ref)
            dm = dm[np.isfinite(dm)]
            agn_rms = float(np.std(dm, ddof=1)) if dm.size >= 3 else None
        rows_n.append({"night": night, "n_frames": len(fr), "n_compact": int(nsrc),
                       "n_comps": len(comps),
                       "comp_rms_median": float(np.median(rms_c)) if rms_c else None,
                       "comp_rms_min": float(np.min(rms_c)) if rms_c else None,
                       "agn_minus_comps_rms": agn_rms,
                       "agn_peak_max": float(np.nanmax(PK[k0])),
                       "agn_saturated": int(np.nanmax(PK[k0]) >= veto)})
    write_table(con, "dw_zeroorder_nights", rows_n)
    med = [r["comp_rms_median"] for r in rows_n if r["comp_rms_median"] is not None]
    verdict = ("paragraph" if med and np.median(med) < ZO_RMS_MAX else "dropped")
    meta(con, zeroorder_utc=utcnow(),
         zeroorder_comp_rms_median=float(np.median(med)) if med else None,
         zeroorder_nights_with_comps=len(med), zeroorder_verdict=verdict)
    for r in rows_n:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print("verdict:", verdict, "median comp rms:", np.median(med) if med else None)
    return 0


# ===========================================================================
# lightcurves (DW-P55) and completeness (DW-P54) — NGC 5238 field stars
# ===========================================================================
LC_BANDS = ("G", "H", "L", "O", "R", "X")
LC_R = (13.0, 18.0)
LC_MIN_COVER = 0.7
LC_EXCLUDE_ARCMIN = 3.0           #: NGC 5238's body
LC_MAG_BINS = ((13.0, 15.0), (15.0, 16.5), (16.5, 18.0))
LC_STARS_PER_BIN = 3
LC_FMAX = 1.0 / 0.05              #: c/d (shortest injected period 0.05 d)


def _lc_matrix(con):
    """Stars x frames calibrated magnitudes for NGC 5238 (science frames,
    the six bands), with per-frame covariates."""
    from macro_dw import dwio
    fr = [dict(r) for r in con.execute(
        "SELECT d.obs_rowid, d.filter, d.readoutm, d.bjd_tdb_mid, d.airmass, d.exptime, "
        "s.fwhm_px, s.sky_adu, s.x_obj, s.y_obj, z.zp, z.apcor_6_to_13p5 "
        "FROM dw_frames d JOIN dw_disposition q USING(obs_rowid) "
        "JOIN dw_solve s USING(obs_rowid) JOIN dw_zp_frames z USING(obs_rowid) "
        "WHERE d.target='NGC5238' AND q.disposition='science' "
        "ORDER BY d.bjd_tdb_mid")]
    fr = [f for f in fr if f["filter"] in LC_BANDS]
    P = load_phot(con, [f["obs_rowid"] for f in fr])
    cat = _cat("NGC5238")
    km = dwio.rkey(dwio.MEAS_APER_ARCSEC)
    sep_deg = np.hypot((cat["ra"] - NGC5238_RADEC[0]) * np.cos(np.radians(NGC5238_RADEC[1])),
                       cat["dec"] - NGC5238_RADEC[1])
    col = {f["obs_rowid"]: j for j, f in enumerate(fr)}
    stars = np.unique(P["star"]).astype(int)
    stars = stars[(cat["r"][stars] > LC_R[0]) & (cat["r"][stars] < LC_R[1])
                  & (sep_deg[stars] * 60 > LC_EXCLUDE_ARCMIN)]
    row = {s_: i for i, s_ in enumerate(stars)}
    caps = {k: v["cap"] for k, v in detector_limits().items()}
    M = np.full((len(stars), len(fr)), np.nan)
    E = np.full_like(M, np.nan)
    for k in range(len(P["star"])):
        i = row.get(int(P["star"][k]))
        if i is None:
            continue
        j = col[int(P["obs_rowid"][k])]
        f, e = P["f" + km][k], P["e" + km][k]
        cap = caps[fr[j]["readoutm"]]
        if f > 0 and P["peak"][k] < cap and P["flag" + km][k] == 0:
            M[i, j] = (-2.5 * np.log10(f / fr[j]["exptime"])
                       + (fr[j]["apcor_6_to_13p5"] or 0) + fr[j]["zp"])
            E[i, j] = 1.0857 * e / f
    keep = np.isfinite(M).mean(axis=1) >= LC_MIN_COVER
    return fr, stars[keep], M[keep], E[keep], cat


def _lc_setup(con):
    fr, stars, M, E, cat = _lc_matrix(con)
    band = np.array([LC_BANDS.index(f["filter"]) for f in fr])
    t = np.array([f["bjd_tdb_mid"] for f in fr])
    # Per-star, per-band offsets removed -> residual matrix for Sys-Rem.
    Rm = M.copy()
    for b in np.unique(band):
        jb = band == b
        Rm[:, jb] -= np.nanmedian(M[:, jb], axis=1)[:, None]
    comps = core.sysrem(Rm, np.nan_to_num(E, nan=1.0) + 0.003)
    xo = np.array([f["x_obj"] for f in fr])
    yo = np.array([f["y_obj"] for f in fr])
    cov = np.column_stack([[f["airmass"] for f in fr], [f["fwhm_px"] for f in fr],
                           xo - np.median(xo), yo - np.median(yo),
                           [f["sky_adu"] for f in fr], *comps])
    cov = np.where(np.isfinite(cov), cov, np.nanmedian(cov, axis=0))
    return fr, stars, M, E, cat, band, t, cov


def cmd_lightcurves(args) -> int:
    """DW-P55: per-star fit of floating per-band offsets + covariates (+ <= 2
    Sys-Rem vectors); chi2_nu of the constant-star model per band with its
    dof, before and after the covariates (standing rule 1)."""
    con = connect()
    fr, stars, M, E, cat, band, t, cov = _lc_setup(con)
    rows = []
    for i in range(len(stars)):
        ok = np.isfinite(M[i])
        y, e, b = M[i][ok], E[i][ok], band[ok]
        e = np.hypot(e, 0.003)            # 3 mmag floor, stated
        X0 = core.nuisance_basis(b, np.empty((len(y), 0)))
        X1 = core.nuisance_basis(b, cov[ok])
        for name, X in (("offsets", X0), ("offsets+covariates", X1)):
            w = 1 / e ** 2
            beta, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
            r = (y - X @ beta) / e
            for bb in np.unique(b):
                m = b == bb
                dof = int(m.sum() - 1 - (X.shape[1] - len(np.unique(b))) * m.sum() / len(b))
                rows.append({"star": int(stars[i]), "r_mag": float(cat["r"][stars[i]]),
                             "model": name, "band": LC_BANDS[bb], "n": int(m.sum()),
                             "dof": max(dof, 1), "chi2nu": float(np.sum(r[m] ** 2) / max(dof, 1)),
                             "rms_mmag": float(1000 * np.std(y[m] - (X @ beta)[m]))})
    write_table(con, "dw_lc_chi2", rows)
    summ = []
    for name in ("offsets", "offsets+covariates"):
        for bb in LC_BANDS:
            v = [r for r in rows if r["model"] == name and r["band"] == bb]
            if v:
                summ.append({"model": name, "band": bb, "n_stars": len(v),
                             "chi2nu_median": float(np.median([r["chi2nu"] for r in v])),
                             "dof_median": float(np.median([r["dof"] for r in v])),
                             "rms_mmag_median": float(np.median([r["rms_mmag"] for r in v]))})
    write_table(con, "dw_lc_summary", summ)
    meta(con, lightcurves_utc=utcnow(), lc_n_stars=len(stars), lc_n_frames=len(fr),
         lc_covariates="airmass,fwhm,x_drift,y_drift,sky,sysrem1,sysrem2")
    for r in summ:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    return 0


def cmd_completeness(args) -> int:
    """DW-P54: sinusoids injected into real NGC 5238 field-star series on
    the real time stamps; recovery by the multiband covariate periodogram;
    90 % amplitude per period per magnitude bin with the signed matched-cell
    bias (standing rule 3)."""
    con = connect()
    fr, stars, M, E, cat, band, t, cov = _lc_setup(con)
    rng = np.random.default_rng(54)
    T = float(np.ptp(t))
    rows, summ = [], []
    rms = np.nanstd(M - np.nanmedian(M, axis=1)[:, None], axis=1)
    for lo, hi in LC_MAG_BINS:
        cand = [i for i in range(len(stars)) if lo <= cat["r"][stars[i]] < hi]
        cand = sorted(cand, key=lambda i: rms[i])
        pick = [cand[k] for k in np.linspace(0, len(cand) - 1,
                                             min(LC_STARS_PER_BIN, len(cand))).astype(int)
                ] if cand else []
        for i in pick:
            ok = np.isfinite(M[i])
            y, e, b, tt = M[i][ok], np.hypot(E[i][ok], 0.003), band[ok], t[ok]
            X = core.nuisance_basis(b, cov[ok])
            f = core.freq_grid(tt, LC_FMAX)
            pg = core.Periodogram(tt, 1 / e ** 2, X, f)
            # The star's own (null) residual series and its threshold.
            resid = (pg.residualise(y[:, None])[:, 0] / pg.sw)
            thr = core.shuffle_threshold(pg, resid, rng)
            own = float(pg.power(y[:, None]).max())
            for P_ in core.INJ_PERIODS_D:
                fi = 1.0 / P_
                phases = rng.uniform(0, 2 * np.pi, core.INJ_PHASES)
                for A in core.INJ_AMPS_MMAG:
                    Y = y[:, None] + (A / 1000) * np.sin(2 * np.pi * fi * tt[:, None]
                                                         + phases[None, :])
                    Pw = pg.power(Y)
                    ok_r = core.recovered(pg, Pw, np.full(Y.shape[1], fi), thr, T)
                    kpk = np.argmax(Pw, axis=0)
                    a_rec = pg.amplitude(Y, kpk) * 1000
                    rows.append({"r_bin": f"{lo}-{hi}", "star": int(stars[i]),
                                 "r_mag": float(cat["r"][stars[i]]), "period_d": float(P_),
                                 "amp_mmag": float(A), "n_trials": int(Y.shape[1]),
                                 "frac": float(ok_r.mean()),
                                 "signed_amp_bias": core.signed_cell_bias(
                                     a_rec, np.full(len(a_rec), A), ok_r),
                                 "star_own_power": own, "fap1_power": thr,
                                 "n_points": int(len(y))})
            print(f"bin {lo}-{hi} star {stars[i]} r={cat['r'][stars[i]]:.2f} "
                  f"N={len(y)} thr={thr:.3f} own={own:.3f}", flush=True)
    write_table(con, "dw_completeness_cells", rows)
    for rb in sorted({r["r_bin"] for r in rows}):
        for P_ in core.INJ_PERIODS_D:
            v = [r for r in rows if r["r_bin"] == rb and r["period_d"] == float(P_)]
            fr_ = np.array([[r["frac"] for r in v if r["amp_mmag"] == A] for A in core.INJ_AMPS_MMAG])
            fr_m = fr_.mean(axis=1)
            a90 = core.amp90(core.INJ_AMPS_MMAG, fr_m)
            # bias in the cell at (or just above) the 90 % amplitude
            k = int(np.argmax(core.INJ_AMPS_MMAG >= (a90 if np.isfinite(a90) else np.inf))) \
                if np.isfinite(a90) else len(core.INJ_AMPS_MMAG) - 1
            bias = [r["signed_amp_bias"] for r in v if r["amp_mmag"] == core.INJ_AMPS_MMAG[k]]
            summ.append({"r_bin": rb, "period_d": float(P_), "amp90_mmag": a90,
                         "cell_amp_mmag": float(core.INJ_AMPS_MMAG[k]),
                         "signed_bias_at_cell": float(np.nanmean(bias)) if bias else None})
    write_table(con, "dw_completeness", summ)
    meta(con, completeness_utc=utcnow())
    for r in summ:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    return 0


# ===========================================================================
# paper, status
# ===========================================================================
def cmd_paper(args) -> int:
    from macro_dw import paper_dw
    radec = {t: object_radec(t) for t in
             [r[0] for r in connect(read_only=True).execute(
                 "SELECT DISTINCT target FROM dw_frames")]}
    for p in paper_dw.build(DB, MANIFEST, radec):
        print(p)
    return 0


def cmd_status(args) -> int:
    con = connect(read_only=True)
    for t in ("dw_frames", "dw_flatprep", "dw_solve", "dw_zp_frames", "dw_disposition",
              "dw_stacks", "dw_halpha", "dw_sersic", "dw_depth", "dw_completeness"):
        try:
            print(t, con.execute(f"SELECT count(*) FROM {t}").fetchone()[0])
        except sqlite3.OperationalError:
            print(t, "-")
    return 0



def _html_table(con, sql: str, digits: int = 4) -> str:
    cur = con.execute(sql)
    cols = [d[0] for d in cur.description]
    out = ["<table><tr>" + "".join(f"<th>{c}</th>" for c in cols) + "</tr>"]
    for r in cur.fetchall():
        cells = []
        for v in r:
            if isinstance(v, float):
                v = f"{v:.{digits}g}"
            cells.append(f"<td>{'' if v is None else v}</td>")
        out.append("<tr>" + "".join(cells) + "</tr>")
    return "\n".join(out + ["</table>"])


def cmd_report(args) -> int:
    """The evidence page docs/DwarfGalaxy_AGN_Survey/dw_paper.html, every
    table drawn from dwarf.sqlite (Question -> Evidence -> Decision)."""
    import shutil
    con = connect(read_only=True)
    m = dict(con.execute("SELECT key, value FROM dw_build_meta").fetchall())
    out = REPO / "docs" / "DwarfGalaxy_AGN_Survey" / "dw_paper.html"
    fig = REPO / "manuscripts" / "DwarfGalaxy_AGN_Survey" / "figures" / "dw_fig4_halpha.png"
    if fig.exists():
        shutil.copy(fig, out.parent / fig.name)
    sec = []
    def S(q, ev, dec):
        sec.append(f"<h2>{q}</h2>{ev}<p><b>Decision.</b> {dec}</p>")
    S("Which frames are science? (DW-P01, DW-P03, DW-P32)",
      _html_table(con, "SELECT disposition, count(*) n FROM dw_disposition GROUP BY 1 ORDER BY 2 DESC")
      + _html_table(con, "SELECT qc_reasons, count(*) n FROM dw_disposition WHERE qc_reasons IS NOT NULL GROUP BY 1"),
      "Only 'science' frames enter stacks; every rejection is pipeline-emitted.")
    S("Is the night-sky superflat flat to 0.3 % of sky? (DW-P11, DW-P12)",
      _html_table(con, "SELECT filter, test, n_build, n_test, rms, rms_noplane FROM dw_flat_validation")
      + _html_table(con, "SELECT * FROM dw_gradient_moon"),
      "L: held-out large-scale residual below 0.3 % (rms column, half splits).")
    S("Zero points and colour terms (DW-P2-ensemble-zp, DW-P21, DW-P37)",
      _html_table(con, "SELECT * FROM dw_zp_colour") + _html_table(con, "SELECT * FROM dw_ps1_check")
      + _html_table(con, "SELECT zmag_flag, count(*) FROM dw_zmag_qc GROUP BY 1"),
      "Every frame is tied to REFCAT2 r with a fitted (g-r) colour term; ZMAG is QC only.")
    S("The NGC 5238 line-flux scale",
      _html_table(con, "SELECT * FROM dw_ha_cal") + _html_table(con, "SELECT * FROM dw_ha_cal_nights"),
      "Calibrates throughput at NGC 5238's Hα wavelength only; band edges stay unknown.")
    S("Hα detections and limits (DW-P36)",
      _html_table(con, "SELECT field, band, dlam_A, n_H, k, k_basis, net, sigma_total, snr, verdict, "
                       "flux_5238scale, flux_lim_3sig, F_pred, pred_over_lim, informative_limit FROM dw_halpha")
      + (f'<p><img src="{fig.name}" width="480"></p>' if fig.exists() else ""),
      f"Venue trigger: {m.get('venue_trigger_n')} in-band objects detected or informative -> {m.get('venue')}.")
    for t, q in (("dw_sersic", "Structure (DW-P35)"), ("dw_depth", "Depth (DW-P34)"),
                 ("dw_completeness", "NGC 5238 completeness (DW-P54)"),
                 ("dw_lc_summary", "Detrending chi2 per band (DW-P55)"),
                 ("dw_zeroorder_nights", "NGC 5548 zero orders (DW-P4y)")):
        try:
            S(q, _html_table(con, f"SELECT * FROM {t}"), "See the manuscript and the ledger note.")
        except sqlite3.OperationalError:
            pass
    html = ("<!doctype html><meta charset='utf-8'><title>Dwarf-Galaxy Hα — evidence</title>"
            "<style>body{font:14px system-ui;max-width:1100px;margin:2em auto;padding:0 1em}"
            "table{border-collapse:collapse;font-size:12px;margin:.5em 0}td,th{border:1px solid #ccc;"
            "padding:2px 6px}</style>"
            f"<h1>Dwarf-Galaxy Hα paper — evidence</h1><p>Generated {utcnow()} by "
            f"pipeline/scripts/run_dw_paper.py report from products/dwarf/dwarf.sqlite "
            f"({m.get('code', DW_CODE_VERSION)}).</p>" + "\n".join(sec))
    out.write_text(html)
    print(out)
    return 0


# ===========================================================================
# broadha  (DW-P4x: one-night broad-Hα triage of NGC 5548, slot '6')
# ===========================================================================
#: Pre-declared in committee/work/dwarf/PREDECLARED.md (DW-P4x rule).
Z5548 = 0.01717
LAM_HA_OBS = 6562.8 * (1 + Z5548)
LAM_O3_OBS = 5006.8 * (1 + Z5548)
CONT_WIN_OBS = ((6420.0, 6480.0), (6900.0, 6980.0))
HA_WIN_OBS = (6480.0, 6900.0)
P4X_MAX_REL_ERR = 0.02
TRACE_HALF, BG_IN, BG_OUT = 4, 8, 14


def _p4x_one(task: dict) -> dict:
    """Trace, extract and self-calibrate one slot-'6' frame; EW(Hα)."""
    from scipy.ndimage import median_filter
    from macro_dw import dwio
    out = {"obs_rowid": task["obs_rowid"], "night": task["night"],
           "jd": task["jd"], "ok": 0}
    try:
        raw, h = dwio.load(task["path"])
        img = raw.astype(float)
        ny, nx = img.shape
        prof = np.median(img[:, 1500:2300], axis=1)
        y0 = int(np.argmax(prof[200:-200])) + 200
        cut = img[y0 - 120:y0 + 121]
        narrow = cut - median_filter(cut, size=(31, 1))
        xs = np.arange(1000, 2800, 10)
        ys = []
        for x in xs:
            col = narrow[:, x - 5:x + 6].mean(axis=1)
            ys.append(np.argmax(col) + y0 - 120)
        ys = np.array(ys, float)
        A = np.c_[np.ones(len(xs)), xs]
        keep = np.ones(len(xs), bool)
        for _ in range(4):
            c, *_ = np.linalg.lstsq(A[keep], ys[keep], rcond=None)
            res = ys - A @ c
            keep = np.abs(res) < max(3 * core.robust_sigma(res[keep]), 1.5)
        out["trace_rms"] = float(np.std(res[keep]))
        out["trace_slope"] = float(c[1])
        cols = np.arange(800, 3000)
        spec = np.zeros(len(cols))
        for k, x in enumerate(cols):
            yc = c[0] + c[1] * x
            yi = int(round(yc))
            ap = img[yi - TRACE_HALF:yi + TRACE_HALF + 1, x].sum()
            bg = np.r_[img[yi - BG_OUT:yi - BG_IN, x], img[yi + BG_IN + 1:yi + BG_OUT + 1, x]]
            spec[k] = ap - (2 * TRACE_HALF + 1) * np.median(bg)
        sm = median_filter(spec, 5)
        cont = median_filter(sm, 151)
        line = sm - cont
        # Hα: strongest excess in the red half; [O III]: strongest excess
        # 350-800 px blueward of it.
        red = (cols > 1900) & (cols < 2600)
        kha = np.nonzero(red)[0][np.argmax(line[red])]
        blue = (cols > cols[kha] - 800) & (cols < cols[kha] - 350)
        ko3 = np.nonzero(blue)[0][np.argmax(line[blue])]
        disp = (LAM_HA_OBS - LAM_O3_OBS) / (cols[kha] - cols[ko3])
        lam = LAM_HA_OBS + disp * (cols - cols[kha])
        win = [(lam >= a) & (lam <= b) for a, b in CONT_WIN_OBS]
        cx = np.r_[lam[win[0]], lam[win[1]]]
        cy = np.r_[sm[win[0]], sm[win[1]]]
        cc = np.polyfit(cx, cy, 1)
        m = (lam >= HA_WIN_OBS[0]) & (lam <= HA_WIN_OBS[1])
        fc = np.polyval(cc, lam[m])
        ew_obs = float(np.sum((sm[m] - fc) / fc) * abs(disp))
        out.update(ok=1, x_ha=int(cols[kha]), x_o3=int(cols[ko3]), disp_A_px=float(disp),
                   ew_rest_A=ew_obs / (1 + Z5548),
                   cont_at_ha=float(np.polyval(cc, LAM_HA_OBS)),
                   peak_ha_snr=float(line[kha] / core.robust_sigma(line)))
    except Exception as ex:                                  # pragma: no cover
        out["error"] = repr(ex)[:200]
    return out


def cmd_broadha(args) -> int:
    from multiprocessing import Pool
    con = connect()
    tasks = [dict(r) for r in con.execute(
        "SELECT d.* FROM dw_frames d JOIN dw_disposition q USING(obs_rowid) "
        "WHERE d.target='NGC 5548' AND d.filter='6' AND q.disposition='spectrum' "
        "ORDER BY d.jd")]
    with Pool(args.workers) as pool:
        res = pool.map(_p4x_one, tasks)
    # A frame passes extraction if its self-calibrated dispersion agrees
    # with the frame-set median to 5 % (a wrong line identification moves
    # it by far more) and Hα stands >= 5 sigma above the local continuum.
    d = np.array([r.get("disp_A_px", np.nan) for r in res])
    dmed = float(np.nanmedian(d))
    for r in res:
        r["pass"] = int(r["ok"] and abs(r["disp_A_px"] / dmed - 1) < 0.05
                        and r["peak_ha_snr"] >= 5)
    write_table(con, "dw_p4x_frames", [{k: r.get(k) for k in
        ("obs_rowid", "night", "jd", "ok", "pass", "trace_rms", "trace_slope",
         "x_ha", "x_o3", "disp_A_px", "ew_rest_A", "cont_at_ha", "peak_ha_snr",
         "error")} for r in res])
    nights = {}
    for r in res:
        if r["pass"]:
            nights.setdefault(r["night"], []).append(r)
    best = sorted(nights, key=lambda n: (-len(nights[n]), n))[0]
    ew = np.array([r["ew_rest_A"] for r in nights[best]])
    rel = float(np.std(ew, ddof=1) / np.sqrt(len(ew)) / np.mean(ew))
    verdict = "paragraph" if rel < P4X_MAX_REL_ERR else "dropped"
    summary = {"night": best, "n_frames": len(ew), "ew_mean_A": float(np.mean(ew)),
               "ew_frame_scatter_A": float(np.std(ew, ddof=1)),
               "ew_frame_scatter_pct": float(100 * np.std(ew, ddof=1) / np.mean(ew)),
               "nightly_mean_rel_err_pct": 100 * rel, "disp_median_A_px": dmed,
               "n_frames_pass_all": int(sum(r["pass"] for r in res)),
               "n_frames_all": len(res), "verdict": verdict}
    write_table(con, "dw_p4x_summary", [summary])
    meta(con, p4x_utc=utcnow(), p4x_verdict=verdict)
    print(summary)
    return 0


# ===========================================================================
# CLI
# ===========================================================================
def main(argv=None) -> int:
    cmds = {k[4:]: v for k, v in globals().items()
            if k.startswith("cmd_") and callable(v)}
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in cmds:
        q = sub.add_parser(name)
        q.add_argument("--workers", type=int, default=10)
        q.add_argument("--limit", type=int, default=0)
        q.add_argument("--refresh", action="store_true")
    args = ap.parse_args(argv)
    return cmds[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
