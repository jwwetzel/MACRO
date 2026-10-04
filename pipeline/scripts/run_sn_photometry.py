#!/usr/bin/env python
"""run_sn_photometry.py — the SN 2023ixf validation-and-limits release.

One resumable CLI for everything downstream of Gate 0 (SYNTHESIS §4
SN2023ixf_LightCurve, strategy §10).  Pure decision logic lives in
``macro_sn.snphot`` (unit-tested in ``pipeline/tests/test_sn_phot.py``);
frame I/O in ``macro_sn.snio``; the paper's figures and numbers in
``macro_sn.paper_sn``.  Everything lands in ONE database:

    products/sn/sn2023ixf.sqlite

and the manifest is opened READ-ONLY.

SUBCOMMANDS, in dependency order (each rebuilds its own tables)
---------------------------------------------------------------
    fetch       external inputs, cached with retrieval dates under
                products/sn/external/: REFCAT2 cone (VizieR J/ApJ/867/105),
                Li et al. 2025 photometry (VizieR J/A+A/703/A168), ZTF alert
                photometry (ALeRCE API), PS1 r effective width (SVO FPS)
    frames      the worklist: every direct-image frame of the Gate 0 freeze
    measure     astrometry + forced photometry per frame (parallel,
                resumable; re-invoke until 'pending 0')        [SN-S4-resolve]
    crosswalk   colour-term regression per filter code + MaxIm->pyscope
                crosswalk                                [SN-S1-filter-crosswalk]
    calibrate   REFCAT2 ensemble calibration, held-out check stars,
                scintillation-aware error model           [SN-S4-ensemble-cal]
    templates   the template table                      [SN-S5-template-table]
    photometry  two-regime SN photometry + overlap test     [SN-S5-photometry]
    peak        the literature peak epoch                  [SN-S5b-peak-epoch]
    excessgate  predicted flash-phase H-alpha excess vs systematic [SN-S6-0]
    flash       (H - '1') differential colour at the flash epochs [SN-S6a]
    residuals   RLMT minus the published world photometry  [SN-S7b]
    variability nightly-mean limits, whole-night bootstrap, injections [SN-S8]
    latetime    late-time stacks, injection-defined 5-sigma limits [SN-S9]
    release     machine-readable tables + README          [SN-S10-release]
    paper       five figures + numbers.tex              [SN-figures, SN-draft]
    status      progress (read-only)

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/run_sn_photometry.py fetch
    $PY pipeline/scripts/run_sn_photometry.py frames
    $PY pipeline/scripts/run_sn_photometry.py measure --workers 10
    ...
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_sn import gate0 as g0                     # noqa: E402
from macro_sn import snphot as sp                    # noqa: E402

REPO = PIPELINE_ROOT.parent
MANIFEST = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
OUT = REPO / "products" / "sn"
DB = OUT / "sn2023ixf.sqlite"
EXT = OUT / "external"

SN_CODE_VERSION = "SN-PHOT v1.0 (2026-10-04)"

#: External sources — each is a public, citable product; the retrieval date
#: is stored beside the cached copy (strategy §9: "cache with retrieval
#: dates").  AAVSO is deliberately absent: its site refuses scripted access.
REFCAT_URL = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"
              "-source=J/ApJ/867/105/refcat2&-c=210.910675+54.31165"
              "&-c.r=30&-c.u=arcmin&-out.max=50000"
              "&-out=RA_ICRS,DE_ICRS,Plx,e_Plx,pmRA,e_pmRA,pmDE,e_pmDE,Gmag,"
              "gmag,e_gmag,rmag,e_rmag,imag,e_imag,rp1,r1,dupvar&rmag=<19")
LI25_URL = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/703/A168/table1.dat"
LI25_README = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/703/A168/ReadMe"
ALERCE_URL = "https://api.alerce.online/ztf/v1/objects/ZTF23aaklqou/lightcurve"
SVO_URL = "https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=PAN-STARRS/PS1.r"


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
    con.execute("CREATE TABLE IF NOT EXISTS sn_build_meta "
                "(key TEXT PRIMARY KEY, value TEXT)")
    for k, v in kw.items():
        con.execute("INSERT OR REPLACE INTO sn_build_meta VALUES (?, ?)",
                    (k, str(v)))


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


# ===========================================================================
# fetch
# ===========================================================================
def _get(url: str, dest: Path) -> int:
    req = urllib.request.Request(url, headers={"User-Agent": "MACRO-RLMT/1.0"})
    with urllib.request.urlopen(req, timeout=300) as r:
        data = r.read()
    dest.write_bytes(data)
    return len(data)


def cmd_fetch(args) -> int:
    EXT.mkdir(parents=True, exist_ok=True)
    got = {}
    for name, url, fn in (("refcat2", REFCAT_URL, "refcat2_m101.tsv"),
                          ("li25", LI25_URL, "li2025_table1.dat"),
                          ("li25_readme", LI25_README, "li2025_ReadMe.txt"),
                          ("ztf_alerce", ALERCE_URL, "ztf_alerce_lc.json"),
                          ("svo_ps1r", SVO_URL, "svo_ps1_r.xml")):
        dest = EXT / fn
        if dest.exists() and not args.refresh:
            print(f"cached  {fn}")
            continue
        n = _get(url, dest)
        got[name] = {"url": url, "file": fn, "bytes": n, "retrieved": utcnow()}
        print(f"fetched {fn} ({n} bytes)")
    log = EXT / "retrieval.json"
    old = json.loads(log.read_text()) if log.exists() else {}
    old.update(got)
    log.write_text(json.dumps(old, indent=1))
    return 0


def load_refcat() -> dict:
    """The cached REFCAT2 cone as numpy columns (row order = star_id)."""
    lines = [ln for ln in (EXT / "refcat2_m101.tsv").read_text().splitlines()
             if ln and not ln.startswith("#")]
    head = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[3:]]
    cols = {}
    for k, name in enumerate(head):
        vals = []
        for r in rows:
            try:
                vals.append(float(r[k]))
            except (ValueError, IndexError):
                vals.append(np.nan)
        cols[name] = np.array(vals)
    cols["star_id"] = np.arange(len(rows))
    return cols


# ===========================================================================
# frames
# ===========================================================================
FRAMES_DDL = """
CREATE TABLE IF NOT EXISTS sn_frames (
    obs_rowid INTEGER PRIMARY KEY, path TEXT, tree TEXT, night TEXT,
    filter TEXT, camera TEXT, readoutm TEXT, exptime REAL, bjd_tdb REAL,
    phase_d REAL, airmass REAL, epoch_role TEXT, band_role TEXT,
    saturation_class TEXT, sn_peak_census REAL, pltsolvd INTEGER,
    focpos REAL, mech_epoch TEXT, gain REAL, clip_adu REAL,
    egain REAL, flipstat TEXT, ccd_temp REAL, instrume TEXT, naxis1 INTEGER,
    naxis2 INTEGER, lin_cap_adu REAL,
    status TEXT DEFAULT 'pending', error TEXT,
    wcs_source TEXT, n_det INTEGER, n_match INTEGER, wcs_rms_arcsec REAL,
    wcs_header TEXT, fwhm_px REAL, bkg_level REAL, bkg_rms REAL,
    calib_recipe TEXT, sn_x REAL, sn_y REAL, sn_flux REAL, sn_err REAL,
    sn_peak REAL, sn_flag INTEGER, measured_at TEXT
);
CREATE TABLE IF NOT EXISTS sn_star_phot (
    obs_rowid INTEGER, star_id INTEGER, x REAL, y REAL,
    flux REAL, err REAL, peak REAL, flag INTEGER,
    PRIMARY KEY (obs_rowid, star_id)
);
"""


def cmd_frames(args) -> int:
    """Copy the worklist out of the Gate 0 freeze.  Every non-dispersed
    frame is measured; which of them may feed which product is decided
    later, by rule, from the columns copied here."""
    m = manifest()
    rows = m.execute("""
        SELECT s.obs_rowid, s.path, s.tree, s.night, s.filter, f.camera,
               s.readoutm, s.exptime, t.bjd_tdb, s.phase_d, f.airmass,
               s.epoch_role, s.band_role, c.saturation_class,
               c.peak_adu, f.pltsolvd, coalesce(f.focpos, f.focuspos),
               me.mech_epoch, f.egain, f.hdr_gain_num, f.flipstat,
               f.ccd_temp, f.instrume, f.naxis1, f.naxis2
        FROM sn_g0_frames s
        JOIN frames f ON f.obs_rowid = s.obs_rowid
        LEFT JOIN frame_times t ON t.obs_rowid = s.obs_rowid
        LEFT JOIN sn_g0_census c ON c.obs_rowid = s.obs_rowid
        LEFT JOIN frame_mech_epoch me ON me.obs_rowid = s.obs_rowid
        WHERE s.dispersion_class != 'dispersed'""").fetchall()
    clips = {r["mode"]: r["clip_adu"] or r["hard_max_adu"] for r in
             m.execute("SELECT mode, clip_adu, hard_max_adu "
                       "FROM s2_ceiling_modes")}
    # Gains and linearity caps are S2 MEASUREMENTS (detector_params):
    # per EGAIN epoch where S2 resolved one (the SN campaign is the 1.054
    # epoch: K = 1.072 +/- 0.015 e-/ADU, cap 1,800 ADU), else per mode.
    dp = {(r[0], r[1]): r[2] for r in m.execute(
        "SELECT era_group, quantity, value FROM detector_params")}
    m.close()
    con = connect()
    con.executescript(FRAMES_DDL)
    n_new = 0
    with con:
        for r in rows:
            mode, egain = r["readoutm"], r[18]
            ekey = f"{r['camera']} {mode} e{egain:.3f}" if egain else None
            gain = (dp.get((ekey, "gain_e_per_adu"))
                    or dp.get((mode, "gain_e_per_adu"))
                    or egain or r[19] or 1.0)
            cap = (dp.get((ekey, "linearity_cap_adu"))
                   or dp.get((mode, "linearity_cap_adu")))
            cur = con.execute(
                "INSERT OR IGNORE INTO sn_frames (obs_rowid, path, tree, "
                "night, filter, camera, readoutm, exptime, bjd_tdb, phase_d,"
                " airmass, epoch_role, band_role, saturation_class, "
                "sn_peak_census, pltsolvd, focpos, mech_epoch, gain, "
                "clip_adu, egain, flipstat, ccd_temp, instrume, naxis1, "
                "naxis2, lin_cap_adu) VALUES "
                "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                tuple(r[:18]) + (float(gain), clips.get(mode), egain,
                                 r["flipstat"], r["ccd_temp"], r["instrume"],
                                 r["naxis1"], r["naxis2"], cap))
            n_new += cur.rowcount
        meta(con, frames_built_at=utcnow(), code_version=SN_CODE_VERSION)
    tot = con.execute("SELECT count(*) FROM sn_frames").fetchone()[0]
    print(f"worklist: {tot} frames ({n_new} new)")
    return 0


# ===========================================================================
# measure  (SN-S4-resolve lives here: every frame gets a refined WCS)
# ===========================================================================
_CAL = None
_CAT = None


def _init_worker():
    global _CAL, _CAT
    from macro_sn import snio
    _CAL = snio.load_calib()
    _CAT = load_refcat()


#: Nominal plate scales (arcsec/px) per camera, used ONLY to project the
#: catalogue for blind matching; the fitted WCS measures the real one.
#: PROVENANCE: the PinPoint CD matrices of solved frames of each camera
#: (AC4040 0.540"/px; QHY600 2x2 0.90"/px) and the iKon's 0.81"/px (TE.F8).
NOMINAL_SCALE = {"AC4040": 0.540, "QHY600": 0.90, "iKon": 0.81}


def _measure(task: dict) -> dict:
    from macro_sn import snio
    t0 = time.time()
    out = {"obs_rowid": task["obs_rowid"], "status": "failed", "stars": None}
    try:
        raw, h = snio.load_data(task["path"])
        img, mask, recipe = snio.calibrate(raw, task["filter"], task["camera"],
                                           _CAL)
        bp = snio.badpix_mask(task, raw.shape)
        if bp is not None:
            mask |= bp
            img[bp] = 0.0
            recipe += "; badpix mask (rlmt_diagnostics.badpix)"
        else:
            recipe += "; no bad-pixel mask for this camera/orientation"

        sub, bkg, obj = snio.detect(img, mask)
        cat = _CAT
        good_cat = np.isfinite(cat["rmag"])
        cra, cdec, cmag = (cat["RA_ICRS"][good_cat], cat["DE_ICRS"][good_cat],
                           cat["rmag"][good_cat])
        w0 = snio.header_wcs(h) if task["pltsolvd"] == 1 else None
        src = "header"
        w, n, rms = (None, 0, None)
        if w0 is not None:
            w, n, rms = snio.refine_wcs(w0, obj, cra, cdec, img.shape)
        if w is None or not sp.wcs_acceptable(n, rms, NOMINAL_SCALE.get(
                task["camera"], 0.54)):
            src = "solved"
            pt = snio.pointing(h)
            if pt is None:
                pt = (g0.SN_RA_DEG, g0.SN_DEC_DEG)
            wb = None
            if task.get("ref_cd") is not None:
                # Known-rotation camera: translation vote first (robust on
                # thin frames), triangle matching only if it fails.
                obj_v = obj
                if len(obj) < 60:
                    # Thin frame (cloud or a 0.5-4 s exposure): a 3-sigma
                    # source list for the vote; the photometry is forced
                    # at catalogue positions and does not use it.
                    import sep
                    obj_v = sep.extract(sub, 3.0, err=bkg.globalrms,
                                        minarea=5, mask=mask)
                wb = snio.offset_vote_solve(obj_v, cra, cdec, task["ref_cd"],
                                            pt[0], pt[1], img.shape)
                if wb is not None:
                    w, n, rms = snio.refine_wcs(wb, obj_v, cra, cdec,
                                                img.shape)
                    src = "solved-vote"
                    if not sp.wcs_acceptable(n, rms, 0.54):
                        # Fixed-CD translation: one free vector, so fewer
                        # stars carry it (sp.WCS_MIN_MATCH_SHIFT).
                        w, n, rms = snio.shift_refine(wb, obj_v, cra, cdec,
                                                      img.shape)
                        src = "solved-shift"
                        if w is not None and sp.shift_acceptable(n, rms, 0.54):
                            pass
                        else:
                            w = None
            if w is None or (src != "solved-shift" and not sp.wcs_acceptable(
                    n, rms, NOMINAL_SCALE.get(task["camera"], 0.54))):
                wb = snio.blind_solve(obj, cra, cdec, cmag, pt[0], pt[1],
                                      NOMINAL_SCALE.get(task["camera"], 0.54),
                                      img.shape)
                src = "solved"
                if wb is not None:
                    w, n, rms = snio.refine_wcs(wb, obj, cra, cdec, img.shape)
        if w is None or not (sp.wcs_acceptable(n, rms, NOMINAL_SCALE.get(
                task["camera"], 0.54)) or (src == "solved-shift"
                                          and sp.shift_acceptable(n, rms, 0.54))):
            out.update(error=f"no acceptable WCS (n={n}, rms={rms})",
                       n_det=len(obj), wcs_source="failed")
            return out
        clip = task["clip_adu"] or 65535.0
        fwhm = snio.fwhm_of(obj, 0.8 * clip)
        if fwhm is None:
            fwhm = float(task.get("hdr_fwhm") or 4.0)
        ny, nx = img.shape
        x, y = w.all_world2pix(cat["RA_ICRS"], cat["DE_ICRS"], 0)
        e = snio.EDGE_PX + snio.ANNULUS_FWHM[1] * fwhm
        ins = np.isfinite(x) & (x > e) & (x < nx - e) & (y > e) & (y < ny - e) \
            & np.isfinite(cat["rmag"]) & (cat["rmag"] < 18.5)
        ids = cat["star_id"][ins]
        fx, fy = x[ins], y[ins]
        flux, err, flag = snio.forced_phot(sub, bkg, mask, fx, fy, fwhm,
                                           task["gain"])
        peak = snio.native_peaks(raw, fx, fy, max(fwhm, 3.0))
        sx, sy = w.all_world2pix([g0.SN_RA_DEG], [g0.SN_DEC_DEG], 0)
        sx, sy = float(sx[0]), float(sy[0])
        if 0 < sx < nx and 0 < sy < ny:
            sf, se, sfl = snio.forced_phot(sub, bkg, mask, np.array([sx]),
                                           np.array([sy]), fwhm, task["gain"])
            spk = snio.native_peaks(raw, [sx], [sy], max(fwhm, 3.0))[0]
            sf, se, sfl = float(sf[0]), float(se[0]), int(sfl[0])
        else:
            sf = se = spk = None
            sfl = -1
        out.update(status="measured", wcs_source=src, n_det=len(obj),
                   n_match=n, wcs_rms_arcsec=rms,
                   wcs_header=w.to_header(relax=True).tostring(),
                   fwhm_px=fwhm, bkg_level=float(bkg.globalback),
                   bkg_rms=float(bkg.globalrms), calib_recipe=recipe,
                   sn_x=sx, sn_y=sy, sn_flux=sf, sn_err=se, sn_peak=spk,
                   sn_flag=sfl,
                   stars=list(zip(ids.tolist(), fx.tolist(), fy.tolist(),
                                  flux.tolist(), err.tolist(),
                                  [None if not np.isfinite(p) else float(p)
                                   for p in peak], flag.tolist())))
    except Exception as exc:                       # recorded, never swallowed
        out.update(error=f"{type(exc).__name__}: {exc}")
    out["measure_s"] = time.time() - t0
    return out


def cmd_measure(args) -> int:
    from multiprocessing import Pool
    con = connect()
    todo = [dict(r) for r in con.execute(
        "SELECT obs_rowid, path, filter, camera, pltsolvd, clip_adu, gain, "
        "instrume, naxis1, naxis2, readoutm, ccd_temp, flipstat "
        "FROM sn_frames WHERE status != 'measured' ORDER BY obs_rowid")]
    if args.limit:
        todo = todo[:args.limit]
    # The campaign camera's measured CD matrix: the median over its own
    # header-solved, accepted frames (one orientation all season).
    from astropy.io.fits import Header
    from astropy.wcs import WCS as _W
    cds = []
    for (hdr,) in con.execute(
            "SELECT wcs_header FROM sn_frames WHERE status='measured' AND "
            "camera='AC4040' AND wcs_source='header' AND epoch_role="
            "'campaign' LIMIT 200"):
        cds.append(_W(Header.fromstring(hdr)).pixel_scale_matrix)
    ref_cd = np.median(np.array(cds), axis=0).tolist() if cds else None
    for t in todo:
        t["ref_cd"] = ref_cd if t["camera"] == "AC4040" else None
    print(f"measuring {len(todo)} frames with {args.workers} workers")
    cols = ("status", "error", "wcs_source", "n_det", "n_match",
            "wcs_rms_arcsec", "wcs_header", "fwhm_px", "bkg_level", "bkg_rms",
            "calib_recipe", "sn_x", "sn_y", "sn_flux", "sn_err", "sn_peak",
            "sn_flag")
    done = 0
    with Pool(args.workers, initializer=_init_worker) as pool:
        for res in pool.imap_unordered(_measure, todo):
            with con:
                con.execute(
                    f"UPDATE sn_frames SET {', '.join(c + '=?' for c in cols)},"
                    f" measured_at=? WHERE obs_rowid=?",
                    tuple(res.get(c) for c in cols) + (utcnow(),
                                                       res["obs_rowid"]))
                con.execute("DELETE FROM sn_star_phot WHERE obs_rowid=?",
                            (res["obs_rowid"],))
                if res["stars"]:
                    con.executemany(
                        "INSERT INTO sn_star_phot VALUES (?,?,?,?,?,?,?,?)",
                        [(res["obs_rowid"],) + tuple(s) for s in res["stars"]])
            done += 1
            if done % 25 == 0 or res["status"] != "measured":
                print(f"  {done}/{len(todo)} {res['obs_rowid']} "
                      f"{res['status']} {res.get('wcs_source')} "
                      f"n={res.get('n_match')} {res.get('error') or ''}",
                      flush=True)
    with con:
        meta(con, measured_at=utcnow())
    return cmd_status(args)


def cmd_status(args) -> int:
    con = connect(read_only=True)
    for r in con.execute("SELECT status, wcs_source, count(*) FROM sn_frames "
                         "GROUP BY 1, 2"):
        print(f"  {r[0]:10s} {str(r[1]):8s} {r[2]}")
    pend = con.execute("SELECT count(*) FROM sn_frames "
                       "WHERE status != 'measured'").fetchone()[0]
    print(f"pending {pend}")
    return 0


COMMANDS = {}

# ===========================================================================
# calibrate + crosswalk  (SN-S4-ensemble-cal, SN-S1-filter-crosswalk)
# ===========================================================================
#: Frame rejection (strategy §4 Step 4): FWHM above this, fewer than
#: MIN_ENSEMBLE ensemble stars, or a transparency more than CLOUD_MAG below
#: the series median (the strategy's 0.5 mag cloud flag, applied to our own
#: zero point rather than to PinPoint's ZMAG).
MAX_FWHM_PX = 6.0
MIN_ENSEMBLE = 8
CLOUD_MAG = 0.5
#: Point screen: S/N below this is not used anywhere in the ensemble.
MIN_SNR = 5.0
#: Ensemble and check stars must have a median per-frame S/N of at least
#: this.  The supernova is measured at S/N 50-300; the zero point and the
#: error model that its errors inherit are set — and tested — on stars in
#: the same well-measured regime, not on faint stars whose extra scatter
#: (galaxy background, faint neighbours) the SN does not share.
ENS_MIN_SNR = 20.0

#: The calibration SETS: which frames share one ensemble.  A set is one
#: (camera, filter code) over one epoch group; the campaign codes are the
#: paper's, the template sets exist for the crosswalk and for template
#: scaling.
CAL_SETS = {
    ("campaign", "AC4040"): ("G", "R", "I", "H", "O", "1"),
    ("template_pre", "AC4040"): ("G", "R", "H", "O"),
    ("template_post", "QHY600"): ("g", "r", "i", "ha"),
    ("template_post", "iKon"): ("g", "r"),
}


#: Frames that are not images by this paper's rule (S2c-dispersed, or
#: S2c-indeterminate without a passed PSF check); filled by cmd_calibrate.
_NOT_IMAGING: set = set()


def not_imaging(con) -> set:
    m = manifest()
    disp = {r[0]: (r[1], r[2]) for r in m.execute(
        "SELECT obs_rowid, dispersion_class, epoch_role FROM sn_g0_census")}
    m.close()
    ok = {r[0]: bool(r[1]) for r in con.execute(
        "SELECT obs_rowid, passes FROM sn_psfcheck")}
    # The indeterminate rule is a CAMPAIGN-photometry rule; template frames
    # are judged by the template table's own star-PSF fits instead.
    return {o for o, (c, role) in disp.items() if c == "dispersed"
            or (c == "indeterminate" and role == "campaign"
                and not ok.get(o, False))}


def load_series(con, role, camera, code):
    """(frames, star ids, inst mag, sig_phot, usable mask) for one set."""
    fr = [dict(r) for r in con.execute("""
        SELECT obs_rowid, exptime, airmass, fwhm_px, night, bjd_tdb,
               lin_cap_adu, clip_adu, tree, phase_d
        FROM sn_frames WHERE status = 'measured' AND epoch_role = ?
          AND camera = ? AND filter = ? AND tree = 'rawimage'
        ORDER BY bjd_tdb""", (role, camera, code))]
    fr = [f for f in fr if f["obs_rowid"] not in _NOT_IMAGING]
    if not fr:
        return None
    idx = {f["obs_rowid"]: j for j, f in enumerate(fr)}
    rows = con.execute(f"""
        SELECT p.obs_rowid, p.star_id, p.flux, p.err, p.peak, p.flag, p.x, p.y
        FROM sn_star_phot p WHERE p.obs_rowid IN
        ({",".join(str(k) for k in idx)})""").fetchall()
    sids = np.array(sorted({r[1] for r in rows}))
    sidx = {s: i for i, s in enumerate(sids)}
    S, F = len(sids), len(fr)
    mag = np.full((S, F), np.nan)
    sig = np.full((S, F), np.nan)
    peak = np.full((S, F), np.nan)
    px = np.full((S, F), np.nan)
    py = np.full((S, F), np.nan)
    cap = np.array([f["lin_cap_adu"] or 0.8 * (f["clip_adu"] or 65535)
                    for f in fr])
    for o, sid, fl, er, pk, flg, x, y in rows:
        i, j = sidx[sid], idx[o]
        peak[i, j] = pk if pk is not None else np.nan
        px[i, j], py[i, j] = x, y
        if fl is None or er is None or fl <= 0 or er <= 0 or flg != 0:
            continue
        if fl / er < MIN_SNR:
            continue
        if pk is not None and pk >= cap[j]:
            continue                  # star-side linearity screen, per point
        mag[i, j] = -2.5 * np.log10(fl)
        sig[i, j] = sp.MAG_PER_REL * er / fl
    return fr, sids, mag, sig, peak, px, py


def _airmass(fr):
    x = np.array([f["airmass"] if f["airmass"] else np.nan for f in fr], float)
    x[~np.isfinite(x)] = np.nanmedian(x) if np.isfinite(x).any() else 1.2
    return x


CAL_DDL = """
DROP TABLE IF EXISTS sn_cal_sets;
CREATE TABLE sn_cal_sets (
    role TEXT, camera TEXT, code TEXT, n_frames INTEGER, n_frames_used INTEGER,
    n_ensemble INTEGER, n_check INTEGER, k2 REAL, scint_C REAL, floor_mag REAL,
    Z REAL, Z_err REAL, cterm REAL, cterm_err REAL, tie_band TEXT,
    tie_rms REAL, n_tie INTEGER, colour_lo REAL, colour_hi REAL,
    ens_chi2nu REAL, ens_dof INTEGER, check_chi2nu REAL, check_dof INTEGER,
    check_verdict TEXT, check_ext_offset REAL, check_ext_rms REAL,
    PRIMARY KEY (role, camera, code));
DROP TABLE IF EXISTS sn_cal_frames;
CREATE TABLE sn_cal_frames (
    obs_rowid INTEGER PRIMARY KEY, role TEXT, camera TEXT, code TEXT,
    used INTEGER, reject_reason TEXT, zp REAL, zp_err REAL, zp_nat REAL,
    n_ens INTEGER, scint_unit REAL, transparency REAL);
DROP TABLE IF EXISTS sn_cal_scint;
CREATE TABLE sn_cal_scint (
    role TEXT, camera TEXT, code TEXT, exptime REAL, n_points INTEGER,
    n_stars INTEGER, rms_obs REAL, rms_model REAL, ratio REAL,
    young_scint_mag REAL, fitted_scint_mag REAL, scint_var_frac REAL);
DROP TABLE IF EXISTS sn_crosswalk;
CREATE TABLE sn_crosswalk (
    role TEXT, camera TEXT, code TEXT, ps1_band TEXT, Z REAL, cterm REAL,
    cterm_err REAL, tie_rms REAL, n INTEGER, identified INTEGER);
"""

#: The PS1 band each code is TIED to once identified; the crosswalk
#: verifies the broadband entries (argmin tie rms) rather than assuming them.
TIE_BAND = {"G": "g", "R": "r", "I": "i", "g": "g", "r": "r", "i": "i",
            "H": "r", "O": "g", "1": "i", "ha": "r"}


def calibrate_set(con, cat, calmask, role, camera, code, fixed_C=None,
                  dflat=True):
    """Solve one set; returns (set row, frame rows, scint rows, xwalk rows)."""
    got = load_series(con, role, camera, code)
    if got is None:
        return None
    fr, sids, mag, sig, peak, px, py = got
    F = len(fr)
    colour = (cat["gmag"] - cat["imag"])[sids]
    is_check = np.array([sp.is_check_star(cat["RA_ICRS"][s], cat["DE_ICRS"][s])
                         for s in sids])
    eligible = calmask[sids]
    X = _airmass(fr)
    t = np.array([f["exptime"] for f in fr], float)
    su = sp.young_scintillation_mag(t, X)
    fwhm = np.array([f["fwhm_px"] or np.nan for f in fr])
    # frame screen 1: seeing (campaign camera only — the template cameras
    # have other pixel scales and are judged by their own median)
    lim = MAX_FWHM_PX if camera == "AC4040" else 1.5 * np.nanmedian(fwhm)
    use_f = np.isfinite(fwhm) & (fwhm <= lim)
    reason = np.where(use_f, "", f"fwhm>{lim:.1f}px").astype(object)
    # stars must be seen on >= 30% of used frames to join the ensemble
    seen = np.sum(np.isfinite(mag[:, use_f]), axis=1)
    med_snr = np.nanmedian(np.where(np.isfinite(sig), sp.MAG_PER_REL / sig,
                                    np.nan), axis=1)
    bright_enough = np.nan_to_num(med_snr) >= ENS_MIN_SNR
    is_comp = eligible & bright_enough & (seen >= max(3, 0.3 * use_f.sum()))
    is_check = is_check & eligible & bright_enough
    m_use = np.where(use_f[None, :], mag, np.nan)
    tie = TIE_BAND[code]
    sol = sp.solve_band(m_use, sig, su, is_comp, is_check, colour,
                        cat[f"{tie}mag"][sids], X, fixed_C=fixed_C)
    # Delta-flat (strategy §5.3) from ENSEMBLE residuals of well-measured
    # points; applied to every star (and later to the SN) by position.
    dcoef = None
    if dflat and camera == "AC4040":
        okd = is_comp[:, None] & ~is_check[:, None] & np.isfinite(
            sol["residual"]) & ~sol["clipped"] & np.isfinite(px)
        dcoef = sp.fit_delta_flat(px[okd], py[okd], sol["residual"][okd],
                                  1.0 / sol["sig"][okd] ** 2)
        corr = sp.eval_delta_flat(dcoef, np.nan_to_num(px, nan=2048.0),
                                  np.nan_to_num(py, nan=2048.0))
        mag = mag - corr.reshape(mag.shape)
        m_use = np.where(use_f[None, :], mag, np.nan)
        sol = sp.solve_band(m_use, sig, su, is_comp, is_check, colour,
                            cat[f"{tie}mag"][sids], X, fixed_C=fixed_C)
    # frame screens 2 and 3: ensemble size and cloud
    n_ens = sol["n_star_used"]
    # m_inst = M + ZP and m_inst = -2.5 log(rate * t), so a frame's
    # exposure-normalised zero point is ZP + 2.5 log t (larger = fainter).
    transp = sol["zp"] + 2.5 * np.log10(t)
    med_t = np.nanmedian(transp[use_f & (n_ens >= MIN_ENSEMBLE)])
    for j in range(F):
        if not use_f[j]:
            continue
        if n_ens[j] < MIN_ENSEMBLE:
            use_f[j], reason[j] = False, f"<{MIN_ENSEMBLE} ensemble stars"
        elif transp[j] > med_t + CLOUD_MAG:
            use_f[j], reason[j] = False, f"cloud (>{CLOUD_MAG} mag)"
    if (~use_f & (reason != "") & ~np.char.startswith(
            reason.astype(str), "fwhm")).any():
        m_use = np.where(use_f[None, :], mag, np.nan)
        sol = sp.solve_band(m_use, sig, su, is_comp, is_check, colour,
                            cat[f"{tie}mag"][sids], X, fixed_C=fixed_C)
        n_ens = sol["n_star_used"]
    # held-out test
    c2, dof, res, e, g_rows, f_cols = sp.check_star_test(
        sol["mcorr"], sol["sig"], np.where(use_f, sol["zp"], np.nan),
        sol["zp_err"], is_check)
    # ensemble chi2 for reference (should be ~1 by construction)
    ok = is_comp[:, None] & ~is_check[:, None] & np.isfinite(sol["residual"]) \
        & ~sol["clipped"] & use_f[None, :]
    ens_c2 = float(np.sum((sol["residual"][ok] / sol["sig"][ok]) ** 2))
    ens_dof = int(ok.sum() - np.any(ok, axis=0).sum() - np.any(ok, axis=1).sum())
    # check stars against the CATALOGUE (external): offset and rms
    Mchk = np.array([np.average(sol["mcorr"][s][np.isfinite(sol["mcorr"][s])
                                               & use_f]
                                - sol["zp"][np.isfinite(sol["mcorr"][s]) & use_f])
                     if (np.isfinite(sol["mcorr"][s]) & use_f).sum() >= 3
                     else np.nan for s in range(len(sids))])
    ext = Mchk - cat[f"{tie}mag"][sids] - sol["Z"] - sol["cterm"] * (
        colour - sp.COLOUR_PIVOT)
    chk_ok = is_check & eligible & np.isfinite(ext)
    tie_rows = sol["tie_rows"]
    tie_res = (sol["mean_mag"] - cat[f"{tie}mag"][sids] - sol["Z"]
               - sol["cterm"] * (colour - sp.COLOUR_PIVOT))[tie_rows]
    lo, hi = (np.percentile(colour[tie_rows], [5, 95]) if len(tie_rows)
              else (np.nan, np.nan))
    setrow = (role, camera, code, F, int(use_f.sum()),
              int((is_comp & ~is_check).sum()), int(len(np.unique(g_rows))),
              sol["k2"], sol["C"], sol["floor"], sol["Z"], sol["Z_err"],
              sol["cterm"], sol["cterm_err"], tie,
              float(np.sqrt(np.mean(tie_res ** 2))) if len(tie_res) else None,
              int(len(tie_rows)), float(lo), float(hi),
              ens_c2 / ens_dof if ens_dof > 0 else None, ens_dof,
              c2, dof, sp.chi2_verdict(c2),
              float(np.median(ext[chk_ok])) if chk_ok.any() else None,
              float(1.4826 * np.median(np.abs(ext[chk_ok]
                                              - np.median(ext[chk_ok]))))
              if chk_ok.any() else None)
    frows = [(fr[j]["obs_rowid"], role, camera, code, int(use_f[j]),
              reason[j] or None, float(sol["zp"][j]), float(sol["zp_err"][j]),
              float(sol["zp"][j] + sol["Z"]), int(n_ens[j]), float(su[j]),
              float(transp[j])) for j in range(F)]
    # scintillation match: bright check residuals grouped by exposure time
    srows = []
    if len(res):
        tt = t[f_cols]
        sphot = sig[g_rows, f_cols]
        for te in np.unique(np.round(tt, 1)):
            m = np.round(tt, 1) == te
            if m.sum() < 5:
                continue
            nst = len(np.unique(g_rows[m]))
            # residuals about a star's mean shrink by (n-1)/n; correct.
            n_per = np.array([np.sum(g_rows == k) for k in g_rows[m]])
            corr = np.sqrt(n_per / np.maximum(n_per - 1, 1))
            rms_o = float(np.sqrt(np.mean((res[m] * corr) ** 2)))
            rms_m = float(np.sqrt(np.mean(e[m] ** 2)))
            # The match statistic is the NORMALISED rms, sqrt(<(r/sigma)^2>):
            # a plain rms ratio is dominated by the largest-error points.
            nrm = float(np.sqrt(np.mean((res[m] * corr / e[m]) ** 2)))
            srows.append((role, camera, code, float(te), int(m.sum()), nst,
                          rms_o, rms_m, nrm,
                          float(np.median(su[f_cols][m])),
                          float(sol["C"] * np.median(su[f_cols][m])),
                          float(np.mean((sol["C"] * su[f_cols][m]) ** 2)
                                / np.mean(e[m] ** 2))))
    # crosswalk: the same ensemble means tied to each PS1 band in turn
    xrows = []
    rms_by = {}
    for b in ("g", "r", "i"):
        okb = is_comp & ~is_check & np.isfinite(sol["mean_mag"]) \
            & np.isfinite(cat[f"{b}mag"][sids])
        if okb.sum() < 8:
            continue
        a, c, cov, kept = sp.robust_wls_line(
            (colour - sp.COLOUR_PIVOT)[okb],
            (sol["mean_mag"] - cat[f"{b}mag"][sids])[okb],
            np.full(okb.sum(), 0.02))
        rr = ((sol["mean_mag"] - cat[f"{b}mag"][sids])[okb]
              - a - c * (colour - sp.COLOUR_PIVOT)[okb])[kept]
        rms_by[b] = float(np.sqrt(np.mean(rr ** 2)))
        xrows.append([role, camera, code, b, a, c, float(np.sqrt(cov[1, 1])),
                      rms_by[b], int(kept.sum()), 0])
    if xrows:
        # Identification = the PS1 band the code responds to most nearly
        # 1:1, i.e. the smallest |colour term|.  (The tie rms cannot
        # discriminate g from i: with c = g - i as the colour, fitting
        # against g with slope c_g is the same fit as against i with slope
        # c_g - 1.)
        best = min(xrows, key=lambda x: abs(x[5]))[3]
        for x in xrows:
            x[-1] = int(x[3] == best)
    okp = is_comp[:, None] & ~is_check[:, None] & np.isfinite(
        sol["residual"]) & ~sol["clipped"] & use_f[None, :]
    nobs = okp.sum(axis=1)
    pooled = (sol["residual"][okp], sig[okp],
              np.broadcast_to(su[None, :], sig.shape)[okp],
              np.broadcast_to(((nobs - 1) / np.maximum(nobs, 1))[:, None],
                              sig.shape)[okp])
    dfl = None
    if dcoef is not None:
        gx, gy = np.meshgrid(np.linspace(200, 3896, 20), np.linspace(200, 3896, 20))
        surf = sp.eval_delta_flat(dcoef, gx.ravel(), gy.ravel())
        dfl = (role, camera, code, *[float(c) for c in dcoef],
               float(1e3 * (surf.max() - surf.min())))
    return setrow, frows, srows, [tuple(x) for x in xrows], pooled, dfl


#: Iterations of the pooled scintillation fit (C converges in 2-3).
POOL_ITER = 3


def cmd_calibrate(args) -> int:
    """Two passes.  (1) The campaign broadband codes G, R, I are solved with
    a common scintillation scale C, fitted POOLED over their ensemble
    residuals (sp.fit_scint_pooled) and iterated; (2) every other set is
    solved with C held at that value (a 32-256 s frame cannot measure C),
    fitting only its own floor."""
    global _NOT_IMAGING
    con = connect()
    _NOT_IMAGING = not_imaging(con)
    cat = load_refcat()
    calmask = sp.calib_star_mask(cat)
    C = 1.5                                   # start: Osborn+15 median
    bb = ("G", "R", "I")
    for it in range(POOL_ITER):
        mats = [calibrate_set(con, cat, calmask, "campaign", "AC4040", c,
                              fixed_C=C) for c in bb]
        r = np.concatenate([m[4][0] for m in mats])
        sph = np.concatenate([m[4][1] for m in mats])
        su = np.concatenate([m[4][2] for m in mats])
        kf = np.concatenate([m[4][3] for m in mats])
        bi = np.concatenate([np.full(len(m[4][0]), k) for k, m in enumerate(mats)])
        C, floors = sp.fit_scint_pooled(r, sph, su, bi,
                                        dof_factor=np.clip(kf, 0.1, 1))
        print(f"pooled pass {it + 1}: C = {C:.3f}, floors (mmag) = "
              f"{', '.join(f'{b} {1e3 * floors[k]:.1f}' for k, b in enumerate(bb))}")
    con.executescript(CAL_DDL + """
        DROP TABLE IF EXISTS sn_cal_dflat;
        CREATE TABLE sn_cal_dflat (role TEXT, camera TEXT, code TEXT,
            c0 REAL, cu REAL, cv REAL, cuu REAL, cuv REAL, cvv REAL,
            ptp_mmag REAL);""")
    for (role, camera), codes in CAL_SETS.items():
        for code in codes:
            try:
                got = calibrate_set(con, cat, calmask, role, camera, code,
                                    fixed_C=C)
            except (np.linalg.LinAlgError, ValueError) as exc:
                # Recorded, not hidden: a set with too few calibrators.
                print(f"{role:13s} {camera:7s} {code:2s} NOT CALIBRATED: "
                      f"{type(exc).__name__}: {exc}")
                continue
            if got is None:
                continue
            setrow, frows, srows, xrows, _pooled, dfl = got
            with con:
                con.execute(f"INSERT INTO sn_cal_sets VALUES "
                            f"({','.join('?' * len(setrow))})", setrow)
                con.executemany("INSERT INTO sn_cal_frames VALUES "
                                "(?,?,?,?,?,?,?,?,?,?,?,?)", frows)
                con.executemany("INSERT INTO sn_cal_scint VALUES "
                                "(?,?,?,?,?,?,?,?,?,?,?,?)", srows)
                con.executemany("INSERT INTO sn_crosswalk VALUES "
                                "(?,?,?,?,?,?,?,?,?,?)", xrows)
                if dfl:
                    con.execute("INSERT INTO sn_cal_dflat VALUES "
                                "(?,?,?,?,?,?,?,?,?,?)", dfl)
            print(f"{role:13s} {camera:7s} {code:2s} frames {setrow[4]}/"
                  f"{setrow[3]} ens {setrow[5]} chk {setrow[6]} "
                  f"c={setrow[12]:+.3f}+-{setrow[13]:.3f} k2={setrow[7]:+.3f} "
                  f"floor={setrow[9]*1e3:.1f}mmag "
                  f"tie_rms={setrow[15] or float('nan'):.3f} "
                  f"check chi2nu={setrow[21]:.2f}/{setrow[22]} "
                  f"[{setrow[23]}] dflat={dfl[-1] if dfl else 0:.1f}mmag")
    with con:
        meta(con, calibrated_at=utcnow(), code_version=SN_CODE_VERSION,
             scint_C=C)
    return 0


COMMANDS["calibrate"] = cmd_calibrate



# ===========================================================================
# photometry  (SN-S5-photometry: the bright, aperture regime + nightly means)
# ===========================================================================
PHOT_DDL = """
DROP TABLE IF EXISTS sn_phot;
CREATE TABLE sn_phot (
    obs_rowid INTEGER PRIMARY KEY, night TEXT, code TEXT, bjd_tdb REAL,
    phase_d REAL, exptime REAL, airmass REAL, census_class TEXT,
    sn_peak_adu REAL, lin_cap_adu REAL, cal_used INTEGER, usable INTEGER,
    exclusion TEXT, m_inst REAL, dflat_mag REAL, m_nat REAL, sig_phot REAL,
    sig_scint REAL, sig_floor REAL, sig_zp REAL, sig_tot REAL,
    colour_gi REAL, colour_in_range INTEGER, m_ps1 REAL, fwhm_px REAL);
DROP TABLE IF EXISTS sn_nightly;
CREATE TABLE sn_nightly (
    night TEXT, code TEXT, n INTEGER, bjd_tdb REAL, phase_d REAL,
    m_nat REAL, err REAL, chi2nu REAL, m_ps1 REAL, colour_gi REAL,
    colour_in_range INTEGER, PRIMARY KEY (night, code));
"""

BROAD = ("G", "R", "I")
PS1_OF = {"G": "g", "R": "r", "I": "i"}


def cmd_photometry(args) -> int:
    con = connect()
    m = manifest()
    census = {r[0]: r[1] for r in m.execute(
        "SELECT obs_rowid, saturation_class FROM sn_g0_census")}
    # S2c's CURRENT verdict, read at photometry time: a broadband frame
    # re-measured as dispersed (SN-G0d) leaves the photometry here even if
    # it was measured before the verdict changed.
    disp = {r[0]: r[1] for r in m.execute(
        "SELECT obs_rowid, dispersion_class FROM sn_g0_census")}
    psf_ok = {r[0]: bool(r[1]) for r in con.execute(
        "SELECT obs_rowid, passes FROM sn_psfcheck")}
    m.close()
    sets = {r["code"]: dict(r) for r in con.execute(
        "SELECT * FROM sn_cal_sets WHERE role='campaign' AND camera='AC4040'")}
    dfl = {r["code"]: [r["c0"], r["cu"], r["cv"], r["cuu"], r["cuv"], r["cvv"]]
           for r in con.execute("SELECT * FROM sn_cal_dflat "
                                "WHERE role='campaign'")}
    # SN frames are photometry only in exposure groups where the held-out
    # check stars VALIDATE the error model (normalised rms inside
    # sp.SCINT_MATCH_WINDOW with >= sp.SCINT_MIN_POINTS points) — the SN-S5
    # acceptance turned into a rule, fixed before any comparison with the
    # published photometry.
    validated = {}
    for r in con.execute("SELECT code, exptime, n_points, ratio FROM "
                         "sn_cal_scint WHERE role = 'campaign'"):
        validated[(r[0], round(r[1], 1))] = (
            r[2] >= sp.SCINT_MIN_POINTS
            and sp.SCINT_MATCH_WINDOW[0] <= r[3] <= sp.SCINT_MATCH_WINDOW[1])
    rows = con.execute("""
        SELECT f.*, c.used, c.zp_nat, c.zp_err, c.scint_unit, c.reject_reason
        FROM sn_frames f LEFT JOIN sn_cal_frames c USING (obs_rowid)
        WHERE f.epoch_role = 'campaign' AND f.filter IN
              ('G','R','I','H','O','1')
        ORDER BY f.bjd_tdb""").fetchall()
    out = []
    for r in rows:
        code = r["filter"]
        st = sets.get(code)
        cls = census.get(r["obs_rowid"])
        excl = []
        if r["tree"] != g0.SCIENCE_TREE:
            excl.append("not science tree")
        if cls not in g0.USABLE_CLASSES:
            excl.append(f"census {cls}")
        if disp.get(r["obs_rowid"]) == "dispersed":
            excl.append("S2c: dispersed")
        elif disp.get(r["obs_rowid"]) == "indeterminate" and not psf_ok.get(
                r["obs_rowid"], False):
            excl.append("S2c indeterminate and own PSF check not passed")
        if r["status"] != "measured":
            excl.append("no WCS")
        elif not r["used"]:
            excl.append(f"calibration: {r['reject_reason'] or 'not in a set'}")
        if r["sn_flux"] is None or (r["sn_flux"] or 0) <= 0:
            excl.append("no SN flux")
        if (r["sn_peak"] is not None and r["lin_cap_adu"]
                and r["sn_peak"] >= r["lin_cap_adu"]):
            excl.append("own peak >= cap")
        if not validated.get((code, round(r["exptime"], 1)), False):
            excl.append("error model not validated at this exposure")
        ok = not excl and st is not None
        rec = dict(obs_rowid=r["obs_rowid"], night=r["night"], code=code,
                   bjd_tdb=r["bjd_tdb"], phase_d=r["phase_d"],
                   exptime=r["exptime"], airmass=r["airmass"],
                   census_class=cls, sn_peak_adu=r["sn_peak"],
                   lin_cap_adu=r["lin_cap_adu"], cal_used=r["used"],
                   usable=int(ok), exclusion="; ".join(excl) or None,
                   fwhm_px=r["fwhm_px"])
        if (r["sn_flux"] or 0) > 0 and r["zp_nat"] is not None and st:
            mi = -2.5 * np.log10(r["sn_flux"])
            dm = float(sp.eval_delta_flat(dfl[code], r["sn_x"], r["sn_y"])) \
                if code in dfl else 0.0
            sph = sp.MAG_PER_REL * r["sn_err"] / r["sn_flux"]
            ssc = st["scint_C"] * r["scint_unit"]
            rec.update(m_inst=mi, dflat_mag=dm, m_nat=mi - dm - r["zp_nat"],
                       sig_phot=sph, sig_scint=ssc, sig_floor=st["floor_mag"],
                       sig_zp=r["zp_err"],
                       sig_tot=float(np.sqrt(sph ** 2 + ssc ** 2
                                             + st["floor_mag"] ** 2
                                             + r["zp_err"] ** 2)))
        out.append(rec)
    # The SN's OWN per-frame error floor.  The star floor (22-39 mmag) is
    # set by S/N 20-100 field stars; the SN is measured at S/N 100-300 near
    # one detector position, and its nights hold 2-6 frames each.  Its floor
    # is therefore fitted by maximum likelihood to its own WITHIN-NIGHT
    # residuals (night means free, dof-corrected) — the strategy's
    # "scatter-based" nightly error, pooled over nights so that two-frame
    # nights do not each carry a two-point variance.  Photon, scintillation
    # and zero-point terms stay explicit; only the floor is replaced.
    sn_floor = {}
    for code in ("G", "R", "I", "H", "O", "1"):
        res, sph, ssc, kf = [], [], [], []
        by = {}
        for rec in out:
            if rec["code"] == code and rec["usable"]:
                by.setdefault(rec["night"], []).append(rec)
        for rs in by.values():
            if len(rs) < 2:
                continue
            v = np.array([x["m_nat"] for x in rs])
            e = np.array([x["sig_tot"] for x in rs])
            mu = np.average(v, weights=1 / e ** 2)
            res += list(v - mu)
            sph += [np.hypot(x["sig_phot"], x["sig_zp"]) for x in rs]
            ssc += [x["sig_scint"] for x in rs]
            kf += [(len(rs) - 1) / len(rs)] * len(rs)
        if len(res) < 10:
            continue
        fl = sp.fit_floor(np.array(res), np.array(sph), np.array(ssc),
                          dof_factor=np.array(kf))
        sn_floor[code] = (fl, len(res), len(by))
        for rec in out:
            if rec["code"] == code and rec.get("m_nat") is not None:
                rec["sig_floor"] = fl
                rec["sig_tot"] = float(np.sqrt(rec["sig_phot"] ** 2
                                               + rec["sig_scint"] ** 2
                                               + rec["sig_zp"] ** 2 + fl ** 2))
    # nightly natural means, then the SN's own (g - i) per night
    nightly = {}
    for code in ("G", "R", "I", "H", "O", "1"):
        by = {}
        for rec in out:
            if rec["code"] == code and rec["usable"]:
                by.setdefault(rec["night"], []).append(rec)
        for night, rs in by.items():
            mu, err, c2, n = sp.weighted_mean([x["m_nat"] for x in rs],
                                              [x["sig_tot"] for x in rs])
            w = 1 / np.array([x["sig_tot"] for x in rs]) ** 2
            nightly[(night, code)] = dict(
                night=night, code=code, n=n, m_nat=mu, err=err, chi2nu=c2,
                bjd_tdb=float(np.average([x["bjd_tdb"] for x in rs], weights=w)),
                phase_d=float(np.average([x["phase_d"] for x in rs], weights=w)))
    cg, ci = sets["G"]["cterm"], sets["I"]["cterm"]
    cn = sorted({k[0] for k in nightly if k[1] == "G"}
                & {k[0] for k in nightly if k[1] == "I"})
    ph = np.array([0.5 * (nightly[(n_, "G")]["phase_d"]
                          + nightly[(n_, "I")]["phase_d"]) for n_ in cn])
    col = np.array([sp.colour_from_natural(nightly[(n_, "G")]["m_nat"],
                                           nightly[(n_, "I")]["m_nat"], cg, ci)
                    for n_ in cn])

    def colour_at(phase):
        return float(np.interp(phase, ph, col))       # clamps at the ends

    for rec in out:
        if rec.get("m_nat") is None or rec["code"] not in BROAD:
            continue
        st = sets[rec["code"]]
        c = colour_at(rec["phase_d"])
        rec.update(colour_gi=c,
                   colour_in_range=int(st["colour_lo"] <= c <= st["colour_hi"]),
                   m_ps1=float(sp.natural_to_ps1(rec["m_nat"], c, st["cterm"],
                                                 st["k2"], rec["airmass"] or 1.2)))
    for (night, code), nt in nightly.items():
        if code in BROAD:
            st = sets[code]
            c = colour_at(nt["phase_d"])
            nt.update(colour_gi=c,
                      colour_in_range=int(st["colour_lo"] <= c <= st["colour_hi"]),
                      m_ps1=float(sp.natural_to_ps1(nt["m_nat"], c, st["cterm"])))
    con.executescript(PHOT_DDL)
    cols = ["obs_rowid", "night", "code", "bjd_tdb", "phase_d", "exptime",
            "airmass", "census_class", "sn_peak_adu", "lin_cap_adu",
            "cal_used", "usable", "exclusion", "m_inst", "dflat_mag", "m_nat",
            "sig_phot", "sig_scint", "sig_floor", "sig_zp", "sig_tot",
            "colour_gi", "colour_in_range", "m_ps1", "fwhm_px"]
    ncols = ["night", "code", "n", "bjd_tdb", "phase_d", "m_nat", "err",
             "chi2nu", "m_ps1", "colour_gi", "colour_in_range"]
    with con:
        con.executemany(f"INSERT INTO sn_phot VALUES ({','.join('?' * len(cols))})",
                        [tuple(r.get(c) for c in cols) for r in out])
        con.executemany(f"INSERT INTO sn_nightly VALUES ({','.join('?' * len(ncols))})",
                        [tuple(n.get(c) for c in ncols) for n in nightly.values()])
        meta(con, photometry_at=utcnow())
        con.execute("DROP TABLE IF EXISTS sn_sn_floor")
        con.execute("CREATE TABLE sn_sn_floor (code TEXT PRIMARY KEY, "
                    "floor_mag REAL, n_points INTEGER, n_nights INTEGER, "
                    "star_floor_mag REAL)")
        con.executemany("INSERT INTO sn_sn_floor VALUES (?,?,?,?,?)",
                        [(k, v[0], v[1], v[2], sets[k]["floor_mag"])
                         for k, v in sn_floor.items()])
    for k, v in sn_floor.items():
        print(f"SN floor {k}: {v[0]*1e3:.1f} mmag from {v[1]} frames on "
              f"{v[2]} nights (star floor {sets[k]['floor_mag']*1e3:.1f})")
    for code in ("G", "R", "I", "H", "1"):
        n_use = sum(1 for r in out if r["code"] == code and r["usable"])
        n_n = sum(1 for k in nightly if k[1] == code)
        print(f"{code}: {n_use} usable frames on {n_n} nights")
    return 0


COMMANDS["photometry"] = cmd_photometry

# ===========================================================================
# overlap  (SN-S5: the template-subtracted regime and the overlap test)
# ===========================================================================
#: Grid half-size for subtraction cutouts (px at 0.54"/px): 601 x 601 px =
#: 5.4' square, enough REFCAT2 stars to fix the science/template flux ratio
#: locally (the template cameras have no flat field, so a whole-frame zero
#: point would carry their vignetting into the ratio).
SUB_HALF = 300
#: Template epochs per broadband code: the 2026-03-21/22 QHY600 'pinwheel
#: galaxy' epoch, the strategy's artifact-rejection workhorse (§3.4); the
#: crosswalk identifies G/R/I with the same PS1 bands as pyscope g/r/i.
TEMPLATE_OF = {"G": ("QHY600", "g"), "R": ("QHY600", "r"),
               "I": ("QHY600", "i")}


def template_nights(con, camera, tcode):
    """The in-focus epochs of a template code (sn_template_table)."""
    return tuple(r[0] for r in con.execute(
        "SELECT night FROM sn_template_table WHERE camera = ? AND filter = ? "
        "AND in_focus = 1 ORDER BY night", (camera, tcode)))


def _grid_image(path, wcs_header, filt, camera):
    """Calibrate (AC4040) or background-subtract (others) one frame and
    resample it onto the SN-centred grid."""
    from astropy.io.fits import Header
    from astropy.wcs import WCS as _W
    from macro_sn import snio
    import sep
    raw, _h = snio.load_data(path)
    img, mask, _ = snio.calibrate(raw, filt, camera, _CAL)
    img = np.ascontiguousarray(img, dtype=np.float32)
    bkg = sep.Background(img, mask=mask, bw=64, bh=64)
    sub = img - bkg.back()
    sub[mask] = np.nan
    w = _W(Header.fromstring(wcs_header))
    grid = snio.sky_grid(g0.SN_RA_DEG, g0.SN_DEC_DEG, SUB_HALF)
    return snio.resample_to_grid(sub, w, grid, SUB_HALF), grid


def _grid_stars(grid, cat, n):
    """REFCAT2 stars on the grid usable for flux ratios (isolated, r<17)."""
    x, y = grid.all_world2pix(cat["RA_ICRS"], cat["DE_ICRS"], 0)
    rp1 = np.nan_to_num(cat["rp1"], nan=99.9)
    ok = (x > 15) & (x < n - 15) & (y > 15) & (y < n - 15) \
        & (cat["rmag"] < 17.5) & (cat["rmag"] > 12.0) & (rp1 >= 5.0)
    return x[ok], y[ok]


def _ratio(a, b, xs, ys, r):
    """Median flux ratio a/b over stars (aperture r, local annulus)."""
    import sep
    fa, _, _ = sep.sum_circle(np.nan_to_num(a).astype(np.float64), xs, ys, r,
                              bkgann=(r + 4, r + 10))
    fb, _, _ = sep.sum_circle(np.nan_to_num(b).astype(np.float64), xs, ys, r,
                              bkgann=(r + 4, r + 10))
    ok = (fa > 0) & (fb > 0)
    if ok.sum() < 3:
        return None, int(ok.sum())
    return float(np.median(fa[ok] / fb[ok])), int(ok.sum())


def build_template(con, cat, code):
    """Median stack of the template epoch on the grid, frames normalised to
    the first by local star ratios.  Returns (stack, fwhm_px, n_frames)."""
    from macro_sn import snio
    cam, tcode = TEMPLATE_OF[code]
    nights = template_nights(con, cam, tcode)
    rows = con.execute(f"""SELECT path, wcs_header, filter, camera FROM sn_frames
        WHERE status='measured' AND camera=? AND filter=? AND night IN
        ({",".join("?" * len(nights))}) ORDER BY obs_rowid""",
                       (cam, tcode, *nights)).fetchall()
    ims = []
    for r in rows:
        g, grid = _grid_image(r[0], r[1], r[2], r[3])
        ims.append(g)
    n = 2 * SUB_HALF + 1
    xs, ys = _grid_stars(grid, cat, n)
    ref = ims[0]
    f0 = snio.grid_fwhm(ref, xs, ys) or 5.0
    scaled = []
    for im in ims:
        k, _n = _ratio(ref, im, xs, ys, 2.0 * f0)
        if k is not None:
            scaled.append(im * k)
    stack = np.nanmedian(np.array(scaled), axis=0)
    return stack, snio.grid_fwhm(stack, xs, ys), len(scaled), (xs, ys)


def subtract_and_measure(sci, f_sci, tpl, f_tpl, stars, r_mult=1.5):
    """PSF-match, scale, subtract; return the two SN fluxes on the matched
    science image: (aperture with annulus, difference image), and the
    matched FWHM, star count and ratio."""
    from macro_sn import snio
    import sep
    xs, ys = stars
    a, b, fm = snio.gauss_match(sci, f_sci, tpl, f_tpl)
    k, nst = _ratio(a, b, xs, ys, 2.0 * fm)
    if k is None:
        return None
    d = snio.plane_background(a - k * b)
    c = float(SUB_HALF)
    r = r_mult * fm
    fa, _, _ = sep.sum_circle(np.nan_to_num(a).astype(np.float64),
                              np.array([c]), np.array([c]), r,
                              bkgann=(4 * fm, 6 * fm))
    fd = snio.aper_sum(d, c, c, r)
    return float(fa[0]), float(fd), fm, nst, k, d


OVERLAP_DDL = """
DROP TABLE IF EXISTS sn_overlap;
CREATE TABLE sn_overlap (obs_rowid INTEGER PRIMARY KEY, code TEXT,
    night TEXT, phase_d REAL, m_nat REAL, f_aper REAL, f_diff REAL,
    dm REAL, fwhm_matched REAL, n_stars INTEGER, ratio REAL);
DROP TABLE IF EXISTS sn_templates;
CREATE TABLE sn_templates (code TEXT PRIMARY KEY, camera TEXT, tcode TEXT,
    nights TEXT, n_frames INTEGER, fwhm_px REAL, n_grid_stars INTEGER);
"""


def _overlap_one(task):
    from macro_sn import snio  # noqa: F401
    try:
        sci, _g = _grid_image(task["path"], task["wcs_header"], task["code"],
                              "AC4040")
        f_sci = __import__("macro_sn.snio", fromlist=["x"]).grid_fwhm(
            sci, *task["stars"]) or task["fwhm_px"]
        got = subtract_and_measure(sci, f_sci, _TPL[task["code"]][0],
                                   _TPL[task["code"]][1], task["stars"])
        if got is None:
            return None
        fa, fd, fm, nst, k, _d = got
        return (task["obs_rowid"], task["code"], task["night"], task["phase_d"],
                task["m_nat"], fa, fd, -2.5 * np.log10(fa / fd) if fa > 0
                and fd > 0 else None, fm, nst, k)
    except Exception as exc:
        print("overlap fail", task["obs_rowid"], exc, flush=True)
        return None


_TPL = {}


def _init_overlap(tpl):
    global _CAL, _TPL
    from macro_sn import snio
    _CAL = snio.load_calib()
    _TPL = tpl


def cmd_overlap(args) -> int:
    from multiprocessing import Pool
    global _CAL
    from macro_sn import snio
    _CAL = snio.load_calib()
    con = connect()
    cat = load_refcat()
    tpl, trow = {}, []
    for code in BROAD:
        stack, f, n, stars = build_template(con, cat, code)
        tpl[code] = (stack, f, stars)
        cam, tcode = TEMPLATE_OF[code]
        nights = template_nights(con, cam, tcode)
        trow.append((code, cam, tcode, ",".join(nights), n, f, len(stars[0])))
        print(f"template {code}: {n} {cam} '{tcode}' frames, FWHM {f:.2f} px, "
              f"{len(stars[0])} grid stars")
    tasks = [dict(obs_rowid=r[0], path=r[1], wcs_header=r[2], code=r[3],
                  night=r[4], phase_d=r[5], m_nat=r[6], fwhm_px=r[7],
                  stars=tpl[r[3]][2])
             for r in con.execute("""SELECT f.obs_rowid, f.path, f.wcs_header,
                    p.code, p.night, p.phase_d, p.m_nat, f.fwhm_px
                FROM sn_phot p JOIN sn_frames f USING (obs_rowid)
                WHERE p.usable = 1 AND p.code IN ('G','R','I')""")]
    with Pool(args.workers, initializer=_init_overlap, initargs=(tpl,)) as pool:
        res = [r for r in pool.imap_unordered(_overlap_one, tasks) if r]
    con.executescript(OVERLAP_DDL)
    with con:
        con.executemany("INSERT INTO sn_overlap VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                        res)
        con.executemany("INSERT INTO sn_templates VALUES (?,?,?,?,?,?,?)", trow)
        meta(con, overlap_at=utcnow())
    for code in BROAD:
        dm = np.array([r[7] for r in res if r[1] == code and r[7] is not None])
        if len(dm):
            print(f"{code}: N={len(dm)} median dm={np.median(dm):+.4f} "
                  f"robust sd={1.4826*np.median(np.abs(dm-np.median(dm))):.4f}")
    return 0


COMMANDS["overlap"] = cmd_overlap

# ===========================================================================
# templates  (SN-S5-template-table, TE.F8)
# ===========================================================================
#: A template epoch is IN FOCUS when its median FWHM (arcsec, from the
#: re-fitted WCS scale) is no more than this factor above the campaign's
#: median FWHM in the same AC4040 filter — or, for the other cameras, above
#: the best epoch of that camera and filter on this field.  Fixed before
#: the table was built.
IN_FOCUS_FACTOR = 1.3

TPL_DDL = """
DROP TABLE IF EXISTS sn_template_table;
CREATE TABLE sn_template_table (
    night TEXT, camera TEXT, filter TEXT, n_frames INTEGER, exptime REAL,
    mech_epoch TEXT, rotation_deg REAL, scale_arcsec REAL, fwhm_px REAL,
    fwhm_arcsec REAL, focpos REAL, focus_ref REAL, focus_offset REAL,
    in_focus INTEGER, role TEXT, PRIMARY KEY (night, camera, filter));
"""


def _star_fwhm(task):
    """Second-moment FWHM (px) at REFCAT2 positions, 12 < r < 15.5,
    isolated, in one frame (41 px stamps, so a defocused PSF fits)."""
    from astropy.io.fits import Header
    from astropy.wcs import WCS as _W
    from macro_sn import snio
    import sep
    raw, _h = snio.load_data(task["path"])
    img, mask, _ = snio.calibrate(raw, task["filter"], task["camera"], _CAL)
    img = np.ascontiguousarray(img, dtype=np.float32)
    sub = img - sep.Background(img, mask=mask, bw=64, bh=64).back()
    w = _W(Header.fromstring(task["wcs_header"]))
    cat = _CAT
    x, y = w.all_world2pix(cat["RA_ICRS"], cat["DE_ICRS"], 0)
    ok = (cat["rmag"] > 12) & (cat["rmag"] < 15.5) & \
        (np.nan_to_num(cat["rp1"], nan=99.9) >= 10) & np.isfinite(x)
    return snio.grid_fwhm(sub, x[ok], y[ok], half=20)


def cmd_templates(args) -> int:
    from astropy.io.fits import Header
    from astropy.wcs import WCS as _W
    con = connect()
    m = manifest()
    rot = {r[0]: r[1] for r in m.execute(
        "SELECT mech_epoch, rotation_deg FROM mech_epoch")}
    # The focus reference per (camera, filter): the median FOCPOS of that
    # camera's SCIENCE frames in the same filter across the whole archive
    # (for AC4040 that is dominated by the campaign itself).
    fref = {}
    for cam, filt in con.execute("SELECT DISTINCT camera, filter FROM sn_frames "
                                 "WHERE epoch_role LIKE 'template%'"):
        v = [r[0] for r in m.execute(
            "SELECT coalesce(focpos, focuspos) FROM frames WHERE camera = ? "
            "AND filter = ? AND is_canonical = 1 AND coalesce(focpos, "
            "focuspos) IS NOT NULL", (cam, filt))]
        fref[(cam, filt)] = float(np.median(v)) if v else None
    m.close()
    camp_fwhm = {}
    for filt, in con.execute("SELECT DISTINCT filter FROM sn_frames "
                             "WHERE epoch_role='campaign'"):
        v = [r[0] * 0.54 for r in con.execute(
            "SELECT fwhm_px FROM sn_frames WHERE epoch_role='campaign' AND "
            "filter=? AND fwhm_px IS NOT NULL AND tree='rawimage'", (filt,))]
        camp_fwhm[filt] = float(np.median(v)) if v else None
    from multiprocessing import Pool
    pool = Pool(args.workers, initializer=_init_worker)
    out = []
    for g in con.execute("""SELECT night, camera, filter, count(*),
                avg(exptime), mech_epoch, avg(focpos), epoch_role
            FROM sn_frames WHERE epoch_role LIKE 'template%'
            GROUP BY night, camera, filter ORDER BY night, filter""").fetchall():
        hs = [r[0] for r in con.execute(
            "SELECT wcs_header FROM sn_frames WHERE night=? AND camera=? AND "
            "filter=? AND status='measured'", g[:3])]
        fw = [v for v in pool.map(_star_fwhm, [dict(r) for r in con.execute(
            "SELECT path, wcs_header, filter, camera FROM sn_frames WHERE "
            "night=? AND camera=? AND filter=? AND status='measured'", g[:3])])
              if v is not None]
        scale = (float(np.median([np.sqrt(abs(np.linalg.det(
            _W(Header.fromstring(h)).pixel_scale_matrix))) * 3600 for h in hs]))
            if hs else None)
        fpx = float(np.median(fw)) if fw else None
        fas = fpx * scale if fpx and scale else None
        ref = fref.get((g[1], g[2]))
        out.append(dict(night=g[0], camera=g[1], filter=g[2], n=g[3],
                        exptime=g[4], mech=g[5], rot=rot.get(g[5]),
                        scale=scale, fpx=fpx, fas=fas, foc=g[6], fref=ref,
                        foff=(g[6] - ref) if (g[6] is not None and ref) else None,
                        role=g[7]))
    pool.close()
    best = {}
    for o in out:
        if o["fas"]:
            k = (o["camera"], o["filter"])
            best[k] = min(best.get(k, 99), o["fas"])
    rows = []
    for o in out:
        if o["camera"] == "AC4040" and camp_fwhm.get(o["filter"]):
            refw = camp_fwhm[o["filter"]]
        else:
            refw = best.get((o["camera"], o["filter"]))
        inf = int(bool(o["fas"] and refw and o["fas"] <= IN_FOCUS_FACTOR * refw))
        rows.append((o["night"], o["camera"], o["filter"], o["n"], o["exptime"],
                     o["mech"], o["rot"], o["scale"], o["fpx"], o["fas"],
                     o["foc"], o["fref"], o["foff"], inf, o["role"]))
    con.executescript(TPL_DDL)
    with con:
        con.executemany("INSERT INTO sn_template_table VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
        meta(con, templates_at=utcnow(), in_focus_factor=IN_FOCUS_FACTOR)
    for r in rows:
        print(f"{r[0]} {r[1]:7s} {r[2]:3s} n={r[3]:3d} {r[4]:6.0f}s "
              f"{r[5]} rot={r[6] if r[6] is None else round(r[6],1)} "
              f"FWHM={r[9] and round(r[9],2)}\" foc={r[10] and round(r[10])} "
              f"off={r[12] and round(r[12])} in_focus={r[13]}")
    return 0


COMMANDS["templates"] = cmd_templates

# ===========================================================================
# literature table + peak epoch  (SN-S5b-peak-epoch)
# ===========================================================================
def load_li25() -> list:
    """Li et al. (2025) photometry, VizieR J/A+A/703/A168 table1 (cached):
    rows of (mjd, phase, band, mag, err, telescope)."""
    out = []
    for ln in (EXT / "li2025_table1.dat").read_text().splitlines():
        if not ln.strip():
            continue
        out.append((float(ln[0:9]), float(ln[10:17]), ln[18].strip(),
                    float(ln[20:26]), float(ln[28:33]), ln[34:].strip()))
    return out


#: Half-width (d) of the window about the brightest published point in
#: which a weighted parabola locates each band's maximum.
PEAK_WINDOW_D = 5.0
PEAK_BOOT = 2000


def cmd_peak(args) -> int:
    from macro_sn import literature as lit
    con = connect()
    li = load_li25()
    rows = []
    rng = np.random.default_rng(20231)
    for b in ("g", "r", "i", "V"):
        d = np.array([(x[1], x[3], x[4]) for x in li if x[2] == b and x[1] < 40])
        k = np.argmin(d[:, 1])
        w = np.abs(d[:, 0] - d[k, 0]) <= PEAK_WINDOW_D
        t, m, e = d[w].T
        def tpk(tt, mm, ee):
            c = np.polyfit(tt, mm, 2, w=1 / ee)
            return -c[1] / (2 * c[0]) if c[0] > 0 else np.nan
        t0 = tpk(t, m, e)
        boots = [tpk(t, m + rng.normal(0, e), e) for _ in range(PEAK_BOOT)]
        boots = np.array(boots)[np.isfinite(boots)]
        lo, hi = np.percentile(boots, [16, 84]) if len(boots) else (np.nan, np.nan)
        # where the window's earliest point already sits at the maximum the
        # parabola is an extrapolation: then the published data say only
        # that the peak is at or before that epoch.
        before = bool(t0 < t.min()) if np.isfinite(t0) else True
        rows.append((b, float(t0), float(lo), float(hi), int(w.sum()),
                     float(d[k, 0]), float(d[k, 1]), int(before),
                     lit.LI25_REF + " (VizieR J/A+A/703/A168), weighted "
                     f"parabola within +/-{PEAK_WINDOW_D:g} d of the "
                     "brightest point"))
    first = {r[0]: r[1] for r in con.execute(
        "SELECT code, min(phase_d) FROM sn_phot WHERE usable = 1 "
        "AND code IN ('G','R','I') GROUP BY code")}
    con.executescript("""
        DROP TABLE IF EXISTS sn_peak;
        CREATE TABLE sn_peak (band TEXT PRIMARY KEY, t_peak_d REAL,
            t_lo REAL, t_hi REAL, n_points INTEGER, t_brightest REAL,
            m_brightest REAL, at_or_before_first INTEGER, basis TEXT);
        DROP TABLE IF EXISTS sn_literature;
        CREATE TABLE sn_literature (key TEXT PRIMARY KEY, value REAL,
            source TEXT);""")
    with con:
        con.executemany("INSERT INTO sn_peak VALUES (?,?,?,?,?,?,?,?,?)", rows)
        lits = [(f"peak_V_{k.split()[0]}", v, k) for k, v in
                lit.PEAK_V_DAYS.items()]
        lits += [(f"first_clean_{b}", v, "this work (sn_phot)")
                 for b, v in first.items()]
        lits += [("t0_mjd_li25", lit.LI25_T0_MJD, lit.LI25_REF),
                 ("f_halpha_max", lit.F_HALPHA_MAX_CGS, lit.F_HALPHA_REF),
                 ("f_halpha_max_phase", lit.F_HALPHA_MAX_PHASE_D, lit.F_HALPHA_REF),
                 ("f_halpha_fade", lit.F_HALPHA_FADE_FACTOR, lit.F_HALPHA_FADE_REF),
                 ("ps1_r_weff", lit.PS1_R_WIDTH_EFF_A, lit.PS1_R_REF)]
        con.executemany("INSERT INTO sn_literature VALUES (?,?,?)", lits)
        meta(con, peak_at=utcnow())
    for r in rows:
        print(f"{r[0]}: t_peak {r[1]:+.2f} d [{r[2]:.2f}, {r[3]:.2f}] "
              f"(n={r[4]}, brightest {r[6]:.3f} at +{r[5]:.2f} d)")
    print("first clean RLMT epoch:", {k: round(v, 2) for k, v in first.items()})
    return 0


COMMANDS["peak"] = cmd_peak


# ===========================================================================
# excessgate  (SN-S6-0-excess-gate)
# ===========================================================================
#: The narrowband systematic the gate judges against, defined before the
#: SN's own (H - '1') colour was formed:
#:   sigma_rep   robust night-to-night rms of the (H - '1') nightly-mean
#:               colour of constant calibration stars with r < NB_REP_RMAX
#:               (the bright end, nearest the SN's S/N);
#:   sigma_loc   the locus-transfer error: the estimator subtracts the
#:               field-star (H - '1') locus AT THE SN's colour, which lies
#:               blueward of every calibration star, so the locus must be
#:               extrapolated; sigma_loc = |linear - quadratic| locus at
#:               the SN colour.  Stellar H-alpha ABSORPTION (EW 2-6 A in
#:               F-K stars, i.e. 3-10% of a 65 A band) is what bends the
#:               locus, and it is absent from the SN's continuum.
#:   sigma_NB = sqrt(sigma_rep^2 + sigma_loc^2).
NB_REP_RMAX = 14.0
GATE_FACTOR = 3.0


def ab_flambda(mag_ab, lam_a):
    """f_lambda (erg s^-1 cm^-2 A^-1) of an AB magnitude at lam_a."""
    fnu = 3631e-23 * 10 ** (-0.4 * mag_ab)
    return fnu * 2.99792458e18 / lam_a ** 2


def narrowband_star_colours(con, cat):
    """Nightly-mean (H - '1') natural colour of each calibration star."""
    cm = sp.calib_star_mask(cat)
    out = {}
    for code in ("H", "1"):
        d = {}
        for n, sid, fl, er, zp, pk, cap in con.execute("""
                SELECT f.night, p.star_id, p.flux, p.err, c.zp_nat, p.peak,
                       f.lin_cap_adu
                FROM sn_star_phot p JOIN sn_frames f USING (obs_rowid)
                JOIN sn_cal_frames c USING (obs_rowid)
                WHERE f.epoch_role = 'campaign' AND f.filter = ? AND c.used = 1
                  AND p.flux > 0 AND p.flag = 0""", (code,)):
            if not cm[sid] or fl / er < 20 or (pk and pk >= cap):
                continue
            d.setdefault((n, sid), []).append(-2.5 * np.log10(fl) - zp)
        out[code] = {k: float(np.mean(v)) for k, v in d.items()}
    by = {}
    for k in set(out["H"]) & set(out["1"]):
        by.setdefault(k[1], []).append(out["H"][k] - out["1"][k])
    return by


def cmd_excessgate(args) -> int:
    from macro_sn import literature as lit
    con = connect()
    cat = load_refcat()
    col = cat["gmag"] - cat["imag"]
    by = narrowband_star_colours(con, cat)
    reps, X, Y = [], [], []
    for sid, v in by.items():
        if len(v) < 4:
            continue
        a = np.array(v)
        X.append(col[sid]); Y.append(a.mean())
        if cat["rmag"][sid] < NB_REP_RMAX:
            reps.extend(a - a.mean())
    reps, X, Y = map(np.array, (reps, X, Y))
    sig_rep = float(1.4826 * np.median(np.abs(reps)))
    c1 = np.polyfit(X, Y, 1)
    c2 = np.polyfit(X, Y, 2)
    # widths from zero-point ratios: K = zp_nat + 2.5 log t per frame;
    # W_X / W_R = 10^(-0.4 (K_X - K_R)), nightly medians, then the median
    # over nights both codes share.
    K = {}
    for code, night, k in con.execute("""
            SELECT c.code, f.night, c.zp_nat + 2.5 * log10(f.exptime) / 1.0
            FROM sn_cal_frames c JOIN sn_frames f USING (obs_rowid)
            WHERE c.role = 'campaign' AND c.used = 1
              AND c.code IN ('R', 'H', '1')"""):
        K.setdefault((code, night), []).append(k)
    widths = {}
    for code in ("H", "1"):
        dk = [np.median(K[(code, n)]) - np.median(K[("R", n)])
              for (c, n) in K if c == code and ("R", n) in K]
        dk = np.array(dk)
        w = lit.PS1_R_WIDTH_EFF_A * 10 ** (-0.4 * dk)
        widths[code] = (float(np.median(w)),
                        float(1.4826 * np.median(np.abs(w - np.median(w)))),
                        len(dk))
    # the SN's continuum and colour at the three flash epochs
    nights = {}
    for night, code, ph, m in con.execute(
            "SELECT night, code, phase_d, m_nat FROM sn_nightly "
            "WHERE code IN ('1', 'H', 'G', 'I') AND phase_d < 7"):
        nights.setdefault(night, {})[code] = (ph, m)
    sets = {r[0]: r[1] for r in con.execute(
        "SELECT code, cterm FROM sn_cal_sets WHERE role='campaign'")}
    gate_rows = []
    any_pass = False
    for tgt in lit.FLASH_EPOCHS_D:
        night = min(nights, key=lambda n: abs(
            nights[n].get("1", (99,))[0] - tgt) if "1" in nights[n] else 99)
        ph, m1 = nights[night]["1"]
        # SN colour: from that night's G and I if both exist, else the
        # nearest night's (the colour enters only via the locus term)
        cn = sorted((n for n in nights if "G" in nights[n] and "I" in nights[n]),
                    key=lambda n: abs(nights[n]["G"][0] - ph))
        if cn:
            gg, ii = nights[cn[0]]["G"][1], nights[cn[0]]["I"][1]
            csn = float(sp.colour_from_natural(gg, ii, sets["G"], sets["I"]))
        else:
            csn = float("nan")
        if not np.isfinite(csn):
            got = con.execute("SELECT colour_gi FROM sn_nightly WHERE code='G' "
                              "ORDER BY phase_d LIMIT 1").fetchone()
            csn = float(got[0])
        sig_loc = float(abs(np.polyval(c1, csn) - np.polyval(c2, csn)))
        sig_nb = float(np.hypot(sig_rep, sig_loc))
        flam = ab_flambda(m1, 6724.0)
        f_up = lit.F_HALPHA_MAX_CGS
        f_mid = lit.F_HALPHA_MAX_CGS * lit.F_HALPHA_FADE_FACTOR ** (
            -max(ph - lit.F_HALPHA_MAX_PHASE_D, 0)
            / (lit.F_HALPHA_FADE_SPAN_D - lit.F_HALPHA_MAX_PHASE_D))
        for W in lit.PLAUSIBLE_WIDTHS_A + (widths["H"][0],):
            eps_up = f_up / (flam * W)
            eps_mid = f_mid / (flam * W)
            dm_up = 2.5 * np.log10(1 + eps_up)
            dm_mid = 2.5 * np.log10(1 + eps_mid)
            passes = dm_up >= GATE_FACTOR * sig_nb
            any_pass |= passes
            gate_rows.append((tgt, night, ph, m1, flam, csn, W, f_up, f_mid,
                              eps_up, eps_mid, dm_up, dm_mid, sig_rep, sig_loc,
                              sig_nb, int(passes),
                              int(dm_up >= GATE_FACTOR * sig_rep)))
    verdict = "PROCEED" if any_pass else "DEMOTE S6a"
    con.executescript("""
        DROP TABLE IF EXISTS sn_excess_gate;
        CREATE TABLE sn_excess_gate (epoch_d REAL, night TEXT, phase_d REAL,
            m_cont REAL, flam_cont REAL, colour_gi REAL, width_a REAL,
            f_line_upper REAL, f_line_central REAL, eps_upper REAL,
            eps_central REAL, dm_upper REAL, dm_central REAL, sig_rep REAL,
            sig_loc REAL, sig_nb REAL, passes INTEGER, passes_rep_only INTEGER);
        DROP TABLE IF EXISTS sn_nb_widths;
        CREATE TABLE sn_nb_widths (code TEXT PRIMARY KEY, width_a REAL,
            width_mad REAL, n_nights INTEGER);""")
    with con:
        con.executemany("INSERT INTO sn_excess_gate VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", gate_rows)
        con.executemany("INSERT INTO sn_nb_widths VALUES (?,?,?,?)",
                        [(k, *v) for k, v in widths.items()])
        meta(con, excess_gate_at=utcnow(), excess_gate_verdict=verdict,
             nb_locus_lin=json.dumps(list(c1)), nb_locus_quad=json.dumps(list(c2)),
             nb_locus_n=len(X), nb_rep_n=len(reps),
             nb_locus_colour_range=json.dumps([float(X.min()), float(X.max())]))
    print(f"widths from ZP ratios: " + ", ".join(
        f"{k} {v[0]:.0f}+-{v[1]:.0f} A ({v[2]} nights)" for k, v in widths.items()))
    print(f"sigma_rep = {sig_rep:.4f} mag ({len(reps)} star-nights, r<{NB_REP_RMAX})")
    for g in gate_rows:
        print(f"+{g[0]:.1f} d ({g[1]}): W={g[6]:.0f} eps_up={g[9]:.3f} "
              f"dm_up={g[11]:.3f} mid={g[12]:.3f} | 3 sig_NB={3*g[15]:.3f} "
              f"(rep {g[13]:.3f}, loc {g[14]:.3f} at g-i={g[5]:+.2f}) "
              f"pass={g[16]} rep-only pass={g[17]}")
    print("GATE:", verdict)
    return 0


COMMANDS["excessgate"] = cmd_excessgate

# ===========================================================================
# residuals  (SN-S7b-residuals-table)
# ===========================================================================
#: Interpolation of the published curve to an RLMT epoch: weighted straight
#: line through the published points within +/- RESID_WINDOW_D, required to
#: have points on BOTH sides of the epoch (no extrapolation).
RESID_WINDOW_D = 3.0


def local_line(t, m, e, t0, win=RESID_WINDOW_D):
    w = np.abs(t - t0) <= win
    if w.sum() < 2 or not (np.any(t[w] < t0) and np.any(t[w] > t0)):
        return None
    a, b, cov = sp.wls_line(t[w] - t0, m[w], 1 / e[w] ** 2)
    # scatter-aware error: formal error inflated by the fit's own chi2_nu
    # when the published points disagree among themselves (never deflated)
    res = m[w] - a - b * (t[w] - t0)
    dof = w.sum() - 2
    c2 = float(np.sum(res ** 2 / e[w] ** 2) / dof) if dof > 0 else 1.0
    return a, float(np.sqrt(cov[0, 0] * max(c2, 1.0))), int(w.sum())


def cmd_residuals(args) -> int:
    from macro_sn import literature as lit
    con = connect()
    li = load_li25()
    L = {b: np.array([(x[1], x[3], x[4]) for x in li if x[2] == b])
         for b in ("g", "r", "i")}
    ours = con.execute("""SELECT night, code, phase_d, m_ps1, err, colour_gi,
                                 colour_in_range FROM sn_nightly
                          WHERE code IN ('G','R','I') ORDER BY phase_d""").fetchall()
    rows = []
    for night, code, ph, m, e, c, inr in ours:
        b = PS1_OF[code]
        got = local_line(L[b][:, 0], L[b][:, 1], L[b][:, 2], ph)
        gg = local_line(L["g"][:, 0], L["g"][:, 1], L["g"][:, 2], ph)
        rr = local_line(L["r"][:, 0], L["r"][:, 1], L["r"][:, 2], ph)
        if got is None or gg is None or rr is None:
            continue
        x = gg[0] - rr[0]
        B0, B1, _ = lit.TONRY12_SDSS_TO_PS1[b]
        li_ps1 = got[0] + B0 + B1 * x
        err = float(np.hypot(e, got[1]))
        rows.append((night, code, ph, m, e, got[0], li_ps1, got[1], got[2],
                     m - li_ps1, err, c, inr, x))
    con.executescript("""
        DROP TABLE IF EXISTS sn_resid;
        CREATE TABLE sn_resid (night TEXT, code TEXT, phase_d REAL,
            rlmt REAL, rlmt_err REAL, li_sdss REAL, li_ps1 REAL, li_err REAL,
            li_n INTEGER, resid REAL, err REAL, colour_gi REAL,
            colour_in_range INTEGER, li_gr_sdss REAL,
            PRIMARY KEY (night, code));
        DROP TABLE IF EXISTS sn_resid_summary;
        CREATE TABLE sn_resid_summary (code TEXT, subset TEXT, n INTEGER,
            offset REAL, offset_err REAL, rms REAL, chi2nu REAL, dof INTEGER,
            colour_lo REAL, colour_hi REAL, PRIMARY KEY (code, subset));
        DROP TABLE IF EXISTS sn_resid_ztf;
        CREATE TABLE sn_resid_ztf (mjd REAL, phase_d REAL, ztf_g REAL,
            ztf_err REAL, rlmt_g REAL, rlmt_err REAL, night TEXT, resid REAL);""")
    summ = []
    sets = {r[0]: (r[1], r[2]) for r in con.execute(
        "SELECT code, colour_lo, colour_hi FROM sn_cal_sets WHERE role='campaign'")}
    for code in BROAD:
        for subset, cond in (("all", lambda r: True),
                             ("colour inside training range", lambda r: r[12] == 1),
                             ("colour outside training range", lambda r: r[12] == 0)):
            rr = [r for r in rows if r[1] == code and cond(r)]
            if len(rr) < 2:
                continue
            v = np.array([r[9] for r in rr]); s_ = np.array([r[10] for r in rr])
            mu, err, c2, n = sp.weighted_mean(v, s_)
            summ.append((code, subset, n, mu, err, float(np.std(v)), c2, n - 1,
                         *sets[code]))
        # the same comparison WITHOUT the Tonry SDSS->PS1 term, so the size
        # of that (stellar-SED) transformation is visible beside the result
        rr = [r for r in rows if r[1] == code]
        if len(rr) >= 2:
            v = np.array([r[3] - r[5] for r in rr])
            s_ = np.array([r[10] for r in rr])
            mu, err, c2, n = sp.weighted_mean(v, s_)
            summ.append((code, "all, published SDSS mags untransformed", n, mu,
                         err, float(np.std(v)), c2, n - 1, *sets[code]))
    # ZTF alert photometry (ALeRCE, difference-image PSF; ZTF g is not PS1 g:
    # a cross-check only, no transformation applied)
    zt = json.loads((EXT / "ztf_alerce_lc.json").read_text())
    zrows = []
    G = {r[0]: (r[1], r[2], r[3]) for r in con.execute(
        "SELECT night, phase_d, m_ps1, err FROM sn_nightly WHERE code='G'")}
    for d in zt.get("detections", []):
        if d.get("fid") != 1 or d.get("magpsf") is None:
            continue
        ph = d["mjd"] - lit.LI25_T0_MJD
        best = min(G.items(), key=lambda kv: abs(kv[1][0] - ph))
        if abs(best[1][0] - ph) <= 1.0:
            zrows.append((d["mjd"], ph, d["magpsf"], d["sigmapsf"], best[1][1],
                          best[1][2], best[0], best[1][1] - d["magpsf"]))
    with con:
        con.executemany("INSERT INTO sn_resid VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        rows)
        con.executemany("INSERT INTO sn_resid_summary VALUES (?,?,?,?,?,?,?,?,?,?)",
                        summ)
        con.executemany("INSERT INTO sn_resid_ztf VALUES (?,?,?,?,?,?,?,?)", zrows)
        meta(con, residuals_at=utcnow())
    for r in summ:
        print(f"{r[0]} [{r[1]}]: N={r[2]} offset={r[3]:+.3f}+-{r[4]:.3f} "
              f"rms={r[5]:.3f} chi2nu={r[6]:.2f}/{r[7]}")
    for z in zrows:
        print(f"ZTF g +{z[1]:.2f} d: ZTF {z[2]:.3f}+-{z[3]:.3f} RLMT {z[4]:.3f} "
              f"-> {z[7]:+.3f}")
    return 0


COMMANDS["residuals"] = cmd_residuals

# ===========================================================================
# variability  (SN-S8-variability-limits)
# ===========================================================================
N_BOOT = 5000
N_INJ = 200
INJ_AMPS = np.geomspace(0.004, 0.40, 18)
VAR_SEED = 2023_0519


def _nightly_from_frames(t, m, s, night):
    """Weighted nightly means (phase, mag, err) from frame arrays."""
    out = []
    for n in np.unique(night):
        k = night == n
        w = 1 / s[k] ** 2
        out.append((np.sum(w * t[k]) / np.sum(w), np.sum(w * m[k]) / np.sum(w),
                    1 / np.sqrt(np.sum(w))))
    a = np.array(out)
    return a[:, 0], a[:, 1], a[:, 2]


def cmd_variability(args) -> int:
    """Night-to-night variability and bump limits (no intra-night
    periodogram — ruled out).  For each band: BIC-chosen spline trend;
    joint trend+bump and trend+sinusoid searches; whole-night bootstrap
    null (standardised nightly residuals resampled onto the real epochs
    and added to the fitted trend); injection-recovery into the per-frame
    data, re-averaged to nightly means, searched exactly as the data were."""
    rng = np.random.default_rng(VAR_SEED)
    con = connect()
    lim_rows, obs_rows = [], []
    for code in BROAD:
        fr = np.array(con.execute(
            "SELECT phase_d, m_nat, sig_tot, night FROM sn_phot "
            "WHERE usable = 1 AND code = ? ORDER BY phase_d", (code,)).fetchall(),
            dtype=object)
        t = fr[:, 0].astype(float); m = fr[:, 1].astype(float)
        s_ = fr[:, 2].astype(float); night = fr[:, 3].astype(str)
        tn, mn, en = _nightly_from_frames(t, m, s_, night)
        ks, bics = sp.choose_knot_spacing(tn, mn, en)
        J = sp.JointSearch(tn, en, ks)
        tr = J.trend(mn)
        r = mn - tr
        z = r / en
        bsnr, bt0, bsig, bamp = J.bump(mn)
        cchi, cP, camp = J.sine_stat(mn)
        nb, nc = np.empty(N_BOOT), np.empty(N_BOOT)
        for i in range(N_BOOT):
            yb = tr + z[rng.integers(0, len(z), len(z))] * en
            nb[i] = J.bump(yb)[0]
            nc[i] = J.sine_stat(yb)[0]
        thr_b = float(np.quantile(nb, 1 - sp.FAP_LEVEL))
        thr_c = float(np.quantile(nc, 1 - sp.FAP_LEVEL))
        fap_b = float((np.sum(nb >= bsnr) + 1) / (N_BOOT + 1))
        fap_c = float((np.sum(nc >= cchi) + 1) / (N_BOOT + 1))
        chi2_trend = float(np.sum(z ** 2) / max(len(r) - (J.H.trace()), 1))
        obs_rows.append((code, len(tn), len(t), chi2_trend, bsnr, bt0, bsig,
                         bamp, thr_b, fap_b, cchi, cP, camp, thr_c, fap_c,
                         float(np.std(r)), ks))
        print(f"{code}: {len(tn)} nights; knots every {ks:g} d; chi2nu about "
              f"trend {chi2_trend:.2f}; bump SNR {bsnr:.2f} (thr {thr_b:.2f}, "
              f"FAP {fap_b:.4f}); sine dchi2 {cchi:.1f} (thr {thr_c:.1f}, "
              f"FAP {fap_c:.4f})")
        shapes = [("bump", w) for w in sp.BUMP_SIGMAS_D] + \
                 [("sine", P) for P in sp.SINE_PERIODS_D]
        for shape, par in shapes:
            rec, bias, rec_cov = [], [], []
            for A in INJ_AMPS:
                hits, bs, hc, nc_ = 0, [], 0, 0
                for _ in range(N_INJ):
                    if shape == "bump":
                        t0 = rng.uniform(tn.min() + 1, tn.max() - 1)
                        inj = -A * np.exp(-0.5 * ((t - t0) / par) ** 2)
                    else:
                        ph0 = rng.uniform(0, 2 * np.pi)
                        inj = A * np.sin(2 * np.pi * t / par + ph0)
                    _t2, mn2, _e2 = _nightly_from_frames(t, m + inj, s_, night)
                    if shape == "bump":
                        det = J.bump(mn2)[0] >= thr_b
                        # signed matched-cell bias: amplitude fitted AT the
                        # injected centre/width, minus the injection
                        afit = -(J.bump_at(mn2, t0, par) - J.bump_at(mn, t0, par))
                    else:
                        det = J.sine_stat(mn2)[0] >= thr_c
                        afit = J.sine_stat(mn2, only=par)[2]
                    hits += int(det)
                    bs.append(afit - A)
                    # a bump is COVERED when a night falls within one sigma
                    # of its centre — the cadence, not the noise, decides
                    # the rest (their fraction is the completeness ceiling)
                    cov = shape != "bump" or np.any(np.abs(tn - t0) <= par)
                    nc_ += int(cov)
                    hc += int(det and cov)
                rec.append(hits / N_INJ)
                rec_cov.append(hc / max(nc_, 1))
                bias.append(float(np.median(bs)))
            rec = np.array(rec)
            rec_cov = np.array(rec_cov)

            def a_at(curve, level):
                ok = np.nonzero(curve >= level)[0]
                if not len(ok):
                    return float("nan"), -1
                k = ok[0]
                return (float(INJ_AMPS[k] if k == 0 else np.interp(
                    level, curve[k - 1:k + 1], INJ_AMPS[k - 1:k + 1])), k)
            a90, k90 = a_at(rec_cov, sp.RECOVERY_LEVEL)
            a50, _ = a_at(rec, 0.5)
            plateau = float(rec[-1])
            # predicted scale (standing rule 2): what pure noise allows with
            # NO trend — (threshold + 1.28) x the matched-filter amplitude
            # error, averaged over centre (bump) / the sine equivalent.
            if shape == "bump":
                den = [np.sum(np.exp(-((tn - c0) / par) ** 2) / en ** 2)
                       for c0 in np.linspace(tn.min() + 1, tn.max() - 1, 50)]
                apred = float(np.mean((thr_b + 1.2816) / np.sqrt(den)))
            else:
                apred = float((np.sqrt(thr_c) + 1.2816)
                              * np.sqrt(2.0 / np.sum(1 / en ** 2)))
            lim_rows.append((code, shape, par, a90, apred,
                             bias[k90] if k90 >= 0 else None,
                             json.dumps([round(x, 3) for x in rec.tolist()]),
                             a50, plateau,
                             json.dumps([round(x, 3) for x in rec_cov.tolist()])))
            print(f"   {shape} {par:>4} d: A90(covered) {a90*1e3:.0f} mmag, "
                  f"A50(all) {a50*1e3:.0f}, ceiling {plateau:.2f}; "
                  f"90% recovery at {a90*1e3:.0f} mmag "
                  f"(noise-only scale {apred*1e3:.0f}; signed bias "
                  f"{(bias[k90] if k90 >= 0 else float('nan'))*1e3:+.1f} mmag)")
    con.executescript("""
        DROP TABLE IF EXISTS sn_var_observed;
        CREATE TABLE sn_var_observed (code TEXT PRIMARY KEY, n_nights INTEGER,
            n_frames INTEGER, chi2nu_trend REAL, bump_snr REAL, bump_t0 REAL,
            bump_sigma REAL, bump_amp REAL, bump_thr REAL, bump_fap REAL,
            sine_chi2 REAL, sine_period REAL, sine_amp REAL, sine_thr REAL,
            sine_fap REAL, resid_rms REAL, knot_spacing_d REAL);
        DROP TABLE IF EXISTS sn_var_limits;
        CREATE TABLE sn_var_limits (code TEXT, shape TEXT, scale_d REAL,
            a90_mag REAL, a_pred_mag REAL, signed_bias_mag REAL,
            recovery_curve TEXT, a50_all_mag REAL, completeness REAL,
            recovery_curve_covered TEXT,
            PRIMARY KEY (code, shape, scale_d));""")
    with con:
        con.executemany("INSERT INTO sn_var_observed VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", obs_rows)
        con.executemany("INSERT INTO sn_var_limits VALUES (?,?,?,?,?,?,?,?,?,?)",
                        lim_rows)
        meta(con, variability_at=utcnow(), inj_amps=json.dumps(INJ_AMPS.tolist()),
             n_boot=N_BOOT, n_inj=N_INJ)
    return 0


COMMANDS["variability"] = cmd_variability

# ===========================================================================
# latetime  (SN-S9-late-time): template-subtracted forced photometry of the
# post-fade stacks — the faint regime — and injection-defined 5-sigma limits
# ===========================================================================
#: The only SN-free reference this archive holds is the pre-explosion
#: AC4040 epoch (2023-05-04); it has G and R but no I, and the template
#: table shows both out of focus (R 4.3"; G beyond the fit range).  The
#: late epochs are subtracted against it after PSF matching to the broader
#: image; i has no SN-free reference and gets no limit.
LATE_EPOCHS = (
    ("+366 d", "iKon", "g", ("2024-05-18",), "G"),
    ("+366 d", "iKon", "r", ("2024-05-18",), "R"),
    ("+1035 d", "QHY600", "g", ("2026-03-21",), "G"),
    ("+1035 d", "QHY600", "r", ("2026-03-21", "2026-03-22"), "R"),
    ("+1045 d", "QHY600", "g", ("2026-04-01",), "G"),
    ("+1045 d", "QHY600", "r", ("2026-04-01",), "R"),
)
#: Random forced apertures that define the empirical noise: this many, in
#: an annulus of 30-90" about the SN (on the arm, same host background).
N_NOISE_AP = 300
NOISE_ANN_ARCSEC = (30.0, 90.0)
INJ_SIGMAS = (3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0)
N_LATE_INJ = 200


def grid_stack(con, cat, camera, code, nights, role):
    from macro_sn import snio
    rows = con.execute(f"""SELECT path, wcs_header, filter, camera FROM sn_frames
        WHERE status='measured' AND camera=? AND filter=? AND epoch_role=?
        AND night IN ({",".join("?" * len(nights))}) ORDER BY obs_rowid""",
                       (camera, code, role, *nights)).fetchall()
    ims = []
    for r in rows:
        g, grid = _grid_image(r[0], r[1], r[2], r[3])
        ims.append(g)
    if not ims:
        return None
    n = 2 * SUB_HALF + 1
    xs, ys = _grid_stars(grid, cat, n)
    f0 = snio.grid_fwhm(ims[0], xs, ys) or 5.0
    scaled = []
    for im in ims:
        k, _n = _ratio(ims[0], im, xs, ys, 2.0 * f0)
        if k is not None:
            scaled.append(im * k)
    st = np.nanmedian(np.array(scaled), axis=0)
    return st, snio.grid_fwhm(st, xs, ys), len(scaled), (xs, ys), grid


def cmd_latetime(args) -> int:
    global _CAL
    from macro_sn import snio
    import sep
    _CAL = snio.load_calib()
    rng = np.random.default_rng(366)
    con = connect()
    cat = load_refcat()
    xw = {(r[0], r[1]): r[2] for r in con.execute(
        "SELECT camera, code, cterm FROM sn_crosswalk WHERE identified = 1 "
        "AND role = 'template_post'")}
    tpl_cache = {}
    out = []
    n = 2 * SUB_HALF + 1
    c = float(SUB_HALF)
    yy, xx = np.mgrid[0:n, 0:n]
    for label, cam, code, nights, tcode in LATE_EPOCHS:
        sci = grid_stack(con, cat, cam, code, nights, "template_post")
        if tcode not in tpl_cache:
            tpl_cache[tcode] = grid_stack(con, cat, "AC4040", tcode,
                                          ("2023-05-04",), "template_pre")
        tpl = tpl_cache[tcode]
        if sci is None or tpl is None:
            out.append((label, cam, code, ",".join(nights), 0, None, None, None,
                        None, None, None, None, None, None, None, None, None,
                        "no frames"))
            continue
        S, fs, ns, stars, grid = sci
        T, ft, nt, _st, _g = tpl
        if ft is None:
            # the pre-explosion G stars are donuts beyond the Gaussian fit
            # range: no defensible PSF match exists
            out.append((label, cam, code, ",".join(nights), ns, fs, None, None,
                        None, None, None, None, None, None, None, None, None,
                        "template PSF unmeasurable (defocused pre-explosion "
                        f"{tcode}); no subtraction"))
            print(f"{label} {cam} {code}: template {tcode} PSF unmeasurable — skipped")
            continue
        got = subtract_and_measure(S, fs, T, ft, stars)
        if got is None:
            continue
        _fa, f_sn, fm, nst, k, D = got
        r_ap = 1.5 * fm
        # empirical noise from random apertures on the arm
        scale = snio.GRID_SCALE_ARCSEC
        rr = rng.uniform(NOISE_ANN_ARCSEC[0] / scale, NOISE_ANN_ARCSEC[1] / scale,
                         N_NOISE_AP)
        th = rng.uniform(0, 2 * np.pi, N_NOISE_AP)
        px, py = c + rr * np.cos(th), c + rr * np.sin(th)
        fr_ = np.array([snio.aper_sum(D, x, y, r_ap) for x, y in zip(px, py)])
        fr_ = fr_[np.isfinite(fr_)]
        sig = float(1.4826 * np.median(np.abs(fr_ - np.median(fr_))))
        # zero point of the PSF-matched science stack from grid stars
        a, _b, _fm = snio.gauss_match(S, fs, T, ft)
        xs, ys = stars
        fl, _e, _f = sep.sum_circle(np.nan_to_num(a).astype(np.float64), xs, ys,
                                    r_ap, bkgann=(4 * fm, 6 * fm))
        ra, dec = grid.all_pix2world(xs, ys, 0)
        idx = [int(np.argmin((cat["RA_ICRS"] - ra_) ** 2 + (cat["DE_ICRS"] - de_) ** 2))
               for ra_, de_ in zip(ra, dec)]
        mcat = cat[f"{code}mag"][idx]
        colr = (cat["gmag"] - cat["imag"])[idx]
        ok = (fl > 0) & np.isfinite(mcat)
        zps = (mcat[ok] + xw.get((cam, code), 0.0) * (colr[ok] - sp.COLOUR_PIVOT)
               + 2.5 * np.log10(fl[ok]))
        zp = float(np.median(zps))
        zp_err = float(1.4826 * np.median(np.abs(zps - zp)) / np.sqrt(len(zps)))
        lim5 = zp - 2.5 * np.log10(5 * sig)
        snr = f_sn / sig
        m_sn = zp - 2.5 * np.log10(f_sn) if f_sn > 0 else None
        # injection: Gaussian sources of the matched PSF injected into the
        # difference image at random arm positions, recovered by the same
        # forced aperture; detection = SNR >= 5
        recs = []
        for kk in INJ_SIGMAS:
            hits = 0
            for _ in range(N_LATE_INJ):
                r0 = rng.uniform(NOISE_ANN_ARCSEC[0] / scale,
                                 NOISE_ANN_ARCSEC[1] / scale)
                t0 = rng.uniform(0, 2 * np.pi)
                x0, y0 = c + r0 * np.cos(t0), c + r0 * np.sin(t0)
                R = int(4 * fm)
                xi, yi = int(round(x0)), int(round(y0))
                Dc = D[yi - R:yi + R + 1, xi - R:xi + R + 1].copy()
                g = snio.gaussian_psf(2 * R + 1, R + (x0 - xi), R + (y0 - yi),
                                      fm, kk * sig)
                f_rec = snio.aper_sum(np.pad(Dc + g, 2), R + 2 + (x0 - xi),
                                      R + 2 + (y0 - yi), r_ap) \
                    if Dc.shape == g.shape else np.nan
                hits += int(np.isfinite(f_rec) and f_rec / sig >= 5.0)
            recs.append(hits / N_LATE_INJ)
        recs = np.array(recs)
        okr = np.nonzero(recs >= 0.9)[0]
        k90 = (float(np.interp(0.9, recs[okr[0] - 1:okr[0] + 1],
                               INJ_SIGMAS[okr[0] - 1:okr[0] + 1]))
               if len(okr) and okr[0] > 0 else (INJ_SIGMAS[0] if len(okr) else np.nan))
        lim90 = zp - 2.5 * np.log10(k90 * sig) if np.isfinite(k90) else None
        out.append((label, cam, code, ",".join(nights), ns, fs, ft, fm, zp,
                    zp_err, int(ok.sum()), sig, f_sn, snr, m_sn, lim5, lim90,
                    json.dumps(dict(zip([str(x) for x in INJ_SIGMAS],
                                        recs.round(3).tolist())))))
        print(f"{label} {cam} {code}: {ns} frames vs {tcode} pre-explosion; "
              f"matched FWHM {fm*scale:.2f}\"; ZP {zp:.2f}+-{zp_err:.2f} "
              f"({int(ok.sum())} stars); SN forced SNR {snr:+.1f}"
              f"{'' if m_sn is None else f' (m={m_sn:.2f})'}; 5-sigma limit "
              f"{lim5:.2f}; 90%-recovery {lim90 if lim90 is None else round(lim90,2)}")
    con.executescript("""
        DROP TABLE IF EXISTS sn_latetime;
        CREATE TABLE sn_latetime (epoch TEXT, camera TEXT, code TEXT,
            nights TEXT, n_frames INTEGER, fwhm_sci_px REAL, fwhm_tpl_px REAL,
            fwhm_matched_px REAL, zp REAL, zp_err REAL, n_zp_stars INTEGER,
            sigma_flux REAL, sn_flux REAL,
            sn_snr REAL, sn_mag REAL, lim5_mag REAL, lim90_mag REAL,
            recovery TEXT);""")
    with con:
        con.executemany("INSERT INTO sn_latetime VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", out)
        meta(con, latetime_at=utcnow())
    return 0


COMMANDS["latetime"] = cmd_latetime

# ===========================================================================
# release  (SN-S10-release): machine-readable tables + README, all from
# the two databases; nothing typed.
# ===========================================================================
REL = OUT / "release"


def _csv(path, header, rows):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(["" if v is None else (f"{v:.6f}" if isinstance(v, float)
                                              else v) for v in r])
    return len(rows)


def cmd_release(args) -> int:
    REL.mkdir(parents=True, exist_ok=True)
    con = connect(read_only=True)
    m = manifest()
    files = {}
    q = ("SELECT p.obs_rowid, f.path, p.night, p.bjd_tdb, p.phase_d, p.code, "
         "p.exptime, p.airmass, p.census_class, p.sn_peak_adu, p.lin_cap_adu, "
         "p.usable, p.exclusion, p.m_nat, p.sig_phot, p.sig_scint, "
         "p.sig_floor, p.sig_zp, p.sig_tot, p.colour_gi, p.colour_in_range, "
         "p.m_ps1 FROM sn_phot p JOIN sn_frames f USING (obs_rowid) "
         "ORDER BY p.bjd_tdb")
    hdr = ["obs_rowid", "archive_path", "night", "bjd_tdb", "phase_d",
           "filter_code", "exptime_s", "airmass", "census_class",
           "sn_peak_adu", "linearity_cap_adu", "usable", "exclusion",
           "mag_natural", "err_photon", "err_scint", "err_floor", "err_zp",
           "err_total", "sn_colour_gi_ps1", "colour_in_training_range",
           "mag_ps1"]
    files["sn2023ixf_rlmt_frames.csv"] = _csv(REL / "sn2023ixf_rlmt_frames.csv",
                                              hdr, con.execute(q).fetchall())
    files["sn2023ixf_rlmt_nightly.csv"] = _csv(
        REL / "sn2023ixf_rlmt_nightly.csv",
        ["night", "filter_code", "n_frames", "bjd_tdb", "phase_d", "mag_natural",
         "err", "chi2nu_within_night", "mag_ps1", "sn_colour_gi_ps1",
         "colour_in_training_range"],
        con.execute("SELECT night, code, n, bjd_tdb, phase_d, m_nat, err, chi2nu,"
                    " m_ps1, colour_gi, colour_in_range FROM sn_nightly "
                    "ORDER BY code, phase_d").fetchall())
    files["saturation_matrix.csv"] = _csv(
        REL / "saturation_matrix.csv",
        ["night", "filter_code", "phase_d", "band_role", "n_frames", "n_measured",
         "n_clean", "n_suspect", "n_rejected", "n_bounded_clean",
         "n_undetermined", "n_spectra", "max_peak_adu", "min_peak_adu",
         "min_exptime", "max_exptime"],
        m.execute("SELECT * FROM sn_g0_matrix ORDER BY night, filter").fetchall())
    files["slot6_frames_not_promoted.csv"] = _csv(
        REL / "slot6_frames_not_promoted.csv",
        ["obs_rowid", "archive_path", "night", "phase_d", "exptime_s",
         "s2c_verdict", "s2c_class"],
        m.execute("SELECT obs_rowid, path, night, phase_d, exptime, "
                  "dispersion_verdict, dispersion_class FROM sn_g0_frames "
                  "WHERE filter = '6' ORDER BY jd").fetchall())
    files["calibration_sets.csv"] = _csv(
        REL / "calibration_sets.csv",
        [d[1] for d in con.execute("PRAGMA table_info(sn_cal_sets)")],
        con.execute("SELECT * FROM sn_cal_sets").fetchall())
    files["filter_crosswalk.csv"] = _csv(
        REL / "filter_crosswalk.csv",
        [d[1] for d in con.execute("PRAGMA table_info(sn_crosswalk)")],
        con.execute("SELECT * FROM sn_crosswalk").fetchall())
    files["template_table.csv"] = _csv(
        REL / "template_table.csv",
        [d[1] for d in con.execute("PRAGMA table_info(sn_template_table)")],
        con.execute("SELECT * FROM sn_template_table").fetchall())
    files["residuals_vs_li2025.csv"] = _csv(
        REL / "residuals_vs_li2025.csv",
        [d[1] for d in con.execute("PRAGMA table_info(sn_resid)")],
        con.execute("SELECT * FROM sn_resid ORDER BY code, phase_d").fetchall())
    files["variability_limits.csv"] = _csv(
        REL / "variability_limits.csv",
        [d[1] for d in con.execute("PRAGMA table_info(sn_var_limits)")],
        con.execute("SELECT * FROM sn_var_limits").fetchall())
    files["late_time.csv"] = _csv(
        REL / "late_time.csv",
        [d[1] for d in con.execute("PRAGMA table_info(sn_latetime)")],
        con.execute("SELECT * FROM sn_latetime").fetchall())
    # README, regenerated from the databases (strategy §4 Step 10)
    tot = m.execute("SELECT count(*), min(night), max(night) FROM sn_g0_frames "
                    "WHERE epoch_role='campaign'").fetchone()
    nn = m.execute("SELECT count(DISTINCT night) FROM sn_g0_frames "
                   "WHERE epoch_role='campaign'").fetchone()[0]
    use = {r[0]: (r[1], r[2]) for r in con.execute(
        "SELECT code, sum(usable), count(DISTINCT CASE WHEN usable=1 THEN night "
        "END) FROM sn_phot GROUP BY code")}
    ver = con.execute("SELECT value FROM sn_build_meta WHERE key='code_version'"
                      ).fetchone()[0]
    lines = [
        "# SN 2023ixf — RLMT +5.4 to +50 d validation and limits release",
        "",
        f"Generated {utcnow()} by `pipeline/scripts/run_sn_photometry.py release` "
        f"({ver}, git {git_commit()[:10]}) from `products/sn/sn2023ixf.sqlite` "
        "and the manifest's Gate 0 tables. Every number below is a query.",
        "",
        f"Campaign: {tot[0]} unique frames (Gate 0 freeze, global "
        f"(basename, jd) dedup, alias-merged) on {nn} nights, night labels "
        f"{tot[1]} to {tot[2]} (local-noon split).",
        "",
        "Usable photometry per filter code (frames / nights):",
    ]
    for k in ("G", "R", "I", "H", "O", "1"):
        if k in use:
            lines.append(f"- {k}: {use[k][0]} / {use[k][1]}")
    lines += ["", "Files:", ""]
    for f, nrow in files.items():
        lines.append(f"- `{f}` — {nrow} rows")
    lines += ["",
              "Magnitudes: `mag_natural` is the RLMT natural system zeroed to "
              "PS1 (REFCAT2) at (g-i)=0.8; `mag_ps1` applies the campaign "
              "colour term with the SN's own (g-i). Errors are photon, "
              "scintillation (Young 1967 law, scale fitted), floor and zero "
              "point in quadrature. Frames flagged `usable=0` carry the "
              "rule that excluded them in `exclusion`.",
              "",
              "The slot-'6' series was triaged and NOT PROMOTED (Gate 0); its "
              "frames are listed as they are.",
              "",
              "Pipeline: this repository (`pipeline/macro_sn/`, "
              "`pipeline/scripts/run_sn_gate0.py`, "
              "`pipeline/scripts/run_sn_photometry.py`). Archive DOI: pending "
              "(Zenodo deposit is an owner action)."]
    (REL / "README.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


COMMANDS["release"] = cmd_release

# ===========================================================================
# psfcheck — the rule for S2c-INDETERMINATE frames
# ===========================================================================
#: S2c could not certify these frames direct or dispersed.  Gate 0 keeps
#: them (exclusion must be earned by a measurement), and this paper keeps
#: them for PHOTOMETRY only if its own measurement shows point sources:
#: at least PSF_MIN_STARS REFCAT2 stars detected within 2.5 px of their
#: catalogue positions under the frame's fitted WCS, with a median
#: elongation (a/b of the detection) no greater than PSF_MAX_ELONG.  A
#: slitless spectrum of a star is a streak of tens of pixels (elongation
#: >> 2); the threshold is fixed before looking at the indeterminate set,
#: and its distribution on the S2c-DIRECT frames is published beside it.
PSF_MIN_STARS = 8
PSF_MAX_ELONG = 1.3
PSF_DIRECT_SAMPLE = 60


def _psf_one(task):
    from astropy.io.fits import Header
    from astropy.wcs import WCS as _W
    from macro_sn import snio
    try:
        raw, _h = snio.load_data(task["path"])
        img, mask, _ = snio.calibrate(raw, task["filter"], "AC4040", _CAL)
        sub, bkg, obj = snio.detect(img, mask)
        w = _W(Header.fromstring(task["wcs_header"]))
        cat = _CAT
        x, y = w.all_world2pix(cat["RA_ICRS"], cat["DE_ICRS"], 0)
        ok = np.isfinite(x) & (cat["rmag"] < 16.5)
        ia, ib = snio.match_xy(x[ok], y[ok], obj["x"], obj["y"], 2.5)
        o = obj[ib]
        good = (o["peak"] < 0.9 * (task["lin_cap_adu"] or 1800)) & (o["flag"] == 0)
        el = o["a"][good] / np.maximum(o["b"][good], 1e-3)
        return (task["obs_rowid"], task["cls"], int(good.sum()),
                float(np.median(el)) if good.sum() else None)
    except Exception:
        return (task["obs_rowid"], task["cls"], 0, None)


def cmd_psfcheck(args) -> int:
    from multiprocessing import Pool
    con = connect()
    m = manifest()
    disp = {r[0]: r[1] for r in m.execute(
        "SELECT obs_rowid, dispersion_class FROM sn_g0_census")}
    m.close()
    rows = [dict(r) for r in con.execute(
        "SELECT obs_rowid, path, filter, wcs_header, lin_cap_adu FROM sn_frames "
        "WHERE status='measured' AND epoch_role='campaign' AND camera='AC4040' "
        "AND filter IN ('G','R','I','H','O','1')")]
    ind = [dict(r, cls="indeterminate") for r in rows
           if disp.get(r["obs_rowid"]) == "indeterminate"]
    dirr = [dict(r, cls="direct") for r in rows if disp.get(r["obs_rowid"]) == "direct"]
    rng = np.random.default_rng(6)
    dirr = [dirr[i] for i in rng.choice(len(dirr), min(PSF_DIRECT_SAMPLE, len(dirr)),
                                        replace=False)]
    with Pool(args.workers, initializer=_init_worker) as pool:
        res = pool.map(_psf_one, ind + dirr)
    con.executescript("""DROP TABLE IF EXISTS sn_psfcheck;
        CREATE TABLE sn_psfcheck (obs_rowid INTEGER PRIMARY KEY, s2c_class TEXT,
            n_stars INTEGER, median_elong REAL, passes INTEGER);""")
    out = [(o, c, n, e, int(n >= PSF_MIN_STARS and e is not None
                             and e <= PSF_MAX_ELONG)) for o, c, n, e in res]
    with con:
        con.executemany("INSERT INTO sn_psfcheck VALUES (?,?,?,?,?)", out)
        meta(con, psfcheck_at=utcnow(), psf_max_elong=PSF_MAX_ELONG,
             psf_min_stars=PSF_MIN_STARS)
    for cls in ("direct", "indeterminate"):
        e = np.array([x[3] for x in out if x[1] == cls and x[3] is not None])
        npass = sum(x[4] for x in out if x[1] == cls)
        ntot = sum(1 for x in out if x[1] == cls)
        print(f"{cls}: {npass}/{ntot} pass; median elongation "
              f"{np.median(e):.3f} (5-95%: {np.percentile(e, 5):.3f}-"
              f"{np.percentile(e, 95):.3f})")
    return 0


COMMANDS["psfcheck"] = cmd_psfcheck

# ===========================================================================
# paper + report
# ===========================================================================
def cmd_paper(args) -> int:
    from macro_sn import paper_sn
    for o in paper_sn.build(DB, MANIFEST):
        print("wrote", o.relative_to(REPO))
    return 0


def cmd_report(args) -> int:
    from macro_sn import report_sn
    print("wrote", report_sn.render(DB, MANIFEST).relative_to(REPO))
    return 0


COMMANDS["paper"] = cmd_paper
COMMANDS["report"] = cmd_report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch")
    p.add_argument("--refresh", action="store_true")
    sub.add_parser("frames")
    p = sub.add_parser("measure")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int, default=0)
    sub.add_parser("status")
    for name in COMMANDS:
        q = sub.add_parser(name)
        q.add_argument("--workers", type=int, default=10)
    args = ap.parse_args(argv)
    fn = {"fetch": cmd_fetch, "frames": cmd_frames, "measure": cmd_measure,
          "status": cmd_status, **COMMANDS}[args.cmd]
    return fn(args)


if __name__ == "__main__":
    sys.exit(main())
