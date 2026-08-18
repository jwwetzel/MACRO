#!/usr/bin/env python
"""G-track validation runner: T CrB grism series, core machinery pass.

Stages (each independently runnable; ``--all`` chains them):

* ``--calibrate``  solve a handful of era-76 Mode0 IMAGING frames from the
                   grism season with the S1 astrometry wrapper, harvest
                   their CD matrices, and adopt the parity-normalized
                   median as the gate's camera model (g_gate_calib).
* ``--plan``       print the frame worklist (roles and counts), no work.
* ``--run``        process frames (gate -> extract x2 -> wavelength ->
                   FITS product + g_extractions row).  Resumable: frames
                   already in g_extractions are skipped; ``--limit N``
                   bounds one batch so no invocation outruns a shell
                   timeout.
* ``--parquet``    bundle every product FITS into one analysis-ready
                   parquet (adds the per-grism fallback wavelength axis
                   for halpha_only frames).
* ``--report``     render docs/pipeline/g_grism.html from the DB.

The archive is READ-ONLY; outputs go to products/grism/ and docs/.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Make pipeline/ importable when run as a script from the repo root.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import sqlite3

from macro_core.astrom import solve_one_frame            # noqa: E402
from macro_grism import db as gdb                        # noqa: E402
from macro_grism import extract as gext                  # noqa: E402
from macro_grism import gate as ggate                    # noqa: E402
from macro_grism import trace as gtrace                  # noqa: E402
from macro_grism import wavelength as gwave              # noqa: E402
from macro_grism.fits_io import GrismLayoutError, load_frame  # noqa: E402

# ---------------------------------------------------------------------------
# Locations (defaults follow the repo layout; every one is overridable).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
DEFAULT_PRODUCTS = REPO_ROOT / "products" / "grism"
ASTROMETRY_CFG = Path.home() / "astrometry-indices" / "astrometry.cfg"

#: The validation sample sizes (the strategy's numbers).
N_TCRB_SAMPLE = 54          # stratified across the 60 nights, both grisms
N_CALIBRATOR = 10           # tet CrB frames (C4 Be-contamination flagged)
N_CD_SOLVES = 6             # imaging frames solved for the camera model

#: tet CrB rows carry this flag everywhere (strategy ruling C4: the
#: calibrator is a Be/shell star — its Halpha region is astrophysically
#: contaminated and must never be used as a featureless reference there).
C4_FLAG = "C4_Be_Halpha"

#: The grism era under validation (T CrB series is entirely era 76 Mode0).
ERA_ID = 76


# ---------------------------------------------------------------------------
# Camera-model calibration
# ---------------------------------------------------------------------------
def pick_calibration_frames(con) -> list[tuple]:
    """Era-76 Mode0 full-frame IMAGING frames spread across the T CrB
    season (2025-02..06): candidates for the CD solve.  Spread is by
    night percentile so pier-side/season drift is sampled, not assumed."""
    rows = con.execute("""
        SELECT path, night, ra_deg, dec_deg FROM frames
        WHERE is_canonical = 1 AND era_id = ? AND tree = 'rawimage'
          AND lower(coalesce(filter,'')) IN ('g','r','i')
          AND imagetyp LIKE 'Light%' AND naxis1 = 4788
          AND exptime BETWEEN 5 AND 300
          AND night BETWEEN '2025-02-15' AND '2025-06-30'
        ORDER BY night""", (ERA_ID,)).fetchall()
    if not rows:
        return []
    # Even spread: one frame at each of N evenly spaced positions.
    idx = np.unique(np.linspace(0, len(rows) - 1, N_CD_SOLVES).astype(int))
    return [rows[i] for i in idx]


def run_calibrate(con, archive: Path, scratch: Path) -> None:
    """Solve the calibration frames, store every result, adopt the
    parity-normalized element-wise median CD."""
    con.execute("DELETE FROM g_gate_calib")     # a re-calibration replaces
    frames = pick_calibration_frames(con)
    from astropy.io import fits as _fits
    cds = []
    for path, night, ra, dec in frames:
        work = scratch / f"solve_{Path(path).stem}"
        res = solve_one_frame(str(archive / path), str(work),
                              str(ASTROMETRY_CFG), "Mode0", 2, ra, dec)
        cd = [None] * 4
        if res["status"] == "solved":
            wcs_files = glob.glob(str(work / "*.wcs"))
            if wcs_files:
                h = _fits.getheader(wcs_files[0])
                cd = [h["CD1_1"], h["CD1_2"], h["CD2_1"], h["CD2_2"]]
                cds.append(cd)
        con.execute("""
            INSERT INTO g_gate_calib (era_id, frame_path, night, status,
                cd1_1, cd1_2, cd2_1, cd2_2, pixscale_arcsec, rotation_deg,
                rms_arcsec, n_matched, adopted)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,0)""",
            (ERA_ID, path, night, res["status"], *cd,
             res["pixscale_arcsec"], res["rotation_deg"],
             res["rms_arcsec"], res["n_matched"]))
        print(f"  calib {night} {res['status']} "
              f"scale={res['pixscale_arcsec']} rot={res['rotation_deg']}")
    if not cds:
        raise SystemExit("no calibration frame solved — gate cannot run")
    # Parity normalization: a meridian flip negates the whole matrix, so
    # flip every solution onto the sign of the first, then take the
    # element-wise median (robust to one bad solve).
    ref_sign = np.sign(cds[0][0])
    normed = [np.array(c) if np.sign(c[0]) == ref_sign else -np.array(c)
              for c in cds]
    adopted = np.median(np.array(normed), axis=0)
    con.execute("""
        INSERT INTO g_gate_calib (era_id, frame_path, night, status,
            cd1_1, cd1_2, cd2_1, cd2_2, adopted)
        VALUES (?, 'ADOPTED-MEDIAN', NULL, 'adopted', ?,?,?,?, 1)""",
        (ERA_ID, *[float(v) for v in adopted]))
    gdb.set_meta(con, "cd_n_solved", str(len(cds)))
    con.commit()
    print(f"adopted CD (median of {len(cds)}): {adopted.tolist()}")


def adopted_cd(con) -> np.ndarray:
    row = con.execute("""
        SELECT cd1_1, cd1_2, cd2_1, cd2_2 FROM g_gate_calib
        WHERE adopted = 1 ORDER BY calib_id DESC LIMIT 1""").fetchone()
    if row is None:
        raise SystemExit("no adopted CD — run --calibrate first")
    return np.array(row, dtype=float).reshape(2, 2)


# ---------------------------------------------------------------------------
# The worklist: which frames, in which role
# ---------------------------------------------------------------------------
def target_reference(con, like: str) -> tuple[float, float]:
    """A target's reference coordinates: the MEDIAN header pointing of its
    canonical grism series.  The median shrugs off the 21 bad pointings
    (they are < 10% of the series), so this stays honest without using
    the outlier flag it is about to help judge."""
    rows = con.execute("""
        SELECT ra_deg, dec_deg FROM frames
        WHERE is_canonical = 1 AND target_best LIKE ?
          AND lower(coalesce(filter,'')) IN ('hrg','lrg') AND era_id = ?
          AND ra_deg IS NOT NULL""", (like, ERA_ID)).fetchall()
    ra = float(np.median([r[0] for r in rows]))
    dec = float(np.median([r[1] for r in rows]))
    return ra, dec


BASE_COLS = ("obs_rowid, path, filter, night, jd, exptime, era_id, "
             "ra_deg, dec_deg, pointing_offset_deg, target_best")


def build_worklist(con) -> list[dict]:
    """The full validation worklist, deterministic, with roles:

    * gate_bad    — ALL 21 known-bad T CrB pointings (the mandated set).
    * gate_good   — a matched good frame per bad one: same grism filter,
                    nearest night, pointing offset < 0.1 deg.
    * tcrb_sample — ~54 good frames stratified across nights, both grisms
                    (alternating by night index so hrg and lrg both cover
                    the season).
    * calibrator  — 10 tet CrB era-76 frames across its exposure strata.
    """
    def rowdicts(sql, params=()):
        cols = [c.strip() for c in BASE_COLS.split(",")]
        return [dict(zip(cols, r)) for r in con.execute(sql, params)]

    bad = rowdicts(f"""
        SELECT {BASE_COLS} FROM frames
        WHERE is_canonical = 1 AND target_best = 'T CrB'
          AND lower(filter) IN ('hrg','lrg') AND pointing_offset_deg > 1
        ORDER BY night, path""")
    good_all = rowdicts(f"""
        SELECT {BASE_COLS} FROM frames
        WHERE is_canonical = 1 AND target_best = 'T CrB'
          AND lower(filter) IN ('hrg','lrg')
          AND pointing_offset_deg < 0.1
        ORDER BY night, path""")
    work, seen = [], set()

    def add(row, role):
        if row["obs_rowid"] in seen:
            return
        seen.add(row["obs_rowid"])
        row = dict(row)
        row["role"] = role
        work.append(row)

    for r in bad:
        add(r, "gate_bad")
    # Matched good frame per bad frame: same filter, nearest night.
    import datetime as dt

    def night_num(n):
        return dt.date.fromisoformat(n).toordinal()

    for b in bad:
        cands = [g for g in good_all if g["filter"] == b["filter"]
                 and g["obs_rowid"] not in seen]
        if not cands:
            continue
        cands.sort(key=lambda g: (abs(night_num(g["night"])
                                      - night_num(b["night"])), g["path"]))
        add(cands[0], "gate_good")
    # Stratified season sample: group good frames by night, walk nights
    # in order, alternate the preferred grism, one frame per night.
    by_night: dict = {}
    for g in good_all:
        by_night.setdefault(g["night"], []).append(g)
    nights = sorted(by_night)
    picked = 0
    for i, night in enumerate(nights):
        if picked >= N_TCRB_SAMPLE:
            break
        pref = "hrg" if i % 2 == 0 else "lrg"
        cands = ([g for g in by_night[night] if g["filter"] == pref
                  and g["obs_rowid"] not in seen]
                 or [g for g in by_night[night]
                     if g["obs_rowid"] not in seen])
        if cands:
            add(sorted(cands, key=lambda g: g["path"])[0], "tcrb_sample")
            picked += 1
    # tet CrB calibrator sample: spread across its exposure strata (the
    # target_best labels bin the series by exptime); shortest exposures
    # first within each stratum — tet CrB is G ~ 4, saturation is the
    # enemy, and short frames are the usable ones.  Frames whose header
    # pointing is > POINTING_TOL_DEG off the tet CrB reference are
    # EXCLUDED here: the first validation pass proved the March-2025
    # header-pointing bug hit the tet CrB series too (150-162 deg off),
    # and a calibrator sample must contain calibrator spectra — the
    # bad-header frames are gate evidence, not calibration material.
    tet_ref = target_reference(con, "tet CrB%")
    strata = con.execute("""
        SELECT DISTINCT target_best FROM frames
        WHERE is_canonical = 1 AND target_best LIKE 'tet CrB%'
          AND era_id = ? AND lower(filter) IN ('hrg','lrg')
        ORDER BY target_best""", (ERA_ID,)).fetchall()
    per = max(1, N_CALIBRATOR // max(1, len(strata)))
    n_cal = 0
    for (stratum,) in strata:
        rows = rowdicts(f"""
            SELECT {BASE_COLS} FROM frames
            WHERE is_canonical = 1 AND target_best = ? AND era_id = ?
              AND lower(filter) IN ('hrg','lrg')
            ORDER BY exptime, night, path""", (stratum, ERA_ID))
        n_stratum = 0
        for r in rows:
            if n_cal >= N_CALIBRATOR or n_stratum >= per:
                break
            if (r["ra_deg"] is None or ggate.angular_offset_deg(
                    r["ra_deg"], r["dec_deg"], *tet_ref)
                    > ggate.POINTING_TOL_DEG):
                continue                       # header-bug frame: skip
            add(r, "calibrator")
            n_cal += 1
            n_stratum += 1
    return work


# ---------------------------------------------------------------------------
# Master darks (the comparison arm's calibration source)
# ---------------------------------------------------------------------------
class DarkLibrary:
    """Era-76 master darks by exposure time, loaded lazily and kept.

    ``nearest(exptime)`` returns (path, exptime, pixels) for the master
    whose exposure is closest in LOG space — 2 s is 'nearer' to 4 s than
    to 0.5 s where dark current is concerned."""

    def __init__(self, con, archive: Path):
        self.archive = archive
        self.rows = con.execute("""
            SELECT exptime, path FROM calib_frames
            WHERE era_id = ? AND kind = 'dark' AND is_master = 1
              AND exptime > 0 ORDER BY exptime""", (ERA_ID,)).fetchall()
        self._cache: dict = {}

    def nearest(self, exptime: float):
        if not self.rows or not exptime or exptime <= 0:
            return None
        best = min(self.rows,
                   key=lambda r: abs(np.log(r[0]) - np.log(exptime)))
        if best[0] not in self._cache:
            data, _, _ = load_frame(str(self.archive / best[1]))
            self._cache[best[0]] = (best[1], data)
        path, data = self._cache[best[0]]
        return path, float(best[0]), data


# ---------------------------------------------------------------------------
# One frame, end to end
# ---------------------------------------------------------------------------
def process_frame(row: dict, cd: np.ndarray, refs: dict, darks: DarkLibrary,
                  archive: Path, spectra_dir: Path, cache_dir: Path) -> list[dict]:
    """Gate + double extraction + wavelength for one frame.  Returns the
    two g_extractions row dicts (flanking, masterdark).  Any failure is
    captured into ``status`` — a bad frame must produce a record, not a
    crashed batch."""
    base = {k: row.get(k) for k in ("obs_rowid", "path", "filter", "night",
                                    "jd", "exptime", "era_id", "role")}
    base["target"] = row["target_best"]
    if row["target_best"].startswith("tet CrB"):
        base["contamination_flag"] = C4_FLAG
    try:
        data, header, layout = load_frame(str(archive / row["path"]))
    except (GrismLayoutError, OSError) as exc:
        base.update(method="flanking", status=f"load_error: {exc}"[:300])
        return [base]
    base["layout"] = layout
    ny, nx = data.shape

    # ---- trace geometry -------------------------------------------------
    xs, ys, amps = gtrace.chunk_peaks(data)
    slope = gtrace.fit_slope(xs, ys, amps)
    _, resid = gtrace.detilted_profile(data, slope)
    u_obs, height = gtrace.main_trace_u(resid)
    base.update(trace_slope=slope, trace_height=height, u_obs=float(u_obs))

    # ---- identity gate --------------------------------------------------
    ra0, dec0 = row.get("ra_deg"), row.get("dec_deg")
    ref = refs["tet" if row["target_best"].startswith("tet") else "tcrb"]
    offset = (ggate.angular_offset_deg(ra0, dec0, *ref)
              if ra0 is not None and dec0 is not None else None)
    stars = np.empty((0, 3))
    if ra0 is not None:
        try:
            stars = ggate.gaia_cone(ra0, dec0, str(cache_dir))
        except Exception as exc:                     # noqa: BLE001
            base.update(method="flanking",
                        status=f"gaia_error: {exc}"[:300])
            return [base]
    preds = (ggate.brightest_prediction(cd, ra0, dec0, stars, slope, ny, nx)
             if len(stars) else {"A": None, "B": None})
    g = ggate.gate_verdict(offset, height, float(u_obs), preds, len(stars))
    base.update(gate_verdict=g.verdict, gate_reason=g.reason,
                pointing_offset_deg=offset, u_pred=g.u_pred_best,
                u_resid_px=g.u_resid_px, gate_parity=g.parity,
                n_gaia=g.n_gaia, brightest_g=g.brightest_g)

    # A frame with no usable trace cannot be extracted; record and stop.
    if height < ggate.MIN_TRACE_HEIGHT_ADU:
        base.update(method="flanking", status="ok_no_trace")
        return [base]

    # ---- extraction, both background methods ---------------------------
    coeffs, n_cent, t_rms = gtrace.fit_trace_centers(data, slope, u_obs)
    base.update(trace_c0=float(coeffs[0]), trace_c1=float(coeffs[1]),
                trace_c2=float(coeffs[2]), trace_rms_px=t_rms,
                trace_n_centroids=n_cent)
    egain = float(header.get("EGAIN", gext.DEFAULT_EGAIN))
    dark = darks.nearest(row.get("exptime"))
    spectra, rows_out = {}, []
    for method in ("flanking", "masterdark"):
        r = dict(base)
        r["method"] = method
        if method == "masterdark":
            if dark is None:
                r["status"] = "no_master_dark"
                rows_out.append(r)
                continue
            r.update(dark_path=dark[0], dark_exptime=dark[1],
                     bg_method="masterdark+flanking")
            spec = gext.extract_spectrum(data, coeffs, egain, dark=dark[2])
        else:
            r["bg_method"] = "flanking"
            spec = gext.extract_spectrum(data, coeffs, egain)
        spectra[method] = spec
        flux = spec["flux"]
        r.update(n_extracted=spec["n_extracted"],
                 n_sat_cols=int((spec["n_sat"] > 0).sum()),
                 peak_flux=float(np.nanmax(flux)),
                 median_flux=float(np.nanmedian(flux)))
        # ---- wavelength anchors (per method: the debt shows up here too)
        w = gwave.solve_wavelength(flux)
        r.update(anchor_status=w["anchor_status"], x_halpha=w["x_halpha"],
                 halpha_snr=w["halpha_snr"],
                 halpha_width_px=w["halpha_width_px"], x_o2b=w["x_o2b"],
                 o2b_snr=w["o2b_snr"], x_o2a=w["x_o2a"],
                 o2a_snr=w["o2a_snr"], disp_a_per_px=w["disp_a_per_px"],
                 disp_source=w["disp_source"])
        if w["x_halpha"] is not None:
            r["snippet_json"] = json.dumps(
                gwave.snippet(flux, w["x_halpha"]))
        r["status"] = "ok"
        rows_out.append(r)

    # ---- the debt number + the FITS product -----------------------------
    if "flanking" in spectra and "masterdark" in spectra:
        debt = gext.median_relative_difference(
            spectra["flanking"]["flux"], spectra["masterdark"]["flux"])
        for r in rows_out:
            r["debt_median_rel_diff"] = debt
    fits_rel = write_spectrum_fits(row, base, rows_out, spectra,
                                   spectra_dir)
    for r in rows_out:
        r["spectrum_fits"] = fits_rel
    return rows_out


def write_spectrum_fits(row, base, rows_out, spectra, spectra_dir: Path):
    """One FITS per frame: a binary-table extension per method with
    x / wavelength (when anchored) / flux / var / background / n_sat.
    Written atomically (tmp + rename)."""
    if not spectra:
        return None
    from astropy.io import fits as _fits
    hdus = [_fits.PrimaryHDU()]
    hdr0 = hdus[0].header
    hdr0["SRCPATH"] = (row["path"], "archive frame (read-only)")
    hdr0["TARGET"] = row["target_best"]
    hdr0["NIGHT"] = row["night"]
    hdr0["JDSTART"] = (row["jd"], "header JD, UTC exposure START")
    hdr0["GATE"] = (base.get("gate_verdict") or "", "identity gate")
    if base.get("contamination_flag"):
        hdr0["C4FLAG"] = (base["contamination_flag"],
                          "Be/shell calibrator: Halpha contaminated")
    for r in rows_out:
        m = r["method"]
        if m not in spectra:
            continue
        spec = spectra[m]
        nx = len(spec["flux"])
        cols = [
            _fits.Column(name="X_PX", format="J",
                         array=np.arange(nx, dtype=np.int32)),
            _fits.Column(name="FLUX_ADU", format="E", array=spec["flux"]),
            _fits.Column(name="VAR_ADU2", format="E", array=spec["var"]),
            _fits.Column(name="BOX_ADU", format="E", array=spec["box"]),
            _fits.Column(name="BG_ADU", format="E", array=spec["bg"]),
            _fits.Column(name="N_SAT", format="I",
                         array=spec["n_sat"].astype(np.int16)),
        ]
        if r.get("x_halpha") is not None and r.get("disp_a_per_px"):
            cols.insert(1, _fits.Column(
                name="WAVE_A", format="E",
                array=gwave.wavelength_axis(nx, r["x_halpha"],
                                            r["disp_a_per_px"])))
        t = _fits.BinTableHDU.from_columns(cols, name=m.upper())
        t.header["BGMETHOD"] = r.get("bg_method") or ""
        t.header["ANCHORS"] = r.get("anchor_status") or ""
        if r.get("x_halpha") is not None:
            t.header["XHALPHA"] = r["x_halpha"]
        if r.get("disp_a_per_px"):
            t.header["DISPAPX"] = (r["disp_a_per_px"],
                                   "A/px, signed, per-frame O2 anchor")
        hdus.append(t)
    name = Path(row["path"]).stem.replace(".fts", "") + "_spec.fits"
    out = spectra_dir / name
    tmp = out.with_suffix(".tmp")
    _fits.HDUList(hdus).writeto(tmp, overwrite=True)
    os.replace(tmp, out)
    return str(Path("spectra") / name)


# ---------------------------------------------------------------------------
# Parquet bundling
# ---------------------------------------------------------------------------
def run_parquet(con, products: Path) -> None:
    """Long-format parquet over every product FITS.  Adds the per-grism
    median dispersion as the fallback wavelength for halpha_only frames
    (recorded as wave_source='grism_median')."""
    import pandas as pd
    from astropy.io import fits as _fits
    # Per-grism median dispersion, computed in Python (SQLite has no
    # median) from the flanking-method rows that carried an O2 anchor.
    disp_rows = con.execute("""
        SELECT filter, disp_a_per_px FROM g_extractions
        WHERE disp_a_per_px IS NOT NULL AND method = 'flanking'""")
    by_filt: dict = {}
    for filt, disp in disp_rows:
        by_filt.setdefault(filt, []).append(disp)
    med = {filt: float(np.median(v)) for filt, v in by_filt.items()}
    rows = con.execute("""
        SELECT obs_rowid, method, target, filter, night, jd, role,
               gate_verdict, contamination_flag, spectrum_fits,
               anchor_status, x_halpha, disp_a_per_px
        FROM g_extractions
        WHERE spectrum_fits IS NOT NULL AND status = 'ok'""").fetchall()
    frames = []
    for (rowid, method, target, filt, night, jd, role, verdict, c4,
         rel, anchor, x_ha, disp) in rows:
        path = products / rel
        if not path.exists():
            continue
        with _fits.open(path) as h:
            if method.upper() not in h:
                continue
            t = h[method.upper()].data
        df = pd.DataFrame({"x_px": t["X_PX"], "flux_adu": t["FLUX_ADU"],
                           "var_adu2": t["VAR_ADU2"]})
        wave_source = None
        if "WAVE_A" in t.dtype.names:
            df["wave_a"] = t["WAVE_A"]
            wave_source = "frame_o2"
        elif x_ha is not None and filt in med:
            df["wave_a"] = gwave.wavelength_axis(len(df), x_ha, med[filt])
            wave_source = "grism_median"
        df["obs_rowid"] = rowid
        df["method"] = method
        df["target"] = target
        df["filter"] = filt
        df["night"] = night
        df["jd_start_utc"] = jd
        df["role"] = role
        df["gate_verdict"] = verdict
        df["contamination_flag"] = c4
        df["wave_source"] = wave_source
        frames.append(df)
    out = products / "spectra_g_validation.parquet"
    tmp = out.with_suffix(".tmp")
    pd.concat(frames, ignore_index=True).to_parquet(tmp, index=False)
    os.replace(tmp, out)
    print(f"parquet: {out} ({sum(len(f) for f in frames):,} rows, "
          f"{len(frames)} spectra)")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    ap.add_argument("--archive", default=str(DEFAULT_ARCHIVE))
    ap.add_argument("--products", default=str(DEFAULT_PRODUCTS))
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--limit", type=int, default=0,
                    help="max frames this batch (0 = no limit)")
    ap.add_argument("--parquet", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args(argv)

    products = Path(args.products)
    spectra_dir = products / "spectra"
    cache_dir = products / "gaia_cache"
    scratch = products / "scratch"
    for d in (spectra_dir, cache_dir, scratch):
        d.mkdir(parents=True, exist_ok=True)

    con = sqlite3.connect(args.manifest)
    gdb.ensure_schema(con)
    gdb.set_meta(con, "code_version", ggate.G_CODE_VERSION)
    con.commit()

    if args.calibrate or args.all:
        run_calibrate(con, Path(args.archive), scratch)

    work = build_worklist(con)
    if args.plan:
        from collections import Counter
        print(Counter(w["role"] for w in work))
        for w in work:
            print(f"  {w['role']:12s} {w['filter']:3s} {w['night']} "
                  f"{w['path']}")

    if args.run or args.all:
        cd = adopted_cd(con)
        refs = {"tcrb": target_reference(con, "T CrB"),
                "tet": target_reference(con, "tet CrB%")}
        gdb.set_meta(con, "ref_tcrb", json.dumps(refs["tcrb"]))
        gdb.set_meta(con, "ref_tet", json.dumps(refs["tet"]))
        darks = DarkLibrary(con, Path(args.archive))
        done = gdb.existing_keys(con)
        todo = [w for w in work if (w["obs_rowid"], "flanking") not in done]
        if args.limit:
            todo = todo[:args.limit]
        print(f"worklist {len(work)}, remaining {len(todo)}")
        t0 = time.time()
        for i, w in enumerate(todo):
            rows = process_frame(w, cd, refs, darks, Path(args.archive),
                                 spectra_dir, cache_dir)
            for r in rows:
                gdb.insert_extraction(con, r)
            con.commit()                     # one frame = one transaction
            v = rows[0].get("gate_verdict")
            print(f"  [{i + 1}/{len(todo)}] {w['role']:11s} "
                  f"{w['filter']:3s} {w['night']} {v} "
                  f"({time.time() - t0:.0f}s)")
        gdb.set_meta(con, "last_run_utc",
                     time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        con.commit()

    if args.parquet or args.all:
        run_parquet(con, products)

    if args.report or args.all:
        from macro_grism.report_g import render_report
        out = render_report(Path(args.manifest))
        print(f"report: {out}")


if __name__ == "__main__":
    main()
