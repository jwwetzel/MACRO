#!/usr/bin/env python
"""CV-R9, CV-R10, CV-R11 — the reduction, instrument and clock checks.

WHAT THIS SCRIPT DOES, AND WHY EACH PIECE IS HERE
--------------------------------------------------
The plan review of 2026-10-03 (``committee/reviews/2026-10-03/``) asked
the CV paper to describe and TEST its reduction, to state the instrument
on measured numbers, and to print its clock residuals.  The statistics of
the revision live in ``run_cv_revision.py``; this script is the three
tasks that are about the data rather than about the estimators.  Every
stage stores its answer in the photometry products database, in ``rv_``
tables and as ``rv_result`` scalars, so the manuscript quotes a macro.

``reduction``  (R9; OA.E5)  One night per instrument era photometered
               TWICE from the same raw pixels, at identical positions and
               apertures: for the eras the observatory server reduced, the
               server product against raw + this programme's own staged
               masters; for the eras this programme reduced itself (no
               server product exists), the production recipe against the
               same recipe WITHOUT the flat, which is the most the flat can
               be doing.  The two differential light curves of the target
               and of every check star are differenced frame by frame and
               the scatter of that difference is set beside the check
               stars' own scatter that night.

``ramp``       (R9; TE.F5)  The flat-field ramp, tested where it would
               bite: the catalogue-tie residual of every tie star against
               its normalised detector column, per tied block, as an
               edge-to-edge change in per cent with a scatter-based error;
               the radial term the tie stage already fitted
               (``cv_cattie_trend``) is carried beside it.

``mech``       (R9; TE, F-3)  The mechanical epoch of every frame, from the
               manifest's ``mech_epoch`` table; the series that span more
               than one state; and, for each of those, the magnitude step
               of every held-out check star between states, which is what
               a fold across states would inherit.

``clock``      (R11; RF.M4, OA.E6)  The absolute-clock verdicts of the
               transit stage (``products/clock/clock_transits.sqlite``)
               copied into the products database so the paper can print
               them; the one eclipse residual the first draft withheld;
               and, for the season with no clock target of its own, the
               offset of its ST LMi epochs from the line through the
               epochs of the clock-verified seasons.

``detector``   (R10; DE.F2/F3/F7)  The measured detector constants the
               instrument section may quote, copied from the manifest's
               ``detector_params``.  Flat-pair gains are used ONLY if the
               detector package has emitted them; otherwise the stage says
               so and stores nothing it would have had to invent.

``all``        reduction, ramp, mech, clock, detector.
``status``     what is in the database.

TABLES WRITTEN (inside products/phot/cv_timeseries.sqlite)
----------------------------------------------------------
``rv_rawred``         per frame and star: flux under both reductions.
``rv_rawred_summary`` per era night: the agreement statistics.
``rv_ramp``           per tied block: x ramp and the radial swing.
``rv_mech_series``    per series: frames per mechanical epoch.
``rv_mech_step``      per multi-state series and check star: the step.
``rv_clock_era``      the transit clock per camera era, copied.
``rv_result``         scalars, stages ``reduction``, ``ramp``, ``mech``,
                      ``clock``, ``detector``.

Databases other than the products database are opened READ-ONLY.  The
archive is never written.

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_checks.py all
    $P pipeline/scripts/run_cv_checks.py reduction --workers 4
"""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "scripts"))

from macro_phot import revision_checks as rc                    # noqa: E402
from macro_phot import photometry as php                        # noqa: E402
from macro_phot import series as sr                             # noqa: E402
# The revision stage's plumbing (connection, the rv_result writer, the
# build stamp) is imported, not copied: one writer, one schema.
import run_cv_revision as rvs                                   # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = rvs.DEFAULT_DB
MANIFEST_DB = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
CLOCK_DB = REPO_ROOT / "products" / "clock" / "clock_transits.sqlite"
ARCHIVE_ROOT = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

CHECKS_CODE_VERSION = ("CV-R checks v1.0 (2026-10-04, R9 reduction + ramp "
                       "+ mechanical epochs, R10 detector constants, R11 "
                       "clock)")

MAX_WORKERS = 6

#: One night per instrument era, the densest of the series that era's
#: science rests on.  ``mode`` says what the second reduction is:
#: ``server_vs_local`` = server product (A) against raw + staged masters
#: (B); ``flat_vs_noflat`` = the production local recipe (A) against the
#: same without its flat (B), for eras the server never reduced.
#: Eras 6, 78 and 79 are not here: era 6 is 15 YZ Cnc frames that enter no
#: measurement, and no master calibration of the QHY600 eras is staged, so
#: no second reduction can be made of them (the stage records that).
REDUCTION_NIGHTS = (
    ("stlmi|e76|g", "2025-02-27", "server_vs_local"),
    ("vvpup|e72|g", "2024-12-03", "server_vs_local"),
    ("stlmi|e7|R", "2024-03-03", "flat_vs_noflat"),
    ("stlmi|e47|z", "2024-04-10", "flat_vs_noflat"),
)

#: Frames per night photometered twice: evenly spaced through the night.
#: 60 gives the frame-to-frame difference scatter to ~9 per cent of itself.
MAX_FRAMES_PER_NIGHT = 60

#: Fewest tie stars a block needs before its ramp is fitted.
MIN_RAMP_STARS = 20

#: A ramp below this, edge to edge, is "flat" for the paper (the plan's
#: acceptance criterion for CV-R9).
RAMP_FLAT_PCT = 0.5

#: Fewest frames a mechanical state needs before a step into it is measured.
MIN_STATE_FRAMES = 10

#: The transit clock's acceptance bar per era, seconds (F-8).
CLOCK_BAR_S = 120.0


def put(con, stage, key, value, unit, fmt, note, text=None,
        origin="measured") -> None:
    """``run_cv_revision.put`` with the provenance wording the paper's
    lint expects: a constant of method says it was set by CV-R and where."""
    if origin == "constant" and not note.lower().startswith("set by"):
        note = f"set by CV-R (run_cv_checks.py, stage {stage}): {note}"
    rvs.put(con, stage, key, value, unit, fmt, note, text=text,
            origin=origin)


# ===========================================================================
# STAGE: reduction — one night per era, reduced twice
# ===========================================================================
def _local_masters(man: sqlite3.Connection, era: int, filt: str,
                   exptime: float, jd: float) -> tuple:
    """Raw + staged masters for ONE frame of a server-reduced era.

    The dark must match the exposure time (the masters are bias-inclusive
    and may not be scaled); where no dark of that length is staged the
    era's master BIAS is used instead, which at -10 C and these exposure
    times differs from a dark by the dark current alone (S2: the Mode0
    master-dark median is the 303 ADU pedestal).  The flat must match the
    filter exactly.  Nearest in time wins, ties by path — the production
    rule (``run_cv_photometry._pick_frame_masters``).
    """
    rows = man.execute("""SELECT role, jd, abs_path, filter, exptime
        FROM stage_cv_timeseries WHERE CAST(era_id AS INT) = ?
        AND role IN ('master_dark', 'master_flat', 'master_bias')""",
                       (era,)).fetchall()
    darks = [(r[1], r[2]) for r in rows if r[0] == "master_dark"
             and sr.dark_exptime_matches(r[4], exptime)]
    bias = [(r[1], r[2]) for r in rows if r[0] == "master_bias"]
    flats = [(r[1], r[2]) for r in rows if r[0] == "master_flat"
             and r[3] == filt]
    d = sr.pick_master(darks, jd)
    kind = "dark"
    if d is None:
        d = sr.pick_master(bias, jd)
        kind = "bias"
    f = sr.pick_master(flats, jd)
    return ((d[1] if d else None), kind if d else None,
            (f[1] if f else None),
            (abs(float(jd) - float(f[0])) if f and f[0] else None))


def _measure_pair(job: dict) -> dict:
    """Photometer ONE frame under two reductions.  Worker; no database.

    Both images get the production treatment — a ``sep`` background model
    subtracted, then circular apertures of the frame's own radius with the
    production sky annulus — at the SAME pixel positions (the production
    detections), so the two flux vectors differ only by the reduction.
    """
    import sep
    from macro_phot import calib as cal
    from macro_phot import extract as ext
    from macro_phot import photometry as ph
    out = {"frame_id": job["frame_id"], "ok": False}
    try:
        raw, _ = ext.read_reduced(Path(job["raw_abs"]))
        if job["mode"] == "server_vs_local":
            a, _ = ext.read_reduced(Path(job["pix_abs"]))
            dark = cal.read_master(Path(job["dark"])) if job["dark"] else None
            flat = cal.read_master(Path(job["flat"])) if job["flat"] else None
            b, recipe_b = cal.apply_masters(raw, dark, flat)
            recipe_a = "server_reduced"
        else:
            dark = cal.read_master(Path(job["dark"])) if job["dark"] else None
            flat = cal.read_master(Path(job["flat"])) if job["flat"] else None
            a, recipe_a = cal.apply_masters(raw, dark, flat)
            b, recipe_b = cal.apply_masters(raw, dark, None)
        x = np.asarray(job["x"], dtype=float)
        y = np.asarray(job["y"], dtype=float)
        r = float(job["aper_px"])
        ann = (ph.SKY_ANNULUS_ARCSEC[0] / ph.APERTURE_RADIUS_ARCSEC * r,
               ph.SKY_ANNULUS_ARCSEC[1] / ph.APERTURE_RADIUS_ARCSEC * r)
        sep.set_extract_pixstack(ext.SEP_PIXSTACK)
        fl = []
        for img in (a, b):
            img = np.ascontiguousarray(img, dtype=np.float32)
            # Dead flat pixels are NaN after local flat-fielding; sep needs
            # finite pixels, and a mask keeps them out of the background.
            mask = ~np.isfinite(img)
            img = np.where(mask, 0.0, img).astype(np.float32)
            bkg = sep.Background(img, mask=mask)
            sub = img - bkg.back()
            f, _, flag = sep.sum_circle(sub, x, y, r, bkgann=ann,
                                        mask=mask)
            f = np.where(flag & sep.APER_HASMASKED, np.nan, f)
            fl.append(f)
        out.update(ok=True, flux_a=fl[0].tolist(), flux_b=fl[1].tolist(),
                   recipe_a=recipe_a, recipe_b=recipe_b)
    except Exception as e:                                   # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"[:200]
    return out


def cmd_reduction(args) -> None:
    """Photometer one night per era twice and compare the light curves."""
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    con.executescript("""
    DROP TABLE IF EXISTS rv_rawred;
    CREATE TABLE rv_rawred (series_key TEXT, night TEXT, frame_id INTEGER,
        bjd_tdb REAL, star_id INTEGER, role TEXT, flux_a REAL, flux_b REAL,
        PRIMARY KEY (series_key, night, frame_id, star_id));
    DROP TABLE IF EXISTS rv_rawred_summary;
    CREATE TABLE rv_rawred_summary (series_key TEXT, night TEXT,
        era_id INTEGER, mode TEXT, recipe_a TEXT, recipe_b TEXT,
        dark_kind TEXT, flat_path TEXT, flat_age_d REAL,
        n_frames INTEGER, n_comp INTEGER, n_check INTEGER,
        target_rms_mag REAL, target_mad_mag REAL, target_offset_mag REAL,
        target_max_abs_mag REAL, check_rms_med_mag REAL,
        check_rms_max_mag REAL, sigma_chk_night_mag REAL,
        sigma_chk_series_mag REAL, agree INTEGER, note TEXT,
        PRIMARY KEY (series_key, night));
    """)
    man = sqlite3.connect(f"file:{MANIFEST_DB}?mode=ro", uri=True)
    rvs.clear_stage(con, "reduction")
    agree_all = []
    for skey, night, mode in REDUCTION_NIGHTS:
        era = int(skey.split("|")[1][1:])
        filt = skey.split("|")[2]
        frames = con.execute("""SELECT frame_id, bjd_tdb, exptime, pixel_path,
            raw_path, master_dark, master_flat, aper_px, jd_header
            FROM cv_frames WHERE series_key=? AND night=? AND
            status='matched' ORDER BY bjd_tdb""", (skey, night)).fetchall()
        if len(frames) > MAX_FRAMES_PER_NIGHT:
            idx = np.unique(np.linspace(0, len(frames) - 1,
                                        MAX_FRAMES_PER_NIGHT).astype(int))
            frames = [frames[i] for i in idx]
        roles = {r[0]: r[1] for r in con.execute(
            "SELECT star_id, role FROM cv_stars WHERE series_key=? AND role "
            "IN ('target','comp','check')", (skey,))}
        jobs, dark_kind, flat_used, flat_age = [], None, None, None
        for fr in frames:
            dets = con.execute("""SELECT star_id, x, y FROM cv_detections
                WHERE frame_id=? AND star_id IS NOT NULL AND
                COALESCE(saturated,0)=0""", (fr["frame_id"],)).fetchall()
            dets = [d for d in dets if d[0] in roles]
            if mode == "server_vs_local":
                dk, dark_kind, fl, flat_age = _local_masters(
                    man, era, filt, fr["exptime"], fr["jd_header"])
                flat_used = fl
            else:
                dk, fl = fr["master_dark"], fr["master_flat"]
                dark_kind = "dark" if dk else None
                flat_used = fl
            jobs.append({"frame_id": fr["frame_id"], "mode": mode,
                         "pix_abs": str(ARCHIVE_ROOT / fr["pixel_path"]),
                         "raw_abs": str(ARCHIVE_ROOT / fr["raw_path"]),
                         "dark": dk, "flat": fl, "aper_px": fr["aper_px"],
                         "sid": [d[0] for d in dets],
                         "x": [d[1] for d in dets],
                         "y": [d[2] for d in dets],
                         "bjd": fr["bjd_tdb"]})
        print(f"  {skey} {night} ({mode}): {len(jobs)} frames, flat "
              f"{Path(flat_used).name if flat_used else 'NONE'}",
              flush=True)
        if flat_used is None:
            con.execute("""INSERT OR REPLACE INTO rv_rawred_summary
                (series_key, night, era_id, mode, n_frames, agree, note)
                VALUES (?,?,?,?,?,?,?)""", (skey, night, era, mode, 0, None,
                "no staged master flat for this filter and era"))
            continue
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            res = list(ex.map(_measure_pair, jobs))
        rows, recipes = [], set()
        for job, r in zip(jobs, res):
            if not r["ok"]:
                print(f"    ! frame {job['frame_id']}: {r.get('error')}")
                continue
            recipes.add((r["recipe_a"], r["recipe_b"]))
            for sid, fa, fb in zip(job["sid"], r["flux_a"], r["flux_b"]):
                rows.append((skey, night, job["frame_id"], job["bjd"], sid,
                             roles[sid], rvs._nz(fa), rvs._nz(fb)))
        con.executemany("INSERT OR REPLACE INTO rv_rawred VALUES "
                        "(?,?,?,?,?,?,?,?)", rows)
        # ---- light curves under both reductions -------------------------
        fids = sorted({r[2] for r in rows})
        sids = sorted({r[4] for r in rows})
        fi = {f: i for i, f in enumerate(fids)}
        si = {s: j for j, s in enumerate(sids)}
        FA = np.full((len(fids), len(sids)), np.nan)
        FB = np.full_like(FA, np.nan)
        for r in rows:
            FA[fi[r[2]], si[r[4]]] = np.nan if r[6] is None else r[6]
            FB[fi[r[2]], si[r[4]]] = np.nan if r[7] is None else r[7]
        comp = [si[s] for s in sids if roles[s] == "comp"]
        # A comparison star enters the ensemble only if it is measured with
        # positive flux on EVERY frame under BOTH reductions, so the two
        # ensembles are the same stars and a missing star cannot masquerade
        # as a reduction difference.
        comp = [j for j in comp if np.all(np.isfinite(FA[:, j]) &
                                          (FA[:, j] > 0) &
                                          np.isfinite(FB[:, j]) &
                                          (FB[:, j] > 0))]
        tgt = [si[s] for s in sids if roles[s] == "target"]
        chk = [si[s] for s in sids if roles[s] == "check"]
        ca, cb = FA[:, comp], FB[:, comp]
        stat_t = {"rms": np.nan, "mad_sigma": np.nan, "offset": np.nan,
                  "max_abs": np.nan}
        if tgt:
            stat_t = rc.compare_light_curves(
                rc.differential_mags(FA[:, tgt[0]], ca),
                rc.differential_mags(FB[:, tgt[0]], cb))
        chk_rms, chk_scatter = [], []
        for j in chk:
            da = rc.differential_mags(FA[:, j], ca)
            db = rc.differential_mags(FB[:, j], cb)
            chk_rms.append(rc.compare_light_curves(da, db)["rms"])
            chk_scatter.append(rc.scatter_about_median(da))
        sig_night = float(np.nanmedian(chk_scatter)) if chk_scatter \
            else np.nan
        sig_series = con.execute("SELECT check_rms_median FROM cv_series "
                                 "WHERE series_key=?", (skey,)).fetchone()[0]
        agree = int(np.isfinite(stat_t["rms"]) and stat_t["rms"] < sig_night)
        agree_all.append(agree)
        ra, rb = sorted(recipes)[0] if recipes else (None, None)
        con.execute("""INSERT OR REPLACE INTO rv_rawred_summary VALUES
            (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
            skey, night, era, mode, ra, rb, dark_kind, flat_used,
            rvs._nz(flat_age), len(fids), len(comp), len(chk),
            rvs._nz(stat_t["rms"]), rvs._nz(stat_t["mad_sigma"]),
            rvs._nz(stat_t["offset"]), rvs._nz(stat_t["max_abs"]),
            rvs._nz(np.nanmedian(chk_rms)) if chk_rms else None,
            rvs._nz(np.nanmax(chk_rms)) if chk_rms else None,
            rvs._nz(sig_night), sig_series, agree,
            "A = " + str(ra) + ", B = " + str(rb)))
        tag = {76: "SevenSix", 72: "SevenTwo", 7: "Seven",
               47: "FourSeven"}[era]
        key = f"rv rawred era {tag.lower()}"
        nt = f"{skey} on {night}, {len(fids)} frames, {mode}"
        put(con, "reduction", f"{key} frames", len(fids), "", "int", nt)
        put(con, "reduction", f"{key} target rms mmag",
                1e3 * stat_t["rms"], "mmag", "f1",
                f"{nt}: RMS of the frame-by-frame difference of the two "
                f"target differential light curves, median removed")
        if chk_rms:
            put(con, "reduction", f"{key} check rms mmag",
                    1e3 * float(np.nanmedian(chk_rms)), "mmag", "f1",
                    f"{nt}: the same for the check stars, median over them")
        put(con, "reduction", f"{key} sigma chk mmag",
                1e3 * sig_night, "mmag", "f0",
                f"{nt}: check-star scatter that night under the production "
                f"reduction, median over the check stars")
        put(con, "reduction", f"{key} offset mmag",
                1e3 * stat_t["offset"], "mmag", "f1",
                f"{nt}: median target difference (a zero point; it cancels "
                f"in differential photometry and is not a disagreement)")
        if flat_age is not None:
            put(con, "reduction", f"{key} flat age d", flat_age, "d",
                    "f0", f"{nt}: age of the staged flat used for B")
        print(f"    target diff rms {1e3 * stat_t['rms']:.1f} mmag vs check "
              f"scatter {1e3 * sig_night:.1f} mmag -> "
              f"{'AGREE' if agree else 'DISAGREE'}", flush=True)
        con.commit()
    put(con, "reduction", "rv rawred nights", len(agree_all), "", "int",
            "era nights reduced twice and compared")
    put(con, "reduction", "rv rawred nights agree", sum(agree_all), "",
            "int", "of those, nights whose two target light curves differ by "
            "less than that night's check-star scatter")
    put(con, "reduction", "rv aperture radius arcsec",
            php.APERTURE_RADIUS_ARCSEC, "arcsec", "f0",
            "fixed sky aperture radius (macro_phot.photometry)",
            origin="constant")
    put(con, "reduction", "rv annulus inner arcsec",
            php.SKY_ANNULUS_ARCSEC[0], "arcsec", "f0",
            "sky annulus inner radius (macro_phot.photometry)",
            origin="constant")
    put(con, "reduction", "rv annulus outer arcsec",
            php.SKY_ANNULUS_ARCSEC[1], "arcsec", "f0",
            "sky annulus outer radius (macro_phot.photometry)",
            origin="constant")
    # The frame census by reduction route (all five targets, measured
    # frames): what the reduction paragraph states.
    for prov, key in (("server_reduced", "server"),
                      ("local_master", "local")):
        n = con.execute("SELECT count(*) FROM cv_frames WHERE status="
                        "'matched' AND provenance=?", (prov,)).fetchone()[0]
        put(con, "reduction", f"rv frames {key} reduced", n, "", "int",
                f"measured frames whose pixels are {prov}")
    lo, hi = con.execute("SELECT min(flat_age_days), max(flat_age_days) FROM "
                         "cv_frames WHERE status='matched' AND provenance="
                         "'local_master' AND flat_age_days IS NOT NULL"
                         ).fetchone()
    put(con, "reduction", "rv local flat age min d", lo, "d", "f0",
            "youngest master flat applied to a locally reduced frame")
    put(con, "reduction", "rv local flat age max d", hi, "d", "f0",
            "oldest master flat applied to a locally reduced frame")
    rvs.stamp(con, "reduction")
    con.commit()
    con.close()


# ===========================================================================
# STAGE: ramp — tie residual against detector column
# ===========================================================================
def cmd_ramp(args) -> None:
    """The flat-field ramp as catalogue-tie residual against detector x."""
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    con.executescript("""
    DROP TABLE IF EXISTS rv_ramp;
    CREATE TABLE rv_ramp (series_key TEXT PRIMARY KEY, catalogue TEXT,
        band TEXT, n INTEGER, naxis1 INTEGER, ramp_pct REAL,
        ramp_err_pct REAL, nsigma REAL, scatter_mag REAL,
        radial_swing_mag REAL, radial_signif REAL);
    """)
    rvs.clear_stage(con, "ramp")
    blocks = con.execute("""SELECT series_key, catalogue, band FROM cv_cattie
        WHERE is_primary=1 AND n_fit >= ? ORDER BY series_key""",
                         (MIN_RAMP_STARS,)).fetchall()
    out = []
    for skey, cat, band in blocks:
        n1 = con.execute("SELECT max(naxis1) FROM cv_frames WHERE "
                         "series_key=?", (skey,)).fetchone()[0]
        pts = con.execute("""SELECT x, resid FROM cv_cattie_star WHERE
            series_key=? AND catalogue=? AND band=? AND in_fit=1""",
                          (skey, cat, band)).fetchall()
        x = np.array([p[0] for p in pts], dtype=float) / float(n1)
        y = np.array([p[1] for p in pts], dtype=float)
        f = rc.ramp_fit(x, y)
        rad = con.execute("""SELECT swing, significance FROM cv_cattie_trend
            WHERE series_key=? AND catalogue=? AND band=? AND
            axis='radius'""", (skey, cat, band)).fetchone()
        out.append((skey, cat, band, f["n"], n1, f["ramp_pct"],
                    f["ramp_err_pct"], f["nsigma"], f["scatter_mag"],
                    rad[0] if rad else None, rad[1] if rad else None))
    con.executemany("INSERT INTO rv_ramp VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    [tuple(rvs._nz(v) if isinstance(v, float) else v
                           for v in r) for r in out])
    ramps = np.array([r[5] for r in out], dtype=float)
    errs = np.array([r[6] for r in out], dtype=float)
    nsig = np.array([r[7] for r in out], dtype=float)
    absr = np.abs(ramps)
    put(con, "ramp", "rv ramp blocks", len(out), "", "int",
            f"primary tied blocks with at least {MIN_RAMP_STARS} tie stars")
    put(con, "ramp", "rv ramp flat bar pct", RAMP_FLAT_PCT, "%", "f1",
            "edge-to-edge ramp called flat (plan acceptance, CV-R9)",
            origin="constant")
    put(con, "ramp", "rv ramp below bar", int(np.sum(absr < RAMP_FLAT_PCT)),
            "", "int", f"blocks whose fitted x ramp is below "
            f"{RAMP_FLAT_PCT} per cent edge to edge")
    put(con, "ramp", "rv ramp consistent with bar",
            int(np.sum(absr - 2 * errs < RAMP_FLAT_PCT)), "", "int",
            f"blocks whose ramp is within 2 sigma of below "
            f"{RAMP_FLAT_PCT} per cent")
    put(con, "ramp", "rv ramp significant", int(np.sum(nsig >= 3)), "",
            "int", "blocks whose x ramp differs from zero at >= 3 sigma")
    put(con, "ramp", "rv ramp median abs pct", float(np.median(absr)),
            "%", "f1", "median |edge-to-edge x ramp| over the blocks")
    k = int(np.argmax(absr))
    put(con, "ramp", "rv ramp max abs pct", absr[k], "%", "f1",
            f"largest |x ramp|, block {out[k][0]}")
    put(con, "ramp", "rv ramp max abs err pct", errs[k], "%", "f1",
            f"its scatter-based error, block {out[k][0]}")
    put(con, "ramp", "rv ramp max block", None, "", "text",
            "block with the largest |x ramp|",
            text=_block_label(out[k][0]))
    for skey, *_rest in out:
        if not skey.startswith("stlmi|"):
            continue
        r = next(o for o in out if o[0] == skey)
        tag = skey.split("|")[1][1:] + skey.split("|")[2]
        tag = tag.replace("76", "SevenSix ").replace("47", "FourSeven ") \
                 .replace("7", "Seven ")
        put(con, "ramp", f"rv ramp stlmi {tag} pct", r[5], "%", "f1",
                f"ST LMi {skey}: edge-to-edge x ramp of the tie residuals")
        put(con, "ramp", f"rv ramp stlmi {tag} err pct", r[6], "%", "f1",
                f"ST LMi {skey}: its scatter-based error")
        if r[9] is not None:
            put(con, "ramp", f"rv ramp stlmi {tag} radial mmag",
                    1e3 * abs(r[9]), "mmag", "f0",
                    f"ST LMi {skey}: |centre-to-corner swing| of the tie "
                    f"residuals fitted against radius (cv_cattie_trend)")
    # Per camera era: the range of the x ramp and of the radial swing over
    # that era's blocks, which is how the photometry section states them.
    era_tag = {"7": "seven", "47": "fourseven", "72": "seventwo",
               "76": "sevensix"}
    for e, tag in era_tag.items():
        sel = [o for o in out if o[0].split("|")[1] == f"e{e}"]
        if not sel:
            continue
        rr = np.array([o[5] for o in sel], dtype=float)
        put(con, "ramp", f"rv ramp era {tag} min pct", float(rr.min()),
                "%", "f1", f"era {e}: smallest signed x ramp over its "
                f"{len(sel)} blocks")
        put(con, "ramp", f"rv ramp era {tag} max pct", float(rr.max()),
                "%", "f1", f"era {e}: largest signed x ramp over its blocks")
        put(con, "ramp", f"rv ramp era {tag} blocks", len(sel), "",
                "int", f"era {e}: tied blocks")
        rs = np.array([abs(o[9]) for o in sel if o[9] is not None],
                      dtype=float)
        if rs.size:
            put(con, "ramp", f"rv radial era {tag} max mmag",
                    1e3 * float(rs.max()), "mmag", "f0",
                    f"era {e}: largest |centre-to-corner swing|")
    rad = np.array([abs(r[9]) for r in out if r[9] is not None], dtype=float)
    rads = np.array([r[10] for r in out if r[9] is not None], dtype=float)
    put(con, "ramp", "rv radial significant", int(np.sum(rads >= 3)), "",
            "int", "blocks whose radial tie-residual trend is >= 3 sigma")
    put(con, "ramp", "rv radial max mmag", 1e3 * float(np.max(rad)),
            "mmag", "f0", "largest |centre-to-corner swing| of the tie "
            "residuals over the blocks")
    put(con, "ramp", "rv radial median mmag", 1e3 * float(np.median(rad)),
            "mmag", "f0", "median |centre-to-corner swing| over the blocks")
    rvs.stamp(con, "ramp")
    con.commit()
    print(f"  {len(out)} blocks; |ramp| median {np.median(absr):.2f}%, "
          f"max {absr[k]:.2f} +/- {errs[k]:.2f}% ({out[k][0]}); "
          f"{int(np.sum(absr < RAMP_FLAT_PCT))} below {RAMP_FLAT_PCT}%; "
          f"{int(np.sum(nsig >= 3))} at >=3 sigma; radial max "
          f"{1e3 * np.max(rad):.0f} mmag")
    con.close()


_TARGET_TEX = {"stlmi": "ST~LMi", "vvpup": "VV~Pup", "euuma": "EU~UMa",
               "anuma": "AN~UMa", "yzcnc": "YZ~Cnc"}


def _block_label(skey: str) -> str:
    t, e, f = skey.split("|")
    return f"{_TARGET_TEX.get(t, t)} era~{e[1:]} ${f}$"


# ===========================================================================
# STAGE: mech — mechanical epoch per series, and the step between states
# ===========================================================================
def cmd_mech(args) -> None:
    """Mechanical epochs per series and check-star steps between them."""
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    con.executescript("""
    DROP TABLE IF EXISTS rv_mech_series;
    CREATE TABLE rv_mech_series (series_key TEXT PRIMARY KEY,
        n_states INTEGER, states_json TEXT);
    DROP TABLE IF EXISTS rv_mech_step;
    CREATE TABLE rv_mech_step (series_key TEXT, star_id INTEGER,
        state_ref TEXT, state TEXT, n_ref INTEGER, n_state INTEGER,
        step_mag REAL, step_err_mag REAL,
        PRIMARY KEY (series_key, star_id, state));
    """)
    rvs.clear_stage(con, "mech")
    man = sqlite3.connect(f"file:{MANIFEST_DB}?mode=ro", uri=True)
    ranges = [(r[0], r[1], r[2]) for r in man.execute(
        "SELECT mech_epoch, first_night, last_night FROM mech_epoch "
        "ORDER BY seq")]
    series = [r[0] for r in con.execute(
        "SELECT DISTINCT series_key FROM cv_frames WHERE status='matched'")]
    multi, worst = [], []
    stl = {}
    for skey in sorted(series):
        fr = con.execute("SELECT frame_id, night FROM cv_frames WHERE "
                         "series_key=? AND status='matched'",
                         (skey,)).fetchall()
        lab = rc.nights_in_ranges([f[1] for f in fr], ranges)
        state_of = {f[0]: s for f, s in zip(fr, lab)}
        counts: dict = {}
        for s in lab:
            counts[s or "none"] = counts.get(s or "none", 0) + 1
        con.execute("INSERT INTO rv_mech_series VALUES (?,?,?)",
                    (skey, len(counts), json.dumps(counts, sort_keys=True)))
        big = [s for s, n in counts.items() if n >= MIN_STATE_FRAMES]
        if len(big) < 2:
            continue
        multi.append(skey)
        ref = max(big, key=lambda s: counts[s])
        steps = []
        for (sid,) in con.execute("SELECT star_id FROM cv_stars WHERE "
                                  "series_key=? AND role='check'", (skey,)):
            pts = con.execute("SELECT frame_id, mag FROM cv_lightcurve WHERE "
                              "series_key=? AND star_id=? AND mag IS NOT "
                              "NULL", (skey, sid)).fetchall()
            by: dict = {}
            for fid, m in pts:
                by.setdefault(state_of.get(fid), []).append(m)
            for s in big:
                if s == ref:
                    continue
                st = rc.state_offset(by.get(ref, []), by.get(s, []))
                con.execute("INSERT INTO rv_mech_step VALUES "
                            "(?,?,?,?,?,?,?,?)",
                            (skey, sid, ref, s, st["n_a"], st["n_b"],
                             rvs._nz(st["step"]), rvs._nz(st["step_err"])))
                if np.isfinite(st["step"]):
                    steps.append(abs(st["step"]))
        if steps:
            worst.append((float(np.median(steps)), skey))
            if skey.startswith("stlmi|"):
                stl[skey] = (float(np.median(steps)), len(big))
    put(con, "mech", "rv mech series", len(series), "", "int",
            "measured series assigned a mechanical epoch per frame")
    put(con, "mech", "rv mech multi state", len(multi), "", "int",
            f"series with >= {MIN_STATE_FRAMES} frames in each of two or "
            f"more mechanical states")
    put(con, "mech", "rv mech multi state list", None, "", "text",
            "the multi-state series, named", text=", ".join(_block_label(s) for s in multi))
    if worst:
        w = max(worst)
        put(con, "mech", "rv mech step max mmag", 1e3 * w[0], "mmag",
                "f0", f"largest median |check-star step| between states, "
                f"over the multi-state series ({w[1]})")
        put(con, "mech", "rv mech step median mmag",
                1e3 * float(np.median([x[0] for x in worst])), "mmag", "f0",
                "median over the multi-state series of the median "
                "|check-star step|")
    for skey, (v, ns) in stl.items():
        f = skey.split("|")[2]
        put(con, "mech", f"rv mech stlmi {f} step mmag", 1e3 * v,
                "mmag", "f0", f"ST LMi {skey}: median |check-star step| "
                f"between its {ns} mechanical states")
    # ST LMi nights by state: which of the timed and colour nights are
    # after the 2025-10 meridian-flip reconfiguration.
    nights = [r[0] for r in con.execute(
        "SELECT DISTINCT night FROM cv_frames WHERE target_key='stlmi' AND "
        "era_id=76 AND status='matched'")]
    lab = rc.nights_in_ranges(nights, ranges)
    states76 = sorted(set(lab))
    put(con, "mech", "rv mech stlmi era seven six states",
            len(states76), "", "int",
            "mechanical states inside ST LMi's era-76 (Mode0) season")
    late = [n for n, s in zip(nights, lab) if s and s.startswith("ASI:2025-10")]
    put(con, "mech", "rv mech stlmi late nights", len(late), "", "int",
            "ST LMi era-76 nights in the post-2025-10 mechanical state")
    rvs.stamp(con, "mech")
    con.commit()
    print(f"  {len(series)} series; {len(multi)} multi-state: "
          f"{', '.join(multi)}")
    for v, s in sorted(worst):
        print(f"    {s}: median |check step| {1e3 * v:.1f} mmag")
    con.close()


# ===========================================================================
# STAGE: clock — the transit clock, printed
# ===========================================================================
def cmd_clock(args) -> None:
    """Copy the per-era clock verdicts and test the untested season."""
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    con.executescript("""
    DROP TABLE IF EXISTS rv_clock_era;
    CREATE TABLE rv_clock_era (clock_era TEXT PRIMARY KEY, label TEXT,
        n_series INTEGER, oc_s REAL, err_s REAL, verdict TEXT);
    """)
    rvs.clear_stage(con, "clock")
    ck = sqlite3.connect(f"file:{CLOCK_DB}?mode=ro", uri=True)
    ck.row_factory = sqlite3.Row
    src = "products/clock/clock_transits.sqlite"
    names = {"A": "A", "B": "B", "C": "C", "D": "D", "E": "E", "F": "F",
             "G": "G"}
    tested = 0
    passed = 0
    for r in ck.execute("SELECT * FROM s3b_era ORDER BY clock_era"):
        con.execute("INSERT INTO rv_clock_era VALUES (?,?,?,?,?,?)",
                    (r["clock_era"], r["label"], r["n_series"], r["oc_s"],
                     r["sig_adopted_s"], r["verdict"]))
        if r["oc_s"] is None:
            continue
        tested += 1
        passed += int(abs(r["oc_s"]) < CLOCK_BAR_S)
        e = names[r["clock_era"]].lower()
        put(con, "clock", f"rv clock era {e} oc s", r["oc_s"], "s", "f0",
                f"transit clock O-C, camera era {r['clock_era']} "
                f"({r['label']}); {src} s3b_era")
        put(con, "clock", f"rv clock era {e} err s", r["sig_adopted_s"],
                "s", "f0", f"its adopted error; {src} s3b_era")
        put(con, "clock", f"rv clock era {e} series", r["n_series"], "",
                "int", f"transit light curves behind it; {src} s3b_era")
    put(con, "clock", "rv clock eras tested", tested, "", "int",
            f"camera eras with at least one timed transit; {src}")
    put(con, "clock", "rv clock eras pass", passed, "", "int",
            f"of those, eras whose central |O-C| < {CLOCK_BAR_S:g} s")
    put(con, "clock", "rv clock bar s", CLOCK_BAR_S, "s", "f0",
            "per-era clock acceptance bar (F-8)", origin="constant")
    leg = {r["key"]: float(r["value"]) for r in ck.execute(
        "SELECT key, value FROM s3b_legacy") if _isnum(r["value"])}
    for k, key, fmt, unit, note in (
            ("agl_oc_s", "rv clock eclipse oc s", "f0", "s",
             "O-C of the one usable eclipse of the detached binary "
             "observed 2024-02-22, against the catalogue ephemeris"),
            ("agl_oc_err_s", "rv clock eclipse err s", "f0", "s",
             "its error"),
            ("agl_ephemeris_part_s", "rv clock eclipse eph s", "f0", "s",
             "the part of that O-C the catalogue ephemeris explains when "
             "the Gaia period is used instead"),
            ("agl_ephemeris_part_err_s", "rv clock eclipse eph err s", "f0",
             "s", "its error"),
            ("agl_implied_minus_gaia_sigma", "rv clock eclipse gaia sigma",
             "f1", "", "the period the eclipse implies, minus the Gaia "
             "period, in units of the combined error"),
            ("agl_cycles", "rv clock eclipse cycles", "int", "",
             "cycles from the catalogue epoch to the eclipse")):
        if k in leg:
            put(con, "clock", key, leg[k], unit, fmt,
                    f"{note}; {src} s3b_legacy")
    for r in ck.execute("SELECT * FROM s3b_bounds"):
        term = r["term"]
        if term.startswith("rolling shutter"):
            cam = ("AC4040" if "AC4040" in r["applies_to"] else
                   "ASI" if "ASI6200" in r["applies_to"] else
                   "QHY" if "QHY600" in r["applies_to"] else None)
            if cam:
                put(con, "clock", f"rv rolling {cam.lower()} s",
                        r["bound_s"], "s", "f2",
                        f"rolling-shutter row-skew bound, {r['applies_to']}: "
                        f"{r['basis']}; {src} s3b_bounds")
        elif term.startswith("mechanical shutter"):
            put(con, "clock", "rv shutter ikon s", r["bound_s"], "s",
                    "f1", f"iKon shutter travel bound (ASSUMED class "
                    f"figure): {r['basis']}; {src} s3b_bounds")
    for r in ck.execute("SELECT * FROM s3b_cv WHERE target_key='euuma'"):
        put(con, "clock", "rv clock euuma offset s", r["cv_offset_s"],
                "s", "f0", f"EU UMa edge offset from its catalogue epoch, "
                f"clock era {r['clock_era']}; {src} s3b_cv")
    # ---- the season with no clock target --------------------------------
    eras = [(r["clock_era"], r["first_night"], r["last_night"]) for r in
            ck.execute("SELECT clock_era, first_night, last_night FROM "
                       "s3b_cv WHERE target_key='stlmi'")]
    verified = {r["clock_era"] for r in ck.execute(
        "SELECT clock_era FROM s3b_era WHERE oc_s IS NOT NULL")}
    for variant, tag in (("v1/pub", "v one"), ("v2/scatter", "v two")):
        ep = con.execute("SELECT night, cycle, oc_s, sigma_s FROM "
                         "rv_oc_epoch WHERE variant=? ORDER BY night",
                         (variant,)).fetchall()
        lab = rc.nights_in_ranges([e[0] for e in ep], eras)
        keep = [i for i, s in enumerate(lab) if s is not None]
        res = rc.interpolated_offset(
            [ep[i][1] for i in keep], [ep[i][2] for i in keep],
            [ep[i][3] for i in keep], [lab[i] in verified for i in keep])
        nv = res["n_verified"]
        put(con, "clock", f"rv clock eraD {tag} offset s",
                res["offset_s"], "s", "f0",
                f"ST LMi ({variant}): weighted mean O-C of the "
                f"{res['n_unverified']} night epochs of clock era D (no "
                f"clock target) about the linear ephemeris through the "
                f"{nv} epochs of clock-verified eras")
        put(con, "clock", f"rv clock eraD {tag} offset err s",
                res["offset_err_s"], "s", "f0", "its error (scatter of the "
                "era-D residuals and the line's uncertainty)")
        put(con, "clock", f"rv clock eraD {tag} nights",
                res["n_unverified"], "", "int", "era-D nights in it")
        put(con, "clock", f"rv clock eraD {tag} verified nights", nv, "",
                "int", "clock-verified nights the line is drawn through")
        print(f"  {variant}: era-D offset {res['offset_s']:+.0f} +/- "
              f"{res['offset_err_s']:.0f} s ({res['n_unverified']} nights "
              f"vs {nv} verified)")
    rvs.stamp(con, "clock")
    con.commit()
    con.close()
    print(f"  {tested} clock eras tested, {passed} pass")


def _isnum(v) -> bool:
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


# ===========================================================================
# STAGE: detector — measured constants for the instrument section
# ===========================================================================
#: Readout modes of the CV photometry -> macro tag.
DETECTOR_MODES = {"High Gain": "hg", "Mode0": "modezero",
                  "1MHz High Sensitivity 16-bit": "onemhz", "Fast": "fast",
                  "High Gain StackPro": "stackpro"}

#: The S2 linearity-cap row that applies to the CV frames of each mode
#: (High Gain is measured per EGAIN epoch; all CV High Gain frames are in
#: the 1.057 epoch).
CAP_ROW = {"High Gain": "AC4040 High Gain e1.057"}


def cmd_detector(args) -> None:
    """Copy the measured detector constants and re-quote the inflation.

    Reads ``detector_params`` and ``s2_linearity_caps`` (manifest, read
    only).  Nothing is stored for a quantity the detector package did not
    measure.  The error-inflation factor of every series was computed with
    the photon term at the NOMINAL gain (``series.NOMINAL_GAIN_E_PER_ADU``);
    with the measured gain K the photon term changes by sqrt(1.057 / K),
    so the inflation factor of a photon-limited series changes by
    sqrt(K / 1.057).  That is re-quoted per mode as the range it can take
    (read-noise- and floor-limited series do not move at all).
    """
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    rvs.clear_stage(con, "detector")
    man = sqlite3.connect(f"file:{MANIFEST_DB}?mode=ro", uri=True)
    dp = {(r[0], r[1]): (r[2], r[3], r[4]) for r in man.execute(
        "SELECT era_group, quantity, value, uncertainty, method FROM "
        "detector_params")}
    caps = {r[0]: r[1:] for r in man.execute(
        "SELECT mode, cap_adu, worst_dev_pct, slope_pct_per_scale, "
        "slope_err, measured_to_fraction FROM s2_linearity_caps")}
    n_gain = 0
    for mode, tag in DETECTOR_MODES.items():
        g = dp.get((mode, "gain_e_per_adu"))
        if g and g[0] is not None:
            n_gain += 1
            put(con, "detector", f"rv det {tag} gain", g[0], "e/ADU",
                    "f3", f"{mode}: measured gain, detector_params ({g[2]})")
            if g[1] is not None:
                put(con, "detector", f"rv det {tag} gain err", g[1],
                        "e/ADU", "f3", f"{mode}: standard error of that quantity")
            fac = math.sqrt(g[0] / sr.NOMINAL_GAIN_E_PER_ADU)
            infl = [r[0] for r in con.execute(
                "SELECT s.chi2_inflation FROM cv_series s WHERE "
                "s.chi2_inflation IS NOT NULL AND s.series_key IN (SELECT "
                "DISTINCT series_key FROM cv_frames WHERE readoutm=?)",
                (mode,))]
            if infl:
                put(con, "detector", f"rv det {tag} infl min",
                        min(infl), "", "f2", f"{mode}: smallest error-"
                        "inflation factor (nominal-gain photon term)")
                put(con, "detector", f"rv det {tag} infl max",
                        max(infl), "", "f2", f"{mode}: largest one")
                put(con, "detector", f"rv det {tag} infl factor", fac,
                        "", "f3", f"{mode}: sqrt(K / 1.057), the change of "
                        "a photon-limited inflation factor with the "
                        "measured gain")
                put(con, "detector", f"rv det {tag} infl meas min",
                        min(infl) * min(1.0, fac), "", "f2",
                        f"{mode}: smallest inflation factor re-quoted with "
                        "the measured gain (bracket over photon- and "
                        "floor-limited series)")
                put(con, "detector", f"rv det {tag} infl meas max",
                        max(infl) * max(1.0, fac), "", "f2",
                        f"{mode}: largest one, same bracket")
        rn = dp.get((mode, "read_noise_e"))
        if rn and rn[0] is not None and "SUPERSEDED" not in (rn[2] or ""):
            put(con, "detector", f"rv det {tag} rn e", rn[0], "e", "f2",
                    f"{mode}: read noise, detector_params ({rn[2]})")
            if rn[1] is not None:
                put(con, "detector", f"rv det {tag} rn err e", rn[1],
                        "e", "f2", f"{mode}: standard error of that quantity")
        c = caps.get(CAP_ROW.get(mode, mode))
        if c and c[0] is not None:
            put(con, "detector", f"rv det {tag} cap adu", c[0], "ADU",
                    "int", f"{mode}: 1 per cent linearity cap that applies "
                    f"to the CV frames (s2_linearity_caps "
                    f"{CAP_ROW.get(mode, mode)})")
            if c[1] is not None:
                put(con, "detector", f"rv det {tag} cap dev pct",
                        abs(c[1]), "%", "f1", f"{mode}: worst |departure "
                        f"from linear| below the cap")
        cm = caps.get(mode)
        if cm and cm[2] is not None:
            put(con, "detector", f"rv det {tag} lin slope pct", cm[2],
                    "%", "f2", f"{mode}: straight-line slope of the "
                    "departure from linear, per cent per full scale, up to "
                    "the cap")
            put(con, "detector", f"rv det {tag} lin slope err pct",
                    cm[3], "%", "f2", f"{mode}: standard error of that quantity")
    sp = dp.get(("High Gain StackPro", "read_noise_adu"))
    hg = dp.get(("High Gain", "read_noise_adu"))
    if sp and hg and sp[0] and hg[0]:
        put(con, "detector", "rv det stackpro rn ratio",
                (sp[0] / hg[0]) ** 2, "", "f1", "StackPro over High Gain "
                "read-noise VARIANCE: the sum of 16 sub-reads")
    # Constants of method the paper states, recorded where they are set.
    put(con, "detector", "rv lin criterion pct", 1.0, "%", "f0",
            "linearity-cap criterion: |departure from linear| <= 1 per cent "
            "(s2_linearity_caps.note)", origin="constant")
    put(con, "detector", "rv ramp clip sigma", rc.RAMP_CLIP_SIGMA, "",
            "f0", "outlier clip of the ramp fit, robust sigmas "
            "(revision_checks.RAMP_CLIP_SIGMA)", origin="constant")
    put(con, "detector", "rv fit v one span cadences", 3.0, "", "f0",
            "half-span of the magnitude estimator's time grid, in cadences "
            "(run_cv_revision.make_fitter: 3.0 * cadence)",
            origin="constant")
    put(con, "detector", "rv band pair trials", rvs.rv.N_BAND_PAIR_TRIALS,
            "", "int", "trials factor of the band-offset tests "
            "(revision_cv.N_BAND_PAIR_TRIALS)", origin="constant")
    put(con, "detector", "rv det measured gains", n_gain, "", "int",
            "readout modes with a measured gain in detector_params")
    put(con, "detector", "rv det nominal gain", sr.NOMINAL_GAIN_E_PER_ADU,
            "e/ADU", "f3", "the gain the photometry's photon term used",
            origin="constant")
    rvs.stamp(con, "detector")
    con.commit()
    con.close()
    print(f"  {n_gain} measured gains in detector_params")


# ===========================================================================
# STAGE: cap — measurements between the measured linearity cap and the veto
# ===========================================================================
#: Readout mode -> the ``detector_params`` group whose ``linearity_cap_adu``
#: applies, per CV era.  High Gain is split by EGAIN epoch in S2 and every
#: CV High Gain frame (era 7) is in the 1.057 epoch (``s2_camera_configs``),
#: whose cap is 2,950 ADU; the StackPro era 6 is the same epoch.
CAP_GROUP_BY_ERA = {7: "AC4040 High Gain e1.057", 6: "High Gain StackPro",
                    47: "1MHz High Sensitivity 16-bit",
                    72: "1MHz High Sensitivity 16-bit", 76: "Mode0",
                    78: "Fast", 79: "Fast"}


def cmd_cap(args) -> None:
    """What the S2 linearity cap withholds that the saturation veto did not.

    The photometry used to withhold a measurement at 0.92 of each mode's
    clip.  The detector package then measured where each mode departs from
    linearity by 1 per cent, and in four modes that cap lies BELOW the
    veto; ``run_cv_photometry.py recap`` now withholds at the lower of the
    two and the chain downstream was re-run.  This stage counts, per
    series and role, the measurements the cap withholds and the old veto
    did not: peak level at or above the threshold now applied
    (``cv_frames.veto_applied_adu``) but below the old veto mapped into the
    same pixels the same way.
    """
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    con.executescript("""
    DROP TABLE IF EXISTS rv_cap;
    CREATE TABLE rv_cap (series_key TEXT PRIMARY KEY, era_id INTEGER,
        cap_raw_adu REAL, old_veto_raw_adu REAL, n_target INTEGER,
        n_target_capped INTEGER, n_comp INTEGER, n_comp_capped INTEGER,
        n_check INTEGER, n_check_capped INTEGER);
    """)
    rvs.clear_stage(con, "cap")
    tot = dict.fromkeys(("t", "tc", "c", "cc", "k", "kc"), 0)
    stl = {}
    for skey, era, mode in con.execute(
            "SELECT DISTINCT series_key, era_id, readoutm FROM cv_frames "
            "WHERE status='matched' ORDER BY series_key").fetchall():
        old = sr.veto_adu(mode)
        new_raw = sr.photometry_veto_adu(mode, era)
        if old is None or new_raw is None or new_raw >= old:
            continue
        # The old veto in measured-pixel units differs from the applied one
        # by (old - new) raw ADU, divided by the flat median for server-
        # reduced frames (veto_in_reduced_adu is linear with slope 1/F).
        tr = None
        prov = con.execute("SELECT provenance FROM cv_frames WHERE "
                           "series_key=? LIMIT 1", (skey,)).fetchone()[0]
        scale = 1.0
        if prov == "server_reduced":
            import run_cv_photometry as s4
            tr = s4._recon_transform(s4.DEFAULT_RECON_DIR, era)
            scale = 1.0 / tr[0] if tr else 1.0
        rows = con.execute("""SELECT s.role,
            sum(CASE WHEN d.peak + f.bkg_adu >= f.veto_applied_adu AND
                          d.peak + f.bkg_adu < f.veto_applied_adu + ?
                     THEN 1 ELSE 0 END), count(*)
            FROM cv_detections d JOIN cv_frames f ON f.frame_id=d.frame_id
            JOIN cv_stars s ON s.series_key=f.series_key AND
                 s.star_id=d.star_id
            WHERE f.series_key=? AND f.status='matched' AND
                  s.role IN ('target','comp','check')
            GROUP BY s.role""", ((old - new_raw) * scale, skey)).fetchall()
        by = {r[0]: (r[1] or 0, r[2]) for r in rows}
        t, c, k = (by.get(x, (0, 0)) for x in ("target", "comp", "check"))
        con.execute("INSERT INTO rv_cap VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (skey, era, new_raw, old, t[1], t[0], c[1], c[0], k[1],
                     k[0]))
        tot["t"] += t[1]; tot["tc"] += t[0]
        tot["c"] += c[1]; tot["cc"] += c[0]
        tot["k"] += k[1]; tot["kc"] += k[0]
        if skey.startswith("stlmi|"):
            stl[skey] = (t[0], c[0], c[1])
        print(f"  {skey:14s} cap {new_raw:6d} < veto {old:6d}: target "
              f"{t[0]}/{t[1]}, comp {c[0]}/{c[1]}, check {k[0]}/{k[1]}")
    nser = con.execute("SELECT count(*) FROM rv_cap").fetchone()[0]
    put(con, "cap", "rv cap series", nser, "", "int",
            "measured series whose mode's linearity cap is below its veto")
    for key, v, note in (
            ("rv cap target capped", tot["tc"], "target measurements the "
             "cap withholds and the old veto did not"),
            ("rv cap target all", tot["t"], "target detections in those "
             "series"),
            ("rv cap comp capped", tot["cc"], "comparison-star measurements "
             "the cap withholds and the old veto did not"),
            ("rv cap comp all", tot["c"], "comparison-star detections in "
             "those series"),
            ("rv cap check capped", tot["kc"], "check-star measurements the "
             "cap withholds and the old veto did not")):
        put(con, "cap", key, v, "", "int", note)
    if tot["c"]:
        put(con, "cap", "rv cap comp capped pct",
                100.0 * tot["cc"] / tot["c"], "%", "f1",
                "the same as a percentage of comparison detections")
    for skey, (tc, cc, ca) in stl.items():
        f = skey.split("|")[2]
        e = skey.split("|")[1][1:]
        tag = {"7": "seven", "47": "fourseven", "76": "sevensix"}.get(e, e)
        put(con, "cap", f"rv cap stlmi {tag} {f} target", tc, "", "int",
                f"ST LMi {skey}: target points the cap withholds")
        put(con, "cap", f"rv cap stlmi {tag} {f} comp", cc, "", "int",
                f"ST LMi {skey}: comparison measurements the cap withholds")
    put(con, "cap", "rv cap hg adu", sr.linearity_cap_adu("High Gain", 7),
            "ADU", "int", "High Gain (EGAIN 1.057 epoch) 1 per cent "
            "linearity cap, detector_params", origin="constant")
    put(con, "cap", "rv cap hg veto adu", sr.veto_adu("High Gain"),
            "ADU", "int", "High Gain saturation veto the photometry used "
            "before the cap", origin="constant")
    rvs.stamp(con, "cap")
    con.commit()
    con.close()


# ===========================================================================
# STAGE: edgeq — the edge-fit chi-squared per band and era (R7)
# ===========================================================================
def cmd_edgeq(args) -> None:
    """Per (era, band): edges, dof range and reduced chi-squared of the fits.

    Reads ``rv_edge`` (the published magnitude estimator, ``v1``), which
    ``run_cv_revision.py edges`` writes, so the paper can state the edge
    model's fit quality per band and era with its degrees of freedom
    (standing rule 1), not only pooled.
    """
    con = rvs.connect(args.db)
    rvs.ensure_tables(con)
    rvs.clear_stage(con, "edgeq")
    tags = {7: "seven", 47: "fourseven", 76: "sevensix"}
    for est, etag in (("v1", "v one"), ("v2", "v two")):
        rows = con.execute("SELECT era_id, band, chi2nu, dof FROM rv_edge "
                           "WHERE estimator=? AND chi2nu IS NOT NULL",
                           (est,)).fetchall()
        groups: dict = {}
        for era, band, x, dof in rows:
            groups.setdefault((era, band), []).append((x, dof))
        for (era, band), v in sorted(groups.items()):
            x = np.array([a for a, _ in v], dtype=float)
            d = [b for _, b in v]
            k = f"rv edgeq {etag} era {tags.get(era, era)} {band}"
            nt = f"{est} edge fits, era {era}, band {band}"
            put(con, "edgeq", f"{k} n", len(v), "", "int", f"{nt}: count")
            put(con, "edgeq", f"{k} dof min", min(d), "", "int",
                f"{nt}: fewest degrees of freedom")
            put(con, "edgeq", f"{k} dof max", max(d), "", "int",
                f"{nt}: most degrees of freedom")
            put(con, "edgeq", f"{k} chinu med", float(np.median(x)), "",
                "f0" if np.median(x) >= 10 else "f1",
                f"{nt}: median reduced chi-squared")
            put(con, "edgeq", f"{k} chinu min", float(x.min()), "",
                "f0" if x.min() >= 10 else "f1",
                f"{nt}: smallest reduced chi-squared")
            put(con, "edgeq", f"{k} chinu max", float(x.max()), "",
                "f0" if x.max() >= 10 else "f1",
                f"{nt}: largest reduced chi-squared")
    rvs.stamp(con, "edgeq")
    con.commit()
    con.close()
    print("  edge-fit chi-squared per era and band stored")


# ===========================================================================
# STAGE: status / all
# ===========================================================================
def cmd_status(args) -> None:
    con = rvs.connect(args.db, read_only=True)
    for st in ("reduction", "ramp", "mech", "clock", "cap", "detector",
               "edgeq"):
        n = con.execute("SELECT count(*) FROM rv_result WHERE stage=?",
                        (st,)).fetchone()[0]
        print(f"  {st:10s} {n:4d} scalars")
    con.close()


def cmd_all(args) -> None:
    for fn in (cmd_reduction, cmd_ramp, cmd_mech, cmd_clock, cmd_cap,
               cmd_detector, cmd_edgeq):
        print(f"\n=== {fn.__name__[4:]} ===", flush=True)
        fn(args)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn, hlp in (
            ("reduction", cmd_reduction, "one night per era, reduced twice"),
            ("ramp", cmd_ramp, "tie residual against detector x"),
            ("mech", cmd_mech, "mechanical epochs and state steps"),
            ("clock", cmd_clock, "transit clock, printed"),
            ("cap", cmd_cap, "measurements between cap and veto"),
            ("detector", cmd_detector, "measured detector constants"),
            ("edgeq", cmd_edgeq, "edge-fit chi2 per era and band"),
            ("status", cmd_status, "what is stored"),
            ("all", cmd_all, "every stage")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("--db", type=Path, default=DEFAULT_DB)
        p.add_argument("--workers", type=int, default=4,
                       choices=range(1, MAX_WORKERS + 1), metavar="N")
        p.set_defaults(func=fn)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
