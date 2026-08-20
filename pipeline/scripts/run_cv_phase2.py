#!/usr/bin/env python
"""CV-S8 — close out Phase 2 photometry: cloud, colour-extinction, cross-era
discipline, and the faint-phase upper limits.

WHAT THIS SCRIPT DOES, AND WHY EACH PIECE IS HERE
--------------------------------------------------
``run_cv_photometry.py`` measured the light curves and ``run_cv_cattie.py``
tied them to a catalogue.  Four Phase-2 tasks were left open, and each one
is a place where a light curve that LOOKS finished is quietly wrong:

``cloud``       The strategy's cloud cut leaned on the per-image ``ZMAG``
                keyword.  Measured coverage: ZMAG is present on 2,506 CV
                frames and absent or zero on every single Sloan-era polar
                frame (VV Pup 0 of 1,353, EU UMa 0 of 208).  So for the
                polars the PRIMARY cloud channel has to be the comparison
                ensemble's own summed flux.  This stage builds that
                channel, CALIBRATES its threshold against the frames that
                do carry an independent ZMAG, and then tests explicitly
                whether the resulting veto preferentially removes the
                target's faint phases — because a veto that does that is
                not a cleaning step, it is a light-curve editor.

``extinction``  The Honeycutt ensemble absorbs everything common to a frame
                into ``ZP_j``, which removes first-order extinction
                exactly.  What survives is the part that depends on a
                star's colour, ``k'' * colour * airmass``.  This stage fits
                it per (era, filter) with an uncertainty and reports the
                answer either way — a coefficient consistent with zero is
                a RESULT, converted into a bound in the error budget, not
                a term to be forced into the reduction.

``crossera``    ST LMi's G/R/I and g/r/i seasons do not overlap in time, so
                the paper compares two within-era analyses and never
                stitches them.  This stage derives the transformation
                between the two natural systems FROM COMPARISON STARS,
                publishes it as data-release metadata, and then VERIFIES
                the discipline in the products: no series mixes eras, and
                no target magnitude was ever colour-transformed.

``forced``      Polars drop below detection in low states.  Dropping those
                epochs censors exactly the faint half of the distribution
                and biases every duty-cycle statistic.  This stage does
                forced aperture photometry at the solved target position on
                every frame where the target was not detected, emits proper
                3-sigma upper limits, and recomputes the state statistics
                with the limits included (Kaplan-Meier) beside the censored
                versions, so the size of the bias is visible.

All the arithmetic lives in ``macro_phot.phase2`` and is unit-tested in
``pipeline/tests/test_phase2.py``.  This script is I/O, staging, parallelism
and bookkeeping.

STAGES (each resumable, each safe to repeat)
--------------------------------------------
    cloud       ensemble-flux-ratio statistic, ZMAG calibration, veto,
                sculpting test
    extinction  second-order colour-extinction fits per (era, filter)
    crossera    transformation metadata + the discipline assertions
    forced      forced photometry and upper limits (parallel, resumable)
    report      docs/CV_TimeSeries/cv_phase2_completion.html + figures
    status      where every stage stands and what it concludes
    all         cloud -> extinction -> crossera -> forced -> report

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_phase2.py cloud
    $P pipeline/scripts/run_cv_phase2.py extinction
    $P pipeline/scripts/run_cv_phase2.py crossera
    $P pipeline/scripts/run_cv_phase2.py forced --workers 6
    $P pipeline/scripts/run_cv_phase2.py report
    $P pipeline/scripts/run_cv_phase2.py status

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite)
--------------------------------------------------------------
``p2_cloud_frame``    per frame: the ensemble flux ratio, its local normal,
                      the relative statistic, the independent ZMAG
                      transmission where one exists, and the veto flag.
``p2_cloud_series``   per series: how many frames the veto removes.
``p2_cloud_night``    per (series, night): the exemplar clear and cloudy
                      nights the threshold was calibrated on.
``p2_cloud_roc``      the calibration itself — false-veto rate and recall
                      at every candidate threshold.
``p2_cloud_bias``     the sculpting test, per series.
``p2_extinction``     one k'' per (era, filter), with its uncertainty, its
                      significance and its error-budget consequence.
``p2_transform``      the published cross-era transformation coefficients.
``p2_discipline``     the assertions, and whether the products satisfy them.
``p2_limits``         one row per forced measurement: position, method,
                      flux, noise, and the limit or the recovered detection.
``p2_limit_series``   per series: how many epochs gained a limit.
``p2_state_stats``    duty-cycle statistics, censored vs limit-aware.
``p2_meta``           build stamps and every constant this run used.

CONCURRENCY
-----------
An astrometry batch and a catalogue-tie re-run may be writing this archive
at the same time.  Every connection sets ``busy_timeout = 300000``, worker
count is hard-capped at 6, transactions are short and per-chunk, and the
manifest is opened READ-ONLY and never rebuilt.
"""

from __future__ import annotations

import argparse
import math
import sqlite3
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import phase2 as p2                            # noqa: E402
from macro_phot import photometry as ph                        # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
DEFAULT_MANIFEST = (REPO_ROOT / "products" / "manifest"
                    / "rlmt-manifest.sqlite")

#: Stamped into ``p2_meta`` and read by the provenance graph.  Bump it when
#: the arithmetic changes, not when a comment does.
PHASE2_CODE_VERSION = "CV-S8 v1.0 (2026-08-19, Phase-2 completion)"

#: Hard cap on worker processes.  This machine is also running an S1
#: astrometry batch and a catalogue-tie re-run; a stage that saturates the
#: disk queue slows both of them and finishes no sooner itself.
MAX_WORKERS = 6

#: Candidate thresholds the ZMAG calibration sweeps.
THRESHOLD_GRID = [round(0.70 + 0.01 * i, 2) for i in range(31)]

BUSY_TIMEOUT_MS = 300_000


# ===========================================================================
# Database plumbing
# ===========================================================================
def connect(path: Path, read_only: bool = False) -> sqlite3.Connection:
    """One connection, always with the long busy timeout.

    ``busy_timeout`` rather than a lock file: two other workflows are
    writing this archive today, and SQLite's own waiter is the only one that
    all three processes agree about.  Five minutes is longer than any single
    transaction any of the three takes.
    """
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    return con


def git_commit() -> str:
    """Short commit of the tree that produced this run, ``-dirty`` when the
    working tree is not clean.  A product stamped with a clean commit that
    was actually built from edited files is a reproducibility claim that is
    false, so the suffix is not optional."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO_ROOT, capture_output=True, text=True,
                             timeout=20)
        sha = out.stdout.strip() or "unknown"
        dirty = subprocess.run(["git", "status", "--porcelain"],
                               cwd=REPO_ROOT, capture_output=True, text=True,
                               timeout=30).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:                                   # noqa: BLE001
        return "unknown"


def ensure_tables(con: sqlite3.Connection) -> None:
    """Create every table this stage owns.  Idempotent by construction."""
    con.executescript("""
    CREATE TABLE IF NOT EXISTS p2_meta (key TEXT PRIMARY KEY, value TEXT);

    CREATE TABLE IF NOT EXISTS p2_cloud_frame (
        series_key TEXT, frame_id INTEGER, night TEXT, bjd_tdb REAL,
        airmass REAL, exptime REAL,
        n_core INTEGER, flux_ratio REAL, local_median REAL, rel_ratio REAL,
        zmag REAL, zmag_transmission REAL, zmag_label TEXT,
        vetoed INTEGER, target_mag REAL, target_detected INTEGER,
        PRIMARY KEY (series_key, frame_id));

    CREATE TABLE IF NOT EXISTS p2_cloud_series (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, n_frames INTEGER, n_nights INTEGER, n_core_median REAL,
        threshold REAL, n_vetoed INTEGER, frac_vetoed REAL,
        rel_ratio_mad REAL, worst_rel_ratio REAL,
        n_zmag INTEGER, note TEXT);

    CREATE TABLE IF NOT EXISTS p2_cloud_night (
        series_key TEXT, night TEXT, n_frames INTEGER,
        zmag_mad REAL, zmag_span REAL, rel_ratio_mad REAL,
        rel_ratio_min REAL, n_vetoed INTEGER, exemplar TEXT,
        PRIMARY KEY (series_key, night));

    CREATE TABLE IF NOT EXISTS p2_cloud_roc (
        threshold REAL PRIMARY KEY, n_clear INTEGER, n_clear_vetoed INTEGER,
        false_veto_rate REAL, n_attenuated INTEGER,
        n_attenuated_vetoed INTEGER, recall REAL, chosen INTEGER);

    CREATE TABLE IF NOT EXISTS p2_cloud_bias (
        series_key TEXT PRIMARY KEY, n INTEGER, n_vetoed INTEGER,
        median_mag_vetoed REAL, median_mag_kept REAL,
        u REAL, z REAL, p_mannwhitney REAL,
        faint_veto_rate REAL, bright_veto_rate REAL,
        prop_diff REAL, prop_z REAL, p_proportion REAL,
        undetected_veto_rate REAL, detected_veto_rate REAL,
        verdict TEXT);

    CREATE TABLE IF NOT EXISTS p2_extinction (
        era_id INTEGER, filter TEXT, era_label TEXT, n_series INTEGER,
        series_keys TEXT, n_points INTEGER, n_stars INTEGER,
        n_frames INTEGER, colour_label TEXT, colour_ref REAL,
        colour_min REAL, colour_max REAL, airmass_ref REAL,
        airmass_min REAL, airmass_max REAL,
        kpp REAL, kpp_err REAL, t_stat REAL, p_value REAL,
        significant INTEGER, chi2nu REAL,
        rms_before_mmag REAL, rms_after_mmag REAL,
        term_p95_mmag REAL, bound_mmag REAL, verdict TEXT, note TEXT,
        n_clipped INTEGER, kpp_err_formal REAL, kpp_err_boot REAL,
        n_boot INTEGER,
        PRIMARY KEY (era_id, filter));

    CREATE TABLE IF NOT EXISTS p2_transform (
        target_key TEXT, era_from INTEGER, band_from TEXT,
        era_to INTEGER, band_to TEXT, kind TEXT, colour_label TEXT,
        n_stars INTEGER, colour_ref REAL, colour_min REAL, colour_max REAL,
        a REAL, a_err REAL, b REAL, b_err REAL,
        rms_mmag REAL, chi2nu REAL, applied_to_targets INTEGER, note TEXT,
        colour_term_from REAL, colour_term_to REAL,
        b_expected REAL, b_expected_err REAL, b_tension_sigma REAL,
        PRIMARY KEY (target_key, era_from, band_from, era_to, band_to));

    CREATE TABLE IF NOT EXISTS p2_discipline (
        check_id TEXT PRIMARY KEY, statement TEXT, n_checked INTEGER,
        n_violation INTEGER, verdict TEXT, detail TEXT);

    CREATE TABLE IF NOT EXISTS p2_limits (
        series_key TEXT, frame_id INTEGER, night TEXT, bjd_tdb REAL,
        x_px REAL, y_px REAL, pos_method TEXT, n_pos_stars INTEGER,
        pos_rms_px REAL, aper_px REAL, exptime REAL, zp REAL,
        flux_adu REAL, flux_err_adu REAL, sky_adu REAL, sky_rms_adu REAL,
        n_pix REAL, n_sky INTEGER, snr REAL,
        outcome TEXT, limit_flux_adu REAL, limit_mag REAL,
        forced_mag REAL, forced_mag_err REAL,
        confidence TEXT, method TEXT, status TEXT, note TEXT,
        PRIMARY KEY (series_key, frame_id));

    CREATE TABLE IF NOT EXISTS p2_limit_series (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, n_matched INTEGER, n_detected INTEGER,
        n_candidates INTEGER, n_forced INTEGER, n_limits INTEGER,
        n_recovered INTEGER, n_failed INTEGER,
        median_limit_mag REAL, faintest_detection REAL,
        pos_method TEXT, note TEXT);

    CREATE TABLE IF NOT EXISTS p2_state_stats (
        series_key TEXT, statistic TEXT, censored_value REAL,
        with_limits_value REAL, delta REAL, note TEXT,
        PRIMARY KEY (series_key, statistic));
    """)
    # Columns added after the first production run.  Declared here rather
    # than by rebuilding the table, because a rebuild of a table another
    # process may be reading is exactly the operation the concurrency rules
    # for this archive forbid.
    _ensure_columns(con, "p2_extinction",
                    {"n_clipped": "INTEGER", "kpp_err_formal": "REAL",
                     "kpp_err_boot": "REAL", "n_boot": "INTEGER"})
    _ensure_columns(con, "p2_transform",
                    {"colour_term_from": "REAL", "colour_term_to": "REAL",
                     "b_expected": "REAL", "b_expected_err": "REAL",
                     "b_tension_sigma": "REAL"})
    con.commit()


def _ensure_columns(con: sqlite3.Connection, table: str,
                    columns: dict[str, str]) -> None:
    """Add any missing column to an existing table.  Idempotent."""
    have = {r[1] for r in con.execute(f"PRAGMA table_info({table})")}
    for name, decl in columns.items():
        if name not in have:
            con.execute(f"ALTER TABLE {table} ADD COLUMN {name} {decl}")


def meta_write(con: sqlite3.Connection, pairs: dict) -> None:
    con.executemany("INSERT OR REPLACE INTO p2_meta VALUES (?,?)",
                    [(k, str(v)) for k, v in pairs.items()])
    con.commit()


def stamp_constants(con: sqlite3.Connection) -> None:
    """Every tunable this run used, written into the product.

    A constant that lives only in source is a constant a reader of the
    database cannot check.  These are the numbers every verdict on the
    report rests on, so they travel with the numbers."""
    meta_write(con, {
        "phase2_code_version": PHASE2_CODE_VERSION,
        "git_commit": git_commit(),
        "cloud_window_half": p2.CLOUD_WINDOW_HALF,
        "cloud_min_frames": p2.CLOUD_MIN_FRAMES,
        "cloud_core_min_frac": p2.CLOUD_CORE_MIN_FRAC,
        "cloud_min_core": p2.CLOUD_MIN_CORE,
        "cloud_threshold_default": p2.CLOUD_THRESHOLD_DEFAULT,
        "cloud_max_false_veto": p2.CLOUD_MAX_FALSE_VETO,
        "zmag_clear_mag": p2.ZMAG_CLEAR_MAG,
        "zmag_atten_mag": p2.ZMAG_ATTEN_MAG,
        "kpp_significance_t": p2.KPP_SIGNIFICANCE_T,
        "limit_sigma": p2.LIMIT_SIGMA,
        "forced_min_stars": p2.FORCED_MIN_STARS,
        "forced_max_rms_px": p2.FORCED_MAX_RMS_PX,
        "forced_detect_snr": p2.FORCED_DETECT_SNR,
        "aperture_radius_arcsec": ph.APERTURE_RADIUS_ARCSEC,
        "sky_annulus_arcsec": str(ph.SKY_ANNULUS_ARCSEC),
        "max_workers": MAX_WORKERS,
    })


def solved_series(con: sqlite3.Connection) -> list[tuple]:
    """The series this stage acts on: solved, with a real ensemble."""
    return con.execute(
        """SELECT series_key, target_key, era_id, filter, target_star_id,
                  n_comp
           FROM cv_series WHERE status='solved' ORDER BY series_key"""
    ).fetchall()


# ===========================================================================
# Stage: cloud  — the ensemble-flux-ratio veto
# ===========================================================================
def _series_night_groups(con, series_key: str) -> list[tuple[str, list]]:
    """Frames of one series, grouped by night, in time order.

    Grouped by NIGHT because "the local normal" is a property of a night's
    sky.  Running a window across a night boundary would compare tonight's
    first frame with a frame from three weeks ago and call the difference
    cloud."""
    rows = con.execute(
        """SELECT frame_id, night, bjd_tdb, airmass, exptime
           FROM cv_frames
           WHERE series_key=? AND status='matched'
           ORDER BY night, bjd_tdb, frame_id""", (series_key,)).fetchall()
    groups: dict[str, list] = {}
    for r in rows:
        groups.setdefault(r[1] or "", []).append(r)
    return sorted(groups.items())


def cmd_cloud(args) -> None:
    """Build the cloud statistic, calibrate the threshold, apply, test."""
    con = connect(args.db)
    ensure_tables(con)
    stamp_constants(con)
    mcon = connect(args.manifest, read_only=True)
    series = solved_series(con)
    print(f"  {len(series)} solved series")

    # The independent channel, loaded once.  A per-frame query against the
    # manifest would be 8,716 round trips across a network drive for a
    # column that fits in a dict.
    zmag_by_frame = {int(r[0]): (float(r[1]) if r[1] is not None else np.nan)
                     for r in mcon.execute(
                         "SELECT obs_rowid, zmag FROM frames "
                         "WHERE zmag IS NOT NULL")}

    con.execute("DELETE FROM p2_cloud_frame")
    con.execute("DELETE FROM p2_cloud_night")
    con.commit()

    # ---- pass 1: the statistic, per series per night ----------------------
    per_series_meta: dict[str, dict] = {}
    for skey, tkey, era, filt, tstar, n_comp in series:
        groups = _series_night_groups(con, skey)
        if not groups:
            continue
        # Every comparison/check measurement of this series, once.
        lc = con.execute(
            """SELECT frame_id, star_id, inst_mag
               FROM cv_lightcurve
               WHERE series_key=? AND role IN ('comp','check')
                 AND saturated=0 AND inst_mag IS NOT NULL""",
            (skey,)).fetchall()
        by_frame: dict[int, list[tuple[int, float]]] = {}
        for fid, sid, im in lc:
            by_frame.setdefault(int(fid), []).append((int(sid), float(im)))
        # The target's own magnitude per frame — needed ONLY for the
        # sculpting test, never for the veto.
        tgt = {int(f): (m if m is not None else float("nan"))
               for f, m in con.execute(
                   """SELECT frame_id, mag FROM cv_lightcurve
                      WHERE series_key=? AND role='target'""", (skey,))}
        rows_out, n_core_list, night_rows = [], [], []
        for night, frames in groups:
            fids = [int(r[0]) for r in frames]
            n_f = len(fids)
            index = {f: i for i, f in enumerate(fids)}
            star_ids, frame_idx, flux_rate = [], [], []
            for f in fids:
                for sid, im in by_frame.get(f, ()):
                    star_ids.append(sid)
                    frame_idx.append(index[f])
                    # inst_mag = -2.5 log10(flux/exptime) + 25, so the flux
                    # RATE is recovered exactly; using the rate rather than
                    # the raw flux is what makes a night that mixes 60 s and
                    # 240 s exposures comparable frame to frame.
                    flux_rate.append(10.0 ** ((p2.INST_MAG_OFFSET - im) / 2.5))
            if not star_ids:
                # A night on which the ensemble measured nothing at all.
                # Its frames are still WRITTEN, with a null statistic and a
                # zero core, because a frame silently missing from the veto
                # product would look to the next reader like a frame that
                # passed the veto.
                for r in frames:
                    fid = int(r[0])
                    tm = _f(tgt.get(fid))
                    rows_out.append((skey, fid, night, r[2], r[3], r[4],
                                     0, None, None, None,
                                     _f(zmag_by_frame.get(fid)), None, None,
                                     0, tm, 1 if tm is not None else 0))
                night_rows.append((skey, night, n_f, None, None, None, None,
                                   0, ""))
                continue
            core = p2.core_ensemble(star_ids, frame_idx)
            ratio, cnt = p2.ensemble_flux_ratio(star_ids, frame_idx,
                                                flux_rate, n_f, core)
            usable = (len(core) >= p2.CLOUD_MIN_CORE
                      and n_f >= p2.CLOUD_MIN_FRAMES)
            local = (p2.running_median(ratio, p2.CLOUD_WINDOW_HALF) if usable
                     else np.full(n_f, np.nan))
            rel = (p2.relative_ratio(ratio) if usable
                   else np.full(n_f, np.nan))
            # Independent evidence: the header zero point, where it exists.
            zm = np.array([zmag_by_frame.get(f, np.nan) for f in fids])
            zt = p2.zmag_transmission(zm)
            zl = p2.label_by_zmag(zt)
            n_core_list.append(len(core))
            for i, r in enumerate(frames):
                fid = int(r[0])
                tm = _f(tgt.get(fid))
                rows_out.append((
                    skey, fid, night, r[2], r[3], r[4],
                    int(cnt[i]), _f(ratio[i]), _f(local[i]), _f(rel[i]),
                    _f(zm[i]), _f(zt[i]), (str(zl[i]) or None),
                    0, tm, 1 if tm is not None else 0))
            finite_z = zm[np.isfinite(zm) & (zm != 0.0)]
            finite_rel = rel[np.isfinite(rel)]
            night_rows.append((
                skey, night, n_f,
                _f(p2.median_abs_deviation(finite_z)),
                _f(float(finite_z.max() - finite_z.min())
                   if finite_z.size >= 2 else np.nan),
                _f(p2.median_abs_deviation(finite_rel)),
                _f(float(finite_rel.min()) if finite_rel.size else np.nan),
                0, ""))
        con.executemany(
            "INSERT OR REPLACE INTO p2_cloud_frame VALUES "
            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows_out)
        con.executemany(
            "INSERT OR REPLACE INTO p2_cloud_night VALUES (?,?,?,?,?,?,?,?,?)",
            night_rows)
        con.commit()
        per_series_meta[skey] = {
            "target_key": tkey, "era_id": era, "filter": filt,
            "n_nights": len(groups),
            "n_core_median": (float(np.median(n_core_list))
                              if n_core_list else float("nan")),
        }
        print(f"    {skey:16s} {len(rows_out):5d} frames, "
              f"{len(groups):3d} nights, core ensemble "
              f"{per_series_meta[skey]['n_core_median']:.0f}")

    # ---- pass 2: calibrate the threshold on the independent evidence ------
    rel, lab = [], []
    for r, l in con.execute(
            "SELECT rel_ratio, zmag_label FROM p2_cloud_frame "
            "WHERE zmag_label IS NOT NULL AND zmag_label<>''"):
        rel.append(r if r is not None else float("nan"))
        lab.append(l)
    roc = p2.roc_table(rel, lab, THRESHOLD_GRID)
    chosen, reason = p2.choose_threshold(roc)
    if chosen is None:
        chosen = p2.CLOUD_THRESHOLD_DEFAULT
        reason = (reason + f"; falling back to the declared default "
                           f"{chosen:.2f}")
    con.execute("DELETE FROM p2_cloud_roc")
    con.executemany(
        "INSERT OR REPLACE INTO p2_cloud_roc VALUES (?,?,?,?,?,?,?,?)",
        [(r["threshold"], r["n_clear"], r["n_clear_vetoed"],
          _f(r["false_veto_rate"]), r["n_attenuated"],
          r["n_attenuated_vetoed"], _f(r["recall"]),
          1 if abs(r["threshold"] - chosen) < 1e-9 else 0) for r in roc])
    con.commit()
    n_clear = roc[0]["n_clear"] if roc else 0
    n_att = roc[0]["n_attenuated"] if roc else 0
    print(f"\n  calibration set: {n_clear:,} independently-CLEAR and "
          f"{n_att:,} independently-ATTENUATED frames")
    print(f"  threshold chosen: {chosen:.2f}")
    print(f"    {reason}")

    # ---- pass 3: apply, summarise, and test for sculpting -----------------
    con.execute("UPDATE p2_cloud_frame SET vetoed = "
                "CASE WHEN rel_ratio IS NOT NULL AND rel_ratio < ? "
                "THEN 1 ELSE 0 END", (chosen,))
    con.execute("""UPDATE p2_cloud_night SET n_vetoed =
                   (SELECT count(*) FROM p2_cloud_frame f
                    WHERE f.series_key=p2_cloud_night.series_key
                      AND f.night=p2_cloud_night.night AND f.vetoed=1)""")
    con.commit()
    _mark_exemplar_nights(con)

    con.execute("DELETE FROM p2_cloud_series")
    con.execute("DELETE FROM p2_cloud_bias")
    for skey, meta in per_series_meta.items():
        row = con.execute(
            """SELECT count(*), sum(vetoed),
                      min(rel_ratio), count(zmag_transmission)
               FROM p2_cloud_frame WHERE series_key=?""", (skey,)).fetchone()
        rels = [r[0] for r in con.execute(
            "SELECT rel_ratio FROM p2_cloud_frame WHERE series_key=? "
            "AND rel_ratio IS NOT NULL", (skey,))]
        n_f, n_v = int(row[0]), int(row[1] or 0)
        con.execute(
            "INSERT OR REPLACE INTO p2_cloud_series VALUES "
            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (skey, meta["target_key"], meta["era_id"], meta["filter"],
             n_f, meta["n_nights"], meta["n_core_median"], chosen,
             n_v, (n_v / n_f) if n_f else None,
             _f(p2.median_abs_deviation(rels)), _f(row[2]), int(row[3] or 0),
             "" if meta["n_core_median"] >= p2.CLOUD_MIN_CORE
             else f"core ensemble below {p2.CLOUD_MIN_CORE}: no veto applied"))
        # --- the sculpting test ---
        mags, vet = [], []
        for m, v in con.execute(
                "SELECT target_mag, vetoed FROM p2_cloud_frame "
                "WHERE series_key=?", (skey,)):
            mags.append(m if m is not None else float("nan"))
            vet.append(bool(v))
        t = p2.sculpting_test(mags, vet)
        # And the sharper question the magnitudes alone cannot answer: is
        # the veto more likely to fire on frames where the target was NOT
        # detected at all?  Those are the faintest epochs of the series by
        # definition, and they carry no magnitude to put in the U test.
        und = con.execute(
            "SELECT sum(CASE WHEN target_mag IS NULL THEN vetoed ELSE 0 END),"
            "       sum(CASE WHEN target_mag IS NULL THEN 1 ELSE 0 END),"
            "       sum(CASE WHEN target_mag IS NOT NULL THEN vetoed ELSE 0 "
            "END), sum(CASE WHEN target_mag IS NOT NULL THEN 1 ELSE 0 END) "
            "FROM p2_cloud_frame WHERE series_key=?", (skey,)).fetchone()
        uv, un, dv, dn = (int(x or 0) for x in und)
        con.execute(
            "INSERT OR REPLACE INTO p2_cloud_bias VALUES "
            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (skey, t["n"], t["n_vetoed"], _f(t["median_mag_vetoed"]),
             _f(t["median_mag_kept"]), _f(t["u"]), _f(t["z"]),
             _f(t["p_mannwhitney"]), _f(t["faint_veto_rate"]),
             _f(t["bright_veto_rate"]), _f(t["prop_diff"]), _f(t["prop_z"]),
             _f(t["p_proportion"]),
             (uv / un) if un else None, (dv / dn) if dn else None,
             t["verdict"]))
    con.commit()
    meta_write(con, {
        "stage_cloud": datetime.now(timezone.utc).isoformat(),
        "cloud_threshold_used": chosen,
        "cloud_threshold_reason": reason,
        "cloud_n_clear": n_clear, "cloud_n_attenuated": n_att,
    })
    tot = con.execute("SELECT count(*), sum(vetoed) FROM "
                      "p2_cloud_frame").fetchone()
    print(f"\n  vetoed {int(tot[1] or 0):,} of {int(tot[0]):,} frames "
          f"({100 * (tot[1] or 0) / max(1, tot[0]):.2f}%)")
    for skey, verdict in con.execute(
            "SELECT series_key, verdict FROM p2_cloud_bias "
            "WHERE verdict LIKE '%VETO EXCESS'"):
        print(f"    !! {skey}: {verdict}")
    con.close()
    mcon.close()


def _mark_exemplar_nights(con: sqlite3.Connection) -> None:
    """Name the one clearest and the one cloudiest night among the nights
    that carry independent ZMAG evidence.

    The task asks the threshold to be calibrated on "a known-clear night vs
    a known-cloudy one".  The ROC does that over the whole labelled set;
    these two rows make the same argument in a form a reader can check by
    eye, and they are chosen by the INDEPENDENT channel (ZMAG spread) so
    that the exemplars are not selected by the statistic being tested."""
    con.execute("UPDATE p2_cloud_night SET exemplar=''")
    rows = con.execute(
        """SELECT series_key, night, zmag_mad, zmag_span, n_frames
           FROM p2_cloud_night
           WHERE zmag_span IS NOT NULL AND n_frames >= 20
           ORDER BY zmag_span""").fetchall()
    if not rows:
        return
    con.execute("UPDATE p2_cloud_night SET exemplar='clear' "
                "WHERE series_key=? AND night=?", (rows[0][0], rows[0][1]))
    con.execute("UPDATE p2_cloud_night SET exemplar='cloudy' "
                "WHERE series_key=? AND night=?", (rows[-1][0], rows[-1][1]))
    con.commit()


def _f(x):
    """A float SQLite will accept, or NULL.  NaN is stored as NULL on
    purpose: SQLite has no NaN, and a silently-coerced 0.0 in a column of
    magnitudes is the kind of number that later gets averaged."""
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


# ===========================================================================
# Stage: extinction — the second-order colour term
# ===========================================================================
def _primary_band(con, series_key: str) -> tuple[str, str]:
    """The primary catalogue band and colour label tied for one series."""
    row = con.execute(
        "SELECT band, colour_label FROM cv_cattie WHERE series_key=? "
        "AND is_primary=1 LIMIT 1", (series_key,)).fetchone()
    return (row[0], row[1]) if row else (None, None)


def cmd_extinction(args) -> None:
    """Fit ``k'' * colour * airmass`` on the comparison ensembles."""
    con = connect(args.db)
    ensure_tables(con)
    stamp_constants(con)
    series = solved_series(con)
    # Group the work by (era, filter) — the grouping the strategy names,
    # because an era is one camera in one configuration and a filter is one
    # bandpass, and k'' is a property of exactly that pair.
    groups: dict[tuple[int, str], list[str]] = {}
    for skey, tkey, era, filt, _t, _n in series:
        groups.setdefault((int(era), filt), []).append(skey)
    con.execute("DELETE FROM p2_extinction")
    con.commit()
    mcon = connect(args.manifest, read_only=True)
    era_label = {r[0]: f"{r[1]} {r[2]}x{r[3]} ({r[4]}..{r[5]})"
                 for r in mcon.execute(
                     "SELECT era_id, readoutm, naxis1, naxis2, first_night, "
                     "last_night FROM eras")}
    mcon.close()
    print(f"  {len(groups)} (era, filter) groups")
    for (era, filt), keys in sorted(groups.items()):
        resid, sig, colour, airmass, star_i, frame_i = [], [], [], [], [], []
        colour_label = None
        n_airmass_reject = 0
        for skey in keys:
            band, clabel = _primary_band(con, skey)
            if band is None:
                continue
            colour_label = colour_label or clabel
            # Star colours come from the catalogue tie: they are CATALOGUE
            # colours of the comparison stars, not colours measured here,
            # so they are independent of the residuals being fitted.  Using
            # a colour derived from our own photometry would correlate the
            # design column with the noise it is being fitted against.
            # `in_fit=1` only: these are the stars the catalogue tie itself
            # trusted, after its own saturation, blend and ambiguity vetoes.
            # A star the tie threw out is a star whose catalogue colour we
            # do not believe, and its colour is the design column here.
            cols = {int(r[0]): float(r[1]) for r in con.execute(
                "SELECT star_id, colour FROM cv_cattie_star "
                "WHERE series_key=? AND band=? AND in_fit=1 "
                "AND colour IS NOT NULL", (skey, band))}
            means = {int(r[0]): float(r[1]) for r in con.execute(
                "SELECT star_id, mean_mag FROM cv_stars WHERE series_key=? "
                "AND mean_mag IS NOT NULL", (skey,))}
            # Airmass, with the physically impossible cards refused.  See
            # macro_phot.phase2.AIRMASS_MAX for the 6,877 that made this
            # filter necessary.
            airm = {}
            for fid, am in con.execute(
                    "SELECT frame_id, airmass FROM cv_frames "
                    "WHERE series_key=? AND status='matched'", (skey,)):
                if am is None:
                    continue
                if p2.AIRMASS_MIN <= float(am) <= p2.AIRMASS_MAX:
                    airm[int(fid)] = float(am)
                else:
                    n_airmass_reject += 1
            cur = con.execute(
                """SELECT star_id, frame_id, mag, inst_mag_err
                   FROM cv_lightcurve
                   WHERE series_key=? AND role IN ('comp','check')
                     AND saturated=0 AND mag IS NOT NULL""", (skey,))
            while True:
                batch = cur.fetchmany(200_000)
                if not batch:
                    break
                for sid, fid, mag, merr in batch:
                    c = cols.get(int(sid))
                    m0 = means.get(int(sid))
                    x = airm.get(int(fid))
                    if c is None or m0 is None or x is None:
                        continue
                    e = float(merr) if merr and merr > 0 else 0.01
                    resid.append(float(mag) - m0)
                    sig.append(e)
                    colour.append(c)
                    airmass.append(x)
                    # Star index must be unique ACROSS series: star 12 of
                    # ST LMi and star 12 of VV Pup are different stars, and
                    # merging them would let the two-way centring subtract
                    # one field's mean from the other field's data.
                    star_i.append(hash((skey, int(sid))) & 0x7FFFFFFF)
                    frame_i.append(int(fid))
        cols_n = 31
        if len(resid) < 200:
            con.execute(
                "INSERT OR REPLACE INTO p2_extinction VALUES "
                "(" + ",".join("?" * cols_n) + ")",
                (era, filt, era_label.get(era, ""), len(keys),
                 ",".join(keys), len(resid), 0, 0, colour_label,
                 None, None, None, None, None, None,
                 None, None, None, None, 0, None, None, None, None, None,
                 "NOT FITTED",
                 "fewer than 200 comparison measurements carrying both a "
                 "catalogue colour and a usable airmass",
                 0, None, None, 0))
            con.commit()
            print(f"    era {era:3d} {filt:2s}  NOT FITTED "
                  f"({len(resid)} points)")
            continue
        fit = p2.fit_kpp(resid, sig, colour, airmass, star_i, frame_i)
        boot = p2.bootstrap_kpp_error(resid, sig, colour, airmass,
                                      star_i, frame_i,
                                      n_boot=int(args.n_boot))
        # The PUBLISHED uncertainty is the larger of the two.  The formal
        # one assumes every point is independent, which it is not; the
        # bootstrap one resamples the stars, which is the unit that is.
        err_pub = fit.kpp_err
        if math.isfinite(boot) and boot > err_pub:
            err_pub = boot
        t_pub = fit.kpp / err_pub if err_pub > 0 else float("nan")
        p_pub = (math.erfc(abs(t_pub) / math.sqrt(2.0))
                 if math.isfinite(t_pub) else float("nan"))
        significant = math.isfinite(t_pub) and abs(t_pub) >= p2.KPP_SIGNIFICANCE_T
        pub = p2.KppFit(fit.n_points, fit.n_stars, fit.n_frames,
                        fit.n_clipped, fit.kpp, err_pub, t_pub, p_pub,
                        significant, fit.chi2nu, fit.rms_before,
                        fit.rms_after, fit.colour_ref, fit.airmass_ref,
                        fit.span_term_p95)
        c_arr, x_arr = np.asarray(colour), np.asarray(airmass)
        c_span = float(np.percentile(c_arr, 95) - np.percentile(c_arr, 5))
        x_span = float(np.percentile(x_arr, 95) - np.percentile(x_arr, 5))
        bound = p2.bound_mmag(pub, c_span, x_span)
        verdict = "SIGNIFICANT" if significant else "CONSISTENT WITH ZERO"
        note = p2.budget_sentence(pub, c_span, x_span)
        if math.isfinite(boot):
            note += (f"; uncertainty published is the larger of the formal "
                     f"{fit.kpp_err:.5f} and the {args.n_boot}-replicate "
                     f"star bootstrap {boot:.5f}")
        else:
            note += ("; the star bootstrap could not run (too few stars), "
                     "so the formal error stands and is optimistic")
        if n_airmass_reject:
            note += (f"; {n_airmass_reject} frames refused for an AIRMASS "
                     f"card outside [{p2.AIRMASS_MIN:g}, "
                     f"{p2.AIRMASS_MAX:g}]")
        if fit.n_clipped:
            note += (f"; {fit.n_clipped:,} points clipped at "
                     f"{p2.KPP_CLIP_SIGMA:g} sigma")
        con.execute(
            "INSERT OR REPLACE INTO p2_extinction VALUES "
            "(" + ",".join("?" * cols_n) + ")",
            (era, filt, era_label.get(era, ""), len(keys), ",".join(keys),
             fit.n_points, fit.n_stars, fit.n_frames, colour_label,
             _f(fit.colour_ref), _f(np.min(c_arr)), _f(np.max(c_arr)),
             _f(fit.airmass_ref), _f(np.min(x_arr)), _f(np.max(x_arr)),
             _f(fit.kpp), _f(err_pub), _f(t_pub), _f(p_pub),
             int(significant), _f(fit.chi2nu),
             _f(1000 * fit.rms_before), _f(1000 * fit.rms_after),
             _f(1000 * fit.span_term_p95), _f(bound), verdict, note,
             fit.n_clipped, _f(fit.kpp_err), _f(boot), int(args.n_boot)))
        con.commit()
        size = (1000 * fit.span_term_p95) if significant else bound
        print(f"    era {era:3d} {filt:2s}  k''={fit.kpp:+.5f} "
              f"+/- {err_pub:.5f} (formal {fit.kpp_err:.5f}, boot "
              f"{boot:.5f})  t={t_pub:+.2f}  {verdict:20s} "
              f"{'effect' if significant else 'bound'} {size:.2f} mmag  "
              f"(n={fit.n_points:,}, {fit.n_stars} stars, "
              f"{fit.n_frames} frames)")
    meta_write(con, {"stage_extinction":
                     datetime.now(timezone.utc).isoformat()})
    con.close()


# ===========================================================================
# Stage: crossera — transformation metadata + the discipline assertions
# ===========================================================================
#: Which (target, era_from, era_to) pairs to transform, and why each pair is
#: here.  The uppercase/lowercase pairs are the science question; the
#: lowercase/lowercase pair is the CONTROL — two eras in the same nominal
#: system, whose transformation must come out at zero if the method works.
CROSS_ERA_PAIRS = (
    ("stlmi", 7, 76, "Johnson-Cousins-labelled -> Sloan-labelled"),
    ("yzcnc", 7, 72, "Johnson-Cousins-labelled -> Sloan-labelled"),
    ("vvpup", 72, 76, "control: Sloan -> Sloan, different camera era"),
)

#: Filter label pairs.  ``G`` and ``g`` are two LABELS whose glass identity
#: is exactly what section 4 of the catalogue-tie page could not settle;
#: this stage measures the difference rather than assuming it away.
BAND_PAIRS = (("G", "g"), ("R", "r"), ("I", "i"))


def cmd_crossera(args) -> None:
    """Derive the transformation coefficients and verify the discipline."""
    con = connect(args.db)
    ensure_tables(con)
    stamp_constants(con)
    con.execute("DELETE FROM p2_transform")
    con.execute("DELETE FROM p2_discipline")
    con.commit()

    # ---- the coefficients, from comparison stars --------------------------
    for tkey, era_a, era_b, kind in CROSS_ERA_PAIRS:
        for band_a, band_b in BAND_PAIRS:
            fa, fb = band_a, band_b
            if era_a in (72, 76):          # the lowercase control pair
                fa = band_a.lower()
            key_a = f"{tkey}|e{era_a}|{fa}"
            key_b = f"{tkey}|e{era_b}|{fb}"
            rows = con.execute(
                """SELECT ma.cat_row, sa.ens_mag, sa.ens_err, sa.colour,
                          sb.ens_mag, sb.ens_err
                   FROM cv_cattie_star sa
                   JOIN cv_cattie ca ON ca.series_key=sa.series_key
                        AND ca.band=sa.band AND ca.catalogue=sa.catalogue
                   JOIN cv_cat_match ma ON ma.catalogue=sa.catalogue
                        AND ma.target_key=ca.target_key
                        AND ma.era_id=ca.era_id AND ma.star_id=sa.star_id
                   JOIN cv_cattie cb ON cb.series_key=?
                        AND cb.catalogue=sa.catalogue AND cb.is_primary=1
                   JOIN cv_cat_match mb ON mb.catalogue=sa.catalogue
                        AND mb.target_key=cb.target_key
                        AND mb.era_id=cb.era_id AND mb.cat_row=ma.cat_row
                   JOIN cv_cattie_star sb ON sb.series_key=cb.series_key
                        AND sb.band=cb.band AND sb.catalogue=cb.catalogue
                        AND sb.star_id=mb.star_id
                   WHERE sa.series_key=? AND ca.is_primary=1
                     AND sa.role='comp' AND sb.role='comp'
                     AND sa.ens_mag IS NOT NULL AND sb.ens_mag IS NOT NULL
                     AND sa.colour IS NOT NULL""",
                (key_b, key_a)).fetchall()
            zp_a = con.execute(
                "SELECT zp, colour_label, colour_term, colour_err "
                "FROM cv_cattie WHERE series_key=? AND is_primary=1",
                (key_a,)).fetchone()
            zp_b = con.execute(
                "SELECT zp, colour_label, colour_term, colour_err "
                "FROM cv_cattie WHERE series_key=? AND is_primary=1",
                (key_b,)).fetchone()
            n_cols = 24
            if not rows or zp_a is None or zp_b is None:
                con.execute(
                    "INSERT OR REPLACE INTO p2_transform VALUES "
                    "(" + ",".join("?" * n_cols) + ")",
                    (tkey, era_a, fa, era_b, fb, kind, None, 0,
                     None, None, None, None, None, None, None, None, None,
                     0, "no common comparison stars carrying a primary tie "
                        "in both eras",
                     None, None, None, None, None))
                con.commit()
                print(f"    {tkey} {fa}->{fb}: no common tie stars")
                continue
            # Both magnitudes are put on their own era's catalogue zero
            # point first.  What is left is the BANDPASS difference, which
            # is the thing a transformation coefficient is supposed to be.
            m_from = [r[1] - zp_a[0] for r in rows]
            m_to = [r[4] - zp_b[0] for r in rows]
            colour = [r[3] for r in rows]
            sig = [math.sqrt((r[2] or 0.01) ** 2 + (r[5] or 0.01) ** 2)
                   for r in rows]
            fit = p2.fit_transform(m_from, m_to, colour, sig)
            # INDEPENDENT PREDICTION of the slope.  The catalogue tie
            # already measured, for each era separately, how far that era's
            # bandpass sits from the catalogue's — that is its colour term
            # k.  If the two eras' natural systems differ only in bandpass,
            # then the slope measured here star-by-star MUST equal
            # k(to) - k(from), which was measured a completely different
            # way (a regression of ensemble magnitude against catalogue
            # magnitude, not a differencing of two ensembles).  The tension
            # between the two is the sharpest available check that this
            # transformation means what it says.
            b_exp = ((zp_b[2] - zp_a[2])
                     if zp_a[2] is not None and zp_b[2] is not None else None)
            b_exp_err = (math.hypot(zp_a[3] or 0.0, zp_b[3] or 0.0)
                         if b_exp is not None else None)
            tension = None
            if b_exp is not None and math.isfinite(fit.b):
                denom = math.hypot(fit.b_err or 0.0, b_exp_err or 0.0)
                tension = (fit.b - b_exp) / denom if denom > 0 else None
            con.execute(
                "INSERT OR REPLACE INTO p2_transform VALUES "
                "(" + ",".join("?" * n_cols) + ")",
                (tkey, era_a, fa, era_b, fb, kind, zp_a[1], fit.n,
                 _f(fit.x_ref), _f(fit.x_min), _f(fit.x_max),
                 _f(fit.a), _f(fit.a_err), _f(fit.b), _f(fit.b_err),
                 _f(1000 * fit.rms), _f(fit.chi2nu), 0,
                 "published as data-release metadata; NEVER applied to a "
                 "target magnitude",
                 _f(zp_a[2]), _f(zp_b[2]), _f(b_exp), _f(b_exp_err),
                 _f(tension)))
            con.commit()
            t_s = f"{tension:+.1f}" if tension is not None else "  -  "
            print(f"    {tkey} {fa}->{fb}: n={fit.n:4d}  "
                  f"a={fit.a:+.4f}+/-{fit.a_err:.4f}  "
                  f"b={fit.b:+.4f}+/-{fit.b_err:.4f}  "
                  f"rms={1000 * fit.rms:.1f} mmag  "
                  f"(tie predicts b={b_exp:+.4f}, tension {t_s} sigma)"
                  if b_exp is not None else
                  f"    {tkey} {fa}->{fb}: n={fit.n:4d}  "
                  f"a={fit.a:+.4f}+/-{fit.a_err:.4f}  "
                  f"b={fit.b:+.4f}+/-{fit.b_err:.4f}  "
                  f"rms={1000 * fit.rms:.1f} mmag")

    # ---- the discipline assertions ---------------------------------------
    _verify_discipline(con)
    meta_write(con, {"stage_crossera": datetime.now(timezone.utc).isoformat()})
    con.close()


def _verify_discipline(con: sqlite3.Connection) -> None:
    """Four assertions about the PRODUCTS, not about our intentions."""
    checks = []

    # 1. No series mixes eras.
    rows = con.execute(
        """SELECT f.series_key, s.era_id, f.era_id
           FROM cv_frames f JOIN cv_series s USING (series_key)
           WHERE f.era_id IS NOT NULL AND f.era_id <> s.era_id""").fetchall()
    n_frames = con.execute("SELECT count(*) FROM cv_frames").fetchone()[0]
    n_bad, ex = p2.verify_no_era_mixing(rows)
    checks.append((
        "no-era-mixing",
        "Every frame in a series belongs to that series' era.",
        n_frames, n_bad,
        "HOLDS" if n_bad == 0 else "VIOLATED", "; ".join(ex)))

    # 2. No series-key parses to an era its own row disagrees with.
    n_key = 0
    bad_key = []
    for skey, era in con.execute("SELECT series_key, era_id FROM cv_series"):
        n_key += 1
        try:
            parsed = int(skey.split("|")[1][1:])
        except (IndexError, ValueError):
            bad_key.append(f"{skey}: unparseable")
            continue
        if parsed != int(era):
            bad_key.append(f"{skey}: key says e{parsed}, row says e{era}")
    checks.append((
        "series-key-honest",
        "Every series key names the era its own row records.",
        n_key, len(bad_key),
        "HOLDS" if not bad_key else "VIOLATED", "; ".join(bad_key[:20])))

    # 3. No target magnitude was transformed.
    has_cal = con.execute(
        "SELECT count(*) FROM pragma_table_info('cv_lightcurve') "
        "WHERE name='cal_mag'").fetchone()[0]
    if has_cal:
        rows = con.execute(
            """SELECT l.series_key, l.cal_mag, l.mag, c.zp
               FROM cv_lightcurve l
               JOIN cv_cattie c ON c.series_key=l.series_key
                    AND c.is_primary=1
               WHERE l.role='target' AND l.cal_mag IS NOT NULL""").fetchall()
        n_bad, ex = p2.verify_no_target_transform(rows)
        checks.append((
            "no-target-transform",
            "Every calibrated target magnitude equals mag - zp exactly: a "
            "zero-point shift and no colour transformation.",
            len(rows), n_bad, "HOLDS" if n_bad == 0 else "VIOLATED",
            "; ".join(ex)))
    else:
        checks.append((
            "no-target-transform",
            "Every calibrated target magnitude equals mag - zp exactly: a "
            "zero-point shift and no colour transformation.",
            0, 0, "NOT APPLICABLE",
            "cv_lightcurve carries no cal_mag column at the time of this "
            "run — the catalogue tie's `apply` stage has not been run "
            "against this database, so there is no calibrated target "
            "magnitude to check.  The assertion is re-run automatically the "
            "next time this stage executes."))

    # 4. The transformation coefficients are metadata and nothing else.
    n_applied = con.execute(
        "SELECT count(*) FROM p2_transform WHERE applied_to_targets=1"
    ).fetchone()[0]
    n_tr = con.execute("SELECT count(*) FROM p2_transform").fetchone()[0]
    checks.append((
        "transform-is-metadata",
        "The cross-era transformation coefficients are published as "
        "metadata and applied to no light-curve point.",
        n_tr, n_applied, "HOLDS" if n_applied == 0 else "VIOLATED",
        "cv_lightcurve.mag and cal_mag are written only by CV-S4 and "
        "CV-S6; this stage writes no column on cv_lightcurve at all."))

    con.executemany(
        "INSERT OR REPLACE INTO p2_discipline VALUES (?,?,?,?,?,?)", checks)
    con.commit()
    for cid, _s, n, nb, verdict, _d in checks:
        print(f"    {cid:24s} {verdict:15s} "
              f"({nb} violations in {n:,} checked)")


# ===========================================================================
# Stage: forced — forced photometry and upper limits
# ===========================================================================
def _target_ref_xy(con, target_key: str, era_id: int,
                   target_star_id) -> tuple[float, float, str, int, float]:
    """Where the target sits on this block's REFERENCE frame.

    Two routes, in order of trust:

    1.  the target is itself a reference star — the ensemble identified it,
        so its pixel position is a measurement;
    2.  the block has no identified target star (EU UMa era 78, whose field
        tie died on an HTTP error), so a plate model is fitted to the
        reference stars that DO carry sky coordinates and evaluated at the
        target's catalogued position.  The fit residual is returned so the
        row can say how well that position is known.
    """
    if target_star_id is not None:
        r = con.execute(
            "SELECT x, y FROM cv_ref_stars WHERE target_key=? AND era_id=? "
            "AND star_id=?", (target_key, era_id, target_star_id)).fetchone()
        if r:
            return float(r[0]), float(r[1]), "ref_star", 1, 0.0
    tie = con.execute(
        "SELECT target_ra, target_dec FROM cv_field_tie WHERE target_key=? "
        "AND era_id=?", (target_key, era_id)).fetchone()
    if not tie or tie[0] is None:
        return (float("nan"), float("nan"), "none", 0, float("nan"))
    stars = con.execute(
        "SELECT ra_deg, dec_deg, x, y FROM cv_ref_stars WHERE target_key=? "
        "AND era_id=? AND ra_deg IS NOT NULL", (target_key, era_id)).fetchall()
    if len(stars) < 8:
        return (float("nan"), float("nan"), "none", len(stars), float("nan"))
    ra = np.array([s[0] for s in stars], dtype=float)
    dec = np.array([s[1] for s in stars], dtype=float)
    xy = np.array([[s[2], s[3]] for s in stars], dtype=float)
    ra0, dec0 = float(np.median(ra)), float(np.median(dec))
    xi, eta = p2.gnomonic_project(ra, dec, ra0, dec0)
    mat, rms = p2.affine_from_pairs(np.column_stack([xi, eta]), xy)
    if not math.isfinite(rms):
        return (float("nan"), float("nan"), "none", len(stars), float("nan"))
    txi, teta = p2.gnomonic_project([tie[0]], [tie[1]], ra0, dec0)
    v = mat @ np.array([txi[0], teta[0], 1.0])
    return float(v[0]), float(v[1]), "plate_model", len(stars), float(rms)


def _forced_worklist(con, skey: str) -> list[int]:
    """Frames of one series that were matched but produced no magnitude."""
    return [int(r[0]) for r in con.execute(
        """SELECT f.frame_id FROM cv_frames f
           WHERE f.series_key=? AND f.status='matched'
             AND f.frame_id NOT IN (
                 SELECT frame_id FROM cv_lightcurve
                 WHERE series_key=? AND role='target' AND mag IS NOT NULL)
           ORDER BY f.frame_id""", (skey, skey))]


def _measure_one(job: dict) -> dict:
    """Worker: read one frame, place the target, measure, bound.

    Runs in a subprocess.  It receives everything it needs as plain data —
    no database handle crosses the process boundary, because a handle that
    did would serialise six workers behind one lock and turn the
    parallelism into a queue.
    """
    from astropy.io import fits
    out = dict(job)
    out.update({"status": "failed", "note": "", "x": None, "y": None,
                "flux": None, "flux_err": None, "sky": None,
                "sky_rms": None, "n_pix": None, "n_sky": None, "snr": None,
                "n_pos": 0, "pos_rms": None})
    try:
        path = Path(job["archive_root"]) / job["pixel_path"]
        if not path.exists():
            out["note"] = f"pixel file missing: {job['pixel_path']}"
            return out
        with fits.open(path, memmap=False) as hdul:
            hdu = (hdul[1] if len(hdul) > 1 and hdul[1].data is not None
                   else hdul[0])
            data = np.ascontiguousarray(hdu.data, dtype=np.float32)
        if job["provenance"] == "local_master":
            from macro_phot import calib as cb
            dark = (cb.read_master(Path(job["master_dark"]))
                    if job.get("master_dark") else None)
            flat = (cb.read_master(Path(job["master_flat"]))
                    if job.get("master_flat") else None)
            data, _recipe = cb.apply_masters(data, dark, flat)
        # --- put the reference grid onto this frame ---
        src = np.array(job["ref_xy"], dtype=float)
        dst = np.array(job["frame_xy"], dtype=float)
        t = p2.similarity_from_pairs(src, dst)
        out["n_pos"] = t.n
        out["pos_rms"] = t.rms_px
        if t.n < p2.FORCED_MIN_STARS or not math.isfinite(t.rms_px):
            out["note"] = (f"transform from {t.n} stars is not usable "
                           f"(need {p2.FORCED_MIN_STARS})")
            return out
        if t.rms_px > p2.FORCED_MAX_RMS_PX:
            out["note"] = (f"transform closes to {t.rms_px:.2f} px, worse "
                           f"than the {p2.FORCED_MAX_RMS_PX} px bar")
            return out
        x, y = p2.apply_similarity(t, job["ref_x"], job["ref_y"])
        out["x"], out["y"] = x, y
        ny, nx = data.shape
        if not (0 <= x < nx and 0 <= y < ny):
            out["note"] = "target position falls outside this frame"
            return out
        r_ap = float(job["aper_px"])
        scale = r_ap / ph.APERTURE_RADIUS_ARCSEC
        r_in = ph.SKY_ANNULUS_ARCSEC[0] * scale
        r_out = ph.SKY_ANNULUS_ARCSEC[1] * scale
        m = p2.forced_aperture(data, x, y, r_ap, r_in, r_out,
                               job.get("egain"))
        out.update({"flux": m.flux, "flux_err": m.flux_err, "sky": m.sky,
                    "sky_rms": m.sky_rms, "n_pix": m.n_pix,
                    "n_sky": m.n_sky, "snr": m.snr, "status": "ok"})
        return out
    except Exception as exc:                              # noqa: BLE001
        out["note"] = f"{type(exc).__name__}: {exc}"[:300]
        return out


def cmd_forced(args) -> None:
    """Forced photometry on every undetected target epoch, then the limits."""
    con = connect(args.db)
    ensure_tables(con)
    stamp_constants(con)
    mcon = connect(args.manifest, read_only=True)
    archive_root = con.execute(
        "SELECT value FROM cv_build_meta WHERE key='archive_root'"
    ).fetchone()[0]
    workers = max(1, min(int(args.workers), MAX_WORKERS))
    series = solved_series(con)

    done = {(r[0], r[1]) for r in con.execute(
        "SELECT series_key, frame_id FROM p2_limits")} \
        if not args.force else set()
    if args.force:
        con.execute("DELETE FROM p2_limits")
        con.commit()

    jobs: list[dict] = []
    pos_note: dict[str, tuple[str, int, float]] = {}
    for skey, tkey, era, filt, tstar, _n in series:
        rx, ry, method, n_pos, pos_rms = _target_ref_xy(con, tkey, era, tstar)
        pos_note[skey] = (method, n_pos, pos_rms)
        if method == "none":
            print(f"    {skey:16s} SKIPPED — no target position: the field "
                  f"tie for {tkey} era {era} left neither an identified "
                  f"target star nor a usable plate model")
            continue
        # The block's reference-star grid, read ONCE per series rather than
        # once per frame: it is the same few hundred rows every time, and
        # re-reading it per frame turned a 20-minute stage into an hour.
        refs = {int(r[0]): (float(r[1]), float(r[2])) for r in
                con.execute("SELECT star_id, x, y FROM cv_ref_stars "
                            "WHERE target_key=? AND era_id=?", (tkey, era))}
        for fid in _forced_worklist(con, skey):
            if (skey, fid) in done:
                continue
            fr = con.execute(
                """SELECT night, bjd_tdb, pixel_path, provenance,
                          master_dark, master_flat, aper_px, exptime, zp
                   FROM cv_frames WHERE frame_id=?""", (fid,)).fetchone()
            if fr is None or not fr[2] or fr[6] is None:
                continue
            dets = con.execute(
                "SELECT star_id, x, y FROM cv_detections "
                "WHERE frame_id=? AND star_id IS NOT NULL", (fid,)).fetchall()
            if len(dets) < p2.FORCED_MIN_STARS:
                continue
            ref_xy, frame_xy = [], []
            for sid, dx, dy in dets:
                rp = refs.get(int(sid))
                if rp is not None:
                    ref_xy.append(list(rp))
                    frame_xy.append([float(dx), float(dy)])
            if len(ref_xy) < p2.FORCED_MIN_STARS:
                continue
            eg = mcon.execute("SELECT egain FROM frames WHERE obs_rowid=?",
                              (fid,)).fetchone()
            jobs.append({
                "series_key": skey, "frame_id": fid, "night": fr[0],
                "bjd_tdb": fr[1], "pixel_path": fr[2], "provenance": fr[3],
                "master_dark": fr[4], "master_flat": fr[5],
                "aper_px": fr[6], "exptime": fr[7], "zp": fr[8],
                "archive_root": archive_root,
                "ref_x": rx, "ref_y": ry, "pos_method": method,
                "ref_xy": ref_xy, "frame_xy": frame_xy,
                "egain": (float(eg[0]) if eg and eg[0] else None),
            })
    print(f"  {len(jobs):,} frames to force-photometer "
          f"({len(done):,} already done), {workers} workers")
    if not jobs:
        _summarise_limits(con)
        con.close()
        mcon.close()
        return

    t0 = time.time()
    n_ok = n_lim = n_rec = n_fail = 0
    buf: list[tuple] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(_measure_one, j) for j in jobs]
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            row = _limit_row(r)
            buf.append(row)
            if r["status"] != "ok":
                n_fail += 1
            else:
                n_ok += 1
                if row[19] == "limit":
                    n_lim += 1
                elif row[19] == "detection":
                    n_rec += 1
            if len(buf) >= 200:
                con.executemany(
                    "INSERT OR REPLACE INTO p2_limits VALUES "
                    "(" + ",".join("?" * 28) + ")", buf)
                con.commit()
                buf.clear()
                rate = i / max(1e-9, time.time() - t0)
                print(f"    {i:,}/{len(jobs):,}  {rate:.1f} frames/s  "
                      f"limits {n_lim:,}  recovered {n_rec:,}  "
                      f"failed {n_fail:,}", flush=True)
    if buf:
        con.executemany("INSERT OR REPLACE INTO p2_limits VALUES "
                        "(" + ",".join("?" * 28) + ")", buf)
        con.commit()
    print(f"  measured {n_ok:,}, of which {n_lim:,} upper limits and "
          f"{n_rec:,} detections source detection had missed; "
          f"{n_fail:,} could not be measured")
    _summarise_limits(con, pos_note)
    meta_write(con, {"stage_forced": datetime.now(timezone.utc).isoformat()})
    con.close()
    mcon.close()


def _limit_row(r: dict) -> tuple:
    """One ``p2_limits`` row from one worker result.

    The OUTCOME column is the honest three-way split this task needs:
    ``limit`` (measured noise, no source), ``detection`` (forced photometry
    found what source detection missed), ``failed`` (we could not measure
    at all, and say so rather than emitting a limit we cannot defend).
    """
    conf = f"{p2.LIMIT_SIGMA:g} sigma, one-sided Gaussian (99.87%)"
    method = ("forced aperture at the transformed reference position; "
              f"limit = {p2.LIMIT_SIGMA:g} x aperture noise "
              "(sky shot + sky-level uncertainty)")
    if r["status"] != "ok":
        return (r["series_key"], r["frame_id"], r["night"], r["bjd_tdb"],
                _f(r["x"]), _f(r["y"]), r["pos_method"], r["n_pos"],
                _f(r["pos_rms"]), _f(r["aper_px"]), _f(r["exptime"]),
                _f(r["zp"]), None, None, None, None, None, None, None,
                "failed", None, None, None, None, conf, method,
                "failed", r["note"])
    lim_flux = p2.limit_flux(r["flux_err"])
    lim_mag = p2.limit_magnitude(lim_flux, r["exptime"], r["zp"])
    snr = r["snr"] if r["snr"] is not None else float("nan")
    detected = math.isfinite(snr) and snr >= p2.FORCED_DETECT_SNR
    fmag = (p2.limit_magnitude(r["flux"], r["exptime"], r["zp"])
            if detected else float("nan"))
    ferr = ((2.5 / math.log(10.0)) / snr) if detected and snr > 0 else None
    return (r["series_key"], r["frame_id"], r["night"], r["bjd_tdb"],
            _f(r["x"]), _f(r["y"]), r["pos_method"], r["n_pos"],
            _f(r["pos_rms"]), _f(r["aper_px"]), _f(r["exptime"]), _f(r["zp"]),
            _f(r["flux"]), _f(r["flux_err"]), _f(r["sky"]), _f(r["sky_rms"]),
            _f(r["n_pix"]), r["n_sky"], _f(snr),
            "detection" if detected else "limit",
            _f(lim_flux), _f(lim_mag), _f(fmag), _f(ferr),
            conf, method, "ok", "")


def _summarise_limits(con: sqlite3.Connection, pos_note: dict | None = None
                      ) -> None:
    """Per-series census and the duty-cycle statistics, twice over."""
    con.execute("DELETE FROM p2_limit_series")
    con.execute("DELETE FROM p2_state_stats")
    for skey, tkey, era, filt, _t, _n in solved_series(con):
        n_matched = con.execute(
            "SELECT count(*) FROM cv_frames WHERE series_key=? AND "
            "status='matched'", (skey,)).fetchone()[0]
        det = [float(r[0]) for r in con.execute(
            "SELECT mag FROM cv_lightcurve WHERE series_key=? AND "
            "role='target' AND mag IS NOT NULL", (skey,))]
        lims = [float(r[0]) for r in con.execute(
            "SELECT limit_mag FROM p2_limits WHERE series_key=? AND "
            "outcome='limit' AND limit_mag IS NOT NULL", (skey,))]
        rec = [float(r[0]) for r in con.execute(
            "SELECT forced_mag FROM p2_limits WHERE series_key=? AND "
            "outcome='detection' AND forced_mag IS NOT NULL", (skey,))]
        n_forced, n_failed = con.execute(
            "SELECT sum(status='ok'), sum(status='failed') FROM p2_limits "
            "WHERE series_key=?", (skey,)).fetchone()
        # Candidates = matched frames that produced no target magnitude.
        # Counted from the products rather than from the worklist so that a
        # resumed run reports the same number as a fresh one.
        n_cand = con.execute(
            """SELECT count(*) FROM cv_frames f WHERE f.series_key=? AND
               f.status='matched' AND f.frame_id NOT IN (
                 SELECT frame_id FROM cv_lightcurve WHERE series_key=?
                 AND role='target' AND mag IS NOT NULL)""",
            (skey, skey)).fetchone()[0]
        pm = (pos_note or {}).get(skey, ("", 0, float("nan")))
        con.execute(
            "INSERT OR REPLACE INTO p2_limit_series VALUES "
            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            (skey, tkey, era, filt, n_matched, len(det), n_cand,
             int(n_forced or 0), len(lims), len(rec), int(n_failed or 0),
             _f(np.median(lims)) if lims else None,
             _f(max(det)) if det else None, pm[0],
             "" if lims or rec else "no undetected epoch could be bounded"))
        # The two versions of every statistic, side by side.
        s = p2.state_statistics(det + rec, lims)
        stats = [
            ("n_epochs_measured", len(det), len(det) + len(rec) + len(lims)),
            ("detected_fraction", 1.0 if det else float("nan"),
             s["detected_fraction_true"]),
            ("median_mag", s["median_censored"], s["median_km"]),
            ("faint_state_fraction", 0.0, s["faint_state_fraction"]),
            ("faintest_constrained_mag", s["faintest_detection"],
             s["median_limit"]),
        ]
        notes = {
            "n_epochs_measured":
                "epochs the light curve reports, before and after limits",
            "detected_fraction":
                "fraction of measured epochs on which the target was "
                "DETECTED; the censored version is 1 by construction, which "
                "is the bias",
            "median_mag":
                "median magnitude: detections only, vs the Kaplan-Meier "
                "median with the limits treated as right-censored (NaN when "
                "the survival curve never reaches 0.5, i.e. when more than "
                "half the epochs are limits and no median is estimable)",
            "faint_state_fraction":
                "fraction of epochs at which the target was below detection",
            "faintest_constrained_mag":
                "faintest DETECTION vs the median UPPER LIMIT — how much "
                "deeper the limits reach than the detections do",
        }
        for name, a, b in stats:
            con.execute(
                "INSERT OR REPLACE INTO p2_state_stats VALUES (?,?,?,?,?,?)",
                (skey, name, _f(a), _f(b),
                 _f((b - a) if (a is not None and b is not None
                                and math.isfinite(float(a))
                                and math.isfinite(float(b))) else None),
                 notes[name]))
    con.commit()


# ===========================================================================
# Stage: report / status
# ===========================================================================
def cmd_report(args) -> None:
    from macro_phot.report_phase2 import render_report
    out = render_report(args.db)
    print(f"  wrote {out}")


def cmd_status(args) -> None:
    """Where every stage stands, and what it currently concludes."""
    con = connect(args.db, read_only=True)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    print(f"database: {args.db}")
    for t in ("p2_cloud_frame", "p2_cloud_series", "p2_cloud_roc",
              "p2_cloud_bias", "p2_extinction", "p2_transform",
              "p2_discipline", "p2_limits", "p2_limit_series",
              "p2_state_stats"):
        n = con.execute(f"SELECT count(*) FROM {t}").fetchone()[0] \
            if t in have else 0
        print(f"  {t:20s} {n:>10,}")
    if "p2_meta" not in have:
        print("  (nothing run yet)")
        con.close()
        return
    meta = dict(con.execute("SELECT key, value FROM p2_meta"))
    print(f"\n  code version: {meta.get('phase2_code_version', '-')}  "
          f"commit {meta.get('git_commit', '-')}")
    if "p2_cloud_series" in have:
        thr = meta.get("cloud_threshold_used", "-")
        tot = con.execute("SELECT count(*), sum(vetoed) FROM "
                          "p2_cloud_frame").fetchone()
        print(f"\n  CLOUD VETO   threshold {thr}  "
              f"vetoed {int(tot[1] or 0):,}/{int(tot[0] or 0):,}")
        print(f"    {meta.get('cloud_threshold_reason', '')}")
        for r in con.execute(
                "SELECT series_key, n_frames, n_vetoed, frac_vetoed "
                "FROM p2_cloud_series ORDER BY series_key"):
            print(f"    {r[0]:16s} {r[2]:4d}/{r[1]:5d} "
                  f"({100 * (r[3] or 0):5.2f}%)")
        faint = con.execute("SELECT count(*) FROM p2_cloud_bias WHERE "
                            "verdict='FAINT-PHASE VETO EXCESS'").fetchone()[0]
        bright = con.execute(
            "SELECT count(*) FROM p2_cloud_bias WHERE "
            "verdict='BRIGHT-PHASE VETO EXCESS'").fetchone()[0]
        print(f"    sculpting test: {faint} series with a FAINT-side excess "
              f"(the dangerous direction), {bright} with a bright-side one")
    if "p2_extinction" in have:
        print("\n  SECOND-ORDER COLOUR EXTINCTION")
        for r in con.execute(
                "SELECT era_id, filter, kpp, kpp_err, t_stat, verdict, "
                "bound_mmag, n_points FROM p2_extinction "
                "ORDER BY era_id, filter"):
            k = f"{r[2]:+.5f}" if r[2] is not None else "    -    "
            e = f"{r[3]:.5f}" if r[3] is not None else "   -   "
            t = f"{r[4]:+.2f}" if r[4] is not None else "  -  "
            b = f"{r[6]:.2f}" if r[6] is not None else "  -  "
            print(f"    era {r[0]:3d} {r[1]:2s}  k''={k} +/- {e}  t={t}  "
                  f"{r[5]:20s} bound {b} mmag  n={r[7]:,}")
    if "p2_transform" in have:
        print("\n  CROSS-ERA TRANSFORMATIONS (metadata; never applied)")
        for r in con.execute(
                "SELECT target_key, band_from, band_to, n_stars, a, a_err, "
                "b, b_err, rms_mmag FROM p2_transform ORDER BY target_key, "
                "band_from"):
            if not r[3]:
                print(f"    {r[0]:6s} {r[1]}->{r[2]}  (no common stars)")
                continue
            print(f"    {r[0]:6s} {r[1]}->{r[2]}  n={r[3]:4d}  "
                  f"a={r[4]:+.4f}+/-{r[5]:.4f}  b={r[6]:+.4f}+/-{r[7]:.4f}  "
                  f"rms={r[8]:.1f} mmag")
    if "p2_discipline" in have:
        print("\n  DISCIPLINE")
        for r in con.execute("SELECT check_id, verdict, n_violation, "
                             "n_checked FROM p2_discipline"):
            print(f"    {r[0]:24s} {r[1]:16s} {r[2]} / {r[3]:,}")
    if "p2_limit_series" in have:
        print("\n  FAINT-PHASE LIMITS")
        for r in con.execute(
                "SELECT series_key, n_matched, n_detected, n_candidates, "
                "n_limits, n_recovered, n_failed, median_limit_mag "
                "FROM p2_limit_series ORDER BY series_key"):
            if not r[3]:
                continue
            ml = f"{r[7]:.2f}" if r[7] is not None else "  -  "
            print(f"    {r[0]:16s} detected {r[2]:4d}/{r[1]:4d}  "
                  f"undetected {r[3]:4d} -> {r[4]:4d} limits, "
                  f"{r[5]:3d} recovered, {r[6]:3d} failed  "
                  f"median limit {ml}")
    con.close()


# ===========================================================================
# CLI
# ===========================================================================
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--db", type=Path, default=DEFAULT_DB)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--workers", type=int, default=MAX_WORKERS,
                   help=f"worker processes for `forced` (hard cap "
                        f"{MAX_WORKERS}; this machine is also running an S1 "
                        f"batch and a catalogue-tie re-run)")
    p.add_argument("--n-boot", type=int, default=24,
                   help="star-bootstrap replicates for the colour-extinction "
                        "uncertainty (default 24; the published error is the "
                        "larger of formal and bootstrap)")
    p.add_argument("--force", action="store_true",
                   help="redo work already recorded (re-measure every "
                        "forced epoch instead of resuming)")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("cloud", "extinction", "crossera", "forced", "report",
                 "status", "all"):
        sub.add_parser(name)
    args = p.parse_args()
    table = {"cloud": cmd_cloud, "extinction": cmd_extinction,
             "crossera": cmd_crossera, "forced": cmd_forced,
             "report": cmd_report, "status": cmd_status}
    if args.cmd == "all":
        for name in ("cloud", "extinction", "crossera", "forced", "report"):
            print(f"\n=== {name} ===", flush=True)
            table[name](args)
        cmd_status(args)
    else:
        table[args.cmd](args)


if __name__ == "__main__":
    main()
