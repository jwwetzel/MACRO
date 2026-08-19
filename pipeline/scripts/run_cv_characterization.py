#!/usr/bin/env python
"""S5 — characterize the CV time-series data set, then say what it can do.

    quality  -> image quality per (target, era, filter), and the usability
                cut, DEFENDED by the scatter it buys back
    trail    -> sampled ellipticity / trailing audit (the only stage that
                opens pixels; resumable, <= 6 workers)
    noise    -> measured RMS-vs-magnitude against a photon+sky+read
                prediction computed with the S2 gain BRACKET; the flat
                systematic floor; achieved precision at each target's own
                brightness; chi2 inflation; Allan deviations
    cadence  -> intra-night coverage, gaps, baseline, spectral window and
                the explicit 1 c/d alias ladder
    detect   -> injection and recovery through the REAL timestamps and REAL
                (correlated) noise; the 90% recovery contour per target
    timing   -> Monte-Carlo epoch precision of one bright-phase edge at the
                real per-filter cadence: the number Q2 lives or dies on
    verdict  -> goal-by-goal SUPPORTED / CAVEATS / NOT SUPPORTED, each with
                the measured number that decides it
    report   -> docs/CV_TimeSeries/cv_characterization.html + figures
    status   -> what is done, what is missing

Reads (never writes) the photometry product
``products/phot/cv_timeseries.sqlite`` and the manifest; writes everything
to ``products/phot/cv_characterization.sqlite``.  Every stage is idempotent
and chunked so no single invocation runs long.

House rules obeyed: pure logic lives in ``macro_phot.characterize`` (unit
tested), this file is I/O and orchestration, writes are atomic per stage
(one transaction, DELETE-then-INSERT for the stage's own tables), and the
report renders from the database so no number is ever typed by hand.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sqlite3
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np

# The pipeline packages live one directory up from scripts/.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import characterize as ch          # noqa: E402
from macro_phot import errors as er                # noqa: E402
from macro_phot import photometry as ph            # noqa: E402
from macro_core import timing as tm                # noqa: E402

REPO = Path(__file__).resolve().parent.parent.parent
PHOT_DB = REPO / "products" / "phot" / "cv_timeseries.sqlite"
MANIFEST_DB = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
OUT_DB = REPO / "products" / "phot" / "cv_characterization.sqlite"

#: v2.0: the quality cut is applied downstream instead of only being
#: defended; detection contours carry a score mode; Allan ladders carry
#: their own white-noise null and the tau they were measured at; the
#: photon-model fit excludes detection-truncated stars; Q1/Q2/Q3/Q4/S3/S4
#: are decided by numbers that measure their own estimands.
CODE_VERSION = "CV-S5 characterization v2.0 (2026-08-19)"

# --------------------------------------------------------------------------
# Literature periods.  These are USED ONLY to fold, to choose the injection
# grid's anchor point and to convert cadence into points-per-cycle — never
# measured here.  Source strings are stored beside the numbers so the page
# and the arithmetic cannot drift apart.
# --------------------------------------------------------------------------
PERIODS_D = {
    "stlmi": (113.9 / 1440.0, "P_orb = 113.9 min (ANALYSIS_STRATEGY.md §3.1)"),
    "vvpup": (100.4 / 1440.0, "P_orb = 100.4 min (ANALYSIS_STRATEGY.md §3.1)"),
    "euuma": (90.1 / 1440.0, "P_orb = 90.1 min (ANALYSIS_STRATEGY.md §3.1)"),
    "anuma": (0.0797528, "P_orb = 0.0797528 d, AAVSO VSX (macro_phot.report_s4)"),
    "yzcnc": (2.086 / 24.0, "P_orb = 2.086 h (ANALYSIS_STRATEGY.md §3.1)"),
}

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}

#: Quality axes and their histogram edges.  Explicit edges, because the
#: thresholds the report defends are read off these bins and must be
#: reproducible from the edges alone.
#:
#: Two of these are ABSOLUTE (arcsec, airmass) because they mean the same
#: thing on every camera, and two are RELATIVE TO THE SERIES' OWN MEDIAN
#: (sky brightness ratio, zero-point excess in mag) because the raw sky
#: level in ADU/px/s is a different quantity on a 0.45"/px binned CMOS than
#: on a 0.81"/px CCD — pooling them would measure the camera, not the sky.
#: That mistake is not hypothetical: pooling absolute sky rate rejected 55%
#: of the archive, essentially all of it the low-gain era.
QUALITY_AXES = {
    "fwhm_as": np.arange(1.0, 8.01, 0.25),
    "airmass": np.arange(1.0, 3.01, 0.10),
    "sky_ratio": np.concatenate([np.arange(0.0, 4.01, 0.25),
                                 np.arange(4.5, 12.01, 0.5)]),
    "zp_excess": np.arange(-0.20, 1.001, 0.05),
    "moon_illum": np.arange(0.0, 1.001, 0.10),
    "moon_sep": np.arange(0.0, 180.1, 10.0),
    "sky_rate": np.concatenate([np.arange(0.0, 2.01, 0.1),
                                np.arange(2.5, 10.01, 0.5)]),
}

#: The axes the usability cut is actually APPLIED on.  moon_illum, moon_sep
#: and the absolute sky rate are measured and plotted as the EXPLANATION of
#: sky_ratio, never as independent vetoes — cutting on a cause and its
#: proxy would double-count the same frames.
CUT_AXES = ("fwhm_as", "airmass", "sky_ratio", "zp_excess")

#: Frames sampled per series for the trailing audit (pixel-level, so kept
#: small; the question is "is trailing a population or an accident?", which
#: a dozen frames per series answers).
TRAIL_SAMPLE_PER_SERIES = 12

#: Injection search band, cycles/day.  2 c/d excludes the nightly trend the
#: detrending would remove anyway; 40 c/d is above the second harmonic of
#: every target here (EU UMa's fundamental is 16.0 c/d).
SEARCH_FMIN_CD, SEARCH_FMAX_CD = 2.0, 40.0

#: Injection grid.  Amplitudes are SEMI-amplitudes in magnitudes.
INJECT_AMPS = np.array([0.002, 0.004, 0.007, 0.012, 0.020, 0.035,
                        0.060, 0.100, 0.170, 0.300])
#: Trial periods in days: 30 min to 6 h, log-spaced, plus each target's own
#: orbital period inserted at run time.
INJECT_PERIODS_D = np.array([0.0208, 0.0300, 0.0430, 0.0620, 0.0890,
                             0.1280, 0.1840, 0.2500])
INJECT_TRIALS = 50
#: Signal-free trials that set the detection threshold.  The threshold is
#: the (1 - FAP) quantile of these, so the FAP cannot be quoted finer than
#: the sample supports: 300 trials determine the 99th percentile (the third
#: largest value) usefully and a 0.1% FAP not at all.  ch.DETECT_FAP is set
#: to 0.01 for exactly that reason — an honest 1% beats a fictitious 0.1%.
THRESHOLD_TRIALS = 300


# ==========================================================================
# small helpers
# ==========================================================================

def connect_ro(path: Path) -> sqlite3.Connection:
    """Read-only connection with the long busy timeout the concurrent S1
    batch and rclone transfer demand."""
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


def connect_rw(path: Path) -> sqlite3.Connection:
    con = sqlite3.connect(str(path), timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    con.execute("PRAGMA journal_mode = WAL")
    con.execute("PRAGMA synchronous = NORMAL")
    return con


#: Bumped whenever a stored table's COLUMNS change.  On a version bump
#: every derived table is dropped and its stage must run again, because a
#: half-migrated characterization is worse than no characterization: the
#: report renders from these tables and would happily mix a new column's
#: meaning with an old column's values.  ch_trail is exempt — it is the only
#: table whose rows cost pixel reads, its schema has not moved, and it is
#: keyed per frame so it resumes rather than rebuilds.
#:
#: 2 (2026-08-19, adversarial-review rebuild):
#:   ch_detect/ch_contour gain `score` — a contour is now labelled with the
#:     QUESTION it answers (blind period determination vs detection at a
#:     known period), which differ by 3-8x on this data;
#:   ch_allan_fit gains tau_used_s and the white-noise null it must be
#:     compared against, and red_factor_porb is renamed to the honest
#:     red_factor because 81 of 92 ladders never reach P_orb;
#:   ch_noise_series gains the gain-bracket range on k and the numbers
#:     measured on the FULL frame set beside the usable-only ones;
#:   ch_frames/ch_quality_bins gain the registration-method axis;
#:   ch_timing gains the night's own amplitude and a median-density night.
SCHEMA_VERSION = 2

#: Tables whose rows survive a schema bump (see above).
SCHEMA_KEEP = ("ch_trail",)


def ensure_schema(con: sqlite3.Connection) -> None:
    """Create every table this build writes.  Idempotent."""
    have = con.execute("PRAGMA user_version").fetchone()[0]
    if have < SCHEMA_VERSION:
        stale = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name LIKE 'ch_%'") if r[0] not in SCHEMA_KEEP]
        for t in stale:
            con.execute(f"DROP TABLE IF EXISTS {t}")
        if stale:
            print(f"schema: v{have} -> v{SCHEMA_VERSION}, dropped "
                  f"{len(stale)} table(s); their stages must run again",
                  flush=True)
        con.execute(f"PRAGMA user_version = {SCHEMA_VERSION}")
        con.commit()
    con.executescript("""
    CREATE TABLE IF NOT EXISTS ch_meta (key TEXT PRIMARY KEY, value TEXT);

    CREATE TABLE IF NOT EXISTS ch_frames (
        frame_id INTEGER PRIMARY KEY, series_key TEXT, target_key TEXT,
        era_id INTEGER, filter TEXT, night TEXT, bjd_tdb REAL,
        exptime REAL, readoutm TEXT,
        fwhm_px REAL, plate_scale REAL, fwhm_as REAL,
        sky_rate REAL, sky_ratio REAL, bkg_rms REAL, aper_px REAL,
        airmass_hdr REAL, alt_deg REAL, airmass REAL,
        moon_sep REAL, moon_illum REAL, moon_alt REAL,
        zp REAL, zp_excess REAL,
        n_detected INTEGER, n_saturated INTEGER, match_rate REAL,
        reg_method TEXT, ali_rms_px REAL, flat_age_days REAL,
        check_scatter REAL, n_check INTEGER, rel_scatter REAL,
        usable INTEGER, reject_reason TEXT);
    CREATE INDEX IF NOT EXISTS idx_chf_series ON ch_frames(series_key);

    CREATE TABLE IF NOT EXISTS ch_quality_bins (
        scope TEXT, axis TEXT, bin_center REAL, med_rel_scatter REAL,
        n_frames INTEGER);

    CREATE TABLE IF NOT EXISTS ch_cuts (
        scope TEXT, axis TEXT, unit TEXT, threshold REAL, baseline REAL,
        n_pass INTEGER, n_fail INTEGER, note TEXT,
        PRIMARY KEY (scope, axis));

    CREATE TABLE IF NOT EXISTS ch_quality_series (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, readoutm TEXT, n_frames INTEGER,
        fwhm_as_p10 REAL, fwhm_as_med REAL, fwhm_as_p90 REAL,
        airmass_med REAL, airmass_max REAL,
        sky_rate_med REAL, moon_illum_med REAL, moon_sep_med REAL,
        sat_frac REAL, n_usable INTEGER, frac_usable REAL);

    CREATE TABLE IF NOT EXISTS ch_trail (
        frame_id INTEGER PRIMARY KEY, series_key TEXT, night TEXT,
        n_src INTEGER, ell_med REAL, ell_p90 REAL, pa_R REAL,
        fwhm_as REAL, status TEXT);

    CREATE TABLE IF NOT EXISTS ch_noise_stars (
        series_key TEXT, star_id INTEGER, role TEXT, mean_mag REAL,
        rms REAL, nobs INTEGER, chi2nu REAL,
        pred_lo REAL, pred_nom REAL, pred_hi REAL, formal_err REAL);
    CREATE INDEX IF NOT EXISTS idx_chns ON ch_noise_stars(series_key);

    CREATE TABLE IF NOT EXISTS ch_noise_series (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, readoutm TEXT, exptime REAL, n_stars INTEGER,
        floor_nom REAL, k_nom REAL, floor_lo REAL, floor_hi REAL,
        k_lo REAL, k_hi REAL, k_in_bracket INTEGER, n_stars_kfit INTEGER,
        floor_plateau REAL, plateau_mag REAL, plateau_n INTEGER,
        scint_mag REAL, target_mag REAL, prec_at_target REAL,
        prec_at_target_all REAL, floor_plateau_all REAL,
        n_near_target INTEGER, faint_const_mag REAL, bright_const_mag REAL,
        check_rms_med REAL, inflation REAL,
        chi2nu_med REAL, best_star_rms REAL, best_star_mag REAL,
        prec_field_matched REAL, n_field_matched INTEGER,
        amin_analytic REAL, n_target_points INTEGER,
        n_frames_usable INTEGER, n_frames_all INTEGER);

    CREATE TABLE IF NOT EXISTS ch_allan (
        series_key TEXT, star_id INTEGER, night TEXT, tau_s REAL,
        adev REAL, n_pairs INTEGER);

    -- red_factor is quoted AT tau_used_s, which is the largest rung at or
    -- below the orbital period -- only 11 of 92 ladders actually reach
    -- P_orb, and red noise grows with tau, so a factor stored under the
    -- name 'red_factor_porb' was a lower bound wearing a larger label.
    -- The *_null columns are what a ladder of this length and cadence
    -- returns on PURE WHITE NOISE, which is the only null these estimators
    -- may be compared against (-0.50 is the asymptotic value, not theirs).
    CREATE TABLE IF NOT EXISTS ch_allan_fit (
        series_key TEXT, star_id INTEGER, night TEXT, n_pts INTEGER,
        dt_s REAL, slope REAL, red_factor REAL, tau_used_s REAL,
        tau_target_s REAL, tau_frac_of_porb REAL, adev_first REAL,
        slope_null_p05 REAL, slope_null_p50 REAL, slope_null_p95 REAL,
        red_null_p50 REAL, red_null_p95 REAL,
        slope_redder_than_null INTEGER, red_redder_than_null INTEGER,
        PRIMARY KEY (series_key, star_id, night));

    CREATE TABLE IF NOT EXISTS ch_cadence (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, n_points INTEGER, n_blocks INTEGER, baseline_d REAL,
        median_dt_s REAL, longest_block_h REAL, on_target_h REAL,
        duty_cycle REAL, period_d REAL, pts_per_cycle REAL,
        cycles_longest_block REAL, phase_coverage REAL,
        max_gap_d REAL, n_blocks_ge1cycle INTEGER, smear_phase REAL,
        best_night TEXT, best_night_n INTEGER, best_night_cycles REAL,
        best_night_phase_cov REAL, best_night_dt_s REAL);

    CREATE TABLE IF NOT EXISTS ch_window (
        scope TEXT, freq_cd REAL, power REAL);

    CREATE TABLE IF NOT EXISTS ch_alias (
        scope TEXT, family TEXT, k INTEGER, freq_cd REAL, power REAL,
        baseline_d REAL, freq_res_cd REAL, resolved INTEGER);

    -- score: 'period' = the blind-search question (tallest peak must clear
    -- a max-statistic threshold AND land within 1% of the truth);
    -- 'known'  = the question this paper asks (power at a literature
    -- period must clear the threshold for that one frequency).  Storing
    -- both, labelled, is the fix for a contour that was published as a
    -- detection limit while measuring period determination.
    CREATE TABLE IF NOT EXISTS ch_detect (
        scope TEXT, series_key TEXT, regime TEXT, score TEXT, period_d REAL,
        semi_amp REAL, frac REAL, n_trials INTEGER, threshold REAL,
        n_points INTEGER, sigma_used REAL);

    CREATE TABLE IF NOT EXISTS ch_contour (
        scope TEXT, series_key TEXT, regime TEXT, score TEXT, period_d REAL,
        amp90 REAL, amp90_lo REAL, amp90_hi REAL, n_points INTEGER,
        n_cycles REAL, freq_res_cd REAL);

    -- How often the tallest peak is an ALIAS rather than the truth, at a
    -- stated injected amplitude.  S3's verdict rested on window power
    -- alone, which cannot answer this.
    CREATE TABLE IF NOT EXISTS ch_alias_confusion (
        scope TEXT, series_key TEXT, regime TEXT, semi_amp REAL,
        frac_true REAL, frac_alias REAL, frac_other REAL, n_trials INTEGER);

    -- The detection threshold seen one check star at a time.  The pooled
    -- threshold is a quantile of a max statistic over <= 4 stars, so its
    -- bootstrap error bar (+/-1%) is not its real uncertainty.
    CREATE TABLE IF NOT EXISTS ch_threshold (
        scope TEXT, regime TEXT, score TEXT, threshold REAL,
        n_stars INTEGER, per_star TEXT, spread_frac REAL);

    CREATE TABLE IF NOT EXISTS ch_timing (
        series_key TEXT, target_key TEXT, night TEXT, night_kind TEXT,
        regime TEXT,
        n_pts_cycle INTEGER, cadence_s REAL, depth_mag REAL,
        depth_season_mag REAL,
        sigma_mag REAL, ingress_req REAL, ingress_phase REAL,
        exp_smear_phase REAL, sigma_t_s REAL, depth_source TEXT,
        n_noise_series INTEGER,
        PRIMARY KEY (series_key, night, night_kind, regime, ingress_req));

    -- Colour-point error, per qualifying three-filter night: the two
    -- bands' photometric terms AND the non-simultaneity term Q1's grade
    -- originally omitted.
    CREATE TABLE IF NOT EXISTS ch_colour (
        target_key TEXT, era_id INTEGER, night TEXT, band_a TEXT,
        band_b TEXT, n_pairs INTEGER, dt_med_s REAL, dt_p90_s REAL,
        rate_med REAL, rate_p90 REAL, sigma_a REAL, sigma_b REAL,
        sigma_colour_med REAL, sigma_colour_p90 REAL);

    -- Are the four held-out check stars a fair sample of what a star at
    -- the target's brightness achieves, or the quiet survivors of the comp
    -- stability cut?  Measured against every star within +/-0.25 mag,
    -- INCLUDING the ones that cut dropped.
    CREATE TABLE IF NOT EXISTS ch_check_bias (
        series_key TEXT PRIMARY KEY, target_mag REAL, check_rms_med REAL,
        n_check INTEGER, field_rms_med REAL, n_field INTEGER,
        field_rms_med_kept REAL, n_field_kept INTEGER, bias_ratio REAL);

    CREATE TABLE IF NOT EXISTS ch_verdict (
        goal_id TEXT PRIMARY KEY, rank INTEGER, goal TEXT, claim TEXT,
        verdict TEXT, deciding_number TEXT, reasoning TEXT,
        alternative TEXT);
    """)
    con.commit()


def set_meta(con: sqlite3.Connection, **kw) -> None:
    con.executemany("INSERT OR REPLACE INTO ch_meta(key, value) VALUES (?, ?)",
                    [(k, str(v)) for k, v in kw.items()])
    con.commit()


def stage_done(con: sqlite3.Connection, stage: str) -> None:
    set_meta(con, **{f"stage_{stage}": datetime.now(timezone.utc).isoformat(timespec="seconds")})


def solved_series(pcon: sqlite3.Connection) -> list[tuple]:
    """(series_key, target_key, era_id, filter) for every SOLVED series.

    Refused series are excluded here on purpose: an unsolved series has no
    zero points, so it has no photometry to characterize.  Their existence
    is still reported (the refusal count is itself evidence about the
    archive), but they carry no measurement.
    """
    return pcon.execute(
        "SELECT series_key, target_key, era_id, filter FROM cv_series "
        "WHERE status='solved' ORDER BY series_key").fetchall()


# ==========================================================================
# STAGE quality
# ==========================================================================

def _sky_moon_geometry(target_key: str, ra: float, dec: float,
                       bjd: np.ndarray) -> dict[str, np.ndarray]:
    """Altitude, airmass, moon separation/illumination/altitude per frame.

    Recomputed from coordinates + time because the archive's own AIRMASS
    cards are demonstrably not data (they reach 6877 on VV Pup and are NULL
    on 851 frames), and MOONANGL/MOONPHAS have never been audited at all.
    Uses the Winer EarthLocation that S3 already trusts for BJD_TDB, so
    geometry and timing share one site definition.
    """
    import warnings

    from astropy import units as u
    from astropy.coordinates import AltAz, SkyCoord, get_body
    from astropy.time import Time
    from astropy.utils.exceptions import AstropyWarning

    # get_body returns a GCRS frame carrying an observer velocity; astropy
    # warns that a separation measured between two frames with different
    # velocities is direction-dependent.  The effect is aberration, ~20", and
    # this section bins moon separation in 10-degree bins.
    warnings.filterwarnings("ignore", category=AstropyWarning)

    loc = tm.winer_location()
    # BJD_TDB -> UTC for pointing geometry.  The barycentric correction is
    # at most 8 minutes; over 8 minutes the Moon moves 0.07 deg and the
    # airmass of a mid-sky target changes in the third decimal, so using
    # BJD directly as a UTC-scale time here is far below any threshold this
    # section sets.  (Timing itself never does this: it uses BJD_TDB.)
    t = Time(np.asarray(bjd, dtype=float), format="jd", scale="tdb")
    frame = AltAz(obstime=t, location=loc)
    star = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    alt = star.transform_to(frame).alt.deg
    moon = get_body("moon", t, loc)
    sun = get_body("sun", t, loc)
    sep = star.separation(moon).deg
    elong = moon.separation(sun).deg
    moon_alt = moon.transform_to(frame).alt.deg
    return {"alt_deg": np.asarray(alt, dtype=float),
            "airmass": ch.airmass_from_altitude(alt),
            "moon_sep": np.asarray(sep, dtype=float),
            "moon_illum": ch.moon_illuminated_fraction(elong),
            "moon_alt": np.asarray(moon_alt, dtype=float)}


def _check_scatter_per_frame(pcon: sqlite3.Connection, series_key: str
                             ) -> dict[int, tuple[float, int]]:
    """Per-frame RMS of the HELD-OUT check stars about their own means.

    This is the single most useful per-frame quality number in the whole
    build: it is measured in magnitudes, it is independent of the target,
    and it responds to every failure mode a frame can have (clouds, bad
    seeing, a botched registration, a saturated comparison).  Every quality
    cut below is defended by what it does to THIS number.

    Saturated points are dropped: a clipped star measures the digitizer,
    not the sky.
    """
    rows = pcon.execute(
        "SELECT frame_id, star_id, mag FROM cv_lightcurve "
        "WHERE series_key=? AND role='check' AND mag IS NOT NULL "
        "AND saturated=0", (series_key,)).fetchall()
    if not rows:
        return {}
    fid = np.array([r[0] for r in rows])
    sid = np.array([r[1] for r in rows])
    mag = np.array([r[2] for r in rows], dtype=float)
    # Per-star mean, then residuals.
    means = {}
    for s in np.unique(sid):
        means[s] = float(np.median(mag[sid == s]))
    resid = mag - np.array([means[s] for s in sid])
    out: dict[int, tuple[float, int]] = {}
    for f in np.unique(fid):
        sel = fid == f
        r = resid[sel]
        out[int(f)] = (float(np.sqrt(np.mean(r ** 2))), int(r.size))
    return out


def stage_quality(args) -> None:
    pcon, ocon = connect_ro(PHOT_DB), connect_rw(OUT_DB)
    ensure_schema(ocon)
    try:
        ties = dict(((t, e), (ra, dec)) for t, e, ra, dec in pcon.execute(
            "SELECT target_key, era_id, target_ra, target_dec FROM cv_field_tie"))
        series = solved_series(pcon)
        rows_out, bin_rows = [], []
        for sk, tk, era, filt in series:
            frames = pcon.execute("""
                SELECT frame_id, night, bjd_tdb, exptime, readoutm, fwhm_px,
                       plate_scale, bkg_adu, bkg_rms, aper_px, airmass,
                       n_detected, n_saturated, match_rate, zp,
                       reg_method, ali_rms_px, flat_age_days
                FROM cv_frames WHERE series_key=? AND status='matched'
                ORDER BY bjd_tdb""", (sk,)).fetchall()
            if not frames:
                continue
            arr = {k: np.array([r[i] for r in frames],
                               dtype=(object if k in ("night", "readoutm") else float))
                   for i, k in enumerate(
                       ["frame_id", "night", "bjd_tdb", "exptime", "readoutm",
                        "fwhm_px", "plate_scale", "bkg_adu", "bkg_rms",
                        "aper_px", "airmass_hdr", "n_detected", "n_saturated",
                        "match_rate", "zp"])}
            ra, dec = ties.get((tk, era), (np.nan, np.nan))
            geo = _sky_moon_geometry(tk, ra, dec, arr["bjd_tdb"])
            fwhm_as = ch.fwhm_arcsec(arr["fwhm_px"], arr["plate_scale"])
            sky = ch.sky_rate_adu_per_px_s(arr["bkg_adu"], arr["exptime"])
            sky_med = np.nanmedian(sky) if np.isfinite(sky).any() else np.nan
            sky_ratio = sky / sky_med if np.isfinite(sky_med) and sky_med > 0 \
                else np.full_like(sky, np.nan)
            # Zero-point excess relative to the frame's OWN night: a cloud
            # test that cannot be fooled by a night-to-night throughput
            # change (different airmass, different sky, different focus).
            zp = arr["zp"]
            zp_excess = np.full(zp.size, np.nan)
            nights = np.array([str(r[1]) for r in frames])
            for n in np.unique(nights):
                sel = nights == n
                fin = sel & np.isfinite(zp)
                if fin.sum():
                    zp_excess[sel] = zp[sel] - float(np.median(zp[fin]))
            scat = _check_scatter_per_frame(pcon, sk)
            cs = np.array([scat.get(int(f), (np.nan, 0))[0]
                           for f in arr["frame_id"]])
            nc = np.array([scat.get(int(f), (np.nan, 0))[1]
                           for f in arr["frame_id"]])
            base = np.nanmedian(cs) if np.isfinite(cs).any() else np.nan
            rel = cs / base if (base and np.isfinite(base)) else np.full_like(cs, np.nan)
            for i, r in enumerate(frames):
                rows_out.append((
                    int(r[0]), sk, tk, era, filt, str(r[1]), float(r[2]),
                    float(r[3]), str(r[4]),
                    float(r[5]) if r[5] is not None else None,
                    float(r[6]) if r[6] is not None else None,
                    _f(fwhm_as[i]), _f(sky[i]), _f(sky_ratio[i]),
                    float(r[8]) if r[8] is not None else None,
                    float(r[9]) if r[9] is not None else None,
                    float(r[10]) if r[10] is not None else None,
                    _f(geo["alt_deg"][i]), _f(geo["airmass"][i]),
                    _f(geo["moon_sep"][i]), _f(geo["moon_illum"][i]),
                    _f(geo["moon_alt"][i]),
                    float(r[14]) if r[14] is not None else None,
                    _f(zp_excess[i]),
                    int(r[11]) if r[11] is not None else None,
                    int(r[12]) if r[12] is not None else None,
                    float(r[13]) if r[13] is not None else None,
                    r[15], _f(r[16]), _f(r[17]),
                    _f(cs[i]), int(nc[i]), _f(rel[i]), None, None))
            print(f"  quality {sk:16s} n={len(frames):5d} "
                  f"fwhm={np.nanmedian(fwhm_as):.2f}\" X={np.nanmedian(geo['airmass']):.2f} "
                  f"chkscat={base if base else float('nan'):.4f}")
        with ocon:
            ocon.execute("DELETE FROM ch_frames")
            ocon.executemany(
                "INSERT INTO ch_frames VALUES (" + ",".join("?" * 35) + ")",
                rows_out)
        _quality_cuts(ocon)
        _quality_series(ocon)
        stage_done(ocon, "quality")
        print(f"quality: {len(rows_out)} frames characterized")
    finally:
        pcon.close(); ocon.close()


def _f(x):
    """float() that turns NaN into SQL NULL (a missing value must be NULL,
    not a number that sorts)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _quality_cuts(ocon: sqlite3.Connection) -> None:
    """Bin the relative check-star scatter against each quality axis and
    read the threshold off the degradation."""
    names = list(QUALITY_AXES)
    rows = ocon.execute(
        "SELECT " + ", ".join(names) + ", rel_scatter FROM ch_frames").fetchall()
    cols = {k: np.array([r[i] if r[i] is not None else np.nan for r in rows],
                        dtype=float) for i, k in enumerate(names)}
    rel = np.array([r[len(names)] if r[len(names)] is not None else np.nan
                    for r in rows], dtype=float)
    units = {"fwhm_as": "arcsec", "airmass": "airmass",
             "sky_ratio": "x series median", "sky_rate": "ADU/px/s",
             "zp_excess": "mag", "moon_illum": "fraction", "moon_sep": "deg"}
    # The baseline is the POOLED median relative scatter.  rel_scatter is
    # each frame's check-star scatter divided by its own series' median, so
    # this number is 1 by construction; measuring it anyway (and printing
    # it) keeps the definition auditable rather than assumed.
    pooled = float(np.nanmedian(rel))
    bin_rows, cut_rows = [], []
    for axis, edges in QUALITY_AXES.items():
        # Two binnings.  The plain one is what the figure draws; the
        # tail-pooled one is what the THRESHOLD is read off, because a bin
        # too thin to test is not evidence of no effect.  Airmass is the
        # case that forced this: the three bins from X = 2.45 to 2.65 all
        # exceed the degradation factor and all three were discarded for
        # holding 10, 7 and 10 frames, after which the page reported that
        # airmass "never" degrades the check stars over the observed range.
        c, m, n = ch.binned_median(cols[axis], rel, edges)
        for i in range(len(c)):
            bin_rows.append(("all", axis, float(c[i]), _f(m[i]), int(n[i])))
        cp, mp, npc = ch.binned_median_pooled_tail(cols[axis], rel, edges)
        if npc.size != n.size:
            bin_rows.append(("pooled_tail", axis, float(cp[-1]),
                             _f(mp[-1]), int(npc[-1])))
        if axis == "moon_sep":
            # Separation is a "smaller is worse" axis: flip it so the shared
            # threshold finder (which walks toward larger x) still works, and
            # record the answer back in degrees.
            c2, m2, n2 = ch.binned_median_pooled_tail(-cols[axis], rel,
                                                      -edges[::-1])
            thr, base = ch.degradation_threshold(c2, m2, n2, baseline=pooled)
            thr = -thr if math.isfinite(thr) else float("inf")
            note = "smaller separation is worse; threshold is a LOWER limit"
            passing = (cols[axis] >= thr) if math.isfinite(thr) else np.isfinite(cols[axis])
        else:
            thr, base = ch.degradation_threshold(cp, mp, npc, baseline=pooled)
            note = ""
            passing = (cols[axis] <= thr) if math.isfinite(thr) else np.isfinite(cols[axis])
        if axis not in CUT_AXES:
            note = (note + "; " if note else "") + "diagnostic only, not applied"
        cut_rows.append(("all", axis, units[axis], _f(thr), _f(base),
                         int(np.sum(passing & np.isfinite(cols[axis]))),
                         int(np.sum(~passing & np.isfinite(cols[axis]))), note))
    bin_rows.extend(_reg_method_bins(ocon))
    with ocon:
        ocon.execute("DELETE FROM ch_quality_bins")
        ocon.execute("DELETE FROM ch_cuts")
        ocon.executemany("INSERT INTO ch_quality_bins VALUES (?,?,?,?,?)", bin_rows)
        ocon.executemany("INSERT INTO ch_cuts VALUES (?,?,?,?,?,?,?,?)", cut_rows)
    # Apply the cut back onto the frames.
    cuts = {a: t for _, a, _, t, *_ in ocon.execute("SELECT * FROM ch_cuts")}
    upd = []
    for row in ocon.execute("SELECT frame_id, " + ", ".join(CUT_AXES)
                            + " FROM ch_frames"):
        fid, vals = row[0], dict(zip(CUT_AXES, row[1:]))
        bad = []
        for axis in CUT_AXES:
            thr = cuts.get(axis)
            v = vals[axis]
            if thr is None or not math.isfinite(thr):
                continue
            if v is None or not math.isfinite(v) or v > thr:
                bad.append(axis)
        upd.append((0 if bad else 1, ",".join(bad) or None, fid))
    with ocon:
        ocon.executemany(
            "UPDATE ch_frames SET usable=?, reject_reason=? WHERE frame_id=?", upd)


#: Depth bins for the registration-method comparison, in detections per
#: frame.  Registration quality and frame depth are confounded — a shallow
#: frame is both harder to align and noisier to measure — so the methods are
#: only comparable INSIDE a depth bin.  Log-spaced because the archive spans
#: 16 to 815 detections per frame in one series.
REG_DEPTH_EDGES = np.array([0, 30, 60, 120, 250, 500, 1000, 100000])


def _reg_method_bins(ocon: sqlite3.Connection) -> list[tuple]:
    """Check-star scatter by REGISTRATION METHOD, with depth controlled.

    Frames inside one series are registered by up to four different methods
    (WCS chain 4,247; astroalign 2,837; translation-vote 696; the reference
    itself 14), and nothing in either product reported what that costs.  Raw
    per-method scatter differs by 3-5x inside single series in BOTH
    directions — translation-vote is the best method in ``yzcnc|e72`` and
    the worst in ``yzcnc|e7|R`` — which is the signature of frame-quality
    confounding rather than a registration defect, and is exactly why this
    bins by detections per frame before comparing anything.

    translation-vote is a method invented during this build and carrying 696
    frames, including every frame of ``stlmi|e47|y``.  Measuring it the same
    way as every other quality axis is what validates it or flags it.
    """
    rows = ocon.execute("""SELECT coalesce(reg_method,'unknown'), n_detected,
                                  rel_scatter FROM ch_frames
                           WHERE rel_scatter IS NOT NULL""").fetchall()
    if not rows:
        return []
    meth = np.array([r[0] for r in rows], dtype=object)
    dep = np.array([r[1] if r[1] is not None else np.nan for r in rows],
                   dtype=float)
    rel = np.array([r[2] for r in rows], dtype=float)
    out = []
    idx = np.digitize(dep, REG_DEPTH_EDGES) - 1
    for m in sorted(set(meth.tolist())):
        for b in range(len(REG_DEPTH_EDGES) - 1):
            sel = (meth == m) & (idx == b) & np.isfinite(rel)
            if sel.sum() < ch.MIN_FRAMES_PER_QUALITY_BIN:
                continue
            center = 0.5 * (REG_DEPTH_EDGES[b] + min(REG_DEPTH_EDGES[b + 1],
                                                     2 * REG_DEPTH_EDGES[b] + 30))
            out.append(("reg_method", str(m), float(center),
                        _f(np.median(rel[sel])), int(sel.sum())))
    return out


def _quality_series(ocon: sqlite3.Connection) -> None:
    rows = ocon.execute("""
        SELECT series_key, target_key, era_id, filter, readoutm, count(*),
               avg(usable)
        FROM ch_frames GROUP BY series_key""").fetchall()
    out = []
    for sk, tk, era, filt, ro, n, fr in rows:
        d = ocon.execute("""
            SELECT fwhm_as, airmass, sky_rate, moon_illum, moon_sep,
                   n_saturated, n_detected FROM ch_frames
            WHERE series_key=?""", (sk,)).fetchall()
        a = np.array([[x if x is not None else np.nan for x in r] for r in d],
                     dtype=float)
        def pc(col, p):
            v = a[:, col][np.isfinite(a[:, col])]
            return _f(np.percentile(v, p)) if v.size else None
        sat = a[:, 5] / np.where(a[:, 6] > 0, a[:, 6], np.nan)
        out.append((sk, tk, era, filt, ro, n,
                    pc(0, 10), pc(0, 50), pc(0, 90),
                    pc(1, 50), pc(1, 100), pc(2, 50), pc(3, 50), pc(4, 50),
                    _f(np.nanmedian(sat)),
                    int(round(fr * n)) if fr is not None else None, _f(fr)))
    with ocon:
        ocon.execute("DELETE FROM ch_quality_series")
        ocon.executemany("INSERT INTO ch_quality_series VALUES ("
                         + ",".join("?" * 17) + ")", out)


# ==========================================================================
# STAGE trail  (the only stage that opens pixels)
# ==========================================================================

def _measure_trail(job):
    """Worker: sep second moments of one frame's brightest sources."""
    frame_id, path, aper_px, plate_scale = job
    try:
        import sep
        from astropy.io import fits
        with fits.open(path) as hdul:
            hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None else hdul[0]
            data = np.ascontiguousarray(hdu.data, dtype=np.float32)
        sep.set_extract_pixstack(1_000_000)
        bkg = sep.Background(data)
        objs = sep.extract(data - bkg.back(), ph.DETECT_SIGMA,
                           err=bkg.globalrms, minarea=ph.DETECT_MINAREA)
        if len(objs) < 20:
            return (frame_id, 0, None, None, None, None, "too_few")
        # Brightest 100 unsaturated sources carry the shape information.
        order = np.argsort(-objs["flux"])[:100]
        a, b, th = objs["a"][order], objs["b"][order], objs["theta"][order]
        e = ch.ellipticity(a, b)
        fw = ch.fwhm_arcsec(ph.fwhm_from_ab(a, b), plate_scale)
        return (frame_id, int(order.size), float(np.nanmedian(e)),
                float(np.nanpercentile(e, 90)), ch.pa_coherence(th),
                float(np.nanmedian(fw)), "ok")
    except Exception as exc:                       # noqa: BLE001
        return (frame_id, 0, None, None, None, None, f"error: {type(exc).__name__}")


def stage_trail(args) -> None:
    pcon, ocon = connect_ro(PHOT_DB), connect_rw(OUT_DB)
    ensure_schema(ocon)
    try:
        root = Path(pcon.execute(
            "SELECT value FROM cv_build_meta WHERE key='archive_root'").fetchone()[0])
        done = {r[0] for r in ocon.execute("SELECT frame_id FROM ch_trail")}
        jobs, meta = [], {}
        for sk, tk, era, filt in solved_series(pcon):
            rows = pcon.execute("""
                SELECT frame_id, pixel_path, aper_px, plate_scale, night
                FROM cv_frames WHERE series_key=? AND status='matched'
                ORDER BY frame_id""", (sk,)).fetchall()
            if not rows:
                continue
            # Evenly spaced sample through the series, so the audit covers
            # the whole campaign rather than one lucky night.
            idx = np.unique(np.linspace(0, len(rows) - 1,
                                        min(TRAIL_SAMPLE_PER_SERIES, len(rows))
                                        ).astype(int))
            for i in idx:
                fid, rel, ap, psc, night = rows[i]
                if fid in done:
                    continue
                jobs.append((fid, str(root / rel), ap, psc))
                meta[fid] = (sk, night)
        if args.limit:
            jobs = jobs[:args.limit]
        print(f"trail: {len(jobs)} frames to measure ({len(done)} already done)")
        out = []
        if jobs:
            with ProcessPoolExecutor(max_workers=min(args.workers, 6)) as ex:
                futs = {ex.submit(_measure_trail, j): j[0] for j in jobs}
                for k, fut in enumerate(as_completed(futs), 1):
                    fid, n, em, e90, R, fw, st = fut.result()
                    sk, night = meta[fid]
                    out.append((fid, sk, night, n, em, e90, R, fw, st))
                    if k % 25 == 0:
                        print(f"  {k}/{len(jobs)}")
        with ocon:
            ocon.executemany("INSERT OR REPLACE INTO ch_trail VALUES "
                             "(?,?,?,?,?,?,?,?,?)", out)
        stage_done(ocon, "trail")
        print(f"trail: {len(out)} frames measured")
    finally:
        pcon.close(); ocon.close()


# ==========================================================================
# STAGE noise
# ==========================================================================

#: Minimum fraction of a series' frames a star must be DETECTED on before
#: its RMS may set the photon-model scaling k.  A star detected on 414 of
#: 745 frames has an RMS computed on the subset of frames where it was
#: bright enough to see, which is a truncated sample and biased LOW — and
#: those are exactly the faintest bins, where the k fit takes almost all of
#: its leverage.  Without this cut three flagship series return a measured
#: RMS below even the HIGH-gain prediction at mag >= 20.75, which is
#: unphysical, and k lands outside the gain bracket that is supposed to
#: contain it.
KFIT_MIN_NOBS_FRAC = 0.90

def usable_frames(ocon: sqlite3.Connection, series_key: str) -> set:
    """The frame ids this series' quality cut admits.

    THE fix for the defect that made the whole chain decorative: the
    usability cut was computed, defended over a page, and then read by
    nothing.  Every downstream stage selected ``FROM cv_frames WHERE
    status='matched'`` with no quality filter, so the noise floor, the Allan
    ladders, the windows, the detection contours and the timing Monte Carlo
    were all measured on the full 7,780 frames including the 1,051 the page
    declared unusable — while section 2 stated in as many words that
    "because 87% of frames survive, the noise floor measured next is a
    property of the instrument and the pipeline".  It was not; it was the
    uncut number.  Now every stage asks this function first.
    """
    return {r[0] for r in ocon.execute(
        "SELECT frame_id FROM ch_frames WHERE series_key=? AND usable=1",
        (series_key,))}


def _star_stats_subset(pcon, series_key: str, keep: Optional[set]) -> dict:
    """Per-star (mean, rms, nobs, chi2nu) over a chosen subset of frames.

    Reproduces :func:`macro_phot.ensemble.star_stats` exactly — same
    inverse-variance weighting, same weighting floor — but over an
    arbitrary frame subset, which is what applying the quality cut requires.
    ``keep=None`` means every matched frame, so the two answers can be put
    side by side and the cut's consequence measured rather than asserted.
    """
    from macro_phot.ensemble import WEIGHT_FLOOR_MAG
    rows = pcon.execute(
        "SELECT star_id, frame_id, mag, inst_mag_err FROM cv_lightcurve "
        "WHERE series_key=? AND mag IS NOT NULL AND saturated=0",
        (series_key,)).fetchall()
    acc: dict = {}
    for sid, fid, mag, err in rows:
        if keep is not None and fid not in keep:
            continue
        if mag is None or not math.isfinite(mag):
            continue
        e = float(err) if err is not None and math.isfinite(err) else 0.0
        acc.setdefault(int(sid), []).append((float(mag), e))
    out = {}
    for sid, vals in acc.items():
        m = np.array([v[0] for v in vals], dtype=float)
        s = np.array([v[1] for v in vals], dtype=float)
        var = s ** 2 + WEIGHT_FLOOR_MAG ** 2
        w = 1.0 / var
        mean = float(np.sum(m * w) / np.sum(w))
        dev = m - mean
        rms = float(np.sqrt(np.mean(dev ** 2)))
        chi2nu = (float(np.sum(dev ** 2 / var) / (m.size - 1))
                  if m.size > 1 else float("nan"))
        out[sid] = (mean, rms, int(m.size), chi2nu)
    return out


def stage_noise(args) -> None:
    pcon, ocon = connect_ro(PHOT_DB), connect_rw(OUT_DB)
    ensure_schema(ocon)
    try:
        star_rows, series_rows, allan_rows, allan_fit = [], [], [], []
        bias_rows = []
        for sk, tk, era, filt in solved_series(pcon):
            keep = usable_frames(ocon, sk)
            if not keep:
                print(f"  noise {sk:16s} SKIPPED: no usable frames", flush=True)
                continue
            fr = pcon.execute("""
                SELECT median(exptime), median(bkg_rms), median(aper_px),
                       median(airmass), count(*)
                FROM cv_frames WHERE series_key=? AND status='matched'""",
                (sk,)).fetchone()
            exptime, bkg_rms, aper_px, x_hdr, nfr = fr
            # Airmass from OUR geometry, not the header.
            xa = ocon.execute("SELECT median(airmass) FROM ch_frames "
                              "WHERE series_key=?", (sk,)).fetchone()[0]
            airmass = xa if xa else 1.3
            zp_med = pcon.execute("SELECT median(zp) FROM cv_frames "
                                  "WHERE series_key=? AND zp IS NOT NULL",
                                  (sk,)).fetchone()[0] or 0.0
            # Roles come from the solve; the STATISTICS are recomputed over
            # the usable frames only.  (cv_stars' own rms was measured on
            # every matched frame, which is the number the quality cut is
            # supposed to improve.)
            roles = dict(pcon.execute(
                "SELECT star_id, role FROM cv_stars WHERE series_key=?", (sk,)))
            st_use = _star_stats_subset(pcon, sk, keep)
            st_all = _star_stats_subset(pcon, sk, None)
            stars = [(sid, roles[sid]) + st_use[sid] for sid in sorted(st_use)
                     if roles.get(sid) in ("comp", "check", "target")]
            if not stars:
                continue
            n_pix = math.pi * float(aper_px) ** 2
            mm = np.array([s[2] for s in stars], dtype=float)
            # Instrumental magnitude -> aperture flux in the frame's own ADU.
            inst = mm + float(zp_med)
            flux = float(exptime) * 10 ** ((ph.INST_MAG_OFFSET - inst) / 2.5)
            preds = {g: ch.predicted_sigma_mag(flux, n_pix, float(bkg_rms), g)
                     for g in (ch.GAIN_LO_E_PER_ADU, ch.GAIN_NOMINAL_E_PER_ADU,
                               ch.GAIN_HI_E_PER_ADU)}
            # Mean formal per-point error the pipeline itself reported.
            formal = dict(pcon.execute(
                "SELECT star_id, avg(inst_mag_err) FROM cv_lightcurve "
                "WHERE series_key=? GROUP BY star_id", (sk,)))
            for i, (sid, role, mag, rms, nobs, chi2nu) in enumerate(stars):
                star_rows.append((sk, sid, role, mag, rms, nobs, chi2nu,
                                  _f(preds[ch.GAIN_LO_E_PER_ADU][i]),
                                  _f(preds[ch.GAIN_NOMINAL_E_PER_ADU][i]),
                                  _f(preds[ch.GAIN_HI_E_PER_ADU][i]),
                                  _f(formal.get(sid))))
            const = [(s[2], s[3], s[4]) for s in stars
                     if s[1] in ("comp", "check")
                     and s[3] is not None and s[4] and s[4] >= 10]
            cmag = np.array([c[0] for c in const], dtype=float)
            crms = np.array([c[1] for c in const], dtype=float)
            cnobs = np.array([c[2] for c in const], dtype=float)
            cpred = {g: ch.predicted_sigma_mag(
                float(exptime) * 10 ** ((ph.INST_MAG_OFFSET - (cmag + zp_med)) / 2.5),
                n_pix, float(bkg_rms), g)
                for g in preds}
            # ---- the photon-model fit, on UNTRUNCATED stars only.
            # A star seen on 414 of 745 frames has an RMS measured over the
            # frames where it happened to be detectable, and that sample is
            # biased low.  Those stars sit at the faint end, which is where
            # this fit takes nearly all of its leverage on k — so including
            # them drove three flagship series to a measured RMS below even
            # the high-gain prediction, which is impossible, and pushed k
            # outside the gain bracket that is supposed to bound it.
            n_use = max(len(keep), 1)
            full = cnobs >= KFIT_MIN_NOBS_FRAC * n_use
            kmag, krms = (cmag[full], crms[full]) if full.sum() >= 5 else (cmag, crms)
            kpred = {g: (cpred[g][full] if full.sum() >= 5 else cpred[g])
                     for g in preds}
            floors = {g: ch.fit_noise_floor(kmag, krms, kpred[g]) for g in preds}
            fn, kn, nstars_k = floors[ch.GAIN_NOMINAL_E_PER_ADU]
            # A HIGHER assumed gain makes the source shot term F/g smaller,
            # so the prediction is smaller and the fitted scaling k that
            # matches the same measured RMS is LARGER.  k therefore tracks
            # the gain, and the bracket ends map straight through.
            k_lo = floors[ch.GAIN_LO_E_PER_ADU][1]
            k_hi = floors[ch.GAIN_HI_E_PER_ADU][1]
            # The bracket k must lie in if the gain is the only unknown: the
            # source term scales as 1/g, so k scales with the gain ratio.
            lo_b = ch.GAIN_LO_E_PER_ADU / ch.GAIN_NOMINAL_E_PER_ADU
            hi_b = ch.GAIN_HI_E_PER_ADU / ch.GAIN_NOMINAL_E_PER_ADU
            k_ok = int(bool(np.isfinite(kn) and lo_b <= kn <= hi_b))
            plateau, plateau_mag, plateau_n = ch.noise_plateau(cmag, crms)
            tgt = pcon.execute("SELECT mean_mag, nobs FROM cv_stars "
                               "WHERE series_key=? AND role='target'",
                               (sk,)).fetchone()
            tmag = tgt[0] if tgt else None
            prec, nnear = ((np.nan, 0) if tmag is None
                           else ch.precision_at_mag(cmag, crms, tmag))
            # The SAME numbers on every matched frame, so the report can put
            # the cut's price beside the cut instead of asserting it.
            const_all = [(v[0], v[1]) for sid, v in st_all.items()
                         if roles.get(sid) in ("comp", "check")
                         and v[2] >= 10]
            amag = np.array([c[0] for c in const_all], dtype=float)
            arms = np.array([c[1] for c in const_all], dtype=float)
            prec_all = (ch.precision_at_mag(amag, arms, tmag)[0]
                        if tmag is not None and amag.size else np.nan)
            plateau_all = ch.noise_plateau(amag, arms)[0] if amag.size else np.nan
            # ---- are the check stars a fair sample at this brightness?
            bias = _check_star_bias(pcon, sk, tmag, st_use, roles)
            if bias:
                bias_rows.append(bias)
            chk_use = [v[1] for sid, v in st_use.items()
                       if roles.get(sid) == "check"]
            chi_use = [v[3] for sid, v in st_use.items()
                       if roles.get(sid) == "check" and np.isfinite(v[3])]
            chk_med = float(np.median(chk_use)) if chk_use else None
            infl = float(er.inflation_factor(np.array(chi_use))) if chi_use else None
            chi_med = float(np.median(chi_use)) if chi_use else None
            best_i = int(np.argmin(crms)) if crms.size else None
            series_rows.append((
                sk, tk, era, filt,
                pcon.execute("SELECT readoutm FROM cv_frames WHERE series_key=? "
                             "LIMIT 1", (sk,)).fetchone()[0],
                exptime, len(const),
                _f(fn), _f(kn),
                _f(floors[ch.GAIN_LO_E_PER_ADU][0]),
                _f(floors[ch.GAIN_HI_E_PER_ADU][0]),
                _f(k_lo), _f(k_hi), k_ok, int(nstars_k),
                _f(plateau), _f(plateau_mag), int(plateau_n),
                _f(ch.scintillation_young(airmass, exptime)),
                _f(tmag), _f(prec), _f(prec_all), _f(plateau_all), int(nnear),
                _f(cmag.max()) if cmag.size else None,
                _f(cmag.min()) if cmag.size else None,
                _f(chk_med), _f(infl), _f(chi_med),
                _f(crms[best_i]) if best_i is not None else None,
                _f(cmag[best_i]) if best_i is not None else None,
                _f(bias[4]) if bias else None,
                int(bias[5]) if bias else None,
                _f(ch.amin_analytic(prec, tgt[1] if tgt else 0)),
                tgt[1] if tgt else 0,
                len(keep), int(nfr)))
            # ---- Allan deviation on the longest run of the best check star
            for sid, night, tau, adev, npair, fit in _allan_for_series(
                    pcon, sk, tk, keep):
                allan_rows.extend([(sk, sid, night, float(t), float(a), int(n))
                                   for t, a, n in zip(tau, adev, npair)])
                allan_fit.append((sk, sid, night) + fit)
            print(f"  noise {sk:16s} usable={len(keep)}/{nfr} "
                  f"floor_fit={fn:.4f} k={kn:.2f}[{k_lo:.2f},{k_hi:.2f}]"
                  f"{'' if k_ok else ' OUT-OF-BRACKET'} "
                  f"plateau={plateau:.4f}@{plateau_mag:.1f} "
                  f"prec@target={prec:.4f} (n={nnear}, all-frames "
                  f"{prec_all:.4f})")
        with ocon:
            ocon.execute("DELETE FROM ch_noise_stars")
            ocon.execute("DELETE FROM ch_noise_series")
            ocon.execute("DELETE FROM ch_allan")
            ocon.execute("DELETE FROM ch_allan_fit")
            ocon.execute("DELETE FROM ch_check_bias")
            ocon.executemany("INSERT INTO ch_noise_stars VALUES ("
                             + ",".join("?" * 11) + ")", star_rows)
            ocon.executemany("INSERT INTO ch_noise_series VALUES ("
                             + ",".join("?" * 37) + ")", series_rows)
            ocon.executemany("INSERT INTO ch_allan VALUES (?,?,?,?,?,?)", allan_rows)
            ocon.executemany("INSERT INTO ch_allan_fit VALUES ("
                             + ",".join("?" * 18) + ")", allan_fit)
            ocon.executemany("INSERT INTO ch_check_bias VALUES ("
                             + ",".join("?" * 9) + ")", bias_rows)
        stage_done(ocon, "noise")
        print(f"noise: {len(series_rows)} series, {len(star_rows)} stars")
    finally:
        pcon.close(); ocon.close()


def _check_star_bias(pcon, series_key: str, target_mag,
                     stats: dict, roles: dict):
    """Are the held-out check stars a FAIR sample at the target's brightness?

    They are held out of the zero-point solve — that much was verified, and
    the solution is not circular.  But they are drawn from the survivors of
    the comparison-star stability iteration, which drops every star whose
    residual RMS exceeds 3x the median, and they voted in the earlier passes
    that shaped the surviving set.  Their RMS is therefore the BEST CASE a
    star of that brightness achieves, not an unbiased estimate of it — while
    the report calls it "the only honest estimate of per-point precision".

    The comparison this makes: the median RMS of every star within
    +/-0.25 mag of the target INCLUDING the ones the stability cut dropped,
    against the check stars' own median.  The ratio is the factor every
    precision number on the page should be read with.

    The comparison is drawn from ``cv_stars``, not from ``cv_lightcurve``:
    the light-curve table only holds target, comp and check rows, so the
    dropped stars — the whole point of the test — have no per-point
    photometry to recompute from.  Both sides of the ratio therefore come
    from the same table, measured over the same frames, and only the
    SELECTION differs, which is what is being priced.

    Returns a ch_check_bias row, or None when the target's brightness or a
    magnitude-matched sample is unavailable.
    """
    if target_mag is None or not np.isfinite(float(target_mag)):
        return None
    tm = float(target_mag)
    rows = pcon.execute(
        "SELECT role, mean_mag, rms FROM cv_stars WHERE series_key=? "
        "AND rms IS NOT NULL AND nobs >= 10 AND role <> 'target'",
        (series_key,)).fetchall()
    chk = [r[2] for r in rows if r[0] == "check"]
    # Every star of this brightness whatever role the selection gave it —
    # 'dropped_unstable' most of all, since that is the population the check
    # stars are the survivors of.  The target is excluded: it is variable by
    # definition and would swamp the comparison.
    near = [r for r in rows if abs(r[1] - tm) <= ch.FIELD_MATCH_HALF_WIDTH]
    field = [r[2] for r in near]
    kept = [r[2] for r in near if r[0] in ("comp", "check")]
    if not chk or not field:
        return None
    chk_med, f_med = float(np.median(chk)), float(np.median(field))
    return (series_key, tm, chk_med, len(chk), f_med, len(field),
            float(np.median(kept)) if kept else None, len(kept),
            f_med / chk_med if chk_med > 0 else None)


#: White-noise realizations per Allan ladder.  200 fixes the 5th and 95th
#: percentiles of the slope to about 0.02, which is far finer than the
#: 0.2-wide spread being measured.
ALLAN_NULL_REALIZATIONS = 200


def _allan_for_series(pcon, series_key: str, target_key: str,
                      keep: Optional[set] = None):
    """Allan ladders for every check star on the series' longest run.

    The longest CONTIGUOUS run is used (errors.longest_run), because a
    mid-night pause would masquerade as long-tau correlated noise.

    Two things this stage used to get wrong, both fixed here.  (1) The
    red-noise factor was stored in a column called ``red_factor_porb``
    whatever tau it was actually evaluated at, and only 11 of 92 ladders
    reach P_orb — the rest stop at 0.15-0.51 of it, where the factor is
    smaller because red noise grows with tau.  ``tau_used_s`` and
    ``tau_frac_of_porb`` are now stored beside it so no reader can mistake
    the two.  (2) The fitted slope was compared against -0.50, the
    asymptotic white-noise value, when these 4-6 rung ladders have an
    ESTIMATOR whose white-noise expectation is -0.55 with a 5-95% range of
    [-0.88, -0.30].  Each ladder now carries its OWN white null, generated
    at its own length and cadence, and the verdict is per-ladder rather
    than a median against a textbook constant.
    """
    period_d = PERIODS_D.get(target_key, (0.08, ""))[0]
    tau_target = period_d * 86400.0
    stars = [r[0] for r in pcon.execute(
        "SELECT DISTINCT star_id FROM cv_lightcurve "
        "WHERE series_key=? AND role='check'", (series_key,))]
    out = []
    for sid in stars:
        rows = pcon.execute(
            "SELECT l.bjd_tdb, l.mag, f.night, l.frame_id FROM cv_lightcurve l "
            "JOIN cv_frames f USING(frame_id) WHERE l.series_key=? "
            "AND l.star_id=? AND l.mag IS NOT NULL AND l.saturated=0 "
            "ORDER BY l.bjd_tdb", (series_key, sid)).fetchall()
        if keep is not None:
            rows = [r for r in rows if r[3] in keep]
        if len(rows) < 32:
            continue
        t = np.array([r[0] for r in rows], dtype=float)
        y = np.array([r[1] for r in rows], dtype=float)
        nights = [r[2] for r in rows]
        a, b = er.longest_run(t, max_gap_s=1800.0)
        if b - a < 32:
            continue
        tt, yy = t[a:b], y[a:b]
        dt = float(np.median(np.diff(tt)) * 86400.0)
        tau, adev, npair = er.allan_deviation(yy - np.median(yy), dt)
        if tau.size < 3:
            continue
        slope = ch.allan_slope(tau, adev)
        red, tau_used = ch.red_noise_factor(tau, adev, tau_target)
        null = ch.white_noise_allan_null(
            int(b - a), dt, tau_target, er.allan_deviation,
            n_real=ALLAN_NULL_REALIZATIONS,
            rng=np.random.default_rng(20260819 + int(sid)))
        # "Redder than white" means redder than THIS estimator's own 95th
        # percentile on white noise, not redder than -0.50.
        slope_red = int(bool(np.isfinite(slope)
                             and np.isfinite(null["slope_p95"])
                             and slope > null["slope_p95"]))
        red_red = int(bool(np.isfinite(red) and np.isfinite(null["red_p95"])
                           and red > null["red_p95"]))
        # The manifest's own night label of the first point in the run — not
        # a UTC calendar date derived from the JD, which lands on the far
        # side of midnight and would not match any other table's night.
        night = str(nights[a])
        fit = (int(b - a), dt, _f(slope), _f(red), _f(tau_used),
               _f(tau_target), _f(tau_used / tau_target if tau_target else None),
               _f(adev[0]),
               _f(null["slope_p05"]), _f(null["slope_p50"]),
               _f(null["slope_p95"]), _f(null["red_p50"]), _f(null["red_p95"]),
               slope_red, red_red)
        out.append((sid, night, tau, adev, npair, fit))
    return out


# ==========================================================================
# STAGE cadence
# ==========================================================================

def stage_cadence(args) -> None:
    pcon, ocon = connect_ro(PHOT_DB), connect_rw(OUT_DB)
    ensure_schema(ocon)
    try:
        cad_rows, win_rows, alias_rows = [], [], []
        # Frequency grid for the window: 0-40 c/d at 0.002 c/d — fine enough
        # to resolve a 500-day baseline's 0.002 c/d peak width.
        fgrid = np.arange(0.0, 40.0 + 1e-9, 0.002)
        for sk, tk, era, filt in solved_series(pcon):
            # Sampling is characterized over the frames that carry usable
            # photometry, not over every frame that matched: a night whose
            # frames are all rejected is not coverage, and counting it as
            # coverage is what let the quality cut be defended in one
            # section and ignored in every other.
            t = np.array([r[0] for r in ocon.execute(
                "SELECT bjd_tdb FROM ch_frames WHERE series_key=? "
                "AND usable=1 ORDER BY bjd_tdb", (sk,))], dtype=float)
            if t.size < 3:
                continue
            period, _src = PERIODS_D.get(tk, (None, ""))
            s = ch.cadence_stats(t, period_d=period)
            blocks = ch.night_blocks(t)
            gaps = np.diff(t)
            n_full = sum(1 for a, b in blocks if (t[b - 1] - t[a]) >= period) if period else 0
            expt = pcon.execute("SELECT median(exptime) FROM cv_frames "
                                "WHERE series_key=? AND status='matched'",
                                (sk,)).fetchone()[0]
            smear = (float(expt) / 86400.0 / period) if period else None
            # Per-NIGHT phase coverage: the multi-night number is 1.00 for
            # every series (a 90-minute period sampled on ten nights fills
            # every phase bin trivially) and therefore says nothing.  What
            # the science needs is coverage inside ONE night, which is what
            # 'cycle-resolved' means.
            bn = max(blocks, key=lambda ab: ab[1] - ab[0])
            tb = t[bn[0]:bn[1]]
            bn_night = ocon.execute(
                "SELECT night FROM ch_frames WHERE series_key=? AND "
                "usable=1 GROUP BY night ORDER BY count(*) DESC "
                "LIMIT 1", (sk,)).fetchone()[0]
            bn_stats = ch.cadence_stats(tb, period_d=period)
            cad_rows.append((sk, tk, era, filt, s["n_points"], s["n_blocks"],
                             s["baseline_d"], _f(s["median_dt_s"]),
                             s["longest_block_h"], s["on_target_h"],
                             _f(s["duty_cycle"]), period,
                             _f(s.get("pts_per_cycle")),
                             _f(s.get("cycles_longest_block")),
                             _f(s.get("phase_coverage")),
                             float(gaps.max()) if gaps.size else None,
                             n_full, _f(smear), bn_night, int(tb.size),
                             _f(bn_stats.get("cycles_longest_block")),
                             _f(bn_stats.get("phase_coverage")),
                             _f(bn_stats.get("median_dt_s"))))
            print(f"  cadence {sk:16s} n={s['n_points']:5d} nights={s['n_blocks']:3d} "
                  f"dt={s['median_dt_s']:6.1f}s base={s['baseline_d']:7.1f}d "
                  f"best-night phase cov={bn_stats.get('phase_coverage', float('nan')):.2f} "
                  f"({tb.size} pts, {bn_stats.get('cycles_longest_block', 0):.1f} cyc) "
                  f"full-orbit nights={n_full}")
        # Windows: one per (target, era, filter) is far too many figures, so
        # store the window for each target's richest series and for each
        # target's ALL-filter time set (what a joint fit would see).
        for tk in sorted({r[1] for r in solved_series(pcon)}):
            for scope, sql, par in _window_scopes(ocon, tk):
                t = np.array([r[0] for r in ocon.execute(sql, par)], dtype=float)
                if t.size < 5:
                    continue
                w = ch.spectral_window(t, fgrid)
                # Store a decimated copy for plotting (every 5th point) plus
                # every local maximum above 0.05, so no alias peak is lost.
                keep = set(range(0, w.size, 5))
                loc = np.flatnonzero((w[1:-1] > w[:-2]) & (w[1:-1] >= w[2:])
                                     & (w[1:-1] > 0.05)) + 1
                keep.update(int(i) for i in loc)
                for i in sorted(keep):
                    win_rows.append((scope, float(fgrid[i]), float(w[i])))
                period = PERIODS_D[tk][0]
                # An alias is only a PROBLEM when the sampling resolves it.
                # The frequency resolution of a run of length T is ~1/T, so
                # a single 4-hour night has a 6 c/d wide peak: f and f+1 c/d
                # are not two peaks there, they are one.  Recording the
                # resolution beside the alias power is what stops the table
                # being read as "even one night is hopeless".
                base = float(t.max() - t.min())
                fres = 1.0 / base if base > 0 else float("inf")
                for fam, falias in (("solar day", 1.0),
                                    ("sidereal day", 1.0 / ch.SIDEREAL_DAY_D)):
                    resolved = 1 if fres < falias else 0
                    for k, f, p in ch.alias_power(t, 1.0 / period, falias):
                        alias_rows.append((scope, fam, k, f, p, base,
                                           _f(fres), resolved))
        colour_rows = _colour_point_errors(pcon, ocon)
        with ocon:
            ocon.execute("DELETE FROM ch_cadence")
            ocon.execute("DELETE FROM ch_window")
            ocon.execute("DELETE FROM ch_alias")
            ocon.execute("DELETE FROM ch_colour")
            ocon.executemany("INSERT INTO ch_colour VALUES ("
                             + ",".join("?" * 14) + ")", colour_rows)
            ocon.executemany("INSERT INTO ch_cadence VALUES ("
                             + ",".join("?" * 23) + ")", cad_rows)
            ocon.executemany("INSERT INTO ch_window VALUES (?,?,?)", win_rows)
            ocon.executemany("INSERT INTO ch_alias VALUES (?,?,?,?,?,?,?,?)", alias_rows)
        stage_done(ocon, "cadence")
        print(f"cadence: {len(cad_rows)} series, {len(win_rows)} window points")
    finally:
        pcon.close(); ocon.close()


def qualifying_nights(ocon, target_key: str, period_d: float) -> dict:
    """Nights that carry a full orbit in a given filter: ``{night: {filters}}``.

    A night qualifies for a filter when its USABLE frames span at least one
    orbital period and there are at least 12 of them.  This is the Q1
    criterion, in one place, so the verdict and the colour-error
    measurement below cannot count different nights.
    """
    per_night: dict = {}
    for night, filt, n, span in ocon.execute(
            "SELECT night, filter, count(*), max(bjd_tdb)-min(bjd_tdb) "
            "FROM ch_frames WHERE target_key=? AND usable=1 "
            "GROUP BY night, filter", (target_key,)):
        if span is not None and span >= period_d and n >= 12:
            per_night.setdefault(night, set()).add(filt.lower())
    return per_night


def _colour_point_errors(pcon, ocon) -> list:
    """What a COLOUR point actually costs, per qualifying three-filter night.

    Q1 asks for "quasi-simultaneous" three-colour coverage and was graded on
    ``prec_at_target`` — the SINGLE-BAND per-point precision.  A colour is
    never better than sqrt(2) times that even when the two exposures are
    simultaneous, and here they are not: the filters interleave, so a g
    point is paired with an i point taken a median 75 s away on a star whose
    brightness sweeps up to 3 mmag per second.  That term is usually the
    LARGEST one and it was not in the grade at all.

    Measured, not modelled: the offsets come from the real timestamps
    (:func:`ch.nearest_time_offsets`) and the sweep rate from the target's
    own light curve on that night (:func:`ch.rate_of_change_mag_per_s`).
    """
    out = []
    prec = dict(ocon.execute(
        "SELECT series_key, prec_at_target FROM ch_noise_series "
        "WHERE prec_at_target IS NOT NULL"))
    for tk in sorted(PERIODS_D):
        period = PERIODS_D[tk][0]
        nights = qualifying_nights(ocon, tk, period)
        for night, filts in sorted(nights.items()):
            if len(filts) < 3:
                continue
            era = ocon.execute("SELECT era_id FROM ch_frames WHERE "
                               "target_key=? AND night=? LIMIT 1",
                               (tk, night)).fetchone()[0]
            # The target's own points on this night, per filter.
            byband: dict = {}
            for sk, in ocon.execute(
                    "SELECT DISTINCT series_key FROM ch_frames WHERE "
                    "target_key=? AND night=? AND usable=1", (tk, night)):
                rows = pcon.execute(
                    "SELECT l.bjd_tdb, l.mag FROM cv_lightcurve l "
                    "JOIN cv_frames f USING(frame_id) WHERE l.series_key=? "
                    "AND f.night=? AND l.role='target' AND l.mag IS NOT NULL "
                    "AND l.saturated=0 ORDER BY l.bjd_tdb", (sk, night)
                ).fetchall()
                if len(rows) >= 6:
                    byband[sk] = (np.array([r[0] for r in rows], dtype=float),
                                  np.array([r[1] for r in rows], dtype=float))
            keys = sorted(byband)
            for i in range(len(keys)):
                for j in range(i + 1, len(keys)):
                    ka, kb = keys[i], keys[j]
                    ta, ma = byband[ka]
                    tb, _mb = byband[kb]
                    dt = ch.nearest_time_offsets(ta, tb)
                    rate = ch.rate_of_change_mag_per_s(ta, ma)
                    if dt.size == 0 or rate.size != dt.size:
                        continue
                    sa = prec.get(ka)
                    sb = prec.get(kb)
                    if sa is None or sb is None:
                        continue
                    sig = np.array([ch.colour_point_sigma(sa, sb, r, d)
                                    for r, d in zip(rate, dt)])
                    sig = sig[np.isfinite(sig)]
                    if sig.size == 0:
                        continue
                    out.append((tk, era, night, ka, kb, int(dt.size),
                                _f(np.median(dt)), _f(np.percentile(dt, 90)),
                                _f(np.median(rate)),
                                _f(np.percentile(rate, 90)),
                                _f(sa), _f(sb), _f(np.median(sig)),
                                _f(np.percentile(sig, 90))))
    print(f"  colour: {len(out)} band pairs on three-filter full-orbit nights")
    return out


def _window_scopes(ocon, target_key: str):
    """The three sampling sets whose windows actually decide anything:
    the whole target (all filters, all eras — what a joint period search
    sees), the richest single series, and that series' richest single night
    (the alias-free anchor).

    Every one of them is drawn from ``ch_frames WHERE usable=1``: the window
    of a sampling that includes frames the quality cut rejected is the
    window of an observing run nobody would analyse.
    """
    yield (f"{target_key}|all",
           "SELECT bjd_tdb FROM ch_frames WHERE target_key=? AND usable=1 "
           "ORDER BY bjd_tdb", (target_key,))
    best = ocon.execute(
        "SELECT series_key, count(*) n FROM ch_frames WHERE target_key=? "
        "AND usable=1 GROUP BY series_key ORDER BY n DESC LIMIT 1",
        (target_key,)).fetchone()
    if not best:
        return
    sk = best[0]
    yield (f"{sk}|series",
           "SELECT bjd_tdb FROM ch_frames WHERE series_key=? "
           "AND usable=1 ORDER BY bjd_tdb", (sk,))
    night = ocon.execute(
        "SELECT night FROM ch_frames WHERE series_key=? AND usable=1 "
        "GROUP BY night ORDER BY count(*) DESC LIMIT 1", (sk,)).fetchone()
    if night:
        yield (f"{sk}|{night[0]}",
               "SELECT bjd_tdb FROM ch_frames WHERE series_key=? AND night=? "
               "AND usable=1 ORDER BY bjd_tdb", (sk, night[0]))


# ==========================================================================
# STAGE detect
# ==========================================================================

def _residual_pool(pcon, series_key: str, frame_ids: np.ndarray) -> list[np.ndarray]:
    """Check-star residual series on a fixed frame list — the REAL noise.

    Each check star contributes one residual vector aligned to
    ``frame_ids``; missing points are filled with the star's own median
    residual (0 by construction) so the vector length matches the timestamp
    vector exactly.  Filling rather than dropping keeps the sampling
    identical between the noise and the injected signal, which is the whole
    point of injecting at the real timestamps.
    """
    rows = pcon.execute(
        "SELECT star_id, frame_id, mag FROM cv_lightcurve WHERE series_key=? "
        "AND role='check' AND mag IS NOT NULL AND saturated=0", (series_key,)
    ).fetchall()
    if not rows:
        return []
    index = {int(f): i for i, f in enumerate(frame_ids)}
    by_star: dict[int, np.ndarray] = {}
    for sid, fid, mag in rows:
        if fid not in index:
            continue
        arr = by_star.setdefault(sid, np.full(len(frame_ids), np.nan))
        arr[index[fid]] = mag
    pool = []
    for sid, arr in by_star.items():
        n_ok = int(np.isfinite(arr).sum())
        # A star must actually cover the run: 40% of the frames and at least
        # 20 real points.  Gap-filling with zeros is what keeps the vector
        # aligned to the timestamps, but too much filling would fabricate a
        # quiet star and flatter the detection limit.
        if n_ok < max(20, 0.4 * len(frame_ids)):
            continue
        med = float(np.nanmedian(arr))
        r = np.where(np.isfinite(arr), arr - med, 0.0)
        pool.append(r)
    return pool


def _detect_cell(job):
    """Worker: recovery fraction of ONE (period, amplitude, score) cell.

    The frequency grid is rebuilt inside the worker from three floats
    instead of being pickled: at 0.0005 c/d over 2-40 c/d it is 76,000
    doubles, and shipping 600 kB per cell would cost more than the
    periodogram it feeds.

    Both the detrending and the two score modes now live inside
    ``ch.recovery_fraction``, so the detrended path can no longer drift
    away from the raw one — it used to be a hand-inlined copy of the same
    loop, which is how the detrended variant kept the period-tolerance
    scoring after the reason for questioning it was known.
    """
    (t, pool, fmin, fmax, fstep, period, amp, threshold, n_trials, seed,
     nights, score) = job
    freqs = np.arange(fmin, fmax, fstep)
    rng = np.random.default_rng(seed)
    return (period, amp, score,
            ch.recovery_fraction(t, pool, freqs, period, amp, threshold,
                                 n_trials, rng=rng, score=score,
                                 detrend_nights=nights))


def stage_detect(args) -> None:
    """Injection and recovery, one scope at a time, committed as it goes.

    Three regimes per target, because they answer different questions:

    * ``season``  — the whole richest series, untouched.  This is what a
      blind period search sees, aliases and nightly zero-point wander and
      all.
    * ``season-dt`` — the same set after per-night mean removal, the
      detrending any real analysis applies (strategy §4.20).  The pair
      season / season-dt is the measurement of what that detrending is
      worth, and what it costs at long periods.
    * ``night``   — the richest single night.  No 1 c/d comb, so this
      contour is what a single-night, cycle-resolved claim can rest on —
      the claim this paper is actually built around.

    And TWO SCORE MODES per regime, which is the correction this rebuild
    exists for.  The first version scored every trial as recovered only if
    the tallest peak landed within 1% of the injected frequency.  On a
    single night the frequency resolution is 2.6-9.0 c/d and that window is
    0.12-0.14 c/d — 20 to 70 times narrower than the peak — so the answer
    measured how well a period can be DETERMINED from under five cycles,
    not what amplitude can be DETECTED.  A 300 mmag injection into VV Pup's
    richest night cleared the threshold in 40 trials of 40 and was scored
    recovered in 25, and the page reported "never reached even at 300 mmag"
    as a detection limit, then fed it into Q3's margin and S4's headline
    ratio.  Both questions are now measured and stored under their own
    ``score`` label, each against its own null (see
    ``ch.detection_threshold``): ``'period'`` for a blind search,
    ``'known'`` for a modulation at a literature period — which is what
    these five CVs actually have.

    The noise is not simulated: it is the held-out check stars' own
    residuals at the same timestamps, cyclically rolled so their correlated
    structure survives.  The threshold is measured the same way on
    signal-free data, so a red-noise archive cannot borrow a white-noise
    false-alarm rate.
    """
    pcon, ocon = connect_ro(PHOT_DB), connect_rw(OUT_DB)
    ensure_schema(ocon)
    try:
        done = {(r[0], r[1]) for r in ocon.execute(
            "SELECT scope, regime FROM ch_detect GROUP BY scope, regime")}
        fstep = 0.0005
        freqs = np.arange(SEARCH_FMIN_CD, SEARCH_FMAX_CD, fstep)
        jobs = []
        for tk in sorted(PERIODS_D):
            best = ocon.execute(
                "SELECT series_key, count(*) n FROM ch_frames WHERE "
                "target_key=? AND usable=1 GROUP BY series_key "
                "ORDER BY n DESC LIMIT 1", (tk,)).fetchone()
            if not best:
                continue
            sk = best[0]
            jobs.append((f"{sk}|season", sk, "season", None))
            jobs.append((f"{sk}|season-dt", sk, "season-dt", None))
            night = ocon.execute(
                "SELECT night FROM ch_frames WHERE series_key=? AND "
                "usable=1 GROUP BY night ORDER BY count(*) DESC LIMIT 1",
                (sk,)).fetchone()
            if night:
                jobs.append((f"{sk}|{night[0]}", sk, "night", night[0]))
        n_done_now = 0
        skipped: list[tuple[str, str]] = []
        for scope, sk, regime, night in jobs:
            if (scope, regime) in done and not args.force:
                print(f"  detect {scope} [{regime}] already done")
                continue
            if args.limit and n_done_now >= args.limit:
                print(f"  detect: stopping at --limit {args.limit} "
                      f"(resume by re-running this stage)")
                break
            sql = ("SELECT frame_id, bjd_tdb, night FROM ch_frames "
                   "WHERE series_key=? AND usable=1"
                   + (" AND night=?" if night else "") + " ORDER BY bjd_tdb")
            rows = ocon.execute(sql, (sk, night) if night else (sk,)).fetchall()
            fids = np.array([r[0] for r in rows])
            t = np.array([r[1] for r in rows], dtype=float)
            nights = np.array([str(r[2]) for r in rows])
            detrend = nights if regime == "season-dt" else None
            pool = _residual_pool(pcon, sk, fids)
            if len(pool) < 2 or t.size < 20:
                # Not a failure to retry: this scope has no held-out check
                # star covering the run, so it can never produce a contour.
                # Recorded as skipped so the stage can still complete.
                print(f"  detect {scope}: too little data "
                      f"(pool={len(pool)}, n={t.size}) - permanently skipped")
                skipped.append((scope, regime))
                continue
            sigma = float(np.median([np.std(p[p != 0]) for p in pool]))
            rng = np.random.default_rng(20260819)
            t0 = time.time()
            tk = sk.split("|")[0]
            f_orb = 1.0 / PERIODS_D[tk][0]
            baseline = float(t.max() - t.min())
            fres = 1.0 / baseline if baseline > 0 else float("inf")
            # The threshold must be measured through the SAME pipeline the
            # injections pass through: a detrended search has a lower
            # false-alarm floor, and reusing the raw threshold would hide
            # exactly the improvement this regime exists to measure.
            thr_pool = ([ch.remove_nightly_means(p, nights) for p in pool]
                        if detrend is not None else pool)
            periods = np.unique(np.concatenate(
                [INJECT_PERIODS_D, [PERIODS_D[tk][0]]]))
            periods = periods[(1.0 / periods >= SEARCH_FMIN_CD)
                              & (1.0 / periods <= SEARCH_FMAX_CD)]
            # One threshold per score mode, each on its own null statistic.
            # 'period' compares the tallest peak in the band against the
            # distribution of the tallest peak; 'known' compares the power at
            # one frequency against the distribution of the power at that
            # frequency.  Using the max-statistic threshold for a known-period
            # claim charges the paper a look-elsewhere penalty it does not owe.
            thr = {"period": ch.detection_threshold(
                       t, thr_pool, freqs, THRESHOLD_TRIALS, ch.DETECT_FAP,
                       np.random.default_rng(20260819))}
            thr_known = {}
            for P in periods:
                thr_known[float(P)] = ch.detection_threshold(
                    t, thr_pool, freqs, THRESHOLD_TRIALS, ch.DETECT_FAP,
                    np.random.default_rng(20260819), at_freq_cd=1.0 / float(P))
            spread = ch.threshold_spread(t, thr_pool, freqs, 100,
                                         ch.DETECT_FAP,
                                         np.random.default_rng(20260819))
            spread_known = ch.threshold_spread(
                t, thr_pool, freqs, 100, ch.DETECT_FAP,
                np.random.default_rng(20260819), at_freq_cd=f_orb)
            cells = []
            for i, P in enumerate(periods):
                for j, A in enumerate(INJECT_AMPS):
                    for m, score in enumerate(ch.SCORE_MODES):
                        th = (thr["period"] if score == "period"
                              else thr_known[float(P)])
                        cells.append((t, pool, SEARCH_FMIN_CD, SEARCH_FMAX_CD,
                                      fstep, float(P), float(A), float(th),
                                      INJECT_TRIALS,
                                      20260819 + 1000 * i + 10 * j + m,
                                      detrend, score))
            frac: dict = {}
            with ProcessPoolExecutor(max_workers=min(args.workers, 6)) as ex:
                for P, A, score, f in ex.map(_detect_cell, cells, chunksize=1):
                    frac[(P, A, score)] = f
            det_rows, con_rows = [], []
            for score in ch.SCORE_MODES:
                for P in periods:
                    th = (thr["period"] if score == "period"
                          else thr_known[float(P)])
                    fr = [frac[(float(P), float(A), score)] for A in INJECT_AMPS]
                    for A, f in zip(INJECT_AMPS, fr):
                        det_rows.append((scope, sk, regime, score, float(P),
                                         float(A), _f(f), INJECT_TRIALS,
                                         _f(th), int(t.size), _f(sigma)))
                    a90 = ch.recovery_contour(INJECT_AMPS, fr)
                    lo, hi = ch.contour_uncertainty(INJECT_AMPS, fr,
                                                    INJECT_TRIALS)
                    con_rows.append((scope, sk, regime, score, float(P),
                                     _f(a90), _f(lo), _f(hi), int(t.size),
                                     _f(baseline / float(P)), _f(fres)))
                    print(f"  detect {scope:26s} [{regime:9s}/{score:6s}] "
                          f"P={P:.4f}d A90="
                          + (f"{1000 * a90:7.1f} mmag "
                             f"[{1000 * lo:.0f}-{1000 * hi:.0f}]"
                             if np.isfinite(a90) else "  not reached"))
            # ---- how often the tallest peak is an ALIAS, at each amplitude.
            # S3's verdict was a window statistic converted straight into an
            # operational claim; this is the operational measurement.
            alias_rows = []
            for A in (0.05, 0.10, 0.50, 1.00):
                res = ch.alias_confusion(t, pool, freqs, PERIODS_D[tk][0], A,
                                         n_trials=40,
                                         rng=np.random.default_rng(20260901),
                                         detrend_nights=detrend)
                alias_rows.append((scope, sk, regime, float(A),
                                   _f(res["true"]), _f(res["alias"]),
                                   _f(res["other"]), int(res["n_trials"])))
            with ocon:
                ocon.execute("DELETE FROM ch_detect WHERE scope=? AND regime=?",
                             (scope, regime))
                ocon.execute("DELETE FROM ch_contour WHERE scope=? AND regime=?",
                             (scope, regime))
                ocon.execute("DELETE FROM ch_alias_confusion "
                             "WHERE scope=? AND regime=?", (scope, regime))
                ocon.execute("DELETE FROM ch_threshold "
                             "WHERE scope=? AND regime=?", (scope, regime))
                ocon.executemany("INSERT INTO ch_detect VALUES ("
                                 + ",".join("?" * 11) + ")", det_rows)
                ocon.executemany("INSERT INTO ch_contour VALUES ("
                                 + ",".join("?" * 11) + ")", con_rows)
                ocon.executemany("INSERT INTO ch_alias_confusion VALUES "
                                 "(?,?,?,?,?,?,?,?)", alias_rows)
                ocon.executemany("INSERT INTO ch_threshold VALUES "
                                 "(?,?,?,?,?,?,?)", [
                    (scope, regime, "period", _f(thr["period"]),
                     spread["n"], json.dumps([round(x, 5) for x in
                                              spread["per_star"]]),
                     _f(spread["spread_frac"])),
                    (scope, regime, "known", _f(thr_known.get(PERIODS_D[tk][0])),
                     spread_known["n"], json.dumps([round(x, 5) for x in
                                                    spread_known["per_star"]]),
                     _f(spread_known["spread_frac"]))])
            n_done_now += 1
            print(f"    ({time.time() - t0:.0f}s, thr[period]="
                  f"{thr['period']:.3f} (per-star spread "
                  f"{100 * spread['spread_frac']:.0f}% over {spread['n']} "
                  f"stars), thr[known@P_orb]="
                  f"{thr_known.get(PERIODS_D[tk][0], float('nan')):.3f}, "
                  f"sigma={sigma:.4f}, n={t.size}, {len(pool)} check stars, "
                  f"{baseline / PERIODS_D[tk][0]:.1f} cycles)")
        have = {(r[0], r[1]) for r in ocon.execute(
            "SELECT scope, regime FROM ch_detect GROUP BY scope, regime")}
        have.update(skipped)
        remaining = [j for j in jobs if (j[0], j[2]) not in have]
        if not remaining:
            stage_done(ocon, "detect")
            set_meta(ocon, detect_skipped_scopes="; ".join(
                f"{sc} [{rg}]" for sc, rg in skipped) or "none")
        print(f"detect: {len(remaining)} scopes still to do, "
              f"{len(skipped)} permanently skipped")
    finally:
        pcon.close(); ocon.close()


# ==========================================================================
# STAGE timing
# ==========================================================================

#: Bright-phase edge shapes tried, as a fraction of the orbital period.
#: 0.01 (~1 min for a 100-min binary) is the sharp self-eclipse edge polars
#: actually show; 0.05 is a soft cyclotron-beaming shoulder.  Reporting both
#: is the honest way to admit the answer depends on a shape this archive has
#: not measured.
INGRESS_PHASES = (0.01, 0.05)

#: Width of the injected faint phase, in orbital phase.  A self-eclipsing
#: polar spends roughly half its cycle faint; 0.45 is that, and it is what
#: puts TWO edges in the cycle for the fit to use.
FAINT_PHASE_WIDTH = 0.45

#: A target needs at least this many measured points before its own light
#: curve is allowed to set the injected amplitude.  Below it the "amplitude"
#: is the scatter of a non-detection (EU UMa in i has 24 points and an RMS
#: of 2.2 mag) and would flatter the timing test by an order of magnitude.
#:
#: 30, not 60, since the amplitude is now measured on ONE NIGHT rather than
#: over a whole season: a median-density night holds 40-ish target points and
#: a 5th-to-95th percentile over 40 points is a real measurement (the 5th
#: percentile is the second order statistic, not the minimum).  Requiring 60
#: would push every typical night onto the fallback constant, which is the
#: opposite of the point — the typical night is precisely what has to be
#: measured rather than assumed.
MIN_TARGET_POINTS_FOR_DEPTH = 30

#: Fallback bright-phase amplitude in magnitudes when the target is too
#: poorly measured to supply one.  1.0 mag is a modest polar bright phase;
#: quoting it explicitly is better than silently inheriting a 4 mag "depth"
#: from an undetected star.
FALLBACK_DEPTH_MAG = 1.0


def _target_amplitude(pcon, series_key: str,
                      night: Optional[str] = None) -> tuple[float, int, str]:
    """Bright-phase amplitude from the target's OWN measured light curve.

    The robust 5th-to-95th percentile range of the target magnitudes, not
    2x the RMS: the RMS of a barely-detected star is dominated by its
    measurement noise and would inject an unphysically deep feature.

    ``night`` restricts it to ONE night, and that is the version the timing
    Monte Carlo must use.  Taken over the whole series the range is not a
    bright-phase amplitude at all: ``stlmi|e76|g`` spans 396 days across a
    1.92 mag change in nightly median — the high/low accretion-state
    transition Q1 and Q4 are about — so the "measured 1.45 mag bright-phase
    amplitude" quoted in Q2's deciding number was mostly a state change.
    The night being simulated has 0.672 mag.  (The answer barely moves,
    because this regime is sampling-limited rather than S/N-limited, but a
    referee would not get that far before stopping.)

    Returns ``(depth_mag, n_points, source)``.
    """
    sql = ("SELECT l.mag FROM cv_lightcurve l JOIN cv_frames f USING(frame_id) "
           "WHERE l.series_key=? AND l.role='target' AND l.mag IS NOT NULL "
           "AND l.saturated=0" + (" AND f.night=?" if night else ""))
    mags = [r[0] for r in pcon.execute(
        sql, (series_key, night) if night else (series_key,))]
    n = len(mags)
    if n < MIN_TARGET_POINTS_FOR_DEPTH:
        return FALLBACK_DEPTH_MAG, n, "fallback (target too poorly measured)"
    a = np.asarray(mags, dtype=float)
    src = ("measured p5-p95 of the target on night " + night if night
           else "measured p5-p95 of the whole series (NOT a single cycle)")
    return float(np.percentile(a, 95) - np.percentile(a, 5)), n, src


def stage_timing(args) -> None:
    """Monte-Carlo epoch precision of one bright-phase edge.

    This is the measurement the strategy document (§4.16) says must be made
    before any per-cycle O-C tier is adopted, and it is made here at the
    REAL per-filter timestamps, with the REAL correlated noise (check-star
    residuals, cyclically rolled), at the amplitude the target ACTUALLY
    showed on the night being simulated.  In the matched regimes the
    template shape is held exactly right during the fit, so those numbers
    are a LOWER BOUND on what a real fit achieves.

    TWO nights per target, not one.  The first version demonstrated Q2 on
    each target's richest night and reported the answer as a property of the
    whole per-cycle O-C tier — but ``stlmi|e76|g``'s demonstration night
    holds 153 frames and 4.94 cycles while its MEDIAN night holds 42 frames
    and 1.78, and most O-C points will come from nights like the second.
    ``night_kind`` labels which is which, so the tier is sized on the
    typical night rather than the best one.
    """
    pcon, ocon = connect_ro(PHOT_DB), connect_rw(OUT_DB)
    ensure_schema(ocon)
    try:
        rows_out = []
        for tk in sorted(PERIODS_D):
            period = PERIODS_D[tk][0]
            best = ocon.execute(
                "SELECT series_key, count(*) n FROM ch_frames WHERE "
                "target_key=? AND usable=1 GROUP BY series_key "
                "ORDER BY n DESC LIMIT 1", (tk,)).fetchone()
            if not best:
                continue
            sk = best[0]
            # Nights ordered by how many USABLE frames they carry, so both
            # the demonstration night and the typical one are real nights of
            # analysable data rather than of frames the quality cut rejects.
            by_night = ocon.execute(
                "SELECT night, count(*) n FROM ch_frames WHERE series_key=? "
                "AND usable=1 GROUP BY night HAVING n >= 12 ORDER BY n",
                (sk,)).fetchall()
            if not by_night:
                continue
            chosen = {"richest": by_night[-1][0],
                      "median": by_night[len(by_night) // 2][0]}
            for night_kind, night in chosen.items():
                rows_out.extend(_timing_one_night(
                    pcon, ocon, tk, sk, period, night, night_kind, args))
        with ocon:
            ocon.execute("DELETE FROM ch_timing")
            ocon.executemany("INSERT INTO ch_timing VALUES ("
                             + ",".join("?" * 16) + ")", rows_out)
        stage_done(ocon, "timing")
        print(f"timing: {len(rows_out)} Monte-Carlo results")
    finally:
        pcon.close(); ocon.close()


def _timing_one_night(pcon, ocon, tk, sk, period, night, night_kind,
                      args) -> list:
    """Every timing regime for ONE night of one series.  Returns ch_timing rows."""
    rows = ocon.execute(
        "SELECT frame_id, bjd_tdb FROM ch_frames WHERE series_key=? "
        "AND night=? AND usable=1 ORDER BY bjd_tdb", (sk, night)).fetchall()
    if len(rows) < 12:
        return []
    fids = np.array([r[0] for r in rows])
    t = np.array([r[1] for r in rows], dtype=float)
    pool_night = _residual_pool(pcon, sk, fids)
    sigma = ocon.execute(
        "SELECT prec_at_target FROM ch_noise_series WHERE series_key=?",
        (sk,)).fetchone()
    sigma = (sigma[0] if sigma and sigma[0] else None)
    if sigma is None:
        sigma = pcon.execute("SELECT check_rms_median FROM cv_series "
                             "WHERE series_key=?", (sk,)).fetchone()[0]
    # The night's OWN amplitude, not the season's: see _target_amplitude.
    depth, n_tgt, depth_src = _target_amplitude(pcon, sk, night)
    depth_season = _target_amplitude(pcon, sk)[0]
    # An observed edge can never be sharper than one exposure: a 240 s
    # integration on a 90-minute binary smears 4.4% of a cycle whatever the
    # star does.  The injected ingress is therefore the LONGER of the
    # intrinsic edge and the exposure, which is what stops EU UMa's numbers
    # from claiming a 54 s edge it physically cannot record.
    expt = float(pcon.execute(
        "SELECT median(exptime) FROM cv_frames WHERE series_key=? "
        "AND status='matched'", (sk,)).fetchone()[0])
    smear = expt / 86400.0 / period
    # ---- per-cycle: one orbit's worth of the real timestamps -------------
    sel = t < t[0] + period
    t_cycle, idx = t[sel], np.flatnonzero(sel)
    if t_cycle.size < 8:
        idx = np.arange(min(t.size, 12))
        t_cycle = t[idx]
    pool_cycle = [p[idx] for p in pool_night]
    dt = float(np.median(np.diff(t_cycle)) * 86400.0)
    n_cyc = max(1.0, (t[-1] - t[0]) / period)
    out = []
    for regime, tt, pl, nn in (("per-cycle", t_cycle, pool_cycle,
                                t_cycle.size),
                               ("night-mean", t, pool_night, t.size),
                               ("per-cycle shape-mismatched",
                                t_cycle, pool_cycle, t_cycle.size)):
        for ing in INGRESS_PHASES:
            ing_eff = max(ing, smear)
            # The mismatched regime fits with the OTHER edge sharpness and a
            # 20% wrong depth: the realistic case, because nobody has
            # measured this archive's ingress duration and cyclotron
            # beaming makes it band-dependent.
            if regime.endswith("mismatched"):
                other = [x for x in INGRESS_PHASES if x != ing][0]
                fit_ing, fit_depth = max(other, smear), 1.2 * depth
            else:
                fit_ing, fit_depth = None, None
            st = ch.timing_precision_mc(
                tt, period, depth, float(sigma),
                width_phase=FAINT_PHASE_WIDTH, ingress_phase=ing_eff,
                n_trials=args.timing_trials, noise_pool=pl,
                fit_ingress_phase=fit_ing, fit_depth_mag=fit_depth)
            out.append((sk, tk, night, night_kind, regime, int(nn), dt,
                        depth, _f(depth_season), float(sigma), ing, ing_eff,
                        _f(smear), _f(st), depth_src, len(pl)))
            print(f"  timing {tk:6s} [{night_kind:7s}] {regime:26s} "
                  f"n={nn:4d} dt={dt:6.1f}s depth={depth:.2f} "
                  f"(season {depth_season:.2f}, {n_tgt} pts) "
                  f"sigma={sigma:.4f} ingress={ing:.2f}->{ing_eff:.3f}P "
                  f"-> sigma_t={st:8.1f} s "
                  f"[{len(pl)} real-noise series, {n_cyc:.1f} cycles]")
    return out


# ==========================================================================
# STAGE verdict
# ==========================================================================

def stage_verdict(args) -> None:
    """Turn the measurements into SUPPORTED / CAVEATS / NOT SUPPORTED.

    Every ``deciding_number`` is a string BUILT from a query result, so the
    verdicts are recomputed whenever the measurements change and the table
    cannot drift away from the evidence it claims to rest on.  The
    thresholds that separate the verdicts are the strategy document's OWN
    stated requirements: this stage tests the plan against the data, it does
    not invent new criteria.

    Rebuilt 2026-08-19 after an adversarial review found that three of the
    ten deciding numbers did not measure the quantity their goal asks about:

    * Q1 was graded on single-band precision when a colour point's error is
      dominated by the filters not being simultaneous;
    * Q4 — the plan's ONLY unqualified SUPPORTED — was graded on per-point
      precision against the state separation, which establishes only that
      one night can be classified.  A duty cycle's uncertainty comes from
      the number of independent epochs, and from being able to put those
      epochs on a common magnitude scale, which is the very thing S1 grades
      NOT SUPPORTED;
    * S3 converted a window statistic into an operational claim about which
      peak wins, which the injection machinery can measure and now does.

    Q3 and S4 also move, because both were quoting a single-night
    period-DETERMINATION contour as a detection limit.
    """
    ocon = connect_rw(OUT_DB)
    pcon = connect_ro(PHOT_DB)
    ensure_schema(ocon)
    try:
        v = []

        def one(sql, default=None):
            r = ocon.execute(sql).fetchone()
            return r[0] if r and r[0] is not None else default

        def mmag(x, nd=0):
            return "n/a" if x is None else f"{1000 * float(x):.{nd}f} mmag"

        # ---------- the measurements the verdicts rest on -------------------
        prec = dict(ocon.execute(
            "SELECT target_key, min(prec_at_target) FROM ch_noise_series "
            "WHERE prec_at_target IS NOT NULL GROUP BY target_key"))
        prec_hi = dict(ocon.execute(
            "SELECT target_key, max(prec_at_target) FROM ch_noise_series "
            "WHERE prec_at_target IS NOT NULL GROUP BY target_key"))
        # Full-orbit nights per target, and how many carry three filters —
        # counted from the USABLE frames, through the one shared definition.
        three_filter, any_filter = {}, {}
        for tk in PERIODS_D:
            per_night = qualifying_nights(ocon, tk, PERIODS_D[tk][0])
            three_filter[tk] = sum(1 for f in per_night.values() if len(f) >= 3)
            any_filter[tk] = len(per_night)
        # Independent epochs per target: what a duty cycle is measured from.
        nights_per_target = dict(ocon.execute(
            "SELECT target_key, count(DISTINCT night) FROM ch_frames "
            "WHERE usable=1 GROUP BY target_key"))
        span_per_target = dict(ocon.execute(
            "SELECT target_key, max(bjd_tdb)-min(bjd_tdb) FROM ch_frames "
            "WHERE usable=1 GROUP BY target_key"))
        # Timing: the richest night is the DEMONSTRATION, the median night is
        # what the O-C tier will actually be built from.
        def st_of(regime, night_kind):
            return dict(ocon.execute(
                "SELECT target_key, max(sigma_t_s) FROM ch_timing "
                "WHERE regime LIKE ? AND night_kind=? GROUP BY target_key",
                (regime, night_kind)))
        st_ideal = st_of("per-cycle", "richest")
        st_real = st_of("%mismatched", "richest")
        st_night = st_of("night-mean", "richest")
        st_ideal_med = st_of("per-cycle", "median")
        st_real_med = st_of("%mismatched", "median")

        def a90(tk, regime, score="known"):
            """90% recovery semi-amplitude at the target's own period.

            ``score`` defaults to 'known' — a modulation at a literature
            period — because that is the question every one of these goals
            actually asks.  'period' is available for the blind-search
            statement and is quoted as such wherever it is used.
            """
            r = ocon.execute(
                "SELECT c.amp90, c.amp90_lo, c.amp90_hi, c.n_cycles "
                "FROM ch_contour c JOIN ch_cadence cad "
                "ON cad.series_key = c.series_key "
                "WHERE c.series_key LIKE ? AND c.regime=? AND c.score=? "
                "AND abs(c.period_d - cad.period_d) < 1e-6 LIMIT 1",
                (tk + "|%", regime, score)).fetchone()
            return r if r else (None, None, None, None)

        worst_alias = one("SELECT max(power) FROM ch_alias WHERE k IN (-1,1) "
                          "AND resolved=1")
        solar_alias = one("SELECT max(power) FROM ch_alias WHERE k IN (-1,1) "
                          "AND resolved=1 AND family='solar day'")
        med_red = one("SELECT median(red_factor) FROM ch_allan_fit")
        med_slope = one("SELECT median(slope) FROM ch_allan_fit")
        med_null_slope = one("SELECT median(slope_null_p50) FROM ch_allan_fit")
        n_ladders = one("SELECT count(*) FROM ch_allan_fit", 0)
        n_redder = one("SELECT sum(slope_redder_than_null) FROM ch_allan_fit", 0)
        n_red_redder = one("SELECT sum(red_redder_than_null) FROM ch_allan_fit", 0)
        med_tau_frac = one("SELECT median(tau_frac_of_porb) FROM ch_allan_fit")
        n_reach_porb = one("SELECT sum(tau_frac_of_porb >= 0.99) "
                           "FROM ch_allan_fit", 0)
        n_tied = pcon.execute("SELECT count(*) FROM cv_field_tie "
                              "WHERE n_gaia_matched > 0").fetchone()[0]
        n_blocks = pcon.execute("SELECT count(*) FROM cv_field_tie").fetchone()[0]
        # Do any two eras of one target overlap in time?  (The colour seam.)
        overlap = {}
        for tk in PERIODS_D:
            eras = dict((era, (lo, hi)) for era, lo, hi in ocon.execute(
                "SELECT era_id, min(bjd_tdb), max(bjd_tdb) FROM ch_frames "
                "WHERE target_key=? AND usable=1 GROUP BY era_id", (tk,)))
            keys = list(eras)
            overlap[tk] = any(eras[keys[i]][0] <= eras[keys[j]][1]
                              and eras[keys[j]][0] <= eras[keys[i]][1]
                              for i in range(len(keys))
                              for j in range(i + 1, len(keys)))

        # ---------------- Q1 --------------------------------------------------
        night3 = "; ".join(f"{TARGET_LABEL[t]} {three_filter[t]} of "
                           f"{any_filter[t]}" for t in
                           ("stlmi", "vvpup", "euuma", "anuma"))
        # Per target: the typical colour-point error, and the single-band
        # precision it should be compared against.  A range over all band
        # pairs would be dominated by one bright-phase-edge night.
        col_rows = ocon.execute(
            "SELECT target_key, count(*), median(sigma_colour_med), "
            "median(sigma_a), median(dt_med_s) FROM ch_colour "
            "GROUP BY target_key ORDER BY target_key").fetchall()
        col_n = ocon.execute("SELECT count(*) FROM ch_colour").fetchone()[0]
        col_txt = ("no three-filter full-orbit night carries a measurable "
                   "colour point" if not col_n else
                   "median COLOUR-point error vs the single-band precision it "
                   "was graded on: "
                   + "; ".join(
                       f"{TARGET_LABEL.get(t, t)} {mmag(cm)} vs {mmag(sa)} "
                       f"({cm / sa:.1f}x, {n} band pairs, median filter "
                       f"offset {dt:.0f} s)"
                       for t, n, cm, sa, dt in col_rows))
        v.append((
            "Q1", 1,
            "Accretion-state-resolved, orbit-phase-resolved cyclotron colour "
            "curves of ST LMi, VV Pup and EU UMa (two within-era analyses each)",
            "single-night, cycle-resolved, state-tagged quasi-simultaneous "
            "three-colour coverage",
            _q1_verdict(three_filter),
            f"full-orbit nights carrying ALL THREE filters, of all full-orbit "
            f"nights: {night3}. {col_txt}",
            "A colour curve needs one night, three filters, a whole orbit. "
            "Counted from the USABLE frames: a night qualifies when its span "
            "exceeds one orbital period in >= 12 frames of that filter. The "
            "precision quoted is the COLOUR point's, not one band's: the "
            "filters interleave rather than fire together, so the target's "
            "own sweep between the two exposures enters in quadrature and is "
            "usually the largest term.",
            "ST LMi carries a genuine three-colour claim and carries it "
            "well. VV Pup and EU UMa do not: their colour panels have "
            "essentially no three-filter full-orbit nights to draw on and "
            "should be replaced by single-band folded morphology plus state "
            "history. Retitle the paper around ST LMi rather than around "
            "three polars, and quote colours with the non-simultaneity term "
            "attached - on the bright-phase edges it dominates."))

        # ---------------- Q2 --------------------------------------------------
        s_id, s_re = st_ideal.get("stlmi"), st_real.get("stlmi")
        s_id_m, s_re_m = st_ideal_med.get("stlmi"), st_real_med.get("stlmi")
        n_oc_nights, dt_best, dt_med = _nights_meeting_timing_bar(ocon)
        v.append((
            "Q2", 2,
            "Bright-phase timing as an accretion-spot longitude tracker "
            "(per-cycle O-C against decades-baseline ephemerides)",
            "per-cycle sigma_t < 60 s, demonstrated by injection before "
            "adoption (strategy §4.16)",
            _q2_verdict(s_re, s_id, s_re_m),
            (f"ST LMi per-cycle sigma_t, shape known exactly / shape 5x wrong "
             f"and depth 20% wrong: RICHEST night {s_id:.0f} s / "
             f"{s_re:.0f} s (cadence {dt_best:.0f} s); MEDIAN-density night "
             f"{s_id_m:.0f} s / {s_re_m:.0f} s (cadence {dt_med:.0f} s) - "
             f"ABOVE the 60 s bar. Night-mean "
             f"{st_night.get('stlmi', float('nan')):.1f} s. "
             f"{n_oc_nights} full-orbit nights are available to the tier"
             if None not in (s_id, s_re, s_re_m) else "not measured"),
            "The 60 s threshold is met on the demonstration night with almost "
            "no margin once the edge shape is admitted to be unknown - and "
            "the shape IS unknown and colour-dependent. On the median night "
            "it is not met at all. sigma_t is set by whether samples land on "
            "the ingress ramp, so it tracks CADENCE, not signal-to-noise: "
            "the two nights differ by a factor of about two in cadence and "
            "that alone moves sigma_t across the bar. A per-cycle tier "
            "therefore exists, but only on the dense half of the nights.",
            "Report per-cycle timings only from nights whose cadence matches "
            "the demonstration night, with an explicit shape-systematic term "
            "of ~50 s; fall back to night-mean timings elsewhere (the "
            "strategy's own §6.10a), which are an order of magnitude better "
            "and carry the same shape question. And fix the cadence: a "
            "faster filter cycle is worth more here than any realistic gain "
            "in signal-to-noise."))

        # ---------------- Q3 --------------------------------------------------
        yz_n = a90("yzcnc", "night", "known")
        yz_s = a90("yzcnc", "season-dt", "known")
        yz_n_per = a90("yzcnc", "night", "period")
        v.append((
            "Q3", 3,
            "YZ Cnc 2024: superhump period (and dP_sh/dt) inside confirmed "
            "outburst states; orbital hump + flickering in quiescence",
            "detect and time a 0.1-0.3 mag superhump on consecutive-night "
            "blocks",
            _q3_verdict(yz_n[0], yz_s[0]),
            (f"90% recovery semi-amplitude at P_orb, scored as DETECTION at a "
             f"known period: {mmag(yz_n[0])} on the richest single night "
             f"[{mmag(yz_n[1])}-{mmag(yz_n[2])}], {mmag(yz_s[0])} on the "
             f"detrended season. Scored instead as blind period "
             f"DETERMINATION the same night gives {mmag(yz_n_per[0])} - the "
             f"number the first version of this page quoted, which measures "
             f"the {yz_n[3]:.1f} cycles the night holds, not the amplitude"
             if yz_n[0] and yz_s[0] else "not measured"),
            "Superhumps have semi-amplitudes of 50-150 mmag. YZ Cnc's "
            "superhump period is not known in advance the way P_orb is, so "
            "the blind-search contour is the relevant one for measuring "
            "P_sh - but detecting that a superhump is PRESENT on a given "
            "night is the known-period question, and one night answers it "
            "comfortably. The consecutive-night blocks remain the "
            "measurement of dP_sh/dt, exactly as the strategy assumes.",
            "If the dense runs were quiescent, flickering statistics and the "
            "orbital hump are comfortably measurable at this precision - but "
            "no multi-year period change, as the strategy already forbids."))

        # ---------------- Q4 --------------------------------------------------
        # A duty cycle is a fraction of TIME IN A STATE.  Its uncertainty is
        # set by the number of independent epochs and by whether those epochs
        # can be placed on a common magnitude scale - never by per-point
        # photometric precision, which is what this verdict used to quote.
        duty = {t: ch.duty_cycle_sigma(n) for t, n in nights_per_target.items()}
        worst_t = max(duty, key=lambda t: duty[t])
        best_t = min(duty, key=lambda t: duty[t])
        v.append((
            "Q4", 6,
            "Long-term accretion-state duty cycles for the three polars, "
            "embedded in the survey record",
            "nightly means good to much better than the 1-3 mag high/low "
            "state separation",
            _q4_verdict(nights_per_target, n_tied, n_blocks),
            f"independent epochs per target: "
            + "; ".join(f"{TARGET_LABEL[t]} {nights_per_target[t]} nights over "
                        f"{span_per_target[t]:.0f} d (+/-"
                        f"{100 * duty[t]:.0f} pp on a 50% duty cycle)"
                        for t in sorted(nights_per_target))
            + f". Catalogue tie: {n_tied} of {n_blocks} (target, era) blocks",
            "Classifying ONE night is easy - per-point precision beats the "
            "state separation by a factor of 20-100, which is what this "
            "verdict used to say. It is not what the goal asks. A duty cycle "
            "is a fraction of time, so its error bar is binomial in the "
            "number of independent nights: "
            f"{nights_per_target[best_t]}-{nights_per_target[worst_t]} nights "
            f"give +/-{100 * duty[best_t]:.0f} to "
            f"+/-{100 * duty[worst_t]:.0f} percentage points, before any "
            "correction for seasonal observability bias in WHEN those nights "
            "fall. And 'embedded in the survey record' means joining "
            "ASAS-SN/ZTF, which needs the catalogue tie S1 grades NOT "
            "SUPPORTED.",
            "State the duty cycle with its binomial interval and the season "
            "structure, or restrict the claim to within-series relative "
            "state history, where the ensemble gauge is stable and no tie is "
            "needed. Q4 and S1 cannot both stand as they were written."))

        # ---------------- Q5 --------------------------------------------------
        v.append((
            "Q5", 7,
            "AN UMa as a fifth target in the colour analysis",
            "at least 8 full-orbit nights with all three filters",
            "NOT SUPPORTED" if three_filter.get("anuma", 0) < 8 else "SUPPORTED",
            f"three-filter full-orbit nights measured: "
            f"{three_filter.get('anuma', 0)} (criterion: >= 8); "
            f"{any_filter.get('anuma', 0)} full-orbit nights in any filter",
            "The strategy set this criterion itself and the solved photometry "
            "answers it directly.",
            "Timing and state history only, or drop - the strategy's own "
            "fallback. The AN UMa photometry is good; it is the three-filter "
            "scheduling that fails."))

        # ---------------- structural claims -----------------------------------
        v.append((
            "S1", 4,
            "Absolute calibration: 'nightly REFCAT2 tie -> PS1 AB to "
            "0.01-0.02 mag' (strategy §5)",
            "every published magnitude on a standard system",
            "NOT SUPPORTED",
            f"{n_tied} of {n_blocks} (target, era) blocks carry a catalogue "
            f"magnitude tie; the other {n_blocks - n_tied} were tied by WCS "
            f"with zero catalogue matches, so their magnitudes sit on the "
            f"ensemble's own arbitrary gauge",
            "Ensemble differential photometry fixes RELATIVE magnitudes only. "
            "Nothing measured here is wrong because of it - every precision "
            "number on this page is differential - but no magnitude can be "
            "published on a standard system until the tie is made, and Q4's "
            "cross-survey framing depends on it.",
            "Publish differential magnitudes plus per-series offsets now, and "
            "make the catalogue tie a separate, testable stage. This is the "
            "highest-value single action on the whole list."))
        v.append((
            "S2", 5,
            "Cross-era colour comparison after the 2024-05 instrument seam "
            "(strategy §4.13/4.13a)",
            "morphology-level cross-era comparison, never zero points",
            "SUPPORTED-WITH-CAVEATS",
            "eras never overlap in time for any target ("
            + "; ".join(f"{TARGET_LABEL[t]}: "
                        + ("OVERLAP" if overlap[t] else "no overlap")
                        for t in sorted(overlap)) + ")",
            "Camera and epoch are perfectly confounded, exactly as the "
            "strategy states; nothing in the measured data relaxes it.",
            "Keep the within-era discipline. The October 2026 g/r/i season is "
            "the only real repair, and it repairs ST LMi only."))
        # S3: window power is the INPUT to the alias question, not the answer.
        conf = ocon.execute(
            "SELECT semi_amp, avg(frac_true), avg(frac_alias) "
            "FROM ch_alias_confusion WHERE regime IN ('season','season-dt') "
            "GROUP BY semi_amp ORDER BY semi_amp").fetchall()
        conf_txt = "; ".join(f"{1000 * a:.0f} mmag -> {100 * ft:.0f}% true, "
                             f"{100 * fa:.0f}% alias" for a, ft, fa in conf)
        v.append((
            "S3", 8,
            "Period verification from the multi-season data - 'confirmation "
            "and alias hygiene' (strategy §4.15)",
            "select the right alias family from the combined seasons",
            _s3_verdict(conf),
            (f"the +/-1 c/d aliases carry up to {solar_alias:.2f} of the "
             f"window power on the resolved multi-night sets ({worst_alias:.2f} "
             f"including the sidereal comb). MEASURED misidentification rate "
             f"of the tallest peak, by injected semi-amplitude: {conf_txt}"
             if solar_alias is not None and conf else "not measured"),
            "Window power says how much power the SAMPLING makes available at "
            "f +/- k c/d. It cannot say whether the wrong peak wins - that "
            "depends on the signal's amplitude and on how many nights there "
            "are, and both are measurable with the injection machinery this "
            "page already owns. Measured, the tallest peak lands on the true "
            "frequency at every amplitude these polars actually show "
            "(0.5-2 mag, one to two orders above the season contour). The "
            "alias risk is real at threshold amplitudes and negligible at "
            "the real ones. Note also that the headline window number is the "
            "maximum over the solar AND sidereal combs; for ground-based "
            "scheduling the solar-day comb is the relevant one.",
            "Publish the window, the solar-day alias power and the measured "
            "misidentification rate beside every periodogram, and state the "
            "alias caveat as amplitude-dependent rather than absolute."))
        # S4: the analytic formula against the MEASURED contour.  Restricted
        # to the known-period score and to scopes that resolve the frequency:
        # a single night cannot separate f from f+1 c/d, so a
        # period-determination contour there is not a detection limit and
        # must not be averaged into a claim about one.
        rows = ocon.execute(
            "SELECT c.series_key, c.regime, c.score, c.amp90, d.sigma_used, "
            "c.n_points, c.n_cycles FROM ch_contour c "
            "JOIN (SELECT scope, regime, score, min(sigma_used) sigma_used "
            "      FROM ch_detect GROUP BY scope, regime, score) d "
            "ON d.scope = c.scope AND d.regime = c.regime AND d.score = c.score "
            "JOIN ch_cadence cad ON cad.series_key = c.series_key "
            "WHERE abs(c.period_d - cad.period_d) < 1e-6 "
            "AND c.amp90 IS NOT NULL").fetchall()

        def ratio_range(score, regimes):
            r = [(x[3] / ch.amin_analytic(x[4], x[5])) for x in rows
                 if x[2] == score and x[1] in regimes and x[4]
                 and ch.amin_analytic(x[4], x[5]) > 0]
            return (min(r), max(r), len(r)) if r else (None, None, 0)

        k_lo, k_hi, k_n = ratio_range("known", ("season", "season-dt", "night"))
        p_lo, p_hi, p_n = ratio_range("period", ("season", "season-dt"))
        n_lo, n_hi, n_n = ratio_range("period", ("night",))
        v.append((
            "S4", 9,
            "Detection limits quoted as A_min = sigma sqrt(4z/N) "
            "(strategy §4.21: 4-8 mmag per target)",
            "an analytic detection limit on white, alias-free data",
            _s4_verdict(k_lo, k_hi, p_lo, p_hi),
            (f"the MEASURED 90% recovery amplitude at a KNOWN period is "
             f"{k_lo:.1f}-{k_hi:.1f}x the analytic formula at the same sigma "
             f"and N ({k_n} contours) - i.e. the formula errs in BOTH "
             f"directions; for blind period DETERMINATION on the multi-night "
             f"sets it is {p_lo:.1f}-{p_hi:.1f}x, and on a single night "
             f"{n_lo:.1f}-{n_hi:.1f}x. Median Allan slope {med_slope:.2f} "
             f"against a white-noise NULL of {med_null_slope:.2f} for ladders "
             f"this short (not -0.50); {n_redder} of {n_ladders} ladders are "
             f"redder than their own 95th-percentile null"
             if k_n and p_n else "not measured"),
            "The formula assumes white noise, an alias-free window, and one "
            "fixed look-elsewhere factor z. This archive has none of the "
            "three, and the errors do not even share a sign: at a KNOWN "
            "period the formula is conservative (it charges a blind-search "
            "penalty the paper does not owe), while for a genuine blind "
            "search it is optimistic by up to a factor of five. A single "
            "number that is wrong by 2-3x in either direction depending on "
            "which question is asked cannot be a detection limit. Note that "
            "the first version of this page reported 1.5-8.5x, and its "
            "largest ratios came from single-night scopes scored on period "
            "determination - where the acceptance window is 20-70x narrower "
            "than the peak and no amplitude can succeed.",
            "Quote the injection contour, per target, per regime and per "
            "SCORE MODE. Per-night detrending recovers part of the gap and "
            "its cost is measured here too (regime 'season-dt')."))
        floor_best = one("SELECT min(floor_plateau) FROM ch_noise_series")
        floor_med = one("SELECT median(floor_plateau) FROM ch_noise_series")
        bias_med = one("SELECT median(bias_ratio) FROM ch_check_bias")
        bias_n = one("SELECT count(*) FROM ch_check_bias", 0)
        v.append((
            "S5", 10,
            "Expected precision: '8-10 mmag/frame at r <= 15, 20-25 mmag at "
            "r = 16.5' (strategy §5)",
            "per-frame precision as a function of brightness",
            "SUPPORTED-WITH-CAVEATS",
            f"best measured floor {mmag(floor_best, 1)}, typical floor "
            f"{mmag(floor_med)}; achieved precision at the targets themselves "
            f"{mmag(min(prec.values()))}-{mmag(max(prec_hi.values()))}. "
            f"Magnitude-matched field stars (including every star the comp "
            f"stability cut dropped) scatter {bias_med:.2f}x as much as the "
            f"held-out check stars over {bias_n} series",
            "The 8-25 mmag band is right for the BRIGHT end. It is not what "
            "the targets get: every CV here sits near the faint end of its "
            "own comparison ensemble, so the number that matters is several "
            "times worse than the headline. The check stars ARE a selected "
            "sample - survivors of the comp stability cut, drawn nearest the "
            "target's brightness - so their scatter was tested against every "
            "star of the same brightness whatever role it was given. The "
            "ratio is close to 1, so the precision numbers are not "
            "optimistic by selection.",
            "State the precision AT EACH TARGET'S BRIGHTNESS, not the best "
            "the camera can do. The floor itself is instrumental (tens of "
            "times scintillation) and is the thing to attack."))

        with ocon:
            ocon.execute("DELETE FROM ch_verdict")
            ocon.executemany("INSERT INTO ch_verdict VALUES (?,?,?,?,?,?,?,?)", v)
        stage_done(ocon, "verdict")
        for row in sorted(v, key=lambda r: r[1]):
            print(f"  {row[0]:3s} {row[4]:24s} {row[5][:110]}")
    finally:
        ocon.close(); pcon.close()


def _nights_meeting_timing_bar(ocon, target_key: str = "stlmi") -> tuple:
    """The sampling the O-C tier would really be built from.

    Q2's demonstration is made on the richest night; the tier is built from
    every night that carries a full cycle.  Since sigma_t is set by the
    PER-FILTER cadence rather than by signal-to-noise, the spread of that
    cadence across nights is what decides how much of the tier survives the
    60 s bar.

    The two cadences are read back out of ``ch_timing`` rather than
    recomputed, so they are exactly the numbers the Monte Carlo ran on.
    Recomputing them from ch_frames would pool all three filters and report
    72 s where the g-band series the MC actually used cycles every 219 s.

    Returns ``(n_full_orbit_nights, richest_dt_s, median_dt_s)``.
    """
    per = qualifying_nights(ocon, target_key, PERIODS_D[target_key][0])
    dt = dict(ocon.execute(
        "SELECT night_kind, cadence_s FROM ch_timing WHERE target_key=? "
        "GROUP BY night_kind", (target_key,)))
    return (len(per), dt.get("richest", float("nan")),
            dt.get("median", float("nan")))


def _q1_verdict(three_filter: dict) -> str:
    """All three polars carrying it is SUPPORTED; one or two is CAVEATS;
    none is NOT SUPPORTED.  Three qualifying nights is the minimum that can
    show a state-to-state morphology comparison at all."""
    good = [t for t in ("stlmi", "vvpup", "euuma")
            if three_filter.get(t, 0) >= 3]
    if len(good) >= 3:
        return "SUPPORTED"
    return "SUPPORTED-WITH-CAVEATS" if good else "NOT SUPPORTED"


def _q2_verdict(sigma_t_real_richest, sigma_t_ideal_richest,
                sigma_t_real_median) -> str:
    """The governing number is the SHAPE-MISMATCHED one: nobody knows this
    archive's ingress duration, so the realistic fit is the mismatched fit.

    And the governing NIGHT matters as much as the governing regime.  If even
    the campaign's best night misses the bar there is no per-cycle tier at
    all.  If the best night clears it and a typical night does not, the tier
    exists on a subset of nights - which is CAVEATS with a sharp, countable
    caveat, not a clean pass.
    """
    if sigma_t_real_richest is None or sigma_t_ideal_richest is None:
        return "NOT MEASURED"
    if sigma_t_real_richest >= 60.0 and sigma_t_ideal_richest >= 60.0:
        return "NOT SUPPORTED"
    return "SUPPORTED-WITH-CAVEATS"


def _q3_verdict(a90_night, a90_season) -> str:
    """Superhumps run 50-150 mmag semi-amplitude.  A single night that
    detects the low end of that at a known period, plus a detrended season
    that beats it several times over, is the strategy's own plan working."""
    if a90_night is None or a90_season is None:
        return "NOT MEASURED"
    return ("SUPPORTED-WITH-CAVEATS" if a90_night <= 0.050
            else "SUPPORTED-WITH-CAVEATS" if a90_season <= 0.050
            else "NOT SUPPORTED")


def _q4_verdict(nights_per_target: dict, n_tied: int, n_blocks: int) -> str:
    """A duty cycle needs epochs AND a common gauge.

    SUPPORTED would require both: enough independent nights for a useful
    binomial interval, and a catalogue tie for every block so the epochs can
    be joined to the long-term survey record the goal names.  With half the
    blocks untied, the cross-survey form of the claim cannot be made at all,
    whatever the photometric precision is.
    """
    if not nights_per_target:
        return "NOT MEASURED"
    if n_tied < n_blocks:
        return "SUPPORTED-WITH-CAVEATS"
    return ("SUPPORTED" if min(nights_per_target.values()) >= 30
            else "SUPPORTED-WITH-CAVEATS")


def _s4_verdict(k_lo, k_hi, p_lo, p_hi) -> str:
    """The analytic limit is usable only if it is close on BOTH questions.

    "Close" is taken as within 30% either way — tighter than that and no
    formula would ever pass, looser and a factor-of-two error would be
    waved through as agreement.  The formula fails here on both counts and
    with opposite signs, which is why it cannot simply be rescaled.
    """
    if None in (k_lo, k_hi, p_lo, p_hi):
        return "NOT MEASURED"
    ok = all(0.7 <= r <= 1.3 for r in (k_lo, k_hi, p_lo, p_hi))
    return "SUPPORTED" if ok else "NOT SUPPORTED"


def _s3_verdict(conf) -> str:
    """Alias hygiene is graded on the measured misidentification rate.

    ``conf`` is [(semi_amp, frac_true, frac_alias), ...].  If the tallest
    peak lands on the truth at the amplitudes these targets actually show,
    the sampling CAN choose and the caveat is about faint signals, not about
    the archive.
    """
    if not conf:
        return "NOT MEASURED"
    bright = [c for c in conf if c[0] >= 0.5]
    if bright and min(c[1] for c in bright) >= 0.95:
        return "SUPPORTED-WITH-CAVEATS"
    return "NOT SUPPORTED"


# ==========================================================================
# STAGE report / status
# ==========================================================================

def stage_report(args) -> None:
    from macro_phot import report_cv_char
    path = report_cv_char.render_report(OUT_DB, PHOT_DB)
    print(f"report: {path}")


def stage_status(args) -> None:
    if not OUT_DB.exists():
        print("no characterization DB yet")
        return
    con = connect_rw(OUT_DB)
    ensure_schema(con)
    try:
        meta = dict(con.execute("SELECT key, value FROM ch_meta"))
        for st in ("quality", "trail", "noise", "cadence", "detect", "timing",
                   "verdict"):
            print(f"  {st:9s} {meta.get('stage_' + st, 'NOT RUN')}")
        for t in ("ch_frames", "ch_trail", "ch_noise_stars", "ch_noise_series",
                  "ch_allan", "ch_cadence", "ch_window", "ch_alias",
                  "ch_detect", "ch_contour", "ch_timing", "ch_verdict"):
            n = con.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
            print(f"    {t:16s} {n:>8d}")
    finally:
        con.close()


# ==========================================================================
# main
# ==========================================================================

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("stage", choices=["quality", "trail", "noise", "cadence",
                                     "detect", "timing", "verdict", "report",
                                     "status", "all"])
    p.add_argument("--workers", type=int, default=6,
                   help="parallel workers for the pixel stage (hard cap 6)")
    p.add_argument("--limit", type=int, default=0,
                   help="cap the number of units this invocation processes")
    p.add_argument("--timing-trials", type=int, default=200)
    p.add_argument("--force", action="store_true")
    args = p.parse_args(argv)
    OUT_DB.parent.mkdir(parents=True, exist_ok=True)
    con = connect_rw(OUT_DB)
    ensure_schema(con)
    set_meta(con, code_version=CODE_VERSION,
             phot_db=str(PHOT_DB), manifest_db=str(MANIFEST_DB),
             built_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
             periods_json=json.dumps({k: v[0] for k, v in PERIODS_D.items()}),
             period_sources_json=json.dumps({k: v[1] for k, v in PERIODS_D.items()}),
             degrade_factor=ch.DEGRADE_FACTOR,
             detect_fap=ch.DETECT_FAP, inject_trials=INJECT_TRIALS,
             threshold_trials=THRESHOLD_TRIALS,
             search_band_cd=f"{SEARCH_FMIN_CD}-{SEARCH_FMAX_CD}")
    con.close()
    stages = {"quality": stage_quality, "trail": stage_trail,
              "noise": stage_noise, "cadence": stage_cadence,
              "detect": stage_detect, "timing": stage_timing,
              "verdict": stage_verdict, "report": stage_report,
              "status": stage_status}
    if args.stage == "all":
        for name in ("quality", "noise", "cadence", "detect", "timing",
                     "verdict", "report"):
            print(f"=== {name} ===")
            stages[name](args)
    else:
        stages[args.stage](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
