#!/usr/bin/env python
"""CV-S9 — Phase 3: the time-series analysis itself.

WHAT THIS SCRIPT DOES, AND WHY EACH PIECE IS HERE
--------------------------------------------------
Phase 1 measured the light curves, Phase 2 decided what may be published
from them.  Phase 3 asks the six questions the project exists for, and each
stage below is a place where a confident wrong answer is easy to produce:

``ephem``       Fetches the published ephemerides from the AAVSO Variable
                Star Index and stores them WITH their provenance and with
                the uncertainty VSX does not publish converted into an
                explicit quoted-precision floor.  Everything downstream
                that says "agrees with the literature" is anchored here, so
                the anchor is a stored payload rather than a number someone
                typed.

``periods``     Per (target, era, filter): a generalised Lomb-Scargle with
                one free constant per night, a phase-dispersion
                minimisation, and THE SPECTRAL WINDOW, stored beside them.
                The window is not decoration.  On every resolved multi-night
                set in this archive the +/-1 c/d sidelobes carry 0.54-0.97
                of the window power, so the periodogram cannot select the
                period and this stage records HOW the family member was
                chosen instead of pretending the question did not arise.

``sigmat``      Injects bright-phase templates at the REAL timestamps of ST
                LMi's densest night with the REAL error model and recovers
                them, over a grid of template-shape and depth errors.  This
                is the measurement that decides whether per-cycle timing is
                publishable against the 60 s threshold, and it exists
                because the formal delta-chi2 = 1 error bar from an edge fit
                on this data comes out at 0-4 s, which is a fiction.

``edges``       Measures bright-phase edge epochs per cycle PER BAND for the
                polars.  Colour-dependent edge timing is a cyclotron result,
                not a nuisance: the emission region's optical depth is
                wavelength dependent, so the band in which the column
                disappears behind the white dwarf's limb first is physics.
                Reported as inter-band differences with uncertainties.

``oc``          Builds the O-C diagram against the published ephemeris and
                analyses the cycle count EXPLICITLY.  Between VSX's epoch
                and our nights there are tens of thousands of cycles, and
                the count is unique only if the published period is good
                enough.  When it is not, this stage says so; an O-C on a
                wrong cycle count is a fabricated result, not a noisy one.

``states``      Classifies each night's accretion state from the calibrated
                magnitudes with a threshold derived from the observed
                bimodality, and computes duty cycles using the Phase-2 upper
                limits so the statistic is not conditioned on the target
                having been bright enough to detect.

``detrend``     Implements joint GP + signal fitting (celerite2, checked
                against a dense pure-numpy reference) and DEMONSTRATES on a
                real series why detrend-then-search is the wrong order, by
                injecting a known signal and recovering it both ways.

All the arithmetic lives in ``macro_phot.phase3`` and is unit-tested in
``pipeline/tests/test_phase3.py``.  This script is I/O, staging, parallelism
and bookkeeping.

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_phase3.py ephem
    $P pipeline/scripts/run_cv_phase3.py periods --workers 6
    $P pipeline/scripts/run_cv_phase3.py sigmat
    $P pipeline/scripts/run_cv_phase3.py edges
    $P pipeline/scripts/run_cv_phase3.py oc
    $P pipeline/scripts/run_cv_phase3.py states
    $P pipeline/scripts/run_cv_phase3.py detrend
    $P pipeline/scripts/run_cv_phase3.py report
    $P pipeline/scripts/run_cv_phase3.py status

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite)
--------------------------------------------------------------
``p3_ephemeris``   the published ephemerides, their source and their
                   quoted-precision uncertainty floor.
``p3_period``      per series: survey peak, orbital-band peak, PDM trough,
                   alias fractions, how the family was chosen, the recovered
                   period with its uncertainty, and the agreement verdict.
``p3_pgram``       decimated periodogram, PDM and SPECTRAL WINDOW traces, so
                   every figure on the page is drawn from the database.
``p3_sigmat``      the sigma_t injection grid and its contour verdict.
``p3_edge``        per (series, cycle): the fitted bright-phase edge epoch.
``p3_band_pair``   inter-band edge-time differences, per night and pooled.
``p3_oc``          per-cycle edge residuals: the INPUTS to the O-C, not
                   publishable epochs -- the injection test does not license
                   a single cycle's edge as one.
``p3_oc_night``    the PUBLISHED timing epochs: one per night per band, the
                   mean of that night's accepted per-cycle edges, with the
                   error bar the injection budget licenses.
``p3_cycle_count`` the cycle-count ambiguity analysis, per target, plus the
                   per-night O-C summary and its refitted ephemeris.
``p3_state_night`` per (series, night): the state classification.
``p3_state_series``thresholds, separability and duty cycles.
``p3_detrend``     the detrend-then-search versus joint-fit demonstration.
``p3_meta``        build stamps and every constant this run used.

CONCURRENCY
-----------
An astrometry batch and a catalogue-tie re-run may be writing this archive
at the same time.  Every connection sets ``busy_timeout = 300000``, worker
count is hard-capped at 6, transactions are short and per-stage, and the
manifest is opened READ-ONLY and never rebuilt.
"""

from __future__ import annotations

import argparse
import json
import math
import sqlite3
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import phase3 as p3                            # noqa: E402
# Kaplan-Meier is Phase 2's, not re-implemented here: the limits these
# nights carry were produced by that module and must be consumed by the
# same estimator that produced its own state statistics, or the two pages
# would disagree about the same epochs.
from macro_phot.phase2 import km_median as p2_km_median        # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
DEFAULT_MANIFEST = (REPO_ROOT / "products" / "manifest"
                    / "rlmt-manifest.sqlite")
CACHE_DIR = REPO_ROOT / "products" / "external" / "vsx"

#: Stamped into ``p3_meta`` and read by the provenance graph.  Bump it when
#: the arithmetic changes, not when a comment does.
PHASE3_CODE_VERSION = "CV-S9 v1.0 (2026-08-20, Phase-3 time series)"

#: Hard cap on worker processes.  This machine is also running an S1
#: astrometry batch and a catalogue-tie re-run; a stage that saturates the
#: disk queue slows both of them and finishes no sooner itself.
MAX_WORKERS = 6

BUSY_TIMEOUT_MS = 300_000
NET_TIMEOUT_S = 45
USER_AGENT = "MACRO-pipeline/1.0 (RLMT archive; research use)"

#: The five CV targets, and the VSX identifier each one is looked up by.
#: The keys are the ``target_key`` values used throughout the products.
VSX_IDENT = {"stlmi": "ST LMi", "anuma": "AN UMa", "vvpup": "VV Pup",
             "euuma": "EU UMa", "yzcnc": "YZ Cnc"}

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}

#: Which targets are POLARS (AM Her stars with a bright phase to time) and
#: which is the dwarf nova.  Read from the VSX variability type when one is
#: available, and used only to decide which stages apply to which target.
POLARS = ("stlmi", "anuma", "vvpup", "euuma")

#: Series with fewer target points than this are not searched.  Twenty is
#: two per PDM bin, below which theta is noise.
MIN_POINTS_FOR_SEARCH = 20

#: Number of points a decimated periodogram trace keeps for the figures.
TRACE_POINTS = 1400

#: The sigma_t grid: template shape errors (multiples of the true edge
#: width) crossed with fractional depth errors.  1.0 / 0.0 is the edge shape
#: known exactly; 5.0 / 0.2 is the strategy's stated worst case.
SHAPE_ERRORS = (1.0, 2.0, 3.0, 5.0)
DEPTH_ERRORS = (0.0, 0.1, 0.2, 0.5)

#: Monte Carlo realizations per sigma_t grid cell.
N_SIGMAT_REAL = 300

#: INJECTED edge widths, as multiples of the width fitted from the fold.
#:
#: This axis exists because the fitted width turned out to be the single
#: most consequential and least known input to the whole test.  A first
#: version of the shape estimator was algebraically pinned at 547 s and
#: returned "CONDITIONAL, sigma_t = 20 s with the shape known"; measuring
#: the width properly gave 29-48 s AT THE FOLD'S SAMPLING FLOOR — i.e. the
#: edge is UNRESOLVED, the number is an upper bound — and the verdict
#: flipped to "NOT PUBLISHABLE, 81-131 s".  A sharper edge is harder to
#: time, not easier, because fewer points land on the ramp and the epoch
#: is bracketed by the 219 s cadence instead of interpolated within it.
#: Since the data cannot say which it is, the grid is run at all three and
#: the verdict has to survive the range rather than a chosen point.
INJECT_WIDTH_FACTORS = (0.2, 1.0, 2.0)

#: The night the sigma_t test runs on: ST LMi's densest, 153 + 141 + 142
#: frames in g, r and i over 9.37 h — 4.94 orbits with each filter sampled
#: every 219 s.  Stored as the LOCAL night label used by ``cv_frames``;
#: the UTC night is one day later, which is the "2025-02-28" the plan names.
SIGMAT_NIGHT = "2025-02-27"
SIGMAT_TARGET = "stlmi"

#: Detrend demonstration: running-median window widths, in units of the
#: orbital period.
DETREND_WINDOWS_P = (0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0)
DETREND_AMPLITUDE = 0.10        # mag, semi-amplitude of the injected signal
DETREND_REALS = 12


# ===========================================================================
# Database plumbing
# ===========================================================================
def connect(path: Path, read_only: bool = False) -> sqlite3.Connection:
    """One connection, always with the long busy timeout.

    ``busy_timeout`` rather than a lock file: two other workflows are
    writing this archive today, and SQLite's own waiter is the only one all
    three processes agree about.  Five minutes is longer than any single
    transaction any of the three takes.
    """
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    return con


def git_commit() -> str:
    """Short commit of the tree that produced this run, ``-dirty`` when the
    working tree is not clean.  A product stamped with a clean commit that
    was actually built from edited files is a false reproducibility claim,
    so the suffix is not optional."""
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


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def ensure_tables(con: sqlite3.Connection) -> None:
    """Create every table this stage owns.  Idempotent by construction."""
    con.executescript("""
    CREATE TABLE IF NOT EXISTS p3_meta (key TEXT PRIMARY KEY, value TEXT);

    CREATE TABLE IF NOT EXISTS p3_ephemeris (
        target_key TEXT PRIMARY KEY, name TEXT, var_type TEXT,
        period_d REAL, period_str TEXT, period_sigma_d REAL,
        sigma_basis TEXT, epoch_bjd REAL, epoch_str TEXT,
        source TEXT, fetched_utc TEXT, payload_path TEXT, note TEXT);

    CREATE TABLE IF NOT EXISTS p3_period (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, n_points INTEGER, n_blocks INTEGER, baseline_d REAL,
        median_cadence_s REAL, phase_coverage REAL,
        f_survey_cd REAL, p_survey_pow REAL, survey_is_orbital INTEGER,
        survey_class TEXT,
        f_ls_cd REAL, p_ls_pow REAL, f_pdm_cd REAL, pdm_theta REAL,
        ls_minus_pub_cd REAL, pdm_minus_pub_cd REAL,
        alias_frac_p1 REAL, alias_frac_m1 REAL, alias_frac_max REAL,
        peak_halfwidth_cd REAL, amplitude_mag REAL, resid_rms_mag REAL,
        period_d REAL, sigma_period_d REAL, sigma_formal_d REAL,
        sigma_boot_d REAL, sigma_resolution_d REAL, n_boot INTEGER,
        sigma_basis TEXT, refine_halfwidth_cd REAL,
        frac_precision REAL, constraint_class TEXT,
        published_d REAL, published_sigma_d REAL,
        deviation_sigma REAL, agrees INTEGER, agree_note TEXT,
        family_code TEXT, family_note TEXT, harmonic_note TEXT,
        detected INTEGER, status TEXT, note TEXT);

    CREATE TABLE IF NOT EXISTS p3_pgram (
        series_key TEXT, panel TEXT, kind TEXT, freq_cd REAL, value REAL);

    CREATE TABLE IF NOT EXISTS p3_sigmat (
        series_key TEXT, night TEXT, inject_factor REAL,
        inject_width_s REAL, shape_error REAL, depth_error REAL,
        sigma_t_s REAL, bias_s REAL, total_error_s REAL, rms_s REAL,
        p95_abs_s REAL, n_ok INTEGER, n_try INTEGER,
        recovered_fraction REAL, passes INTEGER,
        PRIMARY KEY (series_key, inject_factor, shape_error, depth_error));

    CREATE TABLE IF NOT EXISTS p3_sigmat_input (
        series_key TEXT PRIMARY KEY, night TEXT, n_points INTEGER,
        median_cadence_s REAL, median_err_mag REAL, chi2_inflation REAL,
        used_err_mag REAL, depth_mag REAL, edge_width_s REAL,
        edge_width_floor_s REAL,
        bright_width_phase REAL, n_cycles REAL, note TEXT);

    CREATE TABLE IF NOT EXISTS p3_edge (
        series_key TEXT, cycle INTEGER, target_key TEXT, era_id INTEGER,
        filter TEXT, night TEXT, t_edge_bjd REAL, sigma_t_s REAL,
        sigma_t_mc_s REAL, width_s REAL, level_bright REAL,
        level_faint REAL, depth_mag REAL, depth_snr REAL, chi2nu REAL,
        n_points INTEGER, bracket_s REAL, phase REAL,
        accepted INTEGER, reason TEXT,
        PRIMARY KEY (series_key, cycle));

    CREATE TABLE IF NOT EXISTS p3_band_pair (
        target_key TEXT, era_id INTEGER, night TEXT, band_a TEXT,
        band_b TEXT, n_cycles INTEGER, delta_s REAL, sigma_s REAL,
        chi2nu REAL, significant INTEGER, note TEXT,
        PRIMARY KEY (target_key, era_id, night, band_a, band_b));

    CREATE TABLE IF NOT EXISTS p3_oc (
        series_key TEXT, cycle INTEGER, target_key TEXT, filter TEXT,
        night TEXT, t_edge_bjd REAL, sigma_t_s REAL, oc_s REAL,
        oc_sigma_s REAL, phase_offset REAL, count_unique INTEGER,
        PRIMARY KEY (series_key, cycle));

    -- The PUBLISHED timing epochs.  One row per (target, night, band):
    -- the mean of that night's accepted per-cycle edges in that band.
    -- CV-S5's injection test demonstrated that a single cycle's edge does
    -- not reach the 60 s threshold with an error bar of its OWN, so no
    -- per-cycle edge is published carrying its own per-cycle error; p3_oc
    -- holds those as the inputs and this table holds what the paper is
    -- allowed to plot and fit, each row carrying the injection-demonstrated
    -- budget instead.  On nights where only one cycle was timed, n_cycles
    -- is 1 and the epoch IS that cycle -- with the budget error, not the
    -- fit's own.  Bands are NOT pooled: the band-pair stage finds no
    -- significant band-to-band offset, but a wavelength-dependent edge
    -- phase is expected on physical grounds and averaging within a band is
    -- the conservative choice against an effect this data set could not
    -- have detected (see p3_band_pair, 0 of 32 pairs significant).
    -- ``oc_sigma_edge_s`` is the same epoch's error propagated from the
    -- edge fits' OWN Monte-Carlo sigma_t instead of the budget: the
    -- robustness check on the transfer of a one-night budget across bands
    -- and instrument eras.
    CREATE TABLE IF NOT EXISTS p3_oc_night (
        target_key TEXT, night TEXT, filter TEXT, series_key TEXT,
        era_id INTEGER, n_cycles INTEGER, cycle_mean REAL,
        cycle_lo INTEGER, cycle_hi INTEGER, t_mean_bjd REAL,
        oc_s REAL, oc_sigma_s REAL, sigma_random_s REAL,
        sigma_floor_s REAL, within_night_rms_s REAL, meets_threshold INTEGER,
        oc_sigma_edge_s REAL, budget_band TEXT,
        PRIMARY KEY (target_key, night, filter));

    CREATE TABLE IF NOT EXISTS p3_cycle_count (
        target_key TEXT PRIMARY KEY, epoch_bjd REAL, period_d REAL,
        sigma_period_d REAL, sigma_basis TEXT,
        t_first_bjd REAL, t_last_bjd REAL, elapsed_d REAL,
        n_cycles_first REAL, n_cycles_last REAL,
        drift_cycles REAL, unique_count INTEGER,
        sigma_period_max_d REAL, ratio_to_quoted REAL,
        oc_mean_s REAL, oc_rms_s REAL, fitted_period_d REAL,
        fitted_period_sigma_d REAL, fitted_epoch_bjd REAL,
        n_epochs INTEGER, phase_spread REAL, one_feature INTEGER,
        verdict TEXT, note TEXT,
        n_night_epochs INTEGER, n_nights INTEGER, oc_night_rms_s REAL,
        oc_night_wrms_s REAL, oc_night_chi2nu REAL,
        sigma_night_median_s REAL, sigma_night_lo_s REAL,
        sigma_night_hi_s REAL, n_night_at_threshold INTEGER,
        fitted_period_night_d REAL, fitted_period_night_sigma_d REAL,
        fitted_epoch_night_bjd REAL, fit_night_chi2nu REAL,
        n_cycles_span REAL, span_d REAL,
        n_night_single_cycle INTEGER,
        oc_night_chi2nu_edge REAL, sigma_night_edge_median_s REAL,
        quad_coeff_s_per_cycle2 REAL, quad_sigma_s_per_cycle2 REAL,
        pdot REAL, pdot_sigma REAL, pdot_limit3 REAL);

    CREATE TABLE IF NOT EXISTS p3_state_night (
        series_key TEXT, night TEXT, target_key TEXT, era_id INTEGER,
        filter TEXT, n_points INTEGER, phase_coverage REAL,
        median_mag REAL, p10_mag REAL, p90_mag REAL, amplitude_mag REAL,
        n_limits INTEGER, censored INTEGER, gated INTEGER,
        state TEXT, note TEXT,
        PRIMARY KEY (series_key, night));

    CREATE TABLE IF NOT EXISTS p3_state_series (
        series_key TEXT PRIMARY KEY, target_key TEXT, era_id INTEGER,
        filter TEXT, n_nights INTEGER, n_gated INTEGER, n_used INTEGER,
        threshold_mag REAL, threshold_sigma REAL, separability REAL,
        bimodal INTEGER, n_high INTEGER, n_low INTEGER,
        n_intermediate INTEGER, n_censored INTEGER,
        duty_naive REAL, duty_with_limits REAL, duty_bias REAL,
        n_informative_limits INTEGER, n_uninformative INTEGER,
        span_mag REAL, verdict TEXT, note TEXT);

    CREATE TABLE IF NOT EXISTS p3_detrend (
        series_key TEXT, window_periods REAL, window_d REAL,
        amplitude_in REAL, frac_detrend REAL, frac_joint REAL,
        n_detrend INTEGER, n_joint INTEGER, backend TEXT,
        PRIMARY KEY (series_key, window_periods));

    CREATE TABLE IF NOT EXISTS p3_gp_check (
        series_key TEXT PRIMARY KEY, n_points INTEGER,
        amp_celerite REAL, amp_dense REAL, rel_diff REAL,
        loglike_celerite REAL, loglike_dense REAL, ll_abs_diff REAL,
        celerite_eps REAL, verdict TEXT);

    CREATE INDEX IF NOT EXISTS ix_p3_pgram
        ON p3_pgram (series_key, panel, kind);
    CREATE INDEX IF NOT EXISTS ix_p3_edge_target
        ON p3_edge (target_key, night, filter);
    """)
    # ``CREATE TABLE IF NOT EXISTS`` cannot widen a table that already
    # exists, and this archive predates the per-night O-C columns.  Add
    # them one at a time so an archive built by an older revision keeps its
    # rows instead of being rebuilt from scratch.
    have = {r[1] for r in con.execute("PRAGMA table_info(p3_cycle_count)")}
    for col, decl in (("n_night_epochs", "INTEGER"), ("n_nights", "INTEGER"),
                      ("oc_night_rms_s", "REAL"),
                      ("oc_night_wrms_s", "REAL"),
                      ("oc_night_chi2nu", "REAL"),
                      ("sigma_night_median_s", "REAL"),
                      ("sigma_night_lo_s", "REAL"),
                      ("sigma_night_hi_s", "REAL"),
                      ("n_night_at_threshold", "INTEGER"),
                      ("fitted_period_night_d", "REAL"),
                      ("fitted_period_night_sigma_d", "REAL"),
                      ("fitted_epoch_night_bjd", "REAL"),
                      ("fit_night_chi2nu", "REAL"),
                      # Added 2026-08-20 after the second referee report: the
                      # span of the EPOCHS (which is not the count from the
                      # catalogue epoch), the single-cycle epoch count, the
                      # error-transfer robustness check, and the period
                      # derivative the null actually bounds.
                      ("n_cycles_span", "REAL"), ("span_d", "REAL"),
                      ("n_night_single_cycle", "INTEGER"),
                      ("oc_night_chi2nu_edge", "REAL"),
                      ("sigma_night_edge_median_s", "REAL"),
                      ("quad_coeff_s_per_cycle2", "REAL"),
                      ("quad_sigma_s_per_cycle2", "REAL"),
                      ("pdot", "REAL"), ("pdot_sigma", "REAL"),
                      ("pdot_limit3", "REAL")):
        if col not in have:
            con.execute(f"ALTER TABLE p3_cycle_count ADD COLUMN {col} {decl}")
    have_n = {r[1] for r in con.execute("PRAGMA table_info(p3_oc_night)")}
    for col, decl in (("oc_sigma_edge_s", "REAL"),
                      ("budget_band", "TEXT")):
        if col not in have_n:
            con.execute(f"ALTER TABLE p3_oc_night ADD COLUMN {col} {decl}")
    con.commit()


def set_meta(con: sqlite3.Connection, items: dict) -> None:
    con.executemany("INSERT OR REPLACE INTO p3_meta (key, value) "
                    "VALUES (?, ?)",
                    [(k, str(v)) for k, v in items.items()])
    con.commit()


def stamp(con: sqlite3.Connection, stage: str) -> None:
    set_meta(con, {f"stage_{stage}": utcnow(),
                   "phase3_code_version": PHASE3_CODE_VERSION,
                   "git_commit": git_commit()})


def record_stage(key: str) -> None:
    """Tell the provenance graph this stage ran.  Never fatal: a failure to
    record is a bookkeeping problem, and killing a completed six-hour stage
    over one would be worse than the problem."""
    try:
        subprocess.run([sys.executable,
                        str(PIPELINE_ROOT / "scripts" /
                            "check_pipeline_status.py"), "record", key],
                       cwd=REPO_ROOT, capture_output=True, text=True,
                       timeout=180)
    except Exception as exc:                            # noqa: BLE001
        print(f"  ! could not record stage {key}: {exc}")


# ===========================================================================
# Reading the light curves
# ===========================================================================
def load_series(con: sqlite3.Connection, series_key: str) -> dict:
    """Every target point of one series, cloud-vetoed frames removed.

    ``cal_mag`` (the catalogue-tied natural-system magnitude) is the
    magnitude used throughout Phase 3, and rows without one are dropped
    rather than falling back to the instrumental ``mag``: mixing the two
    inside one series would put a zero-point step in the middle of a light
    curve and every period, edge and state statistic downstream would
    inherit it.

    The Phase-2 cloud veto is honoured here by JOINing ``p2_cloud_frame``.
    That table stores a FLAG rather than editing the light curve, which is
    what makes this join a decision this stage takes and can be argued
    with, instead of a fact it inherits.
    """
    rows = con.execute("""
        SELECT l.bjd_tdb, l.cal_mag, l.inst_mag_err, l.frame_id, f.night,
               COALESCE(c.vetoed, 0)
        FROM cv_lightcurve l
        JOIN cv_frames f
          ON f.frame_id = l.frame_id AND f.series_key = l.series_key
        LEFT JOIN p2_cloud_frame c
          ON c.frame_id = l.frame_id AND c.series_key = l.series_key
        WHERE l.series_key = ? AND l.role = 'target'
          AND l.cal_mag IS NOT NULL
        ORDER BY l.bjd_tdb
    """, (series_key,)).fetchall()
    if not rows:
        return {"t": np.array([]), "m": np.array([]), "e": np.array([]),
                "frame_id": np.array([]), "night": [], "n_vetoed": 0}
    keep = [r for r in rows if not r[5]]
    n_vetoed = len(rows) - len(keep)
    if not keep:
        return {"t": np.array([]), "m": np.array([]), "e": np.array([]),
                "frame_id": np.array([]), "night": [], "n_vetoed": n_vetoed}
    t = np.array([r[0] for r in keep], dtype=float)
    m = np.array([r[1] for r in keep], dtype=float)
    e = np.array([r[2] if r[2] is not None else np.nan for r in keep],
                 dtype=float)
    # A missing or absurd error bar becomes the series median rather than a
    # NaN that would silently drop the point out of every weighted fit.
    med = float(np.nanmedian(e)) if np.isfinite(e).any() else 0.02
    e = np.where(np.isfinite(e) & (e > 0), e, med)
    return {"t": t, "m": m, "e": e,
            "frame_id": np.array([r[3] for r in keep], dtype=np.int64),
            "night": [r[4] for r in keep], "n_vetoed": n_vetoed}


def series_rows(con: sqlite3.Connection) -> list[tuple]:
    """(series_key, target_key, era_id, filter) for every solved series."""
    return con.execute(
        "SELECT series_key, target_key, era_id, filter FROM cv_series "
        "WHERE status = 'solved' ORDER BY target_key, era_id, filter"
    ).fetchall()


def inflation_for(con: sqlite3.Connection, series_key: str) -> float:
    """The measured chi2 inflation, or 1.0 when the series has none.

    Every error bar this script uses is multiplied by it.  The
    characterization measured 0.92-3.02 across the sample, and an inflation
    below 1 is NOT applied — an error model that turns out to be slightly
    pessimistic does not license shrinking the bars, because the same
    measurement that says 0.92 has an uncertainty of its own.
    """
    row = con.execute("SELECT inflation FROM cv_error_model WHERE "
                      "series_key = ?", (series_key,)).fetchone()
    if not row or row[0] is None or not np.isfinite(row[0]):
        return 1.0
    return float(max(1.0, row[0]))


# ===========================================================================
# STAGE: ephem — the published ephemerides
# ===========================================================================
def fetch_vsx(ident: str) -> dict:
    """One VSX ``api.object`` record, as a dict.  Raises on failure."""
    url = ("https://www.aavso.org/vsx/index.php?view=api.object&ident="
           + urllib.parse.quote(ident) + "&format=json")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=NET_TIMEOUT_S) as r:
        return json.loads(r.read().decode("utf-8"))


def cmd_ephem(args) -> None:
    """Fetch and store the published ephemerides.

    Cached: the raw JSON payload is written to
    ``products/external/vsx/<key>.json`` and re-read on a later run unless
    ``--force``.  The cache is the point — every "agrees with the published
    value" on the report page traces to a stored payload with a fetch
    timestamp, not to a number in someone's memory.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    con = connect(args.db)
    ensure_tables(con)
    print(f"database: {args.db}")
    n_ok = 0
    for key, ident in VSX_IDENT.items():
        cache = CACHE_DIR / f"{key}.json"
        payload = None
        if cache.exists() and not args.force:
            payload = json.loads(cache.read_text(encoding="utf-8"))
            origin = "cache"
        else:
            try:
                payload = fetch_vsx(ident)
                cache.write_text(json.dumps(payload, indent=1),
                                 encoding="utf-8")
                origin = "network"
            except Exception as exc:                    # noqa: BLE001
                print(f"  {key:6s} FETCH FAILED: {exc}")
                continue
        obj = payload.get("VSXObject", {}) if payload else {}
        p_str = str(obj.get("Period") or "").strip()
        e_str = str(obj.get("Epoch") or "").strip()
        try:
            period = float(p_str)
        except ValueError:
            period = float("nan")
        try:
            epoch = float(e_str)
        except ValueError:
            epoch = float("nan")
        sigma = p3.quoted_precision_sigma(p_str)
        # VSX publishes NO period uncertainty.  Saying so, in the column
        # that carries the uncertainty, is the whole point of this field.
        basis = ("quoted precision of the VSX period string "
                 f"({p_str!r}); VSX publishes no uncertainty, so this is a "
                 "FLOOR and the true published sigma is at least this large")
        note = ""
        if not np.isfinite(epoch):
            note = ("VSX gives no epoch for this star, so no absolute phase, "
                    "no cycle count and no O-C are possible against it")
        con.execute("""
            INSERT OR REPLACE INTO p3_ephemeris
            (target_key, name, var_type, period_d, period_str,
             period_sigma_d, sigma_basis, epoch_bjd, epoch_str, source,
             fetched_utc, payload_path, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (key, str(obj.get("Name") or TARGET_LABEL.get(key, key)),
                     str(obj.get("VariabilityType") or ""),
                     period if np.isfinite(period) else None, p_str,
                     sigma if np.isfinite(sigma) else None, basis,
                     epoch if np.isfinite(epoch) else None, e_str,
                     "AAVSO VSX api.object", utcnow(),
                     str(cache.relative_to(REPO_ROOT)), note))
        n_ok += 1
        print(f"  {key:6s} {obj.get('Name'):10s} {obj.get('VariabilityType'):8s} "
              f"P={p_str:14s} +/-{sigma:.1e}  E={e_str or '(none)'}  [{origin}]")
    con.commit()
    stamp(con, "ephem")
    set_meta(con, {"n_ephemeris": n_ok})
    con.close()
    print(f"  stored {n_ok} ephemerides")


def load_ephemerides(con: sqlite3.Connection) -> dict[str, p3.Ephemeris]:
    """The stored ephemerides, keyed by target."""
    out = {}
    for r in con.execute(
            "SELECT target_key, name, period_d, epoch_bjd, source, "
            "period_str, period_sigma_d, var_type, note FROM p3_ephemeris"):
        if r[2] is None:
            continue
        out[r[0]] = p3.Ephemeris(
            target_key=r[0], name=r[1], period_d=float(r[2]),
            epoch_bjd=float(r[3]) if r[3] is not None else float("nan"),
            source=r[4], period_str=r[5] or "",
            period_sigma_d=float(r[6]) if r[6] is not None else None,
            var_type=r[7] or "", note=r[8] or "")
    return out


# ===========================================================================
# STAGE: periods
# ===========================================================================
def decimate(x, y, n: int = TRACE_POINTS) -> tuple[np.ndarray, np.ndarray]:
    """Keep ``n`` points of a trace WITHOUT losing the peaks.

    Plain slicing would step over a one-bin peak and draw a periodogram
    that is missing its own answer.  This bins the trace and keeps the
    EXTREME point of each bin (the maximum, since every trace stored here
    is one whose peaks matter), so the decimated curve still touches every
    feature the full-resolution one has.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size <= n:
        return x, y
    edges = np.linspace(0, x.size, n + 1).astype(int)
    xs, ys = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        seg = y[a:b]
        if not np.isfinite(seg).any():
            continue
        i = a + int(np.nanargmax(seg))
        xs.append(x[i])
        ys.append(y[i])
    return np.asarray(xs), np.asarray(ys)


def analyse_period(payload: dict) -> dict:
    """One series' period analysis.  Pure enough to run in a worker."""
    t = np.asarray(payload["t"], dtype=float)
    m = np.asarray(payload["m"], dtype=float)
    e = np.asarray(payload["e"], dtype=float)
    sk = payload["series_key"]
    pub_p = payload["published_d"]
    pub_s = payload["published_sigma_d"]
    out: dict = {"series_key": sk, "n_points": int(t.size)}
    if t.size < MIN_POINTS_FOR_SEARCH:
        out["status"] = "too_few_points"
        return out
    blocks = p3.block_index(t)
    n_blocks = int(blocks.max() + 1)
    baseline = float(t.max() - t.min())
    within = np.concatenate([np.diff(t[blocks == b]) for b in range(n_blocks)
                             if (blocks == b).sum() > 1] or [np.array([])])
    cadence_s = float(np.median(within) * 86400.0) if within.size else float("nan")
    out.update(n_blocks=n_blocks, baseline_d=baseline,
               median_cadence_s=cadence_s)
    # --- the survey band, so the page can show what the tallest peak IS ---
    g_survey = p3.frequency_grid(baseline, p3.SURVEY_F_MIN_CD,
                                 p3.SURVEY_F_MAX_CD)
    pw_survey = p3.gls_block_power(t, m, e, g_survey, blocks)
    i = int(np.nanargmax(pw_survey))
    out["f_survey_cd"] = float(g_survey[i])
    out["p_survey_pow"] = float(pw_survey[i])
    out["survey_is_orbital"] = int(p3.ORBITAL_F_MIN_CD <= g_survey[i]
                                   <= p3.ORBITAL_F_MAX_CD)
    # WHAT the tallest peak in the whole search actually is, when it is not
    # the orbit.  Classified rather than lumped together, because the two
    # cases mean opposite things: a peak at 2x the orbital frequency is the
    # orbit, seen through a double-humped light curve that puts more power
    # in its first harmonic than in its fundamental (all three AN UMa
    # series do this, at 2.00x), while a peak below 3 c/d is the accretion
    # state changing on week-to-month timescales and is a different
    # astrophysical signal entirely.  An earlier draft of the report
    # described all of them as "near 2 c/d", which was true of six series
    # and false of five.
    out["survey_class"] = "orbit"
    if not out["survey_is_orbital"]:
        pubf = (1.0 / payload["published_d"]
                if payload["published_d"] and np.isfinite(payload["published_d"])
                else float("nan"))
        ratio = g_survey[i] / pubf if np.isfinite(pubf) and pubf > 0 else np.nan
        if np.isfinite(ratio) and any(abs(ratio - n) < 0.1 for n in (2, 3)):
            out["survey_class"] = f"harmonic {ratio:.2f}x f_orb"
        elif g_survey[i] < 3.0:
            out["survey_class"] = "slow (accretion state)"
        else:
            out["survey_class"] = "other"
    # --- the orbital band ---
    g_orb = p3.frequency_grid(baseline, p3.ORBITAL_F_MIN_CD,
                              p3.ORBITAL_F_MAX_CD)
    pw_orb = p3.gls_block_power(t, m, e, g_orb, blocks)
    theta = p3.pdm_theta(t, m, g_orb, blocks=blocks)
    j = int(np.nanargmax(pw_orb))
    k = int(np.nanargmin(theta))
    out["f_ls_band_cd"] = float(g_orb[j])
    out["p_ls_band_pow"] = float(pw_orb[j])
    out["f_pdm_band_cd"] = float(g_orb[k])
    out["pdm_band_theta"] = float(theta[k])
    # --- the window: the reason this whole stage exists ---
    fr = p3.alias_window_fractions(t)
    out["alias_frac_p1"] = fr.get(1, float("nan"))
    out["alias_frac_m1"] = fr.get(-1, float("nan"))
    out["alias_frac_max"] = float(max(fr.values())) if fr else float("nan")
    f_pub = 1.0 / pub_p if pub_p and np.isfinite(pub_p) and pub_p > 0 else float("nan")
    # Local refinement AT the published frequency.  The half-width must be
    # a few RESOLUTION ELEMENTS (1/T), not a fixed fraction of the alias
    # spacing: on a 400 d baseline the peak is 0.0025 c/d wide, and a
    # +/-0.45 c/d window would happily walk onto a completely different
    # feature 180 resolution elements away and report it as the orbital
    # peak.  Three elements, floored so a very long baseline still gets a
    # window wider than the grid step, and capped below half the alias
    # spacing so it can never cross into a neighbouring family member.
    halfwidth = float(min(0.45, max(3.0 / baseline, 0.002)))
    out["refine_halfwidth_cd"] = halfwidth
    if np.isfinite(f_pub):
        f_ls, p_ls = p3.refine_peak(g_orb, pw_orb, f_pub, halfwidth)
        f_pdm, th = p3.refine_trough(g_orb, theta, f_pub, halfwidth)
        out["f_ls_cd"], out["p_ls_pow"] = f_ls, p_ls
        out["f_pdm_cd"], out["pdm_theta"] = f_pdm, th
        out["ls_minus_pub_cd"] = float(f_ls - f_pub)
        out["pdm_minus_pub_cd"] = float(f_pdm - f_pub)
        out["peak_halfwidth_cd"] = p3.peak_halfwidth(g_orb, pw_orb, f_ls)
        amp = p3.gls_amplitude(t, m, e, f_ls, blocks)
        out["amplitude_mag"] = amp
        # Residual scatter after the best sinusoid + night constants.
        des = p3.sinusoid_design(t - t.mean(), f_ls)
        w = 1.0 / e ** 2
        y0 = p3._project_out_blocks(m[None, :], blocks, w)[0]
        d0 = p3._project_out_blocks(des.T, blocks, w).T
        try:
            coef = np.linalg.solve((d0.T * w) @ d0, (d0.T * w) @ y0)
            resid = y0 - d0 @ coef
            out["resid_rms_mag"] = float(np.std(resid, ddof=1))
        except np.linalg.LinAlgError:
            out["resid_rms_mag"] = float(np.std(m, ddof=1))
        period = 1.0 / f_ls if np.isfinite(f_ls) and f_ls > 0 else float("nan")
        out["period_d"] = period
        out["sigma_formal_d"] = p3.mhb_period_sigma(
            t.size, baseline, out["amplitude_mag"], out["resid_rms_mag"],
            period)
        sig_b, med_b, n_b = p3.night_bootstrap_period(
            t, m, e, blocks, f_pub, min(halfwidth, 0.15),
            n_boot=payload.get("n_boot", p3.N_BOOT))
        out["sigma_boot_d"], out["n_boot"] = sig_b, n_b
        cands = [v for v in (out["sigma_formal_d"], sig_b) if np.isfinite(v)]
        # THE RESOLUTION FLOOR.  The night bootstrap needs at least three
        # nights; with one or two there is nothing to resample and the only
        # surviving estimate is the analytic one, which assumes white noise
        # and no red noise — on a two-night VV Pup set that produced a
        # 5e-5 d error bar and an 18-sigma "disagreement" with a period the
        # data cannot actually distinguish.  When the empirical check is
        # unavailable, fall back to the frequency resolution 1/(2T), which
        # is the width of the thing being centroided.
        floor = float("nan")
        if n_blocks < 3 and np.isfinite(period) and baseline > 0:
            floor = period ** 2 / (2.0 * baseline)
            cands.append(floor)
        out["sigma_resolution_d"] = floor
        out["sigma_period_d"] = float(max(cands)) if cands else float("nan")
        out["sigma_basis"] = (
            "larger of the Montgomery & O'Donoghue analytic sigma and a "
            f"{n_b}-replicate NIGHT bootstrap; the bootstrap resamples "
            "whole nights because the noise is correlated within one and a "
            "point bootstrap would shrink with the exposure count instead "
            "of with the number of independent nights"
            + ("" if n_blocks >= 3 else
               f"; fewer than 3 nights, so the bootstrap is unavailable and "
               f"the 1/(2T) frequency-resolution floor {floor:.2e} d is "
               f"applied in its place"))
        # How much the period is actually CONSTRAINED, which is a different
        # question from whether it agrees.  A single night agrees with
        # everything.
        frac = (out["sigma_period_d"] / period
                if np.isfinite(out["sigma_period_d"]) and period else float("nan"))
        out["frac_precision"] = frac
        out["constraint"] = ("UNKNOWN" if not np.isfinite(frac) else
                             "TIGHT" if frac < 1e-3 else
                             "WEAK" if frac < 1e-2 else "UNINFORMATIVE")
        dev, ok, sentence = p3.agreement(period, out["sigma_period_d"],
                                         pub_p, pub_s)
        out["deviation_sigma"], out["agrees"], out["agree_note"] = dev, ok, sentence
        code, note = p3.classify_family_choice(
            n_blocks, out["alias_frac_max"], out["peak_halfwidth_cd"])
        out["family_code"], out["family_note"] = code, note
        # Is the tallest peak in the band a HARMONIC of the orbit?  AN UMa's
        # is: a double-humped light curve puts more power at 2f than at f,
        # and a search that reported the tallest peak would publish half the
        # period.
        ratio = out["f_ls_band_cd"] / f_pub if f_pub else float("nan")
        harm = ""
        if abs(out["f_ls_band_cd"] - f_pub) > 0.1:
            for n in (2, 3):
                if abs(ratio - n) < 0.02:
                    harm = (f"the tallest peak in the orbital band is at "
                            f"{ratio:.2f} x the orbital frequency — the "
                            f"{n}x harmonic of a non-sinusoidal light curve, "
                            "not a different period")
            if not harm:
                off = out["f_ls_band_cd"] - f_pub
                if abs(abs(off) - 1.0) < 0.1:
                    harm = (f"the tallest peak in the orbital band is "
                            f"{off:+.3f} c/d from the published frequency — "
                            "the +/-1 c/d ALIAS, which is what the window "
                            "predicts and what makes the prior necessary")
        out["harmonic_note"] = harm
        # Detected at all?  A PDM trough that never drops and an LS peak
        # with no power mean the orbit is not in this series, which is a
        # result (YZ Cnc is a dwarf nova, not a polar).
        out["detected"] = int(np.isfinite(th) and th < 0.9
                              and np.isfinite(p_ls) and p_ls > 0.05)
        if not out["detected"]:
            # An agreement test on an undetected signal compares two noise
            # peaks.  Blank it rather than let a wide error bar manufacture
            # a pass: "no modulation found" is the result here.
            out["agrees"] = None
            out["agree_note"] = (
                "No orbital modulation detected in this series (PDM theta "
                f"{th:.3f}, Lomb-Scargle power {p_ls:.3f} at the published "
                "frequency), so there is no recovered period to compare and "
                "the agreement test is not applicable.")
    else:
        out["detected"] = 0
        out["family_code"], out["family_note"] = "NO EPHEMERIS", (
            "No published period for this target, so there is no family "
            "member to select and no agreement to test.")
    out["status"] = "ok"
    # Traces for the figures: the survey panel and a zoom on the orbital
    # peak, with the SPECTRAL WINDOW on the same frequency axis as each.
    traces = []
    fx, fy = decimate(g_survey, pw_survey)
    traces += [("survey", "ls", a, b) for a, b in zip(fx, fy)]
    wx = np.linspace(-3.0, 3.0, 1201)
    ww = p3.spectral_window_power(t, wx)
    traces += [("window", "window", a, b) for a, b in zip(wx, ww)]
    if np.isfinite(f_pub):
        zoom = p3.frequency_grid(baseline, f_pub - 2.5, f_pub + 2.5)
        pz = p3.gls_block_power(t, m, e, zoom, blocks)
        tz = p3.pdm_theta(t, m, zoom, blocks=blocks)
        zx, zy = decimate(zoom, pz)
        traces += [("zoom", "ls", a, b) for a, b in zip(zx, zy)]
        zx2, zy2 = decimate(zoom, -tz)
        traces += [("zoom", "pdm", a, -b) for a, b in zip(zx2, zy2)]
    out["traces"] = traces
    return out


def cmd_periods(args) -> None:
    """Verify the orbital period per (target, era, filter)."""
    con = connect(args.db)
    ensure_tables(con)
    eph = load_ephemerides(con)
    if not eph:
        print("  ! no ephemerides stored — run `ephem` first")
        con.close()
        return
    rows = series_rows(con)
    jobs = []
    for sk, tk, era, filt in rows:
        data = load_series(con, sk)
        e = eph.get(tk)
        jobs.append({
            "series_key": sk, "target_key": tk, "era_id": era,
            "filter": filt, "t": data["t"], "m": data["m"],
            "e": data["e"] * inflation_for(con, sk),
            "published_d": e.period_d if e else float("nan"),
            "published_sigma_d": (e.period_sigma_d if e and e.period_sigma_d
                                  else float("nan")),
            "n_boot": args.n_boot})
    print(f"database: {args.db}")
    print(f"  {len(jobs)} series, {args.workers} workers")
    results = []
    t0 = time.time()
    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=min(args.workers,
                                                 MAX_WORKERS)) as pool:
            futs = {pool.submit(analyse_period, j): j["series_key"]
                    for j in jobs}
            for fut in as_completed(futs):
                results.append((futs[fut], fut.result()))
                print(f"    {futs[fut]:16s} done ({len(results)}/{len(jobs)})",
                      flush=True)
    else:
        for j in jobs:
            results.append((j["series_key"], analyse_period(j)))
            print(f"    {j['series_key']:16s} done", flush=True)
    by_key = {j["series_key"]: j for j in jobs}
    con.execute("DELETE FROM p3_pgram")
    for sk, res in sorted(results):
        j = by_key[sk]
        eobj = eph.get(j["target_key"])
        con.execute("""
            INSERT OR REPLACE INTO p3_period
            (series_key, target_key, era_id, filter, n_points, n_blocks,
             baseline_d, median_cadence_s, phase_coverage,
             f_survey_cd, p_survey_pow, survey_is_orbital, survey_class,
             f_ls_cd, p_ls_pow, f_pdm_cd, pdm_theta,
             ls_minus_pub_cd, pdm_minus_pub_cd,
             alias_frac_p1, alias_frac_m1, alias_frac_max,
             peak_halfwidth_cd, amplitude_mag, resid_rms_mag,
             period_d, sigma_period_d, sigma_formal_d, sigma_boot_d,
             sigma_resolution_d, n_boot, sigma_basis, refine_halfwidth_cd,
             frac_precision, constraint_class,
             published_d, published_sigma_d,
             deviation_sigma, agrees, agree_note, family_code, family_note,
             harmonic_note, detected, status, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
                    ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
            sk, j["target_key"], j["era_id"], j["filter"],
            res.get("n_points"), res.get("n_blocks"), res.get("baseline_d"),
            res.get("median_cadence_s"),
            (p3.phase_coverage(j["t"], eobj.period_d, eobj.epoch_bjd)
             if eobj and np.isfinite(eobj.epoch_bjd) and j["t"].size else None),
            res.get("f_survey_cd"), res.get("p_survey_pow"),
            res.get("survey_is_orbital"), res.get("survey_class"),
            res.get("f_ls_cd"), res.get("p_ls_pow"), res.get("f_pdm_cd"),
            res.get("pdm_theta"), res.get("ls_minus_pub_cd"),
            res.get("pdm_minus_pub_cd"), res.get("alias_frac_p1"),
            res.get("alias_frac_m1"), res.get("alias_frac_max"),
            res.get("peak_halfwidth_cd"), res.get("amplitude_mag"),
            res.get("resid_rms_mag"), res.get("period_d"),
            res.get("sigma_period_d"), res.get("sigma_formal_d"),
            res.get("sigma_boot_d"), res.get("sigma_resolution_d"),
            res.get("n_boot"), res.get("sigma_basis"),
            res.get("refine_halfwidth_cd"), res.get("frac_precision"),
            res.get("constraint"), j["published_d"],
            j["published_sigma_d"], res.get("deviation_sigma"),
            (int(res["agrees"]) if res.get("agrees") is not None else None),
            res.get("agree_note"), res.get("family_code"),
            res.get("family_note"), res.get("harmonic_note"),
            res.get("detected"), res.get("status"),
            f"{res.get('n_points', 0)} points after the Phase-2 cloud veto"))
        con.executemany(
            "INSERT INTO p3_pgram (series_key, panel, kind, freq_cd, value) "
            "VALUES (?,?,?,?,?)",
            [(sk, a, b, float(c), float(d))
             for a, b, c, d in res.get("traces", [])])
    con.commit()
    stamp(con, "periods")
    n_ok = con.execute("SELECT count(*) FROM p3_period WHERE "
                       "status='ok'").fetchone()[0]
    n_agree = con.execute("SELECT count(*) FROM p3_period WHERE "
                          "agrees=1").fetchone()[0]
    set_meta(con, {"n_period_series": n_ok, "n_period_agree": n_agree,
                   "orbital_band_cd": f"{p3.ORBITAL_F_MIN_CD}-"
                                      f"{p3.ORBITAL_F_MAX_CD}",
                   "alias_decidable_max": p3.ALIAS_DECIDABLE_MAX})
    con.close()
    print(f"  {n_ok} series analysed in {time.time() - t0:.0f}s, "
          f"{n_agree} agree with the published period")


# ===========================================================================
# STAGE: sigmat
# ===========================================================================
def measure_bright_phase_shape(t, m, period_d, epoch_d, err=None) -> dict:
    """Depth, bright-phase width and edge width, measured from a fold.

    The injection test has to inject something REAL, and these three
    numbers are what "real" means.

    Depth and bright-phase width come from a 25-bin phase-folded median
    profile: the depth is the 10th-to-90th percentile span of the profile,
    the width is the fraction of bins above the half-depth level.

    THE EDGE WIDTH IS NOT MEASURED FROM THAT PROFILE, and the reason is
    worth stating.  A 25-bin profile of a 0.0791 d orbit has a bin 273 s
    wide, so no profile-based estimator can report an edge sharper than
    about 273 s, and the first version of this function was worse than that:
    it was algebraically pinned at exactly two bins, 547 s, for every series
    it was ever given.  Injecting a 547 s ramp when the true one may be
    150 s makes the recovery far too easy and flatters the whole sigma_t
    contour.

    So the edge width is fitted to the FOLDED POINTS at full time
    resolution instead.  Folding five cycles that are each sampled every
    219 s gives an effective phase sampling several times finer than one
    cycle's cadence, because consecutive cycles sample different phases —
    which is exactly the resolution the profile threw away by binning.
    ``edge_width_floor_d`` is returned beside it so the caller can see when
    the fit has hit the sampling floor and the width is an upper bound.
    """
    ph = p3.phase_of(t, period_d, epoch_d)
    nb = 25
    idx = np.clip((ph * nb).astype(int), 0, nb - 1)
    prof = np.array([np.median(m[idx == i]) if (idx == i).sum() else np.nan
                     for i in range(nb)])
    if np.isfinite(prof).sum() < 8:
        return {}
    bright = float(np.nanpercentile(prof, 10))
    faint = float(np.nanpercentile(prof, 90))
    depth = faint - bright
    half = bright + depth / 2.0
    bright_bins = np.flatnonzero(prof < half)
    if bright_bins.size == 0:
        return {}
    width_phase = bright_bins.size / nb
    # Circular mean of the bright bins, so a bright phase straddling phase
    # zero (ST LMi's does) gets its centre in the right place instead of at
    # 0.5.  The falling edge is half a bright-phase width later.
    centre = float(np.angle(np.exp(2j * np.pi * bright_bins / nb).sum())
                   / (2 * np.pi)) % 1.0
    edge_phase = (centre + width_phase / 2.0) % 1.0
    # --- the edge width, fitted to the folded points ---
    e = (np.asarray(err, dtype=float) if err is not None
         else np.full(t.shape, max(0.01, depth / 50.0)))
    # Signed distance from the edge in DAYS, folded to +/- half a period.
    dt = np.mod(np.asarray(ph) - edge_phase + 0.5, 1.0) - 0.5
    dt = dt * period_d
    sel = np.abs(dt) <= 0.16 * period_d
    med_dt_s = float(np.median(np.diff(np.sort(dt[sel])))) * 86400.0 \
        if sel.sum() > 3 else float("nan")
    floor_d = max(med_dt_s, 20.0) / 86400.0 if np.isfinite(med_dt_s) else 60.0 / 86400.0
    edge_width_d = 2.0 / nb * period_d          # fallback: the old estimate
    edge_basis = "profile-bin fallback (the folded fit did not converge)"
    if sel.sum() >= p3.EDGE_MIN_POINTS:
        widths = np.geomspace(floor_d, 0.30 * period_d, 24)
        fit = p3.fit_edge(dt[sel], m[sel], e[sel],
                          p3.edge_time_grid(0.0, 0.05 * period_d, 201),
                          widths, med_dt_s, min_snr=0.0,
                          max_bracket_cadence=1e9)
        if np.isfinite(fit.width_d) and fit.width_d > 0:
            edge_width_d = float(fit.width_d)
            at_floor = edge_width_d <= widths[0] * 1.001
            edge_basis = (
                f"fitted to {int(sel.sum())} folded points at "
                f"{med_dt_s:.0f} s effective phase sampling"
                + (" — AT THE SAMPLING FLOOR, so this is an upper bound on "
                   "the true ramp" if at_floor else ""))
    return {"depth_mag": depth, "bright_width_phase": width_phase,
            "edge_width_d": edge_width_d, "edge_width_floor_d": floor_d,
            "edge_basis": edge_basis, "edge_phase": edge_phase,
            "profile": prof, "bright": bright, "faint": faint,
            "effective_sampling_s": med_dt_s}


def cmd_sigmat(args) -> None:
    """The sigma_t injection test on ST LMi's densest night."""
    con = connect(args.db)
    ensure_tables(con)
    eph = load_ephemerides(con)
    e_obj = eph.get(SIGMAT_TARGET)
    if not e_obj:
        print("  ! no ephemeris for ST LMi — run `ephem` first")
        con.close()
        return
    print(f"database: {args.db}")
    print(f"  night {SIGMAT_NIGHT} (local; UTC night {SIGMAT_NIGHT[:8]}28), "
          f"target {TARGET_LABEL[SIGMAT_TARGET]}")
    keys = [r[0] for r in con.execute(
        "SELECT DISTINCT series_key FROM cv_frames WHERE target_key=? "
        "AND night=? ORDER BY series_key", (SIGMAT_TARGET, SIGMAT_NIGHT))]
    n_cells = 0
    for sk in keys:
        rows = con.execute("""
            SELECT l.bjd_tdb, l.cal_mag, l.inst_mag_err
            FROM cv_lightcurve l
            JOIN cv_frames f ON f.frame_id = l.frame_id
                            AND f.series_key = l.series_key
            LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                            AND c.series_key = l.series_key
            WHERE l.series_key=? AND l.role='target' AND f.night=?
              AND l.cal_mag IS NOT NULL AND COALESCE(c.vetoed,0)=0
            ORDER BY 1""", (sk, SIGMAT_NIGHT)).fetchall()
        if len(rows) < 30:
            print(f"  {sk:16s} only {len(rows)} points — skipped")
            continue
        a = np.array(rows, dtype=float)
        t, m, e = a[:, 0], a[:, 1], a[:, 2]
        infl = inflation_for(con, sk)
        e = np.where(np.isfinite(e) & (e > 0), e, np.nanmedian(e)) * infl
        cadence = float(np.median(np.diff(t)) * 86400.0)
        shape = measure_bright_phase_shape(t, m, e_obj.period_d,
                                           e_obj.epoch_bjd, err=e)
        if not shape:
            print(f"  {sk:16s} could not measure a bright-phase profile")
            continue
        n_cycles = float((t.max() - t.min()) / e_obj.period_d)
        con.execute("""
            INSERT OR REPLACE INTO p3_sigmat_input
            (series_key, night, n_points, median_cadence_s, median_err_mag,
             chi2_inflation, used_err_mag, depth_mag, edge_width_s,
             edge_width_floor_s, bright_width_phase, n_cycles, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
            sk, SIGMAT_NIGHT, len(t), cadence,
            float(np.median(a[:, 2])), infl, float(np.median(e)),
            shape["depth_mag"], shape["edge_width_d"] * 86400.0,
            shape["edge_width_floor_d"] * 86400.0,
            shape["bright_width_phase"], n_cycles,
            "depth and bright-phase width from this night's own 25-bin "
            "phase-folded median profile; edge width " + shape["edge_basis"]))
        print(f"  {sk:16s} n={len(t)} cadence={cadence:.0f}s "
              f"depth={shape['depth_mag']:.2f} "
              f"edge={shape['edge_width_d'] * 86400:.0f}s "
              f"bright_width={shape['bright_width_phase']:.2f} "
              f"cycles={n_cycles:.2f}")
        grid = []
        nominal = []
        for wf in INJECT_WIDTH_FACTORS:
            inj_w = shape["edge_width_d"] * wf
            for se in SHAPE_ERRORS:
                for de in DEPTH_ERRORS:
                    res = p3.sigma_t_injection(
                        t, e, e_obj.period_d, depth_mag=shape["depth_mag"],
                        edge_width_d=inj_w,
                        bright_width_phase=shape["bright_width_phase"],
                        shape_error=se, depth_error=de,
                        n_real=args.n_real, median_cadence_s=cadence)
                    res["shape_error"], res["depth_error"] = se, de
                    res["inject_factor"] = wf
                    grid.append(res)
                    if wf == 1.0:
                        nominal.append(res)
                    con.execute("""
                        INSERT OR REPLACE INTO p3_sigmat
                        (series_key, night, inject_factor, inject_width_s,
                         shape_error, depth_error, sigma_t_s, bias_s,
                         total_error_s, rms_s, p95_abs_s, n_ok, n_try,
                         recovered_fraction, passes)
                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                        sk, SIGMAT_NIGHT, wf, inj_w * 86400.0, se, de,
                        res["sigma_t_s"], res["bias_s"],
                        res["total_error_s"], res.get("sigma_t_rms_s"),
                        res["p95_abs_s"], res["n_ok"], res["n_try"],
                        res["recovered_fraction"],
                        int(np.isfinite(res["total_error_s"]) and
                            res["total_error_s"] <= p3.SIGMA_T_THRESHOLD_S)))
                    n_cells += 1
                    print(f"      inject {inj_w * 86400:5.0f}s  shape x{se:.0f}"
                          f" depth {de:+4.0%}  sigma_t={res['sigma_t_s']:6.1f}s"
                          f" bias={res['bias_s']:+6.1f}s"
                          f" total={res['total_error_s']:6.1f}s"
                          f" ({res['n_ok']}/{res['n_try']})", flush=True)
            con.commit()
        # The headline verdict is taken on the NOMINAL injected width (the
        # one the fold actually measured); the width sensitivity is reported
        # separately so a reader can see it did not decide the answer.
        verdict, sentence = p3.contour_verdict(nominal)
        wide, wide_sentence = p3.contour_verdict(grid)
        best_by_width = {}
        for r in grid:
            if r["shape_error"] == 1.0 and r["depth_error"] == 0.0:
                best_by_width[r["inject_factor"]] = r["total_error_s"]
        set_meta(con, {f"sigmat_verdict_{sk}": verdict,
                       f"sigmat_sentence_{sk}": sentence,
                       f"sigmat_verdict_allwidths_{sk}": wide,
                       f"sigmat_width_scan_{sk}": "; ".join(
                           f"x{k:g} -> {v:.1f} s"
                           for k, v in sorted(best_by_width.items()))})
        print(f"    -> {verdict}: {sentence}")
        print(f"    -> width sensitivity (shape known exactly): "
              + "; ".join(f"inject x{k:g} -> {v:.1f} s"
                          for k, v in sorted(best_by_width.items())))
        con.commit()
    stamp(con, "sigmat")
    set_meta(con, {"sigmat_threshold_s": p3.SIGMA_T_THRESHOLD_S,
                   "sigmat_n_real": args.n_real,
                   "sigmat_night": SIGMAT_NIGHT,
                   "sigmat_shape_errors": ",".join(str(s) for s in SHAPE_ERRORS),
                   "sigmat_depth_errors": ",".join(str(s) for s in DEPTH_ERRORS)})
    con.close()
    print(f"  {n_cells} grid cells measured")


# ===========================================================================
# STAGE: edges
# ===========================================================================
#: Significance bar for an inter-band edge-time offset, in sigma.  A pooled
#: difference is called a DETECTION only above this; below it the row is a
#: null and must be published as one, with its value and its error bar, not
#: as "the offsets are not uniformly zero" -- which is true of any set of
#: measured differences and carries no inference at all.
BAND_PAIR_SIGMA_BAR = 3.0


def cmd_edges(args) -> None:
    """Bright-phase edge epochs, per cycle, per band, for the polars."""
    con = connect(args.db)
    ensure_tables(con)
    eph = load_ephemerides(con)
    print(f"database: {args.db}")
    n_fit = n_acc = 0
    for sk, tk, era, filt in series_rows(con):
        if tk not in POLARS:
            continue
        e_obj = eph.get(tk)
        if not e_obj or not np.isfinite(e_obj.epoch_bjd):
            continue
        data = load_series(con, sk)
        if data["t"].size < 30:
            continue
        infl = inflation_for(con, sk)
        t, m, e = data["t"], data["m"], data["e"] * infl
        nights = np.array(data["night"])
        # The bright-phase profile tells us WHERE the falling edge is in
        # phase, so the per-cycle fits are seeded from this series' own
        # light curve and not from a number typed into this script.
        shape = measure_bright_phase_shape(t, m, e_obj.period_d,
                                           e_obj.epoch_bjd, err=e)
        if not shape:
            continue
        # The falling edge's phase, from this series' own folded profile.
        edge_phase = shape["edge_phase"]
        within = np.concatenate([np.diff(t[nights == n]) for n in np.unique(nights)
                                 if (nights == n).sum() > 1] or [np.array([])])
        cadence = float(np.median(within) * 86400.0) if within.size else 219.0
        # The Monte Carlo's sigma_t for the matching cell, if it exists.
        # The injection Monte Carlo's total error for the SHAPE-KNOWN,
        # DEPTH-KNOWN cell at the nominal injected width — the most
        # optimistic cell in the grid, and therefore a LOWER bound on the
        # real per-cycle uncertainty.  Used as a floor under the formal bar,
        # never as a replacement for it.
        mc = con.execute(
            "SELECT total_error_s FROM p3_sigmat WHERE series_key=? AND "
            "inject_factor=1.0 AND shape_error=1.0 AND depth_error=0.0",
            (sk,)).fetchone()
        mc_sigma = float(mc[0]) if mc and mc[0] is not None else None
        cyc_all = p3.cycle_number(t, e_obj.period_d, e_obj.epoch_bjd)
        for cyc in np.unique(cyc_all):
            t_guess = (e_obj.epoch_bjd + (cyc + edge_phase) * e_obj.period_d)
            sel = np.abs(t - t_guess) <= p3.EDGE_WINDOW_PHASE * e_obj.period_d
            if sel.sum() < p3.EDGE_MIN_POINTS:
                continue
            tt, mm, ee = t[sel], m[sel], e[sel]
            # Robust pre-clip: a single 4-magnitude outlier (there is one in
            # the r-band dense night) would otherwise set both levels.
            med = float(np.median(mm))
            mad = 1.4826 * float(np.median(np.abs(mm - med)))
            keep = np.abs(mm - med) < 6.0 * max(mad, 0.05)
            tt, mm, ee = tt[keep], mm[keep], ee[keep]
            if tt.size < p3.EDGE_MIN_POINTS:
                continue
            half_span = 3.0 * cadence / 86400.0
            tg = p3.edge_time_grid(t_guess, half_span, 361)
            wg = np.array([0.5, 1.0, 2.0, 4.0]) * shape["edge_width_d"]
            fit = p3.fit_edge(tt, mm, ee, tg, wg, cadence)
            n_fit += 1
            n_acc += int(fit.accepted)
            con.execute("""
                INSERT OR REPLACE INTO p3_edge
                (series_key, cycle, target_key, era_id, filter, night,
                 t_edge_bjd, sigma_t_s, sigma_t_mc_s, width_s,
                 level_bright, level_faint, depth_mag, depth_snr, chi2nu,
                 n_points, bracket_s, phase, accepted, reason)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                sk, int(cyc), tk, era, filt,
                str(nights[np.argmin(np.abs(t - t_guess))]),
                fit.t_edge_d if np.isfinite(fit.t_edge_d) else None,
                fit.sigma_t_s if np.isfinite(fit.sigma_t_s) else None,
                mc_sigma, fit.width_d * 86400.0 if np.isfinite(fit.width_d) else None,
                fit.level_bright, fit.level_faint, fit.depth_mag,
                fit.depth_snr if np.isfinite(fit.depth_snr) else None,
                fit.chi2nu, fit.n_points,
                fit.bracket_s if np.isfinite(fit.bracket_s) else None,
                float(p3.phase_of([fit.t_edge_d], e_obj.period_d,
                                  e_obj.epoch_bjd)[0])
                if np.isfinite(fit.t_edge_d) else None,
                int(fit.accepted), fit.reason))
        con.commit()
        n_here = con.execute("SELECT count(*), sum(accepted) FROM p3_edge "
                             "WHERE series_key=?", (sk,)).fetchone()
        print(f"  {sk:16s} edge phase {edge_phase:.3f}  "
              f"{int(n_here[1] or 0):3d}/{int(n_here[0] or 0):3d} accepted")
    # --- inter-band differences: the cyclotron measurement ---
    con.execute("DELETE FROM p3_band_pair")
    pairs = con.execute("""
        SELECT DISTINCT a.target_key, a.era_id, a.night, a.filter, b.filter
        FROM p3_edge a JOIN p3_edge b
          ON a.target_key=b.target_key AND a.era_id=b.era_id
         AND a.night=b.night AND a.cycle=b.cycle AND a.filter < b.filter
        WHERE a.accepted=1 AND b.accepted=1
        ORDER BY 1,2,3,4,5""").fetchall()
    for tk, era, night, fa, fb in pairs:
        rows = con.execute("""
            SELECT a.t_edge_bjd, a.sigma_t_s, b.t_edge_bjd, b.sigma_t_s,
                   a.sigma_t_mc_s, b.sigma_t_mc_s
            FROM p3_edge a JOIN p3_edge b
              ON a.target_key=b.target_key AND a.era_id=b.era_id
             AND a.night=b.night AND a.cycle=b.cycle
            WHERE a.target_key=? AND a.era_id=? AND a.night=?
              AND a.filter=? AND b.filter=? AND a.accepted=1 AND b.accepted=1
            """, (tk, era, night, fa, fb)).fetchall()
        if not rows:
            continue
        arr = np.array([[r[0], r[1], r[2], r[3]] for r in rows], dtype=float)
        # Error bars: the injection Monte Carlo where it exists, the
        # rescaled formal bar otherwise, and NEVER the smaller of the two.
        mc_a = rows[0][4]
        mc_b = rows[0][5]
        sa = np.maximum(arr[:, 1], mc_a if mc_a else 0.0)
        sb = np.maximum(arr[:, 3], mc_b if mc_b else 0.0)
        d, s, chi2nu = p3.band_difference(arr[:, 0], sa, arr[:, 2], sb)
        con.execute("""
            INSERT OR REPLACE INTO p3_band_pair
            (target_key, era_id, night, band_a, band_b, n_cycles, delta_s,
             sigma_s, chi2nu, significant, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)""", (
            tk, era, night, fa, fb, len(rows), d, s, chi2nu,
            int(np.isfinite(d) and np.isfinite(s) and s > 0
                and abs(d) > BAND_PAIR_SIGMA_BAR * s),
            ("per-cycle paired difference; the error bar is the larger of "
             "the rescaled formal bar and the injection Monte Carlo's "
             "sigma_t, inflated further when the offset is not constant "
             "from cycle to cycle")))
    # --- POOLED across every night: the actual cyclotron measurement ---
    # Most individual nights contribute a single paired cycle, and a single
    # cycle at the sigma_t this cadence supports is a +/-300 s bound that
    # says nothing.  The band offset is a property of the emission region,
    # not of a night, so pooling every paired cycle is the right estimator
    # and it is the only one with the signal-to-noise to detect anything.
    # Stored with night='(pooled)' so it can never be mistaken for a night.
    for tk, era, fa, fb in con.execute("""
            SELECT DISTINCT a.target_key, a.era_id, a.filter, b.filter
            FROM p3_edge a JOIN p3_edge b
              ON a.target_key=b.target_key AND a.era_id=b.era_id
             AND a.cycle=b.cycle AND a.filter < b.filter
            WHERE a.accepted=1 AND b.accepted=1
            ORDER BY 1,2,3,4""").fetchall():
        rows = con.execute("""
            SELECT a.t_edge_bjd, a.sigma_t_s, b.t_edge_bjd, b.sigma_t_s,
                   a.sigma_t_mc_s, b.sigma_t_mc_s
            FROM p3_edge a JOIN p3_edge b
              ON a.target_key=b.target_key AND a.era_id=b.era_id
             AND a.cycle=b.cycle
            WHERE a.target_key=? AND a.era_id=? AND a.filter=? AND b.filter=?
              AND a.accepted=1 AND b.accepted=1""",
                           (tk, era, fa, fb)).fetchall()
        if len(rows) < 2:
            continue
        arr = np.array([[r[0], r[1], r[2], r[3]] for r in rows], dtype=float)
        mc_a = rows[0][4] or 0.0
        mc_b = rows[0][5] or 0.0
        d, s, chi2nu = p3.band_difference(arr[:, 0],
                                          np.maximum(arr[:, 1], mc_a),
                                          arr[:, 2],
                                          np.maximum(arr[:, 3], mc_b))
        con.execute("""
            INSERT OR REPLACE INTO p3_band_pair
            (target_key, era_id, night, band_a, band_b, n_cycles, delta_s,
             sigma_s, chi2nu, significant, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)""", (
            tk, era, "(pooled)", fa, fb, len(rows), d, s, chi2nu,
            int(np.isfinite(d) and np.isfinite(s) and s > 0
                and abs(d) > BAND_PAIR_SIGMA_BAR * s),
            ("every paired cycle in this era, pooled; this is the "
             "publishable inter-band number, and the per-night rows above "
             "are its components rather than independent results")))
    con.commit()
    stamp(con, "edges")
    n_pairs = con.execute("SELECT count(*) FROM p3_band_pair").fetchone()[0]
    n_sig = con.execute("SELECT count(*) FROM p3_band_pair WHERE "
                        "significant=1").fetchone()[0]
    # THE INTER-BAND RESULT, RECORDED AS WHAT IT IS.  ``significant`` is a
    # 3-sigma flag; a stage that stores only the flag lets a reader (and a
    # previous revision of this paper) infer a detection from the mere
    # existence of non-zero differences.  The strongest pooled pair and its
    # significance are therefore stored beside the count, so the null is
    # quotable and cannot be turned into a result by paraphrase.
    top = con.execute("""
        SELECT target_key, era_id, band_a, band_b, n_cycles, delta_s, sigma_s
        FROM p3_band_pair WHERE night='(pooled)' AND sigma_s > 0
        ORDER BY abs(delta_s) / sigma_s DESC LIMIT 1""").fetchone()
    set_meta(con, {"n_edges_fitted": n_fit, "n_edges_accepted": n_acc,
                   "n_band_pairs": n_pairs, "n_band_pairs_significant": n_sig,
                   "band_pair_sigma_bar": BAND_PAIR_SIGMA_BAR,
                   "band_pair_verdict": (
                       f"{n_sig} of {n_pairs} band pairs significant at "
                       f"{BAND_PAIR_SIGMA_BAR:.0f} sigma"
                       + (f"; strongest pooled pair {top[2]}-{top[3]} "
                          f"(era {top[1]}, {top[4]} paired cycles) "
                          f"{top[5]:+.1f} +/- {top[6]:.1f} s = "
                          f"{abs(top[5]) / top[6]:.1f} sigma"
                          if top else ""))})
    con.close()
    print(f"  {n_acc}/{n_fit} edges accepted; {n_sig}/{n_pairs} inter-band "
          f"differences significant at 3 sigma")


# ===========================================================================
# STAGE: oc
# ===========================================================================
#: Circular phase scatter above which pooled edges are timing DIFFERENT
#: features and no O-C may be built from them.
ONE_FEATURE_BAR_CYCLES = 0.05

#: The strategy's per-epoch timing threshold.  CV-S5's injection grid
#: measured a median total error above it, which is why a SINGLE cycle's
#: edge is never published as an epoch and this stage publishes per-night,
#: per-band means instead.
TIMING_THRESHOLD_S = 60.0


def _era_of(series_key: str):
    """``stlmi|e76|g`` -> 76.  ``None`` when the key is not in that shape."""
    parts = str(series_key).split("|")
    if len(parts) >= 2 and parts[1].startswith("e") and parts[1][1:].isdigit():
        return int(parts[1][1:])
    return None


def injection_error_budget(con: sqlite3.Connection) -> dict:
    """The per-epoch timing error the INJECTION test demonstrated, per band.

    The formal and Monte-Carlo error bars attached to an individual edge fit
    are what the fit believes about itself; ``p3_sigmat`` is what recovering
    a synthetic edge from the real timestamps and the real residuals
    actually achieved, and the paper takes the second as governing.  It
    resolves into two pieces that behave differently under averaging:

    * ``sigma_random`` -- the scatter of the recovered epoch, which falls as
      ``1/sqrt(n)`` when n cycles of one night are averaged;
    * ``bias`` -- the systematic offset of the recovered epoch from the
      injected one, which does NOT average down and is therefore an
      irreducible floor under any night mean.

    Returns ``{band_or_'*': (sigma_random_s, bias_floor_s)}``.  The grid was
    run on one night of one target, so the band keys are that target's
    bands; ``'*'`` is the whole-grid fallback.

    THE TRANSFER THIS MAKES, STATED HERE BECAUSE THE PAPER MUST STATE IT
    ---------------------------------------------------------------------
    The grid exists for exactly one night of one target in one readout
    mode.  Looking a band up by name therefore hands a 2025 Mode0 number to
    a 2024 High Gain epoch in the same bandpass slot, and hands the
    whole-grid fallback to any band the grid never saw at all.  That
    transfer is an assumption, not a measurement.  It is made deliberately
    -- it is the only injection-measured budget this programme has -- and
    :func:`night_epochs` records which key each epoch was served from in
    ``p3_oc_night.budget_band`` so the assumption is visible per epoch, and
    the O-C stage recomputes the reduced chi-squared under the edge fits'
    own errors as the check on it.
    """
    out: dict[str, tuple[float, float]] = {}
    rows = con.execute("SELECT series_key, sigma_t_s, bias_s FROM p3_sigmat "
                       "WHERE sigma_t_s IS NOT NULL").fetchall()
    if not rows:
        return out
    by_band: dict[str, list[tuple[float, float]]] = {}
    for r in rows:
        band = str(r[0]).split("|")[-1].lower()
        by_band.setdefault(band, []).append(
            (float(r[1]), abs(float(r[2] or 0.0))))
        by_band.setdefault("*", []).append(
            (float(r[1]), abs(float(r[2] or 0.0))))
    for band, vals in by_band.items():
        out[band] = (float(np.median([v[0] for v in vals])),
                     float(np.median([v[1] for v in vals])))
    return out


def night_epochs(rows, budget: dict) -> list[dict]:
    """Collapse accepted per-cycle edges into one epoch per night per band.

    ``rows`` are the accepted ``p3_edge`` rows of one target as dicts, in
    time order, each carrying ``series_key``, ``cycle``, ``filter``,
    ``night`` and ``t_edge_bjd``.
    Each output row is the mean of the cycles timed in one band on one
    night, with the error bar the injection test licenses:

        sigma_night = sqrt( sigma_random^2 / n  +  bias^2 )

    An unweighted mean is deliberate.  Weighting by the per-edge formal
    error would let the fit's own opinion of itself decide which cycles
    count, and Section 3.1 forbids using a formal error bar at face value
    anywhere in this paper; within a single night the edges are cycles of
    the same star through the same sky, so they are of comparable quality
    by construction.

    ``n`` IS ALLOWED TO BE ONE, AND OFTEN IS.  The rule this stage enforces
    is not "an epoch must average several cycles" -- it is "no per-cycle
    edge is published carrying its own per-cycle error bar".  On a night
    with one accepted edge the epoch is that edge, with the injection
    budget attached instead of the fit's own sigma_t.  ``budget_band``
    records which key of ``budget`` served the row, so an epoch served by
    the whole-grid ``'*'`` fallback (a band the injection grid never saw)
    is distinguishable from one served by its own band, and
    ``oc_sigma_edge_s`` carries the same epoch's error propagated from the
    edge fits themselves as the check on that transfer.
    """
    groups: dict[tuple, list] = {}
    for r in rows:
        groups.setdefault((str(r["night"]), str(r["filter"])), []).append(r)
    out = []
    for (night, filt), g in sorted(groups.items()):
        n = len(g)
        key = filt.lower() if filt.lower() in budget else "*"
        s_rand, s_bias = budget.get(key, (float("nan"), float("nan")))
        sigma = float(np.hypot(s_rand / np.sqrt(n), s_bias))
        cyc = np.array([float(x["cycle"]) for x in g])
        t = np.array([float(x["t_edge_bjd"]) for x in g])
        # The same epoch's error under the edge fits' OWN errors: the
        # Monte-Carlo sigma where the fitter produced one, the rescaled
        # formal bar otherwise, propagated through an unweighted mean.
        es = [max(float(x.get("sigma_t_s") or 0.0),
                  float(x.get("sigma_t_mc_s") or 0.0)) for x in g]
        s_edge = (float(np.sqrt(np.sum(np.square(es))) / n)
                  if all(e > 0 for e in es) else None)
        out.append({
            "night": night, "filter": filt,
            "series_key": str(g[0]["series_key"]),
            "n_cycles": n, "cycle_mean": float(cyc.mean()),
            "cycle_lo": int(cyc.min()), "cycle_hi": int(cyc.max()),
            "t_mean_bjd": float(t.mean()),
            "sigma_random_s": float(s_rand / np.sqrt(n)),
            "sigma_floor_s": float(s_bias),
            "oc_sigma_s": sigma,
            "oc_sigma_edge_s": s_edge,
            "budget_band": key,
            "meets_threshold": int(sigma <= TIMING_THRESHOLD_S),
            "_cycles": cyc, "_times": t,
        })
    return out


def pdot_bound(cycles, oc_s, sigma_s, period_d: float) -> dict:
    """The period derivative the O-C residuals bound, by weighted quadratic.

    A steady ``dP/dt`` puts a term ``0.5 * (dP/dt) * P * E^2`` into an O-C
    curve.  Fitting ``a + b*(E-E0) + c*(E-E0)^2`` to the published epochs
    and reading ``c`` back gives ``dP/dt = 2c / (P * 86400)``, and the
    3-sigma bound ``|dP/dt| + 3 sigma`` is what turns "we report no period
    change" from an absence of evidence into a measured upper limit.

    Returns a dict of ``quad_coeff_s_per_cycle2``, ``quad_sigma_s_per_cycle2``,
    ``pdot``, ``pdot_sigma`` and ``pdot_limit3``, or empty when there are
    too few epochs (fewer than four) to fit a quadratic at all.
    """
    e = np.asarray(cycles, dtype=float)
    y = np.asarray(oc_s, dtype=float)
    s = np.asarray(sigma_s, dtype=float)
    ok = np.isfinite(e) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    if int(ok.sum()) < 4:
        return {}
    e, y, s = e[ok], y[ok], s[ok]
    e0 = float(e.mean())
    x = e - e0
    design = np.vstack([np.ones_like(x), x, x * x]).T
    w = 1.0 / np.square(s)
    xtwx = design.T @ (design * w[:, None])
    try:
        cov = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:                          # pragma: no cover
        return {}
    beta = cov @ (design.T @ (w * y))
    c, sc = float(beta[2]), float(np.sqrt(cov[2, 2]))
    scale = 2.0 / (float(period_d) * 86400.0)
    return {"quad_coeff_s_per_cycle2": c, "quad_sigma_s_per_cycle2": sc,
            "pdot": c * scale, "pdot_sigma": sc * scale,
            "pdot_limit3": (abs(c) + 3.0 * sc) * scale}


def cmd_oc(args) -> None:
    """O-C construction, and the cycle-count analysis that licenses it."""
    con = connect(args.db)
    ensure_tables(con)
    eph = load_ephemerides(con)
    print(f"database: {args.db}")
    con.execute("DELETE FROM p3_oc")
    con.execute("DELETE FROM p3_oc_night")
    con.execute("DELETE FROM p3_cycle_count")
    budget = injection_error_budget(con)
    for tk, e_obj in sorted(eph.items()):
        if not np.isfinite(e_obj.epoch_bjd):
            con.execute("""INSERT OR REPLACE INTO p3_cycle_count
                (target_key, period_d, sigma_period_d, sigma_basis, verdict,
                 note) VALUES (?,?,?,?,?,?)""", (
                tk, e_obj.period_d, e_obj.period_sigma_d,
                "quoted precision", "NO EPHEMERIS EPOCH",
                "VSX publishes a period but no epoch for this star, so "
                "there is no zero point to count cycles from and no O-C is "
                "possible against the catalogue value."))
            print(f"  {tk:6s} NO EPHEMERIS EPOCH")
            continue
        rows = con.execute("""
            SELECT series_key, cycle, filter, night, t_edge_bjd, sigma_t_s,
                   sigma_t_mc_s
            FROM p3_edge WHERE target_key=? AND accepted=1
            ORDER BY t_edge_bjd""", (tk,)).fetchall()
        if not rows:
            con.execute("""INSERT OR REPLACE INTO p3_cycle_count
                (target_key, epoch_bjd, period_d, sigma_period_d,
                 sigma_basis, n_epochs, verdict, note) VALUES (?,?,?,?,?,?,?,?)""",
                        (tk, e_obj.epoch_bjd, e_obj.period_d,
                         e_obj.period_sigma_d, "quoted precision", 0,
                         "NO TIMED EPOCHS",
                         "No bright-phase edge passed the acceptance gate "
                         "for this target, so there is nothing to plot."))
            print(f"  {tk:6s} no accepted edges")
            continue
        # --- is this ONE feature, or several? ---
        # An O-C pools epochs from every filter of every era into a single
        # diagram, and that is only legitimate if they are all timing the
        # SAME edge.  For ST LMi they are: the folded profile puts the
        # falling edge at phase 0.140 in all six series.  For AN UMa the
        # per-series edge phases came out at 0.802, 0.625 and 0.040, because
        # with 4 full-orbit nights out of 11 the folded profile is too poor
        # to locate the edge, and pooling those would produce an O-C of
        # three different features with an rms of a thousand seconds that
        # would read as a period error.  Measure the spread on the circle
        # and refuse to call it an O-C when it is large.
        ph_acc = np.array([r for (r,) in con.execute(
            "SELECT phase FROM p3_edge WHERE target_key=? AND accepted=1 "
            "AND phase IS NOT NULL", (tk,))], dtype=float)
        if ph_acc.size:
            vec = np.exp(2j * np.pi * ph_acc).mean()
            # Circular standard deviation, in phase units.
            r_len = float(np.abs(vec))
            phase_spread = (float(np.sqrt(-2.0 * np.log(r_len)) / (2 * np.pi))
                            if r_len > 1e-12 else 0.5)
        else:
            phase_spread = float("nan")
        one_feature = bool(np.isfinite(phase_spread)
                           and phase_spread < ONE_FEATURE_BAR_CYCLES)
        t_edge = np.array([r[4] for r in rows], dtype=float)
        sig_form = np.array([r[5] if r[5] is not None else np.nan
                             for r in rows], dtype=float)
        sig_mc = np.array([r[6] if r[6] is not None else np.nan
                           for r in rows], dtype=float)
        sigma_s = np.fmax(np.nan_to_num(sig_form, nan=0.0),
                          np.nan_to_num(sig_mc, nan=0.0))
        sigma_s = np.where(sigma_s > 0, sigma_s, np.nanmedian(sigma_s))
        sig_p = (e_obj.period_sigma_d if e_obj.period_sigma_d
                 else float("nan"))
        amb_first = p3.cycle_ambiguity(float(t_edge.min()), e_obj.period_d,
                                       e_obj.epoch_bjd, sig_p,
                                       "VSX quoted precision")
        amb_last = p3.cycle_ambiguity(float(t_edge.max()), e_obj.period_d,
                                      e_obj.epoch_bjd, sig_p,
                                      "VSX quoted precision")
        unique = bool(amb_last.unique)
        cyc = p3.cycle_number(t_edge, e_obj.period_d, e_obj.epoch_bjd)
        oc = p3.oc_seconds(t_edge, cyc, e_obj.period_d, e_obj.epoch_bjd)
        # The bright-phase edge is not at phase zero of the catalogue
        # ephemeris, so the O-C carries a constant offset that is a property
        # of the FEATURE, not of the clock.  Removing its mean is the only
        # honest way to plot it, and the mean is published beside it.
        oc_mean = float(np.mean(oc))
        for (sk, c, filt, night, te, _sf, _sm), o, s in zip(rows, oc, sigma_s):
            if not one_feature:
                # Refused above: these epochs do not time the same feature,
                # so no O-C row is written at all.  Storing them "for
                # reference" is how a refused result gets plotted anyway.
                continue
            con.execute("""INSERT OR REPLACE INTO p3_oc
                (series_key, cycle, target_key, filter, night, t_edge_bjd,
                 sigma_t_s, oc_s, oc_sigma_s, phase_offset, count_unique)
                VALUES (?,?,?,?,?,?,?,?,?,?,?)""", (
                sk, int(c), tk, filt, night, te, float(s),
                float(o - oc_mean), float(s),
                float(oc_mean / (e_obj.period_d * 86400.0)),
                int(unique)))
        fitE, fsE, fitP, fsP, chi2nu = p3.fit_linear_ephemeris(
            cyc, t_edge, sigma_s / 86400.0)

        # --- THE PUBLISHED EPOCHS: one per night per band -----------------
        # Section 4.2 concludes from the injection test that a single
        # cycle's edge does not reach the 60 s threshold on its OWN error
        # bar.  This stage therefore never publishes a per-cycle error; the
        # per-cycle rows above are the INPUTS, and what the paper plots and
        # fits is their per-night, per-band mean carrying the injection
        # budget.  On a night with one accepted edge that mean is that
        # edge -- with the budget error, not the fit's -- and the count of
        # such epochs is published beside the total.
        ne: list[dict] = []
        if one_feature:
            edge_dicts = [
                {"series_key": r[0], "cycle": r[1], "filter": r[2],
                 "night": r[3], "t_edge_bjd": r[4], "sigma_t_s": r[5],
                 "sigma_t_mc_s": r[6]} for r in rows]
            ne = night_epochs(edge_dicts, budget)
            for row in ne:
                o_night = float(np.mean(
                    p3.oc_seconds(row["_times"], row["_cycles"],
                                  e_obj.period_d, e_obj.epoch_bjd)) - oc_mean)
                within = (float(np.std(
                    p3.oc_seconds(row["_times"], row["_cycles"],
                                  e_obj.period_d, e_obj.epoch_bjd), ddof=1))
                    if row["n_cycles"] > 1 else None)
                row["oc_s"] = o_night
                row["within_night_rms_s"] = within
                con.execute("""INSERT OR REPLACE INTO p3_oc_night
                    (target_key, night, filter, series_key, era_id,
                     n_cycles, cycle_mean, cycle_lo, cycle_hi, t_mean_bjd,
                     oc_s, oc_sigma_s, sigma_random_s, sigma_floor_s,
                     within_night_rms_s, meets_threshold,
                     oc_sigma_edge_s, budget_band)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    tk, row["night"], row["filter"], row["series_key"],
                    _era_of(row["series_key"]), row["n_cycles"],
                    row["cycle_mean"], row["cycle_lo"], row["cycle_hi"],
                    row["t_mean_bjd"], o_night, row["oc_sigma_s"],
                    row["sigma_random_s"], row["sigma_floor_s"],
                    within, row["meets_threshold"],
                    row["oc_sigma_edge_s"], row["budget_band"]))
        n_single = n_span = span_d = None
        chi_edge = sg_edge_med = None
        quad: dict = {}
        if ne:
            oc_n = np.array([r["oc_s"] for r in ne], dtype=float)
            sg_n = np.array([r["oc_sigma_s"] for r in ne], dtype=float)
            w = 1.0 / np.square(sg_n)
            n_rms = float(np.sqrt(np.mean(np.square(oc_n))))
            n_wrms = float(np.sqrt(np.sum(w * np.square(oc_n)) / np.sum(w)))
            n_chi2nu = float(np.sum(np.square(oc_n / sg_n)) / oc_n.size)
            fitEn, fsEn, fitPn, fsPn, chi2n = p3.fit_linear_ephemeris(
                np.array([r["cycle_mean"] for r in ne], dtype=float),
                np.array([r["t_mean_bjd"] for r in ne], dtype=float),
                sg_n / 86400.0)
            n_at_thr = int(sum(r["meets_threshold"] for r in ne))
            n_nights = len({r["night"] for r in ne})
            sg_med, sg_lo, sg_hi = (float(np.median(sg_n)), float(sg_n.min()),
                                    float(sg_n.max()))
            # How many epochs are ONE cycle.  The rule is about the error
            # bar, not the estimator, and a reader is entitled to know how
            # often the mean has a single term in it.
            n_single = int(sum(1 for r in ne if r["n_cycles"] == 1))
            # THE SPAN OF THE EPOCHS.  Not the count from the catalogue
            # epoch (n_cycles_last), which is the extrapolation baseline
            # and is 2.5x larger here.  Both matter and they are not
            # interchangeable: the count from the catalogue epoch sets the
            # sensitivity to a period ERROR, this span the leverage on a
            # period CHANGE.
            n_span = float(amb_last.n_cycles - amb_first.n_cycles)
            span_d = float(t_edge.max() - t_edge.min())
            # The error-transfer check: the same chi-squared under the edge
            # fits' own errors instead of the one-night injection budget.
            se = [r["oc_sigma_edge_s"] for r in ne]
            if all(x is not None and np.isfinite(x) and x > 0 for x in se):
                se_a = np.array(se, dtype=float)
                chi_edge = float(np.sum(np.square(oc_n / se_a)) / oc_n.size)
                sg_edge_med = float(np.median(se_a))
            quad = pdot_bound([r["cycle_mean"] for r in ne], oc_n, sg_n,
                              e_obj.period_d)
        else:
            oc_n = np.array([])
            n_rms = n_wrms = n_chi2nu = None
            fitEn = fsEn = fitPn = fsPn = chi2n = None
            n_at_thr = n_nights = None
            sg_med = sg_lo = sg_hi = None

        ratio = (sig_p / amb_last.sigma_period_max_d
                 if np.isfinite(sig_p) and amb_last.sigma_period_max_d > 0
                 else float("nan"))
        if not one_feature:
            verdict = "NOT ONE FEATURE — NO O-C"
            note = (
                f"The {len(rows)} accepted edges for this target scatter over "
                f"{phase_spread:.3f} in orbital phase (circular s.d.), far "
                f"more than the {ONE_FEATURE_BAR_CYCLES:.2f} that would make "
                f"them one feature. "
                f"The folded profile is too sparse here to locate the falling "
                f"edge consistently between filters, so these epochs time "
                f"DIFFERENT things and pooling them into one O-C would "
                f"produce a large scatter that reads as a period error and is "
                f"not one.  The cycle count itself is still unique "
                f"({amb_last.n_cycles:,.0f} cycles, drift "
                f"{amb_last.drift_cycles:.4f}), but there is nothing "
                f"legitimate to plot against it.")
        elif unique:
            verdict = "CYCLE COUNT UNIQUE"
            note = (
                f"{amb_last.n_cycles:,.0f} cycles separate our last timed "
                f"edge from the VSX epoch.  At the quoted-precision floor "
                f"sigma_P = {sig_p:.1e} d the accumulated drift is "
                f"{amb_last.drift_cycles:.3f} cycles, below the half cycle "
                f"at which the integer count stops being determined.  The "
                f"count survives any period uncertainty up to "
                f"{amb_last.sigma_period_max_d:.2e} d, which is "
                f"{1.0 / ratio:,.0f}x the quoted precision — so the "
                f"conclusion holds even if the true published sigma is far "
                f"larger than the rounding of the last digit.  The "
                f"{len(rows)} accepted per-cycle edges are published as "
                f"{len(ne)} per-night, per-band epochs on {n_nights} nights "
                f"(p3_oc_night), spanning {n_span:,.0f} cycles "
                f"({span_d:,.0f} d) between the first and last timed edge — "
                f"which is NOT the {amb_last.n_cycles:,.0f} cycles above, "
                f"that being the extrapolation baseline from the catalogue "
                f"epoch.  The injection test does not license a single "
                f"cycle's edge as an epoch WITH ITS OWN error bar, so none "
                f"is offered as one; on {n_single} of the {len(ne)} epochs "
                f"only one cycle was timed, and those carry the injection "
                f"budget rather than the fit's own sigma_t.")
        else:
            verdict = "CYCLE COUNT NOT UNIQUE"
            note = (
                f"{amb_last.n_cycles:,.0f} cycles separate our last timed "
                f"edge from the VSX epoch, and at sigma_P = {sig_p:.1e} d "
                f"the accumulated drift is {amb_last.drift_cycles:.1f} "
                f"cycles.  The integer count is NOT determined; an O-C "
                f"built on it would be a fabricated result.  A period "
                f"uncertainty of {amb_last.sigma_period_max_d:.2e} d or "
                f"better would be needed.")
        con.execute("""INSERT OR REPLACE INTO p3_cycle_count
            (target_key, epoch_bjd, period_d, sigma_period_d, sigma_basis,
             t_first_bjd, t_last_bjd, elapsed_d, n_cycles_first,
             n_cycles_last, drift_cycles, unique_count, sigma_period_max_d,
             ratio_to_quoted, oc_mean_s, oc_rms_s, fitted_period_d,
             fitted_period_sigma_d, fitted_epoch_bjd, n_epochs,
             phase_spread, one_feature, verdict, note,
             n_night_epochs, n_nights, oc_night_rms_s, oc_night_wrms_s,
             oc_night_chi2nu, sigma_night_median_s, sigma_night_lo_s,
             sigma_night_hi_s, n_night_at_threshold,
             fitted_period_night_d, fitted_period_night_sigma_d,
             fitted_epoch_night_bjd, fit_night_chi2nu,
             n_cycles_span, span_d, n_night_single_cycle,
             oc_night_chi2nu_edge, sigma_night_edge_median_s,
             quad_coeff_s_per_cycle2, quad_sigma_s_per_cycle2,
             pdot, pdot_sigma, pdot_limit3)
             VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
                     ?,?,?,?,?,?,?,?,?,?,?,?,?,
                     ?,?,?,?,?,?,?,?,?,?)""", (
            tk, e_obj.epoch_bjd, e_obj.period_d, sig_p,
            "quoted precision of the VSX period string (VSX publishes no "
            "uncertainty; this is a FLOOR)",
            float(t_edge.min()), float(t_edge.max()), amb_last.elapsed_d,
            amb_first.n_cycles, amb_last.n_cycles, amb_last.drift_cycles,
            int(unique), amb_last.sigma_period_max_d, ratio,
            oc_mean if one_feature else None,
            (float(np.std(oc - oc_mean, ddof=1))
             if oc.size > 1 and one_feature else None),
            fitP if one_feature else None, fsP if one_feature else None,
            fitE if one_feature else None, len(rows),
            phase_spread, int(one_feature), verdict, note,
            len(ne) or None, n_nights, n_rms, n_wrms, n_chi2nu,
            sg_med, sg_lo, sg_hi, n_at_thr, fitPn, fsPn, fitEn, chi2n,
            n_span, span_d, n_single, chi_edge, sg_edge_med,
            quad.get("quad_coeff_s_per_cycle2"),
            quad.get("quad_sigma_s_per_cycle2"),
            quad.get("pdot"), quad.get("pdot_sigma"),
            quad.get("pdot_limit3")))
        print(f"  {tk:6s} {len(rows):4d} accepted edges  n_cycles="
              f"{amb_last.n_cycles:,.0f}  drift={amb_last.drift_cycles:.4f}  "
              f"phase spread {phase_spread:.3f}  {verdict}")
        if ne:
            print(f"         {len(ne)} per-night/band epochs on "
                  f"{n_nights} nights: O-C rms {n_rms:.1f}s, "
                  f"sigma_night {sg_lo:.0f}-{sg_hi:.0f}s (median "
                  f"{sg_med:.0f}s), chi2/nu about zero {n_chi2nu:.2f}, "
                  f"{n_at_thr} at or inside the "
                  f"{TIMING_THRESHOLD_S:.0f} s threshold")
            print(f"         fitted P={fitPn:.8f} +/- {fsPn:.2e} d "
                  f"(VSX {e_obj.period_d:.8f})")
    con.commit()
    stamp(con, "oc")
    n_uniq = con.execute("SELECT count(*) FROM p3_cycle_count WHERE "
                         "unique_count=1").fetchone()[0]
    n_night = con.execute("SELECT count(*) FROM p3_oc_night").fetchone()[0]
    # The provenance of every published error bar, in the release rather
    # than in the prose: the injection grid exists for ONE night of ONE
    # target, and every epoch's sigma is served from it by band slot, with
    # a whole-grid fallback for any band the grid never saw.
    src = con.execute("SELECT DISTINCT series_key, night FROM p3_sigmat "
                      "ORDER BY 1").fetchall()
    n_fb = con.execute("SELECT count(*) FROM p3_oc_night WHERE "
                       "budget_band='*'").fetchone()[0]
    set_meta(con, {"n_cycle_unique": n_uniq,
                   "n_oc_night_epochs": n_night,
                   "one_feature_bar_cycles": ONE_FEATURE_BAR_CYCLES,
                   "timing_threshold_s": TIMING_THRESHOLD_S,
                   "timing_budget_source": "; ".join(
                       f"{r[0]} on {r[1]}" for r in src) or "none",
                   "timing_budget_per_band": "; ".join(
                       f"{b}: sigma_random={v[0]:.2f} s, bias={v[1]:.2f} s"
                       for b, v in sorted(budget.items())),
                   "n_oc_night_budget_fallback": n_fb,
                   "oc_epoch_rule": (
                       "one epoch per night per band, the unweighted mean of "
                       "that night's accepted per-cycle edges, carrying the "
                       "injection-demonstrated budget sigma_random/sqrt(n) "
                       "(+) bias rather than any per-cycle error bar of its "
                       "own; n may be 1, and where it is, the epoch is that "
                       "cycle with the budget attached.  Bands are averaged "
                       "separately as a CONSERVATIVE choice: a "
                       "wavelength-dependent edge phase is expected if the "
                       "cyclotron beaming is wavelength dependent, but "
                       "p3_band_pair does not detect one (0 of 32 pairs "
                       "significant at 3 sigma, largest 1.9 sigma), so "
                       "pooling is guarded against an effect these data "
                       "could not have measured, not against a measured "
                       "one.  The budget itself is measured on a single "
                       "night of a single target and applied by band slot "
                       "across both instrument eras; "
                       "oc_night_chi2nu_edge is the check on that")})
    con.close()


# ===========================================================================
# STAGE: states
# ===========================================================================
def cmd_states(args) -> None:
    """Accretion-state classification and duty cycles, using the limits."""
    con = connect(args.db)
    ensure_tables(con)
    eph = load_ephemerides(con)
    print(f"database: {args.db}")
    con.execute("DELETE FROM p3_state_night")
    con.execute("DELETE FROM p3_state_series")
    have_limits = bool(con.execute(
        "SELECT count(*) FROM sqlite_master WHERE type='table' AND "
        "name='p2_limits'").fetchone()[0])
    for sk, tk, era, filt in series_rows(con):
        e_obj = eph.get(tk)
        data = load_series(con, sk)
        if data["t"].size < 10 or not e_obj:
            continue
        t, m = data["t"], data["m"]
        nights = np.array(data["night"])
        # --- the limits, per night, from Phase 2 ---
        lim: dict[str, list[tuple[float, float]]] = {}
        if have_limits:
            for night, lmag, lt in con.execute(
                    "SELECT night, limit_cal_mag, bjd_tdb FROM p2_limits "
                    "WHERE series_key=? AND outcome='limit' AND "
                    "limit_cal_mag IS NOT NULL", (sk,)):
                lim.setdefault(night, []).append((float(lmag), float(lt)))
        # THE PHASE-COVERAGE GATE APPLIES ONLY WHERE THERE IS AN ORBIT TO
        # GATE AGAINST.  The gate exists because these polars vary by
        # 0.65-1.7 mag around one orbit, so a median over a third of an
        # orbit is a phase measurement wearing a state label.  Where no
        # orbital modulation was detected in section 1 — every YZ Cnc
        # series, a dwarf nova outside superoutburst — there is no such
        # contamination and the gate is not just unnecessary but harmful:
        # applied blindly it discarded all 63 YZ Cnc nights, including the
        # six outburst nights the external record independently confirms.
        det_row = con.execute("SELECT detected FROM p3_period WHERE "
                              "series_key=?", (sk,)).fetchone()
        has_orbit = bool(det_row and det_row[0]) and np.isfinite(e_obj.epoch_bjd)
        per_night = []
        for night in sorted(set(nights.tolist()) | set(lim.keys())):
            sel = nights == night
            n_pts = int(sel.sum())
            lims = lim.get(night, [])
            n_lim = len(lims)
            # Coverage is measured over EVERY epoch of the night, detected
            # or not: a night of pure non-detections still sampled the orbit
            # and its limits still constrain it.
            all_t = np.concatenate([t[sel], np.array([x[1] for x in lims])]) \
                if (n_pts or n_lim) else np.array([])
            cov = (p3.phase_coverage(all_t, e_obj.period_d, e_obj.epoch_bjd)
                   if has_orbit and all_t.size else float("nan"))
            gated = int(has_orbit and np.isfinite(cov)
                        and cov < p3.STATE_MIN_PHASE_COVERAGE)
            gate_note = ("phase coverage below "
                         f"{p3.STATE_MIN_PHASE_COVERAGE:.0%} — this night's "
                         "median is a phase measurement, not a state"
                         ) if gated else ""
            if n_pts == 0:
                if not n_lim:
                    continue
                # Nothing detected all night.  The night's magnitude is the
                # median 3-sigma limit and is a lower bound on faintness.
                per_night.append((
                    night, 0, cov, float(np.median([x[0] for x in lims])),
                    None, None, None, n_lim, 1, gated,
                    ("no detection all night; the night's magnitude is the "
                     "median 3-sigma upper limit, so the star was AT LEAST "
                     "this faint. " + gate_note).strip()))
                continue
            # Detections AND limits: the night's magnitude is the
            # Kaplan-Meier median over both, so the undetected epochs pull
            # it faint instead of being dropped.  The naive median over
            # detections alone is kept in p10/p90 so the size of that pull
            # stays visible per night rather than only in the aggregate.
            det = m[sel]
            note = gate_note
            censored = 0
            if n_lim:
                vals = np.concatenate([det, np.array([x[0] for x in lims])])
                cens = np.concatenate([np.zeros(det.size, bool),
                                       np.ones(n_lim, bool)])
                km = p2_km_median(vals, cens)
                if np.isfinite(km):
                    mag = float(km)
                    note = (f"{n_lim} of {n_pts + n_lim} epochs are upper "
                            f"limits; magnitude is the Kaplan-Meier median "
                            f"({km:.3f}) rather than the detections-only "
                            f"median ({np.median(det):.3f}). " + note).strip()
                else:
                    mag = float(np.median(det))
                    censored = 1
                    note = (f"{n_lim} of {n_pts + n_lim} epochs are upper "
                            "limits and the Kaplan-Meier curve never reaches "
                            "0.5, so the true median is fainter than any "
                            "value this night can name. " + note).strip()
            else:
                mag = float(np.median(det))
            per_night.append((
                night, n_pts, cov, mag,
                float(np.percentile(det, 10)), float(np.percentile(det, 90)),
                float(np.percentile(det, 90) - np.percentile(det, 10)),
                n_lim, censored, gated, note))
        if not per_night:
            continue
        # --- the threshold, from the ungated nights' own bimodality ---
        usable = [r for r in per_night if not r[9] and not r[8]]
        vals = np.array([r[3] for r in usable], dtype=float)
        thr, sep = p3.otsu_threshold(vals) if vals.size >= 4 else (float("nan"),) * 2
        thr_sig = (p3.bootstrap_threshold(vals, n_boot=args.n_boot)
                   if vals.size >= 6 else float("nan"))
        bimodal = int(np.isfinite(sep) and sep > 0.75)
        states = {}
        for r in per_night:
            night = r[0]
            st = ("UNCLASSIFIED" if r[9] else
                  p3.classify_state(r[3], thr, thr_sig))
            if r[8] and st == "HIGH":
                # A limit can never prove a HIGH state: "fainter than X" is
                # not evidence of brightness.
                st = "UNCLASSIFIED"
            states[night] = st
            con.execute("""INSERT OR REPLACE INTO p3_state_night
                (series_key, night, target_key, era_id, filter, n_points,
                 phase_coverage, median_mag, p10_mag, p90_mag,
                 amplitude_mag, n_limits, censored, gated, state, note)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                        (sk, r[0], tk, era, filt, r[1], r[2], r[3], r[4],
                         r[5], r[6], r[7], r[8], r[9], st, r[10]))
        duty = p3.duty_cycle([r[3] for r in per_night if not r[9]],
                             [bool(r[8]) for r in per_night if not r[9]],
                             thr)
        n_high = sum(1 for s in states.values() if s == "HIGH")
        n_low = sum(1 for s in states.values() if s == "LOW")
        n_int = sum(1 for s in states.values() if s == "INTERMEDIATE")
        if not np.isfinite(thr):
            verdict = "NO THRESHOLD"
            note = (f"only {vals.size} ungated nights — too few for a "
                    "threshold, so no state history is published for this "
                    "series")
        elif not bimodal:
            verdict = "UNIMODAL — NO STATE CHANGE RESOLVED"
            note = (f"Otsu separability {sep:.2f} below 0.75: the nightly "
                    "magnitudes are one population, so the threshold cuts a "
                    "single distribution in half and the HIGH/LOW labels "
                    "below are NOT evidence of two accretion states")
        else:
            verdict = "BIMODAL — TWO STATES RESOLVED"
            note = (f"Otsu separability {sep:.2f}: two distinct populations "
                    f"split at {thr:.2f} mag +/- {thr_sig:.2f}")
        con.execute("""INSERT OR REPLACE INTO p3_state_series
            (series_key, target_key, era_id, filter, n_nights, n_gated,
             n_used, threshold_mag, threshold_sigma, separability, bimodal,
             n_high, n_low, n_intermediate, n_censored, duty_naive,
             duty_with_limits, duty_bias, n_informative_limits,
             n_uninformative, span_mag, verdict, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
            sk, tk, era, filt, len(per_night),
            sum(1 for r in per_night if r[9]), int(vals.size),
            thr if np.isfinite(thr) else None,
            thr_sig if np.isfinite(thr_sig) else None,
            sep if np.isfinite(sep) else None, bimodal,
            n_high, n_low, n_int, sum(1 for r in per_night if r[8]),
            duty["naive"] if np.isfinite(duty["naive"]) else None,
            duty["with_limits"] if np.isfinite(duty["with_limits"]) else None,
            duty["bias"] if np.isfinite(duty["bias"]) else None,
            duty["n_informative_limits"], duty["n_uninformative"],
            float(np.ptp(vals)) if vals.size else None, verdict, note))
        print(f"  {sk:16s} {len(per_night):3d} nights ({vals.size} used)  "
              f"thr={thr:.2f}+/-{thr_sig:.2f} sep={sep:.2f}  "
              f"H/L/I={n_high}/{n_low}/{n_int}  "
              f"duty {duty['naive']:.2f} -> {duty['with_limits']:.2f}  "
              f"{verdict}")
    con.commit()
    stamp(con, "states")
    set_meta(con, {"state_min_phase_coverage": p3.STATE_MIN_PHASE_COVERAGE,
                   "state_bimodal_bar": 0.75})
    con.close()


# ===========================================================================
# STAGE: detrend
# ===========================================================================
def cmd_detrend(args) -> None:
    """Joint GP + signal fitting, and why the order matters."""
    con = connect(args.db)
    ensure_tables(con)
    eph = load_ephemerides(con)
    print(f"database: {args.db}")
    # The demonstration runs on the densest night of the series with the
    # most target points, so the sampling is real and the injected signal
    # has several cycles to live in.
    row = con.execute("""
        SELECT f.series_key, f.night, count(*) n FROM cv_frames f
        JOIN cv_series s ON s.series_key = f.series_key
        WHERE s.status='solved' GROUP BY 1,2 ORDER BY n DESC LIMIT 1
        """).fetchone()
    if not row:
        print("  ! no series to demonstrate on")
        con.close()
        return
    sk, night = row[0], row[1]
    tk = sk.split("|")[0]
    e_obj = eph.get(tk)
    rows = con.execute("""
        SELECT l.bjd_tdb, l.cal_mag, l.inst_mag_err FROM cv_lightcurve l
        JOIN cv_frames f ON f.frame_id=l.frame_id AND f.series_key=l.series_key
        LEFT JOIN p2_cloud_frame c ON c.frame_id=l.frame_id
                                  AND c.series_key=l.series_key
        WHERE l.series_key=? AND l.role='target' AND f.night=?
          AND l.cal_mag IS NOT NULL AND COALESCE(c.vetoed,0)=0
        ORDER BY 1""", (sk, night)).fetchall()
    a = np.array(rows, dtype=float)
    if a.shape[0] < 40:
        print(f"  ! {sk} {night} has only {a.shape[0]} points")
        con.close()
        return
    t = a[:, 0]
    infl = inflation_for(con, sk)
    e = np.where(np.isfinite(a[:, 2]) & (a[:, 2] > 0), a[:, 2],
                 np.nanmedian(a[:, 2])) * infl
    period = e_obj.period_d if e_obj else 0.08
    freq = 1.0 / period
    print(f"  demonstrating on {sk} night {night}: {t.size} points, "
          f"{(t.max() - t.min()) * 24:.2f} h, {(t.max() - t.min()) / period:.2f} "
          f"orbits, sigma {np.median(e) * 1000:.1f} mmag")
    # A REAL trend to hide the signal under: a slow smooth drift of the same
    # size as the signal.  Drawn from a GP so it is a fair example of what
    # the joint fit is designed to absorb, and stated so nobody thinks the
    # demonstration was rigged with something the smoother could never track.
    rng = np.random.default_rng(p3.SEED)
    k = p3.matern32_cov(t, 0.15, 0.12)
    k[np.diag_indices_from(k)] += 1e-10
    trend = np.linalg.cholesky(k) @ rng.normal(size=t.size)
    # --- the celerite2 / dense cross-check, on the real sampling ---
    y_check = 16.0 + trend + DETREND_AMPLITUDE * np.sin(
        2 * np.pi * freq * t) + rng.normal(0.0, e)
    fast = p3.joint_gp_fit(t, y_check, e, freq, use_celerite=True)
    slow = p3.joint_gp_fit(t, y_check, e, freq, use_celerite=False)
    rel = (abs(fast["amplitude"] - slow["amplitude"])
           / max(abs(slow["amplitude"]), 1e-12))
    con.execute("""INSERT OR REPLACE INTO p3_gp_check
        (series_key, n_points, amp_celerite, amp_dense, rel_diff,
         loglike_celerite, loglike_dense, ll_abs_diff, celerite_eps, verdict)
        VALUES (?,?,?,?,?,?,?,?,?,?)""", (
        sk, int(t.size), fast["amplitude"], slow["amplitude"], rel,
        fast["loglike"], slow["loglike"],
        abs(fast["loglike"] - slow["loglike"]), p3.CELERITE_MATERN_EPS,
        "AGREE" if rel < 1e-6 else "DISAGREE"))
    print(f"  GP backend check: celerite2 A={fast['amplitude']:.6f} "
          f"dense A={slow['amplitude']:.6f} rel diff {rel:.2e} "
          f"({fast['backend']} vs {slow['backend']}, "
          f"eps={p3.CELERITE_MATERN_EPS:g})")
    con.execute("DELETE FROM p3_detrend WHERE series_key=?", (sk,))
    for wp in DETREND_WINDOWS_P:
        out = p3.detrend_suppression(t, e, freq, DETREND_AMPLITUDE,
                                     window_d=wp * period, trend=trend,
                                     n_real=args.n_real_detrend)
        con.execute("""INSERT OR REPLACE INTO p3_detrend
            (series_key, window_periods, window_d, amplitude_in,
             frac_detrend, frac_joint, n_detrend, n_joint, backend)
            VALUES (?,?,?,?,?,?,?,?,?)""", (
            sk, float(wp), float(wp * period), DETREND_AMPLITUDE,
            out["frac_detrend"], out["frac_joint"], out["n_detrend"],
            out["n_joint"], fast["backend"]))
        print(f"    window {wp:5.1f} P  detrend-then-search recovers "
              f"{out['frac_detrend']:6.1%} of the injected amplitude, "
              f"joint fit {out['frac_joint']:6.1%}")
    con.commit()
    stamp(con, "detrend")
    set_meta(con, {"detrend_series": sk, "detrend_night": night,
                   "detrend_amplitude_mag": DETREND_AMPLITUDE,
                   "detrend_reals": args.n_real_detrend,
                   "celerite_eps": p3.CELERITE_MATERN_EPS})
    con.close()


# ===========================================================================
# STAGE: report
# ===========================================================================
def cmd_report(args) -> None:
    from macro_phot import report_phase3
    path = report_phase3.render_report(args.db)
    print(f"  wrote {path}")
    record_stage("R-CV-S9")


# ===========================================================================
# STAGE: status
# ===========================================================================
def cmd_status(args) -> None:
    """Where every stage stands, and what it currently concludes."""
    con = connect(args.db, read_only=True)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    print(f"database: {args.db}")
    for t in ("p3_ephemeris", "p3_period", "p3_pgram", "p3_sigmat",
              "p3_edge", "p3_band_pair", "p3_oc", "p3_cycle_count",
              "p3_state_night", "p3_state_series", "p3_detrend"):
        n = con.execute(f"SELECT count(*) FROM {t}").fetchone()[0] \
            if t in have else 0
        print(f"  {t:20s} {n:>10,}")
    if "p3_meta" not in have:
        print("  (nothing run yet)")
        con.close()
        return
    meta = dict(con.execute("SELECT key, value FROM p3_meta"))
    print(f"\n  code version: {meta.get('phase3_code_version', '-')}  "
          f"commit {meta.get('git_commit', '-')}")
    if "p3_period" in have:
        print("\n  PERIOD VERIFICATION (recovered vs published)")
        for r in con.execute("""
                SELECT series_key, n_points, n_blocks, period_d,
                       sigma_period_d, published_d, deviation_sigma, agrees,
                       alias_frac_max, family_code, detected, harmonic_note,
                       constraint_class
                FROM p3_period WHERE status='ok'
                ORDER BY target_key, era_id, filter"""):
            if r[3] is None:
                print(f"    {r[0]:16s} no period")
                continue
            mark = "-- " if r[7] is None else ("OK " if r[7] else "XX ")
            print(f"    {mark}{r[0]:16s} n={r[1]:4d} nights={r[2]:3d}  "
                  f"P={r[3]:.7f}+/-{(r[4] or 0):.1e}  "
                  f"pub={r[5]:.7f}  dev={(r[6] or 0):+6.2f}s  "
                  f"alias={r[8]:.2f}  {r[9]:12s} {(r[12] or ''):14s} "
                  f"{'DETECTED' if r[10] else 'not detected'}")
            if r[11]:
                print(f"        ! {r[11]}")
    if "p3_sigmat" in have:
        print("\n  SIGMA_T INJECTION (ST LMi dense night, threshold "
              f"{meta.get('sigmat_threshold_s', '60')} s)")
        for sk in [x[0] for x in con.execute(
                "SELECT DISTINCT series_key FROM p3_sigmat ORDER BY 1")]:
            print(f"    {sk}: {meta.get('sigmat_verdict_' + sk, '-')}")
            print(f"      width scan: {meta.get('sigmat_width_scan_' + sk, '-')}")
            for r in con.execute(
                    "SELECT shape_error, depth_error, sigma_t_s, bias_s, "
                    "total_error_s, passes, inject_width_s FROM p3_sigmat "
                    "WHERE series_key=? AND inject_factor=1.0 "
                    "ORDER BY shape_error, depth_error", (sk,)):
                print(f"      inject {r[6]:4.0f}s shape x{r[0]:.0f} "
                      f"depth {r[1]:+4.0%}  sigma={r[2]:6.1f} "
                      f"bias={r[3]:+6.1f} total={r[4]:6.1f}s  "
                      f"{'PASS' if r[5] else 'FAIL'}")
    if "p3_band_pair" in have:
        print("\n  INTER-BAND EDGE TIMING")
        for r in con.execute(
                "SELECT target_key, era_id, night, band_a, band_b, n_cycles,"
                " delta_s, sigma_s, chi2nu, significant FROM p3_band_pair "
                "ORDER BY target_key, night, band_a"):
            print(f"    {r[0]:6s} e{r[1]:<3d} {r[2]}  {r[3]}-{r[4]}  "
                  f"n={r[5]:2d}  {r[6]:+7.1f} +/- {r[7]:5.1f} s  "
                  f"chi2nu={r[8]:5.2f}  "
                  f"{'SIGNIFICANT' if r[9] else 'consistent with zero'}")
    if "p3_cycle_count" in have:
        print("\n  CYCLE COUNT")
        for r in con.execute(
                "SELECT target_key, n_epochs, n_cycles_last, drift_cycles, "
                "unique_count, verdict, fitted_period_d, "
                "fitted_period_sigma_d, period_d, oc_rms_s "
                "FROM p3_cycle_count ORDER BY target_key"):
            print(f"    {r[0]:6s} {r[5]}")
            if r[2] is not None:
                print(f"      {r[1]} epochs, {r[2]:,.0f} cycles from the VSX "
                      f"epoch, drift {r[3]:.4f} cycles, O-C rms "
                      f"{(r[9] or 0):.1f} s")
                print(f"      fitted P = {r[6]:.8f} +/- {r[7]:.2e} d  "
                      f"(VSX {r[8]:.8f})")
    if "p3_state_series" in have:
        print("\n  ACCRETION STATES (duty cycle = fraction of epochs HIGH)")
        for r in con.execute(
                "SELECT series_key, n_nights, n_used, threshold_mag, "
                "separability, n_high, n_low, n_intermediate, duty_naive, "
                "duty_with_limits, verdict FROM p3_state_series "
                "ORDER BY target_key, era_id, filter"):
            thr = f"{r[3]:.2f}" if r[3] is not None else "  -  "
            sep = f"{r[4]:.2f}" if r[4] is not None else " -  "
            dn = f"{r[8]:.3f}" if r[8] is not None else "  -  "
            dl = f"{r[9]:.3f}" if r[9] is not None else "  -  "
            print(f"    {r[0]:16s} {r[1]:3d} nights  thr={thr} sep={sep}  "
                  f"H/L/I={r[5]}/{r[6]}/{r[7]}  duty {dn} -> {dl}  {r[10]}")
    if "p3_detrend" in have:
        print("\n  DETRENDING DISCIPLINE (fraction of an injected amplitude "
              "recovered)")
        for r in con.execute(
                "SELECT window_periods, frac_detrend, frac_joint FROM "
                "p3_detrend ORDER BY window_periods"):
            print(f"    window {r[0]:5.1f} P   detrend-then-search "
                  f"{r[1]:7.1%}   joint GP+signal {r[2]:7.1%}")
    con.close()


# ===========================================================================
# CLI
# ===========================================================================
def _common(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Every option, on every parser.

    Declared through a shared parent rather than only on the top-level
    parser so that BOTH ``run_cv_phase3.py --workers 6 periods`` and
    ``run_cv_phase3.py periods --workers 6`` work.  Argparse's default
    behaviour accepts only the first, and a stage that silently ran with
    400 bootstrap replicates because the flag landed after the subcommand
    and was rejected is exactly the kind of thing that makes a published
    number unreproducible.
    """
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--workers", type=int, default=MAX_WORKERS,
                        help=f"worker processes for `periods` (hard cap "
                             f"{MAX_WORKERS}; this machine is also running "
                             f"an S1 batch and a catalogue-tie re-run)")
    parser.add_argument("--n-boot", type=int, default=p3.N_BOOT,
                        help="night-bootstrap replicates for the period and "
                             "state-threshold uncertainties")
    parser.add_argument("--n-real", type=int, default=N_SIGMAT_REAL,
                        help="Monte Carlo realizations per sigma_t grid cell")
    parser.add_argument("--n-real-detrend", type=int, default=DETREND_REALS,
                        help="realizations per detrend-demonstration window")
    parser.add_argument("--force", action="store_true",
                        help="refetch the ephemerides instead of using the "
                             "cache")
    return parser


def main() -> None:
    p = _common(argparse.ArgumentParser(description=__doc__.split("\n")[0]))
    parent = _common(argparse.ArgumentParser(add_help=False))
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("ephem", "periods", "sigmat", "edges", "oc", "states",
                 "detrend", "report", "status", "all"):
        sub.add_parser(name, parents=[parent])
    args = p.parse_args()
    table = {"ephem": cmd_ephem, "periods": cmd_periods,
             "sigmat": cmd_sigmat, "edges": cmd_edges, "oc": cmd_oc,
             "states": cmd_states, "detrend": cmd_detrend,
             "report": cmd_report, "status": cmd_status}
    order = ("ephem", "periods", "sigmat", "edges", "oc", "states",
             "detrend", "report")
    if args.cmd == "all":
        for name in order:
            print(f"\n=== {name} ===", flush=True)
            table[name](args)
        record_stage("CV-S9")
        cmd_status(args)
    else:
        table[args.cmd](args)
        if args.cmd in ("periods", "sigmat", "edges", "oc", "states",
                        "detrend"):
            record_stage("CV-S9")


if __name__ == "__main__":
    main()
