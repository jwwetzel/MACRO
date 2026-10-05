#!/usr/bin/env python
"""CV-R1...R8 — the statistics of the CV paper, redone for the committee.

WHAT THIS SCRIPT DOES, AND WHY EACH PIECE IS HERE
--------------------------------------------------
The plan review of 2026-10-03 reopened the CV paper as a major revision
(``committee/reviews/2026-10-03/SYNTHESIS.md`` section 4).  Eight of its
fifteen tasks are statistics, and this script is those eight.  Every stage
below answers a finding by test and stores the answer in the products
database, so the manuscript package that follows can quote a macro rather
than an argument.

``edges``       (R2, R7)  Re-measures every accepted ST LMi edge with two
                further estimators beside the published one, and attaches
                to every edge a real per-edge bootstrap error and its
                grid-quantisation term.  The published v1 estimator fits a
                ramp in MAGNITUDES on a width grid taken from each band's
                own folded profile; ``magc`` fits the same magnitudes on
                one grid common to all bands; ``v2`` fits RELATIVE FLUX on
                that common grid.  The three exist so that "is the band
                offset made by the estimator" is answered by changing the
                estimator and looking.

``inject``      (R2, D3)  Injects ONE achromatic flux edge into the real
                timestamps of every band of every night that contributes
                an epoch — 2024 High Gain nights included — with each
                band's own depth, its own errors and its own real
                residuals, and recovers it with the estimator that band
                actually used.  Stores the SIGNED bias per band and the
                SIGNED differential bias per band pair, with standard
                errors, for every true ramp width tested.  No median of
                absolute values, no transported budget.

``offset``      (R1, D3)  The band-offset question as a paired test:
                same-cycle and same-night differences, scatter-based
                errors, sign / signed-rank / sign-flip permutation tests,
                a night-clustered permutation, both eras separately and
                combined, a three-pair trials factor, a non-parametric
                half-flux crossing as an estimator-independent check, and
                the injection's differential bias set beside the result.

``oc``          (R3)  The O-C refit: per-band constants, an era-offset
                nuisance term, night-level epochs, 2024 dropped, each with
                chi-squared AND its degrees of freedom, and the
                period-derivative bound quoted every way.

``longitude``   (R5)  The timed feature is an accretion-spot edge, so the
                same residuals are a spot LONGITUDE.  Stability in
                degrees, split by accretion state from ``p3_state_night``,
                and the three physical scales a period-derivative bound
                has to be read against.

``superhump``   (R4)  What YZ Cnc's blind-search contours do and do not
                exclude, as counts of run-filters per amplitude.

``yzrefold``    (R4)  Refolds every YZ Cnc scope on the one published
                period that has an error bar (van Paradijs et al. 1994)
                and restates the between-night phase drift under each
                published period error.

``colour``      (R6)  The colour curves as measurements: amplitude, phase
                of the extrema, night-bootstrap errors, per-era
                repeatability, and survival under a 120 s pairing window
                and under an interpolating estimator that is immune to
                non-simultaneity altogether.

``fitquality``  (R7)  The per-edge chi-squared distribution, with degrees
                of freedom, and the three fits the paper should show.

``counts``      (R8)  Staged versus measured frames; what an S/N floor
                does to "measurement"; EU UMa i and r; the FWHM > 5 arcsec
                cut and what it changes.

``figures``     Draws the revision figures from the tables above.

``all``         Every stage, in the order above.

``status``      Prints what is in the database.

All the arithmetic lives in ``macro_phot.revision_cv`` and is unit-tested
in ``pipeline/tests/test_revision_cv.py``.  This script is I/O, staging,
parallelism and bookkeeping.

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_revision.py all
    $P pipeline/scripts/run_cv_revision.py inject --workers 6
    $P pipeline/scripts/run_cv_revision.py offset
    $P pipeline/scripts/run_cv_revision.py figures
    $P pipeline/scripts/run_cv_revision.py status

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite, prefix rv_)
--------------------------------------------------------------------------
``rv_result``        every SCALAR the manuscript may quote: key, value,
                     unit, format, stage and the clause a referee needs.
                     ``numbers_cv`` turns each row into a macro.
``rv_edge``          one row per (edge, estimator): epoch, width, chi2,
                     dof, valley width, grid step, quantisation sigma and
                     the per-edge bootstrap sigma.
``rv_inject_band``   signed injection bias per (night, band, estimator,
                     noise model, true width).
``rv_inject_pair``   signed differential bias per (night, band pair,
                     estimator, noise model, true width a, true width b).
``rv_inject_summary``the same, combined over nights with the weights the
                     real measurement has (pairs or edges per night).
``rv_band_offset``   the paired band-offset tests.
``rv_halfflux``      the non-parametric half-flux crossing per band/era.
``rv_stack``         the cycle-aligned, level-normalised egress profile
                     per band: the picture of the band offset.
``rv_oc_fit``        one row per O-C model variant.
``rv_oc_epoch``      the night-level epochs (N = 17).
``rv_oc_chi2``       chi-squared per band and per era for each variant.
``rv_longitude``     spot-longitude statistics by accretion state.
``rv_superhump``     excluded-amplitude counts per threshold.
``rv_supercycle``    where each YZ Cnc run sits relative to the
                     superoutbursts before and after it.
``rv_yz_refold``     every YZ Cnc scope folded on both published periods.
``rv_colour_curve``  the binned colour curves, per era/colour/method.
``rv_colour_summary``amplitude, extremum phases and errors.
``rv_colour_repeat`` era-to-era repeatability of each curve's shape.
``rv_edge_quality``  the chi-squared distribution summary per band.
``rv_edge_example``  the points of the best, the median and the worst
                     published edge fit, so the figure that shows them is
                     drawn from the database like every other.
``rv_counts``        frame and measurement counts for R8.
``rv_figure``        one row per revision figure drawn.
``rv_meta``          build stamps and every constant this run used.

WHAT THIS SCRIPT DOES NOT TOUCH
-------------------------------
No ``p3_``, ``p4_`` or ``p5_`` table is written.  The published tables are
the record of what the first draft measured, and the revision is defined
against them.  ``manuscripts/CV_TimeSeries/main.tex`` is not edited here
or anywhere in this package.

CONCURRENCY
-----------
Other packages read this database while this runs.  Every connection sets
``busy_timeout = 300000``; transactions are short and per-stage; the
worker count is capped at 6.
"""

from __future__ import annotations

import argparse
import math
import sqlite3
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "scripts"))

from macro_phot import phase3 as p3                             # noqa: E402
from macro_phot import revision_cv as rv                        # noqa: E402
# The published edge stage's own loaders and its own profile measurement
# are IMPORTED, not re-implemented: the revision must see exactly the
# points, the edge phase and the width grid the published fits saw, and
# the only way to guarantee that is to call the same code.
import run_cv_phase3 as s9                                      # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
CHAR_DB = REPO_ROOT / "products" / "phot" / "cv_characterization.sqlite"
PDF_DIR = REPO_ROOT / "manuscripts" / "CV_TimeSeries" / "figures" / "revision"
PNG_DIR = REPO_ROOT / "docs" / "CV_TimeSeries" / "figures" / "cv_revision"

REVISION_CODE_VERSION = ("CV-R v1.0 (2026-10-03, committee major revision: "
                         "R1-R8 statistics)")

MAX_WORKERS = 6
BUSY_TIMEOUT_MS = 300_000

#: The one target with an O-C, and therefore the one this revision's
#: timing stages act on.
TARGET = "stlmi"

#: Instrument eras of the ST LMi timing series and what they are.
ERA_NAME = {7: "2024 High Gain", 47: "2024 1 MHz", 76: "2025-26 Mode0"}

#: Reference nights for the Gaussian-noise cross-check of the injection:
#: the densest night of each camera era.
GAUSS_NIGHTS = ("2025-02-27", "2024-03-03")

#: Amplitude thresholds (semi-amplitude, mag) the superhump statement is
#: tabulated at.  The first is the "floor" the first draft used, for which
#: the literature package found no source; 125 and 150 mmag are the
#: published PEAK superhump semi-amplitudes (smak2010, kato2012; dai2026
#: for YZ Cnc itself in TESS) and are the ones the statement is made at.
SUPERHUMP_THRESHOLDS = (0.050, 0.100, rv.SUPERHUMP_SEMI_PEAK_MAG,
                        rv.SUPERHUMP_SEMI_TESS_MAG, 0.200, 0.250)

#: S/N floors tabulated for "measurement" (R8).  5 is the detection
#: threshold the photometry itself uses (``cv_build_meta.detect_sigma``).
SNR_FLOORS = (3.0, 5.0)

#: Fewest pairs an era must hold before its scatter-based error is trusted
#: as an inverse-variance weight.
MIN_IV_PAIRS = 5

#: Seeing cut tabulated for R8, arcsec.  The photometric aperture radius is
#: 4 arcsec, so at FWHM > 5 arcsec more than a third of a Gaussian star's
#: light falls outside it.
FWHM_CUT_ARCSEC = 5.0


# ===========================================================================
# Plumbing
# ===========================================================================
def connect(path: Path, read_only: bool = False) -> sqlite3.Connection:
    """Open the products database with a long busy timeout."""
    if read_only:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True,
                              timeout=BUSY_TIMEOUT_MS / 1000)
    else:
        con = sqlite3.connect(path, timeout=BUSY_TIMEOUT_MS / 1000)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    con.row_factory = sqlite3.Row
    return con


def git_commit() -> str:
    """Short commit hash, with ``-dirty`` when the tree has local edits."""
    try:
        h = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                           cwd=REPO_ROOT, capture_output=True, text=True,
                           timeout=30).stdout.strip()
        d = subprocess.run(["git", "status", "--porcelain",
                            "--untracked-files=no"], cwd=REPO_ROOT,
                           capture_output=True, text=True,
                           timeout=60).stdout.strip()
        return (h + "-dirty") if d else h
    except Exception:                                       # noqa: BLE001
        return "unknown"


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def ensure_tables(con: sqlite3.Connection) -> None:
    """Create every ``rv_`` table if it is not there yet."""
    con.executescript("""
    CREATE TABLE IF NOT EXISTS rv_meta (key TEXT PRIMARY KEY, value TEXT);
    CREATE TABLE IF NOT EXISTS rv_result (
        key TEXT PRIMARY KEY, value REAL, text TEXT, unit TEXT, fmt TEXT,
        stage TEXT, origin TEXT, note TEXT);
    CREATE TABLE IF NOT EXISTS rv_edge (
        series_key TEXT, cycle INTEGER, estimator TEXT, target_key TEXT,
        era_id INTEGER, filter TEXT, band TEXT, night TEXT,
        t_edge_bjd REAL, oc_raw_s REAL, width_s REAL, level_bright REAL,
        step REAL, snr REAL,
        chi2 REAL, dof INTEGER, chi2nu REAL, n_points INTEGER,
        bracket_s REAL, valley_s REAL, grid_step_s REAL,
        sigma_quant_s REAL, sigma_boot_s REAL, boot_n_ok INTEGER,
        sigma_edge_s REAL, accepted INTEGER, reason TEXT,
        v1_refit_minus_stored_s REAL, n_fwhm_gt_cut INTEGER,
        PRIMARY KEY (series_key, cycle, estimator));
    CREATE TABLE IF NOT EXISTS rv_inject_band (
        night TEXT, era_id INTEGER, series_key TEXT, band TEXT,
        estimator TEXT, noise TEXT, true_width_s REAL, depth_mag REAL,
        n_try INTEGER, n_ok INTEGER, bias_mean_s REAL, bias_se_s REAL,
        bias_median_s REAL, sigma_s REAL, rms_s REAL, pool_basis TEXT,
        pool_n INTEGER, n_real_edges INTEGER,
        PRIMARY KEY (night, band, estimator, noise, true_width_s));
    CREATE TABLE IF NOT EXISTS rv_inject_pair (
        night TEXT, era_id INTEGER, band_a TEXT, band_b TEXT,
        estimator TEXT, noise TEXT, true_width_a_s REAL,
        true_width_b_s REAL, n_ok INTEGER, dbias_mean_s REAL,
        dbias_se_s REAL, dbias_sd_s REAL, n_real_pairs INTEGER,
        PRIMARY KEY (night, band_a, band_b, estimator, noise,
                     true_width_a_s, true_width_b_s));
    CREATE TABLE IF NOT EXISTS rv_inject_summary (
        kind TEXT, era TEXT, estimator TEXT, noise TEXT, band_a TEXT,
        band_b TEXT, true_width_a_s REAL, true_width_b_s REAL,
        n_nights INTEGER, weight_sum REAL, bias_s REAL, bias_se_s REAL,
        sigma_s REAL, matched INTEGER,
        PRIMARY KEY (kind, era, estimator, noise, band_a, band_b,
                     true_width_a_s, true_width_b_s));
    CREATE TABLE IF NOT EXISTS rv_band_offset (
        estimator TEXT, level TEXT, era TEXT, band_a TEXT, band_b TEXT,
        n INTEGER, n_nights INTEGER, mean_s REAL, median_s REAL, sd_s REAL,
        se_s REAL, t REAL, p_t REAL, n_negative INTEGER, p_sign REAL,
        p_wilcoxon REAL, p_perm REAL, perm_basis TEXT, p_perm_cluster REAL,
        p_perm_bonf REAL, p_perm_cluster_bonf REAL, p_wilcoxon_bonf REAL,
        mean_pub_s REAL, sigma_pub_s REAL, chi2_pub REAL, dof_pub INTEGER,
        chi2nu_pub REAL, mean_boot_s REAL, sigma_boot_s REAL,
        chi2_boot REAL, dof_boot INTEGER, chi2nu_boot REAL,
        chi2_het REAL, dof_het INTEGER, p_het REAL, note TEXT,
        PRIMARY KEY (estimator, level, era, band_a, band_b));
    CREATE TABLE IF NOT EXISTS rv_halfflux (
        era TEXT, band TEXT, n_nights INTEGER, n_points INTEGER,
        t_half_s REAL, t_half_err_s REAL, depth_frac REAL, note TEXT,
        PRIMARY KEY (era, band));
    CREATE TABLE IF NOT EXISTS rv_stack (
        scope TEXT, band TEXT, bin INTEGER, dt_s REAL, level REAL,
        err REAL, n INTEGER, PRIMARY KEY (scope, band, bin));
    CREATE TABLE IF NOT EXISTS rv_oc_fit (
        variant TEXT, model TEXT, estimator TEXT, sigma_model TEXT,
        n_epochs INTEGER, n_informative INTEGER, n_params INTEGER,
        chi2 REAL, dof INTEGER, chi2nu REAL, rms_s REAL,
        gamma_s_per_cycle2 REAL, gamma_sigma_budget REAL,
        gamma_sigma_scatter REAL, pdot REAL, pdot_sigma_budget REAL,
        pdot_sigma_scatter REAL, pdot_limit3_budget REAL,
        pdot_limit3_scatter REAL, pdot_nsigma_scatter REAL,
        beta_s_per_cycle REAL, beta_sigma_scatter REAL,
        period_d REAL, period_sigma_d REAL,
        const_g_s REAL, const_g_err_s REAL, const_r_s REAL,
        const_r_err_s REAL, const_i_s REAL, const_i_err_s REAL,
        era_offset_s REAL, era_offset_err_s REAL,
        gamma_sigma_nightboot REAL, pdot_limit3_nightboot REAL,
        singletons TEXT, note TEXT,
        PRIMARY KEY (variant, model));
    CREATE TABLE IF NOT EXISTS rv_oc_epoch (
        variant TEXT, night TEXT, era_id INTEGER, bands TEXT,
        n_bands INTEGER, n_cycles INTEGER, cycle REAL, oc_s REAL,
        sigma_s REAL, within_rms_s REAL,
        PRIMARY KEY (variant, night));
    CREATE TABLE IF NOT EXISTS rv_oc_chi2 (
        variant TEXT, model TEXT, group_kind TEXT, grp TEXT, n INTEGER,
        chi2 REAL, chi2_per_n REAL, rms_s REAL, sigma_median_s REAL,
        PRIMARY KEY (variant, model, group_kind, grp));
    CREATE TABLE IF NOT EXISTS rv_longitude (
        scope TEXT, grp TEXT, n_epochs INTEGER, n_nights INTEGER,
        mean_deg REAL, se_deg REAL, rms_deg REAL, mean_s REAL, rms_s REAL,
        note TEXT, PRIMARY KEY (scope, grp));
    CREATE TABLE IF NOT EXISTS rv_superhump (
        threshold_mag REAL PRIMARY KEY, n_run_filters INTEGER,
        n_with_contour INTEGER, n_excluding INTEGER, note TEXT);
    CREATE TABLE IF NOT EXISTS rv_colour_curve (
        era_id INTEGER, colour TEXT, method TEXT, bin INTEGER, phase REAL,
        median_mag REAL, err_mag REAL, n INTEGER,
        PRIMARY KEY (era_id, colour, method, bin));
    CREATE TABLE IF NOT EXISTS rv_colour_summary (
        era_id INTEGER, colour TEXT, method TEXT, n_pairs INTEGER,
        n_nights INTEGER, n_bins_used INTEGER, amplitude_mag REAL,
        amplitude_err_mag REAL, phase_red REAL, phase_red_err REAL,
        phase_blue REAL, phase_blue_err REAL, mean_mag REAL,
        tie_sys_clip_mag REAL, tie_sys_raw_mag REAL, median_dt_s REAL,
        n_boot_ok INTEGER, note TEXT,
        PRIMARY KEY (era_id, colour, method));
    CREATE TABLE IF NOT EXISTS rv_colour_repeat (
        colour TEXT, method TEXT, era_a INTEGER, era_b INTEGER,
        n_bins INTEGER, rms_diff_mag REAL, chi2 REAL, dof INTEGER,
        chi2nu REAL, r REAL, scale REAL, amp_diff_mag REAL,
        amp_diff_err_mag REAL, dphase_red REAL, dphase_red_err REAL,
        PRIMARY KEY (colour, method, era_a, era_b));
    CREATE TABLE IF NOT EXISTS rv_edge_quality (
        estimator TEXT, grp TEXT, n INTEGER, dof_min INTEGER,
        dof_max INTEGER, chi2nu_min REAL, chi2nu_q25 REAL,
        chi2nu_median REAL, chi2nu_q75 REAL, chi2nu_max REAL,
        n_above_10 INTEGER, n_below_half INTEGER,
        PRIMARY KEY (estimator, grp));
    CREATE TABLE IF NOT EXISTS rv_edge_example (
        label TEXT, series_key TEXT, cycle INTEGER, band TEXT, night TEXT,
        chi2nu REAL, dof INTEGER, t_edge_bjd REAL, width_s REAL,
        level_bright REAL, step REAL, k INTEGER, dt_s REAL, mag REAL,
        err REAL, PRIMARY KEY (label, k));
    CREATE TABLE IF NOT EXISTS rv_counts (
        scope TEXT, item TEXT, n INTEGER, note TEXT,
        PRIMARY KEY (scope, item));
    CREATE TABLE IF NOT EXISTS rv_figure (
        fig_id TEXT PRIMARY KEY, label TEXT, title TEXT, caption TEXT,
        tables_used TEXT, pdf_path TEXT, png_path TEXT, width_in REAL,
        finding TEXT, built_utc TEXT);
    """)
    con.commit()


def set_meta(con: sqlite3.Connection, items: dict) -> None:
    con.executemany("INSERT OR REPLACE INTO rv_meta (key, value) "
                    "VALUES (?, ?)", [(k, str(v)) for k, v in items.items()])
    con.commit()


def stamp(con: sqlite3.Connection, stage: str) -> None:
    set_meta(con, {f"stage_{stage}": utcnow(),
                   "revision_code_version": REVISION_CODE_VERSION,
                   "git_commit": git_commit()})


def put(con: sqlite3.Connection, stage: str, key: str, value, unit: str,
        fmt: str, note: str, text: str = None,
        origin: str = "measured") -> None:
    """Store ONE scalar the manuscript may quote.

    ``fmt`` tells ``numbers_cv`` how to typeset it: ``int``, ``f0`` ...
    ``f3`` (fixed point), ``sci1``/``sci2`` (scientific), ``p`` (a
    p-value), ``pct`` or ``text``.  ``note`` is the clause a referee needs
    and is mandatory — the macro emitter refuses a number without one.

    ``origin`` is ``measured`` (the default: a result of this programme's
    data), ``constant`` (a choice this stage made in code — a window, a
    grid, a threshold) or ``literature`` (an external physical scale).
    The paper promises a reader can separate measurements from everything
    else with one query, so the distinction is recorded where the number
    is made.
    """
    if not note.strip():
        raise ValueError(f"rv_result {key!r} has no note")
    if origin not in ("measured", "constant", "literature"):
        raise ValueError(f"rv_result {key!r}: unknown origin {origin!r}")
    v = None
    if value is not None:
        v = float(value)
        if not math.isfinite(v):
            v = None
    if v is None and fmt != "text":
        # A quantity that does not exist (the standard error of one pair,
        # a bootstrap that could not run) is not stored at all.  A macro
        # that expanded to a missing-value marker would let a sentence be
        # written around a number nobody measured; an undefined macro
        # stops the build instead.
        return
    con.execute("INSERT OR REPLACE INTO rv_result (key, value, text, unit, "
                "fmt, stage, origin, note) VALUES (?,?,?,?,?,?,?,?)",
                (key, v, text, unit, fmt, stage, origin, note))


def clear_stage(con: sqlite3.Connection, stage: str) -> None:
    """Remove a stage's scalars before it rewrites them, so a key that a
    re-run no longer produces cannot survive as a stale macro."""
    con.execute("DELETE FROM rv_result WHERE stage = ?", (stage,))


def _nz(x):
    """NaN / inf -> None for SQLite."""
    if x is None:
        return None
    try:
        return float(x) if math.isfinite(float(x)) else None
    except (TypeError, ValueError):
        return x


# ===========================================================================
# Shared loaders
# ===========================================================================
def series_context(con: sqlite3.Connection, eph) -> dict:
    """Everything the published edge stage derived per ST LMi series.

    For each solved series: the cloud-vetoed light curve with inflated
    errors, the folded-profile edge phase and ramp width, the within-night
    cadence, and from those the v1 width grid — computed by calling the
    published stage's own functions, in the published stage's own order.
    """
    out = {}
    for sk, tk, era, filt in s9.series_rows(con):
        if tk != TARGET:
            continue
        data = s9.load_series(con, sk)
        if data["t"].size < 30:
            continue
        infl = s9.inflation_for(con, sk)
        t, m, e = data["t"], data["m"], data["e"] * infl
        nights = np.array(data["night"])
        shape = s9.measure_bright_phase_shape(t, m, eph.period_d,
                                              eph.epoch_bjd, err=e)
        if not shape:
            continue
        within = np.concatenate(
            [np.diff(t[nights == n]) for n in np.unique(nights)
             if (nights == n).sum() > 1] or [np.array([])])
        cadence = (float(np.median(within) * 86400.0) if within.size
                   else 219.0)
        out[sk] = {
            "series_key": sk, "era": int(era), "filter": filt,
            "band": rv.BAND_SLOT.get(filt, filt.lower()),
            "t": t, "m": m, "e": e, "nights": nights,
            "frame_id": data["frame_id"],
            "edge_phase": float(shape["edge_phase"]),
            "edge_width_d": float(shape["edge_width_d"]),
            "cadence_s": cadence,
            "v1_width_grid_d": np.array([0.5, 1.0, 2.0, 4.0])
            * float(shape["edge_width_d"]),
        }
    return out


def fwhm_by_frame(con: sqlite3.Connection) -> dict:
    """``frame_id -> seeing FWHM in arcsec`` for the target's frames."""
    return {int(r[0]): float(r[1]) for r in con.execute(
        "SELECT frame_id, fwhm_px * plate_scale FROM cv_frames "
        "WHERE target_key = ? AND fwhm_px IS NOT NULL "
        "AND plate_scale IS NOT NULL", (TARGET,))}


# ===========================================================================
# The three estimators, as one closure over a night's fixed quantities
# ===========================================================================
ESTIMATORS = ("v1", "magc", "v2")


def make_fitter(period_d: float, cadence_s: dict, v1_grid_d: dict):
    """``fit_all(band, t, mag, mag_err, t_guess) -> {estimator: t_edge}``.

    ``cadence_s`` and ``v1_grid_d`` are keyed by band slot.  The callable
    applies the published window and clip and then does what the real
    procedure does, in the real order:

    1.  the published v1 fit (magnitudes, that band's width grid, 361 time
        nodes) decides whether an edge is there at all — its three
        acceptance rules are the published gate;
    2.  if it is, ``magc`` (magnitudes, common width grid) and ``v2``
        (relative flux, common width grid) RE-TIME the same points.  They
        do not re-apply the step-S/N rule — whether the edge exists was
        settled in step 1, and a step measured in flux against bright-phase
        flickering is a different number from the same step in magnitudes —
        but they do keep the two rules that are about the epoch itself: not
        in a gap, not on the grid boundary.

    Every estimator returns NaN when its fit is not a measurement.
    """
    v2_grid_d = np.array(rv.V2_WIDTH_GRID_S) / 86400.0
    nan = float("nan")

    def fit_all(b, t, mag, e, t_guess):
        out = {"v1": nan, "magc": nan, "v2": nan}
        tt, mm, ee = rv.edge_window(t, mag, e, t_guess, period_d)
        if tt.size < p3.EDGE_MIN_POINTS:
            return out
        tg1 = p3.edge_time_grid(t_guess, 3.0 * cadence_s[b] / 86400.0, 361)
        f1 = p3.fit_edge(tt, mm, ee, tg1, v1_grid_d[b], cadence_s[b])
        if not f1.accepted:
            return out
        out["v1"] = f1.t_edge_d
        tg2 = rv.v2_time_grid(t_guess, cadence_s[b], period_d)
        fm = rv.fit_edge_profile(tt, mm, ee, tg2, v2_grid_d, cadence_s[b],
                                 min_snr=0.0)
        if fm["accepted"]:
            out["magc"] = fm["t_edge_d"]
        ff, fe, _ = rv.to_relative_flux(mm, ee)
        fv = rv.fit_edge_profile(tt, ff, fe, tg2, v2_grid_d, cadence_s[b],
                                 min_snr=0.0)
        if fv["accepted"]:
            out["v2"] = fv["t_edge_d"]
        return out

    return fit_all


# ===========================================================================
# STAGE: edges — three estimators and a per-edge bootstrap
# ===========================================================================
def cmd_edges(args) -> None:
    """Re-measure every accepted edge three ways; bootstrap each one."""
    con = connect(args.db)
    ensure_tables(con)
    eph = s9.load_ephemerides(con)[TARGET]
    ctx = series_context(con, eph)
    fwhm = fwhm_by_frame(con)
    print(f"database: {args.db}")
    con.execute("DELETE FROM rv_edge")
    clear_stage(con, "edges")
    v2_grid_d = np.array(rv.V2_WIDTH_GRID_S) / 86400.0
    n_done = 0
    max_refit_diff = 0.0
    for sk, c in sorted(ctx.items()):
        stored = con.execute(
            "SELECT cycle, night, t_edge_bjd, sigma_t_s, width_s, chi2nu, "
            "n_points, depth_mag FROM p3_edge WHERE series_key = ? AND "
            "accepted = 1 ORDER BY cycle", (sk,)).fetchall()
        for row in stored:
            cyc = int(row["cycle"])
            t_guess = eph.epoch_bjd + (cyc + c["edge_phase"]) * eph.period_d
            sel = (np.abs(c["t"] - t_guess)
                   <= p3.EDGE_WINDOW_PHASE * eph.period_d)
            tt, mm, ee = rv.edge_window(c["t"], c["m"], c["e"], t_guess,
                                        eph.period_d)
            n_bad_seeing = int(sum(
                1 for fid in c["frame_id"][sel]
                if fwhm.get(int(fid), 0.0) > FWHM_CUT_ARCSEC))
            cad = c["cadence_s"]
            # --- v1: refit with the published code to prove the window
            # reconstruction, then keep the STORED epoch as the record.
            tg1 = p3.edge_time_grid(t_guess, 3.0 * cad / 86400.0, 361)
            f1 = p3.fit_edge(tt, mm, ee, tg1, c["v1_width_grid_d"], cad)
            diff = (f1.t_edge_d - float(row["t_edge_bjd"])) * 86400.0
            max_refit_diff = max(max_refit_diff, abs(diff))
            # The vectorised surrogate on v1's own grids supplies the
            # valley width and the bootstrap; its chi2 surface is the
            # published fit's (the unit tests hold them equal).
            s1 = rv.fit_edge_profile(tt, mm, ee, tg1, c["v1_width_grid_d"],
                                     cad)
            fits = {"v1": (s1, tt, mm, ee, tg1, c["v1_width_grid_d"])}
            # --- magc: magnitudes, common width grid, 2 s time grid.
            tg2 = rv.v2_time_grid(t_guess, cad, eph.period_d)
            fits["magc"] = (rv.fit_edge_profile(tt, mm, ee, tg2, v2_grid_d,
                                                cad, min_snr=0.0),
                            tt, mm, ee, tg2, v2_grid_d)
            # --- v2: relative flux, common width grid, 2 s time grid.
            ff, fe, _ = rv.to_relative_flux(mm, ee)
            fits["v2"] = (rv.fit_edge_profile(tt, ff, fe, tg2, v2_grid_d,
                                              cad, min_snr=0.0),
                          tt, ff, fe, tg2, v2_grid_d)
            for est, (fit, xt, xy, xe, xg, xw) in fits.items():
                def refit(ystar, _xt=xt, _xe=xe, _xg=xg, _xw=xw, _cad=cad):
                    # No acceptance gate on a replicate: see bootstrap_edge.
                    return rv.fit_edge_profile(_xt, ystar, _xe, _xg, _xw,
                                               _cad)["t_edge_d"]
                boot = rv.bootstrap_edge(
                    fit, refit, n_boot=args.n_boot,
                    seed=rv.SEED + 1000 * (cyc % 100000) + len(est))
                sq = rv.quantisation_sigma_s(fit["grid_step_s"])
                # The epoch of record for v1 is the STORED one; its
                # per-edge error adds the valley indifference, because the
                # published argmin lands wherever rounding puts it inside
                # a flat valley.
                t_edge = (float(row["t_edge_bjd"]) if est == "v1"
                          else fit["t_edge_d"])
                valley_sig = (fit["valley_s"] / math.sqrt(12.0)
                              if est == "v1" and np.isfinite(fit["valley_s"])
                              else 0.0)
                sb = boot["sigma_s"]
                s_edge = (float(math.sqrt(sb ** 2 + sq ** 2
                                          + valley_sig ** 2))
                          if np.isfinite(sb) else None)
                oc_raw = (float((t_edge - (eph.epoch_bjd
                                           + cyc * eph.period_d)) * 86400.0)
                          if np.isfinite(t_edge) else None)
                accepted = (1 if est == "v1" else int(fit["accepted"]))
                con.execute("""
                    INSERT OR REPLACE INTO rv_edge
                    (series_key, cycle, estimator, target_key, era_id,
                     filter, band, night, t_edge_bjd, oc_raw_s, width_s,
                     level_bright, step, snr, chi2, dof, chi2nu, n_points,
                     bracket_s,
                     valley_s, grid_step_s, sigma_quant_s, sigma_boot_s,
                     boot_n_ok, sigma_edge_s, accepted, reason,
                     v1_refit_minus_stored_s, n_fwhm_gt_cut)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
                            ?,?,?,?,?)""", (
                    sk, cyc, est, TARGET, c["era"], c["filter"], c["band"],
                    row["night"], _nz(t_edge), oc_raw,
                    _nz(fit["width_d"] * 86400.0),
                    _nz(fit["level_bright"]), _nz(fit["step"]),
                    _nz(fit["snr"]), _nz(fit["chi2"]), fit["dof"],
                    _nz(fit["chi2nu"]), fit["n_points"],
                    _nz(fit["bracket_s"]), _nz(fit["valley_s"]),
                    _nz(fit["grid_step_s"]), _nz(sq), _nz(sb),
                    boot["n_ok"], s_edge, accepted, fit["reason"],
                    diff if est == "v1" else None, n_bad_seeing))
            n_done += 1
        con.commit()
        print(f"  {sk:14s} {len(stored):3d} edges, three estimators each",
              flush=True)
    # --- scalars -----------------------------------------------------------
    q = con.execute
    put(con, "edges", "rv edges total", n_done, "", "int",
        "accepted ST LMi edges re-measured by the revision; equals the "
        "published accepted count because the revision keeps the published "
        "cycle selection and changes only the estimator")
    put(con, "edges", "rv v1 refit max diff s", max_refit_diff, "s", "f3",
        "largest difference between the published edge epoch and the same "
        "fit re-run by the revision on its reconstruction of the window: "
        "the proof that the revision sees the points the published fit saw")
    for est in ("magc", "v2"):
        n_acc = q("SELECT sum(accepted) FROM rv_edge WHERE estimator=?",
                  (est,)).fetchone()[0]
        put(con, "edges", f"rv {est} edges accepted", n_acc, "", "int",
            f"of the published accepted edges, how many the {est} estimator "
            "also accepts under the same three acceptance rules")
    n_valley = q("SELECT count(*) FROM rv_edge WHERE estimator='v1' AND "
                 "valley_s > 0").fetchone()[0]
    put(con, "edges", "rv v1 edges flat valley", n_valley, "", "int",
        "published edges whose chi-squared is flat (to rounding) across "
        "more than one time-grid node at the best width: the fit says only "
        "'somewhere in this interval' and the published epoch is whichever "
        "node rounding favoured")
    vmax = q("SELECT max(valley_s) FROM rv_edge WHERE estimator='v1'"
             ).fetchone()[0]
    put(con, "edges", "rv v1 valley max s", vmax, "s", "f0",
        "widest flat chi-squared valley among the published edges")
    for est in ("v1", "v2"):
        r = q("SELECT min(sigma_edge_s), max(sigma_edge_s) FROM rv_edge "
              "WHERE estimator=? AND accepted=1 AND sigma_edge_s IS NOT "
              "NULL", (est,)).fetchone()
        med = np.median([x[0] for x in q(
            "SELECT sigma_edge_s FROM rv_edge WHERE estimator=? AND "
            "accepted=1 AND sigma_edge_s IS NOT NULL", (est,))])
        put(con, "edges", f"rv {est} sigma edge median s", med, "s", "f0",
            f"median per-edge bootstrap error of the {est} estimator (wild "
            "residual bootstrap, grid quantisation and valley width in "
            "quadrature); a lower bound, calibrated against the scatter of "
            "real band differences in rv_band_offset")
        put(con, "edges", f"rv {est} sigma edge min s", r[0], "s", "f0",
            f"smallest per-edge bootstrap error, {est} estimator")
        put(con, "edges", f"rv {est} sigma edge max s", r[1], "s", "f0",
            f"largest per-edge bootstrap error, {est} estimator")
    n_dist = q("SELECT count(DISTINCT round(sigma_edge_s, 3)) FROM rv_edge "
               "WHERE estimator='v1' AND sigma_edge_s IS NOT NULL"
               ).fetchone()[0]
    put(con, "edges", "rv v1 sigma edge distinct", n_dist, "", "int",
        "distinct per-edge error values among the published edges under "
        "the bootstrap; the first draft's 'Monte-Carlo sigma' took three "
        "values, one per band")
    for era in (7, 76):
        r = q("SELECT min(grid_step_s), max(grid_step_s) FROM rv_edge WHERE "
              "estimator='v1' AND era_id=?", (era,)).fetchone()
        put(con, "edges", f"rv v1 grid step era {era} max s", r[1], "s",
            "f1", f"largest v1 time-grid step in the {ERA_NAME[era]} era "
            "(six cadences over 360 intervals); its quantisation error "
            "step/sqrt(12) is now in every per-edge budget")
    n_dup = q("""SELECT count(*) FROM (SELECT round(oc_raw_s, 3) AS o,
                 count(*) AS k FROM rv_edge WHERE estimator='v1'
                 GROUP BY o HAVING k > 1)""").fetchone()[0]
    put(con, "edges", "rv v1 oc duplicate values", n_dup, "", "int",
        "O-C values (to 1 ms) shared by more than one published edge: the "
        "signature of the quantised time grid")
    set_meta(con, {"edges_n_boot": args.n_boot,
                   "v2_width_grid_s": ",".join(f"{x:.1f}"
                                               for x in rv.V2_WIDTH_GRID_S),
                   "v2_time_step_s": rv.V2_TIME_STEP_S,
                   "valley_rtol": rv.VALLEY_RTOL})
    stamp(con, "edges")
    con.commit()
    con.close()
    print(f"  {n_done} edges; published epochs reproduced to "
          f"{max_refit_diff:.3g} s")


# ===========================================================================
# STAGE: inject — signed bias, per band, per night, real estimator
# ===========================================================================
def _inject_job(payload: dict) -> dict:
    """One night's injection, in a worker process.  Pure: no database."""
    fit_all = make_fitter(payload["period_d"], payload["cadence_s"],
                          payload["v1_grid_d"])
    res = rv.inject_recover(payload["bands"], payload["period_d"],
                            rv.INJECT_WIDTHS_S, fit_all, ESTIMATORS,
                            n_real=payload["n_real"], seed=payload["seed"],
                            noise=payload["noise"])
    res["night"] = payload["night"]
    res["noise"] = payload["noise"]
    return res


def cmd_inject(args) -> None:
    """The injection grid on every night that contributes an epoch."""
    con = connect(args.db)
    ensure_tables(con)
    eph = s9.load_ephemerides(con)[TARGET]
    ctx = series_context(con, eph)
    print(f"database: {args.db}")
    edges = con.execute(
        "SELECT series_key, night, filter, era_id, cycle, depth_mag "
        "FROM p3_edge WHERE target_key = ? AND accepted = 1", (TARGET,)
    ).fetchall()
    by_night: dict[str, dict] = {}
    for r in edges:
        by_night.setdefault(r["night"], {}).setdefault(
            r["series_key"], []).append(r)
    # Residual pool per (series, night); the series' best-sampled night is
    # the fallback for a night too sparse to yield its own.
    pools: dict[tuple, np.ndarray] = {}
    best_pool: dict[str, tuple] = {}
    for sk, c in ctx.items():
        for night in np.unique(c["nights"]):
            sel = c["nights"] == night
            pool = rv.residual_pool(c["t"][sel], c["m"][sel], eph.period_d,
                                    eph.epoch_bjd)
            pools[(sk, str(night))] = pool
            if pool.size > best_pool.get(sk, ("", np.array([])))[1].size:
                best_pool[sk] = (str(night), pool)
    jobs = []
    meta = {}
    for night, per_series in sorted(by_night.items()):
        bands, cad, grid = {}, {}, {}
        info = {}
        for sk, rows in per_series.items():
            if sk not in ctx:
                continue
            c = ctx[sk]
            if c["band"] not in ("g", "r", "i"):
                # The single z epoch has no partner band and no O-C
                # leverage once bands carry their own constants.
                continue
            sel = c["nights"] == night
            pool = pools.get((sk, night), np.array([]))
            basis = "this night's own residuals"
            if pool.size < 40:
                bn, pool = best_pool.get(sk, ("", np.array([])))
                basis = (f"residuals of {bn}, this series' best-sampled "
                         "night (this night is too sparse for its own)")
            b = c["band"]
            bands[b] = {"t": c["t"][sel], "e": c["e"][sel],
                        "depth_mag": float(np.median(
                            [abs(r["depth_mag"]) for r in rows])),
                        "pool": pool}
            cad[b] = c["cadence_s"]
            grid[b] = c["v1_width_grid_d"]
            info[b] = {"series_key": sk, "era": c["era"], "basis": basis,
                       "pool_n": int(pool.size), "n_real": len(rows),
                       "depth": bands[b]["depth_mag"],
                       "cycles": {int(r["cycle"]) for r in rows}}
        if not bands:
            continue
        meta[night] = info
        noises = ["rolled"] + (["gaussian"] if night in GAUSS_NIGHTS else [])
        for noise in noises:
            jobs.append({"night": night, "noise": noise, "bands": bands,
                         "cadence_s": cad, "v1_grid_d": grid,
                         "period_d": eph.period_d, "n_real": args.n_real,
                         "seed": rv.SEED + int(night.replace("-", ""))})
    print(f"  {len(jobs)} night x noise-model jobs, {args.n_real} "
          f"realisations x {len(rv.INJECT_WIDTHS_S)} widths x 3 estimators "
          f"each, {args.workers} workers", flush=True)
    con.execute("DELETE FROM rv_inject_band")
    con.execute("DELETE FROM rv_inject_pair")
    con.commit()
    with ProcessPoolExecutor(max_workers=args.workers) as pool_exec:
        futs = [pool_exec.submit(_inject_job, j) for j in jobs]
        for k, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            night, noise = res["night"], res["noise"]
            info = meta[night]
            for r in res["per_band"]:
                i = info[r["band"]]
                con.execute("""
                    INSERT OR REPLACE INTO rv_inject_band
                    (night, era_id, series_key, band, estimator, noise,
                     true_width_s, depth_mag, n_try, n_ok, bias_mean_s,
                     bias_se_s, bias_median_s, sigma_s, rms_s, pool_basis,
                     pool_n, n_real_edges)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    night, i["era"], i["series_key"], r["band"],
                    r["estimator"], noise, r["true_width_s"], i["depth"],
                    r["n_try"], r["n_ok"], r["bias_mean_s"], r["bias_se_s"],
                    r["bias_median_s"], r["sigma_s"], r["rms_s"],
                    i["basis"], i["pool_n"], i["n_real"]))
            for r in res["per_pair"]:
                ia, ib = info[r["band_a"]], info[r["band_b"]]
                n_pairs = len(ia["cycles"] & ib["cycles"])
                con.execute("""
                    INSERT OR REPLACE INTO rv_inject_pair
                    (night, era_id, band_a, band_b, estimator, noise,
                     true_width_a_s, true_width_b_s, n_ok, dbias_mean_s,
                     dbias_se_s, dbias_sd_s, n_real_pairs)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    night, ia["era"], r["band_a"], r["band_b"],
                    r["estimator"], noise, r["true_width_a_s"],
                    r["true_width_b_s"], r["n_ok"], r["dbias_mean_s"],
                    r["dbias_se_s"], r["dbias_sd_s"], n_pairs))
            con.commit()
            print(f"    [{k}/{len(jobs)}] {night} ({noise})", flush=True)
    summarise_injection(con)
    set_meta(con, {"inject_n_real": args.n_real,
                   "inject_widths_s": ",".join(f"{x:g}"
                                               for x in rv.INJECT_WIDTHS_S),
                   "inject_guess_scatter_s": rv.GUESS_SCATTER_S,
                   "inject_gauss_nights": ",".join(GAUSS_NIGHTS)})
    stamp(con, "inject")
    con.commit()
    con.close()


def matched_widths(con: sqlite3.Connection) -> dict:
    """Per (era, band): the injected width nearest the MEASURED ramp width.

    "Matched cell" needs a definition that is not chosen after seeing the
    bias.  It is: the injection width closest (in log) to the median ramp
    width the v2 flux fit returned for that band's real edges in that era.
    """
    out = {}
    widths = np.array(rv.INJECT_WIDTHS_S)
    for era, band, w in con.execute(
            "SELECT era_id, band, width_s FROM rv_edge WHERE estimator='v2' "
            "AND accepted=1 AND width_s IS NOT NULL"):
        out.setdefault((int(era), band), []).append(float(w))
    return {k: float(widths[np.argmin(np.abs(np.log(widths)
                                             - math.log(np.median(v))))])
            for k, v in out.items()}


def summarise_injection(con: sqlite3.Connection) -> None:
    """Combine the per-night injection over nights, weighted as the real
    measurement is: a band's bias by its real edges per night, a pair's
    differential bias by its real same-cycle pairs per night."""
    con.execute("DELETE FROM rv_inject_summary")
    clear_stage(con, "inject")
    match = matched_widths(con)

    def combine(rows, wkey, vkey, sekey):
        w = np.array([float(r[wkey]) for r in rows])
        v = np.array([float(r[vkey]) for r in rows])
        se = np.array([float(r[sekey]) for r in rows])
        if w.sum() <= 0:
            return None
        ww = w / w.sum()
        return float(np.sum(ww * v)), float(math.sqrt(np.sum((ww * se) ** 2)))

    # --- per band ---------------------------------------------------------
    for est, noise in con.execute("SELECT DISTINCT estimator, noise FROM "
                                  "rv_inject_band").fetchall():
        for era_sel, era_lab in (((7,), "7"), ((76,), "76"),
                                 ((7, 76), "all")):
            ph = ",".join("?" * len(era_sel))
            for band in ("g", "r", "i"):
                for wd in rv.INJECT_WIDTHS_S:
                    rows = con.execute(f"""
                        SELECT bias_mean_s, bias_se_s, sigma_s,
                               n_real_edges FROM rv_inject_band
                        WHERE estimator=? AND noise=? AND band=?
                          AND true_width_s=? AND era_id IN ({ph})
                          AND bias_mean_s IS NOT NULL""",
                                       (est, noise, band, wd, *era_sel)
                                       ).fetchall()
                    if not rows:
                        continue
                    c = combine(rows, "n_real_edges", "bias_mean_s",
                                "bias_se_s")
                    if c is None:
                        continue
                    wts = np.array([r["n_real_edges"] for r in rows], float)
                    sig = float(np.sum(wts * np.array(
                        [r["sigma_s"] for r in rows])) / wts.sum())
                    is_match = int(era_lab != "all" and match.get(
                        (int(era_lab), band)) == wd)
                    con.execute("""
                        INSERT OR REPLACE INTO rv_inject_summary
                        (kind, era, estimator, noise, band_a, band_b,
                         true_width_a_s, true_width_b_s, n_nights,
                         weight_sum, bias_s, bias_se_s, sigma_s, matched)
                        VALUES ('band',?,?,?,?,'',?,?,?,?,?,?,?,?)""", (
                        era_lab, est, noise, band, wd, wd, len(rows),
                        float(wts.sum()), c[0], c[1], sig, is_match))
            # --- per pair ---------------------------------------------------
            for a, b in rv.BAND_PAIRS:
                for wa in rv.INJECT_WIDTHS_S:
                    for wb in rv.INJECT_WIDTHS_S:
                        rows = con.execute(f"""
                            SELECT dbias_mean_s, dbias_se_s, dbias_sd_s,
                                   n_real_pairs FROM rv_inject_pair
                            WHERE estimator=? AND noise=? AND band_a=?
                              AND band_b=? AND true_width_a_s=?
                              AND true_width_b_s=? AND era_id IN ({ph})
                              AND dbias_mean_s IS NOT NULL
                              AND n_real_pairs > 0""",
                                           (est, noise, a, b, wa, wb,
                                            *era_sel)).fetchall()
                        if not rows:
                            continue
                        c = combine(rows, "n_real_pairs", "dbias_mean_s",
                                    "dbias_se_s")
                        if c is None:
                            continue
                        wts = np.array([r["n_real_pairs"] for r in rows],
                                       float)
                        is_match = int(
                            era_lab != "all"
                            and match.get((int(era_lab), a)) == wa
                            and match.get((int(era_lab), b)) == wb)
                        con.execute("""
                            INSERT OR REPLACE INTO rv_inject_summary
                            (kind, era, estimator, noise, band_a, band_b,
                             true_width_a_s, true_width_b_s, n_nights,
                             weight_sum, bias_s, bias_se_s, sigma_s,
                             matched)
                            VALUES ('pair',?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                            era_lab, est, noise, a, b, wa, wb, len(rows),
                            float(wts.sum()), c[0], c[1],
                            float(np.sum(wts * np.array(
                                [r["dbias_sd_s"] for r in rows]))
                                / wts.sum()), is_match))
    con.commit()
    # --- scalars: the signed matched-cell bias, per band and era ----------
    for (era, band), wd in sorted(match.items()):
        put(con, "inject", f"rv matched width era {era} {band} s", wd, "s",
            "int", f"injected flux-ramp width taken as the matched cell for "
            f"band {band} in the {ERA_NAME.get(era, era)} era: the grid "
            "value nearest the median ramp width the v2 flux fit returned "
            "for that band's real edges")
        for est in ("v1", "v2"):
            r = con.execute("""SELECT bias_s, bias_se_s, sigma_s FROM
                rv_inject_summary WHERE kind='band' AND era=? AND
                estimator=? AND noise='rolled' AND band_a=? AND
                true_width_a_s=?""", (str(era), est, band, wd)).fetchone()
            if r is None:
                continue
            put(con, "inject", f"rv bias {est} era {era} {band} s",
                r["bias_s"], "s", "f0",
                f"SIGNED matched-cell bias (recovered minus injected edge "
                f"time) of the {est} estimator for band {band}, "
                f"{ERA_NAME.get(era, era)} era, real residual noise, "
                "combined over that era's nights weighted by real edges "
                "per night; replaces the first draft's median of absolute "
                "biases")
            put(con, "inject", f"rv bias {est} era {era} {band} err s",
                r["bias_se_s"], "s", "f0",
                "standard error of that signed bias over realisations")
            put(con, "inject", f"rv sigma inj {est} era {era} {band} s",
                r["sigma_s"], "s", "f0",
                f"injection scatter (16-84 half range) of one recovered "
                f"edge, {est} estimator, band {band}, "
                f"{ERA_NAME.get(era, era)} era, matched cell")
    # --- the envelope of the signed bias over every true width ------------
    for est in ("v1", "magc", "v2"):
        r = con.execute("""SELECT min(bias_s), max(bias_s) FROM
            rv_inject_summary WHERE kind='band' AND estimator=? AND
            noise='rolled' AND era IN ('7','76')""", (est,)).fetchone()
        put(con, "inject", f"rv bias {est} min s", r[0], "s", "f0",
            f"most negative signed band bias of the {est} estimator over "
            "every band, era and injected width")
        put(con, "inject", f"rv bias {est} max s", r[1], "s", "f0",
            f"most positive signed band bias of the {est} estimator over "
            "every band, era and injected width")
        for a, b in rv.BAND_PAIRS:
            rr = con.execute("""SELECT min(bias_s), max(bias_s),
                max(bias_se_s) FROM rv_inject_summary WHERE kind='pair' AND
                estimator=? AND noise='rolled' AND era='all' AND band_a=?
                AND band_b=?""", (est, a, b)).fetchone()
            put(con, "inject", f"rv dbias {est} {a}{b} min s", rr[0], "s",
                "f0", f"most negative SIGNED differential bias in "
                f"t({a}) - t({b}) of the {est} estimator under an "
                "achromatic injected edge, over all 25 combinations of true "
                f"ramp width in {a} and in {b}, both eras combined with the "
                "real pair weights")
            put(con, "inject", f"rv dbias {est} {a}{b} max s", rr[1], "s",
                "f0", f"most positive signed differential bias in "
                f"t({a}) - t({b}), same grid")
            put(con, "inject", f"rv dbias {est} {a}{b} se max s", rr[2],
                "s", "f0", "largest standard error of those 25 cells")
            eq = con.execute("""SELECT min(bias_s), max(bias_s) FROM
                rv_inject_summary WHERE kind='pair' AND estimator=? AND
                noise='rolled' AND era='all' AND band_a=? AND band_b=? AND
                true_width_a_s = true_width_b_s""", (est, a, b)).fetchone()
            put(con, "inject", f"rv dbias {est} {a}{b} equal min s", eq[0],
                "s", "f0", f"most negative differential bias "
                f"t({a}) - t({b}) when the two bands share one true ramp "
                "width (the strictly achromatic case)")
            put(con, "inject", f"rv dbias {est} {a}{b} equal max s", eq[1],
                "s", "f0", f"most positive differential bias "
                f"t({a}) - t({b}) when the two bands share one true ramp "
                "width")
    n_n = con.execute("SELECT count(DISTINCT night) FROM rv_inject_band "
                      "WHERE noise='rolled'").fetchone()[0]
    put(con, "inject", "rv inject nights", n_n, "", "int",
        "nights on which the injection test was run: every night that "
        "contributes a g, r or i timing epoch, in both camera eras, so no "
        "epoch's bias is transported from another night")
    n7 = con.execute("SELECT count(DISTINCT night) FROM rv_inject_band "
                     "WHERE noise='rolled' AND era_id=7").fetchone()[0]
    put(con, "inject", "rv inject nights era 7", n7, "", "int",
        "of those, nights in the 2024 High Gain era")
    n_own = con.execute("SELECT count(*) FROM (SELECT DISTINCT night, band "
                        "FROM rv_inject_band WHERE noise='rolled' AND "
                        "pool_basis LIKE 'this night%')").fetchone()[0]
    n_all = con.execute("SELECT count(*) FROM (SELECT DISTINCT night, band "
                        "FROM rv_inject_band WHERE noise='rolled')"
                        ).fetchone()[0]
    put(con, "inject", "rv inject night bands", n_all, "", "int",
        "night-band series injected")
    put(con, "inject", "rv inject night bands own noise", n_own, "", "int",
        "of those, the ones dense enough to supply their own residual "
        "pool; the rest borrow the residuals of the same series' "
        "best-sampled night and keep their own timestamps and errors")
    # The analytic part of the magnitude-space bias, for the report.
    depth = {b: np.median([r[0] for r in con.execute(
        "SELECT depth_mag FROM rv_inject_band WHERE band=? AND era_id=76",
        (b,))]) for b in ("g", "r", "i")}
    for b, d in depth.items():
        put(con, "inject", f"rv depth era 76 {b} mag", d, "mag", "f2",
            f"median fitted edge depth in band {b}, 2025-26 Mode0 era")
        put(con, "inject", f"rv mag midpoint frac {b}",
            rv.magnitude_midpoint_shift(d), "", "f3",
            f"fraction of a linear FLUX ramp at which the MAGNITUDE "
            f"midpoint falls for band {b}'s depth: 0.5 would be unbiased, "
            "and the difference between two bands times the ramp width is "
            "the analytic band bias of a magnitude-space fit")
    con.commit()


# ===========================================================================
# STAGE: offset — the band offset as a paired test (R1, D3)
# ===========================================================================
def _edge_table(con: sqlite3.Connection, estimator: str) -> list[dict]:
    """Accepted edges of one estimator, with both error models attached."""
    rows = con.execute("""
        SELECT e.series_key, e.cycle, e.era_id, e.band, e.night,
               e.t_edge_bjd, e.oc_raw_s, e.sigma_edge_s, e.n_fwhm_gt_cut,
               p.sigma_t_s AS formal, p.sigma_t_mc_s AS mc
        FROM rv_edge e JOIN p3_edge p
          ON p.series_key = e.series_key AND p.cycle = e.cycle
        WHERE e.estimator = ? AND e.accepted = 1
          AND e.band IN ('g','r','i') ORDER BY e.t_edge_bjd""",
                       (estimator,)).fetchall()
    out = []
    for r in rows:
        d = dict(r)
        # The first draft's per-edge error: the larger of the rescaled
        # formal bar and the per-series "Monte-Carlo" constant.
        d["sigma_pub_s"] = max(float(r["formal"] or 0.0),
                               float(r["mc"] or 0.0)) or None
        out.append(d)
    return out


def _pairs(edges: list[dict], a: str, b: str, level: str) -> list[dict]:
    """Paired differences ``t_a - t_b``: per cycle, or per night."""
    out = []
    if level == "cycle":
        by = {}
        for e in edges:
            by.setdefault((e["cycle"]), {})[e["band"]] = e
        for cyc, v in sorted(by.items()):
            if a in v and b in v:
                ea, eb = v[a], v[b]
                sp = (math.hypot(ea["sigma_pub_s"], eb["sigma_pub_s"])
                      if ea["sigma_pub_s"] and eb["sigma_pub_s"] else None)
                sb = (math.hypot(ea["sigma_edge_s"], eb["sigma_edge_s"])
                      if ea["sigma_edge_s"] and eb["sigma_edge_s"] else None)
                out.append({"d": (ea["t_edge_bjd"] - eb["t_edge_bjd"])
                            * 86400.0, "night": ea["night"],
                            "era": ea["era_id"], "sigma_pub": sp,
                            "sigma_boot": sb,
                            "bad_seeing": (ea["n_fwhm_gt_cut"] or 0)
                            + (eb["n_fwhm_gt_cut"] or 0)})
    else:
        by = {}
        for e in edges:
            by.setdefault(e["night"], {}).setdefault(e["band"], []).append(e)
        for night, v in sorted(by.items()):
            if a in v and b in v:
                ma = float(np.mean([e["oc_raw_s"] for e in v[a]]))
                mb = float(np.mean([e["oc_raw_s"] for e in v[b]]))

                def prop(es, key):
                    s = [e[key] for e in es]
                    if any(x is None for x in s):
                        return None
                    return math.sqrt(sum(x * x for x in s)) / len(s)
                pa, pb = prop(v[a], "sigma_pub_s"), prop(v[b], "sigma_pub_s")
                ba, bb = (prop(v[a], "sigma_edge_s"),
                          prop(v[b], "sigma_edge_s"))
                out.append({"d": ma - mb, "night": night,
                            "era": v[a][0]["era_id"],
                            "sigma_pub": (math.hypot(pa, pb)
                                          if pa and pb else None),
                            "sigma_boot": (math.hypot(ba, bb)
                                           if ba and bb else None),
                            "bad_seeing": sum((e["n_fwhm_gt_cut"] or 0)
                                              for e in v[a] + v[b])})
    return out


def cmd_offset(args) -> None:
    """Paired band-offset tests, every estimator, both eras, combined."""
    con = connect(args.db)
    ensure_tables(con)
    eph = s9.load_ephemerides(con)[TARGET]
    print(f"database: {args.db}")
    con.execute("DELETE FROM rv_band_offset")
    con.execute("DELETE FROM rv_halfflux")
    clear_stage(con, "offset")
    store = {}
    for est in ("v1", "magc", "v2"):
        edges = _edge_table(con, est)
        for level in ("cycle", "night"):
            for a, b in rv.BAND_PAIRS:
                allp = _pairs(edges, a, b, level)
                per_era = {}
                for era_lab, sel in (("7", [p for p in allp
                                            if p["era"] == 7]),
                                     ("76", [p for p in allp
                                             if p["era"] == 76]),
                                     ("all", allp)):
                    if not sel:
                        continue
                    d = np.array([p["d"] for p in sel])
                    nights = [p["night"] for p in sel]
                    sp = [p["sigma_pub"] for p in sel]
                    sb = [p["sigma_boot"] for p in sel]
                    t_pub = rv.paired_tests(
                        d, clusters=nights,
                        sigma_budget=(np.array(sp, dtype=float)
                                      if all(x for x in sp) else None),
                        n_trials=rv.N_BAND_PAIR_TRIALS)
                    t_boot = rv.paired_tests(
                        d, sigma_budget=(np.array(sb, dtype=float)
                                         if all(x for x in sb) else None))
                    per_era[era_lab] = t_pub
                    store[(est, level, era_lab, a, b)] = (t_pub, t_boot)
                    con.execute("""
                        INSERT OR REPLACE INTO rv_band_offset
                        (estimator, level, era, band_a, band_b, n,
                         n_nights, mean_s, median_s, sd_s, se_s, t, p_t,
                         n_negative, p_sign, p_wilcoxon, p_perm,
                         perm_basis, p_perm_cluster, p_perm_bonf,
                         p_perm_cluster_bonf, p_wilcoxon_bonf, mean_pub_s,
                         sigma_pub_s, chi2_pub, dof_pub, chi2nu_pub,
                         mean_boot_s, sigma_boot_s, chi2_boot, dof_boot,
                         chi2nu_boot, note)
                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
                                ?,?,?,?,?,?,?,?,?,?,?)""", (
                        est, level, era_lab, a, b, t_pub["n"],
                        t_pub["n_clusters"], _nz(t_pub["mean"]),
                        _nz(t_pub["median"]), _nz(t_pub["sd"]),
                        _nz(t_pub["se"]), _nz(t_pub["t"]), _nz(t_pub["p_t"]),
                        t_pub["n_negative"], _nz(t_pub["p_sign"]),
                        _nz(t_pub["p_wilcoxon"]), _nz(t_pub["p_perm"]),
                        t_pub["perm_basis"], _nz(t_pub["p_perm_cluster"]),
                        _nz(t_pub["p_perm_bonf"]),
                        _nz(t_pub["p_perm_cluster_bonf"]),
                        _nz(t_pub["p_wilcoxon_bonf"]),
                        _nz(t_pub["mean_budget"]),
                        _nz(t_pub["sigma_budget"]), _nz(t_pub["chi2"]),
                        t_pub["dof"], _nz(t_pub["chi2nu"]),
                        _nz(t_boot["mean_budget"]),
                        _nz(t_boot["sigma_budget"]), _nz(t_boot["chi2"]),
                        t_boot["dof"], _nz(t_boot["chi2nu"]),
                        ("t(a) - t(b) in seconds, bluer band first; se is "
                         "sd/sqrt(n); *_pub uses the first draft's per-edge "
                         "errors with NO rescaling; *_boot uses the "
                         "per-edge bootstrap errors")))
                # Inverse-variance combination of the two eras, with the
                # heterogeneity test that licenses (or refuses) it.
                # A scatter-based standard error from two or three pairs
                # is itself uncertain by a factor of two, and weighting by
                # it lets a chance agreement dominate; the combination is
                # made only when each era has MIN_IV_PAIRS, and otherwise
                # the pooled 'all' row is the combined-era estimate.
                if "7" in per_era and "76" in per_era and all(
                        np.isfinite(per_era[k]["se"])
                        and per_era[k]["n"] >= MIN_IV_PAIRS
                        for k in ("7", "76")):
                    comb = rv.combine_inverse_variance(
                        [per_era["7"]["mean"], per_era["76"]["mean"]],
                        [per_era["7"]["se"], per_era["76"]["se"]])
                    con.execute("""
                        INSERT OR REPLACE INTO rv_band_offset
                        (estimator, level, era, band_a, band_b, n, mean_s,
                         se_s, t, chi2_het, dof_het, p_het, note)
                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                        est, level, "iv", a, b,
                        per_era["7"]["n"] + per_era["76"]["n"],
                        comb["mean"], comb["se"], comb["mean"] / comb["se"],
                        comb["chi2_het"], comb["dof_het"], comb["p_het"],
                        "inverse-variance combination of the two eras' "
                        "scatter-based estimates; chi2_het tests that the "
                        "eras agree"))
                    store[(est, level, "iv", a, b)] = (comb, None)
    con.commit()

    # --- the non-parametric half-flux crossing ------------------------------
    ctx = series_context(con, eph)
    phi0 = float(np.angle(np.mean([np.exp(2j * np.pi * c["edge_phase"])
                                   for c in ctx.values()
                                   if c["band"] in "gri"]))
                 / (2 * np.pi)) % 1.0
    half_rows = {}
    for era in (7, 76):
        per_band = {}
        for sk, c in ctx.items():
            if c["era"] != era or c["band"] not in ("g", "r", "i"):
                continue
            nights_ok = {}
            for night in np.unique(c["nights"]):
                sel = c["nights"] == night
                if sel.sum() < 15:
                    continue
                m = c["m"][sel]
                f = 10.0 ** (-0.4 * (m - np.percentile(m, 10.0)))
                ph = p3.phase_of(c["t"][sel], eph.period_d, eph.epoch_bjd)
                dt = (np.mod(ph - phi0 + 0.5, 1.0) - 0.5) \
                    * eph.period_d * 86400.0
                nights_ok[str(night)] = (dt, f)
            per_band[c["band"]] = nights_ok
        # Only nights that every band of this era sampled: the stack then
        # averages the same cycles' wander into every band's profile.
        common = sorted(set.intersection(*[set(v) for v in per_band.values()])
                        ) if per_band else []
        if len(common) < 3:
            continue
        rng = np.random.default_rng(rv.SEED + era)

        def crossing(nights, band):
            dt = np.concatenate([per_band[band][n][0] for n in nights])
            f = np.concatenate([per_band[band][n][1] for n in nights])
            return rv.half_level_crossing(dt, f), dt.size
        base = {b: crossing(common, b) for b in per_band}
        boots = {b: [] for b in per_band}
        diffs = {(a, b): [] for a, b in rv.BAND_PAIRS}
        for _ in range(400):
            pick = [common[i] for i in rng.integers(0, len(common),
                                                    len(common))]
            got = {b: crossing(pick, b)[0]["t_half_s"] for b in per_band}
            for b, v in got.items():
                if v is not None:
                    boots[b].append(v)
            for a, b in rv.BAND_PAIRS:
                if got.get(a) is not None and got.get(b) is not None:
                    diffs[(a, b)].append(got[a] - got[b])

        def half(v):
            if len(v) < 50:
                return None
            lo, hi = np.percentile(v, [16.0, 84.0])
            return float((hi - lo) / 2.0)
        for b, (res, npts) in base.items():
            con.execute("""INSERT OR REPLACE INTO rv_halfflux
                (era, band, n_nights, n_points, t_half_s, t_half_err_s,
                 depth_frac, note) VALUES (?,?,?,?,?,?,?,?)""", (
                str(era), b, len(common), int(npts), res["t_half_s"],
                half(boots[b]), res["depth_frac"],
                f"stacked relative flux of nights sampled in every band, "
                f"against time from ephemeris phase {phi0:.3f}; crossing of "
                "the level half-way between the bright and faint plateaux; "
                "error from a 400-replicate bootstrap over nights"))
        for a, b in rv.BAND_PAIRS:
            ta, tb = base[a][0]["t_half_s"], base[b][0]["t_half_s"]
            if ta is None or tb is None:
                continue
            half_rows[(era, a, b)] = (ta - tb, half(diffs[(a, b)]),
                                      len(common))
    con.commit()

    # --- the cycle-aligned stack: the band offset as a picture -------------
    # The night stack above smears every band by the cycle-to-cycle wander
    # of the edge (~80 s).  Aligning each cycle on ONE reference common to
    # its three bands removes the wander without favouring a band: the
    # reference is the mean of that cycle's g, r and i flux-fit epochs, so
    # each band's profile is displaced from it by exactly that band's own
    # offset.  Each window is normalised to its own fitted levels (bright
    # = 1, faint = 0), so bands of different depth share an axis and the
    # half level is 0.5 for all of them.
    con.execute("DELETE FROM rv_stack")
    v2_grid_d = np.array(rv.V2_WIDTH_GRID_S) / 86400.0
    by_cycle: dict[int, dict] = {}
    for r in con.execute("""SELECT series_key, cycle, band, era_id, night,
            t_edge_bjd FROM rv_edge WHERE estimator='v2' AND accepted=1
            AND band IN ('g','r','i')"""):
        by_cycle.setdefault(int(r["cycle"]), {})[r["band"]] = r
    full = {c: v for c, v in by_cycle.items() if len(v) == 3}
    stack_pts: dict[int, dict] = {}
    for cyc, v in sorted(full.items()):
        ref = float(np.mean([v[b]["t_edge_bjd"] for b in ("g", "r", "i")]))
        entry = {"era": int(v["g"]["era_id"]), "night": v["g"]["night"]}
        for b in ("g", "r", "i"):
            c = ctx[v[b]["series_key"]]
            t_guess = eph.epoch_bjd + (cyc + c["edge_phase"]) * eph.period_d
            tt, mm, ee = rv.edge_window(c["t"], c["m"], c["e"], t_guess,
                                        eph.period_d)
            ff, fe, _ = rv.to_relative_flux(mm, ee)
            fit = rv.fit_edge_profile(
                tt, ff, fe, rv.v2_time_grid(t_guess, c["cadence_s"],
                                            eph.period_d),
                v2_grid_d, c["cadence_s"], min_snr=0.0)
            lb, lf = fit["level_bright"], fit["level_faint"]
            entry[b] = ((tt - ref) * 86400.0, (ff - lf) / (lb - lf))
        stack_pts[cyc] = entry
    for scope, sel in (("aligned 76", [c for c in stack_pts
                                       if stack_pts[c]["era"] == 76]),
                       ("aligned all", list(stack_pts))):
        if len(sel) < 5:
            continue

        def cross(cycles, band):
            dt = np.concatenate([stack_pts[c][band][0] for c in cycles])
            lv = np.concatenate([stack_pts[c][band][1] for c in cycles])
            return rv.half_level_crossing(dt, lv, n_bins=20,
                                          half_span_s=900.0), dt.size
        base = {b: cross(sel, b) for b in ("g", "r", "i")}
        rng = np.random.default_rng(rv.SEED + len(sel))
        boots = {b: [] for b in ("g", "r", "i")}
        diffs = {pr: [] for pr in rv.BAND_PAIRS}
        prof = {b: [] for b in ("g", "r", "i")}
        for _ in range(400):
            pick = [sel[i] for i in rng.integers(0, len(sel), len(sel))]
            got = {}
            for b in ("g", "r", "i"):
                res, _ = cross(pick, b)
                got[b] = res["t_half_s"]
                prof[b].append(res["median"])
                if got[b] is not None:
                    boots[b].append(got[b])
            for a, b in rv.BAND_PAIRS:
                if got[a] is not None and got[b] is not None:
                    diffs[(a, b)].append(got[a] - got[b])

        def half(v):
            if len(v) < 50:
                return None
            lo, hi = np.percentile(v, [16.0, 84.0])
            return float((hi - lo) / 2.0)
        for b, (res, npts) in base.items():
            con.execute("""INSERT OR REPLACE INTO rv_halfflux
                (era, band, n_nights, n_points, t_half_s, t_half_err_s,
                 depth_frac, note) VALUES (?,?,?,?,?,?,?,?)""", (
                scope, b, len(sel), int(npts), res["t_half_s"],
                half(boots[b]), None,
                "level-normalised flux of the cycles timed in all three "
                "bands, each cycle aligned on the mean of its three "
                "flux-fit epochs; n_nights holds the number of CYCLES; "
                "error from a 400-replicate bootstrap over cycles"))
            with np.errstate(invalid="ignore"):
                pe = (np.nanpercentile(np.array(prof[b]), 84.0, axis=0)
                      - np.nanpercentile(np.array(prof[b]), 16.0,
                                         axis=0)) / 2.0
            dtb = np.concatenate([stack_pts[c][b][0] for c in sel])
            edges_ = np.linspace(-900.0, 900.0, 21)
            cnt = np.histogram(dtb, bins=edges_)[0]
            for k in range(res["median"].size):
                con.execute("""INSERT OR REPLACE INTO rv_stack (scope, band,
                    bin, dt_s, level, err, n) VALUES (?,?,?,?,?,?,?)""", (
                    scope, b, k, float(res["centres"][k]),
                    _nz(res["median"][k]), _nz(pe[k]), int(cnt[k])))
        for a, b in rv.BAND_PAIRS:
            ta, tb = base[a][0]["t_half_s"], base[b][0]["t_half_s"]
            if ta is None or tb is None:
                continue
            half_rows[(scope, a, b)] = (ta - tb, half(diffs[(a, b)]),
                                        len(sel))
    con.commit()

    # --- scalars -----------------------------------------------------------
    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "offset", key, value, unit, fmt, note,
            origin=origin)

    for est, est_lab in (("v1", "published magnitude-space estimator"),
                         ("magc", "magnitude fit on the common width grid"),
                         ("v2", "flux-space estimator on the common grid")):
        for a, b in rv.BAND_PAIRS:
            pair = f"{a}{b}"
            for era_lab in ("7", "76", "all"):
                rec = store.get((est, "cycle", era_lab, a, b))
                if rec is None:
                    continue
                t = rec[0]
                tag = f"rv off {est} {pair} era {era_lab}"
                S(f"{tag} n", t["n"], "", "int",
                  f"same-cycle pairs timed in both {a} and {b}, {est_lab}, "
                  f"era {era_lab}")
                S(f"{tag} mean s", t["mean"], "s", "f0",
                  f"mean of t({a}) - t({b}) over those pairs")
                S(f"{tag} se s", t["se"], "s", "f0",
                  "its scatter-based standard error, sd/sqrt(n)")
                S(f"{tag} neg", t["n_negative"], "", "int",
                  f"pairs in which the {a} edge precedes the {b} edge")
                if era_lab == "all":
                    S(f"{tag} t", t["t"], "", "f1",
                      "Student t of the mean against zero")
                    S(f"{tag} nights", t["n_clusters"], "", "int",
                      "distinct nights those pairs fall on")
                    S(f"{tag} p perm", t["p_perm"], "", "p",
                      "two-sided sign-flip (label-permutation) p-value of "
                      f"the mean, {t['perm_basis']}")
                    S(f"{tag} p perm bonf", t["p_perm_bonf"], "", "p",
                      "the same times the three-pair trials factor")
                    S(f"{tag} p cluster", t["p_perm_cluster"], "", "p",
                      "sign-flip p-value flipping whole NIGHTS together: "
                      "the test that does not treat cycles of one night as "
                      "independent")
                    S(f"{tag} p cluster bonf", t["p_perm_cluster_bonf"],
                      "", "p",
                      "night-clustered p-value times the three-pair trials "
                      "factor: the headline significance of this pair")
                    S(f"{tag} p wilcoxon", t["p_wilcoxon"], "", "p",
                      "exact two-sided Wilcoxon signed-rank p-value")
                    S(f"{tag} p sign", t["p_sign"], "", "p",
                      "exact two-sided sign-test p-value")
                    S(f"{tag} chi pub", t["chi2nu"], "", "f2",
                      "reduced chi-squared of the pair differences about "
                      "their weighted mean under the FIRST DRAFT's per-edge "
                      "errors, unclipped; far below one means those errors "
                      "were over-stated")
                    S(f"{tag} dof", t["dof"], "", "int",
                      "its degrees of freedom")
                    S(f"{tag} sigma pub s", t["sigma_budget"], "s", "f0",
                      "propagated error of the weighted mean under the "
                      "first draft's per-edge errors, with no max(chi2nu,1)")
                    tb = rec[1]
                    S(f"{tag} chi boot", tb["chi2nu"], "", "f2",
                      "reduced chi-squared of the same differences under "
                      "the per-edge BOOTSTRAP errors: the calibration of "
                      "the bootstrap (1 = calibrated, above 1 = the "
                      "bootstrap under-states the real scatter)")
            rec = store.get((est, "night", "all", a, b))
            if rec is not None:
                t = rec[0]
                tag = f"rv off {est} {pair} night"
                S(f"{tag} n", t["n"], "", "int",
                  f"nights with an epoch in both {a} and {b}, {est_lab}")
                S(f"{tag} mean s", t["mean"], "s", "f0",
                  f"mean over those nights of (night-mean O-C in {a}) minus "
                  f"(night-mean O-C in {b})")
                S(f"{tag} se s", t["se"], "s", "f0",
                  "its scatter-based standard error over nights")
                S(f"{tag} neg", t["n_negative"], "", "int",
                  f"nights on which {a} precedes {b}")
                S(f"{tag} p perm", t["p_perm"], "", "p",
                  f"sign-flip p-value over nights, {t['perm_basis']}")
                S(f"{tag} p perm bonf", t["p_perm_bonf"], "", "p",
                  "the same times the three-pair trials factor")
            rec = store.get((est, "cycle", "iv", a, b))
            if rec is not None:
                c = rec[0]
                tag = f"rv off {est} {pair} iv"
                S(f"{tag} mean s", c["mean"], "s", "f0",
                  "inverse-variance combination of the two eras' "
                  f"scatter-based {a}-{b} offsets")
                S(f"{tag} se s", c["se"], "s", "f0",
                  "its standard error")
                S(f"{tag} p het", c["p_het"], "", "p",
                  "p-value of the heterogeneity chi-squared between the "
                  "two eras (1 dof): a small value would say the eras "
                  "disagree and must not be combined")
    # --- D3: the offset set beside the estimator's differential bias ------
    for est in ("v1", "magc", "v2"):
        for a, b in rv.BAND_PAIRS:
            rec = store.get((est, "cycle", "all", a, b))
            if rec is None:
                continue
            t = rec[0]
            env = con.execute("""SELECT min(bias_s), max(bias_s),
                max(bias_se_s) FROM rv_inject_summary WHERE kind='pair' AND
                estimator=? AND noise='rolled' AND era='all' AND band_a=?
                AND band_b=?""", (est, a, b)).fetchone()
            if env is None or env[0] is None:
                continue
            # The cell of the envelope that explains MOST of the observed
            # offset is the one closest to it.
            closest = env[0] if abs(t["mean"] - env[0]) < abs(
                t["mean"] - env[1]) else env[1]
            resid = t["mean"] - closest
            err = math.hypot(t["se"], env[2] or 0.0)
            tag = f"rv dthree {est} {a}{b}"
            S(f"{tag} resid s", resid, "s", "f0",
              f"observed same-cycle t({a}) - t({b}) minus the injected "
              "differential bias of the cell that explains the most of it "
              "(the most favourable case for 'estimator bias')")
            S(f"{tag} resid err s", err, "s", "f0",
              "its error: the offset's scatter-based standard error and "
              "the injection cell's, in quadrature")
            S(f"{tag} resid nsigma", resid / err if err > 0 else None, "",
              "f1", "that residual in units of its error: what is left of "
              "the offset after the largest estimator bias the achromatic "
              "injection can produce")
    for (era, a, b), (d, e, nn) in sorted(half_rows.items(),
                                          key=lambda kv: str(kv[0])):
        if isinstance(era, str):
            tag = f"rv halfflux {era} {a}{b}"
            S(f"{tag} s", d, "s", "f0",
              f"difference of the half-level crossing times of the "
              f"cycle-aligned, level-normalised {a} and {b} egress stacks "
              f"({era}): read off the stacked points, not off a ramp")
            S(f"{tag} err s", e, "s", "f0",
              f"its error from a bootstrap over the {nn} cycles timed in "
              "all three bands")
            S(f"{tag} cycles", nn, "", "int",
              "cycles timed in all three bands entering that stack")
            continue
        tag = f"rv halfflux era {era} {a}{b}"
        S(f"{tag} s", d, "s", "f0",
          f"difference of the half-flux crossing times of the stacked "
          f"{a} and {b} egress profiles, {ERA_NAME[era]} era: an estimator "
          "with no ramp model, no width grid and no per-cycle fit")
        S(f"{tag} err s", e, "s", "f0",
          f"its error from a bootstrap over the {nn} nights sampled in "
          "every band")
    # --- the literature's prediction beside the measurement (L4) ----------
    # Standing rule 2: every result carries its predicted scale.  Bailey et
    # al. (1985) measured ST LMi's bright phase longer at longer
    # wavelength; interpolated to an optical pair that predicts the SIGN
    # (bluer band first) and the order of magnitude.  It is set beside the
    # measurement, never used to decide D3: that is the injection's job.
    for a, b in rv.BAND_PAIRS:
        pred = rv.predicted_chromatic_shift_s(eph.period_d, a, b)
        S(f"rv off predicted {a}{b} s", pred, "s", "f0",
          f"PREDICTED t({a}) - t({b}) from the literature, not a "
          "measurement: half the white-light-to-J difference in "
          "bright-phase duration of Bailey et al. (1985, bailey1985 "
          "Sect. 7.3), scaled log-linearly in wavelength to this pair; "
          "order of magnitude only, sign firm (the redder band ends later)",
          origin="literature")
        rec = store.get(("v2", "cycle", "all", a, b))
        if rec is not None and np.isfinite(rec[0]["se"]):
            S(f"rv off v2 {a}{b} over predicted", rec[0]["mean"] / pred, "",
              "f1", f"measured flux-space t({a}) - t({b}) divided by that "
              "literature prediction: positive means the measured sign is "
              "the predicted one")
    S("rv off predicted upper s",
      -0.5 * (rv.BRIGHT_DURATION_K - rv.BRIGHT_DURATION_J)
      * eph.period_d * 86400.0, "s", "f0",
      "upper end of the literature scale for a chromatic edge shift, not "
      "a measurement: half the J-to-K difference in bright-phase duration "
      "(bailey1985 Sect. 7.3: 0.31 and 0.39 cycle), the largest "
      "simultaneous band-to-band shift on record for ST LMi",
      origin="literature")
    # FWHM sensitivity of the headline pair (R8).
    edges = _edge_table(con, "v1")
    allp = _pairs(edges, "g", "i", "cycle")
    clean = [p for p in allp if p["bad_seeing"] == 0]
    if clean:
        t = rv.paired_tests(np.array([p["d"] for p in clean]))
        S("rv off v1 gi fwhmcut n", t["n"], "", "int",
          f"same-cycle g-i pairs left when every edge whose fit window "
          f"holds a frame with seeing FWHM above {FWHM_CUT_ARCSEC:g} arcsec "
          "is removed")
        S("rv off v1 gi fwhmcut mean s", t["mean"], "s", "f0",
          "mean g-i edge offset on that seeing-clean subset")
        S("rv off v1 gi fwhmcut se s", t["se"], "s", "f0",
          "its scatter-based standard error")
    stamp(con, "offset")
    con.commit()
    for est in ("v1", "magc", "v2"):
        for a, b in rv.BAND_PAIRS:
            rec = store.get((est, "cycle", "all", a, b))
            if rec:
                t = rec[0]
                print(f"  {est:5s} {a}-{b}: n={t['n']:2d} mean "
                      f"{t['mean']:+7.1f} +/- {t['se']:5.1f} s  "
                      f"neg {t['n_negative']:2d}  p_perm {t['p_perm']:.4f} "
                      f"p_cluster {t['p_perm_cluster']:.4f}  "
                      f"chi2nu(pub) {t['chi2nu']:.2f}")
    con.close()


# ===========================================================================
# STAGE: oc — the O-C refit (R3)
# ===========================================================================
def _epochs(con: sqlite3.Connection, variant: str, eph) -> list[dict]:
    """Per-night per-band epochs for one variant.

    ``v1``        the published ``p3_oc_night`` rows, published errors.
    ``v1corr``    the same with each epoch's SIGNED matched-cell injection
                  bias (that band, that night, v1 estimator) subtracted.
    ``v2``        epochs rebuilt from the flux-space edges in ``rv_edge``.
    """
    out = []
    if variant in ("v1", "v1corr"):
        match = matched_widths(con)
        for r in con.execute("""SELECT night, filter, era_id, n_cycles,
                cycle_mean, oc_s, oc_sigma_s FROM p3_oc_night
                WHERE target_key = ? ORDER BY night, filter""", (TARGET,)):
            band = rv.BAND_SLOT.get(r["filter"], r["filter"].lower())
            oc = float(r["oc_s"])
            if variant == "v1corr":
                wd = match.get((int(r["era_id"]), band))
                b = con.execute("""SELECT bias_mean_s FROM rv_inject_band
                    WHERE night=? AND band=? AND estimator='v1' AND
                    noise='rolled' AND true_width_s=?""",
                                (r["night"], band, wd)).fetchone()
                if b is None or b[0] is None:
                    if band in ("g", "r", "i"):
                        continue
                else:
                    oc -= float(b[0])
            out.append({"night": r["night"], "band": band,
                        "era": int(r["era_id"]),
                        "n_cycles": int(r["n_cycles"]),
                        "cycle": float(r["cycle_mean"]), "oc_s": oc,
                        "sigma_pub": float(r["oc_sigma_s"])})
    else:
        groups = {}
        for r in con.execute("""SELECT night, band, era_id, cycle, oc_raw_s
                FROM rv_edge WHERE estimator = ? AND accepted = 1""",
                             (variant,)):
            groups.setdefault((r["night"], r["band"], int(r["era_id"])),
                              []).append((r["cycle"], r["oc_raw_s"]))
        for (night, band, era), v in sorted(groups.items()):
            out.append({"night": night, "band": band, "era": era,
                        "n_cycles": len(v),
                        "cycle": float(np.mean([x[0] for x in v])),
                        "oc_s": float(np.mean([x[1] for x in v])),
                        "sigma_pub": None})
    # Remove the global mean: the edge is not at phase zero of the
    # catalogue ephemeris, and every model below fits its own constants.
    mean = float(np.mean([e["oc_s"] for e in out]))
    for e in out:
        e["oc_s"] -= mean
    return out


def _scatter_sigmas(ep: list[dict], band_consts: bool, era_term: bool
                    ) -> np.ndarray:
    """Iterated scatter-based epoch errors: sigma_i = s_band / sqrt(n_i).

    Start from unit single-cycle scatter, fit, read each band's implied
    single-cycle scatter off the residuals, refit, five times.  A band
    with fewer than two epochs takes the pooled value.
    """
    cyc = np.array([e["cycle"] for e in ep])
    oc = np.array([e["oc_s"] for e in ep])
    band = [e["band"] for e in ep]
    era = [e["era"] for e in ep] if era_term else None
    n = np.array([e["n_cycles"] for e in ep], dtype=float)
    sig = 100.0 / np.sqrt(n)
    sb: dict = {}
    for _ in range(5):
        fit = rv.fit_oc_model(cyc, oc, sig, band=band if band_consts
                              else None, era=era)
        if fit["resid"] is None:
            break
        sb = rv.scatter_sigma_by_band(fit["resid"], band, n, fit["dof"])
        pooled = float(math.sqrt(np.mean(n * fit["resid"] ** 2)
                                 * len(ep) / fit["dof"]))
        sig = np.array([sb.get(b, pooled) for b in band]) / np.sqrt(n)
    _scatter_sigmas.last = dict(sb)
    return sig


#: The per-band single-cycle scatter of the most recent call, kept so the
#: O-C stage can publish it without changing the function's return type.
_scatter_sigmas.last = {}


def cmd_oc(args) -> None:
    """The O-C refit, every way the committee asked for it."""
    con = connect(args.db)
    ensure_tables(con)
    eph = s9.load_ephemerides(con)[TARGET]
    period_s = eph.period_d * 86400.0
    print(f"database: {args.db}")
    for tab in ("rv_oc_fit", "rv_oc_epoch", "rv_oc_chi2"):
        con.execute(f"DELETE FROM {tab}")
    clear_stage(con, "oc")
    results = {}

    def run(variant, est, sigma_model, model, ep, band_consts, era_term,
            sig, note, night_boot=False):
        cyc = np.array([e["cycle"] for e in ep])
        oc = np.array([e["oc_s"] for e in ep])
        band = [e["band"] for e in ep] if band_consts else None
        era = [e["era"] for e in ep] if era_term else None
        fit = rv.fit_oc_model(cyc, oc, sig, band=band, era=era)
        if fit["coef"] is None:
            con.execute("""INSERT OR REPLACE INTO rv_oc_fit (variant, model,
                estimator, sigma_model, n_epochs, n_params, note)
                VALUES (?,?,?,?,?,?,?)""", (
                variant, model, est, sigma_model, len(ep), fit["n_params"],
                "RANK DEFICIENT: " + note))
            return None
        names = fit["names"]
        sc = math.sqrt(fit["chi2nu"])

        def coef(name):
            if name not in names:
                return None, None
            i = names.index(name)
            return (float(fit["coef"][i]),
                    float(math.sqrt(fit["cov"][i, i])) * sc)
        g, r_, i_ = coef("band:g"), coef("band:r"), coef("band:i")
        eras_fitted = [nm for nm in names if nm.startswith("era:")]
        eo = coef("era:76") if "era:76" in names else (
            coef(eras_fitted[-1]) if eras_fitted else (None, None))
        pdot = rv.pdot_from_quadratic(fit["gamma"], eph.period_d)
        k = 2.0 / period_s
        sb, ss = fit["gamma_sigma_budget"], fit["gamma_sigma_scatter"]
        gboot = None
        if night_boot:
            # Resample NIGHTS with replacement and refit: the error on the
            # quadratic if the unit of independence is the night.
            nights = sorted({e["night"] for e in ep})
            idx = {nt: [j for j, e in enumerate(ep) if e["night"] == nt]
                   for nt in nights}
            rng = np.random.default_rng(rv.SEED + len(model))
            gam = []
            for _ in range(1000):
                pick = rng.integers(0, len(nights), len(nights))
                ii = np.concatenate([idx[nights[j]] for j in pick])
                f2 = rv.fit_oc_model(
                    cyc[ii], oc[ii], sig[ii],
                    band=[band[j] for j in ii] if band else None,
                    era=[era[j] for j in ii] if era else None)
                if f2["coef"] is not None and np.isfinite(f2["gamma"]):
                    gam.append(f2["gamma"])
            if len(gam) >= 200:
                lo, hi = np.percentile(gam, [16.0, 84.0])
                gboot = float((hi - lo) / 2.0)
        con.execute("""
            INSERT OR REPLACE INTO rv_oc_fit
            (variant, model, estimator, sigma_model, n_epochs,
             n_informative, n_params, chi2, dof, chi2nu, rms_s,
             gamma_s_per_cycle2, gamma_sigma_budget, gamma_sigma_scatter,
             pdot, pdot_sigma_budget, pdot_sigma_scatter,
             pdot_limit3_budget, pdot_limit3_scatter, pdot_nsigma_scatter,
             beta_s_per_cycle, beta_sigma_scatter, period_d, period_sigma_d,
             const_g_s, const_g_err_s, const_r_s, const_r_err_s, const_i_s,
             const_i_err_s, era_offset_s, era_offset_err_s,
             gamma_sigma_nightboot, pdot_limit3_nightboot, singletons, note)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
                    ?,?,?,?,?,?,?,?)""", (
            variant, model, est, sigma_model, fit["n"],
            fit["n_informative"], fit["n_params"], fit["chi2"], fit["dof"],
            fit["chi2nu"], fit["rms_s"], fit["gamma"], sb, ss, pdot,
            sb * k, ss * k, (abs(fit["gamma"]) + 3 * sb) * k,
            (abs(fit["gamma"]) + 3 * ss) * k,
            abs(fit["gamma"]) / ss if ss > 0 else None,
            fit["beta"], fit["beta_sigma_scatter"],
            eph.period_d + fit["beta"] / 86400.0,
            fit["beta_sigma_scatter"] / 86400.0,
            g[0], g[1], r_[0], r_[1], i_[0], i_[1], eo[0], eo[1],
            gboot, ((abs(fit["gamma"]) + 3 * gboot) * k if gboot else None),
            "; ".join(fit["singleton_levels"]), note))
        # Per-band and per-era chi-squared (standing rule 1).
        for kind, labels in (("band", [e["band"] for e in ep]),
                             ("era", [str(e["era"]) for e in ep])):
            for row in rv.per_group_chi2(fit["resid"], sig, labels):
                con.execute("""INSERT OR REPLACE INTO rv_oc_chi2
                    (variant, model, group_kind, grp, n, chi2, chi2_per_n,
                     rms_s, sigma_median_s) VALUES (?,?,?,?,?,?,?,?,?)""", (
                    variant, model, kind, row["group"], row["n"],
                    row["chi2"], row["chi2_per_n"], row["rms_s"],
                    row["sigma_median_s"]))
        results[(variant, model)] = {
            "fit": fit, "pdot": pdot, "lim_b": (abs(fit["gamma"])
                                                + 3 * sb) * k,
            "lim_s": (abs(fit["gamma"]) + 3 * ss) * k, "g": g, "r": r_,
            "i": i_, "era": eo, "gboot": gboot,
            "lim_boot": ((abs(fit["gamma"]) + 3 * gboot) * k
                         if gboot else None)}
        print(f"  {variant:12s} {model:16s} n={fit['n']:2d} "
              f"chi2/dof={fit['chi2']:6.1f}/{fit['dof']:2d}  "
              f"Pdot={pdot:+.2e}  "
              f"lim3(budget)={results[(variant, model)]['lim_b']:.2e} "
              f"lim3(scatter)={results[(variant, model)]['lim_s']:.2e}")
        return fit

    for est, variants in (("v1", ("pub", "scatter")),
                          ("v1corr", ("scatter",)),
                          ("magc", ("scatter",)),
                          ("v2", ("scatter",))):
        try:
            ep_all = _epochs(con, est, eph)
        except sqlite3.OperationalError:
            continue
        if len(ep_all) < 8:
            continue
        for sm in variants:
            variant = f"{est}/{sm}"
            for tag, ep in (("", ep_all),
                            (" no2024", [e for e in ep_all
                                         if e["era"] == 76])):
                if len(ep) < 8:
                    continue
                gri = [e for e in ep if e["band"] in ("g", "r", "i")]

                def sig_for(eps, band_consts, era_term):
                    if sm == "pub":
                        return np.array([e["sigma_pub"] for e in eps])
                    return _scatter_sigmas(eps, band_consts, era_term)
                multi_era = len({e["era"] for e in gri}) > 1
                run(variant, est, sm, "pooled" + tag, ep, False, False,
                    sig_for(ep, False, False),
                    "one constant, period correction, quadratic: the first "
                    "draft's model")
                sig_band = sig_for(gri, True, False)
                if sm == "scatter" and not tag:
                    for b_, v_ in sorted(_scatter_sigmas.last.items()):
                        put(con, "oc", f"rv oc {est} cycle scatter {b_} s",
                            v_, "s", "f0",
                            f"single-cycle edge-time scatter of band {b_} "
                            f"implied by the residuals of the {est} "
                            "per-band-constant O-C fit (epoch variance "
                            "s^2/n): the scatter-based replacement for the "
                            "transported injection budget")
                fb = run(variant, est, sm, "band" + tag, gri, True, False,
                         sig_band,
                         "one constant per band; the single z epoch is "
                         "dropped because its own constant would absorb it",
                         night_boot=True)
                if multi_era:
                    run(variant, est, sm, "band+era" + tag, gri, True, True,
                        sig_for(gri, True, True),
                        "per-band constants plus an era offset: the "
                        "quadratic is then constrained by curvature within "
                        "eras only")
                if fb is None:
                    continue
                # Night-level epochs: band constants removed, one epoch
                # per night.
                consts = {nm.split(":")[1]: float(c)
                          for nm, c in zip(fb["names"], fb["coef"])
                          if nm.startswith("band:")}
                sig_b = sig_for(gri, True, False)
                for e, s in zip(gri, sig_b):
                    e["sigma_s"] = float(s)
                nl = rv.night_level_epochs(gri, consts)
                if not tag:
                    for row in nl:
                        con.execute("""INSERT OR REPLACE INTO rv_oc_epoch
                            (variant, night, era_id, bands, n_bands,
                             n_cycles, cycle, oc_s, sigma_s, within_rms_s)
                            VALUES (?,?,?,?,?,?,?,?,?,?)""", (
                            variant, row["night"], row["era"], row["bands"],
                            row["n_bands"], row["n_cycles"], row["cycle"],
                            row["oc_s"], row["sigma_s"],
                            row["within_rms_s"]))
                nl_ep = [{"night": r["night"], "band": "all",
                          "era": r["era"], "n_cycles": r["n_cycles"],
                          "cycle": r["cycle"], "oc_s": r["oc_s"]}
                         for r in nl]
                nl_sig = np.array([r["sigma_s"] for r in nl])
                run(variant, est, sm, "night" + tag, nl_ep, False, False,
                    nl_sig, "one epoch per night after removing the "
                    "per-band constants; the unit of independence is the "
                    "night")
                if multi_era:
                    run(variant, est, sm, "night+era" + tag, nl_ep, False,
                        True, nl_sig, "night-level epochs plus an era "
                        "offset")
    con.commit()

    # --- scalars -----------------------------------------------------------
    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "oc", key, value, unit, fmt, note,
            origin=origin)

    label = {"pooled": "one constant (the first draft's model)",
             "band": "per-band constants",
             "band+era": "per-band constants and an era offset",
             "night": "night-level epochs",
             "night+era": "night-level epochs and an era offset",
             "pooled no2024": "one constant, 2024 dropped",
             "band no2024": "per-band constants, 2024 dropped",
             "night no2024": "night-level epochs, 2024 dropped"}
    for (variant, model), rec in sorted(results.items()):
        fit = rec["fit"]
        vkey = variant.replace("/", " ")
        tag = f"rv oc {vkey} {model}".replace("+", " plus ")
        desc = (f"O-C fit with {label.get(model, model)}, epochs from the "
                f"{variant.split('/')[0]} estimator, "
                + ("the first draft's budget errors" if variant.endswith(
                    "pub") else "scatter-based errors s_band/sqrt(n)"))
        S(f"{tag} n", fit["n"], "", "int", f"epochs in the {desc}")
        S(f"{tag} chisq", fit["chi2"], "", "f1", f"chi-squared of the {desc}")
        S(f"{tag} dof", fit["dof"], "", "int", "its degrees of freedom")
        S(f"{tag} chinu", fit["chi2nu"], "", "f2",
          "its reduced chi-squared, reported as measured (no clipping at "
          "one)")
        S(f"{tag} pdot", rec["pdot"], "", "sci1",
          f"period derivative from the quadratic term of the {desc}")
        S(f"{tag} pdot nsig", abs(fit["gamma"]) / fit["gamma_sigma_scatter"]
          if fit["gamma_sigma_scatter"] > 0 else None, "", "f1",
          "the quadratic term in units of its scatter-based error")
        S(f"{tag} limit budget", rec["lim_b"], "", "sci1",
          "bound |Pdot| + 3 sigma with the error PROPAGATED from the epoch "
          "errors, unscaled")
        S(f"{tag} limit scatter", rec["lim_s"], "", "sci1",
          "bound |Pdot| + 3 sigma with the error scaled by sqrt(chi2_nu) in "
          "whichever direction the residuals demand")
        if rec["lim_boot"] is not None:
            S(f"{tag} limit nightboot", rec["lim_boot"], "", "sci1",
              "bound |Pdot| + 3 sigma with the error from a 1000-replicate "
              "bootstrap over NIGHTS")
        S(f"{tag} rms s", fit["rms_s"], "s", "f0",
          "unweighted rms of the fit residuals")
        if model == "band":
            for b in ("g", "r", "i"):
                c, e = rec[b]
                S(f"{tag} const {b} s", c, "s", "f0",
                  f"fitted O-C constant of band {b} (about the global "
                  "mean): the per-band edge-time offset the O-C itself "
                  "sees")
                S(f"{tag} const {b} err s", e, "s", "f0",
                  "its scatter-scaled error")
            ci, cg = rec["i"], rec["g"]
            if ci[0] is not None and cg[0] is not None:
                ig = fit["names"].index("band:g")
                ii = fit["names"].index("band:i")
                var = (fit["cov"][ig, ig] + fit["cov"][ii, ii]
                       - 2 * fit["cov"][ig, ii]) * fit["chi2nu"]
                S(f"{tag} g minus i s", cg[0] - ci[0], "s", "f0",
                  "difference of the g and i O-C constants: the band offset "
                  "measured from night epochs rather than same-cycle pairs")
                S(f"{tag} g minus i err s", math.sqrt(max(var, 0.0)), "s",
                  "f0", "its error, with the covariance of the two "
                  "constants")
            S(f"{tag} period d", eph.period_d + fit["beta"] / 86400.0, "d",
              "f8", "period refitted with per-band constants (the linear "
              "term of a fit WITHOUT the quadratic would differ slightly; "
              "this is the quadratic fit's)")
            S(f"{tag} period err d", fit["beta_sigma_scatter"] / 86400.0,
              "d", "sci1", "its scatter-scaled error")
        if model in ("band+era", "night+era"):
            c, e = rec["era"]
            S(f"{tag} era offset s", c, "s", "f0",
              "fitted offset of the 2025-26 era relative to the 2024 era, "
              "the nuisance term for the camera change")
            S(f"{tag} era offset err s", e, "s", "f0",
              "its scatter-scaled error")
    # Chi-squared per band under the first draft's budget (referee minor 4).
    for r in con.execute("""SELECT grp, n, chi2, chi2_per_n, rms_s,
            sigma_median_s FROM rv_oc_chi2 WHERE variant='v1/pub' AND
            model='band' AND group_kind='band'"""):
        S(f"rv oc budget chi per n {r['grp']}", r["chi2_per_n"], "", "f2",
          f"chi-squared per epoch of band {r['grp']}'s residuals about its "
          "own constant under the FIRST DRAFT's budget errors: the "
          "per-band calibration of that budget")
        S(f"rv oc budget n {r['grp']}", r["n"], "", "int",
          f"epochs in band {r['grp']}")
        S(f"rv oc budget rms {r['grp']} s", r["rms_s"], "s", "f0",
          f"rms of band {r['grp']}'s residuals")
        S(f"rv oc budget sigma median {r['grp']} s", r["sigma_median_s"],
          "s", "f0", f"median budget error assigned to band {r['grp']}")
    base = results.get(("v1/pub", "pooled"))
    if base:
        for model in ("band", "band+era", "night", "night+era",
                      "band no2024", "night no2024"):
            rec = results.get(("v1/scatter", model))
            if rec is None:
                continue
            S(f"rv oc limit ratio {model}".replace("+", " plus "),
              rec["lim_s"] / base["lim_b"], "", "f2",
              f"ratio of the |Pdot| + 3 sigma bound under "
              f"{label.get(model, model)} (scatter-based errors) to the "
              "first draft's published bound")
    # --- cycle-count uniqueness under the only PUBLISHED sigma_P (L2) -----
    cc = con.execute("SELECT n_cycles_last, elapsed_d, drift_cycles FROM "
                     "p3_cycle_count WHERE target_key=?", (TARGET,)
                     ).fetchone()
    if cc is not None and cc["elapsed_d"] is not None:
        drift = rv.phase_drift_cycles(cc["elapsed_d"], eph.period_d,
                                      rv.CROPPER_PERIOD_SIGMA_D)
        S("rv cycle sigma p cropper d", rv.CROPPER_PERIOD_SIGMA_D, "d",
          "sci1", "the only published uncertainty on ST LMi's period, from "
          "the literature (cropper1986 Sect. 3: 0.07908908 +/- 0.00000008 "
          "d); the VSX period the pipeline folds on has no published error",
          origin="literature")
        S("rv cycle drift cropper", drift, "cycle", "f3",
          "phase drift accumulated from the catalogue epoch to our last "
          "timed edge under Cropper's (1986) published period error: the "
          "cycle count is unique while this stays below 0.5.  Replaces the "
          "first draft's drift, which assumed an error of half the last "
          "printed digit of a catalogue string")
        S("rv cycle drift cropper margin", 0.5 / drift, "", "f0",
          "half a cycle divided by that drift: the factor by which the "
          "integer cycle count is safe under the published period error")
        S("rv cycle drift assumed", cc["drift_cycles"], "cycle", "f4",
          "the first draft's drift under the ASSUMED half-last-digit "
          "error, for comparison")
        pf = results.get(("v1/scatter", "band"))
        if pf is not None:
            dp = (eph.period_d + pf["fit"]["beta"] / 86400.0) - 0.07908908
            sp = math.hypot(pf["fit"]["beta_sigma_scatter"] / 86400.0,
                            rv.CROPPER_PERIOD_SIGMA_D)
            S("rv period minus cropper sigma", dp / sp, "", "f1",
              "our refitted period (per-band constants) minus Cropper's "
              "(1986) 0.07908908 d, in units of the two errors combined")
    n_nights = con.execute("SELECT count(*) FROM rv_oc_epoch WHERE "
                           "variant='v1/scatter'").fetchone()[0]
    S("rv oc night epochs", n_nights, "", "int",
      "night-level epochs: nights with at least one g, r or i epoch (the "
      "seventeenth night of the first draft carries only a z epoch, which "
      "its own band constant absorbs)")
    stamp(con, "oc")
    con.commit()
    con.close()


# ===========================================================================
# STAGE: longitude — the same residuals as a spot longitude (R5)
# ===========================================================================
def cmd_longitude(args) -> None:
    """Spot-longitude stability in degrees, split by accretion state."""
    con = connect(args.db)
    ensure_tables(con)
    eph = s9.load_ephemerides(con)[TARGET]
    print(f"database: {args.db}")
    con.execute("DELETE FROM rv_longitude")
    clear_stage(con, "longitude")
    ep = [e for e in _epochs(con, "v1", eph) if e["band"] in ("g", "r", "i")]
    sig = _scatter_sigmas(ep, True, False)
    cyc = np.array([e["cycle"] for e in ep])
    oc = np.array([e["oc_s"] for e in ep])
    band = [e["band"] for e in ep]
    # Residuals about per-band constants and a LINEAR ephemeris: a
    # longitude is a departure from uniform rotation, so the quadratic is
    # not removed.
    fit = rv.fit_oc_model(cyc, oc, sig, band=band, quadratic=False)
    res = fit["resid"]
    deg = rv.seconds_to_degrees(res, eph.period_d)
    # State of each epoch: that series' own classification of that night.
    state, bright = [], []
    for e in ep:
        r = con.execute("""SELECT s.state, s.median_mag, t.threshold_mag
            FROM p3_state_night s LEFT JOIN p3_state_series t
              ON t.series_key = s.series_key
            WHERE s.target_key = ? AND s.night = ? AND
                  lower(s.filter) = ? AND s.era_id = ?""",
                        (TARGET, e["night"], e["band"], e["era"])).fetchone()
        state.append(r["state"] if r else "UNCLASSIFIED")
        bright.append((float(r["median_mag"]) - float(r["threshold_mag"]))
                      if r and r["median_mag"] is not None
                      and r["threshold_mag"] is not None else np.nan)
    state = np.array(state)
    bright = np.array(bright)
    nights = np.array([e["night"] for e in ep])

    def row(scope, grp, sel, note):
        if sel.sum() == 0:
            return None
        d = deg[sel]
        se = float(d.std(ddof=1) / math.sqrt(d.size)) if d.size > 1 else None
        con.execute("""INSERT OR REPLACE INTO rv_longitude
            (scope, grp, n_epochs, n_nights, mean_deg, se_deg, rms_deg,
             mean_s, rms_s, note) VALUES (?,?,?,?,?,?,?,?,?,?)""", (
            scope, grp, int(sel.sum()), int(np.unique(nights[sel]).size),
            float(d.mean()), se, float(math.sqrt(np.mean(d ** 2))),
            float(res[sel].mean()),
            float(math.sqrt(np.mean(res[sel] ** 2))), note))
        return float(d.mean()), se, int(sel.sum())

    every = np.ones(len(ep), dtype=bool)
    row("all", "all", every, "every g, r, i epoch; residuals about per-band "
        "constants and a linear ephemeris")
    by_state = {}
    for st in ("HIGH", "INTERMEDIATE", "LOW", "UNCLASSIFIED"):
        by_state[st] = row("state", st, state == st,
                           "epochs whose own series classified the night "
                           f"as {st} in p3_state_night")
    for era in (7, 76):
        row("era", str(era), np.array([e["era"] == era for e in ep]),
            f"epochs of the {ERA_NAME[era]} era")

    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "longitude", key, value, unit, fmt, note,
            origin=origin)

    rms_s = float(math.sqrt(np.mean(res ** 2)))
    S("rv long rms s", rms_s, "s", "f0",
      "rms of the per-night per-band edge epochs about per-band constants "
      "and a linear ephemeris")
    S("rv long rms deg", float(rv.seconds_to_degrees(rms_s, eph.period_d)),
      "deg", "f1", "the same as accretion-spot longitude, 360 dt / P")
    S("rv long rms cycles", rms_s / (eph.period_d * 86400.0), "cycle", "f3",
      "the same in cycles")
    nl = con.execute("SELECT oc_s, cycle FROM rv_oc_epoch WHERE "
                     "variant='v1/scatter'").fetchall()
    if len(nl) >= 5:
        c2 = np.array([r["cycle"] for r in nl])
        o2 = np.array([r["oc_s"] for r in nl])
        f2 = rv.fit_oc_model(c2, o2, np.ones(c2.size), quadratic=False)
        rn = float(math.sqrt(np.sum(f2["resid"] ** 2) / f2["dof"]))
        S("rv long night rms s", rn, "s", "f0",
          "scatter of the NIGHT-LEVEL epochs about a linear ephemeris "
          "(dof-corrected): the night-to-night wander of the edge")
        S("rv long night rms deg",
          float(rv.seconds_to_degrees(rn, eph.period_d)), "deg", "f1",
          "the same as spot longitude: the stability of the accretion "
          "region's longitude from night to night")
        S("rv long span yr", float((c2.max() - c2.min()) * eph.period_d
                                   / 365.25), "yr", "f1",
          "time between the first and last night-level epoch")
    for st, got in by_state.items():
        if got is None:
            continue
        S(f"rv long state {st.lower()} n", got[2], "", "int",
          f"epochs classified {st}")
        S(f"rv long state {st.lower()} mean deg", got[0], "deg", "f1",
          f"mean spot-longitude residual of the {st} epochs")
        S(f"rv long state {st.lower()} err deg", got[1], "deg", "f1",
          "its scatter-based standard error")
    hi, lo = state == "HIGH", state == "LOW"
    if hi.sum() >= 3 and lo.sum() >= 3:
        diff = float(deg[hi].mean() - deg[lo].mean())
        err = float(math.hypot(deg[hi].std(ddof=1) / math.sqrt(hi.sum()),
                               deg[lo].std(ddof=1) / math.sqrt(lo.sum())))
        # Label permutation between the two state groups.
        both = np.concatenate([deg[hi], deg[lo]])
        rng = np.random.default_rng(rv.SEED)
        cnt = 0
        n_hi = int(hi.sum())
        for _ in range(20000):
            p = rng.permutation(both)
            if abs(p[:n_hi].mean() - p[n_hi:].mean()) >= abs(diff) - 1e-12:
                cnt += 1
        S("rv long high minus low deg", diff, "deg", "f1",
          "mean spot longitude of HIGH-state epochs minus LOW-state epochs "
          "(positive = the edge is later in the high state)")
        S("rv long high minus low err deg", err, "deg", "f1",
          "its scatter-based error")
        S("rv long high minus low p", (cnt + 1) / 20001.0, "", "p",
          "label-permutation p-value of that difference (20,000 shuffles "
          "of the state labels)")
        # What difference this sample WOULD have detected (standing rule
        # 2): 3 sigma of the difference.
        S("rv long state detectable deg", 3.0 * err, "deg", "f1",
          "the state-dependent longitude shift these epochs would have "
          "shown at 3 sigma: the effect size behind the null")
    ok = np.isfinite(bright)
    if ok.sum() >= 6:
        slope, _, p = rv.slope_permutation_p(bright[ok], deg[ok])
        S("rv long slope deg per mag", slope, "deg/mag", "f1",
          "slope of spot-longitude residual against the night's median "
          "magnitude relative to that series' own state threshold "
          "(positive = later when fainter)")
        S("rv long slope p", p, "", "p",
          "permutation p-value of that slope")
        S("rv long slope n", int(ok.sum()), "", "int",
          "epochs with a state measurement entering the slope")
    # --- the predicted scale of a state-dependent shift (L5) ---------------
    hist_s = rv.historical_edge_shift_s(eph.period_d)
    hist_deg = float(rv.seconds_to_degrees(hist_s, eph.period_d))
    S("rv long historical shift s", hist_s, "s", "f0",
      "PREDICTED SCALE from the literature, not a measurement: how far one "
      "bright-phase edge moved between 1982 and 1985 (cropper1986 Sect. 4: "
      "duration 0.33 to 0.38 cycle; half the change, if symmetric)",
      origin="literature")
    S("rv long historical shift deg", hist_deg, "deg", "f0",
      "the same literature scale as spot longitude: the size of a real "
      "accretion-geometry change in this star's record", origin="literature")
    if hi.sum() >= 3 and lo.sum() >= 3:
        S("rv long historical over detectable", hist_deg / (3.0 * err), "",
          "f1", "the historical 1982-85 edge shift divided by the "
          "state-dependent shift these epochs would show at 3 sigma: above "
          "one means a change of that historical size between high and low "
          "state would have been detected")
        S("rv long historical excluded nsigma",
          (hist_deg - abs(diff)) / err, "", "f1",
          "how many standard errors the measured high-minus-low longitude "
          "difference lies below the historical 1982-85 shift")
    S("rv long rms over historical", float(rv.seconds_to_degrees(
        rms_s, eph.period_d)) / hist_deg, "", "f2",
      "rms spot-longitude scatter of our epochs divided by the historical "
      "1982-85 edge shift")
    # --- the physical scales of a period derivative (standing rule 2) -----
    lim = con.execute("SELECT pdot_limit3_scatter FROM rv_oc_fit WHERE "
                      "variant='v1/scatter' AND model='band'").fetchone()
    lim_era = con.execute("SELECT pdot_limit3_scatter FROM rv_oc_fit WHERE "
                          "variant='v1/scatter' AND model='band+era'"
                          ).fetchone()
    gr = rv.gr_orbital_pdot(eph.period_d)
    S("rv pdot gr scale", gr, "", "sci1",
      f"EXTERNAL SCALE from the literature's quadrupole formula, not a "
      f"measurement: |dP/dt| from gravitational radiation alone at this "
      f"period for point masses of {rv.GR_M1_MSUN:g} and "
      f"{rv.GR_M2_MSUN:g} solar masses", origin="literature")
    for tau in rv.SECULAR_TAU_GYR:
        S(f"rv pdot secular {int(tau)} gyr", rv.secular_pdot(eph.period_d,
                                                             tau), "",
          "sci1", f"EXTERNAL SCALE from the literature (knigge2011): P "
          f"divided by an evolution time of {tau:g} Gyr, the range of "
          "secular orbital evolution below the period gap",
          origin="literature")
    S("rv pdot async scale", rv.ASYNC_POLAR_PDOT, "", "sci2",
      "EXTERNAL SCALE from the literature: spin-period derivative of the "
      "asynchronous polar V1500 Cyg relaxing toward synchronism "
      "(schmidt1995, as quoted by pavlenko2018)", origin="literature")
    lib = rv.libration_pdot(eph.period_d)
    S("rv pdot libration scale", lib, "", "sci1",
      f"EXTERNAL SCALE from the literature: apparent |dP/dt| at the peak "
      f"curvature of a DP Leo-like spot libration "
      f"({rv.LIBRATION_AMP_DEG:g} deg amplitude, "
      f"{rv.LIBRATION_PERIOD_YR:g} yr period; beuermann2014) at ST LMi's "
      "period", origin="literature")
    for tag, row_ in (("", lim), (" era", lim_era)):
        if row_ is None or not row_[0]:
            continue
        which = ("per-band-constant" if not tag else
                 "per-band-constant, era-offset")
        S(f"rv pdot limit{tag} over secular",
          row_[0] / rv.secular_pdot(eph.period_d, rv.SECULAR_TAU_GYR[0]),
          "", "sci1", f"the {which} |Pdot| bound divided by the fastest "
          "secular orbital rate: how far above orbital evolution the bound "
          "sits, i.e. why it is not an orbital result")
        S(f"rv pdot async over limit{tag}", rv.ASYNC_POLAR_PDOT / row_[0],
          "", "f1", f"the V1500 Cyg resynchronisation rate divided by the "
          f"{which} bound: above one means spin evolution at that rate is "
          "excluded")
        S(f"rv pdot limit{tag} over libration", row_[0] / lib, "", "f0",
          f"the {which} bound divided by the DP Leo-like libration scale: "
          "above one means such a libration is NOT constrained")
    for b, lam in rv.BAND_WAVELENGTH_A.items():
        S(f"rv cyclotron harmonic {b}",
          rv.cyclotron_harmonic(lam, rv.STLMI_FIELD_MG), "", "f0",
          f"EXTERNAL SCALE from the literature: cyclotron harmonic number "
          f"at band {b} ({lam:.0f} A) for the {rv.STLMI_FIELD_MG:g} MG "
          "field of ST LMi's main accretion region (campbell2008c)",
          origin="literature")
    stamp(con, "longitude")
    con.commit()
    con.close()
    print(f"  spot longitude rms {rms_s:.0f} s = "
          f"{float(rv.seconds_to_degrees(rms_s, eph.period_d)):.1f} deg over "
          f"{len(ep)} epochs")


# ===========================================================================
# STAGE: superhump — what the YZ Cnc contours exclude (R4)
# ===========================================================================
def cmd_superhump(args) -> None:
    """Excluded-amplitude statement from the blind-search contours."""
    con = connect(args.db)
    ensure_tables(con)
    print(f"database: {args.db}")
    con.execute("DELETE FROM rv_superhump")
    clear_stage(con, "superhump")
    rows = con.execute("SELECT series_key, night, amp90_blind, "
                       "superhump_floor, episode FROM p4_outburst").fetchall()
    summ = rv.excluded_amplitude_summary([r["amp90_blind"] for r in rows],
                                         SUPERHUMP_THRESHOLDS)
    for th, n in sorted(summ["n_excluding"].items()):
        con.execute("INSERT OR REPLACE INTO rv_superhump VALUES (?,?,?,?,?)",
                    (th, summ["n_run_filters"], summ["n_with_contour"], n,
                     "run-filters whose 90% blind-recovery contour is at or "
                     "below this semi-amplitude: the runs in which a "
                     "periodic signal of this size is excluded"))

    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "superhump", key, value, unit, fmt, note,
            origin=origin)

    S("rv sh run filters", summ["n_run_filters"], "", "int",
      "YZ Cnc dense run-filters (night x filter) examined for a superhump")
    S("rv sh with contour", summ["n_with_contour"], "", "int",
      "of those, run-filters for which a 90% blind-recovery contour was "
      "measured; the others carry no sensitivity statement at all")
    S("rv sh contour min mmag", 1000 * summ["contour_min"], "mmag", "f0",
      "smallest 90% blind-recovery semi-amplitude among the run-filters "
      "with a contour: no superhump below this is excluded anywhere")
    S("rv sh contour max mmag", 1000 * summ["contour_max"], "mmag", "f0",
      "largest 90% blind-recovery semi-amplitude: only above this is a "
      "superhump excluded in every run-filter that has a contour")
    S("rv sh contour min full mag", 2 * summ["contour_min"], "mag", "f2",
      "the smallest contour as a FULL (peak-to-peak) amplitude")
    S("rv sh contour max full mag", 2 * summ["contour_max"], "mag", "f2",
      "the largest contour as a full amplitude")
    for th, n in sorted(summ["n_excluding"].items()):
        S(f"rv sh excluding {int(round(1000 * th))} mmag", n, "", "int",
          f"run-filters in which a periodic signal of semi-amplitude "
          f"{1000 * th:.0f} mmag is excluded at 90% (contour at or below "
          "it)")
    S("rv sh peak semi mmag", 1000 * rv.SUPERHUMP_SEMI_PEAK_MAG, "mmag",
      "f0", "from the literature, not a measurement: maximum superhump "
      "semi-amplitude of SU UMa stars at low inclination, half of 0.25 mag "
      "full (smak2010 abstract; kato2012 Sect. 4.7)", origin="literature")
    S("rv sh tess semi mmag", 1000 * rv.SUPERHUMP_SEMI_TESS_MAG, "mmag",
      "f0", "from the literature, not a measurement: YZ Cnc's own superhump "
      "semi-amplitude at the precursor and early plateau in TESS, half of "
      "0.3 mag full (dai2026 Sect. 4.2)", origin="literature")
    n_peak = summ["n_excluding"][rv.SUPERHUMP_SEMI_PEAK_MAG]
    S("rv sh excluding peak", n_peak, "", "int",
      "run-filters whose 90% blind-recovery contour lies at or below the "
      "published peak superhump semi-amplitude (125 mmag): the run-filters "
      "in which a FULLY DEVELOPED superhump is excluded")
    S("rv sh not excluding peak", summ["n_run_filters"] - n_peak, "", "int",
      "run-filters in which even a fully developed superhump is not "
      "excluded: those with no contour and those whose contour is above "
      "125 mmag")
    # --- where each run sits in the supercycle (cv-literature item 3) ----
    # A superhump can outlive its superoutburst (in YZ Cnc late superhumps
    # were seen into the next normal outburst), so "no superoutburst on
    # the night" is weaker than "no superoutburst recently".  The episodes
    # are CV-S7's, from AAVSO and ZTF; nothing is re-classified here.
    con.execute("""CREATE TABLE IF NOT EXISTS rv_supercycle (
        night TEXT PRIMARY KEY, prev_super_end TEXT, days_since REAL,
        next_super_start TEXT, days_until REAL, note TEXT)""")
    con.execute("DELETE FROM rv_supercycle")
    supers = con.execute("""SELECT start_night, end_night FROM
        cv_ext_episode WHERE target='yzcnc' AND kind='SUPEROUTBURST'
        ORDER BY start_night""").fetchall()

    def jd(day):
        return datetime.strptime(day, "%Y-%m-%d").toordinal()
    since, until = [], []
    for night in sorted({r["night"] for r in rows}):
        prev = [x for x in supers if x["end_night"] <= night]
        nxt = [x for x in supers if x["start_night"] > night]
        ds = (jd(night) - jd(prev[-1]["end_night"])) if prev else None
        du = (jd(nxt[0]["start_night"]) - jd(night)) if nxt else None
        con.execute("INSERT OR REPLACE INTO rv_supercycle VALUES "
                    "(?,?,?,?,?,?)", (
                        night, prev[-1]["end_night"] if prev else None, ds,
                        nxt[0]["start_night"] if nxt else None, du,
                        "local night of the dense run against the "
                        "SUPEROUTBURST episodes of cv_ext_episode (AAVSO "
                        "and ZTF, classified by CV-S7)"))
        if ds is not None:
            since.append(ds)
        if du is not None:
            until.append(du)
    if since:
        S("rv sh days since super min", min(since), "d", "int",
          "shortest interval from the end of the preceding superoutburst "
          "(as CV-S7 classified the AAVSO record) to a dense YZ Cnc run")
        S("rv sh days since super max", max(since), "d", "int",
          "longest such interval")
    if until:
        S("rv sh days until super min", min(until), "d", "int",
          "shortest interval from a dense YZ Cnc run to the start of the "
          "next superoutburst")
    nights = sorted({r["night"] for r in rows})
    S("rv sh runs", len(nights), "", "int",
      "distinct YZ Cnc dense runs (nights) behind those run-filters")
    n_with_nights = len({r["night"] for r in rows
                         if r["amp90_blind"] is not None})
    S("rv sh runs with contour", n_with_nights, "", "int",
      "runs with a contour in at least one filter")
    stamp(con, "superhump")
    con.commit()
    con.close()
    print(f"  {summ['n_with_contour']} of {summ['n_run_filters']} "
          f"run-filters carry a contour, "
          f"{1000 * summ['contour_min']:.0f}-"
          f"{1000 * summ['contour_max']:.0f} mmag; excluded at the "
          f"{1000 * rv.SUPERHUMP_SEMI_PEAK_MAG:.0f} mmag published peak in "
          f"{summ['n_excluding'][rv.SUPERHUMP_SEMI_PEAK_MAG]}")


# ===========================================================================
# STAGE: yzrefold — YZ Cnc on the period that has an error bar (R4, L10/L11)
# ===========================================================================
def cmd_yzrefold(args) -> None:
    """Refold every YZ Cnc scope on van Paradijs et al.'s (1994) period.

    The first draft folded YZ Cnc on the VSX period 0.0868 d, whose
    published error (Shafter & Hessman 1988: 0.0002 d) accumulates 0.027
    cycle of phase per day — not the "< 0.01 cycle" the draft quoted, which
    came from an assumed half-last-digit error.  A period thirty times
    more precise exists (0.086924 +/- 0.000007 d).  This stage refits the
    same points, with the same errors and the same harmonic model
    (``final_science.fold_fit``), on both periods, so that what changes is
    the period and nothing else, and restates the between-night drift
    under each published error.

    It writes ``rv_yz_refold`` and touches no ``p4_`` table.  The first
    row of evidence is the reproduction: the fold on the VSX period must
    return the amplitude ``p4_run`` stores, or the two stages are not
    fitting the same points.
    """
    import run_cv_final as s10
    from macro_phot import final_science as fs
    con = connect(args.db)
    ensure_tables(con)
    print(f"database: {args.db}")
    con.execute("""CREATE TABLE IF NOT EXISTS rv_yz_refold (
        scope TEXT PRIMARY KEY, series_key TEXT, kind TEXT, state TEXT,
        nights TEXT, n_points INTEGER, span_d REAL,
        amp_stored_mag REAL, amp_vsx_mag REAL, amp_vsx_sigma REAL,
        phase_vsx REAL, chi2nu_vsx REAL, amp_vp_mag REAL,
        amp_vp_sigma REAL, phase_vp REAL, chi2nu_vp REAL,
        drift_sh_cycles REAL, drift_vp_cycles REAL, slip_cycles REAL,
        amp90_self_mag REAL, amp90_field_mag REAL, detection TEXT)""")
    con.execute("DELETE FROM rv_yz_refold")
    clear_stage(con, "yzrefold")
    p_vsx = float(con.execute("SELECT period_d FROM p3_ephemeris WHERE "
                              "target_key='yzcnc'").fetchone()[0])
    p_vp = rv.YZCNC_VP_PERIOD_D
    infl = {r["series_key"]: (r["chi2_inflation"] or 1.0)
            for r in con.execute("SELECT series_key, chi2_inflation FROM "
                                 "cv_series WHERE target_key='yzcnc'")}
    worst = 0.0
    out = []
    for r in con.execute("SELECT * FROM p4_run WHERE target_key='yzcnc' "
                         "ORDER BY scope").fetchall():
        sk = r["series_key"]
        nights = str(r["nights"]).split("+")
        parts = [s10.load_run(con, sk, n) for n in nights]
        parts = [(n, d) for n, d in zip(nights, parts)
                 if d["t"].size >= s10.FULL_ORBIT_MIN_POINTS]
        if not parts:
            continue
        t = np.concatenate([d["t"] for _n, d in parts])
        m = np.concatenate([d["m"] for _n, d in parts])
        e = np.concatenate([d["e"] for _n, d in parts]) \
            * float(infl.get(sk, 1.0) or 1.0)
        lab = np.concatenate([np.full(d["t"].size, n) for n, d in parts])
        idx = lab if r["kind"] == "block" else None
        f_vsx = fs.fold_fit(t, m, e, p_vsx, s10.YZ_EPOCH_BJD,
                            night_index=idx)
        f_vp = fs.fold_fit(t, m, e, p_vp, s10.YZ_EPOCH_BJD, night_index=idx)
        if r["hump_amp"] is not None:
            worst = max(worst, abs(f_vsx["amp"] - r["hump_amp"]))
        span = float(t.max() - t.min())
        row = {
            "scope": r["scope"], "kind": r["kind"], "state": r["state"],
            "amp_vsx": f_vsx["amp"], "amp_vp": f_vp["amp"],
            "sig_vp": f_vp["amp_sigma"],
            "dphase": rv.circ_diff(f_vp["phase_max"], f_vsx["phase_max"]),
            "drift_sh": rv.phase_drift_cycles(span, p_vsx,
                                              rv.YZCNC_SH_PERIOD_SIGMA_D),
            "drift_vp": rv.phase_drift_cycles(span, p_vp,
                                              rv.YZCNC_VP_PERIOD_SIGMA_D),
            "slip": abs(span / p_vsx - span / p_vp),
            "self": r["amp90_self"], "field": r["amp90_field"],
            "detection": r["detection"]}
        out.append(row)
        con.execute("""INSERT OR REPLACE INTO rv_yz_refold VALUES
            (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
            r["scope"], sk, r["kind"], r["state"], r["nights"], int(t.size),
            span, r["hump_amp"], f_vsx["amp"], f_vsx["amp_sigma"],
            f_vsx["phase_max"], f_vsx["chi2nu"], f_vp["amp"],
            f_vp["amp_sigma"], f_vp["phase_max"], f_vp["chi2nu"],
            row["drift_sh"], row["drift_vp"], row["slip"],
            r["amp90_self"], r["amp90_field"], r["detection"]))

    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "yzrefold", key, value, unit, fmt, note, origin=origin)

    S("rv yz scopes", len(out), "", "int",
      "YZ Cnc scopes (single runs and multi-night blocks, per filter) "
      "refolded")
    S("rv yz refold reproduction mmag", 1000 * worst, "mmag", "f3",
      "largest difference between the hump semi-amplitude p4_run stores "
      "and the same fold repeated here on the same period: the proof the "
      "refold fits the published points")
    S("rv yz vp period d", p_vp, "d", "f6",
      "from the literature, not a measurement: YZ Cnc's orbital period "
      "(vanparadijs1994 Sect. 6), the published value with the smallest "
      "error", origin="literature")
    S("rv yz vp period sigma d", rv.YZCNC_VP_PERIOD_SIGMA_D, "d", "sci1",
      "from the literature: its published error (vanparadijs1994)",
      origin="literature")
    S("rv yz sh period sigma d", rv.YZCNC_SH_PERIOD_SIGMA_D, "d", "sci1",
      "from the literature: the published error of the period VSX "
      "carries (shafter1988, as quoted by later papers)",
      origin="literature")
    S("rv yz day drift sh", rv.phase_drift_cycles(
        1.0, p_vsx, rv.YZCNC_SH_PERIOD_SIGMA_D), "cycle", "f3",
      "phase drift accumulated over one day under the PUBLISHED error of "
      "the period the first draft folded on: the between-night drift the "
      "draft quoted as below 0.01 cycle")
    S("rv yz day drift vp", rv.phase_drift_cycles(
        1.0, p_vp, rv.YZCNC_VP_PERIOD_SIGMA_D), "cycle", "f4",
      "the same over one day under van Paradijs et al.'s period error")
    S("rv yz day slip", abs(1.0 / p_vsx - 1.0 / p_vp), "cycle", "f3",
      "phase slip per day between a fold on 0.0868 d and one on "
      "0.086924 d")
    blocks = [x for x in out if x["kind"] == "block"]
    runs = [x for x in out if x["kind"] != "block"]
    if blocks:
        S("rv yz block drift sh max", max(x["drift_sh"] for x in blocks),
          "cycle", "f3", "largest phase drift across a multi-night block "
          "under the published error of the VSX period")
        S("rv yz block drift vp max", max(x["drift_vp"] for x in blocks),
          "cycle", "f4", "the same under van Paradijs et al.'s error: on "
          "that period the two nights of a block ARE on one phase axis")
        S("rv yz block amp change max mmag", 1000 * max(
            abs(x["amp_vp"] - x["amp_vsx"]) for x in blocks), "mmag", "f1",
          "largest change in a block's fitted hump semi-amplitude when "
          "the fold period is changed to van Paradijs et al.'s")
        S("rv yz block amp vp min mmag",
          1000 * min(x["amp_vp"] for x in blocks), "mmag", "f0",
          "smallest block hump semi-amplitude on van Paradijs et al.'s "
          "period")
        S("rv yz block amp vp max mmag",
          1000 * max(x["amp_vp"] for x in blocks), "mmag", "f0",
          "largest block hump semi-amplitude on van Paradijs et al.'s "
          "period")
    if runs:
        S("rv yz run amp change max mmag", 1000 * max(
            abs(x["amp_vp"] - x["amp_vsx"]) for x in runs), "mmag", "f1",
          "largest change in a single run's fitted hump semi-amplitude "
          "when the fold period is changed: within one night the two "
          "periods differ by a small fraction of a cycle")
    # Night-to-night coherence on the precise period: with 0.001 cycle of
    # drift between the two May nights, a hump locked to the orbit must
    # come back at the same phase.
    may = {}
    for x, r in zip(out, con.execute(
            "SELECT nights, phase_vp, phase_vsx FROM rv_yz_refold ORDER BY "
            "scope").fetchall()):
        if x["kind"] != "block" and str(r["nights"]).startswith("2024-05"):
            may.setdefault(r["nights"], []).append((r["phase_vp"],
                                                    r["phase_vsx"]))
    if len(may) == 2:
        def cmean(v):
            z = np.exp(2j * np.pi * np.asarray(v)).mean()
            return float(np.mod(np.angle(z) / (2 * np.pi), 1.0))
        n1, n2 = sorted(may)
        for k, tag, lab in ((0, "vp", "van Paradijs et al.'s period"),
                            (1, "vsx", "the VSX period")):
            sh = abs(rv.circ_diff(cmean([v[k] for v in may[n2]]),
                                  cmean([v[k] for v in may[n1]])))
            S(f"rv yz night shift {tag}", sh, "cycle", "f2",
              f"shift of the fitted hump phase (circular mean over the "
              f"three filters) from {n1} to {n2}, folded on {lab}")
    # --- against the prior detection (L11) --------------------------------
    prior = rv.YZCNC_PRIOR_HUMP_SEMI_MAG
    S("rv yz prior hump semi mmag", 1000 * prior, "mmag", "f0",
      "from the literature, not a measurement: semi-amplitude of the "
      "orbital hump van Paradijs et al. (1994, Sect. 6) detected in "
      "quiescence, half of 0.5 mag full", origin="literature")
    quiet = [x for x in out if x["state"] == "QUIESCENT"]
    with_self = [x for x in quiet if x["self"] is not None]
    S("rv yz quiescent scopes", len(quiet), "", "int",
      "quiescent YZ Cnc scopes")
    S("rv yz quiescent with contour", len(with_self), "", "int",
      "of those, scopes with a red-noise (self-noise) 90% recovery contour")
    S("rv yz quiescent excluding prior", sum(
        1 for x in with_self if x["self"] <= prior + 1e-12), "", "int",
      "quiescent scopes whose red-noise 90% contour lies at or below the "
      "250 mmag hump of van Paradijs et al.: the scopes in which a hump "
      "of that published size is excluded")
    if quiet:
        S("rv yz quiescent amp vp min mmag",
          1000 * min(x["amp_vp"] for x in quiet), "mmag", "f0",
          "smallest fitted quiescent hump semi-amplitude, van Paradijs "
          "period")
        S("rv yz quiescent amp vp max mmag",
          1000 * max(x["amp_vp"] for x in quiet), "mmag", "f0",
          "largest fitted quiescent hump semi-amplitude, van Paradijs "
          "period")
        S("rv yz prior over ours min", prior / max(x["amp_vp"]
                                                   for x in quiet), "",
          "f1", "the published 250 mmag hump divided by our largest "
          "fitted quiescent semi-amplitude")
    n_det = sum(1 for x in quiet if x["detection"] == "DETECTED")
    S("rv yz quiescent detected", n_det, "", "int",
      "quiescent scopes in which the first draft's own two-test rule "
      "graded the orbital hump DETECTED")
    stamp(con, "yzrefold")
    con.commit()
    con.close()
    print(f"  {len(out)} scopes refolded; stored amplitudes reproduced to "
          f"{1000 * worst:.3g} mmag")


# ===========================================================================
# STAGE: colour — the colour curves as measurements (R6)
# ===========================================================================
def _colour_points(con: sqlite3.Connection, series_key: str,
                   fwhm_cut: bool = False) -> dict:
    """Catalogue-tied, unsaturated, cloud-clean target points of a series,
    exactly as Figure 6 selects them, optionally without poor-seeing
    frames."""
    sql = """
        SELECT l.bjd_tdb, l.cal_mag, f.night,
               f.fwhm_px * f.plate_scale AS fwhm
        FROM cv_lightcurve l
        JOIN cv_frames f ON f.frame_id = l.frame_id
                        AND f.series_key = l.series_key
        LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                                  AND c.series_key = l.series_key
        WHERE l.series_key = ? AND l.role = 'target'
          AND l.cal_mag IS NOT NULL AND l.saturated = 0
          AND COALESCE(c.vetoed, 0) = 0 ORDER BY l.bjd_tdb"""
    rows = con.execute(sql, (series_key,)).fetchall()
    if fwhm_cut:
        rows = [r for r in rows if r["fwhm"] is None
                or r["fwhm"] <= FWHM_CUT_ARCSEC]
    return {"t": np.array([r["bjd_tdb"] for r in rows], dtype=float),
            "m": np.array([r["cal_mag"] for r in rows], dtype=float),
            "night": np.array([r["night"] for r in rows])}


def _nearest_pairs(pa: dict, pb: dict, max_dt_s: float):
    """Nearest-neighbour pairing, as Figure 6 does it, with nights kept."""
    ta, tb = pa["t"], pb["t"]
    if ta.size == 0 or tb.size < 2:
        return np.array([]), np.array([]), np.array([]), np.array([])
    pos = np.clip(np.searchsorted(tb, ta), 1, tb.size - 1)
    left, right = pos - 1, pos
    near = np.where(np.abs(ta - tb[left]) <= np.abs(ta - tb[right]),
                    left, right)
    dt = (ta - tb[near]) * 86400.0
    keep = np.abs(dt) <= max_dt_s
    return (ta[keep], (pa["m"] - pb["m"][near])[keep], dt[keep],
            pa["night"][keep])


def cmd_colour(args) -> None:
    """Amplitude, extremum phases, repeatability and window survival."""
    con = connect(args.db)
    ensure_tables(con)
    eph = s9.load_ephemerides(con)[TARGET]
    print(f"database: {args.db}")
    for tab in ("rv_colour_curve", "rv_colour_summary", "rv_colour_repeat"):
        con.execute(f"DELETE FROM {tab}")
    clear_stage(con, "colour")
    tie = {r["series_key"]: r for r in con.execute(
        "SELECT * FROM cv_cattie WHERE is_primary = 1")}
    filt = {7: {"g": "G", "r": "R", "i": "I"},
            76: {"g": "g", "r": "r", "i": "i"}}
    methods = (
        ("pair600", "nearest exposure within 600 s (the first draft's "
                    "window)"),
        ("pair120", f"nearest exposure within {rv.STRICT_PAIR_WINDOW_S:.0f} "
                    "s (the window the referee asked for)"),
        ("interp", "bluer band linearly interpolated to the redder band's "
                   "exposure times between bracketing exposures no more "
                   "than two cadences apart: no first-order "
                   "non-simultaneity"),
        ("interp_fwhm", "the interpolated estimator with frames of seeing "
                        f"FWHM above {FWHM_CUT_ARCSEC:g} arcsec removed"),
    )
    summ = {}
    for era in (7, 76):
        for a, b in rv.BAND_PAIRS:
            colour = f"{a}-{b}"
            ska = f"{TARGET}|e{era}|{filt[era][a]}"
            skb = f"{TARGET}|e{era}|{filt[era][b]}"
            for method, mnote in methods:
                pa = _colour_points(con, ska, fwhm_cut=method.endswith("fwhm"))
                pb = _colour_points(con, skb, fwhm_cut=method.endswith("fwhm"))
                if pa["t"].size < 20 or pb["t"].size < 20:
                    continue
                if method.startswith("pair"):
                    win = 600.0 if method == "pair600" \
                        else rv.STRICT_PAIR_WINDOW_S
                    t, col, dt, ng = _nearest_pairs(pa, pb, win)
                    med_dt = float(np.median(dt)) if dt.size else None
                else:
                    # Per night, so nothing is interpolated across a gap
                    # between nights; the gap gate is twice the band-a
                    # cadence of that night.
                    ts, cs, ns = [], [], []
                    for night in np.unique(pb["night"]):
                        sa = pa["night"] == night
                        sb_ = pb["night"] == night
                        if sa.sum() < 3:
                            continue
                        cad = float(np.median(np.diff(pa["t"][sa]))
                                    * 86400.0)
                        tt, cc, _ = rv.interpolated_colour(
                            pa["t"][sa], pa["m"][sa], pb["t"][sb_],
                            pb["m"][sb_], 2.0 * cad)
                        ts.append(tt)
                        cs.append(cc)
                        ns.append(np.full(tt.size, night))
                    t = np.concatenate(ts) if ts else np.array([])
                    col = np.concatenate(cs) if cs else np.array([])
                    ng = np.concatenate(ns) if ns else np.array([])
                    med_dt = 0.0
                if t.size < 40:
                    continue
                ph = p3.phase_of(t, eph.period_d, eph.epoch_bjd)
                s = rv.colour_curve_summary(
                    ph, col, ng, seed=rv.SEED + era + len(method))
                ra = tie[ska]["check_rms_clip"] if ska in tie else None
                rb = tie[skb]["check_rms_clip"] if skb in tie else None
                ua = tie[ska]["check_rms"] if ska in tie else None
                ub = tie[skb]["check_rms"] if skb in tie else None
                sys_c = (math.hypot(ra, rb) if ra and rb else None)
                sys_r = (math.hypot(ua, ub) if ua and ub else None)
                con.execute("""
                    INSERT OR REPLACE INTO rv_colour_summary
                    (era_id, colour, method, n_pairs, n_nights,
                     n_bins_used, amplitude_mag, amplitude_err_mag,
                     phase_red, phase_red_err, phase_blue, phase_blue_err,
                     mean_mag, tie_sys_clip_mag, tie_sys_raw_mag,
                     median_dt_s, n_boot_ok, note)
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    era, colour, method, s["n_pairs"], s["n_nights"],
                    s["n_bins_used"], s["amplitude"], s["amplitude_err"],
                    s["phase_red"], s["phase_red_err"], s["phase_blue"],
                    s["phase_blue_err"], s["mean"], sys_c, sys_r, med_dt,
                    s["n_boot_ok"], mnote))
                for k in range(len(s["centres"])):
                    con.execute("""
                        INSERT OR REPLACE INTO rv_colour_curve
                        (era_id, colour, method, bin, phase, median_mag,
                         err_mag, n) VALUES (?,?,?,?,?,?,?,?)""", (
                        era, colour, method, k, float(s["centres"][k]),
                        _nz(s["median"][k]), _nz(s["bin_err"][k]),
                        int(s["count"][k])))
                summ[(era, colour, method)] = s
                if s["amplitude"] is not None:
                    err = s["amplitude_err"]
                    print(f"  era {era:2d} {colour} {method:12s} "
                          f"n={s['n_pairs']:4d} nights={s['n_nights']:2d} "
                          f"amp={s['amplitude']:.3f}"
                          + (f"+/-{err:.3f}" if err else " (no err)")
                          + f" red@{s['phase_red']:.2f} "
                          f"blue@{s['phase_blue']:.2f}")
    # --- repeatability between the two eras --------------------------------
    for a, b in rv.BAND_PAIRS:
        colour = f"{a}-{b}"
        for method, _ in methods:
            sa, sb_ = summ.get((7, colour, method)), summ.get(
                (76, colour, method))
            if not sa or not sb_ or sa["amplitude"] is None \
                    or sb_["amplitude"] is None:
                continue
            rep = rv.curve_repeatability(sa["median"], sa["bin_err"],
                                         sb_["median"], sb_["bin_err"])
            ea, eb = sa["amplitude_err"], sb_["amplitude_err"]
            pa_, pb_ = sa["phase_red_err"], sb_["phase_red_err"]
            con.execute("""
                INSERT OR REPLACE INTO rv_colour_repeat
                (colour, method, era_a, era_b, n_bins, rms_diff_mag, chi2,
                 dof, chi2nu, r, scale, amp_diff_mag, amp_diff_err_mag,
                 dphase_red, dphase_red_err)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                colour, method, 7, 76, rep["n_bins"], _nz(rep["rms_diff"]),
                _nz(rep["chi2"]), rep["dof"], _nz(rep["chi2nu"]),
                _nz(rep["r"]), _nz(rep["scale"]),
                sb_["amplitude"] - sa["amplitude"],
                (math.hypot(ea, eb)
                 if ea is not None and eb is not None else None),
                rv.circ_diff(sb_["phase_red"], sa["phase_red"]),
                (math.hypot(pa_, pb_)
                 if pa_ is not None and pb_ is not None else None)))
    con.commit()

    # --- scalars -----------------------------------------------------------
    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "colour", key, value, unit, fmt, note,
            origin=origin)

    mlab = {"pair600": "six hundred", "pair120": "strict", "interp":
            "interp", "interp_fwhm": "interp seeing"}
    for (era, colour, method), s in sorted(summ.items()):
        if s["amplitude"] is None:
            continue
        ctag = colour.replace("-", "")
        tag = f"rv col {ctag} era {era} {mlab[method]}"
        base = (f"ST LMi {colour} colour curve, {ERA_NAME[era]} era, "
                f"{dict(methods)[method]}")
        S(f"{tag} n", s["n_pairs"], "", "int", f"colour points in the {base}")
        S(f"{tag} nights", s["n_nights"], "", "int",
          "nights contributing those points")
        S(f"{tag} amp mag", s["amplitude"], "mag", "f2",
          f"amplitude of the {base}: mean of the three reddest phase bins "
          f"minus mean of the three bluest, of {rv.COLOUR_BINS} bins")
        S(f"{tag} amp err mag", s["amplitude_err"], "mag", "f2",
          "its error from a bootstrap over nights")
        S(f"{tag} phase red", s["phase_red"], "cycle", "f2",
          "orbital phase (catalogue ephemeris) at which the colour is "
          "reddest")
        S(f"{tag} phase red err", s["phase_red_err"], "cycle", "f2",
          "its night-bootstrap error")
        S(f"{tag} phase blue", s["phase_blue"], "cycle", "f2",
          "orbital phase at which the colour is bluest")
        S(f"{tag} phase blue err", s["phase_blue_err"], "cycle", "f2",
          "its night-bootstrap error")
    for r in con.execute("SELECT * FROM rv_colour_repeat"):
        ctag = r["colour"].replace("-", "")
        tag = f"rv col {ctag} repeat {mlab[r['method']]}"
        S(f"{tag} r", r["r"], "", "f2",
          f"Pearson correlation between the two eras' binned {r['colour']} "
          "curves (each with its own mean removed)")
        S(f"{tag} scale", r["scale"], "", "f2",
          "least-squares amplitude ratio of the 2025-26 curve to the 2024 "
          "curve (1 = same amplitude)")
        S(f"{tag} rms mag", r["rms_diff_mag"], "mag", "f2",
          "rms difference between the two eras' mean-removed curves")
        S(f"{tag} chinu", r["chi2nu"], "", "f2",
          "reduced chi-squared of that difference under the night-bootstrap "
          "bin errors")
        S(f"{tag} dof", r["dof"], "", "int", "its degrees of freedom")
        S(f"{tag} amp diff mag", r["amp_diff_mag"], "mag", "f2",
          "amplitude in 2025-26 minus amplitude in 2024")
        S(f"{tag} amp diff err mag", r["amp_diff_err_mag"], "mag", "f2",
          "night-bootstrap error of that amplitude difference")
        S(f"{tag} dphase", r["dphase_red"], "cycle", "f2",
          "phase of reddest colour in 2025-26 minus that in 2024")
        S(f"{tag} dphase err", r["dphase_red_err"], "cycle", "f2",
          "night-bootstrap error of that phase difference")
    # Tie systematic beside the amplitude (the editor's question to seat 3).
    r = con.execute("""SELECT max(tie_sys_clip_mag), max(tie_sys_raw_mag)
        FROM rv_colour_summary""").fetchone()
    S("rv col tie sys clip max mag", r[0], "mag", "f3",
      "largest catalogue-tie systematic (sigma-clipped check-star "
      "residuals of the two series in quadrature) on any colour zero "
      "point: a SHIFT of a whole curve, so it does not enter an amplitude")
    S("rv col tie sys raw max mag", r[1], "mag", "f3",
      "the same with every held-out check star kept")
    S("rv col strict window s", rv.STRICT_PAIR_WINDOW_S, "s", "int",
      "set by CV-R (revision_cv.STRICT_PAIR_WINDOW_S): the strict pairing "
      "window the referee asked the colour result to survive",
      origin="constant")
    stamp(con, "colour")
    con.commit()
    con.close()


# ===========================================================================
# STAGE: fitquality — the chi-squared the first draft did not print (R7)
# ===========================================================================
def cmd_fitquality(args) -> None:
    """Per-edge chi-squared distribution, with degrees of freedom."""
    con = connect(args.db)
    ensure_tables(con)
    print(f"database: {args.db}")
    con.execute("DELETE FROM rv_edge_quality")
    clear_stage(con, "fitquality")
    rows = con.execute("""SELECT estimator, band, era_id, chi2nu, dof,
        n_points FROM rv_edge WHERE accepted = 1 AND chi2nu IS NOT NULL"""
                       ).fetchall()

    def summarise(est, grp, sel):
        if not sel:
            return None
        c = np.array([r["chi2nu"] for r in sel])
        d = np.array([r["dof"] for r in sel])
        q = np.percentile(c, [25, 50, 75])
        con.execute("""INSERT OR REPLACE INTO rv_edge_quality VALUES
            (?,?,?,?,?,?,?,?,?,?,?,?)""", (
            est, grp, len(sel), int(d.min()), int(d.max()), float(c.min()),
            float(q[0]), float(q[1]), float(q[2]), float(c.max()),
            int(np.sum(c > 10)), int(np.sum(c < 0.5))))
        return {"n": len(sel), "median": float(q[1]), "max": float(c.max()),
                "min": float(c.min()), "above10": int(np.sum(c > 10)),
                "dofmin": int(d.min()), "dofmax": int(d.max())}
    out = {}
    for est in ("v1", "magc", "v2"):
        sel = [r for r in rows if r["estimator"] == est]
        out[est] = summarise(est, "all", sel)
        for b in ("g", "r", "i", "z"):
            summarise(est, f"band {b}", [r for r in sel if r["band"] == b])
        for era in (7, 76):
            summarise(est, f"era {era}", [r for r in sel
                                          if r["era_id"] == era])

    # --- the three fits the paper should show: best, median, worst --------
    # Chosen by rank of the published fit's reduced chi-squared, so the
    # choice is a rule and not a taste.  Their points are stored because
    # the figure must be redrawable from the release alone.
    con.execute("DELETE FROM rv_edge_example")
    eph = s9.load_ephemerides(con)[TARGET]
    ctx = series_context(con, eph)
    ranked = con.execute("""SELECT series_key, cycle, band, night, chi2nu,
        dof, t_edge_bjd, width_s, level_bright, step FROM rv_edge
        WHERE estimator='v1' AND band IN ('g','r','i')
        ORDER BY chi2nu""").fetchall()
    if ranked:
        picks = (("best", ranked[0]), ("median", ranked[len(ranked) // 2]),
                 ("worst", ranked[-1]))
        for label, r in picks:
            c = ctx[r["series_key"]]
            t_guess = (eph.epoch_bjd + (r["cycle"] + c["edge_phase"])
                       * eph.period_d)
            tt, mm, ee = rv.edge_window(c["t"], c["m"], c["e"], t_guess,
                                        eph.period_d)
            for k in range(tt.size):
                con.execute("""INSERT OR REPLACE INTO rv_edge_example
                    VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""", (
                    label, r["series_key"], r["cycle"], r["band"],
                    r["night"], r["chi2nu"], r["dof"], r["t_edge_bjd"],
                    r["width_s"], r["level_bright"], r["step"], k,
                    float((tt[k] - r["t_edge_bjd"]) * 86400.0),
                    float(mm[k]), float(ee[k])))

    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "fitquality", key, value, unit, fmt, note,
            origin=origin)

    for est, lab in (("v1", "published magnitude-space fits"),
                     ("v2", "flux-space fits")):
        s = out.get(est)
        if not s:
            continue
        S(f"rv fit {est} n", s["n"], "", "int", f"accepted edges, {lab}")
        S(f"rv fit {est} chinu median", s["median"], "", "f0",
          f"median reduced chi-squared of the {lab} under the inflated "
          "photometric errors")
        S(f"rv fit {est} chinu min", s["min"], "", "f1",
          f"smallest reduced chi-squared, {lab}")
        S(f"rv fit {est} chinu max", s["max"], "", "f0",
          f"largest reduced chi-squared, {lab}")
        S(f"rv fit {est} above ten", s["above10"], "", "int",
          f"edges with reduced chi-squared above 10, {lab}: the ramp does "
          "not describe flickering, which is why no per-edge formal error "
          "is published")
        S(f"rv fit {est} dof min", s["dofmin"], "", "int",
          "fewest degrees of freedom of any accepted edge fit (points in "
          "the window minus four parameters)")
        S(f"rv fit {est} dof max", s["dofmax"], "", "int",
          "most degrees of freedom of any accepted edge fit")
    S("rv fit window phase", p3.EDGE_WINDOW_PHASE, "cycle", "f2",
      "set by CV-S9 (phase3.EDGE_WINDOW_PHASE): half-width of the fitting "
      "window about the predicted edge",
      origin="constant")
    S("rv fit min points", p3.EDGE_MIN_POINTS, "", "int",
      "set by CV-S9 (phase3.EDGE_MIN_POINTS): fewest points a window may "
      "hold",
      origin="constant")
    S("rv fit min snr", p3.EDGE_MIN_SNR, "", "int",
      "set by CV-S9 (phase3.EDGE_MIN_SNR): smallest step, in units of the "
      "larger of the residual scatter and the median error, for an edge to "
      "be accepted",
      origin="constant")
    S("rv fit max bracket", p3.EDGE_MAX_BRACKET_CADENCE, "", "f1",
      "set by CV-S9 (phase3.EDGE_MAX_BRACKET_CADENCE): largest gap between "
      "the exposures bracketing the fitted epoch, in cadences, for an edge "
      "to be accepted",
      origin="constant")
    S("rv fit v one time nodes", 361, "", "int",
      "set by CV-S9 (run_cv_phase3.cmd_edges): time-grid nodes of the "
      "published fit, spanning three cadences either side of the ephemeris "
      "prediction",
      origin="constant")
    S("rv fit v one width nodes", 4, "", "int",
      "set by CV-S9 (run_cv_phase3.cmd_edges): ramp widths tried by the "
      "published fit, 0.5, 1, 2 and 4 times that series' folded-profile "
      "width",
      origin="constant")
    S("rv fit v two width nodes", len(rv.V2_WIDTH_GRID_S), "", "int",
      "set by CV-R (revision_cv.V2_WIDTH_GRID_S): ramp widths tried by the "
      "common-grid fits",
      origin="constant")
    S("rv fit v two width min s", rv.V2_WIDTH_GRID_S[0], "s", "f0",
      "set by CV-R (revision_cv.V2_WIDTH_GRID_S): narrowest ramp of the "
      "common width grid",
      origin="constant")
    S("rv fit v two width max s", rv.V2_WIDTH_GRID_S[-1], "s", "f0",
      "set by CV-R (revision_cv.V2_WIDTH_GRID_S): widest ramp of the "
      "common width grid",
      origin="constant")
    S("rv fit v two time step s", rv.V2_TIME_STEP_S, "s", "f0",
      "set by CV-R (revision_cv.V2_TIME_STEP_S): time-grid step of the "
      "common-grid fits",
      origin="constant")
    stamp(con, "fitquality")
    con.commit()
    con.close()
    if out.get("v1"):
        print(f"  v1: chi2nu median {out['v1']['median']:.0f}, max "
              f"{out['v1']['max']:.0f}, {out['v1']['above10']} of "
              f"{out['v1']['n']} above 10, dof {out['v1']['dofmin']}-"
              f"{out['v1']['dofmax']}")


# ===========================================================================
# STAGE: counts — what was staged, what was measured (R8)
# ===========================================================================
def cmd_counts(args) -> None:
    """Frame-count language, the S/N floor, EU UMa i/r, the seeing cut."""
    con = connect(args.db)
    ensure_tables(con)
    print(f"database: {args.db}")
    con.execute("DELETE FROM rv_counts")
    clear_stage(con, "counts")

    def C(scope, item, n, note):
        con.execute("INSERT OR REPLACE INTO rv_counts VALUES (?,?,?,?)",
                    (scope, item, int(n) if n is not None else None, note))

    def S(key, value, unit, fmt, note, origin="measured"):
        put(con, "counts", key, value, unit, fmt, note,
            origin=origin)

    # --- staged versus measured -------------------------------------------
    st = {r[0]: (r[1], r[2]) for r in con.execute(
        "SELECT status, count(*), count(DISTINCT night) FROM cv_frames "
        "GROUP BY status")}
    total = sum(v[0] for v in st.values())
    nights_all = con.execute("SELECT count(DISTINCT night) FROM cv_frames"
                             ).fetchone()[0]
    for k, (n, nn) in st.items():
        C("frames", k, n, f"cv_frames.status = {k}; on {nn} nights")
    S("rv frames staged", total, "", "int",
      "light frames STAGED for photometry (every row of cv_frames); the "
      "first draft called this number 'usable'")
    S("rv frames staged nights", nights_all, "", "int",
      "nights carrying a staged frame")
    S("rv frames measured", st.get("matched", (0, 0))[0], "", "int",
      "frames registered to the reference field and MEASURED "
      "(cv_frames.status = 'matched')")
    S("rv frames measured nights", st.get("matched", (0, 0))[1], "", "int",
      "nights carrying a measured frame")
    S("rv frames failed registration", st.get("failed_match", (0, 0))[0],
      "", "int", "staged frames on which registration to the reference "
      "field failed")
    S("rv frames off field", st.get("skipped_pointing", (0, 0))[0], "",
      "int", "staged frames that do not contain the reference field "
      "(pointing more than a degree off)")
    S("rv frames no pixels", st.get("excluded", (0, 0))[0], "", "int",
      "staged frames with no pixel file of the era's provenance")
    S("rv frames measured fraction", st.get("matched", (0, 0))[0] / total,
      "", "pct", "measured frames as a fraction of staged frames")
    # --- the S/N floor on "measurement" -----------------------------------
    # The measured chi2 inflation of each series' photometric errors, as
    # measured: not clipped at one (standing rule 1).
    infl = {r[0]: float(r[1]) if r[1] is not None and r[1] > 0 else 1.0
            for r in con.execute("SELECT series_key, inflation FROM "
                                 "cv_error_model")}
    pts = con.execute("""SELECT series_key, inst_mag_err FROM cv_lightcurve
        WHERE role = 'target' AND cal_mag IS NOT NULL""").fetchall()
    k = 2.5 / math.log(10.0)          # sigma_mag = k / (S/N)
    by_series: dict[str, list] = {}
    for sk, err in pts:
        by_series.setdefault(sk, []).append(
            (float(err) * infl.get(sk, 1.0)) if err is not None else np.nan)
    S("rv points tied", len(pts), "", "int",
      "catalogue-tied target points, the first draft's count of "
      "'measurements'")
    for floor in SNR_FLOORS:
        cut = k / floor
        n_pass = 0
        lost_series = []
        for sk, errs in sorted(by_series.items()):
            e = np.array(errs)
            ok = int(np.sum(np.isfinite(e) & (e <= cut)))
            n_pass += ok
            C(f"snr>={floor:g}", sk, ok,
              f"of {e.size} catalogue-tied target points, those with "
              f"inflated magnitude error at or below {cut:.3f} mag")
            if ok < e.size:
                lost_series.append((sk, e.size - ok, e.size))
        tag = f"rv points snr {int(floor)}"
        S(tag, n_pass, "", "int",
          f"catalogue-tied target points with S/N at or above {floor:g} "
          f"(inflated magnitude error at or below {cut:.3f} mag): what "
          "'measurement' means under that floor")
        S(f"{tag} lost", len(pts) - n_pass, "", "int",
          f"points the S/N >= {floor:g} floor removes")
        S(f"{tag} series touched", len(lost_series), "", "int",
          "series that lose at least one point to that floor")
    for filt in ("i", "r", "g"):
        sk = f"euuma|e76|{filt}"
        e = np.array(by_series.get(sk, []))
        S(f"rv euuma {filt} tied", e.size, "", "int",
          f"EU UMa {filt} catalogue-tied points")
        for floor in SNR_FLOORS:
            S(f"rv euuma {filt} snr {int(floor)}",
              int(np.sum(np.isfinite(e) & (e <= k / floor))), "", "int",
              f"of those, points with S/N at or above {floor:g}")
    # --- does the quoted precision range cover every counted series? ------
    # RF.M6: the abstract's precision range is taken over the series the
    # characterisation could measure a precision for, while the count of
    # "measurements" beside it includes series it could not.  Count them.
    if CHAR_DB.exists():
        ch = connect(CHAR_DB, read_only=True)
        prec = {r[0]: r[1] for r in ch.execute(
            "SELECT series_key, prec_at_target_all FROM ch_noise_series")}
        ch.close()
        k5 = k / 5.0
        no_prec = [sk for sk in by_series if prec.get(sk) is None]
        S("rv points series without precision", len(no_prec), "", "int",
          "series that contribute catalogue-tied target points to the "
          "first draft's measurement count but have NO per-point precision "
          "in ch_noise_series, and so lie outside the precision range "
          "quoted beside that count: " + (", ".join(sorted(no_prec))
                                          or "none"))
        S("rv points in series without precision",
          sum(len(by_series[sk]) for sk in no_prec), "", "int",
          "catalogue-tied points in those series")
        left = [sk for sk in no_prec if int(np.sum(
            np.isfinite(by_series[sk]) & (np.array(by_series[sk]) <= k5)))]
        S("rv points series without precision after floor", len(left), "",
          "int", "of those series, the ones that still contribute a point "
          "once the S/N >= 5 floor is applied: zero means the floor makes "
          "the precision range cover every counted series")
        covered = [1000.0 * prec[sk] for sk in by_series
                   if prec.get(sk) is not None and int(np.sum(
                       np.isfinite(by_series[sk])
                       & (np.array(by_series[sk]) <= k5)))]
        if covered:
            S("rv precision floor min mmag", min(covered), "mmag", "f0",
              "smallest per-point precision (ch_noise_series."
              "prec_at_target_all) among the series that contribute a "
              "point under the S/N >= 5 floor")
            S("rv precision floor max mmag", max(covered), "mmag", "f0",
              "largest per-point precision among those series")
    # --- the seeing cut ----------------------------------------------------
    r = con.execute("""SELECT count(*),
            sum(fwhm_px * plate_scale > ?) FROM cv_frames
        WHERE status = 'matched' AND fwhm_px IS NOT NULL""",
                    (FWHM_CUT_ARCSEC,)).fetchone()
    S("rv fwhm cut arcsec", FWHM_CUT_ARCSEC, "arcsec", "f0",
      "set by CV-R (run_cv_revision.FWHM_CUT_ARCSEC): the seeing cut "
      "tabulated, chosen against the 4 arcsec aperture radius",
      origin="constant")
    S("rv fwhm frames above", r[1], "", "int",
      f"measured frames with seeing FWHM above {FWHM_CUT_ARCSEC:g} arcsec")
    S("rv fwhm frames measured", r[0], "", "int",
      "measured frames with a FWHM")
    for tk, n in con.execute("""SELECT target_key, sum(fwhm_px *
            plate_scale > ?) FROM cv_frames WHERE status = 'matched'
            GROUP BY target_key""", (FWHM_CUT_ARCSEC,)):
        C("fwhm>cut", tk, n or 0, "measured frames above the seeing cut")
    n_pts = con.execute("""SELECT count(*) FROM cv_lightcurve l JOIN
        cv_frames f ON f.frame_id = l.frame_id AND f.series_key =
        l.series_key WHERE l.role = 'target' AND l.cal_mag IS NOT NULL AND
        f.fwhm_px * f.plate_scale > ?""", (FWHM_CUT_ARCSEC,)).fetchone()[0]
    S("rv fwhm points above", n_pts, "", "int",
      "catalogue-tied target points on those frames: what the seeing cut "
      "removes from the measurement count")
    n_edges = con.execute("SELECT count(*) FROM rv_edge WHERE estimator="
                          "'v1' AND n_fwhm_gt_cut > 0").fetchone()[0] \
        if con.execute("SELECT count(*) FROM rv_edge").fetchone()[0] else None
    S("rv fwhm edges touched", n_edges, "", "int",
      "published ST LMi edges whose fit window contains a frame above the "
      "seeing cut")
    mx = con.execute("SELECT max(fwhm_px * plate_scale) FROM cv_frames "
                     "WHERE status = 'matched'").fetchone()[0]
    S("rv fwhm max arcsec", mx, "arcsec", "f1",
      "worst seeing FWHM among measured frames")
    stamp(con, "counts")
    con.commit()
    con.close()
    print(f"  {total} staged, {st.get('matched', (0, 0))[0]} measured; "
          f"{len(pts)} tied points")


# ===========================================================================
# STAGE: figures
# ===========================================================================
def cmd_figures(args) -> None:
    """Draw the revision figures from the ``rv_`` tables."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from macro_phot import figures_cv as fx
    fx.apply_style()
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    PNG_DIR.mkdir(parents=True, exist_ok=True)
    cv = fx.connect_ro(args.db)
    out = connect(args.db)
    ensure_tables(out)
    wanted = args.only or list(fx.REVISION_FIGURE_IDS)
    for fig_id in wanted:
        entry = fx.REVISION_BUILDERS[fig_id]
        try:
            fig, spec = entry["fn"](cv)
        except Exception as exc:                            # noqa: BLE001
            print(f"  ! {fig_id}: FAILED, {type(exc).__name__}: {exc}")
            continue
        stem = f"{spec.fig_id}_{spec.label.split(':')[-1]}"
        pdf, png = PDF_DIR / f"{stem}.pdf", PNG_DIR / f"{stem}.png"
        fig.savefig(pdf, format="pdf")
        fig.savefig(png, format="png", dpi=fx.PNG_DPI)
        plt.close(fig)
        out.execute("""INSERT OR REPLACE INTO rv_figure (fig_id, label,
            title, caption, tables_used, pdf_path, png_path, width_in,
            finding, built_utc) VALUES (?,?,?,?,?,?,?,?,?,?)""", (
            spec.fig_id, spec.label, spec.title, spec.full_caption,
            ", ".join(spec.tables), str(pdf.relative_to(REPO_ROOT)),
            str(png.relative_to(REPO_ROOT)), spec.width_in, spec.note,
            utcnow()))
        out.commit()
        print(f"  + {spec.fig_id}  {spec.title}\n      {png}")
    stamp(out, "figures")
    out.close()


# ===========================================================================
# STAGE: status
# ===========================================================================
def cmd_status(args) -> None:
    con = connect(args.db, read_only=True)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name LIKE "
        "'rv_%'")}
    if not have:
        print("no rv_ tables: run a stage first")
        return
    for t in sorted(have):
        n = con.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
        print(f"  {t:20s} {n:6d} rows")
    print()
    for k, v in con.execute("SELECT key, value FROM rv_meta ORDER BY key"):
        print(f"  {k:28s} {v}")
    con.close()


def cmd_all(args) -> None:
    for fn in (cmd_edges, cmd_inject, cmd_offset, cmd_oc, cmd_longitude,
               cmd_superhump, cmd_yzrefold, cmd_colour, cmd_fitquality,
               cmd_counts,
               cmd_figures):
        print(f"\n=== {fn.__name__[4:]} ===", flush=True)
        fn(args)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--db", type=Path, default=DEFAULT_DB,
                       help="products database (default: %(default)s)")
        p.add_argument("--workers", type=int, default=MAX_WORKERS,
                       choices=range(1, MAX_WORKERS + 1), metavar="N",
                       help=f"worker processes, 1-{MAX_WORKERS}")
        p.add_argument("--n-real", type=int, default=rv.N_INJECT,
                       help="injection realisations per cell")
        p.add_argument("--n-boot", type=int, default=rv.N_EDGE_BOOT,
                       help="bootstrap replicates per edge")
        p.add_argument("--only", nargs="*", default=None,
                       help="figure ids to (re)build")
        return p

    for name, fn, hlp in (
            ("edges", cmd_edges, "three estimators + per-edge bootstrap"),
            ("inject", cmd_inject, "signed injection bias, every night"),
            ("offset", cmd_offset, "paired band-offset tests (R1, D3)"),
            ("oc", cmd_oc, "the O-C refit (R3)"),
            ("longitude", cmd_longitude, "spot longitude by state (R5)"),
            ("superhump", cmd_superhump, "excluded amplitudes (R4)"),
            ("yzrefold", cmd_yzrefold, "YZ Cnc on the 1994 period (R4)"),
            ("colour", cmd_colour, "colour curves as measurements (R6)"),
            ("fitquality", cmd_fitquality, "edge chi2 distribution (R7)"),
            ("counts", cmd_counts, "frame and measurement counts (R8)"),
            ("figures", cmd_figures, "draw the revision figures"),
            ("status", cmd_status, "print what is in the database"),
            ("all", cmd_all, "every stage, in order")):
        common(sub.add_parser(name, help=hlp)).set_defaults(func=fn)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
