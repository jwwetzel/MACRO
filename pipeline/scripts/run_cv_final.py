#!/usr/bin/env python
"""CV-S10 — the two closing science decisions of the CV project.

WHAT THIS SCRIPT DOES, AND WHY EACH PIECE IS HERE
--------------------------------------------------
Phase 3 answered the questions the archive was built to ask.  Two decisions
were left open, and both are DECISIONS -- they change what the manuscript is
allowed to claim, not what the sky did.

``hump``      YZ Cnc on the FALLBACK branch.  CV-S7 established that no
              dense run sits inside a superoutburst, so there is no
              superhump section.  This stage executes the fallback the
              strategy names instead: fold each QUIESCENT dense run on the
              published orbital period and measure the hump's semi-amplitude
              and phase -- with the detection decided by an
              injection-recovery contour measured on THAT run's timestamps
              and THAT run's noise, plus a bootstrap false-alarm threshold
              on the power at the orbital frequency.  Two tests, both
              measured, because either one alone can be fooled.

``flicker``   Flickering amplitude against timescale on the same runs, as a
              structure function.  The photometric noise floor subtracted
              from it is MEASURED, not modelled: the identical statistic
              computed on magnitude-matched field stars seen through the
              same frames.  Those stars carry the same photon noise, the
              same sky, the same ensemble zero-point wander and the same
              atmosphere, and no flickering -- so the quadrature difference
              is the intrinsic variability and nothing else.  YZ Cnc's four
              held-out CHECK stars are about a magnitude BRIGHTER than the
              star at quiescence, which is exactly why they are not used as
              the floor here.

``outburst``  The six dense runs inside NORMAL outbursts, characterised on
              their own terms: amplitude, covered duration, rate of change,
              and the blind-search injection contour that says what a
              coherent modulation would have had to reach on a single night
              to be measurable.  That last number is what makes "no
              superhump period" a measurement rather than an absence.

``gate``      The strategy's §4.19 gate, executed at last: 8 s High Gain
              frames at quiescent V ~ 14.5 might be sky- or read-noise
              dominated, and the fallback was not to be promised until
              somebody showed they are not.  The gate compares the
              precision actually achieved on the quiescent runs against the
              modulation the fallback exists to measure.

``anuma``     AN UMa, per filter.  CV-S5 already graded the three-filter
              colour goal NOT SUPPORTED (4 qualifying nights against a bar
              of 8).  This stage asks the narrower question: filter by
              filter, what can AN UMa support, and what is the deciding
              number?  Each capability is one measured number against one
              stated bar, so a reader who disagrees can disagree with the
              bar in one line instead of with a score nobody can audit.

``verdict``   The headline calls, each with its own deciding number, in the
              same shape CV-S5's ``ch_verdict`` uses.

All the arithmetic lives in ``macro_phot.final_science`` and the injection
machinery is reused unchanged from ``macro_phot.characterize`` -- the same
functions that produced the contours on the CV-S5 page, so the numbers on
the two pages are comparable by construction.  This script is I/O, staging
and bookkeeping.

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_final.py hump
    $P pipeline/scripts/run_cv_final.py flicker
    $P pipeline/scripts/run_cv_final.py outburst
    $P pipeline/scripts/run_cv_final.py gate
    $P pipeline/scripts/run_cv_final.py anuma
    $P pipeline/scripts/run_cv_final.py verdict
    $P pipeline/scripts/run_cv_final.py report
    $P pipeline/scripts/run_cv_final.py all

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite)
--------------------------------------------------------------
``p4_run``       per (dense run, filter): coverage, the folded hump, its
                 injection contour, its false-alarm threshold, the call.
``p4_flicker``   per (run, filter, timescale): the target's structure
                 function, the measured field-star floor, the excess.
``p4_outburst``  per outburst run: amplitude, rate, and the blind-search
                 contour that closes the superhump question.
``p4_gate``      the §4.19 signal-to-noise gate, line by line.
``p4_anuma``     per (filter, capability): measured value, bar, verdict.
``p4_verdict``   the headline decisions with their deciding numbers.
``p4_meta``      build stamps and every constant this run used.

RESUMABILITY
------------
Every stage deletes and rewrites only its own rows, keyed on the scope it
recomputes, so a stage may be re-run at any time.  Nothing here edits
``cv_lightcurve`` or any Phase-1/2/3 table.
"""

from __future__ import annotations

import argparse
import math
import re
import sqlite3
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import characterize as ch          # noqa: E402
from macro_phot import final_science as fs         # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
CHAR_DB = REPO_ROOT / "products" / "phot" / "cv_characterization.sqlite"

#: Stamped into ``p4_meta`` and read by the provenance graph.  Bump it when
#: the arithmetic changes, not when a comment does.
FINAL_CODE_VERSION = "CV-S10 v1.0 (2026-08-20, closing science decisions)"

BUSY_TIMEOUT_MS = 300_000

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}

# --------------------------------------------------------------------------
# YZ Cnc
# --------------------------------------------------------------------------

#: Arbitrary phase zero point, BJD_TDB.  It has to be arbitrary: VSX
#: publishes a period for YZ Cnc and NO epoch, which ``p3_ephemeris`` stores
#: as an explicit note.  Every phase on this page is therefore relative to
#: this constant, and the only phase statements that survive that are the
#: ones taken WITHIN a run -- the shape of the fold, and the offset between
#: two filters observed through the same night.  The ``phase_drift_d``
#: column in ``p4_run`` is the number that proves the point.
YZ_EPOCH_BJD = 2_460_000.0

#: Structure-function timescale bins, seconds.  The lower edge is well below
#: the ~145-215 s per-filter cadence on purpose: the empty short-lag bins
#: are the page's honest statement of what the sampling cannot reach, and
#: deleting them would let a reader assume the range was measured.
TAU_MIN_S, TAU_MAX_S, TAU_BINS = 60.0, 14_400.0, 9

#: Injection amplitudes and search band.  IDENTICAL to CV-S5's, so a
#: contour on this page and a contour on the characterization page are the
#: same statistic and may be compared directly.
INJECT_AMPS = np.array([0.002, 0.004, 0.007, 0.012, 0.020, 0.035,
                        0.060, 0.100, 0.170, 0.300])
SEARCH_FMIN_CD, SEARCH_FMAX_CD = 2.0, 40.0
SEARCH_FSTEP_CD = 0.0005
INJECT_TRIALS = 50
THRESHOLD_TRIALS = 300

#: The signal-to-noise ratio §4.19's gate demands between the modulation the
#: fallback measures and the per-point precision that measures it.  Five:
#: the same bar CV-S9's edge fitter uses to call a step distinguishable from
#: flickering, adopted here rather than invented so that two stages of this
#: pipeline do not hold two different opinions about what "detectable" means.
GATE_RATIO = 5.0

#: The ratio §4.19's floor line demands between the target's own variability
#: and the MEASURED photometric floor at the same timescale.  Two, not five:
#: this line asks whether the frames are noise-DOMINATED, and a signal twice
#: the noise is not dominated by it.  The stricter significance question --
#: is the excess real -- is asked separately, and by a sigma rather than by
#: a ratio, because a ratio has no error bar.
GATE_RATIO_SF = 2.0

# --------------------------------------------------------------------------
# AN UMa
# --------------------------------------------------------------------------

#: The bars every AN UMa capability is graded against.  ARGUE WITH THESE
#: NUMBERS, not with a score: each one is a single threshold, stated with
#: where it comes from, and moving it moves exactly one verdict.
ANUMA_BARS = {
    # A folded morphology panel needs the modulation to be detected at all
    # in that filter, and enough independent full-orbit nights that the fold
    # is not one night wearing a season's clothes.  Three is the smallest
    # number that can show reproducibility (two agreeing nights is a
    # coincidence; three is a pattern).
    "full_orbit_nights": 3,
    # A per-cycle O-C needs epochs.  Five accepted bright-phase edges is the
    # minimum that can show a trend rather than a scatter, and CV-S9's
    # cycle-count stage already refuses to build an O-C from four.
    "accepted_edges": 5,
    # sigma_t against the strategy's own 60 s threshold (§4.16), quoted in
    # the pessimistic shape-mismatched regime because the bright-phase edge
    # shape is unknown and colour dependent.
    "sigma_t_s": 60.0,
    # A duty cycle is a fraction of time, so its error bar is binomial in
    # the number of INDEPENDENT nights.  +/-15 percentage points is the
    # loosest interval that still distinguishes a 50% duty cycle from a 20%
    # or an 80% one; below that the number illustrates rather than measures.
    "duty_halfwidth_pp": 15.0,
    # State history in the relative sense -- which nights were high and
    # which low -- needs the nightly magnitudes to actually separate into
    # two populations.  CV-S9's own bimodality bar.
    "separability": 0.75,
    # The colour goal, restated from CV-S5's Q5 so the two pages cannot
    # drift apart.  It is a property of the SCHEDULE, not of any one filter.
    "three_filter_nights": 8,
}

#: Minimum points in a night before it counts as a full-orbit night in a
#: filter, and the span it must exceed (one orbital period).  Same rule as
#: CV-S5's Q1 census, restated here so this page's counts are reproducible
#: without opening that one.
FULL_ORBIT_MIN_POINTS = 12


# ===========================================================================
# Database plumbing
# ===========================================================================
def connect(path: Path, read_only: bool = False) -> sqlite3.Connection:
    """One connection, always with the long busy timeout.

    Other workflows write this archive; SQLite's own waiter is the only
    lock all of them agree about.
    """
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    con.row_factory = sqlite3.Row
    return con


def git_commit() -> str:
    """Short commit, ``-dirty`` when the tree is not clean.  A product
    stamped with a clean commit but built from edited files is a false
    reproducibility claim, so the suffix is not optional."""
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO_ROOT, capture_output=True, text=True,
                             timeout=20).stdout.strip() or "unknown"
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
    CREATE TABLE IF NOT EXISTS p4_meta (key TEXT PRIMARY KEY, value TEXT);

    CREATE TABLE IF NOT EXISTS p4_run (
        scope TEXT PRIMARY KEY, series_key TEXT, kind TEXT,
        nights TEXT, utc_nights TEXT, target_key TEXT,
        era_id INTEGER, filter TEXT, state TEXT, n_runs INTEGER,
        n_points INTEGER, span_h REAL, cycles REAL, cadence_s REAL,
        median_cal_mag REAL, amp_p5p95 REAL,
        sigma_point REAL, sigma_basis TEXT,
        hump_amp REAL, hump_amp_sigma REAL, hump_amp_harm REAL,
        hump_phase REAL, chi2nu REAL,
        amp90_field REAL, amp90_field_lo REAL, amp90_field_hi REAL,
        thr_field REAL,
        amp90_self REAL, amp90_self_lo REAL, amp90_self_hi REAL,
        thr_self REAL,
        power_forb REAL, n_floor_stars INTEGER, phase_drift_cycles REAL,
        detection TEXT, note TEXT);

    CREATE TABLE IF NOT EXISTS p4_flicker (
        series_key TEXT, night TEXT, state TEXT, filter TEXT,
        tau_s REAL, sf_target REAL, sf_target_sigma REAL, n_pairs INTEGER,
        sf_floor REAL, sf_floor_sigma REAL, n_floor_stars INTEGER,
        sf_excess REAL, excess_sigma REAL, detected INTEGER,
        sigma_formal REAL,
        PRIMARY KEY (series_key, night, tau_s));

    CREATE TABLE IF NOT EXISTS p4_outburst (
        series_key TEXT, night TEXT, utc_night TEXT, filter TEXT,
        n_points INTEGER, span_h REAL, median_cal_mag REAL,
        amp_above_quiescence REAL, amp_p5p95 REAL,
        rate_mag_per_h REAL, rate_sigma REAL, rate_verdict TEXT,
        amp90_blind REAL, amp90_blind_lo REAL, amp90_blind_hi REAL,
        superhump_floor REAL, episode TEXT, structure TEXT, note TEXT,
        PRIMARY KEY (series_key, night));

    CREATE TABLE IF NOT EXISTS p4_gate (
        gate_id TEXT, scope TEXT, quantity TEXT, value REAL, bar REAL,
        ratio REAL, passes INTEGER, unit TEXT, note TEXT,
        PRIMARY KEY (gate_id, scope));

    CREATE TABLE IF NOT EXISTS p4_anuma (
        filter TEXT, capability TEXT, rank INTEGER, measured REAL,
        bar REAL, unit TEXT, verdict TEXT, deciding_number TEXT,
        reasoning TEXT, what_would_change_it TEXT,
        PRIMARY KEY (filter, capability));

    CREATE TABLE IF NOT EXISTS p4_verdict (
        verdict_id TEXT PRIMARY KEY, rank INTEGER, task TEXT,
        question TEXT, verdict TEXT, deciding_number TEXT,
        reasoning TEXT, alternative TEXT);
    """)
    con.commit()


def set_meta(con: sqlite3.Connection, items: dict) -> None:
    con.executemany("INSERT OR REPLACE INTO p4_meta (key, value) VALUES (?,?)",
                    [(k, str(v)) for k, v in items.items()])
    con.commit()


def stamp(con: sqlite3.Connection, stage: str) -> None:
    set_meta(con, {f"stage_{stage}": utcnow(),
                   "final_code_version": FINAL_CODE_VERSION,
                   "git_commit": git_commit()})


def record_stage(key: str) -> None:
    """Tell the provenance graph this stage ran.  Never fatal."""
    try:
        subprocess.run([sys.executable,
                        str(PIPELINE_ROOT / "scripts" /
                            "check_pipeline_status.py"), "record", key],
                       cwd=REPO_ROOT, capture_output=True, text=True,
                       timeout=180)
    except Exception as exc:                            # noqa: BLE001
        print(f"  ! could not record stage {key}: {exc}")


def _f(x):
    """A float SQLite will accept, or NULL.  NaN is not a number a database
    column should carry: a NULL says 'not measured' and a NaN says
    'measured, and the answer was nonsense'."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


# ===========================================================================
# Reading the archive
# ===========================================================================
def ephemeris(con, target_key: str) -> dict:
    """The published ephemeris row, as this stage will use it.

    Read from ``p3_ephemeris`` rather than typed, so the period the fold
    uses is provably the period CV-S9 verified against and the manuscript
    cites.  For YZ Cnc the epoch is NULL and the note says why; this
    function returns it unchanged rather than filling in a default, and the
    caller supplies :data:`YZ_EPOCH_BJD` and says so.
    """
    r = con.execute("SELECT * FROM p3_ephemeris WHERE target_key=?",
                    (target_key,)).fetchone()
    if r is None:
        raise SystemExit(f"no p3_ephemeris row for {target_key}: run "
                         f"run_cv_phase3.py ephem first")
    return dict(r)


def dense_runs(con, target_key: str) -> list[dict]:
    """The dense runs CV-S7 identified, with their states."""
    rows = [dict(r) for r in con.execute(
        "SELECT * FROM cv_ext_verdict WHERE target=?", (target_key,))]
    return fs.select_dense_runs(rows)


def series_of(con, target_key: str) -> list[dict]:
    """Every solved, validated series for a target."""
    return [dict(r) for r in con.execute(
        "SELECT series_key, era_id, filter, chi2_inflation FROM cv_series "
        "WHERE target_key=? AND status='solved' AND target_verdict='measured' "
        "ORDER BY era_id, filter", (target_key,))]


def load_run(con, series_key: str, night: str) -> dict:
    """One night of one series: target points, cloud-vetoed frames removed.

    ``cal_mag`` throughout, exactly as Phase 3 does, and rows without one
    are dropped rather than falling back to the instrumental magnitude:
    mixing the two inside one run would put a zero-point step in the middle
    of a light curve that this stage then measures as variability.
    """
    rows = con.execute("""
        SELECT l.frame_id, l.bjd_tdb, l.cal_mag, l.inst_mag_err
        FROM cv_lightcurve l
        JOIN cv_frames f ON f.frame_id = l.frame_id
                        AND f.series_key = l.series_key
        LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                                  AND c.series_key = l.series_key
        WHERE l.series_key = ? AND f.night = ? AND l.role = 'target'
          AND l.cal_mag IS NOT NULL AND l.saturated = 0
          AND COALESCE(c.vetoed, 0) = 0
        ORDER BY l.bjd_tdb
    """, (series_key, night)).fetchall()
    t = np.array([r["bjd_tdb"] for r in rows], dtype=float)
    m = np.array([r["cal_mag"] for r in rows], dtype=float)
    e = np.array([r["inst_mag_err"] if r["inst_mag_err"] is not None
                  else np.nan for r in rows], dtype=float)
    # A missing or absurd error bar becomes the run's median rather than a
    # NaN that would silently drop the point out of every weighted fit.
    med = float(np.nanmedian(e)) if np.isfinite(e).any() else 0.02
    e = np.where(np.isfinite(e) & (e > 0), e, med)
    return {"frame_id": np.array([r["frame_id"] for r in rows],
                                 dtype=np.int64),
            "t": t, "m": m, "e": e}


def floor_pool(con, series_key: str, night: str, frame_ids, target_mag: float,
               half_width: float = fs.FLOOR_MATCH_HALF_WIDTH
               ) -> tuple[list[np.ndarray], int]:
    """Magnitude-matched field-star residual vectors on the run's frames.

    Returns ``(pool, n_stars)``.  Each vector is one star's magnitudes minus
    its own median, aligned to ``frame_ids``, with missing frames filled
    with zero so the vector length matches the timestamp vector exactly.
    Filling rather than dropping keeps the SAMPLING identical between the
    noise and the signal, which is the entire point of injecting at real
    timestamps -- and it is what makes these same vectors usable both as
    the injection noise pool and as the flickering floor.

    Two selection rules, both stated on the page:

    * within ``half_width`` magnitudes of the TARGET's median on this night,
      so the floor carries the target's photon statistics rather than the
      comparison ensemble's brighter ones;
    * covering at least :data:`macro_phot.final_science.FLOOR_MIN_COVERAGE`
      of the run's frames, because a star seen in a third of the frames has
      a structure function built from different sampling.
    """
    rows = con.execute("""
        SELECT l.star_id, l.frame_id, l.cal_mag
        FROM cv_lightcurve l
        JOIN cv_frames f ON f.frame_id = l.frame_id
                        AND f.series_key = l.series_key
        WHERE l.series_key = ? AND f.night = ? AND l.role IN ('comp','check')
          AND l.cal_mag IS NOT NULL AND l.saturated = 0
    """, (series_key, night)).fetchall()
    if not rows:
        return [], 0
    index = {int(f): i for i, f in enumerate(frame_ids)}
    by_star: dict[int, np.ndarray] = {}
    for r in rows:
        i = index.get(int(r["frame_id"]))
        if i is None:
            continue
        arr = by_star.setdefault(int(r["star_id"]),
                                 np.full(len(frame_ids), np.nan))
        arr[i] = float(r["cal_mag"])
    pool = []
    need = max(3, int(fs.FLOOR_MIN_COVERAGE * len(frame_ids)))
    for arr in by_star.values():
        ok = np.isfinite(arr)
        if int(ok.sum()) < need:
            continue
        med = float(np.median(arr[ok]))
        if abs(med - float(target_mag)) > float(half_width):
            continue
        pool.append(np.where(ok, arr - med, 0.0))
    return pool, len(pool)


# ===========================================================================
# STAGE hump — fold the dense runs and decide whether a hump is there
# ===========================================================================
def _contour(t, pool, period_d, score, freqs, rng_seed=20260820,
             detrend_nights=None):
    """A 90% recovery contour and its binomial error bar, one scope.

    Thin wrapper over the CV-S5 machinery so the two pages' contours are the
    same statistic.  ``score='known'`` asks "is a modulation detectable at a
    frequency the literature already gives?", which is the question the
    orbital hump poses; ``score='period'`` asks "could the period be
    MEASURED from scratch?", which is the question a superhump poses and the
    reason the two are never quoted interchangeably here.
    """
    at = 1.0 / period_d if score == "known" else None
    # The threshold must be measured through the SAME pipeline the
    # injections pass through: a fit that removes one constant per night has
    # a lower false-alarm floor, and reusing the raw threshold would credit
    # the analysis with an improvement it did not earn.
    thr_pool = ([ch.remove_nightly_means(p, detrend_nights) for p in pool]
                if detrend_nights is not None else pool)
    thr = ch.detection_threshold(t, thr_pool, freqs, THRESHOLD_TRIALS,
                                 ch.DETECT_FAP,
                                 np.random.default_rng(rng_seed),
                                 at_freq_cd=at)
    if not np.isfinite(thr):
        return float("nan"), float("nan"), float("nan"), float("nan")
    fr = [ch.recovery_fraction(t, pool, freqs, period_d, float(A), thr,
                               INJECT_TRIALS,
                               rng=np.random.default_rng(rng_seed + 7 * i),
                               score=score, detrend_nights=detrend_nights)
          for i, A in enumerate(INJECT_AMPS)]
    a90 = ch.recovery_contour(INJECT_AMPS, fr)
    lo, hi = ch.contour_uncertainty(INJECT_AMPS, fr, INJECT_TRIALS)
    return a90, lo, hi, thr


def _known_contour(t, pool, nights, period_d, detrend: bool,
                   rng_seed=20260820):
    """A known-period contour whose realizations roll NIGHT BY NIGHT.

    Same statistic as :func:`_contour` with ``score='known'`` -- 90%
    recovery of an injected sinusoid at the published orbital frequency,
    against a 1% false-alarm threshold measured on signal-free realizations
    of the same noise -- and the same amplitude grid, so the numbers are
    comparable with CV-S5's.

    The one difference is the realization: rolls happen within each night
    (:func:`macro_phot.final_science.roll_within_nights`) instead of across
    the whole vector.  On a single-night scope the two are identical; on a
    two-night block the whole-vector roll wraps one night's residuals onto
    the other's timestamps and invents a step at the join, which leaks into
    the orbital band and inflates the threshold.
    """
    f_orb = 1.0 / period_d
    rng = np.random.default_rng(rng_seed)
    pool = [np.asarray(p, dtype=float) for p in pool
            if np.isfinite(p).sum() >= 5]
    if not pool:
        return (float("nan"),) * 4

    def realize():
        base = pool[rng.integers(len(pool))]
        y = fs.roll_within_nights(base, nights, rng)
        return y

    pows = []
    for _ in range(THRESHOLD_TRIALS):
        y = realize()
        if detrend:
            y = ch.remove_nightly_means(y, nights)
        p = ch.ls_power_at(t, y, f_orb)
        if np.isfinite(p):
            pows.append(p)
    if not pows:
        return (float("nan"),) * 4
    thr = float(np.quantile(pows, 1.0 - ch.DETECT_FAP))
    fr = []
    for A in INJECT_AMPS:
        hits = 0
        for _ in range(INJECT_TRIALS):
            y = realize() + ch.inject_sinusoid(t, period_d, float(A),
                                               phase=float(rng.random()))
            if detrend:
                y = ch.remove_nightly_means(y, nights)
            p = ch.ls_power_at(t, y, f_orb)
            if np.isfinite(p) and p > thr:
                hits += 1
        fr.append(hits / float(INJECT_TRIALS))
    a90 = ch.recovery_contour(INJECT_AMPS, fr)
    lo, hi = ch.contour_uncertainty(INJECT_AMPS, fr, INJECT_TRIALS)
    return a90, lo, hi, thr


def _scopes(con, runs, series, period_d: float, sigma_p: float) -> list[dict]:
    """The units a hump is fitted to: single runs, plus foldable blocks.

    A SINGLE RUN is one night in one filter.  A BLOCK is several runs of the
    SAME state, SAME era and SAME filter folded together, and it exists
    because a single 2 h run cannot separate an orbital modulation from a
    trend while two runs a day apart can: a coherent signal adds across the
    gap and red noise does not.

    Two conditions gate a block, and both are arithmetic rather than taste:

    * SAME ERA.  Camera and epoch are perfectly confounded across the
      2024-05 seam (CV-S5's verdict S2), so a fold that crossed it would be
      measuring the instrument.
    * ACCUMULATED PHASE DRIFT within
      :data:`macro_phot.final_science.PHASE_DRIFT_BAR_CYCLES`.  YZ Cnc's
      published period carries only a quoted-precision floor of 5e-5 d, so
      over the 71 days between the February and May quiescent runs the phase
      is uncertain by about half a cycle -- those two may never share a
      phase axis.  Over the ONE day between 2024-05-01 and 2024-05-02 the
      drift is under 0.01 cycles, so those two may.
    """
    out = []
    for s in series:
        sk = s["series_key"]
        # A series belongs to exactly one era, so "same era" is enforced by
        # asking which of the dense nights this series actually observed.
        have = {str(r["night"]) for r in con.execute(
            "SELECT DISTINCT night FROM cv_frames WHERE series_key=?", (sk,))}
        mine = [r for r in runs if str(r["local_night"]) in have]
        for run in mine:
            out.append({"scope": f"{sk}|{run['local_night']}",
                        "kind": "run", "series": s, "runs": [run]})
        # Blocks are built only from QUIESCENT runs: the deliverable of this
        # branch is the quiescent hump, and folding across an outburst would
        # average light curves the star does not repeat.
        grp = [r for r in mine if str(r["state"]) == "QUIESCENT"]
        if len(grp) < 2:
            continue
        span = (max(_night_num(r) for r in grp)
                - min(_night_num(r) for r in grp))
        drift = fs.phase_drift_cycles(span, period_d, sigma_p)
        if drift > fs.PHASE_DRIFT_BAR_CYCLES:
            continue
        nights = "+".join(str(r["local_night"]) for r in grp)
        out.append({"scope": f"{sk}|block:{nights}", "kind": "block",
                    "series": s, "runs": grp, "drift": drift})
    return out


def _night_num(run) -> float:
    """A YYYY-MM-DD night as a day number, for span arithmetic."""
    y, m, d = (int(x) for x in str(run["local_night"]).split("-"))
    return float(datetime(y, m, d, tzinfo=timezone.utc).timestamp()) / 86400.0


def cmd_hump(args) -> None:
    """Fold YZ Cnc's dense runs on the published period and test the hump."""
    con = connect(args.db)
    ensure_tables(con)
    eph = ephemeris(con, "yzcnc")
    period_d = float(eph["period_d"])
    sigma_p = float(eph["period_sigma_d"] or 0.0)
    freqs = np.arange(SEARCH_FMIN_CD, SEARCH_FMAX_CD, SEARCH_FSTEP_CD)
    runs = dense_runs(con, "yzcnc")
    series = series_of(con, "yzcnc")
    infl = {r["series_key"]: (r["chi2_inflation"] or 1.0) for r in series}
    print(f"  YZ Cnc: {len(runs)} dense runs, P = {period_d} d "
          f"(sigma floor {sigma_p:g} d), epoch {YZ_EPOCH_BJD} (ARBITRARY)")
    con.execute("DELETE FROM p4_run WHERE target_key='yzcnc'")
    rows = []
    for sc in _scopes(con, runs, series, period_d, sigma_p):
        s = sc["series"]
        sk = s["series_key"]
        k = float(infl.get(sk, 1.0) or 1.0)
        parts = []
        for run in sc["runs"]:
            d = load_run(con, sk, str(run["local_night"]))
            if d["t"].size >= FULL_ORBIT_MIN_POINTS:
                parts.append((run, d))
        if not parts:
            continue
        t = np.concatenate([d["t"] for _r, d in parts])
        m = np.concatenate([d["m"] for _r, d in parts])
        e = np.concatenate([d["e"] for _r, d in parts])
        fids = np.concatenate([d["frame_id"] for _r, d in parts])
        nightlab = np.concatenate([np.full(d["t"].size,
                                           str(r["local_night"]))
                                   for r, d in parts])
        sig = e * k
        # Nightly constants are fitted JOINTLY with the harmonics for a
        # block (CV-S9's discipline), and a single constant for one run.
        fit = fs.fold_fit(t, m, sig, period_d, YZ_EPOCH_BJD,
                          night_index=(nightlab if sc["kind"] == "block"
                                       else None))
        # Coverage is counted in ORBITS ACTUALLY SAMPLED, not in the
        # baseline: a block spanning a day covers eleven orbits of which it
        # observed two, and only the two can constrain anything.
        cycles = sum(float(d["t"].max() - d["t"].min()) / period_d
                     for _r, d in parts)
        span_h = sum(float(d["t"].max() - d["t"].min()) * 24.0
                     for _r, d in parts)
        cadence = float(np.median(np.diff(np.sort(t)))) * 86400.0
        med_mag = float(np.median(m))
        pool: list[np.ndarray] = []
        n_floor = 0
        for run, d in parts:
            p_, n_ = floor_pool(con, sk, str(run["local_night"]),
                                d["frame_id"], float(np.median(d["m"])))
            n_floor += n_
            # Field-star vectors are per run; for a block they are laid into
            # the block's full frame vector so the injected signal and the
            # noise see identical sampling.
            index = {int(f): i for i, f in enumerate(fids)}
            for v in p_:
                full = np.zeros(fids.size)
                for j, f in enumerate(d["frame_id"]):
                    full[index[int(f)]] = v[j]
                pool.append(full)
        # TWO nulls, and the difference between them is itself a result.
        #
        # The FIELD null is CV-S5's: magnitude-matched field stars through
        # the same frames.  It answers "could the photometry see a hump this
        # size?" and knows nothing about the target.
        #
        # The SELF null is this star's own residuals about the fitted
        # orbital model, cyclically rolled.  It answers the question the
        # paper has to survive -- "could this star's OWN aperiodic
        # flickering have produced a peak this tall at the orbital
        # frequency?" -- and on a flickering dwarf nova it is much the
        # stricter.  The detection call is made against it; the field
        # contour is reported beside it as the instrumental floor, because a
        # reader needs to know which of the two is binding before deciding
        # what more telescope time would buy.
        block = sc["kind"] == "block"
        if len(pool) >= 2:
            a90f, lof, hif, thrf = _known_contour(t, pool, nightlab, period_d,
                                                  block)
        else:
            a90f = lof = hif = thrf = float("nan")
        # The self null is only meaningful where the fold model is NOT
        # degenerate with a trend.  On a sub-cycle scope the five harmonic
        # columns absorb nearly all of the variance, the residuals collapse,
        # and a null built from them would be absurdly permissive (it comes
        # out at 5-8 mmag, below the instrumental floor, which is by itself
        # proof that it is measuring the fit and not the star).  Those
        # scopes publish NULL rather than a flattering number.
        if cycles >= fs.MIN_CYCLES_FOR_DETECTION:
            self_pool = [np.where(np.isfinite(fit["resid"]), fit["resid"],
                                  0.0)]
            a90s, los, his, thrs = _known_contour(t, self_pool, nightlab,
                                                  period_d, block)
        else:
            a90s = los = his = thrs = float("nan")
        power = ch.ls_power_at(t, (ch.remove_nightly_means(m, nightlab)
                                   if block else m), 1.0 / period_d)
        call = fs.detection_call(fit["amp"], a90s, power, thrs, cycles)
        drift = sc.get("drift", 0.0)
        notes = []
        if len(pool) < 2:
            notes.append(f"only {len(pool)} field stars within "
                         f"{fs.FLOOR_MATCH_HALF_WIDTH} mag of the target "
                         f"cover this scope, so the instrumental contour is "
                         f"unmeasured (the self contour, which decides the "
                         f"call, is not)")
        if cycles < fs.MIN_CYCLES_FOR_DETECTION:
            notes.append(f"this scope samples only {cycles:.2f} orbits, "
                         f"below the {fs.MIN_CYCLES_FOR_DETECTION} at which "
                         f"a sinusoid at P_orb stops being degenerate with "
                         f"a trend; the amplitude is reported, the detection "
                         f"is not claimed")
        if sc["kind"] == "block":
            notes.append(f"{len(parts)} runs folded together: accumulated "
                         f"phase drift {drift:.4f} cycles, inside the "
                         f"{fs.PHASE_DRIFT_BAR_CYCLES} bar; one free "
                         f"constant per night fitted jointly with the "
                         f"harmonics")
        rows.append((
            sc["scope"], sk, sc["kind"],
            "+".join(str(r["local_night"]) for r, _d in parts),
            "+".join(str(r["utc_night"]) for r, _d in parts), "yzcnc",
            s["era_id"], s["filter"], str(parts[0][0]["state"]), len(parts),
            int(t.size), _f(span_h), _f(cycles), _f(cadence), _f(med_mag),
            _f(fs.percentile_amplitude(m)), _f(float(np.median(sig))),
            "per-point formal error x measured chi2 inflation "
            f"{k:.2f} (cv_series.chi2_inflation)",
            _f(fit["amp"]), _f(fit["amp_sigma"]), _f(fit["amp_harm"]),
            _f(fit["phase_max"]), _f(fit["chi2nu"]),
            _f(a90f), _f(lof), _f(hif), _f(thrf),
            _f(a90s), _f(los), _f(his), _f(thrs),
            _f(power), int(n_floor), _f(drift), call, "; ".join(notes)))
        print(f"    {sc['scope']:42s} [{str(parts[0][0]['state'])[:9]:9s}] "
              f"n={t.size:3d} cyc={cycles:4.2f} "
              f"A={1000 * fit['amp']:6.1f}+/-{1000 * fit['amp_sigma']:4.1f} "
              f"A90f={1000 * a90f:6.1f} A90s={1000 * a90s:6.1f} mmag  {call}")
    con.executemany(
        "INSERT OR REPLACE INTO p4_run VALUES (" + ",".join("?" * 36) + ")",
        rows)
    set_meta(con, {"yz_period_d": period_d, "yz_period_sigma_d": sigma_p,
                   "yz_epoch_bjd": YZ_EPOCH_BJD,
                   "yz_epoch_basis": "ARBITRARY — VSX publishes no epoch for "
                                     "YZ Cnc (see p3_ephemeris.note)",
                   "n_dense_runs": len(runs),
                   "inject_amps": ",".join(f"{a:g}" for a in INJECT_AMPS),
                   "inject_trials": INJECT_TRIALS,
                   "threshold_trials": THRESHOLD_TRIALS,
                   "detect_fap": ch.DETECT_FAP,
                   "n_harmonics": fs.N_HARMONICS,
                   "min_cycles_for_detection": fs.MIN_CYCLES_FOR_DETECTION,
                   "phase_drift_bar_cycles": fs.PHASE_DRIFT_BAR_CYCLES})
    stamp(con, "hump")
    con.commit()
    con.close()


# ===========================================================================
# STAGE flicker — amplitude against timescale, over a measured floor
# ===========================================================================
def cmd_flicker(args) -> None:
    """Structure functions for every dense run, with a measured noise floor."""
    con = connect(args.db)
    ensure_tables(con)
    eph = ephemeris(con, "yzcnc")
    period_d = float(eph["period_d"])
    edges = fs.log_tau_edges(TAU_MIN_S, TAU_MAX_S, TAU_BINS)
    runs = dense_runs(con, "yzcnc")
    series = series_of(con, "yzcnc")
    infl = {r["series_key"]: (r["chi2_inflation"] or 1.0) for r in series}
    con.execute("DELETE FROM p4_flicker")
    rows = []
    for run in runs:
        night, state = str(run["local_night"]), str(run["state"])
        for s in series:
            sk = s["series_key"]
            d = load_run(con, sk, night)
            if d["t"].size < FULL_ORBIT_MIN_POINTS:
                continue
            med_mag = float(np.median(d["m"]))
            k = float(infl.get(sk, 1.0) or 1.0)
            # The orbital model is REMOVED before the structure function is
            # computed.  A coherent hump has a structure function of its own
            # that rises across exactly the timescales flickering occupies,
            # so leaving it in would let the orbit masquerade as flickering.
            fit = fs.fold_fit(d["t"], d["m"], d["e"] * k, period_d,
                              YZ_EPOCH_BJD)
            resid = fit["resid"]
            centres, sf_t, npair = fs.structure_function(d["t"], resid, edges)
            sig_t = fs.sf_sigma(sf_t, npair, int(d["t"].size))
            pool, n_floor = floor_pool(con, sk, night, d["frame_id"], med_mag)
            curves = [fs.structure_function(d["t"], p, edges)[1] for p in pool]
            if n_floor >= fs.FLOOR_MIN_STARS:
                sf_f, sig_f = fs.median_floor(curves)
            else:
                sf_f = np.full(centres.size, np.nan)
                sig_f = np.full(centres.size, np.nan)
            exc = fs.quadrature_excess(sf_t, sf_f)
            z = fs.excess_significance(sf_t, sig_t, sf_f, sig_f)
            formal = float(np.median(d["e"] * k))
            for i, tau in enumerate(centres):
                det = (1 if np.isfinite(z[i]) and z[i] >= fs.FLICKER_SIGMA_BAR
                       else 0)
                rows.append((sk, night, state, s["filter"], _f(tau),
                             _f(sf_t[i]), _f(sig_t[i]), int(npair[i]),
                             _f(sf_f[i]), _f(sig_f[i]), int(n_floor),
                             _f(exc[i]), _f(z[i]), det, _f(formal)))
            good = np.isfinite(exc)
            print(f"    {sk:14s} {night} [{state[:9]:9s}] floor stars="
                  f"{n_floor:2d}  bins with excess: {int(good.sum())}/"
                  f"{centres.size}  formal={1000 * formal:.1f} mmag")
    con.executemany(
        "INSERT OR REPLACE INTO p4_flicker VALUES (" + ",".join("?" * 15)
        + ")", rows)
    set_meta(con, {"tau_edges_s": ",".join(f"{e:.1f}" for e in edges),
                   "floor_match_half_width": fs.FLOOR_MATCH_HALF_WIDTH,
                   "floor_min_coverage": fs.FLOOR_MIN_COVERAGE,
                   "floor_min_stars": fs.FLOOR_MIN_STARS,
                   "sf_min_pairs": fs.SF_MIN_PAIRS,
                   "flicker_sigma_bar": fs.FLICKER_SIGMA_BAR})
    stamp(con, "flicker")
    con.commit()
    con.close()


# ===========================================================================
# STAGE outburst — the six normal-outburst runs, on their own terms
# ===========================================================================
def cmd_outburst(args) -> None:
    """Characterise the normal-outburst dense runs and close the superhump
    question with a measured blind-search contour."""
    con = connect(args.db)
    ensure_tables(con)
    eph = ephemeris(con, "yzcnc")
    period_d = float(eph["period_d"])
    freqs = np.arange(SEARCH_FMIN_CD, SEARCH_FMAX_CD, SEARCH_FSTEP_CD)
    runs = [r for r in dense_runs(con, "yzcnc")
            if str(r["state"]) == "OUTBURST"]
    series = series_of(con, "yzcnc")
    con.execute("DELETE FROM p4_outburst")
    rows = []
    for run in runs:
        night = str(run["local_night"])
        # ONE filter per run carries the blind-search contour: the richest.
        # The contour is a statement about the SAMPLING, and running it in
        # three filters that share the same night would triple the compute
        # to produce three numbers that differ only by their noise.
        loaded = []
        for s in series:
            d = load_run(con, s["series_key"], night)
            if d["t"].size >= FULL_ORBIT_MIN_POINTS:
                loaded.append((s, d))
        if not loaded:
            continue
        best = max(loaded, key=lambda sd: sd[1]["t"].size)
        for s, d in loaded:
            sk = s["series_key"]
            span_h = float(d["t"].max() - d["t"].min()) * 24.0
            med_mag = float(np.median(d["m"]))
            rate, rate_sig = fs.linear_rate_per_hour(d["t"], d["m"])
            # A rate is only a direction if it is bigger than its own bar.
            if not np.isfinite(rate) or not np.isfinite(rate_sig):
                rverd = "UNMEASURED"
            elif abs(rate) < 3.0 * rate_sig:
                rverd = "FLAT (within 3 sigma)"
            elif rate > 0:
                rverd = "FADING"
            else:
                rverd = "BRIGHTENING"
            a90 = lo = hi = float("nan")
            if s is best[0]:
                pool, n_floor = floor_pool(con, sk, night, d["frame_id"],
                                           med_mag)
                if len(pool) >= 2:
                    a90, lo, hi, _thr = _contour(d["t"], pool, period_d,
                                                 "period", freqs)
            structure = (f"p5-p95 {fs.percentile_amplitude(d['m']):.3f} mag "
                         f"over {span_h:.2f} h at "
                         f"{span_h / 24.0 / period_d:.2f} orbital cycles")
            rows.append((sk, night, str(run["utc_night"]), s["filter"],
                         int(d["t"].size), _f(span_h), _f(med_mag),
                         _f(run["amp"]), _f(fs.percentile_amplitude(d["m"])),
                         _f(rate), _f(rate_sig), rverd,
                         _f(a90), _f(lo), _f(hi),
                         _f(fs.SUPERHUMP_SEMI_AMP_FLOOR),
                         str(run["episode"] or ""), structure,
                         "blind-search contour measured in this filter"
                         if s is best[0] else
                         "contour measured in this run's richest filter"))
            print(f"    {sk:14s} {night} n={d['t'].size:3d} "
                  f"rate={rate:+.4f}+/-{rate_sig:.4f} mag/h {rverd}"
                  + (f"  A90(blind)={1000 * a90:6.1f} mmag"
                     if np.isfinite(a90) else ""))
    con.executemany(
        "INSERT OR REPLACE INTO p4_outburst VALUES (" + ",".join("?" * 19)
        + ")", rows)
    set_meta(con, {"superhump_semi_amp_floor": fs.SUPERHUMP_SEMI_AMP_FLOOR,
                   "n_outburst_runs": len(runs)})
    stamp(con, "outburst")
    con.commit()
    con.close()


# ===========================================================================
# STAGE gate — §4.19's signal-to-noise check, executed
# ===========================================================================
def cmd_gate(args) -> None:
    """Decide §4.19: are 8 s High Gain frames at quiescence good enough?

    The strategy would not promise the quiescent fallback until somebody
    showed that 8 s High Gain frames at quiescent V ~ 14.5 are not sky- or
    read-noise dominated.  "Not dominated" is not measurable as written, so
    it is turned into three arithmetic lines, each answering a different
    half of what the fallback needs:

    ``4.19-floor``    the MEASURED photometric floor at the shortest lag the
                      run samples -- the field-star structure function, not
                      a formal error bar -- against the target's own
                      variability at the same lag.  This is the literal
                      question §4.19 asks.
    ``4.19-flicker``  how many timescale bins clear the significance bar.
                      A floor that is merely acceptable is worth nothing if
                      no bin survives it.
    ``4.19-hump``     the fitted hump against the INSTRUMENTAL contour.
                      This line asks only whether the photometry could see a
                      modulation this size; whether the STAR's own
                      flickering allows it to be called a detection is a
                      separate question and is answered in p4_run, not here.
    """
    con = connect(args.db)
    ensure_tables(con)
    con.execute("DELETE FROM p4_gate")
    rows = []
    runs = [dict(r) for r in con.execute(
        "SELECT * FROM p4_run WHERE target_key='yzcnc' AND state='QUIESCENT' "
        "AND kind='run' ORDER BY nights, filter")]
    if not runs:
        raise SystemExit("no quiescent runs in p4_run: run `hump` first")
    for r in runs:
        sk, night = r["series_key"], r["nights"]
        # The SHORTEST lag this run actually measures, which is a property
        # of the per-filter cadence and not of the exposure time: three
        # filters interleaving at 8 s each still sample the target only
        # every ~190 s.
        f = con.execute(
            "SELECT tau_s, sf_target, sf_floor, excess_sigma, sf_excess "
            "FROM p4_flicker WHERE series_key=? AND night=? "
            "AND sf_target IS NOT NULL AND sf_floor IS NOT NULL "
            "ORDER BY tau_s LIMIT 1", (sk, night)).fetchone()
        if f is None:
            continue
        ratio = (f["sf_target"] / f["sf_floor"]) if f["sf_floor"] else \
            float("nan")
        # The exposure time is the thing §4.19 names, so it is quoted from
        # the frames rather than assumed from the era.
        exp = con.execute("SELECT exptime FROM cv_frames WHERE series_key=? "
                          "AND night=? LIMIT 1", (sk, night)).fetchone()
        exp_s = float(exp[0]) if exp and exp[0] is not None else float("nan")
        rows.append(("4.19-floor", f"{sk}|{night}",
                     "target variability / measured photometric floor, at "
                     "the shortest lag sampled",
                     _f(ratio), GATE_RATIO_SF, _f(ratio),
                     int(np.isfinite(ratio) and ratio >= GATE_RATIO_SF),
                     "ratio",
                     f"at tau = {f['tau_s']:.0f} s the target's structure "
                     f"function is {1000 * f['sf_target']:.0f} mmag and the "
                     f"magnitude-matched field stars' is "
                     f"{1000 * f['sf_floor']:.0f} mmag "
                     f"({r['n_floor_stars']} stars, exposure "
                     f"{exp_s:.0f} s)"))
        n_ok = con.execute(
            "SELECT count(*) FROM p4_flicker WHERE series_key=? AND night=? "
            "AND detected=1", (sk, night)).fetchone()[0]
        n_meas = con.execute(
            "SELECT count(*) FROM p4_flicker WHERE series_key=? AND night=? "
            "AND sf_excess IS NOT NULL", (sk, night)).fetchone()[0]
        rows.append(("4.19-flicker", f"{sk}|{night}",
                     "timescale bins with variability above the floor",
                     float(n_ok), 1.0, float(n_ok), int(n_ok >= 1), "bins",
                     f"{n_ok} of {n_meas} measurable bins clear "
                     f"{fs.FLICKER_SIGMA_BAR:g} sigma on the variance "
                     f"excess"))
        a, c = r["hump_amp"], r["amp90_field"]
        ok = bool(a and c and a >= c)
        # The self contour is deliberately NULL on sub-cycle scopes (see
        # `hump`), so the sentence says that rather than printing a NaN.
        self_txt = (f"{1000 * r['amp90_self']:.0f} mmag"
                    if r["amp90_self"] is not None else
                    "not measurable on this scope, which samples fewer than "
                    f"{fs.MIN_CYCLES_FOR_DETECTION} orbits")
        rows.append(("4.19-hump", f"{sk}|{night}",
                     "fitted hump semi-amplitude / instrumental 90% "
                     "recovery contour",
                     _f((a / c) if (a and c) else float("nan")), 1.0,
                     _f((a / c) if (a and c) else float("nan")), int(ok),
                     "ratio",
                     f"hump {1000 * (a or float('nan')):.0f} mmag against an "
                     f"instrumental contour of "
                     f"{1000 * (c or float('nan')):.0f} mmag; the contour "
                     f"set by the star's OWN flickering is {self_txt} "
                     f"and is reported in p4_run, not here"))
    con.executemany(
        "INSERT OR REPLACE INTO p4_gate VALUES (?,?,?,?,?,?,?,?,?)", rows)
    set_meta(con, {"gate_ratio_sf": GATE_RATIO_SF,
                   "gate_ratio": GATE_RATIO})
    stamp(con, "gate")
    con.commit()
    for gid in ("4.19-floor", "4.19-flicker", "4.19-hump"):
        sel = [r for r in rows if r[0] == gid]
        print(f"    {gid:14s} {sum(1 for r in sel if r[6])}/{len(sel)} "
              f"scopes pass")
    con.close()


# ===========================================================================
# STAGE anuma — the per-filter go/no-go
# ===========================================================================
def full_orbit_nights(con, series_key: str, period_d: float) -> list[str]:
    """Nights in this filter whose span exceeds one orbit in >= 12 points.

    The same rule CV-S5's Q1 census applied, restated here in one place so
    that this page's per-filter counts are reproducible from this file
    alone.  The counts it returns reproduce ``ch_cadence.n_blocks_ge1cycle``
    exactly, which is the cross-check that the rule really is the same one.
    """
    rows = con.execute("""
        SELECT f.night AS night, l.bjd_tdb AS t
        FROM cv_lightcurve l
        JOIN cv_frames f ON f.frame_id = l.frame_id
                        AND f.series_key = l.series_key
        LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                                  AND c.series_key = l.series_key
        WHERE l.series_key = ? AND l.role = 'target'
          AND l.cal_mag IS NOT NULL AND l.saturated = 0
          AND COALESCE(c.vetoed, 0) = 0
        ORDER BY l.bjd_tdb
    """, (series_key,)).fetchall()
    by: dict[str, list[float]] = {}
    for r in rows:
        by.setdefault(str(r["night"]), []).append(float(r["t"]))
    return sorted(n for n, v in by.items()
                  if len(v) >= FULL_ORBIT_MIN_POINTS
                  and (max(v) - min(v)) > period_d)


def cmd_anuma(args) -> None:
    """Grade AN UMa capability by capability, filter by filter."""
    con = connect(args.db)
    ensure_tables(con)
    eph = ephemeris(con, "anuma")
    period_d = float(eph["period_d"])
    series = series_of(con, "anuma")
    con.execute("DELETE FROM p4_anuma")
    per_filter_nights = {}
    for s in series:
        per_filter_nights[s["filter"]] = full_orbit_nights(
            con, s["series_key"], period_d)
    # The colour goal is a property of the SCHEDULE: the nights on which all
    # three filters independently covered a full orbit.
    three = set.intersection(*[set(v) for v in per_filter_nights.values()]) \
        if per_filter_nights else set()
    # The MEASURED per-epoch timing error for this target, from the edges
    # actually fitted.  There is no injection test for AN UMa -- the grid
    # was run on ST LMi's Mode0 night -- so the analytic budget must not
    # stand in for one, and these are the numbers the capability rows quote.
    _es = con.execute(
        "SELECT min(sigma_t_s), max(sigma_t_s), count(*) FROM p3_edge "
        "WHERE target_key='anuma'").fetchone()
    edge_sig_lo, edge_sig_hi, n_edge_all = (_es[0] or float("nan"),
                                            _es[1] or float("nan"),
                                            int(_es[2] or 0))
    _eo = con.execute(
        "SELECT min(sigma_t_s), max(sigma_t_s), count(*) FROM p3_edge "
        "WHERE target_key='anuma' AND accepted=1").fetchone()
    edge_ok_lo, edge_ok_hi, n_edge_ok = (_eo[0] or float("nan"),
                                         _eo[1] or float("nan"),
                                         int(_eo[2] or 0))
    # Seasons actually observed, from the nights themselves — the rate that
    # converts a shortfall in nights into a number of observing seasons.
    all_nights = sorted({n for v in per_filter_nights.values() for n in v})
    seasons = len({_season(n) for n in all_nights}) or 1
    rows = []
    for s in series:
        filt = s["filter"]
        sk = s["series_key"]
        nights = per_filter_nights[filt]
        # ---- 1. folded morphology -------------------------------------
        per = dict(con.execute("SELECT * FROM p3_period WHERE series_key=?",
                               (sk,)).fetchone())
        detected = int(per["detected"] or 0)
        n_full = len(nights)
        # A full-orbit night is only usable for a folded panel if the
        # modulation is detected in THIS filter, so the count that is graded
        # is gated on detection.  Reporting 6 nights against a bar of 3 and
        # then overruling it with a footnote would put a passing number and
        # a failing verdict on the same row.
        usable = n_full if detected else 0
        v = fs.capability_verdict(usable, ANUMA_BARS["full_orbit_nights"])
        need = fs.nights_needed(n_full, ANUMA_BARS["full_orbit_nights"],
                                n_full / seasons if seasons else 0.0)
        rows.append((filt, "folded orbital morphology", 1, float(usable),
                     float(ANUMA_BARS["full_orbit_nights"]), "nights", v,
                     f"{usable} usable full-orbit nights in {filt} "
                     f"({n_full} nights whose SAMPLING covers an orbit by "
                     f"the Section 2.2 census -- >= "
                     f"{FULL_ORBIT_MIN_POINTS} points "
                     f"spanning more than one {period_d:.6f} d orbit -- "
                     f"of which usable counts only those where the "
                     f"modulation is also detected in this filter, which "
                     f"is why the two counts differ); CV-S9 modulation "
                     f"detected = {detected}, PDM theta "
                     f"{per['pdm_theta']:.3f}, LS power at the published "
                     f"frequency {per['p_ls_pow']:.3f}, fitted semi-amplitude "
                     f"{1000 * (per['amplitude_mag'] or 0):.0f} mmag",
                     "A folded panel needs the modulation to be detected in "
                     "THIS filter and enough independent full-orbit nights "
                     "that the fold shows reproducibility rather than one "
                     "night's weather.",
                     "" if v == "SUPPORTED" else
                     (f"{need['shortfall']} more full-orbit {filt} nights "
                      f"(~{need['seasons']:.1f} seasons at this archive's "
                      f"delivered rate)" if detected else
                      "nothing schedulable: the modulation is not detected "
                      "in this filter at all, so more nights of the same "
                      "kind buy depth, not a detection")))
        # ---- 2. bright-phase timing ------------------------------------
        # NOT "per-cycle" timing: §4.2 abolished the per-cycle epoch
        # programme-wide, so grading AN UMa on an estimator the paper says
        # it does not use put a passing standard and a failing rule on the
        # same page.  The question for this target is whether per-NIGHT
        # epochs are constructible, and they are not -- for two independent
        # reasons, both quoted here rather than one.
        edg = con.execute(
            "SELECT count(*) n, sum(accepted) ok FROM p3_edge "
            "WHERE series_key=?", (sk,)).fetchone()
        n_ok = int(edg["ok"] or 0)
        v = fs.capability_verdict(n_ok, ANUMA_BARS["accepted_edges"])
        cyc = con.execute("SELECT verdict, phase_spread FROM p3_cycle_count "
                          "WHERE target_key='anuma'").fetchone()
        # WHY the edges are rejected, from the stored reasons rather than
        # from an assumption about the cadence.  16 of 20 are step-SNR
        # failures; exactly one is a cadence gap.  A faster filter cycle
        # does not fix an S/N-limited edge, and the remedy below says so.
        rej = con.execute(
            "SELECT reason FROM p3_edge WHERE target_key='anuma' AND "
            "accepted=0").fetchall()
        n_snr = sum(1 for r in rej if "step SNR" in str(r["reason"]))
        n_grid = sum(1 for r in rej if "search grid" in str(r["reason"]))
        n_gap = sum(1 for r in rej if " gap," in str(r["reason"]))
        snr_vals = [float(m.group(1)) for m in
                    (re.search(r"step SNR ([0-9.]+)", str(r["reason"]))
                     for r in rej) if m]
        rows.append((filt, "bright-phase timing (O-C)", 2,
                     float(n_ok), float(ANUMA_BARS["accepted_edges"]),
                     "accepted edges", v,
                     f"{n_ok} of {int(edg['n'] or 0)} fitted edges accepted "
                     f"in {filt}. Across the target, {n_snr} of "
                     f"{len(rej)} rejections are step signal-to-noise below "
                     f"5 (measured {min(snr_vals):.1f}-{max(snr_vals):.1f}), "
                     f"{n_grid} are grid-edge epochs and {n_gap} is a "
                     f"cadence gap: the edge is not distinguishable from "
                     f"AN UMa's own flickering on most cycles. Per-NIGHT "
                     f"epochs are unavailable for a second, independent "
                     f"reason: the accepted edges scatter over "
                     f"{float(cyc['phase_spread']):.3f} in orbital phase "
                     f"against the 0.05 one-feature bar, so CV-S9's "
                     f"cycle-count verdict is "
                     f"{cyc['verdict'] if cyc else 'n/a'}",
                     f"The edge itself is only detected on a handful of "
                     f"cycles, so there are almost no epochs to time. The "
                     f"precision is NOT measured for this target: CV-S5's "
                     f"injection test was run on ST LMi's Mode0 night only, "
                     f"and CV-S5's ANALYTIC budget for AN UMa's richest g "
                     f"night (21-35 s) is the estimator §4.16 forbids "
                     f"quoting as achieved. What is measured is the "
                     f"per-epoch error of the fitted edges themselves: "
                     f"{edge_sig_lo:.0f}-{edge_sig_hi:.0f} s over all "
                     f"{n_edge_all} fitted edges, "
                     f"{edge_ok_lo:.0f}-{edge_ok_hi:.0f} s over the "
                     f"{n_edge_ok} accepted ones, every one of them outside "
                     f"the 60 s bar.",
                     "Whatever raises the STEP SIGNAL-TO-NOISE of the edge "
                     "-- longer exposures, fewer filters per cycle so each "
                     "gets more integration, or catching the star in a "
                     "brighter state -- because "
                     f"{n_snr} of {len(rej)} rejections are S/N and only "
                     f"{n_gap} is a cadence gap. A faster filter cycle "
                     "raises the number of sampled edges but not the "
                     "signal-to-noise of any one of them, so on its own it "
                     "would not convert a rejection into an acceptance."),)
        # ---- 3. relative state history ---------------------------------
        st = con.execute("SELECT * FROM p3_state_series WHERE series_key=?",
                         (sk,)).fetchone()
        sep = float(st["separability"]) if st and st["separability"] else \
            float("nan")
        v = fs.capability_verdict(sep, ANUMA_BARS["separability"])
        rows.append((filt, "relative state history",
                     3, _f(sep) or float("nan"),
                     float(ANUMA_BARS["separability"]), "Otsu separability",
                     v,
                     f"Otsu separability {sep:.2f} on {int(st['n_used'])} "
                     f"ungated nights, threshold {st['threshold_mag']:.2f} "
                     f"mag, span {st['span_mag']:.2f} mag; CV-S9 verdict "
                     f"{st['verdict']}",
                     "Classifying one night is easy -- the per-point "
                     "precision beats the state separation many times over. "
                     "What has to be shown is that the nightly magnitudes "
                     "separate into two populations at all.",
                     "" if v == "SUPPORTED" else
                     "more nights spread across states, not more points "
                     "per night"))
        # ---- 4. absolute duty cycle ------------------------------------
        n_used = int(st["n_used"]) if st else 0
        half = 100.0 * ch.duty_cycle_sigma(max(n_used, 1))
        v = fs.capability_verdict(half, ANUMA_BARS["duty_halfwidth_pp"],
                                  higher_is_better=False)
        n_need = int(math.ceil(0.25 / (ANUMA_BARS["duty_halfwidth_pp"]
                                       / 100.0) ** 2))
        rows.append((filt, "absolute duty cycle", 4,
                     _f(half) or float("nan"),
                     float(ANUMA_BARS["duty_halfwidth_pp"]),
                     "percentage points", v,
                     f"{n_used} independent nights give a binomial "
                     f"half-width of +/-{half:.0f} percentage points on a "
                     f"50% duty cycle (duty measured "
                     f"{100 * float(st['duty_with_limits'] or 0):.0f}%)",
                     "A duty cycle is a fraction of TIME, so its error bar "
                     "counts independent nights, not points. Nothing about "
                     "the photometry can shrink it.",
                     f"{max(n_need - n_used, 0)} more independent {filt} "
                     f"nights ({n_need} in total for +/-"
                     f"{ANUMA_BARS['duty_halfwidth_pp']:.0f} pp)"))
        # ---- 5. the colour goal, restated ------------------------------
        v = fs.capability_verdict(len(three),
                                  ANUMA_BARS["three_filter_nights"])
        rows.append((filt, "three-filter colour curves (Q5)", 5,
                     float(len(three)),
                     float(ANUMA_BARS["three_filter_nights"]), "nights", v,
                     f"{len(three)} nights carry a full orbit in ALL THREE "
                     f"filters ({', '.join(sorted(three)) or 'none'}) "
                     f"against the strategy's own bar of "
                     f"{ANUMA_BARS['three_filter_nights']}",
                     "This is a property of the SCHEDULE, not of any one "
                     "filter: it is repeated on every filter's row because "
                     "no single filter can fix it.",
                     f"{ANUMA_BARS['three_filter_nights'] - len(three)} more "
                     f"nights carrying all three filters over a full orbit"))
    con.executemany(
        "INSERT OR REPLACE INTO p4_anuma VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
    set_meta(con, {"anuma_period_d": period_d,
                   "anuma_three_filter_nights": ",".join(sorted(three)),
                   "anuma_seasons": seasons,
                   "anuma_bars": ";".join(f"{k}={v}"
                                          for k, v in ANUMA_BARS.items()),
                   "full_orbit_min_points": FULL_ORBIT_MIN_POINTS})
    stamp(con, "anuma")
    con.commit()
    for f_ in sorted(per_filter_nights):
        print(f"    {f_}: {len(per_filter_nights[f_])} full-orbit nights")
    print(f"    three-filter full-orbit nights: {len(three)} "
          f"({', '.join(sorted(three))})")
    con.close()


def _season(night: str) -> str:
    """Observing season label for a YYYY-MM-DD night.

    A season runs July -> June, because these are northern winter/spring
    targets and a calendar year would cut every AN UMa season in half.
    """
    y, m = int(night[:4]), int(night[5:7])
    return f"{y if m >= 7 else y - 1}/{(y + 1) if m >= 7 else y}"


# ===========================================================================
# STAGE verdict — the headline calls
# ===========================================================================
def cmd_verdict(args) -> None:
    """Write the closing decisions, each with its deciding number.

    Same shape as CV-S5's ``ch_verdict``: question, verdict, THE number that
    decides it, the reasoning, and what the alternative would have been.
    Every number here is a query, so a re-run that moved a measurement moves
    the verdict text with it instead of leaving prose behind.
    """
    con = connect(args.db)
    ensure_tables(con)
    con.execute("DELETE FROM p4_verdict")

    def one(sql, *p):
        r = con.execute(sql, p).fetchone()
        return r[0] if r and r[0] is not None else float("nan")

    def rng(col, where):
        lo = one(f"SELECT min({col}) FROM p4_run WHERE {where}")
        hi = one(f"SELECT max({col}) FROM p4_run WHERE {where}")
        return 1000 * lo, 1000 * hi

    # ---- YZ Cnc: the orbital hump ------------------------------------
    testable = "state='QUIESCENT' AND detection NOT IN ('AMPLITUDE ONLY')"
    n_test = int(one(f"SELECT count(*) FROM p4_run WHERE {testable}"))
    n_det = int(one("SELECT count(*) FROM p4_run WHERE state='QUIESCENT' "
                    "AND detection='DETECTED'"))
    a_lo, a_hi = rng("hump_amp", testable)
    f_lo, f_hi = rng("amp90_field", testable)
    s_lo, s_hi = rng("amp90_self", testable)
    # WHERE THE INSTRUMENTAL CONTOUR IS CLEARED, SCOPE BY SCOPE.  A blanket
    # "the photometry could see a hump this size" is false on any scope
    # whose fitted hump sits BELOW its own instrumental contour, and one of
    # ours does: it is uninformative rather than a non-detection, and the
    # deciding number has to say so or the table contradicts the figure
    # caption and Section 5.4, which already do.
    n_above = int(one(f"SELECT count(*) FROM p4_run WHERE {testable} AND "
                      "hump_amp > amp90_field"))
    below = con.execute(
        f"SELECT nights, filter FROM p4_run WHERE {testable} AND "
        "hump_amp <= amp90_field ORDER BY nights, filter").fetchall()
    below_txt = ("; on the remaining "
                 + (", ".join(f"{r['filter']} run of {r['nights']}"
                              for r in below))
                 + " it does not, so that scope is uninformative rather "
                   "than a non-detection") if below else ""
    # ---- YZ Cnc: flickering ------------------------------------------
    n_fl = int(one("SELECT count(*) FROM p4_flicker WHERE state='QUIESCENT' "
                   "AND detected=1"))
    n_flm = int(one("SELECT count(*) FROM p4_flicker WHERE "
                    "state='QUIESCENT' AND sf_excess IS NOT NULL"))
    fl_lo = 1000 * one("SELECT min(sf_excess) FROM p4_flicker WHERE "
                       "state='QUIESCENT' AND detected=1")
    fl_hi = 1000 * one("SELECT max(sf_excess) FROM p4_flicker WHERE "
                       "state='QUIESCENT' AND detected=1")
    fo_lo = 1000 * one("SELECT min(sf_floor) FROM p4_flicker WHERE "
                       "state='QUIESCENT' AND sf_excess IS NOT NULL")
    fo_hi = 1000 * one("SELECT max(sf_floor) FROM p4_flicker WHERE "
                       "state='QUIESCENT' AND sf_excess IS NOT NULL")
    # The timescales the MEASUREMENT populates, not the edges of the grid it
    # was searched on.  Quoting TAU_MIN_S--TAU_MAX_S here -- an earlier
    # revision did -- advertises the tau_edges_s constant as a measured
    # range, and it disagrees with every other timescale number in the
    # paper because the shortest and longest bins carry no detected excess.
    ta_lo = one("SELECT min(tau_s) FROM p4_flicker WHERE state='QUIESCENT' "
                "AND detected=1")
    ta_hi = one("SELECT max(tau_s) FROM p4_flicker WHERE state='QUIESCENT' "
                "AND detected=1")
    # The two eras' floors, per era rather than typed: the 8 s High Gain
    # frames and the 30 s Sloan-era frames are the comparison §4.19 asked
    # for, and both bounds move if the flicker stage is re-run.
    def _floor_era(prefix):
        lo = one("SELECT min(sf_floor) FROM p4_flicker WHERE "
                 "state='QUIESCENT' AND sf_floor IS NOT NULL AND "
                 "series_key LIKE ?", prefix)
        hi = one("SELECT max(sf_floor) FROM p4_flicker WHERE "
                 "state='QUIESCENT' AND sf_floor IS NOT NULL AND "
                 "series_key LIKE ?", prefix)
        return 1000 * lo, 1000 * hi
    hg_lo, hg_hi = _floor_era("yzcnc|e7|%")
    sl_lo, sl_hi = _floor_era("yzcnc|e72|%")
    # ---- YZ Cnc: outbursts -------------------------------------------
    b_lo = 1000 * one("SELECT min(amp90_blind) FROM p4_outburst WHERE "
                      "amp90_blind IS NOT NULL")
    b_hi = 1000 * one("SELECT max(amp90_blind) FROM p4_outburst WHERE "
                      "amp90_blind IS NOT NULL")
    n_ob = int(one("SELECT count(DISTINCT night) FROM p4_outburst"))
    ob_amp = one("SELECT max(amp_above_quiescence) FROM p4_outburst")
    rise = one("SELECT min(rate_mag_per_h) FROM p4_outburst WHERE "
               "rate_verdict='BRIGHTENING'")
    fade = one("SELECT max(rate_mag_per_h) FROM p4_outburst WHERE "
               "rate_verdict='FADING'")
    n_rate = int(one("SELECT count(*) FROM p4_outburst WHERE rate_verdict "
                     "IN ('BRIGHTENING','FADING')"))
    n_all = int(one("SELECT count(*) FROM p4_outburst"))
    # The blind contour is not measurable on every run-filter: a run too
    # short to place a periodogram maximum carries no contour at all.  The
    # verdict states its coverage, because a range quoted without it reads
    # as a property of all eighteen.
    n_blind = int(one("SELECT count(*) FROM p4_outburst WHERE "
                      "amp90_blind IS NOT NULL"))
    rows = [
        ("YZ-hump", 1, "CV-P3-yzcnc-superhump",
         "Is the quiescent orbital hump detected, and at what amplitude?",
         ("DETECTED" if n_det else
          "NOT DETECTED — PUBLISHED AS AN UPPER LIMIT"),
         f"{n_det} of {n_test} testable quiescent scopes clear both the 90% "
         f"recovery contour and the 1% false-alarm threshold. The folded "
         f"fundamental has semi-amplitude {a_lo:.0f}-{a_hi:.0f} mmag; the "
         f"INSTRUMENTAL contour is {f_lo:.0f}-{f_hi:.0f} mmag and the "
         f"fitted hump exceeds it on {n_above} of the {n_test} testable "
         f"scopes{below_txt}; but the contour set by "
         f"the star's OWN flickering is {s_lo:.0f}-{s_hi:.0f} mmag, so the "
         f"data cannot tell a coherent hump from the flickering",
         "Two nulls, and they disagree, which is the whole result. Against "
         "magnitude-matched field stars the fitted amplitude is comfortably "
         "significant. Against this star's own residuals -- rolled night by "
         "night, so every realization is made of light the telescope "
         "actually recorded -- it is not. YZ Cnc flickers at about ten "
         "times the photometric floor on exactly the timescales an orbital "
         "hump occupies, so the field-star contour alone would have "
         "produced a confident detection of nothing in particular.",
         f"The publishable statement is an upper limit: any COHERENT "
         f"orbital modulation with semi-amplitude above about "
         f"{s_lo:.0f} mmag would have been recovered in nine injections out "
         f"of ten, and none is. Phase zero is arbitrary in any case -- VSX "
         f"publishes no epoch for YZ Cnc -- so only hump SHAPES and "
         f"within-night inter-filter phase offsets were ever available."),
        ("YZ-flicker", 2, "CV-P3-yzcnc-superhump",
         "Is flickering separated from photometric noise, and how?",
         ("MEASURED" if n_fl else "NOT MEASURED"),
         f"{n_fl} of {n_flm} measurable timescale bins on the quiescent "
         f"runs exceed the measured field-star floor by "
         f"{fs.FLICKER_SIGMA_BAR:g} sigma or more on the variance excess. "
         f"Flickering amplitude {fl_lo:.0f}-{fl_hi:.0f} mmag over "
         f"{ta_lo:.0f}-{ta_hi:.0f} s; the floor subtracted is "
         f"{fo_lo:.0f}-{fo_hi:.0f} mmag",
         "The floor is not modelled and not a formal error bar. It is the "
         "SAME structure function computed on field stars within "
         f"{fs.FLOOR_MATCH_HALF_WIDTH} mag of the target through the same "
         "frames, so it carries identical photon noise, sky, ensemble "
         "zero-point wander and atmosphere, and no flickering. YZ Cnc's "
         "four held-out check stars sit about a magnitude BRIGHTER than the "
         "star at quiescence and would have given an optimistic floor.",
         f"The subtraction is done in variance, and bins where the target "
         f"does not exceed the floor are reported as NOT MEASURED rather "
         f"than as a small amplitude. The two eras differ by design: the "
         f"30 s Sloan-era frames reach a {sl_lo:.0f}-{sl_hi:.0f} mmag "
         f"floor, the 8 s High Gain frames only {hg_lo:.0f}-{hg_hi:.0f} "
         f"mmag, which is precisely what §4.19 feared."),
        ("YZ-superhump", 3, "CV-P3-yzcnc-superhump",
         "Can this season measure a superhump period or dP_sh/dt?",
         "NO — AND THAT IS A MEASUREMENT, NOT AN ABSENCE",
         f"the blind-search 90% recovery contour on the outburst dense runs "
         f"is {b_lo:.0f}-{b_hi:.0f} mmag on the {n_blind} of {n_all} "
         f"run-filters long enough to carry one, against superhump "
         f"semi-amplitudes of "
         f"{1000 * fs.SUPERHUMP_SEMI_AMP_FLOOR:.0f} mmag and up; and "
         f"CV-S7 already established that none of the {n_ob} outburst dense "
         f"runs sits inside a superoutburst (brightest {ob_amp:.2f} mag "
         f"above quiescence, against the ~3 mag a superoutburst reaches)",
         "A superhump period is a BLIND period determination -- P_sh is not "
         "known in advance the way P_orb is -- and a single night spanning "
         "about two orbits cannot place a periodogram maximum inside a 1% "
         "frequency window at any amplitude these data reach. Two "
         "independent reasons, either sufficient on its own.",
         "No superhump period, no dP_sh/dt, no O-C of superhump maxima, and "
         "no such promise in the abstract. Figure 11 of the strategy's "
         "figure list becomes the orbital-hump and flickering panel, which "
         "is the alternative the strategy itself names."),
        ("YZ-outburst", 4, "CV-P3-yzcnc-superhump",
         "What do the six normal-outburst dense runs support?",
         "A SEPARATE RESULT, NOT A CONSOLATION",
         f"{n_ob} dense runs inside normal outbursts, three filters at 8 s "
         f"exposures, peaking {ob_amp:.2f} mag above quiescence; "
         f"{n_rate} of {n_all} run-filters show a rate significant at "
         f"3 sigma, from {rise:+.3f} mag/h on the rise to {fade:+.3f} "
         f"mag/h on the decline",
         "Cycle-resolved three-colour coverage of normal-outburst structure "
         "is sampling the sparse survey record cannot provide. The "
         "2024-02-21 run catches a rise in all three filters at once, and "
         "the following nights catch the top and the decline.",
         "Reported as amplitude, covered duration and rate of change per "
         "run. No outburst recurrence time and no disc-instability model "
         "fit: nine dense nights of one season cannot support either."),
    ]
    # ---- AN UMa -------------------------------------------------------
    sup: dict[str, list[tuple[str, str]]] = {}
    for r in con.execute("SELECT filter, capability, verdict FROM p4_anuma "
                         "ORDER BY rank"):
        sup.setdefault(r["filter"], []).append((r["capability"],
                                                r["verdict"]))
    lines = []
    for f_ in sorted(sup):
        ok = [c for c, v in sup[f_] if v == "SUPPORTED"]
        lines.append(f"{f_} — " + ("; ".join(ok) if ok else "nothing"))
    n_sup = sum(1 for f_ in sup for _c, v in sup[f_] if v == "SUPPORTED")
    n_tot = sum(len(v) for v in sup.values())
    three = one("SELECT value FROM p4_meta WHERE "
                "key='anuma_three_filter_nights'")
    n_three = len([x for x in str(three).split(",") if x]) \
        if isinstance(three, str) else 0
    an_sig_lo = one("SELECT min(sigma_t_s) FROM p3_edge WHERE "
                    "target_key='anuma'")
    an_sig_hi = one("SELECT max(sigma_t_s) FROM p3_edge WHERE "
                    "target_key='anuma'")
    an_ok = int(one("SELECT sum(accepted) FROM p3_edge WHERE "
                    "target_key='anuma'"))
    an_fit = int(one("SELECT count(*) FROM p3_edge WHERE "
                     "target_key='anuma'"))
    an_amp = {r[0]: 1000.0 * (r[1] or 0.0) for r in con.execute(
        "SELECT filter, amplitude_mag FROM p3_period WHERE "
        "series_key LIKE 'anuma|%'")}
    rows.append((
        "ANUMA-role", 5, "CV-P4-anuma",
        "What role should AN UMa have in the paper?",
        "REDUCED-SCOPE TARGET",
        f"{n_sup} of {n_tot} graded capabilities are SUPPORTED across the "
        f"three filters. Per filter: " + " | ".join(lines),
        "AN UMa is neither a colour target nor a timing target, and the two "
        "failures have different causes. The colour goal fails on "
        f"SCHEDULING: {n_three} nights carry a full orbit in all three "
        f"filters against a bar of {ANUMA_BARS['three_filter_nights']}. The "
        f"timing goal fails on the EDGE: only {an_ok} of {an_fit} fitted "
        f"bright-phase "
        "edges survive a step signal-to-noise of 5, and CV-S9's cycle-count "
        "stage refuses an O-C on them. It is NOT established that the "
        "precision would otherwise be adequate -- no injection test was run "
        f"for this target, and the measured per-epoch errors of its own "
        f"fitted edges are {an_sig_lo:.0f}-{an_sig_hi:.0f} s, every one "
        f"outside the 60 s bar. What survives is real and costs "
        f"nothing extra: folded morphology in g and r, where the modulation "
        f"is detected at {an_amp.get('g', float('nan')):.0f} and "
        f"{an_amp.get('r', float('nan')):.0f} mmag, and a relative state "
        f"history that is bimodal in all three filters.",
        "Cutting AN UMa entirely would discard a resolved bimodal state "
        "history and a detected orbital modulation already in hand; "
        "promoting it to a full target would require three-filter "
        "full-orbit nights the schedule never delivered. The i band is the "
        "one filter with nothing of its own to offer: the modulation is not "
        "detected there at all (LS power 0.048, PDM theta 0.410, fitted "
        "semi-amplitude 57 mmag), so i appears only inside the state "
        "history."))
    con.executemany(
        "INSERT OR REPLACE INTO p4_verdict VALUES (?,?,?,?,?,?,?,?)", rows)
    stamp(con, "verdict")
    con.commit()
    for r in rows:
        print(f"    [{r[0]:14s}] {r[4]}")
    con.close()


# ===========================================================================
# report / status
# ===========================================================================
def cmd_report(args) -> None:
    from macro_phot.report_final import render_report
    path = render_report(args.db)
    print(f"  wrote {path}")


def cmd_status(args) -> None:
    con = connect(args.db, read_only=True)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table' "
        "AND name LIKE 'p4_%'")}
    print("CV-S10 status")
    for t in ("p4_run", "p4_flicker", "p4_outburst", "p4_gate", "p4_anuma",
              "p4_verdict"):
        n = (con.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
             if t in have else 0)
        print(f"  {t:14s} {n:6d} rows")
    if "p4_run" in have:
        print("\n  YZ Cnc dense runs")
        for r in con.execute(
                "SELECT scope, state, n_points, cycles, hump_amp, "
                "amp90_field, amp90_self, detection FROM p4_run "
                "ORDER BY scope"):
            def _m(v):
                return f"{1000 * v:6.1f}" if v is not None else "   -  "
            print(f"    {r['scope'].replace('yzcnc|', ''):40s} "
                  f"[{r['state'][:9]:9s}] n={r['n_points']:3d} "
                  f"cyc={r['cycles']:4.2f}  A={_m(r['hump_amp'])} "
                  f"A90f={_m(r['amp90_field'])} A90s={_m(r['amp90_self'])} "
                  f"mmag  {r['detection']}")
    if "p4_anuma" in have:
        print("\n  AN UMa, per filter")
        for r in con.execute("SELECT filter, capability, verdict FROM "
                             "p4_anuma ORDER BY filter, rank"):
            print(f"    {r['filter']:>2s}  {r['capability']:52s} "
                  f"{r['verdict']}")
    con.close()


# ===========================================================================
# CLI
# ===========================================================================
def _common(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Every option on every parser, through a shared parent, so that both
    ``run_cv_final.py --db X hump`` and ``run_cv_final.py hump --db X``
    work.  Argparse accepts only the first by default, and a stage that ran
    against the wrong database because the flag landed after the subcommand
    is exactly the kind of thing that makes a published number wrong."""
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--char-db", type=Path, default=CHAR_DB)
    return parser


def main() -> None:
    p = _common(argparse.ArgumentParser(description=__doc__.split("\n")[0]))
    parent = _common(argparse.ArgumentParser(add_help=False))
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("hump", "flicker", "outburst", "gate", "anuma", "verdict",
                 "report", "status", "all"):
        sub.add_parser(name, parents=[parent])
    args = p.parse_args()
    table = {"hump": cmd_hump, "flicker": cmd_flicker,
             "outburst": cmd_outburst, "gate": cmd_gate, "anuma": cmd_anuma,
             "verdict": cmd_verdict, "report": cmd_report,
             "status": cmd_status}
    order = ("hump", "flicker", "outburst", "gate", "anuma", "verdict",
             "report")
    if args.cmd == "all":
        for name in order:
            print(f"\n=== {name} ===", flush=True)
            table[name](args)
        record_stage("CV-S10")
        cmd_status(args)
    else:
        table[args.cmd](args)
        if args.cmd in order and args.cmd != "report":
            record_stage("CV-S10")


if __name__ == "__main__":
    main()
