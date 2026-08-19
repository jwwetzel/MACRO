#!/usr/bin/env python
"""Tie the CV ensembles to a photometric catalogue: relative -> publishable.

WHAT THIS SCRIPT DOES, AND WHY IT IS THE TOP-RANKED ACTION
----------------------------------------------------------
``run_cv_photometry.py`` solved a Honeycutt ensemble for every (target, era,
filter) series and produced 3.18 million light-curve rows.  Every one of
them is DIFFERENTIAL: the ensemble's gauge freedom was fixed by demanding
``mean(ZP) = 0``, which makes the magnitudes internally consistent and
externally meaningless.  The CV characterization graded the strategy's
calibration goal NOT SUPPORTED on that fact -- "7 of 14 (target, era) blocks
carry a catalogue magnitude tie; the other 7 were tied by WCS with zero
catalogue matches" -- and named the repair the highest-value single action
on the whole list.

This is that repair.  For every solved series it fits, on COMPARISON STARS
ONLY,

    m_ensemble - m_catalogue  =  ZP0 + k * (colour - colour_ref)

robustly, with uncertainties on both terms, a stated colour range of
validity, and an achieved accuracy measured on stars the fit never saw.
The published magnitude is ``m_nat = m_ensemble - ZP0``: the telescope's own
natural system, zero-pointed onto the catalogue.  ``k`` is published as
metadata and is NEVER applied to a science target -- CVs are blue, variable
and routinely outside the colour range any transformation was calibrated
over, so transforming them would trade a known-size bandpass error for an
unknown-size extrapolation error.  All of that arithmetic lives in
``macro_phot.cattie`` and is unit-tested; this script is I/O, staging and
bookkeeping.

CATALOGUES
----------
PRIMARY   ATLAS-REFCAT2 (Tonry et al. 2018), via VizieR TAP.  The system the
          strategy itself names, m ~ 19 deep, and it ships its own blend
          metrology (the R1/R10 contamination radii).
SECONDARY Gaia DR3 standardised synthetic photometry
          (``gaiadr3.synthetic_photometry_gspc``), via the ESA archive.
          Shallower (it needs a BP/RP spectrum, so it stops near G = 17.65)
          but it publishes Sloan, PS1-y AND Johnson-Cousins bands for the
          SAME stars from the SAME spectra, which is the only way to
          (a) cover the y filter at all and (b) test whether the uppercase
          'G'/'R'/'I' filter labels of eras 6/7 mean the same glass as the
          lowercase 'g'/'r'/'i' of eras 72/76.
Both are fetched for every field, both are tied, and the DIFFERENCE between
them on the same stars is reported as the systematic floor of this
calibration -- which is a number the error budget needs and which no single
catalogue can supply.

STAGES (each resumable, each safe to repeat)
--------------------------------------------
    fieldfix   re-run the (target, era) field ties that failed, with the
               retry/backoff already built into macro_phot.gaia
    fetch      pull + cache catalogue photometry per FIELD, recording the
               query text and the pull date; never re-fetches silently
    match      match catalogue sources to each block's reference stars,
               with blend metrology and proper-motion propagation
    solve      the robust ZP + colour-term fit, per (series, catalogue, band)
    validate   independent check stars, cross-catalogue difference,
               residual-vs-magnitude and residual-vs-position regressions
    apply      write the calibrated natural-system column on the light curve
    report     docs/CV_TimeSeries/cv_catalogue_tie.html + figures
    status     progress + verdict summary (safe at any time)
    all        fieldfix -> ... -> report, in order

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_cattie.py fetch
    $P pipeline/scripts/run_cv_cattie.py match
    $P pipeline/scripts/run_cv_cattie.py solve
    $P pipeline/scripts/run_cv_cattie.py validate
    $P pipeline/scripts/run_cv_cattie.py apply
    $P pipeline/scripts/run_cv_cattie.py report
    $P pipeline/scripts/run_cv_cattie.py status

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite)
--------------------------------------------------------------
``cv_cat_fetch``     one row per (catalogue, field): the ADQL, the pull
                     date, the row count, the cache file and its sha256.
``cv_cat_match``     one row per (catalogue, target, era, reference star):
                     the matched catalogue source, its separation, its
                     photometry and its blend metrology.
``cv_cattie``        THE DELIVERABLE.  One row per (series, catalogue,
                     band hypothesis): n_tie_stars, ZP +/- err, colour term
                     +/- err, residual RMS, colour range, achieved accuracy
                     on independent checks, and a verdict.
``cv_cattie_star``   per tie star: its role in the fit and its residual
                     (what every figure is drawn from).
``cv_cattie_veto``   the rejection census -- why each candidate was dropped.
``cv_cattie_trend``  residual regressed on magnitude, on x, on y.
``cv_cattie_cross``  ATLAS-REFCAT2 minus Gaia synthetic, on common stars.
``cv_cat_meta``      build stamps and every constant.
and ONE new column on ``cv_lightcurve``: ``cal_mag``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import socket
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import cattie as ct                            # noqa: E402
from macro_phot import gaia as gg                              # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_CACHE = REPO_ROOT / "products" / "phot" / "catalogue_cache"

#: This stage's code version, recorded in cv_cat_meta and used by the
#: provenance DAG to decide whether a downstream page is stale.
#: v1.1 answers adversarial review: the target colour is formed from paired
#: epochs instead of campaign means, the blend veto is taken in the tie band
#: on the WORST contaminant in the aperture instead of in Gaia G on the
#: nearest one, each block's rigid astrometric offset from the catalogue is
#: measured and (where coherent) removed, and the per-star table is cleared
#: before it is rewritten.  Every one of those changes the tie stars or the
#: published colour position, so the version moves.
CATTIE_CODE_VERSION = "CV-S6 catalogue tie v1.1 (2026-08-19)"

#: The manifest and this product are shared with a running rclone transfer
#: and possibly other stages: five minutes of patience beats a spurious
#: 'database is locked' in the middle of a network-bound campaign.
BUSY_TIMEOUT_MS = 300_000

#: Archive endpoints.  VizieR's TAP service rather than its HTTP
#: query_region interface: the latter aborted the connection on a 0.55-deg
#: cone during this build, and TAP answered the same cone in 2 seconds.
VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
REFCAT2_TABLE = '"J/ApJ/867/105/refcat2"'

#: Retry policy for every archive call.  An HTTP 500 from a shared public
#: archive is transient, not a verdict -- this campaign lost two EU UMa
#: field ties to exactly that and they came back on a re-run.  Attempts
#: with a widening pause; only a persistent failure is recorded as one.
#: Six attempts with a widening pause is up to about four minutes of
#: patience per query.  Measured during this build: the ESA archive
#: answered a cone in 30 s one hour and returned HTTP 500 to five
#: consecutive attempts at the same cone the next, then recovered.
NET_RETRIES = 6
NET_PAUSE_S = 15.0

#: Socket timeout, seconds.  astroquery's TAP client sets none, so a
#: server that accepts a connection and then stops answering hangs the
#: whole stage FOREVER -- observed here while the ESA archive was
#: degraded.  A timeout converts that into a retry, which is the only
#: behaviour that lets a resumable stage stay resumable.  Generous,
#: because a healthy-but-loaded archive genuinely took 84 seconds to
#: answer `SELECT TOP 1 source_id ... WHERE source_id = <one id>`.
SOCKET_TIMEOUT_S = 420.0
socket.setdefaulttimeout(SOCKET_TIMEOUT_S)

#: Gaia's synchronous TAP caps every query at 2,000 rows no matter what, so
#: a dense cone is walked in brightness slices (see gaia_cone_sliced).
GAIA_PAGE = 2000

#: Gaia synthetic photometry needs a BP/RP spectrum, which in DR3 means
#: G <~ 17.65.  Asking for fainter stars costs archive time and returns
#: nothing, so the Gaia pull stops here.
GAIA_G_LIMIT = 17.8

#: Reference epochs of the two catalogues (Julian years).  Positions are
#: propagated from these to each block's own observation epoch before
#: matching: 8-10 years at 30 mas/yr is 0.3 arcsec, a quarter of the match
#: tolerance, and high-proper-motion stars are exactly the ones a
#: nearest-neighbour match gets wrong.
CAT_EPOCH_JYR = {"refcat2": 2015.5, "gaia_gspc": 2016.0}

#: The Gaia synthetic bands this stage pulls.  Kept in one list because it
#: is used three times: to build the ADQL, to name the cache columns, and
#: to convert fluxes to magnitude errors.
GSPC_BANDS = ("g_sdss", "r_sdss", "i_sdss", "z_sdss", "y_ps1",
              "b_jkc", "v_jkc", "r_jkc", "i_jkc")


# ===========================================================================
# Small helpers
# ===========================================================================
def connect(db_path: Path) -> sqlite3.Connection:
    """Open the products database patiently; WAL, as every other stage."""
    con = sqlite3.connect(db_path, timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    con.execute("PRAGMA synchronous=NORMAL")
    return con


def git_commit() -> str:
    """HEAD's short hash with an honest '-dirty' marker."""
    try:
        head = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              cwd=REPO_ROOT, capture_output=True, text=True,
                              check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT,
                               capture_output=True, text=True,
                               check=True).stdout.strip()
        return f"{head}-dirty" if dirty else head
    except Exception:
        return ""


def fnum(x):
    """float(x) when finite, else None -- the SQL-safe cast used everywhere."""
    if x is None:
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def meta_write(con: sqlite3.Connection, extra: dict = ()) -> None:
    """Upsert build metadata.  Every constant the report quotes is here, so
    a page can never disagree with the code that produced it."""
    con.execute("""CREATE TABLE IF NOT EXISTS cv_cat_meta
                   (key TEXT PRIMARY KEY, value TEXT)""")
    rows = {
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "cattie_code_version": CATTIE_CODE_VERSION,
        "git_commit": git_commit(),
        "primary_catalogue": ct.PRIMARY_CATALOGUE,
        "cone_radius_deg": str(ct.CONE_RADIUS_DEG),
        "cat_mag_max": str(ct.CAT_MAG_MAX),
        "match_tol_arcsec": str(ct.MATCH_TOL_ARCSEC),
        "ambiguity_factor": str(ct.AMBIGUITY_FACTOR),
        "blend_aperture_arcsec": str(ct.BLEND_APERTURE_ARCSEC),
        "blend_dmag": str(ct.BLEND_DMAG),
        "blend_annulus_arcsec": str(ct.BLEND_ANNULUS_ARCSEC),
        "sat_frac_max": str(ct.SAT_FRAC_MAX),
        "near_veto_frac": str(ct.NEAR_VETO_FRAC),
        "min_tie_stars": str(ct.MIN_TIE_STARS),
        "min_check_stars": str(ct.MIN_CHECK_STARS),
        "holdout_fraction": str(ct.HOLDOUT_FRACTION),
        "holdout_salt": ct.HOLDOUT_SALT,
        "huber_delta": str(ct.HUBER_DELTA),
        "clip_sigma": str(ct.CLIP_SIGMA),
        "accuracy_goal_mag": str(ct.ACCURACY_GOAL_MAG),
        "accuracy_stretch_mag": str(ct.ACCURACY_STRETCH_MAG),
        "trend_sigma": str(ct.TREND_SIGMA),
        "vizier_tap_url": VIZIER_TAP_URL,
        "refcat2_table": REFCAT2_TABLE,
        "gaia_g_limit": str(GAIA_G_LIMIT),
        "net_retries": str(NET_RETRIES),
        **dict(extra),
    }
    con.executemany("INSERT OR REPLACE INTO cv_cat_meta VALUES (?,?)",
                    list(rows.items()))
    con.commit()


#: Substrings that mark an archive error as OURS, not the archive's.  A
#: malformed query fails identically on every attempt, so retrying it four
#: times with a widening pause wastes a minute and, worse, files a code bug
#: under "the archive was down" -- which is exactly the misdiagnosis this
#: stage exists to correct for EU UMa era 80.
PERMANENT_ERROR_MARKS = ("Error 400", "Incorrect ADQL", "Unknown column",
                         "Ambiguous column", "unresolved identifier")


def with_retry(fn, what: str, retries: int | None = None):
    """Run a network call with widening backoff; raise only on persistence.

    Existence proof that this matters: the two EU UMa field ties this stage
    repairs were lost to an archive HTTP 500, and the identical query
    succeeded on a later attempt.  A build that gives up on the first 500
    loses real data to somebody else's load spike.

    The converse matters just as much.  A 400 is a statement about OUR
    query, not about the archive's health, and retrying it dresses a
    syntax error up as an outage.  Those are raised immediately.
    """
    last = None
    n = NET_RETRIES if retries is None else max(1, int(retries))
    for attempt in range(n):
        try:
            return fn()
        except Exception as e:                      # HTTP 500, timeout, reset
            last = e
            if any(mark in str(e) for mark in PERMANENT_ERROR_MARKS):
                print(f"    {what}: {type(e).__name__} is a QUERY error, "
                      f"not an outage — not retrying", flush=True)
                raise
            wait = NET_PAUSE_S * (attempt + 1)
            print(f"    {what}: {type(e).__name__} "
                  f"(attempt {attempt + 1}/{n})"
                  + (f"; retrying in {wait:.0f}s" if attempt < n - 1
                     else ""), flush=True)
            if attempt < n - 1:
                time.sleep(wait)
    raise last


def atomic_write_gz(path: Path, payload: dict) -> str:
    """Write a gzipped-JSON cache atomically; return its sha256.

    Atomic because a half-written cache that LOOKS present is worse than an
    absent one: the next run would skip the fetch and read garbage.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    raw = json.dumps(payload, separators=(",", ":")).encode()
    with gzip.open(tmp, "wb") as fh:
        fh.write(raw)
    tmp.replace(path)
    return hashlib.sha256(raw).hexdigest()


def read_cache(path: Path) -> dict | None:
    """Read a gzipped-JSON cache, or None if it is absent or unreadable."""
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rb") as fh:
            return json.loads(fh.read().decode())
    except Exception:
        return None


def col_list(tab, name: str, dtype=float) -> list:
    """One astropy-table column as a plain Python list, masks -> NaN.

    Necessary, not decorative: a VOTable column with missing values comes
    back MASKED, and ``np.asarray(masked_column)`` silently substitutes the
    FILL VALUE (often 1e20, sometimes 0) for every missing entry.  A zero
    magnitude that looks like a measurement is exactly the kind of number
    that survives every downstream check and ruins a zero point.
    """
    col = tab[name]
    if dtype is int:
        return [int(v) for v in np.asarray(col)]
    try:
        arr = np.asarray(col.filled(np.nan), dtype=float)
    except (AttributeError, TypeError, ValueError):
        arr = np.asarray(col, dtype=float)
    return arr.tolist()


def epoch_jyr(bjd_tdb: float) -> float:
    """Julian year of a BJD_TDB timestamp (for proper-motion propagation)."""
    return 2000.0 + (float(bjd_tdb) - 2451545.0) / 365.25


def propagate(ra, dec, pmra, pmdec, from_jyr: float, to_jyr: float):
    """Move catalogue positions to the observation epoch.

    ``pmra`` is mu_alpha* (already including cos(dec)), which is how both
    catalogues publish it, so the RA increment is pm/cos(dec).  Missing
    proper motions are treated as zero: unknown motion is better modelled as
    none than as an invented value, and a star with no PM in either
    catalogue is usually distant enough for that to be true.
    """
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    pmra = np.nan_to_num(np.asarray(pmra, dtype=float))
    pmdec = np.nan_to_num(np.asarray(pmdec, dtype=float))
    dt = float(to_jyr) - float(from_jyr)
    dd = pmdec * dt / 3.6e6                     # mas/yr * yr -> deg
    dr = pmra * dt / 3.6e6 / np.cos(np.radians(dec))
    return ra + dr, dec + dd


# ===========================================================================
# Stage: fieldfix -- repair the (target, era) ties that failed
# ===========================================================================
def cmd_fieldfix(args) -> None:
    """Re-run the field-tie stage for every (target, era) not already 'ok'.

    Delegated to ``run_cv_photometry.py field`` rather than reimplemented:
    ``cv_field_tie`` and ``cv_ref_stars.ra_deg`` belong to THAT stage, and
    two stages writing one table is how a product acquires two
    incompatible definitions of the same column.  That stage is already
    resumable -- it skips ties whose status is 'ok' -- so this is a
    one-line delegation plus an honest before/after report.

    Why it belongs here at all: a block with no sky positions for its
    reference stars cannot be matched to any catalogue, so the two EU UMa
    eras that failed on archive HTTP 500s are the difference between 208
    measured frames being calibratable and not.
    """
    con = connect(args.db)
    before = {(r[0], r[1]): r[2] for r in con.execute(
        "SELECT target_key, era_id, status FROM cv_field_tie")}
    bad = {k: v for k, v in before.items() if v != "ok"}
    print(f"field ties not 'ok' before: {len(bad)}")
    for (tk, era), st in sorted(bad.items()):
        print(f"    {tk} era{era}: {st}")
    con.close()
    if not bad and not args.force:
        print("  nothing to repair")
        return
    cmd = [sys.executable,
           str(PIPELINE_ROOT / "scripts" / "run_cv_photometry.py"),
           "field"]
    if args.force:
        cmd.append("--force")
    print("  running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=REPO_ROOT, check=False)
    con = connect(args.db)
    after = {(r[0], r[1]): r[2] for r in con.execute(
        "SELECT target_key, era_id, status FROM cv_field_tie")}
    for k in sorted(bad):
        was, now = bad[k], after.get(k, "(gone)")
        verdict = "REPAIRED" if now == "ok" else "still failing"
        print(f"    {k[0]} era{k[1]}: {verdict}\n"
              f"        was: {was}\n        now: {now}")
    meta_write(con, {"stage_fieldfix": datetime.now(timezone.utc).isoformat()})
    con.close()


# ===========================================================================
# Stage: fetch -- catalogue photometry per FIELD, cached
# ===========================================================================
def _catalogues(args) -> tuple[str, ...]:
    """Which catalogues this invocation touches.

    ``--catalogue`` exists for one operational reason and it is worth
    stating: the two archives have wildly different service levels (VizieR
    answered a cone in 2 seconds while the ESA archive was returning HTTP
    500s on the same afternoon), and a stage that could only ever run both
    would hold the primary tie hostage to the cross-check.  The DEFAULT is
    still both, because a run that quietly did half the work would be worse
    than a slow one.
    """
    if getattr(args, "catalogue", None):
        return (args.catalogue,)
    return ct.CATALOGUES


def _field_centres(con) -> list[tuple[str, float, float]]:
    """One cone centre per TARGET (not per era): the fields are the same sky.

    Uses the catalogue coordinates already resolved and recorded by the
    field-tie stage, so the cone is centred on the same position the target
    identification used and no second SIMBAD call can disagree with it.
    """
    rows = con.execute("""SELECT target_key, avg(target_ra), avg(target_dec)
                          FROM cv_field_tie
                          WHERE target_ra IS NOT NULL
                          GROUP BY target_key ORDER BY target_key""").fetchall()
    return [(r[0], float(r[1]), float(r[2])) for r in rows]


#: Page size for the VizieR TAP cone.  Its synchronous service caps a
#: result at 2,000 rows and IGNORES a larger ``MAXREC`` -- measured here:
#: the VV Pup cone holds 25,136 sources and came back as exactly 2,000,
#: with no error and no warning.  A silent truncation is the most dangerous
#: kind of archive failure, because the result LOOKS like an answer, so the
#: cone is paged in brightness slices and the page size is checked.
VIZIER_PAGE = 2000


def refcat2_adql(ra: float, dec: float, radius: float, mag_max: float,
                 mag_min: float = -5.0, page: int = VIZIER_PAGE) -> str:
    """One page of the ATLAS-REFCAT2 cone, verbatim as sent and as stored.

    The column names are quoted case-sensitively on purpose: VizieR's copy
    of this catalogue carries BOTH ``Gmag`` (Gaia G) and ``gmag`` (PS1 g),
    and an unquoted reference to either is rejected as ambiguous.  A build
    that quietly took whichever the server preferred would tie the g band
    to Gaia G and never say so.
    """
    return (
        f'SELECT TOP {page} c.RA_ICRS AS ra, c.DE_ICRS AS dec, '
        'c.pmRA AS pmra, '
        'c.pmDE AS pmdec, c."Gmag" AS gaia_gmag, '
        'c."gmag" AS gmag, c."e_gmag" AS e_gmag, '
        'c."rmag" AS rmag, c."e_rmag" AS e_rmag, '
        'c."imag" AS imag, c."e_imag" AS e_imag, '
        'c."zmag" AS zmag, c."e_zmag" AS e_zmag, '
        'c.dupvar AS dupvar, c.RP1 AS rp1, c.R1 AS r1, c.R10 AS r10 '
        f'FROM {REFCAT2_TABLE} AS c '
        "WHERE 1=CONTAINS(POINT('ICRS',c.RA_ICRS,c.DE_ICRS),"
        f"CIRCLE('ICRS',{ra:.6f},{dec:.6f},{radius:.4f})) "
        f'AND c."rmag" >= {mag_min:.6f} AND c."rmag" < {mag_max:.2f} '
        # ORDER BY names the ALIAS, not the qualified quoted column: this
        # ADQL parser accepts `c."rmag"` in SELECT and WHERE but rejects it
        # after ORDER BY.
        'ORDER BY rmag')


def fetch_refcat2(ra: float, dec: float, radius: float,
                  mag_max: float, retries: int | None = None
                  ) -> tuple[dict, str]:
    """Pull ATLAS-REFCAT2 over one cone, paged; return (columns, queries).

    Paged faintward in r: each page asks for the ``VIZIER_PAGE`` brightest
    sources fainter than the last page's faint end.  The boundary uses
    ``>=`` and the result is de-duplicated on position, because a strict
    ``>`` would silently drop every source sharing the boundary magnitude
    -- the same trap as the Gaia cone, and the same fix.
    """
    from astroquery.utils.tap.core import TapPlus
    tap = TapPlus(url=VIZIER_TAP_URL)
    cols: dict[str, list] = {}
    queries: list[str] = []
    lo = -5.0
    for _page in range(60):
        adql = refcat2_adql(ra, dec, radius, mag_max, lo)
        queries.append(adql)
        tab = with_retry(lambda a=adql: tap.launch_job(a).get_results(),
                         "refcat2 cone page", retries)
        if len(tab) == 0:
            break
        for name in tab.colnames:
            cols.setdefault(name, []).extend(col_list(tab, name))
        if len(tab) < VIZIER_PAGE:
            break
        r = np.asarray(col_list(tab, "rmag"), dtype=float)
        lo = float(np.nanmax(r))
    if cols:
        # De-duplicate on exact position (RA and Dec together are unique in
        # this catalogue; the boundary overlap is the only source of dupes).
        key = list(zip(cols["ra"], cols["dec"]))
        seen: set = set()
        keep = [i for i, k in enumerate(key)
                if not (k in seen or seen.add(k))]
        cols = {n: [v[i] for i in keep] for n, v in cols.items()}
    return cols, "\n".join(queries)


#: Narrowest magnitude window the cone walk will bisect to.  Below this a
#: window that still fills the page is reported as a possible truncation
#: rather than split forever.
GAIA_MIN_WINDOW_MAG = 0.02


def gaia_cone_sliced(ra: float, dec: float, radius: float,
                     g_max: float, retries: int | None = None
                     ) -> tuple[dict, list[str]]:
    """Gaia DR3 positions + G + BP-RP over one cone, in magnitude windows.

    The archive's synchronous service caps EVERY query at 2,000 rows and
    its asynchronous service answered 500 to every attempt during this
    build, so the cone has to be split.  It is split by BISECTION on
    magnitude rather than by an ``ORDER BY`` walk, for a reason that turned
    out to matter while the archive was degraded: ``ORDER BY
    phot_g_mean_mag`` makes the server sort the whole cone before it can
    return the first row, and a sort is exactly what a loaded server kills
    with "canceling statement due to statement timeout".  A bare magnitude
    window needs no sort at all.

    The split is self-correcting and provably complete: any window that
    comes back EXACTLY full is assumed truncated and is bisected, so no
    star can be lost by guessing the window sizes wrong.  A window that
    still fills at ``GAIA_MIN_WINDOW_MAG`` is reported, not silently
    accepted -- a silent truncation is the failure mode this whole routine
    exists to avoid.
    """
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1
    cols: dict[str, list] = {k: [] for k in
                             ("source_id", "ra", "dec", "pmra", "pmdec",
                              "phot_g_mean_mag", "bp_rp")}
    queries: list[str] = []
    todo = [(-5.0, float(g_max))]
    n_query = 0
    while todo and n_query < 60:
        lo, hi = todo.pop()
        adql = (f"SELECT TOP {GAIA_PAGE} source_id, ra, dec, pmra, pmdec, "
                "phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source "
                "WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
                f"CIRCLE('ICRS',{ra:.6f},{dec:.6f},{radius:.4f})) "
                f"AND phot_g_mean_mag >= {lo:.6f} "
                f"AND phot_g_mean_mag < {hi:.6f}")
        if not queries:
            queries.append(adql)
        n_query += 1
        tab = with_retry(lambda a=adql: Gaia.launch_job(a).get_results(),
                         f"gaia cone G in [{lo:.2f}, {hi:.2f})",
                         retries)
        if len(tab) >= GAIA_PAGE and (hi - lo) > GAIA_MIN_WINDOW_MAG:
            mid = 0.5 * (lo + hi)
            todo.extend([(lo, mid), (mid, hi)])
            continue
        if len(tab) >= GAIA_PAGE:
            print(f"      WARNING: G in [{lo:.3f}, {hi:.3f}) still fills "
                  f"the page at the minimum window width — this cone may "
                  f"be truncated", flush=True)
        for k in cols:
            cols[k].extend(col_list(tab, k,
                                    int if k == "source_id" else float))
        print(f"      G in [{lo:.2f}, {hi:.2f}): {len(tab)} "
              f"(total {len(cols['source_id'])})", flush=True)
    # De-duplicate on source_id, keeping first occurrence.
    sid = np.asarray(cols["source_id"], dtype=np.int64)
    if sid.size:
        _, keep = np.unique(sid, return_index=True)
        keep = np.sort(keep)
        for k in cols:
            cols[k] = [cols[k][i] for i in keep]
    return cols, queries


#: HEALPix ranges per query.  ONE.  A first version OR-ed a dozen ranges
#: into each query with an ORDER BY on top, and the archive answered
#: "canceling statement due to statement timeout" every time: an OR of
#: ranges plus a sort is a plan the optimiser will not drive from the
#: index, which is the entire point of using ranges.  A single BETWEEN with
#: no ORDER BY is one contiguous read of one partition.
GSPC_RANGES_PER_QUERY = 1


def fetch_gspc(source_ids: list[int], retries: int | None = None
               ) -> tuple[dict, list[str]]:
    """Standardised synthetic photometry for the sources of one cone.

    NOT queried by cone: a JOIN between ``gaia_source`` and
    ``synthetic_photometry_gspc`` under a CONTAINS predicate exceeded the
    archive's statement timeout on every attempt, synchronous AND
    asynchronous, during this build.  NOT queried by IN-list either, for
    the reason given in ``macro_phot.cattie.hpx_ranges``.  Queried as one
    ``source_id BETWEEN`` range per HEALPix run -- a contiguous read of the
    partition the archive is physically ordered by.

    Magnitudes come with their fluxes so a per-band error can be formed:
    sigma_m = 1.0857 * sigma_F / F.  The table publishes no magnitude
    errors of its own, and a magnitude with no error cannot be weighted.
    """
    from astroquery.gaia import Gaia
    Gaia.ROW_LIMIT = -1
    sel = ["source_id", "c_star"]
    for b in GSPC_BANDS:
        sel += [f"{b}_mag", f"{b}_flux", f"{b}_flux_error", f"{b}_flag"]
    cols: dict[str, list] = {k: [] for k in sel}
    queries: list[str] = []
    wanted = set(int(s) for s in source_ids)
    ranges = ct.hpx_ranges(source_ids)
    print(f"      {len(ranges)} source_id ranges", flush=True)
    for bi, (lo, hi) in enumerate(ranges):
        adql = (f"SELECT TOP {GAIA_PAGE} {', '.join(sel)} "
                "FROM gaiadr3.synthetic_photometry_gspc "
                f"WHERE source_id BETWEEN {lo} AND {hi}")
        if bi == 0:
            queries.append(adql)
        tab = with_retry(lambda a=adql: Gaia.launch_job(a).get_results(),
                         f"gspc range {bi + 1}/{len(ranges)}",
                         retries)
        if len(tab) == 0:
            continue
        sids = col_list(tab, "source_id", int)
        take = [i for i, sd in enumerate(sids) if sd in wanted]
        for k in sel:
            if k in tab.colnames:
                v = col_list(tab, k, int if k == "source_id" else float)
            else:
                v = [float("nan")] * len(tab)
            cols[k].extend([v[i] for i in take])
        if (bi + 1) % 5 == 0 or bi + 1 == len(ranges):
            print(f"      range {bi + 1}/{len(ranges)}: "
                  f"{len(cols['source_id'])} kept", flush=True)
    return cols, queries


def _gspc_to_frame(gaia_cols: dict, gspc_cols: dict) -> dict:
    """Join the cone (positions) to the synthetic photometry (magnitudes).

    Returns ONE column dict shaped exactly like the REFCAT2 one, so that
    every downstream step -- matching, blend metrology, the veto census --
    runs identical code on both catalogues.  That identity is what makes
    the cross-catalogue difference in ``cv_cattie_cross`` a measurement of
    the CATALOGUES rather than of two different pipelines.
    """
    sid_all = np.asarray(gaia_cols["source_id"], dtype=np.int64)
    out: dict[str, list] = {
        "source_id": sid_all.tolist(),
        "ra": list(gaia_cols["ra"]), "dec": list(gaia_cols["dec"]),
        "pmra": list(gaia_cols["pmra"]), "pmdec": list(gaia_cols["pmdec"]),
        "gaia_gmag": list(gaia_cols["phot_g_mean_mag"]),
        "bp_rp": list(gaia_cols["bp_rp"]),
    }
    n = len(sid_all)
    idx = {int(s): i for i, s in enumerate(
        np.asarray(gspc_cols.get("source_id", []), dtype=np.int64))}
    nan = float("nan")
    out["c_star"] = [nan] * n
    for b in GSPC_BANDS:
        out[f"{b}_mag"] = [nan] * n
        out[f"{b}_mag_error"] = [nan] * n
        out[f"{b}_flag"] = [1.0] * n            # 1 = "no synthetic value"
    for i, s in enumerate(sid_all):
        j = idx.get(int(s))
        if j is None:
            continue
        out["c_star"][i] = float(gspc_cols["c_star"][j])
        for b in GSPC_BANDS:
            m = float(gspc_cols[f"{b}_mag"][j])
            f = float(gspc_cols[f"{b}_flux"][j])
            fe = float(gspc_cols[f"{b}_flux_error"][j])
            out[f"{b}_mag"][i] = m
            out[f"{b}_mag_error"][i] = (1.0857 * fe / f
                                        if f and math.isfinite(f) and f > 0
                                        and math.isfinite(fe) else nan)
            # A MISSING flag is treated as a set flag.  "The catalogue did
            # not tell us whether this magnitude is inside its validity
            # range" is not evidence that it is.
            fl = gspc_cols[f"{b}_flag"][j]
            out[f"{b}_flag"][i] = (float(fl)
                                   if fl is not None and math.isfinite(
                                       float(fl)) else 1.0)
    return out


def cmd_fetch(args) -> None:
    """Pull and cache both catalogues for every field, once.

    RESUMABLE AND NEVER SILENT: a cache file that exists is used, and the
    ``cv_cat_fetch`` row records the exact query, the pull date and the
    cache's sha256, so 'where did this magnitude come from' has a written
    answer months later.  ``--force`` re-pulls.
    """
    con = connect(args.db)
    con.execute("""CREATE TABLE IF NOT EXISTS cv_cat_fetch (
        catalogue TEXT, field_key TEXT, ra_deg REAL, dec_deg REAL,
        radius_deg REAL, mag_limit REAL, service TEXT, query TEXT,
        pulled_utc TEXT, n_rows INTEGER, cache_path TEXT, cache_sha256 TEXT,
        note TEXT, PRIMARY KEY (catalogue, field_key))""")
    con.commit()
    fields = _field_centres(con)
    print(f"fields: {len(fields)}")
    # PRIMARY catalogue first, all fields, THEN the secondary.  The order
    # matters operationally: the ESA archive was answering in minutes while
    # VizieR answered in seconds during this build, and the primary tie
    # should not wait behind a cross-check.
    for cat in sorted(_catalogues(args),
                      key=lambda c: (c != ct.PRIMARY_CATALOGUE, c)):
        for tk, ra, dec in fields:
            cache = Path(args.cache) / cat / f"{tk}.json.gz"
            have = con.execute("SELECT n_rows FROM cv_cat_fetch WHERE "
                               "catalogue=? AND field_key=?",
                               (cat, tk)).fetchone()
            if have and cache.exists() and not args.force:
                print(f"  {cat:10s} {tk}: cached ({have[0]} rows) — skipping")
                continue
            print(f"  {cat:10s} {tk}: pulling "
                  f"({ra:.4f}, {dec:.4f}) r={ct.CONE_RADIUS_DEG} deg",
                  flush=True)
            t0 = time.time()
            note = ""
            try:
                if cat == "refcat2":
                    cols, query = fetch_refcat2(ra, dec, ct.CONE_RADIUS_DEG,
                                                ct.CAT_MAG_MAX,
                                                args.retries)
                    service = VIZIER_TAP_URL
                else:
                    gcols, qs = gaia_cone_sliced(ra, dec, ct.CONE_RADIUS_DEG,
                                                 GAIA_G_LIMIT, args.retries)
                    sids = [int(s) for s in gcols["source_id"]]
                    print(f"      cone: {len(sids)} sources; "
                          f"pulling synthetic photometry", flush=True)
                    gspc, qs2 = fetch_gspc(sids, args.retries)
                    cols = _gspc_to_frame(gcols, gspc)
                    query = "\n".join(qs + qs2)
                    service = "https://gea.esac.esa.int/tap-server/tap"
                    n_syn = int(np.sum(np.isfinite(
                        np.asarray(cols["g_sdss_mag"], dtype=float))))
                    note = (f"{n_syn} of {len(sids)} cone sources carry "
                            f"standardised synthetic photometry")
            except Exception as e:
                print(f"      FAILED after "
                      f"{args.retries or NET_RETRIES} attempts: "
                      f"{type(e).__name__}: {e}", flush=True)
                con.execute(
                    "INSERT OR REPLACE INTO cv_cat_fetch VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (cat, tk, ra, dec, ct.CONE_RADIUS_DEG, ct.CAT_MAG_MAX,
                     "", "", datetime.now(timezone.utc).isoformat(), 0, "",
                     "", f"FETCH FAILED: {type(e).__name__}: {e}"))
                con.commit()
                continue
            n = len(cols["ra"])
            payload = {"catalogue": cat, "field_key": tk,
                       "ra_deg": ra, "dec_deg": dec,
                       "radius_deg": ct.CONE_RADIUS_DEG,
                       "mag_limit": ct.CAT_MAG_MAX,
                       "epoch_jyr": CAT_EPOCH_JYR[cat],
                       "pulled_utc": datetime.now(timezone.utc).isoformat(),
                       "query": query, "columns": cols}
            sha = atomic_write_gz(cache, payload)
            con.execute("INSERT OR REPLACE INTO cv_cat_fetch VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        (cat, tk, ra, dec, ct.CONE_RADIUS_DEG, ct.CAT_MAG_MAX,
                         service, query, payload["pulled_utc"], n,
                         str(cache.relative_to(REPO_ROOT)), sha, note))
            con.commit()
            print(f"      {n} rows in {time.time() - t0:.0f}s -> "
                  f"{cache.name}  {note}", flush=True)
    meta_write(con, {"stage_fetch": datetime.now(timezone.utc).isoformat()})
    con.close()


# ===========================================================================
# Stage: match -- catalogue sources onto each block's reference stars
# ===========================================================================
MATCH_SQL = """CREATE TABLE IF NOT EXISTS cv_cat_match (
    catalogue TEXT, target_key TEXT, era_id INTEGER, star_id INTEGER,
    cat_row INTEGER, sep_arcsec REAL, second_sep_arcsec REAL,
    cat_ra REAL, cat_dec REAL, cat_gmag REAL,
    nn_sep_arcsec REAL, nn_dmag REAL,
    bright_sep_arcsec REAL, bright_dmag REAL,
    cat_flag INTEGER, epoch_jyr REAL,
    PRIMARY KEY (catalogue, target_key, era_id, star_id))"""

#: The astrometric zero point of each block, measured against the catalogue.
#: A row exists for EVERY (catalogue, block) that was matched, including the
#: ones where nothing was corrected -- "we measured it and it was 0.2 arcsec"
#: is evidence, and a table that only recorded the corrections would leave a
#: reader unable to tell a clean block from an unexamined one.
ASTROM_SQL = """CREATE TABLE IF NOT EXISTS cv_cat_astrom (
    catalogue TEXT, target_key TEXT, era_id INTEGER,
    method TEXT, n_paired INTEGER,
    dra_arcsec REAL, ddec_arcsec REAL, offset_arcsec REAL,
    scatter_arcsec REAL, applied INTEGER,
    n_match_before INTEGER, n_match_after INTEGER, n_stars INTEGER,
    reason TEXT,
    PRIMARY KEY (catalogue, target_key, era_id))"""


def _block_epoch(con, tk: str, era: int) -> float:
    """Median observation epoch of a (target, era), in Julian years."""
    r = con.execute("""SELECT avg(bjd_tdb) FROM cv_frames
                       WHERE target_key=? AND era_id=? AND bjd_tdb IS NOT NULL
                    """, (tk, era)).fetchone()
    return epoch_jyr(r[0]) if r and r[0] else 2024.0


#: REFCAT2's own contamination radius (arcsec) below which a source is
#: called blended.  ``R1`` is the radius at which the summed flux of the
#: neighbours equals the star's own; inside the 6-arcsec aperture that is a
#: 0.75 mag error, forty times the accuracy goal.  Sources with no
#: contaminating neighbour carry the sentinel 99.9.
REFCAT2_R1_MIN_ARCSEC = ct.BLEND_APERTURE_ARCSEC

#: Gaia's colour-excess consistency statistic ``C*``.  It is the
#: catalogue's own blend detector: BP and RP are measured through 3.5 x 2.1
#: arcsec windows, so a neighbour inflates them relative to G.  0.05 is
#: loose for bright stars and about right at G = 17, which is where these
#: comparison stars live.
GSPC_CSTAR_MAX = 0.05


def _catalogue_flag(cat: str, cols: dict) -> np.ndarray:
    """Per-source 'the catalogue itself distrusts this photometry' flag.

    BOTH catalogues are asked the SAME question -- "did YOU, with better
    angular resolution than this telescope, detect a contaminant or a
    photometry you would not stand behind?" -- and each answers it in its
    own native currency:

    ``refcat2``    ``R1``, the measured radius at which neighbour flux
                   equals the star's own, plus a missing/zero magnitude
                   error (which is how REFCAT2 marks a band with no real
                   measurement behind it).
    ``gaia_gspc``  the per-band standardisation validity flag, plus ``C*``,
                   the BP/RP colour-excess statistic that betrays a
                   neighbour inside Gaia's own 3.5-arcsec BP/RP window.

    Deliberately NOT USED: REFCAT2's ``dupvar``.  It reads 2 for 171 of the
    173 sources in a test cone and 6 for the other two -- i.e. it separates
    nothing at all in this data, so applying it as a veto would delete
    every tie star for a reason that is not a reason.  It is fetched and
    stored so a reader can see that it was looked at and discarded, which
    is what "we checked" has to mean.
    """
    n = len(cols["ra"])
    if cat == "refcat2":
        r1 = np.asarray(cols.get("r1", [99.9] * n), dtype=float)
        blended = np.isfinite(r1) & (r1 < REFCAT2_R1_MIN_ARCSEC)
        return np.where(blended, 1, 0).astype(int)
    fl = np.asarray(cols.get("g_sdss_flag", [1] * n), dtype=float)
    cs = np.asarray(cols.get("c_star", [np.nan] * n), dtype=float)
    bad_cs = np.isfinite(cs) & (np.abs(cs) > GSPC_CSTAR_MAX)
    return np.where((np.isfinite(fl) & (fl != 0)) | bad_cs, 1, 0).astype(int)


def cmd_match(args) -> None:
    """Match each block's reference stars to each catalogue.

    Sky positions for the reference stars come from ``cv_ref_stars.ra_deg``,
    which the field-tie stage filled EITHER from the S1 plate solution (when
    the reference frame has one) OR from the Gaia similarity fit (when it
    does not).  Which route a block took is recorded in
    ``cv_field_tie.method`` and carried into every tie row as
    ``astrom_source``, because a tie standing on a triangle fit deserves to
    be read differently from one standing on a plate solution.

    THE ASTROMETRIC ZERO POINT IS MEASURED FIRST.  Neither upstream route
    has an absolute reference -- a similarity fit onto a Gaia cone can be
    internally perfect and displaced as a whole, and nothing before this
    stage could see that.  So before any photometric match is attempted,
    every block is paired against the catalogue at a deliberately loose
    radius and the rigid part of the disagreement is measured
    (``cattie.rigid_offset``).  It is REMOVED only when it is bigger than
    the photometric tolerance, coherent across stars, and measured on
    enough of them; otherwise the positions are used exactly as written and
    ``cv_cat_astrom.reason`` records the decision.  This stage never writes
    ``cv_ref_stars`` -- the correction lives in the tie's own matching,
    because the field tie's ledger belongs to the field-tie stage.
    """
    con = connect(args.db)
    con.execute(MATCH_SQL)
    con.execute(ASTROM_SQL)
    con.commit()
    ft_method = {(r[0], r[1]): r[2] for r in con.execute(
        "SELECT target_key, era_id, method FROM cv_field_tie")}
    blocks = con.execute("""SELECT r.target_key, r.era_id
                            FROM cv_ref r ORDER BY 1, 2""").fetchall()
    for tk, era in blocks:
        stars = con.execute("""SELECT star_id, ra_deg, dec_deg FROM cv_ref_stars
                               WHERE target_key=? AND era_id=?
                               ORDER BY star_id""", (tk, era)).fetchall()
        have_pos = [s for s in stars if s[1] is not None]
        if not have_pos:
            print(f"  {tk} e{era}: no sky positions for its reference stars "
                  f"— cannot be matched to any catalogue")
            continue
        ep = _block_epoch(con, tk, era)
        sids = np.array([s[0] for s in have_pos], dtype=int)
        sra = np.array([s[1] for s in have_pos], dtype=float)
        sdec = np.array([s[2] for s in have_pos], dtype=float)
        # The astrometric offset is measured on the VETTED stars only --
        # every star any series of this block promoted to comparison or
        # check.  It has to be: a reference list can be mostly spurious
        # detections (EU UMa era 78 carries 2,746 entries and 49 real
        # stars), and pairing noise against a catalogue at a 15-arcsec
        # radius produces a scatter that correctly refuses every offset,
        # including the real one.  Measured on the 49 vetted stars the same
        # block yields a 5.20-arcsec offset with 0.45 arcsec of scatter.
        vet = {r[0] for r in con.execute(
            """SELECT DISTINCT star_id FROM cv_stars
               WHERE target_key=? AND era_id=? AND role IN ('comp','check')""",
            (tk, era))}
        use = np.array([int(s) in vet for s in sids]) if vet \
            else np.ones(len(sids), dtype=bool)
        if not use.any():
            use = np.ones(len(sids), dtype=bool)
        for cat in _catalogues(args):
            done = con.execute("SELECT count(*) FROM cv_cat_match WHERE "
                               "catalogue=? AND target_key=? AND era_id=?",
                               (cat, tk, era)).fetchone()[0]
            if done and not args.force:
                print(f"  {cat:10s} {tk} e{era}: {done} matches — skipping")
                continue
            payload = read_cache(Path(args.cache) / cat / f"{tk}.json.gz")
            if payload is None:
                print(f"  {cat:10s} {tk} e{era}: no cache — run fetch first")
                continue
            cols = payload["columns"]
            cra, cdec = propagate(cols["ra"], cols["dec"],
                                  cols.get("pmra"), cols.get("pmdec"),
                                  payload["epoch_jyr"], ep)
            cgm = np.asarray(cols["gaia_gmag"], dtype=float)
            flag = _catalogue_flag(cat, cols)
            # ---- astrometric zero point of this block, then the match -----
            n_before = int((ct.match_by_sky(sra, sdec, cra, cdec)[0] >= 0).sum())
            off = ct.rigid_offset(sra[use], sdec[use], cra, cdec)
            ura, udec = (ct.apply_offset(sra, sdec, off.dra_arcsec,
                                         off.ddec_arcsec)
                         if off.applied else (sra, sdec))
            idx, sep, sep2 = ct.match_by_sky(ura, udec, cra, cdec)
            hit = idx >= 0
            con.execute(
                "INSERT OR REPLACE INTO cv_cat_astrom VALUES "
                "(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (cat, tk, int(era), ft_method.get((tk, era)), off.n,
                 fnum(off.dra_arcsec), fnum(off.ddec_arcsec),
                 fnum(off.size_arcsec), fnum(off.scatter_arcsec),
                 int(off.applied), n_before, int(hit.sum()), int(use.sum()),
                 off.reason))
            if off.applied:
                print(f"  {cat:10s} {tk} e{era}: ASTROMETRY REFINED — "
                      f"{off.reason}; matches {n_before} -> {int(hit.sum())}",
                      flush=True)
            # Blend CENSUS, in Gaia G: a single band-independent currency in
            # which both catalogues answer the identical question, stored so
            # the two can be compared row for row.  It is NOT the veto.  The
            # veto has to be taken in the band being tied -- what matters is
            # whether the PS1 g magnitude describes the flux our g-band
            # aperture collected -- and ``solve`` recomputes these same
            # metrics per band for that purpose.  Review found the earlier
            # code vetoing in G and fitting in g/r/i, which let 46
            # equal-brightness close pairs (equal in r, two and a half
            # magnitudes apart in G) into the VV Pup fits.
            nm = ct.neighbour_metrics(
                cra[idx[hit]], cdec[idx[hit]], cgm[idx[hit]],
                cra, cdec, cgm, self_index=idx[hit])
            rows = []
            k = 0
            for n, sid in enumerate(sids):
                if not hit[n]:
                    continue
                j = int(idx[n])
                rows.append((cat, tk, int(era), int(sid), j,
                             fnum(sep[n]), fnum(sep2[n]),
                             fnum(cra[j]), fnum(cdec[j]), fnum(cgm[j]),
                             fnum(nm["nn_sep"][k]), fnum(nm["nn_dmag"][k]),
                             fnum(nm["bright_sep"][k]),
                             fnum(nm["bright_dmag"][k]),
                             int(flag[j]), float(ep)))
                k += 1
            con.execute("DELETE FROM cv_cat_match WHERE catalogue=? "
                        "AND target_key=? AND era_id=?", (cat, tk, era))
            con.executemany("INSERT INTO cv_cat_match VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
            con.commit()
            rate = 100.0 * len(rows) / max(len(sids), 1)
            print(f"  {cat:10s} {tk} e{era}: {len(rows)}/{len(sids)} ref "
                  f"stars matched ({rate:.1f}%), epoch {ep:.2f}", flush=True)
    meta_write(con, {"stage_match": datetime.now(timezone.utc).isoformat()})
    con.close()


# ===========================================================================
# Stage: solve -- the robust ZP + colour-term fit
# ===========================================================================
TIE_SQL = """CREATE TABLE IF NOT EXISTS cv_cattie (
    series_key TEXT, target_key TEXT, era_id INTEGER, filter TEXT,
    catalogue TEXT, hypothesis TEXT, band TEXT, band_system TEXT,
    colour_label TEXT, astrom_source TEXT,
    n_candidates INTEGER, n_clean INTEGER, n_fit INTEGER, n_clipped INTEGER,
    n_check INTEGER,
    zp REAL, zp_err REAL, colour_term REAL, colour_err REAL, colour_ref REAL,
    resid_rms REAL, resid_mad REAL, chi2nu REAL,
    colour_min REAL, colour_max REAL, colour_p05 REAL, colour_p95 REAL,
    check_rms REAL, check_rms_clip REAL, check_median REAL,
    check_mad REAL, n_check_outlier INTEGER,
    target_colour REAL, target_colour_source TEXT, colour_position TEXT,
    extrap_err REAL, is_primary INTEGER, verdict TEXT, note TEXT,
    n_colour_pairs INTEGER, colour_scatter REAL,
    PRIMARY KEY (series_key, catalogue, band))"""

#: Columns added after the table first shipped.  ``CREATE TABLE IF NOT
#: EXISTS`` will not add them to a database built by the previous version,
#: and silently writing 39 values into a 41-column row is exactly the kind
#: of failure that produces a plausible wrong number, so the migration is
#: explicit and runs on every invocation.
TIE_ADDED_COLUMNS = (("n_colour_pairs", "INTEGER"), ("colour_scatter", "REAL"))


def _ensure_columns(con, table: str,
                    columns: tuple[tuple[str, str], ...]) -> None:
    """Add missing columns to an existing table (idempotent)."""
    have = {r[1] for r in con.execute(f"PRAGMA table_info({table})")}
    for name, decl in columns:
        if name not in have:
            con.execute(f"ALTER TABLE {table} ADD COLUMN {name} {decl}")
    con.commit()

STAR_SQL = """CREATE TABLE IF NOT EXISTS cv_cattie_star (
    series_key TEXT, catalogue TEXT, band TEXT, star_id INTEGER,
    role TEXT, in_fit INTEGER, ens_mag REAL, ens_err REAL,
    cat_mag REAL, cat_err REAL, colour REAL, delta REAL, resid REAL,
    x REAL, y REAL,
    PRIMARY KEY (series_key, catalogue, band, star_id))"""

VETO_SQL = """CREATE TABLE IF NOT EXISTS cv_cattie_veto (
    series_key TEXT, catalogue TEXT, band TEXT, reason TEXT, n INTEGER,
    PRIMARY KEY (series_key, catalogue, band, reason))"""


def _series_rows(con) -> list[tuple]:
    """Every SOLVED series -- the only ones that have magnitudes to tie."""
    return con.execute("""SELECT series_key, target_key, era_id, filter
                          FROM cv_series WHERE status='solved'
                          ORDER BY target_key, era_id, filter""").fetchall()


def _saturation_census(con, series_key: str) -> dict[int, tuple[float, float]]:
    """Per star: fraction of measurements saturated, and fraction near veto.

    Read from the DETECTIONS, not from the light curve, because the veto is
    a property of the pixels: ``cv_detections.saturated`` was set against
    the frame's own applied veto in the frame's own units, which differ by
    readout mode (High Gain's scale ends at 3,496 ADU while Mode0's ends at
    65,535).  ``peak`` against 0.9x that veto catches the non-linear
    shoulder below the clip, which no flag records.
    """
    out: dict[int, tuple[float, float]] = {}
    q = con.execute("""
        SELECT d.star_id, COUNT(*), SUM(d.saturated),
               SUM(CASE WHEN d.peak >= ? * COALESCE(f.veto_applied_adu,
                                                    f.veto_adu)
                        THEN 1 ELSE 0 END)
        FROM cv_frames f JOIN cv_detections d ON d.frame_id = f.frame_id
        WHERE f.series_key = ? AND d.star_id IS NOT NULL
        GROUP BY d.star_id""", (ct.NEAR_VETO_FRAC, series_key))
    for sid, n, nsat, nnear in q:
        n = max(int(n), 1)
        out[int(sid)] = (float(nsat or 0) / n, float(nnear or 0) / n)
    return out


def _colour_refs(con, tk: str, era: int, cat: str, cols: dict
                 ) -> dict[str, float]:
    """Median catalogue colour of a BLOCK, per colour index.

    Computed once per (target, era, catalogue) and shared by every filter of
    that block.  Two reasons, and the second is the load-bearing one:
    a shared reference makes the g and r zero points directly comparable,
    and :func:`macro_phot.cattie.target_colour_solve` -- the closed-form
    that places the target on the colour axis -- is only valid when both
    bands were centred on the SAME reference colour.
    """
    rows = con.execute("""SELECT cat_row FROM cv_cat_match
                          WHERE catalogue=? AND target_key=? AND era_id=?""",
                       (cat, tk, era)).fetchall()
    idx = np.array([r[0] for r in rows], dtype=int)
    out: dict[str, float] = {}
    if idx.size == 0:
        return out
    seen = set()
    for specs in ct.BAND_CANDIDATES.values():
        for sp in specs:
            if sp.catalogue != cat or sp.colour_label in seen:
                continue
            seen.add(sp.colour_label)
            b = np.asarray(cols.get(sp.colour_blue, []), dtype=float)
            r = np.asarray(cols.get(sp.colour_red, []), dtype=float)
            if b.size == 0 or r.size == 0:
                continue
            c = b[idx] - r[idx]
            c = c[np.isfinite(c)]
            if c.size:
                out[sp.colour_label] = float(np.median(c))
    return out


def _band_blend(con, cache: dict, cat: str, tk: str, era: int,
                cols: dict, mag_col: str) -> dict:
    """Blend metrology for one block IN THE BAND BEING TIED, memoised.

    ``cv_cat_match`` stores the same quantities in Gaia G, which is a good
    census (one currency, both catalogues, comparable row for row) and a bad
    VETO.  Ruling 4 asks whether the catalogue magnitude in the tie band
    describes the flux this telescope's aperture collected, and that
    question is asked in the tie band or not at all: a pair Gaia resolves
    and PS1 does not is two and a half magnitudes apart in G and equal in r,
    so a G-band gate passes it straight into an r-band fit.  Review found 46
    such pairs inside the VV Pup fits.

    Computed once per (catalogue, target, era, magnitude column) over every
    star the block matched, and returned as ``cat_row -> metrics`` so each
    series can index the subset it uses.
    """
    key = (cat, tk, int(era), mag_col)
    if key in cache:
        return cache[key]
    cra = np.asarray(cols["ra"], dtype=float)
    cdec = np.asarray(cols["dec"], dtype=float)
    cm = np.asarray(cols.get(mag_col, []), dtype=float)
    rows = np.array(sorted({r[0] for r in con.execute(
        "SELECT cat_row FROM cv_cat_match WHERE catalogue=? AND target_key=? "
        "AND era_id=?", (cat, tk, int(era)))}), dtype=int)
    if cm.size == 0 or rows.size == 0:
        cache[key] = {}
        return cache[key]
    nm = ct.neighbour_metrics(cra[rows], cdec[rows], cm[rows],
                              cra, cdec, cm, self_index=rows)
    # ``aper_*`` rather than ``nn_*``: the veto is about the WORST
    # contaminant inside the aperture, which is not in general the nearest
    # source.  See cattie.neighbour_metrics for the 55 stars that taught it.
    cache[key] = {int(r): (float(nm["aper_sep"][i]), float(nm["aper_dmag"][i]),
                           float(nm["bright_sep"][i]),
                           float(nm["bright_dmag"][i]))
                  for i, r in enumerate(rows)}
    return cache[key]


def cmd_solve(args) -> None:
    """Fit ZP0 and the colour term for every (series, catalogue, band).

    One row of ``cv_cattie`` per hypothesis, including the ones the data is
    asked to REJECT (the Johnson-Cousins readings of the uppercase filter
    labels).  Storing the losing hypotheses is deliberate: 'we checked, and
    here is by how much it lost' is evidence, while 'we assumed' is not.
    """
    con = connect(args.db)
    for sql in (TIE_SQL, STAR_SQL, VETO_SQL):
        con.execute(sql)
    con.commit()
    _ensure_columns(con, "cv_cattie", TIE_ADDED_COLUMNS)
    # Drop rows for series this run will not revisit.  Same lesson as the
    # per-block delete below: an upsert leaves behind whatever it no longer
    # writes, and a stale row is worse than a missing one because it looks
    # like an answer.
    for tab in ("cv_cattie", "cv_cattie_star", "cv_cattie_veto"):
        con.execute(f"""DELETE FROM {tab} WHERE series_key NOT IN
                        (SELECT series_key FROM cv_series
                         WHERE status='solved')""")
    con.commit()
    astrom = {(r[0], r[1]): (r[2], r[3]) for r in con.execute(
        "SELECT target_key, era_id, method, status FROM cv_field_tie")}
    caches: dict[tuple[str, str], dict] = {}
    crefs: dict[tuple[str, int, str], dict[str, float]] = {}
    blends: dict[tuple, dict] = {}

    for skey, tk, era, filt in _series_rows(con):
        sat = _saturation_census(con, skey)
        # Candidate stars: comparison + check.  Never the target (ruling 2),
        # never the stars the ensemble's stability iteration dropped, never
        # 'field' stars that were measured but never vetted for constancy.
        cand = con.execute("""SELECT star_id, role, mean_mag, rms, nobs
                              FROM cv_stars WHERE series_key=?
                              AND role IN ('comp','check')
                              ORDER BY star_id""", (skey,)).fetchall()
        if not cand:
            continue
        sid = np.array([c[0] for c in cand], dtype=int)
        role = np.array([c[1] for c in cand], dtype=object)
        ens = np.array([c[2] if c[2] is not None else np.nan
                        for c in cand], dtype=float)
        rms = np.array([c[3] if c[3] is not None else np.nan
                        for c in cand], dtype=float)
        nob = np.array([max(int(c[4] or 1), 1) for c in cand], dtype=float)
        # Error on the MEAN magnitude, with the pipeline's own systematic
        # floor: a bright star's formal error is not the accuracy of its
        # mean, because flat-field and PSF systematics do not average down.
        ens_err = np.sqrt(np.square(rms) / nob + 0.005 ** 2)
        satf = np.array([sat.get(int(s), (0.0, 0.0))[0] for s in sid])
        nearf = np.array([sat.get(int(s), (0.0, 0.0))[1] for s in sid])
        xy = {r[0]: (r[1], r[2]) for r in con.execute(
            "SELECT star_id, x, y FROM cv_ref_stars WHERE target_key=? "
            "AND era_id=?", (tk, era))}
        px = np.array([xy.get(int(s), (np.nan, np.nan))[0] for s in sid])
        py = np.array([xy.get(int(s), (np.nan, np.nan))[1] for s in sid])
        held = ct.holdout_mask(sid.tolist(), skey) | (role == "check")

        for cat in ct.CATALOGUES:
            specs = ct.band_candidates(filt, cat)
            if not specs:
                continue
            key = (cat, tk)
            if key not in caches:
                caches[key] = read_cache(
                    Path(args.cache) / cat / f"{tk}.json.gz") or {}
            payload = caches[key]
            if not payload:
                continue
            cols = payload["columns"]
            ckey = (tk, era, cat)
            if ckey not in crefs:
                crefs[ckey] = _colour_refs(con, tk, era, cat, cols)
            cref_map = crefs[ckey]
            m = {r[0]: r[1:] for r in con.execute(
                """SELECT star_id, cat_row, sep_arcsec, second_sep_arcsec,
                          nn_sep_arcsec, nn_dmag, bright_sep_arcsec,
                          bright_dmag, cat_flag
                   FROM cv_cat_match WHERE catalogue=? AND target_key=?
                   AND era_id=?""", (cat, tk, era))}
            if not m:
                continue
            row = np.array([m.get(int(s)) is not None for s in sid])
            cat_row = np.array([m[int(s)][0] if m.get(int(s)) else -1
                                for s in sid], dtype=int)
            g = lambda k: np.array(                       # noqa: E731
                [m[int(s)][k] if m.get(int(s)) else np.nan for s in sid],
                dtype=float)
            sep, sep2 = g(1), g(2)
            cflag = np.array([m[int(s)][7] if m.get(int(s)) else 1
                              for s in sid], dtype=float)
            cflag = np.where(row, cflag, 1.0)

            for sp in specs:
                cm_all = np.asarray(cols.get(sp.mag_col, []), dtype=float)
                if cm_all.size == 0:
                    continue
                ce_all = (np.asarray(cols.get(sp.err_col, []), dtype=float)
                          if sp.err_col and sp.err_col in cols
                          else np.full(cm_all.shape, np.nan))
                cb = np.asarray(cols.get(sp.colour_blue, []), dtype=float)
                cr = np.asarray(cols.get(sp.colour_red, []), dtype=float)
                take = lambda a: np.where(                 # noqa: E731
                    row, a[np.clip(cat_row, 0, len(a) - 1)], np.nan)
                cat_mag = take(cm_all)
                cat_err = take(ce_all)
                colour = take(cb) - take(cr) if cb.size and cr.size \
                    else np.full(len(sid), np.nan)
                # ---- blend metrology IN THIS BAND (ruling 4) -------------
                inf4 = (np.inf, np.inf, np.inf, np.inf)

                def _pick(col):                            # noqa: E306
                    bb = _band_blend(con, blends, cat, tk, era, cols, col)
                    if not len(sid):
                        return (np.array([]),) * 4
                    return tuple(np.array(a) for a in zip(
                        *[bb.get(int(cat_row[n]), inf4) if row[n] else inf4
                          for n in range(len(sid))]))
                nn_sep, nn_dmag, br_sep, br_dmag = _pick(sp.mag_col)
                # The SAME rule read in Gaia G, for the comparison below.
                gn_sep, gn_dmag = _pick("gaia_gmag")[:2]
                keep, census = ct.clean_mask(
                    saturated_frac=satf, near_veto_frac=nearf,
                    blend_sep_arcsec=nn_sep, blend_dmag=nn_dmag,
                    annulus_sep_arcsec=br_sep, annulus_dmag=br_dmag,
                    second_sep_arcsec=sep2, match_sep_arcsec=sep,
                    cat_mag=cat_mag, cat_colour=colour, cat_flag=cflag)
                # How the SAME gate would have read in Gaia G -- the number
                # that says whether taking the veto in the tie band mattered.
                # ``band_only`` is the one that answers it: stars the tie
                # band rejects and Gaia G would have waved through, which is
                # exactly the population review found inside the VV Pup fits.
                g_bad = ((gn_sep < ct.BLEND_APERTURE_ARCSEC)
                         & (gn_dmag < ct.BLEND_DMAG))
                b_bad = ((nn_sep < ct.BLEND_APERTURE_ARCSEC)
                         & (nn_dmag < ct.BLEND_DMAG))
                g_blend, b_blend = int(g_bad.sum()), int(b_bad.sum())
                band_only = int((b_bad & ~g_bad).sum())
                keep &= np.isfinite(ens)
                in_fit = keep & ~held
                in_chk = keep & held
                c_ref = cref_map.get(sp.colour_label, float("nan"))
                if not math.isfinite(c_ref):
                    c_ref = float(np.median(colour[keep])) if keep.any() \
                        else float("nan")
                delta = ens - cat_mag
                sig = np.sqrt(np.square(ens_err)
                              + np.square(np.nan_to_num(cat_err, nan=0.01)))
                fit = ct.robust_line_fit(colour[in_fit], delta[in_fit],
                                         sig[in_fit], x_ref=c_ref)
                cmin, cmax, p05, p95 = ct.colour_range(colour[in_fit])
                nat_chk = ens[in_chk] - fit.zp
                chk = ct.check_accuracy(nat_chk, cat_mag[in_chk],
                                        colour[in_chk], fit.zp, fit.slope,
                                        c_ref)
                # Graded on the CLIPPED check RMS -- see
                # cattie.block_verdict for why, and note that the raw
                # RMS and the outlier count are stored beside it so the
                # clip can never hide what it removed.
                verdict = ct.block_verdict(int(in_fit.sum()), chk["n"],
                                           chk["rms_clip"],
                                           fit.converged)
                method, status = astrom.get((tk, era), (None, None))
                is_primary = int(cat == ct.PRIMARY_CATALOGUE
                                 and sp.hypothesis == "primary")
                if not is_primary and cat != ct.PRIMARY_CATALOGUE \
                        and sp.hypothesis == "primary" \
                        and not ct.band_candidates(filt,
                                                   ct.PRIMARY_CATALOGUE):
                    # The primary catalogue has no analogue for this filter
                    # (the y band).  The secondary is then not a cross-check
                    # but the tie itself, and the row says so.
                    is_primary = 1
                note = ""
                if int(in_fit.sum()) < ct.MIN_TIE_STARS:
                    note = (f"only {int(in_fit.sum())} clean tie stars "
                            f"(minimum {ct.MIN_TIE_STARS})")
                con.execute(
                    "INSERT OR REPLACE INTO cv_cattie VALUES "
                    "(" + ",".join("?" * 41) + ")",
                    (skey, tk, int(era), filt, cat, sp.hypothesis,
                     sp.mag_col, sp.system, sp.colour_label,
                     f"{method or 'none'}/{status or ''}",
                     int(len(sid)), int(keep.sum()), int(in_fit.sum()),
                     int(fit.n_clipped), int(chk["n"]),
                     fnum(fit.zp), fnum(fit.zp_err), fnum(fit.slope),
                     fnum(fit.slope_err), fnum(c_ref),
                     fnum(fit.resid_rms), fnum(fit.resid_mad),
                     fnum(fit.chi2nu),
                     fnum(cmin), fnum(cmax), fnum(p05), fnum(p95),
                     fnum(chk["rms"]), fnum(chk["rms_clip"]),
                     fnum(chk["median"]), fnum(chk["mad"]),
                     int(chk["n_outlier"]),
                     None, None, None, None,
                     is_primary, verdict, note, None, None))
                # DELETE BEFORE INSERT, and the reason is a bug this stage
                # actually shipped.  ``cv_cattie_star`` holds one row per
                # star that SURVIVED the gate, so a star that the gate
                # starts rejecting simply stops being written -- and its row
                # from the previous run, complete with ``in_fit = 1``,
                # stayed behind.  The fit itself was always correct (it runs
                # on arrays, not on this table), but the published per-star
                # table then disagreed with the published fit, and every
                # audit that reads the table -- including the regression
                # test that caught this -- was reading stars the solver had
                # already thrown out.  An upsert cannot express "these rows
                # and no others"; a delete can.
                con.execute("DELETE FROM cv_cattie_star WHERE series_key=? "
                            "AND catalogue=? AND band=?",
                            (skey, cat, sp.mag_col))
                con.execute("DELETE FROM cv_cattie_veto WHERE series_key=? "
                            "AND catalogue=? AND band=?",
                            (skey, cat, sp.mag_col))
                con.executemany(
                    "INSERT OR REPLACE INTO cv_cattie_star VALUES "
                    "(" + ",".join("?" * 15) + ")",
                    [(skey, cat, sp.mag_col, int(sid[n]), str(role[n]),
                      int(in_fit[n]), fnum(ens[n]), fnum(ens_err[n]),
                      fnum(cat_mag[n]), fnum(cat_err[n]), fnum(colour[n]),
                      fnum(delta[n]),
                      fnum(delta[n] - (fit.zp + fit.slope
                                       * (colour[n] - c_ref))),
                      fnum(px[n]), fnum(py[n]))
                     for n in range(len(sid)) if keep[n]])
                con.executemany(
                    "INSERT OR REPLACE INTO cv_cattie_veto VALUES (?,?,?,?,?)",
                    [(skey, cat, sp.mag_col, k, v) for k, v in (
                        ("candidates", census.n_candidates),
                        ("no_cat_mag", census.n_no_cat_mag),
                        ("catalogue_flag", census.n_flagged),
                        ("saturated", census.n_saturated),
                        ("near_veto", census.n_near_veto),
                        ("blend_aperture", census.n_blend_aperture),
                        ("blend_annulus", census.n_blend_annulus),
                        ("ambiguous", census.n_ambiguous),
                        ("clean", census.n_clean),
                        ("held_out", int(in_chk.sum())),
                        ("clipped_by_fit", int(fit.n_clipped)),
                        # Diagnostic pair, not part of the partition: how
                        # many candidates the aperture-blend rule catches in
                        # the TIE BAND, and how many the same rule would
                        # have caught reading Gaia G instead.
                        ("blend_aperture_tieband", b_blend),
                        ("blend_aperture_gband", g_blend),
                        ("blend_aperture_band_only", band_only))])
                con.commit()
                print(f"  {skey:16s} {cat:10s} {sp.mag_col:12s} "
                      f"[{sp.hypothesis:11s}] n_fit={int(in_fit.sum()):4d} "
                      f"ZP={fit.zp:8.4f}+/-{fit.zp_err:.4f} "
                      f"k={fit.slope:+.4f}+/-{fit.slope_err:.4f} "
                      f"rms={fit.resid_rms:.4f} "
                      f"check={chk['rms_clip']:.4f}"
                      f"/{chk['rms']:.4f}({chk['n']}"
                      f"-{chk['n_outlier']})  {verdict}",
                      flush=True)
    meta_write(con, {"stage_solve": datetime.now(timezone.utc).isoformat()})
    con.close()


# ===========================================================================
# Stage: validate -- be adversarial with our own result
# ===========================================================================
TREND_SQL = """CREATE TABLE IF NOT EXISTS cv_cattie_trend (
    series_key TEXT, catalogue TEXT, band TEXT, axis TEXT,
    slope REAL, slope_err REAL, n INTEGER, span REAL,
    significance REAL, swing REAL, significant INTEGER,
    PRIMARY KEY (series_key, catalogue, band, axis))"""

CROSS_SQL = """CREATE TABLE IF NOT EXISTS cv_cattie_cross (
    series_key TEXT, band_a TEXT, band_b TEXT, n_common INTEGER,
    d_zp REAL, d_colour_term REAL,
    star_offset_median REAL, star_offset_rms REAL, star_offset_mad REAL,
    note TEXT, PRIMARY KEY (series_key, band_a, band_b))"""


def cmd_validate(args) -> None:
    """Attack the tie from three sides, and write down what happens.

    1.  TRENDS.  Regress the tie residual on catalogue magnitude (which
        would expose detector non-linearity or a signal-dependent aperture
        correction) and fit a plane in detector (x, y) (which would expose
        a flat-field residual).  The characterization flagged a noise floor
        28 times above scintillation and nobody had looked at (x, y); this
        looks.
    2.  CROSS-CATALOGUE.  Where both catalogues tie the same series, the
        difference between their zero points -- and, more directly, the
        median offset between their magnitudes for the SAME stars -- is
        the systematic floor of any magnitude published here.
    3.  TARGET COLOUR.  State where each science target's colour falls
        relative to the range its own fit interpolated over.
    """
    con = connect(args.db)
    con.execute(TREND_SQL)
    con.execute(CROSS_SQL)
    con.commit()
    _ensure_columns(con, "cv_cattie", TIE_ADDED_COLUMNS)

    # ---- 1. trends ------------------------------------------------------
    ties = con.execute("""SELECT series_key, catalogue, band FROM cv_cattie
                          WHERE n_fit >= ? ORDER BY 1,2,3""",
                       (ct.MIN_TIE_STARS,)).fetchall()
    for skey, cat, band in ties:
        rows = con.execute("""SELECT cat_mag, resid, x, y FROM cv_cattie_star
                              WHERE series_key=? AND catalogue=? AND band=?
                              AND in_fit=1""", (skey, cat, band)).fetchall()
        if len(rows) < 6:
            continue
        cm = np.array([r[0] if r[0] is not None else np.nan for r in rows])
        rs = np.array([r[1] if r[1] is not None else np.nan for r in rows])
        xs = np.array([r[2] if r[2] is not None else np.nan for r in rows])
        ys = np.array([r[3] if r[3] is not None else np.nan for r in rows])
        # RADIUS from the detector centre.  Added after the (x, y) maps were
        # first drawn and showed a bullseye rather than a tilt: positive
        # residuals in the middle of the frame, negative at the corners, in
        # four different fields.  A plane cannot represent that shape, so a
        # plane fit UNDERSTATES it -- which is exactly how a real
        # illumination-correction error hides from a linear test.
        n1, n2 = con.execute(
            """SELECT max(naxis1), max(naxis2) FROM cv_frames
               WHERE series_key=?""", (skey,)).fetchone()
        cx = (n1 or 0) / 2.0 or float(np.nanmean(xs))
        cy = (n2 or 0) / 2.0 or float(np.nanmean(ys))
        rad = np.hypot(xs - cx, ys - cy)
        out = []
        for name, axis in (("cat_mag", cm), ("x", xs), ("y", ys),
                           ("radius", rad)):
            t = ct.residual_trend(name, axis, rs)
            out.append((skey, cat, band, t.name, fnum(t.slope),
                        fnum(t.slope_err), t.n, fnum(t.span),
                        fnum(t.significance), fnum(t.swing),
                        int(t.significant)))
        pf = ct.plane_fit(xs, ys, rs)
        out.append((skey, cat, band, "plane_xy", fnum(pf["cx"]),
                    fnum(pf["cx_err"]), pf["n"],
                    fnum(np.nanmax(xs) - np.nanmin(xs)),
                    fnum(pf["significance"]), fnum(pf["swing"]),
                    int(pf["significant"])))
        con.executemany("INSERT OR REPLACE INTO cv_cattie_trend VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?)", out)
    con.commit()
    n_sig = con.execute("SELECT count(*) FROM cv_cattie_trend "
                        "WHERE significant=1").fetchone()[0]
    n_all = con.execute("SELECT count(*) FROM cv_cattie_trend").fetchone()[0]
    print(f"  trends: {n_sig} of {n_all} regressions significant at "
          f"{ct.TREND_SIGMA:g} sigma")

    # ---- 2. cross-catalogue --------------------------------------------
    con.execute("DELETE FROM cv_cattie_cross")
    pairs = con.execute("""
        SELECT a.series_key, a.band, b.band, a.zp, b.zp,
               a.colour_term, b.colour_term
        FROM cv_cattie a JOIN cv_cattie b ON a.series_key=b.series_key
        WHERE a.catalogue='refcat2' AND b.catalogue='gaia_gspc'
          AND a.hypothesis='primary' AND b.hypothesis='primary'
          AND a.n_fit >= ? AND b.n_fit >= ?""",
        (ct.MIN_TIE_STARS, ct.MIN_TIE_STARS)).fetchall()
    for skey, ba, bb, zpa, zpb, ka, kb in pairs:
        sa = {r[0]: r[1] for r in con.execute(
            "SELECT star_id, cat_mag FROM cv_cattie_star WHERE series_key=? "
            "AND catalogue='refcat2' AND band=?", (skey, ba))}
        sb = {r[0]: r[1] for r in con.execute(
            "SELECT star_id, cat_mag FROM cv_cattie_star WHERE series_key=? "
            "AND catalogue='gaia_gspc' AND band=?", (skey, bb))}
        common = [(sa[k], sb[k]) for k in sa.keys() & sb.keys()
                  if sa[k] is not None and sb[k] is not None]
        if len(common) < 5:
            continue
        d = np.array([a - b for a, b in common], dtype=float)
        con.execute("INSERT OR REPLACE INTO cv_cattie_cross VALUES "
                    "(?,?,?,?,?,?,?,?,?,?)",
                    (skey, f"refcat2:{ba}", f"gaia_gspc:{bb}", len(d),
                     fnum(zpa - zpb) if (zpa is not None and zpb is not None)
                     else None,
                     fnum(ka - kb) if (ka is not None and kb is not None)
                     else None,
                     fnum(np.median(d)),
                     fnum(np.sqrt(np.mean(np.square(d)))),
                     fnum(1.4826 * np.median(np.abs(d - np.median(d)))),
                     "star-by-star offset is REFCAT2 minus Gaia synthetic "
                     "in the nominally same band; it contains the real "
                     "PS1-vs-SDSS system difference as well as either "
                     "catalogue's own error"))
    con.commit()
    print(f"  cross-catalogue: {len(pairs)} series tied by both catalogues")

    # ---- 3. target colours ---------------------------------------------
    _target_colours(con)
    meta_write(con, {"stage_validate": datetime.now(timezone.utc).isoformat()})
    con.close()


def _target_points(con, series_key: str) -> tuple[np.ndarray, np.ndarray]:
    """(time, magnitude) of every finite TARGET measurement in one series."""
    rows = con.execute("""SELECT bjd_tdb, mag FROM cv_lightcurve
                          WHERE series_key=? AND role='target'
                            AND mag IS NOT NULL AND bjd_tdb IS NOT NULL
                          ORDER BY bjd_tdb""", (series_key,)).fetchall()
    return (np.array([r[0] for r in rows], dtype=float),
            np.array([r[1] for r in rows], dtype=float))


def _target_colours(con) -> None:
    """Place each science target on its own fit's colour axis.

    Nothing here writes a magnitude -- ruling 1 stands -- it writes a
    POSITION on the colour axis, which is what decides whether the tie may
    be claimed to apply to the target at all.  The colour is then converted
    from the natural system to the catalogue system by a closed form in the
    two colour terms (``cattie.target_colour_solve``).

    THE COLOUR IS FORMED FROM PAIRED EPOCHS, NOT FROM CAMPAIGN MEANS, and
    the difference is not a refinement.  The first version of this function
    differenced the two filters' ENSEMBLE MEAN magnitudes -- the mean over
    every frame of each series.  For a constant star that is the colour.
    For a polar it is the difference between the target's mean state while
    the blue filter was on the wheel and its mean state while the red one
    was, and those campaigns need not overlap at all: VV Pup era 76's g
    points span 370 days and its r points 55 of them, and the recipe
    published g-r = -1.73 for an object whose colour on the epochs where
    both filters actually observed is +0.04.  That error decided the block's
    stated verdict ("extrapolated", with a 105 mmag extrapolation charge)
    and it was wrong by 1.77 magnitudes.

    So each blue measurement is paired with the nearest red measurement in
    time and only pairs within ``cattie.COLOUR_PAIR_TOL_DAYS`` are kept
    (:func:`macro_phot.cattie.paired_colour`).  Where the two filters never
    sampled the same time -- VV Pup era 72's g and r share no night at all,
    EU UMa era 76's r and i share none either -- the answer is ``unknown``,
    which is the same answer a single-filter era gets and for the same
    reason: this target has no measured colour in this era.  The number of
    surviving pairs and their scatter are stored beside the colour, because
    a CV's colour genuinely moves and one number without its spread would
    hide that.
    """
    # Preferred colour pairs, per colour index: (blue filter, red filter).
    PAIRS = {"g-r": [("g", "r"), ("G", "R")],
             "r-i": [("r", "i"), ("R", "I")],
             "i-z": [("i", "z"), ("I", "z")]}
    rows = con.execute("""SELECT series_key, target_key, era_id, filter,
                                 catalogue, band, colour_label, zp,
                                 colour_term, colour_ref, colour_min,
                                 colour_max, colour_p05, colour_p95,
                                 colour_err
                          FROM cv_cattie WHERE hypothesis='primary'
                          AND n_fit >= ?""", (ct.MIN_TIE_STARS,)).fetchall()
    by = {}
    for r in rows:
        by[(r[1], r[2], r[4], r[3])] = r
    pts: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for r in rows:
        skey, tk, era, filt, cat, band, clab, zp, k, cref = r[:10]
        cmin, cmax, p05, p95, kerr = r[10:15]
        tmag = con.execute("""SELECT mean_mag FROM cv_stars
                              WHERE series_key=? AND role='target'""",
                           (skey,)).fetchone()
        colour = float("nan")
        n_pairs, scatter = 0, float("nan")
        source = "single filter in this era — colour not measurable"
        for blue, red in PAIRS.get(clab, []):
            rb = by.get((tk, era, cat, blue))
            rr = by.get((tk, era, cat, red))
            if rb is None or rr is None:
                continue
            for s in (rb[0], rr[0]):
                if s not in pts:
                    pts[s] = _target_points(con, s)
            tb, mb = pts[rb[0]]
            tr, mr = pts[rr[0]]
            # Natural-system magnitudes: each filter's own ZP0 removed
            # BEFORE the difference, so the colour is on the catalogue's
            # zero point in both bands.
            pc = ct.paired_colour(tb, mb - rb[7], tr, mr - rr[7])
            colour = ct.target_colour_solve(pc.colour, rb[8], rr[8], cref)
            n_pairs, scatter = pc.n_pairs, pc.scatter
            source = (f"{blue}-{red}: {pc.note}"
                      + (f"; natural-system {pc.colour:+.3f} inverted through "
                         f"k({blue})={rb[8]:+.3f}, k({red})={rr[8]:+.3f}"
                         if math.isfinite(pc.colour) else ""))
            break
        pos = ct.colour_position(colour, cmin, cmax, p05, p95)
        # The extrapolation charge carries the colour term's OWN error
        # over the same lever arm: outside the fitted range the model is
        # uncertain both because it may be the wrong shape and because its
        # slope was never well determined out there.
        extrap = ct.colour_extrapolation_error(colour, cmin, cmax, k or 0.0,
                                               kerr or 0.0)
        if tmag is None or tmag[0] is None:
            source = "target undetected in this series"
            pos = "unknown"
            colour, extrap = float("nan"), float("nan")
        con.execute("""UPDATE cv_cattie SET target_colour=?,
                       target_colour_source=?, colour_position=?,
                       extrap_err=?, n_colour_pairs=?, colour_scatter=?
                       WHERE series_key=? AND catalogue=?
                       AND band=?""",
                    (fnum(colour), source, pos, fnum(extrap),
                     int(n_pairs), fnum(scatter), skey, cat, band))
    con.commit()


# ===========================================================================
# Stage: apply -- the calibrated light curve
# ===========================================================================
def cmd_apply(args) -> None:
    """Write ``cv_lightcurve.cal_mag`` = ``mag - ZP0`` for every tied series.

    WHY A COLUMN AND NOT A TABLE.  ``cv_lightcurve`` holds 3.18 million
    rows.  A parallel table would double that on a spinning disk and, worse,
    would make every consumer perform a four-column join to pair a
    calibrated point with its own natural-system parent -- and a join is a
    place where the wrong two rows can meet.  A column cannot be mispaired
    with its own row.  The relative magnitude ``mag`` is untouched beside
    it, so nothing that was measured is lost or overwritten.

    WHY NO COLOUR TERM IS APPLIED.  ``cal_mag`` is the telescope's NATURAL
    system, zero-pointed to the catalogue at the tie stars' median colour.
    Applying ``k`` to the target would be a transformation, which ruling 1
    forbids for exactly these objects.  Applying it to the comparison stars
    but not the target would put two kinds of magnitude in one column.  So:
    one rule for every row, and ``k`` published beside it.

    WHY NO PER-POINT ERROR COLUMN.  ``zp_err`` is a SYSTEMATIC shared by
    every point of the block.  Adding it in quadrature to each point would
    corrupt the differential errors that every precision, period and timing
    result on this data set rests on, and would then be double-counted by
    anyone comparing two points of the same curve.  It stays one number per
    block in ``cv_cattie.zp_err``, to be applied once, to the block.
    """
    con = connect(args.db)
    have = {r[1] for r in con.execute("PRAGMA table_info(cv_lightcurve)")}
    if "cal_mag" not in have:
        con.execute("ALTER TABLE cv_lightcurve ADD COLUMN cal_mag REAL")
        con.commit()
        print("  added column cv_lightcurve.cal_mag")
    ties = con.execute("""SELECT series_key, zp, catalogue, band, verdict
                          FROM cv_cattie WHERE is_primary=1
                          AND verdict LIKE 'TIED%' AND zp IS NOT NULL
                          ORDER BY series_key""").fetchall()
    dupes = con.execute(
        """SELECT series_key, count(*) FROM cv_cattie WHERE is_primary=1
           GROUP BY 1 HAVING count(*) > 1""").fetchall()
    if dupes:
        # INSERT OR REPLACE into the map below would silently pick whichever
        # row came last, i.e. calibrate a light curve against an arbitrary
        # one of two zero points.  Refuse instead.
        raise SystemExit(
            "REFUSING TO APPLY: these series carry more than one PRIMARY "
            "tie, so there is no single zero point to apply — "
            + ", ".join(f"{k} ({n})" for k, n in dupes))
    con.execute("DROP TABLE IF EXISTS temp.zp_map")
    con.execute("CREATE TEMP TABLE zp_map (series_key TEXT PRIMARY KEY, "
                "zp REAL)")
    con.executemany("INSERT OR REPLACE INTO temp.zp_map VALUES (?,?)",
                    [(t[0], t[1]) for t in ties])
    con.commit()
    t0 = time.time()
    con.execute("""UPDATE cv_lightcurve
                   SET cal_mag = mag - (SELECT zp FROM temp.zp_map z
                                        WHERE z.series_key =
                                              cv_lightcurve.series_key)
                   WHERE mag IS NOT NULL
                     AND series_key IN (SELECT series_key FROM temp.zp_map)""")
    con.commit()
    n = con.execute("SELECT count(*) FROM cv_lightcurve "
                    "WHERE cal_mag IS NOT NULL").fetchone()[0]
    total = con.execute("SELECT count(*) FROM cv_lightcurve").fetchone()[0]
    print(f"  calibrated {n:,} of {total:,} light-curve rows "
          f"across {len(ties)} series in {time.time() - t0:.0f}s")
    for skey, zp, cat, band, verdict in ties:
        print(f"    {skey:16s} cal_mag = mag - {zp:.4f}   "
              f"[{cat}:{band}, {verdict}]")
    untied = con.execute("""SELECT series_key FROM cv_series
                            WHERE status='solved' AND series_key NOT IN
                            (SELECT series_key FROM temp.zp_map)
                            ORDER BY 1""").fetchall()
    print(f"  left RELATIVE (cal_mag NULL): {len(untied)} solved series")
    for (s,) in untied:
        print(f"    {s}")
    meta_write(con, {"stage_apply": datetime.now(timezone.utc).isoformat(),
                     "n_calibrated_rows": str(n)})
    con.close()


# ===========================================================================
# Stage: report
# ===========================================================================
def cmd_report(args) -> None:
    """Render the chain-of-evidence page from the database."""
    from macro_phot.report_cattie import render_report
    out = render_report(args.db, Path(args.cache))
    print(f"  wrote {out}")


# ===========================================================================
# Stage: status
# ===========================================================================
def cmd_status(args) -> None:
    """Where the stage stands, and what it currently concludes."""
    con = connect(args.db)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    print(f"database: {args.db}")
    for t in ("cv_cat_fetch", "cv_cat_match", "cv_cat_astrom", "cv_cattie",
              "cv_cattie_star", "cv_cattie_veto", "cv_cattie_trend",
              "cv_cattie_cross"):
        n = con.execute(f"SELECT count(*) FROM {t}").fetchone()[0] \
            if t in have else 0
        print(f"  {t:20s} {n:>10,}")
    if "cv_field_tie" in have:
        bad = con.execute("SELECT count(*) FROM cv_field_tie "
                          "WHERE status!='ok'").fetchone()[0]
        tot = con.execute("SELECT count(*) FROM cv_field_tie").fetchone()[0]
        print(f"  field ties ok:       {tot - bad}/{tot}")
    if "cv_cattie" not in have:
        print("  (no ties yet — run fetch, match, solve)")
        con.close()
        return
    print("\n  per-block verdicts (primary catalogue only):")
    rows = con.execute("""SELECT series_key, catalogue, band, n_fit, zp,
                                 zp_err, colour_term, colour_err, resid_rms,
                                 check_rms_clip, n_check, verdict,
                                 colour_position, n_check_outlier
                          FROM cv_cattie WHERE is_primary=1
                          ORDER BY target_key, era_id, filter""").fetchall()

    def num(x, w=8, nd=4):
        """Fixed-width number, or a dash -- an UNTIED row has NULLs, and a
        status command that crashes on its own worst case is useless
        exactly when it is needed."""
        if x is None or not math.isfinite(float(x)):
            return "-" * w
        return f"{float(x):{w}.{nd}f}"

    for r in rows:
        print(f"    {r[0]:16s} {r[1]:10s} {r[2]:11s} n={r[3] or 0:4d} "
              f"ZP={num(r[4])}+/-{num(r[5], 6)} "
              f"k={num(r[6], 7)}+/-{num(r[7], 6)} "
              f"rms={num(r[8], 6)} check={num(r[9], 6)}"
              f"({r[10] or 0}-{r[13] or 0}) {r[11]:16s} {r[12] or ''}")
    verdicts = [r[11] for r in rows]
    n_blocks = con.execute("SELECT count(*) FROM cv_ref").fetchone()[0]
    n_series = con.execute("SELECT count(*) FROM cv_series "
                           "WHERE status='solved'").fetchone()[0]
    v, deciding = ct.goal_verdict(verdicts, n_series)
    print(f"\n  solved series: {n_series}; (target, era) blocks: {n_blocks}")
    print(f"  CALIBRATION GOAL: {v}")
    print(f"    {deciding}")
    con.close()


# ===========================================================================
# CLI
# ===========================================================================
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--db", type=Path, default=DEFAULT_DB)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--net-timeout", type=float, default=None,
                   help=f"socket timeout in seconds (default "
                        f"{SOCKET_TIMEOUT_S:g}).  Pair it with a low "
                        f"--retries to RECORD a dead archive in a minute "
                        f"rather than wait one out: a stage that cannot "
                        f"write down an outage without first sitting "
                        f"through it is not as resumable as it claims.")
    p.add_argument("--retries", type=int, default=None,
                   help=f"network attempts per query (default "
                        f"{NET_RETRIES}).  Lower it to RECORD an archive "
                        f"outage quickly instead of waiting one out: a "
                        f"recorded failure is evidence, an unattempted "
                        f"catalogue is silence.")
    p.add_argument("--catalogue", choices=sorted(ct.CATALOGUES),
                   help="restrict fetch/match to ONE catalogue (default: "
                        "both; see _catalogues for why this exists)")
    p.add_argument("--force", action="store_true",
                   help="redo work that is already recorded (re-fetch a "
                        "cached cone, re-match a matched block)")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("fieldfix", "fetch", "match", "solve", "validate",
                 "apply", "report", "status", "all"):
        sub.add_parser(name)
    args = p.parse_args()
    if args.net_timeout:
        socket.setdefaulttimeout(float(args.net_timeout))
    table = {"fieldfix": cmd_fieldfix, "fetch": cmd_fetch, "match": cmd_match,
             "solve": cmd_solve, "validate": cmd_validate,
             "apply": cmd_apply, "report": cmd_report, "status": cmd_status}
    if args.cmd == "all":
        for name in ("fieldfix", "fetch", "match", "solve", "validate",
                     "apply", "report"):
            print(f"\n=== {name} ===", flush=True)
            table[name](args)
        cmd_status(args)
    else:
        table[args.cmd](args)


if __name__ == "__main__":
    main()
