#!/usr/bin/env python
"""CV-S7 — the external survey record, and the YZ Cnc branch decision.

WHAT THIS SCRIPT IS FOR
-----------------------
Our CV photometry is DIFFERENTIAL.  A light curve that sits 1.4 magnitudes
above last night's looks exactly like one that does not, until something
outside this telescope supplies an absolute reference.  Two plan tasks turn
on that reference:

``CV-P0-aavso-yzcnc``  THE GATING TASK.  YZ Cnc is an SU UMa dwarf nova.
    Common superhumps occur in SUPEROUTBURST and essentially nowhere else,
    so the strategy's Q3 has two branches -- a superhump period (and maybe
    dP_sh/dt), or an orbital-hump + flickering fallback -- and which one the
    paper takes is a FACT about the accretion state on 2024-02-21 ->
    2024-05-03, not a modelling preference.  Four CV tasks and one whole
    science branch wait on it.

``CV-P0-survey-context``  The long-baseline record for all five targets:
    what the surveys constrain that our nights cannot (state history,
    outburst recurrence), and -- said out loud -- where a survey saturates
    or simply never looked.

SOURCES, AND WHAT EACH ONE COST
-------------------------------
AAVSO   ``https://vsx.aavso.org/index.php?view=api.delim`` -- the AID, by
        AUID, over an explicit JD range.  The richest record by far and the
        only one with coverage inside the decisive window.  Found by reading
        the LCG v2 client's own JavaScript; the documented WebObs HTML route
        ignores its date parameters and would have needed ~500 paged
        requests for one star.
ZTF     IRSA ``nph_light_curves`` cone search.  Public, no credentials.
ASAS-SN The Variable Stars Database ``.json`` export, reached by cone search
        for the star's record UUID.  This is the LEGACY V-band survey; the
        modern g-band Sky Patrol v2 API (``asas-sn.ifa.hawaii.edu``) refuses
        connections from this network, and the failure is RECORDED rather
        than quietly dropped.
ATLAS   Forced photometry needs an account and an API token.  There is none
        here and this script will not create one.  Recorded UNREACHABLE with
        the HTTP status, so the gap is visible instead of being invisible.

THE CIRCULARITY THIS SCRIPT REFUSES TO COMMIT
---------------------------------------------
A large share of the AAVSO rows for YZ Cnc in this window were submitted by
observer ``MALW`` from THIS TELESCOPE ("TAKEN WITH MACRO CONSORTIUM'S ROBERT
L. MUTEL TELESCOPE AT WINER OBSERVATORY").  Classifying our own nights from
our own resubmitted photometry and reporting it as external confirmation
would be circular.  Every row is therefore tagged at parse time
(``macro_phot.external.is_own_observation``), every aggregate is computed
both ways, and each night's verdict records WHICH basis carried it --
independent, bracketed, or own.  The headline verdict is built to survive
deleting every one of our own rows.

STAGES (each resumable, each safe to repeat)
--------------------------------------------
    fetch     pull + cache every source for every target, recording the
              query text, the pull date, the row count and a sha256;
              NEVER re-fetches silently (``--force`` re-pulls)
    classify  measure the quiescent baseline, tag every external night,
              group outburst episodes, and decide the branch
    report    docs/CV_TimeSeries/cv_external_context.html + figures
    status    what is cached, what is classified, the standing verdict
    all       fetch -> classify -> report

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_external.py fetch
    $P pipeline/scripts/run_cv_external.py classify
    $P pipeline/scripts/run_cv_external.py report
    $P pipeline/scripts/run_cv_external.py status

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite)
--------------------------------------------------------------
``cv_external``      THE DELIVERABLE the task names: one row per (target,
                     source) -- n_points, MJD span, bands, cadence, notes.
``cv_ext_fetch``     one row per (target, source): service URL, the exact
                     query, pull date, row count, cache path + sha256, and
                     the failure text when a source was unreachable.
``cv_ext_nightly``   one row per (target, source, UTC night): the median
                     magnitude, how many points, and whether INDEPENDENT
                     observers carried it.
``cv_ext_episode``   outburst episodes with their grade and the reason.
``cv_ext_verdict``   one row per RLMT YZ Cnc night in the window: local
                     night label, UTC date, frame count, state, evidence.
``cv_ext_meta``      build stamps, every constant, and the branch verdict.

Raw responses are cached under ``products/external/<source>/<target>.*.gz``
so a re-run never re-queries a public service, and so the bytes a
classification rests on can be re-read years later.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sqlite3
import ssl
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_phot import external as ex                          # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
DEFAULT_CACHE = REPO_ROOT / "products" / "external"

#: This stage's code version, stamped into cv_ext_meta and read by the
#: provenance DAG to decide whether the page below it has gone stale.
EXTERNAL_CODE_VERSION = "CV-S7 v1.0"

#: How many frames make a run "dense" enough to carry a period analysis.
#: The strategy's dense YZ Cnc blocks are ~100 frames per filter per night
#: cycling three filters; 90 total frames is comfortably below the smallest
#: of them (97) and far above the 1-3 frame snapshot nights, so it separates
#: the two populations without landing on either.
DENSE_RUN_MIN_FRAMES = 90

#: The window the gating task names, as UTC dates.
WINDOW_START_UTC = "2024-02-21"
WINDOW_END_UTC = "2024-05-03"

#: Context padding around that window, in days.  A state ladder needs to see
#: the star's quiescent floor and at least one full outburst to calibrate
#: itself; 80 days each way covers a whole SU UMa supercycle (~130 d for
#: YZ Cnc) so the baseline is measured against real quiescence rather than
#: against whatever the window happened to contain.
#:
#: NOTE this bounds the CLASSIFICATION, not the FETCH.  The two plan tasks
#: want different spans -- the gating one wants a supercycle around the
#: window, the survey-context one wants every year there is -- so the fetch
#: pulls the whole AAVSO record once and the ladder slices the season it
#: needs out of it.  An earlier build bounded the fetch instead and returned
#: 0 AAVSO rows for EU UMa and AN UMa, which would have been published as
#: "AAVSO does not cover these targets" when in fact it holds 875 and 7,615
#: observations of them.
WINDOW_PAD_D = 80

#: The AAVSO pull covers the whole archive.  JD 2400000 is 1858; the upper
#: bound is generous so the query text does not need a clock in it.
AAVSO_JD_MIN = 2400000.0
AAVSO_JD_MAX = 2500000.0

#: The five science targets.  ``auid`` is the AAVSO identifier (the AID key);
#: RA/Dec are the VSX J2000 positions, used for the ZTF and ASAS-SN cone
#: searches.  ``bright_mag`` is each star's bright extreme from VSX, used
#: ONLY to ask whether a survey would have saturated on it -- it is an
#: external catalogue number and the report footer says so.
TARGETS: tuple[dict, ...] = (
    dict(key="stlmi", name="ST LMi", auid="000-BBR-984",
         ra=166.41571, dec=25.10794, vsx_type="AM",
         bright_mag=14.4, faint_mag=18.5),
    dict(key="vvpup", name="VV Pup", auid="000-BBP-388",
         ra=123.77829, dec=-19.05492, vsx_type="AM+ELL",
         bright_mag=13.9, faint_mag=19.6),
    dict(key="euuma", name="EU UMa", auid="000-BBS-438",
         ra=177.48209, dec=28.75203, vsx_type="AM",
         bright_mag=16.45, faint_mag=19.3),
    dict(key="anuma", name="AN UMa", auid="000-BBR-959",
         ra=166.10688, dec=45.05383, vsx_type="AM+E",
         bright_mag=14.5, faint_mag=19.3),
    dict(key="yzcnc", name="YZ Cnc", auid="000-BBP-203",
         ra=122.73604, dec=28.14256, vsx_type="UGSU",
         bright_mag=10.5, faint_mag=16.3),
)
TARGET_BY_KEY = {t["key"]: t for t in TARGETS}

AAVSO_SERVICE = "https://vsx.aavso.org/index.php"
ZTF_SERVICE = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
ASASSN_SERVICE = "https://asas-sn.osu.edu/variables"
ATLAS_SERVICE = "https://fallingstar-data.com/forcedphot/queue/"

#: ZTF cone radius in degrees (~5 arcsec).  Wide enough to catch the
#: object's own oid through proper motion, tight enough to exclude
#: neighbours in these uncrowded high-latitude fields.
ZTF_RADIUS_DEG = 0.0014
#: ASAS-SN lookup radius in ARCMIN.  The Variable Stars Database's ``radius``
#: parameter is arcmin, not degrees: a value of 0.01 that read like a 36
#: arcsec cone is really 0.6 arcsec, and it silently found nothing for three
#: of the five targets while appearing to work on the fourth.  0.5 arcmin is
#: comfortably wider than the 0.17-1.61 arcsec separations these five
#: actually have, and every match is then verified by cross-identification
#: (see macro_phot.external.pick_asassn_match) rather than by proximity.
ASASSN_RADIUS_ARCMIN = 0.5

NET_TIMEOUT_S = 300
NET_RETRIES = 3
USER_AGENT = "MACRO-pipeline/CV-S7 (research; contact james@animal-lamps.com)"


# ===========================================================================
# small I/O helpers
# ===========================================================================

def connect(path: Path, readonly: bool = False) -> sqlite3.Connection:
    """Open the product DB with the long busy timeout the house rule sets.

    An S1 batch solve and an rclone transfer may be running against the same
    spinning disk; 300 s of patience is the difference between waiting and
    failing.
    """
    uri = f"file:{path}?mode=ro" if readonly else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


def _ssl_context() -> ssl.SSLContext:
    """A verifying TLS context, falling back to the certifi bundle.

    This environment's default trust store cannot verify some archive hosts.
    Rather than disable verification globally (which would make every future
    fetch silently unauthenticated), try the platform store first and fall
    back to certifi's bundle only if it is installed.
    """
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return ssl.create_default_context()


def http_get(url: str, retries: int = NET_RETRIES) -> bytes:
    """GET with linear backoff.  Raises the LAST exception on exhaustion.

    Public archives rate-limit and occasionally drop connections; a single
    attempt would turn a transient 503 into a permanent 'source
    unreachable' line in a published report.
    """
    ctx = _ssl_context()
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=NET_TIMEOUT_S,
                                        context=ctx) as r:
                return r.read()
        except Exception as e:                       # noqa: BLE001
            last = e
            if attempt < retries:
                time.sleep(2.0 * attempt)
    raise last                                        # type: ignore[misc]


def atomic_write_gz(path: Path, data: bytes) -> str:
    """Write gzipped bytes atomically; return the sha256 of the RAW bytes.

    The hash is of the uncompressed payload so it does not move when gzip's
    implementation or mtime does -- it identifies what the archive SAID, not
    how we stored it.  Temp-then-rename means a killed process leaves either
    the old cache or the new one, never a half file that would parse into a
    truncated light curve.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, path)
    return hashlib.sha256(data).hexdigest()


def read_cache(path: Path) -> bytes:
    """Read a gzipped cache file back to its raw bytes."""
    with gzip.open(path, "rb") as fh:
        return fh.read()


def git_commit() -> str:
    """Current commit, or '' when git is unavailable or the tree is not one."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO_ROOT, capture_output=True, text=True,
                             timeout=10)
        return out.stdout.strip() if out.returncode == 0 else ""
    except Exception:                                  # noqa: BLE001
        return ""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


# ===========================================================================
# Query builders — one per source.  Kept pure and separate so the exact
# query text stored in cv_ext_fetch is the SAME string that was sent.
# ===========================================================================

def aavso_url(auid: str, jd_from: float, jd_to: float) -> str:
    """The AAVSO AID delimited-data query for one star over one JD range."""
    return (f"{AAVSO_SERVICE}?" + urllib.parse.urlencode({
        "view": "api.delim", "ident": auid,
        "fromjd": f"{jd_from:.5f}", "tojd": f"{jd_to:.5f}",
        "delimiter": "@@@"}))


def ztf_url(ra: float, dec: float, radius: float = ZTF_RADIUS_DEG) -> str:
    """The IRSA ZTF light-curve cone search, all bands, CSV."""
    return (f"{ZTF_SERVICE}?" + urllib.parse.urlencode({
        "POS": f"CIRCLE {ra} {dec} {radius}", "FORMAT": "csv"},
        quote_via=urllib.parse.quote))


def asassn_search_url(ra: float, dec: float,
                      radius_arcmin: float = ASASSN_RADIUS_ARCMIN) -> str:
    """Cone search of the ASAS-SN Variable Stars Database (HTML).

    ``radius`` is in ARCMIN — see the note on ASASSN_RADIUS_ARCMIN.
    """
    return (f"{ASASSN_SERVICE}?" + urllib.parse.urlencode({
        "ra": f"{ra}", "dec": f"{dec}", "radius": f"{radius_arcmin}",
        "commit": "Search"}))


def asassn_json_url(uuid: str) -> str:
    """The V-band light-curve export for one ASAS-SN variable record."""
    return f"{ASASSN_SERVICE}/{uuid}.json"


# ===========================================================================
# Stage: fetch
# ===========================================================================

FETCH_DDL = """CREATE TABLE IF NOT EXISTS cv_ext_fetch (
    target TEXT, source TEXT, service TEXT, query TEXT,
    pulled_utc TEXT, n_rows INTEGER, cache_path TEXT, cache_sha256 TEXT,
    ok INTEGER, note TEXT,
    PRIMARY KEY (target, source))"""


def _record_fetch(con, target, source, service, query, n_rows, cache,
                  sha, ok, note) -> None:
    con.execute("INSERT OR REPLACE INTO cv_ext_fetch VALUES (?,?,?,?,?,?,?,?,?,?)",
                (target, source, service, query, now_utc(), n_rows,
                 cache, sha, 1 if ok else 0, note))
    con.commit()


def _already(con, target, source, cache_root: Path) -> bool:
    """Is this (target, source) already cached AND still on disk?

    Both halves matter: a DB row without its file is a promise the cache
    cannot keep, and re-fetching is then correct.
    """
    row = con.execute("SELECT cache_path, ok FROM cv_ext_fetch "
                      "WHERE target=? AND source=?", (target, source)).fetchone()
    if not row:
        return False
    path, ok = row
    if not ok:
        return False           # a recorded FAILURE is always worth retrying
    return bool(path) and (REPO_ROOT / path).exists()


def cmd_fetch(args) -> None:
    """Pull every source for every target, once, and cache the raw bytes."""
    con = connect(Path(args.db))
    con.execute(FETCH_DDL)
    con.commit()
    cache_root = Path(args.cache)
    sources = args.sources.split(",") if args.sources else \
        ["aavso", "ztf", "asassn", "atlas"]

    jd_from, jd_to = AAVSO_JD_MIN, AAVSO_JD_MAX

    for t in TARGETS:
        for source in sources:
            if _already(con, t["key"], source, cache_root) and not args.force:
                print(f"  {source:8s} {t['name']:8s}: cached — skipping")
                continue
            print(f"  {source:8s} {t['name']:8s}: pulling", flush=True)
            try:
                if source == "aavso":
                    _fetch_aavso(con, t, cache_root, jd_from, jd_to)
                elif source == "ztf":
                    _fetch_ztf(con, t, cache_root)
                elif source == "asassn":
                    _fetch_asassn(con, t, cache_root)
                elif source == "atlas":
                    _fetch_atlas(con, t)
                else:
                    raise ValueError(f"unknown source {source!r}")
            except Exception as e:                     # noqa: BLE001
                print(f"      FAILED: {type(e).__name__}: {e}", flush=True)
                _record_fetch(con, t["key"], source, "", "", 0, "", "",
                              False, f"FETCH FAILED: {type(e).__name__}: {e}")
    con.close()


def _jd_of(utc_date: str) -> float:
    """JD at 00:00 UTC of a ``YYYY-MM-DD`` date."""
    d = datetime.strptime(utc_date, "%Y-%m-%d").replace(tzinfo=timezone.utc)
    return (d - datetime(1858, 11, 17, tzinfo=timezone.utc)).total_seconds() \
        / 86400.0 + ex.JD_MINUS_MJD


def _fetch_aavso(con, t, cache_root: Path, jd_from: float,
                 jd_to: float) -> None:
    """AAVSO AID for one target — the WHOLE record.

    The survey-context task wants every year AAVSO holds (these polars have
    40-year records), and the gating task wants one supercycle around the
    2024 window.  Pulling everything once and slicing at classification time
    serves both from a single cached response, and costs 22 MB and four
    seconds for the largest star.
    """
    url = aavso_url(t["auid"], jd_from, jd_to)
    raw = http_get(url)
    obs = ex.parse_aavso_delim(raw.decode("utf-8", "replace"))
    cache = cache_root / "aavso" / f"{t['key']}.delim.gz"
    sha = atomic_write_gz(cache, raw)
    n_own = sum(1 for o in obs if o.is_own)
    note = (f"JD {jd_from:.1f}-{jd_to:.1f}; {n_own} of {len(obs)} rows are "
            f"RLMT data resubmitted to AAVSO (observer code or comment) and "
            f"are excluded from every independent aggregate")
    _record_fetch(con, t["key"], "aavso", AAVSO_SERVICE, url, len(obs),
                  str(cache.relative_to(REPO_ROOT)), sha, True, note)
    print(f"      {len(obs)} rows ({n_own} ours) -> {cache.name}", flush=True)


def _fetch_ztf(con, t, cache_root: Path) -> None:
    """ZTF DR light curves from IRSA for one target."""
    url = ztf_url(t["ra"], t["dec"])
    raw = http_get(url)
    obs = ex.parse_ztf_csv(raw.decode("utf-8", "replace"))
    cache = cache_root / "ztf" / f"{t['key']}.csv.gz"
    sha = atomic_write_gz(cache, raw)
    at_risk, why = ex.saturation_risk("ztf", t["bright_mag"])
    note = (f"cone r={ZTF_RADIUS_DEG} deg; "
            + ("BRIGHT-END WARNING: " if at_risk else "") + why)
    _record_fetch(con, t["key"], "ztf", ZTF_SERVICE, url, len(obs),
                  str(cache.relative_to(REPO_ROOT)), sha, True, note)
    print(f"      {len(obs)} rows -> {cache.name}", flush=True)


def _fetch_asassn(con, t, cache_root: Path) -> None:
    """ASAS-SN legacy V-band record, via a cone search for its UUID.

    Two requests: the cone search to learn the record id, then the JSON
    export.  Both URLs are stored so the path can be re-walked.
    """
    search_url = asassn_search_url(t["ra"], t["dec"])
    html = http_get(search_url).decode("utf-8", "replace")
    match, why = ex.pick_asassn_match(ex.parse_asassn_search(html), t["name"])
    if match is None:
        _record_fetch(con, t["key"], "asassn", ASASSN_SERVICE, search_url, 0,
                      "", "", False,
                      f"no usable ASAS-SN Variable Stars Database record "
                      f"within {ASASSN_RADIUS_ARCMIN} arcmin of the VSX "
                      f"position: {why}")
        print(f"      no ASAS-SN record found ({why})", flush=True)
        return
    uuid = match["uuid"]
    url = asassn_json_url(uuid)
    raw = http_get(url)
    obs = ex.parse_asassn_json(json.loads(raw.decode("utf-8", "replace")))
    cache = cache_root / "asassn" / f"{t['key']}.json.gz"
    sha = atomic_write_gz(cache, raw)
    note = ("legacy V-band Variable Stars Database export; the modern "
            "g-band Sky Patrol v2 API (asas-sn.ifa.hawaii.edu) refuses "
            "connections from this network, so post-2018 ASAS-SN coverage "
            f"is NOT included here. Match: {why}. "
            f"ASAS-SN name {match['asassn_name']}, type {match['var_type']}, "
            f"mean V {match['mean_vmag']}, amplitude {match['amplitude']}. "
            f"record uuid={uuid}")
    _record_fetch(con, t["key"], "asassn", ASASSN_SERVICE,
                  f"{search_url}\n{url}", len(obs),
                  str(cache.relative_to(REPO_ROOT)), sha, True, note)
    print(f"      {len(obs)} rows -> {cache.name}", flush=True)


def _fetch_atlas(con, t) -> None:
    """Record ATLAS forced photometry as unreachable, with the evidence.

    The service requires a registered account and an API token.  This
    script does not create accounts and does not hold credentials, so the
    correct output is a recorded, dated NEGATIVE -- the HTTP status the
    server actually returned -- rather than an empty row that a later
    reader would mistake for 'ATLAS saw nothing'.
    """
    status = "no response"
    try:
        http_get(ATLAS_SERVICE, retries=1)
        status = "HTTP 200 (unexpected — endpoint normally requires a token)"
    except Exception as e:                             # noqa: BLE001
        status = f"{type(e).__name__}: {e}"
    _record_fetch(con, t["key"], "atlas", ATLAS_SERVICE, ATLAS_SERVICE, 0,
                  "", "", False,
                  f"UNREACHABLE — ATLAS forced photometry requires a "
                  f"registered account and API token; none is configured in "
                  f"this environment and this script does not create one. "
                  f"Probe result: {status}")
    print(f"      ATLAS unreachable ({status})", flush=True)


# ===========================================================================
# Stage: classify
# ===========================================================================

CLASSIFY_DDL = [
    """CREATE TABLE IF NOT EXISTS cv_external (
        target TEXT, source TEXT, n_points INTEGER,
        mjd_min REAL, mjd_max REAL, span_d REAL,
        bands TEXT, n_nights INTEGER, median_gap_d REAL,
        n_independent INTEGER, notes TEXT,
        PRIMARY KEY (target, source))""",
    """CREATE TABLE IF NOT EXISTS cv_ext_nightly (
        target TEXT, source TEXT, utc_night TEXT, mag REAL, n INTEGER,
        n_independent INTEGER, bands TEXT, observers TEXT, spread REAL,
        independent INTEGER, state TEXT, amp REAL,
        PRIMARY KEY (target, source, utc_night))""",
    """CREATE TABLE IF NOT EXISTS cv_ext_episode (
        target TEXT, source TEXT, start_night TEXT, end_night TEXT,
        peak_night TEXT, peak_amp REAL, duration_d REAL,
        plateau_start TEXT, plateau_end TEXT, plateau_d REAL,
        n_nights INTEGER, kind TEXT, why TEXT,
        PRIMARY KEY (target, source, start_night))""",
    # ``mag`` is the magnitude that CARRIED the verdict — the independent
    # nightly median where outside observers covered the night, our own
    # V-equivalent where they did not.  ``own_mag`` is always OUR value,
    # recorded separately and unconditionally, because the two are the same
    # number on independent nights and a figure that plotted ``mag`` against
    # the independent median was therefore plotting a quantity against
    # itself and drawing a perfect 1:1 line with zero scatter.
    """CREATE TABLE IF NOT EXISTS cv_ext_verdict (
        target TEXT, local_night TEXT, utc_night TEXT, n_frames INTEGER,
        filters TEXT, is_dense INTEGER, state TEXT, mag REAL, amp REAL,
        basis TEXT, episode TEXT, evidence TEXT, own_mag REAL,
        PRIMARY KEY (target, local_night))""",
    """CREATE TABLE IF NOT EXISTS cv_ext_meta (key TEXT PRIMARY KEY,
        value TEXT)""",
]


def load_obs(con, target: str, source: str) -> list[ex.Obs]:
    """Re-parse a cached response into normalised observations.

    Reading from the cache rather than the network is what makes
    ``classify`` cheap, repeatable and offline: the classification can be
    re-derived a hundred times without touching a public archive once.
    """
    row = con.execute("SELECT cache_path, ok FROM cv_ext_fetch "
                      "WHERE target=? AND source=?", (target, source)).fetchone()
    if not row or not row[1] or not row[0]:
        return []
    path = REPO_ROOT / row[0]
    if not path.exists():
        return []
    raw = read_cache(path)
    if source == "aavso":
        return ex.parse_aavso_delim(raw.decode("utf-8", "replace"))
    if source == "ztf":
        return ex.parse_ztf_csv(raw.decode("utf-8", "replace"))
    if source == "asassn":
        return ex.parse_asassn_json(json.loads(raw.decode("utf-8", "replace")))
    return []


def our_nights(con, target: str) -> list[tuple[str, int, str]]:
    """Our observing nights for a target: (local night label, frames, filters).

    Read from ``cv_frames``, which is the canonical per-target frame view.
    """
    return [(n, c, f) for n, c, f in con.execute(
        """SELECT night, count(*), group_concat(DISTINCT filter)
           FROM cv_frames WHERE target_key=? GROUP BY night ORDER BY night""",
        (target,))]


def cmd_classify(args) -> None:
    """Measure the ladder, tag every night, grade every episode, decide."""
    con = connect(Path(args.db))
    # Every table below is DERIVED, in full, from the cached raw responses.
    # Dropping them first makes `classify` idempotent rather than
    # accumulative: a re-run after a threshold changes must not leave last
    # run's episodes sitting beside this run's, and INSERT OR REPLACE alone
    # would do exactly that whenever a boundary moved.  The cache — the only
    # thing here that cost a network call — is never touched.
    for name in ("cv_external", "cv_ext_nightly", "cv_ext_episode",
                 "cv_ext_verdict", "cv_ext_meta"):
        con.execute(f"DROP TABLE IF EXISTS {name}")
    for ddl in CLASSIFY_DDL:
        con.execute(ddl)
    con.commit()
    meta: dict[str, str] = {}

    # ---- 1. coverage + nightly points, per (target, source) --------------
    for t in TARGETS:
        for source in ("aavso", "ztf", "asassn"):
            obs = load_obs(con, t["key"], source)
            note = con.execute("SELECT note FROM cv_ext_fetch WHERE target=? "
                               "AND source=?", (t["key"], source)).fetchone()
            cov = ex.coverage(t["key"], source, obs, note[0] if note else "")
            con.execute(
                "INSERT OR REPLACE INTO cv_external VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                (cov.target, cov.source, cov.n_points, cov.mjd_min,
                 cov.mjd_max, cov.span_d, ",".join(cov.bands), cov.n_nights,
                 cov.median_gap_d,
                 sum(1 for o in obs if not o.is_own), cov.notes))
            if not obs:
                continue
            # Nightly points on the V-like ladder.  ZTF/ASAS-SN carry their
            # own bands and are bucketed whole; AAVSO is restricted to the
            # V-like set so the ladder is not polluted by Sloan zero points.
            bands = ex.V_LIKE_BANDS if source == "aavso" else None
            pts = ex.nightly(obs, bands=bands)
            base = ex.quiescent_baseline(pts)
            for p in pts:
                state = (ex.classify_night(p.mag, base) if base is not None
                         else ex.STATE_UNKNOWN)
                con.execute(
                    "INSERT OR REPLACE INTO cv_ext_nightly VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?)",
                    (t["key"], source, p.night, p.mag, p.n, p.n_independent,
                     ",".join(p.bands), ",".join(p.observers[:8]), p.spread,
                     1 if p.independent else 0, state,
                     (base - p.mag) if base is not None else None))
            if base is not None:
                meta[f"baseline_{t['key']}_{source}"] = f"{base:.4f}"
                for ep in ex.find_episodes(pts, base):
                    con.execute(
                        "INSERT OR REPLACE INTO cv_ext_episode VALUES "
                        "(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        (t["key"], source, ep.start, ep.end, ep.peak_night,
                         ep.peak_amp, ep.duration_d, ep.plateau_start,
                         ep.plateau_end, ep.plateau_d, ep.n_nights, ep.kind,
                         ep.why))
        con.commit()

    # ---- 2. the gating question: YZ Cnc, night by night -------------------
    branch, reasoning = _classify_yzcnc(con, meta)
    meta.update({
        "built_utc": now_utc(),
        "external_code_version": EXTERNAL_CODE_VERSION,
        "git_commit": git_commit(),
        "window_start_utc": WINDOW_START_UTC,
        "window_end_utc": WINDOW_END_UTC,
        "window_pad_d": str(WINDOW_PAD_D),
        "dense_run_min_frames": str(DENSE_RUN_MIN_FRAMES),
        "amp_quiescent_max": str(ex.AMP_QUIESCENT_MAX),
        "amp_outburst_min": str(ex.AMP_OUTBURST_MIN),
        "superoutburst_amp_min": str(ex.SUPEROUTBURST_AMP_MIN),
        "superoutburst_days_min": str(ex.SUPEROUTBURST_DAYS_MIN),
        "quiescent_decile": str(ex.QUIESCENT_DECILE),
        "branch": branch,
        "branch_reasoning": reasoning,
    })
    con.executemany("INSERT OR REPLACE INTO cv_ext_meta VALUES (?,?)",
                    sorted(meta.items()))
    con.commit()
    print(f"\nBRANCH: {branch}\n{reasoning}")
    con.close()


def _classify_yzcnc(con, meta: dict) -> tuple[str, str]:
    """Tag every RLMT YZ Cnc night in the window, then pick the branch."""
    all_obs = load_obs(con, "yzcnc", "aavso")
    if not all_obs:
        return "UNDECIDED", "no AAVSO record cached for YZ Cnc"

    # Slice ONE SEASON out of the 30-year record.  A quiescent baseline
    # measured across three decades would average over long-term changes in
    # the star and in the observers' equipment; the ladder that tags these
    # nights has to be the ladder this season stood on.
    jd_lo = _jd_of(WINDOW_START_UTC) - WINDOW_PAD_D
    jd_hi = _jd_of(WINDOW_END_UTC) + WINDOW_PAD_D
    obs = [o for o in all_obs if jd_lo <= o.jd <= jd_hi]
    meta["yzcnc_season_jd_range"] = f"{jd_lo:.1f}-{jd_hi:.1f}"
    meta["yzcnc_season_rows"] = str(len(obs))

    # The ladder is built from INDEPENDENT V-like points only.  That is the
    # whole methodological point: our own resubmitted rows may inform a
    # night's tag, but they may not define the scale that tag is measured
    # against, or the argument becomes circular.
    ind_pts = ex.nightly(obs, bands=ex.V_LIKE_BANDS, independent_only=True)
    base = ex.quiescent_baseline(ind_pts)
    if base is None:
        return "UNDECIDED", "too few independent AAVSO nights to set a baseline"
    meta["baseline_yzcnc_independent"] = f"{base:.4f}"
    episodes = ex.find_episodes(ind_pts, base)

    # Our own AAVSO-submitted photometry, placed on the same ladder by a
    # MEASURED band offset rather than by assuming SR == V.
    own = [o for o in obs if o.is_own]
    offsets: dict[str, float] = {}
    for band in sorted({o.band for o in own}):
        off, n, scat = ex.band_offset(obs, band)
        if off is not None and n >= 3:
            offsets[band] = off
            meta[f"offset_{band}_minus_V"] = f"{off:+.3f} (n={n} nights, " \
                                             f"scatter {scat:.3f} mag)"
    own_pts: dict[str, float] = {}
    if offsets:
        buckets: dict[str, list[float]] = {}
        for o in own:
            if o.band in offsets:
                buckets.setdefault(o.utc_night, []).append(o.mag - offsets[o.band])
        own_pts = {n: ex._median(v) for n, v in buckets.items()}

    ind_by_night = {p.night: p for p in ind_pts}
    rows = []
    for local_night, n_frames, filters in our_nights(con, "yzcnc"):
        utc = ex.night_label_to_utc_date(local_night)
        if not (WINDOW_START_UTC <= utc <= WINDOW_END_UTC):
            continue
        ip = ind_by_night.get(utc)
        episode = ""
        for ep in episodes:
            if ep.start <= utc <= ep.end:
                episode = ep.kind
                break
        if ip is not None:
            state = ex.classify_night(ip.mag, base)
            mag, amp, basis = ip.mag, base - ip.mag, "independent"
            evidence = (f"{ip.n_independent} independent point(s) from "
                        f"{', '.join(ip.observers[:4])} in {','.join(ip.bands)}"
                        f"; median {ip.mag:.2f}")
        elif utc in own_pts:
            mag = own_pts[utc]
            state = ex.classify_night(mag, base)
            amp = base - mag
            basis = "own"
            before, after = ex.bracket(ind_pts, utc)
            excluded, why = ex.superoutburst_excluded_by_brackets(
                before, after, base)
            evidence = (f"no independent observer covered this night; RLMT "
                        f"photometry resubmitted to AAVSO gives V-equivalent "
                        f"{mag:.2f} via measured band offsets. "
                        + ("Superoutburst independently excluded: " + why
                           if excluded else "Bracket test: " + why))
        else:
            before, after = ex.bracket(ind_pts, utc)
            excluded, why = ex.superoutburst_excluded_by_brackets(
                before, after, base)
            state, mag, amp = ex.STATE_UNKNOWN, None, None
            basis = "bracketed" if excluded else "none"
            evidence = why
        rows.append((("yzcnc"), local_night, utc, n_frames, filters,
                     1 if n_frames >= DENSE_RUN_MIN_FRAMES else 0,
                     state, mag, amp, basis, episode, evidence,
                     own_pts.get(utc)))
    con.executemany("INSERT OR REPLACE INTO cv_ext_verdict VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)
    con.commit()

    verdicts = [ex.NightVerdict(utc_night=r[2], local_night=r[1],
                                n_frames=r[3], state=r[6], mag=r[7],
                                amp=r[8], basis=r[9], evidence=r[11],
                                episode=r[10]) for r in rows]
    return ex.branch_recommendation(verdicts, DENSE_RUN_MIN_FRAMES)


# ===========================================================================
# Stage: report / status
# ===========================================================================

def cmd_report(args) -> None:
    from macro_phot import report_external as rp
    path = rp.render_report(Path(args.db))
    print(f"wrote {path}")


def cmd_status(args) -> None:
    """What is cached, what is classified, and the standing verdict."""
    con = connect(Path(args.db), readonly=True)
    have = {r[0] for r in con.execute(
        "SELECT name FROM sqlite_master WHERE type='table'")}
    if "cv_ext_fetch" not in have:
        print("CV-S7: never run (no cv_ext_fetch table)")
        return
    print(f"{'target':8s} {'source':8s} {'rows':>7s} {'ok':>3s}  note")
    for r in con.execute("SELECT target, source, n_rows, ok, note FROM "
                         "cv_ext_fetch ORDER BY target, source"):
        print(f"{r[0]:8s} {r[1]:8s} {r[2]:7d} {r[3]:3d}  {(r[4] or '')[:80]}")
    if "cv_ext_meta" in have:
        meta = dict(con.execute("SELECT key, value FROM cv_ext_meta"))
        print(f"\nbuilt   {meta.get('built_utc','-')}")
        print(f"version {meta.get('external_code_version','-')}")
        print(f"BRANCH  {meta.get('branch','-')}")
        print(f"        {meta.get('branch_reasoning','-')}")
    if "cv_ext_verdict" in have:
        print("\nRLMT YZ Cnc nights in window:")
        # printf('%.2f', NULL) is '0.00' in SQLite, not NULL, so coalesce
        # cannot catch a missing magnitude — a night nobody measured would
        # print as V=0.00, the brightest object in the sky.
        for r in con.execute(
                "SELECT utc_night, n_frames, is_dense, state, "
                "CASE WHEN mag IS NULL THEN '-' ELSE printf('%.2f',mag) END, "
                "basis, episode FROM cv_ext_verdict ORDER BY utc_night"):
            dense = "DENSE" if r[2] else "     "
            print(f"  {r[0]}  {r[1]:4d} fr {dense}  {r[3]:10s} "
                  f"V={r[4]:>6s}  [{r[5]}] {r[6]}")
    con.close()


def cmd_all(args) -> None:
    cmd_fetch(args)
    cmd_classify(args)
    cmd_report(args)


# ===========================================================================
# CLI
# ===========================================================================

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_cv_external.py",
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", default=str(DEFAULT_DB),
                   help="CV products database (default: %(default)s)")
    p.add_argument("--cache", default=str(DEFAULT_CACHE),
                   help="raw-response cache root (default: %(default)s)")
    sub = p.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="pull + cache every source, once")
    f.add_argument("--force", action="store_true",
                   help="re-pull even when a cache entry exists")
    f.add_argument("--sources", default="",
                   help="comma list (aavso,ztf,asassn,atlas); default all")
    f.set_defaults(func=cmd_fetch)

    c = sub.add_parser("classify", help="measure the ladder and decide")
    c.set_defaults(func=cmd_classify)

    r = sub.add_parser("report", help="render the evidence page")
    r.set_defaults(func=cmd_report)

    s = sub.add_parser("status", help="what is cached and what it says")
    s.set_defaults(func=cmd_status)

    a = sub.add_parser("all", help="fetch -> classify -> report")
    a.add_argument("--force", action="store_true")
    a.add_argument("--sources", default="")
    a.set_defaults(func=cmd_all)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
