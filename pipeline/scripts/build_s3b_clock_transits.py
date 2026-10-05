#!/usr/bin/env python
"""Build S3b: an absolute clock check per camera era from archived transits.

WHAT THIS SCRIPT DOES (foundation task F-8 of the 2026-10-03 plan review)
-------------------------------------------------------------------------
S3 established what every header time stamp MEANS and converted it to
BJD_TDB.  It could not say whether the observatory clock was RIGHT: its one
astrophysical check (AG LMi against a survey ephemeris) is uncertain by more
than an hour.  This stage closes that gap with events whose times the
literature predicts to seconds — exoplanet transits and the eclipses of
post-common-envelope binaries that students happened to observe — and
measures (our clock) - (true time) for every camera era that has one.

Stages (``--stage``; default ``all`` runs them in order):

``ephemerides``  load the cached literature ephemerides
                 (``products/clock/ephemerides.json``); with
                 ``--refresh-ephemerides`` re-fetch them first (ExoClock
                 catalogue + one NASA Exoplanet Archive cross-check per
                 planet) and rewrite the cache.
``census``       find, in the manifest (opened READ-ONLY), every time
                 series that points at a clock target; predict the events
                 inside each; admit the ones that are fully covered and
                 predicted well enough.  Every refusal is recorded with its
                 reason.
``photometry``   differential aperture photometry of the admitted series
                 (resumable; one row per frame and per star).
``fit``          build the light curves, fit each event's mid-time blind,
                 bootstrap its error, re-fit under every analysis variant,
                 and run the injection-recovery test.
``summary``      O - C per event, per target, per clock era; the verdicts;
                 the StackPro and raw-vs-reduced pair comparisons; the
                 explanation of S3's -294 s AG LMi residual; what follows
                 for the shared ~1,065 s ST LMi / EU UMa offset.
``report``       figures + ``docs/pipeline/s3b_clock.html``, rendered from
                 the products database alone.

OUTPUTS
-------
* ``products/clock/clock_transits.sqlite``  — every table (``s3b_*``)
* ``products/clock/ephemerides.json``       — the literature inputs, cached
* ``docs/pipeline/s3b_clock.html`` and ``docs/pipeline/figures/s3b/*.png``

READ-ONLY DISCIPLINE
--------------------
The archive, the manifest and the CV products database are opened
read-only.  Nothing outside ``products/clock/`` and the two ``docs`` paths
above is written.

USAGE
-----
    /opt/miniconda3/envs/rlmt-checks/bin/python \\
        pipeline/scripts/build_s3b_clock_transits.py              # all
    ... build_s3b_clock_transits.py --stage fit --workers 8        # one

Every number on the report page is the result of a query against the
products database; none is typed.  The only typed numbers in this file are
INPUTS — literature ephemerides and nominal eclipse shapes — each with its
citation beside it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import sqlite3
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import clock_transits as ct                  # noqa: E402
from macro_core import timing as tm                          # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_CV_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
PRODUCT_DIR = REPO_ROOT / "products" / "clock"
DEFAULT_DB = PRODUCT_DIR / "clock_transits.sqlite"
EPHEMERIS_CACHE = PRODUCT_DIR / "ephemerides.json"
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s3b"
HTML_PATH = DOCS_DIR / "s3b_clock.html"

EXOCLOCK_URL = "https://www.exoclock.space/database/planets_json"
NASA_TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

#: Same definition of "canonical science frame" S3 uses (see
#: build_s3_timing.SCIENCE_WHERE); calibration paths are excluded by S0b's
#: classification, not by IMAGETYP alone.
SCIENCE_WHERE = ("is_canonical = 1 AND "
                 "(imagetyp IS NULL OR imagetyp LIKE 'Light%') AND "
                 "path NOT IN (SELECT path FROM calib_frames)")

#: A series belongs to a clock target when its mean frame centre lies
#: within this radius of the target (every camera's half-field exceeds it).
CONE_RADIUS_ARCMIN = 12.0

#: Dispersed (grism) filters are not photometry.
GRISM_FILTERS = ("hrg", "lrg", "HaGrism", "OGGrism", "6", "W")

#: Quadratic limb-darkening coefficients held fixed in the fits.  A
#: mid-time is insensitive to them by symmetry; the ``ld_low`` / ``ld_high``
#: variants (u1 -+ LD_VARIANT_STEP) measure that insensitivity per event.
LD_PLANET = (0.45, 0.25)        # K/G dwarf host, red optical
LD_VARIANT_STEP = 0.20

#: Source detection for registration (bright stars only).
DETECT_SIGMA = 8.0
DETECT_MINAREA = 9
#: A comparison candidate's reference-frame peak must stay below this
#: fraction of the saturation veto (head-room for seeing changes).
COMP_PEAK_HEADROOM = 0.8
#: Comparison candidates stay this many pixels (plus the annulus) off the
#: frame edge.
EDGE_MARGIN_PX = 40.0
#: The target must be found within this many arcseconds of its catalogue
#: position on the reference frame.
TARGET_FIND_ARCSEC = 20.0
#: How many reference-frame candidates are tried before a series is
#: declared unidentifiable.
MAX_REF_TRIES = 6

# ---------------------------------------------------------------------------
# Eclipsing-binary clock targets (literature inputs, each cited)
# ---------------------------------------------------------------------------
#: Times are BJD_TDB.  ``primary`` is the ephemeris the O - C is quoted
#: against; ``alt`` is an INDEPENDENT one used only to measure how far two
#: published predictions disagree at our epoch.  ``sys_floor_s`` is the
#: published amplitude of the star's own eclipse-timing variations about a
#: linear ephemeris — an astrophysical term no clock test can remove, so
#: it is carried as a systematic on the target's O - C.
#:
#: ``shape`` is only the STARTING template (duration and depth are refitted
#: per event; the mid-time of a symmetric eclipse does not depend on it).
EB_TARGETS = {
    "NSVS 07826147": {
        "aliases": "DD CrB, FBS 1531+381, 2M1533+3759",
        "kind": "eclipse",
        "class": "sdB + dM (HW Vir type), primary eclipse",
        # Sesame/SIMBAD (Gaia DR3), ICRS, fetched 2026-10-03.
        "ra_deg": 233.45601721, "dec_deg": 37.99113761,
        "primary": {
            "t0": 2455611.9265712, "sig_t0": 1.7e-6,
            "period": 0.16177044594, "sig_p": 5e-11,
            "quad": 0.0, "sig_quad": 0.0,
            "ref": ("Basturk et al. 2026, MNRAS (accepted), "
                    "arXiv:2602.09925 — linear ephemeris from 618 primary "
                    "minima through 2025"),
        },
        "alt": {
            "t0": 2455611.926576, "sig_t0": 5e-6,
            "period": 0.161770447, "sig_p": 1e-9,
            "quad": -0.9e-13, "sig_quad": 0.3e-13,
            "ref": ("Pulley et al. 2025, MNRAS 544, 24 "
                    "(arXiv:2507.06748) — quadratic ephemeris, data "
                    "through June 2025"),
        },
        # Pulley et al. 2025: residuals about their ephemeris span about
        # +-15 s through June 2025.
        "sys_floor_s": 15.0,
        # Nominal geometry after For et al. 2010 (ApJ 708, 253):
        # a = 0.98 Rsun, R1 = 0.166, R2 = 0.152, i = 86.6 deg.
        "shape": {"a_over_r": 5.90, "k": 0.916, "b": 0.350,
                  "u1": 0.25, "u2": 0.20, "f1": 0.9,
                  "free": ("a_over_r", "f1"), "baseline_order": 2},
        "min_depth": 0.05,
    },
    "GK Vir": {
        "aliases": "PG 1413+015",
        "kind": "eclipse",
        "class": "white dwarf + dM, total eclipse of the white dwarf",
        "ra_deg": 213.90172013, "dec_deg": 1.28839633,
        "primary": {
            # BMJD(TDB) 42543.33745(3) + 0.3443308470(8) E
            "t0": 2442543.83745, "sig_t0": 3e-5,
            "period": 0.3443308470, "sig_p": 8e-10,
            "quad": 0.0, "sig_quad": 0.0,
            "ref": ("Yates et al. 2026, MNRAS 547, stag358 "
                    "(arXiv:2602.17800), Table 1 — best-fit linear "
                    "ephemeris of the long-term timing programme"),
        },
        "alt": {
            # MJD(BTDB) 42543.3377143(30) + 0.344330838759(92) E
            "t0": 2442543.8377143, "sig_t0": 3.0e-6,
            "period": 0.344330838759, "sig_p": 9.2e-11,
            "quad": 0.0, "sig_quad": 0.0,
            "ref": ("Parsons et al. 2010, MNRAS 407, 2362 "
                    "(arXiv:1005.3958)"),
        },
        # No published amplitude adopted: the disagreement of the two
        # ephemerides at our epoch (computed by the build) IS the
        # systematic for this star.
        "sys_floor_s": 0.0,
        # Nominal geometry after Parsons et al. 2012 (MNRAS 420, 3281):
        # R_WD = 0.0170, R_2 = 0.155 Rsun, a = 1.82 Rsun, i ~ 89.5 deg.
        "shape": {"a_over_r": 107.0, "k": 9.1, "b": 0.0,
                  "u1": 0.30, "u2": 0.20, "f1": 0.8,
                  "free": ("a_over_r", "f1"), "baseline_order": 2},
        "min_depth": 0.05,
    },
}

#: Clock candidates that were examined and REFUSED before any photometry,
#: with the reason — printed on the report so the selection is auditable.
REFUSED_TARGETS = (
    ("V2301 Oph", "eclipsing polar; 2025-04/05 (clock era D)",
     "published ephemerides disagree by minutes at 2025: the quadratic "
     "term of Ramsay & Cropper 2007 (MNRAS 379, 1209), 3.18(62)e-13 d "
     "per cycle^2, alone is uncertain by ~140 s after 162,000 cycles; "
     "frames are 240 s exposures across a ~5 min eclipse"),
    ("NN Ser", "WD + dM; 2026-04-01 (clock era F)",
     "two recent linear ephemerides (Ozdonmez et al. 2023; Yates et al. "
     "2026) differ by ~270 s at 2026 — the star's own period variation "
     "exceeds the 120 s criterion; frames are 180 s exposures"),
    ("QS Vir", "WD + dM; 2023-04, 2026-04/05",
     "eclipse-timing variations of order +-100 s about any linear "
     "ephemeris (Parsons et al. 2010): not a clock at 120 s"),
    ("DE CVn", "WD + dM; 2026-03/04 (clock era F)",
     "only one long-baseline ephemeris found (Yates et al. 2026) and that "
     "paper states the systems deviate significantly from linear; no "
     "second ephemeris to measure the disagreement"),
    ("HAT-P-18 b (2024-05-28), WASP-80 b (2025-06-22)", "transits",
     "no predicted transit falls inside the observed window (the census "
     "table records the nearest event and its margins)"),
)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[S3b] {msg}", flush=True)


def git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        return ""


def open_ro(path: Path) -> sqlite3.Connection:
    """Open a SQLite database strictly read-only (ground rule)."""
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


def open_db(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(path, timeout=300)
    con.execute("PRAGMA journal_mode = WAL")
    return con


def replace_table(con, name: str, create_sql: str, rows, n_cols: int) -> None:
    """Rebuild one derived table atomically (temp name, then swap)."""
    tmp = f"{name}__new"
    con.execute(f"DROP TABLE IF EXISTS {tmp}")
    con.execute(create_sql.format(table=tmp))
    marks = ",".join("?" * n_cols)
    con.executemany(f"INSERT INTO {tmp} VALUES ({marks})", rows)
    con.commit()
    con.execute("BEGIN IMMEDIATE")
    con.execute(f"DROP TABLE IF EXISTS {name}")
    con.execute(f"ALTER TABLE {tmp} RENAME TO {name}")
    con.commit()


def write_meta(con, extra: dict) -> None:
    con.execute("""CREATE TABLE IF NOT EXISTS s3b_build_meta
                   (key TEXT PRIMARY KEY, value TEXT)""")
    base = {"built_utc": datetime.now(timezone.utc).isoformat(),
            "code_version": ct.S3B_CODE_VERSION, "git_commit": git_commit()}
    for k, v in {**base, **extra}.items():
        con.execute("INSERT OR REPLACE INTO s3b_build_meta VALUES (?, ?)",
                    (k, str(v)))
    con.commit()


def curl(url: str, params: dict | None = None, timeout: int = 120) -> bytes:
    """HTTP GET through curl (the archive services reject Python's default
    client but answer curl); raises on a non-zero exit."""
    cmd = ["curl", "-s", "-f", "-m", str(timeout), "-G", url]
    for k, v in (params or {}).items():
        cmd += ["--data-urlencode", f"{k}={v}"]
    out = subprocess.run(cmd, capture_output=True, timeout=timeout + 30)
    if out.returncode != 0:
        raise RuntimeError(f"curl failed ({out.returncode}) for {url}")
    return out.stdout


def sexagesimal_to_deg(ra_hms: str, dec_dms: str) -> tuple[float, float]:
    h, m, s = (float(x) for x in ra_hms.split(":"))
    sign = -1.0 if dec_dms.strip().startswith("-") else 1.0
    d, dm, ds = (abs(float(x)) for x in dec_dms.split(":"))
    return 15.0 * (h + m / 60 + s / 3600), sign * (d + dm / 60 + ds / 3600)


def ang_sep_arcmin(ra1, dec1, ra2, dec2):
    """Great-circle separation (arcmin), vectorised (haversine)."""
    r1, d1, r2, d2 = (np.radians(np.asarray(v, dtype=float))
                      for v in (ra1, dec1, ra2, dec2))
    h = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) \
        * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))) * 60.0


# ---------------------------------------------------------------------------
# Stage: ephemerides
# ---------------------------------------------------------------------------
EXOCLOCK_FIELDS = (
    "name", "star", "ra_j2000", "dec_j2000", "ephem_mid_time",
    "ephem_mid_time_e1", "ephem_mid_time_format", "ephem_period",
    "ephem_period_e1", "ephem_parameters_ref", "inclination",
    "rp_over_rs", "sma_over_rs", "transit_parameters_ref", "depth_r_mmag",
    "duration_hours", "r_mag", "eccentricity")


def refresh_ephemerides(manifest: Path) -> None:
    """Fetch the literature inputs and rewrite the cache file.

    1. The ExoClock catalogue (Kokori et al. 2023, ApJS 265, 4, and its
       living database): one ephemeris per planet in BJD_TDB with errors,
       plus the transit geometry.  Only planets within
       :data:`CONE_RADIUS_ARCMIN` of an archived series are kept.
    2. For each kept planet, the most precise NON-ExoClock ephemeris in
       the NASA Exoplanet Archive ``ps`` table whose time system is
       BJD-TDB: an independent prediction to cross-check against.
    """
    raw = curl(EXOCLOCK_URL)
    cat = json.loads(raw)
    names = list(cat)
    coords = np.array([sexagesimal_to_deg(cat[n]["ra_j2000"],
                                          cat[n]["dec_j2000"])
                       for n in names])
    with open_ro(manifest) as con:
        series = con.execute(f"""
            SELECT avg(ra_deg), avg(dec_deg) FROM frames
            WHERE {SCIENCE_WHERE} AND jd IS NOT NULL AND ra_deg IS NOT NULL
            GROUP BY canonical_target, night, filter, readoutm, era_id, tree
            HAVING count(*) >= ? AND (max(jd) - min(jd)) * 24 >= ?""",
            (ct.MIN_SERIES_FRAMES, ct.MIN_SERIES_SPAN_H)).fetchall()
    sra = np.array([s[0] for s in series])
    sdec = np.array([s[1] for s in series])
    keep = {}
    for n, (ra, dec) in zip(names, coords):
        if (ang_sep_arcmin(sra, sdec, ra, dec) <= CONE_RADIUS_ARCMIN).any():
            p = {k: cat[n].get(k) for k in EXOCLOCK_FIELDS}
            p["ra_deg"], p["dec_deg"] = float(ra), float(dec)
            keep[n] = p
    log(f"ephemerides: {len(keep)} ExoClock planets lie in archived fields")
    alt = {}
    for n in keep:
        nasa = n[:-1] + " " + n[-1]                 # 'WASP-43b' -> 'WASP-43 b'
        q = ("select pl_name,pl_tranmid,pl_tranmiderr1,pl_orbper,"
             "pl_orbpererr1,pl_refname,pl_tsystemref from ps "
             f"where pl_name='{nasa}'")
        try:
            txt = curl(NASA_TAP_URL, {"query": q, "format": "json"}).decode()
            rows = json.loads(txt)
        except Exception as e:                       # noqa: BLE001
            log(f"  NASA archive query failed for {nasa}: {e}")
            continue
        best = None
        for r in rows:
            if None in (r.get("pl_tranmid"), r.get("pl_tranmiderr1"),
                        r.get("pl_orbper"), r.get("pl_orbpererr1")):
                continue
            if (r.get("pl_tsystemref") or "").upper() != "BJD-TDB":
                continue
            ref = _strip_tags(r.get("pl_refname") or "")
            if "Kokori" in ref or "ExoClock" in ref:
                continue                             # must be independent
            if best is None or r["pl_orbpererr1"] < best["pl_orbpererr1"]:
                best = dict(r, ref=ref)
        if best:
            alt[n] = {"t0": best["pl_tranmid"],
                      "sig_t0": best["pl_tranmiderr1"],
                      "period": best["pl_orbper"],
                      "sig_p": best["pl_orbpererr1"], "ref": best["ref"]}
    PRODUCT_DIR.mkdir(parents=True, exist_ok=True)
    payload = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "exoclock_url": EXOCLOCK_URL,
        "exoclock_sha256": hashlib.sha256(raw).hexdigest(),
        "exoclock_n_planets": len(cat),
        "nasa_tap_url": NASA_TAP_URL,
        "planets": keep, "nasa_alt": alt,
    }
    EPHEMERIS_CACHE.write_text(json.dumps(payload, indent=1, sort_keys=True))
    log(f"ephemerides: cache written ({len(alt)} with an independent "
        f"NASA-archive cross-check)")


def _strip_tags(s: str) -> str:
    import re
    return re.sub(r"<[^>]+>", "", s).strip()


def load_targets() -> dict:
    """All clock targets as one dict: name -> description.

    Planets come from the cache file, eclipsing binaries from
    :data:`EB_TARGETS`.  Every target carries ``primary`` (and maybe
    ``alt``) ephemerides in BJD_TDB, coordinates, and an ``EventShape``
    recipe.
    """
    if not EPHEMERIS_CACHE.exists():
        raise SystemExit(f"{EPHEMERIS_CACHE} is missing; run with "
                         "--refresh-ephemerides once (needs network).")
    cache = json.loads(EPHEMERIS_CACHE.read_text())
    targets = {}
    for n, p in cache["planets"].items():
        if (p.get("ephem_mid_time_format") or "") != "BJD_TDB":
            continue
        if (p.get("eccentricity") or 0.0) > 0.1:
            continue        # the circular-orbit template would not be exact
        a_r = float(p["sma_over_rs"])
        b = a_r * math.cos(math.radians(float(p["inclination"])))
        targets[n] = {
            "aliases": p.get("star") or "",
            "kind": "transit", "class": "transiting planet",
            "ra_deg": p["ra_deg"], "dec_deg": p["dec_deg"],
            "primary": {
                "t0": float(p["ephem_mid_time"]),
                "sig_t0": float(p["ephem_mid_time_e1"]),
                "period": float(p["ephem_period"]),
                "sig_p": float(p["ephem_period_e1"]),
                "quad": 0.0, "sig_quad": 0.0,
                "ref": ("ExoClock database (Kokori et al. 2023, ApJS 265, "
                        f"4), ephemeris ref {p.get('ephem_parameters_ref')}"
                        f", retrieved {cache['retrieved_utc'][:10]}"),
            },
            "alt": ({**cache["nasa_alt"][n], "quad": 0.0, "sig_quad": 0.0,
                     "ref": cache["nasa_alt"][n]["ref"]
                     + " (NASA Exoplanet Archive)"}
                    if n in cache["nasa_alt"] else None),
            "sys_floor_s": 0.0,
            "shape": {"a_over_r": a_r, "k": float(p["rp_over_rs"]),
                      "b": float(b), "u1": LD_PLANET[0], "u2": LD_PLANET[1],
                      "f1": 1.0, "free": ("k",), "baseline_order": 1},
            "min_depth": None,
        }
    for n, e in EB_TARGETS.items():
        targets[n] = dict(e)
    return targets


def shape_of(target: dict) -> ct.EventShape:
    s = target["shape"]
    return ct.EventShape(target["primary"]["period"], s["a_over_r"], s["k"],
                         s["b"], s["u1"], s["u2"], s["f1"], s["free"],
                         s["baseline_order"])


def stage_ephemerides(db: sqlite3.Connection) -> dict:
    targets = load_targets()
    cache = json.loads(EPHEMERIS_CACHE.read_text())
    rows = []
    for n, t in sorted(targets.items()):
        for role in ("primary", "alt"):
            e = t.get(role)
            if not e:
                continue
            rows.append((n, role, t["kind"], t["class"], t["ra_deg"],
                         t["dec_deg"], e["t0"], e["sig_t0"], e["period"],
                         e["sig_p"], e.get("quad", 0.0),
                         e.get("sig_quad", 0.0), "BJD_TDB", e["ref"]))
    replace_table(db, "s3b_ephemeris", """CREATE TABLE {table} (
        target TEXT, role TEXT, kind TEXT, class TEXT, ra_deg REAL,
        dec_deg REAL, t0_bjd REAL, sig_t0_d REAL, period_d REAL,
        sig_period_d REAL, quad_d REAL, sig_quad_d REAL, time_system TEXT,
        reference TEXT, PRIMARY KEY (target, role))""", rows, 14)
    replace_table(db, "s3b_refused", """CREATE TABLE {table} (
        target TEXT PRIMARY KEY, what TEXT, reason TEXT)""",
                  list(REFUSED_TARGETS), 3)
    write_meta(db, {"ephemeris_retrieved_utc": cache["retrieved_utc"],
                    "exoclock_sha256": cache["exoclock_sha256"],
                    "exoclock_n_planets": cache["exoclock_n_planets"],
                    "n_targets": len(targets)})
    log(f"ephemerides: {len(targets)} clock targets "
        f"({sum(t['kind'] == 'transit' for t in targets.values())} planets)")
    return targets


# ---------------------------------------------------------------------------
# Stage: census
# ---------------------------------------------------------------------------
def series_id_of(target, night, filt, readoutm, era_id, tree,
                 manifest_target) -> str:
    """Unique key of one photometric series.  The manifest's own target
    name is part of it because one star can be requested twice on one
    night under two names (``WASP 52`` / ``WASP-52``), and those are two
    separate runs."""
    mode = (readoutm or "blank").replace(" ", "")
    req = (manifest_target or "none").replace(" ", "_")
    return f"{target}|{night}|{filt}|{mode}|e{era_id}|{tree}|{req}"


def stage_census(db: sqlite3.Connection, manifest: Path, targets: dict,
                 ephemeris: str) -> None:
    """Find the series, predict the events, admit or refuse each."""
    names = sorted(targets)
    tra = np.array([targets[n]["ra_deg"] for n in names])
    tdec = np.array([targets[n]["dec_deg"] for n in names])
    with open_ro(manifest) as con:
        groups = con.execute(f"""
            SELECT canonical_target, night, filter, readoutm, era_id, tree,
                   count(*), min(jd), max(jd), avg(exptime), avg(ra_deg),
                   avg(dec_deg), max(xbinning)
            FROM frames
            WHERE {SCIENCE_WHERE} AND jd IS NOT NULL AND exptime > 0
              AND ra_deg IS NOT NULL
            GROUP BY canonical_target, night, filter, readoutm, era_id, tree
            HAVING count(*) >= ? AND (max(jd) - min(jd)) * 24 >= ?""",
            (ct.MIN_SERIES_FRAMES, ct.MIN_SERIES_SPAN_H)).fetchall()
    series_rows, event_rows = [], []
    for (ctarget, night, filt, readoutm, era_id, tree, n, jd0, jd1, exp_avg,
         ra, dec, xbin) in groups:
        if filt in GRISM_FILTERS:
            continue
        sep = ang_sep_arcmin(tra, tdec, ra, dec)
        j = int(np.argmin(sep))
        if sep[j] > CONE_RADIUS_ARCMIN:
            continue
        name = names[j]
        t = targets[name]
        eph = t["primary"]
        shape = shape_of(t)
        sid = series_id_of(name, night, filt, readoutm, era_id, tree,
                           ctarget)
        # Window in BJD_TDB at the target's own coordinates.
        (b0, b1), _, _ = tm.bjd_tdb_from_utc(
            np.array([jd0, jd1 + exp_avg / 86400.0]),
            t["ra_deg"], t["dec_deg"], ephemeris=ephemeris)
        t14 = shape.duration_d()
        epochs = ct.events_in_window(b0 - t14, b1 + t14, eph["t0"],
                                     eph["period"])
        if not epochs:                    # record the nearest, as evidence
            epochs = [ct.nearest_epoch(0.5 * (b0 + b1), eph["t0"],
                                       eph["period"])]
        n_adm = 0
        for e in epochs:
            tp, sig = ct.predict_event(eph["t0"], eph["sig_t0"],
                                       eph["period"], eph["sig_p"], e,
                                       eph.get("quad", 0.0),
                                       eph.get("sig_quad", 0.0))
            pre_h = (tp - t14 / 2 - b0) * 24.0
            post_h = (b1 - (tp + t14 / 2)) * 24.0
            sig_s = sig * 86400.0
            reason = ""
            need_h = (ct.MIN_MARGIN_H if t["kind"] == "transit"
                      else ct.MIN_MARGIN_ECLIPSE_H)
            if min(pre_h, post_h) < need_h - ct.CLOCK_SEARCH_PAD_H:
                reason = (f"event not covered even allowing a "
                          f"{ct.CLOCK_SEARCH_PAD_H * 3600:.0f} s clock "
                          f"error (margins {pre_h:+.2f} h / {post_h:+.2f} "
                          f"h; need {need_h} h each side)")
            elif sig_s > ct.EPH_SIGMA_MAX_S:
                reason = (f"ephemeris prediction too loose "
                          f"({sig_s:.0f} s > {ct.EPH_SIGMA_MAX_S:.0f} s)")
            admitted = int(reason == "")
            n_adm += admitted
            event_rows.append((f"{sid}#{e}", sid, name, int(e), tp, sig_s,
                               pre_h, post_h, t14 * 1440.0, admitted,
                               reason))
        series_rows.append((sid, name, t["kind"], ctarget, night, filt,
                            readoutm, int(era_id), tree,
                            ct.clock_era(night, era_id), int(n), jd0, jd1,
                            exp_avg, float(sep[j]), int(xbin or 1),
                            int(n_adm > 0)))
    replace_table(db, "s3b_series", """CREATE TABLE {table} (
        series_id TEXT PRIMARY KEY, target TEXT, kind TEXT,
        manifest_target TEXT, night TEXT, filter TEXT, readoutm TEXT,
        era_id INTEGER, tree TEXT, clock_era TEXT, n_frames INTEGER,
        jd_start REAL, jd_end REAL, exptime_s REAL, offset_arcmin REAL,
        binning INTEGER, admitted INTEGER)""", series_rows, 17)
    replace_table(db, "s3b_census_events", """CREATE TABLE {table} (
        event_id TEXT PRIMARY KEY, series_id TEXT, target TEXT,
        epoch INTEGER, t_pred_bjd REAL, sig_pred_s REAL, pre_h REAL,
        post_h REAL, t14_min REAL, admitted INTEGER, reason TEXT)""",
                  event_rows, 11)
    n_adm_s = sum(r[-1] for r in series_rows)
    n_adm_e = sum(r[9] for r in event_rows)
    write_meta(db, {"census_series": len(series_rows),
                    "census_series_admitted": n_adm_s,
                    "census_events": len(event_rows),
                    "census_events_admitted": n_adm_e,
                    "bjd_ephemeris": ephemeris})
    log(f"census: {len(series_rows)} series on clock targets, "
        f"{n_adm_s} admitted with {n_adm_e} events")


# ---------------------------------------------------------------------------
# Stage: photometry
# ---------------------------------------------------------------------------
#: Header cards that define a frame's sky projection (prefix match).
_WCS_PREFIXES = ("CTYPE", "CRVAL", "CRPIX", "CD1_", "CD2_", "CDELT",
                 "CROTA", "PC1_", "PC2_", "CUNIT", "A_", "B_", "AP_",
                 "BP_", "RADESYS", "EQUINOX", "LONPOLE", "LATPOLE")


def read_frame(path: Path):
    """Pixels (float32) and a PLAIN-DICT header of one archive frame.

    The header is copied card by card into a dict, each read guarded:
    some archive headers carry malformed CONTINUE cards that make astropy
    raise on a whole-header copy, and one unreadable comment card must not
    cost a frame its pixels.  Only cards that can be read are kept.
    """
    from astropy.io import fits
    with fits.open(path) as hdul:
        hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None \
            else hdul[0]
        data = np.ascontiguousarray(hdu.data, dtype=np.float32)
        hdr = {}
        for key in list(hdu.header.keys()):
            if not key or key in ("COMMENT", "HISTORY", "CONTINUE"):
                continue
            try:
                hdr[key] = hdu.header[key]
            except Exception:                        # noqa: BLE001
                continue
    return data, hdr


def wcs_of(hdr: dict, shape: tuple):
    """An astropy WCS built from the projection cards of a header dict
    (None when the frame carries no celestial WCS)."""
    from astropy.io import fits
    from astropy.wcs import WCS
    mini = fits.Header()
    mini["NAXIS"] = 2
    mini["NAXIS1"], mini["NAXIS2"] = int(shape[1]), int(shape[0])
    for key, val in hdr.items():
        if key.startswith(_WCS_PREFIXES):
            try:
                mini[key] = val
            except Exception:                        # noqa: BLE001
                continue
    if "CTYPE1" not in mini:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        wcs = WCS(mini)
    return wcs if wcs.has_celestial else None


def pick_masters(con_manifest, archive: Path, readoutm: str, filt: str,
                 exptime_s: float, jd: float) -> dict:
    """The master dark and master flat that serve one RAW series.

    Same rules as the CV campaign (``run_cv_photometry._pick_frame_masters``),
    read from S0b's ``calib_frames`` inventory: the dark must come from the
    same readout mode and match the exposure time (no scaling: the archived
    masters are bias-inclusive); the flat must come from the same readout
    mode and the same filter, case-sensitively; among equals the one
    nearest in time wins (``macro_phot.series.pick_master``).  Either may be
    absent — the recipe actually applied is recorded per series.

    Why calibrate at all for a mid-time?  Because these telescopes were not
    guided: the stars walk tens of pixels in a night and are re-centred in
    jumps, and on unflattened frames every jump is a 1-2 % step in the
    differential light curve — the size of the transit itself.
    """
    from macro_phot import series as sr
    rows = con_manifest.execute(
        """SELECT kind, jd, path, filter, exptime FROM calib_frames
           WHERE is_master = 1 AND kind IN ('dark', 'flat')
             AND readoutm IS ?""", (readoutm,)).fetchall()
    darks = [(r[1], r[2]) for r in rows if r[0] == "dark"
             and sr.dark_exptime_matches(r[4], exptime_s)]
    flats = [(r[1], r[2]) for r in rows if r[0] == "flat" and r[3] == filt]
    d = sr.pick_master(darks, jd)
    f = sr.pick_master(flats, jd)

    def age(chosen):
        return (abs(jd - chosen[0]) if chosen and chosen[0] is not None
                else None)
    return {"dark": str(archive / d[1]) if d else None,
            "flat": str(archive / f[1]) if f else None,
            "dark_age_d": age(d), "flat_age_d": age(f)}


def calibrate(raw: np.ndarray, dark_path, flat_path):
    """Apply the series' masters (``macro_phot.calib``); a master whose
    shape disagrees with the frame is skipped, and the recipe string says
    exactly what was applied."""
    from macro_phot import calib as cal
    dark = cal.read_master(Path(dark_path)) if dark_path else None
    flat = cal.read_master(Path(flat_path)) if flat_path else None
    if dark is not None and dark.shape != raw.shape:
        dark = None
    if flat is not None and flat.shape != raw.shape:
        flat = None
    return cal.apply_masters(raw, dark, flat)


def detect(data: np.ndarray):
    """Background-subtract and detect bright sources (sep).

    NaN pixels (flat-divided dead pixels in the reduced tree) are replaced
    by the frame median first — sep cannot digest them — which removes
    them from every aperture sum they touch only in the sense that they
    contribute sky; a comparison star sitting on one is caught by the
    stability pruning.
    """
    import sep
    if not np.isfinite(data).all():
        data = np.where(np.isfinite(data), data,
                        np.float32(np.nanmedian(data)))
        data = np.ascontiguousarray(data, dtype=np.float32)
    sep.set_extract_pixstack(1_000_000)
    bkg = sep.Background(data)
    sub = data - bkg.back()
    objs = sep.extract(sub, DETECT_SIGMA, err=bkg.globalrms,
                       minarea=DETECT_MINAREA)
    order = np.argsort(-objs["flux"])
    return data, sub, bkg, objs[order]


def fwhm_of(objs) -> np.ndarray:
    return 2.3548 * np.sqrt((objs["a"] ** 2 + objs["b"] ** 2) / 2.0)


def prepare_series(archive: Path, frames: list, target: dict,
                   veto_adu: float, binning: int, masters: dict) -> dict:
    """Choose the reference frame, find the target on it, pick comparisons.

    Candidates are tried from the middle of the run outwards, plate-solved
    frames first.  On each: detect sources, read the frame's own WCS,
    project the target's catalogue position, and take the nearest
    detection within :data:`TARGET_FIND_ARCSEC`.  Comparison candidates
    are the brightest detections that are clear of the frame edge, clear
    of the saturation veto with head-room, and without a neighbour inside
    their sky annulus brighter than a tenth of themselves.
    """
    from astropy.wcs.utils import proj_plane_pixel_scales
    mid = len(frames) // 2
    order = sorted(range(len(frames)), key=lambda i: abs(i - mid))
    last_reason = "no usable reference frame"
    for i in order[:MAX_REF_TRIES]:
        path = frames[i][0]
        try:
            raw, hdr = read_frame(archive / path)
            data, recipe = calibrate(raw, masters["dark"], masters["flat"])
            data, sub, bkg, objs = detect(data)
        except Exception as e:                       # noqa: BLE001
            last_reason = f"reference unreadable: {type(e).__name__}"
            continue
        if len(objs) < 8:
            last_reason = "fewer than 8 sources on the reference"
            continue
        wcs = wcs_of(hdr, data.shape)
        if wcs is None:
            last_reason = "reference frame has no WCS"
            continue
        scale = float(np.mean(proj_plane_pixel_scales(wcs.celestial))
                      * 3600.0)
        tx, ty = wcs.celestial.all_world2pix(target["ra_deg"],
                                             target["dec_deg"], 0)
        d = np.hypot(objs["x"] - tx, objs["y"] - ty)
        it = int(np.argmin(d))
        if d[it] * scale > TARGET_FIND_ARCSEC:
            last_reason = (f"no detection within {TARGET_FIND_ARCSEC:.0f}"
                           f"\" of the target position")
            continue
        fwhm = float(np.median(fwhm_of(objs[:50])))
        radii = [max(ct.APERTURE_MIN_PX, f * fwhm)
                 for f in ct.APERTURE_FWHM_FACTORS]
        ann = (max(ct.ANNULUS_FWHM_FACTORS[0] * fwhm, radii[-1] + 2.0),
               max(ct.ANNULUS_FWHM_FACTORS[1] * fwhm, radii[-1] + 6.0))
        h, w = data.shape
        margin = EDGE_MARGIN_PX + ann[1]
        ratio = ct.native_peak_ratio(fwhm, binning)
        comps = []
        for j in range(len(objs)):
            if j == it or len(comps) >= ct.N_COMP_CANDIDATES:
                continue
            x, y = float(objs["x"][j]), float(objs["y"][j])
            if min(x, y, w - x, h - y) < margin:
                continue
            # Saturation is judged on the RAW frame (native ADU), never
            # on dark-subtracted, flat-divided values.
            xi, yi = int(round(x)), int(round(y))
            peak = float(raw[max(yi - 2, 0):yi + 3,
                             max(xi - 2, 0):xi + 3].max()) * ratio
            if peak >= COMP_PEAK_HEADROOM * veto_adu:
                continue
            dd = np.hypot(objs["x"] - x, objs["y"] - y)
            near = (dd < ann[1]) & (dd > 0)
            if (objs["flux"][near] > 0.1 * objs["flux"][j]).any():
                continue
            if np.hypot(x - objs["x"][it], y - objs["y"][it]) < ann[1]:
                continue
            comps.append(j)
        if len(comps) < ct.MIN_COMPS:
            last_reason = "fewer than 2 clean comparison stars"
            continue
        idx = [it] + comps
        xy = np.column_stack([objs["x"][idx], objs["y"][idx]]).astype(float)
        sky = wcs.celestial.all_pix2world(xy[:, 0], xy[:, 1], 0)
        return {"ok": True, "ref_path": path, "fwhm_px": fwhm,
                "recipe": recipe,
                "scale_arcsec": scale, "radii": radii, "annulus": ann,
                "star_xy": xy, "star_ra": np.asarray(sky[0]),
                "star_dec": np.asarray(sky[1]),
                "star_flux": objs["flux"][idx].astype(float),
                "ref_bright": np.column_stack(
                    [objs["x"][:800], objs["y"][:800]]).astype(float),
                "target_offset_arcsec": float(d[it] * scale)}
    return {"ok": False, "reason": last_reason}


def measure_one(task: dict) -> dict:
    """Worker: forced multi-aperture photometry of one frame.

    Registration is by astroalign triangle matching of the frame's bright
    detections against the reference frame's (``macro_phot.extract.
    find_series_transform`` — seeded, so re-runs are bit-identical); no
    per-frame WCS is trusted.  Each star's aperture is centred on the
    detection nearest its transformed reference position when one lies
    within the match tolerance, else on the transformed position itself
    (flagged ``matched = 0``).
    """
    import sep
    from macro_phot import extract as ex
    out = {"path": task["path"], "status": "ok", "stars": []}
    try:
        raw, hdr = read_frame(Path(task["archive"]) / task["path"])
        data, _recipe = calibrate(raw, task["dark"], task["flat"])
        data, sub, bkg, objs = detect(data)
    except Exception as e:                           # noqa: BLE001
        out["status"] = f"unreadable: {type(e).__name__}"
        return out
    out.update(bkg_adu=float(bkg.globalback), n_det=int(len(objs)),
               airmass=_float_or_none(hdr.get("AIRMASS")))
    if len(objs) < 5:
        out["status"] = "few_sources"
        return out
    bright = np.column_stack([objs["x"], objs["y"]]).astype(float)
    ref_bright = np.asarray(task["ref_bright"])
    try:
        tf = ex.find_series_transform(bright, ref_bright, seed=task["seed"],
                                      attempts=ex.PRODUCTION_ALIGN_ATTEMPTS)
    except Exception as e:                           # noqa: BLE001
        out["status"] = f"align_failed: {type(e).__name__}"
        return out
    pred = np.asarray(tf.inverse(np.asarray(task["star_xy"])))
    tol = task["tol_px"]
    xs, ys, matched = [], [], []
    for px, py in pred:
        d = np.hypot(objs["x"] - px, objs["y"] - py)
        j = int(np.argmin(d))
        if d[j] <= tol:
            xs.append(float(objs["x"][j]))
            ys.append(float(objs["y"][j]))
            matched.append(1)
        else:
            xs.append(float(px))
            ys.append(float(py))
            matched.append(0)
    xs, ys = np.array(xs), np.array(ys)
    h, w = data.shape
    egain = _float_or_none(hdr.get("EGAIN"))
    gain = egain if egain and egain > 0 else None
    # A star whose sky annulus leaves the frame (drift, a re-pointing) is
    # not measured on this frame: sep refuses an annulus with no pixels,
    # and a truncated annulus is a different sky estimate anyway.
    a_out = float(task["annulus"][1])
    inside = ((xs >= a_out) & (xs < w - a_out)
              & (ys >= a_out) & (ys < h - a_out))
    fl = np.full((len(task["radii"]), len(xs)), np.nan)
    fe = np.full((len(task["radii"]), len(xs)), np.nan)
    if inside.any():
        for ir, r in enumerate(task["radii"]):
            f, e, _flag = sep.sum_circle(sub, xs[inside], ys[inside], r,
                                         err=bkg.globalrms, gain=gain,
                                         bkgann=task["annulus"])
            fl[ir, inside] = f
            fe[ir, inside] = e
    m = fwhm_of(objs)
    frame_fwhm = float(np.median(m[:50]))
    out["fwhm_px"] = frame_fwhm
    rpk = int(math.ceil(task["radii"][0]))
    for s in range(len(xs)):
        xi, yi = int(round(xs[s])), int(round(ys[s]))
        on = bool(inside[s])
        # Peak from the RAW pixels: the saturation veto is a statement
        # about the ADC, in native ADU (standing rule 4).
        peak = float(raw[yi - rpk:yi + rpk + 1,
                         xi - rpk:xi + rpk + 1].max()) if on else None
        out["stars"].append({
            "star": s, "x": float(xs[s]), "y": float(ys[s]),
            "matched": matched[s], "peak": peak, "on_frame": int(on),
            "flux": [float(v) if np.isfinite(v) else None
                     for v in fl[:, s]],
            "err": [float(v) if np.isfinite(v) else None
                    for v in fe[:, s]]})
    return out


def _float_or_none(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None


def veto_for(con_manifest, readoutm: str) -> float:
    """S2 saturation veto (raw ADU) for one readout mode, from the
    manifest's ``detector_params`` (never typed here).  The 5 MHz iKon
    mode has no S2 row and inherits the 1 MHz iKon value; anything else
    unknown falls back to the smallest veto on file (most conservative)."""
    rows = dict(con_manifest.execute(
        "SELECT era_group, value FROM detector_params "
        "WHERE quantity = 'saturation_veto_adu'").fetchall())
    key = (readoutm or "").strip() or "(blank 2026)"
    if key in rows:
        return float(rows[key])
    if "High Sensitivity" in key:
        return float(rows["1MHz High Sensitivity 16-bit"])
    return float(min(rows.values()))


def stage_photometry(db: sqlite3.Connection, manifest: Path, archive: Path,
                     targets: dict, ephemeris: str, workers: int,
                     only: str | None) -> None:
    db.execute("""CREATE TABLE IF NOT EXISTS s3b_phot_series (
        series_id TEXT PRIMARY KEY, status TEXT, ref_path TEXT,
        ref_fwhm_px REAL, scale_arcsec REAL, r0_px REAL, r1_px REAL,
        r2_px REAL, ann_in_px REAL, ann_out_px REAL, n_stars INTEGER,
        veto_adu REAL, target_offset_arcsec REAL, n_frames_listed INTEGER,
        n_withdrawn INTEGER, recipe TEXT, dark_path TEXT, flat_path TEXT,
        dark_age_d REAL, flat_age_d REAL)""")
    db.execute("""CREATE TABLE IF NOT EXISTS s3b_stars (
        series_id TEXT, star INTEGER, ra_deg REAL, dec_deg REAL,
        x_ref REAL, y_ref REAL, ref_flux REAL,
        PRIMARY KEY (series_id, star))""")
    db.execute("""CREATE TABLE IF NOT EXISTS s3b_frames (
        series_id TEXT, path TEXT, jd_utc_start REAL, exptime_s REAL,
        jd_utc_mid REAL, mid_method TEXT, bjd_tdb REAL, bary_ltt_s REAL,
        airmass REAL, fwhm_px REAL, bkg_adu REAL, n_det INTEGER,
        status TEXT, PRIMARY KEY (series_id, path))""")
    db.execute("""CREATE TABLE IF NOT EXISTS s3b_flux (
        series_id TEXT, path TEXT, star INTEGER, x REAL, y REAL,
        matched INTEGER, on_frame INTEGER, peak_adu REAL,
        flux0 REAL, flux1 REAL, flux2 REAL, err0 REAL, err1 REAL,
        err2 REAL, PRIMARY KEY (series_id, path, star))""")
    db.commit()
    series = db.execute("""SELECT series_id, target, night, filter, readoutm,
                                  era_id, tree, manifest_target, binning
                           FROM s3b_series WHERE admitted = 1
                           ORDER BY night, series_id""").fetchall()
    if only:
        series = [s for s in series if only in s[0]]
    con = open_ro(manifest)
    pool = mp.get_context("spawn").Pool(workers) if workers > 1 else None
    try:
        for k, (sid, name, night, filt, readoutm, era_id, tree, mtarget,
                binning) in enumerate(series, 1):
            target = targets[name]
            # Frames of the series; any frame whose time S3 withdrew
            # (raw-vs-reduced stamp disagreement) is excluded up front.
            frames = con.execute(f"""
                SELECT f.path, f.jd, f.exptime FROM frames f
                WHERE {SCIENCE_WHERE} AND f.jd IS NOT NULL AND f.exptime > 0
                  AND f.night = ? AND f.filter = ? AND f.era_id = ?
                  AND f.tree = ? AND f.readoutm IS ?
                  AND f.canonical_target IS ?
                ORDER BY f.jd""",
                (night, filt, era_id, tree, readoutm, mtarget)).fetchall()
            withdrawn = {r[0] for r in con.execute(
                "SELECT path FROM s3_time_outliers").fetchall()}
            n_with = sum(f[0] in withdrawn for f in frames)
            frames = [f for f in frames if f[0] not in withdrawn]
            st = db.execute("SELECT status FROM s3b_phot_series "
                            "WHERE series_id = ?", (sid,)).fetchone()
            if st and st[0] != "ok":
                continue                     # already refused; keep verdict
            veto = veto_for(con, readoutm)
            if not st:
                # Raw series get the era's masters; the server-reduced
                # tree is already dark-subtracted and flat-fielded.
                if tree == "reduced" or not frames:
                    masters = {"dark": None, "flat": None,
                               "dark_age_d": None, "flat_age_d": None}
                else:
                    mid = frames[len(frames) // 2]
                    masters = pick_masters(con, archive, readoutm, filt,
                                           mid[2], mid[1])
                prep = prepare_series(archive, frames, target, veto,
                                      binning, masters)
                if not prep["ok"]:
                    db.execute("INSERT OR REPLACE INTO s3b_phot_series "
                               "(series_id, status, n_frames_listed, "
                               " n_withdrawn) VALUES (?,?,?,?)",
                               (sid, prep["reason"], len(frames), n_with))
                    db.commit()
                    log(f"[{k}/{len(series)}] {sid}: REFUSED — "
                        f"{prep['reason']}")
                    continue
                rel = lambda p: (str(Path(p).relative_to(archive))  # noqa: E731
                                 if p else None)
                db.execute("INSERT OR REPLACE INTO s3b_phot_series VALUES "
                           "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                           (sid, "ok", prep["ref_path"], prep["fwhm_px"],
                            prep["scale_arcsec"], *prep["radii"],
                            *prep["annulus"], len(prep["star_xy"]), veto,
                            prep["target_offset_arcsec"], len(frames),
                            n_with,
                            ("server_reduced" if tree == "reduced"
                             else prep["recipe"]),
                            rel(masters["dark"]) if "dark" in prep["recipe"]
                            else None,
                            rel(masters["flat"]) if "flat" in prep["recipe"]
                            else None,
                            masters["dark_age_d"], masters["flat_age_d"]))
                db.executemany(
                    "INSERT OR REPLACE INTO s3b_stars VALUES (?,?,?,?,?,?,?)",
                    [(sid, s, float(prep["star_ra"][s]),
                      float(prep["star_dec"][s]), float(prep["star_xy"][s, 0]),
                      float(prep["star_xy"][s, 1]),
                      float(prep["star_flux"][s]))
                     for s in range(len(prep["star_xy"]))])
                db.execute("INSERT OR REPLACE INTO s3b_ref_bright VALUES "
                           "(?, ?)", (sid, json.dumps(
                               prep["ref_bright"].round(3).tolist())))
                db.commit()
            ps = db.execute("""SELECT ref_fwhm_px, r0_px, r1_px, r2_px,
                                      ann_in_px, ann_out_px, dark_path,
                                      flat_path
                               FROM s3b_phot_series WHERE series_id = ?""",
                            (sid,)).fetchone()
            star_xy = db.execute("""SELECT x_ref, y_ref FROM s3b_stars
                                    WHERE series_id = ? ORDER BY star""",
                                 (sid,)).fetchall()
            ref_bright = json.loads(db.execute(
                "SELECT xy_json FROM s3b_ref_bright WHERE series_id = ?",
                (sid,)).fetchone()[0])
            done = {r[0] for r in db.execute(
                "SELECT path FROM s3b_frames WHERE series_id = ?", (sid,))}
            todo = [f for f in frames if f[0] not in done]
            if not todo:
                continue
            log(f"[{k}/{len(series)}] {sid}: {len(todo)} frames to measure")
            # Mid-exposure BJD_TDB at the TARGET's coordinates (S3
            # convention 4: frame_times points at the frame centre).
            jd = np.array([f[1] for f in todo])
            exps = np.array([f[2] for f in todo])
            mids = np.array([tm.jd_utc_mid(j, e, readoutm)[0]
                             for j, e in zip(jd, exps)])
            method = tm.mid_method_for(readoutm, float(jd[0]),
                                       float(exps[0]))
            bjd, ltt, _ = tm.bjd_tdb_from_utc(mids, target["ra_deg"],
                                              target["dec_deg"],
                                              ephemeris=ephemeris)
            meta = {f[0]: (f[1], f[2], float(m), float(b), float(lt))
                    for f, m, b, lt in zip(todo, mids, bjd, ltt)}
            tasks = [{"archive": str(archive), "path": f[0],
                      "ref_bright": ref_bright, "star_xy": star_xy,
                      "radii": list(ps[1:4]), "annulus": (ps[4], ps[5]),
                      "dark": str(archive / ps[6]) if ps[6] else None,
                      "flat": str(archive / ps[7]) if ps[7] else None,
                      "tol_px": max(4.0, ps[0]),
                      "seed": int(hashlib.sha1(f[0].encode()).hexdigest()[:8],
                                  16)}
                     for f in todo]
            results = (pool.imap_unordered(measure_one, tasks, chunksize=2)
                       if pool else map(measure_one, tasks))
            for n_done, res in enumerate(results, 1):
                jd0, ex_s, mid, b, lt = meta[res["path"]]
                db.execute("INSERT OR REPLACE INTO s3b_frames VALUES "
                           "(?,?,?,?,?,?,?,?,?,?,?,?,?)",
                           (sid, res["path"], jd0, ex_s, mid, method, b, lt,
                            res.get("airmass"), res.get("fwhm_px"),
                            res.get("bkg_adu"), res.get("n_det"),
                            res["status"]))
                db.executemany(
                    "INSERT OR REPLACE INTO s3b_flux VALUES "
                    "(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    [(sid, res["path"], s["star"], s["x"], s["y"],
                      s["matched"], s["on_frame"], s["peak"], *s["flux"],
                      *s["err"]) for s in res["stars"]])
                if n_done % 50 == 0:
                    db.commit()
            db.commit()
    finally:
        if pool:
            pool.close()
            pool.join()
        con.close()
    write_meta(db, {"photometry_series_ok": db.execute(
        "SELECT count(*) FROM s3b_phot_series WHERE status = 'ok'"
    ).fetchone()[0], "photometry_frames": db.execute(
        "SELECT count(*) FROM s3b_frames").fetchone()[0]})


# ---------------------------------------------------------------------------
# Stage: fit
# ---------------------------------------------------------------------------
def load_series_arrays(db, sid: str) -> dict | None:
    """Photometry of one series as arrays: frames x stars x apertures."""
    frames = db.execute("""SELECT path, bjd_tdb, exptime_s, fwhm_px, airmass
                           FROM s3b_frames
                           WHERE series_id = ? AND status = 'ok'
                           ORDER BY bjd_tdb""", (sid,)).fetchall()
    n_star = db.execute("SELECT count(*) FROM s3b_stars WHERE series_id = ?",
                        (sid,)).fetchone()[0]
    if len(frames) < ct.MIN_FIT_POINTS or n_star < 1 + ct.MIN_COMPS:
        return None
    index = {f[0]: i for i, f in enumerate(frames)}
    flux = np.full((len(frames), n_star, 3), np.nan)
    err = np.full((len(frames), n_star, 3), np.nan)
    peak = np.full((len(frames), n_star), np.nan)
    good = np.zeros((len(frames), n_star), dtype=bool)
    txy = np.full((len(frames), 2), np.nan)
    for (path, star, matched, on_frame, pk, f0, f1, f2, e0, e1, e2, sx,
         sy) in db.execute("""SELECT path, star, matched, on_frame, peak_adu,
                                 flux0, flux1, flux2, err0, err1, err2, x, y
                          FROM s3b_flux WHERE series_id = ?""", (sid,)):
        i = index.get(path)
        if i is None:
            continue
        if star == 0:
            txy[i] = (sx, sy)
        flux[i, star] = [np.nan if v is None else v for v in (f0, f1, f2)]
        err[i, star] = [np.nan if v is None else v for v in (e0, e1, e2)]
        peak[i, star] = pk if pk is not None else np.nan
        # A comparison star must be re-detected on the frame.  The TARGET
        # must not: in a total eclipse it drops below the detection
        # threshold, and discarding exactly those frames would cut the
        # bottom out of the eclipse being timed.  Its aperture then sits
        # at the position the frame-to-reference transform predicts.
        good[i, star] = bool(on_frame) and (bool(matched) or star == 0)
    return {"paths": [f[0] for f in frames],
            "bjd": np.array([f[1] for f in frames]),
            "exptime_s": np.array([f[2] for f in frames]),
            "fwhm": np.array([f[3] if f[3] else np.nan for f in frames]),
            "airmass": np.array([f[4] if f[4] else np.nan for f in frames]),
            "flux": flux, "err": err, "peak": peak, "good": good,
            "txy": txy}


def build_lightcurves(arr: dict, veto_adu: float, binning: int) -> dict:
    """Differential light curve per aperture + which aperture to adopt.

    Saturation (standing rule 4) is judged per frame and per star from the
    raw peak pixel, scaled to NATIVE pixels for the average-binned eras by
    :func:`clock_transits.native_peak_ratio` at that frame's own FWHM.
    """
    ratio = np.array([ct.native_peak_ratio(f if np.isfinite(f) else 3.0,
                                           binning) for f in arr["fwhm"]])
    native_peak = arr["peak"] * ratio[:, None]
    unsat = ~(native_peak >= veto_adu)
    usable = arr["good"] & unsat & np.isfinite(arr["flux"][:, :, 0])
    out = {"apertures": [], "n_target_saturated":
           int((~unsat[:, 0] & arr["good"][:, 0]).sum())}
    for a in range(3):
        fl = arr["flux"][:, :, a]
        ok = usable & (fl > 0)
        keep = ct.select_comparisons(fl[:, 1:], ok[:, 1:])
        if keep.sum() < ct.MIN_COMPS:
            out["apertures"].append(None)
            continue
        lc, use = ct.differential_lightcurve(fl[:, 0], ok[:, 0], fl[:, 1:],
                                             ok[:, 1:], keep)
        sd = ct.successive_difference_rms(lc[use])
        rel = arr["err"][:, 0, a] / np.where(fl[:, 0] > 0, fl[:, 0], np.nan)
        med = np.nanmedian(rel[use]) if use.any() else np.nan
        weight = np.clip(rel / med, 0.5, 3.0)
        out["apertures"].append({"flux": lc, "use": use, "sd_rms": sd,
                                 "err": sd * weight, "keep": keep})
    scores = [a["sd_rms"] if a and np.isfinite(a["sd_rms"]) and
              a["use"].sum() >= ct.MIN_FIT_POINTS else np.inf
              for a in out["apertures"]]
    out["adopted"] = int(np.argmin(scores)) if np.isfinite(min(scores)) \
        else None
    return out


def event_window(kind: str, t: np.ndarray, t_pred: float,
                 shape: ct.EventShape) -> np.ndarray:
    """Which points of a series belong to one event's fit.

    A transit uses the whole run (one event per run).  An eclipse of a
    short-period binary uses a window of half-width min(0.3 P, 3 T14)
    about the PREDICTED time: wide enough that a +-1,065 s displacement
    stays inside it, narrow enough to exclude the secondary eclipse.
    """
    if kind == "transit":
        return np.ones(len(t), dtype=bool)
    half = min(0.3 * shape.period, 3.0 * shape.duration_d())
    return np.abs(t - t_pred) <= half


def analyse_event(job: dict) -> dict:
    """Worker: everything about one event (fit, errors, variants,
    injection).  Pure computation on arrays handed in by the parent."""
    shape = ct.EventShape(**job["shape"])
    res = {"event_id": job["event_id"], "variants": [], "injection": [],
           "points": []}
    ap = job["apertures"][job["adopted"]]
    t, f, e, ex = (np.asarray(ap[k]) for k in ("t", "flux", "err", "ex"))

    def steps_of(a):
        """Re-pointing offsets of one aperture's points (None if none)."""
        st = np.asarray(a["steps"], dtype=float)
        return st if st.size and st.shape[1] > 0 else None

    steps = steps_of(ap)
    res["n_jumps"] = 0 if steps is None else int(steps.shape[1])

    def measure(t, f, e, ex, shape, clip=True, reg=None):
        """Blind fit -> one outlier pass -> refit.  Returns (fit, mask)."""
        fit = ct.fit_event(t, f, e, ex, shape, regressors=reg)
        mask = np.ones(len(t), dtype=bool)
        if fit["status"] not in (ct.STATUS_OK, ct.STATUS_ONE_SIDED):
            return fit, mask
        if clip:
            mask = ct.clip_outliers(f - fit["model"])
            if (~mask).any():
                fit = ct.fit_event(t[mask], f[mask], e[mask], ex[mask],
                                   shape, regressors=None if reg is None
                                   else reg[mask])
        return fit, mask

    fit, mask = measure(t, f, e, ex, shape, reg=steps)
    res["status"] = fit["status"]
    res["n_points"] = int(len(t))
    res["n_clipped"] = int((~mask).sum())
    if fit["status"] not in (ct.STATUS_OK, ct.STATUS_ONE_SIDED):
        return res
    # Depth sanity: is this THE event, on THE star?
    expect = job["depth_expected"]
    if job["min_depth"] is not None:
        if fit["depth"] < job["min_depth"]:
            res["status"] = ct.STATUS_SHALLOW
    elif not (ct.DEPTH_MIN_FRACTION * expect <= fit["depth"]
              <= ct.DEPTH_MAX_FRACTION * expect):
        res["status"] = ct.STATUS_SHALLOW
    tm_, fm, em, exm = t[mask], f[mask], e[mask], ex[mask]
    res.update({k: fit[k] for k in (
        "t0", "k", "a_over_r", "f1", "depth", "chi2", "dof", "rms",
        "n_before", "n_after", "t14_d", "t0_err_formal")})
    res["depth_expected"] = expect
    res["points"] = [(float(a), float(b), float(c), float(d), int(m))
                     for a, b, c, d, m in zip(
                         t, f, e, np.interp(t, tm_, fit["model"]), mask)]
    res["baseline"] = [float(v) for v in np.interp(t, tm_, fit["baseline"])]
    if res["status"] != ct.STATUS_OK:
        return res
    seed = ct.BOOT_SEED + job["seed"]
    sm = None if steps is None else steps[mask]
    boot = ct.bootstrap_t0(tm_, fm, em, exm, shape, fit, seed=seed,
                           regressors=sm)
    bead = ct.prayer_bead_t0(tm_, fm, em, exm, shape, fit, regressors=sm)
    res.update(sig_boot_s=boot["sigma"] * 86400.0,
               boot_bias_s=boot["bias"] * 86400.0,
               boot_n_ok=boot["n_ok"], block_points=boot["block_points"],
               sig_bead_s=bead["sigma"] * 86400.0)

    # ---- analysis variants: how far does each defensible choice move t0?
    def shift_of(vfit):
        return ((vfit["t0"] - fit["t0"]) * 86400.0
                if vfit["status"] == ct.STATUS_OK else None)

    variants = []
    for a_idx, other in enumerate(job["apertures"]):
        if a_idx == job["adopted"] or other is None:
            continue
        vt, vf, ve, vex = (np.asarray(other[k])
                           for k in ("t", "flux", "err", "ex"))
        vfit, _ = measure(vt, vf, ve, vex, shape, reg=steps_of(other))
        variants.append((f"aperture_{a_idx}", vfit))
    alt_order = 2 if shape.baseline_order == 1 else 1
    variants.append((f"baseline_order_{alt_order}", measure(
        t, f, e, ex, shape.copy(baseline_order=alt_order), reg=steps)[0]))
    variants.append(("ld_low", measure(t, f, e, ex, shape.copy(
        u1=max(0.0, shape.u1 - LD_VARIANT_STEP)), reg=steps)[0]))
    variants.append(("ld_high", measure(t, f, e, ex, shape.copy(
        u1=shape.u1 + LD_VARIANT_STEP), reg=steps)[0]))
    if "a_over_r" not in shape.free:
        variants.append(("duration_free", measure(t, f, e, ex, shape.copy(
            free=shape.free + ("a_over_r",)), reg=steps)[0]))
    else:
        variants.append(("k_x1.3", measure(t, f, e, ex, shape.copy(
            k=1.3 * shape.k), reg=steps)[0]))
    variants.append(("no_outlier_clip", measure(t, f, e, ex, shape,
                                                clip=False, reg=steps)[0]))
    if steps is not None:
        variants.append(("no_jump_offsets", measure(t, f, e, ex, shape)[0]))
    # Decorrelation against the two per-frame quantities that drive
    # ground-based systematics.  Not adopted (a regressor that happens to
    # change during ingress can absorb part of the event); recorded so the
    # reader sees how much the mid-time depends on it.
    for name, key in (("decorr_fwhm", "fwhm"), ("decorr_airmass", "airmass")):
        reg = np.asarray(ap[key], dtype=float)
        if np.isfinite(reg).all() and np.ptp(reg) > 0:
            if steps is not None:
                reg = np.column_stack([steps, reg])
            variants.append((name, measure(t, f, e, ex, shape, reg=reg)[0]))
    for name, vfit in variants:
        res["variants"].append((name, vfit["status"], shift_of(vfit)))
    # Model error: the root-mean-square mid-time shift over the variants
    # that returned a fit.  (Not the half-range: one bistable variant —
    # typically the smallest aperture on a drifting star — would then set
    # the error of an otherwise stable measurement all by itself; the rms
    # still counts it, in proportion.)
    shifts = np.array([s for _n, _st, s in res["variants"] if s is not None])
    res["sig_model_s"] = (float(np.sqrt(np.mean(shifts ** 2)))
                          if len(shifts) else 0.0)

    # ---- injection-recovery (standing rules 2 and 3)
    for j, shift in enumerate(ct.INJECT_SHIFTS_S):
        inj = ct.inject_and_recover(tm_, fm, em, exm, shape, fit, shift,
                                    ct.N_INJECT_DRAWS, seed + 1000 + j,
                                    regressors=sm)
        res["injection"].append((shift, inj["n_ok"], inj["n_draws"],
                                 inj["bias_s"], inj["scatter_s"]))
    return res


def anchored_refits(jobs: list, results: list) -> list:
    """Supplementary timing of ONE-SIDED binary eclipses with the shape
    anchored to the same star's complete eclipses.

    A free-shape symmetric template fitted to one flank is degenerate
    (duration against mid-time), which is why such events are refused.
    When the same star has at least :data:`clock_transits.ANCHOR_MIN_EVENTS`
    complete, timed eclipses on this telescope, their median duration and
    depth remove the degeneracy: the flank then fixes the mid-time.  The
    result REPLACES the refused row of that event with status
    ``ok_anchored_shape``; its model error is the rms mid-time shift when
    the anchored duration and depth are moved to the 16th and 84th
    percentiles of the complete eclipses and when the baseline order is
    changed.  The summary reports it and never averages it.
    """
    by_id = {r["event_id"]: r for r in results}
    target_of = lambda eid: eid.split("|", 1)[0]          # noqa: E731
    out = []
    for job in jobs:
        res = by_id.get(job["event_id"])
        if not res or res.get("status") != ct.STATUS_ONE_SIDED:
            continue
        if "a_over_r" not in job["shape"]["free"]:
            continue                                   # eclipses only
        good = [r for r in results if r.get("status") == ct.STATUS_OK
                and target_of(r["event_id"]) == target_of(job["event_id"])]
        if len(good) < ct.ANCHOR_MIN_EVENTS:
            continue
        a_r = np.array([g["a_over_r"] for g in good])
        f1s = np.array([g["f1"] for g in good])
        base = ct.EventShape(**{**job["shape"], "free": (),
                                "baseline_order": 1,
                                "a_over_r": float(np.median(a_r)),
                                "f1": float(np.median(f1s))})
        ap = job["apertures"][job["adopted"]]
        t, f, e, ex = (np.asarray(ap[k]) for k in ("t", "flux", "err", "ex"))
        fit = ct.fit_event(t, f, e, ex, base)
        if fit["status"] not in (ct.STATUS_OK, ct.STATUS_ONE_SIDED):
            continue
        mask = ct.clip_outliers(f - fit["model"])
        if (~mask).any():
            fit = ct.fit_event(t[mask], f[mask], e[mask], ex[mask], base)
            if fit["status"] not in (ct.STATUS_OK, ct.STATUS_ONE_SIDED):
                continue
        tm_, fm, em, exm = t[mask], f[mask], e[mask], ex[mask]
        seed = ct.BOOT_SEED + job["seed"]
        boot = ct.bootstrap_t0(tm_, fm, em, exm, base, fit, seed=seed)
        bead = ct.prayer_bead_t0(tm_, fm, em, exm, base, fit)
        variants = []
        lo_a, hi_a = np.percentile(a_r, [15.865, 84.135])
        lo_f, hi_f = np.percentile(f1s, [15.865, 84.135])
        for name, shp in (
                ("anchor_duration_lo", base.copy(a_over_r=float(hi_a))),
                ("anchor_duration_hi", base.copy(a_over_r=float(lo_a))),
                ("anchor_depth_lo", base.copy(f1=float(lo_f))),
                ("anchor_depth_hi", base.copy(f1=float(hi_f))),
                ("baseline_order_2", base.copy(baseline_order=2))):
            v = ct.fit_event(tm_, fm, em, exm, shp)
            ok = v["status"] in (ct.STATUS_OK, ct.STATUS_ONE_SIDED)
            variants.append((name, v["status"],
                             (v["t0"] - fit["t0"]) * 86400.0 if ok
                             else None))
        shifts = np.array([v[2] for v in variants if v[2] is not None])
        new = {"event_id": job["event_id"], "status": ct.STATUS_ANCHORED,
               "n_points": int(len(t)), "n_clipped": int((~mask).sum()),
               "variants": variants, "injection": [], "n_jumps": 0,
               "depth_expected": res.get("depth_expected"),
               "sig_boot_s": boot["sigma"] * 86400.0,
               "boot_bias_s": boot["bias"] * 86400.0,
               "boot_n_ok": boot["n_ok"],
               "block_points": boot["block_points"],
               "sig_bead_s": bead["sigma"] * 86400.0,
               "sig_model_s": (float(np.sqrt(np.mean(shifts ** 2)))
                               if len(shifts) else 0.0),
               "points": [(float(a), float(b_), float(c), float(d), int(m))
                          for a, b_, c, d, m in zip(
                              t, f, e, np.interp(t, tm_, fit["model"]),
                              mask)],
               "baseline": [float(v) for v in
                            np.interp(t, tm_, fit["baseline"])]}
        new.update({k: fit[k] for k in (
            "t0", "k", "a_over_r", "f1", "depth", "chi2", "dof", "rms",
            "n_before", "n_after", "t14_d", "t0_err_formal")})
        results.remove(res)                  # the refused row is replaced
        out.append(new)
    return out


def stage_fit(db: sqlite3.Connection, manifest: Path, targets: dict,
              workers: int, only: str | None) -> None:
    series = db.execute("""SELECT s.series_id, s.target, s.kind, s.readoutm,
                                  s.binning, p.veto_adu
                           FROM s3b_series s JOIN s3b_phot_series p
                             ON p.series_id = s.series_id
                           WHERE p.status = 'ok'
                           ORDER BY s.night, s.series_id""").fetchall()
    if only:
        series = [s for s in series if only in s[0]]
    jobs, lc_rows, series_rows, dead_events = [], [], [], []
    for sid, name, kind, readoutm, binning, veto in series:
        arr = load_series_arrays(db, sid)
        events = db.execute("""SELECT event_id, epoch, t_pred_bjd
                               FROM s3b_census_events
                               WHERE series_id = ? AND admitted = 1""",
                            (sid,)).fetchall()
        if arr is None:
            series_rows.append((sid, None, 0, 0, None, None, None,
                                "too few measured frames"))
            dead_events += [(ev[0], ct.STATUS_TOO_FEW) for ev in events]
            continue
        lcs = build_lightcurves(arr, veto, binning)
        ad = lcs["adopted"]
        n_meas = int(arr["good"][:, 0].sum())
        if ad is None:
            sat_frac = lcs["n_target_saturated"] / max(n_meas, 1)
            why = (ct.STATUS_SATURATED if sat_frac > 0.5
                   else ct.STATUS_NO_COMPS)
            series_rows.append((sid, None, len(arr["bjd"]),
                                lcs["n_target_saturated"], None, None, None,
                                why))
            dead_events += [(ev[0], why) for ev in events]
            continue
        a = lcs["apertures"][ad]
        series_rows.append((sid, ad, len(arr["bjd"]),
                            lcs["n_target_saturated"], int(a["use"].sum()),
                            int(a["keep"].sum()), float(a["sd_rms"]), "ok"))
        for ai, ap in enumerate(lcs["apertures"]):
            if ap is None:
                continue
            for i in np.flatnonzero(ap["use"]):
                lc_rows.append((sid, ai, arr["paths"][i],
                                float(arr["bjd"][i]),
                                float(arr["exptime_s"][i]),
                                float(ap["flux"][i]), float(ap["err"][i])))
        target = targets[name]
        shape = shape_of(target)
        for event_id, epoch, t_pred in events:
            aps = []
            for ap in lcs["apertures"]:
                if ap is None:
                    aps.append(None)
                    continue
                win = ap["use"] & event_window(kind, arr["bjd"], t_pred,
                                               shape)
                aps.append({"t": arr["bjd"][win], "flux": ap["flux"][win],
                            "err": ap["err"][win],
                            "ex": arr["exptime_s"][win] / 86400.0,
                            "fwhm": arr["fwhm"][win],
                            "airmass": arr["airmass"][win],
                            "steps": ct.jump_steps(arr["txy"][win, 0],
                                                   arr["txy"][win, 1])})
            if aps[ad] is None or len(aps[ad]["t"]) < ct.MIN_FIT_POINTS:
                dead_events.append((event_id, ct.STATUS_TOO_FEW))
                continue
            jobs.append({
                "event_id": event_id, "apertures": aps, "adopted": ad,
                "shape": dict(period=shape.period, a_over_r=shape.a_over_r,
                              k=shape.k, b=shape.b, u1=shape.u1,
                              u2=shape.u2, f1=shape.f1, free=shape.free,
                              baseline_order=shape.baseline_order),
                "depth_expected": shape.nominal_depth(),
                "min_depth": target["min_depth"],
                "seed": int(hashlib.sha1(event_id.encode()).hexdigest()[:6],
                            16)})
    log(f"fit: {len(jobs)} events to analyse "
        f"({len(dead_events)} already refused)")
    if workers > 1 and len(jobs) > 1:
        with mp.get_context("spawn").Pool(workers) as pool:
            results = list(pool.imap_unordered(analyse_event, jobs))
    else:
        results = [analyse_event(j) for j in jobs]
    results += anchored_refits(jobs, results)

    ev_rows, var_rows, inj_rows, pt_rows = [], [], [], []
    for event_id, why in dead_events:
        ev_rows.append((event_id, why) + (None,) * 23)
    for r in results:
        g = r.get
        sig_stat = None
        if g("sig_boot_s") is not None:
            cands = [v for v in (g("sig_boot_s"), g("sig_bead_s"))
                     if v is not None and np.isfinite(v)]
            sig_stat = max(cands) if cands else None
        status = r["status"]
        if status in (ct.STATUS_OK, ct.STATUS_ANCHORED) and sig_stat is None:
            status = "bootstrap_failed"      # no error bar -> no measurement
        ev_rows.append((
            r["event_id"], status, g("n_points"), g("n_clipped"),
            g("t0"), g("sig_boot_s"), g("sig_bead_s"),
            (g("t0_err_formal") * 86400.0
             if g("t0_err_formal") is not None else None),
            sig_stat, g("boot_bias_s"), g("sig_model_s"), g("depth"),
            g("depth_expected"), g("k"), g("a_over_r"), g("f1"),
            g("t14_d") * 1440.0 if g("t14_d") is not None else None,
            g("chi2"), g("dof"), g("rms"), g("n_before"), g("n_after"),
            g("block_points"), g("boot_n_ok"), g("n_jumps")))
        for name, status, shift in r["variants"]:
            var_rows.append((r["event_id"], name, status, shift))
        for row in r["injection"]:
            inj_rows.append((r["event_id"],) + tuple(row))
        base = r.get("baseline") or []
        for i, (t, f, e, m, used) in enumerate(r["points"]):
            pt_rows.append((r["event_id"], t, f, e, m,
                            base[i] if i < len(base) else None, used))
    replace_table(db, "s3b_series_lc", """CREATE TABLE {table} (
        series_id TEXT PRIMARY KEY, aperture_adopted INTEGER,
        n_frames_ok INTEGER, n_target_saturated INTEGER, n_used INTEGER,
        n_comps INTEGER, sd_rms REAL, status TEXT)""", series_rows, 8)
    replace_table(db, "s3b_lightcurve", """CREATE TABLE {table} (
        series_id TEXT, aperture INTEGER, path TEXT, bjd_tdb REAL,
        exptime_s REAL, flux REAL, err REAL,
        PRIMARY KEY (series_id, aperture, path))""", lc_rows, 7)
    replace_table(db, "s3b_fits", """CREATE TABLE {table} (
        event_id TEXT PRIMARY KEY, status TEXT, n_points INTEGER,
        n_clipped INTEGER, t0_bjd REAL, sig_boot_s REAL, sig_bead_s REAL,
        sig_formal_s REAL, sig_stat_s REAL, boot_bias_s REAL,
        sig_model_s REAL, depth REAL, depth_expected REAL, k REAL,
        a_over_r REAL, f1 REAL, t14_min REAL, chi2 REAL, dof INTEGER,
        rms REAL, n_before INTEGER, n_after INTEGER, block_points INTEGER,
        boot_n_ok INTEGER, n_jumps INTEGER)""", ev_rows, 25)
    replace_table(db, "s3b_variants", """CREATE TABLE {table} (
        event_id TEXT, variant TEXT, status TEXT, t0_shift_s REAL,
        PRIMARY KEY (event_id, variant))""", var_rows, 4)
    replace_table(db, "s3b_injection", """CREATE TABLE {table} (
        event_id TEXT, shift_s REAL, n_ok INTEGER, n_draws INTEGER,
        bias_s REAL, scatter_s REAL,
        PRIMARY KEY (event_id, shift_s))""", inj_rows, 6)
    replace_table(db, "s3b_fit_points", """CREATE TABLE {table} (
        event_id TEXT, bjd_tdb REAL, flux REAL, err REAL, model REAL,
        baseline REAL, used INTEGER)""", pt_rows, 7)
    n_ok = sum(r[1] in (ct.STATUS_OK, ct.STATUS_ANCHORED) for r in ev_rows)
    write_meta(db, {"fit_events": len(ev_rows), "fit_events_ok": n_ok,
                    "n_bootstrap": ct.N_BOOTSTRAP,
                    "boot_block_min": ct.BOOT_BLOCK_MIN})
    log(f"fit: {n_ok} of {len(ev_rows)} events timed")


# ---------------------------------------------------------------------------
# Stamp convention: does DATE-OBS mark the START or the MIDDLE of exposure?
# ---------------------------------------------------------------------------
#: Families (readout mode x camera era x acquisition build) in which the
#: convention is tested.  A cell = consecutive same-night frame pairs with
#: exposure e_a followed by e_b; its gap is the 10th percentile of the
#: stamp-to-stamp intervals (the shortest gaps are the machine cycle; the
#: long tail is waiting).  Symmetric differential, for each e_a < e_b:
#:     D = gap(e_a -> e_b) - gap(e_b -> e_a)
#: Every overhead tied to the sequence cancels.  START stamping gives
#: D = e_a - e_b; MID stamping gives D = 0; END stamping gives e_b - e_a.
#: The slope s of D on (e_a - e_b) is therefore 1 / 0 / -1.  Header
#: JD-HELIO cannot make this distinction (MaxIm computes it by ASSUMING a
#: start stamp), which is why S3's probe could not.
CONV_MIN_PAIRS = 3
CONV_MAX_DEXP_S = 300.0
#: ...and at least this large: the cycle gap jitters by a few seconds
#: (download, filter wheel), so a cell whose exposures differ by less than
#: this measures the jitter, not the convention (a 0.1 s vs 0.5 s cell
#: returns "ratios" of +-60).
CONV_MIN_DEXP_S = 8.0


def stamp_convention(con_manifest) -> list[tuple]:
    """Rows of ``s3b_stamp_convention`` (see the constants above)."""
    rows = con_manifest.execute(f"""
        SELECT night, era_id, readoutm, filter, exptime, jd, swcreate
        FROM frames WHERE {SCIENCE_WHERE} AND jd IS NOT NULL
          AND exptime > 0 AND tree = 'rawimage'
        ORDER BY night, jd""").fetchall()
    cells: dict = {}
    for a, b in zip(rows, rows[1:]):
        if a[0] != b[0] or a[1] != b[1] or a[2] != b[2] or a[6] != b[6]:
            continue
        gap = (b[5] - a[5]) * 86400.0
        if gap <= 0 or gap > max(a[4], b[4]) + 180.0:
            continue
        key = (a[2] or "", int(a[1]), a[6] or "", round(a[4], 1),
               round(b[4], 1), a[3] != b[3])
        cells.setdefault(key, []).append(gap)
    fam: dict = {}
    for k, g in cells.items():
        if k[3] >= k[4]:
            continue
        k2 = k[:3] + (k[4], k[3], k[5])
        if k2 not in cells or min(len(g), len(cells[k2])) < CONV_MIN_PAIRS:
            continue
        x = k[3] - k[4]
        if not CONV_MIN_DEXP_S <= abs(x) <= CONV_MAX_DEXP_S:
            continue
        d = float(np.percentile(g, 10) - np.percentile(cells[k2], 10))
        fam.setdefault(k[:3], []).append((x, d))
    out = []
    for (mode, era, sw), v in sorted(fam.items()):
        x = np.array([t[0] for t in v])
        d = np.array([t[1] for t in v])
        r = d / x
        lo, med, hi = np.percentile(r, [25, 50, 75])
        if len(v) < 5:
            verdict = "undetermined (fewer than 5 cells)"
        elif abs(med - 1.0) < 0.25 and lo > 0.5:
            verdict = "START"
        elif abs(med) < 0.25 and hi < 0.5:
            verdict = "MID"
        else:
            verdict = "undetermined (cells disagree)"
        out.append((mode, era, sw, ct.clock_era("", era) if False else
                    None, len(v), float(med), float(lo), float(hi),
                    float(np.mean(r > 0.5)), verdict))
    return out


# ---------------------------------------------------------------------------
# Stage: summary
# ---------------------------------------------------------------------------
def stage_summary(db: sqlite3.Connection, manifest: Path, cv_db: Path,
                  targets: dict) -> None:
    """O - C per event -> per target -> per clock era, and the corollaries."""
    rows = db.execute("""
        SELECT c.event_id, c.series_id, c.target, c.epoch, c.t_pred_bjd,
               c.sig_pred_s, s.night, s.filter, s.readoutm, s.era_id,
               s.tree, s.clock_era, s.kind, f.t0_bjd, f.sig_stat_s,
               f.sig_model_s, f.status
        FROM s3b_census_events c
        JOIN s3b_series s ON s.series_id = c.series_id
        JOIN s3b_fits f ON f.event_id = c.event_id
        WHERE f.status IN (?, ?) ORDER BY f.t0_bjd""",
        (ct.STATUS_OK, ct.STATUS_ANCHORED)).fetchall()
    oc_rows = []
    for (event_id, sid, name, epoch, t_pred, sig_pred, night, filt, readoutm,
         era_id, tree, cera, kind, t0, sig_stat, sig_model, fstatus) in rows:
        t = targets[name]
        oc = (t0 - t_pred) * 86400.0
        alt = t.get("alt")
        alt_diff = None
        if alt:
            e_alt = ct.nearest_epoch(t0, alt["t0"], alt["period"])
            tp_alt, _ = ct.predict_event(alt["t0"], alt["sig_t0"],
                                         alt["period"], alt["sig_p"], e_alt,
                                         alt.get("quad", 0.0),
                                         alt.get("sig_quad", 0.0))
            alt_diff = (tp_alt - t_pred) * 86400.0
        # Astrophysical systematic of the standard itself: the published
        # timing-variation amplitude, or (for eclipsing binaries) the
        # disagreement of two published ephemerides, whichever is larger.
        sys_s = t["sys_floor_s"]
        if kind == "eclipse" and alt_diff is not None:
            sys_s = max(sys_s, abs(alt_diff))
        sig_meas = math.hypot(sig_stat, sig_model)
        if fstatus == ct.STATUS_ANCHORED:
            grade = ct.GRADE_ANCHORED
        elif sig_meas > ct.TIMING_GRADE_MAX_SIGMA_S:
            grade = ct.GRADE_LOW
        else:
            grade = ct.GRADE_TIMING
        oc_rows.append((event_id, sid, name, kind, int(epoch), night, filt,
                        readoutm, int(era_id), tree, cera, t0, t_pred, oc,
                        sig_stat, sig_model, sig_meas, sig_pred, sys_s,
                        math.sqrt(sig_meas ** 2 + sig_pred ** 2
                                  + sys_s ** 2),
                        alt_diff, 0, grade))
    # Primary flag: one row per physical photometric series.  A raw frame
    # and its server-reduced twin are the SAME photons, so only one of a
    # twin pair may enter a mean: the raw one when it was timed, else the
    # reduced one.
    oc_rows = [list(r) for r in oc_rows]
    by_key = {}
    for r in oc_rows:
        by_key.setdefault((r[2], r[4], r[5], r[6], r[7]), []).append(r)
    for members in by_key.values():
        # Only timing-grade rows compete to represent a series; among
        # them the raw one wins, and its reduced twin is labelled.
        cand = [m for m in members if m[22] == ct.GRADE_TIMING]
        raw = [m for m in cand if m[9] != "reduced"]
        chosen = raw if raw else cand[:1]
        for m in cand:
            if any(m is c for c in chosen):
                m[21] = 1
            else:
                m[22] = ct.GRADE_TWIN
    replace_table(db, "s3b_oc", """CREATE TABLE {table} (
        event_id TEXT PRIMARY KEY, series_id TEXT, target TEXT, kind TEXT,
        epoch INTEGER, night TEXT, filter TEXT, readoutm TEXT,
        era_id INTEGER, tree TEXT, clock_era TEXT, t0_bjd REAL,
        t_pred_bjd REAL, oc_s REAL, sig_stat_s REAL, sig_model_s REAL,
        sig_meas_s REAL, sig_eph_s REAL, sig_sys_s REAL, sig_total_s REAL,
        alt_minus_primary_s REAL, is_primary INTEGER, grade TEXT)""",
                  [tuple(r) for r in oc_rows], 23)

    # ---- per (clock era, target), then per clock era -------------------
    prim = [r for r in oc_rows if r[21] == 1]
    tgt_rows, era_rows = [], []
    for cera, label, _f, _l in ct.CLOCK_ERAS:
        in_era = [r for r in prim if r[10] == cera]
        t_oc, t_sig = [], []
        for name in sorted({r[2] for r in in_era}):
            ev = [r for r in in_era if r[2] == name]
            c = ct.combine_oc([r[13] for r in ev], [r[16] for r in ev])
            sig_eph = float(np.mean([r[17] for r in ev]))
            sig_sys = float(max(r[18] for r in ev))
            sig_t = math.sqrt(c["sigma_adopted"] ** 2 + sig_eph ** 2
                              + sig_sys ** 2)
            tgt_rows.append((cera, name, ev[0][3], c["n"],
                             len({(r[4]) for r in ev}), c["wmean"],
                             c["wmean_err"], c["chi2"], c["dof"],
                             c["scatter_err"], c["sigma_adopted"], sig_eph,
                             sig_sys, sig_t, min(r[5] for r in ev),
                             max(r[5] for r in ev)))
            t_oc.append(c["wmean"])
            t_sig.append(sig_t)
        # For information only: the inverse-variance mean of EVERY timed
        # event of the era (low-precision ones included, twins and
        # anchored fits excluded), with total errors and no scatter term.
        every = [r for r in oc_rows if r[10] == cera
                 and r[22] in (ct.GRADE_TIMING, ct.GRADE_LOW)]
        n_low = sum(r[22] == ct.GRADE_LOW for r in every)
        c_all = ct.combine_oc([r[13] for r in every],
                              [r[19] for r in every]) if every else {}
        info = (n_low, c_all.get("wmean"), c_all.get("wmean_err"))
        if not t_oc:
            era_rows.append((cera, label, 0, 0, None, None, None, None,
                             None, None, None,
                             ("NO TIMING-GRADE EVENT IN ERA" if every
                              else "NO CLOCK TARGET IN ERA"),
                             None, None, None, None) + info)
            continue
        c = ct.combine_oc(t_oc, t_sig)
        # Sensitivity to the CV-sized offset: the injected +-1,065 s
        # shifts, pooled over the era's events (signed bias).
        ids = [r[0] for r in in_era]
        marks = ",".join("?" * len(ids))
        inj = db.execute(f"""
            SELECT sum(n_ok), sum(n_draws),
                   sum(bias_s * n_ok) / nullif(sum(n_ok), 0)
            FROM s3b_injection WHERE abs(shift_s) = ?
              AND event_id IN ({marks})""",
            [ct.CV_OFFSET_S] + ids).fetchone()
        era_rows.append((cera, label, len(t_oc), len(in_era), c["wmean"],
                         c["wmean_err"], c["chi2"], c["dof"], c["mean"],
                         c["scatter_err"], c["sigma_adopted"],
                         ct.era_verdict(c["wmean"], c["sigma_adopted"]),
                         3.0 * c["sigma_adopted"], inj[0], inj[1], inj[2])
                        + info)
    replace_table(db, "s3b_target_oc", """CREATE TABLE {table} (
        clock_era TEXT, target TEXT, kind TEXT, n_series INTEGER,
        n_epochs INTEGER, oc_s REAL, oc_formal_err_s REAL, chi2 REAL,
        dof INTEGER, oc_scatter_err_s REAL, sig_meas_adopted_s REAL,
        sig_eph_s REAL, sig_sys_s REAL, sig_total_s REAL,
        first_night TEXT, last_night TEXT,
        PRIMARY KEY (clock_era, target))""", tgt_rows, 16)
    replace_table(db, "s3b_era", """CREATE TABLE {table} (
        clock_era TEXT PRIMARY KEY, label TEXT, n_targets INTEGER,
        n_series INTEGER, oc_s REAL, oc_formal_err_s REAL, chi2 REAL,
        dof INTEGER, oc_unweighted_s REAL, oc_scatter_err_s REAL,
        sig_adopted_s REAL, verdict TEXT, detectable_3sig_s REAL,
        inj1065_n_ok INTEGER, inj1065_n_draws INTEGER,
        inj1065_bias_s REAL, n_low_precision INTEGER, oc_all_s REAL,
        oc_all_err_s REAL)""", era_rows, 19)

    # ---- pair comparisons: StackPro vs plain; raw vs reduced ------------
    pair_rows = []
    for (name, epoch, night, filt, _mode), members in by_key.items():
        raws = [m for m in members if m[9] != "reduced"
                and m[22] != ct.GRADE_ANCHORED]
        reds = [m for m in members if m[9] == "reduced"
                and m[22] != ct.GRADE_ANCHORED]
        for a in raws:
            for b in reds:
                pair_rows.append(("raw_vs_reduced", name, int(epoch), night,
                                  a[0], b[0], b[13] - a[13],
                                  math.hypot(a[16], b[16])))
    by_ev = {}
    for r in (r for r in oc_rows if r[22] != ct.GRADE_ANCHORED
              and r[9] != "reduced"):
        by_ev.setdefault((r[2], r[4], r[5], r[6]), []).append(r)
    for (name, epoch, night, filt), members in by_ev.items():
        sp = [m for m in members if tm.is_stackpro(m[7])]
        pl = [m for m in members if not tm.is_stackpro(m[7])]
        for a in pl:
            for b in sp:
                pair_rows.append(("stackpro_minus_plain", name, int(epoch),
                                  night, a[0], b[0], b[13] - a[13],
                                  math.hypot(a[16], b[16])))
    replace_table(db, "s3b_pairs", """CREATE TABLE {table} (
        kind TEXT, target TEXT, epoch INTEGER, night TEXT,
        event_a TEXT, event_b TEXT, delta_s REAL, sig_s REAL)""",
                  pair_rows, 8)

    # ---- where the time stamps came from (audit, per clock era) ---------
    # Every time on this page is header JD (= DATE-OBS, the exposure
    # start) + EXPTIME/2.  File names are never read for time: from 2024
    # the acquisition scripts embed the SCHEDULED start in the name, which
    # leads the real start by minutes.  Both statements are measured here
    # on the very frames that were timed.
    import re
    name_re = re.compile(r"(\d{4}-\d{2}-\d{2})T(\d{2})-(\d{2})-(\d{2})")
    used = db.execute("""SELECT f.path, f.jd_utc_start, s.clock_era
                         FROM s3b_frames f JOIN s3b_series s
                           USING (series_id)""").fetchall() \
        if db.execute("SELECT count(*) FROM sqlite_master WHERE "
                      "name = 's3b_frames'").fetchone()[0] else []
    audit: dict = {}
    with open_ro(manifest) as con:
        for path, jd_used, cera in used:
            row = con.execute("SELECT date_obs, jd, swcreate FROM frames "
                              "WHERE path = ?", (path,)).fetchone()
            if not row:
                continue
            a = audit.setdefault((cera, row[2] or ""), {"n": 0, "d": [],
                                                        "lead": []})
            a["n"] += 1
            jd_do = tm.parse_date_obs(row[0])
            if jd_do is not None:
                a["d"].append(abs(jd_used - jd_do) * 86400.0)
            m = name_re.search(path.rsplit("/", 1)[-1])
            if m and jd_do is not None:
                jd_name = tm.parse_date_obs(
                    f"{m.group(1)}T{m.group(2)}:{m.group(3)}:{m.group(4)}")
                a["lead"].append((jd_do - jd_name) * 86400.0)
    replace_table(db, "s3b_stamp_audit", """CREATE TABLE {table} (
        clock_era TEXT, swcreate TEXT, n_frames INTEGER,
        max_abs_used_minus_dateobs_s REAL, n_with_name_time INTEGER,
        name_lead_median_s REAL, name_lead_max_s REAL)""",
                  [(k[0], k[1], a["n"],
                    max(a["d"]) if a["d"] else None, len(a["lead"]),
                    float(np.median(a["lead"])) if a["lead"] else None,
                    max(a["lead"]) if a["lead"] else None)
                   for k, a in sorted(audit.items())], 7)

    # ---- start-vs-mid stamping, from cadence -----------------------------
    with open_ro(manifest) as con:
        conv = stamp_convention(con)
    replace_table(db, "s3b_stamp_convention", """CREATE TABLE {table} (
        readoutm TEXT, era_id INTEGER, swcreate TEXT, unused TEXT,
        n_cells INTEGER, slope_median REAL, slope_q25 REAL, slope_q75 REAL,
        frac_cells_start_like REAL, verdict TEXT)""", conv, 10)
    # What a MID stamp would do to the StackPro events timed here: S3's
    # policy adds EXPTIME/2 to every stamp, so a mid-stamped frame's time
    # is late by EXPTIME/2 and its O - C too high by the same amount.
    sp = [(r[0], r[2], r[5], r[13], r[22]) for r in oc_rows
          if tm.is_stackpro(r[7])]
    sp_rows = []
    with open_ro(manifest) as con:
        for eid, name, night, oc, grade in sp:
            sid = eid.split("#")[0]
            e = db.execute("SELECT avg(exptime_s) FROM s3b_frames WHERE "
                           "series_id = ?", (sid,)).fetchone()[0]
            sp_rows.append((eid, name, night, grade, e, oc, oc - e / 2.0))
    replace_table(db, "s3b_stackpro_mid", """CREATE TABLE {table} (
        event_id TEXT PRIMARY KEY, target TEXT, night TEXT, grade TEXT,
        exptime_s REAL, oc_s REAL, oc_if_mid_stamp_s REAL)""", sp_rows, 7)

    # ---- the stated bounds ----------------------------------------------
    replace_table(db, "s3b_bounds", """CREATE TABLE {table} (
        term TEXT, applies_to TEXT, bound_s REAL, sign TEXT, basis TEXT)""",
                  list(ct.TIMING_BOUNDS), 5)

    # ---- S3's AG LMi -294 s residual, explained -------------------------
    legacy = {}
    with open_ro(manifest) as con:
        g = con.execute("""SELECT o_minus_c_s, o_minus_c_err_s,
                                  clock_bound_s, n_points
                           FROM s3_clock_eclipses WHERE tag = 'global'
                             AND status = 'ok'""").fetchone()
        meta = dict(con.execute("SELECT key, value FROM s3_build_meta"))
    era_b = next((r for r in era_rows if r[0] == "B"), None)
    if g and era_b and era_b[4] is not None:
        cycles = float(meta["clock_mean_cycle"])
        p_vsx = float(meta["vsx_period_d"])
        p_gaia = float(meta["gaia_period_d"])
        sig_gaia = float(meta["gaia_period_err_d"])
        # What is left of the AG LMi residual once the measured clock
        # offset of its own era is removed belongs to the star's ephemeris.
        resid = g[0] - era_b[4]
        resid_err = math.hypot(g[1], era_b[10])
        dp = resid / 86400.0 / cycles
        legacy = {
            "agl_oc_s": g[0], "agl_oc_err_s": g[1], "agl_bound_s": g[2],
            "agl_n_points": g[3], "agl_cycles": cycles,
            "agl_nights": meta.get("clock_nights_gated", ""),
            "vsx_period_d": p_vsx, "gaia_period_d": p_gaia,
            "gaia_period_err_d": sig_gaia,
            "gaia_predicted_oc_s": float(meta["gaia_predicted_oc_s"]),
            "gaia_envelope_s": float(meta["gaia_oc_envelope_s"]),
            "vsx_epoch_quant_s": 43.2,
            "era_b_clock_oc_s": era_b[4], "era_b_clock_err_s": era_b[10],
            "agl_ephemeris_part_s": resid,
            "agl_ephemeris_part_err_s": resid_err,
            "agl_implied_period_d": p_vsx + dp,
            "agl_implied_period_err_d": resid_err / 86400.0 / cycles,
            "agl_implied_minus_gaia_sigma":
                (p_vsx + dp - p_gaia) / sig_gaia,
        }
    replace_table(db, "s3b_legacy", """CREATE TABLE {table} (
        key TEXT PRIMARY KEY, value TEXT)""",
                  [(k, str(v)) for k, v in legacy.items()], 2)

    # ---- the shared CV offset, era by era -------------------------------
    cv_rows = []
    if cv_db.exists():
        with open_ro(cv_db) as cv:
            head = {r[0]: r[1:] for r in cv.execute(
                """SELECT target_key, oc_mean_s, oc_rms_s, n_epochs
                   FROM p3_cycle_count WHERE oc_mean_s IS NOT NULL""")}
            nights = cv.execute("""SELECT target_key, night, era_id
                                   FROM p3_oc_night
                                   GROUP BY 1, 2, 3""").fetchall()
        era_by = {r[0]: r for r in era_rows}
        for tkey, (oc_mean, oc_rms, n_ep) in sorted(head.items()):
            per = {}
            for tk, night, era_id in nights:
                if tk == tkey:
                    per.setdefault(ct.clock_era(night, era_id),
                                   []).append(night)
            for cera, nl in sorted(per.items()):
                e = era_by.get(cera)
                cv_rows.append((tkey, cera, len(nl), min(nl), max(nl),
                                oc_mean, oc_rms, n_ep,
                                e[4] if e else None, e[10] if e else None,
                                e[11] if e else None))
    replace_table(db, "s3b_cv", """CREATE TABLE {table} (
        target_key TEXT, clock_era TEXT, n_nights INTEGER,
        first_night TEXT, last_night TEXT, cv_offset_s REAL,
        cv_offset_rms_s REAL, cv_n_epochs INTEGER, clock_oc_s REAL,
        clock_sig_s REAL, clock_verdict TEXT,
        PRIMARY KEY (target_key, clock_era))""", cv_rows, 11)
    write_meta(db, {"summary_events_ok": len(oc_rows),
                    "summary_primary": len(prim),
                    "accept_oc_s": ct.ACCEPT_OC_S,
                    "cv_offset_s": ct.CV_OFFSET_S})
    for r in era_rows:
        if r[4] is not None:
            log(f"era {r[0]} ({r[1]}): O-C = {r[4]:+.1f} +- {r[10]:.1f} s "
                f"from {r[2]} targets / {r[3]} series -> {r[11]}")
        else:
            log(f"era {r[0]} ({r[1]}): {r[11]}")


# ---------------------------------------------------------------------------
# Stage: report (figures + page), rendered from the products database only
# ---------------------------------------------------------------------------
def render_report(db_path: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from macro_core import plotstyle as ps
    from macro_core.report_s0 import _figure, esc, fmt, q, q1, table

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = open_ro(db_path)
    meta = dict(q(con, "SELECT key, value FROM s3b_build_meta"))
    accept = float(meta["accept_oc_s"])
    cvoff = float(meta["cv_offset_s"])

    def f1(x, d=1, signed=False):
        if x is None or (isinstance(x, float) and not math.isfinite(x)):
            return "&mdash;"
        return f"{x:+.{d}f}" if signed else f"{x:.{d}f}"

    # One colour + marker per standard class (colour never alone).
    KIND = {"transit": dict(color=ps.ACCENT, marker="o", label="transit"),
            "eclipse": dict(color=ps.WARN, marker="s",
                            label="binary eclipse")}

    # ---- Figure 1: O - C of every timed event against date -------------
    def fig_oc() -> str:
        rows = q(con, """SELECT t0_bjd, oc_s, sig_total_s, kind, clock_era,
                                is_primary, grade FROM s3b_oc
                         ORDER BY t0_bjd""")
        eras = q(con, """SELECT e.clock_era, e.oc_s, e.sig_adopted_s,
                                min(o.t0_bjd), max(o.t0_bjd)
                         FROM s3b_era e JOIN s3b_oc o
                           ON o.clock_era = e.clock_era
                         WHERE e.oc_s IS NOT NULL GROUP BY 1""")
        with ps.context():
            fig, ax = plt.subplots(figsize=(9.4, 4.6))
            ax.axhspan(-accept, accept, color=ps.tint(ps.GOOD, 0.85),
                       zorder=0)
            ax.axhline(0, color=ps.MUTED, lw=0.8)
            ax.axhline(cvoff, color=ps.BAD, lw=1.2, ls="--")
            ax.annotate(f"the shared CV edge offset, +{cvoff:,.0f} s",
                        (0.01, cvoff), xycoords=("axes fraction", "data"),
                        va="bottom", color=ps.BAD, fontsize=8)
            for kind, kw in KIND.items():
                sel = [r for r in rows if r[3] == kind and r[5]]
                if not sel:
                    continue
                yr = [2000.0 + (r[0] - 2451544.5) / 365.25 for r in sel]
                ax.errorbar(yr, [r[1] for r in sel],
                            yerr=[r[2] for r in sel], ls="none",
                            marker=kw["marker"], ms=5, color=kw["color"],
                            mec=ps.INK, mew=0.4, elinewidth=0.8, zorder=3,
                            label=f"{kw['label']}, timing grade "
                                  f"({len(sel)})")
            for grade, mk, lab in ((ct.GRADE_LOW, "o", "low precision "
                                    "(not in means)"),
                                   (ct.GRADE_ANCHORED, "D", "one-sided, "
                                    "anchored shape (not in means)")):
                sel = [r for r in rows if r[6] == grade]
                if not sel:
                    continue
                yr = [2000.0 + (r[0] - 2451544.5) / 365.25 for r in sel]
                ax.errorbar(yr, [r[1] for r in sel],
                            yerr=[r[2] for r in sel], ls="none", marker=mk,
                            ms=5, mfc="none", mec=ps.MUTED, ecolor=ps.MUTED,
                            elinewidth=0.6, zorder=2,
                            label=f"{lab} ({len(sel)})")
            for cera, oc, sig, t_lo, t_hi in eras:
                y0 = 2000.0 + (t_lo - 2451544.5) / 365.25 - 0.03
                y1 = 2000.0 + (t_hi - 2451544.5) / 365.25 + 0.03
                ax.fill_between([y0, y1], oc - sig, oc + sig,
                                color=ps.INK, alpha=0.25, lw=0)
                ax.plot([y0, y1], [oc, oc], color=ps.INK, lw=1.6)
                ax.annotate(cera, ((y0 + y1) / 2, -accept * 1.9),
                            ha="center", fontsize=9, weight="bold")
            ax.set_yscale("symlog", linthresh=150, linscale=1.6)
            ax.set_ylim(-3000, 20000)
            ticks = [-1000, -300, -120, -60, 0, 60, 120, 300, 1000, 10000]
            ax.set_yticks(ticks)
            ax.set_yticklabels([f"{v:,}" for v in ticks])
            ax.set_xlabel("date (year)")
            ax.set_ylabel("O - C (s)   [linear inside +-150 s]")
            ax.set_title("Observed minus predicted mid-time of every timed "
                         "event; black bars = clock-era means")
            ax.legend(loc="upper right", ncol=2, fontsize=7)
            fig.tight_layout()
            fig.savefig(FIG_DIR / "s3b_oc_timeline.png", dpi=ps.WEB_DPI)
            plt.close(fig)
        return "figures/s3b/s3b_oc_timeline.png"

    # ---- Figure 2..: light-curve galleries, one per clock era ----------
    def fig_gallery(cera: str) -> str | None:
        evs = q(con, """SELECT o.event_id, o.target, o.night, o.filter,
                               o.readoutm, o.tree, o.t0_bjd, o.t_pred_bjd,
                               o.oc_s, o.sig_total_s, o.sig_eph_s
                        FROM s3b_oc o WHERE o.clock_era = ?
                        ORDER BY o.t0_bjd""", (cera,))
        if not evs:
            return None
        ncol = 4
        nrow = int(math.ceil(len(evs) / ncol))
        with ps.context():
            fig, axes = plt.subplots(nrow, ncol, squeeze=False,
                                     figsize=(11.0, 2.3 * nrow + 0.5))
            for ax in axes.ravel():
                ax.set_visible(False)
            for ax, ev in zip(axes.ravel(), evs):
                ax.set_visible(True)
                pts = np.array(q(con, """SELECT bjd_tdb, flux, model,
                                                baseline, used
                                         FROM s3b_fit_points
                                         WHERE event_id = ?
                                         ORDER BY bjd_tdb""", (ev[0],)),
                               dtype=float)
                tref = ev[7]
                x = (pts[:, 0] - tref) * 1440.0
                used = pts[:, 4] > 0
                base = np.where(np.isfinite(pts[:, 3]), pts[:, 3], 1.0)
                ax.plot(x[used], pts[used, 1] / base[used], ls="none",
                        marker=".", ms=3, color=ps.MUTED)
                ax.plot(x[~used], pts[~used, 1] / base[~used], ls="none",
                        marker="x", ms=3, color=ps.BAD)
                ax.plot(x[used], pts[used, 2] / base[used], color=ps.ACCENT,
                        lw=1.2)
                band = ev[10] / 60.0
                ax.axvspan(-band, band, color=ps.tint(ps.GOOD, 0.6),
                           zorder=0)
                ax.axvline(0.0, color=ps.GOOD, lw=1.0)
                ax.axvline((ev[6] - tref) * 1440.0, color=ps.WARN, lw=1.0,
                           ls="--")
                mode = "SP" if tm.is_stackpro(ev[4]) else ""
                ax.set_title(f"{ev[1]} {ev[2]} {ev[3]}{' ' + mode if mode else ''}"
                             f"{' red.' if ev[5] == 'reduced' else ''}\n"
                             f"O-C {ev[8]:+.0f} +- {ev[9]:.0f} s",
                             fontsize=7.5)
                ax.tick_params(labelsize=7)
            for ax in axes[-1]:
                ax.set_xlabel("minutes from predicted mid-time", fontsize=8)
            for ax in axes[:, 0]:
                ax.set_ylabel("relative flux", fontsize=8)
            fig.suptitle(f"Clock era {cera} — {ct.clock_era_label(cera)}: "
                         "detrended light curves, model (blue), predicted "
                         "(green) and fitted (orange dashed) mid-times",
                         fontsize=9)
            fig.tight_layout()
            name = f"s3b_gallery_{cera}.png"
            fig.savefig(FIG_DIR / name, dpi=ps.WEB_DPI)
            plt.close(fig)
        return f"figures/s3b/{name}"

    # ---- Figure: injection-recovery -------------------------------------
    def fig_injection() -> str:
        rows = q(con, """SELECT i.shift_s, i.bias_s, i.scatter_s, i.n_ok,
                                i.n_draws, o.kind
                         FROM s3b_injection i JOIN s3b_oc o
                           ON o.event_id = i.event_id""")
        shifts = sorted({r[0] for r in rows})
        with ps.context():
            fig, (ax, ax2) = plt.subplots(
                1, 2, figsize=(9.4, 3.8),
                gridspec_kw={"width_ratios": [1.6, 1]})
            rng = np.random.default_rng(3)
            for kind, kw in KIND.items():
                sel = [r for r in rows if r[5] == kind and r[3] > 0]
                if not sel:
                    continue
                xs = [shifts.index(r[0]) + rng.uniform(-0.25, 0.25)
                      for r in sel]
                ax.errorbar(xs, [r[1] for r in sel],
                            yerr=[r[2] for r in sel], ls="none",
                            marker=kw["marker"], ms=3.5, color=kw["color"],
                            elinewidth=0.5, alpha=0.8, label=kw["label"])
            ax.axhline(0, color=ps.MUTED, lw=0.8)
            ax.set_xticks(range(len(shifts)))
            ax.set_xticklabels([f"{s:+.0f}" for s in shifts])
            ax.set_xlabel("injected clock offset (s)")
            ax.set_ylabel("recovered - injected (s), signed")
            ax.set_title("Signed recovery error per event")
            ax.legend(loc="upper left")
            frac = []
            for s in shifts:
                sel = [r for r in rows if r[0] == s]
                frac.append(sum(r[3] for r in sel)
                            / max(1, sum(r[4] for r in sel)))
            ax2.bar(range(len(shifts)), frac, color=ps.ACCENT)
            ax2.set_xticks(range(len(shifts)))
            ax2.set_xticklabels([f"{s:+.0f}" for s in shifts], rotation=45)
            ax2.set_ylim(0, 1.05)
            ax2.set_xlabel("injected clock offset (s)")
            ax2.set_ylabel("fraction of draws recovered")
            ax2.set_title("Would it have been seen?")
            fig.tight_layout()
            fig.savefig(FIG_DIR / "s3b_injection.png", dpi=ps.WEB_DPI)
            plt.close(fig)
        return "figures/s3b/s3b_injection.png"

    # ---- Figure: analysis variants --------------------------------------
    def fig_variants() -> str:
        rows = q(con, """SELECT v.variant, v.t0_shift_s FROM s3b_variants v
                         JOIN s3b_oc o ON o.event_id = v.event_id
                         WHERE v.t0_shift_s IS NOT NULL""")
        names = sorted({r[0] for r in rows})
        with ps.context():
            fig, ax = plt.subplots(figsize=(9.4, 3.6))
            rng = np.random.default_rng(5)
            for i, n in enumerate(names):
                ys = np.array([r[1] for r in rows if r[0] == n])
                ax.plot(i + rng.uniform(-0.2, 0.2, len(ys)), ys, ls="none",
                        marker="o", ms=3, color=ps.ACCENT, alpha=0.7)
                ax.plot([i - 0.3, i + 0.3], [np.median(ys)] * 2,
                        color=ps.INK, lw=1.6)
            ax.axhline(0, color=ps.MUTED, lw=0.8)
            ax.axhspan(-accept, accept, color=ps.tint(ps.GOOD, 0.85),
                       zorder=0)
            ax.set_xticks(range(len(names)))
            ax.set_xticklabels(names, rotation=20, ha="right")
            ax.set_ylabel("mid-time shift vs adopted fit (s)")
            ax.set_yscale("symlog", linthresh=60)
            ax.set_title("Every defensible analysis choice, re-fitted: how "
                         "far the mid-time moves (bar = median)")
            fig.tight_layout()
            fig.savefig(FIG_DIR / "s3b_variants.png", dpi=ps.WEB_DPI)
            plt.close(fig)
        return "figures/s3b/s3b_variants.png"

    # ---- Figure: AG LMi explained ---------------------------------------
    def fig_legacy(leg: dict) -> str | None:
        if not leg:
            return None
        cyc = float(leg["agl_cycles"])
        p_vsx, p_gaia = float(leg["vsx_period_d"]), float(leg["gaia_period_d"])
        sig = float(leg["gaia_period_err_d"])
        with ps.context():
            fig, ax = plt.subplots(figsize=(7.4, 3.6))
            p = np.linspace(p_gaia - 1.2 * sig, p_vsx + 1.2 * sig, 200)
            ax.plot((p - p_vsx) * 1e6, (p - p_vsx) * cyc * 86400.0,
                    color=ps.INK, lw=1.2,
                    label="O - C a true period P would produce")
            ax.axvspan((p_gaia - sig - p_vsx) * 1e6,
                       (p_gaia + sig - p_vsx) * 1e6,
                       color=ps.tint(ps.ACCENT, 0.8),
                       label="Gaia DR3 period +- 1 sigma")
            ax.axvline(0, color=ps.MUTED, lw=0.8)
            oc, err = float(leg["agl_oc_s"]), float(leg["agl_oc_err_s"])
            ax.axhspan(oc - err, oc + err, color=ps.tint(ps.WARN, 0.6))
            ax.axhline(oc, color=ps.WARN, lw=1.2,
                       label=f"S3 measurement {oc:+.0f} +- {err:.0f} s")
            ax.set_xlabel("true period minus VSX period (1e-6 d)")
            ax.set_ylabel("AG LMi O - C at cycle "
                          f"{cyc:,.0f} (s)")
            ax.legend(loc="lower right", fontsize=7.5)
            ax.set_title("The -294 s AG LMi residual is a period error of "
                         "two parts in a million")
            fig.tight_layout()
            fig.savefig(FIG_DIR / "s3b_aglmi.png", dpi=ps.WEB_DPI)
            plt.close(fig)
        return "figures/s3b/s3b_aglmi.png"

    # =====================================================================
    # Sections
    # =====================================================================
    era_rows = q(con, """SELECT clock_era, label, n_targets, n_series, oc_s,
                                oc_formal_err_s, chi2, dof, oc_unweighted_s,
                                oc_scatter_err_s, sig_adopted_s, verdict,
                                detectable_3sig_s, inj1065_n_ok,
                                inj1065_n_draws, inj1065_bias_s,
                                n_low_precision, oc_all_s, oc_all_err_s
                         FROM s3b_era ORDER BY clock_era""")
    n_pass = sum(1 for r in era_rows if (r[11] or "").startswith("PASS")
                 and "interval" not in r[11])
    n_meas = sum(1 for r in era_rows if r[4] is not None)

    def chi_txt(chi2, dof):
        if chi2 is None or not dof:
            return "&mdash;"
        nu = chi2 / dof
        flag = " &#9888;" if (nu < 0.5 or nu > 2.0) else ""
        return f"{chi2:.1f}/{dof} = {nu:.2f}{flag}"

    era_tbl = table(
        ["era", "camera / software", "targets", "series",
         "O&minus;C (s)", "formal &plusmn;", "&chi;&sup2;/dof",
         "scatter &plusmn;", "adopted &plusmn;", "3&sigma; reach (s)",
         f"&plusmn;{cvoff:,.0f} s injections recovered", "signed bias (s)",
         "verdict", "low-precision events (not in mean)",
         "information only: mean of ALL events (s)"],
        [[r[0], esc(r[1]), fmt(r[2]), fmt(r[3]), f1(r[4], 1, True),
          f1(r[5]), chi_txt(r[6], r[7]), f1(r[9]), f1(r[10]), f1(r[12], 0),
          (f"{r[13]}/{r[14]}" if r[14] else "&mdash;"), f1(r[15], 1, True),
          esc(r[11]), fmt(r[16]),
          (f"{r[17]:+.0f} &plusmn; {r[18]:.0f}" if r[17] is not None
           else "&mdash;")] for r in era_rows],
        row_classes=[None if (r[11] or "").startswith("PASS") else "warn"
                     for r in era_rows])

    tgt = q(con, """SELECT clock_era, target, kind, n_series, n_epochs,
                           oc_s, oc_formal_err_s, chi2, dof,
                           oc_scatter_err_s, sig_eph_s, sig_sys_s,
                           sig_total_s, first_night, last_night
                    FROM s3b_target_oc ORDER BY clock_era, target""")
    tgt_tbl = table(
        ["era", "standard", "kind", "series", "epochs", "nights",
         "O&minus;C (s)", "formal &plusmn;", "&chi;&sup2;/dof",
         "scatter &plusmn;", "ephemeris &plusmn;", "standard&rsquo;s own "
         "systematic", "total &plusmn;"],
        [[r[0], esc(r[1]), r[2], fmt(r[3]), fmt(r[4]),
          esc(r[13] if r[13] == r[14] else f"{r[13]} &rarr; {r[14]}"),
          f1(r[5], 1, True), f1(r[6]), chi_txt(r[7], r[8]), f1(r[9]),
          f1(r[10]), f1(r[11]), f1(r[12])] for r in tgt])

    ev = q(con, """SELECT o.clock_era, o.target, o.night, o.filter,
                          o.readoutm, o.tree, o.epoch, o.oc_s, o.sig_stat_s,
                          o.sig_model_s, o.sig_eph_s, o.sig_sys_s,
                          o.sig_total_s, f.chi2, f.dof, f.sig_boot_s,
                          f.sig_bead_s, f.sig_formal_s, f.depth,
                          f.depth_expected, f.n_points, f.n_clipped,
                          o.alt_minus_primary_s, o.is_primary,
                          l.n_target_saturated, f.n_jumps, o.grade
                   FROM s3b_oc o JOIN s3b_fits f ON f.event_id = o.event_id
                   JOIN s3b_series_lc l ON l.series_id = o.series_id
                   ORDER BY o.t0_bjd""")
    ev_tbl = table(
        ["era", "standard", "night", "filter", "mode", "tree", "cycle",
         "O&minus;C (s)", "bootstrap", "prayer bead", "formal", "model",
         "ephem.", "total &plusmn;", "&chi;&sup2;/dof", "depth / expected",
         "points (clipped)", "re-pointing offsets",
         "target frames over veto",
         "2nd ephemeris &minus; 1st (s)", "in mean"],
        [[r[0], esc(r[1]), r[2], esc(r[3]), esc(r[4] or "blank"), r[5],
          fmt(r[6]), f1(r[7], 1, True), f1(r[15]), f1(r[16]), f1(r[17]),
          f1(r[9]), f1(r[10]), f1(r[12]), chi_txt(r[13], r[14]),
          f"{r[18]:.4f} / {r[19]:.4f}", f"{r[20]} ({r[21]})", fmt(r[25]),
          fmt(r[24]),
          f1(r[22], 1, True),
          "yes" if r[23] else esc(r[26].replace("_", " "))] for r in ev],
        row_classes=[None if r[23] else "warn" for r in ev])

    # census accounting
    n_series = q1(con, "SELECT count(*) FROM s3b_series")
    n_series_adm = q1(con, "SELECT count(*) FROM s3b_series WHERE admitted=1")
    n_ev = q1(con, "SELECT count(*) FROM s3b_census_events")
    n_ev_adm = q1(con, "SELECT count(*) FROM s3b_census_events "
                       "WHERE admitted = 1")
    n_timed = q1(con, "SELECT count(*) FROM s3b_oc")
    n_prim = q1(con, "SELECT count(*) FROM s3b_oc WHERE is_primary = 1")
    refused = q(con, """SELECT CASE WHEN reason LIKE 'event not covered%'
                                    THEN 'event outside the run even allowing a clock error'
                                    WHEN reason LIKE 'ephemeris%'
                                    THEN 'ephemeris prediction looser than '
                                         || ? || ' s'
                                    ELSE reason END, count(*)
                        FROM s3b_census_events WHERE admitted = 0
                        GROUP BY 1""", (int(ct.EPH_SIGMA_MAX_S),))
    fate = q(con, """SELECT coalesce(f.status, 'photometry refused: '
                                     || p.status), count(*)
                     FROM s3b_census_events c
                     LEFT JOIN s3b_fits f ON f.event_id = c.event_id
                     LEFT JOIN s3b_phot_series p ON p.series_id = c.series_id
                     WHERE c.admitted = 1 GROUP BY 1 ORDER BY 2 DESC""")
    frames_n = q1(con, "SELECT count(*) FROM s3b_frames")
    frames_ok = q1(con, "SELECT count(*) FROM s3b_frames WHERE status='ok'")
    eph_tbl = table(
        ["standard", "role", "T0 (BJD_TDB)", "&plusmn; (d)", "P (d)",
         "&plusmn; (d)", "reference"],
        [[esc(r[0]), r[1], f"{r[2]:.6f}", f"{r[3]:.1e}", f"{r[4]:.9f}",
          f"{r[5]:.1e}", esc(r[6])] for r in q(con, """
            SELECT target, role, t0_bjd, sig_t0_d, period_d, sig_period_d,
                   reference FROM s3b_ephemeris
            WHERE target IN (SELECT target FROM s3b_oc)
            ORDER BY target, role DESC""")])
    recipe_tbl = table(
        ["era", "reduction applied", "series", "flat age (d)",
         "dark age (d)"],
        [[r[0], esc(r[1]), fmt(r[2]),
          (f"{r[3]:.0f}&ndash;{r[4]:.0f}" if r[3] is not None
           else "&mdash;"),
          (f"{r[5]:.0f}&ndash;{r[6]:.0f}" if r[5] is not None
           else "&mdash;")] for r in q(con, """
            SELECT s.clock_era, p.recipe, count(*), min(p.flat_age_d),
                   max(p.flat_age_d), min(p.dark_age_d), max(p.dark_age_d)
            FROM s3b_phot_series p JOIN s3b_series s USING (series_id)
            WHERE p.status = 'ok' GROUP BY 1, 2 ORDER BY 1, 2""")])
    phot_refused = q(con, """SELECT status, count(*) FROM s3b_phot_series
                             WHERE status != 'ok' GROUP BY 1""")
    refused_tbl = table(["candidate", "what", "why it is not a clock"],
                        [[esc(r[0]), esc(r[1]), esc(r[2])] for r in q(
                            con, "SELECT * FROM s3b_refused")])

    pairs = q(con, "SELECT kind, target, night, delta_s, sig_s FROM s3b_pairs "
                   "ORDER BY kind, night")
    sp = [(p[3], p[4]) for p in pairs if p[0] == "stackpro_minus_plain"]
    sp_comb = ct.combine_oc([p[0] for p in sp], [p[1] for p in sp]) \
        if sp else {"n": 0}
    rr = [(p[3], p[4]) for p in pairs if p[0] == "raw_vs_reduced"]
    pair_tbl = table(["comparison", "standard", "night", "&Delta; (s)",
                      "&plusmn; (s)"],
                     [[esc(p[0].replace("_", " ")), esc(p[1]), p[2],
                       f1(p[3], 1, True), f1(p[4])] for p in pairs])
    bounds_tbl = table(["term", "applies to", "bound (s)", "sign", "basis"],
                       [[esc(r[0]), esc(r[1]), f"{r[2]:.2f}", esc(r[3]),
                         esc(r[4])] for r in q(con,
                                               "SELECT * FROM s3b_bounds")])
    stamp_tbl = table(
        ["era", "acquisition software", "frames timed",
         "max |time used &minus; DATE-OBS| (s)", "frames with a time in "
         "the file name", "DATE-OBS minus file-name time: median (s)",
         "max (s)"],
        [[r[0], esc(r[1]), fmt(r[2]), f1(r[3], 3), fmt(r[4]), f1(r[5], 0),
          f1(r[6], 0)] for r in q(con, "SELECT * FROM s3b_stamp_audit")])
    conv_tbl = table(
        ["readout mode", "era", "software", "cells", "median ratio",
         "IQR", "fraction start-like", "verdict"],
        [[esc(r[0]), r[1], esc(r[2]), fmt(r[4]), f"{r[5]:.2f}",
          f"{r[6]:.2f}&ndash;{r[7]:.2f}", f"{r[8]:.2f}", esc(r[9])]
         for r in q(con, "SELECT * FROM s3b_stamp_convention")],
        row_classes=[None if r[0] == "START" else "warn" for r in q(
            con, "SELECT verdict FROM s3b_stamp_convention")])
    spmid_tbl = table(
        ["standard", "night", "grade", "EXPTIME (s)", "O&minus;C as built "
         "(s)", "O&minus;C if mid-stamped (s)"],
        [[esc(r[1]), r[2], esc(r[3]), f1(r[4], 0), f1(r[5], 1, True),
          f1(r[6], 1, True)] for r in q(con, "SELECT * FROM "
                                             "s3b_stackpro_mid")])
    leg = dict(q(con, "SELECT key, value FROM s3b_legacy"))
    cv = q(con, """SELECT target_key, clock_era, n_nights, first_night,
                          last_night, cv_offset_s, cv_offset_rms_s,
                          clock_oc_s, clock_sig_s, clock_verdict
                   FROM s3b_cv ORDER BY target_key, clock_era""")
    cv_tbl = table(
        ["CV target", "clock era", "nights timed", "span",
         "published edge offset (s)", "clock O&minus;C in that era (s)",
         "&plusmn;", "clock verdict", "offset / clock error"],
        [[esc(r[0]), r[1], fmt(r[2]), f"{r[3]} &rarr; {r[4]}",
          f"{r[5]:+.0f} (rms {r[6]:.0f})", f1(r[7], 1, True), f1(r[8]),
          esc(r[9] or "no clock standard in this era"),
          (f"{abs(r[5]) / max(r[8], 1e-9):.0f}&sigma; away"
           if r[8] else "&mdash;")] for r in cv],
        row_classes=[None if r[7] is not None else "warn" for r in cv])

    inj = q(con, """SELECT i.shift_s, sum(i.n_ok), sum(i.n_draws),
                           sum(i.bias_s * i.n_ok) / nullif(sum(i.n_ok), 0),
                           max(abs(i.bias_s))
                    FROM s3b_injection i JOIN s3b_oc o
                      ON o.event_id = i.event_id
                    GROUP BY 1 ORDER BY 1""")
    inj_tbl = table(["injected offset (s)", "draws recovered",
                     "mean signed bias (s)", "largest |bias| of any event (s)"],
                    [[f"{r[0]:+.0f}", f"{r[1]}/{r[2]}", f1(r[3], 1, True),
                      f1(r[4])] for r in inj])
    var = q(con, """SELECT v.variant, count(*), sum(v.t0_shift_s IS NOT NULL),
                           avg(v.t0_shift_s), max(abs(v.t0_shift_s))
                    FROM s3b_variants v JOIN s3b_oc o
                      ON o.event_id = v.event_id GROUP BY 1 ORDER BY 1""")
    var_tbl = table(["variant", "events", "re-fitted ok",
                     "mean signed shift (s)", "largest |shift| (s)"],
                    [[esc(r[0]), fmt(r[1]), fmt(r[2]), f1(r[3], 1, True),
                      f1(r[4])] for r in var])
    chi_all = q(con, """SELECT sum(f.chi2 / f.dof < 0.5),
                               sum(f.chi2 / f.dof > 2.0), count(*),
                               min(f.chi2 / f.dof), max(f.chi2 / f.dof)
                        FROM s3b_fits f JOIN s3b_oc o
                          ON o.event_id = f.event_id""")[0]

    galleries = "".join(
        f'<div class="grid">{_figure(src, f"Clock era {r[0]}: every timed event.  Flux is divided by the fitted baseline; crosses are clipped points; the green band is the literature prediction with its 1-sigma width, the dashed orange line the fitted mid-time.")}</div>'
        for r in era_rows if (src := fig_gallery(r[0])))
    src_oc = fig_oc()
    src_inj = fig_injection()
    src_var = fig_variants()
    src_leg = fig_legacy(leg)

    worst = max((abs(r[4]) + 2 * r[10] for r in era_rows
                 if r[4] is not None), default=float("nan"))
    eras_without = [f"{r[0]} ({esc(r[1])})" for r in era_rows
                    if r[4] is None]
    headline = (f"{n_pass} of {n_meas} measured clock eras pass "
                f"|O&minus;C| &lt; {accept:.0f} s at 2&sigma;; the largest "
                f"|O&minus;C| + 2&sigma; of any measured era is "
                f"{worst:.0f} s")

    legacy_html = ""
    if leg:
        legacy_html = f"""
<section id="aglmi"><h2>5 &middot; S3&rsquo;s AG LMi residual, printed and
explained</h2><div class="q">
<p class="sub"><b>The number (RF.M4).</b>  S3&rsquo;s one eclipsing-binary
check, AG&nbsp;LMi on {esc(leg['agl_nights'])}, gave
O&minus;C&nbsp;=&nbsp;{float(leg['agl_oc_s']):+.0f}&nbsp;&plusmn;&nbsp;{float(leg['agl_oc_err_s']):.0f}&nbsp;s
from {esc(leg['agl_n_points'])} points against the VSX ephemeris, and a
clock bound of {float(leg['agl_bound_s']):,.0f}&nbsp;s.  A five-sigma,
five-minute residual on a timing paper&rsquo;s only external check has to
be explained, not described as &ldquo;weak&rdquo;.</p>
<p class="sub"><b>The explanation.</b>  That eclipse lies
{float(leg['agl_cycles']):,.0f} cycles after the VSX epoch.  VSX quotes the
period as {float(leg['vsx_period_d']):.7f}&nbsp;d with no uncertainty;
Gaia&nbsp;DR3 measures
{float(leg['gaia_period_d']):.7f}&nbsp;&plusmn;&nbsp;{float(leg['gaia_period_err_d']):.1e}&nbsp;d
for the same star, which by itself predicts an O&minus;C of
{float(leg['gaia_predicted_oc_s']):+.0f}&nbsp;s with a
&plusmn;{float(leg['gaia_envelope_s']):,.0f}&nbsp;s envelope.  The clock
of that same camera era (B) is now measured independently on this page:
{float(leg['era_b_clock_oc_s']):+.1f}&nbsp;&plusmn;&nbsp;{float(leg['era_b_clock_err_s']):.1f}&nbsp;s.
Subtracting it leaves
{float(leg['agl_ephemeris_part_s']):+.0f}&nbsp;&plusmn;&nbsp;{float(leg['agl_ephemeris_part_err_s']):.0f}&nbsp;s
that belongs to AG&nbsp;LMi&rsquo;s ephemeris, i.e. a true period of
{float(leg['agl_implied_period_d']):.7f}&nbsp;&plusmn;&nbsp;{float(leg['agl_implied_period_err_d']):.1e}&nbsp;d
&mdash; {abs(float(leg['agl_implied_minus_gaia_sigma'])):.2f}&nbsp;&sigma;
from the Gaia period and between the two catalogue values.  (The VSX epoch
is also quoted to 0.001&nbsp;d, &plusmn;{float(leg['vsx_epoch_quant_s']):.0f}&nbsp;s,
which this conversion ignores; it does not change the conclusion.)  The
residual is a stale survey ephemeris, not a clock error, and the S3
statistical error of {float(leg['agl_oc_err_s']):.0f}&nbsp;s was never the
relevant uncertainty.</p>
{f'<div class="grid">{_figure(src_leg, "The straight line is the O-C any true period would produce at the cycle of the AG LMi eclipse.  The S3 measurement (orange) crosses it inside the Gaia period interval (blue).")}</div>' if src_leg else ''}
</div></section>"""

    cv_have = [r for r in cv if r[7] is not None]
    cv_none = [r for r in cv if r[7] is None]
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S3b &mdash; The Absolute Clock, from archived transits</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S3b &mdash; The Absolute Clock</h1>
  <p>{fmt(n_prim)} independent series ({fmt(n_timed)} with twins) timed
  on literature standards across {n_meas} clock eras &middot; {headline}
  &middot; built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))}, commit
  <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="s3_timing.html">S3, the time axis this checks</a></p>
</header>

<nav>
  <a href="#question">1 Question</a> &middot;
  <a href="#census">2 Census</a> &middot;
  <a href="#events">3 Events</a> &middot;
  <a href="#eras">4 Per-era verdict</a> &middot;
  <a href="#aglmi">5 AG LMi</a> &middot;
  <a href="#cv">6 The 1,065 s offset</a> &middot;
  <a href="#tests">7 Bias tests</a> &middot;
  <a href="#bounds">8 Stated bounds</a> &middot;
  <a href="#limits">9 Limits</a>
</nav>

<section id="question"><h2>1 &middot; Is the observatory clock right?</h2>
<div class="q">
<p class="sub">S3 proved what a header stamp <i>means</i>.  It could not
prove the stamp is <i>true</i>: both header clock cards come from the same
PC, and its one astrophysical check carried a
&plusmn;{float(leg.get('gaia_envelope_s', 0)):,.0f}&nbsp;s ephemeris
envelope.  The plan review (OA.E6, DS.F10, RF.M4) noted that two polars
show the same ~{cvoff:,.0f}&nbsp;s offset from their catalogue
ephemerides &mdash; exactly what a clock error would do.</p>
<p class="sub"><b>Method.</b>  The archive holds transits of well-timed
exoplanets and eclipses of post-common-envelope binaries.  For each fully
covered event: differential aperture photometry on the raw frames;
mid-exposure BJD_TDB recomputed at the star&rsquo;s own coordinates; a
limb-darkened occultation model located by a <i>blind</i> grid over the
whole run (the prediction never tells the fit where to look) and refined by
least squares; a {meta.get('n_bootstrap')}-draw moving-block residual
bootstrap ({meta.get('boot_block_min')}&nbsp;min blocks) and a prayer-bead
permutation, the larger of which is adopted; then O&minus;C against the
literature ephemeris propagated with its own uncertainty.
O&minus;C&nbsp;&gt;&nbsp;0 means our stamps are late.</p>
<div class="grid">{_figure(src_oc, f"Every timed event.  The green band is the acceptance window of +-{accept:.0f} s; the dashed red line is where a clock error equal to the shared CV offset would put the points.  The axis is linear inside +-150 s and logarithmic outside.  Black bars are the per-era means with their adopted 1-sigma errors; letters name the clock eras of section 4.")}</div>
</div></section>

<section id="census"><h2>2 &middot; What the archive offers</h2>
<div class="q">
<p class="sub">Every canonical science series of at least
{ct.MIN_SERIES_FRAMES} frames and {ct.MIN_SERIES_SPAN_H}&nbsp;h whose
field centre lies within {CONE_RADIUS_ARCMIN:.0f}&prime; of a clock
standard was examined: {fmt(n_series)} series, {fmt(n_ev)} predicted
events.  An event is admitted when the prediction puts it inside the
run with {ct.MIN_MARGIN_H}&nbsp;h to spare on each side
({ct.MIN_MARGIN_ECLIPSE_H}&nbsp;h for the deep, sharp binary eclipses)
<i>or misses that by no more than {ct.CLOCK_SEARCH_PAD_H * 3600:.0f}&nbsp;s</i>
&mdash; the selection must not assume the clock it is testing, so
borderline events are measured and the two-sided-coverage decision is
taken after the blind fit &mdash; and its ephemeris predicts
it to better than {ct.EPH_SIGMA_MAX_S:.0f}&nbsp;s:
{fmt(n_ev_adm)} events in {fmt(n_series_adm)} series.  Refusals:
{'; '.join(f"{esc(r[0])}: {r[1]}" for r in refused) or 'none'}.</p>
<p class="sub">Of the admitted events:
{'; '.join(f"<b>{esc(r[0])}</b>: {r[1]}" for r in fate)}.
{fmt(frames_ok)} of {fmt(frames_n)} frames were registered and measured.
A status other than <code>ok</code> is a refusal by a rule fixed before any
O&minus;C was computed (one-sided coverage after the fit, a depth outside
{ct.DEPTH_MIN_FRACTION}&ndash;{ct.DEPTH_MAX_FRACTION}&times; the catalogue
depth, the target over the native-pixel saturation veto on most frames, or
no stable comparison stars).</p>
<h3>How the frames were reduced</h3>
<p class="sub">Raw frames are calibrated with the archive&rsquo;s own master
dark (same readout mode, same exposure time, never scaled) and master flat
(same mode, same filter), the nearest in time of each &mdash; the same rule
the CV campaign uses.  The telescope was not guided for most of these
runs: stars walk tens of pixels and are re-centred in jumps, and without a
flat every jump is a step in the light curve as large as a transit.  The
flats are months old in places (ages below); what they fail to remove is
handled in the fit, where every re-pointing found in the target&rsquo;s
centroid track (a jump of more than {ct.JUMP_PX:.0f}&nbsp;px between
consecutive frames) gets its own offset.  Saturation is judged on the raw
pixels.  {('Series refused at this step: ' + '; '.join(f"{esc(r[0])} ({r[1]})" for r in phot_refused) + '.') if phot_refused else ''}</p>
{recipe_tbl}
<h3>The standards and their ephemerides</h3>
{eph_tbl}
<h3>Candidates refused before any photometry</h3>
{refused_tbl}
</div></section>

<section id="events"><h2>3 &middot; Every timed event</h2><div class="q">
{galleries}
{ev_tbl}
<p class="sub">Errors are in seconds.  <i>bootstrap</i> and <i>prayer
bead</i> are the two red-noise estimates (the larger is the statistical
error); <i>formal</i> is the least-squares value, shown only so the reader
can see how much correlated noise inflates it; <i>model</i> is the
root-mean-square mid-time shift over the analysis variants of section 7.
&chi;&sup2; uses empirical point errors from successive differences, so a
ratio above one measures correlated noise the white-noise estimate cannot
see: {chi_all[1]} of {chi_all[2]} fits have &chi;&sup2;/dof &gt; 2 and
{chi_all[0]} have &chi;&sup2;/dof &lt; 0.5 (range {chi_all[3]:.2f} to
{chi_all[4]:.2f}; flagged &#9888; in the tables).  No error on this page is
rescaled by a &chi;&sup2;: the bootstrap draws from the actual residuals.
The last column says whether a row enters the era means.  <i>twin</i>:
the server-reduced copy of a raw series already counted.  <i>low
precision</i>: the event&rsquo;s own measurement error (statistical and
model in quadrature) exceeds {ct.TIMING_GRADE_MAX_SIGMA_S:.0f}&nbsp;s, so
it cannot test a {accept:.0f}&nbsp;s criterion; the gate looks at the
error, never at the O&minus;C.  <i>anchored</i>: a binary eclipse of
which only one flank was observed, timed with the duration and depth
fixed to the median of the same star&rsquo;s complete eclipses &mdash; a
supplementary check, never averaged.</p>
</div></section>

<section id="eras"><h2>4 &middot; The verdict per clock era</h2>
<div class="q">
<p class="sub">A clock era is one camera under one acquisition-software
build.  Events of one standard are combined first (they share an
ephemeris error, which is then added once, with the standard&rsquo;s own
timing systematic); standards are then combined per era.  Both the formal
and the scatter-based error are shown and the larger is adopted.</p>
{era_tbl}
<h3>Per standard</h3>
{tgt_tbl}
<div class="decision"><b>{headline}.</b>
{('Eras with no usable standard, where the clock is NOT directly verified: ' + '; '.join(eras_without) + '.') if eras_without else 'Every era has a standard.'}</div>
</div></section>
{legacy_html}

<section id="cv"><h2>6 &middot; What this means for the shared
~{cvoff:,.0f} s CV offset</h2><div class="q">
<p class="sub">ST&nbsp;LMi and EU&nbsp;UMa sit about {cvoff:,.0f}&nbsp;s
from their catalogue ephemerides.  If that were the clock, every standard
on this page observed in the same era would sit at +{cvoff:,.0f}&nbsp;s
too (the dashed line of the first figure).</p>
{cv_tbl}
<div class="decision">
{f"In {len(cv_have)} of {len(cv)} (target, era) blocks a clock standard exists in the same era, and in each the clock O&minus;C is tens of seconds or less while the edge offset is over a thousand: <b>the shared offset is not an observatory clock error there</b>." if cv_have else ""}
{f" In {len(cv_none)} block(s) ({'; '.join(f'{esc(r[0])} era {r[1]}, {r[3]} to {r[4]}' for r in cv_none)}) no standard was observed, so the clock is not verified directly; the statement there rests on the same offset being measured in the verified eras on either side." if cv_none else ""}
The offset is therefore a property of the comparison (what feature each
catalogue epoch marks, and how stale each catalogue period is), which is
for the CV analysis to settle; it must not be described as, or corrected
as, a clock term.</div>
</div></section>

<section id="tests"><h2>7 &middot; Would a clock error have been seen,
and can the analysis move the answer?</h2><div class="q">
<p class="sub"><b>Injection.</b>  For every event the fitted model was
removed, re-injected displaced by a known offset, block-resampled noise
added, and the whole blind measurement repeated
({ct.N_INJECT_DRAWS} draws per offset).  The bias is signed, per the
standing rule.</p>
{inj_tbl}
<div class="grid">{_figure(src_inj, "Left: recovered minus injected offset for every event and draw-set (error bar = scatter over draws).  Right: fraction of draws in which the displaced event was still recovered as a two-sided fit; a displaced event that leaves the observing window is, correctly, not recovered.")}</div>
<p class="sub"><b>Variants.</b>  Each event was re-measured with the
other two apertures, the other baseline order, limb darkening moved by
&plusmn;{LD_VARIANT_STEP}, the duration (or the radius ratio) freed, and
without outlier rejection.</p>
{var_tbl}
<div class="grid">{_figure(src_var, "Mid-time shift of every event under every variant, relative to the adopted fit.  The rms shift per event is carried as its model error.")}</div>
<h3>Pairs that should agree</h3>
{pair_tbl if pairs else '<p class="sub">No same-event pairs were timed.</p>'}
<p class="sub">{f"StackPro minus plain High Gain on the same transit: {sp_comb['wmean']:+.1f} &plusmn; {sp_comb['sigma_adopted']:.1f} s from {sp_comb['n']} pair(s) &mdash; a direct measurement of the StackPro mid-time policy that S3 could only bound from cadence." if sp_comb.get('n') else "No transit was timed in both StackPro and plain High Gain, so the StackPro mid-time policy stays at S3&rsquo;s cadence bound."}
{f" Raw versus server-reduced copies of the same frames ({len(rr)} pair(s)) are listed above; they share photons, so their difference tests the reduction, not the clock." if rr else ""}</p>
</div></section>

<section id="bounds"><h2>8 &middot; Timing terms no event here can
measure: stated bounds</h2><div class="q">
<p class="sub">DE.F7 asked for the rolling-shutter and mechanical-shutter
terms to be stated.  They are sub-second, far below what any O&minus;C
above resolves, and are therefore <i>bounds from vendor figures, not
measurements</i>:</p>
{bounds_tbl}
<p class="sub">None of these terms is applied as a correction.  They add
to the S3 error budget for any claim of sub-second absolute timing; they
are negligible against every timing claim currently made by a MACRO
paper.  Also for the record (OA.E6): the two header clock cards DATE-OBS
and TELUT are written by the same PC and are not independent.</p>
<h3>Does DATE-OBS mark the start or the middle of the exposure?</h3>
<p class="sub">S3 read the convention off JD-HELIO, but MaxIm computes
JD-HELIO by <i>assuming</i> a start stamp, so that probe cannot tell.  A
test that can: for consecutive frames with exposures e<sub>a</sub> and
e<sub>b</sub>, compare the shortest stamp-to-stamp gap a&rarr;b with the
gap b&rarr;a.  Every overhead cancels; a start stamp leaves
e<sub>a</sub>&minus;e<sub>b</sub>, a mid stamp leaves zero.  The ratio is
1 for START and 0 for MID, per (readout mode, era, software build):</p>
{conv_tbl}
<p class="sub">Where a family is MID, S3&rsquo;s start + EXPTIME/2 puts
the time late by EXPTIME/2.  The StackPro events timed on this page, with
their O&minus;C under a mid stamp:</p>
{spmid_tbl}
<h3>Where the time stamps came from</h3>
<p class="sub">Every time on this page is the header JD (identical to
DATE-OBS, the exposure start) plus half the exposure.  File names are
never read for time: the scripted eras embed the <i>scheduled</i> start in
the name, which precedes the real start.  Measured on the frames timed
here:</p>
{stamp_tbl}
</div></section>

<section id="limits"><h2>9 &middot; What this does not show</h2>
<div class="q"><ul>
<li>An era&rsquo;s verdict certifies the nights that have a standard.  A
PC clock can step between nights; the evidence against that is the
agreement of standards spread across each era, not a proof.</li>
<li>The light curves are not publication photometry.  Unguided drift,
months-old flats and passing cloud leave correlated noise of several
millimagnitudes; that is why the errors come from a block bootstrap of the
real residuals and why the analysis variants are carried as a model error.
An event whose red noise is large gets a large error and little weight; it
is not removed.</li>
<li>Planet ephemerides come from one living catalogue (ExoClock); the
independent NASA-archive ephemeris of each planet is used only to show
how far a second prediction lies from the first (event table).</li>
<li>The eclipsing binaries have their own timing variations; those are
carried as a systematic per star, not removed.</li>
</ul></div></section>

<footer>Generated by <code>pipeline/scripts/build_s3b_clock_transits.py</code>
from <code>products/clock/clock_transits.sqlite</code> — every number on
this page is the result of a query; none is typed by hand.</footer>
</body></html>"""
    HTML_PATH.write_text(html, encoding="utf-8")
    import re
    for src in re.findall(r'<img src="([^"]+)"', html):
        p = DOCS_DIR / src
        if not p.exists() or p.stat().st_size == 0:
            raise RuntimeError(f"report references missing figure: {src}")
    con.close()
    return HTML_PATH


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
STAGES = ("ephemerides", "census", "photometry", "fit", "summary", "report")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--stage", choices=STAGES + ("all",), default="all")
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--cv-db", type=Path, default=DEFAULT_CV_DB)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--refresh-ephemerides", action="store_true",
                    help="re-fetch the literature ephemerides (network)")
    ap.add_argument("--only", default=None,
                    help="restrict photometry/fit to series ids containing "
                         "this text (debugging)")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    run = STAGES if args.stage == "all" else (args.stage,)
    if args.refresh_ephemerides:
        refresh_ephemerides(args.manifest)
    db = open_db(args.db)
    db.execute("""CREATE TABLE IF NOT EXISTS s3b_ref_bright (
        series_id TEXT PRIMARY KEY, xy_json TEXT)""")
    ephemeris = tm.resolve_ephemeris()
    targets = stage_ephemerides(db)        # cheap; every stage needs them
    if "census" in run:
        stage_census(db, args.manifest, targets, ephemeris)
    if "photometry" in run:
        stage_photometry(db, args.manifest, args.archive, targets,
                         ephemeris, args.workers, args.only)
    if "fit" in run:
        stage_fit(db, args.manifest, targets, args.workers, args.only)
    if "summary" in run:
        stage_summary(db, args.manifest, args.cv_db, targets)
    db.close()
    if "report" in run:
        log(f"report: {render_report(args.db)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
