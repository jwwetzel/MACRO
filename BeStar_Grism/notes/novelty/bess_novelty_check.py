#!/usr/bin/env python
"""BE-S-1a / BE-N1 — the BeSS novelty check for the Be-star grism campaign.

WHAT THIS SCRIPT DOES
---------------------
The committee review of 2026-10-03 (SYNTHESIS U9; journal-editor "Do
BE-S-1a-bess first and alone"; hostile-referee "Step -1 gates before any
pipeline work") ruled that no Be-star pipeline effort may be spent until one
question is answered from OUTSIDE our own data:

    Which of the stars we observed were classical Be stars showing H-alpha
    emission DURING our seasons, as recorded by somebody else?

This script answers that from two external sources and one internal one:

    manifest (read-only)   every grism target, nights per season and grism
    SIMBAD / Sesame        identity, coordinates, spectral type, object type
    BeSS (basebe.obspm.fr) every public spectrum of each target, through the
                           BeSS Simple Spectral Access (SSA) service; the
                           H-alpha spectra that bracket our seasons are
                           downloaded and MEASURED (equivalent width and
                           emission-peak height), never read off a web page.

Every number in ``REPORT.md`` and in the tables beside this file is emitted
by this script from ``novelty.sqlite``.  Nothing is typed.

WHY THE EMISSION STATE IS MEASURED AND NOT LOOKED UP
----------------------------------------------------
BeSS is a catalogue of stars that have EVER been classified Be.  Membership
says nothing about whether the disk existed in 2025: Be stars lose and
rebuild their disks on timescales of months to decades.  The only defensible
"verified active" is a spectrum, taken by an independent observer inside or
next to our season, in which H-alpha stands above the continuum.  So the
script downloads those spectra and applies one fixed measurement to all of
them (see ``measure_halpha``).

SUBCOMMANDS (run in this order; each is idempotent and cached)
--------------------------------------------------------------
    targets    manifest -> nv_nights: grism frames per (name, filter, night)
    resolve    Sesame   -> nv_targets: one row per SIMBAD object; aliases
                           that the manifest keeps apart are merged here
    index      BeSS SSA -> nv_bess: every public spectrum record per target
    fetch      download the H-alpha spectra within +-PAD days of a season
    measure    -> nv_halpha: EW, peak, errors for every downloaded spectrum
    tables     -> tables/*.csv, tables/*.md (the deliverable tables)
    figures    -> figures/*.png
    all        everything above, in order

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    cd "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO"
    $PY BeStar_Grism/notes/novelty/bess_novelty_check.py all

SAFETY
------
The shared manifest is opened strictly READ-ONLY (``mode=ro``).  Everything
this script writes lives under ``BeStar_Grism/notes/novelty/``: its own
database ``novelty.sqlite``, a ``cache/`` of raw BeSS/Sesame responses (so a
re-run is offline and reproducible, and the pull date of every response is
recorded), ``tables/`` and ``figures/``.

CONVENTIONS (stated once, used everywhere)
------------------------------------------
* EW sign: POSITIVE = net absorption, NEGATIVE = net emission (the BeSS /
  Be-star literature convention).
* "Season" = observing year 1 August -> 31 July, labelled e.g. ``2025-26``.
* "Night" = the manifest's local-noon-to-noon night label.
* A target is a SIMBAD object, not a manifest string.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import re
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

#: numpy renamed trapz -> trapezoid in 2.0; use whichever this env has.
_trapz = getattr(np, "trapezoid", None) or np.trapz

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MANIFEST = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
#: Everything this script writes goes under WORK.  It defaults to this
#: directory; ``NOVELTY_WORKDIR`` redirects it (the first full run was made on
#: a local SSD and synced back, because several packages were reading the
#: shared spinning disk at once and each small cache write took seconds).
WORK = Path(os.environ.get("NOVELTY_WORKDIR", HERE))
DB = WORK / "novelty.sqlite"
CACHE = WORK / "cache"
TABLES = WORK / "tables"
FIGURES = WORK / "figures"

sys.path.insert(0, str(REPO / "pipeline"))

# --------------------------------------------------------------------------
# Fixed parameters.  Every threshold the verdict depends on is declared here,
# before any BeSS spectrum was looked at, and is printed into the tables.
# --------------------------------------------------------------------------
GRISM_FILTERS = ("hrg", "lrg", "HaGrism", "OGGrism")   # the Step-0 whitelist
MIN_NIGHTS_POOL = 10       # a non-core target enters the survey at >= this
MIN_NIGHTS_SEASON = 10     # a season counts toward BE-N1 only at >= this
PAD_DAYS = 60              # BeSS spectra within +-PAD of a season "bracket" it
HALPHA = 6562.80           # Angstrom (air)
LINE_WIN = (6543.0, 6583.0)        # EW integration window  (+-20 A)
CONT_WINS = ((6520.0, 6540.0), (6586.0, 6606.0))   # linear-continuum windows
SMOOTH_A = 1.0             # peak is read from the profile smoothed to 1 A,
                           # so a single noisy pixel cannot make an "emission"
PEAK_MIN = 1.05            # emission = smoothed peak >= 5% above continuum ...
PEAK_SNR_MIN = 5.0         # ... AND >= 5 sigma of the smoothed continuum rms
BESS_SSA = "http://basebe.obspm.fr/cgi-bin/ssapBE.pl"
BESS_CONE_DEG = 0.02       # 72 arcsec: BeSS coordinates are the star's own
BESS_T0 = "1990-01-01"     # BeSS holds nothing earlier in practice
SESAME = "https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/S?"
POLITE_S = 0.4             # pause between remote requests

#: The ten stars ``ANALYSIS_STRATEGY.md`` §3.2 calls the "core ten", by the
#: manifest's own ``target_key`` in ``stage_bestar_grism``.  Read from the
#: stage table at run time; this tuple only names Vega as the non-member.
STAGE_NON_SCIENCE = ("vega",)

#: Targets the brief names explicitly, whatever their night count.
ALWAYS_INCLUDE = ("tet CrB", "QQ Gem")

#: Manifest display names Sesame cannot parse as written.  Left side is the
#: manifest name (lower-case, grism/exposure suffix already stripped), right
#: side is what is sent to Sesame.  Nothing else is renamed.
SESAME_NAME = {
    "phecda": "gam UMa",
    "ry scuti": "RY Sct",
    "rho oph c": "rho Oph C",
    "alcyone": "eta Tau",
}

#: Manifest names that are not stars / not part of any Be programme and are
#: skipped without a Sesame call (blank name = frames with no target).
SKIP_NAMES = ("",)


# ==========================================================================
# small utilities
# ==========================================================================
def connect() -> sqlite3.Connection:
    """Open this package's own database (never the shared manifest)."""
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    return con


def manifest_ro() -> sqlite3.Connection:
    """Open the shared manifest READ-ONLY (ground rule)."""
    con = sqlite3.connect(f"file:{MANIFEST}?mode=ro", uri=True, timeout=300)
    con.row_factory = sqlite3.Row
    return con


def season_of(night: str) -> str:
    """Observing-year label for a night 'YYYY-MM-DD' (1 Aug -> 31 Jul)."""
    y, m = int(night[:4]), int(night[5:7])
    y0 = y if m >= 8 else y - 1
    return f"{y0}-{str(y0 + 1)[2:]}"


def strip_suffix(name: str) -> str:
    """Remove the grism/exposure text that leaks from filenames into names.

    The archive's filename parser turns ``QQ Gem hrg 1-2e+02s`` into a target
    of its own.  Everything from the grism token onward is not part of the
    star's name.
    """
    return re.sub(r"\s+(hrg|lrg|HaGrism|OGGrism)\b.*$", "", name or "").strip()


def mjd_of(iso: str) -> float:
    """MJD of an ISO date or datetime string (UTC)."""
    t = dt.datetime.fromisoformat(iso[:19] if len(iso) > 10 else iso)
    return (t - dt.datetime(1858, 11, 17)).total_seconds() / 86400.0


def iso_of(mjd: float) -> str:
    return (dt.datetime(1858, 11, 17) + dt.timedelta(days=mjd)).strftime("%Y-%m-%d")


def cached_get(url: str, stem: str, binary: bool = False, tries: int = 4) -> Path:
    """GET ``url`` once; afterwards serve it from ``cache/``.

    The cache is the audit trail: file mtime is the pull date, and the name
    carries a hash of the URL so two different queries can never collide.
    """
    CACHE.mkdir(exist_ok=True)
    h = hashlib.sha1(url.encode()).hexdigest()[:12]
    path = CACHE / f"{stem}_{h}{'.fits' if binary else '.xml'}"
    if path.exists() and path.stat().st_size > 0:
        return path
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MACRO-RLMT novelty check"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            path.write_bytes(data)
            time.sleep(POLITE_S)
            return path
        except Exception as exc:                       # network: retry, then fail loudly
            last = exc
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"GET failed after {tries} tries: {url}: {last}")


# ==========================================================================
# 1. targets — what did WE observe?
# ==========================================================================
def cmd_targets(_args) -> None:
    """Fold the manifest into grism frames per (name, filter, night, era).

    Selection is the strategy's Step-0 rule verbatim: canonical raw-tree
    Light frames, no read error, grism filter whitelist.
    """
    man = manifest_ro()
    q = f"""
        SELECT canonical_target, target_key, filter, night, era_id, COUNT(*) AS n
        FROM frames
        WHERE filter IN ({','.join('?' * len(GRISM_FILTERS))})
          AND tree = 'rawimage' AND is_canonical = 1
          AND imagetyp LIKE 'Light%' AND error IS NULL
        GROUP BY target_key, filter, night, era_id"""
    rows = man.execute(q, GRISM_FILTERS).fetchall()
    core = [r[0] for r in man.execute(
        "SELECT DISTINCT target_key FROM stage_bestar_grism WHERE role='science'")]
    build = man.execute("SELECT * FROM build_meta").fetchall()
    man.close()

    con = connect()
    con.executescript("""
        DROP TABLE IF EXISTS nv_nights;
        CREATE TABLE nv_nights (manifest_name TEXT, target_key TEXT, name TEXT,
            filter TEXT, night TEXT, season TEXT, era_id INTEGER, n INTEGER,
            staged_core INTEGER);
        DROP TABLE IF EXISTS nv_meta;
        CREATE TABLE nv_meta (key TEXT PRIMARY KEY, value TEXT);""")
    for r in rows:
        name = strip_suffix(r["canonical_target"] or "")
        con.execute("INSERT INTO nv_nights VALUES (?,?,?,?,?,?,?,?,?)",
                    (r["canonical_target"], r["target_key"], name, r["filter"],
                     r["night"], season_of(r["night"]), r["era_id"], r["n"],
                     int(r["target_key"] in core and r["target_key"] not in STAGE_NON_SCIENCE)))
    con.execute("INSERT INTO nv_meta VALUES ('targets_built', ?)",
                (dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",))
    con.execute("INSERT INTO nv_meta VALUES ('manifest_build_meta', ?)",
                (repr([tuple(b) for b in build])[:2000],))
    con.commit()
    n_names = con.execute("SELECT COUNT(DISTINCT lower(name)) FROM nv_nights").fetchone()[0]
    n_fr = con.execute("SELECT SUM(n) FROM nv_nights").fetchone()[0]
    print(f"targets: {n_fr} grism light frames, {n_names} distinct stripped names")


# ==========================================================================
# 2. resolve — who are they?
# ==========================================================================
def sesame(name: str) -> dict | None:
    """Resolve one name through CDS Sesame (SIMBAD only).  None if unknown."""
    url = SESAME + urllib.parse.quote(name)
    path = cached_get(url, "sesame_" + re.sub(r"\W+", "_", name))
    root = ET.parse(path).getroot()
    res = root.find(".//Resolver")
    if res is None or res.find("jradeg") is None:
        return None
    g = lambda tag: (res.find(tag).text.strip() if res.find(tag) is not None
                     and res.find(tag).text else None)
    vmag = None
    for m in res.findall("mag"):
        if m.get("band") == "V" and m.find("v") is not None:
            vmag = float(m.find("v").text)
    return {"main_id": g("oname"), "ra": float(g("jradeg")), "dec": float(g("jdedeg")),
            "otype": g("otype"), "sptype": g("spType"), "vmag": vmag}


def cmd_resolve(_args) -> None:
    """One row per SIMBAD object for every surveyed manifest name.

    Surveyed = the staged core ten, the names in ALWAYS_INCLUDE, and every
    other grism target with >= MIN_NIGHTS_POOL nights.  Two manifest names
    that resolve to the same SIMBAD ``main_id`` are ONE target from here on —
    this is where 'PHECDA' and 'gam UMa', 'Zeta Tau' and 'zet Tau', 'Ry Sct'
    and 'Ry Scuti' are reunited.  To catch aliases whose separate night counts
    each fall below the pool threshold, EVERY name with >= 3 nights is
    resolved, and the threshold is applied to the merged object.
    """
    con = connect()
    names = con.execute("""
        SELECT lower(name) AS lname, MIN(name) AS name, COUNT(DISTINCT night) AS nights,
               MAX(staged_core) AS core
        FROM nv_nights GROUP BY lower(name)""").fetchall()
    con.executescript("""
        DROP TABLE IF EXISTS nv_alias;
        CREATE TABLE nv_alias (lname TEXT PRIMARY KEY, name TEXT, main_id TEXT,
                               nights_alone INTEGER, resolved INTEGER);
        DROP TABLE IF EXISTS nv_targets;
        CREATE TABLE nv_targets (main_id TEXT PRIMARY KEY, label TEXT, ra REAL, dec REAL,
            otype TEXT, sptype TEXT, vmag REAL, staged_core INTEGER, named_in_brief INTEGER,
            nights INTEGER, frames INTEGER, manifest_names TEXT);""")
    always = {a.lower() for a in ALWAYS_INCLUDE}
    info = {}
    for r in names:
        if r["lname"] in SKIP_NAMES:
            continue
        if r["nights"] < 3 and not r["core"] and r["lname"] not in always:
            continue
        s = sesame(SESAME_NAME.get(r["lname"], r["name"]))
        con.execute("INSERT INTO nv_alias VALUES (?,?,?,?,?)",
                    (r["lname"], r["name"], s["main_id"] if s else None, r["nights"], int(bool(s))))
        if s:
            info[s["main_id"]] = s
    con.commit()
    for mid, s in info.items():
        al = con.execute("SELECT lname, name FROM nv_alias WHERE main_id=?", (mid,)).fetchall()
        lnames = [a["lname"] for a in al]
        ph = ",".join("?" * len(lnames))
        nights, frames, core = con.execute(
            f"SELECT COUNT(DISTINCT night), SUM(n), MAX(staged_core) FROM nv_nights "
            f"WHERE lower(name) IN ({ph})", lnames).fetchone()
        brief = int(any(l in always for l in lnames))
        if not (core or brief or nights >= MIN_NIGHTS_POOL):
            continue
        con.execute("INSERT INTO nv_targets VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    (mid, re.sub(r"\s+", " ", re.sub(r"^(\*|V\*|NAME)\s+", "", mid)),
                     s["ra"], s["dec"], s["otype"], s["sptype"], s["vmag"], core, brief,
                     nights, frames, "; ".join(sorted(a["name"] for a in al))))
    con.commit()
    unresolved = con.execute(
        "SELECT name, nights_alone FROM nv_alias WHERE resolved=0 ORDER BY nights_alone DESC").fetchall()
    n = con.execute("SELECT COUNT(*) FROM nv_targets").fetchone()[0]
    print(f"resolve: {n} surveyed SIMBAD objects; unresolved names: "
          f"{[(u['name'], u['nights_alone']) for u in unresolved]}")


def target_nights_sql(con, main_id: str):
    """(lnames, placeholder string) for all manifest names of one object."""
    lnames = [a[0] for a in con.execute("SELECT lname FROM nv_alias WHERE main_id=?", (main_id,))]
    return lnames, ",".join("?" * len(lnames))


# ==========================================================================
# 3. index — what does BeSS hold?
# ==========================================================================
def ssa_query(ra: float, dec: float, t0: str, t1: str, stem: str) -> tuple[str, list[dict]]:
    """One BeSS SSA cone query in a time window -> (status, records)."""
    url = (f"{BESS_SSA}?POS={ra:.5f},{dec:.5f}&SIZE={BESS_CONE_DEG}"
           f"&REQUEST=queryData&TIME={t0}/{t1}")
    path = cached_get(url, f"ssa_{stem}_{t0}_{t1}")
    root = ET.parse(path).getroot()
    status = root.find(".//INFO[@name='QUERY_STATUS']").get("value")
    tab = root.find(".//TABLE")
    fields = [f.get("ID") for f in tab.findall("FIELD")] if tab is not None else []
    recs = [dict(zip(fields, [td.text for td in tr])) for tr in root.iter("TR")]
    return status, recs


def ssa_all(ra: float, dec: float, t0: str, t1: str, stem: str) -> list[dict]:
    """All records in [t0, t1], bisecting the window on the 1000-row OVERFLOW.

    The service returns an EMPTY table with status OVERFLOW when a query
    matches more than 1000 spectra, so an unhandled overflow would read as
    'no spectra' — the exact opposite of the truth for the best-observed
    stars.  Every overflow is therefore split in time until it fits.
    """
    status, recs = ssa_query(ra, dec, t0, t1, stem)
    if status == "OK":
        return recs
    if status != "OVERFLOW":
        raise RuntimeError(f"BeSS SSA status {status} for {stem} {t0}/{t1}")
    a, b = mjd_of(t0), mjd_of(t1)
    if b - a < 2:
        raise RuntimeError(f"BeSS overflow inside a 2-day window for {stem}")
    mid = iso_of((a + b) / 2)
    return ssa_all(ra, dec, t0, mid, stem) + ssa_all(ra, dec, mid, t1, stem)


def cmd_index(args) -> None:
    """Every public BeSS spectrum record of every surveyed target."""
    con = connect()
    con.executescript("""
        CREATE TABLE IF NOT EXISTS nv_bess (main_id TEXT, bess_id TEXT, bess_name TEXT,
            mjd REAL, date TEXT, wl_lo REAL, wl_hi REAL, covers_ha INTEGER,
            acref TEXT, curated TEXT, PRIMARY KEY (main_id, bess_id));
        CREATE TABLE IF NOT EXISTS nv_bess_done (main_id TEXT PRIMARY KEY, pulled TEXT, n INTEGER);""")
    today = args.until
    for t in con.execute("SELECT * FROM nv_targets ORDER BY main_id").fetchall():
        if con.execute("SELECT 1 FROM nv_bess_done WHERE main_id=?", (t["main_id"],)).fetchone():
            continue
        stem = re.sub(r"\W+", "_", t["label"])
        recs = ssa_all(t["ra"], t["dec"], BESS_T0, today, stem)
        seen = set()
        for r in recs:
            bid = r["creatorDID"]
            if bid in seen:            # window edges are inclusive on both sides
                continue
            seen.add(bid)
            lo, hi = float(r["boundsSpectralStart"]) * 1e10, float(r["boundsSpectralStop"]) * 1e10
            covers = int(lo <= CONT_WINS[0][0] and hi >= CONT_WINS[1][1])
            con.execute("INSERT OR REPLACE INTO nv_bess VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (t["main_id"], bid, r["targetName"], float(r["locationTime"]),
                         r["boundsTimeStart"][:19], lo, hi, covers, r["acref"], r["curationDate"]))
        con.execute("INSERT INTO nv_bess_done VALUES (?,?,?)",
                    (t["main_id"], dt.datetime.utcnow().isoformat(timespec="seconds") + "Z", len(seen)))
        con.commit()
        print(f"index: {t['label']:<14} {len(seen):5d} BeSS records", flush=True)


# ==========================================================================
# 4. fetch + measure — was H-alpha in emission?
# ==========================================================================
def season_windows(con, main_id: str) -> list[dict]:
    """Our seasons for one object: label, first/last night (MJD), nights, per-grism nights."""
    lnames, ph = target_nights_sql(con, main_id)
    out = []
    for r in con.execute(
            f"SELECT season, MIN(night) a, MAX(night) b, COUNT(DISTINCT night) nights, SUM(n) frames "
            f"FROM nv_nights WHERE lower(name) IN ({ph}) GROUP BY season ORDER BY season", lnames):
        per = con.execute(
            f"SELECT filter, COUNT(DISTINCT night) FROM nv_nights WHERE lower(name) IN ({ph}) "
            f"AND season=? GROUP BY filter", lnames + [r["season"]]).fetchall()
        out.append({"season": r["season"], "first": r["a"], "last": r["b"],
                    "mjd0": mjd_of(r["a"]), "mjd1": mjd_of(r["b"]) + 1.0,
                    "nights": r["nights"], "frames": r["frames"],
                    "by_grism": {p[0]: p[1] for p in per}})
    return out


def cmd_fetch(_args) -> None:
    """Download every H-alpha-covering BeSS spectrum within +-PAD_DAYS of a season."""
    con = connect()
    n = 0
    for t in con.execute("SELECT * FROM nv_targets ORDER BY main_id").fetchall():
        for w in season_windows(con, t["main_id"]):
            for r in con.execute(
                    "SELECT bess_id, acref FROM nv_bess WHERE main_id=? AND covers_ha=1 "
                    "AND mjd BETWEEN ? AND ?", (t["main_id"], w["mjd0"] - PAD_DAYS, w["mjd1"] + PAD_DAYS)):
                cached_get(r["acref"], "bess_" + r["bess_id"].replace("BeSS:", ""), binary=True)
                n += 1
    print(f"fetch: {n} H-alpha spectra in cache")


def read_bess_spectrum(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Wavelength (A), flux and a few header facts from a BeSS FITS file.

    BeSS serves two layouts: a 1-D primary image with a linear CRVAL1/CDELT1
    axis, and (for most recent uploads) an empty primary plus a binary-table
    extension whose first two columns are wavelength and flux.  Both are
    handled; anything else raises rather than being guessed at.
    """
    from astropy.io import fits
    with fits.open(path) as hdul:
        h0 = hdul[0].header
        meta = {"observer": str(h0.get("OBSERVER", ""))[:60], "inst": str(h0.get("BSS_INST", ""))[:60],
                "resolving_power": h0.get("SPE_RPOW") or h0.get("BSS_ITRP") or h0.get("BSS_ESRP"),
                "date_obs": str(h0.get("DATE-OBS", ""))}
        if hdul[0].data is not None and np.ndim(hdul[0].data) >= 1:
            flux = np.asarray(hdul[0].data, float).ravel()
            wl = h0["CRVAL1"] + h0["CDELT1"] * (np.arange(flux.size) + 1 - h0.get("CRPIX1", 1))
        elif len(hdul) > 1 and hasattr(hdul[1], "columns"):
            d = hdul[1].data
            wl = np.asarray(d.field(0), float)
            flux = np.asarray(d.field(1), float)
        else:
            raise ValueError(f"unrecognised BeSS FITS layout: {path.name}")
    return wl, flux, meta


def measure_halpha(wl: np.ndarray, flux: np.ndarray) -> dict | None:
    """The one H-alpha measurement applied to every BeSS spectrum.

    1. Continuum: straight line fitted to the two flanking windows
       ``CONT_WINS`` (5-95% clipped so a telluric line or cosmic ray cannot
       tilt it).  rms of the residual there = the per-pixel noise.
    2. EW = integral of (1 - F/Fc) over ``LINE_WIN``; positive = absorption.
       Error = per-pixel noise summed in quadrature over the window PLUS the
       continuum-placement term (noise of the continuum mean times the window
       width).  The second term dominates for good spectra and is the honest
       one: EW is a ratio to a continuum somebody has to place.
    3. Peak = maximum of F/Fc inside the window AFTER boxcar smoothing to
       ``SMOOTH_A``.  The maximum of N noisy pixels is biased HIGH (a
       max-statistic), so the unsmoothed peak is never used, and the peak is
       only called emission when it clears both ``PEAK_MIN`` and
       ``PEAK_SNR_MIN`` times the rms of the equally-smoothed continuum.

    Returns None when the spectrum does not cover the windows or is unusable.
    """
    ok = np.isfinite(wl) & np.isfinite(flux)
    wl, flux = wl[ok], flux[ok]
    if wl.size < 20 or wl.min() > CONT_WINS[0][0] or wl.max() < CONT_WINS[1][1]:
        return None
    order = np.argsort(wl)
    wl, flux = wl[order], flux[order]
    in_c = ((wl >= CONT_WINS[0][0]) & (wl <= CONT_WINS[0][1])) | \
           ((wl >= CONT_WINS[1][0]) & (wl <= CONT_WINS[1][1]))
    in_l = (wl >= LINE_WIN[0]) & (wl <= LINE_WIN[1])
    if in_c.sum() < 8 or in_l.sum() < 5:
        return None
    wc, fc = wl[in_c], flux[in_c]
    coef = np.polyfit(wc, fc, 1)
    res = fc - np.polyval(coef, wc)
    lo, hi = np.percentile(res, [5, 95])
    keep = (res >= lo) & (res <= hi)
    if keep.sum() >= 6:
        coef = np.polyfit(wc[keep], fc[keep], 1)
    cont = np.polyval(coef, wl)
    if not np.all(cont[in_l | in_c] > 0):
        return None
    norm = flux / cont
    rms = float(np.std(norm[in_c][keep] if keep.sum() >= 6 else norm[in_c], ddof=2))
    dl = float(np.median(np.diff(wl[in_l])))
    ew = float(_trapz(1.0 - norm[in_l], wl[in_l]))
    width = LINE_WIN[1] - LINE_WIN[0]
    ew_err = float(np.hypot(rms * dl * np.sqrt(in_l.sum()), rms / np.sqrt(max(keep.sum(), 1)) * width))
    k = max(1, int(round(SMOOTH_A / dl)))
    sm = np.convolve(norm, np.ones(k) / k, mode="same")
    edge = (wl > wl.min() + SMOOTH_A) & (wl < wl.max() - SMOOTH_A)
    sm_rms = float(np.std(sm[in_c & edge]))
    core = in_l & edge
    i = int(np.argmax(np.where(core, sm, -np.inf)))
    j = int(np.argmin(np.where(core, sm, np.inf)))
    peak = float(sm[i])
    return {"ew": ew, "ew_err": ew_err, "peak": peak, "peak_wl": float(wl[i]),
            "peak_snr": (peak - 1.0) / sm_rms if sm_rms > 0 else np.nan,
            "min_norm": float(sm[j]), "cont_rms": rms, "dl": dl,
            "emission": int(peak >= PEAK_MIN and sm_rms > 0 and (peak - 1.0) / sm_rms >= PEAK_SNR_MIN)}


def cmd_measure(_args) -> None:
    """Measure every cached H-alpha spectrum; failures are recorded, not dropped."""
    con = connect()
    con.executescript("""
        DROP TABLE IF EXISTS nv_halpha;
        CREATE TABLE nv_halpha (main_id TEXT, bess_id TEXT, mjd REAL, date TEXT,
            ew REAL, ew_err REAL, peak REAL, peak_wl REAL, peak_snr REAL, min_norm REAL,
            cont_rms REAL, dl REAL, emission INTEGER, resolving_power REAL,
            observer TEXT, inst TEXT, status TEXT, PRIMARY KEY (main_id, bess_id));""")
    n_ok = n_bad = 0
    for r in con.execute("SELECT * FROM nv_bess WHERE covers_ha=1").fetchall():
        h = hashlib.sha1(r["acref"].encode()).hexdigest()[:12]
        path = CACHE / f"bess_{r['bess_id'].replace('BeSS:', '')}_{h}.fits"
        if not path.exists():
            continue                       # outside every season window: never fetched
        try:
            wl, fl, meta = read_bess_spectrum(path)
            m = measure_halpha(wl, fl)
            status = "ok" if m else "no_coverage_or_unusable"
        except Exception as exc:
            m, meta, status = None, {}, f"read_error: {type(exc).__name__}"
        m = m or {}
        rp = meta.get("resolving_power")
        try:
            rp = float(rp) if rp not in (None, "") else None
        except (TypeError, ValueError):
            rp = None
        con.execute("INSERT INTO nv_halpha VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (r["main_id"], r["bess_id"], r["mjd"], r["date"], m.get("ew"), m.get("ew_err"),
                     m.get("peak"), m.get("peak_wl"), m.get("peak_snr"), m.get("min_norm"),
                     m.get("cont_rms"), m.get("dl"), m.get("emission"), rp,
                     meta.get("observer"), meta.get("inst"), status))
        n_ok += status == "ok"
        n_bad += status != "ok"
    con.commit()
    print(f"measure: {n_ok} measured, {n_bad} unusable (kept in nv_halpha with status)")


# ==========================================================================
# 5. tables — the deliverable
# ==========================================================================
def classify(n_meas: int, n_em: int) -> str:
    """BeSS-side verdict for one target-season, from measured spectra only.

    EMISSION      every measured bracketing spectrum shows emission above the
                  continuum (>= 2 spectra), i.e. two independent confirmations
    EMISSION(1)   the single bracketing spectrum shows it — one witness
    MIXED         some do, some do not (transition, or resolution-dependent
                  weak emission) — look at the profiles
    NO-EMISSION   none of >= 1 measured spectra shows emission above continuum
                  (does NOT exclude weak emission filling the absorption core)
    UNWITNESSED   BeSS holds no usable H-alpha spectrum within +-PAD_DAYS
    """
    if n_meas == 0:
        return "UNWITNESSED"
    if n_em == n_meas:
        return "EMISSION" if n_meas >= 2 else "EMISSION(1)"
    if n_em == 0:
        return "NO-EMISSION"
    return "MIXED"


def build_summary(con) -> list[dict]:
    """One row per (target, season) with our coverage and the BeSS evidence."""
    rows = []
    for t in con.execute("SELECT * FROM nv_targets ORDER BY staged_core DESC, nights DESC").fetchall():
        tot = con.execute("SELECT COUNT(*), SUM(covers_ha), MIN(date), MAX(date) FROM nv_bess "
                          "WHERE main_id=?", (t["main_id"],)).fetchone()
        for w in season_windows(con, t["main_id"]):
            a, b = w["mjd0"], w["mjd1"]
            q = lambda lo, hi: con.execute(
                "SELECT COUNT(*) FROM nv_bess WHERE main_id=? AND covers_ha=1 AND mjd BETWEEN ? AND ?",
                (t["main_id"], lo, hi)).fetchone()[0]
            before = con.execute("SELECT MAX(mjd) FROM nv_bess WHERE main_id=? AND covers_ha=1 AND mjd<?",
                                 (t["main_id"], a)).fetchone()[0]
            after = con.execute("SELECT MIN(mjd) FROM nv_bess WHERE main_id=? AND covers_ha=1 AND mjd>?",
                                (t["main_id"], b)).fetchone()[0]
            m = con.execute(
                "SELECT ew, ew_err, peak, emission FROM nv_halpha WHERE main_id=? AND status='ok' "
                "AND mjd BETWEEN ? AND ? ORDER BY mjd", (t["main_id"], a - PAD_DAYS, b + PAD_DAYS)).fetchall()
            ew = np.array([x["ew"] for x in m]); pk = np.array([x["peak"] for x in m])
            n_em = int(sum(x["emission"] for x in m))
            rows.append({
                "main_id": t["main_id"], "label": t["label"], "manifest_names": t["manifest_names"],
                "core": t["staged_core"],
                "brief": t["named_in_brief"], "sptype": t["sptype"], "otype": t["otype"], "vmag": t["vmag"],
                "bess_total": tot[0], "bess_ha_total": tot[1] or 0,
                "bess_first": (tot[2] or "")[:10], "bess_last": (tot[3] or "")[:10],
                "season": w["season"], "first": w["first"], "last": w["last"],
                "nights": w["nights"], "frames": w["frames"],
                **{f"n_{g}": w["by_grism"].get(g, 0) for g in GRISM_FILTERS},
                "bess_in_season": q(a, b), "bess_bracket": q(a - PAD_DAYS, b + PAD_DAYS),
                "gap_before_d": None if before is None else round(a - before, 1),
                "gap_after_d": None if after is None else round(after - b, 1),
                "n_meas": len(m), "n_emission": n_em,
                "ew_med": float(np.median(ew)) if len(m) else None,
                "ew_min": float(ew.min()) if len(m) else None,
                "ew_max": float(ew.max()) if len(m) else None,
                "ew_err_med": float(np.median([x["ew_err"] for x in m])) if len(m) else None,
                # observer-to-observer scatter (1.4826 x MAD): the EMPIRICAL error of a
                # BeSS EW.  It is ~10x the formal error because every observer places
                # the continuum, removes (or not) tellurics and resolves the core
                # differently.  A change is real only against THIS number.
                "ew_scatter": (float(1.4826 * np.median(np.abs(ew - np.median(ew))))
                               if len(m) >= 3 else None),
                "peak_med": float(np.median(pk)) if len(m) else None,
                "peak_max": float(pk.max()) if len(m) else None,
                "bess_verdict": classify(len(m), n_em)})
    return rows


def fmt(v, nd=2):
    if v is None:
        return "—"
    return f"{v:.{nd}f}" if isinstance(v, float) else str(v)


def cmd_tables(_args) -> None:
    """Write the deliverable tables (CSV for machines, Markdown for the report)."""
    import csv
    TABLES.mkdir(exist_ok=True)
    con = connect()
    rows = build_summary(con)
    with open(TABLES / "target_season_bess.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    def md_table(sel, path, title):
        lines = [f"<!-- emitted by bess_novelty_check.py tables; do not edit -->", f"### {title}", "",
                 "| Star (SIMBAD) [manifest name] | SpT (SIMBAD) | SIMBAD type | Season | Nights (hrg/lrg/Ha/OG) | Our span | "
                 "BeSS Hα all-time | in season | within ±%d d | gap before / after (d) | measured | "
                 "in emission | EW median [min, max] (Å) | EW scatter (Å) | formal err (Å) | "
                 "peak F/Fc median (max) | BeSS verdict |" % PAD_DAYS,
                 "|" + "---|" * 17]
        for r in sel:
            lines.append(
                f"| {r['label']} [{r['manifest_names']}] | {r['sptype'] or '—'} | {r['otype'] or '—'} | {r['season']} | "
                f"{r['nights']} ({r['n_hrg']}/{r['n_lrg']}/{r['n_HaGrism']}/{r['n_OGGrism']}) | "
                f"{r['first']} → {r['last']} | {r['bess_ha_total']} | {r['bess_in_season']} | "
                f"{r['bess_bracket']} | {fmt(r['gap_before_d'], 0)} / {fmt(r['gap_after_d'], 0)} | "
                f"{r['n_meas']} | {r['n_emission']} | "
                + (f"{fmt(r['ew_med'])} [{fmt(r['ew_min'])}, {fmt(r['ew_max'])}]" if r['n_meas'] else "—")
                + f" | {fmt(r['ew_scatter'])} | {fmt(r['ew_err_med'])}"
                + " | " + (f"{fmt(r['peak_med'])} ({fmt(r['peak_max'])})" if r['n_meas'] else "—")
                + f" | {r['bess_verdict']} |")
        Path(path).write_text("\n".join(lines) + "\n")

    md_table([r for r in rows if r["core"]], TABLES / "core_ten.md",
             "Core ten (staged science targets) — our seasons vs BeSS")
    md_table([r for r in rows if r["brief"] and not r["core"]], TABLES / "brief_named.md",
             "Targets named in the brief (θ CrB, QQ Gem)")
    md_table([r for r in rows if not r["core"] and not r["brief"]], TABLES / "pool.md",
             f"Other grism targets with ≥ {MIN_NIGHTS_POOL} nights (the unadopted pool)")

    # ---- the BE-N1 count: an object is BeSS-verified active when at least one
    # season with >= MIN_NIGHTS_SEASON of OUR nights has verdict EMISSION
    # (>= 2 independent bracketing spectra, all in emission).
    def verified(sel):
        return sorted({r["label"] for r in sel
                       if r["nights"] >= MIN_NIGHTS_SEASON and r["bess_verdict"] == "EMISSION"})

    def single(sel):
        return sorted({r["label"] for r in sel if r["nights"] >= MIN_NIGHTS_SEASON
                       and r["bess_verdict"] in ("EMISSION(1)", "MIXED")} - set(verified(sel)))
    core = [r for r in rows if r["core"]]
    rest = [r for r in rows if not r["core"]]
    lines = ["<!-- emitted by bess_novelty_check.py tables; do not edit -->",
             "### BE-N1 count (rule fixed in the script header before any spectrum was measured)", "",
             f"Verified active = SIMBAD-resolved star with ≥ {MIN_NIGHTS_SEASON} of our grism nights in a season "
             f"AND ≥ 2 measured BeSS Hα spectra within ±{PAD_DAYS} d of that season, ALL with smoothed "
             f"peak F/Fc ≥ {PEAK_MIN} at ≥ {PEAK_SNR_MIN:.0f}σ.", "",
             "| Set | Objects | Verified active (BeSS) | n | One-witness or mixed | n |", "|---|---|---|---|---|---|",
             f"| Core ten | {len({r['label'] for r in core})} | {', '.join(verified(core)) or '—'} | "
             f"{len(verified(core))} | {', '.join(single(core)) or '—'} | {len(single(core))} |",
             f"| Brief-named + pool (≥ {MIN_NIGHTS_POOL} nights) | {len({r['label'] for r in rest})} | "
             f"{', '.join(verified(rest)) or '—'} | {len(verified(rest))} | "
             f"{', '.join(single(rest)) or '—'} | {len(single(rest))} |"]
    (TABLES / "be_n1_count.md").write_text("\n".join(lines) + "\n")

    # ---- where would OUR series add something?  Active seasons ranked by how
    # much denser our night coverage is than BeSS's own in-season H-alpha
    # coverage.  A star BeSS already watches weekly is validation material; a
    # star BeSS saw twice while we saw it on 30 nights is where novelty lives.
    act = [r for r in rows if r["nights"] >= MIN_NIGHTS_SEASON
           and r["bess_verdict"] in ("EMISSION", "EMISSION(1)", "MIXED")]
    act.sort(key=lambda r: -r["nights"] / (r["bess_in_season"] + 1))
    lines = ["<!-- emitted by bess_novelty_check.py tables; do not edit -->",
             "### Active seasons ranked by RLMT nights per BeSS in-season Hα spectrum", "",
             "| Star [manifest name] | SpT | Season | RLMT nights (hrg/lrg) | BeSS Hα in season | "
             "ratio nights/(BeSS+1) | BeSS EW median (Å) | BeSS verdict | core ten? |", "|" + "---|" * 9]
    lines += [f"| {r['label']} [{r['manifest_names']}] | {r['sptype'] or '—'} | {r['season']} | "
              f"{r['nights']} ({r['n_hrg']}/{r['n_lrg']}) | {r['bess_in_season']} | "
              f"{r['nights'] / (r['bess_in_season'] + 1):.1f} | {fmt(r['ew_med'])} | {r['bess_verdict']} | "
              f"{'yes' if r['core'] else 'no'} |" for r in act]
    (TABLES / "novelty_ranking.md").write_text("\n".join(lines) + "\n")

    # ---- did the BeSS EW move across our campaign?  A straight line through
    # every measured spectrum of each core-ten / brief-named star, with the
    # slope error taken from the RESIDUAL SCATTER (never from the formal
    # errors, which are ~10x too small — see ew_scatter above; and no
    # max(chi2,1) rescaling: the scatter is the error model).
    lines = ["<!-- emitted by bess_novelty_check.py tables; do not edit -->",
             "### BeSS Hα EW trend across the campaign (core ten + brief-named; ≥ 5 spectra)", "",
             "| Star [manifest name] | N spectra | Observers | Span | EW mean (Å) | EW sd (Å) | "
             "slope (Å / 100 d) | ± (from residual scatter) | slope / σ | residual sd (Å) |",
             "|" + "---|" * 10]
    for t in con.execute("SELECT * FROM nv_targets WHERE staged_core=1 OR named_in_brief=1 "
                         "ORDER BY staged_core DESC, nights DESC").fetchall():
        m = con.execute("SELECT mjd, ew, observer, date FROM nv_halpha WHERE main_id=? AND status='ok' "
                        "ORDER BY mjd", (t["main_id"],)).fetchall()
        if len(m) < 5:
            continue
        x = np.array([r["mjd"] for r in m]); y = np.array([r["ew"] for r in m])
        coef = np.polyfit(x, y, 1)
        res = y - np.polyval(coef, x)
        sd_res = float(res.std(ddof=2))
        se = sd_res / float(np.sqrt(((x - x.mean()) ** 2).sum()))
        lines.append(f"| {t['label']} [{t['manifest_names']}] | {len(m)} | "
                     f"{len({r['observer'] for r in m})} | {m[0]['date'][:10]} → {m[-1]['date'][:10]} | "
                     f"{y.mean():.2f} | {y.std(ddof=1):.2f} | {coef[0] * 100:+.2f} | {se * 100:.2f} | "
                     f"{coef[0] / se:+.1f} | {sd_res:.2f} |")
    (TABLES / "ew_trends.md").write_text("\n".join(lines) + "\n")

    # ---- alias reunions and unresolved names: evidence for the ledger notes
    al = con.execute("""SELECT t.label, t.manifest_names, t.nights, t.frames FROM nv_targets t
                        WHERE t.manifest_names LIKE '%;%' ORDER BY t.nights DESC""").fetchall()
    un = con.execute("SELECT name, nights_alone FROM nv_alias WHERE resolved=0 "
                     "ORDER BY nights_alone DESC").fetchall()
    lines = ["<!-- emitted by bess_novelty_check.py tables; do not edit -->",
             "### Manifest names reunited by SIMBAD identity", "",
             "| SIMBAD object | Manifest names merged | Nights (merged) | Frames |", "|---|---|---|---|"]
    lines += [f"| {a['label']} | {a['manifest_names']} | {a['nights']} | {a['frames']} |" for a in al]
    lines += ["", "### Names Sesame could not resolve (≥ 3 nights)", "", "| Manifest name | Nights |", "|---|---|"]
    lines += [f"| {u['name']} | {u['nights_alone']} |" for u in un]
    (TABLES / "aliases.md").write_text("\n".join(lines) + "\n")

    st = con.execute("SELECT status, COUNT(*) FROM nv_halpha GROUP BY status").fetchall()
    print("tables: wrote", ", ".join(p.name for p in sorted(TABLES.iterdir())))
    print("        nv_halpha status:", [tuple(s) for s in st])
    print((TABLES / "be_n1_count.md").read_text())


# ==========================================================================
# 6. figures — the plots do the talking
# ==========================================================================
def cmd_figures(_args) -> None:
    """Two figure families.

    ``coverage_<set>.png``  one row per star: our grism nights (ticks) against
        every BeSS H-alpha spectrum (filled = emission above continuum, open
        = none), 2024-08 -> today.  Shows at a glance who is bracketed.
    ``profiles_<star>.png`` the bracketing BeSS H-alpha profiles themselves,
        continuum-normalised, with the EW series beneath and our nights
        marked — the evidence behind each verdict.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from macro_core import plotstyle as ps
    FIGURES.mkdir(exist_ok=True)
    con = connect()
    targets = con.execute("SELECT * FROM nv_targets ORDER BY staged_core DESC, nights DESC").fetchall()
    x0, x1 = mjd_of("2024-08-01"), mjd_of(dt.date.today().isoformat())

    def nights_mjd(mid):
        lnames, ph = target_nights_sql(con, mid)
        return np.array([mjd_of(r[0]) + 0.5 for r in con.execute(
            f"SELECT DISTINCT night FROM nv_nights WHERE lower(name) IN ({ph})", lnames)])

    def coverage(sel, fname, title):
        with ps.context():
            fig, ax = plt.subplots(figsize=(ps.COL_DOUBLE, 0.8 + 0.27 * len(sel)))
            for i, t in enumerate(sel):
                y = len(sel) - 1 - i
                ax.plot(nights_mjd(t["main_id"]), np.full_like(nights_mjd(t["main_id"]), y + 0.18),
                        "|", color=ps.ACCENT, ms=6, mew=0.9)
                allb = np.array([r[0] for r in con.execute(
                    "SELECT mjd FROM nv_bess WHERE main_id=? AND covers_ha=1", (t["main_id"],))])
                ax.plot(allb, np.full_like(allb, y - 0.18), "|", color=ps.WISP, ms=5, mew=0.8)
                for em, kw in ((1, ps.measurement_kw(ps.BAD, "o", size=3.2)),
                               (0, ps.floor_kw(ps.MUTED, "o", size=3.2))):
                    mm = np.array([r[0] for r in con.execute(
                        "SELECT mjd FROM nv_halpha WHERE main_id=? AND status='ok' AND emission=?",
                        (t["main_id"], em))])
                    ax.plot(mm, np.full_like(mm, y - 0.18), **kw)
            ax.set_yticks(range(len(sel)))
            ax.set_yticklabels([t["label"] for t in sel][::-1])
            ax.set_xlim(x0, x1); ax.set_ylim(-0.7, len(sel) - 0.3)
            ticks = [mjd_of(f"{y}-{m:02d}-01") for y in (2024, 2025, 2026) for m in (1, 4, 7, 10)
                     if x0 <= mjd_of(f"{y}-{m:02d}-01") <= x1]
            ax.set_xticks(ticks); ax.set_xticklabels([iso_of(v)[:7] for v in ticks])
            ax.set_xlabel("date (UTC)")
            ax.set_title(title, loc="left")
            from matplotlib.lines import Line2D
            ax.legend(handles=[
                Line2D([], [], ls="none", marker="|", color=ps.ACCENT, label="RLMT grism night"),
                Line2D([], [], ls="none", marker="|", color=ps.WISP, label="BeSS Hα spectrum (not fetched)"),
                ps.measurement_handle("BeSS Hα: emission above continuum", ps.BAD, "o"),
                Line2D([], [], label="BeSS Hα: no emission above continuum",
                       **ps.floor_kw(ps.MUTED, "o", size=4.0))],
                loc="upper center", bbox_to_anchor=(0.5, -0.9 / (0.8 + 0.27 * len(sel)) - 0.02), ncol=2, frameon=False)
            fig.savefig(FIGURES / fname, dpi=ps.WEB_DPI, bbox_inches="tight")
            plt.close(fig)

    core = [t for t in targets if t["staged_core"]]
    rest = [t for t in targets if not t["staged_core"]]
    coverage(core, "coverage_core_ten.png", "Core ten: RLMT grism nights vs BeSS Hα spectra")
    half = (len(rest) + 1) // 2
    coverage(rest[:half], "coverage_pool_a.png", "Brief-named and pool targets (1/2)")
    coverage(rest[half:], "coverage_pool_b.png", "Brief-named and pool targets (2/2)")

    # ---- per-star profile sheets, only where there is something to show
    n_sheets = 0
    for t in targets:
        m = con.execute("SELECT h.*, b.acref FROM nv_halpha h JOIN nv_bess b USING (main_id, bess_id) "
                        "WHERE h.main_id=? AND h.status='ok' ORDER BY h.mjd", (t["main_id"],)).fetchall()
        if not m:
            continue
        with ps.context():
            fig, (a1, a2) = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 3.0),
                                         gridspec_kw={"width_ratios": [1, 1.25]})
            cols = ps.ordinal_colors(len(m))
            for c, r in zip(cols, m):
                h = hashlib.sha1(r["acref"].encode()).hexdigest()[:12]
                wl, fl, _ = read_bess_spectrum(CACHE / f"bess_{r['bess_id'].replace('BeSS:', '')}_{h}.fits")
                o = np.argsort(wl); wl, fl = wl[o], fl[o]
                sel = (wl > 6515) & (wl < 6611) & np.isfinite(fl)
                inc = sel & (((wl >= CONT_WINS[0][0]) & (wl <= CONT_WINS[0][1])) |
                             ((wl >= CONT_WINS[1][0]) & (wl <= CONT_WINS[1][1])))
                cont = np.polyval(np.polyfit(wl[inc], fl[inc], 1), wl[sel])
                a1.plot(wl[sel], fl[sel] / cont, color=c, lw=0.7)
            a1.axhline(1, **ps.reference_kw())
            for w in LINE_WIN:
                a1.axvline(w, color=ps.WISP, lw=0.6)
            a1.set_xlabel("wavelength (Å)"); a1.set_ylabel("F / F$_c$")
            a1.set_title(f"{t['label']}  ({t['sptype'] or '?'}) — {len(m)} BeSS spectra", loc="left")
            mj = np.array([r["mjd"] for r in m]); ew = np.array([r["ew"] for r in m])
            er = np.array([r["ew_err"] for r in m]); em = np.array([r["emission"] for r in m], bool)
            for flag, kw in ((True, ps.measurement_kw(ps.BAD, "o", size=4.0)),
                             (False, ps.floor_kw(ps.MUTED, "o", size=4.0))):
                if (em == flag).any():
                    a2.errorbar(mj[em == flag], ew[em == flag], er[em == flag],
                                ecolor=ps.MUTED, elinewidth=0.6, **kw)
            nm = nights_mjd(t["main_id"])
            lo, hi = a2.get_ylim()
            a2.plot(nm, np.full_like(nm, hi), "|", color=ps.ACCENT, ms=7, mew=0.9, clip_on=False)
            a2.axhline(0, **ps.reference_kw())
            a2.invert_yaxis()
            a2.set_xlabel("MJD"); a2.set_ylabel("Hα EW (Å; negative = emission)")
            a2.set_title("blue ticks: RLMT nights; filled: emission", loc="left")
            fig.tight_layout()
            fig.savefig(FIGURES / f"profiles_{re.sub(r'[^A-Za-z0-9]+', '_', t['label'])}.png", dpi=ps.WEB_DPI)
            plt.close(fig)
            n_sheets += 1
    print(f"figures: 3 coverage panels + {n_sheets} profile sheets in {FIGURES}")


def cmd_all(args) -> None:
    for f in (cmd_targets, cmd_resolve, cmd_index, cmd_fetch, cmd_measure, cmd_tables, cmd_figures):
        f(args)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["targets", "resolve", "index", "fetch", "measure",
                                    "tables", "figures", "all"])
    ap.add_argument("--until", default="2026-10-03",
                    help="last date of the BeSS index (fixed so the cache keys are stable)")
    args = ap.parse_args(argv)
    globals()[f"cmd_{args.cmd}"](args)


if __name__ == "__main__":
    main()
