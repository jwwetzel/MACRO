#!/usr/bin/env python
"""L1 external lookups for the legacy census: VSX and SIMBAD, cached.

WHY THE CENSUS NEEDS ANYTHING EXTERNAL
--------------------------------------
Two census questions cannot be answered from FITS headers:

* **Which legacy targets are eclipsing/contact binaries, and with what
  period?**  The pre-registered gate G1 counts systems with
  "minimum-bearing nights" in ≥ 3 seasons, and a minimum-bearing night is
  defined by the catalogue period P (a run ≥ P/2) or ephemeris.  The
  authority used is the AAVSO International Variable Star Index (VSX), read
  from its VizieR mirror ``B/vsx/vsx``.
* **Where on the sky are the RLMT-era project targets?**  The overlap query
  matches by position as well as by name, so it needs one catalogue
  position per project target.  The names come from the RLMT manifest's
  ``project_counts`` table; the positions from SIMBAD (through the CDS
  Sesame name resolver).

WHAT IT WRITES (into ``products/legacy/legacy.sqlite``)
-------------------------------------------------------
* ``vsx_matches``   one row per queried legacy target: the VSX object it was
                    matched to (or NULLs), type, period, epoch, how it was
                    matched (``name`` / ``position``) and the separation.
* ``rlmt_targets``  one row per RLMT-era project target: project, name,
                    alias key, ICRS position and where the position came from.
* ``external_meta`` query dates, row counts, radii, catalogue names.

The raw VizieR rows are cached in ``products/legacy/external/vsx_cone.csv``
so the match can be re-derived (and audited) without the network; SIMBAD
answers are cached in ``simbad_targets.csv``.  A cached file is NEVER
silently re-fetched: ``--force`` re-pulls.

MATCHING RULE (pure function ``choose_vsx_match``; unit-tested)
---------------------------------------------------------------
VSX objects within ``CONE_ARCSEC`` of the target's requested coordinates
are candidates.  (1) A candidate whose normalized NAME equals the target's
alias key wins (``name``).  (2) Otherwise the nearest candidate within
``POSITION_ARCSEC`` wins (``position``) — the scheduler's OBJRA/OBJDEC are
catalogue coordinates, so a real match sits within a few arcseconds.
(3) Otherwise there is no match.  Fields and galaxies therefore stay
unmatched unless a catalogued variable sits on the requested position, and
the match method and separation are stored so a reader can discount
position-only matches.

USAGE
-----
    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_legacy_external.py            # uses caches
    … build_legacy_external.py --force                       # re-query
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_legacy import census as lc                        # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_DB = REPO_ROOT / "products" / "legacy" / "legacy.sqlite"
DEFAULT_RLMT = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
CACHE_DIR = REPO_ROOT / "products" / "legacy" / "external"

VSX_CATALOG = "B/vsx/vsx"
#: Candidates are fetched within this radius of the requested coordinates.
CONE_ARCSEC = 60.0
#: A position-only match must be this close (see module docstring).
POSITION_ARCSEC = 10.0
#: Targets are queried when they could matter to gate G1 or to a timings
#: table: observed on at least this many nights, or with a run this long.
MIN_NIGHTS_TO_QUERY = 2
MIN_RUN_HOURS_TO_QUERY = 1.0
#: VizieR upload batch size (one HTTP request per batch).
BATCH = 150

#: SIMBAD identifiers for project targets whose manifest name SIMBAD does
#: not resolve as written.
SIMBAD_NAME = {"2023ixf": "SN 2023ixf", "lameri": "lam Eri",
               "phecda": "gam UMa", "phileo": "phi Leo", "tetcrb": "tet CrB"}
#: The SN field is also a target under its host's name.
EXTRA_TARGETS = (("SN2023ixf_LightCurve", "M101", "m101", "M 101"),)


# ---------------------------------------------------------------------------
# Pure matching logic
# ---------------------------------------------------------------------------
def choose_vsx_match(target_key: str, candidates: list[dict],
                     position_arcsec: float = POSITION_ARCSEC):
    """Pick the VSX object for one target.  See the module docstring.

    ``candidates`` are dicts with at least ``Name`` and ``sep_arcsec``.
    Returns ``(candidate, method)`` or ``(None, None)``.
    """
    by_name = [c for c in candidates
               if lc.target_key(c["Name"])[0] == target_key]
    if by_name:
        return min(by_name, key=lambda c: c["sep_arcsec"]), "name"
    near = [c for c in candidates if c["sep_arcsec"] <= position_arcsec]
    if near:
        return min(near, key=lambda c: c["sep_arcsec"]), "position"
    return None, None


# ---------------------------------------------------------------------------
# VSX
# ---------------------------------------------------------------------------
def fetch_vsx(targets: pd.DataFrame) -> pd.DataFrame:
    """Cone-search VSX around every target; return all candidate rows."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    viz = Vizier(columns=["OID", "Name", "Type", "Period", "Epoch", "max",
                          "min", "RAJ2000", "DEJ2000"], row_limit=-1)
    out = []
    for i in range(0, len(targets), BATCH):
        chunk = targets.iloc[i:i + BATCH].reset_index(drop=True)
        coords = SkyCoord(chunk["ra_deg"].to_numpy() * u.deg,
                          chunk["dec_deg"].to_numpy() * u.deg)
        for attempt in range(3):
            try:
                res = viz.query_region(coords, radius=CONE_ARCSEC * u.arcsec,
                                       catalog=VSX_CATALOG)
                break
            except Exception as exc:       # noqa: BLE001 — retried, then raised
                if attempt == 2:
                    raise
                print(f"  VizieR retry after {type(exc).__name__}", flush=True)
                time.sleep(10)
        if len(res):
            t = res[0].to_pandas()
            # _q is the 1-based index of the uploaded position.
            t["target_key"] = chunk["target_key"].to_numpy()[t["_q"] - 1]
            t["q_ra"] = chunk["ra_deg"].to_numpy()[t["_q"] - 1]
            t["q_dec"] = chunk["dec_deg"].to_numpy()[t["_q"] - 1]
            out.append(t)
        print(f"  VSX {min(i + BATCH, len(targets)):,}/{len(targets):,}",
              flush=True)
    if not out:
        return pd.DataFrame()
    cand = pd.concat(out, ignore_index=True)
    d_ra = (cand["RAJ2000"] - cand["q_ra"]) * np.cos(np.radians(cand["q_dec"]))
    d_dec = cand["DEJ2000"] - cand["q_dec"]
    cand["sep_arcsec"] = np.hypot(d_ra, d_dec) * 3600.0
    return cand.drop(columns=[c for c in ("_q",) if c in cand])


def match_vsx(targets: pd.DataFrame, cand: pd.DataFrame) -> pd.DataFrame:
    """One row per queried target with its chosen VSX match (or NULLs)."""
    groups = {k: g.to_dict("records") for k, g in cand.groupby("target_key")} \
        if len(cand) else {}
    rows = []
    for t in targets.itertuples():
        best, method = choose_vsx_match(t.target_key,
                                        groups.get(t.target_key, []))
        rows.append(dict(
            target_key=t.target_key, q_ra_deg=t.ra_deg, q_dec_deg=t.dec_deg,
            n_candidates=len(groups.get(t.target_key, [])),
            vsx_name=best["Name"] if best else None,
            vsx_oid=int(best["OID"]) if best else None,
            vsx_type=best["Type"] if best else None,
            period=_num(best, "Period"), epoch=_num(best, "Epoch"),
            vsx_ra_deg=best["RAJ2000"] if best else None,
            vsx_dec_deg=best["DEJ2000"] if best else None,
            sep_arcsec=round(best["sep_arcsec"], 2) if best else None,
            match_method=method))
    return pd.DataFrame(rows)


def _num(rec, key):
    """A positive float from a candidate record, else None (masked → NaN)."""
    if not rec:
        return None
    v = rec.get(key)
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) and v > 0 else None


# ---------------------------------------------------------------------------
# RLMT-era project targets
# ---------------------------------------------------------------------------
def project_targets(rlmt_path: Path) -> pd.DataFrame:
    """Distinct (project, target, key) from the RLMT manifest, read-only.

    The Dwarf survey's ``__dw_survey__`` sentinel is expanded into its
    field keys (``dw1…``), each with the manifest's median header position
    (those fields have no SIMBAD name to resolve).
    """
    con = sqlite3.connect(f"file:{rlmt_path}?mode=ro", uri=True)
    try:
        pc = pd.read_sql_query(
            "SELECT DISTINCT project, target, target_key FROM project_counts",
            con)
        dw = pd.read_sql_query("""
            SELECT target_key, canonical_target AS target, ra_deg, dec_deg
            FROM frames WHERE is_canonical = 1 AND target_key LIKE 'dw1%'
              AND ra_deg IS NOT NULL""", con)
    finally:
        con.close()
    rows = []
    for r in pc.itertuples():
        if r.target_key == "__dw_survey__":
            for key, g in dw.groupby("target_key"):
                rows.append(dict(project=r.project, target=g["target"].iloc[0],
                                 target_key=key,
                                 ra_deg=float(g["ra_deg"].median()),
                                 dec_deg=float(g["dec_deg"].median()),
                                 position_source="RLMT manifest median"))
        else:
            rows.append(dict(project=r.project, target=r.target,
                             target_key=r.target_key, ra_deg=None,
                             dec_deg=None, position_source=None))
    for project, target, key, _simbad in EXTRA_TARGETS:
        rows.append(dict(project=project, target=target, target_key=key,
                         ra_deg=None, dec_deg=None, position_source=None))
    return pd.DataFrame(rows).drop_duplicates(["project", "target_key"])


def resolve_simbad(tg: pd.DataFrame) -> pd.DataFrame:
    """Fill ICRS positions for named project targets from SIMBAD.

    Goes through the CDS Sesame name resolver (``SkyCoord.from_name``),
    which answers from SIMBAD and is mirrored; the SIMBAD TAP endpoint that
    ``astroquery.simbad`` needs was unreachable from this network on the
    day this was written, and a census must not hang on one host.  A name
    that does not resolve is recorded as such — the target then matches by
    name only, never by a guessed position.
    """
    from astropy.coordinates import SkyCoord
    from astropy.coordinates.name_resolve import NameResolveError
    from astropy.utils.data import conf as data_conf
    extra = {key: simbad for _p, _t, key, simbad in EXTRA_TARGETS}
    tg = tg.copy()
    for i, r in tg.iterrows():
        if pd.notna(r["ra_deg"]):
            continue
        name = extra.get(r["target_key"]) or SIMBAD_NAME.get(
            r["target_key"], r["target"])
        try:
            with data_conf.set_temp("remote_timeout", 30):
                c = SkyCoord.from_name(name)
        except NameResolveError as exc:
            tg.at[i, "position_source"] = f"not resolved: '{name}'"
            print(f"  not resolved: {name} ({str(exc)[:60]})")
            continue
        tg.at[i, "ra_deg"] = float(c.ra.deg)
        tg.at[i, "dec_deg"] = float(c.dec.deg)
        tg.at[i, "position_source"] = f"SIMBAD via CDS Sesame: '{name}'"
        time.sleep(0.2)
    return tg


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--rlmt-manifest", type=Path, default=DEFAULT_RLMT)
    ap.add_argument("--force", action="store_true",
                    help="re-query even when a cache file exists")
    args = ap.parse_args(argv)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    meta = {}

    con = sqlite3.connect(args.db, timeout=120)
    # Write-ahead logging: appends sequentially instead of journalling every
    # replaced page — on the spinning disk this database lives on, replacing
    # the derived tables is otherwise tens of minutes.
    con.execute("PRAGMA journal_mode=WAL")
    try:
        # ---- RLMT-era targets ---------------------------------------------
        simbad_cache = CACHE_DIR / "simbad_targets.csv"
        if simbad_cache.exists() and not args.force:
            tg = pd.read_csv(simbad_cache)
            print(f"rlmt_targets: {len(tg)} from cache")
        else:
            tg = resolve_simbad(project_targets(args.rlmt_manifest))
            tg.to_csv(simbad_cache, index=False)
            meta["simbad_query_utc"] = now
            print(f"rlmt_targets: {len(tg)} resolved")
        tg.to_sql("rlmt_targets", con, index=False, if_exists="replace")

        # ---- VSX ----------------------------------------------------------
        targets = pd.read_sql_query(f"""
            SELECT target_key, ra_deg, dec_deg FROM targets
            WHERE ra_deg IS NOT NULL AND dec_deg IS NOT NULL
              AND (n_nights >= {MIN_NIGHTS_TO_QUERY}
                   OR longest_run_h >= {MIN_RUN_HOURS_TO_QUERY})
            ORDER BY target_key""", con)
        vsx_cache = CACHE_DIR / "vsx_cone.csv"
        if vsx_cache.exists() and not args.force:
            cand = pd.read_csv(vsx_cache)
            missing = set(targets["target_key"]) - set(
                pd.read_csv(CACHE_DIR / "vsx_queried.csv")["target_key"])
            if missing:
                print(f"VSX cache lacks {len(missing)} targets — querying them")
                more = fetch_vsx(targets[targets["target_key"].isin(missing)])
                cand = pd.concat([cand, more], ignore_index=True)
                cand.to_csv(vsx_cache, index=False)
                meta["vsx_query_utc"] = now
        else:
            print(f"querying VSX for {len(targets):,} targets …")
            cand = fetch_vsx(targets)
            cand.to_csv(vsx_cache, index=False)
            meta["vsx_query_utc"] = now
        targets[["target_key"]].to_csv(CACHE_DIR / "vsx_queried.csv",
                                       index=False)
        matches = match_vsx(targets, cand)
        matches.to_sql("vsx_matches", con, index=False, if_exists="replace")
        n_match = int(matches["vsx_name"].notna().sum())
        print(f"vsx_matches: {len(matches):,} targets queried, "
              f"{n_match:,} matched "
              f"({int((matches['match_method'] == 'name').sum()):,} by name)")

        meta.update(vsx_catalog=VSX_CATALOG, vsx_cone_arcsec=CONE_ARCSEC,
                    vsx_position_arcsec=POSITION_ARCSEC,
                    vsx_targets_queried=len(matches),
                    vsx_targets_matched=n_match,
                    vsx_min_nights=MIN_NIGHTS_TO_QUERY,
                    vsx_min_run_hours=MIN_RUN_HOURS_TO_QUERY,
                    external_built_utc=now)
        con.execute("CREATE TABLE IF NOT EXISTS external_meta "
                    "(key TEXT PRIMARY KEY, value TEXT)")
        con.executemany("INSERT OR REPLACE INTO external_meta VALUES (?,?)",
                        [(k, str(v)) for k, v in meta.items()])
        con.commit()
    finally:
        con.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
