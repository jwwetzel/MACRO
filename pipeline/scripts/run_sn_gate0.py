#!/usr/bin/env python
"""Run SN 2023ixf GATE 0 — the three blocking activities, from the pixels.

WHAT THIS SCRIPT DOES
---------------------
``SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md`` §4 Step 0 forbids every
downstream step until three things land.  This script produces all three as
tables in the manifest, so that each is a query rather than a paragraph:

    0a  freeze    the globally deduplicated, alias-merged frame list
    0b  census    the supernova's own peak ADU, measured in every image,
                  screened against S2's MEASURED per-mode ceiling
    0c  triage    the slitless series judged against the pre-registered
                  promotion criterion — using S2c's per-frame dispersion
                  measurement, never a filter label

and then answers, in ``sn_g0_verdict``, the three questions Gate 0 exists to
answer, each with the number that decides it.

THE DEFECT THIS STAGE IS BUILT AROUND
-------------------------------------
Filter slot ``6`` was assumed to be a grism because of its name.  It is
mixed.  On this target S2c measures 61 of its 83 frames as spectra, 3 as
ordinary direct images, and 19 as uncertifiable.  Every rule here reads that
measurement (``frame_dispersion``) and no rule reads a filter string.

SUBCOMMANDS
-----------
    freeze      0a: build sn_g0_frames / sn_g0_dedup / sn_g0_pointing
    census      0b step 1: enqueue every frame for pixel measurement
    measure     0b step 2: measure pending frames.  SAFE TO RE-RUN — a
                killed run loses only the frames in flight, which stay
                pending.  ``--max-seconds N`` returns cleanly after N
                seconds so this can be driven under a command timeout.
    matrix      0b step 3: fold the census into the filter x night
                saturation matrix and the per-band summary
    triage      0c: the grism series assessment
    verdicts    the three Gate 0 answers + the venue posture
    report      render docs/SN2023ixf_LightCurve/sn_gate0.html
    status      progress (read-only)

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/run_sn_gate0.py freeze
    $PY pipeline/scripts/run_sn_gate0.py census
    $PY pipeline/scripts/run_sn_gate0.py measure --workers 6
    $PY pipeline/scripts/run_sn_gate0.py matrix
    $PY pipeline/scripts/run_sn_gate0.py triage
    $PY pipeline/scripts/run_sn_gate0.py verdicts
    $PY pipeline/scripts/run_sn_gate0.py report

CONCURRENCY AND SAFETY
----------------------
The archive is opened strictly READ-ONLY; no FITS file is ever written.  The
manifest is a WAL database other stages write concurrently, so every
connection sets a five-minute busy timeout and every write transaction is
kept short (results flush in small batches).  Only ``sn_g0_*`` tables are
ever written; no existing table is modified.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import math
import sqlite3
import sys
import time
import warnings
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

# Make the pipeline package importable regardless of the working directory.
PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_sn import SN_G0_CODE_VERSION            # noqa: E402
from macro_sn import gate0 as g0                   # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
DEFAULT_ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: Five minutes.  A backstop against a long WAL checkpoint by another
#: writer, not an expected wait.
BUSY_TIMEOUT_MS = 300_000

#: Default worker count.  Capped low on purpose: the archive lives on a
#: spinning drive that a legacy-archive repair transfer may be using at the
#: same time, and six concurrent decompressions already saturate it.
DEFAULT_WORKERS = 6

#: Fallback plate scale, degrees per pixel, used only when a frame carries
#: no CDELT card at all.  PROVENANCE: the CD matrix of this campaign's
#: plate-solved frames (1.5002e-4 deg/px = 0.540"/px at 1x1 binning).
FALLBACK_SCALE_DEG_PX = 1.5e-4


def utcnow() -> str:
    """UTC timestamp, seconds resolution, the format every *_build_meta
    table in this repository already uses."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def connect(path: Path, read_only: bool = False) -> sqlite3.Connection:
    """Open the manifest with the house settings."""
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=BUSY_TIMEOUT_MS / 1000.0)
    con.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MS}")
    con.row_factory = sqlite3.Row
    return con


def record_meta(con: sqlite3.Connection, **kw) -> None:
    """Write the build-meta rows for this stage (key/value, replace)."""
    con.execute("""CREATE TABLE IF NOT EXISTS sn_g0_build_meta (
                       key TEXT PRIMARY KEY, value TEXT)""")
    for k, v in kw.items():
        con.execute("INSERT OR REPLACE INTO sn_g0_build_meta(key, value) "
                    "VALUES (?, ?)", (k, str(v)))


# ===========================================================================
# 0a — FREEZE
# ===========================================================================
FREEZE_DDL = """
DROP TABLE IF EXISTS sn_g0_frames;
CREATE TABLE sn_g0_frames (
    obs_rowid   INTEGER PRIMARY KEY,
    path        TEXT NOT NULL,
    tree        TEXT,
    basename    TEXT,
    night       TEXT,
    jd          REAL,
    phase_d     REAL,
    filter      TEXT,
    exptime     REAL,
    readoutm    TEXT,
    xbinning    REAL,
    naxis1      REAL,
    naxis2      REAL,
    canonical_target TEXT,
    raw_object  TEXT,
    observer    TEXT,
    epoch_role  TEXT,
    band_role   TEXT,
    dispersion_verdict TEXT,
    dispersion_status  TEXT,
    dispersion_class   TEXT,
    is_image    INTEGER,
    pltsolvd    REAL,
    crval1      REAL,
    crval2      REAL,
    objctra     TEXT,
    objctdec    TEXT,
    fwhm        REAL,
    airmass     REAL,
    zmag        REAL,
    dup_copies  INTEGER,
    era_id      INTEGER
);
DROP TABLE IF EXISTS sn_g0_dedup;
CREATE TABLE sn_g0_dedup (
    scope       TEXT,
    n_rows      INTEGER,
    n_canonical INTEGER,
    note        TEXT
);
DROP TABLE IF EXISTS sn_g0_pointing;
CREATE TABLE sn_g0_pointing (
    scope       TEXT PRIMARY KEY,
    n_frames    INTEGER,
    d_ra_arcsec REAL,
    d_dec_arcsec REAL,
    resid_med_px REAL,
    resid_p90_px REAL,
    resid_max_px REAL
);
"""


def cmd_freeze(args) -> int:
    """0a — the manifest freeze.

    Three products, and the second and third exist because the strategy's
    own headline numbers turned out to be claims nobody could re-run:

    * ``sn_g0_frames``   one row per unique frame on this sky, at any epoch,
      alias-merged and globally deduplicated, with the S2c dispersion
      verdict and every derived role attached;
    * ``sn_g0_dedup``    the accounting behind it — catalog rows in, unique
      frames out, per tree;
    * ``sn_g0_pointing`` the MEASURED commanded-versus-true pointing offset,
      which is what makes a search box for the unsolved frames a defensible
      size rather than a guess.
    """
    con = connect(args.manifest)
    try:
        con.executescript(FREEZE_DDL)

        # --- the frames themselves -------------------------------------
        raw = [dict(r) for r in con.execute(g0.FREEZE_SQL)]
        rows = g0.freeze_rows(raw)

        # How many catalog copies does each unique frame stand for?  The
        # dup_group is S0's global key; counting its members is the honest
        # measure of how inflated an un-deduplicated count would be.
        dup_counts = dict(con.execute(f"""
            SELECT dup_group, count(*) FROM frames
            WHERE canonical_target IN ({g0.targets_sql()})
              AND (imagetyp IS NULL OR imagetyp LIKE '%Light%')
            GROUP BY dup_group"""))

        with con:
            con.executemany("""INSERT INTO sn_g0_frames VALUES (
                :obs_rowid,:path,:tree,:basename,:night,:jd,:phase_d,:filter,
                :exptime,:readoutm,:xbinning,:naxis1,:naxis2,
                :canonical_target,:raw_object,:observer,:epoch_role,
                :band_role,:dispersion_verdict,:dispersion_status,
                :dispersion_class,:is_image,:pltsolvd,:crval1,:crval2,
                :objctra,:objctdec,:fwhm,:airmass,:zmag,:dup_copies,:era_id)""",
                [{**r,
                  "raw_object": r.get("object"),
                  "dup_copies": dup_counts.get(r.get("dup_group"), 1)}
                 for r in rows])

        # --- the dedup accounting --------------------------------------
        dedup = [("tree:" + (r["tree"] or "?"), r["n_rows"], r["n_canonical"],
                  "catalog rows on this sky in this tree; canonical = the "
                  "copy S0's global (basename, jd) dedup kept")
                 for r in con.execute(g0.DEDUP_SQL)]
        total_rows = sum(d[1] for d in dedup)
        total_canon = sum(d[2] for d in dedup)
        dedup.append(("ALL TREES", total_rows, total_canon,
                      "the global dedup: every tree at once, which is the "
                      "only scope in which the within-rawimage copies are "
                      "visible"))
        # The campaign slice on its own — the number the paper quotes.
        camp = [r for r in rows if r["epoch_role"] == "campaign"]
        dedup.append(("CAMPAIGN 2023 (unique)", len(camp), len(camp),
                      "unique light frames on this sky between nights "
                      f"{g0.CAMPAIGN_FIRST_NIGHT} and "
                      f"{g0.CAMPAIGN_LAST_NIGHT}, alias-merged"))
        # ...and split by which canonical name carried them, because the
        # strategy's published 1,052 counted only one of the two.
        for name in g0.CAMPAIGN_TARGETS:
            n = sum(1 for r in camp if r["canonical_target"] == name)
            dedup.append((f"CAMPAIGN under canonical_target={name}", n, n,
                          "alias-merged slice; a census taken under only "
                          "one of these names misses the other"))
        with con:
            con.executemany("INSERT INTO sn_g0_dedup VALUES (?,?,?,?)", dedup)

        # --- the measured pointing offset ------------------------------
        # For every plate-solved campaign frame, compare the field centre
        # the telescope was COMMANDED to (OBJCTRA/OBJCTDEC) with the one it
        # actually reached (CRVAL1/CRVAL2).  The median of that difference
        # is the systematic the fallback prediction must correct for; the
        # scatter about the median is what sets the search-box size.
        pt = []
        for r in camp:
            if r.get("pltsolvd") != 1 or r.get("crval1") is None:
                continue
            ra = g0.parse_sexagesimal(r.get("objctra"), is_hours=True)
            dec = g0.parse_sexagesimal(r.get("objctdec"), is_hours=False)
            if ra is None or dec is None:
                continue
            # Convert the RA difference to a true angle on the sky.
            dra = (r["crval1"] - ra) * math.cos(math.radians(r["crval2"])) * 3600.0
            ddec = (r["crval2"] - dec) * 3600.0
            pt.append((dra, ddec))
        if pt:
            med_ra = _median([p[0] for p in pt])
            med_dec = _median([p[1] for p in pt])
            scale_arcsec = FALLBACK_SCALE_DEG_PX * 3600.0
            resid = sorted(math.hypot(p[0] - med_ra, p[1] - med_dec) /
                           scale_arcsec for p in pt)
            with con:
                con.execute("INSERT INTO sn_g0_pointing VALUES (?,?,?,?,?,?,?)",
                            ("campaign", len(pt), med_ra, med_dec,
                             _median(resid), _pct(resid, 90), resid[-1]))
            print(f"pointing: commanded->true median offset "
                  f"({med_ra:+.1f}, {med_dec:+.1f}) arcsec over {len(pt)} "
                  f"solved frames; residual about it "
                  f"{_median(resid):.0f} px median, {_pct(resid, 90):.0f} px "
                  f"90th pct, {resid[-1]:.0f} px max")

        with con:
            record_meta(con, code_version=SN_G0_CODE_VERSION,
                        freeze_built_at=utcnow(),
                        sn_position=g0.SN_POSITION_SOURCE,
                        t0_mjd=g0.T0_MJD, t0_source=g0.T0_SOURCE)

        print(f"froze {len(rows):,} unique frames "
              f"({len(camp):,} in the 2023 campaign) from {total_rows:,} "
              f"catalog rows")
        n_spec = sum(1 for r in camp if not r["is_image"])
        print(f"  of the campaign: {len(camp) - n_spec:,} images, "
              f"{n_spec:,} MEASURED spectra (excluded from the imaging "
              f"census by measurement, not by label)")
        return 0
    finally:
        con.close()


def _median(xs):
    """Median of a non-empty sequence (kept local: numpy is not imported in
    the freeze path, which runs without touching a single pixel)."""
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def _pct(sorted_xs, p):
    """Percentile of an ALREADY SORTED sequence, nearest-rank."""
    if not sorted_xs:
        return None
    k = max(0, min(len(sorted_xs) - 1,
                   int(round(p / 100.0 * (len(sorted_xs) - 1)))))
    return sorted_xs[k]


# ===========================================================================
# 0b — CENSUS: build the queue
# ===========================================================================
CENSUS_DDL = """
CREATE TABLE IF NOT EXISTS sn_g0_census (
    obs_rowid    INTEGER PRIMARY KEY,
    path         TEXT NOT NULL,
    tree         TEXT,
    night        TEXT,
    filter       TEXT,
    exptime      REAL,
    phase_d      REAL,
    epoch_role   TEXT,
    band_role    TEXT,
    dispersion_class TEXT,
    readoutm     TEXT,
    status       TEXT NOT NULL DEFAULT 'pending',
    quality      TEXT,
    pred_x       REAL,
    pred_y       REAL,
    sn_x         REAL,
    sn_y         REAL,
    offset_px    REAL,
    peak_adu     REAL,
    core_max_adu REAL,
    box_max_adu  REAL,
    frame_max_adu REAL,
    sky_adu      REAL,
    isolation_px REAL,
    bound_basis  TEXT,
    n_ge_reject  INTEGER,
    n_ge_clip    INTEGER,
    saturation_class TEXT,
    measure_s    REAL,
    error        TEXT,
    code_version TEXT,
    measured_at  TEXT
);
"""


def cmd_census(args) -> int:
    """0b step 1 — enqueue the frames whose pixels must be read.

    Everything on this sky is queued, at every epoch, including the frames
    S2c measures as spectra: the grism series has its own saturation
    question (a 256 s slitless exposure of an 11th-magnitude supernova can
    clip its own trace), and refusing to look at it would leave Gate 0c
    answering that question from an exposure-time argument instead of from
    the pixels.  The census records what each row IS, so a query can always
    take the imaging-only subset back out.
    """
    con = connect(args.manifest)
    try:
        con.executescript(CENSUS_DDL)
        if args.rebuild:
            with con:
                con.execute("DELETE FROM sn_g0_census")
        done = {r[0] for r in con.execute(
            "SELECT obs_rowid FROM sn_g0_census WHERE status != 'pending'")}
        rows = [dict(r) for r in con.execute("""
            SELECT obs_rowid, path, tree, night, filter, exptime, phase_d,
                   epoch_role, band_role, dispersion_class, readoutm
            FROM sn_g0_frames ORDER BY jd""")]
        new = [r for r in rows if r["obs_rowid"] not in done]
        with con:
            con.executemany("""INSERT OR IGNORE INTO sn_g0_census
                (obs_rowid, path, tree, night, filter, exptime, phase_d,
                 epoch_role, band_role, dispersion_class, readoutm, status)
                VALUES (:obs_rowid,:path,:tree,:night,:filter,:exptime,
                        :phase_d,:epoch_role,:band_role,:dispersion_class,
                        :readoutm,'pending')""", new)
        pending = con.execute("SELECT count(*) FROM sn_g0_census "
                              "WHERE status = 'pending'").fetchone()[0]
        print(f"census queue: {len(rows):,} frames total, "
              f"{len(done):,} already measured, {pending:,} pending")
        return 0
    finally:
        con.close()


# ===========================================================================
# 0b — CENSUS: the measurement itself
# ===========================================================================
def _measure_one(task: dict) -> dict:
    """Open ONE frame and measure the supernova's peak.  Runs in a worker.

    Everything this function needs arrives in ``task`` — it must not touch
    the database, because a process pool sharing one sqlite handle is a
    corruption bug waiting for a slow disk.

    THE MEASUREMENT, STEP BY STEP
    1. open the file read-only and take the last HDU (these are fpacked
       images: HDU 0 is an empty stub, HDU 1 carries the pixels);
    2. if the header holds a celestial WCS, project the supernova's
       catalogue position through it — this is the exact answer;
    3. otherwise project through a north-up tangent plane centred on the
       COMMANDED pointing, corrected by the campaign's measured commanded-
       versus-true offset — this is a seed, not an answer;
    4. take the maximum of a small core stamp and of a large search box;
    5. record where in the core stamp the maximum fell, which is the check
       that decides whether the frame's own WCS is corroborated;
    6. for frames with a WCS, also measure the ISOLATION radius: how far
       away the nearest BRIGHTER pixel is.  That single number is what
       later tells the census how much a box maximum is worth in each band,
       measured rather than assumed.
    """
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS

    t0 = time.time()
    out = {"obs_rowid": task["obs_rowid"], "status": "measured",
           "error": None, "quality": None, "pred_x": None, "pred_y": None,
           "sn_x": None, "sn_y": None, "offset_px": None, "peak_adu": None,
           "core_max_adu": None, "box_max_adu": None, "frame_max_adu": None,
           "sky_adu": None, "isolation_px": None, "bound_basis": None,
           "n_ge_reject": None, "n_ge_clip": None}
    try:
        with warnings.catch_warnings():
            # Astropy warns loudly about the non-standard MaxIm cards in
            # every one of these headers; the warning is true and useless.
            warnings.simplefilter("ignore")
            with fits.open(task["abs_path"], memmap=False) as hdul:
                header = hdul[-1].header
                data = hdul[-1].data
            if data is None or data.ndim != 2:
                out["status"] = "unreadable"
                out["error"] = "no 2-D image data"
                return out
            ny, nx = data.shape

            # --- step 2/3: where is the supernova on this frame? --------
            has_wcs = False
            if header.get("CTYPE1"):
                try:
                    wcs = WCS(header)
                    if wcs.has_celestial:
                        px, py = wcs.all_world2pix(g0.SN_RA_DEG,
                                                   g0.SN_DEC_DEG, 0)
                        px, py = float(px), float(py)
                        has_wcs = True
                except Exception:                        # noqa: BLE001
                    has_wcs = False
            if not has_wcs:
                ra = g0.parse_sexagesimal(header.get("OBJCTRA"), True)
                dec = g0.parse_sexagesimal(header.get("OBJCTDEC"), False)
                if ra is None or dec is None:
                    out["status"] = "unreadable"
                    out["error"] = "no WCS and no pointing cards"
                    return out
                # Apply the campaign's MEASURED commanded->true offset, so
                # the search box is centred on where the telescope actually
                # goes rather than on where it was told to go.
                ra += task["d_ra_deg"] / max(1e-9,
                                             math.cos(math.radians(dec)))
                dec += task["d_dec_deg"]
                scale = abs(float(header.get("CDELT1")
                                  or FALLBACK_SCALE_DEG_PX))
                px, py = g0.gnomonic_pixel(g0.SN_RA_DEG, g0.SN_DEC_DEG,
                                           ra, dec, scale,
                                           nx / 2.0, ny / 2.0)
        out["pred_x"], out["pred_y"] = px, py

        xi, yi = int(round(px)), int(round(py))
        if not (0 <= xi < nx and 0 <= yi < ny):
            # The supernova's position does not fall on this detector.  Not
            # an error — the campaign contains a handful of frames pointed
            # elsewhere on the galaxy — but nothing can be measured.
            out["status"] = "off_frame"
            out["quality"] = "off_frame"
            return out

        def window(half):
            """Clipped slice of `half`-pixel half-width about (xi, yi)."""
            y0, y1 = max(0, yi - half), min(ny, yi + half + 1)
            x0, x1 = max(0, xi - half), min(nx, xi + half + 1)
            return data[y0:y1, x0:x1], x0, y0

        core, cx0, cy0 = window(g0.CORE_HALF_PX)
        box, bx0, by0 = window(g0.BOUND_HALF_PX)
        sky, _, _ = window(g0.SKY_HALF_PX)

        core_max = float(core.max())
        box_max = float(box.max())
        out["core_max_adu"] = core_max
        out["box_max_adu"] = box_max
        out["frame_max_adu"] = float(data.max())
        out["sky_adu"] = float(np.median(sky))

        # --- step 5: did the core stamp's peak land on the prediction? --
        j = np.unravel_index(int(core.argmax()), core.shape)
        sx, sy = j[1] + cx0, j[0] + cy0
        out["sn_x"], out["sn_y"] = float(sx), float(sy)
        offset = math.hypot(sx - px, sy - py)
        out["offset_px"] = offset
        out["quality"] = g0.position_quality(offset, has_wcs)

        # WHICH NUMBER GETS PUBLISHED AS THE PEAK.
        # With a plate solution the supernova is inside the core stamp
        # whether or not it dominates it, so the core maximum is the right
        # number in both WCS cases: a measurement when the peak is centred,
        # a tight upper bound when it is not.  Only a frame with no solution
        # at all has to fall back on the wide box, whose maximum is a much
        # looser bound — and the census labels it as such rather than
        # letting a 501 px box maximum masquerade as a supernova peak.
        out["peak_adu"] = core_max if has_wcs else box_max
        out["bound_basis"] = ("measured" if out["quality"] == "wcs"
                              else ("core" if has_wcs else "box"))

        # --- step 6: isolation, measured only where truth is known ------
        if out["quality"] == "wcs":
            brighter = np.argwhere(box > core_max)
            if len(brighter):
                dy = brighter[:, 0] + by0 - py
                dx = brighter[:, 1] + bx0 - px
                out["isolation_px"] = float(np.min(np.hypot(dx, dy)))
            else:
                # Nothing in a 501 px box outbrightens the supernova: the
                # isolation radius is at least the box half-width, and is
                # recorded as that floor.
                out["isolation_px"] = float(g0.BOUND_HALF_PX)
            out["n_ge_reject"] = int((core >= task["reject_adu"]).sum())
            out["n_ge_clip"] = int((core >= task["clip_adu"]).sum())
    except Exception as exc:                              # noqa: BLE001
        out["status"] = "unreadable"
        out["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        out["measure_s"] = time.time() - t0
    return out


def cmd_measure(args) -> int:
    """0b step 2 — measure the pending frames.  Resumable."""
    con = connect(args.manifest)
    try:
        # The screen comes from S2's MEASURED ceiling.  If it is missing the
        # census refuses to run rather than falling back to a typed number:
        # a screen with nothing behind it is the exact state that had task
        # SN-G0b marked BLOCKED.
        screens = _screens(con)
        pointing = con.execute(
            "SELECT d_ra_arcsec, d_dec_arcsec FROM sn_g0_pointing "
            "WHERE scope = 'campaign'").fetchone()
        d_ra = (pointing["d_ra_arcsec"] / 3600.0) if pointing else 0.0
        d_dec = (pointing["d_dec_arcsec"] / 3600.0) if pointing else 0.0

        rows = con.execute("""
            SELECT c.obs_rowid, c.path, c.readoutm FROM sn_g0_census c
            WHERE c.status = 'pending' ORDER BY c.obs_rowid""").fetchall()
        if args.limit:
            rows = rows[:args.limit]
        if not rows:
            print("nothing pending.")
            return 0

        tasks = []
        for r in rows:
            sc = screens.get(r["readoutm"]) or screens.get("High Gain")
            tasks.append({"obs_rowid": r["obs_rowid"],
                          "abs_path": str(args.archive / r["path"]),
                          "d_ra_deg": d_ra, "d_dec_deg": d_dec,
                          "reject_adu": sc.reject_adu,
                          "clip_adu": sc.clip_adu})

        print(f"measuring {len(tasks):,} frames with {args.workers} workers "
              f"(screen: reject >= {tasks[0]['reject_adu']} ADU)")
        started = time.time()
        done = 0
        batch = []
        with concurrent.futures.ProcessPoolExecutor(args.workers) as pool:
            for res in pool.map(_measure_one, tasks, chunksize=4):
                batch.append(res)
                done += 1
                if len(batch) >= 40:
                    _flush(con, batch, screens)
                    batch = []
                    print(f"  {done:,}/{len(tasks):,} "
                          f"({time.time() - started:.0f}s)", flush=True)
                if args.max_seconds and time.time() - started > args.max_seconds:
                    break
        if batch:
            _flush(con, batch, screens)
        with con:
            record_meta(con, census_measured_at=utcnow(),
                        code_version=SN_G0_CODE_VERSION)
        left = con.execute("SELECT count(*) FROM sn_g0_census "
                           "WHERE status = 'pending'").fetchone()[0]
        print(f"measured {done:,} frames in {time.time() - started:.0f}s; "
              f"{left:,} still pending")
        return 0
    finally:
        con.close()


def _screens(con) -> dict:
    """Per-readout-mode screens, built from S2's measured ceilings."""
    rows = con.execute("SELECT mode, clip_adu, veto_adu FROM "
                       "s2_ceiling_modes").fetchall()
    if not rows:
        raise SystemExit(
            "s2_ceiling_modes is empty: the saturation screen has no "
            "measured clip behind it.  Re-run stage S2 "
            "(`run_s2_campaign.py ceiling` then `params`) first.")
    out = {}
    for r in rows:
        if r["clip_adu"] is None:
            continue
        out[r["mode"]] = g0.screen_for_mode(r["mode"], r["clip_adu"],
                                            r["veto_adu"])
    return out


def _flush(con, batch, screens) -> None:
    """Write one batch of results, classifying each as it goes.

    The saturation class is computed HERE rather than in the worker because
    it depends on the screen, and a screen change must be re-appliable from
    the stored ADU values without re-reading a single pixel — the same
    measure-first-judge-later discipline S2c uses.
    """
    rows = []
    for r in batch:
        sc = screens.get(_mode_of(con, r["obs_rowid"])) or screens["High Gain"]
        cls = None
        if r["status"] == "measured":
            cls = g0.saturation_class(r["peak_adu"], r["quality"], sc)
        rows.append({**r, "saturation_class": cls,
                     "code_version": SN_G0_CODE_VERSION,
                     "measured_at": utcnow()})
    with con:
        con.executemany("""UPDATE sn_g0_census SET
            status=:status, quality=:quality, pred_x=:pred_x, pred_y=:pred_y,
            sn_x=:sn_x, sn_y=:sn_y, offset_px=:offset_px, peak_adu=:peak_adu,
            core_max_adu=:core_max_adu, box_max_adu=:box_max_adu,
            frame_max_adu=:frame_max_adu, sky_adu=:sky_adu,
            isolation_px=:isolation_px, bound_basis=:bound_basis,
            n_ge_reject=:n_ge_reject,
            n_ge_clip=:n_ge_clip, saturation_class=:saturation_class,
            measure_s=:measure_s, error=:error, code_version=:code_version,
            measured_at=:measured_at
            WHERE obs_rowid=:obs_rowid""", rows)


_MODE_CACHE: dict = {}


def _mode_of(con, obs_rowid):
    """Readout mode for one frame, cached (the census asks once per row)."""
    if not _MODE_CACHE:
        for r in con.execute("SELECT obs_rowid, readoutm FROM sn_g0_census"):
            _MODE_CACHE[r["obs_rowid"]] = r["readoutm"]
    return _MODE_CACHE.get(obs_rowid)


# ===========================================================================
# 0b — the saturation matrix
# ===========================================================================
MATRIX_DDL = """
DROP TABLE IF EXISTS sn_g0_matrix;
CREATE TABLE sn_g0_matrix (
    night TEXT, filter TEXT, phase_d REAL, band_role TEXT,
    n_frames INTEGER, n_measured INTEGER,
    n_clean INTEGER, n_suspect INTEGER, n_rejected INTEGER,
    n_bounded_clean INTEGER, n_undetermined INTEGER, n_spectra INTEGER,
    max_peak_adu REAL, min_peak_adu REAL, min_exptime REAL, max_exptime REAL
);
DROP TABLE IF EXISTS sn_g0_bands;
CREATE TABLE sn_g0_bands (
    band_role TEXT, filter TEXT,
    n_frames INTEGER, n_images INTEGER, n_spectra INTEGER,
    n_wcs INTEGER, n_bound INTEGER,
    n_clean INTEGER, n_suspect INTEGER, n_rejected INTEGER,
    n_bounded_clean INTEGER, n_undetermined INTEGER,
    n_usable INTEGER, first_clean_night TEXT, first_clean_phase_d REAL,
    isolation_false_id INTEGER, isolation_tested INTEGER
);
"""


def cmd_matrix(args) -> int:
    """0b step 3 — fold the per-frame census into the two published tables.

    ``sn_g0_matrix`` is the strategy's filter × night saturation matrix.
    ``sn_g0_bands`` is the per-band summary that answers question 1, and it
    carries the isolation calibration beside the counts so that the worth of
    every ``bounded_clean`` row can be read off the same line.
    """
    con = connect(args.manifest)
    try:
        con.executescript(MATRIX_DDL)
        cells = con.execute("""
            SELECT night, filter, band_role,
                   min(phase_d) AS phase_d, count(*) AS n_frames,
                   sum(status = 'measured') AS n_measured,
                   sum(saturation_class = 'clean') AS n_clean,
                   sum(saturation_class = 'suspect') AS n_suspect,
                   sum(saturation_class = 'rejected') AS n_rejected,
                   sum(saturation_class = 'bounded_clean') AS n_bc,
                   sum(saturation_class = 'undetermined') AS n_und,
                   sum(dispersion_class = 'dispersed') AS n_spectra,
                   max(peak_adu) AS max_peak, min(peak_adu) AS min_peak,
                   min(exptime) AS min_exp, max(exptime) AS max_exp
            FROM sn_g0_census WHERE epoch_role = 'campaign'
            GROUP BY night, filter ORDER BY night, filter""").fetchall()
        with con:
            con.executemany("INSERT INTO sn_g0_matrix VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            [(c["night"], c["filter"], c["phase_d"],
                              c["band_role"], c["n_frames"], c["n_measured"],
                              c["n_clean"], c["n_suspect"], c["n_rejected"],
                              c["n_bc"], c["n_und"], c["n_spectra"],
                              c["max_peak"], c["min_peak"],
                              c["min_exp"], c["max_exp"]) for c in cells])

        # --- the per-band summary --------------------------------------
        bands = []
        for r in con.execute("""
                SELECT band_role, filter FROM sn_g0_census
                WHERE epoch_role = 'campaign'
                GROUP BY band_role, filter ORDER BY band_role, filter"""):
            br, filt = r["band_role"], r["filter"]
            rows = con.execute("""
                SELECT * FROM sn_g0_census
                WHERE epoch_role = 'campaign' AND filter = ?
                ORDER BY night, phase_d""", (filt,)).fetchall()
            n_images = sum(1 for x in rows
                           if x["dispersion_class"] != "dispersed")
            counts = {k: sum(1 for x in rows if x["saturation_class"] == k)
                      for k in ("clean", "suspect", "rejected",
                                "bounded_clean", "undetermined")}
            usable = sum(1 for x in rows if g0.is_usable_photometry(dict(x)))
            # First night at which this filter delivers a clean measurement
            # of the supernova — the "true clean start per band" the
            # strategy asks the census to decide.
            first = next((x for x in rows
                          if x["saturation_class"] == "clean"), None)
            iso = [x["isolation_px"] for x in rows]
            wrong, tested = g0.isolation_false_id_rate(iso, g0.BOUND_HALF_PX)
            bands.append((br, filt, len(rows), n_images, len(rows) - n_images,
                          sum(1 for x in rows if x["quality"] == "wcs"),
                          sum(1 for x in rows if x["quality"] == "bound"),
                          counts["clean"], counts["suspect"],
                          counts["rejected"], counts["bounded_clean"],
                          counts["undetermined"], usable,
                          first["night"] if first else None,
                          first["phase_d"] if first else None,
                          wrong, tested))
        with con:
            con.executemany("INSERT INTO sn_g0_bands VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", bands)
            record_meta(con, matrix_built_at=utcnow())
        print(f"matrix: {len(cells)} filter x night cells, "
              f"{len(bands)} band rows")
        for b in bands:
            print(f"  {b[1]:2s} ({b[0]:10s}) n={b[2]:4d} images={b[3]:4d} "
                  f"clean={b[7]:4d} suspect={b[8]:3d} rejected={b[9]:3d} "
                  f"bounded={b[10]:3d} undet={b[11]:3d} USABLE={b[12]:4d}")
        return 0
    finally:
        con.close()


# ===========================================================================
# 0c — the grism triage
# ===========================================================================
TRIAGE_DDL = """
DROP TABLE IF EXISTS sn_g0_triage;
CREATE TABLE sn_g0_triage (
    night TEXT PRIMARY KEY, phase_d REAL,
    n_labelled INTEGER, n_dispersed INTEGER, n_direct INTEGER,
    n_indeterminate INTEGER, min_exptime REAL, max_exptime REAL,
    n_paired_direct INTEGER, paired_filters TEXT,
    max_box_adu REAL, n_at_screen INTEGER, in_flash_window INTEGER
);
DROP TABLE IF EXISTS sn_g0_triage_summary;
CREATE TABLE sn_g0_triage_summary (
    clause TEXT PRIMARY KEY, requirement TEXT, value TEXT, passed INTEGER
);
"""


def cmd_triage(args) -> int:
    """0c — the grism triage, on the measurement rather than the label."""
    con = connect(args.manifest)
    try:
        con.executescript(TRIAGE_DDL)
        screens = _screens(con)
        screen = screens.get("High Gain")

        nights = con.execute("""
            SELECT night, min(phase_d) AS phase_d, count(*) AS n_labelled,
                   sum(dispersion_class = 'dispersed') AS n_disp,
                   sum(dispersion_class = 'direct') AS n_dir,
                   sum(dispersion_class = 'indeterminate') AS n_ind,
                   min(exptime) AS min_exp, max(exptime) AS max_exp,
                   max(box_max_adu) AS max_box,
                   sum(box_max_adu >= ?) AS n_at_screen
            FROM sn_g0_census
            WHERE epoch_role = 'campaign' AND filter = '6'
            GROUP BY night ORDER BY night""",
            (screen.reject_adu,)).fetchall()

        rows = []
        for n in nights:
            # Does this night ALSO hold direct images?  Gate 0c requirement
            # (i): a slitless trace can only be identified if something
            # tells you where the zero order should be, and the strategy
            # insists this be VERIFIED from the manifest rather than
            # assumed.
            paired = con.execute("""
                SELECT filter, count(*) FROM sn_g0_census
                WHERE night = ? AND epoch_role = 'campaign'
                  AND dispersion_class != 'dispersed' AND filter != '6'
                GROUP BY filter ORDER BY filter""", (n["night"],)).fetchall()
            rows.append((n["night"], n["phase_d"], n["n_labelled"],
                         n["n_disp"], n["n_dir"], n["n_ind"],
                         n["min_exp"], n["max_exp"],
                         sum(p[1] for p in paired),
                         ",".join(f"{p[0]}x{p[1]}" for p in paired),
                         n["max_box"], n["n_at_screen"],
                         1 if (n["phase_d"] is not None
                               and n["phase_d"] <= g0.FLASH_PHASE_END_D)
                         else 0))
        with con:
            con.executemany("INSERT INTO sn_g0_triage VALUES "
                            "(?,?,?,?,?,?,?,?,?,?,?,?,?)", rows)

        # --- assemble the series and apply the pre-registered rule ------
        disp_nights = [r for r in rows if r[3] > 0]
        series = g0.GrismSeries(
            n_labelled=sum(r[2] for r in rows),
            n_dispersed=sum(r[3] for r in rows),
            n_direct=sum(r[4] for r in rows),
            n_indeterminate=sum(r[5] for r in rows),
            n_nights=len(disp_nights),
            n_flash_nights=sum(1 for r in disp_nights if r[12]),
            n_nights_with_paired_direct=sum(1 for r in disp_nights if r[8] > 0),
            # The extraction has not been run: this project consumes
            # extracted spectra, it does not own the extractor (strategy
            # §9), and no grism extraction for this target exists in the
            # manifest.  Counted, not assumed — see the query below.
            n_extracted=_n_extracted(con),
            n_contamination_passed=0,
            wavelength_source="",
            n_flats=_n_grism_flats(con))
        decision = g0.grism_promotion(series)
        with con:
            con.executemany(
                "INSERT INTO sn_g0_triage_summary VALUES (?,?,?,?)",
                [(k, v["requirement"], str(v["value"]), int(v["passed"]))
                 for k, v in decision["clauses"].items()])
            record_meta(con, triage_built_at=utcnow(),
                        grism_promoted=int(decision["promoted"]),
                        grism_blocking=",".join(decision["blocking"]),
                        grism_n_dispersed=series.n_dispersed,
                        grism_n_nights=series.n_nights,
                        grism_n_flash_nights=series.n_flash_nights,
                        grism_n_flats=series.n_flats)
        print(f"grism series: {series.n_dispersed} measured spectra over "
              f"{series.n_nights} nights ({series.n_flash_nights} inside the "
              f"flash window), {series.n_direct} measured direct images and "
              f"{series.n_indeterminate} uncertifiable frames wearing the "
              f"same label")
        print(f"promotion: {'YES' if decision['promoted'] else 'NO'} "
              f"— blocking clauses: {', '.join(decision['blocking']) or 'none'}")
        return 0
    finally:
        con.close()


def _n_extracted(con) -> int:
    """How many 1-D spectra of this target has the grism pipeline
    produced?  Counted from ``g_extractions`` if that table exists, so the
    triage's "nothing extracted yet" is a query and not a belief."""
    try:
        return con.execute("""
            SELECT count(*) FROM g_extractions e
            JOIN frames f ON f.obs_rowid = e.obs_rowid
            WHERE f.canonical_target IN ('2023ixf', 'M101')""").fetchone()[0]
    except sqlite3.Error:
        return 0


def _n_grism_flats(con) -> int:
    """Flat frames taken through filter '6', at any epoch.  The strategy
    asserts there are none and makes relative-flux-only policy from it;
    Gate 0 re-checks the assertion instead of inheriting it."""
    try:
        return con.execute(
            "SELECT count(*) FROM frames WHERE filter = '6' "
            "AND imagetyp LIKE '%Flat%'").fetchone()[0]
    except sqlite3.Error:
        return 0


# ===========================================================================
# THE VERDICTS
# ===========================================================================
VERDICT_DDL = """
DROP TABLE IF EXISTS sn_g0_verdict;
CREATE TABLE sn_g0_verdict (
    question_id TEXT PRIMARY KEY, question TEXT, deciding_number TEXT,
    value REAL, verdict TEXT, moved INTEGER, basis TEXT
);
DROP TABLE IF EXISTS sn_g0_astrometry;
CREATE TABLE sn_g0_astrometry (
    scope TEXT PRIMARY KEY, description TEXT,
    k INTEGER, n INTEGER, n_population INTEGER,
    rate_pct REAL, wilson_lo REAL, wilson_hi REAL, verdict TEXT, basis TEXT
);
"""


def cmd_verdicts(args) -> int:
    """Answer the questions Gate 0 exists to answer, plus the venue posture.

    Every row stores the DECIDING NUMBER as a string beside the verdict, so
    a reader never has to reconstruct which quantity a word like "NO-GO" was
    about — and so a re-run that keeps the word while moving the number is
    visible in the table rather than only on the page.
    """
    con = connect(args.manifest)
    try:
        con.executescript(VERDICT_DDL)
        v = []
        astro = []

        # --- Q1: how many frames are genuinely usable? ------------------
        b = con.execute("""
            SELECT sum(n_usable) AS usable, sum(n_frames) AS total,
                   sum(n_spectra) AS spectra, sum(n_rejected) AS rejected,
                   sum(n_suspect) AS suspect, sum(n_undetermined) AS undet,
                   sum(n_clean) AS clean, sum(n_bounded_clean) AS bounded
            FROM sn_g0_bands WHERE band_role = 'broadband'""").fetchone()
        first = con.execute("""
            SELECT min(first_clean_night) AS night,
                   min(first_clean_phase_d) AS phase
            FROM sn_g0_bands WHERE band_role = 'broadband'""").fetchone()
        v.append(("usable-broadband",
                  "How many SN frames are genuinely usable broadband "
                  "photometry, after removing spectra and saturated epochs?",
                  f"{b['usable']} of {b['total']} campaign broadband frames "
                  f"({b['clean']} measured clean + {b['bounded']} clean by "
                  f"bound)",
                  float(b["usable"]), "MEASURED", 0,
                  f"removed: {b['rejected']} over the peak-ADU screen, "
                  f"{b['suspect']} suspect (held pending the linearity "
                  f"curve), {b['undet']} undetermined (no plate solution and "
                  f"a search box that reaches the screen); the first clean "
                  f"broadband frame is night {first['night']} "
                  f"(+{first['phase']:.1f} d)"))

        # --- Q2a: the S1 stratum verdict, recomputed --------------------
        # This is the number the strategy graded NO-GO on, and the one the
        # dispersion repair was expected to move.  It is read here rather
        # than recomputed so that the two pages cannot disagree.
        s1 = con.execute("""
            SELECT s.n_population, s.n_sample, s.n_population_label_gate,
                   sum(e.status = 'solved') AS k
            FROM s1_strata s LEFT JOIN s1_solve_experiment e
              ON e.stratum_id = s.stratum_id
            WHERE s.stratum_id = 'sn_gsense_broadband'""").fetchone()
        k, n = int(s1["k"] or 0), int(s1["n_sample"] or 0)
        lo, hi = g0.wilson_ci(k, n)
        s1_verdict = g0.astrometry_verdict(k, n, s1["n_population"])
        astro.append(("s1-stratum",
                      "S1 stratum sn_gsense_broadband: can the UNSOLVED "
                      "residue of this campaign be batch-solved blind?",
                      k, n, s1["n_population"], g0.rate_pct(k, n), lo, hi,
                      s1_verdict,
                      f"sample of {n} drawn from a population of "
                      f"{s1['n_population']} after the measured-dispersion "
                      f"gate replaced the filter-label gate (the label gate "
                      f"had counted {s1['n_population_label_gate']})"))

        # --- Q2b: what the stratum is actually MADE of ------------------
        # The stratum's name says "broadband"; its population does not.
        # Decomposing it is what turns "the SN campaign is NO-GO" into a
        # statement about which frames.
        comp = con.execute("""
            SELECT f.band_role, count(*) AS n
            FROM frames fr
            JOIN sn_g0_frames f ON f.obs_rowid = fr.obs_rowid
            LEFT JOIN frame_dispersion d ON d.obs_rowid = fr.obs_rowid
            WHERE fr.canonical_target = '2023ixf' AND fr.is_canonical = 1
              AND fr.tree = 'rawimage'
              AND (fr.pltsolvd IS NULL OR fr.pltsolvd != 1)
              AND (d.verdict IS NULL OR d.verdict != 'dispersed')
            GROUP BY 1""").fetchall()
        comp_d = {r["band_role"]: r["n"] for r in comp}

        # --- Q2c: the campaign's ACTUAL astrometric coverage ------------
        # A different and, for a light curve, more useful question: of the
        # campaign's broadband images, how many already carry a plate
        # solution the census has independently CORROBORATED by finding the
        # supernova within 6 px of its catalogue position?  That is a
        # census, not a sample, so the exact rate is judged.
        for role in ("broadband", "narrowband"):
            r = con.execute("""
                SELECT sum(quality = 'wcs') AS k, count(*) AS n
                FROM sn_g0_census
                WHERE epoch_role = 'campaign' AND band_role = ?
                  AND dispersion_class != 'dispersed' AND tree = ?""",
                (role, g0.SCIENCE_TREE)).fetchone()
            kk, nn = int(r["k"] or 0), int(r["n"] or 0)
            lo2, hi2 = g0.wilson_ci(kk, nn)
            astro.append((f"campaign-{role}",
                          f"Campaign {role} images carrying a plate "
                          f"solution the census independently corroborates",
                          kk, nn, nn, g0.rate_pct(kk, nn), lo2, hi2,
                          g0.astrometry_verdict(kk, nn, nn),
                          "corroboration = the brightest pixel in a 21 px "
                          "stamp at the catalogue position of the supernova "
                          f"lands within {g0.MAX_CENTROID_OFFSET_PX:.0f} px "
                          "of it; this is a full census, so the exact rate "
                          "is judged and no sampling interval applies"))

        camp_bb = next(a for a in astro if a[0] == "campaign-broadband")
        v.append(("astrometry",
                  "Does the astrometry verdict move now that the rate is "
                  "computed over real images?",
                  f"{k}/{n} = {g0.rate_pct(k, n):.1f}% on the S1 stratum "
                  f"(unchanged verdict {s1_verdict}); "
                  f"{camp_bb[2]}/{camp_bb[3]} = {camp_bb[5]:.1f}% of "
                  f"campaign broadband images carry a corroborated solution "
                  f"({camp_bb[8]})",
                  float(g0.rate_pct(k, n) or 0.0), s1_verdict, 0,
                  f"the stratum's unsolved population is "
                  f"{comp_d.get('narrowband', 0)} narrowband, "
                  f"{comp_d.get('broadband', 0)} broadband and "
                  f"{comp_d.get('other', 0)} other frames, so its NO-GO is a "
                  f"verdict about long narrowband exposures and not about "
                  f"the gri light curve"))

        # --- Q3: is the dispersed series a flash-phase record? ----------
        meta = dict(con.execute("SELECT key, value FROM sn_g0_build_meta"))
        promoted = meta.get("grism_promoted") == "1"
        n_extracted = _n_extracted(con)
        v.append(("grism",
                  "Is the dispersed series a genuine flash-phase spectral "
                  "record, and does it promote?",
                  f"{meta.get('grism_n_dispersed', '?')} measured spectra on "
                  f"{meta.get('grism_n_nights', '?')} nights, "
                  f"{meta.get('grism_n_flash_nights', '?')} of them inside "
                  f"the flash window; {n_extracted} extracted spectra",
                  float(meta.get("grism_n_flash_nights", 0)),
                  "PROMOTED" if promoted else "NOT PROMOTED", 0,
                  f"blocking clauses: {meta.get('grism_blocking', '')}; "
                  f"{meta.get('grism_n_flats', '?')} filter-'6' flat frames "
                  f"exist in the archive, so the relative-flux-only policy "
                  f"stands"))

        # --- the venue posture ------------------------------------------
        # The narrowband bandpass has not been recovered: task
        # SN-S1-narrowband-curves is BLOCKED on a manufacturer sheet or a
        # physical filter scan, and no transmission curve exists anywhere in
        # this repository or this manifest.
        posture = g0.venue_posture(grism_promoted=promoted,
                                   bandpass_recovered=False)
        v.append(("venue",
                  "Does Gate 0 promote the venue posture from AJ/PASP to "
                  "ApJ, per the strategy's own rule?",
                  f"{n_extracted} contamination-tested flash-phase spectra "
                  f"(the pre-registered rule requires "
                  f"{g0.GRISM_PROMOTION_MIN_NIGHTS}) and 0 recovered "
                  f"narrowband transmission curves (the rule requires 1)",
                  float(n_extracted), posture["posture"],
                  int(posture["moved"]),
                  "the strategy decides the posture in advance: AJ/PASP is "
                  "the base case, ApJ only if Gate 0 promotes the grism or "
                  "recovers the narrowband bandpass; neither happened, so "
                  "the posture stays where it was put"))

        with con:
            con.executemany("INSERT INTO sn_g0_verdict VALUES (?,?,?,?,?,?,?)",
                            v)
            con.executemany("INSERT INTO sn_g0_astrometry VALUES "
                            "(?,?,?,?,?,?,?,?,?,?)", astro)
            record_meta(con, verdicts_built_at=utcnow())
        for row in v:
            print(f"\n[{row[4]}] {row[0]}")
            print(f"  {row[1]}")
            print(f"  deciding number: {row[2]}")
        print()
        for a in astro:
            print(f"  {a[0]:22s} {a[2]:5d}/{a[3]:<5d} = {a[5]:5.1f}%  "
                  f"{a[8]}")
        return 0
    finally:
        con.close()


# ===========================================================================
def cmd_status(args) -> int:
    con = connect(args.manifest, read_only=True)
    try:
        for name in ("sn_g0_frames", "sn_g0_census", "sn_g0_matrix",
                     "sn_g0_bands", "sn_g0_triage", "sn_g0_verdict"):
            try:
                n = con.execute(f"SELECT count(*) FROM {name}").fetchone()[0]
            except sqlite3.Error:
                n = "MISSING"
            print(f"{name:20s} {n}")
        try:
            for r in con.execute("SELECT status, count(*) FROM sn_g0_census "
                                 "GROUP BY status ORDER BY 2 DESC"):
                print(f"  census.{r[0]:12s} {r[1]:,}")
        except sqlite3.Error:
            pass
        return 0
    finally:
        con.close()


def cmd_report(args) -> int:
    from macro_sn import report_gate0
    path = report_gate0.render_report(args.manifest)
    print(f"wrote {path.relative_to(REPO_ROOT)}")
    return 0


# ===========================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("freeze", help="0a: the deduplicated manifest freeze")
    p.set_defaults(func=cmd_freeze)

    p = sub.add_parser("census", help="0b: enqueue frames for measurement")
    p.add_argument("--rebuild", action="store_true",
                   help="discard existing measurements and start over")
    p.set_defaults(func=cmd_census)

    p = sub.add_parser("measure", help="0b: measure pending frames (resumable)")
    p.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--max-seconds", type=float, default=0.0)
    p.set_defaults(func=cmd_measure)

    p = sub.add_parser("matrix", help="0b: the filter x night matrix")
    p.set_defaults(func=cmd_matrix)

    p = sub.add_parser("triage", help="0c: the grism triage")
    p.set_defaults(func=cmd_triage)

    p = sub.add_parser("verdicts", help="the three Gate 0 answers")
    p.set_defaults(func=cmd_verdicts)

    p = sub.add_parser("report", help="render the Gate 0 evidence page")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("status", help="progress (read-only)")
    p.set_defaults(func=cmd_status)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
