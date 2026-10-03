#!/usr/bin/env python
"""DW-N1 / DW-P36-0 — what is already known about every dwarf field, and
whether RLMT's H-alpha frames are deep enough to add anything.

WHAT THIS SCRIPT DOES
---------------------
The committee review of 2026-10-03 (``committee/reviews/2026-10-03``) put a
novelty gate in front of the dwarf-galaxy project (journal-editor, DW-N1:
"per-candidate table of existing H-alpha/spectroscopic status; if most are
already classified, stop") and a depth gate in front of the H-alpha stacking
(physicist, DW-P36-0: "predicted limit versus expected L_Ha per candidate
tabulated; fields where the limit cannot separate dIrr from dSph are
labelled uninformative").  This script produces both tables, and the one
figure that shows them, from four inputs:

    1. the S0 manifest (opened READ-ONLY)        -> which frames exist
    2. ``literature_status.csv`` (beside this file) -> what is published
    3. LVGDB object pages + GALEX AIS cone search -> machine-read, cached
       under ``sources/`` so the run is reproducible offline
    4. the raw H-alpha frames in the archive (READ-ONLY) -> the measured
       sky noise, hence a measured (not asserted) 3-sigma flux limit

Nothing in the output tables is typed: literature values enter only through
``literature_status.csv`` (each with its reference key, resolved in
``references.csv``), and every derived number is computed here.

THE CHAIN OF REASONING, AND WHERE EACH ASSUMPTION ENTERS
--------------------------------------------------------
*Limit.*  For one frame the robust per-pixel noise sigma_pix [ADU] is
measured in a 512 x 512 px box at the candidate's position.  It is measured
twice: on the single frame (which includes fixed-pattern noise from hot
pixels and un-flattened pixel response) and on the DIFFERENCE of two frames
of the same field and night divided by sqrt(2) (which cancels everything
fixed in pixel space and is therefore the random noise that averages down in
a stack).  The single-frame value is reported beside it so the size of the
fixed-pattern term is visible.  Neither is trusted to scale as sqrt(N_pix):
the noise of an aperture-sized sum is measured DIRECTLY as the scatter of
33 x 33 px block sums in the pair difference (``block_sigma``), and the
ratio of that to sigma_pix * 33 — the "correlation factor", 1 for white
noise — multiplies the per-pixel noise of every frame of the field before
the limit is formed.  The statistical 3-sigma limit of an N-frame mean stack in
a circular aperture of radius r is

    F_3sig = 3 * sqrt(sum_i sigma_i^2) / N * sqrt(pi r^2 / s^2) / t * F_1

with s the pixel scale, t the exposure and F_1 the line flux that produces
1 ADU/s.  This is the photon/read-noise floor BEFORE continuum subtraction,
flat-field and sky-gradient systematics; the true limit is worse.

*F_1* needs the filter.  ``ZMAG`` in the frame headers is the magnitude of a
star giving 1 ADU/s (verified below: ~18.3 in H, ~21.8 in R for 512 s and
256 s frames alike).  Treating it as an AB magnitude at 6563 A and the
filter as a top-hat of width W gives F_1 = f_lambda(ZMAG) * W.  **W is not
known** — no transmission curve is on file.  The physicist's memo adopts
W ~ 65 A "from zero-point ratios"; this script uses the same value so that
his numbers can be checked like for like, emits the measured ZP(R)-ZP(H)
that the 65 A rests on, and every limit scales linearly with W.

*In band or not.*  A ~65 A filter centred on rest H-alpha passes the line
only for cz < ~1500 km/s.  Most of these candidates now HAVE published
velocities, so for each one the script states whether H-alpha falls in the
assumed band.  The filter CENTRE is as unknown as its width; the verdicts
flag which candidates sit within one half-width of the band edge.

*Expected flux.*  Where a GALEX FUV magnitude exists the expected H-alpha
flux follows from equating the two SFR calibrations the Local Volume
database itself uses (UNGC eq. 17-18 for FUV; Kaisin & Karachentsev 2013
for H-alpha):

    log F_Ha[erg/s/cm2] = -6.20 - 0.4 * (m_FUV - 1.93 A_B)

which is independent of distance.  Lee et al. (2009) show H-alpha
under-predicts FUV-based SFRs in dwarfs with SFR < ~0.01 Msun/yr, so this is
an UPPER expectation; the verdict therefore requires the prediction to
exceed the limit by a stated margin (``MARGIN``), not merely to touch it.
Where no FUV magnitude exists, a second predictor is used: a galaxy that
formed its stars at a constant rate (P = log(SFR*T0/M*) = 0), from the
published stellar mass and distance.  For the FUV-detected objects both are
computed, and their difference is emitted so the reader can see how good
the fallback is on this very sample.

VERDICT RULES (fixed here, before the table was looked at)
----------------------------------------------------------
``classified``  — a published radial velocity (hence distance and
                  membership) exists.  This is the DW-N1 column.
``ha_verdict``  — one of
    no_ha_data          no H frames in the archive
    out_of_band         published cz puts H-alpha outside the assumed band
    informative         in band (or velocity unknown), no published H-alpha
                        flux, and predicted flux >= MARGIN x our 3-sigma limit
    marginal            as above but limit <= predicted < MARGIN x limit
    below_depth         predicted flux < our 3-sigma limit: a non-detection
                        would say nothing (standing rule 2 of the synthesis)
    no_prediction       in band / unknown velocity but no FUV and no stellar
                        mass to predict from

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    cd "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO"
    $PY DwarfGalaxy_AGN_Survey/notes/novelty/build_dwarf_novelty.py census
    $PY DwarfGalaxy_AGN_Survey/notes/novelty/build_dwarf_novelty.py fetch     # network; refreshes sources/
    $PY DwarfGalaxy_AGN_Survey/notes/novelty/build_dwarf_novelty.py measure   # reads the archive, resumable
    $PY DwarfGalaxy_AGN_Survey/notes/novelty/build_dwarf_novelty.py table     # joins everything, writes out/
    $PY DwarfGalaxy_AGN_Survey/notes/novelty/build_dwarf_novelty.py figure
    $PY DwarfGalaxy_AGN_Survey/notes/novelty/build_dwarf_novelty.py all       # census+measure+table+figure

SAFETY
------
The manifest is opened with ``mode=ro``; the archive is only ever read; the
only files written are under this directory's ``out/`` and ``sources/``.
"""
from __future__ import annotations

import argparse
import csv
import html
import math
import re
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MANIFEST = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
OUT = HERE / "out"
SRC = HERE / "sources"
LIT = HERE / "literature_status.csv"

sys.path.insert(0, str(REPO / "pipeline"))

# ---------------------------------------------------------------- constants
C_KMS = 299792.458
HA_REST_A = 6562.8
#: Assumed H filter: top-hat, centred on rest H-alpha.  NOT MEASURED — both
#: numbers await Cannon's transmission curve.  65 A is the physicist's value.
W_ASSUMED_A = 65.0
BAND_CENTRE_A = HA_REST_A
PIXSCALE_ARCSEC = 0.539      #: S1 astrometry of Dw1446+58 (s1_batch), CDELT = 0.00015 deg
APERTURE_ARCSEC = 10.0       #: the physicist's aperture, kept for like-for-like
BOX_HALF = 256               #: half-size of the noise box [px]
BLOCK = 33                   #: block side [px]; 33^2 = 1089 px ~ the r = 10" aperture (1081 px)
MARGIN = 3.0                 #: predicted/limit ratio required for "informative"
T0_YR = 13.7e9               #: cosmic time used in the P parameter (Kaisin & Karachentsev 2013)
#: log SFR = log F_Ha + 2 log D + 8.98   (Kaisin & Karachentsev 2013, their SFR relation)
KK_HA_CONST = 8.98
#: log SFR = 2.78 - 0.4 m_FUV^c + 2 log D ; m^c = m - 1.93 A_B  (UNGC eq. 17-18)
UNGC_FUV_CONST = 2.78
UNGC_FUV_AB = 1.93
#: Faintest upper limit quoted in the BTA H-alpha table we read (KK2013: log F < -15.25)
BTA_LOGF_LIMIT = -15.25
#: Kennicutt & Evans (2012) L(Ha)->SFR, used ONLY to reproduce the physicist's SFR number
KE12_SFR_PER_LHA = 5.37e-42
MPC_CM = 3.0857e24
#: LVGDB galaxies WITH a published BTA H-alpha flux, used (a) to validate the
#: FUV->H-alpha predictor on real numbers and (b) as possible on-chip flux
#: standards: UGC 9660 lies 16 arcmin from Dw1459+44, UGC 9992 55 arcmin from
#: Dw1533+67 (Nazarova et al. 2025 notes); NGC 5238 is the project's own target.
KNOWN_HA = {"370": ("NGC 5238", "NGC5238"), "837": ("UGC 9660", "Dw1459+44"), "838": ("UGC 9992", "Dw1533+67")}

FIELD_SQL = """
SELECT canonical_target, filter, COUNT(*) AS n_frames, COUNT(DISTINCT night) AS n_nights,
       SUM(exptime) AS exp_s, MIN(night) AS first_night, MAX(night) AS last_night,
       GROUP_CONCAT(DISTINCT readoutm) AS readout
FROM frames
WHERE is_canonical = 1 AND imagetyp LIKE 'Light%' AND error IS NULL
  AND (canonical_target LIKE 'Dw%' OR canonical_target IN ('NGC5238', 'NGC 5548'))
GROUP BY canonical_target, filter
ORDER BY canonical_target, filter
"""


def manifest():
    """Open the shared manifest strictly read-only (ground rule)."""
    con = sqlite3.connect(f"file:{MANIFEST}?mode=ro", uri=True, timeout=300)
    con.row_factory = sqlite3.Row
    return con


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    """Write ``rows`` with a stable column order; every output goes through here."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(x):
    """CSV cell -> float or None ('' and None are missing, never zero)."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------- census
def cmd_census(_args=None) -> list[dict]:
    """Filters, nights and total exposure per field, straight from ``frames``.

    Uses the manifest's own canonical-frame flag (``is_canonical``), so the
    counts are deduplicated exactly as S0 deduplicates them.  NGC 5548 rows
    are listed by header filter string; whether slot '6' frames are images or
    spectra is ``frame_dispersion``'s business and is summarised separately.
    """
    con = manifest()
    rows = [dict(r) for r in con.execute(FIELD_SQL)]
    for r in rows:
        r["exp_h"] = round(r["exp_s"] / 3600.0, 3)
        r["exp_s"] = round(r["exp_s"], 1)
    write_csv(OUT / "field_census.csv", rows)

    # NGC 5548 slot-'6' nights with the S2c pixel verdicts, to place the
    # series on the calendar of the published 2023 campaign.
    q = """SELECT f.night, MIN(f.jd) AS jd_first, MAX(f.jd) AS jd_last, COUNT(*) AS n,
                  SUM(d.verdict = 'dispersed') AS n_dispersed,
                  SUM(d.verdict = 'direct') AS n_direct,
                  SUM(d.verdict NOT IN ('dispersed', 'direct') OR d.verdict IS NULL) AS n_other
           FROM frames f LEFT JOIN frame_dispersion d ON d.obs_rowid = f.obs_rowid
           WHERE f.is_canonical = 1 AND f.canonical_target = 'NGC 5548' AND f.filter = '6'
           GROUP BY f.night ORDER BY f.night"""
    n5548 = [dict(r) for r in con.execute(q)]
    write_csv(OUT / "ngc5548_slot6_nights.csv", n5548)
    con.close()
    print(f"census: {len(rows)} field x filter rows; NGC 5548 slot-6 nights: {len(n5548)}")
    return rows


# -------------------------------------------------------------------- fetch
def _get(url: str, timeout: int = 60) -> bytes:
    return urllib.request.urlopen(url, timeout=timeout).read()


def parse_lvgdb(page: str) -> dict:
    """Pull the handful of fields we use out of one LVGDB object page.

    The page is a fixed-layout table.  After the header cell ``Method`` the
    first seven value cells are, in order, a_Holm, b/a, A_B, m_FUV, B_T,
    m_Ha, K_s (an em-dash means "no value").  Further down, after the cell
    ``lg[SFR],``, come three rows labelled ``UV``, ``Hα`` and ``HI``; the
    H-alpha row carries a flux in erg cm^-2 s^-1 only when the database
    holds an H-alpha measurement.  Its ABSENCE is the fact we record, so the
    parser is tested on pages with (NGC 5238) and without (Dw1735+57) one.
    """
    page = re.sub(r"<style.*?</style>", "", page, flags=re.S)
    text = html.unescape(re.sub(r"<[^>]+>", "|", page))
    flat = [t.strip() for t in text.split("|") if t.strip()]
    out: dict = {}
    pos = re.search(r"(\d\d \d\d \d\d\.\d) \| ([+-]\d\d \d\d \d\d)", " | ".join(flat))
    if pos:
        out["ra_hms"], out["dec_dms"] = pos.group(1).replace(" ", ":"), pos.group(2).replace(" ", ":")
    i = flat.index("Method")
    for name, val in zip(("a_holm", "b_over_a", "A_B", "m_FUV", "B_T", "m_Ha", "K_s"), flat[i + 1:i + 8]):
        out[name] = None if val in ("—", "-") else val
    joined = " | ".join(flat[i:])
    m = re.search(r"\| (-?\d+) \| ([\d.]+) \| \2 \| (?:--> \| )?(\w+) \| A \| 26", joined)
    if m:
        out["V_h"], out["D_Mpc"], out["D_method"] = m.group(1), m.group(2), m.group(3)
    k = flat.index("lg[SFR],")
    iu = flat.index("UV", k)
    ih = flat.index("Hα", iu)
    ii = flat.index("HI", ih)
    sfr = [t for t in flat[iu + 1:ih] if re.fullmatch(r"-?\d+\.\d+", t)]
    out["logSFR_FUV"] = sfr[-1] if sfr else None
    ha = re.search(r"\(([\d.]+)±([\d.]+)\)x10 \| (-?\d+) \| \[erg", " | ".join(flat[ih + 1:ii]))
    out["F_Ha"] = float(ha.group(1)) * 10 ** int(ha.group(3)) if ha else None
    out["e_F_Ha"] = float(ha.group(2)) * 10 ** int(ha.group(3)) if ha else None
    return out


def cmd_fetch(_args=None) -> None:
    """Refresh the two machine-read sources (network).  Safe to skip: the
    cached copies under ``sources/`` are what ``table`` reads."""
    lit = read_csv(LIT)
    ids = sorted({r["lvgdb_id"] for r in lit if r["lvgdb_id"]} | {"370", "837", "838"})
    (SRC / "lvgdb").mkdir(parents=True, exist_ok=True)
    for i in ids:
        url = f"https://www.sao.ru/lv/lvgdb/object.php?name=x&id={i}"        # the id selects the object
        try:
            page = _get(url).decode("koi8-r", "replace")
        except Exception as exc:                                  # keep the cached copy
            print(f"  LVGDB id {i}: {exc} (cached copy kept)")
            continue
        (SRC / "lvgdb" / f"lvgdb_{i}.html").write_text(page)
        time.sleep(0.5)
    rows = []
    for r in lit:
        ra, dec = sexa_to_deg(r["ra_hms"], r["dec_dms"])
        q = urllib.parse.urlencode({
            "-source": "II/335/galex_ais", "-c": f"{ra:.5f} {dec:+.5f}", "-c.rs": "30",
            "-out": "_r,RAJ2000,DEJ2000,FUVmag,e_FUVmag,NUVmag,e_NUVmag", "-sort": "_r"})
        try:
            txt = _get("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?" + q).decode()
        except Exception as exc:
            print(f"  GALEX {r['candidate']}: {exc}")
            return                                                # do not overwrite the cache with a partial
        for line in txt.splitlines():
            p = line.split("\t")
            if len(p) == 7 and re.match(r"\s*\d", p[0]):
                rows.append(dict(candidate=r["candidate"], sep_arcsec=p[0].strip(),
                                 ra=p[1].strip(), dec=p[2].strip(),
                                 FUVmag=p[3].strip(), e_FUVmag=p[4].strip(),
                                 NUVmag=p[5].strip(), e_NUVmag=p[6].strip()))
        time.sleep(0.4)
    write_csv(SRC / "galex_ais_cone30.csv", rows)
    (SRC / "FETCHED.txt").write_text(time.strftime("fetched %Y-%m-%dT%H:%M:%SZ\n", time.gmtime()))
    print(f"fetch: {len(ids)} LVGDB pages, {len(rows)} GALEX AIS rows")


def sexa_to_deg(ra_hms: str, dec_dms: str) -> tuple[float, float]:
    h, m, s = (float(x) for x in ra_hms.split(":"))
    sign = -1.0 if dec_dms.strip().startswith("-") else 1.0
    d, dm, ds = (abs(float(x)) for x in dec_dms.split(":"))
    return 15.0 * (h + m / 60 + s / 3600), sign * (d + dm / 60 + ds / 3600)


# ------------------------------------------------------------------ measure
def robust_sigma(a: np.ndarray, iters: int = 4) -> tuple[float, float]:
    """Median and noise sigma of a sky box, robust to stars and hot pixels.

    Two steps, because the raw frames are INTEGER ADU with sigma of only a
    few ADU: the MAD of such data is itself an integer, so 1.4826*MAD moves
    in steps of ~1.5 ADU (a 20 % quantisation of the answer — found on the
    first run of this script, where every frame returned 7.413 or 8.896).
    The MAD is therefore used only to place the clip; the returned sigma is
    the standard deviation of the pixels inside +-3 sigma_MAD, divided by
    0.98658 (the std of a unit Gaussian truncated at 3 sigma), which is
    continuous.  Returns (median, sigma).
    """
    a = a[np.isfinite(a)].astype(float).ravel()
    med = np.median(a)
    sig = max(1.4826 * np.median(np.abs(a - med)), 0.5)
    keep = np.abs(a - med) <= 3 * sig
    for _ in range(iters):
        if keep.sum() < 100:
            break
        sig = a[keep].std(ddof=1) / 0.98658
        med = np.median(a[keep])
        keep = np.abs(a - med) <= 3 * sig
    return float(med), float(sig)


def block_sigma(diff: np.ndarray, block: int = BLOCK) -> float:
    """Empirical noise of an aperture-sized SUM, from a pair-difference image.

    sigma_pix * sqrt(N) is the aperture noise only if pixels are independent
    and the background is flat.  Neither is guaranteed here (the camera's
    readout modes combine reads; sky gradients and scattered light differ
    between dithered frames), and an optimistic limit is the bias this
    package exists to prevent.  So the box is cut into ``block`` x ``block``
    cells, each cell is summed after removing the box median, and the robust
    scatter of the cell sums is returned, per single frame (hence /sqrt 2).
    Cells holding a star residual are outliers and the MAD ignores them (the
    sums are not integer-quantised at this level, so the MAD is safe here).
    """
    n = (diff.shape[0] // block) * block
    d = diff[:n, :n] - np.median(diff)
    sums = d.reshape(n // block, block, n // block, block).sum(axis=(1, 3)).ravel()
    return float(1.4826 * np.median(np.abs(sums - np.median(sums))) / math.sqrt(2.0))


def load_frame(relpath: str):
    """Read one archive frame (fpack-compressed, READ-ONLY).  Returns
    (data, header) of the first image HDU that has pixels."""
    from astropy.io import fits
    with fits.open(ARCHIVE / relpath, mode="readonly", memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) == 2:
                return np.asarray(hdu.data, dtype=np.float32), hdu.header.copy()
    raise ValueError(f"no image HDU in {relpath}")


def candidate_pixel(header, ra: float, dec: float):
    """Candidate position on the chip from the frame's own WCS, or None.

    Only frames the observatory software plate-solved carry a WCS.  The
    answer matters twice: the noise box is placed there, and a candidate
    that is not on the chip cannot be measured at all.
    """
    import warnings
    from astropy.wcs import WCS
    if "CRVAL1" not in header or "CRPIX1" not in header:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            x, y = WCS(header).all_world2pix([[ra, dec]], 0)[0]
        except Exception:
            return None
    return float(x), float(y)


def cmd_measure(args=None) -> None:
    """Measure sky noise in every H frame of every Dw field (resumable).

    One row per frame in ``out/halpha_frame_noise.csv``; frames already in
    the file are skipped, so a killed run loses only the frame in flight.
    """
    lit = {r["field"]: r for r in read_csv(LIT)}
    con = manifest()
    frames = [dict(r) for r in con.execute(
        """SELECT canonical_target AS field, path, night, jd, exptime, zmag, fwhm, moonphas
           FROM frames WHERE is_canonical = 1 AND imagetyp LIKE 'Light%' AND error IS NULL
             AND canonical_target LIKE 'Dw%' AND filter = 'H'
           ORDER BY canonical_target, jd""")]
    # a WCS-bearing frame per field (any filter) gives the candidate pixel for
    # H frames the observatory did not solve; pointing repeats to ~1 arcmin.
    wcs_donor = {r["field"]: r["path"] for r in con.execute(
        """SELECT canonical_target AS field, MIN(path) AS path FROM frames
           WHERE is_canonical = 1 AND canonical_target LIKE 'Dw%' AND crval1 IS NOT NULL
           GROUP BY canonical_target""")}
    con.close()

    out_path = OUT / "halpha_frame_noise.csv"
    done = {r["path"]: r for r in read_csv(out_path)} if out_path.exists() else {}
    rows = list(done.values())
    fields = ["field", "path", "night", "jd", "exptime", "zmag", "fwhm", "moonphas",
              "x_cand", "y_cand", "pos_basis", "on_chip", "box_x", "box_y", "sky_adu", "sigma_single_adu",
              "pair_with", "sigma_pairdiff_adu", "sigma_block_adu"]
    donor_xy: dict[str, tuple | None] = {}
    field_box: dict[str, tuple[int, int]] = {r["field"]: (int(r["box_x"]), int(r["box_y"])) for r in rows}
    prev = None          # (field, night, path, box_array, (x, y))
    t0 = time.time()
    for fr in frames:
        if fr["path"] in done:
            prev = None
            continue
        ra, dec = sexa_to_deg(lit[fr["field"]]["ra_hms"], lit[fr["field"]]["dec_dms"])
        data, hdr = load_frame(fr["path"])
        ny, nx = data.shape
        xy, basis = candidate_pixel(hdr, ra, dec), "frame WCS"
        if xy is None:
            if fr["field"] not in donor_xy:
                d = wcs_donor.get(fr["field"])
                donor_xy[fr["field"]] = candidate_pixel(load_frame(d)[1], ra, dec) if d else None
            xy, basis = donor_xy[fr["field"]], "WCS of another frame of the field"
        if xy is None:
            xy, basis = (nx / 2, ny / 2), "chip centre (no WCS in any frame of this field)"
        on_chip = 0 <= xy[0] < nx and 0 <= xy[1] < ny
        # One box per FIELD, fixed in pixel space (set by the first frame):
        # the pair difference below needs the identical pixels in both frames,
        # and the candidate moves by only ~1 px between frames of a field.
        if fr["field"] not in field_box:
            field_box[fr["field"]] = (int(min(max(xy[0], BOX_HALF), nx - BOX_HALF)),
                                      int(min(max(xy[1], BOX_HALF), ny - BOX_HALF)))
        cx, cy = field_box[fr["field"]]
        box = data[cy - BOX_HALF:cy + BOX_HALF, cx - BOX_HALF:cx + BOX_HALF]
        sky, s1 = robust_sigma(box)
        row = dict(fr, x_cand=round(xy[0], 1), y_cand=round(xy[1], 1), pos_basis=basis,
                   on_chip=int(on_chip), box_x=cx, box_y=cy, sky_adu=round(sky, 2), sigma_single_adu=round(s1, 3),
                   pair_with="", sigma_pairdiff_adu="", sigma_block_adu="")
        # pair difference with the previous frame of the same field and night,
        # same box in PIXEL space (dither moves the stars, not the hot pixels)
        if prev and prev[0] == fr["field"] and prev[1] == fr["night"] and prev[3].shape == box.shape \
                and prev[4] == (cx, cy):
            _, sd = robust_sigma(box - prev[3])
            row["pair_with"] = prev[2]
            row["sigma_pairdiff_adu"] = round(sd / math.sqrt(2.0), 3)
            row["sigma_block_adu"] = round(block_sigma(box - prev[3]), 2)
        prev = (fr["field"], fr["night"], fr["path"], box.copy(), (cx, cy))
        rows.append(row)
        write_csv(out_path, rows, fields)
        if args and args.max_seconds and time.time() - t0 > args.max_seconds:
            print(f"measure: stopped cleanly after {len(rows)}/{len(frames)} frames")
            return
    print(f"measure: {len(rows)}/{len(frames)} H frames measured")

    # Is a galaxy with a PUBLISHED H-alpha flux on the same chip?  If so the
    # field carries its own end-to-end flux check (filter width included).
    nb = []
    for lv_id, (name, field) in KNOWN_HA.items():
        page = SRC / "lvgdb" / f"lvgdb_{lv_id}.html"
        donor = wcs_donor.get(field)
        if not (page.exists() and donor):
            continue
        lv = parse_lvgdb(page.read_text())
        ra, dec = sexa_to_deg(lv["ra_hms"], lv["dec_dms"])
        data, hdr = load_frame(donor)
        xy = candidate_pixel(hdr, ra, dec)
        ny, nx = data.shape
        nb.append(dict(galaxy=name, field=field, donor_frame=donor, x=round(xy[0], 1), y=round(xy[1], 1),
                       on_chip=int(0 <= xy[0] < nx and 0 <= xy[1] < ny),
                       edge_margin_px=round(min(xy[0], xy[1], nx - xy[0], ny - xy[1]), 0),
                       F_Ha_published=lv["F_Ha"], e_F_Ha=lv["e_F_Ha"]))
    if nb:
        write_csv(OUT / "known_halpha_neighbours.csv", nb)


# -------------------------------------------------------------------- table
def flam_ab(mag: float, lam_a: float = HA_REST_A) -> float:
    """f_lambda [erg/s/cm2/A] of an AB magnitude at wavelength lam_a."""
    fnu = 3.631e-20 * 10 ** (-0.4 * mag)
    return fnu * (C_KMS * 1e13) / lam_a ** 2


def limit_3sigma(sigmas: list[float], exptime: float, zp: float, r_arcsec: float,
                 width_a: float = W_ASSUMED_A) -> float:
    """Statistical 3-sigma line-flux limit of a mean stack (see module doc)."""
    n = len(sigmas)
    npix = math.pi * (r_arcsec / PIXSCALE_ARCSEC) ** 2
    sig_rate = math.sqrt(sum(s * s for s in sigmas)) / n * math.sqrt(npix) / exptime
    return 3.0 * sig_rate * flam_ab(zp) * width_a


def predicted_flux_fuv(m_fuv: float, a_b: float) -> float:
    """Expected OBSERVED H-alpha flux if SFR(Ha) = SFR(FUV)  (distance-free)."""
    return 10 ** (UNGC_FUV_CONST - KK_HA_CONST - 0.4 * (m_fuv - UNGC_FUV_AB * a_b))


def predicted_flux_p0(log_mstar: float, d_mpc: float) -> float:
    """Expected H-alpha flux for constant star formation over T0 (P = 0)."""
    return 10 ** (log_mstar - math.log10(T0_YR) - KK_HA_CONST - 2 * math.log10(d_mpc))


def band_status(v_kms: float | None) -> tuple[str, float | None]:
    """Where redshifted H-alpha falls relative to the ASSUMED top-hat band."""
    if v_kms is None:
        return "velocity unknown", None
    dlam = HA_REST_A * v_kms / C_KMS + (HA_REST_A - BAND_CENTRE_A)
    half = W_ASSUMED_A / 2
    if abs(dlam) <= half:
        return ("in band (edge)" if abs(dlam) > half / 2 else "in band"), dlam
    return ("out of band (within one width)" if abs(dlam) <= 2 * half else "out of band"), dlam


def ha_verdict(n_h: int, band: str, f_lit, f_pred, f_lim) -> str:
    """The pre-declared verdict rules of the module docstring, as code."""
    if n_h == 0:
        return "no_ha_data"
    if band.startswith("out of band"):
        return "out_of_band"
    if f_lit is not None:
        return "already_measured"
    if f_pred is None or f_lim is None:
        return "no_prediction"
    if f_pred >= MARGIN * f_lim:
        return "informative"
    if f_pred >= f_lim:
        return "marginal"
    return "below_depth"


def cmd_table(_args=None) -> list[dict]:
    lit = read_csv(LIT)
    census = read_csv(OUT / "field_census.csv")
    noise = read_csv(OUT / "halpha_frame_noise.csv") if (OUT / "halpha_frame_noise.csv").exists() else []
    galex = read_csv(SRC / "galex_ais_cone30.csv") if (SRC / "galex_ais_cone30.csv").exists() else []
    by_field: dict[str, dict] = {}
    for c in census:
        by_field.setdefault(c["canonical_target"], {})[c["filter"]] = c

    # global H zero point and the R-H difference the 65 A assumption rests on
    con = manifest()
    zp = {f: [r[0] for r in con.execute(
        "SELECT zmag FROM frames WHERE is_canonical=1 AND canonical_target LIKE 'Dw%' "
        "AND filter=? AND zmag IS NOT NULL", (f,))] for f in ("H", "R", "L")}
    con.close()
    zp_h_global = float(np.median(zp["H"]))
    zp_r_global = float(np.median(zp["R"]))

    cf_all = [fnum(n["sigma_block_adu"]) / (BLOCK * fnum(n["sigma_pairdiff_adu"]))
              for n in noise if fnum(n.get("sigma_block_adu"))]
    corr_global = float(np.median(cf_all)) if cf_all else 1.0

    rows = []
    for r in lit:
        f = r["field"]
        cen = by_field.get(f, {})
        row = dict(field=f, candidate=r["candidate"], ra_hms=r["ra_hms"], dec_dms=r["dec_dms"])
        for filt in ("L", "R", "H"):
            c = cen.get(filt)
            row[f"n_{filt}"] = int(c["n_frames"]) if c else 0
            row[f"nights_{filt}"] = int(c["n_nights"]) if c else 0
            row[f"exp_h_{filt}"] = float(c["exp_h"]) if c else 0.0
        row.update(first_listing=r["first_listing"], morph_lit=r["morph_lit"],
                   association=r["association"], v_hel_kms=r["v_hel_kms"], v_method=r["v_method"],
                   v_ref=r["v_ref"], D_Mpc=r["D_Mpc"], D_method=r["D_method"], B_mag=r["B_mag"],
                   spec_emission=r["spec_emission"], spec_ref=r["spec_ref"],
                   HI=("detected S=%s Jy km/s" % r["S_HI_Jykms"]) if r["S_HI_Jykms"] else
                      ("not detected / confused" if r["first_listing"] else ""))

        # ---- LVGDB snapshot (machine-read): A_B, m_FUV, published H-alpha flux
        lv = {}
        if r["lvgdb_id"]:
            p = SRC / "lvgdb" / f"lvgdb_{r['lvgdb_id']}.html"
            if p.exists():
                lv = parse_lvgdb(p.read_text())
        a_b = fnum(lv.get("A_B")) or 0.0
        row["A_B"] = a_b if lv else ""
        row["in_LVGDB"] = int(bool(lv))
        f_lit = lv.get("F_Ha") if lv else None
        row["F_Ha_published"] = f_lit if f_lit is not None else ""
        row["F_Ha_published_src"] = ("LVGDB" if f_lit is not None else
                                     ("none in LVGDB" if lv else "not in LVGDB; none found in refs"))

        # ---- absolute magnitude at the published distance
        d = fnum(r["D_Mpc"]); b = fnum(r["B_mag"]); v = fnum(r["v_hel_kms"])
        row["M_B"] = round(b - 5 * math.log10(d) - 25 - a_b, 2) if (d and b) else ""
        row["classified"] = int(v is not None)

        # ---- FUV magnitude: published table > LVGDB > GALEX AIS nearest <= 10"
        m_fuv, src = fnum(r["m_fuv_pub"]), r["m_fuv_ref"]
        if m_fuv is None and lv.get("m_FUV") and not str(lv["m_FUV"]).startswith(">"):
            m_fuv, src = fnum(lv["m_FUV"]), "LVGDB"
        if m_fuv is None:
            near = [g for g in galex if g["candidate"] == r["candidate"] and g["FUVmag"]
                    and float(g["sep_arcsec"]) <= 10.0]
            if near:
                m_fuv, src = float(near[0]["FUVmag"]), "GALEX AIS (nearest <= 10 arcsec)"
        row["m_FUV"] = m_fuv if m_fuv is not None else ""
        row["m_FUV_src"] = src if m_fuv is not None else ""
        f_fuv = predicted_flux_fuv(m_fuv, a_b) if m_fuv is not None else None
        ms = fnum(r["logMstar"])
        f_p0 = predicted_flux_p0(ms, d) if (ms and d) else None
        row["F_pred_FUV"] = f"{f_fuv:.2e}" if f_fuv else ""
        row["F_pred_P0"] = f"{f_p0:.2e}" if f_p0 else ""
        row["P_FUV"] = round(math.log10(f_fuv / f_p0), 2) if (f_fuv and f_p0) else ""
        f_pred = f_fuv if f_fuv else f_p0
        row["F_pred_basis"] = "FUV" if f_fuv else ("P=0" if f_p0 else "")
        row["L_pred_erg_s"] = f"{f_pred * 4 * math.pi * (d * MPC_CM) ** 2:.2e}" if (f_pred and d) else ""

        # ---- measured limit
        fr = [n for n in noise if n["field"] == f]
        f_lim = None
        if fr:
            sig0 = [fnum(n["sigma_pairdiff_adu"]) or fnum(n["sigma_single_adu"]) for n in fr]
            sig1 = [fnum(n["sigma_single_adu"]) for n in fr]
            # correlation factor: measured block-sum noise / white-noise expectation
            cf = [fnum(n["sigma_block_adu"]) / (BLOCK * fnum(n["sigma_pairdiff_adu"]))
                  for n in fr if fnum(n.get("sigma_block_adu"))]
            corr = float(np.median(cf)) if cf else corr_global
            sig = [s_ * corr for s_ in sig0]
            zps = [fnum(n["zmag"]) for n in fr if fnum(n["zmag"])]
            zp_f = float(np.median(zps)) if zps else zp_h_global
            t = float(np.median([float(n["exptime"]) for n in fr]))
            f_lim = limit_3sigma(sig, t, zp_f, APERTURE_ARCSEC)
            row.update(n_H_measured=len(fr), zp_H=round(zp_f, 2),
                       zp_H_basis="field median" if zps else "survey median (no ZMAG in this field)",
                       sigma_pix_adu=round(float(np.median(sig)), 2),
                       fpn_ratio=round(float(np.median(sig1)) / float(np.median(sig0)), 2),
                       corr_factor=round(corr, 2), n_pairs=len(cf),
                       F_lim_3sig_white=f"{limit_3sigma(sig0, t, zp_f, APERTURE_ARCSEC):.2e}",
                       cand_on_chip=min(int(n["on_chip"]) for n in fr),
                       F_lim_3sig=f"{f_lim:.2e}",
                       L_lim_erg_s=f"{f_lim * 4 * math.pi * (d * MPC_CM) ** 2:.2e}" if d else "",
                       logSFR_lim=round(math.log10(f_lim) + 2 * math.log10(d) + KK_HA_CONST, 2) if d else "",
                       lim_over_BTA=round(f_lim / 10 ** BTA_LOGF_LIMIT, 1))
        band, dlam = band_status(v)
        row["band"] = band
        row["dlam_A"] = round(dlam, 1) if dlam is not None else ""
        row["pred_over_lim"] = round(f_pred / f_lim, 2) if (f_pred and f_lim) else ""
        row["ha_verdict"] = ha_verdict(row["n_H"], band, f_lit, f_pred, f_lim)
        rows.append(row)

    cols = ["field", "candidate", "ra_hms", "dec_dms", "n_L", "nights_L", "exp_h_L", "n_R", "nights_R",
            "exp_h_R", "n_H", "nights_H", "exp_h_H", "first_listing", "morph_lit", "association",
            "v_hel_kms", "v_method", "v_ref", "D_Mpc", "D_method", "B_mag", "A_B", "M_B", "HI",
            "spec_emission", "spec_ref", "in_LVGDB", "F_Ha_published", "F_Ha_published_src",
            "classified", "m_FUV", "m_FUV_src", "F_pred_FUV", "F_pred_P0", "P_FUV", "F_pred_basis",
            "L_pred_erg_s", "n_H_measured", "zp_H", "zp_H_basis", "sigma_pix_adu", "fpn_ratio",
            "corr_factor", "n_pairs", "cand_on_chip", "F_lim_3sig", "F_lim_3sig_white", "L_lim_erg_s", "logSFR_lim",
            "lim_over_BTA", "band", "dlam_A", "pred_over_lim", "ha_verdict"]
    write_csv(OUT / "novelty_table.csv", rows, cols)

    # ---------------- the physicist's arithmetic, re-derived number by number
    lims = [float(r["F_lim_3sig"]) for r in rows if r.get("F_lim_3sig")]
    sig_all = [(fnum(n["sigma_pairdiff_adu"]) or fnum(n["sigma_single_adu"])) * corr_global for n in noise]
    sky = [fnum(n["sky_adu"]) for n in noise]
    vm = [fnum(n["sigma_pairdiff_adu"]) ** 2 / fnum(n["sky_adu"]) for n in noise if fnum(n["sigma_pairdiff_adu"])]
    f1 = flam_ab(18.3) * W_ASSUMED_A
    seven = limit_3sigma([float(np.median(sig_all))] * 7, 512.0, zp_h_global, APERTURE_ARCSEC) if sig_all else None
    l10 = seven * 4 * math.pi * (10 * MPC_CM) ** 2 if seven else None
    dz = zp_r_global - zp_h_global
    check = [
        dict(quantity="ZMAG(H) survey median [mag for 1 ADU/s]", memo="18.3", here=round(zp_h_global, 2),
             note=f"{len(zp['H'])} H frames with ZMAG"),
        dict(quantity="F_1 = flux for 1 ADU/s at ZP 18.3, W 65 A [erg/s/cm2]", memo="8e-15", here=f"{f1:.2e}",
             note="AB f_lambda at 6563 A x 65 A"),
        dict(quantity="3-sigma limit, 7 x 512 s, r = 10 arcsec [erg/s/cm2]", memo="~1e-14",
             here=f"{seven:.2e}" if seven else "",
             note=f"median measured noise incl. correlation factor {corr_global:.2f}; statistical only"),
        dict(quantity="same, as L(Ha) at 10 Mpc [erg/s]", memo="~1e38", here=f"{l10:.2e}" if l10 else "", note=""),
        dict(quantity="same, as SFR at 10 Mpc, Kennicutt & Evans 2012 [Msun/yr]", memo="~6e-4",
             here=f"{l10 * KE12_SFR_PER_LHA:.1e}" if l10 else "", note="the memo's calibration"),
        dict(quantity="same, as SFR at 10 Mpc, LVGDB relation [Msun/yr]", memo="",
             here=f"{seven * 100 * 10 ** KK_HA_CONST:.1e}" if seven else "", note="relation used in the table"),
        dict(quantity="our limit / BTA limit (log F < -15.25)", memo="10-100x shallower",
             here=f"{seven / 10 ** BTA_LOGF_LIMIT:.0f}x" if seven else "",
             note="BTA: 2x600 s on 6 m, faintest upper limit in KK2013 table"),
        dict(quantity="per-field measured limits, min-median-max [erg/s/cm2]", memo="",
             here=(f"{min(lims):.1e} - {float(np.median(lims)):.1e} - {max(lims):.1e}" if lims else ""),
             note="fields have 2-11 H frames, not 7"),
        dict(quantity="block-sum noise / white-noise expectation (33x33 px), median over pairs", memo="",
             here=round(corr_global, 2), note=f"{len(cf_all)} frame pairs; 1.00 = independent pixels, flat background"),
        dict(quantity="raw H sky: variance/mean [ADU], median over pairs", memo="",
             here=round(float(np.median(vm)), 2) if vm else "",
             note=f"median sky {float(np.median(sky)):.0f} ADU incl. any pedestal; header EGAIN 1.054 e/ADU would give >= 0.95 "
                  "for pure sky Poisson noise. For the detector package (D2), not resolved here"),
        dict(quantity="ZP(R) - ZP(H) [mag]", memo="", here=round(dz, 2),
             note=f"=> W_H = W_R x {10 ** (-0.4 * dz):.4f} for equal peak throughput: "
                  f"{1000 * 10 ** (-0.4 * dz):.0f} A if W_R = 1000 A, {1400 * 10 ** (-0.4 * dz):.0f} A if 1400 A"),
    ]
    write_csv(OUT / "physicist_depth_check.csv", check)

    # ---------------- does the FUV predictor reproduce PUBLISHED H-alpha fluxes?
    val = []
    for lv_id, (name, _field) in KNOWN_HA.items():
        page = SRC / "lvgdb" / f"lvgdb_{lv_id}.html"
        if not page.exists():
            continue
        lv = parse_lvgdb(page.read_text())
        pred = predicted_flux_fuv(float(lv["m_FUV"]), float(lv["A_B"]))
        val.append(dict(galaxy=name, m_FUV=lv["m_FUV"], A_B=lv["A_B"], F_Ha_predicted=f"{pred:.2e}",
                        F_Ha_published=f"{lv['F_Ha']:.2e}", log_pred_minus_obs=round(math.log10(pred / lv["F_Ha"]), 2)))
    write_csv(OUT / "predictor_validation.csv", val)

    # ---------------- gate summary
    n = len(rows)
    n_cls = sum(r["classified"] for r in rows)
    n_h = sum(1 for r in rows if r["n_H"] > 0)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["ha_verdict"]] = counts.get(r["ha_verdict"], 0) + 1
    lv_members = [r["candidate"] for r in rows if fnum(r["D_Mpc"]) and fnum(r["D_Mpc"]) < 12]
    summary = [dict(key="fields", value=n),
               dict(key="classified (published velocity)", value=n_cls),
               dict(key="unclassified", value="; ".join(r["candidate"] for r in rows if not r["classified"])),
               dict(key="within Local Volume (D < 12 Mpc)", value="; ".join(lv_members)),
               dict(key="fields with H frames", value=n_h),
               dict(key="published H-alpha flux for any candidate", value=sum(1 for r in rows if r["F_Ha_published"] != "")),
               ] + [dict(key=f"ha_verdict = {k}", value=f"{v}: " + "; ".join(
                   r["candidate"] for r in rows if r["ha_verdict"] == k)) for k, v in sorted(counts.items())]
    write_csv(OUT / "gate_summary.csv", summary)
    write_markdown(rows, check, summary)
    for s in summary:
        print(f"  {s['key']}: {s['value']}")
    return rows


def write_markdown(rows, check, summary) -> None:
    """Human-readable rendering of the three tables (the CSVs are canonical)."""
    def tab(hdr, body):
        return ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)] + \
               ["| " + " | ".join(str(x) for x in b) + " |" for b in body]
    L = ["# Dwarf fields — literature status and H-alpha depth (generated)",
         "", "Generated by `build_dwarf_novelty.py table`; do not edit. Fluxes in erg s^-1 cm^-2.", "",
         "## Census and classification", ""]
    L += tab(["field", "candidate", "L/R/H frames", "H exp [h]", "first listing", "cz [km/s]", "D [Mpc]",
              "M_B", "spectrum", "published F(Ha)"],
             [[r["field"], r["candidate"], f"{r['n_L']}/{r['n_R']}/{r['n_H']}", r["exp_h_H"],
               r["first_listing"], r["v_hel_kms"] or "—", r["D_Mpc"] or "—", r["M_B"] or "—",
               "emission" if r["spec_emission"] else "—", r["F_Ha_published"] or "none"] for r in rows])
    L += ["", "## Depth versus expectation (fields with H frames)", ""]
    L += tab(["candidate", "N_H", "3σ limit", "(white-noise)", "corr", "limit/BTA", "predicted F(Ha)", "basis",
              "pred/limit", "band (assumed 65 Å at 6563)", "verdict"],
             [[r["candidate"], r["n_H"], r.get("F_lim_3sig", ""), r.get("F_lim_3sig_white", ""),
               r.get("corr_factor", ""), r.get("lim_over_BTA", ""),
               r["F_pred_FUV"] or r["F_pred_P0"] or "—", r["F_pred_basis"] or "—", r["pred_over_lim"] or "—",
               r["band"], r["ha_verdict"]] for r in rows if r["n_H"] > 0])
    L += ["", "## The physicist's arithmetic, re-derived", ""]
    L += tab(["quantity", "memo", "here", "note"], [[c["quantity"], c["memo"], c["here"], c["note"]] for c in check])
    if (OUT / "predictor_validation.csv").exists():
        v = read_csv(OUT / "predictor_validation.csv")
        L += ["", "## FUV predictor against published BTA H-alpha fluxes (LVGDB)", ""]
        L += tab(list(v[0].keys()), [list(x.values()) for x in v])
    L += ["", "## Gate summary", ""]
    L += tab(["key", "value"], [[s["key"], s["value"]] for s in summary])
    (OUT / "NOVELTY_TABLE.md").write_text("\n".join(L) + "\n")


# ------------------------------------------------------------------- figure
def cmd_figure(_args=None) -> None:
    """One figure: predicted H-alpha flux against the measured limit, per field.

    Filled markers are predictions (circle = FUV-based, square = P=0
    fallback); the open downward triangle is OUR 3-sigma limit (a floor, in
    the house convention).  Candidates whose published velocity puts
    H-alpha outside the assumed band are drawn in the muted colour: for
    them the comparison is moot.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from macro_core import plotstyle as ps

    rows = [r for r in read_csv(OUT / "novelty_table.csv") if int(r["n_H"]) > 0]
    rows.sort(key=lambda r: -(fnum(r["F_pred_FUV"]) or fnum(r["F_pred_P0"]) or 0))
    x = np.arange(len(rows))
    with ps.context("print"):
        fig, ax = plt.subplots(figsize=(ps.COL_DOUBLE, 3.6))
        for i, r in enumerate(rows):
            out = r["band"].startswith("out of band")
            col = ps.MUTED if out else ps.ACCENT
            if fnum(r["F_pred_FUV"]):
                ax.plot(i, fnum(r["F_pred_FUV"]), **ps.measurement_kw(color=col, marker="o"))
            if fnum(r["F_pred_P0"]):
                ax.plot(i, fnum(r["F_pred_P0"]), **ps.measurement_kw(color=ps.tint(col), marker="s"))
            if fnum(r["F_lim_3sig"]):
                ax.plot(i, fnum(r["F_lim_3sig"]), **ps.floor_kw(color=ps.BAD))
        ax.axhline(10 ** BTA_LOGF_LIMIT, **ps.reference_kw(color=ps.GOOD, style=":"))
        ax.text(len(rows) - 0.5, 10 ** BTA_LOGF_LIMIT * 1.15, "BTA 6 m upper limits (2 x 600 s)",
                ha="right", va="bottom", fontsize=7, color=ps.GOOD)
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r['candidate']}\n{r['v_hel_kms'] or '?'}" for r in rows], rotation=60,
                           ha="right", fontsize=6.5)
        ax.set_xlabel("candidate, with published cz (km s$^{-1}$)")
        ax.set_ylabel(r"H$\alpha$ flux (erg s$^{-1}$ cm$^{-2}$)")
        handles = [ps.measurement_handle("predicted from GALEX FUV", color=ps.ACCENT, marker="o"),
                   ps.measurement_handle("predicted, constant SFR (P = 0)", color=ps.tint(ps.ACCENT), marker="s"),
                   ps.measurement_handle(r"H$\alpha$ outside assumed band", color=ps.MUTED, marker="o"),
                   ps.floor_handle(r"RLMT 3$\sigma$ limit, r = 10$''$ (measured noise)", color=ps.BAD)]
        ax.legend(handles=handles, fontsize=6.5, loc="upper right", frameon=False)
        fig.tight_layout()
        for ext, dpi in ((".png", ps.PNG_DPI), (".pdf", None)):
            fig.savefig((OUT / "fig_dw_halpha_depth").with_suffix(ext), dpi=dpi)
        plt.close(fig)
    print(f"figure: {OUT / 'fig_dw_halpha_depth.png'}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("census", "fetch", "table", "figure", "all"):
        sub.add_parser(name)
    m = sub.add_parser("measure")
    m.add_argument("--max-seconds", type=float, default=0.0,
                   help="return cleanly after this many seconds (resumable)")
    args = ap.parse_args(argv)
    if args.cmd == "all":
        cmd_census(); cmd_measure(argparse.Namespace(max_seconds=0)); cmd_table(); cmd_figure()
    else:
        {"census": cmd_census, "fetch": cmd_fetch, "measure": cmd_measure,
         "table": cmd_table, "figure": cmd_figure}[args.cmd](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
