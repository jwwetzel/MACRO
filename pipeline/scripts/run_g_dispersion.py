#!/usr/bin/env python
"""D1 / G-1 — settle the grism dispersion on hot stars, once per
(grism, mechanical epoch), and write the fixed solution.

WHAT THIS SCRIPT ANSWERS
------------------------
Is the high-resolution grism ~0.47 A/px (observational astronomer, OA.E1)
or ~1.59 A/px (the first-generation code)?  And what is the low-resolution
unit, for which three inconsistent values were on disk (PH.P7)?  The
answer decides whether line-profile and V/R science is in scope for the
T CrB and Be-star projects, so it is settled by measurement on stars
whose spectra cannot be misread — A and B stars with clean continua.

STAGES (each resumable; ``--all`` chains them)
----------------------------------------------
``--census``   pick calibrator frames from the manifest (READ-ONLY): the
               registry of hot stars below, up to ``PER_STAR`` frames per
               (star, grism, mechanical epoch), spread over nights.
``--extract``  reduce them (trace, polynomial sky, optimal extraction)
               into the spectrum cache and ``g_frames``.
``--identify`` identify each frame's absorption features under the adopted
               seed AND under every rival dispersion that was on file
               (same features, same tolerance) -> ``g_line_id``.
``--lines``    measure every line / band edge of the list in every frame,
               seeded by the frame's own identification -> ``g_line_meas``.
``--solve``    one dispersion per (grism, epoch) from STELLAR lines only,
               free zero point per frame -> ``g_dispersion`` and
               ``g_dispersion_resid``; telluric band-edge effective
               wavelengths on that scale -> ``g_telluric_lambda``.
``--figures``  the evidence plots (spectra with line identifications
               under both hypotheses; hypothesis scores; residuals).

Outputs: ``products/grism/grism.sqlite`` (tables above), spectra in
``products/grism/spec1d/``, figures in ``committee/work/grism/figures/``.
The archive and the shared manifest are opened read-only.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_grism import config as gconfig                # noqa: E402
from macro_grism import db as gdb                        # noqa: E402
from macro_grism import linecal as lc                    # noqa: E402
from macro_grism import store as gstore                  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
FIG_DIR = REPO_ROOT / "committee" / "work" / "grism" / "figures"

# ---------------------------------------------------------------------------
# The calibrator registry
# ---------------------------------------------------------------------------
#: Hot stars with frames on disk.  ``like`` matches lower(target_best)
#: (the manifest's exposure-stratum suffixes, e.g. 'PHECDA hrg 0-51s',
#: are covered by the trailing %).  ``spt`` selects the line list: 'A'
#: uses the Balmer lines and O I 7774, 'B' adds He I.  ``note`` records
#: why a star is (or is not) trusted for which purpose.  Radial
#: velocities are NOT needed by the primary fit (a common velocity is
#: absorbed by the per-frame zero point; see linecal.py) and are not
#: listed.
CALIBRATORS = (
    # name,        like,            spt, note
    ("Vega",       "alpha lyr%",    "A", "A0 V; QHY epochs"),
    ("Vega",       "vega%",         "A", "A0 V; Andor epochs"),
    ("Phecda",     "phecda%",       "A", "A0 Ve; v sin i 178 km/s"),
    ("Phecda",     "gam uma%",      "A", "A0 Ve (same star, IAU name)"),
    ("Alphecca",   "alphecca%",     "A", "A0 V + G5 eclipsing"),
    ("Alphecca",   "hip 76267%",    "A", "A0 V (same star)"),
    ("Denebola",   "denebola%",     "A", "A3 V"),
    ("Denebola",   "hip 57632%",    "A", "A3 V (same star)"),
    ("Rasalhague", "hip 86032%",    "A", "A5 III"),
    ("53 Boo",     "53 boo%",       "A", "A5 V"),
    ("109 Vir",    "109 vir%",      "A", "A0 V"),
    ("tet Vir",    "hr 4963%",      "A", "A1 IV; BeStar standard"),
    ("tet Crt",    "hr 4468%",      "A", "B9.5 Vn; BeStar standard"),
    ("xi2 Cet",    "hr 718%",       "A", "B9 III; BeStar standard"),
    ("Mizar",      "mizar%",        "A", "A2 V SB"),
    ("Spica",      "spica%",        "B", "B1 III-IV SB2: K ~120 km/s, "
                                         "absorbed by the zero point"),
    ("Alkaid",     "hip 67301%",    "B", "B3 V"),
    ("eta Hya",    "hr 3454%",      "B", "B3 V; BeStar standard"),
    ("Algol",      "bet per%",      "B", "B8 V eclipsing"),
    ("Rigel",      "rigel%",        "B", "B8 Ia"),
    ("eta Ori",    "eta ori%",      "B", "B1 V + B2e"),
    ("15 Mon",     "15 mon%",       "B", "O7 V"),
    ("zet Oph",    "zeta oph%",     "B", "O9.5 V"),
    ("tet CrB",    "tet crb%",      "B", "B6 Vnne Be/shell — Halpha "
                                         "contaminated (C4); telluric "
                                         "edges and trace only"),
)

#: Frames per (star, grism, mechanical epoch).
PER_STAR = 10

#: Exposure ceiling for a calibrator frame (s).  Longer exposures of
#: these stars saturate; the short strata are the usable ones.
MAX_EXPTIME = 130.0


def census(mcon) -> list[dict]:
    """The calibrator worklist: up to PER_STAR canonical raw frames per
    (star, grism, mechanical epoch), evenly spread across the available
    frames in night order (so one bad night cannot be the whole sample).
    """
    cols = [c.strip() for c in gstore.MANIFEST_COLS.split(",")]
    groups = defaultdict(list)
    # ONE pass over the manifest (the frames table has no index on
    # filter, and 24 LIKE scans of a 600 MB file on a busy spinning disk
    # take minutes); the star patterns are matched here in Python.
    import fnmatch
    rows = mcon.execute(f"""
        SELECT {gstore.MANIFEST_COLS} FROM frames
        WHERE is_canonical = 1 AND tree = 'rawimage'
          AND lower(filter) IN ('hrg','lrg','hagrism','oggrism','hag')
          AND exptime <= ?
        ORDER BY night, path""", (MAX_EXPTIME,)).fetchall()
    patterns = [(name, like.replace("%", "*"), spt)
                for name, like, spt, _note in CALIBRATORS]
    for r in rows:
        d = dict(zip(cols, r))
        target = (d["target_best"] or "").lower()
        for name, pat, spt in patterns:
            if fnmatch.fnmatchcase(target, pat):
                d["star"], d["spt"] = name, spt
                grism = gconfig.grism_unit(d["filter"])
                epoch = gconfig.mech_epoch_id(d["night"])
                groups[(name, grism, epoch)].append(d)
                break
    work = []
    for key in sorted(groups, key=lambda k: tuple(str(v) for v in k)):
        rows = groups[key]
        idx = np.unique(np.linspace(0, len(rows) - 1,
                                    min(PER_STAR, len(rows))).astype(int))
        for i in idx:
            d = dict(rows[i])
            d["sample"] = "calibrator"
            work.append(d)
    return work


def star_spt() -> dict:
    return {name: spt for name, _like, spt, _note in CALIBRATORS}


#: Stars whose stellar lines never enter the primary fit, with the
#: reason.  Their telluric edges and traces are still used.
#:  * tet CrB — Be/shell star (strategy ruling C4).  BeSS shows pure
#:    Halpha absorption throughout this campaign (novelty-be package,
#:    2026-10-03), so it is probably usable; it is kept out until that is
#:    written into the strategy, and used as an INDEPENDENT check instead.
#:  * Rigel   — B8 Ia supergiant: Halpha is a variable wind (P Cygni)
#:    profile and the He I lines are wind-affected; not a wavelength
#:    reference at any resolution.
HALPHA_CONTAMINATED = ("tet CrB", "Rigel")

#: Frame quality cuts for line work: enough signal, and no column where a
#: native pixel may have clipped (standing rule 4).
MIN_SNR = 25.0

#: Heliocentric radial velocities (km/s) of the calibrators whose
#: velocity is constant to a few km/s, with J2000 coordinates (deg) for
#: the barycentric correction.  Used ONLY to put the telluric band edges
#: on an absolute scale (g_telluric_lambda); the dispersion fit itself
#: needs none of this.  Binaries and variables are deliberately absent.
#: VALUES ARE CATALOGUE NUMBERS QUOTED FROM MEMORY (GCRV / SIMBAD) AND
#: MUST BE VERIFIED BEFORE PUBLICATION; +/- 5 km/s (0.11 A) is carried
#: as their systematic uncertainty.
STAR_RV = {
    "Vega":       (279.2347, 38.7837, -13.9),
    "Phecda":     (178.4577, 53.6948, -12.6),
    "Denebola":   (177.2649, 14.5721, -0.2),
    "Rasalhague": (263.7336, 12.5600, 12.6),
    "109 Vir":    (221.5622, 1.8929, -6.1),
    "tet Vir":    (197.4875, -5.5390, -2.9),
    "tet Crt":    (174.1704, -9.8022, 1.0),
    "xi2 Cet":    (37.0398, 8.4601, 11.9),
    "Alkaid":     (206.8852, 49.3133, -10.9),
    "eta Hya":    (130.8062, 3.3987, 21.0),
}
RV_SYS_KMS = 5.0

#: Winer Observatory (Sonoita, AZ): lat, lon (deg), height (m).
SITE = (31.6657, -110.6018, 1515.7)


def usable_flux(spec: dict) -> np.ndarray:
    """The spectrum the line work uses: optimal flux inside the trace
    extent, NaN elsewhere and in any column flagged saturated."""
    f = np.where(spec["inside"], spec["flux"], np.nan)
    return np.where(spec["n_sat"] > 0, np.nan, f)


def calib_rows(con, where: str = "") -> list[dict]:
    cur = con.execute(f"""
        SELECT path, star, grism, mech_epoch, night, jd, exptime,
               snr_median, n_sat_cols, peak_adu, fwhm_px, nx
        FROM g_frames WHERE sample = 'calibrator' AND status = 'ok'
        {where} ORDER BY mech_epoch, grism, star, path""")
    cols = [c[0] for c in cur.description]
    return [dict(zip(cols, r)) for r in cur]


# ---------------------------------------------------------------------------
# Stage: identification under each hypothesis (D1)
# ---------------------------------------------------------------------------
#: A frame's identification is trusted (seeds its precise measurement,
#: votes on the epoch's dispersion direction) with at least this many
#: matched lines, of which at least this many stellar.
ID_MIN_MATCH = 4
ID_MIN_STELLAR = 2

#: The 1.59 A/px reading puts Hbeta on the detector in every hot-star hrg
#: frame.  Half-window (px) over which its depth is measured there.
HBETA_WIN_PX = 12


def pixel_scale(con, path: str) -> float:
    """Seed rescaling for a frame's camera: reference pitch / its pitch."""
    key = con.execute("SELECT detector_key FROM g_frames WHERE path = ?",
                      (path,)).fetchone()[0]
    det = gconfig.detector_table()[key]
    return gconfig.REFERENCE_PIXEL_UM / det.pixel_um


def hbeta_depth_under_159(norm: np.ndarray, x_ha: float, scale: float):
    """Depth of the spectrum where a 1.59 A/px dispersion (either sign)
    would put Hbeta, given Halpha at ``x_ha``: (depth, error) for the
    DEEPER of the two candidate positions — the reading most favourable
    to that hypothesis.  None when neither position is on usable pixels.
    """
    best = None
    off = (lc.LINE_BY_NAME["Hbeta"].wave - lc.HALPHA_A) / 1.59 * scale
    sm = lc.smooth(norm, 2.5)
    noise = lc.depth_noise(norm, 2.5)
    for sign in (1, -1):
        xp = int(round(x_ha + sign * off))
        lo, hi = xp - HBETA_WIN_PX, xp + HBETA_WIN_PX + 1
        if lo < 0 or hi > len(sm) or not np.isfinite(sm[lo:hi]).all():
            continue
        depth = float(1.0 - np.min(sm[lo:hi]))
        if best is None or depth > best:
            best = depth
    return (best, noise) if best is not None else (None, None)


def run_identify(con) -> None:
    """Identify every usable calibrator frame under the adopted seed and
    under each rival dispersion that was on file — same features, same
    tolerance — and store all scores."""
    spt = star_spt()
    con.execute("DELETE FROM g_line_id")
    con.commit()
    rows_out = []
    for r in calib_rows(con):
        if (r["snr_median"] or 0) < MIN_SNR:
            continue
        spec = gstore.load_spec(r["path"])
        if spec is None:
            continue
        norm, _c, _ok = lc.normalize(usable_flux(spec))
        if np.isfinite(norm).sum() < 500:
            continue
        feats = lc.find_features(norm)
        sc = pixel_scale(con, r["path"])
        st = spt[r["star"]]
        idn = lc.identify_lines(feats, r["grism"], st, pixel_scale=sc)
        if idn is None:
            continue
        rivals = {}
        for name, seed in lc.RIVALS[r["grism"]].items():
            alt = lc.identify_lines(feats, r["grism"], st, seed=seed,
                                    pixel_scale=sc)
            rivals[name] = ({"n_match": alt["n_match"],
                             "n_stellar": alt["n_stellar"],
                             "rms_px": alt["rms_px"],
                             "refined": int(alt["refined"])} if alt else
                            {"n_match": 0, "n_stellar": 0, "rms_px": None,
                             "refined": 0})
        hb, hb_err = (hbeta_depth_under_159(norm, idn["x_halpha"], sc)
                      if r["grism"] == "hrg" else (None, None))
        rows_out.append({
            "path": r["path"], "star": r["star"], "grism": r["grism"],
            "mech_epoch": r["mech_epoch"], "spt": st,
            "n_features": len(feats), "n_match": idn["n_match"],
            "n_stellar": idn["n_stellar"], "n_lines": idn["n_lines"],
            "disp_at_ha": idn["disp_at_ha"], "sign": idn["sign"],
            "x_halpha": idn["x_halpha"], "rms_px": idn["rms_px"],
            "coeffs_json": json.dumps(idn["coeffs"]),
            "refined": int(idn["refined"]),
            "matches_json": json.dumps(idn["matches"]),
            "rivals_json": json.dumps(rivals),
            "best_rival_match": max(v["n_match"] for v in rivals.values()),
            "hbeta_depth_159": hb, "hbeta_depth_159_err": hb_err})
    for row in rows_out:                 # one short transaction
        gdb.upsert(con, "g_line_id", row)
    con.commit()
    print(f"identified: {len(rows_out)} frames")


def epoch_signs(con) -> dict:
    """{(grism, epoch): (sign, n_for, n_against)} — the dispersion
    direction of each (grism, epoch), voted by the frames with a trusted
    identification.  A grism bolted in a wheel has ONE direction per
    mechanical epoch; dissenting frames are counted and reported."""
    votes = defaultdict(list)
    for grism, epoch, sign, nm, ns, ref in con.execute(
            "SELECT grism, mech_epoch, sign, n_match, n_stellar, refined "
            "FROM g_line_id"):
        if ref and nm >= ID_MIN_MATCH + 1 and ns >= ID_MIN_STELLAR + 1:
            votes[(grism, epoch)].append(sign)
    out = {}
    for key, v in votes.items():
        pos, neg = sum(1 for s in v if s > 0), sum(1 for s in v if s < 0)
        out[key] = ((1 if pos >= neg else -1), max(pos, neg), min(pos, neg))
    return out


# ---------------------------------------------------------------------------
# Stage: precise line measurement
# ---------------------------------------------------------------------------
#: Fit half-windows in Angstrom per line (converted to px with the local
#: dispersion of the frame's identification).  Balmer lines in A/B stars
#: have Stark wings tens of Angstrom wide; He I, Si II, Ne I and O I are
#: rotationally broadened lines a few Angstrom wide plus the LSF.
HALFWIN_A = {"Hepsilon": 40.0, "Hdelta": 45.0, "Hgamma": 50.0,
             "Hbeta": 55.0, "Halpha": 55.0,
             "HeI5876": 14.0, "SiII6347": 9.0, "SiII6371": 9.0,
             "NeI6402": 9.0, "HeI6678": 14.0, "HeI7065": 14.0,
             "OI7774": 18.0}
MIN_HALFWIN_PX = 9

#: Acceptance of one line measurement: significance, and proximity to
#: where the frame's own identification puts the line.
LINE_MIN_SNR = 6.0
LINE_MAX_SHIFT_PX = 5.0
LINE_MAX_SHIFT_FRAC = 0.004
EDGE_MIN_DEPTH = 0.02

#: A frame re-identified with the epoch's direction imposed, whose
#: identification did not survive the tight re-match, is still used if
#: its seed-level matches scatter by less than this (px) — an A star on
#: the hrg has three features and cannot be "refined" by a quadratic.
REID_MAX_RMS_PX = 4.0


def measure_frame(flux: np.ndarray, var: np.ndarray, coeffs,
                  grism: str, spt: str) -> list[dict]:
    """Every line / edge of the grism's list that lands on the frame,
    measured.  ``coeffs`` (polyval order in lambda - 6562.80) is the
    frame's own IDENTIFICATION polynomial: it says where to look and how
    wide a window is; the measured centre is constrained to it only by
    the acceptance window (``LINE_MAX_SHIFT_PX`` + a fraction of lever).
    """
    out = []
    x0 = float(np.polyval(coeffs, 0.0))
    for name in lc.GRISM_LINES[grism]:
        ln = lc.LINE_BY_NAME[name]
        if not (ln.spt == "*" or any(c in ln.spt for c in spt)):
            continue
        dl = ln.wave - lc.HALPHA_A
        xp = float(np.polyval(coeffs, dl))
        slope = float(np.polyval(np.polyder(coeffs), dl))     # px / A
        if not (0 <= xp < len(flux)) or slope == 0:
            continue
        tol = LINE_MAX_SHIFT_PX + LINE_MAX_SHIFT_FRAC * abs(xp - x0)
        if ln.shape == "line":
            hw = max(MIN_HALFWIN_PX,
                     int(round(HALFWIN_A[ln.name] * abs(slope))))
            m = lc.fit_line(flux, xp, hw, var=var)
            if m is None or m["snr"] < LINE_MIN_SNR or m["amp_frac"] >= 0:
                continue                     # absent, weak, or emission
            if abs(m["x"] - xp) > tol or m["fwhm_px"] > 1.6 * hw:
                continue
            out.append({"line": ln.name, "kind": ln.kind, "shape": ln.shape,
                        "wave_ref": ln.wave, "x": m["x"],
                        "x_err": m["x_err"], "fwhm_px": m["fwhm_px"],
                        "amp_frac": m["amp_frac"], "snr": m["snr"],
                        "rchi2": m["rchi2"]})
        else:
            span = lc.EDGE_FLOOR_SPAN_A[ln.name] * abs(slope)
            m = lc.edge_position(flux, xp, 1 if slope > 0 else -1,
                                 search_px=tol + 4.0,
                                 span_px=max(4.0, span))
            if m is None or m["depth"] < EDGE_MIN_DEPTH:
                continue
            out.append({"line": ln.name, "kind": ln.kind, "shape": ln.shape,
                        "wave_ref": ln.wave, "x": m["x"],
                        "x_err": m["x_err"], "fwhm_px": None,
                        "amp_frac": -m["depth"], "snr": None,
                        "rchi2": None})
    return out


def run_lines(con) -> None:
    """Measure the line list in every calibrator frame whose
    identification is trusted and agrees with its epoch's direction."""
    spt = star_spt()
    signs = epoch_signs(con)
    for key, (sign, n_for, n_against) in sorted(signs.items(),
                                                key=lambda kv: str(kv[0])):
        print(f"  direction {key[1]:9s} {key[0]}: red toward "
              f"{'+x' if sign > 0 else '-x'}  ({n_for} frames for, "
              f"{n_against} against)")
    con.execute("DELETE FROM g_line_meas")
    con.commit()
    cur = con.execute("SELECT path, star, grism, mech_epoch, sign, n_match, "
                      "refined, coeffs_json FROM g_line_id")
    ids = cur.fetchall()
    rows_out = []
    n_reid = n_drop = 0
    for path, star, grism, epoch, sign, nm, refined, cj in ids:
        key = (grism, epoch)
        if key not in signs:
            continue
        coeffs = json.loads(cj)
        spec = gstore.load_spec(path)
        flux = usable_flux(spec)
        norm, _c, _ok = lc.normalize(flux)
        if sign != signs[key][0] or not refined or nm < ID_MIN_MATCH:
            # Too few lines to fix the direction on its own (an A star on
            # the hrg shows only Halpha and two nearly mirror-image
            # telluric bands), or a wrong-direction match: re-identify
            # with the epoch's direction imposed.
            idn = lc.identify_lines(lc.find_features(norm), grism,
                                    spt[star], sign=signs[key][0],
                                    pixel_scale=pixel_scale(con, path))
            n_reid += 1
            if idn is None or idn["n_match"] < 3 or (
                    not idn["refined"] and idn["rms_px"] > REID_MAX_RMS_PX):
                n_drop += 1
                continue
            coeffs = idn["coeffs"]
        fl = np.where(np.isfinite(norm), flux, np.nan)
        for m in measure_frame(fl, spec["var"], coeffs, grism, spt[star]):
            m.update(path=path, star=star, grism=grism, mech_epoch=epoch,
                     used=0)
            rows_out.append(m)
    print(f"  re-identified with the epoch direction: {n_reid} frames "
          f"({n_drop} dropped: no acceptable identification)")
    for m in rows_out:
        gdb.upsert(con, "g_line_meas", m)
    con.commit()
    print(f"line measurements: {len(rows_out)}")


# ---------------------------------------------------------------------------
# Stage: the solution
# ---------------------------------------------------------------------------
#: Acceptance (SYNTHESIS G-1): at least this many distinct lines and an
#: rms residual below this many pixels.
ACCEPT_MIN_LINES = 3
ACCEPT_RMS_PX = 1.0

#: A frame enters the primary fit with at least this many stellar lines.
MIN_STELLAR_PER_FRAME = 2

#: Residual clip of the primary fit (px): a line more than this far from
#: the solution after a first pass is a misidentification, not noise.
RESID_CLIP_PX = 3.0

#: Floor added in quadrature to every line-centre error (px).  The formal
#: Gaussian-fit errors of a high-S/N Balmer line are ~0.02 px, far below
#: what a non-Gaussian, rotationally broadened, blended line can deliver;
#: without a floor the brightest star's Halpha dictates the fit.
LINE_ERR_FLOOR_PX = 0.25

#: A higher polynomial degree is adopted only if its leading coefficient
#: is significant at this many sigma AND it lowers the rms.
DEGREE_SIGMA = 4.0


_BARY_CACHE: dict = {}


def barycentric_kms(ra_deg: float, dec_deg: float, jd: float) -> float:
    """Barycentric radial-velocity correction (km/s) toward a star at
    Winer for a UTC Julian date (astropy; positive = observatory moving
    toward the star)."""
    key = (round(ra_deg, 4), round(dec_deg, 4), round(jd, 3))
    if key in _BARY_CACHE:
        return _BARY_CACHE[key]
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    loc = EarthLocation.from_geodetic(SITE[1] * u.deg, SITE[0] * u.deg,
                                      SITE[2] * u.m)
    sc = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
    v = sc.radial_velocity_correction(
        "barycentric", obstime=Time(jd, format="jd", scale="utc"),
        location=loc)
    _BARY_CACHE[key] = float(v.to(u.km / u.s).value)
    return _BARY_CACHE[key]


def fit_group(stellar: list[dict], x_ref: float) -> dict:
    """The primary fit of one (grism, epoch): stellar lines, wavelength
    as a polynomial in detector position, free constant per frame.

    The DEGREE is chosen first, on all measurements, by significance:
    starting from a straight line, the next degree is adopted while its
    leading coefficient exceeds ``DEGREE_SIGMA`` sigma (and enough
    distinct lines remain to leave one of redundancy).  Only then are
    residuals clipped, at the adopted degree.  (Clipping before choosing
    the degree — as a first draft did — lets a straight line throw away
    the far lines of a curved solution and then look well fitted.)
    Returns the adopted ``solve_dispersion`` dict plus 'use' (the
    measurements kept) and 'n_clipped'.
    """
    n_lines = len({m["line"] for m in stellar})

    def solve(use, degree):
        err = np.hypot([m["x_err"] for m in use], LINE_ERR_FLOOR_PX)
        return lc.solve_dispersion(
            [m["path"] for m in use], [m["wave_ref"] for m in use],
            [m["x"] for m in use], err, degree=degree, x_ref=x_ref)

    degree = 1
    for trial in (2, 3):
        if n_lines < trial + 2:
            break
        sol = solve(stellar, trial)
        if abs(sol["coeffs"][-1]) > DEGREE_SIGMA * sol["coeff_errs"][-1]:
            degree = trial
        else:
            break
    use = list(stellar)
    for _ in range(4):
        sol = solve(use, degree)
        keep = np.abs(sol["resid_px"]) <= RESID_CLIP_PX
        if keep.all():
            break
        use = [m for m, k in zip(use, keep) if k]
        cnt = defaultdict(int)
        for m in use:
            cnt[m["path"]] += 1
        use = [m for m in use if cnt[m["path"]] >= MIN_STELLAR_PER_FRAME]
    sol = solve(use, degree)
    sol["use"] = use
    sol["n_clipped"] = len(stellar) - len(use)
    return sol


#: Epochs whose stellar-only solution is good enough to DEFINE the
#: secondary (telluric) standards: adopted, and at least this many frames
#: of velocity-known stars measuring the edge.
SECONDARY_MIN_FRAMES = 5

#: Error floor (px) of a telluric edge used as a secondary standard: the
#: edge is a half-depth crossing of a band whose depth varies with
#: airmass, not a line centre.
EDGE_ERR_FLOOR_PX = 0.5


def run_solve(con) -> None:
    """Two passes.

    PASS 1 — stellar lines only, per (grism, epoch): the primary
    solution, and on it the effective wavelengths of the telluric band
    edges (stars of known constant velocity only).

    PASS 2 — high-resolution grism only: the O2-gamma and O2-B edge
    wavelengths, averaged over the epochs whose pass-1 solution was
    ADOPTED, become secondary standards and every epoch is re-fitted
    with them added.  This is what lets the A stars contribute (on the
    hrg they show one stellar line, Halpha, but both edges), and it is
    what rescues an epoch with few B-star frames.  The pass-1 numbers
    are kept in ``g_meta`` ('solve_pass1') so the report can show that
    pass 2 did not move the well-determined epochs.  The low-resolution
    grism stays stellar-only: its telluric edges lie redward of every
    stellar line, where the polynomial is an extrapolation and the edge
    wavelength at ~15 A resolution is not known a priori.
    """
    summary1, tel1, _ = _solve_all(con, secondary=None)
    sec = {}
    for feat in ("O2gamma", "O2B"):
        vals = [(t["wave_eff"], t["wave_eff_err"]) for t in tel1
                if t["grism"] == "hrg" and t["feature"] == feat
                and t["n_frames"] >= SECONDARY_MIN_FRAMES
                and summary1.get(("hrg", t["mech_epoch"]),
                                 {}).get("status") == "adopted"]
        if vals:
            w = np.array([1.0 / v[1] ** 2 for v in vals])
            sec[feat] = float(np.sum(w * np.array([v[0] for v in vals]))
                              / w.sum())
    print(f"  secondary standards (hrg): "
          f"{ {k: round(v, 2) for k, v in sec.items()} }")
    summary2, _tel2, writes = _solve_all(con, secondary=sec)
    for t in ("g_dispersion", "g_dispersion_resid", "g_telluric_lambda"):
        con.execute(f"DELETE FROM {t}")
    con.execute("UPDATE g_line_meas SET used = 0")
    used = [(w[1]["path"], w[1]["line"]) for w in writes
            if w[0] == "g_dispersion_resid"]
    for table, row in writes:
        gdb.upsert(con, table, row)
    con.executemany("UPDATE g_line_meas SET used = 1 "
                    "WHERE path = ? AND line = ?", used)
    gdb.set_g_meta(con, "solve_pass1", json.dumps(
        {f"{k[0]}|{k[1]}": v for k, v in summary1.items()}))
    gdb.set_g_meta(con, "secondary_standards", json.dumps(sec))
    con.commit()


def _solve_all(con, secondary):
    """One solving pass over every (grism, epoch).  ``secondary`` maps
    telluric feature -> wavelength to use as an extra line on the hrg
    (None: stellar lines only).  Returns (summary, telluric rows,
    pending table writes); nothing is written here."""
    label = "pass 2 (stellar + telluric)" if secondary else \
        "pass 1 (stellar only)"
    print(f"  --- {label}")
    cur = con.execute("""
        SELECT m.path, m.star, m.grism, m.mech_epoch, m.line, m.kind,
               m.wave_ref, m.x, m.x_err, f.jd, f.nx
        FROM g_line_meas m JOIN g_frames f USING (path)""")
    cols = [c[0] for c in cur.description]
    meas = [dict(zip(cols, r)) for r in cur]
    by_key = defaultdict(list)
    for m in meas:
        by_key[(m["grism"], m["mech_epoch"])].append(m)
    writes, summary, tel_rows = [], {}, []
    for (grism, epoch), rows in sorted(by_key.items(),
                                       key=lambda kv: str(kv[0])):
        stellar = [dict(m) for m in rows if m["kind"] == "stellar"
                   and m["star"] not in HALPHA_CONTAMINATED]
        if secondary and grism == "hrg":
            # Telluric edges of velocity-known stars as extra lines, with
            # the standard's wavelength moved INTO the star's rest frame
            # (the frame constant is defined by the stellar lines).
            for m in rows:
                if (m["kind"] == "telluric" and m["line"] in secondary
                        and m["star"] in STAR_RV):
                    ra, dec, rv = STAR_RV[m["star"]]
                    v = rv - barycentric_kms(ra, dec, m["jd"])
                    e = dict(m)
                    e["wave_ref"] = secondary[m["line"]] / (
                        1.0 + v / 299792.458)
                    e["x_err"] = float(np.hypot(m["x_err"],
                                                EDGE_ERR_FLOOR_PX))
                    stellar.append(e)
        per_frame = defaultdict(list)
        for m in stellar:
            per_frame[m["path"]].append(m)
        stellar = [m for p, ms in per_frame.items()
                   if len(ms) >= MIN_STELLAR_PER_FRAME for m in ms]
        n_lines = len({m["line"] for m in stellar})
        if n_lines < 3 or len(stellar) < 6:
            writes.append(("g_dispersion", {
                "grism": grism, "mech_epoch": epoch, "status": "unsolved",
                "n_meas": len(stellar), "n_lines": n_lines,
                "note": "fewer than three distinct lines"}))
            summary[(grism, epoch)] = {"status": "unsolved"}
            print(f"  {epoch:9s} {grism}: UNSOLVED ({n_lines} lines)")
            continue
        # Reference column: the detector centre (median frame width / 2).
        x_ref = float(np.median([m["nx"] for m in rows])) / 2.0
        sol = fit_group(stellar, x_ref)
        use = sol["use"]
        for m, rp, ra in zip(use, sol["resid_px"], sol["resid_a"]):
            writes.append(("g_dispersion_resid", {
                "grism": grism, "mech_epoch": epoch, "path": m["path"],
                "star": m["star"], "line": m["line"],
                "wave_ref": m["wave_ref"], "x": m["x"], "x_err": m["x_err"],
                "resid_px": float(rp), "resid_a": float(ra)}))
        lines_used = sorted({m["line"] for m in use})
        stars_used = sorted({m["star"] for m in use})
        ok = (len(lines_used) >= ACCEPT_MIN_LINES
              and sol["rms_px"] < ACCEPT_RMS_PX)

        # ---- telluric band edges on the stellar scale ------------------
        const = sol["const"]
        tel = defaultdict(list)
        for m in rows:
            if (m["kind"] != "telluric" or m["path"] not in const
                    or m["star"] not in STAR_RV):
                continue
            ra, dec, rv = STAR_RV[m["star"]]
            # The frame constant was fixed by STELLAR lines, so the scale
            # is in the star's rest frame.  The star recedes from the
            # observer at v = rv - bary; a telluric feature at rest in
            # the observatory reads lambda_scale = lambda_obs / (1 + v/c),
            # i.e. lambda_obs = lambda_scale * (1 + v/c).
            v = rv - barycentric_kms(ra, dec, m["jd"])
            lam = const[m["path"]] + float(
                lc.scale_poly(m["x"], sol["coeffs"], x_ref))
            tel[m["line"]].append((lam * (1.0 + v / 299792.458),
                                   m["star"]))
        o2b, o2b_err = None, None
        for feat, vals in tel.items():
            lam = np.array([v[0] for v in vals])
            med = float(np.median(lam))
            scat = float(1.4826 * np.median(np.abs(lam - med)))
            err = float(np.hypot(scat / np.sqrt(len(lam)),
                                 RV_SYS_KMS / 299792.458 * med))
            writes.append(("g_telluric_lambda", {
                "grism": grism, "mech_epoch": epoch, "feature": feat,
                "wave_head": lc.LINE_BY_NAME[feat].wave,
                "wave_eff": med, "wave_eff_err": err, "scatter_a": scat,
                "n_frames": len(lam),
                "n_stars": len({v[1] for v in vals})}))
            tel_rows.append({"grism": grism, "mech_epoch": epoch,
                             "feature": feat, "wave_eff": med,
                             "wave_eff_err": err, "n_frames": len(lam)})
            if feat == "O2B":
                o2b, o2b_err = med, err
        if secondary and grism == "hrg" and "O2B" in secondary:
            # The check on science frames uses the STANDARD's wavelength,
            # identical for every epoch of this grism.
            o2b, o2b_err = secondary["O2B"], o2b_err
        d_lo, d_hi = (abs(float(lc.local_dispersion(x_ref + dx,
                                                    sol["coeffs"], x_ref)))
                      for dx in (-1000.0, 1000.0))
        writes.append(("g_dispersion", {
            "grism": grism, "mech_epoch": epoch,
            "disp_a_per_px": sol["disp_ref"],
            "disp_err": sol["disp_ref_err"],
            "disp_minus1000": d_lo, "disp_plus1000": d_hi, "x_ref": x_ref,
            "coeffs_json": json.dumps(sol["coeffs"]),
            "coeff_errs_json": json.dumps(sol["coeff_errs"]),
            "degree": sol["degree"],
            "n_frames": sol["n_frames"], "n_stars": len(stars_used),
            "n_lines": len(lines_used), "n_meas": sol["n_meas"],
            "n_clipped": sol["n_clipped"],
            "rms_px": sol["rms_px"], "rms_a": sol["rms_a"],
            "max_abs_resid_px": float(np.max(np.abs(sol["resid_px"]))),
            "chi2": sol["chi2"], "dof": sol["dof"], "rchi2": sol["rchi2"],
            "o2b_wave_eff": o2b, "o2b_wave_err": o2b_err,
            "stars": ", ".join(stars_used), "lines": ", ".join(lines_used),
            "status": "adopted" if ok else "provisional",
            "note": "" if ok else
            f"acceptance not met: {len(lines_used)} lines, rms "
            f"{sol['rms_px']:.2f} px"}))
        summary[(grism, epoch)] = {
            "disp_ref": sol["disp_ref"], "disp_err": sol["disp_ref_err"],
            "rms_px": sol["rms_px"], "degree": sol["degree"],
            "n_frames": sol["n_frames"], "n_lines": len(lines_used),
            "disp_minus1000": d_lo, "disp_plus1000": d_hi,
            "status": "adopted" if ok else "provisional"}
        print(f"  {epoch:9s} {grism}: D(x_ref={x_ref:.0f}) = "
              f"{sol['disp_ref']:+.5f} +/- {sol['disp_ref_err']:.5f} A/px  "
              f"(deg {sol['degree']}; {len(lines_used)} lines, "
              f"{len(stars_used)} stars, {sol['n_frames']} frames, "
              f"{sol['n_clipped']} clipped, rms {sol['rms_px']:.3f} px, "
              f"chi2/dof {sol['chi2']:.1f}/{sol['dof']}; |D| at -1000/"
              f"+1000 px = {d_lo:.3f}/{d_hi:.3f})  "
              f"{'ADOPTED' if ok else 'provisional'}")
    return summary, tel_rows, writes


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--manifest", default=str(gstore.default_manifest()))
    ap.add_argument("--archive", default=str(gstore.DEFAULT_ARCHIVE))
    ap.add_argument("--db", default=str(gconfig.GRISM_DB))
    ap.add_argument("--workers", type=int, default=4)
    for stage in ("census", "extract", "identify", "lines", "solve",
                  "figures", "all"):
        ap.add_argument(f"--{stage}", action="store_true")
    args = ap.parse_args(argv)

    con = gdb.connect_grism(args.db)
    mcon = gdb.connect_manifest_ro(args.manifest)
    work = census(mcon)
    if args.census or args.all:
        from collections import Counter
        c = Counter((w["star"], gconfig.grism_unit(w["filter"]),
                     gconfig.mech_epoch_id(w["night"])) for w in work)
        for k in sorted(c, key=lambda k: (str(k[2]), str(k[1]), k[0])):
            print(f"  {str(k[2]):9s} {str(k[1]):3s} {k[0]:11s} {c[k]:3d}")
        print(f"calibrator worklist: {len(work)} frames")
    if args.extract or args.all:
        n = gstore.reduce_batch(con, work, Path(args.archive),
                                workers=args.workers)
        print(f"reduced {n} frames")
    if args.identify or args.all:
        run_identify(con)
    if args.lines or args.all:
        run_lines(con)
    if args.solve or args.all:
        run_solve(con)
    if args.figures or args.all:
        from macro_grism import figures_g
        for f in figures_g.dispersion_figures(con, FIG_DIR):
            print(f"figure: {f}")
    gdb.set_g_meta(con, "dispersion_code", gstore.REDUCE_VERSION)
    con.commit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
