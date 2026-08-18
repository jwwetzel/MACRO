"""Gaia DR3 field identification: cone query (I/O) + pure projection/match.

Role in S4 (prototype-grade, stated honestly in the report): the ensemble
magnitudes are INSTRUMENTAL.  Gaia supplies three services —

1.  absolute identity for the matched stars (source_id + G magnitude),
2.  an independent plate-scale/orientation check (the similarity transform
    fitted between reference-frame pixels and Gaia tangent-plane
    coordinates measures the scale from the sky, not from header cards),
3.  a single zero-point OFFSET per (target, era, filter) — the median of
    (G - ensemble mag) over comparison stars — that places light curves
    near a recognizable scale WITHOUT claiming absolute calibration
    (no color terms; g/r/i vs broad G differ star by star at the tenths
    level).  Absolute calibration is a later S4 milestone.

The pure pieces (tangent projection, parity handling, offset statistics)
are unit-tested; the ADQL query is the only network call in the package.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

# --------------------------------------------------------------------------
# Tunables.
# --------------------------------------------------------------------------

#: Gaia magnitude limit for the cone (deep enough to cover every ensemble
#: star at these exposure times; shallow enough to keep the match sparse).
GAIA_G_MAX = 19.0

#: Match tolerance between transformed reference stars and Gaia positions,
#: arcsec (generous: proper motions over ~9 years and centroid noise).
GAIA_MATCH_TOL_ARCSEC = 2.0

#: Tolerance for identifying THE TARGET among reference stars, arcsec.
#: Wider than the star-to-star tolerance on purpose: the similarity fit is
#: linear over the whole frame, so residual field distortion reaches a few
#: arcsec at the ~1e-3 level over a half-degree field, and the target ID
#: is a single named-coordinate lookup, not a statistical cross-match.
#: The achieved separation is recorded in ``phot_gaia_tie`` either way —
#: an ID near this bound is visibly weaker evidence than one at 1".
TARGET_ID_TOL_ARCSEC = 4.0

#: Target coordinates fallback (ICRS deg) used ONLY if SIMBAD is
#: unreachable at build time; values from SIMBAD (2026-08-17): AN UMa
#: 11:04:25.68 +45:03:13.9, VV Pup 08:15:06.79 -19:03:17.4.  The DB records
#: which source (live query vs fallback) supplied the numbers.
TARGET_COORDS_FALLBACK: dict[str, tuple[float, float]] = {
    "anuma": (166.10700, 45.05386),
    "vvpup": (123.77829, -19.05483),
}


# --------------------------------------------------------------------------
# Pure geometry
# --------------------------------------------------------------------------

def tangent_project(ra_deg: np.ndarray, dec_deg: np.ndarray,
                    ra0_deg: float, dec0_deg: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Gnomonic (tangent-plane) projection, output in arcsec.

    Standard coordinates (xi, eta) about the tangent point; xi grows toward
    INCREASING RA (east) and eta toward increasing Dec.  Detector x may run
    either way relative to east — the identification tries both parities
    (:func:`parity_candidates`) rather than assuming one.
    """
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    ra0, dec0 = math.radians(ra0_deg), math.radians(dec0_deg)
    cosc = (np.sin(dec0) * np.sin(dec)
            + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec)
           - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(xi) * 3600.0, np.degrees(eta) * 3600.0


def parity_candidates(xi: np.ndarray, eta: np.ndarray
                      ) -> list[tuple[str, np.ndarray]]:
    """The two sky->detector parity hypotheses astroalign must try.

    astroalign fits similarity transforms (rotation + scale + shift) but
    never reflections, so a mirror-flipped detector needs the flip applied
    up front.  Returns [('direct', pts), ('flipped', pts)] with xi negated
    in the flipped hypothesis.
    """
    direct = np.column_stack([xi, eta])
    flipped = np.column_stack([-xi, eta])
    return [("direct", direct), ("flipped", flipped)]


def median_offset(gaia_g: np.ndarray, ens_mag: np.ndarray
                  ) -> tuple[float, float, int]:
    """Median (G - ensemble) offset, its robust scatter, and the star count.

    The scatter is 1.4826 x MAD — with no color terms it is dominated by
    star-to-star color differences, and the report quotes it as exactly
    that: the honest width of a color-termless tie.
    """
    g = np.asarray(gaia_g, dtype=float)
    m = np.asarray(ens_mag, dtype=float)
    ok = np.isfinite(g) & np.isfinite(m)
    if ok.sum() == 0:
        return float("nan"), float("nan"), 0
    d = g[ok] - m[ok]
    med = float(np.median(d))
    mad = float(1.4826 * np.median(np.abs(d - med)))
    return med, mad, int(ok.sum())


# --------------------------------------------------------------------------
# I/O: the cone query (the package's only network call)
# --------------------------------------------------------------------------

def cone_query(ra_deg: float, dec_deg: float, radius_deg: float,
               g_max: float = GAIA_G_MAX) -> dict:
    """Gaia DR3 cone: source_id, ra, dec, G — parallel numpy arrays.

    Synchronous ADQL job.  The Gaia archive caps synchronous queries at
    2,000 rows NO MATTER WHAT — and in the VV Pup field (galactic latitude
    ~ +2 deg) a G < 19 cone holds far more than that, so an uncapped query
    comes back truncated in ARBITRARY order and the field fit starves.
    ``SELECT TOP 2000 ... ORDER BY phot_g_mean_mag`` turns the cap into a
    brightness cut: we always get the BRIGHTEST 2,000 — exactly the stars
    the triangle fit and the ensemble tie need.  Raises on network
    failure; the build script catches and records the outage so the
    pipeline still completes with an instrumental-only zero point.
    """
    from astroquery.gaia import Gaia         # deferred: import cost + net
    Gaia.ROW_LIMIT = -1
    job = Gaia.launch_job(f"""
        SELECT TOP 2000 source_id, ra, dec, phot_g_mean_mag
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),
                           CIRCLE('ICRS', {ra_deg:.6f}, {dec_deg:.6f},
                                  {radius_deg:.6f}))
          AND phot_g_mean_mag < {g_max:.2f}
        ORDER BY phot_g_mean_mag""")
    tab = job.get_results()
    return {
        "source_id": np.asarray(tab["source_id"], dtype=np.int64),
        "ra": np.asarray(tab["ra"], dtype=float),
        "dec": np.asarray(tab["dec"], dtype=float),
        "gmag": np.asarray(tab["phot_g_mean_mag"], dtype=float),
    }


def resolve_target(target_key: str) -> tuple[float, float, str]:
    """Target ICRS coordinates: live SIMBAD if reachable, else the recorded
    fallback.  Returns (ra_deg, dec_deg, source_tag)."""
    name = {"anuma": "AN UMa", "vvpup": "VV Pup"}.get(target_key, target_key)
    try:
        from astroquery.simbad import Simbad
        tab = Simbad.query_object(name)
        # astroquery >= 0.4.8 returns decimal-degree 'ra'/'dec' columns.
        ra = float(tab["ra"][0]); dec = float(tab["dec"][0])
        return ra, dec, "simbad"
    except Exception:
        ra, dec = TARGET_COORDS_FALLBACK[target_key]
        return ra, dec, "fallback_constant"


def identify_reference(ref_xy: np.ndarray, gaia: dict,
                       ra0: float, dec0: float,
                       ref_bright_xy: Optional[np.ndarray] = None,
                       fit_radius_arcsec: Optional[float] = None,
                       seed: Optional[int] = None,
                       ) -> Optional[dict]:
    """Fit the reference frame to Gaia and match every reference star.

    The fit walks the SAME attempt ladder as frame-to-frame matching
    (``macro_phot.extract.ALIGN_ATTEMPTS``) for each parity hypothesis:
    sparse high-latitude fields converge on the first small rung, while
    the dense VV Pup field — where the brightest-2,000 Gaia cone bottoms
    out near G ~ 16 — only locks with the deep 800-star pools.
    ``ref_bright_xy`` must be the UNCLIPPED reference detections in
    flux-descending order: the symmetric-pool lesson from frame matching
    applies here identically (saturated stars present on one side of the
    fit and absent from the other poison the control sets from the top).
    The transform is fitted on bright subsets but APPLIED to every
    reference star, so faint stars still receive identities.
    ``fit_radius_arcsec`` confines the Gaia fit pool to stars near the
    tangent point (pass the frame's inscribed-circle radius): the cone is
    deliberately wider than the frame, and bright off-frame stars would
    otherwise crowd the fit pool with unmatchable triangles.
    ``seed`` (the build script passes the reference frame_id) pins
    astroalign's otherwise-unseedable RANSAC for the whole parity/ladder
    walk (``macro_phot.extract.seeded_ransac``), so the fitted tie — and
    the Gaia identity of every star — reproduces run to run.  Returns
    None when no parity converges.  On success:

    ``{'parity', 'scale_arcsec_per_px', 'rot_deg', 'n_matched',
       'gaia_idx' (per ref star, -1 = unmatched),
       'ref_radec' (per ref star, transformed sky position deg)}``
    """
    import astroalign as aa
    from . import photometry as ph
    from .extract import ALIGN_ATTEMPTS, seeded_ransac
    xi, eta = tangent_project(gaia["ra"], gaia["dec"], ra0, dec0)
    if ref_bright_xy is None:
        ref_bright_xy = ref_xy
    # Footprint-confined, brightness-ordered Gaia fit pool.
    pool = np.arange(len(xi))
    if fit_radius_arcsec is not None:
        inside = np.hypot(xi, eta) <= fit_radius_arcsec
        if inside.sum() >= 75:            # enough stars to fit from
            pool = np.flatnonzero(inside)
    gaia_order = pool[np.argsort(gaia["gmag"][pool])]
    best = None
    with seeded_ransac(seed):
        for parity, sky in parity_candidates(xi, eta):
            tf = None
            for n_pool, n_ctrl in ALIGN_ATTEMPTS:
                try:
                    # source = pixel positions, target = sky (arcsec): the
                    # fitted similarity scale IS the plate scale in "/px.
                    tf, _ = aa.find_transform(ref_bright_xy[:n_pool],
                                              sky[gaia_order[:n_pool]],
                                              max_control_points=n_ctrl)
                    break
                except Exception:
                    continue
            if tf is None:
                continue
            moved = tf(ref_xy)           # every ref star onto the sky plane
            idx = ph.match_one_to_one(sky, moved, GAIA_MATCH_TOL_ARCSEC)
            n = int((idx >= 0).sum())
            if best is None or n > best["n_matched"]:
                best = {"parity": parity,
                        "scale_arcsec_per_px": float(tf.scale),
                        "rot_deg": float(math.degrees(tf.rotation)),
                        "n_matched": n, "gaia_idx": idx, "moved": moved}
    if best is None or best["n_matched"] < 6:
        return None
    # Convert the matched tangent coordinates back to RA/Dec for storage.
    xi_m, eta_m = best["moved"][:, 0], best["moved"][:, 1]
    if best["parity"] == "flipped":
        xi_m = -xi_m
    ra, dec = tangent_deproject(xi_m, eta_m, ra0, dec0)
    best["ref_radec"] = np.column_stack([ra, dec])
    del best["moved"]
    return best


def tangent_deproject(xi_arcsec: np.ndarray, eta_arcsec: np.ndarray,
                      ra0_deg: float, dec0_deg: float
                      ) -> tuple[np.ndarray, np.ndarray]:
    """Inverse gnomonic projection (pure; exact inverse of tangent_project)."""
    xi = np.radians(np.asarray(xi_arcsec, dtype=float) / 3600.0)
    eta = np.radians(np.asarray(eta_arcsec, dtype=float) / 3600.0)
    ra0, dec0 = math.radians(ra0_deg), math.radians(dec0_deg)
    den = np.cos(dec0) - eta * np.sin(dec0)
    ra = ra0 + np.arctan2(xi, den)
    dec = np.arctan((np.sin(dec0) + eta * np.cos(dec0))
                    / np.hypot(xi, den))
    return np.degrees(ra) % 360.0, np.degrees(dec)
