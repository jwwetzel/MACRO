"""Pure S4 photometry logic: plate scales, apertures, matching, magnitudes.

Everything in this module is a pure function — no I/O, no globals mutated —
so every rule the photometry obeys can be unit-tested on hand-built cases
(``pipeline/tests/test_phot.py``).  The extraction worker
(``macro_phot.extract``) and the build script wire these functions to real
pixels.

The binding conventions:

1.  The aperture is a fixed angle on the SKY (arcsec), converted to pixels
    per frame from that frame's own XPIXSZ/FOCALLEN header cards — so the
    same physical aperture is used across camera eras with different plate
    scales (0.45"/px for the binned CMOS era, 0.81"/px for the Andor iKon).
2.  Cross-frame star matching never assumes a WCS (most polar frames are
    unsolved — that is S1's business, not S4's).  Frames are registered to
    a per-(target, era) REFERENCE frame by astroalign triangle matching;
    the pure one-to-one nearest-neighbour assignment lives here.
3.  Instrumental magnitudes are exposure-time-normalized
    (-2.5 log10(flux / exptime) + INST_MAG_OFFSET), so frames of different
    exposure times within one series land on one scale before the ensemble
    solves for the residual per-frame zero point.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Arcseconds per radian / 1000: plate scale ["/px] = RAD2ARCSEC_KILO *
#: pixel_size[um] / focal_length[mm].  (206265 arcsec/rad; the um/mm unit
#: pair contributes the factor 1/1000.)
RAD2ARCSEC_KILO = 206.265

#: Photometric aperture RADIUS on the sky.  Winer seeing is typically
#: 2.5-4" FWHM; a 4" radius (~2.9x a 2.8"-FWHM sigma) captures >99% of a
#: Gaussian PSF while keeping the sky contribution modest.
APERTURE_RADIUS_ARCSEC = 4.0

#: Sky annulus (inner, outer) radii on the sky, arcsec.  Wide enough to
#: dodge the PSF wings, narrow enough to sample LOCAL sky.
SKY_ANNULUS_ARCSEC = (8.0, 12.0)

#: Guard rails on the pixel aperture: below ~2 px aperture photometry is
#: dominated by pixelation; above 25 px something is wrong with the header
#: plate-scale cards and we want a loud failure, not a silent 60-px hole.
APERTURE_MIN_PX, APERTURE_MAX_PX = 2.0, 25.0

#: Source-detection threshold in units of the global background RMS.
DETECT_SIGMA = 5.0

#: Minimum connected pixels above threshold for a detection (sep minarea).
DETECT_MINAREA = 5

#: Detections with a peak above this ADU level are flagged as clipped /
#: non-linear and excluded from the ensemble (both cameras digitize 16 bits;
#: the CMOS reduced frames top out at 65535, the iKon at ~65535 — 55000
#: keeps a safe margin below either full well without vetoing real stars).
PEAK_CLIP_ADU = 55000.0

#: One-to-one match tolerance between transformed frame detections and
#: reference stars, in units of the reference frame's FWHM (a seeing-scaled
#: tolerance survives both plate scales without retuning).
MATCH_TOL_FWHM = 0.8

#: Hard floor on the match tolerance in pixels (undersampled frames).
MATCH_TOL_MIN_PX = 2.0

#: Constant added to instrumental magnitudes so typical values are positive
#: (pure display convention; the ensemble is invariant to it).
INST_MAG_OFFSET = 25.0

#: 2.5 / ln(10): converts a relative flux error to a magnitude error.
MAG_ERR_FACTOR = 2.5 / math.log(10.0)

#: Reference-frame choice: only frames whose FWHM sits at or below this
#: percentile of the series compete (sharp frames resolve close pairs), and
#: among those the one with the most detections wins.
REF_FWHM_PERCENTILE = 75.0

#: Minimum stars astroalign needs for a credible triangle match; frames
#: with fewer detections are marked unmatchable rather than force-fitted.
MIN_STARS_FOR_ALIGN = 8

#: Reference-frame quality control against DOUBLE-IMAGED exposures (guiding
#: jump mid-exposure): every star appears twice, at a constant offset, with
#: near-equal component fluxes.  A doubled reference splits every physical
#: star into two identities and silently poisons the whole series (found
#: the hard way: both VV Pup references of the first build were doubled —
#: 83-89% of their stars had an equal-brightness companion within a few
#: FWHM, vs 0.7% for the clean AN UMa reference).  The detector below asks,
#: per star: is there ANOTHER detection within REF_PAIR_RADIUS_FWHM x the
#: frame FWHM whose flux is within a factor REF_PAIR_FLUX_RATIO?  Genuine
#: close binaries/blends put a few percent of stars in pairs; a doubled
#: frame puts nearly all of them there — the two populations are separated
#: by two orders of magnitude, so the threshold is not delicate.
REF_PAIR_RADIUS_FWHM = 3.0
REF_PAIR_FLUX_RATIO = 2.0
REF_DOUBLED_MAX_FRAC = 0.2


# --------------------------------------------------------------------------
# Plate scale and aperture scaling
# --------------------------------------------------------------------------

def plate_scale_arcsec_per_px(xpixsz_um: Optional[float],
                              focallen_mm: Optional[float]) -> Optional[float]:
    """Plate scale in arcsec/pixel from the header pixel size and focal length.

    XPIXSZ is the BINNED pixel size in microns (the acquisition software
    writes physical_pixel x binning), FOCALLEN the telescope focal length in
    millimetres — so this needs no separate binning correction.  Returns
    ``None`` when either card is missing or non-positive (the caller decides
    whether that frame is usable).

    Cross-check available to auditors: the plate-solved CMOS frames carry a
    CD matrix; |CD1_1| = 0.000125 deg = 0.451"/px agrees with
    206.265 * 7.52 / 3454 = 0.449"/px to 0.5%.
    """
    if xpixsz_um is None or focallen_mm is None:
        return None
    if not (xpixsz_um > 0 and focallen_mm > 0):
        return None
    return RAD2ARCSEC_KILO * float(xpixsz_um) / float(focallen_mm)


def aperture_radius_px(plate_scale: Optional[float],
                       radius_arcsec: float = APERTURE_RADIUS_ARCSEC
                       ) -> Optional[float]:
    """Convert the fixed sky aperture to pixels for one frame's plate scale.

    Clamped to [APERTURE_MIN_PX, APERTURE_MAX_PX]: the clamp bounds protect
    against absurd header values, and hitting a bound is itself recorded by
    the caller (the returned value simply equals the bound).
    Returns ``None`` when the plate scale is unknown.
    """
    if plate_scale is None or plate_scale <= 0:
        return None
    return float(np.clip(radius_arcsec / plate_scale,
                         APERTURE_MIN_PX, APERTURE_MAX_PX))


def fwhm_from_ab(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """FWHM estimate from sep's ellipse semi-axes (Gaussian-equivalent).

    sep's ``a``/``b`` are the second-moment RMS sizes along the ellipse
    axes; for a Gaussian profile FWHM = 2 sqrt(2 ln 2) * sigma, and the
    circularized sigma is sqrt((a^2 + b^2) / 2).
    """
    sigma = np.sqrt((np.asarray(a) ** 2 + np.asarray(b) ** 2) / 2.0)
    return 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma


# --------------------------------------------------------------------------
# Reference-frame selection
# --------------------------------------------------------------------------

def rank_references(stats: Sequence[tuple]) -> list[int]:
    """Preference-ordered candidate reference frames of one (target, era).

    Parameters
    ----------
    stats
        Sequence of ``(frame_id, n_detected, fwhm_px, solved)`` tuples, one
        per successfully extracted frame.  ``solved`` is 1 when the frame
        carries a plate solution (nice-to-have for later sky checks — used
        only as a tie-break, never as a requirement, because whole eras are
        unsolved).

    Returns
    -------
    list[int]
        EVERY eligible ``frame_id``, best first (empty for an empty
        series).  The caller walks this list and takes the first candidate
        that passes quality control (:func:`paired_fraction` — a doubled
        reference must be rejected and the next-best tried, which is why
        this returns a ranking rather than a single winner).

    Rule: among frames whose FWHM sits at or below the series'
    :data:`REF_FWHM_PERCENTILE` percentile (sharp frames split close pairs),
    order by most detections; break ties by solved-first, then by smallest
    frame_id (fully deterministic re-runs).
    """
    usable = [(fid, n, f, s) for fid, n, f, s in stats
              if n is not None and n > 0 and f is not None and f > 0]
    if not usable:
        return []
    cutoff = float(np.percentile([f for _, _, f, _ in usable],
                                 REF_FWHM_PERCENTILE))
    sharp = [(fid, n, f, s) for fid, n, f, s in usable if f <= cutoff]
    # sort() over (n_detected, solved, -frame_id), descending, implements
    # the rule; the negated id makes the SMALLEST id win the final tie.
    sharp.sort(key=lambda t: (t[1], t[3] or 0, -t[0]), reverse=True)
    return [t[0] for t in sharp]


def choose_reference(stats: Sequence[tuple]) -> Optional[int]:
    """The single best reference candidate (see :func:`rank_references`).

    Kept as the one-answer convenience wrapper; callers that can VET a
    candidate (the build script, which sees the detection catalog) must use
    the full ranking so a doubled reference can be rejected and replaced.
    """
    ranking = rank_references(stats)
    return ranking[0] if ranking else None


def paired_fraction(xy: np.ndarray, flux: np.ndarray, radius_px: float,
                    flux_ratio: float = REF_PAIR_FLUX_RATIO) -> float:
    """Fraction of catalog stars with a similar-flux companion nearby.

    The double-image detector for reference QC (constants above).  For each
    star: does any OTHER detection sit within ``radius_px`` with flux
    within a factor ``flux_ratio``?  A double-imaged exposure answers yes
    for nearly every star (each physical star contributes two near-equal
    components at one constant offset); a clean field answers yes only for
    the rare genuine blend.  Pure geometry; the caller supplies
    ``radius_px`` = :data:`REF_PAIR_RADIUS_FWHM` x the frame FWHM.

    Computed in row blocks so a 5,000-star dense-field catalog never
    materializes the full pairwise matrix at once.
    """
    xy = np.asarray(xy, dtype=float)
    flux = np.asarray(flux, dtype=float)
    n = len(xy)
    if n < 2:
        return 0.0
    paired = 0
    block = 512
    for i0 in range(0, n, block):
        blk_xy = xy[i0:i0 + block]
        blk_fl = flux[i0:i0 + block]
        # Distances from this block to EVERY star, self-pairs excluded.
        d = np.hypot(blk_xy[:, 0, None] - xy[None, :, 0],
                     blk_xy[:, 1, None] - xy[None, :, 1])
        for k in range(len(blk_xy)):
            d[k, i0 + k] = np.inf
        # Similar flux: min/max ratio of the pair >= 1/flux_ratio.
        with np.errstate(divide="ignore", invalid="ignore"):
            r = (np.minimum(blk_fl[:, None], flux[None, :])
                 / np.maximum(blk_fl[:, None], flux[None, :]))
        paired += int(((d <= radius_px) & (r >= 1.0 / flux_ratio))
                      .any(axis=1).sum())
    return paired / n


# --------------------------------------------------------------------------
# One-to-one nearest-neighbour matching (post-astroalign)
# --------------------------------------------------------------------------

def match_one_to_one(ref_xy: np.ndarray, xy: np.ndarray,
                     tol_px: float) -> np.ndarray:
    """Assign each detection to at most one reference star, greedily by distance.

    Parameters
    ----------
    ref_xy
        (R, 2) reference-star positions, in reference-frame pixels.
    xy
        (N, 2) detection positions ALREADY TRANSFORMED into the reference
        frame (astroalign supplies the transform; this function is pure
        geometry).
    tol_px
        Maximum accepted match distance in reference pixels.

    Returns
    -------
    np.ndarray
        Length-N integer array: the matched reference-star index per
        detection, or -1 for no match.  Greedy assignment in order of
        increasing pair distance guarantees one-to-one matching (a reference
        star claimed by a closer detection is unavailable to a farther one)
        and is deterministic.
    """
    ref_xy = np.asarray(ref_xy, dtype=float)
    xy = np.asarray(xy, dtype=float)
    out = np.full(len(xy), -1, dtype=int)
    if len(ref_xy) == 0 or len(xy) == 0:
        return out
    # All pairwise distances (fields here hold a few hundred stars, so the
    # (N, R) matrix is tiny; no KD-tree needed for correctness or speed).
    d = np.hypot(xy[:, 0, None] - ref_xy[None, :, 0],
                 xy[:, 1, None] - ref_xy[None, :, 1])
    # Candidate pairs within tolerance, sorted by distance (stable sort on
    # the flat index breaks exact-distance ties deterministically).
    cand = np.argwhere(d <= tol_px)
    order = np.argsort(d[cand[:, 0], cand[:, 1]], kind="stable")
    used_det = np.zeros(len(xy), dtype=bool)
    used_ref = np.zeros(len(ref_xy), dtype=bool)
    for i, j in cand[order]:
        if not used_det[i] and not used_ref[j]:
            out[i] = j
            used_det[i] = used_ref[j] = True
    return out


def match_tolerance_px(ref_fwhm_px: Optional[float]) -> float:
    """Seeing-scaled match tolerance with a pixel floor (pure arithmetic)."""
    if ref_fwhm_px is None or ref_fwhm_px <= 0:
        return MATCH_TOL_MIN_PX
    return max(MATCH_TOL_MIN_PX, MATCH_TOL_FWHM * float(ref_fwhm_px))


# --------------------------------------------------------------------------
# Magnitudes
# --------------------------------------------------------------------------

def instrumental_mag(flux: np.ndarray, exptime: float
                     ) -> np.ndarray:
    """Exposure-normalized instrumental magnitude.

    ``-2.5 log10(flux / exptime) + INST_MAG_OFFSET``; non-positive fluxes
    map to NaN (a star measured at or below zero net counts carries no
    magnitude — the ensemble treats NaN as a missing observation).
    """
    flux = np.asarray(flux, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        m = -2.5 * np.log10(flux / float(exptime)) + INST_MAG_OFFSET
    m = np.asarray(m)
    m[~np.isfinite(m)] = np.nan
    return m


def mag_error(flux: np.ndarray, fluxerr: np.ndarray) -> np.ndarray:
    """Magnitude error from a flux and its error: 1.0857 * fluxerr / flux.

    Non-positive fluxes (or errors) map to NaN, mirroring
    :func:`instrumental_mag`.
    """
    flux = np.asarray(flux, dtype=float)
    fluxerr = np.asarray(fluxerr, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = MAG_ERR_FACTOR * fluxerr / flux
    e = np.asarray(e)
    e[(flux <= 0) | (fluxerr <= 0) | ~np.isfinite(e)] = np.nan
    return e
