"""Minimal star photometry for the S2 linearity and saturation probes.

WHY S2 HAS ITS OWN (SMALL) PHOTOMETRY
-------------------------------------
The science photometry lives in ``macro_phot`` and is tuned for differential
light curves: it vetoes anything near the ceiling, because a clipped star is
useless to a light curve.  The detector probes need the opposite — they
exist to watch what happens to a star AS it approaches the ceiling, so the
bright, nearly-clipped and clipped stars are exactly the ones that must be
measured and kept.  They also need one number ``macro_phot`` never stores:
the RAW peak pixel of each star (background included, no smoothing), which
is the quantity a saturation or linearity cap is written in.

So this module is deliberately small and does three things:

* :func:`measure_stars` — detect sources (``sep``), fixed-aperture fluxes
  at several radii with a local sky annulus, and the raw peak in a small
  box at each star;
* :func:`match_shift` / :func:`match_stars` — match two star lists of the
  same field taken minutes apart (a translation found by offset voting,
  then nearest-neighbour matching);
* :func:`native_peak_factor` — what an AVERAGE-binned peak pixel hides:
  the ratio of the brightest NATIVE pixel to the stored binned pixel for a
  star of given width (standing rule 4: saturation is judged in native
  pixels).

No I/O; inputs are arrays, outputs are plain dicts of arrays.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np
import sep

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Detection threshold in background-RMS units.
DETECT_SIGMA = 6.0

#: Minimum connected pixels for a detection.
DETECT_MINAREA = 6

#: sep background mesh (pixels).  Large against a star, small against the
#: sky gradients of a 37-arcminute field.
BKG_MESH = 128

#: The fixed aperture radii (pixels) every frame is measured at, in ONE
#: pass.  A set of frames compared star by star must use one radius, but
#: which radius is right is only known once every frame's seeing has been
#: measured — and on this archive a second read of a frame costs seconds.
#: So all radii are measured at once and the analysis picks, per set, the
#: smallest radius that is at least :data:`APERTURE_FWHM` times the WORST
#: seeing in the set (:func:`choose_aperture`).
APERTURES_PX = (6.0, 9.0, 12.0, 16.0, 22.0)

#: Aperture radius in units of stellar FWHM.  3 FWHM holds > 99% of a
#: Gaussian and most of a Moffat's wings — which is what makes a flux
#: RATIO between two frames of different seeing safe: the ensemble
#: normalisation removes the common aperture loss, the wide aperture keeps
#: the differential one small.
APERTURE_FWHM = 3.0

#: Sky annulus (inner, outer) in units of the LARGEST aperture radius, so
#: one sky value serves every aperture of a star.
ANNULUS_SCALE = (1.4, 2.0)

#: Half-width of the box in which a star's raw peak pixel is searched.
PEAK_BOX_HALF = 3

#: Matching: offset-vote bin (pixels) and the brightest-N stars that vote.
MATCH_VOTE_BIN_PX = 3.0
MATCH_VOTE_NSTARS = 80

#: Matching tolerance after the shift is removed (pixels).
MATCH_TOL_PX = 4.0

#: Brightest sources kept per frame (bounds the work on crowded fields; the
#: linearity question lives entirely in the bright end).
MAX_SOURCES = 1500


def _annulus_median(data: np.ndarray, x: np.ndarray, y: np.ndarray,
                    r_in: float, r_out: float) -> np.ndarray:
    """Median pixel value in each star's sky annulus.

    One small stamp per star; the annulus footprint (relative to the
    stamp centre) is computed once and re-used, with the star's sub-pixel
    position rounded to the nearest pixel — irrelevant for a median over
    several thousand annulus pixels.
    """
    ny, nx = data.shape
    h = int(math.ceil(r_out))
    yy, xx = np.mgrid[-h:h + 1, -h:h + 1]
    r2 = yy * yy + xx * xx
    foot = (r2 >= r_in * r_in) & (r2 <= r_out * r_out)
    out = np.empty(len(x), dtype=np.float64)
    for i, (xi, yi) in enumerate(zip(x, y)):
        cx, cy = int(round(xi)), int(round(yi))
        y0, y1 = cy - h, cy + h + 1
        x0, x1 = cx - h, cx + h + 1
        if y0 < 0 or x0 < 0 or y1 > ny or x1 > nx:
            # Clipped by the frame edge: use what is inside.
            sub = data[max(y0, 0):min(y1, ny), max(x0, 0):min(x1, nx)]
            f = foot[max(y0, 0) - y0: foot.shape[0] - (y1 - min(y1, ny)),
                     max(x0, 0) - x0: foot.shape[1] - (x1 - min(x1, nx))]
            vals = sub[f]
        else:
            vals = data[y0:y1, x0:x1][foot]
        out[i] = float(np.median(vals)) if vals.size else float("nan")
    return out


def choose_aperture(fwhm_px: float,
                    apertures=APERTURES_PX) -> float:
    """Smallest measured radius >= APERTURE_FWHM x the given FWHM.

    Pass the WORST (largest) frame FWHM of a set.  When even the largest
    radius is too small the largest is returned (and the caller records
    the FWHM, so the shortfall is visible).
    """
    need = APERTURE_FWHM * float(fwhm_px)
    for r in sorted(apertures):
        if r >= need:
            return float(r)
    return float(max(apertures))


def measure_stars(img: np.ndarray, apertures=APERTURES_PX,
                  detect_sigma: float = DETECT_SIGMA) -> dict:
    """Detect stars and measure flux (at every aperture), raw peak, width.

    Parameters
    ----------
    img
        Raw frame (ADU).  Not bias-subtracted: the local sky annulus
        removes any pedestal from the flux, and the PEAK is wanted raw.
    apertures
        Aperture radii in pixels; fluxes are returned for all of them.

    Returns
    -------
    dict with per-star arrays ``x``, ``y``, ``peak_raw`` (ADU, maximum raw
    pixel in the peak box), ``sky`` (ADU/pixel, local annulus level),
    ``fwhm`` (px), ``flag`` (sep's flag); 2-D arrays ``flux`` and
    ``flux_err`` of shape (n_apertures, n_stars) (ADU, sky-subtracted
    aperture sums; the error is the sky-RMS floor times sqrt(area) — a
    FLOOR, not a Poisson error: no gain is assumed here); and scalars
    ``apertures`` (the radii), ``fwhm_med`` and ``bkg_rms``.
    """
    data = np.ascontiguousarray(img, dtype=np.float32)
    ny, nx = data.shape
    aps = tuple(float(a) for a in apertures)
    bkg = sep.Background(data, bw=BKG_MESH, bh=BKG_MESH)
    sub = data - bkg.back()
    rms = float(bkg.globalrms)
    empty = {k: np.array([]) for k in
             ("x", "y", "peak_raw", "sky", "fwhm", "flag")}
    empty.update({"flux": np.empty((len(aps), 0)),
                  "flux_err": np.empty((len(aps), 0)),
                  "apertures": aps, "fwhm_med": None, "bkg_rms": rms})
    try:
        objs = sep.extract(sub, detect_sigma, err=rms,
                           minarea=DETECT_MINAREA)
    except Exception:
        # sep raises when its internal pixel buffer overflows on a frame
        # that is one huge source (a twilight frame, a saturated field).
        return empty
    if len(objs) == 0:
        return empty
    order = np.argsort(objs["flux"])[::-1][:MAX_SOURCES]
    objs = objs[order]
    # FWHM from the second moments (Gaussian-equivalent).
    fwhm = 2.0 * np.sqrt(np.log(2.0) * (objs["a"] ** 2 + objs["b"] ** 2))
    good = np.isfinite(fwhm) & (fwhm > 1.0) & (objs["flag"] == 0)
    fwhm_med = float(np.median(fwhm[good])) if good.any() else None
    r_in, r_out = ANNULUS_SCALE[0] * max(aps), ANNULUS_SCALE[1] * max(aps)
    x, y = objs["x"], objs["y"]
    margin = r_out + 1.0
    inside = ((x > margin) & (x < nx - margin)
              & (y > margin) & (y < ny - margin))
    objs, fwhm, x, y = objs[inside], fwhm[inside], x[inside], y[inside]
    if len(objs) == 0:
        empty["fwhm_med"] = fwhm_med
        return empty
    # Local sky from the annulus, on the RAW frame: the flux is then
    # independent of sep's background model near bright stars.  The sky is
    # the MEDIAN of the annulus pixels, not their mean: a neighbour star
    # inside the annulus would raise a mean and so lower this star's flux
    # (1.3% on a 50-star synthetic field, which is how the unit test
    # caught the mean).
    sky = _annulus_median(data, x, y, r_in, r_out)
    flux = np.empty((len(aps), len(objs)))
    flux_err = np.empty_like(flux)
    for k, r in enumerate(aps):
        raw_sum, _, _ = sep.sum_circle(data, x, y, r)
        area = math.pi * r * r
        flux[k] = raw_sum - sky * area
        # Sky noise in the aperture plus the error of the sky estimate.
        ann_area = math.pi * (r_out ** 2 - r_in ** 2)
        flux_err[k] = rms * math.sqrt(area + area * area / ann_area)
    peak = np.empty(len(objs))
    for i, (xi, yi) in enumerate(zip(x, y)):
        cx, cy = int(round(xi)), int(round(yi))
        box = data[max(cy - PEAK_BOX_HALF, 0): cy + PEAK_BOX_HALF + 1,
                   max(cx - PEAK_BOX_HALF, 0): cx + PEAK_BOX_HALF + 1]
        peak[i] = float(box.max())
    return {"x": np.asarray(x, dtype=np.float64),
            "y": np.asarray(y, dtype=np.float64),
            "flux": flux, "flux_err": flux_err,
            "peak_raw": peak, "sky": np.asarray(sky, dtype=np.float64),
            "fwhm": np.asarray(fwhm, dtype=np.float64),
            "flag": np.asarray(objs["flag"], dtype=np.int64),
            "apertures": aps, "fwhm_med": fwhm_med, "bkg_rms": rms}


def match_shift(ref_xy: np.ndarray, xy: np.ndarray,
                bin_px: float = MATCH_VOTE_BIN_PX,
                n_vote: int = MATCH_VOTE_NSTARS,
                max_shift: float = 600.0) -> Optional[tuple[float, float]]:
    """The translation taking ``xy`` onto ``ref_xy``, by offset voting.

    Both lists must be brightest-first.  Every pair among the brightest
    ``n_vote`` stars of each list casts one vote for its offset; true
    pairs all vote for the same cell, false pairs scatter.  Returns the
    refined (dx, dy) — the mean offset of the votes in the winning cell —
    or None when no cell collects at least 4 votes (too few common stars).
    """
    a = np.asarray(ref_xy, dtype=np.float64)[:n_vote]
    b = np.asarray(xy, dtype=np.float64)[:n_vote]
    if len(a) < 4 or len(b) < 4:
        return None
    dx = (a[:, None, 0] - b[None, :, 0]).ravel()
    dy = (a[:, None, 1] - b[None, :, 1]).ravel()
    keep = (np.abs(dx) <= max_shift) & (np.abs(dy) <= max_shift)
    dx, dy = dx[keep], dy[keep]
    if dx.size == 0:
        return None
    ix = np.floor(dx / bin_px).astype(np.int64)
    iy = np.floor(dy / bin_px).astype(np.int64)
    # Vote in the cell AND its neighbours' union by using two half-shifted
    # grids would be overkill here: the frames are minutes apart and the
    # offset scatter is far below one cell.
    key = ix * 100003 + iy
    vals, counts = np.unique(key, return_counts=True)
    k = int(np.argmax(counts))
    if counts[k] < 4:
        return None
    sel = key == vals[k]
    return float(dx[sel].mean()), float(dy[sel].mean())


def match_stars(ref_xy: np.ndarray, xy: np.ndarray,
                shift: tuple[float, float],
                tol: float = MATCH_TOL_PX) -> np.ndarray:
    """Index into ``xy`` for each reference star, or -1 when unmatched.

    Nearest neighbour after applying ``shift``; a match is kept only when
    it is mutual (the reference star is also the nearest reference to its
    partner) and closer than ``tol`` pixels — so two stars can never claim
    one partner.
    """
    a = np.asarray(ref_xy, dtype=np.float64)
    b = np.asarray(xy, dtype=np.float64) + np.asarray(shift)[None, :]
    out = np.full(len(a), -1, dtype=np.int64)
    if len(a) == 0 or len(b) == 0:
        return out
    d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)
    nearest_b = d2.argmin(axis=1)
    nearest_a = d2.argmin(axis=0)
    for i, j in enumerate(nearest_b):
        if nearest_a[j] == i and d2[i, j] <= tol * tol:
            out[i] = j
    return out


def native_peak_factor(fwhm_binned_px: float, n_bin: int = 2,
                       beta: Optional[float] = None) -> float:
    """Brightest NATIVE pixel over the stored AVERAGE-binned peak pixel.

    An on-camera n x n average stores the mean of n^2 native pixels, so a
    star's brightest native pixel is always above the binned peak value —
    by more the sharper the star.  This returns the worst case over
    sub-pixel centring: the star centred ON a native pixel (which
    maximises the native peak) while that pixel shares its binned
    superpixel with three fainter neighbours.

    The profile is a circular Gaussian of the given binned-pixel FWHM (or
    a Moffat when ``beta`` is given), integrated over pixel areas on a
    fine grid.  The detector-engineer memo's "7-16% for FWHM 6-4 native
    pixels" is this function at 3 and 2 binned pixels.
    """
    fwhm_native = float(fwhm_binned_px) * n_bin
    over = 9                                   # sub-samples per native pixel
    half = int(math.ceil(3.0 * fwhm_native)) + n_bin
    n = (2 * half + 1) * over
    # Native-pixel grid with the star centred on the central native pixel.
    c = (np.arange(n) + 0.5) / over - (half + 0.5)
    xx, yy = np.meshgrid(c, c)
    r2 = xx * xx + yy * yy
    if beta is None:
        sigma = fwhm_native / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        prof = np.exp(-0.5 * r2 / sigma ** 2)
    else:
        alpha = fwhm_native / (2.0 * math.sqrt(2.0 ** (1.0 / beta) - 1.0))
        prof = (1.0 + r2 / alpha ** 2) ** (-beta)
    native = prof.reshape(2 * half + 1, over, 2 * half + 1, over).mean(
        axis=(1, 3))
    peak_native = native[half, half]
    # The central native pixel sits in one corner of its superpixel; by
    # symmetry every corner gives the same binned value.
    binned = native[half:half + n_bin, half:half + n_bin].mean()
    return float(peak_native / binned)
