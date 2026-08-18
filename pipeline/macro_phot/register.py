"""Sky-chained frame registration: S1 plate solutions as a matching route.

The S4 prototype registered every frame with astroalign, because the two
prototype polars were essentially unsolved.  Over the full CV set roughly
half the staged frames now carry an S1 plate solution written as a FITS
header sidecar under ``products/astrom/wcs/``, and for those frames there
is a shorter, sturdier road from a detection to a star identity:

    frame pixels --(frame WCS)--> sky --(reference WCS)--> reference pixels

Then the SAME pure one-to-one nearest-neighbour assignment
(:func:`macro_phot.photometry.match_one_to_one`) that follows an astroalign
transform finishes the job.  The chain costs two matrix evaluations instead
of a RANSAC triangle search, and it cannot starve on a sparse field.

What it can do is be quietly WRONG — a solution fitted on a different pixel
grid, or SIP terms extrapolated past the frame edge — which is why the
build script never trusts a chained match on faith: it measures the match
against :func:`macro_phot.series.wcs_match_ok` and falls back to the
astroalign ladder when the chain does not deliver.

One geometric fact makes the chain legal at all, and it was verified rather
than assumed: the observatory's reduction preserves the pixel grid.  On a
test pair (ST LMi, era 76) the eight brightest stars land within 0.02 pixels
of each other in the raw and reduced frames, so a WCS solved on the RAW
frame describes the REDUCED frame too.  ``pixel_grids_agree`` below is the
reusable form of that check.

Only :func:`load_wcs` touches the filesystem; everything else is arithmetic
on arrays and is unit-tested against hand-built WCS objects.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

#: FITS extension the S1 batch writes its solutions with.
WCS_SUFFIX = ".wcs"

#: Compressed-FITS suffix stripped before the sidecar name is formed: S1
#: names its sidecar after the UNCOMPRESSED file name.
FPACK_SUFFIX = ".fz"


def sidecar_path(wcs_root: Path, frame_relpath: str) -> Path:
    """Where S1 would have written this frame's plate solution.

    ``frame_relpath`` is the archive-relative path recorded in the manifest
    (e.g. ``rawimage/2025-01-09/maw_ST_LMi_i_120s_...fts.fz``).  S1 names
    the sidecar after the file WITHOUT the fpack suffix, so the mapping is
    'drop a trailing .fz, append .wcs' and the directory tree is mirrored.
    Existence is the caller's question — this function only names the file.
    """
    rel = str(frame_relpath)
    if rel.endswith(FPACK_SUFFIX):
        rel = rel[: -len(FPACK_SUFFIX)]
    return Path(wcs_root) / (rel + WCS_SUFFIX)


def load_wcs(path: Path):
    """Read one S1 sidecar into an astropy WCS, or None if unusable.

    Returns None — never raises — for a missing, truncated or
    WCS-less sidecar: a frame without astrometry is an ordinary, expected
    state that must route to astroalign, not an error that stops a batch.
    """
    try:
        from astropy.io import fits
        from astropy.wcs import WCS
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            if "CTYPE1" not in hdr:
                return None
            return WCS(hdr, relax=True)
    except Exception:
        return None


def chain_to_reference(frame_wcs, ref_wcs, xy: np.ndarray) -> np.ndarray:
    """Carry frame-pixel positions into REFERENCE-frame pixels via the sky.

    Parameters
    ----------
    frame_wcs, ref_wcs
        astropy WCS objects for this frame and for the series reference.
    xy
        (N, 2) detection positions in this frame's pixels, in sep's
        0-indexed convention.

    Returns
    -------
    np.ndarray
        (N, 2) positions in the reference frame's 0-indexed pixels, NaN
        wherever the round trip failed (a detection whose sky position
        falls outside the reference solution's valid domain).

    The 0 passed to ``all_pix2world`` / ``all_world2pix`` is the origin
    convention: sep centroids and this pipeline are 0-indexed throughout,
    while FITS itself counts from 1.  Getting that wrong shifts every star
    by exactly one pixel — a shift small enough to survive a generous match
    tolerance and large enough to poison every centroid downstream, so it
    is stated here rather than left to the reader.
    """
    xy = np.asarray(xy, dtype=float)
    if xy.size == 0:
        return np.empty((0, 2), dtype=float)
    try:
        sky = frame_wcs.all_pix2world(xy[:, 0], xy[:, 1], 0)
        out_x, out_y = ref_wcs.all_world2pix(sky[0], sky[1], 0)
    except Exception:
        return np.full(xy.shape, np.nan)
    out = np.column_stack([np.asarray(out_x, dtype=float),
                           np.asarray(out_y, dtype=float)])
    out[~np.isfinite(out).all(axis=1)] = np.nan
    return out


def pixel_grids_agree(xy_a: np.ndarray, xy_b: np.ndarray,
                      tol_px: float = 0.5) -> tuple[bool, float]:
    """Do two detection lists describe the SAME pixel grid?

    Pure geometry, used to justify applying a raw-frame plate solution to
    the reduced version of the same exposure.  Each position in ``xy_a``
    is paired with its nearest neighbour in ``xy_b``; the function returns
    ``(agree, median_separation_px)`` where agreement means the median
    separation is within ``tol_px``.  A reduction that cropped, flipped or
    rebinned the frame produces a median far above any sane tolerance, so
    the test fails loudly instead of subtly.
    """
    a = np.asarray(xy_a, dtype=float)
    b = np.asarray(xy_b, dtype=float)
    if a.size == 0 or b.size == 0:
        return False, float("nan")
    d = np.hypot(a[:, 0, None] - b[None, :, 0],
                 a[:, 1, None] - b[None, :, 1])
    nearest = d.min(axis=1)
    med = float(np.median(nearest))
    return bool(med <= tol_px), med


def sky_of_pixel(wcs, x: float, y: float) -> Optional[tuple[float, float]]:
    """One pixel's ICRS (ra, dec) in degrees, or None if the WCS refuses."""
    try:
        ra, dec = wcs.all_pix2world([float(x)], [float(y)], 0)
        return float(ra[0]), float(dec[0])
    except Exception:
        return None


def pixel_of_sky(wcs, ra_deg: float, dec_deg: float
                 ) -> Optional[tuple[float, float]]:
    """Where a sky position lands in this frame's 0-indexed pixels.

    The build script uses this to find THE TARGET in a plate-solved
    reference frame: the target's catalogue coordinates go in, a pixel
    position comes out, and the nearest reference star to that pixel is
    the target — no Gaia query, no triangle fit, no network.  Returns None
    when the WCS cannot invert (position off the projection's valid sky).
    """
    try:
        x, y = wcs.all_world2pix([float(ra_deg)], [float(dec_deg)], 0)
        if not (np.isfinite(x[0]) and np.isfinite(y[0])):
            return None
        return float(x[0]), float(y[0])
    except Exception:
        return None
