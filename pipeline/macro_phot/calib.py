"""Local master-calibration application for eras with no reduced tree.

Roughly a third of the staged CV frames — ST LMi's whole 2024 season
(era 7), YZ Cnc's main run (era 7), and the small iKon blocks (eras 6, 47)
— have NO server-reduced counterpart anywhere in the archive.  The
provenance rule forbids photometering half a series from the reduced tree
and half from raw pixels, so those eras get ONE uniform local recipe
instead: subtract the era's master dark, divide by the era's master flat.

Why bother, when a local sky annulus already removes an additive pedestal?
Because a raw frame's problem is not its pedestal, it is its HOT PIXELS.
A measured pair (ST LMi, era 76) makes the point without argument: sep
finds 281 "detections" on the raw frame with a median FWHM of 2.4 px, and
161 on the dark-subtracted reduced frame with a median FWHM of 4.7 px.  The
120 extra raw detections are not stars — they are hot pixels sharp enough
to drag the frame's seeing estimate down by a factor of two, and every one
of them is a candidate comparison star the ensemble would have to reject.

This module is thin I/O only.  Which master serves which frame is decided
by the pure rules in :mod:`macro_phot.series` (:func:`pick_master`,
:func:`dark_exptime_matches`); this file reads the chosen files, caches
them, and does the arithmetic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

#: Flat pixels at or below this fraction of the flat's median are treated
#: as dead: dividing by them would manufacture enormous "flux" out of
#: nothing.  Their pixels are set to NaN, which sep's background estimator
#: and the aperture sum both propagate into a rejected measurement rather
#: than a spectacular fake star.
FLAT_MIN_RELATIVE = 0.2

#: Cache of master images already read this process, keyed by path.  Master
#: frames are 16-64 MB each and one era uses only a handful, so caching
#: turns a per-frame re-read into a per-era one; the build script's workers
#: each hold their own cache, which is exactly the intent.
_MASTER_CACHE: dict[str, np.ndarray] = {}


def read_master(path: Path) -> np.ndarray:
    """Read a master calibration frame as native float32 (cached).

    Masters live in the same fpack-tiled or plain-primary FITS shapes as
    science frames, so the extension choice mirrors
    :func:`macro_phot.extract.read_reduced`.
    """
    key = str(path)
    hit = _MASTER_CACHE.get(key)
    if hit is not None:
        return hit
    from astropy.io import fits
    with fits.open(path) as hdul:
        hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None else hdul[0]
        data = np.ascontiguousarray(hdu.data, dtype=np.float32)
    _MASTER_CACHE[key] = data
    return data


def normalized_flat(flat: np.ndarray) -> tuple[np.ndarray, float]:
    """Flat divided by its own median, with dead pixels masked to NaN.

    Returns ``(normalized, median)``.  Normalizing by the median keeps the
    frame's overall count level — and therefore its photon statistics and
    its saturation veto — unchanged; only the pixel-to-pixel throughput
    pattern is removed.
    """
    med = float(np.median(flat))
    if not np.isfinite(med) or med <= 0:
        raise ValueError("master flat has a non-positive median")
    norm = np.asarray(flat, dtype=np.float32) / np.float32(med)
    norm = np.where(norm <= FLAT_MIN_RELATIVE, np.nan, norm)
    return norm.astype(np.float32), med


def apply_masters(data: np.ndarray, dark: Optional[np.ndarray] = None,
                  flat: Optional[np.ndarray] = None
                  ) -> tuple[np.ndarray, str]:
    """Apply the local recipe: ``(data - dark) / normalized_flat``.

    Either master may be absent, and the recipe string returned alongside
    the pixels names exactly what was applied — ``'dark+flat'``,
    ``'dark_only'``, ``'flat_only'`` or ``'none'``.  That string is stored
    per frame in the products database, so a reader can tell at a glance
    whether a given light curve rests on a fully calibrated series or a
    partially calibrated one, without inferring it from era numbers.

    A master whose shape disagrees with the frame is REFUSED rather than
    broadcast or cropped: a mismatched master means the era assignment is
    wrong, and silently trimming one would calibrate the frame with the
    wrong pixels while reporting success.
    """
    out = np.asarray(data, dtype=np.float32)
    steps: list[str] = []
    if dark is not None:
        if dark.shape != out.shape:
            raise ValueError(
                f"master dark shape {dark.shape} != frame {out.shape}")
        out = out - dark
        steps.append("dark")
    if flat is not None:
        if flat.shape != out.shape:
            raise ValueError(
                f"master flat shape {flat.shape} != frame {out.shape}")
        norm, _med = normalized_flat(flat)
        out = out / norm
        steps.append("flat")
    recipe = "+".join(steps) if steps else "none"
    if recipe == "dark":
        recipe = "dark_only"
    elif recipe == "flat":
        recipe = "flat_only"
    return np.ascontiguousarray(out, dtype=np.float32), recipe
