"""Pure S2 saturation logic: is THIS star clipped in THIS frame?

WHY THIS MODULE EXISTS
----------------------
The ceiling stage answers "where does this readout mode's scale end".  The
committee's standing rule 4 asks a different, sharper question and asks it
per frame: *saturation is judged at the target, per frame, in native
pixels*.  A frame-level maximum cannot answer it (the brightest pixel of a
frame is usually a field star or a hot pixel, not the target), and a
binned peak cannot answer it either (an on-camera 2x2 AVERAGE hides a
saturated native pixel behind three unsaturated neighbours).

The observational-astronomer memo (OA.E3) showed what was at stake for
T CrB: opening four frames by hand found the R and I "anchors" flat-topped
at the clip.  This module makes that inspection a function:

* :func:`stamp_peak_stats` — given a small stamp around the target, the
  raw peak, how many pixels sit AT the peak (a flat top is the signature
  of clipping: a real stellar peak is one pixel, a clipped one is a
  plateau), and the star's width from second moments;
* :func:`native_equivalent_peak` — the binned peak scaled to the brightest
  native pixel it can hide (worst-case centring; the factor comes from
  :func:`rlmt_diagnostics.starphot.native_peak_factor`);
* :func:`saturation_verdict` — one of ``clipped`` / ``flat_topped`` /
  ``above_cap`` / ``clean`` from those numbers, with the rule written down
  in one place.

No I/O; the campaign script reads pixels and locates the target.
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

from .starphot import native_peak_factor

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Half-size of the stamp cut around the target position (pixels).
STAMP_HALF = 20

#: After the target is located (WCS or otherwise) the peak is searched
#: within this radius of the predicted position: WCS solutions are good to
#: a pixel, the fallback locator to a few.
RECENTRE_RADIUS_PX = 12

#: Pixels within this fraction of the stamp maximum count as "at the peak".
FLAT_TOP_FRACTION = 0.01

#: A plateau of at least this many pixels at the peak is a flat top.  A
#: critically sampled star has ONE pixel within 1% of its maximum (two
#: when it straddles a pixel boundary); four or more is a clipped core.
FLAT_TOP_MIN_PIXELS = 4

#: A peak within this fraction of the mode's measured clip IS the clip
#: (the High Gain clip is a per-pixel mound a few tens of ADU wide).
CLIP_PROXIMITY_FRACTION = 0.02


def stamp_peak_stats(stamp: np.ndarray, sky: Optional[float] = None) -> dict:
    """Peak, plateau size and width of the star in a stamp.

    Parameters
    ----------
    stamp
        Raw pixels around the target (ADU), target near the centre.
    sky
        Background level to subtract for the width; the stamp's border
        median is used when None.

    Returns
    -------
    dict with ``peak_raw`` (maximum raw pixel), ``peak_y``/``peak_x``
    (its position in the stamp), ``n_at_peak`` (pixels within
    :data:`FLAT_TOP_FRACTION` of the sky-subtracted peak, counted in the
    7x7 box around the peak), ``sky`` and ``fwhm_px`` (Gaussian-equivalent
    FWHM from the sky-subtracted second moments within the half-maximum
    footprint grown by two pixels; NaN when the moments are degenerate).
    """
    a = np.asarray(stamp, dtype=np.float64)
    ny, nx = a.shape
    if sky is None:
        border = np.concatenate([a[0, :], a[-1, :], a[:, 0], a[:, -1]])
        sky = float(np.median(border))
    py, px = np.unravel_index(int(np.argmax(a)), a.shape)
    peak = float(a[py, px])
    amp = peak - sky
    box = a[max(py - 3, 0):py + 4, max(px - 3, 0):px + 4]
    n_at_peak = int(((box - sky) >= (1.0 - FLAT_TOP_FRACTION) * amp).sum()) \
        if amp > 0 else 0
    fwhm = float("nan")
    if amp > 0:
        sub = a - sky
        yy, xx = np.mgrid[0:ny, 0:nx]
        r2 = (yy - py) ** 2 + (xx - px) ** 2
        # Footprint: pixels above 10% of the peak within 3 half-max radii.
        hm = sub >= 0.5 * amp
        r_hm = math.sqrt(max(float(hm.sum()), 1.0) / math.pi)
        foot = (r2 <= (3.0 * r_hm + 2.0) ** 2) & (sub > 0.05 * amp)
        w = np.where(foot, sub, 0.0)
        tot = w.sum()
        if tot > 0:
            cy, cx = (w * yy).sum() / tot, (w * xx).sum() / tot
            var = ((w * ((yy - cy) ** 2 + (xx - cx) ** 2)).sum() / tot) / 2.0
            if var > 0:
                fwhm = 2.0 * math.sqrt(2.0 * math.log(2.0) * var)
    return {"peak_raw": peak, "peak_y": int(py), "peak_x": int(px),
            "n_at_peak": n_at_peak, "sky": float(sky), "fwhm_px": fwhm}


def native_equivalent_peak(peak_raw: float, bias: float, fwhm_px: float,
                           n_bin: int) -> tuple[float, float]:
    """Brightest native pixel an average-binned peak can hide.

    Returns ``(native_peak_adu, factor)``.  For an unbinned frame
    (``n_bin`` 1) the factor is exactly 1.  For an n x n on-camera AVERAGE
    the signal above bias is multiplied by
    :func:`~rlmt_diagnostics.starphot.native_peak_factor` evaluated at the
    star's measured width — the worst case over sub-pixel centring, which
    is the right case for a veto: the question is whether a native pixel
    COULD have clipped.  A missing width falls back to 2 binned pixels
    (sharp seeing; conservative).
    """
    if n_bin <= 1:
        return float(peak_raw), 1.0
    w = fwhm_px if (fwhm_px is not None and np.isfinite(fwhm_px)
                    and fwhm_px > 0.5) else 2.0
    factor = native_peak_factor(w, n_bin=n_bin)
    return float(bias + (peak_raw - bias) * factor), float(factor)


def saturation_verdict(peak_native: float, bias: float, clip: float,
                       cap_fraction: Optional[float],
                       n_at_peak: int) -> str:
    """One verdict per (target, frame).

    In this order:

    * ``clipped``     — the native-equivalent peak is within
      :data:`CLIP_PROXIMITY_FRACTION` of the mode's clip, or above it;
    * ``flat_topped`` — below that, but the core is a plateau of at least
      :data:`FLAT_TOP_MIN_PIXELS` pixels (clipping caught in the act: the
      per-pixel clip on the GSENSE4040 varies by tens of ADU, so a star
      can be clipped a little below the mode's nominal value);
    * ``above_cap``   — unclipped, but above the mode's linearity cap
      (``cap_fraction`` of the usable scale): usable for nothing that
      needs photometric accuracy;
    * ``clean``       — below the cap.

    With no cap for the mode (``cap_fraction`` None) the third verdict is
    never issued and ``clean`` means only "not clipped".
    """
    span = clip - bias
    frac = (peak_native - bias) / span if span > 0 else float("nan")
    if frac >= 1.0 - CLIP_PROXIMITY_FRACTION:
        return "clipped"
    if n_at_peak >= FLAT_TOP_MIN_PIXELS:
        return "flat_topped"
    if cap_fraction is not None and frac > cap_fraction:
        return "above_cap"
    return "clean"


def brightest_central_source(x: np.ndarray, y: np.ndarray, flux: np.ndarray,
                             shape: tuple[int, int],
                             max_offset_fraction: float = 0.25,
                             ) -> Optional[int]:
    """Fallback target locator: the brightest source near the frame centre.

    For a frame with no usable WCS.  A monitoring target is, by
    construction of the observation, near the centre; the brightest
    detection within ``max_offset_fraction`` of the frame's short side of
    the centre is taken.  Returns its index, or None when nothing is
    there.  The campaign records which locator was used for every frame,
    so a reader can separate WCS-located verdicts from these.
    """
    if len(x) == 0:
        return None
    ny, nx = shape
    r = np.hypot(np.asarray(x) - nx / 2.0, np.asarray(y) - ny / 2.0)
    near = r <= max_offset_fraction * min(ny, nx)
    if not near.any():
        return None
    idx = np.flatnonzero(near)
    return int(idx[np.argmax(np.asarray(flux)[idx])])
