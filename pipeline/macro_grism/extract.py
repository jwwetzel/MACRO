"""Horne-style optimal extraction with flanking-band background.

Background policy (the paper's calibration-debt argument lives here):

* The FLANKING BANDS — two strips above and below the trace window —
  measure sky + scattered-light halo + smooth dark current locally, per
  column, and a straight line through the two band medians is subtracted
  across the window.  On a raw frame this doubles as local dark removal:
  any dark structure smooth over ~60 px vanishes with the halo.
* What flanking bands CANNOT remove is pixel-scale dark structure (hot
  pixels inside the aperture).  The recovered era-76 Mode0 master darks
  can.  Every frame is therefore extracted TWICE — once from the raw
  pixels, once from master-dark-subtracted pixels, flanking bands both
  times — and the difference between the two spectra is the calibration
  debt the paper must quote as an error term.

The Horne weighting follows Horne (1986): with a normalized spatial
profile P and per-pixel variance V, the optimal flux estimate per column
is sum(P*(D-B)/V) / sum(P^2/V) — inverse-variance weighting that
down-weights the noisy profile wings (and, with the mask, ignores
saturated pixels instead of eating them).
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .trace import CENTROID_HALFWIN

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report quotes these).
# --------------------------------------------------------------------------

#: Extraction aperture half-width (px).  The trace FWHM measures ~5-9 px
#: on the validation frames; +/-12 px holds > 99% of a 9 px FWHM Gaussian.
APERTURE_HALFWIN = CENTROID_HALFWIN

#: Flanking bands: inner and outer edge offsets from the trace center
#: (px), used on BOTH sides.  Inner edge 30 px keeps the bands out of the
#: trace wings; outer edge 60 px keeps them inside the same patch of halo
#: (the halo varies over ~hundreds of px).
BAND_INNER = 30
BAND_OUTER = 60

#: Spatial profile: median-combined over this many column chunks, then
#: smoothed along the dispersion is NOT done — the profile is one global
#: shape (the traces show no measurable profile change along x at the
#: validation SNR, and a global profile cannot chase noise).
PROFILE_CHUNKS = 16

#: Read noise in electrons for Mode0 (ASI/IMX455 class, gain setting 100).
#: The S2 detector campaign bracketed the CMOS families; Mode0's
#: conversion gain is 0.2467 e-/ADU (header EGAIN) and the read noise a
#: few electrons — 3.5 e- is the conservative bracket top.  The Horne
#: weights change negligibly across 1.5-3.5 e-, and the variance array is
#: labelled with the value used.
READ_NOISE_E = 3.5

#: Fallback conversion gain (e-/ADU) when the header lacks EGAIN.
DEFAULT_EGAIN = 0.2467

#: Saturation threshold (ADU).  Mode0 raws rail at ~16.4k ADU (12-bit ADC
#: scaled x4); pixels at or above this are masked out of the Horne sum and
#: counted per column.
SATURATION_ADU = 16300.0


def flanking_background(column: np.ndarray, yc: float,
                        halfwin: int = APERTURE_HALFWIN,
                        inner: int = BAND_INNER,
                        outer: int = BAND_OUTER):
    """Per-column background from the two flanking bands.

    Returns ``(background_at_window_rows, ok)`` where the background is a
    straight line through (band centers, band medians) evaluated at the
    aperture rows — a LINEAR local model, so a smooth vertical gradient in
    halo/sky/dark is removed exactly (unit-tested).  ``ok`` is False when
    either band leaves the frame; the caller records the column as
    unextractable rather than inventing a background.
    """
    ny = len(column)
    lo, hi = int(round(yc)) - halfwin, int(round(yc)) + halfwin + 1
    b_lo0, b_lo1 = int(round(yc)) - outer, int(round(yc)) - inner
    b_hi0, b_hi1 = int(round(yc)) + inner, int(round(yc)) + outer
    if b_lo0 < 0 or b_hi1 + 1 > ny:
        return None, False
    below = np.median(column[b_lo0:b_lo1])
    above = np.median(column[b_hi0:b_hi1 + 1])
    y_below = (b_lo0 + b_lo1 - 1) / 2.0          # band center rows
    y_above = (b_hi0 + b_hi1) / 2.0
    slope = (above - below) / (y_above - y_below)
    rows = np.arange(lo, hi)
    return below + slope * (rows - y_below), True


def build_profile(cutouts: np.ndarray) -> np.ndarray:
    """Normalized spatial profile from background-subtracted, per-column
    flux-normalized cutouts (shape: n_columns x window).

    Median over columns (cosmic rays lose), floored at zero (a profile is
    a probability density; negative wings are noise), normalized to unit
    sum.  Falls back to a flat profile if everything medians to zero —
    the extraction then degrades to a box sum instead of dividing by 0.
    """
    prof = np.median(cutouts, axis=0)
    prof = np.clip(prof, 0.0, None)
    total = prof.sum()
    if total <= 0:
        return np.full(cutouts.shape[1], 1.0 / cutouts.shape[1])
    return prof / total


def horne_column(window: np.ndarray, background: np.ndarray,
                 profile: np.ndarray, egain: float, read_noise_e: float,
                 sat_adu: float = SATURATION_ADU):
    """Optimal flux for ONE column window (Horne 1986 eq. 8, one pass).

    Variance model per pixel, in ADU^2:  V = (RN/g)^2 + max(D, 0)/g
    with g = e-/ADU — read noise plus photon noise of the TOTAL counts
    (source + background; the background's photons are noise too).
    Saturated pixels are masked.  Returns (flux, variance, n_saturated);
    flux is None when every pixel of the window is masked.
    """
    sat = window >= sat_adu
    good = ~sat
    if not good.any():
        return None, None, int(sat.sum())
    d = window - background                        # net counts, ADU
    var = (read_noise_e / egain) ** 2 + np.clip(window, 0, None) / egain
    p, v = profile[good], var[good]
    denom = float((p * p / v).sum())
    if denom <= 0:
        return None, None, int(sat.sum())
    flux = float((p * d[good] / v).sum() / denom)
    return flux, 1.0 / denom, int(sat.sum())


def extract_spectrum(data: np.ndarray, trace_coeffs: np.ndarray,
                     egain: float = DEFAULT_EGAIN,
                     read_noise_e: float = READ_NOISE_E,
                     dark: Optional[np.ndarray] = None,
                     halfwin: int = APERTURE_HALFWIN) -> dict:
    """Extract the full spectrum along a fitted trace.

    ``dark`` (a master dark, same shape) is subtracted pixel-wise first
    when given — the comparison arm.  Flanking-band background is applied
    in BOTH modes (see the module header for why).

    Returns a dict of 1-D arrays over columns: ``flux`` (optimal, ADU),
    ``var``, ``box`` (plain aperture sum, the cross-check), ``bg``
    (background level at trace center), ``n_sat``, plus scalars
    ``profile`` and ``n_extracted``.  Columns that cannot be extracted
    (windows off frame) carry NaN.
    """
    work = data if dark is None else data - dark
    ny, nx = work.shape
    yc = np.polyval(trace_coeffs, np.arange(nx))
    flux = np.full(nx, np.nan)
    var = np.full(nx, np.nan)
    box = np.full(nx, np.nan)
    bg = np.full(nx, np.nan)
    n_sat = np.zeros(nx, dtype=np.int32)

    # ---- pass 1: background-subtracted cutouts for the spatial profile.
    # Sampled on a chunk grid; each cutout is normalized by its own sum so
    # bright and faint columns vote equally on the SHAPE.
    cut = []
    for x in range(0, nx, max(1, nx // (PROFILE_CHUNKS * 8))):
        b, ok = flanking_background(work[:, x], yc[x], halfwin)
        if not ok:
            continue
        lo = int(round(yc[x])) - halfwin
        w = work[lo:lo + 2 * halfwin + 1, x] - b
        s = w.sum()
        if s > 0:
            cut.append(w / s)
    profile = (build_profile(np.array(cut)) if cut
               else np.full(2 * halfwin + 1, 1.0 / (2 * halfwin + 1)))

    # ---- pass 2: per-column optimal extraction with that profile.
    for x in range(nx):
        b, ok = flanking_background(work[:, x], yc[x], halfwin)
        if not ok:
            continue
        lo = int(round(yc[x])) - halfwin
        w = work[lo:lo + 2 * halfwin + 1, x]
        f, v, ns = horne_column(w, b, profile, egain, read_noise_e)
        n_sat[x] = ns
        if f is None:
            continue
        flux[x] = f
        var[x] = v
        box[x] = float((w - b).sum())
        bg[x] = float(b[halfwin])
    return {"flux": flux, "var": var, "box": box, "bg": bg,
            "n_sat": n_sat, "profile": profile,
            "n_extracted": int(np.isfinite(flux).sum())}


def median_relative_difference(flux_a: np.ndarray,
                               flux_b: np.ndarray) -> Optional[float]:
    """The calibration-debt statistic: median over columns of
    |A - B| / |A|, restricted to columns where both spectra exist and A
    has meaningful signal (|A| above the 25th percentile of |A| — ratio
    against near-zero flux measures nothing but noise).

    A = flanking-only, B = master-dark-subtracted; the paper quotes this
    number as the error bar for season-1's calibration debt.
    """
    both = np.isfinite(flux_a) & np.isfinite(flux_b)
    if both.sum() < 100:
        return None
    a, b = flux_a[both], flux_b[both]
    floor = np.percentile(np.abs(a), 25)
    sig = np.abs(a) >= max(floor, 1e-9)
    if sig.sum() < 50:
        return None
    return float(np.median(np.abs(a[sig] - b[sig]) / np.abs(a[sig])))
