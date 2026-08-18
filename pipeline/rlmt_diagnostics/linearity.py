"""Pure S2 linearity logic: counts-vs-exptime from archival exposure ladders.

THE METHOD
----------
A linear detector doubles its counts when the exposure doubles.  An
exposure *ladder* — the same bright target shot at several exposure times
in one visit (the 2024-05-20 Vega BeStar sequence: 0.0001 s -> 0.1 s) —
turns that statement into a measurement: fit  flux = k * exptime  through
the rung medians and read the per-rung residual.  Rungs that fall LOW at
the top of the ladder are the linearity roll-off (or the hard ceiling);
scatter at the bottom is shutter/exposure-timing error.  No ladder for a
mode is itself a finding and is recorded as one.

Flux extraction must survive a grism (the Vega ladder is spectra, not
points) and needs no WCS: the measure is the background-subtracted sum
over the single brightest ``box`` x ``box`` window in the frame, located by
a coarse box-sum scan.  Pure numpy (integral images), no photometry
package, deterministic.
"""

from __future__ import annotations

import math
import re
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Photometry window (pixels).  Big enough to hold a saturated star's whole
#: profile or a grism trace segment; small enough that the background
#: estimate stays local.
LADDER_BOX = 96

#: Coarse scan stride when locating the brightest window (the exact peak is
#: then refined by centering the box on the brightest coarse cell).
LADDER_SCAN_STRIDE = 32

#: A ladder needs at least this many distinct exposure-time rungs.
MIN_RUNGS = 3

#: ... and at least this many frames per rung to median away transients.
MIN_FRAMES_PER_RUNG = 2

#: Filename exposure-time token: '0p0001s' -> 0.0001, '10s' -> 10.0.  The
#: header EXPTIME rounds sub-millisecond exposures to 0.0 (observed on the
#: Vega ladder), so the filename is the better record at the short end.
_EXPTIME_TOKEN_RE = re.compile(r"(?:^|_)(\d+(?:p\d+)?)s(?:_|$|\.)")


def parse_exptime_token(basename: str) -> Optional[float]:
    """Exposure time encoded in a filename, or None.

    ``kaf_Vega_0p0001s_lrg_0.fts.fz`` -> 0.0001;
    ``mjc_HD_20134_hrg_64s_...`` -> 64.0.  The 'p' is the decimal point
    (filenames cannot carry '.').  Returns None when no token matches —
    callers then fall back to the header EXPTIME.
    """
    m = _EXPTIME_TOKEN_RE.search(basename)
    if not m:
        return None
    return float(m.group(1).replace("p", "."))


def effective_exptime(header_exptime: Optional[float],
                      basename: str) -> Optional[float]:
    """The exposure time to trust for ladder work.

    The filename token wins whenever the header is missing, non-positive,
    or disagrees with the token by more than 20% (the header driver rounds
    0.0001 s to 0.0 — the observed failure).  Otherwise the header stands.
    """
    tok = parse_exptime_token(basename)
    h = header_exptime
    if h is None or not math.isfinite(h) or h <= 0:
        return tok if tok is not None else h
    if tok is not None and tok > 0 and abs(h - tok) > 0.2 * max(h, tok):
        return tok
    return float(h)


def brightest_box_flux(img: np.ndarray, box: int = LADDER_BOX,
                       stride: int = LADDER_SCAN_STRIDE) -> dict:
    """Background-subtracted flux in the brightest box-window of a frame.

    Integral-image box sums on a stride grid find the brightest window;
    the background is the frame median (dominant sky/bias level) times the
    box area.  Returns the flux, the window's corner, the frame's peak
    pixel inside the window, and the frame median — everything the ladder
    fit and the saturation cross-check need.
    """
    a = np.asarray(img, dtype=np.float64)
    ny, nx = a.shape
    b = min(box, ny, nx)
    # Integral image: S[i, j] = sum of a[:i, :j]; box sums in O(1) each.
    S = np.zeros((ny + 1, nx + 1))
    np.cumsum(np.cumsum(a, axis=0), axis=1, out=S[1:, 1:])
    ys = np.arange(0, ny - b + 1, stride)
    xs = np.arange(0, nx - b + 1, stride)
    sums = (S[np.ix_(ys + b, xs + b)] - S[np.ix_(ys, xs + b)]
            - S[np.ix_(ys + b, xs)] + S[np.ix_(ys, xs)])
    iy, ix = np.unravel_index(int(np.argmax(sums)), sums.shape)
    y0, x0 = int(ys[iy]), int(xs[ix])
    med = float(np.median(a))
    window = a[y0:y0 + b, x0:x0 + b]
    return {
        "flux": float(window.sum() - med * b * b),
        "y0": y0, "x0": x0, "box": b,
        "peak_adu": float(window.max()),
        "sky_med": med,
    }


def fit_ladder(exptimes: Sequence[float], fluxes: Sequence[float]) -> Optional[dict]:
    """Fit  flux = k * exptime  and report per-rung residuals.

    ``exptimes``/``fluxes`` are ONE value per rung (the caller medians its
    frames first).  The rate k is the median of flux/exptime over rungs —
    robust, and immune to a single rolled-off top rung dragging the line.
    Residuals are percentages: 100 * (flux / (k * t) - 1).  Returns None
    for fewer than :data:`MIN_RUNGS` usable rungs.
    """
    t = np.asarray(exptimes, dtype=np.float64)
    f = np.asarray(fluxes, dtype=np.float64)
    ok = np.isfinite(t) & np.isfinite(f) & (t > 0)
    t, f = t[ok], f[ok]
    if t.size < MIN_RUNGS:
        return None
    order = np.argsort(t)
    t, f = t[order], f[order]
    rates = f / t
    k = float(np.median(rates))
    if k <= 0:
        return None
    resid_pct = 100.0 * (f / (k * t) - 1.0)
    return {
        "rate_adu_per_s": k,
        "exptimes": t.tolist(),
        "fluxes": f.tolist(),
        "resid_pct": resid_pct.tolist(),
        "max_abs_resid_pct": float(np.max(np.abs(resid_pct))),
        "n_rungs": int(t.size),
    }


def group_ladders(rows: Sequence[tuple],
                  min_rungs: int = MIN_RUNGS,
                  min_frames: int = MIN_FRAMES_PER_RUNG) -> list[dict]:
    """Find exposure ladders in manifest rows.

    ``rows`` = (night, target_key, readoutm, exptime_bin, n_frames) tuples
    (one per group, from SQL).  A ladder = one (night, target, mode) with
    >= ``min_rungs`` distinct exposure bins each holding >= ``min_frames``
    frames.  Returns one dict per ladder, rungs sorted by exposure —
    ready for the campaign script to fetch pixels for.
    """
    from collections import defaultdict
    sets: dict[tuple, dict[float, int]] = defaultdict(dict)
    for night, target, mode, ebin, n in rows:
        if ebin is None or ebin <= 0 or target is None:
            continue
        sets[(night, target, mode)][float(ebin)] = int(n)
    out = []
    for (night, target, mode), rungs in sets.items():
        good = {t: n for t, n in rungs.items() if n >= min_frames}
        if len(good) >= min_rungs:
            out.append({"night": night, "target_key": target, "mode": mode,
                        "rungs": sorted(good), "n_frames": sum(good.values())})
    # Deterministic order: richest ladders first, then by night/target.
    out.sort(key=lambda d: (-len(d["rungs"]), -d["n_frames"],
                            d["night"] or "", d["target_key"]))
    return out
