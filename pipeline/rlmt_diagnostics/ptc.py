"""Pure S2 photon-transfer logic: gain, read noise, StackPro signature.

THE METHOD (classic difference-pair photon transfer, e.g. Janesick)
-------------------------------------------------------------------
Take two frames of the SAME scene under the SAME settings, seconds apart.
Their *difference* cancels everything fixed (scene, dark structure, flat
field, amp glow) and keeps only the random noise, doubled:

    var(a - b) / 2  =  read_noise^2  +  level / gain        [ADU^2]

so a straight line of half-difference-variance against signal level has
slope 1/gain (e-/ADU falls straight out) and intercept read_noise^2.
Repeated darks give the level~0 anchor (read noise); repeated star-field
exposures (the 2023-06-07 Albireo series) sweep the level axis through
their pixel-brightness range, so ONE pair yields many (level, variance)
points by binning pixels on their mean level.

StackPro's fingerprint (the model this campaign TESTED and adopted): a
StackPro frame is the SUM of N_sub sub-exposures.  A sum multiplies the
bias offset by N_sub, the read-noise VARIANCE by N_sub, and the saturation
ceiling by N_sub, while leaving the PTC's Poisson slope UNCHANGED (signal
and signal-variance scale together).  :func:`stackpro_signature` reads
N_sub from those three ratios at once — the measured archive gives 16 on
all three.  (The rival AVERAGING architecture would instead divide every
variance by N_sub and raise the apparent PTC gain by N_sub;
:func:`nsub_estimate` covers that hypothesis and the data rejected it —
the dark-floor variance RATIO came out ~16, not ~1/16.)

All functions are pure numpy; the campaign script owns file I/O and the
pairing of frames, the tests own synthetic Poisson frames with known truth.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Robust sigma-clip threshold for pixels inside one level bin: hot pixels,
#: cosmic rays and stars-on-the-move violate the same-scene assumption and
#: must not inflate the variance.  4 sigma on a Gaussian keeps 99.994% of
#: honest pixels.
CLIP_SIGMA = 4.0

#: Clipping iterations — two passes are enough for the contamination levels
#: seen here (a first pass to find the scale, a second to settle).
CLIP_ITERS = 2

#: Minimum pixels per level bin for a usable (level, variance) point.
MIN_PIXELS_PER_BIN = 2000

#: Pixels above this fraction of the mode's ceiling are excluded from PTC
#: statistics: at the clip the variance collapses artificially (every value
#: is the same code) and would drag the fitted slope down.
PTC_LEVEL_CEILING_FRACTION = 0.85

#: MAD -> sigma conversion for a Gaussian (1 / Phi^-1(3/4)).
MAD_TO_SIGMA = 1.4826

#: Level cap for the LOWER-bound gain fit on star-field (light) pairs.
#: Consecutive star exposures move by a fraction of a pixel between frames
#: (tracking/seeing), which inflates the difference variance wherever the
#: image has gradients — i.e. everywhere bright.  Below this level the
#: pixels are mostly flat sky, where motion adds least, so the fitted
#: slope there gives the LEAST-inflated variance -> a defensible lower
#: bound on the gain (measured var >= true var  =>  1/slope <= gain).
GAIN_LOWER_BOUND_LEVEL_ADU = 300.0


def robust_pair_variance(diff: np.ndarray) -> tuple[float, int]:
    """Half the variance of a difference sample, sigma-clipped.

    Returns ``(var(diff)/2, n_used)``.  The clip is MAD-based (median
    absolute deviation), so a few cosmic-ray or hot-pixel outliers cannot
    seed the estimate they are being tested against.
    """
    d = np.asarray(diff, dtype=np.float64).ravel()
    mask = np.isfinite(d)
    for _ in range(CLIP_ITERS):
        sel = d[mask]
        if sel.size < 2:
            return float("nan"), 0
        med = np.median(sel)
        sigma = MAD_TO_SIGMA * np.median(np.abs(sel - med))
        if sigma == 0:
            break
        mask &= np.abs(d - med) <= CLIP_SIGMA * sigma
    sel = d[mask]
    if sel.size < 2:
        return float("nan"), 0
    return float(np.var(sel, ddof=1) / 2.0), int(sel.size)


def pair_ptc_points(a: np.ndarray, b: np.ndarray,
                    n_bins: int = 12,
                    level_max: Optional[float] = None,
                    ) -> list[dict]:
    """Turn one same-scene frame pair into (level, variance) PTC points.

    Pixels are binned on their mean level ``(a+b)/2`` (log-spaced bins from
    the sky floor up to ``level_max``, default the 99.9th percentile), and
    each populated bin contributes one point: median level, robust
    half-difference variance, pixel count.  Bins with fewer than
    :data:`MIN_PIXELS_PER_BIN` pixels are dropped (their variance estimate
    is noise, not measurement).
    """
    af = np.asarray(a, dtype=np.float64).ravel()
    bf = np.asarray(b, dtype=np.float64).ravel()
    mean = (af + bf) / 2.0
    diff = af - bf
    if level_max is None:
        level_max = float(np.percentile(mean, 99.9))
    lo = float(np.percentile(mean, 1.0))
    if not (level_max > lo > 0):
        # Degenerate frame (constant, or non-positive floor): single bin.
        var, n = robust_pair_variance(diff)
        return ([{"level": float(np.median(mean)), "var": var, "n_pix": n}]
                if n >= MIN_PIXELS_PER_BIN else [])
    edges = np.geomspace(lo, level_max, n_bins + 1)
    out: list[dict] = []
    idx = np.digitize(mean, edges)
    for i in range(1, n_bins + 1):
        sel = idx == i
        if int(sel.sum()) < MIN_PIXELS_PER_BIN:
            continue
        var, n = robust_pair_variance(diff[sel])
        if n < MIN_PIXELS_PER_BIN:
            continue
        out.append({"level": float(np.median(mean[sel])),
                    "var": var, "n_pix": n})
    return out


def fit_ptc(levels: Sequence[float], variances: Sequence[float],
            n_pix: Optional[Sequence[float]] = None) -> Optional[dict]:
    """Weighted straight-line fit  var = rn2 + level/gain.

    Weights are sqrt(n_pix) per point when given (a variance estimated from
    more pixels is better known).  Returns gain (e-/ADU), read noise in ADU
    and electrons, the raw slope/intercept, and simple fit uncertainties
    from the weighted normal equations; ``None`` if fewer than 3 points or
    a non-positive slope (no Poisson signal — darks alone can do this, and
    the caller must then report read noise only).
    """
    x = np.asarray(levels, dtype=np.float64)
    y = np.asarray(variances, dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 3:
        return None
    w = (np.sqrt(np.asarray(n_pix, dtype=np.float64)[ok])
         if n_pix is not None else np.ones_like(x))
    # Weighted least squares on [1, x]: solve the 2x2 normal equations.
    W = w / w.sum()
    xm = float((W * x).sum())
    ym = float((W * y).sum())
    sxx = float((W * (x - xm) ** 2).sum())
    if sxx <= 0:
        return None
    sxy = float((W * (x - xm) * (y - ym)).sum())
    slope = sxy / sxx
    intercept = ym - slope * xm
    if slope <= 0:
        return None
    # Scatter-based uncertainties (assumes the line is the right model).
    resid = y - (intercept + slope * x)
    dof = max(x.size - 2, 1)
    s2 = float((W * resid ** 2).sum()) * x.size / dof
    slope_err = np.sqrt(s2 / (sxx * x.size))
    inter_err = np.sqrt(s2 * (1.0 / x.size + xm ** 2 / (sxx * x.size)))
    gain = 1.0 / slope
    rn_adu = float(np.sqrt(max(intercept, 0.0)))
    return {
        "gain_e_per_adu": float(gain),
        "gain_err": float(slope_err / slope ** 2),   # d(1/s) = ds/s^2
        "read_noise_adu": rn_adu,
        "read_noise_adu_err": (float(inter_err / (2 * rn_adu))
                               if rn_adu > 0 else float(inter_err ** 0.5)),
        "read_noise_e": float(rn_adu * gain),
        "slope": float(slope), "intercept": float(intercept),
        "n_points": int(x.size),
    }


def read_noise_from_dark_points(points: Sequence[tuple[float, float, float]],
                                ) -> Optional[dict]:
    """Read noise from dark-pair PTC points: the variance floor.

    ``points`` = (exptime, level, var) tuples.  At the very bottom of the
    level axis of the SHORTEST dark exposure the pixels are bias-offset +
    read noise + a residual dark-current shot term that is small but NOT
    zero (measured ~1% of the read-noise variance for the archive's 8 s
    shortest darks — quantify it with :func:`dark_shot_fraction`, do not
    assert "none").  Uses the points whose level lies within 5% of the
    group minimum; returns the read noise in ADU (sqrt of the median
    variance), a MAD-based STATISTICAL uncertainty (callers add the
    dark-shot systematic), the offset level itself (the bias pedestal made
    visible) with its floor-group spread, and the counts behind them.
    None when no dark points exist.
    """
    pts = [(e, l, v) for e, l, v in points
           if e is not None and l is not None and v is not None and v >= 0]
    if not pts:
        return None
    e_min = min(p[0] for p in pts)
    grp = [p for p in pts if p[0] == e_min]
    l_min = min(p[1] for p in grp)
    floor = [p[2] for p in grp if p[1] <= 1.05 * l_min]
    floor_levels = [p[1] for p in grp if p[1] <= 1.05 * l_min]
    var_med = float(np.median(floor))
    var_mad = float(1.4826 * np.median(np.abs(np.asarray(floor) - var_med)))
    rn = float(np.sqrt(max(var_med, 0.0)))
    return {"read_noise_adu": rn,
            "read_noise_adu_err": (var_mad / (2 * rn)) if rn > 0 else None,
            "offset_adu": float(l_min),
            # Half-spread of the floor bins' levels: how well "the offset"
            # is one number at all (per-pixel offset structure folded in).
            "offset_adu_err": float((max(floor_levels) - min(floor_levels))
                                    / 2.0),
            "exptime": float(e_min),
            "n_points": len(floor)}


def dark_shot_fraction(points: Sequence[tuple[float, float, float]],
                       ) -> Optional[dict]:
    """How much dark-current shot noise contaminates the shortest-dark floor.

    The variance floor grows linearly with exposure time (dark shot
    variance ~ dark rate x t); comparing the SHORTEST and LONGEST dark
    floors therefore measures the dark-shot term hiding INSIDE the
    shortest floor:

        shot_var(t_short) = (var_long - var_short) * t_short/(t_long - t_short)

    Returns that variance, its fraction of the shortest floor (the honest
    caveat on any "the floor is pure read noise" claim — measured ~1% for
    the archive's 8 s vs 128 s High Gain darks), and the read-noise BIAS in
    ADU it implies (d(sqrt(v)) = dv / (2 sqrt(v))).  None when fewer than
    two distinct exposure times exist or the long floor is not above the
    short one (no measurable growth).
    """
    pts = [(e, l, v) for e, l, v in points
           if e is not None and l is not None and v is not None and v >= 0]
    if not pts:
        return None
    times = sorted({p[0] for p in pts})
    if len(times) < 2:
        return None
    def floor_var(t):
        grp = [p for p in pts if p[0] == t]
        l_min = min(p[1] for p in grp)
        return float(np.median([p[2] for p in grp if p[1] <= 1.05 * l_min]))
    t_s, t_l = times[0], times[-1]
    v_s, v_l = floor_var(t_s), floor_var(t_l)
    if not (v_l > v_s > 0 and t_l > t_s):
        return None
    shot = (v_l - v_s) * t_s / (t_l - t_s)
    return {"shot_var_adu2": float(shot),
            "frac_of_floor": float(shot / v_s),
            "rn_bias_adu": float(shot / (2.0 * np.sqrt(v_s))),
            "t_short": float(t_s), "t_long": float(t_l),
            "var_short": v_s, "var_long": v_l}


def stackpro_signature(dark_points_sp: Sequence[tuple[float, float, float]],
                       dark_points_base: Sequence[tuple[float, float, float]],
                       ceiling_sp: Optional[float] = None,
                       ceiling_base: Optional[float] = None,
                       ) -> Optional[dict]:
    """N_sub from three independent StackPro/HighGain ratios.

    If a StackPro frame is the SUM of N_sub sub-exposures, then relative
    to a single High Gain read: (1) the bias offset multiplies by N_sub,
    (2) the read-noise VARIANCE multiplies by N_sub, and (3) the
    saturation ceiling multiplies by N_sub.  All three are measured
    quantities here (dark-pair PTC floors and the ceiling memo), so the
    consensus integer is triple-checked; the per-ratio misfits expose any
    disagreement.  Returns None when either mode lacks dark points.
    """
    rn_sp = read_noise_from_dark_points(dark_points_sp)
    rn_b = read_noise_from_dark_points(dark_points_base)
    if rn_sp is None or rn_b is None:
        return None
    ratios = {
        "offset_ratio": rn_sp["offset_adu"] / rn_b["offset_adu"],
        "rn_var_ratio": (rn_sp["read_noise_adu"] ** 2
                         / max(rn_b["read_noise_adu"] ** 2, 1e-12)),
    }
    if ceiling_sp and ceiling_base:
        ratios["ceiling_ratio"] = ceiling_sp / ceiling_base
    consensus = int(round(float(np.median(list(ratios.values())))))
    return {**ratios, "nsub": consensus,
            "max_misfit": float(max(abs(v - consensus)
                                    for v in ratios.values()))}


def nsub_estimate(gain_apparent: float, gain_true: float) -> dict:
    """N_sub under the AVERAGING hypothesis (the model the data rejected).

    An AVERAGE of N_sub sub-frames suppresses variance by N_sub, which the
    PTC reads as an apparent gain N_sub times larger than the true single-
    read gain.  Returns the raw ratio and its nearest integer (the physical
    N_sub must be an integer; the distance to it is the sanity check).

    NOTE: the archive's StackPro frames were determined to be SUMS, not
    averages (see :func:`stackpro_signature` and the module docstring) —
    for a sum the PTC slope is UNCHANGED and this estimator reads ~1, so it
    serves as the discriminator between the two architectures, not as the
    adopted N_sub measurement.
    """
    ratio = gain_apparent / gain_true
    return {"ratio": float(ratio), "nsub": int(round(ratio)),
            "misfit": float(abs(ratio - round(ratio)))}


def amp_glow_metric(dark: np.ndarray, edge: int = 256) -> dict:
    """Amp-glow check on one long dark: edge-region excess over the center.

    Amplifier glow lives against the sensor edges (readout electronics
    corners); a spatially uniform dark shows ~zero excess.  Returns the
    median level of the four ``edge``-wide borders, the central median, and
    their difference in ADU — the number the report quotes.
    """
    d = np.asarray(dark, dtype=np.float64)
    ny, nx = d.shape
    center = d[ny // 4: 3 * ny // 4, nx // 4: 3 * nx // 4]
    borders = np.concatenate([
        d[:edge, :].ravel(), d[-edge:, :].ravel(),
        d[edge:-edge, :edge].ravel(), d[edge:-edge, -edge:].ravel()])
    c = float(np.median(center))
    e = float(np.median(borders))
    # Corner hot spots are the classic glow signature: report the hottest
    # corner's median too (an edge median can hide a single glowing corner).
    k = edge
    corners = [d[:k, :k], d[:k, -k:], d[-k:, :k], d[-k:, -k:]]
    hottest = max(float(np.median(cn)) for cn in corners)
    return {"center_med": c, "edge_med": e, "edge_excess": e - c,
            "hottest_corner_med": hottest, "hottest_corner_excess": hottest - c}
