"""Pure S2 empirical-noise logic: measured counts-vs-variance, per mode.

WHY THIS MODULE EXISTS (and why it is not :mod:`rlmt_diagnostics.ptc`)
----------------------------------------------------------------------
:mod:`~rlmt_diagnostics.ptc` answers a *physics* question — what is the
gain, what is the read noise — and to answer it cleanly it needs one very
special night (2023-06-07, repeated darks and repeated Albireo fields, all
of them ``High Gain``).  That night characterises two readout modes and
leaves the other six with no noise statement at all.

The CV time-series project needs a different thing, and needs it for every
mode it actually observed in.  Its per-point error bars are only as good as
its answer to: *at a measured level of L ADU, how much does a pixel of this
detector, in this readout mode, actually fluctuate?*  That question does
not require knowing the gain.  It can be answered by MEASUREMENT ALONE, and
this module is deliberately built to answer it that way:

    take two frames of the same scene seconds apart, bin their pixels on
    level, and record var(a - b)/2 per bin.

The result is a TABLE — level in, variance out — not a formula.  Nothing
here assumes Poisson statistics, a gain value, or a read-noise term; if the
detector does something ugly (a bias wobble, a fixed-pattern residual, a
mode whose variance is flat because the frames are stacked sums), the table
says so and the error model inherits the truth rather than the theory.
:func:`predict_variance` is the interpolator downstream code calls, and it
is honest at the edges: outside the measured level span it says so instead
of extrapolating a curve nobody measured.

The formula-based reading is still available as a CROSS-CHECK — fitting a
straight line to the same points recovers the classic
``var = read_noise^2 + level/gain`` — but that fit is reported next to the
table, never in place of it.

Everything here is a pure function on plain sequences / numpy arrays: no
I/O, no database, no globals.  The campaign script owns frame pairing and
pixel reading; the unit tests own synthetic frames with known truth.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Level bins in a mode's empirical curve.  Log-spaced across the measured
#: span: enough resolution to show curvature, few enough that every bin
#: still pools many pixels from many pairs.
NOISE_CURVE_BINS = 14

#: A curve bin needs points from at least this many distinct frame PAIRS
#: before it is published.  One pair's variance can be wrecked by a single
#: passing cloud, a satellite trail, or a tracking bump; requiring several
#: independent pairs makes the bin a measurement instead of an anecdote.
MIN_PAIRS_PER_BIN = 2

#: ... and at least this many (level, variance) points in total.
MIN_POINTS_PER_BIN = 3

#: A mode needs at least this many usable points before a curve is built
#: at all.  Below it the honest output is "not measured for this mode".
MIN_POINTS_PER_MODE = 12

#: The floor is read from the lowest levels of the measured span: every
#: point within this factor of the minimum measured level.  1.30 keeps the
#: sky-floor bins together without reaching up into the star wings.
FLOOR_LEVEL_FACTOR = 1.30

#: MAD -> sigma conversion for a Gaussian (1 / Phi^-1(3/4)).
MAD_TO_SIGMA = 1.4826

#: Crossover definition: the level at which the measured variance reaches
#: this multiple of the measured floor variance.  x2 is the classic
#: "shot noise now equals read noise" point, but note it is READ OFF THE
#: MEASURED CURVE, not computed from a gain.
CROSSOVER_VAR_MULTIPLE = 2.0

#: A curve must span at least this RATIO of levels (top / bottom) before a
#: log-log slope means anything.  Measured on the archive: the 5 MHz iKon's
#: science pairs all sit on a flat ~6,000 ADU sky, giving a span ratio of
#: 1.05 — over which a straight-line fit happily returned a "slope" of 48,
#: a number with no physical reading whatsoever.  Below this ratio the
#: honest answer is that the mode was measured at ONE level, not across a
#: range, and the slope is refused.
MIN_SPAN_RATIO_FOR_SLOPE = 2.0


def empirical_noise_curve(points: Sequence[tuple],
                          n_bins: int = NOISE_CURVE_BINS,
                          ) -> list[dict]:
    """The measured variance-vs-level table for one readout mode.

    Parameters
    ----------
    points
        ``(level, var, n_pix, pair_id)`` tuples — one per level bin of one
        frame pair, as produced by :func:`rlmt_diagnostics.ptc.pair_ptc_points`
        and tagged with the pair they came from.
    n_bins
        Number of log-spaced level bins spanning the measured range.

    Returns
    -------
    list of dict
        One row per POPULATED bin, in increasing level order, each with:

        * ``level``      — median level of the bin's points (ADU);
        * ``var``        — median half-difference variance (ADU^2);
        * ``var_mad``    — MAD-sigma of those variances: the honest spread
          between independent pairs at the same level, i.e. how repeatable
          the measurement is;
        * ``sigma``      — sqrt(var), the number an error model wants;
        * ``n_points``   — points pooled into the bin;
        * ``n_pairs``    — DISTINCT frame pairs behind them;
        * ``n_pix``      — pixels behind them.

    Bins failing :data:`MIN_PAIRS_PER_BIN` / :data:`MIN_POINTS_PER_BIN` are
    dropped rather than published thin — a variance from one pair is not a
    measurement of the detector, it is a measurement of that night.
    """
    # Keep only physically meaningful points: a positive level (log bins)
    # and a non-negative variance.  NaNs from a degenerate pair vanish here.
    rows = [(float(l), float(v), float(n), str(p))
            for l, v, n, p in points
            if l is not None and v is not None
            and np.isfinite(l) and np.isfinite(v) and l > 0 and v >= 0]
    if len(rows) < MIN_POINTS_PER_MODE:
        return []
    levels = np.array([r[0] for r in rows])
    lo, hi = float(levels.min()), float(levels.max())
    if not hi > lo:
        return []
    # Log-spaced edges: detector noise spans decades of level, so equal
    # RATIOS (not equal differences) are the natural bin width.
    edges = np.geomspace(lo, hi * 1.000001, n_bins + 1)
    idx = np.digitize(levels, edges) - 1
    idx = np.clip(idx, 0, n_bins - 1)
    out: list[dict] = []
    for b in range(n_bins):
        sel = [r for r, k in zip(rows, idx) if k == b]
        if len(sel) < MIN_POINTS_PER_BIN:
            continue
        pairs = {r[3] for r in sel}
        if len(pairs) < MIN_PAIRS_PER_BIN:
            continue
        v = np.array([r[1] for r in sel])
        v_med = float(np.median(v))
        out.append({
            "level": float(np.median([r[0] for r in sel])),
            "var": v_med,
            "var_mad": float(MAD_TO_SIGMA * np.median(np.abs(v - v_med))),
            "sigma": float(np.sqrt(max(v_med, 0.0))),
            "n_points": len(sel),
            "n_pairs": len(pairs),
            "n_pix": float(sum(r[2] for r in sel)),
        })
    return out


def noise_floor(curve: Sequence[dict],
                factor: float = FLOOR_LEVEL_FACTOR) -> Optional[dict]:
    """The measured noise floor at the bottom of a mode's curve.

    Every curve bin whose level lies within ``factor`` of the lowest
    measured level is pooled; the floor is the median of their variances.
    This is the mode's read-plus-dark-plus-whatever-else floor AS MEASURED
    — it is deliberately not called "read noise", because a science frame's
    faintest sky is not a zero-second dark and this module refuses to
    pretend otherwise.  Returns the floor variance, its square root in ADU,
    the level it was measured at, and the bins behind it; None for an empty
    curve.
    """
    rows = [c for c in curve if c["level"] > 0]
    if not rows:
        return None
    l_min = min(c["level"] for c in rows)
    grp = [c for c in rows if c["level"] <= factor * l_min]
    var = float(np.median([c["var"] for c in grp]))
    # Spread ACROSS the pooled bins is the honest uncertainty of a floor
    # read this way (a single bin's var_mad would understate it).
    spread = ([c["var"] for c in grp])
    var_err = (float((max(spread) - min(spread)) / 2.0)
               if len(spread) > 1 else float(grp[0]["var_mad"]))
    sigma = float(np.sqrt(max(var, 0.0)))
    return {
        "floor_var_adu2": var,
        "floor_var_err_adu2": var_err,
        "floor_sigma_adu": sigma,
        # d(sqrt(v)) = dv / (2 sqrt(v)) — propagate the spread into ADU.
        "floor_sigma_err_adu": (var_err / (2.0 * sigma)) if sigma > 0 else None,
        "level_adu": float(np.median([c["level"] for c in grp])),
        "n_bins": len(grp),
        # True when the floor window swallowed the ENTIRE curve: the mode
        # was measured at one level, so this is "the noise at that level",
        # not a floor the curve later rises off.  Anything derived from the
        # floor-vs-rest CONTRAST (the crossover) is circular in that case
        # and callers must decline to publish it.
        "is_whole_curve": len(grp) == len(rows),
    }


def crossover_level(curve: Sequence[dict], floor_var: float,
                    multiple: float = CROSSOVER_VAR_MULTIPLE,
                    ) -> Optional[float]:
    """Level at which the measured variance reaches ``multiple`` x the floor.

    Read off the MEASURED curve by linear interpolation between the two
    bracketing bins — no gain, no Poisson assumption.  Below this level a
    pixel's error is floor-dominated (exposure time buys nothing); above
    it, signal shot noise dominates.  Returns None when the curve never
    reaches the threshold within its measured span (an honest "not
    witnessed" instead of an extrapolation).
    """
    if floor_var <= 0:
        return None
    target = multiple * floor_var
    rows = sorted(curve, key=lambda c: c["level"])
    for a, b in zip(rows, rows[1:]):
        if a["var"] < target <= b["var"]:
            span = b["var"] - a["var"]
            if span <= 0:
                return float(b["level"])
            f = (target - a["var"]) / span
            return float(a["level"] + f * (b["level"] - a["level"]))
    return None


def predict_variance(curve: Sequence[dict], level: float) -> Optional[float]:
    """The empirical error model itself: measured variance at a level.

    Log-log linear interpolation between the measured bins.  Returns None
    OUTSIDE the measured span — the whole point of a measured model is that
    it declines to invent numbers where nobody looked.  Callers that need a
    value beyond the span must either measure more frames or state that
    they extrapolated.
    """
    rows = sorted((c for c in curve if c["level"] > 0 and c["var"] > 0),
                  key=lambda c: c["level"])
    if len(rows) < 2 or level is None or not np.isfinite(level) or level <= 0:
        return None
    if level < rows[0]["level"] or level > rows[-1]["level"]:
        return None
    x = np.log([c["level"] for c in rows])
    y = np.log([c["var"] for c in rows])
    return float(np.exp(np.interp(np.log(level), x, y)))


def level_span_ratio(curve: Sequence[dict]) -> Optional[float]:
    """Top measured level divided by the bottom one.

    How much of a RANGE the mode was actually measured over — the number
    that decides whether the curve is a curve at all or a single point with
    error bars.  None for an empty curve.
    """
    rows = [c["level"] for c in curve if c["level"] > 0]
    if not rows:
        return None
    return float(max(rows) / min(rows))


def curve_shape_index(curve: Sequence[dict],
                      min_span: float = MIN_SPAN_RATIO_FOR_SLOPE,
                      ) -> Optional[float]:
    """Log-log slope of variance against level over the measured span.

    A pure-Poisson detector reads 1.0 (variance proportional to signal); a
    floor-dominated span reads ~0; a span whose variance grows faster than
    the signal (scene motion between the pair, fixed-pattern residual,
    tracking drift) reads >1 and is a WARNING that the pair difference is
    measuring the sky and the mount as well as the detector.  Reported
    beside every curve so a reader can see which regime they are in
    without re-deriving it.

    Returns None when there are fewer than two usable bins OR when the
    measured levels span less than ``min_span`` (see
    :data:`MIN_SPAN_RATIO_FOR_SLOPE`): a slope fitted across a few percent
    of level is a property of the noise in the fit, not of the detector,
    and reporting it would dress a meaningless number as a measurement.
    """
    rows = [c for c in curve if c["level"] > 0 and c["var"] > 0]
    if len(rows) < 2:
        return None
    span = level_span_ratio(rows)
    if span is None or span < min_span:
        return None
    x = np.log([c["level"] for c in rows])
    y = np.log([c["var"] for c in rows])
    xm, ym = x.mean(), y.mean()
    sxx = float(((x - xm) ** 2).sum())
    if sxx <= 0:
        return None
    return float(((x - xm) * (y - ym)).sum() / sxx)
