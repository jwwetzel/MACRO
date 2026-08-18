"""Pure empirical-error-model arithmetic — the S5 seed.

Three validation statistics, computed from CHECK STARS only (the targets
are polars whose orbital modulation is SIGNAL — a target never validates
its own error bars):

1.  Check-star RMS vs magnitude per (era, filter), against the photon-noise
    floor predicted by the measurement errors: where the cloud of constant
    stars sits above the floor, the floor is optimistic.
2.  Reduced chi-square of the constant-star hypothesis -> the error
    INFLATION FACTOR sqrt(median chi2_nu): multiply formal errors by this
    to make a constant star's chi2_nu = 1.
3.  Allan deviation of one long night of one check star: white noise
    integrates down as tau^-1/2; a flattening reveals the correlated-noise
    floor that no amount of averaging removes.

All pure numpy; unit-tested in ``pipeline/tests/test_phot.py`` (white noise
must yield inflation ~1 and an Allan slope of -1/2).
"""

from __future__ import annotations

import math
from typing import Optional

import numpy as np

# --------------------------------------------------------------------------
# Tunables (the report interpolates these).
# --------------------------------------------------------------------------

#: Magnitude-bin width for the RMS-vs-mag floor curves.
RMS_MAG_BIN = 0.5

#: Allan deviation: shortest and longest averaging bins, as counts of
#: consecutive points; the ladder is roughly geometric.
ALLAN_MIN_PTS_PER_BIN = 1
ALLAN_MIN_BINS = 4


def inflation_factor(chi2nu: np.ndarray) -> float:
    """Error inflation factor from check-star reduced chi-squares.

    ``sqrt(median chi2_nu)``: the median resists the one check star that is
    secretly variable; the square root converts a variance ratio into the
    multiplicative factor for the errors themselves.  NaN when no finite
    chi2 is supplied.
    """
    c = np.asarray(chi2nu, dtype=float)
    c = c[np.isfinite(c)]
    if c.size == 0:
        return float("nan")
    return float(math.sqrt(np.median(c)))


def rms_vs_mag_curve(mags: np.ndarray, values: np.ndarray,
                     bin_width: float = RMS_MAG_BIN
                     ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median of ``values`` in magnitude bins (the floor-curve generator).

    Returns ``(bin_center, median_value, n_in_bin)`` for the non-empty
    bins; used both for the check/field-star RMS cloud's running median and
    for the photon-noise floor (feed predicted sigmas as ``values``).
    """
    mags = np.asarray(mags, dtype=float)
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(mags) & np.isfinite(values)
    mags, values = mags[ok], values[ok]
    if mags.size == 0:
        return np.array([]), np.array([]), np.array([])
    lo = math.floor(mags.min() / bin_width) * bin_width
    idx = np.floor((mags - lo) / bin_width).astype(int)
    centers, medians, counts = [], [], []
    for b in np.unique(idx):
        sel = idx == b
        centers.append(lo + (b + 0.5) * bin_width)
        medians.append(float(np.median(values[sel])))
        counts.append(int(sel.sum()))
    return np.array(centers), np.array(medians), np.array(counts)


def allan_deviation(y: np.ndarray, dt_s: float
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Non-overlapping Allan deviation of an (approximately) even series.

    Parameters
    ----------
    y
        Magnitude series of ONE star through ONE night, in observation
        order (finite values only — the caller filters).
    dt_s
        Median cadence of the series in seconds (sets the tau axis).

    Returns
    -------
    (tau_s, adev, n_pairs)
        For each averaging factor m in a geometric ladder: the averaging
        time m*dt_s, the Allan deviation

            sigma_A(m) = sqrt( 0.5 * mean( (ybar_{k+1} - ybar_k)^2 ) )

        over consecutive non-overlapping bins of m points, and the number
        of difference pairs entering the mean.  Ladders stop while at least
        :data:`ALLAN_MIN_BINS` bins remain (fewer differences than that is
        gossip, not statistics).

    For pure white noise of per-point sigma s, sigma_A(m) = s / sqrt(m):
    a -1/2 log-log slope.  Flattening at large tau = correlated noise.
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    n = y.size
    taus, adevs, npairs = [], [], []
    m = ALLAN_MIN_PTS_PER_BIN
    while n // m >= ALLAN_MIN_BINS:
        k = n // m                       # complete bins of m points
        means = y[:k * m].reshape(k, m).mean(axis=1)
        d = np.diff(means)
        taus.append(m * dt_s)
        adevs.append(float(np.sqrt(0.5 * np.mean(d ** 2))))
        npairs.append(len(d))
        # Geometric ladder: 1, 2, 4, 8, ... (doubling keeps bins independent
        # and the points evenly spaced in log tau).
        m *= 2
    return np.array(taus), np.array(adevs), np.array(npairs)


def longest_run(jds: np.ndarray, max_gap_s: float = 1800.0
                ) -> tuple[int, int]:
    """Longest contiguous observing run in a sorted JD array.

    Returns the (start, stop) slice indices of the longest stretch whose
    consecutive gaps never exceed ``max_gap_s`` — the stretch handed to the
    Allan analysis (a mid-night pause would masquerade as long-tau noise).
    """
    jds = np.asarray(jds, dtype=float)
    if jds.size == 0:
        return 0, 0
    gaps = np.diff(jds) * 86400.0
    breaks = np.flatnonzero(gaps > max_gap_s)
    starts = np.concatenate(([0], breaks + 1))
    stops = np.concatenate((breaks + 1, [jds.size]))
    best = int(np.argmax(stops - starts))
    return int(starts[best]), int(stops[best])


def predicted_sigma_with_zp(sig_star: np.ndarray,
                            zp_err: Optional[np.ndarray] = None,
                            floor_mag: float = 0.0) -> np.ndarray:
    """Total predicted per-point sigma: photon + zero-point + floor, in
    quadrature.  This is the denominator the chi2 verdicts use, so the
    inflation factor speaks about everything the formal model claims."""
    v = np.square(np.asarray(sig_star, dtype=float))
    if zp_err is not None:
        v = v + np.square(np.asarray(zp_err, dtype=float))
    return np.sqrt(v + floor_mag ** 2)
