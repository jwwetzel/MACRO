"""be_tscore — the ONE periodogram engine of the Be-star paper, shared by the
injection–recovery run (BE-S-1b) and the real slow-tier search (BE-S11), so
the injections pass through exactly the code the data will.

THE MODEL (why it is a joint fit and not detrend-then-search)
-------------------------------------------------------------
Strategy §4 Step 11 forbids detrend-then-search: a regression fitted first can
absorb a real signal, or manufacture one.  So, at each trial frequency f, the
nightly series y (n nights) is fitted with

    y = N @ beta + a sin(2 pi f t) + b cos(2 pi f t) + noise,

where N is the NUISANCE design: one free offset per mechanical state (the
ruled instrument boundaries — a constant EW offset across a boundary is never
assumed) and the systematics regressors (recomputed airmass, focus position,
CCD temperature, telluric H2O band depth when measured), each standardised.
The power is the fraction of the nuisance-only residual variance that the
sinusoid removes:

    P(f) = 1 - chi2(nuisance + sinusoid) / chi2(nuisance),

computed exactly by projecting the sinusoid columns onto the orthogonal
complement of N (Frisch–Waugh–Lovell), vectorised over all frequencies.  With
N = a single column of ones this is the floating-mean generalised Lomb–Scargle.

SIGNIFICANCE IS GLOBAL
----------------------
``max_power_null`` returns the distribution of the MAXIMUM power over the whole
frequency band (the max-statistic, so the trials factor is inside the FAP —
DS memo, BE-S11), either from Gaussian noise realisations through the same
projection (the injection run, where the noise IS Gaussian) or from a
night-label bootstrap of the real residuals (the real search).

Everything here is pure numpy on arrays; no database, no files.
"""
from __future__ import annotations

import numpy as np


def nuisance_design(states, regressors=None, max_frac: float = 1 / 3):
    """Columns: one indicator per mechanical state present, then each regressor
    standardised.  A regressor with no variance left after the state offsets is
    dropped (it would be degenerate), and the total column count is capped at
    ``max_frac`` of the nights so the nuisance model cannot eat the data.

    Returns (N, names)."""
    states = np.asarray(states)
    cols, names = [], []
    for s in sorted(set(states)):
        cols.append((states == s).astype(float))
        names.append(f"offset[{s}]")
    N = np.column_stack(cols)
    if regressors:
        for name, x in regressors.items():
            x = np.asarray(x, float)
            if np.any(~np.isfinite(x)):
                continue
            # residual of the regressor after the state offsets
            r = x - N @ np.linalg.lstsq(N, x, rcond=None)[0]
            sd = r.std()
            if sd < 1e-9 * (np.abs(x).max() + 1) or sd == 0:
                continue
            if N.shape[1] + 1 > max(len(states) * max_frac, len(set(states))):
                break
            N = np.column_stack([N, (x - x.mean()) / x.std()])
            names.append(name)
    return N, names


def projector(N):
    """Q = I - N (N^T N)^+ N^T — removes the nuisance subspace."""
    pinv = np.linalg.pinv(N.T @ N)
    return np.eye(N.shape[0]) - N @ pinv @ N.T


def power(t, y, Q, freqs):
    """Joint-fit power P(f) for one series (y may be 2-D: n x m realisations).

    Returns (P, a, b): power and the sine/cosine amplitudes at each frequency
    (shapes F or F x m)."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    one = y.ndim == 1
    Y = y[:, None] if one else y                        # n x m
    Yr = Q @ Y                                          # residual after nuisance
    ph = 2 * np.pi * np.outer(t - t.mean(), freqs)      # n x F
    S, Cc = Q @ np.sin(ph), Q @ np.cos(ph)              # residualised columns
    ss, cc, sc = (S * S).sum(0), (Cc * Cc).sum(0), (S * Cc).sum(0)   # F
    sy, cy = S.T @ Yr, Cc.T @ Yr                        # F x m
    det = ss * cc - sc ** 2
    det = np.where(det <= 1e-12 * (ss * cc + 1e-300), np.nan, det)
    a = (cc[:, None] * sy - sc[:, None] * cy) / det[:, None]
    b = (ss[:, None] * cy - sc[:, None] * sy) / det[:, None]
    explained = a * sy + b * cy                         # chi2 reduction
    tot = (Yr * Yr).sum(0)[None, :]
    P = explained / tot
    if one:
        return P[:, 0], a[:, 0], b[:, 0]
    return P, a, b


def max_power_null_gaussian(t, Q, freqs, n_draw: int, rng, chunk: int = 500):
    """Distribution of max_f P(f) under white Gaussian noise, through Q."""
    out = []
    n = len(t)
    for k in range(0, n_draw, chunk):
        m = min(chunk, n_draw - k)
        P, _, _ = power(t, rng.standard_normal((n, m)), Q, freqs)
        out.append(np.nanmax(P, axis=0))
    return np.concatenate(out)


def max_power_null_bootstrap(t, resid, Q, freqs, n_draw: int, rng, chunk: int = 500):
    """Night-label bootstrap: the nuisance-model residuals are permuted among
    nights (time order destroyed, the value distribution kept), and the max
    power over the band is recorded for each draw."""
    out = []
    resid = np.asarray(resid, float)
    for k in range(0, n_draw, chunk):
        m = min(chunk, n_draw - k)
        Y = np.column_stack([rng.permutation(resid) for _ in range(m)])
        P, _, _ = power(t, Y, Q, freqs)
        out.append(np.nanmax(P, axis=0))
    return np.concatenate(out)


def freq_grid(t, f_max: float, oversample: float = 10.0, p_max: float | None = None):
    """Uniform frequency grid from 1/p_max (default: 1/baseline) to f_max with
    spacing 1/(oversample * baseline)."""
    T = float(np.ptp(t))
    f0 = 1.0 / (p_max if p_max else T)
    df = 1.0 / (oversample * T)
    return np.arange(f0, f_max + df, df)
