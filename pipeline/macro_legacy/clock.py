"""Per-season clock audit of the legacy archive — the pure logic (RIG-L1-clock-audit).

WHAT THE AUDIT ASKS
-------------------
The census (``macro_legacy.census``) established from headers what every
legacy ``DATE-OBS`` MEANS: UTC, the start of the exposure — except the
AC4040 under MaxIm 6.28–6.30, where two independent header tests say the
middle.  Headers cannot say whether the acquisition PC's clock was RIGHT.
This module measures that, per observing season, from eclipses of
post-common-envelope binaries (HW Vir-type sdB + dM, and WD + dM) that
students happened to observe.  Their primary eclipses are deep, sharp and
short, so one well-sampled minimum is timed to tens of seconds or better.

THE REFERENCE IS TESS, NOT A LITERATURE EPHEMERIS
-------------------------------------------------
The ledger's criterion is "O − C within 60 s of the literature per season,
or the offset is carried".  For these systems a literature linear
ephemeris is not good to 60 s: every one of them shows eclipse-timing
variations of tens of seconds to minutes (the reason S3b refused NN Ser and
QS Vir as clocks).  The audit therefore predicts each legacy minimum from
eclipse times MEASURED in TESS 2-min photometry of the same star — absolute
BJD_TDB, independent of any observatory clock — fitted with a local
ephemeris that spans or nearly spans the legacy season.  How far a
prediction may be carried outside the TESS span, and what systematic it
then carries, is fixed below (:data:`EXTRAP_MAX_D`,
:func:`reference_prediction`) before any legacy minimum was timed.

Sign convention (the same as S3b): ``O − C = T_observed − T_predicted`` in
seconds, ``T_observed`` read off OUR time axis.  A PC clock running AHEAD
of true time by D seconds gives ``O − C = +D``.

All functions are pure; unit tests in ``pipeline/tests/test_legacy_clock.py``.
"""

from __future__ import annotations

import math
import re
from typing import Optional, Sequence

import numpy as np

SECONDS_PER_DAY = 86400.0

#: The ledger's acceptance criterion for a season (RIG-L1-clock-audit).
CLOCK_ACCEPT_S = 60.0

#: VSX types admitted as clock targets: HW Vir-type sdB + dM and WD + dM
#: binaries (deep, sharp, short primary eclipses).  W UMa-type contact
#: binaries are NOT admitted: their broad minima and secular period changes
#: make a 60 s clock test impossible from an ephemeris years away.
CLOCK_TYPE_RE = re.compile(r"EA/HW|EA/WD|EA\+UV|\bHW\b|/WD", re.IGNORECASE)

#: A legacy minimum is predicted from the TESS ephemeris only if it lies
#: within this many days OUTSIDE the span of the TESS eclipse times.  Inside
#: the span the prediction is an interpolation; beyond this margin the
#: star's own eclipse-timing variations (tens of s per year in this class)
#: make the prediction no clock reference at 60 s, and the event is
#: recorded as out of range rather than timed against a guess.
EXTRAP_MAX_D = 400.0

#: When the TESS times come from a single sector, the ephemeris is a
#: 27-day local one; it may be carried this far only.
EXTRAP_MAX_SINGLE_SECTOR_D = 60.0

#: Phase-folding window about the primary eclipse, in units of the period
#: (excludes the secondary eclipse at phase 0.5 and leaves baseline for the
#: polynomial on both sides of a ~0.1 P primary).
FOLD_HALF_WINDOW_P = 0.3

#: Width of the bins the folded TESS light curve is averaged into (s).
#: Far below the eclipse ingress (several minutes), so binning cannot move
#: the mid-time of a symmetric event; it makes the fit and the bootstrap
#: cheap on ~10^4 points.
FOLD_BIN_S = 30.0

#: Legacy-run coverage required about a predicted minimum, in units of the
#: period on EACH side.  The primary of an HW Vir system lasts ~0.1 P, so
#: 0.15 P each side leaves ingress, egress and some baseline in the run.
RUN_COVER_P = 0.15


def is_clock_type(vsx_type: Optional[str]) -> bool:
    """Is a VSX variability type an admitted clock target?"""
    return bool(vsx_type) and bool(CLOCK_TYPE_RE.search(str(vsx_type)))


def season_label(night: str) -> str:
    """Observing season of a night label: July–June, e.g. '2018/19'.

    Winer's seasons straddle the new year and the monsoon closes the dome
    in July–August, so a July–June season never splits one observing run.
    """
    y, m = int(night[:4]), int(night[5:7])
    start = y if m >= 7 else y - 1
    return f"{start}/{(start + 1) % 100:02d}"


# ---------------------------------------------------------------------------
# TESS eclipse times
# ---------------------------------------------------------------------------
def fold_primary(t: np.ndarray, t_ref: float, period: float,
                 half_window_p: float = FOLD_HALF_WINDOW_P
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Shift every point onto the cycle of ``t_ref`` (the reference minimum).

    Returns ``(t_folded, keep)``: ``t_folded = t − (E − E_ref)·P`` with E the
    nearest cycle of each point, and ``keep`` the points within
    ``half_window_p`` periods of the reference minimum.  The folded times
    are absolute BJD on the reference cycle, so a fit to them returns the
    mid-time of THAT eclipse directly.
    """
    t = np.asarray(t, dtype=float)
    cyc = np.round((t - t_ref) / period)
    tf = t - cyc * period
    return tf, np.abs(tf - t_ref) <= half_window_p * period


def bin_series(t: np.ndarray, f: np.ndarray, width_d: float
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Average points into fixed-width bins.

    Returns ``(t_mean, f_mean, f_err, n)`` per bin with at least 2 points;
    ``f_err`` is the empirical standard error (std / sqrt n), never a
    propagated pipeline error, so red noise in the bin shows up in it.
    """
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    idx = np.floor((t - t.min()) / width_d).astype(int)
    order = np.argsort(idx)
    idx, t, f = idx[order], t[order], f[order]
    edges = np.flatnonzero(np.diff(idx)) + 1
    tm, fm, fe, nn = [], [], [], []
    for a, b in zip(np.r_[0, edges], np.r_[edges, len(idx)]):
        n = b - a
        if n < 2:
            continue
        tm.append(t[a:b].mean())
        fm.append(f[a:b].mean())
        fe.append(f[a:b].std(ddof=1) / math.sqrt(n))
        nn.append(n)
    return np.array(tm), np.array(fm), np.array(fe), np.array(nn)


# ---------------------------------------------------------------------------
# Ephemeris from TESS eclipse times
# ---------------------------------------------------------------------------
def cycle_numbers(times: Sequence[float], t_anchor: float,
                  period: float) -> np.ndarray:
    """Integer cycle of each time relative to an anchor epoch."""
    return np.round((np.asarray(times, dtype=float) - t_anchor) / period
                    ).astype(int)


def fit_ephemeris(cycles: Sequence[int], times: Sequence[float],
                  sigmas: Sequence[float], order: int = 1) -> dict:
    """Weighted least-squares ephemeris ``T(E) = sum_k c_k E^k`` (order 1 or 2).

    Cycles are centred on their weighted mean before fitting so the
    covariance is well conditioned; the returned coefficients refer to the
    centred cycle ``E − e0``.  Reports chi-square and dof; no error is
    rescaled (standing rule 1).
    """
    e = np.asarray(cycles, dtype=float)
    t = np.asarray(times, dtype=float)
    s = np.asarray(sigmas, dtype=float)
    w = 1.0 / s ** 2
    e0 = float(np.sum(w * e) / np.sum(w))
    x = e - e0
    a = np.vander(x, order + 1, increasing=True)
    aw = a * np.sqrt(w)[:, None]
    tw = t * np.sqrt(w)
    cov = np.linalg.inv(aw.T @ aw)
    coef = cov @ aw.T @ tw
    resid = t - a @ coef
    chi2 = float(np.sum(w * resid ** 2))
    return {"order": order, "e0": e0, "coef": coef, "cov": cov,
            "chi2": chi2, "dof": len(t) - (order + 1),
            "resid_s": resid * SECONDS_PER_DAY,
            "period": float(coef[1]),
            "sig_period": float(math.sqrt(cov[1, 1]))}


def predict(fit: dict, cycle: float) -> tuple[float, float]:
    """``(T, sigma_T)`` in days at a cycle, from :func:`fit_ephemeris`."""
    x = float(cycle) - fit["e0"]
    v = np.array([x ** k for k in range(fit["order"] + 1)])
    return float(v @ fit["coef"]), float(math.sqrt(v @ fit["cov"] @ v))


def reference_prediction(cycle: int, t_obs_approx: float,
                         tess_times: Sequence[float],
                         tess_sigmas: Sequence[float],
                         tess_cycles: Sequence[int],
                         n_sectors: int) -> dict:
    """Predict one legacy minimum from the TESS eclipse times.

    * ``range`` — 'interpolated' (inside the TESS span), 'extrapolated'
      (outside but within :data:`EXTRAP_MAX_D`, or
      :data:`EXTRAP_MAX_SINGLE_SECTOR_D` for a single sector) or
      'out_of_range' (no prediction is made).
    * ``t_pred`` / ``sig_stat_s`` — the linear ephemeris and its formal
      error at that cycle.
    * ``sig_sys_s`` — the model systematic: with TESS times in three or
      more sectors, the difference between the linear and the quadratic
      prediction (curvature the data themselves show); with fewer, the
      rms of the linear residuals (the eclipse-timing scatter seen within
      the TESS span), which is all the data can say.

    Nothing here is tuned to the legacy minimum: only its cycle and its
    approximate time (to decide the range) enter.
    """
    t = np.asarray(tess_times, dtype=float)
    lo, hi = float(t.min()), float(t.max())
    gap = max(lo - t_obs_approx, t_obs_approx - hi, 0.0)
    limit = EXTRAP_MAX_SINGLE_SECTOR_D if n_sectors < 2 else EXTRAP_MAX_D
    if gap == 0.0:
        rng = "interpolated"
    elif gap <= limit:
        rng = "extrapolated"
    else:
        return {"range": "out_of_range", "gap_d": gap}
    lin = fit_ephemeris(tess_cycles, t, tess_sigmas, order=1)
    tp, sp = predict(lin, cycle)
    if n_sectors >= 3 and len(t) >= 4:
        quad = fit_ephemeris(tess_cycles, t, tess_sigmas, order=2)
        tq, _ = predict(quad, cycle)
        sys_s = abs(tq - tp) * SECONDS_PER_DAY
        sys_basis = "linear vs quadratic"
    else:
        r = lin["resid_s"]
        sys_s = float(np.sqrt(np.mean(r ** 2))) if len(r) > 2 else 0.0
        sys_basis = "rms of TESS residuals"
    return {"range": rng, "gap_d": gap, "t_pred": tp,
            "sig_stat_s": sp * SECONDS_PER_DAY, "sig_sys_s": sys_s,
            "sys_basis": sys_basis, "period": lin["period"],
            "chi2": lin["chi2"], "dof": lin["dof"]}


def oc_both_conventions(t0_start_conv: float, t_pred: float,
                        exptime_s: float) -> tuple[float, float]:
    """O − C (s) under the start-of-exposure and the mid-exposure readings.

    The light curve is built on mid-exposure times computed as
    ``DATE-OBS + exptime/2`` (the start reading).  If DATE-OBS actually
    marks the middle, every point — and the fitted mid-time — sits
    ``exptime/2`` late, so the mid reading's O − C is smaller by that much.
    """
    oc_start = (t0_start_conv - t_pred) * SECONDS_PER_DAY
    return oc_start, oc_start - 0.5 * exptime_s


def season_verdict(oc_s: float, sigma_s: float) -> str:
    """Verdict against the ledger's 60 s criterion (S3b's wording).

    * ``PASS`` — the whole 2-sigma interval is inside +-60 s.
    * ``OFFSET MEASURED`` — outside 60 s and excluding zero at 3 sigma; the
      offset is carried.
    * ``PASS (central value); 2-sigma interval reaches the criterion``.
    * ``UNRESOLVED`` — anything else.
    """
    if abs(oc_s) + 2.0 * sigma_s < CLOCK_ACCEPT_S:
        return "PASS"
    if abs(oc_s) >= CLOCK_ACCEPT_S and abs(oc_s) > 3.0 * sigma_s:
        return "OFFSET MEASURED"
    if abs(oc_s) < CLOCK_ACCEPT_S:
        return "PASS (central value); 2-sigma interval reaches the criterion"
    return "UNRESOLVED"
