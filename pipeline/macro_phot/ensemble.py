"""Pure Honeycutt-1992 inhomogeneous-ensemble solver.

Honeycutt (1992, PASP 104, 435) models every measured magnitude of star
``i`` on frame ``j`` as

    m_ij = M_i + ZP_j + noise,

where ``M_i`` is the star's ensemble mean magnitude and ``ZP_j`` the frame's
zero-point offset (clouds, airmass, focus — everything that moves ALL stars
on a frame together).  With inverse-variance weights the least-squares
normal equations decouple into two weighted means, so the solution is found
by alternating

    ZP_j = weighted mean over i of (m_ij - M_i)
    M_i  = weighted mean over j of (m_ij - ZP_j)

to convergence.  The model is invariant under M_i -> M_i + c,
ZP_j -> ZP_j - c (the gauge freedom); we fix the gauge by demanding
``mean(ZP) = 0``, so ensemble magnitudes stay on the instrumental scale.

Robustness: after the plain solution converges, residuals beyond
:data:`CLIP_SIGMA` times the global residual RMS are masked and the solve
repeats — cosmic rays and blended measurements lose their vote without any
star or frame being deleted outright.

Everything here is a pure numpy function on masked arrays; the REQUIRED
synthetic-recovery test (inject a known ZP pattern, demand its recovery
within tolerance) lives in ``pipeline/tests/test_phot.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Alternating-solve convergence: stop when no ZP or M moves by more than
#: this many magnitudes between iterations.
CONVERGE_MMAG = 1e-5

#: Iteration cap for one alternating solve (converges in ~10 normally; the
#: cap turns a pathological matrix into a reported non-convergence, not a
#: hang).
MAX_ITER = 200

#: Residuals beyond this many global RMS are masked in the robust re-solve.
CLIP_SIGMA = 4.0

#: Number of robust (clip + re-solve) passes after the first solve.
CLIP_PASSES = 2

#: Error floor added in quadrature to every measurement error before
#: weighting, in magnitudes.  Prevents one bright star with a formally tiny
#: photon error from owning the zero point (flat-field and PSF systematics
#: are never below a few mmag).
WEIGHT_FLOOR_MAG = 0.005

#: Comparison-star candidacy: a star must appear on at least this fraction
#: of the series' frames (ephemeral stars can't anchor a zero point).
COMP_MIN_FRAME_FRAC = 0.5

#: Stability iteration: drop comps whose residual RMS exceeds
#: COMP_RMS_FACTOR x the median comp RMS (and at least COMP_RMS_FLOOR_MAG,
#: so a set of uniformly excellent comps is not decimated by the factor).
COMP_RMS_FACTOR = 2.5
COMP_RMS_FLOOR_MAG = 0.02

#: Stability iteration cap (each pass drops the newly-unstable; in practice
#: this settles in 2-4 passes).
COMP_MAX_PASSES = 10

#: Minimum surviving comparison stars for a trustworthy ensemble.
MIN_COMPS = 5

#: Check stars held OUT of the zero-point solve to validate it: how many,
#: chosen nearest the target's magnitude among stable stars.
N_CHECK_STARS = 4


# --------------------------------------------------------------------------
# Results container
# --------------------------------------------------------------------------

@dataclass
class EnsembleSolution:
    """One converged ensemble solve.

    Attributes
    ----------
    mean_mag : (S,) per-star ensemble mean magnitude (NaN if never seen).
    zp : (F,) per-frame zero point, gauge-fixed to mean 0.
    zp_err : (F,) formal error of each zero point (from the weights).
    n_star_used : (F,) stars contributing to each frame's ZP.
    residual : (S, F) residual m - M - ZP (NaN where unobserved/clipped).
    clipped : (S, F) bool, True where the robust passes masked a point.
    n_iter : total alternating iterations across all passes.
    converged : True if every pass met CONVERGE_MMAG within MAX_ITER.
    """
    mean_mag: np.ndarray
    zp: np.ndarray
    zp_err: np.ndarray
    n_star_used: np.ndarray
    residual: np.ndarray
    clipped: np.ndarray
    n_iter: int
    converged: bool


def _weights(sig: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Inverse-variance weights with the systematic floor; 0 where masked."""
    var = np.square(np.where(np.isfinite(sig), sig, np.inf)) \
        + WEIGHT_FLOOR_MAG ** 2
    w = np.where(mask, 1.0 / var, 0.0)
    return w


def _alternate(mag: np.ndarray, w: np.ndarray
               ) -> tuple[np.ndarray, np.ndarray, int, bool]:
    """One alternating-least-squares solve given fixed weights.

    Returns (M, ZP, n_iter, converged) with the mean-ZP gauge applied.
    """
    S, F = mag.shape
    m0 = np.where(w > 0, mag, 0.0)          # masked entries contribute 0
    # Initialize M with each star's straight weighted mean (ZP = 0 start).
    with np.errstate(invalid="ignore"):
        M = np.where(w.sum(axis=1) > 0,
                     (m0 * w).sum(axis=1) / w.sum(axis=1), np.nan)
    ZP = np.zeros(F)
    converged = False
    it = 0
    for it in range(1, MAX_ITER + 1):
        Mv = np.where(np.isfinite(M), M, 0.0)
        # ZP_j from stars with known M only (a star unseen elsewhere has
        # NaN M and must not vote).  A frame with NO voting stars at all
        # gets NaN, never a fabricated number: 0.0 here once masqueraded
        # as a perfectly ordinary zero point and let measurements on
        # comp-less frames flow into light curves uncorrected.
        wj = np.where(np.isfinite(M)[:, None], w, 0.0)
        denom_f = wj.sum(axis=0)
        with np.errstate(invalid="ignore"):
            ZP_new = np.where(denom_f > 0,
                              ((m0 - Mv[:, None]) * wj).sum(axis=0) / denom_f,
                              np.nan)
        # M_i against the fresh zero points (frames with NaN ZP carry zero
        # weight in m0/w already — Mv-style zeroing keeps the sums finite).
        ZPv = np.where(np.isfinite(ZP_new), ZP_new, 0.0)
        ws = np.where(np.isfinite(ZP_new)[None, :], w, 0.0)
        denom_s = ws.sum(axis=1)
        with np.errstate(invalid="ignore"):
            M_new = np.where(denom_s > 0,
                             ((m0 - ZPv[None, :]) * ws).sum(axis=1) / denom_s,
                             np.nan)
        # Gauge: mean ZP over frames that actually got one = 0; frames
        # without one stay NaN.
        got = denom_f > 0
        if got.any():
            shift = ZP_new[got].mean()
            ZP_new = np.where(got, ZP_new - shift, np.nan)
            M_new = M_new + shift
        # Convergence over the entries defined in BOTH iterations (a frame
        # whose ZP just appeared or vanished cannot report a delta).
        bz = np.isfinite(ZP_new) & np.isfinite(ZP)
        dz = np.max(np.abs(ZP_new[bz] - ZP[bz])) if bz.any() else 0.0
        both = np.isfinite(M_new) & np.isfinite(M)
        dm = np.max(np.abs(M_new[both] - M[both])) if both.any() else 0.0
        M, ZP = M_new, ZP_new
        if max(dz, dm) < CONVERGE_MMAG:
            converged = True
            break
    return M, ZP, it, converged


def solve_ensemble(mag: np.ndarray, sig: np.ndarray,
                   mask: Optional[np.ndarray] = None) -> EnsembleSolution:
    """Robust Honeycutt solve over a (stars x frames) magnitude matrix.

    Parameters
    ----------
    mag, sig
        (S, F) instrumental magnitudes and their errors; NaN = unobserved.
    mask
        Optional (S, F) bool of points allowed to vote (default: every
        finite mag/sig pair).  The comp-selection wrapper passes the comp
        rows only.

    The plain solve runs first; then :data:`CLIP_PASSES` robust passes mask
    residuals beyond :data:`CLIP_SIGMA` global RMS and re-solve.  Clipped
    points keep their residual value in ``residual`` (computed against the
    final solution) but carry ``clipped=True`` and zero weight.
    """
    mag = np.asarray(mag, dtype=float)
    sig = np.asarray(sig, dtype=float)
    base = np.isfinite(mag) & np.isfinite(sig)
    if mask is not None:
        base &= np.asarray(mask, dtype=bool)

    active = base.copy()
    total_iter = 0
    all_conv = True
    M = ZP = None
    for _pass in range(1 + CLIP_PASSES):
        w = _weights(sig, active)
        M, ZP, it, conv = _alternate(mag, w)
        total_iter += it
        all_conv &= conv
        res = mag - M[:, None] - ZP[None, :]
        res_act = res[active & np.isfinite(res)]
        if res_act.size == 0:
            break
        rms = float(np.sqrt(np.mean(res_act ** 2)))
        bad = active & np.isfinite(res) & (np.abs(res) > CLIP_SIGMA * rms)
        if not bad.any():
            break                      # nothing left to clip — done early
        active &= ~bad

    res = mag - M[:, None] - ZP[None, :]
    w = _weights(sig, active)
    denom_f = np.where(np.isfinite(M)[:, None], w, 0.0).sum(axis=0)
    with np.errstate(divide="ignore"):
        zp_err = np.where(denom_f > 0, 1.0 / np.sqrt(denom_f), np.nan)
    return EnsembleSolution(
        mean_mag=M, zp=ZP, zp_err=zp_err,
        n_star_used=(w > 0).sum(axis=0),
        residual=np.where(base, res, np.nan),
        clipped=base & ~active,
        n_iter=total_iter, converged=all_conv)


# --------------------------------------------------------------------------
# Per-star statistics against a FIXED zero-point solution
# --------------------------------------------------------------------------

def star_stats(mag: np.ndarray, sig: np.ndarray, zp: np.ndarray
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate every star against fixed frame zero points.

    Returns ``(mean_mag, rms, nobs, chi2nu)`` per star, where ``chi2nu`` is
    the reduced chi-square of the constant-star hypothesis using the
    measurement errors (with the weighting floor) — the quantity whose
    median over CHECK stars becomes the error inflation factor.
    """
    mag = np.asarray(mag, dtype=float)
    sig = np.asarray(sig, dtype=float)
    ok = np.isfinite(mag) & np.isfinite(sig) & np.isfinite(zp)[None, :]
    corr = np.where(ok, mag - zp[None, :], np.nan)   # ZP-corrected mags
    var = np.square(sig) + WEIGHT_FLOOR_MAG ** 2
    w = np.where(ok, 1.0 / var, 0.0)
    nobs = ok.sum(axis=1)
    # Written with nansum + explicit counts (not nanmean) so all-NaN rows
    # produce NaN quietly instead of a RuntimeWarning per never-seen star.
    with np.errstate(invalid="ignore"):
        mean = np.where(nobs > 0,
                        np.nansum(corr * w, axis=1) / w.sum(axis=1), np.nan)
        dev = corr - mean[:, None]
        rms = np.where(nobs > 0,
                       np.sqrt(np.nansum(np.square(dev), axis=1)
                               / np.maximum(nobs, 1)), np.nan)
        chi2 = np.nansum(np.square(dev) / np.where(ok, var, np.inf), axis=1)
        chi2nu = np.where(nobs > 1, chi2 / np.maximum(nobs - 1, 1), np.nan)
    return mean, rms, nobs, chi2nu


# --------------------------------------------------------------------------
# Comparison-star selection by stability iteration
# --------------------------------------------------------------------------

@dataclass
class CompSelection:
    """Outcome of the stability iteration for one (target, era, filter).

    ``role`` per star: 'comp' (votes in the ZP), 'check' (stable but held
    out to validate), 'dropped_unstable' (was a candidate, failed the RMS
    cut — the polar target lands here by design), 'field' (never a
    candidate: too few frames or excluded up front).
    """
    role: np.ndarray            # (S,) unicode roles
    solution: EnsembleSolution  # final solve on the comp rows only
    n_passes: int
    comp_rms_median: float


def select_comps(mag: np.ndarray, sig: np.ndarray,
                 exclude: Optional[np.ndarray] = None,
                 target_row: Optional[int] = None,
                 n_check: int = N_CHECK_STARS) -> CompSelection:
    """Choose comparison stars by iterative stability, then hold out checks.

    Parameters
    ----------
    mag, sig
        (S, F) instrumental magnitudes / errors, NaN = unobserved.
    exclude
        (S,) bool — stars barred from comp duty regardless of stability
        (the TARGET is always excluded: a polar's orbital modulation must
        never help set the zero point).
    target_row
        Row index of the target, used both to auto-exclude it and to pick
        check stars near its magnitude.
    n_check
        Stable stars to hold out as checks (nearest the target's median
        magnitude — they sample the error model where the science lives).

    The iteration: start from every star seen on >= COMP_MIN_FRAME_FRAC of
    frames (minus exclusions); solve; drop comps whose residual RMS exceeds
    max(COMP_RMS_FACTOR x median comp RMS, COMP_RMS_FLOOR_MAG); repeat
    until stable or COMP_MAX_PASSES.  Variables (the target, any field
    variable) fail the RMS cut and are dropped — that is the design, not a
    failure mode.  After the comp set settles, the n_check stable stars
    closest to the target's magnitude move comp -> check and the final
    solve runs WITHOUT them, so their statistics are honest hold-outs.
    """
    mag = np.asarray(mag, dtype=float)
    sig = np.asarray(sig, dtype=float)
    S, F = mag.shape
    seen = np.isfinite(mag) & np.isfinite(sig)
    candidate = seen.sum(axis=1) >= COMP_MIN_FRAME_FRAC * F
    barred = np.zeros(S, dtype=bool)
    if exclude is not None:
        barred |= np.asarray(exclude, dtype=bool)
    if target_row is not None:
        barred[target_row] = True
    comp = candidate & ~barred

    n_passes = 0
    med = float("nan")
    for n_passes in range(1, COMP_MAX_PASSES + 1):
        if comp.sum() < MIN_COMPS:
            break                       # too few comps — stop shrinking
        sol = solve_ensemble(mag, sig, mask=comp[:, None] & seen)
        # Per-star residual RMS, written to avoid the all-NaN-row warning
        # (a star with no residuals — e.g. never observed — gets NaN RMS).
        n_res = np.isfinite(sol.residual).sum(axis=1)
        sq = np.nansum(np.square(sol.residual), axis=1)
        with np.errstate(invalid="ignore"):
            rms = np.where(n_res > 0, np.sqrt(sq / np.maximum(n_res, 1)),
                           np.nan)
        med = float(np.nanmedian(rms[comp]))
        cut = max(COMP_RMS_FACTOR * med, COMP_RMS_FLOOR_MAG)
        drop = comp & (rms > cut)
        # Never drop below MIN_COMPS: keep the quietest if the cut is too
        # eager (deterministic: sort by rms, keep the smallest).
        if (comp.sum() - drop.sum()) < MIN_COMPS:
            order = np.argsort(np.where(comp, rms, np.inf), kind="stable")
            keep = set(order[:MIN_COMPS].tolist())
            drop = comp.copy()
            drop[list(keep)] = False
        if not drop.any():
            break
        comp &= ~drop

    # ---- hold out check stars nearest the target's magnitude -------------
    role = np.array(["field"] * S, dtype=object)
    role[candidate & ~barred & ~comp] = "dropped_unstable"
    role[comp] = "comp"
    # Per-star median magnitude, all-NaN rows answered with NaN silently
    # (np.nanmedian warns per empty row; never-seen stars are routine here).
    med_mag = np.full(S, np.nan)
    for i in range(S):
        finite = mag[i][np.isfinite(mag[i])]
        if finite.size:
            med_mag[i] = float(np.median(finite))
    if target_row is not None and np.isfinite(med_mag[target_row]):
        anchor = med_mag[target_row]
    else:
        anchor = float(np.nanmedian(med_mag[comp])) if comp.any() else 0.0
    comp_idx = np.flatnonzero(comp)
    # Hold out checks only while MIN_COMPS comps survive.
    n_hold = min(n_check, max(0, len(comp_idx) - MIN_COMPS))
    if n_hold > 0:
        by_dist = comp_idx[np.argsort(np.abs(med_mag[comp_idx] - anchor),
                                      kind="stable")]
        for idx in by_dist[:n_hold]:
            role[idx] = "check"
            comp[idx] = False
    if target_row is not None:
        role[target_row] = "target"

    final = solve_ensemble(mag, sig, mask=comp[:, None] & seen)
    return CompSelection(role=role.astype(str), solution=final,
                         n_passes=n_passes, comp_rms_median=med)
