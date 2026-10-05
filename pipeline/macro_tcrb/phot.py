"""macro_tcrb.phot — Phase B (B-only anchors) and C1 decision logic.

Pure functions; the I/O lives in ``pipeline/scripts/run_tcrb.py``
(stages ``measure``, ``filters``, ``zmag``, ``census``, ``ensemble``,
``errors``, ``flicker``).  Tests: ``pipeline/tests/test_tcrb.py``.

Conventions
-----------
* Instrumental magnitude ``m_inst = -2.5 log10(flux_ADU / t_s)``.
* Zero-point model per frame: ``m_cat - m_inst = zp + k (c - c_ref)``,
  with c a catalogue colour; fitted by Huber-weighted IRLS and a final
  3-sigma clip (:func:`huber_line`).  Unlike ``macro_phot.cattie``'s
  fitter, NO ``max(chi2nu, 1)`` inflation is applied (SYNTHESIS standing
  rule 1): the formal errors are reported with chi2nu and its dof beside
  them, and the EMPIRICAL errors (check-star rms, B5) are what the paper
  quotes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

# ===========================================================================
# Calibration masters (B1)
# ===========================================================================
def pick_master(cands: Sequence[dict], kind: str, mode: str,
                exptime: float, night: str, code: Optional[str] = None
                ) -> Optional[dict]:
    """The master of ``kind`` ('dark' | 'flat') for one frame.

    Mode-matched only (ruling 7: no master crosses readout modes); a flat
    must also match the filter code.  Darks: nearest exposure in log time,
    then nearest night.  Flats: nearest night.  Candidates are dicts with
    keys kind, mode, code, exptime, night, path.  None when no mode-matched
    master exists — the caller then records the frame as uncalibrated for
    that step and says so.
    """
    def nd(n):  # night distance in days (ISO strings)
        from datetime import date
        a = date.fromisoformat(n[:10])
        b = date.fromisoformat(night[:10])
        return abs((a - b).days)
    pool = [c for c in cands if c["kind"] == kind and c["mode"] == mode
            and (kind != "flat" or c["code"] == code)]
    if not pool:
        return None
    if kind == "dark":
        return min(pool, key=lambda c: (round(abs(math.log(
            max(c["exptime"], 1e-3) / max(exptime, 1e-3))), 2), nd(c["night"])))
    return min(pool, key=lambda c: nd(c["night"]))


# ===========================================================================
# Robust zero-point / colour-term fit
# ===========================================================================
@dataclass
class Fit:
    zp: float
    zp_err: float
    k: float
    k_err: float
    c_ref: float
    n: int
    n_clip: int
    chi2: float
    dof: int
    rms: float
    used: np.ndarray


def huber_line(x: np.ndarray, y: np.ndarray, s: np.ndarray,
               c_ref: Optional[float] = None, delta: float = 1.5,
               clip: float = 3.0, n_iter: int = 30,
               fit_slope: bool = True) -> Optional[Fit]:
    """Huber IRLS fit of y = zp + k (x - c_ref); final 3-sigma clip
    (both on the MAD scale of the normalised residuals).

    ``s`` are per-point errors.  Returns None for fewer than 3 usable
    points.  Formal errors come from the weighted normal matrix with NO
    chi2 rescaling; chi2 and dof are returned so the caller reports them.
    """
    x, y, s = (np.asarray(a, float) for a in (x, y, s))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    if ok.sum() < 3:
        return None
    cr = float(np.median(x[ok])) if c_ref is None else float(c_ref)
    X = np.vstack([np.ones_like(x), x - cr]).T if fit_slope else \
        np.ones((x.size, 1))
    w = np.where(ok, 1.0 / s ** 2, 0.0)
    beta = np.zeros(X.shape[1])
    def scale(r):
        # Robust scale of the NORMALISED residuals: the Huber threshold and
        # the clip are set on it, so mis-stated per-point errors cannot
        # make the clip delete good points (or keep bad ones).  It never
        # touches the reported errors.
        rr = r[ok]
        return max(1.4826 * np.median(np.abs(rr - np.median(rr))), 1e-9)

    for _ in range(n_iter):
        W = w.copy()
        r = (y - X @ beta) / s
        a = np.abs(r) / scale(r)
        W *= np.where(a <= delta, 1.0, delta / np.maximum(a, 1e-12))
        W[~ok] = 0
        A = X.T @ (X * W[:, None])
        new = np.linalg.solve(A, X.T @ (W * np.where(ok, y, 0)))
        if np.allclose(new, beta, atol=1e-7):
            beta = new
            break
        beta = new
    r = (y - X @ beta) / s
    keep = ok & (np.abs(r - np.median(r[ok])) <= clip * scale(r))
    if keep.sum() < 3:
        keep = ok
    W = np.where(keep, 1.0 / s ** 2, 0.0)
    A = X.T @ (X * W[:, None])
    beta = np.linalg.solve(A, X.T @ (W * np.where(keep, y, 0)))
    cov = np.linalg.inv(A)
    res = (y - X @ beta)[keep]
    chi2 = float(np.sum((res / s[keep]) ** 2))
    dof = int(keep.sum() - X.shape[1])
    k = float(beta[1]) if fit_slope else 0.0
    ke = float(np.sqrt(cov[1, 1])) if fit_slope else 0.0
    return Fit(float(beta[0]), float(np.sqrt(cov[0, 0])), k, ke, cr,
               int(keep.sum()), int(ok.sum() - keep.sum()), chi2, dof,
               float(np.sqrt(np.mean(res ** 2))), keep)


# ===========================================================================
# Saturation verdict at the target (B0) — standing rule 4
# ===========================================================================
def peak_verdict(peak_native: Optional[float], cap_adu: Optional[float],
                 clip_adu: Optional[float], dispersed: bool) -> str:
    """'dispersed' | 'clipped' | 'above_cap' | 'clean' | 'unknown'.

    ``peak_native`` is the RAW native-pixel peak at the target; ``cap_adu``
    the mode's linearity cap and ``clip_adu`` its ceiling (both raw ADU,
    S2 ``detector_params`` / ``s2_target_verdicts``).
    """
    if dispersed:
        return "dispersed"
    if peak_native is None or cap_adu is None or not np.isfinite(peak_native):
        return "unknown"
    if clip_adu is not None and peak_native >= clip_adu - 1:
        return "clipped"
    if peak_native > cap_adu:
        return "above_cap"
    return "clean"


# ===========================================================================
# Empirical errors (B5) and the flickering limit (C1)
# ===========================================================================
def check_star_rms(resid: np.ndarray) -> float:
    """Robust rms of a check star's residuals about its own mean
    (1.4826 MAD; NaN for fewer than 3 points)."""
    r = np.asarray(resid, float)
    r = r[np.isfinite(r)]
    if r.size < 3:
        return float("nan")
    return float(1.4826 * np.median(np.abs(r - np.median(r))))


def detrend_model(n: int, n_cov: int) -> tuple[int, int]:
    """(polynomial order in time, number of covariates) for a snippet of
    ``n`` frames, keeping >= 3 degrees of freedom: a mean only below 6
    frames, a linear trend from 6, covariates (airmass, FWHM) from 10."""
    if n < 6:
        return 0, 0
    if n < 10:
        return 1, 0
    return 1, min(n_cov, n - 5)


def detrend_rms(t: np.ndarray, m: np.ndarray, X: Optional[np.ndarray] = None,
                order: int = 1) -> tuple[float, int]:
    """rms of m after a joint linear fit: polynomial of ``order`` (0 or 1)
    in time plus the covariates ``X``.  Returns (rms, dof)."""
    t, m = np.asarray(t, float), np.asarray(m, float)
    cols = [np.ones_like(t)] + ([t - t.mean()] if order >= 1 else [])
    if X is not None and X.size:
        for c in np.atleast_2d(X.T if X.ndim == 2 else X[None, :]):
            c = np.asarray(c, float)
            if np.nanstd(c) > 0:
                cols.append((c - c.mean()) / c.std())
    A = np.vstack(cols).T
    coef, *_ = np.linalg.lstsq(A, m, rcond=None)
    r = m - A @ coef
    dof = max(t.size - A.shape[1], 1)
    return float(np.sqrt(np.sum(r ** 2) / dof)), int(dof)


def red_noise(n: int, dt_s: np.ndarray, sigma: float, beta: float,
              rng: np.random.Generator) -> np.ndarray:
    """Power-law (P ~ f^-beta) noise sampled at irregular times.

    Generated on a fine uniform grid by spectral synthesis (Timmer & Koenig
    1995 amplitudes) and interpolated to the sample times; normalised to
    rms ``sigma`` over the samples.  ``dt_s`` are the sample times (s).
    """
    t = np.asarray(dt_s, float) - float(np.min(dt_s))
    span = max(t.max(), 1.0)
    m = 4096
    f = np.fft.rfftfreq(m, d=4 * span / m)[1:]
    amp = f ** (-beta / 2.0)
    ph = rng.normal(size=f.size) + 1j * rng.normal(size=f.size)
    spec = np.concatenate([[0], amp * ph])
    x = np.fft.irfft(spec, n=m)
    grid = np.linspace(0, 4 * span, m)
    off = rng.uniform(0, 3 * span)
    y = np.interp(t + off, grid, x)
    y = y - y.mean()
    sd = y.std()
    return y * (sigma / sd if sd > 0 else 0.0)


def flicker_upper_limit(t_s: np.ndarray, resid: np.ndarray,
                        noise_resid: Sequence[np.ndarray],
                        X: Optional[np.ndarray] = None,
                        betas: Sequence[float] = (1.0, 2.0),
                        amps: Optional[np.ndarray] = None, n_mc: int = 300,
                        level: float = 0.95, seed: int = 1) -> dict:
    """95% upper limit on the flickering rms of one snippet.

    ``resid``: the target's differential magnitudes (comparison-ensemble
    corrected).  ``noise_resid``: check-star residual series of the SAME
    frames (the empirical noise, not a Poisson model).  For each injected
    rms A the injected series is red noise of index beta added to a
    check-star series drawn at random; the detrended rms of each injection
    is compared with the target's observed detrended rms.  The limit is the
    smallest A for which >= 95% of injections exceed the observed rms (an
    rms larger than the observed one would have been seen).  The detection
    statistic is also returned: ``recovered90`` is the injected rms that
    exceeds the pure-noise 95% threshold in 90% of injections — what the
    snippet WOULD have detected (standing rule 2).  When the observed rms
    falls below the noise (few frames), the data-conditioned limit can be
    0; the quoted limit is then the larger of the two (the caller's
    choice, stated in the table).
    """
    rng = np.random.default_rng(seed)
    order, ncov = detrend_model(len(resid), 0 if X is None else
                                np.atleast_2d(X).shape[1])
    X = None if ncov == 0 else np.atleast_2d(X)[:, :ncov]
    obs, dof = detrend_rms(t_s, resid, X, order)
    noise = [np.asarray(n, float) for n in noise_resid
             if np.isfinite(n).all() and len(n) == len(resid)]
    if not noise:
        return {"obs_rms": obs, "dof": dof, "limit": {}, "noise_rms": None,
                "order": order, "n_cov": ncov}
    nrms = np.array([detrend_rms(t_s, n, X, order)[0] for n in noise])
    if amps is None:
        amps = np.linspace(0.0, 0.60, 241)
    # Pure-noise detection threshold: the 95th percentile of the detrended
    # rms of check-star series alone (resampled with random sign flips so
    # the few check stars give a distribution, not three numbers).
    thr = np.percentile([detrend_rms(t_s, noise[rng.integers(len(noise))]
                                     * rng.choice([-1, 1], len(t_s)), X,
                                     order)[0] for _ in range(n_mc)], 95)
    lim, rec = {}, {}
    for beta in betas:
        found = rfound = float("nan")
        for a in amps:
            cnt = det = 0
            for _ in range(n_mc):
                base = noise[rng.integers(len(noise))]
                y = base + red_noise(len(t_s), t_s, a, beta, rng)
                r = detrend_rms(t_s, y, X, order)[0]
                cnt += r > obs
                det += r > thr
            if np.isnan(found) and cnt / n_mc >= level:
                found = float(a)
            if np.isnan(rfound) and det / n_mc >= 0.90:
                rfound = float(a)
            if np.isfinite(found) and np.isfinite(rfound):
                break
        lim[beta], rec[beta] = found, rfound
    return {"obs_rms": obs, "dof": dof, "limit": lim, "recovered90": rec,
            "noise_threshold95": float(thr), "order": order,
            "n_cov": ncov, "noise_rms_median": float(np.median(nrms)),
            "noise_rms_p95": float(np.percentile(nrms, 95))}


# ===========================================================================
# C4 — uncertainties: red-noise inflation and posteriors
# ===========================================================================
def red_noise_beta(series: Sequence[np.ndarray], max_bin: int = 3) -> float:
    """Red-noise inflation factor from binned comparison-star residuals.

    For each check-star series (one night), the rms of means of n
    consecutive points is compared with the white-noise expectation
    sigma_1 / sqrt(n); beta is the median ratio over series and n = 2..max_bin
    (Pont, Zucker & Queloz 2006 time-averaging method).  1 = white.
    """
    ratios = []
    for s in series:
        s = np.asarray(s, float)
        s = s[np.isfinite(s)] - np.nanmean(s)
        if s.size < 4:
            continue
        s1 = s.std(ddof=1)
        for n in range(2, max_bin + 1):
            m = s.size // n
            if m < 2:
                continue
            b = s[:m * n].reshape(m, n).mean(axis=1)
            exp = s1 / np.sqrt(n) * np.sqrt(m / (m - 1))
            if exp > 0:
                ratios.append(b.std(ddof=1) / exp)
    return float(max(np.median(ratios), 1.0)) if ratios else float("nan")


def mean_with_jitter(y: np.ndarray, s: np.ndarray, n_walkers: int = 24,
                     n_steps: int = 3000, seed: int = 11) -> dict:
    """emcee posterior of a mean with an unknown extra scatter (jitter).

    y_i ~ N(mu, s_i^2 + j^2), flat prior on mu, log-uniform-ish prior on
    j >= 0 (uniform in log j over 1e-5..10).  Returns the 16/50/84
    percentiles of mu and j, and the acceptance fraction.
    """
    import emcee
    y, s = np.asarray(y, float), np.asarray(s, float)

    def lnp(th):
        mu, lj = th
        if not -5 < lj < 1:
            return -np.inf
        v = s ** 2 + (10 ** lj) ** 2
        return -0.5 * np.sum((y - mu) ** 2 / v + np.log(v))
    rng = np.random.default_rng(seed)
    p0 = np.c_[np.median(y) + 1e-3 * rng.standard_normal(n_walkers),
               np.log10(np.std(y) + 1e-4) + 0.1 * rng.standard_normal(
                   n_walkers)]
    sm = emcee.EnsembleSampler(n_walkers, 2, lnp)
    sm.random_state = np.random.RandomState(seed).get_state()
    sm.run_mcmc(p0, n_steps, progress=False)
    ch = sm.get_chain(discard=n_steps // 3, flat=True)
    mu = np.percentile(ch[:, 0], [16, 50, 84])
    jj = np.percentile(10 ** ch[:, 1], [16, 50, 84])
    return {"mu": mu.tolist(), "jitter": jj.tolist(),
            "acc": float(np.mean(sm.acceptance_fraction))}
