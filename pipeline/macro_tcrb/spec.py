"""macro_tcrb.spec — Hα equivalent width, LSF degradation, line flux.

Pure functions implementing the measurement definitions PRE-REGISTERED in
``TCrB_Monitoring/ANALYSIS_STRATEGY.md`` §10 (TCRB-A5a block, fixed
2026-10-05T02:40Z; tolerance 03:35Z):

* EW is emission-POSITIVE, in Å, on the frame's wavelength scale.
* Line window: 6562.8 Å ± h, h = max(30 Å, 1.5 × LSF FWHM).
* Pseudo-continuum: a straight line through (median λ, median flux) of
  6470–6520 Å and 6600–6640 Å.
* Sensitivity: every window edge shifted by −10 and +10 Å.
* ARAS spectra are degraded to the RLMT LSF (Gaussian, quadrature
  difference of FWHMs) and measured with the identical windows.

Tests: ``pipeline/tests/test_tcrb.py``.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

HALPHA = 6562.8
BLUE = (6470.0, 6520.0)
RED = (6600.0, 6640.0)
H_MIN = 30.0
H_LSF = 1.5
SHIFT = 10.0


def half_width(lsf_fwhm_a: Optional[float]) -> float:
    """Pre-registered line half-width h (Å)."""
    if lsf_fwhm_a is None or not np.isfinite(lsf_fwhm_a):
        return H_MIN
    return float(max(H_MIN, H_LSF * lsf_fwhm_a))


def equivalent_width(wave: np.ndarray, flux: np.ndarray,
                     var: Optional[np.ndarray] = None, h: float = H_MIN,
                     shift: float = 0.0, centre: float = HALPHA,
                     blue: tuple = BLUE, red: tuple = RED) -> dict:
    """EW (emission positive) with the pre-registered windows.

    ``shift`` moves EVERY window edge by the same amount (the ±10 Å test).
    Returns ew, ew_err (from ``var`` and the continuum-window scatter),
    the continuum at Hα (flux units) and the number of pixels used.
    """
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    v = None if var is None else np.asarray(var, float)
    ok = np.isfinite(w) & np.isfinite(f)
    if v is not None:
        ok &= np.isfinite(v)
    w, f = w[ok], f[ok]
    v = v[ok] if v is not None else None
    if w.size < 10 or w[0] > w[-1]:
        o = np.argsort(w)
        w, f = w[o], f[o]
        v = v[o] if v is not None else None

    def win(lo, hi):
        return (w >= lo + shift) & (w <= hi + shift)
    b, r = win(*blue), win(*red)
    ln = win(centre - h, centre + h)
    nan = float("nan")
    if b.sum() < 3 or r.sum() < 3 or ln.sum() < 3:
        return {"ew": nan, "ew_err": nan, "cont": nan, "n_line": int(ln.sum())}
    xb, yb = np.median(w[b]), np.median(f[b])
    xr, yr = np.median(w[r]), np.median(f[r])
    slope = (yr - yb) / (xr - xb)
    cont = yb + slope * (w - xb)
    if np.any(cont[ln] <= 0):
        return {"ew": nan, "ew_err": nan, "cont": nan, "n_line": int(ln.sum())}
    dl = np.gradient(w)
    ew = float(np.sum((f[ln] / cont[ln] - 1.0) * dl[ln]))
    # Error: pixel noise inside the line window + continuum-level error
    # (median of a window: 1.2533 sigma/sqrt(n), sigma from the window's
    # own robust scatter about a line), propagated to first order.
    var_line = (np.sum((dl[ln] / cont[ln]) ** 2 * v[ln])
                if v is not None else 0.0)

    def med_err(sel):
        x, y = w[sel], f[sel]
        p = np.polyfit(x, y, 1)
        s = 1.4826 * np.median(np.abs(y - np.polyval(p, x)))
        return 1.2533 * s / np.sqrt(sel.sum())
    c_ha = yb + slope * (centre + shift - xb)
    ec = float(np.hypot(med_err(b), med_err(r)) / np.sqrt(2.0))
    # dEW/dC ~ -(EW + width)/C for a uniform continuum shift
    width = float(np.sum(dl[ln]))
    var_cont = ((ew + width) / c_ha * ec) ** 2
    return {"ew": ew, "ew_err": float(np.sqrt(var_line + var_cont)),
            "cont": float(c_ha), "n_line": int(ln.sum())}


def degrade(wave: np.ndarray, flux: np.ndarray, r_native: Optional[float],
            fwhm_target_a: float) -> np.ndarray:
    """Convolve a spectrum to a Gaussian LSF of FWHM ``fwhm_target_a`` (Å).

    The kernel FWHM is the quadrature difference between the target and the
    native resolution at Hα (λ/R); if the native LSF is already broader the
    spectrum is returned unchanged.  Assumes a (near-)linear wavelength
    grid, as every ARAS file is.
    """
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    fn = HALPHA / r_native if r_native and r_native > 0 else 0.0
    if fwhm_target_a <= fn:
        return f.copy()
    k_fwhm = np.sqrt(fwhm_target_a ** 2 - fn ** 2)
    dl = float(np.median(np.diff(w)))
    sig = k_fwhm / 2.3548 / dl
    half = int(np.ceil(5 * sig))
    x = np.arange(-half, half + 1)
    ker = np.exp(-0.5 * (x / sig) ** 2)
    ker /= ker.sum()
    good = np.isfinite(f)
    fz = np.where(good, f, 0.0)
    num = np.convolve(fz, ker, mode="same")
    den = np.convolve(good.astype(float), ker, mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den
    out[den < 0.5] = np.nan
    return out


def line_flux(ew: float, ew_err: float, cont_flux: float,
              cont_err: float) -> tuple[float, float]:
    """F(Hα) = EW × F_λ(continuum at Hα); first-order error."""
    f = ew * cont_flux
    e = float(np.hypot(ew_err * cont_flux, ew * cont_err))
    return float(f), e


def mag_to_flambda(mag: float, zero_flux: float) -> float:
    """F_λ from a magnitude given the band's zero-magnitude flux density."""
    return float(zero_flux * 10 ** (-0.4 * mag))


def orbital_phase(jd: float, t0: float, period: float) -> float:
    """Phase in [0, 1) on the given ephemeris (T0 = giant at max velocity,
    as in Munari et al. 2025)."""
    return float(((jd - t0) / period) % 1.0)


def wavelength_scale(x: np.ndarray, coeffs: list, x_ref: float,
                     x_anchor: float, lam_anchor: float = HALPHA
                     ) -> np.ndarray:
    """λ(x) on the G-1 fixed dispersion with a per-frame zero point.

    G-1 stores λ = c_frame + Σ_j k_j ((x − x_ref)/1000)^j (j = 1..deg);
    the zero point c_frame is set so that λ(x_anchor) = lam_anchor (the
    frame's own Hα centroid, G-1/A3's anchor).  Only the zero point is a
    per-frame quantity (U2).
    """
    def poly(xx):
        u = (np.asarray(xx, float) - x_ref) / 1000.0
        return sum(k * u ** (j + 1) for j, k in enumerate(coeffs))
    c = lam_anchor - poly(x_anchor)
    return c + poly(x)


# ---------------------------------------------------------------------------
# Calibrator (theta CrB) helpers: alignment and feature-pair zero point
# ---------------------------------------------------------------------------
def normalise(flux: np.ndarray, win: int = 301) -> np.ndarray:
    """Flux divided by a running median (continuum-normalised shape)."""
    from scipy.ndimage import median_filter
    f = np.asarray(flux, float)
    good = np.isfinite(f)
    g = np.where(good, f, np.nanmedian(f))
    c = median_filter(g, size=win, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = g / c
    out[~good | (c <= 0)] = np.nan
    return out


def xcorr_lag(ref: np.ndarray, spec: np.ndarray, max_lag: int = 60
              ) -> float:
    """Sub-pixel lag (spec relative to ref) maximising the correlation of
    the two normalised, mean-removed shapes; parabolic peak interpolation."""
    a = np.nan_to_num(np.asarray(ref, float) - np.nanmean(ref))
    b = np.nan_to_num(np.asarray(spec, float) - np.nanmean(spec))
    lags = np.arange(-max_lag, max_lag + 1)
    cc = np.array([np.sum(a[max(0, -l):len(a) - max(0, l)] *
                          b[max(0, l):len(b) - max(0, -l)]) for l in lags])
    i = int(np.argmax(cc))
    if 0 < i < len(cc) - 1:
        y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
        den = y0 - 2 * y1 + y2
        off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    else:
        off = 0.0
    return float(lags[i] + off)


def absorption_minima(norm: np.ndarray, depth: float = 0.05,
                      sep: int = 15) -> np.ndarray:
    """Pixel positions of local minima deeper than ``depth`` below 1."""
    from scipy.signal import find_peaks
    x = np.nan_to_num(1.0 - np.asarray(norm, float))
    pk, _ = find_peaks(x, height=depth, distance=sep)
    return pk


def pair_zero_point(minima: np.ndarray, disp: float, lam_a: float,
                    lam_b: float, tol_px: float = 4.0) -> Optional[float]:
    """The pixel of feature A, from the one pair of minima whose
    separation matches (lam_b - lam_a) / disp within ``tol_px``.

    ``disp`` is signed (Å/px), so the expected pixel offset of B from A is
    (lam_b - lam_a) / disp.  Returns None when no pair, or more than one
    pair, matches (an ambiguous identification is not an identification).
    """
    want = (lam_b - lam_a) / disp
    hits = [(abs((b - a) - want), a) for a in minima for b in minima
            if abs((b - a) - want) <= tol_px]
    if len(hits) != 1:
        return None
    return float(hits[0][1])


def pair_zero_point_poly(norm: np.ndarray, minima: np.ndarray, coeffs: list,
                         x_ref: float, lam_a: float, lam_b: float,
                         tol_px: float = 8.0,
                         valid: Optional[np.ndarray] = None,
                         minima_b: Optional[np.ndarray] = None,
                         prior: Optional[tuple] = None
                         ) -> Optional[float]:
    """Feature-A pixel from an A/B absorption pair on the FULL G-1
    polynomial: for each candidate minimum a (inside ``valid``), anchor
    λ(a) = lam_a, predict where λ = lam_b falls, and require a minimum
    within ``tol_px`` there.  Among the candidates that pass, the pair with
    the largest summed depth wins; ``prior`` = (pixel, half-width) restricts
    A to a window (the same grism's Hα position on the target frames).
    Returns None when no candidate passes."""
    x = np.arange(norm.size, dtype=float)
    best = None
    for a in minima:
        if valid is not None and not valid[a]:
            continue
        lam = wavelength_scale(x, coeffs, x_ref, float(a), lam_a)
        o = np.argsort(lam)
        xb = float(np.interp(lam_b, lam[o], x[o]))
        if not (0 <= xb < norm.size):
            continue
        mb = minima if minima_b is None else minima_b
        near = [b for b in mb if abs(b - xb) <= tol_px and
                (valid is None or valid[b])]
        if not near:
            continue
        b = min(near, key=lambda bb: abs(bb - xb))
        if prior is not None and abs(a - prior[0]) > prior[1]:
            continue
        depth = (1.0 - norm[a]) + (1.0 - norm[b])
        if best is None or depth > best[1]:
            best = (float(a), depth)
    return None if best is None else best[0]
