"""macro_sn.snphot — the decision rules of the SN 2023ixf release.

Pure functions only (numpy/scipy, no file or database access), unit-tested
in ``pipeline/tests/test_sn_phot.py``.  The plumbing is
``pipeline/scripts/run_sn_photometry.py``; the pixels are
``macro_sn.snio``.  Every threshold is a module constant with its
provenance, so the paper's methods section, the tests and the code read the
same number.
"""

from __future__ import annotations

import math
import warnings
import zlib
from typing import Optional

import numpy as np

# ===========================================================================
# 1.  ASTROMETRY ACCEPTANCE  (SN-S4-resolve)
# ===========================================================================

#: A re-fitted WCS is accepted only with at least this many REFCAT2 stars
#: matched within 2.5 px ...
WCS_MIN_MATCH = 12
#: ... and a sky residual no worse than this many pixels (rms).  At the
#: campaign's 0.54"/px that is 0.81", a quarter of the median FWHM: forced
#: apertures of radius 1.5 x FWHM then lose < 0.1% of a star's flux.
WCS_MAX_RMS_PX = 1.5


#: A translation-only solution under the camera's measured CD matrix has
#: two free parameters, so it is accepted on fewer stars — but no fewer than
#: the calibration's own MIN_ENSEMBLE, so a frame solved this way can still
#: be calibrated — and on a tighter residual.
WCS_MIN_MATCH_SHIFT = 8
WCS_MAX_RMS_PX_SHIFT = 1.0


def shift_acceptable(n_match, rms_arcsec, scale_arcsec) -> bool:
    """Fail-closed acceptance of a fixed-CD translation solution."""
    if n_match is None or rms_arcsec is None or not np.isfinite(rms_arcsec):
        return False
    return (n_match >= WCS_MIN_MATCH_SHIFT
            and rms_arcsec <= WCS_MAX_RMS_PX_SHIFT * scale_arcsec)


def wcs_acceptable(n_match: Optional[int], rms_arcsec: Optional[float],
                   scale_arcsec: float) -> bool:
    """Fail-closed acceptance of a fitted plate solution."""
    if n_match is None or rms_arcsec is None or not np.isfinite(rms_arcsec):
        return False
    return n_match >= WCS_MIN_MATCH and rms_arcsec <= WCS_MAX_RMS_PX * scale_arcsec


# ===========================================================================
# 2.  WHICH REFCAT2 STARS MAY CALIBRATE  (SN-S4-ensemble-cal)
# ===========================================================================

#: M101's nucleus (SIMBAD "M 101": 14 03 12.58 +54 20 55.5, J2000) and the
#: radius inside which a REFCAT2 entry must PROVE it is a star.  REFCAT2 is
#: Gaia-selected, and Gaia catalogues the brightest H II knots of M101's
#: disc as point sources; inside 15' (about M101's D25 radius, 14.4') an
#: entry is kept only with a >= 5-sigma proper motion or parallax.
M101_RA, M101_DEC = 210.802417, 54.348750
DISC_RADIUS_DEG = 15.0 / 60.0
STAR_PM_SIGMA = 5.0

#: Calibration-star window — the strategy's §4 Step 4 values: 12 < r < 16
#: for the broadband ensemble (extended to 17 so the 4-8 s frames, where
#: r < 13 stars clip, keep an ensemble), 0.2 < g - r < 1.2, isolated (no
#: neighbour brighter than 0.1 x the star inside REFCAT2's rp1 >= 5").
CAL_R_RANGE = (11.5, 17.0)
CAL_GR_RANGE = (0.2, 1.2)
CAL_RP1_MIN_ARCSEC = 5.0

#: The colour used for colour terms, and its pivot (a typical G star).
COLOUR = ("gmag", "imag")
COLOUR_PIVOT = 0.8

#: One star in CHECK_MODULUS is held out of every fit (deterministically,
#: by a CRC of its REFCAT2 position, so the split never depends on which
#: frames happened to be measured).  These are the held-out check stars of
#: DS's amendment: the chi2_nu claim is made on them alone.
CHECK_MODULUS = 4


def is_check_star(ra: float, dec: float) -> bool:
    """Deterministic held-out assignment from the star's catalogue position."""
    key = f"{ra:.6f},{dec:.6f}".encode()
    return zlib.crc32(key) % CHECK_MODULUS == 0


def ang_sep_deg(ra1, dec1, ra2, dec2):
    """Great-circle separation in degrees (haversine)."""
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    a = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) \
        * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))


def is_stellar(pmra, e_pmra, pmde, e_pmde, plx, e_plx) -> np.ndarray:
    """A >= STAR_PM_SIGMA proper motion or parallax (NaN = not proven)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        pm_sig = np.hypot(pmra / e_pmra, pmde / e_pmde)
        plx_sig = plx / e_plx
    return (np.nan_to_num(pm_sig) >= STAR_PM_SIGMA) | \
        (np.nan_to_num(plx_sig) >= STAR_PM_SIGMA)


def calib_star_mask(cat: dict) -> np.ndarray:
    """REFCAT2 rows eligible to calibrate (before any per-frame screen)."""
    r, g = cat["rmag"], cat["gmag"]
    gr = g - r
    in_disc = ang_sep_deg(cat["RA_ICRS"], cat["DE_ICRS"], M101_RA,
                          M101_DEC) < DISC_RADIUS_DEG
    stellar = is_stellar(cat["pmRA"], cat["e_pmRA"], cat["pmDE"],
                         cat["e_pmDE"], cat["Plx"], cat["e_Plx"])
    dv = np.nan_to_num(cat["dupvar"], nan=2).astype(int)
    ok_flag = (dv % 2 == 0) & (dv < 4)       # not Gaia-variable, not duplicated
    rp1 = np.nan_to_num(cat["rp1"], nan=99.9)
    with np.errstate(invalid="ignore"):
        return ((r > CAL_R_RANGE[0]) & (r < CAL_R_RANGE[1])
                & (gr > CAL_GR_RANGE[0]) & (gr < CAL_GR_RANGE[1])
                & np.isfinite(cat["imag"]) & (rp1 >= CAL_RP1_MIN_ARCSEC)
                & ok_flag & (~in_disc | stellar))


# ===========================================================================
# 3.  THE ERROR MODEL  (SN-S5: an explicit scintillation term)
# ===========================================================================

#: Telescope aperture (m) and site altitude (m): RLMT 0.508 m (FITS APTDIA
#: 508 mm) at Winer Observatory (macro_core.timing.WINER_ALT_M = 1515 m).
APERTURE_M = 0.508
SITE_ALT_M = 1515.0
#: Young (1967) scintillation, in the form of Osborn et al. (2015, MNRAS
#: 452, 1707), eq. 7:  sigma_S = 0.09 D^(-2/3) X^(7/4) exp(-h/H) (2t)^(-1/2)
#: with D in cm, H = 8000 m.  The coefficient 0.09 is Young's; Osborn et al.
#: find the median empirical value ~1.5x higher, which is why the scale C
#: of :func:`fit_error_model` is FITTED on ensemble stars and the a-priori
#: C = 1 is only reported beside it.
YOUNG_COEFF = 0.09
SCALE_HEIGHT_M = 8000.0
MAG_PER_REL = 2.5 / math.log(10.0)


def young_scintillation_mag(exptime_s, airmass) -> np.ndarray:
    """A-priori scintillation noise (mag) of one exposure (Young 1967)."""
    d_cm = APERTURE_M * 100.0
    rel = (YOUNG_COEFF * d_cm ** (-2.0 / 3.0)
           * np.power(np.asarray(airmass, float), 1.75)
           * math.exp(-SITE_ALT_M / SCALE_HEIGHT_M)
           / np.sqrt(2.0 * np.asarray(exptime_s, float)))
    return MAG_PER_REL * rel


def total_sigma(sig_phot, scint_unit, c_scint, floor):
    """Per-measurement error: photon (+sky, read) quadrature scintillation
    (scaled) quadrature a flat-field/position floor."""
    return np.sqrt(np.square(sig_phot) + np.square(c_scint * scint_unit)
                   + floor ** 2)


def fit_error_model(resid, sig_phot, scint_unit, dof_factor=None):
    """Maximum-likelihood (C, floor) for residuals ~ N(0, sigma_total^2).

    ``resid`` are ensemble-star residuals about the converged solution;
    ``dof_factor`` (optional, per point) rescales the variance for the
    fraction of a residual absorbed by fitted parameters (n/(n-1) per star).
    Returns (C, floor_mag).  Fitted on ENSEMBLE stars only; the held-out
    check stars test it.
    """
    from scipy.optimize import minimize
    r = np.asarray(resid, float)
    sp_ = np.asarray(sig_phot, float)
    su = np.asarray(scint_unit, float)
    k = np.ones_like(r) if dof_factor is None else np.asarray(dof_factor, float)

    def nll(p):
        c, f = abs(p[0]), abs(p[1])
        v = (sp_ ** 2 + (c * su) ** 2 + f ** 2) / k
        return 0.5 * np.sum(r ** 2 / v + np.log(v))

    best = minimize(nll, x0=[1.0, 0.005], method="Nelder-Mead",
                    options={"xatol": 1e-5, "fatol": 1e-6, "maxiter": 4000})
    return abs(best.x[0]), abs(best.x[1])


def fit_floor(resid, sig_phot, sig_scint, dof_factor=None):
    """ML floor alone, with the scintillation term held fixed."""
    from scipy.optimize import minimize_scalar
    r = np.asarray(resid, float)
    v0 = np.asarray(sig_phot, float) ** 2 + np.asarray(sig_scint, float) ** 2
    k = np.ones_like(r) if dof_factor is None else np.asarray(dof_factor, float)

    def nll(f):
        v = (v0 + f * f) / k
        return 0.5 * np.sum(r ** 2 / v + np.log(v))

    return float(abs(minimize_scalar(nll, bounds=(0.0, 0.3),
                                     method="bounded").x))


def fit_scint_pooled(resid, sig_phot, scint_unit, band_idx, dof_factor=None):
    """ML scintillation scale C SHARED by several bands, each with its own
    floor.  Scintillation is achromatic to first order (Young's law has no
    wavelength term), and pooling is what gives C leverage: within one
    band nearly every frame has the same exposure time, so C and the floor
    are degenerate there, while across 0.5-32 s they are not.
    Returns (C, {band: floor})."""
    from scipy.optimize import minimize
    r = np.asarray(resid, float)
    sp2 = np.asarray(sig_phot, float) ** 2
    su = np.asarray(scint_unit, float)
    b = np.asarray(band_idx, int)
    nb = int(b.max()) + 1
    k = np.ones_like(r) if dof_factor is None else np.asarray(dof_factor, float)

    def nll(p):
        c = abs(p[0])
        fl = np.abs(np.asarray(p[1:]))[b]
        v = (sp2 + (c * su) ** 2 + fl ** 2) / k
        return 0.5 * np.sum(r ** 2 / v + np.log(v))

    best = minimize(nll, x0=[1.5] + [0.02] * nb, method="Nelder-Mead",
                    options={"xatol": 1e-6, "fatol": 1e-7, "maxiter": 20000})
    return abs(best.x[0]), {i: abs(best.x[1 + i]) for i in range(nb)}


def chi2nu(values, sigmas, groups):
    """Constant-source chi2 about each group's weighted mean.

    Returns (chi2_nu, dof).  dof = sum over groups of (n - 1).  Standing
    rule 1: reported as is — no max(chi2_nu, 1), and < 0.5 is a defect."""
    v = np.asarray(values, float)
    s = np.asarray(sigmas, float)
    gr = np.asarray(groups)
    chi2, dof = 0.0, 0
    for g in np.unique(gr):
        m = gr == g
        if m.sum() < 2:
            continue
        w = 1.0 / s[m] ** 2
        mu = np.sum(w * v[m]) / np.sum(w)
        chi2 += float(np.sum((v[m] - mu) ** 2 * w))
        dof += int(m.sum() - 1)
    return (chi2 / dof if dof else float("nan")), dof


def chi2_verdict(c2nu: float) -> str:
    """Standing rule 1: chi2_nu < 0.5 is as much a defect as > 2."""
    if not np.isfinite(c2nu):
        return "untested"
    if c2nu < 0.5:
        return "DEFECT (errors overestimated)"
    if c2nu > 2.0:
        return "DEFECT (errors underestimated)"
    return "pass"


#: The scintillation term is "matched by the check-star rms" (SN-S5 accept)
#: when, in every exposure-time group with >= SCINT_MIN_POINTS bright-check
#: residuals, observed rms / modelled rms lies inside this window.
SCINT_MATCH_WINDOW = (0.8, 1.25)
SCINT_MIN_POINTS = 30


# ===========================================================================
# 4.  ENSEMBLE CALIBRATION
# ===========================================================================

def wls_line(x, y, w):
    """Weighted straight line y = a + b x; returns (a, b, cov)."""
    x, y, w = map(lambda z: np.asarray(z, float), (x, y, w))
    A = np.c_[np.ones_like(x), x]
    W = A * w[:, None]
    cov = np.linalg.inv(A.T @ W)
    a, b = cov @ (W.T @ y)
    return float(a), float(b), cov


def robust_wls_line(x, y, sig, clip=4.0, passes=3):
    """WLS line with iterative clipping at ``clip`` x the weighted rms."""
    x, y, sig = map(lambda z: np.asarray(z, float), (x, y, sig))
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sig) & (sig > 0)
    for _ in range(passes):
        a, b, cov = wls_line(x[keep], y[keep], 1.0 / sig[keep] ** 2)
        res = y - a - b * x
        s = np.sqrt(np.average(res[keep] ** 2, weights=1 / sig[keep] ** 2))
        new = keep & (np.abs(res) <= clip * max(s, 1e-4))
        if new.sum() == keep.sum():
            break
        keep = new
    return a, b, cov, keep


#: Delta-flat (strategy §5.3): a 2-D quadratic in detector position fitted
#: to ensemble-star residuals when the residual surface exceeds ~5 mmag.
#: Pixel coordinates are normalised by this half-size.
DFLAT_NORM_PX = 2048.0


def dflat_design(x, y):
    u = (np.ravel(np.asarray(x, float)) - DFLAT_NORM_PX) / DFLAT_NORM_PX
    v = (np.ravel(np.asarray(y, float)) - DFLAT_NORM_PX) / DFLAT_NORM_PX
    return np.c_[np.ones_like(u), u, v, u * u, u * v, v * v]


def fit_delta_flat(x, y, res, w):
    """Weighted quadratic surface (mag) of residual vs position.  The
    constant term is returned but never applied (it is degenerate with the
    zero points)."""
    A = dflat_design(x, y)
    W = A * np.asarray(w, float)[:, None]
    coef = np.linalg.lstsq(A.T @ W, W.T @ np.asarray(res, float), rcond=None)[0]
    return coef


def eval_delta_flat(coef, x, y):
    """The surface WITHOUT its constant term, in mag."""
    c = np.array(coef, float).copy()
    c[0] = 0.0
    x = np.asarray(x, float)
    return (dflat_design(x.ravel(), np.asarray(y, float).ravel()) @ c
            ).reshape(x.shape)


def solve_band(mag, sig_phot, scint_unit, is_comp, is_check, colour,
               cat_mag, airmass, n_iter: int = 6, fixed_C=None):
    """The full per-band calibration on a (stars x frames) matrix.

    Model:  m_sf = M_s + ZP_f + k2 (c_s - c0)(X_f - 1) + noise
            M_s - m_cat,s = Z + c_term (c_s - c0)        [absolute tie]
    with noise variance sigma_phot^2 + (C S_f)^2 + floor^2.

    ``M_s`` and ``ZP_f`` come from the proven Honeycutt solver
    (``macro_phot.ensemble.solve_ensemble``) run on ensemble stars only;
    ``k2`` is the second-order colour x airmass extinction (strategy §4
    Step 4) regressed from their residuals; (C, floor) the error model
    fitted to the same residuals.  Check stars never vote.

    Returns a dict of the solution and its diagnostics.
    """
    from macro_phot.ensemble import solve_ensemble
    S, F = mag.shape
    dc = colour - COLOUR_PIVOT
    dX = airmass - 1.0
    C, floor, k2 = (1.0 if fixed_C is None else float(fixed_C)), 0.005, 0.0
    comp_rows = is_comp & ~is_check
    for _ in range(n_iter):
        mcorr = mag - k2 * dc[:, None] * dX[None, :]
        sig = total_sigma(sig_phot, scint_unit[None, :], C, floor)
        sol = solve_ensemble(np.where(comp_rows[:, None], mcorr, np.nan), sig)
        res = mcorr - sol.mean_mag[:, None] - sol.zp[None, :]
        use = comp_rows[:, None] & np.isfinite(res) & ~sol.clipped
        # second-order extinction from the ensemble residuals
        # The colour x airmass regressor, DOUBLE-CENTRED over the points in
        # use: the star means absorb its per-star average and the frame zero
        # points its per-frame average, so only the doubly-centred part can
        # survive in the residuals — regressing on the raw product would
        # shrink every increment (the synthetic-recovery test caught it).
        z = np.where(use, dc[:, None] * dX[None, :], np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN rows
            for _ in range(3):
                z = z - np.nan_to_num(np.nanmean(z, axis=1, keepdims=True))
                z = z - np.nan_to_num(np.nanmean(z, axis=0, keepdims=True))
        xx = z[use]
        if xx.size > 10 and np.ptp(xx) > 0:
            w = 1.0 / sig[use] ** 2
            k2 += float(np.sum(w * xx * res[use]) / np.sum(w * xx * xx))
        nobs = np.sum(use, axis=1)
        kf = np.broadcast_to(((nobs - 1) / np.maximum(nobs, 1))[:, None],
                             res.shape)[use]
        if fixed_C is None:
            C, floor = fit_error_model(res[use], sig_phot[use],
                                       np.broadcast_to(scint_unit[None, :],
                                                       res.shape)[use],
                                       dof_factor=np.clip(kf, 0.1, 1))
        else:
            C = float(fixed_C)
            floor = fit_floor(res[use], sig_phot[use],
                              C * np.broadcast_to(scint_unit[None, :],
                                                  res.shape)[use],
                              dof_factor=np.clip(kf, 0.1, 1))
    mcorr = mag - k2 * dc[:, None] * dX[None, :]
    sig = total_sigma(sig_phot, scint_unit[None, :], C, floor)
    # absolute tie on ensemble star means
    M = sol.mean_mag
    ok = comp_rows & np.isfinite(M) & np.isfinite(cat_mag)
    nob = np.sum(np.isfinite(mcorr) & ~sol.clipped, axis=1)
    m_err = 1.0 / np.sqrt(np.nansum(np.where(np.isfinite(mcorr), 1 / sig ** 2,
                                             0), axis=1))
    tie_sig = np.sqrt(m_err ** 2 + 0.01 ** 2)        # + REFCAT2 floor
    Z, cterm, cov, kept = robust_wls_line(dc[ok], (M - cat_mag)[ok],
                                          tie_sig[ok])
    return {"zp": sol.zp, "zp_err": sol.zp_err, "n_star_used": sol.n_star_used,
            "mean_mag": M, "k2": k2, "C": C, "floor": floor, "Z": Z,
            "cterm": cterm, "cterm_err": float(np.sqrt(cov[1, 1])),
            "Z_err": float(np.sqrt(cov[0, 0])), "mcorr": mcorr, "sig": sig,
            "clipped": sol.clipped, "residual": mcorr - M[:, None]
            - sol.zp[None, :], "tie_rows": np.nonzero(ok)[0][kept],
            "nobs": nob}


def check_star_test(mcorr, sig, zp, zp_err, is_check, min_obs=5):
    """Held-out test: calibrated check-star magnitudes about their own
    weighted means, with the full error model plus each frame's ZP error.
    Returns (chi2_nu, dof, per-point residuals, per-point sigmas, rows)."""
    rows = np.nonzero(is_check)[0]
    vals, sigs, grp, cols = [], [], [], []
    for s in rows:
        ok = np.isfinite(mcorr[s]) & np.isfinite(zp)
        if ok.sum() < min_obs:
            continue
        v = mcorr[s][ok] - zp[ok]
        e = np.sqrt(sig[s][ok] ** 2 + zp_err[ok] ** 2)
        # 5-sigma outlier screen on the star's own median (cosmic rays,
        # satellite trails), counted in the returned dof.
        med = np.median(v)
        good = np.abs(v - med) < 5 * e + 0.05
        vals.append(v[good]); sigs.append(e[good])
        grp.append(np.full(good.sum(), s)); cols.append(np.nonzero(ok)[0][good])
    if not vals:
        return float("nan"), 0, np.array([]), np.array([]), np.array([]), \
            np.array([])
    v, e, g, f = map(np.concatenate, (vals, sigs, grp, cols))
    c2, dof = chi2nu(v, e, g)
    mu = {k: np.average(v[g == k], weights=1 / e[g == k] ** 2)
          for k in np.unique(g)}
    res = v - np.array([mu[k] for k in g])
    return c2, dof, res, e, g, f


def natural_to_ps1(m_nat, colour, cterm, k2=0.0, airmass=1.0):
    """PS1-tied magnitude from the natural-system one."""
    dc = np.asarray(colour, float) - COLOUR_PIVOT
    return m_nat - cterm * dc - k2 * dc * (np.asarray(airmass, float) - 1.0)


def colour_from_natural(g_nat, i_nat, cg, ci):
    """Solve the PS1 (g - i) of a source from its natural G and I:
    col = (G - I + (cg - ci) c0) / (1 + cg - ci)."""
    dc = cg - ci
    return (np.asarray(g_nat) - np.asarray(i_nat) + dc * COLOUR_PIVOT) / (1 + dc)


# ===========================================================================
# 5.  NIGHTLY MEANS
# ===========================================================================

def weighted_mean(v, s):
    """(mean, formal error, chi2_nu about the mean, n)."""
    v, s = np.asarray(v, float), np.asarray(s, float)
    w = 1 / s ** 2
    mu = float(np.sum(w * v) / np.sum(w))
    err = float(1 / math.sqrt(np.sum(w)))
    c2 = float(np.sum(w * (v - mu) ** 2) / (len(v) - 1)) if len(v) > 1 \
        else float("nan")
    return mu, err, c2, len(v)


# ===========================================================================
# 6.  VARIABILITY LIMITS  (SN-S8)
# ===========================================================================

#: Smooth-trend model for the nightly light curve: a cubic smoothing spline
#: in phase with its knot spacing floored at 4 d — the night-to-night
#: analogue of the strategy's GP length-scale floor of >= 2 d (a trend that
#: can bend inside 4 d could absorb the very bumps we are testing for).
TREND_KNOT_SPACING_D = 4.0
#: Detection bar: whole-night-bootstrap false-alarm probability < 0.1%
#: (strategy §4 Step 8) — i.e. the statistic must exceed the 99.9th
#: percentile of its bootstrap null.
FAP_LEVEL = 0.001
#: Recovery fraction that DEFINES a limit (strategy: 90%-recovery curve).
RECOVERY_LEVEL = 0.90
#: Bump widths (Gaussian sigma, d) and sinusoid periods (d) on the grid.
#: Night-to-night only: the intra-night periodogram is dropped (ruling).
BUMP_SIGMAS_D = (1.0, 2.0, 4.0)
SINE_PERIODS_D = (2.5, 5.0, 10.0, 20.0)


def spline_trend(t, y, s, knot_spacing=TREND_KNOT_SPACING_D):
    """Weighted least-squares cubic B-spline with interior knots every
    ``knot_spacing`` days; returns the fitted values at t."""
    from scipy.interpolate import make_lsq_spline
    t = np.asarray(t, float)
    order = np.argsort(t)
    ts, ys, ws = t[order], np.asarray(y)[order], 1 / np.asarray(s)[order] ** 2
    k = 3
    n_int = max(int((ts[-1] - ts[0]) // knot_spacing) - 1, 0)
    inner = np.linspace(ts[0], ts[-1], n_int + 2)[1:-1] if n_int else []
    knots = np.r_[[ts[0]] * (k + 1), inner, [ts[-1]] * (k + 1)]
    # Schoenberg-Whitney needs data between knots; drop knots that fail.
    while True:
        try:
            spl = make_lsq_spline(ts, ys, knots, k=k, w=np.sqrt(ws))
            break
        except Exception:
            if len(inner) == 0:
                raise
            inner = inner[::2] if len(inner) > 1 else []
            knots = np.r_[[ts[0]] * (k + 1), inner, [ts[-1]] * (k + 1)]
    out = np.empty_like(t)
    out[order] = spl(ts)
    return out


#: Candidate knot spacings for the trend; the one with the lowest BIC on
#: the band's own nightly curve is used (never below TREND_KNOT_SPACING_D).
TREND_KNOT_CANDIDATES_D = (4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0)


def choose_knot_spacing(t, y, s, candidates=TREND_KNOT_CANDIDATES_D):
    """BIC-selected knot spacing for :func:`spline_trend`.  Returns
    (spacing, {spacing: BIC})."""
    t = np.asarray(t, float)
    out = {}
    for ks in candidates:
        if ks < TREND_KNOT_SPACING_D:
            continue
        try:
            f = spline_trend(t, y, s, ks)
        except Exception:
            continue
        n_int = max(int((t.max() - t.min()) // ks) - 1, 0)
        k = n_int + 4
        chi2 = float(np.sum(((np.asarray(y) - f) / np.asarray(s)) ** 2))
        out[ks] = chi2 + k * math.log(len(t))
    best = min(out, key=out.get)
    return best, out


# ---------------------------------------------------------------------------
# Joint trend + signal search.  The trend and the candidate signal are
# fitted TOGETHER (strategy §4 Step 7/8: "refit ... signal simultaneously"),
# so a bump the trend could absorb earns no significance — and an injected
# bump the trend absorbs is honestly not recovered.  With the nightly
# epochs and errors fixed, every statistic is a fixed linear (bump) or
# quadratic (sine) form in the data, precomputed once.
# ---------------------------------------------------------------------------

def trend_basis(t, knot_spacing):
    """Clamped cubic B-spline design matrix with interior knots every
    ``knot_spacing`` d (the basis :func:`spline_trend` fits)."""
    from scipy.interpolate import BSpline
    t = np.asarray(t, float)
    k = 3
    n_int = max(int((t.max() - t.min()) // knot_spacing) - 1, 0)
    inner = np.linspace(t.min(), t.max(), n_int + 2)[1:-1] if n_int else []
    knots = np.r_[[t.min()] * (k + 1), inner, [t.max()] * (k + 1)]
    return BSpline.design_matrix(t, knots, k).toarray()


class JointSearch:
    """Precomputed matched filters orthogonal to the trend basis."""

    def __init__(self, t, s, knot_spacing, widths=BUMP_SIGMAS_D,
                 periods=SINE_PERIODS_D):
        self.t = np.asarray(t, float)
        w = 1.0 / np.asarray(s, float) ** 2
        self.w = w
        B = trend_basis(self.t, knot_spacing)
        G = np.linalg.pinv((B.T * w) @ B)
        self.H = B @ G @ (B.T * w)          # trend projector (weighted)
        self.bump_vec, self.bump_meta = [], []
        for sig in widths:
            for t0 in self.t:
                p = np.exp(-0.5 * ((self.t - t0) / sig) ** 2)
                pp = p - self.H @ p
                den = float(np.sum(w * pp * pp))
                if den <= 1e-12:
                    continue
                # amplitude estimate a = (pp^T W y)/den, sigma_a = den^-1/2
                self.bump_vec.append(w * pp / den)
                self.bump_meta.append((t0, sig, den ** -0.5))
        self.bump_vec = np.array(self.bump_vec)
        self.bump_sig = np.array([m[2] for m in self.bump_meta])
        self.sine = []
        for P in periods:
            A = np.c_[np.sin(2 * np.pi * self.t / P), np.cos(2 * np.pi * self.t / P)]
            Ap = A - self.H @ A
            M = (Ap.T * w) @ Ap
            if np.linalg.cond(M) > 1e10:
                continue
            Mi = np.linalg.inv(M)
            self.sine.append((P, (Ap.T * w), Mi))

    def trend(self, y):
        return self.H @ np.asarray(y, float)

    def bump(self, y):
        """(max |SNR|, t0, sigma, amplitude) over the bump grid."""
        a = self.bump_vec @ np.asarray(y, float)
        snr = np.abs(a) / self.bump_sig
        k = int(np.argmax(snr))
        t0, sig, _ = self.bump_meta[k]
        return float(snr[k]), t0, sig, float(a[k])

    def bump_at(self, y, t0, sig):
        """Amplitude fitted at a GIVEN centre and width (matched cell)."""
        p = np.exp(-0.5 * ((self.t - t0) / sig) ** 2)
        pp = p - self.H @ p
        return float(np.sum(self.w * pp * y) / np.sum(self.w * pp * pp))

    def sine_stat(self, y, only=None):
        """(max delta-chi2, period, amplitude) over the period grid."""
        best = (0.0, np.nan, np.nan)
        for P, AtW, Mi in self.sine:
            if only is not None and P != only:
                continue
            v = AtW @ np.asarray(y, float)
            beta = Mi @ v
            chi = float(v @ beta)
            if chi > best[0]:
                best = (chi, P, float(np.hypot(*beta)))
        return best
