"""macro_dw.dwcore — the decision rules of the Dwarf-Galaxy Hα paper.

Pure functions only (numpy/scipy; no file, network or database access),
unit-tested in ``pipeline/tests/test_dw.py``.  Every threshold is a module
constant with its provenance, so the methods section, the tests and the
code read the same number.  The plumbing is
``pipeline/scripts/run_dw_paper.py``; the pixels are ``macro_dw.dwio``.

Sections
--------
1. Frames: night labels, the per-frame disposition (DW-P01) and the QC
   gates (DW-P32, strategy §4 Phase 3.2, thresholds declared there).
2. Flats: the per-frame pedestal/sky regression on the twilight master and
   the large-scale residual metric that decides DW-P11 (< 0.3 % of sky).
3. Depth: Román, Trujillo & Montes (2020) surface-brightness limits and
   Sérsic algebra for the synthetic-dwarf injections (DW-P34/P35).
4. Hα: band position from a velocity, the NGC 5238 end-to-end line-flux
   scale and what it implies for the filter's effective width, and the
   detection rule (DW-P36).
5. Time series: the multiband periodogram with floating per-band offsets
   and covariates fitted SIMULTANEOUSLY with the sinusoid (DW-P55), and the
   injection–recovery bookkeeping with the signed matched-cell bias
   (DW-P54, standing rule 3).
"""
from __future__ import annotations

import math
from datetime import date, timedelta
from typing import Optional, Sequence

import numpy as np

# ===========================================================================
# 0.  Physical and instrumental constants
# ===========================================================================
C_KMS = 299792.458
C_AA_S = 2.99792458e18          #: speed of light in Å/s (f_nu <-> f_lambda)
HA_REST_A = 6562.8
#: Nominal AC4040 plate scale (FOCALLEN 3454 mm, 9 µm pixels).  The paper
#: quotes the MEASURED scale from the refined solutions instead.
PIX_ARCSEC_NOMINAL = 0.5375
#: The filter width assumed by the DW-P36-0 depth table (novelty package);
#: kept only for the in-band/out-of-band labels, which the paper states
#: are conditional on it.
W_ASSUMED_A = 65.0

# ===========================================================================
# 1.  Frames: night labels, disposition, QC gates
# ===========================================================================
#: Strategy header: a night is local noon to noon at Winer (UTC-7); label =
#: calendar date of date(JD - 0.7917).
NIGHT_JD_OFFSET = 0.7917
#: A frame whose header pointing is this far from the field's median
#: pointing did not observe the field (the 2023-03-25 NGC 5548 night sits
#: at 8-11 deg).  Dithers are <= 0.03 deg.
MISPOINT_DEG = 0.25
#: Filters that are dispersers, not bandpasses (S2c verdict for slot '6').
SPECTRAL_FILTERS = frozenset({"6", "HaGrism", "OGGrism", "hrg", "lrg"})

#: QC gates, strategy §4 Phase 3.2, declared 2026-08-16 before any stack:
#: FWHM above 12 px; ensemble ZP more than 0.3 mag below the night's median
#: for that filter; transparency loss above 1 mag against the filter's
#: clear-sky zero point (the 90th percentile of all its frames).
QC_FWHM_MAX_PX = 12.0
QC_ZP_BELOW_NIGHT = 0.3
QC_TRANSPARENCY_LOSS = 1.0
QC_ZP_CLEAR_PCTL = 90.0

#: Plate-solution acceptance (fail-closed): matched REFCAT2 stars and the
#: per-star sky rms.  The per-star rms is dominated by centroid noise of
#: faint stars (a header-started 256 s L frame: 301 stars, 0.56" = 1.0 px);
#: the registration error of the fitted solution is rms/sqrt(n) ~ 0.03".
#: 1.5 px admits that centroid noise and rejects a wrong quad match, whose
#: residuals are tens of pixels.
WCS_MIN_MATCH = 10
WCS_MAX_RMS_PX = 1.5


def night_label(jd: float) -> str:
    """Winer local night label for a UTC Julian date."""
    d = date(2000, 1, 1) + timedelta(days=int(math.floor(jd - NIGHT_JD_OFFSET
                                                        - 2451544.5)))
    return d.isoformat()


def disposition(is_canonical: int, filt: str, offset_deg: Optional[float],
                solved: Optional[bool], qc: Sequence[str]) -> str:
    """One label per staged science frame, in a fixed precedence order:
    duplicate > mispointed > spectrum > unsolved > qc_rejected > science
    (a mispointed spectrum — the 2023-03-25 NGC 5548 night — is reported as
    mispointed: what it failed first is pointing).

    ``solved`` is None for frames never queued for a solve (spectra)."""
    if not is_canonical:
        return "duplicate"
    if offset_deg is not None and offset_deg > MISPOINT_DEG:
        return "mispointed"
    if filt in SPECTRAL_FILTERS:
        return "spectrum"
    if solved is False:
        return "unsolved"
    if qc:
        return "qc_rejected"
    return "science"


def qc_reasons(fwhm_px: Optional[float], zp: Optional[float],
               zp_night_median: Optional[float],
               zp_clear: Optional[float],
               sky_raw_adu: Optional[float] = None,
               lin_cap_adu: Optional[float] = None) -> list[str]:
    """The Phase 3.2 gates that a frame fails (empty list = passes), plus
    standing rule 4 judged at the target: a frame whose raw SKY level
    already exceeds the measured linearity cap (detector_params) puts every
    target pixel in the non-linear regime."""
    out = []
    if sky_raw_adu is not None and lin_cap_adu is not None \
            and sky_raw_adu > lin_cap_adu:
        out.append("nonlinear_sky")
    if fwhm_px is None or not np.isfinite(fwhm_px):
        out.append("no_fwhm")
    elif fwhm_px > QC_FWHM_MAX_PX:
        out.append("fwhm")
    if zp is None or not np.isfinite(zp):
        out.append("no_zp")
        return out
    if zp_night_median is not None and zp < zp_night_median - QC_ZP_BELOW_NIGHT:
        out.append("zp_below_night")
    if zp_clear is not None and zp < zp_clear - QC_TRANSPARENCY_LOSS:
        out.append("transparency")
    return out


def wcs_ok(n_match: Optional[int], rms_arcsec: Optional[float],
           scale_arcsec: float) -> bool:
    """Fail-closed acceptance of a refined plate solution."""
    if n_match is None or rms_arcsec is None or not np.isfinite(rms_arcsec):
        return False
    return n_match >= WCS_MIN_MATCH and rms_arcsec <= WCS_MAX_RMS_PX * scale_arcsec


# ===========================================================================
# 2.  Flats (DW-P11) and gradients (DW-P12)
# ===========================================================================
#: Block size (px) of the per-frame sky maps; 64 px = 35".
FLAT_BLOCK = 64
#: DW-P11 acceptance: held-out residual large-scale structure < 0.3 % of sky.
FLAT_RESID_MAX = 0.003
#: Fringe/periodic-structure threshold of strategy §4 Phase 1.2.
FRINGE_MAX = 0.002


def robust_sigma(a: np.ndarray) -> float:
    """1.4826 x MAD of the finite values (NaN for fewer than 3)."""
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size < 3:
        return float("nan")
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def fit_sky_on_template(block: np.ndarray, tmpl: np.ndarray,
                        clip: float = 3.0, iters: int = 4
                        ) -> tuple[float, float]:
    """Fit ``block = s0 * tmpl + c`` over the finite blocks (clipped LSQ).

    ``tmpl`` is the flat-field shape at block resolution (vignetting
    included), so ``s0`` is the frame's sky in flat-fielded ADU and ``c``
    the additive offset that a pedestal mismatch between the frame and its
    dark master leaves.  Removing ``c`` before normalising is what keeps a
    constant offset from diluting the superflat's contrast."""
    b, t = np.asarray(block, float).ravel(), np.asarray(tmpl, float).ravel()
    m = np.isfinite(b) & np.isfinite(t)
    keep = m.copy()
    s0, c = float("nan"), float("nan")
    for _ in range(iters):
        if keep.sum() < 8:
            break
        A = np.c_[t[keep], np.ones(keep.sum())]
        (s0, c), *_ = np.linalg.lstsq(A, b[keep], rcond=None)
        r = b - (s0 * t + c)
        s = robust_sigma(r[keep])
        keep = m & (np.abs(r) < clip * max(s, 1e-12))
    return float(s0), float(c)


def plane_design(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.c_[np.ones(len(x)), x, y]


def remove_plane(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares plane removed from the finite pixels of a 2-D map.
    Returns (residual map, coefficients [c0, cx, cy] per unit pixel)."""
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    m = np.isfinite(img)
    if m.sum() < 4:
        return img.copy(), np.full(3, np.nan)
    A = plane_design(xx[m] - nx / 2, yy[m] - ny / 2)
    coef, *_ = np.linalg.lstsq(A, img[m], rcond=None)
    out = img - (coef[0] + coef[1] * (xx - nx / 2) + coef[2] * (yy - ny / 2))
    return out, coef


def residual_structure(resid: np.ndarray) -> dict:
    """The DW-P11 metric on a held-out normalised residual map (1 = flat).

    Returns the robust rms of the map about 1 (``rms``), the same after a
    plane is removed (``rms_noplane``: the part the per-frame sky plane of
    Phase 1.2 cannot remove), and the peak-to-valley of the 5th-95th
    percentile.  All as fractions of sky."""
    r = np.asarray(resid, float) - 1.0
    noplane, _ = remove_plane(r)
    f = r[np.isfinite(r)]
    return {"rms": robust_sigma(f),
            "rms_noplane": robust_sigma(noplane[np.isfinite(noplane)]),
            "p95_p05": float(np.percentile(f, 95) - np.percentile(f, 5))
            if f.size else float("nan")}


# ===========================================================================
# 3.  Depth (DW-P34) and structure (DW-P35)
# ===========================================================================
#: Román, Trujillo & Montes (2020, A&A 644, A42): 3 sigma in 10" x 10" boxes.
ROMAN_BOX_ARCSEC = 10.0
ROMAN_NSIG = 3.0
#: Synthetic-dwarf grid (exponential discs, n = 1): central surface
#: brightness and effective radius.  Detection = the disc's flux inside
#: r_e exceeds DETECT_INJ_SIGMA times the empirical noise of r_e apertures
#: (random placements in the same stack, so correlated noise is included).
INJ_MU0 = np.arange(22.0, 27.51, 0.5)
INJ_RE_ARCSEC = np.array([3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0])
INJ_PER_CELL = 20
DETECT_INJ_SIGMA = 5.0


def roman_mu_limit(sigma_box_mean: float, zp: float, pix_arcsec: float,
                   nsig: float = ROMAN_NSIG) -> float:
    """Surface-brightness limit (mag/arcsec^2) from the scatter of mean
    pixel values in boxes, ``sigma_box_mean`` (flux units per pixel, on a
    scale where m = zp - 2.5 log10(flux))."""
    if not (sigma_box_mean > 0):
        return float("nan")
    return float(zp - 2.5 * np.log10(nsig * sigma_box_mean / pix_arcsec ** 2))


def sersic_bn(n: float) -> float:
    """Ciotti & Bertin (1999) asymptotic b_n (accurate to 1e-4 for n>0.36)."""
    return (2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n ** 2)
            + 131 / (1148175 * n ** 3))


def sersic_mu0(mu_e: float, n: float) -> float:
    """Central surface brightness from the one at r_e."""
    return mu_e - 2.5 * sersic_bn(n) / math.log(10)


def sersic_total_mag(mu_e: float, re_arcsec: float, n: float,
                     q: float = 1.0) -> float:
    """Total magnitude of a Sérsic profile (Graham & Driver 2005, eq. 12)."""
    from scipy.special import gamma
    b = sersic_bn(n)
    lum = (2 * math.pi * n * math.exp(b) / b ** (2 * n) * gamma(2 * n)
           * re_arcsec ** 2 * q)
    return mu_e - 2.5 * math.log10(lum)


def exp_disc_image(n_pix: int, x0: float, y0: float, mu0: float,
                   re_arcsec: float, zp: float, pix_arcsec: float,
                   q: float = 1.0, pa_deg: float = 0.0) -> np.ndarray:
    """An exponential disc (Sérsic n = 1) in flux per pixel on a
    ``zp``-scaled image.  r_e = 1.678 h."""
    yy, xx = np.mgrid[0:n_pix, 0:n_pix]
    th = math.radians(pa_deg)
    dx, dy = xx - x0, yy - y0
    u = dx * math.cos(th) + dy * math.sin(th)
    v = -dx * math.sin(th) + dy * math.cos(th)
    r = np.hypot(u, v / q) * pix_arcsec
    h = re_arcsec / sersic_bn(1.0)
    i0 = 10 ** (-0.4 * (mu0 - zp)) * pix_arcsec ** 2
    return i0 * np.exp(-r / h)


def recovery_contour(rec_frac: np.ndarray, mu0: np.ndarray,
                     level: float) -> np.ndarray:
    """For each r_e column, the faintest mu0 at which the recovered
    fraction is still >= ``level`` (linear interpolation; NaN if never)."""
    out = np.full(rec_frac.shape[1], np.nan)
    for j in range(rec_frac.shape[1]):
        f = rec_frac[:, j]
        ok = np.nonzero(f >= level)[0]
        if ok.size == 0:
            continue
        k = ok.max()
        if k == len(mu0) - 1:
            out[j] = mu0[k]
            continue
        f1, f2 = f[k], f[k + 1]
        out[j] = mu0[k] + (mu0[k + 1] - mu0[k]) * (f1 - level) / max(f1 - f2, 1e-9)
    return out


# ===========================================================================
# 4.  Hα (DW-P36) and the NGC 5238 line-flux scale
# ===========================================================================
#: Aperture of the depth table (the physicist's 10"), kept like-for-like.
HA_APER_ARCSEC = 10.0
#: Pre-declared detection rule (2026-10-05, before any stack existed):
#: detection if net/sigma >= 5 with sigma from random apertures on the
#: SAME continuum-subtracted stack (correlated noise and subtraction
#: residuals included); 3 <= net/sigma < 5 is "marginal" and is reported
#: as a non-detection with its limit; the limit is 3 sigma above max(net, 0).
#: Thirteen fields are tested, so a 5-sigma threshold keeps the expected
#: false positives below 1e-5 for Gaussian noise.
HA_DETECT_SIGMA = 5.0
HA_MARGINAL_SIGMA = 3.0


def ha_offset_A(v_kms: Optional[float]) -> Optional[float]:
    """Observed-minus-rest wavelength of Hα (Å) at heliocentric velocity v."""
    if v_kms is None or not np.isfinite(v_kms):
        return None
    return HA_REST_A * v_kms / C_KMS


def flam_per_rate(zp_ab: float, lam_a: float = HA_REST_A) -> float:
    """f_lambda (erg/s/cm^2/Å) of a source giving 1 ADU/s, for an AB zero
    point ``zp_ab`` (m_AB = zp_ab - 2.5 log10 rate)."""
    fnu = 10 ** (-0.4 * (zp_ab + 48.6))
    return fnu * C_AA_S / lam_a ** 2


def line_scale(f_pub: float, rate: float) -> float:
    """Line flux (erg/s/cm^2) per ADU/s from a published total flux and the
    measured continuum-subtracted count rate of the same galaxy."""
    return f_pub / rate


def effective_width(scale_line: float, zp_ab: float,
                    lam_a: float = HA_REST_A) -> float:
    """W_eff = (∫T dλ) / T(λ_line), in Å.

    A continuum of flux density f_λ gives rate f_λ ∫T dλ / K; a line of flux
    F at λ_line gives rate F T(λ_line) / K.  The ratio of the two scales is
    therefore the filter's equivalent width as seen by THAT line — not the
    FWHM, not the centre."""
    return scale_line / flam_per_rate(zp_ab, lam_a)


def ha_verdict(snr: Optional[float], band: str) -> str:
    """Detection label for one object (``band`` = in / out / unknown)."""
    if snr is None or not np.isfinite(snr):
        return "no_measurement"
    if snr >= HA_DETECT_SIGMA:
        return "detected"
    if snr >= HA_MARGINAL_SIGMA:
        return "marginal"
    return "not_detected"


def upper_limit(net: float, sigma: float, nsig: float = 3.0) -> float:
    """Upper limit on a flux: nsig x sigma above max(net, 0)."""
    return max(net, 0.0) + nsig * sigma


# ===========================================================================
# 5.  Time series: multiband periodogram with covariates (DW-P55) and
#     injection–recovery (DW-P54)
# ===========================================================================
#: Period grid of strategy §4 Phase 5.4.
INJ_PERIODS_D = np.geomspace(0.05, 30.0, 16)
INJ_AMPS_MMAG = np.array([2, 5, 10, 20, 50, 100, 200], float)
INJ_PHASES = 20
#: Recovery: the highest peak of the periodogram lies within one
#: resolution element (1/T) of the injected frequency and exceeds the
#: star's own 1 % false-alarm power from shuffled residuals.
FAP_LEVEL = 0.01
N_SHUFFLE = 200
#: At most two Sys-Rem-like components (strategy §4 Phase 5.5).
MAX_SYSREM = 2


def nuisance_basis(band_idx: np.ndarray, covariates: np.ndarray) -> np.ndarray:
    """Design matrix of the nuisance model: one floating offset per band
    plus the frame-level covariates (columns standardised)."""
    bands = np.unique(band_idx)
    X = [(band_idx == b).astype(float) for b in bands]
    cov = np.atleast_2d(np.asarray(covariates, float))
    if cov.size:
        if cov.shape[0] != len(band_idx):
            cov = cov.T
        for j in range(cov.shape[1]):
            c = cov[:, j]
            s = np.std(c)
            if s > 0:
                X.append((c - np.mean(c)) / s)
    return np.column_stack(X)


class Periodogram:
    """Generalised least-squares periodogram with nuisance regressors.

    For each trial frequency f the model is  y = X b + a sin(2πft) +
    c cos(2πft) + noise,  with X the per-band offsets and covariates; the
    power is the fractional chi^2 reduction the sinusoid buys OVER the
    nuisance-only fit.  Because X is the same at every f, the trig columns
    are projected off X once per frequency and the per-trial cost is two
    matrix products — which is what makes thousands of injections cheap.
    Trend and signal are thus fitted simultaneously (never detrend first).
    """

    def __init__(self, t: np.ndarray, w: np.ndarray, X: np.ndarray,
                 freqs: np.ndarray):
        self.t = np.asarray(t, float)
        self.sw = np.sqrt(np.asarray(w, float))
        self.f = np.asarray(freqs, float)
        Xw = X * self.sw[:, None]
        self.Q, _ = np.linalg.qr(Xw)
        ph = 2 * np.pi * np.outer(self.f, self.t)            # F x N
        S = np.sin(ph) * self.sw
        C = np.cos(ph) * self.sw
        S -= (S @ self.Q) @ self.Q.T
        C -= (C @ self.Q) @ self.Q.T
        self.S, self.C = S, C
        ss = np.einsum("ij,ij->i", S, S)
        cc = np.einsum("ij,ij->i", C, C)
        sc = np.einsum("ij,ij->i", S, C)
        det = ss * cc - sc ** 2
        det[det <= 0] = np.inf
        self.ginv = np.stack([cc / det, -sc / det, ss / det], 1)

    def residualise(self, Y: np.ndarray) -> np.ndarray:
        """Weighted data with the nuisance fit removed (N x T)."""
        Yw = np.atleast_2d(Y.T).T * self.sw[:, None]
        return Yw - self.Q @ (self.Q.T @ Yw)

    def power(self, Y: np.ndarray) -> np.ndarray:
        """Fractional chi^2 reduction, F x T, for data columns Y (N x T)."""
        R = self.residualise(Y)
        a = self.S @ R
        b = self.C @ R
        g = self.ginv
        dchi = (g[:, 0:1] * a * a + 2 * g[:, 1:2] * a * b + g[:, 2:3] * b * b)
        chi0 = np.einsum("ij,ij->j", R, R)
        return dchi / np.maximum(chi0, 1e-30)

    def amplitude(self, Y: np.ndarray, k: np.ndarray) -> np.ndarray:
        """Fitted semi-amplitude at frequency index k[j] for column j."""
        R = self.residualise(Y)
        cols = np.arange(R.shape[1])
        a = np.einsum("ij,ji->i", self.S[k], R)
        b = np.einsum("ij,ji->i", self.C[k], R)
        g = self.ginv[k]
        sa = g[:, 0] * a + g[:, 1] * b
        ca = g[:, 1] * a + g[:, 2] * b
        del cols
        return np.hypot(sa, ca)


def freq_grid(t: np.ndarray, fmax: float, oversample: int = 5) -> np.ndarray:
    """Uniform frequency grid from 1/T to ``fmax`` at 1/(oversample T)."""
    T = float(np.ptp(t))
    df = 1.0 / (oversample * T)
    return np.arange(1.0 / T, fmax + df, df)


def shuffle_threshold(pg: Periodogram, resid: np.ndarray, rng,
                      n: int = N_SHUFFLE, level: float = FAP_LEVEL) -> float:
    """Max-power quantile of shuffled residuals (time stamps kept, values
    permuted): the empirical 1 % false-alarm power for this star."""
    Y = np.column_stack([rng.permutation(resid) for _ in range(n)])
    pmax = pg.power(Y).max(axis=0)
    return float(np.quantile(pmax, 1 - level))


def recovered(pg: Periodogram, P: np.ndarray, f_inj: np.ndarray,
              thresh: float, T: float) -> np.ndarray:
    """Boolean per column: peak frequency within 1/T of injection and peak
    power above the false-alarm threshold."""
    k = np.argmax(P, axis=0)
    fpk = pg.f[k]
    pk = P[k, np.arange(P.shape[1])]
    return (np.abs(fpk - f_inj) <= 1.0 / T) & (pk > thresh)


def signed_cell_bias(a_rec: np.ndarray, a_inj: np.ndarray,
                     ok: np.ndarray) -> float:
    """Signed matched-cell bias: mean of (A_rec - A_inj)/A_inj over the
    RECOVERED trials of one cell (standing rule 3: signed, never |.|)."""
    m = ok & np.isfinite(a_rec)
    if not m.any():
        return float("nan")
    return float(np.mean((a_rec[m] - a_inj[m]) / a_inj[m]))


def amp90(amps: np.ndarray, frac: np.ndarray, level: float = 0.9) -> float:
    """Smallest injected amplitude at which recovery reaches ``level``
    (log-linear interpolation; NaN if never)."""
    ok = np.nonzero(frac >= level)[0]
    if ok.size == 0:
        return float("nan")
    k = ok.min()
    if k == 0:
        return float(amps[0])
    f1, f2 = frac[k - 1], frac[k]
    la, lb = np.log(amps[k - 1]), np.log(amps[k])
    return float(np.exp(la + (lb - la) * (level - f1) / max(f2 - f1, 1e-9)))


def sysrem(resid: np.ndarray, err: np.ndarray, n_comp: int = MAX_SYSREM,
           iters: int = 20) -> np.ndarray:
    """Tamuz, Mazeh & Zucker (2005) Sys-Rem on a stars x frames residual
    matrix (NaN = missing).  Returns the frame vectors (n_comp x frames),
    which then enter each star's fit as covariates — they are never
    subtracted on their own."""
    R = np.array(resid, float)
    W = 1.0 / np.square(np.asarray(err, float))
    W[~np.isfinite(R)] = 0.0
    R[~np.isfinite(R)] = 0.0
    comps = []
    for _ in range(n_comp):
        a = np.ones(R.shape[1])
        for _ in range(iters):
            c = (R * W) @ a / np.maximum((W * a ** 2).sum(1), 1e-30)
            a = (R * W).T @ c / np.maximum((W.T * c ** 2).sum(1), 1e-30)
        comps.append(a.copy())
        R = R - np.outer(c, a)
    return np.array(comps)
