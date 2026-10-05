"""Horne-style optimal extraction for slitless grism frames.

WHAT CHANGED ON 2026-10-03, AND WHY (committee findings DE.F1, OA.E7, OA.E8)
-----------------------------------------------------------------------------
The first version of this module carried three hardware constants and one
background model.  The committee showed all four were wrong:

* **Variance.**  It built the Horne variance from header ``EGAIN`` 0.2467
  e-/ADU and 3.5 e- read noise.  Mode0 frames are 2x2 *average*-binned, so
  the effective gain is ~1.0 e-/ADU and the read noise ~2.1 ADU: the read
  term was 201 ADU^2 where ~4.4 is measured, the shot term four times too
  large, and it charged shot noise on the 303 ADU electronic pedestal.
  Every ``halpha_snr`` and every error bar was wrong by about a factor of
  two.  The variance now comes from :class:`macro_grism.config.Detector`
  (one record per camera configuration, fed by the detector package's
  photon-transfer table) and nothing in this file holds a gain.
* **Saturation.**  It masked every pixel above 16,300 ADU as "the Mode0
  rail".  There is no rail: the pile-up near 16.6 kADU is one saturated
  NATIVE hot pixel averaged with three ordinary ones.  Real trace pixels
  between 16.3k and ~56k were being thrown away — on T CrB that is the
  Halpha peak itself.  Saturation is now the detector's native clip
  de-rated by the measured native-vs-binned peaking factor, and isolated
  hot pixels are rejected where they belong: by the optimal-extraction
  outlier test against the spatial profile, and by an optional hot-pixel
  mask.
* **Background.**  It drew a straight line through two flanking bands.  On
  the low-resolution grism the dispersed sky is the image of the field
  stop — a lozenge with a curved top and sharp ends — and a straight line
  through curved sky over-subtracts (OA.E7).  The sky model is now a
  robust per-column POLYNOMIAL through wide flanks (the straight line is
  kept as ``sky="flanking"``, the legacy arm, so the difference between
  the two is a measured systematic and not an argument).

BACKGROUND POLICY
-----------------
* ``sky="flanking"`` — two bands, a straight line (the v1 arm; exact for a
  linear gradient, biased by curvature).
* ``sky="poly"``     — a degree-``SKY_POLY_DEG`` polynomial per column,
  fitted with iterative clipping to every flank row between
  ``SKY_INNER`` and ``SKY_OUTER`` px either side of the trace.  Wide
  flanks see the curvature; clipping removes the instrumental ghost
  traces that run parallel to the main trace inside the flanks.
* a 2-D ``template`` (the sky-lozenge model of ``macro_grism.background``)
  may be subtracted first; either arm then removes what the template left.

How well each arm works is not asserted here: ``run_g_background.py``
extracts an EMPTY aperture offset from the trace (a null test) with every
arm and reports the residual.

THE ESTIMATOR
-------------
Horne (1986), vectorised over columns.  With spatial profile P (unit sum
per column), data D, sky S and variance V:

    f = sum(M P (D - S) / V) / sum(M P^2 / V),   var(f) = 1 / sum(M P^2 / V)

M masks saturated and outlier pixels.  Two details matter for bias:

* the variance is evaluated from the MODEL (f P + S), not from the noisy
  data, after the first pass — data-based weights correlate with the
  noise and bias f low at low signal (Horne 1986, section II.e);
* the profile is allowed to vary slowly along the dispersion
  (``PROFILE_CHUNKS`` chunks, linearly interpolated): a grism's focus is
  chromatic, and one global profile under-weights the wings where the
  trace is broad.  A wrong profile biases the *level* smoothly with
  wavelength; an equivalent width, being a local ratio, is immune to
  first order — and the plain aperture sum ``box`` is carried beside
  ``flux`` so the difference is measured per frame.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .config import Detector
from .trace import CENTROID_HALFWIN

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report quotes these).
# NOTE: no detector constant lives here.  See macro_grism.config.
# --------------------------------------------------------------------------

#: Extraction aperture half-width (px).  The trace FWHM measures ~5-9 px
#: on the validation frames; +/-12 px holds > 99% of a 9 px FWHM Gaussian.
APERTURE_HALFWIN = CENTROID_HALFWIN

#: Legacy flanking bands: inner and outer edge offsets from the trace
#: center (px), used on BOTH sides.
BAND_INNER = 30
BAND_OUTER = 60

#: Polynomial sky: flank rows from SKY_INNER to SKY_OUTER px either side
#: of the trace, degree SKY_POLY_DEG, SKY_CLIP_SIGMA iterative clipping.
#: The inner edge clears the trace wings (FWHM up to ~12 px on defocused
#: nights -> 20 px is > 3.9 sigma); the outer edge spans enough of the
#: field-stop lozenge (~1500 px tall) to see its curvature while staying
#: far inside its sharp ends.
SKY_INNER = 20
SKY_OUTER = 120
SKY_POLY_DEG = 2
SKY_CLIP_SIGMA = 3.0
SKY_CLIP_ITERS = 4

#: The fitted sky is smoothed ALONG the dispersion with a running median
#: of this many columns.  A slitless sky has no sharp spectral features —
#: every sky line is smeared over the full width of the field stop — so
#: the true sky is smooth in x on this scale, while the per-column fit
#: carries independent noise (~1.5 ADU per pixel, fully correlated down
#: the 25 aperture rows, i.e. ~40 ADU per column on a 240 s frame).
#: Smoothing removes that noise from the extracted spectrum instead of
#: adding it to every column.  31 columns is ~1/7 of the sharpest sky
#: structure measured (the ~200-column ends of the lrg lozenge).
SKY_SMOOTH_X = 31

#: Half-height of the rectified cutout around the trace (px).  Must hold
#: the widest sky flank.
RECT_HALF = SKY_OUTER

#: Spatial profile: built in this many chunks along the dispersion and
#: linearly interpolated between chunk centers.
PROFILE_CHUNKS = 16

#: Outlier rejection in the Horne sum: a pixel whose residual from
#: f*P exceeds this many sigma is masked (hot pixels, cosmic rays).  The
#: variance used in the test carries a PROFILE_TOL fractional allowance
#: for profile mismatch, so a bright, slightly mis-modelled core is not
#: mistaken for a cosmic ray.
HORNE_CLIP_SIGMA = 6.0
HORNE_PROFILE_TOL = 0.10
HORNE_ITERS = 3


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
def rectify(data: np.ndarray, yc: np.ndarray,
            half: int = RECT_HALF) -> np.ndarray:
    """Cut a strip of ``2*half+1`` rows that FOLLOWS the trace.

    Row ``half`` of the result is the pixel row nearest the trace center
    in every column (integer shifts — no resampling, so no correlated
    noise and no interpolation of hot pixels into their neighbours).
    Rows that fall off the frame are NaN.  Shape: (2*half+1, nx).
    """
    ny, nx = data.shape
    y0 = np.round(yc).astype(int)
    rows = y0[None, :] + np.arange(-half, half + 1)[:, None]
    inside = (rows >= 0) & (rows < ny)
    cut = data[np.clip(rows, 0, ny - 1),
               np.arange(nx)[None, :]].astype(np.float64)
    cut[~inside] = np.nan
    return cut


# --------------------------------------------------------------------------
# Sky models (all return the sky at EVERY row of the rectified cutout)
# --------------------------------------------------------------------------
def flanking_background(column: np.ndarray, yc: float,
                        halfwin: int = APERTURE_HALFWIN,
                        inner: int = BAND_INNER,
                        outer: int = BAND_OUTER):
    """Per-column background from the two flanking bands (legacy arm,
    scalar form kept for the unit tests and for one-off use).

    Returns ``(background_at_window_rows, ok)`` where the background is a
    straight line through (band centers, band medians) evaluated at the
    aperture rows — a LINEAR local model, so a smooth vertical gradient in
    halo/sky/dark is removed exactly (unit-tested).  ``ok`` is False when
    either band leaves the frame; the caller records the column as
    unextractable rather than inventing a background.
    """
    ny = len(column)
    lo, hi = int(round(yc)) - halfwin, int(round(yc)) + halfwin + 1
    b_lo0, b_lo1 = int(round(yc)) - outer, int(round(yc)) - inner
    b_hi0, b_hi1 = int(round(yc)) + inner, int(round(yc)) + outer
    if b_lo0 < 0 or b_hi1 + 1 > ny:
        return None, False
    below = np.median(column[b_lo0:b_lo1])
    above = np.median(column[b_hi0:b_hi1 + 1])
    y_below = (b_lo0 + b_lo1 - 1) / 2.0          # band center rows
    y_above = (b_hi0 + b_hi1) / 2.0
    slope = (above - below) / (y_above - y_below)
    rows = np.arange(lo, hi)
    return below + slope * (rows - y_below), True


def sky_flanking(cut: np.ndarray, inner: int = BAND_INNER,
                 outer: int = BAND_OUTER) -> np.ndarray:
    """Vectorised legacy arm: a straight line per column through the
    medians of the two flanking bands of a rectified cutout.  Columns
    with an incomplete band are NaN (unextractable, not invented)."""
    half = cut.shape[0] // 2
    dy = np.arange(-half, half + 1, dtype=float)
    lo = cut[half - outer:half - inner]
    hi = cut[half + inner:half + outer + 1]
    with np.errstate(invalid="ignore"), _quiet():
        below = np.nanmedian(lo, axis=0)
        above = np.nanmedian(hi, axis=0)
    bad = np.isnan(lo).any(axis=0) | np.isnan(hi).any(axis=0)
    y_below = -(outer + inner + 1) / 2.0
    y_above = (outer + inner) / 2.0
    slope = (above - below) / (y_above - y_below)
    sky = below[None, :] + slope[None, :] * (dy[:, None] - y_below)
    sky[:, bad] = np.nan
    return sky


def sky_polynomial(cut: np.ndarray, inner: int = SKY_INNER,
                   outer: int = SKY_OUTER, deg: int = SKY_POLY_DEG,
                   clip: float = SKY_CLIP_SIGMA,
                   iters: int = SKY_CLIP_ITERS,
                   exclude: Optional[np.ndarray] = None,
                   smooth_x: int = SKY_SMOOTH_X) -> np.ndarray:
    """Robust polynomial sky: per column, a degree-``deg`` polynomial in
    the cross-dispersion offset, fitted to the flank rows with iterative
    sigma clipping, evaluated at every row, then median-smoothed along
    the dispersion over ``smooth_x`` columns (0/1 disables; see
    ``SKY_SMOOTH_X`` for why smoothing is physically justified).

    Solved for ALL columns at once as a masked linear least-squares
    problem (normal equations per column; the design matrix is shared).
    ``exclude`` is an optional boolean row mask (True = never use this
    row), for known ghost traces.  Columns with fewer than ``4*(deg+1)``
    usable flank pixels are NaN.
    """
    half = cut.shape[0] // 2
    dy = np.arange(-half, half + 1, dtype=float)
    flank = (np.abs(dy) >= inner) & (np.abs(dy) <= outer)
    if exclude is not None:
        flank &= ~exclude
    # Scaled abscissa keeps the normal equations well conditioned.
    t = dy / float(outer)
    basis = np.vstack([t ** k for k in range(deg + 1)])      # (deg+1, R)
    use = flank[:, None] & np.isfinite(cut)                  # (R, nx)
    d = np.where(use, cut, 0.0)
    coef = np.full((deg + 1, cut.shape[1]), np.nan)
    for _ in range(iters + 1):
        w = use.astype(float)
        # Normal matrix A[c] = sum_r w B_i B_j ; rhs b[c] = sum_r w B_i d
        a = np.einsum("ir,jr,rc->cij", basis, basis, w)
        b = np.einsum("ir,rc->ci", basis, w * d)
        enough = use.sum(axis=0) >= 4 * (deg + 1)
        coef[:] = np.nan
        if enough.any():
            coef[:, enough] = np.linalg.solve(
                a[enough], b[enough][..., None])[..., 0].T
        model = np.einsum("ic,ir->rc", coef, basis)
        resid = np.where(use, cut - model, np.nan)
        with np.errstate(invalid="ignore"), _quiet():
            med = np.nanmedian(resid, axis=0)
            sig = 1.4826 * np.nanmedian(np.abs(resid - med), axis=0)
        new_use = use & (np.abs(resid - med) <= clip * np.maximum(sig, 1e-6))
        if (new_use == use).all():
            break
        use = new_use
    if smooth_x and smooth_x > 1:
        from scipy.ndimage import median_filter
        ok = np.isfinite(model)
        filled = np.where(ok, model, 0.0)
        # NaN-aware running median is expensive; the model is NaN only in
        # whole columns (frame edges), so filter the finite block and
        # restore the NaNs afterwards.
        good_cols = ok.all(axis=0)
        if good_cols.sum() > smooth_x:
            sm = median_filter(filled[:, good_cols],
                               size=(1, smooth_x), mode="nearest")
            model = model.copy()
            model[:, good_cols] = sm
    return model


class _quiet:
    """Silence numpy's 'All-NaN slice' RuntimeWarning inside a block: an
    all-NaN column is an expected, handled state here (frame edge), not a
    surprise worth a line of log per column."""

    def __enter__(self):
        import warnings
        self._cm = warnings.catch_warnings()
        self._cm.__enter__()
        warnings.simplefilter("ignore", category=RuntimeWarning)

    def __exit__(self, *exc):
        return self._cm.__exit__(*exc)


# --------------------------------------------------------------------------
# Spatial profile
# --------------------------------------------------------------------------
def build_profile(cutouts: np.ndarray) -> np.ndarray:
    """Normalized spatial profile from background-subtracted, per-column
    flux-normalized cutouts (shape: n_columns x window).

    Median over columns (cosmic rays lose), floored at zero (a profile is
    a probability density; negative wings are noise), normalized to unit
    sum.  Falls back to a flat profile if everything medians to zero —
    the extraction then degrades to a box sum instead of dividing by 0.
    """
    with _quiet():
        prof = np.nanmedian(cutouts, axis=0)
    prof = np.clip(np.nan_to_num(prof), 0.0, None)
    total = prof.sum()
    if total <= 0:
        return np.full(cutouts.shape[1], 1.0 / cutouts.shape[1])
    return prof / total


def profile_map(net: np.ndarray, n_chunks: int = PROFILE_CHUNKS):
    """Slowly varying spatial profile P(row, column) for an aperture
    window ``net`` (rows x nx, sky-subtracted).

    Each chunk's profile is the median of its columns' unit-sum cutouts
    (only columns with positive net flux vote — a profile of noise is a
    random number); profiles are then linearly interpolated between chunk
    centers and re-normalized per column.  Returns (P, chunk_fwhm_px,
    chunk_centers): the per-chunk FWHM of the profile is the trace's
    cross-dispersion width, which for a point source is also the line
    spread function along the dispersion (see ``macro_grism.lsf``).
    """
    n_rows, nx = net.shape
    edges = np.linspace(0, nx, n_chunks + 1).astype(int)
    centers, profs, fwhms = [], [], []
    for i in range(n_chunks):
        seg = net[:, edges[i]:edges[i + 1]]
        tot = np.nansum(seg, axis=0)
        noise = 1.4826 * np.nanmedian(np.abs(seg - np.nanmedian(seg))) \
            if np.isfinite(seg).any() else np.nan
        good = np.isfinite(tot) & (tot > 0)
        if np.isfinite(noise):
            good &= tot > 5.0 * noise * np.sqrt(n_rows)
        if good.sum() < 8:
            continue
        p = build_profile((seg[:, good] / tot[good]).T)
        centers.append(0.5 * (edges[i] + edges[i + 1]))
        profs.append(p)
        fwhms.append(profile_fwhm(p))
    if not profs:
        flat = np.full((n_rows, nx), 1.0 / n_rows)
        return flat, np.array([]), np.array([])
    centers = np.array(centers)
    profs = np.array(profs)                       # (n_good_chunks, rows)
    x = np.arange(nx)
    pmap = np.empty((n_rows, nx))
    for r in range(n_rows):
        pmap[r] = np.interp(x, centers, profs[:, r])
    pmap /= pmap.sum(axis=0, keepdims=True)
    return pmap, np.array(fwhms), centers


def profile_fwhm(prof: np.ndarray) -> float:
    """Full width at half maximum of a 1-D profile, by linear
    interpolation of the half-maximum crossings either side of the peak
    (px).  NaN when a crossing falls outside the window."""
    i = int(np.argmax(prof))
    half = prof[i] / 2.0
    lo = i
    while lo > 0 and prof[lo] > half:
        lo -= 1
    hi = i
    while hi < len(prof) - 1 and prof[hi] > half:
        hi += 1
    if prof[lo] > half or prof[hi] > half:
        return float("nan")
    x_lo = lo + (half - prof[lo]) / (prof[lo + 1] - prof[lo])
    x_hi = hi - (half - prof[hi]) / (prof[hi - 1] - prof[hi])
    return float(x_hi - x_lo)


# --------------------------------------------------------------------------
# The estimator
# --------------------------------------------------------------------------
def horne_column(window: np.ndarray, background: np.ndarray,
                 profile: np.ndarray, detector: Detector):
    """Optimal flux for ONE column window (Horne 1986 eq. 8) — the scalar
    reference implementation the vectorised path is tested against.

    Variance per pixel from ``detector.variance_adu2`` (read noise plus
    photon noise of everything above the pedestal).  Pixels at or above
    ``detector.saturation_cap_adu()`` are masked.  Returns
    (flux, variance, n_saturated); flux is None when every pixel is
    masked.
    """
    sat = window >= detector.saturation_cap_adu()
    good = ~sat & np.isfinite(window)
    if not good.any():
        return None, None, int(sat.sum())
    d = window - background                        # net counts, ADU
    var = detector.variance_adu2(window)
    p, v = profile[good], var[good]
    denom = float((p * p / v).sum())
    if denom <= 0:
        return None, None, int(sat.sum())
    flux = float((p * d[good] / v).sum() / denom)
    return flux, 1.0 / denom, int(sat.sum())


def horne_extract(win: np.ndarray, sky: np.ndarray, pmap: np.ndarray,
                  detector: Detector, hot: Optional[np.ndarray] = None,
                  clip_sigma: float = HORNE_CLIP_SIGMA,
                  iters: int = HORNE_ITERS) -> dict:
    """Vectorised Horne extraction of an aperture window (rows x nx).

    ``win`` holds RAW stored counts (pedestal included — the variance
    model needs them), ``sky`` the sky model at the same pixels (also
    raw-count scale, i.e. including the pedestal), ``pmap`` the spatial
    profile.  ``hot`` is an optional boolean mask of known-bad pixels.

    Pass 1 weights with the data-based variance; later passes use the
    MODEL variance (f*P + sky) and reject outliers against it.  Returns
    flux, var, n_sat (per column), n_rej (per column) and the final mask.
    """
    cap = detector.saturation_cap_adu()
    finite = np.isfinite(win) & np.isfinite(sky)
    # Saturation is judged on pixels that are NOT known-bad: a hot pixel
    # at the clip is a detector defect, masked below, not evidence that
    # the target saturated (standing rule 4).
    sat = finite & (win >= cap)
    if hot is not None:
        sat &= ~hot
    mask = finite & ~sat
    if hot is not None:
        mask &= ~hot
    net = np.where(mask, win - sky, 0.0)
    var = detector.variance_adu2(np.where(finite, win, 0.0))
    flux = np.full(win.shape[1], np.nan)
    fvar = np.full(win.shape[1], np.nan)
    rejected = np.zeros_like(mask)
    for it in range(iters + 1):
        w = np.where(mask, 1.0 / var, 0.0)
        denom = (pmap * pmap * w).sum(axis=0)
        ok = denom > 0
        flux = np.where(ok, (pmap * net * w).sum(axis=0)
                        / np.where(ok, denom, 1.0), np.nan)
        fvar = np.where(ok, 1.0 / np.where(ok, denom, 1.0), np.nan)
        if it == iters:
            break
        # Model variance and outlier test for the next pass.
        f = np.nan_to_num(flux)
        model = f[None, :] * pmap
        var = detector.variance_adu2(model + np.where(finite, sky, 0.0))
        tol = (HORNE_PROFILE_TOL * model) ** 2
        chi2 = np.where(mask, (net - model) ** 2 / (var + tol), 0.0)
        # Reject only the single worst pixel per column per pass (Horne's
        # prescription): a genuine outlier distorts f, which inflates its
        # neighbours' residuals; removing them all at once over-rejects.
        worst = np.argmax(chi2, axis=0)
        cols = np.arange(win.shape[1])
        bad = chi2[worst, cols] > clip_sigma ** 2
        if not bad.any():
            # Converged: one more pass with the model variance in place.
            if it > 0:
                break
            continue
        mask[worst[bad], cols[bad]] = False
        rejected[worst[bad], cols[bad]] = True
        net = np.where(mask, win - sky, 0.0)
    return {"flux": flux, "var": fvar,
            "n_sat": sat.sum(axis=0).astype(np.int32),
            "n_rej": rejected.sum(axis=0).astype(np.int32),
            "mask": mask}


def extract_spectrum(data: np.ndarray, trace_coeffs: np.ndarray,
                     detector: Detector,
                     dark: Optional[np.ndarray] = None,
                     sky: str = "poly",
                     template: Optional[np.ndarray] = None,
                     hot: Optional[np.ndarray] = None,
                     halfwin: int = APERTURE_HALFWIN,
                     row_offset: float = 0.0,
                     null_guard: int = SKY_INNER) -> dict:
    """Extract the full spectrum along a fitted trace.

    Parameters
    ----------
    data : raw frame (stored ADU, pedestal included).
    trace_coeffs : polynomial of the trace center, ``np.polyval`` order.
    detector : the :class:`~macro_grism.config.Detector` for this frame.
    dark : optional master dark (same shape, pedestal included).  When
        given, ``dark - pedestal`` is subtracted first — the comparison
        arm for the calibration-debt statistic.  The pedestal is kept in
        the working frame so the variance model stays valid.
    sky : 'poly' (default) or 'flanking' (legacy straight line).
    template : optional 2-D sky model (stored ADU above pedestal) to
        subtract before the local sky fit — the lozenge template.
    hot : optional boolean bad-pixel mask (frame shape).
    halfwin : aperture half-width in px.
    row_offset : shift the whole aperture this many px in y.  Used ONLY by
        the null test (extract an empty aperture beside the trace).  The
        rows within ``null_guard`` px of the REAL trace are then excluded
        from the sky fit, so the null aperture's sky is interpolated
        across a gap exactly as the real aperture's is.

    Returns a dict of 1-D arrays over columns: ``flux`` (optimal, ADU),
    ``var`` (ADU^2), ``box`` (plain aperture sum), ``box_var``, ``bg``
    (sky level at trace center, pedestal included), ``n_sat``, ``n_rej``,
    ``peak`` (brightest raw pixel in the aperture, known-bad pixels
    excluded — the saturation evidence), ``n_hot`` (known-bad pixels
    masked per column), plus ``fwhm_px`` / ``fwhm_x`` (cross-dispersion profile
    FWHM per chunk) and ``n_extracted``.  Unextractable columns are NaN.
    """
    if sky not in ("poly", "flanking"):
        raise ValueError(f"unknown sky model {sky!r}")
    work = np.asarray(data, dtype=np.float64)
    if dark is not None:
        work = work - (dark - detector.pedestal_adu)
    if template is not None:
        work = work - template
    ny, nx = work.shape
    yc = np.polyval(trace_coeffs, np.arange(nx)) + row_offset
    cut = rectify(work, yc)
    raw_cut = rectify(np.asarray(data, dtype=np.float64), yc)
    hot_cut = rectify(hot.astype(float), yc) > 0.5 if hot is not None \
        else None
    half = cut.shape[0] // 2
    exclude = None
    if row_offset:
        dy = np.arange(-half, half + 1)
        exclude = np.abs(dy + row_offset) <= null_guard
    skymod = (sky_polynomial(cut, exclude=exclude) if sky == "poly"
              else sky_flanking(cut))
    a0, a1 = half - halfwin, half + halfwin + 1
    win, sky_win = cut[a0:a1], skymod[a0:a1]
    net = win - sky_win
    if row_offset:
        # Null aperture: there is no trace to take a profile from, so
        # weight uniformly — the null test is about the SKY model.
        pmap = np.full(win.shape, 1.0 / win.shape[0])
        fwhms, centers = np.array([]), np.array([])
    else:
        pmap, fwhms, centers = profile_map(net)
    res = horne_extract(win, sky_win, pmap, detector,
                        hot=None if hot_cut is None else hot_cut[a0:a1])
    # Plain aperture sum and its variance — the cross-check.
    valid = np.isfinite(net).all(axis=0)
    box = np.where(valid, np.nansum(net, axis=0), np.nan)
    box_var = np.where(valid, np.nansum(
        detector.variance_adu2(np.nan_to_num(raw_cut[a0:a1])), axis=0),
        np.nan)
    with _quiet():
        peak_px = raw_cut[a0:a1]
        if hot_cut is not None:
            peak_px = np.where(hot_cut[a0:a1], np.nan, peak_px)
        peak = np.nanmax(peak_px, axis=0)
    flux = np.where(valid, res["flux"], np.nan)
    return {"flux": flux, "var": np.where(valid, res["var"], np.nan),
            "box": box, "box_var": box_var,
            "bg": np.where(valid, sky_win[halfwin], np.nan),
            "n_sat": res["n_sat"], "n_rej": res["n_rej"], "peak": peak,
            "n_hot": (hot_cut[a0:a1].sum(axis=0).astype(np.int32)
                      if hot_cut is not None
                      else np.zeros(win.shape[1], dtype=np.int32)),
            "profile": pmap, "fwhm_px": fwhms, "fwhm_x": centers,
            "n_extracted": int(np.isfinite(flux).sum())}


#: Row offsets (px) of the empty "null" apertures of the sky-model test,
#: and the column block of the sky-variance measurement.  +/-200 px is
#: beyond the star's halo (invisible against the sky past ~30 px on the
#: 240 s frames) and beyond the nearest instrumental ghosts' rows as seen
#: from the null aperture itself, while still inside the field-stop
#: lozenge whose curvature the sky model has to follow.
NULL_OFFSETS = (-200, 200)
SKYVAR_BLOCK = 100


def null_aperture(data: np.ndarray, trace_coeffs: np.ndarray,
                  detector: Detector, sky: str, offset: int,
                  halfwin: int = APERTURE_HALFWIN) -> np.ndarray:
    """Plain aperture sum of an EMPTY aperture ``offset`` px from the
    trace, with the sky modelled exactly as for a real extraction (same
    flanks, same clipping, same smoothing — centred on the null
    aperture).  A perfect sky model returns zero in every column; the
    median of what it does return, as a fraction of the target's
    continuum, is the sky-method bias (G-4's "no negative continuum").

    If the real trace lies inside the null aperture's sky flanks its
    rows are excluded from the fit.  Returns the per-column sums (NaN
    where the window leaves the frame).
    """
    ny, nx = data.shape
    yc = np.polyval(trace_coeffs, np.arange(nx)) + offset
    cut = rectify(np.asarray(data, dtype=np.float64), yc)
    half = cut.shape[0] // 2
    exclude = None
    if abs(offset) <= SKY_OUTER + SKY_INNER:
        dy = np.arange(-half, half + 1)
        exclude = np.abs(dy + offset) <= SKY_INNER
    model = (sky_polynomial(cut, exclude=exclude) if sky == "poly"
             else sky_flanking(cut))
    a0, a1 = half - halfwin, half + halfwin + 1
    net = cut[a0:a1] - model[a0:a1]
    valid = np.isfinite(net).all(axis=0)
    return np.where(valid, np.nansum(net, axis=0), np.nan)


def clipped_variance(r: np.ndarray, nsig: float = 4.0,
                     iters: int = 5) -> float:
    """Variance of a residual sample after iterative ``nsig`` clipping
    about the median.

    Not the MAD: on dark grism flanks the sky residual has a standard
    deviation of ~2 ADU on INTEGER data, and the MAD of a near-integer
    distribution is quantised (it collapses to ~1 ADU and under-estimates
    a sigma of 2 by half — the variance by four).  The clipped standard
    deviation has no such quantisation; 4-sigma clipping removes the
    ghost traces and hot pixels and biases a Gaussian's variance by only
    0.1%."""
    keep = np.isfinite(r)
    for _ in range(iters):
        med = np.median(r[keep])
        sd = np.std(r[keep])
        new = np.isfinite(r) & (np.abs(r - med) <= nsig * sd)
        if (new == keep).all():
            break
        keep = new
    return float(np.var(r[keep]))


def sky_variance_bins(data: np.ndarray, trace_coeffs: np.ndarray,
                      detector: Detector,
                      block: int = SKYVAR_BLOCK) -> np.ndarray:
    """Measured vs predicted pixel variance of the SKY in the flanks of
    the trace, in column blocks — the G-2 acceptance test.

    For each block of ``block`` columns: the flank pixels' residuals from
    the polynomial sky model give a robust variance (iteratively
    4-sigma-clipped; see :func:`clipped_variance` for why not the
    MAD); the model's median level above the pedestal gives the
    variance the detector record predicts.  Returns an array of rows
    [level_adu, var_measured, var_predicted, n_pixels].

    The polynomial (3 parameters per column, smoothed over 31 columns)
    removes < 2% of the white-noise variance from ~200 flank pixels per
    column, so the measurement is not biased low at the 20% level the
    test is run at.
    """
    ny, nx = data.shape
    yc = np.polyval(trace_coeffs, np.arange(nx))
    cut = rectify(np.asarray(data, dtype=np.float64), yc)
    half = cut.shape[0] // 2
    model = sky_polynomial(cut)
    dy = np.abs(np.arange(-half, half + 1))
    flank = (dy >= SKY_INNER) & (dy <= SKY_OUTER)
    resid = (cut - model)[flank]
    level = model[flank] - detector.pedestal_adu
    rows = []
    for x0 in range(0, nx - block + 1, block):
        r = resid[:, x0:x0 + block].ravel()
        lv = level[:, x0:x0 + block].ravel()
        ok = np.isfinite(r) & np.isfinite(lv)
        if ok.sum() < 0.8 * r.size:
            continue
        r, lv = r[ok], lv[ok]
        var = clipped_variance(r)
        lev = float(np.median(lv))
        pred = float(detector.read_noise_adu ** 2
                     + max(lev, 0.0) / detector.gain_e_per_adu)
        rows.append([lev, float(var), pred, int(ok.sum())])
    return np.array(rows).reshape(-1, 4)


def median_relative_difference(flux_a: np.ndarray,
                               flux_b: np.ndarray) -> Optional[float]:
    """The method-difference statistic: median over columns of
    |A - B| / |A|, restricted to columns where both spectra exist and A
    has meaningful signal (|A| above the 25th percentile of |A| — ratio
    against near-zero flux measures nothing but noise).

    Used for flanking-vs-master-dark (the calibration debt) and for
    straight-line-vs-polynomial sky (the G-4 method difference).
    """
    both = np.isfinite(flux_a) & np.isfinite(flux_b)
    if both.sum() < 100:
        return None
    a, b = flux_a[both], flux_b[both]
    floor = np.percentile(np.abs(a), 25)
    sig = np.abs(a) >= max(floor, 1e-9)
    if sig.sum() < 50:
        return None
    return float(np.median(np.abs(a[sig] - b[sig]) / np.abs(a[sig])))
