"""Pure characterization arithmetic — image quality, noise, cadence, detectability.

This module is the S5 ("what can this data actually do?") counterpart to
``macro_phot.errors``.  Everything here is a pure function of numbers: no
database, no FITS, no matplotlib.  The staged CLI
``pipeline/scripts/run_cv_characterization.py`` supplies the numbers, stores
the answers, and the renderer draws them.  Unit tests live in
``pipeline/tests/test_characterize.py``.

The four questions this file answers, in the order the evidence has to be
built (each answer is the input to the next):

1.  **How good are the images?**  FWHM must be quoted in ARCSEC (three
    plate scales are in play: 0.5375, 0.8062 and 0.4491 "/px), sky level in
    ADU per pixel per second so eras with different exposure times compare,
    airmass recomputed from coordinates and time (the archive's AIRMASS
    headers contain values up to 6877 — they are not data), and moon
    proximity from ephemerides rather than the MOONANGL card.  The cut that
    separates "usable" from "degraded" is then chosen by its measured
    CONSEQUENCE: bin frames by a quality axis, measure the check-star
    scatter in each bin, and cut where the scatter has degraded past a
    stated factor of the best-bin baseline.

2.  **What is the noise floor?**  The measured RMS-vs-magnitude cloud of
    constant stars versus a photon+sky+read prediction computed from the S2
    detector facts.  The gain enters only the SOURCE shot term (the sky and
    read terms come from the per-frame measured background RMS), so the
    honest [0.60, 1.77] e-/ADU bracket produces a prediction BAND, not a
    line.  Where the measured cloud stops following the prediction and goes
    flat, that flat value is the systematic floor; :func:`fit_noise_floor`
    measures it, and :func:`scintillation_young` says whether the
    atmosphere can explain it (on a 0.5 m telescope it cannot).

3.  **What can the sampling resolve?**  The spectral window
    ``|sum exp(-2 pi i f t)|^2 / N^2`` and its alias comb at 1 cycle/day
    (and 1 sidereal day, and the nightly-block harmonics).  A period whose
    aliases carry comparable window power cannot be claimed from this
    sampling alone.

4.  **What amplitude is detectable?**  Injection and recovery through the
    REAL timestamps and REAL noise (check-star residuals, cyclically
    shifted so their correlated structure survives), with the detection
    threshold set by the same machinery run on signal-free data.

House rule followed throughout: a function either returns a number a
student can re-derive from its docstring, or it does not belong here.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Detector facts, quoted from the S2 report (docs/pipeline/s2_detector.html).
# These are MEASURED values; the report interpolates them so the page and
# the arithmetic can never disagree.
# --------------------------------------------------------------------------

#: High Gain single-read noise, ADU (S2 §2: variance floor of the shortest
#: darks, dark-current term removed).
READ_NOISE_ADU_HIGH_GAIN = 3.93

#: The same read noise in electrons at the nominal header gain, with the
#: uncertainty that the GAIN BRACKET (not the ADU statistics) implies.
READ_NOISE_E = 4.15
READ_NOISE_E_ERR = 2.30

#: Gain bracket, e-/ADU.  Lower bound = sky-level light slope, upper bound =
#: dark slope; the header EGAIN sits inside and is adopted as nominal.
#: NEVER quote a single gain as if it were measured.
GAIN_LO_E_PER_ADU = 0.60
GAIN_NOMINAL_E_PER_ADU = 1.057
GAIN_HI_E_PER_ADU = 1.77

#: StackPro is the SUM of this many sub-reads (three independent ratios
#: agree).  Read-noise variance is N_sub x the single read.
STACKPRO_N_SUB = 16

#: Telescope aperture diameter (cm) and site altitude (m) for the
#: scintillation prediction.  RLMT is a 0.5 m at Winer Observatory.
TELESCOPE_DIAMETER_CM = 50.0
SITE_ALTITUDE_M = 1515.0

#: Young's (1967/1993) scintillation coefficient and scale height, the
#: standard form quoted by Dravins et al. 1998 and Osborn et al. 2015.
YOUNG_COEFF = 0.09
ATMOSPHERIC_SCALE_HEIGHT_M = 8000.0

#: 2.5 / ln 10 — relative flux error to magnitude error.
MAG_ERR_FACTOR = 2.5 / math.log(10.0)

# --------------------------------------------------------------------------
# Policy constants for the cuts and grids (single source of truth; the
# report interpolates every one of them).
# --------------------------------------------------------------------------

#: A quality bin counts as DEGRADED when its check-star scatter exceeds the
#: best decile's scatter by this factor.  1.30 = "30% worse than this
#: series can do" — chosen because it is well outside the bin-to-bin
#: sampling noise of a 4-check-star median (which is ~10%) yet still throws
#: away frames whose information content is materially reduced.
DEGRADE_FACTOR = 1.30

#: Minimum frames in a quality bin before its scatter is allowed to set a
#: threshold (a 3-frame bin's median is not a measurement).
MIN_FRAMES_PER_QUALITY_BIN = 15

#: Sidereal day in days — the OTHER alias spacing, 0.27% from 1 c/d; a
#: baseline longer than ~1 yr separates them, shorter does not.
SIDEREAL_DAY_D = 0.9972695787

#: False-alarm probability the detection threshold is set at (bootstrap,
#: not analytic: correlated noise breaks the analytic Baluev/Scargle FAP).
DETECT_FAP = 0.001

#: A recovered period counts as CORRECT when it lands within this
#: fractional distance of the injected one.  1% is far tighter than the
#: 1 c/d alias spacing at these periods, so "correct" really does mean the
#: right peak and not its neighbour.
PERIOD_TOL_FRAC = 0.01

#: Recovery contour level.
RECOVERY_LEVEL = 0.90


# ==========================================================================
# 1 - IMAGE QUALITY
# ==========================================================================

def fwhm_arcsec(fwhm_px, plate_scale_arcsec_per_px) -> np.ndarray:
    """Seeing in ARCSEC from a pixel FWHM and that frame's plate scale.

    The whole point of this one-line function is that it is the ONLY place
    the conversion happens: the archive mixes 0.5375 (High Gain 4096),
    0.8062 (Andor iKon) and 0.4491 "/px (2x2-binned CMOS), so a pixel FWHM
    is not comparable across eras and an arcsec FWHM is.  Non-finite or
    non-positive inputs propagate as NaN rather than as a plausible-looking
    small number.
    """
    f = np.asarray(fwhm_px, dtype=float)
    s = np.asarray(plate_scale_arcsec_per_px, dtype=float)
    out = f * s
    bad = ~np.isfinite(out) | (f <= 0) | (s <= 0)
    out = np.where(bad, np.nan, out)
    return out


def sky_rate_adu_per_px_s(bkg_adu, exptime_s) -> np.ndarray:
    """Sky background as ADU per pixel per SECOND.

    The raw background level is meaningless across an archive whose
    exposures run 8 s to 300 s; dividing by exposure time makes a moonlit
    240 s frame comparable with a dark 8 s one.  (It does NOT make two
    cameras comparable — the gain and pixel solid angle differ — so the
    report always groups this by era.)
    """
    b = np.asarray(bkg_adu, dtype=float)
    t = np.asarray(exptime_s, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = b / t
    return np.where(np.isfinite(out) & (t > 0), out, np.nan)


def airmass_from_altitude(alt_deg) -> np.ndarray:
    """Kasten & Young (1989) airmass from apparent altitude in degrees.

        X = 1 / ( sin(h) + 0.50572 * (h + 6.07995)^-1.6364 )     [h in deg]

    Chosen over the plane-parallel sec(z) because VV Pup never rises above
    ~39 deg altitude from Winer, where sec(z) is already 1% wrong and
    climbing.  Targets below the horizon return NaN — a frame taken with
    the target below the horizon is a pointing bug, not a high airmass, and
    must not be silently rendered as a large number (the archive's own
    AIRMASS headers reach 6877 precisely because nothing refused).
    """
    h = np.asarray(alt_deg, dtype=float)
    with np.errstate(invalid="ignore"):
        denom = (np.sin(np.radians(h))
                 + 0.50572 * np.power(np.clip(h + 6.07995, 1e-9, None),
                                      -1.6364))
        x = 1.0 / denom
    return np.where(np.isfinite(x) & (h > 0), x, np.nan)


def moon_illuminated_fraction(elongation_deg) -> np.ndarray:
    """Illuminated fraction of the lunar disc from the Sun-Moon elongation.

        k = (1 - cos(elongation)) / 2

    New moon (elongation 0) -> 0, full moon (180 deg) -> 1.  Computed from
    ephemeris positions rather than read from the MOONPHAS header, for the
    same reason airmass is recomputed: the headers are unaudited and, where
    they have been audited, wrong.
    """
    e = np.asarray(elongation_deg, dtype=float)
    return (1.0 - np.cos(np.radians(e))) / 2.0


def binned_median(x, y, edges) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median of ``y`` inside bins of ``x`` given explicit bin ``edges``.

    Returns ``(centers, medians, counts)`` for every bin, with NaN medians
    where a bin is empty.  Explicit edges (rather than a bin count) because
    the thresholds the report defends are read off these bins and must be
    reproducible from the edges alone.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    edges = np.asarray(edges, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    centers = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(centers.size, np.nan)
    cnt = np.zeros(centers.size, dtype=int)
    idx = np.digitize(x, edges) - 1
    for b in range(centers.size):
        sel = idx == b
        cnt[b] = int(sel.sum())
        if cnt[b]:
            med[b] = float(np.median(y[sel]))
    return centers, med, cnt


def degradation_threshold(centers, medians, counts,
                          factor: float = DEGRADE_FACTOR,
                          min_count: int = MIN_FRAMES_PER_QUALITY_BIN
                          ) -> tuple[float, float]:
    """Where a quality axis starts costing precision, and the baseline.

    The baseline is the SMALLEST well-populated bin median (what this
    series achieves when the axis is favourable).  Walking outward from the
    baseline bin toward larger ``x``, the threshold is the left edge of the
    first well-populated bin whose median exceeds ``factor * baseline``.
    Returns ``(threshold_x, baseline_y)``; the threshold is ``inf`` when no
    bin ever degrades that far — which is itself a result ("this axis does
    not limit us"), not a failure.

    This is the whole method for defending a cut: the number comes from the
    data's own response, so the only judgement call left in the open is
    ``factor``, which the report prints.
    """
    centers = np.asarray(centers, dtype=float)
    medians = np.asarray(medians, dtype=float)
    counts = np.asarray(counts, dtype=int)
    good = (counts >= min_count) & np.isfinite(medians)
    if not good.any():
        return float("inf"), float("nan")
    base_i = int(np.flatnonzero(good)[np.argmin(medians[good])])
    baseline = float(medians[base_i])
    for i in range(base_i + 1, centers.size):
        if good[i] and medians[i] > factor * baseline:
            return float(centers[i]), baseline
    return float("inf"), baseline


def usable_mask(values: dict[str, np.ndarray],
                thresholds: dict[str, float]) -> np.ndarray:
    """Boolean 'this frame is usable' from per-axis upper thresholds.

    ``values`` maps axis name -> per-frame value, ``thresholds`` maps the
    same names -> the upper limit from :func:`degradation_threshold`.  A
    frame passes when EVERY axis present in both dicts is finite and at or
    below its threshold; a NaN on any axis fails (an unmeasurable frame is
    not a usable frame).  Axes with an infinite threshold pass everything,
    so an axis that never degrades cannot silently veto anything.
    """
    n = len(next(iter(values.values()))) if values else 0
    ok = np.ones(n, dtype=bool)
    for name, thr in thresholds.items():
        if name not in values:
            continue
        v = np.asarray(values[name], dtype=float)
        ok &= np.isfinite(v) & (v <= thr)
    return ok


def ellipticity(a, b) -> np.ndarray:
    """Ellipticity ``1 - b/a`` from sep's semi-axes (0 = round, 1 = a line).

    Trailing from a guiding failure or a wind gust shows up here AND in the
    position-angle coherence of the same detections: real seeing elongation
    is random star-to-star, trailing is not.  :func:`pa_coherence` is the
    companion test.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        e = 1.0 - b / a
    return np.where(np.isfinite(e) & (a > 0), e, np.nan)


def pa_coherence(theta_rad) -> float:
    """How aligned a frame's source elongations are: 0 = random, 1 = one PA.

    Position angles are defined modulo pi (an ellipse has no head or tail),
    so the circular mean is taken on the DOUBLED angle:

        R = | mean( exp(2 i theta) ) |

    A frame of round stars measured by sep produces near-random theta and
    R ~ 1/sqrt(N); a trailed frame drives every star to the same PA and R
    -> 1.  Empty input returns NaN.
    """
    t = np.asarray(theta_rad, dtype=float)
    t = t[np.isfinite(t)]
    if t.size == 0:
        return float("nan")
    return float(np.abs(np.mean(np.exp(2j * t))))


# ==========================================================================
# 2 - NOISE
# ==========================================================================

def predicted_sigma_mag(flux_adu, n_pix_aper: float, bkg_rms_adu,
                        gain_e_per_adu: float,
                        n_sky_pix: Optional[float] = None) -> np.ndarray:
    """Predicted per-point magnitude sigma from first principles.

    The variance of an aperture sum, in ADU^2:

        var = F/g                      source shot noise (F in ADU, g in e-/ADU)
            + n_pix * sigma_bkg^2      sky + read + dark inside the aperture
            + n_pix^2/n_sky * sigma_bkg^2   error of the sky level itself

    ``sigma_bkg`` is the frame's MEASURED background RMS, so the sky, read
    and dark terms are empirical and carry no assumption at all — this is
    why the gain bracket only widens the SOURCE term and the prediction
    stays useful despite an unmeasured gain.  StackPro needs no special
    case for the same reason: its 16-sub-read read noise is already inside
    the measured background RMS.

    The magnitude sigma is ``1.0857 * sqrt(var) / F``.  Non-positive fluxes
    return NaN.
    """
    f = np.asarray(flux_adu, dtype=float)
    s = np.asarray(bkg_rms_adu, dtype=float)
    if n_sky_pix is None:
        # Default annulus: the pipeline's 8-12" annulus around a 4" aperture
        # is (12^2 - 8^2)/4^2 = 5 aperture areas of sky pixels.
        n_sky_pix = 5.0 * n_pix_aper
    with np.errstate(divide="ignore", invalid="ignore"):
        var = (f / float(gain_e_per_adu)
               + n_pix_aper * s ** 2
               + (n_pix_aper ** 2 / float(n_sky_pix)) * s ** 2)
        sig = MAG_ERR_FACTOR * np.sqrt(var) / f
    return np.where(np.isfinite(sig) & (f > 0), sig, np.nan)


def scintillation_young(airmass, exptime_s,
                        diameter_cm: float = TELESCOPE_DIAMETER_CM,
                        altitude_m: float = SITE_ALTITUDE_M) -> np.ndarray:
    """Young's scintillation estimate, in magnitudes (relative-flux ~ mag).

        sigma = 0.09 * D^(-2/3) * X^1.75 * exp(-h/8000) * (2 T)^(-1/2)

    with D in cm, X the airmass, h the site altitude in m and T the
    exposure in seconds.  Quoted because it is the ONE noise term that no
    reduction can remove and that scales with nothing the pipeline
    controls: if the measured systematic floor sits far above this, the
    floor is instrumental (flat field, ensemble solution, colour terms) and
    is therefore in principle fixable.  On a 0.5 m at 1515 m in a 60 s
    exposure at X = 1.3 the answer is ~0.8 mmag — an order of magnitude
    below anything measured here, which is the point.
    """
    x = np.asarray(airmass, dtype=float)
    t = np.asarray(exptime_s, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        sig = (YOUNG_COEFF
               * diameter_cm ** (-2.0 / 3.0)
               * np.power(x, 1.75)
               * math.exp(-altitude_m / ATMOSPHERIC_SCALE_HEIGHT_M)
               * np.power(2.0 * t, -0.5))
    return np.where(np.isfinite(sig) & (x >= 1) & (t > 0), sig, np.nan)


def fit_noise_floor(mag, rms, predicted) -> tuple[float, float, int]:
    """Split a measured RMS-vs-mag cloud into a scaled prediction + a floor.

    Model, fitted in VARIANCE (where the terms add):

        rms^2 = k * predicted^2 + floor^2

    solved by non-negative least squares on the two-column design matrix
    ``[predicted^2, 1]``.  ``k`` is the factor by which the formal photon
    model is wrong (k > 1 means the predicted errors are optimistic even
    before any systematic), and ``floor`` is the magnitude-independent
    scatter that survives at infinite brightness: scintillation + flat-field
    residual + ensemble zero-point error + any real low-level variability of
    the "constant" stars.

    Returns ``(floor_mag, k, n_used)``; NaNs when fewer than three stars are
    usable.  Negative solutions are clipped to zero before the square root,
    so a floor is never reported as an imaginary number.
    """
    m = np.asarray(mag, dtype=float)
    r = np.asarray(rms, dtype=float)
    p = np.asarray(predicted, dtype=float)
    ok = np.isfinite(m) & np.isfinite(r) & np.isfinite(p) & (r > 0) & (p > 0)
    r, p = r[ok], p[ok]
    if r.size < 3:
        return float("nan"), float("nan"), int(r.size)
    A = np.column_stack([p ** 2, np.ones_like(p)])
    coef, *_ = np.linalg.lstsq(A, r ** 2, rcond=None)
    k = max(float(coef[0]), 0.0)
    floor = math.sqrt(max(float(coef[1]), 0.0))
    return floor, k, int(r.size)


def precision_at_mag(mag, rms, target_mag: float,
                     half_width: float = 0.5) -> tuple[float, int]:
    """Achieved per-point precision at a given brightness, measured.

    The median RMS of the constant stars within ``half_width`` magnitudes
    of ``target_mag``.  This is the honest way to state "what precision do
    we get on the target": the target itself is variable and cannot
    validate its own error bars, so its brightness is used only to pick
    which constant stars answer the question.  Returns
    ``(rms_mag, n_stars)``; ``(nan, 0)`` when no star is in range, which
    happens exactly when the target is brighter or fainter than every
    comparison star in its field — a result worth printing.
    """
    m = np.asarray(mag, dtype=float)
    r = np.asarray(rms, dtype=float)
    sel = np.isfinite(m) & np.isfinite(r) & (np.abs(m - target_mag) <= half_width)
    if not sel.any():
        return float("nan"), 0
    return float(np.median(r[sel])), int(sel.sum())


def allan_slope(tau_s, adev) -> float:
    """Log-log slope of an Allan deviation ladder.

    White noise averages down as tau^-1/2, so a slope of -0.5 is the
    licence to average.  A slope shallower than that (toward 0) says
    correlated noise: averaging N points does NOT buy sqrt(N), and any
    error bar computed as sigma/sqrt(N) is a fiction.  Fewer than three
    finite rungs returns NaN.
    """
    t = np.asarray(tau_s, dtype=float)
    a = np.asarray(adev, dtype=float)
    ok = np.isfinite(t) & np.isfinite(a) & (t > 0) & (a > 0)
    if ok.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log10(t[ok]), np.log10(a[ok]), 1)[0])


def red_noise_factor(tau_s, adev, tau_target_s: float) -> float:
    """How much worse than white the noise is at one averaging timescale.

    Extrapolates the FIRST rung (tau_0, sigma_0) forward as white noise,
    ``sigma_white = sigma_0 * sqrt(tau_0 / tau_target)``, and divides the
    measured Allan deviation at (or nearest below) ``tau_target`` by it.
    A factor of 1 means averaging works perfectly; 2 means a mean over that
    timescale is twice as uncertain as the per-point scatter implies.

    ``tau_target`` is chosen by the science: for these polars it is the
    orbital period, because every bright-phase timing and every folded
    colour point is an average over about that long.
    """
    t = np.asarray(tau_s, dtype=float)
    a = np.asarray(adev, dtype=float)
    ok = np.isfinite(t) & np.isfinite(a) & (t > 0) & (a > 0)
    t, a = t[ok], a[ok]
    if t.size < 2:
        return float("nan")
    order = np.argsort(t)
    t, a = t[order], a[order]
    at_or_below = np.flatnonzero(t <= tau_target_s)
    if at_or_below.size == 0:
        return float("nan")
    i = int(at_or_below[-1])
    white = a[0] * math.sqrt(t[0] / t[i])
    if white <= 0:
        return float("nan")
    return float(a[i] / white)


# ==========================================================================
# 3 - CADENCE AND SAMPLING
# ==========================================================================

def night_blocks(times_d, max_gap_h: float = 3.0) -> list[tuple[int, int]]:
    """Split a sorted time array into contiguous observing blocks.

    Returns ``(start, stop)`` slice index pairs.  A gap longer than
    ``max_gap_h`` starts a new block; 3 h is longer than any within-night
    pause in this archive and far shorter than the ~14 h night-to-night
    gap, so blocks are nights (or the separate halves of a night broken by
    weather) without needing a calendar.
    """
    t = np.asarray(times_d, dtype=float)
    if t.size == 0:
        return []
    gaps = np.diff(t) * 24.0
    breaks = np.flatnonzero(gaps > max_gap_h)
    starts = np.concatenate(([0], breaks + 1))
    stops = np.concatenate((breaks + 1, [t.size]))
    return [(int(a), int(b)) for a, b in zip(starts, stops)]


def cadence_stats(times_d, period_d: Optional[float] = None,
                  max_gap_h: float = 3.0) -> dict:
    """Everything the sampling section needs about one series, in one pass.

    Returns a dict with: number of blocks (nights), total baseline in days,
    median within-block sampling interval in seconds, longest block in
    hours, the duty cycle (time actually on target / baseline), and — when
    a period is supplied — points per cycle, cycles in the longest block,
    and the fraction of the orbital phase circle that the whole series
    covers at all (in 20 phase bins).  Phase coverage is the number that
    decides whether a fold is a light curve or a rumour.
    """
    t = np.sort(np.asarray(times_d, dtype=float))
    t = t[np.isfinite(t)]
    out = {"n_points": int(t.size)}
    if t.size == 0:
        return out
    blocks = night_blocks(t, max_gap_h)
    out["n_blocks"] = len(blocks)
    out["baseline_d"] = float(t[-1] - t[0])
    within = np.concatenate([np.diff(t[a:b]) for a, b in blocks
                             if b - a > 1] or [np.array([])])
    out["median_dt_s"] = float(np.median(within) * 86400.0) if within.size else float("nan")
    block_len = np.array([t[b - 1] - t[a] for a, b in blocks])
    out["longest_block_h"] = float(block_len.max() * 24.0)
    out["on_target_h"] = float(block_len.sum() * 24.0)
    out["duty_cycle"] = (float(block_len.sum() / out["baseline_d"])
                         if out["baseline_d"] > 0 else float("nan"))
    if period_d:
        out["period_d"] = float(period_d)
        out["pts_per_cycle"] = (float(period_d * 86400.0 / out["median_dt_s"])
                                if np.isfinite(out["median_dt_s"]) and out["median_dt_s"] > 0
                                else float("nan"))
        out["cycles_longest_block"] = float(block_len.max() / period_d)
        ph = np.mod(t / period_d, 1.0)
        filled = np.unique((ph * 20).astype(int)).size
        out["phase_coverage"] = float(filled / 20.0)
    return out


def spectral_window(times_d, freqs_cd) -> np.ndarray:
    """Spectral window ``W(f) = |sum exp(-2 pi i f t)|^2 / N^2``.

    This is the periodogram the sampling alone would produce for a
    perfectly constant star observed at these instants, normalised so
    ``W(0) = 1``.  Every peak in it is a frequency at which the sampling
    manufactures power; every peak in a real periodogram is the true
    signal CONVOLVED with this.  Reading it before reading any Lomb-Scargle
    is the difference between a period and an alias.
    """
    t = np.asarray(times_d, dtype=float)
    t = t[np.isfinite(t)]
    f = np.asarray(freqs_cd, dtype=float)
    if t.size == 0:
        return np.zeros_like(f)
    t = t - t.mean()
    phase = np.exp(-2j * np.pi * np.outer(f, t))
    return np.abs(phase.sum(axis=1)) ** 2 / t.size ** 2


def alias_ladder(f_true_cd: float, f_alias_cd: float = 1.0,
                 orders: Sequence[int] = (-3, -2, -1, 1, 2, 3)) -> list[float]:
    """The frequencies a true signal is confusable with: ``f +/- k f_alias``.

    Negative results are folded back to their absolute value (a periodogram
    cannot tell +f from -f) and frequencies at or below zero are dropped.
    """
    out = []
    for k in orders:
        f = f_true_cd + k * f_alias_cd
        if abs(f) > 1e-9:
            out.append(abs(float(f)))
    return out


def alias_power(times_d, f_true_cd: float, f_alias_cd: float = 1.0,
                orders: Sequence[int] = (-3, -2, -1, 1, 2, 3)
                ) -> list[tuple[int, float, float]]:
    """Window power carried by each alias of a true frequency.

    Returns ``(order, frequency, window_power)`` per alias, where the window
    is evaluated at the OFFSET ``k * f_alias`` — because the periodogram
    peak an alias produces has height set by ``W(f_peak - f_true)``, not by
    ``W(f_peak)``.  A ratio near 1 means the alias is as tall as the truth
    and the period is not determinable from this sampling alone; the
    published number that decides an alias argument is exactly this.
    """
    offs = np.array([k * f_alias_cd for k in orders], dtype=float)
    w = spectral_window(times_d, offs)
    out = []
    for k, off, p in zip(orders, offs, w):
        f = f_true_cd + off
        if abs(f) > 1e-9:
            out.append((int(k), abs(float(f)), float(p)))
    return out


# ==========================================================================
# 4 - DETECTABILITY
# ==========================================================================

def inject_sinusoid(times_d, period_d: float, semi_amp_mag: float,
                    phase: float = 0.0) -> np.ndarray:
    """A sinusoid of the given SEMI-amplitude, sampled at ``times_d``.

    Semi-amplitude (not peak-to-peak) throughout this module, because that
    is what a Lomb-Scargle amplitude spectrum returns and what the
    A_min = sigma sqrt(4 z / N) formula in the strategy document means.
    """
    t = np.asarray(times_d, dtype=float)
    return semi_amp_mag * np.sin(2 * np.pi * (t / period_d + phase))


def eclipse_template(times_d, period_d: float, depth_mag: float,
                     width_phase: float = 0.10, t0_d: float = 0.0,
                     ingress_phase: float = 0.02) -> np.ndarray:
    """A trapezoidal eclipse/bright-phase edge sampled at ``times_d``.

    ``width_phase`` is the full width at half depth in units of the orbital
    phase, ``ingress_phase`` the duration of each sloped edge.  Returned in
    magnitudes as a POSITIVE dip (fainter), so it is added to a magnitude
    series.  A trapezoid rather than a box because the quantity being
    measured downstream is the epoch of the EDGE, and a box's edge has
    infinite slope, which would flatter the timing test into meaninglessness.
    """
    t = np.asarray(times_d, dtype=float)
    ph = np.mod((t - t0_d) / period_d + 0.5, 1.0) - 0.5
    half = width_phase / 2.0
    ing = max(ingress_phase, 1e-6)
    # Piecewise-linear trapezoid: flat bottom inside +/- (half - ing/2),
    # linear ramps over ing, zero outside +/- (half + ing/2).
    a = np.abs(ph)
    out = np.zeros_like(a)
    inner, outer = half - ing / 2.0, half + ing / 2.0
    out[a <= inner] = 1.0
    ramp = (a > inner) & (a < outer)
    out[ramp] = (outer - a[ramp]) / ing
    return depth_mag * out


def cyclic_noise_realization(residuals, shift: int) -> np.ndarray:
    """Roll a residual series by ``shift`` samples: a noise realization that
    KEEPS its correlated structure.

    Shuffling residuals would whiten them and make every detection limit
    look better than the data can deliver.  Rolling preserves the
    autocorrelation exactly (this is the "prayer-bead" / cyclic-permutation
    trick used for red-noise error bars in transit work) while decorrelating
    the noise from the injected signal's phase, which is all the injection
    test needs.
    """
    r = np.asarray(residuals, dtype=float)
    if r.size == 0:
        return r
    return np.roll(r, int(shift) % r.size)


def ls_peak(times_d, y, freqs_cd) -> tuple[float, float]:
    """Highest Lomb-Scargle peak: ``(frequency, power)``.

    Floating-mean (astropy's default) so a constant offset cannot leak into
    the amplitude, evaluated on the caller's explicit frequency grid so
    the search band is a stated policy rather than an astropy default.
    """
    from astropy.timeseries import LombScargle
    t = np.asarray(times_d, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 5:
        return float("nan"), float("nan")
    power = LombScargle(t, y).power(np.asarray(freqs_cd, dtype=float))
    i = int(np.argmax(power))
    return float(freqs_cd[i]), float(power[i])


def detection_threshold(times_d, residual_pool, freqs_cd,
                        n_trials: int = 200,
                        fap: float = DETECT_FAP,
                        rng: Optional[np.random.Generator] = None) -> float:
    """Bootstrap detection threshold: the ``1 - fap`` quantile of the
    highest periodogram peak found in SIGNAL-FREE data.

    ``residual_pool`` is a list of real check-star residual series measured
    at the same timestamps.  Realizations are cyclic rolls of those series,
    so the threshold inherits the archive's real correlated noise; an
    analytic Baluev FAP computed on this data would be optimistic by the
    same factor the noise is red.
    """
    rng = rng or np.random.default_rng(20260819)
    peaks = []
    pool = [np.asarray(r, dtype=float) for r in residual_pool
            if np.isfinite(r).sum() >= 5]
    if not pool:
        return float("nan")
    for _ in range(n_trials):
        base = pool[rng.integers(len(pool))]
        y = cyclic_noise_realization(base, rng.integers(base.size))
        _, p = ls_peak(times_d, y, freqs_cd)
        if np.isfinite(p):
            peaks.append(p)
    if not peaks:
        return float("nan")
    return float(np.quantile(peaks, 1.0 - fap))


def recovery_fraction(times_d, residual_pool, freqs_cd, period_d: float,
                      semi_amp_mag: float, threshold: float,
                      n_trials: int = 60,
                      period_tol: float = PERIOD_TOL_FRAC,
                      rng: Optional[np.random.Generator] = None) -> float:
    """Fraction of injections recovered at the RIGHT period.

    A trial counts as recovered only when the peak power exceeds
    ``threshold`` AND the peak frequency is within ``period_tol`` of the
    injected one.  Requiring both is what makes the contour a statement
    about measuring a period rather than about noticing that something
    varies — an alias-family peak above threshold is a detection of the
    wrong period, and this pipeline scores it as a failure.
    """
    rng = rng or np.random.default_rng(20260819)
    pool = [np.asarray(r, dtype=float) for r in residual_pool
            if np.isfinite(r).sum() >= 5]
    if not pool or not np.isfinite(threshold):
        return float("nan")
    f_true = 1.0 / period_d
    hits = 0
    for _ in range(n_trials):
        base = pool[rng.integers(len(pool))]
        noise = cyclic_noise_realization(base, rng.integers(base.size))
        sig = inject_sinusoid(times_d, period_d, semi_amp_mag,
                              phase=float(rng.random()))
        f, p = ls_peak(times_d, noise + sig, freqs_cd)
        if np.isfinite(p) and p > threshold and abs(f - f_true) <= period_tol * f_true:
            hits += 1
    return hits / float(n_trials)


def recovery_contour(amps, fracs, level: float = RECOVERY_LEVEL) -> float:
    """Smallest amplitude whose recovery fraction reaches ``level``.

    Linear interpolation between the last amplitude below the level and the
    first at or above it, in log-amplitude (recovery curves are much closer
    to straight in log A).  Returns NaN when the grid never reaches the
    level — the honest answer "not detectable anywhere on this grid", which
    the report prints as such rather than extrapolating.
    """
    a = np.asarray(amps, dtype=float)
    f = np.asarray(fracs, dtype=float)
    ok = np.isfinite(a) & np.isfinite(f)
    a, f = a[ok], f[ok]
    if a.size == 0:
        return float("nan")
    order = np.argsort(a)
    a, f = a[order], f[order]
    above = np.flatnonzero(f >= level)
    if above.size == 0:
        return float("nan")
    i = int(above[0])
    if i == 0:
        return float(a[0])
    a0, a1, f0, f1 = a[i - 1], a[i], f[i - 1], f[i]
    if f1 == f0:
        return float(a1)
    w = (level - f0) / (f1 - f0)
    return float(10 ** (np.log10(a0) + w * (np.log10(a1) - np.log10(a0))))


def timing_precision_mc(times_d, period_d: float, depth_mag: float,
                        sigma_mag: float, width_phase: float = 0.10,
                        ingress_phase: float = 0.02,
                        n_trials: int = 200,
                        rng: Optional[np.random.Generator] = None) -> float:
    """Monte-Carlo epoch uncertainty of ONE eclipse/bright-phase edge, seconds.

    Injects a trapezoidal feature at a random epoch into white noise of
    ``sigma_mag`` at the supplied REAL timestamps of one cycle, then
    recovers the epoch by chi-square minimisation over a fine epoch grid
    with the template shape held fixed (the optimistic case: shape known
    exactly).  The returned robust scatter (1.4826 x MAD) of recovered
    minus true epoch is therefore a LOWER BOUND on what a real fit achieves.

    The reason this is the decisive measurement for bright-phase timing:
    when the sampling interval is long compared with the ingress, the
    answer is set by the sampling, not by the photometry, and no amount of
    signal-to-noise improves it.  A pipeline that assumes otherwise reports
    a timing precision it cannot deliver.
    """
    rng = rng or np.random.default_rng(20260819)
    t = np.asarray(times_d, dtype=float)
    t = t[np.isfinite(t)]
    if t.size < 6 or not np.isfinite(sigma_mag) or sigma_mag <= 0:
        return float("nan")
    span = t.max() - t.min()
    # Epoch search grid: 1/400 of a period is far finer than any plausible
    # answer, so grid granularity never limits the result.
    grid = np.linspace(-0.5 * period_d, 0.5 * period_d, 401)
    errs = []
    for _ in range(n_trials):
        t0 = t.min() + span * float(rng.random())
        y = eclipse_template(t, period_d, depth_mag, width_phase, t0,
                             ingress_phase) + rng.normal(0, sigma_mag, t.size)
        chi2 = np.array([
            np.sum((y - eclipse_template(t, period_d, depth_mag, width_phase,
                                         t0 + g, ingress_phase)) ** 2)
            for g in grid])
        errs.append(grid[int(np.argmin(chi2))])
    e = np.asarray(errs)
    mad = float(np.median(np.abs(e - np.median(e))))
    return 1.4826 * mad * 86400.0


def amin_analytic(sigma_mag: float, n_points: int, z: float = 18.0) -> float:
    """The strategy document's analytic detection limit, for comparison only.

        A_min = sigma * sqrt(4 z / N)

    Reproduced here so the report can put it beside the MEASURED injection
    contour and show the gap.  It assumes white, Gaussian, uncorrelated
    noise and an alias-free window: two assumptions this archive violates,
    which is why the measured contour is the number that governs.
    """
    if n_points <= 0 or not np.isfinite(sigma_mag):
        return float("nan")
    return float(sigma_mag * math.sqrt(4.0 * z / n_points))
