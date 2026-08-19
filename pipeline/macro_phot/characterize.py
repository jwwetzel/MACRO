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
#: 1%, not the 0.1% the strategy document assumes, because the threshold is
#: an empirical quantile of a few hundred signal-free trials and a 0.1%
#: quantile of 300 samples is not a measurement.  A stricter FAP would move
#: the recovery contours by well under the amplitude grid spacing.
DETECT_FAP = 0.01

#: A recovered period counts as CORRECT when it lands within this
#: fractional distance of the injected one.  1% is far tighter than the
#: 1 c/d alias spacing at these periods, so "correct" really does mean the
#: right peak and not its neighbour.
#:
#: READ THE NEXT PARAGRAPH BEFORE QUOTING ANY CONTOUR SCORED WITH THIS.
#: A 1% window is 0.12-0.14 c/d at these periods.  On a SINGLE NIGHT the
#: frequency resolution is 1/T = 2.6-9.0 c/d, so the acceptance window is 20
#: to 70 times narrower than the peak itself and no amplitude, however
#: large, can reliably place the maximum inside it: a 300 mmag injection
#: into VV Pup's richest night exceeds the detection threshold in 40 trials
#: out of 40 while landing in the window in only 25.  A contour scored this
#: way is a PERIOD-DETERMINATION limit, not a detection limit, and the two
#: differ by 3-8x.  :data:`SCORE_MODES` names both and the report must say
#: which one it is quoting.
PERIOD_TOL_FRAC = 0.01

#: The two questions an injection-recovery run can answer, kept apart
#: because this campaign conflated them and three verdicts moved when they
#: were separated:
#:
#: * ``'period'`` — can the sampling MEASURE the period from scratch?  The
#:   highest peak in the search band must clear a threshold set on the same
#:   max statistic AND land within :data:`PERIOD_TOL_FRAC` of the truth.
#:   This is the right question for a blind survey search.
#: * ``'known'`` — can a modulation be DETECTED at a frequency already known
#:   from the literature?  The power at that one frequency must clear a
#:   threshold set on the power at that same frequency in signal-free data.
#:   This is the right question for these five CVs, whose orbital periods
#:   have decades-long published ephemerides, and it is the question the
#:   paper actually asks.
SCORE_MODES = ("period", "known")

#: Recovery contour level.
RECOVERY_LEVEL = 0.90

#: Half-width in magnitudes of the magnitude-matched comparison used to test
#: whether the held-out check stars are a FAIR sample of what a star at the
#: target's brightness achieves, or merely the quiet survivors of the
#: comparison-star stability cut.  +/-0.25 mag is narrow enough that photon
#: noise is nearly constant across the sample and wide enough to hold tens of
#: stars in every field here.
FIELD_MATCH_HALF_WIDTH = 0.25


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


def binned_median_pooled_tail(x, y, edges,
                              min_count: int = MIN_FRAMES_PER_QUALITY_BIN
                              ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """:func:`binned_median`, with the sparse TAIL merged into one wide bin.

    A count rule that silently drops bins converts "too few frames to test"
    into "no effect", and this archive published exactly that mistake:
    airmass was reported to degrade the check stars nowhere "over the range
    actually observed", while the bins from X = 2.45 to 2.65 sit at 1.48,
    1.79 and 1.45 times baseline — all three above the degradation factor,
    all three excluded for holding fewer than 15 frames each.  Pooled they
    hold 27 frames and the response is testable.

    Everything from the first under-populated bin onward is merged into a
    single bin whose median is recomputed over the POOLED SAMPLE (a median
    of medians would be a different and wrong statistic).  Bins before that
    are untouched, so the well-sampled part of the axis reads exactly as
    before.

    The pooled bin's reported x is its LOWER EDGE, not the mean x of its
    frames.  That is deliberate: this bin exists to be a threshold
    candidate, and the only thing it can honestly support is "beyond here
    the data no longer demonstrate their own quality".  Placing it at the
    mean would license everything up to the mean on evidence that covers
    the whole tail.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    edges = np.asarray(edges, dtype=float)
    centers, med, cnt = binned_median(x, y, edges)
    thin = np.flatnonzero((cnt > 0) & (cnt < int(min_count)))
    if thin.size == 0:
        return centers, med, cnt
    # The tail starts at the first sparse bin that has only sparse bins
    # after it — an isolated thin bin in the middle of a well-sampled axis
    # is not a tail and must not swallow everything beyond it.
    start = centers.size
    for i in range(centers.size - 1, -1, -1):
        if cnt[i] == 0 or cnt[i] < int(min_count):
            start = i
        elif cnt[i] >= int(min_count):
            break
    if start >= centers.size:
        return centers, med, cnt
    ok = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[ok], y[ok]
    sel = xs >= edges[start]
    if not sel.any():
        return centers[:start], med[:start], cnt[:start]
    pooled_center = float(edges[start])
    return (np.concatenate([centers[:start], [pooled_center]]),
            np.concatenate([med[:start], [float(np.median(ys[sel]))]]),
            np.concatenate([cnt[:start], [int(sel.sum())]]))


def degradation_threshold(centers, medians, counts,
                          factor: float = DEGRADE_FACTOR,
                          min_count: int = MIN_FRAMES_PER_QUALITY_BIN,
                          baseline: Optional[float] = None,
                          run_length: int = 2,
                          allow_terminal_single: bool = True
                          ) -> tuple[float, float]:
    """Where a quality axis starts costing precision, and the baseline.

    ``baseline`` is what the data achieve when the axis is favourable.
    Supply it (the pooled median of a normalised scatter is the natural
    choice, and equals 1 by construction) or leave it ``None`` to take the
    smallest well-populated bin median.

    The threshold is the centre of the first bin that STARTS a run of
    ``run_length`` consecutive well-populated bins all exceeding
    ``factor * baseline``.  Requiring a run is not fussiness: a median over
    ~20 frames fluctuates by tens of percent, and a single-bin excursion
    would otherwise set an archive-wide cut.  (It did, in the first version
    of this build: a lone bin at 1.4x threw away 69% of the frames.)

    ``allow_terminal_single`` lets the LAST well-populated bin set the
    threshold on its own.  A run of two protects against an interior
    excursion, and an interior excursion by definition has a neighbour to
    confirm it against; the final bin has none, so demanding a run there
    silently converts "the axis ends here" into "no effect".  That is
    exactly what happened to airmass: pooled above X = 2.4 the check stars
    sit at 1.33x baseline over 36 frames, and the page reported that
    airmass "never" degrades them over the observed range.

    Returns ``(threshold_x, baseline_y)``; the threshold is ``inf`` when no
    such run exists — itself a result ("this axis does not limit us"), not
    a failure.

    This is the whole method for defending a cut: the number comes from the
    data's own response, so the only judgement calls left in the open are
    ``factor`` and ``run_length``, which the report prints.
    """
    centers = np.asarray(centers, dtype=float)
    medians = np.asarray(medians, dtype=float)
    counts = np.asarray(counts, dtype=int)
    good = (counts >= min_count) & np.isfinite(medians)
    if not good.any():
        return float("inf"), float("nan") if baseline is None else float(baseline)
    if baseline is None:
        base_i = int(np.flatnonzero(good)[np.argmin(medians[good])])
        baseline = float(medians[base_i])
        start = base_i + 1
    else:
        baseline = float(baseline)
        start = 0
    over = good & (medians > factor * baseline)
    for i in range(start, centers.size - run_length + 1):
        if all(over[i:i + run_length]):
            return float(centers[i]), baseline
    if allow_terminal_single:
        last = np.flatnonzero(good)
        if last.size and over[last[-1]] and last[-1] >= start:
            return float(centers[last[-1]]), baseline
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

    The fit is RELATIVE (every row divided by that star's ``rms^2``), which
    matters more than it looks: an unweighted fit in variance is dominated
    by the faintest stars, whose variance is 100x the bright end's, and it
    returns a "floor" that is really the faint end's photon noise.  That
    error is not hypothetical — the first version of this build reported
    30-60 mmag floors for series whose bright stars sit at 8 mmag.

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
    w = 1.0 / r ** 2
    A = np.column_stack([p ** 2 * w, w])
    coef, *_ = np.linalg.lstsq(A, r ** 2 * w, rcond=None)
    k = max(float(coef[0]), 0.0)
    floor = math.sqrt(max(float(coef[1]), 0.0))
    return floor, k, int(r.size)


def noise_plateau(mag, rms, bin_width: float = 0.5, min_stars: int = 5
                  ) -> tuple[float, float, int]:
    """Model-free version of "where the curve goes flat".

    Bins the constant-star cloud in magnitude, keeps bins with at least
    ``min_stars``, and returns the SMALLEST bin median together with the
    magnitude where it occurs and the number of stars in that bin.  This is
    the floor as a reader sees it on the RMS-vs-magnitude figure: no model,
    no gain assumption, no fit — the best scatter this series achieves on
    any star of any brightness.  Quoted beside the fitted floor so a reader
    can see whether the model and the eye agree.
    """
    m = np.asarray(mag, dtype=float)
    r = np.asarray(rms, dtype=float)
    ok = np.isfinite(m) & np.isfinite(r) & (r > 0)
    m, r = m[ok], r[ok]
    if m.size == 0:
        return float("nan"), float("nan"), 0
    lo = math.floor(m.min() / bin_width) * bin_width
    idx = np.floor((m - lo) / bin_width).astype(int)
    best = (float("inf"), float("nan"), 0)
    for b in np.unique(idx):
        sel = idx == b
        if sel.sum() < min_stars:
            continue
        med = float(np.median(r[sel]))
        if med < best[0]:
            best = (med, lo + (b + 0.5) * bin_width, int(sel.sum()))
    if not math.isfinite(best[0]):
        return float("nan"), float("nan"), 0
    return best


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


def red_noise_factor(tau_s, adev, tau_target_s: float
                     ) -> tuple[float, float]:
    """How much worse than white the noise is, and AT WHICH timescale.

    Extrapolates the FIRST rung (tau_0, sigma_0) forward as white noise,
    ``sigma_white = sigma_0 * sqrt(tau_0 / tau_target)``, and divides the
    measured Allan deviation at (or nearest below) ``tau_target`` by it.
    A factor of 1 means averaging works perfectly; 2 means a mean over that
    timescale is twice as uncertain as the per-point scatter implies.

    Returns ``(factor, tau_used_s)`` — BOTH, always.  It used to return the
    factor alone, and the caller stored it in a column named
    ``red_factor_porb`` on the assumption that the ladder reached the
    orbital period.  Only 11 of 92 ladders do; the other 81 top out at
    0.15-0.51 P_orb, and because red noise GROWS with tau, every one of
    those was a lower bound published under the name of a larger quantity.
    Returning the tau makes that impossible to store silently.
    """
    t = np.asarray(tau_s, dtype=float)
    a = np.asarray(adev, dtype=float)
    ok = np.isfinite(t) & np.isfinite(a) & (t > 0) & (a > 0)
    t, a = t[ok], a[ok]
    if t.size < 2:
        return float("nan"), float("nan")
    order = np.argsort(t)
    t, a = t[order], a[order]
    at_or_below = np.flatnonzero(t <= tau_target_s)
    if at_or_below.size == 0:
        return float("nan"), float("nan")
    i = int(at_or_below[-1])
    white = a[0] * math.sqrt(t[0] / t[i])
    if white <= 0:
        return float("nan"), float(t[i])
    return float(a[i] / white), float(t[i])


def white_noise_allan_null(n_pts: int, dt_s: float, tau_target_s: float,
                           allan_fn, n_real: int = 200,
                           rng: Optional[np.random.Generator] = None
                           ) -> dict:
    """What a LADDER THIS SHORT returns when the noise really is white.

    The campaign compared a median fitted Allan slope of -0.38 against the
    theoretical white value of -0.50 and called the difference evidence of
    correlated noise.  That comparison is not available: these ladders have
    4-6 rungs over about one decade with as few as 3 pairs at the top rung,
    and the ESTIMATOR's expectation on white noise is not -0.50 — it is
    -0.55 with a 5-95% range of [-0.88, -0.30].  A median of -0.38 sits
    inside that range, so the median comparison proves nothing, even though
    a per-ladder comparison against this null does show real red noise in a
    third of them.

    Generates ``n_real`` pure-Gaussian series of ``n_pts`` points at spacing
    ``dt_s``, pushes each through the caller's ``allan_fn(y, dt) ->
    (tau, adev, n_pairs)`` (``macro_phot.errors.allan_deviation``; passed in
    rather than imported so this module keeps its no-heavy-imports rule),
    and returns the percentiles of the slope and of the red-noise factor
    that white noise alone produces.

    Returns ``{'slope_p50', 'slope_p05', 'slope_p95', 'red_p50', 'red_p95',
    'n'}``; NaNs when no realization yields a usable ladder.
    """
    rng = rng or np.random.default_rng(20260819)
    slopes, reds = [], []
    for _ in range(int(n_real)):
        y = rng.normal(0.0, 1.0, int(n_pts))
        tau, adev, _n = allan_fn(y, float(dt_s))
        if np.asarray(tau).size < 3:
            continue
        s = allan_slope(tau, adev)
        r, _tu = red_noise_factor(tau, adev, tau_target_s)
        if np.isfinite(s):
            slopes.append(s)
        if np.isfinite(r):
            reds.append(r)
    if not slopes:
        return {k: float("nan") for k in
                ("slope_p50", "slope_p05", "slope_p95", "red_p50", "red_p95")
                } | {"n": 0}
    return {"slope_p50": float(np.percentile(slopes, 50)),
            "slope_p05": float(np.percentile(slopes, 5)),
            "slope_p95": float(np.percentile(slopes, 95)),
            "red_p50": float(np.percentile(reds, 50)) if reds else float("nan"),
            "red_p95": float(np.percentile(reds, 95)) if reds else float("nan"),
            "n": len(slopes)}


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


def ls_power_at(times_d, y, freq_cd: float) -> float:
    """Lomb-Scargle power at ONE stated frequency.

    The companion to :func:`ls_peak` for the ``'known'`` score mode: when
    the period is already known from decades of published ephemerides, the
    statistic that decides a detection is the power AT that frequency, not
    the height of whatever peak happens to be tallest in the band.  Nothing
    else differs — same floating-mean periodogram, same normalisation — so
    the two statistics are directly comparable against their own nulls.
    """
    from astropy.timeseries import LombScargle
    t = np.asarray(times_d, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 5:
        return float("nan")
    return float(LombScargle(t, y).power(np.array([float(freq_cd)]))[0])


def detection_threshold(times_d, residual_pool, freqs_cd,
                        n_trials: int = 200,
                        fap: float = DETECT_FAP,
                        rng: Optional[np.random.Generator] = None,
                        at_freq_cd: Optional[float] = None) -> float:
    """Bootstrap detection threshold from SIGNAL-FREE data.

    ``residual_pool`` is a list of real check-star residual series measured
    at the same timestamps.  Realizations are cyclic rolls of those series,
    so the threshold inherits the archive's real correlated noise; an
    analytic Baluev FAP computed on this data would be optimistic by the
    same factor the noise is red.

    Two nulls, matching the two questions in :data:`SCORE_MODES`, because a
    test statistic must be compared against ITS OWN null distribution:

    * ``at_freq_cd=None`` — the ``1 - fap`` quantile of the HIGHEST peak
      anywhere in ``freqs_cd``.  The null for a blind search, and it carries
      the whole band's look-elsewhere penalty.
    * ``at_freq_cd=f`` — the ``1 - fap`` quantile of the power at that one
      frequency.  The null for a search at a known period.  It is lower
      than the max-statistic threshold by exactly the look-elsewhere factor
      the known period buys back, which is why using the max threshold for
      a known-period claim understates the data by a large factor.
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
        p = (ls_power_at(times_d, y, at_freq_cd) if at_freq_cd is not None
             else ls_peak(times_d, y, freqs_cd)[1])
        if np.isfinite(p):
            peaks.append(p)
    if not peaks:
        return float("nan")
    return float(np.quantile(peaks, 1.0 - fap))


def threshold_spread(times_d, residual_pool, freqs_cd,
                     n_trials: int = 100,
                     fap: float = DETECT_FAP,
                     rng: Optional[np.random.Generator] = None,
                     at_freq_cd: Optional[float] = None) -> dict:
    """The detection threshold computed from EACH check star separately.

    The pooled threshold is a quantile of a maximum statistic drawn from a
    pool of at most four held-out stars, so it is set by the tail of the
    worst of them rather than by a population.  In ``stlmi|e76|g`` the four
    per-star 99th percentiles are 0.053, 0.062, 0.085 and 0.182: the pooled
    value is essentially the last star alone, and the between-star spread is
    30% of the threshold.  A bootstrap over trials reports +/-1% on that
    number and hides the fact entirely.

    Returns ``{'per_star': [...], 'spread_frac': (max-min)/median, 'n': N}``
    so the report can print the spread beside the threshold instead of a
    precision the pool cannot support.
    """
    rng = rng or np.random.default_rng(20260819)
    per = []
    for r in residual_pool:
        r = np.asarray(r, dtype=float)
        if np.isfinite(r).sum() < 5:
            continue
        thr = detection_threshold(times_d, [r], freqs_cd, n_trials, fap,
                                  rng, at_freq_cd)
        if np.isfinite(thr):
            per.append(float(thr))
    if not per:
        return {"per_star": [], "spread_frac": float("nan"), "n": 0}
    med = float(np.median(per))
    return {"per_star": per,
            "spread_frac": (max(per) - min(per)) / med if med > 0 else float("nan"),
            "n": len(per)}


def classify_peak(f_peak: float, f_true: float, f_alias_cd: float = 1.0,
                  tol_frac: float = PERIOD_TOL_FRAC,
                  max_order: int = 3) -> str:
    """Name what a periodogram peak actually is: truth, alias, or neither.

    Returns ``'true'``, ``'alias'`` or ``'other'``.  The alias test walks
    ``f_true +/- k * f_alias`` for ``k`` up to ``max_order`` and folds
    negative frequencies back (a periodogram cannot tell +f from -f).

    This exists because the campaign converted a WINDOW statistic — the
    +/-1 c/d aliases carry up to 0.92 of the window power — straight into
    the operational claim that "a multi-night periodogram cannot on its own
    choose between the true frequency and its daily aliases".  Window power
    is necessary for that conclusion and nowhere near sufficient: whether
    the wrong peak actually wins depends on amplitude and on the number of
    nights, and both are measurable with the injection machinery already
    built.  This function is what turns each trial into a vote.
    """
    tol = tol_frac * float(f_true)
    if abs(float(f_peak) - float(f_true)) <= tol:
        return "true"
    for k in range(-max_order, max_order + 1):
        if k == 0:
            continue
        f = abs(float(f_true) + k * float(f_alias_cd))
        if f > 1e-9 and abs(float(f_peak) - f) <= tol:
            return "alias"
    return "other"


def recovery_fraction(times_d, residual_pool, freqs_cd, period_d: float,
                      semi_amp_mag: float, threshold: float,
                      n_trials: int = 60,
                      period_tol: float = PERIOD_TOL_FRAC,
                      rng: Optional[np.random.Generator] = None,
                      score: str = "period",
                      detrend_nights=None) -> float:
    """Fraction of injections recovered, under a STATED score mode.

    ``score='period'`` (the blind-search question): the highest peak must
    exceed ``threshold`` and land within ``period_tol`` of the injected
    frequency.  An alias-family peak above threshold is a detection of the
    wrong period and is scored as a failure — which is correct for that
    question, and disastrous when the answer is read as a detection limit
    (see the warning on :data:`PERIOD_TOL_FRAC`).

    ``score='known'`` (the question this paper asks): the power AT the
    injected frequency must exceed ``threshold``, which the caller must
    have computed with ``detection_threshold(..., at_freq_cd=1/period_d)``.
    No frequency tolerance enters, because no frequency is being measured.

    ``detrend_nights``, when given, applies per-night mean removal to the
    injected series exactly as a real analysis would — signal and noise
    together, so the contour pays for the signal power detrending destroys
    as well as banking the red noise it removes.
    """
    if score not in SCORE_MODES:
        raise ValueError(f"score must be one of {SCORE_MODES}, got {score!r}")
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
        y = noise + sig
        if detrend_nights is not None:
            y = remove_nightly_means(y, detrend_nights)
        if score == "known":
            p = ls_power_at(times_d, y, f_true)
            if np.isfinite(p) and p > threshold:
                hits += 1
        else:
            f, p = ls_peak(times_d, y, freqs_cd)
            if (np.isfinite(p) and p > threshold
                    and abs(f - f_true) <= period_tol * f_true):
                hits += 1
    return hits / float(n_trials)


def remove_nightly_means(y, nights) -> np.ndarray:
    """Subtract each night's own mean — the detrending a real analysis does.

    The strategy document requires per-night systematics to be fitted
    jointly with the periodic model (§4.20).  The cheapest honest stand-in
    is removing the nightly mean, and it matters enormously here: night-to-
    night zero-point wander is low-frequency power that the 1 c/d sampling
    comb ALIASES straight into the 2-40 c/d search band, where it raises the
    false-alarm threshold for every frequency at once.  Removing it is not
    free either — it also removes signal power on timescales approaching a
    night — which is exactly why both the raw and the detrended contour are
    computed and published side by side.
    """
    y = np.array(y, dtype=float)
    nights = np.asarray(nights)
    for n in np.unique(nights):
        sel = nights == n
        if sel.sum() >= 2:
            y[sel] -= np.mean(y[sel])
    return y


def alias_confusion(times_d, residual_pool, freqs_cd, period_d: float,
                    semi_amp_mag: float, n_trials: int = 40,
                    f_alias_cd: float = 1.0,
                    rng: Optional[np.random.Generator] = None,
                    detrend_nights=None) -> dict:
    """How often the tallest peak is an ALIAS rather than the true frequency.

    The measurement S3's verdict should have rested on and never made.  A
    window statistic says how much power the sampling makes available at
    ``f +/- k c/d``; it cannot say whether the wrong one wins, because that
    depends on the signal-to-noise the real modulation brings.  Here it is
    measured directly, through the same injections and the same real
    correlated noise: inject at the known period, take the tallest peak in
    the band, and classify it with :func:`classify_peak`.

    Returns ``{'true': f, 'alias': f, 'other': f, 'n_trials': N}`` as
    fractions.  On these data the answer at the real modulation amplitudes
    (0.5-2 mag, one to two orders of magnitude above the season contour) is
    ``true`` in every trial — so the honest S3 statement is amplitude
    dependent, not "cannot choose".
    """
    rng = rng or np.random.default_rng(20260819)
    pool = [np.asarray(r, dtype=float) for r in residual_pool
            if np.isfinite(r).sum() >= 5]
    counts = {"true": 0, "alias": 0, "other": 0}
    if not pool:
        return dict(counts, n_trials=0)
    f_true = 1.0 / period_d
    for _ in range(n_trials):
        base = pool[rng.integers(len(pool))]
        noise = cyclic_noise_realization(base, rng.integers(base.size))
        sig = inject_sinusoid(times_d, period_d, semi_amp_mag,
                              phase=float(rng.random()))
        y = noise + sig
        if detrend_nights is not None:
            y = remove_nightly_means(y, detrend_nights)
        f, _p = ls_peak(times_d, y, freqs_cd)
        if not np.isfinite(f):
            continue
        counts[classify_peak(f, f_true, f_alias_cd)] += 1
    n = max(sum(counts.values()), 1)
    return {k: v / float(n) for k, v in counts.items()} | {"n_trials": n}


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
                        rng: Optional[np.random.Generator] = None,
                        noise_pool: Optional[Sequence[np.ndarray]] = None,
                        fit_ingress_phase: Optional[float] = None,
                        fit_depth_mag: Optional[float] = None) -> float:
    """Monte-Carlo epoch uncertainty of ONE eclipse/bright-phase edge, seconds.

    Injects a trapezoidal feature at a random epoch into noise at the
    supplied REAL timestamps of one cycle, then
    recovers the epoch by chi-square minimisation over a fine epoch grid
    with the template shape held fixed (the optimistic case: shape known
    exactly).  The returned scatter is the 68.3rd percentile of the
    absolute epoch error — the one-sigma half-width, which is what an error
    bar means — and it is a LOWER BOUND on what a real fit achieves.

    Two details that matter, both consequences of coarse sampling.  (a) The
    chi-square is PIECEWISE CONSTANT in the trial epoch: shifting the
    template by less than one sampling interval usually changes no point's
    membership of the eclipse, so whole plateaux of trial epochs tie.  The
    estimator takes the CENTRE of the winning plateau, because ``argmin``
    would take its left edge and manufacture a systematic bias.  (b) For the
    same reason a percentile, not a MAD, measures the scatter: when more
    than half the trials land exactly right the MAD is identically zero and
    would report perfect timing from a sampling-limited measurement.

    The reason this is the decisive measurement for bright-phase timing:
    when the sampling interval is long compared with the ingress, the
    answer is set by the sampling, not by the photometry, and no amount of
    signal-to-noise improves it.  A pipeline that assumes otherwise reports
    a timing precision it cannot deliver.

    ``noise_pool``, when supplied, replaces the Gaussian draw with cyclic
    rolls of REAL check-star residual series measured at these same
    timestamps.  That matters: the Allan analysis shows this archive's noise
    is not white on the orbital timescale, and correlated noise on the
    timescale of the feature being timed is precisely the noise that biases
    an epoch.  Gaussian noise is the optimistic case and is kept only for
    the unit tests, where the right answer is known in closed form.

    ``fit_ingress_phase`` / ``fit_depth_mag`` let the RECOVERY template
    differ from the injected one.  Left at ``None`` (template exactly right)
    the answer is the best case any analysis could ever reach; set to a
    different edge sharpness or amplitude they measure what a wrong assumed
    shape costs.  That is the realistic case here, because the polar ingress
    duration is colour-dependent and this archive has not measured it — so
    the pair of answers brackets the truth instead of asserting one end.
    """
    rng = rng or np.random.default_rng(20260819)
    t = np.asarray(times_d, dtype=float)
    t = t[np.isfinite(t)]
    if t.size < 6 or not np.isfinite(sigma_mag) or sigma_mag <= 0:
        return float("nan")
    span = t.max() - t.min()
    ingress_d = ingress_phase * period_d
    dt_d = float(np.median(np.diff(np.sort(t)))) if t.size > 1 else ingress_d

    fit_ing = ingress_phase if fit_ingress_phase is None else float(fit_ingress_phase)
    fit_depth = depth_mag if fit_depth_mag is None else float(fit_depth_mag)

    def _best_offset(y, t0, grid):
        """Offset minimising chi-square on one grid; plateau CENTRE."""
        # Vectorized template bank: (n_grid, n_points).  The FIT template's
        # shape may deliberately differ from the injected one.
        bank = np.stack([eclipse_template(t, period_d, fit_depth, width_phase,
                                          t0 + g, fit_ing) for g in grid])
        chi2 = np.sum((y[None, :] - bank) ** 2, axis=1)
        tie = chi2 <= chi2.min() + 1e-12 * max(float(chi2.max()), 1.0)
        return float(np.median(grid[tie]))

    # Two-stage search.  Stage 1 sweeps the whole cycle at a quarter of the
    # smaller of (sampling interval, ingress duration) — fine enough that no
    # plateau is missed, coarse enough to be cheap.  Stage 2 refines around
    # the winner 50x finer, so the GRID never sets the answer: without it a
    # densely-sampled, high-S/N case returns exactly zero and would be
    # reported as perfect timing.
    step1 = max(min(dt_d, ingress_d, fit_ing * period_d) / 4.0,
                period_d / 20000.0)
    grid1 = np.arange(-0.5 * period_d, 0.5 * period_d + step1, step1)
    # Stage-2 step: fine enough to resolve the PHOTON-limited edge precision
    # (an edge of duration ``ingress_d`` measured at relative depth
    # sigma/depth locates to ~ingress_d * sigma/depth), capped so the refine
    # grid never exceeds a few thousand points.
    photon_scale = ingress_d * float(sigma_mag) / max(float(depth_mag), 1e-9)
    step2 = max(min(step1 / 50.0, photon_scale / 5.0), 4.0 * step1 / 4000.0)
    pool = [np.asarray(p, dtype=float) for p in (noise_pool or [])
            if np.asarray(p).size == t.size]
    errs = []
    for _ in range(n_trials):
        t0 = t.min() + span * float(rng.random())
        if pool:
            base = pool[rng.integers(len(pool))]
            noise = cyclic_noise_realization(base, rng.integers(base.size))
        else:
            noise = rng.normal(0, sigma_mag, t.size)
        y = eclipse_template(t, period_d, depth_mag, width_phase, t0,
                             ingress_phase) + noise
        g1 = _best_offset(y, t0, grid1)
        grid2 = np.arange(g1 - 2 * step1, g1 + 2 * step1 + step2, step2)
        errs.append(_best_offset(y, t0, grid2))
    e = np.abs(np.asarray(errs))
    return float(np.percentile(e, 68.27)) * 86400.0


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


# ==========================================================================
# 5 - TURNING MEASUREMENTS INTO VERDICTS
#
# Each function here exists because a verdict in the first version of this
# build was decided by a number that did not measure the thing the goal
# asked about.  Keeping them pure and named means the mismatch is visible
# in the verdict table rather than buried in an f-string.
# ==========================================================================

def colour_point_sigma(sigma_a_mag: float, sigma_b_mag: float,
                       rate_mag_per_s: float, dt_s: float) -> float:
    """Uncertainty of ONE colour point from two non-simultaneous bands.

        sigma_colour^2 = sigma_a^2 + sigma_b^2 + (dm/dt * dt)^2

    The first two terms are the per-point precisions the report already
    measures.  The third is the one that was missing, and on these targets
    it dominates: Q1 asks for "quasi-simultaneous" three-colour coverage
    and was graded on the SINGLE-BAND precision, 14-44 mmag.  But g and i
    are taken minutes apart on a star whose brightness sweeps 1-2 mag per
    orbit.  On ST LMi's richest three-filter night the g-to-nearest-i offset
    has a median of 75 s and a 90th percentile of 144 s, while the target's
    own rate of change has a median of 0.30 and a p90 of 1.01 mmag/s — so
    the phase-offset term alone is a median 22 mmag and a p90 of 146 mmag,
    several times the photometric term on the bright-phase edges.

    A colour is never better than sqrt(2) times a single band even when
    simultaneous, which this formula also enforces.
    """
    a = float(sigma_a_mag)
    b = float(sigma_b_mag)
    drift = float(rate_mag_per_s) * float(dt_s)
    if not all(math.isfinite(v) for v in (a, b, drift)):
        return float("nan")
    return math.sqrt(a * a + b * b + drift * drift)


def duty_cycle_sigma(n_epochs: int, p: float = 0.5) -> float:
    """Binomial 1-sigma on a duty cycle measured from ``n_epochs`` nights.

        sigma = sqrt(p (1 - p) / N)

    THE number Q4 should always have been graded on.  A duty cycle is a
    fraction of TIME SPENT in a state; its uncertainty is set by how many
    independent epochs sampled it, not by the per-point photometric
    precision.  Grading it on precision-versus-state-separation (9-77 mmag
    against 1-3 mag) only establishes that a single night can be
    classified, which nobody disputed.  With 13-37 nights per target the
    answer is +/-8 to +/-14 percentage points, and that is before any
    seasonal observability bias in WHEN those nights fall.

    ``p = 0.5`` is the worst case and therefore the honest default.
    """
    n = int(n_epochs)
    if n <= 0:
        return float("nan")
    return math.sqrt(float(p) * (1.0 - float(p)) / n)


def contour_uncertainty(amps, fracs, n_trials: int,
                        level: float = RECOVERY_LEVEL,
                        n_boot: int = 200,
                        rng: Optional[np.random.Generator] = None
                        ) -> tuple[float, float]:
    """A 90% recovery contour's own error bar, by binomial resampling.

    Each grid cell's recovery fraction is ``k`` hits out of ``n_trials``,
    and at a fraction of 0.9 with 50 trials the binomial standard error is
    4.2% — while :func:`recovery_contour` interpolates between two such
    cells and the report printed the answer to 0.1 mmag.  Resampling every
    cell from ``Binomial(n_trials, frac)`` and re-deriving the contour gives
    the 16th-84th percentile band the printed value has to be rounded to.

    Returns ``(lo, hi)``; NaNs when the grid never reaches ``level``.
    """
    rng = rng or np.random.default_rng(20260819)
    a = np.asarray(amps, dtype=float)
    f = np.asarray(fracs, dtype=float)
    ok = np.isfinite(a) & np.isfinite(f)
    a, f = a[ok], f[ok]
    if a.size == 0 or n_trials <= 0:
        return float("nan"), float("nan")
    draws = []
    for _ in range(int(n_boot)):
        fb = rng.binomial(int(n_trials), np.clip(f, 0.0, 1.0)) / float(n_trials)
        c = recovery_contour(a, fb, level)
        if np.isfinite(c):
            draws.append(c)
    if len(draws) < 3:
        return float("nan"), float("nan")
    return (float(np.percentile(draws, 16)),
            float(np.percentile(draws, 84)))


def rate_of_change_mag_per_s(times_d, mags) -> np.ndarray:
    """Per-point |dm/dt| in mag/s from a target's own light curve.

    Symmetric difference where a neighbour exists on each side, one-sided
    at the ends.  Used only to price the non-simultaneity term of a colour
    point (:func:`colour_point_sigma`), so it deliberately measures the
    star's ACTUAL sweep rate on the night rather than a model of it.
    """
    t = np.asarray(times_d, dtype=float)
    m = np.asarray(mags, dtype=float)
    ok = np.isfinite(t) & np.isfinite(m)
    t, m = t[ok], m[ok]
    if t.size < 2:
        return np.array([])
    order = np.argsort(t)
    t, m = t[order], m[order]
    return np.abs(np.gradient(m, t * 86400.0))


def nearest_time_offsets(times_a_d, times_b_d) -> np.ndarray:
    """For each time in A, the gap in SECONDS to the nearest time in B.

    The measured cost of interleaved filters: this is the ``dt`` that
    :func:`colour_point_sigma` multiplies by the target's sweep rate.
    Returns an empty array when either side is empty.
    """
    a = np.sort(np.asarray(times_a_d, dtype=float))
    b = np.sort(np.asarray(times_b_d, dtype=float))
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size == 0 or b.size == 0:
        return np.array([])
    idx = np.searchsorted(b, a)
    lo = np.clip(idx - 1, 0, b.size - 1)
    hi = np.clip(idx, 0, b.size - 1)
    d = np.minimum(np.abs(a - b[lo]), np.abs(a - b[hi]))
    return d * 86400.0
