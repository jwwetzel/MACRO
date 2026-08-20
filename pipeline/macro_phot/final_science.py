"""CV-S10 — the arithmetic behind the two closing science decisions.

Two questions live here, and neither is a measurement of the sky so much as
a decision about what the 2024-2025 archive is allowed to claim.

**YZ Cnc, on the FALLBACK branch.**  CV-S7 measured that no dense run sits
inside a superoutburst, so the superhump branch is closed and the strategy's
fallback -- "orbital hump + flickering statistics, honest but weaker" -- is
what the season supports.  Delivering that fallback properly means three
separate pieces of arithmetic:

* fold each QUIESCENT dense run on the published orbital period and measure
  the hump's semi-amplitude and phase (:func:`fold_fit`), then decide
  whether it is DETECTED against an injection-recovery contour measured on
  that run's own timestamps and that run's own noise -- never by eye;
* separate FLICKERING from photometric noise with a structure function
  (:func:`structure_function`) whose floor is measured, not assumed: the
  same statistic computed on magnitude-matched field stars observed through
  the same frames (:func:`quadrature_excess`);
* characterise the NORMAL-OUTBURST runs on their own terms
  (:func:`linear_rate_per_hour`, :func:`percentile_amplitude`), because six
  dense multi-colour runs inside normal outbursts are a data set, not a
  consolation prize.

**AN UMa, per filter.**  CV-S5 already graded the three-filter colour goal
NOT SUPPORTED on 4 qualifying nights against a >= 8 bar.  That is one goal.
This module asks the narrower question the plan actually poses -- what can
EACH filter support on its own -- and answers it with
:func:`capability_verdict`, which is deliberately dumb: it compares one
measured number against one stated bar and returns a word.  Every judgement
call therefore lives in the CALLER's choice of bar, where it can be argued
with, rather than inside a scoring function nobody can audit.

Everything in this file is a pure function of arrays and floats.  The I/O,
the SQL and the run selection live in ``pipeline/scripts/run_cv_final.py``;
the page lives in :mod:`macro_phot.report_final`.  The tests in
``pipeline/tests/test_final_science.py`` inject known humps, known
flickering and known noise floors and demand them back.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np


# ==========================================================================
# 0.  CONSTANTS — every threshold this stage applies, in one place
# ==========================================================================

#: Harmonics fitted to a folded run.  TWO, not one, and the second is the
#: point: a dwarf nova's quiescent orbital hump is not a sinusoid -- it is a
#: one-sided bump from the bright spot where the stream hits the disc rim --
#: so a pure fundamental systematically under-reads its height.  The
#: SEMI-AMPLITUDE this module quotes is still the FUNDAMENTAL's, because
#: that is the quantity an injection-recovery contour is calibrated in
#: (:func:`macro_phot.characterize.inject_sinusoid` injects a sinusoid), and
#: comparing a peak-to-peak morphology against a sinusoidal contour would be
#: comparing two different numbers.  The harmonic is fitted so it cannot
#: leak into the fundamental, and its own amplitude is reported beside it as
#: the honest statement of how non-sinusoidal the hump is.
N_HARMONICS = 2

#: Half-width, in magnitudes, of the field-star sample that measures the
#: photometric noise floor a flickering amplitude is subtracted against.
#: Same value as :data:`macro_phot.characterize.FIELD_MATCH_HALF_WIDTH`, and
#: for the same reason: wide enough to hold tens of stars in these fields,
#: narrow enough that photon noise barely varies across the sample.
#:
#: This is NOT a cosmetic choice.  YZ Cnc's four held-out CHECK stars sit at
#: cal_mag 13.10-13.26 while the star at quiescence sits at 14.1-14.6, so
#: the check stars are about a magnitude BRIGHTER than the target and their
#: scatter is an optimistic floor.  Matching on brightness is what makes the
#: subtraction defensible.
FLOOR_MATCH_HALF_WIDTH = 0.25

#: A field star must cover this fraction of a run's frames before its
#: structure function may enter the floor.  A star seen in a third of the
#: frames has a structure function built from a different sampling than the
#: target's, and the floor is supposed to differ from the target in exactly
#: one respect: the absence of flickering.
FLOOR_MIN_COVERAGE = 0.50

#: Minimum magnitude-matched stars before a floor is quoted at all.  Below
#: this the median over stars is not a measurement and the excess is
#: reported as UNMEASURED rather than as zero.
FLOOR_MIN_STARS = 5

#: Minimum pairs in a structure-function bin before the bin is used.  A
#: two-pair bin has a 50% error bar on its own variance.
SF_MIN_PAIRS = 20

#: How many sigma of variance excess before flickering is called DETECTED in
#: a timescale bin.  3, on the difference of two variances, each with its
#: own sampling error -- see :func:`excess_significance`.
FLICKER_SIGMA_BAR = 3.0

#: A superhump semi-amplitude, in magnitudes, at the low end of what the
#: literature reports for SU UMa stars (common superhumps run ~0.05-0.15 mag
#: peak-to-peak in the optical).  Used ONLY to state what an outburst run's
#: blind-period contour would have had to reach to measure a superhump
#: period; never as a detection or a claim about YZ Cnc.
SUPERHUMP_SEMI_AMP_FLOOR = 0.050

#: Accumulated phase uncertainty, in cycles, above which two runs may not be
#: put on a common absolute phase axis.  0.1 cycle is the strategy
#: document's own bar for re-deriving a local ephemeris (chair's call 4).
PHASE_DRIFT_BAR_CYCLES = 0.10


# ==========================================================================
# 1.  RUN SELECTION — which nights are "dense runs", and in what state
# ==========================================================================

def select_dense_runs(rows: Sequence[dict], state: Optional[str] = None
                      ) -> list[dict]:
    """The dense runs, optionally restricted to one accretion state.

    ``rows`` are ``cv_ext_verdict`` records as dicts -- the CV-S7 product
    that already decided, per night, whether the night is DENSE and what
    state the star was in, using independent AAVSO photometry wherever it
    exists.  This function does not re-decide either question.  It exists so
    that the selection is one auditable expression rather than a WHERE
    clause repeated in five places, and so a test can prove that changing
    the state filter changes the answer.

    The state labels come from CV-S7 (``QUIESCENT`` / ``OUTBURST`` / ...),
    and the density flag from its own frame-count and span rules.  Both are
    inherited on purpose: a stage that re-derived "is this night dense"
    could disagree with the page that decided the branch.
    """
    out = [dict(r) for r in rows if int(r.get("is_dense") or 0) == 1]
    if state is not None:
        out = [r for r in out if str(r.get("state") or "") == state]
    # Sorted by the LOCAL night, which is the key cv_frames uses.  CV-S7
    # tabulates by UTC night (one day later for these Arizona evenings), and
    # mixing the two conventions is the single easiest way to attribute a
    # quiescent measurement to an outburst night.
    return sorted(out, key=lambda r: str(r.get("local_night") or ""))


# ==========================================================================
# 2.  FOLDING AND THE ORBITAL HUMP
# ==========================================================================

def orbital_phase(times_d, period_d: float, epoch_bjd: float) -> np.ndarray:
    """Orbital phase in [0, 1) from a period and a stated zero point.

    ``epoch_bjd`` is an ARGUMENT and not a constant because YZ Cnc has no
    published epoch at all: the AAVSO VSX record gives a period and leaves
    the epoch blank, which is stored in ``p3_ephemeris`` as an explicit
    note.  Every phase this module produces is therefore relative to a zero
    point the caller chose, and the caller must say so.  What survives that
    limitation is real and useful -- phase DIFFERENCES between filters
    observed through the same night, and the SHAPE of the folded curve --
    and what does not survive is any claim about where the hump sits
    relative to a physical conjunction.
    """
    t = np.asarray(times_d, dtype=float)
    return np.mod((t - float(epoch_bjd)) / float(period_d), 1.0)


def harmonic_design(phase, n_harm: int = N_HARMONICS,
                    night_index=None) -> np.ndarray:
    """Design matrix ``[nightly constants | cos, sin, cos2, sin2, ...]``.

    Linear in its coefficients, which is why the fit below is a one-line
    least-squares solve with an exact covariance rather than an optimiser
    with a convergence story.  A folded light curve is one of the few places
    in this pipeline where the honest model really is linear.

    ``night_index`` turns the single constant into ONE CONSTANT PER NIGHT.
    That is what makes two nights foldable together: YZ Cnc faded 0.28 mag
    between 2024-05-01 and 2024-05-02, and a shared constant would push that
    step into the orbital fit as a spurious hump.  Fitting the nightly
    offsets JOINTLY with the harmonics -- rather than removing nightly means
    first and searching afterwards -- is the discipline CV-S9 measured the
    case for: detrend-then-search eats part of the signal it is looking for,
    and the joint fit does not.

    The returned column order is ``n_nights`` indicator columns followed by
    ``2 * n_harm`` trigonometric columns, so the fundamental is always at
    offsets ``n_nights`` and ``n_nights + 1``.
    """
    ph = np.asarray(phase, dtype=float)
    if night_index is None:
        cols = [np.ones_like(ph)]
    else:
        idx = np.asarray(night_index)
        cols = [(idx == u).astype(float) for u in np.unique(idx)]
    for k in range(1, int(n_harm) + 1):
        cols.append(np.cos(2 * np.pi * k * ph))
        cols.append(np.sin(2 * np.pi * k * ph))
    return np.vstack(cols).T


def weighted_lstsq(design, y, sigma) -> tuple[np.ndarray, np.ndarray, float]:
    """Weighted least squares: ``(coefficients, covariance, chi2 per dof)``.

    Weights are ``1/sigma**2`` with the sigmas the caller supplies, which in
    this pipeline means the per-point error bar already multiplied by the
    series' MEASURED chi-square inflation.  Returning chi2/dof alongside the
    covariance is not decoration: the covariance is only an error bar if the
    model fits, and every caller here rescales by ``sqrt(chi2/dof)`` when it
    exceeds 1 rather than quoting a formal bar the residuals contradict.
    """
    X = np.asarray(design, dtype=float)
    y = np.asarray(y, dtype=float)
    s = np.asarray(sigma, dtype=float)
    ok = np.isfinite(y) & np.isfinite(s) & (s > 0) & np.isfinite(X).all(axis=1)
    X, y, s = X[ok], y[ok], s[ok]
    n, p = X.shape
    if n <= p:
        nan = np.full(p, np.nan)
        return nan, np.full((p, p), np.nan), float("nan")
    w = 1.0 / s
    Xw, yw = X * w[:, None], y * w
    # lstsq rather than a normal-equation inverse: the design matrix is
    # 5 columns of trigonometry on a run that may cover barely one cycle,
    # which is exactly when the normal equations lose their conditioning.
    coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = yw - Xw @ coef
    dof = n - p
    chi2nu = float(resid @ resid / dof)
    try:
        cov = np.linalg.inv(Xw.T @ Xw)
    except np.linalg.LinAlgError:
        cov = np.full((p, p), np.nan)
    return coef, cov, chi2nu


def amplitude_and_phase(a: float, b: float) -> tuple[float, float]:
    """Convert ``a cos + b sin`` into ``(semi-amplitude, phase of maximum)``.

    ``a cos(2 pi ph) + b sin(2 pi ph) = A cos(2 pi (ph - ph_max))`` with
    ``A = hypot(a, b)`` and ``ph_max = atan2(b, a) / 2 pi``.  Phase of
    MAXIMUM in the fitted quantity; every caller here fits MAGNITUDES, so a
    maximum of the fitted curve is the star at its FAINTEST.  The callers
    negate deliberately and say which way round the answer is -- a hump is a
    brightening, and reporting its phase as the phase of maximum magnitude
    would be exactly backwards.
    """
    A = float(math.hypot(float(a), float(b)))
    ph = float(math.atan2(float(b), float(a)) / (2.0 * math.pi)) % 1.0
    return A, ph


def amplitude_sigma(a: float, b: float, cov_aa: float, cov_bb: float,
                    cov_ab: float) -> float:
    """One-sigma error on ``A = hypot(a, b)`` by linear propagation.

    ``dA/da = a/A``, ``dA/db = b/A``, so
    ``var(A) = (a^2 v_aa + b^2 v_bb + 2ab v_ab) / A^2``.  When ``A`` is
    small compared with its own error this expansion is not valid -- the
    amplitude of a non-detection has a Rice distribution, not a Gaussian --
    and the callers say NOT DETECTED from the injection contour instead of
    quoting this number as if it were a measurement.
    """
    A = math.hypot(float(a), float(b))
    if not np.isfinite(A) or A <= 0:
        return float("nan")
    var = (a * a * cov_aa + b * b * cov_bb + 2.0 * a * b * cov_ab) / (A * A)
    return float(math.sqrt(var)) if var > 0 else float("nan")


def fold_fit(times_d, mags, sigma, period_d: float, epoch_bjd: float,
             n_harm: int = N_HARMONICS, night_index=None) -> dict:
    """Fit a folded harmonic model to one run and report the hump.

    Returns a dict with

    ``amp``          fundamental SEMI-amplitude, magnitudes;
    ``amp_sigma``    its propagated error, rescaled by sqrt(chi2/dof) when
                     the fit is worse than the error bars predict;
    ``phase_max``    phase at which the star is BRIGHTEST (the hump peak);
    ``amp_harm``     second-harmonic semi-amplitude -- how non-sinusoidal;
    ``chi2nu``       goodness of fit;
    ``resid``        residuals about the fitted model, in the input order,
                     with NaN where the point was dropped.  This vector is
                     the input to the flickering statistics, so the coherent
                     orbital signal cannot be counted as flickering.

    The residual vector is returned rather than recomputed by the caller so
    that the model subtracted from the flickering analysis is provably the
    model quoted in the hump table.
    """
    t = np.asarray(times_d, dtype=float)
    y = np.asarray(mags, dtype=float)
    s = np.asarray(sigma, dtype=float)
    ph = orbital_phase(t, period_d, epoch_bjd)
    X = harmonic_design(ph, n_harm, night_index)
    n_const = 1 if night_index is None else int(
        np.unique(np.asarray(night_index)).size)
    coef, cov, chi2nu = weighted_lstsq(X, y, s)
    out = {"n": int(np.isfinite(y).sum()), "chi2nu": chi2nu,
           "period_d": float(period_d), "epoch_bjd": float(epoch_bjd),
           "n_const": n_const}
    if not np.isfinite(coef).all():
        out.update({"amp": float("nan"), "amp_sigma": float("nan"),
                    "phase_max": float("nan"), "amp_harm": float("nan"),
                    "mean_mag": float("nan"),
                    "resid": np.full(t.shape, np.nan)})
        return out
    # coef = [nightly constants..., a1, b1, a2, b2, ...]; the fundamental
    # sits immediately after the constants, whatever their number.
    i1 = n_const
    a1, b1 = float(coef[i1]), float(coef[i1 + 1])
    amp, ph_max_mag = amplitude_and_phase(a1, b1)
    # Rescaling the covariance by chi2/dof when chi2/dof > 1 is the standard
    # "the model does not fit, so widen the bar" convention.  It is applied
    # in ONE direction only: a chi2/dof below 1 does not license shrinking an
    # error bar, it means the input errors were generous.
    scale = math.sqrt(max(chi2nu, 1.0)) if np.isfinite(chi2nu) else 1.0
    amp_sig = amplitude_sigma(a1, b1, cov[i1, i1], cov[i1 + 1, i1 + 1],
                              cov[i1, i1 + 1]) * scale
    amp_h = (float(math.hypot(coef[i1 + 2], coef[i1 + 3]))
             if len(coef) >= i1 + 4 else float("nan"))
    # The fit is in MAGNITUDES, so the phase of the model's MAXIMUM is the
    # phase at which the star is FAINTEST.  The hump is a brightening, so
    # its peak sits half a cycle away from that maximum.
    out.update({"amp": amp, "amp_sigma": amp_sig,
                "phase_max": (ph_max_mag + 0.5) % 1.0,
                "amp_harm": amp_h, "mean_mag": float(np.mean(coef[:i1])),
                "resid": y - X @ coef})
    return out


def circular_mean_and_spread(phases) -> tuple[float, float]:
    """Mean phase and circular standard deviation, both in CYCLES.

    Phases live on a circle, so the ordinary mean of 0.98 and 0.02 is 0.5 --
    the point on the circle furthest from both.  The circular mean takes the
    argument of the mean unit vector instead, and the spread is
    ``sqrt(-2 ln R) / 2 pi`` with ``R`` the resultant length: zero when every
    phase agrees, and rising without bound as they scatter.

    This is the statistic that decides whether a modulation REPEATS.  A hump
    fitted independently in three filters through the same night should land
    at the same phase in all three, because it is one physical feature seen
    three ways; a hump fitted on two different nights should also agree, if
    it is the orbit.  When the first agrees and the second does not, the
    modulation is real and is not orbital -- which is a conclusion no
    amplitude on its own can reach.
    """
    ph = np.asarray(phases, dtype=float)
    ph = ph[np.isfinite(ph)]
    if ph.size == 0:
        return float("nan"), float("nan")
    ang = 2.0 * np.pi * ph
    C, S = float(np.mean(np.cos(ang))), float(np.mean(np.sin(ang)))
    R = math.hypot(C, S)
    mean = (math.atan2(S, C) / (2.0 * math.pi)) % 1.0
    if R <= 0.0:
        return mean, float("inf")
    return mean, float(math.sqrt(max(-2.0 * math.log(R), 0.0))
                       / (2.0 * math.pi))


def phase_difference(a: float, b: float) -> float:
    """Signed phase difference ``a - b`` wrapped into [-0.5, +0.5) cycles.

    Half a cycle is the largest difference that exists on a circle, so any
    routine that reports 0.62 is reporting 0.38 badly.
    """
    return float((float(a) - float(b) + 0.5) % 1.0 - 0.5)


def phase_drift_cycles(elapsed_d: float, period_d: float,
                       sigma_period_d: float) -> float:
    """Cycles of accumulated phase uncertainty over an elapsed baseline.

    ``N = elapsed / P`` cycles, each carrying ``sigma_P / P`` of phase
    error, so the accumulated drift is ``N * sigma_P / P``.  This is the
    number that decides whether two runs may share an absolute phase axis.
    For YZ Cnc the published period is quoted to four decimals with no
    uncertainty, so ``sigma_P`` is the quoted-precision FLOOR -- and the
    drift computed from a floor is itself a floor.
    """
    P = float(period_d)
    if P <= 0:
        return float("nan")
    return float(abs(elapsed_d) / P * abs(sigma_period_d) / P)


# ==========================================================================
# 3.  FLICKERING — amplitude against timescale, with a MEASURED floor
# ==========================================================================

def structure_function(times_d, values, tau_edges_s,
                       min_pairs: int = SF_MIN_PAIRS
                       ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """First-order structure function: ``(tau_centre, sigma(tau), n_pairs)``.

    For every pair of points separated by a lag inside a bin, accumulate
    ``(m_i - m_j)^2``; the bin's statistic is ``sqrt(mean / 2)``, which for
    a stationary process is the standard deviation of the variability on
    that timescale.  The factor of two is why: differencing two independent
    samples of the same process doubles the variance.

    A structure function rather than a periodogram because flickering is
    aperiodic, and rather than an Allan deviation because the Allan
    statistic bins in time and therefore mixes the cadence into the answer.
    Here the lag is the axis, and the cadence only limits which lags exist.

    Bins with fewer than ``min_pairs`` pairs return NaN.  They are not
    dropped from the arrays: a report that silently omitted its sparse bins
    would look like it measured a range it never sampled.
    """
    t = np.asarray(times_d, dtype=float) * 86400.0     # seconds
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(t) & np.isfinite(v)
    t, v = t[ok], v[ok]
    edges = np.asarray(tau_edges_s, dtype=float)
    centres = np.sqrt(edges[:-1] * edges[1:])          # log-centre of the bin
    sf = np.full(centres.size, np.nan)
    npair = np.zeros(centres.size, dtype=int)
    if t.size < 2:
        return centres, sf, npair
    # All pairs, once each.  These runs hold at most a few hundred points,
    # so the full O(N^2) upper triangle is ~50k entries -- cheaper than any
    # clever binning, and it cannot silently drop a pair.
    iu, ju = np.triu_indices(t.size, k=1)
    lag = np.abs(t[iu] - t[ju])
    dsq = (v[iu] - v[ju]) ** 2
    idx = np.digitize(lag, edges) - 1
    for b in range(centres.size):
        sel = idx == b
        n = int(sel.sum())
        npair[b] = n
        if n >= int(min_pairs):
            sf[b] = float(np.sqrt(np.mean(dsq[sel]) / 2.0))
    return centres, sf, npair


def sf_sigma(sf, n_pairs, n_points: int) -> np.ndarray:
    """Sampling error on a structure-function point.

    The variance of a variance estimated from ``n`` independent samples is
    ``2 sigma^4 / n``, so the relative error on ``sigma`` is
    ``1 / sqrt(2 n)``.  The catch is that structure-function pairs are NOT
    independent: every point appears in many pairs, so ``n_pairs`` badly
    overstates the information.  The effective count used here is capped at
    the NUMBER OF POINTS in the run, which is the number of independent
    measurements that exist at all.  That is deliberately conservative and
    is stated as such on the page -- a flickering detection that needs the
    optimistic pair count to clear three sigma is not a detection.
    """
    s = np.asarray(sf, dtype=float)
    n = np.minimum(np.asarray(n_pairs, dtype=float), float(max(n_points, 1)))
    with np.errstate(divide="ignore", invalid="ignore"):
        return s / np.sqrt(2.0 * np.maximum(n, 1.0))


def quadrature_excess(sf_total, sf_floor) -> np.ndarray:
    """Flickering amplitude: the floor removed in quadrature.

    ``sigma_flicker = sqrt(sigma_total^2 - sigma_floor^2)``, and NaN where
    the total sits below the floor.  Returning NaN rather than zero is the
    whole discipline of this function: a run whose scatter does not exceed
    the field stars' has not measured a small flickering amplitude, it has
    measured nothing, and a column of zeros would be read as the former.

    The subtraction is only valid if the two variances are independent and
    additive.  They are, by construction: the floor is measured on DIFFERENT
    STARS in the SAME FRAMES, so it carries the same photon statistics, the
    same sky, the same ensemble zero-point wander and the same atmosphere,
    and contains no flickering from the target.
    """
    a = np.asarray(sf_total, dtype=float)
    b = np.asarray(sf_floor, dtype=float)
    var = a * a - b * b
    return np.where(var > 0, np.sqrt(np.maximum(var, 0.0)), np.nan)


def excess_significance(sf_total, sig_total, sf_floor, sig_floor
                        ) -> np.ndarray:
    """How many sigma the VARIANCE excess is, bin by bin.

    The estimand is ``V = sigma_total^2 - sigma_floor^2``, whose error is
    ``sqrt((2 sigma_t s_t)^2 + (2 sigma_f s_f)^2)`` by propagation through
    the squares.  Working in variance rather than in sigma matters: the
    difference of two similar sigmas is a badly behaved statistic near zero,
    while the difference of two variances is the quantity that is actually
    additive and whose null is symmetric about zero.
    """
    a = np.asarray(sf_total, dtype=float)
    sa = np.asarray(sig_total, dtype=float)
    b = np.asarray(sf_floor, dtype=float)
    sb = np.asarray(sig_floor, dtype=float)
    num = a * a - b * b
    den = np.sqrt((2.0 * a * sa) ** 2 + (2.0 * b * sb) ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, np.nan)


def roll_within_nights(values, night_index, rng) -> np.ndarray:
    """Cyclically roll a residual series, EACH NIGHT INDEPENDENTLY.

    The prayer-bead trick -- rolling residuals rather than shuffling them --
    keeps the noise's correlated structure while decorrelating it from any
    signal's phase, and it is what makes a red-noise false-alarm threshold
    possible at all.  Rolling the WHOLE vector across a one-day gap does
    something else, though: it wraps the second night's residuals onto the
    first night's timestamps, manufacturing a step at the join that no
    realization of this noise would ever contain, and the step's power leaks
    into the orbital band and inflates the threshold.

    Rolling night by night keeps every realization made of pieces the
    telescope actually recorded on that night, at that airmass, in that
    seeing.  On a single-night scope this is identical to
    :func:`macro_phot.characterize.cyclic_noise_realization`.
    """
    v = np.asarray(values, dtype=float)
    idx = np.asarray(night_index)
    out = np.array(v, dtype=float)
    for u in np.unique(idx):
        sel = idx == u
        n = int(sel.sum())
        if n > 1:
            out[sel] = np.roll(v[sel], int(rng.integers(n)))
    return out


def median_floor(curves: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Median structure function over the field stars, and its own spread.

    Returns ``(median, sigma_of_the_median)`` where the second is
    ``1.253 * MAD-based sigma / sqrt(n)`` -- the standard error of a median
    for a roughly normal population.  The MEDIAN and not the mean because a
    single variable star that slipped into the magnitude-matched sample
    would drag a mean floor upward and, through the quadrature subtraction,
    ERASE a real flickering detection.  A median cannot be moved that way by
    one star in twenty.
    """
    if not curves:
        return np.array([]), np.array([])
    M = np.vstack([np.asarray(c, dtype=float) for c in curves])
    with np.errstate(invalid="ignore"):
        med = np.nanmedian(M, axis=0)
        mad = np.nanmedian(np.abs(M - med), axis=0)
    n = np.sum(np.isfinite(M), axis=0).astype(float)
    # 1.4826 converts a MAD to a Gaussian sigma; sqrt(pi/2) = 1.2533 converts
    # the standard error of a mean into that of a median.
    sig = 1.4826 * mad * 1.2533 / np.sqrt(np.maximum(n, 1.0))
    return med, np.where(n >= 2, sig, np.nan)


# ==========================================================================
# 4.  OUTBURST RUNS — what a single night of an outburst can say
# ==========================================================================

def linear_rate_per_hour(times_d, mags) -> tuple[float, float]:
    """Straight-line brightness rate through one run: ``(mag/h, sigma)``.

    The simplest statement an outburst run can make that a single nightly
    mean cannot: is the star rising, flat or declining, and how fast.  The
    error is the ordinary least-squares slope error rescaled by the fit's
    own residual scatter, so a run with structure gets a wide bar rather
    than a precise slope through a curve.
    """
    t = np.asarray(times_d, dtype=float)
    y = np.asarray(mags, dtype=float)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 3:
        return float("nan"), float("nan")
    h = (t - t.mean()) * 24.0
    Sxx = float(h @ h)
    if Sxx <= 0:
        return float("nan"), float("nan")
    slope = float(h @ (y - y.mean()) / Sxx)
    resid = y - (y.mean() + slope * h)
    s2 = float(resid @ resid / (t.size - 2))
    return slope, float(math.sqrt(s2 / Sxx))


def percentile_amplitude(mags, lo: float = 5.0, hi: float = 95.0) -> float:
    """Peak-to-peak amplitude as a percentile range, magnitudes.

    p5-p95 rather than max-min because max-min is a statistic of the two
    worst points in the run and grows with the number of points even in pure
    noise.  This is the same definition CV-S9's state classifier uses, so
    the two pages' amplitudes are the same quantity.
    """
    m = np.asarray(mags, dtype=float)
    m = m[np.isfinite(m)]
    if m.size < 5:
        return float("nan")
    return float(np.percentile(m, hi) - np.percentile(m, lo))


# ==========================================================================
# 5.  GATES AND VERDICTS
# ==========================================================================

def snr_gate(sigma_measured_mag: float, amplitude_to_detect_mag: float,
             required_ratio: float = 5.0) -> dict:
    """The strategy's §4.19 gate, as arithmetic.

    §4.19 refuses to promise the quiescent fallback until someone shows that
    8 s High Gain frames at quiescent V ~ 14.5 are not sky- or read-noise
    dominated.  "Not dominated" is not a measurable statement on its own, so
    it is turned into one here: the per-point precision actually achieved on
    the quiescent frames must be at least ``required_ratio`` times smaller
    than the modulation the fallback exists to measure.

    Returns ``{'ratio', 'passes', 'margin_mag'}``.  The caller supplies both
    numbers from products -- the precision from the measured error model,
    the amplitude from the fitted hump -- so the gate cannot be passed by
    choosing a flattering assumption about either.
    """
    s = float(sigma_measured_mag)
    a = float(amplitude_to_detect_mag)
    if not (np.isfinite(s) and np.isfinite(a)) or s <= 0:
        return {"ratio": float("nan"), "passes": False,
                "margin_mag": float("nan")}
    ratio = a / s
    return {"ratio": float(ratio), "passes": bool(ratio >= required_ratio),
            "margin_mag": float(a - required_ratio * s)}


#: Orbital cycles a scope must cover before a fitted modulation at P_orb may
#: be called a DETECTION at all.  Below about one and a half cycles a
#: sinusoid at the orbital period and a smooth nightly trend span nearly the
#: same subspace, so the fit will return an amplitude whatever the star is
#: doing.  1.5 is where the two become separable in the fitted covariance;
#: scopes below it report an AMPLITUDE ONLY and say so.
MIN_CYCLES_FOR_DETECTION = 1.5


def detection_call(amplitude_mag: float, contour_amp90_mag: float,
                   power_at_forb: float, threshold_power: float,
                   cycles: float = float("inf")) -> str:
    """DETECTED / MARGINAL / NOT DETECTED, from two independent tests.

    Both tests must agree before the word DETECTED is used:

    * the measured semi-amplitude must reach the 90% RECOVERY CONTOUR -- the
      amplitude at which nine injections in ten into this scope's own
      timestamps and this scope's own noise come back;
    * the Lomb-Scargle power AT the published orbital frequency must clear a
      1% false-alarm threshold measured on signal-free realisations of the
      same noise.

    They fail differently, which is why both are here.  The contour can be
    beaten by one lucky excursion; the power threshold can be beaten by red
    noise that happens to have power near the orbital frequency.  When they
    disagree the answer is MARGINAL, and a MARGINAL hump is not a hump.

    Ahead of both sits a coverage gate: a scope spanning less than
    :data:`MIN_CYCLES_FOR_DETECTION` orbits returns AMPLITUDE ONLY however
    well it scores, because on such a scope "modulation at P_orb" and
    "trend" are not distinguishable statements about the sky.
    """
    if np.isfinite(cycles) and float(cycles) < MIN_CYCLES_FOR_DETECTION:
        return "AMPLITUDE ONLY"
    a, c = float(amplitude_mag), float(contour_amp90_mag)
    p, thr = float(power_at_forb), float(threshold_power)
    if not np.isfinite(a) or not np.isfinite(c) or not np.isfinite(p) \
            or not np.isfinite(thr):
        return "UNTESTED"
    by_amp, by_power = a >= c, p > thr
    if by_amp and by_power:
        return "DETECTED"
    if by_amp or by_power:
        return "MARGINAL"
    return "NOT DETECTED"


def run_night_label(utc_nights, nights) -> str:
    """The ONE name a dense run is called by in print: its first UTC night.

    ``p4_run.scope`` keys a run by its local OBSERVING night
    (``yzcnc|e7|I|2024-02-20``) while ``p4_run.utc_nights`` records the UTC
    night the frames actually carry (``2024-02-21``).  Both conventions are
    legitimate and the release keeps both, but a paper may use only one:
    Table 4 named the uninformative YZ Cnc run "the I run of 2024-02-20"
    while Figure 11's caption and its own axis label called the same run
    "the 2024-02-21 I run", and a reader comparing the two saw two runs
    where there is one.

    UTC is the convention that survives into print, because every time in
    this paper is a BJD_TDB derived from a UTC timestamp and an observing
    night is a local, site-dependent bookkeeping label.  ``nights`` is the
    fallback for a row that predates ``utc_nights`` being populated.

    Returns the first night of the run; a block scope spanning several
    nights is joined with ``+`` upstream and only its first is named.
    """
    night = utc_nights if utc_nights else nights
    return str(night or "").split("+")[0]


def capability_verdict(measured: float, bar: float,
                       higher_is_better: bool = True) -> str:
    """One measured number against one stated bar -> one word.

    Deliberately trivial.  Every genuinely contestable choice in the AN UMa
    decision is a choice of BAR, and a bar is a number a reader can disagree
    with in one line.  Burying the same choice inside a weighted score would
    make the decision unarguable, which for a go/no-go is the opposite of
    what is wanted.
    """
    m, b = float(measured), float(bar)
    if not np.isfinite(m) or not np.isfinite(b):
        return "UNTESTED"
    ok = (m >= b) if higher_is_better else (m <= b)
    return "SUPPORTED" if ok else "NOT SUPPORTED"


def nights_needed(have: int, bar: int, nights_per_season: float) -> dict:
    """What it would take to move a NOT-SUPPORTED verdict.

    Returns ``{'shortfall', 'seasons'}``.  ``nights_per_season`` is measured
    from the archive (qualifying nights divided by seasons observed), so the
    seasons figure is the archive's own delivered rate rather than an
    optimistic scheduling assumption.  A verdict that cannot say what would
    change it is an opinion.
    """
    short = max(int(bar) - int(have), 0)
    rate = float(nights_per_season)
    seasons = float(short) / rate if rate > 0 else float("inf")
    return {"shortfall": short, "seasons": seasons}


def log_tau_edges(tau_min_s: float, tau_max_s: float, n_bins: int
                  ) -> np.ndarray:
    """Logarithmically spaced structure-function bin edges.

    Log spacing because flickering is broadband: its amplitude is expected
    to vary as a power of the timescale, and a power law is a straight line
    only on a log axis.  Linear bins would put nearly every pair in the last
    bin and leave the short lags -- the ones a 140-200 s cadence can barely
    reach -- with no pairs at all.
    """
    return np.geomspace(float(tau_min_s), float(tau_max_s), int(n_bins) + 1)
