"""CV-R1...R8 — the arithmetic of the committee's major revision.

WHY THIS MODULE EXISTS
----------------------
The plan review of 2026-10-03 (``committee/reviews/2026-10-03/``) reopened
the CV paper.  Its statistics were arithmetically right and inferentially
wrong in four places, and every one of them was the same kind of mistake:
an error bar or a summary that could only ever make a result look weaker.

*   The inter-band edge offset was declared null on a 3-sigma VOTE COUNT
    whose error bars came from a transported budget inflated by
    ``max(chi2_nu, 1)``; the paired differences themselves scatter half as
    much as those bars and their mean is four standard errors from zero.
*   The "bias floor of 5.5 s" was the MEDIAN OF ABSOLUTE VALUES over an
    injection grid whose matched cell held a signed bias of -50 to -80 s.
*   The "per-edge Monte-Carlo sigma" was one constant per series.
*   The O-C pooled bands that differ by ~100 s, counted 36 epochs on 17
    nights as 36 independent data, and had no nuisance term for the camera
    change between the two seasons.

This module holds the replacements, as pure functions of arrays:

1.  **Distribution-free paired tests** (:func:`paired_tests`,
    :func:`signflip_p`, :func:`wilcoxon_signed_rank_p`, :func:`sign_test_p`)
    with SCATTER-BASED errors, a night-clustered variant because cycles of
    one night are not independent, and a Bonferroni trials factor.
2.  **The edge estimators, written down** (:func:`fit_edge_profile`,
    :func:`edge_window`, :func:`to_relative_flux`): the published v1
    estimator is :func:`macro_phot.phase3.fit_edge` in magnitudes on a
    per-band width grid; the v2 robustness estimator is the same ramp
    fitted in FLUX on one width grid common to every band, with a
    symmetric rule for flat chi-squared valleys.  Flux is the space in
    which the half-level point of an edge does not move when a constant
    (unocculted) light source dilutes one band more than another.
3.  **Injection with the real estimator** (:func:`inject_recover`): an
    ACHROMATIC flux ramp is injected at the real timestamps of every band
    of a night, with that band's own depth, errors and real residuals, and
    recovered by the estimator that band actually used.  The output is the
    SIGNED bias per band and the SIGNED differential bias per band pair,
    each with a standard error — never a median of absolute values.
4.  **A per-edge bootstrap** (:func:`bootstrap_edge`): a wild residual
    bootstrap of each real edge fit, which is what "a Monte-Carlo error on
    every edge" has to mean if the phrase is to be used at all.
5.  **The O-C model with its nuisance terms** (:func:`fit_oc_model`,
    :func:`night_level_epochs`, :func:`pdot_from_quadratic`): per-band
    constants, an era offset, and night-level epochs, every variant
    reporting chi-squared WITH its degrees of freedom and both the
    budget-propagated and the scatter-based error on the quadratic term.
6.  **Spot longitude** (:func:`seconds_to_degrees`) and the three physical
    scales a period-derivative bound is to be read against
    (:func:`gr_orbital_pdot`).
7.  **The colour curve as a measurement** (:func:`interpolated_colour`,
    :func:`colour_curve_summary`, :func:`curve_repeatability`): amplitude,
    phases of the extrema, errors from a night bootstrap, and an estimator
    that is immune to the non-simultaneity a pairing window allows.
8.  **Excluded-amplitude statements** (:func:`excluded_amplitude_summary`):
    what a 90 per cent recovery contour does and does not exclude.

Everything here is deterministic given its seed and touches no file, no
database and no network.  ``pipeline/scripts/run_cv_revision.py`` does the
I/O; ``pipeline/tests/test_revision_cv.py`` injects known offsets, known
biases and known colour curves and demands them back.

THE STANDING RULES THIS MODULE IS WRITTEN TO (SYNTHESIS section 5)
-------------------------------------------------------------
* Every chi-squared is returned with its degrees of freedom, and nothing
  in this file contains ``max(chi2_nu, 1)``.  A chi2_nu of 0.3 is reported
  as 0.3: it says the errors are over-stated by 1.8x, which is a defect.
* Every null carries the effect size it would have recovered.
* Injection results are SIGNED, per matched cell.
"""

from __future__ import annotations

import math
from typing import Callable, Optional, Sequence

import numpy as np

from . import phase3 as p3

# ===========================================================================
# 0.  Constants — every choice this revision makes, in one place
# ===========================================================================

#: Canonical band SLOT of every filter label in the ST LMi series.  The
#: 2024 High Gain camera's wheel is labelled ``G/R/I`` and the 2025 Mode0
#: camera's ``g/r/i``.  They are treated as the same three slots when eras
#: are combined, and that is an ASSUMPTION about the glass (stated in the
#: report and testable only as "both eras give the same sign and size"),
#: not a measurement.
BAND_SLOT = {"G": "g", "R": "r", "I": "i", "g": "g", "r": "r", "i": "i",
             "z": "z", "Z": "z", "y": "y", "Y": "y"}

#: The three band pairs, bluer band first, so a NEGATIVE difference always
#: means "the bluer band's edge comes first".
BAND_PAIRS = (("g", "i"), ("g", "r"), ("r", "i"))

#: Trials factor for the band-offset tests: three pairs are examined.  Only
#: two are independent (g-i = g-r + r-i), so multiplying by three is
#: conservative, which is the direction a trials factor should err in.
N_BAND_PAIR_TRIALS = 3

#: v2 estimator: ONE ramp-width grid for every band, in seconds.  Eight
#: geometric steps from well under the 60 s exposure to 600 s, which is
#: longer than any ramp the folded profiles show (0.09 cycle).  The v1
#: estimator took its grid from each series' own folded profile, which
#: handed g a grid ending at 98 s and i one ending at 448 s; a band offset
#: measured with band-dependent model freedom cannot be told from one
#: produced by it, so v2 removes the freedom.
V2_WIDTH_GRID_S = tuple(float(x) for x in np.geomspace(25.0, 600.0, 8))

#: v2 time-grid step, seconds.  Small against every error in the problem,
#: so grid quantisation (step / sqrt(12) = 0.6 s) leaves the budget.
V2_TIME_STEP_S = 2.0

#: Relative tolerance inside which two chi-squared values are THE SAME for
#: the purpose of locating a flat valley.  A ramp that lies wholly between
#: two exposures gives exactly the same model wherever it sits in the gap,
#: so chi-squared is flat there to rounding; ``np.argmin`` then returns
#: whichever node rounding favours.  v2 takes the valley's midpoint.
VALLEY_RTOL = 1e-9

#: Scatter, in seconds, of the ephemeris prediction about the true edge
#: that the injection test gives the estimator as its starting guess.  The
#: real pipeline centres its window and its grid on the ephemeris, and the
#: measured O-C scatter of real edges about it is 84 s; using zero would
#: hand the recovery a better guess than the real fits ever had.
GUESS_SCATTER_S = 90.0

#: True ramp widths injected, seconds, in FLUX.  The real width is not
#: known per band (the fits put it at or below the cadence in g and at one
#: to two cadences in r and i), so the bias is tabulated across the range
#: and the verdict has to hold over it rather than at a chosen point.
INJECT_WIDTHS_S = (30.0, 60.0, 120.0, 240.0, 480.0)

#: Realisations per injection cell and night.  500 gives a standard error
#: on one night's bias of about sigma_t / 22, i.e. 4 s at sigma_t = 90 s,
#: and less once nights are combined — small against the ~100 s offset the
#: test exists to explain or fail to explain.
N_INJECT = 500

#: Bootstrap replicates per real edge.
N_EDGE_BOOT = 400

#: Seed for every random draw in the revision.  Same value as Phase 3.
SEED = p3.SEED

#: Enumerate sign flips exactly up to this many units (2^20 = 1,048,576
#: patterns, about a second); above it, Monte Carlo with N_PERM_MC draws.
EXACT_PERM_MAX_N = 20
N_PERM_MC = 200_000

#: Phase bins of the colour curve.  Twenty bins of a 6,833 s orbit are
#: 342 s wide: one and a half cadences, so every bin of a dense night holds
#: several pairs and the bin is still narrower than the bright phase edges
#: are apart.
COLOUR_BINS = 20

#: Minimum pairs in a colour bin before the bin is used.
COLOUR_MIN_PER_BIN = 5

#: Night-bootstrap replicates for the colour-curve summary.
N_COLOUR_BOOT = 500

#: The strict pairing window the referee asked the colour result to
#: survive, seconds.
STRICT_PAIR_WINDOW_S = 120.0

#: ---------------------------------------------------------------------
#: LITERATURE CONSTANTS.  None of these is a measurement of this
#: programme.  Each was read and verified by the ``cv-literature`` package
#: (``committee/work/cv-literature/literature_scales.md``, which is the
#: record of where in each paper the value stands and how it was read);
#: the BibTeX key of the source is given beside each.  They are carried
#: into ``rv_result`` with ``origin='literature'`` so that the paper's
#: one-query separation of measurements from everything else still works.
#: ---------------------------------------------------------------------

#: Masses for the gravitational-radiation scale: a representative polar
#: below the period gap (knigge2011's donor sequence at 1.9 h).
GR_M1_MSUN = 0.75
GR_M2_MSUN = 0.17

#: Secular orbital evolution as P / tau for tau = 1 and 5 Gyr (knigge2011):
#: the range the literature package adopts as "what the orbit does".
SECULAR_TAU_GYR = (1.0, 5.0)

#: Spin-period derivative of the asynchronous polar V1500 Cyg
#: (schmidt1995, as quoted by pavlenko2018).
ASYNC_POLAR_PDOT = 3.86e-8

#: A DP Leo-like libration of the accretion spot in longitude: amplitude
#: in degrees and period in years (beuermann2014).
LIBRATION_AMP_DEG = 25.0
LIBRATION_PERIOD_YR = 60.0

#: The only PUBLISHED uncertainty on ST LMi's period, in days
#: (cropper1986, Sect. 3: 0.07908908 +/- 0.00000008 d).  The VSX period
#: the pipeline folds on carries no published error at all (it was set
#: "from AAVSO data" by the VSX moderator in 2025), so the first draft's
#: cycle-count drift rested on an ASSUMED half-last-digit error.
CROPPER_PERIOD_SIGMA_D = 8.0e-8

#: Bright-phase duration of ST LMi in 1982 and in 1985, in cycles
#: (cropper1986, Sect. 4).  Half the change is how far ONE edge moved:
#: the historical scale of an accretion-geometry change, and the size a
#: state-dependent longitude shift has to be compared with.
BRIGHT_DURATION_1982 = 0.33
BRIGHT_DURATION_1985 = 0.38

#: Bright-phase duration at half maximum in white light and in J, in
#: cycles, measured simultaneously (bailey1985, Sect. 7.3): the bright
#: phase is LONGER at longer wavelength, so the redder band's edge comes
#: later.  Effective wavelengths in micron for the interpolation.
BRIGHT_DURATION_WHITE = 0.29
BRIGHT_DURATION_J = 0.31
BRIGHT_DURATION_K = 0.39
WAVELENGTH_WHITE_UM = 0.55
WAVELENGTH_J_UM = 1.25

#: Magnetic field of ST LMi's main accretion region, MG (campbell2008c:
#: 12.1 +/- 0.5), and the effective wavelengths of the band slots, A.
STLMI_FIELD_MG = 12.1
BAND_WAVELENGTH_A = {"g": 4810.0, "r": 6170.0, "i": 7520.0}

#: Superhump semi-amplitudes, mag: the maximum at low inclination
#: (smak2010; kato2012: 0.25 mag full) and YZ Cnc's own in TESS at the
#: precursor and early plateau (dai2026: 0.3 mag full).  These replace the
#: first draft's uncited "50 mmag floor" as the amplitudes a superhump
#: search has to be sensitive to.
SUPERHUMP_SEMI_PEAK_MAG = 0.125
SUPERHUMP_SEMI_TESS_MAG = 0.150

#: YZ Cnc: the spectroscopic period with its error (vanparadijs1994,
#: Sect. 6), the error of the period VSX carries (shafter1988, as quoted),
#: and the semi-amplitude of the orbital hump van Paradijs et al. detected
#: in quiescence (0.5 mag full).
YZCNC_VP_PERIOD_D = 0.086924
YZCNC_VP_PERIOD_SIGMA_D = 7.0e-6
YZCNC_SH_PERIOD_SIGMA_D = 2.0e-4
YZCNC_PRIOR_HUMP_SEMI_MAG = 0.250

_G = 6.67430e-11
_C = 2.99792458e8
_MSUN = 1.98847e30


# ===========================================================================
# 1.  Small-sample distributions, without scipy
# ===========================================================================
# The pipeline does not depend on scipy and three p-values are not a reason
# to start.  Each function below is exact or a standard continued fraction,
# and the test file checks all of them against tabulated values.

def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (Lentz)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    d = tiny if abs(d) < tiny else d
    d = 1.0 / d
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + aa / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = tiny if abs(d) < tiny else d
        c = 1.0 + aa / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_front = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                + a * math.log(x) + b * math.log(1.0 - x))
    front = math.exp(ln_front)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def student_t_two_sided_p(t_stat: float, dof: float) -> float:
    """Two-sided p-value of Student's t with ``dof`` degrees of freedom."""
    if not (np.isfinite(t_stat) and dof > 0):
        return float("nan")
    x = dof / (dof + t_stat * t_stat)
    return float(betainc(dof / 2.0, 0.5, x))


def chi2_sf(chi2: float, dof: float) -> float:
    """Upper tail probability of chi-squared: P(X >= chi2 | dof).

    Via the regularised upper incomplete gamma function, series for small
    arguments and a continued fraction for large ones.
    """
    if not (np.isfinite(chi2) and dof > 0):
        return float("nan")
    if chi2 <= 0:
        return 1.0
    a, x = dof / 2.0, chi2 / 2.0
    lg = math.lgamma(a)
    if x < a + 1.0:
        # Series for the lower function P(a, x); return 1 - P.
        term = 1.0 / a
        total = term
        n = a
        for _ in range(2000):
            n += 1.0
            term *= x / n
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return float(max(0.0, 1.0 - total * math.exp(-x + a * math.log(x)
                                                     - lg)))
    # Continued fraction for Q(a, x) (modified Lentz).
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return float(min(1.0, math.exp(-x + a * math.log(x) - lg) * h))


def sign_test_p(n_negative: int, n_nonzero: int) -> float:
    """Exact two-sided binomial sign test.

    ``n_negative`` of ``n_nonzero`` differences are below zero; under the
    null each sign is a fair coin.  Two-sided: twice the smaller tail,
    capped at one.
    """
    n, k = int(n_nonzero), int(n_negative)
    if n <= 0:
        return float("nan")
    k_small = min(k, n - k)
    tail = sum(math.comb(n, j) for j in range(0, k_small + 1)) / 2.0 ** n
    return float(min(1.0, 2.0 * tail))


def wilcoxon_signed_rank_p(d) -> float:
    """Exact two-sided Wilcoxon signed-rank p-value.

    Zeros are dropped (Wilcoxon's own convention); ties in |d| get average
    ranks, and the null distribution is then enumerated for THOSE ranks by
    dynamic programming over doubled ranks, so the p-value is exact with
    ties rather than a normal approximation.  Exact enumeration is cheap
    for every sample this project has (n <= 40).
    """
    x = np.asarray(d, dtype=float)
    x = x[np.isfinite(x) & (x != 0.0)]
    n = x.size
    if n == 0:
        return float("nan")
    order = np.argsort(np.abs(x), kind="mergesort")
    ranks = np.empty(n)
    a = np.abs(x)[order]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and a[j + 1] == a[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    # Doubled ranks are integers even when ties produce half-ranks.
    r2 = np.rint(2.0 * ranks).astype(int)
    total = int(r2.sum())
    w_plus = int(r2[x > 0].sum())
    # counts[s] = number of sign patterns whose positive doubled-rank sum
    # is s.  Python integers, so no overflow at any n.
    counts = [0] * (total + 1)
    counts[0] = 1
    for r in r2:
        for s in range(total, r - 1, -1):
            counts[s] += counts[s - r]
    n_pat = 2 ** n
    lo = min(w_plus, total - w_plus)
    tail = sum(counts[:lo + 1])
    return float(min(1.0, 2.0 * tail / n_pat))


def signflip_p(d, clusters=None, seed: int = SEED) -> tuple[float, str]:
    """Two-sided sign-flip (label-permutation) p-value for a paired mean.

    Under the null "the two bands' edges coincide" each paired difference
    is as likely to carry one sign as the other, so swapping the two band
    LABELS within a pair — flipping the difference's sign — leaves the
    data's distribution unchanged.  The p-value is the fraction of sign
    patterns whose |mean| is at least the observed one.  No error model,
    no Gaussian assumption, no budget.

    ``clusters`` (one label per difference) makes the flip act on whole
    clusters: every pair of one NIGHT flips together.  That is the honest
    version when cycles within a night share a seeing, a sky and an
    accretion state, because it does not credit the test with more
    independent units than there are nights.

    Returns ``(p, basis)`` where ``basis`` says whether the enumeration was
    exact, and over how many patterns.
    """
    x = np.asarray(d, dtype=float)
    ok = np.isfinite(x)
    x = x[ok]
    if x.size == 0:
        return float("nan"), "no data"
    if clusters is not None:
        lab = np.asarray(clusters)[ok]
        _, inv = np.unique(lab, return_inverse=True)
        units = np.bincount(inv, weights=x)
    else:
        units = x
    m = units.size
    obs = abs(float(units.sum()))
    # A tolerance on the comparison, so the identity pattern (and its
    # mirror) always count as "at least as extreme" despite rounding.
    tol = 1e-12 * max(1.0, float(np.sum(np.abs(units))))
    if m <= EXACT_PERM_MAX_N:
        # Enumerate all 2^m patterns in blocks, by bit arithmetic.
        n_pat = 1 << m
        count = 0
        block = 1 << 16
        bits = np.arange(m, dtype=np.int64)
        for start in range(0, n_pat, block):
            idx = np.arange(start, min(start + block, n_pat),
                            dtype=np.int64)
            signs = 1.0 - 2.0 * ((idx[:, None] >> bits) & 1)
            count += int(np.sum(np.abs(signs @ units) >= obs - tol))
        return float(count / n_pat), f"exact, {n_pat} sign patterns"
    rng = np.random.default_rng(seed)
    count = 0
    done = 0
    while done < N_PERM_MC:
        k = min(20_000, N_PERM_MC - done)
        signs = rng.choice((-1.0, 1.0), size=(k, m))
        count += int(np.sum(np.abs(signs @ units) >= obs - tol))
        done += k
    # +1 in both: the observed pattern is itself a draw from the null.
    return float((count + 1) / (N_PERM_MC + 1)), (
        f"Monte Carlo, {N_PERM_MC} sign patterns")


def paired_tests(d, clusters=None, sigma_budget=None,
                 n_trials: int = 1, seed: int = SEED) -> dict:
    """Everything a paired band difference is tested with, in one record.

    ``d`` are the per-unit differences (seconds).  The record carries:

    ``mean``, ``sd``, ``se``
        the sample mean, the sample standard deviation and the
        SCATTER-BASED standard error ``sd / sqrt(n)``.  This is the error
        the committee asked for: it comes from how much the differences
        actually disagree with each other, not from a budget.
    ``t``, ``p_t``
        Student's t on ``n - 1`` degrees of freedom.
    ``n_negative``, ``p_sign``
        the sign test.
    ``p_wilcoxon``, ``p_perm``
        the exact signed-rank test and the sign-flip permutation test on
        the mean (see :func:`signflip_p`).
    ``p_perm_cluster``, ``n_clusters``
        the same permutation test flipping whole clusters (nights).
    ``mean_budget``, ``sigma_budget``, ``chi2``, ``dof``, ``chi2nu``
        when ``sigma_budget`` is given: the inverse-variance weighted mean
        under those per-unit errors, its PROPAGATED error with no
        rescaling in either direction, and the chi-squared of the
        differences about that mean with its degrees of freedom.  A
        chi2_nu far below one says the budget errors are over-stated; it
        is reported, not clipped.
    ``p_*_bonf``
        each p-value multiplied by ``n_trials`` and capped at one.
    """
    x = np.asarray(d, dtype=float)
    ok = np.isfinite(x)
    x = x[ok]
    n = int(x.size)
    nan = float("nan")
    out = {"n": n, "mean": nan, "median": nan, "sd": nan, "se": nan,
           "t": nan, "p_t": nan, "n_negative": 0, "n_positive": 0,
           "p_sign": nan, "p_wilcoxon": nan, "p_perm": nan,
           "perm_basis": "", "p_perm_cluster": nan, "n_clusters": None,
           "mean_budget": nan, "sigma_budget": nan, "chi2": nan,
           "dof": None, "chi2nu": nan, "n_trials": int(n_trials)}
    if n == 0:
        return out
    out["mean"] = float(x.mean())
    out["median"] = float(np.median(x))
    out["n_negative"] = int(np.sum(x < 0))
    out["n_positive"] = int(np.sum(x > 0))
    out["p_sign"] = sign_test_p(out["n_negative"],
                                out["n_negative"] + out["n_positive"])
    if n > 1:
        sd = float(x.std(ddof=1))
        out["sd"] = sd
        out["se"] = sd / math.sqrt(n)
        if sd > 0:
            out["t"] = out["mean"] / out["se"]
            out["p_t"] = student_t_two_sided_p(out["t"], n - 1)
        out["p_wilcoxon"] = wilcoxon_signed_rank_p(x)
    out["p_perm"], out["perm_basis"] = signflip_p(x, seed=seed)
    if clusters is not None:
        lab = np.asarray(clusters)[ok]
        out["n_clusters"] = int(np.unique(lab).size)
        out["p_perm_cluster"], _ = signflip_p(x, clusters=lab, seed=seed)
    if sigma_budget is not None:
        s = np.asarray(sigma_budget, dtype=float)[ok]
        good = np.isfinite(s) & (s > 0)
        if good.sum() >= 1:
            w = 1.0 / s[good] ** 2
            mb = float(np.sum(w * x[good]) / np.sum(w))
            out["mean_budget"] = mb
            out["sigma_budget"] = float(math.sqrt(1.0 / np.sum(w)))
            if good.sum() > 1:
                chi2 = float(np.sum(w * (x[good] - mb) ** 2))
                out["chi2"] = chi2
                out["dof"] = int(good.sum() - 1)
                out["chi2nu"] = chi2 / out["dof"]
    for key in ("p_t", "p_sign", "p_wilcoxon", "p_perm", "p_perm_cluster"):
        v = out[key]
        out[key + "_bonf"] = (float(min(1.0, v * n_trials))
                              if np.isfinite(v) else nan)
    return out


def combine_inverse_variance(means, ses) -> dict:
    """Inverse-variance combination of independent estimates, with the
    heterogeneity chi-squared that says whether combining was legitimate.

    Used to combine the two instrument eras' band offsets.  ``chi2_het`` on
    ``dof_het = k - 1`` tests "the eras measure the same offset"; it is
    returned with its p-value so a combined number is never quoted without
    the evidence that the things combined agree.
    """
    m = np.asarray(means, dtype=float)
    s = np.asarray(ses, dtype=float)
    ok = np.isfinite(m) & np.isfinite(s) & (s > 0)
    m, s = m[ok], s[ok]
    nan = float("nan")
    if m.size == 0:
        return {"mean": nan, "se": nan, "chi2_het": nan, "dof_het": 0,
                "p_het": nan, "k": 0}
    w = 1.0 / s ** 2
    mean = float(np.sum(w * m) / np.sum(w))
    se = float(math.sqrt(1.0 / np.sum(w)))
    chi2 = float(np.sum(w * (m - mean) ** 2))
    dof = int(m.size - 1)
    return {"mean": mean, "se": se, "chi2_het": chi2, "dof_het": dof,
            "p_het": chi2_sf(chi2, dof) if dof > 0 else nan,
            "k": int(m.size)}


# ===========================================================================
# 2.  The edge model, written down
# ===========================================================================
# THE MODEL (this is the text CV-R7 asks the paper to print).
#
# Within a window of +/- 0.17 cycle about the predicted egress the light
# curve y(t) is modelled as two constant levels joined by a linear ramp:
#
#     y(t) = L_b + (L_f - L_b) * clip( (t - t_e)/w + 1/2 , 0, 1 )
#
# with four parameters: the bright level L_b, the faint level L_f, the
# ramp MIDPOINT t_e (the published "edge time") and the ramp width w.  For
# every (t_e, w) on a grid the two levels follow from a weighted linear
# least-squares solve; t_e and w are then chosen by the minimum of
# chi-squared over the grid.  Nothing is optimised iteratively.
#
#   v1 (published):  y is the catalogue-tied MAGNITUDE; w is drawn from
#       four values (0.5, 1, 2, 4) x that series' folded-profile ramp
#       width; t_e from 361 nodes spanning +/- 3 cadences about the
#       ephemeris prediction; argmin on ties.
#   v2 (robustness): y is RELATIVE FLUX; w from the eight-value grid
#       V2_WIDTH_GRID_S common to every band; t_e on a 2 s grid; the
#       midpoint of a flat valley on ties.
#
# The fit is accepted when the step exceeds five times the larger of the
# residual scatter and the median error, the two points bracketing t_e are
# no further apart than 2.5 cadences, and t_e is not on the grid boundary.

def to_relative_flux(mag, mag_err) -> tuple[np.ndarray, np.ndarray, float]:
    """Magnitudes to flux relative to the window's own brightest decile.

    Returns ``(flux, flux_err, m_ref)``.  The reference level only sets the
    unit — the fitted edge time is invariant to it — but choosing the
    bright level makes the fitted step a fractional depth, which is the
    quantity a reader can check against a light curve.
    """
    m = np.asarray(mag, dtype=float)
    e = np.asarray(mag_err, dtype=float)
    m_ref = float(np.nanpercentile(m, 10.0)) if m.size else float("nan")
    f = 10.0 ** (-0.4 * (m - m_ref))
    return f, 0.4 * math.log(10.0) * f * e, m_ref


def edge_window(t, y, e, t_guess_d: float, period_d: float,
                half_phase: float = p3.EDGE_WINDOW_PHASE
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The points one edge is fitted to: the window and the outlier clip.

    Exactly the selection ``run_cv_phase3.cmd_edges`` applies — points
    within ``half_phase`` cycles of the predicted edge, then a six-MAD clip
    about the window median with a 0.05 mag floor on the MAD — lifted here
    so the injection test, the bootstrap and the v2 refit all see the same
    points the published fit saw.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    e = np.asarray(e, dtype=float)
    sel = np.abs(t - t_guess_d) <= half_phase * period_d
    tt, yy, ee = t[sel], y[sel], e[sel]
    if tt.size == 0:
        return tt, yy, ee
    med = float(np.median(yy))
    mad = 1.4826 * float(np.median(np.abs(yy - med)))
    keep = np.abs(yy - med) < 6.0 * max(mad, 0.05)
    return tt[keep], yy[keep], ee[keep]


def v2_time_grid(t_guess_d: float, cadence_s: float, period_d: float,
                 step_s: float = V2_TIME_STEP_S) -> np.ndarray:
    """v2 search grid for the edge epoch, centred on the prediction.

    Half-span: three cadences, as in v1, but never beyond the fitting
    window itself — an edge placed outside the window is not constrained
    by any point in it.  Odd node count so the prediction is a node.
    """
    half = min(3.0 * float(cadence_s), p3.EDGE_WINDOW_PHASE * period_d
               * 86400.0)
    n_half = max(1, int(math.floor(half / step_s)))
    return t_guess_d + step_s / 86400.0 * np.arange(-n_half, n_half + 1)


def fit_edge_profile(times_d, y, dy, t_grid_d, width_grid_d,
                     median_cadence_s: float,
                     min_points: int = p3.EDGE_MIN_POINTS,
                     min_snr: float = p3.EDGE_MIN_SNR,
                     max_bracket_cadence: float = p3.EDGE_MAX_BRACKET_CADENCE,
                     tie: str = "centre") -> dict:
    """The four-parameter ramp fit, vectorised, with the valley disclosed.

    Same model, same chi-squared surface and same acceptance rules as
    :func:`macro_phot.phase3.fit_edge` with free levels (the unit tests
    hold the two surfaces equal to rounding).  Two differences, both
    deliberate:

    * the levels are solved in closed form for every grid node at once, so
      a fit costs a millisecond instead of fifty and a 400-replicate
      bootstrap of every edge is affordable;
    * ``tie='centre'`` returns the MIDPOINT of the contiguous run of time
      nodes whose profiled chi-squared equals the minimum to
      :data:`VALLEY_RTOL`.  When the best ramp lies wholly between two
      exposures chi-squared is flat across the gap and the data say only
      "somewhere in here"; the midpoint is the symmetric estimate, and
      ``valley_s`` records how wide the indifference was.  ``tie='first'``
      returns the first node of the valley instead.

    Returns a dict: ``t_edge_d``, ``width_d``, ``level_bright``,
    ``level_faint``, ``step``, ``snr``, ``chi2``, ``dof``, ``chi2nu``,
    ``n_points``, ``bracket_s``, ``valley_s``, ``grid_step_s``,
    ``accepted``, ``reason``.  ``chi2`` and ``dof`` are both returned
    because a reduced chi-squared without its degrees of freedom cannot be
    judged: these fits have two to eight.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(dy, dtype=float)
    ok = np.isfinite(t) & np.isfinite(yv) & np.isfinite(e) & (e > 0)
    t, yv, e = t[ok], yv[ok], e[ok]
    nan = float("nan")
    tg = np.asarray(t_grid_d, dtype=float)
    wg = np.asarray(width_grid_d, dtype=float)
    out = {"t_edge_d": nan, "width_d": nan, "level_bright": nan,
           "level_faint": nan, "step": nan, "snr": nan, "chi2": nan,
           "dof": None, "chi2nu": nan, "n_points": int(t.size),
           "bracket_s": nan, "valley_s": nan,
           "grid_step_s": (float((tg[1] - tg[0]) * 86400.0)
                           if tg.size > 1 else nan),
           "accepted": False, "reason": ""}
    if t.size < min_points:
        out["reason"] = (f"only {t.size} usable points (need {min_points})")
        return out
    w = 1.0 / e ** 2
    sw = float(w.sum())
    sy = float(np.sum(w * yv))
    syy = float(np.sum(w * yv * yv))
    chi2 = np.full((tg.size, wg.size), np.inf)
    lev_a = np.full((tg.size, wg.size), np.nan)
    lev_b = np.full((tg.size, wg.size), np.nan)
    for j, wid in enumerate(wg):
        # Ramp value of every point for every trial epoch: (n_t, n_pts).
        r = np.clip((t[None, :] - tg[:, None]) / max(float(wid), 1e-12)
                    + 0.5, 0.0, 1.0)
        sr = r @ w
        srr = (r * r) @ w
        sry = r @ (w * yv)
        det = sw * srr - sr * sr
        good = np.isfinite(det) & (np.abs(det) >= 1e-30)
        with np.errstate(divide="ignore", invalid="ignore"):
            b = (sw * sry - sr * sy) / det
            a = (sy - b * sr) / sw
            # chi2 = sum w (y - a - b r)^2, expanded so no (n_t, n_pts)
            # residual array is needed.
            c2 = (syy - 2.0 * a * sy - 2.0 * b * sry + a * a * sw
                  + 2.0 * a * b * sr + b * b * srr)
        chi2[good, j] = np.maximum(c2[good], 0.0)
        lev_a[good, j] = a[good]
        lev_b[good, j] = b[good]
    if not np.isfinite(chi2).any():
        out["reason"] = "no invertible fit on the grid"
        return out
    i0, j0 = np.unravel_index(int(np.argmin(chi2)), chi2.shape)
    chi_min = float(chi2[i0, j0])
    # The valley: contiguous time nodes, at the best width, whose chi2 is
    # the minimum to rounding.  The expanded chi2 above carries rounding of
    # order 1e-12 * syy, so the tolerance is set on syy, not on chi_min
    # (which can be arbitrarily close to zero for a clean edge).
    tol = VALLEY_RTOL * max(syy, 1.0)
    col = chi2[:, j0]
    lo = hi = int(i0)
    while lo > 0 and col[lo - 1] <= chi_min + tol:
        lo -= 1
    while hi < tg.size - 1 and col[hi + 1] <= chi_min + tol:
        hi += 1
    if tie == "centre":
        t_edge = 0.5 * float(tg[lo] + tg[hi])
        i_use = int(round(0.5 * (lo + hi)))
    elif tie == "first":
        t_edge = float(tg[lo])
        i_use = lo
    else:
        raise ValueError(f"tie must be 'centre' or 'first', got {tie!r}")
    width = float(wg[j0])
    bright = float(lev_a[i_use, j0])
    step = float(lev_b[i_use, j0])
    n_free = 4
    dof = max(int(t.size) - n_free, 1)
    r = np.clip((t - t_edge) / max(width, 1e-12) + 0.5, 0.0, 1.0)
    resid = yv - (bright + step * r)
    scatter = float(np.std(resid, ddof=1)) if t.size > 2 else nan
    noise = max(scatter if np.isfinite(scatter) else 0.0,
                float(np.median(e)))
    snr = abs(step) / noise if noise > 0 else float("inf")
    before, after = t[t <= t_edge], t[t >= t_edge]
    bracket_s = (float((after.min() - before.max()) * 86400.0)
                 if before.size and after.size else float("inf"))
    accepted, reason = True, "accepted"
    if not np.isfinite(snr) or snr < min_snr:
        accepted, reason = False, (
            f"step SNR {snr:.1f} below {min_snr:.0f} — not distinguishable "
            "from flickering")
    elif not np.isfinite(bracket_s) or (
            np.isfinite(median_cadence_s) and median_cadence_s > 0 and
            bracket_s > max_bracket_cadence * median_cadence_s):
        accepted, reason = False, (
            f"edge fell in a {bracket_s:.0f} s gap, more than "
            f"{max_bracket_cadence:g}x the {median_cadence_s:.0f} s cadence")
    elif lo == 0 or hi == tg.size - 1:
        accepted, reason = False, (
            "best epoch sits on the edge of the search grid")
    out.update({"t_edge_d": t_edge, "width_d": width,
                "level_bright": bright, "level_faint": bright + step,
                "step": step, "snr": float(snr), "chi2": chi_min,
                "dof": dof, "chi2nu": chi_min / dof,
                "bracket_s": bracket_s,
                "valley_s": float((tg[hi] - tg[lo]) * 86400.0),
                "accepted": accepted, "reason": reason,
                "_resid": resid, "_model": bright + step * r,
                "_t": t, "_e": e})
    return out


def quantisation_sigma_s(grid_step_s: float) -> float:
    """Standard deviation of a value rounded to a grid of this step.

    A uniform distribution of width ``step`` has sigma = step / sqrt(12).
    This is the term CV-R2 adds to the timing budget: the v1 time grid has
    361 nodes over six cadences, so its step is cadence / 60 — 3.65 s at
    the 219 s Mode0 cadence, tens of seconds on the sparse 2024 nights.
    """
    if grid_step_s is None or not np.isfinite(grid_step_s):
        return float("nan")
    return float(abs(grid_step_s) / math.sqrt(12.0))


# ===========================================================================
# 3.  A per-edge bootstrap
# ===========================================================================

def bootstrap_edge(fit: dict, refit: Callable[[np.ndarray], float],
                   n_boot: int = N_EDGE_BOOT, seed: int = SEED) -> dict:
    """Wild residual bootstrap of ONE edge fit.

    ``fit`` is a :func:`fit_edge_profile` record (it carries the best-fit
    model and residuals).  Each replicate rebuilds the window as

        y*_i = model_i + s_i * r_i * sqrt(n / dof)

    with ``s_i`` an independent random sign, and ``refit`` — a callable
    taking the replicate's y-array and returning the fitted edge time in
    days — measures it again.  ``refit`` should NOT apply the acceptance
    gate: the edge being bootstrapped has already been accepted, and
    discarding the replicates that wander furthest (into a gap, below the
    step-S/N bar) would truncate exactly the tail the error bar is for.

    WHY WILD, AND WHY THE INFLATION.  The residuals of these fits are not
    photon noise: chi2_nu runs from 2 to 1,900 because a ramp does not
    describe flickering.  Flipping each residual's sign keeps every
    point's own misfit attached to that point (heteroscedasticity
    preserved) while destroying the particular pattern this cycle
    happened to have.  The factor sqrt(n / dof) undoes the shrinkage a
    four-parameter fit to six-to-twelve points imposes on its own
    residuals; without it the bootstrap is over-confident by up to 1.7x.

    WHAT IT IS NOT.  Sign flips whiten: correlated flickering that drags
    several consecutive points the same way is under-represented, so this
    sigma is a LOWER bound on the per-edge error.  The stage that calls
    this function therefore also tests it — the chi-squared of real
    same-cycle band differences under these sigmas, per band pair, with
    its degrees of freedom — and the paper may quote the per-edge sigma
    only with that calibration beside it.

    Returns ``sigma_s`` (half the 16-84 percentile range), ``rms_s``,
    ``bias_s`` (median replicate minus the fit), ``n_ok``, ``n_boot``.
    """
    nan = float("nan")
    out = {"sigma_s": nan, "rms_s": nan, "bias_s": nan, "n_ok": 0,
           "n_boot": int(n_boot)}
    # ``accepted`` is deliberately not required: the caller decides which
    # edges are bootstrapped (the published set), and one published edge
    # sits in a flat valley the symmetric rule would have rejected.
    if "_resid" not in fit:
        return out
    resid = np.asarray(fit["_resid"], dtype=float)
    model = np.asarray(fit["_model"], dtype=float)
    n = resid.size
    dof = max(int(fit.get("dof") or 1), 1)
    infl = math.sqrt(n / dof)
    rng = np.random.default_rng(seed)
    got = []
    for _ in range(int(n_boot)):
        signs = rng.choice((-1.0, 1.0), size=n)
        te = refit(model + signs * resid * infl)
        if te is not None and np.isfinite(te):
            got.append((te - fit["t_edge_d"]) * 86400.0)
    out["n_ok"] = len(got)
    if len(got) >= 20:
        g = np.asarray(got)
        lo, hi = np.percentile(g, [16.0, 84.0])
        out["sigma_s"] = float((hi - lo) / 2.0)
        out["rms_s"] = float(g.std(ddof=1))
        out["bias_s"] = float(np.median(g))
    return out


# ===========================================================================
# 4.  Injection with the real estimator
# ===========================================================================

def flux_edge_template(times_d, t_edge_d: float, width_d: float,
                       depth_mag: float) -> np.ndarray:
    """An egress in RELATIVE FLUX: bright level 1, linear fall to the faint
    level ``10^(-0.4 depth_mag)`` over ``width_d`` centred on ``t_edge_d``.

    The ramp is linear in flux because that is what an emitting region
    passing behind a limb produces to first order — the visible fraction
    of the region falls, and flux is linear in visible fraction.  The SAME
    ``t_edge_d`` and ``width_d`` in every band is the achromatic null: the
    region's visibility does not depend on wavelength, only how bright it
    is relative to everything else (``depth_mag``) does.
    """
    faint = 10.0 ** (-0.4 * float(depth_mag))
    return 1.0 - (1.0 - faint) * p3.ramp(times_d, t_edge_d, width_d)


def residual_pool(times_d, mag, period_d: float, epoch_d: float,
                  n_bins: int = 20, edge_exclude_phase: float = 0.06
                  ) -> np.ndarray:
    """Real relative-flux residuals of one night about its own fold.

    The noise an injected edge is recovered against must be the star's own:
    flickering is correlated on the timescale of the fit window and a
    Gaussian draw from the photometric errors is a factor of several too
    kind (the real fits' chi2_nu says so).  This returns the night's light
    curve, in flux relative to its bright level, minus a ``n_bins``-bin
    phase-folded median profile, IN TIME ORDER, so a contiguous block of it
    carries the measured short-range correlation.

    Points within ``edge_exclude_phase`` cycles of the profile's two
    steepest gradients are dropped: there the residual is dominated by the
    bin being wider than the edge, which is a property of the profile and
    not of the star, and leaving it in would inject edge-shaped "noise".
    """
    t = np.asarray(times_d, dtype=float)
    m = np.asarray(mag, dtype=float)
    ok = np.isfinite(t) & np.isfinite(m)
    t, m = t[ok], m[ok]
    if t.size < 3 * n_bins:
        return np.array([])
    f = 10.0 ** (-0.4 * (m - np.percentile(m, 10.0)))
    ph = p3.phase_of(t, period_d, epoch_d)
    idx = np.clip((ph * n_bins).astype(int), 0, n_bins - 1)
    prof = np.array([np.median(f[idx == b]) if np.any(idx == b) else np.nan
                     for b in range(n_bins)])
    if not np.isfinite(prof).all():
        # Fill unsampled bins from their neighbours (circularly); a night
        # that leaves a bin empty still has to yield residuals elsewhere.
        good = np.flatnonzero(np.isfinite(prof))
        if good.size < n_bins // 2:
            return np.array([])
        centres = (np.arange(n_bins) + 0.5) / n_bins
        prof = np.interp(centres, np.concatenate([centres[good] - 1,
                                                  centres[good],
                                                  centres[good] + 1]),
                         np.tile(prof[good], 3))
    # Linear interpolation of the profile between bin centres (circular).
    centres = (np.arange(n_bins) + 0.5) / n_bins
    model = np.interp(ph, np.concatenate([[centres[-1] - 1.0], centres,
                                          [centres[0] + 1.0]]),
                      np.concatenate([[prof[-1]], prof, [prof[0]]]))
    resid = f - model
    # The two steepest profile gradients mark the rise and the fall.
    grad = np.abs(np.diff(np.concatenate([prof, prof[:1]])))
    steep = (np.argsort(grad)[-2:] + 1.0) / n_bins     # bin-boundary phases
    keep = np.ones(t.size, dtype=bool)
    for s in steep:
        dphi = np.abs(np.mod(ph - s + 0.5, 1.0) - 0.5)
        keep &= dphi > edge_exclude_phase
    return resid[keep]


def inject_recover(bands: dict, period_d: float, true_widths_s,
                   fit_all: Callable, estimator_names: Sequence[str],
                   n_real: int = N_INJECT,
                   seed: int = SEED, noise: str = "rolled",
                   guess_scatter_s: float = GUESS_SCATTER_S) -> dict:
    """Inject one ACHROMATIC edge into every band of a night; recover each.

    ``bands`` maps a band slot to a dict with ``t`` (BJD of that band's
    real exposures on the night), ``e`` (their real magnitude errors,
    already inflated), ``depth_mag`` (that band's real edge depth) and
    ``pool`` (its real residual pool, relative flux; may be empty).

    ``fit_all`` is a callable
    ``f(band, t, mag, mag_err, t_guess_d) -> {estimator name: t_edge_d}``
    that runs the REAL fits for that band on one synthetic night — the v1
    one with that band's own width grid, and the re-timing estimators —
    and returns NaN for an estimator whose acceptance gate fails.  One
    callable for all estimators, rather than one each, because the real
    procedure is sequential: an edge EXISTS when the published v1 gate
    accepts it, and the other estimators then re-time that same edge.  The
    injection has to reproduce that selection or its recovered fraction and
    its bias are those of a procedure nobody ran.  ``estimator_names``
    lists the keys to tabulate.

    For each of ``n_real`` realisations one true edge time is drawn
    uniformly over the night, and for each true width the same edge is
    written into every band's timestamps (common random numbers across
    widths and bands, so differences between cells are differences of
    estimator and sampling, not of luck).  Noise is a contiguous block of
    the band's real residual pool at a random offset (``noise='rolled'``)
    or Gaussian at the quoted errors (``noise='gaussian'``, the optimistic
    case, kept as a cross-check).

    Returns ``{"per_band": [...], "per_pair": [...]}``:

    ``per_band`` rows carry, per (estimator, band, true width): ``n_try``,
    ``n_ok``, the SIGNED ``bias_mean_s`` with its standard error
    ``bias_se_s``, ``bias_median_s``, ``sigma_s`` (16-84 half range) and
    ``rms_s``.

    ``per_pair`` rows carry, per (estimator, band a, band b, true width in
    a, true width in b): the SIGNED mean of ``(t_a - t_b)`` over the
    realisations both bands recovered, and its standard error.  Because
    the truth is the same instant in both bands, this is the differential
    bias of the band-offset measurement itself — including the interleaved
    filter sequence, the band-dependent depth and the band-dependent width
    grid.  Unequal true widths (a != b) are the "per-band ramp widths"
    case: same midpoint, different shape.
    """
    if noise not in ("rolled", "gaussian"):
        raise ValueError(f"noise must be 'rolled' or 'gaussian', got {noise!r}")
    rng = np.random.default_rng(seed)
    widths = [float(x) for x in true_widths_s]
    names = sorted(bands)
    t_lo = max(float(np.min(bands[b]["t"])) for b in names)
    t_hi = min(float(np.max(bands[b]["t"])) for b in names)
    # recovered[est][band][width] -> array over realisations (NaN = miss)
    rec = {est: {b: {wd: np.full(n_real, np.nan) for wd in widths}
                 for b in names} for est in estimator_names}
    tried = {b: np.zeros(n_real, dtype=bool) for b in names}
    half_win = p3.EDGE_WINDOW_PHASE * period_d
    for k in range(int(n_real)):
        if t_hi - t_lo <= 2 * half_win:
            # A night shorter than one window: put the edge mid-night.
            t_true = 0.5 * (t_lo + t_hi) + rng.uniform(-0.25, 0.25) * half_win
        else:
            t_true = float(rng.uniform(t_lo + half_win, t_hi - half_win))
        t_guess = t_true + rng.normal(0.0, guess_scatter_s) / 86400.0
        for b in names:
            t = np.asarray(bands[b]["t"], dtype=float)
            e = np.asarray(bands[b]["e"], dtype=float)
            sel = np.abs(t - t_guess) <= half_win
            if sel.sum() < p3.EDGE_MIN_POINTS:
                # Draw this band's noise variates anyway?  No: nothing is
                # drawn for a band with no window, and because the draws
                # below are per band from independent child streams the
                # other bands' numbers do not depend on it.
                continue
            tried[b][k] = True
            pool = np.asarray(bands[b].get("pool", []), dtype=float)
            child = np.random.default_rng([seed, k, names.index(b)])
            n_pts = int(t.size)
            if noise == "rolled" and pool.size >= max(20, int(sel.sum())):
                off = int(child.integers(0, pool.size))
                reps = int(math.ceil(n_pts / pool.size)) + 1
                eps = np.tile(np.roll(pool, -off), reps)[:n_pts]
                # A random sign for the whole block: the pool's own mean
                # trend must not bias every realisation the same way.
                eps = eps * (1.0 if child.random() < 0.5 else -1.0)
            else:
                # Gaussian in relative flux at the quoted magnitude errors.
                eps = child.normal(0.0, 1.0, n_pts) * (0.4 * math.log(10.0)
                                                       * e)
            for wd in widths:
                flux = flux_edge_template(t, t_true, wd / 86400.0,
                                          bands[b]["depth_mag"]) + eps
                # A magnitude needs positive flux.  At these depths and
                # noise levels a non-positive draw is rare; it is clipped
                # to a small positive floor rather than dropped, because a
                # dropped point changes the sampling the test is about.
                flux = np.maximum(flux, 1e-3)
                mag = -2.5 * np.log10(flux) + 16.0
                got = fit_all(b, t, mag, e, t_guess)
                for est in estimator_names:
                    te = got.get(est)
                    if te is not None and np.isfinite(te):
                        rec[est][b][wd][k] = (te - t_true) * 86400.0
    per_band, per_pair = [], []
    for est in estimator_names:
        for b in names:
            for wd in widths:
                d = rec[est][b][wd]
                good = d[np.isfinite(d)]
                row = {"estimator": est, "band": b, "true_width_s": wd,
                       "n_try": int(tried[b].sum()), "n_ok": int(good.size),
                       "bias_mean_s": None, "bias_se_s": None,
                       "bias_median_s": None, "sigma_s": None,
                       "rms_s": None}
                if good.size >= 10:
                    lo, hi = np.percentile(good, [16.0, 84.0])
                    row.update({
                        "bias_mean_s": float(good.mean()),
                        "bias_se_s": float(good.std(ddof=1)
                                           / math.sqrt(good.size)),
                        "bias_median_s": float(np.median(good)),
                        "sigma_s": float((hi - lo) / 2.0),
                        "rms_s": float(good.std(ddof=1))})
                per_band.append(row)
        for a, b in BAND_PAIRS:
            if a not in rec[est] or b not in rec[est]:
                continue
            for wa in widths:
                for wb in widths:
                    d = rec[est][a][wa] - rec[est][b][wb]
                    good = d[np.isfinite(d)]
                    row = {"estimator": est, "band_a": a, "band_b": b,
                           "true_width_a_s": wa, "true_width_b_s": wb,
                           "n_ok": int(good.size), "dbias_mean_s": None,
                           "dbias_se_s": None, "dbias_sd_s": None}
                    if good.size >= 10:
                        row.update({
                            "dbias_mean_s": float(good.mean()),
                            "dbias_se_s": float(good.std(ddof=1)
                                                / math.sqrt(good.size)),
                            "dbias_sd_s": float(good.std(ddof=1))})
                    per_pair.append(row)
    return {"per_band": per_band, "per_pair": per_pair}


def magnitude_midpoint_shift(depth_mag: float) -> float:
    """Where, along a LINEAR FLUX ramp, the MAGNITUDE midpoint falls.

    Returns the fraction x in (0.5, 1) of the ramp at which the magnitude
    is half-way between the bright and faint levels.  A flux falling
    linearly from 1 to q = 10^(-0.4 dm) has magnitude midpoint where the
    flux equals the geometric mean sqrt(q), i.e. at

        x = (1 - sqrt(q)) / (1 - q) = 1 / (1 + sqrt(q)).

    For a shallow edge x -> 0.5; for a deep one x -> 1.  A magnitude-space
    ramp fit therefore times a deeper band's edge LATER than a shallower
    band's even when the two are the same event, by (x_a - x_b) of the
    ramp width.  This is the analytic part of the estimator's band bias;
    the injection test measures the rest.
    """
    q = 10.0 ** (-0.4 * float(depth_mag))
    return float(1.0 / (1.0 + math.sqrt(q)))


# ===========================================================================
# 5.  The O-C model with its nuisance terms
# ===========================================================================

def night_level_epochs(rows: Sequence[dict], band_constants: dict,
                       sigma_key: str = "sigma_s") -> list[dict]:
    """Collapse per-night per-band epochs to ONE epoch per night.

    ``rows`` carry ``night``, ``band``, ``cycle``, ``oc_s`` and
    ``sigma_key``.  Each band's constant (its mean edge-time offset, fitted
    elsewhere) is subtracted first, so combining bands does not re-import
    the band offset as scatter; the night's epoch is then the
    inverse-variance weighted mean of what is left.

    The reason this exists: 36 per-night per-band epochs are not 36
    independent data.  Three bands of one night watch the same accretion
    spot through the same cycles, so their residuals move together, and a
    fit that counts them separately claims three times the nights it has.
    Seventeen nights are seventeen measurements.

    Each output row carries ``n_bands``, ``n_cycles``, ``cycle`` (weighted
    mean), ``oc_s``, ``sigma_s`` (the propagated weighted-mean error — the
    fit that consumes these rescales it by the scatter of the nights
    themselves) and ``within_rms_s``, the scatter of the bands about the
    night mean where there are at least two.
    """
    groups: dict[str, list] = {}
    for r in rows:
        groups.setdefault(str(r["night"]), []).append(r)
    out = []
    for night, g in sorted(groups.items()):
        oc = np.array([float(r["oc_s"]) - float(band_constants.get(
            r["band"], 0.0)) for r in g])
        s = np.array([float(r[sigma_key]) for r in g])
        cyc = np.array([float(r["cycle"]) for r in g])
        w = 1.0 / s ** 2
        mean = float(np.sum(w * oc) / np.sum(w))
        out.append({
            "night": night, "n_bands": len(g),
            "n_cycles": int(sum(int(r.get("n_cycles") or 1) for r in g)),
            "bands": "".join(sorted(str(r["band"]) for r in g)),
            "era": g[0].get("era"),
            "cycle": float(np.sum(w * cyc) / np.sum(w)),
            "oc_s": mean,
            "sigma_s": float(math.sqrt(1.0 / np.sum(w))),
            "within_rms_s": (float(np.std(oc, ddof=1)) if len(g) > 1
                             else None),
        })
    return out


def fit_oc_model(cycle, oc_s, sigma_s, band=None, era=None,
                 quadratic: bool = True) -> dict:
    """Weighted least squares of an O-C curve with explicit nuisance terms.

    The model for an epoch at cycle E, in band b, taken in era k, is

        O-C = c_b  +  d_k  +  beta * (E - E0)  +  gamma * (E - E0)^2

    * ``c_b`` — one constant PER BAND when ``band`` is given (one overall
      constant otherwise).  This is the term whose absence let a ~100 s
      band offset masquerade as scatter and, because the band mix changes
      with time, leak into the period and the quadratic.
    * ``d_k`` — an offset for every era after the first when ``era`` is
      given: the nuisance term for "the camera, the filters and the error
      model all changed between 2024 and 2025".  It is degenerate with
      nothing as long as each era spans a range of cycles, but it removes
      the between-era lever arm from gamma — which is the point: with it,
      the quadratic is constrained only by curvature WITHIN eras.
    * ``beta`` — a period correction; ``gamma`` — the period derivative
      (``quadratic=False`` drops it).

    Returns a dict with the coefficients and their covariance, and for the
    quadratic term BOTH errors the committee asked for:

    ``gamma_sigma_budget``  propagated from the supplied ``sigma_s`` with
                            no rescaling in either direction;
    ``gamma_sigma_scatter`` the same multiplied by sqrt(chi2_nu), i.e. the
                            error the residuals' own scatter implies.

    ``chi2``, ``dof`` and ``chi2nu`` are all returned.  A band (or era)
    that contributes a single epoch is fitted exactly by its own constant
    and adds nothing to chi-squared or to gamma; ``n_informative`` counts
    the epochs that are not absorbed that way, and ``singleton_levels``
    names the ones that are.
    """
    c = np.asarray(cycle, dtype=float)
    y = np.asarray(oc_s, dtype=float)
    s = np.asarray(sigma_s, dtype=float)
    n = c.size
    ok = np.isfinite(c) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    if band is not None:
        band = np.asarray([str(b) for b in band])
    if era is not None:
        era = np.asarray([str(k) for k in era])
    if not ok.all():
        c, y, s = c[ok], y[ok], s[ok]
        band = band[ok] if band is not None else None
        era = era[ok] if era is not None else None
        n = c.size
    e0 = float(c.mean()) if n else 0.0
    x = c - e0
    cols, names = [], []
    singletons = []
    if band is not None:
        for b in sorted(set(band)):
            col = (band == b).astype(float)
            cols.append(col)
            names.append(f"band:{b}")
            if col.sum() == 1:
                singletons.append(f"band:{b}")
    else:
        cols.append(np.ones(n))
        names.append("const")
    if era is not None:
        # An era whose epochs are exactly one band's epochs is the same
        # column as that band's constant: it is set aside first, and said
        # so.  Of the eras that remain, the first is the reference level
        # (its offset is zero by definition) and the others get a column.
        free = []
        for k in sorted(set(era)):
            col = (era == k).astype(float)
            if any(np.array_equal(col, cc) for cc in cols):
                singletons.append(f"era:{k} (aliased with a band constant)")
            else:
                free.append((k, col))
        for k, col in free[1:]:
            cols.append(col)
            names.append(f"era:{k}")
    # The cycle axis is scaled to unit half-range before the solve.  In raw
    # cycles the quadratic column is ~1e7 times the constant ones and the
    # normal matrix's condition number reaches 1e15, at which point a
    # perfectly well-posed fit is indistinguishable from a singular one.
    # The coefficients and their covariance are scaled back afterwards.
    xs = float(np.max(np.abs(x))) if n and np.max(np.abs(x)) > 0 else 1.0
    scale = [1.0] * len(cols)
    cols.append(x / xs)
    names.append("beta")
    scale.append(1.0 / xs)
    if quadratic:
        cols.append((x / xs) ** 2)
        names.append("gamma")
        scale.append(1.0 / xs ** 2)
    a = np.column_stack(cols)
    p = a.shape[1]
    scale = np.asarray(scale)
    nan = float("nan")
    out = {"names": names, "n": int(n), "n_params": int(p), "e0": e0,
           "coef": None, "cov": None, "chi2": nan, "dof": None,
           "chi2nu": nan, "gamma": nan, "gamma_sigma_budget": nan,
           "gamma_sigma_scatter": nan, "beta": nan, "beta_sigma_budget": nan,
           "resid": None, "singleton_levels": singletons,
           "n_informative": int(n - sum(1 for z in singletons
                                        if z.startswith("band:"))),
           "rank_deficient": False}
    if n <= p:
        out["rank_deficient"] = True
        return out
    w = 1.0 / s ** 2
    atw = a.T * w
    m = atw @ a
    if np.linalg.matrix_rank(m) < p or np.linalg.cond(m) > 1e12:
        out["rank_deficient"] = True
        return out
    cov = np.linalg.inv(m)
    coef = cov @ (atw @ y)
    resid = y - a @ coef
    # Back to raw cycles: coef -> coef * scale, cov -> D cov D.
    coef = coef * scale
    cov = cov * np.outer(scale, scale)
    chi2 = float(np.sum(w * resid ** 2))
    dof = int(n - p)
    chi2nu = chi2 / dof
    out.update({"coef": coef, "cov": cov, "chi2": chi2, "dof": dof,
                "chi2nu": chi2nu, "resid": resid,
                "wrms_s": float(math.sqrt(np.sum(w * resid ** 2)
                                          / np.sum(w))),
                "rms_s": float(math.sqrt(np.mean(resid ** 2)))})
    ib = names.index("beta")
    out["beta"] = float(coef[ib])
    out["beta_sigma_budget"] = float(math.sqrt(cov[ib, ib]))
    out["beta_sigma_scatter"] = out["beta_sigma_budget"] * math.sqrt(chi2nu)
    if quadratic:
        ig = names.index("gamma")
        out["gamma"] = float(coef[ig])
        out["gamma_sigma_budget"] = float(math.sqrt(cov[ig, ig]))
        out["gamma_sigma_scatter"] = (out["gamma_sigma_budget"]
                                      * math.sqrt(chi2nu))
    return out


def pdot_from_quadratic(gamma_s_per_cycle2: float, period_d: float) -> float:
    """dP/dt (dimensionless) from the quadratic O-C coefficient.

    A steady period derivative puts ``0.5 * Pdot * P * E^2`` into an O-C
    curve (in time units, with P in the same units), so
    ``Pdot = 2 * gamma / P`` with gamma in seconds per cycle squared and P
    in seconds.
    """
    return float(2.0 * gamma_s_per_cycle2 / (period_d * 86400.0))


def per_group_chi2(resid, sigma, groups) -> list[dict]:
    """Chi-squared of fit residuals per group, each with its own count.

    Standing rule 1: a global chi2_nu of 0.9 can be the average of one
    band at 1.3 and another at 0.3, and was.  This returns, per group,
    ``n``, ``chi2``, ``chi2_per_n`` and the residual ``rms`` beside the
    median assigned sigma, so over- and under-stated errors are both
    visible.  ``chi2_per_n`` divides by the number of epochs in the group,
    not by a per-group dof — the fitted parameters are shared between
    groups and cannot be apportioned — and the caller says so.
    """
    r = np.asarray(resid, dtype=float)
    s = np.asarray(sigma, dtype=float)
    g = np.asarray([str(x) for x in groups])
    out = []
    for key in sorted(set(g)):
        sel = g == key
        chi2 = float(np.sum((r[sel] / s[sel]) ** 2))
        out.append({"group": key, "n": int(sel.sum()), "chi2": chi2,
                    "chi2_per_n": chi2 / int(sel.sum()),
                    "rms_s": float(math.sqrt(np.mean(r[sel] ** 2))),
                    "sigma_median_s": float(np.median(s[sel]))})
    return out


def scatter_sigma_by_band(resid, band, n_cycles, dof_total: int) -> dict:
    """Per-band single-cycle scatter implied by the fit residuals.

    Model: an epoch averaging ``n`` cycles in band b has variance
    ``s_b^2 / n``.  Then ``s_b^2 = mean(n * resid^2)`` over that band's
    epochs, corrected for the fitted parameters by the global factor
    ``N / dof_total`` (the parameters are shared, so the correction is
    too).  This is the scatter-based replacement for the transported
    injection budget: it is measured on the epochs it is applied to.
    """
    r = np.asarray(resid, dtype=float)
    b = np.asarray([str(x) for x in band])
    n = np.asarray(n_cycles, dtype=float)
    corr = r.size / max(int(dof_total), 1)
    out = {}
    for key in sorted(set(b)):
        sel = b == key
        if sel.sum() >= 2:
            out[key] = float(math.sqrt(np.mean(n[sel] * r[sel] ** 2) * corr))
    return out


# ===========================================================================
# 6.  Spot longitude and the physical scales of a period derivative
# ===========================================================================

def seconds_to_degrees(seconds, period_d: float):
    """Edge-time offset in seconds to accretion-spot longitude in degrees.

    The bright-phase edge is the accretion region rotating behind the
    white dwarf's limb, so a shift of the edge by dt is a shift of the
    region in longitude by ``360 * dt / P`` (for a synchronised system the
    spin period is the orbital period).  Positive = later = trailing.
    """
    return np.asarray(seconds, dtype=float) * 360.0 / (period_d * 86400.0)


def gr_orbital_pdot(period_d: float, m1_msun: float = GR_M1_MSUN,
                    m2_msun: float = GR_M2_MSUN) -> float:
    """|dP/dt| of two point masses losing angular momentum to gravitational
    radiation, at this period, with no mass transfer.

        Pdot = -(192 pi / 5) (2 pi G Mc / (c^3 P))^(5/3),
        Mc = (m1 m2)^(3/5) / (m1 + m2)^(1/5)   (the chirp mass)

    Returned as a positive number.  It is the FLOOR of what secular
    orbital evolution does below the period gap — mass transfer changes
    the sign and the size by factors of order unity, not by orders of
    magnitude — and exists so the timing bound has a physical scale beside
    it (standing rule 2).
    """
    p = float(period_d) * 86400.0
    m1, m2 = m1_msun * _MSUN, m2_msun * _MSUN
    mc = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2
    x = 2.0 * math.pi * _G * mc / (_C ** 3 * p)
    return float(192.0 * math.pi / 5.0 * x ** (5.0 / 3.0))


def secular_pdot(period_d: float, tau_gyr: float) -> float:
    """P / tau: the period derivative of an orbit evolving on ``tau``."""
    return float(period_d / (tau_gyr * 1.0e9 * 365.25))


def libration_pdot(period_d: float, amp_deg: float = LIBRATION_AMP_DEG,
                   period_yr: float = LIBRATION_PERIOD_YR) -> float:
    """Apparent |dP/dt| at the peak curvature of a sinusoidal libration.

    A spot whose longitude librates as ``A sin(2 pi t / T)`` moves the
    edge by ``(A/360) P sin(2 pi t / T)`` in time, so the O-C curve has
    peak curvature ``(A/360) P (2 pi / T)^2``.  A steady period derivative
    gives ``O-C(t) = 0.5 (Pdot / P) t^2``, curvature ``Pdot / P``.
    Equating the two: ``Pdot = (A/360) (2 pi / T)^2 P^2`` with P and T in
    the same units.  This is the formula ``literature_scales.py`` uses,
    kept identical so the two packages quote one number.
    """
    t_d = period_yr * 365.25
    return float(amp_deg / 360.0 * (2.0 * math.pi / t_d) ** 2 * period_d ** 2)


def historical_edge_shift_s(period_d: float) -> float:
    """How far one bright-phase edge moved between 1982 and 1985, seconds:
    half the change in bright-phase duration (cropper1986), if the change
    was symmetric about the bright-phase centre."""
    return float(0.5 * (BRIGHT_DURATION_1985 - BRIGHT_DURATION_1982)
                 * period_d * 86400.0)


def predicted_chromatic_shift_s(period_d: float, band_a: str,
                                band_b: str) -> float:
    """PREDICTED ``t_a - t_b`` of the bright-phase edge, from the literature.

    Bailey et al. (1985) measured the bright phase 0.02 cycle longer in J
    than in white light; half of that is one edge.  Scaling by
    ``ln(lambda_b / lambda_a) / ln(lambda_J / lambda_white)`` interpolates
    it to an optical band pair.  ORDER OF MAGNITUDE ONLY — a log-linear
    interpolation of one near-infrared measurement to the optical — but its
    SIGN is not in doubt: the redder band ends later, so for a bluer band
    ``a`` the prediction is negative.  Same arithmetic as
    ``literature_scales.py::chromatic_g_i_s``.
    """
    edge = 0.5 * (BRIGHT_DURATION_J - BRIGHT_DURATION_WHITE) \
        * period_d * 86400.0
    la, lb = BAND_WAVELENGTH_A[band_a], BAND_WAVELENGTH_A[band_b]
    return float(-edge * math.log(lb / la)
                 / math.log(WAVELENGTH_J_UM / WAVELENGTH_WHITE_UM))


def phase_drift_cycles(elapsed_d: float, period_d: float,
                       sigma_period_d: float) -> float:
    """Phase uncertainty accumulated over ``elapsed_d`` by a period error:
    ``(elapsed / P) * (sigma_P / P)`` cycles."""
    return float(abs(elapsed_d) / period_d * sigma_period_d / period_d)


def cyclotron_harmonic(wavelength_a: float, field_mg: float) -> float:
    """Cyclotron harmonic number observed at a wavelength, for a field.

    The cyclotron fundamental is at lambda_c = 10,710 A * (100 MG / B) in
    the non-relativistic limit (the cold-plasma value; the thermal shift
    of the harmonics at 10 keV is ~10 per cent and is ignored here because
    only the ORDER of the harmonic is used).  The harmonic number is then
    lambda_c / lambda.
    """
    lam_c = 10710.0 * (100.0 / float(field_mg))
    return float(lam_c / float(wavelength_a))


def slope_permutation_p(x, y, n_perm: int = 20_000, seed: int = SEED
                        ) -> tuple[float, float, float]:
    """Least-squares slope of y on x and its permutation p-value.

    Returns ``(slope, intercept, p)``.  The null is "y does not depend on
    x": permuting y against x destroys any dependence and keeps both
    marginal distributions.  Used for "does the spot longitude move with
    accretion state", where the samples are a dozen nights and no error
    model for the scatter is worth assuming.
    """
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    ok = np.isfinite(xv) & np.isfinite(yv)
    xv, yv = xv[ok], yv[ok]
    if xv.size < 4 or np.ptp(xv) == 0:
        return float("nan"), float("nan"), float("nan")
    xc = xv - xv.mean()
    sxx = float(np.sum(xc * xc))
    slope = float(np.sum(xc * yv) / sxx)
    icpt = float(yv.mean() - slope * xv.mean())
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(int(n_perm)):
        s = float(np.sum(xc * rng.permutation(yv)) / sxx)
        if abs(s) >= abs(slope) - 1e-15:
            count += 1
    return slope, icpt, float((count + 1) / (n_perm + 1))


# ===========================================================================
# 7.  The colour curve as a measurement
# ===========================================================================

def interpolated_colour(t_a, m_a, t_b, m_b, max_gap_s: float
                        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Colour ``a - b`` at band-b times, band a INTERPOLATED to them.

    A pairing window, however narrow, pairs each band-b exposure with the
    band-a exposure taken a fixed fraction of a filter cycle earlier (the
    wheel always turns the same way), so on a steep part of the light
    curve every pair carries the same sign of time offset and the colour
    is biased by slope x offset.  Linear interpolation of band a between
    the exposures immediately before and after the band-b time removes
    that first-order term: the interpolated magnitude is the band-a
    brightness AT the band-b instant if band a varies linearly across its
    own cadence.

    A band-b point is used only when its two band-a neighbours are no
    more than ``max_gap_s`` apart, so nothing is interpolated across a gap
    or between nights.

    Returns ``(t, colour, gap_s)`` with ``gap_s`` the bracketing spacing.
    """
    ta = np.asarray(t_a, dtype=float)
    ma = np.asarray(m_a, dtype=float)
    tb = np.asarray(t_b, dtype=float)
    mb = np.asarray(m_b, dtype=float)
    if ta.size < 2 or tb.size == 0:
        return np.array([]), np.array([]), np.array([])
    order = np.argsort(ta)
    ta, ma = ta[order], ma[order]
    pos = np.searchsorted(ta, tb)
    ok = (pos > 0) & (pos < ta.size)
    pos = np.clip(pos, 1, ta.size - 1)
    t0, t1 = ta[pos - 1], ta[pos]
    gap = (t1 - t0) * 86400.0
    ok &= (gap > 0) & (gap <= float(max_gap_s))
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (tb - t0) / (t1 - t0)
    m_interp = ma[pos - 1] + frac * (ma[pos] - ma[pos - 1])
    return tb[ok], (m_interp - mb)[ok], gap[ok]


def binned_curve(phase, value, n_bins: int = COLOUR_BINS,
                 min_count: int = COLOUR_MIN_PER_BIN
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median of ``value`` in ``n_bins`` equal phase bins.

    Returns ``(centres, medians, counts)``; a bin with fewer than
    ``min_count`` points is NaN, because an under-filled bin is a
    statement about coverage and must not be read as a colour.
    """
    ph = np.asarray(phase, dtype=float)
    v = np.asarray(value, dtype=float)
    ok = np.isfinite(ph) & np.isfinite(v)
    ph, v = ph[ok], v[ok]
    centres = (np.arange(n_bins) + 0.5) / n_bins
    med = np.full(n_bins, np.nan)
    cnt = np.zeros(n_bins, dtype=int)
    if ph.size:
        idx = np.clip((ph * n_bins).astype(int), 0, n_bins - 1)
        for b in range(n_bins):
            sel = v[idx == b]
            cnt[b] = sel.size
            if sel.size >= min_count:
                med[b] = float(np.median(sel))
    return centres, med, cnt


def _curve_stats(centres, med) -> Optional[dict]:
    """Amplitude and extremum phases of one binned curve, or None.

    The amplitude is the difference between the mean of the three reddest
    (largest) and the three bluest (smallest) bin medians — not max minus
    min, which on twenty noisy bins is biased high by selecting the two
    most extreme noise draws.  The extremum PHASES are the circular means
    of those same three bins.
    """
    good = np.isfinite(med)
    if good.sum() < 10:
        return None
    c, m = centres[good], med[good]
    order = np.argsort(m)
    lo, hi = order[:3], order[-3:]

    def cmean(ph):
        z = np.exp(2j * np.pi * ph).mean()
        return float(np.mod(np.angle(z) / (2 * np.pi), 1.0))

    return {"amplitude": float(m[hi].mean() - m[lo].mean()),
            "phase_red": cmean(c[hi]), "phase_blue": cmean(c[lo]),
            "mean": float(m.mean()), "n_bins_used": int(good.sum())}


def circ_diff(a: float, b: float) -> float:
    """Signed circular difference a - b in (-0.5, 0.5] cycles."""
    return float(np.mod(a - b + 0.5, 1.0) - 0.5)


def colour_curve_summary(phase, colour, night, n_bins: int = COLOUR_BINS,
                         min_count: int = COLOUR_MIN_PER_BIN,
                         n_boot: int = N_COLOUR_BOOT, seed: int = SEED
                         ) -> dict:
    """Amplitude and extremum phases of a colour curve, with errors.

    The errors come from a bootstrap over NIGHTS, not over pairs: pairs of
    one night share a zero point, an accretion state and correlated
    flickering, so resampling pairs would return an error bar for a
    sample several times larger than the one that exists.  Each replicate
    draws whole nights with replacement, re-bins and re-measures.

    Returns ``amplitude`` and ``amplitude_err`` (mag), ``phase_red`` /
    ``phase_blue`` with ``*_err`` (cycles; circular), ``n_pairs``,
    ``n_nights``, ``n_bins_used``, ``n_boot_ok`` and the binned curve.
    With fewer than three nights the errors are None — a bootstrap over
    two nights has three distinct outcomes and is not an error bar.
    """
    ph = np.asarray(phase, dtype=float)
    cv = np.asarray(colour, dtype=float)
    ng = np.asarray([str(x) for x in night])
    centres, med, cnt = binned_curve(ph, cv, n_bins, min_count)
    base = _curve_stats(centres, med)
    nights = sorted(set(ng))
    out = {"n_pairs": int(np.isfinite(cv).sum()), "n_nights": len(nights),
           "centres": centres, "median": med, "count": cnt,
           "amplitude": None, "amplitude_err": None, "phase_red": None,
           "phase_red_err": None, "phase_blue": None,
           "phase_blue_err": None, "mean": None, "n_bins_used": 0,
           "n_boot_ok": 0, "bin_err": np.full(n_bins, np.nan)}
    if base is None:
        return out
    out.update(base)
    if len(nights) < 3:
        return out
    rng = np.random.default_rng(seed)
    idx_by_night = {n: np.flatnonzero(ng == n) for n in nights}
    amps, pr, pb, curves = [], [], [], []
    for _ in range(int(n_boot)):
        pick = rng.choice(len(nights), size=len(nights), replace=True)
        idx = np.concatenate([idx_by_night[nights[i]] for i in pick])
        c2, m2, _ = binned_curve(ph[idx], cv[idx], n_bins, min_count)
        st = _curve_stats(c2, m2)
        if st is None:
            continue
        amps.append(st["amplitude"])
        pr.append(circ_diff(st["phase_red"], base["phase_red"]))
        pb.append(circ_diff(st["phase_blue"], base["phase_blue"]))
        curves.append(m2)
    out["n_boot_ok"] = len(amps)
    if len(amps) >= 50:
        def half(v):
            lo, hi = np.percentile(v, [16.0, 84.0])
            return float((hi - lo) / 2.0)
        out["amplitude_err"] = half(amps)
        out["phase_red_err"] = half(pr)
        out["phase_blue_err"] = half(pb)
        cc = np.array(curves)
        with np.errstate(invalid="ignore"):
            lo = np.nanpercentile(cc, 16.0, axis=0)
            hi = np.nanpercentile(cc, 84.0, axis=0)
        out["bin_err"] = (hi - lo) / 2.0
    return out


def curve_repeatability(med_a, err_a, med_b, err_b) -> dict:
    """Do two binned colour curves have the same SHAPE?

    The two curves (same phase bins) are compared after each has its own
    mean removed over the bins both populate, because the zero point of a
    colour carries the catalogue-tie systematic and the two eras used
    different cameras and filters: the claim under test is about shape.

    Returns ``n_bins``, ``rms_diff`` (mag), ``chi2`` / ``dof`` / ``chi2nu``
    of the difference under the bootstrap bin errors, the Pearson
    correlation ``r`` of the two curves and the amplitude ratio from a
    least-squares scale of b onto a (``scale``; 1 means equal amplitude).
    """
    a = np.asarray(med_a, dtype=float)
    b = np.asarray(med_b, dtype=float)
    ea = np.asarray(err_a, dtype=float)
    eb = np.asarray(err_b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    nan = float("nan")
    out = {"n_bins": int(ok.sum()), "rms_diff": nan, "chi2": nan,
           "dof": None, "chi2nu": nan, "r": nan, "scale": nan}
    if ok.sum() < 6:
        return out
    a0 = a[ok] - a[ok].mean()
    b0 = b[ok] - b[ok].mean()
    d = a0 - b0
    out["rms_diff"] = float(math.sqrt(np.mean(d ** 2)))
    if np.std(a0) > 0 and np.std(b0) > 0:
        out["r"] = float(np.corrcoef(a0, b0)[0, 1])
        out["scale"] = float(np.sum(a0 * b0) / np.sum(a0 * a0))
    s2 = ea[ok] ** 2 + eb[ok] ** 2
    good = np.isfinite(s2) & (s2 > 0)
    if good.sum() >= 6:
        out["chi2"] = float(np.sum(d[good] ** 2 / s2[good]))
        out["dof"] = int(good.sum() - 1)
        out["chi2nu"] = out["chi2"] / out["dof"]
    return out


# ===========================================================================
# 8.  What a recovery contour excludes
# ===========================================================================

def excluded_amplitude_summary(contours_mag, thresholds_mag) -> dict:
    """Turn a set of 90 per cent recovery contours into an honest sentence.

    A recovery contour ``a90`` is the SMALLEST semi-amplitude a blind
    search recovers nine times in ten.  A signal with semi-amplitude at or
    above ``a90`` is excluded (at 90 per cent) by a non-detection in that
    run; one below it is NOT.  The previous draft read a contour above the
    literature floor as licence to say superhumps "would have been
    visible", which inverts it.

    ``contours_mag`` is one entry per run-filter, None/NaN where no
    contour could be measured.  For each threshold the function counts
    the run-filters whose contour is at or below it — the runs in which a
    signal of that size IS excluded.

    Returns ``n_run_filters``, ``n_with_contour``, ``contour_min`` /
    ``contour_max`` / ``contour_median`` (mag) and ``n_excluding``, a dict
    threshold -> count.
    """
    c = np.array([float(x) if x is not None else np.nan
                  for x in contours_mag], dtype=float)
    have = c[np.isfinite(c)]
    out = {"n_run_filters": int(c.size), "n_with_contour": int(have.size),
           "contour_min": None, "contour_max": None, "contour_median": None,
           "n_excluding": {}}
    if have.size:
        out["contour_min"] = float(have.min())
        out["contour_max"] = float(have.max())
        out["contour_median"] = float(np.median(have))
    for th in thresholds_mag:
        out["n_excluding"][float(th)] = int(np.sum(have <= float(th) + 1e-12))
    return out


def half_level_crossing(dt_s, flux, n_bins: int = 40,
                        half_span_s: float = 1200.0) -> dict:
    """Non-parametric time of an egress: where the stacked flux profile
    crosses half-way between its bright and faint levels.

    ``dt_s`` are times relative to a COMMON reference (the ephemeris
    prediction of the edge, the same for every band), ``flux`` the
    relative flux.  The points are median-binned; the bright and faint
    levels are the medians of the first and last quarter of the span; the
    crossing is found by linear interpolation between the two bins that
    bracket the half level, searching outward from the steepest descent.

    This is deliberately a different estimator from the ramp fit — no
    width grid, no time grid, no per-cycle fit — so agreement between the
    two on a band offset is agreement between independent reductions of
    the same photons, and disagreement would convict the ramp model.

    Returns ``t_half_s``, ``bright``, ``faint``, ``depth_frac`` and the
    binned profile; ``t_half_s`` is None when no clean crossing exists.
    """
    x = np.asarray(dt_s, dtype=float)
    f = np.asarray(flux, dtype=float)
    ok = np.isfinite(x) & np.isfinite(f) & (np.abs(x) <= half_span_s)
    x, f = x[ok], f[ok]
    edges = np.linspace(-half_span_s, half_span_s, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(n_bins, np.nan)
    idx = np.clip(np.digitize(x, edges) - 1, 0, n_bins - 1)
    for b in range(n_bins):
        sel = f[idx == b]
        if sel.size >= 3:
            med[b] = float(np.median(sel))
    out = {"t_half_s": None, "bright": None, "faint": None,
           "depth_frac": None, "centres": centres, "median": med}
    q = n_bins // 4
    if np.isfinite(med[:q]).sum() < 2 or np.isfinite(med[-q:]).sum() < 2:
        return out
    bright = float(np.nanmedian(med[:q]))
    faint = float(np.nanmedian(med[-q:]))
    if not bright > faint:
        return out
    half = 0.5 * (bright + faint)
    good = np.flatnonzero(np.isfinite(med))
    # First downward crossing of the half level, scanning in time.
    for i, j in zip(good[:-1], good[1:]):
        if med[i] >= half > med[j]:
            frac = (med[i] - half) / (med[i] - med[j])
            out.update({"t_half_s": float(centres[i]
                                          + frac * (centres[j] - centres[i])),
                        "bright": bright, "faint": faint,
                        "depth_frac": float(1.0 - faint / bright)})
            return out
    return out
