"""Pure arithmetic for CV Phase 3 — the time-series analysis itself.

WHY THIS MODULE EXISTS
----------------------
Phase 1 measured the light curves, Phase 2 decided what was allowed to be
published from them.  Phase 3 asks the questions the project was actually
for, and every one of them has a way of producing a confident wrong answer:

1.  **Period.**  A Lomb-Scargle peak is the true signal convolved with the
    spectral window.  On every multi-night set in this archive the +/-1 c/d
    window sidelobes carry 0.54-0.97 of the central power, so the
    periodogram's tallest peak is a *family* of candidates, not a period.
    :func:`gls_block_power`, :func:`pdm_theta` and :func:`alias_family`
    produce the family; :func:`classify_family_choice` records HOW the
    member was chosen, because "we picked the tallest one" and "the
    literature told us which one" are different claims and only one of
    them is defensible here.
2.  **Timing.**  A per-cycle eclipse or bright-phase edge time is only
    publishable if its uncertainty is known, and the formal delta-chi2 = 1
    error bar from an edge fit is a fiction: the trapezoid does not fit the
    flickering, so chi2_nu comes out at 5-400 and the interval it implies is
    absurdly small.  :func:`fit_edge` therefore rescales to chi2_nu = 1
    before quoting anything, and :func:`sigma_t_injection` measures the real
    sigma_t by injecting known edges into the real cadence with the real
    error model and recovering them.
3.  **O-C.**  An O-C diagram needs a cycle count, and a cycle count 20,000
    cycles from the ephemeris zero point is only unique if the published
    period is known well enough that 20,000 * sigma_P stays below half a
    cycle.  :func:`cycle_number` and :func:`cycle_ambiguity` compute that,
    and the second one is allowed to return "NOT UNIQUE" — an O-C on a
    wrong cycle count is a fabricated result, not a noisy one.
4.  **States.**  A high/low classification needs a threshold that came from
    somewhere.  :func:`otsu_threshold` derives one from the observed
    distribution's own bimodality and :func:`bootstrap_threshold` says how
    well determined it is; :func:`duty_cycle` then uses the Phase-2 upper
    limits so the statistic is not computed on the epochs that were bright
    enough to see.
5.  **Detrending.**  Removing a trend before searching removes signal too.
    :func:`running_median_detrend` is the wrong way, :func:`joint_gp_fit`
    (and its pure-numpy reference :func:`matern32_cov` /
    :func:`gp_log_likelihood_dense`) is the right way, and
    :func:`detrend_suppression` measures the difference on an injected
    signal so the claim is a number rather than an opinion.

Everything here is deterministic and pure — no database, no files, no
network.  ``astropy`` and ``celerite2`` are imported lazily inside the two
functions that need them, so importing this module stays cheap and the unit
tests can run without either.  ``pipeline/scripts/run_cv_phase3.py`` does
all of the I/O and calls in here for every number it stores.

WHAT THIS MODULE IS NOT ALLOWED TO ASSUME
------------------------------------------
The characterization and Phase 2 measured what this data can do, and the
constants below are set against those measurements rather than against the
strategy's hopes:

* per-point precision 9-77 mmag by series, chi2 inflation 0.92-3.02;
* +/-1 c/d aliases at 0.54-0.97 of the window power on every resolved
  multi-night set, so :data:`ALIAS_DECIDABLE_MAX` is set where a periodogram
  could in principle decide alone, and no series in this archive reaches it;
* single-night peaks 3-9 c/d wide, so a single night constrains a period to
  about a part in 30 and nothing better;
* three-filter full-orbit nights: ST LMi 20 of 30, AN UMa 4 of 11, VV Pup 1
  of 18, EU UMa 0 — the inter-band timing in section 3 is an ST LMi result;
* ST LMi's dense 2025 nights sample each filter every 219 s (about 18 points
  per 1.898 h orbit), which is what sets the edge-timing floor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

import numpy as np

# ===========================================================================
# 0.  Tunable constants — the single source of truth.  The report reads
#     these rather than repeating them, so changing one changes the page.
# ===========================================================================

#: Frequency band, in cycles per day, that the SURVEY periodogram covers.
#: Deliberately wider than the orbital band so the page can show what the
#: tallest peak in the whole search actually is.  It is usually not the
#: orbit: these targets change accretion state on week-to-month timescales
#: and that variability piles up below about 3 c/d.
SURVEY_F_MIN_CD = 0.5
SURVEY_F_MAX_CD = 40.0

#: The ORBITAL band.  Periods 0.05-0.125 d.  Every published period in this
#: sample (0.0626 d EU UMa to 0.0868 d YZ Cnc) sits inside it with room for
#: a 1 c/d alias on either side, and no harmonic of a nightly or diurnal
#: cadence lands in it.  Stating the band is the point: a search band chosen
#: after looking at the answer is not a search.
ORBITAL_F_MIN_CD = 8.0
ORBITAL_F_MAX_CD = 20.0

#: Oversampling factor for the frequency grid: the grid step is
#: ``1 / (OVERSAMPLE * baseline)``.  Five is the usual minimum for not
#: missing a peak; ten is used here because the peaks that matter are
#: 1/T wide and the alias structure being resolved is 1 c/d away.
OVERSAMPLE = 10

#: Gap, in hours, that starts a new observing block.  Same value as the
#: characterization's :func:`macro_phot.characterize.night_blocks`: longer
#: than any within-night pause in this archive, far shorter than the ~14 h
#: night-to-night gap.
BLOCK_GAP_H = 3.0

#: Number of phase bins and covers for the Stellingwerf PDM statistic.
#: Ten bins over an orbit is about two bins per ingress for these polars,
#: and two covers (bins offset by half a bin width) stops a real edge from
#: falling on a bin boundary and being averaged away.
PDM_BINS = 10
PDM_COVERS = 2

#: Alias orders tested around every candidate frequency.  +/-3 c/d is
#: further than any real ambiguity in this archive and cheap to test.
ALIAS_ORDERS = (-3, -2, -1, 1, 2, 3)

#: Below this window-power fraction in the strongest alias, a multi-night
#: periodogram could plausibly pick its own period.  Set at 0.20 because a
#: sidelobe a fifth as tall as the peak cannot beat it at any realistic
#: signal-to-noise.  NO series in this archive is below it; the constant
#: exists so that "the data cannot choose" is a threshold comparison rather
#: than an assertion.
ALIAS_DECIDABLE_MAX = 0.20

#: A recovered period AGREES with the published one when the two differ by
#: no more than this many combined sigma.  Three, not one: the published
#: value carries a quoted-precision floor rather than a real error bar,
#: so a tight test would be a test of VSX's rounding.
AGREE_SIGMA = 3.0

#: Edge fitting.  The half-width of the phase window a single edge is fitted
#: in, and the minimum number of points that window must contain.  0.17 in
#: phase is about 2.5 orbital cadence steps on each side of the edge at the
#: measured 219 s sampling — enough to pin both levels without reaching the
#: opposite edge of the bright phase.
EDGE_WINDOW_PHASE = 0.17
EDGE_MIN_POINTS = 6

#: An edge fit is only accepted when the two points that bracket the fitted
#: epoch are no further apart than this multiple of the median cadence.  A
#: wider bracket means the edge fell in a gap and the "measurement" is an
#: interpolation between two levels.
EDGE_MAX_BRACKET_CADENCE = 2.5

#: Minimum edge amplitude, in units of the fitted level scatter, for an
#: edge to be called detected.  Five sigma: an edge is a step, and a step
#: at three sigma in a flickering light curve is a flicker.
EDGE_MIN_SNR = 5.0

#: The publishability threshold for per-cycle timing, in seconds, set by
#: the strategy.  Reported against, never quietly moved.
SIGMA_T_THRESHOLD_S = 60.0

#: Minimum fraction of the orbital phase circle a night must cover before
#: its median magnitude may be called that night's accretion state.  These
#: polars vary by 0.65-1.7 mag around one orbit, which is the same size as
#: a state change, so a median over a third of an orbit is a phase
#: measurement wearing a state label.
STATE_MIN_PHASE_COVERAGE = 0.60

#: Number of phase bins the coverage above is measured in.
STATE_PHASE_BINS = 20

#: Bootstrap replicates for the state threshold and for period errors.
#: Resampling is by NIGHT, not by point, because the noise in these series
#: is correlated within a night and a point bootstrap would report an error
#: bar smaller than the data supports.
N_BOOT = 400

#: Random seed for every Monte Carlo in this module.  Stated so a re-run
#: reproduces the published contour exactly.
SEED = 20260820

#: celerite2's ``Matern32Term`` is not the Matern-3/2 function — it is a
#: two-exponential APPROXIMATION to it, and its own docstring warns that it
#: "should be used with care".  ``eps`` controls the quality: measured
#: against the exact kernel the default eps = 0.01 is wrong by 1.8e-5
#: relative, eps = 1e-4 by 1.8e-9, and eps = 1e-6 by 1.8e-13.  1e-6 is used
#: here so the fast path and the dense reference are the SAME kernel to
#: machine precision and the unit test that compares them is a real test
#: rather than a tolerance chosen to make it pass.
CELERITE_MATERN_EPS = 1e-6


# ===========================================================================
# 1.  Ephemerides and cycle counts
# ===========================================================================

@dataclass(frozen=True)
class Ephemeris:
    """A published ephemeris, with the provenance of every field.

    ``period_sigma_d`` is deliberately Optional.  VSX publishes a period and
    an epoch and NO uncertainty, and the difference between "the catalogue
    did not say" and "the uncertainty is the rounding of the last quoted
    digit" is the difference between an honest cycle count and a hopeful
    one.  When it is None the caller must fall back to
    :func:`quoted_precision_sigma` and SAY that is what it did.
    """

    target_key: str
    name: str
    period_d: float
    epoch_bjd: float
    source: str
    period_str: str = ""
    period_sigma_d: Optional[float] = None
    var_type: str = ""
    note: str = ""


def quoted_precision_sigma(value_str: str) -> float:
    """Uncertainty implied by how many digits a catalogue actually printed.

    ``"0.07908912"`` has eight decimal places, so the value stands for
    anything in the last half-digit either way: sigma = 0.5e-8 d.  This is
    a FLOOR, not an error bar — the true published uncertainty is at least
    this and probably larger — and every result derived from it has to say
    so.  Returns NaN for a string with no decimal point, because a period
    quoted as an integer number of days is not a thing this archive has.
    """
    s = (value_str or "").strip()
    if "." not in s:
        return float("nan")
    # Strip an exponent if one is present; the mantissa carries the digits.
    mantissa = s.split("e")[0].split("E")[0]
    decimals = len(mantissa.split(".")[1])
    return 0.5 * 10.0 ** (-decimals)


def phase_of(times_d, period_d: float, epoch_d: float) -> np.ndarray:
    """Orbital phase in [0, 1) for each time, on the given ephemeris."""
    t = np.asarray(times_d, dtype=float)
    return np.mod((t - epoch_d) / period_d, 1.0)


def cycle_number(times_d, period_d: float, epoch_d: float) -> np.ndarray:
    """Integer cycle index of each time, counted from the ephemeris epoch.

    ``floor`` rather than ``round``: the cycle a point belongs to is the one
    whose phase-zero it has already passed, so cycle boundaries line up with
    phase 0 and a "cycle" is a full orbit rather than a window centred on
    one.
    """
    t = np.asarray(times_d, dtype=float)
    return np.floor((t - epoch_d) / period_d).astype(np.int64)


@dataclass(frozen=True)
class CycleAmbiguity:
    """Is the cycle count between an ephemeris epoch and an observation
    unique?  Every field is a published number."""

    n_cycles: float
    elapsed_d: float
    sigma_period_d: float
    drift_cycles: float
    unique: bool
    sigma_period_max_d: float
    basis: str


def cycle_ambiguity(t_obs_d: float, period_d: float, epoch_d: float,
                    sigma_period_d: float, basis: str = "") -> CycleAmbiguity:
    """How badly the period uncertainty smears the cycle count.

    Between the ephemeris zero point and an observation there are
    ``n = (t_obs - E) / P`` cycles.  An error ``sigma_P`` in the period
    moves the predicted time of cycle ``n`` by ``n * sigma_P``, which is
    ``n * sigma_P / P`` CYCLES.  Once that reaches half a cycle the integer
    count is no longer determined and an O-C diagram built on it is showing
    the arithmetic of a guess.

    ``sigma_period_max_d = 0.5 * P / n`` is the largest period uncertainty
    that would still leave the count unique — the number to quote when the
    catalogue publishes no uncertainty at all, because it converts "we do
    not know sigma_P" into "sigma_P would have to be better than THIS".
    """
    elapsed = float(t_obs_d - epoch_d)
    n = elapsed / period_d
    drift = abs(n) * sigma_period_d / period_d
    sig_max = (0.5 * period_d / abs(n)) if abs(n) > 0 else float("inf")
    return CycleAmbiguity(
        n_cycles=float(n), elapsed_d=elapsed,
        sigma_period_d=float(sigma_period_d), drift_cycles=float(drift),
        unique=bool(drift < 0.5), sigma_period_max_d=float(sig_max),
        basis=basis)


def oc_seconds(t_obs_d, cycle, period_d: float, epoch_d: float) -> np.ndarray:
    """Observed minus computed, in seconds, on a linear ephemeris."""
    t = np.asarray(t_obs_d, dtype=float)
    c = np.asarray(cycle, dtype=float)
    return (t - (epoch_d + c * period_d)) * 86400.0


def fit_linear_ephemeris(cycle, t_obs_d, sigma_t_d
                         ) -> tuple[float, float, float, float, float]:
    """Weighted straight line ``t = E + P * cycle``.

    Returns ``(E, sigma_E, P, sigma_P, chi2_nu)``.  The errors are the
    formal ones from the weighted normal equations RESCALED so chi2_nu = 1
    when chi2_nu > 1 — the same discipline the edge fits use, and for the
    same reason: the scatter of real edge times around a linear ephemeris
    is dominated by the star, not by the photon noise, and an error bar
    that ignores that is not an error bar.
    """
    c = np.asarray(cycle, dtype=float)
    t = np.asarray(t_obs_d, dtype=float)
    s = np.asarray(sigma_t_d, dtype=float)
    ok = np.isfinite(c) & np.isfinite(t) & np.isfinite(s) & (s > 0)
    c, t, s = c[ok], t[ok], s[ok]
    if c.size < 3 or np.ptp(c) == 0:
        nan = float("nan")
        return nan, nan, nan, nan, nan
    w = 1.0 / s ** 2
    a = np.column_stack([np.ones_like(c), c])
    atw = a.T * w
    cov = np.linalg.inv(atw @ a)
    p = cov @ (atw @ t)
    resid = t - a @ p
    chi2 = float(np.sum(w * resid ** 2))
    dof = max(c.size - 2, 1)
    chi2nu = chi2 / dof
    scale = math.sqrt(max(chi2nu, 1.0))
    return (float(p[0]), float(math.sqrt(cov[0, 0])) * scale,
            float(p[1]), float(math.sqrt(cov[1, 1])) * scale, float(chi2nu))


# ===========================================================================
# 2.  Periodograms
# ===========================================================================

def frequency_grid(baseline_d: float, f_min: float = SURVEY_F_MIN_CD,
                   f_max: float = SURVEY_F_MAX_CD,
                   oversample: int = OVERSAMPLE) -> np.ndarray:
    """Uniform frequency grid with step ``1 / (oversample * baseline)``.

    Uniform in FREQUENCY, not period: the peak of a periodogram has a
    constant width in frequency (about 1/T) and a period grid therefore
    over-samples the long periods and under-samples the short ones, which
    is how a search misses a peak it was told to look for.
    """
    T = float(baseline_d)
    if not np.isfinite(T) or T <= 0:
        return np.array([])
    df = 1.0 / (float(oversample) * T)
    n = int(math.ceil((f_max - f_min) / df)) + 1
    return f_min + df * np.arange(n)


def block_index(times_d, gap_h: float = BLOCK_GAP_H) -> np.ndarray:
    """Integer observing-block (night) label for each time, 0-based.

    The times must be sorted; the caller always sorts them because the
    light-curve queries do.  A gap longer than ``gap_h`` starts a new block.
    """
    t = np.asarray(times_d, dtype=float)
    if t.size == 0:
        return np.array([], dtype=np.int64)
    gaps = np.diff(t) * 24.0
    return np.concatenate(([0], np.cumsum(gaps > gap_h))).astype(np.int64)


def _project_out_blocks(columns: np.ndarray, blocks: np.ndarray,
                        weights: np.ndarray) -> np.ndarray:
    """Subtract the weighted per-block mean from every column.

    THIS IS THE WHOLE NUISANCE MODEL, and it is worth being precise about
    what it is and is not.  It is ONE free constant per night, fitted
    simultaneously with the sinusoid — the exact generalisation of the
    "floating mean" that an ordinary Lomb-Scargle already fits, from one
    constant for the whole series to one per night.  Because the block
    indicators are orthogonal to each other, fitting them jointly with any
    other columns is algebraically identical to removing each column's
    weighted block mean first, which is what this function does and why the
    search below can stay fast without becoming approximate.

    It is NOT detrending.  A running median or a spline has a free parameter
    every few points and absorbs power at the frequencies being searched; a
    single constant per night can only absorb power at frequencies below
    1/T_night, about 2.5 c/d for these 9 h runs, five times below the
    orbital band.  Section 5 of the report measures that suppression on an
    injected signal instead of asserting it.
    """
    x = np.atleast_2d(columns)
    out = x.copy()
    # np.bincount is the vectorised group-by: one pass for the weight sums,
    # one per column for the weighted sums.
    wsum = np.bincount(blocks, weights=weights)
    wsum[wsum == 0] = 1.0
    for k in range(out.shape[0]):
        num = np.bincount(blocks, weights=weights * out[k])
        out[k] -= (num / wsum)[blocks]
    return out


def gls_block_power(times_d, y, dy, freqs_cd, blocks=None,
                    chunk: int = 512) -> np.ndarray:
    """Generalised Lomb-Scargle with one free constant PER NIGHT.

    At each trial frequency the model is

        y(t) = A cos(2 pi f t) + B sin(2 pi f t) + c_night(t)

    fitted by weighted least squares, and the returned power is the
    fractional chi-squared improvement over the night-constants-only model:

        power = 1 - chi2(model) / chi2(constants only)

    which runs 0 to 1 exactly like astropy's ``normalization='standard'``.
    The night constants are in the NULL model too, so the power is the
    power the periodic term adds and nothing else.

    Why not just call astropy on night-mean-subtracted data: because that
    projects the data but not the cosine and sine columns, and the two stop
    agreeing exactly where it matters — at frequencies where the sampling
    correlates a sinusoid with the block pattern.  This does the projection
    on all three, which is the exact joint fit.

    Evaluated in chunks of ``chunk`` frequencies so a 100,000-point grid on
    a 700-point series never allocates more than a few tens of megabytes.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(dy, dtype=float)
    f = np.asarray(freqs_cd, dtype=float)
    ok = np.isfinite(t) & np.isfinite(yv) & np.isfinite(e) & (e > 0)
    t, yv, e = t[ok], yv[ok], e[ok]
    if t.size < 5 or f.size == 0:
        return np.full(f.shape, np.nan)
    b = (np.asarray(blocks, dtype=np.int64)[ok] if blocks is not None
         else block_index(t))
    # Re-label blocks densely so bincount arrays stay small after masking.
    _, b = np.unique(b, return_inverse=True)
    w = 1.0 / e ** 2
    # Null model: night constants only.  Its residual is the data with the
    # weighted block means removed, and its chi2 is the denominator.
    y0 = _project_out_blocks(yv[None, :], b, w)[0]
    chi2_null = float(np.sum(w * y0 ** 2))
    if chi2_null <= 0:
        return np.full(f.shape, np.nan)
    # Shift the time origin to the weighted mean: this does not change the
    # fit, but it keeps cos/sin arguments small and the normal equations
    # numerically clean over a 400 d baseline.
    t0 = t - float(np.average(t, weights=w))
    power = np.empty(f.size, dtype=float)
    for start in range(0, f.size, chunk):
        fk = f[start:start + chunk]
        ang = 2.0 * np.pi * np.outer(fk, t0)
        cc = _project_out_blocks(np.cos(ang), b, w)
        ss = _project_out_blocks(np.sin(ang), b, w)
        # 2x2 weighted normal equations, solved analytically per frequency.
        scc = np.einsum("ij,j,ij->i", cc, w, cc)
        sss = np.einsum("ij,j,ij->i", ss, w, ss)
        scs = np.einsum("ij,j,ij->i", cc, w, ss)
        syc = np.einsum("ij,j->i", cc, w * y0)
        sys_ = np.einsum("ij,j->i", ss, w * y0)
        det = scc * sss - scs ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            a_hat = (syc * sss - sys_ * scs) / det
            b_hat = (sys_ * scc - syc * scs) / det
            gain = a_hat * syc + b_hat * sys_
        gain = np.where(np.isfinite(gain) & (det > 0), gain, 0.0)
        power[start:start + chunk] = np.clip(gain / chi2_null, 0.0, 1.0)
    return power


def gls_amplitude(times_d, y, dy, freq_cd, blocks=None) -> float:
    """Best-fit semi-amplitude of the sinusoid at ONE frequency, in mag.

    Same model as :func:`gls_block_power`; returned as
    ``sqrt(A^2 + B^2)`` so it is the semi-amplitude a Lomb-Scargle
    amplitude spectrum would report, which is what the detection-limit
    formula in the strategy means.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(dy, dtype=float)
    ok = np.isfinite(t) & np.isfinite(yv) & np.isfinite(e) & (e > 0)
    t, yv, e = t[ok], yv[ok], e[ok]
    if t.size < 5:
        return float("nan")
    b = (np.asarray(blocks, dtype=np.int64)[ok] if blocks is not None
         else block_index(t))
    _, b = np.unique(b, return_inverse=True)
    w = 1.0 / e ** 2
    y0 = _project_out_blocks(yv[None, :], b, w)[0]
    t0 = t - float(np.average(t, weights=w))
    ang = 2.0 * np.pi * float(freq_cd) * t0
    cc = _project_out_blocks(np.cos(ang)[None, :], b, w)[0]
    ss = _project_out_blocks(np.sin(ang)[None, :], b, w)[0]
    scc = float(np.sum(w * cc * cc))
    sss = float(np.sum(w * ss * ss))
    scs = float(np.sum(w * cc * ss))
    syc = float(np.sum(w * cc * y0))
    sys_ = float(np.sum(w * ss * y0))
    det = scc * sss - scs ** 2
    if not np.isfinite(det) or det <= 0:
        return float("nan")
    a_hat = (syc * sss - sys_ * scs) / det
    b_hat = (sys_ * scc - syc * scs) / det
    return float(math.hypot(a_hat, b_hat))


def pdm_theta(times_d, y, freqs_cd, n_bins: int = PDM_BINS,
              n_covers: int = PDM_COVERS, blocks=None) -> np.ndarray:
    """Stellingwerf phase-dispersion minimisation statistic theta(f).

    theta is the pooled within-phase-bin variance divided by the overall
    variance.  A correct period gathers similar magnitudes into the same
    bin, so theta drops well below 1; a wrong one scatters them and theta
    sits at 1.  Unlike Lomb-Scargle it assumes NOTHING about the shape of
    the light curve, which is why it is here: a polar's bright phase is a
    top hat with an edge, and a sinusoid fitted to a top hat puts most of
    its power in harmonics where the fundamental search cannot see it.

    ``blocks`` removes one constant per night first, for the same reason
    and by the same argument as :func:`gls_block_power`.  Overlapping
    covers (bins shifted by ``1/(n_bins*n_covers)``) stop a real edge that
    happens to land on a bin boundary from being averaged into invisibility.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    f = np.asarray(freqs_cd, dtype=float)
    ok = np.isfinite(t) & np.isfinite(yv)
    t, yv = t[ok], yv[ok]
    if t.size < n_bins + 2 or f.size == 0:
        return np.full(f.shape, np.nan)
    if blocks is not None:
        b = np.asarray(blocks, dtype=np.int64)[ok]
        _, b = np.unique(b, return_inverse=True)
        yv = _project_out_blocks(yv[None, :], b, np.ones_like(yv))[0]
    yv = yv - yv.mean()
    total_var = float(np.var(yv, ddof=1))
    if total_var <= 0:
        return np.full(f.shape, np.nan)
    width = 1.0 / n_bins
    starts = np.arange(n_bins * n_covers) * (width / n_covers)
    out = np.empty(f.size, dtype=float)
    for k, fk in enumerate(f):
        ph = np.mod(t * fk, 1.0)
        num = 0.0
        den = 0
        for s0 in starts:
            sel = np.mod(ph - s0, 1.0) < width
            n = int(sel.sum())
            if n > 1:
                num += (n - 1) * float(np.var(yv[sel], ddof=1))
                den += n - 1
        out[k] = (num / den) / total_var if den > 0 else np.nan
    return out


def spectral_window_power(times_d, freqs_cd) -> np.ndarray:
    """``W(f) = |sum exp(-2 pi i f t)|^2 / N^2``, normalised to W(0) = 1.

    Re-exported here (identical to
    :func:`macro_phot.characterize.spectral_window`) so that every figure on
    the Phase-3 page draws its window from the same function that computes
    the alias fractions in its caption.  The window is the periodogram of a
    perfectly constant star observed at these instants: every peak in it is
    a frequency at which the SAMPLING makes power, and a real periodogram is
    the truth convolved with it.
    """
    t = np.asarray(times_d, dtype=float)
    t = t[np.isfinite(t)]
    f = np.asarray(freqs_cd, dtype=float)
    if t.size == 0:
        return np.zeros_like(f)
    t = t - t.mean()
    return np.abs(np.exp(-2j * np.pi * np.outer(f, t)).sum(axis=1)) ** 2 / t.size ** 2


def alias_family(f_cd: float, f_alias_cd: float = 1.0,
                 orders: Sequence[int] = ALIAS_ORDERS) -> list[tuple[int, float]]:
    """The frequencies a signal at ``f_cd`` is confusable with.

    Returns ``(order, frequency)`` pairs for ``f + k * f_alias``, folding
    negative results to their absolute value because a periodogram cannot
    distinguish +f from -f, and dropping anything at or below zero.
    """
    out = []
    for k in orders:
        fk = f_cd + k * f_alias_cd
        if abs(fk) > 1e-9:
            out.append((int(k), abs(float(fk))))
    return out


def alias_window_fractions(times_d, f_alias_cd: float = 1.0,
                           orders: Sequence[int] = ALIAS_ORDERS
                           ) -> dict[int, float]:
    """Window power at each alias OFFSET, i.e. the height of the sidelobe
    an alias would produce relative to the true peak.

    Evaluated at ``k * f_alias``, not at ``f_true + k * f_alias``: the peak
    an alias makes has height ``W(f_peak - f_true)``, and evaluating the
    window at the alias frequency itself would answer a different question.
    A value near 1 means the sidelobe is as tall as the truth and the
    period is not determinable from this sampling alone.
    """
    offs = np.array([k * f_alias_cd for k in orders], dtype=float)
    w = spectral_window_power(times_d, offs)
    return {int(k): float(v) for k, v in zip(orders, w)}


def refine_peak(freqs_cd, power, f_target: float, halfwidth_cd: float
                ) -> tuple[float, float]:
    """Highest point of ``power`` within ``halfwidth_cd`` of ``f_target``,
    with a parabolic interpolation through its two neighbours.

    Returns ``(frequency, power)``.  The interpolation matters: the grid
    step is 1/(10 T) and the peak is 1/T wide, so the grid maximum can be
    up to a twentieth of a peak width from the true maximum, which at a
    400 d baseline is already larger than the period uncertainty being
    quoted.  Returns NaN when the window contains no finite point.
    """
    f = np.asarray(freqs_cd, dtype=float)
    p = np.asarray(power, dtype=float)
    sel = np.isfinite(p) & (np.abs(f - f_target) <= halfwidth_cd)
    if not sel.any():
        return float("nan"), float("nan")
    idx = np.flatnonzero(sel)
    i = int(idx[int(np.argmax(p[idx]))])
    if 0 < i < f.size - 1 and np.isfinite(p[i - 1]) and np.isfinite(p[i + 1]):
        y0, y1, y2 = p[i - 1], p[i], p[i + 1]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            shift = 0.5 * (y0 - y2) / denom          # in grid steps
            shift = float(np.clip(shift, -1.0, 1.0))
            df = f[i + 1] - f[i]
            return float(f[i] + shift * df), float(y1)
    return float(f[i]), float(p[i])


def refine_trough(freqs_cd, theta, f_target: float, halfwidth_cd: float
                  ) -> tuple[float, float]:
    """The PDM equivalent of :func:`refine_peak`: lowest theta near a target.

    Implemented by refining the peak of ``-theta`` so that exactly one
    interpolation lives in this module and the two statistics cannot drift
    apart.
    """
    f, p = refine_peak(freqs_cd, -np.asarray(theta, dtype=float),
                       f_target, halfwidth_cd)
    return f, (-p if np.isfinite(p) else p)


def peak_halfwidth(freqs_cd, power, f_peak: float,
                   frac: float = 0.5) -> float:
    """Half width at ``frac`` of maximum of the peak containing ``f_peak``.

    Walks outward from the peak in both directions until the power drops
    below ``frac`` times its height, and returns half the resulting full
    width.  This is the number that says what a SINGLE NIGHT can do: a
    9 h run gives a peak about 1/0.4 d = 2.7 c/d wide, so the +/-1 c/d
    aliases are not separate features at all — the whole family is one
    blur, and the period is constrained to roughly a part in 5.
    """
    f = np.asarray(freqs_cd, dtype=float)
    p = np.asarray(power, dtype=float)
    if f.size < 3:
        return float("nan")
    i = int(np.argmin(np.abs(f - f_peak)))
    half = frac * p[i]
    lo = i
    while lo > 0 and np.isfinite(p[lo]) and p[lo] > half:
        lo -= 1
    hi = i
    while hi < f.size - 1 and np.isfinite(p[hi]) and p[hi] > half:
        hi += 1
    return float((f[hi] - f[lo]) / 2.0)


def mhb_period_sigma(n_points: int, baseline_d: float, amplitude_mag: float,
                     residual_rms_mag: float, period_d: float) -> float:
    """Analytic period uncertainty, Montgomery & O'Donoghue (1999).

    ``sigma_f = sqrt(6 / N) * sigma_m / (pi * T * A)`` and
    ``sigma_P = P^2 sigma_f``.  It is the Cramer-Rao bound for a sinusoid
    of semi-amplitude ``A`` in white noise of scatter ``sigma_m`` sampled
    ``N`` times over baseline ``T``, and it is quoted here as the OPTIMISTIC
    end of the error bar for exactly that reason: the noise in these series
    is not white (Allan slopes shallower than -0.5 on a third of the
    ladders) and the sampling is not uniform.  The published uncertainty is
    the larger of this and the night bootstrap.
    """
    n = float(n_points)
    T = float(baseline_d)
    a = float(amplitude_mag)
    s = float(residual_rms_mag)
    if n < 4 or T <= 0 or a <= 0 or s <= 0:
        return float("nan")
    sigma_f = math.sqrt(6.0 / n) * s / (math.pi * T * a)
    return float(period_d ** 2 * sigma_f)


def night_bootstrap_period(times_d, y, dy, blocks, f_target: float,
                           halfwidth_cd: float, oversample: int = OVERSAMPLE,
                           n_boot: int = N_BOOT, seed: int = SEED
                           ) -> tuple[float, float, int]:
    """Period uncertainty from resampling NIGHTS with replacement.

    Returns ``(sigma_period_d, median_period_d, n_ok)``.

    Nights, not points.  A point bootstrap treats every measurement as an
    independent draw, which the Allan analysis says they are not, and it
    would return an error bar that shrinks with the number of exposures
    instead of with the number of independent nights.  Each replicate
    re-runs the same local GLS refinement in a window around ``f_target``,
    so what is being propagated is the whole measurement, not a linearised
    version of it.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(dy, dtype=float)
    b = np.asarray(blocks, dtype=np.int64)
    labels = np.unique(b)
    if labels.size < 3:
        return float("nan"), float("nan"), 0
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(b == lab) for lab in labels]
    out = []
    for _ in range(int(n_boot)):
        pick = rng.integers(0, len(groups), len(groups))
        idx = np.concatenate([groups[i] for i in pick])
        # Re-label so two draws of the same night are two independent
        # nuisance constants rather than one shared one; sharing would
        # pretend the resampled series has fewer free parameters than it has.
        newb = np.concatenate([np.full(groups[i].size, j)
                               for j, i in enumerate(pick)])
        tb, yb, eb = t[idx], yv[idx], e[idx]
        order = np.argsort(tb, kind="mergesort")
        tb, yb, eb, nb = tb[order], yb[order], eb[order], newb[order]
        base = float(tb.max() - tb.min())
        if base <= 0:
            continue
        grid = frequency_grid(base, f_target - halfwidth_cd,
                              f_target + halfwidth_cd, oversample)
        if grid.size < 5:
            continue
        p = gls_block_power(tb, yb, eb, grid, nb)
        fr, _ = refine_peak(grid, p, f_target, halfwidth_cd)
        if np.isfinite(fr) and fr > 0:
            out.append(1.0 / fr)
    if len(out) < 10:
        return float("nan"), float("nan"), len(out)
    arr = np.asarray(out, dtype=float)
    return float(np.std(arr, ddof=1)), float(np.median(arr)), int(arr.size)


def classify_family_choice(n_blocks: int, max_alias_fraction: float,
                           peak_halfwidth_cd: float,
                           f_alias_cd: float = 1.0,
                           prior_available: bool = True) -> tuple[str, str]:
    """HOW the alias family member was chosen.  Returns ``(code, sentence)``.

    THE ALIAS RULE IS BINDING and this function is where it is enforced.
    Three outcomes and no fourth:

    ``DATA``
        the sidelobes are below :data:`ALIAS_DECIDABLE_MAX` and the
        periodogram could pick its own member.  No series in this archive
        reaches this; the branch exists so the claim stays falsifiable.
    ``SINGLE-NIGHT``
        one observing block, whose peak is wider than the 1 c/d alias
        spacing.  There is then no alias AMBIGUITY — the family members are
        not separate features — but there is no precision either, and the
        period is only good to about the peak half-width.
    ``PRIOR``
        multi-night, sidelobes tall, family members individually resolved:
        the member has to be chosen by the published ephemeris, and the
        recovered period is then a CONFIRMATION of the literature value at
        the precision of the local peak, not an independent determination.

    When ``prior_available`` is False the PRIOR branch becomes
    ``UNRESOLVED``, which is the honest answer for a target with no
    published ephemeris and tall sidelobes.
    """
    if not np.isfinite(max_alias_fraction):
        return "UNRESOLVED", ("The spectral window could not be evaluated, "
                              "so no statement about aliasing is available.")
    if max_alias_fraction < ALIAS_DECIDABLE_MAX:
        return "DATA", (
            f"The strongest sidelobe carries {max_alias_fraction:.2f} of the "
            f"window power, below the {ALIAS_DECIDABLE_MAX:.2f} at which a "
            "periodogram can distinguish a peak from its own alias, so this "
            "sampling selects the family member on its own.")
    if n_blocks <= 1:
        return "SINGLE-NIGHT", (
            f"One observing block.  The peak is {peak_halfwidth_cd:.2f} c/d "
            f"wide at half maximum against a {f_alias_cd:.0f} c/d alias "
            "spacing, so the family members are not separate features and "
            "there is nothing to choose between; the cost is that the "
            "period is constrained only to about the peak width.")
    if not prior_available:
        return "UNRESOLVED", (
            f"Multi-night sampling with the strongest sidelobe at "
            f"{max_alias_fraction:.2f} of the window power and no published "
            "ephemeris to select a family member.  The period is NOT "
            "resolvable from this data set.")
    return "PRIOR", (
        f"Multi-night sampling with the strongest sidelobe at "
        f"{max_alias_fraction:.2f} of the window power, so the periodogram "
        "cannot choose among the +/-1 c/d family.  The member was selected "
        "by the published ephemeris; what is measured here is agreement "
        "with that value at the precision of the local peak, not an "
        "independent period determination.")


def agreement(period_d: float, sigma_d: float, published_d: float,
              published_sigma_d: float) -> tuple[float, bool, str]:
    """``(deviation_in_sigma, agrees, sentence)`` against a published value.

    The two uncertainties are added in quadrature; when either is missing
    the comparison is made on the one that exists and the sentence says so,
    because a NaN sigma silently dropping out of a quadrature sum is how a
    disagreement gets published as an agreement.
    """
    if not (np.isfinite(period_d) and np.isfinite(published_d)):
        return float("nan"), False, "No recovered period to compare."
    s1 = sigma_d if np.isfinite(sigma_d) else 0.0
    s2 = published_sigma_d if np.isfinite(published_sigma_d) else 0.0
    comb = math.hypot(s1, s2)
    diff = float(period_d - published_d)
    if comb <= 0:
        return float("nan"), False, (
            f"Difference {diff * 86400.0:+.3f} s, but neither value carries "
            "an uncertainty, so no agreement test is possible.")
    dev = diff / comb
    ok = abs(dev) <= AGREE_SIGMA
    return (float(dev), bool(ok),
            f"{'Agrees' if ok else 'DISAGREES'}: recovered minus published "
            f"= {diff * 86400.0:+.4f} s = {dev:+.2f} sigma against the "
            f"{AGREE_SIGMA:.0f} sigma bar.")


# ===========================================================================
# 3.  Edge timing
# ===========================================================================

def ramp(times_d, t_edge_d: float, width_d: float) -> np.ndarray:
    """Unit ramp: 0 before ``t_edge - width/2``, 1 after ``+ width/2``.

    Linear in between.  A ramp and not a step, because the quantity being
    measured is the EPOCH of the transition and a step's epoch is only ever
    localised to the gap between two exposures, whereas a ramp lets a point
    that caught the transition half-way carry information about where in
    that gap it happened.  A polar's bright-phase edge really is a ramp:
    the accretion column takes a finite time to rotate behind the white
    dwarf's limb.
    """
    t = np.asarray(times_d, dtype=float)
    x = (t - float(t_edge_d)) / max(float(width_d), 1e-12) + 0.5
    return np.clip(x, 0.0, 1.0)


@dataclass(frozen=True)
class EdgeFit:
    """One measured bright-phase (or eclipse) edge epoch."""

    t_edge_d: float
    sigma_t_s: float
    width_d: float
    level_bright: float
    level_faint: float
    depth_mag: float
    depth_snr: float
    chi2nu: float
    n_points: int
    bracket_s: float
    accepted: bool
    reason: str


def fit_edge(times_d, y, dy, t_grid_d, width_grid_d,
             median_cadence_s: float,
             min_points: int = EDGE_MIN_POINTS,
             min_snr: float = EDGE_MIN_SNR,
             max_bracket_cadence: float = EDGE_MAX_BRACKET_CADENCE,
             fixed_depth: Optional[float] = None) -> EdgeFit:
    """Fit one edge epoch by profiling over its time and its width.

    The model has four parameters — bright level, faint level, edge epoch,
    edge width — but only two of them are non-linear.  For every
    ``(t_edge, width)`` on the grid the two LEVELS follow from a two-column
    weighted linear solve, so the search is a clean 2-D chi-squared surface
    with no optimiser to get stuck.  Deterministic, and the surface itself
    is the error bar.

    ``fixed_depth`` switches the fit from MEASUREMENT to TEMPLATE MATCHING:
    the step between the two levels is held at the supplied value and only
    the overall level floats, so the design has one free linear parameter
    instead of two.  Real edges are fitted with it None — nobody knows the
    depth in advance.  The sigma_t injection test sets it, because "the
    depth is 20% wrong" is only a meaningful error axis if the recovery is
    actually committed to a depth; with both levels free, a wrong assumed
    depth costs nothing and the test would report a reassuring nothing.

    THE ERROR BAR IS RESCALED, and this is the most important line in the
    function.  The raw fits on the ST LMi dense night return chi2_nu of 5
    to 400 on 6-11 points: the trapezoid does not describe the flickering,
    and a delta-chi2 = 1 interval on a model that misfits by that much
    returns 0-4 s, which would be a published timing precision an order of
    magnitude better than the cadence allows.  Rescaling the errors so
    chi2_nu = 1 before taking the interval converts "the model is wrong" into
    "the error bar is bigger", which is the only honest of the two.  Even
    then the rescaled bar is a FORMAL number; the authority on sigma_t is
    :func:`sigma_t_injection`.

    Returns an :class:`EdgeFit` whose ``accepted`` flag is False, with a
    ``reason``, whenever the fit is not a measurement: too few points, the
    edge landing in a gap, or a step too small to be distinguished from a
    flicker.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(dy, dtype=float)
    ok = np.isfinite(t) & np.isfinite(yv) & np.isfinite(e) & (e > 0)
    t, yv, e = t[ok], yv[ok], e[ok]
    nan = float("nan")
    if t.size < min_points:
        return EdgeFit(nan, nan, nan, nan, nan, nan, nan, nan, int(t.size),
                       nan, False, f"only {t.size} usable points "
                                   f"(need {min_points})")
    w = 1.0 / e ** 2
    tg = np.asarray(t_grid_d, dtype=float)
    wg = np.asarray(width_grid_d, dtype=float)
    chi2 = np.full((tg.size, wg.size), np.inf)
    levels = np.full((tg.size, wg.size, 2), np.nan)
    n_free = 4 if fixed_depth is None else 3
    for i, te in enumerate(tg):
        for j, wid in enumerate(wg):
            r = ramp(t, te, wid)
            if fixed_depth is None:
                # Design: [1, ramp].  Coefficients are (bright level, step).
                a = np.column_stack([np.ones_like(t), r])
                target = yv
            else:
                # Template matching: the step is known, so subtract it and
                # fit only the level the template floats at.
                a = np.ones((t.size, 1))
                target = yv - float(fixed_depth) * r
            atw = a.T * w
            m = atw @ a
            if not np.isfinite(m).all() or abs(np.linalg.det(m)) < 1e-30:
                continue
            p = np.linalg.solve(m, atw @ target)
            res = target - a @ p
            chi2[i, j] = float(np.sum(w * res ** 2))
            levels[i, j] = (p if fixed_depth is None
                            else np.array([p[0], float(fixed_depth)]))
    if not np.isfinite(chi2).any():
        return EdgeFit(nan, nan, nan, nan, nan, nan, nan, nan, int(t.size),
                       nan, False, "no invertible fit on the grid")
    i0, j0 = np.unravel_index(int(np.argmin(chi2)), chi2.shape)
    chi_min = float(chi2[i0, j0])
    dof = max(t.size - n_free, 1)
    chi2nu = chi_min / dof
    # Profile over width, then read the interval off the time axis.  The
    # threshold is delta-chi2 = 1 SCALED by chi2_nu, which is exactly the
    # same as rescaling every error bar by sqrt(chi2_nu) and using 1.
    prof = chi2.min(axis=1)
    thresh = chi_min + max(chi2nu, 1.0)
    inside = np.flatnonzero(prof <= thresh)
    if inside.size >= 2:
        sigma_s = float((tg[inside[-1]] - tg[inside[0]]) / 2.0 * 86400.0)
    else:
        sigma_s = float((tg[1] - tg[0]) * 86400.0) if tg.size > 1 else nan
    t_edge = float(tg[i0])
    width = float(wg[j0])
    bright, step = float(levels[i0, j0, 0]), float(levels[i0, j0, 1])
    faint = bright + step
    # Scatter of the residuals sets the scale a step has to beat — but never
    # a scale SMALLER than the quoted photometric errors.  Without that
    # floor a noiseless (or very well fitted) edge divides by a residual
    # scatter of zero and returns a NaN signal-to-noise, which the gate
    # below then reads as "too small to see" and throws the best
    # measurement in the set away.  A step is never known better than the
    # error bars on the points that define it.
    r = ramp(t, t_edge, width)
    resid = yv - (bright + step * r)
    scatter = float(np.std(resid, ddof=1)) if t.size > 2 else nan
    noise = max(scatter if np.isfinite(scatter) else 0.0,
                float(np.median(e)))
    snr = abs(step) / noise if noise > 0 else float("inf")
    # Bracket: how far apart are the points either side of the fitted epoch?
    before = t[t <= t_edge]
    after = t[t >= t_edge]
    if before.size and after.size:
        bracket_s = float((after.min() - before.max()) * 86400.0)
    else:
        bracket_s = float("inf")
    reason = "accepted"
    accepted = True
    if not np.isfinite(snr) or snr < min_snr:
        accepted, reason = False, (f"step SNR {snr:.1f} below {min_snr:.0f} "
                                   "— not distinguishable from flickering")
    elif not np.isfinite(bracket_s) or (
            np.isfinite(median_cadence_s) and median_cadence_s > 0 and
            bracket_s > max_bracket_cadence * median_cadence_s):
        accepted, reason = False, (
            f"edge fell in a {bracket_s:.0f} s gap, more than "
            f"{max_bracket_cadence:g}x the {median_cadence_s:.0f} s cadence "
            "— the epoch is an interpolation, not a measurement")
    elif i0 in (0, tg.size - 1):
        accepted, reason = False, "best epoch sits on the edge of the search grid"
    return EdgeFit(t_edge_d=t_edge, sigma_t_s=sigma_s, width_d=width,
                   level_bright=bright, level_faint=faint,
                   depth_mag=float(step), depth_snr=float(snr),
                   chi2nu=float(chi2nu), n_points=int(t.size),
                   bracket_s=bracket_s, accepted=accepted, reason=reason)


def edge_time_grid(t_centre_d: float, half_span_d: float, n: int = 481
                   ) -> np.ndarray:
    """Search grid for the edge epoch.  Odd ``n`` so the guess is on it."""
    return np.linspace(t_centre_d - half_span_d, t_centre_d + half_span_d, n)


def band_difference(t_a_d, sig_a_s, t_b_d, sig_b_s
                    ) -> tuple[float, float, float]:
    """Weighted mean inter-band edge-time difference, in seconds.

    Returns ``(delta_s, sigma_s, chi2nu)`` for the per-cycle differences
    ``t_a - t_b``.  Paired per cycle rather than differencing two means:
    the two bands see the same cycle-to-cycle wander of the accretion spot,
    so pairing cancels it and leaves the systematic band offset, which is
    the cyclotron quantity.  ``chi2nu`` is reported because a value far
    above 1 means the offset is not constant from cycle to cycle, which is
    itself the interesting result and must not be hidden inside a smaller
    error bar.
    """
    ta = np.asarray(t_a_d, dtype=float)
    tb = np.asarray(t_b_d, dtype=float)
    sa = np.asarray(sig_a_s, dtype=float)
    sb = np.asarray(sig_b_s, dtype=float)
    ok = (np.isfinite(ta) & np.isfinite(tb) & np.isfinite(sa)
          & np.isfinite(sb) & (sa > 0) & (sb > 0))
    if ok.sum() < 1:
        return float("nan"), float("nan"), float("nan")
    d = (ta[ok] - tb[ok]) * 86400.0
    s = np.hypot(sa[ok], sb[ok])
    w = 1.0 / s ** 2
    mean = float(np.sum(w * d) / np.sum(w))
    var = 1.0 / float(np.sum(w))
    if d.size > 1:
        chi2nu = float(np.sum(w * (d - mean) ** 2) / (d.size - 1))
        # Same discipline as everywhere else: a misfit inflates the bar.
        sigma = math.sqrt(var * max(chi2nu, 1.0))
    else:
        chi2nu = float("nan")
        sigma = math.sqrt(var)
    return mean, sigma, chi2nu


# ===========================================================================
# 4.  The sigma_t injection test
# ===========================================================================

def bright_phase_template(times_d, period_d: float, t0_d: float,
                          depth_mag: float, bright_width_phase: float,
                          edge_width_d: float) -> np.ndarray:
    """A polar's bright phase: a raised plateau with two ramped edges.

    Returned in magnitudes as a NEGATIVE offset over the bright interval
    (brighter = smaller magnitude), so it adds to a faint baseline.  ``t0``
    is the centre of the bright phase.  Both edges get the same ramp width,
    which is the assumption the recovery is later tested against by being
    deliberately violated.
    """
    t = np.asarray(times_d, dtype=float)
    p = float(period_d)
    # Signed phase offset from the bright-phase centre, in DAYS, folded to
    # the nearest cycle so a template evaluated over many cycles repeats.
    dt = np.mod(t - t0_d + p / 2.0, p) - p / 2.0
    half = bright_width_phase * p / 2.0
    w = max(float(edge_width_d), 1e-9)
    # Rising side and falling side, each a clipped linear ramp; their
    # product is the plateau.  Written as two ramps rather than a boolean
    # mask precisely so the edges have finite slope.
    up = np.clip((dt + half) / w + 0.5, 0.0, 1.0)
    down = np.clip((half - dt) / w + 0.5, 0.0, 1.0)
    return -float(depth_mag) * up * down


def sigma_t_injection(times_d, sigma_mag, period_d: float,
                      depth_mag: float, edge_width_d: float,
                      bright_width_phase: float,
                      shape_error: float = 1.0, depth_error: float = 0.0,
                      n_real: int = 200, median_cadence_s: float = 219.0,
                      seed: int = SEED,
                      residuals: Optional[np.ndarray] = None
                      ) -> dict:
    """The measurement that decides whether per-cycle timing is publishable.

    Injects a bright phase of known epoch into the REAL timestamps with the
    REAL per-point error model, then recovers the falling edge with
    :func:`fit_edge` using a template that may be deliberately WRONG:

    ``shape_error``
        multiplies the edge width the recovery assumes relative to the one
        that was injected.  1.0 is the edge shape known exactly; 5.0 is the
        ramp assumed five times longer than it is.
    ``depth_error``
        fractional error in the assumed depth: the recovery template is
        committed to ``depth * (1 + depth_error)``.

    Both errors bite because the recovery is TEMPLATE MATCHING — a single
    fixed width and a fixed depth, with only the edge epoch and an overall
    level free.  That is the whole point.  An earlier version of this
    function let the recovery grid over the width and float both levels;
    it reported that a five-times-wrong shape made no difference, which was
    true of that fit and false of the question, because a fit free to
    re-derive the shape has not been told a wrong one.

    The reason this exists rather than a formal error bar: the fitted
    delta-chi2 = 1 interval on real data comes out at 0-4 s because the
    trapezoid misfits the flickering.  Here the truth is known, so the
    scatter of ``t_recovered - t_true`` IS sigma_t, with no model of the
    star in it.

    ``residuals`` lets the caller pass REAL residuals to be cyclically
    rolled instead of drawing Gaussian noise, which preserves the measured
    correlated structure.  Rolling rather than shuffling, for the same
    reason the Phase-2 injection tests roll: shuffling whitens the noise and
    every timing precision computed on whitened noise is better than the
    data can deliver.

    Returns a dict with ``sigma_t_s`` (the robust 1-sigma, from the 16-84
    percentile half-range so a few catastrophic misses cannot dominate),
    ``bias_s``, ``total_error_s``, ``sigma_t_rms_s``, ``p95_abs_s``,
    ``n_ok``, ``n_try`` and ``recovered_fraction``.

    ``total_error_s = hypot(bias, sigma_t)`` is the number the verdict is
    taken on, and the distinction matters.  A wrong template does not just
    scatter the recovered epochs, it SHIFTS them: assuming a five-times-too-
    wide ramp actually produces a TIGHTER spread than the correct template,
    because a smooth model fits smoothly, while placing every epoch several
    seconds late.  Judging on the scatter alone would have scored the badly
    wrong template as the better measurement.  A bias common to every cycle
    cancels out of an O-C gradient, but it does not cancel out of an O-C
    zero point and it certainly does not cancel between two bands whose
    templates differ — which is exactly the inter-band measurement in
    section 3.
    """
    t = np.sort(np.asarray(times_d, dtype=float))
    sig = np.asarray(sigma_mag, dtype=float)
    if sig.size == 1:
        sig = np.full(t.shape, float(sig))
    ok = np.isfinite(t) & np.isfinite(sig) & (sig > 0)
    t, sig = t[ok], sig[ok]
    nan = float("nan")
    if t.size < EDGE_MIN_POINTS:
        return {"sigma_t_s": nan, "bias_s": nan, "sigma_t_rms_s": nan,
                "p95_abs_s": nan, "n_ok": 0, "n_try": 0,
                "recovered_fraction": nan}
    rng = np.random.default_rng(seed)
    res = (np.asarray(residuals, dtype=float)
           if residuals is not None and np.asarray(residuals).size >= t.size
           else None)
    span = float(t.max() - t.min())
    # Recovery assumes these, and they are wrong on purpose.
    assumed_width = edge_width_d * float(shape_error)
    assumed_depth = depth_mag * (1.0 + float(depth_error))
    deltas = []
    n_try = 0
    for _ in range(int(n_real)):
        # Random bright-phase centre so the edge lands at a random place
        # relative to the sampling; that placement is a real and dominant
        # part of sigma_t at 219 s cadence and must not be averaged away by
        # always injecting at the same phase.
        t0 = float(t.min() + rng.uniform(0.0, period_d))
        model = bright_phase_template(t, period_d, t0, depth_mag,
                                      bright_width_phase, edge_width_d)
        if res is not None:
            noise = np.roll(res[:t.size], int(rng.integers(0, t.size)))
        else:
            noise = rng.normal(0.0, sig)
        yv = 18.0 + model + noise
        # Truth: the falling edge is half a bright-phase width AFTER centre.
        t_true = t0 + bright_width_phase * period_d / 2.0
        while t_true < t.min():
            t_true += period_d
        if t_true > t.max():
            continue
        n_try += 1
        half_span = 3.0 * median_cadence_s / 86400.0
        tg = edge_time_grid(t_true, half_span, 241)
        # ONE width and ONE depth: the recovery is committed to the template
        # it was given, right or wrong.
        wg = np.array([assumed_width])
        sel = np.abs(t - t_true) <= EDGE_WINDOW_PHASE * period_d
        if sel.sum() < EDGE_MIN_POINTS:
            continue
        # The SAME acceptance rule the real pipeline applies, so the
        # recovered fraction below is the fraction of cycles that would
        # actually yield a published epoch, not an idealised one.
        fit = fit_edge(t[sel], yv[sel], sig[sel], tg, wg,
                       median_cadence_s, min_snr=EDGE_MIN_SNR,
                       fixed_depth=assumed_depth)
        if not fit.accepted:
            continue
        deltas.append((fit.t_edge_d - t_true) * 86400.0)
    if len(deltas) < 5:
        return {"sigma_t_s": nan, "bias_s": nan, "total_error_s": nan,
                "sigma_t_rms_s": nan, "p95_abs_s": nan,
                "n_ok": len(deltas), "n_try": n_try,
                "recovered_fraction": (len(deltas) / n_try) if n_try else nan}
    d = np.asarray(deltas, dtype=float)
    lo, hi = np.percentile(d, [16.0, 84.0])
    scatter = float((hi - lo) / 2.0)
    bias = float(np.median(d))
    return {"sigma_t_s": scatter,
            "bias_s": bias,
            "total_error_s": float(math.hypot(bias, scatter)),
            "sigma_t_rms_s": float(np.std(d, ddof=1)),
            "p95_abs_s": float(np.percentile(np.abs(d), 95.0)),
            "n_ok": int(d.size), "n_try": int(n_try),
            "recovered_fraction": float(d.size / n_try) if n_try else nan,
            "span_d": span}


def contour_verdict(grid: Sequence[dict],
                    threshold_s: float = SIGMA_T_THRESHOLD_S,
                    key: str = "total_error_s") -> tuple[str, str]:
    """Turn the sigma_t grid into one publishable sentence.

    ``grid`` is a list of dicts each carrying ``shape_error``,
    ``depth_error`` and ``total_error_s`` (scatter and template-induced bias
    added in quadrature — see :func:`sigma_t_injection` for why the verdict
    is not taken on the scatter alone).  Three verdicts:

    ``PUBLISHABLE``
        every cell of the grid is under the threshold, so the timing holds
        however wrong the assumed template is within the range tested.
    ``CONDITIONAL``
        the exact-shape cells pass and some wrong-shape cells do not: the
        timing is publishable ONLY with a stated template, and the template
        is then part of the result.
    ``NOT PUBLISHABLE``
        even the exact-shape cell misses.
    """
    cells = [g for g in grid if np.isfinite(g.get(key, np.nan))]
    if not cells:
        return "NO RESULT", ("No cell of the grid produced a usable "
                             "recovery, so sigma_t is not measured.")
    worst = max(c[key] for c in cells)
    exact = [c for c in cells
             if abs(c.get("shape_error", 1.0) - 1.0) < 1e-9
             and abs(c.get("depth_error", 0.0)) < 1e-9]
    best = min((c[key] for c in exact), default=float("nan"))
    n_pass = sum(1 for c in cells if c[key] <= threshold_s)
    frac = n_pass / len(cells)
    if worst <= threshold_s:
        return "PUBLISHABLE", (
            f"Every cell of the grid returns sigma_t <= {worst:.1f} s "
            f"against the {threshold_s:.0f} s threshold, including the cells "
            "where the template shape and depth are deliberately wrong.")
    if np.isfinite(best) and best <= threshold_s:
        return "CONDITIONAL", (
            f"With the edge shape known exactly sigma_t = {best:.1f} s, "
            f"inside the {threshold_s:.0f} s threshold; with the shape and "
            f"depth wrong it degrades to {worst:.1f} s.  "
            f"{n_pass} of {len(cells)} grid cells ({frac:.0%}) pass.  The "
            "timing is publishable only with the template stated, and the "
            "template is then part of the result.")
    return "NOT PUBLISHABLE", (
        f"Even with the edge shape known exactly sigma_t = {best:.1f} s, "
        f"outside the {threshold_s:.0f} s threshold.  Per-cycle timing is "
        "not supported by this cadence.")


# ===========================================================================
# 5.  Accretion states
# ===========================================================================

def phase_coverage(times_d, period_d: float, epoch_d: float,
                   n_bins: int = STATE_PHASE_BINS) -> float:
    """Fraction of the orbital phase circle a set of times touches at all.

    The gate on whether a night's median magnitude may be called a STATE.
    These polars vary by 0.65-1.7 mag around one orbit, which is the same
    size as a high-to-low state change, so a median over a third of an orbit
    is a phase measurement wearing a state label.
    """
    t = np.asarray(times_d, dtype=float)
    t = t[np.isfinite(t)]
    if t.size == 0:
        return 0.0
    ph = phase_of(t, period_d, epoch_d)
    return float(np.unique((ph * n_bins).astype(int)).size / n_bins)


def otsu_threshold(values, n_bins: int = 64) -> tuple[float, float]:
    """Threshold that best splits a distribution into two classes.

    Otsu's method: histogram the values, then choose the cut that maximises
    the BETWEEN-class variance ``w0 w1 (mu0 - mu1)^2``.  It is used here
    rather than a round number or an eyeballed gap because it is
    deterministic, it has no free parameters beyond the bin count, and it
    derives the threshold from the observed bimodality instead of importing
    one.  Returns ``(threshold, separability)`` where separability is the
    between-class variance divided by the total variance — a number near 1
    means two genuinely distinct populations, and a number near 0 means the
    method has cut a unimodal distribution in half and the resulting
    "classification" is a fiction.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 4 or np.ptp(v) <= 0:
        return float("nan"), float("nan")
    counts, edges = np.histogram(v, bins=int(n_bins))
    centres = 0.5 * (edges[:-1] + edges[1:])
    total = counts.sum()
    if total == 0:
        return float("nan"), float("nan")
    p = counts / total
    w0 = np.cumsum(p)
    w1 = 1.0 - w0
    mu = np.cumsum(p * centres)
    mu_t = mu[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        mu0 = mu / w0
        mu1 = (mu_t - mu) / w1
        between = w0 * w1 * (mu0 - mu1) ** 2
    between = np.where(np.isfinite(between), between, -np.inf)
    # When the two populations are cleanly separated, EVERY cut inside the
    # empty gap gives identical class weights and means, so the between-class
    # variance is a plateau and ``argmax`` would return its first bin —
    # pressing the threshold right up against the bright population instead
    # of putting it in the middle of the gap.  Take the plateau's midpoint,
    # which is the same answer when the maximum is unique and the sensible
    # one when it is not.
    peak = float(np.max(between))
    plateau = np.flatnonzero(between >= peak - 1e-12 * max(abs(peak), 1.0))
    i = int(round(0.5 * (plateau[0] + plateau[-1])))
    var_total = float(np.var(v))
    sep = float(between[i] / var_total) if var_total > 0 else float("nan")
    return float(centres[i]), sep


def bootstrap_threshold(values, n_boot: int = N_BOOT, n_bins: int = 64,
                        seed: int = SEED) -> float:
    """Standard deviation of the Otsu threshold under resampling.

    The threshold's own uncertainty is what defines the INTERMEDIATE class:
    a night whose magnitude sits within this of the cut has not been
    classified, and saying so is better than assigning it to whichever side
    the histogram binning happened to favour.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 6:
        return float("nan")
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(int(n_boot)):
        s = rng.choice(v, size=v.size, replace=True)
        thr, _sep = otsu_threshold(s, n_bins)
        if np.isfinite(thr):
            out.append(thr)
    if len(out) < 10:
        return float("nan")
    return float(np.std(np.asarray(out), ddof=1))


def classify_state(mag: float, threshold: float, half_width: float
                   ) -> str:
    """HIGH / LOW / INTERMEDIATE for one night's magnitude.

    Brighter (smaller magnitude) than ``threshold - half_width`` is HIGH,
    fainter than ``threshold + half_width`` is LOW, and the band between is
    INTERMEDIATE.  The band is the threshold's own bootstrap uncertainty,
    so "intermediate" means "this data cannot say", not "the star was in a
    physically intermediate state" — a distinction the report has to keep
    making because they are not the same claim.
    """
    if not (np.isfinite(mag) and np.isfinite(threshold)):
        return "UNKNOWN"
    hw = half_width if np.isfinite(half_width) else 0.0
    if mag < threshold - hw:
        return "HIGH"
    if mag > threshold + hw:
        return "LOW"
    return "INTERMEDIATE"


def duty_cycle(mags, censored, threshold: float) -> dict:
    """Fraction of epochs in the HIGH state, censored and uncensored.

    ``censored[i]`` True means ``mags[i]`` is an UPPER LIMIT — the target
    was at least that faint and possibly fainter.  Two numbers come back:

    ``naive``
        detections only.  This is the duty cycle of the epochs at which the
        target was bright enough to be detected, which is close to a
        tautology and is reported to show the size of the bias.
    ``with_limits``
        every censored epoch counted as LOW.  This is exact rather than
        approximate: an upper limit fainter than the threshold PROVES the
        epoch is on the low side, so no survival modelling is needed for
        the duty cycle itself.  Limits that fall on the BRIGHT side of the
        threshold prove nothing, and are counted separately as
        ``n_uninformative`` rather than being silently assigned.
    """
    m = np.asarray(mags, dtype=float)
    c = np.asarray(censored, dtype=bool)
    ok = np.isfinite(m)
    m, c = m[ok], c[ok]
    det = ~c
    n_det = int(det.sum())
    naive = float((m[det] < threshold).sum() / n_det) if n_det else float("nan")
    # A limit is informative only when it is FAINTER than the threshold.
    informative = c & (m > threshold)
    uninformative = c & ~informative
    n_total = int(det.sum() + informative.sum())
    n_high = int((m[det] < threshold).sum())
    with_lim = float(n_high / n_total) if n_total else float("nan")
    return {"naive": naive, "with_limits": with_lim,
            "n_detections": n_det, "n_high": n_high,
            "n_censored": int(c.sum()),
            "n_informative_limits": int(informative.sum()),
            "n_uninformative": int(uninformative.sum()),
            "n_used": n_total,
            "bias": (float(naive - with_lim)
                     if np.isfinite(naive) and np.isfinite(with_lim)
                     else float("nan"))}


# ===========================================================================
# 6.  Detrending discipline: GP + signal, jointly
# ===========================================================================

def matern32_cov(times_d, sigma: float, rho: float) -> np.ndarray:
    """Matern-3/2 covariance matrix, ``sigma^2 (1 + s) exp(-s)``,
    ``s = sqrt(3) |dt| / rho``.

    The pure-numpy reference the celerite2 path is checked against.  Dense
    and O(N^3), so it is only ever used on one night at a time (150 points)
    and in the unit tests — but it is what makes the fast path auditable
    rather than trusted.  Matern-3/2 rather than a squared exponential
    because the trend being modelled is atmospheric and detector drift,
    which is continuous but not smooth, and a squared-exponential kernel
    would insist it is infinitely differentiable and then interpolate
    straight through a real step.
    """
    t = np.asarray(times_d, dtype=float)
    dt = np.abs(t[:, None] - t[None, :])
    s = math.sqrt(3.0) * dt / max(float(rho), 1e-12)
    return float(sigma) ** 2 * (1.0 + s) * np.exp(-s)


def gp_log_likelihood_dense(times_d, resid, yerr, sigma: float, rho: float
                            ) -> float:
    """Gaussian-process log likelihood by dense Cholesky.

    ``-0.5 (r^T K^-1 r + log det K + N log 2pi)`` with
    ``K = matern32 + diag(yerr^2)``.  Returns -inf rather than raising when
    the matrix is not positive definite, so an optimiser can walk into a bad
    corner of parameter space and walk back out.
    """
    r = np.asarray(resid, dtype=float)
    e = np.asarray(yerr, dtype=float)
    k = matern32_cov(times_d, sigma, rho)
    k[np.diag_indices_from(k)] += e ** 2
    try:
        c = np.linalg.cholesky(k)
    except np.linalg.LinAlgError:
        return -np.inf
    alpha = np.linalg.solve(c.T, np.linalg.solve(c, r))
    logdet = 2.0 * float(np.sum(np.log(np.diag(c))))
    return float(-0.5 * (float(r @ alpha) + logdet + r.size * math.log(2 * math.pi)))


def sinusoid_design(times_d, freq_cd: float) -> np.ndarray:
    """``[1, cos(2 pi f t), sin(2 pi f t)]`` — the signal model's columns."""
    t = np.asarray(times_d, dtype=float)
    ang = 2.0 * np.pi * float(freq_cd) * t
    return np.column_stack([np.ones_like(t), np.cos(ang), np.sin(ang)])


def joint_gp_fit(times_d, y, yerr, freq_cd: float,
                 rho_grid_d=None, sigma_grid=None,
                 use_celerite: bool = True) -> dict:
    """Fit ``GP(trend) + sinusoid(freq)`` SIMULTANEOUSLY.

    The signal's three linear coefficients are solved exactly by generalised
    least squares inside the GP covariance at every hyper-parameter grid
    point, so only the two GP hyper-parameters (``sigma``, ``rho``) are
    searched.  A grid rather than an optimiser because a 2-D grid is
    deterministic, reproducible from the published grid definition, and
    cannot report a local minimum as the answer.

    Returns ``amplitude`` (the recovered semi-amplitude, the number the
    demonstration turns on), ``sigma``, ``rho``, ``loglike`` and
    ``backend``.  ``backend`` is ``"celerite2"`` when the fast path ran and
    ``"dense"`` otherwise; the two are checked against each other in the
    unit tests and by the CLI on the real series, because a fast path nobody
    compared is an assumption.

    WHY THIS AND NOT DETREND-THEN-SEARCH.  A trend model fitted to
    ``y`` alone has no way to know which of the variation is trend; it will
    happily absorb the signal, and the search that follows then measures
    what the smoother left behind.  Fitting both at once forces the trend to
    explain only what the sinusoid cannot, which is the entire difference —
    and :func:`detrend_suppression` measures how large that difference is.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(yerr, dtype=float)
    ok = np.isfinite(t) & np.isfinite(yv) & np.isfinite(e) & (e > 0)
    t, yv, e = t[ok], yv[ok], e[ok]
    nan = float("nan")
    if t.size < 10:
        return {"amplitude": nan, "sigma": nan, "rho": nan,
                "loglike": nan, "backend": "none", "n": int(t.size)}
    # REBASE THE TIME AXIS.  These are BJD values around 2.46e6 and the
    # correlation lengths being fitted are ~0.1 d, so the interesting
    # structure lives in the 8th significant figure of the input.  The dense
    # path survives that because it only ever forms differences, but
    # celerite2 keeps the absolute times in its factorisation and at some
    # hyper-parameters the factorisation silently failed — leaving a
    # half-built object whose `apply_inverse` raised on a missing internal
    # array.  With `quiet=True` that failure does not raise from `compute`,
    # so those grid points dropped out of the celerite2 search while
    # remaining in the dense one, and the two backends then optimised over
    # different grids and disagreed by 3.6% in amplitude.  Subtracting the
    # first timestamp fixes the conditioning and does not change the model:
    # the kernel depends only on time differences and the design matrix is
    # centred separately.
    t = t - t[0]
    span = float(t.max() - t.min())
    scatter = float(np.std(yv, ddof=1))
    if rho_grid_d is None:
        # From about one cadence step to about the whole block.  Wider than
        # the signal period at the top end ON PURPOSE: the demonstration is
        # only honest if the GP is ALLOWED to eat the signal and declines to.
        rho_grid_d = np.geomspace(max(span / 200.0, 1e-4), max(span, 1e-3), 12)
    if sigma_grid is None:
        sigma_grid = np.geomspace(max(scatter / 50.0, 1e-5),
                                  max(scatter * 3.0, 1e-4), 12)
    a = sinusoid_design(t - t.mean(), freq_cd)
    backend = "dense"
    gp = None
    if use_celerite:
        try:
            from celerite2 import GaussianProcess, terms      # noqa: PLC0415
            gp = GaussianProcess(
                terms.Matern32Term(sigma=1.0, rho=1.0,
                                   eps=CELERITE_MATERN_EPS), mean=0.0)
            backend = "celerite2"
        except Exception:                                     # noqa: BLE001
            gp = None
            backend = "dense"

    def _solve(sig: float, rho: float):
        """GLS for the three linear coefficients at fixed hyper-parameters,
        plus the marginal log likelihood of the residual."""
        if gp is not None:
            gp.kernel = terms.Matern32Term(sigma=sig, rho=rho,
                                           eps=CELERITE_MATERN_EPS)
            # The WHOLE celerite2 path is guarded, not just `compute`.
            # `quiet=True` makes a failed factorisation return normally with
            # the object left half-built, so the error surfaces later out of
            # `apply_inverse` as an AttributeError on a private array.  A
            # try/except around `compute` alone therefore catches nothing,
            # which is how a grid point can vanish from one backend's search
            # and not the other's.
            try:
                gp.compute(t, yerr=e, quiet=True)
                # K^-1 applied column by column: celerite2's apply_inverse
                # is the whole reason the fast path exists.
                kia = np.column_stack([gp.apply_inverse(a[:, j].copy())
                                       for j in range(a.shape[1])])
                kiy = gp.apply_inverse(yv.copy())
            except Exception:                                 # noqa: BLE001
                return None
            if not (np.isfinite(kia).all() and np.isfinite(kiy).all()):
                return None
        else:
            k = matern32_cov(t, sig, rho)
            k[np.diag_indices_from(k)] += e ** 2
            try:
                c = np.linalg.cholesky(k)
            except np.linalg.LinAlgError:
                return None
            kia = np.linalg.solve(c.T, np.linalg.solve(c, a))
            kiy = np.linalg.solve(c.T, np.linalg.solve(c, yv))
        m = a.T @ kia
        try:
            coef = np.linalg.solve(m, a.T @ kiy)
        except np.linalg.LinAlgError:
            return None
        r = yv - a @ coef
        if gp is not None:
            try:
                ll = float(gp.log_likelihood(r))
            except Exception:                                 # noqa: BLE001
                return None
        else:
            ll = gp_log_likelihood_dense(t, r, e, sig, rho)
        if not np.isfinite(ll):
            return None
        return coef, ll

    best = None
    for sig in np.asarray(sigma_grid, dtype=float):
        for rho in np.asarray(rho_grid_d, dtype=float):
            got = _solve(float(sig), float(rho))
            if got is None:
                continue
            coef, ll = got
            if np.isfinite(ll) and (best is None or ll > best[0]):
                best = (ll, coef, float(sig), float(rho))
    if best is None:
        return {"amplitude": nan, "sigma": nan, "rho": nan,
                "loglike": nan, "backend": backend, "n": int(t.size)}
    ll, coef, sig, rho = best
    return {"amplitude": float(math.hypot(coef[1], coef[2])),
            "sigma": sig, "rho": rho, "loglike": float(ll),
            "backend": backend, "n": int(t.size)}


def running_median_detrend(times_d, y, window_d: float) -> np.ndarray:
    """Subtract a running median of half-width ``window_d / 2``.

    THE WRONG WAY, implemented faithfully so the demonstration is a fair
    fight.  This is the standard "flatten it first" step, and it is standard
    because it works beautifully when the trend timescale is far longer than
    the signal period.  When it is not — and a 0.08 d orbit under a
    0.1 d smoothing window is not — the smoother tracks the signal and
    subtracts it from itself.
    """
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    half = float(window_d) / 2.0
    out = np.empty_like(yv)
    for i in range(yv.size):
        sel = np.abs(t - t[i]) <= half
        out[i] = yv[i] - np.median(yv[sel])
    return out


def detrend_suppression(times_d, yerr, freq_cd: float,
                        amplitude_mag: float, window_d: float,
                        trend: Optional[np.ndarray] = None,
                        seed: int = SEED, n_real: int = 24,
                        use_celerite: bool = True) -> dict:
    """Inject a known sinusoid, then recover it BOTH ways.

    Returns the median recovered semi-amplitude as a fraction of the
    injected one for

    * ``frac_detrend`` — running-median detrend, then a sinusoid fit;
    * ``frac_joint``   — GP and sinusoid fitted together.

    The injected signal is the ONLY periodic content, and both methods see
    exactly the same data, so any difference between the two fractions is
    caused by the ORDER of the operations and by nothing else.  That is the
    whole argument, and it is why the demonstration injects rather than
    arguing from a real light curve whose true amplitude nobody knows.
    """
    t = np.sort(np.asarray(times_d, dtype=float))
    e = np.asarray(yerr, dtype=float)
    if e.size == 1:
        e = np.full(t.shape, float(e))
    rng = np.random.default_rng(seed)
    if trend is None:
        trend = np.zeros_like(t)
    trend = np.asarray(trend, dtype=float)
    fd, fj = [], []
    for _ in range(int(n_real)):
        phase = rng.uniform(0.0, 1.0)
        sig = amplitude_mag * np.sin(2 * np.pi * (freq_cd * t + phase))
        yv = 16.0 + trend + sig + rng.normal(0.0, e)
        # (a) the wrong order
        flat = running_median_detrend(t, yv, window_d)
        a_det = _lsq_amplitude(t, flat, e, freq_cd)
        # (b) the right order
        joint = joint_gp_fit(t, yv, e, freq_cd, use_celerite=use_celerite)
        if np.isfinite(a_det):
            fd.append(a_det / amplitude_mag)
        if np.isfinite(joint["amplitude"]):
            fj.append(joint["amplitude"] / amplitude_mag)
    med = lambda v: float(np.median(v)) if v else float("nan")   # noqa: E731
    return {"frac_detrend": med(fd), "frac_joint": med(fj),
            "n_detrend": len(fd), "n_joint": len(fj),
            "window_d": float(window_d),
            "window_periods": float(window_d * freq_cd),
            "amplitude_in": float(amplitude_mag)}


def _lsq_amplitude(times_d, y, yerr, freq_cd: float) -> float:
    """Weighted least-squares semi-amplitude at a fixed frequency."""
    t = np.asarray(times_d, dtype=float)
    yv = np.asarray(y, dtype=float)
    e = np.asarray(yerr, dtype=float)
    a = sinusoid_design(t - t.mean(), freq_cd)
    w = 1.0 / e ** 2
    atw = a.T * w
    try:
        coef = np.linalg.solve(atw @ a, atw @ yv)
    except np.linalg.LinAlgError:
        return float("nan")
    return float(math.hypot(coef[1], coef[2]))
