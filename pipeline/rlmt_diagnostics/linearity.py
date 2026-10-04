"""Pure S2 linearity logic: counts-vs-exptime from archival exposure ladders.

THE METHOD
----------
A linear detector doubles its counts when the exposure doubles.  An
exposure *ladder* — the same bright target shot at several exposure times
in one visit (the 2024-05-20 Vega BeStar sequence: 0.0001 s -> 0.1 s) —
turns that statement into a measurement: fit  flux = k * exptime  through
the rung medians and read the per-rung residual.  Rungs that fall LOW at
the top of the ladder are the linearity roll-off (or the hard ceiling);
scatter at the bottom is shutter/exposure-timing error.  No ladder for a
mode is itself a finding and is recorded as one.

Flux extraction must survive a grism (the Vega ladder is spectra, not
points) and needs no WCS: the measure is the background-subtracted sum
over the single brightest ``box`` x ``box`` window in the frame, located by
a coarse box-sum scan.  Pure numpy (integral images), no photometry
package, deterministic.
"""

from __future__ import annotations

import math
import re
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Photometry window (pixels).  Big enough to hold a saturated star's whole
#: profile or a grism trace segment; small enough that the background
#: estimate stays local.
LADDER_BOX = 96

#: Coarse scan stride when locating the brightest window (the exact peak is
#: then refined by centering the box on the brightest coarse cell).
LADDER_SCAN_STRIDE = 32

#: A ladder needs at least this many distinct exposure-time rungs.
MIN_RUNGS = 3

#: ... and at least this many frames per rung to median away transients.
MIN_FRAMES_PER_RUNG = 2

#: Filename exposure-time token: '0p0001s' -> 0.0001, '10s' -> 10.0.  The
#: header EXPTIME rounds sub-millisecond exposures to 0.0 (observed on the
#: Vega ladder), so the filename is the better record at the short end.
_EXPTIME_TOKEN_RE = re.compile(r"(?:^|_)(\d+(?:p\d+)?)s(?:_|$|\.)")


def parse_exptime_token(basename: str) -> Optional[float]:
    """Exposure time encoded in a filename, or None.

    ``kaf_Vega_0p0001s_lrg_0.fts.fz`` -> 0.0001;
    ``mjc_HD_20134_hrg_64s_...`` -> 64.0.  The 'p' is the decimal point
    (filenames cannot carry '.').  Returns None when no token matches —
    callers then fall back to the header EXPTIME.
    """
    m = _EXPTIME_TOKEN_RE.search(basename)
    if not m:
        return None
    return float(m.group(1).replace("p", "."))


def effective_exptime(header_exptime: Optional[float],
                      basename: str) -> Optional[float]:
    """The exposure time to trust for ladder work.

    The filename token wins whenever the header is missing, non-positive,
    or disagrees with the token by more than 20% (the header driver rounds
    0.0001 s to 0.0 — the observed failure).  Otherwise the header stands.
    """
    tok = parse_exptime_token(basename)
    h = header_exptime
    if h is None or not math.isfinite(h) or h <= 0:
        return tok if tok is not None else h
    if tok is not None and tok > 0 and abs(h - tok) > 0.2 * max(h, tok):
        return tok
    return float(h)


def brightest_box_flux(img: np.ndarray, box: int = LADDER_BOX,
                       stride: int = LADDER_SCAN_STRIDE) -> dict:
    """Background-subtracted flux in the brightest box-window of a frame.

    Integral-image box sums on a stride grid find the brightest window;
    the background is the frame median (dominant sky/bias level) times the
    box area.  Returns the flux, the window's corner, the frame's peak
    pixel inside the window, and the frame median — everything the ladder
    fit and the saturation cross-check need.
    """
    a = np.asarray(img, dtype=np.float64)
    ny, nx = a.shape
    b = min(box, ny, nx)
    # Integral image: S[i, j] = sum of a[:i, :j]; box sums in O(1) each.
    S = np.zeros((ny + 1, nx + 1))
    np.cumsum(np.cumsum(a, axis=0), axis=1, out=S[1:, 1:])
    ys = np.arange(0, ny - b + 1, stride)
    xs = np.arange(0, nx - b + 1, stride)
    sums = (S[np.ix_(ys + b, xs + b)] - S[np.ix_(ys, xs + b)]
            - S[np.ix_(ys + b, xs)] + S[np.ix_(ys, xs)])
    iy, ix = np.unravel_index(int(np.argmax(sums)), sums.shape)
    y0, x0 = int(ys[iy]), int(xs[ix])
    med = float(np.median(a))
    window = a[y0:y0 + b, x0:x0 + b]
    return {
        "flux": float(window.sum() - med * b * b),
        "y0": y0, "x0": x0, "box": b,
        "peak_adu": float(window.max()),
        "sky_med": med,
    }


def fit_ladder(exptimes: Sequence[float], fluxes: Sequence[float]) -> Optional[dict]:
    """Fit  flux = k * exptime  and report per-rung residuals.

    ``exptimes``/``fluxes`` are ONE value per rung (the caller medians its
    frames first).  The rate k is the median of flux/exptime over rungs —
    robust, and immune to a single rolled-off top rung dragging the line.
    Residuals are percentages: 100 * (flux / (k * t) - 1).  Returns None
    for fewer than :data:`MIN_RUNGS` usable rungs.
    """
    t = np.asarray(exptimes, dtype=np.float64)
    f = np.asarray(fluxes, dtype=np.float64)
    ok = np.isfinite(t) & np.isfinite(f) & (t > 0)
    t, f = t[ok], f[ok]
    if t.size < MIN_RUNGS:
        return None
    order = np.argsort(t)
    t, f = t[order], f[order]
    rates = f / t
    k = float(np.median(rates))
    if k <= 0:
        return None
    resid_pct = 100.0 * (f / (k * t) - 1.0)
    return {
        "rate_adu_per_s": k,
        "exptimes": t.tolist(),
        "fluxes": f.tolist(),
        "resid_pct": resid_pct.tolist(),
        "max_abs_resid_pct": float(np.max(np.abs(resid_pct))),
        "n_rungs": int(t.size),
    }


def fair_ladder_order(candidates: Sequence[tuple],
                      mode_of, quality) -> list[tuple]:
    """Order ladder candidates so every MODE is tried before any repeats.

    Parameters
    ----------
    candidates
        Opaque rows; this function never looks inside them.
    mode_of
        ``row -> mode label``.  Rows sharing a label compete with each
        other and with nobody else.
    quality
        ``row -> sort key``, best first, applied WITHIN a mode.

    Returns
    -------
    list
        Every mode's best candidate, then every mode's second, and so on —
        modes visited in sorted label order so the result is deterministic.

    WHY THIS EXISTS.  Ranking candidates globally by quality is the obvious
    thing and it is wrong here: the archive's frame counts are wildly
    uneven (Mode0 alone offers hundreds of candidate ladders, Low Gain and
    5 MHz offer one or none), so a global ranking spends every processing
    slot on the two richest modes and reports "no archival linearity
    constraint" for the sparse ones — a conclusion about the SCHEDULER
    dressed up as a conclusion about the archive.  Round-robin guarantees
    each mode is tried; a mode that then yields nothing has genuinely
    yielded nothing, and that is a finding worth reporting.
    """
    by_mode: dict[str, list] = {}
    for row in candidates:
        by_mode.setdefault(mode_of(row), []).append(row)
    for rows in by_mode.values():
        rows.sort(key=quality)
    out: list[tuple] = []
    for slot in range(max((len(v) for v in by_mode.values()), default=0)):
        for label in sorted(by_mode):
            if slot < len(by_mode[label]):
                out.append(by_mode[label][slot])
    return out


def group_ladders(rows: Sequence[tuple],
                  min_rungs: int = MIN_RUNGS,
                  min_frames: int = MIN_FRAMES_PER_RUNG) -> list[dict]:
    """Find exposure ladders in manifest rows.

    ``rows`` = (night, target_key, readoutm, exptime_bin, n_frames) tuples
    (one per group, from SQL).  A ladder = one (night, target, mode) with
    >= ``min_rungs`` distinct exposure bins each holding >= ``min_frames``
    frames.  Returns one dict per ladder, rungs sorted by exposure —
    ready for the campaign script to fetch pixels for.
    """
    from collections import defaultdict
    sets: dict[tuple, dict[float, int]] = defaultdict(dict)
    for night, target, mode, ebin, n in rows:
        if ebin is None or ebin <= 0 or target is None:
            continue
        sets[(night, target, mode)][float(ebin)] = int(n)
    out = []
    for (night, target, mode), rungs in sets.items():
        good = {t: n for t, n in rungs.items() if n >= min_frames}
        if len(good) >= min_rungs:
            out.append({"night": night, "target_key": target, "mode": mode,
                        "rungs": sorted(good), "n_frames": sum(good.values())})
    # Deterministic order: richest ladders first, then by night/target.
    out.sort(key=lambda d: (-len(d["rungs"]), -d["n_frames"],
                            d["night"] or "", d["target_key"]))
    return out


# ==========================================================================
# Peak-resolved differential linearity (committee findings DE.F3, F-5)
# ==========================================================================
#
# WHY THE LADDERS ABOVE WERE NOT ENOUGH
# -------------------------------------
# An exposure ladder asks "does flux double when the exposure doubles?" of
# ONE source, and its answer is contaminated by everything else that
# changed between the rungs (transparency, seeing, shutter timing).  The
# review found that the per-mode number S2 v1.2 published from them was the
# BEST ladder per mode — for Mode0 one whose brightest rung peaked at 1% of
# full scale — so it constrained nothing near any saturation veto, while
# four different caps (70%, 80%, 80%, 92%) were in force downstream.
#
# The measurement a cap needs is DIFFERENTIAL and PEAK-RESOLVED: many stars
# of different brightness, measured in the same frames, so that whatever
# changed between frames is common to all of them and cancels in the
# ensemble — leaving, for each star in each frame, a fractional deviation
# from linear response that can be plotted against that star's own peak
# pixel.  A linear detector gives zero at every peak level.  Roll-off
# toward the ceiling shows up as a NEGATIVE deviation growing with peak.
#
# The same arithmetic serves three kinds of data:
#   * an exposure ladder on one field (rate = flux / exptime);
#   * one field in two readout modes at one exposure (StackPro vs High
#     Gain: the 16-read sum never approaches its sub-read clip, so it is
#     the linear reference for the single read);
#   * a time series at fixed exposure, where seeing moves each star's peak
#     up and down from frame to frame (the CV comparison stars).

#: Stars (in frames) whose peak sits below this fraction of full scale
#: define "linear" — the reference regime every deviation is measured
#: against.  0.30 leaves two thirds of the scale to be tested and is far
#: below every cap anyone has proposed (the lowest is 70%).
REF_PEAK_FRACTION = 0.30

#: Peak-fraction bin edges of the published deviation curve.  Finer above
#: 0.8, where the caps in dispute live (0.70, 0.80, 0.92).
PEAK_BIN_EDGES = (0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                  0.85, 0.90, 0.95, 1.00, 1.10)

#: The cap criterion adopted by the committee (F-5): response deviation
#: below this many percent up to the cap.
CAP_CRITERION_PCT = 1.0

#: A bin is "measured" only with at least this many distinct stars (or,
#: for a fixed-exposure series, distinct star-night groups).
MIN_GROUPS_PER_BIN = 3

#: Reference stars need at least this signal-to-noise (flux / floor error).
REF_MIN_SNR = 30.0

#: A frame's factor z is trusted only when at least this many reference
#: stars determine it.  With one reference star (the 32 s Albireo frame in
#: the first pass) the "frame factor" is that star's noise, and every
#: bright star in the frame inherits it as a spurious, peak-independent
#: 5% "deviation".  Frames below the minimum yield no deviations at all.
MIN_REF_PER_FRAME = 10

#: Normalisation iterations (frame factors <-> star base rates): an upper
#: bound, with :data:`NORM_TOL` the stopping rule.  The alternating solve
#: converges SLOWLY when the reference sets of different frames barely
#: overlap (a long exposure's reference stars are the faint ones, a short
#: exposure's are everything) — six iterations, the first version's fixed
#: count, left the frame factors 1% from their fixed point and that 1%
#: appeared as a spurious "non-linearity" in the null-injection test.
NORM_ITERS = 2000

#: Convergence: stop when no frame factor moves by more than this
#: (relative) in one iteration.
NORM_TOL = 1e-9


def peak_fraction(peak_raw, bias: float, ceiling: float):
    """Raw peak as a fraction of the usable scale: (peak - bias)/(ceil - bias)."""
    return (np.asarray(peak_raw, dtype=np.float64) - bias) / (ceiling - bias)


def ensemble_deviation(flux: np.ndarray, flux_err: np.ndarray,
                       exptime: np.ndarray, peak_frac: np.ndarray,
                       ref_max: float = REF_PEAK_FRACTION,
                       ref_min_snr: float = REF_MIN_SNR,
                       frame_is_reference: Optional[np.ndarray] = None,
                       min_ref_per_frame: int = None,
                       ) -> dict:
    """Fractional deviation from linear response, star by star, frame by frame.

    Parameters
    ----------
    flux, flux_err, peak_frac
        Arrays of shape (n_stars, n_frames); NaN where a star was not
        measured in a frame.
    exptime
        Length n_frames.  Use all ones for a fixed-exposure series.
    ref_max
        Measurements with ``peak_frac < ref_max`` (and EXPECTED S/N above
        ``ref_min_snr`` — see the selection note in the code) are the
        linear reference.
    frame_is_reference
        Optional boolean per frame: when given, ONLY measurements in these
        frames may be reference (the StackPro frames, when StackPro is the
        linear reference for High Gain), in addition to the peak test.

    Model
    -----
    ``rate[i, j] = flux[i, j] / exptime[j] = base[i] * z[j] * (1 + dev[i, j])``
    with ``base[i]`` the star's true rate and ``z[j]`` everything common to
    frame j (transparency, aperture loss, exposure-time error, the
    mode-to-mode scale).  ``base`` and ``z`` are solved by alternating
    weighted least squares over REFERENCE measurements only; ``dev`` is
    then evaluated for every measurement.

    Returns
    -------
    dict with ``dev`` and ``dev_err`` (same shape; fractional), ``is_ref``,
    ``base`` (n_stars; NaN for a star with no reference measurement — it
    cannot be anchored and its ``dev`` is NaN), ``z`` and ``z_err``
    (n_frames; the frame factor and its standard error from the reference
    stars' scatter — the common-mode error of every ``dev`` in that
    frame), and ``n_ref`` per frame.

    What this does NOT do: it cannot see a non-linearity that is the same
    at every signal level (a pure scale error is absorbed in ``z``), and
    the reference measurements have median ``dev`` 0 by construction —
    only the measurements ABOVE ``ref_max`` are a test.
    """
    if min_ref_per_frame is None:
        min_ref_per_frame = MIN_REF_PER_FRAME
    f = np.asarray(flux, dtype=np.float64)
    e = np.asarray(flux_err, dtype=np.float64)
    pf = np.asarray(peak_frac, dtype=np.float64)
    t = np.asarray(exptime, dtype=np.float64)
    n_stars, n_frames = f.shape
    with np.errstate(invalid="ignore", divide="ignore"):
        rate = f / t[None, :]
        rate_err = e / t[None, :]
        usable = (np.isfinite(rate) & np.isfinite(rate_err) & (rate_err > 0)
                  & np.isfinite(pf))
        low = usable & (pf < ref_max)
    if frame_is_reference is not None:
        low &= np.asarray(frame_is_reference, dtype=bool)[None, :]
    # SELECTION WITHOUT SELECTION BIAS.  The obvious reference cut — keep a
    # measurement if ITS OWN flux/error exceeds the threshold — is biased:
    # near the threshold a measurement passes only when noise pushed it UP,
    # so faint stars' short exposures enter the reference high, their base
    # rates come out high, and every brighter measurement then reads LOW by
    # 1-2% with no non-linearity present (found by the null-injection test
    # on the first version of this function).  The cut is therefore made
    # on the EXPECTED signal-to-noise, base * z * t / err, which does not
    # know which way the noise went; and base/z are inverse-variance
    # weighted MEANS (weights from the error floor, never from the measured
    # flux), which are unbiased where a median of a skewed selected sample
    # is not.
    w_all = np.where(low, 1.0 / np.square(rate_err), 0.0)
    z = np.ones(n_frames)
    is_ref = low.copy()
    base = np.full(n_stars, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"), \
            _quiet_nan_warnings():
        r0 = np.where(low, rate, 0.0)
        z_prev = z.copy()
        for _ in range(NORM_ITERS):
            w = np.where(is_ref, w_all, 0.0)
            # base_i = sum_j w z r / sum_j w z^2   (least squares for
            # r_ij = base_i z_j with known z).
            zz = np.where(np.isfinite(z), z, 0.0)
            num = (w * zz[None, :] * r0).sum(axis=1)
            den = (w * zz[None, :] ** 2).sum(axis=1)
            base = np.where(den > 0, num / den, np.nan)
            bb = np.where(np.isfinite(base), base, 0.0)
            num = (w * bb[:, None] * r0).sum(axis=0)
            den = (w * bb[:, None] ** 2).sum(axis=0)
            z = np.where(den > 0, num / den, np.nan)
            # The gauge (base -> base*c, z -> z/c) is fixed by the median
            # frame factor being 1, so base is a rate in the units of rate.
            zmed = np.nanmedian(z)
            if np.isfinite(zmed) and zmed > 0:
                z = z / zmed
                base = base * zmed
            expected_snr = (base[:, None] * z[None, :]) / rate_err
            is_ref_new = low & (expected_snr >= ref_min_snr)
            moved = np.nanmax(np.abs(z / z_prev - 1.0)) if _ else np.inf
            same_set = bool(np.array_equal(is_ref_new, is_ref))
            is_ref = is_ref_new
            z_prev = z.copy()
            if same_set and moved < NORM_TOL:
                break
        model = base[:, None] * z[None, :]
        ratio = np.where(is_ref, rate / model, np.nan)
        n_ref = np.sum(np.isfinite(ratio), axis=0)
        mad = MAD_TO_SIGMA_LIN * np.nanmedian(np.abs(ratio - 1.0), axis=0)
        z_err = np.where(n_ref > 1, mad / np.sqrt(np.maximum(n_ref, 1)),
                         np.nan)
        # A frame with too few reference stars has no trustworthy factor.
        z = np.where(n_ref >= min_ref_per_frame, z, np.nan)
        z_err = np.where(n_ref >= min_ref_per_frame, z_err, np.nan)
        # A star with no reference measurement cannot be anchored.
        anchored = np.any(is_ref, axis=1)
        base = np.where(anchored, base, np.nan)
        dev = rate / (base[:, None] * z[None, :]) - 1.0
        dev_err = np.abs(rate_err / (base[:, None] * z[None, :]))
    return {"dev": dev, "dev_err": dev_err, "is_ref": is_ref, "base": base,
            "z": z, "z_err": z_err, "n_ref": n_ref.astype(int)}


#: MAD -> sigma (Gaussian).
MAD_TO_SIGMA_LIN = 1.4826


class _quiet_nan_warnings:
    """Context manager silencing numpy's all-NaN-slice RuntimeWarnings.

    An all-NaN slice (a star never in the reference regime, a frame with no
    reference star) is an expected, handled outcome here — its result is
    NaN and is reported as "not anchored" — so the warning is noise.
    """

    def __enter__(self):
        import warnings
        self._cm = warnings.catch_warnings()
        self._cm.__enter__()
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return self

    def __exit__(self, *exc):
        return self._cm.__exit__(*exc)


def deviation_curve(peak_frac: Sequence[float], dev: Sequence[float],
                    group: Sequence, edges: Sequence[float] = PEAK_BIN_EDGES,
                    min_groups: int = MIN_GROUPS_PER_BIN) -> list[dict]:
    """Bin fractional deviations on peak fraction, with scatter-based errors.

    ``group`` labels the statistically independent unit each point belongs
    to (a star; or a star-night for a long series): points sharing a group
    share that star's base-rate error, so they are first averaged WITHIN
    the group and the bin statistics are then taken ACROSS groups.  Each
    bin returns:

    * ``lo``, ``hi``         — bin edges (peak fraction);
    * ``peak_frac``          — mean peak fraction of the points in it;
    * ``dev_pct``            — median over groups of the group means (%);
    * ``dev_err_pct``        — 1.2533 x MAD-sigma / sqrt(n_groups): the
      standard error of a median, from the measured group-to-group scatter
      (no error bar here is a propagated formal error);
    * ``n_points``, ``n_groups``;
    * ``measured``           — False when fewer than ``min_groups`` groups
      populate the bin (the value is then reported but must not be used).
    """
    pf = np.asarray(peak_frac, dtype=np.float64)
    dv = np.asarray(dev, dtype=np.float64)
    gr = np.asarray(group)
    ok = np.isfinite(pf) & np.isfinite(dv)
    pf, dv, gr = pf[ok], dv[ok], gr[ok]
    out: list[dict] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (pf >= lo) & (pf < hi)
        if not sel.any():
            continue
        means = np.array([dv[sel & (gr == g)].mean()
                          for g in np.unique(gr[sel])])
        med = float(np.median(means))
        if means.size > 1:
            mad = MAD_TO_SIGMA_LIN * float(np.median(np.abs(means - med)))
            # A MAD of zero from a handful of groups is not zero error.
            spread = mad if mad > 0 else float(means.std(ddof=1))
            err = 1.2533 * spread / math.sqrt(means.size)
        else:
            err = float("nan")
        out.append({"lo": float(lo), "hi": float(hi),
                    "peak_frac": float(pf[sel].mean()),
                    "dev_pct": 100.0 * med, "dev_err_pct": 100.0 * err,
                    "n_points": int(sel.sum()), "n_groups": int(means.size),
                    "measured": bool(means.size >= min_groups)})
    return out


def recommend_cap(curve: Sequence[dict],
                  criterion_pct: float = CAP_CRITERION_PCT,
                  ref_max: float = REF_PEAK_FRACTION) -> dict:
    """The linearity cap a deviation curve supports.

    THE RULE (fixed here, before any curve was looked at).  Walk the
    MEASURED bins upward from the reference regime.  The cap is the upper
    edge of the last bin for which this bin and every measured bin below
    it have ``|dev| <= criterion``.  Walking stops at the first measured
    bin that violates the criterion, or at the first gap (an unmeasured
    bin between measured ones is not assumed linear).

    Returns
    -------
    dict
        * ``cap_fraction`` — the cap as a fraction of full scale (None
          when not even the first bin above the reference is measured);
        * ``limited_by``   — ``"nonlinearity"`` (a measured bin above the
          cap violates the criterion), ``"data"`` (the measurements simply
          stop: the cap is "linear as far as it was measured", NOT "linear
          up to here and non-linear beyond") or ``"none"`` (no cap);
        * ``worst_dev_pct`` — the largest |dev| among the bins below the
          cap, signed;
        * ``precision_pct`` — the largest 1-sigma error among those bins:
          when this is not well below the criterion the cap is supported
          only to that precision, and the report says so;
        * ``first_bad``     — the violating bin (dict) or None;
        * ``slope_pct_per_scale`` and ``slope_err`` — weighted straight-
          line slope of dev (%) against peak fraction over the bins up to
          the cap: the "residual-vs-peak slope" of the committee's
          acceptance criterion, in percent per full scale.
    """
    bins = sorted((b for b in curve), key=lambda b: b["lo"])
    cap, worst, prec, first_bad, limited = None, 0.0, 0.0, None, "none"
    used: list[dict] = []
    started = False
    for b in bins:
        if b["hi"] <= ref_max and not started:
            # Reference-regime bins are zero by construction; they anchor
            # the slope fit but cannot fail.
            if b["measured"]:
                used.append(b)
            continue
        started = True
        if not b["measured"]:
            limited = "data"
            break
        if abs(b["dev_pct"]) > criterion_pct:
            first_bad, limited = b, "nonlinearity"
            break
        cap = b["hi"]
        used.append(b)
        if abs(b["dev_pct"]) > abs(worst):
            worst = b["dev_pct"]
        if np.isfinite(b["dev_err_pct"]):
            prec = max(prec, b["dev_err_pct"])
        limited = "data"
    slope, slope_err = None, None
    pts = [b for b in used if np.isfinite(b["dev_err_pct"])
           and b["dev_err_pct"] > 0]
    if len(pts) >= 3:
        x = np.array([b["peak_frac"] for b in pts])
        y = np.array([b["dev_pct"] for b in pts])
        w = 1.0 / np.array([b["dev_err_pct"] for b in pts]) ** 2
        xm = (w * x).sum() / w.sum()
        sxx = (w * (x - xm) ** 2).sum()
        if sxx > 0:
            slope = float((w * (x - xm) * y).sum() / sxx)
            slope_err = float(1.0 / math.sqrt(sxx))
    return {"cap_fraction": cap, "limited_by": limited if cap is not None
            or first_bad is not None else "none",
            "worst_dev_pct": float(worst), "precision_pct": float(prec),
            "first_bad": first_bad, "slope_pct_per_scale": slope,
            "slope_err": slope_err, "n_bins_used": len(used)}


def injection_bias(flux_err: np.ndarray, exptime: np.ndarray,
                   base: np.ndarray, z: np.ndarray,
                   peak_per_flux: np.ndarray, injected,
                   n_trials: int = 200, seed: int = 20261003,
                   edges: Sequence[float] = PEAK_BIN_EDGES,
                   frame_is_reference: Optional[np.ndarray] = None,
                   peak_noise: float = 0.02) -> list[dict]:
    """Signed recovered-minus-injected bias of the deviation curve.

    Standing rule 3: an injection test reports the SIGNED bias per cell,
    never a median of absolute values.  The simulation re-uses the REAL
    geometry of a data set — its stars' base rates, its frame factors, its
    per-measurement error floors, its peak-to-flux ratios (seeing) — and
    replaces only the fluxes, by ``base * z * t * (1 + injected(pf))`` plus
    Gaussian noise of the recorded size.  ``injected`` maps true peak
    fraction to fractional deviation (use ``lambda pf: 0 * pf`` for the
    null test).  The measured peak fraction is perturbed by ``peak_noise``
    (relative) so bin migration is part of what is tested.

    Returns one dict per bin: ``lo``, ``hi``, ``injected_pct`` (mean
    injected deviation of the points that landed in the bin),
    ``recovered_pct`` (mean over trials of the bin's curve value),
    ``bias_pct`` = recovered - injected (SIGNED), ``bias_err_pct`` (its
    standard error over trials) and ``n_trials`` in which the bin was
    measured.
    """
    rng = np.random.default_rng(seed)
    e = np.asarray(flux_err, dtype=np.float64)
    t = np.asarray(exptime, dtype=np.float64)
    b = np.asarray(base, dtype=np.float64)
    zz = np.asarray(z, dtype=np.float64)
    ppf = np.asarray(peak_per_flux, dtype=np.float64)
    true_flux = b[:, None] * zz[None, :] * t[None, :]
    pf_true = true_flux * ppf
    inj = np.asarray(injected(pf_true), dtype=np.float64)
    group = np.repeat(np.arange(e.shape[0])[:, None], e.shape[1], axis=1)
    rec: dict[tuple, list[float]] = {}
    inj_in_bin: dict[tuple, list[float]] = {}
    for _ in range(n_trials):
        flux = true_flux * (1.0 + inj) + rng.normal(0.0, 1.0, e.shape) * e
        pf_meas = pf_true * (1.0 + inj) * (
            1.0 + peak_noise * rng.normal(0.0, 1.0, e.shape))
        res = ensemble_deviation(flux, e, t, pf_meas,
                                 frame_is_reference=frame_is_reference)
        curve = deviation_curve(pf_meas.ravel(), res["dev"].ravel(),
                                group.ravel(), edges=edges)
        for c in curve:
            if not c["measured"]:
                continue
            key = (c["lo"], c["hi"])
            rec.setdefault(key, []).append(c["dev_pct"])
            sel = (np.isfinite(res["dev"]) & (pf_meas >= c["lo"])
                   & (pf_meas < c["hi"]))
            inj_in_bin.setdefault(key, []).append(
                100.0 * float(np.median(inj[sel])))
    out = []
    for key in sorted(rec):
        r = np.array(rec[key])
        i = np.array(inj_in_bin[key])
        d = r - i
        out.append({"lo": key[0], "hi": key[1],
                    "injected_pct": float(i.mean()),
                    "recovered_pct": float(r.mean()),
                    "bias_pct": float(d.mean()),
                    "bias_err_pct": float(d.std(ddof=1) / math.sqrt(d.size))
                    if d.size > 1 else float("nan"),
                    "n_trials": int(d.size)})
    return out
