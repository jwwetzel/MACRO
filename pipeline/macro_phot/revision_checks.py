"""CV-R9 and CV-R11 — the arithmetic of the reduction and clock checks.

WHY THIS MODULE EXISTS
----------------------
The plan review of 2026-10-03 found that the CV paper never described its
reduction (OA.E5), that the server flat-fields leave a sky ramp of a few
per cent across the 2025 Mode0 frames (TE.F5), that one camera "era" spans
several mechanical states of the instrument (TE, F-3), and that the
absolute clock of the season that carries most of ST LMi's timing has no
clock target of its own (F-8, era D).  Each finding is answered by a test,
and the tests need four small pieces of arithmetic, all here, all pure:

1.  **Differential light curves from apertures** (:func:`differential_mags`)
    — target (or check star) flux against the summed flux of the
    comparison stars measured on the same frame.  Used to photometer one
    night per era twice, from two different reductions of the same raw
    pixels, at identical positions and apertures, so that the ONLY thing
    that differs between the two light curves is the reduction.
2.  **Agreement of two light curves** (:func:`compare_light_curves`) — the
    scatter of their frame-by-frame difference after removing its median
    (a constant offset is a zero point and is irrelevant to differential
    photometry), as a plain RMS and as a robust MAD-sigma.
3.  **A ramp in x** (:func:`ramp_fit`) — a straight line through
    catalogue-tie residuals against normalised detector column, reported
    as the edge-to-edge change in per cent with a SCATTER-based error (the
    residuals' own dispersion, not their formal errors, sets the error bar,
    because catalogue errors and unrecognised variables dominate them).
4.  **A step between states** (:func:`state_offset`) and **a clock offset
    interpolated between verified eras** (:func:`interpolated_offset`).

Nothing here touches a file or a database.
``pipeline/scripts/run_cv_checks.py`` does the I/O and
``pipeline/tests/test_revision_checks.py`` checks every function on
synthetic data with a known answer.
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

#: 1.4826 x MAD is the standard deviation of a normal distribution.
MAD_TO_SIGMA = 1.4826

#: Magnitude -> fractional flux for small numbers: d(flux)/flux =
#: 0.4 ln(10) d(mag) = 0.921 d(mag).  Used to quote a ramp in per cent.
MAG_TO_FRAC = 0.4 * math.log(10.0)


def mad_sigma(x: Sequence[float]) -> float:
    """Robust standard deviation, 1.4826 x median absolute deviation."""
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return float(MAD_TO_SIGMA * np.median(np.abs(a - np.median(a))))


def differential_mags(star_flux: np.ndarray, comp_flux: np.ndarray
                      ) -> np.ndarray:
    """``-2.5 log10(star / sum(comps))`` per frame.

    ``star_flux`` has shape (n_frames,); ``comp_flux`` (n_frames, n_comp).
    A comparison star missing on a frame is NaN there and is left out of
    that frame's sum; a frame with no positive comparison flux, or a
    non-positive star flux, returns NaN rather than a number.
    """
    s = np.asarray(star_flux, dtype=float)
    c = np.asarray(comp_flux, dtype=float)
    if c.ndim == 1:
        c = c[:, None]
    tot = np.nansum(np.where(np.isfinite(c) & (c > 0), c, np.nan), axis=1)
    ok = np.isfinite(s) & (s > 0) & np.isfinite(tot) & (tot > 0)
    out = np.full(s.shape, np.nan)
    out[ok] = -2.5 * np.log10(s[ok] / tot[ok])
    return out


def compare_light_curves(a: Sequence[float], b: Sequence[float]) -> dict:
    """How well two light curves of the same frames agree.

    Returns ``n`` (frames where both are finite), ``offset`` (median of
    a - b: a zero point, reported but not counted as disagreement), ``rms``
    and ``mad_sigma`` of the difference about that median, and ``max_abs``
    (the largest single-frame disagreement).
    """
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    d = d[np.isfinite(d)]
    if d.size < 2:
        return {"n": int(d.size), "offset": float("nan"),
                "rms": float("nan"), "mad_sigma": float("nan"),
                "max_abs": float("nan")}
    off = float(np.median(d))
    r = d - off
    return {"n": int(d.size), "offset": off,
            "rms": float(np.sqrt(np.mean(r ** 2))),
            "mad_sigma": mad_sigma(r),
            "max_abs": float(np.max(np.abs(r)))}


def scatter_about_median(x: Sequence[float]) -> float:
    """RMS of a light curve about its own median (NaNs ignored)."""
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size < 2:
        return float("nan")
    return float(np.sqrt(np.mean((a - np.median(a)) ** 2)))


#: Outlier clip for the ramp fit, in robust sigmas about the current line.
#: Tie residuals carry a long tail (blends, variables, mis-matched catalogue
#: entries): in ST LMi's field the RMS is five times the MAD-sigma, and a
#: least-squares line through the tail measures the tail.
RAMP_CLIP_SIGMA = 4.0


def ramp_fit(x_norm: Sequence[float], resid_mag: Sequence[float],
             clip_sigma: float = RAMP_CLIP_SIGMA, max_iter: int = 10
             ) -> dict:
    """Straight line through tie residuals against normalised column.

    ``x_norm`` runs 0 (first column) to 1 (last); ``resid_mag`` is the
    tie residual as ``cv_cattie_star`` stores it, ensemble minus catalogue
    magnitude after the fitted zero point and colour term (any constant is
    absorbed by the intercept), so a POSITIVE residual is a star the
    reduced frames measure too FAINT.
    The slope is the change in magnitude from one edge of the detector to
    the other.  Its error is the ordinary least-squares error with the
    residuals' own scatter as sigma, because the formal catalogue errors
    understate what scatters these points (blends, variables, colour).

    Points further than ``clip_sigma`` robust sigmas from the line are
    clipped, iteratively, before the final fit.

    Returns ``n`` (kept), ``n_clipped``, ``slope_mag`` and
    ``slope_err_mag`` (edge to edge),
    ``ramp_pct`` and ``ramp_err_pct`` (the same as a flux change, per
    cent, POSITIVE when stars on the right-hand side are measured too
    FAINT, i.e. the flat leaves the right-hand side under-corrected),
    ``nsigma`` and ``scatter_mag``.
    """
    x = np.asarray(x_norm, dtype=float)
    y = np.asarray(resid_mag, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n_in = int(x.size)
    nan = float("nan")
    if n_in < 4 or np.ptp(x) <= 0:
        return {"n": n_in, "n_clipped": 0, "slope_mag": nan,
                "slope_err_mag": nan, "ramp_pct": nan, "ramp_err_pct": nan,
                "nsigma": nan, "scatter_mag": nan}
    # Iterative clip about the current line, with a MAD-based sigma so the
    # tail being clipped does not set its own threshold.
    keep = np.ones(n_in, dtype=bool)
    for _ in range(max_iter):
        b0, a0 = np.polyfit(x[keep], y[keep], 1)
        r0 = y - (a0 + b0 * x)
        sig = mad_sigma(r0[keep])
        new = np.abs(r0 - np.median(r0[keep])) <= clip_sigma * sig
        if new.sum() < 4 or np.array_equal(new, keep):
            break
        keep = new
    x, y = x[keep], y[keep]
    n = int(x.size)
    xm = x - x.mean()
    sxx = float(np.sum(xm ** 2))
    b = float(np.sum(xm * (y - y.mean())) / sxx)
    a = float(y.mean())
    res = y - (a + b * xm)
    s = float(np.sqrt(np.sum(res ** 2) / (n - 2)))
    eb = s / math.sqrt(sxx)
    # resid = ensemble - catalogue.  A star measured too faint has an
    # ensemble magnitude too LARGE and therefore a POSITIVE residual, so a
    # positive slope means the right-hand side reads too faint: report that
    # flux deficit as a positive ramp.
    return {"n": n, "n_clipped": n_in - n, "slope_mag": b,
            "slope_err_mag": eb, "ramp_pct": 100.0 * MAG_TO_FRAC * b,
            "ramp_err_pct": 100.0 * MAG_TO_FRAC * eb,
            "nsigma": (abs(b) / eb) if eb > 0 else nan,
            "scatter_mag": s}


def state_offset(mags_a: Sequence[float], mags_b: Sequence[float]) -> dict:
    """Step in one star's magnitude between two instrument states.

    ``mags_a``/``mags_b`` are the star's per-frame magnitudes in the two
    states.  Returns the difference of medians (b minus a) and a standard
    error from each state's robust scatter, so that a step can be compared
    with the scatter of the light curves it would contaminate.
    """
    a = np.asarray(mags_a, dtype=float)
    b = np.asarray(mags_b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size < 3 or b.size < 3:
        return {"n_a": int(a.size), "n_b": int(b.size),
                "step": float("nan"), "step_err": float("nan")}
    se = math.sqrt((1.2533 * mad_sigma(a)) ** 2 / a.size
                   + (1.2533 * mad_sigma(b)) ** 2 / b.size)
    return {"n_a": int(a.size), "n_b": int(b.size),
            "step": float(np.median(b) - np.median(a)), "step_err": se}


def interpolated_offset(cycle: Sequence[float], oc_s: Sequence[float],
                        sigma_s: Sequence[float], verified: Sequence[bool]
                        ) -> dict:
    """Offset of UNVERIFIED epochs from the line through VERIFIED ones.

    A clock offset that is constant within an instrument era shifts every
    epoch of that era by the same amount.  Where the absolute clock of an
    era cannot be tested directly (no clock target was observed in it),
    its offset relative to the eras that WERE tested is still measured by
    the timing series itself: fit a straight line (a linear ephemeris)
    through the epochs of the verified eras, by weighted least squares,
    and take the weighted mean residual of the unverified era's epochs
    about that line.  Its error combines the scatter of those residuals
    and the line's own uncertainty at their mean cycle.

    This assumes the timed feature keeps a constant period across the
    interval, which is exactly what the O-C tests; a period derivative at
    the paper's bound changes the result by far less than its error over
    one year, and the caller states that.

    Returns ``offset_s``, ``offset_err_s``, ``n_verified``,
    ``n_unverified``, ``chi2_verified`` and ``dof_verified``.
    """
    c = np.asarray(cycle, dtype=float)
    y = np.asarray(oc_s, dtype=float)
    s = np.asarray(sigma_s, dtype=float)
    v = np.asarray(verified, dtype=bool)
    nan = float("nan")
    if v.sum() < 3 or (~v).sum() < 1:
        return {"offset_s": nan, "offset_err_s": nan,
                "n_verified": int(v.sum()), "n_unverified": int((~v).sum()),
                "chi2_verified": nan, "dof_verified": 0}
    w = 1.0 / s[v] ** 2
    A = np.vstack([np.ones(v.sum()), c[v]]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    beta = cov @ (A.T @ (w * y[v]))
    chi2 = float(np.sum(w * (y[v] - A @ beta) ** 2))
    dof = int(v.sum() - 2)
    cu, yu, su = c[~v], y[~v], s[~v]
    wu = 1.0 / su ** 2
    resid = yu - (beta[0] + beta[1] * cu)
    off = float(np.sum(wu * resid) / np.sum(wu))
    cbar = float(np.sum(wu * cu) / np.sum(wu))
    g = np.array([1.0, cbar])
    line_var = float(g @ cov @ g)
    err = math.sqrt(1.0 / float(np.sum(wu)) + line_var)
    return {"offset_s": off, "offset_err_s": err,
            "n_verified": int(v.sum()), "n_unverified": int((~v).sum()),
            "chi2_verified": chi2, "dof_verified": dof}


def nights_in_ranges(nights: Sequence[str],
                     ranges: Sequence[tuple[str, str, str]]
                     ) -> list[Optional[str]]:
    """Label each night (ISO date string) with the range that contains it.

    ``ranges`` are ``(label, first_night, last_night)``, inclusive; a night
    in none of them is labelled ``None``.  ISO dates compare correctly as
    strings, which is why no date parsing is needed.
    """
    out: list[Optional[str]] = []
    for n in nights:
        hit = None
        for lab, lo, hi in ranges:
            if lo <= n <= hi:
                hit = lab
                break
        out.append(hit)
    return out
