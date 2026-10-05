"""macro_tcrb.p0 — the Phase-0 gates' decision logic (pure functions).

Each function here is a rule a Phase-0 table is built by; the CLI
(``pipeline/scripts/run_tcrb.py``) does the I/O.  Tests:
``pipeline/tests/test_tcrb.py``.

TCRB-P0-mech-epoch
    A calibration, trace prior or zero-order prediction may be applied to a
    frame only inside the frame's mechanical epoch.  Two epoch vocabularies
    exist: the formal F-3 table (``frame_mech_epoch``; splits at camera,
    flip, coarse rotation AND wheel-map changes) and the grism library's
    working table (``macro_grism.config.MECH_EPOCHS``; merges the
    wheel-relabelling boundaries a grism cannot feel, and adds the 2023
    slot-6 reversal).  A grism-library epoch is safe when every frame it
    holds sits in ONE hardware epoch — :func:`crosswalk_verdicts`.

TCRB-P0-temp-split
    :func:`temp_group` assigns each grism frame to the temperature group its
    dark treatment is stated for; :func:`flank_residual` is the
    flanking-band adequacy statistic.

TCRB-P0-filter-forensics / -zmag-provenance / -shutter-timing
    :func:`colour_slope`, :func:`lambda_eff_from_slope`,
    :func:`zmag_band_match`, :func:`exposure_offset_fit`.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Optional, Sequence

import numpy as np

# ===========================================================================
# Mechanical epochs
# ===========================================================================
#: Formal-boundary causes a grism cannot feel (as in macro_grism.config).
SOFT_CAUSES = {"wheel_map", "rotation_fine"}


def crosswalk_verdicts(pairs: Iterable[tuple[str, str]],
                       formal_cause: dict[str, str]) -> list[dict]:
    """Classify each grism-library epoch against the formal epochs.

    ``pairs``: (grism_epoch, formal_epoch) per frame.  ``formal_cause``:
    formal epoch -> the boundary cause that OPENED it.

    A grism epoch is ``one_to_one`` when its frames sit in one formal epoch;
    ``merged_soft`` when they span several formal epochs separated only by
    soft (wheel-label / fine-rotation) boundaries; ``CROSSES_HARDWARE``
    when any formal epoch after the first was opened by a camera, flip or
    coarse-rotation change — then a solution keyed by that grism epoch
    WOULD be applied across a hardware boundary.
    """
    members: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for g, f in pairs:
        members[g][f] += 1
    out = []
    for g, fs in sorted(members.items()):
        formal = sorted(fs)
        if len(formal) == 1:
            verdict = "one_to_one"
        else:
            # The earliest formal epoch's opening cause is irrelevant (it is
            # the grism epoch's own boundary); every later one must be soft.
            later = formal[1:]
            causes = set()
            for f in later:
                causes |= set((formal_cause.get(f) or "").split(","))
            verdict = ("merged_soft" if causes <= SOFT_CAUSES
                       else "CROSSES_HARDWARE")
        out.append({"grism_epoch": g, "formal_epochs": ",".join(formal),
                    "n_frames": int(sum(fs.values())), "verdict": verdict})
    return out


def master_allowed(science_epoch: Optional[str],
                   master_epoch: Optional[str]) -> bool:
    """A master calibration may be applied only inside its formal epoch."""
    return bool(science_epoch) and science_epoch == master_epoch


# ===========================================================================
# Temperature split
# ===========================================================================
#: Grism-frame temperature groups (CCD-TEMP, deg C).  The ASI was run at a
#: -10 C set point; the 2025 masters are all at -10 C.  "Cold" is within
#: 2.5 C of that; "warm" is anything above -5 C (the regulation failed or
#: was switched off); between is "intermediate".
COLD_MAX_C = -7.5
WARM_MIN_C = -5.0


def temp_group(ccd_temp: Optional[float]) -> str:
    if ccd_temp is None or not np.isfinite(ccd_temp):
        return "unknown"
    if ccd_temp <= COLD_MAX_C:
        return "cold"
    if ccd_temp >= WARM_MIN_C:
        return "warm"
    return "intermediate"


def flank_residual(strip: np.ndarray, half: int, gap: int, width: int,
                   centre: int) -> tuple[float, float]:
    """Flanking-band background residual on a source-free strip.

    ``strip`` is a 2-D cutout (rows = cross-dispersion, columns =
    dispersion) of blank sky.  A pseudo-aperture of ``2*half+1`` rows is
    centred on row ``centre``; the background is the mean of two flanking
    bands of ``width`` rows each, starting ``gap`` rows beyond the aperture
    edge on either side — the same geometry the extraction uses.  Returns
    (mean residual per aperture pixel, its standard error), both in ADU,
    computed column-wise and then averaged over columns.
    """
    a0, a1 = centre - half, centre + half + 1
    lo = slice(a0 - gap - width, a0 - gap)
    hi = slice(a1 + gap, a1 + gap + width)
    ap = strip[a0:a1].mean(axis=0)
    bg = 0.5 * (strip[lo].mean(axis=0) + strip[hi].mean(axis=0))
    d = ap - bg
    d = d[np.isfinite(d)]
    if d.size < 3:
        return float("nan"), float("nan")
    return float(d.mean()), float(d.std(ddof=1) / np.sqrt(d.size))


# ===========================================================================
# Filter forensics: effective wavelength from a colour slope
# ===========================================================================
def colour_slope(m_inst: np.ndarray, m_ref: np.ndarray, colour: np.ndarray,
                 clip: float = 3.0, n_iter: int = 5
                 ) -> tuple[float, float, float, int]:
    """Robust slope of (m_inst - m_ref) against colour.

    Returns (slope, slope_err, intercept, n_used).  A filter bluer than the
    reference band gives a positive slope against BP-RP (red stars are
    fainter in it than in the reference); redder gives negative.  Iterative
    sigma clipping on the MAD scale; the error is the OLS standard error of
    the surviving points.
    """
    y = np.asarray(m_inst, float) - np.asarray(m_ref, float)
    x = np.asarray(colour, float)
    ok = np.isfinite(x) & np.isfinite(y)
    for _ in range(n_iter):
        if ok.sum() < 5:
            return float("nan"), float("nan"), float("nan"), int(ok.sum())
        A = np.vstack([x[ok], np.ones(ok.sum())]).T
        coef, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
        r = y - (coef[0] * x + coef[1])
        mad = 1.4826 * np.median(np.abs(r[ok] - np.median(r[ok])))
        new = ok & (np.abs(r) <= clip * max(mad, 1e-3))
        if new.sum() == ok.sum():
            break
        ok = new
    A = np.vstack([x[ok], np.ones(ok.sum())]).T
    coef, res, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
    n = int(ok.sum())
    r = y[ok] - A @ coef
    s2 = float(r @ r) / max(n - 2, 1)
    cov = s2 * np.linalg.inv(A.T @ A)
    return float(coef[0]), float(np.sqrt(cov[0, 0])), float(coef[1]), n


def lambda_eff_from_slope(slope: float, known: Sequence[tuple[float, float]]
                          ) -> float:
    """Interpolate an effective wavelength from a colour slope.

    ``known``: (slope, lambda_eff) for filters of known identity measured
    the same way.  Slope is monotonic in wavelength over the optical, so a
    linear interpolation (extrapolation at the ends) in slope is used.
    """
    k = sorted(known)
    s = np.array([a for a, _ in k])
    lam = np.array([b for _, b in k])
    if not np.all(np.diff(s) > 0):
        # Monotonic is the premise; if it fails the mapping is not defined.
        return float("nan")
    if slope <= s[0]:
        i = 0
    elif slope >= s[-1]:
        i = len(s) - 2
    else:
        i = int(np.searchsorted(s, slope) - 1)
    f = (slope - s[i]) / (s[i + 1] - s[i])
    return float(lam[i] + f * (lam[i + 1] - lam[i]))


# ===========================================================================
# ZMAG provenance
# ===========================================================================
def zmag_band_match(m_inst_1s: np.ndarray, zmag: float,
                    bands: dict[str, np.ndarray]) -> dict[str, dict]:
    """How well does ``m_inst + ZMAG`` reproduce each catalogue band?

    ``m_inst_1s`` is -2.5 log10(counts / s) per matched star; PinPoint's
    ZMAG is defined so that m_inst + ZMAG = m_catalogue for the band it
    solved against.  For each candidate band, returns the median offset
    (catalogue - implied), its robust scatter and n.  The band PinPoint
    used is the one with the smallest scatter AND an offset consistent
    with an aperture correction (a few tenths at most).
    """
    implied = np.asarray(m_inst_1s, float) + zmag
    out = {}
    for name, cat in bands.items():
        d = np.asarray(cat, float) - implied
        d = d[np.isfinite(d)]
        if d.size < 5:
            out[name] = {"offset": float("nan"), "scatter": float("nan"),
                         "n": int(d.size)}
            continue
        med = float(np.median(d))
        out[name] = {"offset": med,
                     "scatter": float(1.4826 * np.median(np.abs(d - med))),
                     "n": int(d.size)}
    return out


# ===========================================================================
# Short-exposure timing
# ===========================================================================
def exposure_offset_fit(t_nom: np.ndarray, counts: np.ndarray,
                        sigma: np.ndarray) -> dict:
    """Fit counts = rate * (t_nom + delta): the exposure-time offset.

    Linear least squares in (rate, rate*delta).  A shutter/driver that
    delivers a constant extra (or missing) integration delta shows up as a
    non-zero intercept.  Returns rate, delta (s) and their 1-sigma errors
    (delta's from the covariance by first-order propagation), chi2, dof.
    """
    t = np.asarray(t_nom, float)
    c = np.asarray(counts, float)
    s = np.asarray(sigma, float)
    w = 1.0 / s ** 2
    A = np.vstack([t, np.ones_like(t)]).T
    Aw = A * np.sqrt(w)[:, None]
    cw = c * np.sqrt(w)
    coef, *_ = np.linalg.lstsq(Aw, cw, rcond=None)
    cov = np.linalg.inv(Aw.T @ Aw)
    rate, b = coef
    delta = b / rate
    # var(b/r) = (1/r^2) var b + (b^2/r^4) var r - 2 b/r^3 cov
    var_d = (cov[1, 1] / rate ** 2 + b ** 2 * cov[0, 0] / rate ** 4
             - 2 * b * cov[0, 1] / rate ** 3)
    r = (c - A @ coef) / s
    chi2 = float(r @ r)
    dof = int(t.size - 2)
    return {"rate": float(rate), "rate_err": float(np.sqrt(cov[0, 0])),
            "delta_s": float(delta), "delta_err_s": float(np.sqrt(var_d)),
            "chi2": chi2, "dof": dof}
