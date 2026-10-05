"""Wavelength scale of a science frame: FIXED dispersion, per-frame zero
point (committee rulings U2 / G-1; findings OA.E1, PH.P7, TE.F3).

WHAT THIS MODULE USED TO DO, AND WHY IT WAS WITHDRAWN
-----------------------------------------------------
Version 1 solved a wavelength scale per frame from the frame itself: find
Halpha, find a PAIR of telluric dips whose offsets sit in the O2 A/B
lever ratio (3.387), read the dispersion off the pair.  Three reviewers
independently showed the stored results were not measurements of a
physical quantity: accepted T CrB frames carried dispersions from -1.8 to
+2.0 A/px on one grism, bimodal at 0.46 and 1.59 — the ratio of the two
modes being exactly 3.39.  On an M giant, TiO band heads supply decoy
dips at every spacing, and the pair test is satisfied by (TiO 6651, O2-B)
as readily as by (O2-B, O2-A).  A grism at a fixed distance from a
detector has ONE dispersion and ONE sign until someone moves the camera.

WHAT IT DOES NOW
----------------
* The dispersion comes from ``g_dispersion``: one row per (grism,
  mechanical epoch), solved on hot stars by ``run_g_dispersion.py``.
  The scale is a polynomial in DETECTOR position (an optical distortion
  of the grism-camera pair, the same for every star); this module never
  fits a dispersion.
* Per frame it measures exactly one number: the ZERO POINT — the pixel of
  Halpha.  For an emission-line target that is the centroid of the Halpha
  emission; for an absorption-line star, of the absorption core.
* It then performs one CHECK, which is not allowed to change anything:
  the blue edge of the telluric O2 B band must sit at the separation from
  Halpha that the fixed solution predicts.  The measured separation is
  stored per frame; its constancy across a series is the acceptance test
  of the whole scheme (G-1: constant to 1-2%).

A caveat that belongs in every paper using this: a zero point taken from
the target's own Halpha puts the wavelength scale in the TARGET's rest
frame for that line.  Velocities of Halpha itself are therefore zero by
construction; only the separation to a telluric feature (the O2-B check
here) carries velocity information, and it is limited by the edge-
position error stored beside it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from . import config as gconfig
from . import linecal as lc
from .trace import running_median

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report quotes these).
# --------------------------------------------------------------------------

#: Rest wavelength (Angstrom, air) of the zero-point anchor.
HALPHA_A = lc.HALPHA_A

#: Continuum window for the running-median continuum estimate (px) of the
#: emission-peak SEARCH.  Wide enough that Halpha (FWHM up to ~25 px on a
#: defocused night) does not lift its own continuum; narrow enough to
#: follow the TiO band structure of an M giant.
CONT_WIN = 151

#: Emission peak acceptance: minimum SNR of the continuum-subtracted peak
#: against the robust residual noise, and the width band (px, FWHM-like)
#: a sharp emission line may occupy.  Below 3 px is a cosmic-ray hit;
#: above 80 px is molecular band structure, not a line.
PEAK_MIN_SNR = 8.0
PEAK_MIN_WIDTH = 3
PEAK_MAX_WIDTH = 80

#: Half-window (px) of the Gaussian refinement of the Halpha centre, in
#: units of the search FWHM (floor ``FIT_MIN_HALFWIN``).
FIT_HALFWIN_FWHM = 2.5
FIT_MIN_HALFWIN = 14

#: O2-B check: how far from the predicted position the band edge may be
#: looked for — a constant plus a fraction of the predicted separation —
#: and the minimum depth for the edge to count as measured.  The window
#: is deliberately generous (several per cent): the check must be able to
#: FAIL, i.e. to report a separation that disagrees with the solution.
O2B_SEARCH_PX = 8.0
O2B_SEARCH_FRAC = 0.04
O2B_MIN_DEPTH = 0.03

#: Columns to keep either side of Halpha in the stored spectrum snippet.
SNIPPET_HALFWIN = 150
SNIPPET_STRIDE = 3


# --------------------------------------------------------------------------
# The fixed solution
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Solution:
    """One fixed dispersion solution (a row of ``g_dispersion``).

    The scale is a polynomial in DETECTOR position,
    lambda(x) = const + P(x), P(x) = sum_k coeffs[k-1] ((x - x_ref)/1000)^k
    (see ``linecal.solve_dispersion`` for why it is tied to the detector
    and not to the star).  A frame contributes only the constant, fixed
    by the pixel of Halpha.  ``o2b_wave`` is the measured effective
    wavelength of the O2-B blue edge on this scale (None if not
    calibrated).
    """
    grism: str
    mech_epoch: str
    coeffs: tuple
    x_ref: float
    disp_ref: float              # signed A/px at x_ref
    disp_err: float
    rms_px: Optional[float]
    o2b_wave: Optional[float]
    status: str

    @property
    def red_sign(self) -> int:
        return 1 if self.disp_ref > 0 else -1

    def wave_of(self, x, x_halpha: float):
        return lc.wavelength_of(x, x_halpha, self.coeffs, self.x_ref)

    def x_of(self, wave, x_halpha: float):
        return lc.predict_x(wave, x_halpha, self.coeffs, self.x_ref)

    def local_disp(self, x) -> float:
        """|A/px| at a detector position (it varies across the chip: the
        hrg runs from ~0.40 to ~0.51 A/px over the central 2000 px)."""
        return float(abs(lc.local_dispersion(x, self.coeffs, self.x_ref)))


def load_solutions(con) -> dict:
    """{(grism, mech_epoch): Solution} from ``g_dispersion`` (rows with a
    solved dispersion only)."""
    import json
    out = {}
    for row in con.execute("""
            SELECT grism, mech_epoch, coeffs_json, x_ref, disp_a_per_px,
                   disp_err, rms_px, o2b_wave_eff, status
            FROM g_dispersion WHERE coeffs_json IS NOT NULL"""):
        out[(row[0], row[1])] = Solution(
            grism=row[0], mech_epoch=row[1],
            coeffs=tuple(json.loads(row[2])), x_ref=row[3],
            disp_ref=row[4], disp_err=row[5], rms_px=row[6],
            o2b_wave=row[7], status=row[8])
    return out


def solution_for(solutions: dict, grism: Optional[str],
                 night: str) -> Optional[Solution]:
    """The solution that applies to a frame, by grism unit and the
    mechanical epoch of its night — None when that (grism, epoch) has no
    solution.  There is deliberately NO fallback to another epoch: a
    borrowed dispersion is how a re-seated grism goes unnoticed."""
    if grism is None:
        return None
    return solutions.get((grism, gconfig.mech_epoch_id(night)))


# --------------------------------------------------------------------------
# The zero point
# --------------------------------------------------------------------------
def continuum_residual(flux: np.ndarray, win: int = CONT_WIN):
    """(residual, noise): flux minus its running-median continuum, and the
    robust (MAD) noise of that residual.  NaNs pass through as NaN."""
    finite = np.isfinite(flux)
    filled = np.where(finite, flux, np.nanmedian(flux) if finite.any()
                      else 0.0)
    resid = flux - running_median(filled, win)
    ok = resid[np.isfinite(resid)]
    if len(ok) == 0:
        return resid, 1e-9
    noise = float(1.4826 * np.median(np.abs(ok - np.median(ok))))
    return resid, max(noise, 1e-9)


def emission_candidates(flux: np.ndarray, n: int = 3,
                        min_sep: int = 40) -> list[dict]:
    """Up to ``n`` distinct sharp emission features, strongest first
    (each as :func:`find_emission_peak` returns it).

    Why more than one: on the low-resolution grism the target's ZERO-ORDER
    image lies on the trace line, ~2,900 px from Halpha, and is often the
    single brightest compact feature.  Choosing among candidates is the
    caller's job (by the T CrB fingerprint, see the G runner)."""
    resid, noise = continuum_residual(flux)
    resid = np.where(np.isfinite(resid), resid, 0.0)
    out = []
    work = resid.copy()
    for _ in range(n):
        p = _peak_at(work, noise, int(np.argmax(work)))
        if p is None:
            break
        out.append(p)
        lo = max(0, int(p["x"]) - min_sep)
        work[lo:int(p["x"]) + min_sep + 1] = -np.inf
        if not np.isfinite(work).any():
            break
    return out


def find_emission_peak(flux: np.ndarray) -> Optional[dict]:
    """The strongest sharp emission feature: candidate Halpha (search
    only — :func:`halpha_zero_point` refines it).

    Returns {'x': centroid over the half-maximum span, 'snr',
    'width_px'} or None.
    """
    resid, noise = continuum_residual(flux)
    resid = np.where(np.isfinite(resid), resid, 0.0)
    return _peak_at(resid, noise, int(np.argmax(resid)))


def _peak_at(resid: np.ndarray, noise: float, i: int) -> Optional[dict]:
    """The emission feature whose maximum is at ``i``: None when it is
    weaker than PEAK_MIN_SNR or its half-maximum width is outside the
    allowed band (cosmic ray / molecular structure)."""
    height = resid[i]
    if not np.isfinite(height) or height < PEAK_MIN_SNR * noise:
        return None
    half = height / 2.0
    lo = i
    while lo > 0 and resid[lo - 1] > half:
        lo -= 1
    hi = i
    while hi < len(resid) - 1 and resid[hi + 1] > half:
        hi += 1
    width = hi - lo + 1
    if not (PEAK_MIN_WIDTH <= width <= PEAK_MAX_WIDTH):
        return None
    seg = resid[lo:hi + 1]
    x = float((seg * np.arange(lo, hi + 1)).sum() / seg.sum())
    return {"x": x, "snr": float(height / noise), "width_px": int(width)}


def halpha_zero_point(flux: np.ndarray,
                      var: Optional[np.ndarray] = None,
                      peak: Optional[dict] = None) -> Optional[dict]:
    """Zero point from Halpha EMISSION: the given candidate ``peak`` (or
    the strongest sharp emission feature), refined with a Gaussian +
    linear-continuum fit.

    Returns {'x', 'x_err', 'snr', 'fwhm_px', 'amp_frac'} or None.  When
    the Gaussian fit fails or runs away from the search position, the
    search centroid is returned with ``x_err`` = None (flagged, not
    silently trusted).
    """
    if peak is None:
        peak = find_emission_peak(flux)
    if peak is None:
        return None
    hw = max(FIT_MIN_HALFWIN, int(round(FIT_HALFWIN_FWHM * peak["width_px"])))
    fit = lc.fit_line(flux, peak["x"], hw,
                      core_halfwin=max(3, peak["width_px"]), var=var)
    if (fit is None or fit["amp_frac"] <= 0
            or abs(fit["x"] - peak["x"]) > peak["width_px"]):
        return {"x": peak["x"], "x_err": None, "snr": peak["snr"],
                "fwhm_px": float(peak["width_px"]), "amp_frac": None}
    return {"x": fit["x"], "x_err": fit["x_err"], "snr": peak["snr"],
            "fwhm_px": fit["fwhm_px"], "amp_frac": fit["amp_frac"]}


def check_o2b(flux: np.ndarray, x_halpha: float,
              sol: Solution) -> Optional[dict]:
    """Measure the O2-B blue edge where the fixed solution predicts it
    and report where it actually is, in wavelength.

    Returns {'x', 'x_err', 'depth', 'wave' (edge wavelength on this
    frame's scale), 'wave_err', 'dwave' (wave - calibrated edge
    wavelength), 'sep_px' (pixel distance from Halpha), 'sep_pred_px',
    'frac_dev' ((sep - pred)/pred)} or None when the solution has no
    calibrated edge, the edge falls off the spectrum, or no edge deeper
    than ``O2B_MIN_DEPTH`` is found.  NOTHING here feeds back into the
    wavelength scale.

    Because the zero point is the TARGET's Halpha, ``dwave`` contains the
    target's velocity relative to the observer (1 A = 45.7 km/s) as well
    as any scale error; a series' SCATTER in ``dwave`` is the acceptance
    statistic, its mean is not.
    """
    if sol.o2b_wave is None:
        return None
    xp = float(sol.x_of(sol.o2b_wave, x_halpha))
    sep = abs(xp - x_halpha)
    span = lc.EDGE_FLOOR_SPAN_A["O2B"] / sol.local_disp(xp)
    m = lc.edge_position(flux, xp, sol.red_sign,
                         search_px=O2B_SEARCH_PX + O2B_SEARCH_FRAC * sep,
                         span_px=max(4.0, span))
    if m is None or m["depth"] < O2B_MIN_DEPTH:
        return None
    wave = float(sol.wave_of(m["x"], x_halpha))
    meas = abs(m["x"] - x_halpha)
    return {"x": m["x"], "x_err": m["x_err"], "depth": m["depth"],
            "wave": wave, "wave_err": m["x_err"] * sol.local_disp(m["x"]),
            "dwave": wave - sol.o2b_wave,
            "sep_px": meas, "sep_pred_px": sep,
            "frac_dev": (meas - sep) / sep}


def wavelength_axis(nx: int, x_halpha: float, sol: Solution) -> np.ndarray:
    """The wavelength column of a frame under a fixed solution, given the
    pixel of Halpha.  One definition shared by the FITS writer, the EW
    code and the tests."""
    return sol.wave_of(np.arange(nx), x_halpha)


def snippet(flux: np.ndarray, x_center: float,
            halfwin: int = SNIPPET_HALFWIN,
            stride: int = SNIPPET_STRIDE) -> list:
    """The Halpha-region quick-look: [x, flux] pairs (JSON-ready floats,
    NaN -> None) every ``stride`` px within ``halfwin`` of the line."""
    n = len(flux)
    lo = max(0, int(x_center) - halfwin)
    hi = min(n, int(x_center) + halfwin + 1)
    out = []
    for x in range(lo, hi, stride):
        v = flux[x]
        out.append([int(x), float(v) if np.isfinite(v) else None])
    return out
