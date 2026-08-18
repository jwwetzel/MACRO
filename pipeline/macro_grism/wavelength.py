"""Per-frame self-anchored wavelength solution for the grism track.

There is no arc lamp and no absolute flux in season 1 — the wavelength
scale must come from the spectra themselves.  Anchors, in order of trust:

1. **Halpha emission** (6562.8 A).  T CrB is a symbiotic recurrent nova:
   its Halpha emission is the strongest sharp feature on every validation
   frame (hrg peak SNR >> 10).  Its pixel position anchors the zero point
   of the solution frame by frame — which also makes the anchor immune to
   flexure between frames.
2. **The telluric O2 PAIR** (B band 6867.2 A, A band 7593.7 A).  The
   atmosphere imprints both at fixed wavelengths on every spectrum, so
   their pixel offsets from Halpha must sit in the exact ratio
   (7593.7-6562.8)/(6867.2-6562.8) = 3.387 — and that ratio is the
   identification test.  The validation run PROVED a single-dip search is
   ambiguous: the deepest dip near Halpha on the hrg frames (at +650 px)
   was first read as O2-B (implying 0.47 A/px), but every high-SNR frame
   shows a second dip at +193 px, and 654/193 = 3.39 — the +650 dip is
   O2-A, the +193 dip is O2-B, and the true hrg dispersion is ~1.59 A/px.
   A decoy dip has no partner at the exact ratio on the same side, so the
   PAIR is required before any dispersion is claimed; the A band (3.4x
   the lever arm) supplies the number.  The per-(grism) median across
   pair-anchored frames is the adopted fallback, its robust scatter the
   honest uncertainty.
3. Frames where only Halpha is found get the per-grism FALLBACK dispersion
   with anchor_status = 'halpha_only'; frames with neither anchor get
   'none' and NO wavelength column — an unanchored scale is worse than an
   absent one.

Direction convention: on the era-76 frames red runs toward +x for both
grisms (established by the O2 pair landing redward of Halpha on the
validation frames); the search is symmetric and records the sign it found,
so a future grism mounted the other way round is measured, not assumed.

Precision statement (the report quotes it): the anchors are pixel
positions of extrema; their frame-to-frame scatter (in px, converted
through the dispersion) IS the wavelength precision — nothing finer is
claimed.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from .trace import running_median

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report quotes these).
# --------------------------------------------------------------------------

#: Rest wavelengths (Angstrom, air) of the anchors.
HALPHA_A = 6562.8
O2_B_A = 6867.2          # telluric O2 B band head
O2_A_A = 7593.7          # telluric O2 A band head

#: Continuum window for the running-median continuum estimate (px).  Wide
#: enough that Halpha (FWHM ~20 px on hrg) does not lift its own
#: continuum; narrow enough to follow the TiO band structure.
CONT_WIN = 151

#: Emission peak acceptance: minimum SNR of the continuum-subtracted peak
#: against the robust residual noise, and the width band (px, FWHM-like)
#: a sharp emission line may occupy.  Below 3 px is a cosmic-ray hit;
#: above 80 px is molecular band structure, not a line.
PEAK_MIN_SNR = 8.0
PEAK_MIN_WIDTH = 3
PEAK_MAX_WIDTH = 80

#: Telluric dip acceptance: minimum depth SNR for a dip to enter the
#: pair search at all, and the minimum contiguous pixels below half the
#: acceptance threshold (a one-pixel negative excursion is not a band).
DIP_MIN_SNR = 5.0
DIP_MIN_PIXELS = 2

#: The O2 identification ratio: (A-Halpha)/(B-Halpha) separations.  Both
#: dips of a candidate pair must sit at this ratio from Halpha, within
#: O2_RATIO_TOL (fractional) — the measured pairs on the validation
#: frames hit it to 0.1-2%, and no plausible decoy (TiO band heads,
#: artifacts) lands a partner at the exact ratio on the same side.
O2_RATIO = (O2_A_A - HALPHA_A) / (O2_B_A - HALPHA_A)      # = 3.387
O2_RATIO_TOL = 0.06

#: Physical dispersion bracket (|A/px|) a candidate pair must imply.
#: Measured on the validation frames: hrg ~1.59, lrg ~1.9-2.2; the
#: bracket brackets both with margin without allowing absurd solutions.
DISP_BRACKET_A_PER_PX = (0.3, 3.0)

#: Columns to keep either side of Halpha in the stored spectrum snippet
#: (the g_extractions quick-look; full spectra live in the FITS/parquet).
SNIPPET_HALFWIN = 150
SNIPPET_STRIDE = 3


def continuum_residual(flux: np.ndarray, win: int = CONT_WIN):
    """(residual, noise): flux minus its running-median continuum, and the
    robust (MAD) noise of that residual.  NaNs pass through as NaN."""
    filled = np.where(np.isfinite(flux), flux, np.nanmedian(flux))
    resid = flux - running_median(filled, win)
    finite = resid[np.isfinite(resid)]
    noise = float(1.4826 * np.median(np.abs(finite - np.median(finite))))
    return resid, max(noise, 1e-9)


def find_emission_peak(flux: np.ndarray) -> Optional[dict]:
    """The strongest sharp emission feature: candidate Halpha.

    Returns {'x': subpixel position, 'snr', 'width_px'} or None.  The
    subpixel position is the flux-weighted centroid over the peak's
    half-maximum span — robust for asymmetric lines, and its stability
    across frames is exactly what the anchor-precision figure measures.
    """
    resid, noise = continuum_residual(flux)
    resid = np.where(np.isfinite(resid), resid, 0.0)
    i = int(np.argmax(resid))
    height = resid[i]
    if height < PEAK_MIN_SNR * noise:
        return None
    # Half-maximum span around the peak.
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


def significant_dips(flux: np.ndarray, min_snr: float = DIP_MIN_SNR,
                     min_pixels: int = DIP_MIN_PIXELS,
                     gap: int = 5) -> list[dict]:
    """Every significant absorption dip of a spectrum, as
    [{'x', 'snr', 'width_px'}, ...] sorted deepest-first.

    A dip = a run of pixels (gaps under ``gap`` px bridged) whose
    continuum residual drops below ``-min_snr * noise``, at least
    ``min_pixels`` long (a single-pixel excursion is not a band).  'x' is
    the run's minimum (the band head's core), 'width_px' the span below
    half the dip's own depth — the report's honesty check that a claimed
    O2 band has band-like width.
    """
    resid, noise = continuum_residual(flux)
    resid = np.where(np.isfinite(resid), resid, 0.0)
    idx = np.where(resid < -min_snr * noise)[0]
    if len(idx) == 0:
        return []
    groups = np.split(idx, np.where(np.diff(idx) > gap)[0] + 1)
    dips = []
    for g in groups:
        if len(g) < min_pixels:
            continue
        i = int(g[np.argmin(resid[g])])
        depth = -resid[i]
        lo = i
        while lo > 0 and resid[lo - 1] < -depth / 2:
            lo -= 1
        hi = i
        while hi < len(resid) - 1 and resid[hi + 1] < -depth / 2:
            hi += 1
        dips.append({"x": float(i), "snr": float(depth / noise),
                     "width_px": int(hi - lo + 1)})
    dips.sort(key=lambda d: -d["snr"])
    return dips


def find_o2_pair(flux: np.ndarray, x_ref: float,
                 disp_bracket=DISP_BRACKET_A_PER_PX,
                 ratio_tol: float = O2_RATIO_TOL) -> Optional[dict]:
    """The self-validating O2 anchor: a PAIR of dips whose offsets from
    Halpha sit in the exact A:B separation ratio (``O2_RATIO``), on the
    same side, implying a dispersion inside the physical bracket.

    Among valid pairs the highest combined SNR wins.  Returns
    {'x_o2b', 'o2b_snr', 'x_o2a', 'o2a_snr', 'disp_a_per_px'} — the
    dispersion signed (positive = red toward +x) and taken from the A
    band (3.4x the lever arm of B) — or None when no pair qualifies.
    """
    dips = significant_dips(flux)
    sep_a = O2_A_A - HALPHA_A
    best, best_score = None, 0.0
    for d_b in dips:
        xb = d_b["x"] - x_ref
        if xb == 0:
            continue
        for d_a in dips:
            xa = d_a["x"] - x_ref
            if xa == 0 or xb * xa <= 0:          # same side only
                continue
            ratio = xa / xb                      # sign cancels: both same side
            if not (O2_RATIO * (1 - ratio_tol) <= ratio
                    <= O2_RATIO * (1 + ratio_tol)):
                continue
            disp = sep_a / xa                    # signed A/px
            if not (disp_bracket[0] <= abs(disp) <= disp_bracket[1]):
                continue
            score = d_b["snr"] + d_a["snr"]
            if score > best_score:
                best_score = score
                best = {"x_o2b": d_b["x"], "o2b_snr": d_b["snr"],
                        "x_o2a": d_a["x"], "o2a_snr": d_a["snr"],
                        "disp_a_per_px": float(disp)}
    return best


def solve_wavelength(flux: np.ndarray) -> dict:
    """Per-frame anchor hunt.  Returns a dict of everything found:

    ``x_halpha`` / ``halpha_snr`` / ``halpha_width_px`` — the emission
    anchor; ``x_o2b`` / ``x_o2a`` with SNRs (only ever set as a PAIR —
    see ``find_o2_pair``); ``disp_a_per_px`` (signed: positive = red
    toward +x, measured from the O2-A offset); ``disp_source`` in
    {'o2_pair', None}; ``anchor_status`` in {'halpha+o2', 'halpha_only',
    'none'}.
    """
    out = {"x_halpha": None, "halpha_snr": None, "halpha_width_px": None,
           "x_o2b": None, "o2b_snr": None, "x_o2a": None, "o2a_snr": None,
           "disp_a_per_px": None, "disp_source": None,
           "anchor_status": "none"}
    peak = find_emission_peak(flux)
    if peak is None:
        return out
    out.update(x_halpha=peak["x"], halpha_snr=peak["snr"],
               halpha_width_px=peak["width_px"])
    out["anchor_status"] = "halpha_only"
    pair = find_o2_pair(flux, peak["x"])
    if pair is not None:
        out.update(pair)
        out["disp_source"] = "o2_pair"
        out["anchor_status"] = "halpha+o2"
    return out


def wavelength_axis(nx: int, x_halpha: float,
                    disp_a_per_px: float) -> np.ndarray:
    """The wavelength column: lambda(x) = 6562.8 + disp * (x - x_Halpha).
    Pure arithmetic, kept as a function so the FITS writer and the tests
    share one definition."""
    return HALPHA_A + disp_a_per_px * (np.arange(nx) - x_halpha)


def snippet(flux: np.ndarray, x_center: float,
            halfwin: int = SNIPPET_HALFWIN,
            stride: int = SNIPPET_STRIDE) -> list:
    """The Halpha-region quick-look stored in g_extractions: [x, flux]
    pairs (JSON-ready floats, NaN -> None) every ``stride`` px within
    ``halfwin`` of the line.  Small enough for a DB row, detailed enough
    to eyeball a profile without opening the FITS."""
    n = len(flux)
    lo = max(0, int(x_center) - halfwin)
    hi = min(n, int(x_center) + halfwin + 1)
    out = []
    for x in range(lo, hi, stride):
        v = flux[x]
        out.append([int(x), float(v) if np.isfinite(v) else None])
    return out
