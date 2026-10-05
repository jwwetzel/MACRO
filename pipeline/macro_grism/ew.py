"""Halpha equivalent width of an emission-line target — one estimator,
applied identically to every library version.

Step 7 of the committee brief asks how the T CrB equivalent widths and
their errors MOVED when the grism library was corrected.  That question
only has an answer if the estimator is held fixed and only the inputs
change, so the estimator lives here, takes plain arrays, and knows
nothing about which library produced them.

DEFINITION (provisional — see the caveat)
-----------------------------------------
    EW = sum over the line window of (F / C - 1) * dlambda     [A]

emission positive; C is a straight line through the medians of two
pseudo-continuum bands.  Windows (Angstrom, in the wavelength scale that
is passed in):

    line        6562.8 +/- 25
    blue band   6500 - 6535
    red band    6590 - 6625

The bands stop short of the TiO 6651 head redward and the telluric
O2-gamma band blueward.  These windows are a WORKING convention for
comparing libraries; the published convention (Munari's, strategy A.5)
is task TCRB-A5's decision, and a different choice rescales every EW
here by the same few per cent without changing how they moved.

ERROR MODEL
-----------
Two terms, added in quadrature:
* photon/read noise in the line window, from the extraction variance;
* the continuum placement: the standard error of each band's median
  (1.2533 sigma / sqrt(N), sigma from the band's own scatter — which
  includes real TiO structure and so is conservative), propagated
  through EW's dependence on the continuum level, d(EW)/d(lnC) =
  -(EW + W) with W the window width.

A dispersion-free companion, ``ew_px`` (the same sum with dlambda = 1
pixel and the windows converted to pixels with the local dispersion),
is returned so that a change in EW can be split into "the pixels
changed" and "the Angstrom-per-pixel changed".
"""

from __future__ import annotations

from typing import Optional

import numpy as np

HALPHA_A = 6562.80
LINE_HALFWIDTH_A = 25.0
BLUE_BAND_A = (6500.0, 6535.0)
RED_BAND_A = (6590.0, 6625.0)

#: Minimum finite pixels in each continuum band and in the line window.
MIN_BAND_PX = 4
MIN_LINE_PX = 5


def equivalent_width(wave: np.ndarray, flux: np.ndarray,
                     var: Optional[np.ndarray] = None,
                     line_center: float = HALPHA_A) -> Optional[dict]:
    """EW of the line at ``line_center`` on a spectrum with wavelength
    axis ``wave`` (monotonic, either direction).

    Returns None when a band or the line window has too few finite
    pixels, or the continuum is not positive.  Otherwise {'ew_a',
    'ew_err_a', 'ew_px', 'cont' (continuum level at the line, flux
    units), 'cont_snr' (continuum level / per-pixel scatter in the
    bands), 'n_line', 'n_masked' (non-finite pixels inside the window —
    saturated or rejected columns, which bias an emission EW LOW),
    'peak' (max flux in the window)}.
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    ok = np.isfinite(wave) & np.isfinite(flux)
    shift = line_center - HALPHA_A

    def band(lo, hi):
        return ok & (wave >= lo + shift) & (wave <= hi + shift)

    b, r = band(*BLUE_BAND_A), band(*RED_BAND_A)
    in_win = (np.isfinite(wave)
              & (np.abs(wave - line_center) <= LINE_HALFWIDTH_A))
    line = in_win & ok
    if b.sum() < MIN_BAND_PX or r.sum() < MIN_BAND_PX \
            or line.sum() < MIN_LINE_PX:
        return None
    wb, wr = float(np.mean(wave[b])), float(np.mean(wave[r]))
    cb, cr = float(np.median(flux[b])), float(np.median(flux[r]))
    slope = (cr - cb) / (wr - wb)
    cont = cb + slope * (wave - wb)
    c0 = cb + slope * (line_center - wb)
    if c0 <= 0 or np.any(cont[line] <= 0):
        return None
    dl = np.abs(np.gradient(wave))
    ratio = flux / cont - 1.0
    ew = float(np.sum(ratio[line] * dl[line]))
    ew_px = float(np.sum(ratio[line]))
    width = float(np.sum(dl[line]))
    # Continuum placement error.
    sb = 1.4826 * np.median(np.abs(flux[b] - cb))
    sr = 1.4826 * np.median(np.abs(flux[r] - cr))
    se_c = 0.5 * np.hypot(1.2533 * sb / np.sqrt(b.sum()),
                          1.2533 * sr / np.sqrt(r.sum()))
    err_cont = (abs(ew) + width) * se_c / c0
    # Photon / read noise in the window.
    if var is not None:
        v = np.asarray(var, dtype=float)
        good = line & np.isfinite(v)
        err_phot = float(np.sqrt(np.sum(v[good] * (dl[good] / cont[good]) ** 2)))
    else:
        err_phot = float(0.5 * (sb + sr) / c0 * np.sqrt(np.sum(dl[line] ** 2)))
    return {"ew_a": ew, "ew_err_a": float(np.hypot(err_phot, err_cont)),
            "ew_px": ew_px, "cont": float(c0),
            "cont_snr": float(c0 / max(0.5 * (sb + sr), 1e-9)),
            "n_line": int(line.sum()),
            "n_masked": int(in_win.sum() - line.sum()),
            "peak": float(np.max(flux[line]))}
