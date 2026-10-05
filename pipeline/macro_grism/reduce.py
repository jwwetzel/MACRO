"""One frame, pixels to 1-D spectrum — the function every G script calls.

Before this module each runner re-implemented "load, find the trace,
extract" with its own constants, which is how one script came to use a
gain the other did not.  :func:`reduce_frame` is now the single path:

    load (fits_io) -> detector (config) -> trace (trace) -> extent ->
    extraction (extract) -> header cards + quality numbers

and it returns everything as plain numpy arrays and Python scalars so the
caller decides what to store.  It never writes.

The header cards it carries out are the ones the committee asked every
project to regress against (TE.F4, DE.F5): true focuser position
(pyscope's ``FOCPOS`` — MaxIm's ``FOCUSPOS`` is stuck for whole months,
TE.F8), sensor temperature and set-point, flip status, airmass.  The
manifest does not have them yet (finding F-2, owned by ``foundation-s0``);
reading them here from the pixels' own header keeps this track
self-sufficient until it does.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from . import config as gconfig
from . import extract as gext
from . import trace as gtrace
from .fits_io import load_frame

#: Header cards copied into the result (missing cards become None).
HEADER_CARDS = ("DATE-OBS", "JD", "EXPTIME", "FILTER", "INSTRUME",
                "READOUTM", "XBINNING", "GAIN", "OFFSET", "EGAIN",
                "CCD-TEMP", "SET-TEMP", "FOCPOS", "FOCUSPOS", "FOCUSTEM",
                "FLIPSTAT", "TELPIER", "AIRMASS", "OBJCTRA", "OBJCTDEC",
                "RA", "DEC", "OBJECT", "SWCREATE", "FWPOS")

#: Minimum main-trace height (halo-subtracted median-collapse ADU) for a
#: frame to count as containing a spectrum at all.  Good 240 s T CrB
#: frames measure 40-260; a cloud-dead frame measured 4.
MIN_TRACE_HEIGHT_ADU = 15.0


def header_subset(header) -> dict:
    """The regressor cards as a plain dict (floats where numeric)."""
    out = {}
    for card in HEADER_CARDS:
        v = header.get(card)
        if isinstance(v, (int, float, np.integer, np.floating)):
            v = float(v)
        elif v is not None:
            v = str(v).strip()
        out[card] = v
    return out


#: A frame's own dark corners must sit within this many ADU of the
#: tabulated pedestal, else the pedestal is taken from the frame.
PEDESTAL_TOL_ADU = 20.0
CORNER_PX = 200


def frame_pedestal(det, data: np.ndarray, hdr: dict):
    """The detector record with the pedestal THIS frame actually has.

    The tabulated pedestal belongs to one camera setting; the ASI was
    also run at OFFSET 50 (pedestal ~503 ADU instead of ~303 at OFFSET
    30), and a 200 ADU error in the pedestal is 200 ADU of phantom sky
    whose photon noise the variance model then predicts (found by the G-2
    test, 2026-10-05: measured/predicted 0.2 on those frames).  A slitless
    frame's corners receive no first-order sky, so their median IS the
    pedestal (plus dark, negligible at -10 C).  When it differs from the
    table by more than PEDESTAL_TOL_ADU the frame's value is used and the
    provenance says so; otherwise the measured table value stands."""
    from dataclasses import replace
    c = CORNER_PX
    corners = [data[:c, :c], data[:c, -c:], data[-c:, :c], data[-c:, -c:]]
    ped = float(np.median([np.median(q) for q in corners]))
    if abs(ped - det.pedestal_adu) <= PEDESTAL_TOL_ADU:
        return det
    return replace(det, pedestal_adu=ped,
                   provenance=det.provenance + f"; pedestal {ped:.1f} from "
                   f"frame corners (OFFSET card {hdr.get('OFFSET')})")


def hot_pixel_mask(hdr: dict, shape) -> tuple:
    """(mask, note): the detector package's bad-pixel mask for this
    frame's camera, geometry, orientation and sensor temperature
    (``rlmt_diagnostics.badpix``; hot, rail and noisy pixels), or
    (None, 'none') when no product matches — the extraction then relies on
    the Horne outlier test alone, and the row says so."""
    from rlmt_diagnostics import badpix
    ny, nx = shape
    try:
        m = badpix.mask_for_frame(hdr.get("INSTRUME"), nx, ny,
                                  readoutm=hdr.get("READOUTM"),
                                  ccd_temp=hdr.get("CCD-TEMP"),
                                  flipstat=hdr.get("FLIPSTAT"))
    except Exception as exc:                          # noqa: BLE001
        return None, f"error: {type(exc).__name__}"
    if m is None:
        return None, "none"
    return m, f"badpix ({int(m.sum())} px)"


def find_trace(data: np.ndarray) -> dict:
    """Trace geometry of the brightest spectrum on a frame.

    Returns slope, u (detilted row of the main trace), height (its
    halo-subtracted amplitude), the fitted polynomial ``coeffs``, the
    number of centroid blocks used, the fit rms, and ``extent`` — the
    (first, last) column over which the trace actually carries signal
    (None when no block passes).
    """
    xs, ys, amps = gtrace.chunk_peaks(data)
    slope = gtrace.fit_slope(xs, ys, amps)
    _, resid = gtrace.detilted_profile(data, slope)
    u, height = gtrace.main_trace_u(resid)
    coeffs, n_used, rms = gtrace.fit_trace_centers(data, slope, u)
    extent = gtrace.trace_extent(data, coeffs) if rms is not None else None
    return {"slope": float(slope), "u": float(u), "height": float(height),
            "coeffs": coeffs, "n_centroids": int(n_used), "rms_px": rms,
            "extent": extent}


def reduce_frame(path: str, night: Optional[str] = None,
                 sky: str = "poly", dark: Optional[np.ndarray] = None,
                 hot: Optional[np.ndarray] = None,
                 data_header=None, diagnostics: bool = False) -> dict:
    """Reduce one grism frame to a 1-D spectrum.

    ``night`` (local observing night) is only the tie-breaker for
    detector identification on blank pyscope headers.  ``data_header``
    lets a caller that already holds the pixels pass ``(data, header,
    layout)`` instead of re-reading the file.

    Returns a dict with: ``status`` ('ok' | 'no_trace'), ``layout``,
    ``shape``, ``detector`` (the config record used), ``header`` (card
    subset), ``trace`` (see :func:`find_trace`), and — when status is
    'ok' — ``spec`` (the :func:`extract.extract_spectrum` dict) plus the
    quality scalars ``peak_adu`` (brightest raw pixel in the aperture
    inside the trace extent), ``n_sat_cols``, ``snr_median``,
    ``fwhm_px_median`` (cross-dispersion FWHM, median over chunks) and
    ``sky_median_adu`` (sky above pedestal under the trace).

    With ``diagnostics=True`` the pixels are also used, in the same
    pass (a frame costs seconds to read off the archive disk; nothing
    should have to read it twice), for the three tests that need 2-D
    data: ``diag['flux_flanking']`` / ``['box_flanking']`` (the legacy
    straight-line-sky spectrum — the G-4 method difference),
    ``diag['null_<method>_<offset>']`` (empty-aperture sums — the G-4
    null test) and ``diag['skyvar']`` (measured vs predicted sky
    variance — the G-2 test).
    """
    data, header, layout = (data_header if data_header is not None
                            else load_frame(path))
    hdr = header_subset(header)
    det = gconfig.detector_for(hdr.get("INSTRUME"), hdr.get("READOUTM"),
                               night)
    det = frame_pedestal(det, data, hdr)
    hot_note = "given"
    if hot is None:
        hot, hot_note = hot_pixel_mask(hdr, data.shape)
    tr = find_trace(data)
    out = {"status": "ok", "layout": layout, "shape": data.shape,
           "detector": det, "header": hdr, "trace": tr,
           "hot_mask": hot_note,
           "sat_cap_adu": float(det.saturation_cap_adu())}
    if (tr["height"] < MIN_TRACE_HEIGHT_ADU or tr["rms_px"] is None
            or tr["extent"] is None):
        out["status"] = "no_trace"
        return out
    spec = gext.extract_spectrum(data, tr["coeffs"], det, dark=dark,
                                 sky=sky, hot=hot)
    x0, x1 = tr["extent"]
    inside = np.zeros(data.shape[1], dtype=bool)
    inside[x0:x1 + 1] = True
    spec["inside"] = inside
    f, v = spec["flux"], spec["var"]
    sel = inside & np.isfinite(f) & np.isfinite(v) & (v > 0)
    out["spec"] = spec
    out["peak_adu"] = (float(np.nanmax(spec["peak"][inside]))
                       if inside.any() else None)
    out["n_sat_cols"] = int((spec["n_sat"][inside] > 0).sum())
    out["snr_median"] = (float(np.median(f[sel] / np.sqrt(v[sel])))
                         if sel.any() else None)
    fw = spec["fwhm_px"]
    out["fwhm_px_median"] = (float(np.nanmedian(fw))
                             if len(fw) and np.isfinite(fw).any() else None)
    bg = spec["bg"][sel]
    out["sky_median_adu"] = (float(np.median(bg) - det.pedestal_adu)
                             if sel.any() else None)
    if diagnostics:
        legacy = gext.extract_spectrum(data, tr["coeffs"], det, dark=dark,
                                       sky="flanking", hot=hot)
        diag = {"flux_flanking": legacy["flux"],
                "box_flanking": legacy["box"],
                "skyvar": gext.sky_variance_bins(data, tr["coeffs"], det)}
        for method in ("poly", "flanking"):
            for off in gext.NULL_OFFSETS:
                tag = f"null_{method}_{'m' if off < 0 else 'p'}{abs(off)}"
                diag[tag] = gext.null_aperture(data, tr["coeffs"], det,
                                               method, off)
        out["diag"] = diag
    return out
