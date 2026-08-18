"""Thin I/O layer: read one reduced frame, detect sources, measure fluxes.

This is the ONLY module in ``macro_phot`` that touches pixels.  All policy
(thresholds, aperture sizes, clip levels) lives in the pure module
``macro_phot.photometry``; this file just wires sep to those constants and
returns plain dictionaries the build script writes to the database.

Provenance rule (decision 2026-08-18): the S4 prototype reads SERVER-REDUCED
pixels only — the ``reduced/`` tree reached through ``raw_reduced_links`` —
and never mixes reduction provenances within one (target, era) series.
"""

from __future__ import annotations

import itertools
from contextlib import contextmanager
from pathlib import Path
from typing import Optional

import numpy as np
import sep
from astropy.io import fits

from . import photometry as ph
# Sibling package under pipeline/ — the single source of truth for turning
# a FITS header into image dimensions (tile compression included).
from macro_core import fitsgeom

#: Keep only the brightest this-many detections for triangle matching
#: (astroalign cost grows fast with control points; the bright end carries
#: all the geometry).
MAX_ALIGN_STARS = 100

#: The astroalign attempt ladder: (bright-pool size, control points),
#: tried in order until one converges.  Two empirical lessons bought with
#: real frames: (1) sparse high-latitude fields (AN UMa, ~200-400
#: detections) fit fastest with astroalign's default 50 control points;
#: (2) the dense Milky-Way field of VV Pup (b ~ +2 deg, 1,000-5,000
#: detections) needs a DEEPER pool — the brightest-100 sets of a sharp
#: reference and a defocused frame overlap too little there, and only a
#: 300-star pool with 100 control points recovers the transform; (3) the
#: Andor iKon December nights only lock onto their reference with an
#: 800-star pool (dense field + night-to-night depth changes).  A failed
#: attempt costs 0.3-2 s, so the ladder is cheap for the frames that pass
#: rung 1 and decisive for the ones that need the deep rungs.
ALIGN_ATTEMPTS: tuple[tuple[int, int], ...] = ((100, 50), (300, 100),
                                               (800, 200))

#: sep sometimes needs a bigger pixel-buffer for crowded/big frames.
SEP_PIXSTACK = 1_000_000


def read_reduced(path: Path) -> tuple[np.ndarray, dict]:
    """Read one reduced frame: float32 pixels + the header cards S4 needs.

    Reduced files are fpack tiles (CompImageHDU in extension 1) or plain
    primary HDUs; both int16-with-BZERO and float32 payloads occur.  The
    pixel array is cast to native float32 C-order (sep's requirement).
    """
    with fits.open(path) as hdul:
        hdu = hdul[1] if len(hdul) > 1 and hdul[1].data is not None else hdul[0]
        data = np.ascontiguousarray(hdu.data, dtype=np.float32)
        hdr = hdu.header
        # TRUE geometry via the compression-aware resolver, resolved ONCE.
        # A tile-compressed header that astropy could not translate reports
        # NAXIS1 as the BINTABLE row length in bytes (8), not the image
        # width — see macro_core.fitsgeom and the S0e report.
        naxis1, naxis2 = fitsgeom.resolve_geometry_or_none(hdr)
        meta = {
            "exptime": hdr.get("EXPTIME"),
            "jd": hdr.get("JD"),
            "filter": hdr.get("FILTER"),
            "xpixsz": hdr.get("XPIXSZ"),
            "focallen": hdr.get("FOCALLEN"),
            "egain": hdr.get("EGAIN"),
            "pltsolvd": 1 if hdr.get("PLTSOLVD") else 0,
            "naxis1": naxis1,
            "naxis2": naxis2,
            "airmass": hdr.get("AIRMASS"),
            # The frame's own linear sky transform, when it carries one.
            # This is the FALLBACK plate scale for headers whose XPIXSZ /
            # FOCALLEN cards are unusable: the March-2026 'Fast' frames
            # write FOCALLEN = 0.0 while carrying a perfectly good CD
            # matrix, and reading the scale off that matrix rescues a whole
            # observing block that otherwise fails extraction outright.
            # See ``macro_phot.series.resolve_plate_scale``.
            "cd1_1": hdr.get("CD1_1"), "cd1_2": hdr.get("CD1_2"),
            "cd2_1": hdr.get("CD2_1"), "cd2_2": hdr.get("CD2_2"),
            "cdelt1": hdr.get("CDELT1"), "cdelt2": hdr.get("CDELT2"),
        }
    return data, meta


def measure_frame(data: np.ndarray, egain: Optional[float],
                  aper_px: float) -> tuple[dict, dict]:
    """Detect sources and measure aperture fluxes on one reduced frame.

    Parameters
    ----------
    data
        Background-included float32 image (as read).
    egain
        Header e-/ADU gain, used for the shot-noise term of the flux error;
        ``None`` or 0 (the Andor iKon writes EGAIN=0) drops the source shot
        term — the predicted errors are then sky-limited, and the empirical
        inflation factor (section 5 of the report) absorbs the difference.
    aper_px
        Aperture RADIUS in pixels (from the pure layer's sky-fixed scaling).

    Returns
    -------
    (frame_stats, detections)
        ``frame_stats``: n_detected, background level/rms, median FWHM px.
        ``detections``: parallel arrays (x, y, flux, fluxerr, fwhm, peak,
        flag, clipped) for every detection, aperture-measured with a local
        sky annulus.
    """
    sep.set_extract_pixstack(SEP_PIXSTACK)
    bkg = sep.Background(data)
    sub = data - bkg.back()
    objs = sep.extract(sub, ph.DETECT_SIGMA, err=bkg.globalrms,
                       minarea=ph.DETECT_MINAREA)
    x, y = objs["x"], objs["y"]
    # Aperture sum with LOCAL sky annulus (bkgann re-estimates the sky per
    # star inside the annulus, on top of the global model subtraction —
    # belt and braces against background-model residuals near bright stars).
    ann = (ph.SKY_ANNULUS_ARCSEC[0] / ph.APERTURE_RADIUS_ARCSEC * aper_px,
           ph.SKY_ANNULUS_ARCSEC[1] / ph.APERTURE_RADIUS_ARCSEC * aper_px)
    gain = float(egain) if egain else None   # 0/None -> no shot term
    flux, fluxerr, flag = sep.sum_circle(
        sub, x, y, aper_px, err=bkg.globalrms, gain=gain, bkgann=ann)
    fwhm = ph.fwhm_from_ab(objs["a"], objs["b"])
    clipped = objs["peak"] + bkg.globalback > ph.PEAK_CLIP_ADU
    stats = {
        "n_detected": int(len(objs)),
        "bkg_adu": float(bkg.globalback),
        "bkg_rms": float(bkg.globalrms),
        "fwhm_px": float(np.median(fwhm)) if len(objs) else None,
    }
    dets = {
        "x": x, "y": y, "flux": flux, "fluxerr": fluxerr, "fwhm": fwhm,
        "peak": objs["peak"], "flag": flag.astype(int),
        "clipped": clipped.astype(int),
    }
    return stats, dets


def brightest_xy(dets: dict, n_max: int = MAX_ALIGN_STARS) -> np.ndarray:
    """(N, 2) positions of the brightest unclipped detections, for astroalign."""
    order = np.argsort(-np.asarray(dets["flux"]))
    keep = [i for i in order if not dets["clipped"][i]][:n_max]
    return np.column_stack([np.asarray(dets["x"])[keep],
                            np.asarray(dets["y"])[keep]])


@contextmanager
def seeded_ransac(seed: Optional[int]):
    """Make astroalign's RANSAC deterministic inside this context.

    astroalign 2.6.2 shuffles its RANSAC trial order with
    ``np.random.default_rng().shuffle(...)`` — a FRESH generator seeded
    from OS entropy on every call, so neither ``np.random.seed`` nor any
    global generator state can reach it, and two identical calls return
    (slightly) different transforms.  That broke the repo's regenerable-
    products discipline: re-running the match stage flipped a handful of
    star identities per frame.

    The fix routes zero-argument ``np.random.default_rng()`` calls made
    INSIDE the context to generators seeded from ``(seed, call_number)``
    (a valid numpy SeedSequence entropy tuple), so every RANSAC shuffle in
    one attempt ladder draws a deterministic, distinct stream.  Calls WITH
    an explicit seed argument pass through untouched, and the original
    function is always restored on exit.  ``seed=None`` disables the
    patch (library behavior unchanged).  Single-threaded use only — which
    is how the build script runs.
    """
    if seed is None:
        yield
        return
    orig = np.random.default_rng
    calls = itertools.count()
    def _routed(arg=None):              # noqa: ANN001 - mirrors numpy API
        if arg is not None:
            return orig(arg)
        return orig((int(seed) & 0xFFFFFFFF, next(calls)))
    np.random.default_rng = _routed
    try:
        yield
    finally:
        np.random.default_rng = orig


#: The PRODUCTION ladder used by the CV campaign
#: (``pipeline/scripts/run_cv_photometry.py``).  Same three pools as
#: :data:`ALIGN_ATTEMPTS`, but the deepest rung fits with 100 control points
#: instead of 200.  The reason is throughput measured on real frames: for
#: the dense unsolved fields (YZ Cnc, EU UMa) rungs 1 and 2 fail in about
#: two seconds together and rung 3 decides the frame — and rung 3 at 200
#: control points costs ~25 s, which at six workers is thirteen frames a
#: minute and puts a 2,500-frame block three hours away.  astroalign's
#: triangle-matching cost grows steeply with the control-point count while
#: its REACH is set by the pool depth, which is unchanged here; the
#: measured Gaia-tie probe locked the same field at 100 control points.
#: Every accepted transform is still validated the same way afterwards (the
#: one-to-one match rate and the alignment RMS are recorded per frame), so
#: this is a speed-for-reach trade made where reach is not what is scarce.
PRODUCTION_ALIGN_ATTEMPTS: tuple[tuple[int, int], ...] = ((100, 50),
                                                          (300, 100),
                                                          (800, 100))


def find_series_transform(frame_bright: np.ndarray, ref_bright: np.ndarray,
                          seed: Optional[int] = None,
                          attempts: Optional[tuple] = None):
    """Fit frame -> reference with the :data:`ALIGN_ATTEMPTS` ladder.

    Both inputs are flux-descending (brightest first) position arrays at
    least as deep as the ladder's largest pool (shorter arrays simply use
    what they have).  Returns the astroalign transform of the first rung
    that converges; raises the LAST rung's exception when none does, so
    the caller records a real astroalign error, not a wrapper's.

    ``seed`` (the build script passes the frame_id) pins astroalign's
    RANSAC via :func:`seeded_ransac`, making the returned transform — and
    therefore every star identity downstream — reproducible run to run.
    ``attempts`` overrides the ladder itself; the CV campaign passes
    :data:`PRODUCTION_ALIGN_ATTEMPTS`, whose cheaper deepest rung keeps the
    same reach at a fraction of the cost (see that constant's note).
    """
    import astroalign as aa
    last_err: Exception | None = None
    ladder = attempts if attempts is not None else ALIGN_ATTEMPTS
    with seeded_ransac(seed):
        for pool, ctrl in ladder:
            try:
                tf, _ = aa.find_transform(frame_bright[:pool],
                                          ref_bright[:pool],
                                          max_control_points=ctrl)
                return tf
            except Exception as e:      # MaxIterError, ValueError (too few)
                last_err = e
    raise last_err


