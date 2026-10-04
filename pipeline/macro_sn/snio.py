"""macro_sn.snio — the frame I/O layer of the SN 2023ixf photometry.

Everything here touches a file; nothing here decides anything a referee
would argue with.  The decisions (which stars calibrate, how errors add,
what a limit means) live in :mod:`macro_sn.snphot`, which is pure and
unit-tested.  This module turns one archive frame into numbers:

1. **calibrate** — AC4040 High Gain frames: subtract the 0.1 s dark-flat
   master (the bias pedestal, measured 93 ADU by S2, plus a negligible dark)
   and divide by the campaign's one in-window master flat of the same filter
   (2023-06-30, strategy §3.5), normalised to its central median.  Pixels
   whose flat response is under :data:`FLAT_DEAD` are masked.  Frames from
   the other two cameras (the 2024 iKon and 2026 QHY600 template epochs)
   have no in-window master in the archive, so they are only background-
   subtracted; they feed differential (cutout) work, never the zero point.
2. **astrometry** — a frame that carries a PinPoint solution starts from it;
   one that does not is solved blind against REFCAT2 by triangle-invariant
   point matching (astroalign), which is rotation- and scale-free and is run
   in both parities.  Either way the WCS is then RE-FITTED to every matched
   REFCAT2 star (``astropy.wcs.utils.fit_wcs_from_points``) and accepted only
   if :func:`macro_sn.snphot.wcs_acceptable` passes on its matched count and
   sky residual.  The header solution is never trusted unrefined.
3. **forced photometry** — circular apertures of radius
   ``APER_FWHM x FWHM`` at the WCS position of every REFCAT2 star in the
   frame and of the supernova, local sky from a ``4-6 x FWHM`` annulus
   (strategy §4 Step 5), photon + background errors from sep with the
   header gain.  The RAW native-pixel peak inside ``max(FWHM, 3 px)`` of
   each position is recorded beside it, because the saturation screen is
   judged at the target, per frame, in native pixels (standing rule 4).

The archive is opened READ-ONLY; nothing is written next to it.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

warnings.filterwarnings("ignore", module="astropy")

from astropy.io import fits                                  # noqa: E402
from astropy.wcs import WCS, FITSFixedWarning                # noqa: E402
from astropy.coordinates import SkyCoord                     # noqa: E402
import astropy.units as u                                    # noqa: E402

warnings.simplefilter("ignore", FITSFixedWarning)

ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")

#: The campaign camera's calibration masters (all from the one in-window
#: flat epoch, night 2023-06-30; ``calib_frames`` in the manifest).  The
#: 0.1 s dark-flat is used as the bias for every exposure: S2 measures the
#: High Gain pedestal at 93 ADU, and at -18 C the dark current over even a
#: 128 s narrowband exposure is a few ADU — removed anyway by the local sky
#: annulus, which is why no exposure-matched dark is needed for aperture
#: photometry.
BIAS_MASTER = "rawimage/calib/master-dark-flat-01-1x1.fts.fz"
FLAT_MASTERS = {f: f"rawimage/calib/master-flat-{f}-1x1.fts.fz"
                for f in ("G", "R", "I", "H", "O", "1")}
#: Flat-field pixels below this normalised response are masked (vignetted
#: corners and dead columns; 0.17% of the array).
FLAT_DEAD = 0.2

#: Aperture radius and sky annulus in units of the frame's FWHM — the
#: strategy's §4 Step 5 values.
APER_FWHM = 1.5
ANNULUS_FWHM = (4.0, 6.0)

#: Edge margin (px) inside which no forced aperture is placed.
EDGE_PX = 40

#: Detection threshold for the astrometric source list (background sigma).
DETECT_SIGMA = 5.0


@dataclass
class Calib:
    """The AC4040 masters, loaded once per worker process."""
    bias: np.ndarray
    flats: dict
    dead: dict


def load_data(rel_path: str) -> tuple[np.ndarray, fits.Header]:
    """First HDU with data, as float32, plus its header (read-only)."""
    with fits.open(ARCHIVE / rel_path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is not None:
                return np.asarray(hdu.data, dtype=np.float32), hdu.header
    raise ValueError(f"no image data in {rel_path}")


def load_calib() -> Calib:
    """Read the bias and the six normalised flats."""
    bias, _ = load_data(BIAS_MASTER)
    flats, dead = {}, {}
    for f, p in FLAT_MASTERS.items():
        fl, _ = load_data(p)
        fl = fl - bias
        fl /= np.median(fl[1024:3072, 1024:3072])
        d = fl < FLAT_DEAD
        fl[d] = 1.0
        flats[f], dead[f] = fl, d
    return Calib(bias=bias, flats=flats, dead=dead)


def calibrate(raw: np.ndarray, filt: str, camera: str,
              calib: Optional[Calib]) -> tuple[np.ndarray, np.ndarray, str]:
    """Return (calibrated image, bad-pixel mask, recipe name)."""
    if camera == "AC4040" and calib is not None and raw.shape == calib.bias.shape:
        img = raw - calib.bias
        if filt in calib.flats:
            img = img / calib.flats[filt]
            mask = calib.dead[filt].copy()
            recipe = f"bias(0.1s dark-flat)+flat({filt},2023-06-30)"
        else:
            mask = np.zeros(raw.shape, bool)
            recipe = "bias(0.1s dark-flat), no flat for this code"
        img[mask] = 0.0
        return img, mask, recipe
    return raw.copy(), np.zeros(raw.shape, bool), "none (no in-window masters)"


def header_wcs(h: fits.Header) -> Optional[WCS]:
    """The header's celestial WCS if it has one, else None."""
    if not str(h.get("CTYPE1", "")).startswith("RA"):
        return None
    try:
        w = WCS(h).celestial
        return w if w.has_celestial else None
    except Exception:
        return None


def pointing(h: fits.Header) -> Optional[tuple[float, float]]:
    """Commanded pointing (deg) from OBJCTRA/OBJCTDEC or RA/DEC cards."""
    for kr, kd in (("OBJCTRA", "OBJCTDEC"), ("RA", "DEC")):
        r, d = h.get(kr), h.get(kd)
        if r is None or d is None:
            continue
        try:
            if isinstance(r, str) and (":" in r or " " in r.strip()):
                c = SkyCoord(r, d, unit=(u.hourangle, u.deg))
                return c.ra.deg, c.dec.deg
            return float(r), float(d)
        except Exception:
            continue
    return None


def trial_wcs(ra0: float, dec0: float, scale_arcsec: float,
              nx: int, ny: int, parity: int) -> WCS:
    """A north-up TAN WCS at the pointing — the frame the catalogue is
    projected into for blind matching.  ``parity`` = +1 or -1 flips x."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    s = scale_arcsec / 3600.0
    w.wcs.cd = np.array([[-s * parity, 0.0], [0.0, s]])
    return w


def detect(img: np.ndarray, mask: np.ndarray):
    """Background model + 5-sigma source list (sep)."""
    import sep
    img = np.ascontiguousarray(img, dtype=np.float32)
    bkg = sep.Background(img, mask=mask, bw=64, bh=64)
    sub = img - bkg.back()
    sep.set_extract_pixstack(5_000_000)
    sep.set_sub_object_limit(4096)
    try:
        obj = sep.extract(sub, DETECT_SIGMA, err=bkg.globalrms, minarea=5,
                          mask=mask)
    except Exception:
        # A saturated 256 s frame of the galaxy can overflow the deblender;
        # a higher threshold still finds every star the matcher needs.
        obj = sep.extract(sub, 4 * DETECT_SIGMA, err=bkg.globalrms,
                          minarea=5, mask=mask, deblend_nthresh=16)
    return sub, bkg, obj


def match_xy(ax, ay, bx, by, tol):
    """Nearest-neighbour match of points a to points b within ``tol`` px.
    Returns index arrays (ia, ib) of one-to-one pairs."""
    from scipy.spatial import cKDTree
    if len(bx) == 0 or len(ax) == 0:
        return np.array([], int), np.array([], int)
    tree = cKDTree(np.c_[bx, by])
    d, j = tree.query(np.c_[ax, ay], distance_upper_bound=tol)
    ok = np.isfinite(d)
    ia, ib = np.nonzero(ok)[0], j[ok]
    # Enforce one-to-one: keep the closest a for every b.
    order = np.argsort(d[ok])
    seen, keep = set(), []
    for k in order:
        if ib[k] not in seen:
            seen.add(ib[k])
            keep.append(k)
    keep = np.array(sorted(keep), int)
    return ia[keep], ib[keep]


def refine_wcs(w: WCS, obj, cat_ra, cat_dec, shape, sip: bool = True):
    """Re-fit a WCS to all REFCAT2 stars matched under it.

    Two passes: a 15 px match to measure the bulk offset (header solutions
    are good to a few px, blind ones to ~1 px), then a 2.5 px match under
    the shifted solution, then ``fit_wcs_from_points``.  Returns
    ``(wcs, n_matched, rms_arcsec)`` or ``(None, n, None)``.
    """
    from astropy.wcs.utils import fit_wcs_from_points
    ny, nx = shape
    x, y = w.all_world2pix(cat_ra, cat_dec, 0)
    inside = (x > 5) & (x < nx - 5) & (y > 5) & (y < ny - 5)
    if inside.sum() < 8 or len(obj) < 8:
        return None, int(inside.sum()), None
    ox, oy = obj["x"], obj["y"]
    ci = np.nonzero(inside)[0]
    w_cur = w
    for tol in (15.0, 4.0, 2.5):
        x, y = w_cur.all_world2pix(cat_ra[ci], cat_dec[ci], 0)
        ia, ib = match_xy(x, y, ox, oy, tol)
        if len(ia) < 8:
            return None, len(ia), None
        sky = SkyCoord(cat_ra[ci][ia], cat_dec[ci][ia], unit="deg")
        deg = 2 if (sip and len(ia) >= 60) else None
        w_cur = fit_wcs_from_points((ox[ib], oy[ib]), sky,
                                    projection="TAN", sip_degree=deg)
    x, y = w_cur.all_world2pix(cat_ra[ci], cat_dec[ci], 0)
    ia, ib = match_xy(x, y, ox, oy, 2.5)
    if len(ia) < 8:
        return None, len(ia), None
    pred = w_cur.all_pix2world(ox[ib], oy[ib], 0)
    dra = (pred[0] - cat_ra[ci][ia]) * np.cos(np.radians(cat_dec[ci][ia]))
    rms = float(np.sqrt(np.mean((dra ** 2 + (pred[1] - cat_dec[ci][ia]) ** 2)))
                * 3600.0)
    return w_cur, len(ia), rms


def blind_solve(obj, cat_ra, cat_dec, cat_mag, ra0, dec0, scale, shape):
    """Triangle-invariant match of the 60 brightest detections to the 120
    brightest catalogue stars near the pointing, in both parities.  Returns
    a starting WCS (to be refined) or None."""
    import astroalign as aa
    from astropy.wcs.utils import fit_wcs_from_points
    ny, nx = shape
    order = np.argsort(-obj["flux"])[:60]
    src = np.c_[obj["x"][order], obj["y"][order]]
    best = None
    for parity in (+1, -1):
        tw = trial_wcs(ra0, dec0, scale, nx, ny, parity)
        x, y = tw.all_world2pix(cat_ra, cat_dec, 0)
        # Generous window: the commanded pointing can miss by ~250 px.
        pad = 0.25 * max(nx, ny)
        ok = (x > -pad) & (x < nx + pad) & (y > -pad) & (y < ny + pad)
        idx = np.nonzero(ok)[0]
        idx = idx[np.argsort(cat_mag[idx])][:120]
        if len(idx) < 10 or len(src) < 10:
            continue
        dst = np.c_[x[idx], y[idx]]
        try:
            T, (s_pts, d_pts) = aa.find_transform(src, dst,
                                                  max_control_points=60)
        except Exception:
            continue
        if len(s_pts) < 8:
            continue
        # Matched catalogue points back to the sky through the trial WCS.
        ra, dec = tw.all_pix2world(d_pts[:, 0], d_pts[:, 1], 0)
        w0 = fit_wcs_from_points((s_pts[:, 0], s_pts[:, 1]),
                                 SkyCoord(ra, dec, unit="deg"),
                                 projection="TAN")
        if best is None or len(s_pts) > best[1]:
            best = (w0, len(s_pts))
    return None if best is None else best[0]


def fwhm_of(obj, peak_cap: float) -> Optional[float]:
    """Median FWHM (px) of well-measured, unsaturated compact detections."""
    a, b = obj["a"], obj["b"]
    good = (obj["peak"] < peak_cap) & (obj["flux"] > 0) & (a > 0.5) \
        & (b / np.maximum(a, 1e-3) > 0.6) & (obj["flag"] == 0)
    if good.sum() < 5:
        return None
    sel = np.argsort(-obj["flux"][good])[:80]
    f = 2.0 * np.sqrt(np.log(2.0) * (a[good][sel] ** 2 + b[good][sel] ** 2))
    return float(np.median(f))


def native_peaks(raw: np.ndarray, x, y, radius: float) -> np.ndarray:
    """Maximum RAW ADU within ``radius`` px of each position (native
    pixels, before any calibration — the quantity the screen judges)."""
    r = int(np.ceil(radius))
    out = np.full(len(x), np.nan)
    ny, nx = raw.shape
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    disk = xx ** 2 + yy ** 2 <= radius ** 2
    for k, (xi, yi) in enumerate(zip(x, y)):
        cx, cy = int(round(xi)), int(round(yi))
        if cx - r < 0 or cy - r < 0 or cx + r >= nx or cy + r >= ny:
            continue
        out[k] = raw[cy - r:cy + r + 1, cx - r:cx + r + 1][disk].max()
    return out


def forced_phot(sub, bkg, mask, x, y, fwhm, gain):
    """Forced aperture photometry; returns (flux, err, sep_flag)."""
    import sep
    r = APER_FWHM * fwhm
    flux, err, flag = sep.sum_circle(
        np.ascontiguousarray(sub, dtype=np.float32), x, y, r,
        err=bkg.rms(), gain=gain, mask=mask,
        bkgann=(ANNULUS_FWHM[0] * fwhm, ANNULUS_FWHM[1] * fwhm))
    return flux, err, flag


def badpix_mask(task: dict, shape) -> Optional[np.ndarray]:
    """The detector package's bad-pixel mask for this frame, or None when
    no product matches its camera/geometry/orientation (then the frame is
    measured unmasked and the recipe says so)."""
    try:
        from rlmt_diagnostics import badpix
        mk = badpix.mask_for_frame(task.get("instrume") or task.get("camera"),
                                   int(shape[1]), int(shape[0]),
                                   readoutm=task.get("readoutm"),
                                   ccd_temp=task.get("ccd_temp"),
                                   flipstat=task.get("flipstat"))
    except Exception:
        return None
    if mk is None or mk.shape != tuple(shape):
        return None
    return mk


# ===========================================================================
# Template subtraction on a common sky grid (the faint regime, the overlap
# test and the late-time stacks).  Deliberately simple and fully stated:
#   1. every image is resampled (bilinear, flux-conserving by the pixel-area
#      ratio) onto ONE north-up TAN grid centred on the supernova at the
#      campaign scale (0.54"/px), using its own re-fitted WCS;
#   2. each is put on a common flux scale, ZP 25 in its own natural system,
#      from its ensemble zero point (sn_cal_frames.zp_nat);
#   3. stacks are pixel medians of the scaled, resampled frames;
#   4. the sharper image is convolved with a Gaussian to the broader one's
#      FWHM (both FWHMs measured on the grid from REFCAT2 stars);
#   5. difference = science - template; a planar local background is fitted
#      in a 40-80 px annulus and removed;
#   6. forced aperture photometry at the supernova (radius 1.5 x matched
#      FWHM, no further sky), and the same aperture at random positions on
#      the arm defines the empirical noise.
# ===========================================================================
GRID_SCALE_ARCSEC = 0.54
GRID_HALF = 200


def sky_grid(ra0: float, dec0: float, half: int = GRID_HALF,
             scale: float = GRID_SCALE_ARCSEC) -> WCS:
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [half + 1.0, half + 1.0]
    s = scale / 3600.0
    w.wcs.cd = np.array([[-s, 0.0], [0.0, s]])
    return w


def resample_to_grid(img: np.ndarray, w_img: WCS, grid: WCS,
                     half: int = GRID_HALF) -> np.ndarray:
    """Bilinear resampling of ``img`` onto ``grid`` (NaN outside), flux-
    conserving: values are multiplied by (grid pixel area / image pixel
    area)."""
    from scipy.ndimage import map_coordinates
    n = 2 * half + 1
    yy, xx = np.mgrid[0:n, 0:n]
    ra, dec = grid.all_pix2world(xx.ravel(), yy.ravel(), 0)
    xi, yi = w_img.all_world2pix(ra, dec, 0)
    out = map_coordinates(img, [yi, xi], order=1, mode="constant",
                          cval=np.nan).reshape(n, n)
    a_img = abs(np.linalg.det(w_img.pixel_scale_matrix)) * 3600.0 ** 2
    a_grid = (GRID_SCALE_ARCSEC) ** 2
    return out * (a_grid / a_img)


def grid_fwhm(img: np.ndarray, xs, ys, half: int = 10) -> Optional[float]:
    """Median FWHM (px) of circular-Gaussian + constant fits to stamps at
    KNOWN star positions (catalogue, not detections: on a defocused frame a
    detection list is dominated by hot pixels and cosmic rays, and raw
    second moments by noise in the wings).  Needs >= 3 good fits."""
    from scipy.optimize import least_squares
    vals = []
    ny, nx = img.shape
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    for x, y in zip(xs, ys):
        cx, cy = int(round(x)), int(round(y))
        if cx < half + 2 or cy < half + 2 or cx > nx - half - 3 \
                or cy > ny - half - 3:
            continue
        st = img[cy - half:cy + half + 1, cx - half:cx + half + 1]
        if not np.all(np.isfinite(st)):
            continue
        b0 = float(np.median(np.r_[st[0], st[-1], st[:, 0], st[:, -1]]))
        a0 = float(st.max() - b0)
        if a0 <= 0:
            continue

        def resid(p):
            a, dx, dy, sg, b = p
            return (a * np.exp(-0.5 * ((xx - dx) ** 2 + (yy - dy) ** 2)
                               / sg ** 2) + b - st).ravel()

        try:
            r = least_squares(resid, [a0, x - cx, y - cy, 2.0, b0],
                              bounds=([0, -4, -4, 0.3, -np.inf],
                                      [np.inf, 4, 4, half / 2.0, np.inf]))
        except Exception:
            continue
        a, dx, dy, sg, b = r.x
        noise = np.std(r.fun)
        if r.success and a > 10 * noise and 0.3 < sg < half / 2.0 - 0.01:
            vals.append(2.3548 * sg)
    return float(np.median(vals)) if len(vals) >= 3 else None


def gauss_match(a: np.ndarray, fa: float, b: np.ndarray, fb: float):
    """Convolve the sharper of (a, b) to the other's FWHM.  Returns
    (a', b', matched FWHM)."""
    from scipy.ndimage import gaussian_filter
    if fa is None or fb is None:
        return a, b, max(fa or 0, fb or 0) or 4.0
    if fa < fb:
        sig = np.sqrt(max(fb ** 2 - fa ** 2, 0)) / 2.3548
        return np.nan_to_num(gaussian_filter(np.nan_to_num(a), sig)), b, fb
    sig = np.sqrt(max(fa ** 2 - fb ** 2, 0)) / 2.3548
    return a, np.nan_to_num(gaussian_filter(np.nan_to_num(b), sig)), fa


def plane_background(img: np.ndarray, r_in: float = 40, r_out: float = 80):
    """Least-squares plane fitted to the annulus about the grid centre
    (3-sigma clipped); returns the image with it removed."""
    n = img.shape[0]
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]
    rr = np.hypot(xx - c, yy - c)
    m = (rr > r_in) & (rr < r_out) & np.isfinite(img)
    A = np.c_[np.ones(m.sum()), xx[m] - c, yy[m] - c]
    v = img[m]
    keep = np.ones(len(v), bool)
    for _ in range(3):
        coef, *_ = np.linalg.lstsq(A[keep], v[keep], rcond=None)
        res = v - A @ coef
        s = 1.4826 * np.median(np.abs(res[keep]))
        keep = np.abs(res) < 3 * max(s, 1e-9)
    return img - (coef[0] + coef[1] * (xx - c) + coef[2] * (yy - c))


def aper_sum(img: np.ndarray, x: float, y: float, r: float) -> float:
    """Plain circular-aperture sum (sub-pixel by 5x5 supersampling)."""
    n = img.shape[0]
    R = int(np.ceil(r)) + 1
    cx, cy = int(round(x)), int(round(y))
    if cx - R < 0 or cy - R < 0 or cx + R >= n or cy + R >= n:
        return np.nan
    st = img[cy - R:cy + R + 1, cx - R:cx + R + 1]
    sub = (np.arange(5) + 0.5) / 5 - 0.5
    yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
    frac = np.zeros_like(st, dtype=float)
    for dy in sub:
        for dx in sub:
            frac += ((xx + dx - (x - cx)) ** 2 + (yy + dy - (y - cy)) ** 2
                     <= r * r)
    frac /= 25.0
    return float(np.nansum(st * frac))


def gaussian_psf(n: int, x: float, y: float, fwhm: float, flux: float):
    """A unit-normalised circular Gaussian times ``flux`` on an n x n grid."""
    yy, xx = np.mgrid[0:n, 0:n]
    s = fwhm / 2.3548
    g = np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s ** 2)
    return flux * g / (2 * np.pi * s * s)


def offset_vote_solve(obj, cat_ra, cat_dec, ref_cd, ra0, dec0, shape,
                      search_px: float = 700.0, bin_px: float = 3.0):
    """Translation-only solve for a camera whose rotation and scale are
    KNOWN (the campaign's AC4040 sat at one orientation all season, TE).

    The catalogue is projected through the reference CD matrix at the
    commanded pointing; every (detection - star) offset within
    ``search_px`` votes in a 2-D histogram, and the peak bin is the
    pointing error.  Robust with as few as ~6 matched stars, which is what
    the short or thin frames that defeated PinPoint offer.  Returns a
    starting WCS (to be refined) or None when no bin stands out.
    """
    ny, nx = shape
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [nx / 2.0, ny / 2.0]
    w.wcs.cd = np.asarray(ref_cd, float)
    x, y = w.all_world2pix(cat_ra, cat_dec, 0)
    pad = search_px
    ok = (x > -pad) & (x < nx + pad) & (y > -pad) & (y < ny + pad)
    x, y = x[ok], y[ok]
    order = np.argsort(-obj["flux"])[:150]
    ox, oy = obj["x"][order], obj["y"][order]
    if len(ox) < 5 or len(x) < 5:
        return None
    dx = (ox[:, None] - x[None, :]).ravel()
    dy = (oy[:, None] - y[None, :]).ravel()
    m = (np.abs(dx) < search_px) & (np.abs(dy) < search_px)
    nb = int(2 * search_px / bin_px)
    H, xe, ye = np.histogram2d(dx[m], dy[m], bins=nb,
                               range=[[-search_px, search_px]] * 2)
    k = np.unravel_index(np.argmax(H), H.shape)
    peak = H[k]
    bg = np.median(H) + 1.0
    if peak < 6 or peak < 5 * bg:
        return None
    sx = 0.5 * (xe[k[0]] + xe[k[0] + 1])
    sy = 0.5 * (ye[k[1]] + ye[k[1] + 1])
    w.wcs.crpix = [nx / 2.0 + sx, ny / 2.0 + sy]
    return w


def shift_refine(w: WCS, obj, cat_ra, cat_dec, shape, tol: float = 3.0):
    """Translation-only refinement under a FIXED CD matrix: median offset
    of detections matched within ``tol`` px, iterated twice.  Used only
    for the thin frames where a 6-parameter refit has too few stars.
    Returns (wcs, n_matched, rms_arcsec) or (None, n, None)."""
    ny, nx = shape
    w = w.deepcopy()
    n, rms = 0, None
    for _ in range(3):
        x, y = w.all_world2pix(cat_ra, cat_dec, 0)
        ins = (x > 5) & (x < nx - 5) & (y > 5) & (y < ny - 5)
        ia, ib = match_xy(x[ins], y[ins], obj["x"], obj["y"], tol)
        n = len(ia)
        if n < 4:
            return None, n, None
        dx = np.median(obj["x"][ib] - x[ins][ia])
        dy = np.median(obj["y"][ib] - y[ins][ia])
        w.wcs.crpix = [w.wcs.crpix[0] + dx, w.wcs.crpix[1] + dy]
        rx = obj["x"][ib] - x[ins][ia] - dx
        ry = obj["y"][ib] - y[ins][ia] - dy
        scale = np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix))) * 3600.0
        rms = float(np.sqrt(np.mean(rx ** 2 + ry ** 2)) * scale)
    return w, n, rms
