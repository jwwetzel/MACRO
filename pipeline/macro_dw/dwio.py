"""macro_dw.dwio — frame I/O for the Dwarf-Galaxy Hα paper.

Everything here touches a file or pixels; nothing here decides anything a
referee would argue with (the rules are in :mod:`macro_dw.dwcore`).  The
archive is opened READ-ONLY and nothing is written next to it.

Calibration recipe (AC4040, era 2, mechanical epoch AC4040:2023-02-06):

1. **Dark.**  The era's master dark of the frame's readout mode and
   exposure time (``master-dark-{high,spro}-<t>``), which carries the
   pedestal, amp glow and hot pixels.  The masters were taken in April
   2023; the residual pedestal mismatch this leaves (a few to ~20 ADU) is
   measured per frame by :func:`macro_dw.dwcore.fit_sky_on_template` and
   removed before the frame enters a flat, and is otherwise absorbed by
   the per-frame sky.
2. **Flat.**  ``products/dwarf/flats/flat_<f>.fits`` from the ``flats``
   stage: for L, the full-resolution night-sky superflat of the June 2023
   L frames (DW-P11 as ruled); for every other band, the in-window
   2023-06-30 twilight master (pixel-scale structure) times the night-sky
   large-scale correction from that band's own frames.
3. **Astrometry.**  Header (PinPoint) solution when present, else a
   translation vote with the field's reference CD matrix, else a blind
   triangle match — always RE-FITTED to every matched REFCAT2 star with
   SIP-2 (``macro_sn.snio.refine_wcs``) and accepted by
   :func:`macro_dw.dwcore.wcs_ok`.
4. **Photometry.**  sep circular apertures at fixed ANGULAR radii
   (:data:`APER_RADII_ARCSEC`) at the WCS position of every REFCAT2 star,
   local sky from a 20–30" annulus; the growth curve between the radii is
   what turns the 6" measurement aperture into a total flux (DW-P37: the
   declared mode is aperture + curve of growth, no PSF fitting).
"""
from __future__ import annotations

import warnings
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np

warnings.filterwarnings("ignore", module="astropy")

from astropy.io import fits                                  # noqa: E402
from astropy.wcs import WCS, FITSFixedWarning                # noqa: E402

from macro_sn import snio                                    # noqa: E402
from macro_dw import dwcore as core                          # noqa: E402

warnings.simplefilter("ignore", FITSFixedWarning)

ARCHIVE = snio.ARCHIVE
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "products" / "dwarf"
FLATDIR = OUT / "flats"
MASKDIR = OUT / "masks"
BLOCKDIR = OUT / "blocks"
STACKDIR = OUT / "stacks"

#: Twilight masters of 2023-06-30 and the dark-flat of matching exposure
#: (staged pairs in ``stage_dwarfgalaxy_agn_survey``, role master_flat).
TWILIGHT = {
    "L": ("master-flat-L-1x1.fts.fz", "master-dark-flat-01-1x1.fts.fz"),
    "G": ("master-flat-G-1x1.fts.fz", "master-dark-flat-01-1x1.fts.fz"),
    "H": ("master-flat-H-1x1.fts.fz", "master-dark-flat-07-1x1.fts.fz"),
    "R": ("master-flat-R-1x1.fts.fz", "master-dark-flat-013-1x1.fts.fz"),
    "O": ("master-flat-O-1x1.fts.fz", "master-dark-flat-23-1x1.fts.fz"),
    "X": ("master-flat-X-1x1.fts.fz", "master-dark-flat-04-1x1.fts.fz"),
    "W": ("master-flat-W-1x1.fts.fz", "master-dark-flat-04-1x1.fts.fz"),
}
CALIB = "rawimage/calib/"
#: Flat pixels below this normalised response are masked (vignetted
#: corners, dead columns) — the SN pipeline's value, same camera.
FLAT_DEAD = 0.2

#: Photometric aperture radii (arcsec) for the growth curve; sky annulus.
APER_RADII_ARCSEC = (3.0, 4.5, 6.0, 9.0, 13.5)
MEAS_APER_ARCSEC = 6.0


def rkey(r: float) -> str:
    """Column suffix of an aperture radius: 4.5" -> '45'."""
    return f"{int(round(r * 10))}"
SKY_ANNULUS_ARCSEC = (20.0, 30.0)
#: Source-mask growth (px) around the 2-sigma segmentation footprint.
MASK_GROW_PX = 9


def load(rel: str) -> tuple[np.ndarray, fits.Header]:
    return snio.load_data(rel)


def dark_name(readoutm: str, exptime: float) -> str:
    """The era master dark for a readout mode and exposure (nearest of the
    staged ladder 0.125 … 1024 s)."""
    ladder = np.array([0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256,
                       512, 1024], float)
    t = ladder[np.argmin(np.abs(np.log(ladder / max(exptime, 1e-3))))]
    ts = f"{t:g}"
    mode = "spro" if "StackPro" in (readoutm or "") else "high"
    return f"{CALIB}master-dark-{mode}-{ts}-1x1.fts.fz"


@lru_cache(maxsize=8)
def dark(readoutm: str, exptime: float) -> np.ndarray:
    d, _ = load(dark_name(readoutm, exptime))
    return d


@lru_cache(maxsize=8)
def twilight(filt: str) -> tuple[np.ndarray, np.ndarray]:
    """Normalised twilight master (dark-flat subtracted) and its dead mask."""
    fl, df = TWILIGHT[filt]
    a, _ = load(CALIB + fl)
    b, _ = load(CALIB + df)
    t = a - b
    t /= np.median(t[1024:3072, 1024:3072])
    dead = t < FLAT_DEAD
    t[dead] = 1.0
    return t.astype(np.float32), dead


@lru_cache(maxsize=8)
def final_flat(filt: str) -> tuple[np.ndarray, np.ndarray]:
    with fits.open(FLATDIR / f"flat_{filt}.fits") as h:
        f = np.asarray(h[0].data, np.float32)
    dead = ~np.isfinite(f) | (f < FLAT_DEAD)
    f = np.where(dead, 1.0, f).astype(np.float32)
    return f, dead


def dark_subtract(raw: np.ndarray, readoutm: str, exptime: float) -> np.ndarray:
    return raw - dark(readoutm, float(exptime))


def calibrate(raw: np.ndarray, filt: str, readoutm: str, exptime: float
              ) -> tuple[np.ndarray, np.ndarray]:
    """(raw - dark) / final flat; returns (image, bad-pixel mask)."""
    img = dark_subtract(raw, readoutm, exptime)
    f, dead = final_flat(filt)
    img = img / f
    img[dead] = 0.0
    return img.astype(np.float32), dead


def source_mask(img: np.ndarray, base: Optional[np.ndarray] = None,
                thresh: float = 2.0, grow: int = MASK_GROW_PX) -> np.ndarray:
    """Pixels belonging to any source (2-sigma segmentation footprint,
    grown by ``grow`` px), OR'd with ``base``."""
    import sep
    from scipy.ndimage import binary_dilation
    im = np.ascontiguousarray(img, np.float32)
    bm = base if base is not None else np.zeros(im.shape, bool)
    bkg = sep.Background(im, mask=bm, bw=128, bh=128)
    sub = im - bkg.back()
    sep.set_extract_pixstack(20_000_000)
    sep.set_sub_object_limit(8192)
    try:
        _, seg = sep.extract(sub, thresh, err=bkg.globalrms, minarea=8,
                             mask=bm, segmentation_map=True, deblend_cont=1.0)
    except Exception:
        _, seg = sep.extract(sub, 2 * thresh, err=bkg.globalrms, minarea=8,
                             mask=bm, segmentation_map=True, deblend_cont=1.0)
    m = seg > 0
    st = np.ones((3, 3), bool)
    m = binary_dilation(m, structure=st, iterations=grow)
    return m | bm


def block_median(img: np.ndarray, mask: np.ndarray,
                 b: int = core.FLAT_BLOCK, max_masked: float = 0.5
                 ) -> np.ndarray:
    """Masked median in b x b blocks (NaN where > ``max_masked`` masked)."""
    ny, nx = img.shape
    a = np.where(mask, np.nan, img).reshape(ny // b, b, nx // b, b)
    a = a.transpose(0, 2, 1, 3).reshape(ny // b, nx // b, b * b)
    frac = np.isnan(a).mean(axis=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(a, axis=2)
    med[frac > max_masked] = np.nan
    return med


def upsample_blocks(m: np.ndarray, shape=(4096, 4096)) -> np.ndarray:
    """Cubic upsampling of a block map (values at block centres) to pixels."""
    from scipy.ndimage import map_coordinates
    ny, nx = shape
    by, bx = m.shape
    sy, sx = ny / by, nx / bx
    yy = (np.arange(ny) + 0.5) / sy - 0.5
    xx = (np.arange(nx) + 0.5) / sx - 0.5
    # Separable evaluation keeps memory at one full frame.
    from scipy.interpolate import RectBivariateSpline
    spl = RectBivariateSpline(np.arange(by), np.arange(bx), m, kx=3, ky=3)
    del map_coordinates
    return spl(yy, xx).astype(np.float32)


def fill_nan_blocks(m: np.ndarray) -> np.ndarray:
    """Replace NaN blocks by the mean of finite neighbours (iterated)."""
    m = m.copy()
    for _ in range(50):
        bad = ~np.isfinite(m)
        if not bad.any():
            break
        p = np.pad(m, 1, constant_values=np.nan)
        nb = np.stack([p[:-2, 1:-1], p[2:, 1:-1], p[1:-1, :-2], p[1:-1, 2:]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            fill = np.nanmean(nb, axis=0)
        m[bad] = fill[bad]
    return m


# ---------------------------------------------------------------------------
# Astrometry + photometry of one calibrated frame
# ---------------------------------------------------------------------------
def solve(img: np.ndarray, bad: np.ndarray, hdr, cat: dict,
          ref_cd: Optional[np.ndarray]) -> dict:
    """Plate solution refined on REFCAT2.  Returns a dict with ``wcs``
    (or None), ``method``, ``n_match``, ``rms_arcsec``, ``obj`` (sep list)."""
    from scipy.ndimage import median_filter
    sub, bkg, _ = snio.detect(img, bad)
    # Detect on a 3x3-median-filtered copy: isolated hot and RTS pixels (a
    # late-June Hα frame at -14.7 C carries > 10^4 of them) vanish, stars
    # keep their centroids.  ``sub``/``bkg`` stay those of the true frame.
    _, _, obj = snio.detect(median_filter(img, size=3), bad)
    n_raw = len(obj)
    # Hot pixels, RTS pixels and cosmic rays are single-pixel sharp; a
    # warm 512 s Hα frame yields >10^4 of them, and a nearest-neighbour
    # matcher then pairs catalogue stars with noise.  Keep resolved
    # detections only (a, b >= 0.8 px, >= 9 pixels), brightest 2500.
    keep = (obj["a"] >= 0.8) & (obj["b"] >= 0.8) & (obj["npix"] >= 9)
    obj = obj[keep]
    obj = obj[np.argsort(-obj["flux"])[:2500]]
    out = {"wcs": None, "method": None, "n_match": 0, "n_raw": n_raw,
           "rms_arcsec": None, "obj": obj, "sub": sub, "bkg": bkg}
    ra, de = cat["ra"], cat["dec"]
    pt = snio.pointing(hdr)

    def candidates():
        """Starting solutions, cheapest first; generated lazily so the
        blind matcher runs only when everything cheaper has failed."""
        w0 = snio.header_wcs(hdr)
        if w0 is not None:
            yield "header", w0
        if pt is None:
            return
        if ref_cd is not None:
            for name, cd in (("vote", ref_cd), ("vote_flip", -np.asarray(ref_cd))):
                wv = snio.offset_vote_solve(obj, ra, de, cd, pt[0], pt[1],
                                            img.shape)
                if wv is not None:
                    yield name, wv
        wb = snio.blind_solve(obj, ra, de, cat["r"], pt[0], pt[1],
                              core.PIX_ARCSEC_NOMINAL, img.shape)
        if wb is not None:
            yield "blind", wb

    best = None
    for name, w in candidates():
        try:
            wr, n, rms = snio.refine_wcs(w, obj, ra, de, img.shape)
        except Exception:
            continue
        if wr is None:
            continue
        sc = float(np.sqrt(abs(np.linalg.det(wr.pixel_scale_matrix))) * 3600)
        if core.wcs_ok(n, rms, sc):
            best = (name, wr, n, rms, sc)
            break
    if best is not None:
        out.update(method=best[0], wcs=best[1], n_match=best[2],
                   rms_arcsec=best[3], scale=best[4])
    return out


def photometry(sub: np.ndarray, bad: np.ndarray, wcs: WCS, cat: dict,
               gain: float, scale: float, edge: int = 60) -> dict:
    """Forced multi-aperture photometry at every REFCAT2 position on chip.

    ``sub`` is the calibrated frame with sep's mesh background removed —
    only for STAR photometry, where a local annulus sky is subtracted on
    top; the stacks never see a mesh background."""
    import sep
    x, y = wcs.all_world2pix(cat["ra"], cat["dec"], 0)
    ny, nx = sub.shape
    on = (x > edge) & (x < nx - edge) & (y > edge) & (y < ny - edge)
    idx = np.nonzero(on)[0]
    im = np.ascontiguousarray(sub, np.float32)
    rin, rout = (r / scale for r in SKY_ANNULUS_ARCSEC)
    res = {"idx": idx, "x": x[idx], "y": y[idx]}
    for r in APER_RADII_ARCSEC:
        f, e, fl = sep.sum_circle(im, x[idx], y[idx], r / scale, mask=bad,
                                  gain=gain, bkgann=(rin, rout),
                                  err=core.robust_sigma(im[~bad][::97]))
        k = rkey(r)
        res[f"f{k}"], res[f"e{k}"], res[f"flag{k}"] = f, e, fl
    return res


def make_grid(ra0: float, dec0: float, half: int, scale: float) -> WCS:
    return snio.sky_grid(ra0, dec0, half=half, scale=scale)


def resample(img: np.ndarray, w_img: WCS, grid: WCS, half: int,
             scale: float, order: int = 1) -> np.ndarray:
    """Flux-conserving resampling of ``img`` onto ``grid`` (NaN outside)."""
    from scipy.ndimage import map_coordinates
    n = 2 * half + 1
    yy, xx = np.mgrid[0:n, 0:n]
    ra, dec = grid.all_pix2world(xx.ravel(), yy.ravel(), 0)
    # Non-adaptive, divergence detection off: grid points beyond the
    # chip, where a SIP polynomial extrapolates wildly, otherwise drop
    # astropy into a per-point Python loop (minutes per frame).  Those
    # points fall outside the image and become NaN below anyway.
    xi, yi = w_img.all_world2pix(ra, dec, 0, tolerance=1e-3, maxiter=30,
                                 adaptive=False, detect_divergence=False,
                                 quiet=True)
    out = map_coordinates(img, [yi, xi], order=order, mode="constant",
                          cval=np.nan).reshape(n, n)
    a_img = abs(np.linalg.det(w_img.pixel_scale_matrix)) * 3600.0 ** 2
    return (out * (scale ** 2 / a_img)).astype(np.float32)


def wcs_from_text(text: str) -> WCS:
    return WCS(fits.Header.fromstring(text, sep="\n"))


def wcs_to_text(w: WCS) -> str:
    return w.to_header(relax=True).tostring(sep="\n")


# ---------------------------------------------------------------------------
# Stack measurement utilities (DW-P34/P35/P36)
# ---------------------------------------------------------------------------
def read_stack(name: str):
    """(image, weight, WCS, header) of products/dwarf/stacks/<name>.fits."""
    with fits.open(STACKDIR / f"{name}.fits") as h:
        return (np.asarray(h[0].data, float), np.asarray(h[1].data, float),
                WCS(h[0].header), h[0].header)


def stack_mask(img: np.ndarray, thresh: float = 3.0, grow: int = 3
               ) -> np.ndarray:
    """Sources + empty (NaN) pixels on a stack."""
    bad = ~np.isfinite(img)
    return source_mask(np.nan_to_num(img), bad, thresh=thresh, grow=grow)


def aper_net(img: np.ndarray, x: float, y: float, r: float,
             ann: tuple[float, float], mask: Optional[np.ndarray] = None
             ) -> float:
    """Aperture sum minus area x median of a (source-masked) annulus."""
    import sep
    im = np.ascontiguousarray(np.nan_to_num(img), np.float32)
    f, _, _ = sep.sum_circle(im, np.atleast_1d(x), np.atleast_1d(y), r,
                             bkgann=ann, mask=mask)
    return float(f[0])


def random_apertures(img: np.ndarray, mask: np.ndarray, r: float,
                     ann: tuple[float, float], n: int, rng,
                     avoid: Optional[tuple[float, float, float]] = None,
                     edge: float = 60.0) -> np.ndarray:
    """Net sums of ``n`` apertures at random source-free positions (the
    aperture itself may not touch a masked pixel; the annulus may)."""
    import sep
    ny, nx = img.shape
    im = np.ascontiguousarray(np.nan_to_num(img), np.float32)
    good = ~mask
    out = []
    tries = 0
    while len(out) < n and tries < 50 * n:
        tries += 1
        x = rng.uniform(edge + ann[1], nx - edge - ann[1])
        y = rng.uniform(edge + ann[1], ny - edge - ann[1])
        if avoid is not None and np.hypot(x - avoid[0], y - avoid[1]) < avoid[2]:
            continue
        x0, x1 = int(x - r) - 1, int(x + r) + 2
        y0, y1 = int(y - r) - 1, int(y + r) + 2
        yy, xx = np.mgrid[y0:y1, x0:x1]
        disc = (xx - x) ** 2 + (yy - y) ** 2 <= r * r
        if not good[y0:y1, x0:x1][disc].all():
            continue
        f, _, _ = sep.sum_circle(im, np.array([x]), np.array([y]), r,
                                 bkgann=ann, mask=mask)
        out.append(float(f[0]))
    return np.array(out)


def psf_match(a: np.ndarray, fa: float, b: np.ndarray, fb: float):
    """Gaussian-convolve the sharper stack to the broader FWHM (px)."""
    from scipy.ndimage import gaussian_filter
    if fa is None or fb is None or abs(fa - fb) < 0.05:
        return a, b, max(fa or 0, fb or 0)
    if fa < fb:
        s = np.sqrt(fb ** 2 - fa ** 2) / 2.3548
        return _nan_gauss(a, s), b, fb
    s = np.sqrt(fa ** 2 - fb ** 2) / 2.3548
    return a, _nan_gauss(b, s), fa


def _nan_gauss(a: np.ndarray, s: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter
    m = np.isfinite(a)
    num = gaussian_filter(np.where(m, a, 0.0), s)
    den = gaussian_filter(m.astype(float), s)
    out = num / np.maximum(den, 1e-6)
    out[~m] = np.nan
    return out


def stack_fwhm(img: np.ndarray, w: WCS, cat: dict, n: int = 40) -> Optional[float]:
    """FWHM (px) of a stack from Gaussian fits at 14 < r < 17 REFCAT2 stars."""
    x, y = w.all_world2pix(cat["ra"], cat["dec"], 0)
    ny, nx = img.shape
    m = (cat["r"] > 14) & (cat["r"] < 17) & (x > 30) & (x < nx - 30) \
        & (y > 30) & (y < ny - 30)
    idx = np.nonzero(m)[0][:n]
    return snio.grid_fwhm(np.nan_to_num(img), x[idx], y[idx], half=12)
