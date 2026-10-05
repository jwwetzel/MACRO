"""Trace geometry for slitless grism frames.

The Mode0 grism frames (era 76, 4788x3194 bin2) show one dominant
first-order trace: a near-horizontal streak (slope ~ +0.03 px/px for hrg,
~ -0.1 for lrg) riding on a broad scattered-light halo, with a forest of
faint instrumental ghosts parallel to it.  The exploration round proved the
ghosts are *instrumental* — a T CrB frame and a frame of a field 168 deg
away show secondary peaks at the SAME offsets from their main trace — so
nothing here may treat secondary peaks as field stars without independent
evidence (that lesson is written into the gate design).

All decision arithmetic is pure numpy on arrays the tests can synthesize;
the only inputs are pixel arrays and numbers.
"""

from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report quotes these).
# --------------------------------------------------------------------------

#: Column chunks used for the coarse trace-slope fit.  48 chunks of ~100
#: columns each: wide enough that a column median kills cosmic rays, many
#: enough to constrain a straight line well.
SLOPE_CHUNKS = 48

#: Window of the running median filter that models the scattered-light halo
#: in a cross-dispersion profile (px).  Much wider than a trace (FWHM ~
#: 5-10 px) and much narrower than the halo (~1000 px), so it removes the
#: halo and keeps the traces.
HALO_WIN = 301

#: Same idea for the collapsed (detilted) profile, where traces are
#: sharpened by the alignment — a tighter window tracks halo curvature.
PROFILE_HALO_WIN = 151

#: Chunks whose peak amplitude falls below this fraction of the median
#: chunk amplitude are dropped from the slope fit (cloud gaps, frame edges).
MIN_CHUNK_AMP_FRAC = 0.3

#: The detilted profile is built from the central span of columns
#: (fractions of the width) — the trace is guaranteed present there and the
#: frame edges (vignetting, trace run-off) stay out of the median.
COLLAPSE_X_LO = 0.125
COLLAPSE_X_HI = 0.875

#: Column stride of the collapse.  Every 9th column of the central 75% of a
#: 4788-px frame is ~400 samples — the median is stable and the loop cheap.
COLLAPSE_STRIDE = 9

#: Half-width of the window used for per-column flux-weighted trace
#: centroids, and the degree of the polynomial fitted through them.
CENTROID_HALFWIN = 12
TRACE_POLY_DEG = 3

#: Column block of the centroid fit (median-combined) and the search
#: half-windows of the first, wide passes (px).  The low-resolution trace
#: departs from its chord by up to ~35 px, so the first pass must look
#: +/-40 px; the passes then tighten onto the curve.
TRACE_BLOCK = 8
TRACE_PASS_HALFWINS = (40, 20)

#: Largest fraction of centroid blocks one robust-fit pass may clip.
TRACE_MAX_CLIP = 0.2

#: A column whose background-subtracted window sum falls below this many
#: times the local noise is excluded from the centroid fit (no signal — a
#: centroid of noise is a random number).
CENTROID_MIN_SNR = 5.0


def running_median(arr: np.ndarray, win: int) -> np.ndarray:
    """Running median with reflected edges — the halo model.

    scipy's ``median_filter`` in 'reflect' mode, wrapped here so every
    module shares one definition of "the smooth background of a profile".
    """
    from scipy.ndimage import median_filter
    return median_filter(arr, size=win, mode="reflect")


def chunk_peaks(data: np.ndarray, n_chunks: int = SLOPE_CHUNKS,
                halo_win: int = HALO_WIN):
    """Coarse trace samples: for each column chunk, the row of the
    strongest halo-subtracted peak of the chunk's median profile.

    Returns (x_centers, y_peaks, amplitudes) as float arrays.  The column
    median inside a chunk suppresses cosmic rays and single hot columns;
    the halo subtraction stops the broad scattered-light bump from
    outvoting a real trace.
    """
    ny, nx = data.shape
    xs, ys, amps = [], [], []
    for i in range(n_chunks):
        x0, x1 = i * nx // n_chunks, (i + 1) * nx // n_chunks
        prof = np.median(data[:, x0:x1], axis=1)
        resid = prof - running_median(prof, halo_win)
        iy = int(np.argmax(resid))
        xs.append((x0 + x1) / 2.0)
        ys.append(float(iy))
        amps.append(float(resid[iy]))
    return np.array(xs), np.array(ys), np.array(amps)


def fit_slope(xs: np.ndarray, ys: np.ndarray, amps: np.ndarray,
              min_amp_frac: float = MIN_CHUNK_AMP_FRAC) -> float:
    """Straight-line slope through the strong chunk peaks.

    Chunks below ``min_amp_frac`` of the median amplitude are dropped:
    where the trace is weak (clouds, run-off past the frame edge) the
    argmax lands on noise and would drag the fit.  A line (not the deg-2
    curve used later for extraction) is enough here: the slope's only job
    is the detilt shift, where a 10% slope error moves a trace end by
    ~±3 px — well inside the peak-matching tolerance.
    """
    good = amps > min_amp_frac * np.median(amps)
    if good.sum() < 3:
        # Too few trustworthy chunks to fit anything: report a flat trace
        # rather than a line through noise (the gate's height floor will
        # reject such a frame downstream anyway).
        return 0.0
    # Robust line: a single bright compact source off the trace (the
    # low-resolution grism's zero-order image sits on the same detector
    # rows' chunks at the far end) must not tilt it — clip chunks more
    # than 3 MAD-sigma (floor 5 px) from the line and refit.
    keep = good.copy()
    for _ in range(5):
        c = np.polyfit(xs[keep], ys[keep], 1)
        r = ys - np.polyval(c, xs)
        sig = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        new = good & (np.abs(r) <= max(3.0 * sig, 5.0))
        if new.sum() < 3 or (new == keep).all():
            break
        keep = new
    return float(c[0])


def detilted_profile(data: np.ndarray, slope: float,
                     x_lo_frac: float = COLLAPSE_X_LO,
                     x_hi_frac: float = COLLAPSE_X_HI,
                     stride: int = COLLAPSE_STRIDE):
    """Collapse the frame along the dispersion into one cross-dispersion
    profile, after removing the trace tilt.

    Each sampled column is rolled by ``-slope * (x - nx/2)`` so every
    trace becomes horizontal, then the column stack is median-combined
    (median: a star's zero-order blob or a satellite in a few columns
    cannot print through).  Returns ``(profile, halo_subtracted)``.

    The detilt coordinate ("u") this defines — u = y - slope*(x - nx/2) —
    is shared with the gate's prediction side: because every star's trace
    is the SAME translated curve, a star at pixel (x*, y*) produces a
    peak at u = y* - slope*(x* - nx/2), the grism's along-dispersion
    deflection cancelling exactly.
    """
    ny, nx = data.shape
    cols = np.arange(int(nx * x_lo_frac), int(nx * x_hi_frac), stride)
    stack = np.empty((len(cols), ny), dtype=np.float32)
    for k, x in enumerate(cols):
        stack[k] = np.roll(data[:, x], -int(round(slope * (x - nx / 2))))
    prof = np.median(stack, axis=0)
    return prof, prof - running_median(prof, PROFILE_HALO_WIN)


def main_trace_u(halo_subtracted: np.ndarray) -> tuple[int, float]:
    """(u position, height) of the strongest peak of the detilted profile
    — the frame's dominant trace, presumed the pointed target."""
    u = int(np.argmax(halo_subtracted))
    return u, float(halo_subtracted[u])


def fit_trace_centers(data: np.ndarray, slope: float, u_main: float,
                      halfwin: int = CENTROID_HALFWIN,
                      deg: int = TRACE_POLY_DEG):
    """Refine the main trace: block centroids inside a window that follows
    the current trace model, then a polynomial through the good
    centroids — ITERATED, so the window follows the curve it is fitting.

    Returns ``(poly_coeffs, n_used, rms_px)``; ``np.polyval(coeffs, x)``
    is the trace center at column x.

    Why it iterates (2026-10-03 revision).  The first version centroided
    once, inside +/-``halfwin`` px of the coarse straight line.  A real
    trace is curved: on the low-resolution grism it departs from the
    chord by 30+ px toward the ends, so the centroid window slid off the
    trace exactly where the curvature was largest, the polynomial was
    fitted to sky, and the "extraction aperture" then sat beside the
    spectrum (visible as trace peaks 10-35 px off the aperture center in
    the rectified cutouts).  Now pass 1 uses a wide window
    (``TRACE_SEARCH_HALFWIN``) around the straight line and each later
    pass re-centroids in the narrow window around the previous
    polynomial, with sigma clipping; the loop stops when the fit moves
    by less than 0.05 px rms.

    Blocks (``TRACE_BLOCK`` columns, median-combined so one cosmic ray
    cannot move a centroid) with window S/N below ``CENTROID_MIN_SNR``
    are excluded — a centroid of noise is a random number.  Falls back
    to the coarse line when fewer than 10 blocks qualify (``rms_px`` is
    then None, and the caller must treat the trace as unverified).
    """
    ny, nx = data.shape
    x_blk = np.arange(TRACE_BLOCK // 2, nx, TRACE_BLOCK)
    line = np.zeros(deg + 1)
    line[-1] = u_main - slope * (nx / 2)
    line[-2] = slope
    coeffs, n_used, rms = line, 0, None
    passes = TRACE_PASS_HALFWINS + (halfwin,) * 3
    for it, hw in enumerate(passes):
        yc = np.polyval(coeffs, x_blk)
        cx, cy = _block_centroids(data, x_blk, yc, hw)
        if len(cx) < 10:
            break
        # The degree RAMPS UP pass by pass (1, 2, then ``deg``): a cubic
        # fitted to the first, wide-window centroids can bend through
        # the blocks of one half of the trace and run away over the
        # other, and the next window then follows the runaway (seen on
        # the 2025-04-20 T CrB lrg frame, 2026-10-05).
        d_it = min(deg, it + 1)
        # Robust polynomial: fit, clip at 3 sigma (MAD), refit — but
        # never discard more than TRACE_MAX_CLIP of the blocks: clipping
        # that wants to remove more is fitting the wrong curve, not
        # removing outliers.
        keep = np.ones(len(cx), dtype=bool)
        for _ in range(4):
            c = np.polyfit(cx[keep], cy[keep], d_it)
            r = cy - np.polyval(c, cx)
            sig = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
            new = np.abs(r) <= 3.0 * max(sig, 0.05)
            if (new.sum() < max(10, (1 - TRACE_MAX_CLIP) * len(cx))
                    or (new == keep).all()):
                break
            keep = new
        c = np.concatenate([np.zeros(deg - d_it), c])   # common length
        moved = np.sqrt(np.mean((np.polyval(c, x_blk)
                                 - np.polyval(coeffs, x_blk)) ** 2))
        coeffs, n_used = c, int(keep.sum())
        rms = float(np.sqrt(np.mean(r[keep] ** 2)))
        if it >= len(TRACE_PASS_HALFWINS) and moved < 0.05:
            break
    return coeffs, n_used, rms


def _block_centroids(data: np.ndarray, x_blk: np.ndarray, yc: np.ndarray,
                     halfwin: int):
    """Flux-weighted cross-dispersion centroids of column blocks.

    For each block center x (``TRACE_BLOCK`` columns, median-combined),
    the window ``yc +/- halfwin`` is baseline-subtracted (median of the
    window's outer quarter on each side — a local linear-free pedestal
    that the trace core cannot bias) and centroided over its positive
    part.  Returns (x, y) arrays of the blocks that pass the S/N cut.
    """
    ny, nx = data.shape
    cx, cy = [], []
    half_blk = TRACE_BLOCK // 2
    for x, y in zip(x_blk, yc):
        lo, hi = int(round(y)) - halfwin, int(round(y)) + halfwin + 1
        if lo < 0 or hi > ny:
            continue
        seg = data[lo:hi, max(0, x - half_blk):x + half_blk + 1]
        win = np.median(seg, axis=1).astype(np.float64)
        q = max(2, len(win) // 4)
        edge = np.concatenate([win[:q], win[-q:]])
        base = np.median(edge)
        noise = 1.4826 * np.median(np.abs(edge - base)) + 1e-3
        sig = win - base
        if sig.max() < CENTROID_MIN_SNR * noise:
            continue                             # no believable signal
        w = np.clip(sig - 2.0 * noise, 0, None)   # core only: wings and
        if w.sum() <= 0:                         # sky slope cannot pull
            continue
        cx.append(float(x))
        cy.append(lo + float((w * np.arange(len(win))).sum() / w.sum()))
    return np.array(cx), np.array(cy)


def trace_extent(data: np.ndarray, coeffs: np.ndarray,
                 halfwin: int = CENTROID_HALFWIN,
                 min_snr: float = CENTROID_MIN_SNR):
    """(x_first, x_last): the column range over which the trace actually
    carries signal.

    A slitless trace does not span the detector: the low-resolution grism
    puts the whole 4000-9500 A spectrum on ~2700 of 4788 columns, and
    past its ends the polynomial is an EXTRAPOLATION through empty sky.
    Everything downstream (line search, sky-method statistics) must stay
    inside this range.  Measured from the same block centroids as the
    fit: the first and last block whose window passes the S/N cut and
    whose centroid lies within 2 px of the polynomial.
    """
    ny, nx = data.shape
    x_blk = np.arange(TRACE_BLOCK // 2, nx, TRACE_BLOCK)
    cx, cy = _block_centroids(data, x_blk, np.polyval(coeffs, x_blk),
                              halfwin)
    if len(cx) == 0:
        return None
    on = np.abs(cy - np.polyval(coeffs, cx)) <= 2.0
    if not on.any():
        return None
    return int(cx[on].min()), int(cx[on].max())


def profile_noise(halo_subtracted: np.ndarray, exclude_u: int,
                  exclude_halfwin: int = 200) -> float:
    """Robust noise of a detilted profile, measured away from the main
    trace (± ``exclude_halfwin`` px around it is masked so the trace and
    its ghost forest cannot inflate their own detection threshold)."""
    mask = np.ones(len(halo_subtracted), dtype=bool)
    lo = max(0, exclude_u - exclude_halfwin)
    hi = min(len(halo_subtracted), exclude_u + exclude_halfwin)
    mask[lo:hi] = False
    r = halo_subtracted[mask]
    return float(1.4826 * np.median(np.abs(r - np.median(r))))
