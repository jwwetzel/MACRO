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
TRACE_POLY_DEG = 2

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
    return float(np.polyfit(xs[good], ys[good], 1)[0])


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
    """Refine the main trace: per-column flux-weighted centroid inside a
    window that follows the coarse (slope, u) line, then a deg-2
    polynomial through the good centroids.

    Returns ``(poly_coeffs, n_used, rms_px)``; ``np.polyval(coeffs, x)``
    is the trace center at column x.  Columns with window SNR below
    ``CENTROID_MIN_SNR`` are excluded (an emission-line star's continuum
    can vanish into the noise between features; fitting centroids there
    would chase noise).  Falls back to the coarse line when fewer than 10
    columns qualify.
    """
    ny, nx = data.shape
    xs = np.arange(0, nx, 4)                    # every 4th column: plenty
    cx, cy = [], []
    for x in xs:
        yc = u_main + slope * (x - nx / 2)      # coarse center here
        lo, hi = int(yc) - halfwin, int(yc) + halfwin + 1
        if lo < 0 or hi > ny:
            continue
        win = data[lo:hi, x].astype(np.float64)
        base = np.median(win)                    # local pedestal
        sig = win - base
        noise = 1.4826 * np.median(np.abs(sig - np.median(sig))) + 1e-3
        if sig.sum() < CENTROID_MIN_SNR * noise * np.sqrt(len(win)):
            continue                             # no believable signal
        w = np.clip(sig, 0, None)
        if w.sum() <= 0:
            continue
        cx.append(x)
        cy.append(lo + float((w * np.arange(len(win))).sum() / w.sum()))
    if len(cx) < 10:
        # Coarse fallback: the straight line as a degree-``deg`` poly.
        coeffs = np.zeros(deg + 1)
        coeffs[-1] = u_main - slope * (nx / 2)
        coeffs[-2] = slope
        return coeffs, len(cx), None
    cx, cy = np.array(cx, float), np.array(cy)
    coeffs = np.polyfit(cx, cy, deg)
    rms = float(np.sqrt(np.mean((np.polyval(coeffs, cx) - cy) ** 2)))
    return coeffs, len(cx), rms


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
