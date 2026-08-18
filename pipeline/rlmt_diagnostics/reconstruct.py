"""Pure S2 master-reconstruction logic: recover D and F from raw/reduced pairs.

THE EXPERIMENT (designed 2026-08-18)
------------------------------------
The archive's ``reduced/`` tree was produced by a pipeline we never audited.
But S0b linked 106,925 reduced frames back to their raw parents, and the
reduced headers describe the arithmetic (pyscope, "Calibration mode: ccd"):

    reduced = (raw - D) / F + PEDESTAL

with D the master dark(+bias) actually subtracted, F the flat actually
divided (often "flat correction NOT performed", i.e. F = 1), and PEDESTAL a
constant (1000 ADU) added back so the integer file has no negative pixels.
Rearranged per pixel, with r = reduced - PEDESTAL:

    raw = F * r + D

— a straight line, one per pixel, whose slope IS the flat value and whose
intercept IS the dark value at that pixel.  Many (raw, reduced) pairs of
DIFFERENT scenes sweep r through a range of values, so a least-squares line
per pixel *reconstructs the calibration frames the pipeline used* without
ever having seen them.  Era 47 (Andor iKon) is the control: the archive
holds its actual master bias/dark, so the reconstruction can be graded
against ground truth (residual map + RMS).

Robustness: the reduction ran astro-scrappy cosmic-ray cleaning, which
REPLACES isolated pixels in the reduced frame — those pixels break the
linear model for that one pair, so the fit iteratively rejects >k-sigma
residuals per pixel (they are few and uncorrelated across pairs).

All functions are pure numpy on stacked arrays; the campaign script owns
file I/O, pair selection and the pixel-grid subsampling bookkeeping.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Sigma-clip threshold for per-pixel residuals across pairs (cosmic-ray
#: replacements, satellite trails, stars that moved between pointings).
RECON_CLIP_SIGMA = 4.0

#: Clipping iterations for the per-pixel fit.
RECON_CLIP_ITERS = 3

#: Minimum surviving pairs for a pixel's (F, D) to be reported.
RECON_MIN_PAIRS = 8

#: Pixels whose raw value exceeds this fraction of the era's ceiling in a
#: given pair are masked for that pair — a clipped raw pixel no longer sits
#: on the line.
RECON_SAT_FRACTION = 0.95

#: Subsampled pixel grid: one central square + four corner squares, sizes in
#: pixels, corners inset by the margin (stay off the physical edge rows the
#: readout electronics distort).
RECON_CENTER_SIZE = 512
RECON_CORNER_SIZE = 128
RECON_EDGE_MARGIN = 32


def sample_regions(ny: int, nx: int,
                   center: int = RECON_CENTER_SIZE,
                   corner: int = RECON_CORNER_SIZE,
                   margin: int = RECON_EDGE_MARGIN,
                   ) -> list[tuple[str, slice, slice]]:
    """The five postage-stamp regions the experiment fits, as slices.

    Central ``center`` x ``center`` square plus four ``corner`` squares
    inset ``margin`` pixels from each edge.  Sizes shrink automatically for
    small sensors (never overlapping, never off-image); the region list is
    deterministic so re-runs and the report agree on geometry.
    """
    c = min(center, ny // 2, nx // 2)
    k = min(corner, ny // 4, nx // 4)
    m = min(margin, max(0, ny // 8), max(0, nx // 8))
    cy0, cx0 = (ny - c) // 2, (nx - c) // 2
    return [
        ("center", slice(cy0, cy0 + c), slice(cx0, cx0 + c)),
        ("corner_tl", slice(m, m + k), slice(m, m + k)),
        ("corner_tr", slice(m, m + k), slice(nx - m - k, nx - m)),
        ("corner_bl", slice(ny - m - k, ny - m), slice(m, m + k)),
        ("corner_br", slice(ny - m - k, ny - m), slice(nx - m - k, nx - m)),
    ]


#: Largest raw-minus-reduced shape difference (pixels, per axis) treated as
#: a reduction-pipeline crop worth aligning; bigger differences mean a
#: resampled/stacked product that is off-model for the experiment.
MAX_CROP_PX = 64

#: Patch size for the crop-offset search (big enough to contain real scene
#: structure, small enough that the 2-D offset scan stays instant).
CROP_PATCH = 256


def find_crop_offset(raw: np.ndarray, reduced: np.ndarray,
                     patch: int = CROP_PATCH) -> Optional[dict]:
    """Locate a cropped reduced frame inside its raw parent.

    The 2026 pipeline writes reduced frames a few rows/columns SMALLER than
    the raw sensor readout (overscan/edge trim: observed 4800x3211 ->
    4787x3193).  The crop offset is not recorded anywhere, so it is
    measured: a central ``patch`` of the reduced frame is slid over every
    feasible integer offset (the shape difference caps the search) and the
    offset minimizing the median absolute pixel difference wins — valid
    because reduction preserves pixel identity up to a smooth D/F transform,
    so the true alignment is dramatically better-matched than any shifted
    one.  Returns ``{"dy", "dx", "mad"}`` (reduced[0,0] sits at
    raw[dy, dx]), or None when the shapes do not describe a small crop.
    """
    ry, rx = raw.shape
    ny, nx = reduced.shape
    oy, ox = ry - ny, rx - nx
    if not (0 <= oy <= MAX_CROP_PX and 0 <= ox <= MAX_CROP_PX
            and (oy, ox) != (0, 0)):
        return None
    p = min(patch, ny // 2, nx // 2)
    cy, cx = (ny - p) // 2, (nx - p) // 2
    ref = np.asarray(reduced[cy:cy + p, cx:cx + p], dtype=np.float64)
    ref = ref - np.median(ref)          # remove DC so D cannot bias the MAD
    best = None
    for dy in range(oy + 1):
        for dx in range(ox + 1):
            win = np.asarray(raw[cy + dy:cy + dy + p, cx + dx:cx + dx + p],
                             dtype=np.float64)
            mad = float(np.median(np.abs((win - np.median(win)) - ref)))
            if best is None or mad < best["mad"]:
                best = {"dy": dy, "dx": dx, "mad": mad}
    return best


def fit_pixel_lines(reduced: np.ndarray, raw: np.ndarray,
                    sat_adu: Optional[float] = None) -> dict:
    """Per-pixel robust straight-line fit  raw = F * reduced + D.

    Parameters
    ----------
    reduced, raw
        Float arrays of shape ``(n_pairs, n_pixels)``: the SAME pixel across
        many pairs.  ``reduced`` must already be pedestal-subtracted.
    sat_adu
        Raw values >= ``RECON_SAT_FRACTION * sat_adu`` are masked per pair
        (clipped pixels are off the line by construction).

    Returns
    -------
    dict of 1-D arrays over pixels: ``F`` (slope), ``D`` (intercept),
    ``n_used`` (surviving pairs), ``rms`` (residual RMS, ADU).  Pixels with
    fewer than :data:`RECON_MIN_PAIRS` survivors get NaN F/D.

    Vectorized: the normal equations are computed with masked sums over the
    pair axis; :data:`RECON_CLIP_ITERS` rounds of residual clipping at
    :data:`RECON_CLIP_SIGMA` * (per-pixel residual RMS) remove cosmic-ray
    replacements and moving sources.
    """
    r = np.asarray(reduced, dtype=np.float64)
    y = np.asarray(raw, dtype=np.float64)
    if r.shape != y.shape or r.ndim != 2:
        raise ValueError("reduced and raw must share shape (n_pairs, n_pix)")
    mask = np.isfinite(r) & np.isfinite(y)
    if sat_adu is not None:
        mask &= y < RECON_SAT_FRACTION * sat_adu

    F = np.full(r.shape[1], np.nan)
    D = np.full(r.shape[1], np.nan)
    rms = np.full(r.shape[1], np.nan)

    # ---- Robust seeding pass --------------------------------------------
    # A cosmic-ray replacement moves a sample along the *x* axis (the
    # reduced value was rewritten), which makes it a HIGH-LEVERAGE outlier:
    # a least-squares line bends toward it and its own residual comes out
    # small, so residual clipping alone cannot find it.  Seed with a
    # median-of-slopes estimator instead (slope through each sample from
    # the pixel's median point; take the median), which one leverage point
    # cannot move, then clip against THAT line before any least squares.
    rn = np.where(mask, r, np.nan)
    yn = np.where(mask, y, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        xmed = np.nanmedian(rn, axis=0)
        ymed = np.nanmedian(yn, axis=0)
        dx = rn - xmed
        s = (yn - ymed) / dx
        # Samples sitting at the median x carry no slope information (and
        # divide by ~0); exclude them from the slope median only.
        s[np.abs(dx) < 1e-6] = np.nan
        slope0 = np.nanmedian(s, axis=0)
        inter0 = ymed - slope0 * xmed
        res0 = np.abs(yn - (slope0 * rn + inter0))
        mad0 = 1.4826 * np.nanmedian(res0, axis=0)
        keep = res0 <= RECON_CLIP_SIGMA * np.maximum(mad0, 1e-9)
    seeded = np.isfinite(slope0)
    # Only apply the seed clip where the seed line existed; elsewhere keep
    # the original mask (the LSQ pass below handles those pixels).
    mask = np.where(seeded[None, :], mask & keep, mask)

    for it in range(RECON_CLIP_ITERS + 1):
        n = mask.sum(axis=0)                          # pairs per pixel
        ok = n >= RECON_MIN_PAIRS
        if not ok.any():
            break
        w = mask.astype(np.float64)
        sw = np.maximum(n, 1)
        rm = (w * r).sum(axis=0) / sw                 # mean reduced
        ym = (w * y).sum(axis=0) / sw                 # mean raw
        rc = (r - rm) * w                             # centered, masked
        yc = (y - ym) * w
        srr = (rc * rc).sum(axis=0)
        sry = (rc * yc).sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            slope = np.where(srr > 0, sry / srr, np.nan)
        inter = ym - slope * rm
        resid = (y - (slope * r + inter)) * w
        with np.errstate(invalid="ignore"):
            rr = np.sqrt((resid ** 2).sum(axis=0) / np.maximum(n - 2, 1))
        F = np.where(ok, slope, np.nan)
        D = np.where(ok, inter, np.nan)
        rms = np.where(ok, rr, np.nan)
        if it == RECON_CLIP_ITERS:
            break
        # Clip: drop pair samples whose residual exceeds k * a ROBUST
        # per-pixel scale.  The scale is MAD-based (median |residual| over
        # the surviving pairs, Gaussian-calibrated by 1.4826), because a
        # plain RMS is inflated by the very outlier being tested — a single
        # large cosmic-ray replacement among ~20 pairs can hide inside an
        # RMS-derived threshold.  Zero-scale pixels clip nothing (floor).
        absres = np.where(mask, np.abs(y - (slope * r + inter)), np.nan)
        with np.errstate(invalid="ignore"):
            mad_scale = 1.4826 * np.nanmedian(absres, axis=0)
        thr = RECON_CLIP_SIGMA * np.maximum(mad_scale, 1e-9)
        newmask = mask & (np.abs(y - (slope * r + inter)) <= thr)
        if newmask.sum() == mask.sum():
            break                                     # converged: no change
        mask = newmask

    n_final = mask.sum(axis=0)
    return {"F": F, "D": D, "n_used": n_final.astype(np.int32), "rms": rms}


def summarize_reconstruction(F: np.ndarray, D: np.ndarray,
                             rms: np.ndarray) -> dict:
    """Headline statistics for one era's reconstruction (report + DB row).

    Medians and robust scatter (MAD-sigma) of the recovered flat and dark
    over the fitted pixels, the fraction of pixels that produced a fit, and
    the median residual RMS.  A flat median near 1.0 with small scatter =
    "flat correction NOT performed" confirmed from pixels, not comments.
    """
    def med_mad(x):
        v = x[np.isfinite(x)]
        if v.size == 0:
            return float("nan"), float("nan")
        med = float(np.median(v))
        return med, float(1.4826 * np.median(np.abs(v - med)))

    f_med, f_mad = med_mad(F)
    d_med, d_mad = med_mad(D)
    return {
        "flat_median": f_med, "flat_mad_sigma": f_mad,
        "dark_median": d_med, "dark_mad_sigma": d_mad,
        "fit_fraction": float(np.isfinite(F).mean()),
        "rms_median": float(np.nanmedian(rms)) if np.isfinite(rms).any()
                      else float("nan"),
    }


def flat_dark_correlation(F: np.ndarray, D: np.ndarray) -> float:
    """Pearson correlation of the per-pixel F and D estimates.

    THE DEGENERACY DIAGNOSTIC: the per-pixel line  raw = F*r + D  only
    separates slope from intercept when the pairs sweep r through a real
    range of levels.  With poor level diversity the estimates trade off
    against each other (a too-high F is compensated by a too-low D), which
    shows up as a strong NEGATIVE F-D correlation and a wide per-pixel F
    scatter — era 47's fits (36 near-constant-sky pairs) measure corr
    ~ -0.4 with 9% F scatter, so its per-pixel flat values are degeneracy
    noise and only the MEDIAN F carries meaning.  Well-conditioned eras
    (76/80, ~1% scatter) measure |corr| ~ 0.  Storing this number per era
    lets a reader tell which regime each verdict lives in.  Returns NaN
    when fewer than 10 pixels have finite F and D.
    """
    f = np.asarray(F, dtype=np.float64).ravel()
    d = np.asarray(D, dtype=np.float64).ravel()
    ok = np.isfinite(f) & np.isfinite(d)
    if ok.sum() < 10:
        return float("nan")
    return float(np.corrcoef(f[ok], d[ok])[0, 1])


def recon_verdict(flat_median: Optional[float], pedestal: Optional[float],
                  truth_rms: Optional[float]) -> str:
    """One-line verdict for an era's reconstruction table row.

    CAREFUL about what an identity fit proves: F ~ 1 with D ~ 0 and no
    pedestal demonstrates ONLY that the two trees hold the same pixels —
    no dark/flat applied BETWEEN them.  It says nothing about whether that
    shared content is raw or calibrated (era 79's review lesson: both of
    its trees turned out to hold the SAME already-calibrated products, so
    the earlier "reduced = uncalibrated raw copy" reading had the direction
    exactly backwards).  Which tree is miscast is a manifest/header
    question the caller must answer separately; this function no longer
    pretends the fit answers it.
    """
    if flat_median is None:
        return "unfittable (no readable same-geometry pairs)"
    if abs(flat_median - 1.0) < 0.005 and (pedestal or 0) == 0:
        return ("IDENTITY: both trees hold the same pixels "
                "(no relative calibration between them)")
    parts = ["dark subtracted"]
    parts.append("flat applied" if abs(flat_median - 1.0) > 0.02
                 else "no flat (F = 1)")
    if truth_rms is not None:
        parts.append(f"matches archived master to {truth_rms:.1f} ADU RMS")
    return ", ".join(parts)


def residual_vs_truth(D_recon: np.ndarray, truth: np.ndarray) -> dict:
    """Grade a reconstructed dark against the archive's actual master.

    The reconstruction and the master may sit on different DC offsets (the
    master carries the bias level; D absorbs whatever constant the pipeline
    used), so the comparison removes the MEDIAN offset first and then
    reports the structural agreement: robust RMS of (D_recon - truth -
    offset), plus the offset itself (a number worth reporting — it is the
    pedestal/bias bookkeeping of the unaudited pipeline made visible).
    """
    d = np.asarray(D_recon, dtype=np.float64)
    t = np.asarray(truth, dtype=np.float64)
    ok = np.isfinite(d) & np.isfinite(t)
    if not ok.any():
        return {"offset": float("nan"), "resid_rms": float("nan"),
                "resid_mad_sigma": float("nan"), "n_pix": 0}
    resid = d[ok] - t[ok]
    offset = float(np.median(resid))
    r = resid - offset
    return {
        "offset": offset,
        "resid_rms": float(np.sqrt(np.mean(r ** 2))),
        "resid_mad_sigma": float(1.4826 * np.median(np.abs(r))),
        "n_pix": int(ok.sum()),
    }
