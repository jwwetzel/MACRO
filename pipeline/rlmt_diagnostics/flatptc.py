"""Pure S2 flat-pair photon transfer: the gain measurement downstream uses.

WHY THIS MODULE EXISTS (and why it is not :mod:`rlmt_diagnostics.ptc`)
----------------------------------------------------------------------
:mod:`~rlmt_diagnostics.ptc` was built around one night of darks and star
fields.  Darks carry almost no Poisson signal, and star fields move between
exposures, so its two slopes could only BRACKET the gain ([0.60, 1.77]
e-/ADU on High Gain) and it reported read noise 0.0 from four railed fits.
The 2026-10-03 committee review (DE.F1/F2/F4, DS.F7) pointed out that the
archive holds what a photon-transfer curve actually needs and what S2 never
read: FLAT FIELDS — a uniformly illuminated detector, where two consecutive
frames differ by nothing except noise.  This module is the textbook
flat-pair PTC (Janesick 2001), written so that every known bias of the
estimator is either removed analytically or measured and carried:

* **Illumination drift between the two frames** (a twilight sky fades by
  several percent between consecutive exposures; a panel lamp flickers).
  The plain difference ``a - b`` then contains the flat-field pattern times
  the level change.  We difference ``a - r*b`` with ``r`` the measured
  level ratio, tile by tile, and propagate ``r`` through the variance
  EXACTLY (see :func:`effective_signal_and_variance`) instead of pretending
  it is 1.
* **Fixed pattern** (pixel response non-uniformity, dust, vignetting)
  cancels in the scaled difference because both frames see the same
  pattern; what survives of it would grow as signal SQUARED, so
  :func:`fit_ptc_line` also reports the quadratic coefficient as the
  diagnostic that it did cancel.
* **Sigma clipping** removes cosmic rays and stars in twilight flats, and
  in doing so trims the tails of the honest Gaussian too.  The variance of
  a Gaussian truncated at +/-k sigma is smaller than sigma^2 by a known
  factor; :data:`CLIP_VARIANCE_FACTOR` restores it (0.11% at 4 sigma — small,
  but it is a bias with a known sign, and it is removed, not ignored).
* **Bias level**.  Signal is level minus bias, so a bias error of db moves
  a fixed-intercept gain by db/S.  The adopted estimator therefore fits the
  flat points JOINTLY with bias/dark-pair points at zero signal, and the
  free-intercept slope (which does not care what the bias is) is reported
  beside it as the cross-check.
* **Pixel-to-pixel correlation** (inter-pixel capacitance, or a readout
  that mixes neighbours) lowers the per-pixel variance and so inflates a
  PTC gain.  :func:`diff_autocorrelation` measures the nearest-neighbour
  correlation of the difference image; the campaign carries the implied
  gain correction as a stated systematic.

THE BINNING QUESTION (committee disagreement D2)
------------------------------------------------
An IMX455 read out "2x2" can hand back either the SUM or the AVERAGE of
four native pixels.  The header EGAIN (0.2467 e-/ADU for the ASI at gain
100) is the NATIVE-pixel gain either way.  For a superpixel collecting N
electrons in total:

* SUM:      level = N/g,        variance = N/g^2      ->  K_eff = g
* AVERAGE:  level = N/(4g),     variance = N/(16 g^2) ->  K_eff = 4g

so the flat-pair PTC of binned frames answers the question by itself:
:func:`binning_verdict` compares the measured K_eff with g and 4g.  The
independent hot-pixel check (one saturated native pixel averaged with three
normal ones lands at (65535 + 3*pedestal)/4) lives in
:mod:`rlmt_diagnostics.badpix`.

Everything here is a pure function on numpy arrays / plain sequences: no
I/O, no database, no globals.  The campaign script owns file reading and
pairing; the unit tests own synthetic Poisson flats with known truth,
including the cases this estimator must survive (drifting level, vignetted
flat, stars, averaged binning).
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: Tile edge in pixels.  A tile is the unit over which the level ratio r is
#: measured and the difference variance taken: small enough that vignetting
#: and twilight gradients are flat across it, large enough (4,096 pixels)
#: that one tile's variance is known to sqrt(2/4095) = 2.2%.
TILE_PX = 64

#: Fraction of each frame edge excluded before tiling.  Edges carry amp
#: glow, vignetting roll-off and (on the reduced 2026 frames) crop seams.
EDGE_FRACTION = 0.06

#: Robust clip threshold (in MAD-sigmas) for pixels inside one tile.
CLIP_SIGMA = 4.0

#: Clip iterations (find the scale, then settle).
CLIP_ITERS = 2

#: MAD -> sigma for a Gaussian.
MAD_TO_SIGMA = 1.4826


def _truncated_gaussian_variance_factor(k: float) -> float:
    """Variance of a unit Gaussian truncated at +/-k, i.e. what a k-sigma
    clip leaves behind:  1 - 2 k phi(k) / (2 Phi(k) - 1)."""
    phi = math.exp(-0.5 * k * k) / math.sqrt(2.0 * math.pi)
    big_phi = 0.5 * (1.0 + math.erf(k / math.sqrt(2.0)))
    return 1.0 - 2.0 * k * phi / (2.0 * big_phi - 1.0)


#: A k-sigma clip of honest Gaussian noise UNDERESTIMATES its variance by
#: this factor (0.99893 at 4 sigma).  Every clipped variance in this module
#: is divided by it, so the clip costs outliers and nothing else.
CLIP_VARIANCE_FACTOR = _truncated_gaussian_variance_factor(CLIP_SIGMA)

#: A tile is usable only if this fraction of its pixels survives masking
#: and clipping.  A tile that loses a third of its pixels is a tile with a
#: star, a dust mote that moved, or a bad column in it.
MIN_TILE_FRACTION = 0.70

#: Flat pairs whose level ratio leaves this window are refused: beyond it
#: the "same scene" assumption has failed (a different lamp setting, a
#: cloud, twilight moving too fast for the exposure).
RATIO_WINDOW = (0.6, 1.0 / 0.6)

#: Level bins per pair.  A flat is not flat — vignetting spreads its tiles
#: over +/-10-20% in level — and that spread is free leverage on the PTC
#: slope, independent of any bias-level assumption.  Tiles are grouped into
#: this many equal-population level bins and each bin becomes one
#: (signal, variance) point.
BINS_PER_PAIR = 4

#: A level bin needs at least this many tiles to yield a point.
MIN_TILES_PER_BIN = 12

#: Tile-variance outliers beyond this many MAD-sigmas of their bin are
#: dropped before the bin mean (a star that moved between twilight flats
#: owns its tile's variance).  Symmetric, so it does not bias the mean of a
#: symmetric distribution; the fraction dropped is recorded.
TILE_OUTLIER_SIGMA = 4.0

#: Bootstrap resamples for the gain uncertainty (pairs are the resampling
#: unit — points inside a pair share its bias error and its ratio r).
N_BOOTSTRAP = 2000

#: Fixed seed: the bootstrap is part of a published number and must
#: regenerate bit-for-bit.
BOOTSTRAP_SEED = 20261003

#: Acceptance criterion from the committee synthesis (F-4): gain known to
#: this relative precision.
GAIN_TARGET_REL_ERR = 0.03

#: Binning verdict: the measured K_eff/g must sit within this relative
#: distance of 1 (sum) or n_native (average) to be called.
BINNING_VERDICT_TOL = 0.15


# --------------------------------------------------------------------------
# Tiling
# --------------------------------------------------------------------------
def tile_stack(img: np.ndarray, tile: int = TILE_PX,
               edge_fraction: float = EDGE_FRACTION) -> np.ndarray:
    """Cut an image into non-overlapping tiles: (n_tiles, tile*tile).

    The frame is trimmed by ``edge_fraction`` on every side, then by
    whatever remainder does not fill a whole tile (taken off the high-index
    edges).  Returns a float64 array with one ROW per tile — the layout
    every per-tile statistic below vectorises over.
    """
    a = np.asarray(img, dtype=np.float64)
    ny, nx = a.shape
    ey, ex = int(round(ny * edge_fraction)), int(round(nx * edge_fraction))
    a = a[ey:ny - ey, ex:nx - ex]
    ty, tx = a.shape[0] // tile, a.shape[1] // tile
    if ty == 0 or tx == 0:
        return np.empty((0, tile * tile))
    a = a[:ty * tile, :tx * tile]
    return (a.reshape(ty, tile, tx, tile).swapaxes(1, 2)
             .reshape(ty * tx, tile * tile))


def effective_signal_and_variance(var_d: np.ndarray, r: np.ndarray,
                                  mean_b: np.ndarray
                                  ) -> tuple[np.ndarray, np.ndarray]:
    """Propagate the level ratio r through a scaled-difference variance.

    With per-pixel variance  sigma^2(S) = RN^2 + S/K  and  S_a = r S_b:

        var(a - r b) = RN^2 (1 + r^2) + (r + r^2) S_b / K

    Dividing by (1 + r^2) puts this back in the single-frame form
    ``V = RN^2 + S_eff / K`` with

        V     = var(a - r b) / (1 + r^2)
        S_eff = S_b * r (1 + r) / (1 + r^2)

    For r = 1 this is the textbook ``var(a - b)/2`` at signal S.  Returns
    ``(S_eff, V)``.
    """
    r = np.asarray(r, dtype=np.float64)
    denom = 1.0 + r * r
    return (np.asarray(mean_b, dtype=np.float64) * r * (1.0 + r) / denom,
            np.asarray(var_d, dtype=np.float64) / denom)


def tile_pair_stats(a: np.ndarray, b: np.ndarray, bias=0.0,
                    valid: Optional[np.ndarray] = None,
                    tile: int = TILE_PX,
                    edge_fraction: float = EDGE_FRACTION) -> dict:
    """Per-tile photon-transfer statistics of one same-scene frame pair.

    Parameters
    ----------
    a, b
        The two frames (raw ADU, same shape).
    bias
        Bias level to subtract before anything else: a scalar or a frame
        (master bias / dark of the same exposure).  Pass 0 for the
        bias-agnostic use (sky pairs fitted in raw level).
    valid
        Optional boolean frame, True where a pixel may be used (bad-pixel
        and star masks).  Masked pixels are simply absent from every
        statistic.

    Returns
    -------
    dict of 1-D arrays, one entry per USABLE tile:

    * ``signal``  — S_eff, the Poisson-weighted bias-subtracted level (ADU);
    * ``var``     — V, the clip-corrected single-frame variance (ADU^2);
    * ``ratio``   — the tile's level ratio r = <a>/<b>;
    * ``n_pix``   — pixels behind the tile's variance;

    plus scalars ``n_tiles`` (all tiles) and ``n_used``.

    Method, per tile: r from the tile medians; d = a - r b; pixels beyond
    :data:`CLIP_SIGMA` MAD-sigmas of d are dropped (twice); r is then
    re-measured as the ratio of MEANS over the surviving pixels (a mean,
    not a median, because Poisson signal is a mean) and d recomputed; the
    variance of the survivors is divided by :data:`CLIP_VARIANCE_FACTOR`
    and propagated through :func:`effective_signal_and_variance`.
    """
    A = tile_stack(np.asarray(a, dtype=np.float64) - bias, tile, edge_fraction)
    B = tile_stack(np.asarray(b, dtype=np.float64) - bias, tile, edge_fraction)
    n_tiles = A.shape[0]
    empty = {"signal": np.array([]), "var": np.array([]),
             "ratio": np.array([]), "n_pix": np.array([]),
             "n_tiles": int(n_tiles), "n_used": 0}
    if n_tiles == 0:
        return empty
    keep = np.isfinite(A) & np.isfinite(B)
    if valid is not None:
        keep &= tile_stack(np.asarray(valid, dtype=np.float64),
                           tile, edge_fraction) > 0.5
    A = np.where(keep, A, np.nan)
    B = np.where(keep, B, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        medB = np.nanmedian(B, axis=1)
        r = np.nanmedian(A, axis=1) / medB
        for it in range(CLIP_ITERS + 1):
            d = A - r[:, None] * B
            if it == CLIP_ITERS:
                break
            med = np.nanmedian(d, axis=1)
            sig = MAD_TO_SIGMA * np.nanmedian(np.abs(d - med[:, None]), axis=1)
            bad = np.abs(d - med[:, None]) > CLIP_SIGMA * sig[:, None]
            # A zero MAD (constant tile) clips nothing: sig = 0 would
            # otherwise reject every pixel that is not exactly the median.
            bad &= (sig[:, None] > 0)
            A = np.where(bad, np.nan, A)
            B = np.where(bad, np.nan, B)
            # Ratio of MEANS over the survivors: the Poisson-correct level.
            r = np.nanmean(A, axis=1) / np.nanmean(B, axis=1)
        n_pix = np.sum(np.isfinite(d), axis=1)
        var_d = np.nanvar(d, axis=1, ddof=1) / CLIP_VARIANCE_FACTOR
        mean_b = np.nanmean(B, axis=1)
    sig_eff, var = effective_signal_and_variance(var_d, r, mean_b)
    ok = (np.isfinite(sig_eff) & np.isfinite(var) & np.isfinite(r)
          & (n_pix >= MIN_TILE_FRACTION * tile * tile)
          & (r > RATIO_WINDOW[0]) & (r < RATIO_WINDOW[1]))
    return {"signal": sig_eff[ok], "var": var[ok], "ratio": r[ok],
            "n_pix": n_pix[ok].astype(np.float64),
            "n_tiles": int(n_tiles), "n_used": int(ok.sum())}


def pair_points(stats: dict, n_bins: int = BINS_PER_PAIR,
                min_tiles: int = MIN_TILES_PER_BIN) -> list[dict]:
    """Collapse one pair's tiles into a few (signal, variance) points.

    Tiles are sorted on signal and cut into ``n_bins`` equal-population
    bins.  Inside each bin, tiles whose variance sits more than
    :data:`TILE_OUTLIER_SIGMA` MAD-sigmas from the bin median are dropped
    (a star that moved, a cosmic-ray shower), and the survivors give:

    * ``signal``   — mean signal of the bin (ADU);
    * ``var``      — mean tile variance (ADU^2);
    * ``var_err``  — standard error of that mean, from the tile-to-tile
      scatter (an EMPIRICAL error: it knows nothing about Gaussianity);
    * ``n_tiles``  — tiles used; ``n_dropped`` — tiles rejected as outliers;
    * ``ratio``    — mean level ratio of the bin's tiles.

    The variance of a bin whose tiles span a range of signal is the mean
    of a LINEAR function of signal, so pairing it with the mean signal is
    exact, not an approximation.  The ERROR, however, must not count that
    in-bin trend as scatter (it would inflate every error bar and push
    chi2/dof below 1 — a defect by standing rule 1): the pair's own
    variance-vs-signal trend is removed (a Theil-Sen-like median slope
    over the pair's tiles) before the outlier screen and the standard
    error are computed.
    """
    s = np.asarray(stats["signal"], dtype=np.float64)
    v = np.asarray(stats["var"], dtype=np.float64)
    r = np.asarray(stats["ratio"], dtype=np.float64)
    if s.size < min_tiles:
        return []
    n_bins = max(1, min(n_bins, s.size // min_tiles))
    order = np.argsort(s)
    # In-pair trend: slope between the medians of the lower and upper
    # signal halves (robust; zero when the pair has no level spread).
    half = s.size // 2
    lo, hi = order[:half], order[half:]
    ds = np.median(s[hi]) - np.median(s[lo])
    trend = ((np.median(v[hi]) - np.median(v[lo])) / ds) if ds > 0 else 0.0
    resid = v - trend * (s - np.median(s))
    out: list[dict] = []
    for chunk in np.array_split(order, n_bins):
        vv, ss, rr, ee = v[chunk], s[chunk], r[chunk], resid[chunk]
        med = np.median(ee)
        mad = MAD_TO_SIGMA * np.median(np.abs(ee - med))
        good = (np.abs(ee - med) <= TILE_OUTLIER_SIGMA * mad
                if mad > 0 else np.ones(vv.size, dtype=bool))
        if int(good.sum()) < min_tiles:
            continue
        out.append({
            "signal": float(ss[good].mean()),
            "var": float(vv[good].mean()),
            "var_err": float(ee[good].std(ddof=1) / math.sqrt(good.sum())),
            "n_tiles": int(good.sum()),
            "n_dropped": int((~good).sum()),
            "ratio": float(rr[good].mean()),
        })
    return out


def pair_read_noise(a: np.ndarray, b: np.ndarray,
                    valid: Optional[np.ndarray] = None,
                    tile: int = TILE_PX,
                    edge_fraction: float = EDGE_FRACTION) -> Optional[dict]:
    """Read noise from one bias (or shortest-dark) pair.

    Two zero-signal frames differ only by read noise, so the tile variance
    of ``a - b`` halved IS the read-noise variance.  No level ratio is
    applied (the level is the bias pedestal, and scaling a pedestal would
    inject its structure).  Returns the variance, its standard error from
    the tile scatter, the read noise in ADU with propagated error, the
    pedestal level and the tile count; None when no tile survives.
    """
    A = tile_stack(a, tile, edge_fraction)
    B = tile_stack(b, tile, edge_fraction)
    if A.shape[0] == 0:
        return None
    keep = np.isfinite(A) & np.isfinite(B)
    if valid is not None:
        keep &= tile_stack(np.asarray(valid, dtype=np.float64),
                           tile, edge_fraction) > 0.5
    d = np.where(keep, A - B, np.nan)
    with np.errstate(invalid="ignore"):
        for _ in range(CLIP_ITERS):
            med = np.nanmedian(d, axis=1)
            sig = MAD_TO_SIGMA * np.nanmedian(np.abs(d - med[:, None]), axis=1)
            bad = ((np.abs(d - med[:, None]) > CLIP_SIGMA * sig[:, None])
                   & (sig[:, None] > 0))
            d = np.where(bad, np.nan, d)
        n_pix = np.sum(np.isfinite(d), axis=1)
        var = np.nanvar(d, axis=1, ddof=1) / 2.0 / CLIP_VARIANCE_FACTOR
    ok = np.isfinite(var) & (n_pix >= MIN_TILE_FRACTION * tile * tile)
    if int(ok.sum()) < MIN_TILES_PER_BIN:
        return None
    var = var[ok]
    med = np.median(var)
    mad = MAD_TO_SIGMA * np.median(np.abs(var - med))
    good = (np.abs(var - med) <= TILE_OUTLIER_SIGMA * mad
            if mad > 0 else np.ones(var.size, dtype=bool))
    v = float(var[good].mean())
    v_err = float(var[good].std(ddof=1) / math.sqrt(good.sum()))
    rn = math.sqrt(max(v, 0.0))
    return {"var": v, "var_err": v_err, "read_noise_adu": rn,
            "read_noise_adu_err": (v_err / (2.0 * rn)) if rn > 0 else None,
            "level": float(np.nanmedian(np.where(keep, (A + B) / 2.0,
                                                 np.nan))),
            "n_tiles": int(good.sum())}


def diff_autocorrelation(a: np.ndarray, b: np.ndarray, bias=0.0,
                         tile: int = TILE_PX,
                         edge_fraction: float = EDGE_FRACTION) -> dict:
    """Nearest-neighbour correlation of a flat pair's difference image.

    Independent pixels give 0.  Inter-pixel capacitance, charge diffusion
    or a readout that mixes neighbours gives a POSITIVE correlation, and a
    per-pixel PTC then overstates the gain by about ``1/(1 - 4 rho)`` for
    an isotropic coupling (variance removed from each pixel reappears as
    covariance with its four neighbours).  Returns the lag-(0,1) and
    lag-(1,0) correlation coefficients, each the median over tiles of the
    tile's own lagged correlation of the scaled difference, and their
    standard errors from the tile scatter.  Tiles are screened with the
    same clip as the variance (outliers replaced by the tile median so a
    cosmic ray cannot fake a correlation).
    """
    af = np.asarray(a, dtype=np.float64) - bias
    bf = np.asarray(b, dtype=np.float64) - bias
    ny, nx = af.shape
    ey, ex = int(round(ny * edge_fraction)), int(round(nx * edge_fraction))
    af, bf = af[ey:ny - ey, ex:nx - ex], bf[ey:ny - ey, ex:nx - ex]
    ty, tx = af.shape[0] // tile, af.shape[1] // tile
    rho_x, rho_y = [], []
    for iy in range(ty):
        for ix in range(tx):
            sa = af[iy * tile:(iy + 1) * tile, ix * tile:(ix + 1) * tile]
            sb = bf[iy * tile:(iy + 1) * tile, ix * tile:(ix + 1) * tile]
            mb = np.median(sb)
            if not np.isfinite(mb) or mb == 0:
                continue
            d = sa - (np.median(sa) / mb) * sb
            med = np.median(d)
            sig = MAD_TO_SIGMA * np.median(np.abs(d - med))
            if not sig > 0:
                continue
            d = np.where(np.abs(d - med) > CLIP_SIGMA * sig, med, d) - med
            v = float((d * d).mean())
            rho_x.append(float((d[:, 1:] * d[:, :-1]).mean() / v))
            rho_y.append(float((d[1:, :] * d[:-1, :]).mean() / v))
    if len(rho_x) < MIN_TILES_PER_BIN:
        return {"rho_x": None, "rho_y": None, "rho_x_err": None,
                "rho_y_err": None, "n_tiles": len(rho_x)}
    rx, ry = np.array(rho_x), np.array(rho_y)
    n = rx.size
    return {"rho_x": float(rx.mean()), "rho_y": float(ry.mean()),
            "rho_x_err": float(rx.std(ddof=1) / math.sqrt(n)),
            "rho_y_err": float(ry.std(ddof=1) / math.sqrt(n)),
            "n_tiles": int(n)}


# --------------------------------------------------------------------------
# The fit
# --------------------------------------------------------------------------
def _wls(x: np.ndarray, y: np.ndarray, sigma: np.ndarray, order: int
         ) -> Optional[tuple[np.ndarray, np.ndarray, float]]:
    """Weighted polynomial least squares with honest sigmas.

    Returns (coefficients low->high order, covariance, chi2) or None when
    the design matrix is singular.  ``sigma`` are 1-sigma errors on y, so
    chi2 is a real chi-square and the covariance is NOT rescaled by it
    (standing rule 1: chi2/dof is reported, never absorbed).
    """
    w = 1.0 / np.square(sigma)
    X = np.vander(x, order + 1, increasing=True)
    XtW = X.T * w
    try:
        cov = np.linalg.inv(XtW @ X)
    except np.linalg.LinAlgError:
        return None
    beta = cov @ (XtW @ y)
    chi2 = float((w * (y - X @ beta) ** 2).sum())
    return beta, cov, chi2


def fit_ptc_line(points: Sequence[dict],
                 n_boot: int = N_BOOTSTRAP,
                 seed: int = BOOTSTRAP_SEED) -> Optional[dict]:
    """Fit  V = RN^2 + S/K  to photon-transfer points, with honest errors.

    Parameters
    ----------
    points
        Dicts with ``signal``, ``var``, ``var_err`` and ``pair_id``.
        Flat-pair points and (optionally) bias-pair points at signal 0 go
        in together: the bias pairs then pin the intercept with a real
        measurement instead of an extrapolation.

    Returns
    -------
    dict or None
        * ``gain`` (e-/ADU) = 1/slope and ``gain_err_formal`` from the
          weighted-least-squares covariance (meaningful only if the
          per-point errors are right — which ``chi2``/``dof`` tests);
        * ``gain_err_boot`` — the standard deviation of the gain over
          :data:`N_BOOTSTRAP` resamples of whole PAIRS (cluster bootstrap:
          points from one pair share its ratio and bias error), and
          ``gain_boot_bias`` — bootstrap mean minus point estimate, the
          estimator's measured small-sample bias (SIGNED);
        * ``intercept`` (ADU^2), ``intercept_err``, ``read_noise_adu`` and
          its error — the fit's own read noise, to be compared with the
          bias-pair value, never substituted for it;
        * ``chi2``, ``dof``, ``chi2nu`` — reported as measured.  A value
          far from 1 in EITHER direction is a defect of the error model
          and is stated, not hidden;
        * ``quad_coeff``, ``quad_coeff_err``, ``quad_z`` — the S^2 term of
          a three-parameter refit: residual fixed pattern or non-linearity
          shows up here.  ``gain_quad`` is the gain that refit implies at
          zero signal (the sensitivity of the answer to the model);
        * ``n_points``, ``n_pairs``, ``signal_lo``, ``signal_hi``.

        None for fewer than 3 points, fewer than 2 pairs, or a
        non-positive slope.
    """
    pts = [p for p in points
           if p.get("var_err") and p["var_err"] > 0
           and np.isfinite(p["signal"]) and np.isfinite(p["var"])]
    if len(pts) < 3:
        return None
    x = np.array([p["signal"] for p in pts], dtype=np.float64)
    y = np.array([p["var"] for p in pts], dtype=np.float64)
    s = np.array([p["var_err"] for p in pts], dtype=np.float64)
    ids = np.array([str(p["pair_id"]) for p in pts])
    uniq = np.unique(ids)
    if uniq.size < 2 or np.ptp(x) <= 0:
        return None
    lin = _wls(x, y, s, 1)
    if lin is None or lin[0][1] <= 0:
        return None
    beta, cov, chi2 = lin
    slope, inter = float(beta[1]), float(beta[0])
    dof = int(x.size - 2)
    gain = 1.0 / slope
    out = {
        "gain": gain,
        "gain_err_formal": float(math.sqrt(cov[1, 1]) / slope ** 2),
        "slope": slope, "intercept": inter,
        "intercept_err": float(math.sqrt(cov[0, 0])),
        "chi2": chi2, "dof": dof,
        "chi2nu": (chi2 / dof) if dof > 0 else None,
        "n_points": int(x.size), "n_pairs": int(uniq.size),
        "signal_lo": float(x.min()), "signal_hi": float(x.max()),
    }
    rn = math.sqrt(inter) if inter > 0 else None
    out["read_noise_adu"] = rn
    out["read_noise_adu_err"] = (out["intercept_err"] / (2.0 * rn)
                                 if rn else None)
    # Quadratic sensitivity: does an S^2 term want to exist?
    quad = _wls(x, y, s, 2) if x.size >= 4 else None
    if quad is not None and quad[0][1] > 0:
        qb, qcov, _ = quad
        qerr = float(math.sqrt(max(qcov[2, 2], 0.0)))
        out.update({"quad_coeff": float(qb[2]), "quad_coeff_err": qerr,
                    "quad_z": (float(qb[2]) / qerr) if qerr > 0 else None,
                    "gain_quad": float(1.0 / qb[1])})
    else:
        out.update({"quad_coeff": None, "quad_coeff_err": None,
                    "quad_z": None, "gain_quad": None})
    # Cluster bootstrap over pairs.
    rng = np.random.default_rng(seed)
    members = [np.flatnonzero(ids == u) for u in uniq]
    gains = []
    for _ in range(n_boot):
        pick = rng.integers(0, uniq.size, uniq.size)
        idx = np.concatenate([members[k] for k in pick])
        if np.ptp(x[idx]) <= 0:
            continue
        fit = _wls(x[idx], y[idx], s[idx], 1)
        if fit is None or fit[0][1] <= 0:
            continue
        gains.append(1.0 / fit[0][1])
    if len(gains) >= max(50, n_boot // 4):
        g = np.array(gains)
        out["gain_err_boot"] = float(g.std(ddof=1))
        out["gain_boot_bias"] = float(g.mean() - gain)
        out["n_boot_ok"] = int(g.size)
    else:
        out["gain_err_boot"] = None
        out["gain_boot_bias"] = None
        out["n_boot_ok"] = len(gains)
    return out


def fit_ptc_within_groups(points: Sequence[dict],
                          n_boot: int = N_BOOTSTRAP,
                          seed: int = BOOTSTRAP_SEED) -> Optional[dict]:
    """Fit  V = a_g + S/K  with a FREE additive term per flat group.

    WHY A SECOND ESTIMATOR.  The single-line fit assumes that a flat pair
    differs from a bias pair only by Poisson noise.  Sky flats break that:
    they contain faint stars, and between two exposures the telescope has
    moved, so each star is in the difference image twice (once positive,
    once negative).  The clip removes the bright cores; the faint ones
    remain as an ADDITIVE variance that does not scale with the twilight
    level.  A line pinned at the bias-pair variance then tilts to absorb
    it — on the AC4040 twilight flats that tilt is -0.9% in gain with a
    chi-square of 500 for 82 degrees of freedom, which is how it was
    found.

    This estimator is immune to any contamination that is constant within
    a group (one twilight sequence through one filter): each group gets
    its own intercept ``a_g`` and only the WITHIN-group variation of
    variance with signal — the sky brightening from frame to frame, and
    the vignetting spread inside each pair — determines the common slope.
    It is equally immune to a bias-level error (a shift of every S in the
    group).  The price: a group whose pairs all sit at one level (panel
    flats at a fixed lamp setting) carries no information here.

    Parameters
    ----------
    points
        Dicts with ``signal``, ``var``, ``var_err``, ``pair_id`` and
        ``group``.

    Returns
    -------
    dict or None
        The same gain/chi2/bootstrap keys as :func:`fit_ptc_line`
        (``dof`` = points - groups - 1), plus ``groups``: a list of
        ``{"group", "intercept", "intercept_err", "n_points",
        "signal_lo", "signal_hi"}`` — the per-group additive terms, to be
        compared with the read-noise variance (their excess over it IS the
        contamination).  None when fewer than 3 points carry within-group
        leverage or the slope is non-positive.
    """
    pts = [p for p in points
           if p.get("var_err") and p["var_err"] > 0
           and np.isfinite(p["signal"]) and np.isfinite(p["var"])]
    if len(pts) < 4:
        return None
    x = np.array([p["signal"] for p in pts], dtype=np.float64)
    y = np.array([p["var"] for p in pts], dtype=np.float64)
    sig = np.array([p["var_err"] for p in pts], dtype=np.float64)
    ids = np.array([str(p["pair_id"]) for p in pts])
    grp = np.array([str(p["group"]) for p in pts])

    def solve(sel: np.ndarray) -> Optional[tuple]:
        """Weighted within-group slope on a subset; None if no leverage."""
        w = 1.0 / sig[sel] ** 2
        xs, ys, gs = x[sel], y[sel], grp[sel]
        xt = np.empty_like(xs)
        yt = np.empty_like(ys)
        n_groups = 0
        for g in np.unique(gs):
            m = gs == g
            wm = w[m].sum()
            xt[m] = xs[m] - (w[m] * xs[m]).sum() / wm
            yt[m] = ys[m] - (w[m] * ys[m]).sum() / wm
            n_groups += 1
        sxx = float((w * xt * xt).sum())
        if sxx <= 0:
            return None
        slope = float((w * xt * yt).sum() / sxx)
        chi2 = float((w * (yt - slope * xt) ** 2).sum())
        return slope, sxx, chi2, n_groups

    full = np.ones(x.size, dtype=bool)
    res = solve(full)
    if res is None or res[0] <= 0:
        return None
    slope, sxx, chi2, n_groups = res
    dof = int(x.size - n_groups - 1)
    if dof < 1:
        return None
    gain = 1.0 / slope
    out = {"gain": gain,
           "gain_err_formal": float(1.0 / math.sqrt(sxx) / slope ** 2),
           "slope": slope, "intercept": None, "intercept_err": None,
           "read_noise_adu": None, "read_noise_adu_err": None,
           "chi2": chi2, "dof": dof, "chi2nu": chi2 / dof,
           "quad_coeff": None, "quad_coeff_err": None, "quad_z": None,
           "gain_quad": None,
           "n_points": int(x.size), "n_pairs": int(np.unique(ids).size),
           "signal_lo": float(x.min()), "signal_hi": float(x.max())}
    # Per-group additive terms.
    groups = []
    w = 1.0 / sig ** 2
    for g in np.unique(grp):
        m = grp == g
        wm = w[m].sum()
        xb = float((w[m] * x[m]).sum() / wm)
        yb = float((w[m] * y[m]).sum() / wm)
        a_g = yb - slope * xb
        # Error: the group mean's own error and the slope's lever arm.
        a_err = math.sqrt(1.0 / wm + (xb ** 2) / sxx)
        groups.append({"group": g, "intercept": a_g, "intercept_err": a_err,
                       "n_points": int(m.sum()),
                       "signal_lo": float(x[m].min()),
                       "signal_hi": float(x[m].max())})
    out["groups"] = groups
    # Cluster bootstrap over pairs (group labels travel with their pairs).
    rng = np.random.default_rng(seed)
    uniq = np.unique(ids)
    members = [np.flatnonzero(ids == u) for u in uniq]
    gains = []
    for _ in range(n_boot):
        pick = rng.integers(0, uniq.size, uniq.size)
        idx = np.concatenate([members[k] for k in pick])
        r = solve(idx)
        if r is None or r[0] <= 0:
            continue
        gains.append(1.0 / r[0])
    if len(gains) >= max(50, n_boot // 4):
        g_arr = np.array(gains)
        out["gain_err_boot"] = float(g_arr.std(ddof=1))
        out["gain_boot_bias"] = float(g_arr.mean() - gain)
        out["n_boot_ok"] = int(g_arr.size)
    else:
        out["gain_err_boot"], out["gain_boot_bias"] = None, None
        out["n_boot_ok"] = len(gains)
    return out


def additive_excess(groups: Sequence[dict], zero_var: float,
                    zero_var_err: float = 0.0) -> Optional[dict]:
    """How much additive variance the flat groups carry above read noise.

    ``groups`` is the list from :func:`fit_ptc_within_groups`;
    ``zero_var`` the read-noise variance from zero-signal pairs.  Returns
    the mean excess ``a_g - zero_var`` over groups, its standard error
    from the group-to-group scatter (combined with the zero-level error),
    the significance ``z``, and the group count.  A significant positive
    excess says the flats contain something a bias pair does not (stars,
    scattered light), and that any fit pinned at the bias variance is
    biased low in gain.  None for no groups.
    """
    a = np.array([g["intercept"] for g in groups], dtype=np.float64)
    if a.size == 0:
        return None
    ex = a - zero_var
    mean = float(ex.mean())
    if a.size > 1:
        err = float(ex.std(ddof=1) / math.sqrt(a.size))
    else:
        err = float(groups[0]["intercept_err"])
    err = float(math.hypot(err, zero_var_err))
    return {"excess": mean, "excess_err": err,
            "z": (mean / err) if err > 0 else None, "n_groups": int(a.size)}


def adopted_gain_error(fit: dict) -> float:
    """The statistical error adopted for a fitted gain.

    The LARGER of the formal and the pair-bootstrap error: the formal one
    is right only if the per-point errors are, the bootstrap one only if
    there are enough pairs to resample.  Taking the larger cannot
    understate either failure.
    """
    cands = [e for e in (fit.get("gain_err_formal"), fit.get("gain_err_boot"))
             if e is not None and np.isfinite(e)]
    return float(max(cands))


def full_scale_electrons(clip_adu: float, bias_adu: float, gain: float,
                         gain_err: float = 0.0, bias_err: float = 0.0
                         ) -> tuple[float, float]:
    """Effective full scale in electrons: (clip - bias) * K, with error.

    This is the largest signal the mode can RECORD, not necessarily the
    pixel full well: on High Gain the 12-bit-consistent clip arrives long
    before the GSENSE4040 well fills.  Error = quadrature of the gain and
    bias-level contributions.
    """
    span = float(clip_adu) - float(bias_adu)
    fs = span * gain
    return fs, float(math.hypot(span * gain_err, gain * bias_err))


def binning_verdict(k_eff: float, k_eff_err: float, native_gain: float,
                    n_native: int = 4,
                    tol: float = BINNING_VERDICT_TOL) -> dict:
    """Is an on-camera 2x2 bin a SUM or an AVERAGE of native pixels?

    ``native_gain`` is the unbinned e-/ADU (the header EGAIN, itself
    checked against an unbinned PTC when one exists).  A sum leaves the
    PTC gain at ``native_gain``; an average multiplies it by ``n_native``
    (see the module docstring).  Returns the measured ratio
    ``k_eff / native_gain`` with its error, the verdict (``"average"``,
    ``"sum"`` or ``"neither"``) and, for each hypothesis, how many of the
    measurement's own sigmas away it lies.
    """
    ratio = k_eff / native_gain
    ratio_err = k_eff_err / native_gain
    z_sum = (ratio - 1.0) / ratio_err if ratio_err > 0 else None
    z_avg = (ratio - n_native) / ratio_err if ratio_err > 0 else None
    if abs(ratio - n_native) <= tol * n_native:
        verdict = "average"
    elif abs(ratio - 1.0) <= tol:
        verdict = "sum"
    else:
        verdict = "neither"
    return {"ratio": float(ratio), "ratio_err": float(ratio_err),
            "verdict": verdict, "z_from_sum": z_sum, "z_from_average": z_avg,
            "n_native": int(n_native)}


def star_free_mask(mean_img: np.ndarray, tile: int = TILE_PX,
                   nsigma: float = 3.0, grow: int = 2) -> np.ndarray:
    """Boolean mask of sky pixels in a science frame (True = sky).

    For SKY-pair photon transfer (configurations with no flats on disk):
    the frame is cut into tiles, each tile's median and MAD-sigma define
    its sky, and pixels more than ``nsigma`` above it are flagged as
    sources and grown by ``grow`` pixels to cover their wings.

    The cut is made on the pair's MEAN image.  For two frames of equal
    noise the mean and the difference are statistically independent, so
    selecting pixels on the mean does not truncate — and therefore does
    not bias — the distribution of the difference whose variance is the
    measurement.
    """
    m = np.asarray(mean_img, dtype=np.float64)
    ny, nx = m.shape
    sky = np.ones_like(m, dtype=bool)
    for y0 in range(0, ny, tile):
        for x0 in range(0, nx, tile):
            sub = m[y0:y0 + tile, x0:x0 + tile]
            med = np.median(sub)
            sig = MAD_TO_SIGMA * np.median(np.abs(sub - med))
            sky[y0:y0 + tile, x0:x0 + tile] = sub <= med + nsigma * max(sig,
                                                                         1e-9)
    if grow > 0:
        src = ~sky
        grown = src.copy()
        for dy in range(-grow, grow + 1):
            for dx in range(-grow, grow + 1):
                if dy == 0 and dx == 0:
                    continue
                grown |= np.roll(np.roll(src, dy, axis=0), dx, axis=1)
        sky = ~grown
    return sky
