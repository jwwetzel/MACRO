"""Line identification and the fixed grism dispersion (findings D1, G-1).

THE QUESTION THIS MODULE SETTLES
--------------------------------
Two readings of the same high-resolution-grism spectra were on file.  The
first-generation code (``wavelength.py`` v1) read the two telluric dips
redward of Halpha as the O2 B and A bands and adopted ~1.59 A/px.  The
observational astronomer (OA.E1) read the same dips as TiO 6651 and O2-B
and derived ~0.47 A/px — a factor 3.4, which is the difference between
R ~ 200 and R ~ 2,000 and decides whether line-profile science is
possible at all.  The two readings differ by exactly the O2 A/B lever
ratio (3.387), which is why a two-dip pattern on a COOL star cannot tell
them apart: an M giant's TiO band heads supply decoys at every spacing.

A HOT star can.  An A or B star has a clean continuum with a handful of
known stellar lines (the Balmer series, He I), and the atmosphere stamps
the same telluric bands on it.  Under the wrong dispersion the predicted
positions of those lines fall on featureless continuum.

So this module does two things, in this order:

1. :func: — a HYPOTHESIS TEST.  Given a candidate
   dispersion (the adopted seed, or one of the rival values that were on
   file), it tries every deep feature as Halpha with both signs and
   counts how many of the a-priori lines land on a real feature.  It is
   run once per rival, on the same features, with the same tolerance.
   The truth matches nearly every line the star can show; a wrong
   hypothesis matches its anchor plus chance.
2. :func: — the precise solution.  Line centres are
   measured individually (:func:, Gaussian + local linear
   continuum; :func: for band heads) and fitted with ONE
   dispersion polynomial per (grism, mechanical epoch), with a free zero
   point per frame — the only thing a slitless frame is allowed to have
   of its own (TE.F3, PH.P7).

A NOTE ON WHAT FAILED (kept because it is evidence).  Three fully blind
statistics were tried and withdrawn on 2026-10-03: a count of "confirmed"
lines over a dispersion grid (noise dips at the trace ends confirmed
anything), a template cosine similarity (it PREFERRED ~1.65 A/px on the
high-resolution grism, for the same reason version 1 did), and an
unseeded pair-assignment search (a tolerance loose enough to survive the
grism's 20% non-linearity matched 7 of 11 lines for almost any scale).  The cause is
physical and worth knowing: the hrg spectrum of every hot star carries
two telluric features almost symmetric about Halpha — the O2 B band head
at +304 A and the O2 gamma band head at -286 A — so any statistic that
rewards "a band at the right distance" is satisfied on either side and at
either scale.  What breaks the symmetry is the set of narrow STELLAR
lines (He I 5876/6678/7065, Si II 6347/6371, Ne I 6402), whose spacings
no wrong scale reproduces, and the absence of Hbeta.  So the question
is put the other way round: each candidate value is given its best
chance and scored on matched lines, at a tolerance that cannot confuse
values a factor 3.4 apart.

WAVELENGTH REFERENCES, AND WHAT IS PRIMARY
------------------------------------------
* **Primary: stellar lines.**  Laboratory air wavelengths; a star's
  radial velocity shifts all of them by the same factor (1 + v/c), which
  the per-frame zero point absorbs to first order and which changes the
  fitted dispersion by v/c <= 4e-4 (0.04%) even for Spica's 120 km/s
  orbit — forty times below the 1-2% acceptance.
* **Secondary: telluric band edges.**  The O2 bands are not lines at this
  resolution; what is measurable is the half-depth point of the sharp
  blue edge of each band (the R-branch head).  Its EFFECTIVE wavelength
  depends on the line spread function, so it is not typed from a table:
  it is *measured* on the stellar scale (``g_telluric_lambda``) and then
  used as a secondary standard.  The band-head wavelengths below are the
  a-priori values the blind scan uses and the measured values are checked
  against.

Everything here is pure numpy on 1-D arrays; no file or pixel access.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .trace import running_median

# --------------------------------------------------------------------------
# The line list (air wavelengths, Angstrom)
# --------------------------------------------------------------------------
HALPHA_A = 6562.80


@dataclass(frozen=True)
class Line:
    """One reference feature.

    ``kind``: 'stellar' (a line that moves with the star), 'telluric'
    (fixed in the observatory frame).  ``shape``: 'line' (Gaussian-like
    absorption) or 'edge' (band head: sharp blue edge, shaded to the
    red).  ``spt``: spectral classes in which the stellar line is strong
    enough to measure ('B', 'A'; telluric features use '*').
    """
    name: str
    wave: float
    kind: str
    shape: str
    spt: str
    weight: float = 1.0      # prior depth relative to Halpha (blind scan)


#: A-priori line list.  Stellar wavelengths are laboratory air values
#: (NIST ASD).  Telluric entries are the R-branch band HEADS (the
#: conventional 6276.6 / 6867.2 / 7593.7 A) — a-priori only: the
#: effective half-depth wavelengths at this instrument's resolution are
#: measured and stored in ``g_telluric_lambda``.
LINES = (
    # Balmer series
    Line("Hepsilon", 3970.07, "stellar", "line", "AB", 0.8),
    Line("Hdelta", 4101.73, "stellar", "line", "AB", 1.0),
    Line("Hgamma", 4340.46, "stellar", "line", "AB", 1.0),
    Line("Hbeta", 4861.32, "stellar", "line", "AB", 1.0),
    Line("Halpha", HALPHA_A, "stellar", "line", "AB", 1.0),
    # He I, Ne I (B stars) and Si II (late B / early A)
    Line("HeI5876", 5875.62, "stellar", "line", "B", 0.25),
    Line("SiII6347", 6347.11, "stellar", "line", "B", 0.10),
    Line("SiII6371", 6371.37, "stellar", "line", "B", 0.06),
    Line("NeI6402", 6402.25, "stellar", "line", "B", 0.06),
    Line("HeI6678", 6678.15, "stellar", "line", "B", 0.25),
    Line("HeI7065", 7065.19, "stellar", "line", "B", 0.15),
    # O I triplet (blend centroid; A stars)
    Line("OI7774", 7773.4, "stellar", "line", "A", 0.3),
    # Telluric band heads
    Line("O2gamma", 6276.6, "telluric", "edge", "*", 0.08),
    Line("O2B", 6867.2, "telluric", "edge", "*", 0.35),
    Line("O2A", 7593.7, "telluric", "edge", "*", 0.8),
)

#: Lines each grism is asked about.  The low-resolution unit cannot
#: resolve the weak He I / Si II / Ne I lines (depth < 2% at its ~15 A
#: resolution).  Hbeta IS offered to the high-resolution unit although it
#: lies outside that unit's true range: a hypothesis that needs Hbeta on
#: the chip (1.59 A/px does) must be free to claim it, and to fail.
GRISM_LINES = {
    "hrg": ("Hbeta", "HeI5876", "O2gamma", "SiII6347", "SiII6371",
            "NeI6402", "Halpha", "HeI6678", "O2B", "HeI7065", "O2A",
            "OI7774"),
    "lrg": ("Hepsilon", "Hdelta", "Hgamma", "Hbeta", "Halpha", "O2B",
            "O2A"),
}

#: Where the MINIMUM of a band sits redward of its head at this
#: resolution (A) — used only by the blind identification, which works
#: from feature minima (R-branch maxima; a-priori, +/-2 A, absorbed by
#: the identification tolerance).
BAND_MIN_OFFSET_A = {"O2gamma": 2.5, "O2B": 3.0, "O2A": 6.0}

LINE_BY_NAME = {ln.name: ln for ln in LINES}

#: Full extent (A) of the O2 bands redward of their heads (R + P
#: branches) — the blind-scan template draws a band over this span.
BAND_SPAN_A = {"O2gamma": 40.0, "O2B": 70.0, "O2A": 80.0}

#: Span (A) redward of a band head within which the band FLOOR is sought
#: when measuring the edge: the R-branch, between the head and the band
#: origin (O2-B: 6867 -> 6884; O2-A: 7594 -> 7621; gamma: 6277 -> 6288).
EDGE_FLOOR_SPAN_A = {"O2gamma": 12.0, "O2B": 16.0, "O2A": 22.0}

#: Continuum window (px) of the running-median normalisation.  Wide
#: enough that the Stark wings of an A-star Balmer line do not lift
#: their own continuum on the high-resolution grism.
CONT_WIN_PX = 301

#: Dispersion bracket (|A/px|) the blind identification may return.  It
#: contains both disputed hrg values (0.47, 1.59), the three lrg values
#: on file (1.1, 1.9-2.2, 3.2) and a factor of two of margin.
SCAN_DISP_MIN = 0.20
SCAN_DISP_MAX = 4.50

#: Blind identification: feature detection floor, and how many of the
#: most prominent features enter the search.
ID_MIN_PROMINENCE = 0.012
ID_DETECT_SIGMA = 5.0
ID_MAX_FEATURES = 14

#: Matching tolerance of the seeded identification (px + fraction of the
#: pixel distance from the Halpha anchor) and of its quadratic
#: refinement; how many of the deepest features are tried as Halpha.
ID_TOL_PX = 4.0
ID_TOL_FRAC = 0.025
ID_REFINE_TOL_PX = 2.5
ID_REFINE_TOL_FRAC = 0.002
ID_N_ANCHORS = 6

#: Dispersion HYPOTHESES for :func:, as (a, b) with
#: x - x_Halpha = a (lambda - 6562.8) + b (lambda - 6562.8)^2 in the
#: red-toward-positive sense.  "adopted" entries are SEEDS read by eye
#: off two frames (eta Hya hrg 2025-12-17, Phecda lrg 2025-01-23); they
#: only tell the matcher where to look — the adopted solution is the
#: least-squares fit of solve_dispersion to the measured centres.
#: The other entries are the values that were on file before 2026-10-03
#: and are tested, not used.
SEEDS = {"hrg": (2.19, -2.6e-4), "lrg": (0.424, -1.14e-5)}
RIVALS = {"hrg": {"1.59 A/px (v1 code)": (1.0 / 1.59, 0.0)},
          "lrg": {"1.1 A/px (stored mode)": (1.0 / 1.1, 0.0),
                  "3.2 A/px (S2c factor 2)": (1.0 / 3.2, 0.0)}}

#: Search window of the seeded precise measurement (runner).
SCAN_WIN_PX = 4.0
SCAN_WIN_FRAC = 0.03


# --------------------------------------------------------------------------
# Continuum normalisation
# --------------------------------------------------------------------------
def normalize(flux: np.ndarray, win: int = CONT_WIN_PX,
              floor_frac: float = 0.15):
    """(normalised, continuum, usable): flux divided by a running-median
    continuum.

    The median is a pseudo-continuum that sits slightly BELOW the true
    one inside wide lines; that depresses depths but not centres, which
    are all this module uses.  ``usable`` is False where the continuum
    falls below ``floor_frac`` of its maximum (the ends of the trace:
    ratios against nearly nothing measure noise).
    """
    finite = np.isfinite(flux)
    if not finite.any():
        nan = np.full(len(flux), np.nan)
        return nan, nan, finite
    # Fill gaps and the regions beyond the trace ends with the NEAREST
    # valid value before filtering.  (Filling with the global median, as
    # v1 did, makes the running median near a trace end a mixture of
    # real spectrum and a constant, and the normalised spectrum then
    # ramps away from 1 over half a window — a fake broad "absorption"
    # at each end that the feature finder duly reported.)
    idx = np.where(finite)[0]
    nearest = idx[np.clip(np.searchsorted(idx, np.arange(len(flux))), 0,
                          len(idx) - 1)]
    left = idx[np.clip(np.searchsorted(idx, np.arange(len(flux)),
                                       side="right") - 1, 0, len(idx) - 1)]
    pick = np.where(np.abs(nearest - np.arange(len(flux)))
                    <= np.abs(left - np.arange(len(flux))), nearest, left)
    fill = flux[pick]
    cont = running_median(fill, win)
    usable = finite & (cont > floor_frac * np.nanmax(cont[finite]))
    # The half-window at each end of the trace has a one-sided continuum
    # estimate; it is not used.
    usable[:idx[0] + win // 4] = False
    usable[idx[-1] - win // 4 + 1:] = False
    with np.errstate(divide="ignore", invalid="ignore"):
        norm = np.where(usable, flux / cont, np.nan)
    return norm, cont, usable


def smooth(y: np.ndarray, sigma_px: float) -> np.ndarray:
    """Gaussian smoothing that ignores NaNs (normalised convolution)."""
    from scipy.ndimage import gaussian_filter1d
    ok = np.isfinite(y)
    num = gaussian_filter1d(np.where(ok, y, 0.0), sigma_px, mode="nearest")
    den = gaussian_filter1d(ok.astype(float), sigma_px, mode="nearest")
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0.3, num / den, np.nan)


# --------------------------------------------------------------------------
# Feature finding
# --------------------------------------------------------------------------
def absorption_candidates(norm: np.ndarray, n: int = 4,
                          sigma_px: float = 3.0,
                          min_sep_px: int = 40) -> list[dict]:
    """The ``n`` deepest distinct absorption features of a normalised
    spectrum, as [{'x', 'depth'}, ...] deepest first.

    Depth is measured on the Gaussian-smoothed spectrum (single-pixel
    noise spikes cannot win); two candidates closer than ``min_sep_px``
    are one feature.
    """
    d = 1.0 - smooth(norm, sigma_px)
    d = np.where(np.isfinite(d), d, -np.inf)
    out = []
    work = d.copy()
    for _ in range(n):
        i = int(np.argmax(work))
        if not np.isfinite(work[i]) or work[i] <= 0:
            break
        out.append({"x": float(i), "depth": float(work[i])})
        lo, hi = max(0, i - min_sep_px), min(len(work), i + min_sep_px + 1)
        work[lo:hi] = -np.inf
    return out


def depth_noise(norm: np.ndarray, sigma_px: float = 3.0) -> float:
    """Noise of the SMOOTHED depth spectrum (robust): the MAD of the
    difference between the smoothed spectrum and a wider smoothing of
    itself — what a featureless stretch scatters by at feature scale."""
    d = smooth(norm, sigma_px) - smooth(norm, 6.0 * sigma_px)
    d = d[np.isfinite(d)]
    if len(d) < 50:
        return float("nan")
    return float(1.4826 * np.median(np.abs(d - np.median(d))))


def find_features(norm: np.ndarray, sigma_px: float = 2.5,
                  max_n: int = ID_MAX_FEATURES) -> list[dict]:
    """Absorption features of a normalised spectrum for the blind
    identification: local minima whose PROMINENCE exceeds both
    ``ID_MIN_PROMINENCE`` and ``ID_DETECT_SIGMA`` x the smoothed noise.
    Returns [{'x', 'depth'}, ...], the ``max_n`` most prominent, sorted
    by position.  Prominence (not raw depth) is what makes a line on the
    shoulder of a band, or on a sloping continuum, count as a feature.
    """
    from scipy.signal import find_peaks
    sm = smooth(norm, sigma_px)
    d = np.where(np.isfinite(sm), 1.0 - sm, 0.0)
    noise = depth_noise(norm, sigma_px)
    floor = max(ID_MIN_PROMINENCE,
                ID_DETECT_SIGMA * (noise if np.isfinite(noise) else 0.0))
    pk, prop = find_peaks(d, prominence=floor, width=2)
    order = np.argsort(prop["prominences"])[::-1][:max_n]
    out = [{"x": float(pk[i]), "depth": float(prop["prominences"][i])}
           for i in order]
    out.sort(key=lambda f: f["x"])
    return out


def _id_wave(ln: Line) -> float:
    """Wavelength at which a feature MINIMUM is expected for a line."""
    return ln.wave + BAND_MIN_OFFSET_A.get(ln.name, 0.0)


def _count_matches(xs: np.ndarray, waves: np.ndarray, coeffs,
                   pivot_px: float, tol0: float, tol_frac: float):
    """Greedy one-to-one matching of features to lines under a model
    x(lambda) = polyval(coeffs, lambda - HALPHA_A).  Tolerance grows with
    pixel distance from ``pivot_px``.  Returns [(i_feature, i_line,
    residual_px), ...]."""
    pred = np.polyval(coeffs, waves - HALPHA_A)
    pairs = []
    for m, xp in enumerate(pred):
        j = int(np.argmin(np.abs(xs - xp)))
        r = xs[j] - xp
        if abs(r) <= tol0 + tol_frac * abs(xp - pivot_px):
            pairs.append((j, m, float(r)))
    # One feature may serve one line only: keep the closest claim.
    best = {}
    for j, m, r in pairs:
        if j not in best or abs(r) < abs(best[j][2]):
            best[j] = (j, m, r)
    return sorted(best.values(), key=lambda t: t[1])


def identify_lines(features: Sequence[dict], grism: str, spt: str,
                   seed=None, pixel_scale: float = 1.0,
                   sign: Optional[int] = None) -> Optional[dict]:
    """Identify the features of a hot-star spectrum under ONE stated
    hypothesis about the dispersion, and say how well it does.

    ``seed`` = (a, b): the hypothesis, as px/A at Halpha and px/A^2, in
    the red-toward-positive sense (default: ``SEEDS[grism]``).  Nothing
    else is assumed: every one of the ``ID_N_ANCHORS`` deepest features
    is tried as Halpha, with both signs of the dispersion.  For each
    (anchor, sign) the list's lines are predicted and matched one-to-one
    to features within ``ID_TOL_PX`` + ``ID_TOL_FRAC`` x lever arm; the
    best (anchor, sign) is refined with a free quadratic through its
    matches and re-counted at ``ID_REFINE_TOL_PX``.

    This is a HYPOTHESIS TEST, not a fit: call it once per candidate
    dispersion and compare the scores.  A true hypothesis matches nearly
    every line the star can show, with sub-pixel residuals after the
    refinement; a false one matches the anchor plus chance.  The
    tolerance (2.5% of lever arm) is what lets one seed serve epochs
    whose camera spacing differs by ~1.5%, and it is far too tight for
    one disputed value to be mistaken for another (they differ by 240%).

    Returns None when fewer than 3 features are available, else
    {'n_match', 'n_stellar', 'n_lines' (candidates on offer), 'matches':
    [(line, x, resid_px)], 'coeffs' (polyval order in lambda - 6562.80,
    signed), 'disp_at_ha' (signed A/px), 'x_halpha', 'rms_px', 'sign',
    'refined' (True when >= 4 lines survived the tight re-match)}.
    Ranking: refined first, then most matches, most stellar matches,
    smallest rms.
    """
    a0, b0 = SEEDS[grism] if seed is None else seed
    a0, b0 = a0 * pixel_scale, b0 * pixel_scale
    cand = [LINE_BY_NAME[n] for n in GRISM_LINES[grism]
            if LINE_BY_NAME[n].spt == "*"
            or any(c in LINE_BY_NAME[n].spt for c in spt)]
    if len(features) < 3:
        return None
    xs = np.array([f["x"] for f in features])
    waves = np.array([_id_wave(ln) for ln in cand])
    stellar = np.array([ln.kind == "stellar" for ln in cand])
    i_ha = [i for i, ln in enumerate(cand) if ln.name == "Halpha"][0]
    anchors = sorted(features, key=lambda f: -f["depth"])[:ID_N_ANCHORS]
    best, best_key = None, None
    for anc in anchors:
        for sgn in ((1, -1) if sign is None else (sign,)):
            coeffs = [sgn * b0, sgn * a0, anc["x"]]
            mt = _count_matches(xs, waves, coeffs, anc["x"],
                                ID_TOL_PX, ID_TOL_FRAC)
            if not any(t[1] == i_ha for t in mt):
                continue
            refined = False
            for _ in range(2):          # free quadratic, tight re-match
                if len(mt) < 4:
                    break
                c = np.polyfit(waves[[t[1] for t in mt]] - HALPHA_A,
                               xs[[t[0] for t in mt]], 2)
                mt2 = _count_matches(xs, waves, c, anc["x"],
                                     ID_REFINE_TOL_PX, ID_REFINE_TOL_FRAC)
                if len(mt2) < 4:
                    break
                mt, coeffs, refined = mt2, list(c), True
            resid = np.array([t[2] for t in mt])
            rms = float(np.sqrt(np.mean(resid ** 2)))
            # A solution that survived the tight re-match outranks any
            # that did not: five lines within +/-2.5% of lever arm are
            # weaker evidence than four lines within 2.5 px.
            key = (int(refined), len(mt),
                   int(sum(stellar[t[1]] for t in mt)), -rms)
            if best_key is None or key > best_key:
                best_key = key
                best = {"n_match": len(mt), "n_lines": len(cand),
                        "n_stellar": key[2], "refined": refined,
                        "matches": [(cand[t[1]].name, float(xs[t[0]]), t[2])
                                    for t in mt],
                        "coeffs": [float(v) for v in coeffs],
                        "disp_at_ha": float(1.0 / coeffs[-2]),
                        "x_halpha": float(coeffs[-1]),
                        "rms_px": rms, "sign": sgn}
    return best


# --------------------------------------------------------------------------
# Precise feature measurement
# --------------------------------------------------------------------------
def fit_line(flux: np.ndarray, x_guess: float, halfwin: int,
             core_halfwin: Optional[int] = None,
             var: Optional[np.ndarray] = None) -> Optional[dict]:
    """Centre of one absorption or emission line: Gaussian + straight-line
    continuum, least squares, on RAW flux (not the normalised spectrum —
    the continuum slope is fitted, not assumed flat).

    Works for either sign (the amplitude is free).  Returns {'x',
    'x_err', 'sigma_px', 'fwhm_px', 'amp_frac' (amplitude / local
    continuum; negative = absorption), 'snr', 'rchi2'} or None when the
    window is unusable or the fit fails or lands outside the window.

    ``x_err`` is the formal least-squares error scaled by sqrt(rchi2)
    when rchi2 > 1 (a Gaussian is not the true shape of a rotationally
    broadened or Stark-broadened line; the scaled error absorbs that).
    """
    from scipy.optimize import least_squares
    lo = max(0, int(round(x_guess)) - halfwin)
    hi = min(len(flux), int(round(x_guess)) + halfwin + 1)
    x = np.arange(lo, hi, dtype=float)
    y = flux[lo:hi].astype(float)
    good = np.isfinite(y)
    if good.sum() < 12:
        return None
    x, y = x[good], y[good]
    if var is not None:
        e = np.sqrt(np.clip(var[lo:hi][good], 1e-12, None))
    else:
        # Empirical noise from the pixel-to-pixel differences (the line
        # and continuum are smooth on the 1-px scale).
        e = np.full(len(y), max(1.4826 * np.median(np.abs(np.diff(y)))
                                / np.sqrt(2.0), 1e-9))
    cont0 = np.median(np.concatenate([y[:5], y[-5:]]))
    chw = core_halfwin or max(3, halfwin // 3)
    core = np.abs(x - x_guess) <= chw
    if not core.any():
        return None
    dev = y - cont0
    j = int(np.argmax(np.abs(np.where(core, dev, 0.0))))
    p0 = [dev[j], x[j], max(2.0, chw / 2.0), cont0, 0.0]

    def model(p):
        return (p[0] * np.exp(-0.5 * ((x - p[1]) / p[2]) ** 2)
                + p[3] + p[4] * (x - x_guess))

    try:
        sol = least_squares(lambda p: (model(p) - y) / e, p0,
                            bounds=([-np.inf, lo, 0.7, -np.inf, -np.inf],
                                    [np.inf, hi, halfwin, np.inf, np.inf]))
    except Exception:                               # noqa: BLE001
        return None
    p = sol.x
    dof = max(1, len(y) - 5)
    rchi2 = float((sol.fun ** 2).sum() / dof)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac)
    except np.linalg.LinAlgError:
        return None
    x_err = float(np.sqrt(max(cov[1, 1], 0.0)) * np.sqrt(max(rchi2, 1.0)))
    amp_err = float(np.sqrt(max(cov[0, 0], 0.0)) * np.sqrt(max(rchi2, 1.0)))
    if not (lo + 2 <= p[1] <= hi - 3) or p[3] <= 0:
        return None
    return {"x": float(p[1]), "x_err": x_err, "sigma_px": float(p[2]),
            "fwhm_px": float(2.3548 * p[2]),
            "amp_frac": float(p[0] / p[3]),
            "snr": float(abs(p[0]) / max(amp_err, 1e-12)),
            "rchi2": rchi2}


def edge_position(flux: np.ndarray, x_head: float, red_sign: int,
                  search_px: int, span_px: float,
                  sigma_px: float = 1.5) -> Optional[dict]:
    """Half-depth position of a telluric band's BLUE edge.

    A band head is a cliff: continuum on the blue side, absorption
    beginning abruptly and shading to the red.  The measured quantity is
    where the (lightly smoothed) spectrum crosses HALF of the drop from
    the blue-side continuum to the band floor.

    * blue-side continuum: a straight line fitted over a window of
      ``2*span_px`` blueward of the search range and extrapolated across
      the band (the local continuum slope is real — the red-end response
      falls steeply — and a flat estimate would bias the crossing);
    * band floor: the minimum of flux/continuum within ``span_px``
      redward of the steepest descent;
    * the crossing is located by linear interpolation between pixels.

    ``red_sign`` is +1 when wavelength increases with x, -1 otherwise.
    ``search_px`` is how far either side of ``x_head`` the cliff may sit.
    Returns {'x', 'depth', 'x_err'} or None.  ``x_err`` is the noise-
    propagated crossing error: sigma(flux ratio) / |slope at crossing|.
    """
    n = len(flux)
    sm = smooth(flux.astype(float), sigma_px)
    # Work in a frame where red is toward +index.
    if red_sign < 0:
        sm = sm[::-1]
        raw = flux[::-1]
        x_h = (n - 1) - x_head
    else:
        raw = flux
        x_h = x_head
    span = int(round(span_px))
    b0 = int(round(x_h - search_px - 2 * span))
    b1 = int(round(x_h - search_px))
    r1 = int(round(x_h + search_px + span))
    if b0 < 0 or r1 >= n or b1 - b0 < 8:
        return None
    xb = np.arange(b0, b1)
    yb = sm[b0:b1]
    okb = np.isfinite(yb)
    if okb.sum() < 8:
        return None
    c = np.polyfit(xb[okb], yb[okb], 1)
    xs = np.arange(b1, r1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sm[b1:r1] / np.polyval(c, xs)
    if not np.isfinite(ratio).all():
        return None
    # Steepest descent inside the search range = the cliff.
    grad = np.gradient(ratio)
    in_search = np.abs(xs - x_h) <= search_px
    if not in_search.any():
        return None
    j_cliff = int(np.argmin(np.where(in_search, grad, np.inf)))
    floor_seg = ratio[j_cliff:min(len(ratio), j_cliff + span + 1)]
    floor = float(np.min(floor_seg))
    top = float(np.median(ratio[max(0, j_cliff - span):max(1, j_cliff - 2)]))
    depth = top - floor
    if depth <= 0:
        return None
    half = top - 0.5 * depth
    # Walk blueward from the floor to the first pixel above half-depth.
    j_floor = j_cliff + int(np.argmin(floor_seg))
    j = j_floor
    while j > 0 and ratio[j] < half:
        j -= 1
    if ratio[j] < half or j + 1 >= len(ratio):
        return None
    slope = ratio[j + 1] - ratio[j]
    if slope >= 0:
        return None
    x_cross = xs[j] + (half - ratio[j]) / slope
    # Noise on the ratio from the blue-side scatter of the RAW spectrum.
    rb = raw[b0:b1] / np.polyval(c, xb)
    noise = 1.4826 * np.nanmedian(np.abs(rb - np.nanmedian(rb)))
    # Smoothing reduces single-pixel noise by ~ (2 sqrt(pi) sigma)^-1/2.
    noise_sm = noise / np.sqrt(max(1.0, 2.0 * np.sqrt(np.pi) * sigma_px))
    x_err = float(noise_sm / abs(slope))
    if red_sign < 0:
        x_cross = (n - 1) - x_cross
    return {"x": float(x_cross), "depth": float(depth), "x_err": x_err}


# --------------------------------------------------------------------------
# The fixed-dispersion solution
# --------------------------------------------------------------------------
def solve_dispersion(frame_ids: Sequence, waves: Sequence[float],
                     xs: Sequence[float], x_errs: Sequence[float],
                     degree: int = 2, x_ref: float = 0.0) -> dict:
    """ONE wavelength polynomial in DETECTOR coordinates for many frames;
    one constant per frame.

    Model, per measurement i in frame f, with t = (x_i - x_ref)/1000:

        lambda_i = c_f + k1 t + k2 t^2 [+ k3 t^3]

    WHY DETECTOR COORDINATES (2026-10-03).  The first form of this fit
    was x(lambda) with a per-frame pixel shift, i.e. it assumed the
    non-linearity of the scale travels WITH the star.  The T CrB series
    falsified that: the pixel distance from Halpha to the O2-B edge is
    not constant but changes by ~5% when the star sits 500 px from its
    usual place, at very nearly the rate expected if the SAME
    non-linearity were instead fixed to the detector.  A scale that is
    stretched as a function of where on the chip the light lands is an
    optical distortion of the grism-camera combination, and it is the
    same for every star; so the polynomial is a function of x, and the
    only per-frame freedom is an additive constant in WAVELENGTH (which
    star position, and the target's radial velocity, set).

    The measurement errors are in x; they are converted to wavelength
    errors with the local slope and the fit is iterated once.  Weighted
    linear least squares; parameter errors are formal errors scaled by
    sqrt(chi2_nu) when chi2_nu > 1, and chi2_nu is REPORTED as it comes
    (standing rule 1: never clipped at 1).

    Returns {'coeffs': [k1, k2, (k3)] (A per kilo-pixel^j), 'coeff_errs',
    'x_ref', 'degree', 'disp_ref' (signed A/px at x_ref), 'disp_ref_err',
    'const': {frame: c_f}, 'resid_a', 'resid_px' (input order),
    'rms_px', 'rms_a', 'chi2', 'dof', 'rchi2', 'n_frames', 'n_meas'}.
    """
    frame_ids = list(frame_ids)
    waves = np.asarray(waves, dtype=float)
    xs = np.asarray(xs, dtype=float)
    x_errs = np.asarray(x_errs, dtype=float)
    frames = sorted(set(frame_ids), key=str)
    col = {f: k for k, f in enumerate(frames)}
    nf = len(frames)
    a_mat = np.zeros((len(xs), nf + degree))
    for i, f in enumerate(frame_ids):
        a_mat[i, col[f]] = 1.0
    t = (xs - x_ref) / 1000.0
    for k in range(1, degree + 1):
        a_mat[:, nf + k - 1] = t ** k
    slope = np.full(len(xs), 1.0)          # |A/px|, refined below
    for _ in range(2):
        w = 1.0 / np.clip(x_errs * slope, 1e-6, None)
        sol, *_ = np.linalg.lstsq(a_mat * w[:, None], waves * w, rcond=None)
        coeffs = [float(v) for v in sol[nf:]]
        slope = np.abs(local_dispersion(xs, coeffs, x_ref))
    resid = waves - a_mat @ sol
    dof = len(xs) - a_mat.shape[1]
    chi2 = float(((resid * w) ** 2).sum())
    rchi2 = chi2 / dof if dof > 0 else float("nan")
    cov = np.linalg.pinv((a_mat * w[:, None]).T @ (a_mat * w[:, None]))
    scale = np.sqrt(rchi2) if dof > 0 and rchi2 > 1 else 1.0
    errs = [float(np.sqrt(cov[nf + k, nf + k]) * scale)
            for k in range(degree)]
    resid_px = resid / slope
    return {"coeffs": coeffs, "coeff_errs": errs, "degree": degree,
            "x_ref": float(x_ref),
            "disp_ref": coeffs[0] / 1000.0, "disp_ref_err": errs[0] / 1000.0,
            "const": {f: float(sol[col[f]]) for f in frames},
            "resid_a": resid, "resid_px": resid_px,
            "rms_px": float(np.sqrt(np.mean(resid_px ** 2))),
            "rms_a": float(np.sqrt(np.mean(resid ** 2))),
            "chi2": chi2, "dof": int(dof), "rchi2": float(rchi2),
            "n_frames": nf, "n_meas": int(len(xs))}


def scale_poly(x, coeffs: Sequence[float], x_ref: float):
    """P(x) = sum_k coeffs[k-1] t^k, t = (x - x_ref)/1000: the wavelength
    scale up to an additive constant (Angstrom)."""
    t = (np.asarray(x, dtype=float) - x_ref) / 1000.0
    out = np.zeros_like(t)
    for k, c in enumerate(coeffs, start=1):
        out = out + c * t ** k
    return out


def local_dispersion(x, coeffs: Sequence[float], x_ref: float):
    """Signed A per px at detector position x (dP/dx)."""
    t = (np.asarray(x, dtype=float) - x_ref) / 1000.0
    out = np.zeros_like(t)
    for k, c in enumerate(coeffs, start=1):
        out = out + k * c * t ** (k - 1)
    return out / 1000.0


def wavelength_of(x, x_anchor: float, coeffs: Sequence[float],
                  x_ref: float, anchor_wave: float = HALPHA_A):
    """Wavelength at detector position x when the feature of wavelength
    ``anchor_wave`` (Halpha by default) sits at ``x_anchor``:
    lambda = anchor_wave + P(x) - P(x_anchor)."""
    return (anchor_wave + scale_poly(x, coeffs, x_ref)
            - scale_poly(x_anchor, coeffs, x_ref))


def predict_x(wave, x_anchor: float, coeffs: Sequence[float], x_ref: float,
              anchor_wave: float = HALPHA_A, iters: int = 8):
    """Detector position of a wavelength (inverse of
    :func:`wavelength_of`), by Newton iteration from the linear estimate
    (the scale is monotonic over the detector)."""
    wave = np.asarray(wave, dtype=float)
    d0 = float(local_dispersion(x_anchor, coeffs, x_ref))
    x = x_anchor + (wave - anchor_wave) / d0
    for _ in range(iters):
        f = wavelength_of(x, x_anchor, coeffs, x_ref, anchor_wave) - wave
        x = x - f / local_dispersion(x, coeffs, x_ref)
    return x
