"""The identity gate for slitless frames.

STATUS (2026-10-05): the v1 gate below (Gaia field verification behind a
header-pointing test) is SUPERSEDED by the pixel fingerprint gate at the
end of this module (findings OA.E2, TE.F2: the header must never be the
sole reason a frame is rejected).  The v1 functions are kept because the
August validation rows in ``g_extractions`` were produced by them and the
report compares the two.

v1 — Gaia DR3 field verification:

A grism frame has no WCS and no point sources, so "is this really the
claimed field?" must be answered from the dispersed content itself.  The
exploration round (documented in the G report) established three facts that
fix the design:

1. **Ghosts, not stars.**  The secondary peaks of the cross-dispersion
   profile sit at the same offsets from the main trace on a T CrB frame
   and on a frame pointed 168 deg away — they are instrumental ghosts of
   the bright star.  A naive "match all trace positions to Gaia" gate
   would match ghosts and lie.  Only the MAIN trace position is treated as
   astrometric evidence.
2. **The camera grid is known without the grism.**  Plate-solving Mode0
   *imaging* frames from the grism season gives the era's CD matrix
   (0.4508"/px, rotation ~0.3 deg) — so a Gaia star's pixel position on a
   grism frame is predictable from the header pointing alone, up to the
   telescope's meridian-flip parity (a German equatorial rotates the field
   180 deg across the pier; the pier side is not in the headers, so BOTH
   parities are tried and the better one kept).
3. **The detilt coordinate cancels the grism.**  In u = y - slope*(x-nx/2)
   the grism's along-dispersion deflection drops out (see trace.py), so
   the predicted u of the brightest Gaia star in the cone can be compared
   directly with the observed main-trace u.

The gate then asks two independent questions and requires BOTH yeses:

* **header-vs-target** — is the header pointing within
  ``POINTING_TOL_DEG`` of the claimed target's reference coordinates?
  (The S0 manifest's pointing-outlier law; the 21 known-bad T CrB frames
  fail here by 50-168 deg.)
* **content-vs-header** — does the frame's main trace sit within
  ``U_TOL_PX`` (best parity) of where Gaia says the brightest star of the
  header's field should fall?  This is the arm that checks the PIXELS:
  measured on the good T CrB series the best-parity residual stays within
  ~250 px, while the known-bad frames miss by 700-1200 px — their content
  does not match even their own headers (the mount kept imaging T CrB
  while the header recorded a bogus pointing).

Everything that decides is pure; the Vizier cone search is the one impure
function and it caches to disk so re-runs are offline and deterministic.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report quotes these).
# --------------------------------------------------------------------------

#: Gate code version, recorded into g_build_meta.
G_CODE_VERSION = "G v2.1 (2026-10-05)"

#: header-vs-target tolerance.  The manifest's own pointing-outlier law
#: (manifest.POINTING_OUTLIER_DEG): the good T CrB series scatters at
#: ~0.013 deg, the bad frames start at 50 deg — 1 deg splits the two
#: populations with orders of magnitude to spare on both sides.
POINTING_TOL_DEG = 1.0

#: content-vs-header tolerance on the main-trace u residual (px, best
#: parity).  Calibrated on the good T CrB series: best-parity residuals
#: stay within ~250 px (the header pointing itself is only good to a few
#: arcmin); the known-bad frames miss by >= 700 px.  400 px (~3 arcmin)
#: sits between the populations with margin on both sides.
U_TOL_PX = 400.0

#: Minimum main-trace height (halo-subtracted median-collapse ADU) for the
#: frame to count as containing a spectrum at all.  Good 240 s T CrB
#: frames measure 40-260; one cloud-dead frame measured 4.  A frame
#: without a believable trace cannot pass any identity check.
MIN_TRACE_HEIGHT_ADU = 15.0

#: Gaia cone radius.  Half the frame diagonal is 0.36 deg at 0.4508"/px
#: (4788 x 3194 bin2); 0.45 deg covers every star whose trace can reach
#: the collapse window even after a worst-case in-tolerance pointing error.
CONE_RADIUS_DEG = 0.45

#: Gaia G limit for the cone.  The gate only ever uses the BRIGHTEST star
#: of the field; 14 keeps the cone small while guaranteeing the brightest
#: star of any field the archive points at is inside it.
CONE_G_LIMIT = 14.0

#: Pixel margin inside which a predicted star counts as on-frame.  Slack
#: (one U_TOL) is added around the physical frame so an in-tolerance
#: pointing error cannot push the true brightest star out of the model.
ONFRAME_MARGIN_PX = 400.0

#: Mode0 bin2 plate scale (arcsec/px) — S2 registry fact, used only for
#: reporting; the pixel mapping itself uses the solved CD matrix.
PLATE_SCALE_ARCSEC = 0.4508


@dataclass(frozen=True)
class GateResult:
    """One frame's identity-gate outcome, ready for the DB row."""
    verdict: str                 # ACCEPT / REJECT
    reason: str                  # 'ok' or the first failed check
    pointing_offset_deg: Optional[float]
    trace_height: Optional[float]
    u_obs: Optional[float]
    u_pred_best: Optional[float]
    u_resid_px: Optional[float]  # obs - pred, best parity
    parity: Optional[str]        # 'A' (as solved) or 'B' (meridian flip)
    n_gaia: int
    brightest_g: Optional[float]


# --------------------------------------------------------------------------
# Geometry: sky -> pixel -> detilt coordinate, both parities
# --------------------------------------------------------------------------

def sky_to_pixel_offsets(cd_matrix: np.ndarray, ra0: float, dec0: float,
                         ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """Pixel offsets (dx, dy) from the frame center for sky positions,
    using the era's solved CD matrix (deg/px).

    Small-field gnomonic approximation: the intermediate world coordinates
    are (Δα·cosδ, Δδ) in degrees — over a 0.5 deg field the TAN projection
    correction is < 0.1 px, far below the gate tolerance.  Returns an
    (N, 2) array.
    """
    cdi = np.linalg.inv(cd_matrix)
    dra = (np.asarray(ra) - ra0) * np.cos(np.radians(np.asarray(dec)))
    ddec = np.asarray(dec) - dec0
    return (cdi @ np.vstack([dra, ddec])).T


def predicted_u(cd_matrix: np.ndarray, parity: str, ra0: float,
                dec0: float, ra: np.ndarray, dec: np.ndarray,
                slope: float, ny: int, nx: int) -> np.ndarray:
    """Predicted detilt coordinate u for stars, one parity.

    Parity 'A' uses the CD matrix as solved from the calibration imaging
    frames; parity 'B' is the meridian flip (both axes negated).  The
    star's trace peaks at u = y_center + dy - slope*dx — the same detilt
    formula the observed profile uses, so the two sides are comparable
    (see trace.detilted_profile for why the grism deflection cancels).
    """
    xy = sky_to_pixel_offsets(cd_matrix, ra0, dec0, ra, dec)
    if parity == "B":
        xy = -xy
    return ny / 2.0 + xy[:, 1] - slope * xy[:, 0]


def on_frame(cd_matrix: np.ndarray, parity: str, ra0: float, dec0: float,
             ra: np.ndarray, dec: np.ndarray, ny: int, nx: int,
             margin: float = ONFRAME_MARGIN_PX) -> np.ndarray:
    """Boolean mask: which stars land on (or within ``margin`` of) the
    frame for this parity — the candidate set for 'brightest star'."""
    xy = sky_to_pixel_offsets(cd_matrix, ra0, dec0, ra, dec)
    if parity == "B":
        xy = -xy
    return ((np.abs(xy[:, 0]) < nx / 2 + margin)
            & (np.abs(xy[:, 1]) < ny / 2 + margin))


def brightest_prediction(cd_matrix: np.ndarray, ra0: float, dec0: float,
                         stars: np.ndarray, slope: float,
                         ny: int, nx: int):
    """For each parity: (u of the brightest on-frame Gaia star, its G).

    ``stars`` is an (N, 3) array of (ra, dec, Gmag).  Returns a dict
    {'A': (u, g) | None, 'B': ...} — None when no star lands on frame for
    that parity (possible at the survey's sparse high-latitude pointings).
    """
    out = {}
    for parity in ("A", "B"):
        mask = on_frame(cd_matrix, parity, ra0, dec0,
                        stars[:, 0], stars[:, 1], ny, nx)
        if not mask.any():
            out[parity] = None
            continue
        sub = stars[mask]
        i0 = int(np.argmin(sub[:, 2]))         # smallest G = brightest
        u = predicted_u(cd_matrix, parity, ra0, dec0,
                        sub[i0:i0 + 1, 0], sub[i0:i0 + 1, 1],
                        slope, ny, nx)[0]
        out[parity] = (float(u), float(sub[i0, 2]))
    return out


# --------------------------------------------------------------------------
# The verdict (pure — the full truth table is unit-tested)
# --------------------------------------------------------------------------

def gate_verdict(pointing_offset_deg: Optional[float],
                 trace_height: Optional[float],
                 u_obs: Optional[float],
                 predictions: dict,
                 n_gaia: int,
                 pointing_tol: float = POINTING_TOL_DEG,
                 u_tol: float = U_TOL_PX,
                 min_height: float = MIN_TRACE_HEIGHT_ADU) -> GateResult:
    """Combine the checks into one verdict.  Checks run in evidence order
    and the FIRST failure names the reason; every number that was computed
    is kept on the result even when a later check never ran (the DB row is
    the forensic record).

    Fail-closed policy: a frame that cannot show its evidence (no pointing
    in the header, no Gaia stars, no trace) is REJECTED, never waved
    through.
    """
    # Best-parity prediction, computed up front so it lands in the record
    # regardless of which check fails.
    best = None                  # (resid, u_pred, parity)
    if u_obs is not None:
        for parity in ("A", "B"):
            p = predictions.get(parity)
            if p is None:
                continue
            resid = u_obs - p[0]
            if best is None or abs(resid) < abs(best[0]):
                best = (resid, p[0], parity, p[1])
    resid, u_pred, parity, gmag = best if best else (None, None, None, None)

    def result(verdict, reason):
        return GateResult(verdict=verdict, reason=reason,
                          pointing_offset_deg=pointing_offset_deg,
                          trace_height=trace_height, u_obs=u_obs,
                          u_pred_best=u_pred, u_resid_px=resid,
                          parity=parity, n_gaia=n_gaia, brightest_g=gmag)

    if pointing_offset_deg is None:
        return result("REJECT", "no_header_pointing")
    if pointing_offset_deg > pointing_tol:
        return result("REJECT", "header_off_target")
    if trace_height is None or trace_height < min_height:
        return result("REJECT", "no_trace")
    if n_gaia == 0:
        return result("REJECT", "no_gaia_catalog")
    if best is None:
        return result("REJECT", "no_onframe_star")
    if abs(resid) > u_tol:
        return result("REJECT", "field_mismatch")
    return result("ACCEPT", "ok")


# --------------------------------------------------------------------------
# The impure edge: cached Gaia DR3 cone search (Vizier)
# --------------------------------------------------------------------------

def cone_cache_key(ra: float, dec: float, radius: float,
                   g_limit: float) -> str:
    """Deterministic cache filename for one cone.  Coordinates are rounded
    to 0.01 deg (~36") — pointings inside the same 36" share a cone, which
    is exactly right: the cone radius carries 400 px of slack."""
    return (f"gaia_{ra:.2f}_{dec:+.2f}_r{radius:.2f}"
            f"_g{g_limit:.1f}.json").replace("+", "p").replace("-", "m")


#: Empty-cone retries.  The first validation run PROVED Vizier sometimes
#: returns an empty table list transiently (the T CrB cone at 239.88+25.91
#: came back empty while its 36-arcsec neighbour cone held 83 stars) —
#: and a cached empty then poisons every later run as 'no_gaia_catalog'.
#: An empty answer is therefore only believed (and only cached) after it
#: repeats EMPTY_RETRIES times; a non-empty answer is believed at once.
EMPTY_RETRIES = 3
EMPTY_RETRY_SLEEP_S = 5.0


#: Vizier servers, tried in order.  During the validation run the
#: Strasbourg host went through a phase of answering every cone with an
#: EMPTY table list (no error!) while the CfA mirror answered correctly —
#: so a mirror is not an optimization, it is required for correctness.
VIZIER_SERVERS = ("vizier.cds.unistra.fr", "vizier.cfa.harvard.edu")


def _vizier_cone_query(ra: float, dec: float, radius: float,
                       g_limit: float) -> list:
    """The one real network call: one Vizier cone, as [[ra, dec, G], ...].

    Vizier (catalog I/355/gaiadr3) rather than the ESA TAP endpoint: the
    exploration round measured ~7 s per Vizier cone against a TAP call
    that ran past two minutes.  Servers are tried in ``VIZIER_SERVERS``
    order; the first NON-EMPTY answer wins (an empty answer from one
    server is checked against the next before it is believed); an
    exception moves to the next server, and only the last server's
    exception propagates.
    """
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    last_exc, answered = None, False
    for server in VIZIER_SERVERS:
        try:
            v = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag"],
                       column_filters={"Gmag": f"<{g_limit}"},
                       row_limit=50000)
            v.VIZIER_SERVER = server
            tables = v.query_region(
                SkyCoord(ra * u.deg, dec * u.deg),
                radius=radius * u.deg, catalog="I/355/gaiadr3")
        except Exception as exc:                 # noqa: BLE001
            last_exc = exc
            continue
        answered = True
        rows = []
        if tables:                               # empty list = no rows
            t = tables[0]
            rows = [[float(a), float(b), float(g)] for a, b, g in
                    zip(t["RA_ICRS"], t["DE_ICRS"], t["Gmag"])
                    if np.isfinite(g)]
        if rows:
            return rows
    if not answered and last_exc is not None:
        raise last_exc                           # no server ever answered
    return []                                    # answered, and empty


def gaia_cone(ra: float, dec: float, cache_dir: str,
              radius: float = CONE_RADIUS_DEG,
              g_limit: float = CONE_G_LIMIT,
              _query=None, _sleep=None) -> np.ndarray:
    """Gaia DR3 stars (ra, dec, Gmag) within ``radius`` of a pointing,
    as an (N, 3) array — served from the JSON cache when present.

    A NON-EMPTY cache entry is always trusted.  An EMPTY cache entry is
    NOT trusted (see ``EMPTY_RETRIES``): it is re-queried, and replaced
    the moment the query returns stars — so the transient-failure empties
    of an earlier run heal themselves on the next pass.  A fresh empty
    answer is cached only after ``EMPTY_RETRIES`` consecutive empty
    replies: a pointing with no G<14 star IS a fact worth caching, but
    only once it has been asked firmly enough to rule out a hiccup.

    ``_query`` / ``_sleep`` are test injection points for the network
    call and the retry pause (the tests exercise the retry ladder without
    astroquery or wall-clock).
    """
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    cache = Path(cache_dir) / cone_cache_key(ra, dec, radius, g_limit)
    if cache.exists():
        cached = json.loads(cache.read_text()) or []
        if cached:                               # non-empty: trusted
            return np.array(cached).reshape(-1, 3)
    query = _query or _vizier_cone_query
    sleep = _sleep or (lambda s: __import__("time").sleep(s))
    rows: list = []
    for attempt in range(EMPTY_RETRIES):
        try:
            rows = query(ra, dec, radius, g_limit)
        except Exception:                        # noqa: BLE001
            # A raised query (Vizier read timeout — observed in the
            # validation run) gets the same ladder as an empty reply;
            # only the LAST attempt's exception propagates to the
            # caller, which records it as a gaia_error row.
            if attempt == EMPTY_RETRIES - 1:
                raise
            sleep(EMPTY_RETRY_SLEEP_S)
            continue
        if rows:
            break                                # stars found: believe it
        if attempt < EMPTY_RETRIES - 1:
            sleep(EMPTY_RETRY_SLEEP_S)           # hiccup? ask again
    # Atomic-ish cache write: temp then rename, same directory.
    tmp = cache.with_suffix(".tmp")
    tmp.write_text(json.dumps(rows))
    os.replace(tmp, cache)
    return np.array(rows or []).reshape(-1, 3)


def angular_offset_deg(ra1, dec1, ra2, dec2) -> float:
    """Great-circle separation in degrees (haversine — exact at any
    separation; the bad pointings reach 168 deg where small-angle math
    would lie)."""
    p1, p2 = math.radians(dec1), math.radians(dec2)
    dl = math.radians(ra2 - ra1)
    a = (math.sin((p2 - p1) / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2)
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(a))))


# ==========================================================================
# v2 (2026-10-05): the PIXEL identity gate — a spectral fingerprint
# ==========================================================================
# WHY THE v1 GATE IS RETIRED (findings OA.E2, TE.F2, RF)
# -----------------------------------------------------
# v1 asked the header first: a frame whose RA/Dec card was > 1 deg from
# the target was rejected before a pixel was looked at.  The 21 T CrB
# frames it rejected that way carry a stale pointing card (RA 05h29, a
# position below the horizon at the time) while the pixels show T CrB's
# spectrum.  Its pixel arm compared one coordinate of one star under two
# meridian-flip parities the mount never used.  A header can be wrong in
# either direction; a spectrum cannot pretend to be T CrB's.
#
# WHAT T CrB LOOKS LIKE THROUGH THE GRISM
# ---------------------------------------
# A symbiotic: an M4 III continuum (TiO bands — the 7054 A head is a
# 20-40% step on the hrg) with strong, narrow Halpha emission.  Few stars
# on the archive's grism target lists show both, and none shows them at
# T CrB's relative strengths.  The gate therefore asks, from the pixels
# alone:
#
#   1. is there a trace at all, and a sharp emission line to anchor on?
#   2. with that line as Halpha and the FIXED dispersion of the frame's
#      (grism, mechanical epoch), does the response-free 6300-7400 A
#      spectrum correlate with T CrB's template (median of the OTHER
#      nights' T CrB frames)?
#   3. is the TiO 7054 step present?
#
# Thresholds come from T CrB's own frames only (see FP_COMPLETENESS); the
# false-accept rate on non-T CrB frames of the same grisms, camera state
# and exposure regime is then MEASURED and published.  The header pointing
# is recorded beside every verdict and never consulted.

#: Fingerprint window (A) and the wavelength step of the common grid.
FP_WAVE = (6300.0, 7400.0)
FP_STEP_A = {"hrg": 1.0, "lrg": 4.0}

#: Halpha itself is excluded from the correlation (+/- this many A): the
#: anchor line would otherwise correlate with itself on ANY emission-line
#: star, and the test is about the continuum's fingerprint.
FP_HALPHA_GUARD_A = 30.0

#: TiO 7054 band: flux just redward of the head over flux just blueward.
TIO_BLUE_A = (6990.0, 7045.0)
TIO_RED_A = (7060.0, 7120.0)

#: Pre-registered thresholds (2026-10-05, before scoring the controls).
#: THEY FAILED and are kept as the record: on the raw (response-shaped)
#: spectra every hot star's continuum correlates with T CrB's at r ~ 0.9,
#: because the grism's response dominates the shape, and the response
#: alone lowers the red/blue ratio across 7054 A to ~0.8 on the hrg; 51 of
#: 160 non-T CrB frames (all hrg Be stars) passed.  The report prints the
#: pre-registered verdicts beside the adopted ones.
FP_MIN_R_PREREG = 0.80
TIO_MAX_STEP_PREREG = 0.90

#: The adopted gate (v2.1) removes the response before correlating: each
#: spectrum is divided by its own running median over FP_HIGHPASS_A, so
#: only structure narrower than the response — band heads and lines —
#: enters the correlation.  Its two thresholds are DERIVED FROM T CrB
#: ALONE (the value that keeps FP_COMPLETENESS of the target's frames,
#: two-fold cross-validated by night parity), so the control frames never
#: influence them and the false-accept rate measured on the controls is
#: an out-of-sample number.
FP_HIGHPASS_A = 120.0
FP_COMPLETENESS = 0.95

GATE_ACCEPT, GATE_REJECT = "ACCEPT", "REJECT"


def resample(wave: np.ndarray, flux: np.ndarray, grid: np.ndarray):
    """Flux on a wavelength grid by linear interpolation (NaN outside the
    finite, sorted part of the input)."""
    ok = np.isfinite(wave) & np.isfinite(flux)
    if ok.sum() < 10:
        return np.full(len(grid), np.nan)
    w, f = wave[ok], flux[ok]
    order = np.argsort(w)
    w, f = w[order], f[order]
    out = np.interp(grid, w, f, left=np.nan, right=np.nan)
    return out


def normalise(f: np.ndarray) -> np.ndarray:
    """Divide by the median (a fingerprint is about SHAPE, not level)."""
    med = np.nanmedian(f)
    return f / med if np.isfinite(med) and med > 0 else f * np.nan


def tio_step(grid: np.ndarray, f: np.ndarray) -> Optional[float]:
    """Median flux redward of the TiO 7054 head over blueward."""
    b = (grid >= TIO_BLUE_A[0]) & (grid <= TIO_BLUE_A[1]) & np.isfinite(f)
    r = (grid >= TIO_RED_A[0]) & (grid <= TIO_RED_A[1]) & np.isfinite(f)
    if b.sum() < 5 or r.sum() < 5:
        return None
    blue = float(np.median(f[b]))
    return float(np.median(f[r]) / blue) if blue > 0 else None


def highpass(grid: np.ndarray, f: np.ndarray,
             win_a: float = FP_HIGHPASS_A) -> np.ndarray:
    """f / running-median(f) - 1 over a window of ``win_a`` Angstrom:
    the spectrum with the instrument response (and the star's own
    broad continuum) divided out."""
    from scipy.ndimage import median_filter
    step = float(np.median(np.diff(grid)))
    win = max(3, int(round(win_a / step)) | 1)
    ok = np.isfinite(f)
    if ok.sum() < win:
        return np.full(len(f), np.nan)
    filled = np.interp(np.arange(len(f)), np.where(ok)[0], f[ok])
    cont = median_filter(filled, size=win, mode="nearest")
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(ok & (cont > 0), filled / cont - 1.0, np.nan)
    return out


def derive_thresholds(r_values, step_values,
                      completeness: float = FP_COMPLETENESS):
    """(r_min, step_max) that keep ``completeness`` of the TARGET's own
    frames: the (1 - completeness) quantile of r and the completeness
    quantile of the TiO step.  Nothing but target frames goes in."""
    r = np.asarray([v for v in r_values if v is not None], dtype=float)
    st = np.asarray([v for v in step_values if v is not None], dtype=float)
    return (float(np.quantile(r, 1.0 - completeness)),
            float(np.quantile(st, completeness)))


def fingerprint_r(grid: np.ndarray, f: np.ndarray,
                  template: np.ndarray) -> Optional[float]:
    """Pearson correlation of a normalised spectrum with the template on
    the common grid, Halpha excluded; None with < 50% overlap."""
    use = (np.isfinite(f) & np.isfinite(template)
           & (np.abs(grid - 6562.8) > FP_HALPHA_GUARD_A))
    if use.sum() < 0.5 * len(grid):
        return None
    a, b = f[use], template[use]
    if np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def fingerprint_verdict(has_trace: bool, halpha_snr: Optional[float],
                        r: Optional[float], step: Optional[float],
                        r_min: float = FP_MIN_R_PREREG,
                        step_max: float = TIO_MAX_STEP_PREREG
                        ) -> tuple[str, str]:
    """The pixel gate's verdict and its reason (first failed check).
    Nothing here takes a header value."""
    if not has_trace:
        return GATE_REJECT, "no_trace"
    if halpha_snr is None:
        return GATE_REJECT, "no_emission_line"
    if r is None:
        return GATE_REJECT, "no_wavelength_overlap"
    if r < r_min:
        return GATE_REJECT, "fingerprint_mismatch"
    if step is None:
        return GATE_REJECT, "tio_window_off_spectrum"
    if step > step_max:
        return GATE_REJECT, "no_tio_step"
    return GATE_ACCEPT, "ok"
