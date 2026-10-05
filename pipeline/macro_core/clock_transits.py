"""Pure S3b logic: an ABSOLUTE clock check from archived transits and eclipses.

WHY THIS MODULE EXISTS
----------------------
S3 (``macro_core.timing``) proves what instant a header stamp refers to and
converts it to BJD_TDB to a millisecond.  What S3 could not do is say whether
the acquisition PC's clock was RIGHT: its one external check (AG LMi against
a survey ephemeris) carries a +-4,100 s ephemeris envelope, and the committee
review of 2026-10-03 (OA.E6, DS.F10, RF.M4) pointed out that a shared
~1,065 s offset in two polar ephemerides is exactly the size a clock error
would have to be.  The archive already contains the cure: students observed
dozens of exoplanet transits and post-common-envelope eclipsing binaries
whose mid-times are predicted to SECONDS by the literature.  Timing those
events with our own stamps measures (observatory clock) - (true time)
directly, per camera era.

Sign convention used everywhere: ``O - C = T_observed - T_predicted`` in
seconds, with T_observed read off OUR time axis.  A clock that runs AHEAD of
true time by D seconds stamps every event D seconds late, so it produces
``O - C = +D``.

WHAT IS HERE (all pure, unit-tested in ``pipeline/tests/test_clock_transits.py``)
----------------------------------------------------------------------------------
* a limb-darkened occultation model that is exact for ANY radius ratio —
  planet in front of star (k ~ 0.1), M dwarf in front of sdB (k ~ 1), M dwarf
  in front of white dwarf (k ~ 10) — by one-dimensional quadrature over the
  occulted disc.  No external transit package is needed, and the same
  symmetric template serves every clock target;
* the light-curve assembly rules (comparison-star pruning, the native-pixel
  saturation veto demanded by standing rule 4, empirical point errors);
* the mid-time fit: a blind grid over the whole observing window followed by
  a non-linear refinement, so a 1,000 s clock error would be FOUND, not
  assumed away by starting at the prediction;
* its error: a moving-block residual bootstrap (preserves red noise) and a
  prayer-bead cross-check;
* ephemeris propagation with its uncertainty, and the per-era combination
  (formal and scatter-based errors side by side; no chi-square rescaling);
* the stated bounds on rolling-shutter and mechanical-shutter timing terms
  (DE.F7) and the clock-era map.

The I/O — manifest queries, pixels, the products database, figures and the
report page — lives in ``pipeline/scripts/build_s3b_clock_transits.py``.
"""

from __future__ import annotations

import functools
import math
from typing import Optional, Sequence

import numpy as np

#: Recorded in ``s3b_build_meta``; bump when database content would change.
S3B_CODE_VERSION = "S3b v1.0 (2026-10-03)"

SECONDS_PER_DAY = 86400.0

# --------------------------------------------------------------------------
# Policy constants (single source of truth; the report interpolates these)
# --------------------------------------------------------------------------

#: Acceptance criterion of foundation task F-8 (SYNTHESIS section 3):
#: |O - C| < 120 s per era, or the offset is measured and carried.
ACCEPT_OC_S = 120.0

#: The shared CV edge offset whose origin this stage has to settle
#: (ST LMi 1,071 +- 97 s, EU UMa 1,060 +- 84 s; OA.E6).  Used ONLY as an
#: injected test size — never as a prior on any fit.
CV_OFFSET_S = 1065.0

#: A literature ephemeris may serve as a clock only if its propagated
#: 1-sigma prediction error at our epoch is at most this (half the
#: acceptance criterion): a looser one cannot certify 120 s.
EPH_SIGMA_MAX_S = 60.0

#: Census admission: out-of-event data required before first contact and
#: after last contact, in hours.  A one-sided event is refused for the same
#: reason S3's AG LMi gate refuses it — a symmetric template fitted to one
#: flank returns the flank's slope, not the midpoint.
MIN_MARGIN_H = 0.15
#: The same margin for a binary eclipse.  An HW Vir or white-dwarf eclipse
#: is ten to thirty times deeper than a transit and its ingress and egress
#: are sharp, so its midpoint is pinned by the two flanks themselves and
#: needs only enough out-of-eclipse data to anchor the baseline; observers
#: also schedule ~35-minute visits around a ~24-minute eclipse.  Three
#: minutes each side (plus the post-fit MIN_SIDE_POINTS gate) is the rule.
MIN_MARGIN_ECLIPSE_H = 0.05
#: THE SELECTION MUST NOT ASSUME THE ANSWER.  The margins above are
#: computed from the PREDICTED event time on OUR time axis; if our clock
#: were wrong by ~1,000 s, an event that the prediction puts half outside
#: the run could in truth be fully inside it (and vice versa).  A census
#: that admits only events predicted to be fully covered would therefore
#: silently discard exactly the series that could reveal a large clock
#: error.  So the census admits every event whose predicted margins fall
#: short by no more than this pad (twice the CV-sized offset), and the
#: decision "is the event two-sided?" is taken AFTER the blind fit, on the
#: fitted mid-time (MIN_SIDE_POINTS).
CLOCK_SEARCH_PAD_H = 2.0 * CV_OFFSET_S / 3600.0

#: Census admission: minimum frames and minimum run length of a series.
MIN_SERIES_FRAMES = 20
MIN_SERIES_SPAN_H = 0.5

#: After the fit: minimum number of points before first and after last
#: contact (same number S3 uses for its coverage gate).
MIN_SIDE_POINTS = 3

#: Aperture radii, in units of the reference frame's FWHM.  All three are
#: measured; the one with the smallest successive-difference scatter is
#: adopted BEFORE any fit, and the other two become robustness variants.
APERTURE_FWHM_FACTORS = (1.5, 2.25, 3.0)
#: Sky annulus (inner, outer) in units of the reference FWHM.
ANNULUS_FWHM_FACTORS = (4.0, 6.0)
#: Smallest aperture radius in pixels (pixelation floor).
APERTURE_MIN_PX = 3.0

#: Comparison stars: how many candidates are measured, how few may remain.
N_COMP_CANDIDATES = 10
MIN_COMPS = 2
#: A comparison star is dropped when its own differential scatter exceeds
#: this multiple of the median comparison scatter (iterated).
COMP_RMS_FACTOR = 2.5
#: A star must be measured (and unsaturated) on this fraction of frames.
COMP_MIN_FRAME_FRAC = 0.9

#: Frames whose ensemble flux falls below this fraction of the series
#: median are cloud/tracking casualties and are dropped.
ENSEMBLE_MIN_FRACTION = 0.5

#: A re-pointing: the target's centroid moves by more than this many
#: pixels between consecutive frames.  Unguided drift here is ~0.1 px per
#: frame, so a jump of several pixels is the mount being re-centred, and
#: every re-centring puts the star on different pixels (different flat-
#: field error, different hot pixels) — a step in the differential light
#: curve whose time is known from the centroids, independently of the
#: photometry.
JUMP_PX = 5.0
#: A segment between two jumps shorter than this many points is not given
#: its own offset (it could not constrain one).
JUMP_MIN_SEGMENT = 8
#: At most this many offsets per fit (the largest jumps are kept).
JUMP_MAX_STEPS = 6

#: Outlier rejection after the first fit: |residual| > this many robust
#: sigma (1.4826 MAD).  One pass, then one refit.
CLIP_SIGMA = 5.0

#: Residual bootstrap: number of resamples and block length in minutes.
#: The block must be long enough to carry the correlated (seeing, guiding,
#: flat-field) noise of ground-based differential photometry from one
#: resample to the next; ten minutes is the ingress time-scale of the
#: transits timed here, which is the time-scale that moves a mid-time.
N_BOOTSTRAP = 300
BOOT_BLOCK_MIN = 10.0
#: Deterministic seed base (each event adds its own id).
BOOT_SEED = 20261003

#: Injection-recovery: the signed time shifts (s) re-injected into each
#: event's own residuals, and the number of noise draws per shift.
INJECT_SHIFTS_S = (-CV_OFFSET_S, -300.0, -ACCEPT_OC_S, 0.0,
                   ACCEPT_OC_S, 300.0, CV_OFFSET_S)
N_INJECT_DRAWS = 8

#: Gauss-Legendre nodes of the occultation integral (see
#: :func:`occultation_depth`; accuracy is pinned by a unit test).
N_QUAD = 48

#: Exposure-time smearing: the model is averaged over this many sub-samples
#: when an exposure is longer than SMEAR_MIN_EXPTIME_S.
SMEAR_MIN_EXPTIME_S = 20.0
N_SMEAR = 5

#: Event status strings written to ``s3b_events.status``.
STATUS_OK = "ok"
STATUS_ONE_SIDED = "one_sided_coverage"
STATUS_TOO_FEW = "too_few_points"
STATUS_NO_DIP = "no_dip_found"
STATUS_SATURATED = "target_saturated"
STATUS_NO_COMPS = "no_comparison_stars"
STATUS_SHALLOW = "depth_inconsistent"
#: A one-sided eclipse re-fitted with its duration and depth FIXED to the
#: median of the same star's complete eclipses (see the build script's
#: ``anchored_refit``).  Supplementary: reported, never in an era mean.
STATUS_ANCHORED = "ok_anchored_shape"
#: Minimum number of complete, timed eclipses of a star before its shape
#: may anchor a one-sided one.
ANCHOR_MIN_EVENTS = 3

#: Precision gate for era means.  An event whose own measurement error
#: (statistical and model, in quadrature) exceeds the acceptance criterion
#: cannot say anything about that criterion; averaging several such
#: events mostly averages their unmodelled systematics.  They are
#: recorded and plotted as "low precision" and kept out of the means.
#: The gate looks at the ERROR only, never at the O - C value.
TIMING_GRADE_MAX_SIGMA_S = ACCEPT_OC_S
GRADE_TIMING = "timing"
GRADE_LOW = "low_precision"
GRADE_TWIN = "twin"
GRADE_ANCHORED = "anchored"

#: A fitted event must reproduce at least this fraction (and at most the
#: inverse) of the catalogue depth to be accepted as THE event rather than
#: a noise dip or a wrong star.
DEPTH_MIN_FRACTION = 0.4
DEPTH_MAX_FRACTION = 2.5

#: Minimum points in a fit.
MIN_FIT_POINTS = 12


# --------------------------------------------------------------------------
# Stated timing bounds that no event in the archive can measure (DE.F7)
# --------------------------------------------------------------------------

#: Each row: (term, applies to, bound in seconds, sign, basis).  These are
#: STATED BOUNDS from vendor figures, not measurements; the measured O - C
#: of each era bounds their sum empirically at a much coarser level, and
#: the report prints both.  ``sign`` says which way the true mid-exposure
#: lies relative to ``start + EXPTIME/2``.
TIMING_BOUNDS = (
    ("rolling shutter (row-sequential start of exposure)",
     "SBIG Aluma AC4040 (GSENSE4040), eras 1-44", 0.75, "late, 0..bound",
     "vendor full-frame download time 0.75 s is an upper limit on the "
     "top-to-bottom row skew; the target row lies somewhere inside it"),
    ("rolling shutter (row-sequential start of exposure)",
     "ZWO ASI6200MM (IMX455), eras 75-77, 84", 0.32, "late, 0..bound",
     "vendor full-resolution frame rate 3.19 fps -> one frame scan "
     "<= 1/3.19 s"),
    ("rolling shutter (row-sequential start of exposure)",
     "QHY600 (IMX455), eras 78-82", 0.40, "late, 0..bound",
     "vendor 16-bit full-frame rate 2.5 fps -> one frame scan <= 0.40 s"),
    ("mechanical shutter travel (open and close)",
     "Andor iKon-L 936, eras 45-74", 0.10, "either, |x| <= bound",
     "ASSUMED: the vendor documents a 45 mm iris shutter but the open/"
     "close time was not found in the documentation consulted; 0.1 s is "
     "a generous class figure for such shutters and is flagged for "
     "confirmation.  An iris shutter also exposes the field centre "
     "longer than the corners, which shifts no mid-time"),
    ("StackPro sub-read dead time (already in S3)",
     "AC4040 'High Gain StackPro'", 5.675, "late, 0..bound",
     "macro_core.timing.STACKPRO_MID_WORST_CASE_S (cadence bound); this "
     "stage MEASURES it wherever one transit was observed in StackPro "
     "and plain High Gain on the same night"),
)


# --------------------------------------------------------------------------
# Clock eras
# --------------------------------------------------------------------------

#: A "clock era" is a stretch of nights that shares one camera AND one
#: acquisition-software build string (SWCREATE).  The software string is
#: the best available proxy for "the same acquisition PC and clock
#: discipline"; the camera decides the readout timing terms.  Boundaries
#: are (first night, last night) inclusive, in the manifest's ``night``
#: convention.  The label is what the report prints.
CLOCK_ERAS = (
    ("A", "AC4040, 2023 season (MaxIm 6.30)", "2023-02-01", "2023-07-31"),
    ("B", "AC4040, 2023-24 season (MaxIm 6.30/6.40)",
     "2023-10-01", "2024-03-31"),
    ("C", "Andor iKon-L, 2024 (MaxIm 6.40)", "2024-04-01", "2024-12-12"),
    ("D", "ASI6200, 2024-12 to 2025-06 (MaxIm 6.40)",
     "2024-12-13", "2025-08-31"),
    ("E", "ASI6200, 2025-10 to 2026-03 (MaxIm 6.30)",
     "2025-09-01", "2026-03-20"),
    ("F", "QHY600, 2026-03 to 2026-06 (MaxIm 6.40)",
     "2026-03-21", "2026-06-27"),
    ("G", "QHY600, pyscope (2026-06-28 on)", "2026-06-28", "2099-12-31"),
)

#: The iKon and the ASI overlapped for ten nights in December 2024; the
#: manifest era id (not the night) says which camera took a frame.
IKON_ERA_IDS = range(45, 75)
ASI_ERA_IDS = (75, 76, 77, 84)


def clock_era(night: str, era_id: Optional[int] = None) -> str:
    """Clock-era code (``'A'``..``'G'``) of one night, or ``'?'``.

    ``era_id`` breaks the December-2024 overlap, when the iKon and the ASI
    were both taking frames: an iKon era id always maps to ``'C'`` and an
    ASI era id observed before September 2025 to ``'D'``.
    """
    if era_id is not None:
        if int(era_id) in IKON_ERA_IDS:
            return "C"
        if int(era_id) in ASI_ERA_IDS and night < "2025-09-01":
            return "D"
    for code, _label, first, last in CLOCK_ERAS:
        if first <= night <= last:
            return code
    return "?"


def clock_era_label(code: str) -> str:
    """Human label of a clock-era code."""
    for c, label, _f, _l in CLOCK_ERAS:
        if c == code:
            return label
    return "unknown era"


# --------------------------------------------------------------------------
# The occultation model
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=8)
def _quadrature(n_nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes and weights mapped onto [0, pi] (cached: the
    fit evaluates the model tens of thousands of times)."""
    x, wq = np.polynomial.legendre.leggauss(n_nodes)
    return 0.5 * math.pi * (x + 1.0), 0.5 * math.pi * wq


def occultation_depth(z, k: float, u1: float = 0.0, u2: float = 0.0,
                      n_nodes: int = N_QUAD) -> np.ndarray:
    """Fraction of a limb-darkened disc's flux hidden by an opaque disc.

    The occulted star has unit radius and quadratic limb darkening
    ``I(mu) = 1 - u1 (1 - mu) - u2 (1 - mu)^2`` with ``mu = sqrt(1 - r^2)``.
    The occulter has radius ``k`` (same units) and its centre lies at
    projected distance ``z`` from the star's centre.  ``k`` may be smaller
    than one (a planet), about one (the M dwarf in an HW Vir system) or
    much larger (an M dwarf covering a white dwarf): nothing below assumes
    a small occulter.

    Method.  Split the stellar disc into thin rings of radius ``r``.  The
    fraction of a ring that lies inside the occulter is elementary geometry
    (the angle subtended by the chord where the two circles cross):

    * ring entirely inside the occulter (``r + z <= k``)      -> 1
    * ring entirely outside (``|r - z| >= k``)                -> 0
    * otherwise  ``arccos((r^2 + z^2 - k^2) / (2 r z)) / pi``

    The hidden flux is the intensity-weighted integral of that fraction
    over the disc, done in two pieces so that it is both accurate and a
    SMOOTH function of ``z`` (a least-squares fit differentiates it):

    * the fully covered core ``r < k - z`` (present only when ``k > z``)
      is integrated in closed form — with ``s = r^2`` the intensity is
      ``c0 + c1 sqrt(1 - s) + c2 (1 - s)``, whose integral is elementary;
    * the partially covered zone ``|z - k| < r < min(z + k, 1)`` is
      integrated by Gauss-Legendre quadrature in the angle ``phi`` of the
      substitution ``r = r_a + (r_b - r_a)(1 - cos phi)/2``, which piles
      nodes onto both ends of the zone where the integrand has square-root
      behaviour (the chord closing, and the stellar limb).

    With the default 48 nodes the uniform-disc case agrees with the
    analytic two-circle overlap to < 2e-6 of the stellar flux for every
    radius ratio used here (unit test).

    Returns an array shaped like ``z`` with values in [0, 1].
    """
    z_in = np.abs(np.asarray(z, dtype=float))
    zz = z_in.ravel()
    c0, c1, c2 = 1.0 - u1 - u2, u1 + 2.0 * u2, -u2

    def cumulative(s):
        """Integral of the intensity over disc area fraction 0..s."""
        s = np.clip(s, 0.0, 1.0)
        return (c0 * s + c1 * (2.0 / 3.0) * (1.0 - (1.0 - s) ** 1.5)
                + c2 * (s - 0.5 * s * s))

    norm = cumulative(np.array(1.0))
    # --- fully covered core -------------------------------------------
    r_core = np.clip(k - zz, 0.0, 1.0)
    hidden = cumulative(r_core ** 2)
    # --- partially covered zone ---------------------------------------
    r_a = np.minimum(np.abs(zz - k), 1.0)
    r_b = np.minimum(zz + k, 1.0)
    width = np.clip(r_b - r_a, 0.0, None)
    phi, wphi = _quadrature(n_nodes)                # nodes on [0, pi]
    r = r_a[:, None] + width[:, None] * 0.5 * (1.0 - np.cos(phi))[None, :]
    drdphi = width[:, None] * 0.5 * np.sin(phi)[None, :]
    mu = np.sqrt(np.clip(1.0 - r * r, 0.0, None))
    inten = c0 + c1 * mu + c2 * mu * mu
    zc = np.maximum(zz, 1e-12)[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosang = (r * r + zc * zc - k * k) / (2.0 * r * zc)
    frac = np.arccos(np.clip(np.nan_to_num(cosang, nan=1.0), -1.0, 1.0)) \
        / math.pi
    hidden = hidden + np.sum(inten * frac * 2.0 * r * drdphi
                             * wphi[None, :], axis=1)
    out = np.clip(hidden / norm, 0.0, 1.0)
    return out.reshape(z_in.shape)


def uniform_overlap_depth(z, k: float) -> np.ndarray:
    """Analytic hidden fraction of a UNIFORM unit disc (two-circle overlap).

    The textbook lens-area formula; used by the unit tests as the
    independent truth :func:`occultation_depth` must reproduce at
    ``u1 = u2 = 0``.
    """
    z = np.abs(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    full = z <= abs(1.0 - k)
    out[full] = min(k * k, 1.0)
    part = (~full) & (z < 1.0 + k)
    zp = z[part]
    a1 = np.arccos(np.clip((zp * zp + 1.0 - k * k) / (2.0 * zp), -1, 1))
    a2 = np.arccos(np.clip((zp * zp + k * k - 1.0) / (2.0 * zp * k), -1, 1))
    lens = a1 + k * k * a2 - 0.5 * np.sqrt(np.clip(
        (-zp + 1 + k) * (zp + 1 - k) * (zp - 1 + k) * (zp + 1 + k), 0, None))
    out[part] = lens / math.pi
    return out


def projected_separation(t, t0: float, period: float, a_over_r: float,
                         b: float) -> np.ndarray:
    """Sky-projected centre separation for a circular orbit, in units of
    the occulted star's radius.

    ``z^2 = (a/R)^2 sin^2(phi) + b^2 cos^2(phi)``, ``phi = 2 pi (t - t0)/P``
    — exact for a circular orbit with impact parameter ``b = (a/R) cos i``.
    On the far side of the orbit (``cos(phi) < 0``) the occulter is BEHIND
    the star; the separation is returned as +inf there so that half of the
    orbit contributes no occultation.  The function is even in
    ``t - t0``, which is the only property the mid-time fit relies on.
    """
    phi = 2.0 * math.pi * (np.asarray(t, dtype=float) - t0) / period
    z = np.sqrt((a_over_r * np.sin(phi)) ** 2 + (b * np.cos(phi)) ** 2)
    return np.where(np.cos(phi) > 0.0, z, np.inf)


def event_duration_d(period: float, a_over_r: float, k: float,
                     b: float) -> float:
    """First-to-last-contact duration (days) of the circular-orbit event.

    ``T14 = (P/pi) asin( sqrt((1+k)^2 - b^2) / sqrt((a/R)^2 - b^2) )``;
    returns 0 for a geometry that never touches (``b >= 1 + k``).
    """
    num = (1.0 + k) ** 2 - b * b
    den = a_over_r ** 2 - b * b
    if num <= 0 or den <= 0:
        return 0.0
    return period / math.pi * math.asin(min(1.0, math.sqrt(num / den)))


def event_model(t, exptime_d, t0: float, period: float, a_over_r: float,
                k: float, b: float, u1: float, u2: float,
                f1: float = 1.0) -> np.ndarray:
    """Normalised flux of one occultation event, exposure-averaged.

    ``flux = 1 - f1 * occultation_depth(z(t))``.  ``f1`` is the fraction of
    the system's light that comes from the occulted star (1 for a
    transiting planet; < 1 for a binary whose companion also shines).  Each
    exposure longer than :data:`SMEAR_MIN_EXPTIME_S` is averaged over
    :data:`N_SMEAR` sub-samples spanning it, because a 120-300 s exposure
    smears a one-minute white-dwarf ingress and an unsmeared model would
    put the mid-time where the SHARP curve fits best.
    """
    t = np.asarray(t, dtype=float)
    ex = np.broadcast_to(np.asarray(exptime_d, dtype=float), t.shape)
    if np.nanmax(ex) * SECONDS_PER_DAY < SMEAR_MIN_EXPTIME_S:
        z = projected_separation(t, t0, period, a_over_r, b)
        return 1.0 - f1 * _depth_finite(z, k, u1, u2)
    offs = (np.arange(N_SMEAR) + 0.5) / N_SMEAR - 0.5      # mid-point rule
    tt = t[:, None] + ex[:, None] * offs[None, :]
    z = projected_separation(tt, t0, period, a_over_r, b)
    d = _depth_finite(z.ravel(), k, u1, u2).reshape(z.shape)
    return 1.0 - f1 * d.mean(axis=1)


def _depth_finite(z: np.ndarray, k: float, u1: float, u2: float
                  ) -> np.ndarray:
    """:func:`occultation_depth` evaluated only where the discs can touch
    (skips the out-of-event bulk of a light curve, which is most of it)."""
    z = np.asarray(z, dtype=float)
    out = np.zeros(z.shape)
    touch = z < 1.0 + k
    if touch.any():
        out[touch] = occultation_depth(z[touch], k, u1, u2)
    return out


# --------------------------------------------------------------------------
# Ephemeris propagation
# --------------------------------------------------------------------------

def nearest_epoch(t: float, t0: float, period: float) -> int:
    """Integer cycle number of the event nearest to time ``t``."""
    return int(round((t - t0) / period))


def predict_event(t0: float, sig_t0: float, period: float, sig_p: float,
                  epoch: int, quad: float = 0.0, sig_quad: float = 0.0
                  ) -> tuple[float, float]:
    """Predicted event time and its 1-sigma uncertainty (days) at ``epoch``.

    ``T = t0 + P E + quad E^2``; the uncertainty adds the three terms in
    quadrature, i.e. it treats the published parameters as uncorrelated.
    That is the correct reading for an ephemeris quoted at its
    zero-covariance epoch (the ExoClock convention) and CONSERVATIVE-OR-
    EQUAL otherwise only when the published covariance is negative; where
    it matters the build cross-checks every prediction against a second,
    independent ephemeris instead of trusting one error bar.
    """
    e = float(epoch)
    t = t0 + period * e + quad * e * e
    sig = math.sqrt(sig_t0 ** 2 + (e * sig_p) ** 2 + (e * e * sig_quad) ** 2)
    return t, sig


def events_in_window(t_start: float, t_end: float, t0: float, period: float
                     ) -> list[int]:
    """Cycle numbers of every predicted event inside ``[t_start, t_end]``."""
    e_lo = int(math.ceil((t_start - t0) / period))
    e_hi = int(math.floor((t_end - t0) / period))
    return list(range(e_lo, e_hi + 1))


# --------------------------------------------------------------------------
# Light-curve assembly
# --------------------------------------------------------------------------

def native_peak_ratio(fwhm_px: float, binning: int) -> float:
    """Worst-case (brightest NATIVE pixel) / (binned-average pixel) ratio.

    Standing rule 4: saturation is judged at the target, per frame, in
    native pixels.  The 2x2-binned eras store the AVERAGE of four native
    pixels, so a clipped native pixel hides behind three unclipped
    neighbours (OA.E8).  For a Gaussian PSF of FWHM ``fwhm_px`` BINNED
    pixels whose peak sits on one native pixel, the brightest native pixel
    exceeds the mean of its bin-block by the factor returned here; a binned
    peak times this factor is the native peak to compare with the ADC
    ceiling.  ``binning = 1`` returns exactly 1.
    """
    if binning <= 1 or not fwhm_px or fwhm_px <= 0:
        return 1.0
    sigma_native = fwhm_px * binning / 2.3548
    offs = np.arange(binning, dtype=float)
    g = np.exp(-(offs[:, None] ** 2 + offs[None, :] ** 2)
               / (2.0 * sigma_native ** 2))
    return float(1.0 / g.mean())


def successive_difference_rms(y) -> float:
    """White-noise scatter from point-to-point differences (robust).

    ``1.4826 MAD(diff(y)) / sqrt(2)``: insensitive to a slow trend or to
    the event itself, so it can rank apertures and set point errors BEFORE
    any model is fitted (no forking path through the fit).
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 4:
        return float("nan")
    d = np.diff(y)
    return float(1.4826 * np.median(np.abs(d - np.median(d))) / math.sqrt(2))


def jump_steps(x, y, jump_px: float = JUMP_PX,
               min_segment: int = JUMP_MIN_SEGMENT,
               max_steps: int = JUMP_MAX_STEPS) -> np.ndarray:
    """Step-function regressors at the target's re-pointing jumps.

    Given the target's centroid track, returns an ``(n, m)`` array with one
    column per accepted jump: 0 before the jump, 1 from it onwards.  A jump
    is a frame-to-frame centroid displacement larger than ``jump_px``; a
    jump that would leave a segment shorter than ``min_segment`` points on
    either side of it (counted to the neighbouring accepted jump or the
    series end) is ignored, and if more than ``max_steps`` remain the
    largest displacements are kept.  ``m`` may be zero.

    These columns enter the baseline of the ADOPTED fit: the cause of the
    step is physical and its epoch is fixed by astrometry, so modelling it
    is not a choice made by looking at the light curve.  The fit without
    them is kept as a variant.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 2 * min_segment:
        return np.zeros((n, 0))
    d = np.hypot(np.diff(x), np.diff(y))
    cand = [int(i) + 1 for i in np.flatnonzero(d > jump_px)]
    cand.sort(key=lambda i: -d[i - 1])                 # largest first
    kept: list[int] = []
    for i in cand:
        edges = sorted(kept + [0, n])
        left = max(e for e in edges if e <= i)
        right = min(e for e in edges if e > i)
        if i - left >= min_segment and right - i >= min_segment:
            kept.append(i)
        if len(kept) >= max_steps:
            break
    kept.sort()
    out = np.zeros((n, len(kept)))
    for j, i in enumerate(kept):
        out[i:, j] = 1.0
    return out


def select_comparisons(flux: np.ndarray, usable: np.ndarray,
                       min_frac: float = COMP_MIN_FRAME_FRAC,
                       rms_factor: float = COMP_RMS_FACTOR,
                       min_comps: int = MIN_COMPS) -> np.ndarray:
    """Choose the comparison ensemble by iterative stability pruning.

    Parameters
    ----------
    flux
        ``(n_frames, n_comp)`` comparison-star fluxes.
    usable
        Same shape, True where the star was measured and unsaturated.

    Returns the boolean keep-mask over the comparison stars.  Rule: a star
    must be usable on ``min_frac`` of the frames; then, repeatedly, each
    remaining star's light curve is formed against the SUM OF THE OTHERS
    and the worst star is dropped while its successive-difference scatter
    exceeds ``rms_factor`` x the median scatter.  The target never enters,
    so the selection cannot learn anything about the event being timed.
    """
    flux = np.asarray(flux, dtype=float)
    usable = np.asarray(usable, dtype=bool)
    n_comp = flux.shape[1]
    keep = usable.mean(axis=0) >= min_frac
    while keep.sum() > max(min_comps, 2):
        idx = np.flatnonzero(keep)
        rms = np.empty(len(idx))
        for j, i in enumerate(idx):
            others = [o for o in idx if o != i]
            ok = usable[:, idx].all(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = flux[ok, i] / flux[ok][:, others].sum(axis=1)
            rms[j] = successive_difference_rms(ratio / np.nanmedian(ratio))
        med = float(np.nanmedian(rms))
        worst = int(np.nanargmax(rms))
        if not np.isfinite(med) or rms[worst] <= rms_factor * med:
            break
        keep[idx[worst]] = False
    if keep.sum() < min_comps:
        keep[:] = False
    assert keep.shape == (n_comp,)
    return keep


def differential_lightcurve(target_flux: np.ndarray, target_ok: np.ndarray,
                            comp_flux: np.ndarray, comp_ok: np.ndarray,
                            keep: np.ndarray) -> tuple[np.ndarray,
                                                       np.ndarray]:
    """Normalised target/ensemble flux ratio and the per-frame use mask.

    A frame is used when the target and EVERY kept comparison star were
    measured and unsaturated on it, and the ensemble flux is at least
    :data:`ENSEMBLE_MIN_FRACTION` of its series median (clouds).  Requiring
    every kept star on every used frame keeps the ensemble's composition
    fixed, so no frame-to-frame change of membership can masquerade as a
    step in the light curve.  The ratio is divided by its median over the
    used frames.
    """
    ens = comp_flux[:, keep].sum(axis=1)
    use = target_ok & comp_ok[:, keep].all(axis=1) & (target_flux > 0) \
        & (ens > 0)
    if use.any():
        use &= ens >= ENSEMBLE_MIN_FRACTION * np.median(ens[use])
    ratio = np.full(len(target_flux), np.nan)
    ratio[use] = target_flux[use] / ens[use]
    if use.any():
        ratio /= np.nanmedian(ratio[use])
    return ratio, use


# --------------------------------------------------------------------------
# The mid-time fit
# --------------------------------------------------------------------------

class EventShape:
    """The fixed and free shape parameters of one clock target's event.

    ``free`` names which of ``('k', 'a_over_r', 'f1')`` the fit may move;
    the mid-time and the polynomial baseline are always free.  A transit
    frees ``k`` (the depth) and keeps the literature duration; an eclipsing
    binary keeps ``k`` and frees ``a_over_r`` (the duration) and ``f1``
    (the depth), because for a deep eclipse those two are what the data
    constrain.  None of these choices can move the midpoint of a symmetric
    event; the build proves that empirically by re-fitting with the other
    choices and recording the shifts (``s3b_variants``).
    """

    def __init__(self, period: float, a_over_r: float, k: float, b: float,
                 u1: float, u2: float, f1: float = 1.0,
                 free: Sequence[str] = ("k",), baseline_order: int = 1):
        self.period = float(period)
        self.a_over_r = float(a_over_r)
        self.k = float(k)
        self.b = float(b)
        self.u1 = float(u1)
        self.u2 = float(u2)
        self.f1 = float(f1)
        self.free = tuple(free)
        self.baseline_order = int(baseline_order)

    def copy(self, **changes) -> "EventShape":
        kw = dict(period=self.period, a_over_r=self.a_over_r, k=self.k,
                  b=self.b, u1=self.u1, u2=self.u2, f1=self.f1,
                  free=self.free, baseline_order=self.baseline_order)
        kw.update(changes)
        return EventShape(**kw)

    def duration_d(self) -> float:
        return event_duration_d(self.period, self.a_over_r, self.k, self.b)

    def nominal_depth(self) -> float:
        """Mid-event depth of the nominal shape (fraction of total flux)."""
        return float(self.f1 * occultation_depth(
            np.array([self.b]), self.k, self.u1, self.u2)[0])


def _baseline_matrix(x: np.ndarray, order: int,
                     regressors: Optional[np.ndarray] = None) -> np.ndarray:
    """Design matrix of the out-of-event baseline.

    Columns: a polynomial in time of the given order, plus (optionally)
    one column per external regressor — per-frame quantities such as
    seeing FWHM or airmass, mean-subtracted and scaled to unit scatter so
    the least-squares problem stays well conditioned.  The regressors are
    NOT part of the adopted fit; they exist so the build can ask "does the
    mid-time move if the light curve is decorrelated against seeing?" and
    carry the answer as a model error.
    """
    base = np.vander(x, order + 1, increasing=True)
    if regressors is None:
        return base
    r = np.asarray(regressors, dtype=float)
    if r.ndim == 1:
        r = r[:, None]
    r = r - np.nanmean(r, axis=0)
    sd = np.nanstd(r, axis=0)
    r = np.nan_to_num(r / np.where(sd > 0, sd, 1.0))
    return np.column_stack([base, r])


def _solve_baseline(flux, w, occ, design):
    """Weighted linear least squares for the baseline (``design`` columns)
    that multiplies the occultation curve ``occ``; returns (coeffs,
    model)."""
    a = design * occ[:, None]
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(a * sw[:, None], flux * sw, rcond=None)
    return coef, a @ coef


def _unpack(shape: EventShape, theta: np.ndarray) -> dict:
    """Map the non-linear parameter vector onto named shape values."""
    vals = {"k": shape.k, "a_over_r": shape.a_over_r, "f1": shape.f1}
    for name, v in zip(shape.free, theta[1:]):
        vals[name] = float(v)
    return vals


def _model_for(theta, t, ex, shape: EventShape, flux, w, design):
    v = _unpack(shape, theta)
    occ = event_model(t, ex, float(theta[0]), shape.period, v["a_over_r"],
                      v["k"], shape.b, shape.u1, shape.u2, v["f1"])
    coef, model = _solve_baseline(flux, w, occ, design)
    return occ, coef, model


def grid_scan_t0(t, flux, err, exptime_d, shape: EventShape,
                 n_per_duration: int = 24, regressors=None
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Chi-square of the best fixed-shape fit on a grid of trial mid-times.

    The grid spans the WHOLE data window (first to last point) in steps of
    ``duration / n_per_duration``.  This is the blind part of the fit: the
    literature prediction is never used to choose where to look, so an
    event displaced by a clock error of any size that still leaves it in
    the window is found where it IS.  At each trial mid-time the nominal
    shape is held and only the polynomial baseline and a depth scale are
    solved (linear least squares).  Trials whose depth scale is not
    positive (a brightening) get chi-square = +inf.
    """
    t = np.asarray(t, dtype=float)
    flux = np.asarray(flux, dtype=float)
    w = 1.0 / np.asarray(err, dtype=float) ** 2
    dur = max(shape.duration_d(), 4.0 * float(np.median(np.diff(t))))
    step = dur / n_per_duration
    grid = np.arange(t.min(), t.max() + step, step)
    x = t - t.mean()
    base = _baseline_matrix(x, shape.baseline_order, regressors)
    sw = np.sqrt(w)
    chi2 = np.full(len(grid), np.inf)
    for i, t0 in enumerate(grid):
        occ = event_model(t, exptime_d, t0, shape.period, shape.a_over_r,
                          shape.k, shape.b, shape.u1, shape.u2, shape.f1)
        dip = 1.0 - occ                       # >= 0 inside the event
        if dip.max() <= 0:
            continue
        a = np.column_stack([base, -dip])     # flux ~ poly - s * dip
        coef, *_ = np.linalg.lstsq(a * sw[:, None], flux * sw, rcond=None)
        if coef[-1] <= 0:
            continue
        chi2[i] = float(np.sum(w * (flux - a @ coef) ** 2))
    return grid, chi2


def fit_event(t, flux, err, exptime_d, shape: EventShape,
              t0_start: Optional[float] = None, regressors=None) -> dict:
    """Fit one event's mid-time (and its free shape parameters).

    Parameters
    ----------
    t, flux, err, exptime_d
        Mid-exposure BJD_TDB, normalised flux, 1-sigma errors, and the
        exposure times in DAYS (scalar or per point).
    shape
        The :class:`EventShape` (fixed values, which are free, baseline
        order).
    t0_start
        Skip the blind grid and start the refinement here — used ONLY by
        the bootstrap and prayer-bead resamples of an already-located
        event, never for the measurement itself.
    regressors
        Optional ``(n, m)`` array of per-point external quantities added
        to the baseline (see :func:`_baseline_matrix`); variants only.

    Returns a dict with ``status`` and, when ``ok``: ``t0``, the fitted
    ``k`` / ``a_over_r`` / ``f1``, ``depth`` (mid-event fractional depth),
    ``chi2``, ``dof``, ``rms`` (unweighted residual scatter), ``model``,
    ``baseline`` (the polynomial alone), ``n_before`` / ``n_after`` (points
    outside first/last contact), ``t14_d`` and the formal ``t0_err_formal``
    from the Jacobian (reported, never adopted: see the bootstrap).
    """
    from scipy.optimize import least_squares

    t = np.asarray(t, dtype=float)
    flux = np.asarray(flux, dtype=float)
    err = np.asarray(err, dtype=float)
    ex = np.broadcast_to(np.asarray(exptime_d, dtype=float), t.shape)
    n = len(t)
    w = 1.0 / err ** 2
    x = t - t.mean()
    design = _baseline_matrix(x, shape.baseline_order, regressors)
    n_par = 1 + len(shape.free) + design.shape[1]
    if n < max(MIN_FIT_POINTS, n_par + 3):
        return {"status": STATUS_TOO_FEW, "n_points": n}

    if t0_start is None:
        grid, chi2 = grid_scan_t0(t, flux, err, ex, shape,
                                  regressors=regressors)
        if not np.isfinite(chi2).any():
            return {"status": STATUS_NO_DIP, "n_points": n}
        t0_start = float(grid[int(np.argmin(chi2))])

    # The optimiser works in SCALED variables: the mid-time as an offset
    # from the mean time in units of duration/20, each shape parameter as
    # a ratio to its nominal value.  This is not cosmetic.  A raw BJD is
    # ~2.46e6, and a finite-difference Jacobian steps a parameter by
    # ~1.5e-8 of its own magnitude — 53 MINUTES for a BJD — which destroys
    # the derivative with respect to the one parameter this whole stage
    # exists to measure.  In scaled variables the step is ~0.2 ms.
    dur = max(shape.duration_d(), 1e-4)
    t_ref = float(t.mean())
    t_unit = dur / 20.0
    nominal = [getattr(shape, name) for name in shape.free]

    # The model is evaluated on times RELATIVE to t_ref for the same
    # reason: at BJD 2.46e6 one unit in the last place of a double is
    # 40 microseconds, and a derivative taken across a few of those is
    # mostly rounding noise.
    tau = t - t_ref

    def to_theta(p):
        return np.array([p[0] * t_unit]
                        + [pj * nj for pj, nj in zip(p[1:], nominal)])

    p0 = [(t0_start - t_ref) / t_unit] + [1.0] * len(nominal)
    lo = [(t.min() - dur - t_ref) / t_unit] \
        + [_BOUNDS[name][0] for name in shape.free]
    hi = [(t.max() + dur - t_ref) / t_unit] \
        + [_BOUNDS[name][1] for name in shape.free]
    # f1 is a light fraction: it can never exceed 1.
    for j, name in enumerate(shape.free):
        if name == "f1":
            hi[1 + j] = min(hi[1 + j], 1.0 / nominal[j])
            p0[1 + j] = min(p0[1 + j], hi[1 + j])
    p0 = list(np.clip(p0, lo, hi))

    def resid(p):
        _occ, _coef, model = _model_for(to_theta(p), tau, ex, shape, flux,
                                        w, design)
        return (flux - model) / err

    # diff_step 1e-3 of a scaled unit = duration/20,000 (a fraction of a
    # second): large against rounding, small against any structure in the
    # light curve.  Tolerances are tightened because the chi-square
    # surface of a noisy transit is shallow near its minimum and the
    # default stopping rule leaves the mid-time tens of seconds short.
    sol = least_squares(resid, p0, bounds=(lo, hi), method="trf",
                        diff_step=1e-3, ftol=1e-12, xtol=1e-12, gtol=1e-12,
                        max_nfev=400)
    theta = to_theta(sol.x)
    occ, coef, model = _model_for(theta, tau, ex, shape, flux, w, design)
    theta = theta.copy()
    theta[0] += t_ref                      # back to absolute BJD_TDB
    v = _unpack(shape, theta)
    res = flux - model
    chi2 = float(np.sum(w * res ** 2))
    dof = n - n_par
    t14 = event_duration_d(shape.period, v["a_over_r"], v["k"], shape.b)
    t0 = float(theta[0])
    n_before = int(np.sum(t < t0 - t14 / 2.0))
    n_after = int(np.sum(t > t0 + t14 / 2.0))
    depth = float(v["f1"] * occultation_depth(
        np.array([shape.b]), v["k"], shape.u1, shape.u2)[0])
    # Formal error from the Jacobian (J^T J)^-1 — printed beside the
    # bootstrap so a reader can see how much red noise inflates it.
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac)
        t0_err_formal = float(math.sqrt(max(cov[0, 0], 0.0)) * t_unit)
    except np.linalg.LinAlgError:
        t0_err_formal = float("nan")
    status = STATUS_OK
    if min(n_before, n_after) < MIN_SIDE_POINTS:
        status = STATUS_ONE_SIDED
    return {"status": status, "n_points": n, "t0": t0, "k": v["k"],
            "a_over_r": v["a_over_r"], "f1": v["f1"], "depth": depth,
            "chi2": chi2, "dof": dof, "rms": float(np.std(res)),
            "model": model, "baseline": design @ coef,
            "n_before": n_before, "n_after": n_after, "t14_d": t14,
            "t0_err_formal": t0_err_formal, "theta": theta}


#: Multiplicative bounds on the free shape parameters, relative to their
#: nominal values (wide: they exist to keep the optimiser out of
#: unphysical corners, not to constrain the answer).
_BOUNDS = {"k": (0.2, 3.0), "a_over_r": (0.4, 2.5), "f1": (0.05, 20.0)}


def clip_outliers(resid: np.ndarray, n_sigma: float = CLIP_SIGMA
                  ) -> np.ndarray:
    """Boolean keep-mask: |resid - median| <= n_sigma x 1.4826 MAD."""
    resid = np.asarray(resid, dtype=float)
    med = np.median(resid)
    mad = 1.4826 * np.median(np.abs(resid - med))
    if not np.isfinite(mad) or mad <= 0:
        return np.ones(len(resid), dtype=bool)
    return np.abs(resid - med) <= n_sigma * mad


def block_length_points(t: np.ndarray, block_min: float = BOOT_BLOCK_MIN
                        ) -> int:
    """Number of points in one bootstrap block of ``block_min`` minutes."""
    cadence = float(np.median(np.diff(np.asarray(t, dtype=float))))
    if not np.isfinite(cadence) or cadence <= 0:
        return 1
    return int(max(1, round(block_min / 1440.0 / cadence)))


def resample_blocks(resid: np.ndarray, block: int,
                    rng: np.random.Generator) -> np.ndarray:
    """One moving-block bootstrap resample of a residual series.

    Blocks of ``block`` consecutive residuals are drawn with replacement
    from all overlapping positions and concatenated to the original
    length, so correlations shorter than a block survive the resampling.
    """
    n = len(resid)
    block = int(max(1, min(block, n)))
    n_blocks = int(math.ceil(n / block))
    starts = rng.integers(0, n - block + 1, size=n_blocks)
    out = np.concatenate([resid[s:s + block] for s in starts])
    return out[:n]


def bootstrap_t0(t, flux, err, exptime_d, shape: EventShape, fit: dict,
                 n_boot: int = N_BOOTSTRAP, seed: int = BOOT_SEED,
                 block_min: float = BOOT_BLOCK_MIN, regressors=None) -> dict:
    """Moving-block residual bootstrap of the mid-time.

    Synthetic light curves = best-fit model + block-resampled residuals;
    each is re-fitted from the best-fit mid-time (the event has already
    been located blind; the resamples measure how far noise moves it).

    Returns ``sigma`` (half the 16-84 percentile range, days), ``bias``
    (median of the resampled mid-times minus the fit, days; signed) and
    ``n_ok``.
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(t, dtype=float)
    res = np.asarray(flux, dtype=float) - fit["model"]
    block = block_length_points(t, block_min)
    start_shape = shape.copy(**{name: fit[name] for name in shape.free})
    t0s = []
    for _ in range(n_boot):
        synth = fit["model"] + resample_blocks(res, block, rng)
        f = fit_event(t, synth, err, exptime_d, start_shape,
                      t0_start=fit["t0"], regressors=regressors)
        if f["status"] in (STATUS_OK, STATUS_ONE_SIDED):
            t0s.append(f["t0"])
    if len(t0s) < max(20, n_boot // 4):
        return {"sigma": float("nan"), "bias": float("nan"),
                "n_ok": len(t0s), "block_points": block}
    t0s = np.asarray(t0s)
    lo, med, hi = np.percentile(t0s, [15.865, 50.0, 84.135])
    return {"sigma": float((hi - lo) / 2.0),
            "bias": float(med - fit["t0"]), "n_ok": len(t0s),
            "block_points": block}


def prayer_bead_t0(t, flux, err, exptime_d, shape: EventShape, fit: dict,
                   max_shifts: int = 60, regressors=None) -> dict:
    """Residual-permutation ("prayer bead") scatter of the mid-time.

    The residual series is cyclically shifted by up to ``max_shifts``
    evenly spaced offsets and re-added to the model; the standard
    deviation of the re-fitted mid-times is the alternative red-noise
    error printed beside the bootstrap.  The build ADOPTS the larger of
    the two (a rule fixed here, before any number was seen).
    """
    t = np.asarray(t, dtype=float)
    res = np.asarray(flux, dtype=float) - fit["model"]
    n = len(res)
    shifts = np.unique(np.linspace(1, n - 1, min(max_shifts, n - 1)
                                   ).astype(int))
    start_shape = shape.copy(**{name: fit[name] for name in shape.free})
    t0s = []
    for s in shifts:
        synth = fit["model"] + np.roll(res, int(s))
        f = fit_event(t, synth, err, exptime_d, start_shape,
                      t0_start=fit["t0"], regressors=regressors)
        if f["status"] in (STATUS_OK, STATUS_ONE_SIDED):
            t0s.append(f["t0"])
    if len(t0s) < 5:
        return {"sigma": float("nan"), "n_ok": len(t0s)}
    return {"sigma": float(np.std(t0s)), "n_ok": len(t0s)}


def inject_and_recover(t, flux, err, exptime_d, shape: EventShape,
                       fit: dict, shift_s: float, n_draws: int,
                       seed: int, regressors=None) -> dict:
    """Would a clock error of ``shift_s`` seconds have been recovered?

    The fitted event is REMOVED from the data (leaving baseline x
    residuals), the same event is re-injected displaced by ``shift_s``,
    the residuals are block-resampled, and the full measurement — blind
    grid included — is repeated.  Returns the mean SIGNED recovery error
    (recovered - injected, seconds; standing rule 3 forbids a median of
    absolute values), its scatter, and how many draws returned a usable
    two-sided fit.  A shift that pushes the event out of the observing
    window is reported as unrecoverable (``n_ok = 0``), which is itself
    the answer to "could this series have seen it?".
    """
    rng = np.random.default_rng(seed)
    t = np.asarray(t, dtype=float)
    ex = np.broadcast_to(np.asarray(exptime_d, dtype=float), t.shape)
    res = np.asarray(flux, dtype=float) - fit["model"]
    block = block_length_points(t)
    occ_new = event_model(t, ex, fit["t0"] + shift_s / SECONDS_PER_DAY,
                          shape.period, fit["a_over_r"], fit["k"], shape.b,
                          shape.u1, shape.u2, fit["f1"])
    clean = fit["baseline"] * occ_new
    errs = []
    for _ in range(n_draws):
        synth = clean + resample_blocks(res, block, rng)
        f = fit_event(t, synth, err, ex, shape, regressors=regressors)
        if f["status"] != STATUS_OK:
            continue
        rec = (f["t0"] - fit["t0"]) * SECONDS_PER_DAY
        errs.append(rec - shift_s)
    if not errs:
        return {"shift_s": shift_s, "n_ok": 0, "n_draws": n_draws,
                "bias_s": float("nan"), "scatter_s": float("nan")}
    return {"shift_s": shift_s, "n_ok": len(errs), "n_draws": n_draws,
            "bias_s": float(np.mean(errs)),
            "scatter_s": float(np.std(errs)) if len(errs) > 1 else 0.0}


# --------------------------------------------------------------------------
# Combining O - C values
# --------------------------------------------------------------------------

def combine_oc(oc_s: Sequence[float], sig_s: Sequence[float]) -> dict:
    """Combine independent O - C values; formal AND scatter-based errors.

    Returns the inverse-variance weighted mean and its formal error, the
    chi-square of the values about that mean with its degrees of freedom,
    and — independently of every quoted error bar — the plain mean and the
    standard error of the mean from the observed scatter (``N >= 2``).

    No error is rescaled by ``sqrt(chi2_nu)`` and no ``max(chi2_nu, 1)``
    is applied (standing rule 1): both numbers are reported, and
    ``sigma_adopted`` is simply the LARGER of the formal and the
    scatter-based error, so an over-dispersed set cannot hide behind
    optimistic error bars and an under-dispersed one cannot shrink its
    quoted error below what its own bars support.
    """
    oc = np.asarray(oc_s, dtype=float)
    sig = np.asarray(sig_s, dtype=float)
    n = len(oc)
    if n == 0:
        return {"n": 0}
    w = 1.0 / sig ** 2
    mean = float(np.sum(w * oc) / np.sum(w))
    err = float(1.0 / math.sqrt(np.sum(w)))
    chi2 = float(np.sum(w * (oc - mean) ** 2))
    out = {"n": n, "wmean": mean, "wmean_err": err, "chi2": chi2,
           "dof": n - 1, "mean": float(np.mean(oc)),
           "scatter_err": None, "sigma_adopted": err}
    if n >= 2:
        sem = float(np.std(oc, ddof=1) / math.sqrt(n))
        out["scatter_err"] = sem
        out["sigma_adopted"] = max(err, sem)
    return out


def era_verdict(oc_s: float, sigma_s: float,
                accept_s: float = ACCEPT_OC_S) -> str:
    """Verdict string for one clock era against the F-8 criterion.

    * ``PASS`` — the whole 2-sigma interval lies inside +-accept: the
      clock is right to better than the criterion.
    * ``OFFSET MEASURED`` — the interval excludes zero at 3 sigma AND the
      value is outside the criterion: a real offset, to be carried.
    * ``PASS (central value); 2-sigma interval reaches the criterion`` —
      the estimate is inside but its error is too large to certify it.
    * ``UNRESOLVED`` — everything else.
    """
    if abs(oc_s) + 2.0 * sigma_s < accept_s:
        return "PASS"
    if abs(oc_s) >= accept_s and abs(oc_s) > 3.0 * sigma_s:
        return "OFFSET MEASURED"
    if abs(oc_s) < accept_s:
        return "PASS (central value); 2-sigma interval reaches the criterion"
    return "UNRESOLVED"
