"""Pure S3 timing logic: the shared time axis every MACRO paper depends on.

This module answers, per canonical science frame, ONE question: *at what
Barycentric Dynamical Time did the midpoint of this exposure occur at the
solar-system barycenter?* — i.e. it computes mid-exposure BJD_TDB from the
header's UTC exposure start, from scratch.  Header JD-HELIO is NEVER used
(ROADMAP S3 rule); the audit in ``build_s3_timing.py`` shows why it must
not be (it is a heliocentric-only, mid-exposure MaxIm value with its own
approximations — evidence in ``docs/pipeline/s3_timing.html``).

Everything here is either a *pure function* (no I/O, unit-tested on
hand-built cases) or a *thin, deterministic wrapper* around astropy's time
machinery (vectorized, same inputs -> same outputs).  The build script
(``pipeline/scripts/build_s3_timing.py``) wires these functions to the
manifest; the report renderer (``macro_core.report_s3``) reads only the
database tables the build wrote.

The timing conventions implemented here, with their evidence:

1.  Header JD == DATE-OBS == UTC exposure START.  Verified over all
    198,289 canonical frames with both cards: |JD - DATE-OBS| <= 2 ms
    except a single reduced-tree frame whose JD the reduction pipeline
    re-stamped (the known S0b ``stem_jd_drift`` behavior).  START (not
    mid/end) is proven by the header's own JD-HELIO: across eras,
    JD-HELIO - JD - (our heliocentric correction) == EXPTIME/2 to ~0.1 s
    — the acquisition software computed its heliocentric stamp at
    start + EXPTIME/2, so the base stamp is the start.
2.  Mid-exposure = start + EXPTIME/2 for every readout family, including
    StackPro (see :data:`STACKPRO_DEADTIME_BOUND_S` for the worst case).
3.  BJD_TDB = (UTC start + EXPTIME/2) -> TDB scale -> + barycentric light
    travel time toward the frame's sky position, computed with astropy
    ``Time.light_travel_time`` at the Winer EarthLocation with a JPL DE
    ephemeris.  Sub-second end-to-end (see the error budget below).
4.  The sky position used is the FRAME CENTER (manifest ra_deg/dec_deg).
    The barycentric correction changes by <= ~1.3 s from frame center to
    a corner (half-diagonal ~26 arcmin x 499 s/rad); a paper that needs
    sub-second absolute times for an off-center object must recompute
    with :func:`bjd_tdb_from_utc` at the object's own coordinates — the
    ``frame_times`` table stores every input needed to do that.

Error budget of convention 3 (each term, worst case):
  * DE440s barycenter geometry ............ < 1 ms
  * UTC->TDB (leap seconds + 32.184 s + periodic) ... exact to < 1 ms
  * Earth-rotation observer term (Winer vs geocenter)  < 21.3 ms total;
    ignoring UT1-UTC (< 0.9 s of rotation) perturbs it by < 1.5 us
  * header DATE-OBS stamping (software clock write) ... ~10 ms class
  * StackPro mid-time policy .............. < 0.12 s (bounded below)
  * frame-center vs object position ....... < 1.3 s at a frame corner,
    ~0 for the targeted object near center (recompute to remove)
"""

from __future__ import annotations

import math
import re
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Site and policy constants (single source of truth — the report and the
# build script interpolate these; changing one here changes everything).
# --------------------------------------------------------------------------

#: Winer Observatory (ROADMAP S3; header SITELAT/SITELONG agree: 31 39 56,
#: -110 36 07).  Geodetic degrees, meters above the WGS84 ellipsoid — the
#: ~30 m geoid/ellipsoid ambiguity moves the barycentric correction by
#: ~0.1 us and is irrelevant at our precision target.
WINER_LAT_DEG = 31.6656
WINER_LON_DEG = -110.6018
WINER_ALT_M = 1515.0

#: Number of summed sub-reads in a StackPro frame (S2 detector campaign:
#: three independent PTC ratios all give N_sub = 16).
N_SUB_STACKPRO = 16

#: Empirical bound on the TOTAL internal dead time of one StackPro frame
#: (all 15 inter-sub-read boundaries together), in seconds.  Derived from
#: the archive itself (build stage ``cadence``): in back-to-back StackPro
#: series the minimum observed gap between consecutive exposure starts
#: minus EXPTIME is 0.24 s (0.5 s series, n=83 gaps) — the frame's true
#: wall-clock span therefore cannot exceed EXPTIME + 0.24 s, so the sum of
#: any inter-sub gaps is <= 0.24 s.  The sub-read readout time is a fixed
#: property of the sensor, independent of exposure length, so the bound
#: from the shortest series applies to every StackPro exposure.
STACKPRO_DEADTIME_BOUND_S = 0.24

#: Worst-case error of the StackPro mid-time policy (start + EXPTIME/2),
#: in seconds: if the whole dead-time bound sat between sub-reads, the true
#: photon-weighted midpoint would shift late by at most half the bound.
#: Independent corroboration: MaxIm's own JD-HELIO on a 1024 s StackPro
#: frame equals start + 512 s + heliocentric correction to 0.11 s — the
#: acquisition software also treats EXPTIME as the contiguous total span.
STACKPRO_MID_WORST_CASE_S = STACKPRO_DEADTIME_BOUND_S / 2.0

#: Ephemeris preference order.  ``de440s`` (JPL, 1849-2150, ~10 MB) is the
#: precision choice; ``builtin`` (ERFA analytic, ~km-level Earth position
#: = sub-ms timing) is the offline fallback.  The build records which one
#: actually ran in ``s3_build_meta`` and in every row's ``bjd_method``.
EPHEMERIS_PREFERENCE: tuple[str, ...] = ("de440s", "builtin")

#: Mid-time method identifiers written into ``frame_times.mid_method``.
MID_PLAIN = "start_plus_half_exptime"
MID_STACKPRO = "stackpro_sum_midpoint_half_exptime"
MID_NO_JD = "no_jd"
MID_EXPTIME_NONPOS = "exptime_nonpos_start_used"

#: The code-version string recorded in ``s3_build_meta`` and in every
#: ``frame_times`` row (kept here rather than ``macro_core/__init__`` so S3
#: work does not touch files a concurrent stage may be editing).
S3_CODE_VERSION = "S3 v1.0 (2026-08-18)"


# --------------------------------------------------------------------------
# DATE-OBS parsing (pure)
# --------------------------------------------------------------------------

#: Fractional-seconds tail of an ISO timestamp.  Headers write 0-3 fraction
#: digits ('...T03:53:57.60'); Python 3.10's ``fromisoformat`` accepts only
#: exactly 3 or 6, so :func:`parse_date_obs` normalizes the tail first.
_FRAC_RE = re.compile(r"\.(\d+)$")

#: JD of the Unix epoch — same constant the S0 manifest uses; duplicated
#: here (with the same value, asserted in tests) to keep this module
#: importable on its own.
UNIX_EPOCH_JD = 2440587.5


def parse_date_obs(date_obs: Optional[str]) -> Optional[float]:
    """Parse a header DATE-OBS string into a UTC Julian Date.

    Handles every format observed in the archive's 198k headers:
    ``YYYY-MM-DDThh:mm:ss`` with 0, 1, 2, or 3 fractional-second digits.
    Returns ``None`` for missing/blank/unparseable input rather than
    raising — the caller records the failure, it never crashes the build.

    Pure arithmetic via the Unix epoch (no astronomy library): DATE-OBS
    carries millisecond precision at best, far above where UTC subtleties
    (leap seconds inside an exposure) could bite a comparison against the
    header JD, which is what this function exists for.
    """
    if date_obs is None or not str(date_obs).strip():
        return None
    text = str(date_obs).strip()
    # Normalize the fractional tail to exactly 6 digits for fromisoformat.
    m = _FRAC_RE.search(text)
    if m:
        text = text[: m.start()] + "." + (m.group(1) + "000000")[:6]
    from datetime import datetime, timezone
    try:
        moment = datetime.fromisoformat(text).replace(tzinfo=timezone.utc)
    except ValueError:
        return None
    return moment.timestamp() / 86400.0 + UNIX_EPOCH_JD


# --------------------------------------------------------------------------
# Mid-exposure policy (pure)
# --------------------------------------------------------------------------

def is_stackpro(readoutm: Optional[str]) -> bool:
    """True when the readout mode is a StackPro (on-camera sum) family.

    Matched on the substring 'StackPro' in READOUTM — the era registry
    contains exactly one such family ('High Gain StackPro', 3 pinned eras);
    a substring match keeps any future 'Low Gain StackPro' frames on the
    correct policy automatically.
    """
    return "stackpro" in (readoutm or "").lower()


def mid_method_for(readoutm: Optional[str], jd: Optional[float],
                   exptime_s: Optional[float]) -> str:
    """Return the ``frame_times.mid_method`` identifier for one frame.

    * no usable JD           -> :data:`MID_NO_JD` (no time axis at all);
    * EXPTIME missing or <=0 -> :data:`MID_EXPTIME_NONPOS` (the start is
      used as the mid; the S0 ``exptime_nonpos`` QC flag already marks
      these frames, 156 canonical science cases);
    * StackPro readout       -> :data:`MID_STACKPRO` (same arithmetic,
      distinct label so the policy's worst case is queryable);
    * everything else        -> :data:`MID_PLAIN`.
    """
    if jd is None or (isinstance(jd, float) and math.isnan(jd)):
        return MID_NO_JD
    if exptime_s is None or (isinstance(exptime_s, float)
                             and math.isnan(exptime_s)) or exptime_s <= 0:
        return MID_EXPTIME_NONPOS
    return MID_STACKPRO if is_stackpro(readoutm) else MID_PLAIN


def jd_utc_mid(jd_start: Optional[float], exptime_s: Optional[float],
               readoutm: Optional[str] = None
               ) -> tuple[Optional[float], str]:
    """Mid-exposure UTC JD for one frame, plus the method identifier.

    Every readout family uses start + EXPTIME/2:

    * **Plain frames** — DATE-OBS/JD is the UTC exposure start (evidence:
      module docstring, convention 1) and EXPTIME the shutter-open span,
      so the midpoint is exact up to the header's own stamping precision.
    * **StackPro frames** — a StackPro frame is the on-camera SUM of
      :data:`N_SUB_STACKPRO` = 16 sub-reads (S2).  ASSUMPTION, stated:
      the sub-reads are contiguous and together span EXPTIME, so the
      photon-weighted midpoint of the sum is start + EXPTIME/2.  The
      cadence of back-to-back StackPro series bounds any violation: total
      internal dead time <= :data:`STACKPRO_DEADTIME_BOUND_S` = 0.24 s,
      hence worst-case mid-time error
      <= :data:`STACKPRO_MID_WORST_CASE_S` = 0.12 s (late).

    EXPTIME <= 0 (a header pathology, never a real exposure) falls back to
    the start instant, labeled :data:`MID_EXPTIME_NONPOS`.
    """
    method = mid_method_for(readoutm, jd_start, exptime_s)
    if method == MID_NO_JD:
        return None, method
    if method == MID_EXPTIME_NONPOS:
        return float(jd_start), method
    return float(jd_start) + float(exptime_s) / 2.0 / 86400.0, method


def worst_case_mid_error_s(readoutm: Optional[str]) -> float:
    """Worst-case mid-time policy error (seconds) for a readout family.

    Plain families: 0 by construction (the policy IS the definition of the
    midpoint; residual uncertainty is the header stamping precision, which
    belongs to the DATE-OBS card, not to the policy).  StackPro: the
    empirically bounded :data:`STACKPRO_MID_WORST_CASE_S`.
    """
    return STACKPRO_MID_WORST_CASE_S if is_stackpro(readoutm) else 0.0


# --------------------------------------------------------------------------
# BJD_TDB / HJD (thin deterministic astropy wrappers, vectorized)
# --------------------------------------------------------------------------

def winer_location():
    """The Winer Observatory ``EarthLocation`` (astropy), from the module
    constants.  A function (not a module-level object) so importing this
    module never touches astropy's lazy machinery."""
    from astropy import units as u
    from astropy.coordinates import EarthLocation
    return EarthLocation(lat=WINER_LAT_DEG * u.deg,
                         lon=WINER_LON_DEG * u.deg,
                         height=WINER_ALT_M * u.m)


def resolve_ephemeris(preference: Sequence[str] = EPHEMERIS_PREFERENCE
                      ) -> str:
    """First ephemeris in ``preference`` that actually loads on this host.

    ``de440s`` needs the jplephem package and a one-time kernel download
    (cached by astropy); ``builtin`` always works.  The chosen name is
    recorded in the build metadata so a re-run on another host cannot
    silently change precision class without leaving a trace.
    """
    from astropy.coordinates import solar_system_ephemeris
    for name in preference:
        try:
            with solar_system_ephemeris.set(name):
                pass
            return name
        except Exception:
            continue
    return "builtin"


def bjd_tdb_from_utc(jd_utc, ra_deg, dec_deg, ephemeris: str = "de440s",
                     location=None):
    """Convert UTC JD(s) at Winer to mid-frame BJD_TDB toward given sky
    position(s).  The core S3 conversion — everything else is bookkeeping.

    Parameters
    ----------
    jd_utc, ra_deg, dec_deg
        Scalars or equal-length arrays: the (mid-exposure) UTC Julian Date
        and the ICRS position the light-travel correction points at.
    ephemeris
        Solar-system ephemeris name for astropy (``de440s``/``builtin``).
    location
        ``EarthLocation`` override; default = Winer.  Pass ``None``
        explicitly-with-geocenter semantics is NOT supported: tests that
        want the geocenter pass ``location=EarthLocation.from_geocentric
        (0,0,0, unit='m')`` — an explicit choice, never an accident.

    Returns
    -------
    (bjd_tdb, ltt_s, tdb_minus_utc_s)
        Arrays (or scalars, matching the input shape): the barycentric
        JD in the TDB scale; the light-travel term in seconds; and the
        TDB-minus-UTC scale offset in seconds (leap seconds + 32.184 s +
        periodic terms) — stored per frame so any row can be undone.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord, solar_system_ephemeris
    from astropy.time import Time
    loc = location if location is not None else winer_location()
    t = Time(np.asarray(jd_utc, dtype=float), format="jd", scale="utc",
             location=loc)
    target = SkyCoord(ra=np.asarray(ra_deg, dtype=float) * u.deg,
                      dec=np.asarray(dec_deg, dtype=float) * u.deg)
    with solar_system_ephemeris.set(ephemeris):
        ltt = t.light_travel_time(target, kind="barycentric")
    tdb = t.tdb
    bjd = tdb.jd + ltt.jd            # JD arithmetic: adds the ltt days
    tdb_minus_utc = (tdb.jd - t.jd) * 86400.0
    return bjd, ltt.sec, tdb_minus_utc


def hjd_utc_from_utc(jd_utc, ra_deg, dec_deg, ephemeris: str = "de440s",
                     location=None):
    """Heliocentric JD in the UTC scale — the header JD-HELIO convention.

    Exists ONLY for the S3 audit (comparing header JD-HELIO against an
    independent computation) and for reading old literature ephemerides
    quoted in HJD.  Science time stamps use :func:`bjd_tdb_from_utc`;
    HJD_UTC differs from BJD_TDB by up to ~73 s (scale offset) plus up to
    ~4 s (helio vs bary geometry) and has no place in a paper's time axis.
    """
    from astropy import units as u
    from astropy.coordinates import SkyCoord, solar_system_ephemeris
    from astropy.time import Time
    loc = location if location is not None else winer_location()
    t = Time(np.asarray(jd_utc, dtype=float), format="jd", scale="utc",
             location=loc)
    target = SkyCoord(ra=np.asarray(ra_deg, dtype=float) * u.deg,
                      dec=np.asarray(dec_deg, dtype=float) * u.deg)
    with solar_system_ephemeris.set(ephemeris):
        ltt = t.light_travel_time(target, kind="heliocentric")
    return t.jd + ltt.jd, ltt.sec


# --------------------------------------------------------------------------
# Eclipse-timing helpers for the clock validation (pure)
# --------------------------------------------------------------------------

def fold_phase(jd, epoch: float, period: float):
    """Orbital phase in [-0.5, +0.5) with the eclipse ephemeris at 0.

    Vectorized; plain numpy.  ``epoch`` and ``jd`` must be on the SAME
    time standard (the clock stage compares HJD-based VSX epochs against
    our own HJD_UTC mid-times — like against like).
    """
    ph = (np.asarray(jd, dtype=float) - epoch) / period % 1.0
    return np.where(ph >= 0.5, ph - 1.0, ph)


def fit_eclipse_offset(phases, dmags, errs=None,
                       ph0_grid=None, width_grid=None
                       ) -> dict:
    """Fit a symmetric eclipse dip to phase-folded differential photometry.

    Model: ``dmag(ph) = depth * exp(-(ph - ph0)^2 / (2 w^2))`` on a
    per-call zero baseline (the caller subtracts each night/filter's
    out-of-eclipse median first).  A Gaussian is not the true shape of an
    EA minimum, but for a ROUGH mid-time from sparse snapshots only the
    SYMMETRY matters: a symmetric template fit to a symmetric dip finds
    the center without modeling limb darkening or contact points.

    Deterministic grid search over (ph0, w) with the depth solved by
    weighted linear least squares at each grid point; the returned
    ``ph0_err`` is the half-width of the ``chi2 <= chi2_min + 1`` interval
    along the ph0 axis (profile over w), floored at one grid step.

    Returns a dict: ph0, ph0_err, depth, width, chi2_min, n_points.
    """
    ph = np.asarray(phases, dtype=float)
    dm = np.asarray(dmags, dtype=float)
    w_ = (np.ones_like(dm) if errs is None
          else 1.0 / np.maximum(np.asarray(errs, dtype=float), 1e-4) ** 2)
    if ph0_grid is None:
        ph0_grid = np.arange(-0.06, 0.0601, 0.0005)
    if width_grid is None:
        width_grid = np.arange(0.005, 0.0451, 0.0025)
    best = None
    chi2_by_ph0 = np.full(len(ph0_grid), np.inf)
    for i, ph0 in enumerate(ph0_grid):
        for w in width_grid:
            model = np.exp(-((ph - ph0) ** 2) / (2.0 * w ** 2))
            denom = float(np.sum(w_ * model * model))
            if denom <= 0:
                continue
            depth = float(np.sum(w_ * model * dm)) / denom
            if depth <= 0:           # an eclipse dims the star; reject
                continue             # inverted "brightenings" outright
            chi2 = float(np.sum(w_ * (dm - depth * model) ** 2))
            if chi2 < chi2_by_ph0[i]:
                chi2_by_ph0[i] = chi2
            if best is None or chi2 < best["chi2_min"]:
                best = {"ph0": float(ph0), "depth": depth,
                        "width": float(w), "chi2_min": chi2,
                        "n_points": int(len(ph))}
    if best is None:
        return {"ph0": None, "ph0_err": None, "depth": None, "width": None,
                "chi2_min": None, "n_points": int(len(ph))}
    # Profile-likelihood interval on ph0: all grid ph0 whose best-over-w
    # chi2 lies within 1 of the global minimum.
    ok = chi2_by_ph0 <= best["chi2_min"] + 1.0
    step = float(ph0_grid[1] - ph0_grid[0])
    ph0_err = max((ph0_grid[ok].max() - ph0_grid[ok].min()) / 2.0, step) \
        if ok.any() else step
    best["ph0_err"] = float(ph0_err)
    return best
