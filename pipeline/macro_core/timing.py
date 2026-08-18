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
    re-stamped (the known S0b ``stem_jd_drift`` behavior).  That test has
    a blind spot it cannot see past — it catches a re-stamp that moved
    ONE card, and the reduction pipeline sometimes moved BOTH — so the
    build ALSO cross-checks every exposure against its other copy via
    S0b's ``raw_reduced_links``; 46 disagree by more than
    :data:`JD_SIBLING_DISAGREE_S` and the reduced copy's BJD is withdrawn
    (see ``s3_time_outliers``).  START (not
    mid/end) is proven by the header's own JD-HELIO: across eras,
    JD-HELIO - JD - (our heliocentric correction) == EXPTIME/2 to ~0.1 s
    — the acquisition software computed its heliocentric stamp at
    start + EXPTIME/2, so the base stamp is the start.  THE EXCEPTION,
    stated because it is not proven: the 2026 ``pyscope`` eras write no
    JD-HELIO card at all and copy DATE-OBS verbatim into TELUT, so
    NEITHER probe exists for them; their rows carry
    :data:`START_UNVERIFIED` in ``frame_times.start_evidence``.
2.  Mid-exposure = start + EXPTIME/2 for every readout family, including
    StackPro (see :data:`STACKPRO_DEADTIME_BOUND_S` for the worst case,
    which is SECONDS, not milliseconds — read that constant before using
    a StackPro time for sub-second work).
3.  BJD_TDB = (UTC start + EXPTIME/2) -> TDB scale -> + barycentric light
    travel time toward the frame's sky position, computed with astropy
    ``Time.light_travel_time`` at the Winer EarthLocation with a JPL DE
    ephemeris.  The conversion itself is sub-millisecond; the INPUT
    stamp carries the rest of the budget below.
4.  The sky position used is the FRAME CENTER (manifest ra_deg/dec_deg).
    The barycentric correction changes by up to ~3.8 s from frame center
    to a corner on the widest field in the archive (half-diagonal ~26
    arcmin x :data:`LIGHT_TIME_PER_RAD_S`; the exact per-era ceiling is
    measured by the header-audit stage and printed in the report); a
    paper that needs sub-second absolute times for an off-center object
    must recompute with :func:`bjd_tdb_from_utc` at the object's own
    coordinates — the ``frame_times`` table stores every input needed to
    do that.

Error budget of convention 3 (each term, worst case):
  * DE440s barycenter geometry ............ < 1 ms
  * UTC->TDB (leap seconds + 32.184 s + periodic) ... exact to < 1 ms
  * Earth-rotation observer term (Winer vs geocenter)  < 21.3 ms total;
    ignoring UT1-UTC (< 0.9 s of rotation) perturbs it by < 1.5 us
  * header DATE-OBS stamping (software clock write) ... ~10 ms class for
    every MaxIm family (they write milliseconds), but 1 SECOND for the
    2026 pyscope eras, which stamp whole seconds only — up to 0.5 s of
    systematic if pyscope truncates rather than rounds
  * StackPro mid-time policy .............. <= STACKPRO_MID_WORST_CASE_S
    (seconds; bounded from cadence, see below)
  * frame-center vs object position ....... up to ~3.8 s at a corner of
    the widest field, ~0 for the targeted object near center (recompute
    at the object's coordinates to remove it entirely)
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
#: (all 15 inter-sub-read boundaries together), in seconds.  Measured by
#: the build stage ``cadence`` and re-checked against this constant on
#: every run (a mismatch is logged loudly).
#:
#: HOW IT IS MEASURED, and why not the obvious way.  The obvious estimator
#: — the smallest (gap - EXPTIME) over all back-to-back StackPro pairs in
#: the archive — is an EXTREME ORDER STATISTIC over ~15,000 gaps, so a
#: single bad time stamp sets it.  It did: the three smallest "overheads"
#: in the whole archive (0.24 s, 0.64 s, 0.71 s) are all one frame landing
#: < 1.2 s from its neighbour inside a series whose own measured cadence
#: is 11-13 s.  A camera that demonstrably needs ~12 s per frame cannot
#: deliver two exposures 0.74 s apart; those pairs are time-stamp defects
#: (they are listed in ``s3_cadence_outliers``), not cadence.  The same
#: estimator returns NEGATIVE overheads (gap < EXPTIME, physically
#: impossible) in several plain High Gain cells — proof the raw minimum is
#: not a measurement of anything.
#:
#: What is used instead is robust by construction: for each same-config
#: run we take the MEDIAN gap (immune to single bad stamps), keep only
#: runs whose cadence is regular (:data:`CADENCE_REGULAR_SPREAD`), and
#: take the smallest (median gap - EXPTIME) over all such StackPro runs.
#: The physical statement is exact: a camera that repeatedly delivers a
#: frame every ``median gap`` seconds must fit exposure + internal dead
#: time + readout/save inside that cycle, so the internal dead time is
#: <= median gap - EXPTIME.  Threshold-insensitive: perturbing the
#: regularity cut moves this number by < 1 s, where the raw-minimum
#: estimator swings from 2.7 s to 8.9 s.
STACKPRO_DEADTIME_BOUND_S = 11.35

#: Worst-case error of the StackPro mid-time policy (start + EXPTIME/2),
#: in seconds: if the whole dead-time bound sat between sub-reads, the true
#: photon-weighted midpoint would shift late by at most half the bound.
#: NOTE this is SECONDS.  It is an honest ceiling from cadence alone, and
#: cadence is weak evidence: the gap it measures is dominated by readout
#: and file save, not by anything happening between sub-reads.  A camera
#: manual, or a lab measurement of the sub-read readout time, would
#: collapse it by orders of magnitude — until then no StackPro frame may
#: be quoted at sub-second absolute accuracy.
#:
#: What does NOT corroborate it: MaxIm's own JD-HELIO on a 1024 s StackPro
#: frame equals start + 512 s + heliocentric correction to 0.13 s.  That
#: only shows MaxIm ALSO assumes EXPTIME is the contiguous total span —
#: it constrains the software's convention, not the detector's physics.
STACKPRO_MID_WORST_CASE_S = STACKPRO_DEADTIME_BOUND_S / 2.0

# --------------------------------------------------------------------------
# Cadence-measurement policy (stage ``cadence``; the report interpolates
# these so the page always states the cuts it actually applied).
# --------------------------------------------------------------------------

#: Minimum frames for a (mode, target, night, exptime) run to count as a
#: series at all.
CADENCE_MIN_RUN = 5

#: Gaps at or beyond this multiple of a series' own median gap are pauses
#: (clouds, refocus, a slew), not cadence.
CADENCE_GAP_CEILING = 3.0

#: Gaps BELOW this fraction of a series' own median gap cannot be genuine
#: back-to-back pairs — the camera has just demonstrated, over the rest of
#: the same series, that it cannot cycle that fast.  They are discarded
#: from the overhead statistics and listed in ``s3_cadence_outliers``.
CADENCE_MIN_GAP_FRACTION = 0.5

#: A series counts as REGULAR (and so contributes its median gap to the
#: dead-time bound) when it has at least this many in-cadence gaps ...
CADENCE_REGULAR_MIN_GAPS = 4

#: ... and their interquartile spread is at most this fraction of the
#: median.  Both cuts exist so a ragged run (filter changes, dithers,
#: guiding losses) cannot masquerade as a machine cycle time.
CADENCE_REGULAR_SPREAD = 0.15

# --------------------------------------------------------------------------
# Clock-validation policy (stage ``clock``; ONE definition, imported by
# both the build script and the report so the figure and the fit cannot
# drift apart).
# --------------------------------------------------------------------------

#: |phase| beyond this is out of eclipse: these points set each
#: (readout config, filter) group's zero level.
CLOCK_OOE_PHASE = 0.08

#: Points inside this |phase| enter the dip fit.
CLOCK_FIT_PHASE = 0.12

#: A per-night eclipse fit is only a MEASUREMENT if the night's sampled
#: phases bracket the minimum.  A one-sided arc fits the eclipse's flank
#: (or the baseline's slope) and returns a confident, meaningless centre —
#: so a night must carry at least this many points strictly below phase 0
#: AND at least this many strictly above it, or its row is written with
#: status :data:`CLOCK_STATUS_ONE_SIDED` and a NULL O-C.
CLOCK_MIN_SIDE_POINTS = 3

#: Minimum points for a fit to be attempted at all (below this a row is
#: still written, with status :data:`CLOCK_STATUS_TOO_FEW`, so a
#: configured night can never vanish from the table without a trace).
CLOCK_MIN_FIT_POINTS = 8

#: Status strings written into ``s3_clock_eclipses.status``.
CLOCK_STATUS_OK = "ok"
CLOCK_STATUS_ONE_SIDED = "one_sided_coverage"
CLOCK_STATUS_TOO_FEW = "too_few_points"
CLOCK_STATUS_NO_DIP = "no_dip_found"

#: Light-travel time for one radian of angular offset at 1 au, seconds:
#: the small-angle conversion behind the frame-center caveat (an object
#: theta radians from the frame center differs from the stored
#: frame-center BJD by up to theta * this).  1 au = 499.00478 light
#: seconds (IAU 2012), and the barycentric baseline is ~1 au.
LIGHT_TIME_PER_RAD_S = 499.00478

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

#: ``frame_times.bjd_method`` values for the rows that have NO BJD, plus
#: the one that means "we had a BJD and withdrew it".
BJD_NO_JD = "no_jd"
BJD_NO_COORDS = "no_coords"
BJD_JD_DISAGREES = "jd_disagrees_with_raw_parent"

#: ``frame_times.start_evidence`` values.  Convention 1 (header JD is the
#: exposure START) is PROVEN per readout family by the header audit, and
#: the proof does not cover every family — this column says, on the row,
#: which case a consumer is looking at.
START_VERIFIED = "verified_jd_helio_or_telut"
START_UNVERIFIED = "unverified_no_helio_no_telut"

#: Above this |raw JD - reduced JD| disagreement (seconds) between two
#: copies of the SAME exposure (S0b ``raw_reduced_links``), the reduced
#: copy's stamp is not trustworthy and S3 withdraws its BJD.  Chosen to
#: sit far above the ~1 ms header write quantization and far below the
#: smallest real re-stamp (the observed disagreements are either < 1 s or
#: > 20 s — nothing lands near this line).
JD_SIBLING_DISAGREE_S = 10.0

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
      cadence of back-to-back StackPro series bounds any violation:
      total internal dead time <= :data:`STACKPRO_DEADTIME_BOUND_S`,
      hence worst-case mid-time error
      <= :data:`STACKPRO_MID_WORST_CASE_S` (late).  That ceiling is
      SECONDS wide — cadence cannot separate sub-read dead time from
      readout and save — so a StackPro frame is not a sub-second
      absolute time stamp until better evidence exists.

    EXPTIME <= 0 (a header pathology, never a real exposure) falls back to
    the start instant, labeled :data:`MID_EXPTIME_NONPOS`.
    """
    method = mid_method_for(readoutm, jd_start, exptime_s)
    if method == MID_NO_JD:
        return None, method
    if method == MID_EXPTIME_NONPOS:
        return float(jd_start), method
    return float(jd_start) + float(exptime_s) / 2.0 / 86400.0, method


def series_cadence(jds_sorted: Sequence[float], exptime_s: float,
                   min_gap_fraction: Optional[float] = None,
                   gap_ceiling: Optional[float] = None,
                   regular_min_gaps: Optional[int] = None,
                   regular_spread: Optional[float] = None) -> dict:
    """Robust cadence statistics for ONE same-config run of frames.

    ``jds_sorted`` are the exposure-start Julian Dates of a run of frames
    taken with one (readout mode, target, night, EXPTIME); they must be
    sorted ascending.  Returns a dict:

    * ``median_gap_s``  — median of the positive inter-start gaps.  A
      MEDIAN, not a minimum: one duplicated or out-of-order time stamp
      cannot move it, and it is the machine's actual cycle time.
    * ``kept_s``        — the gaps that are genuine back-to-back cadence:
      at least :data:`CADENCE_MIN_GAP_FRACTION` and less than
      :data:`CADENCE_GAP_CEILING` times the median.
    * ``short_idx``     — indices (into the gap array, i.e. gap i runs
      from frame i to frame i+1) of the gaps discarded as impossibly
      short.  These are time-stamp defects and get named individually in
      ``s3_cadence_outliers``.
    * ``regular``       — True when the kept gaps are numerous and tight
      enough (:data:`CADENCE_REGULAR_MIN_GAPS`,
      :data:`CADENCE_REGULAR_SPREAD`) that ``median_gap_s`` may be read
      as a machine cycle time.
    * ``spread``        — interquartile range of the kept gaps divided by
      the median (the regularity measure itself, so a caller can report
      how tight a series was rather than only the yes/no).
    * ``overhead_s``    — ``median_gap_s - exptime_s`` for a regular
      series, else None: the quantity that bounds everything the camera
      does between one exposure start and the next, internal sub-read
      dead time included.

    The four cut parameters default to the module constants; passing them
    explicitly exists so the build can SWEEP them and record how much the
    resulting bound actually moves, instead of the page asserting that it
    is insensitive.

    Pure numpy; no I/O.  Returns the same dict shape (with Nones) for a
    run too short to say anything about, so callers never branch on
    exceptions.
    """
    frac = (CADENCE_MIN_GAP_FRACTION if min_gap_fraction is None
            else min_gap_fraction)
    ceiling = CADENCE_GAP_CEILING if gap_ceiling is None else gap_ceiling
    min_gaps = (CADENCE_REGULAR_MIN_GAPS if regular_min_gaps is None
                else regular_min_gaps)
    spread_cut = (CADENCE_REGULAR_SPREAD if regular_spread is None
                  else regular_spread)
    jds = np.asarray(jds_sorted, dtype=float)
    empty = {"median_gap_s": None, "kept_s": np.array([]), "short_idx": [],
             "regular": False, "spread": None, "overhead_s": None,
             "n_gaps": 0}
    if jds.size < 2:
        return empty
    gaps = np.diff(jds) * 86400.0          # inter-start gaps, seconds
    positive = gaps > 0                    # a zero/negative gap is a dup
    if not positive.any():
        return empty
    median_gap = float(np.median(gaps[positive]))
    # Genuine cadence sits in a band around the series' own median: too
    # long = a pause between blocks, too short = an impossible cycle.
    too_short = positive & (gaps < frac * median_gap)
    in_band = positive & ~too_short & (gaps < ceiling * median_gap)
    kept = gaps[in_band]
    spread = None
    regular = False
    if kept.size >= min_gaps and median_gap > 0:
        spread = float(np.percentile(kept, 75) - np.percentile(kept, 25)) \
            / median_gap
        regular = spread <= spread_cut
    return {"median_gap_s": median_gap,
            "kept_s": kept,
            "short_idx": [int(i) for i in np.flatnonzero(too_short)],
            "regular": regular,
            "spread": spread,
            "overhead_s": (median_gap - float(exptime_s)) if regular else None,
            "n_gaps": int(in_band.sum())}


#: Plate scales outside this range cannot belong to any instrument that
#: could sit on a 3.45 m focal length telescope; a candidate outside it is
#: a broken header card, not a measurement.
PIXEL_SCALE_RANGE_ARCSEC = (0.01, 3.0)


def pixel_scale_arcsec(naxis1, naxis2, cdelt1=None, cdelt2=None,
                       cd1_1=None, cd1_2=None, xpixsz=None, focallen=None,
                       secpix1=None) -> tuple[Optional[float], str]:
    """Plate scale in arcsec/pixel from whatever cards a header carries,
    plus the NAME of the card it came from.

    Preference order, and why it is NOT "WCS first":

    1. ``SECPIX1`` — the configured plate scale, written by the
       acquisition software as a property of the instrument;
    2. ``206.265 * XPIXSZ[um] / FOCALLEN[mm]`` — the optics.  XPIXSZ is
       already the BINNED pixel in this archive (7.52 um at 2x2 for a
       3.76 um sensor), so no binning factor is applied.  FOCALLEN is not
       unit-safe here (the 2026 pyscope headers write 3.454, i.e. METRES,
       and some 2026 MaxIm headers write 0.0), which the range check
       below catches;
    3. the WCS CD matrix (``sqrt(CD1_1^2 + CD1_2^2)``, degrees/pixel);
    4. ``CDELT1`` (degrees/pixel).

    The WCS ranks LAST because a header WCS is not necessarily a solved
    one: a 0.009 s Mode0 frame in this archive carries a CD matrix
    claiming 3.08 arcsec/px on a camera whose optics say 0.449, with no
    PLTSOLVD card to back it — a leftover from a failed solve.  The
    instrument's own scale cannot be wrong in that way.

    Returns ``(None, 'unknown')`` when nothing usable is present or every
    candidate is outside :data:`PIXEL_SCALE_RANGE_ARCSEC`.
    """
    lo, hi = PIXEL_SCALE_RANGE_ARCSEC
    candidates: list[tuple[Optional[float], str]] = []
    if secpix1:
        candidates.append((abs(float(secpix1)), "SECPIX1"))
    if xpixsz and focallen:
        # 206.265 arcsec/rad x (um / mm) = 206.265 * um/mm arcsec/px.
        candidates.append((206.264806 * float(xpixsz) / float(focallen),
                           "XPIXSZ/FOCALLEN"))
    if cd1_1 is not None and cd1_2 is not None:
        candidates.append((math.hypot(float(cd1_1), float(cd1_2)) * 3600.0,
                           "CD matrix"))
    if cdelt1 is not None and float(cdelt1) != 0.0:
        candidates.append((abs(float(cdelt1)) * 3600.0, "CDELT1"))
    for scale, source in candidates:
        if scale is not None and lo <= scale <= hi:
            return float(scale), source
    return None, "unknown"


def field_corner_light_time_s(naxis1, naxis2, pixel_scale_arcsec_
                              ) -> Optional[float]:
    """Worst-case BJD difference (seconds) between a frame's CENTER and
    its most distant CORNER — the size of the caveat that travels with
    every ``frame_times`` row.

    An object theta radians away from the position the barycentric
    correction pointed at differs from the stored stamp by at most
    ``theta * LIGHT_TIME_PER_RAD_S`` (the full value is reached when the
    offset lies along the Earth-barycenter direction).  theta here is the
    frame's half-diagonal.  Returns None if the geometry is unknown.
    """
    if not naxis1 or not naxis2 or not pixel_scale_arcsec_:
        return None
    half_diag_arcsec = 0.5 * math.hypot(float(naxis1), float(naxis2)) \
        * float(pixel_scale_arcsec_)
    half_diag_rad = half_diag_arcsec / 206264.806
    return half_diag_rad * LIGHT_TIME_PER_RAD_S


def phase_coverage(phases, width: float = 0.0) -> tuple[int, int]:
    """Count fit points strictly below and strictly above phase zero.

    The pair a coverage gate needs: an eclipse mid-time is only measured
    when the sampling BRACKETS the minimum.  ``width`` optionally ignores
    a band around zero (points inside the dip itself constrain the depth,
    not the symmetry), so a caller can demand real shoulders on each
    side.  Pure; returns (n_before, n_after).
    """
    ph = np.asarray(phases, dtype=float)
    return int(np.sum(ph < -abs(width))), int(np.sum(ph > abs(width)))


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
