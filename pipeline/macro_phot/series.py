"""Pure production-scale series logic for the CV time-series campaign.

The S4 prototype (``build_s4_photometry.py``) proved the photometry core on
two polars whose pixels were ALL server-reduced and whose star matching was
ALL astroalign.  Running the same core over the five staged CV targets adds
four questions the prototype never had to answer, and every one of them is
a RULE rather than a computation.  Those rules live here, as pure functions
on plain values, so a student auditor can read them without a telescope and
``pipeline/tests/test_series.py`` can test them without a pixel:

1.  **Which pixels may a series use?**  ``stage_cv_timeseries`` stages RAW
    frames.  Some eras have a complete server-reduced counterpart, some
    have none at all, and a few are partial.  The provenance rule
    (:func:`choose_provenance`) picks ONE provenance for a whole
    (target, era) — never a mixture inside one series — and names what it
    dropped.
2.  **Where is the saturation ceiling?**  The S2 campaign measured a
    DIFFERENT clip for every readout mode (High Gain digitizes 12 bits and
    clips at 3,496 ADU; Mode0/Fast fill 16 bits; StackPro sums 16
    sub-reads and clips near 56,062).  A single global threshold — the
    prototype's 55,000 ADU — vetoes nothing at all on High Gain, which is
    exactly the era where saturation bites hardest.  :func:`veto_adu` and
    :func:`veto_in_reduced_adu` carry the measured per-mode numbers, the
    latter mapping the RAW threshold through a server reduction whose
    (dark, flat, pedestal) the S2 reconstruction measured.
3.  **How is a frame registered?**  About half the staged frames now carry
    a plate solution from S1.  Where BOTH a frame and its series reference
    are solved, the star match can go through the sky instead of through
    astroalign's triangles — faster and immune to triangle starvation —
    but only if the resulting match is credible (:func:`wcs_match_ok`),
    otherwise the astroalign ladder still runs.
4.  **Is the series photometrable at all?**  An image can be too small to
    hold one aperture plus its sky annulus, and :func:`geometry_verdict`
    says so once, in one place, instead of leaving mysterious extraction
    failures in a log.  It is judged against RESOLVED pixel dimensions
    only.  RETIRED PREMISE, recorded because the campaign nearly acted on
    it: EU UMa's era-80 frames were believed to be 8-pixel-wide readout
    STRIPS narrower than a 4-arcsec aperture, and all 207 were excluded on
    that basis.  They are nothing of the kind — the files hold 4,800 x
    3,211 raw / 4,787 x 3,193 reduced images, and the 8 is the BINTABLE
    ROW LENGTH of a tile-compressed header that the S0 scan copied without
    translating.  ``macro_core.fitsgeom`` exists to translate it, and
    geometry is now judged only against dimensions that resolver returns.

Nothing in this module does I/O, and nothing here imports anything heavier
than numpy — the same house standard as ``macro_phot.photometry``.
"""

from __future__ import annotations

import math
from typing import Iterable, Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Series identity
# --------------------------------------------------------------------------

#: Field separator inside a series key.  Chosen because no target key, era
#: label, or filter name in the archive contains it, so the key round-trips.
SERIES_SEP = "|"


def series_key(target_key: str, era_id: int, filt: Optional[str]) -> str:
    """The canonical name of one (target, era, filter) series.

    Every photometry row in the products database carries this string, so
    a reader can never accidentally pool two eras (different camera, plate
    scale and zero point) or two filters into one light curve.  A missing
    FILTER card becomes the literal ``'none'`` rather than SQL NULL: NULL
    silently loses rows to every ``=`` comparison, and an unfiltered frame
    is a real observing state that deserves a visible name.
    """
    f = "none" if filt is None or str(filt).strip() == "" else str(filt)
    return f"{target_key}{SERIES_SEP}e{int(era_id)}{SERIES_SEP}{f}"


def parse_series_key(key: str) -> tuple[str, int, str]:
    """Inverse of :func:`series_key` — ``(target_key, era_id, filter)``."""
    tk, era, filt = key.split(SERIES_SEP)
    return tk, int(era.lstrip("e")), filt


# --------------------------------------------------------------------------
# S2 detector facts: the per-mode saturation ceiling and veto
# --------------------------------------------------------------------------

#: Adopted per-mode ceilings and saturation vetoes in RAW ADU, MEASURED by
#: the S2 detector campaign and published in ``docs/pipeline/s2_detector.html``
#: (the campaign's own tables; the veto is
#: ``floor(0.92 * ceiling / 100) * 100`` — see
#: ``rlmt_diagnostics.ceiling.veto_threshold``).  Copied here as the S4
#: consumer's single source of truth because the S2 tables were carried in
#: the manifest database that S0's rebuild replaced; the report remains the
#: evidence.  Modes S2 sampled but never saw saturate ('Low Gain',
#: '5MHz High Sensitivity 16-bit') are ABSENT on purpose: an unmeasured
#: ceiling must read as "unknown", not as a fabricated 65,535.
S2_MODE_CEILING_ADU: dict[str, int] = {
    "Mode0": 65535,
    "Fast": 65534,
    "(blank 2026)": 65534,
    "1MHz High Sensitivity 16-bit": 64674,
    "High Gain": 3496,
    "High Gain StackPro": 56062,
}

S2_MODE_VETO_ADU: dict[str, int] = {
    "Mode0": 60200,
    "Fast": 60200,
    "(blank 2026)": 60200,
    "1MHz High Sensitivity 16-bit": 59500,
    "High Gain": 3200,
    "High Gain StackPro": 51500,
}

#: Nominal gain used for the photon-shot term of every flux error, e-/ADU.
#: S2 bracketed the true gain at [0.60, 1.77] e-/ADU and the header card
#: EGAIN reads 1.057 on the CMOS eras; the header value is demonstrably
#: unusable on others (Mode0 writes 0.247, the iKon writes 0 or nothing),
#: so ONE nominal value is applied everywhere and the bracket is recorded
#: beside it.  The consequence is honest and bounded: a gain wrong by the
#: full bracket scales the predicted photon error by at most sqrt(1.77 /
#: 0.60) = 1.7, and the empirical inflation factor measured from held-out
#: check stars (macro_phot.errors) absorbs exactly that kind of constant
#: mis-scaling.  Never quote the nominal as if it were measured.
NOMINAL_GAIN_E_PER_ADU = 1.057
GAIN_BRACKET_E_PER_ADU = (0.60, 1.77)


def veto_adu(readoutm: Optional[str]) -> Optional[int]:
    """The S2 saturation veto for one readout mode, or None if unmeasured.

    None is a real answer: it means S2 never saw this mode saturate, so
    the caller must record "no veto applied" instead of inventing one.
    """
    if readoutm is None:
        return None
    return S2_MODE_VETO_ADU.get(str(readoutm).strip())


def veto_in_reduced_adu(veto_raw: Optional[float], flat_med: float,
                        dark_med: float, pedestal: float) -> Optional[float]:
    """Map a RAW-ADU saturation veto into SERVER-REDUCED ADU.

    The S2 reconstruction experiment fitted the server pipeline as

        reduced = (raw - D) / F + pedestal

    per era (``products/detector/recon/eraNN.npz``: ``D`` the dark image,
    ``F`` the flat image, ``pedestal`` the additive constant).  A star that
    saturated the DETECTOR at ``veto_raw`` therefore appears in the reduced
    frame at this level, so the veto travels with the pixels instead of
    being re-guessed.  Medians of D and F stand in for the images: the
    pixel-to-pixel spread of a flat is a few percent, far below the 8%
    margin the 0.92 veto fraction already leaves under the hard clip.

    Returns None when the mode has no measured veto (None in, None out).
    """
    if veto_raw is None:
        return None
    if not (flat_med and flat_med > 0):
        return None
    return (float(veto_raw) - float(dark_med)) / float(flat_med) \
        + float(pedestal)


def applied_veto_adu(veto_adu_raw: Optional[float],
                     dark_median_adu: Optional[float]) -> Optional[float]:
    """The veto level that actually applies to the pixels being measured.

    A locally-calibrated frame has had its master dark SUBTRACTED before
    anything is measured, so a threshold defined on raw pixels sits that
    much too high on the calibrated ones.  (The flat is median-normalised
    and shifts no levels, which is why it does not appear here.)  Server-
    reduced frames pass ``dark_median_adu=None``: their veto was already
    mapped through the reduction by :func:`veto_in_reduced_adu`.

    This one line used to live inline in the extraction worker, and the
    number it produced was never written down — the database recorded the
    RAW veto while the saturation flags had been computed with this one.
    9,377 era-7 detections could not be reproduced from the product as a
    result.  A derived threshold that decides a published statistic has to
    be a named function with a column of its own.
    """
    if veto_adu_raw is None:
        return None
    if dark_median_adu is None:
        return float(veto_adu_raw)
    return float(veto_adu_raw) - float(dark_median_adu)


#: A recorded image dimension below this many pixels cannot be a real
#: science image from any camera in this archive (the smallest is 1,024 px
#: on its short axis), so it is a header artifact — almost always the
#: BINTABLE row length of a tile-compressed frame, which is 8 on EU UMa's
#: era 80.  Frames flagged here get their geometry re-resolved from the
#: file rather than trusted, because the alternative is what already
#: happened once: 207 good frames excluded as "8-pixel readout strips".
MIN_PLAUSIBLE_AXIS_PX = 64


def geometry_is_implausible(naxis1: Optional[float],
                            naxis2: Optional[float],
                            min_axis: int = MIN_PLAUSIBLE_AXIS_PX) -> bool:
    """True when a recorded (naxis1, naxis2) cannot describe a real image.

    Deliberately NOT a photometry decision — :func:`geometry_verdict` makes
    those.  This is a provenance test on the NUMBER: it asks whether the
    dimension is worth believing at all, so a caller knows when to go back
    to the pixels.  Missing values count as implausible, because "unknown"
    and "8" both mean the same thing here: go and look.
    """
    if naxis1 is None or naxis2 is None:
        return True
    try:
        short = min(float(naxis1), float(naxis2))
    except (TypeError, ValueError):
        return True
    return not math.isfinite(short) or short < float(min_axis)


def saturated_mask(peak_adu: np.ndarray, background_adu: float,
                   veto: Optional[float]) -> np.ndarray:
    """Which detections sit at or above the saturation veto.

    sep reports ``peak`` on the BACKGROUND-SUBTRACTED image, while the S2
    veto is a level on the image as digitized — so the background is added
    back before the comparison.  With no veto for this mode (None) nothing
    is flagged, and the caller records that no veto was available: an
    unmeasured ceiling must not silently become "nothing is saturated".
    """
    peak = np.asarray(peak_adu, dtype=float)
    if veto is None:
        return np.zeros(peak.shape, dtype=bool)
    return (peak + float(background_adu)) >= float(veto)


# --------------------------------------------------------------------------
# Provenance: which pixels a (target, era) series is allowed to use
# --------------------------------------------------------------------------

#: A (target, era) uses server-reduced pixels when at least this fraction of
#: its staged science frames HAVE a reduced counterpart.  Between the
#: threshold and 100% the unlinked minority is DROPPED, not raw-substituted:
#: the provenance rule forbids two reduction histories inside one series,
#: and a dropped frame is counted in the selection ledger where a reader
#: can see it.
PROVENANCE_MIN_LINK_FRAC = 0.5


def choose_provenance(n_staged: int, n_linked: int,
                      has_master_calib: bool,
                      min_frac: float = PROVENANCE_MIN_LINK_FRAC
                      ) -> tuple[str, str]:
    """Pick ONE pixel provenance for a whole (target, era).

    Returns ``(provenance, reason)`` where provenance is one of:

    ``'server_reduced'``
        The observatory's own reduced tree.  Used when the reduced
        counterpart covers at least ``min_frac`` of the staged frames;
        frames without one are excluded from the series.
    ``'local_master'``
        No usable reduced tree, but era-matched MASTER calibrations are
        staged: this build subtracts the master dark and divides by the
        master flat itself, uniformly across the era.  One recipe, one
        provenance, applied to every frame of the series.
    ``'raw'``
        Neither.  Photometry proceeds on raw pixels and the series is
        flagged: a raw frame still carries dark current and hot pixels
        that a local sky annulus only partly hides, and its detection list
        is measurably contaminated (a test pair showed 281 raw detections
        collapse to 161 real stars once dark-subtracted).

    The reason string is stored verbatim in the products database so the
    decision is auditable without re-deriving it.
    """
    if n_staged <= 0:
        return "none", "no staged frames"
    frac = n_linked / float(n_staged)
    if frac >= min_frac:
        dropped = n_staged - n_linked
        return "server_reduced", (
            f"{n_linked}/{n_staged} staged frames have a reduced counterpart "
            f"({frac:.0%} >= {min_frac:.0%}); {dropped} unlinked frame(s) "
            f"excluded rather than mixed with raw pixels")
    if has_master_calib:
        return "local_master", (
            f"only {n_linked}/{n_staged} staged frames are reduced "
            f"({frac:.0%} < {min_frac:.0%}); era-matched master "
            f"calibrations applied locally to every frame instead")
    return "raw", (
        f"only {n_linked}/{n_staged} staged frames are reduced "
        f"({frac:.0%}) and no era-matched master calibration is staged; "
        f"raw pixels used and the series flagged")


# --------------------------------------------------------------------------
# Plate scale: header cards first, the frame's own sky transform second
# --------------------------------------------------------------------------

def plate_scale_from_cd(cd1_1, cd1_2, cd2_1, cd2_2) -> Optional[float]:
    """Plate scale in arcsec/pixel from a FITS CD matrix.

    The CD matrix maps pixel offsets to degrees on the sky, so the area of
    one pixel on the sky is |det CD| square degrees and the linear scale is
    its square root.  Taking the determinant rather than |CD1_1| is what
    makes this correct for a ROTATED frame, where CD1_1 alone shrinks with
    the cosine of the rotation angle and would understate the scale — by up
    to 100% at 90 degrees.  Missing off-diagonal terms are treated as zero
    (the common CDELT-style header).  Returns None when the matrix is
    absent or degenerate.
    """
    if cd1_1 is None and cd2_2 is None:
        return None
    a = float(cd1_1 or 0.0)
    b = float(cd1_2 or 0.0)
    c = float(cd2_1 or 0.0)
    d = float(cd2_2 or 0.0)
    det = abs(a * d - b * c)
    if not (det > 0):
        return None
    return math.sqrt(det) * 3600.0


def resolve_plate_scale(meta: dict) -> tuple[Optional[float], str]:
    """This frame's plate scale and the EVIDENCE it rests on.

    Order of preference, with the basis string recorded per frame:

    1. ``'header_optics'`` — XPIXSZ (binned pixel size, um) and FOCALLEN
       (mm).  The primary route: it is independent of any astrometric
       solution and available on nearly every frame.
    2. ``'header_cd'`` — the frame's own CD matrix, i.e. a scale MEASURED
       against the sky rather than computed from optics.  This exists for
       exactly the reason a fallback must: the March-2026 'Fast' frames
       write ``FOCALLEN = 0.0``, which is not a small focal length but a
       missing card, and rejecting them for 'no plate scale' discarded a
       74-frame EU UMa observing block that was carrying a perfectly good
       0.4509"/px CD matrix all along.
    3. ``'header_cdelt'`` — the older CDELT1/CDELT2 form of the same thing.

    Returns ``(scale_or_None, basis)``; the basis is ``'none'`` when no
    route worked, so the failure is named rather than inferred.
    """
    from . import photometry as _ph          # local: avoids a cycle at import
    s = _ph.plate_scale_arcsec_per_px(meta.get("xpixsz"), meta.get("focallen"))
    if s is not None:
        return s, "header_optics"
    s = plate_scale_from_cd(meta.get("cd1_1"), meta.get("cd1_2"),
                            meta.get("cd2_1"), meta.get("cd2_2"))
    if s is not None:
        return s, "header_cd"
    d1, d2 = meta.get("cdelt1"), meta.get("cdelt2")
    if d1 and d2:
        s = math.sqrt(abs(float(d1) * float(d2))) * 3600.0
        if s > 0:
            return s, "header_cdelt"
    return None, "none"


# --------------------------------------------------------------------------
# Geometry: can this frame be aperture-photometered at all?
# --------------------------------------------------------------------------

#: The shorter image axis must span at least this many aperture RADII for
#: aperture photometry to mean anything: the sky annulus alone reaches
#: SKY_ANNULUS_ARCSEC[1] / APERTURE_RADIUS_ARCSEC = 3 radii from a star's
#: centre, so an image narrower than 6 radii cannot hold one complete
#: star-plus-sky measurement across its short axis.
MIN_AXIS_APERTURE_RADII = 6.0


def geometry_verdict(naxis1: Optional[float], naxis2: Optional[float],
                     aper_px: Optional[float]) -> tuple[bool, str]:
    """Is a frame of this shape photometrable with this aperture?

    Returns ``(ok, reason)``.  The rule guards one thing only: an image
    whose SHORT axis cannot hold a photometric aperture plus its sky
    annulus, where every "measurement" would really be the image's own
    edge.  ``naxis1``/``naxis2`` must be RESOLVED dimensions — what
    ``macro_core.fitsgeom`` returns for the frame's real pixels — never a
    raw NAXIS card copied out of a tile-compressed header, where NAXIS1 is
    the BINTABLE row length in bytes and has nothing to do with the image.
    Feeding this function that number is how the campaign came within one
    commit of throwing away 207 perfectly good EU UMa frames.
    """
    if naxis1 is None or naxis2 is None:
        return False, "image dimensions unknown"
    if aper_px is None or aper_px <= 0:
        return False, "no plate scale, so no aperture in pixels"
    short = min(float(naxis1), float(naxis2))
    need = MIN_AXIS_APERTURE_RADII * float(aper_px)
    if short < need:
        return False, (
            f"short axis {short:.0f} px < {MIN_AXIS_APERTURE_RADII:g} x "
            f"aperture radius {aper_px:.1f} px = {need:.0f} px: the "
            f"photometric aperture does not fit inside the frame")
    return True, "aperture and sky annulus fit inside the frame"


# --------------------------------------------------------------------------
# Registration: sky-chained (WCS) vs triangle-matched (astroalign)
# --------------------------------------------------------------------------

#: A WCS-chained match must recover at least this fraction of the stars the
#: smaller of (frame detections, reference stars) could possibly share, and
#: at least this many stars outright, before it is believed.  A plate
#: solution can be right about the sky and still be useless here — solved
#: on a different pixel grid, or with SIP terms fitted off the frame edge —
#: and a silently bad chain would hand the ensemble a scrambled star
#: identity map.  Failing the test is not fatal: the astroalign ladder runs
#: instead and the fallback is recorded per frame.
WCS_MIN_MATCH_FRAC = 0.30
WCS_MIN_MATCH_ABS = 8


def registration_method(frame_has_wcs: bool, ref_has_wcs: bool) -> str:
    """The registration route this frame is entitled to try first.

    Chaining through the sky needs BOTH ends solved — the frame to leave
    its pixels and the reference to receive them.  Everything else goes to
    astroalign, which needs no astrometry at all.
    """
    return "wcs" if (frame_has_wcs and ref_has_wcs) else "astroalign"


def wcs_match_ok(n_match: int, n_det: int, n_ref: int,
                 min_frac: float = WCS_MIN_MATCH_FRAC,
                 min_abs: int = WCS_MIN_MATCH_ABS) -> bool:
    """Is a WCS-chained star match credible enough to keep?

    The denominator is ``min(n_det, n_ref)`` — the most stars the two
    catalogs could share — so a shallow frame against a deep reference is
    judged on what it could have matched, not on the reference's depth.
    """
    if n_match < min_abs:
        return False
    possible = min(int(n_det), int(n_ref))
    if possible <= 0:
        return False
    return (n_match / float(possible)) >= min_frac


#: How deep a PLATE-SOLVED reference candidate must be, as a fraction of
#: the deepest candidate's detection count, before its plate solution is
#: worth preferring.  A solved reference buys three things — the cheap
#: sky-chained route for every solved frame in the era, a target
#: identification straight from catalogue coordinates with no network
#: query, and a real sky position for every star in the catalog — but the
#: reference's own star list is the ceiling on what any frame can match,
#: so a solved-but-shallow frame would quietly cap the whole series.  0.6
#: buys the astrometry whenever it is nearly free and refuses to pay for
#: it with a third of the field.
REF_SOLVED_MIN_DEPTH_FRAC = 0.6


def order_reference_candidates(ranking: Sequence[int],
                               n_detected: dict, solved: dict,
                               min_frac: float = REF_SOLVED_MIN_DEPTH_FRAC
                               ) -> list[int]:
    """Re-order reference candidates to prefer a plate-solved one — if it is
    deep enough to deserve the promotion.

    ``ranking`` is the quality order from
    :func:`macro_phot.photometry.rank_references` (best first).  Solved
    candidates whose detection count reaches ``min_frac`` of the deepest
    candidate's move to the front, keeping their relative order; every
    other candidate follows in unchanged rank order, so nothing is ever
    dropped — the caller still walks the whole list until one passes the
    double-image quality control.
    """
    cand = list(ranking)
    if not cand:
        return []
    deepest = max(float(n_detected.get(c) or 0) for c in cand)
    floor = min_frac * deepest
    promoted = [c for c in cand
                if solved.get(c) and float(n_detected.get(c) or 0) >= floor]
    rest = [c for c in cand if c not in set(promoted)]
    return promoted + rest


#: How many times the CHANCE-coincidence expectation a star match must
#: exceed before a proposed registration is believed, and the floor on the
#: absolute count.  See :func:`matches_beat_chance`.
CHANCE_MATCH_FACTOR = 10.0
CHANCE_MATCH_MIN_ABS = 8


def expected_chance_matches(n_det: int, n_ref: int, tol_px: float,
                            area_px2: float) -> float:
    """How many star pairs a WRONG registration would match by luck.

    Scatter ``n_ref`` reference stars at random over an image of
    ``area_px2`` pixels; the chance that a given detection lands within
    ``tol_px`` of one of them is ``n_ref * pi * tol^2 / area``.  Multiply by
    ``n_det`` detections for the expected count.  Concretely, for a 4,787 x
    3,193 EU UMa frame with 1,839 reference stars and a 2-pixel tolerance,
    a completely wrong shift is expected to match 0.6 pairs.
    """
    if area_px2 <= 0 or tol_px <= 0:
        return float("inf")
    per_det = float(n_ref) * math.pi * float(tol_px) ** 2 / float(area_px2)
    return float(n_det) * per_det


def matches_beat_chance(n_match: int, n_det: int, n_ref: int, tol_px: float,
                        area_px2: float,
                        factor: float = CHANCE_MATCH_FACTOR,
                        min_abs: int = CHANCE_MATCH_MIN_ABS) -> bool:
    """Is a proposed registration supported by more matches than luck allows?

    The fraction-based gate (:func:`wcs_match_ok`) asks a frame to match a
    THIRD of everything it detected, which is the right question for a deep
    frame and the wrong one for a cloudy one: a 240-second EU UMa exposure
    through cirrus reports 400 detections of which perhaps 60 are stars, so
    demanding 120 matches rejects a registration that is in fact certain.
    This gate asks the question that actually matters — could luck have
    produced this? — and for a wrong shift the answer is a fraction of one
    pair, so even a dozen real matches clears it by an order of magnitude.
    The absolute floor keeps a nearly-empty frame from qualifying on two
    coincidences in a sparse field.
    """
    if n_match < min_abs:
        return False
    chance = expected_chance_matches(n_det, n_ref, tol_px, area_px2)
    return float(n_match) >= factor * max(chance, 0.1)


def match_rate(n_match: int, n_det: int, n_ref: int) -> float:
    """Matched fraction on the same denominator :func:`wcs_match_ok` uses.

    Reported per frame and averaged per series as the campaign's
    'match rate' diagnostic.  Returns NaN when nothing could match, which
    a mean over frames must therefore skip rather than count as zero.
    """
    possible = min(int(n_det), int(n_ref))
    if possible <= 0:
        return float("nan")
    return float(n_match) / float(possible)


# --------------------------------------------------------------------------
# Master-calibration selection
# --------------------------------------------------------------------------

#: A master dark serves a science frame only when their exposure times
#: agree within this relative tolerance (with the absolute floor below).
#: The archive's masters sit on the power-of-two ladder the acquisition
#: software uses, and eras 6/7 science frames sit on exactly the same
#: ladder (8/16/32/64/128 s), so exact matches are the norm.  Scaling a
#: mismatched dark is REFUSED rather than approximated: these masters
#: include the bias, and scaling a bias-inclusive dark by an exposure ratio
#: corrupts the bias along with the dark current.
DARK_EXPTIME_REL_TOL = 0.005
DARK_EXPTIME_ABS_TOL = 0.02


def dark_exptime_matches(dark_exptime: Optional[float],
                         sci_exptime: Optional[float]) -> bool:
    """Does this master dark's exposure time serve this science frame?"""
    if dark_exptime is None or sci_exptime is None:
        return False
    d, s = float(dark_exptime), float(sci_exptime)
    if math.isnan(d) or math.isnan(s):
        return False
    return abs(d - s) <= max(DARK_EXPTIME_ABS_TOL, DARK_EXPTIME_REL_TOL * s)


def pick_master(candidates: Sequence[tuple], sci_jd: Optional[float]
                ) -> Optional[tuple]:
    """Choose one master from several equally eligible ones.

    ``candidates`` are ``(jd, path, ...)`` tuples already filtered for
    kind, era, filter and (for darks) exposure time.  The winner is the
    master nearest in time to the science frame, ties broken by path
    string — a total order, so the same frame calibrates identically on
    every re-run.  A candidate with no JD sorts last but remains usable
    (several archived masters carry no date at all).
    """
    if not candidates:
        return None
    def key(c):
        jd = c[0]
        far = float("inf") if (jd is None or sci_jd is None) else abs(
            float(jd) - float(sci_jd))
        return (far, str(c[1]))
    return sorted(candidates, key=key)[0]


# --------------------------------------------------------------------------
# Series-level admission and diagnostics
# --------------------------------------------------------------------------

#: A series needs at least this many successfully matched frames before an
#: ensemble solve is attempted.  Below it the Honeycutt solve is not wrong,
#: it is empty: the comparison-star stability iteration needs stars seen on
#: at least COMP_MIN_FRAME_FRAC (half) of the frames, and "half of six" is
#: not a statistic anybody should publish a zero point from.
MIN_SERIES_FRAMES = 12

#: The campaign's contract with the next analyst: every solved series must
#: hold out at least this many CHECK stars, which never vote in the zero
#: point and therefore give an honest error estimate.  Series that cannot
#: are solved anyway and flagged — the measurement still exists, but its
#: error bars are unvalidated and must be labelled so.
MIN_CHECK_STARS = 3


def series_admission(n_matched_frames: int,
                     min_frames: int = MIN_SERIES_FRAMES
                     ) -> tuple[bool, str]:
    """May this series be handed to the ensemble solver? With the reason."""
    if n_matched_frames <= 0:
        return False, "no frames survived extraction and matching"
    if n_matched_frames < min_frames:
        return False, (
            f"only {n_matched_frames} matched frame(s) < {min_frames}: too "
            f"few for a comparison-star stability iteration")
    return True, f"{n_matched_frames} matched frames"


#: A target with fewer than this many finite magnitudes is not a light
#: curve.  12 matches MIN_SERIES_FRAMES: a series admitted on 12 frames
#: whose target was measured on all of them is the smallest thing this
#: campaign is willing to call a measurement.
MIN_TARGET_POINTS = 12

#: Fractions of the APPLIED saturation veto at which a target measurement
#: is counted as "near the ceiling".  A hard clip is not where a detector
#: stops being linear — response rolls over below it — so a binary
#: saturated/not flag can read 0 while the brightest points of the series
#: are being quietly compressed.  0.90 and 0.85 bracket the region where
#: that matters; both are reported so a reader can see the shape of the
#: approach rather than one threshold's verdict.
NEAR_VETO_FRACTIONS = (0.85, 0.90)


def target_series_verdict(n_rows: int, n_measured: int,
                          min_points: int = MIN_TARGET_POINTS) -> str:
    """One word on whether this series delivers a target light curve.

    ``n_rows`` is how many light-curve rows the target has; ``n_measured``
    how many of them carry a finite magnitude.  The two differ whenever a
    frame got no ensemble zero point or the target's flux was non-positive,
    and the difference is not cosmetic: ``stlmi|e47|y`` was recorded as
    ``solved`` with a ``validated`` check verdict and 10 "target points",
    every one of which was NULL — the series has no y-band light curve at
    all, while every summary of the run presented one.

    Returns ``'measured'`` (a usable target light curve), ``'sparse'``
    (some magnitudes, too few to be one), or ``'undetected'`` (no finite
    magnitude anywhere, or no target row at all).
    """
    if n_rows <= 0 or n_measured <= 0:
        return "undetected"
    return "measured" if n_measured >= min_points else "sparse"


def near_veto_counts(level_adu, veto_applied: Optional[float],
                     fractions: Sequence[float] = NEAR_VETO_FRACTIONS
                     ) -> dict[float, int]:
    """How many measurements sit within a stated fraction of the ceiling.

    ``level_adu`` is each measurement's peak PLUS its frame background —
    the same quantity :func:`saturated_mask` compares, i.e. the level on
    the image as digitized — and ``veto_applied`` is the threshold that was
    actually used (see :func:`applied_veto_adu`), never the raw one.
    Returns ``{fraction: count}``, empty when there is no veto to measure
    against, because "no ceiling was measured for this mode" must not
    render as "nothing was near the ceiling".
    """
    if veto_applied is None or not math.isfinite(float(veto_applied)) \
            or float(veto_applied) <= 0:
        return {}
    lv = np.asarray(level_adu, dtype=float)
    lv = lv[np.isfinite(lv)]
    return {float(f): int((lv >= float(f) * float(veto_applied)).sum())
            for f in fractions}


def check_star_verdict(n_check: int, min_check: int = MIN_CHECK_STARS
                       ) -> str:
    """One word on whether the held-out check stars validate this series."""
    if n_check >= min_check:
        return "validated"
    if n_check > 0:
        return "weak"
    return "unvalidated"


def mean_ignoring_nan(values: Iterable[float]) -> float:
    """Mean of the finite entries; NaN when there are none.

    Used for the per-series match rate, whose per-frame value is NaN for
    any frame that could not have matched anything.  Written out rather
    than reaching for ``np.nanmean`` so an all-NaN series answers NaN
    quietly instead of emitting a RuntimeWarning per empty series.
    """
    v = np.asarray(list(values), dtype=float)
    v = v[np.isfinite(v)]
    return float(v.mean()) if v.size else float("nan")
