"""macro_sn.gate0 — the decision rules of SN 2023ixf Gate 0.

Pure functions only.  Nothing here opens a file, a database or a socket; the
plumbing lives in ``pipeline/scripts/run_sn_gate0.py`` and the rendering in
``macro_sn.report_gate0``.  That split is what makes every rule below
testable without an archive drive attached, and it is why the thresholds are
module constants rather than literals buried in a query.

WHAT GATE 0 IS
--------------
``SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md`` §4 Step 0 names three blocking
activities and forbids all downstream work until they land:

* **0a  manifest freeze**   — one row per unique frame, deduplicated
  GLOBALLY (including within ``rawimage/``, whose July directories hold
  wholesale copies of earlier nights) and with the target's name aliases
  merged, so that a frame labelled ``NGC5457`` on the discovery night is
  counted once and counted at all.
* **0b  saturation census** — the peak ADU of the supernova ITSELF measured
  from the pixels in every direct-imaging frame of the campaign, screened
  against the detector's measured ceiling, and published as a filter × night
  matrix.  This decides the true clean start per band.
* **0c  grism triage**      — whether the slitless series is a usable
  spectral record, judged against pre-registered promotion criteria.

THE ONE RULE THAT CHANGED
-------------------------
The strategy was written believing filter slot ``6`` was a grism because of
its LABEL.  Stage S2c has since measured dispersion frame by frame, and the
slot is MIXED: on this target 61 of its 83 frames carry measured traces, 3
are measured direct images, and 19 could not be certified either way.  The
same label-vs-measurement defect had contaminated the S1 astrometry
experiment, whose published rate and failure taxonomy for this campaign were
computed over a denominator containing spectra.

So every gate in this module reads the MEASUREMENT
(``frame_dispersion.verdict``) and never a filter string.  The rule is
stated once, in :func:`is_direct_image`, and is deliberately identical to
the one S1 v1.2 adopted in ``macro_core.astrom.is_measured_spectrum``:
exclusion must be EARNED by a positive measurement.  A frame is removed from
the imaging census only when S2c measured it dispersed; ``direct``,
``indeterminate`` and never-measured frames all stay in, because "we did not
look" is not evidence of a spectrum.

WHY A BOX MAXIMUM IS A LEGITIMATE ANSWER
----------------------------------------
Two thirds of the campaign's frames carry a plate solution in their header,
and for those the supernova's pixel is known exactly.  The rest do not, and
the commanded pointing misses the true field centre by ~200 px.  Rather than
guess, the census reports what it can actually defend: over a search box
large enough to contain the supernova wherever it fell, the box MAXIMUM is
an upper bound on the supernova's peak.  A frame whose entire box sits under
the screen is therefore provably unsaturated — a one-sided but genuine
measurement — while a frame whose box exceeds it is recorded as
``undetermined`` rather than as a saturated epoch it may not be.  See
:func:`saturation_class`.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

# The astrometry module owns two things this stage must not re-implement:
# the dispersion-class vocabulary (so "dispersed" means the same string in
# both stages) and the GO / CAUTION / NO-GO arithmetic (so the Gate 0
# astrometry answer is comparable, threshold for threshold, with the S1
# stratum verdict it is being asked about).  Re-deriving either here would
# create two definitions that could drift apart silently, which is the exact
# failure this whole correction exists to repair.
from macro_core import astrom


# ===========================================================================
# 1.  CONSTANTS — every one of them carries its provenance in the comment
# ===========================================================================

#: The supernova's J2000 position, in degrees.
#:
#: PROVENANCE: the discovery position of SN 2023ixf reported by Itagaki and
#: circulated as TNS AT 2023ixf — RA 14:03:38.562, Dec +54:18:41.94 (J2000).
#: It is a LITERATURE CONSTANT, not a measurement made in this repository,
#: and it is declared here rather than typed into a query so that the page,
#: the census and the tests all read the same three numbers.  The census
#: validates it in passing: in every plate-solved frame the brightest pixel
#: within a 21 px stamp of this position lands within ~1 px of it, which a
#: wrong position could not do.
SN_RA_DEG = 210.910675
SN_DEC_DEG = 54.311650
SN_POSITION_SOURCE = ("Itagaki discovery position, TNS AT 2023ixf "
                      "(14:03:38.562 +54:18:41.94, J2000)")

#: Adopted explosion epoch, MJD.  PROVENANCE: Hosseinzadeh+23 (ApJL 953 L16)
#: and Li+24 (Nature 627 754) bracket first light at MJD 60082.75–60082.83;
#: the strategy adopts that range and this module takes its midpoint.
#:
#: PHASE IS DISPLAY-ONLY.  No Gate 0 verdict is a function of phase: the
#: census screens on measured ADU and the triage screens on measured
#: dispersion.  Phase appears on the page so a reader can see WHEN a night
#: sits, and its uncertainty (±0.04 d) is smaller than the nightly cadence,
#: so no row could move between nights if the midpoint were wrong.
T0_MJD = 60082.79
T0_MJD_LOW = 60082.75
T0_MJD_HIGH = 60082.83
T0_SOURCE = "Hosseinzadeh+23 (ApJL 953 L16); Li+24 (Nature 627 754)"

#: MJD = JD - this.  The IAU definition, not a fitted quantity.
MJD_OFFSET = 2400000.5

#: Local-noon night labels bounding the 2023 campaign.  The manifest's
#: ``night`` column splits at local noon, so the label runs one day behind
#: the UT date of a post-midnight exposure: night ``2023-05-19`` is the UT
#: 2023-05-20 discovery-week epoch the strategy calls "+1.6 d".
CAMPAIGN_FIRST_NIGHT = "2023-05-19"
CAMPAIGN_LAST_NIGHT = "2023-07-07"

#: The pre-explosion template epoch (night label; UT 2023-05-05).
TEMPLATE_PRE_NIGHT = "2023-05-04"

#: Filter codes, grouped by the role they play in this paper.  These are
#: MaxIm single-character wheel codes whose physical identity is still
#: unproven (that is Step 1's job); the grouping here is by role only —
#: which product a frame could feed — and never by an assumed bandpass.
BROADBAND_FILTERS = ("G", "R", "I")
NARROWBAND_FILTERS = ("H", "O", "1")
#: ``6`` is the mixed slot S2c measures per frame; ``X`` and ``L`` are the
#: luminance/unknown pair the strategy drops from photometry outright.
OTHER_FILTERS = ("6", "X", "L")

#: Fractions of the MEASURED per-mode clip at which a peak-pixel measurement
#: stops being trustworthy.
#:
#: The strategy (§3.3) states the screen in raw ADU — "reject above 2,800,
#: flag 2,400–2,800" — but it derived both numbers by eye from an ASSUMED
#: ~3,500 ADU clip, because the S2 tables that hold the measured ceiling had
#: been destroyed when it was written.  S2 has since been rebuilt and
#: measures the High Gain clip at 3,496 ADU.  Storing the fractions rather
#: than the ADU values means the screen follows the measurement: applied to
#: the measured clip these give 2,797 and 2,398, reproducing the strategy's
#: hand numbers to better than 0.2% — which is itself the check that the
#: strategy's assumed ceiling was right.
REJECT_CLIP_FRACTION = 0.80          # 2,800 / 3,500
SUSPECT_CLIP_FRACTION = 2400.0 / 3500.0

#: Half-width in pixels of the small stamp whose maximum is taken as the
#: supernova's peak when the position is known exactly.  At the campaign's
#: median FWHM of 3.5–4.0 px a 21 × 21 stamp contains the whole PSF with
#: room for the ~1 px WCS residual, and is far too small to reach a
#: neighbouring star.
CORE_HALF_PX = 10

#: Half-width of the large box whose maximum is the UPPER BOUND used when
#: the position is not known exactly.  It must exceed the worst commanded-
#: versus-true pointing residual in the campaign, which the freeze measures
#: from the plate-solved frames themselves (~240 px about the campaign
#: median).  250 px is that residual rounded up.
BOUND_HALF_PX = 250

#: Half-width of the annulus box used to estimate the local sky pedestal.
SKY_HALF_PX = 60

#: A frame is only claimed as a per-frame supernova measurement when the
#: brightest pixel found in the core stamp sits within this many pixels of
#: the predicted position.  Beyond it the stamp maximum is some other
#: source, and the row is downgraded to an upper bound.
MAX_CENTROID_OFFSET_PX = 6.0

#: Flash-phase window in days after explosion.  PROVENANCE: the flash
#: (highly ionised, narrow, electron-scattering-winged) features in
#: SN 2023ixf are reported gone by ~day 7–8 (Jacobson-Galán+23 ApJL 954 L42;
#: Bostroem+23; Hiramatsu+23).  Used to COUNT nights, never to judge one.
FLASH_PHASE_END_D = 8.0

#: The strategy's own pre-registered grism promotion criterion (§4 Step 0c):
#: Hα emission with electron-scattering wings at SNR > 10 on at least this
#: many nights, AND a passed offset-trace contamination test.  Both clauses
#: are required; this constant is only the first one's threshold.
GRISM_PROMOTION_MIN_NIGHTS = 3

#: The two venue postures the strategy allows, and the rule between them
#: (§1, §2): AJ/PASP is the decided base case; ApJ is taken ONLY if Gate 0
#: promotes the grism or recovers the narrowband bandpass.
VENUE_BASE = "AJ/PASP"
VENUE_UPSIDE = "ApJ"


# ===========================================================================
# 2.  COORDINATES — turning a header into a pixel
# ===========================================================================

def parse_sexagesimal(text, is_hours: bool) -> Optional[float]:
    """Parse ``'14 03 38.56'`` / ``'+54:18:41.9'`` into degrees.

    ``is_hours`` multiplies by 15 (an RA written in hours).  Returns None
    for anything unparseable, because a missing pointing card is a normal
    state in this archive and must not raise in the middle of a 1,100-frame
    campaign pass.

    The sign is taken from the STRING, not from the degrees field: a
    declination of ``-00 30 00`` has a zero degrees field whose float
    carries no sign, and reading the sign off the number would silently
    flip it to +00°30'.
    """
    if text is None:
        return None
    s = str(text).strip()
    if not s:
        return None
    # Split on whitespace or colons — MaxIm writes both, sometimes mixed.
    parts = re.split(r"[\s:]+", s)
    try:
        # Take up to three fields and pad the missing ones with zero, so
        # '14 03' and '14' parse as 14h03m00s and 14h00m00s.
        vals = [float(p) for p in parts[:3]]
    except ValueError:
        return None
    while len(vals) < 3:
        vals.append(0.0)
    sign = -1.0 if s.lstrip().startswith("-") else 1.0
    deg = abs(vals[0]) + vals[1] / 60.0 + vals[2] / 3600.0
    return sign * deg * (15.0 if is_hours else 1.0)


def gnomonic_pixel(ra_deg: float, dec_deg: float,
                   ra0_deg: float, dec0_deg: float,
                   scale_deg_px: float,
                   crpix_x: float, crpix_y: float) -> tuple[float, float]:
    """Where does (ra, dec) fall on a north-up tangent-plane frame centred
    on (ra0, dec0)?  Returns 0-indexed (x, y) pixels.

    This is the FALLBACK projection, used only when a frame carries no plate
    solution of its own.  It assumes north up, east left and a square pixel
    — all three true of this camera to within the ~1° rotation the solved
    frames show, which over the 250 px search box costs at most ~4 px.  It
    is never used to claim a position: it seeds a search whose success is
    then checked against the found peak (:func:`position_quality`).
    """
    # Standard gnomonic (TAN) projection.  Everything in radians.
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    ra0 = math.radians(ra0_deg)
    dec0 = math.radians(dec0_deg)
    cos_c = (math.sin(dec0) * math.sin(dec) +
             math.cos(dec0) * math.cos(dec) * math.cos(ra - ra0))
    if cos_c <= 0:
        # More than 90° from the tangent point: the projection diverges and
        # the answer would be meaningless rather than merely imprecise.
        raise ValueError("position is on the far side of the tangent plane")
    xi = math.cos(dec) * math.sin(ra - ra0) / cos_c
    eta = ((math.cos(dec0) * math.sin(dec) -
            math.sin(dec0) * math.cos(dec) * math.cos(ra - ra0)) / cos_c)
    # Radians on the tangent plane -> degrees -> pixels.  East is left, so
    # the xi axis runs the opposite way to the pixel x axis.
    x = crpix_x - math.degrees(xi) / scale_deg_px
    y = crpix_y + math.degrees(eta) / scale_deg_px
    return (x, y)


def position_quality(offset_px: Optional[float],
                     has_wcs: bool,
                     max_offset_px: float = MAX_CENTROID_OFFSET_PX) -> str:
    """Classify how well this frame's supernova position is known.

    ``wcs``       the frame carried its own plate solution AND the brightest
                  pixel in the core stamp landed within ``max_offset_px`` of
                  the predicted position: the measurement is of the
                  supernova, and the frame's WCS is corroborated by it.
    ``wcs_off``   the frame carried a plate solution but no well-centred
                  peak sits where the supernova should be.  Three things
                  produce this and the census does not try to tell them
                  apart: the supernova is too faint to dominate its own
                  stamp, the point-spread function is so broad that the
                  blob's maximum wanders (the pre-explosion template epoch,
                  whose 14–16 px FWHM the strategy already flagged), or the
                  solve is wrong.  In all three the CORE-STAMP maximum is
                  still a valid upper bound on the supernova's peak — the
                  supernova is inside that stamp — so that is what such a
                  row reports.
    ``bound``     no plate solution at all: the supernova's position is
                  known only to the pointing residual, so only the wide
                  search box's maximum is defensible, and it is a much
                  looser bound.

    A frame is never promoted from ``bound`` on the strength of finding
    *something* bright in the box, because in the narrowband bands a field
    star outshines the supernova often enough to make that inference wrong
    about half the time — a rate the census measures rather than assumes.
    """
    if not has_wcs:
        return "bound"
    if offset_px is None:
        return "wcs_off"
    return "wcs" if offset_px <= max_offset_px else "wcs_off"


# ===========================================================================
# 3.  THE DISPERSION GATE — measurement, never label
# ===========================================================================

def dispersion_class(verdict) -> str:
    """The S2c per-frame verdict, normalised.  Delegates to
    ``macro_core.astrom`` so the two stages cannot drift apart."""
    return astrom.dispersion_class(verdict)


def is_direct_image(verdict) -> bool:
    """THE imaging gate: is this frame an IMAGE at all?

    True for everything except a frame S2c MEASURED as dispersed.  The
    asymmetry is deliberate and is the whole correction: excluding a frame
    requires a positive measurement that it is a spectrum, while including
    one requires only that no such measurement exists.  Under the retired
    label rule this campaign's 83 slot-``6`` frames were all called spectra,
    which deleted 3 measured direct images unseen; and 19 frames nobody
    could certify were being treated as though somebody had.
    """
    return not astrom.is_measured_spectrum(verdict)


def band_role(filt) -> str:
    """Which product could this frame feed: ``broadband`` / ``narrowband`` /
    ``other``.  A grouping by ROLE, not by an assumed bandpass — the filter
    codes' physical identity is Step 1's job and is not assumed here."""
    f = (filt or "").strip()
    if f in BROADBAND_FILTERS:
        return "broadband"
    if f in NARROWBAND_FILTERS:
        return "narrowband"
    return "other"


def epoch_role(night: str) -> str:
    """Which epoch of the project a night belongs to.

    ``campaign``       the 2023 monitoring run Gate 0 censuses.
    ``template_pre``   the single pre-explosion epoch (UT 2023-05-05).
    ``template_post``  everything after the campaign: the post-fade
                       template epochs and the late-time anchors.
    """
    if night == TEMPLATE_PRE_NIGHT:
        return "template_pre"
    if CAMPAIGN_FIRST_NIGHT <= night <= CAMPAIGN_LAST_NIGHT:
        return "campaign"
    if night < CAMPAIGN_FIRST_NIGHT:
        return "template_pre"
    return "template_post"


def phase_days(jd: Optional[float], t0_mjd: float = T0_MJD) -> Optional[float]:
    """Days since the adopted explosion epoch.  DISPLAY ONLY (see T0_MJD)."""
    if jd is None:
        return None
    return (jd - MJD_OFFSET) - t0_mjd


# ===========================================================================
# 4.  THE SATURATION SCREEN
# ===========================================================================

@dataclass(frozen=True)
class Screen:
    """The three ADU levels that judge a peak, all derived from the MEASURED
    per-mode ceiling rather than typed.

    ``clip_adu``     where S2 measured this readout mode's ADU scale to end.
    ``veto_adu``     S2's own adopted saturation veto for the mode.
    ``reject_adu``   the strategy's 80%-of-clip screen, applied to the
                     measured clip.
    ``suspect_adu``  the strategy's lower flag level, likewise.
    """

    mode: str
    clip_adu: int
    veto_adu: Optional[int]
    reject_adu: int
    suspect_adu: int


def screen_for_mode(mode: str, clip_adu, veto_adu,
                    reject_fraction: float = REJECT_CLIP_FRACTION,
                    suspect_fraction: float = SUSPECT_CLIP_FRACTION) -> Screen:
    """Build the screen for one readout mode from S2's measured numbers.

    Raises when the clip is missing: a census that silently fell back to a
    typed 2,800 would publish a screen with nothing behind it, which is the
    precise state the strategy flagged as blocking (task SN-G0b was BLOCKED
    on exactly this missing measurement).
    """
    if clip_adu is None:
        raise ValueError(f"no measured clip for readout mode {mode!r}; "
                         f"S2 must supply s2_ceiling_modes.clip_adu before "
                         f"a saturation screen can be applied")
    clip = int(clip_adu)
    return Screen(mode=mode, clip_adu=clip,
                  veto_adu=None if veto_adu is None else int(veto_adu),
                  reject_adu=int(round(reject_fraction * clip)),
                  suspect_adu=int(round(suspect_fraction * clip)))


def saturation_class(peak_adu: Optional[float], quality: str,
                     screen: Screen) -> str:
    """Judge one frame's supernova peak against the screen.

    The five outcomes, and why five rather than three:

    ``clean``          peak measured on the supernova, below the flag level.
    ``suspect``        measured, between the flag level and the reject
                       screen — usable only once the linearity curve exists.
    ``rejected``       measured, at or above the reject screen: the pixel is
                       in the non-linear or clipped regime and no photometry
                       may be taken from it.
    ``bounded_clean``  the position was NOT known, but the whole search box
                       stays under the flag level — so wherever the
                       supernova fell in it, it was unsaturated.  A genuine
                       one-sided measurement.
    ``undetermined``   the position was not known and the box does reach the
                       screen.  Something in that box is bright; the census
                       declines to say it was the supernova.

    The last class is the honest one.  Folding it into ``rejected`` would
    manufacture saturated epochs out of field stars, and folding it into
    ``clean`` would licence photometry on frames nobody has checked.
    """
    if peak_adu is None:
        return "undetermined"
    if quality == "wcs":
        # Position known and corroborated: this IS the supernova's peak.
        if peak_adu >= screen.reject_adu:
            return "rejected"
        if peak_adu >= screen.suspect_adu:
            return "suspect"
        return "clean"
    # Position not established.  The number in hand is a box maximum, which
    # can only ever be an UPPER bound on the supernova's peak, so only the
    # "under the screen" direction can be concluded.
    if peak_adu < screen.suspect_adu:
        return "bounded_clean"
    return "undetermined"


#: The classes from which broadband photometry may be taken.  ``suspect``
#: is deliberately NOT here: the strategy makes it conditional on a
#: linearity curve that Step 2 has not yet produced.
USABLE_CLASSES = ("clean", "bounded_clean")

#: The archive tree that holds the campaign's science exposures.  The
#: strategy's own path-canonicalization rule (§4 Step 0a) names it: "the
#: ``rawimage`` copy wins".
#:
#: The clause exists because the census found what the rule was for.  Two
#: canonical frames on this sky live under ``mjc/misc/neg10_test/`` and are
#: named ``*.raw.fts.fz`` — detector engineering products from a negative-
#: temperature test that happen to carry the campaign's target name and
#: fall inside the campaign's dates.  One of them reads 65,535 ADU at the
#: supernova's position, a value the 12-bit High Gain channel cannot
#: physically produce, which is how they were caught.  Without this clause
#: the other one would have been counted as a clean broadband epoch.
SCIENCE_TREE = "rawimage"


def is_usable_photometry(row: dict) -> bool:
    """The composite Gate 0 test for question 1 — is this frame genuinely
    usable broadband photometry?

    Five clauses, each of which has to be true, and each of which is a
    measurement or a stated manifest rule rather than an assumption:

    1. it belongs to the campaign (not a template or late-time epoch);
    2. it comes from the science tree, not an engineering directory;
    3. S2c did not measure it as a spectrum;
    4. its filter is one of the three broadband codes;
    5. the pixel census puts the supernova's peak below the flag level,
       either directly or by a bound.
    """
    return (row.get("epoch_role") == "campaign"
            and row.get("tree") == SCIENCE_TREE
            and is_direct_image(row.get("dispersion_verdict"))
            and band_role(row.get("filter")) == "broadband"
            and row.get("saturation_class") in USABLE_CLASSES)


# ===========================================================================
# 5.  ARITHMETIC — rates, verdicts, intervals
# ===========================================================================

def wilson_ci(k: int, n: int) -> tuple[float, float]:
    """95% Wilson interval; delegates to ``macro_core.astrom``."""
    return astrom.wilson_ci(k, n)


def astrometry_verdict(k: int, n: int,
                       n_population: Optional[int] = None) -> str:
    """GO / CAUTION / NO-GO on an astrometric success rate.

    Delegates to ``macro_core.astrom.verdict_for`` ON PURPOSE.  The question
    Gate 0 is asked — "does the astrometry verdict move?" — is only
    answerable if the Gate 0 number and the S1 number are judged by the SAME
    thresholds.  A second copy of 0.80 / 0.50 here would let the two answers
    diverge without either page noticing.
    """
    return astrom.verdict_for(k, n, n_population)


def rate_pct(k: int, n: int) -> Optional[float]:
    """k/n as a percentage, or None for an empty denominator."""
    if not n:
        return None
    return 100.0 * k / n


def isolation_false_id_rate(isolation_px: Sequence[Optional[float]],
                            radius_px: float) -> tuple[int, int]:
    """How often would "brightest pixel in the box" NOT be the supernova?

    ``isolation_px`` is, per plate-solved frame, the distance from the
    supernova to the nearest pixel BRIGHTER than it.  A search box of
    half-width ``radius_px`` picks the wrong source exactly when that
    distance is STRICTLY smaller than the box.  Returns (n_wrong, n_tested).

    The comparison is strict on purpose.  When nothing in the search box
    outbrightens the supernova the census cannot measure how far away the
    next brighter pixel is, so it records the box half-width itself as a
    FLOOR ("at least this isolated").  A non-strict comparison would read
    every one of those floors as a confusion and turn the cleanest possible
    result — the supernova is the brightest thing for 250 px in every
    direction — into a 100% false-identification rate.

    This is the calibration that lets the census say what a box maximum is
    worth in each band, measured on the frames where the answer is known,
    instead of assuming the supernova is always the brightest thing near it.
    """
    tested = [d for d in isolation_px if d is not None]
    wrong = [d for d in tested if d < radius_px]
    return (len(wrong), len(tested))


# ===========================================================================
# 6.  GATE 0c — the grism triage
# ===========================================================================

@dataclass(frozen=True)
class GrismSeries:
    """What the measurement says about the slitless series.

    Every field is a count or a flag derived from the manifest; none is an
    opinion.  ``n_extracted`` and ``n_contamination_passed`` exist so the
    promotion rule has somewhere to read the two things Gate 0c has NOT
    done, rather than leaving their absence implicit.
    """

    n_labelled: int                  # frames carrying the slot-'6' label
    n_dispersed: int                 # ...that S2c measures as spectra
    n_direct: int                    # ...that S2c measures as IMAGES
    n_indeterminate: int             # ...S2c could not certify
    n_nights: int                    # nights holding >=1 measured spectrum
    n_flash_nights: int              # ...of those, at or before day 8
    n_nights_with_paired_direct: int  # nights that also hold direct images
    n_extracted: int                 # 1-D spectra extracted so far
    n_contamination_passed: int      # ...that passed the offset-trace test
    wavelength_source: str           # named calibration source, or ''
    n_flats: int                     # filter-6 flat frames in the archive


def grism_promotion(series: GrismSeries,
                    min_nights: int = GRISM_PROMOTION_MIN_NIGHTS) -> dict:
    """Apply the strategy's OWN pre-registered promotion criterion.

    §4 Step 0c: promote to a headline figure and a dedicated section only on
    "Hα emission with electron-scattering wings at SNR > 10 across >= 3
    nights AND a passed offset-trace contamination test".  Otherwise: one
    appendix paragraph and no further effort.

    Returns the verdict together with every clause's state, so the page can
    show WHICH clause failed rather than only that one did.  The rule is
    applied exactly as pre-registered: a series that is real, well sampled
    and never extracted does not promote, because the criterion is about
    measured spectra and not about the frames that could yield them.
    """
    clauses = {
        # Clause 1 — enough nights of spectra to be a series at all.  This
        # one is answerable from the manifest today.
        "nights": {
            "requirement": f"measured spectra on >= {min_nights} nights",
            "value": series.n_nights,
            "passed": series.n_nights >= min_nights,
        },
        # Clause 2 — the flash phase actually covered.  Not part of the
        # pre-registered bar, but the reason the product was proposed, so
        # it is reported beside it and never allowed to substitute for it.
        "flash_nights": {
            "requirement": "flash-phase nights covered (context, not a bar)",
            "value": series.n_flash_nights,
            "passed": series.n_flash_nights >= 1,
        },
        # Clause 3 — the extraction itself.
        "extracted": {
            "requirement": f"Hα at SNR > 10 on >= {min_nights} nights",
            "value": series.n_extracted,
            "passed": series.n_extracted >= min_nights,
        },
        # Clause 4 — the contamination test, mandatory and independent.
        "contamination": {
            "requirement": "offset-trace contamination test passed",
            "value": series.n_contamination_passed,
            "passed": series.n_contamination_passed >= min_nights,
        },
        # Clause 5 — a named, dated wavelength solution or a declared
        # self-calibration.
        "wavelength": {
            "requirement": "named/dated wavelength calibration source",
            "value": series.wavelength_source or "none identified",
            "passed": bool(series.wavelength_source),
        },
    }
    # The BAR is clauses 3, 4 and 5 — the ones about spectra.  Clause 1 is
    # necessary but not sufficient, and clause 2 is context.
    bar = ("extracted", "contamination", "wavelength")
    promoted = (clauses["nights"]["passed"]
                and all(clauses[c]["passed"] for c in bar))
    return {
        "clauses": clauses,
        "bar": bar,
        "promoted": promoted,
        "blocking": [c for c in ("nights", *bar) if not clauses[c]["passed"]],
    }


def venue_posture(grism_promoted: bool, bandpass_recovered: bool) -> dict:
    """The strategy's own venue rule, applied to Gate 0's outcome.

    §1 and §2 decide the posture in advance and in one sentence: AJ/PASP is
    the base case; ApJ is taken ONLY if Gate 0 promotes the grism or
    recovers the narrowband bandpass.  Both conditions are disjunctive, both
    are about things Gate 0 either did or did not achieve, and neither is a
    judgement call — which is why the rule was written down before the
    evidence arrived.
    """
    promote = bool(grism_promoted) or bool(bandpass_recovered)
    return {
        "posture": VENUE_UPSIDE if promote else VENUE_BASE,
        "moved": promote,
        "grism_promoted": bool(grism_promoted),
        "bandpass_recovered": bool(bandpass_recovered),
    }


# ===========================================================================
# 7.  SQL the plumbing runs — kept here so the row scope is reviewable
# ===========================================================================

#: The campaign's canonical-target set.  Both names are the SAME SKY: S0's
#: alias merge maps ``NGC5457``, ``M101``, ``ngc 5457`` and ``pinwheel
#: galaxy`` onto canonical target ``M101``, while frames whose OBJECT card
#: already named the transient map onto ``2023ixf``.  A census that took
#: only the second would silently drop the discovery-week frames, which is
#: exactly the error the strategy's own headline count (1,052) contains.
CAMPAIGN_TARGETS = ("2023ixf", "M101")


def targets_sql() -> str:
    """SQL IN-list of the campaign's canonical target names."""
    return ", ".join("'" + t.replace("'", "''") + "'"
                     for t in CAMPAIGN_TARGETS)


#: One row per UNIQUE frame on this sky, at any epoch, with the S2c verdict
#: attached.  ``is_canonical = 1`` is S0's global dedup flag — one row per
#: (basename, jd) group across every tree, including within ``rawimage``.
FREEZE_SQL = f"""
SELECT  f.obs_rowid, f.path, f.tree, f.basename, f.night, f.jd,
        f.filter, f.exptime, f.readoutm, f.xbinning, f.naxis1, f.naxis2,
        f.canonical_target, f.object, f.observer, f.imagetyp,
        f.pltsolvd, f.crval1, f.crval2, f.objctra, f.objctdec,
        f.fwhm, f.airmass, f.zmag, f.dup_group, f.era_id,
        d.verdict AS dispersion_verdict,
        d.status  AS dispersion_status
FROM    frames f
LEFT JOIN frame_dispersion d ON d.obs_rowid = f.obs_rowid
WHERE   f.canonical_target IN ({targets_sql()})
  AND   f.is_canonical = 1
  AND   (f.imagetyp IS NULL OR f.imagetyp LIKE '%Light%')
ORDER BY f.jd
"""

#: Every catalog row on this sky, canonical or not — the denominator of the
#: dedup accounting.  Gate 0a's whole point is that the difference between
#: this count and the one above is 3x, and that the duplicates live INSIDE
#: ``rawimage`` as well as across trees.
DEDUP_SQL = f"""
SELECT  f.tree,
        count(*)                          AS n_rows,
        sum(CASE WHEN f.is_canonical = 1 THEN 1 ELSE 0 END) AS n_canonical
FROM    frames f
WHERE   f.canonical_target IN ({targets_sql()})
  AND   (f.imagetyp IS NULL OR f.imagetyp LIKE '%Light%')
GROUP BY f.tree
ORDER BY n_rows DESC
"""


def freeze_rows(raw_rows: Iterable[dict]) -> list[dict]:
    """Decorate the frozen manifest rows with every derived column.

    Pure: takes dicts, returns dicts.  Everything the rest of Gate 0 keys on
    — epoch role, band role, dispersion class, whether the frame is an image
    at all, phase — is computed HERE and once, so no downstream query can
    apply a different rule to the same frame.
    """
    out = []
    for r in raw_rows:
        night = r.get("night") or ""
        row = dict(r)
        row["epoch_role"] = epoch_role(night)
        row["band_role"] = band_role(r.get("filter"))
        row["dispersion_class"] = dispersion_class(r.get("dispersion_verdict"))
        row["is_image"] = 1 if is_direct_image(r.get("dispersion_verdict")) else 0
        row["phase_d"] = phase_days(r.get("jd"))
        out.append(row)
    return out
