"""Per-frame, empirical filter identity: is this frame DISPERSED or DIRECT?

WHY THIS MODULE EXISTS
----------------------
The RLMT filter wheel carries compact grism spectrometers alongside ordinary
photometric filters, and the FILTER header card that names them changed
vocabulary three times in three years:

    2023-02 -> 2024-03   slot number ``6``
    2024-04 -> 2025-03   ``HaGrism`` / ``OGGrism``
    2025-01 -> 2026-06   ``hrg`` / ``lrg``

Worse, the slot-number epoch is not trustworthy at all: ``6`` reads as a
dispersed spectrum on some targets and as an ordinary direct image on
others.  Two committee panels read the same card and reached opposite
conclusions (one called it "grism", one called it "luminance") because the
card alone genuinely cannot settle it.  Every downstream project — T CrB,
the BeStar campaign, SN 2023ixf, the Dwarf/AGN survey — needs to know which
of its frames are spectra, and a label nobody can verify is not an answer.

So we stop reading the label and measure the PIXELS.

THE PHYSICS BEING DETECTED
--------------------------
A grism is a slitless disperser: it spreads every source in the field into a
line along ONE fixed axis set by the grating rulings.  Three consequences
follow, and this module measures all three:

1. ELONGATION.  A dispersed source is enormously eccentric.  Direct frames
   put every star at the seeing disc's axis ratio, a/b ~= 1.0-1.3 even when
   the star saturates.  Measured traces run a/b = 30-150.  There is no
   overlap, and the gap is more than an order of magnitude wide.

2. A SHARED AXIS.  Every trace on a grism frame points the SAME way,
   because they are all ruled by the same grating.  This is the signature
   that separates a grism from the two things that also elongate sources:
   a cosmic-ray hit (one streak, random angle) and a wind-shaken or
   badly-tracked exposure (all sources elongated, but at the drift angle
   and — crucially — only mildly, a/b of a few).  Position angle is AXIAL,
   not directional: a trace at +90 deg and one at -90 deg lie on the same
   axis, so all circular statistics here run on DOUBLED angles, where that
   wrap closes correctly.

3. A MIXED POPULATION.  A grism field is not uniformly streaky.  The bright
   target throws a long first-order trace, the undispersed zeroth order sits
   beside it as a round blob, and faint field stars are detected only at
   their trace peak and also look round.  A frame therefore reads as
   dispersed on the strength of its BRIGHT sources; demanding that most
   sources be streaky would reject real spectra.

That third point is why the classifier does not use a plain median over all
detections.  The measured median a/b over the ten brightest sources of a
real ``hrg`` frame can be as low as 2.5 purely because seven of those ten
are faint round field stars sharing the frame with two 900-pixel traces.
The median is still recorded (it is the statistic the survey's first pass
quoted, and the report plots it), but the VERDICT is decided by counting
trace-like sources and testing whether they share an axis.

THE SECOND AXIS: DISPERSION STRENGTH
------------------------------------
Two different grisms flew.  The high-dispersion H-alpha unit (``hrg`` /
``HaGrism``) spreads a given star roughly twice as far as the broad-spectrum
low-dispersion unit (``lrg`` / ``OGGrism``).  ``classify_strength`` reads the
trace aspect ratio and reports ``high`` / ``low`` / ``ambiguous``.  Aspect
ratio is used rather than raw trace length in pixels because the archive
spans three cameras and two binnings, and a dimensionless ratio survives
that where a pixel count does not.

The honest caveat, which the report states in full: trace aspect is
confounded by source brightness and exposure time, because a brighter star's
trace is detected further into its faint wings before it drops below the
extraction threshold.  The two grism populations therefore OVERLAP, and
``classify_strength`` returns ``ambiguous`` across that overlap instead of
guessing.  The calibrated numbers are in the report.

WHAT IS PURE HERE
-----------------
Everything that decides anything.  ``axial_stats``, ``select_bright``,
``trace_flags``, ``summarize_sources``, ``classify_frame`` and
``classify_strength`` see nothing but plain numpy arrays and scalars, so the
unit tests drive them with hand-built source lists whose truth is known —
including the cases that MUST NOT read as a grism (a round star field, a
lone cosmic-ray streak, a trailed exposure, an empty frame).  Only
``extract_sources`` and ``measure_file`` touch pixels or disk.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Optional, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Version stamp recorded into every frame_dispersion row.  Bump whenever a
# change here would alter a stored verdict, so a later reader can tell which
# rules produced the numbers in front of them.
# ---------------------------------------------------------------------------
DISPERSION_CODE_VERSION = "S2c v1.0 (2026-08-18)"

# ---------------------------------------------------------------------------
# Tunable constants — single source of truth; the report interpolates these
# rather than repeating them, so page and code can never disagree.
# ---------------------------------------------------------------------------

#: Detection threshold in units of the background RMS.  5 sigma is the
#: standard faint-source floor; the escalation ladder below raises it when a
#: crowded or nebulous frame overflows the extractor's pixel stack.
DETECT_SIGMA = 5.0

#: Threshold ladder tried in order.  A deep field (or a reduced frame with a
#: flattened, near-zero background) can put more than the extractor's
#: 300k-pixel object budget above 5 sigma; rather than lose the frame we
#: re-extract brighter.  The threshold actually used is recorded per frame.
DETECT_SIGMA_LADDER = (5.0, 15.0, 50.0, 150.0)

#: Background mesh (px).  64 px is small enough to track the vignetting and
#: sky gradients of a 4800-px frame, and large enough that a 9-px-wide,
#: 900-px-long trace never fills a mesh cell and gets absorbed into "sky".
BACK_SIZE = 64

#: Background median-filter size in mesh cells — smooths the mesh so one
#: bright box does not punch a hole in the sky model.
BACK_FILTER = 3

#: Minimum connected pixels for the extractor to call something a source.
EXTRACT_MINAREA = 9

#: Sources smaller than this are discarded before any statistic is formed.
#: Cosmic-ray hits and hot pixels live below it; a real trace is thousands
#: of pixels and a real star at 5 sigma is comfortably above it.
MIN_NPIX = 20

#: A source with a degenerate minor axis cannot give a meaningful ratio.
MIN_MINOR_AXIS_PX = 0.5

#: How many of the brightest sources the frame statistics are formed over.
#: Ten is enough to be robust against one weird detection and few enough
#: that a bright target's trace is not drowned by field stars.
BRIGHT_N = 10

#: A source counts as a TRACE when it is both very eccentric AND long in
#: absolute pixels.  Both gates matter: the ratio alone would admit a
#: 3-px-long sliver with a sub-pixel minor axis, and the length alone would
#: admit a big fat galaxy.
TRACE_MIN_AB = 5.0
TRACE_MIN_A_PX = 15.0

#: How tightly the traces on one frame must share an axis, in degrees of
#: axial circular scatter.  Real grism frames measure a few degrees; a
#: random-angle collection scatters near the 52-deg ceiling of the statistic.
TRACE_MAX_PA_SCATTER_DEG = 20.0

#: Two traces sharing an axis is already a strong statement.  One trace
#: cannot be checked for a shared axis at all, so a solo trace must clear a
#: much higher eccentricity bar before it is believed — this is the case of
#: a bright spectrophotometric standard alone in the field, which is a large
#: and legitimate part of the archive.
MIN_TRACES_FOR_AXIS = 2
SOLO_TRACE_MIN_AB = 20.0

#: ...and the solo rule additionally requires a SPARSE field.  This gate was
#: added after a 60-s luminance frame of M57 — 1,278 detected sources, every
#: bright one round — was called "dispersed" on the strength of a single
#: 1,095-px streak across it.  That streak was a satellite, and the giveaway
#: is the 1,278 round stars it flew past: a grating disperses EVERYTHING, so
#: a rich field with exactly one smear is the one thing a grism cannot
#: produce.  Measured sparse grism frames (focus sweeps, bright standards)
#: carry 1-7 usable sources; the crowded false positive carried 1,278.  The
#: gate sits far above the former and far below the latter.
SOLO_MAX_SOURCES = 60

#: A frame with fewer usable sources than this is not evidence of anything;
#: it returns ``indeterminate`` rather than a guess.
MIN_SOURCES = 1

#: A frame with NO trace-like source reads as direct only if its bright
#: sources really are round.  Above this median ratio something is elongating
#: the field (trailing, wind, focus astigmatism) and the frame is called
#: indeterminate instead of direct — we decline to certify it as clean
#: photometry when we can see that something is wrong with it.
DIRECT_MAX_MEDIAN_AB = 3.0

#: Dispersion-strength split.  The measured quantity is TRACE LENGTH AS A
#: FRACTION OF FRAME WIDTH — the trace's semi-major sigma divided by the
#: decoded frame width — not the aspect ratio.
#:
#: Aspect ratio was the obvious first choice and it does not work.  On the
#: 2025-01-23 focus sweep, where the SAME star was shot through both grisms
#: on the same night with the same camera, the two units measured a/b = 83
#: (high) against a/b = 61 (low): overlapping distributions, useless as a
#: split.  The same frames measured trace length fractions of 0.206 against
#: 0.060 — a clean factor of 3.4.  The reason is that a/b divides by the
#: minor axis, which is the seeing width, so every change in focus or
#: atmosphere feeds straight into the statistic; length alone does not care.
#: Dividing by frame width (rather than using raw pixels) is what lets one
#: pair of constants span three cameras and two binnings.
STRENGTH_LOW_MAX_FRAC = 0.11
STRENGTH_HIGH_MIN_FRAC = 0.15

#: Verdict vocabulary — the only three strings ``classify_frame`` may emit.
VERDICT_DISPERSED = "dispersed"
VERDICT_DIRECT = "direct"
VERDICT_INDETERMINATE = "indeterminate"

#: Known-truth label sets used to calibrate and to score the classifier.
#: These are the labels whose meaning nobody disputes.  Slot '6' and 'W' are
#: deliberately absent: they are the QUESTION, not the calibration.
KNOWN_DISPERSED_FILTERS = ("hrg", "lrg", "HaGrism", "OGGrism", "HaG")
KNOWN_DIRECT_FILTERS = ("g", "r", "i", "V", "R", "I", "B", "L")

#: Which known-dispersed labels are the high-dispersion H-alpha unit and
#: which are the low-dispersion broad-spectrum one.
HIGH_DISPERSION_FILTERS = ("hrg", "HaGrism", "HaG")
LOW_DISPERSION_FILTERS = ("lrg", "OGGrism")


# ---------------------------------------------------------------------------
# Pure geometry: axial (mod-180) circular statistics
# ---------------------------------------------------------------------------
def axial_stats(theta_deg: Sequence[float]) -> tuple[Optional[float],
                                                     Optional[float]]:
    """Mean axis and scatter of a set of ORIENTATIONS, in degrees.

    An orientation is an axis, not an arrow: 0 deg and 180 deg are the same
    line, and so are +89 and -89 (they differ by 2 deg, not 178).  Averaging
    such angles arithmetically is simply wrong — the mean of {+89, -89} is 0,
    perpendicular to both inputs.

    The fix is the standard doubling trick: map each angle t to 2t, where the
    180-deg ambiguity becomes a full 360-deg turn and ordinary circular
    statistics apply; average the unit vectors there; halve the result back.

    Returns ``(pa_deg, scatter_deg)`` with pa in [0, 180) and scatter in
    [0, 90], or ``(None, None)`` for an empty input.  Scatter is the circular
    standard deviation sqrt(-2 ln R) of the doubled angles, halved back into
    orientation space; R is the resultant length, so perfectly aligned axes
    give 0 and uniformly random ones approach the 90-deg ceiling.
    """
    arr = np.asarray(theta_deg, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return (None, None)
    # Double the angles, average as unit vectors on the circle.
    doubled = np.radians(2.0 * arr)
    c = float(np.mean(np.cos(doubled)))
    s = float(np.mean(np.sin(doubled)))
    resultant = math.hypot(c, s)
    # Halve the mean direction back, and fold into [0, 180).
    pa = (math.degrees(math.atan2(s, c)) / 2.0) % 180.0
    if resultant <= 1e-12:
        # Perfectly opposed axes: the mean direction is undefined, and the
        # scatter is maximal.  Report the ceiling rather than an infinity.
        return (pa, 90.0)
    scatter = math.degrees(math.sqrt(max(-2.0 * math.log(resultant), 0.0))) / 2.0
    return (pa, min(scatter, 90.0))


# ---------------------------------------------------------------------------
# Pure selection: which detections are usable, and which are traces
# ---------------------------------------------------------------------------
def usable_mask(npix: Sequence[float], b: Sequence[float],
                min_npix: int = MIN_NPIX,
                min_minor: float = MIN_MINOR_AXIS_PX) -> np.ndarray:
    """Boolean mask of detections worth measuring at all.

    Two cuts, each rejecting a specific known impostor: ``npix`` rejects
    cosmic-ray hits and hot pixels (a few pixels each), and ``b`` rejects
    degenerate slivers whose minor axis would make a/b meaningless or
    infinite.
    """
    npix_a = np.asarray(npix, dtype=float)
    b_a = np.asarray(b, dtype=float)
    return (npix_a >= min_npix) & (b_a >= min_minor) & np.isfinite(b_a)


def select_bright(flux: Sequence[float], keep: int = BRIGHT_N) -> np.ndarray:
    """Indices of the ``keep`` brightest entries, brightest first.

    Brightness ordering is the whole game: a grism disperses every source,
    but only the bright ones carry enough flux per resolution element to be
    detected along their whole trace.  Faint stars in the SAME frame are
    detected as round stubs.  Ranking by flux puts the evidence first.
    """
    f = np.asarray(flux, dtype=float)
    if f.size == 0:
        return np.empty(0, dtype=int)
    order = np.argsort(f)[::-1]
    return order[:keep]


def trace_flags(a: Sequence[float], b: Sequence[float],
                min_ab: float = TRACE_MIN_AB,
                min_a: float = TRACE_MIN_A_PX) -> np.ndarray:
    """Boolean mask marking detections shaped like a dispersed trace.

    BOTH gates are required, and each one blocks a different false positive:
    the ratio gate alone would accept a 3-px sliver with a 0.5-px minor axis
    (a hot-pixel pair), and the absolute-length gate alone would accept a
    resolved galaxy or a defocused blob.  A dispersed spectrum is long AND
    thin, so it must clear both.
    """
    a_a = np.asarray(a, dtype=float)
    b_a = np.asarray(b, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(b_a > 0, a_a / b_a, 0.0)
    return (ratio >= min_ab) & (a_a >= min_a)


# ---------------------------------------------------------------------------
# The per-frame measurement record
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FrameShape:
    """Everything measured from one frame's pixels, before any judgement.

    Kept separate from the verdict on purpose: the measurement is expensive
    (it decompresses and extracts a 15-megapixel frame) while the verdict is
    free.  Storing the numbers means the thresholds can be recalibrated later
    over the whole archive without re-reading a single byte of pixel data.
    """
    n_detected: int          # raw extractor detections, before any cut
    n_sources: int           # detections surviving the usability cuts
    n_bright: int            # how many entered the bright-set statistics
    median_ab: Optional[float]      # median a/b over the bright set
    max_ab: Optional[float]         # largest a/b in the bright set
    pa_median: Optional[float]      # bright-set mean axis, deg [0,180)
    pa_scatter: Optional[float]     # bright-set axial scatter, deg
    n_trace: int                    # bright sources shaped like traces
    trace_frac: Optional[float]     # n_trace / n_bright
    trace_ab: Optional[float]       # median a/b over the traces
    trace_a_px: Optional[float]     # median semi-major sigma of traces, px
    trace_pa: Optional[float]       # trace mean axis, deg
    trace_pa_scatter: Optional[float]   # trace axial scatter, deg
    detect_sigma: float             # threshold that actually extracted
    height: int = 0                 # frame rows, as really decoded
    width: int = 0                  # frame columns, as really decoded

    def as_dict(self) -> dict:
        """Plain dict for the database writer."""
        return asdict(self)


#: An honest empty measurement, for a frame with nothing usable in it.
def empty_shape(detect_sigma: float = DETECT_SIGMA, n_detected: int = 0,
                height: int = 0, width: int = 0) -> FrameShape:
    """The measurement record of a frame that yielded no usable source."""
    return FrameShape(
        n_detected=n_detected, n_sources=0, n_bright=0,
        median_ab=None, max_ab=None, pa_median=None, pa_scatter=None,
        n_trace=0, trace_frac=None, trace_ab=None, trace_a_px=None,
        trace_pa=None, trace_pa_scatter=None,
        detect_sigma=detect_sigma, height=height, width=width)


def summarize_sources(a: Sequence[float], b: Sequence[float],
                      theta_deg: Sequence[float], flux: Sequence[float],
                      npix: Sequence[float],
                      detect_sigma: float = DETECT_SIGMA,
                      bright_n: int = BRIGHT_N,
                      height: int = 0, width: int = 0) -> FrameShape:
    """Reduce one frame's source list to the :class:`FrameShape` record.

    Pure: it sees five parallel arrays of source properties and nothing else,
    so the tests can feed it a hand-written grism field, a hand-written star
    field, a single cosmic ray, or nothing at all.

    ``theta_deg`` is the major-axis orientation in degrees; the extractor
    reports it in radians over [-pi/2, +pi/2], and the caller converts.
    """
    a_a = np.asarray(a, dtype=float)
    b_a = np.asarray(b, dtype=float)
    th_a = np.asarray(theta_deg, dtype=float)
    fl_a = np.asarray(flux, dtype=float)
    np_a = np.asarray(npix, dtype=float)
    n_detected = int(a_a.size)

    # Step 1 — drop the detections that cannot carry a shape measurement.
    keep = usable_mask(np_a, b_a)
    a_a, b_a, th_a, fl_a, np_a = (a_a[keep], b_a[keep], th_a[keep],
                                  fl_a[keep], np_a[keep])
    if a_a.size == 0:
        return empty_shape(detect_sigma, n_detected, height, width)

    # Step 2 — keep the brightest few; that is where a grism writes itself.
    idx = select_bright(fl_a, bright_n)
    a_b, b_b, th_b = a_a[idx], b_a[idx], th_a[idx]
    ratio = a_b / b_b

    # Step 3 — the bright-set summary (the statistic the first survey pass
    # quoted, kept for continuity and for the calibration scatter plot).
    pa_med, pa_scat = axial_stats(th_b)

    # Step 4 — the trace sub-population and its shared-axis test, which is
    # what the verdict actually rests on.
    tmask = trace_flags(a_b, b_b)
    n_trace = int(np.count_nonzero(tmask))
    if n_trace:
        t_pa, t_scat = axial_stats(th_b[tmask])
        trace_ab = float(np.median(ratio[tmask]))
        trace_a = float(np.median(a_b[tmask]))
    else:
        t_pa, t_scat, trace_ab, trace_a = None, None, None, None

    return FrameShape(
        n_detected=n_detected,
        n_sources=int(a_a.size),
        n_bright=int(a_b.size),
        median_ab=float(np.median(ratio)),
        max_ab=float(np.max(ratio)),
        pa_median=pa_med,
        pa_scatter=pa_scat,
        n_trace=n_trace,
        trace_frac=float(n_trace) / float(a_b.size),
        trace_ab=trace_ab,
        trace_a_px=trace_a,
        trace_pa=t_pa,
        trace_pa_scatter=t_scat,
        detect_sigma=detect_sigma,
        height=height, width=width)


# ---------------------------------------------------------------------------
# Pure judgement
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Verdict:
    """A classification plus the reason for it, in plain words.

    The reason string is not decoration.  A verdict that disagrees with a
    header label will be argued about, and the argument goes better when the
    row itself says which rule fired.
    """
    verdict: str
    strength_class: str
    reason: str


def classify_frame(shape: FrameShape) -> Verdict:
    """Decide dispersed / direct / indeterminate from a measured frame.

    The rules, in the order they are tried, each with the case it exists for:

    1. **Nothing to go on.**  Fewer than ``MIN_SOURCES`` usable detections —
       a cloud, a closed dome, a badly under-exposed frame.  There is no
       evidence either way, so we say so: ``indeterminate``.

    2. **Two or more traces sharing an axis** -> ``dispersed``.  This is the
       grism signature proper.  Cosmic rays fail it (random angles); a
       trailed exposure fails the trace gates themselves, which demand a/b
       >= 5 AND a >= 15 px — drift smears stars by a few pixels, not by
       fifteen at fifteen-to-one.

    3. **Two or more traces NOT sharing an axis** -> ``indeterminate``.
       Something is streaking this frame, but not a grating.  Satellite
       trails and cosmic-ray showers land here, and so would a genuinely
       broken exposure.  Refusing to call it is the honest answer.

    4. **One extreme trace** -> ``dispersed``.  The shared-axis test needs
       two traces and there is only one, so the eccentricity bar is raised to
       ``SOLO_TRACE_MIN_AB``.  This rule carries the bright-standard frames
       (Vega, Spica) where the target is the only thing in the field.

    5. **One modest trace** -> ``indeterminate``.  A single a/b = 6 object
       could be a faint spectrum, a satellite, or an edge-on galaxy.  One
       object is not a population.

    6. **No traces, round bright sources** -> ``direct``.  Ordinary imaging.

    7. **No traces, but the field is elongated** -> ``indeterminate``.  Not a
       grism, but not clean photometry either; the frame is trailed or
       astigmatic and the caller should know before using it.
    """
    if shape.n_sources < MIN_SOURCES:
        return Verdict(VERDICT_INDETERMINATE, "n/a",
                       "no usable sources extracted")

    if shape.n_trace >= MIN_TRACES_FOR_AXIS:
        scat = shape.trace_pa_scatter
        if scat is not None and scat <= TRACE_MAX_PA_SCATTER_DEG:
            return Verdict(
                VERDICT_DISPERSED, strength_of(shape),
                f"{shape.n_trace} traces share an axis "
                f"(PA scatter {scat:.1f} deg)")
        return Verdict(
            VERDICT_INDETERMINATE, "n/a",
            f"{shape.n_trace} elongated sources but no common axis "
            f"(PA scatter {scat:.1f} deg)" if scat is not None else
            f"{shape.n_trace} elongated sources, axis undefined")

    if shape.n_trace == 1:
        # Only one object is streaked.  With no second trace to corroborate
        # the axis, two things must hold: the streak must be extreme, AND
        # the field must be sparse.  The sparsity gate is what rejects a
        # satellite crossing a rich star field — see SOLO_MAX_SOURCES.
        solo_ab = shape.trace_ab if shape.trace_ab is not None else 0.0
        if solo_ab < SOLO_TRACE_MIN_AB:
            return Verdict(VERDICT_INDETERMINATE, "n/a",
                           f"one modest elongated source (a/b {solo_ab:.1f})")
        if shape.n_sources > SOLO_MAX_SOURCES:
            return Verdict(
                VERDICT_INDETERMINATE, "n/a",
                f"one streak (a/b {solo_ab:.0f}) in a rich field of "
                f"{shape.n_sources} round sources — satellite or artefact, "
                "not a grating")
        return Verdict(
            VERDICT_DISPERSED, strength_of(shape),
            f"solo trace, a/b {solo_ab:.0f}, sparse field "
            f"({shape.n_sources} sources) — isolated bright target")

    # No trace-like source at all.
    med = shape.median_ab if shape.median_ab is not None else 0.0
    if med <= DIRECT_MAX_MEDIAN_AB:
        return Verdict(VERDICT_DIRECT, "n/a",
                       f"no traces; bright sources round (median a/b "
                       f"{med:.2f})")
    return Verdict(VERDICT_INDETERMINATE, "n/a",
                   f"no traces but field elongated (median a/b {med:.2f}) "
                   "— trailed or astigmatic")


def trace_length_fraction(trace_a_px: Optional[float],
                          width: int) -> Optional[float]:
    """Trace semi-major sigma as a fraction of the frame's decoded width.

    Normalising by width is what makes one threshold work across the
    archive's three cameras and two binnings: a 2x-binned GSENSE frame is
    4,788 px wide and an iKon frame is 4,096, and the same physical
    dispersion lands on different pixel counts in each.  The frame width
    used here is the width of the array that actually DECODED — never the
    NAXIS1 header card, which in a tile-compressed file describes the
    compressed table's row length in bytes and not the image at all.
    """
    if trace_a_px is None or not width:
        return None
    return float(trace_a_px) / float(width)


def classify_strength(trace_a_px: Optional[float], width: int) -> str:
    """Split a dispersed frame into the high- or low-dispersion grism.

    The high-dispersion H-alpha unit spreads a given star roughly three
    times further than the low-dispersion broad-spectrum unit, and trace
    LENGTH is what carries that signal cleanly (see the constants above for
    why aspect ratio does not).  Frames landing between the two calibrated
    bounds are reported ``ambiguous`` rather than assigned: length still
    grows somewhat with source brightness and exposure, because a brighter
    trace stays above the extraction threshold further into its wings, and
    that spread is real rather than something a threshold can wish away.
    """
    frac = trace_length_fraction(trace_a_px, width)
    if frac is None:
        return "n/a"
    if frac >= STRENGTH_HIGH_MIN_FRAC:
        return "high"
    if frac <= STRENGTH_LOW_MAX_FRAC:
        return "low"
    return "ambiguous"


def strength_of(shape: FrameShape) -> str:
    """Convenience: strength class straight from a measured frame."""
    return classify_strength(shape.trace_a_px, shape.width)


def expected_verdict(filter_name: Optional[str]) -> Optional[str]:
    """The verdict a KNOWN label implies, or None when the label is the
    question rather than the answer.

    Used only for scoring the classifier against ground truth; nothing in
    the production path consults a label.
    """
    if not filter_name:
        return None
    if filter_name in KNOWN_DISPERSED_FILTERS:
        return VERDICT_DISPERSED
    if filter_name in KNOWN_DIRECT_FILTERS:
        return VERDICT_DIRECT
    return None


def expected_strength(filter_name: Optional[str]) -> Optional[str]:
    """The grism unit a known dispersed label implies, else None."""
    if not filter_name:
        return None
    if filter_name in HIGH_DISPERSION_FILTERS:
        return "high"
    if filter_name in LOW_DISPERSION_FILTERS:
        return "low"
    return None


# ---------------------------------------------------------------------------
# Impure edge: pixels in, FrameShape out
# ---------------------------------------------------------------------------
def extract_sources(data: np.ndarray,
                    sigma_ladder: Sequence[float] = DETECT_SIGMA_LADDER,
                    bright_n: int = BRIGHT_N) -> FrameShape:
    """Background-subtract, extract, and summarize one 2-D image array.

    The threshold ladder exists because the extractor carries a fixed budget
    of "pixels currently above threshold" (300k by default, raised here).  A
    nebulous or very crowded frame can blow that budget at 5 sigma; instead
    of losing the frame we step the threshold up and record which rung
    succeeded, so a reader can see that the frame was measured shallowly.
    """
    import sep

    # Raise the extractor's object-pixel budget once per process.  A single
    # grism trace can occupy a quarter-million pixels on its own, which the
    # stock 300k budget cannot hold alongside anything else.
    sep.set_extract_pixstack(int(3e6))

    # sep needs a contiguous native-byte-order float array.
    img = np.ascontiguousarray(data, dtype=np.float32)
    height, width = (img.shape if img.ndim == 2 else (0, 0))
    if img.ndim != 2 or img.size == 0:
        return empty_shape(sigma_ladder[0], 0, height, width)

    bkg = sep.Background(img, bw=BACK_SIZE, bh=BACK_SIZE,
                         fw=BACK_FILTER, fh=BACK_FILTER)
    img = img - bkg
    rms = float(bkg.globalrms)
    if not np.isfinite(rms) or rms <= 0:
        # A perfectly flat frame (all-zero, or a synthetic) has no noise
        # scale to threshold against; nothing can be measured from it.
        return empty_shape(sigma_ladder[0], 0, height, width)

    sources = None
    used_sigma = float(sigma_ladder[0])
    for sigma in sigma_ladder:
        try:
            # deblend_cont=1.0 disables deblending on purpose: a dispersed
            # trace is exactly the kind of extended, multi-peaked object the
            # deblender would shatter into fragments, destroying the very
            # elongation we came to measure.
            sources = sep.extract(img, float(sigma), err=rms,
                                  minarea=EXTRACT_MINAREA, deblend_cont=1.0)
            used_sigma = float(sigma)
            break
        except Exception:
            # Pixel-stack overflow (or any extractor refusal): step brighter.
            continue
    if sources is None or len(sources) == 0:
        return empty_shape(used_sigma, 0 if sources is None else len(sources),
                           height, width)

    return summarize_sources(
        a=sources["a"], b=sources["b"],
        theta_deg=np.degrees(sources["theta"]),
        flux=sources["flux"], npix=sources["npix"],
        detect_sigma=used_sigma, bright_n=bright_n,
        height=height, width=width)


def measure_file(path: str) -> FrameShape:
    """Open one archive frame and measure it.

    Packaging (plain / fpack / repackaged) is resolved by the shared grism
    reader, which inspects the HDU list instead of trusting a convention.
    That reader also side-steps the archive's tile-compression trap: in a
    tile-compressed FITS the on-disk NAXIS1/NAXIS2 describe the compressed
    BINTABLE, not the image, and only the Z-prefixed keywords give the real
    dimensions.  Reading pixels through astropy makes the question moot —
    the shape recorded here is the shape of the array that actually decoded.
    """
    import sys
    from pathlib import Path

    # The grism reader lives in the sibling package; make it importable
    # whether this module was imported from the repo root or from pipeline/.
    pkg_root = Path(__file__).resolve().parent.parent
    if str(pkg_root) not in sys.path:
        sys.path.insert(0, str(pkg_root))
    from macro_grism.fits_io import load_frame

    data, _header, _layout = load_frame(path)
    return extract_sources(data)
