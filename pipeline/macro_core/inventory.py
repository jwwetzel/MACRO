"""Pure S0b inventory logic: raw<->reduced linkage and calibration census.

Everything in this module is a *pure function* or a *data table*: no I/O, no
globals mutated, no hidden state — the same contract as ``macro_core.manifest``.
The build script (``pipeline/scripts/build_s0b_inventory.py``) wires these
functions to the S0 manifest database; the unit tests
(``pipeline/tests/test_inventory.py``) exercise every function on hand-built
cases, including the cases that MUST NOT work (a science frame classified as a
flat, a stem match jumping to a different night).

S0b answers two questions S0 deliberately left open:

1.  **Which reduced/ file came from which raw frame?**  S0 proved that
    (basename, JD) dedup cannot see the ``reduced/`` tree's *renamed* copies
    (the 29,737 (target, JD) collision pairs it handed us).  Exploration of
    the manifest shows the rename is almost always mechanical — the pipeline
    appends ``_calibrated`` (rarely ``_cal`` or ``_wcs``) before the
    extension and keeps the header JD — plus a small population whose JD was
    *rewritten* during reduction (15–70 s drifts, same filename).  The match
    ladder below encodes exactly those observations, strongest evidence
    first; anything the ladder cannot place is an orphan and is REPORTED,
    never hidden.

2.  **Which calibration frames exist for which camera era?**  Era identity
    comes from the S0 ``frames.era_id`` column (keyed on READOUTM, geometry,
    binning, EGAIN by ``manifest.era_key`` — S0b joins it, never recomputes
    it).  Calibration *kind* is normalized here from IMAGETYP where the
    header is explicit, and from a short list of observed filename
    conventions where it is not (master files written with IMAGETYP =
    'Light Frame'; iKon/grism twilight-flat series).

3.  **What physically changed on the telescope, and when?**  (Added by the
    plan review of 2026-10-03, findings F-3 / TE.F1 / DE.F5.)  An era is a
    HEADER history: it sees a readout-mode string change and misses a
    camera rotated by 180 degrees.  The ``mech_epoch`` layer built by the
    last section of this module is the HARDWARE history beneath the eras —
    camera identity, software flip state, filter-wheel map, and measured
    sky rotation — and it bounds which calibration frame may serve which
    science frame (:func:`calib_valid_for`).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from statistics import median
from typing import Iterable, Optional, Sequence

from . import manifest as _mf

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these,
# so changing a value here changes both the pipeline and its documentation).
# --------------------------------------------------------------------------

#: Two JDs closer than this are "the same instant" for matching purposes.
#: 1e-7 day = 8.6 ms — far below the 1 s granularity of the filename
#: timestamps, far above float noise on a JD read back from SQLite.
JD_EQUAL_TOL_DAYS = 1e-7

#: Dark-vs-science exposure-time match tolerance.  Header EXPTIME is a
#: driver-reported float ("16 s" is stored as 15.9999628067017), so exact
#: equality is wrong; the next *genuinely different* exposure in the archive
#: ladder differs by >= 25% (…, 8, 16, 32, …).  A relative window of 0.5%
#: (with a 0.02 s absolute floor for sub-second exposures) sits three orders
#: of magnitude above header fuzz and fifty below the ladder spacing, so it
#: can neither split a real match nor bridge two real exposure settings.
DARK_MATCH_REL_TOL = 0.005
DARK_MATCH_ABS_TOL = 0.02

#: Acquisition specs for the October shopping list: the minimum number of
#: RAW calibration frames per requirement before we call it covered.
#: Standard CCD-reduction practice (median-combine needs enough frames that
#: the master's noise is negligible against a single science frame):
#: >= 20 zero-second frames for a bias, >= 15 darks per exposure time,
#: >= 10 twilight flats per filter.
SPEC_N_BIAS = 20
SPEC_N_DARK = 15
SPEC_N_FLAT = 10

#: Header IMAGETYP values that explicitly declare a calibration kind.
#: ('FLAT' is the abbreviated form some master files carry.)
CALIB_KIND_OF_IMAGETYP: dict[str, str] = {
    "Bias Frame": "bias",
    "Dark Frame": "dark",
    "Flat Field": "flat",
    "FLAT": "flat",
}

#: Processing suffixes the reduction pipeline appends to a raw basename
#: (observed in the manifest: 27,122 ``_calibrated``, 48 ``_cal``, a handful
#: of ``_wcs``).  Lowercase; compared case-insensitively.
#: One vocabulary, defined once: S0's dedup and S0b's link ladder must agree
#: on what a processing suffix is, so this IS the manifest's set.
REDUCED_SUFFIX_TOKENS: frozenset[str] = _mf.PROCESSING_TOKENS

#: Filename patterns that mark a flat series written with IMAGETYP =
#: 'Light Frame' (or blank): the iKon / grism twilight-flat convention
#: ``<filter>.flat<n>…`` ('ha.flat3', 'OGG.flat2', 'HRG.flat6_2024-05-15…')
#: and the 1MHz commissioning series ('g.Flat-light-1MHz.20s.2').
_FLAT_SERIES_RE = re.compile(r"^[A-Za-z]{1,3}\.flat\d", re.IGNORECASE)
_FLAT_LIGHT_RE = re.compile(r"flat-light", re.IGNORECASE)

#: Science FILTER strings that collide with the calibration-kind vocabulary.
#: The archive contains a real case: era 76 holds 3 grism science exposures
#: of HD 6343 (2025-01-19, filenames say hrg/lrg, IMAGETYP = 'Light Frame')
#: whose FILTER card reads 'dark' — a filter-wheel/header glitch, not a real
#: filter.  A "flat in filter dark" is not an acquirable item, so filters on
#: this list are excluded from flat REQUIREMENTS on the shopping list; the
#: coverage matrix still records their cells (nothing hidden), and the
#: report states the exclusion with query-derived counts.
CALIB_VOCAB_FILTERS = frozenset({"bias", "dark", "flat"})

#: Target-key prefix that identifies the Dwarf survey family — the SAME rule
#: the S0 build uses for the ``__dw_survey__`` reconciliation selector
#: (``build_s0_manifest.build_project_counts``); duplicated as a named
#: constant so the two stages can never silently diverge (the S0b tests
#: assert the S0 build script still uses this prefix).
DW_SURVEY_PREFIX = "dw1"


# --------------------------------------------------------------------------
# Basename surgery: compression, extensions, processing suffixes
# --------------------------------------------------------------------------

#: ``strip_compression`` and ``frame_stem`` moved to ``macro_core.manifest``
#: on 2026-10-03, because S0's dedup now needs them (finding F-1) and S0 must
#: not import its own downstream stage.  They are re-exported here under
#: their original names so every existing caller keeps working.
strip_compression = _mf.strip_compression
frame_stem = _mf.frame_stem


def reduced_stem(basename: str) -> str:
    """Return the raw-frame stem a reduced filename points back to.

    Strips (repeatedly, innermost last) the processing suffixes the
    reduction pipeline appends: ``_calibrated``, ``_cal``, ``_wcs``, and a
    short copy counter that immediately FOLLOWS such a suffix
    (``…_calibrated_1``).  A bare trailing number is never stripped on its
    own — filenames like ``Vega_0p1s_hrg_7`` end in a legitimate frame
    index, and eating it would fabricate matches.

    ``mpg_NGC_7619_g_180s_2026-06-30T10-30-00_calibrated_1.fts`` →
    ``mpg_NGC_7619_g_180s_2026-06-30T10-30-00``.
    """
    tokens = frame_stem(basename).split("_")
    while len(tokens) > 1:
        last = tokens[-1].lower()
        if last in REDUCED_SUFFIX_TOKENS:
            tokens.pop()                      # the suffix itself
            continue
        if (last.isdigit() and len(last) <= 2 and len(tokens) > 2
                and tokens[-2].lower() in REDUCED_SUFFIX_TOKENS):
            tokens.pop()                      # copy counter AFTER a suffix
            continue
        break
    return "_".join(tokens)


# --------------------------------------------------------------------------
# Calibration-kind normalization
# --------------------------------------------------------------------------

def is_master(basename: str) -> bool:
    """True for stacked calibration *products* (``master…`` filenames).

    Masters are usable calibrations but cannot be re-derived or counted as
    raw acquisition frames, so the coverage matrix tallies them separately
    and the shopping-list specs count raw frames only.
    """
    return basename.lower().startswith("master")


def calib_kind(imagetyp: Optional[str], basename: str) -> Optional[str]:
    """Normalize one frame's calibration kind: ``bias``/``dark``/``flat``.

    Returns ``None`` for science (and anything unrecognized — S0b never
    guesses).  Precedence:

    1.  An explicit header IMAGETYP wins (:data:`CALIB_KIND_OF_IMAGETYP`);
        'Bias Frame' / 'Dark Frame' / 'Flat Field' / 'FLAT' cover 3,834
        catalog rows.
    2.  ``master…`` filenames: the archive holds 99 master flats written
        with IMAGETYP = 'Light Frame' (e.g. ``calib/2024_June/
        master_flat_g_1x1_Readout2_3s``).  Kind comes from the name, with
        **dark checked before flat** so a flat-field *dark*
        (``master-dark-flat-16-1x1``, ``master_flatdark_…``) is a dark, not
        a flat.
    3.  The iKon/grism twilight-flat series (``ha.flat3``,
        ``HRG.flat6_2024-05-15…``, ``g.Flat-light-1MHz…``), also written as
        'Light Frame'.

    Everything else — including the 'Fringe field N' *sky* exposures, which
    are Light frames of a fringe-calibration field, not detector
    calibrations — stays ``None``.
    """
    it = (imagetyp or "").strip()
    if it in CALIB_KIND_OF_IMAGETYP:
        return CALIB_KIND_OF_IMAGETYP[it]
    low = basename.lower()
    if is_master(basename):
        # Order matters: 'master-dark-flat-…' is a DARK for flat exposures.
        if "bias" in low:
            return "bias"
        if "dark" in low:
            return "dark"
        if "flat" in low:
            return "flat"
        return None
    if _FLAT_SERIES_RE.match(low) or _FLAT_LIGHT_RE.search(low):
        return "flat"
    return None


def is_science(imagetyp: Optional[str], kind: Optional[str]) -> bool:
    """True when a frame counts as science for the coverage matrix.

    Science = not classified as a calibration, AND the header either says
    'Light Frame' or says nothing at all.  The blank-IMAGETYP clause is a
    deliberate, load-bearing decision: the 2026-06-28 → 2026-07-02 nights
    (the newest frames in the archive, the CURRENT camera) were written
    without an IMAGETYP card, and their filenames/filters are plainly
    science.  Excluding them would hide precisely the era whose calibration
    gaps the October run must fill.
    """
    if kind is not None:
        return False
    it = (imagetyp or "").strip()
    return it == "" or it.startswith("Light")


def is_calib_vocab_filter(filter_name: Optional[str]) -> bool:
    """True when a science FILTER string collides with the calibration
    vocabulary (:data:`CALIB_VOCAB_FILTERS`) — a header glitch, never a
    physical filter.

    Such strings must not spawn flat requirements ("flat dark x >=10" is a
    physically meaningless acquisition item that could confuse the October
    ops request if read literally).  Comparison is case-insensitive and
    whitespace-stripped; ``None`` (no filter card) is NOT a collision — a
    blank filter is a separate, honest fact.
    """
    if filter_name is None:
        return False
    return str(filter_name).strip().lower() in CALIB_VOCAB_FILTERS


# --------------------------------------------------------------------------
# Dark exposure-time matching and binning
# --------------------------------------------------------------------------

def exptime_bin(exptime: Optional[float]) -> Optional[float]:
    """Canonical exposure-time bin: the value rounded to 3 significant digits.

    Groups the driver's float fuzz onto one label (15.9999628 → 16.0,
    0.09998 → 0.1, 2.0000159 → 2.0) while keeping every genuinely distinct
    exposure setting in the archive separate (adjacent settings differ by
    >= 25%).  Missing exposure times bin to ``None``; non-positive ones to
    0.0 (bias-like).
    """
    if exptime is None or (isinstance(exptime, float) and math.isnan(exptime)):
        return None
    x = float(exptime)
    if x <= 0:
        return 0.0
    digits = 2 - math.floor(math.log10(abs(x)))
    return round(x, digits)


def dark_matches(dark_exptime: Optional[float],
                 science_exptime: Optional[float]) -> bool:
    """True when a dark's exposure time serves a science exposure time.

    Exact-match policy with the documented tolerance
    (:data:`DARK_MATCH_REL_TOL` relative, :data:`DARK_MATCH_ABS_TOL`
    absolute floor): |dark − science| <= max(0.02 s, 0.5% · science).
    """
    if dark_exptime is None or science_exptime is None:
        return False
    d, s = float(dark_exptime), float(science_exptime)
    if math.isnan(d) or math.isnan(s):
        return False
    return abs(d - s) <= max(DARK_MATCH_ABS_TOL, DARK_MATCH_REL_TOL * abs(s))


# --------------------------------------------------------------------------
# The raw<->reduced match ladder
# --------------------------------------------------------------------------
# Exact copies (same basename AND same JD) are already grouped by S0's
# dup_group and never reach this ladder — the build links them directly with
# method 'same_basename_jd'.  The ladder below places every OTHER reduced
# frame, strongest evidence first:
#
#   stem_jd            same JD and the filename is the raw one plus a known
#                      processing suffix — the rename S0 predicted.
#   stem_jd_drift      same stem, same NIGHT, but the JD was rewritten
#                      during reduction (observed 15–70 s drifts).  The
#                      night gate stops a drift match from jumping between
#                      two different visits of the same field.
#   target_jd          same (target, JD) — S0's collision-pair evidence —
#                      when the filename gives nothing (fully renamed).
#   target_jd_ambiguous  same (target, JD) matches SEVERAL raw frames
#                      (burst sequences sharing a start second); every
#                      candidate pair is recorded, none invented.
#   (orphan)           the ladder returns [] and the build records the
#                      reduced frame with no raw parent — stacks, class
#                      products, frames whose raw copy never reached the
#                      archive.  Characterized in the report, never hidden.

def link_reduced(stem: str, jd: Optional[float], night: Optional[str],
                 target_key: Optional[str],
                 raw_by_stem: dict[str, Sequence[tuple]],
                 raw_by_target_jd: dict[tuple, Sequence[tuple]],
                 ) -> list[tuple[int, str, Optional[float]]]:
    """Place one reduced frame on the match ladder.

    Parameters
    ----------
    stem
        :func:`reduced_stem` of the reduced file's basename.
    jd, night, target_key
        The reduced frame's own header JD, S0 night label, and alias key.
    raw_by_stem
        ``frame_stem`` → sequence of ``(raw_id, jd, night)`` over canonical
        non-reduced frames.
    raw_by_target_jd
        ``(target_key, round(jd, 7))`` → sequence of ``(raw_id, jd, night)``
        over the same frames.

    Returns
    -------
    list of (raw_id, match_method, jd_drift_seconds)
        One entry per matched raw parent (several only for
        ``target_jd_ambiguous``); empty list = orphan.  Selection within a
        rung is deterministic: smallest raw_id wins ties.
    """
    cands = raw_by_stem.get(stem, ())

    # Rung 1: stem + identical JD — the mechanical '_calibrated' rename.
    if jd is not None:
        exact = [c for c in cands
                 if c[1] is not None and abs(c[1] - jd) <= JD_EQUAL_TOL_DAYS]
        if exact:
            best = min(exact, key=lambda c: c[0])
            return [(best[0], "stem_jd", 0.0)]

    # Rung 2: stem + same night, JD rewritten by the reduction pipeline.
    if night is not None:
        same_night = [c for c in cands if c[2] == night and c[1] is not None
                      and jd is not None]
        if same_night:
            best = min(same_night, key=lambda c: (abs(c[1] - jd), c[0]))
            drift_s = (jd - best[1]) * 86400.0
            return [(best[0], "stem_jd_drift", drift_s)]

    # Rung 3: (target, JD) — the S0 collision-pair evidence.
    if target_key is not None and jd is not None:
        tj = raw_by_target_jd.get((target_key, round(jd, 7)), ())
        if len(tj) == 1:
            return [(tj[0][0], "target_jd", 0.0)]
        if len(tj) > 1:
            return [(c[0], "target_jd_ambiguous", 0.0)
                    for c in sorted(tj, key=lambda c: c[0])]

    # Off the ladder: orphan.
    return []


# --------------------------------------------------------------------------
# Requirement arithmetic (shared by coverage matrix and shopping list)
# --------------------------------------------------------------------------

def coverage_status(n_raw_calib: int, n_master: int, spec_n: int) -> str:
    """Classify one coverage cell.

    * ``ok``          — enough RAW frames to (re)build a master to spec.
    * ``partial``     — some raw frames, fewer than spec.
    * ``master_only`` — no raw frames, but a stacked master product exists
      (usable, not re-derivable).
    * ``missing``     — nothing at all.
    """
    if n_raw_calib >= spec_n:
        return "ok"
    if n_raw_calib > 0:
        return "partial"
    if n_master > 0:
        return "master_only"
    return "missing"


def gap_spec(need_kind: str, req_key: Optional[str], n_raw_calib: int,
             spec_n: int) -> str:
    """Human-readable acquisition spec for one shopping-list row.

    e.g. ``dark 240s x >=15 (have 0)``; ``flat g x >=10 (have 3)``;
    ``bias x >=20 (have 0)``.
    """
    middle = f" {req_key}" if req_key else ""
    return f"{need_kind}{middle} x >={spec_n} (have {n_raw_calib})"


def fmt_exptime(x: Optional[float]) -> str:
    """Compact exposure-time label for specs and report tables (240 → '240s',
    0.1 → '0.1s', None → '?s')."""
    if x is None:
        return "?s"
    return f"{x:g}s"


def projects_of_target(target_key: Optional[str],
                       project_of_key: dict[str, frozenset[str]],
                       dw_projects: frozenset[str] = frozenset(),
                       ) -> frozenset[str]:
    """Projects that claim one target key, per the S0 project_counts lists.

    ``project_of_key`` maps explicit target keys to project-name sets;
    ``dw_projects`` is the set claiming the Dwarf survey family, applied to
    every key with the :data:`DW_SURVEY_PREFIX` prefix (the S0 selector).
    """
    if target_key is None:
        return frozenset()
    out = set(project_of_key.get(target_key, frozenset()))
    if target_key.startswith(DW_SURVEY_PREFIX):
        out |= dw_projects
    return frozenset(out)


# ==========================================================================
# Mechanical epochs — the hardware history beneath the eras  (F-3, TE.F1)
# ==========================================================================
# THE DEFECT.  The era registry keys on four header values (READOUTM, NAXIS,
# XBINNING, EGAIN).  The telescope has had four cameras and a dozen
# mechanical states, and the registry sees neither number: era 76 runs
# straight through the 2025 monsoon shutdown, across which the sky rotated
# by 179.3 degrees on the detector.  A flat field, a dust-donut map, a
# grism trace prior — all are valid for one orientation only — and nothing
# in the manifest could say which side of the boundary a frame was on.
#
# THE EVIDENCE USED, strongest first:
#
#   camera      which camera took the frame (manifest.camera_id) — exact.
#   flipstat    MaxIm's software flip state, from the header — exact.
#   wheel map   the filter wheel's slot list (FWALLNAM) — exact.
#   rotation    the position angle S1 measured on plate-solved frames.
#
# THE ROTATION TRAP, and how it is avoided.  The measured position angle is
# NOT a pure hardware number: it carries a term that depends on where the
# telescope points (polar misalignment turns the field by an angle that
# grows as 1/cos(dec) and varies with hour angle).  Measured on the ASI
# camera's first season, with nothing touched: +0.03 deg at dec +56,
# +0.45 at dec +25, +0.70 at dec -19 — a 0.7-degree spread from pointing
# alone.  A rule "nightly median moved by more than 0.3 deg" therefore
# fires on nearly every change of target list (it finds 30 "re-seats" in
# 2025 where the hardware log has none).  Two tests replace it:
#
#   FINE    compare the SAME TARGET across nights.  Pointing cancels; what
#           is left is the hardware.  A step is declared when the median
#           same-target difference exceeds MECH_ROT_STEP_DEG.
#   COARSE  compare nightly medians, but only against a threshold
#           (MECH_ROT_COARSE_DEG) set well above the measured pointing
#           spread.  It needs no shared target, so it still sees a flipped
#           or grossly rotated camera on a night of new fields.
#
# A step between 0.3 and 1.5 degrees on a night that shares no target with
# the nights before it is INVISIBLE to both tests.  That is stated, not
# hidden: every epoch records how many of its night-to-night transitions
# each test could actually examine.

#: A same-target rotation step larger than this starts a new mechanical
#: epoch.  The value is the chair's (SYNTHESIS F-3, from TE.F1); its
#: false-alarm rate on this archive is MEASURED by the null distribution
#: :func:`segment_mech_epochs` returns, not assumed.
MECH_ROT_STEP_DEG = 0.3

#: A nightly-median rotation change larger than this starts a new epoch
#: even with no shared target.  Must exceed the pointing-induced spread
#: (0.7 deg measured, see above) with margin.
MECH_ROT_COARSE_DEG = 1.5

#: Rotation samples beyond this |declination| are ignored: the pointing
#: term scales as 1/cos(dec) and the one dec +85 field in the archive
#: (NGC 188) sits 1.9 deg from its night's other fields.
MECH_ROT_MAX_ABS_DEC = 65.0

#: Fewest plate-solved frames for a night (or a target on a night) to count
#: as rotation evidence.  One or two solves are an anecdote; a wrong
#: astrometric match is rare but not absent.
MECH_ROT_MIN_FRAMES = 3

#: Two nights' solves of "the same target" must also agree on the sky to
#: this tolerance (great-circle degrees).  Target names like 'stars' or
#: 'guide field' are re-used for different fields; the name alone is not
#: evidence that the pointing term cancels.
MECH_TARGET_MATCH_DEG = 0.5

#: Dark/bias set-point tolerance, deg C (DE.F5: "set-point +/- 2 C").
SET_TEMP_TOL_C = 2.0

#: Slot names that are two SPELLINGS of one element.  The wheel map is
#: written by whichever program took the frame, and the two programs name
#: the grisms differently: MaxIm's wheel labels say 'OGGrism' / 'HaGrism',
#: pyscope's say 'lrg' / 'hrg' (the same pairing staging.GRISM_ALL records).
#: In January 2025 the two alternated night by night on an untouched wheel;
#: without this table every such night would read as a reloaded wheel.
#: Applied after case-folding.  Keep it to names PROVEN to be one element.
WHEEL_SLOT_SYNONYMS: dict[str, str] = {
    "oggrism": "lrg",
    "hagrism": "hrg",
}


def wheel_compare_key(wheel_map: Optional[str]) -> Optional[tuple[str, ...]]:
    """A wheel map reduced to what a PHYSICAL change would alter.

    ``wheel_map`` is the ``|``-joined slot list S0 stores in
    ``frames.fwallnam``.  The key is the slot tuple, case-folded and passed
    through :data:`WHEEL_SLOT_SYNONYMS`, so a map re-spelled by different
    software compares equal and a slot holding a different element does
    not.  ``None`` for a missing map.
    """
    if wheel_map is None or not str(wheel_map).strip():
        return None
    return tuple(WHEEL_SLOT_SYNONYMS.get(tok.strip().lower(),
                                         tok.strip().lower())
                 for tok in str(wheel_map).split("|"))


#: Boundary causes, as written into ``mech_epoch.boundary_cause``.
CAUSE_CAMERA = "camera"                  # first appearance of this camera
CAUSE_CAMERA_SWAP = "camera_swap"        # another camera was mounted between
CAUSE_FLIPSTAT = "flipstat"              # software flip state changed
CAUSE_WHEEL = "wheel_map"                # a wheel slot's content changed
CAUSE_ROT_FINE = "rotation_fine"         # same-target step > 0.3 deg
CAUSE_ROT_COARSE = "rotation_coarse"     # nightly-median step > 1.5 deg

#: Causes that change WHERE A PIXEL LANDS IN THE FILE.  A new camera is a
#: new sensor; a flip-state change mirrors the saved array.  Rotating or
#: re-seating a camera, or reloading the wheel, moves the sky and the dust
#: but leaves every hot pixel and the bias structure at the same array
#: coordinates.  Darks and biases are therefore bounded by these causes
#: only; flats by every cause.  (TE.F1: "a hot-pixel map showing the same
#: pixel coordinates before and after … would leave darks valid (flats
#: still not)".)
DETECTOR_CAUSES: frozenset[str] = frozenset(
    {CAUSE_CAMERA, CAUSE_CAMERA_SWAP, CAUSE_FLIPSTAT})


def circ_diff(a: float, b: float) -> float:
    """Signed smallest difference ``a - b`` between two angles, degrees,
    in (-180, 180].  359.9 vs 0.1 is -0.2, not 359.8."""
    d = (float(a) - float(b) + 180.0) % 360.0 - 180.0
    return 180.0 if d == -180.0 else d


def circular_median(angles: Iterable[float]) -> tuple[float, float]:
    """Robust centre and scatter of a set of angles: ``(median, MAD)``.

    A plain median fails across the 0/360 seam (the ASI's first season
    sits at +0.4 deg with some frames at 359.9).  The angles are first
    referred to their circular MEAN direction, where the seam is as far
    away as it can be; the median and the median absolute deviation are
    taken there and the median is folded back into [0, 360).

    The MAD is returned unscaled (multiply by 1.4826 for a Gaussian sigma).
    """
    vals = [float(a) for a in angles]
    if not vals:
        raise ValueError("circular_median needs at least one angle")
    sx = sum(math.cos(math.radians(a)) for a in vals)
    sy = sum(math.sin(math.radians(a)) for a in vals)
    centre = math.degrees(math.atan2(sy, sx))
    offs = [circ_diff(a, centre) for a in vals]
    med = median(offs)
    mad = median(abs(o - med) for o in offs)
    return (centre + med) % 360.0, mad


@dataclass(frozen=True)
class RotSample:
    """One plate-solved frame's contribution: where it pointed and the
    position angle S1 measured there."""
    target: Optional[str]
    rot_deg: float
    ra_deg: float
    dec_deg: float


@dataclass(frozen=True)
class NightState:
    """Everything known about one camera on one night.

    ``flipstat`` / ``wheel_map`` are the night's modal header values, or
    ``None`` when no frame of the night carries the card (the whole AC4040
    era has no wheel map).  ``rot`` holds the night's plate-solved frames.

    ``wheel_maps`` holds the night's modal map FOR EACH SLOT COUNT seen
    that night.  The telescope reports its wheels in more than one way —
    in November 2024 single nights carry both a 16-slot 'Dual Wheels' map
    and a 5-slot 'FLI' map — and a 16-slot list cannot be compared with a
    5-slot one: they are different descriptions, not different states.
    Maps are therefore only ever compared with the last map OF THE SAME
    LENGTH.  When ``wheel_maps`` is empty, ``wheel_map`` alone is used.
    """
    camera: str
    night: str
    n_frames: int
    flipstat: Optional[str] = None
    wheel_map: Optional[str] = None
    swcreate: Optional[str] = None
    telpier: Optional[str] = None
    rot: tuple[RotSample, ...] = ()
    wheel_maps: tuple[str, ...] = ()


@dataclass(frozen=True)
class _TargetRot:
    """A target's rotation on one night: median, standard error, position."""
    rot: float
    se: float
    n: int
    ra: float
    dec: float


def _usable(samples: Iterable[RotSample]) -> list[RotSample]:
    """Rotation samples inside the declination limit."""
    return [s for s in samples if abs(s.dec_deg) <= MECH_ROT_MAX_ABS_DEC]


def night_rotation(state: NightState) -> Optional[tuple[float, float, int]]:
    """The night's rotation level: ``(median, MAD, n)`` over every usable
    solve, or ``None`` with fewer than :data:`MECH_ROT_MIN_FRAMES`."""
    use = _usable(state.rot)
    if len(use) < MECH_ROT_MIN_FRAMES:
        return None
    med, mad = circular_median(s.rot_deg for s in use)
    return med, mad, len(use)


def target_rotations(state: NightState) -> dict[str, _TargetRot]:
    """Per-target rotation on one night, for the same-target (fine) test.

    Only NAMED targets with at least :data:`MECH_ROT_MIN_FRAMES` usable
    solves qualify.  The standard error is 1.4826 MAD / sqrt(n) — the
    scatter of the solves themselves, which is all one night can know.
    """
    groups: dict[str, list[RotSample]] = {}
    for s in _usable(state.rot):
        if s.target:
            groups.setdefault(s.target, []).append(s)
    out: dict[str, _TargetRot] = {}
    for tgt, grp in groups.items():
        if len(grp) < MECH_ROT_MIN_FRAMES:
            continue
        med, mad = circular_median(g.rot_deg for g in grp)
        ra, dec = _mf.median_radec([g.ra_deg for g in grp],
                                   [g.dec_deg for g in grp])
        out[tgt] = _TargetRot(rot=med, se=1.4826 * mad / math.sqrt(len(grp)),
                              n=len(grp), ra=ra, dec=dec)
    return out


def mech_epoch_id(camera: str, first_night: str) -> str:
    """Identifier of a mechanical epoch: ``<camera>:<first night>``.

    A NAME, not a counter, on purpose.  Era ids are a pinned registry
    because a counter renumbers when a configuration is discovered
    mid-timeline; an id made of the camera and the night the state began is
    stable under every later ingest and says what it is.
    """
    return f"{camera}:{first_night}"


def segment_mech_epochs(states: Sequence[NightState]
                        ) -> tuple[list[dict], list[dict], list[dict]]:
    """Cut the telescope's history into mechanical epochs.

    Parameters
    ----------
    states
        One :class:`NightState` per (camera, night) that holds any frame,
        in any order.

    Returns
    -------
    epochs : list of dict
        One row per mechanical epoch, chronological — the ``mech_epoch``
        table.  ``boundary_cause`` says why the epoch BEGAN; ``step_deg`` is
        the measured rotation step into it (same-target where a shared
        target existed, nightly-median otherwise) with ``step_err_deg``;
        ``n_transitions`` / ``n_tested_fine`` / ``n_tested_coarse`` say how
        many of the night-to-night transitions INSIDE the epoch each
        rotation test was able to examine.
    nights : list of dict
        One row per (camera, night) — the ``night_mech_epoch`` table.
        ``certain = 0`` marks a night that carries no rotation evidence and
        lies between the last evidenced night of one epoch and a boundary
        found by rotation alone: the hardware changed somewhere in that
        gap and the night cannot be placed on either side.  It is filed
        under the EARLIER epoch and no flat may be matched to it.
    null : list of dict
        Every same-target, night-to-night rotation difference measured
        INSIDE an epoch.  This is the empirical distribution of the fine
        statistic when nothing changed — the false-alarm test of
        :data:`MECH_ROT_STEP_DEG`.

    The walk is per camera, in night order.  A new epoch begins when

    * the camera is seen for the first time, or returns after another
      camera was mounted in between (``camera`` / ``camera_swap``);
    * the software flip state changes (``flipstat``);
    * a filter-wheel slot's content changes (``wheel_map``): the night's
      map differs, after :func:`wheel_compare_key`, from the last map of
      the SAME slot count.  A map of a different length is a different
      description of the wheels and is never compared slot by slot — but
      the FIRST appearance, on a camera, of a slot count it has never
      reported before is itself a boundary: the wheel went from 7 named
      slots to 9 on 2026-06-28, and no slot-by-slot comparison could see
      that.  A slot count that merely RE-appears (the 16-slot and 5-slot
      descriptions alternated through November 2024) is not;
    * the same-target rotation step exceeds :data:`MECH_ROT_STEP_DEG`
      (``rotation_fine``), or the nightly median moves by more than
      :data:`MECH_ROT_COARSE_DEG` (``rotation_coarse``).

    A header value of ``None`` (card absent that night) never triggers a
    boundary and never overwrites the epoch's last known value.
    """
    by_cam: dict[str, list[NightState]] = {}
    for st in states:
        by_cam.setdefault(st.camera, []).append(st)
    # Every (night, camera) pair, to detect a different camera in between.
    cam_nights = sorted((st.night, st.camera) for st in states)

    def _other_camera_between(cam: str, lo: str, hi: str) -> bool:
        return any(lo < n < hi and c != cam for n, c in cam_nights)

    epochs: list[dict] = []
    nights: list[dict] = []
    null: list[dict] = []

    for cam in sorted(by_cam):
        seq = sorted(by_cam[cam], key=lambda s: s.night)
        cur: Optional[dict] = None           # the open epoch's working state
        seen_lens: set[int] = set()          # slot counts this camera has shown
        det_first: Optional[str] = None      # first night of the detector run
        prev_night: Optional[str] = None

        def _open(st: NightState, causes: list[str], step, step_err,
                  step_n, step_basis) -> dict:
            return {
                "camera": cam, "first_night": st.night, "causes": causes,
                "step": step, "step_err": step_err, "step_n": step_n,
                "step_basis": step_basis,
                "flipstat": None, "wheel_map": None,
                "wheel_by_len": {},      # slot count -> last compare key
                "levels": [],            # nightly median rotations
                "targets": {},           # target -> list of _TargetRot
                "night_rows": [],        # indices into ``nights``
                "n_frames": 0, "n_rot_frames": 0,
                "n_transitions": 0, "n_fine": 0, "n_coarse": 0,
                "swcreate": {}, "telpier": {},
            }

        def _close(ep: dict, last_night: str) -> None:
            level = circular_median(ep["levels"]) if ep["levels"] else None
            modal = lambda d: max(sorted(d), key=d.get) if d else None
            epochs.append({
                "mech_epoch": mech_epoch_id(cam, ep["first_night"]),
                "camera": cam,
                "first_night": ep["first_night"], "last_night": last_night,
                "n_nights": len(ep["night_rows"]),
                "n_frames": ep["n_frames"],
                "boundary_cause": ",".join(ep["causes"]),
                "detector_epoch": mech_epoch_id(cam, ep["det_first"]),
                "rotation_deg": None if level is None else level[0],
                "rotation_mad_deg": None if level is None else level[1],
                "n_rot_nights": len(ep["levels"]),
                "n_rot_frames": ep["n_rot_frames"],
                "step_deg": ep["step"], "step_err_deg": ep["step_err"],
                "step_n_targets": ep["step_n"],
                "step_basis": ep["step_basis"],
                "flipstat": ep["flipstat"], "wheel_map": ep["wheel_map"],
                "swcreate": modal(ep["swcreate"]),
                "telpier": modal(ep["telpier"]),
                "n_transitions": ep["n_transitions"],
                "n_tested_fine": ep["n_fine"],
                "n_tested_coarse": ep["n_coarse"],
                "n_gap_nights": sum(1 for i in ep["night_rows"]
                                    if not nights[i]["certain"]),
            })

        for st in seq:
            nrot = night_rotation(st)
            trot = target_rotations(st)
            causes: list[str] = []
            step = step_err = step_basis = None
            step_n = 0
            tested_fine = tested_coarse = False
            deltas: list[tuple[str, float, float, float]] = []

            if cur is None:
                causes.append(CAUSE_CAMERA)
            else:
                if _other_camera_between(cam, prev_night, st.night):
                    causes.append(CAUSE_CAMERA_SWAP)
                if (st.flipstat is not None and cur["flipstat"] is not None
                        and st.flipstat != cur["flipstat"]):
                    causes.append(CAUSE_FLIPSTAT)
                for wm in (st.wheel_maps or (st.wheel_map,)):
                    key = wheel_compare_key(wm)
                    if key is None:
                        continue
                    last = cur["wheel_by_len"].get(len(key))
                    if (last is not None and last != key) or (
                            seen_lens and len(key) not in seen_lens):
                        causes.append(CAUSE_WHEEL)
                        break
                # ---- fine test: the same target, before and after ---------
                for tgt, now in trot.items():
                    hist = cur["targets"].get(tgt)
                    if not hist:
                        continue
                    ra0, dec0 = _mf.median_radec([h.ra for h in hist],
                                                 [h.dec for h in hist])
                    if _mf.angular_separation_deg(
                            now.ra, now.dec, ra0, dec0) > MECH_TARGET_MATCH_DEG:
                        continue        # same name, different field
                    ref, ref_mad = circular_median(h.rot for h in hist)
                    ref_se = (1.4826 * ref_mad / math.sqrt(len(hist))
                              if len(hist) > 1 else hist[0].se)
                    deltas.append((tgt, circ_diff(now.rot, ref),
                                   math.hypot(now.se, ref_se), now.dec))
                if deltas:
                    tested_fine = True
                    ds = [d[1] for d in deltas]
                    step = median(ds)
                    step_n = len(ds)
                    scatter = (1.4826 * median(abs(d - step) for d in ds)
                               / math.sqrt(len(ds))) if len(ds) > 1 else 0.0
                    formal = math.sqrt(sum(d[2] ** 2 for d in deltas)) \
                        / len(deltas)
                    step_err = max(scatter, formal)
                    step_basis = "same_target"
                    if abs(step) > MECH_ROT_STEP_DEG:
                        causes.append(CAUSE_ROT_FINE)
                # ---- coarse test: nightly medians, generous threshold -----
                if nrot is not None and cur["levels"]:
                    tested_coarse = True
                    level, level_mad = circular_median(cur["levels"])
                    d = circ_diff(nrot[0], level)
                    if step is None:
                        step, step_basis = d, "nightly_median"
                        step_err = math.hypot(
                            1.4826 * nrot[1] / math.sqrt(nrot[2]),
                            1.4826 * level_mad
                            / math.sqrt(len(cur["levels"])))
                    if abs(d) > MECH_ROT_COARSE_DEG:
                        causes.append(CAUSE_ROT_COARSE)

            if causes:
                if cur is not None:
                    # A boundary found by ROTATION ALONE is dated only to
                    # "after the last night that showed the old angle": the
                    # trailing nights of the old epoch with no rotation
                    # evidence could belong to either side.
                    if all(c in (CAUSE_ROT_FINE, CAUSE_ROT_COARSE)
                           for c in causes):
                        for i in reversed(cur["night_rows"]):
                            if nights[i]["n_rot"] >= MECH_ROT_MIN_FRAMES:
                                break
                            nights[i]["certain"] = 0
                            nights[i]["basis"] = "gap"
                    _close(cur, prev_night)
                if det_first is None or any(c in DETECTOR_CAUSES
                                            for c in causes):
                    det_first = st.night
                cur = _open(st, causes, step, step_err, step_n, step_basis)
                cur["det_first"] = det_first
            else:
                cur["n_transitions"] += 1
                cur["n_fine"] += tested_fine
                cur["n_coarse"] += tested_coarse
                # Inside an epoch every same-target difference is a draw
                # from the "nothing changed" distribution.
                for tgt, d, se, dec in deltas:
                    null.append({"camera": cam, "night": st.night,
                                 "mech_epoch": mech_epoch_id(
                                     cam, cur["first_night"]),
                                 "target": tgt, "delta_deg": d,
                                 "delta_se_deg": se, "dec_deg": dec})

            # ---- fold this night into the open epoch -----------------------
            if st.flipstat is not None:
                cur["flipstat"] = st.flipstat
            if st.wheel_map is not None:
                cur["wheel_map"] = st.wheel_map
            for wm in (st.wheel_maps or (st.wheel_map,)):
                key = wheel_compare_key(wm)
                if key is not None:
                    cur["wheel_by_len"][len(key)] = key
                    seen_lens.add(len(key))
            if nrot is not None:
                cur["levels"].append(nrot[0])
                cur["n_rot_frames"] += nrot[2]
            for tgt, tr in trot.items():
                cur["targets"].setdefault(tgt, []).append(tr)
            cur["n_frames"] += st.n_frames
            for bucket, val in (("swcreate", st.swcreate),
                                ("telpier", st.telpier)):
                if val:
                    cur[bucket][val] = cur[bucket].get(val, 0) + st.n_frames
            cur["night_rows"].append(len(nights))
            nights.append({
                "camera": cam, "night": st.night,
                "mech_epoch": mech_epoch_id(cam, cur["first_night"]),
                "detector_epoch": mech_epoch_id(cam, cur["det_first"]),
                "certain": 1,
                "basis": "rotation" if nrot is not None else "header",
                "n_frames": st.n_frames,
                "rot_median_deg": None if nrot is None else nrot[0],
                "rot_mad_deg": None if nrot is None else nrot[1],
                "n_rot": 0 if nrot is None else nrot[2],
                "flipstat": st.flipstat, "wheel_map": st.wheel_map,
            })
            prev_night = st.night
        if cur is not None:
            _close(cur, prev_night)

    epochs.sort(key=lambda e: (e["first_night"], e["camera"]))
    for i, ep in enumerate(epochs, 1):
        ep["seq"] = i
    return epochs, nights, null


# --------------------------------------------------------------------------
# Calibration validity across mechanical boundaries
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class EpochTag:
    """Where one frame sits in the hardware history."""
    mech_epoch: Optional[str]
    detector_epoch: Optional[str]
    certain: bool = True


def calib_valid_for(kind: str, science: EpochTag, calib: EpochTag) -> bool:
    """May a calibration frame of ``kind`` be applied to a science frame?

    THE RULE THAT F-3 EXISTS TO ENFORCE: **no calibration is applied across
    a boundary that invalidates it.**

    * a ``flat`` records the illumination pattern — vignetting, dust,
      filter — in detector coordinates.  Any re-seat, rotation, flip or
      wheel reload moves it.  A flat is valid only inside the science
      frame's own ``mech_epoch``, and only when BOTH frames are placed with
      certainty (an uncertain-gap night cannot be shown not to cross).
    * a ``bias`` or ``dark`` records the sensor.  It survives a rotation
      or a wheel reload and dies with a camera change or a flip-state
      change: valid inside the same ``detector_epoch``.

    A frame with no epoch at all (unknown camera) is valid for nothing and
    nothing is valid for it: unknown is never read as "same".
    """
    if kind == "flat":
        return (science.mech_epoch is not None
                and science.mech_epoch == calib.mech_epoch
                and bool(science.certain) and bool(calib.certain))
    if kind in ("bias", "dark"):
        return (science.detector_epoch is not None
                and science.detector_epoch == calib.detector_epoch)
    raise ValueError(f"unknown calibration kind {kind!r}")


#: Outcomes of :func:`settings_match`.
SETTINGS_MATCH = "match"
SETTINGS_MISMATCH = "mismatch"
SETTINGS_UNVERIFIED = "unverified"


def settings_match(sci_gain: Optional[str], sci_offset: Optional[float],
                   sci_temp: Optional[float], cal_gain: Optional[str],
                   cal_offset: Optional[float], cal_temp: Optional[float],
                   ) -> str:
    """Do a dark/bias frame's camera settings match a science frame's?

    DE.F5: calibration matching must key on gain setting, offset setting
    and cooler set-point (within :data:`SET_TEMP_TOL_C`) — era 76 pools
    Mode0 at -10 C with Mode0 at 0 C, and a -10 C dark does not describe a
    0 C frame.

    Three outcomes, because "the header does not say" is not "the same":

    * ``mismatch``   — a setting known on BOTH frames differs;
    * ``unverified`` — nothing contradicts, but a setting the science
      frame records is missing from the calibration frame (stacked masters
      often lose their cards).  Counted separately; never as a match;
    * ``match``      — every setting the science frame records is present
      on the calibration frame and agrees.

    A setting the SCIENCE frame does not record constrains nothing (the
    AC4040 headers carry no GAIN or OFFSET card at all; its gain is in the
    era key).
    """
    def _blank(v) -> bool:
        return v is None or (isinstance(v, float) and math.isnan(v)) \
            or (isinstance(v, str) and not v.strip())

    unverified = False
    for sci, cal, numeric_tol in ((sci_gain, cal_gain, None),
                                  (sci_offset, cal_offset, 0.0),
                                  (sci_temp, cal_temp, SET_TEMP_TOL_C)):
        if _blank(sci):
            continue
        if _blank(cal):
            unverified = True
            continue
        if numeric_tol is None:
            if str(sci).strip() != str(cal).strip():
                return SETTINGS_MISMATCH
        elif abs(float(sci) - float(cal)) > numeric_tol:
            return SETTINGS_MISMATCH
    return SETTINGS_UNVERIFIED if unverified else SETTINGS_MATCH
