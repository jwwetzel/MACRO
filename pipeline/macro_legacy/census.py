"""Legacy-archive census — the pure decision logic (L0, L1, L2).

Every function here is side-effect free: no file is opened and no database
is touched.  ``pipeline/scripts/build_legacy_census.py`` feeds these
functions the ``scan`` table and writes what they return; the unit tests in
``pipeline/tests/test_legacy_census.py`` feed them hand-built cases.

The definitions in this module are NOT free choices.  The ones that decide
the go/no-go (canonical frame, night, series, run, season, calibrated
night, minimum-bearing night, time-convention pass, gates G0–G3) were
written down in ``Legacy_Rigel/notes/GO_NOGO_PREREGISTERED.md`` before any
header was scanned, and each constant below cites the clause it implements.
Changing one after the census has been seen is forbidden by §5 of that
note; a definition that proves unworkable is reported as a deviation, with
the result under both readings.

Sections
--------
1.  Paths and the transfer-manifest reconciliation
2.  Time: DATE-OBS parsing, night labels, path dates, sidereal check
3.  Camera identity and the mechanical timeline
4.  Frame kind, filters, targets
5.  The dedup rule (canonical frames)
6.  Runs and seasons
7.  Calibration availability and flat pairs
8.  Eclipsing systems and minimum-bearing nights
9.  The header-time convention audit
10. The pre-registered gates and the decision rule
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from typing import Iterable, Mapping, Optional, Sequence

from macro_core import manifest as s0   # the shared S0 name/night rules

# ===========================================================================
# 1.  PATHS AND THE TRANSFER-MANIFEST RECONCILIATION
# ===========================================================================

#: Suffixes that mark an fpack-compressed file.
_FZ = ".fz"


def _missing(x) -> bool:
    """True for None and for the float NaN pandas substitutes for it."""
    return x is None or (isinstance(x, float) and math.isnan(x))


def logical_path(path: str) -> str:
    """The archive-relative name of the EXPOSURE a file holds.

    ``2015/day004/frm00400.fts.fz`` and ``2015/day004/frm00400.fts`` are the
    same exposure in two encodings; the transfer manifest lists the second
    form.  Stripping a trailing ``.fz`` puts disk and manifest in one
    namespace.
    """
    return path[:-len(_FZ)] if path.lower().endswith(_FZ) else path


#: Reconciliation statuses.  The first is a census row; every other one is a
#: NAMED EXCLUSION (pre-registration §2, "files-on-disk = rows + named
#: exclusions").
REC_SCANNED = "scanned"                       # .fz on disk, header read
REC_UNREADABLE = "fz_unreadable"              # .fz on disk, header unreadable
REC_UNCOMPRESSED_ONLY = "uncompressed_only"   # no .fz twin: never opened
REC_NOT_ON_DISK = "manifest_not_on_disk"      # listed by the crawl, absent
REC_NOT_FITS = "not_fits"                     # a non-FITS file on disk
#: Order used wherever the statuses are tabulated.
REC_ORDER = (REC_SCANNED, REC_UNREADABLE, REC_UNCOMPRESSED_ONLY,
             REC_NOT_ON_DISK, REC_NOT_FITS)


def reconcile_one(in_manifest: bool, has_fz: bool, has_fits: bool,
                  fz_error: bool, is_fits_name: bool = True) -> str:
    """Status of ONE logical path, given where it was found.

    Parameters are plain facts: is the path in the transfer manifest; is a
    ``.fz`` on disk; is an uncompressed copy on disk; did the ``.fz`` header
    fail to read; does the name look like a FITS file at all.

    The rule, in priority order:

    * a readable ``.fz`` is a census row, whatever else exists (an
      uncompressed twin beside it is the dedupe job's unfinished business,
      not a second frame);
    * an unreadable ``.fz`` is a named exclusion — but it still counts as a
      scanned file (it has a row, with ``error`` set);
    * an uncompressed file with no ``.fz`` twin is a named exclusion: the
      scan opens ``.fz`` only (these are the frames fpack could not
      compress — damaged or truncated downloads);
    * a manifest path with nothing on disk is a named exclusion;
    * anything else on disk (logs, thumbnails) is ``not_fits``.
    """
    if has_fz:
        return REC_UNREADABLE if fz_error else REC_SCANNED
    if has_fits:
        return REC_UNCOMPRESSED_ONLY
    if in_manifest:
        return REC_NOT_ON_DISK
    if not is_fits_name:
        return REC_NOT_FITS
    # Not in the manifest, not on disk: cannot be reached from either list.
    raise ValueError("a path that is neither in the manifest nor on disk "
                     "cannot be reconciled")


def strip_first_component(rel_path: str) -> str:
    """What the first (BAD) manifest did to every path: drop the year.

    ``2015/day004/frm00400.fts`` → ``day004/frm00400.fts``.  The collision
    audit re-applies this to the good manifest to count how many distinct
    exposures the bad crawl would have written onto one another.
    """
    head, sep, tail = rel_path.partition("/")
    return tail if sep else head


# ===========================================================================
# 2.  TIME
# ===========================================================================

#: The archive spans 2015–2022 (plus a 2023 calibration folder).  A header
#: date outside this window is not an observation date; it is a clock that
#: was never set (the 1970-01-01 stamps of the AC4040 commissioning nights).
VALID_YEAR_MIN = 2014
VALID_YEAR_MAX = 2023

#: JD of the Unix epoch, and of J2000.0 (for the sidereal-time check).
UNIX_EPOCH_JD = 2440587.5
J2000_JD = 2451545.0

#: Winer Observatory east longitude in degrees (header LONGITUD −110:36:06).
SITE_LONGITUDE_DEG = -(110.0 + 36.0 / 60.0 + 6.0 / 3600.0)

_DATE_OBS_RE = re.compile(
    r"^(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2}):(\d{2})(\.\d+)?$")


def parse_date_obs(text: Optional[str]) -> Optional[float]:
    """DATE-OBS string → JD (UTC), or None if absent / malformed / unset.

    Accepts ``YYYY-MM-DDThh:mm:ss[.fff]``.  Returns None for a date outside
    :data:`VALID_YEAR_MIN`–:data:`VALID_YEAR_MAX`: such a frame has "no
    usable DATE-OBS" and is a named exclusion of every science count
    (pre-registration §2), not a frame observed in 1970.
    """
    if _missing(text) or not text:
        return None
    m = _DATE_OBS_RE.match(str(text).strip())
    if not m:
        return None
    y, mo, d, h, mi, s = (int(g) for g in m.groups()[:6])
    frac = float(m.group(7)) if m.group(7) else 0.0
    if not VALID_YEAR_MIN <= y <= VALID_YEAR_MAX:
        return None
    try:
        moment = datetime(y, mo, d, h, mi, s, tzinfo=timezone.utc)
    except ValueError:
        return None
    return (moment.timestamp() + frac) / 86400.0 + UNIX_EPOCH_JD


def date_obs_decimals(text: Optional[str]) -> Optional[int]:
    """Number of fractional-second digits DATE-OBS carries (0 = whole seconds).

    The timestamp resolution is part of the timing budget: a whole-second
    stamp is uncertain by up to 1 s before any clock question is asked.
    """
    if _missing(text) or not text:
        return None
    m = _DATE_OBS_RE.match(str(text).strip())
    if not m:
        return None
    return len(m.group(7)) - 1 if m.group(7) else 0


def night_label(jd: Optional[float]) -> Optional[str]:
    """Local-noon-to-noon night label (pre-registration §2, "Night").

    Delegates to the S0 rule so that a legacy night and an RLMT night are
    the same thing: the calendar date of (JD − 19 h).
    """
    return s0.night_label(jd)


_DAY_DIR_RE = re.compile(r"(?:^|/)(\d{4})/day(\d{3})/")
_ARCHIVAR_RE = re.compile(r"(?:^|/)archivar/(\d{4})/[A-Za-z0-9_]{3}(\d{3})")
_OLD_RE = re.compile(r"(?:^|/)old/day(\d{3})/")


def path_date(path: str) -> tuple[Optional[int], Optional[int]]:
    """``(year, day-of-year)`` encoded in an archive path, where present.

    Three layouts exist: ``YYYY/dayDDD/…``; ``archivar/YYYY/pppDDDss…``
    (the day is characters 4–6 of the file name); and ``old/dayDDD/…``
    (day known, year not — returned as ``(None, DDD)``).

    The path date is the ONLY date some frames have (unset camera clock),
    and for all others it is an independent check on DATE-OBS: the
    directory was named by the scheduler, not by the camera PC.
    """
    m = _DAY_DIR_RE.search(path)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = _ARCHIVAR_RE.search(path)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = _OLD_RE.search(path)
    if m:
        return None, int(m.group(1))
    return None, None


def doy_of_jd(jd: float, shift_days: float = 0.0) -> tuple[int, int]:
    """``(year, day-of-year)`` of a JD, optionally shifted (UT → local)."""
    moment = (datetime(1970, 1, 1, tzinfo=timezone.utc)
              + timedelta(days=jd - UNIX_EPOCH_JD + shift_days))
    return moment.year, moment.timetuple().tm_yday


def path_date_offset_days(path: str, jd: Optional[float]) -> Optional[int]:
    """Whole days between the path's day-of-year and DATE-OBS's UT date.

    0 means the directory is named for the UT date of the exposure; any
    other value on a whole population is a convention to document, and a
    scattered non-zero value is a misfiled or mis-stamped frame.  The
    result is wrapped into −182…+182 so a New-Year crossing reads as ±1.
    Returns None when either side is unknown.
    """
    year, doy = path_date(path)
    if doy is None or _missing(jd):
        return None
    y, d = doy_of_jd(jd)
    if year is not None and year != y:
        # Different calendar years: express as a signed day count.
        delta = (date(y, 1, 1) + timedelta(days=d - 1)
                 - (date(year, 1, 1) + timedelta(days=doy - 1))).days
        return delta
    diff = d - doy
    if diff > 182:
        diff -= 365
    elif diff < -182:
        diff += 365
    return diff


def parse_sexagesimal(text: Optional[str]) -> Optional[float]:
    """``'H:M:S'`` / ``'D M S'`` → decimal (same unit as the first field).

    Tolerates the Talon style with leading blanks and signs
    (``' -0:27:47.7'``).  Returns None when the text is not three numeric
    fields.
    """
    if _missing(text):
        return None
    parts = str(text).replace(":", " ").split()
    if len(parts) != 3:
        return None
    try:
        a, b, c = (float(p) for p in parts)
    except ValueError:
        return None
    value = abs(a) + b / 60.0 + c / 3600.0
    return -value if parts[0].lstrip().startswith("-") else value


def gmst_hours(jd_ut: float) -> float:
    """Greenwich mean sidereal time (hours) at a UT Julian date.

    IAU 1982 expression, good to ~0.1 s over these decades — far finer than
    the 1-s header LST it is compared with.  UT1−UTC (< 0.9 s) is ignored
    and is part of the stated tolerance of the sidereal check.
    """
    d = jd_ut - J2000_JD
    t = d / 36525.0
    gmst = (280.46061837 + 360.98564736629 * d
            + 0.000387933 * t * t - t * t * t / 38710000.0)
    return (gmst % 360.0) / 15.0


def lst_residual_seconds(jd_ut: Optional[float], lst_text: Optional[str],
                         longitude_deg: float = SITE_LONGITUDE_DEG
                         ) -> Optional[float]:
    """Header LST minus the LST implied by DATE-OBS, in seconds of time.

    This is the one header-internal check that can tell UTC from local
    time and a start stamp from a grossly different instant: the mount
    computed LST from its own clock and the site longitude, the camera
    software wrote DATE-OBS.  A residual near zero says DATE-OBS is UTC and
    refers to (about) the instant the LST was sampled.  It CANNOT certify
    the absolute clock when both cards derive from the same PC clock —
    that needs a measured event (an eclipse minimum against a literature
    ephemeris), which is photometry and outside a header census.

    The result is wrapped into ±12 h and expressed in SI seconds (sidereal
    seconds ÷ 1.0027379).
    """
    lst = parse_sexagesimal(lst_text)
    if _missing(jd_ut) or lst is None:
        return None
    expected = (gmst_hours(jd_ut) + longitude_deg / 15.0) % 24.0
    diff = (lst - expected + 12.0) % 24.0 - 12.0
    return diff * 3600.0 / 1.00273790935


# ===========================================================================
# 3.  CAMERA IDENTITY AND THE MECHANICAL TIMELINE
# ===========================================================================

#: INSTRUME text → short camera name.  First match wins.  Every pattern is a
#: string that actually occurs in the archive's headers (the census build
#: fails loudly on an INSTRUME value none of these recognises, so a camera
#: can never be silently folded into "unknown").
CAMERA_RULES: tuple[tuple[str, str], ...] = (
    (r"PL\s*16803", "FLI PL16803"),
    (r"Apogee.*F47|F47", "Apogee F47"),
    (r"Aspen.*CG42|CG42", "Andor Aspen CG42"),
    (r"IKON.*936|iKon", "Andor iKon-L 936"),
    (r"STXL.*6303|STXL", "SBIG STXL-6303"),
    (r"6303", "SBIG 6303e"),
    (r"AC4040|Aluma", "SBIG Aluma AC4040"),
    # 'DL Imaging' is the Diffraction Limited driver name MaxIm reports for
    # the Aluma AC4040 from February 2023 (4096 px, 9 µm, StackPro readout
    # modes that only that camera has).
    (r"^DL Imaging$", "SBIG Aluma AC4040"),
    # 'Andor Tech' is a vendor, not a model: 2048 px of 13.5 µm is the
    # CCD42-40, which both Andor cameras in this archive carry.
    (r"^Andor Tech$", "CCD42-40 camera (model not in header)"),
)
_CAMERA_RES = tuple((re.compile(p, re.IGNORECASE), name)
                    for p, name in CAMERA_RULES)

#: Label for frames whose header carries no INSTRUME at all.
CAMERA_UNLABELLED = "(no INSTRUME)"


def camera_name(instrume: Optional[str]) -> Optional[str]:
    """Short camera name for an INSTRUME value.

    Returns :data:`CAMERA_UNLABELLED` for a missing/blank card and None for
    a non-blank value that no rule recognises (the build treats None as an
    error to be resolved by adding a rule, never by guessing).
    """
    if _missing(instrume) or not str(instrume).strip():
        return CAMERA_UNLABELLED
    for rx, name in _CAMERA_RES:
        if rx.search(str(instrume)):
            return name
    return None


#: Native (unbinned) pixel pitch in µm of each camera's sensor.  This is the
#: check on INSTRUME: MaxIm writes the pixel size it reads from the camera
#: driver (XPIXSZ, binned), whereas INSTRUME is a label the scheduler
#: copies from its configuration and which can go stale across a camera
#: swap.  A frame whose pixel pitch contradicts its label was NOT taken
#: with the labelled camera.
CAMERA_PIXEL_UM: dict[str, float] = {
    "FLI PL16803": 9.0,          # KAF-16803, 4096 x 4096
    "Apogee F47": 13.0,          # e2v CCD47-10, 1024 x 1024
    "SBIG 6303e": 9.0,           # KAF-6303E, 3072 x 2048
    "SBIG STXL-6303": 9.0,       # KAF-6303E, 3072 x 2048
    "Andor Aspen CG42": 13.5,    # e2v CCD42-40, 2048 x 2048
    "Andor iKon-L 936": 13.5,    # e2v CCD42-40, 2048 x 2048
    "SBIG Aluma AC4040": 9.0,    # GSENSE4040, 4096 x 4096
    "CCD42-40 camera (model not in header)": 13.5,
}
#: Largest native format (x pixels) of each sensor: a frame wider than this
#: cannot have come from it either.
CAMERA_NATIVE_X: dict[str, int] = {
    "FLI PL16803": 4096, "Apogee F47": 1024, "SBIG 6303e": 3072,
    "SBIG STXL-6303": 3072, "Andor Aspen CG42": 2048,
    "Andor iKon-L 936": 2048, "SBIG Aluma AC4040": 4096,
    "CCD42-40 camera (model not in header)": 2048,
}
#: Name given to frames identified by their sensor signature alone (label
#: missing or contradicted), keyed by (native pixel µm, native x ≤ limit).
#: Two cameras share each of these sensors, so the MODEL cannot be
#: recovered from the header — the name says so instead of guessing.
SENSOR_ONLY_NAMES: tuple[tuple[float, int, str], ...] = (
    (13.0, 1024, "Apogee F47"),            # unique sensor in this archive
    (9.0, 3072, "KAF-6303 camera (model not in header)"),
    (13.5, 2048, "CCD42-40 camera (model not in header)"),
    (9.0, 4096, "4096-px 9-µm camera (model not in header)"),
)
#: Pixel-pitch tolerance (µm) when comparing header and sensor.
PIXEL_TOL_UM = 0.26


def camera_identity(instrume: Optional[str], pixel_um: Optional[float],
                    native_x: Optional[float]) -> tuple[Optional[str], str]:
    """``(camera, basis)`` — the camera a frame was really taken with.

    ``pixel_um`` is the NATIVE pixel pitch (header XPIXSZ ÷ binning) and
    ``native_x`` the frame width × binning; either may be None.

    * ``basis = 'label'``            — INSTRUME recognised and not
      contradicted by the sensor signature (or no signature to test with).
    * ``basis = 'sensor (label stale)'`` — INSTRUME names a camera whose
      sensor cannot produce this pixel pitch / frame width; the frame is
      attributed by sensor signature instead.
    * ``basis = 'sensor (no label)'``  — INSTRUME blank; attributed by
      signature.
    * ``basis = 'unidentified'``     — no label and no usable signature.

    Returns ``(None, 'unrecognised label')`` for a non-blank INSTRUME no
    rule knows — the build stops on those rather than guess.
    """
    label = camera_name(instrume)
    if label is None:
        return None, "unrecognised label"
    has_pix = pixel_um is not None and not math.isnan(pixel_um) and pixel_um > 0
    has_x = native_x is not None and not math.isnan(native_x)
    if label != CAMERA_UNLABELLED:
        contradicted = (
            (has_pix and abs(pixel_um - CAMERA_PIXEL_UM[label]) > PIXEL_TOL_UM)
            or (has_x and native_x > CAMERA_NATIVE_X[label]))
        if not contradicted:
            return label, "label"
        basis = "sensor (label stale)"
    else:
        basis = "sensor (no label)"
    if has_pix:
        for pitch, max_x, name in SENSOR_ONLY_NAMES:
            if abs(pixel_um - pitch) <= PIXEL_TOL_UM and \
                    (not has_x or native_x <= max_x):
                return name, basis
    return CAMERA_UNLABELLED, "unidentified"


def binning_of(xbinning: Optional[float], xfactor: Optional[float]
               ) -> Optional[int]:
    """Binning factor from XBINNING (MaxIm) or XFACTOR (Talon/scheduler)."""
    for v in (xbinning, xfactor):
        if v is not None and not (isinstance(v, float) and math.isnan(v)) \
                and v >= 1:
            return int(round(v))
    return None


#: A rotation step larger than this between consecutive nightly medians
#: starts a new mechanical epoch (telescope engineer F1: "rotation step
#: > 0.3°").
ROTATION_STEP_DEG = 0.3
#: A night needs at least this many plate-solved frames for its median
#: rotation to count (single solves scatter by more than the step).
ROTATION_MIN_FRAMES = 5


def fold_rotation(crota_deg: float) -> float:
    """Fold a position angle into [−90, +90): a pier flip is not a re-mount.

    A German-equatorial flip rotates the field by exactly 180°; that is an
    operating state, recorded separately (FLIPSTAT), not a mechanical
    change.  Folding modulo 180° removes it so that only a genuine camera
    re-seating (a few tenths of a degree to a few degrees) shows as a step.
    """
    return (crota_deg + 90.0) % 180.0 - 90.0


def _rot_delta(a: float, b: float) -> float:
    """Signed difference a − b on the 180°-periodic circle."""
    return (a - b + 90.0) % 180.0 - 90.0


def rotation_epochs(nightly: Sequence[tuple[str, float]],
                    step_deg: float = ROTATION_STEP_DEG
                    ) -> list[tuple[str, str, float, int, int]]:
    """Segment nightly median rotations into epochs at steps > ``step_deg``.

    ``nightly`` is ``[(night, folded median rotation), …]`` sorted by night.

    Pass 1 — a new segment starts when a night's rotation differs from the
    RUNNING MEDIAN of the current segment by more than the step (comparison
    with the segment, not the previous night, so slow drift within
    tolerance does not chain into a false break).

    Pass 2 — a ONE-night segment whose two neighbours agree with each other
    within the step is an *excursion* (one night of bad plate solutions, or
    a camera put on and taken off again), not a re-mount: the three are
    merged and the night is counted in ``n_excursions``.  A one-night
    segment between two DIFFERENT neighbours is kept: that is where the
    change happened.

    Returns ``[(first night, last night, median rotation, n nights,
    n excursion nights), …]``.
    """
    segs: list[dict] = []
    for night, rot in nightly:
        if segs and abs(_rot_delta(rot, _median(segs[-1]["rots"]))) <= step_deg:
            segs[-1]["last"] = night
            segs[-1]["rots"].append(rot)
        else:
            segs.append(dict(first=night, last=night, rots=[rot], exc=0))
    merged = True
    while merged:
        merged = False
        for i in range(1, len(segs) - 1):
            a, mid, b = segs[i - 1], segs[i], segs[i + 1]
            if len(mid["rots"]) == 1 and abs(_rot_delta(
                    _median(a["rots"]), _median(b["rots"]))) <= step_deg:
                a["rots"] += b["rots"]
                a["last"] = b["last"]
                a["exc"] += mid["exc"] + b["exc"] + 1
                del segs[i:i + 2]
                merged = True
                break
    return [(g["first"], g["last"], _median(g["rots"]), len(g["rots"]),
             g["exc"]) for g in segs]


def _median(values: Sequence[float]) -> float:
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


# ===========================================================================
# 4.  FRAME KIND, FILTERS, TARGETS
# ===========================================================================

KIND_LIGHT = "light"
KIND_BIAS = "bias"
KIND_DARK = "dark"
KIND_FLAT = "flat"
KIND_FOCUS = "focus"
KIND_TEST = "test"

#: File-name prefixes (first three characters) the scheduler reserves for
#: its own housekeeping frames rather than an observer's programme.
FOCUS_PREFIXES = frozenset({"foc"})

#: OBJECT values that name an engineering activity, not a sky target.  Keys
#: are compared after :func:`macro_core.manifest.normalize_target`.
_TEST_OBJECT_RE = re.compile(
    r"^(test\d*|tests|focus\d*|focusing|pointing\w*|point\d+|flat\d*|flats|"
    r"skyflat\d*|skyflats|twilight\d*|dark\d*|darks|bias\d*|zenith|park|"
    r"sync\d*|tpoint\w*|mapping\d*|collimation\d*|unknown|none|object)$")
# NOTE the anchors and the digit-only tails: 'Dark Doodad' or 'Flaming Star'
# normalise to 'darkdoodad' / 'flamingstar' and must stay science targets.


def frame_kind(imagetyp: Optional[str], basename: str,
               object_name: Optional[str], exptime: Optional[float]) -> str:
    """Classify one frame as light / bias / dark / flat / focus / test.

    Evidence, strongest first:

    1. the file-name prefix ``foc`` — the scheduler's autofocus frames (they
       carry ``IMAGETYP = Light Frame`` and a star's name, so the header
       alone would count them as science);
    2. IMAGETYP, where the software writes one (MaxIm: ``Light Frame``,
       ``Bias Frame``, ``Dark Frame``, ``Flat Field``);
    3. where IMAGETYP is absent (Talon, 2015): OBJECT text naming a
       calibration kind, then a zero exposure as bias;
    4. an OBJECT that names an engineering activity (``test``,
       ``Pointing4W`` …) makes a light frame a ``test`` frame.

    The must-NOT case this ordering protects: a science frame is never
    promoted to a calibration kind by its FILTER or target name alone.
    """
    stem = basename.lower()
    if stem[:3] in FOCUS_PREFIXES:
        return KIND_FOCUS
    it = "" if _missing(imagetyp) else str(imagetyp).strip().lower()
    obj_key = target_key(object_name)[0] or ""
    if it:
        if "bias" in it:
            return KIND_BIAS
        if "dark" in it:
            return KIND_DARK
        if "flat" in it:
            return KIND_FLAT
    else:
        if obj_key.startswith("bias"):
            return KIND_BIAS
        if obj_key.startswith("dark"):
            return KIND_DARK
        if obj_key.startswith(("flat", "skyflat", "twilight", "domeflat")):
            return KIND_FLAT
        if not _missing(exptime) and exptime == 0 and not obj_key:
            return KIND_BIAS
    if _TEST_OBJECT_RE.match(obj_key):
        return KIND_TEST
    return KIND_LIGHT


@dataclass(frozen=True)
class FilterInfo:
    """A normalized filter: short band name, photometric system, tieability."""
    band: str        #: series key, e.g. 'V', 'g', 'Red', 'Ha', 'clear'
    system: str      #: 'Johnson-Cousins', 'Sloan', 'RGB', 'narrow', 'none', …
    tieable: bool    #: can be tied to an RLMT-era band (pre-reg. G2 list)


#: Talon (2015, the Rigel) wrote ONE LETTER for a Johnson–Cousins + Hα
#: wheel.  Only under Talon does a bare letter name a photometric band.
TALON_LETTERS: dict[str, FilterInfo] = {
    "B": FilterInfo("B", "Johnson-Cousins", True),
    "V": FilterInfo("V", "Johnson-Cousins", True),
    "R": FilterInfo("R", "Johnson-Cousins", True),
    "I": FilterInfo("I", "Johnson-Cousins", True),
    "N": FilterInfo("clear", "none", False),
    "H": FilterInfo("Ha", "narrow", False),
}

#: MaxIm wrote ``'<slot> - <name>'`` in many spellings.  The NAME part
#: (lower-cased, punctuation folded) is matched against these patterns in
#: order; first match wins.  THE DISTINCTION THAT MATTERS: a named colour
#: (``Red``, ``Blue``, ``Green``, ``Visual``) does not say which glass was in
#: the slot, so it is its own band and is NOT tieable under gate G2, whose
#: list is B, V, R, I, g, r, i "and their Johnson/Sloan counterparts".  Only
#: a name that says Sloan or Johnson–Cousins ("J-C", "JC") is tieable.
FILTER_NAME_RULES: tuple[tuple[str, FilterInfo], ...] = (
    (r"grism|spectrom|beam ?s?plitter", FilterInfo("grism", "dispersive", False)),
    (r"sloan ?_?'?g'?$", FilterInfo("g", "Sloan", True)),
    (r"sloan ?_?'?r'?$", FilterInfo("r", "Sloan", True)),
    (r"sloan ?_?'?i'?$", FilterInfo("i", "Sloan", True)),
    (r"sloan ?_?'?z'?$", FilterInfo("z", "Sloan", False)),
    (r"^(j-?c|johnson)[ _]*b(lue)?$", FilterInfo("B", "Johnson-Cousins", True)),
    (r"^(j-?c|johnson)[ _]*v(isual)?$", FilterInfo("V", "Johnson-Cousins", True)),
    (r"^(j-?c|johnson)[ _]*r(ed)?$", FilterInfo("R", "Johnson-Cousins", True)),
    (r"^(j-?c|johnson)[ _]*i$", FilterInfo("I", "Johnson-Cousins", True)),
    (r"visual$", FilterInfo("Visual", "named colour (glass not stated)", False)),
    (r"red$", FilterInfo("Red", "named colour (glass not stated)", False)),
    (r"^green$", FilterInfo("Green", "named colour (glass not stated)", False)),
    (r"^blue$", FilterInfo("Blue", "named colour (glass not stated)", False)),
    (r"^luminances?$", FilterInfo("Lum", "broadband", False)),
    (r"^h(ydrogen)? ?alpha$", FilterInfo("Ha", "narrow", False)),
    (r"^(oiii|oxygen ?(iii)?)$", FilterInfo("OIII", "narrow", False)),
    (r"^sulphur ?ii$|^sii$", FilterInfo("SII", "narrow", False)),
    (r"^longpass ?(ir)?$", FilterInfo("W", "longpass", False)),
    (r"^none$", FilterInfo("clear", "none", False)),
)
_FILTER_NAME_RES = tuple((re.compile(p), info) for p, info in FILTER_NAME_RULES)
_SLOT_NAME_RE = re.compile(r"^(\w)\s*-\s*(.+)$")


def filter_info(raw: Optional[str], talon: bool = False) -> FilterInfo:
    """Normalize a header FILTER value.

    ``talon`` says the frame was written by Talon (the 2015 Rigel), the
    only software whose single letters name photometric bands.

    * ``'<slot> - <name>'`` (MaxIm): the name is matched against
      :data:`FILTER_NAME_RULES`.
    * a bare token under Talon: :data:`TALON_LETTERS`.
    * a bare token under MaxIm (``'R'``, ``'6'``): a wheel SLOT whose glass
      the header does not name — band ``'slot R'``, never tieable.

    Unknown strings are NOT guessed into a band: they become a band named
    by the cleaned raw text with system ``'unrecognised'``, so they form
    their own series and are visible in the filter census.  A missing card
    is band ``'(none)'``.
    """
    if _missing(raw) or not str(raw).strip():
        return FilterInfo("(none)", "missing", False)
    text = " ".join(str(raw).split())
    m = _SLOT_NAME_RE.match(text)
    if m:
        name = " ".join(m.group(2).lower().replace("_", " ").split())
        for rx, info in _FILTER_NAME_RES:
            if rx.search(name):
                return info
        return FilterInfo(text, "unrecognised", False)
    if talon and text.upper() in TALON_LETTERS:
        return TALON_LETTERS[text.upper()]
    if len(text) == 1:
        return FilterInfo(f"slot {text}", "slot letter (glass not stated)",
                          False)
    return FilterInfo(text, "unrecognised", False)


def target_key(object_name: Optional[str]) -> tuple[Optional[str], Optional[str]]:
    """``(alias-group key, cleaned display name)`` for an OBJECT value.

    Exactly the S0 rule chain (``macro_core.manifest.normalize_target``):
    case, spaces, underscores and letter-hyphens folded, genitives and the
    explicit synonym table applied.  Using the same function is what makes
    "T CrB in the legacy archive" and "T CrB in the RLMT manifest" the same
    key for the overlap query.
    """
    n = s0.normalize_target(None if _missing(object_name) else object_name)
    return n.key, n.cleaned


# ===========================================================================
# 5.  THE DEDUP RULE (CANONICAL FRAMES)
# ===========================================================================

def dup_key(camera: Optional[str], date_obs: Optional[str],
            exptime: Optional[float], filter_band: str,
            naxis1: Optional[int], naxis2: Optional[int],
            row_id) -> tuple:
    """Header identity key of one exposure (pre-registration §2).

    ``(camera, DATE-OBS as recorded, exposure, filter, geometry)``.  Two
    files with the same key are copies of one exposure.  A frame with no
    usable DATE-OBS cannot be identified this way and gets a key unique to
    itself — it is never merged with anything, and it is excluded from
    science counts by the "no usable DATE-OBS" clause, not by dedup.
    """
    if parse_date_obs(date_obs) is None:
        return ("__nodate__", row_id)
    exp = None if _missing(exptime) else round(float(exptime), 3)
    return (camera, date_obs.strip(), exp, filter_band, naxis1, naxis2)


def canonical_rank(calstat: Optional[str], path: str) -> tuple:
    """Sort key choosing the canonical member of a duplicate group.

    Pre-registered order: "a raw file over a reduced/calibrated twin, then
    the lexicographically smallest path".  In this archive a calibrated
    file is one whose header carries a non-blank ``CALSTAT`` (the scheduler
    stamps ``B``/``D``/``F`` when it applies bias/dark/flat in place).
    """
    stamped = not _missing(calstat) and bool(str(calstat).strip())
    return (1 if stamped else 0, path)


# ===========================================================================
# 6.  RUNS AND SEASONS
# ===========================================================================

#: Pre-registration §2, "Run": no gap longer than 30 min.
RUN_GAP_DAYS = 30.0 / 1440.0
#: Float slack on that comparison (half a second): stamps are whole seconds,
#: so a gap of exactly 30:00 occurs and must not split on a rounding bit.
_RUN_GAP_EPS_DAYS = 0.5 / 86400.0
#: Pre-registration §2, "Season": nights clustered by gaps > 120 d.
SEASON_GAP_DAYS = 120


def split_runs(jds: Sequence[float], gap_days: float = RUN_GAP_DAYS
               ) -> list[tuple[float, float, int]]:
    """Split one night's exposure-start JDs of one series into runs.

    Returns ``[(first start, last start, n frames), …]``.  A run's length is
    last start − first start (pre-registered), so a single frame is a run
    of length zero.  Input need not be sorted.
    """
    runs: list[list] = []
    for jd in sorted(jds):
        if runs and jd - runs[-1][1] <= gap_days + _RUN_GAP_EPS_DAYS:
            runs[-1][1] = jd
            runs[-1][2] += 1
        else:
            runs.append([jd, jd, 1])
    return [(a, b, n) for a, b, n in runs]


def night_ordinal(night: str) -> int:
    """Proleptic-Gregorian ordinal of a ``YYYY-MM-DD`` night label."""
    return date.fromisoformat(night).toordinal()


def assign_seasons(nights: Iterable[str],
                   gap_days: int = SEASON_GAP_DAYS) -> dict[str, int]:
    """Map each night to a season index (1, 2, …) for ONE target.

    Seasons are clusters of that target's nights separated by gaps longer
    than ``gap_days`` — a definition that needs no coordinates and is the
    same for circumpolar and ecliptic targets.
    """
    out: dict[str, int] = {}
    season, previous = 0, None
    for night in sorted(set(nights)):
        o = night_ordinal(night)
        if previous is None or o - previous > gap_days:
            season += 1
        out[night] = season
        previous = o
    return out


# ===========================================================================
# 7.  CALIBRATION AVAILABILITY AND FLAT PAIRS
# ===========================================================================

#: Pre-registration §2, "Calibrated night": within ±30 nights.
CALIB_WINDOW_NIGHTS = 30


def within_window(night_ord: int, sorted_ords: Sequence[int],
                  window: int = CALIB_WINDOW_NIGHTS) -> bool:
    """Is any element of ``sorted_ords`` within ±``window`` of ``night_ord``?"""
    import bisect
    i = bisect.bisect_left(sorted_ords, night_ord - window)
    return i < len(sorted_ords) and sorted_ords[i] <= night_ord + window


def is_calibrated_night(night_ord: int, flat_ords: Sequence[int],
                        zero_ords: Sequence[int],
                        window: int = CALIB_WINDOW_NIGHTS) -> bool:
    """The pre-registered availability test for one series night.

    ``flat_ords`` — sorted ordinals of nights with a flat in the SAME
    camera, binning and filter.  ``zero_ords`` — sorted ordinals of nights
    with a bias OR dark on the same camera and binning.  Both must have a
    member within the window.  This says the frames EXIST; it does not say
    they are good, matched in temperature, or unsaturated.
    """
    return (within_window(night_ord, flat_ords, window)
            and within_window(night_ord, zero_ords, window))


#: Two flats count as a pair when their exposure times agree to this
#: relative tolerance (same illumination level is then plausible; the
#: photon-transfer measurement itself checks the levels on pixels).
FLAT_PAIR_REL_TOL = 0.005


def count_flat_pairs(exptimes: Sequence[float],
                     rel_tol: float = FLAT_PAIR_REL_TOL) -> tuple[int, int]:
    """``(n pairs, n distinct exposure levels with a pair)`` in one group.

    The group is one camera × binning × filter × night.  Exposures are
    clustered by relative tolerance; a cluster of k flats yields ⌊k/2⌋
    disjoint pairs.  The number of distinct levels is what a photon-transfer
    curve needs (one level gives one point; a gain needs a slope).
    """
    pairs = levels = 0
    cluster: list[float] = []
    for e in sorted(x for x in exptimes if x is not None and x > 0):
        if cluster and e - cluster[0] > rel_tol * cluster[0]:
            pairs += len(cluster) // 2
            levels += len(cluster) >= 2
            cluster = []
        cluster.append(e)
    pairs += len(cluster) // 2
    levels += len(cluster) >= 2
    return pairs, levels


# ===========================================================================
# 8.  ECLIPSING SYSTEMS AND MINIMUM-BEARING NIGHTS
# ===========================================================================

#: VSX type tokens that are NOT stellar eclipses although they start with
#: 'E': planetary transits, ellipsoidal variables, EX Lupi-type eruptives.
_NOT_ECLIPSING = frozenset({"EP", "ELL", "EXOR"})


def is_eclipsing_type(vsx_type: Optional[str]) -> bool:
    """Does a VSX variability type describe an eclipsing/contact binary?

    VSX types combine tokens with ``+``, ``/`` and ``|`` (``EW/KW``,
    ``EA+BY``, ``UGSU+E``).  A system is eclipsing when any token is ``E``,
    ``EA``, ``EB``, ``EW`` or ``EC`` (the GCVS eclipse classes), optionally
    with an uncertainty colon.  Planet transits (``EP``) and ellipsoidal
    variables (``ELL``) are excluded: the gate is about stellar period
    change read from eclipse minima.

    This token rule was fixed in code before gate G1 was evaluated; the
    report gives G1 under the wider reading (``EP`` included) as well.
    """
    if not vsx_type:
        return False
    for token in re.split(r"[+/|]", vsx_type.upper()):
        token = token.strip().rstrip(":")
        if token in _NOT_ECLIPSING:
            continue
        if token in ("E", "EA", "EB", "EW", "EC"):
            return True
    return False


#: Pre-registration §2: without a period, the guaranteed clause is run ≥ 3 h.
NO_PERIOD_RUN_DAYS = 3.0 / 24.0
#: Pre-registration §2: a predicted minimum needs ≥ 30 min of data each side.
MINIMUM_MARGIN_DAYS = 30.0 / 1440.0


def run_guarantees_minimum(run_days: float, period: Optional[float]) -> bool:
    """Clause (a): a run of length ≥ P/2 contains a minimum whatever the phase.

    Minima of an eclipsing binary recur every P/2 (primary, secondary), so
    any window at least that long holds one.  With no catalogue period the
    pre-registered fallback is a run of at least three hours.
    """
    need = period / 2.0 if period and period > 0 else NO_PERIOD_RUN_DAYS
    return run_days >= need


def predicted_minima_in_run(start: float, end: float, epoch: Optional[float],
                            period: Optional[float],
                            margin_days: float = MINIMUM_MARGIN_DAYS
                            ) -> list[float]:
    """Clause (b): predicted minima inside a run with the required margin.

    ``start``/``end`` and ``epoch`` must be on one time scale (the build
    converts run limits to HJD before calling, since catalogue epochs are
    heliocentric).  Minima are predicted at ``epoch + k·P/2`` for integer k
    — primary and secondary alike.  Returns the predicted times that lie in
    ``[start + margin, end − margin]``.

    The prediction inherits the catalogue ephemeris error (a period wrong by
    1e-6 d is ~30 min after 2×10⁴ cycles), which is why this clause is
    published separately from clause (a) and never alone decides a claim.
    """
    if not epoch or not period or period <= 0:
        return []
    lo, hi = start + margin_days, end - margin_days
    if hi < lo:
        return []
    half = period / 2.0
    k0 = math.ceil((lo - epoch) / half)
    out = []
    t = epoch + k0 * half
    while t <= hi:
        out.append(t)
        t += half
    return out


# ===========================================================================
# 9.  THE HEADER-TIME CONVENTION AUDIT
# ===========================================================================

#: Pre-registration §2, "Time-convention pass": independent card within 2 s
#: on ≥ 99% of the frames that carry both.
TIME_AGREE_SECONDS = 2.0
TIME_AGREE_FRACTION = 0.99
#: Slack allowed when testing whether consecutive exposures overlap: one
#: whole-second timestamp quantum on each stamp.
OVERLAP_SLACK_SECONDS = 2.0


def card_agreement(diffs_seconds: Sequence[float],
                   tol: float = TIME_AGREE_SECONDS) -> tuple[int, int, float]:
    """``(n, n within tol, fraction)`` for a list of card differences."""
    n = len(diffs_seconds)
    good = sum(1 for d in diffs_seconds if abs(d) < tol)
    return n, good, (good / n if n else float("nan"))


def stamp_hypothesis_violations(frames: Sequence[tuple[float, float]],
                                slack_s: float = OVERLAP_SLACK_SECONDS
                                ) -> dict[str, int]:
    """Test what DATE-OBS marks, from exposure overlaps alone.

    ``frames`` is ``[(jd, exptime_s), …]`` for consecutive exposures of ONE
    camera on one night, in time order.  A camera cannot expose two frames
    at once, so:

    * if the stamp is the exposure **start**, then
      ``stamp[i+1] − stamp[i] ≥ exptime[i]``;
    * if the stamp is the exposure **end**, then
      ``stamp[i+1] − stamp[i] ≥ exptime[i+1]``;
    * if it is the **middle**, the gap is ≥ the mean of the two.

    With equal exposures the three are indistinguishable, so only pairs
    with ``exptime[i] ≠ exptime[i+1]`` are counted as informative.  The
    hypothesis with zero (or by far the fewest) violations is the
    convention the software used.  Returned keys: ``pairs`` (all consecutive
    pairs), ``informative``, ``viol_start``, ``viol_mid``, ``viol_end``,
    and ``viol_start_all`` (start-hypothesis violations over ALL pairs —
    the "overlapping exposures" count used by the monotonicity clause).
    """
    out = dict(pairs=0, informative=0, viol_start=0, viol_mid=0, viol_end=0,
               viol_start_all=0)
    for (jd0, e0), (jd1, e1) in zip(frames, frames[1:]):
        if e0 is None or e1 is None:
            continue
        gap = (jd1 - jd0) * 86400.0
        out["pairs"] += 1
        if gap < e0 - slack_s:
            out["viol_start_all"] += 1
        if abs(e0 - e1) <= slack_s:
            continue
        out["informative"] += 1
        out["viol_start"] += gap < e0 - slack_s
        out["viol_end"] += gap < e1 - slack_s
        out["viol_mid"] += gap < 0.5 * (e0 + e1) - slack_s
    return out


def time_convention_pass(n_both: int, n_agree: int, n_duplicate_stamps: int,
                         n_overlaps: int,
                         frac_required: float = TIME_AGREE_FRACTION) -> bool:
    """The pre-registered pass for one camera/software epoch.

    Requires (i) an independent time card present on at least one frame and
    agreeing with DATE-OBS within 2 s on ≥ 99% of the frames that carry
    both, (ii) no duplicate timestamps within a run and (iii) no
    non-monotonic (overlapping) timestamps within a run.  An epoch with no
    second card at all cannot pass — "cannot be checked" is not "agrees".
    """
    if n_both <= 0:
        return False
    return (n_agree / n_both >= frac_required
            and n_duplicate_stamps == 0 and n_overlaps == 0)


# ===========================================================================
# 10. THE PRE-REGISTERED GATES AND THE DECISION RULE
# ===========================================================================

#: G1 — pre-registration §4.
G1_MIN_SYSTEMS = 3
G1_MIN_SEASONS = 3
#: G2 — pre-registration §4.
G2_MIN_NIGHTS = 10
G2_MIN_SEASONS = 2
G2_TIMING_RUN_DAYS = 1.0 / 24.0
#: G3 — pre-registration §4.
G3_MIN_CALIBRATED_NIGHTS = 30

OUTCOME_NOT_DECIDABLE = "NOT DECIDABLE — fix the census"
OUTCOME_RELEASE_ONLY = "NO-GO — archive data release only"
OUTCOME_TRANSFER = "TRANSFER — overlap frames handed to the owning project"
OUTCOME_GO_CANDIDATE = "GO-CANDIDATE"


def g1_passes(seasons_per_system: Sequence[int]) -> bool:
    """G1: ≥ 3 eclipsing systems each with minimum-bearing nights in ≥ 3 seasons."""
    return sum(1 for s in seasons_per_system
               if s >= G1_MIN_SEASONS) >= G1_MIN_SYSTEMS


def g2_series_passes(n_nights: int, n_seasons: int, tieable: bool,
                     is_timing_target: bool, longest_run_days: float) -> bool:
    """G2 for ONE series of an RLMT-era target.

    ≥ 10 nights over ≥ 2 seasons in a tieable filter; a timing target (the
    polars, YZ Cnc) additionally needs one run ≥ 1 h.
    """
    if not tieable or n_nights < G2_MIN_NIGHTS or n_seasons < G2_MIN_SEASONS:
        return False
    if is_timing_target and longest_run_days < G2_TIMING_RUN_DAYS:
        return False
    return True


def g3_series_passes(n_calibrated_nights: int) -> bool:
    """G3 for ONE series: ≥ 30 calibrated nights."""
    return n_calibrated_nights >= G3_MIN_CALIBRATED_NIGHTS


def decide(g0: bool, g1: bool, g2: bool, g3: bool,
           g3_question_named: bool = False) -> list[str]:
    """Apply the pre-registered decision rule; return every outcome that applies.

    1. G0 fails → NOT DECIDABLE (alone, whatever else holds).
    2. G0 holds, none of G1–G3 → release only.
    3. G2 → TRANSFER (not a sixth project).
    4. G1, or G3 with a nameable question and reader → GO-CANDIDATE.
       G3 without a nameable question falls back to rule 2 — unless rule 3
       or G1 already applies.

    Rules 3 and 4 can both apply; both are then returned.
    """
    if not g0:
        return [OUTCOME_NOT_DECIDABLE]
    outcomes = []
    if g2:
        outcomes.append(OUTCOME_TRANSFER)
    if g1 or (g3 and g3_question_named):
        outcomes.append(OUTCOME_GO_CANDIDATE)
    return outcomes or [OUTCOME_RELEASE_ONLY]


# ===========================================================================
# 11. THE LEGACY FILENAME CONVENTION (RIG-L0-filename-parser)
# ===========================================================================

#: The scheduler's name for a frame: three lower-case letters, the
#: three-digit day of year, two hexadecimal digits, then the FITS suffix —
#: ``gbm10028.fts`` = request code ``gbm``, day 100, sequence 0x28.  Seen on
#: every camera and both acquisition systems; ``foc`` is the autofocus
#: request.
_SCHED_NAME_RE = re.compile(
    r"^([a-z]{3})(\d{3})([0-9a-f]{2})\.(?:fts|fit|fits)(?:\.fz)?$")
CONV_SCHEDULER = "scheduler (ppp DDD ss)"
CONV_OTHER = "other"


@dataclass(frozen=True)
class LegacyName:
    """What a legacy file name says.  Fields are None when absent."""
    convention: str
    request: Optional[str]      #: three-letter request/user code
    doy: Optional[int]          #: day of year (UT date of the night's start)
    seq: Optional[int]          #: sequence number within the day (hex)


def parse_legacy_filename(basename: str) -> LegacyName:
    """Decode one legacy file name.

    The convention carries a REQUEST CODE, a DAY OF YEAR and a SEQUENCE
    NUMBER — and nothing else.  Target, filter and exposure are not in the
    name (unlike the RLMT-era ``user_target_filter_exp_date`` names), so
    they are read from the header for every frame; the parser's job is to
    recover what the name does hold and let the build check it against the
    header (day of year against DATE-OBS, request code against OBSERVER).
    """
    m = _SCHED_NAME_RE.match(basename.lower())
    if not m:
        return LegacyName(CONV_OTHER, None, None, None)
    doy = int(m.group(2))
    if not 1 <= doy <= 366:
        return LegacyName(CONV_OTHER, None, None, None)
    return LegacyName(CONV_SCHEDULER, m.group(1), doy, int(m.group(3), 16))


# ===========================================================================
# 12. MULTI-ARCHIVE KEYING (RIG-L0-multi-archive)
# ===========================================================================

#: Value of the ``archive_root`` column for every legacy row (the RLMT
#: manifest's rows live under ``rlmt-archive``).
ARCHIVE_ROOT = "legacy-archive"


def legacy_era_key(camera: Optional[str], focallen_mm: Optional[float],
                   telescop: Optional[str]) -> tuple[str, str]:
    """``(camera, optics)`` — the legacy era key.

    Keyed on CAMERA and FOCAL LENGTH, not on a telescope name: the 0.5 m
    was called 'Gemini', 'Iowa Robotic Telescope' and 'Robert L. Mutel
    Telescope' over the years without any optical change, and two names for
    one configuration must not split it.  Talon headers carry no FOCALLEN;
    for them the TELESCOP text stands in, labelled as such, so the Rigel
    configuration can never share a key with the 0.5 m.
    """
    cam = camera or CAMERA_UNLABELLED
    if not _missing(focallen_mm) and focallen_mm and focallen_mm > 0:
        return cam, f"f = {int(round(float(focallen_mm)))} mm"
    tel = "" if _missing(telescop) else str(telescop).strip()
    return cam, f"telescope '{tel or 'unknown'}' (no FOCALLEN)"


# ===========================================================================
# 13. TRUNCATED-AT-SOURCE FRAMES (RIG-L0-dedup-reconcile)
# ===========================================================================

TRUNC_TWIN_ADOPTED = "truncated_twin_adopted"
TRUNC_TWIN_MISMATCH = "truncated_twin_mismatch"
TRUNC_LOST = "truncated_lost"
TRUNC_UNREADABLE = "truncated_header_unreadable"


def truncated_disposition(trunc_date_obs: Optional[str],
                          twin_date_obs: Sequence[Optional[str]]) -> str:
    """What happens to one frame that is truncated at the source.

    A twin is an intact ``.fz`` of the SAME BASENAME elsewhere in the
    archive (another year or tree).  The twin stands in for the truncated
    frame only if it carries the identical ``DATE-OBS`` — same name alone is
    not identity (the year-stripped crawl showed that names repeat across
    years).  No twin with the same stamp ⇒ the exposure is lost; a twin
    whose stamp differs is a different exposure and the truncated one is
    still lost (recorded separately so the mismatch is visible).
    """
    if _missing(trunc_date_obs) or not str(trunc_date_obs).strip():
        return TRUNC_UNREADABLE
    stamp = str(trunc_date_obs).strip()
    twins = [str(t).strip() for t in twin_date_obs if not _missing(t) and t]
    if stamp in twins:
        return TRUNC_TWIN_ADOPTED
    if twins:
        return TRUNC_TWIN_MISMATCH
    return TRUNC_LOST
