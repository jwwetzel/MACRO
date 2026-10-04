"""Pure S0 manifest logic for the MACRO/RLMT archive.

Everything in this module is a *pure function* or a *data table*: no I/O, no
globals mutated, no hidden state.  The build script
(``pipeline/scripts/build_s0_manifest.py``) wires these functions to the
observation catalog; the unit tests (``pipeline/tests/test_manifest.py``)
exercise every function on hand-built cases, including the cases that MUST
NOT work (e.g. two distinct dwarf fields merging).

The binding conventions implemented here come from ROADMAP.md section 1:

1.  Dedup is GLOBAL on (exposure, jd) across and within trees.  Until
    2026-10-03 "exposure" meant the basename; it now means the basename's
    EXPOSURE ROOT — the name with its packaging (``.fts`` vs ``.fts.fz``)
    and its processing suffixes (``_calibrated``, ``_cal``, ``_wcs``)
    removed — because a reduced or re-packaged copy of a frame is the same
    exposure, not a second one (finding F-1 / TE.F6; see
    :func:`exposure_name`).
2.  ``rawimage`` is the canonical tree, with a documented, data-driven
    exceptions mechanism (see ``TREE_PRIORITY_EXCEPTIONS``).  Within a
    duplicate group the LEAST-PROCESSED copy outranks tree priority: a raw
    frame always beats its own calibrated derivative.
3.  Era/camera assignment keys on (READOUTM, NAXIS geometry, XBINNING,
    EGAIN) — never on filter name or date.  An era is a HEADER history;
    the hardware history (camera, orientation, wheel map) is the
    ``mech_epoch`` table S0b builds beneath it (finding F-3 / TE.F1).
4.  Night label = local-noon-to-noon calendar date of (JD − 0.7917);
    Winer Observatory sits at UTC−7, so local noon is 19:00 UT.
5.  Header JD is stored as-is (UTC exposure start).  BJD_TDB is stage S3's
    job; nothing here fabricates mid-exposure times.
6.  Alias merges are allowed only via a pure normalization rule, or via the
    explicit synonym table gated by a 0.2-degree coordinate cone.  Fuzzy
    string similarity is never used.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from statistics import median
from typing import Iterable, Optional, Sequence

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these,
# so changing a value here changes both the pipeline and its documentation).
# --------------------------------------------------------------------------

#: Shift subtracted from header JD before taking the calendar date, so that
#: the date rolls over at local noon (Winer is UTC−7 → noon = 19:00 UT, and
#: 19:00 UT is JD fraction 0.29167; 0.29167 + 0.5 = 0.7917 places the date
#: boundary exactly there).  Adopted portfolio-wide (Dwarf convention).
NIGHT_SHIFT_DAYS = 0.7917

#: JD of the Unix epoch 1970-01-01T00:00:00 UTC; used to turn a JD into a
#: Python datetime without any astronomy dependency (keeps this module pure).
UNIX_EPOCH_JD = 2440587.5

#: Radius of the coordinate-cone gate: two distinct raw names may be merged
#: by the synonym table only if their median plate-solved coordinates agree
#: within this angle.  ROADMAP S0 spec value.
CONE_RADIUS_DEG = 0.2

#: A frame whose pointing offset from its target's reference coordinates
#: exceeds this angle is flagged as a pointing outlier (T CrB ground truth:
#: 21/247 grism frames beyond 1 degree).
POINTING_OUTLIER_DEG = 1.0

#: Minimum number of plate-solved canonical frames a target needs before we
#: trust their median coordinates as the target's reference position.
MIN_SOLVED_FOR_REFERENCE = 3

#: Header airmass outside [1, 10] (or the -999 sentinel, which is < 1) is
#: physically impossible and gets flagged — never recomputed here (S3's job).
AIRMASS_MIN = 1.0
AIRMASS_MAX = 10.0

#: Decimal places kept when rounding EGAIN for the era key.  Header EGAIN is
#: a driver-reported float (e.g. 1.05697000026703); three decimals separate
#: every genuinely distinct gain in the catalog while absorbing float fuzz.
EGAIN_DECIMALS = 3


# --------------------------------------------------------------------------
# Tree policy: which copy of a duplicated frame is canonical.
# --------------------------------------------------------------------------

#: Default tree preference, highest priority first.  ``rawimage`` is the
#: canonical pixel source (ROADMAP convention 1).  ``reduced`` is *last*
#: because it is unaudited and demonstrably wrong at least once (BeStar
#: Phecda pair).  Trees absent from this list rank after all listed trees.
DEFAULT_TREE_PRIORITY: tuple[str, ...] = (
    "rawimage",     # canonical pixels, by-night directories
    "macalester",   # school mirror; superset for some targets (see below)
    "external",     # consortium mirror; held exactly one copy of each SN frame
    "iKon",         # Andor iKon camera nights (CV VV Pup 2024-11/12)
    "grism",        # small dedicated grism tree
    "hagrism",
    "coe", "iowa", "knox", "augustana", "mjc",
    "calib", "cmos_tests", "AASPrep", "ASTR1070",
    "classes", "daily", "forcolin", "latency", "other",
    "failed", "tmp_filter_quarantine",
    "reduced",      # unaudited — always the last resort
)

#: Documented per-target exceptions to the default tree policy, keyed by the
#: *normalized target key* (see :func:`normalize_target`).  Each value is a
#: tree-priority tuple that replaces the default for that target's frames.
#:
#: NGC 5548: the ``macalester`` tree holds 143 unique frames vs 126 under
#: ``rawimage`` (the Dwarf/AGN strategy documents macalester as the superset),
#: so macalester is preferred to keep the campaign's canonical paths in one
#: tree.  The unique-frame COUNT is unaffected — dedup is global — only the
#: chosen canonical path changes.
TREE_PRIORITY_EXCEPTIONS: dict[str, tuple[str, ...]] = {
    "ngc5548": ("macalester", "rawimage") + tuple(
        t for t in DEFAULT_TREE_PRIORITY if t not in ("macalester", "rawimage")
    ),
}


def tree_rank(tree: Optional[str],
              priority: Sequence[str] = DEFAULT_TREE_PRIORITY) -> int:
    """Return the sort rank of ``tree`` under a priority list (lower wins).

    Unknown or missing trees rank after every listed tree, so a frame in an
    unrecognized directory can still be canonical when it is the only copy,
    but never outranks a listed tree.
    """
    try:
        return list(priority).index(tree)
    except ValueError:
        return len(priority)


def choose_canonical(members: Sequence[tuple[str, str]],
                     target_key: Optional[str] = None) -> int:
    """Pick the canonical member of one duplicate group.

    Parameters
    ----------
    members
        Sequence of ``(tree, path)`` pairs — every catalog row that shares
        this group's (basename, jd) identity.
    target_key
        Normalized target key of the group, used to look up a documented
        tree-policy exception (e.g. NGC 5548 prefers ``macalester``).

    Returns
    -------
    int
        Index into ``members`` of the canonical row.  Selection is fully
        deterministic: fewest processing suffixes first (a raw frame beats
        its ``_calibrated`` twin), then best tree rank, then
        lexicographically smallest path (so within one tree the earliest
        night directory wins — the SN 2023ixf wholesale-copied July
        directories lose to the original May/June directories).
    """
    # Look up the effective tree priority: the documented exception if this
    # target has one, the archive-wide default otherwise.
    priority = TREE_PRIORITY_EXCEPTIONS.get(target_key or "", DEFAULT_TREE_PRIORITY)
    # Decorate each member with (processing depth, rank, path, index) and
    # take the minimum — tuple comparison gives us exactly the deterministic
    # ordering we want.  Processing depth comes FIRST (F-1): since duplicate
    # groups now hold a raw frame together with its reduced derivatives, the
    # canonical member must be the one whose pixels nobody has touched, even
    # when a derivative sits in a better-ranked tree (the archive holds 138
    # ``_calibrated`` files inside ``rawimage`` itself).
    best = min(
        (exposure_name(basename_of(path)).n_processing,
         tree_rank(tree, priority), path, idx)
        for idx, (tree, path) in enumerate(members)
    )
    return best[3]


# --------------------------------------------------------------------------
# Basenames and duplicate identity
# --------------------------------------------------------------------------

def basename_of(path: str) -> str:
    """Return the filename component of an archive path.

    The catalog stores POSIX-style relative paths (``tree/.../file.fts``);
    the basename is everything after the last slash.  Paths without a slash
    are already basenames.
    """
    return path.rsplit("/", 1)[-1]


# --------------------------------------------------------------------------
# Exposure identity: one exposure, many filenames (finding F-1, TE.F6)
# --------------------------------------------------------------------------
# THE DEFECT.  S0 v1.0 deduplicated on (basename, jd).  That sees an exact
# copy, and nothing else.  The observatory's own reduction pipeline writes a
# calibrated copy of every frame as ``<raw name>_calibrated.fts.fz`` with the
# header JD untouched; pyscope writes ``<raw name>_wcs`` after a plate solve;
# the Calibrations tree holds uncompressed ``x.fts`` beside ``x.fts.fz``.
# None of those share a basename with their parent, so every one of them was
# counted as a SECOND canonical exposure: 25,442 era-79 and 1,681 era-82
# reduced rows sat at ``is_canonical = 1`` beside the raw frames they were
# made from, and every count of 2026 data was inflated by them.
#
# THE RULE.  Two files are the same exposure when their header JD is
# IDENTICAL and one name is the other plus packaging and processing
# suffixes.  Both halves are required.  The JD half is what makes the rule
# safe: a name like ``xyz_cal_M42`` has a 'cal' token that is not a suffix
# at all, but it can only ever merge with a file that ALSO shares its JD to
# the last bit — which is the twin signature itself.

#: FITS filename extensions, and the compression suffixes that may follow.
FITS_EXTENSIONS: tuple[str, ...] = (".fts", ".fit", ".fits")
COMPRESSION_EXTENSIONS: tuple[str, ...] = (".fz", ".gz")

#: Underscore-delimited tokens the reduction software appends to a raw
#: stem.  Census of the 330,865-row catalog (2026-10-03): ``_calibrated``
#: 27,400 rows, ``_wcs`` 5,160, ``_cal`` 450.  Compared case-insensitively.
#: ``red`` is deliberately absent: it is a FILTER name in this archive
#: (``…_ngc3169_red_60s…``), not a reduction suffix.
PROCESSING_TOKENS: frozenset[str] = frozenset({"calibrated", "cal", "wcs"})

#: A copy counter is at most this many digits (``_1`` .. ``_99``).
_COUNTER_MAX_DIGITS = 2


def strip_compression(basename: str) -> str:
    """Drop a trailing compression suffix (``.fz``/``.gz``) if present.

    The archive stores most files fpack-compressed (``x.fts.fz``) but also
    holds plain ``x.fts`` copies of the same frame; both must reduce to the
    same identity.
    """
    low = basename.lower()
    for ext in COMPRESSION_EXTENSIONS:
        if low.endswith(ext):
            return basename[: -len(ext)]
    return basename


def frame_stem(basename: str) -> str:
    """Return the extension-free identity of a FITS filename.

    ``mlw_V426_Oph_g_5s_2026-06-27T05-40-49.fts.fz`` and the uncompressed
    ``….fts`` both become ``mlw_V426_Oph_g_5s_2026-06-27T05-40-49``.
    A name with no recognized extension is returned unchanged.
    """
    s = strip_compression(basename)
    low = s.lower()
    for ext in FITS_EXTENSIONS:
        if low.endswith(ext):
            return s[: -len(ext)]
    return s


@dataclass(frozen=True)
class ExposureName:
    """One filename, taken apart into the pieces dedup reasons about.

    Attributes
    ----------
    stem
        The extension-free basename (:func:`frame_stem`).
    base
        ``stem`` cut immediately BEFORE its first processing token — the
        name of the raw frame this file claims to derive from.  Equal to
        ``stem`` when the name carries no processing token.
    base_no_counter
        ``base`` with one trailing copy counter removed, or ``None`` when
        ``base`` does not end in one.  Needed for the form
        ``…T05-16-00_1_wcs_calibrated``, where pyscope put a ``_1`` counter
        BETWEEN the raw name and the suffix.  It is only ever used when a
        raw frame of exactly that name exists at the same JD (see
        :func:`exposure_root`), because a trailing number is far more often
        a legitimate frame index (``…_r_1_wcs`` is a copy of ``…_r_1``).
    n_processing
        How many processing tokens the stem carries; 0 = an untouched raw
        name.  This is the "processing depth" the canonical choice sorts on.
    """
    stem: str
    base: str
    base_no_counter: Optional[str]
    n_processing: int


def exposure_name(basename: str) -> ExposureName:
    """Take a FITS basename apart (see :class:`ExposureName`).

    The first token is never treated as a processing token: a name cannot
    BEGIN with a suffix, and observer codes are three lower-case letters
    that could collide (``cal_M13_…``).

    ``mpg_M87_g_90s_2026-06-29T05-16-00_1_wcs_calibrated.fts.fz`` →
    stem ``…T05-16-00_1_wcs_calibrated``, base ``…T05-16-00_1``,
    base_no_counter ``…T05-16-00``, n_processing 2.
    """
    stem = frame_stem(basename)
    tokens = stem.split("_")
    hits = [i for i, tok in enumerate(tokens)
            if i > 0 and tok.lower() in PROCESSING_TOKENS]
    if not hits:
        return ExposureName(stem=stem, base=stem, base_no_counter=None,
                            n_processing=0)
    head = tokens[:hits[0]]
    no_counter: Optional[str] = None
    if (len(head) > 1 and head[-1].isdigit()
            and len(head[-1]) <= _COUNTER_MAX_DIGITS):
        no_counter = "_".join(head[:-1])
    return ExposureName(stem=stem, base="_".join(head),
                        base_no_counter=no_counter, n_processing=len(hits))


#: How a duplicate group member relates to the exposure's root name.
ROOT_IS_SELF = "self"                      # an unprocessed name: its own root
ROOT_BY_SUFFIX = "processing_suffix"       # root = name before the suffix
ROOT_BY_COUNTER = "suffix_counter"         # ... after also dropping a counter
ROOT_ORPHAN = "derivative_no_parent"       # no raw frame of that name exists


def exposure_root(name: ExposureName, has_raw) -> tuple[str, str]:
    """Resolve one file to the ROOT name of the exposure it belongs to.

    Parameters
    ----------
    name
        The :class:`ExposureName` of the file.
    has_raw
        Callable ``stem -> bool``: does an UNPROCESSED file of exactly this
        stem exist AT THE SAME JD as this file?  The build supplies a set
        lookup; the unit tests supply a lambda.

    Returns
    -------
    (root, how)
        ``root`` is the stem every copy of this exposure shares; ``how`` is
        one of the ``ROOT_*`` constants and is the audit trail.

    The ladder, strongest evidence first:

    1.  an unprocessed name is its own root;
    2.  a processed name whose ``base`` exists as a raw frame at the same JD
        roots there — the mechanical ``_calibrated`` rename;
    3.  otherwise, if dropping a copy counter from ``base`` reaches a raw
        frame at the same JD, it roots there;
    4.  otherwise the raw parent is not in the archive: the file roots at
        its own ``base``, so that two derivatives of one missing parent
        (``x_wcs`` and ``x_wcs_cal``) still collapse to one exposure.
    """
    if name.n_processing == 0:
        return name.stem, ROOT_IS_SELF
    if has_raw(name.base):
        return name.base, ROOT_BY_SUFFIX
    if name.base_no_counter is not None and has_raw(name.base_no_counter):
        return name.base_no_counter, ROOT_BY_COUNTER
    return name.base, ROOT_ORPHAN


#: Why a NON-canonical row is a duplicate of its group's canonical row.
DUP_CANONICAL = "canonical"                # the row IS the canonical member
DUP_SAME_BASENAME = "same_basename_jd"     # an exact copy (the v1.0 rule)
DUP_PACKAGING = "packaging"                # same stem, different extension
DUP_PROCESSING = "processing_suffix"       # a reduced/solved derivative
DUP_PROCESSING_COUNTER = "suffix_counter"  # ... reached through a counter


def dup_basis(basename: str, canonical_basename: str, is_canonical: bool,
              how: str = ROOT_IS_SELF,
              canonical_how: str = ROOT_IS_SELF) -> str:
    """Name the evidence that makes one row a duplicate of its group head.

    Every non-canonical row in ``frames`` carries the answer, so a reader
    can count exactly how many frames each dedup rule removed and audit the
    weakest ones (``suffix_counter``) by hand.  ``how``/``canonical_how``
    are the :func:`exposure_root` outcomes of the row and of the head.
    """
    if is_canonical:
        return DUP_CANONICAL
    if basename == canonical_basename:
        return DUP_SAME_BASENAME
    if frame_stem(basename) == frame_stem(canonical_basename):
        return DUP_PACKAGING
    if ROOT_BY_COUNTER in (how, canonical_how):
        return DUP_PROCESSING_COUNTER
    return DUP_PROCESSING


def dup_key(basename: str, jd: Optional[float],
            row_id: Optional[int] = None) -> tuple:
    """Return the global duplicate-group key for a frame.

    Duplicate identity is (exposure root, jd).  Callers pass the frame's
    EXPOSURE ROOT (:func:`exposure_root`) as ``basename`` — for an
    unprocessed, uniquely packaged file that is simply its stem, so the
    original property holds: the same exposure written to two places keeps
    its name and its header JD, while two genuinely different exposures
    never share both.  Frames with no JD (the handful of unreadable files)
    can never be proven duplicates, so each becomes its own singleton
    group, keyed by ``row_id`` to keep them distinct.
    """
    if jd is None or (isinstance(jd, float) and math.isnan(jd)):
        # Unreadable header → no JD → not mergeable with anything.
        return ("__nojd__", row_id)
    return (basename, float(jd))


# --------------------------------------------------------------------------
# Night labels (local-noon-to-noon)
# --------------------------------------------------------------------------

def night_label(jd: Optional[float]) -> Optional[str]:
    """Return the local-noon-to-noon night label for a header JD.

    The night is the calendar date of (JD − 0.7917): subtracting 0.7917 days
    moves the date rollover to 19:00 UT = 12:00 local at Winer (UTC−7), so
    every frame of one observing night — dusk through dawn — shares one
    label, and the label equals the local *evening* date.

    Returns ``None`` for a missing JD (unreadable files).  The JD itself is
    stored untouched elsewhere: it remains the UTC exposure START, and any
    barycentric/mid-exposure correction is stage S3's job, not S0's.
    """
    if jd is None or (isinstance(jd, float) and math.isnan(jd)):
        return None
    # Convert the shifted JD to a UTC datetime via the Unix epoch — pure
    # arithmetic, no astronomy library needed at header precision (~ms).
    days_since_epoch = (jd - NIGHT_SHIFT_DAYS) - UNIX_EPOCH_JD
    moment = datetime(1970, 1, 1, tzinfo=timezone.utc) + timedelta(days=days_since_epoch)
    return moment.date().isoformat()


# --------------------------------------------------------------------------
# Era keying (READOUTM, NAXIS geometry, XBINNING, EGAIN — never filter/date)
# --------------------------------------------------------------------------

def era_key(readoutm: Optional[str], naxis1: Optional[float],
            naxis2: Optional[float], xbinning: Optional[float],
            egain: Optional[float]) -> tuple:
    """Return the hashable camera-era key for one frame.

    The key is exactly the tuple the roadmap mandates — readout mode string,
    detector geometry, binning, and gain — normalized so that header float
    fuzz cannot split one physical configuration into several eras:

    * ``readoutm`` is stripped of surrounding whitespace; missing → ``""``.
    * NAXIS/XBINNING are integers in the headers but floats in the catalog;
      they are cast back to int (missing → ``None``).
    * EGAIN is rounded to :data:`EGAIN_DECIMALS` places (missing → ``None``);
      the two GSENSE gains 1.057 and 1.054 remain distinct at 3 decimals,
      which is intentional — they mark different driver epochs.

    Filter names and dates never enter the key (BeStar lesson: mislabeled
    hrg/lrg frames, era-C repackaging).
    """
    def _num(x, cast):
        # Missing values arrive as None or NaN depending on the reader.
        if x is None or (isinstance(x, float) and math.isnan(x)):
            return None
        return cast(x)

    return (
        (readoutm or "").strip(),
        _num(naxis1, int),
        _num(naxis2, int),
        _num(xbinning, int),
        _num(egain, lambda v: round(float(v), EGAIN_DECIMALS)),
    )


# --------------------------------------------------------------------------
# Target-name normalization (data-driven, auditable)
# --------------------------------------------------------------------------
# Every rule below was derived from patterns actually observed in the
# catalog's 3,014 distinct target_best values; each rule's name appears in
# the aliases table so a reviewer can trace exactly why two names merged.

#: Leading junk prefixes seen in header OBJECT values ('* tet CrB',
#: 'V* KV UMa') — SIMBAD-style object-type markers, not part of the name.
_PREFIX_JUNK_RE = re.compile(r"^(?:\*\s+|v\*\s+)", re.IGNORECASE)

#: Filename-derived observer/date prefixes leaking into target_best, e.g.
#: 'mjcMay01 yzcnc' (observer initials + month + day).  Pattern: 2–4 lower
#: letters, capitalized 3-letter month, 2 digits — as one leading token.
_OBSERVER_PREFIX_RE = re.compile(
    r"^[a-z]{2,4}(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\d{2}$"
)

#: Exposure tokens leaked from filenames: '0-25s', '2-4s', '0p001s', '60s'
#: (the '-' and 'p' both stand for a decimal point in filename encoding).
_EXPOSURE_TOKEN_RE = re.compile(r"^\d+(?:[-p]\d+)?s$", re.IGNORECASE)

#: Filter / grism tokens that ride along with a leaked exposure token
#: ('PHECDA lrg 0-25s', 'Vega 0p001s g 2').  Derived from the catalog's
#: filter column.  A trailing token from this set is stripped ONLY as part
#: of an exposure-leakage tail (see :func:`_strip_exposure_leakage`) — never
#: on its own, so real names ending in a short word are safe.
_FILTER_TOKENS = frozenset(
    {"hrg", "lrg", "ha", "og", "hagrism", "oggrism",
     "lum", "red", "green", "blue", "empty"}
    | set("griubvlwohxyz") | {"6", "1"}
)

#: Bare frame-index tokens that trail some leaked names ('Vega 0p1s hrg 7').
_INDEX_TOKEN_RE = re.compile(r"^\d{1,2}$")

#: Series suffixes seen on CV filename-parsed names ('ST-LMi-y-series',
#: 'STLMi-z-series') — campaign bookkeeping, not part of the target name.
_SERIES_SUFFIX_RE = re.compile(r"[-_\s](?:[a-z])[-_\s]series\s*$", re.IGNORECASE)

#: Constellation-genitive spelling variants observed in the catalog
#: ('RR Lyrae' vs 'RR Lyr').  Applied to the final token only.
GENITIVE_MAP: dict[str, str] = {
    "lyrae": "lyr",
}

#: Explicit known-synonym table, applied to the *normalized key*.  Entries
#: merge two names that no pure string rule can relate; every entry MUST
#: pass the coordinate-cone gate at build time (checked, recorded, and the
#: merge is refused if the check fails).  Keep this list short and famous.
SYNONYM_TABLE: dict[str, str] = {
    # Vega's Bayer designation vs its proper name (BeStar standard series).
    "alphalyr": "vega",
    # M101 vs its NGC number (SN 2023ixf host field; frames labeled either).
    "ngc5457": "m101",
    # M101's common name.  The observer typed 'pinwheel galaxy' on the
    # 2026-03-21/22 nights: 245 catalog rows, 140 of them canonical rawimage
    # science (35 each g/r/i at 30 s + 35 ha at 90 s).  No string rule can
    # relate a common name to a Messier number, so it belongs here; the cone
    # gate passes at 0.031 deg against the natives of 'm101' (measured
    # 2026-08-18).  Until this merge landed, the deepest post-fade M101-field
    # epoch in the archive was invisible to the SN project — the very
    # material its §3.4 template plan says it lacks.
    "pinwheelgalaxy": "m101",
    # SN 2023ixf with a sequence counter fused to the name.  Two frames on
    # 2023-11-28 (the strategy's first named post-fade template epoch) carry
    # OBJECT='2023ixf1'/'2023ixf2'; both sit within 0.05 deg of the SN, and
    # the '2023ixf1' cone gate passes at 0.023 deg (measured 2026-08-18).
    # Listed explicitly rather than handled by a "strip a trailing digit"
    # rule: such a rule would have to know which stems are real target keys,
    # which is corpus knowledge this pure normalizer deliberately lacks.
    "2023ixf1": "2023ixf",
    "2023ixf2": "2023ixf",
}


@dataclass(frozen=True)
class NormalizedName:
    """Result of normalizing one raw target_best value.

    Attributes
    ----------
    key
        The canonical grouping key (lowercase, separator-free).  Two raw
        names alias the same target exactly when their keys are equal.
    cleaned
        The human-readable cleaned name (junk stripped, whitespace
        collapsed) — used to pick a display name for the group.
    rules
        Tuple of rule names that actually fired, in application order —
        the audit trail written into the aliases table.
    pre_synonym_key
        The key as it stood BEFORE the synonym table fired (equal to
        ``key`` when no synonym applied).  The build uses it to run the
        coordinate-cone gate on synonym merges and to undo a merge that
        fails the gate.
    """
    key: Optional[str]
    cleaned: Optional[str]
    rules: tuple[str, ...] = field(default_factory=tuple)
    pre_synonym_key: Optional[str] = None


def _strip_exposure_leakage(tokens: list[str]) -> tuple[list[str], bool]:
    """Drop a trailing exposure-leakage tail from a token list.

    Filename-parsed targets leak trailing tokens in two observed orders:
    ``<name> <filter> <exp>`` (``tet CrB hrg 2-4s``) and
    ``<name> <exp> <filter> [<index>]`` (``Vega 0p001s lrg 5``).  We remove
    the longest trailing run made up ONLY of exposure / filter / index
    tokens, and only when that run contains at least one exposure token —
    the exposure token is the smoking gun; filter-like or digit-like tokens
    alone never trigger a strip (protects 'ZTFJ… 5', 'ksi UMa B', 'AG LMi').
    At least one head token must survive.
    """
    def _droppable(tok: str) -> bool:
        low = tok.lower()
        return bool(_EXPOSURE_TOKEN_RE.match(low)) \
            or low in _FILTER_TOKENS \
            or bool(_INDEX_TOKEN_RE.match(low))

    # Walk candidate split points left to right; the first split whose whole
    # tail is droppable AND contains an exposure token wins (longest tail).
    for start in range(1, len(tokens)):
        tail = tokens[start:]
        if all(_droppable(t) for t in tail) and \
           any(_EXPOSURE_TOKEN_RE.match(t.lower()) for t in tail):
            return tokens[:start], True
    return tokens, False


def normalize_target(raw: Optional[str]) -> NormalizedName:
    """Normalize one raw ``target_best`` value into an alias-group key.

    The pipeline of pure string rules (each recorded when it fires):

    1.  ``blank``            — None/empty input → key ``None``.
    2.  ``whitespace``       — strip ends, collapse internal runs.
    3.  ``prefix_junk``      — drop leading ``* `` / ``V* `` markers.
    4.  ``observer_prefix``  — drop a leading ``mjcMay01``-style token.
    5.  ``series_suffix``    — drop trailing ``-y-series`` / ``-z-series``.
    6.  ``exposure_tokens``  — drop a trailing exposure-leakage tail
        (``PHECDA lrg 0-25s`` → ``PHECDA``).
    7.  ``genitive``         — canonical constellation spelling
        (``RR Lyrae`` → ``RR Lyr``).
    8.  ``casefold``         — lowercase (``T CrB`` ≡ ``t crb``).
    9.  ``separators``       — drop underscores and spaces entirely, and
        drop hyphens *only when followed by a letter* (``ST-LMi`` →
        ``stlmi``) — a hyphen before a digit is kept because it may be a
        declination sign (``ZTFJ082835 05-052702 1``).
    10. ``synonym``          — explicit table lookup on the finished key
        (``alphalyr`` → ``vega``); the build gates every such merge with
        the 0.2-degree coordinate cone before accepting it.

    Fuzzy string similarity is deliberately absent: a merge happens through
    a named rule above or not at all.
    """
    rules: list[str] = []

    # Rule 1: blank / missing names form no alias group at all.
    if raw is None or not str(raw).strip():
        return NormalizedName(key=None, cleaned=None, rules=("blank",),
                              pre_synonym_key=None)
    text = str(raw)

    # Rule 2: whitespace hygiene ('  eta UMa', 'RR  LYR').
    collapsed = " ".join(text.split())
    if collapsed != text:
        rules.append("whitespace")
    text = collapsed

    # Rule 3: SIMBAD-style leading markers ('* tet CrB', 'V* KV UMa').
    stripped = _PREFIX_JUNK_RE.sub("", text)
    if stripped != text:
        rules.append("prefix_junk")
    text = stripped

    # Rule 4: observer/date filename prefixes ('mjcMay01 yzcnc').
    tokens = text.split(" ")
    if len(tokens) >= 2 and _OBSERVER_PREFIX_RE.match(tokens[0]):
        tokens = tokens[1:]
        rules.append("observer_prefix")

    # Rule 5: campaign series suffixes ('ST-LMi-y-series').
    joined = " ".join(tokens)
    unsuffixed = _SERIES_SUFFIX_RE.sub("", joined)
    if unsuffixed != joined:
        rules.append("series_suffix")
    tokens = unsuffixed.split(" ")

    # Rule 6: exposure-token leakage ('PHECDA lrg 0-25s', 'Vega 0p1s hrg 7').
    tokens, stripped_exp = _strip_exposure_leakage(tokens)
    if stripped_exp:
        rules.append("exposure_tokens")

    # Rule 7: constellation genitive variants ('RR Lyrae' → 'RR Lyr').
    if tokens and tokens[-1].lower() in GENITIVE_MAP:
        tokens = tokens[:-1] + [GENITIVE_MAP[tokens[-1].lower()]]
        rules.append("genitive")

    cleaned = " ".join(tokens)

    # Rule 8: case folding.
    lowered = cleaned.lower()
    if lowered != cleaned:
        rules.append("casefold")

    # Rule 9: separator unification.  Spaces/underscores vanish; hyphens
    # vanish only before letters (keep declination signs before digits).
    key = re.sub(r"-(?=[a-z])", "", lowered)
    key = key.replace("_", "").replace(" ", "")
    if key != lowered:
        rules.append("separators")

    # Rule 10: explicit synonym table (cone-gated by the build).
    pre_synonym_key = key
    if key in SYNONYM_TABLE:
        key = SYNONYM_TABLE[key]
        rules.append("synonym")

    return NormalizedName(key=key, cleaned=cleaned, rules=tuple(rules),
                          pre_synonym_key=pre_synonym_key)


# --------------------------------------------------------------------------
# Hardware-state columns from the header re-scrape (finding F-2, DE.F5)
# --------------------------------------------------------------------------
# ``rescan_geometry.py hdr-run`` stores the hardware cards of every file AS
# THE HEADER SPELLS THEM (text; NULL = no such card, '' = blank card).  The
# functions below turn that text into the typed ``frames`` columns.  Typing
# is done here rather than in the scraper because it needs judgement the
# scraper must not exercise: ``GAIN`` is the string '4x' on the Andor iKon
# (a preamp setting) and the number 100 on the ASI (a gain register), so it
# is kept twice — once verbatim, once numeric where a number exists.

#: Normalized camera identity, from INSTRUME.  Four cameras have been on
#: the telescope (TE.F1); INSTRUME names the DRIVER, which changed twice
#: for the QHY600 without the camera changing (MaxIm ASCOM capture ->
#: pyscope native), so three strings map to one camera.
CAMERA_OF_INSTRUME: dict[str, str] = {
    "DL Imaging": "AC4040",                      # SBIG Aluma AC4040 (GSENSE)
    "Andor CCD/EMCCD (SDK2)": "iKon",            # Andor iKon-L 936
    "ASI Camera (1)": "ASI",                     # ZWO ASI (IMX455)
    "QHYCCD-Cameras-Capture": "QHY600",          # QHY600 via MaxIm/ASCOM
    "QHY600Pro": "QHY600",                       # QHY600 via pyscope
}

#: CAMNAME prefixes that identify a camera when INSTRUME is absent.  The
#: 854 pyscope-native frames of 2026-06-28/29 carry no INSTRUME card at all;
#: their CAMNAME reads 'QHY600MPCIE-<serial>'.  MaxIm-era frames all say
#: CAMNAME = 'ASCOM', which identifies nothing and is not listed.
CAMERA_OF_CAMNAME_PREFIX: tuple[tuple[str, str], ...] = (
    ("QHY600", "QHY600"),
)


def camera_id(instrume: Optional[str], camname: Optional[str]) -> Optional[str]:
    """Normalized camera identity of a frame, or ``None`` when unknown.

    INSTRUME decides when it is a known string; otherwise a CAMNAME prefix
    may.  An UNKNOWN string returns ``None`` rather than a guess — a frame
    with no provable camera gets no mechanical epoch, and the report counts
    such frames instead of filing them under the nearest camera.
    """
    key = (instrume or "").strip()
    if key in CAMERA_OF_INSTRUME:
        return CAMERA_OF_INSTRUME[key]
    name = (camname or "").strip()
    for prefix, cam in CAMERA_OF_CAMNAME_PREFIX:
        if name.startswith(prefix):
            return cam
    return None


def card_float(text) -> Optional[float]:
    """Numeric value of a re-scraped card, or ``None``.

    ``None`` for an absent card, a blank card, and a card whose value is
    not a number (the iKon's ``GAIN = '4x'``) — never 0.0, which is a real
    set-point and a real cooler power.
    """
    if text is None:
        return None
    t = str(text).strip()
    if not t:
        return None
    try:
        v = float(t)
    except ValueError:
        return None
    return None if math.isnan(v) else v


_WHEEL_SLOT_RE = re.compile(r"'([^']*)'")


def wheel_slots(fwallnam: Optional[str]) -> Optional[tuple[str, ...]]:
    """The filter-wheel slot names, in slot order, from a FWALLNAM value.

    The card holds a Python tuple/list literal written by pyscope —
    ``('g', 'lrg', 'r', 'i', 'ha', 'hrg', 'empty')`` — so the slot names
    are simply its quoted items.  Case is PRESERVED: the 2024-11-08 change
    from 'HaGrism' to 'hrg' is a reloaded configuration, and whether a
    rename is also a physical reload is for the reader of ``mech_epoch`` to
    judge with the boundary in front of them, not for a case-fold to hide.

    Returns ``None`` for an absent or blank card (the whole AC4040 era:
    MaxIm wrote no wheel map), and for a value with no quoted item.
    """
    if fwallnam is None or not str(fwallnam).strip():
        return None
    slots = tuple(_WHEEL_SLOT_RE.findall(str(fwallnam)))
    return slots or None


def wheel_map_key(fwallnam: Optional[str]) -> Optional[str]:
    """Canonical one-line spelling of a wheel map: slots joined by ``|``.

    This is what ``frames.fwallnam`` stores and what the mechanical-epoch
    segmentation compares; two headers describe the same wheel exactly
    when their keys are equal.
    """
    slots = wheel_slots(fwallnam)
    return None if slots is None else "|".join(slots)


#: The typed hardware columns S0 adds to ``frames``, with the meaning of
#: each (rendered into the report and the data README — one source).
HARDWARE_FRAME_COLUMNS: tuple[tuple[str, str], ...] = (
    ("camera", "normalized camera identity (AC4040 / iKon / ASI / QHY600)"),
    ("hdr_gain", "GAIN card exactly as written ('4x' on the iKon, 100 on "
                 "the ASI, 56 on the QHY600)"),
    ("hdr_gain_num", "GAIN as a number where it is one, else NULL"),
    ("hdr_offset", "OFFSET card (camera bias-offset setting), numeric"),
    ("set_temp", "SET-TEMP: cooler set-point, deg C"),
    ("coolpowr", "COOLPOWR: cooler power, percent"),
    ("focpos", "FOCPOS: pyscope's focuser position — the true one; "
               "frames.focuspos is MaxIm's FOCUSPOS, which sticks (TE.F8)"),
    ("flipstat", "FLIPSTAT: MaxIm's software flip state; '' = card present "
                 "and blank, NULL = no card"),
    ("telpier", "TELPIER: mount pier side"),
    ("fwpos", "FWPOS: filter-wheel slot index of this exposure"),
    ("fwallnam", "FWALLNAM: the wheel's slot map, slots joined by '|'"),
    ("camname", "CAMNAME: camera name (identifies the pyscope-native QHY)"),
    ("fwname", "FWNAME: filter-wheel hardware name"),
    ("xpixsz", "XPIXSZ: binned pixel pitch, microns"),
    ("hdr_scanned", "1 = this file's header was re-scraped; 0 = it was not "
                    "(every hardware column is then NULL)"),
)


def hardware_columns(card: dict) -> dict:
    """Typed ``frames`` columns for one row of the header re-scrape.

    ``card`` maps the re-scrape's column names (``h_gain``, ``h_set_temp``,
    …, plus ``h_scanned``) to their stored text; a row that was never
    re-scraped passes an empty dict and gets NULLs with ``hdr_scanned = 0``.
    Text columns keep the blank/absent distinction ('' vs None); numeric
    columns are ``None`` for both.
    """
    def _text(key):
        v = card.get(key)
        return None if v is None else str(v).strip()

    fwpos = card_float(card.get("h_fwpos"))
    return {
        "camera": camera_id(card.get("h_instrume"), card.get("h_camname")),
        "hdr_gain": _text("h_gain"),
        "hdr_gain_num": card_float(card.get("h_gain")),
        "hdr_offset": card_float(card.get("h_offset")),
        "set_temp": card_float(card.get("h_set_temp")),
        "coolpowr": card_float(card.get("h_coolpowr")),
        "focpos": card_float(card.get("h_focpos")),
        "flipstat": _text("h_flipstat"),
        "telpier": _text("h_telpier"),
        "fwpos": None if fwpos is None else int(fwpos),
        "fwallnam": wheel_map_key(card.get("h_fwallnam")),
        "camname": _text("h_camname"),
        "fwname": _text("h_fwname"),
        "xpixsz": card_float(card.get("h_xpixsz")),
        "hdr_scanned": 1 if card.get("h_scanned") else 0,
    }


# --------------------------------------------------------------------------
# Pointing math
# --------------------------------------------------------------------------

def angular_separation_deg(ra1: float, dec1: float,
                           ra2: float, dec2: float) -> float:
    """Great-circle separation between two sky positions, in degrees.

    Uses the Vincenty formula, which is numerically stable at every
    separation (haversine loses precision near 180 degrees; the plain
    spherical cosine law loses it near 0).  Inputs and output in degrees.
    """
    lam1, phi1 = math.radians(ra1), math.radians(dec1)
    lam2, phi2 = math.radians(ra2), math.radians(dec2)
    dlam = lam2 - lam1
    # Vincenty numerator: length of the cross-product of the unit vectors.
    num = math.hypot(
        math.cos(phi2) * math.sin(dlam),
        math.cos(phi1) * math.sin(phi2)
        - math.sin(phi1) * math.cos(phi2) * math.cos(dlam),
    )
    # Vincenty denominator: dot product of the unit vectors.
    den = math.sin(phi1) * math.sin(phi2) \
        + math.cos(phi1) * math.cos(phi2) * math.cos(dlam)
    return math.degrees(math.atan2(num, den))


def median_radec(ras: Iterable[float],
                 decs: Iterable[float]) -> tuple[float, float]:
    """Median sky position of a set of coordinates, RA-wrap-aware.

    Declination is a plain median.  Right ascension needs care when the
    samples straddle 0h/24h (e.g. 359.9 and 0.1 must average to 0.0, not
    180.0): if the raw RA spread exceeds 180 degrees we lift the low RAs by
    360, take the median, and fold back into [0, 360).
    """
    ra_list = [float(r) for r in ras]
    dec_list = [float(d) for d in decs]
    if not ra_list:
        raise ValueError("median_radec needs at least one coordinate pair")
    if max(ra_list) - min(ra_list) > 180.0:
        # The set straddles the 0/360 seam — unwrap before taking a median.
        ra_list = [r + 360.0 if r < 180.0 else r for r in ra_list]
    ra_med = median(ra_list) % 360.0
    return ra_med, median(dec_list)


#: How a target's reference position was derived.  Recorded per target so
#: the pointing audit says which evidence it rests on instead of leaving a
#: silent NULL behind.
POINTING_REF_SOLVED = "plate_solved"     # median of plate-solved frames
POINTING_REF_HEADER = "header_median"    # median of header coordinates
POINTING_REF_UNSUPPORTED = "header_unsupported"   # headers disagree: refused
POINTING_REF_NONE = "none"               # too few coordinates: no reference

#: Fraction of a target's own frames that must lie within
#: :data:`POINTING_OUTLIER_DEG` of a HEADER-derived reference before that
#: reference is accepted.  See :func:`pointing_reference` for why.
POINTING_REF_SUPPORT_FRACTION = 0.5


def pointing_reference(ras: Iterable[Optional[float]],
                       decs: Iterable[Optional[float]],
                       pltsolvd: Iterable[object],
                       min_frames: int = MIN_SOLVED_FOR_REFERENCE,
                       ) -> tuple[Optional[float], Optional[float], str]:
    """Reference position of ONE target, with the basis that produced it.

    THE NGC 5548 LESSON (2026-08-18 review).  The original rule was "median
    of the target's plate-solved canonical frames, if at least
    ``MIN_SOLVED_FOR_REFERENCE`` exist".  NGC 5548 has **zero** plate-solved
    frames out of 279, so it got no reference, so every offset was NaN, so
    ``pointing_gt1deg`` could never fire — on precisely the target whose
    strategy document rules an entire night "unusable, pointing-failure" and
    whose staging rule promises that "its pointing flag, not its absence, is
    the record".  The flagging machinery was alive (21k frames flagged
    archive-wide) and blind exactly where it was needed.

    So: prefer plate-solved coordinates, and fall back to the median of the
    target's HEADER coordinates when too few solves exist.  The fallback is
    weaker evidence — header RA/Dec is where the mount believed it was
    pointing — but a median over many frames is robust to a MINORITY of
    mispointed ones, which is the failure mode being detected: NGC 5548's
    269 on-target rows median cleanly at Dec +25.1 and the 10 rows of the
    parked-mount night stand out at ~8 deg.

    THE SUPPORT TEST (why the fallback can also refuse).  A plain median is
    NOT robust when a target's headers split into two clusters of similar
    size: the median then lands between them and flags everything, including
    the frames that were on target.  Measured on the real archive: Mizar has
    30 canonical frames — 10 with no coordinates, 10 at its true position,
    and 10 carrying a stale RA/Dec card also seen on Denebola, Gamma Vir and
    five other targets.  A bare median put the reference between the two
    clusters and flagged all 20.  So a header-derived reference is accepted
    only when at least :data:`POINTING_REF_SUPPORT_FRACTION` of the target's
    coordinate-bearing frames lie within :data:`POINTING_OUTLIER_DEG` of it;
    otherwise the basis is ``header_unsupported`` and NO offsets are
    computed — the honest answer when the headers contradict each other and
    nothing external can say which cluster is right.  The test is applied
    only to the fallback: plate-solved coordinates are trusted evidence and
    keep their original, unchanged behaviour.

    Returns ``(ra, dec, basis)``; ``(None, None, ...)`` when no reference is
    established, with the basis saying which way it failed.
    """
    solved_ra: list[float] = []
    solved_dec: list[float] = []
    all_ra: list[float] = []
    all_dec: list[float] = []
    for ra, dec, solved in zip(ras, decs, pltsolvd):
        # Skip missing coordinates, including NaN (which fails == itself).
        if ra is None or dec is None or ra != ra or dec != dec:
            continue
        all_ra.append(float(ra))
        all_dec.append(float(dec))
        if solved == 1 or solved is True:
            solved_ra.append(float(ra))
            solved_dec.append(float(dec))
    if len(solved_ra) >= min_frames:
        ra0, dec0 = median_radec(solved_ra, solved_dec)
        return ra0, dec0, POINTING_REF_SOLVED
    if len(all_ra) >= min_frames:
        ra0, dec0 = median_radec(all_ra, all_dec)
        # The support test: how many of this target's own frames actually
        # sit at the position their median claims?
        near = sum(1 for ra, dec in zip(all_ra, all_dec)
                   if angular_separation_deg(ra, dec, ra0, dec0)
                   <= POINTING_OUTLIER_DEG)
        if near >= POINTING_REF_SUPPORT_FRACTION * len(all_ra):
            return ra0, dec0, POINTING_REF_HEADER
        return None, None, POINTING_REF_UNSUPPORTED
    return None, None, POINTING_REF_NONE


# --------------------------------------------------------------------------
# QC flags
# --------------------------------------------------------------------------

def qc_flags(error: Optional[str], exptime: Optional[float],
             airmass: Optional[float], jd: Optional[float],
             ra_deg: Optional[float], target_key: Optional[str],
             pointing_offset_deg: Optional[float]) -> str:
    """Compute the comma-joined QC-flag string for one frame.

    Flags mark; they never delete — every catalog row keeps its manifest
    row, and downstream stages decide what each flag disqualifies.

    Tokens (alphabetical by cause, stable order in the output):

    * ``header_error``      — the cataloger could not read the file at all.
    * ``exptime_nonpos``    — EXPTIME ≤ 0 (fine for a bias, fatal for a
      light frame; the flag records the fact, the consumer applies context).
    * ``airmass_garbage``   — header airmass < 1, > 10, or the −999
      sentinel.  Flagged only: recomputation is S3's job (ROADMAP conv. 7).
    * ``no_jd``             — no header JD → no night label, no dedup proof.
    * ``no_coords``         — no usable RA/Dec → no pointing validation.
    * ``blank_target``      — no target name → excluded from alias groups.
    * ``pointing_gt1deg``   — pointing offset beyond
      :data:`POINTING_OUTLIER_DEG` (the T CrB grism failure mode).
    """
    def _missing(x) -> bool:
        return x is None or (isinstance(x, float) and math.isnan(x))

    flags: list[str] = []
    if error is not None and str(error).strip():
        flags.append("header_error")
    if not _missing(exptime) and float(exptime) <= 0:
        flags.append("exptime_nonpos")
    if not _missing(airmass) and \
       (float(airmass) < AIRMASS_MIN or float(airmass) > AIRMASS_MAX):
        flags.append("airmass_garbage")
    if _missing(jd):
        flags.append("no_jd")
    if _missing(ra_deg):
        flags.append("no_coords")
    if target_key is None:
        flags.append("blank_target")
    if not _missing(pointing_offset_deg) and \
       float(pointing_offset_deg) > POINTING_OUTLIER_DEG:
        flags.append("pointing_gt1deg")
    return ",".join(flags)


# --------------------------------------------------------------------------
# Strategy-claimed reference values (Output C / report section 7)
# --------------------------------------------------------------------------
# Source: each project's ANALYSIS_STRATEGY.md, rev 2026-08-16, read from the
# repo on 2026-08-17.  These are the numbers the five panels PUBLISHED in
# their inventory tables; S0 reproduces each one from the manifest and the
# report shows the diff.
#
# Metric definitions (implemented in the build script — any metric name not
# listed here raises at build time, so this comment cannot go stale
# silently).  Every strategy's own accounting rule is "the canonical tree
# only" (rawimage — or, for NGC 5548, macalester per the documented tree
# exception), so the like-for-like manifest number counts canonical frames
# from the target's PRIMARY tree; the build also records the fully global
# canonical count in a separate column so nothing is hidden:
#
#   rows_all_trees : raw catalog rows, all trees, Light frames (no dedup) —
#                    what a panel quotes when it counts rows.
#   unique_light   : canonical (deduplicated) error-free Light frames.
#   grism_light    : unique_light restricted to filters hrg/lrg.
#   grism4_light   : unique_light restricted to hrg/lrg/HaGrism/OGGrism.
#
# Each entry: (project, target_key, metric,
#              claimed_frames, claimed_nights, where_the_claim_lives).
STRATEGY_CLAIMS: tuple[tuple, ...] = (
    # --- TCrB_Monitoring/ANALYSIS_STRATEGY.md section 3 -------------------
    ("TCrB_Monitoring", "tcrb", "rows_all_trees", 862, 86,
     "sec.3: '862 T CrB light-frame rows' across trees; 86 distinct nights"),
    ("TCrB_Monitoring", "tcrb", "unique_light", 471, 86,
     "sec.3: '471 unique light frames = 414 T CrB + 57 t crb' — counted "
     "rawimage ROWS; global (basename, jd) dedup collapses 69 "
     "within-rawimage duplicate copies the panel's row count kept"),
    ("TCrB_Monitoring", "tcrb", "grism_light", 247, 60,
     "sec.1: '247 grism spectra on 60 nights' (lrg 147 + hrg 100)"),
    ("TCrB_Monitoring", "tetcrb", "grism_light", 403, 40,
     "sec.3: theta CrB calibrator grism series '403 unique rawimage "
     "frames, 40 nights' (hrg + lrg)"),
    # --- CV_TimeSeries/ANALYSIS_STRATEGY.md section 3.1 -------------------
    ("CV_TimeSeries", "stlmi", "unique_light", 3157, 39,
     "sec.3.1 table: ST LMi 3,157 raw light frames / 39 nights (the iKon "
     "tree holds a further 109 unique frames outside rawimage)"),
    ("CV_TimeSeries", "yzcnc", "unique_light", 1920, 26,
     "sec.3.1 table: YZ Cnc 1,920 frames / 26 nights"),
    ("CV_TimeSeries", "vvpup", "unique_light", 1353, 29,
     "sec.3.1 table: VV Pup 1,353 frames / 29 nights"),
    ("CV_TimeSeries", "euuma", "unique_light", 993, 32,
     "sec.3.1 table: EU UMa 993 frames / 32 nights (reduced/ holds 257 "
     "RENAMED copies whose JDs all collide with rawimage frames)"),
    ("CV_TimeSeries", "anuma", "unique_light", 1279, 14,
     "sec.2 Q5: AN UMa 1,279 raw light frames, 14 nights"),
    # --- SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md section 3.1 ------------
    ("SN2023ixf_LightCurve", "2023ixf", "unique_light", 1052, 35,
     "sec.3.1: 1,052 unique light frames (3,018 rawimage rows); the '35 "
     "nights' includes the two saturated first epochs labeled NGC5457 "
     "(May 20) and M101 (May 21), which carry other target names"),
    # --- BeStar_Grism/ANALYSIS_STRATEGY.md section 3.2 --------------------
    ("BeStar_Grism", "spica", "grism4_light", 728, 29,
     "sec.3.2 table: Spica 728 raw frames / 29 nights"),
    ("BeStar_Grism", "phecda", "grism4_light", 333, 40,
     "sec.3.2 table: Phecda 333 / 40"),
    ("BeStar_Grism", "phileo", "grism4_light", 261, 39,
     "sec.3.2 table: phi Leo 261 / 39"),
    ("BeStar_Grism", "lameri", "grism4_light", 252, 47,
     "sec.3.2 table: lambda Eri 252 / 47"),
    ("BeStar_Grism", "vega", "grism4_light", 447, 19,
     "sec.3.2 table: Vega era C 367/18 + era A ladder 80/1 (Alpha Lyr "
     "alias, cone-gated synonym merge)"),
    # --- DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md section 3 ------------
    ("DwarfGalaxy_AGN_Survey", "ngc5548", "unique_light", 133, 15,
     "sec.3: 133 unique on-target frames / 15 nights AFTER excluding the "
     "10-frame mispointed night 2023-03-25; the manifest keeps all 143 "
     "unique frames / 16 nights (mispointing is a pointing flag, not a "
     "dedup fact; primary tree macalester per the documented exception)"),
    ("DwarfGalaxy_AGN_Survey", "ngc5238", "unique_light", 528, 21,
     "sec.3: NGC 5238 ~528 unique frames, 21 nights"),
    ("DwarfGalaxy_AGN_Survey", "__dw_survey__", "unique_light", 499, None,
     "sec.3: Dw survey unique totals L 323 + Halpha 89 + R 87 = 499 "
     "(19 Dw* fields combined)"),
)


def primary_tree(target_key: Optional[str]) -> str:
    """Return the canonical tree for a target under the tree policy.

    ``rawimage`` for everyone except the documented exceptions
    (NGC 5548 → ``macalester``).  This is the tree each strategy's own
    accounting rule counted, so reconciliation compares like with like.
    """
    return TREE_PRIORITY_EXCEPTIONS.get(target_key or "",
                                        DEFAULT_TREE_PRIORITY)[0]

#: Filters counted as the classic slitless-grism pair (T CrB spectra).
GRISM_FILTERS: frozenset[str] = frozenset({"hrg", "lrg"})

#: Filters counted as grism-family for the BeStar selection (its section
#: 3.1 rule: hrg, lrg, HaGrism, OGGrism).
GRISM4_FILTERS: frozenset[str] = frozenset({"hrg", "lrg", "HaGrism", "OGGrism"})
