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
"""

from __future__ import annotations

import math
import re
from typing import Iterable, Optional, Sequence

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
REDUCED_SUFFIX_TOKENS: frozenset[str] = frozenset({"calibrated", "cal", "wcs"})

#: FITS extensions (compression suffix handled separately).
_FITS_EXTENSIONS = (".fts", ".fit", ".fits")
_COMPRESSION_EXTENSIONS = (".fz", ".gz")

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

def strip_compression(basename: str) -> str:
    """Drop a trailing compression suffix (``.fz``/``.gz``) if present.

    The archive stores most files fpack-compressed (``x.fts.fz``) but the
    reduced tree also holds plain ``x.fts`` copies of the same frame; both
    must reduce to the same identity.
    """
    low = basename.lower()
    for ext in _COMPRESSION_EXTENSIONS:
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
    for ext in _FITS_EXTENSIONS:
        if low.endswith(ext):
            return s[: -len(ext)]
    return s


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
