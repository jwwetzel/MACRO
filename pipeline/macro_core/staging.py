"""Pure S0c staging-manifest logic: the per-project working set, without copies.

WHY THIS MODULE EXISTS (the no-copy law)
----------------------------------------
S0 exists because copies proliferated: the archive's 330,865 catalog rows
collapse to 198,294 canonical frames — 132,571 rows (40%) are duplicate
copies of frames that already existed somewhere else (wholesale re-copied
night directories, consortium mirrors, renamed reduced-tree files).  Copying
science frames into per-project working directories would restart exactly
that disease.  The committee's ruling, implemented here:

    **No project ever copies a frame.  The staging manifest IS the working
    set**: a script-emitted, provenance-complete list of frames per project;
    every downstream stage reads the immutable archive directly *through*
    the manifest (relative path + absolute path per row).

Everything in this module is a pure function or a data table — no I/O, no
globals mutated.  The build script (``pipeline/scripts/build_s0c_staging.py``)
wires it to the manifest database; the unit tests
(``pipeline/tests/test_staging.py``) exercise every rule, including the
rules that MUST reject (reduced-tree rows, non-canonical rows, grism frames
in a photometric selection, ...).

WHAT A STAGING MANIFEST CONTAINS
--------------------------------
For each of the five projects:

* **science rows** — the canonical Light frames the project's published
  selection rule claims (tree policy already handled upstream by S0's
  ``is_canonical``; the reduced tree is excluded outright because its
  canonical rows are renamed copies — the S0b lesson).
* **calibration rows** — for every camera era the project's science
  touches, ALL of that era's calibration frames from the S0b census
  (bias/dark/flat, raw and recovered ``Calibrations/`` masters alike), with
  an explicit ``match_basis`` recording *why* each row is present.  Era is
  the only match key used at staging time (kind/exposure/filter narrowing
  is each stage's job, with the S0b coverage matrix as its guide) — staging
  over-includes deliberately so no stage has to go back to the census.

INTEGRITY SURROGATE (why size, not a checksum)
----------------------------------------------
Each row carries ``size_bytes`` from the S0 catalog scan — a cheap tripwire
against truncated or replaced files at read time.  A real content checksum would
mean re-hashing the 3.3 TiB archive; that is a separate, deliberate decision
(and a multi-day I/O job), not something to smuggle into a staging build.
:data:`CHECKSUM_NOTE` states this in every README the build writes.
"""

from __future__ import annotations

import posixpath
import re
from dataclasses import dataclass, field
from typing import Optional

from . import inventory as inv
from . import manifest as mf

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Root of the immutable archive.  Every relative ``path`` in the manifest
#: resolves against this directory; nothing under it is ever written.
DEFAULT_ARCHIVE_ROOT = "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive"

#: The one-line justification for the size-only integrity surrogate,
#: interpolated into every ``data/README.md`` and the report.
CHECKSUM_NOTE = (
    "size_bytes is an integrity SURROGATE, not a checksum: it comes from the "
    "S0 catalog scan and catches truncation/replacement at read time.  A "
    "content hash would require re-reading the full 3.3 TiB archive — that "
    "is a separate archive-custody decision, not part of a staging build."
)

#: Column order of every stage table and stage_manifest.csv — one schema for
#: science and calibration rows, so a stage can stream the whole file.
STAGE_CSV_COLUMNS: tuple[str, ...] = (
    "path",                 # archive-relative POSIX path (the identity)
    "abs_path",             # archive_root + path (what a stage open()s)
    "role",                 # 'science' | 'bias' | 'dark' | 'flat' | 'master_*'
    "match_basis",          # why the row is here (see MATCH_BASIS_*)
    "tree",                 # top-level archive tree the canonical copy lives in
    "era_id",               # S0 pinned camera-era registry id
    "night",                # local-noon-to-noon night label
    "jd",                   # header JD — UTC exposure START (S3 owns BJD_TDB)
    "filter",               # header/filename filter string, as cataloged
    "exptime",              # header EXPTIME, seconds
    "canonical_target",     # S0 alias-merged display name (science rows)
    "target_key",           # S0 normalized target key (science rows)
    "dup_group",            # S0 global duplicate-group id
    "qc_flags",             # S0 comma-joined QC flags (flags mark, never drop)
    "pointing_offset_deg",  # offset from the target's reference position
    "size_bytes",           # integrity surrogate — see CHECKSUM_NOTE
    "obs_rowid",            # catalog row id: the join key back to obs/frames
    "stage_build_id",       # which S0c build emitted this row
)

#: match_basis for science rows: the row is present because the project's
#: published selection rule (encoded in PROJECT_SELECTIONS) claims it.
MATCH_BASIS_SCIENCE = "selection_rule"

#: match_basis for calibration rows: the frame's S0 era_id equals an era the
#: project's science frames touch.  Masters are flagged by their role
#: (``master_bias``/``master_dark``/``master_flat``), not by a basis change.
MATCH_BASIS_CALIB = "era_exact"

#: match_basis for a frame that carries NO target name but whose coordinates
#: fall inside a staged target's cone (see :func:`cone_match`).  Deliberately
#: distinct from MATCH_BASIS_SCIENCE: a coordinate match is a *candidate*
#: identification, not the project's published selection rule.
MATCH_BASIS_CONE = "cone_candidate"

#: Role of an adjudicated science frame.
ROLE_SCIENCE = "science"

#: Role of a cone-matched, name-less frame.  It is IN the working set (so a
#: project's Step-0 identity work happens inside the manifest instead of
#: reaching back into ``frames`` and bypassing S0c provenance) but it is NOT
#: ``science``: every downstream stage must opt in explicitly.
ROLE_SCIENCE_UNRESOLVED = "science_unresolved"

#: The two roles that carry a target identity.  Everything else in a stage
#: table is calibration.  Consumers must use these tuples (or the SQL
#: fragments below) rather than testing ``role = 'science'`` by hand — that
#: test silently reclassified cone candidates as calibration.
SCIENCE_ROLES: tuple[str, ...] = (ROLE_SCIENCE, ROLE_SCIENCE_UNRESOLVED)


def _role_sql(roles: tuple[str, ...], negate: bool = False) -> str:
    """SQL membership test over ``roles`` — generated, never hand-typed.

    Role names are bare identifiers-in-quotes with no apostrophes (asserted
    below), so the literal list is safe to interpolate into a WHERE clause.
    """
    assert all("'" not in r for r in roles), "role name with a quote"
    listed = ", ".join(f"'{r}'" for r in roles)
    return f"role {'NOT ' if negate else ''}IN ({listed})"


#: ``role IN ('science','science_unresolved')`` — every identity-bearing row.
SQL_SCIENCE_ROLES = _role_sql(SCIENCE_ROLES)

#: ``role NOT IN (...)`` — every calibration row (raw frames and masters).
SQL_CALIB_ROLES = _role_sql(SCIENCE_ROLES, negate=True)


# ---------------------------------------------------------------------------
# The five project selections — reviewable data, not code.
# ---------------------------------------------------------------------------

#: The 19 Dwarf-survey fields, verified against the manifest's canonical
#: target keys (DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md section 3:
#: "Memo 1 said 21 fields; verified 19").  Encoded explicitly — a prefix
#: rule could silently absorb a future mislabeled 'dw...' name.
DW_FIELDS: tuple[str, ...] = (
    "dw1403+49", "dw1409+51", "dw1418+46", "dw1441+51", "dw1446+58",
    "dw1459+44", "dw1533+67", "dw1539+45", "dw1558+67", "dw1559+46",
    "dw1608+40", "dw1615+54", "dw1617+46", "dw1633+69", "dw1643+07",
    "dw1645+46", "dw1709+74", "dw1721+71", "dw1735+57",
)

#: EVERY spelling of "this frame is dispersed light" that the archive uses.
#:
#: THE 'HaG' LESSON (2026-08-18 review).  The archive names its grisms five
#: different ways across three acquisition systems: ``hrg``/``lrg`` (pyscope
#: filenames), ``HaGrism``/``OGGrism`` (MaxIm wheel labels) and ``HaG`` (the
#: Andor/iKon tree's own label, 192 canonical rows).  A photometric selection
#: that enumerates only four of them leaks dispersed spectra into an aperture
#: photometry set: on 2024-04-16 the observer shot sequence-paired triples per
#: target — ``ST LMi-0001_hires`` (FILTER='HaG', 120 s), ``…_lowres``
#: (FILTER='OGGrism', 60 s), ``…_r`` (FILTER='r', 20 s) — and the four-name
#: blacklist excluded the ``_lowres`` twin while staging its ``_hires`` twin
#: as photometry.  The same stage table simultaneously carried
#: ``HaG.flat1..5`` as role='flat', i.e. the pipeline already knew HaG was a
#: filter needing its own flats.  ONE constant now feeds every rule below, so
#: a sixth spelling is a one-line change in one place.
GRISM_ALL: frozenset[str] = frozenset(
    {"hrg", "lrg", "HaGrism", "OGGrism", "HaG"})

#: Non-grism filter strings the CV strategy's canonical accounting rule also
#: excludes: a wheel slot with no bandpass ('6'), an open position ('empty'),
#: and the unresolved 'W' code (CV_TimeSeries/ANALYSIS_STRATEGY.md section 3:
#: "photometric filters only (exclude grisms, `empty`, `W`, `6`)").
CV_NON_PHOTOMETRIC_FILTERS: frozenset[str] = frozenset({"empty", "W", "6"})

#: The CV exclusion set, DERIVED so it can never fall behind GRISM_ALL again.
CV_EXCLUDED_FILTERS: frozenset[str] = GRISM_ALL | CV_NON_PHOTOMETRIC_FILTERS

#: The BeStar grism filter whitelist (BeStar_Grism/ANALYSIS_STRATEGY.md
#: Step 0: "explicit whitelist hrg/lrg/HaGrism/OGGrism; `lrgblue`
#: logged-and-excluded").  Same set as manifest.GRISM4_FILTERS.
#:
#: NOTE this is deliberately NARROWER than :data:`GRISM_ALL`: 'HaG' is a real
#: grism spelling, but it belongs to the Andor/iKon and ``grism/`` trees and
#: Step 0's published whitelist names four strings, not five.  Widening it
#: would silently change a published inventory; verified harmless either way
#: — no BeStar target has a single canonical 'HaG' frame (regression-tested
#: against the live manifest in test_staging.py).  The tests also assert this
#: set is a SUBSET of GRISM_ALL, so the two can only ever drift narrower.
BESTAR_GRISM_FILTERS: frozenset[str] = frozenset(
    {"hrg", "lrg", "HaGrism", "OGGrism"})

#: T CrB filters excluded from the science set.  TCrB_Monitoring/
#: ANALYSIS_STRATEGY.md section 3 rules on this code directly: "`H` (6,
#: presumed Halpha, all single-epoch 2024-03-13 — **excluded from science
#: regardless of P0-2 mapping**; filter table only)".  The filter-forensics
#: table the strategy still wants is built from S0's ``frames`` (which keeps
#: every frame), not from the working set — staging honours the published
#: science exclusion.
TCRB_EXCLUDED_FILTERS: frozenset[str] = frozenset({"H"})

#: Radius of the cone-candidate match for name-less frames, in degrees.
#:
#: Chosen, not inherited: S0's synonym gate uses 0.2° to decide whether two
#: NAMES are the same object, which is a different question.  Here we ask
#: whether a frame with no name at all was pointed AT a staged target, so the
#: scale is the pointing error of a frame that is otherwise on target.  0.25°
#: sits above the observed on-target scatter (the three matches in the live
#: archive land at 0.026°) and far below the nearest confusion case (the
#: 2025-01-23 focus frames sit 1.4° from λ Eri and are correctly refused).
#: Widening it to 2° would sweep in 228 Rigel focus frames — measured, and
#: the reason this number is small and stated here rather than guessed.
CONE_CANDIDATE_RADIUS_DEG: float = 0.25


@dataclass(frozen=True)
class ProjectSelection:
    """One project's science-frame selection rule, as reviewable data.

    Attributes
    ----------
    project
        The project directory name (``<repo>/<project>/data/`` receives the
        CSV) — matches the ``project`` column of S0's project_counts table.
    targets
        Normalized S0 target keys (``frames.target_key``) the project claims.
    filter_whitelist
        When set, a science frame is kept ONLY if its filter is in this set
        (BeStar: grism filters only).  ``None`` = no whitelist.
    filter_blacklist
        Filters explicitly excluded (CV: grisms/empty/W/6).  Applied after
        the whitelist; a blank filter never matches a blacklist entry.
    pending_alias_targets
        TRANSITIONAL keys only: target keys that the S0 alias fix already
        committed to :data:`macro_core.manifest.SYNONYM_TABLE` will fold into
        ``targets`` at the next S0 rebuild.  Listing them here keeps the
        working set correct BEFORE that rebuild; afterwards the keys no
        longer exist in ``frames`` and the entry is a harmless no-op.  A test
        asserts every entry resolves through SYNONYM_TABLE into ``targets``,
        so this list cannot drift away from the alias fix it mirrors.
    cone_radius_deg
        When set, name-less frames (blank ``target_key``) that pass every
        other gate and land within this many degrees of a staged target's
        reference position are emitted as ``science_unresolved`` rows with
        ``match_basis='cone_candidate'``.  ``None`` disables cone matching.
    rule
        One human sentence stating the selection — printed in the README,
        the report, and nowhere restated by hand.
    source
        Where the target list/rule was published (strategy doc section).
    """
    project: str
    targets: tuple[str, ...]
    rule: str
    source: str
    filter_whitelist: Optional[frozenset[str]] = None
    filter_blacklist: frozenset[str] = field(default_factory=frozenset)
    pending_alias_targets: tuple[str, ...] = ()
    cone_radius_deg: Optional[float] = None

    @property
    def all_target_keys(self) -> frozenset[str]:
        """Every key gate 5 accepts: the published list + pending aliases."""
        return frozenset(self.targets) | frozenset(self.pending_alias_targets)


#: The five staging selections.  Target lists come from
#: ``manifest.STRATEGY_CLAIMS`` plus each project's ANALYSIS_STRATEGY.md
#: (rev 2026-08-16) — the docstring of each entry's ``source`` says which.
PROJECT_SELECTIONS: tuple[ProjectSelection, ...] = (
    ProjectSelection(
        project="TCrB_Monitoring",
        targets=("tcrb", "tetcrb"),
        filter_blacklist=TCRB_EXCLUDED_FILTERS,
        rule=("Canonical error-free Light frames of T CrB and the θ CrB "
              "calibrator in every filter EXCEPT 'H' — the 2025 grism "
              "series, the 2023–2024 imaging anchors, and the calibrator "
              "series are one working set; the six single-epoch 2024-03-13 "
              "'H' frames are excluded from science by §3's explicit "
              "ruling (they remain visible in S0's frames table, which is "
              "where the filter-forensics table is built)."),
        source=("TCrB_Monitoring/ANALYSIS_STRATEGY.md §3 (T CrB 471 unique "
                "rawimage light frames — 402 after global dedup — + θ CrB "
                "412-frame grism calibrator series; 'H' excluded from "
                "science regardless of P0-2 mapping); STRATEGY_CLAIMS "
                "tcrb/tetcrb rows."),
    ),
    ProjectSelection(
        project="CV_TimeSeries",
        targets=("stlmi", "yzcnc", "vvpup", "euuma", "anuma"),
        filter_blacklist=CV_EXCLUDED_FILTERS,
        rule=("Canonical error-free Light frames of the five CVs in "
              "photometric filters only (ALL FIVE grism spellings — hrg, "
              "lrg, HaGrism, OGGrism, HaG — plus 'empty', 'W', '6' "
              "excluded) — the strategy's canonical accounting rule, "
              "widened from rawimage-only to all canonical trees so the "
              "iKon-tree VV Pup/YZ Cnc/ST LMi frames stage too; the "
              "widening imported the iKon tree's filter vocabulary, which "
              "is why the exclusion set is derived from GRISM_ALL."),
        source=("CV_TimeSeries/ANALYSIS_STRATEGY.md §3 canonical rule + "
                "§3.1 per-target table; STRATEGY_CLAIMS CV rows."),
    ),
    ProjectSelection(
        project="SN2023ixf_LightCurve",
        targets=("2023ixf", "m101"),
        pending_alias_targets=("pinwheelgalaxy", "2023ixf1", "2023ixf2"),
        rule=("Canonical error-free Light frames labeled 2023ixf, plus ALL "
              "canonical M101/NGC5457/'Pinwheel Galaxy' field frames — the "
              "saturated first epochs (2023-05-20/21), the pre-explosion "
              "template (2023-05-05) and every post-fade template epoch "
              "carry the host's name, not the SN's, and two of them carry a "
              "sequence digit fused to the SN's name ('2023ixf1/2')."),
        source=("SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md §3.1 (campaign "
                "start resolution) + §3.4 (templates); STRATEGY_CLAIMS "
                "2023ixf row.  ngc5457, 'pinwheel galaxy', 2023ixf1 and "
                "2023ixf2 are cone-gated S0 synonyms (SYNONYM_TABLE)."),
    ),
    ProjectSelection(
        project="BeStar_Grism",
        targets=("lameri", "5cnc", "69ori", "53boo", "hd70340",  # science
                 "phecda", "phileo", "spica",                    # refs/short tier
                 "hr3454", "hr4963", "vega"),                    # standards
        filter_whitelist=BESTAR_GRISM_FILTERS,
        cone_radius_deg=CONE_CANDIDATE_RADIUS_DEG,
        rule=("Canonical error-free Light frames of the core-ten targets "
              "plus Vega (Alpha Lyr synonym-merged, era-A ladder included) "
              "in the grism filter whitelist hrg/lrg/HaGrism/OGGrism; "
              "'lrgblue' and direct-imaging filters excluded by the "
              "whitelist.  T CrB is deliberately absent ('not this paper').  "
              "Name-less grism frames within "
              f"{CONE_CANDIDATE_RADIUS_DEG:g}° of a staged target also enter "
              "the working set, as role='science_unresolved' rows Step 0 "
              "must adjudicate before any of them counts as science."),
        source=("BeStar_Grism/ANALYSIS_STRATEGY.md §3.2 inventory table + "
                "Step 0 filter whitelist and its blank-target_best cone "
                "match; STRATEGY_CLAIMS BeStar rows."),
    ),
    ProjectSelection(
        project="DwarfGalaxy_AGN_Survey",
        targets=("ngc5548", "ngc5238") + DW_FIELDS,
        rule=("Canonical error-free Light frames of NGC 5548 (slot-'6' "
              "campaign + the 2024-06-05 revisit), NGC 5238, and the 19 "
              "verified Dw survey fields, all filters.  The mispointed "
              "2023-03-25 NGC 5548 night stays IN the manifest — its "
              "pointing flag, not its absence, is the record (S0 derives "
              "NGC 5548's reference position from header coordinates "
              "because none of its 279 frames is plate-solved, so the "
              "~8° offset flags as pointing_gt1deg like any other)."),
        source=("DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md §3 (19 fields "
                "verified; NGC 5548 143-frame call); STRATEGY_CLAIMS Dwarf "
                "rows (__dw_survey__ sentinel expanded to DW_FIELDS)."),
    ),
)


def selection_for(project: str) -> ProjectSelection:
    """Look up one project's selection by name (KeyError if unknown)."""
    for sel in PROJECT_SELECTIONS:
        if sel.project == project:
            return sel
    raise KeyError(f"no staging selection for project {project!r}")


# ---------------------------------------------------------------------------
# Pure predicates and small helpers
# ---------------------------------------------------------------------------

def stage_table_name(project: str) -> str:
    """SQLite table name for a project's staging manifest.

    ``stage_`` + the project name lowercased with every non-alphanumeric
    run collapsed to ``_`` — deterministic, collision-free across the five
    projects (asserted in the tests), and always a bare SQL identifier so
    it can be interpolated into DDL without quoting.
    """
    slug = re.sub(r"[^a-z0-9]+", "_", project.lower()).strip("_")
    return f"stage_{slug}"


def passes_frame_gates(sel: ProjectSelection,
                       imagetyp: Optional[str],
                       error: Optional[str],
                       is_canonical: object,
                       tree: Optional[str],
                       filter_name: Optional[str],
                       basename: str = "") -> bool:
    """The target-independent gates every staged frame must pass.

    Split out of :func:`is_staged_science` so the cone-candidate path
    (:func:`is_cone_candidate`) applies the SAME frame-quality and filter
    rules — a name-less frame that would not have been science had it been
    named must not sneak in through the coordinate door.

    1.  **canonical** — S0's global (basename, jd) dedup already chose the
        one true copy; everything else is a duplicate.
    2.  **not the reduced tree** — its canonical rows are renamed copies of
        rawimage frames (S0b linkage evidence); staging them would re-stage
        the same exposure twice under two names.
    3.  **error-free** — a frame the cataloger could not read has no
        trustworthy metadata to stage.
    4.  **science, not calibration** — mirror of S0b's science universe:
        :func:`macro_core.inventory.is_science` (Light frames plus the
        blank-IMAGETYP 2026 nights), with the calibration kind computed
        from the same header/basename rules.
    5.  **the project's filter rule** — whitelist first (if any), then
        blacklist.  A blank filter passes a blacklist (nothing to match)
        but fails a whitelist (not provably a listed filter).
    """
    # Gate 1 — canonical only (accept 1/True; None/0/NaN all fail).
    if not (is_canonical == 1 or is_canonical is True):
        return False
    # Gate 2 — the reduced tree never stages.
    if tree == "reduced":
        return False
    # Gate 3 — error-free.
    if error is not None and str(error).strip():
        return False
    # Gate 4 — science under the shared S0b definition.
    kind = inv.calib_kind(imagetyp, basename or "")
    if not inv.is_science(imagetyp, kind):
        return False
    # Gate 5 — filter whitelist / blacklist.
    if sel.filter_whitelist is not None:
        if filter_name not in sel.filter_whitelist:
            return False
    if filter_name is not None and filter_name in sel.filter_blacklist:
        return False
    return True


def is_staged_science(sel: ProjectSelection,
                      target_key: Optional[str],
                      imagetyp: Optional[str],
                      error: Optional[str],
                      is_canonical: object,
                      tree: Optional[str],
                      filter_name: Optional[str],
                      basename: str = "") -> bool:
    """Does one manifest frame belong to ``sel``'s science set?

    :func:`passes_frame_gates` plus the target gate: an exact match against
    :attr:`ProjectSelection.all_target_keys` (the published target list plus
    any transitional pending-alias keys).  Exact, never prefix: the alias
    merging already happened in S0, so one key IS one target, and LIKE-prefix
    matching is the documented way frames get silently dropped.
    """
    if not passes_frame_gates(sel, imagetyp, error, is_canonical, tree,
                              filter_name, basename):
        return False
    return target_key in sel.all_target_keys


def _blank(value: Optional[str]) -> bool:
    """True when a catalog string is missing or whitespace-only."""
    return value is None or not str(value).strip()


def is_cone_candidate(sel: ProjectSelection,
                      target_key: Optional[str],
                      imagetyp: Optional[str],
                      error: Optional[str],
                      is_canonical: object,
                      tree: Optional[str],
                      filter_name: Optional[str],
                      ra_deg: Optional[float],
                      dec_deg: Optional[float],
                      basename: str = "") -> bool:
    """Is this frame eligible for the cone match (before running it)?

    WHY THIS EXISTS.  BeStar_Grism/ANALYSIS_STRATEGY.md §3.2 records 671
    raw-tree grism-light rows with blank ``target_best`` and requires Step 0
    to resolve them by coordinate.  Under an exact-target-key gate those rows
    could never enter any stage table, so Step 0 would have had to query
    ``frames`` directly — the provenance bypass S0c exists to prevent.  The
    eligibility test is: same frame gates as science, NO target name at all,
    and usable coordinates.  Whether the frame then matches is
    :func:`cone_match`'s job.
    """
    if sel.cone_radius_deg is None:
        return False
    if not _blank(target_key):
        return False            # named frames go through gate 5, not here
    if ra_deg is None or dec_deg is None:
        return False            # a cone match needs a position
    if ra_deg != ra_deg or dec_deg != dec_deg:   # NaN != NaN
        return False
    return passes_frame_gates(sel, imagetyp, error, is_canonical, tree,
                              filter_name, basename)


def cone_match(ra_deg: float, dec_deg: float,
               refs: dict[str, tuple[float, float]],
               radius_deg: float) -> Optional[tuple[str, float]]:
    """Nearest reference position within ``radius_deg``, or ``None``.

    ``refs`` maps target key → (ra, dec) in degrees.  Returns
    ``(target_key, separation_deg)`` for the SINGLE nearest match — ties are
    impossible in practice and the nearest rule makes the outcome
    deterministic regardless of dict order.  Separation uses S0's own
    :func:`macro_core.manifest.angular_separation_deg`, so a cone here means
    exactly what a cone means in the alias gate and the pointing audit.
    """
    best: Optional[tuple[str, float]] = None
    for key in sorted(refs):                       # sorted = deterministic
        ra0, dec0 = refs[key]
        sep = mf.angular_separation_deg(ra_deg, dec_deg, ra0, dec0)
        if sep <= radius_deg and (best is None or sep < best[1]):
            best = (key, sep)
    return best


def role_of_calib(kind: str, is_master: object) -> str:
    """Staging role of a calibration frame: ``bias`` .. ``master_flat``.

    ``kind`` is the S0b-normalized kind (bias/dark/flat); masters — the
    recovered ``Calibrations/`` and ``calib/`` master products — get a
    ``master_`` prefix so a stage can prefer or exclude them with one
    string test.  An unknown kind raises: staging never guesses.
    """
    if kind not in ("bias", "dark", "flat"):
        raise ValueError(f"unknown calibration kind {kind!r}")
    return f"master_{kind}" if (is_master == 1 or is_master is True) else kind


def abs_archive_path(archive_root: str, rel_path: str) -> str:
    """Absolute path of an archive-relative manifest path.

    POSIX join, no normalization tricks: the catalog's relative paths are
    already clean, and the archive root may contain spaces (it does), which
    is why every consumer must quote this value in shell contexts.
    """
    return posixpath.join(archive_root, rel_path)


def farm_link_name(night: Optional[str], basename: str) -> str:
    """Symlink filename inside the optional browsable farm.

    Prefixed with the night label so that identically named frames from
    different nights (a real archive pattern) cannot collide inside one
    role directory.  Frames with no night label group under ``no-night``.
    """
    return f"{night or 'no-night'}_{basename}"


def science_row(sel: ProjectSelection, frame: dict, archive_root: str,
                build_id: str) -> dict:
    """Build one science staging row (STAGE_CSV_COLUMNS order) from a
    manifest ``frames`` record (column-name → value mapping)."""
    return {
        "path": frame["path"],
        "abs_path": abs_archive_path(archive_root, frame["path"]),
        "role": "science",
        "match_basis": MATCH_BASIS_SCIENCE,
        "tree": frame.get("tree"),
        "era_id": frame.get("era_id"),
        "night": frame.get("night"),
        "jd": frame.get("jd"),
        "filter": frame.get("filter"),
        "exptime": frame.get("exptime"),
        "canonical_target": frame.get("canonical_target"),
        "target_key": frame.get("target_key"),
        "dup_group": frame.get("dup_group"),
        "qc_flags": frame.get("qc_flags"),
        "pointing_offset_deg": frame.get("pointing_offset_deg"),
        "size_bytes": frame.get("size"),
        "obs_rowid": frame.get("obs_rowid"),
        "stage_build_id": build_id,
    }


def cone_candidate_row(frame: dict, matched_key: str, sep_deg: float,
                       archive_root: str, build_id: str) -> dict:
    """Build one ``science_unresolved`` row for a cone-matched frame.

    Same schema as every other row, with three deliberate differences:

    * ``role`` is ``science_unresolved`` — never ``science``;
    * ``match_basis`` is ``cone_candidate`` — never ``selection_rule``;
    * ``target_key`` carries the CANDIDATE key and
      ``pointing_offset_deg`` carries the measured separation from that
      target's reference position, so Step 0 can rank and adjudicate the
      candidates without recomputing anything.

    ``canonical_target`` stays NULL: the frame has no name in the archive,
    and inventing one is precisely the silent promotion this row prevents.
    """
    return {
        "path": frame["path"],
        "abs_path": abs_archive_path(archive_root, frame["path"]),
        "role": ROLE_SCIENCE_UNRESOLVED,
        "match_basis": MATCH_BASIS_CONE,
        "tree": frame.get("tree"),
        "era_id": frame.get("era_id"),
        "night": frame.get("night"),
        "jd": frame.get("jd"),
        "filter": frame.get("filter"),
        "exptime": frame.get("exptime"),
        "canonical_target": None,
        "target_key": matched_key,
        "dup_group": frame.get("dup_group"),
        "qc_flags": frame.get("qc_flags"),
        "pointing_offset_deg": sep_deg,
        "size_bytes": frame.get("size"),
        "obs_rowid": frame.get("obs_rowid"),
        "stage_build_id": build_id,
    }


def calib_row(calib: dict, archive_root: str, build_id: str) -> dict:
    """Build one calibration staging row from an S0b ``calib_frames``
    record (plus its ``size`` joined from ``frames``).

    Calibration frames have no target identity or pointing — those columns
    stay NULL so the one shared schema serves both row types.
    """
    return {
        "path": calib["path"],
        "abs_path": abs_archive_path(archive_root, calib["path"]),
        "role": role_of_calib(calib["kind"], calib.get("is_master")),
        "match_basis": MATCH_BASIS_CALIB,
        "tree": calib.get("tree"),
        "era_id": calib.get("era_id"),
        "night": calib.get("night"),
        "jd": calib.get("jd"),
        "filter": calib.get("filter"),
        "exptime": calib.get("exptime"),
        "canonical_target": None,
        "target_key": None,
        "dup_group": None,
        "qc_flags": None,
        "pointing_offset_deg": None,
        "size_bytes": calib.get("size"),
        "obs_rowid": calib.get("obs_rowid"),
        "stage_build_id": build_id,
    }


# ---------------------------------------------------------------------------
# Published-inventory reconciliation: strategy docs vs the working set
# ---------------------------------------------------------------------------
#
# WHY THIS TABLE EXISTS.  Three separate numbers in three strategy documents
# were tree-doubled or stale at the 2026-08-18 review — "30 frames/band" for
# the SN template stacks (canonical reality 15), "20 frames" for the NGC 5548
# revisit (10), "403 frames / 40 nights" for the theta CrB series (412 / 42,
# because two more nights arrived after the doc was written).  Every one of
# them is quoted in a manuscript.  Staging was right in all three cases; the
# prose drifted, and nothing made the drift visible.
#
# So each doc-published inventory number is encoded here NEXT TO the query
# that reproduces it from the working set, and the S0c report renders the
# diff.  A doc that drifts (or an archive that grows) now shows up as a
# highlighted row on the evidence page instead of at referee time.

@dataclass(frozen=True)
class StageClaim:
    """One inventory number a strategy document publishes.

    Attributes
    ----------
    project
        Which project's stage table the claim is measured against.
    label
        What the document says it is counting, in its own words.
    claimed_frames, claimed_nights
        The published numbers.  ``claimed_nights=None`` when the document
        states no night count for this row.
    where
        A SQL WHERE fragment over that project's stage table.  The renderer
        ALWAYS conjoins ``role = 'science'`` itself, so a claim can never
        accidentally count calibration or cone-candidate rows.
    source
        Document and section the claim comes from.
    """
    project: str
    label: str
    claimed_frames: int
    claimed_nights: Optional[int]
    where: str
    source: str


def assert_safe_where(where: str) -> str:
    """Reject a claim fragment that is not a bare read-only expression.

    The fragments are repo-authored data, not user input, but they ARE
    interpolated into SQL, so the constraint is enforced rather than trusted:
    no statement separator, no comment marker, and no attached-statement
    keyword.  Returns the fragment so it can be used inline.
    """
    lowered = where.lower()
    assert ";" not in where, f"claim fragment with a statement separator: {where}"
    assert "--" not in where, f"claim fragment with a comment marker: {where}"
    for word in ("drop", "delete", "insert", "update", "attach", "pragma"):
        assert word not in lowered.split(), \
            f"claim fragment with the keyword {word!r}: {where}"
    return where


#: Every inventory number the five strategy documents publish that the stage
#: tables can reproduce.  Values are the CORRECTED ones (the 2026-08-18 doc
#: edits); a future drift in either direction lights up the report.
STAGE_CLAIMS: tuple[StageClaim, ...] = (
    # --- CV: the §3.1 per-target table, under §3's canonical rule ---------
    # NOTE the rule excludes non-photometric filters, so these are the
    # photometric counts (3,150 not 3,157 etc.) — the doc's own rule applied
    # to the doc's own table, which is where three of its five rows drifted.
    StageClaim("CV_TimeSeries", "ST LMi — rawimage photometric light frames",
               3150, 39, "tree = 'rawimage' AND target_key = 'stlmi'",
               "CV_TimeSeries/ANALYSIS_STRATEGY.md §3.1"),
    StageClaim("CV_TimeSeries", "YZ Cnc — rawimage photometric light frames",
               1915, 26, "tree = 'rawimage' AND target_key = 'yzcnc'",
               "CV_TimeSeries/ANALYSIS_STRATEGY.md §3.1"),
    StageClaim("CV_TimeSeries", "VV Pup — rawimage photometric light frames",
               1277, 28, "tree = 'rawimage' AND target_key = 'vvpup'",
               "CV_TimeSeries/ANALYSIS_STRATEGY.md §3.1"),
    StageClaim("CV_TimeSeries", "EU UMa — rawimage photometric light frames",
               993, 32, "tree = 'rawimage' AND target_key = 'euuma'",
               "CV_TimeSeries/ANALYSIS_STRATEGY.md §3.1"),
    StageClaim("CV_TimeSeries", "AN UMa — rawimage photometric light frames",
               1279, 14, "tree = 'rawimage' AND target_key = 'anuma'",
               "CV_TimeSeries/ANALYSIS_STRATEGY.md §3.1"),
    # --- T CrB: the two grism series --------------------------------------
    StageClaim("TCrB_Monitoring", "T CrB grism series (the paper's spine)",
               247, 60,
               "target_key = 'tcrb' AND \"filter\" IN ('hrg', 'lrg')",
               "TCrB_Monitoring/ANALYSIS_STRATEGY.md §1/§3"),
    StageClaim("TCrB_Monitoring", "θ CrB grism calibrator series",
               412, 42,
               "target_key = 'tetcrb' AND \"filter\" IN ('hrg', 'lrg')",
               "TCrB_Monitoring/ANALYSIS_STRATEGY.md §3"),
    # --- SN: the primary broadband template epoch -------------------------
    StageClaim("SN2023ixf_LightCurve", "2024-05-19 deep g template stack",
               15, 1, "night = '2024-05-18' AND \"filter\" = 'g'",
               "SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md §3.4"),
    StageClaim("SN2023ixf_LightCurve", "2024-05-19 deep r template stack",
               15, 1, "night = '2024-05-18' AND \"filter\" = 'r'",
               "SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md §3.4"),
    StageClaim("SN2023ixf_LightCurve",
               "2026-03-21/22 M101-field epoch ('pinwheel galaxy')",
               140, 2, "night IN ('2026-03-21', '2026-03-22')",
               "SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md §3.4 "
               "(added 2026-08-18: the deep pre-2026-04 template epoch)"),
    # --- Dwarf: the NGC 5548 campaign and its revisit ----------------------
    StageClaim("DwarfGalaxy_AGN_Survey", "NGC 5548 slot-'6' campaign (2023)",
               143, 16,
               "target_key = 'ngc5548' AND night < '2024-01-01'",
               "DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md §3"),
    StageClaim("DwarfGalaxy_AGN_Survey", "NGC 5548 2024-06-05 grism revisit",
               10, 1,
               "target_key = 'ngc5548' AND night > '2024-01-01'",
               "DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md §3"),
    # --- BeStar: the §3.2 per-target inventory table ----------------------
    StageClaim("BeStar_Grism", "Spica", 728, 29, "target_key = 'spica'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "HR 3454 (η Hya)", 471, 34,
               "target_key = 'hr3454'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "HR 4963 (θ Vir)", 368, 33,
               "target_key = 'hr4963'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "Phecda", 333, 40, "target_key = 'phecda'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "53 Boo", 330, 23, "target_key = '53boo'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "69 Ori", 288, 51, "target_key = '69ori'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "φ Leo", 261, 39, "target_key = 'phileo'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "λ Eri", 252, 47, "target_key = 'lameri'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "5 Cnc", 194, 45, "target_key = '5cnc'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "HD 70340", 181, 30,
               "target_key = 'hd70340'",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    # Vega splits into three instrument tiers, and the doc's table stops at
    # era C — the era-83 tier arrived with a later ingest.
    StageClaim("BeStar_Grism", "Vega — era-A exposure ladder (2024-05-20)",
               80, 1, "target_key = 'vega' AND era_id = 72",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "Vega — era-C standard series", 367, 18,
               "target_key = 'vega' AND era_id IN (78, 80)",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2"),
    StageClaim("BeStar_Grism", "Vega — era-83 tier (post-doc ingest)",
               30, 3, "target_key = 'vega' AND era_id = 83",
               "BeStar_Grism/ANALYSIS_STRATEGY.md §3.2 "
               "(added 2026-08-18: nights the doc's table pre-dates)"),
)


def claims_for(project: str) -> tuple[StageClaim, ...]:
    """Every published inventory claim measured against ``project``."""
    return tuple(c for c in STAGE_CLAIMS if c.project == project)
