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

#: Filters the CV strategy's canonical accounting rule EXCLUDES
#: (CV_TimeSeries/ANALYSIS_STRATEGY.md section 3: "photometric filters only
#: (exclude grisms, `empty`, `W`, `6`)").
CV_EXCLUDED_FILTERS: frozenset[str] = frozenset(
    {"hrg", "lrg", "HaGrism", "OGGrism", "empty", "W", "6"})

#: The BeStar grism filter whitelist (BeStar_Grism/ANALYSIS_STRATEGY.md
#: Step 0: "explicit whitelist hrg/lrg/HaGrism/OGGrism; `lrgblue`
#: logged-and-excluded").  Same set as manifest.GRISM4_FILTERS; restated
#: here as the selection's own datum so a change to either is a visible diff.
BESTAR_GRISM_FILTERS: frozenset[str] = frozenset(
    {"hrg", "lrg", "HaGrism", "OGGrism"})


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


#: The five staging selections.  Target lists come from
#: ``manifest.STRATEGY_CLAIMS`` plus each project's ANALYSIS_STRATEGY.md
#: (rev 2026-08-16) — the docstring of each entry's ``source`` says which.
PROJECT_SELECTIONS: tuple[ProjectSelection, ...] = (
    ProjectSelection(
        project="TCrB_Monitoring",
        targets=("tcrb", "tetcrb"),
        rule=("Canonical error-free Light frames of T CrB and the θ CrB "
              "calibrator, all filters — the 2025 grism series, the "
              "2023–2024 imaging anchors, and the calibrator series are one "
              "working set."),
        source=("TCrB_Monitoring/ANALYSIS_STRATEGY.md §3 (T CrB 471 unique "
                "rawimage light frames — 402 after global dedup — + θ CrB "
                "403-frame grism calibrator series); STRATEGY_CLAIMS "
                "tcrb/tetcrb rows."),
    ),
    ProjectSelection(
        project="CV_TimeSeries",
        targets=("stlmi", "yzcnc", "vvpup", "euuma", "anuma"),
        filter_blacklist=CV_EXCLUDED_FILTERS,
        rule=("Canonical error-free Light frames of the five CVs in "
              "photometric filters only (grisms, 'empty', 'W', '6' "
              "excluded) — the strategy's canonical accounting rule, "
              "widened from rawimage-only to all canonical trees so the "
              "iKon-tree VV Pup/YZ Cnc/ST LMi frames stage too."),
        source=("CV_TimeSeries/ANALYSIS_STRATEGY.md §3 canonical rule + "
                "§3.1 per-target table; STRATEGY_CLAIMS CV rows."),
    ),
    ProjectSelection(
        project="SN2023ixf_LightCurve",
        targets=("2023ixf", "m101"),
        rule=("Canonical error-free Light frames labeled 2023ixf, plus ALL "
              "canonical M101/NGC5457 field frames — the saturated first "
              "epochs (2023-05-20/21), the pre-explosion template "
              "(2023-05-05) and every post-fade template epoch carry the "
              "host's name, not the SN's."),
        source=("SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md §3.1 (campaign "
                "start resolution) + §3.4 (templates); STRATEGY_CLAIMS "
                "2023ixf row.  ngc5457 is a cone-gated S0 synonym of m101."),
    ),
    ProjectSelection(
        project="BeStar_Grism",
        targets=("lameri", "5cnc", "69ori", "53boo", "hd70340",  # science
                 "phecda", "phileo", "spica",                    # refs/short tier
                 "hr3454", "hr4963", "vega"),                    # standards
        filter_whitelist=BESTAR_GRISM_FILTERS,
        rule=("Canonical error-free Light frames of the core-ten targets "
              "plus Vega (Alpha Lyr synonym-merged, era-A ladder included) "
              "in the grism filter whitelist hrg/lrg/HaGrism/OGGrism; "
              "'lrgblue' and direct-imaging filters excluded by the "
              "whitelist.  T CrB is deliberately absent ('not this paper')."),
        source=("BeStar_Grism/ANALYSIS_STRATEGY.md §3.2 inventory table + "
                "Step 0 filter whitelist; STRATEGY_CLAIMS BeStar rows."),
    ),
    ProjectSelection(
        project="DwarfGalaxy_AGN_Survey",
        targets=("ngc5548", "ngc5238") + DW_FIELDS,
        rule=("Canonical error-free Light frames of NGC 5548 (slot-'6' "
              "campaign + the 2024-06-05 revisit), NGC 5238, and the 19 "
              "verified Dw survey fields, all filters.  The mispointed "
              "2023-03-25 NGC 5548 night stays IN the manifest — its "
              "pointing flag, not its absence, is the record."),
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


def is_staged_science(sel: ProjectSelection,
                      target_key: Optional[str],
                      imagetyp: Optional[str],
                      error: Optional[str],
                      is_canonical: object,
                      tree: Optional[str],
                      filter_name: Optional[str],
                      basename: str = "") -> bool:
    """Does one manifest frame belong to ``sel``'s science set?

    The gates, in order (all must pass):

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
    5.  **the project's target list** (exact key match — the alias merging
        already happened in S0, so one key IS the target).
    6.  **the project's filter rule** — whitelist first (if any), then
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
    # Gate 5 — the project's target list.
    if target_key not in sel.targets:
        return False
    # Gate 6 — filter whitelist / blacklist.
    if sel.filter_whitelist is not None:
        if filter_name not in sel.filter_whitelist:
            return False
    if filter_name is not None and filter_name in sel.filter_blacklist:
        return False
    return True


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
