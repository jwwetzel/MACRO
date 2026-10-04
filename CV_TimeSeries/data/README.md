# CV_TimeSeries — staging manifest (`stage_manifest.csv`)

**THE NO-COPY LAW.** No frame is ever copied into this directory. S0 exists
because copies proliferated (132k duplicate rows in the archive catalog);
this manifest **is** the working set. Every pipeline stage reads the
immutable archive directly through the paths below — the archive is
read-only, always.

**This file is regenerable, not precious** (`*/data/` is gitignored). Run
this from anywhere — the path is absolute and quoted because the repo path
contains spaces:

    /opt/miniconda3/envs/rlmt-checks/bin/python \
        "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO/pipeline/scripts/build_s0c_staging.py"

**Selection rule (science rows).** Canonical error-free Light frames of the five CVs in photometric filters only (ALL FIVE grism spellings — hrg, lrg, HaGrism, OGGrism, HaG — plus 'empty', 'W', '6' excluded) — the strategy's canonical accounting rule, widened from rawimage-only to all canonical trees so the iKon-tree VV Pup/YZ Cnc/ST LMi frames stage too; the widening imported the iKon tree's filter vocabulary, which is why the exclusion set is derived from GRISM_ALL.
Source: CV_TimeSeries/ANALYSIS_STRATEGY.md §3 canonical rule + §3.1 per-target table; STRATEGY_CLAIMS CV rows.

**Calibration rows.** For every camera era the science frames touch, that
era's calibration frames from the S0b census are included (raw frames and
recovered `Calibrations/` masters alike) **wherever they may be applied to
this science**, `match_basis = 'era_epoch_exact'`. Staging
deliberately over-includes on kind/exposure/filter; each stage narrows those
with the S0b coverage matrix as its guide.

**The boundary rule.** An era is a header history; the hardware history is
the `mech_epoch` column. A flat is staged only if some science frame of this
project shares its `mech_epoch` (the camera was not re-seated, rotated,
flipped or re-wheeled in between); a dark or bias only if one shares its
`detector_epoch` (same camera, same flip state). **Match per frame**: a
calibration row may be applied to a science row only when their epoch
columns agree — `macro_core.inventory.calib_valid_for` is the rule.

**This build (S0c v1.1 (2026-10-03) @ 2026-10-03T23:54Z):** 8,641 science rows +
0 cone-candidate rows + 4,665 calibration rows.

## Columns

| column | meaning |
|---|---|
| `path` | archive-relative POSIX path — the frame's identity |
| `abs_path` | absolute archive path (QUOTE IT: the root has spaces) |
| `role` | `science`, `science_unresolved` (cone candidate — NOT science until a project adjudicates it), `bias`/`dark`/`flat`, or `master_*` products |
| `match_basis` | `selection_rule` (science: the rule below), `cone_candidate` (no target name; matched by coordinates) or `era_epoch_exact` (calibration: same S0 era as this project's science) |
| `tree` | top-level archive tree holding the canonical copy |
| `era_id` | S0 pinned camera-era registry id |
| `night` | local-noon-to-noon night label |
| `jd` | header JD = **UTC exposure START** (BJD_TDB is stage S3's job — never use this for timing) |
| `filter` | cataloged filter string |
| `exptime` | header EXPTIME (s) |
| `canonical_target` | S0 alias-merged display name (science rows) |
| `target_key` | S0 normalized target key (science rows) |
| `dup_group` | S0 global duplicate-group id |
| `qc_flags` | S0 QC flags — flags mark, they never delete |
| `pointing_offset_deg` | offset from the target's reference position |
| `size_bytes` | integrity surrogate (see note below) |
| `obs_rowid` | catalog/manifest join key |
| `stage_build_id` | S0c build that emitted the row |
| `mech_epoch` | S0b mechanical epoch `<camera>:<first night>` — a FLAT may be applied only to science in the same one |
| `detector_epoch` | S0b detector epoch — a DARK or BIAS may be applied only to science in the same one |
| `epoch_certain` | 1 = the night is placed in its epoch with certainty; 0 = it lies in the gap before a rotation-only boundary (no flat is valid) |

**Integrity note.** size_bytes is an integrity SURROGATE, not a checksum: it comes from the S0 catalog scan and catches truncation/replacement at read time.  A content hash would require re-reading the full 3.3 TiB archive — that is a separate archive-custody decision, not part of a staging build.

## Optional symlink farm (`frames/`)

`build_s0c_staging.py --symlink-farm` materializes
`frames/<role>/<night>_<basename>` symlinks into the archive for humans who
want a browsable view. The farm is **disposable** (delete and regenerate at
will) and **Dropbox does not sync symlink targets** — it is a local
convenience, never a transport mechanism.
