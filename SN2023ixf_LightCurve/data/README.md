# SN2023ixf_LightCurve — staging manifest (`stage_manifest.csv`)

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

**Selection rule (science rows).** Canonical error-free Light frames labeled 2023ixf, plus ALL canonical M101/NGC5457/'Pinwheel Galaxy' field frames — the saturated first epochs (2023-05-20/21), the pre-explosion template (2023-05-05) and every post-fade template epoch carry the host's name, not the SN's, and two of them carry a sequence digit fused to the SN's name ('2023ixf1/2').
Source: SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md §3.1 (campaign start resolution) + §3.4 (templates); STRATEGY_CLAIMS 2023ixf row.  ngc5457, 'pinwheel galaxy', 2023ixf1 and 2023ixf2 are cone-gated S0 synonyms (SYNONYM_TABLE).

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

**This build (S0c v1.1 (2026-10-03) @ 2026-10-03T23:54Z):** 1,349 science rows +
0 cone-candidate rows + 472 calibration rows.

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
