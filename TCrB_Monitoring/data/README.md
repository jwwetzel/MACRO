# TCrB_Monitoring — staging manifest (`stage_manifest.csv`)

**THE NO-COPY LAW.** No frame is ever copied into this directory. S0 exists
because copies proliferated (132k duplicate rows in the archive catalog);
this manifest **is** the working set. Every pipeline stage reads the
immutable archive directly through the paths below — the archive is
read-only, always.

**This file is regenerable, not precious** (`*/data/` is gitignored):

    /opt/miniconda3/envs/rlmt-checks/bin/python \
        pipeline/scripts/build_s0c_staging.py

**Selection rule (science rows).** Canonical error-free Light frames of T CrB and the θ CrB calibrator, all filters — the 2025 grism series, the 2023–2024 imaging anchors, and the calibrator series are one working set.
Source: TCrB_Monitoring/ANALYSIS_STRATEGY.md §3 (T CrB 471 unique rawimage light frames — 402 after global dedup — + θ CrB 403-frame grism calibrator series); STRATEGY_CLAIMS tcrb/tetcrb rows.

**Calibration rows.** For every camera era the science frames touch, ALL of
that era's calibration frames from the S0b census are included (raw frames
and recovered `Calibrations/` masters alike), `match_basis =
'era_exact'`. Staging deliberately over-includes; each stage
narrows by kind/exposure/filter with the S0b coverage matrix as its guide.

**This build (S0c v1.0 (2026-08-18) @ 2026-08-18T15:25Z):** 826 science rows +
1,472 calibration rows.

## Columns

| column | meaning |
|---|---|
| `path` | archive-relative POSIX path — the frame's identity |
| `abs_path` | absolute archive path (QUOTE IT: the root has spaces) |
| `role` | `science`, `bias`/`dark`/`flat`, or `master_*` products |
| `match_basis` | `selection_rule` (science: the rule below) or `era_exact` (calibration: same S0 era as this project's science) |
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

**Integrity note.** size_bytes is an integrity SURROGATE, not a checksum: it comes from the S0 catalog scan and catches truncation/replacement at read time.  A content hash would require re-reading the full 3.3 TiB archive — that is a separate archive-custody decision, not part of a staging build.

## Optional symlink farm (`frames/`)

`build_s0c_staging.py --symlink-farm` materializes
`frames/<role>/<night>_<basename>` symlinks into the archive for humans who
want a browsable view. The farm is **disposable** (delete and regenerate at
will) and **Dropbox does not sync symlink targets** — it is a local
convenience, never a transport mechanism.
