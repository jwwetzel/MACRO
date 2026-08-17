# MACRO

Analysis pipeline and projects for the MACRO consortium.

The [MACRO Consortium](https://macroconsortium.org/) (Augustana College, Coe College,
Knox College, Macalester College, University of Iowa) operates the 0.5 m Robert L.
Mutel Telescope (RLMT) at Winer Observatory, Sonoita AZ. This repository holds the
consortium's shared analysis pipeline, the five archival research projects built on
the 2023–2026 RLMT archive (3.34 TiB, ~329k FITS), and the transparent
chain-of-evidence reports that document every analytical decision.

## Principles

1. **Every number is script-emitted.** No hand-transcribed counts; each paper ships a
   SQL/code appendix that regenerates its tables.
2. **Every decision is documented.** Each pipeline stage and project maintains an HTML
   chain-of-evidence report (`docs/`): question → evidence plot → decision → consequence.
3. **Self-contained, readable code.** Each script runs standalone, states its inputs
   and outputs, and is commented line-by-line so any consortium member can audit it.
4. **No untested statistics.** Detection claims ride on injection–recovery; error bars
   are χ²-validated; periods survive alias tests or they don't ship.

## Layout

| Path | Contents |
|---|---|
| `ROADMAP.md` | Cross-project sequencing, shared-pipeline design, risk register |
| `pipeline/` | `macro_core` and sibling packages (manifest, astrometry, timing, photometry, time series, grism) |
| `docs/` | GitHub Pages site: chain-of-evidence reports for the pipeline and each project |
| `TCrB_Monitoring/` | Recurrent nova T CrB pre-eruption baseline (strategy + scripts + figures) |
| `CV_TimeSeries/` | Polars & dwarf novae time-series photometry |
| `SN2023ixf_LightCurve/` | SN 2023ixf early multi-band campaign |
| `BeStar_Grism/` | Be-star Hα slitless-grism monitoring |
| `DwarfGalaxy_AGN_Survey/` | Dwarf-candidate deep imaging + NGC 5548 monitoring |

Each project directory carries its committee-reviewed `ANALYSIS_STRATEGY.md`
(post-internal-referee) and grows `scripts/`, `figures/`, and `notes/` as analysis
proceeds. Manuscripts are deliberately not tracked here.

## Data

The RLMT archive and the observation catalog (`rlmt-catalog.sqlite`, one row per FITS
with header metadata) live outside the repo on local storage. The pipeline treats the
catalog as read-only ground truth for inventory and reads pixels only from the archive's
`rawimage/` tree. Products land under `products/` (untracked, regenerable).

## Environment

```bash
conda env create -f environment.yml   # env name: rlmt-checks
```
