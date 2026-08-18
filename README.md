# MACRO

Analysis pipeline and projects for the MACRO consortium.

The [MACRO Consortium](https://macroconsortium.org/) (Augustana College, Coe College,
Knox College, Macalester College, University of Iowa) operates the 0.5 m Robert L.
Mutel Telescope (RLMT) at Winer Observatory, Sonoita AZ. This repository holds the
consortium's shared analysis pipeline, the five archival research projects built on
the 2023–2026 RLMT archive (3.34 TiB, ~329k FITS), and the transparent
chain-of-evidence reports that document every analytical decision.

## Evidence site

Every analytical decision is published as a chain-of-evidence report at
**[jwwetzel.github.io/MACRO](https://jwwetzel.github.io/MACRO/)** — question, plot, decision,
consequence, with every number emitted by the pipeline rather than typed by hand.

| Project | What it is | Evidence |
|---|---|---|
| T CrB pre-eruption monitoring | Hα series through the 2025 dip recovery + photometric anchors | [report](docs/TCrB_Monitoring/index.html) |
| Cataclysmic-variable time series | cycle-resolved light curves of 3 polars + YZ Cnc superhumps | [report](docs/CV_TimeSeries/index.html) |
| SN 2023ixf early light curve | gri/narrowband +5.4→+50 d; possible flash-phase grism series | [report](docs/SN2023ixf_LightCurve/index.html) |
| Be-star grism campaign | ~3-day-cadence Hα equivalent-width monitoring | [report](docs/BeStar_Grism/index.html) |
| Dwarf galaxies + NGC 5548 | deep Hα imaging + band-integrated AGN light curve | [report](docs/DwarfGalaxy_AGN_Survey/index.html) |
| Rigel-era legacy archive *(candidate)* | pre-MACRO 2015–2023 contact-binary program, 210k frames | [report](docs/Legacy_Rigel/index.html) |

Shared pipeline evidence (manifest, calibration inventory, astrometry, detector, timing,
photometry, staging) is indexed on the [hub page](docs/index.html).

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
