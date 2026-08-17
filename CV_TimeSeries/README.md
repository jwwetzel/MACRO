# Cataclysmic-Variable Time Series

Archival light curves of the RLMT's deeply covered polars and dwarf novae: orbital light
curves, eclipse timings, and outburst statistics.

| Target | Frames | Nights | Span |
|--------|-------:|-------:|------|
| ST LMi | 3,102 | 37 | 2024-01 → 2026-02 |
| VV Pup | 1,253 | 27 | 2024-11 → 2025-11 |
| YZ Cnc | 1,462+ | 22+ | 2024-01 → 2024-03 |
| EU UMa | 919 | 30 | 2025-01 → 2026-06 |
| AG LMi | 968 | 36 | 2023-03 → 2026-01 |
| QQ Gem | 105 | 24 | 2024-12 → 2026-04 |

## Data
- Archive: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/`
- Catalog: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite` (table `obs`)
- Python: conda env `rlmt-checks` (astropy 6.1.7)

## Starter query
```sql
SELECT target_best, path, date_obs, jd, filter, exptime, fwhm, zmag
FROM obs
WHERE target_best IN ('ST LMi','VV Pup','YZ Cnc','EU UMa','AG LMi','QQ Gem')
  AND imagetyp LIKE 'Light%' AND error IS NULL
ORDER BY target_best, jd;
```

## Notes
- Photometric pipeline is shared with `../TCrB_Monitoring` — build once, reuse here.
- High-cadence intra-night sequences (e.g. ST LMi) are suited to eclipse/orbital work;
  sparse multi-year coverage (AG LMi) to outburst statistics.
