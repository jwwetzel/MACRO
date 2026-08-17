# T CrB Pre-Eruption Monitoring

Quiescent photometric baseline of the recurrent nova **T CrB** ahead of its anticipated
eruption. RLMT coverage: **414 frames over 83 nights, 2023-05-12 → 2025-06-24** (check for
newer nights as the archive syncs). Headers carry per-image photometric zero-points (`ZMAG`)
and seeing (`FWHM`), enabling a calibrated multi-band light curve of the pre-eruption state —
scientifically valuable the moment the nova goes off.

## Data
- Archive: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/` (raw: `rawimage/<night>/`, pipeline: `reduced/`)
- Catalog: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite` (table `obs`, one row per FITS)
- Python: conda env `rlmt-checks` (astropy 6.1.7)

## Starter query
```sql
SELECT path, date_obs, filter, exptime, airmass, fwhm, zmag
FROM obs
WHERE target_best = 'T CrB' AND imagetyp LIKE 'Light%' AND error IS NULL
ORDER BY date_obs;
```

## Notes
- Cross-check against AAVSO photometry for validation.
- `target_best` merges the OBJECT header with pyscope filename parsing; also try
  `LIKE '%CrB%'` to catch naming variants.
