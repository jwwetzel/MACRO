# SN 2023ixf Early Light Curve

Densely sampled multi-band early light curve of **SN 2023ixf in M101**: **3,018 frames over
33 nights, 2023-05-22 → 2023-07-07**, beginning within days of discovery (2023-05-19).

## Data
- Archive: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/`
- Catalog: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite` (table `obs`)
- Python: conda env `rlmt-checks` (astropy 6.1.7)

## Starter query
```sql
SELECT path, date_obs, jd, filter, exptime, airmass, fwhm, zmag
FROM obs
WHERE (target_best LIKE '%2023ixf%' OR target_best LIKE '%M101%')
  AND imagetyp LIKE 'Light%' AND error IS NULL
ORDER BY jd;
```

## Notes
- Galaxy-background subtraction matters here — template imaging of M101 may exist in the
  archive from before discovery (search `target_best LIKE '%M101%' AND jd < 2460084`).
- Compare against published 2023ixf photometry to gauge calibration quality; potential
  angle: combining with the community's dense coverage or student-led reproduction study.
