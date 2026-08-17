# Be-Star Grism Campaign

Hα emission variability in bright Be stars from the RLMT's slitless grism program:
**~23k low-res (`lrg`) and high-res (`hrg`) grism frames** on a deliberate multi-month
cadence, 2024 → 2026. Core targets: **Phecda (γ UMa), φ Leo, 53 Boo, 69 Ori, λ Eri,
Spica, HR 3454, HR 4963, 5 Cnc, HD 70340**.

## Data
- Archive: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/` (also see `grism/`, `hagrism/` trees)
- Catalog: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite` (table `obs`)
- Python: conda env `rlmt-checks` (astropy 6.1.7)

## Starter query
```sql
SELECT target_best, path, date_obs, filter, exptime
FROM obs
WHERE (filter IN ('hrg','lrg','HaGrism','OGGrism') OR fn_filter IN ('hrg','lrg'))
  AND imagetyp LIKE 'Light%' AND error IS NULL
ORDER BY target_best, date_obs;
```

## Notes
- Filename parsing quirk: decimal exposures (e.g. `0-25s`) leak into `fn_target`, so
  Phecda appears as several variants (`PHECDA`, `PHECDA lrg 0-25s`, …) — match with
  `fn_target LIKE 'PHECDA%'` etc.
- Extraction pipeline: trace + collapse spectra, wavelength-calibrate on known lines,
  then track Hα equivalent width vs. time. Thesis-grade coherent dataset.
