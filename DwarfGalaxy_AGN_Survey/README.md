# Dwarf-Galaxy Survey + AGN Monitoring

Two related deep-imaging programs from the 2023 season:

1. **Dwarf-galaxy candidate fields** (`Dw####+##`, e.g. Dw1403+49, Dw1409+51, Dw1418+46,
   Dw1533+67): ~13–16 nights each, Jun–Jul 2023 — low-surface-brightness candidate
   confirmation via deep stacking.
2. **NGC 5548 AGN monitoring**: 116 frames over 13 nights (2023-03 → 2023-04) —
   continuum reverberation-style sampling; NGC 5238 also has 21 nights.

## Data
- Archive: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/`
- Catalog: `/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite` (table `obs`)
- Python: conda env `rlmt-checks` (astropy 6.1.7)

## Starter query
```sql
SELECT target_best, COUNT(*) frames, COUNT(DISTINCT substr(date_obs,1,10)) nights
FROM obs
WHERE (target_best LIKE 'Dw1%' OR target_best LIKE 'NGC 55%' OR target_best LIKE 'NGC55%')
  AND imagetyp LIKE 'Light%' AND error IS NULL
GROUP BY target_best ORDER BY nights DESC;
```

## Notes
- Dwarf fields need careful flat-fielding/fringe correction for LSB work — see the
  archive's `fringes/` tree and the consortium's `defringe` pipeline on GitHub.
- Check whether the Dw fields tie to a published survey (Mutel/collaborators) before
  scoping the science goal.
