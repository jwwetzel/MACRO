<!-- emitted by BeStar_Grism/scripts/be_step0.py; do not edit -->
### Header re-scrape audit (F-2 cards, re-opened where NULL)

| camera | card | column | n_frames | n_null_manifest | n_manifest_null_but_card_present | n_filled_from_header | n_null_card_present_after |
|---|---|---|---|---|---|---|---|
| ASI | CCD-TEMP | ccd_temp | 3497 | 0 | 0 | 0 | 0 |
| ASI | COOLPOWR | coolpowr | 3497 | 0 | 0 | 0 | 0 |
| ASI | FLIPSTAT | flipstat | 3497 | 0 | 0 | 0 | 0 |
| ASI | FOCPOS | focpos | 3497 | 0 | 0 | 0 | 0 |
| ASI | FWALLNAM | fwallnam | 3497 | 0 | 0 | 0 | 0 |
| ASI | FWPOS | fwpos | 3497 | 0 | 0 | 0 | 0 |
| ASI | GAIN | hdr_gain | 3497 | 0 | 0 | 0 | 0 |
| ASI | OFFSET | hdr_offset | 3497 | 0 | 0 | 0 | 0 |
| ASI | SET-TEMP | set_temp | 3497 | 0 | 0 | 0 | 0 |
| ASI | TELPIER | telpier | 3497 | 0 | 0 | 0 | 0 |
| QHY600 | CCD-TEMP | ccd_temp | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | COOLPOWR | coolpowr | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | FLIPSTAT | flipstat | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | FOCPOS | focpos | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | FWALLNAM | fwallnam | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | FWPOS | fwpos | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | GAIN | hdr_gain | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | OFFSET | hdr_offset | 1236 | 5 | 5 | 5 | 0 |
| QHY600 | SET-TEMP | set_temp | 1236 | 0 | 0 | 0 | 0 |
| QHY600 | TELPIER | telpier | 1236 | 0 | 0 | 0 | 0 |
| iKon | CCD-TEMP | ccd_temp | 89 | 0 | 0 | 0 | 0 |
| iKon | COOLPOWR | coolpowr | 89 | 0 | 0 | 0 | 0 |
| iKon | FLIPSTAT | flipstat | 89 | 0 | 0 | 0 | 0 |
| iKon | FOCPOS | focpos | 89 | 0 | 0 | 0 | 0 |
| iKon | FWALLNAM | fwallnam | 89 | 0 | 0 | 0 | 0 |
| iKon | FWPOS | fwpos | 89 | 0 | 0 | 0 | 0 |
| iKon | GAIN | hdr_gain | 89 | 0 | 0 | 0 | 0 |
| iKon | OFFSET | hdr_offset | 89 | 89 | 0 | 0 | 0 |
| iKon | SET-TEMP | set_temp | 89 | 0 | 0 | 0 | 0 |
| iKon | TELPIER | telpier | 89 | 0 | 0 | 0 | 0 |

Acceptance: `n_null_card_present_after` = 0 everywhere — a NULL is left in `be_frames` only where the card is absent from every HDU of the file. Rows with `n_manifest_null_but_card_present` > 0 are manifest (F-2) defects, filled here from the header and reported to the S0 owner.

Built 2026-10-05T05:32:36+00:00 at commit b89290c from manifest S0 v1.1 (2026-10-03) 2026-10-03T22:11:55.291316+00:00 and S2c S2c v1.2 (2026-10-04) (read-only).
