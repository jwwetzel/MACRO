<!-- emitted by BeStar_Grism/scripts/be_calib.py; do not edit -->
### Calibration from the archive, per mechanical state (BE-S2)

| mech_state | bias | dark | flat | grism_flats | std_nights | dark_subtraction | flat_field | error_floor_label |
|---|---|---|---|---|---|---|---|---|
| S1 Andor | 1 | 1 | 48 | 13 | 0.00 | master-dark arm available + flanking band | none (no hrg/lrg flat; BE-X1 dither test deferred) | LOWER-BOUND-ONLY (intra-night) |
| S2 ASI pre-monsoon | 102 | 733 | 56 | 0 | 0.00 | master-dark arm available + flanking band | none (no hrg/lrg flat; BE-X1 dither test deferred) | LOWER-BOUND-ONLY (intra-night) |
| S3 ASI post-monsoon (flipped) | 4 | 116 | 24 | 0 | 28.00 | master-dark arm available + flanking band | none (no hrg/lrg flat; BE-X1 dither test deferred) | standards night-to-night (BE-S10) |
| S4 QHY night 1 | 0 | 0 | 0 | 0 | 0.00 | flanking band only | none (no hrg/lrg flat; BE-X1 dither test deferred) | LOWER-BOUND-ONLY (intra-night) |
| S5 QHY re-seated | 0 | 0 | 0 | 0 | 34.00 | flanking band only | none (no hrg/lrg flat; BE-X1 dither test deferred) | standards night-to-night (BE-S10) |

#### Per F-3 epoch (archive `calib_frames`, any filter)

| mech_state | mech_epoch | bias | dark | flat | grism_flats | gain_config | gain_e_per_adu | gain_err | gain_status |
|---|---|---|---|---|---|---|---|---|---|
| S1 Andor | iKon:2024-05-12 | 0 | 0 | 23 | 13 | iKon 1MHz 4x | 0.940 | 0.007 | measured |
| S1 Andor | iKon:2024-10-03 | 1 | 1 | 8 | 0 | iKon 5MHz | — | — | UNMEASURED |
| S1 Andor | iKon:2024-10-25 | 0 | 0 | 0 | 0 | iKon 5MHz | — | — | UNMEASURED |
| S1 Andor | iKon:2024-11-01 | 0 | 0 | 17 | 0 | iKon 1MHz 4x | 0.940 | 0.007 | measured |
| S2 ASI pre-monsoon | ASI:2024-12-16 | 102 | 733 | 56 | 0 | ASI Mode0 2x2 | 1.046 | 0.015 | measured |
| S3 ASI post-monsoon (flipped) | ASI:2025-10-11 | 4 | 116 | 24 | 0 | ASI Mode0 2x2 | 1.046 | 0.015 | measured |
| S4 QHY night 1 | QHY600:2026-03-21 | 0 | 0 | 0 | 0 | QHY600 Fast 2x2 | 1.288 | 0.060 | measured |
| S5 QHY re-seated | QHY600:2026-03-22 | 0 | 0 | 0 | 0 | QHY600 Fast 2x2 | 1.288 | 0.060 | measured |
| S5 QHY re-seated | QHY600:2026-04-01 | 0 | 0 | 0 | 0 | QHY600 Fast 2x2 | 1.288 | 0.060 | measured |
| S5 QHY re-seated | QHY600:2026-04-14 | 0 | 0 | 0 | 0 | QHY600 Fast 2x2 | 1.288 | 0.060 | measured |
| S5 QHY re-seated | QHY600:2026-04-23 | 0 | 0 | 0 | 0 | QHY600 Fast 2x2 | 1.288 | 0.060 | measured |
| S5 QHY re-seated | QHY600:2026-06-28 | 0 | 0 | 0 | 0 | QHY600 Fast 2x2 | 1.288 | 0.060 | measured |

Method (stated in the paper): flanking-band background everywhere, master-dark arm where the state has darks; no flat-field division (no hrg/lrg flat exists for any camera; the Andor's HaGrism/OGGrism flats are not applied, for one method across states — the dither test that would bound the pixel-response residual is BE-X1, 2027 backlog); gain/read noise/saturation from the measured `s2_camera_configs`/`s2_linearity_caps`. The QHY600 now mounted has no bias, dark or flat in the archive; its acquisition is BE-X3-calib-acquisition (2027 backlog). Era B (ASI) before the standards epoch carries LOWER-BOUND-ONLY errors; ASI after the monsoon and the QHY take their floor from the standards (BE-S10).
