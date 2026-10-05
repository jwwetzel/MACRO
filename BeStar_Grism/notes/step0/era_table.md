<!-- emitted by BeStar_Grism/scripts/be_step0.py; do not edit -->
### Instrument table — the five mechanical states and their F-3 epochs

| mech_state | first_night | last_night | be_frames | gain_configs |
|---|---|---|---|---|
| S1 Andor | 2024-05-12 | 2024-11-07 | 89 | iKon 1MHz 4x,iKon 5MHz |
| S2 ASI pre-monsoon | 2024-12-16 | 2025-06-23 | 1479 | ASI Mode0 2x2 |
| S3 ASI post-monsoon (flipped) | 2025-10-11 | 2026-03-16 | 2018 | ASI Mode0 2x2,ASI Mode0 2x2 e0.780 |
| S4 QHY night 1 | 2026-03-21 | 2026-03-21 | 0 | QHY600 Fast 2x2 |
| S5 QHY re-seated | 2026-03-22 | 2026-07-02 | 1236 | QHY600 Fast 2x2 |

| mech_state | mech_epoch | boundary_cause | first_night | last_night | first_night_be | last_night_be | n_nights | n_frames | grisms | gain_config | gain_e_per_adu | gain_err | gain_status | read_noise_e | clip_adu |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S1 Andor | iKon:2024-05-12 | rotation_fine | 2024-05-12 | 2024-05-27 | 2024-05-19 | 2024-05-19 | 1 | 80 | HaGrism,OGGrism | iKon 1MHz 4x | 0.940 | 0.007 | measured | 6.945 | 64674.000 |
| S1 Andor | iKon:2024-10-03 | rotation_coarse | 2024-10-03 | 2024-10-24 | 2024-10-11 | 2024-10-11 | 1 | 3 | HaGrism | iKon 5MHz | — | — | UNMEASURED | — | — |
| S1 Andor | iKon:2024-10-25 | rotation_fine | 2024-10-25 | 2024-10-31 | 2024-10-31 | 2024-10-31 | 1 | 3 | HaGrism | iKon 5MHz | — | — | UNMEASURED | — | — |
| S1 Andor | iKon:2024-11-01 | rotation_fine | 2024-11-01 | 2024-11-07 | 2024-11-01 | 2024-11-01 | 1 | 3 | HaGrism | iKon 1MHz 4x | 0.940 | 0.007 | measured | 6.945 | 64674.000 |
| S2 ASI pre-monsoon | ASI:2024-12-16 | wheel_map | 2024-12-16 | 2025-06-23 | 2024-12-16 | 2025-06-13 | 71 | 1479 | HaGrism,OGGrism,hrg,lrg | ASI Mode0 2x2 | 1.046 | 0.015 | measured | 2.393 | 65535.000 |
| S3 ASI post-monsoon (flipped) | ASI:2025-10-11 | flipstat | 2025-10-11 | 2026-03-16 | 2025-10-28 | 2026-03-16 | 82 | 2008 | hrg,lrg | ASI Mode0 2x2 | 1.046 | 0.015 | measured | 2.393 | 65535.000 |
| S3 ASI post-monsoon (flipped) | ASI:2025-10-11 | flipstat | 2025-10-11 | 2026-03-16 | 2025-11-30 | 2025-11-30 | 1 | 10 | hrg | ASI Mode0 2x2 e0.780 | — | — | UNMEASURED | — | 65535.000 |
| S4 QHY night 1 | QHY600:2026-03-21 | camera | 2026-03-21 | 2026-03-21 | — | — | 0 | 0 |  | QHY600 Fast 2x2 | 1.288 | 0.060 | measured | — | 65534.000 |
| S5 QHY re-seated | QHY600:2026-03-22 | flipstat,wheel_map | 2026-03-22 | 2026-03-26 | 2026-03-23 | 2026-03-26 | 4 | 129 | hrg,lrg | QHY600 Fast 2x2 | 1.288 | 0.060 | measured | — | 65534.000 |
| S5 QHY re-seated | QHY600:2026-04-01 | wheel_map | 2026-04-01 | 2026-04-11 | 2026-04-01 | 2026-04-11 | 4 | 117 | hrg,lrg | QHY600 Fast 2x2 | 1.288 | 0.060 | measured | — | 65534.000 |
| S5 QHY re-seated | QHY600:2026-04-14 | wheel_map | 2026-04-14 | 2026-04-21 | 2026-04-15 | 2026-04-17 | 2 | 206 | hrg,lrg | QHY600 Fast 2x2 | 1.288 | 0.060 | measured | — | 65534.000 |
| S5 QHY re-seated | QHY600:2026-04-23 | wheel_map | 2026-04-23 | 2026-06-27 | 2026-04-24 | 2026-06-24 | 25 | 779 | hrg,lrg | QHY600 Fast 2x2 | 1.288 | 0.060 | measured | — | 65534.000 |
| S5 QHY re-seated | QHY600:2026-06-28 | wheel_map | 2026-06-28 | 2026-07-02 | 2026-06-28 | 2026-06-28 | 1 | 5 | hrg | QHY600 Fast 2x2 | 1.288 | 0.060 | measured | — | 65534.000 |

No-crossing test: frames without a mechanical epoch = 0; (grism, epoch) calibration keys spanning two states = 0; epochs whose Be frames fall outside the epoch bounds = 0. Gains are `s2_camera_configs` measurements (F-4); rows marked UNMEASURED carry no gain and are not used for a noise model.

Built 2026-10-05T05:32:36+00:00 at commit b89290c from manifest S0 v1.1 (2026-10-03) 2026-10-03T22:11:55.291316+00:00 and S2c S2c v1.2 (2026-10-04) (read-only).
