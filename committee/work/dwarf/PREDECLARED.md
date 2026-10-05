# DwarfGalaxy_AGN_Survey — pre-declared rules (written 2026-10-05T02:42Z, before any stack or Hα measurement exists)

Per Seat 6 (committee/reviews/2026-10-05-seat6-dwarf-abstract.md, rev. 2 decision, edit 3).

- **In-band sample** (\NinBand): objects with Hα frames AND a published velocity putting Hα inside the
  assumed 65 Å band at 6562.8 Å (novelty_table.csv band = "in band"): Dw1459+44, Dw1533+67, Dw1558+67,
  Dw1559+46. Three have published emission spectra (Dw1459+44 SDSS; Dw1558+67, Dw1559+46 BTA).
- **Detection** (dwcore.HA_DETECT_SIGMA): net/σ ≥ 5 in the 10″ aperture, σ from random apertures on the same
  continuum-subtracted stack; 3 ≤ net/σ < 5 = marginal, reported as a limit. Limit = max(net,0) + 3σ.
- **Informative limit**: the 3σ limit is below the predicted flux by at least the DW-N1 margin, i.e.
  F_pred ≥ 3 × F_lim (F_pred = FUV-based where GALEX detects the object, else constant-SFR P=0, exactly as
  in novelty_table.csv; margin 3 because the FUV predictor over-predicts by up to +0.33 dex on the
  validation set).
- **Venue trigger**: AJ (short paper, Figs 2–5) if ≥ 2 in-band objects are detected or have informative
  limits; otherwise RNAAS with one figure (Fig 4). Evaluated once, by script, from the DW-P36 table.
- **NGC 5238** is the flux calibrator: its Hα result is morphological only; no NGC 5238 flux is quoted as a result.

## Addendum (2026-10-05T03:37Z, still before any Hα aperture is measured)
- **Continuum scale.** H and R stacks are both tied to REFCAT2 r at g−r = 0.6, so the continuum scale is
  k = 10^(−0.4 (c_H − c_R)(gr_obj − 0.6)) with c the fitted colour terms and gr_obj inferred from the object's
  own L−R colour in the same aperture (gr_obj = 0.6 + (m_L − m_R)/(c_L − c_R)); pivot colour where L−R is
  not measured at S/N ≥ 10. The continuum-scale uncertainty (σ_gr propagated) enters σ_total in quadrature
  with the random-aperture σ.
- **Detection** requires net/σ_total ≥ 5 (σ_total ⊇ σ_rand). Limits use σ_total.
- **NGC 5238 aperture**: the radius where the net Hα curve of growth rises by < 2 % over the next 10″,
  evaluated on the all-night stack; the same aperture on every per-night stack gives the scatter.

## DW-P4x rule (2026-10-05T04:14Z, before any NGC 5548 spectrum is extracted)
- No grism-library solution exists for the AC4040 slot '6' (G-1 epochs are ASI/QHY/ANDOR only), so the
  triage self-calibrates each frame's wavelength scale on NGC 5548's own lines: [O III] λ5007 and Hα at
  z = 0.01717 (two points → linear dispersion per frame).
- EW(Hα) = the whole Hα+[N II] complex (unresolved at this LSF) above a straight continuum through
  observed windows 6420–6480 Å and 6900–6980 Å (the second avoids the O₂ B band), rest frame.
- Night: the slot-'6' night with the most frames passing extraction (ties → earliest).
- **Earns a paragraph** iff the nightly-mean EW uncertainty (per-frame scatter/√N) is < 2 %, a third of the
  2023 Hβ variability amplitude (F_var = 5.9 %, Xi et al. 2025); otherwise DW-P4x is dropped.

## Evaluation (2026-10-05T05:27Z, by run_dw_paper.py halpha; dw_build_meta venue_trigger_n)
- In-band objects detected or with informative limits: 1 (Dw1459+44 detected, net/σ = 16.1). Dw1558+67's
  prediction is 2.95× its limit (threshold 3): not informative. Dw1533+67 1.11×, Dw1559+46 0.93×.
- Trigger (≥ 2 → AJ) returns **RNAAS with Fig 4 only**. Figs 2, 3, 5 are not built.
