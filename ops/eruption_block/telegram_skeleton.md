# T CrB eruption — first-report skeleton (ATel / RNAAS) — DRAFT

> **Status: DRAFT skeleton. Every `[…]` is a placeholder to be filled by a script from the
> pipeline databases at the time; nothing here is a measurement.** Written 2026-10-03 in
> response to committee finding ED `TCRB-N2-eruption-contingency`. The compiled LaTeX
> version, with macros emitted from the equivalent-width database, belongs with the T CrB
> manuscript and is not built here (see `committee/work/ops/REPORT.md`).

**Title:** Slitless Hα spectroscopy of T CrB [N_HOURS] hours after eruption, with a
[N_BASELINE]-epoch pre-eruption baseline from the same instrument

**Authors:** [per the consortium eruption author list — `AUTHORSHIP.md` §6.4]

**Text (ATel: ≤ 4,000 characters):**

We report slitless grism spectroscopy of the recurrent nova T CrB obtained with the 0.5 m
Robert L. Mutel Telescope (Winer Observatory, Arizona) beginning [UT_FIRST], [N_HOURS] h
after the eruption reported by [DISCOVERY_REF]. Spectra were taken through two grisms
(dispersions [DISP_LRG] and [DISP_HRG] Å/pixel; resolving power at Hα [R_LRG] and
[R_HRG], measured from [LSF_SOURCE]) at [N_EPOCHS] epochs between [UT_FIRST] and
[UT_LAST].

The Hα emission equivalent width was [EW_FIRST] ± [EW_FIRST_ERR] Å at [UT_FIRST] and
[EW_LAST] ± [EW_LAST_ERR] Å at [UT_LAST]; the line FWHM was [FWHM_KMS] ± [FWHM_ERR]
km/s [state whether resolved: instrumental FWHM is [LSF_KMS] km/s].
[One sentence on other lines only if the wavelength solution of that night passes its
acceptance test.]

For comparison, the same instrument recorded T CrB on [N_BASELINE_NIGHTS] nights between
[BASELINE_START] and [BASELINE_END], during the recovery from the 2023–24 pre-eruption
dip. Over that interval the Hα equivalent width ranged from [EW_BASE_MIN] to
[EW_BASE_MAX] Å (median [EW_BASE_MED] Å, night-to-night floor [EW_FLOOR] Å from the
comparison star θ CrB). [If the 2026 QHY600 frames are used: the 2026 measurements were
made with a different camera; the offset between the two, from standards observed with
both, is [SPLICE_OFFSET] ± [SPLICE_ERR] Å.]

Equivalent widths are measured in fixed windows after degrading to a common resolution;
uncertainties are empirical. The continuum itself brightened by [DELTA_V] mag, so the
equivalent width and the line flux are reported separately: F(Hα) = [FLUX] ±
[FLUX_ERR] erg cm⁻² s⁻¹ using [PHOTOMETRY_SOURCE] for the continuum. Reduced spectra and
the measurement table are available at [DOI_OR_URL]. Observations continue at [CADENCE].

**Rules for filling this in**

1. No number is typed. Each placeholder maps to a macro emitted from the grism database.
2. An equivalent width is never quoted without the continuum statement (standing rule 6).
3. Line widths are quoted only with the measured instrumental width beside them.
4. If the 2026 camera's splice offset has not been measured, say so and quote 2026 and
   2025 values separately.
5. The author list is fixed in advance under `AUTHORSHIP.md` §6.4, not on the night.
