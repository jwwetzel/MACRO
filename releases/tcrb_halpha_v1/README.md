# T CrB — RLMT slitless-grism Halpha release v1

Generated 2026-10-05T06:42:46Z by `pipeline/scripts/run_tcrb.py release` (git b89290c) from `products/tcrb/tcrb.sqlite` and the shared grism library. Every number is a query.

- `spectra/` — 181 reduced 1-D spectra (FITS binary tables: wavelength on the fixed G-1 dispersion with a per-frame Halpha zero point; optimal flux, its variance; boxcar flux). Index: `tcrb_rlmt_spectra_index.csv`.
- `aras_format/` — the same 181 spectra on a linear grid, normalised at 6470-6520 A, with BSS keywords, for submission to the ARAS T CrB database.
- `tcrb_rlmt_nightly_ew.csv` — 83 nightly EWs (per grism; error = statistical (+) floor).
- `tcrb_rlmt_line_flux.csv` — line flux = EW x AAVSO Rc continuum, orbital phase (Munari et al. 2025 ephemeris).
- `aras_cross_validation.csv`, `b_anchors_2024.csv`, `filter_codes.csv`, `flickering_limits.csv`.

EW convention (pre-registered, ANALYSIS_STRATEGY.md §10): emission positive; Halpha 6562.8 A +- max(30 A, 1.5 LSF FWHM); pseudo-continuum through the medians of 6470-6520 and 6600-6640 A.

Archive DOI: pending (the Zenodo deposit and the ARAS submission are owner actions).
