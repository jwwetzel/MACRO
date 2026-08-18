"""macro_grism — the MACRO G-track slitless-spectroscopy core.

Modules (each header explains its evidence base; the design decisions all
trace to the T CrB grism strategy and the exploration documented in
``docs/pipeline/g_grism.html``):

* ``fits_io``     — HDU-resolution FITS reader (plain / fpack / repackaged).
* ``trace``       — trace geometry: slope fit, detilted cross-dispersion
                    profile, main-trace location, per-column centers.
* ``gate``        — the identity gate: Gaia DR3 field verification of a
                    slitless frame with no WCS.
* ``extract``     — Horne-style optimal extraction with flanking-band
                    background, plus the master-dark comparison arm.
* ``wavelength``  — per-frame self-anchored wavelength solution (Halpha +
                    telluric O2), per-grism dispersion aggregation.
* ``db``          — the g_* manifest tables (new tables only; existing
                    tables are never modified).
* ``report_g``    — the chain-of-evidence HTML report (reads ONLY the DB).
"""
