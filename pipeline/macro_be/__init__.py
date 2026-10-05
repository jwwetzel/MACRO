"""macro_be — the Be-star grism paper's code version (provenance stage ``BE``).

The Be-star pipeline itself lives beside its strategy, in
``BeStar_Grism/scripts/`` (be_step0 -> be_calib -> be_external -> be_cadence ->
be_extract -> be_measure -> be_series -> be_injection -> be_figures ->
be_numbers), and writes ``BeStar_Grism/products/bestar.sqlite`` and
``BeStar_Grism/products/be_grism.sqlite``.  This package exists only so the
provenance graph can read one code-version string; bump it whenever a change
would alter a published number.
"""
BE_CODE_VERSION = "BE v1.0 (2026-10-05)"
