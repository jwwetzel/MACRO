"""The grism track's manifest tables (NEW tables only — S0/S0b/S1 tables
are never modified; the manifest stays append-only across stages).

* ``g_gate_calib``   — the era's camera-orientation calibration: one row
                       per solved imaging frame, plus the adopted CD.
* ``g_extractions``  — one row per (frame, background method): the gate
                       verdict, trace statistics, wavelength anchors, the
                       Halpha snippet, and the flanking-vs-dark debt.
* ``g_build_meta``   — code version, timestamps, sample definition.

Writes go through one transaction per batch; a crashed run leaves the
last complete batch, and the runner resumes by skipping rows that exist
(the (obs_rowid, method) uniqueness is the resume key).

2026-10-03 (committee wave 0).  The calibration-library tables below
(``G2_SCHEMA``) are written to the grism track's OWN database,
``products/grism/grism.sqlite``, because the shared manifest is read-only
for every package but ``foundation-s0`` this wave; the integrated rebuild
lands them in the manifest unchanged.  Every table is keyed on the archive
``path`` (stable across an S0 rebuild) and carries ``obs_rowid`` only as a
convenience copy of the manifest key at the time of the run.

* ``g_frames``          — one row per frame this track has reduced: the
                          header cards the committee asked for (true
                          focus, sensor temperature, flip), mechanical
                          epoch, detector key, trace geometry and the
                          1-D quality numbers.
* ``g_line_id``         — D1: line identification under the adopted and
                          the rival dispersion hypotheses, per calibrator.
* ``g_line_meas``       — every measured line / band-edge centre.
* ``g_dispersion``      — G-1: THE solution, one row per (grism, epoch).
* ``g_dispersion_resid``— its per-line residuals (the acceptance test).
* ``g_telluric_lambda`` — measured effective wavelengths of the telluric
                          band edges on the stellar scale.
* ``g_zero_point``      — per-frame zero point and the Halpha-O2B
                          separation check for science frames.
* ``g_variance_check``  — G-2: predicted vs measured sky variance.
* ``g_identity``        — G-3: the pixel-based identity gate.
* ``g_null_sky``        — G-4: empty-aperture residual per sky method.
* ``g_lsf``             — G-5: measured line-spread function.
* ``g_ew``              — step 7: equivalent widths, old and new library.
* ``g_meta``            — key/value provenance.
"""

from __future__ import annotations

import sqlite3

G_SCHEMA = """
CREATE TABLE IF NOT EXISTS g_gate_calib (
    calib_id     INTEGER PRIMARY KEY,
    era_id       INTEGER,
    frame_path   TEXT,           -- the imaging frame that was solved
    night        TEXT,
    status       TEXT,           -- solved / unsolved / bad_solve
    cd1_1 REAL, cd1_2 REAL, cd2_1 REAL, cd2_2 REAL,   -- deg/px
    pixscale_arcsec REAL,
    rotation_deg REAL,
    rms_arcsec   REAL,
    n_matched    INTEGER,
    adopted      INTEGER DEFAULT 0    -- 1 = the CD the gate uses
);
CREATE TABLE IF NOT EXISTS g_extractions (
    obs_rowid    INTEGER NOT NULL,   -- frames.obs_rowid
    method       TEXT NOT NULL,      -- 'flanking' | 'masterdark'
    path         TEXT,
    target       TEXT,
    filter       TEXT,               -- hrg | lrg
    night        TEXT,
    jd           REAL,               -- header JD, UTC exposure START (S3
                                     -- owns BJD; nothing here converts)
    exptime      REAL,
    era_id       INTEGER,
    role         TEXT,               -- 'tcrb_sample' | 'gate_bad' |
                                     -- 'gate_good' | 'calibrator'
    layout       TEXT,               -- FITS packaging that was resolved
    -- identity gate ------------------------------------------------------
    gate_verdict TEXT,               -- ACCEPT | REJECT
    gate_reason  TEXT,
    pointing_offset_deg REAL,
    u_obs REAL, u_pred REAL, u_resid_px REAL,
    gate_parity  TEXT,
    n_gaia       INTEGER,
    brightest_g  REAL,
    -- trace ---------------------------------------------------------------
    trace_height REAL,
    trace_slope  REAL,
    trace_c0 REAL, trace_c1 REAL, trace_c2 REAL,   -- centers poly (deg 2)
    trace_rms_px REAL,
    trace_n_centroids INTEGER,
    -- extraction ----------------------------------------------------------
    bg_method    TEXT,               -- 'flanking' | 'masterdark+flanking'
    dark_path    TEXT,               -- master used (masterdark rows only)
    dark_exptime REAL,
    n_extracted  INTEGER,
    n_sat_cols   INTEGER,            -- columns with any saturated pixel
    peak_flux    REAL,               -- max optimal flux (ADU)
    median_flux  REAL,
    -- wavelength ----------------------------------------------------------
    anchor_status TEXT,              -- 'halpha+o2' | 'halpha_only' | 'none'
    x_halpha REAL, halpha_snr REAL, halpha_width_px INTEGER,
    x_o2b REAL, o2b_snr REAL,
    x_o2a REAL, o2a_snr REAL,
    disp_a_per_px REAL,              -- per-frame measurement (signed)
    disp_source  TEXT,
    -- products ------------------------------------------------------------
    snippet_json TEXT,               -- [x, flux] pairs around Halpha
    spectrum_fits TEXT,              -- products-relative output path
    contamination_flag TEXT,         -- 'C4_Be_Halpha' on tet CrB rows
    debt_median_rel_diff REAL,       -- flanking vs masterdark (both rows
                                     -- of a frame carry the same number)
    status       TEXT,               -- 'ok' | error text
    UNIQUE (obs_rowid, method)
);
CREATE TABLE IF NOT EXISTS g_build_meta (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""


G2_SCHEMA = """
CREATE TABLE IF NOT EXISTS g_frames (
    path         TEXT PRIMARY KEY,   -- archive-relative (stable key)
    obs_rowid    INTEGER,            -- manifest key at run time
    sample       TEXT,               -- calibrator | tcrb | tcrb_slot |
                                     -- ngc5548 | identity_control | ...
    target       TEXT,               -- manifest target_best
    star         TEXT,               -- calibrator registry name (if any)
    filter       TEXT,               -- FILTER card as written
    grism        TEXT,               -- hrg | lrg | NULL (slot 6 / W)
    night        TEXT,
    mech_epoch   TEXT,               -- macro_grism.config working table
    jd           REAL,               -- header JD, UTC exposure START
    exptime      REAL,
    detector_key TEXT,
    instrume TEXT, readoutm TEXT, flipstat TEXT, telpier TEXT,
    ccd_temp REAL, set_temp REAL,
    focpos REAL,                     -- pyscope FOCPOS (true focuser)
    focuspos REAL,                   -- MaxIm FOCUSPOS (can be stuck)
    airmass REAL,
    hdr_ra TEXT, hdr_dec TEXT,       -- OBJCTRA / OBJCTDEC as written
    pointing_offset_deg REAL,        -- manifest value (header vs target)
    ny INTEGER, nx INTEGER, layout TEXT,
    -- trace -----------------------------------------------------------
    status       TEXT,               -- ok | no_trace | error text
    trace_height REAL, trace_slope REAL, trace_u REAL,
    trace_rms_px REAL, trace_n INTEGER,
    trace_x0 INTEGER, trace_x1 INTEGER,     -- signal extent (columns)
    trace_coeffs TEXT,               -- JSON polynomial, polyval order
    -- extraction quality --------------------------------------------------
    peak_adu REAL,                   -- brightest raw pixel in aperture
    n_sat_cols INTEGER,              -- native-saturation-possible columns
    snr_median REAL,
    fwhm_px REAL,                    -- cross-dispersion FWHM (median)
    sky_adu REAL,                    -- sky above pedestal under trace
    spec_file TEXT,                  -- products-relative .npz
    code_version TEXT,
    hot_mask TEXT,                   -- badpix product used, or 'none'
    sat_cap_adu REAL,                -- saturation cap applied (stored ADU)
    pedestal_adu REAL,               -- pedestal used (table or frame)
    offset REAL                      -- OFFSET card
);
CREATE TABLE IF NOT EXISTS g_line_id (
    path TEXT PRIMARY KEY,
    star TEXT, grism TEXT, mech_epoch TEXT, spt TEXT,
    n_features INTEGER,
    -- adopted-seed hypothesis ------------------------------------------
    n_match INTEGER, n_stellar INTEGER, n_lines INTEGER,
    disp_at_ha REAL,                 -- signed A/px from the ID refinement
    sign INTEGER, x_halpha REAL, rms_px REAL,
    coeffs_json TEXT,                -- polyval order, lambda - 6562.80
    refined INTEGER,                 -- 1 = >=4 lines survived tight re-match
    matches_json TEXT,               -- [[line, x, resid_px], ...]
    -- rival hypotheses (same features, same tolerance) ------------------
    rivals_json TEXT,                -- {name: {n_match, n_stellar, rms}}
    best_rival_match INTEGER,
    hbeta_depth_159 REAL,            -- hrg: depth where 1.59 A/px puts
    hbeta_depth_159_err REAL         --   Hbeta (best of both signs)
);
CREATE TABLE IF NOT EXISTS g_line_meas (
    path TEXT NOT NULL, line TEXT NOT NULL,
    star TEXT, grism TEXT, mech_epoch TEXT,
    kind TEXT, shape TEXT,           -- stellar|telluric ; line|edge
    wave_ref REAL,                   -- a-priori wavelength (air, A)
    x REAL, x_err REAL,              -- measured centre / edge (px)
    fwhm_px REAL, amp_frac REAL, snr REAL, rchi2 REAL,
    used INTEGER,                    -- 1 = entered the primary fit
    PRIMARY KEY (path, line)
);
CREATE TABLE IF NOT EXISTS g_dispersion (
    grism TEXT NOT NULL, mech_epoch TEXT NOT NULL,
    disp_a_per_px REAL,              -- at x_ref; signed: + = red to +x
    disp_err REAL,
    disp_minus1000 REAL,             -- |A/px| at x_ref - 1000 px
    disp_plus1000 REAL,              -- |A/px| at x_ref + 1000 px
    x_ref REAL,                      -- detector reference column
    coeffs_json TEXT,                -- [k1, k2, (k3)]: lambda = c_frame +
                                     --   sum k_j ((x - x_ref)/1000)^j
    coeff_errs_json TEXT,
    degree INTEGER,
    n_frames INTEGER, n_stars INTEGER, n_lines INTEGER, n_meas INTEGER,
    n_clipped INTEGER,
    rms_px REAL, rms_a REAL, max_abs_resid_px REAL,
    chi2 REAL, dof INTEGER, rchi2 REAL,
    o2b_wave_eff REAL, o2b_wave_err REAL,   -- O2-B edge on this scale
    stars TEXT, lines TEXT,
    status TEXT,                     -- adopted | provisional | unsolved
    note TEXT,
    PRIMARY KEY (grism, mech_epoch)
);
CREATE TABLE IF NOT EXISTS g_dispersion_resid (
    grism TEXT, mech_epoch TEXT, path TEXT, star TEXT, line TEXT,
    wave_ref REAL, x REAL, x_err REAL, resid_px REAL, resid_a REAL,
    PRIMARY KEY (path, line)
);
CREATE TABLE IF NOT EXISTS g_telluric_lambda (
    grism TEXT NOT NULL, mech_epoch TEXT NOT NULL, feature TEXT NOT NULL,
    wave_head REAL,                  -- a-priori band-head wavelength
    wave_eff REAL, wave_eff_err REAL,       -- measured on stellar scale
    scatter_a REAL, n_frames INTEGER, n_stars INTEGER,
    PRIMARY KEY (grism, mech_epoch, feature)
);
CREATE TABLE IF NOT EXISTS g_zero_point (
    path TEXT PRIMARY KEY,
    sample TEXT, target TEXT, grism TEXT, mech_epoch TEXT, night TEXT,
    anchor TEXT,                     -- halpha_em | none
    x_halpha REAL, x_halpha_err REAL, halpha_snr REAL,
    halpha_fwhm_px REAL, halpha_amp_frac REAL,
    disp_at_halpha REAL,             -- |A/px| at the Halpha pixel
    x_o2b REAL, x_o2b_err REAL, o2b_depth REAL,
    o2b_wave REAL, o2b_wave_err REAL,       -- edge wavelength, this frame
    o2b_dwave REAL,                  -- minus the calibrated edge wavelength
    sep_px REAL,                     -- |x_o2b - x_halpha|
    sep_pred_px REAL,                -- detector-coordinate solution
    sep_frac_dev REAL,               -- (sep - pred) / pred
    status TEXT
);
CREATE TABLE IF NOT EXISTS g_variance_check (
    path TEXT PRIMARY KEY,
    sample TEXT, grism TEXT, mech_epoch TEXT, detector_key TEXT,
    n_bins INTEGER,
    level_lo REAL, level_hi REAL,    -- sky above pedestal (ADU) spanned
    var_ratio_median REAL,           -- measured / predicted variance
    var_ratio_lo REAL, var_ratio_hi REAL,   -- 16th / 84th pct over bins
    ptc_gain REAL, ptc_gain_err REAL,       -- e-/ADU from var vs level
    ptc_rn_adu REAL, ptc_rn_err REAL,
    bins_json TEXT                   -- [[level, var_meas, var_pred, n]..]
);
CREATE TABLE IF NOT EXISTS g_identity (
    path TEXT PRIMARY KEY,
    sample TEXT, claimed TEXT, grism TEXT, mech_epoch TEXT, night TEXT,
    pointing_offset_deg REAL,
    old_verdict TEXT, old_reason TEXT,      -- v1 gate (header-first)
    cc_peak REAL,                    -- fingerprint correlation (high-passed)
    cc_raw REAL,                     -- same on raw spectra (pre-registered)
    verdict_prereg TEXT,             -- verdict under the pre-registered
                                     -- thresholds (kept as the record)
    cc_lag_px REAL,
    cc_second REAL,                  -- runner-up peak (ambiguity)
    halpha_em_snr REAL, halpha_ew_px REAL,
    x_halpha REAL,                   -- the Halpha hypothesis chosen
    n_candidates INTEGER,            -- emission features considered
    tio_step REAL,                   -- flux ratio across TiO 7054 head
    cont_adu REAL,                   -- continuum level (ADU/column)
    verdict TEXT,                    -- ACCEPT | REJECT
    reason TEXT,
    truth TEXT                       -- target | non_target (controls)
);
CREATE TABLE IF NOT EXISTS g_null_sky (
    path TEXT NOT NULL, method TEXT NOT NULL, row_offset INTEGER NOT NULL,
    sample TEXT, grism TEXT, mech_epoch TEXT,
    null_median_adu REAL,            -- median empty-aperture sum
    null_p16 REAL, null_p84 REAL,
    cont_median_adu REAL,            -- target continuum, same columns
    null_frac REAL,                  -- null_median / cont_median
    n_cols INTEGER,
    PRIMARY KEY (path, method, row_offset)
);
CREATE TABLE IF NOT EXISTS g_lsf (
    path TEXT PRIMARY KEY,
    sample TEXT, target TEXT, grism TEXT, mech_epoch TEXT, night TEXT,
    focpos REAL, focus_offset REAL,  -- focpos minus nightly g-band focus
    ccd_temp REAL,
    fwhm_xd_px REAL,                 -- cross-dispersion FWHM at Halpha
    fwhm_xd_min_px REAL,             -- best-focus chunk on the frame
    edge_width_px REAL,              -- O2B blue-edge 16-84% width
    line_fwhm_px REAL,               -- narrowest unresolved line
    lsf_fwhm_px REAL,                -- adopted (see lsf.py)
    lsf_fwhm_a REAL, lsf_fwhm_kms REAL, resolving_power REAL,
    off_nominal INTEGER
);
CREATE TABLE IF NOT EXISTS g_ew (
    path TEXT NOT NULL, library TEXT NOT NULL,   -- 'v1' | 'v2'
    grism TEXT, night TEXT, jd REAL, mech_epoch TEXT,
    gate TEXT,
    disp_a_per_px REAL, disp_source TEXT,
    x_halpha REAL,
    ew_a REAL, ew_err_a REAL,        -- emission positive
    ew_box_a REAL,                   -- same, plain aperture sum
    ew_px REAL,                      -- dispersion-free EW in pixels
    cont_adu REAL, cont_snr REAL,
    line_peak_adu REAL, n_sat_cols INTEGER, n_masked_line INTEGER,
    sky_method TEXT, status TEXT,
    PRIMARY KEY (path, library)
);
CREATE TABLE IF NOT EXISTS g_method_diff (
    path TEXT PRIMARY KEY, grism TEXT,
    neg_frac_poly REAL,              -- columns < -3 sigma, polynomial sky
    neg_frac_flanking REAL,          -- same, straight-line (v1) sky
    sky_method_diff REAL,            -- median |poly - flanking| / poly
    box_opt_diff REAL                -- median |box - optimal| / optimal
);
CREATE TABLE IF NOT EXISTS g_saturation (
    path TEXT PRIMARY KEY, grism TEXT, night TEXT,
    peak_adu REAL,                   -- brightest good pixel in aperture
    sat_cap_adu REAL,                -- measured cap applied
    n_sat_cols INTEGER,              -- columns with a pixel >= cap
    n_sat_cols_halpha INTEGER,       -- ... within the EW window
    hot_mask TEXT, masked_frac REAL, -- bad-pixel fraction of the aperture
    verdict TEXT                     -- clean | flag | discard
);
CREATE TABLE IF NOT EXISTS g_meta (
    key TEXT PRIMARY KEY,
    value TEXT
);
"""


def ensure_schema(con: sqlite3.Connection) -> None:
    """Create the g_* tables when absent (idempotent — IF NOT EXISTS)."""
    con.executescript(G_SCHEMA)


def connect_grism(path=None) -> sqlite3.Connection:
    """Open (creating if needed) the grism track's own database and make
    sure the calibration-library tables exist.  WAL is NOT enabled: the
    file lives on a Dropbox volume, where -wal/-shm sidecars are the
    documented way to lose a transaction."""
    from .config import GRISM_DB
    target = path or GRISM_DB
    import os
    os.makedirs(os.path.dirname(str(target)), exist_ok=True)
    con = sqlite3.connect(str(target), timeout=600)
    con.executescript(G2_SCHEMA)
    _add_missing_columns(con)
    return con


def _add_missing_columns(con: sqlite3.Connection) -> None:
    """Bring tables created by an older schema up to date (ADD COLUMN
    only; nothing is dropped or rewritten)."""
    import re
    for block in re.findall(r"CREATE TABLE IF NOT EXISTS (\w+) \((.*?)\n\);",
                            G2_SCHEMA, flags=re.S):
        table, body = block
        have = {r[1] for r in con.execute(f"PRAGMA table_info({table})")}
        for line in body.splitlines():
            m = re.match(r"\s*(\w+)\s+(TEXT|REAL|INTEGER|BLOB)\b", line)
            if m and m.group(1) not in have and m.group(1) != "PRIMARY":
                con.execute(f"ALTER TABLE {table} ADD COLUMN "
                            f"{m.group(1)} {m.group(2)}")
    con.commit()


def connect_manifest_ro(path) -> sqlite3.Connection:
    """The shared manifest, READ-ONLY (ground rule: only ``foundation-s0``
    writes it this wave).  ``mode=ro`` makes a stray INSERT an error
    instead of a race."""
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=120)


def upsert(con: sqlite3.Connection, table: str, row: dict) -> None:
    """INSERT OR REPLACE one row from a plain dict; keys that are not
    columns of ``table`` are an error (a misspelt column must not vanish
    silently), missing columns become NULL."""
    cols = [r[1] for r in con.execute(f"PRAGMA table_info({table})")]
    unknown = set(row) - set(cols)
    if unknown:
        raise KeyError(f"{table}: unknown column(s) {sorted(unknown)}")
    con.execute(
        f"INSERT OR REPLACE INTO {table} ({','.join(cols)}) "
        f"VALUES ({','.join('?' * len(cols))})",
        [row.get(c) for c in cols])


def set_g_meta(con: sqlite3.Connection, key: str, value) -> None:
    con.execute("INSERT OR REPLACE INTO g_meta (key, value) VALUES (?, ?)",
                (key, str(value)))


def existing_keys(con: sqlite3.Connection) -> set:
    """(obs_rowid, method) pairs already recorded — the resume key set."""
    try:
        return set(con.execute(
            "SELECT obs_rowid, method FROM g_extractions"))
    except sqlite3.OperationalError:
        return set()


def insert_extraction(con: sqlite3.Connection, row: dict) -> None:
    """Insert one g_extractions row from a plain dict (missing keys become
    NULL).  ``INSERT OR REPLACE`` on the (obs_rowid, method) key makes a
    deliberate re-run of a frame overwrite its old record instead of
    stacking duplicates."""
    cols = [r[1] for r in con.execute(
        "PRAGMA table_info(g_extractions)")]
    vals = [row.get(c) for c in cols]
    con.execute(
        f"INSERT OR REPLACE INTO g_extractions ({','.join(cols)}) "
        f"VALUES ({','.join('?' * len(cols))})", vals)


def set_meta(con: sqlite3.Connection, key: str, value: str) -> None:
    con.execute("INSERT OR REPLACE INTO g_build_meta (key, value) "
                "VALUES (?, ?)", (key, str(value)))
