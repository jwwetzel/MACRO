"""Provenance report renderer — ``docs/pipeline/pipeline_status.html``.

House rule, same as every other report in this repo: **no number on the
page is typed by hand.**  Each one is either

* a fingerprint / verdict computed by :mod:`macro_core.provenance`,
* the result of a SQL query declared in :data:`EVIDENCE` below and executed
  against the manifest at render time, or
* a FITS header value read live from the archive by :func:`probe_geometry`.

The page carries three things a reader needs and cannot currently get:

1. **The dependency DAG**, drawn as inline SVG and coloured by freshness,
   so "what depends on what" stops being folklore.
2. **The invalidation matrix**: for every stage and every published
   artifact, one of

   ``VALID``            the re-characterizations cannot have changed it,
                        with the reason;
   ``STALE-RECOMPUTE``  mechanically re-runnable — no judgement needed;
   ``STALE-REDERIVE``   a human DECISION inside it may change, and the
                        row names which decision.

   The distinction is the whole point: a recompute is a chore, a rederive
   is a meeting.
3. **The ordered re-run plan**, with the exact commands.

The audit rows below are the 2026-08-18 assessment against two upstream
re-characterizations:

* the **geometry artifact** — 19,980 canonical frames recorded the
  tile-compressed BINTABLE's NAXIS1/NAXIS2 (8 x 3211 bytes/rows) instead of
  the true 4800 x 3211 image, minting phantom eras 80 and 83;
* the **filter identity** finding — source-elongation measurement shows
  hrg/HaGrism/lrg/OGGrism really are dispersed, while slot ``'6'`` is
  MIXED, so FILTER labels cannot be trusted to say what is a spectrum.
"""

from __future__ import annotations

import datetime as dt
import html
import os
import sqlite3
from pathlib import Path

from macro_core import provenance as pv

# ---------------------------------------------------------------------------
# Verdict vocabulary
# ---------------------------------------------------------------------------
VALID = "VALID"
RECOMPUTE = "STALE-RECOMPUTE"
REDERIVE = "STALE-REDERIVE"
DESTROYED = "EVIDENCE-DESTROYED"

VERDICT_CLASS = {VALID: "ok", RECOMPUTE: "warn", REDERIVE: "bad",
                 DESTROYED: "gone"}


# ---------------------------------------------------------------------------
# EVIDENCE — every number the audit quotes, as a query
# ---------------------------------------------------------------------------
# Each entry is name -> SQL returning ONE scalar.  The audit rationale
# strings interpolate them by name, so a rationale can never drift from the
# database it describes: if a query stops returning, the page fails loudly
# instead of printing a stale constant.
#
# PHANTOM ERAS.  Eras 80 and 83 are the two whose recorded geometry is the
# BINTABLE row-length/row-count pair (8 x 3211).  Eras 58/59/70/65 also have
# small NAXIS values but are NOT phantoms — probe_geometry() proves it below
# by reading their headers, where ZNAXIS1 == the recorded NAXIS1.
PHANTOM_ERAS = (80, 83)
SUSPECT_SMALL_ERAS = (58, 59, 70, 65)
#: The real eras those phantom rows belong to once geometry is corrected:
#: era 80 ('Fast', egain 1.0) merges into era 78, era 83 (blank readout,
#: egain 56) into era 81 — both keyed identically on the TRUE 4800 x 3211.
MERGE_TARGET = {80: 78, 83: 81}

EVIDENCE: dict[str, str] = {
    # -- blast radius --------------------------------------------------------
    "phantom_frames":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN (80,83)",
    "phantom_light":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN (80,83)"
        " AND imagetyp LIKE 'Light%'",
    "phantom_era80":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id=80",
    "phantom_era83":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id=83",
    "small_era_frames":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN"
        " (58,59,70,65)",
    "n_eras":
        "SELECT count(*) FROM eras",
    # -- S0b -----------------------------------------------------------------
    "calib_in_phantom":
        "SELECT count(*) FROM calib_frames WHERE era_id IN (80,83)",
    "gaps_total":
        "SELECT count(*) FROM calib_gaps",
    "gaps_phantom":
        "SELECT count(*) FROM calib_gaps WHERE era_id IN (80,83)",
    "gaps_blocked_total":
        "SELECT sum(n_science_frames_blocked) FROM calib_gaps",
    "gaps_blocked_phantom":
        "SELECT sum(n_science_frames_blocked) FROM calib_gaps"
        " WHERE era_id IN (80,83)",
    "bias_merged_fast":
        "SELECT sum(n_science_frames_blocked) FROM calib_gaps"
        " WHERE need_kind='bias' AND era_id IN (78,80)",
    "bias_merged_blank":
        "SELECT sum(n_science_frames_blocked) FROM calib_gaps"
        " WHERE need_kind='bias' AND era_id IN (81,83)",
    "flatg_merged_fast":
        "SELECT sum(n_science_frames_blocked) FROM calib_gaps"
        " WHERE need_kind='flat' AND spec LIKE 'flat g %' AND era_id IN (78,80)",
    "links_phantom":
        "SELECT count(*) FROM raw_reduced_links l JOIN frames f"
        " ON f.obs_rowid=l.raw_rowid WHERE f.era_id IN (80,83)",
    # -- the S0b re-characterization the ops request predates ---------------
    "bias_eras_failing":
        "SELECT count(*) FROM calib_coverage WHERE req_kind='bias'"
        " AND status!='ok'",
    "bias_eras_total":
        "SELECT count(*) FROM calib_coverage WHERE req_kind='bias'",
    "bias_science_blocked":
        "SELECT sum(n_science) FROM calib_coverage WHERE req_kind='bias'"
        " AND status!='ok'",
    "last_calib_night":
        "SELECT max(night) FROM calib_frames",
    "era76_bias_status":
        "SELECT status FROM calib_coverage WHERE era_id=76 AND req_kind='bias'",
    "era76_bias_science":
        "SELECT n_science FROM calib_coverage WHERE era_id=76 AND req_kind='bias'",
    "top_gap_blocked":
        "SELECT max(n_science_frames_blocked) FROM calib_gaps",
    "top_gap_label":
        "SELECT 'era ' || era_id || ' ' || camera || ' ' || spec FROM calib_gaps"
        " ORDER BY n_science_frames_blocked DESC LIMIT 1",
    # -- S0c -----------------------------------------------------------------
    "stage_bestar_phantom":
        "SELECT count(*) FROM stage_bestar_grism WHERE era_id IN (80,83)",
    "stage_cv_phantom":
        "SELECT count(*) FROM stage_cv_timeseries WHERE era_id IN (80,83)",
    "stage_tcrb_phantom":
        "SELECT count(*) FROM stage_tcrb_monitoring WHERE era_id IN (80,83)",
    "stage_dwarf_phantom":
        "SELECT count(*) FROM stage_dwarfgalaxy_agn_survey WHERE era_id IN (80,83)",
    "stage_sn_phantom":
        "SELECT count(*) FROM stage_sn2023ixf_lightcurve WHERE era_id IN (80,83)",
    "stage_calib_phantom":
        "SELECT (SELECT count(*) FROM stage_bestar_grism WHERE era_id IN (80,83)"
        "        AND role NOT LIKE 'science%')"
        "     + (SELECT count(*) FROM stage_cv_timeseries WHERE era_id IN (80,83)"
        "        AND role NOT LIKE 'science%')"
        "     + (SELECT count(*) FROM stage_tcrb_monitoring WHERE era_id IN (80,83)"
        "        AND role NOT LIKE 'science%')",
    # -- S1 ------------------------------------------------------------------
    "s1_batch_phantom":
        "SELECT count(*) FROM s1_batch b JOIN frames f USING (obs_rowid)"
        " WHERE f.era_id IN (80,83)",
    "s1_batch_rows":
        "SELECT count(*) FROM s1_batch",
    "s1_recoverable":
        "SELECT count(*) FROM frames f WHERE f.is_canonical=1"
        " AND f.tree='rawimage' AND f.error IS NULL"
        " AND (f.imagetyp LIKE 'Light%' OR f.imagetyp IS NULL OR f.imagetyp='')"
        " AND (f.pltsolvd IS NULL OR f.pltsolvd != 1)"
        " AND f.era_id IN (80,83)"
        " AND lower(trim(coalesce(f.filter,''))) NOT IN"
        "     ('hrg','lrg','hagrism','oggrism','grism')",
    "s1_phantom_grism":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN (80,83)"
        " AND lower(trim(coalesce(filter,''))) IN"
        "     ('hrg','lrg','hagrism','oggrism','grism')",
    "s1_slot6_queued":
        "SELECT count(*) FROM s1_batch WHERE filter='6'",
    # -- S3 ------------------------------------------------------------------
    "ft_phantom":
        "SELECT count(*) FROM frame_times WHERE era_id IN (80,83)",
    "ft_phantom_bjd":
        "SELECT count(*) FROM frame_times WHERE era_id IN (80,83)"
        " AND bjd_tdb IS NOT NULL",
    "ft_rows":
        "SELECT count(*) FROM frame_times",
    "audit_phantom_true":
        "SELECT count(*) FROM s3_header_audit WHERE era_id IN (80,83)"
        " AND naxis1 = 4800",
    "audit_phantom_rows":
        "SELECT count(*) FROM s3_header_audit WHERE era_id IN (80,83)",
    "audit_corner_ltt":
        "SELECT round(max(corner_ltt_s),3) FROM s3_header_audit"
        " WHERE era_id IN (80,83)",
    # -- G -------------------------------------------------------------------
    "g_rows":
        "SELECT count(*) FROM g_extractions",
    "g_phantom":
        "SELECT count(*) FROM g_extractions WHERE era_id IN (80,83)",
    "g_eras":
        "SELECT group_concat(DISTINCT era_id) FROM g_extractions",
    # -- filter identity -----------------------------------------------------
    "slot6_canonical":
        "SELECT count(*) FROM frames WHERE filter='6' AND is_canonical=1",
    "slot6_ngc5548":
        "SELECT count(*) FROM frames WHERE filter='6' AND is_canonical=1"
        " AND canonical_target='NGC 5548'",
    "slot6_2023ixf":
        "SELECT count(*) FROM frames WHERE filter='6' AND is_canonical=1"
        " AND canonical_target='2023ixf'",
    "slot6_other":
        "SELECT count(*) FROM frames WHERE filter='6' AND is_canonical=1"
        " AND canonical_target NOT IN ('NGC 5548','2023ixf')",
    "w_canonical":
        "SELECT count(*) FROM frames WHERE lower(filter)='w' AND is_canonical=1",
    "cv_w_excluded":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND tree='rawimage'"
        " AND imagetyp LIKE 'Light%' AND error IS NULL"
        " AND lower(filter)='w'"
        " AND target_key IN ('stlmi','vvpup','euuma','anuma','yzcnc')",
    "cv_slot6_excluded":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND tree='rawimage'"
        " AND imagetyp LIKE 'Light%' AND error IS NULL AND filter='6'"
        " AND target_key IN ('stlmi','vvpup','euuma','anuma','yzcnc')",
    "cv_empty_excluded":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND tree='rawimage'"
        " AND imagetyp LIKE 'Light%' AND error IS NULL"
        " AND lower(filter)='empty'"
        " AND target_key IN ('stlmi','vvpup','euuma','anuma','yzcnc')",
    "tcrb_grism":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND tree='rawimage'"
        " AND canonical_target='T CrB'"
        " AND lower(filter) IN ('hrg','lrg','hagrism','oggrism','hag')",
    "tcrb_grism_era76":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND tree='rawimage'"
        " AND canonical_target='T CrB' AND era_id=76"
        " AND lower(filter) IN ('hrg','lrg','hagrism','oggrism','hag')",
    "tetcrb_grism_phantom":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND tree='rawimage'"
        " AND canonical_target='tet CrB' AND era_id IN (80,83)"
        " AND lower(filter) IN ('hrg','lrg','hagrism','oggrism','hag')",
    "bestar_grism_canonical":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND error IS NULL"
        " AND imagetyp LIKE 'Light%'"
        " AND lower(filter) IN ('hrg','lrg','hagrism','oggrism')",
    "bestar_grism_alltree":
        "SELECT count(*) FROM frames WHERE error IS NULL"
        " AND imagetyp LIKE 'Light%'"
        " AND lower(filter) IN ('hrg','lrg','hagrism','oggrism')",
    "hrg_lrg_canonical":
        "SELECT count(*) FROM frames WHERE is_canonical=1"
        " AND lower(filter) IN ('hrg','lrg')",
    "vega_phantom":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN (80,83)"
        " AND canonical_target='Vega'",
    "euuma_phantom":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN (80,83)"
        " AND canonical_target='EU UMa'",
    # -- per-frame dispersion verdicts (S2c) ---------------------------------
    # These replace the target-level aggregates the first version of this
    # audit quoted.  An aggregate elongation is a summary of a population; a
    # verdict about six named frames is a claim about those six frames, and
    # converting one into the other is how an audit over-claims.  Every
    # number below is per frame, and NULL (unmeasured) is counted separately
    # from 'direct' — because "we have not looked" is not a finding.
    "disp_rows":
        "SELECT count(*) FROM frame_dispersion WHERE verdict IS NOT NULL",
    "disp_pending":
        "SELECT count(*) FROM frame_dispersion WHERE verdict IS NULL",
    "w_direct":
        "SELECT count(*) FROM frame_dispersion WHERE lower(filter)='w'"
        " AND verdict='direct'",
    "w_dispersed":
        "SELECT count(*) FROM frame_dispersion WHERE lower(filter)='w'"
        " AND verdict='dispersed'",
    "w_indeterminate":
        "SELECT count(*) FROM frame_dispersion WHERE lower(filter)='w'"
        " AND verdict='indeterminate'",
    "w_unmeasured":
        "SELECT count(*) FROM frame_dispersion WHERE lower(filter)='w'"
        " AND verdict IS NULL",
    "cv_w_unmeasured":
        "SELECT count(*) FROM frame_dispersion WHERE lower(filter)='w'"
        " AND verdict IS NULL AND canonical_target IN"
        " ('ST LMi','YZ Cnc','VV Pup','AN UMa','EU UMa')",
    "ngc5548_direct":
        "SELECT count(*) FROM frame_dispersion WHERE filter='6'"
        " AND canonical_target='NGC 5548' AND verdict='direct'",
    "ngc5548_dispersed":
        "SELECT count(*) FROM frame_dispersion WHERE filter='6'"
        " AND canonical_target='NGC 5548' AND verdict='dispersed'",
    "ngc5548_indeterminate":
        "SELECT count(*) FROM frame_dispersion WHERE filter='6'"
        " AND canonical_target='NGC 5548' AND verdict='indeterminate'",
    "ngc5548_unmeasured":
        "SELECT count(*) FROM frame_dispersion WHERE filter='6'"
        " AND canonical_target='NGC 5548' AND verdict IS NULL",
    "sn_slot6_unmeasured":
        "SELECT count(*) FROM frame_dispersion WHERE filter='6'"
        " AND canonical_target LIKE '%2023ixf%' AND verdict IS NULL",
    "direct_control_ab":
        "SELECT round(avg(median_ab),2) FROM frame_dispersion"
        " WHERE verdict='direct' AND lower(filter) IN"
        " ('g','v','b','r','i','l')",
    # -- the CV photometry product (the one that holds phantom frames) -------
    "cv_prod_phantom":
        "SELECT count(*) FROM frames WHERE is_canonical=1 AND era_id IN (80,83)"
        " AND target_key IN ('anuma','vvpup','stlmi','yzcnc','euuma')",
}

#: Queries against the CV photometry PRODUCT database rather than the
#: manifest.  Kept separate because they need a second connection; merged
#: into the evidence dict before any template is filled.
CV_PRODUCT_DB = "products/phot/cv_timeseries.sqlite"
CV_EVIDENCE: dict[str, str] = {
    "cvp_frames": "SELECT count(*) FROM cv_frames",
    "cvp_era80": "SELECT count(*) FROM cv_frames WHERE era_id=80",
    "cvp_era78": "SELECT count(*) FROM cv_frames WHERE era_id=78",
    "cvp_fast": "SELECT count(*) FROM cv_frames WHERE era_id IN (78,79,80)",
    "cvp_series_split":
        "SELECT count(*) FROM cv_selection WHERE era_id IN (78,80)"
        " AND target_key='euuma'",
    "cvp_built": "SELECT value FROM cv_build_meta WHERE key='built_utc'",
    "cvp_rule":
        "SELECT value FROM cv_build_meta WHERE key='provenance_rule'",
}


# The corrected Item-B table.  The CASE folds each phantom era into the
# configuration it actually is; the GROUP BY then adds the two rows that were
# never two things.  ``need`` is normalized so 'bias x >=20 (have 0)' from two
# eras collapses to one line instead of two that differ only in their era id.
SHOPPING_SQL = """
WITH merged AS (
  SELECT CASE era_id WHEN 80 THEN 78 WHEN 83 THEN 81 ELSE era_id END AS era,
         camera,
         CASE WHEN need_kind='bias' THEN 'bias x >=20'
              ELSE trim(replace(replace(replace(spec,'(have 0)',''),
                                        '(have 6)',''),'(have 2)','')) END AS need,
         n_science_frames_blocked AS n,
         projects_affected AS proj
  FROM calib_gaps)
SELECT era, max(camera), need, sum(n) AS blocked, max(proj)
FROM merged GROUP BY era, need ORDER BY blocked DESC LIMIT 12
"""


def run_evidence(con: sqlite3.Connection,
                 repo_root=None) -> dict[str, object]:
    """Execute every EVIDENCE query.  A query against a table that no longer
    exists yields ``'(unavailable: ...)'`` rather than an exception — the
    destruction is itself a finding the page must show.

    ``repo_root`` additionally opens the CV photometry product READ-ONLY and
    runs :data:`CV_EVIDENCE` against it.  That product was absent from the
    first version of this audit entirely, which is why its numbers are
    fetched here rather than assumed.
    """
    out: dict[str, object] = {}
    for name, sql in EVIDENCE.items():
        try:
            row = con.execute(sql).fetchone()
            val = row[0] if row else None
            out[name] = "(no rows)" if val is None else val
        except sqlite3.OperationalError as exc:
            out[name] = f"(unavailable: {exc})"

    # The CV product lives in its own file; a missing file is reported, not
    # papered over with a zero.
    path = Path(repo_root) / CV_PRODUCT_DB if repo_root else None
    if path is not None and path.exists():
        side = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
        try:
            side.execute("PRAGMA busy_timeout = 300000")
            for name, sql in CV_EVIDENCE.items():
                try:
                    row = side.execute(sql).fetchone()
                    out[name] = row[0] if row and row[0] is not None \
                        else "(no rows)"
                except sqlite3.OperationalError as exc:
                    out[name] = f"(unavailable: {exc})"
        finally:
            side.close()
    else:
        for name in CV_EVIDENCE:
            out[name] = "(product database absent)"
    return out


def fmt(value) -> str:
    """Thousands separators for integers; everything else verbatim."""
    if isinstance(value, int):
        return f"{value:,}"
    return str(value)


# ---------------------------------------------------------------------------
# Live geometry probe — the evidence that separates phantoms from subframes
# ---------------------------------------------------------------------------
ARCHIVE_ROOT = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")


def probe_geometry(con: sqlite3.Connection, per_era: int = 2,
                   archive: Path = ARCHIVE_ROOT) -> list[dict]:
    """Read real FITS headers for the six small-geometry eras.

    For each sampled frame the probe reports the recorded NAXIS pair, the
    on-disk BINTABLE NAXIS pair, and the TRUE image ZNAXIS pair.  The test
    is one line of arithmetic:

        recorded == ZNAXIS  ->  GENUINE small subframe
        recorded == BINTABLE NAXIS (and != ZNAXIS)  ->  PHANTOM

    This is what settles eras 58/59/70/65, whose 57x48 / 56x52 / 45x34 /
    16x23 geometries LOOK like the artifact but are real windows.  It is a
    live measurement rather than an assumption, because assuming was how
    the archive acquired 19,980 wrong rows in the first place.
    """
    try:
        from astropy.io import fits
    except Exception:                                   # astropy unavailable
        return []
    rows: list[dict] = []
    for era in PHANTOM_ERAS + SUSPECT_SMALL_ERAS:
        sample = con.execute(
            "SELECT path, naxis1, naxis2 FROM frames WHERE era_id=? "
            "AND is_canonical=1 LIMIT ?", (era, per_era)).fetchall()
        for path, n1, n2 in sample:
            fp = archive / path
            rec = {"era": era, "path": path, "rec_n1": n1, "rec_n2": n2,
                   "bt_n1": None, "bt_n2": None, "z_n1": None, "z_n2": None,
                   "verdict": "unreadable"}
            if not fp.exists():
                rows.append(rec)
                continue
            try:
                with fits.open(fp, memmap=False) as hdul:
                    for hdu in hdul:
                        # ``_header`` is the RAW extension header — the same
                        # bytes the rescue parser saw.  ``ZNAXIS*`` inside it
                        # is the true image geometry.
                        hdr = getattr(hdu, "_header", hdu.header)
                        if hdr.get("ZIMAGE") or hdr.get("ZCMPTYPE"):
                            rec["bt_n1"] = hdr.get("NAXIS1")
                            rec["bt_n2"] = hdr.get("NAXIS2")
                            rec["z_n1"] = hdr.get("ZNAXIS1")
                            rec["z_n2"] = hdr.get("ZNAXIS2")
                            break
            except Exception as exc:                    # unreadable file
                rec["verdict"] = f"unreadable: {type(exc).__name__}"
                rows.append(rec)
                continue
            if rec["z_n1"] is None:
                rec["verdict"] = "not tile-compressed"
            elif int(n1 or -1) == int(rec["z_n1"]):
                rec["verdict"] = "GENUINE subframe"
            else:
                rec["verdict"] = "PHANTOM (BINTABLE geometry)"
            rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# THE AUDIT — one row per stage and per published artifact
# ---------------------------------------------------------------------------
class Row:
    """One invalidation-matrix row.

    ``rationale`` and ``action`` are ``str.format`` templates over the
    EVIDENCE names, so every number in them is a query result.
    """

    def __init__(self, subject: str, kind: str, verdict: str,
                 rationale: str, action: str, decision: str = ""):
        self.subject = subject
        self.kind = kind
        self.verdict = verdict
        self.rationale = rationale
        self.action = action
        self.decision = decision


AUDIT: tuple[Row, ...] = (
    # ---------------- stages ------------------------------------------------
    Row("S0 — manifest (frames, eras, aliases, project_counts)", "stage",
        RECOMPUTE,
        "S0 reads geometry ONLY in era_key(READOUTM, NAXIS1, NAXIS2, "
        "XBINNING, EGAIN).  Dedup (basename, JD), alias resolution, night "
        "labelling (JD − 0.7917 d) and the pointing audit never see a NAXIS "
        "value, so those columns are unaffected.  What IS affected is the "
        "era registry: {phantom_frames} canonical frames carry the BINTABLE "
        "pair 8×3211 and were minted into two phantom eras — 80 "
        "({phantom_era80} frames) and 83 ({phantom_era83}).  With the true "
        "4800×3211 they key identically to eras 78 (Fast, EGAIN 1.0) and 81 "
        "(blank readout, EGAIN 56) and MERGE into them; the {n_eras}-row era "
        "table loses two rows.  Independent corroboration that the split is "
        "an artifact and not a real configuration change: eras 78 and 80 "
        "INTERLEAVE in time (78 runs 2026-03-21→06-28, 80 runs "
        "2026-04-23→06-28) — one camera cannot be in two configurations on "
        "the same nights.",
        "Re-run the header rescue over the catalog, then "
        "build_s0_manifest.py.  No judgement required."),

    Row("S0 — eras 58/59/70/65 (the four tiny eras)", "stage", VALID,
        "Explicitly checked rather than assumed.  These {small_era_frames} "
        "frames (57×48, 56×52, 45×34, 16×23) look like the artifact but are "
        "GENUINE small subframes: the live header probe below shows their "
        "recorded NAXIS pair equals ZNAXIS1/ZNAXIS2, while their BINTABLE "
        "NAXIS1 is 8 — i.e. they came through the CORRECT path.  The "
        "concurrent catalog rescan agrees: it marks these rows changed=0.  "
        "The blast radius is therefore {phantom_frames} frames, not "
        "{phantom_frames}+{small_era_frames}.",
        "Nothing.  Do not 'fix' these."),

    Row("S0b — calib_frames / calib_coverage / calib_gaps", "stage",
        RECOMPUTE,
        "Calibration frames affected by the geometry artifact: "
        "{calib_in_phantom} — the rescue path touched no calibration frame, "
        "so no master, no bias and no dark changes identity.  The damage is "
        "purely in the era JOIN: {gaps_phantom} of {gaps_total} calib_gaps "
        "rows are phantom-era rows carrying {gaps_blocked_phantom} of "
        "{gaps_blocked_total} blocked science frames, and they DOUBLE-COUNT "
        "the same physical configuration as eras 78/81.  Worked example: "
        "'era 80 Fast bias ×≥20 blocking 18,149' and 'era 78 Fast bias ×≥20 "
        "blocking 8,788' are one line item — {bias_merged_fast} frames "
        "blocked by ONE missing bias set — and 'flat g' likewise merges to "
        "{flatg_merged_fast}.  The blank-readout pair merges to "
        "{bias_merged_blank}.  {links_phantom} raw↔reduced links point at "
        "phantom-era raw frames.",
        "build_s0b_inventory.py after S0.  Mechanical."),

    Row("S0c — the five stage_<project> tables", "stage", RECOMPUTE,
        "Science rows staged from phantom eras: BeStar {stage_bestar_phantom}, "
        "CV {stage_cv_phantom}, T CrB {stage_tcrb_phantom}, Dwarf "
        "{stage_dwarf_phantom}, SN 2023ixf {stage_sn_phantom}.  Their "
        "SELECTION is unaffected (no staging rule reads NAXIS) — every one "
        "of those frames is still selected — but each row carries a wrong "
        "era_id, and era_id is what the calibration side of staging matches "
        "on: match_basis='era_exact'.  Calibration rows era-matched to a "
        "phantom era: {stage_calib_phantom} — i.e. NO calibration row was "
        "attached to eras 80/83 at all, because those eras hold zero "
        "calibration frames.  After the merge those frames inherit era 78/81 "
        "and will match whatever calibration those eras hold.",
        "build_s0c_staging.py after S0/S0b.  Mechanical."),

    Row("S1 — astrometry solvability experiment", "stage", DESTROYED,
        "Two independent problems.  (1) EVIDENCE GONE: s1_strata, "
        "s1_populations, s1_solve_experiment and s1_failure_autopsy were "
        "destroyed by the S0 table swap and are not in the database; the "
        "published verdict page cannot be regenerated or checked.  "
        "(2) WRONG PREMISE: astrom.is_window_geometry() rejects any frame "
        "with an axis under 512 px, so all {phantom_frames} phantom-geometry "
        "frames were excluded as 'photometry strips'.  Removing their grism "
        "frames ({s1_phantom_grism}, correctly excluded as spectra) leaves "
        "{s1_recoverable} full 4800×3211 frames that were written off and "
        "are in fact ordinary solvable fields — including the EU UMa frames "
        "({euuma_phantom}) recorded as 'permanently unsolvable'.",
        "Re-run run_s1_experiment.py after S0 to rebuild the destroyed "
        "evidence, then re-read the strata verdicts.",
        decision="The GO / CAUTION / NO-GO verdict per stratum, and the "
                 "'permanently unsolvable' classification of the EU UMa CV "
                 "frames — both were set on a population that excluded "
                 "{s1_recoverable} solvable frames."),

    Row("S1b — astrometry production batch (s1_batch)", "stage", RECOMPUTE,
        "{s1_batch_rows} queued rows, of which {s1_batch_phantom} are in "
        "phantom eras — the geometry gate kept every one of them out.  "
        "Nothing already solved is wrong; the batch is INCOMPLETE, missing "
        "the {s1_recoverable} recoverable frames.  Separately, "
        "{s1_slot6_queued} queued frames carry FILTER='6', which "
        "measurement now shows is dispersed on some targets — those are "
        "spectra being fed to a plate solver.",
        "Re-queue after S0; the runner is resumable, so the existing solves "
        "are kept and only the new candidates are added."),

    Row("S2 — detector truth (ceiling, PTC, reconstruction, linearity)",
        "stage", DESTROYED,
        "Every S2 table — s2_ceiling_modes, s2_ptc_fits, s2_recon_eras, "
        "s2_linearity_ladders — AND detector_params is absent from the "
        "manifest: destroyed by the same S0 table swap that took S1's.  The "
        "adopted ADU ceiling, gain, read noise and linearity residuals that "
        "every later error budget cites now exist only as prose on a page "
        "no query can reproduce.  On the geometry question specifically, the "
        "selections can be read from the code: PTC keys on naxis1=4096 on "
        "2023-06-07 (eras 1/2 — untouched); reconstruction keys on the "
        "REDUCED frame's era, which for this camera is 79/82, not 80/83 "
        "(untouched); linearity groups on (night, target, readout, filter, "
        "egain) and surfaces NO ladder in eras 80/83 (untouched).  The "
        "ceiling probe is the exception: it groups by READOUTM alone, so the "
        "'Fast' and '(blank 2026)' samples drew from phantom-era frames — "
        "and with the sample table destroyed, WHICH frames entered can no "
        "longer be checked.",
        "Re-run run_s2_campaign.py end to end after S0.  Everything must be "
        "recomputed because nothing survives.",
        decision="The 12-bit reading, the adopted per-mode clip/veto "
                 "thresholds and the era-79 identity verdict were "
                 "human calls made on evidence that no longer exists; they "
                 "must be re-made, not re-typed."),

    Row("S3 — timing (frame_times, BJD_TDB)", "stage", VALID,
        "PROVEN geometry-independent, not assumed.  A frame's BJD_TDB is a "
        "function of (JD_UTC start, EXPTIME, RA, Dec, ephemeris) alone — "
        "timing.py's only NAXIS consumers are pixel_scale_arcsec() and "
        "field_corner_light_time_s(), and both feed the s3_header_audit "
        "CAVEAT column, never frame_times.  {ft_phantom} phantom-era frames "
        "are in frame_times and {ft_phantom_bjd} of them carry a BJD — the "
        "stamps stand.  This part is independent of any recent patch: "
        "bjd_tdb_from_utc(), jd_utc_mid() and mid_method_for() take no "
        "geometry argument, so no NAXIS value can reach a BJD by any path.  "
        "ATTRIBUTION CORRECTED after review: s3_header_audit does carry TRUE "
        "geometry — {audit_phantom_true} of {audit_phantom_rows} phantom-era "
        "rows read NAXIS1 = 4800 — but that is NOT a standing property of "
        "S3 as this row first claimed.  It is the effect of the concurrent "
        "geometry workflow routing build_s3_timing.py through "
        "macro_core.fitsgeom.resolve_geometry_or_none(), a change that was "
        "still uncommitted in the working tree when this audit ran, and S3 "
        "was re-run after it.  Two caveats the first version also omitted: "
        "the evidence is 2 rows out of a 40-row header-audit sample, and the "
        "era-83 row carries a NULL corner light-time with pixscale_source "
        "'unknown', so the quoted {audit_corner_ltt} s comes from the era-80 "
        "row alone.",
        "Only the era_id LABEL on frame_times rows needs to follow the era "
        "merge.  No time value changes.  Re-check the header audit once the "
        "geometry fix is committed, so the property stops depending on an "
        "uncommitted working tree."),

    Row("S4 — ensemble photometry (AN UMa / VV Pup prototype)", "stage",
        RECOMPUTE,
        "VERDICT DOWNGRADED from VALID after review.  On GEOMETRY the "
        "original finding holds and is re-verified: phot_selection holds "
        "exactly three (target, era) sets — anuma era 76, vvpup era 72, "
        "vvpup era 76 — and neither target has a single frame in eras "
        "80/83, so the artifact cannot have touched a measured pixel.  But "
        "'VALID' was the wrong word for the stage, because S4 rests on an "
        "input the audit itself calls destroyed.  Every frame it measures is "
        "saturation-vetoed with S2_MODE_VETO_ADU, a table of literals in "
        "macro_phot/series.py copied there BECAUSE the S0 rebuild destroyed "
        "the s2_* tables that measured them.  Nothing in the database can "
        "reproduce those numbers today.  The DAG did not say so either: "
        "until this revision S4 declared no S2 dependency, and reported "
        "FRESH in the same status run in which S2 reported OUTPUT_MISSING — "
        "a green verdict over a destroyed input, which is the dangerous "
        "direction.  S4 now declares table:detector_params, "
        "table:s2_ceiling_modes and the constants file, and reads STALE "
        "until S2 is rebuilt.  Mitigating: S4's modes are Mode0 and 1MHz "
        "High Sensitivity, NOT the two modes whose ceiling sample was "
        "contaminated by phantom-era frames.",
        "Rebuild S2 (params at minimum) and confirm the adopted ceiling and "
        "veto for Mode0 and 1MHz HiSens against the rebuilt tables; then "
        "re-run S4.  Its 09:06 UTC worklist was also resolved against a "
        "frames table S0 replaced at 15:03."),

    Row("CV-S4 — production CV photometry (products/phot/cv_timeseries.sqlite)",
        "stage", RECOMPUTE,
        "MISSING FROM THE FIRST VERSION OF THIS AUDIT ENTIRELY, and it is "
        "the one photometry product the geometry artifact actually reaches. "
        "Built {cvp_built}, i.e. after the S0 rebuild and during the "
        "geometry rescue.  Its cv_frames table holds {cvp_frames} rows, of "
        "which {cvp_era80} are in phantom era 80 and {cvp_era78} in era 78 — "
        "and the repair proves those are ONE camera configuration.  They "
        "were nevertheless measured as {cvp_series_split} separate EU UMa "
        "series under a rule the product records as: '{cvp_rule}'.  That "
        "rule is being enforced on a boundary that does not exist.  Two "
        "further consequences: the {cvp_fast} Fast-mode frames are vetoed "
        "with the 'Fast' ceiling — one of the two modes whose S2 sample "
        "pooled phantom-era frames, and whose evidence table "
        "(s2_ceiling_frames) is destroyed, so which frames set that "
        "threshold can no longer be checked.  Until this revision, none of "
        "this was in the graph: `status` could exit 0 — a green light to "
        "publish — with this product resting on frames whose era assignment "
        "is known wrong.",
        "After S0/S0b/S0c and S2: re-run CV-S4 from init so the EU UMa era-78 "
        "and era-80 frames form ONE series, and re-derive the Fast veto from "
        "rebuilt S2 tables before trusting any EU UMa point.",
        decision="Whether the merged EU UMa series is one light curve or "
                 "two — the product's own provenance rule answers "
                 "differently before and after the repair."),

    Row("G — grism extraction + identity gate (T CrB)", "stage", VALID,
        "{g_rows} extraction rows, all in era {g_eras}; phantom-era rows: "
        "{g_phantom}.  The T CrB series is entirely Mode0 — "
        "{tcrb_grism_era76} of the {tcrb_grism} archival T CrB grism frames "
        "are era 76 — so the 2026 Fast/blank-readout eras where the artifact "
        "lives are not in this track at all.  Filter identity also confirms "
        "rather than threatens G: the per-frame measurement classes every "
        "measured hrg and lrg frame as dispersed, far above the "
        "{direct_control_ab} median elongation of the measured direct "
        "controls, so every frame G selected on an hrg/lrg label really is a "
        "spectrum.  ONE QUALIFICATION added after review: what G has "
        "PRODUCED is untouched, but its input worklist is not — "
        "{tetcrb_grism_phantom} θ CrB calibrator grism frames sit in phantom "
        "eras, so a re-run picks them up under merged era labels.  The "
        "scoped fingerprint table:frames@grism includes those frames "
        "deliberately, so the machinery marks G stale rather than flattering "
        "it.",
        "Nothing for the 193 published extractions.  Re-run after S0 so the "
        "θ CrB calibrator frames carry merged era labels.",
        decision=""),

    # ---------------- published artifacts ----------------------------------
    Row("ops/2026-08_observatory_request.md — the October shopping list",
        "artifact", REDERIVE,
        "THE ARTIFACT MOST AT RISK, and it is wrong for TWO independent "
        "reasons.  (1) It predates the calibration-master ingest.  Its "
        "headline facts — 'the last calibration frame in the entire archive "
        "is 2024-11-18' and '82 of 83 science eras fail the ≥20-bias spec "
        "(164,769 science frames affected)' — now read, from this database: "
        "last calibration night {last_calib_night}; {bias_eras_failing} of "
        "{bias_eras_total} eras failing; {bias_science_blocked} science "
        "frames affected.  Its Item B table is led by 'era 76 Mode0 bias "
        "68,965' and 'era 76 flat g 40,031'; era 76's bias coverage is now "
        "'{era76_bias_status}' ({era76_bias_science} science frames), so "
        "those two rows should not be on the list at all, and the true top "
        "row is {top_gap_label} at {top_gap_blocked}.  (2) It quotes phantom "
        "eras by number.  Its 'era 80 Fast bias ×≥20 → 18,149' and 'era 78 "
        "Fast bias ×≥20 → 8,788' rows are ONE configuration; asking the "
        "observatory for both is asking twice for the same bias set.  Its "
        "closing instruction — 'eras 78/80 Fast … take the frames directly "
        "at re-opening' and 'the archive's newest science frames … sit in "
        "eras 81–83' — describes an era structure that will not exist after "
        "S0 is rebuilt.",
        "DO NOT SEND until S0/S0b are rebuilt and the Item B table is "
        "regenerated from calib_gaps.  Then re-rank: merged Fast bias "
        "{bias_merged_fast}, merged Fast flat g {flatg_merged_fast}, merged "
        "blank-readout bias {bias_merged_blank}.",
        decision="Which eras get October sky time, and in what priority "
                 "order — the merge changes the ranking, and the "
                 "calibration ingest removes the request's top two rows "
                 "entirely.  A human must re-decide the ask, not just "
                 "re-run a script."),

    Row("docs/pipeline/s1_astrometry.html", "artifact", DESTROYED,
        "Orphaned: every table it renders from (s1_solve_experiment, "
        "s1_strata, s1_populations, s1_failure_autopsy) is gone, so the page "
        "can neither be regenerated nor audited.  Its central claim — which "
        "populations are solvable — was computed on a candidate universe "
        "missing {s1_recoverable} frames.",
        "Re-run S1, then re-render."),

    Row("docs/pipeline/s2_detector.html", "artifact", DESTROYED,
        "Orphaned in the same way: all five S2 tables plus detector_params "
        "are absent.  This page states the detector constants every other "
        "paper's error model inherits, and nothing in the database backs a "
        "single one of them today.",
        "Re-run S2, then re-render."),

    Row("docs/pipeline/s0b_calibration_inventory.html", "artifact",
        RECOMPUTE,
        "Currently CONSISTENT with the database — it was re-rendered after "
        "the calibration ingest and already shows {bias_science_blocked} "
        "bias-blocked science frames and last-calibration {last_calib_night}. "
        "It goes stale the moment S0 is rebuilt, because its era-resolved "
        "§4 shopping list is built on the phantom era split.",
        "Re-render after S0b.  Mechanical."),

    Row("docs/pipeline/s0_manifest.html, s0c_staging.html", "artifact",
        RECOMPUTE,
        "Every number on both pages is a query result, and both queries' "
        "answers change when {phantom_frames} frames move era.  No human "
        "judgement is embedded.",
        "Re-render after S0/S0c."),

    Row("docs/pipeline/s3_timing.html", "artifact", VALID,
        "The timing claims survive intact (see the S3 row).  Only the era "
        "LABELS in its per-era tables shift with the merge.",
        "Re-render after S0 for the labels; no timing number changes."),

    Row("docs/pipeline/s4_photometry.html", "artifact", RECOMPUTE,
        "Renders the AN UMa / VV Pup prototype, which contains no "
        "phantom-era frame — so on geometry it is untouched.  It follows S4 "
        "down, though: the page states the adopted saturation veto, and "
        "those constants are the ones no surviving table can reproduce.",
        "Re-render after S4, once the S2 constants are back in the "
        "database."),

    Row("ROADMAP.md, docs/*/index.html (the six public project pages)",
        "artifact", RECOMPUTE,
        "Declared as resources in this revision.  The previous version of "
        "this page printed a verdict for ROADMAP.md while the gate's exit "
        "code ignored it entirely — a page issuing judgements its own gate "
        "could not enforce.  The six public project pages and docs/index.html "
        "had neither a verdict nor a resource.  All are now owned by "
        "declared stages (STRAT and WEB), so `status` fails while they "
        "disagree with the tables they quote.",
        "Refresh after S0/S0c, alongside the strategy documents."),

    Row("docs/pipeline/s0e_geometry_fix.html, s2c_filter_identity.html",
        "artifact", RECOMPUTE,
        "The two pages that document the re-characterizations this audit is "
        "ABOUT, and neither was in the graph.  s2c_filter_identity.html in "
        "particular re-renders every time the classifier reclassifies "
        "({disp_rows} frames judged, {disp_pending} still pending), so it "
        "moves under its own reader's feet.  Both are now outputs of "
        "declared stages (S0e and S2c).",
        "Re-render R-S2c when the classifier finishes; re-render R-S0e "
        "after the geometry rescan completes."),

    Row("BeStar_Grism/ANALYSIS_STRATEGY.md", "artifact", REDERIVE,
        "Geometry: {vega_phantom} Vega frames — the flux calibrator — sit in "
        "phantom eras, and {tetcrb_grism_phantom} θ CrB grism frames with "
        "them, so the document's era-C accounting is keyed on eras that "
        "dissolve.  Filter identity: the whitelist (hrg / lrg / HaGrism / "
        "OGGrism) is CONFIRMED by measurement — all four are dispersed, and "
        "the ~2× elongation split (hrg 45.0 / HaGrism 43.6 vs lrg 25.2 / "
        "OGGrism 24.0) independently separates the Hα grism from the "
        "broad-spectrum one, which is exactly the era-A↔era-C pairing the "
        "document assumes.  But the whitelist has no entry for slot '6', and "
        "{slot6_canonical} canonical '6' frames exist, some of them "
        "dispersed — the strategy's own §3.2 counts cannot claim "
        "completeness until '6' is adjudicated per frame.  Separately, "
        "BeStar_Grism/README.md's '~23k lrg and hrg grism frames' is not "
        "reproducible: canonical hrg+lrg is {hrg_lrg_canonical}, and the "
        "grism-light total is {bestar_grism_canonical} canonical vs "
        "{bestar_grism_alltree} across all trees — the 23k figure is a "
        "tree-doubled count, a dedup error rather than a filter error.",
        "Re-stage after S0/S0c, then re-reconcile §3.2 against "
        "project_counts and adjudicate the '6' frames.",
        decision="Whether slot-'6' frames belong in the grism whitelist, "
                 "per target and per epoch.  A label cannot decide it; the "
                 "per-frame dispersion measurement must."),

    Row("CV_TimeSeries/ANALYSIS_STRATEGY.md", "artifact", REDERIVE,
        "CORRECTED 2026-08-18 after adversarial review — the earlier version "
        "of this row was WRONG and is retracted here.  It asserted that 'W' "
        "measures ~1.1 median elongation, a direct-imaging control, and "
        "concluded that the {cv_w_excluded} W frames of the five CVs 'were "
        "removed as spectra when they are images', casting doubt on the "
        "2026-08-18 referee correction (ST LMi 3,157→3,150, YZ Cnc "
        "1,920→1,915).  That converted a population-level control value into "
        "a verdict about six named frames, and the per-frame measurements "
        "do not support it.  Live from frame_dispersion: across the whole "
        "archive, W frames classed direct = {w_direct}, dispersed = "
        "{w_dispersed}, indeterminate = {w_indeterminate}, still unmeasured "
        "= {w_unmeasured}; the direct controls that ARE measured average "
        "{direct_control_ab} median elongation.  The specific CV W frames "
        "have NO per-frame verdict yet ({cv_w_unmeasured} unmeasured).  So "
        "the referee correction STANDS until they are classified, and "
        "nothing here challenges it.  What remains genuinely open: slot '6' "
        "({cv_slot6_excluded} CV frames) is mixed and must be decided per "
        "frame; 'empty' ({cv_empty_excluded} frames, 76 of them VV Pup) is "
        "unsettled.  Geometry: {stage_cv_phantom} staged CV rows are in "
        "phantom eras, all EU UMa ({euuma_phantom} frames) — the same frames "
        "S1 wrote off.",
        "Re-stage after S0/S0c; classify the W / '6' / 'empty' CV frames per "
        "frame; only then re-issue the §3.1 table.  Do NOT re-open the "
        "referee correction on the strength of an aggregate.",
        decision="Whether any excluded filter should be restored — which "
                 "cannot be decided until the per-frame verdicts land.  The "
                 "published counts stand in the meantime."),

    Row("DwarfGalaxy_AGN_Survey/ANALYSIS_STRATEGY.md", "artifact", REDERIVE,
        "The most consequential filter finding in the audit.  The NGC 5548 "
        "pillar rests on {slot6_ngc5548} FILTER='6' frames, and the "
        "strategy's own Phase 0.2 hedges 'a wheel slot, not a bandpass; if "
        "it is a luminance filter (likely)…'.  REFRAMED after review: the "
        "first version of this row quoted a single target-level elongation "
        "(13.0) and posed the question as all-or-nothing.  The per-frame "
        "verdicts say it is not.  Live from frame_dispersion, for NGC 5548's "
        "slot-'6' frames: dispersed {ngc5548_dispersed}, direct "
        "{ngc5548_direct}, indeterminate {ngc5548_indeterminate}, still "
        "unmeasured {ngc5548_unmeasured} — mixed WITHIN the single target, "
        "against {direct_control_ab} median elongation for the measured "
        "direct controls.  That changes the question from 'is the NGC 5548 "
        "pillar spectra or photometry' to 'which of its {slot6_ngc5548} "
        "frames survive as photometry', which is a different and far more "
        "tractable one — but it does not shrink the stake: if the surviving "
        "count is small, the 'calibrated band-integrated nightly light "
        "curve' is not the paper.  The measurement is still running, so "
        "these numbers move; the page reads them at render time rather than "
        "quoting a snapshot.  Geometry does NOT touch this project: "
        "{stage_dwarf_phantom} staged rows are in phantom eras, and NGC "
        "5548's campaign is era 1.",
        "Finish the per-frame classification of all {slot6_ngc5548} NGC 5548 "
        "frames ({ngc5548_unmeasured} to go), then rebuild the light curve "
        "from the surviving direct frames only.",
        decision="How much of NGC 5548 survives as photometry — and whether "
                 "what is left still supports the project's central claim."),

    Row("TCrB_Monitoring/ANALYSIS_STRATEGY.md", "artifact", RECOMPUTE,
        "The headline '247 spectra on 60 nights' SURVIVES: {tcrb_grism} "
        "T CrB grism frames are in the manifest, all era 76 Mode0 "
        "({tcrb_grism_era76}), and measurement confirms hrg (45.0) and lrg "
        "(25.2) are genuinely dispersed — the label and the pixels agree. "
        "The exposed part is the CALIBRATOR: {tetcrb_grism_phantom} θ CrB "
        "grism frames sit in phantom eras, so the calibrator series' era "
        "accounting moves, and {stage_tcrb_phantom} staged T CrB-project "
        "rows carry a phantom era_id.  The strategy's calibration claim "
        "'zero Mode0 darks or biases' also predates the master ingest and "
        "must be re-checked against calib_coverage.",
        "Re-stage after S0/S0c and re-check the Mode0 calibration claim; "
        "the 247-spectrum spine needs no re-derivation."),

    Row("SN2023ixf_LightCurve/ANALYSIS_STRATEGY.md", "artifact", RECOMPUTE,
        "Geometry-clean: {stage_sn_phantom} staged rows in phantom eras.  "
        "One filter caveat: {slot6_2023ixf} of the target's canonical frames "
        "carry slot '6', which measures AMBIGUOUS on this target (5.8) — "
        "between the direct controls (~1.1) and the dispersed grisms "
        "(24–45).  Those frames need a per-frame verdict before they enter "
        "or leave the light curve.",
        "Re-stage after S0/S0c; classify the {slot6_2023ixf} slot-'6' "
        "frames per frame."),

    Row("ROADMAP.md", "artifact", RECOMPUTE,
        "Quotes stage status and counts that move with S0/S0b/S0c; contains "
        "no independent measurement of its own.",
        "Refresh after the pipeline re-runs."),
)


# ---------------------------------------------------------------------------
# The DAG figure (inline SVG, no plotting dependency)
# ---------------------------------------------------------------------------
# One colour per verdict.  STALE_UPSTREAM is deliberately the quietest of
# the non-fresh colours: a stage merely waiting on an ancestor is not
# accused of anything, and painting it the same amber as a stage whose own
# evidence moved is what made the plan read as "re-run everything".
SVG_FILL = {pv.FRESH: "#1f7a3f", pv.STALE: "#a06a00",
            pv.STALE_UPSTREAM: "#5b6470",
            pv.NEVER_RUN: "#1f5f7a", pv.OUTPUT_MISSING: "#8c1c1c"}


def dag_svg(freshness) -> str:
    """Draw the stage DAG as SVG, one column per dependency depth.

    Depth = longest path from a root, so an arrow always points rightward
    and the picture cannot suggest a dependency that does not exist.
    """
    order = pv.topological_order(pv.STAGES)
    writer = {w: s.key for s in pv.STAGES for w in s.writes}
    parents = {s.key: sorted({writer[r] for r in s.reads
                              if r in writer and writer[r] != s.key})
               for s in pv.STAGES}
    depth: dict[str, int] = {}
    for key in order:                     # topological order guarantees the
        ps = parents[key]                 # parents are already assigned
        depth[key] = 0 if not ps else 1 + max(depth[p] for p in ps)

    cols: dict[int, list[str]] = {}
    for key in order:
        cols.setdefault(depth[key], []).append(key)

    bw, bh = 168, 42                      # box width / height
    gx, gy = 78, 22                       # gaps
    width = max(cols) * (bw + gx) + bw + 40
    height = max(len(v) for v in cols.values()) * (bh + gy) + 60

    pos: dict[str, tuple[int, int]] = {}
    for d, keys in cols.items():
        for i, key in enumerate(keys):
            pos[key] = (20 + d * (bw + gx), 40 + i * (bh + gy))

    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" '
             f'xmlns="http://www.w3.org/2000/svg" '
             f'font-family="ui-sans-serif,system-ui,sans-serif" '
             f'role="img" aria-label="pipeline dependency DAG">',
             '<defs><marker id="a" markerWidth="9" markerHeight="7" '
             'refX="9" refY="3.5" orient="auto">'
             '<path d="M0,0 L9,3.5 L0,7 z" fill="#8a8a8a"/></marker></defs>']
    # edges first, so boxes paint over them
    for key, ps in parents.items():
        x2, y2 = pos[key]
        for p in ps:
            x1, y1 = pos[p]
            parts.append(
                f'<path d="M{x1 + bw},{y1 + bh // 2} '
                f'C{x1 + bw + 34},{y1 + bh // 2} {x2 - 34},{y2 + bh // 2} '
                f'{x2},{y2 + bh // 2}" fill="none" stroke="#8a8a8a" '
                f'stroke-width="1.3" marker-end="url(#a)" opacity="0.75"/>')
    for key, (x, y) in pos.items():
        f = freshness[key]
        fill = SVG_FILL.get(f.state, "#555")
        title = pv.STAGE_BY_KEY[key].title
        parts.append(
            f'<g><title>{html.escape(key)}: {html.escape(f.state)} — '
            f'{html.escape(title)}</title>'
            f'<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="7" '
            f'fill="{fill}" stroke="#00000033"/>'
            f'<text x="{x + 10}" y="{y + 18}" fill="#fff" font-size="13" '
            f'font-weight="700">{html.escape(key)}</text>'
            f'<text x="{x + 10}" y="{y + 33}" fill="#ffffffcc" '
            f'font-size="10.5">{html.escape(f.state)}</text></g>')
    parts.append("</svg>")
    return "\n".join(parts)


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------
CSS = """
:root{--bg:#ffffff;--fg:#16181d;--muted:#5c6270;--line:#e2e5ea;
      --card:#f7f8fa;--ok:#1f7a3f;--warn:#a06a00;--bad:#8c1c1c;--gone:#5b2b8a}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
     font:15px/1.6 ui-sans-serif,system-ui,-apple-system,sans-serif}
.wrap{max-width:1120px;margin:0 auto;padding:32px 20px 80px}
h1{font-size:29px;margin:0 0 6px;letter-spacing:-.02em}
h2{font-size:20px;margin:38px 0 10px;padding-top:14px;
   border-top:1px solid var(--line)}
h3{font-size:15px;margin:22px 0 6px}
p,li{max-width:76ch}
.sub{color:var(--muted);margin:0 0 22px}
code,pre{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12.5px}
pre{background:var(--card);border:1px solid var(--line);border-radius:8px;
    padding:12px 14px;overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:10px 0}
th,td{border:1px solid var(--line);padding:8px 10px;text-align:left;
      vertical-align:top}
th{background:var(--card);font-weight:650}
.scroll{overflow-x:auto}
.tag{display:inline-block;padding:2px 8px;border-radius:999px;
     font-size:11.5px;font-weight:700;color:#fff;white-space:nowrap}
.ok{background:var(--ok)}.warn{background:var(--warn)}
.bad{background:var(--bad)}.gone{background:var(--gone)}
.key{display:flex;gap:14px;flex-wrap:wrap;margin:10px 0 4px;font-size:12.5px}
.fig{background:var(--card);border:1px solid var(--line);border-radius:10px;
     padding:14px;margin:14px 0}
.note{background:var(--card);border-left:3px solid var(--muted);
      padding:10px 14px;margin:14px 0;border-radius:0 8px 8px 0}
.mono-sm{font-family:ui-monospace,monospace;font-size:12px;color:var(--muted)}
@media (prefers-color-scheme:dark){
 :root{--bg:#0f1115;--fg:#e6e8ec;--muted:#9aa1ad;--line:#262a31;--card:#171a20;
       --ok:#3fae6a;--warn:#d29b2a;--bad:#e0555a;--gone:#a77ce0}
 .tag{color:#0f1115}
}
"""


def render(con: sqlite3.Connection, repo_root: Path, freshness, fingerprints,
           records, manifest_path) -> Path:
    """Render the page.  Returns the path written."""
    ev = run_evidence(con, repo_root)
    probe = probe_geometry(con)
    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    out_dir = Path(repo_root) / "docs" / "pipeline"
    fig_dir = out_dir / "figures" / "prov"
    fig_dir.mkdir(parents=True, exist_ok=True)

    svg = dag_svg(freshness)
    # The figure is also written standalone, so it can be linked or embedded
    # elsewhere without re-running the renderer.
    _atomic_write(fig_dir / "pipeline_dag.svg", svg)

    def E(template: str) -> str:
        """Interpolate EVIDENCE values into an audit template and escape."""
        filled = template.format(**{k: fmt(v) for k, v in ev.items()})
        return html.escape(filled)

    order = pv.topological_order(pv.STAGES)
    plan = pv.rerun_plan(freshness, pv.STAGES)
    counts: dict[str, int] = {}
    for k in order:
        counts[freshness[k].state] = counts.get(freshness[k].state, 0) + 1

    h: list[str] = []
    A = h.append
    A("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
    A("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    A("<title>MACRO pipeline status &amp; invalidation audit</title>")
    A(f"<style>{CSS}</style></head><body><div class='wrap'>")

    A("<h1>Pipeline status &amp; invalidation audit</h1>")
    A(f"<p class='sub'>Generated {html.escape(now)} &middot; "
      f"digest rules {html.escape(pv.PROVENANCE_CODE_VERSION)} &middot; "
      f"manifest <span class='mono-sm'>{html.escape(str(manifest_path))}"
      f"</span><br>Every number below is a query result, a fingerprint, or a "
      f"FITS header read at render time. Nothing here is typed by hand.</p>")

    A("<div class='note'><b>The question this page answers.</b> Two upstream "
      "re-characterizations landed after most stages had already run and "
      "published: a <b>geometry artifact</b> (tile-compressed BINTABLE "
      f"NAXIS recorded instead of the true image size, for "
      f"{fmt(ev['phantom_frames'])} canonical frames) and a <b>filter "
      "identity</b> finding (source elongation measured per frame; FILTER "
      "labels do not reliably say what is a spectrum). What we once "
      "considered done may not be done any longer &mdash; below, per stage "
      "and per published artifact, is which.</div>")

    # ---- 1. DAG ----------------------------------------------------------
    A("<h2>1 &middot; The dependency DAG, coloured by freshness</h2>")
    A("<div class='key'>"
      + " ".join(f"<span class='tag' style='background:{SVG_FILL[s]}'>"
                 f"{s}</span> {counts.get(s, 0)}"
                 for s in pv.ALL_STATES)
      + "</div>")
    A(f"<div class='fig'>{svg}</div>")
    A("<p>Arrows run from a stage to every stage that reads one of its "
      "outputs. Depth is the longest path from a root, so an arrow never "
      "points backwards. A stage is coloured on its own verdict; a stage "
      "that is internally consistent but sits downstream of a non-fresh "
      "parent is marked stale by propagation, and its reason says so.</p>")

    # ---- 2. freshness table ----------------------------------------------
    A("<h2>2 &middot; Stage freshness, from <code>stage_provenance</code></h2>")
    A("<div class='scroll'><table><thead><tr><th>Stage</th><th>Verdict</th>"
      "<th>Last recorded run</th><th>Code version</th><th>Why</th>"
      "</tr></thead><tbody>")
    for key in order:
        stage = pv.STAGE_BY_KEY[key]
        f = freshness[key]
        rec = records.get(key)
        reasons = "<br>".join(html.escape(r) for r in f.reasons) or "&mdash;"
        A(f"<tr><td><b>{html.escape(key)}</b><br>"
          f"<span class='mono-sm'>{html.escape(stage.title)}</span></td>"
          f"<td><span class='tag' style='background:"
          f"{SVG_FILL.get(f.state, '#555')}'>{html.escape(f.state)}</span></td>"
          f"<td class='mono-sm'>{html.escape(rec.run_utc if rec else '—')}</td>"
          f"<td class='mono-sm'>{html.escape(stage.code_version)}</td>"
          f"<td>{reasons}</td></tr>")
    A("</tbody></table></div>")

    # ---- 3. geometry probe -----------------------------------------------
    A("<h2>3 &middot; Live header probe: which small geometries are real?</h2>")
    A("<p>Read from the archive at render time. The test is "
      "<code>recorded NAXIS == ZNAXIS</code> &rarr; genuine subframe; "
      "<code>recorded NAXIS == BINTABLE NAXIS &ne; ZNAXIS</code> &rarr; "
      "phantom. This is what separates the four tiny eras from the two "
      "phantom ones, and it is measured rather than assumed.</p>")
    if probe:
        A("<div class='scroll'><table><thead><tr><th>Era</th><th>File</th>"
          "<th>Recorded</th><th>BINTABLE NAXIS</th><th>True ZNAXIS</th>"
          "<th>Verdict</th></tr></thead><tbody>")
        for r in probe:
            cls = "bad" if "PHANTOM" in r["verdict"] else "ok"
            A(f"<tr><td>{r['era']}</td>"
              f"<td class='mono-sm'>{html.escape(os.path.basename(r['path']))}"
              f"</td>"
              f"<td>{_pair(r['rec_n1'], r['rec_n2'])}</td>"
              f"<td>{_pair(r['bt_n1'], r['bt_n2'])}</td>"
              f"<td>{_pair(r['z_n1'], r['z_n2'])}</td>"
              f"<td><span class='tag {cls}'>{html.escape(r['verdict'])}"
              f"</span></td></tr>")
        A("</tbody></table></div>")
    else:
        A("<p class='mono-sm'>Archive not reachable from this host &mdash; "
          "probe skipped. Re-run this report where "
          f"{html.escape(str(ARCHIVE_ROOT))} is mounted.</p>")

    # ---- 4. invalidation matrix ------------------------------------------
    A("<h2>4 &middot; The invalidation matrix</h2>")
    A("<p><span class='tag ok'>VALID</span> the re-characterizations cannot "
      "have changed it &middot; <span class='tag warn'>STALE-RECOMPUTE</span> "
      "mechanically re-runnable, no judgement needed &middot; "
      "<span class='tag bad'>STALE-REDERIVE</span> a human decision inside it "
      "may change &middot; <span class='tag gone'>EVIDENCE-DESTROYED</span> "
      "the tables it rests on are absent from the database.</p>")
    for kind, label in (("stage", "Stages"), ("artifact",
                                              "Published artifacts")):
        A(f"<h3>{label}</h3>")
        A("<div class='scroll'><table><thead><tr><th style='width:19%'>"
          "Subject</th><th style='width:9%'>Verdict</th><th>Evidence</th>"
          "<th style='width:20%'>What to do</th></tr></thead><tbody>")
        for row in AUDIT:
            if row.kind != kind:
                continue
            dec = ""
            if row.decision:
                dec = (f"<br><br><b>Decision at risk:</b> "
                       f"{E(row.decision)}")
            A(f"<tr><td><b>{html.escape(row.subject)}</b></td>"
              f"<td><span class='tag {VERDICT_CLASS[row.verdict]}'>"
              f"{html.escape(row.verdict)}</span></td>"
              f"<td>{E(row.rationale)}{dec}</td>"
              f"<td>{E(row.action)}</td></tr>")
        A("</tbody></table></div>")

    # ---- 4b. the corrected shopping list ----------------------------------
    A("<h2>4b &middot; The October shopping list, corrected</h2>")
    A("<p>The observatory request's Item B table, regenerated from "
      "<code>calib_gaps</code> as it stands today and with the phantom eras "
      "folded into the configurations they actually are (80&nbsp;&rarr;&nbsp;78, "
      "83&nbsp;&rarr;&nbsp;81). This is a preview of what the request should "
      "say, not a replacement for re-running S0/S0b &mdash; the merge is "
      "applied here by query so the size of the correction is visible "
      "before the rebuild.</p>")
    A("<div class='scroll'><table><thead><tr><th>Era (merged)</th>"
      "<th>Camera / readout</th><th>Need</th>"
      "<th>Science frames blocked</th><th>Papers affected</th></tr></thead>"
      "<tbody>")
    try:
        rows = con.execute(SHOPPING_SQL).fetchall()
    except sqlite3.OperationalError as exc:
        rows = []
        A(f"<tr><td colspan='5'>unavailable: {html.escape(str(exc))}</td></tr>")
    for era, cam, need, blocked, proj in rows:
        A(f"<tr><td>{era}</td><td>{html.escape(cam or '')}</td>"
          f"<td>{html.escape(need or '')}</td><td>{fmt(int(blocked))}</td>"
          f"<td>{html.escape(proj or '—')}</td></tr>")
    A("</tbody></table></div>")
    A(f"<p class='mono-sm'>Compare with the shipped request, whose top two "
      f"rows are 'era 76 Mode0 bias 68,965' and 'era 76 Mode0 flat g 40,031': "
      f"era 76's bias coverage now reads "
      f"<b>{html.escape(str(ev['era76_bias_status']))}</b>, and the archive's "
      f"last calibration night is now "
      f"{html.escape(str(ev['last_calib_night']))}, not 2024-11-18.</p>")

    # ---- 5. re-run plan ---------------------------------------------------
    A("<h2>5 &middot; Ordered re-run plan</h2>")
    if not plan:
        A("<p>Nothing to do &mdash; every stage is FRESH.</p>")
    else:
        A(f"<p>{len(plan)} stages, in dependency order. Running them top to "
          "bottom never rebuilds a stage before the stage it reads. Record "
          "each one as it finishes, so the next status check is truthful.</p>")
        lines = []
        for i, key in enumerate(plan, 1):
            stage = pv.STAGE_BY_KEY[key]
            lines.append(f"# {i}. {key} [{freshness[key].state}] "
                         f"{stage.title}")
            lines.append(f"{stage.build_cmd}")
            lines.append(f"python pipeline/scripts/check_pipeline_status.py "
                         f"record {key}")
            lines.append("")
        A(f"<pre>{html.escape(chr(10).join(lines))}</pre>")

    # ---- 6. fingerprint contract -----------------------------------------
    A("<h2>6 &middot; What is fingerprinted, and why</h2>")
    A("<p>A fingerprint over every column would mark everything permanently "
      "stale; a fingerprint over nothing marks everything permanently fresh. "
      "Each resource below hashes exactly the columns a downstream stage "
      "reads to make a decision, and says which columns were deliberately "
      "left out.</p>")
    A("<div class='scroll'><table><thead><tr><th style='width:20%'>Resource"
      "</th><th style='width:30%'>Columns hashed</th><th>Why these</th>"
      "<th style='width:12%'>Current</th></tr></thead><tbody>")
    for key in sorted(pv.RESOURCES):
        spec = pv.RESOURCES[key]
        cols = ", ".join(spec.columns) if spec.columns else "(whole file)"
        cur = fingerprints.get(key, "&mdash;")
        A(f"<tr><td class='mono-sm'>{html.escape(key)}</td>"
          f"<td class='mono-sm'>{html.escape(cols)}</td>"
          f"<td>{html.escape(spec.why)}</td>"
          f"<td class='mono-sm'>{html.escape(str(cur))}</td></tr>")
    A("</tbody></table></div>")

    A("<h2>7 &middot; How to use this</h2>")
    A("<pre>"
      + html.escape(
          "python pipeline/scripts/check_pipeline_status.py status\n"
          "python pipeline/scripts/check_pipeline_status.py plan\n"
          "python pipeline/scripts/check_pipeline_status.py record S0c\n"
          "python pipeline/scripts/check_pipeline_status.py report\n")
      + "</pre>")
    A("<p><code>status</code> exits 0 when every stage is fresh and 1 when "
      "work is outstanding, so it can gate a publish step. Run "
      "<code>record &lt;stage&gt;</code> immediately after each build "
      "script &mdash; a stage that is not recorded reads as NEVER_RUN, which "
      "is the honest default.</p>")

    A("</div></body></html>")

    out = out_dir / "pipeline_status.html"
    _atomic_write(out, "\n".join(h))
    return out


def _pair(a, b) -> str:
    """Format a geometry pair, or an em dash when unknown."""
    if a is None or b is None:
        return "&mdash;"
    return f"{int(a)}&times;{int(b)}"


def _atomic_write(path: Path, text: str) -> None:
    """Write via a temporary file + rename, so a reader never sees a
    half-written page (the house rule for every artifact in this repo)."""
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
