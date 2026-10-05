#!/usr/bin/env python
"""be_step0 — Step 0 of the RE-DRAWN Be-star paper: sample, master frame table,
dispositions, header re-scrape audit, the instrument (mechanical-state) table,
and the per-night sampling table that decides the short-tier admission.

WHY IT EXISTS
-------------
The novelty gate (BE-N1-gate) failed for the ten stars of ANALYSIS_STRATEGY §3.2
and passed for the campaign as observed; the chair re-drew the paper to the
BeSS-verified active Be stars with >= 10 RLMT grism nights in a season
(strategy §10, *Sample re-draw*).  The S0c staging table ``stage_bestar_grism``
holds the old ten, so Step 0 is re-issued here for the new sample, reading the
shared manifest strictly READ-ONLY and writing only this project's database
``BeStar_Grism/products/bestar.sqlite``.

TABLES WRITTEN (all regenerable; every number in the notes is emitted from them)
------------------------------------------------------------------------------
be_sample        one row per star: role (science / standard / null), SIMBAD id,
                 spectral type, the season(s) that verified it, BeSS EW median.
be_frames        one row per canonical raw grism Light frame of those stars:
                 night, season, BJD_TDB, recomputed airmass, mechanical epoch
                 and ruled mechanical state, gain configuration, header cards,
                 S2c verdict and the DISPOSITION that decides extraction.
be_header_audit  per (camera, card): frames, NULLs in the manifest, and how
                 many of those NULL frames carry the card in the FITS header
                 when it is re-opened (must be zero for BE-S0-header-rescrape).
be_era_table     one row per mechanical epoch (F-3) touched by the sample,
                 grouped into the five ruled mechanical states, with the
                 measured gain from ``detector_params`` (or "unmeasured").
be_nights        per (star, night, grism): frames, BJD span — the sampling.
be_short_tier    per star: long nights (> 2 h span) and short-tier admission.
be_meta          build provenance.

RULES APPLIED (each is a ruling, cited)
--------------------------------------
* Sample: taken verbatim from the novelty tables (chair, 2026-10-04).  QQ Gem is
  HD 46264 and is one of the 19; it is not counted twice.
* Standards: eta Hya (HR 3454), tet Vir (HR 4963), Vega.  Null-test stars: every
  star BeSS shows WITHOUT emission in a standards-epoch season with >= 10 of
  our nights (strategy §10; detection rule only from 2025-12-05, SYNTHESIS §4).
* Five mechanical states (SYNTHESIS §4 Era table; TE.F1): Andor; ASI before the
  2025 monsoon; ASI after it (flipped); QHY night 1; QHY with the grisms
  re-seated.  The finer F-3 epochs are kept as rows beneath them, because a
  calibration must not cross ANY F-3 boundary.
* Dispositions (SYNTHESIS §4 "Disposition column"): S2c dispersed -> extract;
  indeterminate -> extract only if the pixel identity gate (G-3) finds the
  target's spectrum, else reject; direct -> excluded (not a spectrum);
  EXPTIME <= 0 -> excluded.

USAGE
-----
    cd "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO"
    /opt/miniconda3/envs/rlmt-checks/bin/python BeStar_Grism/scripts/be_step0.py all
Subcommands: build (tables), rescrape (header audit, opens FITS headers), notes
(Markdown tables), all.
"""
from __future__ import annotations

import argparse
import datetime as dt
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import be_common as C  # noqa: E402

sys.path.insert(0, str(C.REPO / "pipeline"))

SCRIPT = "BeStar_Grism/scripts/be_step0.py"
OUT = C.NOTES / "step0"

#: Header cards of the F-2 re-scrape, as (manifest column, FITS keyword).
CARDS = (("ccd_temp", "CCD-TEMP"), ("set_temp", "SET-TEMP"), ("coolpowr", "COOLPOWR"),
         ("hdr_gain", "GAIN"), ("hdr_offset", "OFFSET"), ("focpos", "FOCPOS"),
         ("flipstat", "FLIPSTAT"), ("telpier", "TELPIER"), ("fwpos", "FWPOS"),
         ("fwallnam", "FWALLNAM"))

#: SIMBAD ids of the standards (strategy §3.2 roles, unchanged by the re-draw).
STANDARDS = ("* eta Hya", "* tet Vir", "* alf Lyr")


def strip_suffix(name: str) -> str:
    """Remove grism/exposure text leaking from filenames into target names
    (identical to bess_novelty_check.strip_suffix, which built the aliases)."""
    return re.sub(r"\s+(hrg|lrg|HaGrism|OGGrism)\b.*$", "", name or "").strip()


# ===========================================================================
# 1. the sample
# ===========================================================================
def load_sample() -> pd.DataFrame:
    """Science, standard and null stars, from the novelty tables (verbatim)."""
    tss = pd.read_csv(C.NOVELTY_TABLES / "target_season_bess.csv")
    nv = C.novelty_ro()
    tg = pd.read_sql("SELECT main_id, label, ra, dec, otype, sptype FROM nv_targets", nv)
    act = tss[(tss.bess_verdict == "EMISSION") & (tss.nights >= C.MIN_NIGHTS_SEASON)]
    # Null-test stars: no emission, a standards-epoch season (2025-26), >= 10 nights.
    null = tss[(tss.bess_verdict == "NO-EMISSION") & (tss.nights >= C.MIN_NIGHTS_SEASON)
               & (tss.season == C.season_of(C.STANDARDS_EPOCH))]
    rows = []
    for mid, g in act.groupby("main_id"):
        rows.append(dict(main_id=mid, role="science",
                         verified_seasons=",".join(sorted(g.season)),
                         bess_ew_med=float(g.ew_med.median()),
                         manifest_names=g.manifest_names.iloc[0]))
    for mid, g in null.groupby("main_id"):
        if mid in act.main_id.values:
            continue
        rows.append(dict(main_id=mid, role="null", verified_seasons=",".join(sorted(g.season)),
                         bess_ew_med=float(g.ew_med.median()),
                         manifest_names=g.manifest_names.iloc[0]))
    for mid in STANDARDS:
        rows.append(dict(main_id=mid, role="standard", verified_seasons="",
                         bess_ew_med=np.nan,
                         manifest_names=tg.loc[tg.main_id == mid, "label"].iloc[0]))
    s = pd.DataFrame(rows).merge(tg, on="main_id", how="left")
    return s


# ===========================================================================
# 2. the master frame table
# ===========================================================================
def ruled_state(camera: str, mech_epoch: str) -> str:
    """The five mechanical states of the ruling, from the F-3 epoch id."""
    if camera == "iKon":
        return "S1 Andor"
    if camera == "ASI":
        return "S2 ASI pre-monsoon" if mech_epoch < "ASI:2025-07" else "S3 ASI post-monsoon (flipped)"
    if camera == "QHY600":
        return "S4 QHY night 1" if mech_epoch == "QHY600:2026-03-21" else "S5 QHY re-seated"
    return "unknown"


def gain_config(camera, readoutm, xbin, egain, swcreate) -> str:
    """The ``s2_camera_configs`` key a frame's gain is measured under."""
    if camera == "ASI":
        return "ASI Mode0 2x2 e0.780" if egain and abs(egain - 0.78) < 0.01 else (
            "ASI Mode0 2x2" if xbin == 2 else "ASI Mode0 1x1")
    if camera == "QHY600":
        return "QHY600 Fast 2x2" if (readoutm or "") == "Fast" else "QHY600 pyscope 2x2"
    if camera == "iKon":
        return "iKon 1MHz 4x" if (readoutm or "").startswith("1MHz") else "iKon 5MHz"
    return "unknown"


def airmass(ra, dec, jd_utc) -> np.ndarray:
    """Recomputed airmass (astropy AltAz at Winer); the header column is
    corrupt (strategy §3.4).  sec z of the true altitude; NaN below 3 deg."""
    from astropy import units as u
    from astropy.coordinates import AltAz, SkyCoord
    from astropy.time import Time
    from macro_core.timing import winer_location
    c = SkyCoord(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg)
    t = Time(np.asarray(jd_utc), format="jd", scale="utc")
    alt = c.transform_to(AltAz(obstime=t, location=winer_location())).alt.deg
    with np.errstate(invalid="ignore", divide="ignore"):
        # Pickering (2002) airmass, good to the horizon.
        h = alt
        x = 1.0 / np.sin(np.radians(h + 244.0 / (165.0 + 47.0 * np.abs(h) ** 1.1)))
    return np.where(alt > 3.0, x, np.nan)


def cmd_build(_a) -> None:
    sample = load_sample()
    nv = C.novelty_ro()
    alias = pd.read_sql("SELECT lname, main_id FROM nv_alias", nv)
    man = C.manifest_ro()
    ph = ",".join("?" * len(C.GRISM_FILTERS))
    cols = ", ".join(f"f.{c}" for c, _ in CARDS)
    fr = pd.read_sql(f"""
        SELECT f.obs_rowid, f.path, f.canonical_target, f.target_key, f.filter, f.night,
               f.era_id, f.exptime, f.jd, f.camera, f.readoutm, f.xbinning, f.egain,
               f.swcreate, f.naxis1, f.naxis2, f.pointing_offset_deg, f.focuspos, {cols},
               m.mech_epoch, m.detector_epoch, m.epoch_certain,
               t.jd_utc_mid, t.bjd_tdb,
               d.verdict AS s2c_verdict, d.strength_class AS s2c_strength,
               d.verdict_basis AS s2c_basis
        FROM frames f
        LEFT JOIN frame_mech_epoch m USING (obs_rowid)
        LEFT JOIN frame_times t ON t.path = f.path
        LEFT JOIN frame_dispersion d ON d.obs_rowid = f.obs_rowid
        WHERE f.filter IN ({ph}) AND f.tree = 'rawimage' AND f.is_canonical = 1
          AND f.imagetyp LIKE 'Light%' AND f.error IS NULL""", man, params=C.GRISM_FILTERS)
    build = dict(man.execute("SELECT key, value FROM build_meta").fetchall())
    s2c = dict(man.execute("SELECT key, value FROM s2c_build_meta").fetchall())
    man.close()

    fr["lname"] = fr.canonical_target.map(lambda s: strip_suffix(s).lower())
    fr = fr.merge(alias, on="lname", how="left")
    fr = fr[fr.main_id.isin(sample.main_id)].merge(
        sample[["main_id", "role", "label", "ra", "dec"]], on="main_id", how="left")
    fr["season"] = fr.night.map(C.season_of)
    fr["standards_epoch"] = (fr.night >= C.STANDARDS_EPOCH).astype(int)
    fr["mech_state"] = [ruled_state(c, m or "") for c, m in zip(fr.camera, fr.mech_epoch)]
    fr["gain_config"] = [gain_config(*r) for r in
                         fr[["camera", "readoutm", "xbinning", "egain", "swcreate"]].itertuples(index=False)]
    jd = fr.jd_utc_mid.fillna(fr.jd + fr.exptime / 86400 / 2)
    fr["airmass_calc"] = airmass(fr.ra, fr.dec, jd)

    # --- disposition (SYNTHESIS §4 Disposition) ---------------------------
    disp, why = [], []
    for v, e in zip(fr.s2c_verdict, fr.exptime):
        if not (e and e > 0):
            disp.append("exclude"); why.append("EXPTIME <= 0")
        elif v == "dispersed":
            disp.append("extract"); why.append("S2c measured dispersed")
        elif v == "indeterminate":
            disp.append("extract_if_identity")
            why.append("S2c indeterminate: kept only if the G-3 pixel identity gate finds the target spectrum")
        elif v == "direct":
            disp.append("exclude"); why.append("S2c measured direct image: no spectrum")
        else:
            disp.append("exclude"); why.append("no S2c verdict")
    fr["disposition"], fr["disposition_reason"] = disp, why
    fr = fr.drop(columns=["lname"])

    con = C.be_db()
    sample.to_sql("be_sample", con, if_exists="replace", index=False)
    fr.to_sql("be_frames", con, if_exists="replace", index=False)
    con.execute("CREATE UNIQUE INDEX IF NOT EXISTS be_frames_id ON be_frames(obs_rowid)")

    # --- the per-night sampling and the short-tier admission --------------
    n = (fr[fr.disposition != "exclude"]
         .groupby(["main_id", "label", "role", "night", "season", "filter", "mech_epoch"])
         .agg(n_frames=("obs_rowid", "size"), bjd_first=("bjd_tdb", "min"),
              bjd_last=("bjd_tdb", "max"), airmass_med=("airmass_calc", "median"))
         .reset_index())
    n["span_min"] = (n.bjd_last - n.bjd_first) * 1440.0
    n.to_sql("be_nights", con, if_exists="replace", index=False)
    per_night = (fr[fr.disposition != "exclude"].groupby(["main_id", "label", "role", "night"])
                 .agg(b0=("bjd_tdb", "min"), b1=("bjd_tdb", "max")).reset_index())
    per_night["span_min"] = (per_night.b1 - per_night.b0) * 1440.0
    st = (per_night.groupby(["main_id", "label", "role"])
          .agg(nights=("night", "size"),
               nights_std_epoch=("night", lambda x: int((x >= C.STANDARDS_EPOCH).sum())),
               long_nights=("span_min", lambda x: int((x > C.LONG_NIGHT_MIN).sum())),
               max_span_min=("span_min", "max"), median_span_min=("span_min", "median"))
          .reset_index())
    st["short_tier"] = (st.long_nights >= C.MIN_LONG_NIGHTS).astype(int)
    st.to_sql("be_short_tier", con, if_exists="replace", index=False)

    # --- the instrument table ---------------------------------------------
    man = C.manifest_ro()
    me = pd.read_sql("SELECT * FROM mech_epoch", man)
    cfg = pd.read_sql("SELECT config, gain_e_per_adu, gain_err, gain_basis, read_noise_e, "
                      "clip_adu, bias_adu FROM s2_camera_configs", man)
    man.close()
    use = (fr.groupby(["mech_epoch", "mech_state", "gain_config"])
           .agg(n_frames=("obs_rowid", "size"), n_nights=("night", "nunique"),
                first_night_be=("night", "min"), last_night_be=("night", "max"),
                grisms=("filter", lambda x: ",".join(sorted(set(x)))))
           .reset_index())
    era = use.merge(me[["mech_epoch", "camera", "first_night", "last_night", "boundary_cause",
                        "rotation_deg", "flipstat", "wheel_map"]], on="mech_epoch", how="left")
    era = era.merge(cfg, left_on="gain_config", right_on="config", how="left").drop(columns="config")
    # A ruled state with no Be frame (QHY night 1) still gets its row, so the
    # table has one row per mechanical state of the ruling.
    for st_name, ep in (("S4 QHY night 1", "QHY600:2026-03-21"),):
        if st_name not in set(era.mech_state):
            m = me[me.mech_epoch == ep].iloc[0]
            era = pd.concat([era, pd.DataFrame([dict(
                mech_epoch=ep, mech_state=st_name, gain_config="QHY600 Fast 2x2", n_frames=0,
                n_nights=0, grisms="", camera=m.camera, first_night=m.first_night,
                last_night=m.last_night, boundary_cause=m.boundary_cause, rotation_deg=m.rotation_deg,
                flipstat=m.flipstat, wheel_map=m.wheel_map)]).merge(
                    cfg, left_on="gain_config", right_on="config", how="left").drop(columns="config")],
                ignore_index=True)
    era["gain_status"] = np.where(era.gain_e_per_adu.notna(), "measured", "UNMEASURED")
    era = era.sort_values("first_night")
    era.to_sql("be_era_table", con, if_exists="replace", index=False)

    # --- the no-crossing test: each (grism, mech_epoch) calibration key maps to
    # exactly one F-3 epoch AND one ruled state, and every frame has an epoch.
    no_epoch = int(fr.mech_epoch.isna().sum())
    cross = int(fr.groupby(["filter", "mech_epoch"]).mech_state.nunique().gt(1).sum())
    span_ok = int(((era.first_night_be < era.first_night) | (era.last_night_be > era.last_night)).sum())

    meta = {
        "built_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": subprocess.run(["git", "-C", str(C.REPO), "rev-parse", "--short", "HEAD"],
                                     capture_output=True, text=True).stdout.strip(),
        "manifest_build": f"{build.get('code_version')} {build.get('built_utc')}",
        "s2c_build": s2c.get("code_version", ""),
        "n_frames": len(fr), "n_frames_science": int((fr.role == "science").sum()),
        "frames_without_mech_epoch": no_epoch,
        "calib_keys_crossing_a_state": cross,
        "epochs_where_be_frames_fall_outside_epoch_bounds": span_ok,
    }
    con.execute("DROP TABLE IF EXISTS be_meta")
    con.execute("CREATE TABLE be_meta (key TEXT PRIMARY KEY, value TEXT)")
    con.executemany("INSERT INTO be_meta VALUES (?,?)", [(k, str(v)) for k, v in meta.items()])
    con.commit()
    print(meta)


# ===========================================================================
# 3. header re-scrape audit
# ===========================================================================
def cmd_rescrape(_a) -> None:
    """For every manifest NULL in the F-2 card columns, re-open the FITS header
    (all HDUs — the era-C files keep the image in an extension) and record
    whether the card is really absent.  A NULL where the card exists is a
    defect; the value read is written into be_frames so nothing downstream
    depends on the defect."""
    from astropy.io import fits
    con = C.be_db()
    fr = pd.read_sql("SELECT * FROM be_frames", con)
    rows, fills = [], []
    for col, key in CARDS:
        for cam, g in fr.groupby("camera"):
            nulls = g[g[col].isna()]
            present = 0
            for r in nulls.itertuples():
                try:
                    with fits.open(C.ARCHIVE_ROOT / r.path, memmap=False) as hl:
                        val = next((h.header[key] for h in hl if key in h.header), None)
                except Exception as exc:              # unreadable: count it, never hide it
                    val = f"ERROR {exc}"
                if val is not None:
                    present += 1
                    fills.append((str(val), r.obs_rowid, col))
            rows.append(dict(camera=cam, card=key, column=col, n_frames=len(g),
                             n_null_manifest=len(nulls), n_manifest_null_but_card_present=present,
                             n_filled_from_header=present))
    audit = pd.DataFrame(rows)
    audit.to_sql("be_header_audit", con, if_exists="replace", index=False)
    for val, oid, col in fills:
        con.execute(f"UPDATE be_frames SET {col}=? WHERE obs_rowid=?", (val, oid))
    con.commit()
    # After the fill, the residual: NULLs left in be_frames whose card exists.
    audit["n_null_card_present_after"] = [
        int(con.execute(f"SELECT COUNT(*) FROM be_frames WHERE camera=? AND {c} IS NULL",
                        (cam,)).fetchone()[0]) - (n - f)
        for cam, c, n, f in zip(audit.camera, audit.column, audit.n_null_manifest,
                                audit.n_manifest_null_but_card_present)]
    audit.to_sql("be_header_audit", con, if_exists="replace", index=False)
    con.commit()
    print(audit[audit.n_null_manifest > 0].to_string())


# ===========================================================================
# 4. emitted notes
# ===========================================================================
def cmd_notes(_a) -> None:
    con = C.be_db()
    meta = dict(con.execute("SELECT key, value FROM be_meta").fetchall())
    foot = (f"\nBuilt {meta['built_utc']} at commit {meta['git_commit']} from manifest "
            f"{meta['manifest_build']} and S2c {meta['s2c_build']} (read-only).\n")

    s = pd.read_sql("""
        SELECT s.role, s.label AS star, s.main_id AS simbad, s.sptype, s.verified_seasons,
               s.bess_ew_med AS bess_ew_med_A,
               COUNT(f.obs_rowid) AS frames,
               COUNT(DISTINCT f.night) AS nights,
               SUM(f.standards_epoch) AS frames_std_epoch,
               SUM(f.filter='hrg') AS hrg, SUM(f.filter='lrg') AS lrg,
               SUM(f.filter IN ('HaGrism','OGGrism')) AS Ha_OG,
               SUM(f.disposition='extract') AS extract,
               SUM(f.disposition='extract_if_identity') AS if_identity,
               SUM(f.disposition='exclude') AS excluded
        FROM be_sample s LEFT JOIN be_frames f USING (main_id)
        GROUP BY s.main_id ORDER BY s.role DESC, nights DESC""", con)
    C.write_md(OUT / "sample.md", "Re-drawn sample (chair's ruling on BE-N1-gate) — frames and dispositions",
               C.md_table(s) + "\nSample rule and list: `notes/novelty/tables/be_n1_count.md`. QQ Gem = HD 46264 "
               "(science row) — one of the 19, not an addition.\n" + foot, SCRIPT)

    d = pd.read_sql("""SELECT role, s2c_verdict, disposition, disposition_reason, COUNT(*) AS frames,
                       COUNT(DISTINCT main_id) AS stars FROM be_frames
                       GROUP BY 1,2,3,4 ORDER BY 1,2""", con)
    qq = pd.read_sql("""SELECT label AS star, main_id AS simbad, role, COUNT(*) AS frames,
                        COUNT(DISTINCT night) AS nights, MIN(night) AS first, MAX(night) AS last,
                        SUM(disposition='extract') AS extract FROM be_frames
                        WHERE main_id='HD  46264'""", con)
    C.write_md(OUT / "dispositions.md", "Dispositions — every frame of the re-drawn set, and QQ Gem",
               C.md_table(d) + "\n**QQ Gem disposition: SAMPLE MEMBER** (verified active in both seasons, "
               "novelty report §4).\n\n" + C.md_table(qq) + foot, SCRIPT)

    e = pd.read_sql("""SELECT mech_state, mech_epoch, boundary_cause, first_night, last_night,
                       first_night_be, last_night_be, n_nights, n_frames, grisms, gain_config,
                       gain_e_per_adu, gain_err, gain_status, read_noise_e, clip_adu
                       FROM be_era_table ORDER BY first_night""", con)
    st = pd.read_sql("""SELECT mech_state, MIN(first_night) AS first_night, MAX(last_night) AS last_night,
                        SUM(n_frames) AS be_frames, GROUP_CONCAT(DISTINCT gain_config) AS gain_configs
                        FROM be_era_table GROUP BY mech_state ORDER BY mech_state""", con)
    C.write_md(OUT / "era_table.md", "Instrument table — the five mechanical states and their F-3 epochs",
               C.md_table(st) + "\n" + C.md_table(e, 3) +
               f"\nNo-crossing test: frames without a mechanical epoch = {meta['frames_without_mech_epoch']}; "
               f"(grism, epoch) calibration keys spanning two states = {meta['calib_keys_crossing_a_state']}; "
               f"epochs whose Be frames fall outside the epoch bounds = "
               f"{meta['epochs_where_be_frames_fall_outside_epoch_bounds']}. Gains are `s2_camera_configs` "
               "measurements (F-4); rows marked UNMEASURED carry no gain and are not used for a noise model.\n"
               + foot, SCRIPT)

    if con.execute("SELECT name FROM sqlite_master WHERE name='be_header_audit'").fetchone():
        h = pd.read_sql("SELECT * FROM be_header_audit ORDER BY camera, card", con)
        C.write_md(OUT / "header_audit.md", "Header re-scrape audit (F-2 cards, re-opened where NULL)",
                   C.md_table(h) + "\nAcceptance: `n_null_card_present_after` = 0 everywhere — a NULL is left in "
                   "`be_frames` only where the card is absent from every HDU of the file. Rows with "
                   "`n_manifest_null_but_card_present` > 0 are manifest (F-2) defects, filled here from the header "
                   "and reported to the S0 owner.\n" + foot, SCRIPT)

    t = pd.read_sql("""SELECT role, label AS star, nights, nights_std_epoch, long_nights,
                       max_span_min, median_span_min, short_tier FROM be_short_tier
                       ORDER BY role DESC, short_tier DESC, long_nights DESC""", con)
    C.write_md(OUT / "short_tier.md", "Short-tier admission (≥ 3 nights of > 2 h span; strategy §3.3)",
               C.md_table(t, 1) + foot, SCRIPT)
    print("notes written to", OUT)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cmd", choices=["build", "rescrape", "notes", "all"])
    a = ap.parse_args(argv)
    if a.cmd in ("build", "all"):
        cmd_build(a)
    if a.cmd in ("rescrape", "all"):
        cmd_rescrape(a)
    if a.cmd in ("notes", "all"):
        cmd_notes(a)


if __name__ == "__main__":
    main()
