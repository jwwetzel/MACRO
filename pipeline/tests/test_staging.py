"""Unit tests for macro_core.staging — every pure S0c function and datum.

Each test class mirrors one rule the staging encodes.  The must-NOT cases
(a reduced-tree row staging, a grism frame entering a photometric set, a
duplicate copy staging beside its canonical twin) matter as much as the
must cases: the staging manifest is what every downstream stage trusts, so
a silent over- or under-selection poisons every project at once.

The build script's DataFrame path is also under test (the S0 lesson: bugs
ship in the wiring, not just the pure layer).

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import sys
from pathlib import Path

import pytest

# Make the package importable regardless of pytest's working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import pandas as pd

from macro_core import staging as stg
from macro_core.manifest import STRATEGY_CLAIMS

import build_s0c_staging as build


# ---------------------------------------------------------------------------
# The selections as data: complete, consistent, reviewable
# ---------------------------------------------------------------------------
class TestSelectionsData:
    def test_five_projects_with_unique_names_and_tables(self):
        names = [s.project for s in stg.PROJECT_SELECTIONS]
        assert len(names) == 5
        assert len(set(names)) == 5
        tables = [stg.stage_table_name(n) for n in names]
        assert len(set(tables)) == 5

    def test_table_names_are_bare_sql_identifiers(self):
        # Interpolated into DDL without quoting — must never need quoting.
        for s in stg.PROJECT_SELECTIONS:
            t = stg.stage_table_name(s.project)
            assert t.startswith("stage_")
            assert t.replace("_", "").isalnum()
            assert t == t.lower()

    def test_every_strategy_claim_target_is_staged_by_its_project(self):
        """The reconciliation guarantee: every target a strategy document
        claims (STRATEGY_CLAIMS) is selected by that same project's staging
        rule — the __dw_survey__ sentinel via the explicit DW_FIELDS list."""
        for project, tkey, *_ in STRATEGY_CLAIMS:
            sel = stg.selection_for(project)
            if tkey == "__dw_survey__":
                assert set(stg.DW_FIELDS) <= set(sel.targets)
            else:
                assert tkey in sel.targets, (project, tkey)

    def test_dw_fields_are_nineteen_verified_keys(self):
        # The strategy's verified count ("Memo 1 said 21 fields; verified
        # 19") — and every key is a normalized dwHHMM+DD form.
        assert len(stg.DW_FIELDS) == 19
        assert len(set(stg.DW_FIELDS)) == 19
        assert all(k.startswith("dw1") and "+" in k for k in stg.DW_FIELDS)

    def test_no_target_claimed_by_two_projects(self):
        # A frame may only ever be claimed through one project's science
        # rule — cross-project sharing happens by reading the OTHER
        # project's manifest, never by double-listing.
        seen: dict[str, str] = {}
        for s in stg.PROJECT_SELECTIONS:
            for t in s.targets:
                assert t not in seen, (t, seen.get(t), s.project)
                seen[t] = s.project

    def test_bestar_excludes_tcrb_by_design(self):
        # "T CrB … not this paper" — the BeStar strategy's own ruling.
        assert "tcrb" not in stg.selection_for("BeStar_Grism").targets

    def test_selection_for_unknown_project_raises(self):
        with pytest.raises(KeyError):
            stg.selection_for("Nonexistent_Project")

    def test_cv_blacklist_matches_strategy_rule(self):
        # CV §3: "exclude grisms, `empty`, `W`, `6`".
        assert stg.CV_EXCLUDED_FILTERS == frozenset(
            {"hrg", "lrg", "HaGrism", "OGGrism", "empty", "W", "6"})

    def test_bestar_whitelist_matches_step0_rule(self):
        assert stg.BESTAR_GRISM_FILTERS == frozenset(
            {"hrg", "lrg", "HaGrism", "OGGrism"})


# ---------------------------------------------------------------------------
# The science predicate, gate by gate
# ---------------------------------------------------------------------------
def _ok(sel, **over):
    """A frame that passes every gate for ``sel``; override to break one."""
    base = dict(target_key=sel.targets[0], imagetyp="Light Frame",
                error=None, is_canonical=1, tree="rawimage",
                filter_name=(next(iter(sel.filter_whitelist))
                             if sel.filter_whitelist else "R"),
                basename="frame_0001.fts.fz")
    base.update(over)
    return stg.is_staged_science(sel, **base)


class TestSciencePredicate:
    sel = stg.selection_for("TCrB_Monitoring")

    def test_accepts_canonical_light_frame_of_listed_target(self):
        assert _ok(self.sel)

    def test_rejects_noncanonical_duplicate_copy(self):
        # THE core law: the duplicate copy never stages beside its twin.
        assert not _ok(self.sel, is_canonical=0)
        assert not _ok(self.sel, is_canonical=None)
        assert not _ok(self.sel, is_canonical=float("nan"))

    def test_rejects_reduced_tree(self):
        # S0b lesson: reduced canonical rows are renamed copies.
        assert not _ok(self.sel, tree="reduced")

    def test_rejects_header_error_frames(self):
        assert not _ok(self.sel, error="header unreadable")
        assert _ok(self.sel, error="   ")   # whitespace-only = no error

    def test_rejects_unlisted_target_and_blank_target(self):
        assert not _ok(self.sel, target_key="vega")
        assert not _ok(self.sel, target_key=None)

    def test_rejects_calibration_frames(self):
        assert not _ok(self.sel, imagetyp="Dark Frame")
        assert not _ok(self.sel, imagetyp="Flat Field")
        # A master flat written with IMAGETYP='Light Frame' (real archive
        # pattern) is caught by its basename.
        assert not _ok(self.sel, basename="master_flat_g_1x1.fts.fz")

    def test_accepts_blank_imagetyp_2026_nights(self):
        # The load-bearing S0b decision: the newest camera wrote no
        # IMAGETYP card; excluding it would hide the current instrument.
        assert _ok(self.sel, imagetyp=None)
        assert _ok(self.sel, imagetyp="   ")

    def test_cv_blacklist_excludes_grism_and_junk_filters(self):
        cv = stg.selection_for("CV_TimeSeries")
        assert _ok(cv, filter_name="g")
        assert _ok(cv, filter_name=None)          # blank passes a blacklist
        for f in ("hrg", "lrg", "HaGrism", "OGGrism", "empty", "W", "6"):
            assert not _ok(cv, filter_name=f), f

    def test_bestar_whitelist_admits_only_grism_filters(self):
        be = stg.selection_for("BeStar_Grism")
        for f in ("hrg", "lrg", "HaGrism", "OGGrism"):
            assert _ok(be, filter_name=f), f
        assert not _ok(be, filter_name="lrgblue")  # logged-and-excluded
        assert not _ok(be, filter_name="V")
        assert not _ok(be, filter_name=None)       # blank fails a whitelist

    def test_dwarf_selection_takes_all_filters(self):
        dw = stg.selection_for("DwarfGalaxy_AGN_Survey")
        for f in ("6", "L", "R", "Ha", None):
            assert _ok(dw, target_key="ngc5548", filter_name=f), f


# ---------------------------------------------------------------------------
# Roles, paths, farm links
# ---------------------------------------------------------------------------
class TestRolesAndPaths:
    def test_role_of_calib_raw_and_master(self):
        assert stg.role_of_calib("bias", 0) == "bias"
        assert stg.role_of_calib("dark", 1) == "master_dark"
        assert stg.role_of_calib("flat", True) == "master_flat"
        assert stg.role_of_calib("flat", None) == "flat"

    def test_role_of_calib_never_guesses(self):
        with pytest.raises(ValueError):
            stg.role_of_calib("fringe", 0)

    def test_abs_path_preserves_spaces_and_layout(self):
        # The real archive root contains spaces — the join must not mangle
        # them (and consumers must quote; the README says so).
        p = stg.abs_archive_path(
            "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive",
            "rawimage/2025-02-21/tcrb_lrg_240s_0.fts.fz")
        assert p == ("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive/"
                     "rawimage/2025-02-21/tcrb_lrg_240s_0.fts.fz")

    def test_farm_link_name_disambiguates_by_night(self):
        # Identically named frames from different nights must not collide
        # inside one role directory.
        a = stg.farm_link_name("2023-05-21", "xwg1411c.fts.fz")
        b = stg.farm_link_name("2023-07-07", "xwg1411c.fts.fz")
        assert a != b
        assert stg.farm_link_name(None, "x.fts").startswith("no-night_")


# ---------------------------------------------------------------------------
# Row builders: one shared schema, provenance-complete
# ---------------------------------------------------------------------------
class TestRowBuilders:
    frame = dict(obs_rowid=42, path="rawimage/2025-02-21/t.fts.fz",
                 tree="rawimage", size=8_800_000, jd=2460727.9,
                 night="2025-02-21", target_key="tcrb",
                 canonical_target="T CrB", filter="lrg", exptime=240.0,
                 era_id=76, dup_group=7, qc_flags="",
                 pointing_offset_deg=0.02)
    calib = dict(obs_rowid=99, path="Calibrations/master_dark_240s.fts",
                 tree="Calibrations", size=36_000_000, jd=None, night=None,
                 era_id=76, exptime=240.0, filter=None, kind="dark",
                 is_master=1)

    def test_science_row_columns_and_content(self):
        sel = stg.selection_for("TCrB_Monitoring")
        row = stg.science_row(sel, self.frame, "/arch", "S0c test")
        assert tuple(row) == stg.STAGE_CSV_COLUMNS
        assert row["role"] == "science"
        assert row["match_basis"] == stg.MATCH_BASIS_SCIENCE
        assert row["abs_path"] == "/arch/rawimage/2025-02-21/t.fts.fz"
        assert row["size_bytes"] == 8_800_000
        assert row["era_id"] == 76
        assert row["stage_build_id"] == "S0c test"

    def test_calib_row_master_role_and_null_target(self):
        row = stg.calib_row(self.calib, "/arch", "S0c test")
        assert tuple(row) == stg.STAGE_CSV_COLUMNS
        assert row["role"] == "master_dark"
        assert row["match_basis"] == stg.MATCH_BASIS_CALIB
        # Calibration frames have no target identity or pointing.
        assert row["canonical_target"] is None
        assert row["target_key"] is None
        assert row["pointing_offset_deg"] is None


# ---------------------------------------------------------------------------
# The build script's DataFrame wiring (tiny synthetic manifest)
# ---------------------------------------------------------------------------
def _frames_df() -> pd.DataFrame:
    """Six frames: 2 stageable T CrB, 1 duplicate copy, 1 reduced,
    1 wrong target, 1 Vega grism frame (era 80)."""
    cols = build.FRAME_COLUMNS
    rows = [
        # canonical T CrB lrg frame, era 76  -> stages
        dict(obs_rowid=1, path="rawimage/n1/a.fts", tree="rawimage",
             basename="a.fts", size=100, jd=2460001.5, night="2023-01-01",
             target_key="tcrb", canonical_target="T CrB",
             imagetyp="Light Frame", filter="lrg", exptime=240.0,
             era_id=76, is_canonical=1, dup_group=1, error=None,
             qc_flags="", pointing_offset_deg=0.1),
        # canonical T CrB imaging frame, era 2  -> stages
        dict(obs_rowid=2, path="rawimage/n2/b.fts", tree="rawimage",
             basename="b.fts", size=200, jd=2460002.5, night="2023-01-02",
             target_key="tcrb", canonical_target="T CrB",
             imagetyp="Light Frame", filter="R", exptime=8.0,
             era_id=2, is_canonical=1, dup_group=2, error=None,
             qc_flags="pointing_gt1deg", pointing_offset_deg=1.4),
        # duplicate copy of frame 1 (mirror tree)  -> must NOT stage
        dict(obs_rowid=3, path="external/n1/a.fts", tree="external",
             basename="a.fts", size=100, jd=2460001.5, night="2023-01-01",
             target_key="tcrb", canonical_target="T CrB",
             imagetyp="Light Frame", filter="lrg", exptime=240.0,
             era_id=76, is_canonical=0, dup_group=1, error=None,
             qc_flags="", pointing_offset_deg=0.1),
        # reduced-tree canonical row  -> must NOT stage
        dict(obs_rowid=4, path="reduced/n1/a_calibrated.fts", tree="reduced",
             basename="a_calibrated.fts", size=90, jd=2460001.6,
             night="2023-01-01", target_key="tcrb",
             canonical_target="T CrB", imagetyp="Light Frame",
             filter="lrg", exptime=240.0, era_id=79, is_canonical=1,
             dup_group=3, error=None, qc_flags="",
             pointing_offset_deg=None),
        # canonical frame of a target no project claims  -> no stage
        dict(obs_rowid=5, path="rawimage/n3/c.fts", tree="rawimage",
             basename="c.fts", size=300, jd=2460003.5, night="2023-01-03",
             target_key="mars", canonical_target="Mars",
             imagetyp="Light Frame", filter="R", exptime=1.0,
             era_id=2, is_canonical=1, dup_group=4, error=None,
             qc_flags="", pointing_offset_deg=None),
        # canonical Vega grism frame, era 80  -> stages for BeStar only
        dict(obs_rowid=6, path="grism/n4/v.fts", tree="grism",
             basename="v.fts", size=400, jd=2460004.5, night="2023-01-04",
             target_key="vega", canonical_target="Vega",
             imagetyp="Light Frame", filter="hrg", exptime=0.1,
             era_id=80, is_canonical=1, dup_group=5, error=None,
             qc_flags="", pointing_offset_deg=0.05),
    ]
    return pd.DataFrame(rows, columns=cols)


def _calib_df() -> pd.DataFrame:
    rows = [
        # era 76: one raw dark + one recovered master dark
        dict(obs_rowid=10, path="calib/d76.fts", tree="calib",
             night="2023-01-01", jd=2460001.2, era_id=76, exptime=240.0,
             filter=None, kind="dark", is_master=0, size=50),
        dict(obs_rowid=11, path="Calibrations/md76.fts", tree="Calibrations",
             night=None, jd=None, era_id=76, exptime=240.0, filter=None,
             kind="dark", is_master=1, size=60),
        # era 2: one bias
        dict(obs_rowid=12, path="calib/b2.fts", tree="calib",
             night="2023-01-02", jd=2460002.2, era_id=2, exptime=0.0,
             filter=None, kind="bias", is_master=0, size=40),
        # era 99: calib for an era NO project's science touches
        dict(obs_rowid=13, path="calib/f99.fts", tree="calib",
             night="2023-01-05", jd=2460005.2, era_id=99, exptime=1.0,
             filter="g", kind="flat", is_master=0, size=30),
    ]
    return pd.DataFrame(rows)


class TestBuildProjectStage:
    def test_tcrb_stage_contents(self):
        sel = stg.selection_for("TCrB_Monitoring")
        df = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                       "/arch", "bid")
        # Science: exactly the two canonical T CrB frames — never the
        # duplicate copy, the reduced row, or the unclaimed target.
        sci = df[df["role"] == "science"]
        assert sorted(sci["obs_rowid"]) == [1, 2]
        # Calibration: eras 76 and 2 attach (raw dark + master + bias);
        # era 99 (untouched) must NOT ride along.
        cal = df[df["role"] != "science"]
        assert sorted(cal["obs_rowid"]) == [10, 11, 12]
        assert set(cal["match_basis"]) == {stg.MATCH_BASIS_CALIB}
        assert (cal[cal["obs_rowid"] == 11]["role"] == "master_dark").all()
        # One schema throughout, in the documented column order.
        assert tuple(df.columns) == stg.STAGE_CSV_COLUMNS
        # QC flags travel with the row — flagged, never dropped.
        flagged = sci[sci["obs_rowid"] == 2]["qc_flags"].iloc[0]
        assert flagged == "pointing_gt1deg"

    def test_science_rows_sort_before_calibration(self):
        sel = stg.selection_for("TCrB_Monitoring")
        df = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                       "/arch", "bid")
        roles = list(df["role"])
        first_cal = roles.index(next(r for r in roles if r != "science"))
        assert all(r == "science" for r in roles[:first_cal])

    def test_bestar_stage_gets_vega_grism_only(self):
        sel = stg.selection_for("BeStar_Grism")
        df = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                       "/arch", "bid")
        sci = df[df["role"] == "science"]
        assert sorted(sci["obs_rowid"]) == [6]
        # Era 80 has no calib rows in the toy census -> zero calib rows,
        # and that absence is visible, not an error.
        assert len(df[df["role"] != "science"]) == 0

    def test_build_is_deterministic(self):
        sel = stg.selection_for("TCrB_Monitoring")
        a = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                      "/arch", "bid")
        b = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                      "/arch", "bid")
        pd.testing.assert_frame_equal(a, b)


class TestBuildSafety:
    def test_stage_tables_never_collide_with_protected_tables(self):
        for s in stg.PROJECT_SELECTIONS:
            assert stg.stage_table_name(s.project) \
                not in build.PROTECTED_TABLES
        assert not (set(build.S0C_FIXED_TABLES) & build.PROTECTED_TABLES)

    def test_readme_states_law_columns_and_regeneration(self):
        sel = stg.selection_for("CV_TimeSeries")
        text = build.readme_text(sel, 100, 50, "bid")
        assert "NO-COPY LAW" in text
        assert "build_s0c_staging.py" in text
        assert "Dropbox does not sync symlink targets" in text
        for col in stg.STAGE_CSV_COLUMNS:
            assert f"`{col}`" in text
        # The checksum decision is stated, not hidden.
        assert "integrity SURROGATE" in stg.CHECKSUM_NOTE
        assert stg.CHECKSUM_NOTE in text
