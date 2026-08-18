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
        # CV §3: "exclude grisms, `empty`, `W`, `6`" — where "grisms" means
        # EVERY spelling the archive uses, including the iKon tree's 'HaG'.
        assert stg.CV_EXCLUDED_FILTERS == frozenset(
            {"hrg", "lrg", "HaGrism", "OGGrism", "HaG", "empty", "W", "6"})

    def test_cv_blacklist_is_derived_from_grism_all(self):
        # THE 'HaG' REGRESSION.  The blacklist enumerated four of the five
        # grism spellings, so ten dispersed iKon frames staged as photometry.
        # Deriving the set is what makes a sixth spelling a one-line fix.
        assert stg.GRISM_ALL <= stg.CV_EXCLUDED_FILTERS
        assert stg.CV_EXCLUDED_FILTERS == \
            stg.GRISM_ALL | stg.CV_NON_PHOTOMETRIC_FILTERS

    def test_grism_all_covers_every_known_archive_spelling(self):
        # The five strings verified present in frames.filter on 2026-08-18.
        # A new acquisition system adding a sixth breaks this test in the
        # live-manifest check below, not silently in a science set.
        assert stg.GRISM_ALL == frozenset(
            {"hrg", "lrg", "HaGrism", "OGGrism", "HaG"})

    def test_bestar_whitelist_matches_step0_rule(self):
        assert stg.BESTAR_GRISM_FILTERS == frozenset(
            {"hrg", "lrg", "HaGrism", "OGGrism"})

    def test_bestar_whitelist_is_a_subset_of_all_grisms(self):
        # BeStar's published Step-0 whitelist names four strings; it may be
        # narrower than GRISM_ALL, but it must never contain a string that
        # is not a grism at all.
        assert stg.BESTAR_GRISM_FILTERS <= stg.GRISM_ALL

    def test_tcrb_excludes_the_H_filter_from_science(self):
        # TCrB §3 rules on this code directly: "'H' (6, presumed Halpha …
        # excluded from science regardless of P0-2 mapping; filter table
        # only)".  Six frames staged as ordinary science before this.
        sel = stg.selection_for("TCrB_Monitoring")
        assert "H" in sel.filter_blacklist
        assert not _ok(sel, filter_name="H")
        assert _ok(sel, filter_name="hrg")     # the series itself survives

    def test_pending_alias_targets_resolve_into_published_targets(self):
        """THE INVARIANT that ties the transitional keys to the S0 fix.

        Every ``pending_alias_targets`` entry must be a key the committed
        SYNONYM_TABLE folds into one of the project's own ``targets``.  So
        the transitional list cannot drift away from the alias fix it
        mirrors, and it cannot smuggle in a target the project never
        published.
        """
        from macro_core.manifest import SYNONYM_TABLE
        for sel in stg.PROJECT_SELECTIONS:
            for key in sel.pending_alias_targets:
                assert key in SYNONYM_TABLE, (sel.project, key)
                assert SYNONYM_TABLE[key] in sel.targets, (sel.project, key)

    def test_sn_stages_the_pinwheel_and_sequence_digit_names(self):
        # The two S0 alias misses, both of which deleted named template
        # material from the SN working set (140 frames and 2 frames).
        sn = stg.selection_for("SN2023ixf_LightCurve")
        for key in ("pinwheelgalaxy", "2023ixf1", "2023ixf2"):
            assert key in sn.all_target_keys, key
            assert _ok(sn, target_key=key), key
        # …and after the S0 rebuild the merged key stages just the same.
        assert _ok(sn, target_key="m101")
        assert _ok(sn, target_key="2023ixf")

    def test_science_roles_and_their_sql_agree(self):
        # The report and the build must classify rows the same way; the
        # fragments are generated from the tuple, never typed twice.
        assert stg.SCIENCE_ROLES == ("science", "science_unresolved")
        assert stg.SQL_SCIENCE_ROLES == \
            "role IN ('science', 'science_unresolved')"
        assert stg.SQL_CALIB_ROLES == \
            "role NOT IN ('science', 'science_unresolved')"

    def test_only_bestar_enables_the_cone_clause(self):
        enabled = {s.project for s in stg.PROJECT_SELECTIONS
                   if s.cone_radius_deg is not None}
        assert enabled == {"BeStar_Grism"}
        assert stg.selection_for("BeStar_Grism").cone_radius_deg == \
            stg.CONE_CANDIDATE_RADIUS_DEG

    def test_every_stage_claim_targets_a_real_project(self):
        names = {s.project for s in stg.PROJECT_SELECTIONS}
        assert {c.project for c in stg.STAGE_CLAIMS} <= names
        for c in stg.STAGE_CLAIMS:
            assert c.claimed_frames > 0
            assert c.claimed_nights is None or c.claimed_nights > 0
            stg.assert_safe_where(c.where)      # raises on anything unsafe

    def test_claim_fragments_reject_statement_injection(self):
        for bad in ("1=1; DROP TABLE frames", "1=1 -- comment",
                    "1=1 AND (SELECT 1) ; delete FROM frames"):
            with pytest.raises(AssertionError):
                stg.assert_safe_where(bad)


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

    def test_cv_rejects_every_grism_spelling(self):
        """THE 'HaG' REGRESSION, at the predicate.

        Ten iKon-tree frames staged as CV photometry because the blacklist
        named four grism strings and the archive uses five.  The paired
        exposure sequence is the proof: '…_hires' (HaG, 120 s) and
        '…_lowres' (OGGrism, 60 s) are the same target, the same night, the
        same dispersed observation — one was excluded, its twin was not.
        """
        cv = stg.selection_for("CV_TimeSeries")
        for f in sorted(stg.GRISM_ALL):
            assert not _ok(cv, target_key="stlmi", filter_name=f,
                           basename="ST LMi-0001_hires.fts.fz"), f
        # The photometric twin from the same sequence still stages.
        assert _ok(cv, target_key="stlmi", filter_name="r",
                   basename="ST LMi-0001_r.fts.fz")


# ---------------------------------------------------------------------------
# The cone clause: name-less frames enter the working set, but not as science
# ---------------------------------------------------------------------------
def _cone_ok(sel, **over):
    """A name-less frame eligible for ``sel``'s cone match; override to break."""
    base = dict(target_key=None, imagetyp="Light Frame", error=None,
                is_canonical=1, tree="rawimage", filter_name="hrg",
                ra_deg=279.23, dec_deg=38.78, basename="IU2_x_hrg_5s_0.fts")
    base.update(over)
    return stg.is_cone_candidate(sel, **base)


class TestConeCandidates:
    be = stg.selection_for("BeStar_Grism")
    cv = stg.selection_for("CV_TimeSeries")

    def test_eligible_nameless_grism_frame(self):
        assert _cone_ok(self.be)

    def test_named_frames_never_take_the_cone_path(self):
        # A frame WITH a name is adjudicated by gate 5, not by coordinates —
        # otherwise a mislabeled frame could be silently re-identified.
        assert not _cone_ok(self.be, target_key="vega")
        assert not _cone_ok(self.be, target_key="tcrb")
        # Whitespace-only names count as name-less (the catalog writes '').
        assert _cone_ok(self.be, target_key="")
        assert _cone_ok(self.be, target_key="   ")

    def test_cone_path_applies_the_same_frame_gates_as_science(self):
        # A name-less frame that would not have been science had it been
        # named must not get in through the coordinate door.
        assert not _cone_ok(self.be, is_canonical=0)
        assert not _cone_ok(self.be, tree="reduced")
        assert not _cone_ok(self.be, error="header unreadable")
        assert not _cone_ok(self.be, imagetyp="Dark Frame")
        assert not _cone_ok(self.be, filter_name="V")     # off the whitelist
        assert not _cone_ok(self.be, basename="master_flat_g.fts")

    def test_cone_needs_usable_coordinates(self):
        assert not _cone_ok(self.be, ra_deg=None)
        assert not _cone_ok(self.be, dec_deg=None)
        assert not _cone_ok(self.be, ra_deg=float("nan"))
        assert not _cone_ok(self.be, dec_deg=float("nan"))

    def test_projects_without_the_clause_never_cone_match(self):
        assert self.cv.cone_radius_deg is None
        assert not _cone_ok(self.cv, filter_name="g")

    def test_cone_match_picks_the_nearest_within_radius(self):
        refs = {"vega": (279.235, 38.784), "spica": (201.298, -11.154)}
        hit = stg.cone_match(279.23, 38.78, refs, 0.25)
        assert hit is not None and hit[0] == "vega"
        assert hit[1] < 0.01

    def test_cone_match_refuses_outside_the_radius(self):
        # The measured confusion case: the 2025-01-23 focus frames sit on
        # Rigel, 1.4 deg from λ Eri, and must NOT be claimed as λ Eri.
        refs = {"lameri": (77.275, -8.736)}
        assert stg.cone_match(78.633, -8.206, refs, 0.25) is None
        assert stg.cone_match(78.633, -8.206, refs, 2.0) is not None

    def test_cone_match_is_deterministic_between_equal_candidates(self):
        # Two references at the same separation: the sorted-key scan makes
        # the winner independent of dict insertion order.
        a = {"aaa": (10.0, 0.0), "zzz": (350.0, 0.0)}
        b = {"zzz": (350.0, 0.0), "aaa": (10.0, 0.0)}
        assert stg.cone_match(0.0, 0.0, a, 20.0) == \
            stg.cone_match(0.0, 0.0, b, 20.0)

    def test_cone_row_is_never_mistakable_for_science(self):
        frame = dict(obs_rowid=7, path="rawimage/n/IU2_Vega_hrg_5s_0.fts.fz",
                     tree="rawimage", size=1000, jd=2460999.5,
                     night="2025-10-20", filter="hrg", exptime=5.0,
                     era_id=80, dup_group=11, qc_flags="")
        row = stg.cone_candidate_row(frame, "vega", 0.0259, "/arch", "bid")
        assert tuple(row) == stg.STAGE_CSV_COLUMNS
        assert row["role"] == stg.ROLE_SCIENCE_UNRESOLVED != stg.ROLE_SCIENCE
        assert row["match_basis"] == stg.MATCH_BASIS_CONE
        # The candidate identity is recorded, the display name is NOT
        # invented, and the separation is carried for Step 0 to rank on.
        assert row["target_key"] == "vega"
        assert row["canonical_target"] is None
        assert row["pointing_offset_deg"] == 0.0259


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
             qc_flags="", pointing_offset_deg=0.05,
             ra_deg=279.235, dec_deg=38.784),
        # two more Vega grism frames so the target reaches the 3-frame
        # minimum for a reference position (the cone can then form)
        dict(obs_rowid=7, path="grism/n4/v2.fts", tree="grism",
             basename="v2.fts", size=400, jd=2460004.6, night="2023-01-04",
             target_key="vega", canonical_target="Vega",
             imagetyp="Light Frame", filter="hrg", exptime=0.1,
             era_id=80, is_canonical=1, dup_group=6, error=None,
             qc_flags="", pointing_offset_deg=0.05,
             ra_deg=279.230, dec_deg=38.780),
        dict(obs_rowid=8, path="grism/n4/v3.fts", tree="grism",
             basename="v3.fts", size=400, jd=2460004.7, night="2023-01-04",
             target_key="vega", canonical_target="Vega",
             imagetyp="Light Frame", filter="lrg", exptime=0.1,
             era_id=80, is_canonical=1, dup_group=7, error=None,
             qc_flags="", pointing_offset_deg=0.05,
             ra_deg=279.240, dec_deg=38.788),
        # NAME-LESS grism frame on Vega's position -> cone candidate, and
        # only for the selection that enables the clause.
        dict(obs_rowid=9, path="rawimage/n5/IU2_Vega_hrg_5s_0.fts.fz",
             tree="rawimage", basename="IU2_Vega_hrg_5s_0.fts.fz", size=400,
             jd=2460005.5, night="2023-01-05", target_key=None,
             canonical_target=None, imagetyp="Light Frame", filter="hrg",
             exptime=5.0, era_id=80, is_canonical=1, dup_group=8,
             error=None, qc_flags="", pointing_offset_deg=None,
             ra_deg=279.236, dec_deg=38.785),
        # NAME-LESS grism frame 1.4 deg away (the Rigel focus-frame case)
        # -> eligible, but outside the cone: must NOT stage.
        dict(obs_rowid=10, path="grism/n5/focus_10115_5s_HighRes_0.fts.fz",
             tree="grism", basename="focus_10115_5s_HighRes_0.fts.fz",
             size=400, jd=2460005.6, night="2023-01-05", target_key=None,
             canonical_target=None, imagetyp="Light Frame", filter="hrg",
             exptime=5.0, era_id=80, is_canonical=1, dup_group=9,
             error=None, qc_flags="", pointing_offset_deg=None,
             ra_deg=280.635, dec_deg=38.785),
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
        assert sorted(sci["obs_rowid"]) == [6, 7, 8]
        # Era 80 has no calib rows in the toy census -> zero calib rows,
        # and that absence is visible, not an error.
        assert len(df[~df["role"].isin(stg.SCIENCE_ROLES)]) == 0

    def test_bestar_cone_stages_the_nameless_frame_on_target(self):
        """The 671-blank-target finding, at the wiring.

        A name-less grism frame pointed at a staged target now enters the
        working set — as ``science_unresolved``, so Step 0 adjudicates it
        INSIDE the manifest instead of querying ``frames`` behind S0c's
        back.  The frame 1.4 deg away (the Rigel focus-frame case) stays
        out: the cone is a pointing test, not a catch-all.
        """
        sel = stg.selection_for("BeStar_Grism")
        df = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                       "/arch", "bid")
        cone = df[df["role"] == stg.ROLE_SCIENCE_UNRESOLVED]
        assert sorted(cone["obs_rowid"]) == [9]
        assert set(cone["match_basis"]) == {stg.MATCH_BASIS_CONE}
        assert set(cone["target_key"]) == {"vega"}
        # It is NOT science, and no science row was invented for it.
        assert 9 not in set(df[df["role"] == "science"]["obs_rowid"])

    def test_projects_without_the_cone_clause_emit_no_cone_rows(self):
        for name in ("TCrB_Monitoring", "CV_TimeSeries",
                     "SN2023ixf_LightCurve", "DwarfGalaxy_AGN_Survey"):
            df = build.build_project_stage(
                stg.selection_for(name), _frames_df(), _calib_df(),
                "/arch", "bid")
            assert (df["role"] == stg.ROLE_SCIENCE_UNRESOLVED).sum() == 0

    def test_cone_rows_sort_between_science_and_calibration(self):
        sel = stg.selection_for("BeStar_Grism")
        df = build.build_project_stage(sel, _frames_df(), _calib_df(),
                                       "/arch", "bid")
        rank = [0 if r == "science" else
                1 if r == stg.ROLE_SCIENCE_UNRESOLVED else 2
                for r in df["role"]]
        assert rank == sorted(rank)

    def test_reference_positions_need_enough_frames(self):
        # Two frames are not a position: a target with fewer than the
        # minimum simply has no cone and claims nothing.
        sci = pd.DataFrame({"target_key": ["vega", "vega", "spica"],
                            "ra_deg": [279.23, 279.24, 201.3],
                            "dec_deg": [38.78, 38.79, -11.15]})
        assert build.reference_positions(sci) == {}
        assert set(build.reference_positions(sci, min_frames=2)) == {"vega"}

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

    def test_readme_regeneration_command_is_runnable_from_anywhere(self):
        """The command a student auditor copies must actually run.

        The first version printed a repo-RELATIVE script path with no
        working-directory instruction, so running it from the directory the
        README lives in gave file-not-found.  The path must be absolute, and
        quoted, because the repo path contains spaces.
        """
        text = build.readme_text(stg.selection_for("TCrB_Monitoring"),
                                 1, 1, "bid")
        script = build.PIPELINE_ROOT / "scripts" / "build_s0c_staging.py"
        assert script.is_absolute() and script.exists()
        assert f'"{script}"' in text          # absolute AND quoted
        # No bare relative invocation survives anywhere in the README.
        assert "\n        pipeline/scripts/" not in text

    def test_readme_documents_the_cone_role_only_where_it_runs(self):
        cone_text = build.readme_text(stg.selection_for("BeStar_Grism"),
                                      1, 1, "bid", 3)
        plain_text = build.readme_text(stg.selection_for("CV_TimeSeries"),
                                       1, 1, "bid", 0)
        assert stg.MATCH_BASIS_CONE in cone_text
        assert "Cone candidates" in cone_text
        assert "Cone candidates" not in plain_text

    def test_git_commit_marks_a_dirty_tree(self, monkeypatch):
        """THE 17ef904 LESSON: a bare HEAD hash lies about what ran."""
        import subprocess as sp

        class _R:
            def __init__(self, out):
                self.stdout = out

        def fake(cmd, **kw):
            return _R("abc1234\n" if "rev-parse" in cmd else " M staging.py\n")
        monkeypatch.setattr(sp, "run", fake)
        assert build._git_commit() == "abc1234-dirty"

        def clean(cmd, **kw):
            return _R("abc1234\n" if "rev-parse" in cmd else "")
        monkeypatch.setattr(sp, "run", clean)
        assert build._git_commit() == "abc1234"

    def test_staging_source_hash_identifies_the_rule_file(self):
        import hashlib
        source = (build.PIPELINE_ROOT / "macro_core" / "staging.py").read_bytes()
        assert build._staging_source_hash() == \
            hashlib.sha256(source).hexdigest()[:12]


# ---------------------------------------------------------------------------
# The live manifest: rules checked against reality, not only against fixtures
# ---------------------------------------------------------------------------
#
# WHY THIS CLASS EXISTS.  Every S0c rule was unit-tested against synthetic
# frames and every unit test passed — while ten dispersed spectra sat in the
# CV photometry set, 142 M101-field frames sat outside the SN working set,
# and six explicitly-excluded 'H' frames sat in the T CrB science set.  Pure
# tests cannot see any of that: they only know the vocabulary the fixtures
# use.  The reviewer found all three by running SQL against the real
# manifest, so the guard belongs here too.
#
# These tests are READ-ONLY (mode=ro URI, generous busy timeout so they
# coexist with a live S1 batch) and skip cleanly when the manifest or its
# stage tables are absent — a fresh clone still runs the suite.

import sqlite3                                                  # noqa: E402

MANIFEST_PATH = (Path(__file__).resolve().parent.parent.parent
                 / "products" / "manifest" / "rlmt-manifest.sqlite")


def _live_con():
    """Read-only connection to the live manifest, or None when absent."""
    if not MANIFEST_PATH.exists():
        return None
    con = sqlite3.connect(f"file:{MANIFEST_PATH}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


def _has_table(con, name: str) -> bool:
    return bool(con.execute(
        "SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?",
        (name,)).fetchone()[0])


@pytest.fixture(scope="module")
def live():
    con = _live_con()
    if con is None:
        pytest.skip(f"live manifest not present: {MANIFEST_PATH}")
    tables = [stg.stage_table_name(s.project) for s in stg.PROJECT_SELECTIONS]
    missing = [t for t in tables if not _has_table(con, t)]
    if missing:
        con.close()
        pytest.skip(f"stage tables not built yet: {missing}")
    yield con
    con.close()


class TestLiveManifest:
    def test_no_dispersed_frame_is_staged_as_photometric_science(self, live):
        """THE 'HaG' REGRESSION, against the archive's real vocabulary.

        The bug was a PARTIAL exclusion: four of five grism spellings
        blacklisted, so the fifth stayed in a photometric science set.  So
        the invariant is not "no project stages grisms" (T CrB and BeStar
        are grism papers) but "a project that excludes any grism spelling
        stages none of them".
        """
        listed = ", ".join(f"'{f}'" for f in sorted(stg.GRISM_ALL))
        checked = 0
        for sel in stg.PROJECT_SELECTIONS:
            if not (stg.GRISM_ALL & sel.filter_blacklist):
                continue           # this project never claimed to exclude them
            checked += 1
            assert stg.GRISM_ALL <= sel.filter_blacklist, \
                f"{sel.project} excludes only part of GRISM_ALL"
            tbl = stg.stage_table_name(sel.project)
            n = live.execute(
                f"SELECT count(*) FROM {tbl} WHERE role = 'science' "
                f'AND "filter" IN ({listed})').fetchone()[0]
            assert n == 0, f"{sel.project} stages {n} dispersed frames"
        assert checked, "no project exercises the grism exclusion any more"

    def test_no_filter_string_escapes_the_grism_vocabulary(self, live):
        """A sixth grism spelling must break a test, not a science set.

        Paired ``_hires``/``_lowres`` basenames mark a dispersed exposure
        sequence.  Every FILTER value on such a frame is either a known
        grism string or one of the two documented non-grism strays (a '6'
        wheel slot on a MaxIm-labeled low-res frame, and an iKon focus
        image with no FILTER card) — enumerated so a NEW stray fails here.
        """
        rows = live.execute(
            'SELECT DISTINCT "filter" FROM frames '
            "WHERE basename LIKE '%\\_hires%' ESCAPE '\\' "
            "   OR basename LIKE '%\\_lowres%' ESCAPE '\\'").fetchall()
        seen = {r[0] for r in rows}
        known_strays = {None, "6"}
        assert seen - known_strays <= stg.GRISM_ALL, \
            f"unknown filter on a dispersed frame: {seen - known_strays}"

    def test_bestar_loses_nothing_by_the_narrower_whitelist(self, live):
        """'HaG' is a grism, but no BeStar target has one — so Step 0's
        published four-string whitelist costs the paper nothing."""
        be = stg.selection_for("BeStar_Grism")
        keys = ", ".join(f"'{k}'" for k in sorted(be.all_target_keys))
        n = live.execute(
            f'SELECT count(*) FROM frames WHERE "filter" = \'HaG\' '
            f"AND is_canonical = 1 AND target_key IN ({keys})").fetchone()[0]
        assert n == 0

    def test_tcrb_H_frames_are_not_science(self, live):
        n = live.execute(
            "SELECT count(*) FROM stage_tcrb_monitoring "
            "WHERE role = 'science' AND \"filter\" = 'H'").fetchone()[0]
        assert n == 0

    def test_sn_stages_the_deep_M101_field_epoch(self, live):
        """The 'pinwheel galaxy' finding, stated so it survives the rebuild.

        Before the alias fix these 140 frames had target_key
        'pinwheelgalaxy'; after it they have 'm101'.  Either way they must
        be IN the SN working set, because they are the deepest post-fade
        M101-field material in the archive and the template plan turns on
        exactly that depth.
        """
        rows = live.execute(
            "SELECT \"filter\", count(*) FROM stage_sn2023ixf_lightcurve "
            "WHERE role = 'science' "
            "AND night IN ('2026-03-21', '2026-03-22') "
            'GROUP BY "filter" ORDER BY "filter"').fetchall()
        assert dict(rows) == {"g": 35, "ha": 35, "i": 35, "r": 35}

    def test_sn_stages_the_2023_11_28_post_fade_template_epoch(self, live):
        """The '2023ixf1'/'2023ixf2' finding, likewise rebuild-proof.

        A trailing sequence digit fused to the object name deleted a NAMED
        template epoch from the working set.  The frames are identified
        here by position and night, not by the key that was the bug.
        """
        rows = live.execute(
            "SELECT count(*) FROM stage_sn2023ixf_lightcurve "
            "WHERE role = 'science' AND night = '2023-11-27' "
            "AND \"filter\" = 'L'").fetchone()
        assert rows[0] == 2

    def test_every_published_inventory_claim_reconciles(self, live):
        """The doc-drift guard.  A strategy number that stops matching the
        working set — because prose drifted OR because the archive grew —
        fails here and on the S0c report, not at referee time."""
        from macro_core import report_s0c
        bad = []
        for claim in stg.STAGE_CLAIMS:
            n_frames, n_nights = report_s0c.measure_claim(live, claim)
            if n_frames != claim.claimed_frames or (
                    claim.claimed_nights is not None
                    and n_nights != claim.claimed_nights):
                bad.append((claim.project, claim.label,
                            claim.claimed_frames, n_frames,
                            claim.claimed_nights, n_nights))
        assert not bad, f"published claims no longer reconcile: {bad}"

    def test_mispointed_ngc5548_night_is_flagged_once_s0_is_rebuilt(self, live):
        """The NGC 5548 pointing finding.

        The fix lives in S0 (``manifest.pointing_reference``), which cannot
        be re-run while the S1 batch holds the manifest.  This test
        activates itself the moment S0 IS rebuilt — the presence of the
        ``pointing_ref_basis`` column is the signal — so the finding cannot
        be quietly forgotten between waves.
        """
        cols = {r[1] for r in live.execute("PRAGMA table_info(frames)")}
        if "pointing_ref_basis" not in cols:
            pytest.skip("S0 not yet rebuilt with the pointing fallback")
        n_off = live.execute(
            "SELECT count(*) FROM frames WHERE target_key = 'ngc5548' "
            "AND pointing_offset_deg IS NOT NULL").fetchone()[0]
        assert n_off > 0, "NGC 5548 still has no pointing offsets"
        n_flagged = live.execute(
            "SELECT count(*) FROM stage_dwarfgalaxy_agn_survey "
            "WHERE role = 'science' AND target_key = 'ngc5548' "
            "AND night = '2023-03-25' "
            "AND qc_flags LIKE '%pointing_gt1deg%'").fetchone()[0]
        assert n_flagged == 10, \
            "the mispointed 2023-03-25 night is still unflagged"
