"""Unit tests for macro_core.inventory — every pure S0b function.

Each test class mirrors one decision the inventory encodes; the must-NOT
cases (a science frame classified as a flat, a stem match jumping nights, a
tolerance bridging two real exposure settings) are as important as the must
cases, because a silent false link or false calibration poisons the October
shopping list.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import sys
from pathlib import Path

import pytest

# Make the package importable regardless of pytest's working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
# The build script's DataFrame paths are also under test (the S0 era_id .map
# bug shipped exactly because only the pure layer was tested).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import pandas as pd

from macro_core import inventory as inv

import build_s0b_inventory as build


# ---------------------------------------------------------------------------
# Basename surgery
# ---------------------------------------------------------------------------
class TestStems:
    def test_compression_and_extension_strip(self):
        assert inv.frame_stem("a_b.fts.fz") == "a_b"
        assert inv.frame_stem("a_b.fts") == "a_b"
        assert inv.frame_stem("a_b.FIT") == "a_b"
        assert inv.frame_stem("a_b.fits.gz") == "a_b"

    def test_no_extension_returned_unchanged(self):
        assert inv.frame_stem("weird_name") == "weird_name"

    def test_calibrated_suffix_stripped(self):
        # The dominant observed rename: raw name + '_calibrated'.
        assert inv.reduced_stem(
            "mlw_V426_Oph_g_5s_2026-06-27T05-40-49_calibrated.fts.fz"
        ) == "mlw_V426_Oph_g_5s_2026-06-27T05-40-49"

    def test_cal_and_wcs_suffixes_stripped(self):
        assert inv.reduced_stem("1070_M13_10s_Ha_0_cal.fts.fz") \
            == "1070_M13_10s_Ha_0"
        assert inv.reduced_stem("knh_ngc3169_green_60s_x_wcs.fts.fz") \
            == "knh_ngc3169_green_60s_x"

    def test_copy_counter_after_suffix_stripped(self):
        # '…_calibrated_1' — counter rides ON the suffix, both go.
        assert inv.reduced_stem(
            "mpg_NGC_7619_g_180s_2026-06-30T10-30-00_calibrated_1.fts"
        ) == "mpg_NGC_7619_g_180s_2026-06-30T10-30-00"

    def test_bare_trailing_number_is_NOT_stripped(self):
        # MUST NOT: a bare frame index is part of the identity
        # (BeStar ladder files end in one).
        assert inv.reduced_stem("Vega_0p1s_hrg_7.fts.fz") == "Vega_0p1s_hrg_7"

    def test_unknown_suffix_is_NOT_stripped(self):
        # '_test2' is not a known processing suffix; guessing is banned.
        assert inv.reduced_stem("x_2026-06-30T10-30-00_test2.fts") \
            == "x_2026-06-30T10-30-00_test2"

    def test_plain_raw_name_passes_through(self):
        assert inv.reduced_stem("jos_5_Cnc_hrg_300s_2025-11-07T09-21-05.fts.fz") \
            == "jos_5_Cnc_hrg_300s_2025-11-07T09-21-05"


# ---------------------------------------------------------------------------
# Calibration-kind normalization
# ---------------------------------------------------------------------------
class TestCalibKind:
    def test_explicit_imagetyp_wins(self):
        assert inv.calib_kind("Bias Frame", "anything.fts") == "bias"
        assert inv.calib_kind("Dark Frame", "anything.fts") == "dark"
        assert inv.calib_kind("Flat Field", "anything.fts") == "flat"
        assert inv.calib_kind("FLAT", "anything.fts") == "flat"

    def test_master_flat_labeled_light_frame(self):
        # 99 real archive files: master flats written as 'Light Frame'.
        assert inv.calib_kind(
            "Light Frame", "master_flat_g_1x1_Readout2_3s.fts.fz") == "flat"

    def test_flatdark_is_a_dark_not_a_flat(self):
        # Order matters: a flat-field DARK contains both words.
        assert inv.calib_kind(
            "Dark Frame", "master_flatdark_1x1_HighGain_0-7s.fts.fz") == "dark"
        assert inv.calib_kind(
            "Light Frame", "master-dark-flat-16-1x1.fts.fz") == "dark"

    def test_ikon_flat_series(self):
        assert inv.calib_kind("Light Frame", "ha.flat3.fts.fz") == "flat"
        assert inv.calib_kind("Light Frame", "OGG.flat2.fts.fz") == "flat"
        assert inv.calib_kind(
            "Light Frame",
            "HRG.flat6_2024-05-15T12-02-02.fts.fz") == "flat"
        assert inv.calib_kind(
            "Light Frame", "g.Flat-light-1MHz.20s.2.fts.fz") == "flat"

    def test_science_frame_is_NOT_calibration(self):
        # MUST NOT: ordinary science names never classify.
        assert inv.calib_kind(
            "Light Frame",
            "mjc_V426_Oph_g_5s_2026-07-01T07-41-29.fts.fz") is None
        assert inv.calib_kind(None,
                              "mpg_M87_r_90s_2026-06-29T05-23-40.fts.fz") is None

    def test_fringe_field_stays_science(self):
        # Fringe frames are sky exposures of a fringe field — science.
        assert inv.calib_kind(
            "Light Frame",
            "irm_Fringe_field_24_z_120s_2024-10-05T09-22-35.fts.fz") is None

    def test_master_flag(self):
        assert inv.is_master("master_flat_g_1x1.fts.fz")
        assert inv.is_master("master-dark-high-8-1x1.fts.fz")
        assert not inv.is_master("ha.flat3.fts.fz")


# ---------------------------------------------------------------------------
# Science selection
# ---------------------------------------------------------------------------
class TestIsScience:
    def test_light_frame_is_science(self):
        assert inv.is_science("Light Frame", None)

    def test_blank_imagetyp_is_science(self):
        # The 2026-06/07 nights lack IMAGETYP entirely — they are the
        # CURRENT camera and must not vanish from the shopping list.
        assert inv.is_science(None, None)
        assert inv.is_science("", None)

    def test_calibration_kind_disqualifies(self):
        # MUST NOT: a master flat labeled 'Light Frame' is not science.
        assert not inv.is_science("Light Frame", "flat")

    def test_explicit_calib_imagetyp_disqualifies(self):
        assert not inv.is_science("Dark Frame", "dark")


# ---------------------------------------------------------------------------
# Exposure-time binning and dark matching
# ---------------------------------------------------------------------------
class TestExptime:
    def test_bins_absorb_driver_float_fuzz(self):
        assert inv.exptime_bin(15.9999628067017) == 16.0
        assert inv.exptime_bin(0.09998016059399) == 0.1
        assert inv.exptime_bin(2.0000159740448) == 2.0
        assert inv.exptime_bin(240.0) == 240.0

    def test_bins_keep_real_settings_apart(self):
        # MUST NOT: adjacent ladder settings never share a bin.
        assert inv.exptime_bin(8.0) != inv.exptime_bin(16.0)
        assert inv.exptime_bin(0.125) != inv.exptime_bin(0.25)

    def test_missing_and_nonpositive(self):
        assert inv.exptime_bin(None) is None
        assert inv.exptime_bin(float("nan")) is None
        assert inv.exptime_bin(-1.0) == 0.0

    def test_dark_matches_within_tolerance(self):
        assert inv.dark_matches(15.9999628067017, 16.0)
        assert inv.dark_matches(240.00003, 240.0)
        # Sub-second: absolute floor carries the match.
        assert inv.dark_matches(0.09998016059399, 0.1)

    def test_dark_does_NOT_match_across_ladder_steps(self):
        # MUST NOT: the tolerance can never bridge two real settings.
        assert not inv.dark_matches(8.0, 16.0)
        assert not inv.dark_matches(120.0, 240.0)
        assert not inv.dark_matches(0.25, 0.5)

    def test_missing_never_matches(self):
        assert not inv.dark_matches(None, 240.0)
        assert not inv.dark_matches(240.0, None)


# ---------------------------------------------------------------------------
# The match ladder
# ---------------------------------------------------------------------------
class TestLinkReduced:
    # A tiny raw-side world: two frames of one target on one night, plus a
    # second visit of the same field on another night.
    RAW_BY_STEM = {
        "a_T_g_10s_T1": [(1, 2460000.10, "2023-05-31")],
        "a_T_g_10s_T2": [(2, 2460000.20, "2023-05-31")],
        "a_T_g_10s_T9": [(9, 2460030.10, "2023-06-30")],
    }
    RAW_BY_TJD = {
        ("t", round(2460000.10, 7)): [(1, 2460000.10, "2023-05-31")],
        ("t", round(2460000.20, 7)): [(2, 2460000.20, "2023-05-31")],
        # A burst second: two raw frames share (target, JD).
        ("b", round(2460001.30, 7)): [(5, 2460001.30, "2023-06-01"),
                                      (6, 2460001.30, "2023-06-01")],
    }

    def link(self, stem, jd, night, tkey):
        return inv.link_reduced(stem, jd, night, tkey,
                                self.RAW_BY_STEM, self.RAW_BY_TJD)

    def test_stem_with_identical_jd(self):
        out = self.link("a_T_g_10s_T1", 2460000.10, "2023-05-31", "t")
        assert out == [(1, "stem_jd", 0.0)]

    def test_stem_with_rewritten_jd_same_night(self):
        # 30 s drift, same night: the observed reduction-pipeline rewrite.
        out = self.link("a_T_g_10s_T1", 2460000.10 - 30 / 86400.0,
                        "2023-05-31", "t")
        assert len(out) == 1
        rid, method, drift = out[0]
        assert (rid, method) == (1, "stem_jd_drift")
        # JD arithmetic near 2.46e6 has ~4e-5 s of float granularity, so the
        # drift is exact only to that level — the tolerance reflects it.
        assert drift == pytest.approx(-30.0, abs=1e-3)

    def test_stem_drift_does_NOT_jump_nights(self):
        # MUST NOT: same stem exists on 2023-06-30 (raw_id 9), but a frame
        # labeled a different night may not drift-match it; with no target
        # either, it must fall off the ladder entirely.
        out = self.link("a_T_g_10s_T9", 2460031.10, "2023-07-01", None)
        assert out == []

    def test_target_jd_when_stem_unknown(self):
        out = self.link("totally_renamed_file", 2460000.20, "2023-05-31", "t")
        assert out == [(2, "target_jd", 0.0)]

    def test_target_jd_ambiguous_records_every_candidate(self):
        out = self.link("renamed", 2460001.30, "2023-06-01", "b")
        assert [(r, m) for r, m, _ in out] == \
            [(5, "target_jd_ambiguous"), (6, "target_jd_ambiguous")]

    def test_orphan_when_nothing_matches(self):
        assert self.link("stack_of_everything", 2460500.5, "2024-10-01",
                         "unknowntarget") == []

    def test_no_jd_no_stem_match_is_orphan(self):
        assert self.link("a_T_g_10s_T1", None, None, None) == []


# ---------------------------------------------------------------------------
# Coverage arithmetic
# ---------------------------------------------------------------------------
class TestCoverageStatus:
    def test_ok_at_spec(self):
        assert inv.coverage_status(15, 0, 15) == "ok"
        assert inv.coverage_status(40, 2, 15) == "ok"

    def test_partial_below_spec(self):
        assert inv.coverage_status(3, 0, 15) == "partial"

    def test_master_only(self):
        assert inv.coverage_status(0, 1, 15) == "master_only"

    def test_missing(self):
        assert inv.coverage_status(0, 0, 15) == "missing"

    def test_gap_spec_strings(self):
        assert inv.gap_spec("dark", "240s", 0, 15) == "dark 240s x >=15 (have 0)"
        assert inv.gap_spec("bias", None, 3, 20) == "bias x >=20 (have 3)"
        assert inv.gap_spec("flat", "g", 2, 10) == "flat g x >=10 (have 2)"

    def test_fmt_exptime(self):
        assert inv.fmt_exptime(240.0) == "240s"
        assert inv.fmt_exptime(0.1) == "0.1s"
        assert inv.fmt_exptime(None) == "?s"


class TestCalibVocabFilter:
    # Regression (adversarial review, 2026-08-17): era 76 holds 3 grism
    # science frames whose FILTER card reads 'dark' — a header glitch that
    # spawned a nonsense 'flat dark x >=10' shopping-list row.
    def test_vocab_collisions_detected(self):
        assert inv.is_calib_vocab_filter("dark")
        assert inv.is_calib_vocab_filter("Dark")       # case-insensitive
        assert inv.is_calib_vocab_filter(" bias ")     # whitespace-stripped
        assert inv.is_calib_vocab_filter("flat")

    def test_real_filters_are_NOT_collisions(self):
        # MUST NOT: physical filters and grism labels always pass through.
        for f in ("g", "r", "i", "B", "ha", "hrg", "lrg", "lum", "6"):
            assert not inv.is_calib_vocab_filter(f)

    def test_missing_and_blank_are_NOT_collisions(self):
        # A blank filter is a separate honest fact, not a glitch.
        assert not inv.is_calib_vocab_filter(None)
        assert not inv.is_calib_vocab_filter("")
        assert not inv.is_calib_vocab_filter("(blank)")


class TestProjectsOfTarget:
    P = {"tcrb": frozenset({"TCrB_Monitoring"}),
         "stlmi": frozenset({"CV_TimeSeries"})}
    DW = frozenset({"DwarfGalaxy_AGN_Survey"})

    def test_explicit_target(self):
        assert inv.projects_of_target("tcrb", self.P, self.DW) == \
            frozenset({"TCrB_Monitoring"})

    def test_dw_prefix_rule(self):
        assert inv.projects_of_target("dw1403+49", self.P, self.DW) == \
            frozenset({"DwarfGalaxy_AGN_Survey"})

    def test_unknown_target_maps_to_nothing(self):
        assert inv.projects_of_target("m87", self.P, self.DW) == frozenset()
        assert inv.projects_of_target(None, self.P, self.DW) == frozenset()

    def test_prefix_matches_s0_build_selector(self):
        # The S0 build script's __dw_survey__ selector and this module's
        # prefix constant must be the same rule, forever.
        import inspect
        import build_s0_manifest as s0build
        src = inspect.getsource(s0build.build_project_counts)
        assert f'startswith("{inv.DW_SURVEY_PREFIX}")' in src


# ---------------------------------------------------------------------------
# The build script's DataFrame paths (toy archive, end to end)
# ---------------------------------------------------------------------------
def _toy_frames() -> pd.DataFrame:
    """A six-frame toy archive exercising every linkage and coverage path.

    raw canonical science (era 1) .. rowid 1
    its exact reduced copy         .. rowid 2  (same basename, same JD)
    its renamed reduced copy       .. rowid 3  ('_calibrated', same JD)
    reduced orphan                 .. rowid 4
    raw dark, matching exptime     .. rowid 5
    master flat labeled Light      .. rowid 6
    """
    base = dict(canonical_target=None, error=None, readoutm="Fast",
                camtemp=-10.0, ccd_temp=-10.0)
    rows = [
        dict(obs_rowid=1, path="rawimage/n1/a_T_g_10s_T1.fts.fz",
             tree="rawimage", basename="a_T_g_10s_T1.fts.fz",
             jd=2460000.10, night="2023-05-31", target_key="t",
             imagetyp="Light Frame", filter="g", exptime=10.0, era_id=1,
             is_canonical=1, dup_group=100, **base),
        dict(obs_rowid=2, path="reduced/n1/a_T_g_10s_T1.fts.fz",
             tree="reduced", basename="a_T_g_10s_T1.fts.fz",
             jd=2460000.10, night="2023-05-31", target_key="t",
             imagetyp="Light Frame", filter="g", exptime=10.0, era_id=1,
             is_canonical=0, dup_group=100, **base),
        dict(obs_rowid=3, path="reduced/n1/a_T_g_10s_T1_calibrated.fts.fz",
             tree="reduced", basename="a_T_g_10s_T1_calibrated.fts.fz",
             jd=2460000.10, night="2023-05-31", target_key="t",
             imagetyp="Light Frame", filter="g", exptime=10.0, era_id=1,
             is_canonical=1, dup_group=101, **base),
        dict(obs_rowid=4, path="reduced/n1/deep_stack_of_T.fts.fz",
             tree="reduced", basename="deep_stack_of_T.fts.fz",
             jd=2460000.90, night="2023-05-31", target_key=None,
             imagetyp="Light Frame", filter="g", exptime=100.0, era_id=1,
             is_canonical=1, dup_group=102, **base),
        dict(obs_rowid=5, path="rawimage/n1/dark_10s_001.fts.fz",
             tree="rawimage", basename="dark_10s_001.fts.fz",
             jd=2460000.60, night="2023-05-31", target_key=None,
             imagetyp="Dark Frame", filter=None, exptime=10.0000159,
             era_id=1, is_canonical=1, dup_group=103, **base),
        dict(obs_rowid=6, path="calib/master_flat_g_1x1.fts.fz",
             tree="calib", basename="master_flat_g_1x1.fts.fz",
             jd=2460000.55, night="2023-05-31", target_key=None,
             imagetyp="Light Frame", filter="g", exptime=3.0, era_id=1,
             is_canonical=1, dup_group=104, **base),
    ]
    df = pd.DataFrame(rows)
    df["calib_kind"] = [inv.calib_kind(it, bn)
                        for it, bn in zip(df["imagetyp"], df["basename"])]
    return df


class TestBuildLinks:
    def test_toy_archive_links_every_population(self):
        links = build.build_links(_toy_frames())
        by_reduced = {r.reduced_rowid: r for r in links.itertuples()}
        # Exact copy links through the dup_group.
        assert by_reduced[2].match_method == "same_basename_jd"
        assert by_reduced[2].raw_rowid == 1
        # Renamed copy links through the stem ladder.
        assert by_reduced[3].match_method == "stem_jd"
        assert by_reduced[3].raw_rowid == 1
        assert by_reduced[3].raw_path == "rawimage/n1/a_T_g_10s_T1.fts.fz"
        # The stack is an orphan with NULL raw columns — recorded, not lost.
        assert by_reduced[4].match_method == "orphan"
        assert pd.isna(by_reduced[4].raw_rowid)
        # Nothing else invented: exactly the three reduced rows appear.
        assert sorted(by_reduced) == [2, 3, 4]

    def test_agrees_with_pure_ladder_member_by_member(self):
        # The vectorized script path must agree with the pure function —
        # the S0 lesson (era_id .map bug) applied to S0b.
        df = _toy_frames()
        links = build.build_links(df)
        raw = df[(df.tree != "reduced") & (df.is_canonical == 1)]
        raw_by_stem, raw_by_tjd = {}, {}
        for rid, bn, jd, night, tk in zip(raw.obs_rowid, raw.basename,
                                          raw.jd, raw.night, raw.target_key):
            raw_by_stem.setdefault(inv.frame_stem(bn), []).append(
                (int(rid), float(jd), night))
            if tk is not None:
                raw_by_tjd.setdefault((tk, round(float(jd), 7)), []).append(
                    (int(rid), float(jd), night))
        row3 = df[df.obs_rowid == 3].iloc[0]
        pure = inv.link_reduced(inv.reduced_stem(row3.basename),
                                float(row3.jd), row3.night, row3.target_key,
                                raw_by_stem, raw_by_tjd)
        script = links[links.reduced_rowid == 3].iloc[0]
        assert pure == [(script.raw_rowid, script.match_method,
                         script.jd_drift_s)]


class TestBuildCalibAndCoverage:
    def test_calib_frames_kinds_and_masters(self):
        calib = build.build_calib_frames(_toy_frames())
        kinds = dict(zip(calib["obs_rowid"], calib["kind"]))
        assert kinds == {5: "dark", 6: "flat"}
        masters = dict(zip(calib["obs_rowid"], calib["is_master"]))
        assert masters == {5: 0, 6: 1}

    def test_science_selection_excludes_reduced_and_calib(self):
        sci = build.select_science(_toy_frames())
        # Only the raw canonical science frame: the reduced copies and both
        # calibration frames are out.
        assert list(sci["obs_rowid"]) == [1]

    def test_coverage_and_gaps_on_the_toy_era(self):
        cov, gaps = build.build_coverage(
            _toy_frames(),
            {"t": frozenset({"ToyProject"})}, frozenset())
        cell = {(r.req_kind, r.req_key): r for r in cov.itertuples()}
        # Bias: none at all -> missing, and a gap row naming the project.
        assert cell[("bias", None)].status == "missing"
        # Dark at 10s: one raw dark (fuzzy exptime) matches, below spec.
        dark = cell[("dark", "10s")]
        assert (dark.n_calib_raw, dark.status) == (1, "partial")
        # No bias in the era -> the scaled-dark note must be ABSENT (None).
        assert dark.scaled_dark_ok is None or pd.isna(dark.scaled_dark_ok)
        # Flat g: master only — usable but not to spec.
        flat = cell[("flat", "g")]
        assert (flat.n_calib_raw, flat.n_calib_master,
                flat.status) == (0, 1, "master_only")
        # Every gap row names the blocked project via project_counts logic.
        assert set(gaps["projects_affected"]) == {"ToyProject"}
        assert set(gaps["need_kind"]) == {"bias", "dark", "flat"}
        # Ranking: descending by blocked science frames.
        blocked = list(gaps["n_science_frames_blocked"])
        assert blocked == sorted(blocked, reverse=True)

    def test_glitch_filter_stays_in_matrix_but_off_the_shopping_list(self):
        # Regression (adversarial review, 2026-08-17): a science frame whose
        # FILTER header collides with the calibration vocabulary ('dark')
        # keeps its coverage cell — the matrix hides nothing — but must NOT
        # emit a 'flat dark x >=N' shopping-list row: that is not an
        # acquirable item and could mislead the October ops request.
        df = _toy_frames()
        glitch = df[df["obs_rowid"] == 1].copy()
        glitch["obs_rowid"] = 7
        glitch["path"] = "rawimage/n1/mjc_HD_6343_hrg_83s.fts.fz"
        glitch["basename"] = "mjc_HD_6343_hrg_83s.fts.fz"
        glitch["filter"] = "dark"          # the header glitch under test
        glitch["dup_group"] = 105
        df = pd.concat([df, glitch], ignore_index=True)
        cov, gaps = build.build_coverage(
            df, {"t": frozenset({"ToyProject"})}, frozenset())
        cell = {(r.req_kind, r.req_key): r for r in cov.itertuples()}
        # The matrix cell exists and tells the truth (no flats for it).
        assert cell[("flat", "dark")].status == "missing"
        # ...but no shopping-list row asks anyone to acquire it.
        assert not any(s.startswith("flat dark ") for s in gaps["spec"])
        # MUST-still: real filters keep gapping exactly as before.
        assert any(s.startswith("flat g ") for s in gaps["spec"])
        assert any(s.startswith("bias ") for s in gaps["spec"])
