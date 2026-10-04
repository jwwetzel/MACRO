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


# ===========================================================================
# Mechanical epochs (finding F-3 / TE.F1, plan review 2026-10-03)
# ===========================================================================
# Every fixture below is a hand-built list of NightState objects, so the
# exact evidence each rule sees is on the page.  The rotation values are
# modelled on the real archive: the ASI camera sat at +0.45 deg before the
# 2025 monsoon and at 179.77 deg after it; a field at dec -19 reads +0.25
# deg higher than one at dec +25 with nothing touched (the pointing term).

def _rot(target, rot, ra=150.0, dec=25.0, n=5):
    """``n`` identical plate solves of one target."""
    return tuple(inv.RotSample(target, rot, ra, dec) for _ in range(n))


def _night(night, rot=(), cam="ASI", flip="Flip/Mirror", wheel="g|r|i",
           maps=()):
    return inv.NightState(camera=cam, night=night, n_frames=10,
                          flipstat=flip, wheel_map=wheel, rot=tuple(rot),
                          wheel_maps=tuple(maps))


def _causes(epochs):
    return [e["boundary_cause"] for e in epochs]


class TestCircularStatistics:
    def test_circ_diff_takes_the_short_way_round(self):
        assert inv.circ_diff(0.1, 359.9) == pytest.approx(0.2)
        assert inv.circ_diff(359.9, 0.1) == pytest.approx(-0.2)
        assert inv.circ_diff(179.75, 0.43) == pytest.approx(179.32)

    def test_median_across_the_seam(self):
        # A plain median of these is 0.2 by luck of ordering; add 359.x
        # values and it would be ~180.  The circular one is 0.0.
        med, mad = inv.circular_median([359.8, 359.9, 0.0, 0.1, 0.2])
        assert inv.circ_diff(med, 0.0) == pytest.approx(0.0, abs=1e-9)
        assert mad == pytest.approx(0.1)

    def test_empty_input_raises(self):
        with pytest.raises(ValueError):
            inv.circular_median([])


class TestWheelCompareKey:
    def test_software_respellings_are_one_wheel(self):
        # MaxIm's labels vs pyscope's, same seven slots (January 2025).
        maxim = "g|OGGrism|r|i|HaGrism|Ha|Dark"
        pyscope = "g|lrg|r|i|hrg|ha|dark"
        assert inv.wheel_compare_key(maxim) == inv.wheel_compare_key(pyscope)

    def test_a_swapped_slot_is_a_different_wheel(self):
        assert inv.wheel_compare_key("lum|red|green|blue|empty") != \
            inv.wheel_compare_key("lum|red|green|blue|g")

    def test_missing_map_has_no_key(self):
        assert inv.wheel_compare_key(None) is None
        assert inv.wheel_compare_key("  ") is None


class TestSegmentMechEpochs:
    def test_first_night_opens_an_epoch_named_for_camera_and_night(self):
        epochs, nights, _ = inv.segment_mech_epochs(
            [_night("2024-12-14", _rot("a", 0.45))])
        assert len(epochs) == 1
        assert epochs[0]["mech_epoch"] == "ASI:2024-12-14"
        assert epochs[0]["boundary_cause"] == inv.CAUSE_CAMERA
        assert nights[0]["mech_epoch"] == "ASI:2024-12-14"

    def test_a_quiet_season_is_one_epoch(self):
        states = [_night(f"2025-02-{d:02d}", _rot("stlmi", 0.45 + 0.01 * d))
                  for d in range(1, 9)]
        epochs, _, null = inv.segment_mech_epochs(states)
        assert len(epochs) == 1
        assert epochs[0]["first_night"] == "2025-02-01"
        assert epochs[0]["last_night"] == "2025-02-08"
        assert epochs[0]["n_transitions"] == 7
        assert epochs[0]["n_tested_fine"] == 7
        # Every within-epoch same-target difference lands in the null.
        assert len(null) == 7
        assert all(abs(x["delta_deg"]) < inv.MECH_ROT_STEP_DEG for x in null)

    def test_the_pointing_trap_does_not_fire(self):
        """THE CASE THE NAIVE RULE GETS WRONG.  Nothing was touched; the
        target list changed from a dec +25 field to a dec -19 field and the
        nightly median moved by 0.25 deg, then to a dec +56 field (-0.44).
        A 'nightly median moved > 0.3 deg' rule would cut here.  With no
        shared target the fine test cannot run, and 0.44 is far below the
        coarse threshold — so: one epoch, and the transitions are recorded
        as NOT fine-tested."""
        states = [
            _night("2025-03-01", _rot("stlmi", 0.45, dec=25.0)),
            _night("2025-03-02", _rot("vvpup", 0.70, ra=123.0, dec=-19.0)),
            _night("2025-03-03", _rot("ngc5906", 0.01, ra=229.0, dec=56.0)),
        ]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert len(epochs) == 1
        assert epochs[0]["n_tested_fine"] == 0
        assert epochs[0]["n_tested_coarse"] == 2

    def test_a_same_target_step_cuts_and_reports_its_size(self):
        states = [
            _night("2024-10-28", _rot("v826aur", 181.56, dec=45.0)),
            _night("2024-10-31", _rot("v826aur", 181.50, dec=45.0)),
            _night("2024-11-04", _rot("v826aur", 180.77, dec=45.0)),
        ]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert _causes(epochs) == [inv.CAUSE_CAMERA, inv.CAUSE_ROT_FINE]
        step = epochs[1]["step_deg"]
        assert step == pytest.approx(180.77 - 181.53, abs=1e-6)
        assert epochs[1]["step_basis"] == "same_target"
        assert epochs[1]["step_n_targets"] == 1

    def test_the_threshold_is_where_the_rule_says_it_is(self):
        # A hair under the step threshold: one epoch.  A hair over: two.
        # (Not tested AT the threshold: 10.0 + 0.3 - 10.0 is not 0.3 in
        # binary floating point, and a test must not depend on which side
        # the rounding falls.)
        for delta, n_epochs in ((inv.MECH_ROT_STEP_DEG - 0.01, 1),
                                (inv.MECH_ROT_STEP_DEG + 0.01, 2)):
            states = [_night("2025-01-01", _rot("a", 10.0)),
                      _night("2025-01-02", _rot("a", 10.0 + delta))]
            assert len(inv.segment_mech_epochs(states)[0]) == n_epochs

    def test_same_name_different_field_is_not_a_shared_target(self):
        # 'stars' is re-used for different fields; 20 deg apart on the sky
        # the pointing term does not cancel and the pair proves nothing.
        states = [_night("2025-01-01", _rot("stars", 10.0, ra=100.0)),
                  _night("2025-01-02", _rot("stars", 10.9, ra=120.0))]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert len(epochs) == 1
        assert epochs[0]["n_tested_fine"] == 0

    def test_the_monsoon_flip_is_found_without_a_shared_target(self):
        states = [
            _night("2025-06-23", _rot("rhooph", 0.45, dec=-23.0)),
            _night("2025-10-15", _rot("sao054471", 179.77, dec=36.0),
                   flip=""),
        ]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert len(epochs) == 2
        assert inv.CAUSE_FLIPSTAT in epochs[1]["boundary_cause"]
        assert inv.CAUSE_ROT_COARSE in epochs[1]["boundary_cause"]
        assert epochs[1]["step_deg"] == pytest.approx(179.32, abs=1e-6)
        assert epochs[1]["step_basis"] == "nightly_median"
        # A flip-state change is a DETECTOR boundary: darks do not cross it.
        assert epochs[0]["detector_epoch"] != epochs[1]["detector_epoch"]

    def test_a_reseat_keeps_the_detector_epoch(self):
        states = [_night("2024-11-30", _rot("m1", 204.2)),
                  _night("2024-12-20", _rot("m1", 190.0))]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert len(epochs) == 2
        assert epochs[0]["detector_epoch"] == epochs[1]["detector_epoch"]

    def test_high_declination_solves_are_ignored(self):
        # NGC 188 (dec +85) sits ~1.9 deg from its night's other fields with
        # nothing touched; beyond the limit a solve is not rotation evidence.
        states = [_night("2025-10-30", _rot("m31", 179.72, dec=41.0)),
                  _night("2025-10-31", _rot("ngc188", 177.60, dec=85.0))]
        epochs, nights, _ = inv.segment_mech_epochs(states)
        assert len(epochs) == 1
        assert nights[1]["n_rot"] == 0

    def test_too_few_solves_are_not_evidence(self):
        states = [_night("2025-01-01", _rot("a", 10.0)),
                  _night("2025-01-02", _rot("a", 40.0, n=2))]
        assert len(inv.segment_mech_epochs(states)[0]) == 1

    def test_a_swapped_wheel_slot_cuts(self):
        states = [_night("2024-12-02", wheel="lum|red|green|blue|empty"),
                  _night("2024-12-04", wheel="lum|red|green|blue|g")]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert _causes(epochs) == [inv.CAUSE_CAMERA, inv.CAUSE_WHEEL]
        # A wheel reload does not move a hot pixel.
        assert epochs[0]["detector_epoch"] == epochs[1]["detector_epoch"]

    def test_a_respelled_wheel_does_not_cut(self):
        states = [_night("2025-01-20", wheel="g|lrg|r|i|hrg|ha|dark"),
                  _night("2025-01-21", wheel="g|OGGrism|r|i|HaGrism|Ha|Dark"),
                  _night("2025-01-22", wheel="g|lrg|r|i|hrg|ha|dark")]
        assert len(inv.segment_mech_epochs(states)[0]) == 1

    def test_maps_of_different_lengths_are_never_compared(self):
        # November 2024: the same nights report a 16-slot 'Dual Wheels' map
        # and a 5-slot 'FLI' map.  They are two descriptions, not two states.
        dual = "|".join(f"s{i}" for i in range(16))
        states = [
            _night("2024-11-08", wheel=dual, maps=(dual,)),
            _night("2024-11-09", wheel="lum|red|green|blue|empty",
                   maps=("lum|red|green|blue|empty", dual)),
            # ... and a later change in the 5-slot family alone DOES cut.
            _night("2024-11-16", wheel="lum|red|green|blue|g",
                   maps=("lum|red|green|blue|g", dual)),
        ]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert [e["first_night"] for e in epochs] == [
            "2024-11-08",       # the camera's first night
            "2024-11-09",       # a slot count never reported before appears
            "2024-11-16"]       # a slot of the 5-slot wheel changes

    def test_a_new_slot_count_cuts_once_and_its_return_does_not(self):
        # 2026-06-28: the wheel went from 7 named slots to 9.  No slot-by-
        # slot comparison can see that; the new COUNT is the evidence.
        seven = "g|lrg|r|i|ha|hrg|empty"
        nine = seven + "|lum|sii"
        states = [_night("2026-06-27", wheel=seven),
                  _night("2026-06-28", wheel=nine),
                  _night("2026-06-29", wheel=seven),   # seen before: no cut
                  _night("2026-06-30", wheel=nine)]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert _causes(epochs) == [inv.CAUSE_CAMERA, inv.CAUSE_WHEEL]
        assert epochs[1]["first_night"] == "2026-06-28"

    def test_a_missing_card_never_cuts_and_never_forgets(self):
        states = [_night("2025-01-01", flip="Flip/Mirror"),
                  _night("2025-01-02", flip=None, wheel=None),
                  _night("2025-01-03", flip="Flip/Mirror"),
                  _night("2025-01-04", flip="")]
        epochs, _, _ = inv.segment_mech_epochs(states)
        assert [e["first_night"] for e in epochs] == ["2025-01-01",
                                                      "2025-01-04"]
        assert epochs[0]["flipstat"] == "Flip/Mirror"
        assert epochs[1]["flipstat"] == ""

    def test_another_camera_in_between_is_a_remount(self):
        states = [_night("2024-01-01", cam="AC4040"),
                  _night("2024-02-01", cam="iKon"),
                  _night("2024-03-01", cam="AC4040")]
        epochs, _, _ = inv.segment_mech_epochs(states)
        got = {e["mech_epoch"]: e["boundary_cause"] for e in epochs}
        assert got == {"AC4040:2024-01-01": inv.CAUSE_CAMERA,
                       "iKon:2024-02-01": inv.CAUSE_CAMERA,
                       "AC4040:2024-03-01": inv.CAUSE_CAMERA_SWAP}
        # A remounted camera is a new detector epoch: nothing proves its
        # flip state or readout survived the bench.
        det = {e["mech_epoch"]: e["detector_epoch"] for e in epochs}
        assert det["AC4040:2024-03-01"] == "AC4040:2024-03-01"
        assert [e["seq"] for e in epochs] == [1, 2, 3]

    def test_nights_between_a_rotation_only_boundary_are_uncertain(self):
        """The hardware moved somewhere between the last night that showed
        the old angle and the first that showed the new one.  The nights
        in between carry no rotation: they are filed under the earlier
        epoch and flagged, so no flat can be matched to them."""
        states = [
            _night("2025-01-01", _rot("a", 10.0)),
            _night("2025-01-02"),                 # no solves
            _night("2025-01-03"),                 # no solves
            _night("2025-01-04", _rot("a", 12.0)),
        ]
        epochs, nights, _ = inv.segment_mech_epochs(states)
        assert len(epochs) == 2
        by_night = {n["night"]: n for n in nights}
        assert by_night["2025-01-01"]["certain"] == 1
        assert by_night["2025-01-02"]["certain"] == 0
        assert by_night["2025-01-03"]["certain"] == 0
        assert by_night["2025-01-02"]["mech_epoch"] == "ASI:2025-01-01"
        assert by_night["2025-01-04"]["certain"] == 1
        assert epochs[0]["n_gap_nights"] == 2

    def test_a_header_dated_boundary_leaves_no_uncertain_nights(self):
        # The flip card dates the change exactly; nothing is ambiguous.
        states = [_night("2025-01-01", _rot("a", 10.0)),
                  _night("2025-01-02"),
                  _night("2025-01-04", _rot("b", 190.0, ra=10.0), flip="")]
        _, nights, _ = inv.segment_mech_epochs(states)
        assert all(n["certain"] == 1 for n in nights)

    def test_input_order_does_not_matter(self):
        states = [_night(f"2025-02-{d:02d}", _rot("a", 10.0 + (d > 4)))
                  for d in range(1, 9)]
        fwd = inv.segment_mech_epochs(states)[0]
        rev = inv.segment_mech_epochs(list(reversed(states)))[0]
        assert [e["mech_epoch"] for e in fwd] == [e["mech_epoch"]
                                                  for e in rev]


class TestCalibValidFor:
    SCI = inv.EpochTag("ASI:2024-12-13", "ASI:2024-12-13", True)

    def test_flat_needs_the_same_mechanical_epoch(self):
        assert inv.calib_valid_for("flat", self.SCI, self.SCI)
        other = inv.EpochTag("ASI:2025-10-15", "ASI:2024-12-13", True)
        assert not inv.calib_valid_for("flat", self.SCI, other)

    def test_flat_needs_certainty_on_both_sides(self):
        gap = inv.EpochTag("ASI:2024-12-13", "ASI:2024-12-13", False)
        assert not inv.calib_valid_for("flat", gap, self.SCI)
        assert not inv.calib_valid_for("flat", self.SCI, gap)

    def test_dark_and_bias_follow_the_detector_epoch(self):
        reseat = inv.EpochTag("ASI:2025-03-01", "ASI:2024-12-13", True)
        flipped = inv.EpochTag("ASI:2025-10-15", "ASI:2025-10-15", True)
        for kind in ("dark", "bias"):
            assert inv.calib_valid_for(kind, self.SCI, reseat)
            assert not inv.calib_valid_for(kind, self.SCI, flipped)

    def test_no_epoch_is_never_the_same_epoch(self):
        nowhere = inv.EpochTag(None, None, True)
        for kind in ("flat", "dark", "bias"):
            assert not inv.calib_valid_for(kind, nowhere, nowhere)

    def test_unknown_kind_raises(self):
        with pytest.raises(ValueError):
            inv.calib_valid_for("fringe", self.SCI, self.SCI)


class TestSettingsMatch:
    """DE.F5: era 76 pools Mode0 at -10 C with Mode0 at 0 C."""

    def test_same_settings_match(self):
        assert inv.settings_match("100", 30.0, -10.0, "100", 30.0, -10.0) \
            == inv.SETTINGS_MATCH

    def test_set_point_within_tolerance_matches(self):
        assert inv.settings_match("100", 30.0, -10.0, "100", 30.0, -8.5) \
            == inv.SETTINGS_MATCH

    def test_a_ten_degree_set_point_difference_does_not(self):
        # 80 of the 247 T CrB 240 s grism frames were taken at 0 to -2 C;
        # the four 240 s master darks are all at -10 C.
        assert inv.settings_match("100", 30.0, 0.0, "100", 30.0, -10.0) \
            == inv.SETTINGS_MISMATCH

    def test_gain_and_offset_must_agree(self):
        assert inv.settings_match("100", 30.0, -10.0, "56", 30.0, -10.0) \
            == inv.SETTINGS_MISMATCH
        assert inv.settings_match("100", 30.0, -10.0, "100", 76.0, -10.0) \
            == inv.SETTINGS_MISMATCH

    def test_a_calibration_frame_without_the_card_is_unverified(self):
        assert inv.settings_match("100", 30.0, -10.0, None, None, None) \
            == inv.SETTINGS_UNVERIFIED
        # ... but a KNOWN disagreement still wins over an unknown.
        assert inv.settings_match("100", 30.0, -10.0, None, None, 0.0) \
            == inv.SETTINGS_MISMATCH

    def test_what_the_science_frame_does_not_record_constrains_nothing(self):
        # AC4040 headers carry no GAIN/OFFSET card at all.
        assert inv.settings_match(None, None, -20.0, "x", 5.0, -20.0) \
            == inv.SETTINGS_MATCH
        assert inv.settings_match("", None, None, None, None, None) \
            == inv.SETTINGS_MATCH


# ---------------------------------------------------------------------------
# Build-script wiring for the mechanical layer (tiny synthetic manifest)
# ---------------------------------------------------------------------------
def _mech_frames() -> pd.DataFrame:
    """Six science frames either side of a flip, with calibration on both
    sides, all in ONE header era (the era-76 situation in miniature)."""
    def row(i, night, **over):
        base = dict(
            obs_rowid=i, path=f"rawimage/{night}/f{i}.fts",
            tree="rawimage", basename=f"f{i}.fts", jd=2460000.0 + i,
            night=night, target_key="tcrb", canonical_target="T CrB",
            imagetyp="Light Frame", filter="g", exptime=240.0, era_id=76,
            is_canonical=1, dup_group=i, error=None, readoutm="Mode0",
            camtemp=None, ccd_temp=None, camera="ASI",
            flipstat="Flip/Mirror", fwallnam="g|r|i", swcreate=None,
            telpier=None, hdr_gain="100", hdr_offset=30.0, set_temp=-10.0,
            dup_basis="canonical")
        base.update(over)
        return base
    rows = [
        row(1, "2025-05-01"), row(2, "2025-05-01"), row(3, "2025-05-02"),
        # after the flip; frame 6 also at a different set-point
        row(4, "2025-10-15", flipstat=""), row(5, "2025-10-15", flipstat=""),
        row(6, "2025-10-16", flipstat="", set_temp=0.0),
        # calibration BEFORE the flip: a flat, a matching dark, a bias
        row(10, "2025-05-01", imagetyp="Flat Field", exptime=2.0),
        row(11, "2025-05-01", imagetyp="Dark Frame", filter=None),
        row(12, "2025-05-01", imagetyp="Bias Frame", filter=None,
            exptime=0.0),
    ]
    df = pd.DataFrame(rows)
    df["calib_kind"] = [inv.calib_kind(it, bn) for it, bn
                        in zip(df["imagetyp"], df["basename"])]
    return df


class TestMechWiring:
    def _with_epochs(self):
        df = _mech_frames()
        mech = build.build_mech_tables(
            df, pd.DataFrame(columns=["obs_rowid", "target_key",
                                      "rotation_deg", "solved_ra",
                                      "solved_dec"]))
        return build.attach_epochs(df, mech["frame_mech_epoch"]), mech

    def test_night_states_take_the_modal_header_value(self):
        df = _mech_frames()
        df.loc[df["obs_rowid"] == 2, "flipstat"] = None     # card absent
        states = build.build_night_states(df, pd.DataFrame())
        by = {s.night: s for s in states}
        assert by["2025-05-01"].flipstat == "Flip/Mirror"
        assert by["2025-10-15"].flipstat == ""
        assert by["2025-05-01"].wheel_maps == ("g|r|i",)

    def test_the_flip_splits_one_era_into_two_epochs(self):
        df, mech = self._with_epochs()
        ep = mech["mech_epoch"]
        assert list(ep["mech_epoch"]) == ["ASI:2025-05-01", "ASI:2025-10-15"]
        assert inv.CAUSE_FLIPSTAT in ep["boundary_cause"].iloc[1]
        assert set(df[df["obs_rowid"] <= 3]["mech_epoch"]) == \
            {"ASI:2025-05-01"}
        assert set(df[df["obs_rowid"].isin([4, 5, 6])]["mech_epoch"]) == \
            {"ASI:2025-10-15"}

    def test_coverage_never_counts_a_calibration_across_the_flip(self):
        """THE F-3 RULE AT THE COUNTING LAYER.  The era-level matrix says
        era 76 has a g flat, a 240 s dark and a bias.  Per mechanical
        epoch: the pre-flip science has them; the post-flip science, in the
        same era, has NONE — the frames exist and may not be applied."""
        df, _ = self._with_epochs()
        cov, gaps = build.build_coverage_mech(df, {}, frozenset())
        pre = cov[cov["mech_epoch"] == "ASI:2025-05-01"]
        post = cov[cov["mech_epoch"] == "ASI:2025-10-15"]
        assert pre["n_calib_raw"].sum() == 3
        assert post["n_calib_raw"].sum() == 0
        assert post["n_calib_master"].sum() == 0
        assert set(post["status"]) == {"missing"}
        # The era-level matrix, for contrast, counts all three for the
        # whole era — which is exactly the defect.
        era_cov, _ = build.build_coverage(df, {}, frozenset())
        assert era_cov["n_calib_raw"].sum() == 3
        assert set(gaps["mech_epoch"]) == {"ASI:2025-05-01",
                                           "ASI:2025-10-15"}

    def test_set_point_classes_are_separate_requirements(self):
        df, _ = self._with_epochs()
        cov, _ = build.build_coverage_mech(df, {}, frozenset())
        post_bias = cov[(cov["mech_epoch"] == "ASI:2025-10-15")
                        & (cov["req_kind"] == "bias")]
        assert sorted(post_bias["req_key"]) == [
            "+0C, gain 100, offset 30", "-10C, gain 100, offset 30"]
        assert list(post_bias.sort_values("req_key")["n_science"]) == [1, 2]

    def test_a_dark_at_the_wrong_set_point_is_not_counted(self):
        df, _ = self._with_epochs()
        # Move the pre-flip science to 0 C: the -10 C dark no longer serves.
        df.loc[df["obs_rowid"].isin([1, 2, 3]), "set_temp"] = 0.0
        cov, _ = build.build_coverage_mech(df, {}, frozenset())
        dark = cov[(cov["mech_epoch"] == "ASI:2025-05-01")
                   & (cov["req_kind"] == "dark")]
        assert dark["n_calib_raw"].sum() == 0

    def test_a_calibration_with_unrecorded_settings_is_unverified(self):
        df, _ = self._with_epochs()
        df.loc[df["obs_rowid"] == 11, ["hdr_gain", "hdr_offset",
                                       "set_temp"]] = [None, None, None]
        cov, _ = build.build_coverage_mech(df, {}, frozenset())
        dark = cov[(cov["mech_epoch"] == "ASI:2025-05-01")
                   & (cov["req_kind"] == "dark")].iloc[0]
        assert dark["n_calib_raw"] == 0 and dark["n_unverified"] == 1

    def test_calib_census_carries_the_epoch_columns_last(self):
        df, _ = self._with_epochs()
        calib = build.build_calib_frames(df)
        assert list(calib.columns)[-7:] == [
            "camera", "mech_epoch", "detector_epoch", "epoch_certain",
            "hdr_gain", "hdr_offset", "set_temp"]
        assert set(calib["mech_epoch"]) == {"ASI:2025-05-01"}

    def test_frames_with_no_camera_get_no_epoch(self):
        df = _mech_frames()
        df.loc[df["obs_rowid"] == 3, "camera"] = None
        mech = build.build_mech_tables(df, pd.DataFrame())
        assert 3 not in set(mech["frame_mech_epoch"]["obs_rowid"])
        out = build.attach_epochs(df, mech["frame_mech_epoch"])
        row = out[out["obs_rowid"] == 3].iloc[0]
        assert pd.isna(row["mech_epoch"]) and row["epoch_certain"] == 0

    def test_null_summary_states_the_false_alarm_numbers(self):
        null = pd.DataFrame({"delta_deg": [0.01, -0.02, 0.0, 0.35, -0.01]})
        got = build.null_summary(null)
        assert got["mech_null_n"] == "5"
        assert got["mech_null_frac_above_step"] == "0.20000"
        assert got["mech_null_max_abs_deg"] == "0.3500"
        assert build.null_summary(null.iloc[0:0]) == {"mech_null_n": "0"}

    def test_injected_steps_are_recovered_where_targets_repeat(
            self, monkeypatch):
        """The power test, on a fixture where the answer is known: one
        target every night, so every transition is fine-testable.  A 1 deg
        step must always be found on the night it happens; a 0.2 deg step
        (below the rule) must never be."""
        monkeypatch.setattr(build, "INJECTION_TRIALS", 20)
        monkeypatch.setattr(build, "INJECTION_STEPS_DEG", (0.2, 1.0))
        states = [_night(f"2025-02-{d:02d}", _rot("stlmi", 0.45))
                  for d in range(1, 15)]
        epochs, _, _ = inv.segment_mech_epochs(states)
        inj = build.injection_test(states, pd.DataFrame(epochs))
        got = {r["step_deg"]: r for _, r in inj.iterrows()}
        assert got[1.0]["frac_detected"] == 1.0
        assert got[1.0]["frac_exact"] == 1.0
        assert got[1.0]["median_delay_nights"] == 0.0
        assert got[0.2]["frac_detected"] == 0.0

    def test_injection_into_an_epoch_with_no_repeats_is_blind(
            self, monkeypatch):
        """...and the honest counterpart: a different field every night
        gives the fine test nothing, so a 1 deg step (below the coarse
        threshold) is invisible — which the table must SAY, not hide."""
        monkeypatch.setattr(build, "INJECTION_TRIALS", 20)
        monkeypatch.setattr(build, "INJECTION_STEPS_DEG", (1.0, 5.0))
        states = [_night(f"2025-02-{d:02d}",
                         _rot(f"field{d}", 0.45, ra=10.0 * d))
                  for d in range(1, 15)]
        epochs, _, _ = inv.segment_mech_epochs(states)
        inj = build.injection_test(states, pd.DataFrame(epochs))
        got = {r["step_deg"]: r for _, r in inj.iterrows()}
        assert got[1.0]["frac_detected"] == 0.0
        assert got[5.0]["frac_detected"] == 1.0      # coarse test sees it


class TestLinkMethodNaming:
    def test_s0_deduplicated_twins_keep_the_published_vocabulary(self):
        """S2's reconstruction experiment selects links by method name.  A
        ``_calibrated`` twin is ``stem_jd`` whether S0's dedup grouped it
        (v1.1) or this stage's ladder found it (v1.0)."""
        rows = [
            dict(obs_rowid=1, path="rawimage/n/a.fts.fz", tree="rawimage",
                 basename="a.fts.fz", jd=2460000.5, night="n",
                 target_key="t", dup_group=1, is_canonical=1,
                 dup_basis="canonical"),
            dict(obs_rowid=2, path="reduced/n/a.fts.fz", tree="reduced",
                 basename="a.fts.fz", jd=2460000.5, night="n",
                 target_key="t", dup_group=1, is_canonical=0,
                 dup_basis="same_basename_jd"),
            dict(obs_rowid=3, path="reduced/n/a_calibrated.fts.fz",
                 tree="reduced", basename="a_calibrated.fts.fz",
                 jd=2460000.5, night="n", target_key="t", dup_group=1,
                 is_canonical=0, dup_basis="processing_suffix"),
            dict(obs_rowid=4, path="reduced/n/a_1_wcs.fts.fz",
                 tree="reduced", basename="a_1_wcs.fts.fz", jd=2460000.5,
                 night="n", target_key="t", dup_group=1, is_canonical=0,
                 dup_basis="suffix_counter"),
        ]
        links = build.build_links(pd.DataFrame(rows))
        got = dict(zip(links["reduced_rowid"], links["match_method"]))
        assert got == {2: "same_basename_jd", 3: "stem_jd",
                       4: "stem_jd_counter"}
        assert set(links["raw_rowid"]) == {1}
