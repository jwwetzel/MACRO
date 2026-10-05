"""Unit tests for macro_legacy — the legacy-archive census rules.

Each test class mirrors one decision the census encodes.  The must-NOT
cases (a focus frame counted as science, two cameras pooled into one
series, a pier flip read as a re-mount, a 1970 clock read as an observation,
an unevaluated gate treated as a pass) matter as much as the must cases:
every one of them would silently move the pre-registered go/no-go.

The gate thresholds are asserted against the numbers written in
``Legacy_Rigel/notes/GO_NOGO_PREREGISTERED.md`` — if someone edits a
constant, the test that quotes the note fails.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest \
        pipeline/tests/test_legacy_census.py -q
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

# Make the packages importable regardless of pytest's working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
# The build scripts' DataFrame paths are also under test (the S0 era_id
# lesson: a bug shipped because only the pure layer was tested).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from macro_legacy import census as lc
from macro_legacy import scan as lscan

import build_legacy_census as build
import build_legacy_external as external
import build_legacy_scan as scanner


# ---------------------------------------------------------------------------
# Scan: header -> row
# ---------------------------------------------------------------------------
class TestScanRow:
    def test_column_names(self):
        assert lscan.col("CCD-TEMP") == "ccd_temp"
        assert lscan.col("DATE-OBS") == "date_obs"

    def test_strings_stripped_numbers_floated(self):
        row = lscan.row_from_header({
            "INSTRUME": "  Andor IKON L 936 ", "EXPTIME": "90.0",
            "NAXIS1": 2048, "NAXIS2": 2048, "BITPIX": 16})
        assert row["instrume"] == "Andor IKON L 936"
        assert row["exptime"] == 90.0
        assert (row["naxis1"], row["naxis2"], row["bitpix"]) == (2048, 2048, 16.0)

    def test_missing_card_is_null_not_zero(self):
        row = lscan.row_from_header({"NAXIS1": 10, "NAXIS2": 10})
        assert row["egain"] is None and row["imagetyp"] is None

    def test_blank_string_is_kept_as_empty(self):
        # "the card exists and is empty" is information (blank INSTRUME).
        row = lscan.row_from_header({"NAXIS1": 1, "NAXIS2": 1, "INSTRUME": " "})
        assert row["instrume"] == ""

    def test_logical_is_not_a_number(self):
        row = lscan.row_from_header({"NAXIS1": 1, "NAXIS2": 1, "GAIN": True})
        assert row["gain"] is None

    def test_compressed_geometry_uses_znaxis(self):
        # The S0e trap: a BINTABLE's NAXIS1/2 are row bytes and row count.
        row = lscan.row_from_header({
            "ZIMAGE": True, "ZCMPTYPE": "RICE_1", "NAXIS1": 8, "NAXIS2": 3211,
            "ZNAXIS1": 4800, "ZNAXIS2": 3211, "BITPIX": 8, "ZBITPIX": 16})
        assert (row["naxis1"], row["naxis2"]) == (4800, 3211)
        assert row["bitpix"] == 16.0

    def test_unresolvable_geometry_is_an_error_not_a_guess(self):
        row = lscan.row_from_header({"ZIMAGE": True, "NAXIS1": 8, "NAXIS2": 9})
        assert row["naxis1"] is None and "GeometryError" in row["error"]

    def test_keyset_signature_ignores_container_cards(self):
        a, text = lscan.keyset_signature(["JD", "FILTER", "COMMENT", "ZIMAGE",
                                          "CHECKSUM", "NAXIS1"])
        b, _ = lscan.keyset_signature(["FILTER", "JD", "HISTORY", "TFORM1"])
        assert a == b and text == "FILTER,JD"

    def test_keyset_signature_distinguishes_writers(self):
        talon, _ = lscan.keyset_signature(["JD", "RAWHENC", "CAMTEMP"])
        maxim, _ = lscan.keyset_signature(["JD", "SWCREATE", "CCD-TEMP"])
        assert talon != maxim

    def test_scan_one_never_raises_on_a_missing_file(self, tmp_path):
        row = lscan.scan_one((str(tmp_path), "nope.fts.fz"))
        assert row["path"] == "nope.fts.fz"
        assert row["error"].startswith("FileNotFoundError")

    def test_scan_one_reads_a_real_compressed_file(self, tmp_path):
        fits = pytest.importorskip("astropy.io.fits")
        import numpy as np
        hdr = fits.Header()
        hdr["INSTRUME"] = "SBIG Aluma AC4040"
        hdr["DATE-OBS"] = "2022-03-01T05:00:00"
        hdr["EXPTIME"] = 30.0
        hdr["FILTER"] = "G - Sloan_g"
        hdr["COMMENT"] = "programme note"
        hdu = fits.CompImageHDU(np.zeros((20, 30), dtype="int16"), header=hdr)
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(tmp_path / "a.fts.fz")
        row = lscan.scan_one((str(tmp_path), "a.fts.fz"))
        assert row["error"] is None
        assert (row["naxis1"], row["naxis2"]) == (30, 20)
        assert row["instrume"] == "SBIG Aluma AC4040"
        assert row["exptime"] == 30.0 and row["n_hdu"] == 2
        assert "programme note" in row["comment"]


class TestScanScript:
    def test_file_kind(self):
        assert scanner.file_kind("frm00400.fts.fz") == "fz"
        assert scanner.file_kind("frm00400.FTS") == "fits"
        assert scanner.file_kind("notes.txt") == "other"

    def test_walk_needs_no_stat_and_classifies(self, tmp_path):
        (tmp_path / "2015" / "day004").mkdir(parents=True)
        (tmp_path / "2015" / "day004" / "a.fts.fz").write_bytes(b"x")
        (tmp_path / "2015" / "day004" / "a.fts").write_bytes(b"x")
        (tmp_path / "readme.txt").write_text("x")
        got = dict(scanner.walk_archive(tmp_path))
        assert got == {"2015/day004/a.fts.fz": "fz",
                       "2015/day004/a.fts": "fits", "readme.txt": "other"}

    def test_todo_skips_clean_rows_and_retries_errors(self, tmp_path):
        con = sqlite3.connect(tmp_path / "t.sqlite")
        scanner.ensure_scan_tables(con)
        scanner.write_disk_files(con, [("a.fz", "fz"), ("b.fz", "fz"),
                                       ("c.fz", "fz"), ("d.fts", "fits")])
        con.execute("INSERT INTO scan (path) VALUES ('a.fz')")
        con.execute("INSERT INTO scan (path, error) VALUES ('b.fz', 'boom')")
        assert scanner.todo_paths(con) == ["b.fz", "c.fz"]

    def test_seed_paths_append_fz_and_skip_done(self, tmp_path):
        con = sqlite3.connect(tmp_path / "t.sqlite")
        scanner.ensure_scan_tables(con)
        con.execute("INSERT INTO scan (path) VALUES ('2015/day004/a.fts.fz')")
        csv = tmp_path / "m.csv"
        csv.write_text("file_id,rel_path,name\n1,2015/day004/a.fts,a.fts\n"
                       "2,2015/day004/b.fts,b.fts\n")
        assert scanner.seed_paths(con, csv) == ["2015/day004/b.fts.fz"]


# ---------------------------------------------------------------------------
# 1. Paths and reconciliation
# ---------------------------------------------------------------------------
class TestReconciliation:
    def test_logical_path(self):
        assert lc.logical_path("2015/day004/a.fts.fz") == "2015/day004/a.fts"
        assert lc.logical_path("2015/day004/a.fts") == "2015/day004/a.fts"

    def test_readable_fz_is_a_row_even_with_a_twin(self):
        assert lc.reconcile_one(True, True, True, False) == lc.REC_SCANNED

    def test_unreadable_fz_is_named(self):
        assert lc.reconcile_one(True, True, False, True) == lc.REC_UNREADABLE

    def test_uncompressed_only_is_named(self):
        assert lc.reconcile_one(True, False, True, False) \
            == lc.REC_UNCOMPRESSED_ONLY

    def test_manifest_path_absent_from_disk_is_named(self):
        assert lc.reconcile_one(True, False, False, False) \
            == lc.REC_NOT_ON_DISK

    def test_non_fits_file_is_named(self):
        assert lc.reconcile_one(False, False, False, False,
                                is_fits_name=False) == lc.REC_NOT_FITS

    def test_unreachable_path_is_refused(self):
        with pytest.raises(ValueError):
            lc.reconcile_one(False, False, False, False, is_fits_name=True)

    def test_strip_first_component_reproduces_the_bad_manifest(self):
        assert lc.strip_first_component("2015/day004/frm00400.fts") \
            == "day004/frm00400.fts"
        assert lc.strip_first_component("2017/day004/frm00400.fts") \
            == "day004/frm00400.fts"          # … the collision
        assert lc.strip_first_component("solo.fts") == "solo.fts"


# ---------------------------------------------------------------------------
# 2. Time
# ---------------------------------------------------------------------------
class TestTime:
    def test_parse_talon_and_maxim_stamps(self):
        # Talon sample: DATE-OBS 2015-01-04T12:10:40.558 <-> JD 2457027.007413871
        jd = lc.parse_date_obs("2015-01-04T12:10:40.558")
        assert abs(jd - 2457027.007413871) * 86400 < 0.01
        # MaxIm sample: 2017-04-10T04:40:25 <-> JD 2457853.694733796
        jd = lc.parse_date_obs("2017-04-10T04:40:25")
        assert abs(jd - 2457853.694733796) * 86400 < 0.01

    def test_unset_clock_is_not_an_observation(self):
        assert lc.parse_date_obs("1970-01-01T00:00:02.670") is None

    @pytest.mark.parametrize("bad", [None, "", "2017-13-40T00:00:00",
                                     "10/04/2017", "2017-04-10"])
    def test_malformed_is_none(self, bad):
        assert lc.parse_date_obs(bad) is None

    def test_decimals(self):
        assert lc.date_obs_decimals("2017-04-10T04:40:25") == 0
        assert lc.date_obs_decimals("2015-01-04T12:10:40.558") == 3
        assert lc.date_obs_decimals(None) is None

    def test_night_is_local_evening_date(self):
        # 04:40 UT on 10 April is the night of 9 April at Winer (UTC-7).
        assert lc.night_label(lc.parse_date_obs("2017-04-10T04:40:25")) \
            == "2017-04-09"
        # 02:00 UT (19:00 local, the rollover) belongs to the new night.
        assert lc.night_label(lc.parse_date_obs("2017-04-10T19:30:00")) \
            == "2017-04-10"

    def test_path_date_layouts(self):
        assert lc.path_date("2015/day004/frm00400.fts.fz") == (2015, 4)
        assert lc.path_date("archivar/2017/foc0050f.fts.fz") == (2017, 5)
        assert lc.path_date("old/day336/gac33615.fts.fz") == (None, 336)
        assert lc.path_date("stray.fts.fz") == (None, None)

    def test_path_offset_zero_when_directory_is_ut_date(self):
        jd = lc.parse_date_obs("2015-01-04T12:10:40.558")
        assert lc.path_date_offset_days("2015/day004/frm0040f.fts.fz", jd) == 0
        assert lc.path_date_offset_days("2015/day005/x.fts.fz", jd) == -1

    def test_path_offset_year_end_and_unknowns(self):
        jd = lc.parse_date_obs("2016-01-01T03:00:00")
        assert lc.path_date_offset_days("2015/day365/x.fts.fz", jd) == 1
        assert lc.path_date_offset_days("old/day001/x.fts.fz", jd) == 0
        assert lc.path_date_offset_days("stray.fts.fz", jd) is None
        assert lc.path_date_offset_days("2015/day004/x.fts.fz", None) is None

    def test_sexagesimal(self):
        assert lc.parse_sexagesimal("12:45:08.09") == pytest.approx(12.752247)
        assert lc.parse_sexagesimal(" -0:27:47.7") == pytest.approx(-0.463250)
        assert lc.parse_sexagesimal("08 22 35") == pytest.approx(8.376389)
        assert lc.parse_sexagesimal("garbage") is None

    def test_lst_residual_on_real_headers(self):
        # Talon frame: JD 2457027.007413871, LST 11:43:25 (1-s card).
        r = lc.lst_residual_seconds(2457027.007413871, "11:43:25")
        assert abs(r) < 2.0
        # MaxIm frame: JD 2457853.694733796, LST 10:32:26.19.
        r = lc.lst_residual_seconds(2457853.694733796, "10:32:26.19")
        assert abs(r) < 2.0

    def test_lst_residual_exposes_local_time(self):
        # The same frame stamped in local time (UTC-7) is 7 sidereal-ish
        # hours off — the failure this check exists to catch.
        r = lc.lst_residual_seconds(2457853.694733796 - 7 / 24, "10:32:26.19")
        assert abs(abs(r) - 7 * 3600) < 120

    def test_lst_residual_none_when_unknown(self):
        assert lc.lst_residual_seconds(None, "10:00:00") is None
        assert lc.lst_residual_seconds(2457853.5, None) is None


# ---------------------------------------------------------------------------
# 3. Camera identity and mechanical timeline
# ---------------------------------------------------------------------------
class TestCamera:
    @pytest.mark.parametrize("text, name", [
        ("FLI ProLine PL16803 Rev 1.01", "FLI PL16803"),
        ("Apogee F47 1Kx1K 13u pixel", "Apogee F47"),
        ("SBIG 6303e 3072x2047 9u pixels", "SBIG 6303e"),
        ("SBIG6303e", "SBIG 6303e"),
        ("SBIG STXL-6303 3 CCD Camera", "SBIG STXL-6303"),
        ("Andor Aspen  CG42", "Andor Aspen CG42"),
        ("Andor IKON L 936", "Andor iKon-L 936"),
        ("SBIG Aluma AC4040", "SBIG Aluma AC4040"),
        ("SBIG STXL-6303e", "SBIG STXL-6303"),
        ("IKON L936", "Andor iKon-L 936"),
        ("DL Imaging", "SBIG Aluma AC4040"),
        ("Andor Tech", "CCD42-40 camera (model not in header)"),
    ])
    def test_every_observed_instrume_is_recognised(self, text, name):
        assert lc.camera_name(text) == name

    def test_stxl_is_not_folded_into_6303e(self):
        # Same sensor, different camera: rule order must keep them apart.
        assert lc.camera_name("SBIG STXL-6303 3 CCD Camera") != "SBIG 6303e"

    def test_blank_is_unlabelled_unknown_is_none(self):
        assert lc.camera_name("") == lc.CAMERA_UNLABELLED
        assert lc.camera_name(None) == lc.CAMERA_UNLABELLED
        assert lc.camera_name("QHY600") is None

    def test_label_accepted_when_sensor_agrees(self):
        assert lc.camera_identity("Apogee F47 1Kx1K 13u pixel", 13.0, 1024) \
            == ("Apogee F47", "label")

    def test_label_accepted_without_a_signature(self):
        # Talon headers carry no pixel size: nothing contradicts the label.
        assert lc.camera_identity("FLI ProLine PL16803 Rev 1.01", None, 4096) \
            == ("FLI PL16803", "label")

    def test_stale_label_loses_to_the_sensor(self):
        # INSTRUME says F47 (1024 px, 13 µm); the frame is 3072 px of 9 µm.
        cam, basis = lc.camera_identity("Apogee F47 1Kx1K 13u pixel", 9.0, 3072)
        assert basis == "sensor (label stale)"
        assert cam == "KAF-6303 camera (model not in header)"

    def test_stale_label_caught_by_width_alone(self):
        cam, basis = lc.camera_identity("Apogee F47 1Kx1K 13u pixel", None, 3072)
        assert basis == "unidentified" and cam == lc.CAMERA_UNLABELLED

    def test_6303_label_on_a_ccd42_frame_is_stale(self):
        # Observed 2016-10 and 2018-02: label not updated at a camera swap.
        assert lc.camera_identity("SBIG 6303e", 13.5, 2048) == (
            "CCD42-40 camera (model not in header)", "sensor (label stale)")

    def test_no_label_attributed_by_sensor(self):
        assert lc.camera_identity("", 9.0, 3072) == (
            "KAF-6303 camera (model not in header)", "sensor (no label)")

    def test_subframe_does_not_contradict(self):
        # A 512-px subframe of the F47 is still the F47.
        assert lc.camera_identity("Apogee F47 1Kx1K 13u pixel", 13.0, 512)[1] \
            == "label"

    def test_unrecognised_label_is_not_guessed(self):
        assert lc.camera_identity("QHY600", 3.76, 9576) \
            == (None, "unrecognised label")

    def test_binning_prefers_xbinning(self):
        assert lc.binning_of(2.0, 1.0) == 2
        assert lc.binning_of(None, 2.0) == 2
        assert lc.binning_of(None, None) is None
        assert lc.binning_of(float("nan"), 2.0) == 2

    def test_fold_rotation_removes_pier_flip(self):
        assert lc.fold_rotation(181.8) == pytest.approx(1.8)
        assert lc.fold_rotation(-178.2) == pytest.approx(1.8)
        assert lc.fold_rotation(1.8) == pytest.approx(1.8)

    def test_rotation_step_starts_an_epoch(self):
        nightly = [("2016-01-01", 1.0), ("2016-01-02", 1.1),
                   ("2016-01-03", 3.6), ("2016-01-04", 3.7)]
        ep = lc.rotation_epochs(nightly)
        assert [(e[0], e[1], e[3]) for e in ep] == [
            ("2016-01-01", "2016-01-02", 2), ("2016-01-03", "2016-01-04", 2)]

    def test_drift_within_tolerance_is_one_epoch(self):
        nightly = [(f"2016-01-{d:02d}", 1.0 + 0.02 * d) for d in range(1, 11)]
        assert len(lc.rotation_epochs(nightly)) == 1

    def test_one_night_excursion_is_counted_not_promoted(self):
        nightly = [("a", 1.0), ("b", 1.05), ("c", -2.5), ("d", 1.0), ("e", 0.95)]
        ep = lc.rotation_epochs(nightly)
        assert len(ep) == 1 and ep[0][3] == 4 and ep[0][4] == 1

    def test_one_night_between_different_epochs_is_kept(self):
        nightly = [("a", 1.0), ("b", 1.0), ("c", 2.0), ("d", 3.6), ("e", 3.6)]
        assert len(lc.rotation_epochs(nightly)) == 3

    def test_rotation_wraps_at_90(self):
        # +89.9 and -89.9 are 0.2 deg apart on the 180-deg circle.
        assert len(lc.rotation_epochs([("a", 89.9), ("b", -89.9)])) == 1


# ---------------------------------------------------------------------------
# 4. Frame kind, filters, targets
# ---------------------------------------------------------------------------
class TestKind:
    def test_focus_prefix_beats_light_frame(self):
        # The must-NOT: an autofocus frame carries 'Light Frame' and a star.
        assert lc.frame_kind("Light Frame", "foc10028.fts.fz", "Vega", 5.0) \
            == lc.KIND_FOCUS

    def test_maxim_vocabulary(self):
        assert lc.frame_kind("Light Frame", "gbm10028.fts.fz", "BH_Lyn", 90) \
            == lc.KIND_LIGHT
        assert lc.frame_kind("Bias Frame", "x.fts.fz", "", 0) == lc.KIND_BIAS
        assert lc.frame_kind("Dark Frame", "x.fts.fz", "", 60) == lc.KIND_DARK
        assert lc.frame_kind("Flat Field", "x.fts.fz", "", 3) == lc.KIND_FLAT

    def test_talon_without_imagetyp(self):
        assert lc.frame_kind(None, "frm0040f.fts.fz", "NGC 4666", 60) \
            == lc.KIND_LIGHT
        assert lc.frame_kind(None, "frm0040f.fts.fz", "flat V", 3) \
            == lc.KIND_FLAT
        assert lc.frame_kind(None, "frm0040f.fts.fz", "", 0) == lc.KIND_BIAS

    def test_engineering_objects_are_tests(self):
        assert lc.frame_kind("Light Frame", "gai1.fts.fz", "test", 5) \
            == lc.KIND_TEST
        assert lc.frame_kind("Light Frame", "gai1.fts.fz", "Pointing4W", 5) \
            == lc.KIND_TEST

    def test_filter_or_name_never_makes_a_calibration(self):
        # A science target is not a flat because of how it is spelled.
        assert lc.frame_kind("Light Frame", "x.fts.fz", "Flaming Star", 60) \
            == lc.KIND_LIGHT
        assert lc.frame_kind("Light Frame", "x.fts.fz", "Dark Doodad", 60) \
            == lc.KIND_LIGHT


class TestFilters:
    def test_talon_letters_are_johnson(self):
        f = lc.filter_info("V", talon=True)
        assert (f.band, f.system, f.tieable) == ("V", "Johnson-Cousins", True)

    def test_bare_letter_under_maxim_is_only_a_slot(self):
        # The must-NOT: MaxIm's 'R' is a wheel slot, not Johnson R.
        f = lc.filter_info("R", talon=False)
        assert f.band == "slot R" and f.tieable is False
        assert lc.filter_info("6").band == "slot 6"

    @pytest.mark.parametrize("raw, band", [
        ("G - Sloan_g", "g"), ("G - Sloan 'g'", "g"), ("G - Sloang", "g"),
        ("G - Sloan g", "g"), ("R - Sloan_r", "r"), ("R - Sloan 'r'", "r"),
        ("I - Sloan_i", "i"), ("I - Sloan 'i'", "i"),
        ("B - J-C blue", "B"), ("B - JC Blue", "B"), ("B -J-C blue", "B")])
    def test_named_standard_systems_are_tieable(self, raw, band):
        f = lc.filter_info(raw)
        assert (f.band, f.tieable) == (band, True)

    @pytest.mark.parametrize("raw, band", [
        ("R - Red", "Red"), ("X -Astrodon_ Red", "Red"), ("B - Blue", "Blue"),
        ("G - Green", "Green"), ("V - Visual", "Visual"),
        ("V - Astrodon  Visual", "Visual")])
    def test_named_colour_is_not_a_standard_band(self, raw, band):
        # 'R - Red' does not say which glass: never tieable, and never the
        # same series as Johnson R.
        f = lc.filter_info(raw)
        assert (f.band, f.tieable) == (band, False)

    @pytest.mark.parametrize("raw", [
        "G - GRISM", "6 - Grism600", "T - Transmission grism",
        "3- Grism300", "A - Grism300a", "S - spectrometer beamplitter"])
    def test_grisms_are_not_filters(self, raw):
        f = lc.filter_info(raw)
        assert f.system == "dispersive" and not f.tieable

    def test_green_is_not_sloan_g(self):
        assert lc.filter_info("G - Green").band != lc.filter_info(
            "G - Sloan_g").band

    def test_narrow_and_other_spellings_merge(self):
        assert {lc.filter_info(x).band for x in (
            "H - Halpha", "H - HAlpha", "H - Hydrogen Alpha")} == {"Ha"}
        assert {lc.filter_info(x).band for x in (
            "O - OIII", "O - Oxygen", "O - OxygenIII", "O - Oxygen III")} \
            == {"OIII"}
        assert {lc.filter_info(x).band for x in (
            "W-  Longpass", "W - LongpassIR", "W - Longpass IR")} == {"W"}
        assert {lc.filter_info(x).band for x in (
            "L - Luminance", "L- luminance", "L - Luminances")} == {"Lum"}

    def test_unknown_is_its_own_band(self):
        f = lc.filter_info("Q - Mystery")
        assert f.system == "unrecognised" and f.band == "Q - Mystery"
        assert f.tieable is False

    def test_missing(self):
        assert lc.filter_info(None).system == "missing"
        assert lc.filter_info("  ").band == "(none)"

    def test_every_tieable_band_is_in_the_preregistered_list(self):
        allowed = {"B", "V", "R", "I", "g", "r", "i"}
        rules = [i for _p, i in lc.FILTER_NAME_RULES] \
            + list(lc.TALON_LETTERS.values())
        assert {f.band for f in rules if f.tieable} <= allowed


class TestTargets:
    def test_alias_merge_matches_s0(self):
        assert lc.target_key("T CrB")[0] == lc.target_key("t_crb")[0] == "tcrb"
        assert lc.target_key("ST-LMi")[0] == "stlmi"
        assert lc.target_key("RW COM")[0] == lc.target_key("RW Com")[0]

    def test_blank_has_no_key(self):
        assert lc.target_key("")[0] is None and lc.target_key(None)[0] is None


# ---------------------------------------------------------------------------
# 5. Dedup
# ---------------------------------------------------------------------------
class TestDedup:
    def test_same_header_identity_is_one_exposure(self):
        a = lc.dup_key("cam", "2017-04-10T04:40:25", 90.0, "g", 2048, 2048, 1)
        b = lc.dup_key("cam", "2017-04-10T04:40:25", 90.0, "g", 2048, 2048, 2)
        assert a == b

    @pytest.mark.parametrize("change", [
        dict(camera="other"), dict(date_obs="2017-04-10T04:40:26"),
        dict(exptime=60.0), dict(band="r"), dict(naxis1=1024)])
    def test_any_differing_field_is_a_different_exposure(self, change):
        base = dict(camera="cam", date_obs="2017-04-10T04:40:25", exptime=90.0,
                    band="g", naxis1=2048, naxis2=2048)
        k0 = lc.dup_key(base["camera"], base["date_obs"], base["exptime"],
                        base["band"], base["naxis1"], base["naxis2"], 1)
        alt = {**base, **change}
        k1 = lc.dup_key(alt["camera"], alt["date_obs"], alt["exptime"],
                        alt["band"], alt["naxis1"], alt["naxis2"], 2)
        assert k0 != k1

    def test_unset_clock_frames_never_merge(self):
        a = lc.dup_key("cam", "1970-01-01T00:00:02.670", 5.0, "W", 4096, 4096, 1)
        b = lc.dup_key("cam", "1970-01-01T00:00:02.670", 5.0, "W", 4096, 4096, 2)
        assert a != b

    def test_raw_beats_calibrated_then_smallest_path(self):
        members = [("BDF", "2019/day100/a.fts.fz"), ("", "old/day100/a.fts.fz"),
                   (None, "2019/day100/b.fts.fz")]
        best = min(members, key=lambda m: lc.canonical_rank(*m))
        assert best == (None, "2019/day100/b.fts.fz")


# ---------------------------------------------------------------------------
# 6. Runs and seasons
# ---------------------------------------------------------------------------
class TestRunsSeasons:
    def test_preregistered_constants(self):
        assert lc.RUN_GAP_DAYS == pytest.approx(30 / 1440)
        assert lc.SEASON_GAP_DAYS == 120

    def test_gap_over_30_min_splits(self):
        t0 = 2457000.5
        jds = [t0, t0 + 10 / 1440, t0 + 20 / 1440,          # run 1: 20 min
               t0 + 60 / 1440, t0 + 65 / 1440]              # 40-min gap
        runs = lc.split_runs(jds)
        assert [r[2] for r in runs] == [3, 2]
        assert (runs[0][1] - runs[0][0]) * 1440 == pytest.approx(20)

    def test_gap_of_exactly_30_min_does_not_split(self):
        t0 = 2457000.5
        assert len(lc.split_runs([t0, t0 + lc.RUN_GAP_DAYS])) == 1

    def test_single_frame_is_a_zero_length_run(self):
        assert lc.split_runs([2457000.5]) == [(2457000.5, 2457000.5, 1)]

    def test_unsorted_input(self):
        t0 = 2457000.5
        assert lc.split_runs([t0 + 0.01, t0, t0 + 0.005])[0][2] == 3

    def test_seasons_split_at_gaps_over_120_days(self):
        s = lc.assign_seasons(["2016-03-01", "2016-04-15", "2016-08-20",
                               "2017-03-01"])
        assert s["2016-03-01"] == s["2016-04-15"] == 1
        assert s["2016-08-20"] == 2          # 127-day gap
        assert s["2017-03-01"] == 3

    def test_gap_of_exactly_120_days_is_one_season(self):
        s = lc.assign_seasons(["2016-01-01", "2016-04-30"])
        assert set(s.values()) == {1}


# ---------------------------------------------------------------------------
# 7. Calibration availability and flat pairs
# ---------------------------------------------------------------------------
class TestCalibration:
    def test_window_is_30_nights_inclusive(self):
        assert lc.CALIB_WINDOW_NIGHTS == 30
        assert lc.within_window(1000, [970])
        assert lc.within_window(1000, [1030])
        assert not lc.within_window(1000, [969, 1031])
        assert not lc.within_window(1000, [])

    def test_calibrated_needs_flat_and_zero(self):
        assert lc.is_calibrated_night(1000, [995], [1010])
        assert not lc.is_calibrated_night(1000, [995], [])
        assert not lc.is_calibrated_night(1000, [], [1010])

    def test_flat_pairs_cluster_by_exposure(self):
        # 3 s x4 -> 2 pairs; 5 s x2 -> 1 pair; 7 s x1 -> none.
        assert lc.count_flat_pairs([3.0, 3.0, 3.001, 3.0, 5.0, 5.0, 7.0]) \
            == (3, 2)

    def test_flat_pairs_ignore_null_and_zero(self):
        assert lc.count_flat_pairs([None, 0.0, 2.0]) == (0, 0)
        assert lc.count_flat_pairs([]) == (0, 0)


# ---------------------------------------------------------------------------
# 8. Eclipsers
# ---------------------------------------------------------------------------
class TestEclipsers:
    @pytest.mark.parametrize("vtype", ["EW", "EW/KW", "EA+BY", "EB:", "E",
                                       "UGSU+E", "EA/SD|EB"])
    def test_eclipsing_types(self, vtype):
        assert lc.is_eclipsing_type(vtype)

    @pytest.mark.parametrize("vtype", ["EP", "ELL", "ROT", "UGSU", "AM",
                                       "EXOR", "DSCT", None, ""])
    def test_not_eclipsing_types(self, vtype):
        assert not lc.is_eclipsing_type(vtype)

    def test_run_of_half_a_period_guarantees_a_minimum(self):
        assert lc.run_guarantees_minimum(0.17, 0.3336)      # W UMa: P/2=0.167
        assert not lc.run_guarantees_minimum(0.16, 0.3336)

    def test_no_period_needs_three_hours(self):
        assert lc.run_guarantees_minimum(3.0 / 24, None)
        assert not lc.run_guarantees_minimum(2.9 / 24, None)

    def test_predicted_minimum_needs_margin_each_side(self):
        epoch, period = 2457000.0, 0.4
        # Minima at epoch + k*0.2.  Run brackets the one at +0.2 by 45 min.
        lo, hi = epoch + 0.2 - 45 / 1440, epoch + 0.2 + 45 / 1440
        assert lc.predicted_minima_in_run(lo, hi, epoch, period) \
            == pytest.approx([epoch + 0.2])
        # Only 20 min after the minimum: margin not met.
        hi_short = epoch + 0.2 + 20 / 1440
        assert lc.predicted_minima_in_run(lo, hi_short, epoch, period) == []

    def test_predicted_counts_primary_and_secondary(self):
        epoch, period = 2457000.0, 0.4
        got = lc.predicted_minima_in_run(epoch - 0.05, epoch + 0.45, epoch,
                                         period)
        assert got == pytest.approx([epoch, epoch + 0.2, epoch + 0.4])

    def test_no_ephemeris_no_prediction(self):
        assert lc.predicted_minima_in_run(1.0, 2.0, None, 0.4) == []
        assert lc.predicted_minima_in_run(1.0, 2.0, 0.5, None) == []


# ---------------------------------------------------------------------------
# 9. Time-convention audit
# ---------------------------------------------------------------------------
class TestTimeAudit:
    def test_card_agreement(self):
        n, good, frac = lc.card_agreement([0.0, 0.5, -1.9, 2.1, 30.0])
        assert (n, good) == (5, 3) and frac == pytest.approx(0.6)

    def _sequence(self, stamp):
        """Alternating 60 s / 10 s exposures, 5 s overhead, stamped at
        start / mid / end of each exposure."""
        frames, t = [], 0.0
        for i in range(40):
            e = 60.0 if i % 2 == 0 else 10.0
            mark = {"start": t, "mid": t + e / 2, "end": t + e}[stamp]
            frames.append((2457000.5 + mark / 86400.0, e))
            t += e + 5.0
        return frames

    def test_start_stamps_violate_only_the_other_hypotheses(self):
        v = lc.stamp_hypothesis_violations(self._sequence("start"))
        assert v["informative"] == 39
        assert v["viol_start"] == 0 and v["viol_end"] > 0 and v["viol_mid"] > 0

    def test_end_stamps_violate_the_start_hypothesis(self):
        v = lc.stamp_hypothesis_violations(self._sequence("end"))
        assert v["viol_end"] == 0 and v["viol_start"] > 0

    def test_equal_exposures_are_uninformative(self):
        frames = [(2457000.5 + i * 35 / 86400.0, 30.0) for i in range(20)]
        v = lc.stamp_hypothesis_violations(frames)
        assert v["pairs"] == 19 and v["informative"] == 0

    def test_overlapping_exposures_are_counted(self):
        frames = [(2457000.5, 60.0), (2457000.5 + 20 / 86400.0, 60.0)]
        assert lc.stamp_hypothesis_violations(frames)["viol_start_all"] == 1

    def test_pass_requires_a_second_card(self):
        assert not lc.time_convention_pass(0, 0, 0, 0)
        assert lc.time_convention_pass(1000, 995, 0, 0)
        assert not lc.time_convention_pass(1000, 980, 0, 0)     # 98% < 99%
        assert not lc.time_convention_pass(1000, 1000, 1, 0)    # duplicate
        assert not lc.time_convention_pass(1000, 1000, 0, 1)    # overlap

    def test_convention_verdict(self):
        from types import SimpleNamespace as NS
        v = build.convention_verdict
        assert v(NS(n_informative_pairs=2000, viol_start=0, viol_mid=200,
                    viol_end=400)) == build.CONV_START
        assert v(NS(n_informative_pairs=2000, viol_start=400, viol_mid=200,
                    viol_end=0)) == build.CONV_END
        # The AC4040 signature: start AND end violated, middle never.
        assert v(NS(n_informative_pairs=2540, viol_start=580, viol_mid=0,
                    viol_end=690)) == build.CONV_MID
        assert v(NS(n_informative_pairs=0, viol_start=0, viol_mid=0,
                    viol_end=0)).startswith("not identifiable")
        assert v(NS(n_informative_pairs=50, viol_start=0, viol_mid=0,
                    viol_end=0)).startswith("not identifiable")
        assert v(NS(n_informative_pairs=2000, viol_start=300, viol_mid=300,
                    viol_end=300)).startswith("inconsistent")

    def test_mid_stamps_are_recognised_end_to_end(self):
        v = lc.stamp_hypothesis_violations(self._sequence("mid"))
        assert v["viol_mid"] == 0 and v["viol_start"] > 0 and v["viol_end"] > 0

    @pytest.mark.parametrize("stamp, slope", [("start", 0.0), ("mid", -0.5),
                                              ("end", -1.0)])
    def test_lst_slope_separates_the_conventions(self, stamp, slope):
        # LST is sampled as the exposure begins; DATE-OBS is the stamp.
        import numpy as np
        rng = np.random.default_rng(1)
        e = rng.choice([16.0, 64.0, 256.0], 400)
        frac = {"start": 0.0, "mid": 0.5, "end": 1.0}[stamp]
        resid = -frac * e - 2.0 + rng.normal(0, 0.3, 400)
        resid[:5] = 15000.0                       # stale LST cards
        got, icpt, n = build.lst_exposure_slope(e, resid)
        assert got == pytest.approx(slope, abs=0.01)
        assert icpt == pytest.approx(-2.0, abs=0.2) and n >= 390

    def test_lst_slope_needs_more_than_one_exposure_time(self):
        import math
        s, _i, _n = build.lst_exposure_slope([60.0] * 50, [0.0] * 50)
        assert math.isnan(s)


# ---------------------------------------------------------------------------
# 10. Gates and decision rule — pinned to the pre-registration note
# ---------------------------------------------------------------------------
class TestGates:
    def test_thresholds_match_the_preregistration(self):
        note = (Path(__file__).resolve().parents[2] / "Legacy_Rigel" / "notes"
                / "GO_NOGO_PREREGISTERED.md").read_text()
        assert "≥ 3 eclipsing/contact systems each with minimum-bearing " \
               "nights in\n  ≥ 3 seasons" in note
        assert (lc.G1_MIN_SYSTEMS, lc.G1_MIN_SEASONS) == (3, 3)
        assert "≥ 10 nights in one series, spread over\n  ≥ 2 seasons" in note
        assert (lc.G2_MIN_NIGHTS, lc.G2_MIN_SEASONS) == (10, 2)
        assert "run ≥ 1 h is also required" in note
        assert lc.G2_TIMING_RUN_DAYS == pytest.approx(1 / 24)
        assert "≥ 30 calibrated nights in one series" in note
        assert lc.G3_MIN_CALIBRATED_NIGHTS == 30
        assert "no gap longer than\n  30 min" in note
        assert "gaps > 120 d" in note and "±30 nights" in note
        assert "< 2 s on ≥ 99%" in note
        assert (lc.TIME_AGREE_SECONDS, lc.TIME_AGREE_FRACTION) == (2.0, 0.99)

    def test_g1(self):
        assert lc.g1_passes([3, 4, 5])
        assert not lc.g1_passes([3, 4, 2, 2, 1])
        assert not lc.g1_passes([])

    def test_g2_series(self):
        ok = dict(n_nights=10, n_seasons=2, tieable=True,
                  is_timing_target=False, longest_run_days=0.0)
        assert lc.g2_series_passes(**ok)
        assert not lc.g2_series_passes(**{**ok, "n_nights": 9})
        assert not lc.g2_series_passes(**{**ok, "n_seasons": 1})
        assert not lc.g2_series_passes(**{**ok, "tieable": False})

    def test_g2_timing_target_needs_an_hour(self):
        base = dict(n_nights=12, n_seasons=3, tieable=True,
                    is_timing_target=True)
        assert not lc.g2_series_passes(**base, longest_run_days=0.5 / 24)
        assert lc.g2_series_passes(**base, longest_run_days=1.0 / 24)

    def test_g3_series(self):
        assert lc.g3_series_passes(30) and not lc.g3_series_passes(29)

    def test_g0_failure_overrides_everything(self):
        assert lc.decide(False, True, True, True, True) \
            == [lc.OUTCOME_NOT_DECIDABLE]

    def test_nothing_passes_is_release_only(self):
        assert lc.decide(True, False, False, False) \
            == [lc.OUTCOME_RELEASE_ONLY]

    def test_g2_alone_is_a_transfer_not_a_go(self):
        out = lc.decide(True, False, True, False)
        assert out == [lc.OUTCOME_TRANSFER]
        assert lc.OUTCOME_GO_CANDIDATE not in out

    def test_g1_is_a_go_candidate(self):
        assert lc.decide(True, True, False, False) \
            == [lc.OUTCOME_GO_CANDIDATE]

    def test_g3_without_a_named_question_is_not_a_go(self):
        assert lc.decide(True, False, False, True, g3_question_named=False) \
            == [lc.OUTCOME_RELEASE_ONLY]
        assert lc.decide(True, False, False, True, g3_question_named=True) \
            == [lc.OUTCOME_GO_CANDIDATE]

    def test_transfer_and_go_can_both_apply(self):
        assert lc.decide(True, True, True, False) \
            == [lc.OUTCOME_TRANSFER, lc.OUTCOME_GO_CANDIDATE]

    def test_g3_named_questions_is_empty_by_design(self):
        # Naming a question for a G3 target is a human decision; the
        # shipped table must not pre-empt it.
        assert build.G3_NAMED_QUESTIONS == {}


# ---------------------------------------------------------------------------
# The build script's DataFrame paths, end to end on a synthetic archive
# ---------------------------------------------------------------------------
def _scan_row(path, date_obs, **kw):
    """One synthetic scan row with every column the build reads."""
    row = dict.fromkeys(build.FRAME_SCAN_COLUMNS)
    row.update(path=path, date_obs=date_obs, telescop="Iowa Robotic Telescope",
               instrume="Andor IKON L 936", swcreate="MaxIm DL Version 6.11 X",
               imagetyp="Light Frame", filter="G - Sloan_g", object="W UMa",
               exptime=30.0, xbinning=1.0, xpixsz=13.5, naxis1=2048,
               naxis2=2048, crota2=1.0, lst="10:00:00",
               objra="09:43:45.5", objdec="55:57:09")
    row.update(kw)
    jd = lc.parse_date_obs(date_obs)
    if row["jd"] is None:
        row["jd"] = jd
    return row


@pytest.fixture()
def synthetic_frames():
    rows = []
    # Three seasons of W UMa, each with one 5-hour night at 2-min cadence
    # (P/2 = 4.0 h -> minimum guaranteed) and one short night.
    for year in (2017, 2018, 2019):
        for k in range(150):
            m = 3 * 60 + 2 * k
            rows.append(_scan_row(
                f"{year}/day060/gaa060{k:02x}.fts.fz",
                f"{year}-03-01T{m // 60:02d}:{m % 60:02d}:00"))
        for k in range(5):
            rows.append(_scan_row(
                f"{year}/day070/gaa070{k:02x}.fts.fz",
                f"{year}-03-11T04:{2 * k:02d}:00"))
    # A duplicate copy of one frame in another tree, calibrated in place.
    rows.append(_scan_row("old/day060/gaa06000.fts.fz", "2017-03-01T03:00:00",
                          calstat="BDF"))
    # A focus frame, an unset-clock frame, an unreadable file.
    rows.append(_scan_row("2017/day060/foc06001.fts.fz", "2017-03-01T02:00:00",
                          object="Vega"))
    rows.append(_scan_row("2017/day061/gaa06100.fts.fz",
                          "1970-01-01T00:00:02.670", jd=2440587.5))
    bad = dict.fromkeys(build.FRAME_SCAN_COLUMNS)
    bad.update(path="2017/day062/bad.fts.fz", error="OSError: truncated")
    rows.append(bad)
    # Flats and a bias on the 2017 long night.
    for k in range(4):
        rows.append(_scan_row(f"2017/day060/flt060{k:02x}.fts.fz",
                              f"2017-03-01T01:{k:02d}:00",
                              imagetyp="Flat Field", object="", exptime=3.0))
    rows.append(_scan_row("2017/day060/bia06000.fts.fz", "2017-03-01T01:30:00",
                          imagetyp="Bias Frame", object="", exptime=0.0))
    return build.derive_frames(pd.DataFrame(rows))


class TestBuildFrames:
    def test_every_row_kept_and_dispositioned(self, synthetic_frames):
        f = synthetic_frames
        counts = f["exclusion"].fillna("science").value_counts().to_dict()
        assert counts == {"science": 465, "duplicate_copy": 1,
                          "not_science:focus": 1, "no_usable_date_obs": 1,
                          "header_unreadable": 1, "not_science:flat": 4,
                          "not_science:bias": 1}
        assert len(f) == sum(counts.values())

    def test_raw_copy_is_canonical_over_calibrated_twin(self, synthetic_frames):
        f = synthetic_frames
        grp = f[f["date_obs"] == "2017-03-01T03:00:00"]
        grp = grp[grp["kind"] == lc.KIND_LIGHT]
        assert len(grp) == 2
        canon = grp[grp["is_canonical"] == 1]["path"].iloc[0]
        assert canon == "2017/day060/gaa06000.fts.fz"

    def test_series_runs_seasons_and_calibration(self, synthetic_frames):
        f = synthetic_frames
        sci = f[f["is_science"] == 1]
        seasons = {k: lc.assign_seasons(g["night"])
                   for k, g in sci.groupby("target_key")}
        flats, zeros = build.calibration_index(f)
        runs = build.build_runs(sci)
        series = build.build_series(sci, runs, flats, zeros, seasons)
        assert len(series) == 1
        s = series.iloc[0]
        assert (s["n_nights"], s["n_seasons"]) == (6, 3)
        assert s["longest_run_h"] == pytest.approx(149 * 2 / 60)
        assert s["n_nights_run_3h"] == 3
        # Flats + bias exist only around the 2017 nights (within 30 nights).
        assert s["n_calibrated_nights"] == 2
        assert s["tieable"] == 1

    def test_flat_pairs_table(self, synthetic_frames):
        p = build.build_flat_pairs(synthetic_frames)
        assert len(p) == 1 and p["n_pairs"].iloc[0] == 2

    def test_stints_and_time_audit(self, synthetic_frames):
        f = synthetic_frames
        stints = build.build_stints(f)
        assert list(stints["camera"]) == ["Andor iKon-L 936"]
        audit = build.build_time_audit(f)
        a = audit.iloc[0]
        assert a["n_no_usable_date"] == 1
        assert a["frac_jd_agree"] == 1.0 and a["card_pass"] == 1
        # 30-s science, 3-s flats and 0-s bias give unequal pairs; all
        # stamps are starts, so the start hypothesis is never violated.
        assert a["viol_start"] == 0

    def test_camera_table_reports_identity_basis(self, synthetic_frames):
        cams = build.build_cameras(synthetic_frames)
        assert cams["camera"].iloc[0] == "Andor iKon-L 936"
        assert "label" in cams["identity_basis"].iloc[0]
        assert cams["n_label_stale"].iloc[0] == 0


class TestBuildReconciliation:
    def _db(self, tmp_path):
        con = sqlite3.connect(tmp_path / "t.sqlite")
        scanner.ensure_scan_tables(con)
        scanner.write_disk_files(con, [
            ("2015/day004/a.fts.fz", "fz"), ("2015/day004/a.fts", "fits"),
            ("2015/day004/b.fts.fz", "fz"), ("2015/day004/c.fts", "fits"),
            ("2015/day004/extra.fts.fz", "fz"), ("log.txt", "other")])
        con.executemany("INSERT INTO scan (path, error) VALUES (?, ?)", [
            ("2015/day004/a.fts.fz", None), ("2015/day004/b.fts.fz", "boom"),
            ("2015/day004/extra.fts.fz", None)])
        con.commit()
        return con

    def test_both_identities_close(self, tmp_path):
        con = self._db(tmp_path)
        manifest = [("1", "2015/day004/a.fts"), ("2", "2015/day004/b.fts"),
                    ("3", "2015/day004/c.fts"), ("4", "2015/day004/gone.fts")]
        rec = build.build_reconciliation(con, manifest)
        status = dict(zip(rec["logical_path"], rec["status"]))
        assert status == {
            "2015/day004/a.fts": lc.REC_SCANNED,
            "2015/day004/b.fts": lc.REC_UNREADABLE,
            "2015/day004/c.fts": lc.REC_UNCOMPRESSED_ONLY,
            "2015/day004/gone.fts": lc.REC_NOT_ON_DISK,
            "2015/day004/extra.fts": lc.REC_SCANNED,
            "log.txt": lc.REC_NOT_FITS}
        s = dict(zip(*[build.reconciliation_summary(rec, con)[c]
                       for c in ("quantity", "n")]))
        assert s["FILE IDENTITY residual (must be 0)"] == 0
        assert s["MANIFEST IDENTITY residual (must be 0)"] == 0
        assert s["scan rows minus .fz files on disk (must be 0)"] == 0
        assert s["  uncompressed twin of a scanned .fz (named exclusion)"] == 1
        assert s["on disk but not in the manifest (FITS)"] == 1

    def test_collision_audit_proves_distinct_exposures(self):
        good = [("1", "2015/day004/a.fts"), ("2", "2017/day004/a.fts"),
                ("3", "2016/day009/z.fts")]
        bad = [("1", "day004/a.fts"), ("2", "day004/a.fts"),
               ("3", "day009/z.fts")]
        frames = pd.DataFrame({
            "logical_path": ["2015/day004/a.fts", "2017/day004/a.fts",
                             "2016/day009/z.fts"],
            "date_obs": ["2015-01-04T05:00:00", "2017-01-04T05:00:00",
                         "2016-01-09T05:00:00"],
            "error": [None, None, None]})
        audit, meta = build.build_collision_audit(good, bad, frames)
        assert meta["collision_groups"] == 1
        assert meta["collision_paths_at_risk"] == 1
        assert meta["collision_ids_only_in_good"] == 0
        assert meta["collision_bad_distinct_paths"] == 2
        assert meta["collision_groups_all_distinct"] == 1
        assert audit["n_distinct_date_obs"].iloc[0] == 2


class TestSoftwareFamily:
    def test_maxim_version_without_serial(self):
        assert build.software_family("MaxIm DL Version 6.08 150819 2UT5A",
                                     "Gemini") == "MaxIm DL 6.08"

    def test_talon_inferred_from_rigel(self):
        assert build.software_family(None, "Rigel System") \
            == "Talon (no SWCREATE)"

    def test_unknown(self):
        assert build.software_family(None, "Gemini") == "(no SWCREATE)"


class TestVsxMatch:
    def test_name_match_wins_over_a_nearer_stranger(self):
        cands = [dict(Name="Gaia DR3 1", sep_arcsec=1.0),
                 dict(Name="W UMa", sep_arcsec=4.0)]
        best, how = external.choose_vsx_match("wuma", cands)
        assert best["Name"] == "W UMa" and how == "name"

    def test_position_match_needs_to_be_close(self):
        near = [dict(Name="ZTF J1", sep_arcsec=8.0)]
        far = [dict(Name="ZTF J1", sep_arcsec=25.0)]
        assert external.choose_vsx_match("field7", near)[1] == "position"
        assert external.choose_vsx_match("field7", far) == (None, None)

    def test_no_candidates(self):
        assert external.choose_vsx_match("m31", []) == (None, None)
