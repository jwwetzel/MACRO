"""Unit tests for macro_core.manifest — every pure S0 function.

Each test class mirrors one decision the manifest encodes; the must-NOT
cases (distinct targets staying separate, non-duplicates staying apart) are
as important as the must cases, because a silent false merge poisons every
downstream count.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import math
import sqlite3
import sys
from contextlib import closing
from pathlib import Path

import pytest

# Make the package importable regardless of pytest's working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
# The build script itself is also under test (its vectorized pandas paths
# must agree with the pure functions — the era_id .map bug shipped exactly
# because only the pure layer was tested).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import pandas as pd

from macro_core import manifest as m

import build_s0_manifest as build


# ---------------------------------------------------------------------------
# Night labels: local-noon-to-noon boundary (Winer, UTC−7 → 19:00 UT)
# ---------------------------------------------------------------------------
class TestNightLabel:
    # JD of 2023-06-01T12:00:00 UT (JD integers fall at UT noon).
    JD_NOON_UT = 2460097.0

    def test_before_local_noon_belongs_to_previous_night(self):
        # 18:00 UT = 11:00 local — still the *previous* observing night.
        jd = self.JD_NOON_UT + 6.0 / 24.0
        assert m.night_label(jd) == "2023-05-31"

    def test_after_local_noon_starts_the_new_night(self):
        # 20:00 UT = 13:00 local — the new night has begun.
        jd = self.JD_NOON_UT + 8.0 / 24.0
        assert m.night_label(jd) == "2023-06-01"

    def test_evening_and_dawn_share_one_label(self):
        # A real observing sequence: 02:00 UT (19:00 local, evening) through
        # 12:00 UT (05:00 local, dawn) of the next UTC day must be ONE night
        # labeled with the local evening date.
        evening = 2460098.0 - 10.0 / 24.0     # 2023-06-02T02:00 UT
        dawn = 2460098.0                      # 2023-06-02T12:00 UT
        assert m.night_label(evening) == "2023-06-01"
        assert m.night_label(dawn) == "2023-06-01"

    def test_boundary_is_at_the_configured_shift(self):
        # Immediately either side of the exact boundary JD fraction: the
        # label must flip.  Boundary sits at (integer + 0.7917) − 0.5 in
        # shifted space; test one second on each side.
        boundary = self.JD_NOON_UT + (m.NIGHT_SHIFT_DAYS - 0.5)
        second = 1.0 / 86400.0
        assert m.night_label(boundary - second) == "2023-05-31"
        assert m.night_label(boundary + second) == "2023-06-01"

    def test_missing_jd_gives_no_label(self):
        assert m.night_label(None) is None
        assert m.night_label(float("nan")) is None


# ---------------------------------------------------------------------------
# Alias normalization: every rule, and the merges that MUST NOT happen
# ---------------------------------------------------------------------------
class TestNormalizeTarget:
    def key(self, raw):
        return m.normalize_target(raw).key

    # --- case folding / whitespace ---------------------------------------
    def test_case_variants_merge(self):
        assert self.key("T CrB") == self.key("t crb") == "tcrb"
        assert self.key("ST LMi") == self.key("ST LMI") == self.key("st lmi")

    def test_whitespace_and_leading_spaces(self):
        assert self.key("  eta UMa") == self.key("eta UMa")
        assert self.key("RR  LYR") == self.key("RR LYR")

    # --- separator unification -------------------------------------------
    def test_space_and_hyphen_variants_merge(self):
        assert self.key("YZCnc") == self.key("YZ Cnc") == "yzcnc"
        assert self.key("5cnc") == self.key("5 Cnc") == "5cnc"
        assert self.key("NUUMa") == self.key("NU UMa")
        assert self.key("ngc5548") == self.key("NGC 5548") == "ngc5548"
        assert self.key("lameri") == self.key("lam eri") == self.key("Lam Eri")

    def test_hyphen_before_digit_is_preserved(self):
        # The hyphen in a ZTF name is a declination SIGN — collapsing it
        # would alias a southern object onto a northern one.
        assert "-" in self.key("ZTFJ082835 05-052702 1")

    # --- exposure-token leakage ------------------------------------------
    def test_filter_then_exposure_leakage_stripped(self):
        assert self.key("PHECDA lrg 0-25s") == self.key("PHECDA") == "phecda"
        assert self.key("tet CrB hrg 2-4s") == self.key("tet CrB")
        assert self.key("W Uma g 0-12s") == self.key("W Uma")
        assert self.key("10 Cas r 0-5s") == self.key("10 Cas")

    def test_exposure_then_filter_then_index_leakage_stripped(self):
        # The Vega ladder order: '<name> <exp> <filter> [<frame index>]'.
        assert self.key("Vega 0p001s lrg 5") == "vega"
        assert self.key("Vega 0p1s hrg") == "vega"
        assert self.key("Vega 0p0001s lrg 7") == "vega"

    def test_rules_recorded_for_audit(self):
        info = m.normalize_target("PHECDA lrg 0-25s")
        assert "exposure_tokens" in info.rules

    def test_bare_trailing_digit_without_exposure_is_kept(self):
        # 'ZTFJ104433 60+292809 5': the trailing 5 is a decimal of the
        # declination, and there is NO exposure token — nothing may strip.
        assert self.key("ZTFJ104433 60+292809 5").endswith("5")

    def test_short_final_words_survive(self):
        # Names whose last token merely LOOKS like a filter letter must not
        # lose it: no exposure token, no strip.
        assert self.key("ksi UMa B") == "ksiumab"
        assert self.key("AG LMi") == "aglmi"

    # --- prefixes and suffixes -------------------------------------------
    def test_simbad_prefixes_stripped(self):
        assert self.key("* tet CrB") == self.key("tet CrB") == "tetcrb"
        assert self.key("V* KV UMa") == self.key("KV UMa")

    def test_observer_date_prefix_stripped(self):
        assert self.key("mjcMay01 yzcnc") == "yzcnc"
        assert self.key("mjcApr30 yzcnc") == "yzcnc"
        assert self.key("mjcMay02 TIC 48217457") == self.key("TIC 48217457")

    def test_series_suffix_stripped(self):
        assert self.key("ST-LMi-y-series") == "stlmi"
        assert self.key("STLMi-z-series") == "stlmi"

    def test_genitive_normalized(self):
        assert self.key("RR Lyrae") == self.key("RR Lyr") == "rrlyr"

    # --- synonym table (cone-gated at build time) ------------------------
    def test_synonym_maps_and_records_pre_key(self):
        info = m.normalize_target("Alpha Lyr")
        assert info.key == "vega"
        assert info.pre_synonym_key == "alphalyr"
        assert "synonym" in info.rules

    def test_pinwheel_galaxy_merges_into_m101(self):
        """The observer's common name for M101, 245 catalog rows.

        140 of them are canonical rawimage science on 2026-03-21/22 — the
        deepest post-fade M101-field epoch in the archive, and invisible to
        the SN project until this merge.  No string rule can relate a common
        name to a Messier number, so it is a synonym-table entry; the build
        cone-gates it (measured 0.031 deg from the natives of 'm101').
        """
        info = m.normalize_target("Pinwheel Galaxy")
        assert info.key == "m101"
        assert info.pre_synonym_key == "pinwheelgalaxy"
        assert "synonym" in info.rules
        assert self.key("pinwheel galaxy") == self.key("M101")

    def test_sequence_digit_sn_names_merge_into_the_sn(self):
        """'2023ixf1'/'2023ixf2': a trailing sequence counter fused to the
        object name deleted a NAMED template epoch (2023-11-28, 64 s L)
        from the SN working set.  Both sit within 0.05 deg of the SN."""
        for raw in ("2023ixf1", "2023ixf2"):
            info = m.normalize_target(raw)
            assert info.key == "2023ixf", raw
            assert "synonym" in info.rules, raw

    # --- the merges that must NOT happen ---------------------------------
    def test_distinct_dwarf_fields_stay_separate(self):
        # The spec's named regression: two different survey fields.
        assert self.key("Dw1403+49") != self.key("Dw1409+51")

    def test_sn_and_calibration_names_stay_themselves(self):
        assert self.key("2023ixf") == "2023ixf"
        assert self.key("2023ixf") != self.key("M101")

    def test_sequence_digit_merge_is_explicit_not_a_digit_rule(self):
        """The counter merge is a TABLE entry, not a "strip trailing digits"
        rule — such a rule would need to know which stems are real target
        keys, which this pure normalizer deliberately cannot know, and it
        would happily eat the survey-field and HR-number names below."""
        assert self.key("Dw1403+49") == "dw1403+49"
        assert self.key("HR 3454") == "hr3454"
        assert self.key("2023ixf3") == "2023ixf3"   # not in the table

    def test_tcrb_and_thetacrb_stay_separate(self):
        # Target and its spectrophotometric calibrator: 8 degrees apart.
        assert self.key("T CrB") != self.key("tet CrB")

    def test_blank_names_form_no_group(self):
        for raw in (None, "", "   "):
            info = m.normalize_target(raw)
            assert info.key is None
            assert info.rules == ("blank",)


# ---------------------------------------------------------------------------
# Duplicate identity and canonical selection (tree policy)
# ---------------------------------------------------------------------------
class TestDedup:
    def test_same_basename_and_jd_is_one_group(self):
        assert m.dup_key("a.fts", 2460100.5) == m.dup_key("a.fts", 2460100.5)

    def test_different_jd_is_a_different_group(self):
        assert m.dup_key("a.fts", 2460100.5) != m.dup_key("a.fts", 2460101.5)

    def test_missing_jd_never_merges(self):
        # Unreadable files cannot be proven duplicates: rowid keeps them
        # apart even with identical basenames.
        assert m.dup_key("a.fts", None, row_id=1) != \
               m.dup_key("a.fts", None, row_id=2)
        assert m.dup_key("a.fts", float("nan"), row_id=3) != \
               m.dup_key("a.fts", float("nan"), row_id=4)

    def test_rawimage_beats_every_other_tree(self):
        members = [("reduced", "reduced/x/a.fts"),
                   ("external", "external/x/a.fts"),
                   ("rawimage", "rawimage/2023-05-21/a.fts")]
        assert m.choose_canonical(members) == 2

    def test_reduced_is_the_last_resort(self):
        members = [("reduced", "reduced/x/a.fts"),
                   ("macalester", "macalester/x/a.fts")]
        assert m.choose_canonical(members) == 1

    def test_within_tree_earliest_path_wins(self):
        # The SN 2023ixf failure mode: the same frame wholesale-copied into
        # a later night directory of the SAME tree.  The original (earlier,
        # lexicographically smaller) night directory must win.
        members = [("rawimage", "rawimage/2023-07-07/xwg1411c.fts.fz"),
                   ("rawimage", "rawimage/2023-05-21/xwg1411c.fts.fz")]
        assert m.choose_canonical(members) == 1

    def test_tree_policy_exception_ngc5548(self):
        # Documented exception: NGC 5548's canonical copy lives under
        # macalester (the superset tree).  Same members, no exception key →
        # rawimage; with the NGC 5548 key → macalester.
        members = [("rawimage", "rawimage/2023-03-24/n5548a.fts"),
                   ("macalester", "macalester/x/n5548a.fts")]
        assert m.choose_canonical(members) == 0
        assert m.choose_canonical(members, target_key="ngc5548") == 1

    def test_unknown_tree_still_canonical_when_alone(self):
        assert m.choose_canonical([("mystery", "mystery/a.fts")]) == 0

    def test_basename_extraction(self):
        assert m.basename_of("rawimage/2023-05-21/xwg1411c.fts.fz") == \
            "xwg1411c.fts.fz"
        assert m.basename_of("bare.fts") == "bare.fts"


# ---------------------------------------------------------------------------
# Era keys: stability against header fuzz, sensitivity to real changes
# ---------------------------------------------------------------------------
class TestEraKey:
    def test_float_fuzz_and_whitespace_do_not_split_an_era(self):
        a = m.era_key("High Gain", 4096.0, 4096.0, 1.0, 1.05697000026703)
        b = m.era_key("High Gain ", 4096, 4096, 1, 1.0570)
        assert a == b

    def test_distinct_gains_are_distinct_eras(self):
        # The two GSENSE driver epochs differ in the 3rd decimal: real.
        a = m.era_key("High Gain", 4096, 4096, 1, 1.05697000026703)
        b = m.era_key("High Gain", 4096, 4096, 1, 1.05371475219727)
        assert a != b

    def test_geometry_and_binning_split_eras(self):
        base = m.era_key("Mode0", 4788, 3194, 2, 0.2467)
        assert base != m.era_key("Mode0", 9576, 6388, 1, 0.2467)
        assert base != m.era_key("Mode0", 4788, 3194, 1, 0.2467)

    def test_missing_values_are_a_key_not_a_crash(self):
        # The pyscope-transition rows with blank READOUTM/EGAIN still form
        # a (single, well-defined) era.
        a = m.era_key(None, 8.0, 3211.0, 2.0, None)
        b = m.era_key("", 8, 3211, 2, float("nan"))
        assert a == b

    def test_filter_and_date_are_not_in_the_key(self):
        # By construction the signature takes neither — this test documents
        # the roadmap rule for future editors.
        import inspect
        params = inspect.signature(m.era_key).parameters
        assert "filter" not in params and "date_obs" not in params


# ---------------------------------------------------------------------------
# Pointing math on known angles
# ---------------------------------------------------------------------------
class TestPointingMath:
    def test_zero_separation(self):
        assert m.angular_separation_deg(123.4, -56.7, 123.4, -56.7) == \
            pytest.approx(0.0, abs=1e-12)

    def test_quarter_circle_on_equator_and_meridian(self):
        assert m.angular_separation_deg(0, 0, 90, 0) == pytest.approx(90.0)
        assert m.angular_separation_deg(0, 0, 0, 90) == pytest.approx(90.0)

    def test_ra_wraparound(self):
        # 0.2 degrees apart across the 0h seam — NOT 359.8 degrees.
        assert m.angular_separation_deg(359.9, 0, 0.1, 0) == \
            pytest.approx(0.2, abs=1e-9)

    def test_pole_to_pole(self):
        assert m.angular_separation_deg(10, 90, 200, -90) == \
            pytest.approx(180.0)

    def test_high_declination_ra_offsets_shrink(self):
        # 1 degree of RA at dec 60 is ~0.5 degrees on the sky.
        sep = m.angular_separation_deg(10, 60, 11, 60)
        assert sep == pytest.approx(math.degrees(
            math.acos(math.sin(math.radians(60))**2
                      + math.cos(math.radians(60))**2
                      * math.cos(math.radians(1)))), rel=1e-9)

    def test_median_radec_plain(self):
        ra, dec = m.median_radec([10, 11, 12], [5, 6, 7])
        assert (ra, dec) == (11, 6)

    def test_median_radec_wraps(self):
        # Samples straddling 0h must average near 0, never near 180.
        ra, dec = m.median_radec([359.9, 0.1, 0.0], [1, 2, 3])
        assert min(ra, 360 - ra) < 0.2
        assert dec == 2

    def test_median_radec_empty_raises(self):
        with pytest.raises(ValueError):
            m.median_radec([], [])


# ---------------------------------------------------------------------------
# The pointing reference, and the target it used to be blind on
# ---------------------------------------------------------------------------
class TestPointingReference:
    def test_prefers_plate_solved_coordinates(self):
        # Three solved frames at dec 25.1 and two unsolved strays at 33.3:
        # the reference must come from the solves alone.
        ra0, dec0, basis = m.pointing_reference(
            [214.50, 214.51, 214.52, 209.43, 222.52],
            [25.12, 25.13, 25.14, 33.35, 33.33],
            [1, 1, 1, None, None])
        assert basis == m.POINTING_REF_SOLVED
        assert dec0 == pytest.approx(25.13)
        assert ra0 == pytest.approx(214.51)

    def test_falls_back_to_header_median_when_nothing_is_solved(self):
        """THE NGC 5548 REGRESSION, on the real shape of that target.

        279 frames, ZERO plate-solved, so the old solved-only rule produced
        no reference, no offsets and no QC flags — on precisely the target
        whose strategy rules a whole night 'unusable, pointing-failure'.
        The mispointed night is the mount parked and not tracking: dec
        +33.3 while the target sits at +25.1, RA marching at the sidereal
        rate.  A median over the frames is robust to that minority.
        """
        on_target_ra = [214.50] * 9
        on_target_dec = [25.12] * 9
        # the parked-mount night: RA marches, dec is 8 degrees off
        bad_ra = [209.43, 213.0, 218.0, 222.52]
        bad_dec = [33.35, 33.34, 33.34, 33.33]
        ra0, dec0, basis = m.pointing_reference(
            on_target_ra + bad_ra, on_target_dec + bad_dec,
            [None] * 13)
        assert basis == m.POINTING_REF_HEADER
        assert dec0 == pytest.approx(25.12)
        # …and every mispointed frame now clears the outlier threshold.
        for ra, dec in zip(bad_ra, bad_dec):
            off = m.angular_separation_deg(ra, dec, ra0, dec0)
            assert off > m.POINTING_OUTLIER_DEG
            assert "pointing_gt1deg" in m.qc_flags(
                None, 256.0, 1.07, 2460000.5, ra, "ngc5548", off)

    def test_header_reference_is_refused_when_headers_disagree(self):
        """THE MIZAR CASE — the fallback's own failure mode, closed.

        30 canonical frames: 10 with no coordinates, 10 at Mizar's true
        position, 10 carrying a stale RA/Dec card (seen on seven other
        targets too).  A bare median lands between the two clusters and
        flags BOTH — including the frames that were on target.  With no
        external truth to break the tie, the honest answer is no reference
        at all, recorded as such.
        """
        real = [(201.04, 54.96)] * 10
        stale = [(76.37, 52.84)] * 10
        ras = [p[0] for p in real + stale] + [None] * 10
        decs = [p[1] for p in real + stale] + [None] * 10
        ra0, dec0, basis = m.pointing_reference(ras, decs, [None] * 30)
        assert basis == m.POINTING_REF_UNSUPPORTED
        assert (ra0, dec0) == (None, None)

    def test_a_clear_majority_still_carries_the_reference(self):
        # NGC 5548's real shape: 269 on target, 10 parked — the support
        # test must NOT refuse this one.
        ras = [214.50] * 269 + [209.43] * 10
        decs = [25.12] * 269 + [33.35] * 10
        ra0, dec0, basis = m.pointing_reference(ras, decs, [None] * 279)
        assert basis == m.POINTING_REF_HEADER
        assert dec0 == pytest.approx(25.12)

    def test_support_test_does_not_touch_plate_solved_references(self):
        # Trusted evidence keeps its original behaviour even when the
        # solved frames themselves scatter beyond the outlier threshold.
        ra0, dec0, basis = m.pointing_reference(
            [10.0, 40.0, 200.0], [5.0, -5.0, 40.0], [1, 1, 1])
        assert basis == m.POINTING_REF_SOLVED
        assert ra0 is not None

    def test_no_reference_when_too_few_coordinates(self):
        ra0, dec0, basis = m.pointing_reference(
            [10.0, 10.1], [5.0, 5.1], [1, 1])
        assert (ra0, dec0, basis) == (None, None, m.POINTING_REF_NONE)

    def test_missing_and_nan_coordinates_are_skipped(self):
        nan = float("nan")
        ra0, dec0, basis = m.pointing_reference(
            [10.0, None, nan, 10.2, 10.1], [5.0, 5.0, nan, 5.2, 5.1],
            [None, None, None, None, None])
        assert basis == m.POINTING_REF_HEADER
        assert (ra0, dec0) == pytest.approx((10.1, 5.1))

    def test_solved_minority_still_wins_when_it_meets_the_minimum(self):
        # Evidence quality beats sample size: 3 solves outrank 50 headers.
        ra0, dec0, basis = m.pointing_reference(
            [10.0, 10.0, 10.0] + [200.0] * 50,
            [5.0, 5.0, 5.0] + [-40.0] * 50,
            [1, 1, 1] + [None] * 50)
        assert basis == m.POINTING_REF_SOLVED
        assert (ra0, dec0) == (10.0, 5.0)


# ---------------------------------------------------------------------------
# QC flags
# ---------------------------------------------------------------------------
class TestQcFlags:
    def clean(self, **over):
        """A frame with nothing wrong, overridable per test."""
        base = dict(error=None, exptime=60.0, airmass=1.2, jd=2460100.5,
                    ra_deg=150.0, target_key="tcrb",
                    pointing_offset_deg=0.01)
        base.update(over)
        return m.qc_flags(**base)

    def test_clean_frame_has_no_flags(self):
        assert self.clean() == ""

    def test_header_error(self):
        assert "header_error" in self.clean(error="OSError: Empty FITS")

    def test_exptime_nonpos(self):
        assert "exptime_nonpos" in self.clean(exptime=0.0)
        assert "exptime_nonpos" in self.clean(exptime=-1.0)
        assert "exptime_nonpos" not in self.clean(exptime=0.085)

    def test_airmass_garbage_low_high_and_sentinel(self):
        assert "airmass_garbage" in self.clean(airmass=0.99)
        assert "airmass_garbage" in self.clean(airmass=10.5)
        assert "airmass_garbage" in self.clean(airmass=-999.0)
        assert "airmass_garbage" not in self.clean(airmass=1.0)
        assert "airmass_garbage" not in self.clean(airmass=None)

    def test_missing_jd_and_coords_and_target(self):
        assert "no_jd" in self.clean(jd=None)
        assert "no_coords" in self.clean(ra_deg=None)
        assert "blank_target" in self.clean(target_key=None)

    def test_pointing_outlier_flag(self):
        assert "pointing_gt1deg" in self.clean(pointing_offset_deg=1.5)
        assert "pointing_gt1deg" not in self.clean(pointing_offset_deg=0.9)
        assert "pointing_gt1deg" not in self.clean(pointing_offset_deg=None)

    def test_flags_join_with_commas(self):
        flags = self.clean(exptime=0.0, airmass=-999.0)
        assert flags == "exptime_nonpos,airmass_garbage"


# ---------------------------------------------------------------------------
# Build-script integration: the vectorized pandas paths must agree with the
# tested pure functions.  These tests exist because the first shipped build
# had a bug (era_id NULL for every missing-EGAIN frame) that no pure-function
# test could catch — pandas' .map(dict) silently corrupts lookups of tuple
# keys containing None.  Everything here runs build_frames/resolve_aliases
# on small synthetic DataFrames, no catalog required.
# ---------------------------------------------------------------------------
def frame_row(**over):
    """One synthetic catalog row with nothing wrong; override per test."""
    base = dict(
        obs_rowid=1, path="rawimage/2024-01-01/a.fts", tree="rawimage",
        target_best="T CrB", jd=2460310.6, error=None,
        readoutm="High Gain", naxis1=4096.0, naxis2=4096.0, xbinning=1.0,
        egain=1.057, pltsolvd=0, ra_deg=None, dec_deg=None,
        exptime=60.0, airmass=1.2, imagetyp="Light Frame", filter="V",
    )
    base.update(over)
    return base


class TestBuildFramesPointingAudit:
    """Regression: the pointing audit must not go blind on an unsolved target.

    Bugs ship in the wiring, not just the pure layer — the pure fallback can
    be perfect while section 3f still hands it only the plate-solved rows.
    This runs the real ``build_frames`` on an NGC-5548-shaped fixture: many
    on-target frames plus one parked-mount night, and NOTHING plate-solved.
    """

    def _fixture(self):
        rows = []
        oid = 0
        # 9 on-target frames across three nights, no plate solves at all.
        for night, jd0 in (("2023-03-23", 2460027.6), ("2023-03-24", 2460028.6),
                           ("2023-03-27", 2460031.6)):
            for k in range(3):
                oid += 1
                rows.append(frame_row(
                    obs_rowid=oid, path=f"macalester/mrf/{night}/f{oid}.fts",
                    tree="macalester", target_best="NGC 5548",
                    jd=jd0 + 0.01 * k, pltsolvd=None,
                    ra_deg=214.50 + 0.001 * k, dec_deg=25.12,
                    exptime=256.0, filter="6"))
        # The parked-mount night: dec 8 degrees off, RA marching sidereally.
        for k, (ra, dec) in enumerate(
                [(209.43, 33.35), (213.00, 33.34), (222.52, 33.33)]):
            oid += 1
            rows.append(frame_row(
                obs_rowid=oid, path=f"macalester/mrf/2023-03-25/f{oid}.fts",
                tree="macalester", target_best="NGC 5548",
                jd=2460029.6 + 0.02 * k, pltsolvd=None,
                ra_deg=ra, dec_deg=dec, exptime=256.0, filter="6"))
        return pd.DataFrame(rows)

    def _built(self):
        return build.build_frames(self._fixture(),
                                  {"NGC 5548": "ngc5548"},
                                  {"ngc5548": "NGC 5548"})[0]

    def test_unsolved_target_still_gets_offsets(self):
        frames = self._built()
        # Before the fix: 100% NULL, because no frame was plate-solved.
        assert frames["pointing_offset_deg"].notna().all()
        assert (frames["pointing_ref_basis"] == "header_median").all()

    def test_the_parked_mount_night_is_flagged(self):
        frames = self._built()
        bad = frames[frames["path"].str.contains("2023-03-25")]
        good = frames[~frames["path"].str.contains("2023-03-25")]
        assert len(bad) == 3 and len(good) == 9
        assert (bad["pointing_offset_deg"] > m.POINTING_OUTLIER_DEG).all()
        assert bad["qc_flags"].str.contains("pointing_gt1deg").all()
        # …and the on-target frames are NOT swept up with them.
        assert (good["pointing_offset_deg"] < 0.1).all()
        assert not good["qc_flags"].str.contains("pointing_gt1deg").any()

    def test_solved_targets_keep_the_plate_solved_basis(self):
        df = self._fixture()
        df.loc[df.index[:9], "pltsolvd"] = 1
        frames = build.build_frames(df, {"NGC 5548": "ngc5548"},
                                    {"ngc5548": "NGC 5548"})[0]
        assert (frames["pointing_ref_basis"] == "plate_solved").all()
        # The solved subset is the on-target one, so the verdict is the same.
        bad = frames[frames["path"].str.contains("2023-03-25")]
        assert bad["qc_flags"].str.contains("pointing_gt1deg").all()


class TestBuildFramesEraAssignment:
    """Regression: eras whose key contains None (blank EGAIN/READOUTM) must
    still assign era_id to their frames — 29 of 83 eras were orphaned in
    the first shipped manifest."""

    def _fixture(self):
        rows = [
            # Three frames of a missing-EGAIN configuration (the Andor iKon
            # commissioning epoch's failure mode: egain absent from headers).
            frame_row(obs_rowid=1, path="rawimage/2024-04-04/i1.fts",
                      jd=2460405.7, readoutm="1MHz High Sensitivity 16-bit",
                      naxis1=2048.0, naxis2=2048.0, egain=float("nan")),
            frame_row(obs_rowid=2, path="rawimage/2024-04-05/i2.fts",
                      jd=2460406.7, readoutm="1MHz High Sensitivity 16-bit",
                      naxis1=2048.0, naxis2=2048.0, egain=None),
            # A fully-keyed era for contrast.
            frame_row(obs_rowid=3, path="rawimage/2024-01-01/g1.fts",
                      jd=2460310.6),
            # Blank READOUTM as well — still a valid (single) era.
            frame_row(obs_rowid=4, path="rawimage/2024-02-01/p1.fts",
                      jd=2460341.6, readoutm=None, naxis1=8.0,
                      naxis2=3211.0, xbinning=2.0, egain=None),
            # An unreadable-header row: the ONLY case allowed era_id NULL.
            frame_row(obs_rowid=5, path="rawimage/2024-03-01/bad.fts",
                      jd=None, error="OSError: Empty FITS",
                      readoutm=None, naxis1=None, naxis2=None,
                      xbinning=None, egain=None),
        ]
        return pd.DataFrame(rows)

    def test_every_error_free_frame_gets_an_era(self):
        frames, eras = build.build_frames(
            self._fixture(), {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        ok = frames["error"].isna()
        # The build-time contract: only header-error rows may lack era_id.
        assert not (frames["era_id"].isna() & ok).any()
        assert frames.loc[~ok, "era_id"].isna().all()

    def test_eras_table_and_frames_column_agree(self):
        frames, eras = build.build_frames(
            self._fixture(), {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        ok = frames["error"].isna()
        # No orphan eras: every declared era_id is carried by >=1 frame,
        # and vice versa (the shipped bug orphaned 29 of 83 eras).
        declared = set(eras["era_id"])
        carried = set(frames.loc[ok, "era_id"].dropna().astype(int))
        assert declared == carried
        # And the per-era frame counts reconcile exactly.
        for _, era in eras.iterrows():
            n_in_frames = int((frames["era_id"] == era["era_id"]).sum())
            assert n_in_frames == era["n_frames"]


class TestBuildFramesCanonicalSelection:
    """The vectorized dedup ranking must equal choose_canonical()
    member-by-member — single source of truth for the tree policy."""

    def test_vectorized_selection_matches_choose_canonical(self):
        rows = [
            # Group A: cross-tree — rawimage must beat external and reduced.
            frame_row(obs_rowid=1, path="reduced/x/a.fts", tree="reduced",
                      jd=2460100.5),
            frame_row(obs_rowid=2, path="rawimage/2023-05-21/a.fts",
                      tree="rawimage", jd=2460100.5),
            frame_row(obs_rowid=3, path="external/x/a.fts", tree="external",
                      jd=2460100.5),
            # Group B: within-tree tie — earliest path (original night) wins.
            frame_row(obs_rowid=4, path="rawimage/2023-07-07/b.fts",
                      jd=2460101.5),
            frame_row(obs_rowid=5, path="rawimage/2023-05-21/b.fts",
                      jd=2460101.5),
            # Group C: the documented NGC 5548 exception — macalester wins.
            frame_row(obs_rowid=6, path="rawimage/2023-03-24/c.fts",
                      jd=2460102.5, target_best="NGC 5548"),
            frame_row(obs_rowid=7, path="macalester/x/c.fts",
                      tree="macalester", jd=2460102.5,
                      target_best="NGC 5548"),
            # Group D: unknown tree, sole copy — still canonical.
            frame_row(obs_rowid=8, path="mystery/d.fts", tree="mystery",
                      jd=2460103.5),
        ]
        df = pd.DataFrame(rows)
        key_of_raw = {"T CrB": "tcrb", "NGC 5548": "ngc5548"}
        display = {"tcrb": "T CrB", "ngc5548": "NGC 5548"}
        frames, _ = build.build_frames(df, key_of_raw, display)
        # Exactly one canonical row per group, and it is the row the pure,
        # unit-tested choose_canonical() picks.
        for _, grp in frames.groupby("dup_group"):
            assert int(grp["is_canonical"].sum()) == 1
            members = list(zip(grp["tree"], grp["path"]))
            want_label = grp.index[
                m.choose_canonical(members,
                                   target_key=grp["target_key"].iloc[0])]
            got_label = grp.index[grp["is_canonical"] == 1][0]
            assert got_label == want_label


class TestSynonymDisplayDirection:
    """Regression: the SYNONYM_TABLE documents merge DIRECTION, so the
    displayed canonical name must be the destination's native name even
    when the merged-in name has more catalog rows ('Alpha Lyr' at 794 rows
    displayed as the canonical name of the vega group in the first build)."""

    def _fixture(self):
        rows = []
        # 10 'Alpha Lyr' rows — the bigger population, plate-solved at Vega.
        for i in range(10):
            rows.append(frame_row(
                obs_rowid=1 + i, path=f"rawimage/x/al{i}.fts",
                target_best="Alpha Lyr", jd=2460200.5 + i,
                pltsolvd=1, ra_deg=279.23, dec_deg=38.78))
        # 2 native 'Vega ...' ladder rows — fewer, also solved at Vega.
        for i in range(2):
            rows.append(frame_row(
                obs_rowid=100 + i, path=f"rawimage/x/v{i}.fts",
                target_best="Vega 0p001s lrg 5", jd=2460300.5 + i,
                pltsolvd=1, ra_deg=279.24, dec_deg=38.78))
        return pd.DataFrame(rows)

    def test_native_name_wins_the_display_vote(self):
        aliases, key_of_raw, display_of_key = build.resolve_aliases(
            self._fixture())
        # The merge itself: alphalyr -> vega (cone-gated, coords agree).
        assert key_of_raw["Alpha Lyr"] == "vega"
        assert key_of_raw["Vega 0p001s lrg 5"] == "vega"
        # The display name follows the documented arrow: 'Vega', never
        # 'Alpha Lyr', regardless of row counts.
        assert display_of_key["vega"] == "Vega"
        shown = set(aliases.loc[aliases["target_key"] == "vega",
                                "canonical_target"])
        assert shown == {"Vega"}


class TestStrategyClaimsTable:
    """Guard the claims table's shape: a malformed entry (wrong arity, or a
    metric the build script does not implement) must fail HERE, not at
    2 a.m. inside a full catalog build."""

    def test_every_claim_is_well_formed(self):
        implemented = {"rows_all_trees", "unique_light",
                       "grism_light", "grism4_light"}
        for claim in m.STRATEGY_CLAIMS:
            assert len(claim) == 6, f"claim arity drifted: {claim!r}"
            project, tkey, metric, cf, cn, source = claim
            assert metric in implemented, f"unimplemented metric: {metric}"
            assert isinstance(project, str) and isinstance(tkey, str)
            assert cf is None or isinstance(cf, int)
            assert cn is None or isinstance(cn, int)
            assert isinstance(source, str) and source


class TestEraIdRegistry:
    """Regression (2026-08-18): era ids are a registry — a rebuild that sees
    NEW camera configurations must never renumber ids already published.
    Trigger: the Calibrations/ recovery added mid-timeline configurations
    that a pure first-on-sky ordering would have spliced into, silently
    re-pointing every published "era N" reference at a different camera."""

    def _first_build_fixture(self):
        # Two configurations: A (oldest on sky) and B (newer) -> ids 1, 2.
        return pd.DataFrame([
            frame_row(obs_rowid=1, path="rawimage/2023-06-08/a1.fts",
                      jd=2460103.8, readoutm="High Gain"),
            frame_row(obs_rowid=2, path="rawimage/2024-06-08/b1.fts",
                      jd=2460469.8, readoutm="Mode0",
                      naxis1=4788.0, naxis2=3194.0, xbinning=2.0,
                      egain=0.2467),
        ])

    def _second_build_fixture(self):
        # Same two configurations PLUS a new one (C) whose first JD falls
        # BETWEEN A and B — exactly the splice case.
        df = self._first_build_fixture()
        extra = pd.DataFrame([
            frame_row(obs_rowid=3, path="Calibrations/masters/c1.fts",
                      jd=2460300.5, readoutm="Mode0",
                      naxis1=9576.0, naxis2=6388.0, xbinning=1.0,
                      egain=0.2467),
        ])
        return pd.concat([df, extra], ignore_index=True)

    def _registry_of(self, eras):
        # Rebuild {era_key: era_id} exactly the way load_prior_era_ids does.
        return {m.era_key(r.readoutm, r.naxis1, r.naxis2, r.xbinning,
                          r.egain): int(r.era_id)
                for r in eras.itertuples()}

    def test_new_mid_timeline_key_never_renumbers_published_ids(self):
        # First build: no registry -> ids by first-on-sky order.
        _, eras1 = build.build_frames(
            self._first_build_fixture(), {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        reg1 = self._registry_of(eras1)
        # Second build: pass the first build's registry, add the splicer.
        frames2, eras2 = build.build_frames(
            self._second_build_fixture(), {"T CrB": "tcrb"},
            {"tcrb": "T CrB"}, prior_era_ids=reg1)
        reg2 = self._registry_of(eras2)
        # Every previously-published id survives unchanged ...
        for key, eid in reg1.items():
            assert reg2[key] == eid, "published era id was renumbered"
        # ... and the newcomer appends AFTER the existing maximum instead of
        # splicing into the middle, despite its earlier first-on-sky JD.
        new_ids = set(reg2.values()) - set(reg1.values())
        assert new_ids == {max(reg1.values()) + 1}

    def test_without_registry_ordering_is_first_on_sky(self):
        # First-build behaviour is untouched: oldest configuration is era 1.
        frames, eras = build.build_frames(
            self._first_build_fixture(), {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        oldest = frames.loc[frames["jd"].idxmin(), "era_id"]
        assert oldest == 1


class TestSiblingTablePreservation:
    """Regression (2026-08-18): an S0 rebuild must NOT destroy tables added
    by downstream stages.  S0 swaps a freshly built temp file over the live
    manifest, and that whole-file swap wiped s1_strata, s1_solve_experiment,
    s1_failure_autopsy and detector_params during the Calibrations ingest —
    the accepted evidence behind the astrometry verdict and detector memo."""

    def _tiny_s0_frames(self):
        return pd.DataFrame([
            frame_row(obs_rowid=1, path="rawimage/2024-01-01/a.fts",
                      jd=2460310.6),
        ])

    def _write_first_manifest(self, out):
        frames, eras = build.build_frames(
            self._tiny_s0_frames(), {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        aliases = pd.DataFrame([{"target_best": "T CrB", "target_key": "tcrb",
                                 "canonical_target": "T CrB", "n_frames": 1,
                                 "method": "identity",
                                 "cone_check_passed": None}])
        counts = pd.DataFrame([{"project": "TCrB", "target_key": "tcrb",
                                "metric": "unique_light", "claim_frames": 1,
                                "claim_nights": 1, "manifest_frames": 1,
                                "manifest_nights": 1, "source": "test"}])
        build.write_manifest(out, frames, aliases, eras, counts, Path("cat.db"))

    def test_downstream_tables_survive_a_rebuild(self, tmp_path):
        out = tmp_path / "rlmt-manifest.sqlite"
        self._write_first_manifest(out)
        # A downstream stage adds its own evidence table, with an index.
        with closing(sqlite3.connect(out)) as con:
            con.execute("CREATE TABLE detector_params "
                        "(era_group TEXT, quantity TEXT, value REAL)")
            con.execute("INSERT INTO detector_params VALUES "
                        "('High Gain', 'ceiling_adu', 3496.0)")
            con.execute("CREATE INDEX ix_dp ON detector_params(era_group)")
            con.commit()
        # S0 rebuilds over the top of it.
        self._write_first_manifest(out)
        with closing(sqlite3.connect(out)) as con:
            rows = con.execute("SELECT era_group, quantity, value "
                               "FROM detector_params").fetchall()
            names = {r[0] for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='index'")}
            carried = con.execute(
                "SELECT value FROM build_meta WHERE key='carried_tables'"
            ).fetchone()[0]
        assert rows == [("High Gain", "ceiling_adu", 3496.0)], \
            "downstream evidence table was destroyed by the S0 rebuild"
        assert "ix_dp" in names, "carried table lost its index"
        assert "detector_params" in carried, \
            "build_meta must record what was carried (staleness is auditable)"

    def test_s0_owned_tables_are_rebuilt_not_carried(self, tmp_path):
        out = tmp_path / "rlmt-manifest.sqlite"
        self._write_first_manifest(out)
        # Poison an S0-owned table: the rebuild must overwrite it, never
        # carry it forward (otherwise stale frames rows would accumulate).
        with closing(sqlite3.connect(out)) as con:
            con.execute("INSERT INTO frames (path) VALUES ('POISON')")
            con.commit()
        self._write_first_manifest(out)
        with closing(sqlite3.connect(out)) as con:
            n_poison = con.execute(
                "SELECT COUNT(*) FROM frames WHERE path='POISON'").fetchone()[0]
        assert n_poison == 0


class TestStaleWalRemovedOnSwap:
    """Regression (2026-08-19): the atomic swap left the previous database's
    write-ahead log beside the new file, and SQLite replayed that log against
    the new database — reporting 'malformed database schema' on a manifest
    that was internally perfect.  Rebuilding must clear the sidecars.

    Since 2026-10-03 the default rebuild is an in-place transaction, which
    cannot leave a stale log at all (the log belongs to the same file); the
    scenario is kept as a test because it is the one that actually happened,
    and because the ``--fresh-file`` path still replaces the file."""

    def _write(self, out):
        frames, eras = build.build_frames(
            pd.DataFrame([frame_row(obs_rowid=1,
                                    path="rawimage/2024-01-01/a.fts",
                                    jd=2460310.6)]),
            {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        aliases = pd.DataFrame([{"target_best": "T CrB", "target_key": "tcrb",
                                 "canonical_target": "T CrB", "n_frames": 1,
                                 "method": "identity",
                                 "cone_check_passed": None}])
        counts = pd.DataFrame([{"project": "T", "target_key": "tcrb",
                                "metric": "unique_light", "claim_frames": 1,
                                "claim_nights": 1, "manifest_frames": 1,
                                "manifest_nights": 1, "source": "test"}])
        build.write_manifest(out, frames, aliases, eras, counts, Path("cat.db"))

    def test_rebuild_clears_stale_wal_and_database_stays_readable(self, tmp_path):
        out = tmp_path / "rlmt-manifest.sqlite"
        self._write(out)
        # Simulate a sibling stage that had the manifest open in WAL mode and
        # was STILL HOLDING IT when the rebuild ran — which is exactly what
        # happened (several workflows shared the manifest).  The connection is
        # deliberately left open: closing it cleanly checkpoints the WAL away,
        # and it is the un-checkpointed sidecar that does the damage.
        sibling = sqlite3.connect(out)
        sibling.execute("PRAGMA journal_mode=WAL")
        sibling.execute("CREATE TABLE sibling (x)")
        sibling.execute("INSERT INTO sibling VALUES (1)")
        sibling.commit()
        assert Path(str(out) + "-wal").exists(), "fixture needs a live WAL"

        try:
            self._write(out)                  # the rebuild that used to corrupt
        finally:
            sibling.close()

        assert not Path(str(out) + "-wal").exists(), "stale WAL survived the swap"
        assert not Path(str(out) + "-shm").exists(), "stale SHM survived the swap"
        with closing(sqlite3.connect(out)) as con:
            assert con.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
            # And the rebuild's own content is readable, which is the symptom
            # that was lost: the old failure raised on the first real query.
            assert con.execute("SELECT count(*) FROM frames").fetchone()[0] == 1


# ---------------------------------------------------------------------------
# F-1: one exposure, many filenames (plan review 2026-10-03, TE.F6)
# ---------------------------------------------------------------------------
class TestExposureName:
    """Every name here is a pattern taken from the real catalog."""

    def test_a_raw_name_is_its_own_base(self):
        n = m.exposure_name("knh_Berkeley_39_g_120s_2026-03-22T03-50-31.fts.fz")
        assert n.stem == "knh_Berkeley_39_g_120s_2026-03-22T03-50-31"
        assert n.base == n.stem
        assert n.base_no_counter is None
        assert n.n_processing == 0

    def test_calibrated_suffix_points_at_the_raw_name(self):
        n = m.exposure_name(
            "knh_Berkeley_39_g_120s_2026-03-22T03-50-31_calibrated.fts.fz")
        assert n.base == "knh_Berkeley_39_g_120s_2026-03-22T03-50-31"
        assert n.n_processing == 1

    def test_stacked_suffixes_count_their_depth(self):
        n = m.exposure_name("mel_M11_g_10s_2026-06-29T09-26-37_wcs_calibrated.fts.fz")
        assert n.base == "mel_M11_g_10s_2026-06-29T09-26-37"
        assert n.n_processing == 2

    def test_packaging_does_not_change_the_stem(self):
        assert m.exposure_name("master_flat_g.fts").stem == \
            m.exposure_name("master_flat_g.fts.fz").stem == "master_flat_g"
        assert m.frame_stem("x.fits.fz") == m.frame_stem("x.fit") == "x"

    def test_counter_before_the_suffix_is_offered_not_assumed(self):
        n = m.exposure_name(
            "mpg_M87_g_90s_2026-06-29T05-16-00_1_wcs_calibrated.fts.fz")
        assert n.base == "mpg_M87_g_90s_2026-06-29T05-16-00_1"
        assert n.base_no_counter == "mpg_M87_g_90s_2026-06-29T05-16-00"

    def test_trailing_text_after_the_suffix_is_dropped_with_it(self):
        # 'BeStars_…_r_0_wcs_2024-05-10T08-06-32' is a solved copy of
        # 'BeStars_…_r_0'; 'ias_…_cal_Nov19' a class's calibrated copy.
        assert m.exposure_name(
            "BeStars_HD138749_0p08s_r_0_wcs_2024-05-10T08-06-32.fts.fz"
        ).base == "BeStars_HD138749_0p08s_r_0"
        assert m.exposure_name(
            "ias_M16_blue_60s_2024-11-15T01-11-05_cal_Nov19.fts.fz"
        ).base == "ias_M16_blue_60s_2024-11-15T01-11-05"

    def test_a_leading_token_is_never_a_suffix(self):
        # An observer code can be any three letters, 'cal' included.
        n = m.exposure_name("cal_M13_g_10s_2024-05-01T03-00-00.fts.fz")
        assert n.n_processing == 0

    def test_red_is_a_filter_not_a_reduction_suffix(self):
        n = m.exposure_name("knh_knh001_ngc3169_red_60s_2024-04-30T03-16-15.fts")
        assert n.n_processing == 0

    def test_longer_words_containing_a_token_are_untouched(self):
        n = m.exposure_name("grismCalibration_BD-00_3227_180s_hrg_0.fts.fz")
        assert n.n_processing == 0


class TestExposureRoot:
    RAW = "x_g_10s_2026-06-29T05-16-00"

    def test_a_raw_frame_is_its_own_root(self):
        name = m.exposure_name(self.RAW + ".fts.fz")
        assert m.exposure_root(name, lambda s: False) == \
            (self.RAW, m.ROOT_IS_SELF)

    def test_a_twin_roots_at_its_raw_parent(self):
        name = m.exposure_name(self.RAW + "_calibrated.fts.fz")
        assert m.exposure_root(name, lambda s: s == self.RAW) == \
            (self.RAW, m.ROOT_BY_SUFFIX)

    def test_the_counter_is_dropped_only_when_that_parent_exists(self):
        name = m.exposure_name(self.RAW + "_1_wcs.fts.fz")
        # The parent WITH the counter exists: it is a frame index.
        assert m.exposure_root(name, lambda s: s == self.RAW + "_1") == \
            (self.RAW + "_1", m.ROOT_BY_SUFFIX)
        # Only the parent WITHOUT it exists: it was a copy counter.
        assert m.exposure_root(name, lambda s: s == self.RAW) == \
            (self.RAW, m.ROOT_BY_COUNTER)

    def test_an_orphan_derivative_roots_at_its_own_base(self):
        # 2,054 '_wcs' frames of 2024-04/05 are the ONLY copy: pyscope
        # wrote the solved file and nothing else.  They stay canonical.
        name = m.exposure_name(self.RAW + "_wcs.fts.fz")
        assert m.exposure_root(name, lambda s: False) == \
            (self.RAW, m.ROOT_ORPHAN)

    def test_two_derivatives_of_one_missing_parent_share_a_root(self):
        a = m.exposure_root(m.exposure_name(self.RAW + "_wcs.fts.fz"),
                            lambda s: False)[0]
        b = m.exposure_root(m.exposure_name(self.RAW + "_wcs_cal.fts.fz"),
                            lambda s: False)[0]
        assert a == b == self.RAW


class TestDupBasis:
    def test_every_kind_of_duplicate_is_named(self):
        raw = "a_g_10s_T.fts.fz"
        assert m.dup_basis(raw, raw, True) == m.DUP_CANONICAL
        assert m.dup_basis(raw, raw, False) == m.DUP_SAME_BASENAME
        assert m.dup_basis("a_g_10s_T.fts", raw, False) == m.DUP_PACKAGING
        assert m.dup_basis("a_g_10s_T_calibrated.fts.fz", raw, False,
                           m.ROOT_BY_SUFFIX) == m.DUP_PROCESSING
        assert m.dup_basis("a_g_10s_T_1_wcs.fts.fz", raw, False,
                           m.ROOT_BY_COUNTER) == m.DUP_PROCESSING_COUNTER


class TestCanonicalPrefersTheLeastProcessed:
    def test_raw_beats_its_calibrated_twin_in_the_same_tree(self):
        # 138 '_calibrated' files live inside rawimage itself.
        members = [("rawimage", "rawimage/n/a_calibrated.fts.fz"),
                   ("rawimage", "rawimage/n/a.fts.fz")]
        assert m.choose_canonical(members) == 1

    def test_processing_depth_outranks_tree_priority(self):
        # An unprocessed copy in the last-resort tree still beats a
        # processed copy in the best one: untouched pixels come first.
        members = [("rawimage", "rawimage/n/a_calibrated.fts.fz"),
                   ("reduced", "reduced/n/a.fts")]
        assert m.choose_canonical(members) == 1

    def test_among_derivatives_the_shallowest_wins(self):
        members = [("rawimage", "rawimage/n/a_wcs_cal.fts.fz"),
                   ("rawimage", "rawimage/n/a_wcs.fts.fz")]
        assert m.choose_canonical(members) == 1


class TestBuildFramesTwinDedup:
    """F-1 at the wiring: the vectorized build must group twins, pick the
    raw member, and label every duplicate with its evidence."""

    T = "x_g_120s_2026-03-22T03-50-31"

    def _fixture(self):
        t, qhy = self.T, dict(readoutm="Fast", xbinning=2.0, egain=1.0)
        rows = [
            # The raw QHY frame (4800x3211, era 'raw') ...
            frame_row(obs_rowid=1, path=f"rawimage/2026-03-22/{t}.fts.fz",
                      jd=2461121.66, naxis1=4800.0, naxis2=3211.0, **qhy),
            # ... its reduced twin: overscan trimmed, so a DIFFERENT era.
            frame_row(obs_rowid=2, tree="reduced",
                      path=f"reduced/2026-03-22/{t}_calibrated.fts.fz",
                      jd=2461121.66, naxis1=4787.0, naxis2=3193.0, **qhy),
            # ... and an exact copy of the raw frame in another tree.
            frame_row(obs_rowid=3, tree="macalester",
                      path=f"macalester/x/{t}.fts.fz",
                      jd=2461121.66, naxis1=4800.0, naxis2=3211.0, **qhy),
            # A second exposure with the SAME name pattern, different JD:
            # its twin must not attach to the first frame.
            frame_row(obs_rowid=4, path="rawimage/2026-03-22/y.fts.fz",
                      jd=2461121.70, naxis1=4800.0, naxis2=3211.0, **qhy),
            frame_row(obs_rowid=5, tree="reduced",
                      path="reduced/2026-03-22/y_calibrated.fts.fz",
                      jd=2461121.70, naxis1=4787.0, naxis2=3193.0, **qhy),
            # A twin whose JD was REWRITTEN in reduction: not provably the
            # same exposure, stays its own group (S0b's ladder links it).
            frame_row(obs_rowid=6, tree="reduced",
                      path="reduced/2026-03-22/y_calibrated_drift.fts",
                      jd=2461121.7004, naxis1=4787.0, naxis2=3193.0, **qhy),
            # A solved-only frame with no raw parent: canonical on its own.
            frame_row(obs_rowid=7, path="rawimage/2024-04-29/z_wcs.fts.fz",
                      jd=2460429.68),
            # Packaging twin: same stem, compressed and not.
            frame_row(obs_rowid=8, tree="calib", path="calib/m/mf.fts.fz",
                      jd=2460600.5),
            frame_row(obs_rowid=9, tree="Calibrations",
                      path="Calibrations/m/mf.fts", jd=2460600.5),
        ]
        return pd.DataFrame(rows)

    def _build(self):
        return build.build_frames(self._fixture(), {"T CrB": "tcrb"},
                                  {"tcrb": "T CrB"})

    def test_twins_join_their_raw_parent_and_the_raw_frame_wins(self):
        frames, _ = self._build()
        f = frames.set_index("obs_rowid")
        assert f.loc[1, "dup_group"] == f.loc[2, "dup_group"] \
            == f.loc[3, "dup_group"]
        assert list(f.loc[[1, 2, 3], "is_canonical"]) == [1, 0, 0]
        assert f.loc[2, "dup_basis"] == m.DUP_PROCESSING
        assert f.loc[3, "dup_basis"] == m.DUP_SAME_BASENAME
        assert f.loc[1, "dup_basis"] == m.DUP_CANONICAL

    def test_a_different_jd_never_merges(self):
        f = self._build()[0].set_index("obs_rowid")
        assert f.loc[4, "dup_group"] == f.loc[5, "dup_group"]
        assert f.loc[4, "dup_group"] != f.loc[1, "dup_group"]
        # The JD-drifted twin is NOT merged: identity needs the same JD.
        assert f.loc[6, "dup_group"] != f.loc[4, "dup_group"]
        assert f.loc[6, "is_canonical"] == 1

    def test_a_derivative_with_no_parent_stays_canonical(self):
        f = self._build()[0].set_index("obs_rowid")
        assert f.loc[7, "is_canonical"] == 1
        assert f.loc[7, "dup_basis"] == m.DUP_CANONICAL

    def test_packaging_twins_merge(self):
        f = self._build()[0].set_index("obs_rowid")
        assert f.loc[8, "dup_group"] == f.loc[9, "dup_group"]
        assert f.loc[8, "is_canonical"] == 1          # 'calib' is listed
        assert f.loc[9, "dup_basis"] == m.DUP_PACKAGING

    def test_exactly_one_canonical_row_per_group(self):
        frames, _ = self._build()
        assert (frames.groupby("dup_group")["is_canonical"].sum() == 1).all()

    def test_vectorized_choice_still_equals_the_pure_function(self):
        frames, _ = self._build()
        for _, grp in frames.groupby("dup_group"):
            members = list(zip(grp["tree"], grp["path"]))
            want = grp.index[m.choose_canonical(
                members, target_key=grp["target_key"].iloc[0])]
            assert grp.index[grp["is_canonical"] == 1][0] == want

    def test_the_reduced_era_becomes_an_alias_of_the_raw_era(self):
        """Acceptance line of F-1: 'eras 79/82 become reduced aliases'.

        The trimmed-geometry era holds rows 2, 5 and 6.  Rows 2 and 5 are
        copies of raw-era frames; row 6 (JD drift) is not — so with the
        test's 2-of-3 the era is NOT an alias at the 95% threshold, and
        becomes one when the drifted row is removed."""
        _, eras = self._build()
        trimmed = eras[eras["naxis1"] == 4787].iloc[0]
        assert trimmed["n_rows"] == 3 and trimmed["n_canonical"] == 1
        assert pd.isna(trimmed["alias_of_era"])

        df = self._fixture()
        frames, eras = build.build_frames(df[df["obs_rowid"] != 6],
                                          {"T CrB": "tcrb"},
                                          {"tcrb": "T CrB"})
        trimmed = eras[eras["naxis1"] == 4787].iloc[0]
        raw = eras[eras["naxis1"] == 4800].iloc[0]
        assert trimmed["n_canonical"] == 0
        assert trimmed["alias_of_era"] == raw["era_id"]
        assert trimmed["alias_fraction"] == 1.0
        # An era of real exposures is nobody's alias.
        assert pd.isna(raw["alias_of_era"])
        # And the alias keeps its registry id: pinned, never reused.
        assert trimmed["era_id"] != raw["era_id"]


# ---------------------------------------------------------------------------
# F-2: typed hardware-state columns from the header re-scrape
# ---------------------------------------------------------------------------
class TestHardwareColumns:
    def test_camera_identity_from_the_driver_string(self):
        assert m.camera_id("DL Imaging", None) == "AC4040"
        assert m.camera_id("Andor CCD/EMCCD (SDK2)", "ASCOM") == "iKon"
        assert m.camera_id("ASI Camera (1)", "ASCOM") == "ASI"
        # Three driver strings, one camera.
        assert m.camera_id("QHYCCD-Cameras-Capture", "ASCOM") == "QHY600"
        assert m.camera_id("QHY600Pro", None) == "QHY600"

    def test_camname_identifies_the_frames_with_no_instrume(self):
        # 854 pyscope-native frames of 2026-06-28/29 have no INSTRUME card.
        assert m.camera_id(None, "QHY600MPCIE-ec89ab488ac5aca73") == "QHY600"

    def test_an_unknown_camera_is_none_not_a_guess(self):
        assert m.camera_id("Some New Driver", "ASCOM") is None
        assert m.camera_id(None, None) is None

    def test_card_float_keeps_zero_and_refuses_non_numbers(self):
        assert m.card_float("-0.0") == 0.0          # a real set-point
        assert m.card_float("  -10.000000000000000 ") == -10.0
        assert m.card_float("4x") is None           # the iKon's preamp GAIN
        assert m.card_float("") is None and m.card_float(None) is None

    def test_wheel_map_from_both_header_dialects(self):
        maxim = "('g', 'OGGrism', 'r')"
        pyscope = "['g', 'lrg', 'r']"
        assert m.wheel_slots(maxim) == ("g", "OGGrism", "r")
        assert m.wheel_map_key(pyscope) == "g|lrg|r"
        # Case is preserved: a rename is evidence, not noise to fold away.
        assert m.wheel_map_key("('Lum',)") == "Lum"
        assert m.wheel_map_key(None) is None and m.wheel_map_key("") is None

    def test_typed_row_keeps_blank_apart_from_absent(self):
        row = m.hardware_columns({
            "h_instrume": "ASI Camera (1)", "h_camname": "ASCOM",
            "h_gain": "100", "h_offset": "30", "h_set_temp": "-10.0",
            "h_coolpowr": "28", "h_focpos": "9980.31", "h_flipstat": "",
            "h_fwpos": "0", "h_fwallnam": "('g', 'r')", "h_scanned": 1})
        assert row["camera"] == "ASI"
        assert row["hdr_gain"] == "100" and row["hdr_gain_num"] == 100.0
        assert row["hdr_offset"] == 30.0 and row["set_temp"] == -10.0
        assert row["flipstat"] == ""        # card present, blank
        assert row["telpier"] is None       # card absent
        assert row["fwpos"] == 0 and row["fwallnam"] == "g|r"
        assert row["hdr_scanned"] == 1

    def test_an_unscanned_row_is_all_null_and_says_so(self):
        row = m.hardware_columns({})
        assert row["hdr_scanned"] == 0
        assert all(v is None for k, v in row.items() if k != "hdr_scanned")
        assert set(row) == {c for c, _ in m.HARDWARE_FRAME_COLUMNS}


class TestAttachHardwareColumns:
    def _frames(self, **hdr):
        rows = [
            frame_row(obs_rowid=1, path="rawimage/n/a.fts.fz", jd=2460310.6,
                      instrume=None, swcreate=None, ccd_temp=None, **hdr),
            # an exact copy that was NOT re-scraped
            frame_row(obs_rowid=2, tree="macalester",
                      path="macalester/x/a.fts.fz", jd=2460310.6,
                      instrume=None, swcreate=None, ccd_temp=None),
            # a reduced derivative that was NOT re-scraped
            frame_row(obs_rowid=3, tree="reduced",
                      path="reduced/n/a_calibrated.fts.fz", jd=2460310.6,
                      instrume=None, swcreate=None, ccd_temp=None),
        ]
        return pd.DataFrame(rows)

    HDR = dict(h_instrume="ASI Camera (1)", h_swcreate="MaxIm DL 6.40",
               h_ccd_temp="-9.9", h_set_temp="-10", h_gain="100",
               h_flipstat="Flip/Mirror", h_scanned=1)

    def test_nulls_are_filled_where_the_card_exists(self):
        frames, _ = build.build_frames(self._frames(**self.HDR),
                                       {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        f = frames.set_index("obs_rowid")
        assert f.loc[1, "instrume"] == "ASI Camera (1)"
        assert f.loc[1, "swcreate"] == "MaxIm DL 6.40"
        assert f.loc[1, "ccd_temp"] == -9.9
        assert f.loc[1, "camera"] == "ASI" and f.loc[1, "set_temp"] == -10.0
        assert f.loc[1, "hdr_scanned"] == build.HDR_SCANNED
        # The raw h_* columns do not leak into frames.
        assert not [c for c in frames.columns if c.startswith("h_")]

    def test_an_exact_copy_inherits_and_a_derivative_does_not(self):
        frames, _ = build.build_frames(self._frames(**self.HDR),
                                       {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        f = frames.set_index("obs_rowid")
        # Same bytes, same header: inherited, and labelled as inherited.
        assert f.loc[2, "hdr_scanned"] == build.HDR_INHERITED
        assert f.loc[2, "camera"] == "ASI" and f.loc[2, "set_temp"] == -10.0
        # A calibrated derivative has its own header: never inherited.
        assert f.loc[3, "hdr_scanned"] == build.HDR_NOT_SCANNED
        assert pd.isna(f.loc[3, "set_temp"])

    def test_a_catalog_value_is_never_overwritten(self):
        df = self._frames(**self.HDR)
        df.loc[0, "instrume"] = "DL Imaging"        # the original scan's
        frames, _ = build.build_frames(df, {"T CrB": "tcrb"},
                                       {"tcrb": "T CrB"})
        assert frames.set_index("obs_rowid").loc[1, "instrume"] == \
            "DL Imaging"

    def test_without_a_rescrape_the_columns_exist_and_are_null(self):
        frames, _ = build.build_frames(self._frames(), {"T CrB": "tcrb"},
                                       {"tcrb": "T CrB"})
        assert set(frames["hdr_scanned"]) == {build.HDR_NOT_SCANNED}
        assert frames["camera"].isna().all()


# ---------------------------------------------------------------------------
# The in-place swap (2026-10-03): siblings untouched, writers not lost
# ---------------------------------------------------------------------------
class TestInPlaceSwap:
    def _tables(self, n=1):
        frames, eras = build.build_frames(
            pd.DataFrame([frame_row(obs_rowid=i,
                                    path=f"rawimage/2024-01-01/a{i}.fts",
                                    jd=2460310.6 + i) for i in range(n)]),
            {"T CrB": "tcrb"}, {"tcrb": "T CrB"})
        aliases = pd.DataFrame([{"raw_name": "T CrB", "target_key": "tcrb"}])
        counts = pd.DataFrame([{"project": "T", "target_key": "tcrb"}])
        return frames, aliases, eras, counts

    def _write(self, out, n=1, **kw):
        frames, aliases, eras, counts = self._tables(n)
        build.write_manifest(out, frames, aliases, eras, counts,
                             Path("cat.db"), **kw)

    def _meta(self, out):
        with closing(sqlite3.connect(out)) as con:
            return dict(con.execute("SELECT key, value FROM build_meta"))

    def test_first_build_writes_a_file_and_rebuild_swaps_in_place(
            self, tmp_path):
        out = tmp_path / "m.sqlite"
        self._write(out)
        assert self._meta(out)["swap"] == "fresh_file_replace"
        inode = out.stat().st_ino
        self._write(out, n=3)
        assert self._meta(out)["swap"] == "in_place_transaction"
        # THE POINT: the file is the same file.  A process that had it open
        # is still attached to the live manifest, not to an orphaned inode.
        assert out.stat().st_ino == inode
        with closing(sqlite3.connect(out)) as con:
            assert con.execute("SELECT count(*) FROM frames").fetchone()[0] == 3
            idx = {r[0] for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='index'")}
            left = [r[0] for r in con.execute(
                "SELECT name FROM sqlite_master WHERE name LIKE '%_s0_tmp'")]
        assert {"ix_frames_dup", "ix_frames_era", "ix_frames_rowid"} <= idx
        assert left == [], "temporary tables survived the swap"

    def test_a_sibling_writer_holding_the_file_open_loses_nothing(
            self, tmp_path):
        """The defect the file replace still had: a stage that kept the
        manifest open across an S0 rebuild went on writing to an inode
        with no name.  With the transaction swap its later write lands in
        the live manifest."""
        out = tmp_path / "m.sqlite"
        self._write(out)
        sibling = sqlite3.connect(out)
        sibling.execute("PRAGMA journal_mode=WAL")
        sibling.execute("CREATE TABLE detector_params (k TEXT, v REAL)")
        sibling.execute("INSERT INTO detector_params VALUES ('gain', 1.0)")
        sibling.commit()
        try:
            self._write(out, n=2)                 # S0 rebuilds underneath
            # The sibling writes AFTER the rebuild, on its old connection.
            sibling.execute("INSERT INTO detector_params VALUES ('rn', 3.5)")
            sibling.commit()
        finally:
            sibling.close()
        with closing(sqlite3.connect(out)) as con:
            assert con.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
            got = dict(con.execute("SELECT k, v FROM detector_params"))
            n = con.execute("SELECT count(*) FROM frames").fetchone()[0]
        assert got == {"gain": 1.0, "rn": 3.5}
        assert n == 2

    def test_a_reader_mid_rebuild_sees_a_whole_manifest(self, tmp_path):
        out = tmp_path / "m.sqlite"
        self._write(out)
        with closing(sqlite3.connect(out)) as con:
            con.execute("PRAGMA journal_mode=WAL")
        reader = sqlite3.connect(f"file:{out}?mode=ro", uri=True)
        try:
            reader.execute("BEGIN")
            before = reader.execute("SELECT count(*) FROM frames").fetchone()[0]
            self._write(out, n=5)
            # Same snapshot: still the old manifest, complete.
            assert reader.execute(
                "SELECT count(*) FROM frames").fetchone()[0] == before == 1
            reader.execute("COMMIT")
            assert reader.execute(
                "SELECT count(*) FROM frames").fetchone()[0] == 5
        finally:
            reader.close()

    def test_fresh_file_still_carries_siblings(self, tmp_path):
        out = tmp_path / "m.sqlite"
        self._write(out)
        with closing(sqlite3.connect(out)) as con:
            con.execute("CREATE TABLE s1_batch (obs_rowid INTEGER)")
            con.execute("INSERT INTO s1_batch VALUES (7)")
            con.commit()
        self._write(out, n=2, fresh_file=True)
        meta = self._meta(out)
        assert meta["swap"] == "fresh_file_replace"
        assert "s1_batch" in meta["carried_tables"]
        with closing(sqlite3.connect(out)) as con:
            assert con.execute("SELECT * FROM s1_batch").fetchall() == [(7,)]
