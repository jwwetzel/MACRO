"""Unit tests for macro_core.manifest — every pure S0 function.

Each test class mirrors one decision the manifest encodes; the must-NOT
cases (distinct targets staying separate, non-duplicates staying apart) are
as important as the must cases, because a silent false merge poisons every
downstream count.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import math
import sys
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

    # --- the merges that must NOT happen ---------------------------------
    def test_distinct_dwarf_fields_stay_separate(self):
        # The spec's named regression: two different survey fields.
        assert self.key("Dw1403+49") != self.key("Dw1409+51")

    def test_sn_and_calibration_names_stay_themselves(self):
        assert self.key("2023ixf") == "2023ixf"
        assert self.key("2023ixf") != self.key("M101")

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
