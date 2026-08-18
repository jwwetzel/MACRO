"""Unit tests for the CV production layer: series rules, WCS chaining, calib.

The photometry core is already covered by ``test_phot.py``.  What is new in
the production campaign is a set of RULES — which pixels a series may read,
where its saturation ceiling sits, whether its frames can be photometered
at all, and how a plate solution is allowed to replace a triangle match.
Each of them was written because getting it wrong produces a plausible
number rather than an error, so each of them is tested against the case
that would have slipped through:

* a global saturation threshold that vetoes nothing on High Gain,
* an 8-pixel readout strip accepted as a normal frame,
* a series quietly built from half reduced and half raw pixels,
* a WCS chain believed on the strength of six accidental matches,
* a master dark scaled onto the wrong exposure time.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import calib as cal
from macro_phot import photometry as ph
from macro_phot import register as reg
from macro_phot import series as sr


# ---------------------------------------------------------------------------
# Series identity
# ---------------------------------------------------------------------------
class TestSeriesKey:
    def test_round_trip(self):
        k = sr.series_key("stlmi", 76, "g")
        assert k == "stlmi|e76|g"
        assert sr.parse_series_key(k) == ("stlmi", 76, "g")

    def test_case_is_preserved(self):
        # Era 7 writes 'G' and era 76 writes 'g'.  They are different
        # filters on different cameras; folding them together would merge
        # two zero-point systems into one light curve.
        assert sr.series_key("stlmi", 7, "G") != sr.series_key("stlmi", 7, "g")

    def test_missing_filter_gets_a_visible_name(self):
        # NULL loses rows to every SQL '=' comparison; 'none' does not.
        assert sr.series_key("yzcnc", 6, None).endswith("|none")
        assert sr.series_key("yzcnc", 6, "").endswith("|none")

    def test_float_era_id_is_accepted(self):
        # The staging table stores era_id as REAL; the key must still be
        # an integer label so it matches keys built from the eras table.
        assert sr.series_key("vvpup", 72.0, "r") == "vvpup|e72|r"


# ---------------------------------------------------------------------------
# S2 saturation vetoes, per readout mode
# ---------------------------------------------------------------------------
class TestSaturationVeto:
    def test_high_gain_veto_is_twelve_bit(self):
        # THE bug this replaces: the prototype's single 55,000 ADU
        # threshold is 16x above High Gain's entire dynamic range, so it
        # vetoed nothing at all on the 2024 season.
        assert sr.veto_adu("High Gain") == 3200
        assert sr.veto_adu("High Gain") < 55000

    def test_modes_differ_by_orders_of_magnitude(self):
        assert sr.veto_adu("Mode0") == 60200
        assert sr.veto_adu("1MHz High Sensitivity 16-bit") == 59500
        assert sr.veto_adu("High Gain StackPro") == 51500

    def test_stackpro_ceiling_is_about_sixteen_single_reads(self):
        # StackPro frames are SUMS of 16 sub-reads, so their ceiling should
        # land near 16x the single-read High Gain clip.  This is a
        # consistency check on the two independently measured numbers, not
        # a definition: 16 * 3496 = 55,936 vs the measured 56,062.
        implied = 16 * sr.S2_MODE_CEILING_ADU["High Gain"]
        assert abs(implied - sr.S2_MODE_CEILING_ADU["High Gain StackPro"]) < 500

    def test_unmeasured_mode_answers_none_not_a_guess(self):
        # S2 sampled Low Gain and never saw it saturate.  'Unknown' must
        # not silently become 65,535.
        assert sr.veto_adu("Low Gain") is None
        assert sr.veto_adu("5MHz High Sensitivity 16-bit") is None
        assert sr.veto_adu(None) is None

    def test_veto_maps_through_a_measured_reduction(self):
        # Era 76: F = 1.1086, D = 302.8 ADU, pedestal = 1000 (S2 recon).
        v = sr.veto_in_reduced_adu(60200, 1.1086, 302.8, 1000.0)
        assert 54000 < v < 56000
        # Identity reduction leaves the threshold where it was.
        assert sr.veto_in_reduced_adu(60200, 1.0, 0.0, 0.0) == 60200

    def test_veto_map_is_none_safe(self):
        assert sr.veto_in_reduced_adu(None, 1.1, 300.0, 1000.0) is None
        assert sr.veto_in_reduced_adu(60200, 0.0, 300.0, 1000.0) is None

    def test_saturated_mask_adds_the_background_back(self):
        # sep reports peaks on the BACKGROUND-SUBTRACTED image; the S2 veto
        # is a level on the image as digitized.  A star 3,000 ADU above a
        # 300 ADU sky sits at 3,300 and IS saturated on High Gain.
        peaks = np.array([2000.0, 3000.0, 100.0])
        m = sr.saturated_mask(peaks, 300.0, sr.veto_adu("High Gain"))
        assert list(m) == [False, True, False]

    def test_no_veto_flags_nothing(self):
        m = sr.saturated_mask(np.array([1e9]), 0.0, None)
        assert m.tolist() == [False]


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------
class TestProvenance:
    def test_full_reduced_coverage(self):
        prov, why = sr.choose_provenance(1279, 1279, has_master_calib=True)
        assert prov == "server_reduced"
        assert "1279/1279" in why

    def test_partial_coverage_drops_the_minority_it_does_not_mix(self):
        # VV Pup era 76: 721 of 795 staged frames are reduced.  The 74
        # without a counterpart must be EXCLUDED, never filled in with raw
        # pixels — that would put two reduction histories in one series.
        prov, why = sr.choose_provenance(795, 721, has_master_calib=True)
        assert prov == "server_reduced"
        assert "74 unlinked" in why

    def test_no_reduced_tree_uses_local_masters(self):
        # ST LMi era 7: 1,246 staged, none reduced, masters staged.
        prov, why = sr.choose_provenance(1246, 0, has_master_calib=True)
        assert prov == "local_master"

    def test_no_reduced_tree_and_no_masters_is_raw_and_says_so(self):
        prov, why = sr.choose_provenance(100, 0, has_master_calib=False)
        assert prov == "raw"
        assert "flagged" in why

    def test_empty_series(self):
        prov, _ = sr.choose_provenance(0, 0, has_master_calib=True)
        assert prov == "none"


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
class TestPlateScaleResolution:
    def test_header_optics_wins_when_present(self):
        # Era 76: XPIXSZ 7.52 um (binned), FOCALLEN 3454 mm -> 0.449 "/px.
        s, basis = sr.resolve_plate_scale(
            {"xpixsz": 7.52, "focallen": 3454.0, "cd1_1": -0.000125,
             "cd2_2": 0.000125})
        assert basis == "header_optics"
        assert s == pytest.approx(0.449, abs=0.002)

    def test_cd_matrix_rescues_a_zero_focal_length(self):
        # THE bug this fixes: the March-2026 'Fast' frames write
        # FOCALLEN = 0.0, which is a MISSING card, not a short lens.  All
        # 74 frames of one EU UMa block failed extraction outright while
        # carrying a perfectly good CD matrix.
        s, basis = sr.resolve_plate_scale(
            {"xpixsz": 7.52, "focallen": 0.0,
             "cd1_1": -0.0001252494762232, "cd1_2": 0.0,
             "cd2_1": 0.0, "cd2_2": 0.0001252503528998})
        assert basis == "header_cd"
        assert s == pytest.approx(0.4509, abs=0.001)

    def test_rotation_does_not_shrink_the_scale(self):
        # A 45-degree rotated CD matrix has |CD1_1| = scale/sqrt(2); using
        # that alone would understate the scale by 30%.  The determinant
        # is rotation-invariant, which is why it is what we take.
        import math as _m
        d, r = 0.000125, _m.radians(45.0)
        s = sr.plate_scale_from_cd(-d * _m.cos(r), d * _m.sin(r),
                                   d * _m.sin(r), d * _m.cos(r))
        assert s == pytest.approx(0.45, abs=0.005)

    def test_cdelt_is_the_last_resort(self):
        s, basis = sr.resolve_plate_scale(
            {"xpixsz": None, "focallen": None,
             "cdelt1": -0.000125, "cdelt2": 0.000125})
        assert basis == "header_cdelt"
        assert s == pytest.approx(0.45, abs=0.005)

    def test_no_route_names_the_failure(self):
        s, basis = sr.resolve_plate_scale({"xpixsz": None, "focallen": None})
        assert s is None and basis == "none"

    def test_degenerate_matrix_is_refused(self):
        assert sr.plate_scale_from_cd(0.0, 0.0, 0.0, 0.0) is None
        assert sr.plate_scale_from_cd(None, None, None, None) is None
        # Singular (both rows parallel) — no area, no scale.
        assert sr.plate_scale_from_cd(1e-4, 1e-4, 1e-4, 1e-4) is None


class TestGeometry:
    def test_readout_strip_is_refused(self):
        # EU UMa era 80: 8 x 3,211 pixel strips against a ~9 px aperture
        # radius.  The aperture is wider than the image.
        ok, why = sr.geometry_verdict(8, 3211, 8.91)
        assert not ok
        assert "8 px" in why or "short axis 8" in why

    def test_normal_frame_passes(self):
        ok, _ = sr.geometry_verdict(4788, 3194, 8.91)
        assert ok

    def test_unknown_dimensions_refuse_rather_than_assume(self):
        assert not sr.geometry_verdict(None, 3194, 8.91)[0]
        assert not sr.geometry_verdict(4788, 3194, None)[0]

    def test_boundary_is_the_stated_multiple(self):
        aper = 10.0
        need = sr.MIN_AXIS_APERTURE_RADII * aper
        assert sr.geometry_verdict(need, 4000, aper)[0]
        assert not sr.geometry_verdict(need - 1, 4000, aper)[0]


# ---------------------------------------------------------------------------
# Registration policy
# ---------------------------------------------------------------------------
class TestRegistrationPolicy:
    def test_both_ends_must_be_solved(self):
        assert sr.registration_method(True, True) == "wcs"
        assert sr.registration_method(True, False) == "astroalign"
        assert sr.registration_method(False, True) == "astroalign"
        assert sr.registration_method(False, False) == "astroalign"

    def test_a_handful_of_matches_is_not_a_registration(self):
        # 6 matches out of 200 possible is what a WRONG solution produces;
        # believing it would scramble every star identity in the frame.
        assert not sr.wcs_match_ok(6, 200, 300)

    def test_a_good_chain_is_accepted(self):
        assert sr.wcs_match_ok(150, 200, 300)

    def test_denominator_is_the_smaller_catalog(self):
        # A shallow 40-star frame against a 4,000-star reference is judged
        # on the 40 it could have matched, not on the reference's depth.
        assert sr.wcs_match_ok(30, 40, 4000)

    def test_empty_catalogs_fail_closed(self):
        assert not sr.wcs_match_ok(0, 0, 100)
        assert np.isnan(sr.match_rate(0, 0, 100))

    def test_match_rate_matches_the_gate(self):
        assert sr.match_rate(150, 200, 300) == pytest.approx(0.75)


class TestReferencePreference:
    def test_a_deep_solved_candidate_is_promoted(self):
        ranking = [10, 11, 12]
        n = {10: 500, 11: 480, 12: 460}
        solved = {10: 0, 11: 1, 12: 0}
        assert sr.order_reference_candidates(ranking, n, solved)[0] == 11

    def test_a_shallow_solved_candidate_is_not_promoted(self):
        # Astrometry is worth having, but not at the price of two thirds of
        # the star catalog: the reference's depth caps every frame's match.
        ranking = [10, 11, 12]
        n = {10: 500, 11: 100, 12: 460}
        solved = {10: 0, 11: 1, 12: 0}
        assert sr.order_reference_candidates(ranking, n, solved)[0] == 10

    def test_nothing_is_ever_dropped(self):
        # The caller walks the list until one candidate passes double-image
        # QC, so demotion must never mean deletion.
        ranking = [10, 11, 12]
        n = {10: 500, 11: 100, 12: 460}
        solved = {10: 0, 11: 1, 12: 0}
        out = sr.order_reference_candidates(ranking, n, solved)
        assert sorted(out) == sorted(ranking)

    def test_relative_order_survives_promotion(self):
        ranking = [10, 11, 12, 13]
        n = {10: 500, 11: 490, 12: 480, 13: 470}
        solved = {10: 0, 11: 1, 12: 0, 13: 1}
        assert sr.order_reference_candidates(ranking, n, solved) == \
            [11, 13, 10, 12]

    def test_empty_and_unsolved_cases(self):
        assert sr.order_reference_candidates([], {}, {}) == []
        assert sr.order_reference_candidates([1, 2], {1: 5, 2: 4},
                                             {1: 0, 2: 0}) == [1, 2]


# ---------------------------------------------------------------------------
# Master-calibration selection
# ---------------------------------------------------------------------------
class TestMasterSelection:
    def test_exact_exposure_match_required(self):
        assert sr.dark_exptime_matches(64.0, 64.0)
        assert sr.dark_exptime_matches(63.99996, 64.0)
        # 32 s dark for a 64 s exposure: refused, never scaled.  These
        # masters include the bias, and scaling one by an exposure ratio
        # would scale the bias with it.
        assert not sr.dark_exptime_matches(32.0, 64.0)

    def test_missing_values_never_match(self):
        assert not sr.dark_exptime_matches(None, 64.0)
        assert not sr.dark_exptime_matches(64.0, None)
        assert not sr.dark_exptime_matches(float("nan"), 64.0)

    def test_nearest_in_time_wins(self):
        cands = [(2460000.0, "far.fts"), (2460300.0, "near.fts")]
        assert sr.pick_master(cands, 2460310.0)[1] == "near.fts"

    def test_ties_break_on_path_so_reruns_agree(self):
        cands = [(2460300.0, "b.fts"), (2460300.0, "a.fts")]
        assert sr.pick_master(cands, 2460300.0)[1] == "a.fts"

    def test_undated_master_is_usable_but_sorts_last(self):
        cands = [(None, "undated.fts"), (2460300.0, "dated.fts")]
        assert sr.pick_master(cands, 2460300.0)[1] == "dated.fts"
        assert sr.pick_master([(None, "undated.fts")], 2460300.0) is not None

    def test_no_candidates(self):
        assert sr.pick_master([], 2460300.0) is None


# ---------------------------------------------------------------------------
# Series admission
# ---------------------------------------------------------------------------
class TestAdmission:
    def test_single_frame_series_is_refused(self):
        # EU UMa era 79 stages exactly one frame.
        ok, why = sr.series_admission(1)
        assert not ok
        assert "1 matched frame" in why

    def test_healthy_series_admitted(self):
        assert sr.series_admission(500)[0]

    def test_check_star_verdicts(self):
        assert sr.check_star_verdict(4) == "validated"
        assert sr.check_star_verdict(3) == "validated"
        assert sr.check_star_verdict(1) == "weak"
        assert sr.check_star_verdict(0) == "unvalidated"

    def test_mean_ignoring_nan(self):
        assert sr.mean_ignoring_nan([1.0, np.nan, 3.0]) == pytest.approx(2.0)
        assert np.isnan(sr.mean_ignoring_nan([np.nan, np.nan]))
        assert np.isnan(sr.mean_ignoring_nan([]))


# ---------------------------------------------------------------------------
# WCS sidecars and the sky chain
# ---------------------------------------------------------------------------
def _fake_wcs(crval, crpix=(100.0, 100.0), scale_deg=0.000125):
    """A minimal TAN WCS, built in memory (no file, no network)."""
    from astropy.wcs import WCS
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = list(crval)
    w.wcs.crpix = list(crpix)
    w.wcs.cdelt = [-scale_deg, scale_deg]
    return w


class TestSidecarNaming:
    def test_fpack_suffix_is_dropped(self):
        p = reg.sidecar_path(Path("/w"), "rawimage/2025-01-09/f.fts.fz")
        assert p == Path("/w/rawimage/2025-01-09/f.fts.wcs")

    def test_uncompressed_name_is_left_alone(self):
        p = reg.sidecar_path(Path("/w"), "rawimage/2025-01-09/f.fts")
        assert p == Path("/w/rawimage/2025-01-09/f.fts.wcs")

    def test_missing_sidecar_loads_as_none_not_an_exception(self):
        assert reg.load_wcs(Path("/nonexistent/nothing.wcs")) is None


class TestSkyChain:
    def test_identical_solutions_are_the_identity(self):
        w = _fake_wcs((166.4, 25.1))
        xy = np.array([[100.0, 100.0], [250.0, 40.0], [3.0, 900.0]])
        out = reg.chain_to_reference(w, w, xy)
        assert np.allclose(out, xy, atol=1e-6)

    def test_a_shifted_reference_shifts_the_pixels(self):
        # Move the reference's CRPIX by 50 px: every star must land 50 px
        # away, in the same direction, on the reference grid.
        a = _fake_wcs((166.4, 25.1), crpix=(100.0, 100.0))
        b = _fake_wcs((166.4, 25.1), crpix=(150.0, 100.0))
        xy = np.array([[100.0, 100.0], [200.0, 300.0]])
        out = reg.chain_to_reference(a, b, xy)
        assert np.allclose(out[:, 0] - xy[:, 0], 50.0, atol=1e-3)
        assert np.allclose(out[:, 1] - xy[:, 1], 0.0, atol=1e-3)

    def test_chained_positions_match_one_to_one(self):
        # The chain hands its output to the SAME pure matcher astroalign
        # feeds, so the two routes are interchangeable downstream.
        a = _fake_wcs((166.4, 25.1), crpix=(100.0, 100.0))
        b = _fake_wcs((166.4, 25.1), crpix=(103.0, 97.0))
        ref_xy = np.array([[100.0, 100.0], [400.0, 220.0], [50.0, 700.0]])
        # Invent frame pixels that correspond to those reference stars.
        back = reg.chain_to_reference(b, a, ref_xy)
        moved = reg.chain_to_reference(a, b, back)
        idx = ph.match_one_to_one(ref_xy, moved, 1.0)
        assert list(idx) == [0, 1, 2]

    def test_empty_input(self):
        w = _fake_wcs((166.4, 25.1))
        assert reg.chain_to_reference(w, w, np.empty((0, 2))).shape == (0, 2)

    def test_pixel_of_sky_inverts_sky_of_pixel(self):
        w = _fake_wcs((166.4, 25.1))
        sky = reg.sky_of_pixel(w, 137.0, 88.0)
        assert sky is not None
        back = reg.pixel_of_sky(w, *sky)
        assert back == pytest.approx((137.0, 88.0), abs=1e-6)


class TestPixelGrids:
    def test_identical_grids_agree(self):
        xy = np.array([[10.0, 10.0], [200.0, 300.0], [55.0, 4.0]])
        ok, med = reg.pixel_grids_agree(xy, xy)
        assert ok and med == pytest.approx(0.0)

    def test_a_shifted_grid_is_caught(self):
        # This is the test that licensed applying a RAW-frame plate
        # solution to the REDUCED version of the same exposure: on the real
        # pair the median separation was 0.02 px.  A reduction that cropped
        # the frame would look like this instead.
        xy = np.array([[10.0, 10.0], [200.0, 300.0], [55.0, 4.0]])
        ok, med = reg.pixel_grids_agree(xy, xy + 40.0)
        assert not ok and med > 5.0

    def test_empty_catalog_does_not_claim_agreement(self):
        ok, med = reg.pixel_grids_agree(np.empty((0, 2)), np.array([[1.0, 1.0]]))
        assert not ok and np.isnan(med)


# ---------------------------------------------------------------------------
# Local master calibration
# ---------------------------------------------------------------------------
class TestCalibration:
    def test_flat_normalization_preserves_the_count_level(self):
        flat = np.full((8, 8), 4000.0, dtype=np.float32)
        norm, med = cal.normalized_flat(flat)
        assert med == pytest.approx(4000.0)
        assert np.allclose(norm, 1.0)

    def test_dead_flat_pixels_become_nan_not_giant_flux(self):
        flat = np.full((4, 4), 1000.0, dtype=np.float32)
        flat[0, 0] = 1.0                       # a dead pixel
        norm, _ = cal.normalized_flat(flat)
        assert np.isnan(norm[0, 0])
        data = np.full((4, 4), 500.0, dtype=np.float32)
        out, recipe = cal.apply_masters(data, None, flat)
        assert np.isnan(out[0, 0])             # not 500,000 counts of "star"
        assert recipe == "flat_only"

    def test_recipe_names_exactly_what_was_applied(self):
        data = np.full((4, 4), 500.0, dtype=np.float32)
        dark = np.full((4, 4), 100.0, dtype=np.float32)
        flat = np.full((4, 4), 1000.0, dtype=np.float32)
        assert cal.apply_masters(data, dark, flat)[1] == "dark+flat"
        assert cal.apply_masters(data, dark, None)[1] == "dark_only"
        assert cal.apply_masters(data, None, flat)[1] == "flat_only"
        assert cal.apply_masters(data, None, None)[1] == "none"

    def test_dark_subtraction_is_arithmetic_not_approximation(self):
        data = np.full((4, 4), 500.0, dtype=np.float32)
        dark = np.full((4, 4), 120.0, dtype=np.float32)
        out, _ = cal.apply_masters(data, dark, None)
        assert np.allclose(out, 380.0)

    def test_mismatched_master_is_refused(self):
        # A master of the wrong shape means the era assignment is wrong.
        # Cropping or broadcasting it would calibrate with the wrong
        # pixels while reporting success.
        data = np.zeros((4, 4), dtype=np.float32)
        with pytest.raises(ValueError):
            cal.apply_masters(data, np.zeros((8, 8), dtype=np.float32), None)
        with pytest.raises(ValueError):
            cal.apply_masters(data, None, np.ones((8, 8), dtype=np.float32))

    def test_degenerate_flat_is_refused(self):
        with pytest.raises(ValueError):
            cal.normalized_flat(np.zeros((4, 4), dtype=np.float32))


# ---------------------------------------------------------------------------
# The gain policy
# ---------------------------------------------------------------------------
class TestGainPolicy:
    def test_nominal_sits_inside_the_measured_bracket(self):
        lo, hi = sr.GAIN_BRACKET_E_PER_ADU
        assert lo < sr.NOMINAL_GAIN_E_PER_ADU < hi

    def test_bracket_bounds_the_error_scaling(self):
        # The whole justification for one nominal gain: even a worst-case
        # error scales the predicted photon sigma by well under 2, which
        # the empirical inflation factor absorbs.
        lo, hi = sr.GAIN_BRACKET_E_PER_ADU
        assert np.sqrt(hi / lo) < 2.0


# ---------------------------------------------------------------------------
# Translation voting (the dither-aware registration route)
# ---------------------------------------------------------------------------
class TestTranslationVote:
    def _field(self, n=60, seed=3):
        rng = np.random.default_rng(seed)
        return rng.uniform(0, 4000, size=(n, 2))

    def test_recovers_a_pure_dither(self):
        ref = self._field()
        shift = np.array([-161.0, 38.7])       # a real measured EU UMa dither
        frame = ref - shift
        got = reg.vote_translation(frame, ref)
        assert got is not None
        assert np.allclose(got, shift, atol=0.5)

    def test_works_when_the_frame_is_far_shallower(self):
        # THE case astroalign cannot do: 60 frame stars against an
        # 1,839-star reference.  The shift is still one modal offset.
        ref = self._field(n=1839, seed=7)
        shift = np.array([120.0, -75.0])
        frame = ref[:60] - shift
        got = reg.vote_translation(frame, ref)
        assert got is not None and np.allclose(got, shift, atol=1.0)

    def test_survives_noise_and_non_common_stars(self):
        rng = np.random.default_rng(11)
        ref = self._field(n=200, seed=5)
        shift = np.array([40.0, -12.0])
        common = ref[:120] - shift + rng.normal(0, 0.3, size=(120, 2))
        spurious = rng.uniform(0, 4000, size=(80, 2))   # frame-only junk
        frame = np.vstack([common, spurious])
        got = reg.vote_translation(frame, ref)
        assert got is not None and np.allclose(got, shift, atol=1.0)

    def test_two_unrelated_fields_produce_no_proposal(self):
        # A frame pointed somewhere else must NOT be handed a translation:
        # the offsets scatter and no bin can win.
        a = self._field(n=150, seed=1)
        b = self._field(n=150, seed=2)
        assert reg.vote_translation(a, b) is None

    def test_too_few_stars_declines(self):
        assert reg.vote_translation(np.zeros((2, 2)), np.zeros((10, 2))) is None
        assert reg.vote_translation(np.empty((0, 2)), np.zeros((10, 2))) is None

    def test_offset_on_a_bin_edge_is_still_recovered(self):
        # A shift landing exactly on a bin boundary splits its votes; the
        # refinement pass exists precisely so that does not bias the answer.
        ref = self._field(n=120, seed=9)
        shift = np.array([reg.VOTE_BIN_PX * 10.0, reg.VOTE_BIN_PX * 4.0])
        got = reg.vote_translation(ref - shift, ref)
        assert got is not None and np.allclose(got, shift, atol=0.5)


# ---------------------------------------------------------------------------
# The chance-coincidence gate
# ---------------------------------------------------------------------------
class TestChanceGate:
    # A real EU UMa era-80 frame: 4,787 x 3,193 px, 1,839 reference stars,
    # 2 px tolerance.
    AREA = 4787.0 * 3193.0

    def test_a_wrong_shift_matches_almost_nothing(self):
        chance = sr.expected_chance_matches(400, 1839, 2.0, self.AREA)
        assert chance < 1.0

    def test_a_dozen_matches_on_a_cloudy_frame_are_believed(self):
        # THE case the fraction gate wrongly rejected: 400 detections of
        # which only ~60 are stars, so a third of 400 is unreachable — but
        # 19 matches beat a 0.6-pair expectation by thirty times.
        assert not sr.wcs_match_ok(19, 400, 1839)
        assert sr.matches_beat_chance(19, 400, 1839, 2.0, self.AREA)

    def test_the_absolute_floor_still_applies(self):
        # Two coincidences in a sparse field are not a registration, no
        # matter how small the chance expectation is.
        assert not sr.matches_beat_chance(3, 50, 60, 2.0, self.AREA)

    def test_a_crowded_field_demands_more(self):
        # Pack 40,000 stars into a small frame and chance matches become
        # common, so the same 19 matches no longer prove anything.
        small = 500.0 * 500.0
        assert not sr.matches_beat_chance(19, 400, 40000, 2.0, small)

    def test_degenerate_inputs_fail_closed(self):
        assert not sr.matches_beat_chance(100, 400, 1839, 2.0, 0.0)
        assert not sr.matches_beat_chance(100, 400, 1839, 0.0, self.AREA)
