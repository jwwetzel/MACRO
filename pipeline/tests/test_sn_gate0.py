"""Unit tests for macro_sn.gate0 — the SN 2023ixf Gate 0 pure logic.

Same philosophy as test_astrom / test_dispersion: every decision function is
exercised on hand-built cases, and the cases that must NOT pass are tested as
hard as the ones that must.  The rules under test are the ones a referee
would attack:

* the imaging gate must exclude ONLY on a positive dispersion measurement,
  because the retired filter-label gate is exactly what contaminated the
  published rate this stage was built to repair;
* a box maximum must never be allowed to masquerade as a supernova peak;
* the isolation calibration must not read its own "nothing brighter within
  the box" floor as a confusion;
* the venue rule must be the strategy's pre-registered disjunction and
  nothing more generous.
"""

from __future__ import annotations

import math

import pytest

from macro_sn import gate0 as g0
from macro_sn.gate0 import (
    GrismSeries, Screen, band_role, dispersion_class, epoch_role,
    gnomonic_pixel, grism_promotion, is_direct_image, is_usable_photometry,
    isolation_false_id_rate, parse_sexagesimal, phase_days, position_quality,
    rate_pct, saturation_class, screen_for_mode, venue_posture)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def screen() -> Screen:
    """The campaign's real screen: High Gain, S2-measured clip 3,496 ADU."""
    return screen_for_mode("High Gain", clip_adu=3496, veto_adu=3200)


def frame(**over) -> dict:
    """A template usable-photometry row: campaign, rawimage, G, clean."""
    base = dict(epoch_role="campaign", tree="rawimage", filter="G",
                dispersion_verdict=None, saturation_class="clean")
    base.update(over)
    return base


# ---------------------------------------------------------------------------
# the screen: derived from the measurement, never typed
# ---------------------------------------------------------------------------
class TestScreen:

    def test_screen_reproduces_the_strategys_hand_numbers(self):
        """The strategy typed 2,800 / 2,400 ADU from an ASSUMED ~3,500 clip.
        Applying the same fractions to S2's MEASURED 3,496 must land on
        those numbers — that agreement is the check that the assumption was
        sound, and it is the reason the fractions are stored rather than the
        ADU values."""
        s = screen()
        assert s.reject_adu == 2797            # strategy said 2,800
        assert s.suspect_adu == 2397           # strategy said 2,400
        assert abs(s.reject_adu - 2800) / 2800 < 0.002
        assert abs(s.suspect_adu - 2400) / 2400 < 0.002

    def test_screen_refuses_to_invent_a_clip(self):
        """A census that silently fell back to a typed screen would publish
        a threshold with nothing behind it — which is the precise state that
        had task SN-G0b marked BLOCKED."""
        with pytest.raises(ValueError, match="no measured clip"):
            screen_for_mode("High Gain", clip_adu=None, veto_adu=3200)

    def test_a_different_mode_gets_a_different_screen(self):
        """The post-fade templates were taken in Fast/1MHz modes whose clip
        is 16-bit; applying the High Gain screen to them would reject every
        frame."""
        fast = screen_for_mode("Fast", clip_adu=65534, veto_adu=60200)
        assert fast.reject_adu > 50000
        assert fast.reject_adu != screen().reject_adu


# ---------------------------------------------------------------------------
# the imaging gate — the correction this stage carries
# ---------------------------------------------------------------------------
class TestImagingGate:

    def test_only_a_measured_spectrum_is_excluded(self):
        assert not is_direct_image("dispersed")
        assert is_direct_image("direct")
        assert is_direct_image("indeterminate")
        assert is_direct_image(None)
        assert is_direct_image("")

    def test_a_grism_label_does_not_exclude_a_measured_image(self):
        """Slot '6' is a MIXED slot: on this target three of its frames are
        measured direct images.  The retired label rule deleted them unseen;
        the measured rule must keep them."""
        assert is_direct_image("direct") is True

    def test_indeterminate_is_kept_not_dropped(self):
        """Exclusion has to be EARNED by a positive measurement.  S2c's
        commonest indeterminate reason is 'no usable sources extracted',
        which describes a blank image at least as often as a spectrum."""
        assert is_direct_image("indeterminate") is True

    def test_dispersion_class_normalises_the_absence_of_a_verdict(self):
        assert dispersion_class(None) == "unmeasured"
        assert dispersion_class("  ") == "unmeasured"
        assert dispersion_class("DISPERSED") == "dispersed"

    def test_an_unknown_verdict_is_reported_verbatim(self):
        """A new S2c verdict must show up in the census table, not disappear
        into a known bucket."""
        assert dispersion_class("smeared") == "smeared"


# ---------------------------------------------------------------------------
# saturation classification
# ---------------------------------------------------------------------------
class TestSaturationClass:

    def test_measured_peaks_get_the_three_measured_classes(self):
        s = screen()
        assert saturation_class(1000.0, "wcs", s) == "clean"
        assert saturation_class(2500.0, "wcs", s) == "suspect"
        assert saturation_class(3300.0, "wcs", s) == "rejected"

    def test_thresholds_are_inclusive_at_the_bottom(self):
        s = screen()
        assert saturation_class(s.suspect_adu, "wcs", s) == "suspect"
        assert saturation_class(s.suspect_adu - 1, "wcs", s) == "clean"
        assert saturation_class(s.reject_adu, "wcs", s) == "rejected"

    def test_a_bound_can_only_conclude_in_one_direction(self):
        """This is the whole point of the fifth class.  Under the box the
        supernova is provably unsaturated wherever it fell; over it, nothing
        follows — the bright thing may be a field star."""
        s = screen()
        assert saturation_class(500.0, "bound", s) == "bounded_clean"
        assert saturation_class(3400.0, "bound", s) == "undetermined"
        assert saturation_class(3400.0, "wcs_off", s) == "undetermined"

    def test_a_bound_is_never_called_saturated(self):
        """Folding 'undetermined' into 'rejected' would manufacture
        saturated epochs out of field stars."""
        s = screen()
        for q in ("bound", "wcs_off"):
            assert saturation_class(9999.0, q, s) != "rejected"

    def test_no_measurement_is_undetermined(self):
        assert saturation_class(None, "wcs", screen()) == "undetermined"


# ---------------------------------------------------------------------------
# position quality
# ---------------------------------------------------------------------------
class TestPositionQuality:

    def test_a_centred_peak_on_a_solved_frame_is_a_measurement(self):
        assert position_quality(0.4, has_wcs=True) == "wcs"

    def test_an_offset_peak_on_a_solved_frame_is_only_a_bound(self):
        """The pre-explosion template epoch has a 14-16 px FWHM, so its
        blob maxima wander well off the catalogue position."""
        assert position_quality(11.0, has_wcs=True) == "wcs_off"

    def test_no_solution_is_always_a_bound(self):
        """Never promoted on the strength of finding something bright: in
        narrowband a field star outshines the supernova often enough to make
        that inference wrong."""
        assert position_quality(0.0, has_wcs=False) == "bound"
        assert position_quality(None, has_wcs=False) == "bound"


# ---------------------------------------------------------------------------
# the composite usability test
# ---------------------------------------------------------------------------
class TestUsablePhotometry:

    def test_the_template_row_passes(self):
        assert is_usable_photometry(frame()) is True

    def test_a_measured_spectrum_is_not_photometry(self):
        assert not is_usable_photometry(frame(dispersion_verdict="dispersed"))

    def test_a_narrowband_frame_is_not_broadband_photometry(self):
        assert not is_usable_photometry(frame(filter="H"))

    def test_a_template_epoch_is_not_campaign_photometry(self):
        assert not is_usable_photometry(frame(epoch_role="template_pre"))

    def test_a_suspect_frame_is_held_back(self):
        """'suspect' is conditional on a linearity curve Step 2 has not
        produced, so it may not be counted as usable today."""
        assert not is_usable_photometry(frame(saturation_class="suspect"))

    def test_a_bounded_clean_frame_counts(self):
        assert is_usable_photometry(frame(saturation_class="bounded_clean"))

    def test_an_engineering_frame_is_excluded_by_tree(self):
        """Two canonical frames on this sky live under mjc/misc/neg10_test/
        and are detector tests.  One reads 65,535 ADU at the supernova's
        position in a 12-bit channel, which is how they were caught."""
        assert not is_usable_photometry(frame(tree="mjc"))


# ---------------------------------------------------------------------------
# isolation calibration
# ---------------------------------------------------------------------------
class TestIsolation:

    def test_the_floor_is_not_a_confusion(self):
        """When nothing in the box outbrightens the supernova the census
        records the box half-width as a FLOOR.  Counting those as
        confusions would turn the cleanest possible result into a 100%
        false-identification rate."""
        wrong, tested = isolation_false_id_rate([250.0, 250.0, 250.0], 250.0)
        assert (wrong, tested) == (0, 3)

    def test_a_genuinely_closer_star_is_a_confusion(self):
        wrong, tested = isolation_false_id_rate([12.0, 250.0, 30.0], 250.0)
        assert (wrong, tested) == (2, 3)

    def test_unmeasured_rows_are_not_counted_either_way(self):
        wrong, tested = isolation_false_id_rate([None, None, 10.0], 250.0)
        assert (wrong, tested) == (1, 1)

    def test_a_tighter_radius_is_a_weaker_claim(self):
        iso = [12.0, 40.0, 250.0]
        assert isolation_false_id_rate(iso, 20.0)[0] == 1
        assert isolation_false_id_rate(iso, 100.0)[0] == 2


# ---------------------------------------------------------------------------
# coordinates
# ---------------------------------------------------------------------------
class TestCoordinates:

    def test_sexagesimal_hours_and_degrees(self):
        assert parse_sexagesimal("14 03 38.562", True) == pytest.approx(
            210.910675, abs=1e-6)
        assert parse_sexagesimal("+54:18:41.94", False) == pytest.approx(
            54.311650, abs=1e-6)

    def test_a_negative_declination_inside_the_first_degree(self):
        """The sign has to come from the STRING: -00 30 00 has a zero
        degrees field whose float carries no sign at all."""
        assert parse_sexagesimal("-00 30 00", False) == pytest.approx(-0.5)

    def test_missing_fields_pad_with_zero(self):
        assert parse_sexagesimal("14 03", True) == pytest.approx(14.05 * 15)

    def test_unparseable_pointing_returns_none_rather_than_raising(self):
        """A missing pointing card is a normal state in this archive and
        must not abort a 1,461-frame pass."""
        assert parse_sexagesimal(None, True) is None
        assert parse_sexagesimal("", True) is None
        assert parse_sexagesimal("no such thing", True) is None

    def test_the_tangent_point_lands_on_the_reference_pixel(self):
        x, y = gnomonic_pixel(210.0, 54.0, 210.0, 54.0, 1.5e-4, 2048, 2048)
        assert x == pytest.approx(2048.0)
        assert y == pytest.approx(2048.0)

    def test_east_is_left_and_north_is_up(self):
        """A source at higher RA must fall at LOWER x; a source at higher
        declination at HIGHER y.  Getting this backwards would put every
        search box on the wrong side of the field."""
        x, _ = gnomonic_pixel(210.01, 54.0, 210.0, 54.0, 1.5e-4, 2048, 2048)
        _, y = gnomonic_pixel(210.0, 54.01, 210.0, 54.0, 1.5e-4, 2048, 2048)
        assert x < 2048.0
        assert y > 2048.0

    def test_the_scale_is_honoured(self):
        """0.01 deg north at 1.5e-4 deg/px is 66.7 px."""
        _, y = gnomonic_pixel(210.0, 54.01, 210.0, 54.0, 1.5e-4, 2048, 2048)
        assert y - 2048.0 == pytest.approx(0.01 / 1.5e-4, rel=1e-3)

    def test_the_far_side_of_the_sky_raises(self):
        with pytest.raises(ValueError):
            gnomonic_pixel(30.0, -54.0, 210.0, 54.0, 1.5e-4, 2048, 2048)


# ---------------------------------------------------------------------------
# roles and phase
# ---------------------------------------------------------------------------
class TestRoles:

    def test_band_roles(self):
        assert band_role("G") == "broadband"
        assert band_role("H") == "narrowband"
        assert band_role("6") == "other"
        assert band_role(None) == "other"

    def test_epoch_roles(self):
        assert epoch_role("2023-05-04") == "template_pre"
        assert epoch_role("2023-05-19") == "campaign"
        assert epoch_role("2023-07-06") == "campaign"
        assert epoch_role("2024-05-18") == "template_post"

    def test_phase_is_measured_from_the_adopted_explosion_epoch(self):
        """Night 2023-05-23 (UT 05-24) is the strategy's '+5.4 d' clean
        start; the census must agree with that arithmetic."""
        jd = g0.MJD_OFFSET + g0.T0_MJD + 5.4
        assert phase_days(jd) == pytest.approx(5.4)

    def test_phase_of_a_missing_time_is_none_not_zero(self):
        assert phase_days(None) is None

    def test_the_adopted_t0_sits_inside_its_published_bracket(self):
        assert g0.T0_MJD_LOW <= g0.T0_MJD <= g0.T0_MJD_HIGH


# ---------------------------------------------------------------------------
# rates
# ---------------------------------------------------------------------------
class TestRates:

    def test_rate_pct(self):
        assert rate_pct(18, 48) == pytest.approx(37.5)

    def test_an_empty_denominator_is_none_not_zero(self):
        """A band with no frames has no rate; reporting 0% would read as a
        total failure rather than as an absence of data."""
        assert rate_pct(0, 0) is None

    def test_the_verdict_thresholds_are_s1s_own(self):
        """Gate 0's astrometry answer is only comparable with the S1 stratum
        verdict if both are judged by the same bar.  A second copy of the
        thresholds here is exactly what must not exist."""
        from macro_core import astrom
        assert g0.astrometry_verdict(48, 48, 48) == astrom.verdict_for(
            48, 48, 48)
        assert g0.astrometry_verdict(18, 48) == "NO-GO"
        assert g0.astrometry_verdict(553, 630, 630) == "GO"


# ---------------------------------------------------------------------------
# Gate 0c — the pre-registered promotion rule
# ---------------------------------------------------------------------------
def series(**over) -> GrismSeries:
    base = dict(n_labelled=83, n_dispersed=61, n_direct=3, n_indeterminate=19,
                n_nights=13, n_flash_nights=3, n_nights_with_paired_direct=13,
                n_extracted=0, n_contamination_passed=0,
                wavelength_source="", n_flats=0)
    base.update(over)
    return GrismSeries(**base)


class TestGrismPromotion:

    def test_the_real_series_does_not_promote(self):
        """61 measured spectra over 13 nights is a real series and still
        fails the bar, because the bar is about EXTRACTED, contamination-
        tested, wavelength-calibrated spectra and none exist."""
        d = grism_promotion(series())
        assert d["promoted"] is False
        assert set(d["blocking"]) == {"extracted", "contamination",
                                      "wavelength"}

    def test_frames_alone_never_promote(self):
        """The failure mode this rule exists to prevent: promoting a product
        because the frames for it exist."""
        d = grism_promotion(series(n_nights=40, n_flash_nights=10))
        assert d["promoted"] is False

    def test_a_complete_triage_promotes(self):
        d = grism_promotion(series(n_extracted=4, n_contamination_passed=3,
                                   wavelength_source="2023-06 arc lamp"))
        assert d["promoted"] is True
        assert d["blocking"] == []

    def test_a_failed_contamination_test_blocks_everything_else(self):
        """The offset-trace test is mandatory and independent: ambient
        H II-region Halpha must not be promotable."""
        d = grism_promotion(series(n_extracted=5, n_contamination_passed=0,
                                   wavelength_source="self-calibrated"))
        assert d["promoted"] is False
        assert "contamination" in d["blocking"]

    def test_too_few_nights_blocks_even_a_perfect_extraction(self):
        d = grism_promotion(series(n_nights=2, n_extracted=5,
                                   n_contamination_passed=5,
                                   wavelength_source="arc"))
        assert d["promoted"] is False
        assert "nights" in d["blocking"]

    def test_flash_coverage_is_context_and_never_a_substitute(self):
        """Covering the flash phase is why the product was proposed, not a
        reason to publish it."""
        d = grism_promotion(series(n_flash_nights=10))
        assert d["clauses"]["flash_nights"]["passed"] is True
        assert d["promoted"] is False


# ---------------------------------------------------------------------------
# the venue rule
# ---------------------------------------------------------------------------
class TestVenue:

    def test_the_base_case_holds_when_neither_condition_is_met(self):
        v = venue_posture(grism_promoted=False, bandpass_recovered=False)
        assert v["posture"] == g0.VENUE_BASE
        assert v["moved"] is False

    def test_either_condition_alone_promotes(self):
        assert venue_posture(True, False)["posture"] == g0.VENUE_UPSIDE
        assert venue_posture(False, True)["posture"] == g0.VENUE_UPSIDE

    def test_the_rule_is_a_disjunction_and_nothing_more_generous(self):
        """The strategy pre-registered exactly two routes to ApJ.  A third
        one appearing here would be the drafting-stage rationalisation the
        posture was decided in advance to prevent."""
        assert venue_posture(False, False)["moved"] is False
