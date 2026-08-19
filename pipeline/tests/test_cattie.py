"""Unit tests for ``macro_phot.cattie`` — the catalogue-tie arithmetic.

Every test states the fact about PHOTOMETRY it is protecting, not just the
fact about the code.  Synthetic inputs whose right answer is known in closed
form, and at least one test per function that would fail if the function
were quietly replaced by something plausible-but-wrong.

The four committee rulings are each defended by name below:
  1. never transform the target        -> test_target_colour_* (placement
                                          only; no magnitude is produced)
  2. solve on comparison stars         -> enforced by the caller; the split
                                          machinery is test_holdout_*
  3. ZP AND colour term, with errors
     and a stated valid range          -> test_robust_line_fit_*,
                                          test_colour_range_*
  4. saturation/blend veto catalogue
     stars too                         -> test_clean_mask_*
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from macro_phot import cattie as ct


# ==========================================================================
# 1 - band mapping: which catalogue column does a FILTER label mean?
# ==========================================================================

def test_every_campaign_filter_has_a_band_or_is_declared_untieable():
    """The six filter labels this archive actually uses must each either
    map to a catalogue band or be absent from the map ON PURPOSE — a label
    that silently falls through would be tied to nothing and reported as
    tied to something."""
    for f in ("g", "r", "i", "z", "y", "G", "R", "I"):
        assert ct.band_candidates(f), f"filter {f!r} has no band hypothesis"


def test_y_band_is_served_only_by_the_secondary_catalogue():
    """ATLAS-REFCAT2 stops at z.  If this ever returns a REFCAT2 spec for
    y, some column has been mistaken for a y magnitude."""
    assert ct.band_candidates("y", "refcat2") == ()
    assert ct.primary_band("y", "gaia_gspc").mag_col == "y_ps1_mag"


def test_uppercase_labels_carry_a_johnson_alternative_to_be_tested():
    """Eras 6/7 write 'G','R','I' where later eras write 'g','r','i'.  The
    primary hypothesis is that the glass is the same; the ALTERNATIVE must
    exist so the data can reject it rather than the analyst assume it."""
    alts = ct.band_candidates("R", hypothesis="alternative")
    assert any(b.mag_col == "r_jkc_mag" for b in alts)
    assert ct.primary_band("R", "refcat2").mag_col == "rmag"


def test_colour_index_is_the_conventional_neighbour_pair():
    """g and r are fitted against g-r, i against r-i, z against i-z: the
    colour has to straddle the band or the fit has no leverage on it."""
    assert ct.primary_band("g", "refcat2").colour_label == "g-r"
    assert ct.primary_band("i", "refcat2").colour_label == "r-i"
    assert ct.primary_band("z", "refcat2").colour_label == "i-z"


# ==========================================================================
# 2 - the train / holdout split
# ==========================================================================

def test_holdout_mask_is_deterministic():
    """A published accuracy number that changes when you re-run the script
    is not a measurement."""
    ids = list(range(200))
    a = ct.holdout_mask(ids, "stlmi|e7|G")
    b = ct.holdout_mask(ids, "stlmi|e7|G")
    assert np.array_equal(a, b)


def test_holdout_mask_hits_about_the_requested_fraction():
    ids = list(range(4000))
    m = ct.holdout_mask(ids, "vvpup|e76|g", fraction=0.25)
    assert 0.22 < m.mean() < 0.28


def test_holdout_mask_differs_between_series():
    """A star unlucky in g must not be systematically unlucky in r as well,
    or the 'independent' check stars are one sample wearing three hats."""
    ids = list(range(500))
    g = ct.holdout_mask(ids, "stlmi|e76|g")
    r = ct.holdout_mask(ids, "stlmi|e76|r")
    assert not np.array_equal(g, r)
    # ... and the overlap must look like chance, not like a shared rule.
    assert 0.02 < float(np.mean(g & r)) < 0.12


def test_holdout_mask_zero_fraction_holds_nothing_out():
    assert not ct.holdout_mask(list(range(50)), "k", fraction=0.0).any()


# ==========================================================================
# 3 - the cleanliness gate (ruling 4)
# ==========================================================================

def _clean_kwargs(n: int) -> dict:
    """n perfectly clean candidate stars — the baseline every test spoils."""
    return dict(
        saturated_frac=np.zeros(n), near_veto_frac=np.zeros(n),
        blend_sep_arcsec=np.full(n, 30.0), blend_dmag=np.full(n, 5.0),
        annulus_sep_arcsec=np.full(n, 30.0), annulus_dmag=np.full(n, 5.0),
        second_sep_arcsec=np.full(n, 20.0), match_sep_arcsec=np.full(n, 0.3),
        cat_mag=np.full(n, 16.0), cat_colour=np.full(n, 0.5),
        cat_flag=np.zeros(n))


def test_clean_mask_keeps_clean_stars():
    keep, census = ct.clean_mask(**_clean_kwargs(10))
    assert keep.all()
    assert census.n_clean == 10


def test_clean_mask_rejects_any_saturation_at_all():
    """A comparison star that clips on the best-seeing frames has a mean
    magnitude biased faint by an unrecoverable amount.  Tolerance is zero
    because there are always more candidates than the fit needs."""
    kw = _clean_kwargs(5)
    kw["saturated_frac"] = np.array([0.0, 0.001, 0.0, 0.5, 0.0])
    keep, census = ct.clean_mask(**kw)
    assert list(keep) == [True, False, True, False, True]
    assert census.n_saturated == 2


def test_clean_mask_rejects_the_non_linear_shoulder_below_the_clip():
    """High Gain's scale ends at 3,496 ADU; the detector is already
    non-linear at 0.9 of the applied veto, and nothing flags that."""
    kw = _clean_kwargs(3)
    kw["near_veto_frac"] = np.array([0.0, 0.02, 0.0])
    keep, _ = ct.clean_mask(**kw)
    assert list(keep) == [True, False, True]


def test_clean_mask_rejects_a_bright_neighbour_inside_the_aperture():
    """4-arcsec aperture: a neighbour at 3 arcsec that is only 1 mag fainter
    puts 40% of somebody else's light into this star's flux."""
    kw = _clean_kwargs(3)
    kw["blend_sep_arcsec"] = np.array([30.0, 3.0, 3.0])
    kw["blend_dmag"] = np.array([5.0, 1.0, 6.0])
    keep, census = ct.clean_mask(**kw)
    # star 1 blended; star 2 has a neighbour equally close but 6 mag fainter
    # (0.4% of the flux), which is tolerated on purpose.
    assert list(keep) == [True, False, True]
    assert census.n_blend_aperture == 1


def test_clean_mask_rejects_a_brighter_neighbour_in_the_sky_annulus():
    """The 8-12 arcsec background annulus: a brighter star inside it biases
    the subtracted sky, which biases the flux, which biases the tie."""
    kw = _clean_kwargs(2)
    kw["annulus_sep_arcsec"] = np.array([30.0, 10.0])
    kw["annulus_dmag"] = np.array([5.0, -2.0])
    keep, census = ct.clean_mask(**kw)
    assert list(keep) == [True, False]
    assert census.n_blend_annulus == 1


def test_clean_mask_rejects_an_ambiguous_identification():
    """Two catalogue sources at comparable distance means we do not know
    which star we measured; a coin-flip identity is worse than no tie star."""
    kw = _clean_kwargs(2)
    kw["match_sep_arcsec"] = np.array([0.3, 0.5])
    kw["second_sep_arcsec"] = np.array([20.0, 0.8])   # 0.8 < 2.5 * 0.5
    keep, census = ct.clean_mask(**kw)
    assert list(keep) == [True, False]
    assert census.n_ambiguous == 1


def test_clean_mask_census_partitions_every_candidate_exactly_once():
    """The arithmetic that proves no star was silently lost: every
    rejection reason plus the survivors must add up to the candidates."""
    kw = _clean_kwargs(8)
    kw["saturated_frac"][0] = 1.0
    kw["near_veto_frac"][1] = 1.0
    kw["blend_sep_arcsec"][2] = 2.0
    kw["blend_dmag"][2] = 0.0
    kw["annulus_sep_arcsec"][3] = 9.0
    kw["annulus_dmag"][3] = -1.0
    kw["second_sep_arcsec"][4] = 0.4
    kw["cat_mag"][5] = np.nan
    kw["cat_flag"][6] = 1
    _, c = ct.clean_mask(**kw)
    total = (c.n_saturated + c.n_near_veto + c.n_blend_aperture
             + c.n_blend_annulus + c.n_ambiguous + c.n_no_cat_mag
             + c.n_flagged + c.n_clean)
    assert total == c.n_candidates == 8
    assert c.n_clean == 1


def test_clean_mask_applies_the_same_rules_whatever_the_catalogue():
    """There is deliberately no catalogue argument.  If ATLAS-REFCAT2 and
    Gaia were vetoed by different rules, the cross-catalogue difference
    would measure the rules, not the catalogues."""
    import inspect
    params = set(inspect.signature(ct.clean_mask).parameters)
    assert "catalogue" not in params and "cat" not in params


# ==========================================================================
# 4 - matching and blend metrology
# ==========================================================================

def test_match_by_sky_finds_the_nearest_source_within_tolerance():
    star_ra, star_dec = [10.0], [20.0]
    # 0.5 arcsec east, and 10 arcsec east
    cra = [10.0 + 0.5 / 3600 / math.cos(math.radians(20)),
           10.0 + 10.0 / 3600 / math.cos(math.radians(20))]
    cdec = [20.0, 20.0]
    idx, sep, sep2 = ct.match_by_sky(star_ra, star_dec, cra, cdec)
    assert idx[0] == 0
    assert sep[0] == pytest.approx(0.5, abs=0.02)
    assert sep2[0] == pytest.approx(10.0, abs=0.05)


def test_match_by_sky_refuses_a_source_beyond_tolerance():
    """1.2 arcsec is 1.5-2.7 pixels at these plate scales.  An identity may
    be probable; a PHOTOMETRIC tie must be certain."""
    idx, _, _ = ct.match_by_sky([10.0], [20.0], [10.001], [20.0])
    assert idx[0] == -1


def test_match_by_sky_survives_empty_inputs():
    idx, sep, sep2 = ct.match_by_sky([], [], [1.0], [2.0])
    assert len(idx) == 0
    idx, _, _ = ct.match_by_sky([1.0], [2.0], [], [])
    assert idx[0] == -1


def test_neighbour_metrics_measures_separation_and_contrast():
    """One star at the origin with a companion 4 arcsec away, 2 mag
    fainter, and a bright star 20 arcsec away."""
    cosd = math.cos(math.radians(0.0))
    cra = [0.0, 4.0 / 3600 / cosd, 20.0 / 3600 / cosd]
    cdec = [0.0, 0.0, 0.0]
    cmag = [15.0, 17.0, 12.0]
    out = ct.neighbour_metrics([0.0], [0.0], [15.0], cra, cdec, cmag,
                               self_index=np.array([0]))
    assert out["nn_sep"][0] == pytest.approx(4.0, abs=0.05)
    assert out["nn_dmag"][0] == pytest.approx(2.0, abs=1e-6)
    assert out["bright_sep"][0] == pytest.approx(20.0, abs=0.2)
    assert out["bright_dmag"][0] == pytest.approx(-3.0, abs=1e-6)


def test_neighbour_metrics_excludes_the_star_from_its_own_search():
    """Without ``self_index`` every star is its own nearest neighbour at
    zero separation and the blend veto deletes the entire catalogue."""
    out = ct.neighbour_metrics([0.0], [0.0], [15.0], [0.0], [0.0], [15.0],
                               self_index=np.array([0]))
    assert not np.isfinite(out["nn_sep"][0])


# ==========================================================================
# 5 - the robust line fit (ruling 3)
# ==========================================================================

def test_robust_line_fit_recovers_a_planted_zero_point_and_colour_term():
    """The REQUIRED synthetic-recovery test: plant ZP and k, demand them
    back.  This is the test that would fail if the fit were quietly
    replaced by a plain median offset."""
    rng = np.random.default_rng(7)
    colour = rng.uniform(0.2, 1.4, 300)
    truth_zp, truth_k, c_ref = -3.25, 0.117, 0.8
    y = truth_zp + truth_k * (colour - c_ref) + rng.normal(0, 0.01, 300)
    fit = ct.robust_line_fit(colour, y, np.full(300, 0.01), x_ref=c_ref)
    assert fit.zp == pytest.approx(truth_zp, abs=0.002)
    assert fit.slope == pytest.approx(truth_k, abs=0.006)
    assert fit.resid_rms == pytest.approx(0.01, rel=0.2)
    assert fit.converged


def test_robust_line_fit_is_not_dragged_by_a_variable_star():
    """A comparison sample always contains a few variables and blends.  An
    ordinary least-squares fit lets one of them tilt a PUBLISHED colour
    term; Huber weighting plus a 4-sigma clip must not."""
    rng = np.random.default_rng(11)
    colour = rng.uniform(0.2, 1.4, 120)
    y = -3.0 + 0.10 * (colour - 0.8) + rng.normal(0, 0.008, 120)
    y[0] += 1.5                          # one grossly wrong star
    y[1] -= 0.9                          # and one more
    fit = ct.robust_line_fit(colour, y, np.full(120, 0.008), x_ref=0.8)
    assert fit.zp == pytest.approx(-3.0, abs=0.004)
    assert fit.slope == pytest.approx(0.10, abs=0.01)
    assert fit.n_clipped >= 2


def test_robust_line_fit_centres_on_the_median_colour_by_default():
    x = np.array([0.0, 1.0, 2.0, 3.0, 10.0])
    fit = ct.robust_line_fit(x, np.zeros(5))
    assert fit.x_ref == pytest.approx(2.0)


def test_centring_decorrelates_the_zero_point_from_the_colour_term():
    """Why x_ref exists.  Fit the SAME data centred and uncentred: the
    zero-point error must be smaller when centred, because an uncentred
    intercept inherits the slope's uncertainty over the whole lever arm."""
    rng = np.random.default_rng(3)
    colour = rng.uniform(0.6, 1.4, 80)          # far from zero
    y = -3.0 + 0.3 * (colour - 1.0) + rng.normal(0, 0.02, 80)
    centred = ct.robust_line_fit(colour, y, np.full(80, 0.02))
    at_zero = ct.robust_line_fit(colour, y, np.full(80, 0.02), x_ref=0.0)
    assert centred.zp_err < at_zero.zp_err
    # ... and the two describe the same line.
    assert (at_zero.zp + at_zero.slope * centred.x_ref) == \
        pytest.approx(centred.zp, abs=1e-6)


def test_robust_line_fit_inflates_errors_when_the_scatter_exceeds_the_model():
    """Catalogue errors understate the real star-to-star scatter (bandpass
    mismatch beyond a linear term is a per-star effect).  Quoting the
    formal error of an obviously under-dispersed model would understate the
    zero point by exactly the factor the fit itself measured."""
    rng = np.random.default_rng(5)
    x = rng.uniform(0, 1, 200)
    y = rng.normal(0, 0.05, 200)                 # true scatter 0.05
    honest = ct.robust_line_fit(x, y, np.full(200, 0.05))
    optimistic = ct.robust_line_fit(x, y, np.full(200, 0.005))  # claimed 10x
    assert optimistic.chi2nu > 10
    assert optimistic.zp_err == pytest.approx(honest.zp_err, rel=0.35)


def test_robust_line_fit_returns_nan_rather_than_a_guess_when_starved():
    fit = ct.robust_line_fit([1.0, 2.0], [0.0, 1.0])
    assert math.isnan(fit.zp) and not fit.converged


def test_robust_line_fit_tolerates_missing_errors():
    """A star with no usable error must cost it no more and no less than an
    average star's vote — not silently drop out of the sample."""
    x = np.linspace(0, 1, 30)
    y = 2.0 + 0.5 * (x - 0.5)
    e = np.full(30, 0.01)
    e[3] = np.nan
    e[7] = 0.0
    fit = ct.robust_line_fit(x, y, e, x_ref=0.5)
    assert fit.n_used == 30
    assert fit.zp == pytest.approx(2.0, abs=1e-6)


# ==========================================================================
# 6 - colour range, and where a target sits in it (ruling 3)
# ==========================================================================

def test_colour_range_reports_span_and_core():
    c = list(np.linspace(0.0, 1.0, 101))
    cmin, cmax, p05, p95 = ct.colour_range(c)
    assert (cmin, cmax) == (0.0, 1.0)
    assert p05 == pytest.approx(0.05, abs=1e-6)
    assert p95 == pytest.approx(0.95, abs=1e-6)


def test_colour_position_names_the_three_regimes_and_the_unknown():
    assert ct.colour_position(0.5, 0.0, 1.0, 0.1, 0.9) == "inside-core"
    assert ct.colour_position(0.95, 0.0, 1.0, 0.1, 0.9) == "inside-span"
    assert ct.colour_position(1.4, 0.0, 1.0, 0.1, 0.9) == "extrapolated"
    assert ct.colour_position(float("nan"), 0.0, 1.0, 0.1, 0.9) == "unknown"


def test_extrapolation_error_is_zero_inside_and_grows_outside():
    """A CV bluer than every tie star is not calibrated by that fit, and the
    size of the lie is the colour term times how far outside it sits."""
    assert ct.colour_extrapolation_error(0.5, 0.0, 1.0, 0.2, 0.01) == 0.0
    e = ct.colour_extrapolation_error(-0.3, 0.0, 1.0, 0.2, 0.01)
    assert e == pytest.approx(0.3 * 0.2 + 0.3 * 0.01, abs=1e-9)


def test_target_colour_solve_inverts_the_two_band_relation_exactly():
    """Round trip: take a catalogue colour C, produce the natural-system
    colour the telescope would measure, and demand C back.  This function
    places the target on the colour axis; it never produces a magnitude,
    which is what keeps ruling 1 intact."""
    c_true, k_b, k_r, c_ref = 0.35, 0.12, -0.05, 0.80
    nat = c_true + (k_b - k_r) * (c_true - c_ref)
    assert ct.target_colour_solve(nat, k_b, k_r, c_ref) == \
        pytest.approx(c_true, abs=1e-12)


def test_target_colour_solve_refuses_a_collapsed_colour_baseline():
    """If the two colour terms differ by exactly -1 the instrument carries
    no colour information at all; NaN is the honest answer, not a number."""
    assert math.isnan(ct.target_colour_solve(0.4, -0.5, 0.5, 0.8))


def test_target_colour_solve_is_identity_when_both_terms_vanish():
    assert ct.target_colour_solve(0.42, 0.0, 0.0, 0.9) == pytest.approx(0.42)


# ==========================================================================
# 7 - the adversarial checks (ruling 5)
# ==========================================================================

def test_residual_trend_recovers_a_planted_non_linearity():
    """If the detector were non-linear, the tie residual would slope
    against catalogue magnitude.  Plant 5 mmag per magnitude and demand
    it back, with the SWING that a reader actually cares about."""
    rng = np.random.default_rng(2)
    mag = rng.uniform(13.0, 18.0, 200)
    resid = 0.005 * (mag - 15.0) + rng.normal(0, 0.002, 200)
    t = ct.residual_trend("cat_mag", mag, resid)
    assert t.slope == pytest.approx(0.005, abs=0.0008)
    assert t.significant
    assert t.swing == pytest.approx(0.005 * t.span, rel=0.2)


def test_residual_trend_calls_pure_noise_insignificant():
    """The test has to be able to say NO, or its yes means nothing."""
    rng = np.random.default_rng(4)
    x = rng.uniform(13, 18, 200)
    t = ct.residual_trend("cat_mag", x, rng.normal(0, 0.01, 200))
    assert not t.significant


def test_residual_trend_reports_rather_than_guesses_on_a_tiny_sample():
    t = ct.residual_trend("x", [1.0, 2.0], [0.0, 0.1])
    assert not t.significant and math.isnan(t.slope)


def test_plane_fit_recovers_a_planted_flat_field_tilt():
    """The flat-field question in its natural form: is the tie residual
    tilted across the detector?  Plant a tilt and demand the corner-to-
    corner swing back."""
    rng = np.random.default_rng(9)
    x = rng.uniform(0, 4000, 300)
    y = rng.uniform(0, 4000, 300)
    z = 0.01 + 2e-6 * x - 1e-6 * y + rng.normal(0, 0.003, 300)
    out = ct.plane_fit(x, y, z)
    assert out["cx"] == pytest.approx(2e-6, rel=0.25)
    assert out["cy"] == pytest.approx(-1e-6, rel=0.4)
    assert out["significant"]
    assert out["swing"] == pytest.approx(3e-6 * 4000, rel=0.35)


def test_plane_fit_calls_a_flat_detector_flat():
    rng = np.random.default_rng(10)
    out = ct.plane_fit(rng.uniform(0, 4000, 300), rng.uniform(0, 4000, 300),
                       rng.normal(0, 0.01, 300))
    assert not out["significant"]


def test_check_accuracy_applies_the_colour_term_to_check_stars():
    """Held-out stars are not the science target: their colours are known,
    so the FULL transformation may be applied to them.  That is the point —
    it tests the calibration the way a user would apply it."""
    colour = np.array([0.2, 0.8, 1.4])
    c_ref, k, zp = 0.8, 0.25, -3.0
    cat = np.array([15.0, 16.0, 17.0])
    nat = cat + k * (colour - c_ref)          # perfect data
    out = ct.check_accuracy(nat, cat, colour, zp, k, c_ref)
    assert out["n"] == 3
    assert out["rms"] == pytest.approx(0.0, abs=1e-12)


def test_check_accuracy_exposes_a_zero_point_bias():
    """The MEDIAN residual is the half that does not average away over
    stars, and is therefore the more dangerous half."""
    colour = np.array([0.5, 0.6, 0.7, 0.8])
    cat = np.array([15.0, 15.5, 16.0, 16.5])
    out = ct.check_accuracy(cat + 0.037, cat, colour, 0.0, 0.0, 0.6)
    assert out["median"] == pytest.approx(0.037, abs=1e-9)
    assert out["rms"] == pytest.approx(0.037, abs=1e-9)


def test_check_accuracy_on_an_empty_sample_is_not_zero_error():
    """No check stars must never read as perfect accuracy."""
    out = ct.check_accuracy([], [], [], 0.0, 0.0, 0.0)
    assert out["n"] == 0 and math.isnan(out["rms"])
    assert math.isnan(out["rms_clip"])


def test_check_accuracy_separates_the_typical_star_from_the_outliers():
    """Measured on this archive: 91% of held-out stars land inside 0.05 mag
    while two or three per block land 0.4-1.3 mag out.  The raw RMS is
    dominated by those few and would report the calibration as ten times
    worse than it is; the clipped RMS is the typical star.  BOTH must come
    back, because they answer different questions and neither alone is
    honest."""
    colour = np.full(40, 0.6)
    cat = np.linspace(15.0, 18.0, 40)
    nat = cat + 0.012 * np.sin(np.arange(40))     # ~9 mmag of real scatter
    nat[0] += 1.30                                # a variable
    nat[1] -= 0.75                                # a blend
    out = ct.check_accuracy(nat, cat, colour, 0.0, 0.0, 0.6)
    assert out["n"] == 40
    assert out["n_outlier"] == 2
    assert out["outlier_frac"] == pytest.approx(0.05)
    assert out["rms"] > 0.15                      # dominated by the two
    assert out["rms_clip"] < 0.02                 # the typical star
    assert out["rms_clip"] < out["rms"] / 5


def test_check_accuracy_clips_nothing_when_nothing_is_an_outlier():
    """The clip must be able to do nothing, or 'we clipped' would always be
    an excuse rather than a measurement."""
    rng = np.random.default_rng(21)
    cat = np.linspace(15.0, 18.0, 60)
    out = ct.check_accuracy(cat + rng.normal(0, 0.01, 60), cat,
                            np.full(60, 0.6), 0.0, 0.0, 0.6)
    assert out["n_outlier"] == 0
    assert out["rms_clip"] == pytest.approx(out["rms"])


# ==========================================================================
# 8 - verdicts, graded against the strategy's own numbers
# ==========================================================================

def test_block_verdict_is_graded_on_the_clipped_check_rms():
    """The fit clipped at CLIP_SIGMA; grading it with an UNCLIPPED check
    would measure the clipping rather than the calibration.  This test
    pins which of the two numbers the verdict consumes."""
    out = ct.check_accuracy(
        np.concatenate([np.full(30, 15.005), [16.9]]),
        np.full(31, 15.0), np.full(31, 0.6), 0.0, 0.0, 0.6)
    assert out["rms"] > ct.ACCURACY_GOAL_MAG        # raw fails the goal
    assert out["rms_clip"] <= ct.ACCURACY_GOAL_MAG  # clipped passes it
    assert ct.block_verdict(50, out["n"], out["rms_clip"], True) \
        == "TIED-STRETCH"
    assert ct.block_verdict(50, out["n"], out["rms"], True) \
        == "TIED-ABOVE-GOAL"


def test_block_verdict_grades_against_the_strategy_band():
    n = ct.MIN_TIE_STARS + 5
    c = ct.MIN_CHECK_STARS + 2
    assert ct.block_verdict(n, c, 0.008, True) == "TIED-STRETCH"
    assert ct.block_verdict(n, c, 0.015, True) == "TIED-GOAL"
    assert ct.block_verdict(n, c, 0.055, True) == "TIED-ABOVE-GOAL"


def test_block_verdict_refuses_to_tie_on_too_few_stars():
    """Below the minimum the honest answer is UNTIED.  A fabricated tie is
    worse than an honest gauge."""
    assert ct.block_verdict(ct.MIN_TIE_STARS - 1, 10, 0.005, True) == "UNTIED"
    assert ct.block_verdict(50, 10, 0.005, False) == "UNTIED"


def test_block_verdict_marks_an_unverifiable_tie_as_such():
    """A zero point with nothing independent behind it is the weakest
    passing grade, and it must not be able to masquerade as the best."""
    assert ct.block_verdict(50, 1, 0.001, True) == "TIED-UNVERIFIED"
    assert ct.block_verdict(50, 10, float("nan"), True) == "TIED-UNVERIFIED"


def test_goal_verdict_requires_coverage_and_accuracy_together():
    """The strategy's goal has two halves — every block tied, and tied to
    0.01-0.02 mag.  SUPPORTED needs both."""
    v, _ = ct.goal_verdict(["TIED-GOAL"] * 4, 4)
    assert v == "SUPPORTED"
    v, _ = ct.goal_verdict(["TIED-GOAL", "TIED-ABOVE-GOAL"], 2)
    assert v == "SUPPORTED-WITH-CAVEATS"
    v, _ = ct.goal_verdict(["TIED-GOAL", "UNTIED"], 2)
    assert v == "PARTIALLY SUPPORTED"
    v, _ = ct.goal_verdict(["UNTIED", "UNTIED"], 2)
    assert v == "NOT SUPPORTED"


def test_goal_verdict_states_the_deciding_number():
    """A verdict without its deciding number is an opinion."""
    _, deciding = ct.goal_verdict(["TIED-GOAL", "UNTIED", "TIED-STRETCH"], 3)
    assert "2 of 3" in deciding and "3" in deciding


# ==========================================================================
# 9 - archive access: HEALPix ranges instead of an IN-list
# ==========================================================================

def test_hpx_ranges_covers_every_source_it_was_given():
    """The one property that must never break: a range set that omits a
    source silently drops a tie star, and nothing downstream would notice."""
    rng = np.random.default_rng(1)
    ids = (4019898211471807872 + rng.integers(0, 5_000_000_000_000, 400))
    ranges = ct.hpx_ranges(ids)
    for s in ids:
        assert any(lo <= s <= hi for lo, hi in ranges), s


def test_hpx_ranges_coalesces_adjacent_healpix_pixels():
    """Sky-adjacent stars are id-adjacent, so a cone must collapse to a few
    runs — that collapse is the whole point of using ranges."""
    shift = 35 + 2 * (12 - ct.GAIA_HPX_LEVEL)
    base = 4019898211471807872 >> shift
    ids = [(base + k) << shift for k in range(30)]
    assert len(ct.hpx_ranges(ids)) == 1


def test_hpx_ranges_splits_when_the_sky_is_not_contiguous():
    shift = 35 + 2 * (12 - ct.GAIA_HPX_LEVEL)
    base = 4019898211471807872 >> shift
    ids = [base << shift, (base + 1) << shift, (base + 500) << shift]
    assert len(ct.hpx_ranges(ids)) == 2


def test_hpx_ranges_is_empty_for_no_sources():
    assert ct.hpx_ranges([]) == []


# ==========================================================================
# 10 - the external check: does the answer land where the literature says?
# ==========================================================================

def test_literature_check_accepts_a_target_inside_its_published_range():
    """The one validation that does not use the fitting catalogue, so the
    one that could catch a sign error or a mis-scaled gauge."""
    lo, hi, _ = ct.TARGET_V_RANGE["stlmi"]
    verdict, why = ct.literature_check("stlmi", 0.5 * (lo + hi))
    assert verdict == "inside"
    assert "AAVSO" in why


def test_literature_check_would_catch_an_inverted_zero_point():
    """A zero point applied with the wrong SIGN moves a magnitude by twice
    the zero point — several magnitudes here.  This is the test that would
    notice, and it must not be lenient enough to miss it."""
    verdict, _ = ct.literature_check("yzcnc", 10.3 - 2 * 2.86)
    assert verdict == "outside"


def test_literature_check_grades_a_near_miss_as_near_not_as_failure():
    """A polar in a deep low state sits below a range catalogued from its
    brighter epochs, and the bands measured here are not V.  A tenth or a
    magnitude outside is a caveat; ten magnitudes is a bug."""
    lo, _, _ = ct.TARGET_V_RANGE["euuma"]
    verdict, _ = ct.literature_check("euuma", lo - 0.5)
    assert verdict == "near"


def test_literature_check_is_honest_about_an_unknown_target():
    assert ct.literature_check("newthing", 15.0)[0] == "unknown"
    assert ct.literature_check("stlmi", float("nan"))[0] == "unknown"


def test_literature_ranges_are_broad_enough_to_be_a_bound_not_a_fit():
    """If any published range were narrow, someone would eventually be
    tempted to tune a zero point to it.  They are all wide."""
    for key, (lo, hi, src) in ct.TARGET_V_RANGE.items():
        assert hi - lo >= 2.0, key
        assert src, key


# ==========================================================================
# 12 - REGRESSIONS from adversarial review (2026-08-19)
#
# Each test below reproduces, in closed form, a defect that a reviewer found
# by EXECUTING the stage against the real products.  They are grouped here
# rather than scattered so that a future reader can see what the page was
# once wrong about, and what stops it being wrong about that again.
# ==========================================================================

# --- 12.1  the target colour must come from PAIRED EPOCHS ------------------
#
# The first implementation formed a variable target's colour by differencing
# the two filters' ENSEMBLE MEAN magnitudes.  For VV Pup era 76 the g
# campaign ran 370 days and the r campaign 55 of them, so the "colour" was
# the difference between two different accretion states: -1.73 published
# against +0.04 measured on shared epochs.  A 1.77 mag error that decided
# the block's stated verdict.

def test_paired_colour_recovers_a_variable_targets_true_colour():
    """A target that brightens by two magnitudes between campaigns still has
    a constant colour.  Paired epochs must find it; a difference of campaign
    means must not, and this test asserts BOTH so the distinction cannot be
    optimised away."""
    # Interleaved observations: blue at t, red 2 minutes later, all night.
    t = np.arange(0.0, 0.3, 0.01)
    state = np.where(t < 0.15, 0.0, -2.0)      # the object brightens midway
    blue = 17.0 + state
    red = blue - 0.40                          # colour is +0.40 throughout
    pc = ct.paired_colour(t, blue, t + 0.0014, red)
    assert pc.n_pairs == len(t)
    assert pc.colour == pytest.approx(0.40, abs=1e-9)
    assert pc.scatter == pytest.approx(0.0, abs=1e-9)
    # The recipe that was wrong: unequal sampling of the two states makes a
    # difference of means say something else entirely.
    naive = float(np.mean(blue[:20])) - float(np.mean(red[20:]))
    assert abs(naive - 0.40) > 1.0


def test_paired_colour_refuses_filters_that_never_observed_together():
    """VV Pup era 72's g and r share no night.  There is no colour to
    report, and reporting one anyway is how a campaign-mean difference gets
    mistaken for a measurement."""
    t_blue = np.linspace(0.0, 0.2, 20)
    t_red = np.linspace(5.0, 5.2, 20)          # five days later
    pc = ct.paired_colour(t_blue, np.full(20, 17.0),
                          t_red, np.full(20, 16.5))
    assert pc.n_pairs == 0
    assert not math.isfinite(pc.colour)
    assert "did not sample the same time" in pc.note


def test_paired_colour_needs_more_than_a_couple_of_orbital_phases():
    """One or two pairs of a strongly modulated object is one or two orbital
    phases, not a colour.  Below MIN_COLOUR_PAIRS the answer is unknown."""
    t_blue = np.array([0.0, 1.0, 2.0, 3.0])
    t_red = np.array([0.0, 1.0])               # only two epochs coincide
    pc = ct.paired_colour(t_blue, np.full(4, 17.0), t_red, np.full(2, 16.5))
    assert pc.n_pairs == 2 < ct.MIN_COLOUR_PAIRS
    assert not math.isfinite(pc.colour)


def test_paired_colour_reports_the_scatter_of_a_genuinely_variable_colour():
    """A polar's colour really does move.  The median is the answer the
    position test needs; the scatter must come back with it, or one number
    reads as more definite than the data is."""
    t = np.arange(0.0, 0.1, 0.005)
    swing = 0.2 * np.sin(2 * np.pi * t / 0.07)
    pc = ct.paired_colour(t, 17.0 + swing, t, np.full(len(t), 16.6))
    assert pc.colour == pytest.approx(0.4 + float(np.median(swing)), abs=1e-9)
    assert pc.scatter > 0.05


def test_paired_colour_tolerance_is_a_real_gate_not_a_formality():
    """Widening the tolerance must admit pairs and tightening it must
    exclude them; a tolerance that did nothing would let the campaign-mean
    error back in through the side door."""
    t_blue = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    t_red = t_blue + 0.02                       # 29 minutes apart
    kw = dict(t_blue=t_blue, m_blue=np.full(6, 17.0),
              t_red=t_red, m_red=np.full(6, 16.5))
    assert ct.paired_colour(**kw, tol_days=0.005).n_pairs == 0
    assert ct.paired_colour(**kw, tol_days=0.05).n_pairs == 6


def test_paired_colour_pairs_each_point_with_its_nearest_in_time():
    """Not the first, not the last: the NEAREST.  A pairing that walked the
    red series in order would silently drift out of step."""
    t_blue = np.array([1.0, 2.0, 3.0])
    t_red = np.array([3.001, 2.001, 1.001])     # deliberately unsorted
    m_red = np.array([13.0, 12.0, 11.0])        # the value at each red epoch
    pc = ct.paired_colour(t_blue, np.array([21.0, 22.0, 23.0]),
                          t_red, m_red, tol_days=0.01, min_pairs=1)
    assert pc.n_pairs == 3
    # Correct pairing gives 21-11, 22-12, 23-13 = 10 every time.
    assert pc.colour == pytest.approx(10.0)
    assert pc.scatter == pytest.approx(0.0, abs=1e-9)


# --- 12.2  a rigid astrometric offset is measurable, and must be measured --
#
# EU UMa era 78's reference solution placed every star 5.2 arcsec west of
# where the catalogue does, coherently, with 0.4 arcsec of scatter.  Not one
# of its 49 comparison stars matched inside the 1.2 arcsec photometric
# tolerance, and the block was reported untieable for "bad astrometry".

def _grid(n=40, ra0=180.0, dec0=30.0, span=0.2):
    """A small square field of catalogue sources."""
    rng = np.random.default_rng(7)
    return (ra0 + span * (rng.random(n) - 0.5) / math.cos(math.radians(dec0)),
            dec0 + span * (rng.random(n) - 0.5))


def test_rigid_offset_measures_and_removes_a_whole_field_translation():
    """The defect in closed form: displace every star by the same vector and
    the offset must come back, be applied, and restore every match."""
    cra, cdec = _grid()
    sra, sdec = ct.apply_offset(cra, cdec, +5.20, -0.05)   # our positions
    before = ct.match_by_sky(sra, sdec, cra, cdec)[0]
    assert (before >= 0).sum() == 0, "the defect must actually break matching"
    off = ct.rigid_offset(sra, sdec, cra, cdec)
    assert off.applied
    assert off.dra_arcsec == pytest.approx(-5.20, abs=0.02)
    assert off.ddec_arcsec == pytest.approx(+0.05, abs=0.02)
    assert off.scatter_arcsec < 0.1
    ura, udec = ct.apply_offset(sra, sdec, off.dra_arcsec, off.ddec_arcsec)
    assert (ct.match_by_sky(ura, udec, cra, cdec)[0] >= 0).sum() == len(cra)


def test_rigid_offset_refuses_an_incoherent_pairing():
    """A block whose positions are individually WRONG, not collectively
    displaced, has no single translation to remove.  Applying a median of
    noise would bend the astrometry onto the catalogue it is about to be
    tied to, which is the one thing this stage must never do."""
    cra, cdec = _grid(60)
    rng = np.random.default_rng(3)
    sra = cra + (rng.random(60) - 0.5) * 20 / 3600.0
    sdec = cdec + (rng.random(60) - 0.5) * 20 / 3600.0
    off = ct.rigid_offset(sra, sdec, cra, cdec)
    assert not off.applied
    assert "scatter" in off.reason


def test_rigid_offset_leaves_a_small_disagreement_alone():
    """Sub-arcsecond disagreement is what a good plate solution looks like.
    The match tolerance already absorbs it, and 'correcting' it would be
    fitting the cross-match's own noise."""
    cra, cdec = _grid()
    sra, sdec = ct.apply_offset(cra, cdec, 0.25, 0.10)
    off = ct.rigid_offset(sra, sdec, cra, cdec)
    assert not off.applied
    assert off.size_arcsec == pytest.approx(0.269, abs=0.02)
    assert "below the" in off.reason


def test_rigid_offset_refuses_to_measure_from_a_handful_of_stars():
    """A median of three pairings is not a measurement of anything."""
    cra, cdec = _grid(4)
    sra, sdec = ct.apply_offset(cra, cdec, 5.0, 0.0)
    off = ct.rigid_offset(sra, sdec, cra, cdec)
    assert not off.applied and off.n < ct.ASTROM_REFINE_MIN_STARS


def test_rigid_offset_search_radius_exceeds_the_photometric_tolerance():
    """The gate that made the defect invisible.  An offset bigger than the
    match tolerance is the only kind worth finding, so the search radius
    must be wider than the tolerance by a wide margin — or the test cannot
    see the thing it exists to see."""
    assert ct.ASTROM_LOOSE_TOL_ARCSEC >= 5 * ct.MATCH_TOL_ARCSEC
    assert ct.ASTROM_REFINE_MIN_ARCSEC < ct.MATCH_TOL_ARCSEC


def test_apply_offset_is_invertible_and_scales_ra_by_cos_dec():
    """A great-circle offset in RA is not a coordinate offset in RA.  At
    dec 60 the coordinate step is twice the angle."""
    ra, dec = np.array([100.0]), np.array([60.0])
    out_ra, out_dec = ct.apply_offset(ra, dec, 3600.0, 0.0)
    assert out_ra[0] - 100.0 == pytest.approx(2.0, rel=1e-3)
    back = ct.apply_offset(out_ra, out_dec, -3600.0, 0.0)
    assert back[0][0] == pytest.approx(100.0, abs=1e-9)


# --- 12.3  the blend veto belongs in the BAND BEING TIED -------------------
#
# The veto was computed from Gaia G while the fit ran against PS1 g/r/i.
# Where Gaia resolves a pair that PS1 does not, the G contrast is large and
# the r contrast is nil, so the pair passed the gate and entered the fit.

def test_blend_veto_read_in_the_wrong_band_passes_a_real_blend():
    """The exact shape of the VV Pup pairs review found: a neighbour 0.8
    arcsec away, 2.52 mag fainter in Gaia G and equal in r.  Read in G it
    passes; read in r it must not."""
    star = dict(saturated_frac=np.zeros(1), near_veto_frac=np.zeros(1),
                annulus_sep_arcsec=np.array([99.0]),
                annulus_dmag=np.array([np.inf]),
                second_sep_arcsec=np.array([99.0]),
                match_sep_arcsec=np.array([0.2]),
                cat_mag=np.array([16.238]), cat_colour=np.array([0.45]),
                cat_flag=np.zeros(1))
    passes_in_g, _ = ct.clean_mask(blend_sep_arcsec=np.array([0.80]),
                                   blend_dmag=np.array([+2.522]), **star)
    fails_in_r, census = ct.clean_mask(blend_sep_arcsec=np.array([0.80]),
                                       blend_dmag=np.array([-0.006]), **star)
    assert passes_in_g[0], "the old gate really did let this through"
    assert not fails_in_r[0], "the tie-band gate must catch it"
    assert census.n_blend_aperture == 1


def test_neighbour_metrics_answers_in_whatever_band_it_is_given():
    """The function is band-agnostic by construction — which is precisely
    why the CALLER has to pass the tie band's magnitudes and not a
    convenient common currency."""
    ra = np.array([180.0, 180.0002])            # 0.62 arcsec apart at dec 0
    dec = np.array([0.0, 0.0])
    g_like = np.array([16.3, 18.9])             # resolved pair, large dmag
    r_like = np.array([16.2, 16.2])             # unresolved in this band
    in_g = ct.neighbour_metrics(ra[:1], dec[:1], g_like[:1], ra, dec, g_like,
                                self_index=np.array([0]))
    in_r = ct.neighbour_metrics(ra[:1], dec[:1], r_like[:1], ra, dec, r_like,
                                self_index=np.array([0]))
    assert in_g["nn_sep"][0] == pytest.approx(in_r["nn_sep"][0])
    assert in_g["nn_dmag"][0] > ct.BLEND_DMAG        # waved through
    assert in_r["nn_dmag"][0] < ct.BLEND_DMAG        # caught


# --- 12.4  numbers that were typed by hand -------------------------------

def test_galactic_latitude_is_computed_not_remembered():
    """The page asserted VV Pup sits at galactic latitude +2.  It sits at
    +8.7.  Three reference positions, to a hundredth of a degree."""
    from macro_phot.report_cattie import _galactic
    for (ra, dec), (l_ref, b_ref) in (
            ((123.77833733311, -19.05493491834), (239.65, +8.71)),
            ((122.73604659258, +28.14257240475), (194.08, +28.81)),
            ((166.10689840347, +45.05387227348), (165.83, +62.15))):
        l, b = _galactic(ra, dec)
        assert l == pytest.approx(l_ref, abs=0.01)
        assert b == pytest.approx(b_ref, abs=0.01)


def test_saturation_fraction_is_not_a_module_constant():
    """The 7.2% saturated fraction was prose, not a measurement, and the
    database disagreed with it under one reading of the population.  If a
    literal ever reappears in ``cattie`` as a saturation FRACTION, this
    fails: the number belongs to a query, with its population named."""
    src = open(ct.__file__, encoding="utf-8").read()
    assert "7.2% of its" not in src


def test_the_only_hand_entered_numbers_are_the_declared_ones():
    """``cattie`` declares TARGET_V_RANGE as its only hand-entered values.
    That claim is worth a test because the page's footer repeats it."""
    for key, (lo, hi, src) in ct.TARGET_V_RANGE.items():
        assert "AAVSO" in src, key


# --- 12.5  the aperture veto must look at the WORST contaminant -----------
#
# Found while writing the regression test for 12.3: after the BAND was
# fixed, 55 blended stars were still inside the VV Pup fits, because the
# veto tested the NEAREST neighbour rather than the brightest one inside the
# aperture.  Aperture photometry sums every source in the circle.

def test_aperture_veto_sees_past_a_faint_nearest_neighbour():
    """A faint source at 1 arcsec and an EQUAL-brightness source at 3 arcsec
    are both inside a 6-arcsec aperture.  The nearest-neighbour test reports
    the harmless one and the star survives; the aperture test must report
    the one that actually corrupts the magnitude."""
    # dec 0, so 1 deg of RA is 1 deg on the sky.
    ra = np.array([10.0, 10.0 + 1.0 / 3600, 10.0 + 3.0 / 3600])
    dec = np.zeros(3)
    mag = np.array([16.0, 22.0, 16.05])     # near+faint, far+equal
    nm = ct.neighbour_metrics(ra[:1], dec[:1], mag[:1], ra, dec, mag,
                              self_index=np.array([0]))
    assert nm["nn_sep"][0] == pytest.approx(1.0, abs=0.01)
    assert nm["nn_dmag"][0] == pytest.approx(6.0, abs=0.01)   # looks clean
    assert nm["aper_sep"][0] == pytest.approx(3.0, abs=0.01)
    assert nm["aper_dmag"][0] == pytest.approx(0.05, abs=0.01)  # is not
    keep, census = ct.clean_mask(
        saturated_frac=np.zeros(1), near_veto_frac=np.zeros(1),
        blend_sep_arcsec=nm["aper_sep"][:1], blend_dmag=nm["aper_dmag"][:1],
        annulus_sep_arcsec=np.array([99.0]), annulus_dmag=np.array([np.inf]),
        second_sep_arcsec=np.array([99.0]), match_sep_arcsec=np.array([0.2]),
        cat_mag=np.array([16.0]), cat_colour=np.array([0.5]),
        cat_flag=np.zeros(1))
    assert not keep[0] and census.n_blend_aperture == 1


def test_aperture_veto_is_silent_when_the_aperture_is_empty():
    """A star alone in its aperture must not be charged for a neighbour
    that happens to be the nearest source in the whole field."""
    ra = np.array([10.0, 10.02])            # 72 arcsec away
    dec = np.zeros(2)
    mag = np.array([16.0, 15.0])
    nm = ct.neighbour_metrics(ra[:1], dec[:1], mag[:1], ra, dec, mag,
                              self_index=np.array([0]))
    assert not np.isfinite(nm["aper_sep"][0])
    assert not np.isfinite(nm["aper_dmag"][0])
    assert np.isfinite(nm["nn_sep"][0])     # the nearest source still exists


def test_aperture_radius_is_the_one_the_veto_uses():
    """If these two ever diverge, the metric answers about one circle and
    the gate decides about another."""
    ra = np.array([10.0, 10.0 + (ct.BLEND_APERTURE_ARCSEC - 0.1) / 3600])
    nm = ct.neighbour_metrics(ra[:1], np.zeros(1), np.array([16.0]),
                              ra, np.zeros(2), np.array([16.0, 16.0]),
                              self_index=np.array([0]))
    assert np.isfinite(nm["aper_dmag"][0])
    ra_far = np.array([10.0, 10.0 + (ct.BLEND_APERTURE_ARCSEC + 0.1) / 3600])
    nm2 = ct.neighbour_metrics(ra_far[:1], np.zeros(1), np.array([16.0]),
                               ra_far, np.zeros(2), np.array([16.0, 16.0]),
                               self_index=np.array([0]))
    assert not np.isfinite(nm2["aper_dmag"][0])
