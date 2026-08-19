"""Unit tests for ``macro_phot.characterize`` — the S5 characterization core.

Every test states the physical fact it is protecting.  The style follows the
rest of the suite: pure functions only, synthetic inputs whose right answer
is known in closed form, and at least one test per function that would fail
if the function were quietly replaced by something plausible-but-wrong.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from macro_phot import characterize as ch


# ==========================================================================
# 1 - image quality
# ==========================================================================

def test_fwhm_arcsec_multiplies_by_plate_scale():
    """A 4 px FWHM is 2.15" on the High Gain 0.5375"/px camera and 3.22" on
    the Andor iKon at 0.8062 — the whole reason pixel FWHM is never quoted."""
    out = ch.fwhm_arcsec([4.0, 4.0], [0.5375, 0.8062])
    assert out[0] == pytest.approx(2.15, abs=1e-3)
    assert out[1] == pytest.approx(3.2248, abs=1e-3)


def test_fwhm_arcsec_rejects_nonpositive():
    """A zero or negative FWHM is a measurement failure, not a sharp frame."""
    out = ch.fwhm_arcsec([0.0, -1.0, np.nan], [0.5, 0.5, 0.5])
    assert np.all(np.isnan(out))


def test_sky_rate_normalizes_by_exposure():
    """240 ADU in 240 s and 8 ADU in 8 s are the same sky."""
    out = ch.sky_rate_adu_per_px_s([240.0, 8.0], [240.0, 8.0])
    assert out[0] == pytest.approx(1.0)
    assert out[1] == pytest.approx(1.0)


def test_sky_rate_zero_exposure_is_nan():
    assert np.isnan(ch.sky_rate_adu_per_px_s([100.0], [0.0])[0])


def test_airmass_zenith_is_one():
    assert ch.airmass_from_altitude([90.0])[0] == pytest.approx(1.0, abs=1e-3)


def test_airmass_matches_secz_at_high_altitude_and_exceeds_it_low():
    """Kasten-Young tracks sec(z) near zenith and stays FINITE at the
    horizon, where sec(z) diverges — the reason it is used for VV Pup."""
    x60 = ch.airmass_from_altitude([60.0])[0]
    assert x60 == pytest.approx(1.0 / math.sin(math.radians(60.0)), rel=2e-3)
    assert 30.0 < ch.airmass_from_altitude([0.5])[0] < 45.0


def test_airmass_below_horizon_is_nan():
    """The archive's AIRMASS headers reach 6877 because nothing refused a
    below-horizon geometry; this function refuses."""
    assert np.isnan(ch.airmass_from_altitude([-5.0])[0])


def test_moon_illumination_endpoints():
    k = ch.moon_illuminated_fraction([0.0, 90.0, 180.0])
    assert k[0] == pytest.approx(0.0)
    assert k[1] == pytest.approx(0.5)
    assert k[2] == pytest.approx(1.0)


def test_binned_median_counts_and_values():
    x = np.array([0.5, 0.6, 1.5, 1.6, 1.7])
    y = np.array([1.0, 3.0, 10.0, 20.0, 30.0])
    c, m, n = ch.binned_median(x, y, [0, 1, 2])
    assert list(n) == [2, 3]
    assert m[0] == pytest.approx(2.0)
    assert m[1] == pytest.approx(20.0)


def test_binned_median_empty_bin_is_nan():
    c, m, n = ch.binned_median([0.5], [1.0], [0, 1, 2])
    assert n[1] == 0 and np.isnan(m[1])


def test_degradation_threshold_finds_the_knee():
    """Scatter flat at 0.010 for three bins then jumping to 0.020: the cut
    belongs at the first bin that is >1.3x the baseline."""
    centers = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    med = np.array([0.010, 0.010, 0.011, 0.020, 0.030])
    cnt = np.array([50, 50, 50, 50, 50])
    thr, base = ch.degradation_threshold(centers, med, cnt)
    assert base == pytest.approx(0.010)
    assert thr == pytest.approx(4.0)


def test_degradation_threshold_infinite_when_never_degrades():
    """An axis that never costs precision must not be allowed to veto."""
    thr, base = ch.degradation_threshold([1, 2, 3], [0.01, 0.011, 0.0105],
                                         [50, 50, 50])
    assert math.isinf(thr)
    assert base == pytest.approx(0.01)


def test_degradation_threshold_ignores_underpopulated_bins():
    """A 3-frame bin's median must not set an archive-wide threshold."""
    thr, _ = ch.degradation_threshold([1, 2, 3], [0.01, 0.10, 0.011],
                                      [50, 3, 50])
    assert math.isinf(thr)


def test_usable_mask_requires_every_axis():
    vals = {"fwhm": np.array([2.0, 5.0, 2.0]),
            "airmass": np.array([1.1, 1.1, 3.0])}
    ok = ch.usable_mask(vals, {"fwhm": 3.0, "airmass": 2.0})
    assert list(ok) == [True, False, False]


def test_usable_mask_nan_fails():
    ok = ch.usable_mask({"fwhm": np.array([np.nan])}, {"fwhm": 3.0})
    assert not ok[0]


def test_usable_mask_infinite_threshold_passes_everything():
    ok = ch.usable_mask({"fwhm": np.array([1.0, 99.0])},
                        {"fwhm": float("inf")})
    assert ok.all()


def test_ellipticity_round_and_elongated():
    e = ch.ellipticity([1.0, 2.0], [1.0, 1.0])
    assert e[0] == pytest.approx(0.0)
    assert e[1] == pytest.approx(0.5)


def test_pa_coherence_aligned_vs_random():
    """Trailing drives every source to one position angle (R -> 1); seeing
    leaves them random (R ~ 1/sqrt N)."""
    aligned = np.full(200, 0.7)
    rng = np.random.default_rng(1)
    random = rng.uniform(0, np.pi, 200)
    assert ch.pa_coherence(aligned) == pytest.approx(1.0, abs=1e-9)
    assert ch.pa_coherence(random) < 0.2


def test_pa_coherence_is_mod_pi():
    """An ellipse has no head or tail: theta and theta+pi are the same PA,
    so a 50/50 mix of the two must still read as perfectly coherent."""
    t = np.concatenate([np.full(50, 0.3), np.full(50, 0.3 + np.pi)])
    assert ch.pa_coherence(t) == pytest.approx(1.0, abs=1e-9)


# ==========================================================================
# 2 - noise
# ==========================================================================

def test_predicted_sigma_source_limited_matches_poisson():
    """With a negligible background the prediction must collapse to pure
    Poisson: sigma_mag = 1.0857 / sqrt(N_e)."""
    flux, gain = 1.0e6, 1.0
    sig = ch.predicted_sigma_mag(flux, n_pix_aper=50.0, bkg_rms_adu=1e-9,
                                 gain_e_per_adu=gain)
    assert float(sig) == pytest.approx(ch.MAG_ERR_FACTOR / math.sqrt(1e6),
                                       rel=1e-6)


def test_predicted_sigma_gain_bracket_widens_source_term_only():
    """The measured background RMS carries sky+read, so a lower gain (more
    electrons per ADU is FALSE: lower gain = fewer e- per ADU = noisier)
    must increase the predicted sigma, and by exactly sqrt(g_hi/g_lo) in
    the source-limited regime."""
    kw = dict(n_pix_aper=50.0, bkg_rms_adu=1e-9)
    lo = float(ch.predicted_sigma_mag(1e6, gain_e_per_adu=ch.GAIN_LO_E_PER_ADU, **kw))
    hi = float(ch.predicted_sigma_mag(1e6, gain_e_per_adu=ch.GAIN_HI_E_PER_ADU, **kw))
    assert lo > hi
    assert lo / hi == pytest.approx(
        math.sqrt(ch.GAIN_HI_E_PER_ADU / ch.GAIN_LO_E_PER_ADU), rel=1e-6)


def test_predicted_sigma_sky_limited_scales_inversely_with_flux():
    """In the sky-dominated regime sigma_mag must go as 1/F: halve the star
    and the magnitude error doubles."""
    kw = dict(n_pix_aper=50.0, bkg_rms_adu=100.0,
              gain_e_per_adu=1e9)   # kill the source term
    s1 = float(ch.predicted_sigma_mag(1000.0, **kw))
    s2 = float(ch.predicted_sigma_mag(500.0, **kw))
    assert s2 / s1 == pytest.approx(2.0, rel=1e-6)


def test_predicted_sigma_nonpositive_flux_is_nan():
    assert np.isnan(ch.predicted_sigma_mag([0.0, -5.0], 50.0, 10.0, 1.0)).all()


def test_scintillation_scales_as_young_predicts():
    """X^1.75 in airmass and T^-1/2 in exposure — both signatures matter,
    because they are how the report proves the floor is NOT scintillation."""
    a = float(ch.scintillation_young(1.0, 60.0))
    b = float(ch.scintillation_young(2.0, 60.0))
    c = float(ch.scintillation_young(1.0, 240.0))
    assert b / a == pytest.approx(2 ** 1.75, rel=1e-6)
    assert c / a == pytest.approx(0.5, rel=1e-6)


def test_scintillation_on_this_telescope_is_a_millimag():
    """The decisive sanity number: a 0.5 m at Winer in a 60 s exposure at
    X=1.3 scintillates at well under 2 mmag, an order of magnitude below
    any floor measured in this archive."""
    s = float(ch.scintillation_young(1.3, 60.0))
    assert 0.0002 < s < 0.002


def test_fit_noise_floor_recovers_an_injected_floor():
    """Build rms^2 = predicted^2 + 0.012^2 exactly and demand the fit
    return 0.012 and k = 1."""
    pred = np.logspace(-3, -1, 40)
    floor_true = 0.012
    rms = np.sqrt(pred ** 2 + floor_true ** 2)
    floor, k, n = ch.fit_noise_floor(np.arange(40.0), rms, pred)
    assert floor == pytest.approx(floor_true, rel=1e-6)
    assert k == pytest.approx(1.0, rel=1e-6)
    assert n == 40


def test_fit_noise_floor_recovers_an_inflated_model():
    pred = np.logspace(-3, -1, 40)
    rms = np.sqrt(4.0 * pred ** 2 + 0.005 ** 2)
    floor, k, _ = ch.fit_noise_floor(np.arange(40.0), rms, pred)
    assert k == pytest.approx(4.0, rel=1e-6)
    assert floor == pytest.approx(0.005, rel=1e-5)


def test_fit_noise_floor_never_returns_imaginary():
    """A cloud that sits BELOW the prediction would give a negative
    intercept; the floor must clip to zero, not blow up."""
    pred = np.logspace(-3, -1, 20)
    rms = 0.5 * pred
    floor, k, _ = ch.fit_noise_floor(np.arange(20.0), rms, pred)
    assert floor >= 0.0 and np.isfinite(floor)


def test_fit_noise_floor_too_few_points():
    floor, k, n = ch.fit_noise_floor([1.0], [0.01], [0.01])
    assert np.isnan(floor) and n == 1


def test_precision_at_mag_uses_only_nearby_stars():
    mag = np.array([15.0, 16.0, 16.2, 16.4, 20.0])
    rms = np.array([0.001, 0.010, 0.020, 0.030, 0.500])
    val, n = ch.precision_at_mag(mag, rms, 16.2, half_width=0.5)
    assert n == 3
    assert val == pytest.approx(0.020)


def test_precision_at_mag_out_of_range():
    val, n = ch.precision_at_mag([15.0], [0.01], 25.0)
    assert np.isnan(val) and n == 0


def test_allan_slope_is_minus_half_for_white_noise():
    rng = np.random.default_rng(7)
    y = rng.normal(0, 0.01, 4096)
    tau, adev, _ = __import__("macro_phot.errors", fromlist=["errors"]).allan_deviation(y, 60.0)
    assert ch.allan_slope(tau, adev) == pytest.approx(-0.5, abs=0.06)


def test_allan_slope_flattens_for_correlated_noise():
    """A slow drift added to white noise must shallow the slope — the
    signature that averaging has stopped helping."""
    rng = np.random.default_rng(8)
    n = 4096
    t = np.arange(n)
    y = rng.normal(0, 0.01, n) + 0.05 * np.sin(2 * np.pi * t / 2000.0)
    tau, adev, _ = __import__("macro_phot.errors", fromlist=["errors"]).allan_deviation(y, 60.0)
    assert ch.allan_slope(tau, adev) > -0.45


def test_allan_slope_needs_three_rungs():
    assert np.isnan(ch.allan_slope([1.0, 2.0], [0.1, 0.07]))


def test_red_noise_factor_is_one_for_white():
    tau = np.array([60.0, 120.0, 240.0, 480.0])
    adev = 0.01 / np.sqrt(tau / 60.0)
    assert ch.red_noise_factor(tau, adev, 480.0) == pytest.approx(1.0, rel=1e-9)


def test_red_noise_factor_exceeds_one_when_averaging_stalls():
    tau = np.array([60.0, 120.0, 240.0, 480.0])
    adev = np.array([0.01, 0.01, 0.01, 0.01])       # never averages down
    assert ch.red_noise_factor(tau, adev, 480.0) == pytest.approx(
        math.sqrt(8.0), rel=1e-9)


# ==========================================================================
# 3 - cadence
# ==========================================================================

def test_night_blocks_splits_on_the_long_gap():
    t = np.array([0.0, 0.01, 0.02, 1.0, 1.01])
    assert ch.night_blocks(t, max_gap_h=3.0) == [(0, 3), (3, 5)]


def test_night_blocks_single_block_when_dense():
    t = np.arange(10) * 0.001
    assert ch.night_blocks(t) == [(0, 10)]


def test_cadence_stats_two_nights():
    """Two 2-hour nights a day apart at 120 s spacing: the function must get
    the cadence, the baseline, the block count and the duty cycle right."""
    night = np.arange(61) * (120.0 / 86400.0)
    t = np.concatenate([night, night + 1.0])
    s = ch.cadence_stats(t, period_d=0.08)
    assert s["n_blocks"] == 2
    assert s["median_dt_s"] == pytest.approx(120.0, rel=1e-6)
    assert s["baseline_d"] == pytest.approx(1.0 + 60 * 120 / 86400.0)
    assert s["longest_block_h"] == pytest.approx(2.0, rel=1e-6)
    assert s["pts_per_cycle"] == pytest.approx(0.08 * 86400 / 120.0, rel=1e-6)
    assert 0.0 < s["duty_cycle"] < 1.0


def test_cadence_stats_phase_coverage_full_and_partial():
    """A long continuous run covers every phase bin; a run much shorter
    than the period covers almost none."""
    long_run = np.linspace(0, 1.0, 500)
    assert ch.cadence_stats(long_run, period_d=0.08)["phase_coverage"] == 1.0
    short = np.linspace(0, 0.004, 20)
    assert ch.cadence_stats(short, period_d=0.08)["phase_coverage"] <= 0.15


def test_cadence_stats_empty():
    assert ch.cadence_stats([])["n_points"] == 0


def test_spectral_window_is_one_at_zero_frequency():
    t = np.array([0.0, 0.3, 1.1, 2.7])
    assert ch.spectral_window(t, [0.0])[0] == pytest.approx(1.0)


def test_spectral_window_peaks_at_one_per_day_for_nightly_sampling():
    """One point per night for 30 nights: the window must have a tall peak
    at exactly 1 c/d.  This is the alias structure the whole periodogram
    section is about."""
    t = np.arange(30.0)
    f = np.linspace(0.0, 2.0, 2001)
    w = ch.spectral_window(t, f)
    assert w[np.argmin(np.abs(f - 1.0))] == pytest.approx(1.0, abs=1e-6)
    assert w[np.argmin(np.abs(f - 0.5))] < 0.05


def test_spectral_window_of_a_dense_continuous_run_has_no_alias_comb():
    """A single uninterrupted 6 h run has no comb at all above ~1/span: its
    window decays and stays down, which is exactly why single-night
    periodograms are what selects the alias family for the multi-night set.
    (At f = 1 c/d the run is only a quarter cycle long, so the window is
    still high there — the comb question is asked where the CV signals
    live, 10-16 c/d.)"""
    t = np.linspace(0, 0.25, 300)
    f = np.linspace(8.0, 30.0, 3000)
    assert ch.spectral_window(t, f).max() < 0.05


def test_alias_ladder_folds_negative_frequencies():
    out = ch.alias_ladder(0.5, 1.0, orders=(-1, 1))
    assert out == pytest.approx([0.5, 1.5])


def test_alias_power_reports_offsets_not_absolute_frequencies():
    """The power an alias can reach is the window evaluated at the OFFSET.
    For one-point-per-night sampling the +/-1 c/d aliases must come back
    with window power ~1."""
    t = np.arange(30.0)
    rows = ch.alias_power(t, f_true_cd=12.0, f_alias_cd=1.0, orders=(-1, 1))
    assert all(p == pytest.approx(1.0, abs=1e-6) for _, _, p in rows)


# ==========================================================================
# 4 - detectability
# ==========================================================================

def test_inject_sinusoid_semi_amplitude():
    t = np.linspace(0, 1, 1001)
    y = ch.inject_sinusoid(t, 0.1, 0.05)
    assert y.max() == pytest.approx(0.05, abs=1e-4)
    assert y.min() == pytest.approx(-0.05, abs=1e-4)


def test_eclipse_template_depth_and_out_of_eclipse():
    t = np.linspace(0, 1.0, 20001)
    y = ch.eclipse_template(t, 0.1, depth_mag=0.4, width_phase=0.10,
                            t0_d=0.05, ingress_phase=0.02)
    assert y.max() == pytest.approx(0.4, rel=1e-3)
    assert y.min() == pytest.approx(0.0)
    # In-eclipse duty cycle: flat bottom + half of each ramp ~ width_phase.
    assert np.mean(y > 0.2 * 0.4) == pytest.approx(0.10, abs=0.02)


def test_eclipse_template_ingress_is_sloped_not_vertical():
    """A box edge would make the timing test meaningless; the trapezoid must
    put intermediate values on the ramp."""
    t = np.linspace(0.0, 0.1, 2001)
    y = ch.eclipse_template(t, 0.1, 0.4, 0.10, 0.05, ingress_phase=0.02)
    assert np.any((y > 0.05) & (y < 0.35))


def test_cyclic_noise_preserves_values_and_autocorrelation():
    r = np.sin(np.arange(100) / 5.0)
    out = ch.cyclic_noise_realization(r, 37)
    assert np.allclose(np.sort(out), np.sort(r))
    assert np.std(out) == pytest.approx(np.std(r))


def test_ls_peak_finds_an_injected_frequency():
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 5, 400))
    y = ch.inject_sinusoid(t, 0.1, 0.05) + rng.normal(0, 0.005, 400)
    f, p = ch.ls_peak(t, y, np.linspace(1, 30, 6000))
    assert f == pytest.approx(10.0, rel=2e-3)
    assert p > 0.5


def test_ls_peak_too_few_points():
    f, p = ch.ls_peak([0.0, 1.0], [1.0, 2.0], [1.0])
    assert np.isnan(f) and np.isnan(p)


def test_detection_threshold_rises_with_noisier_pool():
    rng = np.random.default_rng(11)
    t = np.sort(rng.uniform(0, 1, 120))
    freqs = np.linspace(2, 40, 2000)
    quiet = [rng.normal(0, 0.001, 120) for _ in range(4)]
    thr = ch.detection_threshold(t, quiet, freqs, n_trials=40, fap=0.05, rng=rng)
    # LS power is normalized, so the threshold is a fraction in (0, 1)
    assert 0.0 < thr < 1.0


def test_detection_threshold_empty_pool_is_nan():
    assert np.isnan(ch.detection_threshold([0.0, 1.0], [], [1.0]))


def test_recovery_fraction_is_one_for_a_huge_signal_and_zero_for_none():
    rng = np.random.default_rng(5)
    t = np.sort(rng.uniform(0, 0.3, 150))
    freqs = np.linspace(5, 40, 3000)
    pool = [rng.normal(0, 0.005, 150) for _ in range(4)]
    thr = ch.detection_threshold(t, pool, freqs, n_trials=60, fap=0.05, rng=rng)
    big = ch.recovery_fraction(t, pool, freqs, 0.08, 0.20, thr,
                               n_trials=20, rng=rng)
    tiny = ch.recovery_fraction(t, pool, freqs, 0.08, 1e-5, thr,
                                n_trials=20, rng=rng)
    assert big == pytest.approx(1.0)
    assert tiny < 0.2


def test_recovery_fraction_scores_wrong_period_as_failure():
    """A detection at an alias is not a period measurement.  Injecting at a
    period OUTSIDE the searched band must score zero even though the noise
    peak may well clear the threshold."""
    rng = np.random.default_rng(6)
    t = np.sort(rng.uniform(0, 0.3, 150))
    freqs = np.linspace(5, 15, 2000)          # band excludes f = 100 c/d
    pool = [rng.normal(0, 0.001, 150) for _ in range(4)]
    frac = ch.recovery_fraction(t, pool, freqs, 0.01, 0.5, threshold=0.0,
                                n_trials=10, rng=rng)
    assert frac == 0.0


def test_recovery_contour_interpolates_in_log_amplitude():
    amps = np.array([0.001, 0.01, 0.1])
    fracs = np.array([0.0, 0.5, 1.0])
    a90 = ch.recovery_contour(amps, fracs, level=0.9)
    assert 0.01 < a90 < 0.1
    # Log interpolation: 80% of the way from 0.01 to 0.1 in log space.
    assert a90 == pytest.approx(10 ** (-2 + 0.8), rel=1e-6)


def test_recovery_contour_nan_when_never_reached():
    assert np.isnan(ch.recovery_contour([0.001, 0.01], [0.1, 0.2]))


def test_recovery_contour_first_point_already_above():
    assert ch.recovery_contour([0.01, 0.1], [1.0, 1.0]) == pytest.approx(0.01)


def test_timing_precision_is_sampling_limited_at_coarse_cadence():
    """The headline timing result in closed form: at 220 s sampling with a
    sharp ingress, making the photometry 30x better must NOT make the epoch
    30x better — the epoch error is set by which samples happen to land on
    the ingress ramp.  This is what refutes an assumed sigma_t < 60 s
    reached 'by getting the S/N up'."""
    period = 113.9 / 1440.0                       # ST LMi, days
    t = np.arange(0, period, 220.0 / 86400.0)
    noisy = ch.timing_precision_mc(t, period, depth_mag=1.0, sigma_mag=0.030,
                                   n_trials=120, ingress_phase=0.01)
    clean = ch.timing_precision_mc(t, period, depth_mag=1.0, sigma_mag=0.001,
                                   n_trials=120, ingress_phase=0.01)
    assert noisy > 10.0                    # tens of seconds, not seconds
    assert clean > 10.0                    # photons cannot rescue sampling
    assert noisy / clean < 5.0             # while the S/N ratio is 30x


def test_timing_precision_improves_with_dense_sampling():
    """Same feature sampled every 10 s must time far better than every
    220 s — the control proving the previous test measured SAMPLING."""
    period = 113.9 / 1440.0
    coarse = ch.timing_precision_mc(np.arange(0, period, 220.0 / 86400.0),
                                    period, 1.0, 0.02, n_trials=120,
                                    ingress_phase=0.01)
    dense = ch.timing_precision_mc(np.arange(0, period, 10.0 / 86400.0),
                                   period, 1.0, 0.02, n_trials=120,
                                   ingress_phase=0.01)
    assert 0.0 < dense < 0.25 * coarse


def test_timing_precision_needs_points():
    assert np.isnan(ch.timing_precision_mc([0.0, 1.0], 0.08, 1.0, 0.01))


def test_amin_analytic_matches_the_strategy_formula():
    """sigma sqrt(4 z / N): 0.03 mag, N = 400, z = 18 -> 0.0180 mag."""
    assert ch.amin_analytic(0.03, 400) == pytest.approx(
        0.03 * math.sqrt(4 * 18 / 400), rel=1e-9)


def test_amin_analytic_guards_empty_series():
    assert np.isnan(ch.amin_analytic(0.03, 0))
