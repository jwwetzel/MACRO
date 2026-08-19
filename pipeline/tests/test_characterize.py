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
    factor, tau_used = ch.red_noise_factor(tau, adev, 480.0)
    assert factor == pytest.approx(1.0, rel=1e-9)
    assert tau_used == pytest.approx(480.0)


def test_red_noise_factor_exceeds_one_when_averaging_stalls():
    tau = np.array([60.0, 120.0, 240.0, 480.0])
    adev = np.array([0.01, 0.01, 0.01, 0.01])       # never averages down
    assert ch.red_noise_factor(tau, adev, 480.0)[0] == pytest.approx(
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


# ==========================================================================
# Additions made when the production build exposed the failure modes below
# ==========================================================================

def test_fit_noise_floor_is_not_dragged_by_the_faint_end():
    """The regression that this test locks down: an UNWEIGHTED fit in
    variance is dominated by the faintest stars, whose variance is 100x the
    bright end's, and returns their photon noise as a 'floor'.  Build a
    cloud whose true floor is 8 mmag but whose faint end reaches 300 mmag,
    and demand the fit still find 8."""
    pred = np.logspace(-3, np.log10(0.3), 200)
    rms = np.sqrt(pred ** 2 + 0.008 ** 2)
    floor, k, _ = ch.fit_noise_floor(np.arange(200.0), rms, pred)
    assert floor == pytest.approx(0.008, rel=0.02)
    assert k == pytest.approx(1.0, rel=0.02)


def test_noise_plateau_finds_the_flat_bottom():
    """Model-free floor: the smallest well-populated magnitude bin median."""
    mag = np.concatenate([np.full(20, 15.2), np.full(20, 18.2),
                          np.full(20, 21.2)])
    rms = np.concatenate([np.full(20, 0.012), np.full(20, 0.030),
                          np.full(20, 0.200)])
    floor, at_mag, n = ch.noise_plateau(mag, rms)
    assert floor == pytest.approx(0.012)
    assert at_mag == pytest.approx(15.25, abs=0.3)
    assert n == 20


def test_noise_plateau_ignores_sparse_bins():
    """A two-star bin cannot define the floor of a series."""
    mag = np.concatenate([np.full(2, 14.2), np.full(20, 18.2)])
    rms = np.concatenate([np.full(2, 0.001), np.full(20, 0.030)])
    floor, _, n = ch.noise_plateau(mag, rms, min_stars=5)
    assert floor == pytest.approx(0.030)
    assert n == 20


def test_noise_plateau_empty():
    assert np.isnan(ch.noise_plateau([], [])[0])


def test_degradation_threshold_run_length_rejects_a_lone_spike():
    """One noisy bin must not set an archive-wide cut (it did, in the first
    version of this build: a lone 1.4x bin threw away 69% of the frames)."""
    centers = np.arange(1.0, 9.0)
    med = np.array([1.0, 1.0, 1.6, 1.0, 1.0, 1.0, 1.0, 1.0])
    cnt = np.full(8, 50)
    thr, base = ch.degradation_threshold(centers, med, cnt, baseline=1.0)
    assert math.isinf(thr)
    assert base == pytest.approx(1.0)


def test_degradation_threshold_run_length_accepts_a_sustained_rise():
    centers = np.arange(1.0, 9.0)
    med = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.6, 1.7])
    thr, _ = ch.degradation_threshold(centers, med, np.full(8, 50), baseline=1.0)
    assert thr == pytest.approx(7.0)


def test_degradation_threshold_supplied_baseline_beats_the_minimum_bin():
    """With an explicit baseline the walk starts at bin 0, so a rise before
    the minimum bin is not skipped."""
    centers = np.array([1.0, 2.0, 3.0, 4.0])
    med = np.array([2.0, 2.0, 0.5, 0.5])
    thr, base = ch.degradation_threshold(centers, med, np.full(4, 50),
                                         baseline=1.0)
    assert base == pytest.approx(1.0)
    assert thr == pytest.approx(1.0)


def test_timing_precision_uses_a_supplied_noise_pool():
    """A pool of REAL residuals must be used in place of the Gaussian draw:
    feeding an all-zero pool has to give a far better answer than feeding
    noisy Gaussians of the same nominal sigma."""
    period = 113.9 / 1440.0
    t = np.arange(0, period, 120.0 / 86400.0)
    quiet = ch.timing_precision_mc(t, period, 1.0, 0.05, n_trials=60,
                                   ingress_phase=0.02,
                                   noise_pool=[np.zeros(t.size)] * 3)
    loud = ch.timing_precision_mc(t, period, 1.0, 0.05, n_trials=60,
                                  ingress_phase=0.02)
    assert quiet < loud


def test_timing_precision_penalises_a_mismatched_template():
    """Fitting with the wrong edge sharpness and depth must cost precision -
    the difference between the optimistic bound and the realistic case."""
    period = 113.9 / 1440.0
    t = np.arange(0, period, 220.0 / 86400.0)
    right = ch.timing_precision_mc(t, period, 1.0, 0.03, width_phase=0.45,
                                   ingress_phase=0.05, n_trials=150)
    wrong = ch.timing_precision_mc(t, period, 1.0, 0.03, width_phase=0.45,
                                   ingress_phase=0.05, n_trials=150,
                                   fit_ingress_phase=0.01,
                                   fit_depth_mag=1.2)
    assert wrong > right


def test_timing_precision_ignores_a_wrong_length_pool():
    """A residual vector that does not match the timestamp vector is not a
    noise realization for these timestamps and must be dropped, not
    broadcast."""
    period = 113.9 / 1440.0
    t = np.arange(0, period, 120.0 / 86400.0)
    out = ch.timing_precision_mc(t, period, 1.0, 0.02, n_trials=30,
                                 noise_pool=[np.zeros(t.size + 5)])
    assert np.isfinite(out) and out > 0


# ==========================================================================
# 6 - the adversarial-review regressions (2026-08-19)
#
# Each test below protects a specific published number that was wrong.  The
# comment names the number, so a future reader can see what was at stake
# rather than only what is asserted.
# ==========================================================================

class TestScoreModesAreDifferentQuestions:
    """A single-night contour scored on period tolerance is not a detection
    limit.  VV Pup's richest night cleared the threshold in 40 trials of 40
    at 300 mmag and was scored recovered in 25, and "never reached even at
    300 mmag" was published as a detection limit."""

    def _night(self, n=68, dt_s=200.0):
        # One real-ish night: 68 points at 200 s = 3.8 h, ~1.6 cycles of a
        # 100 min period, which is VV Pup's actual situation.
        return np.arange(n) * dt_s / 86400.0

    def test_the_acceptance_window_is_narrower_than_the_peak(self):
        """The arithmetic that makes 'period' scoring impossible on one
        night: 1% of f_orb against the 1/T frequency resolution."""
        t = self._night()
        period_d = 100.4 / 1440.0
        f_orb = 1.0 / period_d
        window_cd = ch.PERIOD_TOL_FRAC * f_orb          # acceptance half-width
        resolution_cd = 1.0 / (t.max() - t.min())       # peak half-width
        assert resolution_cd / window_cd > 20

    def _correlated_pool(self, t, seed=7):
        """Four random-walk residual series.

        White noise will NOT reproduce this defect: with white noise on a
        regular grid a high-amplitude peak's frequency is determined far
        better than 1/T and both score modes succeed.  The archive's
        check-star residuals are red (section 2 measures it), and it is red
        noise that makes the peak position wander across an acceptance
        window narrower than the peak.
        """
        rng = np.random.default_rng(seed)
        pool = [np.cumsum(rng.normal(0, 0.006, t.size)) for _ in range(4)]
        return [p - p.mean() for p in pool]

    def test_known_score_detects_what_period_score_misses(self):
        """Same injections, same noise: scoring at the known frequency
        recovers signals that period-scoring throws away.

        This is the published defect in one assertion.  At 300 mmag - the
        amplitude at which VV Pup's night was reported as "never reached" -
        the signal is detected in every trial; it is the FREQUENCY that
        cannot be pinned down from 1.6 cycles.
        """
        t = self._night()
        period_d = 100.4 / 1440.0
        freqs = np.arange(2.0, 40.0, 0.002)
        pool = self._correlated_pool(t)
        thr_p = ch.detection_threshold(t, pool, freqs, 100, 0.01,
                                       np.random.default_rng(1))
        thr_k = ch.detection_threshold(t, pool, freqs, 100, 0.01,
                                       np.random.default_rng(1),
                                       at_freq_cd=1.0 / period_d)
        got_p = ch.recovery_fraction(t, pool, freqs, period_d, 0.30, thr_p,
                                     40, rng=np.random.default_rng(2),
                                     score="period")
        got_k = ch.recovery_fraction(t, pool, freqs, period_d, 0.30, thr_k,
                                     40, rng=np.random.default_rng(2),
                                     score="known")
        assert got_k == pytest.approx(1.0)
        assert got_p < 0.95
        # And well below the contour level, so 'period' scoring would report
        # this amplitude as NOT recovered while the signal is plainly there.
        got_p_low = ch.recovery_fraction(t, pool, freqs, period_d, 0.10,
                                         thr_p, 40,
                                         rng=np.random.default_rng(2),
                                         score="period")
        got_k_low = ch.recovery_fraction(t, pool, freqs, period_d, 0.10,
                                         thr_k, 40,
                                         rng=np.random.default_rng(2),
                                         score="known")
        assert got_p_low < ch.RECOVERY_LEVEL < got_k_low

    def test_the_known_threshold_is_lower_than_the_max_threshold(self):
        """A known period buys back the look-elsewhere penalty.  Charging it
        anyway is what made the known-period limits 3-8x too pessimistic."""
        rng = np.random.default_rng(11)
        t = self._night()
        freqs = np.arange(2.0, 40.0, 0.002)
        pool = [rng.normal(0, 0.02, t.size) for _ in range(4)]
        thr_p = ch.detection_threshold(t, pool, freqs, 80, 0.01,
                                       np.random.default_rng(3))
        thr_k = ch.detection_threshold(t, pool, freqs, 80, 0.01,
                                       np.random.default_rng(3),
                                       at_freq_cd=14.3)
        assert thr_k < thr_p

    def test_an_unknown_score_mode_is_refused(self):
        with pytest.raises(ValueError):
            ch.recovery_fraction([1.0] * 10, [np.zeros(10)], [1.0], 0.1,
                                 0.01, 0.5, score="whatever")


class TestPeakClassification:
    def test_the_truth_the_alias_and_neither(self):
        f_orb = 1.0 / (113.9 / 1440.0)          # ST LMi, 12.64 c/d
        assert ch.classify_peak(f_orb, f_orb) == "true"
        assert ch.classify_peak(f_orb + 1.0, f_orb) == "alias"
        assert ch.classify_peak(f_orb - 2.0, f_orb) == "alias"
        assert ch.classify_peak(f_orb + 0.4, f_orb) == "other"

    def test_negative_aliases_fold_back(self):
        """A periodogram cannot tell +f from -f, so f - k c/d below zero is
        still a confusable frequency at its absolute value."""
        assert ch.classify_peak(1.5, 1.5, f_alias_cd=3.0) == "true"
        assert ch.classify_peak(1.5, 1.5, f_alias_cd=3.0, max_order=1) == "true"
        # f_true = 2, alias spacing 3 -> |2 - 3| = 1 is a real alias.
        assert ch.classify_peak(1.0, 2.0, f_alias_cd=3.0) == "alias"


class TestRedNoiseLabelling:
    """81 of 92 ladders never reached P_orb, and every one of them stored its
    factor in a column called red_factor_porb."""

    def test_the_tau_actually_used_is_returned(self):
        tau = np.array([60.0, 120.0, 240.0])
        adev = 0.01 / np.sqrt(tau / 60.0)
        factor, tau_used = ch.red_noise_factor(tau, adev, 6800.0)
        # The ladder stops at 240 s; the target is 6,800 s (a 113 min orbit).
        assert tau_used == pytest.approx(240.0)
        assert tau_used < 6800.0
        assert factor == pytest.approx(1.0, rel=1e-9)

    def test_a_ladder_entirely_above_the_target_answers_nan(self):
        factor, tau_used = ch.red_noise_factor([500.0, 1000.0],
                                               [0.01, 0.008], 100.0)
        assert math.isnan(factor) and math.isnan(tau_used)

    def test_the_white_null_is_not_minus_one_half(self):
        """THE comparison that was wrong: a 4-6 rung ladder's estimator does
        not have expectation -0.50, so -0.38 is not evidence against white
        noise on its own."""
        from macro_phot import errors as er
        null = ch.white_noise_allan_null(120, 200.0, 6800.0,
                                         er.allan_deviation, n_real=120,
                                         rng=np.random.default_rng(5))
        assert null["n"] > 100
        # The estimator's own median sits well below -0.50 ...
        assert null["slope_p50"] < -0.50
        # ... and its 5-95 band is wide enough to contain -0.38 comfortably.
        assert null["slope_p05"] < -0.38 < null["slope_p95"]

    def test_the_white_null_red_factor_is_not_one(self):
        """A short ladder returns a red-noise factor above 1 on pure white
        noise often enough that '1.5x worse than white' needs the null."""
        from macro_phot import errors as er
        null = ch.white_noise_allan_null(120, 200.0, 6800.0,
                                         er.allan_deviation, n_real=120,
                                         rng=np.random.default_rng(6))
        assert null["red_p95"] > 1.2


class TestPooledTailBinning:
    """Airmass: the bins from X = 2.45 to 2.65 all exceed the degradation
    factor, all three hold fewer than 15 frames, and the page reported that
    airmass never degrades the check stars over the observed range."""

    def _axis(self, n_tail=20):
        # 200 good frames below 2.4, then a degraded tail spread so thinly
        # across the top of the axis that no single bin can carry it - the
        # real airmass situation, where 36 frames sit in ten bins.
        x = np.concatenate([np.linspace(1.0, 2.3, 200),
                            np.linspace(2.45, 2.95, n_tail)])
        y = np.concatenate([np.full(200, 1.0), np.full(n_tail, 1.6)])
        return x, y, np.arange(1.0, 3.01, 0.1)

    def test_the_thin_tail_is_pooled_into_one_testable_bin(self):
        x, y, edges = self._axis()
        c, m, n = ch.binned_median_pooled_tail(x, y, edges, min_count=15)
        assert n[-1] == 20                   # the whole tail, in one bin
        assert m[-1] == pytest.approx(1.6)
        # And the pooled bin is placed at its LOWER EDGE, not its mean x:
        # it is a threshold candidate, and the threshold is where the
        # demonstrable-quality region ends.
        assert c[-1] < np.mean(x[x >= c[-1]])

    def test_plain_binning_discards_that_tail_entirely(self):
        x, y, edges = self._axis()
        _c, _m, n = ch.binned_median(x, y, edges)
        thin = n[(n > 0) & (n < 15)]
        assert thin.size >= 4               # every tail bin is under-populated

    def test_the_pooled_tail_can_set_a_threshold_alone(self):
        """A well-populated FINAL bin has no neighbour to confirm it, so
        requiring a run of two silently turns 'the axis ends here' into 'no
        effect'."""
        x, y, edges = self._axis()
        c, m, n = ch.binned_median_pooled_tail(x, y, edges, min_count=15)
        thr, _base = ch.degradation_threshold(c, m, n, baseline=1.0)
        assert math.isfinite(thr)
        assert thr == pytest.approx(2.4, abs=0.11)

    def test_a_tail_too_thin_even_pooled_stays_untestable(self):
        """Honesty in the other direction: if pooling still cannot reach the
        count rule, the answer is 'no threshold' - which the report must
        render as 'not testable above here', never as 'no effect'."""
        x, y, edges = self._axis(n_tail=5)
        c, m, n = ch.binned_median_pooled_tail(x, y, edges, min_count=15)
        thr, _base = ch.degradation_threshold(c, m, n, baseline=1.0)
        assert math.isinf(thr)

    def test_an_interior_spike_still_needs_a_run(self):
        """The run rule is not abandoned — only the terminal bin is exempt."""
        centers = np.arange(1.0, 6.0)
        med = np.array([1.0, 1.0, 1.9, 1.0, 1.0])
        thr, _ = ch.degradation_threshold(centers, med, np.full(5, 50),
                                          baseline=1.0)
        assert math.isinf(thr)


class TestColourPointError:
    """Q1 was graded on single-band precision.  A colour point costs at least
    sqrt(2) of that, and here the non-simultaneity term dominates."""

    def test_simultaneous_colour_is_root_two_worse(self):
        got = ch.colour_point_sigma(0.03, 0.03, rate_mag_per_s=0.0, dt_s=0.0)
        assert got == pytest.approx(0.03 * math.sqrt(2))

    def test_the_offset_term_can_dominate(self):
        """ST LMi 2025-02-27: median g-to-i offset 75 s, p90 sweep rate
        1.01 mmag/s -> 76 mmag from the offset alone, against ~30 mmag of
        photometry."""
        got = ch.colour_point_sigma(0.029, 0.030, rate_mag_per_s=0.00101,
                                    dt_s=75.0)
        assert got > 0.075
        assert got > 1.7 * ch.colour_point_sigma(0.029, 0.030, 0.0, 0.0)

    def test_nearest_time_offsets_finds_the_closer_neighbour(self):
        a = np.array([0.0, 1.0])
        b = np.array([-0.4, 0.9])
        got = ch.nearest_time_offsets(a, b)
        assert got[0] == pytest.approx(0.4 * 86400.0)
        assert got[1] == pytest.approx(0.1 * 86400.0)

    def test_rate_of_change_is_per_second(self):
        # 0.1 mag over 100 s.
        t = np.array([0.0, 100.0, 200.0]) / 86400.0
        m = np.array([0.0, 0.1, 0.2])
        got = ch.rate_of_change_mag_per_s(t, m)
        assert np.allclose(got, 1e-3)


class TestDutyCycleUncertainty:
    """Q4's deciding number was per-point precision against the state
    separation, which measures whether ONE night can be classified."""

    def test_thirteen_nights_cannot_beat_fourteen_points(self):
        # AN UMa: 13 usable nights.
        assert ch.duty_cycle_sigma(13) == pytest.approx(0.1387, abs=1e-3)

    def test_more_nights_help_as_root_n(self):
        assert (ch.duty_cycle_sigma(13) / ch.duty_cycle_sigma(52)
                == pytest.approx(2.0, rel=1e-6))

    def test_no_epochs_is_not_a_measurement(self):
        assert math.isnan(ch.duty_cycle_sigma(0))


class TestContourUncertainty:
    """A90 was printed to 0.1 mmag from 50 trials per cell."""

    def test_the_band_brackets_the_point_estimate(self):
        amps = np.array([0.002, 0.005, 0.01, 0.02, 0.05])
        fracs = np.array([0.0, 0.1, 0.5, 0.92, 1.0])
        a90 = ch.recovery_contour(amps, fracs)
        lo, hi = ch.contour_uncertainty(amps, fracs, n_trials=50,
                                        rng=np.random.default_rng(4))
        assert lo <= a90 <= hi
        # 50 trials cannot support 0.1 mmag: the band is a real fraction of
        # the value.
        assert (hi - lo) / a90 > 0.05

    def test_a_grid_that_never_reaches_the_level_answers_nan(self):
        lo, hi = ch.contour_uncertainty([0.01, 0.02], [0.1, 0.2], 50)
        assert math.isnan(lo) and math.isnan(hi)
