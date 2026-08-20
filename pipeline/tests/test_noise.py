"""Unit tests for the S2 empirical-noise logic (rlmt_diagnostics.noise).

Same discipline as the sibling S2 tests: every case builds its own tiny
synthetic input with KNOWN truth and checks the pure function recovers it —
including the cases that MUST NOT produce an answer (a thin bin, a single
frame pair, a level outside the measured span).
"""

from __future__ import annotations

import numpy as np
import pytest

from rlmt_diagnostics import noise


def synth_points(gain=1.0, read_noise_adu=4.0, levels=None, pairs=4,
                 n_pix=5000, seed=0):
    """(level, var, n_pix, pair_id) points from a KNOWN noise law.

    The truth is var = read_noise^2 + level/gain — the classic photon
    transfer line.  Each pair contributes its own noisy realisation, so the
    curve builder sees exactly what the archive gives it: several
    independent measurements per level.
    """
    rng = np.random.default_rng(seed)
    levels = levels if levels is not None else np.geomspace(20, 20000, 10)
    out = []
    for p in range(pairs):
        for lv in levels:
            truth = read_noise_adu ** 2 + lv / gain
            # A variance estimated from n_pix pixels scatters by
            # sqrt(2/n_pix) fractionally — emulate that, nothing more.
            v = truth * (1.0 + rng.normal(0, np.sqrt(2.0 / n_pix)))
            out.append((float(lv), float(v), n_pix, f"pair{p}"))
    return out


class TestEmpiricalNoiseCurve:
    def test_recovers_the_known_noise_law(self):
        pts = synth_points(gain=2.0, read_noise_adu=5.0)
        curve = noise.empirical_noise_curve(pts)
        assert curve, "a well-sampled mode must produce a curve"
        # Every published bin must sit on the truth line to a fraction of
        # a percent — the curve IS the measurement, so it must not distort.
        for c in curve:
            truth = 25.0 + c["level"] / 2.0
            assert c["var"] == pytest.approx(truth, rel=0.05)
        # Levels must come out sorted: downstream interpolation relies on it.
        assert [c["level"] for c in curve] == sorted(c["level"] for c in curve)

    def test_bins_record_their_independent_pair_count(self):
        curve = noise.empirical_noise_curve(synth_points(pairs=5))
        assert curve
        assert all(c["n_pairs"] >= noise.MIN_PAIRS_PER_BIN for c in curve)

    def test_single_pair_publishes_nothing(self):
        """One pair is an anecdote about one night, not a detector curve."""
        pts = synth_points(pairs=1, levels=np.geomspace(20, 20000, 40))
        assert noise.empirical_noise_curve(pts) == []

    def test_too_few_points_for_a_mode(self):
        pts = synth_points(pairs=2, levels=np.array([100.0, 200.0]))
        assert noise.empirical_noise_curve(pts) == []

    def test_ignores_nonfinite_and_nonpositive_points(self):
        pts = synth_points()
        pts += [(0.0, 5.0, 100, "junk"), (float("nan"), 5.0, 100, "junk"),
                (100.0, float("nan"), 100, "junk")]
        curve = noise.empirical_noise_curve(pts)
        assert curve
        assert all(np.isfinite(c["level"]) and np.isfinite(c["var"])
                   for c in curve)

    def test_sigma_is_the_root_of_the_variance(self):
        curve = noise.empirical_noise_curve(synth_points())
        assert all(c["sigma"] == pytest.approx(np.sqrt(c["var"]))
                   for c in curve)


class TestNoiseFloor:
    def test_reads_the_floor_at_the_bottom_of_the_curve(self):
        # Levels start well below the read-noise knee so the lowest bins
        # really are floor-dominated.
        pts = synth_points(gain=1.0, read_noise_adu=6.0,
                           levels=np.geomspace(1.0, 20000, 12))
        curve = noise.empirical_noise_curve(pts)
        fl = noise.noise_floor(curve)
        assert fl is not None
        # truth at the lowest level ~ 36 + 1 = 37 ADU^2 -> sigma ~ 6.08
        assert fl["floor_sigma_adu"] == pytest.approx(6.08, rel=0.10)
        assert fl["n_bins"] >= 1

    def test_empty_curve_has_no_floor(self):
        assert noise.noise_floor([]) is None


class TestCrossover:
    def test_finds_where_variance_doubles_the_floor(self):
        # gain 1, RN 6 -> floor var 36; 2x floor = 72 -> level ~ 36 ADU.
        pts = synth_points(gain=1.0, read_noise_adu=6.0,
                           levels=np.geomspace(1.0, 20000, 20))
        curve = noise.empirical_noise_curve(pts)
        fl = noise.noise_floor(curve)
        x = noise.crossover_level(curve, fl["floor_var_adu2"])
        assert x is not None
        assert x == pytest.approx(37.0, rel=0.35)

    def test_returns_none_when_never_witnessed(self):
        """A curve that never climbs off its floor must say so."""
        curve = [{"level": 10.0, "var": 36.0}, {"level": 12.0, "var": 36.2}]
        assert noise.crossover_level(curve, 36.0) is None

    def test_zero_floor_is_refused(self):
        assert noise.crossover_level([{"level": 1.0, "var": 1.0}], 0.0) is None


class TestPredictVariance:
    def test_interpolates_inside_the_measured_span(self):
        pts = synth_points(gain=2.0, read_noise_adu=5.0)
        curve = noise.empirical_noise_curve(pts)
        lv = 1000.0
        got = noise.predict_variance(curve, lv)
        assert got == pytest.approx(25.0 + lv / 2.0, rel=0.08)

    def test_refuses_to_extrapolate_below_the_span(self):
        curve = noise.empirical_noise_curve(synth_points())
        assert noise.predict_variance(curve, 0.001) is None

    def test_refuses_to_extrapolate_above_the_span(self):
        curve = noise.empirical_noise_curve(synth_points())
        assert noise.predict_variance(curve, 1e9) is None

    def test_degenerate_inputs(self):
        assert noise.predict_variance([], 100.0) is None
        assert noise.predict_variance([{"level": 1.0, "var": 1.0}], 1.0) is None
        curve = noise.empirical_noise_curve(synth_points())
        assert noise.predict_variance(curve, -5.0) is None


class TestCurveShapeIndex:
    def test_poisson_span_reads_one(self):
        """Far above the floor, variance tracks signal: log-log slope 1."""
        pts = synth_points(gain=1.0, read_noise_adu=0.01,
                           levels=np.geomspace(1000, 50000, 12))
        curve = noise.empirical_noise_curve(pts)
        assert noise.curve_shape_index(curve) == pytest.approx(1.0, abs=0.05)

    def test_floor_dominated_span_reads_near_zero(self):
        """Far below the knee, variance is flat: log-log slope ~ 0."""
        pts = synth_points(gain=1.0, read_noise_adu=100.0,
                           levels=np.geomspace(1.0, 50.0, 12))
        curve = noise.empirical_noise_curve(pts)
        assert noise.curve_shape_index(curve) == pytest.approx(0.0, abs=0.05)

    def test_too_short_a_curve(self):
        assert noise.curve_shape_index([{"level": 5.0, "var": 5.0}]) is None


class TestNarrowSpanIsRefused:
    """A curve measured over a few percent of level is not a curve.

    The archive's 5 MHz iKon science pairs all sit on a flat ~6,000 ADU
    sky: the measured levels span a ratio of 1.05, over which a straight
    line fit returned a "log-log slope" of 48.  That number describes the
    fit's own noise, not the detector, and must not be published.
    """

    def _flat_curve(self):
        # Ten bins spread over 5% of level — the 5 MHz situation.
        pts = []
        for p in range(4):
            for lv in np.linspace(6000, 6300, 10):
                pts.append((float(lv), 500.0 + (lv - 6000) * 0.3, 4000,
                            f"pair{p}"))
        return noise.empirical_noise_curve(pts)

    def test_slope_refused_over_a_narrow_span(self):
        curve = self._flat_curve()
        assert curve, "the curve itself is still published"
        assert noise.level_span_ratio(curve) < noise.MIN_SPAN_RATIO_FOR_SLOPE
        assert noise.curve_shape_index(curve) is None

    def test_slope_allowed_once_the_span_is_wide(self):
        curve = noise.empirical_noise_curve(
            synth_points(gain=1.0, read_noise_adu=0.01,
                         levels=np.geomspace(1000, 50000, 12)))
        assert noise.level_span_ratio(curve) >= noise.MIN_SPAN_RATIO_FOR_SLOPE
        assert noise.curve_shape_index(curve) == pytest.approx(1.0, abs=0.05)

    def test_floor_reports_when_it_swallowed_the_whole_curve(self):
        fl = noise.noise_floor(self._flat_curve())
        assert fl is not None
        assert fl["is_whole_curve"] is True

    def test_floor_knows_when_the_curve_rises_off_it(self):
        curve = noise.empirical_noise_curve(
            synth_points(gain=1.0, read_noise_adu=6.0,
                         levels=np.geomspace(1.0, 20000, 12)))
        fl = noise.noise_floor(curve)
        assert fl["is_whole_curve"] is False

    def test_span_ratio_of_an_empty_curve(self):
        assert noise.level_span_ratio([]) is None
