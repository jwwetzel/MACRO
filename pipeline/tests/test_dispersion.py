"""Unit tests for the S2c per-frame dispersion classifier.

Every test builds a synthetic source list with KNOWN truth — a grism field,
an ordinary star field, a lone cosmic ray, a trailed exposure, an empty
frame — and checks the pure functions read it the way a human looking at the
picture would.  The cases that MUST NOT read as a grism get as much
attention as the ones that must: a classifier that says "spectrum" to a
satellite trail would quietly poison every downstream light curve.

Shape conventions used throughout, matching the extractor:
    a        semi-major second-moment sigma, px
    b        semi-minor, px
    theta    major-axis orientation in DEGREES over [-90, +90]
    flux     total counts (only its ORDER matters to the code)
    npix     connected pixels above threshold
"""

from __future__ import annotations

import numpy as np
import pytest

from rlmt_diagnostics import dispersion as dsp


# ---------------------------------------------------------------------------
# Synthetic frame builders — each returns the five parallel arrays
# ---------------------------------------------------------------------------
def _sources(specs):
    """Turn a list of (a, b, theta_deg, flux, npix) tuples into arrays."""
    arr = np.array(specs, dtype=float)
    return dict(a=arr[:, 0], b=arr[:, 1], theta_deg=arr[:, 2],
                flux=arr[:, 3], npix=arr[:, 4])


def grism_field():
    """A real grism frame's shape: two long traces on a common axis, the
    round zeroth order beside them, and faint round field stars."""
    return _sources([
        (905.0, 8.8, 1.3, 1.0e7, 163004),    # first-order trace of target
        (731.0, 5.0, 1.4, 2.1e6, 55578),     # trace of a second star
        (19.8, 18.0, -9.8, 8.0e5, 5385),     # undispersed zeroth order
        (8.9, 3.6, 1.3, 3.3e3, 235),         # faint star, trace unresolved
        (5.2, 2.0, 4.9, 1.4e3, 99),
        (3.3, 1.8, 9.2, 1.0e3, 69),
    ])


def star_field():
    """An ordinary direct image: every bright source at the seeing disc."""
    rng = np.random.default_rng(7)
    specs = []
    for i in range(12):
        specs.append((2.9, 2.7, float(rng.uniform(-90, 90)),
                      1.0e7 / (i + 1), 600 - 20 * i))
    return _sources(specs)


def trailed_field():
    """A wind-shaken exposure: EVERY source elongated, same drift angle, but
    only mildly — a/b of about three, and only a few pixels long."""
    specs = [(6.0, 2.0, 33.0 + 0.4 * i, 1.0e7 / (i + 1), 500 - 15 * i)
             for i in range(10)]
    return _sources(specs)


# ---------------------------------------------------------------------------
# axial_stats — the mod-180 circular statistics
# ---------------------------------------------------------------------------
class TestAxialStats:
    def test_plus_and_minus_ninety_are_the_same_axis(self):
        """The case that breaks naive averaging: +89 and -89 are 2 deg
        apart, not 178.  Their mean axis must be ~90, not ~0."""
        pa, scat = dsp.axial_stats([89.0, -89.0])
        assert pa == pytest.approx(90.0, abs=0.5)
        assert scat < 2.0

    def test_aligned_axes_have_zero_scatter(self):
        pa, scat = dsp.axial_stats([30.0, 30.0, 30.0])
        assert pa == pytest.approx(30.0, abs=1e-6)
        assert scat == pytest.approx(0.0, abs=1e-6)

    def test_random_axes_scatter_widely(self):
        rng = np.random.default_rng(3)
        _pa, scat = dsp.axial_stats(rng.uniform(-90, 90, size=400))
        assert scat > 30.0

    def test_perpendicular_pair_is_maximally_scattered(self):
        """0 and 90 deg cancel exactly under doubling; the mean direction is
        undefined and the code must report the ceiling, not crash."""
        _pa, scat = dsp.axial_stats([0.0, 90.0])
        assert scat == pytest.approx(90.0)

    def test_empty_returns_none(self):
        assert dsp.axial_stats([]) == (None, None)

    def test_pa_folded_into_zero_to_one_eighty(self):
        pa, _ = dsp.axial_stats([-45.0, -45.0])
        assert 0.0 <= pa < 180.0
        assert pa == pytest.approx(135.0, abs=1e-6)


# ---------------------------------------------------------------------------
# selection helpers
# ---------------------------------------------------------------------------
class TestSelection:
    def test_usable_mask_rejects_cosmics_and_slivers(self):
        mask = dsp.usable_mask(npix=[500, 5, 500], b=[2.0, 2.0, 0.1])
        assert list(mask) == [True, False, False]

    def test_select_bright_orders_by_flux(self):
        idx = dsp.select_bright([1.0, 50.0, 10.0], keep=2)
        assert list(idx) == [1, 2]

    def test_select_bright_handles_empty(self):
        assert dsp.select_bright([]).size == 0

    def test_trace_flags_need_both_ratio_and_length(self):
        # long+thin -> trace; short sliver -> not; long+fat galaxy -> not
        flags = dsp.trace_flags(a=[900.0, 3.0, 40.0], b=[9.0, 0.5, 20.0])
        assert list(flags) == [True, False, False]


# ---------------------------------------------------------------------------
# summarize_sources — the measurement record
# ---------------------------------------------------------------------------
class TestSummarize:
    def test_grism_field_measures_two_aligned_traces(self):
        s = dsp.summarize_sources(**grism_field())
        assert s.n_trace == 2
        assert s.trace_pa_scatter < 1.0          # the two traces agree
        assert s.trace_ab > 100.0
        # The MEDIAN over the bright set is low despite the obvious spectrum:
        # this is exactly why the verdict must not rest on the median.
        assert s.median_ab < 5.0

    def test_star_field_has_no_traces(self):
        s = dsp.summarize_sources(**star_field())
        assert s.n_trace == 0
        assert s.median_ab == pytest.approx(2.9 / 2.7, rel=1e-6)

    def test_cosmic_ray_is_cut_before_measurement(self):
        """A 12-pixel streak with a razor minor axis is a cosmic ray, and it
        must not survive the usability cut into any statistic."""
        s = dsp.summarize_sources(**_sources([
            (2.9, 2.7, 10.0, 1.0e6, 600),
            (30.0, 0.3, 45.0, 5.0e5, 12),        # the cosmic ray
        ]))
        assert s.n_sources == 1
        assert s.n_trace == 0

    def test_empty_frame_returns_empty_shape(self):
        s = dsp.summarize_sources(a=[], b=[], theta_deg=[], flux=[], npix=[])
        assert s.n_sources == 0 and s.median_ab is None and s.n_trace == 0

    def test_bright_set_is_capped(self):
        specs = [(2.9, 2.7, 0.0, float(100 - i), 300) for i in range(40)]
        s = dsp.summarize_sources(**_sources(specs), bright_n=10)
        assert s.n_sources == 40 and s.n_bright == 10


# ---------------------------------------------------------------------------
# classify_frame — the verdict, and the impostors it must reject
# ---------------------------------------------------------------------------
class TestClassify:
    def test_grism_field_reads_dispersed(self):
        v = dsp.classify_frame(dsp.summarize_sources(**grism_field()))
        assert v.verdict == dsp.VERDICT_DISPERSED
        assert "share an axis" in v.reason

    def test_star_field_reads_direct(self):
        v = dsp.classify_frame(dsp.summarize_sources(**star_field()))
        assert v.verdict == dsp.VERDICT_DIRECT

    def test_trailed_field_is_not_called_a_grism(self):
        """Wind smear elongates everything on a common axis — the shared-axis
        test alone would be fooled.  The absolute-length and ratio gates are
        what save it, so this frame must NOT read as dispersed."""
        v = dsp.classify_frame(dsp.summarize_sources(**trailed_field()))
        assert v.verdict != dsp.VERDICT_DISPERSED

    def test_satellite_trails_at_random_angles_are_indeterminate(self):
        """Two long streaks that do NOT share an axis are not a grating."""
        s = dsp.summarize_sources(**_sources([
            (400.0, 2.0, 10.0, 9.0e6, 40000),
            (380.0, 2.0, 75.0, 8.0e6, 38000),
            (2.9, 2.7, 3.0, 1.0e5, 300),
        ]))
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "common axis" in v.reason

    def test_isolated_bright_standard_reads_dispersed_on_one_trace(self):
        """Vega through the grism: one enormous trace, three faint round
        stubs.  There is no second trace to corroborate the axis, so the
        solo rule must carry it."""
        s = dsp.summarize_sources(**_sources([
            (557.0, 6.7, 7.5, 2.1e8, 247246),
            (13.3, 2.9, 89.0, 4.1e3, 275),
            (5.4, 3.6, 29.1, 2.4e3, 177),
            (4.1, 2.3, 88.6, 1.5e3, 98),
        ]))
        assert dsp.classify_frame(s).verdict == dsp.VERDICT_DISPERSED

    def test_satellite_across_a_rich_star_field_is_not_a_grism(self):
        """The real false positive that forced the sparsity gate: a 60-s
        luminance frame of M57 with 1,278 round stars and ONE 1,095-px
        streak across it.  A grating disperses EVERYTHING, so a rich field
        with exactly one smear is the one thing a grism cannot make."""
        specs = [(1095.0, 9.6, 74.1, 5.0e7, 40000)]          # the satellite
        specs += [(2.9, 2.7, float(10 * i % 90), 1.0e7 / (i + 1), 500)
                  for i in range(1, 1278)]                    # the star field
        s = dsp.summarize_sources(**_sources(specs))
        assert s.n_trace == 1 and s.n_sources > dsp.SOLO_MAX_SOURCES
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "rich field" in v.reason

    def test_one_modest_streak_is_indeterminate(self):
        """A single a/b = 6 object could be an edge-on galaxy or a satellite;
        one object is not a population, so we decline to rule."""
        s = dsp.summarize_sources(**_sources([
            (18.0, 3.0, 20.0, 1.0e6, 900),
            (2.9, 2.7, 40.0, 1.0e5, 300),
        ]))
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE

    def test_empty_frame_is_indeterminate_not_direct(self):
        """An empty frame is unmeasured, not photometry.  Calling it direct
        would silently certify a cloudy frame as good imaging."""
        s = dsp.summarize_sources(a=[], b=[], theta_deg=[], flux=[], npix=[])
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "no usable sources" in v.reason

    def test_badly_elongated_field_is_not_certified_as_direct(self):
        specs = [(9.0, 2.0, 30.0, 1.0e6 / (i + 1), 400) for i in range(8)]
        v = dsp.classify_frame(dsp.summarize_sources(**_sources(specs)))
        assert v.verdict == dsp.VERDICT_INDETERMINATE


# ---------------------------------------------------------------------------
# classify_strength — the second axis
# ---------------------------------------------------------------------------
class TestStrength:
    def test_high_and_low_split_at_the_calibrated_bounds(self):
        # The 2025-01-23 focus sweep's measured values, same star both ways.
        assert dsp.classify_strength(988.0, 4788) == "high"    # frac 0.206
        assert dsp.classify_strength(288.0, 4788) == "low"     # frac 0.060

    def test_overlap_band_returns_ambiguous_not_a_guess(self):
        mid = 0.5 * (dsp.STRENGTH_LOW_MAX_FRAC + dsp.STRENGTH_HIGH_MIN_FRAC)
        assert dsp.classify_strength(mid * 4788, 4788) == "ambiguous"

    def test_missing_measurement_is_not_a_class(self):
        assert dsp.classify_strength(None, 4788) == "n/a"
        assert dsp.classify_strength(500.0, 0) == "n/a"

    def test_width_normalisation_crosses_camera_formats(self):
        """The same physical dispersion on a 4,096-px iKon frame and a
        4,788-px GSENSE frame must land in the SAME class — that is the
        whole reason the statistic is a fraction and not a pixel count."""
        frac = 0.20
        assert (dsp.classify_strength(frac * 4096, 4096)
                == dsp.classify_strength(frac * 4788, 4788) == "high")

    def test_aspect_ratio_would_have_failed_here(self):
        """Documents WHY length replaced a/b: on the focus sweep the two
        grisms measured a/b 83 vs 61 — indistinguishable — while their
        length fractions differed by a factor of 3.4."""
        assert abs(988.0 / 4788) / abs(288.0 / 4788) > 3.0

    def test_direct_verdicts_carry_no_strength_class(self):
        v = dsp.classify_frame(dsp.summarize_sources(**star_field()))
        assert v.strength_class == "n/a"


# ---------------------------------------------------------------------------
# label bookkeeping (used only for scoring, never in the production path)
# ---------------------------------------------------------------------------
class TestLabels:
    def test_known_labels_map_to_expected_verdicts(self):
        assert dsp.expected_verdict("hrg") == dsp.VERDICT_DISPERSED
        assert dsp.expected_verdict("g") == dsp.VERDICT_DIRECT

    def test_the_disputed_labels_have_no_expectation(self):
        """Slot '6' and 'W' are the QUESTION.  Giving them an expected answer
        would smuggle the assumption we are trying to test into the scoring."""
        assert dsp.expected_verdict("6") is None
        assert dsp.expected_verdict("W") is None
        assert dsp.expected_verdict(None) is None

    def test_grism_units_map_to_strength_classes(self):
        assert dsp.expected_strength("hrg") == "high"
        assert dsp.expected_strength("HaGrism") == "high"
        assert dsp.expected_strength("lrg") == "low"
        assert dsp.expected_strength("OGGrism") == "low"
        assert dsp.expected_strength("g") is None


# ---------------------------------------------------------------------------
# the impure edge, exercised on synthetic pixels
# ---------------------------------------------------------------------------
class TestExtractSources:
    def _blank(self, shape=(300, 300), rng_seed=1):
        rng = np.random.default_rng(rng_seed)
        return rng.normal(100.0, 5.0, size=shape).astype(np.float32)

    def _add_gaussian(self, img, x, y, amp, sx, sy, angle_deg=0.0):
        """Paint one elliptical Gaussian at (x, y) rotated by angle."""
        h, w = img.shape
        yy, xx = np.mgrid[0:h, 0:w]
        t = np.radians(angle_deg)
        dx, dy = xx - x, yy - y
        u = dx * np.cos(t) + dy * np.sin(t)
        v = -dx * np.sin(t) + dy * np.cos(t)
        img += amp * np.exp(-0.5 * ((u / sx) ** 2 + (v / sy) ** 2))
        return img

    def test_synthetic_star_field_reads_direct(self):
        img = self._blank()
        for x, y in [(50, 60), (150, 80), (220, 200), (90, 240)]:
            self._add_gaussian(img, x, y, 3000.0, 3.0, 3.0)
        shape = dsp.extract_sources(img)
        assert shape.n_sources >= 4
        assert dsp.classify_frame(shape).verdict == dsp.VERDICT_DIRECT

    def test_synthetic_grism_field_reads_dispersed(self):
        """Two long thin Gaussians on a common axis plus round field stars —
        the picture a grism actually makes."""
        img = self._blank()
        self._add_gaussian(img, 150, 100, 4000.0, 60.0, 2.5, 20.0)
        self._add_gaussian(img, 150, 200, 3000.0, 55.0, 2.5, 20.0)
        for x, y in [(40, 40), (260, 60)]:
            self._add_gaussian(img, x, y, 1200.0, 3.0, 3.0)
        shape = dsp.extract_sources(img)
        assert shape.n_trace >= 2
        v = dsp.classify_frame(shape)
        assert v.verdict == dsp.VERDICT_DISPERSED
        # And the measured axis must be the one we painted.
        assert shape.trace_pa == pytest.approx(20.0, abs=6.0)

    def test_flat_frame_yields_no_measurement(self):
        """An all-constant array has no noise scale to threshold against;
        the code must return an empty record rather than divide by zero."""
        shape = dsp.extract_sources(np.full((100, 100), 42.0, dtype=np.float32))
        assert shape.n_sources == 0
        assert dsp.classify_frame(shape).verdict == dsp.VERDICT_INDETERMINATE

    def test_non_2d_input_is_refused_cleanly(self):
        shape = dsp.extract_sources(np.zeros((10,), dtype=np.float32))
        assert shape.n_sources == 0
