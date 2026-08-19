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


def mildly_trailed_field():
    """A wind-shaken exposure: EVERY source elongated on a common axis, but
    only mildly — a/b of about three, and only a few pixels long.  This one
    is genuinely below the trace gates, and is kept to pin that they work."""
    specs = [(6.0, 2.0, 33.0 + 0.4 * i, 1.0e7 / (i + 1), 500 - 15 * i)
             for i in range(10)]
    return _sources(specs)


def badly_trailed_field(pa: float = 2.2):
    """A REAL mount slip, built from measured numbers rather than wishful
    ones — the YZ Cnc 30-s i-band frames of 2024-05-02, which are labelled
    direct imaging and which the classifier called dispersed.

    This fixture exists to state an uncomfortable truth rather than hide it.
    An earlier version of this file synthesised its "trailed" frame at a/b
    3.0 and 6 px long, comfortably below TRACE_MIN_AB (5.0) and
    TRACE_MIN_A_PX (15.0), and then asserted the classifier did not call it a
    grism — a result guaranteed by arithmetic, which read as coverage while
    testing nothing.  The archive's actual trailed frames measure a/b near 90
    with traces 140 px long and a mutual PA scatter of 0.3 deg: they clear
    both trace gates and the parallelism test with room to spare, and they
    are geometrically indistinguishable from dispersion.

    What separates them is WHERE they point.  A mount slip runs along the
    drift, which has no reason to coincide with the grating axis; the default
    PA here is the measured 2.2 deg of the real frames, which does coincide,
    and so this fixture at its default IS still called dispersed.  That is
    the honest state of the classifier and the test below asserts it.
    """
    specs = [(141.7, 1.6, pa + 0.3 * ((i % 3) - 1), 1.0e7 / (i + 1),
              4000 - 100 * i) for i in range(10)]
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
        assert "share the grating axis" in v.reason

    def test_star_field_reads_direct(self):
        v = dsp.classify_frame(dsp.summarize_sources(**star_field()))
        assert v.verdict == dsp.VERDICT_DIRECT

    def test_mild_trailing_is_below_the_trace_gates(self):
        """Wind smear of a few pixels at a/b 3 does not reach TRACE_MIN_AB or
        TRACE_MIN_A_PX, so no trace is registered at all and the frame cannot
        be called dispersed.  This is the EASY half of trailing."""
        s = dsp.summarize_sources(**mildly_trailed_field())
        assert s.n_trace == 0
        assert dsp.classify_frame(s).verdict != dsp.VERDICT_DISPERSED

    def test_severe_trailing_on_the_grating_axis_is_still_a_false_positive(self):
        """The HARD half, asserted rather than wished away.

        A real mount slip (measured: a/b 90, 142 px, scatter 0.3 deg) clears
        every morphological gate this module has.  When the drift happens to
        run along the grating axis, nothing in the classifier can tell it
        from a spectrum, and it IS called dispersed.  This test pins the known
        residual error so that a future change which claims to fix trailing
        has to change this assertion deliberately."""
        s = dsp.summarize_sources(**badly_trailed_field(pa=2.2))
        assert s.n_trace >= 2 and s.trace_ab > dsp.TRACE_MIN_AB
        assert s.trace_a_px > dsp.TRACE_MIN_A_PX
        assert s.trace_pa_scatter < dsp.TRACE_MAX_PA_SCATTER_DEG
        assert dsp.classify_frame(s).verdict == dsp.VERDICT_DISPERSED

    def test_severe_trailing_off_the_grating_axis_is_rejected(self):
        """...but the same slip pointing anywhere else IS caught, because the
        grating's axis is fixed in detector coordinates and a mount's is not.
        This is the half of the trailing population the axis gate recovers."""
        s = dsp.summarize_sources(**badly_trailed_field(pa=47.0))
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "off the grating axis" in v.reason


class TestGratingAxis:
    """The gate that asks not just 'are the traces parallel?' but 'are they
    parallel to the GRATING?' — the archive's largest false-positive source
    was detector columns and bleed trails, which run at PA 90."""

    def test_offset_is_axial_not_directional(self):
        """PA 179 is 1 deg from an axis at 0, not 179 deg from it."""
        assert dsp.grating_axis_offset(179.0) == pytest.approx(1.0)
        assert dsp.grating_axis_offset(1.0) == pytest.approx(1.0)
        assert dsp.grating_axis_offset(-1.0) == pytest.approx(1.0)
        assert dsp.grating_axis_offset(180.0) == pytest.approx(0.0)

    def test_ninety_degrees_is_the_maximum_possible_offset(self):
        assert dsp.grating_axis_offset(90.0) == pytest.approx(90.0)

    def test_missing_angle_is_none_not_zero(self):
        """An unmeasured angle must not be silently read as 'aligned'."""
        assert dsp.grating_axis_offset(None) is None
        assert dsp.grating_axis_offset(float("nan")) is None

    def test_detector_columns_do_not_forge_a_grism(self):
        """The real mechanism: a field of column defects and bleed trails at
        PA 90.  They are PERFECTLY parallel — more parallel than any real
        spectrum — so a parallelism-only rule called them dispersed."""
        s = dsp.summarize_sources(**_sources([
            (523.9, 0.8, 90.0, 5.0e6, 30000),
            (410.0, 0.9, 90.0, 4.0e6, 25000),
            (380.0, 0.8, 89.9, 3.0e6, 22000),
            (2.9, 2.7, 12.0, 1.0e5, 300),
        ]))
        assert s.trace_pa_scatter < dsp.TRACE_MAX_PA_SCATTER_DEG
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "off the grating axis" in v.reason

    def test_solo_streak_off_axis_is_rejected(self):
        """The solo branch needs the axis test too — a single bleed trail in
        a sparse field otherwise reads as an isolated bright standard."""
        s = dsp.summarize_sources(**_sources([
            (1364.7, 4.9, 90.0, 5.0e7, 60000),
            (2.9, 2.7, 10.0, 1.0e5, 300),
            (2.8, 2.7, 40.0, 9.0e4, 290),
        ]))
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "off the grating axis" in v.reason

    def test_a_real_grism_frame_still_passes(self):
        """The gate must cost the true positives essentially nothing: the
        labelled grism populations sit at PA 0/180 by construction."""
        v = dsp.classify_frame(dsp.summarize_sources(**grism_field()))
        assert v.verdict == dsp.VERDICT_DISPERSED
        assert "grating axis" in v.reason


class TestDirectNeedsEvidence:
    """``direct`` is a certificate that the frame is safe for aperture
    photometry, not a shrug — so it needs a population behind it."""

    def test_one_round_source_cannot_certify_direct(self):
        """The false-negative channel: a faint spectrum that simply did not
        register leaves one round star, and the old rule voted 'clean'."""
        s = dsp.summarize_sources(**_sources([(2.9, 2.7, 12.0, 1.0e5, 300)]))
        v = dsp.classify_frame(s)
        assert v.verdict == dsp.VERDICT_INDETERMINATE
        assert "too little evidence" in v.reason

    def test_three_round_sources_are_enough(self):
        s = dsp.summarize_sources(**_sources([
            (2.9, 2.7, 12.0, 1.0e5, 300),
            (2.8, 2.7, 40.0, 9.0e4, 290),
            (3.0, 2.6, 71.0, 8.0e4, 280)]))
        assert dsp.classify_frame(s).verdict == dsp.VERDICT_DIRECT

    def test_the_dispersed_side_stays_reachable_on_one_source(self):
        """The asymmetry is deliberate: a bright standard alone in its field
        is a legitimate grism frame and must still be callable."""
        s = dsp.summarize_sources(**_sources([
            (905.0, 8.8, 1.3, 1.0e7, 163004)]))
        assert dsp.classify_frame(s).verdict == dsp.VERDICT_DISPERSED

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

    def test_defocused_frame_with_scattered_angles_is_rejected(self):
        """The control-sample false positive that forced the PA gate down
        from 20 deg to 5: a 2-second r-band exposure so defocused its blobs
        sat at 86, -73 and -43 deg.  Two of them cleared the trace gates,
        and their 15.1-deg axial scatter squeaked under the old bound.
        Real gratings hold their traces parallel to ~0.15 deg."""
        s = dsp.summarize_sources(**_sources([
            (38.9, 6.28, -72.7, 5.2e4, 662),
            (45.7, 4.88, -43.1, 2.7e4, 546),
            (15.2, 3.43, 86.5, 2.4e5, 914),
        ]))
        assert s.n_trace >= 2
        assert s.trace_pa_scatter > dsp.TRACE_MAX_PA_SCATTER_DEG
        assert dsp.classify_frame(s).verdict == dsp.VERDICT_INDETERMINATE

    def test_real_grism_alignment_clears_the_gate_with_room(self):
        """The measured median alignment of the labelled grism populations
        is 0.15 deg — two orders of magnitude inside the gate."""
        s = dsp.summarize_sources(**_sources([
            (905.0, 8.8, 1.30, 1.0e7, 163004),
            (731.0, 5.0, 1.45, 2.1e6, 55578),
        ]))
        assert s.trace_pa_scatter < 0.5
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
    def test_a_long_trace_names_the_high_dispersion_unit(self):
        # 988 px on a 4,788-px frame = 0.206, well past the 0.15 bound.
        assert dsp.classify_strength(988.0, 4788) == "high"

    def test_a_short_trace_is_ambiguous_and_never_called_low(self):
        """THE central honesty constraint of this axis.  A short trace is
        made by BOTH grisms — the H-alpha unit on a faint target looks like
        the broad unit on a bright one — so 'low' is never asserted.  A
        'low' call measured only 61% pure against a 47% base rate."""
        for a_px in (288.0, 150.0, 40.0, 1.0):
            assert dsp.classify_strength(a_px, 4788) == "ambiguous"
        assert "low" not in {dsp.classify_strength(a, 4788)
                             for a in (10.0, 100.0, 300.0, 600.0)}

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

    def test_threshold_is_the_calibrated_one(self):
        w = 4788
        eps = 1e-6
        assert dsp.classify_strength(
            (dsp.STRENGTH_HIGH_MIN_FRAC + eps) * w, w) == "high"
        assert dsp.classify_strength(
            (dsp.STRENGTH_HIGH_MIN_FRAC - 0.01) * w, w) == "ambiguous"

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
        """Two long thin Gaussians on the GRATING's axis plus round field
        stars — the picture a grism actually makes.

        The traces are painted at PA 0 rather than at an arbitrary angle,
        because that is where the archive's 18,312 labelled grism frames
        measure: the disperser is fixed in detector coordinates, so a
        synthetic 'grism field' at 20 deg was never a grism field."""
        img = self._blank()
        self._add_gaussian(img, 150, 100, 4000.0, 60.0, 2.5, 0.0)
        self._add_gaussian(img, 150, 200, 3000.0, 55.0, 2.5, 0.0)
        for x, y in [(40, 40), (260, 60)]:
            self._add_gaussian(img, x, y, 1200.0, 3.0, 3.0)
        shape = dsp.extract_sources(img)
        assert shape.n_trace >= 2
        v = dsp.classify_frame(shape)
        assert v.verdict == dsp.VERDICT_DISPERSED
        # And the measured axis must be the one we painted.
        assert dsp.grating_axis_offset(shape.trace_pa) < 6.0

    def test_synthetic_streaks_off_the_grating_axis_are_rejected(self):
        """The same picture rotated 20 deg is NOT a grism — it is drift or a
        pair of satellite trails, and end-to-end from pixels the classifier
        must say so."""
        img = self._blank()
        self._add_gaussian(img, 150, 100, 4000.0, 60.0, 2.5, 20.0)
        self._add_gaussian(img, 150, 200, 3000.0, 55.0, 2.5, 20.0)
        for x, y in [(40, 40), (260, 60)]:
            self._add_gaussian(img, x, y, 1200.0, 3.0, 3.0)
        shape = dsp.extract_sources(img)
        assert shape.n_trace >= 2
        assert dsp.classify_frame(shape).verdict == dsp.VERDICT_INDETERMINATE

    def test_flat_frame_yields_no_measurement(self):
        """An all-constant array has no noise scale to threshold against;
        the code must return an empty record rather than divide by zero."""
        shape = dsp.extract_sources(np.full((100, 100), 42.0, dtype=np.float32))
        assert shape.n_sources == 0
        assert dsp.classify_frame(shape).verdict == dsp.VERDICT_INDETERMINATE

    def test_non_2d_input_is_refused_cleanly(self):
        shape = dsp.extract_sources(np.zeros((10,), dtype=np.float32))
        assert shape.n_sources == 0
