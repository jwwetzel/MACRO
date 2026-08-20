"""Unit tests for the S2 pure logic (rlmt_diagnostics package).

Every test builds its own tiny synthetic case with KNOWN truth — synthetic
Poisson frames with a chosen gain, a fabricated histogram with a chosen
clip, a fabricated dark/flat pair — and checks the pure functions recover
it.  The cases that MUST NOT work are tested too: a smooth histogram must
not read as a clip, a hot pixel must not fake one, a ladder with too few
rungs must return None.
"""

from __future__ import annotations

import numpy as np
import pytest

from rlmt_diagnostics import ceiling as ceil
from rlmt_diagnostics import linearity as lin
from rlmt_diagnostics import ptc
from rlmt_diagnostics import reconstruct as rec


# ---------------------------------------------------------------------------
# ceiling
# ---------------------------------------------------------------------------
class TestFindClip:
    def _tail_hist(self, top=4000, floor=30):
        """A smooth falling bright tail with no clip."""
        h = np.zeros(5000, dtype=np.int64)
        adu = np.arange(top)
        h[:top] = (1e6 * np.exp(-adu / 300.0)).astype(np.int64) + floor
        return h

    def test_detects_pileup_spike(self):
        h = self._tail_hist(top=3526)
        h[3526] = 200_000                       # the pileup spike
        out = ceil.find_clip(h)
        assert out is not None
        assert out["clip_adu"] == 3526
        assert out["ratio"] >= ceil.PILEUP_MIN_RATIO

    def test_smooth_tail_is_not_a_clip(self):
        assert ceil.find_clip(self._tail_hist()) is None

    def test_lone_hot_pixel_is_not_a_clip(self):
        h = self._tail_hist(top=3000)
        h[4200] = 5                             # a hot pixel's few samples
        assert ceil.find_clip(h) is None

    def test_low_adu_spike_ignored(self):
        h = np.zeros(1000, dtype=np.int64)
        h[10] = 10_000_000                      # bias peak near zero
        assert ceil.find_clip(h) is None

    def test_empty_histogram(self):
        assert ceil.find_clip(np.zeros(10, dtype=np.int64)) is None

    def test_detects_cmos_pileup_mound(self):
        """The GSENSE4040 shape: per-pixel clip spread makes a ~60-bin-wide
        terminal mound (measured High Gain peak ~3,496), not one spike."""
        h = self._tail_hist(top=3450, floor=50)
        adu = np.arange(3440, 3560)
        h[3440:3560] += (700_000 * np.exp(-0.5 * ((adu - 3496) / 20.0) ** 2)
                         ).astype(np.int64)
        out = ceil.find_clip(h)
        assert out is not None
        assert abs(out["clip_adu"] - 3496) <= 2

    def test_truncated_sky_slope_rejected(self):
        """The 5 MHz iKon false positive: a steep sky peak whose dense range
        ends mid-slope must NOT read as a clip (mass continues above)."""
        h = np.zeros(20000, dtype=np.int64)
        adu = np.arange(20000)
        # Sky peak at 5,930 with broad wings; nothing special at the top.
        h += (2_000_000 * np.exp(-0.5 * ((adu - 5930) / 60.0) ** 2)
              ).astype(np.int64)
        h[6000:12000] += np.linspace(900, 100, 6000).astype(np.int64)
        assert ceil.find_clip(h) is None

    def test_clip_above_dead_tail(self):
        """A clip can terminate a tail that has fallen to zero counts."""
        h = np.zeros(4001, dtype=np.int64)
        h[100:2000] = 1000
        h[4000] = 50_000                        # spike far beyond the tail
        out = ceil.find_clip(h)
        assert out is not None and out["clip_adu"] == 4000


class TestVetoAndBits:
    def test_veto_matches_cv_adoption(self):
        # The measured High Gain clip must reproduce the CV team's 3,200.
        assert ceil.veto_threshold(3526) == 3200

    def test_veto_rounds_down(self):
        assert ceil.veto_threshold(1000) == 900   # 920 -> floor to 900

    def test_bit_depth_12(self):
        out = ceil.bit_depth_reading(3526)
        assert out["bits"] == 12
        assert out["adc_full_scale"] == 4096
        assert out["unused_codes"] == 4096 - 1 - 3526

    def test_bit_depth_16(self):
        assert ceil.bit_depth_reading(41143)["bits"] == 16

    def test_prior_comparison_exact(self):
        # Regression (adversarial review): the report must COMPUTE its
        # agreement-with-prior sentence, never hand-assert it.  Equal
        # integers — and only equal integers — earn "exactly".
        assert ceil.prior_comparison(3200, ceil.PRIOR_CV_VETO_ADU) \
            == "reproduces it exactly"

    def test_prior_comparison_discrepancy_is_stated(self):
        out = ceil.prior_comparison(3100, ceil.PRIOR_CV_VETO_ADU)
        assert "exactly" not in out
        assert "-100" in out                     # the signed discrepancy


class TestHistHelpers:
    def test_merge_hist_pads(self):
        out = ceil.merge_hist(np.array([1, 2]), np.array([1, 1, 5]))
        assert out.tolist() == [2, 3, 5]

    def test_frame_top_stats(self):
        img = np.zeros((10, 10), dtype=np.uint16)
        img[0, :3] = 3526
        st = ceil.frame_top_stats(img, veto_adu=3200)
        assert st["max_adu"] == 3526 and st["n_at_max"] == 3
        assert st["n_ge_veto"] == 3

    def test_frame_max_cluster_stackpro_shape(self):
        # 60 frames piling at ~56,000 + 60 unsaturated frames scattered.
        rng = np.random.default_rng(7)
        maxes = ([56000 + int(x) for x in rng.normal(0, 250, 60)]
                 + list(rng.uniform(900, 40000, 60)))
        out = ceil.frame_max_cluster(maxes)
        assert out is not None
        assert abs(out["clip_adu"] - 56000) < 600
        assert out["cluster_frac"] >= 0.3
        # Review addition: a cluster ceiling carries its members' spread as
        # the uncertainty (never NULL for a measured cluster).
        assert 0.0 <= out["mad_adu"] < 600

    def test_position_diversity_hot_pixel(self):
        # One stable hot pixel (with jitter) -> diversity ~ 1/n.
        pos = [(1423 + i % 3, 1522 - i % 2) for i in range(20)]
        assert ceil.position_diversity(pos) <= 0.1

    def test_position_diversity_true_ceiling(self):
        rng = np.random.default_rng(9)
        pos = [(int(y), int(x)) for y, x in rng.uniform(0, 4000, (20, 2))]
        assert ceil.position_diversity(pos) >= 0.8

    def test_position_diversity_empty(self):
        assert ceil.position_diversity([]) == 0.0

    def test_frame_max_cluster_rejects_scatter(self):
        rng = np.random.default_rng(8)
        assert ceil.frame_max_cluster(rng.uniform(1000, 60000, 200)) is None

    def test_mode_group_blank(self):
        assert ceil.mode_group("") == "(blank 2026)"
        assert ceil.mode_group(None) == "(blank 2026)"
        assert ceil.mode_group(" High Gain ") == "High Gain"


# ---------------------------------------------------------------------------
# ptc
# ---------------------------------------------------------------------------
class TestPtc:
    def _pair(self, gain=1.05, rn_adu=3.0, seed=1):
        """Two synthetic same-scene frames with known gain and read noise."""
        rng = np.random.default_rng(seed)
        # Scene: levels sweeping 50..2500 ADU (a star field's pixel range).
        level_adu = np.geomspace(50, 2500, 400_000)
        e = level_adu * gain                     # expected electrons
        a = rng.poisson(e) / gain + rng.normal(0, rn_adu, e.size)
        b = rng.poisson(e) / gain + rng.normal(0, rn_adu, e.size)
        return a, b

    def test_recovers_gain_and_read_noise(self):
        a, b = self._pair(gain=1.05, rn_adu=3.0)
        pts = ptc.pair_ptc_points(a, b, n_bins=10)
        assert len(pts) >= 5
        fit = ptc.fit_ptc([p["level"] for p in pts],
                          [p["var"] for p in pts],
                          [p["n_pix"] for p in pts])
        assert fit is not None
        assert fit["gain_e_per_adu"] == pytest.approx(1.05, rel=0.05)
        assert fit["read_noise_adu"] == pytest.approx(3.0, abs=1.0)

    def test_stackpro_variance_suppression(self):
        # Averaging N_sub=4 sub-frames quarters the variance -> apparent
        # gain 4x the true gain -> nsub_estimate reads 4.
        a, b = self._pair(gain=4 * 1.05, rn_adu=1.5, seed=2)
        pts = ptc.pair_ptc_points(a, b, n_bins=10)
        fit = ptc.fit_ptc([p["level"] for p in pts],
                          [p["var"] for p in pts])
        ns = ptc.nsub_estimate(fit["gain_e_per_adu"], 1.05)
        assert ns["nsub"] == 4 and ns["misfit"] < 0.3

    def test_robust_pair_variance_ignores_cosmic_rays(self):
        rng = np.random.default_rng(3)
        d = rng.normal(0, 2.0, 100_000)
        d[:20] = 5000.0                          # cosmic-ray hits
        var, n = ptc.robust_pair_variance(d)
        assert var == pytest.approx(2.0 ** 2 / 2, rel=0.05)
        assert n < d.size                        # the hits were clipped

    def test_fit_ptc_rejects_degenerate(self):
        assert ptc.fit_ptc([100.0, 100.0], [5.0, 5.0]) is None
        # Negative slope (variance falling with level) is not a PTC.
        assert ptc.fit_ptc([10, 100, 1000], [50, 30, 5]) is None

    def test_read_noise_from_dark_points(self):
        # Shortest darks: variance floor 15.4 ADU^2 at the offset level 94.
        pts = [(8.0, 94.0, 15.4), (8.0, 95.0, 15.8), (8.0, 170.0, 40.0),
               (64.0, 94.0, 39.0), (64.0, 1000.0, 170.0)]
        rn = ptc.read_noise_from_dark_points(pts)
        assert rn is not None
        assert rn["read_noise_adu"] == pytest.approx(np.sqrt(15.6), rel=0.02)
        assert rn["offset_adu"] == 94.0 and rn["exptime"] == 8.0

    def test_read_noise_empty(self):
        assert ptc.read_noise_from_dark_points([]) is None

    def test_read_noise_reports_offset_spread(self):
        # Review addition: the bias offset carries the floor bins' level
        # half-spread as its uncertainty ((95 - 94) / 2 here).
        pts = [(8.0, 94.0, 15.4), (8.0, 95.0, 15.8), (8.0, 170.0, 40.0)]
        rn = ptc.read_noise_from_dark_points(pts)
        assert rn["offset_adu_err"] == pytest.approx(0.5)

    def test_dark_shot_fraction_measures_floor_growth(self):
        # The archive shape the review flagged: 8 s floor 15.4 ADU^2,
        # 128 s floor 17.6 — the 8 s floor hides (17.6-15.4)*8/120 ADU^2
        # of dark shot variance (~1% of the floor, NOT zero).
        pts = [(8.0, 94.0, 15.4), (8.0, 95.0, 15.5),
               (128.0, 96.0, 17.6), (128.0, 97.0, 17.7)]
        out = ptc.dark_shot_fraction(pts)
        assert out is not None
        assert out["shot_var_adu2"] == pytest.approx(2.2 * 8 / 120, rel=0.05)
        assert out["frac_of_floor"] == pytest.approx(0.0095, abs=0.002)
        # The implied read-noise bias: dv / (2 sqrt(v)) ~ 0.019 ADU.
        assert out["rn_bias_adu"] == pytest.approx(0.0187, abs=0.003)

    def test_dark_shot_fraction_needs_two_exptimes(self):
        assert ptc.dark_shot_fraction([(8.0, 94.0, 15.4)]) is None
        # No growth (long floor below short) -> no measurable term.
        assert ptc.dark_shot_fraction([(8.0, 94.0, 15.4),
                                       (128.0, 94.0, 15.0)]) is None

    def test_stackpro_signature_reads_16(self):
        # The measured archive shape: offset x16, RN variance x16,
        # ceiling x16 -> N_sub consensus 16.
        hg = [(8.0, 94.0, 15.4), (8.0, 96.0, 15.5)]
        sp = [(32.0, 1504.0, 246.0), (32.0, 1520.0, 250.0)]
        sig = ptc.stackpro_signature(sp, hg, ceiling_sp=56062,
                                     ceiling_base=3496)
        assert sig is not None
        assert sig["nsub"] == 16
        assert sig["max_misfit"] < 1.0

    def test_amp_glow_metric_flags_hot_corner(self):
        img = np.full((1024, 1024), 100.0)
        img[-200:, -200:] += 40.0                # glowing corner
        m = ptc.amp_glow_metric(img, edge=128)
        assert m["hottest_corner_excess"] == pytest.approx(40.0, abs=1.0)
        assert m["center_med"] == pytest.approx(100.0)

    def test_pair_points_exclude_ceiling(self):
        a, b = self._pair()
        pts = ptc.pair_ptc_points(a, b, level_max=500.0)
        assert all(p["level"] <= 500.0 for p in pts)


# ---------------------------------------------------------------------------
# reconstruct
# ---------------------------------------------------------------------------
class TestReconstruct:
    def _pairs(self, n_pairs=20, n_pix=500, seed=4):
        """Synthetic raw/reduced pairs with a known dark and flat."""
        rng = np.random.default_rng(seed)
        D = 90.0 + 10.0 * rng.random(n_pix)       # per-pixel dark
        F = 0.9 + 0.2 * rng.random(n_pix)         # per-pixel flat
        sky = rng.uniform(50, 1500, size=(n_pairs, 1))   # scene sweep
        red = sky + rng.normal(0, 2.0, (n_pairs, n_pix))  # reduced (ped-free)
        raw = F * red + D
        return red, raw, D, F

    def test_recovers_dark_and_flat(self):
        red, raw, D, F = self._pairs()
        fit = rec.fit_pixel_lines(red, raw)
        assert np.nanmedian(np.abs(fit["F"] - F)) < 0.01
        assert np.nanmedian(np.abs(fit["D"] - D)) < 5.0

    def test_cosmic_ray_replacements_are_clipped(self):
        red, raw, D, F = self._pairs()
        red = red.copy()
        red[3, ::7] += 5000.0                     # astro-scrappy rewrites
        fit = rec.fit_pixel_lines(red, raw)
        assert np.nanmedian(np.abs(fit["F"] - F)) < 0.02
        # The poisoned pair was dropped for the affected pixels.
        assert fit["n_used"][::7].max() < red.shape[0]

    def test_too_few_pairs_gives_nan(self):
        red, raw, _, _ = self._pairs(n_pairs=rec.RECON_MIN_PAIRS - 1)
        fit = rec.fit_pixel_lines(red, raw)
        assert np.isnan(fit["F"]).all()

    def test_saturated_raw_masked(self):
        red, raw, D, F = self._pairs()
        fit = rec.fit_pixel_lines(red, raw, sat_adu=1000.0)
        # Pixels keep their fit (enough unsaturated pairs) and stay accurate.
        assert np.nanmedian(np.abs(fit["F"] - F)) < 0.02

    def test_sample_regions_geometry(self):
        regs = rec.sample_regions(4096, 4096)
        names = [r[0] for r in regs]
        assert names == ["center", "corner_tl", "corner_tr",
                         "corner_bl", "corner_br"]
        ys, xs = regs[0][1], regs[0][2]
        assert ys.stop - ys.start == rec.RECON_CENTER_SIZE
        for _nm, ys, xs in regs:
            assert 0 <= ys.start < ys.stop <= 4096
            assert 0 <= xs.start < xs.stop <= 4096

    def test_sample_regions_small_sensor(self):
        for _nm, ys, xs in rec.sample_regions(100, 80):
            assert 0 <= ys.start < ys.stop <= 100
            assert 0 <= xs.start < xs.stop <= 80

    def test_find_crop_offset(self):
        rng = np.random.default_rng(11)
        raw = rng.uniform(90, 4000, (600, 700))
        red = raw[18:18 + 582, 13:13 + 687] - 90.0   # crop + dark shift
        off = rec.find_crop_offset(raw, red)
        assert off is not None
        assert (off["dy"], off["dx"]) == (18, 13)

    def test_find_crop_offset_rejects_same_shape_and_big_diff(self):
        a = np.zeros((100, 100))
        assert rec.find_crop_offset(a, a) is None
        assert rec.find_crop_offset(np.zeros((400, 400)),
                                    np.zeros((100, 100))) is None

    def test_residual_vs_truth_removes_offset(self):
        truth = np.linspace(80, 120, 1000)
        out = rec.residual_vs_truth(truth + 1000.0, truth)  # pedestal shift
        assert out["offset"] == pytest.approx(1000.0)
        assert out["resid_rms"] == pytest.approx(0.0, abs=1e-9)

    def test_recon_verdict_identity_claims_only_sameness(self):
        # REGRESSION for the shipped era-79 blocker: an identity fit
        # (F ~ 1, D ~ 0, no pedestal) proves the two trees hold the same
        # pixels — it must NOT be read as "reduced = uncalibrated raw
        # copy" (era 79's shared content turned out to be CALIBRATED; the
        # old wording had the direction exactly backwards).
        v = rec.recon_verdict(1.00002, 0.0, None)
        assert "no relative calibration" in v
        low = v.lower()
        assert "raw copy" not in low
        assert "uncalibrated" not in low
        assert "byte" not in low                 # no byte-identity claim

    def test_recon_verdict_flat_and_truth(self):
        v = rec.recon_verdict(1.067, 1000.0, 2.04)
        assert "flat applied" in v and "2.0 ADU RMS" in v
        assert "no flat" in rec.recon_verdict(0.9998, 1000.0, None)
        assert rec.recon_verdict(None, None, None).startswith("unfittable")

    def test_flat_dark_correlation_flags_degeneracy(self):
        rng = np.random.default_rng(21)
        F = 1.0 + 0.09 * rng.standard_normal(4000)
        # The era-47 regime: D trades off against F (dD/dF ~ -33 ADU).
        D = 300.0 - 33.0 * (F - 1.0) + 0.5 * rng.standard_normal(4000)
        assert rec.flat_dark_correlation(F, D) < -0.9
        # Independent F and D (well-conditioned era): |corr| ~ 0.
        D_ind = 300.0 + 0.5 * rng.standard_normal(4000)
        assert abs(rec.flat_dark_correlation(F, D_ind)) < 0.1

    def test_flat_dark_correlation_handles_nans(self):
        F = np.full(50, np.nan)
        assert np.isnan(rec.flat_dark_correlation(F, F))

    def test_summarize_flags_unit_flat(self):
        F = np.full(100, 1.0)
        D = np.full(100, 95.0)
        s = rec.summarize_reconstruction(F, D, np.full(100, 2.0))
        assert s["flat_median"] == pytest.approx(1.0)
        assert s["dark_median"] == pytest.approx(95.0)
        assert s["fit_fraction"] == 1.0


# ---------------------------------------------------------------------------
# linearity
# ---------------------------------------------------------------------------
class TestLinearity:
    def test_parse_exptime_token(self):
        assert lin.parse_exptime_token("kaf_Vega_0p0001s_lrg_0.fts.fz") \
            == pytest.approx(0.0001)
        assert lin.parse_exptime_token("mjc_HD_20134_hrg_64s_x.fts") == 64.0
        assert lin.parse_exptime_token("no_token_here.fts") is None

    def test_effective_exptime_prefers_token_when_header_rounds_to_zero(self):
        assert lin.effective_exptime(0.0, "kaf_Vega_0p0001s_lrg_0.fts.fz") \
            == pytest.approx(0.0001)

    def test_effective_exptime_keeps_matching_header(self):
        # Header 63.9999 vs token 64: agreement -> header stands.
        assert lin.effective_exptime(63.9999, "x_64s_y.fts") \
            == pytest.approx(63.9999)

    def test_brightest_box_finds_the_star(self):
        img = np.full((512, 512), 100.0)
        img[300:310, 200:210] = 5000.0            # the star
        ph = lin.brightest_box_flux(img, box=96, stride=32)
        assert ph["y0"] <= 300 < ph["y0"] + 96
        assert ph["x0"] <= 200 < ph["x0"] + 96
        assert ph["flux"] == pytest.approx(100 * (5000 - 100), rel=0.01)
        assert ph["peak_adu"] == 5000.0

    def test_fit_ladder_flags_rolloff(self):
        t = [0.001, 0.01, 0.1, 1.0]
        f = [1000 * x for x in t[:3]] + [1000 * 1.0 * 0.7]   # top rung -30%
        fit = lin.fit_ladder(t, f)
        assert fit is not None
        assert fit["resid_pct"][-1] == pytest.approx(-30.0, abs=1.0)
        assert fit["max_abs_resid_pct"] == pytest.approx(30.0, abs=1.0)

    def test_fit_ladder_needs_three_rungs(self):
        assert lin.fit_ladder([1, 2], [10, 20]) is None

    def test_group_ladders(self):
        rows = [
            ("2024-05-20", "vega", "1MHz", 0.0001, 32),
            ("2024-05-20", "vega", "1MHz", 0.001, 32),
            ("2024-05-20", "vega", "1MHz", 0.01, 32),
            ("2024-05-20", "vega", "1MHz", 0.1, 32),
            ("2024-05-20", "other", "1MHz", 1.0, 5),     # one rung: no ladder
            ("2024-05-21", "x", "Fast", 1.0, 2),
            ("2024-05-21", "x", "Fast", 2.0, 1),         # rung too thin
            ("2024-05-21", "x", "Fast", 4.0, 2),
        ]
        out = lin.group_ladders(rows)
        assert len(out) == 1
        assert out[0]["target_key"] == "vega" and len(out[0]["rungs"]) == 4


# ---------------------------------------------------------------------------
# linearity: fair scheduling across readout modes
# ---------------------------------------------------------------------------
class TestFairLadderOrder:
    """The scheduler must try every mode before it repeats any mode.

    The archive is wildly uneven — hundreds of Mode0 candidates against one
    5 MHz candidate — so a globally-ranked list never reaches the sparse
    modes, and the campaign then reports "no archival linearity constraint"
    for modes it simply never looked at.
    """

    def _rows(self):
        # (mode, quality) pairs; lower quality value = better.
        return [("Mode0", i) for i in range(5)] + \
               [("High Gain", i) for i in range(3)] + \
               [("5MHz", 0)]

    def test_every_mode_appears_before_any_mode_repeats(self):
        out = lin.fair_ladder_order(self._rows(),
                                    mode_of=lambda r: r[0],
                                    quality=lambda r: r[1])
        first_three = [r[0] for r in out[:3]]
        assert sorted(first_three) == ["5MHz", "High Gain", "Mode0"]

    def test_best_candidate_of_each_mode_comes_first(self):
        out = lin.fair_ladder_order(self._rows(),
                                    mode_of=lambda r: r[0],
                                    quality=lambda r: r[1])
        assert all(r[1] == 0 for r in out[:3])

    def test_sparse_mode_is_reached_within_one_slot(self):
        """The failure this guards: 227 Mode0 candidates burying 5 MHz."""
        rows = [("Mode0", i) for i in range(227)] + [("5MHz", 0)]
        out = lin.fair_ladder_order(rows, mode_of=lambda r: r[0],
                                    quality=lambda r: r[1])
        assert ("5MHz", 0) in out[:2]

    def test_keeps_every_candidate_and_is_deterministic(self):
        rows = self._rows()
        a = lin.fair_ladder_order(rows, lambda r: r[0], lambda r: r[1])
        b = lin.fair_ladder_order(rows, lambda r: r[0], lambda r: r[1])
        assert a == b
        assert sorted(a) == sorted(rows)

    def test_empty_input(self):
        assert lin.fair_ladder_order([], lambda r: r[0], lambda r: r[1]) == []
