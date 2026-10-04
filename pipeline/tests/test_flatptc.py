"""Unit tests for the S2 flat-pair photon transfer (rlmt_diagnostics.flatptc).

Same discipline as the sibling S2 tests: every case builds synthetic frames
with KNOWN gain and read noise and checks the pure functions recover them —
including the situations the estimator was written to survive (a level that
drifts between the two frames, a vignetted flat with pixel-response
non-uniformity, an additive contaminant that a bias pair does not have, an
on-camera average) and the ones in which it must decline to answer.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from rlmt_diagnostics import flatptc as fp

SHAPE = (640, 640)


def synth_flat(rng, level, gain=1.06, read_noise=3.9, bias=95.0,
               pattern=None, shape=SHAPE):
    """One synthetic flat: Poisson electrons / gain + Gaussian read noise.

    ``level`` is the mean signal in ADU above bias; ``pattern`` is the
    illumination x response map (1.0 = flat).  Returned in integer ADU,
    like a real frame.
    """
    pat = np.ones(shape) if pattern is None else pattern
    electrons = rng.poisson(level * gain * pat)
    return np.round(electrons / gain + rng.normal(0.0, read_noise, shape)
                    + bias)


def vignette(shape=SHAPE, depth=0.25, prnu=0.01, seed=7):
    """A radially vignetted illumination with pixel-to-pixel non-uniformity."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    cy, cx = shape[0] / 2.0, shape[1] / 2.0
    r2 = ((yy - cy) ** 2 + (xx - cx) ** 2) / (cy * cy + cx * cx)
    return (1.0 - depth * r2) * (1.0 + prnu * rng.standard_normal(shape))


def ptc_points(levels_ratios, gain=1.06, read_noise=3.9, bias=95.0,
               pattern=None, seed=0, group="g", tile=32):
    """Flat-pair points for a list of (level, ratio) pairs."""
    rng = np.random.default_rng(seed)
    pts = []
    for i, (level, ratio) in enumerate(levels_ratios):
        a = synth_flat(rng, level, gain, read_noise, bias, pattern)
        b = synth_flat(rng, level / ratio, gain, read_noise, bias, pattern)
        st = fp.tile_pair_stats(a, b, bias, tile=tile)
        for p in fp.pair_points(st):
            pts.append({**p, "pair_id": f"f{i}", "group": group})
    return pts


def zero_points(n=3, read_noise=3.9, bias=95.0, seed=100, tile=32):
    rng = np.random.default_rng(seed)
    pts = []
    for j in range(n):
        a = np.round(rng.normal(bias, read_noise, SHAPE))
        b = np.round(rng.normal(bias, read_noise, SHAPE))
        rn = fp.pair_read_noise(a, b, tile=tile)
        pts.append({"signal": 0.0, "var": rn["var"],
                    "var_err": rn["var_err"], "pair_id": f"z{j}",
                    "group": "zero"})
    return pts


class TestClipCorrection:
    def test_truncated_gaussian_factor_matches_simulation(self):
        """The clip-bias correction is a closed form; check it against a
        brute-force truncated Gaussian so a typo cannot hide in it."""
        rng = np.random.default_rng(1)
        x = rng.standard_normal(4_000_000)
        kept = x[np.abs(x) <= fp.CLIP_SIGMA]
        assert kept.var() == pytest.approx(fp.CLIP_VARIANCE_FACTOR, abs=1e-3)
        assert 0.998 < fp.CLIP_VARIANCE_FACTOR < 1.0

    def test_a_three_sigma_clip_would_cost_much_more(self):
        """Why the correction matters: it is 0.1% at 4 sigma but 2.7% at 3."""
        assert fp._truncated_gaussian_variance_factor(3.0) == pytest.approx(
            0.9733, abs=1e-3)


class TestEffectiveSignal:
    def test_reduces_to_the_textbook_form_at_ratio_one(self):
        s, v = fp.effective_signal_and_variance(np.array([200.0]),
                                                np.array([1.0]),
                                                np.array([1000.0]))
        assert s[0] == pytest.approx(1000.0)
        assert v[0] == pytest.approx(100.0)          # var(a-b)/2

    def test_ratio_propagation_is_exact(self):
        """var(a - r b) for known RN, K, S_b and r must map back onto the
        single-frame law V = RN^2 + S_eff/K with no approximation."""
        rn2, k, sb, r = 16.0, 1.3, 800.0, 0.75
        var_d = rn2 * (1 + r * r) + (r + r * r) * sb / k
        s_eff, v = fp.effective_signal_and_variance(var_d, r, sb)
        assert v == pytest.approx(rn2 + s_eff / k, rel=1e-12)


class TestTilePairStats:
    def test_recovers_variance_of_a_plain_pair(self):
        rng = np.random.default_rng(2)
        a = synth_flat(rng, 1500.0)
        b = synth_flat(rng, 1500.0)
        st = fp.tile_pair_stats(a, b, 95.0, tile=32)
        truth = 3.9 ** 2 + 1.0 / 12.0 + 1500.0 / 1.06
        assert st["n_used"] == st["n_tiles"] > 50
        assert np.mean(st["var"]) == pytest.approx(truth, rel=0.01)
        assert np.mean(st["signal"]) == pytest.approx(1500.0, rel=0.005)

    def test_fixed_pattern_cancels_despite_level_drift(self):
        """A vignetted, non-uniform flat whose level changed by 25% between
        the two frames: the plain difference would carry the pattern (PRNU
        1% of 3000 ADU = 30 ADU rms, against 53 ADU of noise); the scaled
        difference must not."""
        rng = np.random.default_rng(3)
        pat = vignette()
        a = synth_flat(rng, 3000.0, pattern=pat)
        b = synth_flat(rng, 2400.0, pattern=pat)
        st = fp.tile_pair_stats(a, b, 95.0, tile=32)
        resid = st["var"] - (3.9 ** 2 + 1 / 12.0 + st["signal"] / 1.06)
        assert abs(np.mean(resid)) / np.mean(st["var"]) < 0.01
        assert np.median(st["ratio"]) == pytest.approx(1.25, rel=0.005)

    def test_cosmic_rays_do_not_inflate_the_variance(self):
        rng = np.random.default_rng(4)
        a = synth_flat(rng, 800.0)
        b = synth_flat(rng, 800.0)
        hits = rng.integers(0, SHAPE[0], (300, 2))
        a[hits[:, 0], hits[:, 1]] += 5000.0
        st = fp.tile_pair_stats(a, b, 95.0, tile=32)
        truth = 3.9 ** 2 + 1 / 12.0 + 800.0 / 1.06
        assert np.mean(st["var"]) == pytest.approx(truth, rel=0.015)

    def test_masked_pixels_are_absent(self):
        rng = np.random.default_rng(5)
        a = synth_flat(rng, 800.0)
        b = synth_flat(rng, 800.0)
        a[:, ::7] += 400.0                    # a bad-column comb
        valid = np.ones(SHAPE, dtype=bool)
        valid[:, ::7] = False
        st = fp.tile_pair_stats(a, b, 95.0, valid=valid, tile=32)
        truth = 3.9 ** 2 + 1 / 12.0 + 800.0 / 1.06
        assert np.mean(st["var"]) == pytest.approx(truth, rel=0.015)

    def test_a_ratio_outside_the_window_is_refused(self):
        rng = np.random.default_rng(6)
        a = synth_flat(rng, 3000.0)
        b = synth_flat(rng, 900.0)            # ratio 3.3 — not the same scene
        assert fp.tile_pair_stats(a, b, 95.0, tile=32)["n_used"] == 0

    def test_frame_too_small_for_a_tile(self):
        st = fp.tile_pair_stats(np.ones((20, 20)), np.ones((20, 20)), 0.0)
        assert st["n_tiles"] == 0 and st["n_used"] == 0


class TestPairPoints:
    def test_errors_do_not_count_the_level_trend_as_scatter(self):
        """A vignetted pair spreads its tiles over a range of signal; the
        error bars must reflect noise about that trend, or chi2/dof falls
        far below 1 (a defect by standing rule 1)."""
        chi2nu = []
        for seed in range(6):
            pts = ptc_points([(400, 1.03), (900, 0.95), (1800, 1.08),
                              (2600, 1.0)], pattern=vignette(), seed=seed)
            fit = fp.fit_ptc_line(pts + zero_points(seed=200 + seed),
                                  n_boot=100)
            chi2nu.append(fit["chi2nu"])
        assert 0.6 < np.mean(chi2nu) < 1.6

    def test_too_few_tiles_gives_no_points(self):
        assert fp.pair_points({"signal": np.ones(5), "var": np.ones(5),
                               "ratio": np.ones(5)}) == []


class TestReadNoise:
    def test_bias_pair_read_noise(self):
        rng = np.random.default_rng(8)
        a = np.round(rng.normal(300.0, 2.3, SHAPE))
        b = np.round(rng.normal(300.0, 2.3, SHAPE))
        rn = fp.pair_read_noise(a, b, tile=32)
        truth = math.sqrt(2.3 ** 2 + 1.0 / 12.0)     # + quantisation
        assert rn["read_noise_adu"] == pytest.approx(truth, rel=0.01)
        assert rn["level"] == pytest.approx(300.0, abs=0.6)
        assert abs(rn["read_noise_adu"] - truth) < 5 * rn["read_noise_adu_err"]

    def test_hot_pixels_cancel(self):
        rng = np.random.default_rng(9)
        hot = np.zeros(SHAPE)
        hot[rng.integers(0, SHAPE[0], 2000),
            rng.integers(0, SHAPE[1], 2000)] = 3000.0
        a = np.round(rng.normal(300.0, 2.3, SHAPE)) + hot
        b = np.round(rng.normal(300.0, 2.3, SHAPE)) + hot
        rn = fp.pair_read_noise(a, b, tile=32)
        assert rn["read_noise_adu"] == pytest.approx(
            math.sqrt(2.3 ** 2 + 1 / 12.0), rel=0.01)


class TestFitLine:
    def test_recovers_gain_without_bias(self):
        """The estimator's own bias, measured: mean over seeds within a
        tenth of a percent of truth, and inside the quoted errors."""
        gains, errs = [], []
        for seed in range(8):
            pts = ptc_points([(300, 1.05), (600, 0.93), (1200, 1.1),
                              (2400, 1.0), (3000, 0.9)],
                             pattern=vignette(), seed=seed)
            fit = fp.fit_ptc_line(pts + zero_points(seed=300 + seed),
                                  n_boot=200)
            gains.append(fit["gain"])
            errs.append(fp.adopted_gain_error(fit))
        gains = np.array(gains)
        assert gains.mean() == pytest.approx(1.06, rel=2e-3)
        # Quoted errors must be of the size of the actual scatter.
        assert 0.3 < np.mean(errs) / gains.std(ddof=1) < 3.0

    def test_intercept_is_the_read_noise(self):
        pts = ptc_points([(300, 1.0), (1200, 1.0), (2400, 1.0)], seed=11)
        fit = fp.fit_ptc_line(pts + zero_points(seed=12), n_boot=100)
        assert fit["read_noise_adu"] == pytest.approx(
            math.sqrt(3.9 ** 2 + 1 / 12.0), rel=0.02)

    def test_quadratic_term_flags_uncancelled_pattern(self):
        """Feed the fit variances that carry an S^2 term (what residual
        fixed pattern or a non-linear response would do): the diagnostic
        must light up while a clean curve leaves it dark."""
        clean = [{"signal": s, "var": 15.0 + s / 1.06, "var_err": 0.5,
                  "pair_id": f"p{i}"}
                 for i, s in enumerate(np.linspace(200, 3000, 12))]
        dirty = [{**p, "var": p["var"] + 2e-5 * p["signal"] ** 2}
                 for p in clean]
        assert abs(fp.fit_ptc_line(clean, n_boot=60)["quad_z"] or 0) < 3
        assert fp.fit_ptc_line(dirty, n_boot=60)["quad_z"] > 5

    def test_chi2_is_reported_not_absorbed(self):
        """Errors ten times too small must show as chi2/dof ~ 100, and the
        formal gain error must NOT be rescaled to hide it."""
        rng = np.random.default_rng(13)
        pts = [{"signal": s, "var": 15.0 + s / 1.06 + rng.normal(0, 5.0),
                "var_err": 0.5, "pair_id": f"p{i}"}
               for i, s in enumerate(np.linspace(200, 3000, 30))]
        fit = fp.fit_ptc_line(pts, n_boot=100)
        assert fit["chi2nu"] > 30
        # ... and the bootstrap, which resamples the real scatter, is what
        # carries the honest error.
        assert fit["gain_err_boot"] > 3 * fit["gain_err_formal"]
        assert fp.adopted_gain_error(fit) == fit["gain_err_boot"]

    def test_refuses_degenerate_input(self):
        one_pair = [{"signal": s, "var": 15 + s, "var_err": 1.0,
                     "pair_id": "only"} for s in (100, 200, 300)]
        assert fp.fit_ptc_line(one_pair) is None
        assert fp.fit_ptc_line(one_pair[:2]) is None
        falling = [{"signal": s, "var": 5000 - s, "var_err": 1.0,
                    "pair_id": f"p{i}"} for i, s in enumerate((100, 500, 900))]
        assert fp.fit_ptc_line(falling) is None


class TestWithinGroups:
    """The estimator for sky flats: immune to an additive contaminant."""

    def _contaminated(self, seed=0):
        # Two "filters", each a brightening twilight sequence, each with
        # its own additive variance (faint stars in the difference image).
        pts = []
        for g, (extra, lv0) in enumerate([(3.0, 350.0), (6.0, 380.0)]):
            seq = ptc_points([(lv0 * 1.07 ** k, 0.985) for k in range(6)],
                             pattern=vignette(), seed=seed + 10 * g,
                             group=f"filter{g}")
            for p in seq:
                p["var"] += extra
                p["pair_id"] = f"g{g}-{p['pair_id']}"
            pts += seq
        return pts

    def test_pinned_fit_is_biased_by_an_additive_term(self):
        """The defect that motivated the estimator (found on the AC4040
        twilight flats): pinning the line at the bias variance tilts it."""
        pts = self._contaminated()
        pinned = fp.fit_ptc_line(pts + zero_points(seed=41), n_boot=100)
        within = fp.fit_ptc_within_groups(pts, n_boot=100)
        assert pinned["gain"] < 1.06 * 0.995        # biased LOW
        # ... and the misfit shows in chi2 relative to the model that
        # allows the additive term.
        assert pinned["chi2nu"] > within["chi2nu"]

    def test_within_group_fit_is_not(self):
        gains = [fp.fit_ptc_within_groups(self._contaminated(seed=s),
                                          n_boot=100)["gain"]
                 for s in range(6)]
        assert np.mean(gains) == pytest.approx(1.06, rel=0.01)

    def test_group_intercepts_recover_the_additive_terms(self):
        fit = fp.fit_ptc_within_groups(self._contaminated(seed=3), n_boot=60)
        rn2 = 3.9 ** 2 + 1 / 12.0
        got = {g["group"]: g["intercept"] - rn2 for g in fit["groups"]}
        assert got["filter0"] == pytest.approx(3.0, abs=2.5)
        assert got["filter1"] == pytest.approx(6.0, abs=2.5)
        ex = fp.additive_excess(fit["groups"], rn2, 0.05)
        assert ex["n_groups"] == 2 and ex["excess"] > 0

    def test_no_within_group_leverage_gives_no_answer(self):
        """Panel flats at one lamp setting: every pair at the same level,
        every point alone in its bin — nothing to regress on."""
        pts = [{"signal": 30000.0, "var": 28000.0, "var_err": 30.0,
                "pair_id": f"p{i}", "group": "panel"} for i in range(6)]
        assert fp.fit_ptc_within_groups(pts) is None

    def test_additive_excess_none_for_no_groups(self):
        assert fp.additive_excess([], 15.0) is None


class TestAutocorrelation:
    def test_independent_pixels_read_zero(self):
        rng = np.random.default_rng(20)
        a, b = synth_flat(rng, 1500.0), synth_flat(rng, 1500.0)
        ac = fp.diff_autocorrelation(a, b, 95.0, tile=64)
        assert abs(ac["rho_x"]) < 5 * ac["rho_x_err"] + 0.003
        assert abs(ac["rho_y"]) < 5 * ac["rho_y_err"] + 0.003

    def test_a_smoothed_frame_reads_positive(self):
        """Mix each pixel with its right-hand neighbour (what charge
        sharing does): the x correlation must appear, the y one must not."""
        rng = np.random.default_rng(21)
        def frame():
            f = synth_flat(rng, 1500.0)
            return 0.8 * f + 0.2 * np.roll(f, 1, axis=1)
        ac = fp.diff_autocorrelation(frame(), frame(), 95.0, tile=64)
        assert ac["rho_x"] > 0.15
        assert abs(ac["rho_y"]) < 0.02


class TestBinning:
    """Committee disagreement D2: is the 2x2 an average or a sum?"""

    @staticmethod
    def _binned_pair(rng, mode, native_gain=0.2467, level_native=6000.0):
        def frame():
            e = rng.poisson(level_native * native_gain, (640, 640))
            adu = e / native_gain + rng.normal(0, 4.0, (640, 640)) + 300.0
            blocks = adu.reshape(320, 2, 320, 2)
            return np.round(blocks.mean(axis=(1, 3)) if mode == "average"
                            else blocks.sum(axis=(1, 3)))
        return frame(), frame()

    @pytest.mark.parametrize("mode,expect", [("average", "average"),
                                             ("sum", "sum")])
    def test_ptc_gain_identifies_the_architecture(self, mode, expect):
        rng = np.random.default_rng(30)
        bias = 300.0 if mode == "average" else 1200.0
        a, b = self._binned_pair(rng, mode)
        st = fp.tile_pair_stats(a, b, bias, tile=32, edge_fraction=0.0)
        rn2_native = 16.0
        rn2 = rn2_native / 4.0 if mode == "average" else rn2_native * 4.0
        k_eff = np.mean(st["signal"]) / (np.mean(st["var"]) - rn2 - 1 / 12.0)
        v = fp.binning_verdict(k_eff, 0.02 * k_eff, 0.2467)
        assert v["verdict"] == expect
        assert v["ratio"] == pytest.approx(4.0 if mode == "average" else 1.0,
                                           rel=0.05)

    def test_a_gain_that_is_neither_says_so(self):
        assert fp.binning_verdict(0.6, 0.01, 0.2467)["verdict"] == "neither"


class TestFullScale:
    def test_full_scale_and_error(self):
        fs, err = fp.full_scale_electrons(3511.0, 94.0, 1.068, 0.003, 1.0)
        assert fs == pytest.approx((3511 - 94) * 1.068)
        assert err == pytest.approx(math.hypot(3417 * 0.003, 1.068 * 1.0))


class TestStarMask:
    def test_masks_sources_and_keeps_sky(self):
        rng = np.random.default_rng(40)
        img = rng.normal(500.0, 10.0, (256, 256))
        img[100:104, 60:64] += 400.0
        sky = fp.star_free_mask(img, tile=64)
        assert not sky[101, 61]
        assert not sky[99, 61]                  # grown by the dilation
        assert sky.mean() > 0.95

    def test_sky_pair_variance_survives_masked_stars(self):
        """Stars that move between two science frames must not reach the
        variance once the mask is applied."""
        rng = np.random.default_rng(41)
        truth = 3.9 ** 2 + 1 / 12.0 + 400.0 / 1.06
        a, b = synth_flat(rng, 400.0), synth_flat(rng, 400.0)
        for y, x in rng.integers(20, SHAPE[0] - 20, (120, 2)):
            a[y:y + 3, x:x + 3] += 900.0
            b[y:y + 3, x + 1:x + 4] += 900.0   # shifted by a pixel
        valid = fp.star_free_mask((a + b) / 2.0, tile=32)
        st = fp.tile_pair_stats(a, b, 95.0, valid=valid, tile=32)
        assert np.mean(st["var"]) == pytest.approx(truth, rel=0.03)
