"""Unit tests for the hot-pixel flip test (pipeline/scripts/check_mech_flip.py).

The script answers one question about a mechanical boundary: did the SENSOR
turn, or only the saved ARRAY?  Everything that decides the answer is a pure
function over boolean maps, so these tests build synthetic sensors — a known
set of hot pixels on a noisy frame — and check that each transform is told
apart from the others and from chance.  No archive, no manifest.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import check_mech_flip as flip                                # noqa: E402

SHAPE = (120, 160)          # deliberately not square: mirrors must differ


def hot_set(n=60, seed=1):
    """``n`` distinct hot-pixel coordinates, away from the array centre of
    symmetry so that no transform maps the set onto itself."""
    rng = np.random.default_rng(seed)
    flat = rng.choice(SHAPE[0] * SHAPE[1], size=n, replace=False)
    return np.unravel_index(flat, SHAPE)


def frame(hot, seed, level=300.0, noise=5.0, excess=400.0, stars=0):
    """One synthetic exposure: flat level + noise + hot pixels (+ stars)."""
    rng = np.random.default_rng(seed)
    img = level + noise * rng.standard_normal(SHAPE)
    img[hot] += excess
    for _ in range(stars):
        # A star: a 5x5 Gaussian blob at a random place — wider than one
        # pixel, and never in the same place in two frames.
        y, x = rng.integers(4, SHAPE[0] - 4), rng.integers(4, SHAPE[1] - 4)
        yy, xx = np.mgrid[-2:3, -2:3]
        img[y - 2:y + 3, x - 2:x + 3] += 3000.0 * np.exp(
            -(yy ** 2 + xx ** 2) / 4.0)
    return img.astype(np.float32)


def maps_of(hot, seeds, **kw):
    """(strict, loose) hot-pixel maps from a stack of synthetic frames."""
    zs = [flip.residual_in_sigma(frame(hot, s, **kw)) for s in seeds]
    return (flip.persistent_mask([z > flip.CANDIDATE_SIGMA for z in zs],
                                 flip.PERSISTENCE),
            flip.persistent_mask([z > flip.LOOSE_SIGMA for z in zs],
                                 flip.LOOSE_PERSISTENCE))


class TestCandidates:
    def test_hot_pixels_are_found_and_noise_is_not(self):
        hot = hot_set()
        mask = flip.candidate_mask(frame(hot, seed=3))
        assert mask[hot].all()
        # 8 sigma on 19,200 Gaussian pixels: essentially nothing else.
        assert mask.sum() - len(hot[0]) <= 1

    def test_a_star_is_mostly_not_a_hot_pixel_and_never_persists(self):
        hot = hot_set()
        masks = [flip.candidate_mask(frame(hot, seed=s, stars=6))
                 for s in range(8)]
        strict = flip.persistent_mask(masks)
        # Whatever a star core leaves in one frame is gone from the stack.
        assert strict[hot].all()
        assert strict.sum() == len(hot[0])

    def test_a_constant_image_has_no_outliers(self):
        assert flip.residual_in_sigma(np.full(SHAPE, 303.0)) is None
        assert not flip.candidate_mask(np.full(SHAPE, 303.0)).any()

    def test_a_master_with_zero_mad_still_yields_its_hot_pixels(self):
        """Regression: integer median-combined master darks have a residual
        that is exactly zero in most pixels, so MAD = 0 — and the first
        version found NO hot pixels in six of twelve post-flip masters."""
        hot = hot_set()
        rng = np.random.default_rng(0)
        img = np.full(SHAPE, 303.0, dtype=np.float32)
        # 30% of pixels one count high: MAD is 0, the image is not constant.
        img[rng.random(SHAPE) < 0.3] += 1.0
        img[hot] += 500.0
        z = flip.residual_in_sigma(img)
        assert z is not None
        assert (z > flip.CANDIDATE_SIGMA)[hot].all()

    def test_persistence_needs_frames_of_one_geometry(self):
        with pytest.raises(ValueError):
            flip.persistent_mask([])
        with pytest.raises(ValueError):
            flip.persistent_mask([np.zeros((4, 4), bool),
                                  np.zeros((4, 5), bool)])


class TestTransformsAreToldApart:
    """One sensor, saved four ways.  Each saved orientation must be
    recognised as itself and as nothing else."""

    @pytest.mark.parametrize("saved_as, verdict", [
        ("identity", "physical_rotation"),   # hot pixels stay put
        ("rot180", "software_flip"),         # the saved array was rotated
        ("mirror_lr", "inconclusive"),       # neither story fits
    ])
    def test_verdict(self, saved_as, verdict):
        hot = hot_set()
        sb, lb = maps_of(hot, seeds=range(10, 18))
        sa, la = maps_of(hot, seeds=range(30, 38))
        fn = flip.TRANSFORMS[saved_as]
        sa, la = fn(sa), fn(la)              # what the files would hold
        rows = flip.compare_maps(sb, lb, sa, la, n_shifts=10)
        by = {r["transform"]: r["frac"] for r in rows}
        assert by[saved_as] == pytest.approx(1.0)
        for other in set(flip.TRANSFORMS) - {saved_as}:
            assert by[other] < 0.1
        assert by["shift_control_max"] < 0.1
        assert flip.verdict_of(rows, ceiling=1.0) == verdict

    def test_two_different_sensors_match_under_nothing(self):
        sb, lb = maps_of(hot_set(seed=1), seeds=range(10, 18))
        sa, la = maps_of(hot_set(seed=2), seeds=range(30, 38))
        rows = flip.compare_maps(sb, lb, sa, la, n_shifts=10)
        assert flip.verdict_of(rows, ceiling=1.0) == "inconclusive"
        assert max(r["frac"] for r in rows) < 0.1


class TestContainment:
    def test_a_warmer_noisier_side_does_not_break_the_match(self):
        """WHY CONTAINMENT.  The same sensor, but the 'after' side is taken
        with a brighter, noisier background (as the real post-shutdown
        frames were, at 0 C against -10 C), so its faint hot pixels fall
        below 8 sigma.  Strict-vs-strict coincidence drops; strict-in-loose
        containment does not."""
        ys, xs = hot_set(n=80)
        rng = np.random.default_rng(7)
        excess = rng.uniform(60.0, 400.0, size=80)     # a range of strengths
        hot = (ys, xs)
        sb, lb = maps_of(hot, seeds=range(10, 18), excess=excess, noise=5.0)
        sa, la = maps_of(hot, seeds=range(30, 38), excess=excess, noise=9.0)
        strict_coincidence = np.count_nonzero(sb & sa) / max(
            1, min(sb.sum(), sa.sum()))
        found, asked = flip.containment(sb, lb, sa, la)
        assert sa.sum() < sb.sum()                    # the sides DO differ
        assert found / asked > 0.95
        assert found / asked >= strict_coincidence

    def test_containment_is_symmetric_and_bounded(self):
        a = np.zeros(SHAPE, bool)
        b = np.zeros(SHAPE, bool)
        a[5, 5] = a[6, 6] = True
        b[5, 5] = True
        assert flip.containment(a, a, b, b) == flip.containment(b, b, a, a)
        found, asked = flip.containment(a, a, b, b)
        assert (found, asked) == (2, 3)

    def test_empty_maps_give_zero_not_a_division_error(self):
        z = np.zeros(SHAPE, bool)
        rows = flip.compare_maps(z, z, z, z, n_shifts=3)
        assert all(r["frac"] == 0.0 for r in rows)
        assert flip.verdict_of(rows, ceiling=0.0) == "inconclusive"


class TestVerdictRule:
    def _rows(self, **frac):
        base = dict(identity=0.0, rot180=0.0, mirror_lr=0.0, mirror_ud=0.0,
                    shift_control_mean=0.0, shift_control_max=0.0)
        base.update(frac)
        return [{"transform": k, "frac": v} for k, v in base.items()]

    def test_the_real_dark_arm_numbers_read_as_a_software_flip(self):
        # ASI 2025 monsoon boundary, dark arm, as measured 2026-10-03.
        rows = self._rows(identity=0.010, rot180=0.996, mirror_lr=0.010,
                          mirror_ud=0.010, shift_control_max=0.011)
        assert flip.verdict_of(rows, ceiling=1.0) == "software_flip"

    def test_a_match_no_better_than_chance_decides_nothing(self):
        rows = self._rows(rot180=0.05, shift_control_max=0.04)
        assert flip.verdict_of(rows, ceiling=1.0) == "inconclusive"

    def test_a_match_far_below_the_ceiling_decides_nothing(self):
        rows = self._rows(rot180=0.2, shift_control_max=0.001)
        assert flip.verdict_of(rows, ceiling=0.98) == "inconclusive"

    def test_two_transforms_both_matching_decide_nothing(self):
        rows = self._rows(identity=0.9, rot180=0.9)
        assert flip.verdict_of(rows, ceiling=1.0) == "inconclusive"
