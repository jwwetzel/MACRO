"""Unit tests for the S2 bad-pixel masks (rlmt_diagnostics.badpix) and the
camera/configuration key (rlmt_diagnostics.camera).

The mask logic is tested on synthetic frame samples with planted defects:
hot pixels, the three "rails" an on-camera 2x2 average makes of saturated
native pixels, random-telegraph-signal pixels — and stars, which fire the
same per-frame test and must NOT end up in a mask because they never
repeat at one pixel across pointings.
"""

from __future__ import annotations

import numpy as np
import pytest

from rlmt_diagnostics import badpix as bp
from rlmt_diagnostics import camera as cam

SHAPE = (200, 240)
PED = 303.0


def sky_frame(rng, sky=40.0, n_stars=30):
    """A science frame: pedestal + sky + read noise + stars at random
    positions (different in every call — a different pointing)."""
    img = PED + sky + rng.normal(0.0, 3.0, SHAPE)
    yy, xx = np.mgrid[0:SHAPE[0], 0:SHAPE[1]]
    for _ in range(n_stars):
        y, x = rng.uniform(5, SHAPE[0] - 5), rng.uniform(5, SHAPE[1] - 5)
        img += rng.uniform(300, 20000) * np.exp(
            -((yy - y) ** 2 + (xx - x) ** 2) / (2 * 1.3 ** 2))
    return img


def plant(img, hot, rails):
    """Add the planted defects to one frame (in place) and return it."""
    for (y, x), level in hot.items():
        img[y, x] += level
    for (y, x), k in rails.items():
        local = img[y, x]
        img[y, x] = (k * bp.NATIVE_FULL_SCALE + (4 - k) * local) / 4.0
    return img


HOT = {(20, 30): 400.0, (50, 60): 900.0, (120, 200): 250.0}
RAILS = {(80, 90): 1, (81, 150): 2, (150, 40): 3}


class TestRailLevels:
    def test_levels_for_the_asi_pedestal(self):
        """OA.E8: hot pixels pile at (65535 + 3 x 304)/4 ~ 16.6 kADU."""
        lv = bp.rail_levels(304.0)
        assert lv[0] == pytest.approx(16611.75)
        assert lv[1] == pytest.approx(32919.5)
        assert lv[2] == pytest.approx(49227.25)

    def test_rail_index_uses_the_local_level(self):
        img = np.full((30, 30), 5000.0)             # bright sky
        img[10, 10] = (bp.NATIVE_FULL_SCALE + 3 * 5000.0) / 4.0
        img[20, 20] = bp.rail_levels(304.0)[0]      # rail for a DARK sky
        _exc, med = bp.neighbour_excess(img)
        k = bp.rail_index(img, med)
        assert k[10, 10] == 1
        assert k[20, 20] == 0       # not a rail at this sky level
        assert k.sum() == 1

    def test_full_saturation_is_not_a_rail(self):
        img = np.full((30, 30), 400.0)
        img[10, 10] = bp.NATIVE_FULL_SCALE
        _exc, med = bp.neighbour_excess(img)
        assert bp.rail_index(img, med)[10, 10] == 0


class TestNeighbourExcess:
    def test_isolates_a_single_pixel_spike(self):
        rng = np.random.default_rng(0)
        img = 300.0 + rng.normal(0, 3, SHAPE)
        img += np.linspace(0, 500, SHAPE[1])[None, :]      # a sky gradient
        img[100, 100] += 200.0
        exc, _med = bp.neighbour_excess(img)
        assert exc[100, 100] == pytest.approx(200.0, abs=15.0)
        assert bp.robust_sigma(exc) == pytest.approx(3.0 * np.sqrt(1.5),
                                                     rel=0.25)

    def test_two_adjacent_bad_pixels_are_both_found(self):
        img = np.full((40, 40), 300.0)
        img[20, 20] += 500.0
        img[20, 21] += 500.0
        exc, _ = bp.neighbour_excess(img)
        assert exc[20, 20] > 400 and exc[20, 21] > 400


class TestPersistence:
    def _accumulate(self, n_frames=24, n_native=4, seed=1):
        rng = np.random.default_rng(seed)
        acc = bp.PersistenceAccumulator(SHAPE, n_native)
        for _ in range(n_frames):
            acc.add(plant(sky_frame(rng), HOT, RAILS if n_native > 1 else {}))
        return acc

    def test_planted_defects_are_masked_and_stars_are_not(self):
        acc = self._accumulate()
        mask = acc.mask()
        for pos in HOT:
            assert mask[pos] & bp.BAD_HOT
        for pos in RAILS:
            assert mask[pos] & bp.BAD_RAIL
        # 24 frames x 30 stars fired the per-frame test hundreds of times;
        # none of them may persist.
        n_bad = int((mask > 0).sum())
        assert n_bad == len(HOT) + len(RAILS)

    def test_rail_index_is_recovered_per_pixel(self):
        acc = self._accumulate()
        rk = acc.rail_k()
        for pos, k in RAILS.items():
            assert rk[pos] == k

    def test_unbinned_camera_has_no_rail_test(self):
        acc = self._accumulate(n_native=1)
        assert not (acc.mask() & bp.BAD_RAIL).any()

    def test_a_pixel_bad_in_a_minority_of_frames_is_not_masked(self):
        rng = np.random.default_rng(2)
        acc = bp.PersistenceAccumulator(SHAPE, 1)
        for i in range(20):
            img = sky_frame(rng)
            if i < 6:                       # fires in 30% of frames
                img[33, 44] += 800.0
            acc.add(img)
        assert acc.mask()[33, 44] == 0
        assert acc.hot[33, 44] == 6         # ... but the count is kept

    def test_shape_mismatch_raises(self):
        acc = bp.PersistenceAccumulator(SHAPE, 1)
        with pytest.raises(ValueError):
            acc.add(np.zeros((10, 10)))


class TestTemporalNoise:
    def test_rts_pixel_is_found_and_quiet_hot_pixel_is_not(self):
        rng = np.random.default_rng(3)
        stack = PED + rng.normal(0, 2.0, (12,) + SHAPE)
        stack[:, 10, 10] += 5000.0                      # hot but steady
        stack[::2, 60, 70] += 60.0                      # telegraph jumps
        noisy, typical = bp.temporal_noise_mask(stack)
        assert typical == pytest.approx(2.0, rel=0.25)
        assert noisy[60, 70]
        assert not noisy[10, 10]
        assert noisy.sum() < 0.001 * noisy.size + 2

    def test_too_few_frames_raises(self):
        with pytest.raises(ValueError):
            bp.temporal_noise_mask(np.zeros((3, 8, 8)))


class TestAgreement:
    def test_counts_and_fractions(self):
        a = np.zeros((10, 10), bool)
        b = np.zeros((10, 10), bool)
        a[0, :4] = True
        b[0, 2:8] = True
        ag = bp.mask_agreement(a, b)
        assert (ag["n_a"], ag["n_b"], ag["n_both"]) == (4, 6, 2)
        assert ag["a_in_b"] == pytest.approx(0.5)
        assert ag["b_in_a"] == pytest.approx(2 / 6)
        assert bp.mask_agreement(np.zeros((2, 2), bool), b[:2, :2])["a_in_b"] \
            is None


class TestProductApi:
    def _write(self, tmp_path, key, shape, bad):
        mask = np.zeros(shape, np.uint8)
        for (y, x), flag in bad.items():
            mask[y, x] = flag
        return bp.save_mask(tmp_path, key, mask, np.zeros(shape, np.uint16),
                            None, {"mask_key": key})

    def test_round_trip_and_flag_selection(self, tmp_path):
        key = bp.mask_key("AC4040", 60, 40)
        self._write(tmp_path, key, (40, 60),
                    {(1, 2): bp.BAD_HOT, (3, 4): bp.BAD_NOISY})
        m = bp.load_mask("AC4040", (40, 60), directory=tmp_path)
        assert m.dtype == bool and m.sum() == 2
        only_hot = bp.load_mask("AC4040", (40, 60), flags=bp.BAD_HOT,
                                directory=tmp_path)
        assert only_hot[1, 2] and not only_hot[3, 4]
        assert bp.available_masks(tmp_path) == [key]
        assert bp.mask_metadata(key, tmp_path)["mask_key"] == key

    def test_unknown_camera_returns_none_not_an_empty_mask(self, tmp_path):
        assert bp.load_mask("QHY600", (40, 60), directory=tmp_path) is None
        self._write(tmp_path, bp.mask_key("AC4040", 60, 40), (40, 60), {})
        assert bp.load_mask("QHY600", (40, 60), directory=tmp_path) is None
        # Right camera, wrong geometry, no crop given: also None.
        assert bp.load_mask("AC4040", (30, 50), directory=tmp_path) is None

    def test_nearest_temperature_group_is_chosen(self, tmp_path):
        self._write(tmp_path, bp.mask_key("ASI", 60, 40, -10), (40, 60),
                    {(1, 1): bp.BAD_HOT})
        self._write(tmp_path, bp.mask_key("ASI", 60, 40, 0), (40, 60),
                    {(1, 1): bp.BAD_HOT, (2, 2): bp.BAD_HOT})
        assert bp.load_mask("ASI", (40, 60), temp_c=-9.0,
                            directory=tmp_path).sum() == 1
        assert bp.load_mask("ASI", (40, 60), temp_c=-1.0,
                            directory=tmp_path).sum() == 2
        # No temperature given: the coldest group.
        assert bp.load_mask("ASI", (40, 60), directory=tmp_path).sum() == 1

    def test_cropped_frame_gets_the_cropped_mask(self, tmp_path):
        self._write(tmp_path, bp.mask_key("QHY600", 60, 40), (40, 60),
                    {(20, 30): bp.BAD_RAIL})
        m = bp.load_mask("QHY600", (22, 47), crop=(18, 13),
                         directory=tmp_path)
        assert m.shape == (22, 47)
        assert m[2, 17] and m.sum() == 1

    def test_lookup_by_header_values(self, tmp_path):
        self._write(tmp_path, bp.mask_key("ASI", 60, 40, -10), (40, 60),
                    {(5, 5): bp.BAD_HOT})
        m = bp.mask_for_frame("ASI Camera (1)", 60, 40, "Mode0", -10.1,
                              directory=tmp_path)
        assert m is not None and m[5, 5]

    def test_sample_evenly(self):
        assert bp.sample_evenly(list(range(5)), 10) == list(range(5))
        s = bp.sample_evenly(list(range(100)), 5)
        assert s[0] == 0 and s[-1] == 99 and len(s) == 5


class TestCameraKey:
    """Four cameras, not one camera in six modes (DE memo, opening table)."""

    @pytest.mark.parametrize("instrume,mode,xb,n1,egain,expect", [
        ("DL Imaging", "High Gain", 1, 4096, 1.0537, "AC4040 High Gain e1.054"),
        ("DL Imaging", "High Gain", 1, 4096, 1.05697, "AC4040 High Gain e1.057"),
        ("DL Imaging", "High Gain StackPro", 1, 4096, 1.0537,
         "AC4040 StackPro e1.054"),
        ("DL Imaging", "Low Gain", 1, 4096, 22.6, "AC4040 Low Gain"),
        ("Andor CCD/EMCCD (SDK2)", "1MHz High Sensitivity 16-bit", 1, 2048,
         None, "iKon 1MHz"),
        ("Andor CCD/EMCCD (SDK2)", "5MHz High Sensitivity 16-bit", 1, 2048,
         0.0, "iKon 5MHz"),
        ("ASI Camera (1)", "Mode0", 2, 4788, 0.24666, "ASI Mode0 2x2"),
        ("ASI Camera (1)", "Mode0", 1, 9576, 0.24666, "ASI Mode0 1x1"),
        ("ASI Camera (1)", "Mode0", 2, 4788, 0.78, "ASI Mode0 2x2 e0.780"),
        ("QHYCCD-Cameras-Capture", "Fast", 2, 4800, 1.0, "QHY600 Fast 2x2"),
        ("QHY600Pro", None, 2, 4800, 56.0, "QHY600 pyscope 2x2"),
        (None, None, 2, 4800, 56.0, "QHY600 pyscope 2x2"),   # blank INSTRUME
    ])
    def test_config_key(self, instrume, mode, xb, n1, egain, expect):
        assert cam.config_key(instrume, mode, xb, n1, egain) == expect

    def test_mode0_and_fast_are_different_cameras(self):
        assert cam.camera_of("ASI Camera (1)") == "ASI"
        assert cam.camera_of("QHYCCD-Cameras-Capture") == "QHY600"

    def test_unknown_stays_unknown(self):
        assert cam.camera_of("Some Other Camera", 1234, None) == "unknown"
        assert cam.config_key("Some Other Camera", None, 1, 1234).startswith(
            "unknown")

    @pytest.mark.parametrize("config,mode", [
        ("AC4040 High Gain e1.054", "High Gain"),
        ("AC4040 High Gain e1.057", "High Gain"),
        ("AC4040 StackPro e1.057", "High Gain StackPro"),
        ("AC4040 Low Gain", "Low Gain"),
        ("iKon 1MHz", "1MHz High Sensitivity 16-bit"),
        ("ASI Mode0 2x2", "Mode0"),
        ("ASI Mode0 1x1", "Mode0"),
        ("QHY600 Fast 2x2", "Fast"),
        ("QHY600 pyscope 2x2", "(blank 2026)"),
    ])
    def test_config_maps_back_to_its_s2_mode(self, config, mode):
        assert cam.mode_label_of(config) == mode

    def test_native_pixel_count(self):
        assert cam.n_native_per_pixel("ASI Mode0 2x2") == 4
        assert cam.n_native_per_pixel("ASI Mode0 2x2 e0.780") == 4
        assert cam.n_native_per_pixel("ASI Mode0 1x1") == 1
        assert cam.n_native_per_pixel("AC4040 High Gain e1.057") == 1
        assert cam.camera_name("QHY600 Fast 2x2") == "QHY600"
