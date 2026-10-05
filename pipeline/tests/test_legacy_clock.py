"""Unit tests for macro_legacy.clock — the legacy per-season clock audit.

The must-NOT cases matter most here: a prediction carried years beyond the
TESS span must be REFUSED, not timed against a guess; the two stamp
readings must differ by exactly half an exposure and in the right
direction; a fold must return the eclipse on the reference cycle, not on
some other one.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest \
        pipeline/tests/test_legacy_clock.py -q
"""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_legacy import clock as lk


class TestTargets:
    @pytest.mark.parametrize("t", ["EA/HW", "EA/HW+V361HYA", "EA/WD", "EA+UV"])
    def test_pceb_types_admitted(self, t):
        assert lk.is_clock_type(t)

    @pytest.mark.parametrize("t", ["EW", "EW/KW", "EA", "EB", "EP", None, ""])
    def test_contact_and_plain_eclipsers_refused(self, t):
        assert not lk.is_clock_type(t)

    def test_season_runs_july_to_june(self):
        assert lk.season_label("2018-12-31") == "2018/19"
        assert lk.season_label("2019-06-30") == "2018/19"
        assert lk.season_label("2019-07-01") == "2019/20"
        assert lk.season_label("2000-01-05") == "1999/00"


class TestFold:
    def test_fold_lands_every_eclipse_on_the_reference_cycle(self):
        p, t_ref = 0.1, 2459000.0
        t = t_ref + np.array([-5 * p, 0.0, 3 * p, 7 * p + 0.01])
        tf, keep = lk.fold_primary(t, t_ref, p)
        assert np.allclose(tf, [t_ref, t_ref, t_ref, t_ref + 0.01])
        assert keep.all()

    def test_secondary_is_outside_the_window(self):
        p, t_ref = 0.1, 2459000.0
        tf, keep = lk.fold_primary(np.array([t_ref + 0.5 * p]), t_ref, p)
        assert not keep[0]

    def test_bin_series_uses_empirical_errors(self):
        rng = np.random.default_rng(1)
        t = np.repeat(np.arange(10) * 1.0, 25) + rng.uniform(0, 0.5, 250)
        t[0] = 0.0          # bins are anchored at the earliest point
        f = 1.0 + rng.normal(0, 0.01, 250)
        tb, fb, eb, n = lk.bin_series(t, f, 1.0)
        assert len(tb) == 10 and (n == 25).all()
        assert np.allclose(eb, 0.01 / 5, rtol=0.5)

    def test_single_point_bins_are_dropped(self):
        tb, _fb, _eb, _n = lk.bin_series(np.array([0.0, 5.0, 5.1]),
                                         np.ones(3), 1.0)
        assert len(tb) == 1


class TestEphemeris:
    def _data(self, p=0.1, t0=2459000.0, quad=0.0, sig_s=2.0, seed=3):
        rng = np.random.default_rng(seed)
        e = np.array([0, 150, 3000, 3150, 6000, 6150])
        t = t0 + p * e + quad * e ** 2 + rng.normal(0, sig_s / 86400, len(e))
        return e, t, np.full(len(e), sig_s / 86400)

    def test_linear_fit_recovers_period(self):
        e, t, s = self._data()
        fit = lk.fit_ephemeris(e, t, s, order=1)
        assert fit["period"] == pytest.approx(0.1, abs=5 * fit["sig_period"])
        assert fit["dof"] == 4

    def test_predict_matches_truth_between_points(self):
        e, t, s = self._data(sig_s=0.5)
        fit = lk.fit_ephemeris(e, t, s)
        tp, sp = lk.predict(fit, 4500)
        assert abs(tp - (2459000.0 + 450.0)) * 86400 < 5 * sp * 86400 + 1

    def test_cycle_numbers(self):
        c = lk.cycle_numbers([2459000.0, 2459000.3001], 2459000.0, 0.1)
        assert list(c) == [0, 3]

    def test_interpolation_and_range_limits(self):
        e, t, s = self._data()
        inside = lk.reference_prediction(4500, t[0] + 450, t, s, e, 3)
        assert inside["range"] == "interpolated"
        near = lk.reference_prediction(-1000, t[0] - 100, t, s, e, 3)
        assert near["range"] == "extrapolated"
        far = lk.reference_prediction(-9000, t[0] - 900, t, s, e, 3)
        assert far["range"] == "out_of_range" and "t_pred" not in far

    def test_single_sector_reference_is_kept_short(self):
        e, t, s = np.array([0, 150]), np.array([2459000.0, 2459015.0]), \
            np.full(2, 2 / 86400)
        ok = lk.reference_prediction(500, 2459050.0, t, s, e, 1)
        assert ok["range"] == "extrapolated"
        no = lk.reference_prediction(1500, 2459150.0, t, s, e, 1)
        assert no["range"] == "out_of_range"

    def test_curvature_shows_up_as_the_systematic(self):
        e, t, s = self._data(quad=2e-9, sig_s=0.5)
        p = lk.reference_prediction(9000, t[0] + 900, t, s, e, 3)
        assert p["sys_basis"] == "linear vs quadratic"
        # 2e-9 d cycle^-2 bends the prediction by seconds-to-minutes at 9000.
        assert p["sig_sys_s"] > 5.0


class TestOC:
    def test_mid_reading_is_half_an_exposure_earlier(self):
        t_pred = 2459000.5
        t0 = t_pred + 30.0 / 86400
        s, m = lk.oc_both_conventions(t0, t_pred, exptime_s=64.0)
        assert s == pytest.approx(30.0)
        assert m == pytest.approx(-2.0)

    def test_verdicts(self):
        assert lk.season_verdict(10.0, 5.0) == "PASS"
        assert lk.season_verdict(200.0, 10.0) == "OFFSET MEASURED"
        assert lk.season_verdict(40.0, 15.0).startswith("PASS (central")
        assert lk.season_verdict(70.0, 40.0) == "UNRESOLVED"

    def test_acceptance_is_the_ledger_value(self):
        assert lk.CLOCK_ACCEPT_S == 60.0
