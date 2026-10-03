"""Tests for pipeline/scripts/ops_exposure.py — exposures for the QHY600.

The script's numbers end up in a request sent to the observatory, so the
arithmetic is pinned here on synthetic frames whose answer is known by
construction, and on closed forms:

* the blind and along-trace peak statistics recover an injected trace level
  and an injected emission line, and ignore isolated hot pixels (the
  IMX455's 16.5 kADU "rail" pixels) and a zero-order image;
* a trace polynomial that does not belong to the frame is caught by the
  contrast figure rather than silently measured;
* the Gaussian peak fraction integrates to the right total;
* the zero-point, scintillation and exposure formulas against hand values;
* the camera bridge takes the upper quartile, rejects twilight and clipped
  frames, and refuses to adopt a ratio from fewer than three frames.

No archive access: nothing here opens a FITS file or the manifest.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import ops_exposure as ox                        # noqa: E402


def synthetic_frame(level=5000.0, line=0.0, x_line=600, slope=0.02,
                    y0=150.0, shape=(300, 1200), sky=300.0, seed=1):
    """A sky pedestal plus one tilted trace of Gaussian cross-section
    (sigma 2 px) and constant peak ``level``, with an optional emission
    line of extra peak ``line`` (sigma 5 px along the trace)."""
    rng = np.random.default_rng(seed)
    ny, nx = shape
    y, x = np.mgrid[0:ny, 0:nx]
    yc = y0 + slope * x
    along = level + line * np.exp(-0.5 * ((x - x_line) / 5.0) ** 2)
    on = (x > 100) & (x < 1100)
    img = sky + np.where(on, along, 0.0) * np.exp(-0.5 * ((y - yc) / 2.0) ** 2)
    return (img + rng.normal(0, 3.0, shape)).astype(np.float32), (0.0, slope, y0)


class TestBlindPeak:
    def test_recovers_continuum_level(self):
        img, _ = synthetic_frame(level=5000.0)
        s = ox.trace_peak_blind(img)
        assert s["bg"] == pytest.approx(300.0, abs=2.0)
        # The 3x3 median of a sigma = 2 px profile sits a few % under peak.
        assert s["cont_p90"] == pytest.approx(5000.0, rel=0.06)

    def test_hot_pixels_do_not_count(self):
        img, _ = synthetic_frame(level=5000.0)
        rng = np.random.default_rng(7)
        ys, xs = rng.integers(0, 300, 200), rng.integers(0, 1200, 200)
        img[ys, xs] = 16_600.0               # the average-binning rail
        s = ox.trace_peak_blind(img)
        assert s["peak_max"] < 6000.0
        assert s["cont_p90"] == pytest.approx(5000.0, rel=0.06)

    def test_zero_order_does_not_move_the_continuum(self):
        img, _ = synthetic_frame(level=5000.0)
        img[40:50, 30:40] += 40_000.0        # a compact, much brighter blob
        s = ox.trace_peak_blind(img)
        assert s["peak_max"] > 30_000.0      # it IS the frame maximum ...
        assert s["cont_p90"] == pytest.approx(5000.0, rel=0.06)   # ... only


class TestAlongTrace:
    def test_line_and_continuum(self):
        img, coeffs = synthetic_frame(level=400.0, line=1200.0, x_line=600)
        s = ox.trace_peak_along(img, coeffs, 600.0)
        assert s["cont_p90"] == pytest.approx(400.0, rel=0.10)
        # A line narrow in BOTH directions loses ~12 % to the 3x3 median
        # (documented in ``despike``); it must not lose more.
        assert 0.85 * 1600.0 < s["line_peak"] <= 1600.0 * 1.02
        assert s["contrast"] > 10.0
        assert s["peak_max"] == s["line_peak"]

    def test_sky_lozenge_is_not_counted_as_starlight(self):
        """OA.E4/E7: the slitless sky is a bright lozenge under the trace,
        not the frame median.  900 ADU of it must not move a 400 ADU
        continuum."""
        img, coeffs = synthetic_frame(level=400.0, line=1200.0, x_line=600,
                                      shape=(400, 1200))
        y = np.arange(img.shape[0])[:, None]
        # Rows 70-250 of 400: under the trace and both sky bands, but less
        # than half the frame, so the frame median stays at the pedestal.
        img += (900.0 * (np.abs(y - 160) < 90)).astype(np.float32)
        s = ox.trace_peak_along(img, coeffs, 600.0)
        assert s["cont_p90"] == pytest.approx(400.0, rel=0.10)
        assert s["sky_local"] == pytest.approx(900.0, rel=0.02)
        assert s["contrast"] > 10.0

    def test_brighter_field_star_elsewhere_is_ignored(self):
        img, coeffs = synthetic_frame(level=400.0, line=1200.0)
        img[280:290, 500:520] += 30_000.0    # off the trace and sky bands
        s = ox.trace_peak_along(img, coeffs, 600.0)
        assert s["line_peak"] < 2000.0

    def test_wrong_trace_has_no_contrast(self):
        img, _ = synthetic_frame(level=400.0)
        s = ox.trace_peak_along(img, (0.0, 0.0, 70.0), 600.0)  # empty sky
        assert s["contrast"] < 3.0


class TestFormulas:
    def test_gaussian_peak_fraction_integrates_to_one(self):
        fwhm = 6.0
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        y, x = np.mgrid[-40:41, -40:41]
        g = np.exp(-0.5 * (x ** 2 + y ** 2) / sigma ** 2)
        assert ox.gaussian_peak_fraction(fwhm) == pytest.approx(
            g.max() / g.sum(), rel=1e-3)

    def test_rate_from_zmag(self):
        assert ox.rate_from_zmag(23.0, 23.0) == pytest.approx(1.0)
        assert ox.rate_from_zmag(18.0, 23.0) == pytest.approx(100.0)

    def test_scintillation_hand_value(self):
        # 50.8^(-2/3) = 0.07290; exp(-1515.7/8000) = 0.8274;
        # 0.09 * 0.07290 * 0.8274 / sqrt(2) = 3.839e-3 -> 4.17 mmag.
        assert ox.scintillation_mmag(1.0, 1.0, 50.8, 1515.7) \
            == pytest.approx(1085.7 * 3.839e-3, rel=2e-3)

    def test_scintillation_scalings(self):
        a = ox.scintillation_mmag(1.0, 1.0, 50.8, 1515.7)
        assert ox.scintillation_mmag(4.0, 1.0, 50.8, 1515.7) \
            == pytest.approx(a / 2.0)
        assert ox.scintillation_mmag(1.0, 2.0, 50.8, 1515.7) \
            == pytest.approx(a * 2 ** 1.75)

    def test_exposure_for_peak(self):
        assert ox.exposure_for_peak(2500.0, 25_000.0) == pytest.approx(10.0)

    def test_robust(self):
        med, q75, sig, n = ox.robust([1, 2, 3, 4, 100, np.nan])
        assert (med, n) == (3.0, 5)
        assert q75 == 4.0
        assert sig == pytest.approx(1.4826)
        assert ox.robust([])[3] == 0

    def test_evenly_is_deterministic_and_spans(self):
        rows = list(range(100))
        pick = ox.evenly(rows, 8)
        assert len(pick) == 8 and pick[0] == 0 and pick[-1] == 99
        assert ox.evenly(rows, 8) == pick
        assert ox.evenly([1, 2], 8) == [1, 2]


def _bridge_row(camera, cont, bg=300.0, peak=None, target="η Hya",
                filt="hrg"):
    return {"kind": "bridge", "target": target, "filter": filt,
            "camera": camera, "cont_p90": cont, "bg": bg,
            "peak_max": cont * 1.1 if peak is None else peak, "error": "",
            "line_peak": np.nan, "contrast": np.nan, "exptime": 60.0}


class TestBridge:
    def test_ratio_is_upper_quartile_over_upper_quartile(self):
        rows = [_bridge_row("ASI Mode0", v) for v in (8000, 9000, 10000, 5000)]
        rows += [_bridge_row("QHY MaxIm", v) for v in (4000, 4500, 5000, 2500)]
        b = ox.bridge_table(rows)[0]
        assert b["ASI Mode0"]["n"] == 4
        assert b["QHY MaxIm"]["ratio"] == pytest.approx(0.5, rel=1e-6)

    def test_twilight_and_clipped_frames_rejected(self):
        rows = [_bridge_row("ASI Mode0", 9000) for _ in range(4)]
        rows += [_bridge_row("QHY MaxIm", 4500) for _ in range(4)]
        rows.append(_bridge_row("QHY MaxIm", 64_000, bg=300.0, peak=64_800))
        rows.append(_bridge_row("QHY MaxIm", 20_000, bg=5000.0))
        b = ox.bridge_table(rows)[0]
        assert b["QHY MaxIm"]["n"] == 4
        assert b["QHY MaxIm"]["n_rejected"] == 2
        assert b["QHY MaxIm"]["ratio"] == pytest.approx(0.5)

    def test_no_star_frames_rejected(self):
        rows = [_bridge_row("ASI Mode0", 9000) for _ in range(4)]
        rows += [_bridge_row("QHY MaxIm", v) for v in (4500, 4500, 4500, 30)]
        b = ox.bridge_table(rows)[0]
        assert (b["QHY MaxIm"]["n"], b["QHY MaxIm"]["n_rejected"]) == (3, 1)

    def test_adopted_ratio_needs_three_frames(self):
        rows = [_bridge_row("ASI Mode0", 9000) for _ in range(4)]
        rows += [_bridge_row("QHY MaxIm", 4500) for _ in range(2)]
        bridge = ox.bridge_table(rows)
        assert ox.adopted_ratio(bridge, "hrg", "QHY MaxIm")[3] == 0

    def test_adopted_ratio_is_median_over_stars(self):
        rows = []
        for star, r in (("η Hya", 0.7), ("θ Vir", 0.8), ("θ Crt", 1.2)):
            rows += [_bridge_row("ASI Mode0", 10_000, target=star)] * 3
            rows += [_bridge_row("QHY MaxIm", 10_000 * r, target=star)] * 3
        med, lo, hi, n = ox.adopted_ratio(ox.bridge_table(rows), "hrg",
                                          "QHY MaxIm")
        assert (n, med) == (3, pytest.approx(0.8))
        assert (lo, hi) == (pytest.approx(0.7), pytest.approx(1.2))


def test_tcrb_levels_drop_low_contrast_frames():
    def row(cont, line, contrast):
        return {"kind": "tcrb", "filter": "hrg", "error": "",
                "cont_p90": cont, "line_peak": line, "peak_max": line,
                "sky_local": 50.0, "contrast": contrast, "exptime": 240.0}
    rows = [row(400, 1500, 20), row(420, 1700, 15), row(90, 95, 1.1)]
    t = ox.tcrb_levels(rows, "hrg")
    assert (t["n"], t["n_all"]) == (2, 3)
    assert t["line_max"] == 1700
    assert t["cont"][0] == pytest.approx(410.0)


class TestLadder:
    def test_brackets_both_ends_by_one_step(self):
        rungs = ox.ladder(1.0, 16.0, step=4.0, t_min=0.01, t_max=1000.0)
        assert rungs == [0.25, 1.0, 4.0, 16.0, 64.0]

    def test_clipped_to_limits(self):
        rungs = ox.ladder(0.01, 100.0, step=4.0, t_min=0.05, t_max=240.0)
        assert rungs[0] == 0.05 and rungs[-1] == 240.0
        assert all(b > a for a, b in zip(rungs, rungs[1:]))

    def test_faint_regime_ends_at_the_ceiling(self):
        # lo = 500 / 4 = 125 -> 120 (2 s.f., round-half-even); hi is
        # clipped to the ceiling.
        assert ox.ladder(500.0, 5000.0) == [120.0, ox.LADDER_MAX_S]

    def test_round_sig(self):
        assert ox.round_sig(0.0123) == 0.012
        assert ox.round_sig(187.3) == 190.0


def test_eruption_plan_scales_with_brightness(tmp_path):
    rate = {"hrg": 10.0, "lrg": 20.0, "g": 1000.0, "r": 2000.0, "i": 2500.0}
    plan = ox.eruption_plan(rate)
    assert [r["tag"] for r in plan] == ["1_peak", "2_early_decline",
                                        "3_late_decline"]
    # Brighter regime -> shorter exposures, rung for rung.
    assert plan[0]["ladders"]["hrg"][0] < plan[1]["ladders"]["hrg"][0] \
        < plan[2]["ladders"]["hrg"][0]
    # Twice the rate -> lrg no longer than hrg.
    assert plan[1]["ladders"]["lrg"][0] <= plan[1]["ladders"]["hrg"][0]
    # No direct imaging until the late decline.
    assert "g" not in plan[0]["ladders"] and "g" in plan[2]["ladders"]
    # V = 4 on hrg: 25,000 / (10 * 10**(0.4*6)) = 9.95 s; the peak regime's
    # faint end, bracketed by x4, must reach it.
    assert plan[0]["ladders"]["hrg"][-1] == pytest.approx(40.0, rel=0.05)
    paths = ox.write_eruption_block(plan, tmp_path)
    assert len(paths) == 4
    sch = (tmp_path / "tcrb_eruption_1_peak.sch").read_text()
    assert sch.startswith("# DRAFT") and "block start" in sch
    assert "filter g " not in sch
    n_csv = len((tmp_path / "tcrb_eruption_exposures.csv")
                .read_text().strip().splitlines()) - 1
    assert n_csv == sum(len(v) for r in plan for v in r["ladders"].values())
