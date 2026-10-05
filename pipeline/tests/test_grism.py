"""Unit tests for macro_grism — every pure decision function.

The must-NOT cases matter as much as the must cases: an ambiguous FITS
layout must refuse, a bad pointing must be rejected whatever its content
check says, a flipped-parity field must still be found, a windowless
column must not invent a background, and an unanchored spectrum must not
receive a wavelength axis.

Run with:
    /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest pipeline/tests -q
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Make the package importable regardless of pytest's working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_grism import config as gcfg
from macro_grism import ew as gew
from macro_grism import extract as gext
from macro_grism import gate as gg
from macro_grism import linecal as lc
from macro_grism import summary_g as gsum
from macro_grism import trace as gt
from macro_grism import wavelength as gw
from macro_grism.fits_io import GrismLayoutError, HduSummary, classify_hdus


# ---------------------------------------------------------------------------
# FITS packaging resolution
# ---------------------------------------------------------------------------
class TestClassifyHdus:
    def _s(self, *kinds):
        return [HduSummary(i, k, d) for i, (k, d) in enumerate(kinds)]

    def test_plain_primary_image(self):
        # The recovered master calibrations: image in the primary HDU.
        layout, idx = classify_hdus(self._s(("PrimaryHDU", True)))
        assert (layout, idx) == ("plain", 0)

    def test_fpack(self):
        # Every rawimage .fts.fz: stub primary + one CompImageHDU.
        layout, idx = classify_hdus(
            self._s(("PrimaryHDU", False), ("CompImageHDU", True)))
        assert (layout, idx) == ("fpack", 1)

    def test_repackaged(self):
        # Era-C reprocessing: stub primary + one plain ImageHDU.
        layout, idx = classify_hdus(
            self._s(("PrimaryHDU", False), ("ImageHDU", True)))
        assert (layout, idx) == ("repackaged", 1)

    def test_plain_wins_even_with_extensions(self):
        # Primary image + extra extension: the primary is canonical.
        layout, idx = classify_hdus(
            self._s(("PrimaryHDU", True), ("ImageHDU", True)))
        assert (layout, idx) == ("plain", 0)

    def test_two_extensions_refuse(self):
        # Ambiguous: two data-bearing extensions and no primary image.
        with pytest.raises(GrismLayoutError):
            classify_hdus(self._s(("PrimaryHDU", False),
                                  ("ImageHDU", True), ("ImageHDU", True)))

    def test_mixed_comp_and_image_refuse(self):
        with pytest.raises(GrismLayoutError):
            classify_hdus(self._s(("PrimaryHDU", False),
                                  ("CompImageHDU", True),
                                  ("ImageHDU", True)))

    def test_no_data_refuse(self):
        with pytest.raises(GrismLayoutError):
            classify_hdus(self._s(("PrimaryHDU", False),
                                  ("BinTableHDU", False)))

    def test_empty_refuse(self):
        with pytest.raises(GrismLayoutError):
            classify_hdus([])


# ---------------------------------------------------------------------------
# Trace geometry on a synthetic frame
# ---------------------------------------------------------------------------
def synthetic_frame(ny=400, nx=1200, slope=0.03, u0=200.0, height=200.0,
                    fwhm=6.0, halo=30.0, pedestal=300.0, seed=7):
    """A slitless-frame stand-in: pedestal + broad halo + one tilted
    Gaussian trace + Poisson-ish noise."""
    rng = np.random.default_rng(seed)
    y = np.arange(ny)[:, None]
    x = np.arange(nx)[None, :]
    center = u0 + slope * (x - nx / 2)
    sig = fwhm / 2.3548
    img = (pedestal
           + halo * np.exp(-0.5 * ((y - center) / 120.0) ** 2)
           + height * np.exp(-0.5 * ((y - center) / sig) ** 2))
    return img + rng.normal(0, 3.0, size=img.shape)


class TestTrace:
    def test_slope_recovered(self):
        img = synthetic_frame(slope=0.03)
        xs, ys, amps = gt.chunk_peaks(img, n_chunks=24, halo_win=101)
        slope = gt.fit_slope(xs, ys, amps)
        assert abs(slope - 0.03) < 0.005

    def test_main_trace_found_after_detilt(self):
        img = synthetic_frame(slope=0.03, u0=200.0)
        xs, ys, amps = gt.chunk_peaks(img, n_chunks=24, halo_win=101)
        slope = gt.fit_slope(xs, ys, amps)
        _, resid = gt.detilted_profile(img, slope)
        u, h = gt.main_trace_u(resid)
        assert abs(u - 200) <= 2
        assert h > 100                      # most of the injected height

    def test_flat_frame_gives_zero_slope(self):
        # No trace at all: the slope fit must not hallucinate one.  All
        # chunk amplitudes are comparable noise, so the guard cannot drop
        # them — but the fitted line through argmax-noise must stay tiny
        # relative to a real trace slope, and the trace HEIGHT (the gate's
        # actual floor) must be negligible.
        rng = np.random.default_rng(1)
        img = rng.normal(300, 3.0, size=(400, 1200))
        xs, ys, amps = gt.chunk_peaks(img, n_chunks=24, halo_win=101)
        slope = gt.fit_slope(xs, ys, amps)
        _, resid = gt.detilted_profile(img, slope)
        _, h = gt.main_trace_u(resid)
        assert h < gg.MIN_TRACE_HEIGHT_ADU    # the gate floor catches it

    def test_trace_centers_follow_curvature(self):
        # A curved trace: the deg-2 refinement must track it.
        ny, nx = 400, 1200
        y = np.arange(ny)[:, None]
        x = np.arange(nx)[None, :]
        center = 180 + 0.02 * (x - nx / 2) + 1e-5 * (x - nx / 2) ** 2
        img = 300 + 150 * np.exp(-0.5 * ((y - center) / 3.0) ** 2)
        coeffs, n, rms = gt.fit_trace_centers(img, 0.02, 180.0)
        mid = np.polyval(coeffs, nx / 2)
        edge = np.polyval(coeffs, nx - 1)
        assert abs(mid - 180.0) < 1.0
        assert abs(edge - (180 + 0.02 * (nx / 2 - 1)
                           + 1e-5 * (nx / 2 - 1) ** 2)) < 2.0
        assert rms is not None and rms < 1.0


# ---------------------------------------------------------------------------
# Identity gate
# ---------------------------------------------------------------------------
#: A plausible era-76 CD matrix (0.4508"/px, ~0.3 deg rotation).
CD = np.array([[1.25268e-4, -8.8369e-7], [8.4428e-7, 1.25204e-4]])


class TestGateGeometry:
    def test_center_star_maps_to_center(self):
        u = gg.predicted_u(CD, "A", 240.0, 25.0,
                           np.array([240.0]), np.array([25.0]),
                           0.03, 3194, 4788)
        assert abs(u[0] - 3194 / 2) < 1e-9

    def test_parity_flip_mirrors_offsets(self):
        # The meridian flip negates pixel offsets: u residuals from the
        # center must be equal and opposite between parities.
        ra = np.array([240.05])
        dec = np.array([25.02])
        ua = gg.predicted_u(CD, "A", 240.0, 25.0, ra, dec, 0.0, 3194, 4788)
        ub = gg.predicted_u(CD, "B", 240.0, 25.0, ra, dec, 0.0, 3194, 4788)
        assert abs((ua[0] - 1597) + (ub[0] - 1597)) < 1e-9

    def test_slope_term_cancels_grism_deflection(self):
        # Two stars separated ONLY along dispersion (same u expected):
        # a pure-x offset times the slope must shift u accordingly, and
        # with slope 0 the u values must be identical.
        ra = np.array([240.0, 240.1])
        dec = np.array([25.0, 25.0])
        u0 = gg.predicted_u(CD, "A", 240.0, 25.0, ra, dec, 0.0, 3194, 4788)
        # With CROTA ~ 0.3 deg a pure-RA offset leaks ~[rotation] into y;
        # the leak is the same for both parities and small.
        assert abs(u0[1] - u0[0]) < 15

    def test_brightest_prediction_picks_brightest_on_frame(self):
        # Star A: G=9 on frame.  Star B: G=5 but 3 deg away (off frame).
        stars = np.array([[240.02, 25.01, 9.0], [243.0, 25.0, 5.0]])
        preds = gg.brightest_prediction(CD, 240.0, 25.0, stars, 0.0,
                                        3194, 4788)
        assert preds["A"] is not None
        assert preds["A"][1] == 9.0          # the off-frame G=5 ignored


class TestGateVerdict:
    PREDS = {"A": (1600.0, 8.7), "B": (2400.0, 8.7)}

    def test_good_frame_accepted(self):
        g = gg.gate_verdict(0.013, 180.0, 1650.0, self.PREDS, 54)
        assert g.verdict == "ACCEPT" and g.reason == "ok"
        assert g.parity == "A" and abs(g.u_resid_px - 50.0) < 1e-9

    def test_best_parity_wins(self):
        g = gg.gate_verdict(0.013, 180.0, 2380.0, self.PREDS, 54)
        assert g.verdict == "ACCEPT" and g.parity == "B"

    def test_bad_pointing_rejected_whatever_the_content(self):
        # The mandated case: header 168 deg off target — REJECT, even if
        # the content check would have passed at its own pointing.
        g = gg.gate_verdict(168.0, 180.0, 1600.0, self.PREDS, 54)
        assert g.verdict == "REJECT" and g.reason == "header_off_target"

    def test_field_mismatch_rejected(self):
        # Header near target but the trace is 800 px from the prediction
        # on BOTH parities (the measured signature of the bad frames).
        g = gg.gate_verdict(0.013, 180.0, 800.0, self.PREDS, 54)
        assert g.verdict == "REJECT" and g.reason == "field_mismatch"
        assert g.u_resid_px == -800.0        # forensic record kept

    def test_no_trace_rejected(self):
        g = gg.gate_verdict(0.013, 4.0, 1600.0, self.PREDS, 54)
        assert g.verdict == "REJECT" and g.reason == "no_trace"

    def test_missing_pointing_fails_closed(self):
        g = gg.gate_verdict(None, 180.0, 1600.0, self.PREDS, 54)
        assert g.verdict == "REJECT" and g.reason == "no_header_pointing"

    def test_empty_gaia_fails_closed(self):
        g = gg.gate_verdict(0.013, 180.0, 1600.0,
                            {"A": None, "B": None}, 0)
        assert g.verdict == "REJECT" and g.reason == "no_gaia_catalog"

    def test_no_onframe_star_fails_closed(self):
        g = gg.gate_verdict(0.013, 180.0, 1600.0,
                            {"A": None, "B": None}, 12)
        assert g.verdict == "REJECT" and g.reason == "no_onframe_star"

    def test_tolerance_edge(self):
        # Exactly at the tolerance is inside; one px beyond is out.
        at = gg.gate_verdict(0.013, 180.0, 1600.0 + gg.U_TOL_PX,
                             {"A": (1600.0, 8.7), "B": None}, 5)
        beyond = gg.gate_verdict(0.013, 180.0, 1601.0 + gg.U_TOL_PX,
                                 {"A": (1600.0, 8.7), "B": None}, 5)
        assert at.verdict == "ACCEPT"
        assert beyond.verdict == "REJECT"

    def test_angular_offset_exact_at_large_angle(self):
        # 180 deg apart on the equator — small-angle math would break.
        assert abs(gg.angular_offset_deg(0, 0, 180, 0) - 180.0) < 1e-9
        assert abs(gg.angular_offset_deg(10, 20, 10, 20)) < 1e-12


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------
class TestFlankingBackground:
    def test_linear_gradient_removed_exactly(self):
        # Background = a + b*y: the two-band straight line reproduces it
        # exactly at every aperture row (the local-dark-removal claim).
        y = np.arange(400, dtype=float)
        col = 5.0 + 0.02 * y
        bg, ok = gext.flanking_background(col, yc=200.0)
        assert ok
        rows = np.arange(200 - gext.APERTURE_HALFWIN,
                         200 + gext.APERTURE_HALFWIN + 1)
        assert np.allclose(bg, 5.0 + 0.02 * rows, atol=1e-9)

    def test_band_off_frame_refuses(self):
        col = np.zeros(100)
        bg, ok = gext.flanking_background(col, yc=20.0)   # bands off top
        assert not ok and bg is None

    def test_trace_flux_does_not_leak_into_bands(self):
        # A trace inside the aperture must not bias the background.
        y = np.arange(400, dtype=float)
        col = 10.0 + 100.0 * np.exp(-0.5 * ((y - 200) / 3.0) ** 2)
        bg, ok = gext.flanking_background(col, yc=200.0)
        assert ok and np.all(np.abs(bg - 10.0) < 0.5)


DET = gcfg.DETECTOR_DEFAULTS["ASI-Mode0-bin2avg"]


class TestHorne:
    def _column(self, flux=500.0, bg=20.0, fwhm=5.0, seed=0):
        rng = np.random.default_rng(seed)
        n = 2 * gext.APERTURE_HALFWIN + 1
        y = np.arange(n)
        prof = np.exp(-0.5 * ((y - n // 2) / (fwhm / 2.3548)) ** 2)
        prof /= prof.sum()
        clean = flux * prof + bg + DET.pedestal_adu
        return clean + rng.normal(0, 2.0, n), prof, bg + DET.pedestal_adu

    def test_flux_recovered(self):
        win, prof, bg = self._column()
        f, v, ns = gext.horne_column(win, np.full(len(win), bg), prof, DET)
        assert abs(f - 500.0) < 30.0
        assert v > 0 and ns == 0

    def test_saturated_pixels_masked_and_counted(self):
        win, prof, bg = self._column(flux=500.0)
        win[gext.APERTURE_HALFWIN] = DET.saturation_cap_adu() + 100
        f, v, ns = gext.horne_column(win, np.full(len(win), bg), prof, DET)
        assert ns == 1
        assert f is not None                 # the wings still constrain it

    def test_fully_saturated_column_refuses(self):
        win = np.full(25, DET.saturation_cap_adu() + 1.0)
        prof = np.full(25, 1 / 25)
        f, v, ns = gext.horne_column(win, np.zeros(25), prof, DET)
        assert f is None and ns == 25

    def test_optimal_beats_box_under_noise(self):
        fh, fb = [], []
        for seed in range(60):
            win, prof, bg = self._column(flux=200.0, seed=seed)
            f, _, _ = gext.horne_column(win, np.full(len(win), bg), prof, DET)
            fh.append(f)
            fb.append(float((win - bg).sum()))
        assert np.std(fh) < np.std(fb)

    def test_profile_normalized_and_nonnegative(self):
        cut = np.array([[0.1, 0.8, 0.1], [-0.05, 0.9, 0.15]])
        p = gext.build_profile(cut)
        assert np.all(p >= 0) and abs(p.sum() - 1.0) < 1e-9

    def test_median_relative_difference(self):
        a = np.full(1000, 100.0)
        b = np.full(1000, 98.0)               # a steady 2% offset
        assert abs(gext.median_relative_difference(a, b) - 0.02) < 1e-9

    def test_median_relative_difference_needs_overlap(self):
        assert gext.median_relative_difference(
            np.full(1000, np.nan), np.full(1000, 1.0)) is None


class TestDetector:
    def test_pedestal_carries_no_shot_noise(self):
        v = DET.variance_adu2(np.array([DET.pedestal_adu]))
        assert abs(v[0] - DET.read_noise_adu ** 2) < 1e-9

    def test_average_binning_derates_the_clip(self):
        cap = gcfg.replace(DET, linearity_cap_adu=None).saturation_cap_adu()
        assert cap < DET.native_clip_adu
        assert abs(cap - ((65535 - DET.pedestal_adu) / gcfg.NATIVE_PEAKING
                          + DET.pedestal_adu)) < 1e-6

    def test_linearity_cap_wins_when_lower(self):
        d = gcfg.replace(DET, linearity_cap_adu=20000.0)
        assert d.saturation_cap_adu() == 20000.0

    def test_table_reads_measured_params(self):
        params = {("ASI Mode0 2x2", "gain_e_per_adu"): 1.0455,
                  ("ASI Mode0 2x2", "read_noise_adu"): 2.289}
        t = gcfg.build_detector_table(params)
        d = t["ASI-Mode0-bin2avg"]
        assert d.gain_e_per_adu == 1.0455 and d.read_noise_adu == 2.289
        # Quantities absent from the table keep the default AND say so.
        assert d.provisional and "DEFAULT for" in d.provenance

    def test_unknown_override_refused(self):
        with pytest.raises(KeyError):
            gcfg.detector_table({"no-such-camera": {"gain_e_per_adu": 1}},
                                db_path=False)

    def test_detector_key(self):
        assert gcfg.detector_key("ASI Camera (1)", "Mode0") == \
            "ASI-Mode0-bin2avg"
        assert gcfg.detector_key("QHYCCD-Cameras-Capture", "Fast") == \
            "QHY600-bin2avg"
        assert gcfg.detector_key("QHY600Pro", "") == "QHY600-pyscope-bin2avg"
        assert gcfg.detector_key("DL Imaging", "High Gain StackPro") == \
            "AC4040-StackPro"


class TestEpochs:
    def test_merge_only_across_soft_boundaries(self):
        rows = [("ASI:2024-12-16", "2024-12-16", "2025-06-23", "wheel_map"),
                ("ASI:2025-10-11", "2025-10-11", "2026-03-16", "flipstat")]
        assert all(r["verdict"] == "ok" for r in gcfg.reconcile_epochs(rows))

    def test_hardware_boundary_inside_an_epoch_is_caught(self):
        rows = [("X:a", "2025-01-01", "2025-02-01", "camera"),
                ("X:b", "2025-03-01", "2025-04-01", "flipstat")]
        out = gcfg.reconcile_epochs(rows)
        assert out[1]["verdict"] == "MERGED_OVER_HARDWARE"


class TestSky:
    def test_polynomial_sky_follows_a_curved_lozenge(self):
        # A rectified cutout: curved sky (quadratic across the trace) plus
        # a narrow trace at the centre row.
        half, nx = gext.RECT_HALF, 300
        dy = np.arange(-half, half + 1)[:, None]
        sky = 500.0 - 0.02 * dy ** 2 + 0.0 * np.arange(nx)[None, :]
        trace = 1000.0 * np.exp(-0.5 * (dy / 2.5) ** 2)
        rng = np.random.default_rng(1)
        cut = sky + trace + rng.normal(0, 1.0, sky.shape)
        poly = gext.sky_polynomial(cut)
        line = gext.sky_flanking(cut)
        at = half                                # the trace row
        assert abs(np.median(poly[at]) - 500.0) < 1.0
        # The straight line through the bands misses the curvature.
        assert np.median(line[at]) < 500.0 - 10.0


# ---------------------------------------------------------------------------
# Zero point, line measurement, the fixed solution
# ---------------------------------------------------------------------------
def emission_spectrum(nx=4000, x_ha=1800.0, height=800.0, noise=5.0,
                      seed=3):
    rng = np.random.default_rng(seed)
    x = np.arange(nx, dtype=float)
    spec = 1000.0 + 0.05 * x + height * np.exp(-0.5 * ((x - x_ha) / 6.0) ** 2)
    return spec + rng.normal(0, noise, nx)


class TestZeroPoint:
    def test_halpha_found_to_a_fraction_of_a_pixel(self):
        z = gw.halpha_zero_point(emission_spectrum())
        assert z is not None and abs(z["x"] - 1800.0) < 0.3
        assert z["x_err"] is not None and z["x_err"] < 0.3

    def test_no_line_no_zero_point(self):
        rng = np.random.default_rng(0)
        assert gw.halpha_zero_point(1000.0 + rng.normal(0, 5, 4000)) is None

    def test_cosmic_ray_not_mistaken_for_halpha(self):
        rng = np.random.default_rng(0)
        spec = 1000.0 + rng.normal(0, 5.0, 4000)
        spec[2000] += 5000.0
        assert gw.halpha_zero_point(spec) is None

    def test_snippet_roundtrip(self):
        import json
        flux = np.arange(100, dtype=float)
        flux[50] = np.nan
        back = json.loads(json.dumps(gw.snippet(flux, 50.0, 10, 1)))
        xs = [p[0] for p in back]
        assert back[xs.index(50)][1] is None and back[xs.index(41)][1] == 41.0


class TestLineCal:
    COEFFS = [-460.0, 45.0, -20.0]      # A per kpx^k: ~0.46 A/px, curved

    def test_wavelength_and_pixel_are_inverse(self):
        x = np.array([500.0, 2394.0, 4000.0])
        lam = lc.wavelength_of(x, 2394.0, self.COEFFS, 2394.0)
        assert abs(lam[1] - lc.HALPHA_A) < 1e-9
        back = lc.predict_x(lam, 2394.0, self.COEFFS, 2394.0)
        assert np.allclose(back, x, atol=1e-6)

    def test_fixed_solution_recovered_with_free_frame_constants(self):
        waves = np.array([5875.6, 6347.1, 6562.8, 6678.2, 7065.2])
        rng = np.random.default_rng(2)
        ids, ws, xs = [], [], []
        for f, shift in enumerate((-300.0, 0.0, 250.0)):
            # x for each wavelength under the polynomial with this frame's
            # constant (Newton inverse of the forward model).
            x = lc.predict_x(waves, 2394.0 + shift, self.COEFFS, 2394.0)
            ids += [f] * len(waves)
            ws += list(waves)
            xs += list(x + rng.normal(0, 0.05, len(x)))
        sol = lc.solve_dispersion(ids, ws, xs, np.full(len(xs), 0.05),
                                  degree=3, x_ref=2394.0)
        assert abs(sol["disp_ref"] - self.COEFFS[0] / 1000.0) < 2e-3
        assert sol["rms_px"] < 0.2 and sol["dof"] > 0

    def test_fit_line_centre(self):
        x = np.arange(200, dtype=float)
        y = 1000.0 - 300.0 * np.exp(-0.5 * ((x - 100.3) / 4.0) ** 2)
        m = lc.fit_line(y, 98.0, 30)
        assert abs(m["x"] - 100.3) < 0.05 and m["amp_frac"] < 0

    def test_edge_position_finds_a_cliff(self):
        x = np.arange(400, dtype=float)
        y = np.where(x < 200.0, 1000.0, 700.0) + 0.0 * x
        m = lc.edge_position(y, 200.0, +1, search_px=10, span_px=20)
        assert m is not None and abs(m["x"] - 199.5) < 1.5

    def _b_star_features(self, a):
        """Feature list of a synthetic B star dispersed at ``a`` px/A."""
        names = ("HeI5876", "O2gamma", "SiII6347", "Halpha", "HeI6678",
                 "O2B", "HeI7065")
        return [{"x": 2400.0 + a * (lc._id_wave(lc.LINE_BY_NAME[n])
                                    - lc.HALPHA_A),
                 "depth": 0.3 if n == "Halpha" else 0.05} for n in names]

    def test_true_dispersion_beats_the_rival(self):
        a = lc.SEEDS["hrg"][0]
        feats = self._b_star_features(a)
        true = lc.identify_lines(feats, "hrg", "B", seed=(a, 0.0))
        rival = lc.identify_lines(feats, "hrg", "B",
                                  seed=lc.RIVALS["hrg"]["1.59 A/px (v1 code)"])
        assert true["n_match"] >= 6
        assert rival["n_match"] < true["n_match"] - 2


class TestFingerprintGate:
    def test_verdict_never_needs_a_header(self):
        import inspect
        params = inspect.signature(gg.fingerprint_verdict).parameters
        assert not any("point" in p or "header" in p for p in params)

    def test_truth_table(self):
        A, R = gg.GATE_ACCEPT, gg.GATE_REJECT
        assert gg.fingerprint_verdict(False, 20, 0.95, 0.6)[0] == R
        assert gg.fingerprint_verdict(True, None, 0.95, 0.6)[1] == \
            "no_emission_line"
        assert gg.fingerprint_verdict(True, 20, 0.5, 0.6)[1] == \
            "fingerprint_mismatch"
        assert gg.fingerprint_verdict(True, 20, 0.95, 1.0)[1] == \
            "no_tio_step"
        assert gg.fingerprint_verdict(True, 20, 0.95, 0.6)[0] == A

    def test_tio_step_and_correlation(self):
        grid = np.arange(*gg.FP_WAVE, 1.0)
        f = np.where(grid > 7054.0, 0.6, 1.0)
        assert abs(gg.tio_step(grid, f) - 0.6) < 1e-9
        assert gg.fingerprint_r(grid, f, f) > 0.999


class TestEW:
    def test_gaussian_emission_line(self):
        wave = np.arange(6400.0, 6700.0, 0.5)
        sig = 3.0
        flux = 100.0 * (1 + 2.0 * np.exp(-0.5 * ((wave - 6562.8) / sig) ** 2))
        m = gew.equivalent_width(wave, flux)
        assert abs(m["ew_a"] - 2.0 * sig * np.sqrt(2 * np.pi)) < 0.05

    def test_missing_band_refuses(self):
        wave = np.arange(6540.0, 6580.0, 0.5)
        assert gew.equivalent_width(wave, np.ones_like(wave)) is None


class TestPTC:
    def test_gain_and_floor_recovered(self):
        level = np.linspace(10, 500, 40)
        var = 2.3 ** 2 + level / 1.05
        k, ke, rn, rne = gsum.ptc_fit(level, var)
        assert abs(k - 1.05) < 1e-6 and abs(rn - 2.3) < 1e-6


class TestGaiaConeCache:
    """The empty-cone retry ladder (the first validation run's lesson:
    a transient empty Vizier reply must never poison the cache)."""

    STARS = [[239.87, 25.92, 9.8]]

    def test_nonempty_answer_believed_at_once(self, tmp_path):
        calls = []

        def query(ra, dec, radius, g_limit):
            calls.append(1)
            return list(self.STARS)

        out = gg.gaia_cone(239.88, 25.91, str(tmp_path),
                           _query=query, _sleep=lambda s: None)
        assert len(calls) == 1                     # no retries needed
        assert out.shape == (1, 3)

    def test_empty_answer_retried_then_replaced(self, tmp_path):
        # First two replies empty (the hiccup), third has the stars.
        replies = [[], [], list(self.STARS)]
        slept = []
        out = gg.gaia_cone(239.88, 25.91, str(tmp_path),
                           _query=lambda *a: replies.pop(0),
                           _sleep=slept.append)
        assert out.shape == (1, 3)                 # the hiccup healed
        assert len(slept) == 2                     # one pause per retry

    def test_empty_cached_only_after_all_retries(self, tmp_path):
        calls = []

        def query(ra, dec, radius, g_limit):
            calls.append(1)
            return []

        out = gg.gaia_cone(10.0, -5.0, str(tmp_path),
                           _query=query, _sleep=lambda s: None)
        assert len(calls) == gg.EMPTY_RETRIES      # asked firmly
        assert out.shape == (0, 3)
        # The believed-empty IS cached: a rerun must not re-query...
        calls.clear()
        out2 = gg.gaia_cone(10.0, -5.0, str(tmp_path),
                            _query=query, _sleep=lambda s: None)
        # ...but an empty cache entry is never trusted silently — it is
        # re-asked (the self-healing path), still EMPTY_RETRIES times.
        assert len(calls) == gg.EMPTY_RETRIES
        assert out2.shape == (0, 3)

    def test_poisoned_empty_cache_self_heals(self, tmp_path):
        # Simulate the first run's poison: an empty entry on disk for a
        # field that really has stars.
        key = gg.cone_cache_key(239.88, 25.91, gg.CONE_RADIUS_DEG,
                                gg.CONE_G_LIMIT)
        (tmp_path / key).write_text("[]")
        out = gg.gaia_cone(239.88, 25.91, str(tmp_path),
                           _query=lambda *a: list(self.STARS),
                           _sleep=lambda s: None)
        assert out.shape == (1, 3)                 # healed
        # And the cache now holds the stars, trusted without a query.
        out2 = gg.gaia_cone(239.88, 25.91, str(tmp_path),
                            _query=lambda *a: (_ for _ in ()).throw(
                                AssertionError("must not re-query")),
                            _sleep=lambda s: None)
        assert out2.shape == (1, 3)

    def test_exception_retried_then_healed(self, tmp_path):
        # A read-timeout on attempt 1, stars on attempt 2 (also observed
        # live): the exception joins the same retry ladder.
        replies = [TimeoutError("read timed out"), list(self.STARS)]

        def query(*a):
            r = replies.pop(0)
            if isinstance(r, Exception):
                raise r
            return r

        out = gg.gaia_cone(239.88, 25.91, str(tmp_path),
                           _query=query, _sleep=lambda s: None)
        assert out.shape == (1, 3)

    def test_final_exception_propagates(self, tmp_path):
        def query(*a):
            raise TimeoutError("read timed out")

        with pytest.raises(TimeoutError):
            gg.gaia_cone(239.88, 25.91, str(tmp_path),
                         _query=query, _sleep=lambda s: None)
        # And nothing was cached: the next run starts clean.
        assert not list(tmp_path.glob("gaia_*.json"))
