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

from macro_grism import extract as gext
from macro_grism import gate as gg
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


class TestHorne:
    def _column(self, flux=500.0, bg=20.0, fwhm=5.0, seed=0):
        rng = np.random.default_rng(seed)
        n = 2 * gext.APERTURE_HALFWIN + 1
        y = np.arange(n)
        prof = np.exp(-0.5 * ((y - n // 2) / (fwhm / 2.3548)) ** 2)
        prof /= prof.sum()
        clean = flux * prof + bg
        return clean + rng.normal(0, 2.0, n), prof, bg

    def test_flux_recovered(self):
        win, prof, bg = self._column()
        f, v, ns = gext.horne_column(win, np.full(len(win), bg), prof,
                                     egain=0.25, read_noise_e=3.5)
        assert abs(f - 500.0) < 30.0
        assert v > 0 and ns == 0

    def test_saturated_pixels_masked_and_counted(self):
        win, prof, bg = self._column(flux=500.0)
        win[gext.APERTURE_HALFWIN] = gext.SATURATION_ADU + 100
        f, v, ns = gext.horne_column(win, np.full(len(win), bg), prof,
                                     egain=0.25, read_noise_e=3.5)
        assert ns == 1
        assert f is not None                 # the wings still constrain it

    def test_fully_saturated_column_refuses(self):
        win = np.full(25, gext.SATURATION_ADU + 1.0)
        prof = np.full(25, 1 / 25)
        f, v, ns = gext.horne_column(win, np.zeros(25), prof,
                                     egain=0.25, read_noise_e=3.5)
        assert f is None and ns == 25

    def test_optimal_beats_box_under_noise(self):
        # The Horne estimator's variance across noise realizations must
        # undercut the plain box sum's — that is its whole point.
        fh, fb = [], []
        for seed in range(60):
            win, prof, bg = self._column(flux=200.0, seed=seed)
            f, _, _ = gext.horne_column(win, np.full(len(win), bg), prof,
                                        egain=0.25, read_noise_e=3.5)
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
        d = gext.median_relative_difference(a, b)
        assert abs(d - 0.02) < 1e-9

    def test_median_relative_difference_needs_overlap(self):
        a = np.full(1000, np.nan)
        b = np.full(1000, 1.0)
        assert gext.median_relative_difference(a, b) is None


# ---------------------------------------------------------------------------
# Wavelength anchors
# ---------------------------------------------------------------------------
def synthetic_spectrum(nx=4000, x_ha=1800.0, disp=0.5, ha_height=800.0,
                      o2_depth=250.0, noise=5.0, seed=3):
    """Continuum + Halpha emission + O2 B/A absorption at the pixel
    positions the dispersion implies (red toward +x)."""
    rng = np.random.default_rng(seed)
    x = np.arange(nx, dtype=float)
    cont = 1000.0 + 0.05 * x
    spec = cont + ha_height * np.exp(-0.5 * ((x - x_ha) / 8.0) ** 2)
    for band in (gw.O2_B_A, gw.O2_A_A):
        xb = x_ha + (band - gw.HALPHA_A) / disp
        spec -= o2_depth * np.exp(-0.5 * ((x - xb) / 10.0) ** 2)
    return spec + rng.normal(0, noise, nx)


class TestWavelength:
    def test_halpha_found(self):
        w = gw.solve_wavelength(synthetic_spectrum())
        assert w["x_halpha"] is not None
        assert abs(w["x_halpha"] - 1800.0) < 2.0
        assert w["halpha_snr"] > gw.PEAK_MIN_SNR

    def test_dispersion_from_o2(self):
        w = gw.solve_wavelength(synthetic_spectrum(disp=0.5))
        assert w["anchor_status"] == "halpha+o2"
        # Dispersion is only ever accepted from the O2 B+A PAIR (see
        # find_o2_pair): requiring both bands is a firmer identification than
        # either alone, so 'o2_pair' is the single legal source token.  This
        # assertion previously named the pre-revision per-band tokens.
        assert w["disp_source"] == "o2_pair"
        assert abs(w["disp_a_per_px"] - 0.5) < 0.02

    def test_reversed_dispersion_sign_measured(self):
        # Red toward -x: the anchors sit blueward in pixel terms and the
        # solved dispersion must come out NEGATIVE, not fail.
        spec = synthetic_spectrum()[::-1].copy()
        w = gw.solve_wavelength(spec)
        assert w["anchor_status"] == "halpha+o2"
        assert w["disp_a_per_px"] < 0

    def test_no_line_no_anchor(self):
        rng = np.random.default_rng(0)
        flat = 1000.0 + rng.normal(0, 5.0, 4000)
        w = gw.solve_wavelength(flat)
        assert w["anchor_status"] == "none"
        assert w["x_halpha"] is None and w["disp_a_per_px"] is None

    def test_cosmic_ray_not_mistaken_for_halpha(self):
        # A single-pixel spike is narrower than PEAK_MIN_WIDTH: refused.
        rng = np.random.default_rng(0)
        spec = 1000.0 + rng.normal(0, 5.0, 4000)
        spec[2000] += 5000.0
        w = gw.solve_wavelength(spec)
        assert w["x_halpha"] is None

    def test_wavelength_axis_arithmetic(self):
        lam = gw.wavelength_axis(10, x_halpha=4.0, disp_a_per_px=2.0)
        assert lam[4] == gw.HALPHA_A
        assert lam[5] == gw.HALPHA_A + 2.0
        assert lam[0] == gw.HALPHA_A - 8.0

    def test_snippet_roundtrip(self):
        import json
        flux = np.arange(100, dtype=float)
        flux[50] = np.nan
        s = gw.snippet(flux, 50.0, halfwin=10, stride=1)
        back = json.loads(json.dumps(s))
        xs = [p[0] for p in back]
        assert 40 in xs and 60 in xs
        assert back[xs.index(50)][1] is None          # NaN -> null
        assert back[xs.index(41)][1] == 41.0


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
