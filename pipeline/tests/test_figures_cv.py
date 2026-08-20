"""CV-S11 tests: the pure arithmetic behind the manuscript's figures and
its ``numbers.tex``.

WHAT IS TESTED HERE, AND WHY THESE THINGS
------------------------------------------
Two kinds of failure would put a false statement into a published paper,
and both have actually happened during this stage's development, so both
are pinned:

1. **A pure function that quietly does the wrong arithmetic.**  Folding,
   phase binning, quasi-simultaneous pairing and robust limits are the
   operations every figure rests on.  A fold that is off by a constant
   still looks like a fold.
2. **A formatting or lookup slip that silently deletes a result.**  The
   state palette is keyed in lower case and the classifier writes upper
   case; when that mismatch existed, every accretion state in Figures 5,
   7 and 8 drew in the unknown-state grey, and the figures still looked
   fine.  Likewise a LaTeX macro name containing a digit is rejected by
   TeX, and a caption carrying a bare underscore ends the build.

The tests below use synthetic arrays, not the products database: they are
about whether the arithmetic is right, not about what this season's data
happen to say.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_phot import figures_cv as fx        # noqa: E402
from macro_phot import numbers_cv as nx        # noqa: E402


# ===========================================================================
# fold_phase
# ===========================================================================
def test_fold_phase_puts_the_epoch_at_zero():
    """A time equal to the epoch folds to phase zero, not to phase one."""
    assert fx.fold_phase([2460000.0], 0.1, 2460000.0)[0] == pytest.approx(0.0)


def test_fold_phase_is_periodic_and_in_range():
    t = 2460000.0 + np.array([0.0, 0.025, 0.05, 0.075, 0.1, 0.125, 1.7])
    ph = fx.fold_phase(t, 0.1, 2460000.0)
    assert np.all((ph >= 0.0) & (ph < 1.0))
    # One whole period later must land on the same phase.
    assert ph[4] == pytest.approx(ph[0], abs=1e-9)
    assert ph[5] == pytest.approx(ph[1], abs=1e-9)


def test_fold_phase_handles_times_before_the_epoch():
    """Negative elapsed time must wrap forward, not produce a negative phase.

    Several targets here were observed before their catalogue epoch, so
    this is the normal case and not an edge case.
    """
    ph = fx.fold_phase([2459999.975], 0.1, 2460000.0)
    assert 0.0 <= ph[0] < 1.0
    assert ph[0] == pytest.approx(0.75, abs=1e-9)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), None])
def test_fold_phase_refuses_an_impossible_period(bad):
    with pytest.raises(ValueError):
        fx.fold_phase([2460000.0], bad, 2460000.0)


# ===========================================================================
# mad_sigma and robust_ylim
# ===========================================================================
def test_mad_sigma_recovers_a_gaussian_width():
    rng = np.random.default_rng(11)
    x = rng.normal(5.0, 0.2, 20_000)
    assert fx.mad_sigma(x) == pytest.approx(0.2, rel=0.05)


def test_mad_sigma_ignores_a_wild_outlier():
    """The reason this is used everywhere instead of an RMS."""
    x = np.concatenate([np.full(200, 1.0), np.full(200, 1.2), [1e6]])
    assert fx.mad_sigma(x) < 1.0
    assert np.std(x) > 1e4


def test_mad_sigma_on_empty_is_nan_not_zero():
    """Zero scatter is a claim; no measurement is not."""
    assert math.isnan(fx.mad_sigma([]))


def test_robust_ylim_is_not_dragged_by_one_bad_point():
    y = np.concatenate([np.linspace(16.0, 16.5, 500), [20.0]])
    lo, hi = fx.robust_ylim(y)
    assert hi < 17.5, "one 20th-magnitude frame must not set the axis"
    assert lo < 16.0 and hi > 16.5, "the real data must still fit"


def test_robust_ylim_survives_degenerate_input():
    lo, hi = fx.robust_ylim([3.0, 3.0, 3.0])
    assert hi > lo


# ===========================================================================
# phase_bin
# ===========================================================================
def test_phase_bin_recovers_a_known_sinusoid():
    rng = np.random.default_rng(3)
    ph = rng.uniform(0, 1, 6000)
    y = 0.3 * np.sin(2 * np.pi * ph)
    c, m, s, n = fx.phase_bin(ph, y, n_bins=20)
    assert n.sum() == 6000
    good = np.isfinite(m)
    assert good.all()
    np.testing.assert_allclose(m[good], 0.3 * np.sin(2 * np.pi * c[good]),
                               atol=0.02)


def test_phase_bin_leaves_empty_bins_as_nan():
    """An empty phase bin is a statement about coverage.

    Filling it in would let a reader mistake sampling for shape, which is
    the single most common way a folded light curve lies.
    """
    ph = np.linspace(0.0, 0.45, 200)          # nothing above phase 0.5
    c, m, s, n = fx.phase_bin(ph, np.ones_like(ph), n_bins=10)
    assert np.isnan(m[-1])
    assert n[-1] == 0


def test_phase_bin_respects_min_count():
    ph = np.array([0.05, 0.55, 0.56, 0.57])
    c, m, s, n = fx.phase_bin(ph, np.ones(4), n_bins=2, min_count=3)
    assert n[0] == 1 and n[1] == 3
    assert np.isnan(m[0]), "a one-point bin is not a median"
    assert not np.isnan(m[1]), "a three-point bin meets the bar"


# ===========================================================================
# pair_quasi_simultaneous
# ===========================================================================
def test_pairing_matches_the_nearest_point_and_subtracts():
    ta = np.array([0.0, 1.0])
    tb = np.array([0.0, 1.0])
    t, colour, dt = fx.pair_quasi_simultaneous(
        ta, np.array([10.0, 11.0]), tb, np.array([9.0, 9.5]),
        max_dt_s=1e9)
    np.testing.assert_allclose(colour, [1.0, 1.5])
    np.testing.assert_allclose(dt, [0.0, 0.0])


def test_pairing_rejects_pairs_wider_than_the_gate():
    """The gate is the whole point: a pair half an orbit apart is not a
    colour, it is two different states subtracted."""
    ta = np.array([0.0])
    tb = np.array([0.02])                      # 1728 s away
    t, colour, dt = fx.pair_quasi_simultaneous(
        ta, np.array([10.0]), tb, np.array([9.0]), max_dt_s=600.0)
    assert t.size == 0


def test_pairing_on_empty_input_returns_empty():
    t, c, d = fx.pair_quasi_simultaneous([], [], [1.0], [1.0])
    assert t.size == c.size == d.size == 0


def test_pairing_is_order_independent_in_b():
    """Band B need not arrive sorted; the matcher sorts it itself."""
    ta = np.array([0.0, 0.001])
    tb = np.array([0.001, 0.0])
    mb = np.array([2.0, 1.0])
    _, colour, _ = fx.pair_quasi_simultaneous(
        ta, np.array([10.0, 20.0]), tb, mb, max_dt_s=1e9)
    np.testing.assert_allclose(colour, [9.0, 18.0])


# ===========================================================================
# state normalisation -- the bug that silently greyed out three figures
# ===========================================================================
@pytest.mark.parametrize("raw,expect", [
    ("HIGH", "high"), ("high", "high"), ("Low", "low"),
    ("INTERMEDIATE", "intermediate"), ("UNCLASSIFIED", "unclassified"),
    ("UNKNOWN", "unknown"), (None, "unknown"), ("", "unknown"),
    ("something new", "unknown"),
])
def test_normalise_state(raw, expect):
    assert fx.normalise_state(raw) == expect


def test_every_normalised_state_has_a_colour():
    """The invariant the greyed-out-figures bug violated."""
    for raw in ("HIGH", "LOW", "INTERMEDIATE", "UNCLASSIFIED", "UNKNOWN",
                None, "anything at all"):
        assert fx.normalise_state(raw) in fx.STATE_COLOR


# ===========================================================================
# date axis
# ===========================================================================
def test_night_to_ordinal_is_zero_at_the_epoch():
    assert fx.night_to_ordinal(["2024-01-01"])[0] == 0.0
    assert fx.night_to_ordinal(["2024-01-02"])[0] == 1.0


def test_night_to_ordinal_returns_nan_for_junk():
    out = fx.night_to_ordinal(["not-a-date", None, "2024-06-01"])
    assert math.isnan(out[0]) and math.isnan(out[1])
    assert out[2] > 0


@pytest.mark.parametrize("span_days", [40, 200, 800, 4000, 20000])
def test_year_ticks_never_floods_the_axis(span_days):
    """The bug this pins produced a hundred overlapping labels on a
    forty-year AAVSO baseline."""
    pos, lab = fx.year_ticks(0.0, float(span_days), max_ticks=10)
    assert len(pos) == len(lab)
    assert len(pos) <= 12, f"{len(pos)} ticks over {span_days} d"
    assert all(0.0 - 1 <= p <= span_days + 1 for p in pos)


# ===========================================================================
# series keys and labels
# ===========================================================================
def test_series_parts_and_label():
    assert fx.series_parts("stlmi|e76|g") == ("stlmi", 76, "g")
    assert fx.series_label("stlmi|e76|g") == "ST LMi Mode0 g"


def test_every_filter_has_a_colour_and_a_marker():
    for f in ("G", "g", "R", "r", "I", "i", "z", "y"):
        assert f in fx.BAND_COLOR
        assert f in fx.BAND_MARKER


def test_matching_bands_share_a_colour_across_eras():
    """G and g are the same part of the spectrum through two bandpasses.

    Giving them different colours would let a reader read an instrument
    change as a physical one; the bandpass difference is stated in the
    caption instead.
    """
    assert fx.BAND_COLOR["G"] == fx.BAND_COLOR["g"]
    assert fx.BAND_COLOR["R"] == fx.BAND_COLOR["r"]
    assert fx.BAND_COLOR["I"] == fx.BAND_COLOR["i"]


def test_every_readout_mode_has_a_short_label():
    for mode in fx.MODE_MARKER:
        assert mode in fx.MODE_SHORT


# ===========================================================================
# FigureSpec captions -- these strings are pasted into LaTeX
# ===========================================================================
def test_caption_escapes_table_names():
    """A bare ``cv_frames`` in a \\caption is a subscript outside maths
    mode and ends the tectonic run."""
    spec = fx.FigureSpec(fig_id="figXX", label="fig:x", title="t",
                         caption="Body.", tables=("cv_frames", "p3_period"))
    cap = spec.full_caption
    assert "\\texttt{cv\\_frames}" in cap
    assert "cv_frames" not in cap.replace("cv\\_frames", "")


def test_substitute_caption_states_the_substitution():
    spec = fx.FigureSpec(fig_id="figXX", label="fig:x", title="t",
                         caption="Body.", substitute=True,
                         substitute_reason="the nights do not exist")
    assert "SUBSTITUTE" in spec.full_caption
    assert "the nights do not exist." in spec.full_caption


def test_every_registered_builder_is_callable_and_declares_its_databases():
    assert len(fx.FIGURE_IDS) == 13, "the strategy specifies thirteen"
    for fig_id, entry in fx.BUILDERS.items():
        assert callable(entry["fn"])
        assert entry["needs"], f"{fig_id} declares no database"
        assert set(entry["needs"]) <= {"cv", "ch", "man"}


def test_aastex_widths_are_the_only_two_allowed():
    assert fx.COL_SINGLE == 3.5
    assert fx.COL_DOUBLE == 7.1


# ===========================================================================
# numbers_cv: macro names
# ===========================================================================
def test_macro_name_is_letters_only():
    """LaTeX control sequences may not contain digits."""
    assert nx.tex_macro_name("stlmi_e7_nights") == "NumStlmiESevenNights"
    assert nx.tex_macro_name("fig01", prefix="Cap") == "CapFigZeroOne"
    for key in ("a1", "2024 frames", "st lmi", "p4_run-scopes"):
        assert nx.tex_macro_name(key).isalpha()


def test_macro_names_are_distinct_for_distinct_keys():
    keys = ["st lmi nights", "st lmi frames", "vv pup nights",
            "sigma t ideal s", "sigma t injected median s"]
    names = [nx.tex_macro_name(k) for k in keys]
    assert len(set(names)) == len(names)


# ===========================================================================
# numbers_cv: formatting
# ===========================================================================
def test_fmt_int_uses_a_latex_thin_space():
    assert nx.fmt_int(8716) == "8\\,716"
    assert nx.fmt_int(3174618) == "3\\,174\\,618"
    assert nx.fmt_int(42) == "42"
    assert nx.fmt_int(-1500) == "-1\\,500"


def test_fmt_int_on_nothing_is_none():
    assert nx.fmt_int(None) is None
    assert nx.fmt_int("not a number") is None


def test_fmt_float_rejects_nan_and_inf():
    """A NaN printed into a paper as 'nan' is a typesetting bug; a NaN
    that becomes the missing-value marker is a visible warning."""
    assert nx.fmt_float(float("nan")) is None
    assert nx.fmt_float(float("inf")) is None
    assert nx.fmt_float(1.239, 2) == "1.24"


def test_fmt_range_and_sci():
    assert nx.fmt_range(9, 77, 0) == "9--77"
    assert nx.fmt_range(0.6, 1.77, 2) == "0.60--1.77"
    assert nx.fmt_range(None, 5) is None
    assert nx.fmt_sci(5e-9, 1) == "5.0 \\times 10^{-9}"


def test_fmt_percent():
    assert nx.fmt_percent(0.266, 1) == "26.6"
    assert nx.fmt_percent(None) is None


# ===========================================================================
# numbers_cv: rendering
# ===========================================================================
def test_render_defines_the_missing_marker_first():
    tex = nx.render_tex([nx.Number("a b", "1")])
    assert "\\providecommand{\\NumMissing}" in tex
    assert tex.index("NumMissing") < tex.index("NumAB")


def test_render_emits_the_marker_for_an_unmeasured_value():
    """A macro that silently vanished would let a sentence lose its number
    and still compile."""
    tex = nx.render_tex([nx.Number("missing thing", None)])
    assert "\\newcommand{\\NumMissingThing}{\\NumMissing}" in tex


def test_render_refuses_duplicate_macros():
    with pytest.raises(ValueError):
        nx.render_tex([nx.Number("a b", "1"), nx.Number("a  b", "2")])


def test_render_carries_the_unit_and_source_as_a_comment():
    tex = nx.render_tex([nx.Number("hg ceiling adu", "3\\,496", "ADU",
                                   "detector_params")])
    line = [l for l in tex.splitlines() if "HgCeilingAdu" in l][0]
    assert line.split("%")[0].strip().endswith("{3\\,496}")
    assert "ADU" in line and "detector_params" in line


def test_number_body_falls_back_to_the_marker():
    assert nx.Number("k", None).body == "\\NumMissing"
    assert nx.Number("k", "7").body == "7"


# ===========================================================================
# numbers_cv: LaTeX escaping in table cells
# ===========================================================================
@pytest.mark.parametrize("raw,expect", [
    ("cv_series", "cv\\_series"),
    ("50%", "50\\%"),
    ("a & b", "a \\& b"),
    ("#1", "\\#1"),
    ("", ""),
    (None, ""),
])
def test_table_cell_escaping(raw, expect):
    assert nx._esc(raw) == expect


def test_table_cell_escaping_survives_a_real_verdict_string():
    """Verdict strings from p4_verdict carry underscores and per-cent signs
    and go straight into a table cell."""
    raw = "0 of 6 scopes clear the 90% contour; dP_sh/dt not measurable"
    out = nx._esc(raw)
    assert "\\%" in out and "\\_" in out
    assert "%" not in out.replace("\\%", "")
    assert "_" not in out.replace("\\_", "")
