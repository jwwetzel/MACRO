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
import re
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


# ===========================================================================
# pdot_envelope_seconds -- the curve Figure 9(a) draws
# ===========================================================================
def test_envelope_removes_a_constant_and_a_linear_term():
    """The envelope is the quadratic AFTER an ephemeris fit absorbs what it
    can, so a weighted straight line through it must be flat.  Figure 9 used
    to draw a bare re-centred parabola and call it this."""
    e = np.linspace(0.0, 1000.0, 40)
    s = np.full_like(e, 10.0)
    env = fx.pdot_envelope_seconds(e, e, s, 1e-8, 0.08)
    w = 1.0 / s ** 2
    design = np.vstack([np.ones_like(e), e]).T
    beta = np.linalg.solve(design.T @ (design * w[:, None]),
                           design.T @ (w * env))
    assert abs(beta[0]) < 1e-6 and abs(beta[1]) < 1e-9, (
        "a weighted constant-plus-linear fit to the envelope is not zero, "
        "so the terms the ephemeris fit absorbs have not been removed")


def test_envelope_still_carries_the_curvature():
    """Removing the linear terms must not remove the signal: the second
    difference is what a Pdot puts in, and it must survive."""
    e = np.linspace(0.0, 1000.0, 41)
    s = np.full_like(e, 10.0)
    env = fx.pdot_envelope_seconds(e, e, s, 1e-8, 0.08)
    expect = 0.5 * 1e-8 * 0.08 * 86400.0          # a2, seconds per cycle^2
    step = e[1] - e[0]
    second = np.diff(env, 2) / step ** 2
    assert np.allclose(second, 2.0 * expect, rtol=1e-6)


def test_envelope_scales_linearly_with_the_bound():
    e = np.linspace(0.0, 500.0, 20)
    s = np.full_like(e, 5.0)
    a = fx.pdot_envelope_seconds(e, e, s, 1e-9, 0.08)
    b = fx.pdot_envelope_seconds(e, e, s, 3e-9, 0.08)
    assert np.allclose(b, 3.0 * a, rtol=1e-9)


def test_envelope_falls_back_when_there_is_nothing_to_project_out():
    """Fewer than three usable epochs cannot define a projection; the
    function must return a curve rather than raise inside a figure build."""
    e = np.array([0.0, 10.0])
    out = fx.pdot_envelope_seconds(e, e, np.array([np.nan, np.nan]),
                                   1e-9, 0.08)
    assert np.all(np.isfinite(out))


def test_envelope_report_counts_the_epochs_outside():
    """The caption generator may only state containment when this is zero,
    and Figure 9's caption asserted containment when it was 27."""
    rep = fx._envelope_report([10.0, -300.0, 5.0], [8.0, 8.0, 8.0],
                              [50.0, 50.0, 50.0])
    assert rep["n"] == 3
    assert rep["n_outside"] == 1
    assert rep["n_under_error"] == 0
    assert rep["env_max"] == pytest.approx(50.0)


def test_envelope_report_gives_the_two_ends_separately():
    """Referee 4, minor.  Figure 9's caption said the envelope "reaches only
    293 s at the ends of the baseline"; it reaches 293 s at the late end and
    49 s at the early one, because the curve is a parabola with a weighted
    linear fit removed and is strongly asymmetric.  The report has to be
    able to tell a caption the two ends apart, and has to order them by
    cycle rather than by the order the query happened to return."""
    rep = fx._envelope_report([0.0, 0.0, 0.0], [8.0, 8.0, 8.0],
                              [293.0, -76.0, 49.0],
                              cycles=[21868.5, 16000.0, 13181.0])
    assert rep["env_first"] == pytest.approx(49.0)
    assert rep["env_last"] == pytest.approx(293.0)
    assert rep["env_max"] == pytest.approx(293.0)
    # Without cycles it must still answer, in array order.
    plain = fx._envelope_report([0.0], [8.0], [49.0])
    assert plain["env_first"] == plain["env_last"] == pytest.approx(49.0)


# ===========================================================================
# A scope's name carries every night it folds
# ===========================================================================
def test_a_two_night_block_is_not_labelled_by_one_night():
    """Referee 4, minor.  Three of Figure 11's six rows are two-night blocks
    and were labelled '2024-05-02', which is one of their nights."""
    from macro_phot import final_science as fs
    assert fs.run_scope_label("2024-05-02+2024-05-03", None) == \
        "2024-05-02+03"
    assert fs.run_scope_label("2024-02-21", "2024-02-20") == "2024-02-21"
    assert fs.is_multi_night_scope("2024-05-02+2024-05-03", None)
    assert not fs.is_multi_night_scope("2024-02-21", "2024-02-20")


def test_a_block_across_a_month_boundary_keeps_the_second_month():
    """The day alone is ambiguous when the two nights are in different
    months, so the label keeps the whole date there."""
    from macro_phot import final_science as fs
    assert fs.run_scope_label("2024-04-30+2024-05-01", None) == \
        "2024-04-30+2024-05-01"


# ===========================================================================
# dynamic_range_ratios -- one comparison, one pair of numbers
# ===========================================================================
def test_dynamic_range_ratios_span_the_sixteen_bit_modes():
    """§2.1 said 'nearly twenty' and Figure 3's caption 'a factor of 16'
    for the same comparison; both come from here now."""
    params = {("High Gain", "ceiling_adu"): 3496.0,
              ("High Gain", "adc_bits"): 12.0,
              ("Mode0", "ceiling_adu"): 65535.0,
              ("Mode0", "adc_bits"): 16.0,
              ("High Gain StackPro", "ceiling_adu"): 56062.0,
              ("High Gain StackPro", "adc_bits"): 16.0}
    lo, hi = nx.dynamic_range_ratios(params)
    assert lo == pytest.approx(56062.0 / 3496.0, rel=1e-9)
    assert hi == pytest.approx(65535.0 / 3496.0, rel=1e-9)
    assert lo < hi, "the modes span a range and the paper must quote one"


def test_dynamic_range_ratios_ignore_twelve_bit_modes():
    params = {("High Gain", "ceiling_adu"): 3496.0,
              ("High Gain", "adc_bits"): 12.0,
              ("Low Gain", "ceiling_adu"): 3000.0,
              ("Low Gain", "adc_bits"): 12.0,
              ("Mode0", "ceiling_adu"): 65535.0,
              ("Mode0", "adc_bits"): 16.0}
    lo, hi = nx.dynamic_range_ratios(params)
    assert lo == hi == pytest.approx(65535.0 / 3496.0, rel=1e-9)


def test_dynamic_range_ratios_on_missing_input():
    assert nx.dynamic_range_ratios({}) == (None, None)
    assert nx.dynamic_range_ratios(
        {("High Gain", "ceiling_adu"): 3496.0}) == (None, None)


# ===========================================================================
# The scope clause and the panel list, shared between renderers
# ===========================================================================
def test_the_precision_scope_clause_retracts_the_hold_out_rule():
    """Figure 2's caption and §3.1 read this one string; it must not
    reintroduce the rule §3.1 exists to retract."""
    clause = nx.PRECISION_SCOPE_CLAUSE
    assert "not held-out statistics" in clause
    assert "held-out check stars of a tied solve" not in clause
    # It is pasted straight into a \caption, so a subscript OUTSIDE maths
    # mode would end the tectonic run.  Strip the maths spans, then look.
    outside_maths = re.sub(r"\$[^$]*\$", "", clause).replace("\\_", "")
    assert "_" not in outside_maths, (
        "the scope clause carries a subscript outside maths mode and would "
        "break the build when pasted into a caption")
    assert clause.count("$") % 2 == 0
    assert clause.count("{") == clause.count("}")


def test_the_colour_panel_pairs_are_the_four_figure_six_draws():
    assert len(nx.COLOUR_PANEL_PAIRS) == 4
    assert all(len(p) == 3 for p in nx.COLOUR_PANEL_PAIRS)
    assert {p[0] for p in nx.COLOUR_PANEL_PAIRS} == {7, 76}


def test_fmt_range_collapses_a_degenerate_range():
    """Table 2's caption printed '4--4 in every block', which reads as a
    spread these data do not have."""
    assert nx.fmt_range(4, 4, 0) == "4"
    assert nx.fmt_range(4, 5, 0) == "4--5"
    assert nx.fmt_range(0.10, 0.104, 2) == "0.10", (
        "two values that round to the same printed number are the same "
        "printed number")


# ===========================================================================
# numbers_cv: which released database a macro came from  (referee 5)
# ===========================================================================
def _fake_release():
    """Three in-memory stand-ins for the three released databases."""
    import sqlite3
    cv = sqlite3.connect(":memory:")
    cv.execute("CREATE TABLE cv_series (a)")
    ch = sqlite3.connect(":memory:")
    ch.execute("CREATE TABLE ch_contour (a)")
    man = sqlite3.connect(":memory:")
    man.execute("CREATE TABLE detector_params (a)")
    return cv, ch, man


def test_resolve_databases_records_the_file_each_table_lives_in():
    """§7 names three databases and says which macros come from which; the
    mapping is measured here rather than asserted in prose."""
    cv, ch, man = _fake_release()
    out = nx.resolve_databases(
        [nx.Number("a", "1", source="cv_series", note="n"),
         nx.Number("b", "2", source="ch_contour", note="n"),
         nx.Number("c", "3", source="detector_params", note="n")],
        cv, ch, man)
    assert [n.db for n in out] == ["cv_timeseries.sqlite",
                                   "cv_characterization.sqlite",
                                   "rlmt-manifest.sqlite"]
    assert {n.kind for n in out} == {"measured"}


def test_an_external_constant_resolves_to_no_database_and_is_flagged():
    cv, ch, man = _fake_release()
    out = nx.resolve_databases(
        [nx.Number("a", "60", source="ANALYSIS_STRATEGY §4", note="n"),
         nx.Number("b", "0.079", source="VSX catalogue", note="n"),
         nx.Number("c", "50", source="literature: superhumps", note="n"),
         nx.Number("d", "12", source="CV-S10 constant", note="n")],
        cv, ch, man)
    assert {n.kind for n in out} == {"external"}
    assert {n.db for n in out} == {""}


def test_a_source_that_names_nothing_stops_the_build():
    """The failure mode this guard exists for: a mistyped table name would
    otherwise be silently reclassified as 'not a measurement'."""
    cv, ch, man = _fake_release()
    with pytest.raises(ValueError, match="names no table"):
        nx.resolve_databases(
            [nx.Number("a", "1", source="cv_seires", note="n")],
            cv, ch, man)


def test_a_table_in_two_databases_is_refused_as_ambiguous():
    cv, ch, man = _fake_release()
    ch.execute("CREATE TABLE cv_series (a)")
    with pytest.raises(ValueError, match="ambiguous"):
        nx.resolve_databases(
            [nx.Number("a", "1", source="cv_series", note="n")],
            cv, ch, man)


def test_a_macro_with_no_note_is_refused():
    """72 macros once recorded a table name and nothing else, in a paper
    whose §7 promises a note for every value."""
    cv, ch, man = _fake_release()
    with pytest.raises(ValueError, match="carries no note"):
        nx.resolve_databases(
            [nx.Number("a", "1", source="cv_series", note="   ")],
            cv, ch, man)


def test_the_rendered_macro_line_names_the_database_too():
    n = nx.Number("hg ceiling adu", "3\\,496", "ADU", "detector_params",
                  "note", "rlmt-manifest.sqlite", "measured")
    line = [l for l in nx.render_tex([n]).splitlines()
            if "HgCeilingAdu" in l][0]
    assert "detector_params in rlmt-manifest.sqlite" in line
    ext = nx.Number("sigma t threshold s", "60", "s",
                    "ANALYSIS_STRATEGY §4", "note", "", "external")
    line = [l for l in nx.render_tex([ext]).splitlines()
            if "SigmaTThresholdS" in l][0]
    assert "external constant" in line
