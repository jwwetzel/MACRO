"""The manuscript's own invariants, checked against the emitted files.

``test_figures_cv.py`` protects the pure functions and ``test_cv_products.py``
the products.  Neither can catch the class of defect three rounds of review
have now found in this paper: a macro that is arithmetically fine but is
quoted in a sentence it does not support, a caption that survives a
retraction made in the body, or two generators that name the same thing two
ways.  Those live in the seam between the emitters and ``main.tex``, and
this file is where that seam is tested.

Each test names the printed claim that was wrong.  All of them skip when the
manuscript has not been built, because ``manuscripts/`` is not in git.
"""

from __future__ import annotations

import math
import re
import sqlite3
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MANUSCRIPT = REPO_ROOT / "manuscripts" / "CV_TimeSeries"
PHOT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"


def _text(name: str) -> str:
    p = MANUSCRIPT / name
    if not p.exists():
        pytest.skip(f"{name} not emitted in this checkout")
    return p.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def numbers() -> dict:
    """``{macro name without backslash: body}`` from ``numbers.tex``."""
    out = {}
    for m in re.finditer(r"\\newcommand\{\\(Num[A-Za-z]+)\}\{(.*)\}\s*(?:%.*)?$",
                         _text("numbers.tex"), re.M):
        out[m.group(1)] = m.group(2)
    assert out, "numbers.tex parsed to nothing"
    return out


@pytest.fixture(scope="module")
def captions() -> dict:
    out = {}
    for m in re.finditer(r"\\newcommand\{\\(Cap[A-Za-z]+)\}\{(.*)\}$",
                         _text("captions.tex"), re.M):
        out[m.group(1)] = m.group(2)
    assert out, "captions.tex parsed to nothing"
    return out


@pytest.fixture(scope="module")
def body() -> str:
    return _text("main.tex")


@pytest.fixture(scope="module")
def phot():
    if not PHOT_DB.exists():
        pytest.skip("cv_timeseries.sqlite not built in this checkout")
    con = sqlite3.connect(f"file:{PHOT_DB}?mode=ro", uri=True, timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    yield con
    con.close()


def _num(numbers: dict, name: str) -> float:
    """A macro body as a float, undoing this module's LaTeX formatting."""
    raw = numbers[name]
    raw = raw.replace("\\,", "").replace("$", "").strip()
    m = re.fullmatch(r"([-0-9.]+)\s*\\times\s*10\^\{(-?\d+)\}", raw)
    if m:
        return float(m.group(1)) * 10.0 ** int(m.group(2))
    return float(raw)


# ---------------------------------------------------------------------------
# The headline bound, in its three forms
# ---------------------------------------------------------------------------
class TestThePeriodChangeBoundIsSelfConsistent:
    """Referee 3's blocker.  The null was printed three ways --- a
    dimensionless |Pdot|, a rate in s/yr, and a timescale P/|Pdot| --- and the
    third disagreed with the other two by a factor 12.6, because the emitter
    computed 1/|Pdot| and dropped the period from the numerator.  Any reader
    with the abstract alone could do the arithmetic.  These make the three
    forms unable to disagree again."""

    def test_timescale_times_rate_is_the_orbital_period(self, numbers, phot):
        """P/|Pdot| in yr multiplied by |Pdot| in s/yr must return P in s.
        This is the identity the published numbers violated."""
        timescale_yr = _num(numbers, "NumStLmiPdotTimescaleYr")
        rate_s_per_yr = _num(numbers, "NumStLmiPdotLimitSPerYr")
        period_d = phot.execute("SELECT period_d FROM p3_cycle_count WHERE "
                                "target_key='stlmi'").fetchone()[0]
        period_s = float(period_d) * 86400.0
        # The tolerance is set by the PRINTED precision, not by the
        # arithmetic: the rate is published to two decimals (0.11 from
        # 0.1140, a half-ulp of 4.4 per cent) and the timescale to one
        # significant figure (0.8 per cent).  Eight per cent covers both
        # with room, and still catches the factor-12.6 error this test
        # exists for by two orders of magnitude.
        assert timescale_yr * rate_s_per_yr == pytest.approx(period_s,
                                                             rel=0.08), (
            f"P/|Pdot| = {timescale_yr:.3g} yr times "
            f"{rate_s_per_yr:.3g} s/yr is "
            f"{timescale_yr * rate_s_per_yr:.4g} s, but the orbital period "
            f"is {period_s:.1f} s: the three forms of the bound disagree")

    def test_the_rate_is_the_dimensionless_bound_in_seconds_per_year(
            self, numbers):
        """The s/yr figure must be the dimensionless bound times a year."""
        dimensionless = _num(numbers, "NumStLmiPdotLimit")
        rate = _num(numbers, "NumStLmiPdotLimitSPerYr")
        period_s = 6833.3
        assert rate == pytest.approx(dimensionless * 365.25 * 86400.0,
                                     rel=0.06)
        # And the timescale follows from those two alone.
        assert period_s / rate == pytest.approx(
            _num(numbers, "NumStLmiPdotTimescaleYr"), rel=0.05)

    def test_the_paper_says_which_three_sigma_convention_it_uses(
            self, body, numbers):
        """``pdot_limit3`` is |Pdot_fit| + 3 sigma, not 3 sigma; the two
        differ by half again here and a reader cannot tell which is meant
        unless the paper says."""
        assert "NumStLmiPdotThreeSigma" in numbers
        assert "\\dot{P}_{\\rm fit}| + 3\\sigma" in body, (
            "§5.1 does not name the convention behind its '3-sigma bound'")


# ---------------------------------------------------------------------------
# The constant subtracted from the O-C
# ---------------------------------------------------------------------------
class TestTheOCOffsetIsDisclosed:
    """Referee 3, major.  CV-S9 subtracts the mean per-cycle O-C --- 1,071 s,
    0.157 cycles --- before writing p3_oc.  The subtraction is right; it
    reached no word of the manuscript, which said in three places that the
    residuals were 'against the catalogue ephemeris', and a reader took that
    to mean the edge falls at the catalogue epoch's phase."""

    def test_the_offset_is_emitted_as_a_macro(self, numbers, phot):
        stored = phot.execute("SELECT oc_mean_s, period_d FROM p3_cycle_count "
                              "WHERE target_key='stlmi'").fetchone()
        assert _num(numbers, "NumStLmiOcOffsetS") == pytest.approx(
            float(stored[0]), rel=0.01)
        assert _num(numbers, "NumStLmiOcOffsetCycles") == pytest.approx(
            float(stored[0]) / (float(stored[1]) * 86400.0), abs=0.001)

    def test_the_body_states_what_the_residuals_are_measured_against(
            self, body):
        """The disclosure must be in §4.3, where the O-C is defined, and not
        only in a caption."""
        assert "\\NumStLmiOcOffsetS" in body
        assert "What the residuals are measured against" in body

    def test_figure_nine_does_not_claim_the_residuals_are_against_the_epoch(
            self, captions):
        cap = captions["CapFigZeroNine"]
        assert "residuals against the catalogue PERIOD" in cap
        assert "does not fall at the catalogue epoch's phase zero" in cap

    def test_the_reduced_chisq_loses_the_absorbed_constant(self, numbers,
                                                           phot):
        """One parameter was estimated from these edges, so nu = N - 1.
        The paper quoted chi-squared per epoch and called it reduced."""
        rows = phot.execute("SELECT oc_s, oc_sigma_s FROM p3_oc_night "
                            "WHERE target_key='stlmi'").fetchall()
        chi2 = sum((o / s) ** 2 for o, s in rows)
        dof = len(rows) - 1
        assert int(_num(numbers, "NumStLmiOcDof")) == dof
        assert _num(numbers, "NumStLmiOcChisq") == pytest.approx(chi2 / dof,
                                                                 abs=0.005)


# ---------------------------------------------------------------------------
# The band-offset bound and its quantifier
# ---------------------------------------------------------------------------
class TestTheBandOffsetBoundMatchesItsQuantifier:
    """Referee 3, major.  '134 s bounds ANY band-to-band offset' quoted the
    tightest of five pooled pairs; the weakest allows 313 s, so a 250 s
    offset in G-I is consistent with these data and the paper said it was
    excluded.  Figure 9(b) printed all five honestly and the text did not."""

    def _pooled(self, phot):
        return phot.execute(
            "SELECT band_a, band_b, delta_s, sigma_s FROM p3_band_pair "
            "WHERE target_key='stlmi' AND lower(night) LIKE '%pooled%' "
            "AND sigma_s IS NOT NULL").fetchall()

    def test_the_weakest_bound_is_the_max_over_pooled_pairs(self, numbers,
                                                            phot):
        bounds = [abs(d) + 2.0 * s for _, _, d, s in self._pooled(phot)]
        assert _num(numbers, "NumBandOffsetWeakestBoundS") == pytest.approx(
            max(bounds), abs=1.0)
        assert _num(numbers, "NumBandOffsetBoundS") == pytest.approx(
            min(bounds), abs=1.0)
        assert max(bounds) > min(bounds), (
            "the pooled pairs no longer differ; if that is real the two "
            "macros may be merged, but not before")

    def test_every_universal_claim_uses_the_weakest_bound(self, body):
        """A sentence quantifying over the pairs --- 'any', 'every' --- may
        carry only the weakest bound.  Every use of the TIGHTEST macro must
        therefore sit in a clause that names it as the tightest, or as its
        own pair; an unqualified one is the defect this test exists for."""
        for m in re.finditer(r"\\NumBandOffsetBoundS", body):
            window = body[max(0, m.start() - 320):m.end() + 200]
            window = window.replace("\n", " ")
            named = ("tightest" in window
                     or "\\NumBandOffsetBoundPair" in window)
            assert named, (
                "\\NumBandOffsetBoundS is quoted without naming it as the "
                f"tightest pair, near: ...{window!r}")
        # And the weakest macro must actually be used, or the fix is
        # cosmetic: the paper would simply have stopped quantifying.
        assert "\\NumBandOffsetWeakestBoundS" in body

    def test_the_caption_gives_both_ends(self, captions, phot):
        cap = captions["CapFigZeroNine"]
        bounds = [abs(d) + 2.0 * s for _, _, d, s in self._pooled(phot)]
        assert f"{max(bounds):.0f}~s" in cap
        assert f"{min(bounds):.0f}~s" in cap
        assert "in EVERY pair" in cap


# ---------------------------------------------------------------------------
# Captions that must not outlive a retraction
# ---------------------------------------------------------------------------
class TestCaptionsAgreeWithTheBody:
    """Referee 3, major.  §3.1 was rewritten to retract the held-out rule for
    the headline precision; Figure 2's caption --- the caption of the panel
    those very numbers are plotted on --- went on asserting it, and asserted
    it of the noise floor too, which is fitted over the same mixed
    population.  A caption surviving a retraction is how a paper ends up
    claiming two things."""

    def test_figure_two_does_not_reassert_the_retracted_rule(self, captions):
        cap = captions["CapFigZeroTwo"]
        assert "permits only on the held-out check stars" not in cap, (
            "Figure 2's caption still carries the rule §3.1 retracts")
        assert "LOCAL FITS over the ensemble and check stars together" in cap

    def test_the_scope_clause_is_one_string_in_both_places(self, captions,
                                                           numbers, body):
        """The clause is emitted once and used twice; if it is ever forked,
        the two can drift and this fails."""
        clause = numbers["NumPrecisionScopeClause"]
        assert clause in captions["CapFigZeroTwo"], (
            "Figure 2's caption no longer uses the emitted scope clause")
        assert "\\NumPrecisionScopeClause" in body, (
            "§3.1 no longer uses the emitted scope clause")

    def test_figure_nine_never_asserts_the_epochs_lie_inside_the_envelope(
            self, captions):
        """27 of 36 epochs lie outside the drawn envelope, which is a signal
        shape and not an error band.  The caption claimed containment."""
        cap = captions["CapFigZeroNine"]
        assert "epochs sit inside it" not in cap
        assert "not an error band" in cap

    def test_figure_nines_envelope_counts_match_the_epochs(self, captions,
                                                           phot):
        """Whatever the caption says about the envelope must be recomputable
        from the same rows the panel drew."""
        from macro_phot.figures_cv import pdot_envelope_seconds
        rows = phot.execute(
            "SELECT cycle_mean, oc_s, oc_sigma_s FROM p3_oc_night "
            "WHERE target_key='stlmi' ORDER BY cycle_mean").fetchall()
        lim, per_d = phot.execute(
            "SELECT pdot_limit3, period_d FROM p3_cycle_count "
            "WHERE target_key='stlmi'").fetchone()
        e = [r[0] for r in rows]
        s = [r[2] for r in rows]
        env = pdot_envelope_seconds(e, e, s, float(lim), float(per_d))
        outside = sum(1 for r, v in zip(rows, env) if abs(r[1]) > abs(v))
        assert f"{outside} of the {len(rows)} epochs scatter further" in \
            captions["CapFigZeroNine"], (
            f"the caption does not state the {outside} epochs that lie "
            f"outside the envelope the code draws")

    def test_figure_twelve_does_not_call_the_gap_a_sensitivity_cost(
            self, captions, body):
        """For four of five series the detrended contour is BETTER than the
        raw one, so the gap has the opposite sign from the asserted cost."""
        cap = captions["CapFigOneTwo"]
        assert "the detrending costs in sensitivity" not in cap
        assert "NOT a sensitivity cost" in cap
        assert "sensitivity\ncost of detrending" not in body
        assert "shows the sensitivity" not in body


# ---------------------------------------------------------------------------
# One run, one name
# ---------------------------------------------------------------------------
class TestARunIsNamedTheSameWayEverywhere:
    """Referee 3, minor.  Table 4 named the uninformative YZ Cnc run by its
    local observing night (2024-02-20) and Figure 11's caption and axis by
    its UTC night (2024-02-21).  Same run, two dates, nothing explaining
    that two conventions were in play."""

    def test_the_verdict_names_the_run_the_figure_labels(self, phot,
                                                          captions):
        from macro_phot import final_science as fs
        row = phot.execute(
            "SELECT nights, utc_nights, filter FROM p4_run WHERE "
            "state='QUIESCENT' AND detection NOT IN ('AMPLITUDE ONLY') "
            "AND hump_amp <= amp90_field").fetchone()
        if row is None:
            pytest.skip("no scope below its instrumental contour")
        night = fs.run_night_label(row["utc_nights"] if hasattr(row, "keys")
                                   else row[1], row[0])
        verdict = phot.execute("SELECT deciding_number FROM p4_verdict "
                               "WHERE verdict_id='YZ-hump'").fetchone()[0]
        assert f"run of {night}" in verdict, (
            f"Table 4's hump row does not name the run by its UTC night "
            f"{night}")
        assert f"the {night} $" in captions["CapFigOneOne"], (
            f"Figure 11's caption does not name the same run {night}")

    def test_the_helper_prefers_utc_and_falls_back(self):
        from macro_phot import final_science as fs
        assert fs.run_night_label("2024-02-21", "2024-02-20") == "2024-02-21"
        assert fs.run_night_label(None, "2024-02-20") == "2024-02-20"
        assert fs.run_night_label("2024-05-02+2024-05-03", None) == "2024-05-02"


# ---------------------------------------------------------------------------
# Two populations, two names
# ---------------------------------------------------------------------------
class TestTheTwoCheckStarPopulationsAreDistinguished:
    """Referee 3, major.  §3.1 listed the catalogue-tie accuracy among the
    statistics carried by the four ensemble check stars.  It is measured on
    15--513 catalogue stars withheld from the TIE fit, and Table 2 prints a
    statistic of each in adjacent columns under one name."""

    def test_the_counts_really_do_differ(self, phot):
        solve = {r[0] for r in phot.execute(
            "SELECT DISTINCT n_check FROM cv_series WHERE n_target_points>0")}
        tie_lo, tie_hi = phot.execute(
            "SELECT min(n_check), max(n_check) FROM cv_cattie WHERE "
            "is_primary=1 AND verdict LIKE 'TIED%'").fetchone()
        assert solve == {4}
        assert not (tie_lo == tie_hi == 4), (
            "the two populations now coincide; if that is real the paper "
            "may stop distinguishing them, but not before")

    def test_the_body_does_not_attribute_the_tie_accuracy_to_solve_stars(
            self, body):
        # The retracted sentence listed three statistics on the four stars.
        assert ("the catalogue-tie accuracy of Section~\\ref{sec:tie}, and "
                "the error-bar\ninflation factor" not in body)
        assert "solve check stars" in body
        assert "tie check stars" in body

    def test_table_two_says_its_two_columns_count_different_stars(self):
        tables = _text("tables.tex")
        assert "TWO COLUMNS OF THIS TABLE COUNT DIFFERENT HELD-OUT STARS" \
            in tables


# ---------------------------------------------------------------------------
# The remaining minors
# ---------------------------------------------------------------------------
class TestTheMinorInconsistencies:

    def test_the_dynamic_range_ratio_is_one_measured_range(self, captions,
                                                            body, numbers):
        """§2.1 said 'nearly twenty', Figure 3 said 'a factor of 16'; the
        second is the nominal bit-depth ratio and not a measurement."""
        assert "factor of nearly twenty" not in body
        assert "\\NumDynamicRangeRatioRange" in body
        lo = _num(numbers, "NumDynamicRangeRatioMin")
        hi = _num(numbers, "NumDynamicRangeRatioMax")
        assert f"{lo:.1f}--{hi:.1f}" in captions["CapFigZeroThree"], (
            "Figure 3's caption no longer quotes the same measured range "
            "§2.1 does")

    def test_table_three_carries_both_anuma_timing_reasons(self, phot):
        """§5.3 claimed Table 3 carried both; it carried the edge count."""
        caps = {r[0] for r in phot.execute(
            "SELECT DISTINCT capability FROM p4_anuma")}
        assert any("one feature" in c for c in caps), (
            f"p4_anuma has no one-feature row: {sorted(caps)}")
        row = phot.execute(
            "SELECT measured, bar FROM p4_anuma WHERE capability LIKE "
            "'%one feature%' LIMIT 1").fetchone()
        spread, bar = phot.execute(
            "SELECT phase_spread, (SELECT value FROM p3_meta WHERE "
            "key='one_feature_bar_cycles') FROM p3_cycle_count "
            "WHERE target_key='anuma'").fetchone()
        assert row[0] == pytest.approx(float(spread), rel=1e-6)
        assert row[1] == pytest.approx(float(bar), rel=1e-6)

    def test_the_tie_worst_case_is_labelled_as_what_it_is_the_max_of(
            self, numbers):
        """'the worst case' was the largest unclipped VALUE; the largest
        SENSITIVITY to the clipping choice is a different block."""
        assert numbers["NumTieWorstBlock"] != \
            numbers["NumTieMostSensitiveBlock"]
        assert _num(numbers, "NumTieMostSensitiveRatio") >= \
            _num(numbers, "NumTieWorstRatio")

    def test_a_per_block_sentence_uses_that_blocks_own_check_count(
            self, body, numbers):
        """§3.2 said 'removing one star of 15--513', printing the whole-survey
        range where one block's own count belonged."""
        assert "of \\NumTieCheckStarsRange{}\nmoves the number" not in body
        assert "\\NumTieWorstCheckStars" in body
        assert int(_num(numbers, "NumTieWorstCheckStars")) > 0

    def test_the_two_untied_blocks_are_named_and_are_different(self, numbers,
                                                                body):
        """'the one block with no usable tie' meant EU UMa's Fast block in
        §3.2 and ST LMi's y block in §7."""
        assert numbers["NumTieUntiedBlock"] != \
            numbers["NumUntiedBlockNoTieStage"]
        assert "\\NumUntiedBlockNoTieStage" in body
        assert "\\NumTieUntiedBlock" in body

    def test_the_stlmi_y_block_is_not_called_undetected(self, phot):
        """Table 2's caption called it 'the target is undetected'; the star
        is detected on every frame and has no zero point."""
        tables = _text("tables.tex")
        rows = phot.execute(
            "SELECT n_target_rows FROM cv_series WHERE "
            "series_key='stlmi|e47|y'").fetchone()
        if rows is None:
            pytest.skip("stlmi|e47|y absent from this build")
        assert f"{rows[0]} instrumental target detections" in tables
        assert "2 solved series in which the target is undetected" not in \
            tables

    def test_figure_six_tie_bars_are_stated_in_the_text(self, body, numbers):
        """The Mode0 g-r panel carries a 200 mmag unclipped bar and §3.3
        left it to the picture."""
        assert "\\NumTieBarRangeUnclippedMmag" in body
        assert "--" in numbers["NumTieBarRangeUnclippedMmag"]


# ---------------------------------------------------------------------------
# The law the file header states
# ---------------------------------------------------------------------------
class TestNoNumberIsTyped:

    def test_every_macro_the_body_uses_is_emitted(self, body, numbers,
                                                   captions):
        used = set(re.findall(r"\\(Num[A-Za-z]+)", body))
        used |= set(re.findall(r"\\(Num[A-Za-z]+)",
                               "".join(captions.values())))
        used.discard("NumMissing")
        missing = sorted(used - set(numbers))
        assert not missing, (
            f"main.tex references macros numbers.tex does not define: "
            f"{missing}")

    def test_no_emitted_value_is_the_missing_marker(self, numbers):
        blank = sorted(k for k, v in numbers.items()
                       if "NumMissing" in v)
        assert not blank, f"unmeasured macros reached the paper: {blank}"

    def test_the_build_log_is_free_of_warnings(self):
        """The paper's thesis is that the build is clean, so the log is an
        artefact of the paper and is checked like one.  Only warnings the
        TeX engine raises about THIS document count; tectonic's own notice
        about a non-UTF-8 byte in the bundled lineno.sty is a property of
        the aastex701 class's dependencies and is not ours to fix."""
        log = MANUSCRIPT / "main.log"
        if not log.exists():
            pytest.skip("main.log not present; run tectonic first")
        text = log.read_text(errors="replace")
        bad = [ln for ln in text.splitlines()
               if re.search(r"Overfull|Underfull|LaTeX Warning|"
                            r"Token not allowed|Citation .* undefined|"
                            r"Reference .* undefined", ln)]
        assert not bad, "the build is not warning-free:\n" + "\n".join(bad)
