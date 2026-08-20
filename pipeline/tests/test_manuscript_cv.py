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


# ---------------------------------------------------------------------------
# Referee 4
# ---------------------------------------------------------------------------
class TestUntiedTargetRowsAreSplitByWhetherTheStarWasSeen:
    """Referee 4, major.  §7 --- the section whose stated purpose is that any
    statement above can be checked with one query --- called all 245 untied
    rows "target detections", and said the 235 in the two EU UMa blocks were
    points "the fitted colour relation could not place on the standard
    system".  All 235 have inst_mag, inst_mag_err and mag NULL and match an
    aperture measurement of negative flux: they are non-detections, and a
    colour relation cannot fail to place a magnitude that was never
    measured.  The split now happens at the emitter, on inst_mag."""

    @staticmethod
    def _untied_series(phot) -> list:
        """The blocks that have untied target rows at all.

        Read from ``cv_series`` rather than by scanning ``cv_lightcurve``:
        a series' untied rows are exactly ``n_target_rows -
        n_target_points``, the counts below check that the two agree, and
        naming the blocks lets every query here ride the light curve's
        ``(series_key, ...)`` index instead of scanning a 1.2 GB table.
        """
        return [r[0] for r in phot.execute(
            "SELECT series_key FROM cv_series "
            "WHERE n_target_rows > n_target_points ORDER BY series_key")]

    def test_the_two_counts_are_the_instrumental_magnitude_split(self,
                                                                 numbers,
                                                                 phot):
        keys = self._untied_series(phot)
        marks = ",".join("?" * len(keys))
        det, non, total = phot.execute(
            f"SELECT sum(inst_mag IS NOT NULL), sum(inst_mag IS NULL), "
            f"count(*) FROM cv_lightcurve WHERE series_key IN ({marks}) "
            f"AND role='target' AND cal_mag IS NULL", keys).fetchone()
        from_series = phot.execute(
            "SELECT sum(n_target_rows - n_target_points) FROM cv_series"
        ).fetchone()[0]
        assert total == from_series, (
            "cv_series and cv_lightcurve disagree about how many target "
            "rows carry no catalogue-tied magnitude")
        assert _num(numbers, "NumPointsUntiedDetections") == det
        assert _num(numbers, "NumPointsUntiedNonDetections") == non
        assert _num(numbers, "NumPointsUntied") == total

    def test_the_non_detections_really_have_no_flux(self, numbers, phot):
        """The evidence for calling them non-detections, from the apertures
        themselves rather than from the absence of a magnitude.

        ``+d.star_id`` keeps the planner off ``idx_det_star``, which for a
        target star spans every frame in the archive; the join rides the
        detections' (frame_id, det_id) key instead.
        """
        keys = self._untied_series(phot)
        marks = ",".join("?" * len(keys))
        n, neg = phot.execute(f"""
            SELECT count(*), sum(COALESCE(d.flux, 0) < 0)
            FROM cv_lightcurve l
            LEFT JOIN cv_detections d ON d.frame_id = l.frame_id
                                     AND +d.star_id = l.star_id
            WHERE l.series_key IN ({marks})
              AND l.role='target' AND l.cal_mag IS NULL
              AND l.inst_mag IS NULL""", keys).fetchone()
        assert neg == n, (
            f"only {neg} of {n} untied rows without an instrumental "
            f"magnitude have negative aperture flux; the paper's reason for "
            f"calling them non-detections no longer holds")
        assert _num(numbers, "NumPointsUntiedNonDetectionNegFlux") == neg

    def test_no_macro_calls_the_undetected_rows_detections(self, phot):
        """The note is in the release, so a reader querying p5_number finds
        it as surely as a reader of §7 finds the prose.  A count of rows
        with no instrumental magnitude may not be described as detections,
        and may not be blamed on the colour relation."""
        keys = self._untied_series(phot)
        marks = ",".join("?" * len(keys))
        non = phot.execute(
            f"SELECT count(*) FROM cv_lightcurve WHERE series_key IN "
            f"({marks}) AND role='target' AND cal_mag IS NULL "
            f"AND inst_mag IS NULL", keys).fetchone()[0]
        bad = []
        for name, body, note in phot.execute(
                "SELECT macro, value_tex, COALESCE(note,'') FROM p5_number"):
            try:
                v = int(str(body).replace("\\,", "").strip())
            except (TypeError, ValueError):
                continue
            if v != non:
                continue
            low = str(note).lower()
            if "colour relation could not place" in low:
                bad.append((name, "blames the colour relation"))
            # "detections" alone, not inside "non-detections"
            if re.search(r"(?<!non-)(?<!non )detections", low) and \
                    "non-detection" not in low:
                bad.append((name, "calls them detections"))
        assert not bad, f"p5_number misdescribes the non-detections: {bad}"

    def test_the_body_says_non_detection_and_names_the_blocks(self, body,
                                                              numbers):
        assert "\\NumPointsUntiedNonDetections" in body
        assert "\\NumPointsUntiedDetections" in body
        assert "\\NumPointsUntiedInTiedBlocks" not in body, (
            "§7 still uses the tie-verdict split, which is not what "
            "distinguishes these rows")
        assert "fitted colour\n  relation could not place" not in body
        assert "$i$ and $r$" in numbers["NumPointsUntiedNonDetectionBlocks"]


class TestTheAnUmaEdgeErrorsAreNotAllOutsideTheThreshold:
    """Referee 4, major.  §5.3 and Conclusion 6 said AN UMa's fitted-edge
    errors "all lie outside the 60 s threshold" in sentences printing 57 s
    as the range's lower end.  One fitted edge is inside, and was rejected
    on step signal-to-noise rather than on precision."""

    def _range(self, numbers, name):
        lo, hi = numbers[name].replace("\\,", "").split("--")
        return float(lo), float(hi)

    def test_the_fitted_range_really_does_straddle_the_threshold(
            self, numbers, phot):
        lo, hi = self._range(numbers, "NumAnUmaEdgeSigmaRangeS")
        thr = _num(numbers, "NumSigmaTThresholdS")
        n, outside = phot.execute(
            "SELECT count(*), sum(sigma_t_s >= ?) FROM p3_edge "
            "WHERE target_key='anuma'", (thr,)).fetchone()
        assert _num(numbers, "NumAnUmaEdgesOutsideThreshold") == outside
        assert _num(numbers, "NumAnUmaEdgesInsideThreshold") == n - outside
        # The accepted set is the one the universal claim is true of.
        alo, _ = self._range(numbers, "NumAnUmaEdgeSigmaAcceptedRangeS")
        assert alo >= thr, "the accepted edges no longer clear the threshold"
        assert lo < thr, (
            "the fitted range no longer straddles the threshold; if that is "
            "real the prose may simplify, but the test must be retired "
            "deliberately and not by drift")

    def test_no_sentence_quantifies_universally_over_the_fitted_range(
            self, body, numbers):
        """Every sentence quoting the FITTED range against the threshold is
        checked: if it universally quantifies, the range's lower end must be
        outside the threshold, and it is not."""
        lo, _ = self._range(numbers, "NumAnUmaEdgeSigmaRangeS")
        thr = _num(numbers, "NumSigmaTThresholdS")
        if lo >= thr:
            pytest.skip("the fitted range clears the threshold outright")
        universal = re.compile(r"\b(every|all|each|none|whole)\b", re.I)
        # A universal claim is legitimate only when the sentence also names
        # the set it is universal over: the ACCEPTED range, or the count of
        # fitted edges that are outside.
        qualified = ("NumAnUmaEdgeSigmaAcceptedRangeS",
                     "NumAnUmaEdgesOutsideThreshold")
        offenders = []
        for sentence in re.split(r"(?<=[.;])\s", body):
            if "NumAnUmaEdgeSigmaRangeS" not in sentence:
                continue
            if "NumSigmaTThresholdS" not in sentence:
                continue
            if not universal.search(sentence):
                continue
            if not any(q in sentence for q in qualified):
                offenders.append(" ".join(sentence.split())[:200])
        assert not offenders, (
            f"a sentence quantifies universally over the "
            f"{lo:.0f} s--... fitted range against a {thr:.0f} s threshold "
            f"its own lower end is inside, without naming the subset the "
            f"claim is true of: {offenders}")

    def test_the_products_do_not_carry_the_false_claim_either(self, phot,
                                                              numbers):
        """The same sentence was generated into p4_anuma and p4_verdict, so
        a prose-only fix would leave the release contradicting the paper."""
        thr = _num(numbers, "NumSigmaTThresholdS")
        lo = phot.execute("SELECT min(sigma_t_s) FROM p3_edge WHERE "
                          "target_key='anuma'").fetchone()[0]
        if lo >= thr:
            pytest.skip("no fitted edge inside the threshold any more")
        texts = [r[0] for r in phot.execute(
            "SELECT COALESCE(reasoning,'')||' '||COALESCE(deciding_number,'')"
            " FROM p4_anuma")]
        texts += [r[0] for r in phot.execute(
            "SELECT COALESCE(reasoning,'') FROM p4_verdict "
            "WHERE verdict_id='ANUMA-role'")]
        bad = [t[:160] for t in texts
               if re.search(r"every one of them outside", t, re.I)]
        assert not bad, f"the release still carries the false claim: {bad}"


class TestTheTransferCheckIsDescribedAsWhatItIs:
    """Referee 4, major.  §4.2 justified transporting the one-night injection
    budget by recomputing the chi-squared "from the edge fits' OWN
    Monte-Carlo errors, which are measured on each night in its own era".
    No 2024 edge has a Monte-Carlo error at all, and those epochs are
    exactly the ones whose budget was transported: on the epochs the check
    exists to validate, it falls back to the formal bar."""

    def test_the_epochs_without_a_monte_carlo_error_are_the_transported_ones(
            self, numbers, phot):
        edges = phot.execute(
            "SELECT night, filter, sigma_t_mc_s FROM p3_edge "
            "WHERE target_key='stlmi' AND accepted=1").fetchall()
        by_ep: dict = {}
        for night, filt, mc in edges:
            by_ep.setdefault((night, filt), []).append(mc)
        ocn = phot.execute(
            "SELECT night, filter, era_id FROM p3_oc_night "
            "WHERE target_key='stlmi'").fetchall()
        formal = [(n, f, e) for n, f, e in ocn
                  if not all(m is not None for m in by_ep.get((n, f), [None]))]
        transported = [(n, f, e) for n, f, e in ocn if e != 76]
        assert sorted(formal) == sorted(transported), (
            "the epochs with no Monte-Carlo error are no longer exactly the "
            "epochs whose budget was transported; §4.2's wording assumes "
            "they coincide")
        assert _num(numbers, "NumStLmiEpochsEdgeFormal") == len(formal)
        assert _num(numbers, "NumStLmiEpochsEdgeMonteCarlo") == \
            len(ocn) - len(formal)
        assert _num(numbers, "NumStLmiEpochsEdgeFormal") == \
            _num(numbers, "NumStLmiEpochsTransferredBudget")

    def test_the_body_no_longer_claims_a_per_era_measurement(self, body):
        assert "measured on each night in its own\nera" not in body
        assert "which are measured on each night in its own era" not in body
        assert "\\NumStLmiEpochsEdgeFormal" in body, (
            "§4.2 does not say how many epochs fall back to the formal bar")
        assert "\\NumStLmiOcChisqEdgeMonteCarlo" in body, (
            "§4.2 does not quote the part of the check that is what it "
            "claims to be")

    def test_limitation_three_no_longer_calls_it_the_only_reason(self, body):
        assert "which is the only\n  reason we let it stand" not in body
        assert "weaker\n  check than its name suggests" in body

    def test_the_emitters_note_says_which_errors_the_check_uses(self, phot):
        note = phot.execute(
            "SELECT note FROM p5_number WHERE macro='NumStLmiOcChisqEdge'"
        ).fetchone()[0]
        assert "rescaled formal bar" in note
        assert "measured on every night in its own era" not in note


class TestThePerEraChiSquaredCountsCloseOnEveryEpoch:
    """Referee 4, minor.  The split gave 8 High Gain and 27 Mode0 epochs
    against a set of 36, leaving the single 1MHz HS 16-bit z epoch out."""

    def test_every_era_gets_a_macro_and_they_sum_to_the_epoch_count(
            self, numbers, phot):
        eras = phot.execute(
            "SELECT era_id, count(*) FROM p3_oc_night WHERE "
            "target_key='stlmi' GROUP BY 1").fetchall()
        tags = {7: "HighGain", 47: "OneMhzHsZ", 76: "ModeZero"}
        total = 0
        for era, n in eras:
            tag = tags.get(era)
            assert tag, f"era {era} has no macro tag; the sum will not close"
            assert _num(numbers, f"NumStLmiEpochs{tag}") == n
            assert f"NumStLmiOcChisq{tag}" in numbers
            total += n
        assert total == _num(numbers, "NumStLmiOcEpochs")

    def test_the_body_quotes_all_three(self, body):
        for tag in ("HighGain", "ModeZero", "OneMhzHsZ"):
            assert f"\\NumStLmiEpochs{tag}" in body
            assert f"\\NumStLmiOcChisq{tag}" in body


class TestTheRemainingReferee4Minors:

    def test_figure_nines_envelope_names_both_ends_separately(self,
                                                              captions):
        """The envelope reaches 293 s at one end and 49 s at the other; the
        caption said "only 293 s at the ends" (plural)."""
        cap = captions["CapFigZeroNine"]
        assert "at the ends of the baseline" not in cap
        assert "at the late end of the baseline" in cap
        assert "at the\nearly end" in cap or "at the early end" in cap

    def test_the_timed_offset_and_the_fold_phase_share_a_bin(self, numbers,
                                                             phot):
        """§4.3 called 0.157 cycles "the same number" as a 0.16--0.19 range
        that does not contain it.  What is true is that they fall in one bin
        of the fold, and the paper now claims exactly that."""
        bins = int(_num(numbers, "NumFoldProfileBins"))
        offset = _num(numbers, "NumStLmiOcOffsetCycles")
        fall_lo = float(numbers["NumStLmiFallPhaseRange"].split("--")[0])
        assert int(offset * bins) == int(fall_lo * bins), (
            f"{offset} and {fall_lo} are no longer in one 1/{bins} bin")
        lo, hi = (float(x) for x in
                  numbers["NumStLmiOcOffsetBinRange"].split("--"))
        assert lo <= offset < hi
        assert "which is the same\nnumber arrived at" not in _text("main.tex")

    def test_the_order_of_magnitude_names_what_it_is_against(self, numbers,
                                                             body):
        """0.152 fails the 0.05 bar by 3.0, not by an order of magnitude;
        the order of magnitude is against ST LMi's 0.014."""
        bar_ratio = _num(numbers, "NumAnUmaPhaseSpreadOverBar")
        st_ratio = _num(numbers, "NumAnUmaPhaseSpreadOverStLmi")
        assert bar_ratio < 10 <= st_ratio, (
            "the two comparisons no longer differ in order of magnitude; "
            "the sentence may be simplified deliberately")
        assert "fail by an\norder of magnitude" not in body
        assert "\\NumAnUmaPhaseSpreadOverBar" in body
        assert "\\NumAnUmaPhaseSpreadOverStLmi" in body

    def test_the_colour_census_partitions_the_tied_blocks(self, numbers,
                                                          phot):
        """12 + 14 = 26 over a population the sentence fixed at 25.  The
        26th was the untied block, which has no colour term at all."""
        tied = _num(numbers, "NumTieTied")
        assert (_num(numbers, "NumTieColourInside")
                + _num(numbers, "NumTieColourUnsafe") == tied)
        assert (_num(numbers, "NumTieExtrapolated")
                + _num(numbers, "NumTieColourUnknown")
                == _num(numbers, "NumTieColourUnsafe"))
        db_tied = phot.execute(
            "SELECT count(*) FROM cv_cattie WHERE is_primary=1 "
            "AND verdict LIKE 'TIED%'").fetchone()[0]
        assert tied == db_tied
        inside = phot.execute(
            "SELECT count(*) FROM cv_cattie WHERE is_primary=1 AND verdict "
            "LIKE 'TIED%' AND colour_position IN "
            "('inside-span','inside-core')").fetchone()[0]
        assert _num(numbers, "NumTieColourInside") == inside

    def test_both_tie_median_ratios_are_against_the_same_goal(self, numbers):
        """"between about one and a half and six times the goal" compared
        the clipped median with the goal's upper bound and the unclipped one
        with its lower bound."""
        goal_hi = _num(numbers, "NumTieGoalHiMmag")
        clip = _num(numbers, "NumTieMedianAccuracyMmag")
        raw = _num(numbers, "NumTieMedianAccuracyUnclippedMmag")
        assert _num(numbers, "NumTieMedianRatioToGoal") == \
            pytest.approx(clip / goal_hi, abs=0.06)
        assert _num(numbers, "NumTieMedianRatioUnclippedToGoal") == \
            pytest.approx(raw / goal_hi, abs=0.06)

    def test_the_body_does_not_say_one_and_a_half_times_the_goal(self, body):
        assert "between about one and a half and six times the" not in body
        assert "\\NumTieMedianRatioToGoal" in body

    def test_figure_eleven_does_not_call_a_two_night_block_a_run(
            self, captions, phot):
        """Three of the six rows fold two nights; the caption called all six
        "quiescent dense runs" and the axis named the blocks by one night."""
        n_block = phot.execute(
            "SELECT count(*) FROM p4_run WHERE upper(state)='QUIESCENT' "
            "AND hump_amp IS NOT NULL AND amp90_self IS NOT NULL "
            "AND kind='block'").fetchone()[0]
        cap = captions["CapFigOneOne"]
        assert "For each quiescent dense run," not in cap
        if n_block:
            assert f"{n_block} two-night blocks" in cap
            assert "each row labelled with every night it folds" in cap

    def test_the_body_says_how_the_six_scopes_split(self, body, numbers,
                                                    phot):
        assert "\\NumHumpScopesBlocks" in body
        assert "\\NumHumpScopesRuns" in body
        assert (_num(numbers, "NumHumpScopesRuns")
                + _num(numbers, "NumHumpScopesBlocks")
                == _num(numbers, "NumHumpScopesTested"))
