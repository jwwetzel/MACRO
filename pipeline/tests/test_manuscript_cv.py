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


# ---------------------------------------------------------------------------
# Referee 5
# ---------------------------------------------------------------------------
CHAR_DB = REPO_ROOT / "products" / "phot" / "cv_characterization.sqlite"
MANIFEST_DB = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

#: The three databases §7 says the release comprises, by the file name the
#: emitter records in ``p5_number.db``.
RELEASE_DBS = {
    "cv_timeseries.sqlite": PHOT_DB,
    "cv_characterization.sqlite": CHAR_DB,
    "rlmt-manifest.sqlite": MANIFEST_DB,
}


@pytest.fixture(scope="module")
def release_tables() -> dict:
    """``{database file name: {table and view names}}`` for all three."""
    out = {}
    for name, path in RELEASE_DBS.items():
        if not path.exists():
            pytest.skip(f"{name} not built in this checkout")
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
        con.execute("PRAGMA busy_timeout = 300000")
        out[name] = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type IN ('table','view')")}
        con.close()
    return out


@pytest.fixture(scope="module")
def macro_rows(phot) -> list:
    """Every p5_number row as a dict, or a skip if the stage has not run."""
    cols = {r[1] for r in phot.execute("PRAGMA table_info(p5_number)")}
    if not {"db", "kind"} <= cols:
        pytest.skip("p5_number predates the db/kind columns; re-run "
                    "run_cv_paper.py numbers")
    return [dict(zip(("macro", "value", "unit", "source", "note", "db",
                      "kind"), r))
            for r in phot.execute("SELECT macro, value_tex, unit, source, "
                                  "note, db, kind FROM p5_number")]


class TestTheReleaseIsTheDatabasesTheSectionNames:
    """Referee 5, major.  §7 opened "released as a single SQLite database"
    and closed by claiming that running the emitter "against the released
    database" reproduces every number.  Thirty-one macros --- the whole of
    Table 1, the abstract's own precision range, the injection contours,
    the check-bias ratio --- resolve in neither that database nor anything
    it references, and five of the thirteen figures are drawn from more
    than one database.  §7 now names all three and says which macros come
    from which, and the mapping is measured rather than asserted."""

    def test_every_macro_resolves_in_the_database_it_records(
            self, macro_rows, release_tables):
        """The claim under test is exactly the reader's: take a macro, open
        the database §7 sends you to, find the table it names."""
        unresolved = []
        for r in macro_rows:
            if r["kind"] == "external":
                continue
            tables = release_tables.get(r["db"])
            assert tables is not None, (
                f"{r['macro']} records database {r['db']!r}, which §7 does "
                f"not name")
            if r["source"] not in tables:
                unresolved.append(f"{r['macro']}: {r['source']} not in "
                                  f"{r['db']}")
        assert not unresolved, (
            "macros whose source table is not in the database they claim: "
            + "; ".join(unresolved))

    def test_no_external_constant_hides_in_a_products_table(
            self, macro_rows, release_tables):
        """The other half of the same invariant: a value flagged external
        must not be a query against anything released, or the flag is a
        label rather than a fact."""
        bad = [r["macro"] for r in macro_rows if r["kind"] == "external"
               and any(r["source"] in t for t in release_tables.values())]
        assert not bad, (
            f"macros flagged external whose source IS a released table: "
            f"{bad}")

    def test_section_seven_names_all_three_and_no_longer_says_one(
            self, body):
        for name in ("cv\\_timeseries.sqlite", "cv\\_characterization.sqlite",
                     "rlmt-manifest.sqlite"):
            assert name in body, f"§7 does not name {name}"
        assert "released as a single SQLite database" not in body
        assert "against the released database" not in body

    def test_the_census_macros_partition_the_macro_table(self, numbers,
                                                          macro_rows):
        """§7's own counts are emitted from p5_number, so they cannot drift
        from it: the four parts must sum to the whole."""
        total = _num(numbers, "NumMacrosTotal")
        parts = {
            "NumMacrosPhot": "cv_timeseries.sqlite",
            "NumMacrosCharacterisation": "cv_characterization.sqlite",
            "NumMacrosManifest": "rlmt-manifest.sqlite",
        }
        assert total == len(macro_rows)
        for macro, db in parts.items():
            assert _num(numbers, macro) == sum(
                1 for r in macro_rows if r["db"] == db), macro
        assert _num(numbers, "NumMacrosExternal") == sum(
            1 for r in macro_rows if r["kind"] == "external")
        assert (sum(_num(numbers, m) for m in parts)
                + _num(numbers, "NumMacrosExternal") == total)

    def test_every_table_a_figure_was_drawn_from_is_in_the_release(
            self, phot, release_tables):
        """The same claim for the figures.  §7 says running the generator
        against the release redraws every one of them, so every table a
        figure recorded reading has to be in a database §7 names."""
        rows = phot.execute(
            "SELECT fig_id, COALESCE(tables_used,'') FROM p5_figure"
        ).fetchall()
        assert rows, "no figures recorded; run run_cv_paper.py figures"
        stray = []
        for fig_id, used in rows:
            for table in [t.strip() for t in used.split(",") if t.strip()]:
                if not any(table in t for t in release_tables.values()):
                    stray.append(f"{fig_id}: {table}")
        assert not stray, (
            f"figures drawn from tables that are in none of the released "
            f"databases: {stray}")

    def test_the_figure_census_macros_match_the_figure_table(self, numbers,
                                                              phot,
                                                              release_tables):
        rows = phot.execute(
            "SELECT COALESCE(tables_used,'') FROM p5_figure").fetchall()
        dbs = []
        for (used,) in rows:
            tables = [t.strip() for t in used.split(",") if t.strip()]
            dbs.append({name for name, have in release_tables.items()
                        for t in tables if t in have})
        assert _num(numbers, "NumFiguresTotal") == len(rows)
        assert _num(numbers, "NumFiguresBeyondPhot") == sum(
            1 for d in dbs if d - {"cv_timeseries.sqlite"})
        assert _num(numbers, "NumFiguresMultiDatabase") == sum(
            1 for d in dbs if len(d) > 1)
        assert _num(numbers, "NumFiguresBeyondPhot") > 0, (
            "if no figure needs another database any more, §7's sentence "
            "about them should be retired deliberately")

    def test_the_headline_numbers_are_attributed_to_their_own_database(
            self, macro_rows):
        """The two the referee named: the abstract's per-point precision is
        not in the photometry database, and neither is any of Table 1."""
        by_macro = {r["macro"]: r for r in macro_rows}
        assert by_macro["NumPrecisionRangeMmag"]["db"] == \
            "cv_characterization.sqlite"
        table_one = [r for r in macro_rows
                     if r["source"] in ("detector_params",
                                        "s2_ceiling_modes")]
        assert table_one, "Table 1's macros have vanished from p5_number"
        assert {r["db"] for r in table_one} == {"rlmt-manifest.sqlite"}


class TestExternalConstantsAreSeparableWithOneQuery:
    """Referee 5, major.  §1 and Conclusion 10 promised that every external
    constant carries "the strategy section or the citation as its source,
    never a products table", so filtering on the source field separates
    constants from measurements.  Six macros broke it, among them the 60 s
    timing threshold quoted in the abstract, the 0.05 one-feature bar and
    the 50 mmag literature superhump floor: all three were sourced to a
    stage-meta or products table, so a reader doing what §1 describes was
    left believing they were measurements of this programme."""

    #: The constants the referee named, plus the catalogue periods.  Each
    #: must be flagged external; the test does not care which origin string
    #: it carries, only that it is not a products table.
    NAMED = ("NumSigmaTThresholdS", "NumOneFeatureBar",
             "NumStateSeparabilityBar", "NumFullOrbitMinPoints",
             "NumSuperhumpFloorMmag", "NumStLmiPeriodD", "NumStLmiPeriodMin",
             "NumVvPupPeriodD", "NumYzCncPeriodD", "NumEuUmaPeriodD",
             "NumAnUmaPeriodD")

    def test_the_named_constants_are_flagged_external(self, macro_rows):
        by_macro = {r["macro"]: r for r in macro_rows}
        wrong = [m for m in self.NAMED
                 if by_macro[m]["kind"] != "external"]
        assert not wrong, (
            f"still recorded as measurements of this programme: {wrong}")

    def test_a_note_that_names_an_outside_origin_is_never_a_measurement(
            self, macro_rows):
        """The lint that closes this at the root.  A note saying the value
        was set in advance, fixed by a stage, or published elsewhere
        describes something we did not measure, and such a value may not
        carry a products table as its source."""
        outside = re.compile(
            r"set in advance|set by CV-S\d|published VSX|"
            r"published superhump|taken from the literature", re.I)
        bad = [r["macro"] for r in macro_rows
               if r["kind"] == "measured" and outside.search(r["note"] or "")]
        assert not bad, (
            f"macros whose own note says they came from outside this "
            f"programme, sourced to a products table: {bad}")

    def test_every_external_constant_says_where_it_came_from(self,
                                                             macro_rows):
        """The inverse: a flag with no explanation is not provenance."""
        origin = re.compile(
            r"set in advance|set by CV-S\d|published|literature|catalogue",
            re.I)
        silent = [r["macro"] for r in macro_rows if r["kind"] == "external"
                  and not origin.search(r["note"] or "")]
        assert not silent, (
            f"external constants whose note does not name an origin: "
            f"{silent}")

    def test_the_promise_in_the_introduction_matches_the_table(self, body,
                                                               numbers,
                                                               macro_rows):
        n_external = sum(1 for r in macro_rows if r["kind"] == "external")
        assert _num(numbers, "NumMacrosExternal") == n_external
        assert "\\NumMacrosExternal" in body, (
            "§1 promises a reader can separate the constants but does not "
            "say how many there are to find")
        assert "never with a products table" in " ".join(body.split()), (
            "§1 no longer states the rule the flag enforces")


class TestARangeIsQuotedOverThePopulationItsSentenceIsAbout:
    """Referee 5, minor.  §4.5 printed a separability range of 0.65--0.99
    over the 11 series that "separate into two populations", in a sentence
    that had just set the bar for separating at 0.75.  The range was over
    all 15 GRADED series; over the 11 it is 0.76--0.99.  The generic lint
    below is the AN UMa range/threshold test made general, because that one
    walked a single macro by name and so could not see this."""

    def _range(self, numbers, name):
        raw = numbers[name].replace("\\,", "")
        lo, hi = raw.split("--") if "--" in raw else raw.split(" to ")
        return float(lo), float(hi)

    def test_the_bimodal_range_clears_the_bar_and_matches_the_database(
            self, numbers, phot):
        lo, hi = self._range(numbers, "NumStateSeparabilityBimodalRange")
        bar = _num(numbers, "NumStateSeparabilityBar")
        db_lo, db_hi = phot.execute(
            "SELECT min(separability), max(separability) FROM "
            "p3_state_series WHERE bimodal=1").fetchone()
        assert lo == pytest.approx(db_lo, abs=0.005)
        assert hi == pytest.approx(db_hi, abs=0.005)
        assert lo >= bar, (
            "a series counted as bimodal now sits below the bar; that would "
            "be a defect in the classifier, not in the sentence")

    def test_the_graded_range_is_still_emitted_and_is_the_wider_one(
            self, numbers, phot):
        glo, ghi = self._range(numbers, "NumStateSeparabilityRange")
        blo, _ = self._range(numbers, "NumStateSeparabilityBimodalRange")
        db_lo = phot.execute(
            "SELECT min(separability) FROM p3_state_series WHERE "
            "separability IS NOT NULL").fetchone()[0]
        assert glo == pytest.approx(db_lo, abs=0.005)
        assert glo <= blo

    def test_the_sentence_about_the_bimodal_series_quotes_their_own_range(
            self, body):
        sentences = [s for s in re.split(r"(?<=[.;])\s", body)
                     if "\\NumStateSeriesBimodal" in s]
        assert sentences, "§4.5's bimodal sentence has gone"
        for s in sentences:
            if "\\NumStateSeparability" not in s:
                continue
            assert "\\NumStateSeparabilityBimodalRange" in s, (
                "the sentence about the series that separate quotes a "
                "separability range that is not theirs: "
                + " ".join(s.split())[:200])

    def test_no_claim_straddles_a_bar_it_prints_without_saying_over_what(
            self, body, macro_rows, numbers):
        """The general form of this defect, and of the AN UMa one before
        it.  A range that STRADDLES a bar printed beside it is ambiguous by
        construction --- part of what it covers clears the bar and part does
        not --- so the claim must name the population the range is over.

        Scope is a sentence and the one before it, because that is how the
        bar reached this sentence: §4.5 set the 0.75 bar in one sentence and
        quoted the 0.65--0.99 range in the next.  Comparisons are made only
        between macros in the same unit, so a cycle count beside a threshold
        in seconds is not treated as a comparison.  A claim escapes either
        by naming the wider set in words, or by printing the counts on each
        side of the bar, which is what the paragraphs that legitimately
        straddle one already do."""
        unit = {r["macro"]: (r["unit"] or "") for r in macro_rows}
        rng = {k: v for k, v in numbers.items()
               if re.fullmatch(r"[-0-9.]+--[-0-9.]+", v.replace("\\,", ""))}
        bars = {k: v for k, v in numbers.items()
                if re.search(r"(Bar|Threshold[A-Za-z]*)$", k)
                and re.fullmatch(r"[-0-9.]+", v.replace("\\,", ""))}
        scope = re.compile(r"\bgraded\b|\bfitted\b|\battempted\b|"
                           r"\bfull grid\b|\ball \\Num", re.I)
        split_macro = re.compile(
            r"\\Num[A-Za-z]*(AtThreshold|BelowThreshold|Outside[A-Za-z]*"
            r"|Inside[A-Za-z]*)\b")
        sentences = re.split(r"(?<=[.;])\s", body)
        offenders = []
        for i, sentence in enumerate(sentences):
            window = " ".join(sentences[max(0, i - 1):i + 1])
            for rk in rng:
                if f"\\{rk}" not in sentence:
                    continue
                for bk in bars:
                    if f"\\{bk}" not in window:
                        continue
                    if unit.get(rk, "") != unit.get(bk, ""):
                        continue
                    lo, hi = (float(x) for x in
                              rng[rk].replace("\\,", "").split("--"))
                    bar = float(bars[bk].replace("\\,", ""))
                    if not lo < bar <= hi:
                        continue
                    if scope.search(window) or split_macro.search(window):
                        continue
                    offenders.append(
                        f"{rk} ({rng[rk]}) straddles {bk} ({bars[bk]}): "
                        + " ".join(window.split())[-170:])
        assert not offenders, (
            "a claim quotes a range that straddles a bar printed beside it "
            "without naming the population the range is over: "
            + " | ".join(offenders))


class TestTheMacroTableRecordsWhatSectionSevenSaysItRecords:
    """Referee 5, minor.  §7 and Conclusion 10 said the macro table records
    "the query that produced" each value.  p5_number has no query column and
    never had one --- source is a bare table name, note is prose --- and 72
    of its 293 rows carried no note either, so those rows recorded a table
    name and nothing else."""

    def test_every_macro_carries_a_note(self, macro_rows):
        silent = sorted(r["macro"] for r in macro_rows
                        if not (r["note"] or "").strip())
        assert not silent, (
            f"{len(silent)} macros record a table name and nothing else: "
            f"{silent[:12]}")

    def test_the_paper_does_not_claim_a_query_the_table_does_not_hold(
            self, body, phot):
        cols = {r[1] for r in phot.execute("PRAGMA table_info(p5_number)")}
        if "query" in cols:
            pytest.skip("p5_number now holds the SQL; the claim may return")
        flat = " ".join(body.split())
        assert "the query that produced" not in flat
        assert "recording the query and" not in flat
        assert "a note on how the value was derived" in flat, (
            "§7 no longer says what the macro table does record")


class TestFigureSevenNamesTheEraTheWayTheTablesDo:
    """Referee 5, minor.  Figure 7's caption called VV Pup's 2024 era
    "iKon" --- the camera model, a label that appears nowhere in Table 1,
    Table 2 or Figure 1's legend --- so a reader could not map it onto any
    row of the table beside it.  Table 2's column was headed "Camera" while
    its cells were readout modes, which is the same slippage the other
    way."""

    ERA_LABEL_72 = "1MHz High Sensitivity 16-bit"

    def test_the_caption_uses_the_label_the_tables_use(self, captions):
        cap = captions["CapFigZeroSeven"]
        assert "iKon" not in cap
        assert self.ERA_LABEL_72 in cap
        assert "split by READOUT MODE" in cap

    def test_that_label_really_is_in_the_instrument_table(self):
        tables = _text("tables.tex")
        assert self.ERA_LABEL_72 in tables, (
            "the caption names an era label Table 1 does not carry")

    def test_the_series_table_column_is_named_for_what_its_cells_hold(self):
        tables = _text("tables.tex")
        assert "\\colhead{Readout mode}" in tables
        assert "\\colhead{Camera}" not in tables


# ---------------------------------------------------------------------------
# Referee 6
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def manifest():
    if not MANIFEST_DB.exists():
        pytest.skip("rlmt-manifest.sqlite not built in this checkout")
    con = sqlite3.connect(f"file:{MANIFEST_DB}?mode=ro", uri=True, timeout=300)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA busy_timeout = 300000")
    yield con
    con.close()


@pytest.fixture(scope="module")
def clock_rows(manifest) -> list:
    """Exactly the rows Figure 13(b) draws, through the figure's own helper."""
    from macro_phot.figures_cv import clock_disagreement_rows
    rows = [dict(r) for r in manifest.execute(
        "SELECT * FROM s3_dateobs_audit WHERE n_frames >= 50 "
        "ORDER BY n_frames DESC")]
    if not rows:
        pytest.skip("s3_dateobs_audit is empty in this checkout")
    return clock_disagreement_rows(rows)


class TestFigureThirteenDrawsAFloorAsAFloor:
    """Referee 6, minor.  Panel (b)'s abscissa is the measured worst
    |TELUT - DATE-OBS|, but a mode whose worst frame disagrees by less than
    its own header stamp resolution is drawn AT that resolution --- the
    conservative choice, and the right one.  The defect was that the floor
    carried the SAME marker as a measurement: the 2026 blank mode stamps to
    1 s while its frames agree to 4e-05 s, so a filled marker sat a decade
    right of the panel's own 100 ms alarm line with no red count, appearing
    to contradict the §2.3 claim the panel is cited for.  The caption
    compounded it by calling the dotted guide "the timestamp resolution"
    when it is the FINEST across modes and every bar starts at its own."""

    def test_a_marker_past_the_alarm_line_is_a_floor_or_a_counted_frame(
            self, clock_rows, captions):
        """The defect, stated as an invariant.  A marker right of the 100 ms
        line either carries frames that really exceed it, or it is a floor
        --- and if it is a floor the caption must say so by name, because
        nothing else on the page distinguishes it."""
        from macro_phot.figures_cv import CLOCK_ALARM_S
        cap = captions["CapFigOneThree"]
        for d in clock_rows:
            if d["worst_s"] <= CLOCK_ALARM_S:
                continue
            if d["n_gt_100ms"]:
                continue
            assert d["floored"], (
                f"{d['label']} is drawn past the alarm line with no frame "
                f"exceeding it and is not marked as a floor")
            assert d["label"] in cap, (
                f"{d['label']} is drawn past the alarm line as a floor and "
                f"the caption does not name it")

    def test_the_caption_names_every_floor_limited_mode(self, clock_rows,
                                                         captions):
        """Generated from the audit, not typed: a mode that becomes
        floor-limited (or stops being one) moves the caption with it."""
        cap = captions["CapFigOneThree"]
        floors = [d for d in clock_rows if d["floored"]]
        assert f"{len(floors)} of the {len(clock_rows)} modes drawn are " \
               f"floor-limited" in cap, (
            f"the caption does not state the {len(floors)} floor-limited "
            f"modes of {len(clock_rows)} the code draws")
        for d in floors:
            assert d["label"] in cap, (
                f"floor-limited mode {d['label']} is unnamed in the caption")
        for d in clock_rows:
            if d["floored"]:
                continue
            assert d["worst_s"] == d["measured_s"], (
                f"{d['label']} is drawn somewhere other than its "
                f"measurement without being called a floor")

    def test_the_caption_does_not_call_the_dotted_line_the_resolution(
            self, captions):
        """It is the minimum across modes; each bar starts at its own."""
        cap = captions["CapFigOneThree"]
        assert "running from the timestamp resolution (dotted)" not in cap
        assert "Each bar starts at THAT mode's own timestamp resolution" in cap
        assert "the FINEST resolution across the" in cap
        assert "An OPEN marker is a FLOOR and not a measurement" in cap

    def test_the_floor_rule_is_the_one_the_panel_and_caption_share(self):
        """The pure helper both use, on rows chosen to sit either side of
        the rule, so a change to one cannot leave the other behind."""
        from macro_phot.figures_cv import (clock_disagreement_rows,
                                           clock_floor_clause)
        rows = [
            {"readoutm": "Mode0", "stamp_resolution_s": 0.001,
             "max_abs_s": 211.2, "n_gt_100ms": 1},
            {"readoutm": "(blank)", "stamp_resolution_s": 1.0,
             "max_abs_s": 4.0e-05, "n_gt_100ms": 0},
        ]
        got = clock_disagreement_rows(rows)
        assert [d["floored"] for d in got] == [False, True]
        assert got[0]["worst_s"] == pytest.approx(211.2)
        assert got[1]["worst_s"] == pytest.approx(1.0)
        clause = clock_floor_clause(got)
        assert "1 of the 2 modes drawn are floor-limited" in clause
        assert "(blank)" in clause and "Mode0" not in clause
        assert "not a disagreement" in clause
        assert clock_floor_clause([got[0]]).startswith("No mode is "
                                                       "floor-limited")

    def test_the_body_claim_the_panel_supports_is_still_true(self, body,
                                                              clock_rows):
        """§2.3 says the two cards agree to the timestamp resolution in every
        mode bar one counted frame.  The panel may not be read against it."""
        assert "agree\nto the timestamp resolution in every readout mode" in \
            body or "agree to the timestamp resolution in every readout " \
                    "mode" in " ".join(body.split())
        counted = [d for d in clock_rows if d["n_gt_100ms"]]
        assert len(counted) <= 1, (
            f"more than one mode now exceeds the alarm line: "
            f"{[d['label'] for d in counted]}; §2.3's single-frame sentence "
            f"no longer describes the audit")


class TestThePromiseIsWhatTheReleaseActuallyDelivers:
    """Referee 6, minor.  §1 and §7 promised, without qualification, that
    any measured statement could be checked "in one query"; §7's own closing
    paragraph then says a value computed in code over the rows of more than
    one query has no single query to record, and several macros are exactly
    that.  The promise that IS one query --- separating the external
    constants from the measurements on the kind flag --- is unaffected."""

    def test_the_unqualified_promises_are_gone(self, body):
        flat = " ".join(body.split())
        assert "data release in one query" not in flat
        assert "can be checked with one query" not in flat
        assert "check any measured statement below" not in flat

    def test_what_replaced_them_says_what_takes_one_query_and_what_does_not(
            self, body):
        flat = " ".join(body.split())
        assert "not that one query settles everything" in flat
        assert "Not every one of them is then a single query" in flat
        assert "re-running the released emitter" in flat
        assert "trace every one that is to a named table in a named " \
               "database" in flat

    def test_the_closing_retraction_still_stands(self, body):
        """The sentence that was accurate all along, and which the two
        promises above now agree with rather than contradict."""
        flat = " ".join(body.split())
        assert "has no single query to record" in flat

    def test_no_sentence_promises_one_query_for_any_statement(self, body):
        """The recurrence guard.  Both retracted sentences had the same
        shape: a universal over the paper's statements, and "one query" as
        what settles it.  A sentence may still say one query settles
        something --- the external flag really does --- and may still say
        what needs more than one; what it may not do is quantify over every
        statement and promise a single query for all of them."""
        flat = " ".join(body.split())
        # "more than one query" is the honest half of the pair and must not
        # be read as a promise of one.
        flat = flat.replace("more than one query", "several queries")
        universal = re.compile(r"\b(any|every|each)\s+(measured\s+)?"
                               r"(statement|number|value|claim)")
        qualifier = re.compile(
            r"flag|column of (?:that|the) table|not that one query|"
            r"[Nn]ot every one of them|no single query to record")
        offenders = [s for s in re.split(r"(?<=[.;])\s", flat)
                     if "one query" in s and universal.search(s)
                     and not qualifier.search(s)]
        assert not offenders, (
            "a sentence promises one query for every statement the paper "
            "makes, which §7's own closing paragraph denies: "
            + " | ".join(offenders))


class TestNoConstantIsTypedWhereAMacroHoldsIt:
    """Referee 6, minor.  Six constants were typed into the prose although a
    macro or a table cell already held them: the 40-bin fold, AN UMa's bars
    of 8 nights and 15 percentage points, the 600 s colour-pairing window,
    the 150 s half-exposure and the 90--125 minute period span.  Every one
    was correct on the day; the objection is that a re-run moves the macro
    and leaves the sentence behind, which is the single failure mode this
    paper's apparatus exists to prevent.

    The lint below is the general form.  A prose sentence may not carry a
    numeric literal that equals an emitted macro's value when that macro's
    UNIT is named in the same sentence --- which is precisely how each of
    the six read, and which no ordinary English number ("one query", "three
    databases", "five targets") can trip, because those carry no unit."""

    #: The unit strings ``p5_number`` uses, and how each is written in prose.
    #: A macro whose unit is not here is not checked, because there is no
    #: reliable way to see it in a sentence.
    UNIT_IN_PROSE = {
        "s": [r"\bs\b", r"\bseconds?\b"],
        "min": [r"\bmin\b", r"\bminutes?\b"],
        "mmag": [r"\bmmag\b"],
        "mag": [r"\bmag\b"],
        "ADU": [r"\bADU\b"],
        "cycles": [r"\bcycles?\b"],
        "nights": [r"\bnights?\b"],
        "bins": [r"\bbins?\b"],
        "per cent": [r"per\s+cent"],
        "percentage points": [r"percentage\s+points"],
        "d": [r"\bd\b"],
        "yr": [r"\byr\b"],
        "mag/h": [r"mag/h"],
    }

    #: Spelled-out numbers count as typed constants.  "a bar of eight" is
    #: the same defect as "a bar of 8" and was the harder of the two to see.
    WORD_NUMBER = {"one": "1", "two": "2", "three": "3", "four": "4",
                   "five": "5", "six": "6", "seven": "7", "eight": "8",
                   "nine": "9", "ten": "10", "eleven": "11", "twelve": "12",
                   "thirteen": "13", "fourteen": "14", "fifteen": "15",
                   "sixteen": "16", "twenty": "20"}

    LITERAL = re.compile(
        r"(?<![\\A-Za-z0-9.])(\d+(?:\.\d+)?(?:--\d+(?:\.\d+)?)?)(?![0-9])")

    @classmethod
    def _literals(cls, sentence: str) -> set:
        # Math delimiters are typography, not content: "$90$--$125$" is the
        # range "90--125" typed, and reading it as two separate numbers is
        # how this one hid.
        s = sentence.replace("$", "")
        out = {m.group(1) for m in cls.LITERAL.finditer(s)}
        out |= {cls.WORD_NUMBER[m.group(1)] for m in
                re.finditer(r"\b(" + "|".join(cls.WORD_NUMBER) + r")\b", s)}
        return out

    @staticmethod
    def _prose(body: str) -> list:
        text = "\n".join(re.sub(r"(?<!\\)%.*$", "", ln)
                         for ln in body.splitlines())
        if "\\begin{document}" in text:
            text = text[text.index("\\begin{document}"):]
        return re.split(r"(?<=[.;:])\s", text)

    def test_the_six_constants_are_macros_now(self, body, numbers):
        for macro in ("NumFoldProfileBins", "NumAnUmaColourNightsBar",
                      "NumAnUmaDutyHalfwidthBar", "NumColourPairWindowS",
                      "NumPolarHalfExposureS", "NumOrbitalPeriodRangeMin"):
            assert macro in numbers, f"{macro} is not emitted"
            assert f"\\{macro}" in body, (
                f"{macro} is emitted but the prose still types its value")
        assert "in the 40-bin folded profile" not in body
        assert "against a bar of eight" not in body
        assert "percentage points against a bar of 15" not in body
        assert "within 600~s of each other" not in body
        assert "up to 150~s for our longest" not in body
        assert "$90$--$125$~minute" not in body

    def test_those_macros_agree_with_the_tables_and_the_code(self, numbers,
                                                              phot):
        """The point of moving them: each now comes from the place that
        already held it."""
        from macro_phot.numbers_cv import COLOUR_PAIR_WINDOW_S, FOLD_BINS
        assert _num(numbers, "NumFoldProfileBins") == FOLD_BINS
        assert _num(numbers, "NumColourPairWindowS") == COLOUR_PAIR_WINDOW_S
        tables = _text("tables.tex")
        for macro, capability in (
                ("NumAnUmaColourNightsBar",
                 "three-filter colour curves (Q5)"),
                ("NumAnUmaDutyHalfwidthBar", "absolute duty cycle")):
            bar = phot.execute("SELECT bar FROM p4_anuma WHERE capability=? "
                               "LIMIT 1", (capability,)).fetchone()
            if bar is None:
                pytest.skip(f"p4_anuma has no {capability} row")
            assert _num(numbers, macro) == pytest.approx(float(bar[0]))
            assert f"& {float(bar[0]):.2f} &" in tables, (
                f"Table 3 no longer renders the bar {macro} quotes")
        span = numbers["NumOrbitalPeriodRangeMin"].split("--")
        db_lo, db_hi = phot.execute(
            "SELECT min(period_d), max(period_d) FROM p3_ephemeris "
            "WHERE period_d IS NOT NULL").fetchone()
        assert float(span[0]) == pytest.approx(1440.0 * db_lo, abs=0.5)
        assert float(span[-1]) == pytest.approx(1440.0 * db_hi, abs=0.5)

    def test_no_prose_sentence_types_a_value_a_macro_already_holds(
            self, body, macro_rows):
        """The lint that would catch the next one.  Scope is the sentence,
        because that is where a value and its unit sit together."""
        held = [(r["macro"], (r["value"] or "").replace("\\,", "").strip(),
                 (r["unit"] or "").strip())
                for r in macro_rows]
        held = [(m, v, u) for m, v, u in held
                if v and u in self.UNIT_IN_PROSE]
        offenders = []
        for sentence in self._prose(body):
            lits = self._literals(sentence)
            if not lits:
                continue
            for macro, value, unit in held:
                if value not in lits:
                    continue
                if not any(re.search(p, sentence)
                           for p in self.UNIT_IN_PROSE[unit]):
                    continue
                offenders.append(
                    f"{macro} ({value} {unit}) typed in: "
                    + " ".join(sentence.split())[:150])
        assert not offenders, (
            "the prose types a value an emitted macro already holds, so a "
            "re-run would move the macro and leave the sentence behind: "
            + " | ".join(sorted(set(offenders))))
