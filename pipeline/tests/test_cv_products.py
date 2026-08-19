"""Invariants of the built CV products, checked against the products themselves.

The unit tests in ``test_series.py`` and ``test_characterize.py`` protect the
ARITHMETIC.  The defects an adversarial review found on 2026-08-19 were not
arithmetic defects: every function did what its docstring said, and the
product still asserted things that were not true, because a later stage
overwrote an earlier stage's column, or a column was written under a name
that no longer described it, or a cut was computed and then read by nothing.
No pure-function test can catch that class.  These do.

Each test names the published number that was wrong.  All of them skip when
the products are absent, so a fresh checkout still runs green — the products
are 1.2 GB and are not in git.
"""

from __future__ import annotations

import math
import sqlite3
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PHOT_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
CHAR_DB = REPO_ROOT / "products" / "phot" / "cv_characterization.sqlite"


def _con(path: Path) -> sqlite3.Connection:
    """Read-only, with the long busy timeout the concurrent S1 batch needs."""
    if not path.exists():
        pytest.skip(f"{path.name} not built in this checkout")
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


@pytest.fixture(scope="module")
def phot():
    con = _con(PHOT_DB)
    yield con
    con.close()


@pytest.fixture(scope="module")
def char():
    con = _con(CHAR_DB)
    yield con
    con.close()


def _has(con, table: str) -> bool:
    return bool(con.execute("SELECT count(*) FROM sqlite_master WHERE "
                            "type='table' AND name=?", (table,)).fetchone()[0])


# ---------------------------------------------------------------------------
# The photometry product
# ---------------------------------------------------------------------------
class TestSaturationIsReproducible:
    """'7.2% of all High Gain detections are saturated' has to be
    recomputable from the threshold the product records.

    It was not.  Extraction lowered the veto by the master dark's median and
    compared pixels against 3,104-3,105 ADU; a later `init` wrote the RAW
    3,200 back over the column.  10,371 stored flags disagreed with the only
    threshold the product still held.
    """

    def _mismatches(self, phot, column: str) -> int:
        return phot.execute(f"""
            SELECT count(*) FROM cv_detections d JOIN cv_frames f
              USING(frame_id)
            WHERE f.provenance <> 'server_reduced' AND d.peak IS NOT NULL
              AND {column} IS NOT NULL
              AND d.saturated <> (CASE WHEN d.peak + coalesce(f.bkg_adu, 0.0)
                                       >= {column} THEN 1 ELSE 0 END)"""
        ).fetchone()[0]

    def test_the_applied_veto_reproduces_every_stored_flag(self, phot):
        assert self._mismatches(phot, "f.veto_applied_adu") == 0

    def test_the_raw_veto_does_not(self, phot):
        """The other half of the proof: if BOTH thresholds reproduced the
        flags, the applied column would be telling us nothing."""
        assert self._mismatches(phot, "f.veto_adu") > 0

    def test_locally_calibrated_frames_carry_the_dark_shift(self, phot):
        rows = phot.execute(
            "SELECT veto_adu, veto_applied_adu, dark_median_adu FROM cv_frames "
            "WHERE provenance='local_master' AND master_dark IS NOT NULL "
            "AND veto_adu IS NOT NULL").fetchall()
        assert rows
        for raw, applied, dark in rows:
            assert applied == pytest.approx(raw - dark, abs=1e-6)

    def test_server_reduced_frames_are_not_shifted_twice(self, phot):
        """Their veto was already mapped through the reduction; subtracting a
        local dark as well would move the ceiling a second time."""
        n = phot.execute(
            "SELECT count(*) FROM cv_frames WHERE provenance='server_reduced' "
            "AND veto_applied_adu IS NOT NULL "
            "AND abs(veto_applied_adu - veto_adu) > 1e-6").fetchone()[0]
        assert n == 0


class TestGeometryIsResolvedNotAsserted:
    """EU UMa's 207 era-80 frames read 8 x 3,211 in the product for the whole
    production run — the BINTABLE row length of a tile-compressed header —
    after the campaign had already established that the files hold
    4,800 x 3,211 raw / 4,787 x 3,193 reduced images and had un-excluded them
    on exactly that basis."""

    def test_no_frame_claims_an_impossible_image_size(self, phot):
        from macro_phot import series as sr
        bad = [r for r in phot.execute(
            "SELECT frame_id, naxis1, naxis2 FROM cv_frames")
            if sr.geometry_is_implausible(r[1], r[2])]
        assert not bad, f"{len(bad)} frames carry an impossible geometry"

    def test_era_80_is_recorded_at_its_real_size(self, phot):
        rows = phot.execute(
            "SELECT DISTINCT naxis1, naxis2, geom_basis FROM cv_frames "
            "WHERE era_id = 80").fetchall()
        if not rows:
            pytest.skip("era 80 not present in this product")
        for nx, ny, basis in rows:
            assert nx > 1000 and ny > 1000
            assert basis == "resolved"


class TestSeriesSummariesMeanWhatTheySay:
    """cv_series.n_target_points counted light-curve ROWS, including rows
    whose magnitude is NULL, and overstated three series by up to 8.7x.  And
    stlmi|e47|y was recorded as solved, check-'validated', with 10 target
    points — all ten NULL, with no target row in cv_stars at all."""

    def test_target_points_counts_measurements(self, phot):
        for sk, n in phot.execute(
                "SELECT series_key, n_target_points FROM cv_series"):
            real = phot.execute(
                "SELECT count(*) FROM cv_lightcurve WHERE series_key=? "
                "AND role='target' AND mag IS NOT NULL", (sk,)).fetchone()[0]
            assert n == real, f"{sk}: {n} recorded vs {real} measured"

    def test_target_rows_is_kept_separately(self, phot):
        """The row count is still available — it just is not the headline."""
        over = phot.execute(
            "SELECT count(*) FROM cv_series "
            "WHERE n_target_rows < n_target_points").fetchone()[0]
        assert over == 0

    def test_a_series_without_a_target_says_so(self, phot):
        for sk, verdict, n in phot.execute(
                "SELECT series_key, target_verdict, n_target_points "
                "FROM cv_series WHERE status='solved'"):
            if n == 0:
                assert verdict == "undetected", \
                    f"{sk} has no target magnitudes but reads {verdict!r}"

    def test_the_near_ceiling_count_exists_where_the_ceiling_is_close(self, phot):
        """n_target_saturated is 0 in every series, which reads as 'saturation
        never touched the targets'.  YZ Cnc's era-7 target reaches 98.8% of
        the applied veto on the block that carries the superhump amplitude."""
        row = phot.execute(
            "SELECT sum(n_target_saturated), sum(n_target_near_veto_090) "
            "FROM cv_series WHERE status='solved'").fetchone()
        assert row[0] == 0                    # the flag really is zero ...
        assert row[1] and row[1] > 0          # ... and it is not the whole story

    def test_calibration_age_is_recorded(self, phot):
        """Era 7 is flat-fielded with masters 100-150 days old, and that
        appeared nowhere in the product while the report named flat-field
        residual as the leading suspect for the noise floor."""
        n, oldest = phot.execute(
            "SELECT count(*), max(flat_age_days) FROM cv_frames "
            "WHERE flat_age_days IS NOT NULL").fetchone()
        assert n > 0
        assert oldest > 90


# ---------------------------------------------------------------------------
# The characterization product
# ---------------------------------------------------------------------------
class TestTheQualityCutIsActuallyApplied:
    """The cut was computed, defended over a page, and read by nothing: every
    later stage selected `FROM cv_frames WHERE status='matched'` with no
    quality filter, while section 2 stated that the noise floor was therefore
    'a property of the instrument and the pipeline'."""

    def test_the_noise_stage_used_only_usable_frames(self, char):
        for sk, used, total in char.execute(
                "SELECT series_key, n_frames_usable, n_frames_all "
                "FROM ch_noise_series"):
            usable = char.execute(
                "SELECT count(*) FROM ch_frames WHERE series_key=? AND "
                "usable=1", (sk,)).fetchone()[0]
            assert used == usable, f"{sk}: noise used {used} of {usable} usable"
            assert used <= total

    def test_the_cut_is_not_a_no_op(self, char):
        """If every frame passed, applying the cut would prove nothing."""
        used = char.execute("SELECT sum(n_frames_usable), sum(n_frames_all) "
                            "FROM ch_noise_series").fetchone()
        assert used[0] < used[1]

    def test_the_cadence_stage_used_only_usable_frames(self, char):
        for sk, n in char.execute(
                "SELECT series_key, n_points FROM ch_cadence"):
            usable = char.execute(
                "SELECT count(*) FROM ch_frames WHERE series_key=? AND "
                "usable=1", (sk,)).fetchone()[0]
            assert n == usable

    def test_airmass_carries_a_measured_threshold(self, char):
        """'Airmass gets no threshold at all - over the range actually
        observed it never degrades the check stars' was published while three
        consecutive bins above X = 2.45 sat above the degradation factor,
        excluded one by one for holding under 15 frames each."""
        thr = char.execute("SELECT threshold FROM ch_cuts "
                           "WHERE axis='airmass'").fetchone()[0]
        assert thr is not None and math.isfinite(thr)


class TestContoursAreLabelledByQuestion:
    """A single-night contour scored on 1% period tolerance was published as a
    detection limit and fed into Q3's margin and S4's headline ratio."""

    def test_both_score_modes_are_stored(self, char):
        from macro_phot import characterize as ch
        got = {r[0] for r in char.execute(
            "SELECT DISTINCT score FROM ch_contour")}
        assert got == set(ch.SCORE_MODES)

    def test_the_two_scores_disagree_on_a_single_night(self, char):
        """If they agreed, the distinction would be bookkeeping.  They do not:
        on one night the period-determination limit is several times the
        detection limit."""
        rows = char.execute("""
            SELECT k.series_key, k.amp90, p.amp90 FROM ch_contour k
            JOIN ch_contour p ON p.scope = k.scope AND p.regime = k.regime
                 AND abs(p.period_d - k.period_d) < 1e-9 AND p.score='period'
            JOIN ch_cadence c ON c.series_key = k.series_key
            WHERE k.score='known' AND k.regime='night'
              AND abs(k.period_d - c.period_d) < 1e-6
              AND k.amp90 IS NOT NULL AND p.amp90 IS NOT NULL""").fetchall()
        assert rows
        assert all(p > 2.0 * k for _sk, k, p in rows)

    def test_every_contour_carries_an_uncertainty(self, char):
        """A90 was printed to 0.1 mmag from 50 trials per cell."""
        n_bad = char.execute(
            "SELECT count(*) FROM ch_contour WHERE amp90 IS NOT NULL "
            "AND (amp90_lo IS NULL OR amp90_hi IS NULL)").fetchone()[0]
        assert n_bad == 0

    def test_the_threshold_spread_over_stars_is_recorded(self, char):
        """The pooled threshold is a quantile of a max statistic over at most
        four held-out stars, so its bootstrap error bar is not its real
        uncertainty."""
        if not _has(char, "ch_threshold"):
            pytest.skip("threshold spread not built")
        n = char.execute("SELECT count(*) FROM ch_threshold "
                         "WHERE spread_frac IS NOT NULL").fetchone()[0]
        assert n > 0


class TestRedNoiseIsLabelledByTau:
    """The factor was stored as `red_factor_porb` whatever tau it was measured
    at, and most ladders never reach P_orb."""

    def test_the_tau_used_is_never_above_the_target(self, char):
        n = char.execute(
            "SELECT count(*) FROM ch_allan_fit WHERE tau_used_s IS NOT NULL "
            "AND tau_used_s > tau_target_s * 1.000001").fetchone()[0]
        assert n == 0

    def test_the_label_would_have_been_wrong(self, char):
        """Most ladders stop well short of the orbital period, so a column
        named after P_orb was a lower bound wearing a larger name."""
        short = char.execute(
            "SELECT count(*) FROM ch_allan_fit "
            "WHERE tau_frac_of_porb < 0.99").fetchone()[0]
        total = char.execute("SELECT count(*) FROM ch_allan_fit").fetchone()[0]
        assert total > 0
        assert short > total / 2

    def test_every_ladder_carries_its_own_white_null(self, char):
        """-0.50 is the asymptotic white value, not this estimator's."""
        n_bad = char.execute(
            "SELECT count(*) FROM ch_allan_fit "
            "WHERE slope IS NOT NULL AND slope_null_p95 IS NULL").fetchone()[0]
        assert n_bad == 0
        med_null = char.execute(
            "SELECT median(slope_null_p50) FROM ch_allan_fit").fetchone()[0]
        assert med_null < -0.50


class TestVerdictsRestOnTheirOwnEstimands:
    def test_every_goal_has_a_verdict_and_a_number(self, char):
        rows = char.execute(
            "SELECT goal_id, verdict, deciding_number FROM ch_verdict").fetchall()
        assert len(rows) >= 10
        for gid, verdict, num in rows:
            assert verdict in ("SUPPORTED", "SUPPORTED-WITH-CAVEATS",
                               "NOT SUPPORTED", "NOT MEASURED"), gid
            assert num and num != "not measured", gid

    def test_q4_is_decided_by_epochs_not_by_per_point_precision(self, char):
        """Q4 asks for a duty cycle — a fraction of TIME in a state — and was
        graded on 9-77 mmag per point against a 1-3 mag state separation,
        which only establishes that one night can be classified."""
        num = char.execute("SELECT deciding_number FROM ch_verdict "
                           "WHERE goal_id='Q4'").fetchone()[0]
        assert "nights" in num
        assert "pp on a 50% duty cycle" in num

    def test_q4_and_s1_do_not_contradict_each_other(self, char):
        """Q4's 'embedded in the survey record' is the operation S1 grades NOT
        SUPPORTED.  They cannot both stand as originally written."""
        v = dict(char.execute("SELECT goal_id, verdict FROM ch_verdict"))
        if v.get("S1") == "NOT SUPPORTED":
            assert v.get("Q4") != "SUPPORTED"

    def test_q1_is_decided_by_a_colour_error(self, char):
        """A colour point costs at least sqrt(2) of a single band, and here
        the non-simultaneity term usually dominates."""
        num = char.execute("SELECT deciding_number FROM ch_verdict "
                           "WHERE goal_id='Q1'").fetchone()[0]
        assert "COLOUR-point error" in num

    def test_s3_is_decided_by_a_measured_misidentification_rate(self, char):
        """Window power is the input to the alias question, not the answer."""
        num = char.execute("SELECT deciding_number FROM ch_verdict "
                           "WHERE goal_id='S3'").fetchone()[0]
        assert "misidentification rate" in num
        assert char.execute(
            "SELECT count(*) FROM ch_alias_confusion").fetchone()[0] > 0
