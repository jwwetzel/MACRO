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


# ---------------------------------------------------------------------------
# The catalogue tie (CV-S6), after the 2026-08-19 review
#
# These check the PRODUCT, not the arithmetic.  Every defect below was found
# by running the stage against the real database and comparing what it
# published against what the tables underneath it say.
# ---------------------------------------------------------------------------
class TestCatalogueTie:

    def _skip(self, phot):
        if not _has(phot, "cv_cattie"):
            pytest.skip("catalogue tie not built in this checkout")

    def test_the_tie_table_carries_the_columns_the_report_reads(self, phot):
        """The report indexes ``cv_cattie`` by NAME, and the table grew two
        columns after it first shipped.  A database built by the older code
        must be migrated, not silently read short."""
        self._skip(phot)
        from macro_phot.report_cattie import TIE_COLS
        have = {r[1] for r in phot.execute("PRAGMA table_info(cv_cattie)")}
        for name in TIE_COLS:
            assert name in have, f"cv_cattie is missing {name!r}"

    def test_a_target_colour_exists_only_where_epochs_were_shared(self, phot):
        """THE defect: a variable target's colour was formed by differencing
        two campaign means, so VV Pup era 76 was published at g-r = -1.73
        when its shared-epoch colour is +0.04.  A colour may now be quoted
        only where enough near-simultaneous pairs exist to support it."""
        self._skip(phot)
        from macro_phot import cattie as ct
        # UNTIED blocks are excluded deliberately.  'unknown' means "a fit
        # exists here but the target's colour could not be measured"; an
        # UNTIED block has no fit at all, so every column including the
        # position is NULL.  Collapsing the two would throw away the more
        # informative distinction — and the case is real: when the retired
        # era 80 merged into era 78, EU UMa's combined series lost its tie.
        rows = phot.execute(
            "SELECT series_key, target_colour, colour_position, "
            "n_colour_pairs FROM cv_cattie "
            "WHERE is_primary=1 AND verdict != 'UNTIED'").fetchall()
        assert rows, "no primary tie rows"
        for skey, colour, pos, npair in rows:
            if colour is not None and math.isfinite(colour):
                assert (npair or 0) >= ct.MIN_COLOUR_PAIRS, (
                    f"{skey} quotes a colour from {npair} epoch pairs")
                assert pos != "unknown", skey
            else:
                assert pos == "unknown", (
                    f"{skey} has no colour but is placed as {pos!r}")

    def test_an_unplaceable_target_is_charged_no_extrapolation(self, phot):
        """'unknown' must not be quietly costed as if it were inside the
        range: an unmeasured colour bounds nothing."""
        self._skip(phot)
        for skey, extrap in phot.execute(
                "SELECT series_key, extrap_err FROM cv_cattie "
                "WHERE is_primary=1 AND colour_position='unknown'"):
            assert extrap is None, f"{skey} carries an extrapolation charge"

    def test_every_matched_block_has_a_measured_astrometric_offset(self, phot):
        """A block whose displacement from the catalogue was never measured
        is indistinguishable, in the product, from one measured and found
        clean.  Both must be rows."""
        self._skip(phot)
        if not _has(phot, "cv_cat_astrom"):
            pytest.skip("astrometry census not built")
        blocks = {(r[0], r[1]) for r in phot.execute(
            "SELECT DISTINCT target_key, era_id FROM cv_cat_match "
            "WHERE catalogue='refcat2'")}
        seen = {(r[0], r[1]) for r in phot.execute(
            "SELECT target_key, era_id FROM cv_cat_astrom "
            "WHERE catalogue='refcat2'")}
        assert blocks <= seen, f"unmeasured blocks: {sorted(blocks - seen)}"

    def test_a_removed_offset_actually_recovered_matches(self, phot):
        """EU UMa era 78 sat 5.2 arcsec from the catalogue, coherently, and
        was reported untieable for 'bad astrometry'.  Any offset this stage
        removes must be justified by the matches it brings back — a
        correction that changed nothing was noise."""
        self._skip(phot)
        if not _has(phot, "cv_cat_astrom"):
            pytest.skip("astrometry census not built")
        rows = phot.execute(
            "SELECT target_key, era_id, n_match_before, n_match_after, "
            "scatter_arcsec, offset_arcsec FROM cv_cat_astrom "
            "WHERE applied=1").fetchall()
        for tk, era, before, after, scat, size in rows:
            assert after > before, f"{tk} e{era}: offset changed nothing"
            assert size > scat, f"{tk} e{era}: offset is smaller than its own scatter"

    def test_no_fitted_tie_star_is_blended_in_its_own_band(self, phot):
        """Ruling 4, checked in the band that was actually tied rather than
        in the census currency.  46 equal-brightness close pairs survived
        into the VV Pup fits while the veto was read in Gaia G."""
        self._skip(phot)
        import gzip
        import json
        import numpy as np
        from macro_phot import cattie as ct
        # VV Pup: the crowded field the defect lived in, and the only one
        # where the two readings disagree often enough to be a test.
        cache = (REPO_ROOT / "products" / "phot" / "catalogue_cache"
                 / "refcat2" / "vvpup.json.gz")
        if not cache.exists():
            pytest.skip("catalogue cache not present in this checkout")
        with gzip.open(cache, "rt") as fh:
            cols = json.load(fh)["columns"]
        cra = np.asarray(cols["ra"], dtype=float)
        cdec = np.asarray(cols["dec"], dtype=float)
        cosd = float(np.cos(np.radians(np.median(cdec))))
        offenders = 0
        for band in ("gmag", "rmag", "imag"):
            cm = np.asarray(cols[band], dtype=float)
            rows = phot.execute("""
                SELECT DISTINCT m.cat_row FROM cv_cattie_star s
                JOIN cv_cattie c ON c.series_key=s.series_key
                 AND c.catalogue=s.catalogue AND c.band=s.band
                JOIN cv_cat_match m ON m.catalogue=s.catalogue
                 AND m.target_key=c.target_key AND m.era_id=c.era_id
                 AND m.star_id=s.star_id
                WHERE c.is_primary=1 AND s.in_fit=1
                  AND c.target_key='vvpup' AND c.band=?""", (band,)).fetchall()
            for (j,) in rows:
                d = np.hypot((cra - cra[j]) * cosd, cdec - cdec[j]) * 3600.0
                d[j] = np.inf
                near = d < ct.BLEND_APERTURE_ARCSEC
                if near.any() and float(np.nanmin(cm[near] - cm[j])) \
                        < ct.BLEND_DMAG:
                    offenders += 1
        assert offenders == 0, (
            f"{offenders} fitted VV Pup tie stars have a catalogue neighbour "
            f"inside {ct.BLEND_APERTURE_ARCSEC:g}\" and within "
            f"{ct.BLEND_DMAG:g} mag IN THE TIE BAND — the veto is being "
            f"taken in the wrong band again")
        # ...and the census must record that the gate was read in the tie
        # band at all, so a silent revert cannot pass the check above by
        # accident on a quieter build.
        assert phot.execute(
            "SELECT count(*) FROM cv_cattie_veto "
            "WHERE reason='blend_aperture_tieband'").fetchone()[0] > 0

    def test_the_calibrated_column_never_carries_a_colour_transformation(
            self, phot):
        """Ruling 1.  cal_mag must be exactly mag - ZP0 for every row of
        every tied series; a colour term anywhere in that column would be a
        transformation applied to blue, variable targets."""
        self._skip(phot)
        for skey, zp in phot.execute(
                "SELECT series_key, zp FROM cv_cattie WHERE is_primary=1 "
                "AND verdict LIKE 'TIED%'"):
            worst = phot.execute(
                "SELECT max(abs(cal_mag - (mag - ?))) FROM cv_lightcurve "
                "WHERE series_key=? AND cal_mag IS NOT NULL", (zp, skey)
            ).fetchone()[0]
            if worst is not None:
                assert worst < 1e-6, f"{skey}: cal_mag is not mag - ZP0"

    def test_the_per_star_table_agrees_with_the_fit_it_describes(self, phot):
        """A shipped bug: ``cv_cattie_star`` was upserted and never cleared,
        so a star the gate STARTED rejecting kept its previous row — with
        in_fit=1 — while the fit beside it had correctly dropped it.  The
        fit was right and the table published something else, which is the
        worse of the two failures because audits read the table."""
        self._skip(phot)
        for skey, cat, band, n_clean, n_fit in phot.execute(
                "SELECT series_key, catalogue, band, n_clean, n_fit "
                "FROM cv_cattie"):
            rows, fitted = phot.execute(
                "SELECT count(*), sum(in_fit) FROM cv_cattie_star "
                "WHERE series_key=? AND catalogue=? AND band=?",
                (skey, cat, band)).fetchone()
            assert rows == n_clean, (
                f"{skey}/{band}: {rows} per-star rows against n_clean="
                f"{n_clean}")
            assert (fitted or 0) == n_fit, (
                f"{skey}/{band}: {fitted} rows marked in_fit against n_fit="
                f"{n_fit}")

    def test_no_orphan_tie_rows_survive_a_rerun(self, phot):
        """Every tie row must belong to a series that is currently solved."""
        self._skip(phot)
        for tab in ("cv_cattie", "cv_cattie_star", "cv_cattie_veto"):
            n = phot.execute(
                f"""SELECT count(*) FROM {tab} WHERE series_key NOT IN
                    (SELECT series_key FROM cv_series WHERE status='solved')"""
            ).fetchone()[0]
            assert n == 0, f"{tab} carries {n} orphan rows"


# ---------------------------------------------------------------------------
# The second referee report, 2026-08-20
#
# Every test below names a statement the manuscript made that the products
# contradicted.  None of them is an arithmetic defect: the stages computed
# what their docstrings said, and the PAPER said something else, either
# because a claim was inferred from the existence of a measurement rather
# than from its significance, or because a hand-written summary string
# survived a re-run of the numbers underneath it.  These pin the products
# so the same paraphrase cannot come back.
# ---------------------------------------------------------------------------
class TestTheInterBandOffsetIsANonDetection:
    """Blocker 1.  The abstract and Conclusion 5 of the previous revision
    asserted that a bright-phase edge epoch is band dependent.  Every row of
    ``p3_band_pair`` carries ``significant=0`` and no pooled pair reaches
    even 2 sigma.  "Not uniformly zero" is true of any set of measured
    differences; it is not a detection."""

    def _skip(self, phot):
        if not _has(phot, "p3_band_pair"):
            pytest.skip("phase 3 not built in this checkout")

    def test_no_band_pair_is_significant(self, phot):
        self._skip(phot)
        rows = phot.execute("SELECT count(*), sum(significant) "
                            "FROM p3_band_pair").fetchone()
        assert rows[0] > 0
        assert (rows[1] or 0) == 0, (
            "a band pair became significant: the paper publishes this as a "
            "NON-DETECTION and its abstract, Section 5.1, Conclusion 5 and "
            "Figure 9's caption all have to be rewritten if it is one")

    def test_no_pooled_pair_reaches_two_sigma(self, phot):
        """The bound the paper publishes is 2 sigma, so the null has to hold
        at 2 sigma and not only at the 3 sigma acceptance bar."""
        self._skip(phot)
        for r in phot.execute(
                "SELECT band_a, band_b, era_id, delta_s, sigma_s FROM "
                "p3_band_pair WHERE night='(pooled)' AND sigma_s > 0"):
            assert abs(r[3]) / r[4] < 2.0, (
                f"{r[0]}-{r[1]} in era {r[2]} reaches "
                f"{abs(r[3]) / r[4]:.1f} sigma")

    def test_the_epoch_rule_does_not_call_the_offset_measured(self, phot):
        """``p3_meta.oc_epoch_rule`` is released text and justified averaging
        within a band by "a band-dependent edge epoch" it had not measured."""
        self._skip(phot)
        rule = phot.execute("SELECT value FROM p3_meta WHERE "
                            "key='oc_epoch_rule'").fetchone()
        if rule is None:
            pytest.skip("oc stage predates the epoch-rule note")
        low = str(rule[0]).lower()
        assert "does not detect" in low or "not detect" in low
        assert "conservative" in low


class TestPublishedEpochsAndTheirSpan:
    """Majors 2, 3 and 7.  The rule the paper states, the span it quotes and
    the provenance of the error bars it plots."""

    def _skip(self, phot):
        if not _has(phot, "p3_oc_night"):
            pytest.skip("phase 3 not built in this checkout")

    def test_single_cycle_epochs_exist_and_are_counted(self, phot):
        """"No single orbital cycle's edge is published as a timing epoch"
        was false: the mean of one cycle IS that cycle.  The rule is about
        the error bar, and the count of one-cycle epochs must be published
        beside it."""
        self._skip(phot)
        n_single, n_all = phot.execute(
            "SELECT sum(n_cycles = 1), count(*) FROM p3_oc_night "
            "WHERE target_key='stlmi'").fetchone()
        stored = phot.execute(
            "SELECT n_night_single_cycle FROM p3_cycle_count "
            "WHERE target_key='stlmi'").fetchone()[0]
        assert n_single > 0, (
            "if no epoch is a single cycle any more, the paper's disclosure "
            "of 18 such epochs is stale and must be regenerated")
        assert stored == n_single, (
            f"p3_cycle_count says {stored} single-cycle epochs, p3_oc_night "
            f"has {n_single} of {n_all}")

    def test_every_epoch_carries_the_budget_and_not_its_own_sigma(self, phot):
        """The rule that IS true: an epoch's error bar is the injection
        budget, never the edge fit's own sigma_t."""
        self._skip(phot)
        for r in phot.execute(
                "SELECT target_key, night, filter, oc_sigma_s, "
                "sigma_random_s, sigma_floor_s FROM p3_oc_night"):
            want = math.hypot(r[4], r[5])
            assert abs(r[3] - want) < 1e-6, f"{r[0]} {r[1]} {r[2]}"

    def test_the_epoch_span_is_not_the_catalogue_baseline(self, phot):
        """The abstract quoted 21,869 cycles as the span of the epochs; that
        is the count from the CATALOGUE epoch and is 2.5x larger."""
        self._skip(phot)
        c_first, c_last, span, span_d, t_first, t_last = phot.execute(
            "SELECT n_cycles_first, n_cycles_last, n_cycles_span, span_d, "
            "t_first_bjd, t_last_bjd FROM p3_cycle_count "
            "WHERE target_key='stlmi'").fetchone()
        assert span is not None, "span column not populated"
        assert abs(span - (c_last - c_first)) < 1e-6
        assert abs(span_d - (t_last - t_first)) < 1e-6
        assert span < 0.6 * c_last, (
            "the two counts have converged; check that the paper still "
            "distinguishes them")

    def test_the_error_budget_comes_from_one_night_and_says_so(self, phot):
        """Every published error bar is transported from a single night of a
        single target.  The release has to record that, and the epochs that
        cross an instrument seam have to be countable."""
        self._skip(phot)
        nights = phot.execute(
            "SELECT count(DISTINCT night), count(DISTINCT series_key) "
            "FROM p3_sigmat").fetchone()
        assert nights[0] == 1, (
            "the injection grid now covers more than one night; §4.2's "
            "disclosure of the transfer needs revisiting")
        src = phot.execute("SELECT value FROM p3_meta WHERE "
                           "key='timing_budget_source'").fetchone()
        assert src is not None and str(src[0]).strip(), (
            "the provenance of the timing budget is not in the release")
        # The band a row was served from must be recorded, so an epoch on
        # the whole-grid fallback is distinguishable from one on its band.
        n_null = phot.execute("SELECT count(*) FROM p3_oc_night WHERE "
                              "budget_band IS NULL").fetchone()[0]
        assert n_null == 0, "budget_band is not populated for every epoch"

    def test_the_null_survives_the_edge_fits_own_errors(self, phot):
        """The check on the transfer.  If the two chi-squareds diverge, the
        flagship null rests on the transported budget and the paper may not
        say the residuals are the demonstrated error."""
        self._skip(phot)
        chi_budget, chi_edge = phot.execute(
            "SELECT oc_night_chi2nu, oc_night_chi2nu_edge FROM "
            "p3_cycle_count WHERE target_key='stlmi'").fetchone()
        assert chi_edge is not None
        assert chi_budget < 2.0 and chi_edge < 2.0
        assert abs(chi_budget - chi_edge) < 0.5, (
            "the transported budget and the edge fits' own errors no longer "
            "agree; the transfer is doing work and must be defended")

    def test_the_period_change_null_carries_a_bound(self, phot):
        """A null is publishable when it excludes something.  The quadratic
        must be insignificant AND its 3-sigma bound must be stored."""
        self._skip(phot)
        coeff, sigma, limit = phot.execute(
            "SELECT quad_coeff_s_per_cycle2, quad_sigma_s_per_cycle2, "
            "pdot_limit3 FROM p3_cycle_count "
            "WHERE target_key='stlmi'").fetchone()
        assert limit is not None and limit > 0
        assert abs(coeff) < 3.0 * sigma, (
            "the quadratic term is now significant: this is a period-change "
            "DETECTION and the paper's headline null is wrong")


class TestVerdictStringsReproduceFromTheirSources:
    """Major 8.  ``p4_verdict.deciding_number`` is printed verbatim as
    Table 4.  Twice now a hand-written summary in it has outlived the
    numbers underneath.  These recompute the claims from the source rows."""

    def _skip(self, phot):
        if not _has(phot, "p4_verdict"):
            pytest.skip("phase 4 not built in this checkout")

    def test_the_hump_row_does_not_make_a_blanket_contour_claim(self, phot):
        """One of the six testable scopes has a fitted hump BELOW its own
        instrumental contour, so "the photometry could see a hump this size"
        is false there and the row may not say it unqualified."""
        self._skip(phot)
        num = phot.execute("SELECT deciding_number FROM p4_verdict "
                           "WHERE verdict_id='YZ-hump'").fetchone()[0]
        n_test, n_above = phot.execute(
            "SELECT count(*), sum(hump_amp > amp90_field) FROM p4_run "
            "WHERE state='QUIESCENT' AND detection NOT IN ('AMPLITUDE ONLY')"
        ).fetchone()
        assert f"{n_above} of the {n_test} testable scopes" in num, (
            "Table 4's hump row does not carry the scope-by-scope count "
            f"({n_above} of {n_test}) that Figure 11 and §5.4 use")
        assert "so the photometry could see a hump this size" not in num
        if n_above < n_test:
            assert "uninformative" in num

    def test_the_anuma_row_grades_the_estimator_the_paper_uses(self, phot):
        """§4.2 abolished per-cycle timing programme-wide; Table 3 went on
        grading AN UMa on "per-cycle bright-phase timing (O-C)"."""
        self._skip(phot)
        caps = {r[0] for r in phot.execute(
            "SELECT DISTINCT capability FROM p4_anuma")}
        assert not any(c.startswith("per-cycle") for c in caps), (
            f"p4_anuma still grades a per-cycle capability: {sorted(caps)}")

    def test_the_anuma_timing_remedy_follows_the_stored_reasons(self, phot):
        """§5.3 and §6.3 blamed a slow filter cycle; 16 of 20 rejections are
        step S/N and exactly one is a cadence gap.  A faster cycle does not
        fix an S/N-limited edge, and the remedy has to say so."""
        self._skip(phot)
        if not _has(phot, "p3_edge"):
            pytest.skip("phase 3 not built in this checkout")
        reasons = [str(r[0]) for r in phot.execute(
            "SELECT reason FROM p3_edge WHERE target_key='anuma' "
            "AND accepted=0")]
        n_snr = sum(1 for t in reasons if "step SNR" in t)
        n_gap = sum(1 for t in reasons if " gap," in t)
        assert n_snr > n_gap, (
            "the rejection mix has changed; §5.3's diagnosis and §6.3's "
            "remedy are derived from it and must be regenerated")
        row = phot.execute(
            "SELECT deciding_number, what_would_change_it FROM p4_anuma "
            "WHERE rank=2 AND filter='g'").fetchone()
        assert f"{n_snr} of {len(reasons)} rejections" in row[0]
        assert "SIGNAL-TO-NOISE" in row[1].upper()
        assert "faster filter cycle" in row[1].lower(), (
            "the remedy must still name the faster filter cycle in order "
            "to say why it is NOT the first thing to try")


class TestTheTieAccuracyIsQuotedWithItsChoice:
    """Major 5.  The paper's largest systematic is a sigma-clipped check-star
    RMS and the unclipped value is 2.6x worse.  Both have to be available,
    and the conclusion has to hold either way."""

    def _skip(self, phot):
        if not _has(phot, "cv_cattie"):
            pytest.skip("catalogue tie not built in this checkout")

    def test_both_statistics_exist_and_differ(self, phot):
        self._skip(phot)
        rows = phot.execute(
            "SELECT check_rms, check_rms_clip FROM cv_cattie "
            "WHERE is_primary=1 AND verdict LIKE 'TIED%'").fetchall()
        assert rows
        assert all(r[0] is not None and r[1] is not None for r in rows)
        assert all(r[0] >= r[1] - 1e-9 for r in rows), (
            "a clipped RMS exceeds its unclipped parent, which cannot happen")

    def test_the_calibration_goal_is_missed_on_either_statistic(self, phot):
        """The paper says the goal is missed.  If clipping ever brought the
        median inside 10--20 mmag, the sentence would depend on the choice."""
        self._skip(phot)
        rows = phot.execute(
            "SELECT check_rms, check_rms_clip FROM cv_cattie "
            "WHERE is_primary=1 AND verdict LIKE 'TIED%' "
            "ORDER BY check_rms_clip").fetchall()
        med_clip = sorted(r[1] for r in rows)[len(rows) // 2]
        med_raw = sorted(r[0] for r in rows)[len(rows) // 2]
        assert med_clip > 0.020 and med_raw > 0.020, (
            "the tie now meets the 10--20 mmag goal on one statistic; §3.2, "
            "the abstract and Conclusion 9 all say it does not")
