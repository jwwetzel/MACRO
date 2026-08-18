"""Unit tests for macro_core.batch — the S1b production-batch pure logic.

Same philosophy as test_astrom: every decision function is exercised on
hand-built cases, INCLUDING the cases that must NOT work (a residue frame
sneaking into the queue, a terminal status silently overwritten, a
starved frame reaching the solver through the QC gate).
"""

from __future__ import annotations

import sqlite3

import pytest

from macro_core import astrom, batch
from macro_core.batch import (
    POLICY_BY_STRATUM, PRIOR_MEDIAN_S, S1_VERDICTS, STRATUM_POLICY,
    allowed_transition, build_queue_rows, eta_seconds, qc_pregate,
    saturation_adu_for, scratch_wcs_name, sidecar_rel_path,
    stratum_median_s)


# ---------------------------------------------------------------------------
# helpers: candidate rows in astrom.BASE_COLS shape
# ---------------------------------------------------------------------------
def row(obs_rowid=1, **over) -> dict:
    """A CV-polar Mode0 bin2 Sloan candidate (stratum
    cv_mode0_sloan_short) unless overridden."""
    base = dict(obs_rowid=obs_rowid, path=f"rawimage/x/{obs_rowid}.fts.fz",
                target_key="stlmi", canonical_target="ST LMi",
                readoutm="Mode0", xbinning=2, filter="g", exptime=60.0,
                naxis1=4788, naxis2=3194, ra_deg=151.0, dec_deg=25.0,
                night="2025-03-01")
    base.update(over)
    return base


def metrics(**over) -> dict:
    """image_metrics-shaped dict for a healthy star field."""
    base = dict(n_sources=120, n_psf_sources=80, median_elongation=1.1,
                median_a_px=2.0, bright_median_a_px=2.5,
                saturated_fraction=0.001, bkg_rms=12.0)
    base.update(over)
    return base


# ---------------------------------------------------------------------------
# the policy table itself: it must agree with the experiment's design
# ---------------------------------------------------------------------------
class TestPolicy:
    def test_policy_covers_every_experiment_stratum_exactly_once(self):
        # The batch must have an opinion about every stratum S1 measured
        # — no more, no less (a policy for a stratum that does not exist
        # would queue nothing silently).
        experiment = {s.stratum_id for s in astrom.STRATA}
        policy = {p.stratum_id for p in STRATUM_POLICY}
        assert policy == experiment
        assert len(STRATUM_POLICY) == len(policy)   # no duplicates

    def test_verdicts_and_priors_cover_every_stratum(self):
        for p in STRATUM_POLICY:
            assert p.stratum_id in S1_VERDICTS
            assert p.stratum_id in PRIOR_MEDIAN_S

    def test_gating_follows_the_verdicts_mechanically(self):
        # The rule the module documents: GO strata run ungated, every
        # below-GO stratum runs behind the QC pre-gate.
        for p in STRATUM_POLICY:
            assert p.qc_gated == (S1_VERDICTS[p.stratum_id] != "GO"), \
                p.stratum_id

    def test_priorities_are_a_permutation_in_population_order(self):
        # Priorities must be unique and 1..N (ties would make the queue
        # order depend on rowid alone across populations).
        prios = [p.priority for p in STRATUM_POLICY]
        assert sorted(prios) == list(range(1, len(prios) + 1))
        # Population blocks appear in the mandated order:
        # CV polars -> dwarf -> SN -> facility backlog.
        seen = []
        for p in sorted(STRATUM_POLICY, key=lambda p: p.priority):
            if p.population not in seen:
                seen.append(p.population)
        assert seen == list(batch.POPULATIONS)

    def test_sn_is_gated(self):
        # The task's hard requirement: SN frames ONLY behind the gate.
        assert POLICY_BY_STRATUM["sn_gsense_broadband"].qc_gated


# ---------------------------------------------------------------------------
# queue construction
# ---------------------------------------------------------------------------
class TestBuildQueue:
    def test_orders_by_priority_then_rowid(self):
        cands = [
            row(obs_rowid=5, target_key="rhooph", exptime=5.0),  # backlog
            row(obs_rowid=3),                                    # CV short
            row(obs_rowid=9, target_key="2023ixf",
                readoutm="High Gain", xbinning=1),               # SN
            row(obs_rowid=1),                                    # CV short
        ]
        q = build_queue_rows(cands)
        assert [x["obs_rowid"] for x in q] == [1, 3, 9, 5]
        assert [x["population"] for x in q] == [
            "cv_polars", "cv_polars", "sn2023ixf", "facility"]

    def test_residue_and_unsolvable_are_excluded(self):
        cands = [
            row(obs_rowid=1, filter="hrg"),          # grism: not solvable
            row(obs_rowid=2, naxis1=8),              # window geometry
            row(obs_rowid=3, target_key="weird",
                readoutm="Odd Readout"),             # residue: no stratum
            row(obs_rowid=4),                        # the only keeper
        ]
        q = build_queue_rows(cands)
        assert [x["obs_rowid"] for x in q] == [4]

    def test_rows_carry_policy_fields_and_do_not_mutate_input(self):
        cand = row(obs_rowid=7)
        before = dict(cand)
        q = build_queue_rows([cand])
        assert q[0]["stratum_id"] == "cv_mode0_sloan_short"
        assert q[0]["qc_gated"] == 0 and q[0]["priority"] == 1
        assert cand == before                # caller's row untouched

    def test_gated_flag_set_for_below_go_strata(self):
        q = build_queue_rows([row(obs_rowid=1, target_key="2023ixf",
                                  readoutm="High Gain", xbinning=1)])
        assert q[0]["qc_gated"] == 1


# ---------------------------------------------------------------------------
# the QC pre-gate
# ---------------------------------------------------------------------------
class TestQcPregate:
    def test_healthy_field_passes(self):
        passed, diag = qc_pregate(metrics())
        assert passed and diag == "stars present (QC pass)"

    def test_starved_frame_skips(self):
        # Fewer PSF-shaped sources than the autopsy minimum: blank/cloud.
        passed, diag = qc_pregate(metrics(n_psf_sources=3))
        assert not passed and diag.startswith("starved")

    def test_hot_pixel_blank_frame_skips(self):
        # The round-4 autopsy case: enough fake 'PSF' pairs, but the
        # BRIGHTEST detections are sub-star spikes — blank in truth.
        passed, diag = qc_pregate(metrics(bright_median_a_px=0.7))
        assert not passed and diag.startswith("starved")

    def test_defocused_frame_skips(self):
        passed, diag = qc_pregate(metrics(bright_median_a_px=12.0))
        assert not passed and diag.startswith("defocused")

    def test_trailing_frame_skips(self):
        passed, diag = qc_pregate(metrics(median_elongation=2.4))
        assert not passed and diag.startswith("trailing")

    def test_saturated_frame_skips(self):
        passed, diag = qc_pregate(metrics(saturated_fraction=0.06))
        assert not passed and diag.startswith("saturated")

    def test_unreadable_metrics_skip(self):
        # n_sources None = pixels never measured: fail closed.
        passed, diag = qc_pregate({"n_sources": None})
        assert not passed and diag == "unreadable"

    def test_saturation_rail_per_family(self):
        # High Gain (incl. StackPro) uses the GSENSE clip; others 65k.
        assert saturation_adu_for("High Gain") == 3500.0
        assert saturation_adu_for("High Gain StackPro") == 3500.0
        assert saturation_adu_for("Mode0") == 65000.0
        assert saturation_adu_for(None) == 65000.0


# ---------------------------------------------------------------------------
# sidecar naming
# ---------------------------------------------------------------------------
class TestSidecarNaming:
    def test_scratch_name_strips_fz_then_extension(self):
        # Mirrors solve_one_frame: funpack 'x.fts.fz' -> 'x.fts', then
        # solve-field writes 'x.wcs'.
        assert scratch_wcs_name("a/b/x.fts.fz") == "x.wcs"
        assert scratch_wcs_name("a/b/x.fts") == "x.wcs"
        assert scratch_wcs_name("a/b/x.fits.fz") == "x.wcs"
        assert scratch_wcs_name("noext") == "noext.wcs"

    def test_sidecar_mirrors_archive_tree(self):
        assert sidecar_rel_path("rawimage/2024/x.fts.fz") \
            == "rawimage/2024/x.fts.wcs"
        assert sidecar_rel_path("rawimage/2024/x.fts") \
            == "rawimage/2024/x.fts.wcs"

    def test_two_frames_never_collide(self):
        # Same basename in different nights: directories keep them apart.
        a = sidecar_rel_path("raw/n1/x.fts.fz")
        b = sidecar_rel_path("raw/n2/x.fts.fz")
        assert a != b


# ---------------------------------------------------------------------------
# status transitions
# ---------------------------------------------------------------------------
class TestTransitions:
    def test_pending_may_finish_any_terminal_way(self):
        for t in ("solved", "bad_solve", "failed", "skipped_qc"):
            assert allowed_transition("pending", t)

    def test_terminal_states_are_immutable_except_requeue(self):
        for old in ("solved", "bad_solve", "failed", "skipped_qc"):
            for new in ("solved", "bad_solve", "failed", "skipped_qc"):
                assert not allowed_transition(old, new)
            assert allowed_transition(old, "pending")   # explicit requeue

    def test_pending_to_pending_is_not_a_transition(self):
        assert not allowed_transition("pending", "pending")


# ---------------------------------------------------------------------------
# ETA arithmetic
# ---------------------------------------------------------------------------
class TestEta:
    def test_prior_used_until_enough_live_frames(self):
        sid = "cv_mode0_sloan_short"
        assert stratum_median_s(sid, 9.9, 3) == PRIOR_MEDIAN_S[sid]
        assert stratum_median_s(sid, 9.9,
                                batch.MIN_FRAMES_FOR_LIVE_MEDIAN) == 9.9

    def test_unknown_stratum_gets_conservative_default(self):
        assert stratum_median_s("nope", None, 0) == 6.0

    def test_eta_weights_each_stratum_by_its_own_median(self):
        pend = {"a": 100, "b": 50}
        med = {"a": 2.0, "b": 4.0}
        # (100*2 + 50*4) / 10 workers = 40 s.
        assert eta_seconds(pend, med, 10) == pytest.approx(40.0)

    def test_eta_zero_workers_is_zero_not_crash(self):
        assert eta_seconds({"a": 10}, {"a": 1.0}, 0) == 0.0

    def test_full_backlog_eta_matches_the_report_ballpark(self):
        # The report projected ~4 h for the 38k stratified backlog at
        # 10 workers; the same numbers through eta_seconds must land in
        # that ballpark (this pins the priors against typos).
        backlog = {"cv_mode0_sloan_short": 2341, "cv_mode0_sloan_long": 1663,
                   "cv_ikon_sloan": 558, "cv_gsense_misc": 40,
                   "sn_gsense_broadband": 368, "dwarf_gsense_deep": 247,
                   "mode0_backlog_short": 9232, "mode0_backlog_long": 8478,
                   "fast_fullframe": 2560, "ikon_backlog": 12605}
        hours = eta_seconds(backlog, PRIOR_MEDIAN_S, 10) / 3600.0
        assert 3.0 < hours < 4.5


class TestExecuteResilient:
    """Regression (2026-08-18): a multi-hour batch died with 24,662 frames
    still pending because a sibling stage held the manifest write lock past
    the connection's 120 s timeout.  Contention must be retried; real errors
    must still surface immediately."""

    def _script(self):
        import importlib.util, sys
        from pathlib import Path
        p = (Path(__file__).resolve().parent.parent / "scripts"
             / "run_s1_batch.py")
        spec = importlib.util.spec_from_file_location("run_s1_batch", p)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["run_s1_batch"] = mod
        spec.loader.exec_module(mod)
        return mod

    def test_retries_lock_errors_then_succeeds(self, monkeypatch):
        mod = self._script()
        monkeypatch.setattr(mod.time, "sleep", lambda *_: None)  # no real wait

        class FlakyCon:
            def __init__(self):
                self.calls = 0

            def execute(self, sql, params):
                self.calls += 1
                if self.calls < 3:
                    raise sqlite3.OperationalError("database is locked")
                return None

        con = FlakyCon()
        mod.execute_resilient(con, "UPDATE t SET a=?", [1])
        assert con.calls == 3, "should have retried until the lock cleared"

    def test_real_errors_are_not_retried(self, monkeypatch):
        mod = self._script()
        monkeypatch.setattr(mod.time, "sleep", lambda *_: None)

        class BrokenCon:
            def __init__(self):
                self.calls = 0

            def execute(self, sql, params):
                self.calls += 1
                raise sqlite3.OperationalError("no such column: nonsense")

        con = BrokenCon()
        with pytest.raises(sqlite3.OperationalError):
            mod.execute_resilient(con, "UPDATE t SET nonsense=?", [1])
        assert con.calls == 1, "a genuine SQL error must fail immediately"
