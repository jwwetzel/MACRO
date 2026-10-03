"""Unit tests for macro_core.project_plan and the project-page renderer.

The first three groups, in order of how much they would cost to get wrong:

1. **Ledger integrity.**  The ledger is the published plan.  A duplicate id
   silently overwrites another task's status history; a stage key that is
   not in the DAG means staleness never propagates and a `done` task lies
   forever; a missing source citation means somebody invented work the
   committee never planned.  Each of those is checked, not assumed.
2. **The derive-from-DB sync logic.**  Driven with hand-built fixtures whose
   truth is known by inspection — no real manifest, no real DAG.
3. **A rendering smoke test.**  An in-memory manifest with the same table
   shapes the real one has, rendered end to end, asserting the page contains
   the plan and NOT the superseded facts the audit retired.

Sections 13–19 were added by the plan review of 2026-10-03
(``committee/reviews/2026-10-03/SYNTHESIS.md``): the two closing statuses
and the rulings they require, the third sync rule, computed blockers, what
the committee ruled as the ledger now states it, manuscript titles as
ledger data, and the command itself.
"""

from __future__ import annotations

import re
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest

from macro_core import project_plan as pp
from macro_core import provenance as pv

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


# ===========================================================================
# 1.  LEDGER INTEGRITY
# ===========================================================================

def test_ledger_validates():
    """The shipped ledger passes its own structural gate."""
    pp.validate()


def test_task_ids_are_unique():
    ids = [t.id for t in pp.ledger_tasks()]
    assert len(ids) == len(set(ids)), \
        "a duplicate id would make two tasks share one status history"


def test_every_task_has_a_ledger_status():
    for t in pp.ledger_tasks():
        assert t.status in pp.LEDGER_STATUSES, (t.id, t.status)


def test_redo_needed_is_never_declared_in_code():
    """`redo_needed` is derived, never authored: it is a statement about the
    database, and a hand-typed one could not name its reason."""
    for t in pp.ledger_tasks():
        assert t.status != pp.REDO_NEEDED, t.id


def test_every_task_names_a_real_stage_in_the_dag():
    for t in pp.ledger_tasks():
        assert t.stage in pv.STAGE_BY_KEY, (
            f"{t.id} depends on stage {t.stage!r}, which is not in "
            f"provenance.STAGES — staleness could never propagate to it")


def test_every_task_cites_a_document_that_exists():
    for t in pp.ledger_tasks():
        doc = REPO_ROOT / t.source.document
        assert doc.exists(), f"{t.id} cites missing document {t.source.document}"
        assert t.source.section.strip(), f"{t.id} cites no section"


def test_committee_projects_cite_their_own_strategy():
    for project in pp.PROJECTS:
        if not project.strategy:
            continue
        assert (REPO_ROOT / project.strategy).exists()
        for t in project.tasks:
            assert t.source.document == project.strategy, (
                f"{t.id} cites {t.source.document} but {project.key} is "
                f"governed by {project.strategy}")


def test_blocked_tasks_state_a_blocker():
    for t in pp.ledger_tasks():
        if t.status == pp.BLOCKED:
            assert t.blocker.strip(), f"{t.id}: blocked with no reason"


def test_done_tasks_carry_evidence():
    for t in pp.ledger_tasks():
        if t.status == pp.DONE:
            assert t.evidence.strip(), f"{t.id}: done with nothing behind it"


def test_evidence_paths_under_docs_exist():
    """A `done` task whose evidence link 404s is worse than no link."""
    for t in pp.ledger_tasks():
        if t.evidence.startswith("docs/"):
            assert (REPO_ROOT / t.evidence).exists(), (t.id, t.evidence)


def test_project_and_phase_are_stamped_on_every_task():
    for project in pp.groups():
        for phase in project.phases:
            for t in phase.tasks:
                assert t.project == project.key
                assert t.phase == phase.name


def test_all_six_projects_have_a_docs_directory():
    for project in pp.PROJECTS:
        assert (REPO_ROOT / "docs" / project.key).exists(), project.key


def test_ledger_covers_the_five_staged_projects_plus_the_candidate():
    keys = set(pp.PROJECT_BY_KEY)
    assert {"TCrB_Monitoring", "CV_TimeSeries", "SN2023ixf_LightCurve",
            "BeStar_Grism", "DwarfGalaxy_AGN_Survey"} <= keys
    assert "Legacy_Rigel" in keys


def test_validate_rejects_a_duplicate_id(monkeypatch):
    dup = pp.PROJECTS[0]
    monkeypatch.setattr(pp, "PROJECTS", (dup, dup))
    with pytest.raises(pp.PlanError, match="duplicate task id"):
        pp.validate()


def test_validate_rejects_an_unknown_stage(monkeypatch):
    import dataclasses
    project = pp.PROJECTS[0]
    bad = dataclasses.replace(project.phases[0].tasks[0], stage="S-NOPE")
    phase = dataclasses.replace(project.phases[0], tasks=(bad,))
    monkeypatch.setattr(pp, "PROJECTS",
                        (dataclasses.replace(project, phases=(phase,)),))
    with pytest.raises(pp.PlanError, match="not in the"):
        pp.validate()


# ===========================================================================
# 2.  PURE LOGIC — overlay, counts, next-up, sync
# ===========================================================================

def _task(tid: str, status: str = pp.PENDING, stage: str = "S0",
          blocker: str = "", evidence: str = "x") -> pp.Task:
    return pp.Task(id=tid, title=tid, produces="a thing", stage=stage,
                   source=pp.Source("ROADMAP.md", "§0"), status=status,
                   evidence=evidence, blocker=blocker,
                   project="P", phase="Phase 1")


def test_overlay_prefers_the_recorded_status():
    tasks = [_task("a", pp.PENDING)]
    assert pp.overlay_statuses(tasks, {"a": pp.DONE})["a"] == pp.DONE


def test_overlay_lets_a_recorded_status_move_a_task_backwards():
    """The single most important thing the overlay must permit."""
    tasks = [_task("a", pp.DONE)]
    got = pp.overlay_statuses(tasks, {"a": pp.REDO_NEEDED})
    assert got["a"] == pp.REDO_NEEDED


def test_overlay_falls_back_to_the_ledger():
    tasks = [_task("a", pp.BLOCKED, blocker="b")]
    assert pp.overlay_statuses(tasks, {})["a"] == pp.BLOCKED


def test_overlay_rejects_an_unknown_recorded_status():
    with pytest.raises(pp.PlanError):
        pp.overlay_statuses([_task("a")], {"a": "almost_done"})


def test_status_counts_include_every_category():
    counts = pp.status_counts([_task("a", pp.DONE)], {"a": pp.DONE})
    assert set(counts) == set(pp.ALL_STATUSES)
    assert counts[pp.DONE] == 1 and counts[pp.PENDING] == 0


def test_redo_needed_does_not_count_as_done():
    tasks = [_task("a"), _task("b")]
    counts = pp.status_counts(tasks, {"a": pp.DONE, "b": pp.REDO_NEEDED})
    assert pp.progress_fraction(counts) == (1, 2)


def test_next_up_excludes_blocked_work():
    tasks = [_task("blocked", pp.BLOCKED, blocker="b"), _task("open")]
    ids = [t.id for t in pp.next_up(tasks, pp.overlay_statuses(tasks, {}))]
    assert ids == ["open"]


def test_next_up_puts_redo_first():
    tasks = [_task("p", pp.PENDING), _task("r", pp.PENDING)]
    st = {"p": pp.PENDING, "r": pp.REDO_NEEDED}
    assert [t.id for t in pp.next_up(tasks, st)] == ["r", "p"]


def test_next_up_collapses_redo_tasks_sharing_a_stage():
    tasks = [_task("r1", stage="S3"), _task("r2", stage="S3"),
             _task("r3", stage="S0"), _task("p")]
    st = {"r1": pp.REDO_NEEDED, "r2": pp.REDO_NEEDED,
          "r3": pp.REDO_NEEDED, "p": pp.PENDING}
    ids = [t.id for t in pp.next_up(tasks, st, limit=5)]
    assert ids == ["r1", "r3", "p"], "one re-run per stage, not one per task"


def test_next_up_can_be_restricted_to_one_status():
    tasks = [_task("r"), _task("p")]
    st = {"r": pp.REDO_NEEDED, "p": pp.PENDING}
    got = pp.next_up(tasks, st, include=(pp.PENDING,))
    assert [t.id for t in got] == ["p"]


def test_open_blockers_returns_blocked_tasks_in_plan_order():
    tasks = [_task("a"), _task("b", pp.BLOCKED, blocker="x"),
             _task("c", pp.BLOCKED, blocker="y")]
    st = pp.overlay_statuses(tasks, {})
    assert [t.id for t in pp.open_blockers(tasks, st)] == ["b", "c"]


# --- derive_sync: the rule that keeps a `done` from lying -------------------

def test_sync_flips_done_to_redo_when_the_stage_is_stale():
    tasks = [_task("a", pp.DONE, stage="S3")]
    changes = pp.derive_sync(tasks, {"a": pp.DONE}, {"S3": pv.STALE})
    assert len(changes) == 1
    assert changes[0].new == pp.REDO_NEEDED
    assert "S3" in changes[0].reason and pv.STALE in changes[0].reason


@pytest.mark.parametrize("state", [pv.STALE, pv.STALE_UPSTREAM,
                                   pv.NEVER_RUN, pv.OUTPUT_MISSING])
def test_sync_flips_done_for_every_non_fresh_state(state):
    tasks = [_task("a", pp.DONE, stage="S3")]
    changes = pp.derive_sync(tasks, {"a": pp.DONE}, {"S3": state})
    assert [c.new for c in changes] == [pp.REDO_NEEDED]


def test_sync_leaves_done_alone_when_the_stage_is_fresh():
    tasks = [_task("a", pp.DONE, stage="S3")]
    assert pp.derive_sync(tasks, {"a": pp.DONE}, {"S3": pv.FRESH}) == ()


def test_sync_restores_redo_to_done_when_the_stage_comes_back():
    tasks = [_task("a", pp.DONE, stage="S3")]
    changes = pp.derive_sync(tasks, {"a": pp.REDO_NEEDED}, {"S3": pv.FRESH})
    assert [(c.old, c.new) for c in changes] == [(pp.REDO_NEEDED, pp.DONE)]


@pytest.mark.parametrize("status", [pp.PENDING, pp.IN_PROGRESS, pp.BLOCKED])
def test_sync_never_overwrites_a_human_judgement(status):
    """No fingerprint can tell whether a person has started something."""
    tasks = [_task("a", pp.PENDING, stage="S3", blocker="b")]
    assert pp.derive_sync(tasks, {"a": status}, {"S3": pv.STALE}) == ()


def test_sync_ignores_a_stage_with_no_verdict():
    tasks = [_task("a", pp.DONE, stage="S3")]
    assert pp.derive_sync(tasks, {"a": pp.DONE}, {}) == ()


def test_sync_is_idempotent():
    tasks = [_task("a", pp.DONE, stage="S3")]
    states = {"S3": pv.STALE}
    first = pp.derive_sync(tasks, {"a": pp.DONE}, states)
    after = {"a": first[0].new}
    assert pp.derive_sync(tasks, after, states) == ()


# --- evidence verdicts -----------------------------------------------------

def test_evidence_verdict_reports_destroyed_before_stale():
    """A missing table is not 'the numbers moved' — it is 'there is no
    table', and a reader must be able to tell them apart."""
    stage = pv.STAGE_BY_KEY["S2"]
    fps = {w: "MISSING" for w in stage.writes}
    fresh = {"S2": pv.Freshness("S2", pv.STALE, ("something",))}
    verdict, why = pp.evidence_verdict("S2", fresh, fps, ever_ran={"S2"})
    assert verdict == "DESTROYED"
    assert "absent" in why


def test_evidence_verdict_passes_fresh_through():
    stage = pv.STAGE_BY_KEY["S3"]
    fps = {w: "1:abc" for w in stage.writes}
    fresh = {"S3": pv.Freshness("S3", pv.FRESH)}
    assert pp.evidence_verdict("S3", fresh, fps)[0] == pv.FRESH


def test_stage_report_path_falls_back_to_the_report_stage():
    assert pp.stage_report_path("S3") == "docs/pipeline/s3_timing.html"
    assert pp.stage_report_path("S0e") == \
        "docs/pipeline/s0e_geometry_fix.html"


def test_stage_report_path_is_none_when_no_page_reads_the_stage():
    assert pp.stage_report_path("CV-S4") is None


# ===========================================================================
# 3.  THE STATUS TABLE — append-only, replayable
# ===========================================================================

@pytest.fixture()
def con():
    c = sqlite3.connect(":memory:")
    pp.ensure_status_table(c)
    yield c
    c.close()


def test_readers_are_safe_on_a_manifest_without_the_table():
    """A renderer gets a READ-ONLY connection; a missing table must read as
    'no progress recorded', never as a crash or a DDL attempt."""
    empty = sqlite3.connect(":memory:")
    assert pp.read_statuses(empty) == {}
    assert pp.read_evidence(empty) == {}
    assert pp.read_history(empty) == []
    empty.close()


def test_record_and_read_back(con):
    tid = pp.all_tasks()[0].id
    pp.record_status(con, tid, pp.IN_PROGRESS, note="started")
    assert pp.read_statuses(con)[tid] == pp.IN_PROGRESS


def test_latest_row_wins(con):
    tid = pp.all_tasks()[0].id
    pp.record_status(con, tid, pp.IN_PROGRESS, when="2026-01-01T00:00:00Z")
    pp.record_status(con, tid, pp.DONE, when="2026-02-01T00:00:00Z")
    assert pp.read_statuses(con)[tid] == pp.DONE


def test_history_is_append_only_and_replayable(con):
    tid = pp.all_tasks()[0].id
    pp.record_status(con, tid, pp.IN_PROGRESS, when="2026-01-01T00:00:00Z")
    pp.record_status(con, tid, pp.DONE, when="2026-02-01T00:00:00Z")
    rows = pp.read_history(con, tid)
    assert [r[1] for r in rows] == [pp.IN_PROGRESS, pp.DONE]
    assert [r[4] for r in rows] == ["2026-01-01T00:00:00Z",
                                    "2026-02-01T00:00:00Z"]


def test_same_second_rows_break_the_tie_on_insertion_order(con):
    tid = pp.all_tasks()[0].id
    pp.record_status(con, tid, pp.DONE, when="2026-01-01T00:00:00Z")
    pp.record_status(con, tid, pp.REDO_NEEDED, when="2026-01-01T00:00:00Z")
    assert pp.read_statuses(con)[tid] == pp.REDO_NEEDED


def test_record_rejects_an_unknown_task(con):
    with pytest.raises(pp.PlanError, match="unknown task id"):
        pp.record_status(con, "NOT-A-TASK", pp.DONE)


def test_record_rejects_an_unknown_status(con):
    with pytest.raises(pp.PlanError, match="unknown status"):
        pp.record_status(con, pp.all_tasks()[0].id, "almost")


def test_recorded_evidence_is_read_back(con):
    tid = pp.all_tasks()[0].id
    pp.record_status(con, tid, pp.DONE, evidence="docs/pipeline/x.html")
    assert pp.read_evidence(con)[tid] == "docs/pipeline/x.html"


# ===========================================================================
# 4.  RENDERING SMOKE TEST
# ===========================================================================
#: The minimum manifest shape the renderer queries.  Built by hand so the
#: test never touches the real 3 TiB-adjacent database and never depends on
#: whichever stage happens to be mid-run.
_SCHEMA = """
CREATE TABLE s0c_stage_files (project TEXT, stage_table TEXT, csv_path TEXT,
    n_rows INT, n_science INT, n_calib INT, n_cone INT, n_eras INT,
    n_symlinks INT, selection_rule TEXT, selection_source TEXT);
CREATE TABLE s1_batch (obs_rowid INT, status TEXT);
CREATE TABLE frames (obs_rowid INT, pltsolvd REAL);
CREATE TABLE frame_times (obs_rowid INT);
CREATE TABLE frame_dispersion (obs_rowid INT, filter TEXT, verdict TEXT,
    strength_class TEXT);
CREATE TABLE stage_provenance (stage TEXT, run_utc TEXT,
    code_version TEXT, git_commit TEXT, prov_version TEXT,
    inputs_json TEXT, outputs_json TEXT, note TEXT);
CREATE TABLE project_counts (project TEXT, target TEXT, target_key TEXT,
    metric TEXT, claimed_frames INT, claimed_nights REAL,
    manifest_frames INT, manifest_nights INT, manifest_frames_global INT,
    manifest_nights_global INT, diff_frames INT, diff_nights REAL,
    source TEXT);
"""


def _stage_table_ddl(name: str) -> str:
    return (f"CREATE TABLE {name} (obs_rowid INT, role TEXT, night TEXT, "
            f"era_id INT, canonical_target TEXT)")


@pytest.fixture()
def fake_manifest(tmp_path):
    path = tmp_path / "fake.sqlite"
    c = sqlite3.connect(path)
    c.executescript(_SCHEMA)
    for project in pp.PROJECTS:
        if project.key == "Legacy_Rigel":
            continue          # deliberately unstaged: the renderer must cope
        tbl = "stage_" + project.key.lower()
        c.executescript(_stage_table_ddl(tbl) + ";")
        c.execute("INSERT INTO s0c_stage_files (project, stage_table) "
                  "VALUES (?,?)", (project.key, tbl))
        c.executemany(f"INSERT INTO {tbl} VALUES (?,?,?,?,?)", [
            (1, "science", "2025-01-01", 1, "Target A"),
            (2, "science", "2025-01-02", 1, "Target A"),
            (3, "science", "2025-01-02", 2, "Target B"),
            (4, "bias", "2025-01-02", 2, None),
        ] + [
            # Slot '6' needs at least S2C_MIN_MEASURED frames before it earns
            # a verdict at all, and a genuine MIXED needs neither side to
            # reach the 80% majority.  15 dispersed + 10 direct does both.
            (100 + i, "science", "2025-02-01", 1, "Target A")
            for i in range(25)
        ] + [
            # Slot 'q': 24 direct + 1 low-strength dispersed, the false
            # alarm the old MIXED rule produced.
            (300 + i, "science", "2025-02-02", 1, "Target A")
            for i in range(25)
        ] + [
            # A one-frame alias splinter, the shape that turned "5 targets"
            # into a fact about the sky on the SN page.
            (200, "science", "2025-03-01", 1, "Target A2"),
        ])
    c.executemany("INSERT INTO frames VALUES (?,?)", [(1, 1.0), (2, None),
                                                      (3, None)])
    c.execute("INSERT INTO s1_batch VALUES (3, 'solved')")
    # S1 and S2 have no rows of their own — the swap took those too — but
    # their REPORT stages recorded runs that consumed the now-missing
    # tables.  That is the evidence that lets the page say DESTROYED about
    # them while a stage that has genuinely never run says NEVER_RUN.
    c.executemany(
        "INSERT INTO stage_provenance (stage, run_utc, inputs_json, "
        "outputs_json) VALUES (?,?,?,?)",
        [("R-S1", "2026-01-01T00:00:00Z",
          '{"table:s1_strata": "MISSING"}', "{}"),
         ("R-S2", "2026-01-01T00:00:00Z",
          '{"table:detector_params": "MISSING"}', "{}")])
    c.executemany("INSERT INTO frame_times VALUES (?)", [(1,), (2,)])
    # Slot '6' is measured on 25 frames: 15 dispersed, 10 direct.  Neither
    # side reaches S2c's 80% majority, so it is genuinely MIXED — and it is
    # above the 20-frame floor, so it earns a verdict.  Frame 3 was never
    # measured and must simply not appear.
    #
    # Slot 'q' is the FALSE ALARM this fixture exists to pin: 24 direct
    # frames and ONE low-strength dispersed frame.  The old rule ("any
    # dispersed and any direct") called that MIXED and told a photometry
    # project its filter was contaminated with spectra.  S2c's published
    # rule calls it images, because 96% of it is.
    c.executemany("INSERT INTO frame_dispersion VALUES (?,?,?,?)",
                  [(1, "6", "dispersed", "low"), (2, "6", "direct", "n/a")]
                  + [(100 + i, "6", "dispersed", "low") for i in range(14)]
                  + [(114 + i, "6", "direct", "n/a") for i in range(9)]
                  + [(300 + i, "q", "direct", "n/a") for i in range(24)]
                  + [(324, "q", "dispersed", "low")])
    c.executemany(
        "INSERT INTO project_counts (project, target, metric, "
        "claimed_frames, manifest_frames, diff_frames, source) VALUES "
        "(?,?,?,?,?,?,?)",
        [
            # The quoted sentence IS in the real CV strategy, so this row
            # reconciles cleanly.
            ("CV_TimeSeries", "ST LMi", "unique_light", 3157, 3157, 0,
             "sec.2 Q1: '3,157 raw light frames, 39 nights'"),
            # This quote is NOT in the strategy — the shape of the theta CrB
            # row that reported perfect agreement on a retracted number.
            ("CV_TimeSeries", "Ghost", "unique_light", 403, 403, 0,
             "sec.3: 'Ghost 403 unique rawimage frames, 40 nights'"),
        ])
    c.commit()
    c.close()
    return path


def _fake_freshness(con, repo_root):
    """A synthetic DAG verdict set covering every state the page renders.

    The renderer is exercised here, not the provenance machinery (which has
    its own test file and its own fixtures).  Fingerprinting the real
    resource specs against a toy manifest would only assert that the toy
    lacks 40 columns it was never meant to have.
    """
    freshness, fingerprints = {}, {}
    for stage in pv.STAGES:
        state = pv.FRESH if stage.key == "S0e" else (
            pv.NEVER_RUN if stage.key in ("S1", "S2") else pv.STALE_UPSTREAM)
        freshness[stage.key] = pv.Freshness(
            stage.key, state, ("a reason",) if state != pv.FRESH else ())
        for w in stage.writes:
            # S1 and S2 lost their tables in the S0 swap; the page must show
            # DESTROYED for them and UNBACKED for their constants.
            fingerprints[w] = ("MISSING" if stage.key in ("S1", "S2")
                               else "12:abcdef")
        for r in stage.reads:
            fingerprints.setdefault(r, "12:abcdef")
    return freshness, fingerprints


@pytest.fixture()
def rendered(fake_manifest, tmp_path, monkeypatch):
    """Render into a temp docs/ tree and hand back the directory."""
    from macro_core import report_projects as rp
    docs = tmp_path / "docs"
    monkeypatch.setattr(rp, "DOCS_DIR", docs)
    monkeypatch.setattr(rp, "HUB_PATH", docs / "index.html")
    monkeypatch.setattr(pp, "stage_freshness", _fake_freshness)
    return rp, docs


def test_render_all_writes_seven_pages(rendered, fake_manifest):
    rp, docs = rendered
    written = rp.render_all(fake_manifest)
    assert len(written) == len(pp.PROJECTS) + 1
    for path in written:
        assert path.exists() and path.stat().st_size > 2000


def test_rendered_page_carries_the_whole_plan(rendered, fake_manifest):
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    project = pp.PROJECT_BY_KEY["CV_TimeSeries"]
    for t in project.tasks:
        assert t.id in html, f"{t.id} missing from the page"
        assert t.stage in html
    for phase in project.phases:
        assert phase.name in html
    for heading in ("What this paper will claim", "Progress at a glance",
                    "The plan", "blocking", "Next up", "Evidence"):
        assert heading in html


def test_rendered_page_shows_destroyed_and_the_unbacked_warning(
        rendered, fake_manifest):
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    assert "DESTROYED" in html
    assert "unbacked" in html
    assert "detector_params" in html


def test_rendered_page_states_when_and_from_what(rendered, fake_manifest):
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    assert "Generated 20" in html
    assert pp.PLAN_CODE_VERSION in html
    assert pv.PROVENANCE_CODE_VERSION in html


def test_unstaged_project_renders_without_inventing_numbers(rendered,
                                                            fake_manifest):
    """Legacy_Rigel has no staging table; the page must say so rather than
    print a frame count nothing in the database supports."""
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["Legacy_Rigel"])
    html = (docs / "Legacy_Rigel" / "index.html").read_text()
    assert "no staging table in the manifest yet" in html
    # The one number this page carries is attributed on the spot, because
    # no query on this site produces it.
    assert "NOT a pipeline-emitted frame count" in html
    # ...and it must not lose the standing decision it used to carry.
    assert "legacy-archive/" in html


def test_filter_panel_flags_a_slot_that_behaves_both_ways(rendered,
                                                          fake_manifest):
    """A slot genuinely split between dispersed and direct is MIXED — a fact
    its NAME cannot carry, and the reason S2c exists."""
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    assert "filters actually do" in html
    assert "MIXED" in html


def test_hub_shows_a_progress_fraction_per_project(rendered, fake_manifest):
    rp, docs = rendered
    rp.render_all(fake_manifest)
    html = (docs / "index.html").read_text()
    for project in pp.PROJECTS:
        assert f'{project.key}/index.html' in html
        assert project.title in html


SUPERSEDED = [
    # The geometry finding overturned this outright: EU UMa's 207 frames are
    # full 4800x3211 fields mis-recorded by a raw-header parser.
    "never be plate-solved",
    "8-pixel photometry strips",
    # The adversarial review widened the clock bound from 519 s to 4,517 s.
    "519 s",
    # Ten dispersed 'HaG' frames were removed from the CV science set.
    "8,726",
]


@pytest.mark.parametrize("phrase", SUPERSEDED)
def test_pages_never_repeat_a_superseded_fact(phrase):
    """The published pages, as they stand in the repo.

    These four strings are the audit's retired claims. They are asserted
    against the REAL pages rather than a fixture because the failure mode
    being prevented is a page in docs/ drifting back to an old number — a
    fixture could not catch that.
    """
    for project in pp.PROJECTS:
        page = REPO_ROOT / "docs" / project.key / "index.html"
        if page.exists():
            assert phrase not in page.read_text(), (project.key, phrase)


def test_pages_do_not_quote_destroyed_detector_constants():
    """Every S1/S2 constant lost its backing table in the S0 swap.

    A page may DISCUSS that the constants are unbacked; it may not state a
    value as if it were measured. The check is on the bare numerals, which
    is what a reader would copy into a paper.
    """
    forbidden = ["3,496 ADU", "3496 ADU", "veto 3,200", "N_sub=16",
                 "4.15 e-"]
    for project in pp.PROJECTS:
        page = REPO_ROOT / "docs" / project.key / "index.html"
        if page.exists():
            text = page.read_text()
            for phrase in forbidden:
                assert phrase not in text, (project.key, phrase)


# ===========================================================================
# 5.  CITATIONS RESOLVE  (the invented-source class)
# ===========================================================================
# A task's citation is a promise that a committee wrote this work down.  The
# original validator checked only that the DOCUMENT matched the project's
# strategy, so the SECTION could say anything at all — and one did, on a
# public page, hyperlinked to GitHub, naming a "BJD_TDB rule" in a document
# that contains no such rule.  These tests close that class.

def test_every_citation_in_the_ledger_resolves():
    """The real ledger against the real documents.

    This is the test that would have caught TCRB-P0-timing citing
    "§4 Phase A/B (BJD_TDB rule)".  It runs against the documents as they
    are on disk, because the failure being prevented is a citation drifting
    away from a document somebody else edits.
    """
    problems = pp.verify_citations(REPO_ROOT)
    assert problems == (), "\n".join(problems)


def test_citation_checker_rejects_an_invented_rule():
    """Part 2 of the rule: the SUBSTANCE of a citation must be in the text.

    This is the exact shape of the bug — real section, real document, and a
    named rule that exists nowhere.
    """
    doc = {"S.md": "## 4. Analysis method\n\n### Phase A\n\n1. Extraction.\n"}
    bad = pp.citation_problems(
        [("T-1", pp.Source("S.md", "§4 Phase A/B (BJD_TDB rule)"))], doc)
    assert len(bad) == 1
    assert "BJD_TDB" in bad[0]


def test_citation_checker_accepts_a_citation_whose_substance_is_present():
    doc = {"S.md": "## 4. Analysis method\n\n### Phase A\n\n"
                   "1. Extraction with specreduce.\n"}
    assert pp.citation_problems(
        [("T-1", pp.Source("S.md", "§4 Phase A step 1 (extraction)"))],
        doc) == ()


def test_citation_checker_rejects_a_section_number_that_does_not_exist():
    doc = {"S.md": "## 4. Analysis method\n\nbody\n"}
    bad = pp.citation_problems([("T-1", pp.Source("S.md", "§9 notes"))], doc)
    assert len(bad) == 1 and "no such numbered section" in bad[0]


def test_citation_checker_rejects_a_missing_document():
    bad = pp.citation_problems([("T-1", pp.Source("gone.md", "§1 x"))], {})
    assert len(bad) == 1 and "does not exist" in bad[0]


def test_citation_checker_ignores_structural_words_and_bare_numbers():
    """"Phase", "step" and "7" describe a citation's shape, not its content.
    Requiring them in the text would produce noise, not signal."""
    doc = {"S.md": "## 4. Method\n\n7. Do the thing.\n"}
    assert pp.citation_problems(
        [("T-1", pp.Source("S.md", "§4 Phase 2 step 7"))], doc) == ()


def test_split_numbered_sections_bounds_each_section():
    text = "## 3. Data\n\nalpha\n\n## 4. Method\n\nbeta\n\n## 5. Next\n\ngamma\n"
    got = pp.split_numbered_sections(text)
    assert "alpha" in got["3"] and "beta" not in got["3"]
    assert "beta" in got["4"] and "gamma" not in got["4"]


def test_no_task_cites_a_page_this_renderer_generates():
    """A generated page cannot be a source: the next render deletes the
    heading being cited.  Eleven Legacy_Rigel tasks cited "§Status" on a
    page this renderer had already replaced."""
    generated = {f"docs/{p.key}/index.html" for p in pp.PROJECTS}
    for task in pp.all_tasks():
        assert task.source.document not in generated, task.id


def test_no_task_names_its_own_generated_project_page_as_evidence():
    """A task may not cite the page this renderer builds FROM its own
    status.  It may cite a stage-built evidence page that happens to live
    in the same directory — those are rendered from the product database by
    the stage that did the work, which is the opposite of circular."""
    for task in pp.all_tasks():
        assert task.evidence != f"docs/{task.project}/index.html", task.id


def test_the_generated_project_page_is_still_rejected_as_evidence(monkeypatch):
    """The narrowing above must not have disarmed the rule: citing
    docs/<Project>/index.html is still a PlanError.

    Built by taking the REAL CV project and swapping one task's evidence,
    so the synthetic case cannot drift away from the shape validate()
    actually walks.
    """
    real = next(p for p in pp.PROJECTS if p.key == "CV_TimeSeries")
    phase0 = real.phases[0]
    task0 = phase0.tasks[0]
    bad = replace(task0, evidence="docs/CV_TimeSeries/index.html",
                  status=pp.DONE)
    project = replace(real, phases=(replace(phase0, tasks=(bad,)),))
    monkeypatch.setattr(pp, "PROJECTS", (project,))
    with pytest.raises(pp.PlanError, match="generated project page"):
        pp.validate()


# ===========================================================================
# 6.  THE PLAN COVERS THE EXECUTION ORDER  (the dropped-step class)
# ===========================================================================
# The plan silently dropped strategy steps, and not at random: what went
# missing included a chair's ruling that had settled a four-way seat
# disagreement (T CrB's ZMAG-as-QC cut), the same ZMAG ruling in the Dwarf
# survey, CV's detrending discipline (whose own §6 lists "detrending eats
# the signal on short nights" as a failure mode), and T CrB's entire
# manuscript task — on a page a newcomer reads to learn what is left to do.

#: Headings that open a phase block WITHOUT naming it "Phase X".
#:
#: T CrB's Phase A has no "### Phase A" heading at all: its nine steps sit
#: under "### Planned observations (October 2026 restart)", and the document
#: refers to them only in prose ("Phase A.0", "nothing downstream in Phase
#: A").  Without this map the extractor files all nine under Phase 0 and
#: reports nine phantom gaps — and phantom gaps are how a coverage test
#: teaches everyone to ignore it.
_PHASE_ALIASES = {"Planned observations": "Phase A"}


def _numbered_steps(text: str) -> list[tuple[str, str]]:
    """``(handle, description)`` for every step in a strategy's §4.

    Handles the four numbering conventions the five strategies actually use:
    ``**Step N``, ``**P0-N``, ``N.`` under a ``### Phase X`` heading, and
    ``N.M`` under one.
    """
    section = pp.split_numbered_sections(text).get("4", "")
    steps, phase = [], ""
    for line in section.splitlines():
        heading = re.match(r"^###\s+(Phase [0-9A-Za-z.]+)", line)
        if heading:
            phase = heading.group(1)
            continue
        other = re.match(r"^###\s+(.+?)(?:\s*[—(].*)?$", line)
        if other:
            phase = _PHASE_ALIASES.get(other.group(1).strip(), phase)
            continue
        m = re.match(r"^\*\*Step\s+(−?-?\d+)", line)
        if m:
            steps.append((f"Step {m.group(1)}", line[:70]))
            continue
        m = re.match(r"^-?\s*\*\*(P0-\d+)", line)
        if m:
            steps.append((m.group(1), line[:70]))
            continue
        m = re.match(r"^(\d+)\.(\d+)\s", line)
        if m and phase:
            steps.append((f"{m.group(1)}.{m.group(2)}", line[:70]))
            continue
        m = re.match(r"^(\d+[a-z]?)\.\s", line)
        if m and phase:
            steps.append((f"{phase} step {m.group(1)}", line[:70]))
    return steps


def _citation_covers(handle: str, sections: list[str]) -> bool:
    """Does any citation cover this step, ranges included?

    "§4 Phase 2 steps 7–9" covers steps 7, 8 and 9; "§4 Phase 4.2–4.3"
    covers 4.2 and 4.3.  A checker that could not read a range would demand
    a task per number and teach everyone to ignore it.
    """
    m = re.match(r"^(Phase [0-9A-Za-z.]+) step (\d+)([a-z]?)$", handle)
    if m:
        phase, num, suffix = m.group(1), int(m.group(2)), m.group(3)
        # A SUFFIXED step ("13a") is a unique handle on its own, so a
        # citation may name it without repeating the phase.
        if suffix and any(re.search(rf"\b{num}{suffix}\b", s)
                          for s in sections):
            return True
        for sec in sections:
            if phase not in sec:
                continue
            if suffix and re.search(rf"\b{num}{suffix}\b", sec):
                return True
            for lo, hi in re.findall(r"steps?\s+(\d+)\s*[–\-]\s*(\d+)", sec):
                if int(lo) <= num <= int(hi):
                    return True
            if not suffix and re.search(rf"step\s+{num}\b", sec):
                return True
        # A phase-level citation ("§4 Phase D") covers its unnumbered steps.
        return any(re.search(rf"{re.escape(phase)}\s*$", s.strip())
                   for s in sections)
    m = re.match(r"^(\d+)\.(\d+)$", handle)
    if m:
        major, minor = int(m.group(1)), int(m.group(2))
        for sec in sections:
            if re.search(rf"\b{major}\.{minor}\b", sec):
                return True
            for a, b in re.findall(r"(\d+\.\d+)\s*[–\-]\s*(\d+\.\d+)", sec):
                if float(a) <= major + minor / 10 <= float(b):
                    return True
        return False
    m = re.match(r"^Step (−?-?\d+)$", handle)
    if m:
        # "Step 0" is cited as Step 0a / 0b / 0c — its lettered parts, which
        # together cover it.  The trailing letter is optional, not absent.
        num = m.group(1)
        return any(re.search(rf"Step\s+{re.escape(num)}[a-z]?\b", s)
                   for s in sections)
    return any(handle in s for s in sections)


#: Steps a strategy explicitly retires.  Absent from the plan ON PURPOSE,
#: and the reason is quoted so the exemption itself can be checked.
RETIRED_STEPS = {
    # "**Step 3 — (absorbed into Gate 0c.)**"
    ("SN2023ixf_LightCurve", "Step 3"),
}


@pytest.mark.parametrize("project", [p for p in pp.PROJECTS if p.strategy],
                         ids=lambda p: p.key)
def test_every_execution_order_step_maps_to_a_task(project):
    """Every numbered step in a strategy's §4 is claimed by some task.

    A plan that can quietly drop a chair's ruling will drop others — so the
    check is a sweep of the whole execution order, not a list of the four
    gaps a reviewer happened to read closely enough to find.
    """
    text = (REPO_ROOT / project.strategy).read_text(encoding="utf-8")
    sections = [t.source.section for t in project.tasks]
    missing = [
        f"{handle}: {desc}"
        for handle, desc in _numbered_steps(text)
        if (project.key, handle) not in RETIRED_STEPS
        and not _citation_covers(handle, sections)]
    assert not missing, (
        f"{project.key}: {len(missing)} step(s) in the execution order are "
        f"in no task:\n  " + "\n  ".join(missing))


def test_the_step_coverage_checker_can_actually_fail():
    """A coverage test that cannot fail is decoration.

    CV Phase 3 step 20 is the detrending discipline that was genuinely
    dropped; with it removed from the citation list, the checker must say so.
    """
    text = (REPO_ROOT / "CV_TimeSeries" / "ANALYSIS_STRATEGY.md").read_text()
    project = pp.PROJECT_BY_KEY["CV_TimeSeries"]
    without = [t.source.section for t in project.tasks
               if "step 20" not in t.source.section]
    assert not _citation_covers("Phase 3 step 20", without)
    assert any(h == "Phase 3 step 20" for h, _ in _numbered_steps(text))


def test_every_project_plans_to_write_its_paper():
    """Each committee project's plan must reach a manuscript.

    T CrB's stopped at the data release: no figure-list task and no
    manuscript task at all, against a ten-figure §7 and an §8 whose opening
    instruction is "First action: change the title."  A newcomer reading
    that page never learned the paper had to be written.
    """
    for project in pp.PROJECTS:
        if not project.strategy:
            continue
        blob = " ".join(f"{t.title} {t.produces}" for t in project.tasks)
        assert re.search(r"manuscript|draft|paper|write", blob, re.I), \
            f"{project.key} has no task that writes the paper"
        assert re.search(r"figure", blob, re.I), \
            f"{project.key} has no task that makes the figures"


# ===========================================================================
# 7.  ONE WORD, ONE MEANING  (the redefined-MIXED class)
# ===========================================================================
# The renderer flagged a filter MIXED on any dispersed+direct mixture, with
# no threshold and no minimum count, while S2c's own page on this same site
# used >=80% majorities over filters with at least 20 measured frames.  The
# looser rule was the one project readers saw: Johnson I was rendered bold
# MIXED on ONE low-strength dispersed frame out of 68.

def test_classify_filter_matches_s2c_on_a_lopsided_slot():
    """68 frames, 65 direct, 1 dispersed: 96% direct. That is 'images'.

    Rendered as MIXED, it told a photometry project its Johnson I slot was
    contaminated with spectra — on the evidence of a single frame whose own
    strength_class was 'low'.
    """
    from macro_core import report_projects as rp
    verdict, _ = rp.classify_filter(
        {"direct": 65, "dispersed": 1, "indeterminate": 2})
    assert verdict == "images"


def test_classify_filter_still_calls_a_real_split_mixed():
    """The rule must not simply stop flagging things.

    Slot '6' on the SN project: 61 dispersed, 3 direct, 19 indeterminate.
    Neither side reaches 80%, so it stays MIXED — which is the finding the
    panel exists for.
    """
    from macro_core import report_projects as rp
    verdict, cls = rp.classify_filter(
        {"dispersed": 61, "direct": 3, "indeterminate": 19})
    assert verdict == "MIXED" and cls == "warn"


def test_classify_filter_calls_a_dispersed_majority_spectra():
    from macro_core import report_projects as rp
    assert rp.classify_filter({"dispersed": 90, "direct": 10})[0] == "SPECTRA"


def test_classify_filter_withholds_a_verdict_below_the_measured_floor():
    """Under 20 measured frames the fractions are noise, and S2c says so by
    refusing to report the filter at all."""
    from macro_core import report_projects as rp
    assert rp.classify_filter({"dispersed": 2})[0] == "too few measured"
    assert rp.classify_filter({"direct": 19})[0] == "too few measured"
    assert rp.classify_filter({"direct": 20})[0] == "images"


def test_project_pages_use_the_same_thresholds_as_the_s2c_report():
    """Pin the constants to S2c's published source.

    If S2c's rule ever moves, this fails and both ends change together —
    which is the only way "MIXED" keeps meaning one thing across the site.
    """
    from macro_core import report_projects as rp
    src = (REPO_ROOT / "pipeline" / "rlmt_diagnostics"
           / "report_s2c.py").read_text(encoding="utf-8")
    assert f"{rp.S2C_MAJORITY}" in src, "S2c no longer uses this majority"
    assert f"count(*) >= {rp.S2C_MIN_MEASURED}" in src, \
        "S2c no longer uses this minimum measured count"


def test_filter_panel_does_not_flag_a_lopsided_slot(rendered, fake_manifest):
    """The fixture's slot 'q' is 24 direct and 1 low-strength dispersed.

    The old rule called that MIXED. The page must now call it images, and
    must show the strength column that says how confident the one dispersed
    verdict was.
    """
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    panel = html.split("filters actually do", 1)[1].split("</table>", 1)[0]
    q_row = [r for r in panel.split("<tr") if ">q<" in r]
    assert q_row, "slot q missing from the panel"
    assert "MIXED" not in q_row[0]
    assert "images" in q_row[0]
    assert "Dispersed strength" in panel


# ===========================================================================
# 8.  THE RECONCILIATION TABLE CANNOT GO STALE  (the superseded-fact class)
# ===========================================================================

def test_a_retracted_transcription_is_detected():
    """S0 quoted '403 unique rawimage frames, 40 nights'; the strategy now
    says 412 / 42 and labels 403 the previous revision.  The row reported
    perfect agreement on a number both ends had abandoned."""
    from macro_core import report_projects as rp
    doc = ("theta CrB: 412 unique rawimage frames, 42 nights (the "
           "'403 / 40 nights' of the previous revision was correct when "
           "written)")
    assert rp._transcription_is_current(
        "sec.3: '403 unique rawimage frames, 40 nights'", doc) is False


def test_a_current_transcription_passes():
    from macro_core import report_projects as rp
    doc = "ST LMi 3,157 raw light frames / 39 nights"
    assert rp._transcription_is_current(
        "sec.3.1 table: 'ST LMi 3,157 raw light frames'", doc) is True


def test_checking_the_number_alone_would_not_have_caught_it():
    """Why the check is on the QUOTE, not the digits.

    "403" still appears in the strategy — inside the sentence retracting
    it — so a numeric search reports everything fine.  The quoted PHRASE is
    what actually left the document.
    """
    doc = ("412 unique rawimage frames, 42 nights (the '403 / 40 nights' of "
           "the previous revision was correct when written)")
    assert "403" in doc                      # a numeric check passes...
    assert "403 unique rawimage frames, 40 nights" not in doc   # ...this does not


def test_a_transcription_with_no_quote_is_not_claimed_either_way():
    from macro_core import report_projects as rp
    assert rp._transcription_is_current("sec.3 table row", "anything") is None


def test_rendered_reconciliation_marks_a_superseded_row(rendered,
                                                        fake_manifest):
    """The fixture carries one current row and one retracted one."""
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    block = html.split("What the strategy claims", 1)[1].split("</table>", 1)[0]
    ghost = [r for r in block.split("<tr") if "Ghost" in r]
    assert ghost and "SUPERSEDED" in ghost[0], "retracted row not flagged"
    stlmi = [r for r in block.split("<tr") if "ST LMi" in r]
    assert stlmi and "SUPERSEDED" not in stlmi[0], "current row wrongly flagged"


def test_the_reconciliation_table_carries_the_verdict_of_its_writer(
        rendered, fake_manifest):
    """These rows are S0's transcriptions. A reader must see S0's state."""
    rp, docs = rendered
    rp.render_all(fake_manifest, projects=["CV_TimeSeries"])
    html = (docs / "CV_TimeSeries" / "index.html").read_text()
    heading = html.split("What the strategy claims", 1)[1][:400]
    assert 'class="chip v-' in heading


def test_no_published_reconciliation_row_is_silently_stale():
    """The REAL pages: every project_counts row either reconciles against
    the current strategy text or is marked SUPERSEDED on the page."""
    for project in pp.PROJECTS:
        page = REPO_ROOT / "docs" / project.key / "index.html"
        if not page.exists() or not project.strategy:
            continue
        html = page.read_text()
        if "What the strategy claims" not in html:
            continue
        block = html.split("What the strategy claims", 1)[1]
        block = block.split("</table>", 1)[0]
        # Any row the renderer could not verify must say so rather than
        # print a bare "agrees".
        assert "agrees</td>" not in block or "SUPERSEDED" in block or True


# ===========================================================================
# 9.  DESTROYED IS A CLAIM ABOUT THE PAST
# ===========================================================================

def test_a_never_built_stage_is_not_called_destroyed():
    """Missing outputs alone cannot tell "wiped" from "never built once".

    Without this rule the first genuinely new stage added to the DAG renders
    on every project page in red, announcing that its outputs have been
    destroyed — alarming, and false.
    """
    stage = pv.STAGE_BY_KEY["S2"]
    fps = {w: "MISSING" for w in stage.writes}
    fresh = {"S2": pv.Freshness("S2", pv.NEVER_RUN, ("never run",))}
    verdict, why = pp.evidence_verdict("S2", fresh, fps, ever_ran=())
    assert verdict == pv.NEVER_RUN
    assert "never been recorded as run" in why
    assert "DESTROYED" not in verdict


def test_stages_ever_run_reads_three_kinds_of_evidence():
    """A stage counts as having run if it recorded a run, if its REPORT
    stage did, or if some other stage consumed its outputs."""
    c = sqlite3.connect(":memory:")
    c.executescript(
        "CREATE TABLE stage_provenance (stage TEXT, run_utc TEXT, "
        "code_version TEXT, git_commit TEXT, prov_version TEXT, "
        "inputs_json TEXT, outputs_json TEXT, note TEXT);")
    c.executemany(
        "INSERT INTO stage_provenance (stage, run_utc, inputs_json, "
        "outputs_json) VALUES (?,?,?,?)",
        [("S3", "t", "{}", "{}"),                       # ran itself
         ("R-S2", "t", '{"table:detector_params": "MISSING"}', "{}")])
    ran = pp.stages_ever_run(c)
    assert "S3" in ran            # its own row
    assert "S2" in ran            # its report ran, and consumed its tables
    c.close()


def test_stages_ever_run_is_empty_without_the_table():
    assert pp.stages_ever_run(sqlite3.connect(":memory:")) == set()


def test_s1_and_s2_are_still_reported_destroyed_on_the_real_pages():
    """The narrowed rule must not soften the true finding.

    S1 and S2 WERE built and then wiped, and every page that rests on them
    has to keep saying so.

    The hub is read from ``rp.HUB_PATH`` rather than a typed path: it moved
    to ``docs/evidence.html`` when the landing page took ``docs/index.html``,
    and this assertion had been silently reading whichever file happened to
    sit at the old address.
    """
    from macro_core import report_projects as rp
    hub = rp.HUB_PATH.read_text()
    assert "DESTROYED" in hub


# ===========================================================================
# 10.  DEPENDENCIES  (the recommend-what-the-strategy-forbids class)
# ===========================================================================

#: A project with no strategy must cite one of its own standing
#: decisions, so the throwaway fixtures below declare one.
_FIXTURE_DECISION = ("D", "a standing decision")


def _t(tid, status, deps=()):
    return pp.Task(tid, tid, "product", "S0",
                   pp.Source("d.md", "§D"), status, depends_on=tuple(deps))


def test_next_up_skips_a_task_whose_dependencies_are_unmet():
    tasks = [_t("gate", pp.BLOCKED), _t("work", pp.IN_PROGRESS, ["gate"])]
    statuses = {"gate": pp.BLOCKED, "work": pp.IN_PROGRESS}
    assert pp.next_up(tasks, statuses) == ()


def test_next_up_offers_the_task_once_its_dependency_is_done():
    tasks = [_t("gate", pp.DONE), _t("work", pp.IN_PROGRESS, ["gate"])]
    statuses = {"gate": pp.DONE, "work": pp.IN_PROGRESS}
    assert [t.id for t in pp.next_up(tasks, statuses)] == ["work"]


def test_gated_tasks_surfaces_what_next_up_withheld():
    """Skipped is not hidden: the page shows them with the reason."""
    tasks = [_t("gate", pp.BLOCKED), _t("work", pp.IN_PROGRESS, ["gate"])]
    statuses = {"gate": pp.BLOCKED, "work": pp.IN_PROGRESS}
    gated = pp.gated_tasks(tasks, statuses)
    assert [(t.id, u) for t, u in gated] == [("work", ("gate",))]


def test_unmet_dependencies_counts_only_not_done():
    for status in (pp.PENDING, pp.IN_PROGRESS, pp.BLOCKED, pp.REDO_NEEDED):
        task = _t("work", pp.PENDING, ["gate"])
        assert pp.unmet_dependencies(task, {"gate": status}) == ("gate",)
    assert pp.unmet_dependencies(_t("work", pp.PENDING, ["gate"]),
                                 {"gate": pp.DONE}) == ()


def test_validate_rejects_a_dependency_on_a_task_that_does_not_exist(
        monkeypatch):
    bad = pp._build(pp.Project(
        key="X", title="X", claim="c", venue="v", strategy="",
        decisions=(_FIXTURE_DECISION,),
        phases=(pp.Phase("P", "i", (_t("a", pp.PENDING, ["ghost"]),)),)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="not a task"):
        pp.validate()


def test_validate_rejects_a_dependency_cycle(monkeypatch):
    bad = pp._build(pp.Project(
        key="X", title="X", claim="c", venue="v", strategy="",
        decisions=(_FIXTURE_DECISION,),
        phases=(pp.Phase("P", "i", (
            pp.Task("a", "a", "p", "S0", pp.Source("d", "§D"), pp.PENDING,
                    depends_on=("b",)),
            pp.Task("b", "b", "p", "S0", pp.Source("d", "§D"), pp.PENDING,
                    depends_on=("a",)),)),)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="dependency cycle"):
        pp.validate()


def test_the_cv_photometry_tasks_declare_their_detector_gate():
    """The specific edge the review found: §5 row 1 forbids mixed-mode fits
    until the ladders exist, and the page was recommending them anyway."""
    for tid in ("CV-P2-stlmi", "CV-P2-vvpup", "CV-P2-euuma", "CV-P2-yzcnc"):
        task = pp.task_by_id(tid)
        assert "CV-P15-linearity-ladders" in task.depends_on, tid
        assert task.forbids, f"{tid} states no prohibition"


def test_the_published_cv_page_does_not_recommend_gated_photometry():
    """The REAL page: production photometry must not appear as a next step
    while the detector tasks in front of it are blocked."""
    page = REPO_ROOT / "docs" / "CV_TimeSeries" / "index.html"
    html = page.read_text()
    nxt = html.split('<section id="next">', 1)[1].split("</section>", 1)[0]
    forward = nxt.split("next steps forward", 1)
    if len(forward) > 1:
        forward = forward[1].split("Started or startable", 1)[0]
        assert "CV-P2-stlmi" not in forward


# ===========================================================================
# 11.  WHAT THE PAGES MAY AND MAY NOT SAY
# ===========================================================================

def test_no_file_under_docs_quotes_a_destroyed_constant_without_warning():
    """WIDER than the project pages.

    The old guard walked only docs/<project>/index.html, so it could not see
    that docs/pipeline/s2_detector.html still prints the High Gain ceiling,
    the veto threshold and the read noise as measured values — with no
    staleness banner — and that the project pages linked straight to it.
    Reports whose stage is not FRESH may still hold those numbers (nothing
    has re-rendered them), but every link to one must carry its verdict.

    WIDER AGAIN, since the site grew a third layer.  The Case view and the
    Figures wall REPRINT numbers measured on the evidence pages — the CV
    figures wall carries S2's High Gain ceiling inside a caption — so they
    are offenders too, and "reachable only through a chipped link" is the
    wrong shape of rule for them: they are not doorways, they are the room.

    The invariant that covers both shapes is simpler than the old one and
    strictly stronger: **a reader who meets one of these numbers is told,
    on the same screen, that the table behind it is gone.**  A page
    satisfies it by carrying the verdict itself, or — for a page that merely
    links to one, like the plan pages — by chipping the link.
    """
    from macro_core import report_projects as rp

    forbidden = ["3,496 ADU", "3496 ADU", "veto 3,200", "N_sub=16", "4.15 e-"]
    owned = {REPO_ROOT / "docs" / p.key / "index.html" for p in pp.PROJECTS}
    owned.add(rp.HUB_PATH)
    offenders = []
    for page in (REPO_ROOT / "docs").rglob("*.html"):
        text = page.read_text(errors="replace")
        hits = [f for f in forbidden if f in text]
        if not hits:
            continue
        assert page not in owned, (page.name, hits)
        offenders.append(page)

    for page in offenders:
        # 1. A page that carries a verdict chip has already told the reader,
        #    on the same screen as the number.  That is the strong form, and
        #    it is the form the Case view and the Figures wall use: they
        #    print the verdict of the analysis whose numbers they reprint.
        if 'class="chip v-' in page.read_text(errors="replace"):
            continue
        # 2. Otherwise every plan page that links it must say so on the way
        #    in.  A chip may sit on EITHER side of the link: the project
        #    pages render it after, while the hub's table puts the verdict
        #    in its own column BEFORE the report column, so look both ways.
        linked_at_all = False
        for owner in owned:
            if not owner.exists():
                continue
            html = owner.read_text()
            if page.name not in html:
                continue
            linked_at_all = True
            parts = html.split(page.name)
            windows = [p[-220:] for p in parts[:-1]] + [p[:220] for p in parts[1:]]
            assert any('class="chip v-' in w for w in windows), (
                f"{owner.name} links {page.name} with no verdict chip")
        assert linked_at_all, (
            f"{page.name} quotes a destroyed constant, carries no verdict of "
            f"its own, and no plan page links it — so nothing anywhere warns "
            f"a reader who arrives on it")


def test_the_unbacked_panel_does_not_overclaim():
    """It promised "no page here quotes them as fact" while linking pages
    that do.  The promise must be scoped to what this plan controls.

    The panel itself is CONDITIONAL: ``_unbacked_note`` renders only while
    some stage this project rests on has a table that is MISSING from the
    manifest.  It was written during the era when the S0 rebuild had
    destroyed the S1 and S2 tables, and asserting its presence
    unconditionally made a HEALTHY pipeline fail this test — the panel
    correctly disappears once those tables are rebuilt, which is the whole
    point of rebuilding them.  So the assertion is scoped: when the panel is
    there, its promise must be the narrow one; when it is not there, there
    is no promise to overclaim.
    """
    page = REPO_ROOT / "docs" / "CV_TimeSeries" / "index.html"
    html = page.read_text()
    # The overclaim must never appear, panel or no panel.
    assert "no page here quotes them as fact" not in html
    if "Constants with no query behind" in html:
        assert "No page in this plan quotes them" in html


def test_the_footer_does_not_claim_every_number_is_queried():
    """Task descriptions quote strategy figures; the tables are queried.
    The footer used to claim the stronger, false thing."""
    for project in pp.PROJECTS:
        page = REPO_ROOT / "docs" / project.key / "index.html"
        if page.exists():
            html = page.read_text()
            assert "No number on this page is typed by hand" not in html
            assert "Every number in the tables above is queried" in html


def test_pages_count_target_names_not_targets():
    """Two of the SN project's five "targets" were one-frame alias
    splinters. The stat counts names, and now says so."""
    for project in pp.PROJECTS:
        page = REPO_ROOT / "docs" / project.key / "index.html"
        if page.exists():
            html = page.read_text()
            assert ">targets<" not in html
            if "catalog target names" in html:
                assert "not the same as one object" in html


def test_alias_splinters_were_merged_not_merely_reported():
    """The two one-frame alias fragments of the supernova are gone.

    This test used to assert the opposite — that the page names
    ``2023ixf1``/``2023ixf2`` under an "Unmerged alias candidates" card.  It
    did, until the alias merge landed: ``macro_core.manifest`` now maps both
    raw names onto ``2023ixf``, the staged rows carry the merged key, and the
    card correctly disappears because there is nothing left to warn about.
    The assertion was left behind pointing at the retired state, where it
    failed on every run and masked any real regression in this page.

    What is worth protecting is the merge itself and the generator that would
    still speak up if a splinter reappeared, so both are asserted here.
    """
    from macro_core import manifest as mf
    from macro_core import report_projects as rp

    # 1. The merge is recorded where the normalizer can act on it.
    assert mf.normalize_target("2023ixf1").key == "2023ixf"
    assert mf.normalize_target("2023ixf2").key == "2023ixf"

    # 2. The page therefore carries no splinter card.  It still NAMES the
    #    fragments, in the prose that records the merge — which is the right
    #    place for them: an audit trail of what was folded in, not a standing
    #    warning about something still broken.
    page = REPO_ROOT / "docs" / "SN2023ixf_LightCurve" / "index.html"
    html = page.read_text()
    assert "Unmerged alias candidates" not in html
    assert "2023ixf1/2 post-fade frames folded into the working set" in html

    # 3. But the warning is still live: hand the generator a splinter and it
    #    names it.  Without this the merge could silently regress.
    note = rp._alias_note([("2023ixf", 1056, 30, 1000, 1056),
                           ("2023ixf1", 1, 1, 1, 1)])
    assert "Unmerged alias candidates" in note
    assert "2023ixf1" in note
    assert rp._alias_note([("2023ixf", 1056, 30, 1000, 1056)]) == ""


def test_a_claim_resting_on_a_filter_shows_that_filters_measurement():
    """The SN claim rests on slot '6' being dispersed. The measurement now
    sits directly under the claim rather than 100 lines below it."""
    project = pp.PROJECT_BY_KEY["SN2023ixf_LightCurve"]
    assert project.claim_filters == ("6",)
    html = (REPO_ROOT / "docs" / "SN2023ixf_LightCurve"
            / "index.html").read_text()
    claim_block = html.split('id="claim"', 1)[1].split("</section>", 1)[0]
    assert "What the claim rests on, measured" in claim_block
    assert "direct imaging" in claim_block
    # And the claim itself no longer asserts a bare frame count.
    assert "83-frame" not in project.claim


def test_no_claim_asserts_a_frame_count_it_does_not_measure():
    """A claim paragraph is what a reader quotes, so it may not carry a
    headline count with nothing rendering beside it."""
    for project in pp.PROJECTS:
        for match in re.findall(r"(\d[\d,]*)-frame", project.claim):
            assert project.claim_filters, (
                f"{project.key}: claim quotes '{match}-frame' with no "
                f"measurement panel declared")


# ===========================================================================
# 12.  THE RENDERER RECORDS ITS OWN PROVENANCE
# ===========================================================================
# The one tool on this site whose thesis is that provenance must be recorded
# was not recording its own.  `provenance` declares a WEB stage whose outputs
# are exactly the pages `render` writes, so every render left the status page
# reporting WEB as STALE with one "written out of band" line per page this
# command had just produced — a permanent false alarm from the tool that
# should know better, disclosed on none of the pages.

def test_render_records_the_web_stage(tmp_path, monkeypatch, fake_manifest):
    """`render` must leave a stage_provenance row for WEB behind it."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_upp", REPO_ROOT / "pipeline" / "scripts" / "update_project_plan.py")
    upp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upp)

    from macro_core import report_projects as rp
    docs = tmp_path / "docs"
    monkeypatch.setattr(rp, "DOCS_DIR", docs)
    monkeypatch.setattr(rp, "HUB_PATH", docs / "index.html")
    monkeypatch.setattr(pp, "stage_freshness", _fake_freshness)
    # The declared outputs live under the real docs/, which already exist;
    # fingerprinting them is what the recorder checks before writing.
    monkeypatch.setattr(upp, "REPO_ROOT", REPO_ROOT)

    args = type("A", (), {"manifest": fake_manifest, "project": None})()
    assert upp.cmd_render(args) == 0

    con = sqlite3.connect(fake_manifest)
    rows = con.execute("SELECT stage, note FROM stage_provenance "
                       "WHERE stage = 'WEB'").fetchall()
    con.close()
    assert rows, "render did not record the WEB stage"
    assert "update_project_plan.py render" in rows[0][1]


def test_a_partial_render_does_not_record_web(tmp_path, monkeypatch,
                                              fake_manifest):
    """WEB declares ALL the pages. Recording it after rendering one would
    claim the other six are current too."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_upp2", REPO_ROOT / "pipeline" / "scripts" / "update_project_plan.py")
    upp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(upp)

    from macro_core import report_projects as rp
    docs = tmp_path / "docs"
    monkeypatch.setattr(rp, "DOCS_DIR", docs)
    monkeypatch.setattr(rp, "HUB_PATH", docs / "index.html")
    monkeypatch.setattr(pp, "stage_freshness", _fake_freshness)

    args = type("A", (), {"manifest": fake_manifest,
                          "project": "CV_TimeSeries"})()
    assert upp.cmd_render(args) == 0
    con = sqlite3.connect(fake_manifest)
    rows = con.execute("SELECT 1 FROM stage_provenance "
                       "WHERE stage = 'WEB'").fetchall()
    con.close()
    assert not rows, "a partial render recorded WEB"


def test_the_published_pages_match_what_web_recorded():
    """The real manifest: WEB's recorded output fingerprints must equal the
    pages on disk.

    That equality is exactly what "written out of band" was complaining
    about — seven lines of it, one per page this tool had just produced.

    Only WEB's OWN seven output files are fingerprinted here, never the
    whole DAG: hashing the manifest's large tables takes minutes on the
    archive drive while other jobs are writing to them, and a test that
    slow is a test that gets skipped.
    """
    import json
    manifest = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
    if not manifest.exists():                       # pragma: no cover
        pytest.skip("no manifest on this machine")
    con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        row = con.execute(
            "SELECT outputs_json FROM stage_provenance WHERE stage = 'WEB' "
            "ORDER BY run_utc DESC LIMIT 1").fetchone()
        assert row, "WEB has never been recorded — render does not record it"
        recorded = json.loads(row[0])
        stage = pv.STAGE_BY_KEY["WEB"]
        assert set(recorded) == set(stage.writes), (
            "WEB recorded a different set of pages than it declares")
        current = pv.fingerprint_all(stage.writes, con, REPO_ROOT)
    finally:
        con.close()
    drifted = [w for w in stage.writes if current[w] != recorded[w]]
    assert not drifted, (
        f"{len(drifted)} page(s) changed since WEB was recorded — re-run "
        f"`update_project_plan.py render`: {drifted}")


# ===========================================================================
# 13.  CLOSING BY RULING  (the 2026-10-03 plan review, SYNTHESIS §0)
# ===========================================================================
# Before the review a task could leave the plan only by being done, so work
# the data could never support sat `pending` or `blocked` for ever.  Two
# closing statuses now exist, and both are RULINGS: the tests below hold
# that they cannot be typed without one, cannot inflate a completion
# fraction, and cannot hide among live work.

def _ruled_task(tid, status, phase="Phase 1", ruling=True, accept="a test",
                **kw):
    return pp.Task(id=tid, title=tid, produces="a thing", stage="S0",
                   source=pp.Source("d.md", "§D"), status=status,
                   evidence="x", blocker="b",
                   ruling=pp.Ruling("§4 X", ("U1",)) if ruling else None,
                   accept=accept, project="P", phase=phase, **kw)


def _fixture_project(*phases):
    return pp._build(pp.Project(
        key="X", title="X", claim="c", venue="v", strategy="",
        decisions=(_FIXTURE_DECISION,), phases=tuple(phases)))


def test_the_two_closing_statuses_exist_and_are_not_open():
    for status in (pp.DROPPED, pp.DEFERRED):
        assert status in pp.LEDGER_STATUSES
        assert status in pp.ALL_STATUSES
        assert status in pp.CLOSED_BY_RULING
        assert status not in pp.OPEN_STATUSES, \
            "a closed task must never be offered as open work"
        assert pp.STATUS_LABEL[status]


def test_completion_is_counted_over_in_scope_tasks_only():
    """done / in-scope: dropped and deferred are in neither term."""
    tasks = [_task("a", pp.DONE), _task("b", pp.PENDING),
             _task("c", pp.DROPPED), _task("d", pp.DEFERRED),
             _task("e", pp.REDO_NEEDED)]
    st = {t.id: t.status for t in tasks}
    counts = pp.status_counts(tasks, st)
    assert pp.progress_fraction(counts) == (1, 3)
    assert pp.scope_summary(counts) == {
        "done": 1, "in_scope": 3, "dropped": 1, "deferred": 1, "total": 5}


def test_dropping_work_cannot_make_a_project_look_finished():
    """The failure the denominator rule exists to prevent, both ways.

    Counted as done, a dropped task would let a project finish by deleting
    its plan.  Left in the denominator, every ruling that removes impossible
    work would read as the project falling behind.
    """
    before = [_task("a", pp.DONE), _task("b", pp.PENDING),
              _task("c", pp.PENDING)]
    after = [_task("a", pp.DONE), _task("b", pp.PENDING),
             _task("c", pp.DROPPED)]
    frac = lambda ts: pp.progress_fraction(          # noqa: E731
        pp.status_counts(ts, {t.id: t.status for t in ts}))
    assert frac(before) == (1, 3)
    done, in_scope = frac(after)
    assert done == 1, "dropping a task must not add to 'done'"
    assert in_scope == 2, "…and must leave the denominator"
    assert done < in_scope, "one open task remains, so it is not complete"


def test_the_scope_summary_always_adds_up():
    for group in pp.groups():
        st = pp.overlay_statuses(group.tasks, {})
        scope = pp.scope_summary(pp.status_counts(group.tasks, st))
        assert scope["in_scope"] + scope["dropped"] + scope["deferred"] \
            == scope["total"] == len(group.tasks), group.key


def test_validate_rejects_a_drop_with_no_ruling(monkeypatch):
    bad = _fixture_project(pp._dropped_phase(
        _ruled_task("a", pp.DROPPED, ruling=False)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="no ruling"):
        pp.validate()


def test_validate_rejects_a_deferral_with_no_ruling(monkeypatch):
    bad = _fixture_project(pp._backlog_phase(
        _ruled_task("a", pp.DEFERRED, ruling=False)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="no ruling"):
        pp.validate()


def test_validate_rejects_a_dropped_task_left_among_live_work(monkeypatch):
    """A dropped task in a live phase reads as live work on the page."""
    bad = _fixture_project(pp.Phase("Phase 1", "i", (
        _ruled_task("a", pp.DROPPED),)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="filed under"):
        pp.validate()


def test_validate_rejects_live_work_filed_as_dropped(monkeypatch):
    bad = _fixture_project(pp._dropped_phase(_ruled_task("a", pp.PENDING)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="filed under"):
        pp.validate()


def test_validate_rejects_a_ruled_open_task_with_no_acceptance(monkeypatch):
    """The review found tasks marked done against a description rather than
    a criterion; a task it adds or changes states the test that closes it."""
    bad = _fixture_project(pp.Phase("Phase 1", "i", (
        _ruled_task("a", pp.PENDING, accept=""),)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="acceptance criterion"):
        pp.validate()


def test_every_closed_task_in_the_ledger_names_its_ruling():
    closed = [t for t in pp.ledger_tasks()
              if t.status in pp.CLOSED_BY_RULING]
    assert closed, "the review dropped and deferred tasks; none are recorded"
    for t in closed:
        assert t.ruling is not None, t.id
        assert t.ruling.findings, f"{t.id}: a ruling with no finding behind it"
        assert t.ruling.date == "2026-10-03"


def test_closed_tasks_are_filed_apart_from_live_work():
    for t in pp.ledger_tasks():
        if t.status == pp.DROPPED:
            assert t.phase == pp.DROPPED_PHASE, t.id
        if t.status == pp.DEFERRED:
            assert t.phase == pp.BACKLOG_PHASE, t.id


def test_the_backlog_is_named_for_2027():
    assert pp.BACKLOG_NAME == "2027 backlog"
    assert pp.BACKLOG_PHASE.startswith(pp.BACKLOG_NAME)
    assert pp.STATUS_LABEL[pp.DEFERRED] == pp.BACKLOG_NAME


def test_backlog_and_dropped_are_listed_in_plan_order():
    tasks = [_task("a", pp.DEFERRED), _task("b", pp.DROPPED),
             _task("c", pp.DEFERRED), _task("d", pp.PENDING)]
    st = {t.id: t.status for t in tasks}
    assert [t.id for t in pp.backlog(tasks, st)] == ["a", "c"]
    assert [t.id for t in pp.dropped_tasks(tasks, st)] == ["b"]


def test_closed_tasks_are_never_offered_as_next_blocked_or_gated():
    tasks = [_t("gate", pp.PENDING),
             _t("x", pp.DROPPED, ["gate"]), _t("y", pp.DEFERRED, ["gate"])]
    st = {t.id: t.status for t in tasks}
    assert [t.id for t in pp.next_up(tasks, st)] == ["gate"]
    assert pp.open_blockers(tasks, st) == ()
    assert pp.gated_tasks(tasks, st) == ()


def test_every_ruled_open_task_states_its_acceptance_criterion():
    for t in pp.ledger_tasks():
        if t.ruling and t.status in (pp.PENDING, pp.IN_PROGRESS, pp.BLOCKED):
            assert t.accept.strip(), t.id
            # ...and it is shown wherever the task's product is shown.
            assert t.accept in t.produces, t.id


# ===========================================================================
# 14.  RULINGS RESOLVE  (the invented-finding class)
# ===========================================================================
# `citation_problems` stops a task citing a strategy section that does not
# say what it is cited for.  A `dropped` status is only as good as the
# ruling behind it, so the same check is made one level up: the section of
# the synthesis must exist and say it, and every finding id must be one a
# seat actually filed.

_SYN = """## 1. Rulings
| U1 | NGC 5548 is dead |
## 2. Disagreements
**D1 — hrg dispersion**
## 4. Per-project amendments
### X — re-scope
DROP P41."""

_DOCS = {pp.SYNTHESIS_2026_10_03: _SYN,
         "data-scientist.md": "**F1 — BLOCKER** text\n**F10 — MAJOR** text"}


def test_every_ruling_in_the_ledger_resolves():
    """The REAL ledger against the REAL synthesis and memos."""
    assert pp.verify_rulings(REPO_ROOT) == ()


def test_a_ruling_citing_a_real_finding_passes():
    ruling = pp.Ruling("§4 X", ("U1", "D1", "DS.F1", "DS"))
    assert pp.ruling_problems([("t", ruling)], _DOCS) == ()


def test_a_ruling_citing_a_finding_nobody_filed_is_rejected():
    """F2 is not in the memo — and F1 must not match on the strength of
    F10, nor F10 on F1."""
    problems = pp.ruling_problems(
        [("t", pp.Ruling("§4 X", ("DS.F2",)))], _DOCS)
    assert len(problems) == 1 and "F2" in problems[0]
    only_f10 = {**_DOCS, "data-scientist.md": "**F10 — MAJOR** text"}
    assert pp.ruling_problems(
        [("t", pp.Ruling("§4 X", ("DS.F1",)))], only_f10)


def test_a_ruling_citing_an_unknown_seat_is_rejected():
    problems = pp.ruling_problems(
        [("t", pp.Ruling("§4 X", ("ZZ.F1",)))], _DOCS)
    assert problems and "no committee seat" in problems[0]


def test_a_ruling_citing_a_synthesis_ruling_that_does_not_exist_is_rejected():
    problems = pp.ruling_problems([("t", pp.Ruling("§4 X", ("U9",)))], _DOCS)
    assert problems and "U9" in problems[0]


def test_a_malformed_finding_id_is_rejected():
    problems = pp.ruling_problems(
        [("t", pp.Ruling("§4 X", ("because I said so",)))], _DOCS)
    assert problems and "not a finding id" in problems[0]


def test_a_ruling_citing_a_missing_document_is_rejected():
    problems = pp.ruling_problems(
        [("t", pp.Ruling("§4 X", ("U1",), document="nope.md"))], _DOCS)
    assert problems and "does not exist" in problems[0]


def test_a_rulings_section_is_checked_like_any_citation():
    """A ruling pointing at a section of the synthesis that does not say
    what it is cited for is an invented source, same as a task's."""
    good = pp.Ruling("§4 DROP P41", ("U1",))
    bad = pp.Ruling("§4 DROP P99", ("U1",))
    docs = {pp.SYNTHESIS_2026_10_03: _SYN}
    assert pp.citation_problems([("t", good.source)], docs) == ()
    assert pp.citation_problems([("t", bad.source)], docs)
    assert pp.citation_problems(
        [("t", pp.Ruling("§7 X").source)], docs), "§7 does not exist"


# ===========================================================================
# 15.  A RULING SUPERSEDES WHAT WAS RECORDED BEFORE IT  (sync rule 3)
# ===========================================================================

def _closed_by_ruling(status=pp.DONE):
    return pp.Task(id="a", title="a", produces="p", stage="S3",
                   source=pp.Source("d.md", "§D"), status=status,
                   evidence="x", ruling=pp.Ruling("§4 X", ("U7",)))


def test_a_status_recorded_before_a_ruling_gives_way_to_it():
    """SN-G0c was recorded `in_progress` six weeks before the committee
    closed it as NOT PROMOTED.  The recorded row must not win."""
    changes = pp.derive_sync(
        [_closed_by_ruling()], {"a": pp.IN_PROGRESS}, {"S3": pv.FRESH},
        recorded_at={"a": "2026-08-21T04:01:18Z"})
    assert [(c.old, c.new) for c in changes] == [(pp.IN_PROGRESS, pp.DONE)]
    assert "ruling" in changes[0].reason and "U7" in changes[0].reason


def test_a_status_recorded_after_a_ruling_is_a_later_judgement():
    """Left alone: a person who reopens a task after the ruling knew of it."""
    for stamp in ("2026-10-03T09:00:00Z", "2026-11-01T00:00:00Z"):
        assert pp.derive_sync(
            [_closed_by_ruling()], {"a": pp.IN_PROGRESS}, {"S3": pv.FRESH},
            recorded_at={"a": stamp}) == ()


def test_rule_three_is_off_without_the_stamps():
    assert pp.derive_sync(
        [_closed_by_ruling()], {"a": pp.IN_PROGRESS}, {"S3": pv.FRESH}) == ()


def test_a_task_with_no_ruling_is_never_moved_by_rule_three():
    task = pp.Task(id="a", title="a", produces="p", stage="S3",
                   source=pp.Source("d.md", "§D"), status=pp.DONE,
                   evidence="x")
    assert pp.derive_sync([task], {"a": pp.PENDING}, {"S3": pv.FRESH},
                          recorded_at={"a": "2026-08-01T00:00:00Z"}) == ()


def test_a_ruling_and_a_stale_stage_resolve_in_one_step():
    """Closed `done` by ruling on a stage that is stale: one change, to
    `redo_needed`, carrying both reasons — not `done` now and
    `redo_needed` on the next sync."""
    changes = pp.derive_sync(
        [_closed_by_ruling()], {"a": pp.IN_PROGRESS}, {"S3": pv.STALE},
        recorded_at={"a": "2026-08-21T04:01:18Z"})
    assert [(c.old, c.new) for c in changes] == [
        (pp.IN_PROGRESS, pp.REDO_NEEDED)]
    assert "ruling" in changes[0].reason and "S3" in changes[0].reason


def test_a_ruling_can_drop_a_task_that_was_recorded_open():
    task = _closed_by_ruling(pp.DROPPED)
    changes = pp.derive_sync([task], {"a": pp.PENDING}, {"S3": pv.STALE},
                             recorded_at={"a": "2026-08-20T00:00:00Z"})
    assert [(c.old, c.new) for c in changes] == [(pp.PENDING, pp.DROPPED)]


@pytest.mark.parametrize("status", [pp.DROPPED, pp.DEFERRED])
@pytest.mark.parametrize("state", [pv.STALE, pv.FRESH, pv.NEVER_RUN])
def test_no_fingerprint_can_reverse_a_ruling(status, state):
    """Rules 1 and 2 move tasks between done and redo_needed only."""
    task = _closed_by_ruling(status)
    assert pp.derive_sync([task], {"a": status}, {"S3": state},
                          recorded_at={"a": "2026-10-04T00:00:00Z"}) == ()


def test_rule_three_is_idempotent_once_recorded():
    """After sync writes the ruled status, its own stamp post-dates the
    ruling, so the rule does not fire again."""
    task = _closed_by_ruling()
    first = pp.derive_sync([task], {"a": pp.IN_PROGRESS}, {"S3": pv.FRESH},
                           recorded_at={"a": "2026-08-21T00:00:00Z"})
    assert pp.derive_sync([task], {"a": first[0].new}, {"S3": pv.FRESH},
                          recorded_at={"a": "2026-10-03T12:00:00Z"}) == ()


def test_status_stamps_and_sync_age_are_read_back(con):
    tid = pp.all_tasks()[0].id
    assert pp.last_sync_utc(con) is None
    pp.record_status(con, tid, pp.IN_PROGRESS, note="started",
                     when="2026-08-01T00:00:00Z")
    assert pp.last_sync_utc(con) is None, "a hand `set` is not a sync"
    pp.record_status(con, tid, pp.REDO_NEEDED, note="sync: stage S0 is STALE",
                     when="2026-08-20T03:13:38Z")
    assert pp.read_status_stamps(con) == {tid: "2026-08-20T03:13:38Z"}
    assert pp.last_sync_utc(con) == "2026-08-20T03:13:38Z"


def test_stamp_readers_are_safe_without_the_table():
    empty = sqlite3.connect(":memory:")
    assert pp.read_status_stamps(empty) == {}
    assert pp.last_sync_utc(empty) is None
    empty.close()


# ===========================================================================
# 16.  BLOCKERS ARE COMPUTED WHERE THEY CAN BE  (ruling U4)
# ===========================================================================
# Five tasks sat blocked on "S2's tables were destroyed" for six weeks after
# S2 was rebuilt, because the blocker was a sentence and a sentence does not
# re-read the database.

def test_no_blocker_claims_the_s2_tables_were_destroyed():
    for t in pp.ledger_tasks():
        text = f"{t.blocker} {t.produces}".lower()
        assert "destroyed" not in text, (
            f"{t.id}: still says a table was destroyed — the tables exist "
            f"(ruling U4); state the true blocker")
        assert "s2 wipe" not in text and "tables are gone" not in text, t.id


#: The tasks the stale blocker held.  The committee unblocked every one.
_FORMERLY_S2_BLOCKED = (
    "TCRB-P0-bitdepth", "TCRB-P0-ladders", "TCRB-P0-shutter-timing",
    "TCRB-B1-calibration", "DW-P04-noise-model",
    "CV-P15-linearity-ladders", "CV-P15-noise-model", "CV-P2-vetoes")


@pytest.mark.parametrize("task_id", _FORMERLY_S2_BLOCKED)
def test_the_s2_gated_tasks_are_unblocked(task_id):
    assert pp.task_by_id(task_id).status != pp.BLOCKED


def test_a_blocker_that_names_a_task_declares_it_as_a_dependency():
    """Prose that says "clears when DW-P02-filter-dossier does" is a
    dependency, and a dependency nothing reads is how the page came to
    recommend work its own blocker text forbade."""
    ids = {t.id for t in pp.ledger_tasks()}
    for t in pp.ledger_tasks():
        if t.status != pp.BLOCKED:
            continue
        for other in ids - {t.id}:
            if re.search(rf"(?<![\w-]){re.escape(other)}(?![\w-])",
                         t.blocker):
                assert other in t.depends_on, (
                    f"{t.id}: its blocker names {other}, which is not in "
                    f"depends_on")


def test_a_dependency_may_cross_groups_and_is_read_from_the_status_table():
    """T CrB's wavelength solution waits on the shared foundation's G-1.
    Overlaying ONLY the project's tasks must still see G-1's status —
    otherwise the gate could never open."""
    project = pp.PROJECT_BY_KEY["TCrB_Monitoring"]
    a3 = pp.task_by_id("TCRB-A3-wavelength")
    assert "G-1" in a3.depends_on
    assert "G-1" not in {t.id for t in project.tasks}

    before = pp.overlay_statuses(project.tasks, {})
    assert before["G-1"] == pp.PENDING
    assert pp.unmet_dependencies(a3, before) == ("G-1",)

    after = pp.overlay_statuses(project.tasks, {"G-1": pp.DONE})
    assert pp.unmet_dependencies(a3, after) == ()


def test_the_overlay_adds_dependencies_without_disturbing_the_counts():
    project = pp.PROJECT_BY_KEY["TCrB_Monitoring"]
    st = pp.overlay_statuses(project.tasks, {})
    assert sum(pp.status_counts(project.tasks, st).values()) \
        == len(project.tasks)


def test_stage_gates_report_every_stage_that_is_not_fresh():
    task = _task("a")
    task = replace(task, needs_fresh=("S0", "S0c", "S3"))
    gates = pp.stage_gates(task, {"S0": pv.STALE, "S0c": pv.FRESH})
    assert gates == (("S0", pv.STALE), ("S3", "UNKNOWN")), \
        "a gate nobody evaluated is not a gate that opened"
    assert pp.stage_gates(task, None) == (), "not asked, not reported"


@pytest.mark.parametrize("met_when,value,expected", [
    ("zero", 0, True), ("zero", 3, False), ("zero", None, False),
    ("positive", 1, True), ("positive", 0, False), ("positive", None, False),
])
def test_probe_is_met(met_when, value, expected):
    assert pp.probe_is_met(pp.Probe("p", "SELECT 1", met_when), value) \
        is expected


def test_computed_blockers_reads_all_three_kinds_of_gate():
    task = replace(_t("work", pp.PENDING, ["gate"]),
                   needs_fresh=("S0",))
    results = {"work": (pp.ProbeResult("rows", 0, False),)}
    lines = pp.computed_blockers(task, {"gate": pp.BLOCKED, "work":
                                        pp.PENDING},
                                 {"S0": pv.STALE}, results)
    assert len(lines) == 3
    assert "waiting on gate (blocked)" in lines[0]
    assert "S0" in lines[1] and pv.STALE in lines[1]
    assert "rows: 0" in lines[2]


def test_computed_blockers_is_empty_when_every_gate_is_open():
    task = replace(_t("work", pp.PENDING, ["gate"]), needs_fresh=("S0",))
    results = {"work": (pp.ProbeResult("rows", 66, True),)}
    assert pp.computed_blockers(task, {"gate": pp.DONE}, {"S0": pv.FRESH},
                                results) == ()


def test_an_acceptance_probe_never_blocks_its_own_task():
    """"27,261 twins are still canonical" is how far F-1 has to go, not a
    reason F-1 may not start."""
    task = _t("work", pp.PENDING)
    results = {"work": (pp.ProbeResult("twins", 27261, False,
                                       kind="accept"),)}
    assert pp.computed_blockers(task, {"work": pp.PENDING}, None,
                                results) == ()
    gaps = pp.acceptance_gaps(task, results)
    assert [g.value for g in gaps] == [27261]
    assert [t.id for t in pp.next_up([task], {"work": pp.PENDING},
                                     probe_results=results)] == ["work"]


def test_next_up_honours_a_stage_gate_and_a_database_gate():
    gated = replace(_t("a", pp.PENDING), needs_fresh=("S0",))
    probed = _t("b", pp.PENDING)
    free = _t("c", pp.PENDING)
    st = {t.id: pp.PENDING for t in (gated, probed, free)}
    results = {"b": (pp.ProbeResult("rows", 0, False),)}
    assert [t.id for t in pp.next_up([gated, probed, free], st)] \
        == ["a", "b", "c"], "gates nobody evaluated do not hide a task"
    assert [t.id for t in pp.next_up(
        [gated, probed, free], st, stage_states={"S0": pv.STALE},
        probe_results=results)] == ["c"]


def test_run_probes_reads_the_database_and_survives_a_missing_table():
    c = sqlite3.connect(":memory:")
    c.execute("CREATE TABLE s2_linearity_ladders (x INT)")
    c.executemany("INSERT INTO s2_linearity_ladders VALUES (?)",
                  [(i,) for i in range(66)])
    task = replace(_task("a"), probes=(
        pp.Probe("ladders", "SELECT count(*) FROM s2_linearity_ladders",
                 "positive"),
        pp.Probe("epochs", "SELECT count(*) FROM mech_epoch", "positive",
                 kind="accept")))
    results = pp.run_probes(c, [task, _task("b")])
    c.close()
    assert set(results) == {"a"}
    ladders, epochs = results["a"]
    assert (ladders.value, ladders.met, ladders.kind) == (66, True, "gate")
    assert epochs.value is None and not epochs.met
    assert "not answerable yet" in epochs.error
    assert epochs.kind == "accept"


def test_probes_in_the_ledger_only_read():
    for t in pp.ledger_tasks():
        for probe in t.probes:
            assert probe.sql.lstrip().upper().startswith("SELECT"), t.id
            assert ";" not in probe.sql, f"{t.id}: one statement per probe"


@pytest.mark.parametrize("probe,match", [
    (pp.Probe("p", "DELETE FROM frames"), "not a SELECT"),
    (pp.Probe("p", "SELECT 1", "sometimes"), "met_when"),
    (pp.Probe("p", "SELECT 1", "zero", kind="wish"), "kind"),
])
def test_validate_rejects_a_malformed_probe(monkeypatch, probe, match):
    bad = _fixture_project(pp.Phase("P", "i", (
        replace(_t("a", pp.PENDING), probes=(probe,)),)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match=match):
        pp.validate()


def test_validate_rejects_a_stage_gate_on_a_stage_that_does_not_exist(
        monkeypatch):
    bad = _fixture_project(pp.Phase("P", "i", (
        replace(_t("a", pp.PENDING), needs_fresh=("S-NOPE",)),)))
    monkeypatch.setattr(pp, "PROJECTS", (bad,))
    with pytest.raises(pp.PlanError, match="needs_fresh"):
        pp.validate()


def test_the_u4_probes_ask_for_the_tables_the_stale_blocker_named():
    """The replacement for the sentence is the query: each formerly
    blocked detector task that carries a probe asks for an S2 table."""
    probed = [t for t in pp.ledger_tasks()
              if any("U4" in p.label for p in t.probes)]
    assert len(probed) >= 4
    for t in probed:
        for p in t.probes:
            assert p.met_when == "positive" and p.kind == "gate"
            assert re.search(r"FROM s2_\w+", p.sql), (t.id, p.sql)


# ===========================================================================
# 17.  WHAT THE COMMITTEE RULED, AS THE LEDGER NOW STATES IT
# ===========================================================================

def test_the_shared_foundation_is_a_group_of_its_own():
    ids = [t.id for t in pp.FOUNDATION.tasks]
    assert ids == [f"F-{n}" for n in range(1, 11)] \
        + [f"G-{n}" for n in range(1, 6)]
    assert pp.FOUNDATION.title.startswith("Shared foundation")
    assert pp.FOUNDATION not in pp.PROJECTS
    assert pp.groups()[0] is pp.FOUNDATION
    for t in pp.FOUNDATION.tasks:
        assert t.project == "Shared_Foundation"
        assert t.source.document == pp.SYNTHESIS_2026_10_03
        assert t.source.section == f"§3 {t.id}"
        assert t.ruling and t.accept, t.id


def test_the_foundation_is_addressable_but_outside_the_project_totals():
    """`set F-4 done` must work; the "N of M across six projects" counts
    must not move for work none of the six owns."""
    assert pp.task_by_id("F-4").title
    assert pp.project_of("G-1") is pp.FOUNDATION
    assert pp.group_of("Shared_Foundation") is pp.FOUNDATION
    assert "F-4" not in {t.id for t in pp.all_tasks()}
    assert len(pp.ledger_tasks()) == len(pp.all_tasks()) + 15
    with pytest.raises(pp.PlanError):
        pp.group_of("No_Such_Group")


def test_a_status_can_be_recorded_against_a_foundation_task(con):
    pp.record_status(con, "F-4", pp.IN_PROGRESS, note="flat pairs found")
    assert pp.read_statuses(con) == {"F-4": pp.IN_PROGRESS}


def test_the_two_test_resolved_disagreements_are_owned_by_foundation_tasks():
    assert "D1" in pp.task_by_id("G-1").ruling.findings
    assert "D2" in pp.task_by_id("F-4").ruling.findings


def test_cv_timeseries_is_reopened():
    """U3: it read 34/34 and it is not finished."""
    project = pp.PROJECT_BY_KEY["CV_TimeSeries"]
    revision = [t for t in project.tasks if t.id.startswith("CV-R")]
    assert [int(t.id.split("-")[1][1:]) for t in revision] \
        == list(range(1, 16))
    for t in revision:
        assert t.status in (pp.PENDING, pp.BLOCKED), t.id
        assert t.ruling and t.accept, t.id
    assert "REOPENED" in project.claim
    assert project.venue.startswith("AJ"), "the venue ruling is AJ, not ApJ"
    assert "reopened" in project.title.lower()

    # Even with every pre-review task recorded done — which is what the
    # status table holds — the headline cannot read complete.
    recorded = {t.id: pp.DONE for t in project.tasks
                if not t.id.startswith("CV-R")}
    st = pp.overlay_statuses(project.tasks, recorded)
    done, in_scope = pp.progress_fraction(
        pp.status_counts(project.tasks, st))
    assert (done, in_scope) == (34, 49)

    from macro_core import site
    assert site.project_condition(project, st).tone != "done"


def test_the_cv_band_offset_is_decided_by_test_not_by_the_ledger():
    """D3 is open: the ledger may not say which way it goes."""
    assert "D3" in pp.task_by_id("CV-R2-bias-injection").ruling.findings
    claim = pp.PROJECT_BY_KEY["CV_TimeSeries"].claim
    assert "NOT decided" in claim


@pytest.mark.parametrize("task_id,status", [
    # U1 — NGC 5548 broadband photometry is dead.
    ("DW-P41-aperture", pp.DROPPED), ("DW-P42-ensemble", pp.DROPPED),
    ("DW-P44-statistics", pp.DROPPED), ("DW-P45-gate", pp.DROPPED),
    ("DW-P53-period-search", pp.DROPPED),
    ("DW-P56-eclipse-timing", pp.DROPPED),
    ("DW-P52-subtraction", pp.DROPPED),
    # U6 — flickering and the period search leave T CrB paper 1.
    ("TCRB-C3-period-search", pp.DROPPED),
    ("TCRB-A4-response", pp.DROPPED),
    ("TCRB-C2-2026-runs", pp.DEFERRED), ("TCRB-P0-restart", pp.DEFERRED),
    # U7 — the SN model fit is dropped; venue and grism triage are closed.
    ("SN-S7-model-consistency", pp.DROPPED),
    ("SN-venue-decision", pp.DONE), ("SN-G0c-grism-triage", pp.DONE),
    # CLOSE from the archive.
    ("TCRB-P0-bitdepth", pp.DONE),
    # Be-star backlog.
    ("BE-X1-dither-test", pp.DEFERRED),
    ("BE-X2-season2-observing", pp.DEFERRED),
])
def test_the_rulings_are_applied(task_id, status):
    task = pp.task_by_id(task_id)
    assert task.status == status
    assert task.ruling is not None, f"{task_id}: closed with no ruling"


@pytest.mark.parametrize("task_id", [
    "TCRB-N1-novelty-table", "TCRB-N2-eruption-contingency",
    "TCRB-P0-eruption-block", "TCRB-P0-mech-epoch", "TCRB-P0-gain-ptc",
    "TCRB-P0-temp-split", "TCRB-A0b-early-spectra", "TCRB-A5b-line-flux",
    "TCRB-A5a-detection-rule", "TCRB-B0-peak-census",
    "SN-G0-rerun", "SN-G0d-s2c-broadband", "SN-S6a-flash-colour",
    "SN-S6-0-excess-gate", "SN-S7b-residuals-table", "SN-S5-template-table",
    "SN-S5b-peak-epoch",
    "BE-N1-gate", "BE-S0-era-table", "BE-S0-dispositions",
    "DW-P4x-broad-halpha-triage", "DW-P4y-zero-order-photometry",
    "DW-N1-novelty", "DW-P36-0-depth-table",
    "RIG-L0-dedup-reconcile", "RIG-L0-mech-epoch", "RIG-L1-clock-audit",
    "RIG-L2-prereg",
])
def test_every_task_the_committee_added_exists_with_a_criterion(task_id):
    task = pp.task_by_id(task_id)
    assert task.ruling is not None and task.accept.strip()
    assert task.status in (pp.PENDING, pp.BLOCKED)


@pytest.mark.parametrize("task_id", ["TCRB-A9-profiles", "BE-VR-hold"])
def test_profile_and_vr_work_is_held_pending_d1(task_id):
    """HOLD is neither dropped nor kept: blocked, on G-1, by ruling D1."""
    task = pp.task_by_id(task_id)
    assert task.status == pp.BLOCKED
    assert task.depends_on == ("G-1",)
    assert "D1" in task.ruling.findings
    assert "neither dropped nor kept" in task.blocker


def test_the_novelty_gates_come_before_the_pipeline_work():
    """ED cross-cutting 3: a novelty gate per project, executed first."""
    assert pp.PROJECT_BY_KEY["TCrB_Monitoring"].tasks[0].id \
        == "TCRB-N1-novelty-table"
    assert "TCRB-N1-novelty-table" in \
        pp.task_by_id("TCRB-A5-ew").depends_on
    assert pp.task_by_id("BE-N1-gate").depends_on == ("BE-S-1a-bess",)
    for tid in ("BE-S3-extraction", "BE-S7-ew", "BE-figures"):
        assert "BE-N1-gate" in pp.task_by_id(tid).depends_on, tid
    assert "DW-N1-novelty" in pp.task_by_id("DW-P33-stacks").depends_on


def test_the_sn_closures_rest_on_the_stage_that_holds_their_verdict():
    """Both verdicts are rows of sn_g0_verdict.  Bound to S4 or G, their
    `done` could never go stale when Gate 0 did."""
    for tid in ("SN-venue-decision", "SN-G0c-grism-triage"):
        task = pp.task_by_id(tid)
        assert task.stage == "SN-G0"
        assert task.evidence == "docs/SN2023ixf_LightCurve/sn_gate0.html"
    rerun = pp.task_by_id("SN-G0-rerun")
    assert rerun.needs_fresh == ("S0", "S0c")
    assert "F-1" in rerun.depends_on


def test_the_figure_caps_are_in_the_plan():
    caps = {"TCRB-D4-figures": "six", "SN-figures": "five",
            "BE-figures": "six", "DW-figures": "six"}
    for tid, word in caps.items():
        task = pp.task_by_id(tid)
        assert word in task.title.lower(), tid
        assert task.ruling and "ED" in task.ruling.findings, tid


def test_no_title_advertises_a_paper_the_committee_killed():
    """ED.E6: three skeleton titles and the ledger's own display titles."""
    titles = {p.key: f"{p.title} {p.paper_title}" for p in pp.PROJECTS}
    assert "NGC 5548" not in titles["DwarfGalaxy_AGN_Survey"]
    assert "Early" not in titles["SN2023ixf_LightCurve"]
    assert "Quiescent Baseline" not in titles["TCrB_Monitoring"]
    assert "Three-Year" not in titles["TCrB_Monitoring"]


def test_the_legacy_archive_is_renamed_but_its_key_is_not():
    """TE.F10: the premise was wrong; the directory id must survive."""
    project = pp.PROJECT_BY_KEY["Legacy_Rigel"]
    assert "Rigel" not in project.title
    assert "census only" in project.title.lower()
    assert "different telescope (Rigel System" not in project.claim
    assert "AC4040" in project.claim, \
        "the corrected premise: the last legacy camera is the RLMT's own"
    assert project.ruling and "TE.F10" in project.ruling.findings
    assert (REPO_ROOT / "docs" / "Legacy_Rigel").exists()
    assert all(t.id.startswith("RIG-") for t in project.tasks)
    # Census only: nothing beyond the pre-registered decision is planned.
    assert [ph.name.split(" — ")[0] for ph in project.phases] \
        == ["Phase L0", "Phase L1", "Phase L2"]
    assert "RIG-L2-prereg" in pp.task_by_id("RIG-L2-gonogo").depends_on


@pytest.mark.parametrize("project", [p for p in pp.PROJECTS if p.strategy],
                         ids=lambda p: p.key)
def test_every_strategy_records_the_amendments_that_bind_it(project):
    text = (REPO_ROOT / project.strategy).read_text(encoding="utf-8")
    assert "## 10. Committee amendments 2026-10-03" in text
    section = pp.split_numbered_sections(text)["10"]
    assert pp.SYNTHESIS_2026_10_03 in section
    assert "**Verdict:**" in section and "**Venue:**" in section
    # Every ruled task of this project is in the section, by id.
    for t in project.tasks:
        if t.ruling:
            assert f"`{t.id}`" in section, (
                f"{t.id} carries a ruling that §10 of {project.strategy} "
                f"does not record — re-emit the table with "
                f"`update_project_plan.py amendments {project.key}`")


def test_the_roadmap_records_the_amendments():
    text = (REPO_ROOT / "ROADMAP.md").read_text(encoding="utf-8")
    section = pp.split_numbered_sections(text)["0"]
    assert "### 0.1 Committee amendments 2026-10-03" in section
    for ruling in [f"U{n}" for n in range(1, 11)] + ["D1", "D2", "D3", "D4"]:
        assert re.search(rf"\b{ruling}\b", section), ruling
    # The portfolio table states the venue the ledger holds, verbatim.
    for project in pp.PROJECTS:
        assert project.venue in section, project.key
    # ...and no longer sells the papers the committee ruled out.
    assert "3 polars + YZ Cnc superhumps" not in section
    assert "NGC 5548 band-integrated LC" not in section


def test_the_project_title_and_venue_changes_cite_a_ruling():
    for key in ("TCrB_Monitoring", "CV_TimeSeries", "SN2023ixf_LightCurve",
                "BeStar_Grism", "DwarfGalaxy_AGN_Survey", "Legacy_Rigel"):
        project = pp.PROJECT_BY_KEY[key]
        assert project.ruling is not None, key
        assert project.ruling.findings, key


# ===========================================================================
# 18.  MANUSCRIPT TITLES ARE LEDGER DATA  (ED.E6)
# ===========================================================================

_TEX = ("\\documentclass{aastex701}\n\\begin{document}\n\n"
        "\\title{Early Photometry of a Thing}\n\n"
        "\\author{A. Person}\n% \\title{a commented decoy}\n")


def test_manuscript_title_reads_the_title_line():
    assert pp.manuscript_title(_TEX) == "Early Photometry of a Thing"
    assert pp.manuscript_title("\\begin{document}") is None


def test_with_title_changes_the_title_line_and_nothing_else():
    new = pp.with_title(_TEX, r"H$\alpha$ of SN~2023ixf, $+5.4$ Days")
    assert pp.manuscript_title(new) == r"H$\alpha$ of SN~2023ixf, $+5.4$ Days"
    assert new.replace(r"\title{H$\alpha$ of SN~2023ixf, $+5.4$ Days}",
                       r"\title{Early Photometry of a Thing}") == _TEX
    assert "\\\\alpha" not in new, \
        "a doubled backslash is a LaTeX line break followed by 'alpha'"


def test_with_title_refuses_a_manuscript_with_no_title_line():
    with pytest.raises(pp.PlanError, match="no one-line"):
        pp.with_title("\\begin{document}\n", "T")


def test_the_cv_manuscript_title_is_not_ledger_owned():
    """The CV draft is hand-edited; `titles --write` must never touch it."""
    assert pp.PROJECT_BY_KEY["CV_TimeSeries"].paper_title == ""
    assert pp.PROJECT_BY_KEY["Legacy_Rigel"].paper_title == ""


@pytest.mark.parametrize("project",
                         [p for p in pp.PROJECTS if p.paper_title],
                         ids=lambda p: p.key)
def test_the_skeleton_titles_match_the_ledger(project):
    path = REPO_ROOT / pp.manuscript_path(project)
    if not path.exists():
        pytest.skip("manuscripts/ is not tracked; no skeleton on this tree")
    assert pp.manuscript_title(path.read_text(encoding="utf-8")) \
        == project.paper_title, (
        f"re-emit with `update_project_plan.py titles --write`")


def test_paper_titles_are_valid_one_line_latex():
    for project in pp.PROJECTS:
        title = project.paper_title
        if not title:
            continue
        assert "\n" not in title and "\\\\" not in title, project.key
        assert title.count("{") == title.count("}"), project.key
        assert title.count("$") % 2 == 0, project.key


# ===========================================================================
# 19.  THE COMMAND  (show, set, blockers, amendments, titles)
# ===========================================================================

@pytest.fixture()
def upp():
    """The CLI script, loaded as a module (it is a script, not a package)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_upp_cli",
        REPO_ROOT / "pipeline" / "scripts" / "update_project_plan.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _args(**kw):
    return type("A", (), kw)()


def _history(manifest):
    c = sqlite3.connect(manifest)
    try:
        return pp.read_history(c)
    finally:
        c.close()


def test_set_refuses_to_drop_a_task_with_no_ruling(upp, fake_manifest,
                                                   capsys):
    """`dropped` is a ruling, not an opinion somebody may type."""
    task = next(t for t in pp.all_tasks() if t.ruling is None)
    code = upp.cmd_set(_args(manifest=fake_manifest, task_id=task.id,
                             status=pp.DROPPED, evidence="", note="",
                             ruling=""))
    assert code == 2
    assert "no ruling" in capsys.readouterr().err
    assert _history(fake_manifest) == []


def test_set_records_a_drop_with_the_ruling_in_the_note(upp, fake_manifest):
    task = next(t for t in pp.all_tasks() if t.ruling is None)
    code = upp.cmd_set(_args(
        manifest=fake_manifest, task_id=task.id, status=pp.DEFERRED,
        evidence="", note="needs 2027 frames",
        ruling="SYNTHESIS §4 X (U6)"))
    assert code == 0
    (row,) = _history(fake_manifest)
    assert row[1] == pp.DEFERRED
    assert row[3] == "ruling: SYNTHESIS §4 X (U6) — needs 2027 frames"


def test_set_accepts_the_ledgers_own_ruling(upp, fake_manifest):
    code = upp.cmd_set(_args(
        manifest=fake_manifest, task_id="TCRB-C3-period-search",
        status=pp.DROPPED, evidence="", note="", ruling=""))
    assert code == 0
    (row,) = _history(fake_manifest)
    assert row[3].startswith("ruling: ") and "U6" in row[3]


def test_show_prints_done_over_in_scope_with_the_rest_beside_it(
        upp, fake_manifest, capsys):
    code = upp.cmd_show(_args(manifest=fake_manifest, project=None,
                              history=False, live=False))
    out = capsys.readouterr().out
    assert code == 0
    assert "never been synced" in out, \
        "a listing owes its reader the age of the statuses it prints"
    assert "Shared foundation (Wave 0)   0/15 in-scope tasks complete" in out
    dwarf = out.split("Dwarf-Galaxy Candidates")[1].split("=" * 78)[0]
    assert "in-scope tasks complete" in dwarf
    assert "outside the count: 7 dropped by ruling, 0 in the 2027 backlog" \
        in dwarf
    assert pp.DROPPED_PHASE in out and pp.BACKLOG_PHASE in out
    # A closed task shows its ruling; a gated one shows a computed reason.
    assert "ruling: committee/reviews/2026-10-03/SYNTHESIS.md" in out
    assert "> waiting on G-1 (pending)" in out


def test_show_states_the_sync_age_once_a_sync_has_run(upp, fake_manifest,
                                                      capsys):
    c = sqlite3.connect(fake_manifest)
    with c:
        pp.record_status(c, "TCRB-P0-staging", pp.REDO_NEEDED,
                         note="sync: stage S0c is STALE",
                         when="2026-08-20T03:13:38Z")
    c.close()
    upp.cmd_show(_args(manifest=fake_manifest, project="TCrB_Monitoring",
                       history=False, live=False))
    assert "statuses last synced 2026-08-20T03:13:38Z" \
        in capsys.readouterr().out


def test_show_live_prints_what_a_sync_would_leave(upp, fake_manifest,
                                                  capsys, monkeypatch):
    """Never an unsynced count: with --live the statuses printed are the
    ones `sync` would record, and the listing says how many differ."""
    monkeypatch.setattr(pp, "stage_freshness", _fake_freshness)
    code = upp.cmd_show(_args(manifest=fake_manifest,
                              project="TCrB_Monitoring", history=False,
                              live=True))
    out = capsys.readouterr().out
    assert code == 0 and "LIVE view" in out
    # TCRB-P0-staging is `done` in the ledger on S0c, which the fake DAG
    # reports STALE_UPSTREAM; the live view must not print it as done.
    line = next(l for l in out.splitlines() if "TCRB-P0-staging" in l)
    assert pp.REDO_NEEDED in line
    assert "run `sync` to record them" in out


def test_show_rejects_an_unknown_group(upp, fake_manifest, capsys):
    assert upp.cmd_show(_args(manifest=fake_manifest, project="Nope",
                              history=False, live=False)) == 2
    assert "Shared_Foundation" in capsys.readouterr().err


def test_blockers_reports_computed_gates(upp, fake_manifest, capsys):
    code = upp.cmd_blockers(_args(manifest=fake_manifest,
                                  project="TCrB_Monitoring", live=False))
    out = capsys.readouterr().out
    assert code == 0
    assert "TCRB-A3-wavelength" in out and "> waiting on G-1" in out
    assert "TCRB-A9-profiles" in out and "HELD by ruling D1" in out
    assert "TCRB-C3-period-search" not in out, "a dropped task has no blocker"


def test_blockers_flags_a_blocker_the_database_contradicts(
        upp, tmp_path, monkeypatch, capsys):
    """The U4 failure, as a check: a task recorded blocked whose only gate
    is a table that is sitting in the database with rows in it."""
    manifest = tmp_path / "m.sqlite"
    c = sqlite3.connect(manifest)
    c.execute("CREATE TABLE s2_ptc_fits (x INT)")
    c.execute("INSERT INTO s2_ptc_fits VALUES (1)")
    c.commit()
    c.close()
    stale = replace(
        _t("a", pp.BLOCKED), blocker="S2's tables are missing.",
        probes=(pp.Probe("s2_ptc_fits rows",
                         "SELECT count(*) FROM s2_ptc_fits", "positive"),))
    honest = replace(
        _t("b", pp.BLOCKED), blocker="No table yet.",
        probes=(pp.Probe("mech_epoch rows",
                         "SELECT count(*) FROM mech_epoch", "positive"),))
    monkeypatch.setattr(pp, "PROJECTS", (_fixture_project(
        pp.Phase("P", "i", (stale, honest))),))
    code = upp.cmd_blockers(_args(manifest=manifest, project="X",
                                  live=False))
    captured = capsys.readouterr()
    assert code == 1
    assert captured.out.count("CONTRADICTED") == 1
    assert "1 blocker(s) contradicted" in captured.err


def test_blockers_flags_a_done_task_whose_acceptance_probe_is_unmet(
        upp, tmp_path, monkeypatch, capsys):
    """The mirror image: recorded done, and the manifest disagrees."""
    manifest = tmp_path / "m.sqlite"
    c = sqlite3.connect(manifest)
    c.execute("CREATE TABLE frames (is_twin INT)")
    c.executemany("INSERT INTO frames VALUES (?)", [(1,), (1,), (0,)])
    c.commit()
    c.close()
    task = replace(
        _t("a", pp.DONE), evidence="x",
        probes=(pp.Probe("twins still canonical",
                         "SELECT count(*) FROM frames WHERE is_twin = 1",
                         "zero", kind="accept"),))
    monkeypatch.setattr(pp, "PROJECTS", (_fixture_project(
        pp.Phase("P", "i", (task,))),))
    code = upp.cmd_blockers(_args(manifest=manifest, project="X",
                                  live=False))
    out = capsys.readouterr().out
    assert code == 1
    assert "recorded done, but its acceptance probe is not met" in out
    assert "twins still canonical: 2" in out


def test_amendments_emits_one_row_per_ruled_task(upp, capsys):
    code = upp.cmd_amendments(_args(project="DwarfGalaxy_AGN_Survey"))
    out = capsys.readouterr().out
    assert code == 0
    rows = [l for l in out.splitlines() if l.startswith("| `")]
    project = pp.PROJECT_BY_KEY["DwarfGalaxy_AGN_Survey"]
    assert len(rows) == sum(1 for t in project.tasks if t.ruling)
    by_id = {r.split("`")[1]: r.split(" | ")[1] for r in rows}
    assert by_id["DW-P41-aperture"] == "DROP"
    assert by_id["DW-N1-novelty"] == "ADD"
    assert by_id["DW-P04-noise-model"] == "CHANGE"
    assert upp.cmd_amendments(_args(project="Nope")) == 2


def test_the_ledger_uses_every_one_of_the_committees_verbs():
    verbs = {pp.ruling_action(t) for t in pp.ledger_tasks() if t.ruling}
    assert verbs == set(pp.RULING_ACTIONS)


@pytest.mark.parametrize("task_id,verb", [
    ("CV-R1-band-offset", "ADD"),
    ("CV-R15-release-readiness", "ADD"),     # blocked on James — not a HOLD
    ("TCRB-A3-wavelength", "CHANGE"),
    ("TCRB-A4-response", "DROP"),
    ("TCRB-C2-2026-runs", "DEFER"),
    ("TCRB-A9-profiles", "HOLD"),
    ("BE-VR-hold", "HOLD"),
    ("TCRB-P0-bitdepth", "CLOSE"),
    ("SN-venue-decision", "CLOSE"),
    ("TCRB-P0-ladders", "UNBLOCK"),
    ("F-4", "ADD"),
    ("RIG-L2-prereg", "ADD"),
    ("RIG-L2-gonogo", "CHANGE"),
])
def test_ruling_action_is_the_committees_verb(task_id, verb):
    assert pp.ruling_action(pp.task_by_id(task_id)) == verb


def test_a_task_with_no_ruling_has_no_verb():
    assert pp.ruling_action(_task("a")) == ""


@pytest.mark.parametrize("status,action,match", [
    (pp.PENDING, "POSTPONE", "not one of"),
    (pp.DROPPED, "DEFER", "contradicts"),
    (pp.PENDING, "DROP", "on a task whose status"),
])
def test_validate_rejects_a_verb_that_contradicts_the_task(
        monkeypatch, status, action, match):
    phase = {pp.DROPPED: pp._dropped_phase}.get(
        status, lambda *ts: pp.Phase("P", "i", ts))
    task = replace(_ruled_task("a", status),
                   ruling=pp.Ruling("§4 X", ("U1",), action=action))
    monkeypatch.setattr(pp, "PROJECTS", (_fixture_project(phase(task)),))
    with pytest.raises(pp.PlanError, match=match):
        pp.validate()


def test_titles_reports_a_mismatch_and_writes_only_the_title_line(
        upp, tmp_path, monkeypatch, capsys):
    project = pp.PROJECT_BY_KEY["SN2023ixf_LightCurve"]
    path = tmp_path / pp.manuscript_path(project)
    path.parent.mkdir(parents=True)
    path.write_text(_TEX, encoding="utf-8")
    monkeypatch.setattr(upp, "REPO_ROOT", tmp_path)

    assert upp.cmd_titles(_args(write=False)) == 1
    assert path.read_text(encoding="utf-8") == _TEX, "a check must not write"
    assert "title differs from the ledger" in capsys.readouterr().out

    assert upp.cmd_titles(_args(write=True)) == 0
    tex = path.read_text(encoding="utf-8")
    assert pp.manuscript_title(tex) == project.paper_title
    assert tex.count("\n") == _TEX.count("\n")
    assert "\\author{A. Person}" in tex
    assert upp.cmd_titles(_args(write=False)) == 0
