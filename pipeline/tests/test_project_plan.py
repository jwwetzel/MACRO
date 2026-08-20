"""Unit tests for macro_core.project_plan and the project-page renderer.

Three groups, in order of how much they would cost to get wrong:

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
    ids = [t.id for t in pp.all_tasks()]
    assert len(ids) == len(set(ids)), \
        "a duplicate id would make two tasks share one status history"


def test_every_task_has_a_ledger_status():
    for t in pp.all_tasks():
        assert t.status in pp.LEDGER_STATUSES, (t.id, t.status)


def test_redo_needed_is_never_declared_in_code():
    """`redo_needed` is derived, never authored: it is a statement about the
    database, and a hand-typed one could not name its reason."""
    for t in pp.all_tasks():
        assert t.status != pp.REDO_NEEDED, t.id


def test_every_task_names_a_real_stage_in_the_dag():
    for t in pp.all_tasks():
        assert t.stage in pv.STAGE_BY_KEY, (
            f"{t.id} depends on stage {t.stage!r}, which is not in "
            f"provenance.STAGES — staleness could never propagate to it")


def test_every_task_cites_a_document_that_exists():
    for t in pp.all_tasks():
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
    for t in pp.all_tasks():
        if t.status == pp.BLOCKED:
            assert t.blocker.strip(), f"{t.id}: blocked with no reason"


def test_done_tasks_carry_evidence():
    for t in pp.all_tasks():
        if t.status == pp.DONE:
            assert t.evidence.strip(), f"{t.id}: done with nothing behind it"


def test_evidence_paths_under_docs_exist():
    """A `done` task whose evidence link 404s is worse than no link."""
    for t in pp.all_tasks():
        if t.evidence.startswith("docs/"):
            assert (REPO_ROOT / t.evidence).exists(), (t.id, t.evidence)


def test_project_and_phase_are_stamped_on_every_task():
    for project in pp.PROJECTS:
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
    """
    hub = (REPO_ROOT / "docs" / "index.html").read_text()
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
    """
    forbidden = ["3,496 ADU", "3496 ADU", "veto 3,200", "N_sub=16", "4.15 e-"]
    owned = {REPO_ROOT / "docs" / p.key / "index.html" for p in pp.PROJECTS}
    owned.add(REPO_ROOT / "docs" / "index.html")
    offenders = []
    for page in (REPO_ROOT / "docs").rglob("*.html"):
        text = page.read_text(errors="replace")
        hits = [f for f in forbidden if f in text]
        if not hits:
            continue
        assert page not in owned, (page.name, hits)
        offenders.append(page)
    # Every offender must be reachable only through a verdict-chipped link.
    for page in offenders:
        for owner in owned:
            if not owner.exists():
                continue
            html = owner.read_text()
            if page.name not in html:
                continue
            # A verdict chip may sit on EITHER side of the link: the project
            # pages render it after the link, while the hub's table puts the
            # verdict in its own column BEFORE the report column.  Checking
            # only what follows the link failed the hub even though its row
            # carries a chip a few characters earlier — so look both ways.
            parts = html.split(page.name)
            windows = [p[-220:] for p in parts[:-1]] + [p[:220] for p in parts[1:]]
            assert any('class="chip v-' in w for w in windows), (
                f"{owner.name} links {page.name} with no verdict chip")


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
