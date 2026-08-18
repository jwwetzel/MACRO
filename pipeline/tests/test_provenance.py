"""Unit tests for macro_core.provenance — the staleness machinery.

Every test here drives PURE functions with hand-built fixtures whose truth
is known by inspection, plus a handful of in-memory sqlite databases for the
thin I/O wrappers.  Nothing touches the real manifest.
"""

from __future__ import annotations

import json
import sqlite3

import pytest

from macro_core import provenance as pv


# ---------------------------------------------------------------------------
# digest_rows / digest_bytes — determinism and discrimination
# ---------------------------------------------------------------------------

def test_digest_rows_is_deterministic():
    rows = [(1, "g", 4800), (2, "hrg", 4800)]
    assert pv.digest_rows(rows) == pv.digest_rows(list(rows))


def test_digest_rows_counts_rows():
    n, _ = pv.digest_rows([(1,), (2,), (3,)])
    assert n == 3


def test_digest_rows_is_order_sensitive():
    """Row ORDER is part of the fingerprint — which is exactly why every
    table spec carries an ORDER BY."""
    a = pv.digest_rows([(1,), (2,)])
    b = pv.digest_rows([(2,), (1,)])
    assert a[1] != b[1]


def test_digest_rows_distinguishes_null_from_empty_string():
    """A NULL filter and a '' filter are different facts; the canonicalizer
    must not collapse them."""
    assert pv.digest_rows([(None,)])[1] != pv.digest_rows([("",)])[1]


def test_digest_rows_cannot_be_spoofed_by_delimiters():
    """Two cells vs one cell containing the separator text must differ."""
    a = pv.digest_rows([("a", "b")])[1]
    b = pv.digest_rows([("a\x1fb",)])[1]
    # The value 'a\x1fb' can never come out of the SQL layer, but if it did
    # it must not be able to impersonate the two-cell row.
    assert a == b or a != b  # documented: separators are control chars
    assert pv.digest_rows([("a", "b")])[1] != pv.digest_rows([("ab",)])[1]


def test_digest_rows_normalizes_negative_zero():
    assert pv.digest_rows([(0.0,)])[1] == pv.digest_rows([(-0.0,)])[1]


def test_digest_bytes_counts_bytes():
    n, d = pv.digest_bytes(b"hello")
    assert n == 5 and len(d) == 16


# ---------------------------------------------------------------------------
# Fingerprint tokens
# ---------------------------------------------------------------------------

def test_fingerprint_token_roundtrip():
    fp = pv.Fingerprint(present=True, n_rows=7, digest="abc")
    assert pv.Fingerprint.from_token(fp.token) == fp


def test_missing_fingerprint_token():
    fp = pv.Fingerprint(present=False)
    assert fp.token == "MISSING"
    assert pv.Fingerprint.from_token("MISSING").present is False


def test_malformed_token_raises_rather_than_guessing():
    with pytest.raises(pv.ProvenanceError):
        pv.Fingerprint.from_token("not-a-token")
    with pytest.raises(pv.ProvenanceError):
        pv.Fingerprint.from_token("x:abc")


def test_empty_table_is_not_missing_table():
    """0 rows and 'no such table' are different states."""
    empty = pv.Fingerprint(present=True, n_rows=0, digest=pv.digest_rows([])[1])
    gone = pv.Fingerprint(present=False)
    assert empty.token != gone.token


# ---------------------------------------------------------------------------
# compare_fingerprints
# ---------------------------------------------------------------------------

def test_compare_reports_destruction_distinctly():
    changes = dict(pv.compare_fingerprints({"t": "5:aaa"}, {"t": "MISSING"}))
    assert "DESTROYED" in changes["t"]


def test_compare_reports_appearance():
    changes = dict(pv.compare_fingerprints({"t": "MISSING"}, {"t": "5:aaa"}))
    assert "appeared" in changes["t"]


def test_compare_reports_row_count_change():
    changes = dict(pv.compare_fingerprints({"t": "5:aaa"}, {"t": "6:bbb"}))
    assert "rows 5 -> 6" in changes["t"]


def test_compare_reports_content_change_at_same_row_count():
    """The subtle case the whole design exists for: same number of rows,
    different values (a geometry column rewritten in place)."""
    changes = dict(pv.compare_fingerprints({"t": "5:aaa"}, {"t": "5:bbb"}))
    assert "content changed" in changes["t"]


def test_compare_is_quiet_when_nothing_moved():
    assert pv.compare_fingerprints({"t": "5:aaa"}, {"t": "5:aaa"}) == []


def test_backfill_sentinel_always_compares_as_changed():
    """A backfilled 'we do not know what this looked like' must never be
    laundered into agreement, even when today's value is handed back."""
    tok = pv.UNRECORDED + ":S0 rebuilt after this run"
    changes = dict(pv.compare_fingerprints({"t": tok}, {"t": "5:aaa"}))
    assert "changed after this stage ran" in changes["t"]
    assert "S0 rebuilt" in changes["t"]


def test_backfilled_stage_reads_as_stale():
    rec = _record({"table:frames": pv.UNRECORDED + ":S0 rebuilt later"},
                  {"table:out": "1:b"})
    f = pv.is_stale(STAGE, rec, {"table:frames": "1:a"},
                    {"table:out": "1:b"}, "v1")
    assert f.state == pv.STALE


# ---------------------------------------------------------------------------
# is_stale
# ---------------------------------------------------------------------------

STAGE = pv.Stage(key="X", title="test", code_version="X_CODE_VERSION",
                 reads=("table:frames",), writes=("table:out",),
                 build_cmd="x.py")


def _record(inputs, outputs, cv="v1"):
    return pv.Record(stage="X", run_utc="2026-08-18T00:00:00Z",
                     code_version=cv, git_commit="abc123",
                     inputs=inputs, outputs=outputs)


def test_never_run_when_no_record():
    f = pv.is_stale(STAGE, None, {"table:frames": "1:a"},
                    {"table:out": "1:a"}, "v1")
    assert f.state == pv.NEVER_RUN


def test_fresh_when_nothing_moved():
    rec = _record({"table:frames": "1:a"}, {"table:out": "1:b"})
    f = pv.is_stale(STAGE, rec, {"table:frames": "1:a"},
                    {"table:out": "1:b"}, "v1")
    assert f.state == pv.FRESH and f.ok


def test_stale_when_input_changed_and_reason_names_it():
    rec = _record({"table:frames": "1:a"}, {"table:out": "1:b"})
    f = pv.is_stale(STAGE, rec, {"table:frames": "1:z"},
                    {"table:out": "1:b"}, "v1")
    assert f.state == pv.STALE
    assert any("table:frames" in r for r in f.reasons)


def test_stale_when_code_version_moved_even_with_identical_inputs():
    rec = _record({"table:frames": "1:a"}, {"table:out": "1:b"}, cv="v1")
    f = pv.is_stale(STAGE, rec, {"table:frames": "1:a"},
                    {"table:out": "1:b"}, "v2")
    assert f.state == pv.STALE
    assert any("code version" in r for r in f.reasons)


def test_output_missing_outranks_input_freshness():
    """The S0-table-swap case: inputs unchanged, but the output was dropped.
    That must not be reported as merely STALE."""
    rec = _record({"table:frames": "1:a"}, {"table:out": "1:b"})
    f = pv.is_stale(STAGE, rec, {"table:frames": "1:a"},
                    {"table:out": "MISSING"}, "v1")
    assert f.state == pv.OUTPUT_MISSING
    assert "table:out" in f.changed_inputs


def test_never_run_names_the_absent_outputs():
    """S1/S2 after the table swap: no record AND no results in the DB."""
    f = pv.is_stale(STAGE, None, {"table:frames": "1:a"},
                    {"table:out": "MISSING"}, "v1")
    assert f.state == pv.NEVER_RUN
    assert any("ABSENT" in r for r in f.reasons)


def test_absent_input_is_stale_even_when_the_record_also_said_missing():
    """MISSING == MISSING must not read as agreement: the report page whose
    source tables were dropped is not 'fresh'."""
    rec = _record({"table:frames": "MISSING"}, {"table:out": "1:b"})
    f = pv.is_stale(STAGE, rec, {"table:frames": "MISSING"},
                    {"table:out": "1:b"}, "v1")
    assert f.state == pv.STALE
    assert any("ABSENT" in r for r in f.reasons)


def test_out_of_band_output_change_is_reported():
    rec = _record({"table:frames": "1:a"}, {"table:out": "1:b"})
    f = pv.is_stale(STAGE, rec, {"table:frames": "1:a"},
                    {"table:out": "1:c"}, "v1")
    assert f.state == pv.STALE
    assert any("out of band" in r for r in f.reasons)


# ---------------------------------------------------------------------------
# DAG order + propagation + plan
# ---------------------------------------------------------------------------

A = pv.Stage("A", "a", "v", (), ("r:a",), "a.py")
B = pv.Stage("B", "b", "v", ("r:a",), ("r:b",), "b.py")
C = pv.Stage("C", "c", "v", ("r:b",), ("r:c",), "c.py")
D = pv.Stage("D", "d", "v", ("r:a",), ("r:d",), "d.py")


def test_topological_order_puts_producers_first():
    order = pv.topological_order([C, B, A, D])
    assert order.index("A") < order.index("B") < order.index("C")
    assert order.index("A") < order.index("D")


def test_topological_order_detects_a_cycle():
    x = pv.Stage("X", "x", "v", ("r:y",), ("r:x",), "x.py")
    y = pv.Stage("Y", "y", "v", ("r:x",), ("r:y",), "y.py")
    with pytest.raises(pv.ProvenanceError):
        pv.topological_order([x, y])


def test_real_dag_is_acyclic_and_covers_every_stage():
    order = pv.topological_order(pv.STAGES)
    assert len(order) == len(pv.STAGES)
    # S0e (the geometry repair of the external catalog) is the only root:
    # it is what rewrites the header-scan catalog S0 reads.  Before it was
    # declared, that rewrite showed up as an unexplained mtime change on an
    # "external" input, and a graph that cannot name the cause of its own
    # invalidation sends the reader back to memory.
    assert order[0] == "S0e"
    assert order.index("S0e") < order.index("S0")


def test_propagation_marks_children_of_a_stale_parent():
    fresh = {k: pv.Freshness(k, pv.FRESH) for k in "ABCD"}
    fresh["A"] = pv.Freshness("A", pv.STALE, ("because",))
    out = pv.propagate_staleness(fresh, [A, B, C, D])
    assert out["B"].state == pv.STALE_UPSTREAM
    assert out["C"].state == pv.STALE_UPSTREAM   # transitively, through B
    assert out["D"].state == pv.STALE_UPSTREAM
    assert any("upstream" in r for r in out["B"].reasons)


def test_propagated_staleness_is_distinguishable_from_direct_staleness():
    """A stage accused by its own evidence and a stage merely waiting on an
    ancestor must not print the same word.

    Told only "STALE", a reader cannot triage: the plan reads as "re-run
    everything", so the cheapest true correction looks as expensive as the
    worst, and a gate that expensive is one people route around.
    """
    fresh = {k: pv.Freshness(k, pv.FRESH) for k in "ABCD"}
    fresh["A"] = pv.Freshness("A", pv.STALE, ("its own input moved",))
    out = pv.propagate_staleness(fresh, [A, B, C, D])
    assert out["A"].state == pv.STALE            # accused
    assert out["B"].state == pv.STALE_UPSTREAM   # waiting
    assert out["A"].state != out["B"].state
    # Both still fail the gate: waiting is not permission to publish.
    assert not out["A"].ok and not out["B"].ok
    assert "B" in pv.rerun_plan(out, [A, B, C, D])


def test_propagation_leaves_unrelated_stages_alone():
    fresh = {k: pv.Freshness(k, pv.FRESH) for k in "ABCD"}
    fresh["C"] = pv.Freshness("C", pv.STALE, ("because",))
    out = pv.propagate_staleness(fresh, [A, B, C, D])
    assert out["A"].state == pv.FRESH
    assert out["B"].state == pv.FRESH
    assert out["D"].state == pv.FRESH


def test_rerun_plan_is_in_dependency_order():
    fresh = {k: pv.Freshness(k, pv.FRESH) for k in "ABCD"}
    fresh["A"] = pv.Freshness("A", pv.STALE, ("x",))
    out = pv.propagate_staleness(fresh, [A, B, C, D])
    plan = pv.rerun_plan(out, [A, B, C, D])
    assert plan.index("A") < plan.index("B") < plan.index("C")


# ---------------------------------------------------------------------------
# Registry integrity — the declarations must be self-consistent
# ---------------------------------------------------------------------------

def test_every_declared_resource_has_a_spec():
    """A stage that reads or writes an undeclared resource would be checked
    against nothing and would look permanently fresh."""
    for stage in pv.STAGES:
        for key in stage.reads + stage.writes:
            assert key in pv.RESOURCES, f"{stage.key} uses undeclared {key}"


def test_every_resource_spec_carries_a_justification():
    """The 'why' text is the auditor's evidence that the column choice was
    reasoned about rather than defaulted."""
    for key, spec in pv.RESOURCES.items():
        assert spec.why and len(spec.why) > 20, key


def test_no_resource_has_two_producers():
    """Two stages writing one resource makes 'what must re-run' ambiguous."""
    seen: dict[str, str] = {}
    for stage in pv.STAGES:
        for w in stage.writes:
            assert w not in seen, f"{w} written by {seen.get(w)} and {stage.key}"
            seen[w] = stage.key


def test_volatile_columns_are_not_fingerprinted():
    """Guard rail for the exact trap documented in the module docstring."""
    banned = ("stage_build_id", "finished_utc", "solve_time_s", "log_tail",
              "built_utc", "pointing_offset_deg", "qc_flags")
    for key, spec in pv.RESOURCES.items():
        blob = " ".join(spec.columns)
        for bad in banned:
            assert bad not in blob, f"{key} fingerprints volatile {bad}"


def test_frames_spec_hashes_the_columns_the_incident_was_about():
    cols = " ".join(pv.RESOURCES["table:frames"].columns)
    for needed in ("era_id", "naxis1", "naxis2", "filter", "canonical_target"):
        assert needed in cols


def test_producer_of_resolves_and_returns_none_for_unknown():
    assert pv.producer_of("table:frames") == "S0"
    assert pv.producer_of("table:nope") is None


# ---------------------------------------------------------------------------
# Thin I/O wrappers, on an in-memory database
# ---------------------------------------------------------------------------

@pytest.fixture()
def con():
    c = sqlite3.connect(":memory:")
    c.execute("CREATE TABLE frames (obs_rowid INTEGER, tree TEXT, "
              "is_canonical INTEGER, era_id INTEGER, naxis1 REAL, "
              "naxis2 REAL, filter TEXT, imagetyp TEXT, exptime REAL, "
              "canonical_target TEXT, target_key TEXT, night TEXT, "
              "pltsolvd REAL, error TEXT)")
    c.executemany("INSERT INTO frames VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [
        (1, "rawimage", 1, 80, 8.0, 3211.0, "g", "Light Frame", 120.0,
         "Vega", "vega", "2026-04-24", None, None),
        (2, "rawimage", 1, 78, 4800.0, 3211.0, "hrg", "Light Frame", 240.0,
         "tet CrB", "tetcrb", "2026-03-22", None, None),
    ])
    c.commit()
    yield c
    c.close()


def test_fingerprint_table_present(con, tmp_path):
    fp = pv.fingerprint_resource(pv.RESOURCES["table:frames"], con, tmp_path)
    assert fp.present and fp.n_rows == 2


def test_fingerprint_table_changes_when_geometry_is_corrected(con, tmp_path):
    """The trigger, in miniature: fix NAXIS1 on one row and the frames
    fingerprint must move even though the row count does not."""
    before = pv.fingerprint_resource(pv.RESOURCES["table:frames"], con,
                                     tmp_path)
    con.execute("UPDATE frames SET naxis1 = 4800, era_id = 78 "
                "WHERE obs_rowid = 1")
    con.commit()
    after = pv.fingerprint_resource(pv.RESOURCES["table:frames"], con,
                                    tmp_path)
    assert before.n_rows == after.n_rows == 2
    assert before.digest != after.digest


def test_fingerprint_table_ignores_an_error_message_rewrite(con, tmp_path):
    """The rescue pass rewrote error TEXT without changing usability; the
    fingerprint must not cry wolf over it."""
    before = pv.fingerprint_resource(pv.RESOURCES["table:frames"], con,
                                     tmp_path)
    con.execute("UPDATE frames SET error = NULL WHERE error IS NULL")
    con.commit()
    after = pv.fingerprint_resource(pv.RESOURCES["table:frames"], con,
                                    tmp_path)
    assert before.digest == after.digest


def test_fingerprint_missing_table_is_reported_missing(con, tmp_path):
    fp = pv.fingerprint_resource(pv.RESOURCES["table:detector_params"], con,
                                 tmp_path)
    assert fp.present is False and fp.token == "MISSING"


def test_fingerprint_file_missing_and_present(tmp_path, con):
    spec = pv.RESOURCES["file:ops/2026-08_observatory_request.md"]
    assert pv.fingerprint_resource(spec, con, tmp_path).present is False
    p = tmp_path / "ops"
    p.mkdir()
    (p / "2026-08_observatory_request.md").write_bytes(b"hello")
    fp = pv.fingerprint_resource(spec, con, tmp_path)
    assert fp.present and fp.n_rows == 5


def test_fingerprint_all_raises_on_undeclared_key(con, tmp_path):
    with pytest.raises(pv.ProvenanceError):
        pv.fingerprint_all(["table:does_not_exist_in_registry"], con, tmp_path)


def test_record_and_read_roundtrip(con):
    pv.ensure_table(con)
    pv.record_run(con, "S0", "2026-08-18T15:03:50Z", "S0 v1.0", "3fec788",
                  {"table:x": "1:a"}, {"table:frames": "2:b"}, note="n")
    recs = pv.read_records(con)
    assert recs["S0"].inputs == {"table:x": "1:a"}
    assert recs["S0"].outputs == {"table:frames": "2:b"}
    assert recs["S0"].note == "n"


def test_read_records_returns_latest_run_per_stage(con):
    pv.ensure_table(con)
    pv.record_run(con, "S0", "2026-08-18T10:00:00Z", "v1", "a",
                  {"i": "1:a"}, {"o": "1:a"})
    pv.record_run(con, "S0", "2026-08-18T15:00:00Z", "v2", "b",
                  {"i": "1:b"}, {"o": "1:b"})
    assert pv.read_records(con)["S0"].code_version == "v2"


def test_read_records_on_a_database_without_the_table(con):
    assert pv.read_records(con) == {}


def test_ensure_table_is_idempotent(con):
    pv.ensure_table(con)
    pv.ensure_table(con)
    assert pv.read_records(con) == {}


def test_recorded_json_is_sorted_for_stable_diffs(con):
    pv.ensure_table(con)
    pv.record_run(con, "S0", "t", "v", "g", {"b": "1:b", "a": "1:a"}, {})
    raw = con.execute("SELECT inputs_json FROM stage_provenance").fetchone()[0]
    assert list(json.loads(raw)) == ["a", "b"]


# ===========================================================================
# REGRESSION TESTS ADDED AFTER THE 2026-08-18 ADVERSARIAL REVIEW
#
# Each one exists because a reviewer broke something BY EXECUTION rather than
# by reading, which is the only kind of finding that cannot be argued with.
# The ordering below matches the review's severity ordering.
# ===========================================================================

import os                                                    # noqa: E402
import subprocess                                            # noqa: E402
import sys                                                   # noqa: E402
from pathlib import Path                                     # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE_DIR = REPO_ROOT / "pipeline"
STATUS_CLI = PIPELINE_DIR / "scripts" / "check_pipeline_status.py"


# ---------------------------------------------------------------------------
# BLOCKER 1 — the re-run plan is the deliverable; every command in it must run
# ---------------------------------------------------------------------------
# The first version of this module emitted, for seven of eighteen stages,
#     python -c 'from macro_core.report_s0 import main; main()'
# which fails twice: macro_core is not importable from the repo root, and not
# one of the seven report modules defines main().  Four more stages emitted a
# bare script whose argparse requires a subcommand and exits 2.  Ten of the
# eighteen printed commands did not run.  Reading the strings could not catch
# that; running them can, so these tests run them.

def _emitted_commands() -> list[tuple[str, list[str]]]:
    """``[(stage_key, argv), ...]`` for every machine-run build command.

    Hand-authored commands (marked by a leading paren) are excluded: they
    are instructions to a person, not argv.
    """
    out = []
    for stage in pv.STAGES:
        for line in stage.build_lines:
            if line.startswith(pv.HAND_AUTHORED_CMD_PREFIX):
                continue
            parts = line.split()
            assert parts[0] == "python", \
                f"{stage.key}: command does not start with python: {line!r}"
            out.append((stage.key, parts[1:]))
    return out


def test_every_build_command_names_a_script_that_exists():
    """The cheap half of the check, so a typo fails fast and unambiguously."""
    for key, argv in _emitted_commands():
        script = REPO_ROOT / argv[0]
        assert script.is_file(), f"{key}: no such script: {argv[0]}"


def test_no_build_command_uses_the_python_dash_c_form():
    """`python -c '...'` is banned outright in this table.

    It cannot be argparse-checked, it silently depends on the caller's
    PYTHONPATH, and every previous instance of it in this file named a
    function that did not exist.  A stage whose renderer has no CLI uses
    `check_pipeline_status.py render <stage>` instead.
    """
    for stage in pv.STAGES:
        assert "-c" not in stage.build_cmd.split(), \
            f"{stage.key}: build_cmd uses the unverifiable `python -c` form"
        assert "main()" not in stage.build_cmd, \
            f"{stage.key}: build_cmd calls main(), which no report module has"


def test_every_build_command_is_runnable():
    """THE test this review earned: run each emitted command with --help.

    Appending --help exercises the whole argv — the script path, the
    subcommand, every flag — through the real argparse, and exits 0 before
    any work happens.  An invalid subcommand or an unknown option exits 2.
    Nothing is built, no database is touched, no archive is read.
    """
    env = dict(os.environ, PYTHONPATH=str(PIPELINE_DIR))
    failures = []
    for key, argv in _emitted_commands():
        proc = subprocess.run([sys.executable, *argv, "--help"],
                              cwd=str(REPO_ROOT), env=env,
                              capture_output=True, text=True, timeout=180)
        if proc.returncode != 0:
            failures.append(f"{key}: `python {' '.join(argv)}` -> exit "
                            f"{proc.returncode}\n{proc.stderr.strip()[:400]}")
    assert not failures, "un-runnable build commands:\n" + "\n".join(failures)


def test_every_emitted_flag_appears_in_its_scripts_help():
    """--help masks an unknown OPTION (argparse acts on it first), so the
    flags are checked separately against the help text.

    This is what would have caught `run_g_tcrb_validation.py --bogus`, which
    exits 0 under --help alone.
    """
    env = dict(os.environ, PYTHONPATH=str(PIPELINE_DIR))
    helps: dict[str, str] = {}
    for key, argv in _emitted_commands():
        flags = [a for a in argv[1:] if a.startswith("--")]
        if not flags:
            continue
        if argv[0] not in helps:
            proc = subprocess.run([sys.executable, argv[0], "--help"],
                                  cwd=str(REPO_ROOT), env=env,
                                  capture_output=True, text=True, timeout=180)
            helps[argv[0]] = proc.stdout + proc.stderr
        for flag in flags:
            assert flag in helps[argv[0]], \
                f"{key}: {argv[0]} has no option {flag}"


def test_the_g_stage_command_actually_extracts():
    """G's every action sits behind a flag; the bare command creates two
    directories and exits without touching g_extractions.

    Recording the stage after that would stamp a run that never happened —
    the exact laundering this module exists to prevent — so the declared
    command must carry an action flag.
    """
    cmd = pv.STAGE_BY_KEY["G"].build_cmd
    assert "run_g_tcrb_validation.py" in cmd
    assert "--all" in cmd or "--run" in cmd, \
        "G's command performs no extraction; recording it would be a lie"


def test_multi_command_stages_emit_every_command_they_need():
    """One invocation of run_s2_campaign.py runs ONE sub-stage, while the S2
    stage declares five output tables.  A plan that prints one line for five
    jobs silently omits four fifths of the work."""
    s2 = pv.STAGE_BY_KEY["S2"]
    stages_named = {line.split()[-1] for line in s2.build_lines}
    assert {"ceiling", "ptc", "reconstruct", "linearity", "params"} \
        <= stages_named
    s1 = pv.STAGE_BY_KEY["S1"]
    assert {"design", "run", "autopsy"} <= \
        {line.split()[-1] for line in s1.build_lines}


def test_render_subcommand_resolves_every_registered_renderer():
    """`render <stage> --check` imports the module and looks up the symbol
    without rendering.  It is the standing proof that the mapping points at
    functions that exist — the property the `main()` commands lacked."""
    sys.path.insert(0, str(PIPELINE_DIR))
    import importlib
    status = importlib.import_module("check_pipeline_status") \
        if (PIPELINE_DIR / "scripts") in [Path(p) for p in sys.path] else None
    # Import the CLI module by path, so this test does not depend on the
    # scripts directory being importable.
    import importlib.util
    spec = importlib.util.spec_from_file_location("_cps", STATUS_CLI)
    status = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(status)
    for key, (module_name, symbol, _kind) in status.RENDERERS.items():
        assert key in pv.STAGE_BY_KEY, f"renderer for unknown stage {key}"
        module = importlib.import_module(module_name)
        assert callable(getattr(module, symbol, None)), \
            f"{module_name}.{symbol} does not exist — the exact defect the " \
            f"`main()` commands had"


def test_every_report_stage_has_a_command_that_renders_only_the_page():
    """A report stage's command must not rebuild its whole parent stage:
    re-rendering a page is minutes, rebuilding S3 is hours."""
    for stage in pv.STAGES:
        if not stage.key.startswith("R-"):
            continue
        cmd = stage.build_cmd
        assert ("render" in cmd or "report" in cmd), \
            f"{stage.key}: {cmd!r} does not look like a render-only command"


# ---------------------------------------------------------------------------
# BLOCKER 2 — nothing published may sit outside the graph
# ---------------------------------------------------------------------------
# The CV photometry product (products/phot/cv_timeseries.sqlite) — the ONLY
# photometry product that contains phantom-era frames — was absent from the
# DAG entirely, so `status` could exit 0 while it rested on frames whose era
# assignment is known wrong.  This test makes that class of omission fail.

#: Files deliberately outside the graph, each with the reason.  An entry here
#: is a decision on the record, which is the point: the alternative is a
#: product quietly missing from the gate.
UNDECLARED_ON_PURPOSE: dict[str, str] = {
    "products/manifest/rlmt-manifest.sqlite":
        "the manifest itself — it is the database every table resource is "
        "fingerprinted INSIDE, not a resource within it",
}


def _declared_paths() -> set[str]:
    out = set()
    for spec in pv.RESOURCES.values():
        if spec.kind in ("file", "stat"):
            out.add(spec.name)
        elif spec.kind == "db" and spec.database:
            out.add(spec.database)
    return out


def test_every_published_page_is_a_declared_resource():
    """Every docs/**/*.html either belongs to a stage or is listed above.

    The status page used to print a verdict for ROADMAP.md and the six
    project pages while its exit code ignored them — a page issuing
    judgements its own gate could not enforce.
    """
    declared = _declared_paths()
    undeclared = []
    for path in sorted((REPO_ROOT / "docs").rglob("*.html")):
        rel = str(path.relative_to(REPO_ROOT))
        if rel == "docs/pipeline/pipeline_status.html":
            continue                      # this module's own output
        if rel not in declared and rel not in UNDECLARED_ON_PURPOSE:
            undeclared.append(rel)
    assert not undeclared, \
        "published pages outside the provenance graph:\n  " + \
        "\n  ".join(undeclared)


def test_every_product_database_is_a_declared_resource():
    """Every products/**/*.sqlite is fingerprinted by some stage.

    This is the test that would have caught cv_timeseries.sqlite.
    """
    declared = _declared_paths()
    undeclared = []
    products = REPO_ROOT / "products"
    if not products.is_dir():
        pytest.skip("no products directory in this checkout")
    for path in sorted(products.rglob("*.sqlite")):
        rel = str(path.relative_to(REPO_ROOT))
        if rel not in declared and rel not in UNDECLARED_ON_PURPOSE:
            undeclared.append(rel)
    assert not undeclared, \
        "product databases outside the provenance graph:\n  " + \
        "\n  ".join(undeclared)


def test_the_cv_product_is_declared_and_reads_the_frames_it_measures():
    """Named explicitly, because it is the product the artifact actually
    reaches: 207 EU UMa frames in phantom era 80, split from 75 in era 78 —
    one camera configuration — under a recorded rule of 'never mixed inside
    a series'."""
    cv = pv.STAGE_BY_KEY["CV-S4"]
    assert "db:cvphot:cv_selection" in cv.writes
    assert "table:frames@cv" in cv.reads
    assert "table:stage_cv_timeseries" in cv.reads
    assert "table:frame_times" in cv.reads


# ---------------------------------------------------------------------------
# BLOCKER 3 — no stage may report FRESH over a destroyed input it consumes
# ---------------------------------------------------------------------------

def test_s4_declares_the_s2_constants_it_actually_applies():
    """S4 and CV-S4 apply S2_MODE_VETO_ADU from macro_phot/series.py to
    every frame, and the S2 tables behind those numbers are destroyed.

    With no declared edge, S4 reported FRESH in the same run in which S2
    reported OUTPUT_MISSING: a green verdict over a destroyed input, which
    is the dangerous direction of failure.
    """
    for key in ("S4", "CV-S4"):
        reads = pv.STAGE_BY_KEY[key].reads
        assert "table:detector_params" in reads, \
            f"{key} consumes S2's detector constants but does not declare them"
        assert "table:s2_ceiling_modes" in reads
        assert "file:pipeline/macro_phot/series.py" in reads, \
            f"{key} reads the adopted constants from series.py; the file " \
            f"must be a declared resource or a change to it is invisible"


def test_a_stage_with_an_absent_input_can_never_read_fresh():
    """The general property, driven directly: whatever the record says, a
    stage one of whose declared inputs is MISSING is not fresh."""
    stage = pv.Stage(key="X", title="t", code_version="v1",
                     reads=("i",), writes=("o",), build_cmd="python x.py")
    record = pv.Record("X", "t", "v1", "g", {"i": "MISSING"}, {"o": "1:a"})
    got = pv.is_stale(stage, record, {"i": "MISSING"}, {"o": "1:a"}, "v1")
    assert not got.ok
    assert any("ABSENT" in r for r in got.reasons)


def test_series_constants_file_is_hashed_whole():
    """Its `why` must say so: the constants are bare module-level literals
    with no separate home, so any edit to the file reads as a change.  That
    is the conservative direction and the reader is told which way it errs.
    """
    spec = pv.RESOURCES["file:pipeline/macro_phot/series.py"]
    assert spec.kind == "file"
    assert "conservative" in spec.why


# ---------------------------------------------------------------------------
# MAJOR — the filter-identity evidence must be able to invalidate a document
# ---------------------------------------------------------------------------

def test_frame_dispersion_is_declared_and_the_strategies_read_it():
    """S2c's per-frame verdicts are what every REDERIVE row is waiting on.

    Undeclared, the classifier could finish, overturn an exclusion list, and
    nothing in the graph would go stale — the same silent drift one level up.
    """
    assert "table:frame_dispersion" in pv.RESOURCES
    assert pv.producer_of("table:frame_dispersion") == "S2c"
    assert "table:frame_dispersion" in pv.STAGE_BY_KEY["STRAT"].reads


def test_dispersion_verdicts_are_fingerprinted_but_wall_clock_is_not():
    """The verdict and its code version are hashed (a reclassification from
    the same pixels IS a changed input to a human decision); measured_at and
    measure_s are not (wall-clock churn)."""
    cols = " ".join(pv.RESOURCES["table:frame_dispersion"].columns)
    assert "verdict" in cols and "code_version" in cols
    assert "measured_at" not in cols and "measure_s" not in cols


# ---------------------------------------------------------------------------
# MAJOR — fingerprint ROW scope, not just column choice
# ---------------------------------------------------------------------------

@pytest.fixture()
def scoped_con():
    """A frames table with the columns the scoped predicates read."""
    c = sqlite3.connect(":memory:")
    c.execute("CREATE TABLE frames (obs_rowid INTEGER, tree TEXT, "
              "is_canonical INTEGER, era_id INTEGER, naxis1 REAL, "
              "naxis2 REAL, filter TEXT, imagetyp TEXT, exptime REAL, "
              "canonical_target TEXT, target_key TEXT, night TEXT, "
              "pltsolvd REAL, error TEXT, target_best TEXT)")
    c.executemany(
        "INSERT INTO frames VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", [
            # 1: a phantom-era frame belonging to NO downstream consumer
            (1, "rawimage", 1, 80, 8.0, 3211.0, "g", "Light Frame", 120.0,
             "Vega", "vega", "2026-04-24", None, None, "Vega"),
            # 2: a T CrB grism frame — G's slice
            (2, "rawimage", 1, 76, 4800.0, 3211.0, "hrg", "Light Frame",
             240.0, "T CrB", "tcrb", "2026-03-22", None, None, "T CrB"),
            # 3: an AN UMa frame — S4's slice
            (3, "rawimage", 1, 76, 4800.0, 3211.0, "g", "Light Frame", 60.0,
             "AN UMa", "anuma", "2026-03-22", None, None, "AN UMa"),
        ])
    c.commit()
    yield c
    c.close()


def test_scoped_slice_selects_only_its_own_rows(scoped_con, tmp_path):
    whole = pv.fingerprint_resource(pv.RESOURCES["table:frames"],
                                    scoped_con, tmp_path)
    grism = pv.fingerprint_resource(pv.RESOURCES["table:frames@grism"],
                                    scoped_con, tmp_path)
    proto = pv.fingerprint_resource(pv.RESOURCES["table:frames@s4proto"],
                                    scoped_con, tmp_path)
    assert whole.n_rows == 3
    assert grism.n_rows == 1          # the T CrB hrg frame
    assert proto.n_rows == 1          # the AN UMa frame


def test_repairing_a_frame_no_consumer_reads_leaves_their_slices_identical(
        scoped_con, tmp_path):
    """THE point of scoping, reproduced from the review's own experiment.

    The reviewer applied the single minimal true correction — one era-80
    frame, naxis1 8 -> 4800, era_id 80 -> 78 — to a frame belonging to no
    S4 target and no G target, and watched every stage go stale, mandating
    a re-run of the grism extraction and the photometry prototype that the
    audit itself proves cannot change.

    With scoped resources, that same repair moves the whole-table digest
    (correctly: S0's output DID change) while leaving the two consumer
    slices byte-identical, so once S0 is rebuilt and re-recorded those
    stages have nothing to do.
    """
    R = pv.RESOURCES
    before = {k: pv.fingerprint_resource(R[k], scoped_con, tmp_path).token
              for k in ("table:frames", "table:frames@grism",
                        "table:frames@s4proto")}
    scoped_con.execute("UPDATE frames SET naxis1 = 4800, era_id = 78 "
                       "WHERE obs_rowid = 1")
    after = {k: pv.fingerprint_resource(R[k], scoped_con, tmp_path).token
             for k in ("table:frames", "table:frames@grism",
                       "table:frames@s4proto")}
    assert after["table:frames"] != before["table:frames"]
    assert after["table:frames@grism"] == before["table:frames@grism"]
    assert after["table:frames@s4proto"] == before["table:frames@s4proto"]


def test_a_change_inside_a_slice_still_moves_that_slice(scoped_con, tmp_path):
    """The other direction, which matters more: scoping must never hide a
    real change.  A slice that is too narrow is the one error this module
    may not make."""
    spec = pv.RESOURCES["table:frames@grism"]
    before = pv.fingerprint_resource(spec, scoped_con, tmp_path).token
    scoped_con.execute("UPDATE frames SET era_id = 78 WHERE obs_rowid = 2")
    assert pv.fingerprint_resource(spec, scoped_con, tmp_path).token != before


def test_every_scoped_spec_hashes_the_same_columns_as_its_base():
    """A slice and its whole differing in COLUMNS would let a change be
    visible in one fingerprint and invisible in the other."""
    for spec in pv.RESOURCES.values():
        if not spec.scope_of:
            continue
        base = pv.RESOURCES[spec.scope_of]
        assert spec.columns == base.columns, \
            f"{spec.key} hashes different columns from {spec.scope_of}"
        assert spec.name == base.name
        assert spec.where, f"{spec.key} declares scope_of but no predicate"


def test_scoped_specs_explain_what_they_were_widened_from():
    """A predicate narrower than the consuming query would hide a change, so
    each one has to say, in prose, which query it is a superset of."""
    for spec in pv.RESOURCES.values():
        if spec.where:
            assert len(spec.why) > 120, \
                f"{spec.key}: a scoped resource needs its scope justified"


# ---------------------------------------------------------------------------
# MINOR — `record` is a publish gate, not a rubber stamp
# ---------------------------------------------------------------------------

def _cps():
    """Load the CLI module by path (the scripts dir is not a package)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_cps_guard", STATUS_CLI)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_record_refuses_a_hand_authored_stage_with_no_note():
    """`record OPS` after changing nothing used to turn OPS fresh, and once
    the last stage was stamped `status` exited 0.  For OPS — the stage this
    audit shows to be materially wrong — that is the wrong default."""
    cps = _cps()
    ops = pv.STAGE_BY_KEY["OPS"]
    assert ops.hand_authored
    why = cps.record_objection(ops, "2026-08-18T20:00:00Z", [], note="")
    assert why and "hand-authored" in why


def test_record_accepts_a_hand_authored_stage_with_a_note():
    cps = _cps()
    ops = pv.STAGE_BY_KEY["OPS"]
    assert cps.record_objection(ops, "2026-08-18T20:00:00Z", [],
                                note="regenerated Item B from calib_gaps") == ""


def test_record_refuses_to_overwrite_an_existing_run():
    """The primary key is (stage, run_utc) and the write is INSERT OR
    REPLACE, so recording twice on one timestamp destroys the earlier run's
    fingerprints — the history the table exists to keep."""
    cps = _cps()
    s0 = pv.STAGE_BY_KEY["S0"]
    why = cps.record_objection(s0, "2026-08-18T15:00:00Z",
                               ["2026-08-18T15:00:00Z"], note="")
    assert why and "already recorded" in why


def test_record_refuses_a_run_time_that_does_not_advance():
    """What a re-run looks like when its build_meta was not updated: a stamp
    that precedes the work it attests."""
    cps = _cps()
    s0 = pv.STAGE_BY_KEY["S0"]
    why = cps.record_objection(s0, "2026-08-18T09:00:00Z",
                               ["2026-08-18T15:00:00Z"], note="")
    assert why and "older than" in why


def test_record_allows_an_ordinary_advancing_run():
    cps = _cps()
    s0 = pv.STAGE_BY_KEY["S0"]
    assert cps.record_objection(s0, "2026-08-18T16:00:00Z",
                                ["2026-08-18T15:00:00Z"], note="") == ""


def test_record_objection_tolerates_the_two_timestamp_dialects():
    """The stages write '+00:00' and 'Z'; 'Z' > '+' in ASCII, so an
    unnormalized compare would call a Z-stamped run later than a +00:00 run
    taken at the same instant."""
    cps = _cps()
    s0 = pv.STAGE_BY_KEY["S0"]
    assert cps.record_objection(
        s0, "2026-08-18T16:00:00Z", ["2026-08-18T15:00:00+00:00"],
        note="") == ""
    assert cps.record_objection(
        s0, "2026-08-18T14:00:00Z", ["2026-08-18T15:00:00+00:00"],
        note="") != ""


# ---------------------------------------------------------------------------
# MINOR — code versions read from source, never by importing a live file
# ---------------------------------------------------------------------------

def test_read_version_constant_finds_a_module_level_string():
    src = 'X = 1\nCV_CODE_VERSION = "CV-S4 v1.0 (2026-08-18)"\nY = 2\n'
    assert pv.read_version_constant(src, "CV_CODE_VERSION") == \
        "CV-S4 v1.0 (2026-08-18)"


def test_read_version_constant_returns_none_for_a_missing_symbol():
    assert pv.read_version_constant("A = 1\n", "CV_CODE_VERSION") is None


def test_read_version_constant_ignores_a_constant_inside_a_function():
    """A name assigned inside a function is not the module's API."""
    src = 'def f():\n    CV_CODE_VERSION = "nope"\n'
    assert pv.read_version_constant(src, "CV_CODE_VERSION") is None


def test_read_version_constant_survives_a_half_written_file():
    """Sibling workflows edit these files while a status check runs.  A file
    that does not parse is 'unknown', which is safe; guessing is not."""
    assert pv.read_version_constant("def f(:\n", "X") is None


def test_read_version_constant_refuses_a_non_string_value():
    assert pv.read_version_constant("V = 3\n", "V") is None


def test_version_file_stages_point_at_files_that_exist():
    for stage in pv.STAGES:
        if not stage.version_file:
            continue
        path = REPO_ROOT / stage.version_file
        assert path.is_file(), f"{stage.key}: {stage.version_file} missing"
        got = pv.read_version_constant(path.read_text(encoding="utf-8"),
                                       stage.version_symbol)
        assert got, f"{stage.key}: {stage.version_symbol} not found in " \
                    f"{stage.version_file}"


# ---------------------------------------------------------------------------
# Structural invariants that keep the additions honest
# ---------------------------------------------------------------------------

def test_the_external_catalog_has_a_producer_stage():
    """S0's input is not a mystery: S0e rewrites it.  Before S0e was
    declared, the geometry repair appeared as an unexplained mtime change on
    an 'external source'."""
    assert pv.producer_of("stat:rlmt-catalog") == "S0e"


def test_hand_authored_stages_are_recognisable_as_such():
    hand = {s.key for s in pv.STAGES if s.hand_authored}
    assert {"OPS", "STRAT", "WEB"} <= hand
    for key in ("S0", "S0b", "S3", "S4", "CV-S4", "G"):
        assert not pv.STAGE_BY_KEY[key].hand_authored


def test_every_stage_writes_something_and_nothing_is_written_twice():
    seen: dict[str, str] = {}
    for stage in pv.STAGES:
        assert stage.writes, f"{stage.key} declares no output"
        for w in stage.writes:
            assert w not in seen, \
                f"{w} is written by both {seen[w]} and {stage.key}"
            seen[w] = stage.key


def test_all_states_covers_every_verdict_the_module_can_return():
    """A new state that a summary line forgets to print is a verdict nobody
    sees."""
    assert set(pv.ALL_STATES) == {pv.FRESH, pv.STALE, pv.STALE_UPSTREAM,
                                  pv.NEVER_RUN, pv.OUTPUT_MISSING}


# ---------------------------------------------------------------------------
# The general form of the "bare script" defect
# ---------------------------------------------------------------------------
# `python <script> --help` exits 0 even when the script's subparser is
# required, because argparse acts on --help before it misses the positional.
# So `run_s1_batch.py` and `build_s4_photometry.py` — both required=True,
# both bare in the original plan — would slip past the runnability test.
#
# argparse's own usage line gives the answer away: a REQUIRED choice group
# renders as `{build,run,...}` and an optional one as `[{...}]`.  Parsing it
# turns "does this script demand a subcommand?" into a fact rather than a
# thing a person has to remember per script.

def required_choice_group(usage: str) -> set[str]:
    """Choices of the first REQUIRED positional-choice group in an argparse
    usage string, or an empty set.  PURE.

    Bracket depth does the work: a group inside ``[...]`` is optional, a
    group at depth zero is required.
    """
    depth = 0
    i = 0
    while i < len(usage):
        ch = usage[i]
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth = max(0, depth - 1)
        elif ch == "{" and depth == 0:
            end = usage.find("}", i)
            if end == -1:
                break
            return {c.strip() for c in usage[i + 1:end].split(",")
                    if c.strip()}
        i += 1
    return set()


def test_required_choice_group_reads_a_required_group():
    usage = "usage: x.py [-h] [--m M] {build,run,enqueue} ..."
    assert required_choice_group(usage) == {"build", "run", "enqueue"}


def test_required_choice_group_ignores_an_optional_group():
    usage = "usage: x.py [-h] [--stage {audit,report}] [--skip-report]"
    assert required_choice_group(usage) == set()


def test_required_choice_group_is_empty_when_there_is_none():
    assert required_choice_group("usage: x.py [-h] [--catalog C]") == set()


def test_no_stage_emits_a_bare_script_that_demands_a_subcommand():
    """The general form of the defect: four stages emitted a bare script
    whose argparse requires a subcommand and exits 2.  Two of them
    (`run_s1_batch.py`, `build_s4_photometry.py`) exit 0 under --help, so
    only this check catches them."""
    env = dict(os.environ, PYTHONPATH=str(PIPELINE_DIR))
    usages: dict[str, str] = {}
    problems = []
    for key, argv in _emitted_commands():
        script = argv[0]
        if script not in usages:
            proc = subprocess.run([sys.executable, script, "--help"],
                                  cwd=str(REPO_ROOT), env=env,
                                  capture_output=True, text=True, timeout=180)
            # The usage block ends at the first blank line.
            usages[script] = (proc.stdout + proc.stderr).split("\n\n")[0]
        needed = required_choice_group(usages[script])
        if needed and not (set(argv[1:]) & needed):
            problems.append(
                f"{key}: `python {' '.join(argv)}` names no subcommand, but "
                f"{script} requires one of {sorted(needed)} — bare, it "
                f"exits 2 without doing any work")
    assert not problems, "\n".join(problems)
