#!/usr/bin/env python3
"""update_project_plan.py — the command you run as work completes.

THE WORKING RHYTHM THIS EXISTS FOR

    publish the plan  ->  do the work  ->  update the page the moment a step
    completes  ->  plan the next step  ->  repeat

The third arrow is the one that always breaks, because updating a page by
hand is a chore nobody does at the moment of completion — and a page updated
later is a page written from memory.  So the update is one command:

    python pipeline/scripts/update_project_plan.py set CV-P2-stlmi done \\
        --evidence docs/pipeline/s4_photometry.html --note 'ensemble solved'
    python pipeline/scripts/update_project_plan.py render

SUBCOMMANDS
-----------
``show [project]``  the plan and current status, as text: ``done`` of the
                    IN-SCOPE tasks, the dropped count beside it, and the
                    2027 backlog listed apart.  It states WHEN the statuses
                    were last synced — a count printed from statuses last
                    reconciled six weeks ago is a count about six weeks ago.
                    ``--live`` evaluates the DAG and the database probes now
                    (about a minute) and prints the statuses a sync would
                    leave, so the count is never an unsynced one.
                    ``--history`` replays the recorded status changes.
``set <task-id> <status> [--evidence X] [--note Y] [--ruling Z]``
                    record a status change, stamped with the UTC time.
                    Append-only: nothing is overwritten, so "when did this
                    become true?" always has an answer.  ``dropped`` and
                    ``deferred`` are RULINGS: they are refused unless the
                    ledger's task carries one or ``--ruling`` cites it.
``sync``            recompute the statuses the DATABASE or a dated RULING
                    has already decided: a ``done`` task whose stage is no
                    longer FRESH flips to ``redo_needed`` (and back when it
                    returns); a status recorded BEFORE a committee ruling on
                    its task gives way to the status that ruling set.
                    ``--dry-run`` prints without writing.
``blockers [project]``
                    every blocked or gated task with what holds it —
                    dependencies read from the status table, stage gates
                    from the DAG, probes from the manifest — beside whatever
                    prose it still carries.  Exits non-zero if a blocker's
                    prose is contradicted by its own probe.
``amendments <project>``
                    the markdown table of every task a committee ruling
                    added, changed or closed — emitted from the ledger, for
                    the strategy's "Committee amendments" section.
``titles [--write]``
                    compare each manuscript skeleton's ``\\title{}`` with the
                    title the ledger holds; exit non-zero on a mismatch.
                    ``--write`` re-emits the title line and nothing else.
``render``          regenerate every project page and the hub, then
                    RECORD the WEB provenance stage it just satisfied.

WRITE DISCIPLINE
----------------
Only ``project_plan_status`` is ever written, only by ``set`` and ``sync``,
inside one transaction, with ``busy_timeout = 300000`` so a concurrent solve
or build waits rather than fails.  ``show`` and ``render`` open the manifest
READ-ONLY.
"""

from __future__ import annotations

import argparse
import sqlite3
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent
REPO_ROOT = PIPELINE_DIR.parent
sys.path.insert(0, str(PIPELINE_DIR))

from macro_core import project_plan as pp                    # noqa: E402
from macro_core import provenance as pv                      # noqa: E402

DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

_COLOR = {
    pp.DONE: "\033[32m",
    pp.IN_PROGRESS: "\033[36m",
    pp.REDO_NEEDED: "\033[33m",
    pp.BLOCKED: "\033[31m",
    pp.PENDING: "\033[90m",
    pp.DROPPED: "\033[2m",
    pp.DEFERRED: "\033[35m",
}
_RESET = "\033[0m"


def _paint(status: str, text: str) -> str:
    if not sys.stdout.isatty():
        return text
    return f"{_COLOR.get(status, '')}{text}{_RESET}"


def open_manifest(path: Path, read_only: bool) -> sqlite3.Connection:
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=300.0)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


# ---------------------------------------------------------------------------
# show
# ---------------------------------------------------------------------------
def _live_state(con):
    """``(stage_states, probe_results)`` — the DAG verdicts and every
    probe's number, read now.  The expensive half of ``--live``."""
    freshness, _ = pp.stage_freshness(con, REPO_ROOT)
    states = {k: f.state for k, f in freshness.items()}
    return states, pp.run_probes(con, pp.ledger_tasks())


def _sync_age_line(con) -> str:
    """One line saying how old the statuses about to be printed are."""
    last = pp.last_sync_utc(con)
    if last is None:
        return ("statuses have never been synced against the DAG — run "
                "`sync`, or `show --live` for the synced view.")
    return (f"statuses last synced {last}. A stage that went stale since "
            f"is not reflected below — `show --live` evaluates it now.")


def _print_group(project, statuses, states, probes) -> None:
    """One group's plan: the in-scope count, then every phase."""
    counts = pp.status_counts(project.tasks, statuses)
    scope = pp.scope_summary(counts)
    print()
    print("=" * 78)
    print(f"{project.title}   {scope['done']}/{scope['in_scope']} "
          f"in-scope tasks complete")
    summary = "  ".join(f"{pp.STATUS_LABEL[s]}={counts[s]}"
                        for s in pp.ALL_STATUSES
                        if counts[s] and s not in pp.CLOSED_BY_RULING)
    print(f"  {summary}")
    # The two counts the fraction leaves out, always printed with it.
    print(f"  outside the count: {scope['dropped']} dropped by ruling, "
          f"{scope['deferred']} in the {pp.BACKLOG_NAME}")
    if project.ruling:
        print(f"  ruling: {project.ruling}")
    print(f"  venue: {project.venue}")
    print("=" * 78)
    for phase in project.phases:
        print(f"\n  {phase.name}")
        for t in phase.tasks:
            s = statuses[t.id]
            print(f"    {_paint(s, f'[{s:12}]')} {t.id:28} {t.title}")
            pad = f"{'':17}   "
            if s in pp.CLOSED_BY_RULING:
                print(f"{pad}ruling: {t.ruling or '(recorded by `set`)'}")
                continue
            if s == pp.BLOCKED and t.blocker:
                print(f"{pad}! {t.blocker}")
            if s in pp.OPEN_STATUSES or s == pp.BLOCKED:
                for line in pp.computed_blockers(t, statuses, states, probes):
                    print(f"{pad}> {line}")
            for gap in pp.acceptance_gaps(t, probes):
                flag = ("CONTRADICTED — recorded done, but"
                        if s == pp.DONE else "to go:")
                print(f"{pad}~ {flag} {gap}")
    nxt = pp.next_up(project.tasks, statuses, limit=3,
                     stage_states=states, probe_results=probes)
    if nxt:
        print("\n  NEXT UP:")
        for t in nxt:
            print(f"    - {t.id}: {t.title}")


def cmd_show(args) -> int:
    pp.validate()
    con = open_manifest(args.manifest, read_only=True)
    try:
        if args.history:
            rows = pp.read_history(con, args.project if args.project and
                                   args.project in
                                   {t.id for t in pp.ledger_tasks()} else None)
            if not rows:
                print("no status changes recorded yet.")
                return 0
            print(f"{'when':21} {'task':28} {'status':13} note")
            for task_id, status, evidence, note, when in rows:
                extra = note or evidence or ""
                print(f"{when:21} {task_id:28} {status:13} {extra}")
            return 0

        known = {g.key: g for g in pp.groups()}
        keys = [args.project] if args.project else list(known)
        for key in keys:
            if key not in known:
                print(f"unknown project: {key!r}; known: "
                      f"{', '.join(known)}", file=sys.stderr)
                return 2

        recorded = pp.read_statuses(con)
        tasks = pp.ledger_tasks()
        statuses = pp.overlay_statuses(tasks, recorded)
        states = probes = None
        if args.live:
            # The synced view: what `sync` would leave, without writing it.
            # `show` holds a read-only connection, so it cannot record the
            # changes — it applies them in memory and says how many there
            # are, which is the honest alternative to printing a count it
            # knows to be stale.
            states, probes = _live_state(con)
            changes = pp.derive_sync(tasks, statuses, states,
                                     pp.read_status_stamps(con))
            for c in changes:
                statuses[c.task_id] = c.new
            print(f"LIVE view: DAG and probes evaluated now. "
                  f"{len(changes)} recorded status(es) differ from what "
                  f"the database and the rulings say"
                  + (" — run `sync` to record them." if changes else "."))
        else:
            print(_sync_age_line(con))
        for key in keys:
            _print_group(known[key], statuses, states, probes)
        return 0
    finally:
        con.close()


# ---------------------------------------------------------------------------
# blockers
# ---------------------------------------------------------------------------
def cmd_blockers(args) -> int:
    """Every blocked or gated task, with what holds it — computed.

    The audit this exists for: five tasks sat blocked on "S2's tables were
    destroyed" for six weeks after the tables were rebuilt, because the
    blocker was a sentence and a sentence does not re-read the database.
    Here each gate is evaluated when it is shown, and a blocked task whose
    probes are all MET while it has no other gate is reported as
    CONTRADICTED — its prose says "blocked" and the database says "open".
    So is the mirror image: a task recorded ``done`` whose acceptance probe
    the manifest does not satisfy.
    """
    pp.validate()
    con = open_manifest(args.manifest, read_only=True)
    try:
        known = {g.key: g for g in pp.groups()}
        keys = [args.project] if args.project else list(known)
        for key in keys:
            if key not in known:
                print(f"unknown project: {key!r}; known: "
                      f"{', '.join(known)}", file=sys.stderr)
                return 2
        recorded = pp.read_statuses(con)
        tasks = pp.ledger_tasks()
        statuses = pp.overlay_statuses(tasks, recorded)
        states = None
        if args.live:
            states, probes = _live_state(con)
        else:
            # Probes are cheap (one scalar query each); the DAG is not.
            probes = pp.run_probes(con, tasks)
            print("stage gates not evaluated (pass --live for the DAG "
                  "verdicts; it takes about a minute).")

        contradicted = 0
        for key in keys:
            group = known[key]
            lines = []
            for t in group.tasks:
                s = statuses[t.id]
                if s in pp.CLOSED_BY_RULING:
                    continue
                gaps = pp.acceptance_gaps(t, probes)
                if s == pp.DONE or s == pp.REDO_NEEDED:
                    # A finished task is only of interest here if the
                    # database disagrees that it is finished.
                    if s == pp.DONE and gaps:
                        contradicted += len(gaps)
                        lines.append(f"  {_paint(s, f'[{s:11}]')} {t.id}: "
                                     f"{t.title}")
                        for gap in gaps:
                            lines.append(
                                f"      CONTRADICTED: recorded done, but "
                                f"its acceptance probe is not met — {gap}")
                    continue
                gates = pp.computed_blockers(t, statuses, states, probes)
                gate_results = [r for r in probes.get(t.id, ())
                                if r.kind == "gate"]
                if s != pp.BLOCKED and not gates and not gaps:
                    continue
                lines.append(f"  {_paint(s, f'[{s:11}]')} {t.id}: {t.title}")
                for g in gates:
                    lines.append(f"      > {g}")
                for r in gate_results:
                    if r.met:
                        lines.append(f"      = gate met — {r}")
                for gap in gaps:
                    lines.append(f"      ~ to go: {gap}")
                if s == pp.BLOCKED:
                    lines.append(f"      ! {t.blocker}")
                    if (gate_results and all(r.met for r in gate_results)
                            and not gates):
                        contradicted += 1
                        lines.append(
                            "      CONTRADICTED: every database gate on "
                            "this task is met and nothing else holds it, "
                            "yet it is recorded as blocked.")
            if lines:
                print(f"\n{group.title}")
                print("\n".join(lines))
        if contradicted:
            print(f"\n{contradicted} blocker(s) contradicted by the "
                  f"database.", file=sys.stderr)
            return 1
        return 0
    finally:
        con.close()


# ---------------------------------------------------------------------------
# amendments
# ---------------------------------------------------------------------------
def cmd_amendments(args) -> int:
    """Print the markdown table of one group's ruled tasks.

    This is what goes under "Committee amendments" in a strategy document.
    It is emitted, not typed, for the reason every other table on this site
    is: the ids in a strategy section that tasks then CITE must be the ids
    the ledger actually holds, and a transcription is where they diverge.
    """
    pp.validate()
    try:
        group = pp.group_of(args.project)
    except pp.PlanError as exc:
        print(exc, file=sys.stderr)
        return 2
    print("| Task | Action | Ruling | What it now is | Accept / why |")
    print("|---|---|---|---|---|")
    for t in group.tasks:
        if t.ruling is None:
            continue
        action = pp.ruling_action(t)
        cite = t.ruling.section + (
            f" ({', '.join(t.ruling.findings)})" if t.ruling.findings else "")
        # Closed tasks state their reason in `produces`; open ones state
        # the test that closes them.
        produces = t.produces.split(" Accept: ")[0]
        why = t.accept or produces
        def cell(text):
            return text.replace("|", "\\|")
        print(f"| `{t.id}` | {action} | {cell(cite)} | {cell(t.title)} | "
              f"{cell(why)} |")
    return 0


# ---------------------------------------------------------------------------
# titles
# ---------------------------------------------------------------------------
def cmd_titles(args) -> int:
    """Hold each manuscript skeleton to the title the ledger carries.

    A title is the first claim a reader meets, and three skeletons were
    making claims their strategies had already withdrawn (ED.E6).  Projects
    with an empty ``paper_title`` are skipped — the CV manuscript is a
    hand-edited draft this command must never write to.
    """
    pp.validate()
    mismatched = 0
    for project in pp.PROJECTS:
        if not project.paper_title:
            continue
        path = REPO_ROOT / pp.manuscript_path(project)
        if not path.exists():
            # manuscripts/ is not tracked, so a clean checkout has none.
            print(f"{project.key}: no manuscript at "
                  f"{pp.manuscript_path(project)} — skipped")
            continue
        tex = path.read_text(encoding="utf-8")
        current = pp.manuscript_title(tex)
        if current == project.paper_title:
            print(f"{project.key}: title matches the ledger")
            continue
        mismatched += 1
        print(f"{project.key}: title differs from the ledger")
        print(f"    manuscript: {current}")
        print(f"    ledger:     {project.paper_title}")
        if args.write:
            path.write_text(pp.with_title(tex, project.paper_title),
                            encoding="utf-8")
            print(f"    rewrote the \\title line of "
                  f"{pp.manuscript_path(project)}")
    if mismatched and not args.write:
        print(f"\n{mismatched} title(s) differ; re-emit with "
              f"`titles --write`.", file=sys.stderr)
        return 1
    return 0


# ---------------------------------------------------------------------------
# set
# ---------------------------------------------------------------------------
def cmd_set(args) -> int:
    pp.validate()
    try:
        task = pp.task_by_id(args.task_id)
    except pp.PlanError as exc:
        print(f"{exc}\n\nKnown ids for a project: "
              f"update_project_plan.py show <Project>", file=sys.stderr)
        return 2
    if args.status not in pp.ALL_STATUSES:
        print(f"unknown status {args.status!r}; expected one of "
              f"{', '.join(pp.ALL_STATUSES)}", file=sys.stderr)
        return 2
    if args.status == pp.DONE and not args.evidence and not task.evidence:
        print("refusing to mark done with no evidence: pass --evidence "
              "<report page or product path>.\nA 'done' with nothing behind "
              "it is exactly the claim this machinery exists to prevent.",
              file=sys.stderr)
        return 2
    note = args.note or ""
    if args.status in pp.CLOSED_BY_RULING:
        # `dropped` and `deferred` are rulings, not opinions (SYNTHESIS §0).
        # The ledger's own ruling suffices; otherwise the citation must be
        # given here, and it is stored in the note so the history carries
        # the reason beside the status.
        cited = args.ruling or (str(task.ruling) if task.ruling else "")
        if not cited:
            print(f"refusing to mark {args.status} with no ruling: a task "
                  f"leaves the plan only by being done or by a committee "
                  f"ruling that is cited.\nPass --ruling '<document> "
                  f"<section> (<finding ids>)', or add the ruling to the "
                  f"task in macro_core.project_plan.", file=sys.stderr)
            return 2
        note = f"ruling: {cited}" + (f" — {note}" if note else "")

    con = open_manifest(args.manifest, read_only=False)
    try:
        with con:
            stamp = pp.record_status(
                con, args.task_id, args.status,
                evidence=args.evidence or task.evidence,
                note=note)
        print(f"{args.task_id} -> {args.status}  ({stamp})")
        print(f"  {task.project} / {task.phase} / {task.title}")
        if args.evidence:
            print(f"  evidence: {args.evidence}")
        print("\nRe-render the pages:  python "
              "pipeline/scripts/update_project_plan.py render")
        return 0
    finally:
        con.close()


# ---------------------------------------------------------------------------
# sync
# ---------------------------------------------------------------------------
def cmd_sync(args) -> int:
    pp.validate()
    con = open_manifest(args.manifest, read_only=args.dry_run)
    try:
        freshness, _ = pp.stage_freshness(con, REPO_ROOT)
        states = {k: f.state for k, f in freshness.items()}
        recorded = pp.read_statuses(con)
        # The WHOLE ledger — the shared foundation's tasks rest on stages
        # too — and with the stamps, so a status recorded before a
        # committee ruling gives way to it (sync rule 3).
        tasks = pp.ledger_tasks()
        statuses = pp.overlay_statuses(tasks, recorded)
        changes = pp.derive_sync(tasks, statuses, states,
                                 pp.read_status_stamps(con))

        summary = {}
        for key, f in freshness.items():
            summary[f.state] = summary.get(f.state, 0) + 1
        print("stage verdicts: " + "  ".join(
            f"{s}={summary[s]}" for s in pv.ALL_STATES if s in summary))

        if not changes:
            print("nothing to sync — every recorded status already agrees "
                  "with the database.")
            return 0

        print(f"\n{len(changes)} status change(s) the database or a dated "
              f"ruling has already decided:")
        for c in changes:
            task = pp.task_by_id(c.task_id)
            print(f"  {c.task_id:28} {c.old} -> {_paint(c.new, c.new)}")
            print(f"{'':4}   {task.project} / {task.title}")
            print(f"{'':4}   {c.reason}")

        if args.dry_run:
            print("\n--dry-run: nothing written.")
            return 1

        stamp = pp.utcnow()
        with con:
            for c in changes:
                pp.record_status(
                    con, c.task_id, c.new,
                    evidence=pp.task_by_id(c.task_id).evidence,
                    note=f"sync: {c.reason}", when=stamp)
        print(f"\nrecorded {len(changes)} change(s) at {stamp}.")
        print("Re-render the pages:  python "
              "pipeline/scripts/update_project_plan.py render")
        return 0
    finally:
        con.close()


# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------
def cmd_render(args) -> int:
    from macro_core import report_projects as rp
    from macro_core import site

    written = rp.render_all(args.manifest,
                            projects=[args.project] if args.project else None)
    # A page a renderer has just rewritten has lost its chrome, so the site
    # is reassembled in the same breath.  Doing it here rather than asking a
    # person to remember a second command is the same reasoning that made
    # `set` and `render` one workflow in the first place: the step nobody
    # runs is the step that has to be automatic.
    #
    # `rp.DOCS_DIR` is read rather than assumed so a render into a temporary
    # tree (which the tests do) assembles that tree and never touches the
    # real docs/.
    site_written = site.build_site(manifest=args.manifest,
                                   docs_dir=rp.DOCS_DIR,
                                   repo_root=REPO_ROOT)
    # The chrome pass rewrites the project pages it just wrapped, so report
    # each path once, in the order it was first produced.
    seen = {p.resolve() for p in written}
    written = list(written) + [p for p in site_written
                               if p.resolve() not in seen]
    for path in written:
        # Repo-relative when it is under the repo; absolute otherwise, so a
        # render into a temporary tree reports rather than raises.
        try:
            shown = path.relative_to(REPO_ROOT)
        except ValueError:
            shown = path
        print(f"wrote {shown}")
    if args.project:
        print("\npartial render — WEB not recorded (it declares ALL the "
              "pages, and recording it now would claim the others are "
              "current too). Re-render everything to record it.")
        return 0
    return _record_web(args.manifest, written)


def _record_web(manifest: Path, written) -> int:
    """Record the WEB stage this command just satisfied.

    A tool whose whole thesis is that provenance must be recorded was not
    recording its own.  ``provenance`` declares a WEB stage whose outputs
    are exactly the pages ``render`` writes, so every render left
    ``check_pipeline_status.py status`` reporting WEB as STALE with one
    "written out of band" line per page this command had just produced — a
    permanent false alarm, generated by the one tool that should know
    better, and disclosed on none of the pages.
    """
    from macro_core import project_plan as pp
    from macro_core import provenance as pv

    stage = pv.STAGE_BY_KEY["WEB"]
    con = sqlite3.connect(str(manifest), timeout=300.0)
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        keys = sorted(set(stage.reads) | set(stage.writes))
        prints = pv.fingerprint_all(keys, con, REPO_ROOT)
        missing = [w for w in stage.writes if prints.get(w) == "MISSING"]
        if missing:
            # Never record a run that did not produce what it promised.
            print(f"\nWEB not recorded: declared output(s) still missing: "
                  f"{', '.join(sorted(missing))}")
            return 1
        run_utc = pp.utcnow()
        if run_utc in pv.recorded_run_times(con, "WEB"):
            print(f"\nWEB already recorded at {run_utc}; nothing written.")
            return 0
        # The code version recorded must be the one the DAG will compare
        # against, which for a stage declared hand_authored is its own
        # literal.  Writing the plan version here instead would make every
        # render report WEB as STALE on a code-version change — trading the
        # old false alarm for a new one.  The real generator version goes in
        # the note, which is what the note is for; when WEB is redeclared as
        # a generated stage with a version_file, this becomes a one-line
        # change and the note stops being the only record.
        pv.record_run(
            con, "WEB", run_utc,
            stage.code_version, _git_commit(),
            {k: prints[k] for k in stage.reads},
            {k: prints[k] for k in stage.writes},
            note=f"rendered {len(written)} page(s) by "
                 f"update_project_plan.py render "
                 f"({pp.PLAN_CODE_VERSION}, macro_core.report_projects)")
        print(f"\nrecorded stage WEB at {run_utc} ({len(written)} pages).")
        return 0
    finally:
        con.close()


def _git_commit() -> str:
    """Short commit of the working tree, with ``-dirty`` when it is.

    Same rule as ``check_pipeline_status.git_commit``: pages built from a
    dirty tree cannot be reproduced from the commit id alone, and the
    record has to say so.
    """
    def _run(*args: str) -> subprocess.CompletedProcess:
        return subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                              capture_output=True, text=True, timeout=20)
    try:
        out = _run("rev-parse", "--short", "HEAD")
        if out.returncode != 0:
            return "(unknown)"
        sha = out.stdout.strip()
        dirty = _run("status", "--porcelain")
        return f"{sha}-dirty" if dirty.stdout.strip() else sha
    except (OSError, subprocess.SubprocessError):   # pragma: no cover
        return "(unknown)"


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("show", help="print the plan and current status")
    s.add_argument("project", nargs="?",
                   help="one project key (or a task id with --history)")
    s.add_argument("--history", action="store_true",
                   help="replay the recorded status changes instead")
    s.add_argument("--live", action="store_true",
                   help="evaluate the DAG and the probes now and print the "
                        "synced view (about a minute)")
    s.set_defaults(func=cmd_show)

    s = sub.add_parser("blockers",
                       help="every blocked or gated task, with what holds "
                            "it — computed")
    s.add_argument("project", nargs="?")
    s.add_argument("--live", action="store_true",
                   help="also evaluate stage gates against the DAG")
    s.set_defaults(func=cmd_blockers)

    s = sub.add_parser("amendments",
                       help="markdown table of a group's ruled tasks")
    s.add_argument("project")
    s.set_defaults(func=cmd_amendments)

    s = sub.add_parser("titles",
                       help="check manuscript titles against the ledger")
    s.add_argument("--write", action="store_true",
                   help="re-emit the \\title line of each skeleton")
    s.set_defaults(func=cmd_titles)

    s = sub.add_parser("set", help="record a status change")
    s.add_argument("task_id")
    s.add_argument("status", choices=list(pp.ALL_STATUSES))
    s.add_argument("--evidence", default="",
                   help="repo-relative report page or product path")
    s.add_argument("--note", default="")
    s.add_argument("--ruling", default="",
                   help="the committee ruling behind a dropped/deferred "
                        "status, when the ledger's task carries none")
    s.set_defaults(func=cmd_set)

    s = sub.add_parser("sync", help="recompute statuses the DB has decided")
    s.add_argument("--dry-run", action="store_true")
    s.set_defaults(func=cmd_sync)

    s = sub.add_parser("render", help="regenerate the project pages + hub")
    s.add_argument("project", nargs="?", help="render only this project")
    s.set_defaults(func=cmd_render)

    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
