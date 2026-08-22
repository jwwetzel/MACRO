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
``show [project]``  the plan and current status, as text.  ``--history``
                    replays the recorded status changes instead.
``set <task-id> <status> [--evidence X] [--note Y]``
                    record a status change, stamped with the UTC time.
                    Append-only: nothing is overwritten, so "when did this
                    become true?" always has an answer.
``sync``            recompute the statuses the DATABASE has already decided:
                    a task marked ``done`` whose stage is no longer FRESH
                    flips to ``redo_needed``, with the reason recorded in
                    the note; one whose stage came back FRESH flips back.
                    ``--dry-run`` prints without writing.
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
def cmd_show(args) -> int:
    pp.validate()
    con = open_manifest(args.manifest, read_only=True)
    try:
        if args.history:
            rows = pp.read_history(con, args.project if args.project and
                                   args.project in
                                   {t.id for t in pp.all_tasks()} else None)
            if not rows:
                print("no status changes recorded yet.")
                return 0
            print(f"{'when':21} {'task':28} {'status':13} note")
            for task_id, status, evidence, note, when in rows:
                extra = note or evidence or ""
                print(f"{when:21} {task_id:28} {status:13} {extra}")
            return 0

        recorded = pp.read_statuses(con)
        keys = ([args.project] if args.project
                else [p.key for p in pp.PROJECTS])
        for key in keys:
            if key not in pp.PROJECT_BY_KEY:
                print(f"unknown project: {key!r}; known: "
                      f"{', '.join(pp.PROJECT_BY_KEY)}", file=sys.stderr)
                return 2
            project = pp.PROJECT_BY_KEY[key]
            statuses = pp.overlay_statuses(project.tasks, recorded)
            counts = pp.status_counts(project.tasks, statuses)
            done, total = pp.progress_fraction(counts)
            print()
            print("=" * 78)
            print(f"{project.title}   {done}/{total} tasks complete")
            summary = "  ".join(f"{pp.STATUS_LABEL[s]}={counts[s]}"
                                for s in pp.ALL_STATUSES if counts[s])
            print(f"  {summary}")
            print("=" * 78)
            for phase in project.phases:
                print(f"\n  {phase.name}")
                for t in phase.tasks:
                    s = statuses[t.id]
                    print(f"    {_paint(s, f'[{s:12}]')} {t.id:28} "
                          f"{t.title}")
                    if s == pp.BLOCKED and t.blocker:
                        print(f"{'':17}   ! {t.blocker}")
            nxt = pp.next_up(project.tasks, statuses, limit=3)
            if nxt:
                print("\n  NEXT UP:")
                for t in nxt:
                    print(f"    - {t.id}: {t.title}")
        return 0
    finally:
        con.close()


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

    con = open_manifest(args.manifest, read_only=False)
    try:
        with con:
            stamp = pp.record_status(
                con, args.task_id, args.status,
                evidence=args.evidence or task.evidence,
                note=args.note or "")
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
        tasks = pp.all_tasks()
        statuses = pp.overlay_statuses(tasks, recorded)
        changes = pp.derive_sync(tasks, statuses, states)

        summary = {}
        for key, f in freshness.items():
            summary[f.state] = summary.get(f.state, 0) + 1
        print("stage verdicts: " + "  ".join(
            f"{s}={summary[s]}" for s in pv.ALL_STATES if s in summary))

        if not changes:
            print("nothing to sync — every recorded status already agrees "
                  "with the database.")
            return 0

        print(f"\n{len(changes)} status change(s) the database has already "
              f"decided:")
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
    s.set_defaults(func=cmd_show)

    s = sub.add_parser("set", help="record a status change")
    s.add_argument("task_id")
    s.add_argument("status", choices=list(pp.ALL_STATUSES))
    s.add_argument("--evidence", default="",
                   help="repo-relative report page or product path")
    s.add_argument("--note", default="")
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
