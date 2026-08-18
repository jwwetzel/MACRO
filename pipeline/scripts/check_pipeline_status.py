#!/usr/bin/env python3
"""check_pipeline_status.py — is anything we already computed still true?

THE COMMAND TO RUN BEFORE TRUSTING ANY RESULT.

    python pipeline/scripts/check_pipeline_status.py status
    python pipeline/scripts/check_pipeline_status.py status --verbose
    python pipeline/scripts/check_pipeline_status.py backfill
    python pipeline/scripts/check_pipeline_status.py record S0c
    python pipeline/scripts/check_pipeline_status.py plan

Prints the stage DAG with a freshness verdict per stage, the specific
reason each stale stage is stale, and the ORDERED list of stages that must
re-run.  Every judgement comes from ``macro_core.provenance`` — this file
only walks the graph, formats, and (for ``record``/``backfill``) writes
rows into the manifest's ``stage_provenance`` table.

SUBCOMMANDS
-----------
``status``    fingerprint every declared resource, compare against the
              recorded provenance, print the DAG + verdicts + re-run plan.
``plan``      just the ordered re-run plan with the exact commands.
``backfill``  seed ``stage_provenance`` for stages that already ran, using
              each stage's OWN ``*_build_meta`` timestamp as the run time.
              Inputs whose producer was rebuilt AFTER that timestamp are
              recorded with the ``UNRECORDED`` sentinel, so the system
              starts out truthful (stale where it is stale) instead of
              declaring an unexamined pipeline healthy.
``record``    record the CURRENT state of one stage as a fresh run.  Run it
              immediately after a stage's build script finishes.

WRITE DISCIPLINE
----------------
Only ``stage_provenance`` is ever written, only by ``backfill``/``record``,
inside one transaction, with ``busy_timeout = 300000`` so a concurrent
solve or build waits rather than fails.  ``status`` and ``plan`` open the
manifest READ-ONLY: a status check can never perturb what it inspects.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sqlite3
import subprocess
import sys
from pathlib import Path

# pipeline/scripts/check_pipeline_status.py -> pipeline/ on sys.path, so the
# packages import the same way they do from every other build script.
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent
REPO_ROOT = PIPELINE_DIR.parent
sys.path.insert(0, str(PIPELINE_DIR))

from macro_core import provenance as pv                      # noqa: E402

DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

#: Terminal colours, suppressed when stdout is not a TTY (so the output
#: pipes into a file or a log cleanly).
_COLOR = {
    pv.FRESH: "\033[32m",            # green
    pv.STALE: "\033[33m",            # amber
    pv.STALE_UPSTREAM: "\033[90m",   # grey — waiting, not accused
    pv.NEVER_RUN: "\033[36m",        # cyan
    pv.OUTPUT_MISSING: "\033[31m",   # red
}
_RESET = "\033[0m"

#: Report stages whose renderer has no command-line entry point of its own.
#: ``render <stage>`` calls ``render_report`` directly for these, which is
#: what makes the re-run plan's report steps executable without adding a
#: ``main()`` to a stage module another workflow is currently editing.
#: ``arg`` says what the renderer takes: the manifest, or the S4 product.
RENDERERS: dict[str, tuple[str, str, str]] = {
    "R-S0":  ("macro_core.report_s0",  "render_report", "manifest"),
    "R-S0b": ("macro_core.report_s0b", "render_report", "manifest"),
    "R-S0c": ("macro_core.report_s0c", "render_report", "manifest"),
}


def _paint(state: str, text: str) -> str:
    """Colour a verdict when the terminal can show it."""
    if not sys.stdout.isatty():
        return text
    return f"{_COLOR.get(state, '')}{text}{_RESET}"


def utcnow() -> str:
    """ISO-8601 UTC stamp, seconds resolution."""
    return dt.datetime.now(dt.timezone.utc).replace(
        microsecond=0).isoformat().replace("+00:00", "Z")


def git_commit() -> str:
    """Short git commit of the working tree, '(unknown)' outside a repo.

    ``-dirty`` is appended when the tree has uncommitted changes, because a
    stage built from a dirty tree cannot be reproduced from the commit id
    alone and the record must say so.
    """
    try:
        out = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse",
                              "--short", "HEAD"], capture_output=True,
                             text=True, timeout=20)
        if out.returncode != 0:
            return "(unknown)"
        sha = out.stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO_ROOT), "status",
                                "--porcelain"], capture_output=True,
                               text=True, timeout=30)
        return sha + ("-dirty" if dirty.stdout.strip() else "")
    except Exception:                                  # git absent / timeout
        return "(unknown)"


def open_manifest(path: Path, read_only: bool) -> sqlite3.Connection:
    """Open the manifest.  READ-ONLY unless a write subcommand asked for
    otherwise; always with a five-minute busy timeout, because sibling
    stages (an S1 batch solve, an S3 build) hold their own connections."""
    uri = f"file:{path}?mode=ro" if read_only else f"file:{path}"
    con = sqlite3.connect(uri, uri=True, timeout=300.0)
    con.execute("PRAGMA busy_timeout = 300000")
    return con


# ---------------------------------------------------------------------------
# Reading each stage's own build metadata — the evidence backfill rests on
# ---------------------------------------------------------------------------
def stage_meta(con: sqlite3.Connection, stage: pv.Stage) -> dict[str, str]:
    """``{key: value}`` from a stage's ``*_build_meta`` table, or ``{}``.

    These tables are the stages' own self-reports: when they ran, which code
    version, which commit.  They are the only surviving record of the runs
    that happened before this module existed, so backfill uses them rather
    than inventing a timestamp.
    """
    if not stage.meta_table:
        return {}
    row = con.execute("SELECT 1 FROM sqlite_master WHERE type='table' "
                      "AND name=?", (stage.meta_table,)).fetchone()
    if row is None:
        return {}
    return {k: v for k, v in con.execute(
        f"SELECT key, value FROM {stage.meta_table}")}


#: Stages whose build metadata lives in their OWN product database rather
#: than in the manifest: (repo-relative db path, meta table).
PRODUCT_META: dict[str, tuple[str, str]] = {
    "S4": ("products/phot/anuma_vvpup_prototype.sqlite", "s4_build_meta"),
    "CV-S4": ("products/phot/cv_timeseries.sqlite", "cv_build_meta"),
}


def product_meta(stage_key: str) -> dict[str, str]:
    """``{key: value}`` from a product database's build-meta table, or {}.

    S4 and CV-S4 write their own sqlite files and keep their self-report
    there.  Opened READ-ONLY: reading a product must never be able to alter
    one.
    """
    entry = PRODUCT_META.get(stage_key)
    if entry is None:
        return {}
    rel, table = entry
    path = REPO_ROOT / rel
    if not path.exists():
        return {}
    side = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
    try:
        side.execute("PRAGMA busy_timeout = 300000")
        return {k: v for k, v in side.execute(
            f"SELECT key, value FROM {table}")}
    except sqlite3.Error:
        return {}
    finally:
        side.close()


def artifact_mtime(stage: pv.Stage) -> str:
    """For hand-authored / rendered file stages: the newest write time of
    the files they produce, used as the run time when no build_meta exists.

    A file's mtime is weaker evidence than a build_meta stamp — it is the
    last time somebody SAVED the file, not the last time its numbers were
    re-derived — and the printed report says so wherever it is used.
    """
    stamps = []
    for key in stage.writes:
        spec = pv.RESOURCES[key]
        if spec.kind != "file":
            continue
        p = REPO_ROOT / spec.name
        if p.exists():
            stamps.append(p.stat().st_mtime)
    if not stamps:
        return ""
    return dt.datetime.fromtimestamp(max(stamps), dt.timezone.utc).replace(
        microsecond=0).isoformat().replace("+00:00", "Z")


def _iso(stamp: str) -> str:
    """Normalize an ISO-8601 UTC stamp for lexicographic comparison.

    The stages write two dialects — ``...+00:00`` (datetime.isoformat) and
    ``...Z`` (hand-written) — and ``'Z' > '+'`` in ASCII, so an unnormalized
    string compare would report a Z-stamped run as later than a +00:00 run
    taken at the same instant.  Trailing sub-second digits are also dropped:
    they are precision, not ordering information, at this granularity.
    """
    s = stamp.strip().replace("+00:00", "").rstrip("Z")
    return s.split(".")[0]


def resource_mtime(spec: pv.ResourceSpec) -> str:
    """Last-modified stamp of a resource that lives in a FILE, or ''.

    Used for inputs no stage in the DAG produces — above all the external
    header-scan catalog, whose rewrite by the geometry rescue is exactly the
    event the backfill has to notice.  Table resources return '' (a table
    has no mtime; its producer stage's build_meta carries that evidence).
    """
    if spec.kind == "table":
        return ""
    name = spec.database or spec.name
    path = name if os.path.isabs(name) else str(REPO_ROOT / name)
    if not os.path.exists(path):
        return ""
    return dt.datetime.fromtimestamp(
        os.stat(path).st_mtime, dt.timezone.utc).replace(
            microsecond=0).isoformat().replace("+00:00", "Z")


def stage_run_time(con: sqlite3.Connection, stage: pv.Stage) -> tuple[str, str]:
    """``(run_utc, evidence)`` for a stage that already ran, or ``('','')``.

    Preference order, strongest evidence first:
      1. the stage's own build_meta timestamp in the manifest;
      2. for S4 / CV-S4, the same key inside their product database;
      3. the newest mtime of the files it writes (rendered/hand-authored).

    Three spellings of "when did this run" are accepted, because the stages
    genuinely use three (``built_utc``, ``last_run_utc``, and S2c's
    ``built_at``).  Guessing one and silently returning '' for the others
    would make a stage that HAS run look like one that never did.
    """
    meta = stage_meta(con, stage)
    source = stage.meta_table or ""
    if not meta:
        meta = product_meta(stage.key)
        if meta:
            source = PRODUCT_META[stage.key][1]
    for key in ("built_utc", "last_run_utc", "built_at"):
        if meta.get(key):
            return meta[key], f"{source or 'build_meta'}.{key}"
    mtime = artifact_mtime(stage)
    if mtime:
        return mtime, "newest mtime of the files it writes (weak evidence)"
    return "", ""


def stage_code_version_meta(con: sqlite3.Connection,
                            stage: pv.Stage) -> dict[str, str]:
    """The stage's build_meta, wherever it lives (manifest or product)."""
    meta = stage_meta(con, stage)
    return meta or product_meta(stage.key)


def stage_code_version(con: sqlite3.Connection, stage: pv.Stage) -> str:
    """The code version a stage RECORDED at its last run (not today's).

    CV-S4 stamps ``cv_code_version`` (its own) alongside ``s4_code_version``
    (the library it borrows); the stage's own version is the one that
    decides its freshness.
    """
    meta = stage_code_version_meta(con, stage)
    if stage.key == "CV-S4":
        return meta.get("cv_code_version", "")
    return meta.get("code_version", "")


def current_code_version(stage: pv.Stage) -> str:
    """Today's value of the stage's version constant, read live.

    Three sources, in this order:

    * hand-authored stages (OPS, STRAT, WEB) and the external S0e tool have
      no constant; they return their literal marker so the version test is
      a no-op for them;
    * a stage with ``version_file`` has its constant PARSED out of that
      file's source — see :func:`pv.read_version_constant` for why parsing
      beats importing when sibling workflows are editing the file;
    * everything else imports the constant from the package that owns it.
    """
    if stage.hand_authored:
        return stage.code_version
    if stage.version_file:
        path = REPO_ROOT / stage.version_file
        if not path.exists():
            return f"({stage.version_symbol}: source file missing)"
        got = pv.read_version_constant(
            path.read_text(encoding="utf-8", errors="replace"),
            stage.version_symbol)
        return got or f"({stage.version_symbol}: unreadable)"
    versions = pv._code_versions()
    # A report stage ("R-S0b") is versioned by the code of the stage whose
    # tables it renders ("S0b") — the renderer and the builder ship together,
    # so a builder version bump must also invalidate its page.
    lookup = stage.key[2:] if stage.key.startswith("R-") else stage.key
    return versions.get(lookup, stage.code_version)


# ---------------------------------------------------------------------------
# The freshness computation, shared by every subcommand
# ---------------------------------------------------------------------------
def evaluate(con: sqlite3.Connection) -> tuple[dict, dict, dict]:
    """Return ``(freshness, fingerprints, records)`` for the whole DAG.

    ``fingerprints`` maps every declared resource key to its current token —
    computed ONCE and shared, so a resource that appears as an input to four
    stages is hashed once, not four times.
    """
    all_keys = sorted({k for s in pv.STAGES for k in s.reads + s.writes})
    fingerprints = pv.fingerprint_all(all_keys, con, REPO_ROOT)
    records = pv.read_records(con)
    raw: dict[str, pv.Freshness] = {}
    for stage in pv.STAGES:
        inputs = {k: fingerprints[k] for k in stage.reads}
        outputs = {k: fingerprints[k] for k in stage.writes}
        raw[stage.key] = pv.is_stale(stage, records.get(stage.key), inputs,
                                     outputs, current_code_version(stage))
    return pv.propagate_staleness(raw, pv.STAGES), fingerprints, records


# ---------------------------------------------------------------------------
# Subcommand: status
# ---------------------------------------------------------------------------
def cmd_status(args) -> int:
    con = open_manifest(args.manifest, read_only=True)
    freshness, fps, records = evaluate(con)
    order = pv.topological_order(pv.STAGES)
    writer = {w: s.key for s in pv.STAGES for w in s.writes}

    print("=" * 78)
    print("MACRO pipeline status — " + utcnow())
    print(f"manifest : {args.manifest}")
    print(f"repo     : {REPO_ROOT}  @ {git_commit()}")
    print(f"digests  : {pv.PROVENANCE_CODE_VERSION}")
    print("=" * 78)

    counts: dict[str, int] = {}
    for key in order:
        stage = pv.STAGE_BY_KEY[key]
        f = freshness[key]
        counts[f.state] = counts.get(f.state, 0) + 1
        parents = sorted({writer[r] for r in stage.reads
                          if r in writer and writer[r] != key})
        dep = ("<- " + ", ".join(parents)) if parents else "<- (external catalog)"
        rec = records.get(key)
        when = rec.run_utc if rec else "never recorded"
        print()
        print(f"{_paint(f.state, f'[{f.state:14}]')} {key:6} {stage.title}")
        print(f"{'':17} {dep}")
        print(f"{'':17} last run: {when}   code: {current_code_version(stage)}")
        for reason in f.reasons:
            print(f"{'':17}   ! {reason}")
        if args.verbose:
            for w in stage.writes:
                print(f"{'':17}   out {w} = {fps[w]}")

    print()
    print("-" * 78)
    print("SUMMARY: " + "  ".join(
        f"{_paint(k, k)}={v}" for k, v in sorted(counts.items())))
    plan = pv.rerun_plan(freshness, pv.STAGES)
    print()
    print_plan(plan, freshness)
    con.close()
    # Exit status is machine-readable: 0 = everything fresh, 1 = work to do.
    return 0 if not plan else 1


def print_plan(plan, freshness) -> None:
    """The ordered re-run plan, with the exact command(s) per stage.

    A stage may need SEVERAL commands (S1's design/run/autopsy, S2's six
    sub-stages); each is printed on its own line, because a plan that
    silently omits two thirds of the work is worse than no plan.
    """
    if not plan:
        print("RE-RUN PLAN: nothing to do — every stage is FRESH.")
        return
    print(f"ORDERED RE-RUN PLAN ({len(plan)} stages, dependency order):")
    for i, key in enumerate(plan, 1):
        stage = pv.STAGE_BY_KEY[key]
        f = freshness[key]
        tag = "" if f.state != pv.STALE_UPSTREAM else \
            "   (waiting on an ancestor — its own inputs still match)"
        print(f"  {i:2}. {key:6} [{f.state}] {stage.title}{tag}")
        for line in stage.build_lines:
            print(f"      $ {line}")
        if not stage.hand_authored:
            print(f"      $ python pipeline/scripts/check_pipeline_status.py "
                  f"record {key}")
        else:
            print(f"      $ python pipeline/scripts/check_pipeline_status.py "
                  f"record {key} --note '<what you re-derived>'")


def cmd_plan(args) -> int:
    con = open_manifest(args.manifest, read_only=True)
    freshness, _, _ = evaluate(con)
    print_plan(pv.rerun_plan(freshness, pv.STAGES), freshness)
    con.close()
    return 0


# ---------------------------------------------------------------------------
# Subcommand: backfill
# ---------------------------------------------------------------------------
def cmd_backfill(args) -> int:
    """Seed stage_provenance from the evidence the repo already holds.

    For each stage that shows evidence of having run:

    * ``run_utc``  = its own build_meta timestamp (or artifact mtime);
    * ``outputs``  = TODAY's output fingerprints.  Honest: these are the
      bytes a reader would find right now, which is what any verdict about
      them has to be about;
    * ``inputs``   = today's input fingerprints, EXCEPT for inputs whose
      producing stage ran LATER than this stage did.  Those are recorded
      with the ``UNRECORDED`` sentinel and the reason, which makes the
      stage read STALE from the first status run — the whole point of
      backfilling rather than declaring victory.
    """
    con = open_manifest(args.manifest, read_only=False)
    pv.ensure_table(con)
    all_keys = sorted({k for s in pv.STAGES for k in s.reads + s.writes})
    fps = pv.fingerprint_all(all_keys, con, REPO_ROOT)
    writer = {w: s.key for s in pv.STAGES for w in s.writes}

    # First pass: when did each stage run?
    run_times: dict[str, tuple[str, str]] = {}
    for stage in pv.STAGES:
        run_times[stage.key] = stage_run_time(con, stage)

    n_written = 0
    for stage in pv.STAGES:
        run_utc, evidence = run_times[stage.key]
        produced_anything = any(fps[w] != "MISSING" for w in stage.writes)
        if not run_utc and not produced_anything:
            print(f"  skip   {stage.key:6} no evidence it ever ran")
            continue
        if not run_utc:
            run_utc = "(unknown)"
            evidence = "outputs exist but no timestamp survives"

        inputs: dict[str, str] = {}
        for r in stage.reads:
            src = writer.get(r)
            if src:
                src_time = run_times.get(src, ("", ""))[0]
                src_label = src
            else:
                # No stage in the DAG produces this input — it comes from
                # outside (the external header-scan catalog).  Its file
                # mtime is the best available "when did it last change".
                src_time = resource_mtime(pv.RESOURCES[r])
                src_label = "external source"
            # String comparison is valid here: every timestamp is ISO-8601
            # UTC, which sorts lexicographically in time order.  Timestamps
            # are normalized to a common suffix first so that a stored
            # '+00:00' and a rendered 'Z' compare on their digits, not on
            # their punctuation.
            if src_time and run_utc != "(unknown)" and \
                    _iso(src_time) > _iso(run_utc):
                inputs[r] = (f"{pv.UNRECORDED}:{src_label} changed "
                             f"{src_time}, after this stage ran {run_utc}")
            else:
                inputs[r] = fps[r]
        outputs = {w: fps[w] for w in stage.writes}
        note = f"backfilled {utcnow()}; run time evidence: {evidence}"
        pv.record_run(con, stage.key, run_utc,
                      stage_code_version(con, stage)
                      or current_code_version(stage),
                      git_commit(), inputs, outputs, note)
        n_written += 1
        flagged = sum(1 for v in inputs.values() if v.startswith(pv.UNRECORDED))
        gone = sum(1 for v in outputs.values() if v == "MISSING")
        print(f"  record {stage.key:6} run_utc={run_utc}  "
              f"inputs={len(inputs)} ({flagged} unrecorded)  "
              f"outputs={len(outputs)} ({gone} missing)")
    con.commit()
    con.close()
    print(f"\nbackfilled {n_written} stages into {pv.PROVENANCE_TABLE}.")
    print("Run `status` next — stages whose inputs moved will now say so.")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: record
# ---------------------------------------------------------------------------
def record_objection(stage: pv.Stage, run_utc: str, prior: list[str],
                     note: str) -> str:
    """Why this ``record`` should be refused, or ``''`` to allow it.  PURE.

    ``record`` writes the row that makes a stage read FRESH, and a FRESH
    stage is a publish gate held open.  Before these checks it was an
    unguarded rubber stamp: typing ``record OPS`` while changing nothing
    turned OPS fresh, and once the last stage was stamped ``status`` exited
    0.  For OPS — the stage this audit has just shown to be materially
    wrong — that is precisely the wrong default.

    Three objections, each answering a different way the stamp can lie:

    1. a HAND-AUTHORED stage recorded with no note.  Nothing in the
       filesystem can distinguish "I re-derived these numbers" from "I
       opened the file", so a person has to say which, in writing;
    2. a run time that ALREADY EXISTS in the table.  The primary key is
       (stage, run_utc) and the write is INSERT OR REPLACE, so this would
       overwrite the earlier run rather than append one;
    3. a run time that does not ADVANCE past the newest recorded run.  That
       is what a re-run looks like when its build_meta was not updated: the
       stage would be stamped fresh under a timestamp older than work
       already recorded — a stamp preceding the thing it attests.

    Each is overridable with ``--force``, which is recorded in the note, so
    the override is itself evidence rather than a silent exception.
    """
    if stage.hand_authored and not note.strip():
        return (f"{stage.key} is hand-authored: nothing on disk can show "
                f"that its numbers were re-derived rather than merely "
                f"re-saved.  Re-run with --note '<what you re-derived>' "
                f"(or --force to record without one).")
    if run_utc in prior:
        return (f"a run at {run_utc} is already recorded for {stage.key}; "
                f"recording again would REPLACE it and lose that run's "
                f"fingerprints.  Pass --run-utc with the new run's time, or "
                f"--force to overwrite deliberately.")
    if prior and _iso(run_utc) < _iso(max(prior)):
        return (f"run time {run_utc} is older than the newest recorded run "
                f"{max(prior)}.  This is what a re-run whose build_meta was "
                f"not updated looks like; the stamp would precede the work "
                f"it attests.  Pass --run-utc explicitly, or --force.")
    return ""


def cmd_record(args) -> int:
    """Record one stage's CURRENT state as a completed run.

    Run this straight after the stage's build script: the inputs it saw are
    the inputs on disk right now, and the outputs it produced are the ones
    in the database right now.  See :func:`record_objection` for what this
    refuses to record, and why.
    """
    stage = pv.STAGE_BY_KEY.get(args.stage)
    if stage is None:
        print(f"unknown stage {args.stage!r}; known: "
              f"{', '.join(pv.STAGE_BY_KEY)}", file=sys.stderr)
        return 2
    con = open_manifest(args.manifest, read_only=False)
    run_utc = args.run_utc or stage_run_time(con, stage)[0] or utcnow()
    prior = pv.recorded_run_times(con, stage.key)
    objection = record_objection(stage, run_utc, prior, args.note or "")
    if objection and not args.force:
        con.close()
        print(f"REFUSED: {objection}", file=sys.stderr)
        return 3
    inputs = pv.fingerprint_all(stage.reads, con, REPO_ROOT)
    outputs = pv.fingerprint_all(stage.writes, con, REPO_ROOT)
    note = args.note or f"recorded by check_pipeline_status at {utcnow()}"
    if objection and args.force:
        # The override is kept IN the record.  A forced stamp that looks
        # like an ordinary one would be the worst of both worlds.
        note = f"{note} [--force overrode: {objection}]"
    pv.record_run(con, stage.key, run_utc,
                  stage_code_version(con, stage) or current_code_version(stage),
                  git_commit(), inputs, outputs, note=note)
    con.close()
    print(f"recorded {stage.key} @ {run_utc}: "
          f"{len(inputs)} inputs, {len(outputs)} outputs")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: render
# ---------------------------------------------------------------------------
def cmd_render(args) -> int:
    """Re-render one report stage's published page.

    Exists because three of the seven report modules have no CLI: their
    entry point is ``render_report(path)`` and nothing calls it except the
    build script that also rebuilds the whole stage.  The re-run plan needs
    a command that renders the page and ONLY the page, so it lives here —
    inside this module's own files — rather than as a ``main()`` bolted on
    to stage code a sibling workflow is editing.

    ``--check`` resolves the renderer and returns without rendering, which
    is what the regression test uses to prove every mapping still points at
    a function that exists.
    """
    entry = RENDERERS.get(args.stage)
    if entry is None:
        known = ", ".join(sorted(RENDERERS))
        print(f"no renderer registered for {args.stage!r}; known: {known}.\n"
              f"Other report stages have their own CLI — see `plan`.",
              file=sys.stderr)
        return 2
    module_name, symbol, arg_kind = entry
    import importlib
    module = importlib.import_module(module_name)
    fn = getattr(module, symbol, None)
    if fn is None:
        print(f"{module_name} has no {symbol}()", file=sys.stderr)
        return 2
    target = args.manifest if arg_kind == "manifest" else \
        REPO_ROOT / "products" / "phot" / "anuma_vvpup_prototype.sqlite"
    if args.check:
        print(f"{args.stage}: {module_name}.{symbol} resolves; "
              f"would render from {target}")
        return 0
    print(f"wrote {fn(target)}")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: report
# ---------------------------------------------------------------------------
def cmd_report(args) -> int:
    """Render docs/pipeline/pipeline_status.html from the database."""
    from macro_core import report_provenance
    con = open_manifest(args.manifest, read_only=True)
    freshness, fps, records = evaluate(con)
    out = report_provenance.render(con, REPO_ROOT, freshness, fps, records,
                                   manifest_path=args.manifest)
    con.close()
    print(f"wrote {out}")
    return 0


# ---------------------------------------------------------------------------
# Subcommand: resources (the fingerprint contract, printed)
# ---------------------------------------------------------------------------
def cmd_resources(args) -> int:
    """Print every resource and the justification for its column choice —
    so a reader can audit the fingerprint design without reading code."""
    for key in sorted(pv.RESOURCES):
        spec = pv.RESOURCES[key]
        print(f"\n{key}  [{spec.kind}]")
        if spec.columns:
            print(f"  columns : {', '.join(spec.columns)}")
            print(f"  order   : {spec.order_by}")
        print(f"  why     : {spec.why}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                    help="manifest database (default: %(default)s)")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("status", help="DAG + freshness + re-run plan")
    p.add_argument("--verbose", action="store_true",
                   help="also print each output's fingerprint")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("plan", help="just the ordered re-run plan")
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("backfill", help="seed provenance for past runs")
    p.set_defaults(func=cmd_backfill)

    p = sub.add_parser("record", help="record one stage's current state")
    p.add_argument("stage")
    p.add_argument("--run-utc", default="")
    p.add_argument("--note", default="",
                   help="what this run re-derived; REQUIRED for the "
                        "hand-authored stages (OPS, STRAT, WEB)")
    p.add_argument("--force", action="store_true",
                   help="record despite an objection (the objection is "
                        "written into the record)")
    p.set_defaults(func=cmd_record)

    p = sub.add_parser("render", help="re-render one report stage's page")
    p.add_argument("stage")
    p.add_argument("--check", action="store_true",
                   help="resolve the renderer and exit without rendering")
    p.set_defaults(func=cmd_render)

    p = sub.add_parser("report", help="render the status page")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("resources", help="print the fingerprint contract")
    p.set_defaults(func=cmd_resources)

    args = ap.parse_args(argv)
    if not getattr(args, "func", None):
        args = ap.parse_args((argv or []) + ["status"])
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
