#!/usr/bin/env python
"""CV-S11 — the manuscript's figure set and its numbers, both by script.

WHAT THIS SCRIPT DOES, AND WHY EACH PIECE IS HERE
--------------------------------------------------
The CV project's last two tasks are a thirteen-figure set and a manuscript
draft.  Both are governed by the same law the web reports already obey: a
published number is the result of a query, never a keystroke.  This script
is how that law reaches a LaTeX document.

``figures``   Draws every figure in §7 of ``ANALYSIS_STRATEGY.md`` from the
              products database.  Each one is written TWICE from the same
              in-memory figure: a vector PDF into
              ``manuscripts/CV_TimeSeries/figures/`` for the paper, and a
              raster PNG into ``docs/CV_TimeSeries/figures/cv_paper/`` for
              the web page.  Writing both from one ``Figure`` object is not
              a convenience: it is the only arrangement in which the plot
              in the paper and the plot on the website cannot disagree.

``numbers``   Emits ``manuscripts/CV_TimeSeries/numbers.tex``, a file of
              ``\\newcommand`` definitions, one per value the prose is
              allowed to state, plus ``captions.tex``, one macro per
              figure caption, plus ``tables.tex``, the four measured
              tables.  ``main.tex`` inputs all three.  After this runs
              there is no number, no caption and no table cell in the
              manuscript source that a person typed.  Also available under
              the name ``report``, which is what the provenance graph
              calls a stage that renders a published artefact.

``manifest``  Prints what was built and, for the figures that could not be
              built as specified, the substitution and its reason.

``all``       figures, then numbers, then manifest.

THE SUBSTITUTIONS, STATED ONCE
-------------------------------
Three of the strategy's thirteen figures rest on observations that do not
exist, and one was written as a conditional whose other branch is the one
that happened.  Rather than drawing a plausible-looking version, each is
built as an honest substitute whose reason is stored in ``p5_figure``,
printed by ``manifest``, and carried into the caption that reaches the
paper.  See ``macro_phot.figures_cv``'s module docstring for the four.

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_cv_paper.py figures
    $P pipeline/scripts/run_cv_paper.py figures --only fig05 fig06
    $P pipeline/scripts/run_cv_paper.py numbers    # == 'report'
    $P pipeline/scripts/run_cv_paper.py manifest
    $P pipeline/scripts/run_cv_paper.py all

WHICH DATABASES THIS READS, AND WHICH IT WRITES
------------------------------------------------
It READS three, because this paper's numbers come from three: ``CV_DB``
(the photometry products), ``CH_DB`` (the characterisation products -- the
noise model, the timing budget, the check-star bias, the injection
contours) and ``MAN_DB`` (the frame manifest -- Table 1's detector
constants and per-mode ceilings).  All three are part of the release, and
§7 of the manuscript names all three; ``numbers_cv.resolve_databases``
records which one every macro came from so that description is emitted
rather than asserted.  It WRITES only ``CV_DB``.

TABLES WRITTEN (all inside products/phot/cv_timeseries.sqlite)
--------------------------------------------------------------
``p5_figure``   one row per figure: id, LaTeX label, title, full caption,
                the tables it was drawn from, both output paths, its width
                in inches, and -- when it is a substitute -- what was
                planned and why it could not be drawn.
``p5_number``   one row per macro: the macro name, the formatted LaTeX
                body, the unit, the source table, the RELEASED DATABASE
                that table is in, whether it is a measurement or an
                external constant, and the clause a referee needs about
                how it was derived.
``p5_meta``     build stamps and every constant this run used.

RESUMABILITY
------------
``figures`` deletes and rewrites only the rows for the figures it drew, so
``--only`` is safe.  ``numbers`` rewrites ``p5_number`` whole, because a
macro that disappeared from the collector must disappear from the file as
well or the paper keeps compiling against a value nothing produces.
Nothing here writes to any Phase-1, -2, -3 or -4 table.
"""

from __future__ import annotations

import argparse
import os
import sqlite3
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

import matplotlib                                      # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                        # noqa: E402

from macro_phot import figures_cv as fx                # noqa: E402
from macro_phot import numbers_cv as nx                # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
CV_DB = REPO_ROOT / "products" / "phot" / "cv_timeseries.sqlite"
CH_DB = REPO_ROOT / "products" / "phot" / "cv_characterization.sqlite"
MAN_DB = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

MANUSCRIPT_DIR = REPO_ROOT / "manuscripts" / "CV_TimeSeries"
PDF_DIR = MANUSCRIPT_DIR / "figures"
PNG_DIR = REPO_ROOT / "docs" / "CV_TimeSeries" / "figures" / "cv_paper"
NUMBERS_TEX = MANUSCRIPT_DIR / "numbers.tex"
CAPTIONS_TEX = MANUSCRIPT_DIR / "captions.tex"
TABLES_TEX = MANUSCRIPT_DIR / "tables.tex"

#: Stamped into ``p5_meta`` and read by the provenance graph.  Bump it when
#: a figure's ARITHMETIC or a macro's definition changes, not when a
#: comment or a colour does.
PAPER_CODE_VERSION = ("CV-S11 v1.1 (2026-08-20, manuscript figures + "
                      "numbers; p5_number carries db + kind)")

BUSY_TIMEOUT_MS = 300_000


# ===========================================================================
# Database plumbing
# ===========================================================================
def connect(path: Path, read_only: bool = False) -> sqlite3.Connection:
    """One connection, always with the long busy timeout.

    Other workflows write these archives; SQLite's own waiter is the only
    thing that keeps a concurrent Phase-3 re-run from turning this build
    into a 'database is locked' traceback halfway through figure nine.
    """
    if read_only:
        con = sqlite3.connect(f"file:{path}?mode=ro", uri=True, timeout=300)
    else:
        con = sqlite3.connect(str(path), timeout=300)
    con.row_factory = sqlite3.Row
    con.execute(f"PRAGMA busy_timeout={BUSY_TIMEOUT_MS}")
    return con


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def git_commit() -> str:
    """The commit this build came from, or a marker saying it is unknown."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO_ROOT, capture_output=True, text=True,
                             timeout=30)
        return out.stdout.strip() or "unknown"
    except Exception:                                       # noqa: BLE001
        return "unknown"


def ensure_tables(con: sqlite3.Connection) -> None:
    """Create this stage's three tables if they are not there yet."""
    con.executescript("""
        CREATE TABLE IF NOT EXISTS p5_figure (
            fig_id TEXT PRIMARY KEY,
            label TEXT, title TEXT, caption TEXT, tables_used TEXT,
            pdf_path TEXT, png_path TEXT, width_in REAL,
            substitute INTEGER, substitute_reason TEXT, note TEXT,
            built_utc TEXT
        );
        CREATE TABLE IF NOT EXISTS p5_number (
            macro TEXT PRIMARY KEY,
            key TEXT, value_tex TEXT, unit TEXT, source TEXT, note TEXT,
            db TEXT, kind TEXT
        );
        CREATE TABLE IF NOT EXISTS p5_meta (key TEXT PRIMARY KEY,
                                            value TEXT);
    """)
    # ``db`` and ``kind`` arrived after the first build wrote this table, and
    # CREATE TABLE IF NOT EXISTS will not add a column to a table that is
    # already there.  Without the migration an old products database keeps a
    # six-column p5_number and the insert below fails at build time.
    have = {r[1] for r in con.execute("PRAGMA table_info(p5_number)")}
    for col in ("db", "kind"):
        if col not in have:
            con.execute(f"ALTER TABLE p5_number ADD COLUMN {col} TEXT")
    con.commit()


def set_meta(con: sqlite3.Connection, pairs: dict) -> None:
    con.executemany("INSERT INTO p5_meta(key,value) VALUES(?,?) "
                    "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    [(k, str(v)) for k, v in pairs.items()])
    con.commit()


def record_stage(key: str) -> None:
    """Tell the provenance graph this stage ran.  Never fatal."""
    try:
        subprocess.run([sys.executable,
                        str(PIPELINE_ROOT / "scripts" /
                            "check_pipeline_status.py"), "record", key],
                       cwd=REPO_ROOT, capture_output=True, text=True,
                       timeout=180)
    except Exception as exc:                                # noqa: BLE001
        print(f"  ! could not record stage {key}: {exc}")


def write_atomic(path: Path, text: str) -> None:
    """Write through a temp file in the same directory, then rename.

    A half-written ``numbers.tex`` is worse than none: LaTeX would compile
    it and the paper would carry whichever macros happened to land.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


# ===========================================================================
# figures
# ===========================================================================
def cmd_figures(args) -> None:
    fx.apply_style()
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    PNG_DIR.mkdir(parents=True, exist_ok=True)

    cv = connect(CV_DB, read_only=True)
    ch = connect(CH_DB, read_only=True) if CH_DB.exists() else None
    man = connect(MAN_DB, read_only=True) if MAN_DB.exists() else None
    out = connect(CV_DB)
    ensure_tables(out)

    wanted = args.only or list(fx.FIGURE_IDS)
    unknown = [w for w in wanted if w not in fx.BUILDERS]
    if unknown:
        raise SystemExit(f"unknown figure id(s): {', '.join(unknown)}; "
                         f"known: {', '.join(fx.FIGURE_IDS)}")

    handles = {"cv": cv, "ch": ch, "man": man}
    built, skipped = 0, []
    for fig_id in wanted:
        entry = fx.BUILDERS[fig_id]
        missing = [n for n in entry["needs"] if handles[n] is None]
        if missing:
            skipped.append((fig_id, f"missing database(s): "
                                    f"{', '.join(missing)}"))
            print(f"  - {fig_id}: SKIPPED, {skipped[-1][1]}")
            continue
        cons = [handles[n] for n in entry["needs"]]
        try:
            fig, spec = entry["fn"](*cons)
        except Exception as exc:                            # noqa: BLE001
            skipped.append((fig_id, f"{type(exc).__name__}: {exc}"))
            print(f"  ! {fig_id}: FAILED, {skipped[-1][1]}")
            continue

        pdf = PDF_DIR / f"{spec.fig_id}_{spec.label.split(':')[-1]}.pdf"
        png = PNG_DIR / f"{spec.fig_id}_{spec.label.split(':')[-1]}.png"
        # Both from the SAME figure object, in the same call: the paper's
        # plot and the website's plot cannot drift apart.
        fig.savefig(pdf, format="pdf")
        fig.savefig(png, format="png", dpi=fx.PNG_DPI)
        plt.close(fig)

        out.execute("DELETE FROM p5_figure WHERE fig_id=?", (spec.fig_id,))
        out.execute("""INSERT INTO p5_figure(fig_id,label,title,caption,
                       tables_used,pdf_path,png_path,width_in,substitute,
                       substitute_reason,note,built_utc)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (spec.fig_id, spec.label, spec.title, spec.full_caption,
                     ", ".join(spec.tables),
                     str(pdf.relative_to(REPO_ROOT)),
                     str(png.relative_to(REPO_ROOT)), spec.width_in,
                     1 if spec.substitute else 0, spec.substitute_reason,
                     spec.note, utcnow()))
        out.commit()
        mark = " [SUBSTITUTE]" if spec.substitute else ""
        print(f"  + {spec.fig_id}  {spec.title}{mark}")
        print(f"      pdf {pdf.relative_to(REPO_ROOT)}")
        print(f"      png {png.relative_to(REPO_ROOT)}")
        built += 1

    set_meta(out, {"stage_figures": utcnow(),
                   "paper_code_version": PAPER_CODE_VERSION,
                   "git_commit": git_commit(),
                   "figures_built": built,
                   "png_dpi": fx.PNG_DPI,
                   "col_single_in": fx.COL_SINGLE,
                   "col_double_in": fx.COL_DOUBLE})
    print(f"\n{built} figure(s) written; {len(skipped)} not built.")
    for fid, why in skipped:
        print(f"  ! {fid}: {why}")
    record_stage("CV-S11")


# ===========================================================================
# numbers
# ===========================================================================
def cmd_numbers(args) -> None:
    cv = connect(CV_DB, read_only=True)
    ch = connect(CH_DB, read_only=True)
    man = connect(MAN_DB, read_only=True)
    out = connect(CV_DB)
    ensure_tables(out)

    numbers = nx.collect(cv, ch, man)
    stamp = (f"built {utcnow()} from {CV_DB.name}, {CH_DB.name} and "
             f"{MAN_DB.name} at commit {git_commit()} by "
             f"{PAPER_CODE_VERSION}")
    write_atomic(NUMBERS_TEX, nx.render_tex(numbers, stamp=stamp))

    out.execute("DELETE FROM p5_number")
    out.executemany("""INSERT INTO p5_number(macro,key,value_tex,unit,
                       source,note,db,kind) VALUES(?,?,?,?,?,?,?,?)""",
                    [(n.macro, n.key, n.body, n.unit, n.source, n.note,
                      n.db, n.kind) for n in numbers])
    out.commit()

    missing = [n for n in numbers if n.value is None]
    print(f"  + {NUMBERS_TEX.relative_to(REPO_ROOT)}: "
          f"{len(numbers)} macros, {len(missing)} unmeasured")
    for n in missing:
        print(f"      ! \\{n.macro} has no value in the database")

    # The captions the figure builders returned, as macros, so the paper's
    # \caption{} bodies are emitted by the same script that drew the panels.
    figs = out.execute("""SELECT fig_id,label,title,caption FROM p5_figure
                          ORDER BY fig_id""").fetchall()
    lines = [
        "%% captions.tex -- GENERATED FILE.  DO NOT EDIT.",
        "%% One macro per figure caption, emitted by the same script that "
        "drew the figure,",
        "%% so that a caption cannot describe a panel the code no longer "
        "produces.",
        f"%% {stamp}",
        "",
    ]
    for r in figs:
        macro = nx.tex_macro_name(r["fig_id"], prefix="Cap")
        lines.append(f"\\newcommand{{\\{macro}}}{{{r['caption']}}}")
    lines.append("")
    write_atomic(CAPTIONS_TEX, "\n".join(lines))
    print(f"  + {CAPTIONS_TEX.relative_to(REPO_ROOT)}: "
          f"{len(figs)} caption macros")

    write_atomic(TABLES_TEX, nx.render_tables(cv, man, stamp=stamp))
    print(f"  + {TABLES_TEX.relative_to(REPO_ROOT)}: measured tables")

    set_meta(out, {"stage_numbers": utcnow(),
                   "paper_code_version": PAPER_CODE_VERSION,
                   "git_commit": git_commit(),
                   "numbers_emitted": len(numbers),
                   "numbers_unmeasured": len(missing)})
    record_stage("R-CV-S11")


# ===========================================================================
# manifest
# ===========================================================================
def cmd_manifest(args) -> None:
    con = connect(CV_DB, read_only=True)
    if not con.execute("SELECT count(*) FROM sqlite_master WHERE type='table'"
                       " AND name='p5_figure'").fetchone()[0]:
        raise SystemExit("no p5_figure table: run 'figures' first")
    rows = con.execute("SELECT * FROM p5_figure ORDER BY fig_id").fetchall()
    print(f"{len(rows)} figure(s) in the products database\n")
    for r in rows:
        mark = "SUBSTITUTE" if r["substitute"] else "as specified"
        print(f"{r['fig_id']}  {r['title']}   [{mark}]")
        print(f"    label   \\ref{{{r['label']}}}")
        print(f"    tables  {r['tables_used']}")
        print(f"    pdf     {r['pdf_path']}")
        print(f"    png     {r['png_path']}")
        if r["substitute"]:
            print(f"    why     {r['substitute_reason']}")
        if r["note"]:
            print(f"    note    {r['note']}")
        print()
    n = con.execute("SELECT count(*) FROM p5_number").fetchone()[0] \
        if con.execute("SELECT count(*) FROM sqlite_master WHERE "
                       "type='table' AND name='p5_number'").fetchone()[0] \
        else 0
    print(f"{n} macro(s) in p5_number")


def cmd_all(args) -> None:
    cmd_figures(args)
    cmd_numbers(args)
    cmd_manifest(args)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("figures", help="draw every manuscript figure")
    p.add_argument("--only", nargs="*", default=None,
                   help="figure ids to (re)build, e.g. fig05 fig06")
    p.set_defaults(func=cmd_figures)

    p = sub.add_parser("numbers",
                       help="emit numbers.tex, captions.tex and tables.tex")
    p.set_defaults(func=cmd_numbers)

    # 'report' is the same action under the name the provenance graph uses
    # for every stage that renders a published artefact rather than
    # measuring something.  Emitting the manuscript's macros IS that stage
    # for this paper, and the graph's own test refuses a report stage whose
    # command does not look like one.
    p = sub.add_parser("report", help="alias for 'numbers'")
    p.set_defaults(func=cmd_numbers)

    p = sub.add_parser("manifest", help="print what was built")
    p.set_defaults(func=cmd_manifest)

    p = sub.add_parser("all", help="figures, numbers, manifest")
    p.add_argument("--only", nargs="*", default=None)
    p.set_defaults(func=cmd_all)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
