"""be_common — paths, fixed parameters and small helpers shared by the Be-star scripts.

WHY THIS MODULE EXISTS
----------------------
Every Be-star script (``be_step0.py``, ``be_injection.py``, ``be_external.py``,
``be_series.py`` …) reads the same three things: the shared manifest (READ-ONLY),
the novelty database that defines the sample (the chair's ruling on BE-N1-gate,
2026-10-04, takes the sample from it verbatim), and this project's own database
``BeStar_Grism/products/bestar.sqlite``.  Declaring the paths and the fixed
parameters once means no two scripts can disagree about what "the standards
epoch" or "a long night" is.

Nothing here touches the archive or writes anything; it is constants and
connection helpers only.
"""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROJECT = REPO / "BeStar_Grism"
MANIFEST = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
NOVELTY_DB = PROJECT / "notes" / "novelty" / "novelty.sqlite"
NOVELTY_TABLES = PROJECT / "notes" / "novelty" / "tables"
#: This project's database.  Gitignored (``products/``, ``*.sqlite``); every
#: table in it is regenerable by the scripts beside this file.
BE_DB = Path(os.environ.get("BE_DB", PROJECT / "products" / "bestar.sqlite"))
#: Tracked, script-emitted Markdown tables (the evidence the ledger links).
NOTES = PROJECT / "notes"
ARCHIVE_ROOT = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
PY = "/opt/miniconda3/envs/rlmt-checks/bin/python"

# ---------------------------------------------------------------------------
# Fixed parameters (declared before any series was looked at)
# ---------------------------------------------------------------------------
GRISM_FILTERS = ("hrg", "lrg", "HaGrism", "OGGrism")   # Step-0 whitelist
#: Detection claims exist only from this night on (SYNTHESIS §4, BE-S10).
STANDARDS_EPOCH = "2025-12-05"
#: Sample rule of the novelty gate (bess_novelty_check.py header).
MIN_NIGHTS_SEASON = 10
#: Short-tier admission (strategy §3.3): >= 3 nights each spanning > 2 h.
LONG_NIGHT_MIN = 120.0
MIN_LONG_NIGHTS = 3


def manifest_ro() -> sqlite3.Connection:
    """The shared manifest, strictly read-only (ground rule)."""
    con = sqlite3.connect(f"file:{MANIFEST}?mode=ro", uri=True, timeout=600)
    con.row_factory = sqlite3.Row
    return con


def novelty_ro() -> sqlite3.Connection:
    """The novelty database that defines the sample, read-only."""
    con = sqlite3.connect(f"file:{NOVELTY_DB}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def be_db() -> sqlite3.Connection:
    """This project's own database (read-write)."""
    BE_DB.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(BE_DB, timeout=600)
    con.row_factory = sqlite3.Row
    return con


def season_of(night: str) -> str:
    """Observing-year label for a night 'YYYY-MM-DD' (1 Aug -> 31 Jul)."""
    y, m = int(night[:4]), int(night[5:7])
    y0 = y if m >= 8 else y - 1
    return f"{y0}-{str(y0 + 1)[2:]}"


def md_table(df, floatfmt: int = 2) -> str:
    """A pandas frame as a GitHub Markdown table (no tabulate dependency)."""
    cols = list(df.columns)

    def cell(v):
        if v is None:
            return "—"
        if isinstance(v, float):
            if v != v:                       # NaN
                return "—"
            return f"{v:.{floatfmt}f}"
        return str(v)

    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for row in df.itertuples(index=False):
        lines.append("| " + " | ".join(cell(v) for v in row) + " |")
    return "\n".join(lines) + "\n"


def write_md(path: Path, title: str, body: str, script: str) -> None:
    """Write an emitted Markdown table with its do-not-edit header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"<!-- emitted by {script}; do not edit -->\n### {title}\n\n{body}",
                    encoding="utf-8")
