"""macro_tcrb.db — the T CrB project database and the read-only inputs.

One SQLite file, ``products/tcrb/tcrb.sqlite``, owned by this project.  The
manifest and the grism library are SHARED and are only ever opened
read-only here (``file:...?mode=ro``): GROUND_RULES forbid this project to
write any table it does not own, and a read-only URI makes that a property
of the connection rather than of the programmer's care.

Every stage writes its tables through :func:`write_table`, which replaces
the table wholesale (stages are idempotent: re-running one reproduces its
table exactly) and records the stage, the row count, the code version and
the UTC time in ``tcrb_meta``.
"""

from __future__ import annotations

import json
import sqlite3
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

from macro_tcrb import TCRB_CODE_VERSION

REPO = Path(__file__).resolve().parents[2]
TCRB_DB = REPO / "products" / "tcrb" / "tcrb.sqlite"
MANIFEST_DB = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
GRISM_DB = REPO / "products" / "grism" / "grism.sqlite"
ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
EXTERNAL = REPO / "TCrB_Monitoring" / "external_data"
NOVELTY = REPO / "TCrB_Monitoring" / "notes" / "novelty"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:                                  # pragma: no cover
        return "unknown"


def connect() -> sqlite3.Connection:
    """Read-write connection to the project database (created on demand)."""
    TCRB_DB.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(TCRB_DB)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("""CREATE TABLE IF NOT EXISTS tcrb_meta (
        stage TEXT PRIMARY KEY, tables TEXT, n_rows TEXT, code_version TEXT,
        git_commit TEXT, run_utc TEXT, extra TEXT)""")
    return con


def manifest_ro() -> sqlite3.Connection:
    """Read-only connection to the shared manifest."""
    return sqlite3.connect(f"file:{MANIFEST_DB}?mode=ro", uri=True)


def grism_ro() -> sqlite3.Connection:
    """Read-only connection to the shared grism library."""
    return sqlite3.connect(f"file:{GRISM_DB}?mode=ro", uri=True)


def write_table(con: sqlite3.Connection, name: str, columns: Sequence[str],
                rows: Iterable[Sequence], types: dict | None = None) -> int:
    """Replace table ``name`` with ``rows``; return the row count.

    ``types`` maps a column to its SQLite affinity (default: no affinity,
    which stores Python ints/floats/str faithfully).
    """
    types = types or {}
    rows = [tuple(r) for r in rows]
    con.execute(f'DROP TABLE IF EXISTS "{name}"')
    cols = ", ".join(f'"{c}" {types.get(c, "")}'.strip() for c in columns)
    con.execute(f'CREATE TABLE "{name}" ({cols})')
    if rows:
        ph = ", ".join("?" for _ in columns)
        con.executemany(f'INSERT INTO "{name}" VALUES ({ph})', rows)
    con.commit()
    return len(rows)


def record_stage(con: sqlite3.Connection, stage: str,
                 counts: dict[str, int], extra: dict | None = None) -> None:
    con.execute(
        "INSERT OR REPLACE INTO tcrb_meta VALUES (?,?,?,?,?,?,?)",
        (stage, ",".join(counts), json.dumps(counts), TCRB_CODE_VERSION,
         git_commit(), utc_now(), json.dumps(extra or {}, default=str)))
    con.commit()


def q(con: sqlite3.Connection, sql: str, *args) -> list[tuple]:
    return con.execute(sql, args).fetchall()


def q1(con: sqlite3.Connection, sql: str, *args):
    r = con.execute(sql, args).fetchone()
    return None if r is None else r[0]
