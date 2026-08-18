#!/usr/bin/env python
"""Parse pyscope-style filenames into fn_user / fn_target / fn_filter columns.

Pattern (token order varies):  <user>_<target...>_<filter>_<exp>s_<datetime>[_n].fts.fz
Also fills `target_best` = OBJECT header if present, else filename target.
"""
import os
import re
import sqlite3

DB = "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite"
EXP_RE = re.compile(r"^\d+(\.\d+)?s$")
DT_RE = re.compile(r"^20\d\d-\d\d-\d\d")
USER_RE = re.compile(r"^[a-z]{2,4}\d*$")

db = sqlite3.connect(DB)
filters = {r[0] for r in db.execute(
    "SELECT DISTINCT filter FROM obs WHERE filter IS NOT NULL AND filter != ''")}
filters |= {"halpha", "grism", "red", "green", "clear"}

for c in ("fn_user", "fn_target", "fn_filter", "target_best"):
    try:
        db.execute(f'ALTER TABLE obs ADD COLUMN "{c}"')
    except sqlite3.OperationalError:
        pass  # already added


def parse(basename):
    stem = re.sub(r"\.(fts|fits|fit)(\.fz)?$", "", basename, flags=re.I)
    tokens = stem.split("_")
    if len(tokens) < 3 or not USER_RE.match(tokens[0]):
        return None, None, None
    user = tokens[0]
    exp_i = next((i for i, t in enumerate(tokens) if EXP_RE.match(t)), None)
    dt_i = next((i for i, t in enumerate(tokens) if DT_RE.match(t)), None)
    end = min(x for x in (exp_i, dt_i, len(tokens)) if x is not None)
    filt = None
    if end - 1 > 1 and tokens[end - 1] in filters:
        filt = tokens[end - 1]
        end -= 1
    elif exp_i is not None and exp_i + 1 < len(tokens) and tokens[exp_i + 1] in filters:
        filt = tokens[exp_i + 1]
    target = " ".join(tokens[1:end]) or None
    return user, target, filt


rows = db.execute("SELECT path FROM obs WHERE error IS NULL").fetchall()
batch = []
for (path,) in rows:
    u, t, f = parse(os.path.basename(path))
    batch.append((u, t, f, path))
db.executemany(
    "UPDATE obs SET fn_user=?, fn_target=?, fn_filter=? WHERE path=?", batch)
db.execute("""UPDATE obs SET target_best =
    COALESCE(NULLIF(TRIM(object),''), fn_target)""")
db.commit()

n = db.execute("SELECT COUNT(*) FROM obs WHERE fn_target IS NOT NULL").fetchone()[0]
tb = db.execute("""SELECT COUNT(*) FROM obs WHERE target_best IS NOT NULL
    AND imagetyp LIKE 'Light%'""").fetchone()[0]
print(f"fn_target parsed: {n}; light frames with target_best: {tb}")
db.execute("CREATE INDEX IF NOT EXISTS idx_target ON obs(target_best)")
db.execute("CREATE INDEX IF NOT EXISTS idx_date ON obs(date_obs)")
db.commit()
