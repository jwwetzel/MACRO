#!/usr/bin/env python
"""Build a SQLite observation catalog from RLMT archive FITS headers.

Scans every .fz/.fts/.fit/.fits under the archive root, extracts key header
values (merged primary + first extension), and writes one row per file to
rlmt-catalog.sqlite. Corrupt/unreadable files land in the same table with a
non-null `error` column. Safe to re-run: it skips paths already cataloged.
"""
import os
import sys
import sqlite3
import warnings
from concurrent.futures import ProcessPoolExecutor

warnings.filterwarnings("ignore")

ROOT = "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive"
DB = "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite"
EXTS = (".fz", ".fts", ".fit", ".fits")
WORKERS = 16

STR_KEYS = [
    "DATE-OBS", "FILTER", "IMAGETYP", "OBJECT", "OBSERVER", "READOUTM",
    "CAMSN", "TELESCOP", "INSTRUME", "OBJCTRA", "OBJCTDEC", "RA", "DEC",
    "NOTES", "COMMENT", "ORIGIN", "SWCREATE",
]
NUM_KEYS = [
    "NAXIS1", "NAXIS2", "XBINNING", "YBINNING", "EXPTIME", "JD", "AIRMASS",
    "FWHM", "ZMAG", "MOONANGL", "MOONPHAS", "FOCUSPOS", "FOCUSTEM",
    "CAMTEMP", "CCD-TEMP", "EGAIN", "PRIORITY", "CRVAL1", "CRVAL2",
]
BOOL_KEYS = ["PLTSOLVD"]


def col(k):
    return k.lower().replace("-", "_")


def parse_sex(val, is_ra):
    """Sexagesimal 'H M S' / 'H:M:S' -> degrees (RA in hours -> deg)."""
    try:
        parts = str(val).replace(":", " ").split()
        d = abs(float(parts[0])) + float(parts[1]) / 60 + float(parts[2]) / 3600
        if str(parts[0]).strip().startswith("-"):
            d = -d
        return d * 15.0 if is_ra else d
    except Exception:
        return None


def scan_one(path):
    from astropy.io import fits
    rel = os.path.relpath(path, ROOT)
    row = {"path": rel, "tree": rel.split(os.sep)[0], "error": None,
           "size": None, "ra_deg": None, "dec_deg": None}
    for k in STR_KEYS + NUM_KEYS + BOOL_KEYS:
        row[col(k)] = None
    try:
        row["size"] = os.path.getsize(path)
        with fits.open(path, memmap=False, ignore_missing_simple=True) as h:
            hdr = fits.Header()
            for hdu in h[:2]:
                hdr.update(hdu.header)
        for k in STR_KEYS:
            if k in hdr:
                row[col(k)] = str(hdr[k]).strip()
        for k in NUM_KEYS:
            if k in hdr:
                try:
                    row[col(k)] = float(hdr[k])
                except (TypeError, ValueError):
                    pass
        for k in BOOL_KEYS:
            if k in hdr:
                row[col(k)] = int(bool(hdr[k]))
        # best RA/Dec in degrees: plate solution first, else pointing
        if row.get("crval1") is not None and row.get("crval2") is not None:
            row["ra_deg"], row["dec_deg"] = row["crval1"], row["crval2"]
        elif row.get("ra") and row.get("dec"):
            row["ra_deg"] = parse_sex(row["ra"], True)
            row["dec_deg"] = parse_sex(row["dec"], False)
        elif row.get("objctra") and row.get("objctdec"):
            row["ra_deg"] = parse_sex(row["objctra"], True)
            row["dec_deg"] = parse_sex(row["objctdec"], False)
    except Exception as e:
        row["error"] = f"{type(e).__name__}: {e}"[:300]
    return row


def main():
    cols = (["path", "tree", "size", "ra_deg", "dec_deg", "error"]
            + [col(k) for k in STR_KEYS + NUM_KEYS + BOOL_KEYS])
    db = sqlite3.connect(DB)
    db.execute("CREATE TABLE IF NOT EXISTS obs (%s, PRIMARY KEY(path))"
               % ", ".join(f'"{c}"' for c in cols))
    done = {r[0] for r in db.execute("SELECT path FROM obs")}
    print(f"already cataloged: {len(done)}", flush=True)

    todo = []
    for dirpath, _dirs, files in os.walk(ROOT):
        for f in files:
            if f.lower().endswith(EXTS):
                p = os.path.join(dirpath, f)
                if os.path.relpath(p, ROOT) not in done:
                    todo.append(p)
    print(f"files to scan: {len(todo)}", flush=True)

    # Name the columns explicitly: the live table has since grown extra
    # enrichment columns (fn_user/fn_target/fn_filter/target_best, added by
    # enrich_filenames.py via ALTER TABLE) that this scanner doesn't supply.
    # A bare VALUES list would need every column and fails; named columns
    # let the enrichment fields default to NULL until that script fills them.
    ins = "INSERT OR REPLACE INTO obs (%s) VALUES (%s)" % (
        ", ".join(f'"{c}"' for c in cols), ",".join("?" * len(cols)))
    batch, n = [], 0
    with ProcessPoolExecutor(WORKERS) as ex:
        for row in ex.map(scan_one, todo, chunksize=64):
            batch.append([row.get(c) for c in cols])
            n += 1
            if len(batch) >= 2000:
                db.executemany(ins, batch)
                db.commit()
                batch = []
            if n % 10000 == 0:
                print(f"scanned {n}/{len(todo)}", flush=True)
    if batch:
        db.executemany(ins, batch)
        db.commit()
    errs = db.execute("SELECT COUNT(*) FROM obs WHERE error IS NOT NULL").fetchone()[0]
    total = db.execute("SELECT COUNT(*) FROM obs").fetchone()[0]
    print(f"DONE: {total} rows, {errs} errors", flush=True)
    db.close()


if __name__ == "__main__":
    main()
