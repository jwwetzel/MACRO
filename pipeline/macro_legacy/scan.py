"""L0 header scan — how one legacy FITS file becomes one database row.

This is the only module of ``macro_legacy`` that opens FITS files, and it
opens HEADERS only: no pixel is decompressed (astropy loads HDUs lazily, and
only ``.header`` is ever touched), so the scan is cheap and cannot be
affected by damaged data units.

WHAT IS READ, AND WHY THESE CARDS
---------------------------------
The card lists below are the union of what the committee asked for
(detector engineer: INSTRUME, camera serial, READOUTM, GAIN/EGAIN, SET-TEMP,
CCD-TEMP, binning, geometry, SWCREATE; telescope engineer: telescope, focal
length, flip state; referee: every time card, for the convention audit) and
what the first look at eight sample headers showed is actually written by
the two acquisition systems in the archive:

* **Talon** (2015, the Rigel): no IMAGETYP, no SWCREATE, no JD-HELIO;
  binning in XFACTOR/YFACTOR, temperature in CAMTEMP, a fractional-second
  DATE-OBS and a JD card.
* **MaxIm DL** (late 2015 onward, the 0.5 m): SBIG-extension cards —
  IMAGETYP, XBINNING, SET-TEMP/CCD-TEMP, FOCALLEN, SWCREATE, plus the
  scheduler's own XFACTOR, CAMTEMP, CALSTAT and CALSTART.

Because a census that silently lacks a card cannot say whether the card was
absent or merely not looked for, every file also records the SIGNATURE of
its full keyword set (``keyset_sig``); the distinct keyword sets are stored
once in the ``keysets`` table.  "Does any legacy header carry GAIN?" is then
a query, not a re-scan.

TILE COMPRESSION
----------------
Every scanned file is an fpack'd ``.fz``: an empty primary HDU followed by a
compressed-image extension.  The cards of both HDUs are merged (extension
wins), and the true image geometry is resolved by
``macro_core.fitsgeom.resolve_geometry`` — never read from NAXIS1/NAXIS2 of
a BINTABLE (the S0e lesson).
"""

from __future__ import annotations

import hashlib
import os
import warnings
from typing import Mapping, Optional

from macro_core import fitsgeom

# ---------------------------------------------------------------------------
# The cards.  Column name = keyword lower-cased with '-' -> '_'.
# ---------------------------------------------------------------------------
#: String-valued cards, stored stripped.
STR_KEYS = (
    # --- identity of the instrument and the software that wrote the file
    "TELESCOP", "INSTRUME", "DETECTOR", "CAMERA", "CAMSN", "SERIALNO",
    "READOUTM", "SWCREATE", "SWMODIFY", "SWSERIAL", "SWOWNER", "ORIGIN",
    "FLIPSTAT", "PIERSIDE", "TELPIER", "ROWORDER", "BAYERPAT", "HCOMSTAT",
    "FWALLNAM",
    # --- what the frame is
    "IMAGETYP", "FILTER", "OBJECT", "OBSERVER", "NOTES",
    "CALSTAT", "CALSTART",
    # --- every time card either acquisition system is known to write
    "DATE-OBS", "TIME-OBS", "UT", "DATE", "DATE-AVG", "DATE-END", "TIMESYS",
    "LST", "HA",
    # --- pointing (sexagesimal strings in both systems)
    "RA", "DEC", "OBJRA", "OBJDEC", "OBJCTRA", "OBJCTDEC",
)
#: Numeric cards, stored as REAL (unparseable values become NULL).
NUM_KEYS = (
    # --- detector configuration
    "XBINNING", "YBINNING", "XFACTOR", "YFACTOR", "OFFSET1", "OFFSET2",
    "XORGSUBF", "YORGSUBF", "XPIXSZ", "YPIXSZ", "GAIN", "EGAIN", "OFFSET",
    "RDNOISE", "PEDESTAL", "BZERO", "BSCALE",
    "SET-TEMP", "CCD-TEMP", "CAMTEMP", "COOLPOWR",
    # --- optics and mechanics
    "FOCALLEN", "APTDIA", "FOCUSPOS", "FOCPOS", "FOCUSTEM", "FWPOS",
    # --- exposure and time
    "EXPTIME", "EXPOSURE", "DARKTIME", "JD", "JD-OBS", "JD-HELIO", "HJD",
    "BJD", "MJD-OBS", "MJD",
    # --- pointing / astrometry / image quality as written at the telescope
    "AIRMASS", "CRVAL1", "CRVAL2", "CROTA2", "CDELT1", "CDELT2",
    "FWHMH", "FWHMV", "MOONANGL", "MOONPHAS", "PRIORITY",
)
#: Maximum stored length of a free-text value (COMMENT, NOTES, error).
TEXT_LIMIT = 300


def col(keyword: str) -> str:
    """Database column name for a FITS keyword (``'CCD-TEMP' -> 'ccd_temp'``)."""
    return keyword.lower().replace("-", "_")


#: Every column of the ``scan`` table, in order.  File facts first, then the
#: cards, then the derived per-file facts.
SCAN_COLUMNS = (
    ["path", "size", "mtime", "error", "n_hdu", "n_bad_cards"]
    + [col(k) for k in STR_KEYS]
    + [col(k) for k in NUM_KEYS]
    + ["naxis1", "naxis2", "bitpix", "comment", "keyset_sig", "keyset"]
)


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested without any FITS file)
# ---------------------------------------------------------------------------
def keyset_signature(keywords) -> tuple[str, str]:
    """Return ``(signature, canonical keyword list)`` for a header.

    The signature is the first 12 hex digits of the SHA-1 of the sorted,
    de-duplicated keyword names joined by commas.  COMMENT/HISTORY/blank
    cards are excluded (they repeat and say nothing about the writer), as
    are the tile-compression bookkeeping cards (``Z*``, ``TFORM*`` …) and
    checksums, which describe our fpack run, not the camera software.
    """
    keep = sorted({k for k in keywords if k and not _is_bookkeeping(k)})
    text = ",".join(keep)
    return hashlib.sha1(text.encode()).hexdigest()[:12], text


_BOOKKEEPING_EXACT = frozenset({
    "COMMENT", "HISTORY", "SIMPLE", "EXTEND", "XTENSION", "PCOUNT",
    "GCOUNT", "TFIELDS", "CHECKSUM", "DATASUM", "EXTNAME", "END",
    "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
})
_BOOKKEEPING_PREFIX = ("Z", "TTYPE", "TFORM", "TUNIT")


def _is_bookkeeping(keyword: str) -> bool:
    """True for cards that describe the container, not the observation."""
    return (keyword in _BOOKKEEPING_EXACT
            or keyword.startswith(_BOOKKEEPING_PREFIX))


def row_from_header(hdr: Mapping, comments: list[str] | None = None) -> dict:
    """Turn one merged header (any mapping) into the card part of a row.

    Pure: takes a mapping, returns a dict keyed by column name.  Strings are
    stripped (an all-blank string is kept as ``""`` — "the card exists and
    is empty" is information the camera-identity census needs, distinct from
    NULL = "no such card").  Numbers go through ``float``; a logical or an
    unparseable value becomes NULL rather than 0.

    Geometry and bit depth come from the compression-aware resolver.  A
    header whose geometry cannot be resolved yields NULL geometry and an
    ``error`` entry in the returned dict; the caller decides what to do.
    """
    row: dict = {}
    for k in STR_KEYS:
        row[col(k)] = str(hdr[k]).strip()[:TEXT_LIMIT] if k in hdr else None
    for k in NUM_KEYS:
        row[col(k)] = _as_float(hdr[k]) if k in hdr else None
    try:
        row["naxis1"], row["naxis2"] = fitsgeom.resolve_geometry(hdr)
    except fitsgeom.GeometryError as exc:
        row["naxis1"] = row["naxis2"] = None
        row["error"] = f"GeometryError: {exc}"[:TEXT_LIMIT]
    # ZBITPIX is the image's bit depth in a compressed extension; BITPIX
    # there describes the binary table (always 8).
    bitpix = hdr.get("ZBITPIX", hdr.get("BITPIX")) \
        if fitsgeom.is_compressed_header(hdr) else hdr.get("BITPIX")
    row["bitpix"] = _as_float(bitpix)
    text = " | ".join(c.strip() for c in (comments or []) if c and c.strip())
    row["comment"] = text[:TEXT_LIMIT] or None
    row["keyset_sig"], row["keyset"] = keyset_signature(hdr.keys())
    return row


def _as_float(value) -> Optional[float]:
    """``float(value)`` or None; logicals are refused (True is not 1.0)."""
    if value is None or isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# File access
# ---------------------------------------------------------------------------
#: FITS-standard free-text comments fpack and cfitsio add to every file; not
#: information about the observation, so not stored.
_BOILERPLATE = ("FITS (Flexible Image Transport System)",
                "and Astrophysics', volume 376")


def scan_one(args: tuple[str, str]) -> dict:
    """Scan one file.  ``args = (archive_root, relative_path)``.

    Never raises: any failure (file vanished under a concurrent job,
    truncated header, unreadable card) is recorded in the row's ``error``
    column, so files-on-disk always equals rows.  A single malformed card
    costs that card, not the frame (the CONTINUE-card lesson of the RLMT
    catalog scanner): cards are copied one by one inside a guard.
    """
    from astropy.io import fits            # imported in the worker process
    warnings.filterwarnings("ignore")
    root, rel = args
    path = os.path.join(root, rel)
    row = dict.fromkeys(SCAN_COLUMNS)
    row["path"] = rel
    try:
        st = os.stat(path)
        row["size"], row["mtime"] = st.st_size, st.st_mtime
        merged: dict = {}
        comments: list[str] = []
        n_bad = 0
        with fits.open(path, memmap=False, lazy_load_hdus=True,
                       ignore_missing_simple=True) as hdul:
            # Primary + first extension only: that is the whole frame in an
            # fpack'd file, and touching further HDUs would force reads.
            n_hdu = 0
            for i in (0, 1):
                try:
                    header = hdul[i].header
                except IndexError:
                    break
                n_hdu += 1
                for card in header.cards:
                    try:
                        key, value = card.keyword, card.value
                    except Exception:
                        n_bad += 1
                        continue
                    if key == "COMMENT":
                        text = str(value)
                        if not any(b in text for b in _BOILERPLATE):
                            comments.append(text)
                    elif key not in ("HISTORY", ""):
                        merged[key] = value
        row.update(row_from_header(merged, comments))
        row["n_hdu"], row["n_bad_cards"] = n_hdu, n_bad
    except Exception as exc:               # noqa: BLE001 — recorded, not hidden
        row["error"] = f"{type(exc).__name__}: {exc}"[:TEXT_LIMIT]
    return row
