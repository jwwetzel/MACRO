"""Image geometry from FITS headers — the one place that knows about
tile compression.

WHY THIS MODULE EXISTS (the S0e geometry artifact, 2026-08-18)
--------------------------------------------------------------
A tile-compressed FITS file (``fpack``/``funpack``, the ``.fts.fz`` files
that make up the RLMT raw archive) does NOT store the image as an image.
It stores it as a **binary table** whose rows hold compressed tiles.  That
substitution has a trap in it, and the trap cost this project 18,381
wrongly-excluded frames:

* the on-disk header of that BINTABLE is a *table* header.  Its ``NAXIS1``
  is the **width of one table row in bytes** and its ``NAXIS2`` is the
  **number of table rows**.  For the RLMT's 4800x3211 detector those come
  out as ``NAXIS1 = 8`` (one 1PB variable-length-array descriptor: two
  4-byte ints) and ``NAXIS2 = 3211`` (one row per tile — the tiles are
  full rows, ``ZTILE1=4800, ZTILE2=1``).
* the TRUE image dimensions live in ``ZNAXIS1`` / ``ZNAXIS2``, alongside
  ``ZIMAGE = T`` and ``ZCMPTYPE`` which mark the extension as a compressed
  image in the first place.

Read the table header naively and a full 4800x3211 field looks like an
"8x3211 pixel sub-frame strip".  That is precisely the phantom geometry
that reached the catalog: 19,980 rows, a phantom camera era, and an
astrometry gate (``astrom.is_window_geometry``) that threw the frames away
as too narrow to plate-solve.  They are full frames.  Nothing about them
was ever narrow.

``astropy`` translates ``Z*`` for you when it can build a ``CompImageHDU``
— so code that reads ``hdu.header`` is already safe.  The rows that went
wrong are the ones where astropy *could not* finish reading the header
(these files carry a malformed ``CONTINUE`` card in ``FWALLNAM`` that makes
``Header.update`` raise ``VerifyError``) and some fallback path read the
raw table header instead, taking ``NAXIS1``/``NAXIS2`` at face value.

So this module offers two things, both pure and both unit-tested:

* :func:`resolve_geometry` — given any header-like mapping, return the
  TRUE image geometry, preferring ``ZNAXIS*`` whenever the header says it
  is a compressed image, and never silently guessing.
* :func:`parse_card_block` — a raw 80-character FITS card parser for the
  rescue path, tolerant of the malformed ``CONTINUE`` cards astropy
  rejects, so a frame astropy cannot fully parse still yields its
  geometry instead of a phantom.

HOUSE RULE ENFORCED HERE: a header that cannot be understood raises
:class:`GeometryError`.  It never returns a plausible-looking number.  The
whole incident happened because a wrong answer was easier to produce than
an honest failure.
"""

from __future__ import annotations

import re
from typing import Iterable, Mapping, Optional

__all__ = [
    "GeometryError",
    "is_compressed_header",
    "resolve_geometry",
    "parse_card_block",
    "geometry_from_card_block",
    "CARD_LEN",
]

#: Every FITS header card is exactly 80 bytes.  Headers are padded to
#: 2880-byte blocks (36 cards).  Both numbers are fixed by the standard.
CARD_LEN = 80

#: Keywords that mark an extension as a *tile-compressed image* rather than
#: a genuine binary table.  ``ZIMAGE = T`` is the formal flag from the FITS
#: tiled-image convention; ``ZCMPTYPE`` (RICE_1, GZIP_1, PLIO_1, ...) names
#: the algorithm and is present in every real fpack output.  Either one is
#: enough to distrust NAXIS1/NAXIS2 — we do not require both, because a
#: truncated rescue-path read may only have reached one of them.
_ZIMAGE_KEY = "ZIMAGE"
_ZCMPTYPE_KEY = "ZCMPTYPE"

#: A card image looks like  KEYWORD = value / comment.  We only need the
#: keyword (cols 1-8) and the value field, and only for a handful of
#: numeric/logical keywords, so the grammar stays deliberately small.
_CARD_RE = re.compile(
    r"^(?P<key>[A-Z0-9_\-]{1,8})\s*=\s*(?P<val>.*?)(?:\s*/(?![^']*')[^/]*)?$"
)

#: Integer value field, e.g. "                4800".
_INT_RE = re.compile(r"^[+-]?\d+$")


class GeometryError(ValueError):
    """Raised when a header cannot be resolved to an honest geometry.

    Deliberately a hard error, not a ``None`` return: the S0e incident was
    caused by a fallback that produced a wrong-but-plausible number when it
    should have refused.  Callers that legitimately tolerate unknown
    geometry (the catalog scanner records the message in its ``error``
    column) must catch this explicitly and say so.
    """


def is_compressed_header(hdr: Mapping) -> bool:
    """True when ``hdr`` describes a tile-compressed image extension.

    ``hdr`` is any mapping with FITS keywords as keys — an
    ``astropy.io.fits.Header``, or the plain dict that
    :func:`parse_card_block` returns.

    The test is the presence of the compression markers, NOT the value of
    ``XTENSION``: a rescue-path parse may have skipped ``XTENSION``, and a
    genuine uncompressed image never carries ``ZIMAGE``/``ZCMPTYPE`` at
    all, so presence alone is both necessary and sufficient here.

    ``ZIMAGE`` is checked for truthiness because different readers hand it
    over as ``True``, ``'T'``, or ``1`` depending on how the card was
    parsed; ``ZCMPTYPE`` merely has to exist.
    """
    if _ZCMPTYPE_KEY in hdr:
        return True
    if _ZIMAGE_KEY in hdr:
        v = hdr[_ZIMAGE_KEY]
        # 'F' is the FITS logical false; treat it as "not compressed".
        return not (v is False or v == "F" or v == 0)
    return False


def _as_int(value, keyword: str) -> int:
    """Coerce one header value to int, or raise :class:`GeometryError`.

    Catalog values arrive as floats (SQLite REAL), astropy hands over ints,
    and a raw card parse hands over strings.  All three must land on the
    same integer, and anything else must fail loudly rather than default.
    """
    if value is None:
        raise GeometryError(f"{keyword} is missing")
    if isinstance(value, bool):
        # bool is a subclass of int; a logical here means a corrupt header.
        raise GeometryError(f"{keyword} is a logical, not a dimension")
    try:
        f = float(value)
    except (TypeError, ValueError):
        raise GeometryError(f"{keyword} is not numeric: {value!r}") from None
    if f != int(f):
        raise GeometryError(f"{keyword} is not integral: {value!r}")
    n = int(f)
    if n <= 0:
        raise GeometryError(f"{keyword} is not positive: {n}")
    return n


def resolve_geometry(hdr: Mapping) -> tuple[int, int]:
    """Return the TRUE ``(naxis1, naxis2)`` image geometry from ``hdr``.

    The whole point of this function, in one sentence: **when the header
    says the extension is a compressed image, the geometry is ZNAXIS1 /
    ZNAXIS2, and NAXIS1 / NAXIS2 are table bookkeeping that must be
    ignored.**

    Three cases, and no fourth:

    * compressed (``ZIMAGE``/``ZCMPTYPE`` present) → ``ZNAXIS1``,
      ``ZNAXIS2``.  If those are absent or unusable the header is
      self-contradictory and we raise — falling back to NAXIS here is
      exactly the bug this module exists to prevent.
    * uncompressed → ``NAXIS1``, ``NAXIS2``.
    * neither usable → :class:`GeometryError`.

    Raises :class:`GeometryError` rather than returning a sentinel.
    """
    if is_compressed_header(hdr):
        if "ZNAXIS1" not in hdr or "ZNAXIS2" not in hdr:
            raise GeometryError(
                "header is marked as a compressed image "
                f"({_ZIMAGE_KEY}/{_ZCMPTYPE_KEY} present) but carries no "
                "ZNAXIS1/ZNAXIS2 — refusing to fall back to the BINTABLE "
                "NAXIS values, which are row-bytes and row-count")
        return (_as_int(hdr["ZNAXIS1"], "ZNAXIS1"),
                _as_int(hdr["ZNAXIS2"], "ZNAXIS2"))
    if "NAXIS1" not in hdr or "NAXIS2" not in hdr:
        raise GeometryError("header carries neither ZNAXIS* nor NAXIS* "
                            "image dimensions")
    return (_as_int(hdr["NAXIS1"], "NAXIS1"),
            _as_int(hdr["NAXIS2"], "NAXIS2"))


def _parse_value(raw: str):
    """Turn one card's value field into a Python value.

    Only the shapes this module needs are recognised — quoted strings,
    logicals, integers — and anything else comes back as the stripped raw
    text.  Callers only ever ask for geometry keywords and the compression
    markers, so a loose tail here is harmless; :func:`_as_int` is the gate
    that decides whether a value is acceptable as a dimension.
    """
    raw = raw.strip()
    if not raw:
        return None
    if raw.startswith("'"):
        # Quoted string: take up to the closing quote, un-double '' escapes.
        end = raw.find("'", 1)
        while end != -1 and end + 1 < len(raw) and raw[end + 1] == "'":
            end = raw.find("'", end + 2)
        body = raw[1:end] if end != -1 else raw[1:]
        return body.replace("''", "'").strip()
    if raw in ("T", "F"):
        return raw == "T"
    if _INT_RE.match(raw):
        return int(raw)
    return raw


def parse_card_block(block: bytes | str) -> dict:
    """Parse a raw 80-character FITS card block into ``{keyword: value}``.

    This is the RESCUE path: it exists for files whose headers astropy
    refuses to finish (the RLMT archive's malformed ``CONTINUE`` cards,
    where a ``CONTINUE`` follows a non-string value and
    ``Header.update`` raises ``VerifyError``).  Rather than lose the frame
    — or, far worse, fall back to the BINTABLE ``NAXIS`` values and invent
    an 8-pixel sub-frame — we read the cards ourselves.

    Parsing rules, kept minimal on purpose:

    * the block is split into fixed 80-byte cards; a trailing partial card
      is a malformed header and raises :class:`GeometryError`;
    * ``END`` terminates the header;
    * ``CONTINUE``, ``COMMENT``, ``HISTORY`` and blank cards are SKIPPED —
      they never carry geometry, and skipping them is exactly what makes
      this parser survive where astropy stops;
    * the FIRST occurrence of a keyword wins, matching the FITS convention
      and keeping the result stable if a header repeats a card.

    Raises :class:`GeometryError` on a block that is not card-structured.
    """
    if isinstance(block, bytes):
        # FITS headers are ASCII by standard; be strict so that a binary
        # blob handed here fails loudly instead of yielding mojibake.
        try:
            text = block.decode("ascii")
        except UnicodeDecodeError as e:
            raise GeometryError(f"header block is not ASCII: {e}") from None
    else:
        text = block

    if not text:
        raise GeometryError("empty header block")
    if len(text) % CARD_LEN:
        raise GeometryError(
            f"header block is {len(text)} bytes, not a whole number of "
            f"{CARD_LEN}-byte cards — refusing to guess card boundaries")

    out: dict = {}
    for i in range(0, len(text), CARD_LEN):
        card = text[i:i + CARD_LEN]
        key = card[:8].strip()
        if key == "END":
            break
        if not key or key in ("CONTINUE", "COMMENT", "HISTORY"):
            continue
        m = _CARD_RE.match(card.strip())
        if not m or m.group("key") != key:
            # A card with no '=' (or a keyword we cannot read) is not fatal
            # on its own — skip it and let resolve_geometry decide whether
            # what survived is enough.
            continue
        if key not in out:                       # first occurrence wins
            out[key] = _parse_value(m.group("val"))
    if not out:
        raise GeometryError("header block contained no parsable cards")
    return out


def resolve_geometry_or_none(hdr: Mapping) -> tuple[Optional[int],
                                                    Optional[int]]:
    """Like :func:`resolve_geometry`, but returns ``(None, None)`` instead
    of raising when the header cannot be understood.

    This is NOT the silent fallback that caused the S0e artifact, and the
    difference matters.  The bug returned a *plausible wrong number*
    (8 x 3211) that flowed downstream as fact.  This returns ``None``,
    which every consumer already treats as "geometry unknown" — the
    astrometry gate refuses to promise a solvable field for it, and the
    timing code declines to compute a pixel scale from it.  Unknown is a
    safe answer; wrong is not.

    Use it in bulk scanners that must not die on one bad file; use
    :func:`resolve_geometry` anywhere a missing answer should stop the run.
    """
    try:
        return resolve_geometry(hdr)
    except GeometryError:
        return (None, None)


def geometry_from_card_block(block: bytes | str) -> tuple[int, int]:
    """Convenience: :func:`parse_card_block` then :func:`resolve_geometry`.

    The rescue path in one call, with the same loud-failure contract.
    """
    return resolve_geometry(parse_card_block(block))
