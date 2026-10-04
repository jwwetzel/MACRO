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
* :func:`read_hardware_cards` — the same tolerant card reader pointed at
  the cards that record the HARDWARE state of a frame (gain/offset
  setting, cooler, focuser, flip state, wheel map).  Added for the F-2
  header re-scrape; see the section at the bottom of this file.

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
    "HARDWARE_CARDS",
    "parse_hardware_cards",
    "read_hardware_cards",
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


# ===========================================================================
# Hardware-state cards (finding F-2, plan review 2026-10-03)
# ===========================================================================
# WHY THIS LIVES HERE.  The original catalog scan read headers through
# astropy and kept a fixed list of cards.  Two things went wrong with that,
# and both are the same story as the geometry artifact above:
#
# * the cards that say WHAT HARDWARE TOOK THE FRAME — gain and offset
#   setting, cooler set-point and power, the true focuser position, the
#   software flip state, the pier side, the filter-wheel slot and the wheel's
#   whole slot map — were never in the list, so the manifest could not tell a
#   re-seated camera from an untouched one (TE.F1, DE.F5);
# * ~20k files carry a malformed ``CONTINUE`` card *inside ``FWALLNAM``
#   itself*, so the one card that records the wheel map is exactly the card
#   astropy refuses to parse.
#
# So the re-scrape reads the raw 80-byte cards with the tolerant parser this
# module already owns.  It is also an order of magnitude cheaper than
# ``fits.open``: no HDU objects, no table machinery, one or two short reads
# per file — which matters when 330k files sit on one spinning disk.

#: The header cards the re-scrape collects, in the column order of the
#: ``hdr_rescrape`` table.  The first twelve are the list ruled by the
#: chair's synthesis (F-2); the last four are IDENTITY helpers the
#: mechanical-epoch table needs and that cost nothing to read in the same
#: pass: ``CAMNAME`` names the camera on the pyscope-native frames whose
#: ``INSTRUME`` is absent, ``XPIXSZ`` is the binned pixel pitch (a physical
#: camera property), ``FWNAME`` names the wheel hardware, and ``FOCUSPOS`` is
#: MaxIm's focuser card — kept beside pyscope's ``FOCPOS`` precisely because
#: TE.F8 showed it sticks, and a stuck card can only be demonstrated by
#: storing both.
HARDWARE_CARDS: tuple[str, ...] = (
    "INSTRUME", "GAIN", "OFFSET", "SET-TEMP", "CCD-TEMP", "COOLPOWR",
    "FOCPOS", "FLIPSTAT", "TELPIER", "FWPOS", "FWALLNAM", "SWCREATE",
    "CAMNAME", "XPIXSZ", "FWNAME", "FOCUSPOS",
)

#: A FITS header is written in 2880-byte blocks of 36 cards.
BLOCK_LEN = 2880

#: Upper bound on the header bytes read from one HDU.  The longest real
#: header in the archive is ~620 cards (50 kB); 400 blocks (1.15 MB) is far
#: beyond any header and stops a corrupt file with no ``END`` card from being
#: read to its last byte.
MAX_HEADER_BLOCKS = 400


def _card_string_value(field: str) -> tuple[str, bool]:
    """Extract a quoted string from one card's value field.

    Returns ``(text, continues)``; ``continues`` is True when the string
    ends in the long-string continuation marker ``&``.

    Two dialects must both survive:

    * the standard one, where an embedded quote is doubled —
      ``'(''g'', ''OGGrism'', &'``;
    * pyscope's native header (2026-06-28 onward), which writes the wheel
      map with RAW embedded quotes — ``'('g', 'lrg', 'r')' / Filter`` — an
      illegal card under the standard grammar, which would read it as the
      one-character string ``(``.

    The standard parse is tried first and accepted only when what follows
    the closing quote is empty or a comment.  Otherwise the card is the
    second dialect and the string is everything between the first quote and
    the last quote that precedes the comment separator.
    """
    body = field.lstrip()
    # --- standard grammar: scan to the first undoubled quote --------------
    i, out = 1, []
    while i < len(body):
        ch = body[i]
        if ch == "'":
            if i + 1 < len(body) and body[i + 1] == "'":
                out.append("'")
                i += 2
                continue
            break
        out.append(ch)
        i += 1
    rest = body[i + 1:].strip()
    if i < len(body) and (not rest or rest.startswith("/")):
        text = "".join(out).rstrip()
    else:
        # --- raw-quote dialect: last quote before the comment separator ---
        cut = body.rfind(" /")
        head = body if cut == -1 else body[:cut]
        end = head.rfind("'")
        text = (head[1:end] if end > 0 else head[1:]).rstrip()
    continues = text.endswith("&")
    return (text[:-1] if continues else text), continues


def parse_hardware_cards(block: bytes | str,
                         wanted: Iterable[str] = HARDWARE_CARDS) -> dict:
    """Pull the wanted cards out of a raw header block → ``{KEY: text}``.

    Every value comes back as STRIPPED TEXT, exactly as the header spells
    it; typing is the manifest's job, because one card means different
    things on different cameras (``GAIN`` is ``'4x'`` on the Andor iKon and
    ``100`` on the ASI) and coercing here would destroy that.

    Three states are kept apart, because the acceptance test for the
    re-scrape ("no nulls where the card exists") turns on the difference:

    * key ABSENT from the result — the header has no such card;
    * value ``''`` — the card exists and is blank (``OFFSET  =  / Image
      offset``; the post-monsoon ``FLIPSTAT= '        '``, whose blankness
      is itself the measurement);
    * anything else — the value.

    ``CONTINUE`` cards are FOLLOWED here (unlike :func:`parse_card_block`,
    which skips them): the wheel map is a long string and its tail lives in
    them.  The first occurrence of a keyword wins.

    Raises :class:`GeometryError` on a block that is not card-structured.
    """
    text = block.decode("ascii", "replace") if isinstance(block, bytes) \
        else block
    if not text or len(text) % CARD_LEN:
        raise GeometryError(
            f"header block is {len(text)} bytes, not a whole number of "
            f"{CARD_LEN}-byte cards")
    want = frozenset(wanted)
    out: dict[str, str] = {}
    open_key: Optional[str] = None       # a string still being continued
    for i in range(0, len(text), CARD_LEN):
        card = text[i:i + CARD_LEN]
        key = card[:8].strip()
        if key == "END":
            break
        if key == "CONTINUE":
            if open_key is not None:
                more, cont = _card_string_value(card[8:]) \
                    if "'" in card[8:] else ("", False)
                out[open_key] += more
                if not cont:
                    open_key = None
            continue
        open_key = None
        if key not in want or key in out or card[8:10] != "= ":
            continue
        field = card[10:]
        if field.lstrip().startswith("'"):
            value, cont = _card_string_value(field)
            out[key] = value
            if cont:
                open_key = key
        else:
            # Unquoted: number, logical, or nothing.  Cut the comment.
            out[key] = field.split("/", 1)[0].strip()
    # A continued string may carry trailing pad from its last segment.
    return {k: v.strip() for k, v in out.items()}


def _read_one_header(fh) -> bytes:
    """Read one HDU header (through its ``END`` card) from an open binary
    file positioned at a header start.  Returns ``b''`` at end of file."""
    chunks: list[bytes] = []
    for _ in range(MAX_HEADER_BLOCKS):
        blk = fh.read(BLOCK_LEN)
        if len(blk) < BLOCK_LEN:
            if not chunks and not blk:
                return b""
            raise GeometryError("file ends inside a header block")
        chunks.append(blk)
        # END is a card of its own: 'END' + 77 spaces at an 80-byte boundary.
        for j in range(0, BLOCK_LEN, CARD_LEN):
            if blk[j:j + 8] == b"END     ":
                return b"".join(chunks)
    raise GeometryError(f"no END card within {MAX_HEADER_BLOCKS} blocks")


def read_hardware_cards(path: str,
                        wanted: Iterable[str] = HARDWARE_CARDS) -> dict:
    """Read the wanted hardware cards of one FITS file from disk.

    Reads the primary header and — only when the primary holds no image
    (``NAXIS = 0``, the layout of every fpack ``.fz`` file, whose real
    header sits in the first extension) — the first extension header,
    which follows immediately because a dataless primary has no data
    blocks.  Nothing beyond those header blocks is read, and the file is
    opened read-only.

    ``.gz`` files are read through :mod:`gzip`; everything else is read
    raw.  Cards from the primary win over the extension's, matching the
    first-occurrence rule.

    Raises :class:`GeometryError` for a file with no parsable header, and
    lets ``OSError`` through for a file that cannot be opened — the caller
    records either in its ``error`` column.
    """
    import gzip
    opener = gzip.open if path.lower().endswith(".gz") else open
    with opener(path, "rb") as fh:
        primary = _read_one_header(fh)
        if not primary:
            raise GeometryError("empty file")
        cards = parse_hardware_cards(primary, wanted)
        head = parse_card_block(primary)
        if head.get("NAXIS") == 0:
            ext = _read_one_header(fh)
            if ext:
                for k, v in parse_hardware_cards(ext, wanted).items():
                    cards.setdefault(k, v)
    return cards
