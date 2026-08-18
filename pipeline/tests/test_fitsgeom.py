"""Unit tests for macro_core.fitsgeom — the tile-compression geometry trap.

Every test builds a SYNTHETIC 80-character card block, so the suite never
touches the archive and a student can read the exact bytes under test.

The three mandated cases:
  * compressed   → ZNAXIS wins, the BINTABLE NAXIS values are ignored
  * uncompressed → NAXIS wins
  * malformed    → fails loudly (GeometryError), never silently
"""

import pytest

from macro_core import fitsgeom as fg


# ---------------------------------------------------------------------------
# Synthetic card-block helpers
# ---------------------------------------------------------------------------
def card(keyword: str, value: str = None, comment: str = None) -> str:
    """Build one exactly-80-character FITS card image."""
    if value is None:                       # COMMENT/HISTORY/END style
        return f"{keyword:<80}"[:80]
    body = f"{keyword:<8}= {value:>20}"
    if comment:
        body += f" / {comment}"
    return f"{body:<80}"[:80]


def block(*cards: str) -> str:
    """Join cards and append END, as a real header block does."""
    return "".join(cards) + card("END")


#: The real RLMT compressed header, reduced to the cards that matter.
#: NAXIS1=8 is the table ROW LENGTH IN BYTES and NAXIS2=3211 the ROW COUNT;
#: the true image is 4800x3211 and lives in ZNAXIS1/ZNAXIS2.
COMPRESSED = block(
    card("XTENSION", "'BINTABLE'", "binary table extension"),
    card("BITPIX", "8", "8-bit bytes"),
    card("NAXIS", "2", "2-dimensional binary table"),
    card("NAXIS1", "8", "width of table in bytes"),
    card("NAXIS2", "3211", "number of rows in table"),
    card("PCOUNT", "15051851"),
    card("TFIELDS", "1"),
    card("ZIMAGE", "T", "extension contains compressed image"),
    card("ZBITPIX", "16"),
    card("ZNAXIS", "2"),
    card("ZNAXIS1", "4800"),
    card("ZNAXIS2", "3211"),
    card("ZTILE1", "4800"),
    card("ZTILE2", "1"),
    card("ZCMPTYPE", "'RICE_1  '", "compression algorithm"),
    card("FILTER", "'g       '"),
)

#: A genuinely uncompressed image extension — no Z* keywords anywhere.
UNCOMPRESSED = block(
    card("SIMPLE", "T", "conforms to FITS standard"),
    card("BITPIX", "16", "array data type"),
    card("NAXIS", "2"),
    card("NAXIS1", "4788"),
    card("NAXIS2", "3194"),
    card("FILTER", "'V       '"),
)


# ---------------------------------------------------------------------------
# Case 1 — compressed: ZNAXIS wins
# ---------------------------------------------------------------------------
def test_compressed_block_resolves_to_znaxis():
    """The headline regression: 8x3211 must never be reported again."""
    assert fg.geometry_from_card_block(COMPRESSED) == (4800, 3211)


def test_compressed_block_is_detected_as_compressed():
    hdr = fg.parse_card_block(COMPRESSED)
    assert fg.is_compressed_header(hdr) is True
    # And the table bookkeeping is still parsed — we ignore it, not lose it.
    assert hdr["NAXIS1"] == 8 and hdr["NAXIS2"] == 3211


def test_zcmptype_alone_is_enough_to_distrust_naxis():
    """A truncated rescue read that caught ZCMPTYPE but not ZIMAGE must
    still refuse the table NAXIS values."""
    hdr = {"NAXIS1": 8, "NAXIS2": 3211, "ZCMPTYPE": "RICE_1",
           "ZNAXIS1": 4800, "ZNAXIS2": 3211}
    assert fg.resolve_geometry(hdr) == (4800, 3211)


def test_zimage_false_is_treated_as_uncompressed():
    """ZIMAGE = F is a real (if odd) header: NAXIS is then authoritative."""
    hdr = {"ZIMAGE": False, "NAXIS1": 2048, "NAXIS2": 2048}
    assert fg.resolve_geometry(hdr) == (2048, 2048)


def test_catalog_style_float_values_resolve():
    """SQLite hands geometry back as REAL; ints must come out the far side."""
    hdr = {"ZIMAGE": True, "ZNAXIS1": 4800.0, "ZNAXIS2": 3211.0,
           "NAXIS1": 8.0, "NAXIS2": 3211.0}
    assert fg.resolve_geometry(hdr) == (4800, 3211)


# ---------------------------------------------------------------------------
# Case 2 — uncompressed: NAXIS wins
# ---------------------------------------------------------------------------
def test_uncompressed_block_resolves_to_naxis():
    assert fg.geometry_from_card_block(UNCOMPRESSED) == (4788, 3194)


def test_uncompressed_header_is_not_flagged_compressed():
    assert fg.is_compressed_header(fg.parse_card_block(UNCOMPRESSED)) is False


def test_genuinely_small_uncompressed_frame_survives():
    """A REAL sub-frame window must still read as small — the fix must not
    paper over genuine window geometry (the 105x97 guider population)."""
    hdr = {"NAXIS1": 105, "NAXIS2": 97}
    assert fg.resolve_geometry(hdr) == (105, 97)


# ---------------------------------------------------------------------------
# Case 3 — malformed: fails loudly, never silently
# ---------------------------------------------------------------------------
def test_compressed_without_znaxis_raises_not_falls_back():
    """THE critical negative test.  A compressed header missing ZNAXIS must
    raise — silently returning (8, 3211) is the original bug."""
    hdr = {"ZIMAGE": True, "ZCMPTYPE": "RICE_1", "NAXIS1": 8, "NAXIS2": 3211}
    with pytest.raises(fg.GeometryError, match="refusing to fall back"):
        fg.resolve_geometry(hdr)


def test_ragged_card_block_raises():
    """A block that is not a whole number of 80-byte cards is corrupt."""
    with pytest.raises(fg.GeometryError, match="whole number"):
        fg.parse_card_block(COMPRESSED[:-7])


def test_empty_block_raises():
    with pytest.raises(fg.GeometryError, match="empty header block"):
        fg.parse_card_block("")


def test_block_with_no_parsable_cards_raises():
    with pytest.raises(fg.GeometryError, match="no parsable cards"):
        fg.parse_card_block(card("COMMENT nothing here") + card("END"))


def test_missing_dimensions_raise():
    with pytest.raises(fg.GeometryError, match="neither ZNAXIS"):
        fg.resolve_geometry({"BITPIX": 16})


@pytest.mark.parametrize("bad", [0, -4800, "abc", None, 12.5, True])
def test_non_dimension_values_raise(bad):
    """Zero, negative, non-numeric, missing, fractional and logical values
    are all refused rather than coerced into a plausible number."""
    with pytest.raises(fg.GeometryError):
        fg.resolve_geometry({"ZIMAGE": True, "ZNAXIS1": bad, "ZNAXIS2": 3211})


def test_non_ascii_block_raises():
    with pytest.raises(fg.GeometryError, match="not ASCII"):
        fg.parse_card_block(b"\xff" * 80)


# ---------------------------------------------------------------------------
# The rescue path's reason for existing: malformed CONTINUE cards
# ---------------------------------------------------------------------------
def test_malformed_continue_cards_are_skipped_not_fatal():
    """The archive's FWALLNAM cards are followed by CONTINUE cards with
    non-string values, which makes astropy's Header.update raise.  The raw
    parser must step over them and still deliver the geometry."""
    hostile = block(
        card("XTENSION", "'BINTABLE'"),
        card("NAXIS1", "8"),
        card("NAXIS2", "3211"),
        card("ZIMAGE", "T"),
        card("ZCMPTYPE", "'RICE_1  '"),
        card("FWALLNAM", "'Lum&'"),
        card("CONTINUE", "12345"),          # non-string value: astropy dies
        card("CONTINUE", "678"),
        card("ZNAXIS1", "4800"),
        card("ZNAXIS2", "3211"),
    )
    assert fg.geometry_from_card_block(hostile) == (4800, 3211)


def test_first_occurrence_of_a_keyword_wins():
    dupe = block(card("NAXIS1", "4788"), card("NAXIS2", "3194"),
                 card("NAXIS1", "999"))
    assert fg.parse_card_block(dupe)["NAXIS1"] == 4788


# ---------------------------------------------------------------------------
# The scanner's other half: merging cards without losing a frame to one
# malformed card.
#
# The archive's real trigger is a CONTINUE card following a non-string
# FWALLNAM value, which makes astropy raise VerifyError from deep inside
# Card.value.  That exact byte sequence is awkward to synthesize through
# astropy's own constructors (they normalize it away), and pinning a test to
# astropy's internal parsing quirk would make the suite fragile.  So the
# guard is tested against its CONTRACT instead: given a card that raises
# when read, skip that card and keep the frame.  A stub card reproduces
# that condition exactly and cannot drift with astropy versions.
# ---------------------------------------------------------------------------
import importlib.util          # noqa: E402
import sys                     # noqa: E402
from pathlib import Path       # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))


def _build_catalog():
    """Load build_catalog.py by path (it is a CLI script, not a module)."""
    p = (Path(__file__).resolve().parent.parent / "scripts"
         / "build_catalog.py")
    spec = importlib.util.spec_from_file_location("build_catalog", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _GoodCard:
    """A card that reads normally."""

    def __init__(self, keyword, value, comment=""):
        self.keyword, self.value, self.comment = keyword, value, comment


class _PoisonCard:
    """A card that raises when its value is read — what a malformed
    CONTINUE card does inside astropy."""

    keyword = "FWALLNAM"
    comment = ""

    @property
    def value(self):
        from astropy.io.fits.verify import VerifyError
        raise VerifyError("CONTINUE cards must have string values.")


class _FakeHeader:
    def __init__(self, cards):
        self.cards = cards


def test_merge_cards_tolerant_survives_a_hostile_card():
    """One unreadable card must cost that card only — the frame's other
    keywords, and above all its geometry, must still arrive."""
    from astropy.io import fits
    bc = _build_catalog()

    src = _FakeHeader([
        _GoodCard("NAXIS1", 4800),
        _GoodCard("NAXIS2", 3211),
        _GoodCard("FILTER", "g"),
        _PoisonCard(),                    # the frame-killer
        _GoodCard("EGAIN", 1.0),          # must still be reached
    ])
    dest = fits.Header()
    n_bad = bc.merge_cards_tolerant(dest, src)

    assert n_bad == 1                     # exactly the hostile card
    assert dest["NAXIS1"] == 4800         # and the frame survived intact
    assert dest["NAXIS2"] == 3211
    assert dest["FILTER"] == "g"
    assert dest["EGAIN"] == 1.0
    assert "FWALLNAM" not in dest


def test_merge_cards_tolerant_reports_a_clean_header_as_clean():
    """No false positives: a header with no bad cards reports zero skips."""
    from astropy.io import fits
    bc = _build_catalog()
    src = _FakeHeader([_GoodCard("NAXIS1", 4788), _GoodCard("NAXIS2", 3194)])
    dest = fits.Header()
    assert bc.merge_cards_tolerant(dest, src) == 0
    assert (dest["NAXIS1"], dest["NAXIS2"]) == (4788, 3194)


def test_scanner_geometry_comes_from_the_resolver():
    """End-to-end on the merged header: a scanner that merged a RAW table
    header (astropy could not build a CompImageHDU) must still record the
    true image size, not the table's row bookkeeping."""
    from astropy.io import fits
    merged = fits.Header()
    merged["NAXIS1"] = 8
    merged["NAXIS2"] = 3211
    merged["ZIMAGE"] = True
    merged["ZCMPTYPE"] = "RICE_1"
    merged["ZNAXIS1"] = 4800
    merged["ZNAXIS2"] = 3211
    assert fg.resolve_geometry(merged) == (4800, 3211)
