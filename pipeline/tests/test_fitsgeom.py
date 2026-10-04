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


def test_genuinely_small_plain_image_survives():
    """A REAL sub-frame window in a PLAIN (uncompressed) image must still
    read as small — the fix must not paper over genuine window geometry."""
    hdr = {"NAXIS1": 105, "NAXIS2": 97}
    assert fg.resolve_geometry(hdr) == (105, 97)


def test_genuinely_small_COMPRESSED_frame_survives():
    """The case the real control group actually exercises, and the one that
    must never regress.

    All 91 rows the re-scan left unchanged are tile-compressed ``.fts.fz``
    files — Andor iKon focus and guide windows.  Their headers therefore
    carry BOTH a BINTABLE ``NAXIS1`` (a row length in bytes, small) AND a
    genuinely small ``ZNAXIS1``, and the resolver has to tell those two
    small numbers apart: read the wrong one and a 45x34 focus window becomes
    an 8x34 phantom, which is the very artifact this module exists to undo.

    Testing this against a plain uncompressed header instead would prove
    only the trivial case — the same slip that once described the control
    group as 'uncompressed' in the report prose."""
    hdr = {"ZIMAGE": True, "ZCMPTYPE": "RICE_1",
           "NAXIS1": 8, "NAXIS2": 34,          # the BINTABLE's own shape
           "ZNAXIS1": 45, "ZNAXIS2": 34}       # the real image
    assert fg.is_compressed_header(hdr) is True
    assert fg.resolve_geometry(hdr) == (45, 34)


def test_compressed_small_and_compressed_phantom_are_distinguished():
    """Side by side: the SAME container and the same machinery must give a
    small answer for a genuine window and a full-frame answer for a phantom.
    That contrast is what makes the 91-row control group evidence at all."""
    window = {"ZIMAGE": True, "ZCMPTYPE": "RICE_1", "NAXIS1": 8,
              "NAXIS2": 48, "ZNAXIS1": 57, "ZNAXIS2": 48}
    phantom = {"ZIMAGE": True, "ZCMPTYPE": "RICE_1", "NAXIS1": 8,
               "NAXIS2": 3211, "ZNAXIS1": 4800, "ZNAXIS2": 3211}
    assert fg.resolve_geometry(window) == (57, 48)
    assert fg.resolve_geometry(phantom) == (4800, 3211)


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


# ===========================================================================
# Hardware-state cards (finding F-2, plan review 2026-10-03)
# ===========================================================================
# Same discipline as above: every header is a synthetic card block, built
# byte for byte from the dialects the real archive contains.

def raw_card(text: str) -> str:
    """One card image from literal text, padded to exactly 80 characters —
    for the malformed cards the tidy ``card()`` helper cannot spell."""
    assert len(text) <= 80, text
    return f"{text:<80}"


def pad_block(cards: str) -> bytes:
    """Pad a card string to a whole number of 2880-byte FITS blocks."""
    n = -(-len(cards) // fg.BLOCK_LEN) * fg.BLOCK_LEN
    return f"{cards:<{n}}".encode("ascii")


#: MaxIm/ASI header (2024-12): standard long-string wheel map, doubled quotes.
ASI_HEADER = block(
    raw_card("XTENSION= 'BINTABLE'           / binary table extension"),
    raw_card("SET-TEMP=  -10.000000000000000 /CCD temperature setpoint in C"),
    raw_card("CCD-TEMP=  -10.000000000000000 /CCD temperature at start"),
    raw_card("XPIXSZ  =   7.5199999999999996 /Pixel Width in microns"),
    raw_card("GAIN    =                  100"),
    raw_card("OFFSET  =                   30"),
    raw_card("SWCREATE= 'MaxIm DL Version 6.40 231201 0HTQT' /Name of software"),
    raw_card("FOCUSPOS=                 9980 /Focuser position in steps"),
    raw_card("INSTRUME= 'ASI Camera (1)'"),
    raw_card("FLIPSTAT= 'Flip/Mirror'"),
    raw_card("COOLPOWR=                   28 / Cooler power in percent"),
    raw_card("CAMNAME = 'ASCOM   '           / Name of camera"),
    raw_card("TELPIER = 'pierEast'           / Telescope pier side"),
    raw_card("FWPOS   =                    0 / Filter wheel position"),
    raw_card("FWNAME  = 'Dual Wheels'        / Filter wheel name"),
    raw_card("FWALLNAM= '(''g'', ''OGGrism'', ''r'', ''i'', ''HaGrism'', &'"),
    raw_card("CONTINUE  '''z'', ''ha'')&'"),
    raw_card("CONTINUE  '' / Filter wheel names"),
    raw_card("FOCPOS  =     9980.31127732742 / Focuser position"),
)

#: pyscope-native header (2026-06-29): the wheel map is written with RAW
#: embedded quotes — an illegal card that astropy reads as the string '('.
PYSCOPE_HEADER = block(
    raw_card("SWCREATE= 'pyscope'            / Software used to create file"),
    raw_card("GAIN    =                   56 / Electronic gain"),
    raw_card("OFFSET  = ''                   / Image offset"),
    raw_card("COOLPOWR=   7.8431372549019605 / Cooler power in percent"),
    raw_card("SET-TEMP=                 -0.0 / Camera temperature setpoint [C]"),
    raw_card("CAMNAME = 'QHY600MPCIE-ec89ab488ac5aca73' / Name of camera"),
    raw_card("FWALLNAM= '('g', 'lrg', 'r', 'i', 'ha', 'hrg', 'empty')' / Filter"),
    raw_card("FOCPOS  =     8127.22652775466 / Focuser position"),
)


class TestParseHardwareCards:
    def test_standard_header_values_come_back_as_written(self):
        got = fg.parse_hardware_cards(ASI_HEADER)
        assert got["INSTRUME"] == "ASI Camera (1)"
        assert got["GAIN"] == "100" and got["OFFSET"] == "30"
        assert got["SET-TEMP"] == "-10.000000000000000"
        assert got["FLIPSTAT"] == "Flip/Mirror"
        assert got["TELPIER"] == "pierEast" and got["FWPOS"] == "0"
        assert got["FOCPOS"] == "9980.31127732742"
        assert got["FOCUSPOS"] == "9980"

    def test_the_wheel_map_is_followed_through_its_continue_cards(self):
        # The card astropy cannot parse is the one that records the wheel.
        got = fg.parse_hardware_cards(ASI_HEADER)
        assert got["FWALLNAM"] == \
            "('g', 'OGGrism', 'r', 'i', 'HaGrism', 'z', 'ha')"

    def test_pyscope_raw_quote_dialect(self):
        got = fg.parse_hardware_cards(PYSCOPE_HEADER)
        assert got["FWALLNAM"] == \
            "('g', 'lrg', 'r', 'i', 'ha', 'hrg', 'empty')"
        assert got["CAMNAME"] == "QHY600MPCIE-ec89ab488ac5aca73"
        assert got["SWCREATE"] == "pyscope"

    def test_blank_card_is_empty_string_absent_card_is_absent(self):
        """The distinction the acceptance test turns on."""
        blanks = block(
            raw_card("OFFSET  =  / Image offset"),
            raw_card("FLIPSTAT= '        '"),
            raw_card("GAIN    = '4x      '"),
        )
        got = fg.parse_hardware_cards(blanks)
        assert got["OFFSET"] == ""          # present, no value
        assert got["FLIPSTAT"] == ""        # present, blank — a measurement
        assert got["GAIN"] == "4x"          # the iKon's preamp setting
        assert "TELPIER" not in got         # absent
        assert "OFFSET" in fg.parse_hardware_cards(PYSCOPE_HEADER)
        assert fg.parse_hardware_cards(PYSCOPE_HEADER)["OFFSET"] == ""

    def test_first_occurrence_wins_and_unwanted_cards_are_ignored(self):
        twice = block(raw_card("GAIN    =                  100"),
                      raw_card("GAIN    =                   56"),
                      raw_card("OBJECT  = 'T CrB'"))
        assert fg.parse_hardware_cards(twice) == {"GAIN": "100"}

    def test_a_continue_with_no_open_string_is_harmless(self):
        stray = block(raw_card("GAIN    =                  100"),
                      raw_card("CONTINUE  'orphan' / nothing to continue"))
        assert fg.parse_hardware_cards(stray) == {"GAIN": "100"}

    def test_not_card_structured_raises(self):
        with pytest.raises(fg.GeometryError):
            fg.parse_hardware_cards("GAIN = 100")
        with pytest.raises(fg.GeometryError):
            fg.parse_hardware_cards("")


class TestReadHardwareCards:
    def _fpack_like(self, tmp_path, ext_header):
        """A file laid out like an fpack output: a dataless primary HDU,
        then the extension that carries the real header."""
        primary = block(card("SIMPLE", "T"), card("BITPIX", "16"),
                        card("NAXIS", "0"), card("EXTEND", "T"))
        path = tmp_path / "f.fts.fz"
        path.write_bytes(pad_block(primary) + pad_block(ext_header)
                         + b"\0" * 2880)
        return str(path)

    def test_reads_the_extension_of_a_dataless_primary(self, tmp_path):
        got = fg.read_hardware_cards(self._fpack_like(tmp_path, ASI_HEADER))
        assert got["INSTRUME"] == "ASI Camera (1)"
        assert got["FWALLNAM"].startswith("('g', 'OGGrism'")

    def test_plain_file_is_read_from_its_primary(self, tmp_path):
        primary = block(card("SIMPLE", "T"), card("BITPIX", "16"),
                        card("NAXIS", "2"), card("NAXIS1", "4"),
                        card("NAXIS2", "4"),
                        raw_card("INSTRUME= 'Andor CCD/EMCCD (SDK2)'"),
                        raw_card("GAIN    = '4x      '"))
        path = tmp_path / "m.fts"
        path.write_bytes(pad_block(primary) + b"\0" * 2880)
        got = fg.read_hardware_cards(str(path))
        assert got == {"INSTRUME": "Andor CCD/EMCCD (SDK2)", "GAIN": "4x"}

    def test_a_header_spanning_several_blocks(self, tmp_path):
        # 40 filler cards push the wanted card into the second block.
        filler = "".join(raw_card(f"FILL{i:04d}=                    1")
                         for i in range(40))
        ext = filler + block(raw_card("TELPIER = 'pierWest'"))
        got = fg.read_hardware_cards(self._fpack_like(tmp_path, ext))
        assert got["TELPIER"] == "pierWest"

    def test_gzip_files_are_read_through_gzip(self, tmp_path):
        import gzip
        primary = block(card("SIMPLE", "T"), card("NAXIS", "2"),
                        raw_card("FLIPSTAT= 'Flip/Mirror'"))
        path = tmp_path / "g.fts.gz"
        with gzip.open(path, "wb") as fh:
            fh.write(pad_block(primary))
        assert fg.read_hardware_cards(str(path)) == \
            {"FLIPSTAT": "Flip/Mirror"}

    def test_truncated_and_empty_files_fail_loudly(self, tmp_path):
        empty = tmp_path / "e.fts"
        empty.write_bytes(b"")
        with pytest.raises(fg.GeometryError):
            fg.read_hardware_cards(str(empty))
        cut = tmp_path / "c.fts"
        cut.write_bytes(b"SIMPLE  =                    T" + b" " * 100)
        with pytest.raises(fg.GeometryError):
            fg.read_hardware_cards(str(cut))

    def test_a_header_with_no_end_card_is_refused_not_read_forever(
            self, tmp_path, monkeypatch):
        monkeypatch.setattr(fg, "MAX_HEADER_BLOCKS", 3)
        path = tmp_path / "n.fts"
        path.write_bytes(pad_block(card("SIMPLE", "T")) * 5)
        with pytest.raises(fg.GeometryError, match="no END card"):
            fg.read_hardware_cards(str(path))


# ---------------------------------------------------------------------------
# The re-scrape driver (rescan_geometry.py hdr-*) on a hand-built catalog
# ---------------------------------------------------------------------------
import sqlite3                                               # noqa: E402
import sys                                                   # noqa: E402
from pathlib import Path                                     # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import rescan_geometry as rg                                 # noqa: E402


@pytest.fixture
def toy_archive(tmp_path, monkeypatch):
    """A three-file archive and its catalog, with the script pointed at it.

    ``a`` is a healthy ASI frame the original scan read fully; ``b`` is one
    whose INSTRUME/SWCREATE the original scan LOST (NULL in ``obs``) though
    the cards are in the file; ``c`` does not exist on disk and was already
    an error row in the original scan.
    """
    root = tmp_path / "archive"
    (root / "rawimage").mkdir(parents=True)
    primary = block(card("SIMPLE", "T"), card("BITPIX", "16"),
                    card("NAXIS", "0"), card("EXTEND", "T"))
    for name in ("a", "b"):
        (root / "rawimage" / f"{name}.fts.fz").write_bytes(
            pad_block(primary) + pad_block(ASI_HEADER))
    db = tmp_path / "catalog.sqlite"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE obs (path TEXT PRIMARY KEY, tree TEXT, "
                "error TEXT, instrume TEXT, swcreate TEXT, ccd_temp REAL, "
                "focuspos REAL)")
    con.executemany("INSERT INTO obs VALUES (?,?,?,?,?,?,?)", [
        ("rawimage/a.fts.fz", "rawimage", None, "ASI Camera (1)",
         "MaxIm DL Version 6.40 231201 0HTQT", -10.0, 9980.0),
        ("rawimage/b.fts.fz", "rawimage", None, None, None, None, None),
        ("reduced/c.fts.fz", "reduced", "OSError: Empty or corrupt FITS file",
         None, None, None, None),
    ])
    con.commit()
    con.close()
    monkeypatch.setattr(rg, "ROOT", str(root))
    monkeypatch.setattr(rg, "DB", str(db))
    return db


class _Args:
    workers, limit, retry_errors = 2, None, False


class TestHeaderRescrapeDriver:
    def test_column_names_avoid_the_sql_keyword(self):
        assert rg.hdr_column("OFFSET") == "h_offset"
        assert rg.hdr_column("SET-TEMP") == "h_set_temp"
        assert len(rg.HDR_COLUMNS) == len(fg.HARDWARE_CARDS)

    def test_run_is_additive_resumable_and_never_touches_obs(
            self, toy_archive):
        con = sqlite3.connect(toy_archive)
        before = con.execute("SELECT * FROM obs ORDER BY path").fetchall()
        con.close()
        assert rg.cmd_hdr_run(_Args()) == 0
        con = sqlite3.connect(toy_archive)
        assert con.execute("SELECT * FROM obs ORDER BY path").fetchall() \
            == before, "the re-scrape must not write one byte of obs"
        rows = dict(con.execute(
            "SELECT path, error IS NULL FROM hdr_rescrape"))
        assert rows == {"rawimage/a.fts.fz": 1, "rawimage/b.fts.fz": 1,
                        "reduced/c.fts.fz": 0}
        got = con.execute(
            "SELECT h_gain, h_offset, h_set_temp, h_flipstat, h_fwallnam, "
            "n_cards FROM hdr_rescrape WHERE path = 'rawimage/b.fts.fz'"
        ).fetchone()
        assert got[:4] == ("100", "30", "-10.000000000000000", "Flip/Mirror")
        assert got[4].startswith("('g', 'OGGrism'")
        assert got[5] == 16                    # every wanted card was there
        # Resumable: a second run finds nothing left to read.
        assert rg.hdr_todo(con) == []
        # ... unless asked to retry the unreadable file.
        assert rg.hdr_todo(con, retry_errors=True) == ["reduced/c.fts.fz"]
        con.close()

    def test_read_order_is_rawimage_first_reduced_last(self, toy_archive):
        con = sqlite3.connect(toy_archive)
        con.execute("INSERT INTO obs (path, tree) VALUES "
                    "('iKon/z.fts.fz', 'iKon')")
        rg.ensure_hdr_table(con)
        assert rg.hdr_todo(con) == [
            "rawimage/a.fts.fz", "rawimage/b.fts.fz", "iKon/z.fts.fz",
            "reduced/c.fts.fz"]
        con.close()

    def test_verify_counts_the_control_group_and_the_recovered_nulls(
            self, toy_archive):
        rg.cmd_hdr_run(_Args())
        con = sqlite3.connect(toy_archive)
        c = rg.hdr_verify_counts(con)
        con.close()
        assert (c["n_obs"], c["n_scanned"], c["n_unreadable"]) == (3, 3, 1)
        # The unreadable file was unreadable to the original scan too.
        assert c["n_unreadable_new"] == 0
        # Control group: file 'a', the four cards both scans read, agree.
        assert c["control"]["h_instrume"] == (1, 0)
        assert c["control"]["h_ccd_temp"] == (1, 0)
        assert c["control"]["h_focuspos"] == (1, 0)
        # F-2 acceptance, measured: file 'b' had NULLs where cards exist.
        assert c["recovered"]["h_instrume"] == 1
        assert c["recovered"]["h_swcreate"] == 1
        assert c["recovered"]["h_ccd_temp"] == 1
        # value / blank / absent over the two readable rows.
        assert c["cards"]["h_gain"] == (2, 0, 0)
        assert rg.cmd_hdr_verify(_Args()) == 0

    def test_verify_fails_when_the_two_scans_disagree(self, toy_archive):
        rg.cmd_hdr_run(_Args())
        con = sqlite3.connect(toy_archive)
        con.execute("UPDATE obs SET instrume = 'Some Other Camera' "
                    "WHERE path = 'rawimage/a.fts.fz'")
        con.commit()
        assert rg.hdr_verify_counts(con)["control"]["h_instrume"] == (1, 1)
        con.close()
        assert rg.cmd_hdr_verify(_Args()) == 1

    def test_verify_fails_when_rows_are_missing(self, toy_archive):
        _Args.limit = 1
        try:
            rg.cmd_hdr_run(_Args())
        finally:
            _Args.limit = None
        assert rg.cmd_hdr_verify(_Args()) == 1
