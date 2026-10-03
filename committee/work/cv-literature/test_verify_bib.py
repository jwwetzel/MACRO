"""Offline tests for verify_bib.py -- the parser and the comparison logic, on synthetic input whose
right answer is known by construction.  No network: registry lookups are not exercised here (the live
run is ``python verify_bib.py``, whose table is ``bib_verification.md``).

Run:  /opt/miniconda3/envs/rlmt-checks/bin/python -m pytest committee/work/cv-literature/test_verify_bib.py -q
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_bib as vb  # noqa: E402

BIB = r'''
%% a comment, with braces { that must not confuse the parser
@ARTICLE{cropper1986,
  author  = {{Cropper}, M.},
  title   = "{Polarimetry of ST LMi (CW1103+254)}",
  journal = {\mnras},
  year    = 1986,
  volume  = 222,
  pages   = {853--870},
  doi     = {10.1093/mnras/222.4.853}
}

@MISC{collab,
  author = {{Gaia Collaboration} and {G{\"a}nsicke}, B.~T. and others},
  title  = "{The {G}aia mission}",
  year   = 2016
}
'''

#: A registry record as crossref() would reduce it, matching the first entry above.
REC = {"registry": "Crossref", "first_family": "cropper", "years": [1986], "volume": "222",
       "pages": ["853"], "title": "Polarimetry of ST LMi (CW1103 + 254)", "container": "MNRAS"}


def test_parser_reads_both_value_styles_and_bare_numbers():
    entries = vb.parse_bib(BIB)
    assert [e["key"] for e in entries] == ["cropper1986", "collab"]
    e = entries[0]
    assert e["type"] == "ARTICLE" and e["year"] == "1986" and e["volume"] == "222"
    assert e["title"] == "{Polarimetry of ST LMi (CW1103+254)}"
    assert e["doi"] == "10.1093/mnras/222.4.853"


def test_first_family_skips_collaborations_and_folds_accents():
    entries = vb.parse_bib(BIB)
    assert vb.first_family(entries[0]["author"]) == "cropper"
    assert vb.first_family(entries[1]["author"]) == "gansicke"


def test_matching_record_has_no_disagreement():
    assert vb.compare(vb.parse_bib(BIB)[0], REC) == []


def test_each_wrong_field_is_named():
    entry = vb.parse_bib(BIB)[0]
    for field, wrong, word in (("first_family", "kato", "first author"), ("years", [1994], "year"),
                               ("volume", "66", "volume"), ("pages", ["30"], "first page"),
                               ("title", "Thermal and non-thermal X-rays from a supernova remnant", "title overlap")):
        bad = vb.compare(entry, {**REC, field: wrong})
        assert len(bad) == 1 and bad[0].startswith(word), (field, bad)


def test_a_doi_pointing_at_another_paper_is_caught():
    """The failure this script exists for: a recalled DOI that resolves, but to somebody else's paper."""
    entry = {"key": "kato2014", "author": "{{Kato}, T. and others}", "year": "2014", "volume": "66", "pages": "30",
             "title": "Survey of period variations of superhumps in SU UMa-type dwarf novae. V."}
    other = {"registry": "Crossref", "first_family": "yamauchi", "years": [2014], "volume": "66", "pages": ["2"],
             "title": "Thermal and non-thermal X-rays from the Galactic supernova remnant G348.5+0.1", "container": "PASJ"}
    bad = vb.compare(entry, other)
    assert any(b.startswith("first author") for b in bad) and any(b.startswith("title overlap") for b in bad)


def test_year_tolerates_print_versus_online_but_not_more():
    entry = vb.parse_bib(BIB)[0]
    assert vb.compare(entry, {**REC, "years": [1985, 1987]}) == []
    assert vb.compare(entry, {**REC, "years": [1984]}) != []
