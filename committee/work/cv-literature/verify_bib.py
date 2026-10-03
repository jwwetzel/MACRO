#!/usr/bin/env python
"""Verify every entry of ``manuscripts/CV_TimeSeries/references.bib`` against its registry of record.

WHY THIS SCRIPT EXISTS
----------------------
``references.bib`` is hand-maintained (a citation is an editorial act, not a measurement), and the
committee's rule for this package is that no entry may be invented or mis-transcribed.  A hand-typed
bibliography cannot promise that; a check against the registries can.  This script is that check, and it
is the only source of the verification table in ``bib_verification.md`` -- nobody types a "verified".

(The rule earns its keep: while drafting, two DOIs recalled from memory resolved to unrelated papers, and
this check is what caught them.)

WHAT IS CHECKED, PER ENTRY (first method that applies)
------------------------------------------------------
1. ``doi``     -> Crossref ``/works/<doi>`` (DataCite for 10.5281 Zenodo DOIs).  The registry record must
                  agree with the entry on: first-author family name, year (within one, because print and
                  online years differ), volume (when both give one), first page or article number (when
                  both give one), and title (token overlap >= 0.6).
2. ``eprint``  -> arXiv API ``id_list``.  First-author family name and title must agree.
3. ``adsurl``  -> the ADS full-text scan ``https://adsabs.harvard.edu/pdf/<bibcode>`` must return a PDF.
                  A non-existent bibcode returns 403 there, which is what makes this a test.
4. ``url``     -> the page must load and contain the string registered in ``URL_MUST_CONTAIN``.
5. Entries no registry can answer for are listed in ``MANUAL`` with the evidence a human checked.
   Zenodo concept DOIs (``CONCEPT_DOI``) are checked for resolution and first author only.

Every disagreement is printed with the field that disagreed; a single MISMATCH or UNRESOLVED makes the
exit status non-zero.

USAGE
-----
    cd <repo>
    /opt/miniconda3/envs/rlmt-checks/bin/python committee/work/cv-literature/verify_bib.py
    # optional: --bib PATH  --out PATH  --refresh (ignore the cache)

Network access is required on the first run.  Registry answers are cached (reduced to the compared
fields) in ``bib_verification_cache.json`` beside this script, so a re-run is fast and reproducible and
the evidence a verdict rested on stays inspectable.  Only the standard library is used.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_BIB = REPO / "manuscripts" / "CV_TimeSeries" / "references.bib"
DEFAULT_OUT = HERE / "bib_verification.md"
CACHE_PATH = HERE / "bib_verification_cache.json"

#: Polite-pool contact for Crossref (their etiquette asks for one).
MAILTO = "jwetzel@coe.edu"
USER_AGENT = f"MACRO-cv-literature-verify/1.0 (mailto:{MAILTO})"

#: Title agreement threshold: fraction of the entry's title tokens found in the registry title.
#: 0.6 tolerates LaTeX markup, registry subtitles held in a separate field and ALL-CAPS records,
#: while still rejecting a DOI that points at a different paper (those score near zero).
TITLE_OVERLAP_MIN = 0.6

#: Entries whose ``url`` is the evidence: the page must contain this string.
URL_MUST_CONTAIN = {
    # Konkoly Observatory's own archive of its Mitteilungen; ADS holds the bibcode but serves no scan.
    "walker1965": "Photoelectric Observations of VV Puppis",
}

#: Entries whose DOI is a Zenodo *concept* DOI.  A concept DOI always resolves to the newest release,
#: so the registry's year and title describe a later version than the one cited and cannot be compared;
#: only resolution and the first author are checked.  (The fix is a version DOI -- see the entry's note.)
CONCEPT_DOI = {
    "photutils": ("Zenodo concept DOI: resolves to the latest release, so year/title are not comparable. "
                  "NOTE the pipeline does not import photutils at all (photometry is `sep`); the citation "
                  "should leave main.tex -- see citation_map.md."),
}

#: Entries no registry can verify, with what was checked by hand and when.
MANUAL = {
    "aavso_aid": ("AAVSO's own citation template for the International Database (AAVSO Data Usage "
                  "Guidelines dated 2025-10-01, https://www.aavso.org/data-usage-guidelines); aavso.org "
                  "refuses scripted requests (Cloudflare challenge), so this cannot be fetched here. "
                  "The year must match the year the data were downloaded -- chair/James to confirm."),
}


# --------------------------------------------------------------------------------------------------
# BibTeX parsing -- deliberately minimal: this file's own conventions, not the whole grammar.
# --------------------------------------------------------------------------------------------------
def parse_bib(text: str) -> list[dict]:
    """Return one dict per entry: ``{"type", "key", <field>: <value>, ...}``.

    Handles the two value styles the file uses -- ``{...}`` (brace-balanced), ``"{...}"`` -- and bare
    numbers.  Lines starting with ``%`` outside entries are comments and are skipped by construction,
    because only text following an at-sign is read.
    """
    entries = []
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", text):
        etype, key = m.group(1).upper(), m.group(2)
        i, depth, start = m.end(), 1, m.end()
        while i < len(text) and depth:            # find the entry's closing brace
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            i += 1
        body = text[start:i - 1]
        entry = {"type": etype, "key": key}
        j = 0
        while True:
            fm = re.compile(r"\s*([A-Za-z]+)\s*=\s*").match(body, j)
            if not fm:
                break
            name, j = fm.group(1).lower(), fm.end()
            if body[j] in '{"':
                opener = body[j]
                k, d = j + 1, 1
                if opener == "{":
                    while d:
                        d += {"{": 1, "}": -1}.get(body[k], 0)
                        k += 1
                else:                              # "..." with braces allowed inside
                    while body[k] != '"' or d != 1:
                        d += {"{": 1, "}": -1}.get(body[k], 0)
                        k += 1
                    k += 1
                value = body[j + 1:k - 1]
                j = k
            else:
                vm = re.compile(r"[^,\s]+").match(body, j)
                value, j = vm.group(0), vm.end()
            entry[name] = value.strip()
            cm = re.compile(r"\s*,").match(body, j)
            if not cm:
                break
            j = cm.end()
        entries.append(entry)
    return entries


# --------------------------------------------------------------------------------------------------
# Normalisation helpers
# --------------------------------------------------------------------------------------------------
_TEX_ACCENT = re.compile(r"\\[`'^\"~=.Hcvuk]\s*\{?\\?([A-Za-z])\}?")


def fold(s: str) -> str:
    """Lower-case ASCII skeleton of a string: LaTeX accents, markup and punctuation removed."""
    s = _TEX_ACCENT.sub(r"\1", s or "")
    s = s.replace(r"\i", "i").replace(r"\l", "l").replace(r"\o", "o").replace(r"\ss", "ss")
    s = re.sub(r"<[^>]+>", " ", s)                       # registry HTML tags
    s = re.sub(r"\\[A-Za-z]+", " ", s)                   # remaining TeX commands
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def first_family(author_field: str) -> str:
    """Family name of the first *person* in a BibTeX author field (collaboration names skipped).

    The family name is everything before the first comma at brace depth zero, so names that carry
    braces of their own -- ``{G{\\"a}nsicke}, B.~T.`` -- are read whole and then folded.
    """
    for person in re.split(r"\s+and\s+", author_field or ""):
        if "collaboration" in person.lower():
            continue
        depth, family = 0, person
        for i, ch in enumerate(person):
            depth += {"{": 1, "}": -1}.get(ch, 0)
            if ch == "," and depth == 0:
                family = person[:i]
                break
        if fold(family):
            return fold(family).replace(" ", "")
    return ""


def title_overlap(entry_title: str, registry_title: str) -> float:
    """Fraction of the entry's title tokens (length > 2) present in the registry title."""
    want = [t for t in fold(entry_title).split() if len(t) > 2]
    have = set(fold(registry_title).split())
    return sum(t in have for t in want) / len(want) if want else 0.0


def first_page(pages: str) -> str:
    return re.split(r"-+", (pages or "").strip())[0].strip().lower()


# --------------------------------------------------------------------------------------------------
# Registry lookups (each returns a reduced record, cached)
# --------------------------------------------------------------------------------------------------
def _get(url: str, timeout: float = 60.0, tries: int = 3) -> tuple[int, bytes]:
    """HTTP GET with retries; returns (status, body).  Status 0 means no answer at all."""
    last = 0
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            last = e.code
            if e.code in (403, 404):               # a definite "no": do not hammer the server
                return e.code, b""
        except Exception:                           # timeouts, resets, 5xx gateways
            last = 0
        time.sleep(2.0 * (attempt + 1))
    return last, b""


def crossref(doi: str) -> dict | None:
    status, body = _get("https://api.crossref.org/works/" + urllib.parse.quote(doi) + f"?mailto={MAILTO}")
    if status != 200:
        return None
    m = json.loads(body)["message"]
    people = [a for a in m.get("author", []) if "family" in a]
    years = sorted({(m.get(k) or {}).get("date-parts", [[None]])[0][0]
                    for k in ("published-print", "published-online", "issued")} - {None})
    title = " ".join((m.get("title") or [""]) + (m.get("subtitle") or []))
    return {"registry": "Crossref", "first_family": fold(people[0]["family"]) if people else "",
            "years": years, "volume": m.get("volume") or "",
            "pages": [p for p in (first_page(m.get("page", "")), (m.get("article-number") or "").lower()) if p],
            "title": re.sub(r"\s+", " ", title), "container": (m.get("container-title") or [""])[0]}


def datacite(doi: str) -> dict | None:
    status, body = _get("https://api.datacite.org/dois/" + urllib.parse.quote(doi))
    if status != 200:
        return None
    a = json.loads(body)["data"]["attributes"]
    creators = a.get("creators") or [{}]
    fam = creators[0].get("familyName") or (creators[0].get("name") or "").split(",")[0]
    return {"registry": "DataCite", "first_family": fold(fam), "years": [a.get("publicationYear")],
            "volume": "", "pages": [], "title": (a.get("titles") or [{}])[0].get("title", ""),
            "container": a.get("publisher") if isinstance(a.get("publisher"), str) else "Zenodo"}


def arxiv(eprint: str) -> dict | None:
    status, body = _get("https://export.arxiv.org/api/query?id_list=" + urllib.parse.quote(eprint))
    if status != 200 or b"<entry>" not in body:
        return None
    e = body.decode("utf-8", "replace").split("<entry>")[1]
    if "<title>Error</title>" in e:
        return None
    title = re.search(r"<title>(.*?)</title>", e, re.S).group(1)
    names = re.findall(r"<name>(.*?)</name>", e)
    return {"registry": "arXiv", "first_family": fold(names[0].split()[-1]) if names else "",
            "years": [int(re.search(r"<published>(\d{4})", e).group(1))], "volume": "", "pages": [],
            "title": re.sub(r"\s+", " ", title).strip(), "container": "arXiv"}


def ads_scan(adsurl: str) -> dict | None:
    """True existence test of a bibcode: ADS serves its scan as a PDF (a bogus bibcode gives 403)."""
    bibcode = urllib.parse.unquote(adsurl.rstrip("/").split("/abs/")[-1].split("/")[0])
    status, body = _get("https://adsabs.harvard.edu/pdf/" + urllib.parse.quote(bibcode), timeout=120, tries=4)
    if status == 200 and body[:4] == b"%PDF":
        return {"registry": "ADS scan", "bibcode": bibcode, "bytes": len(body)}
    return None


def url_contains(url: str, needle: str) -> dict | None:
    status, body = _get(url)
    if status == 200 and fold(needle) in fold(body.decode("utf-8", "replace")):
        return {"registry": "publisher page", "url": url, "found": needle}
    return None


# --------------------------------------------------------------------------------------------------
# The verdict
# --------------------------------------------------------------------------------------------------
def compare(entry: dict, rec: dict) -> list[str]:
    """Return the list of disagreements between a bib entry and a registry record (empty = agree)."""
    bad = []
    fam = first_family(entry.get("author", ""))
    # Registry family names may carry particles/ordering differently ("de Miguel" vs "Miguel"):
    # accept containment either way, reject anything else.
    reg = rec["first_family"].replace(" ", "")
    if fam and reg and not (fam in reg or reg in fam):
        bad.append(f"first author: bib '{fam}' vs registry '{rec['first_family']}'")
    year = int(re.sub(r"\D", "", entry.get("year", "0")) or 0)
    if rec["years"] and all(abs(year - y) > 1 for y in rec["years"] if y):
        bad.append(f"year: bib {year} vs registry {rec['years']}")
    if entry.get("volume") and rec["volume"] and entry["volume"] != rec["volume"]:
        bad.append(f"volume: bib {entry['volume']} vs registry {rec['volume']}")
    page = first_page(entry.get("pages", ""))
    if page and rec["pages"] and not any(page.lstrip("0") == p.lstrip("0") or page in p or p in page
                                         for p in rec["pages"]):
        bad.append(f"first page: bib {page} vs registry {rec['pages']}")
    ov = title_overlap(entry.get("title", ""), rec["title"])
    if ov < TITLE_OVERLAP_MIN:
        bad.append(f"title overlap {ov:.2f} < {TITLE_OVERLAP_MIN}: registry '{rec['title'][:80]}'")
    return bad


def verify(entry: dict, cache: dict, refresh: bool) -> tuple[str, str, str]:
    """Return (status, method, detail) for one entry.  status in VERIFIED/MISMATCH/UNRESOLVED/MANUAL."""
    key = entry["key"]
    if key in MANUAL:
        return "MANUAL", "by hand", MANUAL[key]

    def cached(tag: str, ident: str, fn):
        ck = f"{tag}:{ident}"
        if refresh or ck not in cache:
            cache[ck] = fn()
            time.sleep(0.15)
        return cache[ck]

    if entry.get("doi"):
        doi = entry["doi"]
        rec = cached("doi", doi, lambda: (datacite(doi) if doi.startswith("10.5281/") else crossref(doi)))
        if rec is None:
            return "UNRESOLVED", "doi", f"{doi} did not resolve at the registry"
        bad = compare(entry, rec)
        if key in CONCEPT_DOI:
            bad = [x for x in bad if x.startswith("first author")]
            return (("MISMATCH", rec["registry"], "; ".join(bad)) if bad
                    else ("VERIFIED", rec["registry"], f"{doi} -- {CONCEPT_DOI[key]}"))
        return ("MISMATCH", rec["registry"], "; ".join(bad)) if bad else ("VERIFIED", rec["registry"], doi)
    if entry.get("eprint"):
        rec = cached("arxiv", entry["eprint"], lambda: arxiv(entry["eprint"]))
        if rec is None:
            return "UNRESOLVED", "arXiv", f"arXiv:{entry['eprint']} not found"
        bad = [b for b in compare(entry, rec) if not b.startswith("year")]   # journal year != posting year
        return ("MISMATCH", "arXiv", "; ".join(bad)) if bad else ("VERIFIED", "arXiv", "arXiv:" + entry["eprint"])
    if key in URL_MUST_CONTAIN and entry.get("url"):
        rec = cached("url", entry["url"], lambda: url_contains(entry["url"], URL_MUST_CONTAIN[key]))
        return (("VERIFIED", "publisher page", entry["url"]) if rec
                else ("UNRESOLVED", "url", f"{entry['url']} did not show '{URL_MUST_CONTAIN[key]}'"))
    if entry.get("adsurl"):
        rec = cached("ads", entry["adsurl"], lambda: ads_scan(entry["adsurl"]))
        return (("VERIFIED", "ADS scan", rec["bibcode"]) if rec
                else ("UNRESOLVED", "ADS scan", f"no scan served for {entry['adsurl']}"))
    return "UNRESOLVED", "none", "entry carries no doi, eprint, adsurl or registered url"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bib", type=Path, default=DEFAULT_BIB)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--refresh", action="store_true", help="ignore cached registry answers")
    args = ap.parse_args(argv)

    entries = parse_bib(args.bib.read_text(encoding="utf-8"))
    keys = [e["key"] for e in entries]
    dup = sorted({k for k in keys if keys.count(k) > 1})
    cache = json.loads(CACHE_PATH.read_text()) if CACHE_PATH.exists() else {}

    rows, counts = [], {}
    for e in entries:
        status, method, detail = verify(e, cache, args.refresh)
        counts[status] = counts.get(status, 0) + 1
        rows.append((e["key"], e["type"], e.get("year", ""), status, method, detail))
        print(f"{status:10s} {e['key']:22s} {method:14s} {detail[:110]}")
    # Failed lookups are not cached: a gateway timeout today must not read as "unresolved" tomorrow.
    CACHE_PATH.write_text(json.dumps({k: v for k, v in cache.items() if v is not None}, indent=1, sort_keys=True))

    lines = [f"# Bibliography verification -- `{args.bib.relative_to(REPO)}`", "",
             f"Emitted by `committee/work/cv-literature/verify_bib.py` on "
             f"{_dt.datetime.now(_dt.timezone.utc):%Y-%m-%d %H:%M} UTC. Do not edit by hand.", "",
             f"**{len(entries)} entries**: " + ", ".join(f"{n} {s}" for s, n in sorted(counts.items()))
             + (f". **Duplicate keys: {dup}**" if dup else ". No duplicate keys."), "",
             "Method: DOI -> Crossref/DataCite record must match first author, year (+-1), volume, first "
             "page/article number and title; eprint -> arXiv record must match first author and title; "
             "ADS scan -> the bibcode's scan must be served as a PDF (a non-existent bibcode returns 403); "
             "publisher page -> page must contain the stated title; MANUAL -> evidence stated.", "",
             "| key | type | year | status | method | evidence / disagreement |", "|---|---|---|---|---|---|"]
    lines += [f"| `{k}` | {t} | {y} | {s} | {m} | {d.replace('|', '/')} |" for k, t, y, s, m, d in rows]
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n{len(entries)} entries: {counts}; duplicates: {dup or 'none'}; table -> {args.out}")
    return 1 if (dup or counts.get("MISMATCH") or counts.get("UNRESOLVED")) else 0


if __name__ == "__main__":
    sys.exit(main())
