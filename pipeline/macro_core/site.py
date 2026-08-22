"""The site assembler — three navigation layers over pages that already exist.

WHY THIS MODULE EXISTS
----------------------
Every page on this site was already generated from data, every number on it
already had a query behind it, and every page but one was already linked from
somewhere.  The site was still unreadable, and the reviewer who was asked to
circulate it said why:

    "frankly quite intimidating, and quite a load of stuff to crawl through
    ... random links that take you places with no clear bread crumbs or
    provenance or going forward or backward or any relation to anything."

That is not a connectivity defect.  Measured before this module was written,
exactly one page (``pipeline/pipeline_status.html``) was unreachable and two
had no link back.  The defect is that **a reader has no sense of place**: no
breadcrumb, no project identity, no position in an argument, no way to tell
whether the page in front of them is the whole story or one panel of it.

So this module adds no facts.  It adds *position*, in three layers:

  LAYER 1  a project row on every page — the shared pipeline plus the six
           projects, each carrying its own live condition, so the row is a
           status board and not just a switch.
  LAYER 2  a view row inside each area — The Case, Figures, Draft Paper,
           Plan & Status, Evidence Detail.  Only views that EXIST are linked;
           a view that does not exist yet says so, in place, rather than
           linking a stub.
  LAYER 3  a sticky left rail carrying THE QUESTIONS of the current view, in
           order, as jump links, with the answers in the main pane.

plus a breadcrumb and previous/next links on every page.

WHERE EVERYTHING COMES FROM (nothing here is hand-maintained)
-------------------------------------------------------------
* the AREAS and their order          ``project_plan.PROJECTS`` (+ the shared
                                     pipeline, which every project stands on)
* each project's condition           ``project_plan.progress_fraction``
* the pipeline's condition           ``provenance`` freshness of its stages
* which page belongs to which area   the page's own path, and
                                     ``project_plan.stage_report_path`` for
                                     the pipeline stages a project's tasks
                                     declare a dependence on
* the order of pages inside an area  ``provenance.topological_order`` — the
                                     build order of the DAG, so the evidence
                                     chain is read in the order it was made
* THE QUESTIONS                      harvested out of the generated evidence
                                     pages' own Socratic sections (the
                                     ``Question / Evidence / Decision /
                                     Consequence`` blocks the report
                                     renderers already emit), plus every
                                     not-yet-done task in the plan ledger,
                                     which is a question nobody has answered
* the deciding number per question   the ``<span class="tag">`` the section
                                     head already carries
* the figures and their captions     the ``<figure>`` blocks of the pages
                                     themselves, plus the ``p5_figure`` table
                                     for the manuscript figures, which belong
                                     to no page
* whether a draft paper exists       a compiled ``manuscripts/<key>/main.pdf``
* the draft paper's shape            its own ``\\section{}`` headings

THE ONE HAND-STATED TABLE is :data:`SHORT_LABEL` — the two-or-three-word name
each project wears in the top row.  A project's ``title`` is a sentence
("T CrB Pre-Eruption Monitoring"); a navigation row needs "T CrB".  Every
rule tried for deriving one from the other failed on at least one of the six.
:func:`validate` therefore requires an entry for every project key, so adding
a project fails loudly here instead of drifting silently on the page.

HOW EXISTING PAGES KEEP WORKING
-------------------------------
The thirteen report renderers are not touched and not forked.  Each writes
its page exactly as before; :func:`wrap_page` then re-reads that page and
puts the chrome around its body, between marker comments.  Re-wrapping an
already-wrapped page recovers the original body from between those markers
first, so the operation is idempotent and a renderer may re-run at any time.
``update_project_plan.py render`` calls :func:`build_site` after it renders,
so the normal working rhythm rebuilds the chrome without anybody remembering.
"""

from __future__ import annotations

import html as _html
import re
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

from . import project_plan as pp
from . import provenance as pv

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"

SITE_CODE_VERSION = "SITE v1.0 (2026-08-21)"

#: The area key of the shared pipeline — the one area that is not a project.
PIPELINE = "pipeline"

#: Its title.  Named the way the row names it.
PIPELINE_TITLE = "Data & Instrument Characterization"

#: THE ONE HAND-STATED TABLE.  See the module docstring.  Guarded by
#: :func:`validate`, which every build runs first.
SHORT_LABEL: dict[str, str] = {
    "TCrB_Monitoring":        "T CrB",
    "CV_TimeSeries":          "CV Time Series",
    "SN2023ixf_LightCurve":   "SN 2023ixf",
    "BeStar_Grism":           "Be-Star",
    "DwarfGalaxy_AGN_Survey": "Dwarf Galaxies",
    "Legacy_Rigel":           "Legacy Rigel",
}

#: Markers that fence the original page body inside a wrapped page, so the
#: wrap can be undone and re-applied without ever eating the content.
CONTENT_BEGIN = "<!--MACRO-SITE:CONTENT-->"
CONTENT_END = "<!--/MACRO-SITE:CONTENT-->"

#: How long a rail entry may be.  The rail is 258 px wide, which is about
#: 34 characters a line, so this is roughly two lines and a bit — long enough
#: to carry the sense of a question, short enough that twenty of them are one
#: glance rather than one page.  The full text is always in the pane.
RAIL_CHARS = 74

#: Files under ``docs/`` this module writes itself.  They are never treated
#: as harvestable evidence pages — a Case page listing its own questions as
#: evidence would be a hall of mirrors.
GENERATED_NAMES = frozenset({"case.html", "figures.html", "paper.html",
                             "evidence.html"})

#: The two views that exist only while a condition holds: a Figures wall
#: needs figures, a Draft Paper page needs prose.  ``provenance``'s WEB stage
#: declares every page this module ALWAYS writes, and deliberately not these
#: two — a declared output that is legitimately absent reads as
#: OUTPUT_MISSING forever, which would turn the one stage covering the whole
#: public site into a permanent red light.  Their gate is
#: ``build_site.py --check``, which rebuilds and compares rather than
#: fingerprinting, and is the stronger of the two checks anyway.
CONDITIONAL_VIEWS = frozenset({"figures.html", "paper.html"})

#: Pages a renderer other than this one owns but that are still ours to wrap.
#: ``docs/index.html`` is the landing (written here); ``docs/evidence.html``
#: is the full evidence index (written by ``report_projects``).
LANDING = "index.html"
EVIDENCE_INDEX = "evidence.html"


class SiteError(RuntimeError):
    """Raised when the site cannot be assembled honestly."""


# ===========================================================================
# 1.  SMALL HTML HELPERS
# ===========================================================================

def esc(s) -> str:
    """Escape text for HTML.  Same rule as the report renderers'."""
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


_TAG_RE = re.compile(r"<[^>]+>")
_WS_RE = re.compile(r"\s+")


def strip_tags(fragment: str) -> str:
    """Plain text of an HTML fragment, entities resolved and space collapsed.

    Used only for things that must be TEXT — a rail label, a ``<title>``, a
    breadcrumb.  Never for the main pane, which keeps the renderer's markup.
    """
    return _WS_RE.sub(" ", _html.unescape(_TAG_RE.sub(" ", fragment))).strip()


def clip(text: str, limit: int) -> str:
    """Trim to ``limit`` characters on a word boundary, with an ellipsis.

    The rail has to be scannable in fifteen seconds, which means its entries
    have to fit on one or two lines.  The full text is always one click away
    in the pane, so clipping the label costs a reader nothing.
    """
    text = text.strip()
    if len(text) <= limit:
        return text
    cut = text[:limit].rsplit(" ", 1)[0].rstrip(" ,;:.—-")
    return f"{cut}…"


_DASHES = ("—", "–", "-")


def within(area_title: str, page_title: str) -> str:
    """A page's title with the area's own name taken off the front.

    Every CV report is titled "Cataclysmic-Variable Time Series — the
    catalogue tie", which is right on the page itself and wrong six times
    over in a list on the CV project's own Figures view, where the reader
    already knows which project they are in. Stripping it is what turns a
    column of near-identical strings into six distinguishable names.
    """
    for dash in _DASHES:
        prefix = f"{area_title} {dash} "
        if page_title.startswith(prefix):
            rest = page_title[len(prefix):]
            return rest[0].upper() + rest[1:] if rest else page_title
    return page_title


def rel_href(from_rel: str, to_rel: str) -> str:
    """A docs-relative link from one docs-relative page to another.

    Both arguments are paths under ``docs/`` with forward slashes.  Computed
    rather than typed because the same chrome is emitted onto pages at two
    different depths, and a hand-written ``../`` is the classic way a
    generated site develops dead links only the deep pages show.
    """
    up = [".."] * (from_rel.count("/"))
    return "/".join(up + [to_rel]) if up else to_rel


# ===========================================================================
# 2.  WHAT WE HARVEST
# ===========================================================================

@dataclass(frozen=True)
class Question:
    """One Socratic unit, lifted out of a page that already argued it.

    ``ask`` is the question the section poses, ``deciding`` the number the
    section head already displays as its verdict, ``decision`` the callout
    the page draws underneath.  All three are the renderer's own markup; this
    module writes none of them.
    """

    page: str            # docs-relative path of the page it lives on
    anchor: str          # id to jump to on that page
    number: int          # 1-based position within the page
    heading: str         # the section's <h2>, as HTML
    ask: str             # the question prose, as HTML ("" when the heading
                         # is itself the question)
    deciding: str        # the section head's tag — the deciding number
    decision: str        # the <div class="decision"> callout, as HTML
    fig_src: str         # docs-relative figure path ("" for none)
    fig_caption: str

    @property
    def label(self) -> str:
        """What the rail shows: the question if there is one, else the
        heading.  A heading like "Pointing validation" tells a reader less
        than "Do the coordinates in the header point where the telescope was
        looking?" — and the pages carry both."""
        text = strip_tags(self.ask) or strip_tags(self.heading)
        return clip(_drop_leading_number(text), RAIL_CHARS)

    @property
    def title(self) -> str:
        return strip_tags(self.heading)


@dataclass(frozen=True)
class Figure:
    """One generated plot, with the caption its own page gave it."""

    src: str             # docs-relative
    caption: str         # HTML
    page: str            # docs-relative page it belongs to ("" if none)
    anchor: str
    category: str        # the figure directory, e.g. "s3" or "cv_paper"


@dataclass(frozen=True)
class PageDoc:
    """A harvested evidence page."""

    rel: str
    title: str
    lead: str
    questions: tuple[Question, ...]
    figures: tuple[Figure, ...]


_NUM_PREFIX_RE = re.compile(r"^\s*\d+(?:\.\d+)?\s*[·.\u00b7:—-]?\s+")


def _drop_leading_number(text: str) -> str:
    """"3 · Era structure" -> "Era structure".

    The numbers are the page's own ordering, and the rail supplies its own;
    printing both gives "3 3 · Era structure", which is how a generated
    navigation starts to look machine-made.
    """
    return _NUM_PREFIX_RE.sub("", text)


# --- the parser -----------------------------------------------------------
# Deliberately regex over the generated markup rather than an HTML parser.
# These pages are emitted by thirteen renderers that all share one house
# shape (`<section id>` / `.bhead` / `<h3>Question</h3>` / `.decision`), and
# the shape IS the contract this module depends on.  A tolerant parser would
# quietly return nothing when a renderer drifted off the shape; these
# patterns return nothing loudly, and `test_site.py` asserts a minimum
# harvest per page so the drift fails a test instead of emptying a rail.

_BODY_RE = re.compile(r"<body[^>]*>(.*)</body>", re.S | re.I)
_H1_RE = re.compile(r"<h1[^>]*>(.*?)</h1>", re.S | re.I)
_HEADER_P_RE = re.compile(r"<header[^>]*>.*?<p[^>]*>(.*?)</p>", re.S | re.I)
_SUB_P_RE = re.compile(r"<p class=['\"]sub['\"][^>]*>(.*?)</p>", re.S | re.I)
_SECTION_SPLIT_RE = re.compile(r"(?=<section\b)", re.I)
_SECTION_ID_RE = re.compile(r"<section\b[^>]*\bid=['\"]([^'\"]+)['\"]", re.I)
_H2_RE = re.compile(r"<h2[^>]*>(.*?)</h2>", re.S | re.I)
_TAG_SPAN_RE = re.compile(r"<span class=['\"]tag['\"][^>]*>(.*?)</span>",
                          re.S | re.I)
_QUESTION_RE = re.compile(
    r"<h3[^>]*>\s*Question\s*</h3>\s*<p[^>]*>(.*?)</p>", re.S | re.I)
_DECISION_RE = re.compile(
    r"<div class=['\"]decision['\"][^>]*>(.*?)</div>", re.S | re.I)
_IMG_RE = re.compile(r"<img[^>]*\bsrc=['\"]([^'\"]+)['\"][^>]*>", re.I)
_FIGURE_RE = re.compile(
    r"<figure\b[^>]*>(.*?)</figure>", re.S | re.I)
_FIGCAPTION_RE = re.compile(r"<figcaption[^>]*>(.*?)</figcaption>", re.S | re.I)
_ANCHOR_SPAN = '<span class="qanchor" id="{}"></span>'
_ANCHOR_SPAN_RE = re.compile(
    r"<span class=['\"]qanchor['\"] id=['\"]([^'\"]+)['\"]></span>", re.I)


def page_body(html: str) -> str:
    """The inside of ``<body>``, or the whole document if it has no body."""
    match = _BODY_RE.search(html)
    return match.group(1) if match else html


def original_content(body: str) -> str:
    """The renderer's own markup, recovered from a wrapped page.

    This is what makes :func:`wrap_page` idempotent.  Without it, a second
    build would nest chrome inside chrome and the third would nest it again.
    """
    if CONTENT_BEGIN in body and CONTENT_END in body:
        return body.split(CONTENT_BEGIN, 1)[1].rsplit(CONTENT_END, 1)[0]
    return body


def assign_anchors(content: str) -> str:
    """Give every section a stable id, inventing one only where none exists.

    Two pages (``cv_external_context`` and ``pipeline_status``) emit sections
    with no id at all, so their headings could not be linked to — which is
    precisely the "random links that take you places" problem seen from the
    other end: there was nowhere to land.  Ids are assigned by POSITION, so
    re-running the underlying renderer and re-assigning produces the same
    ids; ids the renderer supplies itself are never touched.
    """
    if "<section" in content.lower():
        out, n = [], 0
        for chunk in _SECTION_SPLIT_RE.split(content):
            if chunk.lstrip().lower().startswith("<section"):
                n += 1
                if not _SECTION_ID_RE.search(chunk):
                    chunk = re.sub(r"<section\b", f'<section id="sec-{n}"',
                                   chunk, count=1, flags=re.I)
            out.append(chunk)
        return "".join(out)
    # No sections at all: anchor the h2s in place.
    if _ANCHOR_SPAN_RE.search(content):
        return content                      # already anchored; idempotent
    counter = {"n": 0}

    def _anchor(match: re.Match) -> str:
        counter["n"] += 1
        return _ANCHOR_SPAN.format(f"sec-{counter['n']}") + match.group(0)

    return re.sub(r"<h2\b[^>]*>", _anchor, content, flags=re.I)


def _section_chunks(content: str) -> list[tuple[str, str]]:
    """``(anchor, html)`` per section, for either page shape."""
    if "<section" in content.lower():
        chunks = []
        for chunk in _SECTION_SPLIT_RE.split(content):
            if not chunk.lstrip().lower().startswith("<section"):
                continue
            match = _SECTION_ID_RE.search(chunk)
            chunks.append((match.group(1) if match else "", chunk))
        return chunks
    # h2-anchored shape: split on the anchor spans assign_anchors inserted.
    parts = _ANCHOR_SPAN_RE.split(content)
    chunks = []
    for i in range(1, len(parts) - 1, 2):
        chunks.append((parts[i], parts[i + 1]))
    return chunks


def harvest_page(path: Path, rel: str) -> PageDoc:
    """Read one generated page and lift its questions and figures out of it.

    The page is *not* modified here.  :func:`wrap_page` writes the anchored
    content back; this function anchors a copy so the harvest and the wrap
    agree on every id.
    """
    raw = path.read_text(encoding="utf-8")
    content = assign_anchors(original_content(page_body(raw)))
    h1 = _H1_RE.search(content)
    title = strip_tags(h1.group(1)) if h1 else rel
    lead_match = _HEADER_P_RE.search(content) or _SUB_P_RE.search(content)
    lead = strip_tags(lead_match.group(1)) if lead_match else ""

    questions: list[Question] = []
    figures: list[Figure] = []
    for number, (anchor, chunk) in enumerate(_section_chunks(content), start=1):
        h2 = _H2_RE.search(chunk)
        tag = _TAG_SPAN_RE.search(chunk.split("</div>", 1)[0]) \
            if 'class="bhead"' in chunk else None
        if tag is None:
            tag = _TAG_SPAN_RE.search(chunk)
        ask = _QUESTION_RE.search(chunk)
        decision = _DECISION_RE.search(chunk)
        fig_src, fig_cap = _first_figure(chunk, rel)
        if h2 is None and ask is None:
            continue                        # not a Socratic unit; skip it
        questions.append(Question(
            page=rel, anchor=anchor, number=number,
            heading=(h2.group(1).strip() if h2 else ""),
            ask=(ask.group(1).strip() if ask else ""),
            deciding=(tag.group(1).strip() if tag else ""),
            decision=(decision.group(1).strip() if decision else ""),
            fig_src=fig_src, fig_caption=fig_cap))
        for src, caption in _all_figures(chunk, rel):
            figures.append(Figure(src=src, caption=caption, page=rel,
                                  anchor=anchor,
                                  category=_category_of(src)))
    # Figures outside any section still belong to the page.
    seen = {f.src for f in figures}
    for src, caption in _all_figures(content, rel):
        if src not in seen:
            figures.append(Figure(src=src, caption=caption, page=rel,
                                  anchor="", category=_category_of(src)))
            seen.add(src)
    return PageDoc(rel=rel, title=title, lead=lead,
                   questions=tuple(questions), figures=tuple(figures))


def _resolve(rel_page: str, src: str) -> str:
    """A page-relative ``<img src>`` turned into a docs-relative path."""
    base = Path(rel_page).parent
    return str((base / src).as_posix()).replace("/./", "/")


def _normalise(path: str) -> str:
    parts: list[str] = []
    for piece in path.split("/"):
        if piece in ("", "."):
            continue
        if piece == ".." and parts:
            parts.pop()
        else:
            parts.append(piece)
    return "/".join(parts)


def _first_figure(chunk: str, rel_page: str) -> tuple[str, str]:
    figs = _all_figures(chunk, rel_page)
    return figs[0] if figs else ("", "")


def _all_figures(chunk: str, rel_page: str) -> list[tuple[str, str]]:
    out = []
    for block in _FIGURE_RE.findall(chunk):
        img = _IMG_RE.search(block)
        if not img:
            continue
        caption = _FIGCAPTION_RE.search(block)
        out.append((_normalise(_resolve(rel_page, img.group(1))),
                    caption.group(1).strip() if caption else ""))
    return out


def _category_of(src: str) -> str:
    """The figure directory, which is the category the pages already use."""
    parts = src.split("/")
    return parts[-2] if len(parts) >= 2 else "figures"


# ===========================================================================
# 3.  THE AREAS — layer 1
# ===========================================================================

@dataclass(frozen=True)
class Area:
    """One entry in the top row: the shared pipeline, or one project."""

    key: str
    label: str                  # what the row shows
    title: str                  # the full name
    is_pipeline: bool
    project: Optional[pp.Project] = None

    @property
    def dir(self) -> str:
        return PIPELINE if self.is_pipeline else self.key

    @property
    def home(self) -> str:
        """The page the top row points at — always the Case."""
        return (f"{PIPELINE}/index.html" if self.is_pipeline
                else f"{self.key}/case.html")


def areas() -> tuple[Area, ...]:
    """The top row, in order: the shared pipeline, then the projects.

    The project order is ``PROJECTS``' own order, which is the order the plan
    ledger declares them in.  Nothing re-sorts it here — a navigation row
    with its own opinion about ordering is a second source of truth.
    """
    out = [Area(key=PIPELINE, label="Data & Instrument",
                title=PIPELINE_TITLE, is_pipeline=True)]
    for project in pp.PROJECTS:
        out.append(Area(key=project.key,
                        label=SHORT_LABEL[project.key],
                        title=project.title, is_pipeline=False,
                        project=project))
    return tuple(out)


def validate() -> None:
    """Every project must have a row label.  Run before every build."""
    missing = [p.key for p in pp.PROJECTS if p.key not in SHORT_LABEL]
    if missing:
        raise SiteError(
            f"no SHORT_LABEL for {', '.join(missing)} — the top navigation "
            f"row cannot name a project it has no short name for. Add one to "
            f"macro_core.site.SHORT_LABEL.")
    extra = [k for k in SHORT_LABEL if k not in pp.PROJECT_BY_KEY]
    if extra:
        raise SiteError(
            f"SHORT_LABEL names {', '.join(extra)}, which the plan ledger "
            f"does not declare — a row entry with no project behind it.")


# ===========================================================================
# 4.  WHICH PAGES BELONG WHERE
# ===========================================================================

def evidence_pages(docs_dir: Path) -> dict[str, list[str]]:
    """``area key -> docs-relative evidence pages``, in DAG build order.

    Membership comes from the path (``docs/pipeline/*`` is the shared
    pipeline; ``docs/<ProjectKey>/*`` is that project).  ORDER comes from
    ``provenance.topological_order``, so the chain reads in the order it was
    built rather than the order the filesystem happens to return.
    """
    order: dict[str, int] = {}
    for i, stage_key in enumerate(pv.topological_order(pv.STAGES)):
        rel = pp.stage_report_path(stage_key)
        if rel and rel.startswith("docs/"):
            order.setdefault(rel[len("docs/"):], i)

    found: dict[str, list[str]] = {}
    for area in areas():
        directory = docs_dir / area.dir
        if not directory.is_dir():
            found[area.key] = []
            continue
        pages = []
        for path in sorted(directory.glob("*.html")):
            name = path.name
            if name in GENERATED_NAMES:
                continue
            if name == "index.html" and not area.is_pipeline:
                continue        # the project's Plan & Status tab, not evidence
            if area.is_pipeline and name == "index.html":
                continue        # the pipeline's Case, written here
            pages.append(f"{area.dir}/{name}")
        pages.sort(key=lambda rel: (order.get(rel, 10_000), rel))
        found[area.key] = pages
    return found


def inherited_pages(project: pp.Project) -> list[str]:
    """The shared-pipeline pages a project's own tasks rest on, in DAG order.

    Derived from ``Task.stage``: the edge that already makes staleness
    propagate into the plan is the same edge that says which pipeline
    evidence this project depends on.  No second list.
    """
    order = {k: i for i, k in enumerate(pv.topological_order(pv.STAGES))}
    wanted: dict[str, int] = {}
    for task in project.tasks:
        rel = pp.stage_report_path(task.stage)
        if rel and rel.startswith("docs/" + PIPELINE + "/"):
            wanted.setdefault(rel[len("docs/"):], order.get(task.stage, 9999))
    return [rel for rel, _ in sorted(wanted.items(), key=lambda kv: kv[1])]


# ===========================================================================
# 5.  THE CHROME
# ===========================================================================

@dataclass
class Condition:
    """What an area's top-row entry says about itself right now."""

    text: str
    tone: str = ""          # a status class, or "" for plain


@dataclass
class Chrome:
    """Everything the three layers need for one page."""

    area: Area
    view: str                       # which layer-2 tab is current
    page_rel: str                   # docs-relative path of the page
    crumbs: Sequence[tuple[str, str]]     # (label, href) — last has href ""
    rail: str                       # rendered rail HTML
    conditions: Mapping[str, Condition]
    tabs: Sequence["Tab"]
    prev: Optional[tuple[str, str]] = None      # (label, href)
    next: Optional[tuple[str, str]] = None


@dataclass(frozen=True)
class Tab:
    label: str
    rel: str                # docs-relative target, "" when the view is absent
    why_absent: str = ""


def render_topbar(chrome: Chrome) -> str:
    items = [
        f'<a class="brand" href="{rel_href(chrome.page_rel, LANDING)}">'
        f'<span class="dot"></span>MACRO</a>']
    for area in areas():
        current = ' aria-current="page"' if area.key == chrome.area.key else ""
        shared = " shared" if area.is_pipeline else ""
        cond = chrome.conditions.get(area.key, Condition(""))
        items.append(
            f'<a class="navitem{shared}"{current} '
            f'href="{rel_href(chrome.page_rel, area.home)}" '
            f'title="{esc(area.title)}">'
            f'<span class="n">{esc(area.label)}</span>'
            f'<span class="m">{esc(cond.text)}</span></a>')
    return ('<div class="topbar"><nav class="topbar-inner" '
            'aria-label="Projects">' + "".join(items) + "</nav></div>")


def render_tabbar(chrome: Chrome) -> str:
    out = []
    for tab in chrome.tabs:
        if not tab.rel:
            out.append(f'<span class="tab absent" title="{esc(tab.why_absent)}">'
                       f'{esc(tab.label)} — none yet</span>')
            continue
        current = ' aria-current="page"' if tab.label == chrome.view else ""
        out.append(f'<a class="tab"{current} '
                   f'href="{rel_href(chrome.page_rel, tab.rel)}">'
                   f'{esc(tab.label)}</a>')
    return ('<div class="tabbar"><nav class="tabbar-inner" '
            f'aria-label="{esc(chrome.area.title)} views">'
            + "".join(out) + "</nav></div>")


def render_crumbs(chrome: Chrome) -> str:
    items = []
    for label, href in chrome.crumbs:
        if href:
            items.append(f'<li><a href="{rel_href(chrome.page_rel, href)}">'
                         f'{esc(label)}</a></li>')
        else:
            items.append(f'<li aria-current="page">{esc(label)}</li>')
    return ('<div class="crumbs"><nav aria-label="Breadcrumb"><ol>'
            + "".join(items) + "</ol></nav></div>")


def render_prevnext(chrome: Chrome) -> str:
    if not chrome.prev and not chrome.next:
        return ""
    parts = []
    if chrome.prev:
        label, href = chrome.prev
        parts.append(f'<a class="pv" href="{rel_href(chrome.page_rel, href)}">'
                     f'<span class="dir">&larr; Previous</span>{esc(label)}</a>')
    if chrome.next:
        label, href = chrome.next
        parts.append(f'<a class="nx" href="{rel_href(chrome.page_rel, href)}">'
                     f'<span class="dir">Next &rarr;</span>{esc(label)}</a>')
    return f'<div class="prevnext">{"".join(parts)}</div>'


def render_rail(groups: Sequence[tuple[str, Sequence[tuple[str, str, str]]]],
                back: Optional[tuple[str, str]] = None,
                note: str = "", page_rel: str = "") -> str:
    """The sticky left rail.

    ``groups`` is ``(group heading, [(number, label, href)])``.  The heading
    may be ``""`` for an ungrouped list.  A rail with one entry is not worth
    the column, so an empty rail renders as an empty string and the layout
    collapses to a single column.
    """
    if not any(entries for _, entries in groups):
        return ""
    out = ['<aside class="rail" aria-label="Questions in this view">']
    if back:
        label, href = back
        out.append(f'<a class="railback" href="{rel_href(page_rel, href)}">'
                   f'&larr; {esc(label)}</a>')
    out.append("<h2>The questions</h2>")
    if note:
        out.append(f'<p class="railnote">{note}</p>')
    for heading, entries in groups:
        if not entries:
            continue
        if heading:
            out.append(f'<div class="railgroup">{esc(heading)}</div>')
        out.append("<ol>")
        for number, label, href in entries:
            out.append(f'<li><a href="{href}">'
                       f'<span class="qn">{esc(number)}</span> {esc(label)}'
                       f"</a></li>")
        out.append("</ol>")
    out.append("</aside>")
    return "".join(out)


def shell(chrome: Chrome, title: str, body: str, footer: str,
          head_extra: str = "") -> str:
    """One complete page: head, three layers, pane, footer."""
    css = rel_href(chrome.page_rel, "assets/macro.css")
    return (
        f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">\n'
        f'<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        f"<title>{esc(title)}</title>\n"
        f'<link rel="stylesheet" href="{css}">\n{head_extra}'
        f"</head><body>\n"
        f"{render_topbar(chrome)}\n{render_tabbar(chrome)}\n"
        f"{render_crumbs(chrome)}\n"
        f'<div class="layout">\n{chrome.rail}\n<div class="pane">\n'
        f"{body}\n{render_prevnext(chrome)}\n</div>\n</div>\n"
        f"{footer}\n</body></html>\n")


# ===========================================================================
# 6.  CONDITIONS — what each top-row entry says about itself
# ===========================================================================

def project_condition(project: pp.Project,
                      statuses: Mapping[str, str]) -> Condition:
    counts = pp.status_counts(project.tasks, statuses)
    done, total = pp.progress_fraction(counts)
    tone = "done" if done == total else (
        "blocked" if counts.get(pp.BLOCKED) and not counts.get(pp.IN_PROGRESS)
        else "in_progress")
    return Condition(f"{done}/{total} done", tone)


def analysis_stage_of(page_rel: str) -> str:
    """The stage whose ANALYSIS a report page publishes, or "".

    Two stages point at every report: the analysis (``S2``) and the review
    that renders it (``R-S2``).  Which one a view should quote depends on
    what the view is claiming.

    Everywhere on this site, the claim is about the NUMBERS — "can I trust
    what this page says?" — and that is the analysis stage.  The difference
    is not cosmetic: S2's detector tables were wiped, so ``S2`` reads
    DESTROYED while ``R-S2`` reads merely STALE, and a page that showed the
    renderer's verdict beside S2's High Gain ceiling would be telling a
    reader the number was out of date when the truth is that nothing backs
    it at all.

    Preferring the non-``R-`` key is not a naming convention doing real
    work by accident: ``project_plan.stage_report_path`` is the same map the
    Plan pages read, and this picks the same side of it they do.
    """
    keys = [s.key for s in pv.STAGES
            if pp.stage_report_path(s.key) == f"docs/{page_rel}"]
    analysis = [k for k in keys if not k.startswith("R-")]
    if analysis:
        return analysis[0]
    return keys[0] if keys else ""


def pipeline_stages(docs_dir: Path) -> list[tuple[str, str]]:
    """``(stage key, page)`` for each page the shared pipeline publishes.

    ONE stage per page.  Counting through ``stage_report_path`` directly
    would double-count: every pipeline stage has a paired ``R-*`` review
    stage naming the same report, which made the top row announce "0/17
    stages fresh" over ten pages.  A denominator a reader cannot count on
    the page is a number that teaches them to distrust the rest.
    """
    out = []
    directory = docs_dir / PIPELINE
    for path in sorted(directory.glob("*.html")) if directory.is_dir() else []:
        if path.name in GENERATED_NAMES or path.name == "index.html":
            continue
        stage = analysis_stage_of(f"{PIPELINE}/{path.name}")
        if stage:
            out.append((stage, f"{PIPELINE}/{path.name}"))
    return out


def pipeline_condition(freshness: Mapping[str, pv.Freshness],
                       stages: Sequence[tuple[str, str]]) -> Condition:
    """How much of the shared pipeline is currently trustworthy.

    Only the stages that publish a page under ``docs/pipeline/`` are counted:
    those are the ones this area actually shows, and a fraction that counted
    stages a reader cannot open would be a number with nothing behind it.
    """
    if not stages:
        return Condition("the shared foundation")
    fresh = sum(1 for key, _ in stages
                if freshness.get(key) and freshness[key].state == pv.FRESH)
    tone = "done" if fresh == len(stages) else "redo_needed"
    return Condition(f"{fresh}/{len(stages)} stages fresh", tone)


# ===========================================================================
# 7.  THE BUILDER
# ===========================================================================

class Build:
    """One whole site build.  Everything expensive happens once, here."""

    def __init__(self, docs_dir: Path, repo_root: Path,
                 manifest: Optional[Path] = None):
        validate()
        self.docs = docs_dir
        self.repo = repo_root
        self.generated = pp.utcnow()
        self.areas = areas()
        self.area_by_key = {a.key: a for a in self.areas}

        self.freshness: dict[str, pv.Freshness] = {}
        self.fingerprints: dict[str, str] = {}
        self.statuses: dict[str, str] = {}
        self.ever_ran: set[str] = set()
        if manifest and Path(manifest).exists():
            con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True,
                                  timeout=300.0)
            con.execute("PRAGMA busy_timeout = 300000")
            try:
                self.freshness, self.fingerprints = pp.stage_freshness(
                    con, repo_root)
                self.statuses = pp.read_statuses(con)
                self.ever_ran = pp.stages_ever_run(con)
            except sqlite3.Error:
                # A build must never be blocked by a locked archive drive.
                # It degrades to the ledger's own declared statuses, and the
                # footer says which of the two the page is showing.
                self.statuses = {}
            finally:
                con.close()
        self.live = bool(self.freshness)

        self.pages: dict[str, PageDoc] = {}
        self.by_area = evidence_pages(docs_dir)
        for rels in self.by_area.values():
            for rel in rels:
                path = docs_dir / rel
                if path.exists():
                    self.pages[rel] = harvest_page(path, rel)

        self.pipeline_stages = pipeline_stages(docs_dir)
        self.conditions: dict[str, Condition] = {}
        for area in self.areas:
            if area.is_pipeline:
                self.conditions[area.key] = pipeline_condition(
                    self.freshness, self.pipeline_stages)
            else:
                overlay = pp.overlay_statuses(area.project.tasks,
                                             self.statuses)
                self.conditions[area.key] = project_condition(area.project,
                                                              overlay)
        self.paper_figs = _paper_figures(repo_root, docs_dir)

    # -- per-area facts ---------------------------------------------------
    def questions_of(self, area: Area) -> list[Question]:
        out: list[Question] = []
        for rel in self.by_area.get(area.key, []):
            doc = self.pages.get(rel)
            if doc:
                out.extend(doc.questions)
        return out

    def inherited_questions(self, area: Area) -> list[tuple[str, list[Question]]]:
        if area.is_pipeline or area.project is None:
            return []
        out = []
        for rel in inherited_pages(area.project):
            doc = self.pages.get(rel)
            if doc and doc.questions:
                out.append((rel, list(doc.questions)))
        return out

    def open_tasks(self, area: Area) -> list[pp.Task]:
        if area.project is None:
            return []
        overlay = pp.overlay_statuses(area.project.tasks, self.statuses)
        return [t for t in area.project.tasks if overlay[t.id] != pp.DONE]

    def status_of(self, task: pp.Task) -> str:
        return pp.overlay_statuses([task], self.statuses)[task.id]

    def figures_of(self, area: Area) -> list[Figure]:
        seen, out = set(), []
        for rel in self.by_area.get(area.key, []):
            doc = self.pages.get(rel)
            if not doc:
                continue
            for fig in doc.figures:
                if fig.src not in seen:
                    seen.add(fig.src)
                    out.append(fig)
        for fig in self.paper_figs:
            if fig.src.startswith(area.dir + "/") and fig.src not in seen:
                seen.add(fig.src)
                out.append(fig)
        return out

    def paper_pdf(self, area: Area) -> str:
        """Repo-relative path of a real, compiled draft, or "" if none.

        TWO tests, and the second one earns its keep.  A ``main.tex`` on its
        own is not a draft: five of the six projects carry the same 1.2 kB
        AASTeX skeleton.  A compiled ``main.pdf`` is not enough either —
        T CrB has one, produced by running LaTeX over that skeleton, and a
        "Draft Paper" tab opening on a title page followed by six headings
        and six ``% TODO`` comments is precisely the stub-link this rebuild
        exists to remove.

        So the test is whether any section has been WRITTEN — see
        :func:`tex_sections`.  Structure alone does not distinguish them
        (the skeleton declares all six standard headings); prose does.  This
        is measured from the manuscript itself rather than kept as a list of
        "which projects have a paper", which would go stale the first day
        somebody starts writing the second one.
        """
        if area.project is None:
            return ""
        folder = self.repo / "manuscripts" / area.key
        pdf, tex = folder / "main.pdf", folder / "main.tex"
        if not (pdf.exists() and tex.exists()):
            return ""
        sections = tex_sections(tex.read_text(encoding="utf-8",
                                              errors="replace"))
        return (f"manuscripts/{area.key}/main.pdf"
                if any(s.written for s in sections) else "")

    # -- the layer-2 row --------------------------------------------------
    def tabs_for(self, area: Area) -> list[Tab]:
        directory = area.dir
        tabs: list[Tab] = []
        if area.is_pipeline:
            tabs.append(Tab("The Case", f"{directory}/index.html"))
        else:
            tabs.append(Tab("The Case", f"{directory}/case.html"))
        figs = self.figures_of(area)
        tabs.append(Tab("Figures", f"{directory}/figures.html" if figs else "",
                        "no generated plots belong to this area yet"))
        if area.is_pipeline:
            tabs.append(Tab("Freshness & DAG",
                            f"{directory}/pipeline_status.html"))
        else:
            pdf = self.paper_pdf(area)
            tabs.append(Tab("Draft Paper", f"{directory}/paper.html" if pdf
                            else "",
                            "no compiled draft yet — only the AASTeX "
                            "skeleton; see Plan & Status for the writing "
                            "task and its blocker"))
            tabs.append(Tab("Plan & Status", f"{directory}/index.html"))
        tabs.append(Tab("Evidence Detail", f"{directory}/evidence.html"))
        return tabs

    def chrome(self, area: Area, view: str, page_rel: str, rail: str,
               crumbs: Sequence[tuple[str, str]],
               prev=None, nxt=None) -> Chrome:
        return Chrome(area=area, view=view, page_rel=page_rel, crumbs=crumbs,
                      rail=rail, conditions=self.conditions,
                      tabs=self.tabs_for(area), prev=prev, next=nxt)


def _paper_figures(repo_root: Path, docs_dir: Path) -> list[Figure]:
    """The manuscript figures, with the captions the manuscript gives them.

    These thirteen plots belong to no evidence page — they are built straight
    into the paper — so nothing on the site showed them.  ``p5_figure`` holds
    one row per figure with its label, title and full caption, which is the
    same table the LaTeX pastes from, so the Figures view and the paper
    cannot disagree about what a figure says.
    """
    db = repo_root / "products" / "phot" / "cv_timeseries.sqlite"
    if not db.exists():
        return []
    try:
        con = sqlite3.connect(f"file:{db}?mode=ro", uri=True, timeout=30.0)
    except sqlite3.Error:                                # pragma: no cover
        return []
    try:
        rows = con.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name='p5_figure'").fetchall()
        if not rows:
            return []
        out = []
        for fig_id, title, caption, png in con.execute(
                "SELECT fig_id, title, caption, png_path FROM p5_figure "
                "ORDER BY fig_id"):
            if not png:
                continue
            path = Path(png)
            try:
                rel = path.resolve().relative_to(docs_dir.resolve()).as_posix()
            except (ValueError, OSError):
                # png_path is recorded as an absolute build path; fall back
                # to matching the file name inside the CV figure directory.
                candidates = list(docs_dir.rglob(f"cv_paper/{path.name}"))
                if not candidates:
                    continue
                rel = candidates[0].relative_to(docs_dir).as_posix()
            out.append(Figure(
                src=rel,
                caption=f"<b>{esc(title)}</b> — {esc(caption)}",
                page="", anchor="", category="cv_paper"))
        return out
    except sqlite3.Error:                                # pragma: no cover
        return []
    finally:
        con.close()


# ===========================================================================
# 8.  THE VIEWS
# ===========================================================================

def _q_id(index: int) -> str:
    return f"q{index}"


def render_case(build: Build, area: Area) -> Path:
    """THE CENTREPIECE.

    For every question the area has already answered: the question, the
    deciding number, the figure if there is one, the decision it produced,
    and a link to the page carrying the full working.  Then the questions
    the plan has posed and nobody has answered yet.
    """
    rel = area.home
    own = build.questions_of(area)
    inherited = build.inherited_questions(area)
    open_tasks = build.open_tasks(area)

    # The rail is grouped by the analysis step each question came from.  A
    # flat list of 54 entries is not scannable in fifteen seconds — it is the
    # same wall of undifferentiated links the old site was, moved into a
    # narrower column.  Ten headed groups of five are.
    rail_groups: list[tuple[str, list[tuple[str, str, str]]]] = []
    index = 0
    for page_rel in build.by_area.get(area.key, []):
        doc = build.pages.get(page_rel)
        if not doc or not doc.questions:
            continue
        entries = []
        for q in doc.questions:
            entries.append((str(index + 1), q.label, f"#{_q_id(index)}"))
            index += 1
        rail_groups.append((within(area.title, doc.title), entries))

    inh_entries = []
    n = 0
    for page_rel, questions in inherited:
        for q in questions:
            n += 1
            inh_entries.append(
                (f"i{n}", q.label,
                 rel_href(rel, page_rel) + (f"#{q.anchor}" if q.anchor else "")))
    rail_groups.append(("Rests on the shared pipeline", inh_entries))
    rail_groups.append(("Still open", [
        (f"o{i + 1}", clip(t.title, RAIL_CHARS), f"#open-{esc(t.id)}")
        for i, t in enumerate(open_tasks)]))

    rail = render_rail(
        rail_groups,
        note=("Every entry is a question this project had to answer. "
              "Click the one you doubt."),
        page_rel=rel)

    body: list[str] = []
    subject = PIPELINE_TITLE if area.is_pipeline else area.title
    body.append(f"<header><h1>{esc(subject)} — the case</h1>")
    if area.is_pipeline:
        body.append(
            "<p>The instrument and the archive, argued one question at a "
            "time. Every project on this site stands on the answers below; "
            "when one of them stops being true, the projects that rest on it "
            "say so on their own pages.</p></header>")
    else:
        body.append(f"<p>{esc(area.project.claim)}</p></header>")
        body.append(f'<div class="decision"><b>Where it is going:</b> '
                    f"{esc(area.project.venue)}</div>")

    # The lead counts only what is actually here.  "0 inherited from the
    # shared pipeline, 0 still open in the plan" on the pipeline's own page
    # is three numbers where one is true, and two zeroes a reader has to
    # work out are structural rather than a report of nothing done.
    tally = [f"<b>{len(own)}</b> question(s) answered on "
             + ("these pages" if area.is_pipeline
                else "this project's own evidence pages")]
    if inh_entries:
        tally.append(f"<b>{len(inh_entries)}</b> inherited from the shared "
                     f"pipeline")
    if open_tasks:
        tally.append(f"<b>{len(open_tasks)}</b> still open in the plan")
    body.append(f'<p class="lead">{", ".join(tally)}.</p>')

    if own:
        body.append('<h2 id="answered">The argument so far</h2>')
    for i, q in enumerate(own):
        body.append(_question_block(build, area, q, i, rel, own))
    if not own:
        body.append(
            '<div class="absentnote"><b>No evidence page of its own yet.</b> '
            "Everything this project has established so far it established "
            "through the shared pipeline — those questions are listed below, "
            "and the work still to do is at the bottom of this page.</div>")

    if inherited:
        body.append('<h2 id="rests-on">What this rests on</h2>')
        body.append('<p class="lead">These questions were answered once, for '
                    "the whole archive, and this project inherits their "
                    "answers through the tasks that name their stages.</p>")
        for page_rel, questions in inherited:
            doc = build.pages[page_rel]
            # THE VERDICT TRAVELS WITH THE LINK.  A pipeline report whose
            # stage is no longer FRESH still holds the numbers it printed
            # when it ran — S2's detector constants are the standing case —
            # and a project that links one without saying so has quietly
            # re-asserted a measurement whose table no longer exists.  The
            # rule already governed the Plan pages; the Case page is a new
            # way into the same reports, so it carries the same chip.
            body.append(f'<div class="phase"><h3><a href="'
                        f'{rel_href(rel, page_rel)}">'
                        f"{esc(within(area.title, doc.title))}</a> "
                        f"{_page_verdict(build, page_rel)}</h3>")
            for q in questions:
                href = rel_href(rel, page_rel) + (f"#{q.anchor}"
                                                  if q.anchor else "")
                deciding = (f' <span class="tag">{q.deciding}</span>'
                            if q.deciding else "")
                body.append(f'<p class="sub"><a href="{href}">'
                            f"{esc(q.label)}</a>{deciding}</p>")
            body.append("</div>")

    if open_tasks:
        body.append('<h2 id="open">Still open</h2>')
        body.append('<p class="lead">Each of these is a question the '
                    "committee strategy put on the plan and nobody has "
                    "answered yet. The citation is where it came from; the "
                    "blocker, where there is one, is what would clear it."
                    "</p>")
        for task in open_tasks:
            body.append(_open_block(build, area, task, rel))

    footer = _footer(build, area, rel)
    crumbs = [("Home", LANDING), (area.title, area.home), ("The Case", "")]
    chrome = build.chrome(area, "The Case", rel, rail, crumbs)
    path = build.docs / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(shell(chrome, f"MACRO — {subject}: the case",
                          "\n".join(body), footer), encoding="utf-8")
    return path


def _question_block(build: Build, area: Area, q: Question, index: int,
                    rel: str, siblings: Sequence[Question]) -> str:
    doc = build.pages[q.page]
    href = rel_href(rel, q.page) + (f"#{q.anchor}" if q.anchor else "")
    # The heading loses the source page's own section number: this page
    # supplies one, and "Q1  1 · The artifact" is how a generated navigation
    # starts to look machine-made.
    heading = esc(_drop_leading_number(q.title)) if q.title else esc(q.label)
    # THE VERDICT SITS BESIDE THE NUMBER IT QUALIFIES.  This view reprints a
    # deciding number that was measured on another page, and that page's
    # stage may since have been re-run, invalidated, or had its tables wiped
    # outright — S2's detector constants are the standing example, still
    # printed as fact by a report whose backing tables no longer exist.  The
    # rule the project pages already follow is that no such number is shown
    # without its verdict; a new way into the same numbers inherits the rule.
    out = [f'<div class="qblock" id="{_q_id(index)}">',
           f'<div class="qhead"><span class="qnum">Q{index + 1}</span>'
           f"<h2>{heading}</h2>{_page_verdict(build, q.page)}</div>"]
    if q.ask:
        out.append(f'<p class="qask">{q.ask}</p>')
    if q.deciding:
        out.append(f'<div class="answer"><span class="lbl">The deciding '
                   f'number</span><div class="deciding">{q.deciding}</div>'
                   f"</div>")
    if q.fig_src:
        out.append(f'<figure class="qfig"><a href="{rel_href(rel, q.fig_src)}">'
                   f'<img loading="lazy" src="{rel_href(rel, q.fig_src)}" '
                   f'alt="{esc(strip_tags(q.fig_caption) or q.title)}"></a>')
        if q.fig_caption:
            out.append(f"<figcaption>{q.fig_caption}</figcaption>")
        out.append("</figure>")
    if q.decision:
        out.append(f'<div class="answer"><span class="lbl">The decision it '
                   f'produced</span><div class="decision">{q.decision}</div>'
                   f"</div>")
    nav = []
    if index:
        nav.append(f'<a href="#{_q_id(index - 1)}">&larr; previous question</a>')
    if index + 1 < len(siblings):
        nav.append(f'<a href="#{_q_id(index + 1)}">next question &rarr;</a>')
    out.append(f'<div class="qfoot"><span class="full">'
               f'<a href="{href}">Full working: {esc(doc.title)}</a></span>'
               + "".join(f"<span>{link}</span>" for link in nav) + "</div>")
    out.append("</div>")
    return "".join(out)


def _open_block(build: Build, area: Area, task: pp.Task, rel: str) -> str:
    status = build.status_of(task)
    out = [f'<div class="qblock" id="open-{esc(task.id)}">',
           f'<div class="qhead"><span class="qnum">OPEN</span>'
           f"<h2>{esc(task.title)}</h2>"
           f'<span class="chip {esc(status)}">{esc(pp.STATUS_LABEL[status])}'
           f"</span></div>",
           f'<p class="qask">What would settle it: {esc(task.produces)}</p>']
    if task.blocker:
        out.append(f'<div class="blockcard"><b>Blocked:</b> '
                   f"{esc(task.blocker)}</div>")
    if task.forbids:
        out.append(f'<div class="openq"><b>The strategy forbids running this '
                   f"first:</b> {esc(task.forbids)}</div>")
    if task.depends_on:
        out.append(f'<p class="src">Waits on: '
                   f"{esc(', '.join(task.depends_on))}</p>")
    plan_rel = f"{area.dir}/index.html"
    out.append(f'<div class="qfoot"><span class="src">'
               f"{esc(task.phase)} &middot; rests on stage "
               f"<code>{esc(task.stage)}</code> &middot; derived from "
               f"{esc(str(task.source))}</span>"
               f'<span><a href="{rel_href(rel, plan_rel)}#plan">'
               f"See it in the plan</a></span></div>")
    out.append("</div>")
    return "".join(out)


# --- the figures view -----------------------------------------------------

def render_figures(build: Build, area: Area) -> Optional[Path]:
    figs = build.figures_of(area)
    if not figs:
        return None
    rel = f"{area.dir}/figures.html"
    groups: dict[str, list[Figure]] = {}
    for fig in figs:
        groups.setdefault(fig.category, []).append(fig)

    # The category's name is the TITLE OF THE PAGE that owns its figures —
    # derived, and far more useful than the directory slug. The manuscript
    # figures own no page, so they are named for the manuscript.
    def label(category: str, members: Sequence[Figure]) -> tuple[str, str]:
        pages = [m.page for m in members if m.page]
        if pages:
            doc = build.pages.get(pages[0])
            if doc:
                return within(area.title, doc.title), pages[0]
        return "Draft-paper figures (built straight into the manuscript)", ""

    ordered = []
    page_order = {rel_: i for i, rel_ in
                  enumerate(build.by_area.get(area.key, []))}
    for category, members in groups.items():
        name, page = label(category, members)
        ordered.append((page_order.get(page, 9_999), category, name, page,
                        members))
    ordered.sort(key=lambda row: (row[0], row[1]))

    rail = render_rail([("", [(str(len(members)), name, f"#cat-{category}")
                              for _, category, name, _, members in ordered])],
                       note="Every generated plot in this area, grouped by "
                            "the analysis step that produced it.",
                       page_rel=rel)

    body = [f"<header><h1>{esc(area.title)} — figures</h1>",
            f"<p>{len(figs)} generated plot(s) in "
            f"{len(ordered)} group(s). Every one is written by a script in "
            f"this repository and carries the caption its own page gave it. "
            f"Click a plot for the full-size PNG.</p></header>"]
    for _, category, name, page, members in ordered:
        # A caption is a quoted measurement like any other — the CV wall
        # reprints S2's High Gain ceiling inside one — so a group carries the
        # verdict of the page whose figures it is showing.
        body.append(f'<div class="figwall" id="cat-{esc(category)}">'
                    f"<h2>{esc(name)} {_page_verdict(build, page)}</h2>"
                    if page else
                    f'<div class="figwall" id="cat-{esc(category)}">'
                    f"<h2>{esc(name)}</h2>")
        if page:
            body.append(f'<p class="catnote">{len(members)} figure(s) &middot; '
                        f'<a href="{rel_href(rel, page)}">read the analysis '
                        f"they belong to</a></p>")
        else:
            body.append(f'<p class="catnote">{len(members)} figure(s) &middot; '
                        f"captions come from the manuscript's own "
                        f"<code>p5_figure</code> table, so this page and the "
                        f"paper cannot disagree.</p>")
        body.append('<div class="figgrid">')
        for fig in members:
            src = rel_href(rel, fig.src)
            body.append(f'<figure><a href="{src}">'
                        f'<img loading="lazy" src="{src}" '
                        f'alt="{esc(strip_tags(fig.caption) or fig.src)}">'
                        f"</a>")
            if fig.caption:
                body.append(f"<figcaption>{fig.caption}</figcaption>")
            if fig.page:
                anchor = f"#{fig.anchor}" if fig.anchor else ""
                doc = build.pages.get(fig.page)
                where = within(area.title, doc.title) if doc else fig.page
                body.append(f'<div class="figsrc"><a href="'
                            f'{rel_href(rel, fig.page)}{anchor}">in context: '
                            f"{esc(where)}</a></div>")
            body.append("</figure>")
        body.append("</div></div>")

    crumbs = [("Home", LANDING), (area.title, area.home), ("Figures", "")]
    chrome = build.chrome(area, "Figures", rel, rail, crumbs)
    path = build.docs / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(shell(chrome, f"MACRO — {area.title}: figures",
                          "\n".join(body), _footer(build, area, rel)),
                    encoding="utf-8")
    return path


# --- the draft-paper view -------------------------------------------------

_TEX_SECTION_RE = re.compile(r"^\\section\{(.+?)\}", re.M)
_TEX_TITLE_RE = re.compile(r"\\title\{(.*?)\}\s*\n\s*\n", re.S)
_TEX_COMMENT_RE = re.compile(r"(?<!\\)%.*")
#: Where the body of a manuscript stops.  Everything from here on is
#: apparatus — acknowledgments, the bibliography call, the closing tag — and
#: it all sits after the last ``\section``, so leaving it in would credit the
#: final section with prose nobody wrote.  That is not hypothetical: it is
#: what made the T CrB skeleton's Conclusions measure 372 characters when its
#: five siblings measured under 25.
_TEX_ENDMATTER_RE = re.compile(
    r"\\begin\{acknowledgments\}|\\bibliography\b|\\end\{document\}")
#: Control sequences, environments and their braces.  A section made
#: entirely of ``\begin{figure} \includegraphics \end{figure}`` contains no
#: prose, and the question this measurement answers is whether anybody has
#: written anything yet.
_TEX_COMMAND_RE = re.compile(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?|[{}\\&$_^~]")

#: How much prose a section needs before it counts as written.  Calibrated
#: against the two things that actually exist: the AASTeX skeleton five
#: projects share, whose sections hold nothing but a ``% TODO`` comment and
#: measure 0 characters once comments are stripped, and the CV draft, whose
#: shortest real section runs to several thousand.  There is nothing in
#: between, so the threshold is not a fine judgement — it only has to sit
#: inside a very wide gap, and it is stated here so a reader can see what
#: "written" is being taken to mean.
WRITTEN_CHARS = 200


@dataclass(frozen=True)
class TexSection:
    """One ``\\section`` of a manuscript, and whether it has been written."""

    name: str
    is_appendix: bool
    chars: int

    @property
    def written(self) -> bool:
        return self.chars >= WRITTEN_CHARS


def _tex_clean(text: str) -> str:
    """LaTeX heading text, made readable.  Macros become their own names."""
    text = re.sub(r"\\rlmt\b", "RLMT", text)
    text = re.sub(r"\\[a-zA-Z]+\*?\s*", "", text)
    return text.replace("~", " ").replace("{", "").replace("}", "").strip()


def tex_sections(tex: str) -> list[TexSection]:
    """The manuscript's sections, in order, each with its prose length.

    Three things come out before the count, and each of them was putting a
    skeleton on the wrong side of the line: comments (a skeleton section is
    not empty on disk — it holds ``% TODO``), the end matter after the last
    section, and the control sequences themselves.  What is left is prose.
    """
    body, marker, appendix = tex.partition("\\appendix")
    out: list[TexSection] = []
    for chunk, is_app in ((body, False), (appendix, bool(marker))):
        if not chunk:
            continue
        # Per chunk, not once over the whole file: AASTeX puts the
        # acknowledgments BEFORE ``\appendix`` and the bibliography after
        # it, so a single cut at the first end-matter marker would delete
        # the appendix along with the apparatus.
        end_matter = _TEX_ENDMATTER_RE.search(chunk)
        if end_matter:
            chunk = chunk[:end_matter.start()]
        matches = list(_TEX_SECTION_RE.finditer(chunk))
        for i, match in enumerate(matches):
            start = match.end()
            end = matches[i + 1].start() if i + 1 < len(matches) else len(chunk)
            prose = _TEX_COMMENT_RE.sub("", chunk[start:end])
            prose = _TEX_COMMAND_RE.sub(" ", prose)
            out.append(TexSection(_tex_clean(match.group(1)), is_app,
                                  len(" ".join(prose.split()))))
    return out


def render_paper(build: Build, area: Area) -> Optional[Path]:
    """The draft-paper view.

    It shows the manuscript's own shape and says, per section, whether there
    is prose there yet — because "there is a draft" and "the draft is
    finished" are different sentences, and a page that links a PDF without
    distinguishing them invites a collaborator to open a title page.
    """
    pdf_rel = build.paper_pdf(area)
    if not pdf_rel:
        return None
    rel = f"{area.dir}/paper.html"
    tex_path = build.repo / "manuscripts" / area.key / "main.tex"
    tex = tex_path.read_text(encoding="utf-8", errors="replace") \
        if tex_path.exists() else ""
    title_match = _TEX_TITLE_RE.search(tex)
    title = _tex_clean(title_match.group(1)) if title_match else area.title
    sections = tex_sections(tex)
    written = [s for s in sections if s.written]

    pdf_href = "../" * (rel.count("/") + 1) + pdf_rel
    figs = [f for f in build.paper_figs if f.src.startswith(area.dir + "/")]

    rail = render_rail(
        [("Sections", [(str(i + 1), s.name, f"#sec-{i + 1}")
                       for i, s in enumerate(sections) if not s.is_appendix]),
         ("Appendix", [(str(i + 1), s.name, f"#sec-{i + 1}")
                       for i, s in enumerate(sections) if s.is_appendix])],
        note="The draft's own section order, read out of "
             "<code>main.tex</code>.",
        page_rel=rel)

    body = [f"<header><h1>{esc(title)}</h1>",
            f"<p>The draft for {esc(area.title)}, in "
            f"<code>manuscripts/{esc(area.key)}/</code> &middot; "
            f"{len(written)} of {len(sections)} sections carry prose &middot; "
            f"{len(figs)} figure(s) &middot; {esc(area.project.venue)}"
            f"</p></header>",
            '<div class="doors">'
            f'<a class="door" href="{pdf_href}"><span class="step">Read</span>'
            f"<h3>The compiled draft (PDF)</h3>"
            f"<p>The current build of the manuscript, exactly as LaTeX "
            f"produced it.</p></a>"
            f'<a class="door" href="'
            f"https://github.com/jwwetzel/MACRO/blob/main/manuscripts/"
            f'{esc(area.key)}/main.tex"><span class="step">Source</span>'
            f"<h3>main.tex</h3><p>Every number in the text is a macro "
            f"expanded from <code>numbers.tex</code>, which is generated "
            f"from the products database — so the paper and this site quote "
            f"the same measurement or neither does.</p></a>"
            f'<a class="door" href="{rel_href(rel, area.dir + "/figures.html")}">'
            f'<span class="step">Look</span><h3>The figures</h3>'
            f"<p>All {len(figs)} manuscript figures with the captions the "
            f"paper gives them, beside the analysis plots.</p></a></div>",
            '<h2 id="shape">The shape of the argument</h2>',
            '<p class="lead">Read out of the manuscript itself, so this list '
            "cannot fall behind the draft. The figure beside each heading is "
            "its prose length with LaTeX comments and control sequences "
            "removed — a measurement, not a verdict: a short appendix that "
            "is mostly figures is short on purpose.</p>"]
    for i, section in enumerate(sections, start=1):
        kind = "Appendix" if section.is_appendix else f"§{i}"
        # The chip states what was measured.  An earlier version printed
        # "written" / "outline only" from the same threshold that decides
        # whether a DRAFT exists at all — a threshold calibrated on a gap
        # between 15 characters and 4,262, and therefore meaningless applied
        # to a 156-character appendix that is three figures and a sentence.
        # A number a reader can judge beats a label that judges for them.
        chip = ('<span class="chip pending">no prose yet</span>'
                if section.chars < 50 else
                f'<span class="chip done">{section.chars:,} characters'
                f"</span>")
        body.append(f'<div class="qblock" id="sec-{i}">'
                    f'<div class="qhead"><span class="qnum">{kind}</span>'
                    f"<h2>{esc(section.name)}</h2>{chip}</div></div>")

    crumbs = [("Home", LANDING), (area.title, area.home), ("Draft Paper", "")]
    chrome = build.chrome(area, "Draft Paper", rel, rail, crumbs)
    path = build.docs / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(shell(chrome, f"MACRO — {area.title}: draft paper",
                          "\n".join(body), _footer(build, area, rel)),
                    encoding="utf-8")
    return path


# --- the evidence index ---------------------------------------------------

def render_evidence(build: Build, area: Area) -> Path:
    rel = f"{area.dir}/evidence.html"
    pages = build.by_area.get(area.key, [])
    rail = render_rail(
        [("", [(str(i + 1),
                within(area.title,
                       build.pages[p].title) if p in build.pages else p,
                f"#pg-{i + 1}") for i, p in enumerate(pages)])],
        note="Each page below is one analysis step, in the order the "
             "pipeline built it.",
        page_rel=rel)
    body = [f"<header><h1>{esc(area.title)} — evidence detail</h1>",
            f"<p>The full working behind every answer on the Case page: "
            f"{len(pages)} report(s), in build order. These are the long "
            f"pages — every table, every intermediate plot, every caveat. "
            f"Nothing here is a summary.</p></header>"]
    for i, page_rel in enumerate(pages, start=1):
        doc = build.pages.get(page_rel)
        if not doc:
            continue
        verdict = _page_verdict(build, page_rel)
        body.append(f'<div class="qblock" id="pg-{i}">'
                    f'<div class="qhead"><span class="qnum">{i}</span>'
                    f'<h2><a href="{rel_href(rel, page_rel)}">'
                    f"{esc(within(area.title, doc.title))}</a></h2>"
                    f"{verdict}</div>")
        if doc.lead:
            body.append(f'<p class="qask">{esc(clip(doc.lead, 320))}</p>')
        body.append(f'<div class="qfoot"><span>{len(doc.questions)} '
                    f"question(s) &middot; {len(doc.figures)} figure(s)"
                    f"</span><span><code>docs/{esc(page_rel)}</code></span>"
                    f"</div></div>")
    if not pages:
        body.append('<div class="absentnote">This area has no evidence page '
                    "of its own yet. What it rests on is listed on its Case "
                    "page, under <b>What this rests on</b>.</div>")
    crumbs = [("Home", LANDING), (area.title, area.home),
              ("Evidence Detail", "")]
    chrome = build.chrome(area, "Evidence Detail", rel, rail, crumbs)
    path = build.docs / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(shell(chrome, f"MACRO — {area.title}: evidence",
                          "\n".join(body), _footer(build, area, rel)),
                    encoding="utf-8")
    return path


#: Verdict word -> chip class.  Same map the project pages use, and it must
#: stay the same map: DESTROYED gets its own class because a reader has to be
#: able to tell "the numbers moved" from "there is no table any more".
VERDICT_CLASS = {
    pv.FRESH: "v-fresh", pv.STALE: "v-stale", pv.STALE_UPSTREAM: "v-wait",
    pv.NEVER_RUN: "v-never", pv.OUTPUT_MISSING: "v-gone",
    "DESTROYED": "v-gone", "UNKNOWN": "v-wait",
}


def _page_verdict(build: Build, page_rel: str) -> str:
    """The live verdict of the analysis a page publishes, as a chip.

    The verdict comes from ``project_plan.evidence_verdict`` rather than the
    raw freshness state, so this says DESTROYED where the project pages say
    DESTROYED.  Two names for the same condition, on two views of the same
    report, is how a reader learns to trust neither.
    """
    stage_key = analysis_stage_of(page_rel)
    if not stage_key or not build.freshness:
        return ""
    verdict, _why = pp.evidence_verdict(stage_key, build.freshness,
                                        build.fingerprints, build.ever_ran)
    css = VERDICT_CLASS.get(verdict, "v-wait")
    return (f'<span class="chip {css}">{esc(stage_key)}: '
            f"{esc(verdict)}</span>")


# ===========================================================================
# 9.  THE LANDING PAGE
# ===========================================================================

def render_landing(build: Build) -> Path:
    """The front door.

    James's fear is that the site is intimidating, so the job of this page is
    to orient a colleague in thirty seconds and then get out of the way. It
    answers, above the fold and in this order: what this is, how big it is,
    and what to read first — three doors, not sixteen stages. The sixteen
    stages are still there, one click down, for the reader who wants them.
    """
    rel = LANDING
    pipeline_area = build.area_by_key[PIPELINE]
    n_questions = sum(len(d.questions) for d in build.pages.values())
    n_figures = len({f.src for d in build.pages.values() for f in d.figures}
                    | {f.src for f in build.paper_figs})
    n_tasks = len(pp.all_tasks())
    n_pages = len(build.pages)
    drafts = [a for a in build.areas if build.paper_pdf(a)]

    rail = render_rail([("", [
        ("1", "What this is", "#what"),
        ("2", "Where to start", "#start"),
        ("3", "The shared pipeline", "#pipeline"),
        ("4", "The projects", "#projects"),
        ("5", "How to read a page", "#howto")])],
        note="", page_rel=rel)

    cards = []
    for area in build.areas:
        if area.is_pipeline:
            continue
        cond = build.conditions[area.key]
        tabs = [t for t in build.tabs_for(area) if t.rel]
        links = " ".join(
            f'<a href="{rel_href(rel, t.rel)}">{esc(t.label)}</a>'
            for t in tabs)
        n_own = len(build.by_area.get(area.key, []))
        cards.append(
            f'<div class="card"><h3><a href="{rel_href(rel, area.home)}">'
            f"{esc(area.title)}</a></h3>"
            f'<div class="cardmeta">'
            f'<span class="chip {esc(cond.tone)}">{esc(cond.text)}</span>'
            f"<span>{n_own} evidence page(s)</span></div>"
            f'<p class="claim">{esc(clip(area.project.claim, 260))}</p>'
            f'<div class="cardlinks">{links}</div></div>')

    pipe_cond = build.conditions[PIPELINE]
    body = [
        '<div class="hero" id="what">',
        '<p class="kicker">MACRO Consortium &middot; Robert L. Mutel '
        "Telescope</p>",
        "<h1>Six papers, one telescope, and the evidence for every "
        "decision.</h1>",
        "<p>This site is the working record of the MACRO archive: 3.34 TiB "
        "of imaging taken with the 0.5 m Robert L. Mutel Telescope at Winer "
        "Observatory between 2023 and 2026, and the six analyses being built "
        "on it. Every page is generated from the databases it describes, so "
        "no number on this site was typed by a person, and a page that has "
        "gone out of date says so rather than going quiet.</p>",
        f'<div class="factstrip">'
        f'<div><span class="n">6</span><span class="k">projects</span></div>'
        f'<div><span class="n">{n_tasks}</span>'
        f'<span class="k">planned tasks</span></div>'
        f'<div><span class="n">{n_questions}</span>'
        f'<span class="k">questions answered</span></div>'
        f'<div><span class="n">{n_figures}</span>'
        f'<span class="k">generated figures</span></div>'
        f'<div><span class="n">{n_pages}</span>'
        f'<span class="k">evidence reports</span></div>'
        f'<div><span class="n">{len(drafts)}</span>'
        f'<span class="k">draft paper{"" if len(drafts) == 1 else "s"}'
        f"</span></div></div>",
        "</div>",

        '<h2 id="start">Where to start</h2>',
        '<p class="lead">Three doors, depending on why you came. Nobody '
        "needs to read this site in order.</p>",
        '<div class="doors">',
        f'<a class="door" href="{rel_href(rel, pipeline_area.home)}">'
        '<span class="step">Start here</span>'
        "<h3>Can I trust the data?</h3>"
        "<p>The instrument and the archive, argued one question at a time: "
        "what the timestamps mean, what the detector does, which frames have "
        "astrometry, and which FILTER names are actually spectra. Every "
        "project stands on these answers.</p></a>",
    ]
    if drafts:
        area = drafts[0]
        body.append(
            f'<a class="door" href="{rel_href(rel, area.dir + "/paper.html")}">'
            '<span class="step">The science</span>'
            f"<h3>Read the draft paper</h3>"
            f"<p>{esc(SHORT_LABEL.get(area.key, area.title))} is the one "
            f"analysis that has reached a compiled manuscript. It is the "
            f"shortest route to what this archive can actually claim.</p></a>")
    body += [
        f'<a class="door" href="{rel_href(rel, "evidence.html")}">'
        '<span class="step">Everything</span>'
        "<h3>The full evidence index</h3>"
        "<p>Every report, every pipeline stage, each with its current "
        "freshness verdict. The long way round, for the reader who wants "
        "to audit rather than read.</p></a>",
        "</div>",

        '<h2 id="pipeline">The shared pipeline</h2>',
        f'<p class="lead">One characterization of the instrument and the '
        f"archive, done once, that all six projects inherit. It is currently "
        f'<span class="chip {esc(pipe_cond.tone)}">{esc(pipe_cond.text)}'
        f"</span> — and where a stage has gone stale, every project resting "
        f"on it shows the same verdict rather than quietly keeping its old "
        f"answer.</p>",
        f'<div class="doors">'
        f'<a class="door" href="{rel_href(rel, pipeline_area.home)}">'
        f'<span class="step">Case</span><h3>The argument</h3>'
        f"<p>Question by question, in build order.</p></a>"
        f'<a class="door" href="{rel_href(rel, PIPELINE + "/figures.html")}">'
        f'<span class="step">Figures</span><h3>The plots</h3>'
        f"<p>Every diagnostic figure the pipeline produced.</p></a>"
        f'<a class="door" href="'
        f'{rel_href(rel, PIPELINE + "/pipeline_status.html")}">'
        f'<span class="step">Status</span><h3>Freshness &amp; the DAG</h3>'
        f"<p>What is still true, what has been invalidated, and by what."
        f"</p></a></div>",

        '<h2 id="projects">The projects</h2>',
        '<p class="lead">Each opens on the same shape: the case, the '
        "figures, the draft where there is one, the plan, and the full "
        "evidence underneath.</p>",
        f'<div class="cards">{"".join(cards)}</div>',

        '<h2 id="howto">How to read a page here</h2>',
        '<p class="lead">Every page on this site wears the same three '
        "layers, so you always know where you are.</p>",
        "<ol>",
        "<li><b>The top row</b> names the shared pipeline and the six "
        "projects, and each entry carries its own condition — progress for a "
        "project, freshness for the pipeline.</li>",
        "<li><b>The second row</b> is the five views inside one project: the "
        "case, the figures, the draft paper, the plan, and the evidence "
        "detail. A view that does not exist yet says so instead of linking "
        "an empty page.</li>",
        "<li><b>The left rail</b> lists the questions of whatever you are "
        "looking at, in order. Scan it; click the one you doubt.</li>",
        "</ol>",
        '<p class="lead">Nothing on this site is hand-maintained. The '
        "navigation is derived from the plan ledger, the provenance DAG and "
        "the pages' own section headings, which is why it cannot drift away "
        "from what is actually here.</p>",
    ]

    # The landing belongs to no area, so it wears no layer-2 row and nothing
    # in layer 1 is marked current — the reader is above all of it.
    chrome = Chrome(
        area=Area(key="", label="", title="MACRO", is_pipeline=False),
        view="", page_rel=rel, crumbs=[("Home", "")], rail=rail,
        conditions=build.conditions, tabs=[])
    path = build.docs / rel
    path.write_text(
        shell(chrome, "MACRO — the Robert L. Mutel Telescope archive",
              "\n".join(body), _footer(build, None, rel)),
        encoding="utf-8")
    return path


# ===========================================================================
# 10.  WRAPPING THE PAGES THAT ALREADY EXIST
# ===========================================================================

_CSS_LINK_RE = re.compile(r"<link[^>]+macro\.css", re.I)
_HEAD_RE = re.compile(r"(<head[^>]*>)", re.I)


def wrap_page(build: Build, area: Area, page_rel: str,
              siblings: Sequence[str]) -> Path:
    """Put the site chrome around a page a report renderer wrote.

    The renderer's markup is untouched — it is re-emitted verbatim between
    the content markers. What changes is everything around it: the three
    layers, a breadcrumb, and previous/next through this area's evidence
    chain, so the page finally says where it sits.
    """
    path = build.docs / page_rel
    raw = path.read_text(encoding="utf-8")
    head = raw.split("<body", 1)[0]
    content = assign_anchors(original_content(page_body(raw)))
    doc = build.pages.get(page_rel) or harvest_page(path, page_rel)

    rail = render_rail(
        [("", [(str(q.number), q.label,
                f"#{q.anchor}" if q.anchor else "#") for q in doc.questions])],
        back=("Back to the case", area.home),
        note="The questions this page answers, in order.",
        page_rel=page_rel)

    index = siblings.index(page_rel) if page_rel in siblings else -1
    prev = nxt = None
    if index > 0:
        p = siblings[index - 1]
        prev = (build.pages[p].title if p in build.pages else p, p)
    if 0 <= index < len(siblings) - 1:
        n = siblings[index + 1]
        nxt = (build.pages[n].title if n in build.pages else n, n)

    view = ("Freshness & DAG" if page_rel.endswith("pipeline_status.html")
            else "Plan & Status" if page_rel.endswith("index.html")
            else "Evidence Detail")
    crumbs = [("Home", LANDING), (area.title, area.home),
              (view, f"{area.dir}/evidence.html"
               if view == "Evidence Detail" else ""),
              (clip(doc.title, 60), "")]
    crumbs = [c for c in crumbs if c[0]]
    chrome = build.chrome(area, view, page_rel, rail, crumbs, prev, nxt)

    if not _CSS_LINK_RE.search(head):
        css = rel_href(page_rel, "assets/macro.css")
        head = _HEAD_RE.sub(
            lambda m: m.group(1) + f'\n<link rel="stylesheet" href="{css}">',
            head, count=1)

    body = (f"{render_topbar(chrome)}\n{render_tabbar(chrome)}\n"
            f"{render_crumbs(chrome)}\n"
            f'<div class="layout">\n{chrome.rail}\n<div class="pane">\n'
            f"{CONTENT_BEGIN}{content}{CONTENT_END}\n"
            f"{render_prevnext(chrome)}\n"
            f"{_footer(build, area, page_rel)}\n</div>\n</div>\n")
    path.write_text(f"{head}<body>\n{body}</body></html>\n", encoding="utf-8")
    return path


def _footer(build: Build, area: Optional[Area], page_rel: str) -> str:
    """The chrome's own footer.

    It carries NO build timestamp, deliberately.  Every report page already
    stamps the moment its own numbers were queried, which is the stamp that
    means something; a second one saying when the navigation was last
    reassembled would only be a promise that the furniture is recent.  It
    would also make every page differ on every build, which costs the two
    things a generated site most needs: ``build_site.py --check`` telling the
    truth, and a git diff that shows what actually moved.
    """
    live = ("live provenance verdicts and the recorded plan status"
            if build.live else
            "the plan ledger's declared statuses — the manifest was not "
            "readable when this was built, so no live verdict is shown")
    return (
        f'<footer><p>The navigation, the questions and the figure groupings '
        f"on this page are generated by <code>macro_core.site</code> "
        f"({esc(SITE_CODE_VERSION)}) from the plan ledger "
        f"(<code>macro_core.project_plan</code>, {esc(pp.PLAN_CODE_VERSION)}), "
        f"the provenance DAG (<code>macro_core.provenance</code>, "
        f"{esc(pv.PROVENANCE_CODE_VERSION)}), the products databases, and the "
        f"evidence pages' own section headings — never hand-maintained. "
        f"Built from {live}.</p>"
        f'<p>Rebuild with <code>python pipeline/scripts/build_site.py</code>. '
        f'Source: <a href="https://github.com/jwwetzel/MACRO">'
        f"github.com/jwwetzel/MACRO</a>.</p></footer>")


# ===========================================================================
# 11.  ENTRY POINT
# ===========================================================================

DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"


def build_site(manifest: Optional[Path] = None,
               docs_dir: Optional[Path] = None,
               repo_root: Optional[Path] = None) -> list[Path]:
    """Assemble the whole site.  Returns every path written, in order.

    Safe to run at any time and any number of times: the generated views are
    overwritten from scratch, and the wrap of an existing page recovers that
    page's original body first, so running this twice produces the same
    bytes as running it once.
    """
    docs = Path(docs_dir) if docs_dir else DOCS_DIR
    repo = Path(repo_root) if repo_root else REPO_ROOT
    build = Build(docs, repo, Path(manifest) if manifest else DEFAULT_MANIFEST)

    written: list[Path] = [render_landing(build)]
    for area in build.areas:
        written.append(render_case(build, area))
        figures = render_figures(build, area)
        if figures:
            written.append(figures)
        paper = render_paper(build, area)
        if paper:
            written.append(paper)
        written.append(render_evidence(build, area))

    # Wrap everything a renderer wrote, including the project Plan pages and
    # the full evidence index, so no page on the site is without chrome.
    for area in build.areas:
        siblings = build.by_area.get(area.key, [])
        for page_rel in siblings:
            written.append(wrap_page(build, area, page_rel, siblings))
        if not area.is_pipeline:
            plan = f"{area.key}/index.html"
            if (docs / plan).exists():
                written.append(wrap_page(build, area, plan, siblings))
    index = docs / EVIDENCE_INDEX
    if index.exists():
        written.append(_wrap_evidence_index(build, EVIDENCE_INDEX))

    # Only when this really is the repo's own docs/ tree: a build into a
    # temporary copy (the tests, and ``--check``) must not rewrite the repo
    # root's redirect stub out from under it.
    if docs.resolve() == (repo / "docs").resolve():
        written.append(render_root_redirect(build))
    _remove_stale_views(build, written)
    return written


def render_root_redirect(build: Build) -> Path:
    """The repo-root stub GitHub Pages actually serves.

    Pages serves this repository from its ROOT, so ``/MACRO/`` lands here
    rather than on ``docs/index.html``.  The stub redirects, and carries a
    working link list for a reader who arrives with meta-refresh disabled.

    It is generated for the same reason everything else here is: the version
    it replaces listed the six projects by hand, pointing every one of them
    at its Plan page, and it had already fallen behind — the front door of
    the site, hand-typed, naming a set of links that no longer matched the
    site's own shape.
    """
    rows = "".join(
        f'<tr><td>{esc(area.title)}</td>'
        f'<td><a href="docs/{esc(area.home)}">the case</a></td></tr>'
        for area in build.areas)
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MACRO — the Robert L. Mutel Telescope archive</title>
<!-- GENERATED by macro_core.site ({esc(SITE_CODE_VERSION)}); edit that,
     not this.  GitHub Pages serves this repository from its ROOT, so this
     file is what a visitor to jwwetzel.github.io/MACRO/ gets.  The site
     itself lives in docs/; redirect there, and keep working links below for
     a reader with meta-refresh off.  If Pages is later pointed at /docs,
     this file simply stops being served. -->
<meta http-equiv="refresh" content="0; url=docs/index.html">
<link rel="canonical" href="docs/index.html">
<link rel="stylesheet" href="docs/assets/macro.css">
</head><body>
<div class="layout"><div class="pane">
<header>
  <h1>MACRO Consortium — the analysis record</h1>
  <p>Redirecting to the front page&hellip; if nothing happens,
     <a href="docs/index.html">open it here</a>.</p>
</header>
<section><div class="stage">
<table class="data">
<tr><th>Area</th><th>Open at</th></tr>
{rows}
</table>
<p class="sub">Everything else &mdash; the figures, the plans, the draft
paper, the full evidence &mdash; is one row of tabs inside any of these.
Code: <a href="https://github.com/jwwetzel/MACRO">github.com/jwwetzel/MACRO</a></p>
</div></section>
</div></div>
</body></html>
"""
    path = build.repo / "index.html"
    path.write_text(html, encoding="utf-8")
    return path


def _remove_stale_views(build: Build, written: Sequence[Path]) -> None:
    """Delete a generated view that this build no longer produces.

    A view exists because a condition holds — a Figures wall because there
    are figures, a Draft Paper page because prose has been written.  When
    the condition stops holding, the page has to go, or the site keeps a
    door open onto a room that no longer exists.  This is not hypothetical:
    the first build here published a T CrB draft-paper page because a PDF
    existed, and tightening the test to "has anybody written anything yet"
    would have left that page on disk, still linked from nothing, still
    served — the exact orphan class this rebuild set out to close.

    Only ever removes names this module writes.  Nothing a report renderer
    produced can be deleted from here.
    """
    keep = {p.resolve() for p in written}
    for area in build.areas:
        directory = build.docs / area.dir
        if not directory.is_dir():
            continue
        for name in GENERATED_NAMES:
            path = directory / name
            if path.exists() and path.resolve() not in keep:
                path.unlink()


def _wrap_evidence_index(build: Build, page_rel: str) -> Path:
    """The site-wide evidence index belongs to no single area.

    It wears layer 1 and a breadcrumb but no layer-2 row — there is no
    project whose views it is one of, and inventing one would be the kind of
    almost-right navigation this rebuild exists to remove.
    """
    path = build.docs / page_rel
    raw = path.read_text(encoding="utf-8")
    head = raw.split("<body", 1)[0]
    content = assign_anchors(original_content(page_body(raw)))
    doc = harvest_page(path, page_rel)
    rail = render_rail(
        [("", [(str(q.number), q.label, f"#{q.anchor}" if q.anchor else "#")
               for q in doc.questions])],
        back=("Back to the front page", LANDING),
        note="", page_rel=page_rel)
    chrome = Chrome(
        area=Area(key="", label="", title="All evidence", is_pipeline=False),
        view="", page_rel=page_rel,
        crumbs=[("Home", LANDING), ("All evidence", "")],
        rail=rail, conditions=build.conditions, tabs=[])
    body = (f"{render_topbar(chrome)}\n{render_crumbs(chrome)}\n"
            f'<div class="layout">\n{chrome.rail}\n<div class="pane">\n'
            f"{CONTENT_BEGIN}{content}{CONTENT_END}\n"
            f"{_footer(build, None, page_rel)}\n</div>\n</div>\n")
    path.write_text(f"{head}<body>\n{body}</body></html>\n", encoding="utf-8")
    return path
