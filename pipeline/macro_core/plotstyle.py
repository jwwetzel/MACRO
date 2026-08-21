"""macro_core.plotstyle — the ONE definition of the house figure style.

WHY THIS MODULE EXISTS
----------------------
Until this module existed there were two house styles in the repository and
neither knew about the other.  ``macro_core.report_s0`` carried a ``DARK``
rcParams dict — figures drawn light-on-dark so they would sit invisibly on
the site's dark page background — and every other report renderer imported
it, adding its own private pastel constants (``GOOD = "#9fd8ae"``,
``BAD = "#f0a3a3"``, a ``FILTER_COLOR`` map per file) that were tuned by eye
against that dark ground.  Meanwhile ``macro_phot.figures_cv`` had been
written to a genuinely publication-grade standard for AASTeX: white ground,
Okabe--Ito palette, marker shape as a mandatory second channel, units on
every axis, two fixed column widths.  Eighty-five figures obeyed the first
convention and thirteen obeyed the second, and a reader moving between them
could not tell whether a colour change meant something.

This module is the second convention, extracted so it has exactly one
definition.  ``figures_cv`` imports its palette and rcParams from here
rather than declaring them; every report renderer imports from here too.
There is no dark style any more.

THE RULES, AND WHY EACH ONE IS A RULE
-------------------------------------
**White ground.**  A figure is a printable object.  A dark-background PNG
cannot be put in a journal, cannot be pasted into a talk with a light
template, and prints as a black rectangle.  The web page can be dark or
light; the figure is neither, it is white.

**Okabe--Ito, and colour is never the only channel.**  The eight-hue
Okabe--Ito set stays distinguishable under deuteranopia, protanopia and
tritanopia.  That is necessary and not sufficient: a figure also gets
photocopied, printed in greyscale, and read on a projector that crushes
saturation.  So every categorical series that is given a colour is also
given a distinct MARKER or LINE STYLE.  :func:`series` hands out both at
once so a caller cannot take the colour and forget the marker, and it pairs
them as a Latin square, so eight hues carry SIXTY-FOUR distinguishable
series.  That is why nothing here reaches for ``tab20`` when a panel has
twelve lines on it: ``tab20`` is neither colour-blind safe nor greyscale
safe, and it was the second palette hiding in this repository after the
dark one was removed.  A set of series with a NATURAL ORDER — nights in a
run, rungs of a ladder — takes :func:`ordinal_colors` instead, because
there the reader should see a sequence and not consult a legend.

**Yellow is not in the cycle.**  ``#F0E442`` is part of Okabe--Ito and is
close to invisible as a line or a small marker on white.  It is kept in
:data:`OKABE_ITO` because the set is quoted as a set, and excluded from
:data:`CYCLE`; use it only as a large filled area with an edge.

**A floor is not a measurement.**  A detection limit, an upper limit, a
censored night, a resolution-limited timing residual: these are statements
about what the data COULD have seen, and drawing them with the same marker
as a measurement invites a reader to fit a trend through them.  The house
convention, learned the expensive way in ``figures_cv`` (see the Figure 13
audit and the Figure 8 censored nights), is an OPEN downward triangle with
no fill: :func:`floor_kw`.  Measurements are filled markers with a thin
dark edge: :func:`measurement_kw`.

**Units on axes.**  Not enforceable in code, stated here because this is
where the house rules live.  "amplitude" is not an axis label; "amplitude
(mmag)" is.

**Two widths, two profiles.**  ``COL_SINGLE`` (3.5 in) and ``COL_DOUBLE``
(7.1 in) are the AASTeX column widths; a figure drawn at any other width is
rescaled by the journal and its type stops matching the body text.  The
``print`` profile is what the manuscript figures use.  The ``web`` profile
is the same palette and the same ink at a slightly larger type size and a
lower raster resolution, for the report pages, whose figures are read on a
screen at a metre rather than in a column at forty centimetres.  The two
profiles differ ONLY in type size, dpi and whether the grid is on by
default; every colour is shared.

USAGE
-----
::

    from macro_core import plotstyle as ps

    with ps.context():                       # web profile
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.plot(x, y, **ps.measurement_kw(ps.ACCENT))
        ax.plot(x, lim, **ps.floor_kw(ps.MUTED))
        ax.set_xlabel("exposure time (s)")
        fig.savefig(path, dpi=ps.WEB_DPI)

    ps.apply("print")                        # manuscript figures
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Sequence

import matplotlib
matplotlib.use("Agg")                        # headless: we only write files
import matplotlib.pyplot as plt              # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D          # noqa: E402

# ===========================================================================
# Page geometry
# ===========================================================================

#: AASTeX two-column body: one column is 3.5 in, the full text block 7.1 in.
COL_SINGLE = 3.5
COL_DOUBLE = 7.1

#: Raster resolution for the manuscript's PNG companions.  The PDF that
#: LaTeX includes is vector and unaffected by this.
PNG_DPI = 200

#: Raster resolution for the report pages.  Above the 110 the site spec
#: requires, below PNG_DPI because these are 7-to-9-inch canvases and the
#: page has ninety-odd of them.
WEB_DPI = 130

# ===========================================================================
# Palette
# ===========================================================================

#: Okabe--Ito: eight hues chosen to stay distinguishable under the three
#: common forms of colour-blindness, plus the neutral grey the set is
#: conventionally quoted with.  Named, not indexed, so a reader of the code
#: can see which physical thing each colour means.
OKABE_ITO = {
    "black":     "#000000",
    "orange":    "#E69F00",
    "skyblue":   "#56B4E9",
    "green":     "#009E73",
    "yellow":    "#F0E442",
    "blue":      "#0072B2",
    "vermilion": "#D55E00",
    "purple":    "#CC79A7",
    "grey":      "#999999",
}

#: The ground every figure is drawn on, and the ink it is drawn in.  ``INK``
#: is not pure black: at 8 pt a pure-black serif on pure white shimmers, and
#: #1a1a1a prints identically while reading calmer on screen.
PAPER = "#ffffff"
INK = "#1a1a1a"
#: Spine and tick colour — a shade lighter than the type so the frame does
#: not compete with the data.
RULE = "#4d4d4d"
#: Grid lines.  Faint enough to read through, dark enough to survive a
#: 130-dpi raster.  One value, because three different near-whites across
#: three files is precisely the drift this module exists to stop.
GRID = "#e3e3e3"

#: Semantic roles.  Report code should name the ROLE, not the hue: what a
#: reader must learn is that orange means "look at this", not that orange
#: means "#E69F00".
ACCENT = OKABE_ITO["blue"]        #: the primary data series
WARN = OKABE_ITO["orange"]        #: outliers, disagreements, thresholds
GOOD = OKABE_ITO["green"]         #: confirmations, passes, repaired rows
BAD = OKABE_ITO["vermilion"]      #: contradictions, failures, vetoes
OTHER = OKABE_ITO["purple"]       #: a fourth category when three are taken
SECOND = OKABE_ITO["skyblue"]     #: a second shade of the primary role
#: De-emphasised but still readable: reference lines, annotations, context
#: series a reader should see past rather than at.
MUTED = "#666666"
#: Genuinely backgrounded: rejected points, "neither" categories.
FAINT = OKABE_ITO["grey"]
#: Drawn BEHIND the data and not meant to be read as a series at all: the
#: hairline joining two markers of the same row, the cloud of points a panel
#: is showing the absence of structure in.  Lighter than FAINT on purpose,
#: and one value rather than the five near-greys this replaced.
WISP = "#c9c9c9"

#: The categorical cycle, in the order a figure should consume it.  Yellow
#: is deliberately absent (see the module docstring).
CYCLE = (ACCENT, BAD, GOOD, WARN, OTHER, SECOND, INK, FAINT)

#: Marker shapes, in the order a figure should consume them.  Paired with
#: :data:`CYCLE` by :func:`series` so colour is never load-bearing alone.
MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")

#: Line styles, for the cases where a series is a line and not points.
LINESTYLES = ("-", "--", ":", "-.", (0, (3, 1, 1, 1)), (0, (5, 1)),
              (0, (1, 1)), (0, (4, 2, 1, 2)))

#: The marker that means "this is a limit, not a measurement".  Drawn open.
FLOOR_MARKER = "v"

# ---------------------------------------------------------------------------
# The recurring categorical maps, defined once for the whole repository.
# ---------------------------------------------------------------------------

#: Filter -> colour.  Both the 2024 Johnson-ish G/R/I set and the Sloan
#: g/r/i set map to the SAME three hues on purpose: the era-seam argument
#: this archive keeps making is that these are different bandpasses
#: measuring the same three parts of the spectrum, and a reader comparing
#: two era panels should see blue against blue.  The bandpass difference is
#: stated in the axis label and the caption, never smuggled in as a colour.
BAND_COLOR = {
    "G": OKABE_ITO["blue"],      "g": OKABE_ITO["blue"],
    "R": OKABE_ITO["vermilion"], "r": OKABE_ITO["vermilion"],
    "I": OKABE_ITO["green"],     "i": OKABE_ITO["green"],
    "z": OKABE_ITO["purple"],    "y": OKABE_ITO["orange"],
    "V": OKABE_ITO["skyblue"],   "B": OKABE_ITO["purple"],
    "C": OKABE_ITO["black"],     "": OKABE_ITO["grey"],
}

#: Filter -> marker.  The second channel that makes a band-coloured figure
#: survive a greyscale print.
BAND_MARKER = {
    "G": "o", "g": "o",
    "R": "s", "r": "s",
    "I": "^", "i": "^",
    "z": "D", "y": "v",
    "V": "P", "B": "X", "C": "*", "": "o",
}

#: Readout mode -> marker.  Camera and epoch are confounded across most of
#: this archive, so the camera has to be visible in any figure that spans a
#: seam — otherwise a reader attributes an instrument step to the sky.
MODE_MARKER = {
    "High Gain": "o",
    "High Gain StackPro": "P",
    "1MHz High Sensitivity 16-bit": "s",
    "Mode0": "^",
    "Fast": "D",
}

#: Calibration kind -> colour.  Three roles, three separable hues.
KIND_COLOR = {"bias": GOOD, "dark": ACCENT, "flat": WARN}


# ===========================================================================
# rcParams
# ===========================================================================

#: Everything both profiles share: the ground, the ink, the frame, the
#: fonts, and the two settings (``pdf.fonttype``/``ps.fonttype`` = 42) that
#: embed TrueType rather than Type-3.  Type-3 is the matplotlib default and
#: several publishers reject it outright.
_BASE = {
    "figure.facecolor": PAPER,
    "savefig.facecolor": PAPER,
    "axes.facecolor": PAPER,
    "axes.edgecolor": RULE,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": RULE,
    "ytick.color": RULE,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "grid.color": GRID,
    "grid.alpha": 1.0,
    # The grid is a reading aid, not data.  Without this matplotlib
    # paints it OVER the bars and the markers, and a pale grid line
    # across a filled bar reads as a boundary in the data.
    "axes.axisbelow": True,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
    "legend.frameon": False,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "axes.prop_cycle": plt.cycler(color=list(CYCLE)),
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
}

#: Manuscript figures: 8 pt type against an AASTeX caption, no grid (the
#: journal column is narrow and a grid there is noise), hairline rules.
_PRINT = {
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "axes.linewidth": 0.7,
    "grid.linewidth": 0.4,
    "lines.linewidth": 1.0,
    "lines.markersize": 3.0,
    "axes.grid": False,
    "figure.dpi": 110,
}

#: Report-page figures: read on a screen on a canvas two to three times the
#: area, so the type goes up and the grid comes on.  Same ink, same hues.
_WEB = {
    "font.size": 10,
    "axes.labelsize": 10,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.linewidth": 0.9,
    "grid.linewidth": 0.6,
    "lines.linewidth": 1.3,
    "lines.markersize": 4.0,
    "axes.grid": True,
    "figure.dpi": WEB_DPI,
}

PROFILES = {"print": _PRINT, "web": _WEB}


def rc(profile: str = "web", **overrides) -> dict:
    """The rcParams dict for one profile, with optional local overrides.

    Returned rather than applied so a caller can hand it to
    ``plt.rc_context`` itself, which is what the report renderers do.
    """
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}; "
                         f"expected one of {sorted(PROFILES)}")
    out = dict(_BASE)
    out.update(PROFILES[profile])
    out.update(overrides)
    return out


def apply(profile: str = "web", **overrides) -> None:
    """Set the house style globally for the rest of the process.

    Used by the manuscript figure driver, which draws thirteen figures in
    one run and wants the style set once.  Report renderers use
    :func:`context` instead, so that importing a renderer never mutates a
    caller's matplotlib state.
    """
    plt.rcParams.update(rc(profile, **overrides))


@contextmanager
def context(profile: str = "web", **overrides):
    """Draw inside the house style, and leave matplotlib as it was found."""
    with plt.rc_context(rc(profile, **overrides)):
        yield


#: Backwards-compatible aliases.  ``figures_cv.apply_style`` was the print
#: profile before this module existed and is still called by
#: ``pipeline/scripts/run_cv_paper.py``; ``STYLE`` is the dict the report
#: renderers pass to ``plt.rc_context``.  Both are this module's output.
STYLE = rc("web")
PRINT_STYLE = rc("print")


# ===========================================================================
# Named helpers for the recurring cases
# ===========================================================================
def _pair(i: int, second: Sequence) -> int:
    """Index into a second channel so that (colour, channel) never repeats.

    A Latin square, not two independent cycles.  Within each block of eight
    series the colour index is ``i % 8`` and this is ``(i % 8 + i // 8) % 8``
    — a permutation, so every series in a block has BOTH a unique hue and a
    unique shape.  Across blocks the pair advances, giving 64 combinations
    before anything repeats.  Two independent ``% 8`` cycles would instead
    lock hue to shape and repeat after eight, which is how a twenty-series
    Allan plot ends up needing ``tab20`` and losing colour-blind safety.
    """
    n = len(second)
    return (i % len(CYCLE) + i // len(CYCLE)) % n


def series(i: int) -> dict:
    """Colour AND marker for the ``i``-th categorical series.

    Handing both out together is the point: a caller that asks for a colour
    gets a shape with it, so no figure can end up relying on hue alone.
    """
    return {"color": CYCLE[i % len(CYCLE)],
            "marker": MARKERS[_pair(i, MARKERS)]}


def line_series(i: int) -> dict:
    """Colour AND line style for the ``i``-th series drawn as a line.

    The pairing is what lets a figure carry more series than the palette has
    hues without reaching for a rainbow map: eight colours times eight dash
    patterns, all of them colour-blind safe and all of them surviving a
    greyscale print.
    """
    return {"color": CYCLE[i % len(CYCLE)],
            "linestyle": LINESTYLES[_pair(i, LINESTYLES)]}


def ordinal_colors(n: int):
    """``n`` colours along :data:`SEQ_CMAP` for an ORDERED set of series.

    Nights in a run, rungs of a ladder — categories with a natural order,
    where a reader should see the sequence rather than look each one up in a
    legend.  Sampling starts at 0.3 rather than 0.0 because the pale end of
    the ramp is the "empty cell" colour and a line drawn in it is invisible.
    """
    import numpy as _np
    if int(n) <= 0:
        return []
    return [SEQ_CMAP(x) for x in _np.linspace(0.30, 1.0, int(n))]


def measurement_kw(color: str = ACCENT, marker: str = "o",
                   size: float | None = None, **kw) -> dict:
    """Plot kwargs for a MEASURED point: filled, with a thin dark edge.

    The dark edge is what keeps a cloud of light-hued markers legible where
    they overlap, and what distinguishes a measurement at a glance from the
    open :func:`floor_kw` marker of a limit.
    """
    out = {"marker": marker, "linestyle": "none",
           "markerfacecolor": color, "markeredgecolor": INK,
           "markeredgewidth": 0.35}
    if size is not None:
        out["markersize"] = size
    out.update(kw)
    return out


def floor_kw(color: str = MUTED, marker: str = FLOOR_MARKER,
             size: float | None = None, **kw) -> dict:
    """Plot kwargs for a FLOOR: an open downward triangle, no fill.

    A floor is an upper limit, a detection threshold, a resolution-limited
    residual — a statement about what could have been seen, not about what
    was.  Drawn hollow so that a reader never fits a trend through it.
    """
    out = {"marker": marker, "linestyle": "none",
           "markerfacecolor": "none", "markeredgecolor": color,
           "markeredgewidth": 0.9}
    if size is not None:
        out["markersize"] = size
    out.update(kw)
    return out


def floor_handle(label: str, color: str = MUTED) -> Line2D:
    """A legend entry for the floor marker, drawn the way the data are."""
    return Line2D([], [], label=label, **floor_kw(color, size=4.0))


def measurement_handle(label: str, color: str = ACCENT,
                       marker: str = "o") -> Line2D:
    """A legend entry for a measured series, drawn the way the data are."""
    return Line2D([], [], label=label, **measurement_kw(color, marker,
                                                        size=4.0))


def reference_kw(color: str = MUTED, style: str = "--", **kw) -> dict:
    """Plot kwargs for a reference line: a threshold, a zero, an identity.

    Muted and dashed, because it is the thing the data are compared AGAINST
    and must not read as another series.
    """
    out = {"color": color, "linestyle": style, "linewidth": 1.0}
    out.update(kw)
    return out


def tint(color: str, amount: float = 0.55) -> str:
    """A paler version of a house colour, mixed ``amount`` toward the paper.

    For the "same thing, derived" case ONLY — a master calibration frame
    beside the raw frames it was built from, the box of a distribution whose
    whiskers are drawn in the full hue.  A tint is not a new category: two
    unrelated series must differ by hue and marker, not by saturation, which
    is the axis a greyscale print destroys first.
    """
    a = float(min(max(amount, 0.0), 1.0))
    r, g, b = matplotlib.colors.to_rgb(color)
    pr, pg, pb = matplotlib.colors.to_rgb(PAPER)
    return matplotlib.colors.to_hex((r + (pr - r) * a,
                                     g + (pg - g) * a,
                                     b + (pb - b) * a))


def band_color(filt) -> str:
    """Colour for a filter name, falling back to grey for an unknown one."""
    return BAND_COLOR.get(str(filt or ""), FAINT)


def band_marker(filt) -> str:
    """Marker for a filter name, falling back to a circle."""
    return BAND_MARKER.get(str(filt or ""), "o")


# ---------------------------------------------------------------------------
# Colormaps.  Two, both readable on white, both ending dark enough that a
# full-scale cell is unambiguous at 130 dpi.
# ---------------------------------------------------------------------------
#: The ONE sequential ramp: near-white -> Okabe--Ito blue -> deep blue,
#: for any "how much of this is there" heatmap.  It starts a shade off the
#: paper rather than at pure white so that an EMPTY cell is still visibly a
#: cell: a zero drawn in #ffffff on a #ffffff canvas erases the grid the
#: reader is counting rows in.  It ends dark enough that a full-scale cell
#: is unambiguous at 130 dpi, and it is monotonic in lightness, so it also
#: survives the greyscale print that a rainbow map does not.
SEQ_CMAP = LinearSegmentedColormap.from_list(
    "macro_seq", ["#f4f8fb", "#9ecae1", ACCENT, "#02405f"])

#: Vermilion -> pale -> blue.  For a signed quantity where zero is
#: meaningful.  Vermilion and blue rather than red and green, because red
#: against green is the one pair a deuteranope cannot separate.
DIV_CMAP = LinearSegmentedColormap.from_list(
    "macro_div", [BAD, "#f2c9ae", "#f4f4f4", "#9ecae1", ACCENT])

#: PICTURES OF THINGS are not heatmaps of quantities, and they keep a
#: perceptual astronomical ramp rather than the house blue: a 512-pixel
#: stamp of a dark frame is read for its texture, not for a value looked up
#: against a key.  Named here so that the one place it is allowed is
#: written down, and so a change reaches every stamp at once.
IMAGE_CMAP = "magma"
IMAGE_FLAT_CMAP = "viridis"


def text_on(color: str) -> str:
    """Legible type colour for a label written ON a solid ``color``.

    Chosen by relative luminance (WCAG), not by eye: Okabe--Ito orange is
    light enough that white type on it fails contrast while black passes,
    and Okabe--Ito blue is the other way round.  Deciding that once, here,
    is what stops a badge somewhere on the site being unreadable because
    its palette entry changed and its text colour did not.
    """
    r, g, b = matplotlib.colors.to_rgb(color)

    def lin(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4

    lum = 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)
    # Contrast against white vs against near-black, the WCAG ratio both ways.
    return PAPER if (1.05 / (lum + 0.05)) >= ((lum + 0.05) / 0.06) else INK


def ink_on(fraction: float, threshold: float = 0.55) -> str:
    """Text colour for a label written ON a :data:`SEQ_CMAP` cell.

    ``fraction`` is the cell's position in 0..1 along the colormap.  Past
    ``threshold`` the cell is dark enough that only white type is legible;
    below it, only dark type is.  Written once here so that every heatmap in
    the repository makes the same decision.
    """
    return PAPER if float(fraction) >= threshold else INK
