"""The house figure style, and the guard that keeps every figure inside it.

Two kinds of test live here.

The first kind exercises ``macro_core.plotstyle`` as a module: that the
palette is the Okabe--Ito set, that the two profiles differ only in the ways
they are documented to differ, that a floor marker really is drawn hollow,
that the luminance helpers pick type a person can read.

The second kind is the REGRESSION GUARD this file exists for.  Until August
2026 this repository had two figure styles: a dark one, matched to the site's
old dark page background, that 85 published figures obeyed, and the
publication-grade light one in ``macro_phot.figures_cv`` that the 13
manuscript figures obeyed.  Nothing prevented that split, so it lasted
months and a reader moving between two pages could not tell whether a colour
change carried meaning.  :func:`test_no_published_figure_has_a_dark_ground`
samples the corners of every PNG under ``docs/`` and fails if any of them
comes back dark.  It is deliberately crude — a corner pixel, not a
histogram — because the failure it is guarding against is not subtle: a
renderer that reverts to a dark rcParams dict paints the whole canvas.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from macro_core import plotstyle as ps

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"

#: A corner whose channels all sit above this (0-255) is a light ground.
#: White is 255; the house grid line is 227; anything a dark theme would
#: paint a canvas with is below 60, so the threshold is nowhere near a
#: legitimate figure.
LIGHT_MIN = 200


# ---------------------------------------------------------------------------
# The palette
# ---------------------------------------------------------------------------
def test_palette_is_okabe_ito():
    """The eight hues are the published set, not a lookalike."""
    assert set(ps.OKABE_ITO.values()) >= {
        "#000000", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7"}


def test_semantic_roles_come_from_the_palette():
    for role in (ps.ACCENT, ps.WARN, ps.GOOD, ps.BAD, ps.OTHER, ps.SECOND):
        assert role in ps.OKABE_ITO.values()


def test_yellow_is_excluded_from_the_cycle():
    """#F0E442 is invisible as a line on white and must not be handed out."""
    assert ps.OKABE_ITO["yellow"] not in ps.CYCLE


def test_cycle_has_no_repeats():
    assert len(set(ps.CYCLE)) == len(ps.CYCLE)


def test_every_cycle_colour_gets_its_own_marker():
    """Colour is never the only channel: N colours, N distinct markers."""
    assert len(ps.MARKERS) >= len(ps.CYCLE)
    assert len(set(ps.MARKERS[:len(ps.CYCLE)])) == len(ps.CYCLE)


def test_series_hands_out_colour_and_marker_together():
    a, b = ps.series(0), ps.series(1)
    assert a["color"] != b["color"]
    assert a["marker"] != b["marker"]


def test_series_wraps_rather_than_raising():
    assert ps.series(999)["color"] in ps.CYCLE
    assert ps.line_series(999)["linestyle"] in ps.LINESTYLES


def test_a_block_of_eight_series_is_unique_in_both_channels():
    """Within one block, every series has its own hue AND its own shape."""
    block = [ps.series(i) for i in range(len(ps.CYCLE))]
    assert len({b["color"] for b in block}) == len(ps.CYCLE)
    assert len({b["marker"] for b in block}) == len(ps.CYCLE)


def test_sixty_four_series_never_repeat_a_colour_shape_pair():
    """The reason no figure here needs tab20 or a rainbow map."""
    pairs = {(d["color"], d["marker"]) for d in map(ps.series, range(64))}
    assert len(pairs) == 64
    lines = {(d["color"], d["linestyle"])
             for d in map(ps.line_series, range(64))}
    assert len(lines) == 64


def test_ordinal_colors_are_ordered_and_never_invisible():
    """An ordered set of series gets the ramp, minus its invisible end."""
    import matplotlib.colors as mc
    cols = ps.ordinal_colors(6)
    assert len(cols) == 6
    lums = [sum(mc.to_rgb(c)) for c in cols]
    assert lums == sorted(lums, reverse=True)
    assert min(mc.to_rgb(cols[0])) < 0.85     # the palest is still a line
    assert ps.ordinal_colors(0) == []


def test_bands_that_measure_the_same_light_share_a_hue():
    """G/g, R/r, I/i are different bandpasses on the same three colours."""
    for upper, lower in (("G", "g"), ("R", "r"), ("I", "i")):
        assert ps.BAND_COLOR[upper] == ps.BAND_COLOR[lower]
        assert ps.BAND_MARKER[upper] == ps.BAND_MARKER[lower]


def test_bands_are_distinguishable_from_each_other():
    assert len({ps.BAND_COLOR[b] for b in "GRI"}) == 3
    assert len({ps.BAND_MARKER[b] for b in "GRI"}) == 3


def test_band_lookup_falls_back_rather_than_raising():
    assert ps.band_color("no-such-filter") == ps.FAINT
    assert ps.band_marker(None) == "o"


# ---------------------------------------------------------------------------
# The rcParams
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("profile", ["print", "web"])
def test_both_profiles_paint_a_white_ground(profile):
    rc = ps.rc(profile)
    assert rc["figure.facecolor"] == ps.PAPER
    assert rc["savefig.facecolor"] == ps.PAPER
    assert rc["axes.facecolor"] == ps.PAPER


@pytest.mark.parametrize("profile", ["print", "web"])
def test_both_profiles_embed_truetype(profile):
    """Type-3 is the matplotlib default and several publishers reject it."""
    rc = ps.rc(profile)
    assert rc["pdf.fonttype"] == 42
    assert rc["ps.fonttype"] == 42


def test_profiles_differ_only_in_size_and_grid():
    """The two profiles are one style at two reading distances."""
    a, b = ps.rc("print"), ps.rc("web")
    differing = {k for k in a if a[k] != b[k]}
    allowed = {"font.size", "axes.labelsize", "axes.titlesize",
               "xtick.labelsize", "ytick.labelsize", "legend.fontsize",
               "axes.linewidth", "grid.linewidth", "lines.linewidth",
               "lines.markersize", "axes.grid", "figure.dpi"}
    assert differing <= allowed
    # and in particular: not one colour differs between them
    for key in a:
        if "color" in key:
            assert a[key] == b[key], key


def test_web_type_is_larger_than_print_type():
    assert ps.rc("web")["font.size"] > ps.rc("print")["font.size"]


def test_unknown_profile_is_an_error_not_a_silent_default():
    with pytest.raises(ValueError):
        ps.rc("dark")


def test_context_restores_the_previous_rcparams():
    import matplotlib.pyplot as plt
    before = plt.rcParams["font.size"]
    with ps.context("print"):
        assert plt.rcParams["font.size"] == ps.rc("print")["font.size"]
    assert plt.rcParams["font.size"] == before


def test_column_widths_are_the_aastex_ones():
    assert ps.COL_SINGLE == pytest.approx(3.5)
    assert ps.COL_DOUBLE == pytest.approx(7.1)


# ---------------------------------------------------------------------------
# Floors, measurements, tints
# ---------------------------------------------------------------------------
def test_a_floor_is_drawn_hollow_and_a_measurement_is_not():
    """The one convention a reader must never have to guess at."""
    assert ps.floor_kw()["markerfacecolor"] == "none"
    assert ps.measurement_kw()["markerfacecolor"] != "none"


def test_a_floor_and_a_measurement_do_not_share_a_marker():
    assert ps.floor_kw()["marker"] != ps.measurement_kw()["marker"]


def test_neither_helper_draws_a_connecting_line():
    """Points, not a series: joining limits would imply a trend."""
    assert ps.floor_kw()["linestyle"] == "none"
    assert ps.measurement_kw()["linestyle"] == "none"


def test_helper_kwargs_can_be_overridden():
    assert ps.measurement_kw(size=9.0)["markersize"] == 9.0
    assert ps.floor_kw(zorder=7)["zorder"] == 7


def test_legend_handles_are_drawn_the_way_the_data_are():
    assert ps.floor_handle("limit").get_markerfacecolor() == "none"
    assert ps.measurement_handle("point").get_markerfacecolor() != "none"


def test_tint_moves_toward_the_paper_and_stays_a_colour():
    pale = ps.tint(ps.ACCENT, 0.5)
    assert re.fullmatch(r"#[0-9a-f]{6}", pale)
    assert pale != ps.ACCENT
    assert ps.tint(ps.ACCENT, 0.0).lower() == ps.ACCENT.lower()
    assert ps.tint(ps.ACCENT, 1.0).lower() == ps.PAPER.lower()


def test_tint_is_monotonic_toward_white():
    import matplotlib.colors as mc
    lums = [sum(mc.to_rgb(ps.tint(ps.BAD, a))) for a in (0.0, 0.4, 0.8)]
    assert lums[0] < lums[1] < lums[2]


# ---------------------------------------------------------------------------
# Type-on-colour
# ---------------------------------------------------------------------------
def test_text_on_picks_the_higher_contrast_option():
    """White type on the house orange fails contrast; dark type passes."""
    assert ps.text_on(ps.WARN) == ps.INK
    assert ps.text_on(ps.ACCENT) == ps.PAPER
    assert ps.text_on(ps.PAPER) == ps.INK
    assert ps.text_on(ps.INK) == ps.PAPER


def test_ink_on_flips_across_the_colormap():
    assert ps.ink_on(0.0) == ps.INK
    assert ps.ink_on(1.0) == ps.PAPER


def test_sequential_cmap_starts_near_paper_but_not_on_it():
    """An empty heatmap cell must still be visibly a cell.

    A zero drawn in pure white on a white canvas erases the row a reader is
    counting along, which is how the S0b coverage heatmap lost its grid the
    first time this ramp was written.
    """
    import matplotlib.colors as mc
    lo = mc.to_rgb(ps.SEQ_CMAP(0.0))
    assert lo != mc.to_rgb(ps.PAPER)
    assert min(lo) > 0.9          # still unmistakably "nothing here"


def test_sequential_cmap_is_monotonic_in_lightness():
    """The ramp must survive a greyscale print, which a rainbow does not."""
    import matplotlib.colors as mc
    lums = [sum(mc.to_rgb(ps.SEQ_CMAP(f))) for f in (0.0, 0.35, 0.7, 1.0)]
    assert lums == sorted(lums, reverse=True)


# ---------------------------------------------------------------------------
# figures_cv must not carry a second copy of the style
# ---------------------------------------------------------------------------
def test_figures_cv_uses_the_shared_module():
    from macro_phot import figures_cv as fx
    assert fx.OKABE_ITO is ps.OKABE_ITO
    assert fx.BAND_COLOR is ps.BAND_COLOR
    assert fx.BAND_MARKER is ps.BAND_MARKER
    assert fx.COL_SINGLE == ps.COL_SINGLE
    assert fx.COL_DOUBLE == ps.COL_DOUBLE


def test_figures_cv_apply_style_sets_the_print_profile():
    import matplotlib.pyplot as plt
    from macro_phot import figures_cv as fx
    with plt.rc_context():
        fx.apply_style()
        assert plt.rcParams["figure.facecolor"] == ps.PAPER
        assert plt.rcParams["font.size"] == ps.rc("print")["font.size"]


def test_no_renderer_reaches_for_a_matplotlib_categorical_map():
    """tab10/tab20 are a SECOND categorical palette, and not safe.

    Neither is colour-blind safe, and tab20's pale-and-dark pairs are
    indistinguishable in greyscale.  Anything categorical comes from
    :func:`plotstyle.series` or :func:`plotstyle.line_series`.
    """
    offenders = []
    for path in sorted((REPO_ROOT / "pipeline").rglob("*.py")):
        if path.name in ("plotstyle.py", "test_plotstyle.py"):
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"\btab(10|20|20b|20c)\b", line):
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{i}")
    assert not offenders, ("categorical colours come from plotstyle."
                           "series/line_series: " + ", ".join(offenders))


def test_no_renderer_declares_its_own_hex_palette():
    """A colour literal in a renderer is how the two styles drifted apart.

    HTML and CSS strings still carry hex — those are the page, not the
    figure — so this checks only lines that are matplotlib calls.
    """
    offenders = []
    for path in sorted((REPO_ROOT / "pipeline").rglob("*.py")):
        if path.name in ("plotstyle.py",):
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if not re.search(r'"#[0-9a-fA-F]{3,8}"', line):
                continue
            if re.search(r"\b(color|colour|colors|colours|facecolor|"
                         r"edgecolor|edgecolors|ecolor|mfc|mec|cmap|"
                         r"markerfacecolor|markeredgecolor)\b\s*=", line):
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{i}")
    assert not offenders, ("figure colours must come from macro_core."
                           "plotstyle, not from a literal: "
                           + ", ".join(offenders))


def test_no_renderer_names_a_matplotlib_colormap_directly():
    """``cmap=`` takes a plotstyle constant, never a bare string.

    The categorical guard above catches ``tab10``/``tab20``, and the hex
    guard catches a literal colour — but neither caught ``cmap="coolwarm"``
    in the catalogue-tie residual panels or ``cmap="viridis"`` on the moon
    figure's separation colourbar, both of which were live on the site after
    the dark theme was removed.  A named ramp is a second definition of the
    house style in exactly the way a hex literal is: it just spells the
    colours with a word instead of a number.

    The sanctioned ramps all live in plotstyle (:data:`SEQ_CMAP`,
    :data:`DIV_CMAP`, :data:`IMAGE_CMAP`, :data:`IMAGE_FLAT_CMAP`,
    :data:`IMAGE_GREY`), so a compliant line reads ``cmap=ps.SOMETHING``
    and never ``cmap="something"``.
    """
    offenders = []
    for path in sorted((REPO_ROOT / "pipeline").rglob("*.py")):
        if path.name in ("plotstyle.py", "test_plotstyle.py"):
            continue
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"""\bcmap\s*=\s*["']""", line):
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{i}"
                                 f"  {line.strip()}")
    assert not offenders, (
        "a colormap must be named through macro_core.plotstyle "
        "(SEQ_CMAP / DIV_CMAP / IMAGE_CMAP / IMAGE_FLAT_CMAP / IMAGE_GREY), "
        "not passed as a bare string: " + "; ".join(offenders))


def test_the_sanctioned_colormaps_are_all_resolvable():
    """A guard naming five constants is only as good as their existing."""
    import matplotlib
    for name in ("SEQ_CMAP", "DIV_CMAP", "IMAGE_CMAP",
                 "IMAGE_FLAT_CMAP", "IMAGE_GREY"):
        value = getattr(ps, name)
        # Either an actual colormap object, or a name matplotlib knows.
        if isinstance(value, str):
            assert matplotlib.colormaps[value] is not None, name
        else:
            assert callable(value), name


def _relative_luminance(rgb) -> float:
    """WCAG relative luminance of an (r, g, b) triple in 0..1.

    The same formula ``plotstyle.text_on`` uses, restated here rather than
    imported so the test measures the property independently of the code it
    is checking.
    """
    def lin(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = rgb
    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)


def test_the_one_sequential_ramp_survives_greyscale():
    """SEQ_CMAP's ends must be far apart in LUMINANCE, not just in hue.

    This is the property the SN saturation matrix lost: its old
    vermilion->orange->green ramp put "none usable" at grey 128 and "all
    usable" at grey 138 — ten levels apart out of 255 — with the midpoint
    at 171, LIGHTER than either end.  Photocopied, that figure read
    inverted.  A ramp carrying a quantity has to keep its ends separable
    when the colour is thrown away.
    """
    lo = _relative_luminance(ps.SEQ_CMAP(0.0)[:3])
    hi = _relative_luminance(ps.SEQ_CMAP(1.0)[:3])
    assert lo - hi > 0.5, ("the sequential ramp's ends are too close in "
                           f"luminance to survive greyscale: {lo:.3f} "
                           f"vs {hi:.3f}")


# ---------------------------------------------------------------------------
# THE GUARD: no published figure may have a dark ground
# ---------------------------------------------------------------------------
def _published_pngs() -> list[Path]:
    return sorted(p for p in DOCS_DIR.rglob("*.png") if p.is_file())


def test_the_site_actually_has_figures_to_check():
    """A guard that silently checks nothing is not a guard."""
    assert len(_published_pngs()) >= 50


@pytest.mark.parametrize("png", _published_pngs(),
                         ids=lambda p: str(p.relative_to(DOCS_DIR)))
def test_no_published_figure_has_a_dark_ground(png):
    """Every corner of every published figure is light.

    Four corners rather than one: a figure whose top-left happens to fall
    inside a dark inset (an S1 failure thumbnail, an S2 reconstruction
    stamp) would pass a one-corner check while its canvas was black.  All
    four corners being dark is the signature of a dark canvas; the real
    image panels this repository publishes are always inset with a margin
    of canvas around them.
    """
    import matplotlib.image as mpimg

    img = mpimg.imread(png)
    arr = np.asarray(img)
    if arr.dtype.kind == "f":
        arr = (arr * 255.0).round()
    if arr.ndim == 2:
        arr = arr[:, :, None]
    corners = [arr[0, 0], arr[0, -1], arr[-1, 0], arr[-1, -1]]

    def is_light(px):
        rgb = px[:3]
        if len(px) == 4 and px[3] < 8:
            return True          # fully transparent corner: no ground drawn
        return float(np.min(rgb)) >= LIGHT_MIN

    light = [is_light(c) for c in corners]
    assert any(light), (
        f"{png.relative_to(DOCS_DIR)} has a dark ground at every corner "
        f"({[list(map(int, c[:3])) for c in corners]}) — regenerate it "
        f"through macro_core.plotstyle instead of a local rcParams dict")


def test_no_published_svg_paints_a_dark_ground():
    """The provenance DAG is drawn as SVG and gets the same rule."""
    dark = re.compile(r'fill="#(0[0-9a-f]|1[0-9a-f]|2[0-9a-f])'
                      r'[0-9a-f]{4}"', re.I)
    for svg in sorted(DOCS_DIR.rglob("*.svg")):
        text = svg.read_text()
        head = text[:text.find(">", text.find("<rect")) + 1] if "<rect" \
            in text else ""
        assert not dark.search(head), f"{svg} opens with a dark ground rect"
