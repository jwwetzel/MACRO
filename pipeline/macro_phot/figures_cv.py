"""CV-S11 — the manuscript figure set, emitted from the products database.

WHY THIS MODULE EXISTS
----------------------
The web reports under ``docs/CV_TimeSeries/`` already obey one law: every
number on a published page is the result of a query executed by the script
that writes the page.  Nothing is typed.  This module applies that same law
to the PAPER.  Each figure in ``CV_TimeSeries/ANALYSIS_STRATEGY.md`` §7 is
built here by a function that opens the products database, asks it for the
rows, and draws them.  There is no hand-made plot anywhere in the
manuscript, and re-running this module after a pipeline change redraws
every panel from whatever the database now says.

WHAT "PUBLICATION QUALITY" MEANS HERE, CONCRETELY
-------------------------------------------------
* **Vector PDF for the manuscript, raster PNG for the web page.**  The PDF
  is what LaTeX includes; the PNG is what the browser shows.  Both come out
  of the same ``Figure`` object in the same call, so they can never drift.
* **AASTeX-friendly widths.**  ``COL_SINGLE`` (3.5 in) fits one column of a
  two-column ApJ page; ``COL_DOUBLE`` (7.1 in) spans both.  A figure drawn
  at any other width is resized by the journal and its fonts stop matching
  the body text, so those two numbers are the only widths used.
* **Colour-blind-safe palette.**  The Okabe--Ito eight-colour set, which is
  distinguishable under deuteranopia, protanopia and tritanopia.  Colour is
  never the ONLY channel: every series that is coloured is also given a
  distinct marker or line style, so the figure survives a greyscale print.
* **Axis labels carry units.**  Every axis label in this module states the
  quantity and its unit.  A bare "amplitude" is not a label.
* **Captions are data, not decoration.**  Each builder returns a caption
  that names what the panel shows AND which database table each series was
  read from.  Those captions are written to ``p5_figure`` and pasted into
  the manuscript by the numbers emitter, so a reader can go from a panel to
  the rows behind it without opening the code.

WHAT CANNOT BE DRAWN, AND WHY THAT IS STATED RATHER THAN FAKED
---------------------------------------------------------------
Three of the strategy's thirteen figures rest on observations that do not
exist.  Each of those is drawn as an HONEST SUBSTITUTE and the substitution
is recorded in ``p5_figure.substitute_reason``, printed in the caption, and
carried into the manuscript:

``fig06``  The cyclotron colour--phase diagram was specified for ST LMi
           *and* VV Pup, per camera.  VV Pup has **1 three-filter
           full-orbit night out of 18**, and its two cameras are fully
           confounded with epoch, so a per-camera colour--phase panel would
           be one night wearing a season's clothes.  EU UMa has **0 of 25**.
           Substitute: ST LMi only, two era panels, with the missing
           targets' night counts printed inside the figure.
``fig07``  The VV Pup / EU UMa folded curves were specified as three-filter
           panels.  They are drawn in the bands that actually have
           full-orbit coverage, and EU UMa's merged 2026 Fast series is
           excluded entirely -- it has five comparison stars, ZERO check
           stars, and no catalogue tie, so nothing from it may be shown
           with the confidence a validated series earns.
``fig09``  The O--C diagram was specified for ST LMi, VV Pup and EU UMa.
           CV-S9's cycle-count stage graded VV Pup and AN UMa **NOT ONE
           FEATURE -- NO O-C** (their accepted edges scatter over 0.13 and
           0.15 in orbital phase, against the 0.05 that would make them one
           timeable feature), and EU UMa has two accepted edges.  Substitute:
           ST LMi's two eras per band, plus the band-offset panel that the
           strategy asked not to be averaged away.
``fig11``  The strategy made this figure conditional: superhump periodogram
           and Kato-style O--C **if** a dense run fell in a superoutburst,
           "else orbital-hump/flickering statistics".  CV-S7 measured that
           no dense run did.  This is therefore the specified fallback, not
           a substitute chosen here.

EVERY OTHER FIGURE IS THE FIGURE THE STRATEGY ASKED FOR.

STRUCTURE OF THIS MODULE
------------------------
Pure functions first -- folding, binding, binning, pairing, formatting --
each of which takes arrays and returns arrays and is exercised directly by
``pipeline/tests/test_figures_cv.py``.  Then the database readers.  Then
one builder per figure, each returning ``(Figure, FigureSpec)``.  The CLI
in ``pipeline/scripts/run_cv_paper.py`` is I/O and bookkeeping only.
"""

from __future__ import annotations

import math
import sqlite3
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt              # noqa: E402
from matplotlib.figure import Figure         # noqa: E402
from matplotlib.lines import Line2D          # noqa: E402

# The macro emitter, for the strings and lists a caption and the prose must
# SHARE rather than keep in step by hand: the scope of the hold-out rule
# (§3.1 and Figure 2's caption), and which colour pairs Figure 6 draws
# (§3.3's tie bars).  numbers_cv imports nothing from here, so this is not
# a cycle.
from . import numbers_cv as _nx              # noqa: E402
from . import final_science as _fs           # noqa: E402

# ===========================================================================
# Page geometry and style
# ===========================================================================

#: AASTeX two-column body: one column is 3.5 in, the full text block 7.1 in.
#: Drawing at any other width means the journal rescales the file and the
#: figure's fonts stop matching the caption's.  These are the only widths.
COL_SINGLE = 3.5
COL_DOUBLE = 7.1

#: Raster resolution for the web copy.  The PDF is vector and unaffected.
PNG_DPI = 200

#: Okabe--Ito: eight hues chosen to stay distinguishable under the three
#: common forms of colour-blindness.  Named, not indexed, so a reader of the
#: code can see which physical thing each colour means.
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

#: Filter -> colour.  Both the 2024 Johnson-ish G/R/I set and the Sloan
#: g/r/i set map to the SAME three hues on purpose: the paper's whole
#: era-seam argument is that these are different bandpasses measuring the
#: same three parts of the spectrum, and a reader comparing two era panels
#: should see blue against blue.  The bandpass difference is stated in the
#: axis label and the caption, never smuggled in as a colour change.
BAND_COLOR = {
    "G": OKABE_ITO["blue"],      "g": OKABE_ITO["blue"],
    "R": OKABE_ITO["vermilion"], "r": OKABE_ITO["vermilion"],
    "I": OKABE_ITO["green"],     "i": OKABE_ITO["green"],
    "z": OKABE_ITO["purple"],    "y": OKABE_ITO["orange"],
}

#: Filter -> marker.  The second channel that makes the figure survive a
#: greyscale print, and the reason colour alone is never load-bearing.
BAND_MARKER = {
    "G": "o", "g": "o",
    "R": "s", "r": "s",
    "I": "^", "i": "^",
    "z": "D", "y": "v",
}

#: Readout mode -> marker.  Camera and epoch are confounded for VV Pup and
#: for the ST LMi era seam, so the camera has to be visible in any figure
#: that spans the seam -- otherwise a reader attributes an instrument step
#: to the star.
MODE_MARKER = {
    "High Gain": "o",
    "High Gain StackPro": "P",
    "1MHz High Sensitivity 16-bit": "s",
    "Mode0": "^",
    "Fast": "D",
}

#: Readout mode -> the abbreviation a cramped axis can carry.  The full
#: strings are in the instrument table; these exist so a three-panel figure
#: does not have to choose between legible labels and legible data.
MODE_SHORT = {
    "High Gain": "High Gain",
    "High Gain StackPro": "HG StackPro",
    "1MHz High Sensitivity 16-bit": "1MHz HS",
    "5MHz High Sensitivity 16-bit": "5MHz HS",
    "3MHz High Sensitivity 16-bit": "3MHz HS",
    "Mode0": "Mode0",
    "Fast": "Fast",
    "Low Gain": "Low Gain",
    "(blank)": "(blank)",
    "HDR": "HDR",
    "None": "(none)",
}

#: Accretion state -> colour.  High and low states are the paper's main
#: categorical split and they get the two most separable hues in the set.
STATE_COLOR = {
    "high": OKABE_ITO["vermilion"],
    "low": OKABE_ITO["blue"],
    "intermediate": OKABE_ITO["orange"],
    "censored": OKABE_ITO["grey"],
    "unclassified": OKABE_ITO["grey"],
    "unknown": "#c9c9c9",
}

#: Target key -> the name a reader knows.
TARGET_LABEL = {
    "stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
    "anuma": "AN UMa", "yzcnc": "YZ Cnc",
}

#: Era id -> the camera/readout label the instrument table uses.  Read from
#: the manifest in principle; pinned here because these seven ids are the
#: only ones this project's series use and the figure legend must be stable
#: across re-runs even when the manifest gains eras for other projects.
ERA_LABEL = {
    6: "High Gain StackPro", 7: "High Gain",
    47: "1MHz High Sensitivity 16-bit", 72: "1MHz High Sensitivity 16-bit",
    76: "Mode0", 78: "Fast", 79: "Fast",
}


def apply_style() -> None:
    """Set the rcParams every figure in this module is drawn under.

    ``pdf.fonttype = 42`` embeds TrueType rather than Type-3, which is what
    lets a journal's production system re-flow and search the text in the
    figure.  Type-3 is the matplotlib default and several publishers reject
    it outright.
    """
    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 8.5,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "axes.linewidth": 0.7,
        "grid.linewidth": 0.4,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "figure.dpi": 110,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    })


# ===========================================================================
# The record each builder returns
# ===========================================================================
@dataclass
class FigureSpec:
    """Everything about a figure except its pixels.

    ``caption`` is the sentence that goes under the figure in the paper and
    on the web page.  ``tables`` is the list of database tables the panels
    were read from, and it is appended to the caption verbatim -- that is
    the promise this project makes about every published figure, and it has
    to be machine-checkable, not a habit.
    """

    fig_id: str                       # "fig01"
    label: str                        # LaTeX \label target
    title: str                        # short human title
    caption: str                      # full caption, no table list
    tables: tuple[str, ...] = ()      # source tables, appended to caption
    width_in: float = COL_DOUBLE
    substitute: bool = False
    substitute_reason: str = ""
    note: str = ""

    @property
    def full_caption(self) -> str:
        """Caption plus the provenance clause every figure here carries.

        The table names are wrapped in ``\\texttt`` with their underscores
        escaped, because this string is pasted straight into a LaTeX
        ``\\caption``: a bare ``cv_frames`` there is a subscript outside
        maths mode and ends the ``tectonic`` run.
        """
        parts = [self.caption.strip()]
        if self.substitute:
            reason = self.substitute_reason.strip().rstrip(".")
            parts.append("SUBSTITUTE FOR THE PLANNED FIGURE: "
                         + reason + ".")
        if self.tables:
            pretty = [f"\\texttt{{{t.replace('_', chr(92) + '_')}}}"
                      for t in self.tables]
            parts.append("Drawn from " + ", ".join(pretty) + ".")
        return " ".join(p for p in parts if p)


# ===========================================================================
# Pure functions -- no database, no matplotlib, no files
# ===========================================================================
def fold_phase(t_bjd, period_d: float, epoch_bjd: float):
    """Orbital phase in [0, 1) for times ``t_bjd``.

    The one arithmetic operation the whole polar analysis rests on, written
    once so that no figure can fold on a slightly different convention than
    another.  ``epoch_bjd`` is the zero point; for YZ Cnc, whose catalogue
    entry has NO epoch, the caller passes an arbitrary constant and the
    caption says so, because a phase measured from an arbitrary zero is
    meaningful WITHIN a run and meaningless between runs.
    """
    t = np.asarray(t_bjd, dtype=float)
    if not (period_d and math.isfinite(period_d) and period_d > 0):
        raise ValueError(f"period must be a positive finite number, got "
                         f"{period_d!r}")
    return np.mod((t - float(epoch_bjd)) / float(period_d), 1.0)


def mad_sigma(x) -> float:
    """Median absolute deviation rescaled to a Gaussian sigma.

    Used everywhere a scatter is quoted, because a single cloud-affected
    point moves an RMS and does not move this.
    """
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    return 1.4826 * float(np.median(np.abs(a - np.median(a))))


def phase_bin(phase, y, n_bins: int = 40, min_count: int = 3):
    """Median-bin ``y`` against ``phase``, returning empty bins as NaN.

    Returns ``(centres, medians, sigmas, counts)``.  A bin with fewer than
    ``min_count`` points comes back NaN rather than as a one-point "median":
    an empty phase bin is a statement about coverage, and filling it in
    would let a reader read sampling as shape.
    """
    ph = np.asarray(phase, dtype=float)
    yy = np.asarray(y, dtype=float)
    ok = np.isfinite(ph) & np.isfinite(yy)
    ph, yy = ph[ok], yy[ok]
    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    med = np.full(n_bins, np.nan)
    sig = np.full(n_bins, np.nan)
    cnt = np.zeros(n_bins, dtype=int)
    if ph.size:
        idx = np.clip(np.digitize(ph, edges) - 1, 0, n_bins - 1)
        for b in range(n_bins):
            sel = yy[idx == b]
            cnt[b] = sel.size
            if sel.size >= min_count:
                med[b] = float(np.median(sel))
                sig[b] = mad_sigma(sel) / math.sqrt(sel.size)
    return centres, med, sig, cnt


def pair_quasi_simultaneous(t_a, m_a, t_b, m_b, max_dt_s: float = 600.0):
    """Match band-B points onto band-A times, keeping only close pairs.

    Returns ``(t_matched, colour, dt_s)`` where ``colour = m_a - m_b``.

    A colour is only a colour if the two magnitudes were measured at nearly
    the same moment.  These stars vary by tenths of a magnitude within one
    92-minute orbit, so a pair separated by half an orbit is not a colour,
    it is two different states subtracted.  ``max_dt_s`` is the gate, and
    ``dt_s`` is returned so a figure can show how well it was met rather
    than merely assert that it was.
    """
    ta = np.asarray(t_a, dtype=float)
    ma = np.asarray(m_a, dtype=float)
    tb = np.asarray(t_b, dtype=float)
    mb = np.asarray(m_b, dtype=float)
    if ta.size == 0 or tb.size == 0:
        return (np.array([]), np.array([]), np.array([]))
    order = np.argsort(tb)
    tb, mb = tb[order], mb[order]
    # Nearest neighbour in time, found by insertion point rather than an
    # O(N*M) distance matrix: some of these nights carry 200 points a band.
    pos = np.clip(np.searchsorted(tb, ta), 1, tb.size - 1)
    left, right = pos - 1, pos
    take_left = np.abs(ta - tb[left]) <= np.abs(ta - tb[right])
    nearest = np.where(take_left, left, right)
    dt_d = np.abs(ta - tb[nearest])
    keep = dt_d * 86400.0 <= float(max_dt_s)
    return (ta[keep], ma[keep] - mb[nearest][keep], dt_d[keep] * 86400.0)


def dpdt_envelope(cycles, sigma_t_s: float):
    """The |dP/dt| a timing campaign of this precision could detect.

    A quadratic term in an O--C curve reaches ``0.5 * (dP/dt) * P * E^2``
    seconds after ``E`` cycles.  Setting that equal to the demonstrated
    per-epoch timing error and solving gives the envelope drawn on the O--C
    figure: below the curve the star's period could be changing and this
    data set would never know.  Returned in SECONDS of O--C so it can be
    plotted directly on the same axis as the residuals.
    """
    e = np.asarray(cycles, dtype=float)
    return float(sigma_t_s) * np.ones_like(e)


def quadratic_oc_seconds(cycles, dpdt_dimensionless: float, period_d: float):
    """O--C in seconds produced by a steady period change ``dP/dt``.

    ``dpdt_dimensionless`` is dP/dt with P in the same units as dt, i.e.
    dimensionless.  This is the curve a reader compares the scatter against
    to see what rate of period change the data would have caught.
    """
    e = np.asarray(cycles, dtype=float)
    return 0.5 * float(dpdt_dimensionless) * float(period_d) * 86400.0 * e**2


def pdot_envelope_seconds(cycles_at, cycles_fit, sigma_fit,
                          dpdt_dimensionless: float, period_d: float):
    """The O--C a steady ``dP/dt`` LEAVES BEHIND after an ephemeris refit.

    A period derivative puts ``0.5 (dP/dt) P E^2`` into an O--C, but a
    constant and a linear term in ``E`` are degenerate with the epoch and
    the period, so a fit that solves for those absorbs part of any real
    curvature.  What a data set can detect is therefore not the raw
    parabola but its residual after that absorption, and drawing the raw
    parabola overstates the signal at one end of the baseline and
    understates it at the other.

    ``cycles_fit`` and ``sigma_fit`` are the epochs and their errors, whose
    inverse-variance weights define the projection -- the same weights
    :func:`run_cv_phase3.pdot_bound` uses.  ``cycles_at`` is where the
    curve is wanted, which for a smooth plot is a dense grid.

    Returns the signed residual curve in seconds; a caller drawing a
    two-sided envelope plots plus and minus its absolute value.
    """
    ea = np.asarray(cycles_at, dtype=float)
    ef = np.asarray(cycles_fit, dtype=float)
    sf = np.asarray(sigma_fit, dtype=float)
    ok = np.isfinite(ef) & np.isfinite(sf) & (sf > 0)
    if int(ok.sum()) < 3:
        return quadratic_oc_seconds(ea - float(np.mean(ef)),
                                    dpdt_dimensionless, period_d)
    ef, sf = ef[ok], sf[ok]
    w = 1.0 / np.square(sf)
    sig_fit = quadratic_oc_seconds(ef, dpdt_dimensionless, period_d)
    design = np.vstack([np.ones_like(ef), ef]).T
    beta = np.linalg.solve(design.T @ (design * w[:, None]),
                           design.T @ (w * sig_fit))
    return (quadratic_oc_seconds(ea, dpdt_dimensionless, period_d)
            - np.vstack([np.ones_like(ea), ea]).T @ beta)


def _envelope_report(oc_s, sigma_s, envelope_s) -> dict:
    """How the epochs actually stand against the drawn envelope.

    Returns the counts and sizes a caption may state, so that the caption
    describes the curve the code drew rather than the curve its author
    remembered.  ``n_outside`` is the number of epochs whose residual
    exceeds the envelope; a caption may only claim containment when it is
    zero, and :func:`fig09_oc` does not claim it at all.
    """
    y = np.abs(np.asarray(oc_s, dtype=float))
    s = np.asarray(sigma_s, dtype=float)
    e = np.abs(np.asarray(envelope_s, dtype=float))
    return {
        "n": int(y.size),
        "n_outside": int(np.sum(y > e)),
        "n_under_error": int(np.sum(e < s)),
        "env_max": float(e.max()) if e.size else float("nan"),
        "sigma_median": float(np.median(s)) if s.size else float("nan"),
    }


def robust_ylim(y, pad_frac: float = 0.08, k: float = 6.0):
    """Axis limits set by the bulk of the data, not by its worst point.

    Returns ``(lo, hi)`` spanning the median plus/minus ``k`` robust sigmas,
    widened by ``pad_frac``, and clipped to the actual data range so the
    axis never claims space where nothing was measured.
    """
    a = np.asarray(y, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return (0.0, 1.0)
    if a.size < 4:
        lo, hi = float(np.min(a)), float(np.max(a))
    else:
        med = float(np.median(a))
        s = mad_sigma(a)
        if not math.isfinite(s) or s <= 0:
            s = float(np.std(a)) or 1.0
        lo = max(float(np.min(a)), med - k * s)
        hi = min(float(np.max(a)), med + k * s)
    if hi <= lo:
        lo, hi = lo - 0.5, hi + 0.5
    pad = pad_frac * (hi - lo)
    return (lo - pad, hi + pad)


def series_parts(series_key: str) -> tuple[str, int, str]:
    """``stlmi|e76|g`` -> ``("stlmi", 76, "g")``."""
    tgt, era, filt = series_key.split("|")
    return tgt, int(era.lstrip("e")), filt


def series_label(series_key: str) -> str:
    """``stlmi|e76|g`` -> ``ST LMi Mode0 g``, the label a legend carries."""
    tgt, era, filt = series_parts(series_key)
    return f"{TARGET_LABEL.get(tgt, tgt)} {ERA_LABEL.get(era, f'e{era}')} {filt}"


def night_to_ordinal(nights: Sequence[str]) -> np.ndarray:
    """ISO night strings -> days since 2024-01-01, for a linear date axis.

    A plain ``datetime`` axis would work, but this project's figures are
    also compared against BJD-indexed panels, and a single float day number
    keeps the two commensurate without a timezone anywhere in the code.
    """
    import datetime as _dt
    base = _dt.date(2024, 1, 1).toordinal()
    out = np.full(len(nights), np.nan)
    for i, n in enumerate(nights):
        try:
            y, m, d = (int(v) for v in str(n).split("-"))
            out[i] = _dt.date(y, m, d).toordinal() - base
        except Exception:                                   # noqa: BLE001
            continue
    return out


def year_ticks(day_lo: float, day_hi: float, max_ticks: int = 10):
    """Tick positions and labels for the day-since-2024-01-01 axis.

    The spacing ADAPTS to the span.  A panel showing four months of one
    season wants month labels; a panel showing thirty years of AAVSO
    coverage wants decades.  Choosing one fixed spacing produces either an
    unreadable comb of overlapping labels or a bare axis, and both have
    happened here, so the choice is made from the span and capped at
    ``max_ticks`` labels.
    """
    import datetime as _dt
    base = _dt.date(2024, 1, 1).toordinal()
    span_y = max(1e-6, (day_hi - day_lo) / 365.25)
    pos, lab = [], []
    if span_y <= 1.6:
        # Months, thinned so at most max_ticks land on the axis.
        step = max(1, int(math.ceil(span_y * 12.0 / max_ticks)))
        y0 = _dt.date.fromordinal(base + int(math.floor(day_lo))).year
        y1 = _dt.date.fromordinal(base + int(math.ceil(day_hi))).year
        k = 0
        for y in range(y0, y1 + 1):
            for m in range(1, 13):
                d = _dt.date(y, m, 1).toordinal() - base
                if day_lo - 1 <= d <= day_hi + 1:
                    if k % step == 0:
                        pos.append(d)
                        lab.append(f"{y}-{m:02d}")
                    k += 1
        return pos, lab
    # Whole years, on a 1/2/5/10 ladder so the labels stay round numbers.
    for step in (1, 2, 5, 10, 20, 50):
        if span_y / step <= max_ticks:
            break
    y0 = _dt.date.fromordinal(base + int(math.floor(day_lo))).year
    y1 = _dt.date.fromordinal(base + int(math.ceil(day_hi))).year
    for y in range(y0 - (y0 % step), y1 + 1, step):
        d = _dt.date(max(y, 1), 1, 1).toordinal() - base
        if day_lo - 1 <= d <= day_hi + 1:
            pos.append(d)
            lab.append(str(y))
    return pos, lab


def normalise_state(state) -> str:
    """``'HIGH'``, ``'high'``, ``None`` -> a key ``STATE_COLOR`` knows.

    The state classifier writes upper case and the palette is keyed in
    lower; without this the whole of Figures 5, 7 and 8 draws in the
    unknown-state grey, which silently deletes the paper's main categorical
    result from its own figures.  That happened once; this function is why
    it cannot happen again, and ``test_figures_cv.py`` pins it.
    """
    s = str(state or "unknown").strip().lower()
    return s if s in STATE_COLOR else ("unknown" if s != "" else "unknown")


# ===========================================================================
# Database readers.  Every one of these is the ONLY route to its numbers.
# ===========================================================================
def connect_ro(path: Path) -> sqlite3.Connection:
    """Read-only connection.  These databases are products, not scratch."""
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def table_exists(con: sqlite3.Connection, name: str) -> bool:
    return bool(con.execute(
        "SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?",
        (name,)).fetchone()[0])


def read_ephemeris(con) -> dict[str, dict]:
    """Published ephemerides, keyed by target.

    Read rather than typed so that the period a figure folds on is provably
    the period CV-S9 verified against and the manuscript cites.  YZ Cnc's
    ``epoch_bjd`` is NULL here and stays NULL: the caller must decide, in
    the open, what to do about a star with no published epoch.
    """
    return {r["target_key"]: dict(r)
            for r in con.execute("SELECT * FROM p3_ephemeris")}


def read_target_points(con, series_key: str, night: Optional[str] = None):
    """Target light curve for one series, cloud-vetoed frames removed.

    ``cal_mag`` throughout -- the catalogue-tied magnitude -- exactly as
    Phase 3 and CV-S10 use, and rows without one are DROPPED rather than
    falling back to the instrumental magnitude.  Mixing the two inside one
    curve would put a zero-point step in the middle of a light curve that
    the figure then shows as variability.
    """
    sql = """
        SELECT l.frame_id, l.bjd_tdb, l.cal_mag, l.inst_mag_err,
               f.night, f.airmass, f.exptime, f.readoutm
        FROM cv_lightcurve l
        JOIN cv_frames f ON f.frame_id = l.frame_id
                        AND f.series_key = l.series_key
        LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                                  AND c.series_key = l.series_key
        WHERE l.series_key = ? AND l.role = 'target'
          AND l.cal_mag IS NOT NULL AND l.saturated = 0
          AND COALESCE(c.vetoed, 0) = 0
    """
    args: list = [series_key]
    if night is not None:
        sql += " AND f.night = ?"
        args.append(night)
    sql += " ORDER BY l.bjd_tdb"
    rows = con.execute(sql, args).fetchall()
    return {
        "bjd": np.array([r["bjd_tdb"] for r in rows], dtype=float),
        "mag": np.array([r["cal_mag"] for r in rows], dtype=float),
        "err": np.array([(r["inst_mag_err"] if r["inst_mag_err"] else np.nan)
                         for r in rows], dtype=float),
        "night": [r["night"] for r in rows],
    }


def read_rows(con, sql: str, args: Sequence = ()) -> list[dict]:
    """Every ad-hoc query in this module goes through here, as dicts."""
    return [dict(r) for r in con.execute(sql, tuple(args))]


# ===========================================================================
# Figure builders.  One function per strategy figure.
# ===========================================================================
def _empty_panel(ax, message: str) -> None:
    """What a panel says when the observations behind it do not exist.

    Never a blank box and never a plausible-looking curve drawn from
    something else: the sentence that says which rows were missing.
    """
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center",
            va="center", fontsize=6.5, color=OKABE_ITO["grey"], wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])


def fig01_coverage(cv, ext_targets=("stlmi", "vvpup", "euuma",
                                    "anuma", "yzcnc")):
    """Figure 1 -- coverage and cadence map, all five targets."""
    fig, axes = plt.subplots(2, 1, figsize=(COL_DOUBLE, 5.0), sharex=True,
                             gridspec_kw={"height_ratios": [2.1, 1.0]})
    ax, ax2 = axes

    rows = read_rows(cv, """
        SELECT target_key, era_id, filter, night, count(*) n,
               (max(bjd_tdb) - min(bjd_tdb)) * 24.0 span_h
        FROM cv_frames
        WHERE status IS NULL OR status NOT IN ('excluded')
        GROUP BY target_key, era_id, filter, night
    """)
    ext = read_rows(cv, """
        SELECT target, source, utc_night FROM cv_ext_nightly
        WHERE independent = 1
    """) if table_exists(cv, "cv_ext_nightly") else []

    ycat = [t for t in ext_targets]
    ypos = {t: i for i, t in enumerate(ycat)}

    # The window is the RLMT CAMPAIGN, not the survey baseline.  AAVSO
    # reaches back to 1980 for AN UMa; drawn on one axis with the RLMT
    # nights, that baseline compresses every observation this paper reports
    # into two millimetres of ink.  The decades-long record is Figure 8's
    # subject and it is shown there; this panel is about cadence.
    rlmt_days = night_to_ordinal([r["night"] for r in rows])
    rlmt_days = rlmt_days[np.isfinite(rlmt_days)]
    if rlmt_days.size:
        pad = 0.06 * max(30.0, rlmt_days.max() - rlmt_days.min())
        x_lo, x_hi = rlmt_days.min() - pad, rlmt_days.max() + pad
    else:
        x_lo, x_hi = -30.0, 30.0

    # Survey epochs first, underneath, in pale grey: the point of the panel
    # is what RLMT adds ON TOP of the sparse survey record.
    n_ext_in, n_ext_all = 0, 0
    for r in ext:
        y = ypos.get(r["target"])
        if y is None:
            continue
        n_ext_all += 1
        x = night_to_ordinal([r["utc_night"]])[0]
        if math.isfinite(x) and x_lo <= x <= x_hi:
            n_ext_in += 1
            ax.plot([x], [y + 0.30], marker="|", ms=3.5, lw=0,
                    color="#c4c4c4", zorder=1)

    seen_band, seen_mode = set(), set()
    for r in rows:
        y = ypos.get(r["target_key"])
        if y is None:
            continue
        x = night_to_ordinal([r["night"]])[0]
        if not math.isfinite(x):
            continue
        f = r["filter"]
        mode = ERA_LABEL.get(r["era_id"], "")
        ax.plot([x], [y], marker=MODE_MARKER.get(mode, "o"),
                ms=3.4, lw=0, mfc=BAND_COLOR.get(f, OKABE_ITO["grey"]),
                mec="none", alpha=0.85, zorder=3)
        seen_band.add(f)
        seen_mode.add(mode)

    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([TARGET_LABEL[t] for t in ycat])
    # Headroom above the top target so the two legends sit on empty axis
    # rather than on YZ Cnc's nights.
    ax.set_ylim(-0.7, len(ycat) + 0.9)
    ax.grid(axis="x", color="#eeeeee", zorder=0)
    ax.set_ylabel("target")

    band_handles = [Line2D([], [], marker="o", lw=0, ms=4,
                           mfc=BAND_COLOR[f], mec="none", label=f)
                    for f in ("G", "g", "R", "r", "I", "i", "z", "y")
                    if f in seen_band]
    mode_handles = [Line2D([], [], marker=MODE_MARKER[m], lw=0, ms=4,
                           mfc="none", mec="k", mew=0.7, label=m)
                    for m in MODE_MARKER if m in seen_mode]
    grey = Line2D([], [], marker="|", lw=0, ms=5, color="#bbbbbb",
                  label="survey epoch (ZTF/ASAS-SN/AAVSO)")
    leg1 = ax.legend(handles=band_handles, loc="upper left", ncol=4,
                     title="filter", title_fontsize=6.5,
                     bbox_to_anchor=(0.005, 0.99))
    ax.add_artist(leg1)
    ax.legend(handles=mode_handles + [grey], loc="upper right", ncol=1,
              title="readout mode", title_fontsize=6.5,
              bbox_to_anchor=(0.998, 0.99))

    # Lower panel: orbital cycles covered per night.  This is the number
    # that says what RLMT adds that a survey cannot: a survey epoch is a
    # point, and these are continuous runs several cycles long.
    per = read_rows(cv, "SELECT target_key, period_d FROM p3_ephemeris")
    porb = {r["target_key"]: r["period_d"] for r in per}
    for r in rows:
        t = r["target_key"]
        p = porb.get(t)
        if not p or r["span_h"] is None:
            continue
        x = night_to_ordinal([r["night"]])[0]
        cyc = float(r["span_h"]) / 24.0 / float(p)
        if math.isfinite(x) and cyc > 0:
            ax2.plot([x], [cyc], marker=BAND_MARKER.get(r["filter"], "o"),
                     ms=2.8, lw=0, mfc=BAND_COLOR.get(r["filter"],
                                                      OKABE_ITO["grey"]),
                     mec="none", alpha=0.8)
    ax2.axhline(1.0, color=OKABE_ITO["black"], lw=0.8, ls="--")
    ax2.text(0.005, 1.05, "one full orbit", transform=ax2.get_yaxis_transform(),
             fontsize=6, va="bottom", color=OKABE_ITO["black"])
    ax2.set_yscale("log")
    ax2.set_ylabel("orbital cycles covered\nper night per filter")
    ax2.grid(color="#eeeeee")

    ax.set_xlim(x_lo, x_hi)
    ax2.set_xlim(x_lo, x_hi)
    pos, lab = year_ticks(x_lo, x_hi, max_ticks=9)
    ax2.set_xticks(pos)
    ax2.set_xticklabels(lab, rotation=45, ha="right")
    ax2.set_xlabel("UTC night")
    ax.text(0.5, -0.035,
            f"axis limited to the RLMT campaign; {n_ext_all - n_ext_in:,} "
            f"earlier survey epochs fall outside it and appear in Fig. 8",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.6,
            color=OKABE_ITO["grey"])

    pmin = [1440.0 * float(r["period_d"]) for r in per
            if r["period_d"] and r["target_key"] in ypos]
    p_lo, p_hi = (min(pmin), max(pmin)) if pmin else (float("nan"),) * 2

    spec = FigureSpec(
        fig_id="fig01", label="fig:coverage",
        title="Coverage and cadence map",
        caption=(
            "(a) Every RLMT night of every target, coloured by filter and "
            "marked by readout mode, with independent survey epochs "
            "(ZTF, ASAS-SN, AAVSO) underlaid as pale ticks. The camera "
            "changes are visible as marker changes because camera and epoch "
            "are confounded for VV Pup and across the 2024-05 ST LMi seam, "
            "and an instrument step must never be readable as a stellar "
            "one. (b) Orbital cycles covered per night per filter, on a "
            "logarithmic axis, against the one-cycle line: this is what "
            "these data add to the sparse survey record, which samples "
            f"these {p_lo:.0f}--{p_hi:.0f} minute binaries one point at a "
            "time. This is a "
            "census of what was OBSERVED, so EU~UMa's 2026 Fast-mode "
            "nights appear here even though no measurement in this paper "
            "uses them (Section~\\ref{sec:vvpupeuuma}); a coverage map "
            "that omitted observed nights would be a different claim."),
        tables=("cv_frames", "cv_ext_nightly", "p3_ephemeris"),
        width_in=COL_DOUBLE)
    return fig, spec


def fig02_rms_vs_mag(ch, cv, man):
    """Figure 2 -- RMS against magnitude, per camera, model overplotted.

    ONLY series that carry a catalogue tie appear.  The abscissa is a
    catalogue-tied magnitude, so a series with no tie has no position on
    it; and the annotated systematic floor is a precision statement, which
    §3.1 allows only on held-out check stars of a converged, tied solve.
    EU UMa's 2026 Fast-mode block is neither -- five comparison stars, no
    check stars, no tie and no target detection -- and an earlier revision
    of this figure gave it a panel of its own, whose 1--2 mmag floor was
    the smallest number in the figure and rested on nothing the paper is
    allowed to rest a precision claim on.
    """
    tied = {r["series_key"] for r in read_rows(cv, """
        SELECT series_key FROM cv_cattie
        WHERE is_primary = 1 AND verdict LIKE 'TIED%'
    """)}
    series = [r for r in read_rows(ch, """
        SELECT series_key, target_key, era_id, filter, readoutm, exptime,
               floor_nom, k_nom, floor_lo, floor_hi, k_lo, k_hi,
               target_mag, prec_at_target, check_rms_med, n_stars,
               best_star_mag, best_star_rms
        FROM ch_noise_series WHERE n_stars >= 20
        ORDER BY era_id, filter
    """) if r["series_key"] in tied]
    modes: dict[str, list[dict]] = {}
    for r in series:
        modes.setdefault(r["readoutm"] or "unknown", []).append(r)
    order = [m for m in ("High Gain", "1MHz High Sensitivity 16-bit",
                         "Mode0", "Fast") if m in modes]
    n = max(1, len(order))
    ncol = min(2, n)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(COL_DOUBLE, 2.5 * nrow),
                             squeeze=False, sharex=True, sharey=True)
    flat = [a for row in axes for a in row]

    for ax, mode in zip(flat, order):
        floors_lo, floors_hi = [], []
        for r in modes[mode]:
            stars = read_rows(ch, """
                SELECT mean_mag, rms, pred_lo, pred_nom, pred_hi
                FROM ch_noise_stars WHERE series_key = ?
                  AND rms IS NOT NULL AND mean_mag IS NOT NULL
            """, (r["series_key"],))
            if not stars:
                continue
            f = r["filter"]
            col = BAND_COLOR.get(f, OKABE_ITO["grey"])
            m = np.array([s["mean_mag"] for s in stars])
            v = np.array([s["rms"] for s in stars]) * 1000.0
            ax.plot(m, v, marker=BAND_MARKER.get(f, "o"), lw=0, ms=1.5,
                    mfc=col, mec="none", alpha=0.25, zorder=2)
            # THE MODEL, drawn as the quantity it actually predicts: the
            # per-star photon-and-read term added IN QUADRATURE to the
            # series' own measured systematic floor.  Plotting the photon
            # term alone -- the first version of this figure did -- puts
            # the model an order of magnitude below every point and invites
            # the reader to conclude the photometry is broken.  It is not:
            # the floor is real, it is measured, and it is the thing the
            # paper's precision claims are about.
            o = np.argsort(m)
            fl_lo = float(r["floor_lo"] or r["floor_nom"] or 0.0)
            fl_hi = float(r["floor_hi"] or r["floor_nom"] or 0.0)
            plo = np.array([s["pred_lo"] for s in stars], dtype=float)[o]
            phi = np.array([s["pred_hi"] for s in stars], dtype=float)[o]
            good = np.isfinite(plo) & np.isfinite(phi)
            if good.sum() > 3:
                ax.fill_between(
                    m[o][good],
                    1000.0 * np.hypot(plo[good], fl_lo),
                    1000.0 * np.hypot(phi[good], fl_hi),
                    color=col, alpha=0.30, lw=0, zorder=3)
            floors_lo.append(1000.0 * fl_lo)
            floors_hi.append(1000.0 * fl_hi)
            # The star that actually did best, and the precision reached at
            # the CV's own magnitude: two points, not a claim about a curve.
            if r["best_star_mag"] and r["best_star_rms"]:
                ax.plot([r["best_star_mag"]], [1000.0 * r["best_star_rms"]],
                        marker="*", ms=6, mfc=col, mec="k", mew=0.4,
                        zorder=6)
            if r["target_mag"] and r["prec_at_target"]:
                ax.plot([r["target_mag"]], [1000.0 * r["prec_at_target"]],
                        marker="P", ms=4.5, mfc=col, mec="k", mew=0.4,
                        zorder=6)
        if floors_lo:
            ax.axhspan(min(floors_lo), max(floors_hi), color="#000000",
                       alpha=0.055, lw=0, zorder=1)
            ax.text(0.985, 0.05,
                    f"measured systematic floor "
                    f"{min(floors_lo):.0f}--{max(floors_hi):.0f} mmag",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=5.6, color=OKABE_ITO["grey"])
        ax.set_yscale("log")
        ax.set_ylim(0.5, 400)
        ax.set_title(mode, loc="left")
        ax.grid(color="#eeeeee")
    for ax in flat[len(order):]:
        ax.set_visible(False)
    for row in axes:
        row[0].set_ylabel("per-point scatter (mmag)")
    for ax in axes[-1]:
        ax.set_xlabel("catalogue-tied magnitude (mag)")
    handles = [Line2D([], [], marker=BAND_MARKER.get(f, "o"), lw=0, ms=4,
                      mfc=BAND_COLOR[f], mec="none", label=f)
               for f in ("G", "g", "R", "r", "I", "i", "z")
               if any(s["filter"] == f for s in series)]
    handles += [
        Line2D([], [], marker="*", lw=0, ms=6, mfc="w", mec="k",
               label="best star"),
        Line2D([], [], marker="P", lw=0, ms=4.5, mfc="w", mec="k",
               label="at the CV's magnitude")]
    flat[0].legend(handles=handles, loc="upper left", ncol=3, fontsize=5.6)

    # The gain bracket the noise band is drawn from, read rather than typed.
    gb = read_rows(man, """SELECT quantity, value FROM detector_params
                           WHERE era_group='High Gain' AND quantity IN
                           ('gain_lower_bound_e_per_adu',
                            'gain_upper_bound_e_per_adu')""")
    gv = {r["quantity"]: float(r["value"]) for r in gb}
    g_lo = gv.get("gain_lower_bound_e_per_adu", float("nan"))
    g_hi = gv.get("gain_upper_bound_e_per_adu", float("nan"))

    spec = FigureSpec(
        fig_id="fig02", label="fig:rmsmag",
        title="Per-point scatter against magnitude, per camera",
        caption=(
            "Measured scatter of every ensemble and check star against its "
            "catalogue-tied magnitude, one panel per readout mode. The "
            "shaded curve is the noise model: the per-star photon and read "
            "term added in quadrature to that series' own MEASURED "
            f"systematic floor, drawn as a band rather than a line because "
            f"the detector gain is bracketed at {g_lo:.2f}--{g_hi:.2f} "
            "e$^{-}$/ADU and a "
            "single curve would claim a calibration this instrument does "
            "not yet have. The horizontal grey band is that floor, which is "
            "what limits every bright star. Stars mark the best-performing "
            "star in each series and crosses the precision reached at the "
            "CV's own magnitude, which is where the model has to be "
            "believed. The scatter of individual field stars well above the "
            "band is the field's own variability and blending, not the "
            "instrument's. "
            # THE SCOPE OF THE HOLD-OUT RULE, FROM THE SAME STRING §3.1
            # USES.  This caption used to assert that a precision statement
            # "is permitted only on the held-out check stars of a tied
            # solve" -- a rule §3.1 had already been rewritten to retract,
            # on the panel whose crosses are the counter-example, and false
            # of the annotated floor as well, which is fitted over
            # comparison AND check stars exactly as the crosses are.
            + _nx.PRECISION_SCOPE_CLAUSE + ". "
            "Only catalogue-tied series appear, because the abscissa is a "
            "tied magnitude: EU~UMa's untied 2026 "
            "Fast-mode block therefore has no panel here, and the Fast "
            "readout mode no measured floor "
            "(Section~\\ref{sec:vvpupeuuma})."),
        tables=("ch_noise_stars", "ch_noise_series", "cv_cattie",
                "detector_params"),
        width_in=COL_DOUBLE)
    return fig, spec


def fig03_linearity(man):
    """Figure 3 -- linearity ladders and the adopted saturation vetoes."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_DOUBLE, 2.8),
                             gridspec_kw={"width_ratios": [1.45, 1.0],
                                          "wspace": 0.30})
    ax, ax2 = axes

    ladders = read_rows(man, """
        SELECT ladder_id, mode, night, n_rungs FROM s2_linearity_ladders
        WHERE n_rungs >= 3
    """)
    ceil = {r["mode"]: dict(r) for r in read_rows(
        man, "SELECT * FROM s2_ceiling_modes")}
    mode_color = {
        "High Gain": OKABE_ITO["blue"],
        "High Gain StackPro": OKABE_ITO["purple"],
        "1MHz High Sensitivity 16-bit": OKABE_ITO["green"],
        "Mode0": OKABE_ITO["vermilion"],
        "Fast": OKABE_ITO["orange"],
    }
    # One point per RUNG, plus a binned median curve per mode.  Joining
    # the rungs of each ladder -- the first version did -- draws sixty-six
    # zigzag polylines on top of each other and hides the only thing the
    # panel is for: where each mode stops being linear.
    per_mode: dict[str, list[tuple[float, float]]] = {}
    for L in ladders:
        mode = L["mode"]
        if mode not in mode_color:
            continue
        for r in read_rows(man, """
                SELECT peak_med, resid_pct FROM s2_linearity_rungs
                WHERE ladder_id = ?""", (L["ladder_id"],)):
            pk, rs = r["peak_med"], r["resid_pct"]
            if pk and rs is not None and math.isfinite(float(pk)) \
                    and math.isfinite(float(rs)) and float(pk) > 0:
                per_mode.setdefault(mode, []).append((float(pk), float(rs)))
    for mode, pts in per_mode.items():
        col = mode_color[mode]
        pk = np.array([q[0] for q in pts])
        rs = np.array([q[1] for q in pts])
        ax.plot(pk, np.abs(rs), marker=MODE_MARKER.get(mode, "o"), lw=0,
                ms=1.9, mfc=col, mec="none", alpha=0.30, zorder=2)
        # Binned median of |departure|: the mode's typical behaviour, which
        # is what a veto threshold is chosen against.
        edges = np.logspace(math.log10(max(pk.min(), 1.0)),
                            math.log10(pk.max()), 8)
        cx, cy = [], []
        for a_, b_ in zip(edges[:-1], edges[1:]):
            sel = np.abs(rs[(pk >= a_) & (pk < b_)])
            if sel.size >= 5:
                cx.append(math.sqrt(a_ * b_))
                cy.append(float(np.median(sel)))
        if len(cx) >= 3:
            ax.plot(cx, cy, "-", lw=1.4, color=col, zorder=4, label=mode)
    # Each mode's adopted saturation veto, in that mode's colour.  These are
    # the thresholds the photometry actually applies, so they belong on the
    # panel that justifies them.
    for mode, col in mode_color.items():
        c = ceil.get(mode)
        if c and c.get("veto_adu"):
            ax.axvline(float(c["veto_adu"]), lw=0.9, ls="--", color=col,
                       alpha=0.85, zorder=3)
    ax.axhline(1.0, color="k", lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.02, 2e3)
    ax.set_xlabel("median peak pixel level (ADU)")
    ax.set_ylabel("|departure| from a linear\nflux--exposure law (per cent)")
    ax.grid(color="#eeeeee")
    ax.legend(loc="upper left", fontsize=5.6)

    # Right panel: the three ADU numbers per mode that the pipeline uses.
    labels, hard, clip, veto, cols = [], [], [], [], []
    for mode, col in mode_color.items():
        c = ceil.get(mode)
        if not c or not c.get("clip_adu"):
            continue
        labels.append(MODE_SHORT.get(mode, mode))
        hard.append(float(c["hard_max_adu"]))
        clip.append(float(c["clip_adu"]))
        veto.append(float(c["veto_adu"]))
        cols.append(col)
    y = np.arange(len(labels))
    for i in range(len(labels)):
        ax2.plot([veto[i], clip[i]], [y[i], y[i]], lw=3.0, color=cols[i],
                 alpha=0.35, solid_capstyle="butt")
        ax2.plot([clip[i]], [y[i]], marker="|", ms=9, color=cols[i], mew=1.4)
        ax2.plot([veto[i]], [y[i]], marker="<", ms=4.5, color=cols[i],
                 mec="none")
        ax2.plot([hard[i]], [y[i]], marker="x", ms=4, color="k", mew=0.8)
    ax2.set_yticks(y)
    ax2.set_yticklabels(labels, fontsize=6.5)
    ax2.set_ylim(-0.6, len(labels) - 0.4)
    ax2.set_xscale("log")
    ax2.set_xlim(2e3, 1.2e5)
    ax2.set_xlabel("pixel level (ADU)")
    ax2.grid(axis="x", color="#eeeeee")
    ax2.legend(handles=[
        Line2D([], [], marker="|", lw=0, ms=8, color="k",
               label="measured clip (ceiling)"),
        Line2D([], [], marker="<", lw=0, ms=5, color="k",
               label="adopted veto (0.92 ceiling)"),
        Line2D([], [], marker="x", lw=0, ms=5, color="k",
               label="hard observed maximum")],
        loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1,
        fontsize=5.6)

    # The two numbers the caption's punchline rests on, read not typed.
    _dp = {(r["era_group"], r["quantity"]): float(r["value"])
           for r in read_rows(man, "SELECT era_group, quantity, value FROM "
                                   "detector_params")}
    hg_bits = _dp.get(("High Gain", "adc_bits"), float("nan"))
    hg_clip = _dp.get(("High Gain", "ceiling_adu"), float("nan"))
    # Through the shared helper, and as a RANGE.  This caption used to
    # print the minimum ratio alone as "a factor of 16 below every 16-bit
    # mode" while §2.1 printed the maximum as "nearly twenty" -- two
    # numbers for one comparison, and the smaller of them indistinguishable
    # from the nominal bit-depth ratio it is not.
    r_lo, r_hi = _nx.dynamic_range_ratios(_dp)

    spec = FigureSpec(
        fig_id="fig03", label="fig:linearity",
        title="Linearity ladders and adopted saturation vetoes",
        caption=(
            "(a) Twilight and dome exposure ladders per readout mode: the "
            "absolute departure of the measured flux from a linear "
            "flux--exposure law against the median peak pixel level. Faint "
            "points are individual rungs; heavy curves are the binned "
            "median per mode, which is what a veto threshold is chosen "
            "against. Dashed vertical lines are each mode's adopted "
            "saturation veto, the dotted horizontal line one per cent. (b) The three levels the pipeline uses "
            "per mode: the hard observed maximum, the measured clip that "
            f"defines the ceiling, and the veto, which is 0.92 of the "
            f"ceiling rounded down to the nearest 100 ADU. "
            f"High Gain's {hg_bits:.0f}-bit ceiling of {hg_clip:,.0f} ADU "
            f"sits a MEASURED factor of {r_lo:.1f}--{r_hi:.1f} below the "
            "16-bit modes --- the spread is real, because their ceilings "
            "are measured pileup clips and not the same number --- which "
            "is why its vetoes cannot be shared with the Sloan-era data. "
            "The bit depths would suggest a flat factor of 16; that is a "
            "nominal ratio and not one of this figure's measurements."),
        tables=("s2_linearity_ladders", "s2_linearity_rungs",
                "s2_ceiling_modes"),
        width_in=COL_DOUBLE)
    return fig, spec


def fig04_periodograms(cv, ch, picks=None):
    """Figure 4 -- spectral windows and periodograms with the alias caveat."""
    if picks is None:
        picks = ["stlmi|e76|g", "anuma|e76|g", "vvpup|e76|g",
                 "euuma|e76|g", "yzcnc|e7|G"]
    picks = [p for p in picks if read_rows(
        cv, "SELECT 1 FROM p3_period WHERE series_key=? AND status='ok'",
        (p,))]
    n = max(1, len(picks))
    fig, axes = plt.subplots(n, 2, figsize=(COL_DOUBLE, 1.25 * n + 0.6),
                             squeeze=False, sharex="col")

    eph = read_ephemeris(cv)
    for i, sk in enumerate(picks):
        tgt, era, filt = series_parts(sk)
        per = read_rows(cv, "SELECT * FROM p3_period WHERE series_key=?",
                        (sk,))[0]
        col = BAND_COLOR.get(filt, OKABE_ITO["grey"])

        # Left: the window function -- what the sampling alone would produce
        # against a constant star.  It is the reason no period here is a
        # determination.
        axw = axes[i][0]
        win = read_rows(cv, """SELECT freq_cd, value FROM p3_pgram
                               WHERE series_key=? AND panel='window'
                               ORDER BY freq_cd""", (sk,))
        if win:
            fw = np.array([r["freq_cd"] for r in win])
            pw = np.array([r["value"] for r in win])
            axw.plot(fw, pw, lw=0.7, color=OKABE_ITO["black"])
            axw.fill_between(fw, 0, pw, color=OKABE_ITO["grey"], alpha=0.25,
                             lw=0)
            for k in (-1, 1):
                axw.axvline(k, lw=0.6, ls=":", color=OKABE_ITO["vermilion"])
        else:
            _empty_panel(axw, "no window computed")
        axw.set_ylabel(f"{TARGET_LABEL[tgt]} {filt}\nwindow power",
                       fontsize=6.2)
        axw.set_ylim(0, 1.05)

        # Right: the LS periodogram over the survey band, with the PUBLISHED
        # frequency marked and its +/-1 c/d alias pair.  The alias fraction
        # is printed, because it is the number that forbids a claim.
        axp = axes[i][1]
        pg = read_rows(cv, """SELECT freq_cd, value FROM p3_pgram
                              WHERE series_key=? AND panel='survey'
                                AND kind='ls' ORDER BY freq_cd""", (sk,))
        if pg:
            fp = np.array([r["freq_cd"] for r in pg])
            pp = np.array([r["value"] for r in pg])
            axp.plot(fp, pp, lw=0.6, color=col)
            f_pub = 1.0 / float(eph[tgt]["period_d"])
            axp.axvline(f_pub, lw=0.7, color=OKABE_ITO["black"])
            for k in (-1, 1):
                axp.axvline(f_pub + k, lw=0.6, ls="--",
                            color=OKABE_ITO["vermilion"])
            af = per["alias_frac_max"]
            if af is not None:
                axp.text(0.985, 0.90, f"alias power {float(af):.2f}",
                         transform=axp.transAxes, ha="right", va="top",
                         fontsize=6, color=OKABE_ITO["vermilion"])
            axp.text(0.985, 0.66,
                     f"{per['constraint_class'] or '—'} / "
                     f"{per['family_code'] or '—'}",
                     transform=axp.transAxes, ha="right", va="top",
                     fontsize=6, color=OKABE_ITO["grey"])
        else:
            _empty_panel(axp, "no periodogram computed")
        axp.set_ylabel("LS power", fontsize=6.2)
        for a in (axw, axp):
            a.grid(color="#f2f2f2")
            a.tick_params(labelsize=6)

    axes[-1][0].set_xlabel("frequency offset from peak (cycles d$^{-1}$)")
    axes[-1][1].set_xlabel("frequency (cycles d$^{-1}$)")
    axes[0][0].set_title("spectral window", fontsize=7, loc="left")
    axes[0][1].set_title("Lomb--Scargle periodogram", fontsize=7, loc="left")

    # The caption's alias figure, computed in the same call that drew the
    # panels.  Typing it here -- an earlier revision said 0.92 -- lets the
    # caption understate the paper's own worst case while §4.1 states it
    # correctly, so the figure and the text disagree about the number that
    # justifies refusing every period claim.
    af = [r["alias_frac_max"] for r in read_rows(
        cv, "SELECT alias_frac_max FROM p3_period WHERE status='ok' "
            "AND alias_frac_max IS NOT NULL")]
    a_lo, a_hi = (min(af), max(af)) if af else (float("nan"), float("nan"))
    n_ok = len(af)

    spec = FigureSpec(
        fig_id="fig04", label="fig:periodograms",
        title="Spectral windows and periodograms",
        caption=(
            "For one representative series per target: (left) the spectral "
            "window, the periodogram the sampling alone would produce "
            "against a constant star, with the $\\pm1$ cycle d$^{-1}$ "
            "positions dotted; (right) the Lomb--Scargle periodogram over "
            "the survey band, with the published orbital frequency solid "
            "and its $\\pm1$ cycle d$^{-1}$ aliases dashed. The annotated "
            "alias power is the fraction of the peak carried by the "
            f"strongest alias. Over the {n_ok} series searched it spans "
            f"{a_lo:.2f}--{a_hi:.2f}, which is why "
            "no multi-night period in this paper is presented as a "
            "determination: each is a confirmation of the catalogue value "
            "within a stated alias family, and the family is named."),
        tables=("p3_pgram", "p3_period", "p3_ephemeris"),
        width_in=COL_DOUBLE)
    return fig, spec


def _fold_panel(ax, cv, series_key, eph_row, states, n_bins=40,
                show_points=True):
    """One folded-light-curve panel, split by accretion state.

    Shared by figures 5 and 7 so that ST LMi's fold and VV Pup's fold are
    the same statistic drawn the same way -- a reader comparing them is
    comparing stars, not plotting conventions.
    """
    pts = read_target_points(cv, series_key)
    if pts["bjd"].size == 0:
        _empty_panel(ax, "no catalogue-tied target points")
        return 0
    period = float(eph_row["period_d"])
    epoch = eph_row["epoch_bjd"]
    if epoch is None:
        _empty_panel(ax, "no published epoch: fold not defined")
        return 0
    ph = fold_phase(pts["bjd"], period, float(epoch))
    night_state = {r["night"]: normalise_state(r["state"]) for r in states}
    groups: dict[str, list[int]] = {}
    for i, n in enumerate(pts["night"]):
        groups.setdefault(night_state.get(n, "unknown"), []).append(i)
    # Draw the unclassified nights first and underneath: they are context,
    # not result, and a legend that leads with them buries the two states
    # the paper is about.
    order = sorted(groups, key=lambda s: (s in ("unknown", "unclassified"),
                                          s))
    for state in order:
        idx = groups[state]
        col = STATE_COLOR.get(state, OKABE_ITO["grey"])
        sel = np.array(idx, dtype=int)
        if show_points:
            ax.plot(np.concatenate([ph[sel], ph[sel] + 1.0]),
                    np.concatenate([pts["mag"][sel], pts["mag"][sel]]),
                    marker=".", lw=0, ms=1.2, color=col, alpha=0.30,
                    zorder=2)
        c, m, s, k = phase_bin(ph[sel], pts["mag"][sel], n_bins=n_bins)
        good = np.isfinite(m)
        if good.sum() >= 3:
            ax.errorbar(np.concatenate([c[good], c[good] + 1.0]),
                        np.concatenate([m[good], m[good]]),
                        yerr=np.concatenate([s[good], s[good]]),
                        fmt="o", ms=2.4, lw=0.7, capsize=0, color=col,
                        mec="none", zorder=4,
                        label=f"{state} ({len(idx)} pts)")
    # Limits set by the bulk, not by the one frame that caught a satellite:
    # a single 20th-magnitude outlier otherwise flattens the entire fold
    # into the top centimetre of the panel.
    lo, hi = robust_ylim(pts["mag"], k=5.0)
    ax.set_ylim(hi, lo)                       # inverted: bright at the top
    ax.set_xlim(0, 2)
    ax.grid(color="#f2f2f2")
    return pts["bjd"].size


def fig05_stlmi_folds(cv):
    """Figure 5 -- ST LMi folded light curves, one column per era."""
    eras = [(7, ["G", "R", "I"], "High Gain (2024 Jan--Mar)"),
            (76, ["g", "r", "i"], "Mode0 (2025 Jan--2026 Feb)")]
    fig, axes = plt.subplots(3, 2, figsize=(COL_DOUBLE, 5.4),
                             sharex=True, squeeze=False)
    eph = read_ephemeris(cv)["stlmi"]
    for c, (era, bands, title) in enumerate(eras):
        for r, band in enumerate(bands):
            ax = axes[r][c]
            sk = f"stlmi|e{era}|{band}"
            states = read_rows(cv, """SELECT night, state FROM p3_state_night
                                      WHERE series_key=?""", (sk,))
            n = _fold_panel(ax, cv, sk, eph, states)
            ax.text(0.015, 0.06, band, transform=ax.transAxes, fontsize=8,
                    color=BAND_COLOR.get(band, "k"), fontweight="bold")
            if n:
                ax.legend(loc="upper right", fontsize=5.8, ncol=1)
            if r == 0:
                ax.set_title(title, fontsize=7.5, loc="left")
            if c == 0:
                ax.set_ylabel("catalogue-tied\nmagnitude (mag)", fontsize=6.5)
    for ax in axes[-1]:
        ax.set_xlabel("orbital phase (cycles, repeated)")

    spec = FigureSpec(
        fig_id="fig05", label="fig:stlmifold",
        title="ST LMi folded light curves, per era",
        caption=(
            "ST LMi folded on the catalogue ephemeris, shown twice over in "
            "phase, one column per instrument era and one row per band. "
            "Small points are individual exposures; large points are "
            "median phase bins with a median-absolute-deviation error on "
            "the median, and empty bins are left empty rather than "
            "interpolated. Colour separates the accretion states classified "
            "in \\texttt{p3\\_state\\_night}, and that state palette is "
            "the one thing the two columns do share. The two columns are "
            "NOT combined and share no MAGNITUDE axis, which is this "
            "figure's ordinate: the 2024 G/R/I and 2025 "
            "g/r/i seasons do not overlap in time and were taken through "
            "different cameras and different bandpasses, so the comparison "
            "the paper makes is between morphologies, never between "
            "magnitudes."),
        tables=("cv_lightcurve", "cv_frames", "p2_cloud_frame",
                "p3_ephemeris", "p3_state_night"),
        width_in=COL_DOUBLE)
    return fig, spec


def fig06_colour_phase(cv, max_dt_s=600.0):
    """Figure 6 -- colour against orbital phase.  SUBSTITUTE: ST LMi only."""
    eph = read_ephemeris(cv)["stlmi"]
    period, epoch = float(eph["period_d"]), float(eph["epoch_bjd"])
    # The panel list lives in numbers_cv.COLOUR_PANEL_PAIRS, because §3.3
    # quotes the range of tie bars THESE panels carry and the two must be
    # computed over the same four pairs.
    _pp: dict[int, list] = {}
    for _e, _a, _b in _nx.COLOUR_PANEL_PAIRS:
        _pp.setdefault(_e, []).append((_a, _b))
    eras = [(e, _pp[e][0], _pp[e][1],
             f"{ERA_LABEL.get(e, e)} ({'2024' if e == 7 else '2025'})")
            for e in sorted(_pp)]
    fig, axes = plt.subplots(2, 2, figsize=(COL_DOUBLE, 4.2), sharex=True,
                             squeeze=False)

    # The systematic floor on any colour drawn here: the two series' own
    # catalogue-tie residuals added in quadrature.  It is a SHIFT of the
    # whole curve, not a scatter, so it is drawn as a bar rather than as an
    # error on each point -- confusing the two is how a 25 mmag calibration
    # gets read as a 25 mmag measurement.
    tie = {r["series_key"]: r for r in read_rows(
        cv, "SELECT * FROM cv_cattie WHERE is_primary=1")}

    for c, (era, pair_a, pair_b, title) in enumerate(eras):
        for r, (ba, bb) in enumerate((pair_a, pair_b)):
            ax = axes[r][c]
            ska, skb = f"stlmi|e{era}|{ba}", f"stlmi|e{era}|{bb}"
            pa = read_target_points(cv, ska)
            pb = read_target_points(cv, skb)
            t, colr, dt = pair_quasi_simultaneous(
                pa["bjd"], pa["mag"], pb["bjd"], pb["mag"], max_dt_s)
            if t.size < 10:
                _empty_panel(ax, f"{ba}-{bb}: {t.size} quasi-simultaneous "
                                 f"pairs within {max_dt_s:.0f} s")
                continue
            ph = fold_phase(t, period, epoch)
            ax.plot(np.concatenate([ph, ph + 1]),
                    np.concatenate([colr, colr]), marker=".", lw=0, ms=1.4,
                    color=OKABE_ITO["grey"], alpha=0.35)
            cb, mb, sb, kb = phase_bin(ph, colr, n_bins=30)
            good = np.isfinite(mb)
            ax.errorbar(np.concatenate([cb[good], cb[good] + 1]),
                        np.concatenate([mb[good], mb[good]]),
                        yerr=np.concatenate([sb[good], sb[good]]),
                        fmt="o", ms=2.6, lw=0.8, capsize=0,
                        color=BAND_COLOR.get(ba, "k"), mec="none")
            # TWO BARS, because the tie accuracy depends on whether outlier
            # check stars are clipped and the paper may not make that
            # choice invisibly.  Inner solid bar: the sigma-clipped
            # residuals.  Outer pale bar: every held-out star kept.  The
            # colour zero point is uncertain by somewhere between them, and
            # a single bar drawn at the smaller one understates it.
            ra = tie.get(ska, {}).get("check_rms_clip") or 0.0
            rb = tie.get(skb, {}).get("check_rms_clip") or 0.0
            ua = tie.get(ska, {}).get("check_rms") or 0.0
            ub = tie.get(skb, {}).get("check_rms") or 0.0
            sysbar = math.hypot(float(ra), float(rb))
            sysbar_raw = math.hypot(float(ua), float(ub))
            lo, hi = robust_ylim(colr)
            y0 = lo + 0.14 * (hi - lo)
            if sysbar_raw > sysbar:
                ax.errorbar([0.09], [y0], yerr=[sysbar_raw], fmt="none",
                            ecolor=OKABE_ITO["grey"], lw=2.6, alpha=0.45,
                            capsize=2.6)
            ax.errorbar([0.09], [y0], yerr=[sysbar], fmt="none",
                        ecolor=OKABE_ITO["black"], lw=1.1, capsize=2.0)
            ax.text(0.13, y0,
                    f"tie systematic\n{1000 * sysbar:.0f} mmag clipped\n"
                    f"{1000 * sysbar_raw:.0f} mmag unclipped",
                    fontsize=5.6, va="center")
            ax.set_ylim(lo, hi)
            ax.set_ylabel(f"${ba}-{bb}$ (mag)", fontsize=7)
            ax.text(0.985, 0.93, f"{t.size} pairs, $\\Delta t<${max_dt_s:.0f} s",
                    transform=ax.transAxes, ha="right", va="top", fontsize=5.8,
                    color=OKABE_ITO["grey"])
            ax.grid(color="#f2f2f2")
            if r == 0:
                ax.set_title(title, fontsize=7.5, loc="left")
    for ax in axes[-1]:
        ax.set_xlabel("orbital phase (cycles, repeated)")
        ax.set_xlim(0, 2)

    # The census that justifies the substitution, printed INSIDE the figure
    # and computed by the SAME function the manuscript's census macros use,
    # so the figure cannot state one coverage and §2.2 another.
    _mp = _nx.one(cv, "SELECT value FROM p4_meta WHERE "
                      "key='full_orbit_min_points'")
    _mp = int(_mp) if _mp is not None else _nx.FULL_ORBIT_MIN_POINTS_DEFAULT
    _cen = {t: _nx.coverage_census(cv, t, _mp) for t in ("vvpup", "euuma")}
    _vv, _eu = _cen["vvpup"], _cen["euuma"]
    fig.text(0.5, -0.02,
             f"VV Pup and EU UMa panels are not drawn: VV Pup has "
             f"{len(_vv['three'])} of {len(_vv['any'])} three-filter "
             f"full-orbit nights (and its two cameras are fully "
             f"confounded with epoch); EU UMa has {len(_eu['three'])} of "
             f"{len(_eu['any'])}.",
             ha="center", fontsize=6.2, color=OKABE_ITO["vermilion"])
    _p_min = 1440.0 * float(read_rows(
        cv, "SELECT period_d FROM p3_ephemeris WHERE target_key='stlmi'"
    )[0]["period_d"])

    spec = FigureSpec(
        fig_id="fig06", label="fig:colourphase",
        title="Colour against orbital phase, ST LMi, per era",
        caption=(
            "Quasi-simultaneous colours of ST LMi against orbital phase, "
            "one column per instrument era. A pair enters only if the two "
            f"exposures are within {max_dt_s:.0f} s of each other, because "
            f"these stars move by tenths of a magnitude within one "
            f"{_p_min:.0f}-minute orbit and "
            "a pair separated by half an orbit is not a colour. Small "
            "points are individual pairs, large points are median phase "
            "bins. The bars in each panel are the catalogue-tie "
            "systematic, the two series' check-star residuals added in "
            "quadrature, drawn twice: the black bar from the "
            "SIGMA-CLIPPED residuals and the pale bar behind it from the "
            "same residuals with every held-out star kept, because that "
            "choice moves the zero-point uncertainty by a factor of a few "
            "and Section~\\ref{sec:tie} declines to make it silently. It "
            "shifts the whole curve and does not scatter it, "
            "so colour SHAPE is measured here at the per-point precision "
            "while colour ZERO POINT carries that bar. No colour-dependent "
            "extinction term is applied anywhere in this paper; the 3$\\sigma$ "
            "bound on that term is tabulated instead."),
        tables=("cv_lightcurve", "cv_frames", "cv_cattie", "p3_ephemeris"),
        width_in=COL_DOUBLE, substitute=True,
        substitute_reason=(
            f"the planned figure paired ST LMi with per-camera VV Pup "
            f"panels; VV Pup has {len(_vv['three'])} three-filter "
            f"full-orbit nights of {len(_vv['any'])} and EU UMa "
            f"{len(_eu['three'])} of {len(_eu['any'])}, so there is no "
            f"night from which those panels could be built"))
    return fig, spec


def fig07_vvpup_euuma_folds(cv):
    """Figure 7 -- VV Pup (per camera) and EU UMa folded curves. SUBSTITUTE."""
    picks = [("vvpup|e72|g", "VV Pup — 1MHz HS (2024-11/12)"),
             ("vvpup|e76|g", "VV Pup — Mode0 (2025)"),
             ("vvpup|e76|r", "VV Pup — Mode0 (2025)"),
             ("euuma|e76|g", "EU UMa — Mode0 (2025)")]
    fig, axes = plt.subplots(2, 2, figsize=(COL_DOUBLE, 4.2), sharex=True,
                             squeeze=False)
    eph = read_ephemeris(cv)
    flat = [a for row in axes for a in row]
    for ax, (sk, title) in zip(flat, picks):
        tgt, era, band = series_parts(sk)
        states = read_rows(cv, """SELECT night, state FROM p3_state_night
                                  WHERE series_key=?""", (sk,))
        n = _fold_panel(ax, cv, sk, eph[tgt], states, n_bins=30)
        ax.set_title(f"{title}  {band}", fontsize=7, loc="left")
        ax.set_ylabel("catalogue-tied\nmagnitude (mag)", fontsize=6.5)
        if n:
            ax.legend(loc="lower right", fontsize=5.6)
    for ax in axes[-1]:
        ax.set_xlabel("orbital phase (cycles, repeated)")

    fig.text(0.5, -0.02,
             "EU UMa's merged 2026 Fast series is excluded from this and "
             "every figure: five comparison stars, ZERO check stars, and no "
             "catalogue tie — it cannot carry a validated magnitude.",
             ha="center", fontsize=6.2, color=OKABE_ITO["vermilion"])

    spec = FigureSpec(
        fig_id="fig07", label="fig:vvpupeuumafold",
        title="VV Pup and EU UMa folded light curves",
        caption=(
            "Folded light curves of VV Pup and EU UMa on their catalogue "
            "ephemerides, drawn exactly as Figure~\\ref{fig:stlmifold}. VV "
            "Pup is split by CAMERA and never combined across the two: its "
            "iKon (2024 Nov--Dec) and Mode0 (2025) epochs are fully "
            "confounded with instrument, so a joint fold would put a "
            "zero-point step inside a phase curve. EU UMa is shown only in "
            "the Mode0 $g$ series that has a catalogue tie."),
        tables=("cv_lightcurve", "cv_frames", "p2_cloud_frame",
                "p3_ephemeris", "p3_state_night"),
        width_in=COL_DOUBLE, substitute=True,
        substitute_reason=(
            "the planned figure was a three-filter panel per target; only "
            "the bands with full-orbit coverage are shown, and EU UMa's "
            "2026 Fast series carries no target photometry to fold at "
            "all -- it is untied, holds out no check stars, and the star "
            "is not detected in it"))
    return fig, spec


def fig08_state_history(cv, first_night="2018-01-01"):
    """Figure 8 -- long-term accretion-state histories over the survey record.

    ``first_night`` defaults to the start of ZTF.  AAVSO visual estimates
    for these stars reach back forty years; drawn on the same axis they
    compress the modern, calibrated record into a sliver and their own
    scatter dominates the panel.  The window is stated in the caption and
    the number of earlier epochs is printed in the figure, so the choice is
    visible rather than silent.
    """
    targets = ["stlmi", "vvpup", "euuma", "anuma", "yzcnc"]
    # 7.0 in, not 7.6.  At \textwidth in a two-column AASTeX figure* the
    # panel is reproduced at its full 7.1 in width, so its height goes onto
    # the page unscaled; at 7.6 in this float plus its caption overflowed
    # the text block by 31.9 pt and LaTeX warned on every build.  A figure
    # that does not fit its page is a layout defect, and the fix belongs
    # in the code that draws it rather than in a \resizebox in the source.
    fig, axes = plt.subplots(len(targets), 1, figsize=(COL_DOUBLE, 7.0),
                             sharex=True, squeeze=False)
    src_color = {"ztf": "#9ab4d0", "asassn": "#d0b49a", "aavso": "#b7c9a8"}
    x_min = float(night_to_ordinal([first_night])[0])
    n_before = 0

    for ax, tgt in zip([a for row in axes for a in row], targets):
        ext = read_rows(cv, """SELECT source, utc_night, mag FROM
                               cv_ext_nightly WHERE target=? AND mag IS NOT
                               NULL AND independent=1""", (tgt,))
        yvals: list[float] = []
        for src in ("aavso", "asassn", "ztf"):
            sel = [r for r in ext if r["source"] == src]
            if not sel:
                continue
            x = night_to_ordinal([r["utc_night"] for r in sel])
            y = np.array([r["mag"] for r in sel], dtype=float)
            keep = np.isfinite(x) & (x >= x_min)
            n_before += int((np.isfinite(x) & (x < x_min)).sum())
            ax.plot(x[keep], y[keep], ".", ms=1.5, color=src_color[src],
                    alpha=0.75, zorder=1, label=f"{src.upper()} nightly")
            yvals.extend(y[keep].tolist())

        nights = read_rows(cv, """
            SELECT n.series_key, n.night, n.median_mag, n.state, n.censored,
                   s.filter
            FROM p3_state_night n JOIN cv_series s
                 ON s.series_key = n.series_key
            WHERE n.target_key = ? AND n.median_mag IS NOT NULL""", (tgt,))
        for r in nights:
            x = night_to_ordinal([r["night"]])[0]
            if not math.isfinite(x) or x < x_min:
                continue
            col = STATE_COLOR[normalise_state(r["state"])]
            mk = BAND_MARKER.get(r["filter"], "o")
            yvals.append(float(r["median_mag"]))
            if r["censored"]:
                ax.plot([x], [r["median_mag"]], marker="v", ms=4.0, lw=0,
                        mfc="none", mec=col, mew=0.9, zorder=5)
            else:
                ax.plot([x], [r["median_mag"]], marker=mk, ms=3.6, lw=0,
                        mfc=col, mec="k", mew=0.35, zorder=5)
        thr = read_rows(cv, """SELECT threshold_mag, separability, verdict
                               FROM p3_state_series WHERE target_key=?
                               AND threshold_mag IS NOT NULL""", (tgt,))
        for t in thr:
            ax.axhline(float(t["threshold_mag"]), lw=0.6, ls="--",
                       color=OKABE_ITO["black"], alpha=0.5)
        lim = read_rows(cv, """SELECT median_limit_cal_mag, n_limits
                               FROM p2_limit_series WHERE target_key=?
                               AND median_limit_cal_mag IS NOT NULL""",
                        (tgt,))
        if lim:
            med = float(np.median([r["median_limit_cal_mag"] for r in lim]))
            ax.axhline(med, lw=0.7, ls=":", color=OKABE_ITO["vermilion"])
            ax.text(0.995, 0.06, f"median faint limit {med:.1f} mag "
                                 f"— shallower than the median detection, "
                                 f"so faint fractions are UPPER BOUNDS",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=5.2, color=OKABE_ITO["vermilion"])
        lo, hi = robust_ylim(yvals, k=4.0)
        ax.set_ylim(hi, lo)                   # inverted: bright at the top
        ax.set_ylabel(f"{TARGET_LABEL[tgt]}\nmagnitude", fontsize=6.5)
        ax.grid(color="#f4f4f4")
        ax.tick_params(labelsize=6)
    ax_last = axes[-1][0]
    lo, hi = ax_last.get_xlim()
    pos, lab = year_ticks(max(lo, x_min), hi, max_ticks=10)
    ax_last.set_xlim(max(lo, x_min), hi)
    ax_last.set_xticks(pos)
    ax_last.set_xticklabels(lab, fontsize=6)
    ax_last.set_xlabel("UTC night")
    axes[0][0].text(0.5, 1.62,
                    f"axis begins {first_night} (start of ZTF); "
                    f"{n_before:,} earlier survey epochs, mostly visual "
                    f"AAVSO estimates, lie off the left of every panel",
                    transform=axes[0][0].transAxes, ha="center", va="bottom",
                    fontsize=5.6, color=OKABE_ITO["grey"])
    # EVERY class the panels actually draw.  p3_state_night holds five, and
    # a legend that documents two of them leaves the grey and yellow
    # markers that dominate the ST LMi and EU UMa panels unexplained --
    # while Figure 5's legend, drawn from the same palette, lists them.
    counts = {str(r["state"]).lower(): r["n"] for r in read_rows(
        cv, "SELECT state, count(*) n FROM p3_state_night GROUP BY state")}
    # The caption's shallow-limit tally, by the same test §3.5 applies: a
    # limit is informative only if it goes deeper than the star's own
    # median detection.  Typed into the caption, this is a number that
    # cannot follow a re-run of the forced-photometry stage.
    lim_rows = read_rows(cv, """SELECT series_key, median_limit_cal_mag
                                FROM p2_limit_series
                                WHERE median_limit_cal_mag IS NOT NULL""")
    n_limit, n_shallow = len(lim_rows), 0
    for r in lim_rows:
        mags = sorted(x["cal_mag"] for x in read_rows(
            cv, "SELECT cal_mag FROM cv_lightcurve WHERE series_key=? AND "
                "role='target' AND cal_mag IS NOT NULL", (r["series_key"],)))
        if mags and r["median_limit_cal_mag"] < mags[len(mags) // 2]:
            n_shallow += 1
    handles = [
        Line2D([], [], marker="o", lw=0, ms=4,
               mfc=STATE_COLOR[normalise_state(s)], mec="none",
               label=f"RLMT night, {lab} ({counts.get(s, 0)})")
        for s, lab in (("high", "high state"), ("low", "low state"),
                       ("intermediate", "inside the threshold's own "
                                       "uncertainty"),
                       ("unclassified", "coverage gate not cleared"),
                       ("unknown", "series has no measurable threshold"))
        if counts.get(s)]
    handles += [
        Line2D([], [], marker="v", lw=0, ms=4, mfc="none",
               mec=OKABE_ITO["grey"], label="censored night (upper limit)"),
        Line2D([], [], ls="--", color="k", lw=0.8, label="state threshold"),
        Line2D([], [], ls=":", color=OKABE_ITO["vermilion"], lw=0.9,
               label="median faint limit"),
    ]
    # Above the top panel, not inside it: eight entries laid over ST LMi's
    # nights hid the very markers they were there to explain.
    axes[0][0].legend(handles=handles, loc="lower center",
                      bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=5.4,
                      columnspacing=1.0, handletextpad=0.4)

    spec = FigureSpec(
        fig_id="fig08", label="fig:states",
        title="Long-term accretion-state histories",
        caption=(
            "Nightly RLMT medians (large symbols, coloured by classified "
            "accretion state, marker by filter) over the independent survey "
            "record (small pale points: ZTF, ASAS-SN, AAVSO nightly means, "
            "with any RLMT data resubmitted to AAVSO removed). The legend "
            "lists every class the panels draw, with the number of nights "
            "in each: two of the five are states, and the other three are "
            "ways of saying the classification could not be made. Dashed "
            "lines "
            "are the per-series state thresholds; open triangles are nights "
            "on which the star was not detected and only an upper limit "
            f"exists. The dotted red line is the median faint limit: for "
            f"{n_shallow} of {n_limit} series it lies SHALLOWER than the "
            "median "
            "detection, so every faint-state fraction in this paper is an "
            "upper bound on a low-state duty cycle and not a measurement "
            "of one."),
        tables=("p3_state_night", "p3_state_series", "cv_ext_nightly",
                "p2_limit_series", "cv_series"),
        width_in=COL_DOUBLE)
    return fig, spec


def fig09_oc(cv):
    """Figure 9 -- ST LMi O--C and band-dependent edge offsets.  SUBSTITUTE."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_DOUBLE, 3.4),
                             gridspec_kw={"width_ratios": [1.75, 1.0],
                                          "wspace": 0.40})
    ax, ax2 = axes

    # THE PUBLISHED EPOCHS: one per night per band, from p3_oc_night.  The
    # per-cycle residuals in p3_oc are the inputs -- CV-S5's injection test
    # does not license a single cycle's edge as an epoch, and §4.2 forbids
    # publishing one -- so they are drawn as pale background points, and
    # the epochs the paper fits are the filled ones with error bars.
    raw = read_rows(cv, """SELECT * FROM p3_oc WHERE target_key='stlmi'
                           ORDER BY cycle""")
    for r in raw:
        band = series_parts(r["series_key"])[2]
        ax.plot([r["cycle"]], [r["oc_s"]], marker=".", ms=2.0, lw=0,
                color=BAND_COLOR.get(band, OKABE_ITO["grey"]), alpha=0.30,
                zorder=2)

    oc = read_rows(cv, """SELECT * FROM p3_oc_night WHERE target_key='stlmi'
                          ORDER BY cycle_mean""")
    by_series: dict[str, list[dict]] = {}
    for r in oc:
        by_series.setdefault(r["series_key"], []).append(r)
    for sk, rows in sorted(by_series.items()):
        _, era, band = series_parts(sk)
        e = np.array([r["cycle_mean"] for r in rows], dtype=float)
        y = np.array([r["oc_s"] for r in rows], dtype=float)
        s = np.array([(r["oc_sigma_s"] or np.nan) for r in rows], dtype=float)
        ax.errorbar(e, y, yerr=s, fmt=BAND_MARKER.get(band, "o"), ms=3.4,
                    lw=0, elinewidth=0.7, capsize=0,
                    color=BAND_COLOR.get(band, OKABE_ITO["grey"]),
                    mfc=("none" if era == 7 else
                         BAND_COLOR.get(band, OKABE_ITO["grey"])),
                    mec=BAND_COLOR.get(band, OKABE_ITO["grey"]), mew=0.8,
                    zorder=4,
                    label=f"{ERA_LABEL.get(era, era)} {band}")
    ax.axhline(0.0, color="k", lw=0.7)

    # The two bands are the two ERROR SCALES this diagram lives between:
    # what a single cycle's edge achieved in the injection test (outer),
    # and what a per-night mean of several cycles achieves (inner).  The
    # published epochs are the second; the first is drawn to show why they
    # are not the first.
    cc = read_rows(cv, "SELECT * FROM p3_cycle_count WHERE "
                       "target_key='stlmi'")
    sig = read_rows(cv, """SELECT sigma_t_s, total_error_s FROM p3_sigmat""")
    thr = 60.0
    if sig and cc:
        s_cycle = float(np.nanmedian([r["total_error_s"] for r in sig]))
        s_night = float(cc[0]["sigma_night_median_s"])
        ax.axhspan(-s_cycle, s_cycle, color=OKABE_ITO["grey"], alpha=0.14,
                   lw=0, zorder=0)
        ax.axhspan(-s_night, s_night, color=OKABE_ITO["grey"], alpha=0.22,
                   lw=0, zorder=0)
        for sgn in (-1, 1):
            ax.axhline(sgn * thr, lw=0.7, ls="-.",
                       color=OKABE_ITO["vermilion"])
        # Reduced chi-squared on nu = N - 1: one constant, the edge's own
        # phase offset from the catalogue epoch, was fitted out of these
        # same edges (see the caption), so the denominator is not N.
        _nu = max(len(oc) - 1, 1)
        _chi2 = sum((r["oc_s"] / r["oc_sigma_s"]) ** 2 for r in oc
                    if r["oc_sigma_s"])
        ax.text(0.012, 0.025,
                f"per-CYCLE {s_cycle:.0f} s (outer band): fails {thr:.0f} s\n"
                f"per-NIGHT {s_night:.0f} s (inner band); rms "
                f"{float(cc[0]['oc_night_rms_s']):.0f} s, "
                f"$\\chi^2/\\nu$ = {_chi2 / _nu:.2f} ($\\nu$ = {_nu})",
                transform=ax.transAxes, fontsize=5.6, va="bottom",
                color=OKABE_ITO["black"])
        # THE BOUND THE NULL BUYS.  A null is worth what it excludes, so
        # draw the O-C a period derivative at the 3-sigma limit would have
        # left BEHIND in these residuals: the quadratic a2*E^2 with the
        # constant and linear terms removed the way the fit removes them --
        # by weighted least squares under the fit's own 1/sigma^2 weights,
        # since those two terms are degenerate with the epoch and the
        # period and a real Pdot would be partly absorbed into both.
        #
        # It was drawn as a bare re-centred parabola, which is NOT what the
        # fit absorbs, and the caption then called it a region the epochs
        # "sit inside".  They do not and cannot: the envelope is a SIGNAL
        # SHAPE, and the epochs scatter about zero by their own 84 s error,
        # which is larger than the envelope over most of the baseline --
        # that being precisely why no such curvature is detectable.  The
        # caption is now generated from the drawn curve (see
        # ``_envelope_report``) and refuses to assert containment.
        lim = cc[0]["pdot_limit3"] if "pdot_limit3" in cc[0].keys() else None
        per_d = cc[0]["period_d"]
        env_report = None
        if lim and per_d and oc:
            e_ep = np.array([r["cycle_mean"] for r in oc], dtype=float)
            s_ep = np.array([r["oc_sigma_s"] for r in oc], dtype=float)
            y_ep = np.array([r["oc_s"] for r in oc], dtype=float)
            e_ax = np.linspace(e_ep.min(), e_ep.max(), 400)
            env_ax = pdot_envelope_seconds(e_ax, e_ep, s_ep, float(lim),
                                           float(per_d))
            env_ep = pdot_envelope_seconds(e_ep, e_ep, s_ep, float(lim),
                                           float(per_d))
            for sgn in (-1.0, 1.0):
                ax.plot(e_ax, sgn * np.abs(env_ax), lw=0.8, ls=":",
                        color=OKABE_ITO["green"], zorder=1)
            ax.plot([], [], lw=0.8, ls=":", color=OKABE_ITO["green"],
                    label=f"$|\\dot{{P}}|$ = {float(lim):.1e} (3$\\sigma$ "
                          f"limit)")
            env_report = _envelope_report(y_ep, s_ep, env_ep)
    ax.set_xlabel("cycle number since the catalogue epoch")
    ax.set_ylabel("O$-$C (s)")
    ax.grid(color="#f2f2f2")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=4,
              fontsize=5.6, columnspacing=0.9, handletextpad=0.35)

    # Right: the band-to-band edge offsets, which the strategy demanded not
    # be averaged away.  The result is a NULL, and the panel has to show it
    # as one: every pooled offset is printed with its error bar and its
    # significance, so a reader can see that the largest reaches 1.9 sigma
    # against a 3 sigma bar rather than having to infer a detection from
    # the mere existence of non-zero differences.
    #
    # Only pairs WITH a pooled measurement get a row.  A pair seen on one
    # night has no pooled estimate (the stage needs two paired cycles), and
    # an earlier version of this panel drew that row anyway, with a lone
    # pale marker and no diamond and no explanation.
    bp = read_rows(cv, """SELECT * FROM p3_band_pair WHERE target_key='stlmi'
                          ORDER BY era_id, band_a, band_b, night""")
    groups: dict[tuple, list[dict]] = {}
    for r in bp:
        groups.setdefault((r["era_id"], r["band_a"], r["band_b"]),
                          []).append(r)

    def _pooled_of(rows_):
        return [r for r in rows_ if "pooled" in str(r["night"]).lower()]

    keys = sorted(k for k in groups if _pooled_of(groups[k]))
    dropped = sorted(set(groups) - set(keys))
    n_sig = sum(1 for r in bp if r["significant"])
    for i, key in enumerate(keys):
        era, ba, bb = key
        rows_ = groups[key]
        singles = [r for r in rows_ if "pooled" not in str(r["night"]).lower()]
        pooled = _pooled_of(rows_)
        col = BAND_COLOR.get(ba, OKABE_ITO["grey"])
        if singles:
            jit = np.linspace(-0.22, 0.22, len(singles))
            ax2.errorbar([r["delta_s"] for r in singles], i + jit,
                         xerr=[r["sigma_s"] for r in singles], fmt="o",
                         ms=2.0, lw=0, elinewidth=0.5, capsize=0,
                         color=col, ecolor="#d5d5d5", alpha=0.85, mec="none")
        p0 = pooled[0]
        ax2.errorbar([p0["delta_s"]], [i], xerr=[p0["sigma_s"]], fmt="D",
                     ms=3.4, lw=0, elinewidth=1.1, capsize=1.8,
                     color=OKABE_ITO["black"], zorder=5)
        # THE NUMBER, ON THE PANEL.  Panel (a) already annotates its own
        # error scales; the deciding number of panel (b) is a significance
        # and it belongs in print, not in a database a reader has to open.
        ax2.annotate(
            f"{p0['delta_s']:+.0f}$\\pm${p0['sigma_s']:.0f} s "
            f"({abs(p0['delta_s']) / p0['sigma_s']:.1f}$\\sigma$, "
            f"n={int(p0['n_cycles'])})",
            xy=(p0["delta_s"], i), xytext=(0, 9.0),
            textcoords="offset points", ha="center", fontsize=4.9,
            color=OKABE_ITO["black"], zorder=6)
    if keys:
        ax2.axvline(0.0, color="k", lw=0.7)
        ax2.set_yticks(np.arange(len(keys)))
        ax2.set_yticklabels(
            [f"{ba}$-${bb}  {ERA_LABEL.get(era, era)[:9]}"
             for era, ba, bb in keys], fontsize=5.8)
        ax2.set_ylim(-0.6, len(keys) - 0.15)
        ax2.set_xlabel("band-to-band edge offset (s)")
        ax2.grid(axis="x", color="#f2f2f2")
        ax2.set_title(f"{n_sig} of {len(bp)} significant at 3$\\sigma$",
                      fontsize=6.0, color=OKABE_ITO["vermilion"], pad=3)
        ax2.legend(handles=[
            Line2D([], [], marker="o", lw=0, ms=3,
                   color=OKABE_ITO["grey"], label="one night"),
            Line2D([], [], marker="D", lw=0, ms=3.4,
                   color=OKABE_ITO["black"], label="pooled")],
            loc="lower left", fontsize=5.4)
    else:
        _empty_panel(ax2, "no paired band edges")

    fig.text(0.5, -0.03,
             "AN UMa and VV Pup have NO O$-$C: CV-S9 graded their accepted "
             "edges 'not one feature' (phase scatter 0.15 and 0.13 against "
             "a 0.05 bar). EU UMa has two accepted edges.",
             ha="center", fontsize=6.0, color=OKABE_ITO["vermilion"])

    # The caption's numbers come from the same rows the panel drew, so a
    # re-run that moved a measurement moves the caption with it.
    _cc0 = cc[0] if cc else {}
    _n_single = sum(1 for r in oc if r["n_cycles"] == 1)
    _n_xfer = sum(1 for r in oc if r["era_id"] != 76)
    _pool = [r for r in bp if "pooled" in str(r["night"]).lower()
             and r["sigma_s"]]
    _top = (max(_pool, key=lambda r: abs(r["delta_s"]) / r["sigma_s"])
            if _pool else None)
    # BOTH ENDS of the pooled bound.  The caption used to quote the
    # tightest beside the words "any such offset", which is a universal
    # quantifier carrying the most favourable of five numbers; the panel
    # itself prints all five honestly, so the caption disagreed with the
    # figure it describes.
    def _b2(r):
        return abs(r["delta_s"]) + 2 * r["sigma_s"]
    _bnd = min(_pool, key=_b2) if _pool else None
    _weak = max(_pool, key=_b2) if _pool else None
    # WHAT THESE RESIDUALS ARE MEASURED AGAINST, SAID IN THE CAPTION.
    # CV-S9 subtracts the mean per-cycle O-C before writing p3_oc: the
    # bright-phase edge is not at phase zero of the VSX ephemeris, so the
    # raw residual carries a constant that belongs to the feature.  The
    # subtraction is right and it was invisible -- the caption opened
    # "residuals against the catalogue ephemeris", which a reader takes to
    # mean the edge falls at the catalogue epoch's phase.  It does not.
    _off_s = _cc0["oc_mean_s"] if "oc_mean_s" in _cc0.keys() else None
    _per_s = (float(_cc0["period_d"]) * 86400.0
              if "period_d" in _cc0.keys() and _cc0["period_d"] else None)
    _offset_txt = ""
    if _off_s and _per_s:
        _offset_txt = (
            f"The timed bright-phase edge does not fall at the catalogue "
            f"epoch's phase zero: it sits {float(_off_s):,.0f}~s "
            f"({float(_off_s) / _per_s:.3f} of a cycle) after it, that "
            f"constant is a property of the FEATURE and of the "
            f"catalogue's choice of fiducial rather than of the clock, and "
            f"it is removed before anything here is plotted or fitted. "
            f"What this panel therefore tests is whether that interval is "
            f"CONSTANT, never whether it is zero.")

    # WHAT THE DOTTED ENVELOPE IS, FROM THE CURVE THE CODE DREW.  The
    # previous caption ended "and the epochs sit inside it"; 23 of the 36
    # lie outside, because the envelope is a signal shape and not an error
    # band.  ``_envelope_report`` measures the relation and this clause
    # states it; there is no branch here that can assert containment.
    _env_txt = ""
    if env_report:
        _env_txt = (
            "The dotted envelope is not an error band and the epochs are "
            "not expected to lie within it: it is the O$-$C a steady "
            "period derivative at this data set's 3$\\sigma$ upper bound "
            "would still have left behind after the constant and linear "
            "terms are absorbed as the ephemeris fit absorbs them, i.e. "
            "the CURVATURE the epochs would have had to follow. It reaches "
            f"only {env_report['env_max']:.0f}~s at the ends of the "
            f"baseline and is smaller than a single epoch's error bar on "
            f"{env_report['n_under_error']} of the {env_report['n']} "
            f"epochs, against a median epoch error of "
            f"{env_report['sigma_median']:.0f}~s; "
            f"{env_report['n_outside']} of the {env_report['n']} epochs "
            "scatter further from zero than the envelope does. That the "
            "scatter is larger than the signal, and shows no curvature, is "
            "exactly why no such $\\dot{P}$ can be distinguished.")

    _drop_txt = ""
    if dropped:
        _drop_txt = (
            " Band pairs seen on a single night carry no pooled estimate "
            "and are not drawn ("
            + ", ".join(f"${ba}-{bb}$ in {ERA_LABEL.get(era, era)}"
                        for era, ba, bb in dropped) + ").")
    _null_txt = ""
    if _top is not None and _bnd is not None:
        _null_txt = (
            f" Every pooled offset is consistent with zero: {n_sig} of "
            f"{len(bp)} pairs are significant at 3$\\sigma$, the strongest "
            f"pooled pair is ${_top['band_a']}-{_top['band_b']}$ at "
            f"{_top['delta_s']:+.0f}$\\pm${_top['sigma_s']:.0f}~s "
            f"({abs(_top['delta_s']) / _top['sigma_s']:.1f}$\\sigma$ over "
            f"{int(_top['n_cycles'])} paired cycles). What the pairs bound "
            f"at 2$\\sigma$ spans a factor of "
            f"{_b2(_weak) / _b2(_bnd):.1f}: the tightest, "
            f"${_bnd['band_a']}-{_bnd['band_b']}$ in "
            f"{ERA_LABEL.get(_bnd['era_id'], _bnd['era_id'])}, excludes "
            f"offsets above {_b2(_bnd):.0f}~s, while the weakest, "
            f"${_weak['band_a']}-{_weak['band_b']}$ in "
            f"{ERA_LABEL.get(_weak['era_id'], _weak['era_id'])}, excludes "
            f"only those above {_b2(_weak):.0f}~s --- so an offset present "
            f"in EVERY pair is bounded at {_b2(_weak):.0f}~s and not at "
            f"{_b2(_bnd):.0f}~s. This panel is a NON-DETECTION and is "
            "reported as one.")
    spec = FigureSpec(
        fig_id="fig09", label="fig:oc",
        title="ST LMi O$-$C and the inter-band edge-offset null",
        caption=(
            "(a) Bright-phase timing residuals against the catalogue "
            "PERIOD. " + _offset_txt +
            " Symbols with error bars are the PUBLISHED epochs: "
            "one per night per band, each the mean of that night's "
            "accepted per-cycle edges, open for the 2024 High Gain era and "
            "filled for the 2025 Mode0 era. Pale dots behind them are the "
            "per-cycle edges those means are made of; none is published "
            "carrying its own per-cycle error bar, because the injection "
            "test of Section~\\ref{sec:timing} showed that one cycle's "
            "edge does not reach the 60~s threshold (dash-dot lines). On "
            f"{_n_single} of the {len(oc)} epochs only one cycle was timed, "
            "so the mean is that cycle with the injection budget attached "
            "instead of the fit's own error. That budget is measured on a "
            "SINGLE Mode0 night of ST~LMi and served out by band slot, so "
            f"the {_n_xfer} epochs from the other instrument eras carry an "
            "error bar measured in an era they were not observed in; "
            "Section~\\ref{sec:timing} gives the reduced chi-squared "
            "recomputed under the edge fits' own errors as the check on "
            "that. The outer shaded band is the per-cycle error, the inner "
            "the median error of a per-night epoch. " + _env_txt + " (b) "
            "Band-to-band "
            "offsets of the same edge on the same cycle, shown rather than "
            "averaged away, one row per band pair with the pooled estimate "
            "as the black diamond and its value, error and significance "
            "printed above it." + _null_txt + _drop_txt
            + " Bands are nonetheless averaged separately throughout, as "
            "the conservative choice against a wavelength-dependent edge "
            "phase that the cyclotron picture predicts and these data are "
            "not precise enough to detect or exclude at the size it would "
            "have."),
        tables=("p3_oc", "p3_oc_night", "p3_band_pair", "p3_sigmat",
                "p3_cycle_count"),
        width_in=COL_DOUBLE, substitute=True,
        substitute_reason=(
            "the planned figure carried VV Pup and EU UMa panels; CV-S9's "
            "cycle-count stage graded VV Pup and AN UMa 'NOT ONE FEATURE "
            "-- NO O-C' and EU UMa has two accepted edges, so only ST LMi "
            "has a timeable feature to plot"))
    return fig, spec


def fig10_yzcnc_season(cv):
    """Figure 10 -- YZ Cnc's 2024 season with independent outburst states."""
    fig = plt.figure(figsize=(COL_DOUBLE, 4.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.5, 1.0], hspace=0.42,
                          wspace=0.28)
    ax = fig.add_subplot(gs[0, :])

    aavso = read_rows(cv, """SELECT utc_night, mag, state FROM cv_ext_nightly
                             WHERE target='yzcnc' AND source='aavso'
                               AND mag IS NOT NULL AND independent=1""")
    xa = night_to_ordinal([r["utc_night"] for r in aavso])
    ya = np.array([r["mag"] for r in aavso], dtype=float)
    order = np.argsort(xa)
    ax.plot(xa[order], ya[order], "-", lw=0.5, color="#cccccc", zorder=1)
    ax.plot(xa, ya, ".", ms=2.0, color=OKABE_ITO["grey"], alpha=0.8,
            zorder=2, label="AAVSO nightly (independent)")

    for era, bands in ((7, "GRI"), (72, "gri")):
        for band in bands:
            sk = f"yzcnc|e{era}|{band}"
            pts = read_target_points(cv, sk)
            if pts["bjd"].size == 0:
                continue
            x = night_to_ordinal(pts["night"])
            ax.plot(x, pts["mag"], marker=BAND_MARKER.get(band, "o"), lw=0,
                    ms=1.8, mfc=BAND_COLOR.get(band, "k"), mec="none",
                    alpha=0.75, zorder=3, label=f"RLMT {band} (e{era})")

    ob = read_rows(cv, """SELECT utc_night, episode, structure,
                                 amp_above_quiescence FROM p4_outburst""")
    shaded = set()
    for r in ob:
        x = night_to_ordinal([r["utc_night"]])[0]
        if math.isfinite(x) and r["utc_night"] not in shaded:
            ax.axvspan(x - 0.5, x + 0.5, color=OKABE_ITO["orange"],
                       alpha=0.18, lw=0, zorder=0)
            shaded.add(r["utc_night"])
    # The panel is the SEASON, not YZ Cnc's sixty-year AAVSO record.  With
    # the full baseline on the axis the 2024 campaign is a single vertical
    # smear and the outburst structure the figure exists to show is
    # invisible.  Limits come from the RLMT points themselves.
    rl_days = night_to_ordinal([r["night"] for r in read_rows(
        cv, "SELECT DISTINCT night FROM cv_frames WHERE target_key='yzcnc'")])
    rl_days = rl_days[np.isfinite(rl_days)]
    if rl_days.size:
        pad = 0.10 * max(20.0, rl_days.max() - rl_days.min())
        ax.set_xlim(rl_days.min() - pad, rl_days.max() + pad)
    ax.invert_yaxis()
    ax.set_ylabel("magnitude (mag)")
    ax.set_xlabel("UTC night")
    ax.grid(color="#f4f4f4")
    lo, hi = ax.get_xlim()
    ylo, yhi = robust_ylim(
        [y for x, y in zip(xa, ya) if lo <= x <= hi] + [10.0], k=5.0)
    ax.set_ylim(yhi, ylo)
    pos, lab = year_ticks(lo, hi, max_ticks=8)
    ax.set_xticks(pos)
    ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=6)
    ax.set_xlabel("")                    # the tick labels already say it
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.24), ncol=4,
              fontsize=5.6, columnspacing=1.0, handletextpad=0.4)
    ax.text(0.5, 1.02, "shaded: nights CV-S7 classified as normal outburst "
                       "from independent AAVSO photometry",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=5.6,
            color=OKABE_ITO["orange"])

    # Three dense-run insets, chosen by point count: what a survey epoch
    # cannot show.
    runs = read_rows(cv, """SELECT scope, series_key, nights, utc_nights,
                                   state, n_points, span_h
                            FROM p4_run WHERE n_points >= 40
                            ORDER BY n_points DESC""")
    picked, seen_night = [], set()
    for r in runs:
        key = str(r["utc_nights"])
        if key in seen_night:
            continue
        picked.append(r)
        seen_night.add(key)
        if len(picked) == 3:
            break
    for j, r in enumerate(picked):
        axi = fig.add_subplot(gs[1, j])
        night = str(r["nights"]).split("+")[0]
        pts = read_target_points(cv, r["series_key"], night=night)
        if pts["bjd"].size == 0:
            _empty_panel(axi, "no points")
            continue
        _, _, band = series_parts(r["series_key"])
        t0 = float(np.min(pts["bjd"]))
        axi.plot((pts["bjd"] - t0) * 24.0, pts["mag"], marker=".", lw=0,
                 ms=2.0, color=BAND_COLOR.get(band, "k"))
        axi.invert_yaxis()
        axi.set_xlabel("hours from run start", fontsize=6)
        axi.set_ylabel("mag", fontsize=6)
        axi.set_title(f"{str(r['utc_nights']).split('+')[0]} {band} "
                      f"({r['state']})", fontsize=6, loc="left")
        axi.tick_params(labelsize=5.5)
        axi.grid(color="#f4f4f4")

    # The run census and the peak amplitude, from the tables that measured
    # them, so the caption cannot outlive a re-run of CV-S10.
    n_dense = read_rows(cv, "SELECT count(DISTINCT nights) n FROM p4_run "
                            "WHERE kind='run'")[0]["n"]
    n_quiet = read_rows(cv, "SELECT count(DISTINCT nights) n FROM p4_run "
                            "WHERE kind='run' AND upper(state)='QUIESCENT'"
                        )[0]["n"]
    n_burst = read_rows(cv, "SELECT count(DISTINCT night) n FROM p4_outburst"
                        )[0]["n"]
    ob_amp = read_rows(cv, "SELECT max(amp_above_quiescence) a FROM "
                           "p4_outburst")[0]["a"] or float("nan")

    spec = FigureSpec(
        fig_id="fig10", label="fig:yzcncseason",
        title="YZ Cnc season overview with dense-run insets",
        caption=(
            "(top) YZ Cnc through the 2024 season: RLMT points in three "
            "bands over the independent AAVSO nightly record, with nights "
            "classified as normal outburst shaded. That classification is "
            "made from AAVSO photometry alone and survives deleting every "
            "RLMT row from the AAVSO archive, so it is an INPUT to this "
            "paper's YZ Cnc branch and not a conclusion of it. (bottom) "
            f"Three dense runs at full cadence: {n_dense} such runs exist, "
            f"{n_quiet} quiescent and {n_burst} inside normal outbursts, "
            f"and the brightest reaches {ob_amp:.2f} mag above quiescence "
            "against the roughly 3 mag a superoutburst attains."),
        tables=("cv_lightcurve", "cv_frames", "cv_ext_nightly",
                "p4_outburst", "p4_run"),
        width_in=COL_DOUBLE)
    return fig, spec


def _scope_label(run) -> str:
    """``yzcnc|e7|I|2024-02-20`` -> ``the 2024-02-21 $I$ run``.

    The night comes from :func:`final_science.run_night_label`, which is
    also what ``run_cv_final``'s verdict strings use, so Table 4 and this
    caption cannot name the same run by two different dates.
    """
    if run is None:
        return "none"
    band = series_parts(run["series_key"])[2]
    return (f"the {_fs.run_night_label(run['utc_nights'], run['nights'])} "
            f"${band}$ run")


def fig11_yzcnc_fallback(cv):
    """Figure 11 -- the fallback the strategy names: hump limits + flickering."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_DOUBLE, 3.1))
    ax, ax2 = axes

    # QUIESCENT only.  The hump is an ORBITAL modulation of the accretion
    # stream, and folding an outburst run on the orbital period measures
    # the outburst's slope wearing the orbit's clothes.  CV-S10 tested six
    # quiescent scopes and those are the six that belong on this panel.
    runs = read_rows(cv, """SELECT * FROM p4_run
                            WHERE hump_amp IS NOT NULL
                              AND amp90_self IS NOT NULL
                              AND upper(state) = 'QUIESCENT'
                            ORDER BY scope""")
    y = np.arange(len(runs))
    for i, r in enumerate(runs):
        _, _, band = series_parts(r["series_key"])
        col = BAND_COLOR.get(band, OKABE_ITO["grey"])
        amp = 1000.0 * float(r["hump_amp"])
        f90 = 1000.0 * float(r["amp90_field"] or np.nan)
        s90 = 1000.0 * float(r["amp90_self"] or np.nan)
        ax.plot([f90, s90], [i, i], "-", lw=0.8, color="#dddddd", zorder=1)
        ax.plot([amp], [i], "o", ms=3.4, color=col, mec="none", zorder=4)
        ax.plot([f90], [i], "|", ms=7, color=OKABE_ITO["green"], mew=1.2,
                zorder=3)
        ax.plot([s90], [i], "|", ms=7, color=OKABE_ITO["vermilion"], mew=1.2,
                zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{_fs.run_night_label(r['utc_nights'], r['nights'])} "
                        f"{series_parts(r['series_key'])[2]}"
                        for r in runs], fontsize=5.6)
    ax.set_xscale("log")
    ax.set_xlabel("semi-amplitude (mmag)")
    ax.set_ylim(-0.8, max(0.6, len(runs) - 0.4))
    ax.grid(axis="x", color="#f2f2f2")
    ax.legend(handles=[
        Line2D([], [], marker="o", lw=0, ms=4, color="k",
               label="fitted hump semi-amplitude"),
        Line2D([], [], marker="|", lw=0, ms=6, color=OKABE_ITO["green"],
               mew=1.2, label="90% recovery vs field stars"),
        Line2D([], [], marker="|", lw=0, ms=6, color=OKABE_ITO["vermilion"],
               mew=1.2, label="90% recovery vs the star's own residuals")],
        loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=1,
        fontsize=5.6)
    ax.set_title("(a) the two nulls, quiescent runs only", fontsize=7,
                 loc="left")

    fl = read_rows(cv, """SELECT * FROM p4_flicker
                          WHERE sf_target IS NOT NULL ORDER BY tau_s""")
    by_night: dict[str, list[dict]] = {}
    for r in fl:
        by_night.setdefault(f"{r['night']}|{r['filter']}", []).append(r)
    for key, rows in sorted(by_night.items()):
        band = key.split("|")[-1]
        tau = np.array([r["tau_s"] for r in rows], dtype=float)
        exc = np.array([(r["sf_excess"] if r["sf_excess"] is not None
                         else np.nan) for r in rows], dtype=float) * 1000.0
        det = np.array([bool(r["detected"]) for r in rows])
        ax2.plot(tau[det], exc[det], "-o", lw=0.6, ms=2.0,
                 color=BAND_COLOR.get(band, OKABE_ITO["grey"]), alpha=0.7,
                 mec="none")
        ax2.plot(tau[~det], exc[~det], "x", ms=2.6,
                 color=BAND_COLOR.get(band, OKABE_ITO["grey"]), alpha=0.5,
                 mew=0.7)
    floor = np.array([r["sf_floor"] for r in fl if r["sf_floor"]],
                     dtype=float) * 1000.0
    tauf = np.array([r["tau_s"] for r in fl if r["sf_floor"]], dtype=float)
    if floor.size:
        o = np.argsort(tauf)
        ax2.plot(tauf[o], floor[o], ".", ms=1.4, color="#bbbbbb", zorder=0,
                 label="measured floor (field stars)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("timescale $\\tau$ (s)")
    ax2.set_ylabel("structure-function amplitude (mmag)")
    ax2.grid(color="#f2f2f2")
    ax2.legend(loc="upper left", fontsize=5.4)
    ax2.set_title("(b) flickering over a measured floor", fontsize=7,
                  loc="left")

    # How many scopes the hump actually clears the INSTRUMENTAL contour on.
    # The panel already draws this correctly -- the I-band High Gain row
    # has its filled marker to the LEFT of the green tick -- and an earlier
    # revision's caption said "everywhere", contradicting its own figure
    # and erasing the one scope on which a hump of the fitted size could
    # not have been seen.
    above = [r for r in runs if r["hump_amp"] is not None
             and r["amp90_field"] is not None
             and r["hump_amp"] > r["amp90_field"]]
    n_above = len(above)
    below = next((r for r in runs if r not in above), None)

    spec = FigureSpec(
        fig_id="fig11", label="fig:yzcncfallback",
        title="YZ Cnc: orbital-hump limits and flickering",
        caption=(
            "(a) For each quiescent dense run, the fitted orbital-hump "
            "semi-amplitude against TWO injection--recovery contours "
            "measured on that run's own timestamps: against "
            "magnitude-matched field stars seen through the same frames, "
            "and against the star's own night-rolled residuals, which carry "
            f"its flickering. The hump exceeds the instrumental contour on "
            f"{n_above} of the {len(runs)} scopes and the red-noise "
            f"contour on none, so it is reported as an upper limit and not "
            f"a detection. On the remaining scope "
            f"({_scope_label(below)}) the fitted hump sits BELOW the "
            f"instrumental contour, so that run could not have shown a "
            f"hump of the fitted size at all and is uninformative rather "
            "than a non-detection. (b) "
            "Flickering as a structure function against timescale, with the "
            "photometric floor measured on the same field stars and "
            "subtracted in quadrature; crosses are bins that did not clear "
            "3$\\sigma$ and are reported as NOT MEASURED, never as zero."),
        tables=("p4_run", "p4_flicker"),
        width_in=COL_DOUBLE,
        note=("This is the strategy's own conditional branch, not a "
              "substitution: §7 item 11 specifies the superhump analysis "
              "'else orbital-hump/flickering statistics', and CV-S7 "
              "measured that no dense run falls inside a superoutburst."))
    return fig, spec


def _detrend_gap_direction(by_key) -> tuple:
    """``(median ratio, series where detrending is worse, series compared)``.

    Pairs each series' ``season`` contour against its ``season-dt`` one at
    the injected periods both carry, and returns the median of the
    detrended-over-raw amplitude ratio across series.  Below one means
    detrending RECOVERS a smaller amplitude, which is better sensitivity
    and the opposite of the "cost" Figure 12(a)'s caption used to assert.
    """
    raw: dict[str, dict] = {}
    dtd: dict[str, dict] = {}
    for (scope, score), rr in by_key.items():
        if score != "known":
            continue
        regime = rr[0]["regime"]
        if regime not in ("season", "season-dt"):
            continue
        book = raw if regime == "season" else dtd
        book[rr[0]["series_key"]] = {r["period_d"]: r["amp90"] for r in rr
                                     if r["amp90"]}
    meds = []
    for sk in sorted(set(raw) & set(dtd)):
        shared = sorted(set(raw[sk]) & set(dtd[sk]))
        if shared:
            meds.append(float(np.median([dtd[sk][p] / raw[sk][p]
                                         for p in shared])))
    if not meds:
        return (float("nan"), 0, 0)
    return (float(np.median(meds)), int(sum(m > 1.0 for m in meds)),
            len(meds))


def fig12_injection(ch, cv):
    """Figure 12 -- 90 per cent recovery contours at the real timestamps."""
    fig, axes = plt.subplots(1, 2, figsize=(COL_DOUBLE, 2.9), sharey=True)
    ax, ax2 = axes

    # ch_contour carries TWO scores per scope, and they answer different
    # questions: ``known`` is recovery when the period is known in advance,
    # ``period`` is recovery from a blind search.  The first version of
    # this figure plotted both under one label, which drew each contour as
    # a sawtooth between the two and made the whole panel meaningless.
    rows = read_rows(ch, """SELECT * FROM ch_contour ORDER BY scope, score,
                            period_d""")
    by_key: dict[tuple, list[dict]] = {}
    for r in rows:
        by_key.setdefault((r["scope"], r["score"]), []).append(r)
    tgt_color = {"stlmi": OKABE_ITO["blue"], "vvpup": OKABE_ITO["vermilion"],
                 "euuma": OKABE_ITO["green"], "anuma": OKABE_ITO["purple"],
                 "yzcnc": OKABE_ITO["orange"]}

    def _draw(target_ax, want_regime, want_score, style, label_fmt,
              shade=True):
        for (scope, score), rr in sorted(by_key.items()):
            if rr[0]["regime"] != want_regime or score != want_score:
                continue
            tgt = rr[0]["series_key"].split("|")[0]
            p = np.array([r["period_d"] for r in rr], dtype=float)
            a = np.array([r["amp90"] for r in rr], dtype=float) * 1000.0
            lo = np.array([(r["amp90_lo"] or np.nan)
                           for r in rr], dtype=float) * 1000.0
            hi = np.array([(r["amp90_hi"] or np.nan)
                           for r in rr], dtype=float) * 1000.0
            o = np.argsort(p)
            col = tgt_color.get(tgt, "k")
            target_ax.plot(p[o] * 24.0, a[o], style, lw=1.0, color=col,
                           label=label_fmt.format(
                               t=TARGET_LABEL.get(tgt, tgt)))
            good = np.isfinite(lo[o]) & np.isfinite(hi[o])
            if shade and good.sum() > 3:
                target_ax.fill_between(p[o][good] * 24.0, lo[o][good],
                                       hi[o][good], color=col, alpha=0.10,
                                       lw=0)

    # (a) THE TRADE detrending makes, not its cost.  The panel was labelled
    # "the cost of detrending" and the caption told a reader that the gap
    # between a pair of curves IS that cost.  Measured, the dashed curves
    # mostly sit BELOW the solid ones -- a smaller amplitude recovered 90
    # per cent of the time, i.e. better sensitivity -- because the filter
    # removes red noise along with signal.  The 25--123 per cent injected-
    # signal loss that forbids filtering is a different measurement
    # (p3_detrend), and the two point opposite ways.  The direction is
    # measured here and stated, never left to the eye.
    _draw(ax, "season", "known", "-", "{t} (raw)")
    _draw(ax, "season-dt", "known", "--", "{t} (detrended)", shade=False)
    _dt_med, _dt_worse, _dt_n = _detrend_gap_direction(by_key)
    # (b) what a single night can do, and what a BLIND search costs on it.
    _draw(ax2, "night", "known", "-", "{t} (period known)")
    _draw(ax2, "night", "period", "--", "{t} (blind search)", shade=False)

    for a_ in (ax, ax2):
        a_.set_xscale("log")
        a_.set_yscale("log")
        a_.set_xlabel("injected period (h)")
        a_.grid(color="#f2f2f2")
        a_.legend(fontsize=5.0, loc="upper left", ncol=2,
                  columnspacing=0.8, handletextpad=0.4, handlelength=1.6)
    ax.set_ylabel("semi-amplitude at 90% recovery (mmag)")
    ax.set_title("(a) whole season: raw against detrended", fontsize=7,
                 loc="left")
    ax2.set_title("(b) one night: the cost of a blind search", fontsize=7,
                  loc="left")

    ob = read_rows(cv, """SELECT amp90_blind, superhump_floor FROM
                          p4_outburst WHERE amp90_blind IS NOT NULL""")
    if ob:
        blind = 1000.0 * np.array([r["amp90_blind"] for r in ob], dtype=float)
        ax2.axhspan(float(np.min(blind)), float(np.max(blind)),
                    color=OKABE_ITO["grey"], alpha=0.18, lw=0, zorder=0)
        sh = [r["superhump_floor"] for r in ob if r["superhump_floor"]]
        if sh:
            ax2.axhline(1000.0 * float(np.min(sh)), lw=0.8, ls="-.",
                        color=OKABE_ITO["black"])
            ax2.text(0.02, 0.06, "superhump semi-amplitudes start here",
                     transform=ax2.transAxes, fontsize=5.4)

    spec = FigureSpec(
        fig_id="fig12", label="fig:injection",
        title="Injection--recovery contours at the real timestamps",
        caption=(
            "Semi-amplitude at which an injected sinusoid is recovered 90 "
            "per cent of the time, as a function of injected period, "
            "computed by injecting into the ACTUAL timestamps and the "
            "actual residuals of each series rather than into a simulated "
            "cadence. (a) Whole-season injections with the period known in "
            "advance: solid for the raw series, dashed after detrending. "
            "The gap between a pair of curves is NOT a sensitivity cost, "
            "and this panel is the measurement that says so: the detrended "
            f"contour is the LOWER of the two --- a smaller amplitude "
            f"recovered 90 per cent of the time --- for "
            f"{_dt_n - _dt_worse} of the {_dt_n} series, at a median "
            f"detrended/raw ratio of {_dt_med:.2f}, because the filter "
            "removes red noise along with signal. What forbids filtering "
            "first is a different measurement, the 25--123 per cent of an "
            "injected signal a naive filter destroys at these periods "
            "(Section~\\ref{sec:analysis}); the two point opposite ways, "
            "and the joint fit is what avoids having to choose. (b) "
            "Single-night "
            "injections, solid when the period is known and dashed for a "
            "blind search over the same band; the gap is the price of not "
            "knowing where to look. The grey band is the blind-search "
            "contour measured on the YZ Cnc normal-outburst runs and the "
            "dash-dot line the lower edge of published superhump "
            "semi-amplitudes: the band lies above that line, which is what "
            "turns `no superhump period' from an absence into a "
            "measurement."),
        tables=("ch_contour", "p4_outburst"),
        width_in=COL_DOUBLE)
    return fig, spec


def fig13_timing_audit(man):
    """Figure 13 (appendix) -- the header-time audit and the clock validator."""
    fig, axes = plt.subplots(1, 3, figsize=(COL_DOUBLE, 2.7),
                             gridspec_kw={"width_ratios": [1.15, 1.30, 1.0],
                                          "wspace": 0.45})
    ax, ax2, ax3 = axes

    # (a) The whole question in one histogram.  If the header JD were a
    # MID-exposure time the ratio would pile up at 0.5; if it were the
    # exposure START it piles up at 0.  Plotting the raw difference against
    # exposure time -- the first version of this panel -- puts every point
    # on the zero line and lets the half-exposure guide dominate the axis,
    # which shows the answer far less clearly than counting does.
    rows = read_rows(man, """SELECT jd_minus_dateobs_s, exptime_s
                             FROM s3_header_audit
                             WHERE jd_minus_dateobs_s IS NOT NULL
                               AND exptime_s > 0""")
    if rows:
        ratio = np.array([r["jd_minus_dateobs_s"] / r["exptime_s"]
                          for r in rows], dtype=float)
        ratio = ratio[np.isfinite(ratio)]
        ax.hist(ratio, bins=np.linspace(-0.1, 0.6, 36),
                color=OKABE_ITO["blue"], alpha=0.85)
        ax.axvline(0.0, color=OKABE_ITO["black"], lw=0.9)
        ax.axvline(0.5, ls="--", lw=0.9, color=OKABE_ITO["vermilion"])
        ax.text(0.03, 0.97, "exposure\nSTART", transform=ax.transAxes,
                fontsize=5.6, va="top")
        ax.text(0.5, 0.55, "mid-exposure:\nwhat it is NOT", rotation=90,
                transform=ax.transAxes, fontsize=5.6, va="center",
                ha="right", color=OKABE_ITO["vermilion"])
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_xlabel("(JD$_\\mathrm{hdr}-$DATE-OBS) / $t_\\mathrm{exp}$")
        ax.set_ylabel(f"frames (of {len(rows):,} audited)")
    else:
        _empty_panel(ax, "no header audit rows")
    ax.grid(color="#f2f2f2")
    ax.set_title("(a) what the header JD is", fontsize=7, loc="left")

    # (b) The clock cards agree to the stamp resolution in every mode.  The
    # informative number is therefore not the median -- it is exactly zero
    # everywhere -- but the WORST frame, and how many exceed 100 ms.
    aud = read_rows(man, """SELECT * FROM s3_dateobs_audit
                            WHERE n_frames >= 50 ORDER BY n_frames DESC""")
    if aud:
        y = np.arange(len(aud))
        # Each mode's bar starts at ITS OWN timestamp resolution: the
        # 2026 blank-mode headers stamp to 1 s and every other mode to 1 ms,
        # so one shared floor would draw eight modes as if they were coarser
        # than they are.
        res = float(np.nanmin([r["stamp_resolution_s"] or 1e-3
                               for r in aud]))
        for i, r in enumerate(aud):
            own = float(r["stamp_resolution_s"] or 1e-3)
            worst = max(float(r["max_abs_s"] or 0.0), own)
            ax2.plot([own, worst], [i, i], "-", lw=2.2,
                     color=OKABE_ITO["skyblue"], alpha=0.7,
                     solid_capstyle="butt")
            ax2.plot([own], [i], "|", ms=5, color=OKABE_ITO["black"],
                     mew=0.8)
            ax2.plot([worst], [i], "o", ms=3.2, color=OKABE_ITO["blue"],
                     mec="none")
            if r["n_gt_100ms"]:
                ax2.text(worst * 1.6, i, f"{int(r['n_gt_100ms'])}",
                         fontsize=5.2, va="center",
                         color=OKABE_ITO["vermilion"])
        ax2.axvline(res, ls=":", lw=0.9, color=OKABE_ITO["black"])
        ax2.text(res * 1.4, len(aud) - 1.6, "stamp\nresolution",
                 fontsize=5.0, va="center", color=OKABE_ITO["grey"])
        ax2.axvline(0.1, ls="--", lw=0.8, color=OKABE_ITO["vermilion"])
        ax2.set_yticks(y)
        # Heavily abbreviated: this panel sits between two others, and the
        # full readout-mode strings are long enough to overrun both of its
        # neighbours whichever side they are placed on.  The mapping is in
        # MODE_SHORT and the instrument table gives the full names.
        ax2.set_yticklabels([MODE_SHORT.get(r["readoutm"],
                                            str(r["readoutm"])[:9])
                             for r in aud], fontsize=5.4)
        ax2.set_xscale("log")
        ax2.set_xlim(res * 0.4, 1e4)
        ax2.set_xlabel("worst |TELUT $-$ DATE-OBS| (s)")
    else:
        _empty_panel(ax2, "no DATE-OBS audit rows")
    ax2.grid(axis="x", color="#f2f2f2")
    ax2.set_title("(b) worst clock disagreement per mode", fontsize=7,
                  loc="left")

    ec = read_rows(man, "SELECT * FROM s3_clock_eclipses ORDER BY tag")
    ok = [r for r in ec if r["o_minus_c_s"] is not None]
    bad = [r for r in ec if r["o_minus_c_s"] is None]
    # 'global' is the pooled row.  When only one dated epoch was usable the
    # pooled row IS that epoch, and plotting both draws one measurement
    # twice -- which reads as two independent confirmations.
    dated = [r for r in ok if str(r["tag"]).lower() != "global"]
    if len(dated) == 1 and len(ok) == 2:
        ok = dated
    if ok:
        y = np.arange(len(ok))
        v = np.array([r["o_minus_c_s"] for r in ok], dtype=float)
        e = np.array([(r["o_minus_c_err_s"] or np.nan) for r in ok],
                     dtype=float)
        b = [r["clock_bound_s"] for r in ok if r["clock_bound_s"]]
        if b:
            bound = float(np.max(b))
            ax3.axvspan(-bound, bound, color=OKABE_ITO["grey"], alpha=0.16,
                        lw=0)
            ax3.text(0.03, 0.95,
                     f"bound $\\pm${bound:.0f} s\n(weak: it is the "
                     f"validator's,\nnot the photometry's)",
                     transform=ax3.transAxes, fontsize=5.2, va="top")
        ax3.errorbar(v, y, xerr=e, fmt="o", ms=3.4, lw=0, elinewidth=1.0,
                     capsize=1.8, color=OKABE_ITO["green"])
        ax3.axvline(0.0, color="k", lw=0.7)
        ax3.set_yticks(y)
        ax3.set_yticklabels([str(r["tag"]) for r in ok], fontsize=5.4)
        ax3.set_ylim(-0.8, max(0.8, len(ok) - 0.2))
        ax3.set_xlabel("eclipse O$-$C (s)")
        if bad:
            ax3.text(0.5, -0.28,
                     f"{len(bad)} further epoch(s) not usable: "
                     + "; ".join(str(r["status"]) for r in bad),
                     transform=ax3.transAxes, ha="center", va="top",
                     fontsize=5.0, color=OKABE_ITO["grey"])
    else:
        _empty_panel(ax3, "no clock-validation eclipses with a usable O-C")
    ax3.grid(axis="x", color="#f2f2f2")
    ax3.set_title("(c) the clock validator", fontsize=7, loc="left")

    spec = FigureSpec(
        fig_id="fig13", label="fig:timingaudit",
        title="Timing audit (appendix)",
        caption=(
            "(a) The header JD card minus DATE-OBS, divided by the exposure "
            "time, for every audited frame. The distribution piles up at "
            "zero and not at one half, which is the measurement that the "
            "header JD is exposure START in UTC and must never be used as a "
            "mid-exposure time. Every time in this paper is a "
            "BJD$_\\mathrm{TDB}$ recomputed from DATE-OBS, the exposure "
            "time and the target coordinates. (b) The WORST disagreement "
            "between the two clock cards, TELUT and DATE-OBS, in each "
            "readout mode, on a logarithmic axis running from the "
            "timestamp resolution (dotted) to the 100 ms line (dashed); "
            "the red count is how many frames exceed 100 ms. Medians are "
            "not plotted because they are exactly zero in every mode. (c) "
            "The independent clock validator: eclipse timings of a detached "
            "eclipsing binary observed with the same instrument. Its bound "
            "is weak, and it is quoted as a bound rather than as a "
            "correction for exactly that reason: the timing evidence this "
            "paper relies on is panels (a) and (b)."),
        tables=("s3_header_audit", "s3_dateobs_audit", "s3_clock_eclipses"),
        width_in=COL_DOUBLE)
    return fig, spec


#: The registry the CLI iterates.  ``needs`` names which databases the
#: builder wants, so a missing product produces a clear message instead of
#: a traceback three frames deep.
BUILDERS: dict[str, dict] = {
    "fig01": {"fn": fig01_coverage, "needs": ("cv",)},
    "fig02": {"fn": fig02_rms_vs_mag, "needs": ("ch", "cv", "man")},
    "fig03": {"fn": fig03_linearity, "needs": ("man",)},
    "fig04": {"fn": fig04_periodograms, "needs": ("cv", "ch")},
    "fig05": {"fn": fig05_stlmi_folds, "needs": ("cv",)},
    "fig06": {"fn": fig06_colour_phase, "needs": ("cv",)},
    "fig07": {"fn": fig07_vvpup_euuma_folds, "needs": ("cv",)},
    "fig08": {"fn": fig08_state_history, "needs": ("cv",)},
    "fig09": {"fn": fig09_oc, "needs": ("cv",)},
    "fig10": {"fn": fig10_yzcnc_season, "needs": ("cv",)},
    "fig11": {"fn": fig11_yzcnc_fallback, "needs": ("cv",)},
    "fig12": {"fn": fig12_injection, "needs": ("ch", "cv")},
    "fig13": {"fn": fig13_timing_audit, "needs": ("man",)},
}

FIGURE_IDS: tuple[str, ...] = tuple(sorted(BUILDERS))
