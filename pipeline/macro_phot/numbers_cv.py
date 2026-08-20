"""CV-S11 — ``numbers.tex``: every value in the manuscript, emitted by query.

THE LAW THIS MODULE ENFORCES
----------------------------
The web reports under ``docs/CV_TimeSeries/`` obey one rule: no number on a
published page was typed by a person.  Each is the result of a query the
report script executed.  A manuscript is a published page with a different
renderer, so it obeys the same rule -- and this module is how.

It reads the products databases and writes a single LaTeX file of
``\\newcommand`` definitions.  ``main.tex`` does ``\\input{numbers}`` and
then writes ``\\NumStLMiNights`` where a person would otherwise have typed
``20``.  The consequence is the point: when a pipeline stage is re-run and
a number moves, the paper moves with it on the next ``tectonic`` build, and
a number that becomes unmeasurable becomes a compile-time error rather than
a stale claim in print.

WHY A MACRO FILE AND NOT A TABLE-GENERATOR
-------------------------------------------
Tables are easy to emit and easy to keep honest.  PROSE is where numbers rot
-- "we obtained 3,157 frames" survives three re-reductions without anyone
noticing it now says something false.  Macros put the prose under the same
discipline as the tables.

HOW A NUMBER GETS IN HERE
--------------------------
Through :func:`collect`, which is a list of ``Number`` records.  Each one
names the macro, the SQL that produced it, the unit, and the table it came
from.  ``p5_number`` stores all four, so a referee asking "where does 25
mmag come from?" is answered by one query against the product, not by
reading the paper's source.

FORMATTING RULES
----------------
* Integers get a thin space as the thousands separator (``3\\,157``), which
  is the AAS style and does not break across a line.
* A quantity with an error bar is emitted as TWO macros (``...``, ``...Err``)
  and never as a pre-formatted ``$a \\pm b$`` string: the paper decides how
  to typeset it, this module decides what it is.
* A missing value emits the macro anyway, expanding to
  ``\\NumMissing`` -- which ``main.tex`` defines as a visible marker.  A
  macro that silently vanishes would let a sentence lose its number and
  still compile.
"""

from __future__ import annotations

import cmath
import math
import re
import sqlite3
from dataclasses import dataclass
from typing import Any, Optional, Sequence

# ===========================================================================
# Pure formatting -- unit-testable without a database
# ===========================================================================

#: Digit -> word, because a LaTeX control sequence may contain LETTERS ONLY.
#: ``\Num2024Frames`` is not a macro name; ``\NumTwoZeroTwoFourFrames`` is.
_DIGIT_WORD = {"0": "Zero", "1": "One", "2": "Two", "3": "Three",
               "4": "Four", "5": "Five", "6": "Six", "7": "Seven",
               "8": "Eight", "9": "Nine"}


def tex_macro_name(key: str, prefix: str = "Num") -> str:
    """``stlmi_e7_nights`` -> ``NumStlmiESevenNights``.

    A LaTeX control sequence is letters only, so digits become their English
    names and every other character is a word separator.  The mapping is
    deterministic and injective on the keys this module uses, which is what
    lets the same key be looked up from a test.
    """
    parts = [p for p in re.split(r"[^A-Za-z0-9]+", str(key)) if p]
    out = []
    for p in parts:
        chunk = "".join(_DIGIT_WORD[c] if c.isdigit() else c for c in p)
        out.append(chunk[:1].upper() + chunk[1:])
    name = prefix + "".join(out)
    if not name.isalpha():
        raise ValueError(f"macro name {name!r} from key {key!r} is not "
                         f"letters-only; LaTeX will reject it")
    return name


def fmt_int(value: Any) -> Optional[str]:
    """Integer with a LaTeX thin space every three digits, or None.

    ``\\,`` rather than a comma: AAS style, and it cannot be mistaken for a
    decimal point by a reader whose locale uses one.
    """
    if value is None:
        return None
    try:
        n = int(round(float(value)))
    except (TypeError, ValueError):
        return None
    s = f"{abs(n):,}".replace(",", "\\,")
    return ("-" if n < 0 else "") + s


def fmt_float(value: Any, nd: int = 2) -> Optional[str]:
    """Fixed-point float, or None for anything that is not a finite number.

    NaN and infinity return None on purpose.  A NaN printed into a paper as
    "nan" is a typesetting bug; a NaN that becomes the missing-value marker
    is a visible reminder that something was not measured.
    """
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return f"{v:.{nd}f}"


def fmt_sci(value: Any, nd: int = 2) -> Optional[str]:
    """LaTeX scientific notation, e.g. ``5.0 \\times 10^{-9}``."""
    if value is None:
        return None
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v) or v == 0.0:
        return fmt_float(v, nd)
    exp = int(math.floor(math.log10(abs(v))))
    mant = v / (10.0 ** exp)
    return f"{mant:.{nd}f} \\times 10^{{{exp}}}"


def fmt_range(lo: Any, hi: Any, nd: int = 0, dash: str = "--") -> Optional[str]:
    """``9--77``: the honest shape for a spread this project measures a lot.

    A range is one macro, not two, because the paper never wants to split
    them and a split invites one half being updated without the other.

    A degenerate range collapses to the value: every block holding out four
    check stars must print "4" and not "4--4", which reads as a spread the
    data do not have.
    """
    a, b = fmt_float(lo, nd), fmt_float(hi, nd)
    if a is None or b is None:
        return None
    return a if a == b else f"{a}{dash}{b}"


def fmt_percent(value: Any, nd: int = 0) -> Optional[str]:
    """A fraction in [0, 1] rendered as a per-cent number without the sign.

    The ``\\%`` is left to the prose: some sentences want "31 per cent".
    """
    v = fmt_float(None if value is None else 100.0 * float(value), nd)
    return v


@dataclass(frozen=True)
class Number:
    """One value the manuscript is allowed to state.

    ``key`` is the stable name; ``value`` the already-formatted LaTeX body;
    ``unit`` the unit it is in (empty for counts and fractions); ``source``
    the database table it came from; ``note`` the one clause a referee needs
    to know about how it was derived.  All four go into ``p5_number``.
    """

    key: str
    value: Optional[str]
    unit: str = ""
    source: str = ""
    note: str = ""

    @property
    def macro(self) -> str:
        return tex_macro_name(self.key)

    @property
    def body(self) -> str:
        return self.value if self.value is not None else "\\NumMissing"


def render_tex(numbers: Sequence[Number], stamp: str = "") -> str:
    """The ``numbers.tex`` body: one ``\\newcommand`` per Number.

    ``\\NumMissing`` is defined first, and loudly: a value the database
    could not supply appears in the PDF as a marker a proof-reader cannot
    miss, rather than as a gap or a stale constant.
    """
    seen: dict[str, str] = {}
    lines = [
        "%% numbers.tex -- GENERATED FILE.  DO NOT EDIT.",
        "%% Emitted by pipeline/macro_phot/numbers_cv.py from the CV "
        "products database.",
        "%% Every value in the manuscript prose comes from here, so that no "
        "number is typed",
        "%% by hand and a pipeline re-run propagates into the paper on the "
        "next build.",
    ]
    if stamp:
        lines.append(f"%% {stamp}")
    lines += [
        "",
        "%% What a value the database could not supply looks like in print.",
        "\\providecommand{\\NumMissing}{\\textbf{[NUMBER MISSING]}}",
        "",
    ]
    for n in numbers:
        macro = n.macro
        if macro in seen:
            raise ValueError(f"duplicate macro {macro} from keys "
                             f"{seen[macro]!r} and {n.key!r}")
        seen[macro] = n.key
        tail = f"  % {n.unit}" if n.unit else ""
        if n.source:
            tail += f" [{n.source}]" if tail else f"  % [{n.source}]"
        lines.append(f"\\newcommand{{\\{macro}}}{{{n.body}}}{tail}")
    lines.append("")
    return "\n".join(lines)


# ===========================================================================
# Small query helpers
# ===========================================================================
def one(con: sqlite3.Connection, sql: str, args: Sequence = ()) -> Any:
    """First column of the first row, or None."""
    r = con.execute(sql, tuple(args)).fetchone()
    return None if r is None else r[0]


def rows(con: sqlite3.Connection, sql: str, args: Sequence = ()) -> list[dict]:
    cur = con.execute(sql, tuple(args))
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


def _minmax(values, scale: float = 1.0):
    """Finite min and max of an iterable, scaled.  ``(None, None)`` if empty."""
    v = [float(x) * scale for x in values
         if x is not None and math.isfinite(float(x))]
    return (min(v), max(v)) if v else (None, None)


def _median(values, scale: float = 1.0):
    """Finite median of an iterable, scaled.  ``None`` if empty.

    ``statistics.median`` averages the two central values of an even-length
    sample.  Taking ``sorted(v)[len(v)//2]`` instead -- an earlier revision
    did -- makes the macro file and the figure disagree by a second on a
    144-cell injection grid, which is exactly the kind of drift this module
    exists to prevent.
    """
    v = sorted(float(x) * scale for x in values
               if x is not None and math.isfinite(float(x)))
    if not v:
        return None
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


#: The qualifying rule for a "full-orbit night", fixed by CV-S10 and stored
#: in ``p4_meta.full_orbit_min_points``: at least this many cloud-clean,
#: unsaturated, catalogue-tied target points on that night in that series,
#: spanning more than one orbital period.  It is quoted here so that the
#: census in §2.2 is reproducible from the release rather than asserted.
FULL_ORBIT_MIN_POINTS_DEFAULT = 12

#: Instrument-era and target labels, shared by the macro collector and the
#: table renderers.  Defined up here rather than beside the renderers
#: because ``collect`` names an era in a macro's note and nothing should
#: depend on the order two module-level dicts happen to appear in.
ERA_LABEL = {6: "High Gain StackPro", 7: "High Gain",
             47: "1MHz HS 16-bit", 72: "1MHz HS 16-bit",
             76: "Mode0", 78: "Fast", 79: "Fast"}

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}

def pdot_timescale_yr(period_d: Any, pdot: Any) -> Optional[float]:
    """``P / |dP/dt|`` in years, with the PERIOD in the numerator.

    This existed inline in the macro table and was written as
    ``1 / (pdot * 365.25)``, which drops ``P`` entirely.  ``pdot`` is
    dimensionless (days per day), so that expression is not a time at all;
    it printed 7.6e5 yr in the abstract and a conclusion where the paper's
    own P and Pdot give 6.0e4 yr, a factor 12.6, in the paper whose thesis
    is that a script-emitted number cannot go stale.

    It is a named function rather than an expression so that the regression
    test in ``test_phase3.py`` can hold the emitter itself to the identity
    --- halving the period halves the timescale --- instead of skipping.
    """
    try:
        p, d = float(period_d), float(pdot)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(p) and math.isfinite(d)) or d == 0.0:
        return None
    return p / abs(d) / 365.25


def dynamic_range_ratios(params: dict) -> tuple:
    """``(min, max)`` of each 16-bit mode's ceiling over High Gain's.

    ``params`` maps ``(era_group, quantity)`` to a float, as
    ``detector_params`` supplies it.  Returns ``(None, None)`` when either
    side is missing.

    THE REASON THIS IS A FUNCTION.  §2.1 called the 2024/2025 dynamic-range
    step "nearly twenty" and Figure 3's caption "a factor of 16" -- the
    first being Mode0 over High Gain (18.7) and the second High Gain
    StackPro over High Gain (16.0), which also happens to equal the NOMINAL
    bit-depth ratio 2^16/2^12 and so read as though it had been assumed
    rather than measured.  Both are correct and neither is the whole
    answer, because "every 16-bit mode" is not one ceiling: the measured
    16-bit clips run 56,062 to 65,535 ADU.  The paper now quotes the RANGE,
    from here, in both places.
    """
    hg = params.get(("High Gain", "ceiling_adu"))
    sixteen = [v for (g, q), v in params.items()
               if q == "ceiling_adu" and params.get((g, "adc_bits")) == 16.0]
    if not hg or not sixteen:
        return (None, None)
    return (min(sixteen) / hg, max(sixteen) / hg)


#: WHAT THE HOLD-OUT RULE COVERS, WRITTEN ONCE.
#:
#: §3.1 was rewritten to retract the claim that the headline precision is a
#: held-out measurement: it is a local fit of the scatter--magnitude
#: relation, over ensemble AND check stars together, evaluated at the CV's
#: magnitude.  Figure 2's caption -- the caption of the very panel those
#: numbers are plotted on -- went on asserting the retracted rule for
#: another revision, and asserted it of the noise floor as well, which is
#: fitted over the same mixed population.  The clause therefore lives here,
#: is emitted as ``\NumPrecisionScopeClause`` for §3.1, and is imported by
#: :mod:`macro_phot.figures_cv` for the caption.  One string, two renderers.
PRECISION_SCOPE_CLAUSE = (
    "The crosses at the CV's own magnitude and the noise floor annotated on "
    "each panel are both LOCAL FITS over the ensemble and check stars "
    "together, not held-out statistics: the held-out quantities in this "
    "paper are the check-star scatter $\\sigma_{\\rm chk}$, the "
    "catalogue-tie accuracy and the error-bar inflation factor")

#: The four panels of Figure 6, as ``(era_id, band_a, band_b)``.  Defined
#: here and imported by :mod:`macro_phot.figures_cv` so that the tie bars
#: §3.3 quotes in the text are computed over exactly the panels the figure
#: draws.  A referee found the text describing a 25--64 mmag choice while
#: the headline colour panel carried a 200 mmag bar; that can only be
#: prevented by the text and the panel reading one list.
COLOUR_PANEL_PAIRS = ((7, "G", "R"), (7, "R", "I"),
                      (76, "g", "r"), (76, "r", "i"))

#: Which filters count towards a "three-filter" night.  The 2024 High Gain
#: era wrote G/R/I and the 2025 Sloan era g/r/i for the same three
#: bandpass slots, and a night lies entirely inside one era, so the two
#: spellings are the same slot and are folded together here.
COLOUR_BANDS = ("g", "r", "i")


def full_orbit_nights(cv: sqlite3.Connection, series_key: str,
                      period_d: float, min_points: int) -> set:
    """The nights on which ONE series covered a full orbit.

    Exactly CV-S10's rule (``run_cv_final.full_orbit_nights``), repeated
    here so the manuscript's census is a query against the released
    database rather than four numbers typed into a source file.  Cloud-
    vetoed frames, saturated detections and untied points are excluded,
    because a night whose orbit is covered only by frames the photometry
    threw away did not cover an orbit.
    """
    by: dict[str, list[float]] = {}
    for night, t in cv.execute("""
            SELECT f.night, l.bjd_tdb
            FROM cv_lightcurve l
            JOIN cv_frames f ON f.frame_id = l.frame_id
                            AND f.series_key = l.series_key
            LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                                      AND c.series_key = l.series_key
            WHERE l.series_key = ? AND l.role = 'target'
              AND l.cal_mag IS NOT NULL AND l.saturated = 0
              AND COALESCE(c.vetoed, 0) = 0""", (series_key,)):
        by.setdefault(str(night), []).append(float(t))
    return {n for n, v in by.items()
            if len(v) >= min_points and (max(v) - min(v)) > period_d}


def _folded_profile(cv: sqlite3.Connection, series_key: str,
                    target_key: str, n_bins: int = 40, min_count: int = 3):
    """The binned folded light curve of one series: ``(centres, medians)``.

    The same statistic Figure 5 draws, computed here so that the one
    qualitative sentence about the fold's SHAPE carries a number a reader
    can check instead of an impression.  ``None`` where the fold is not
    defined or too sparse to bin.
    """
    eph = cv.execute("SELECT period_d, epoch_bjd FROM p3_ephemeris "
                     "WHERE target_key=?", (target_key,)).fetchone()
    if eph is None or eph[0] is None or eph[1] is None:
        return None
    period, epoch = float(eph[0]), float(eph[1])
    pts = cv.execute("""
        SELECT l.bjd_tdb, l.cal_mag FROM cv_lightcurve l
        LEFT JOIN p2_cloud_frame c ON c.frame_id = l.frame_id
                                  AND c.series_key = l.series_key
        WHERE l.series_key = ? AND l.role='target'
          AND l.cal_mag IS NOT NULL AND l.saturated = 0
          AND COALESCE(c.vetoed, 0) = 0""", (series_key,)).fetchall()
    if len(pts) < n_bins * min_count:
        return None
    bins: list[list[float]] = [[] for _ in range(n_bins)]
    for t, m in pts:
        ph = ((float(t) - epoch) / period) % 1.0
        bins[min(int(ph * n_bins), n_bins - 1)].append(float(m))
    med = [(sorted(b)[len(b) // 2] if len(b) >= min_count else None)
           for b in bins]
    if sum(x is not None for x in med) < n_bins - 4:
        return None
    # A gap would make the gradient of the neighbouring bins meaningless,
    # so fill the few empty bins by linear interpolation on the circle.
    for i, x in enumerate(med):
        if x is None:
            prev = next((med[(i - k) % n_bins] for k in range(1, n_bins)
                         if med[(i - k) % n_bins] is not None), None)
            nxt = next((med[(i + k) % n_bins] for k in range(1, n_bins)
                        if med[(i + k) % n_bins] is not None), None)
            med[i] = (0.5 * (prev + nxt) if prev is not None
                      and nxt is not None else (prev or nxt))
    centres = [(i + 0.5) / n_bins for i in range(n_bins)]
    return centres, med


def coverage_census(cv: sqlite3.Connection, target_key: str,
                    min_points: int) -> dict:
    """The three-filter full-orbit night census for one target.

    Returns the set of nights covering a full orbit in ANY filter, the set
    covering one in all three of ``COLOUR_BANDS``, and the per-band counts.
    §2.2 calls this the census that decides which figures can exist, so it
    is emitted as a table (``\\NumXxxThreeFilterNights`` and friends) and
    the qualifying nights themselves are listed in ``p5_number``'s note.
    """
    p = one(cv, "SELECT period_d FROM p3_ephemeris WHERE target_key=?",
            (target_key,))
    if p is None:
        return {"any": set(), "three": set(), "per_band": {}}
    per_band: dict[str, set] = {}
    for sk, filt in cv.execute(
            "SELECT series_key, filter FROM cv_series WHERE target_key=? "
            "AND status='solved'", (target_key,)):
        per_band.setdefault(str(filt).lower(), set()).update(
            full_orbit_nights(cv, sk, float(p), min_points))
    have = [b for b in COLOUR_BANDS if b in per_band]
    three = (set.intersection(*[per_band[b] for b in have])
             if len(have) == len(COLOUR_BANDS) else set())
    any_night: set = set()
    for v in per_band.values():
        any_night |= v
    return {"any": any_night, "three": three, "per_band": per_band}


# ===========================================================================
# The collector.  One section per manuscript section, in the paper's order.
# ===========================================================================
def collect(cv: sqlite3.Connection, ch: sqlite3.Connection,
            man: sqlite3.Connection) -> list[Number]:
    """Every number the manuscript prose is allowed to state.

    Grouped in the order the paper uses them, so that a section of the paper
    and a block of this function can be read side by side.
    """
    N: list[Number] = []

    def add(key, value, unit="", source="", note=""):
        N.append(Number(key, value, unit, source, note))

    # -- §2 Observations: the census ------------------------------------
    add("targets", fmt_int(one(cv, "SELECT count(DISTINCT target_key) "
                               "FROM cv_series")),
        source="cv_series", note="targets with at least one staged series")
    add("series total", fmt_int(one(cv, "SELECT count(*) FROM cv_series")),
        source="cv_series")
    add("series solved", fmt_int(one(cv, "SELECT count(*) FROM cv_series "
                                     "WHERE status='solved'")),
        source="cv_series", note="series with a converged ensemble solution")
    add("frames total", fmt_int(one(cv, "SELECT count(*) FROM cv_frames")),
        source="cv_frames", note="staged light frames, raw-tree alias-merged")
    add("nights total", fmt_int(one(cv, "SELECT count(DISTINCT night) "
                                    "FROM cv_frames")), source="cv_frames")
    # CATALOGUE-TIED target measurements.  The bare ``role='target'`` count
    # includes detections in series whose tie did not converge, and the
    # per-target macros below already exclude those, so quoting the bare
    # count as the total made the five per-target numbers fail to add up.
    add("points total", fmt_int(one(cv, "SELECT count(*) FROM cv_lightcurve "
                                    "WHERE role='target' AND cal_mag IS NOT "
                                    "NULL")),
        source="cv_lightcurve",
        note="catalogue-tied target measurements, all series; the five "
             "per-target counts sum to exactly this")
    # UNTIED TARGET POINTS, SPLIT BY WHY.  An earlier revision described all
    # of them as sitting in series "whose tie did not converge", which is
    # true of only one of the three blocks: the other two ARE tied and
    # appear in Table 2, and their untied rows are individual detections
    # that the converged solution could not place on the standard system.
    untied = rows(cv, """
        SELECT l.series_key AS sk, count(*) AS n,
               COALESCE(c.verdict, 'no tie') AS verdict
        FROM cv_lightcurve l
        LEFT JOIN cv_cattie c ON c.series_key = l.series_key
                             AND c.is_primary = 1
        WHERE l.role='target' AND l.cal_mag IS NULL
        GROUP BY 1, 3 ORDER BY n DESC""")
    n_in_tied = sum(r["n"] for r in untied
                    if str(r["verdict"]).startswith("TIED"))
    n_in_untied = sum(r["n"] for r in untied
                      if not str(r["verdict"]).startswith("TIED"))
    add("points untied", fmt_int(sum(r["n"] for r in untied)),
        source="cv_lightcurve",
        note="target detections carrying no catalogue-tied magnitude: "
             + "; ".join(f"{r['sk']} {r['n']} ({r['verdict']})"
                         for r in untied))
    add("points untied in tied blocks", fmt_int(n_in_tied),
        source="cv_lightcurve",
        note="detections inside blocks whose tie DID converge and which "
             "appear in Table 2: individual points the fitted colour "
             "relation could not place on the standard system, not the "
             "output of a failed solve")
    add("points untied in untied blocks", fmt_int(n_in_untied),
        source="cv_lightcurve",
        note="detections in the one block with no usable tie at all")
    add("lightcurve rows", fmt_int(one(cv, "SELECT count(*) FROM "
                                       "cv_lightcurve")),
        source="cv_lightcurve",
        note="all roles: target, comparison and check")

    for tgt, tag in (("stlmi", "st lmi"), ("vvpup", "vv pup"),
                     ("euuma", "eu uma"), ("anuma", "an uma"),
                     ("yzcnc", "yz cnc")):
        add(f"{tag} frames", fmt_int(one(
            cv, "SELECT count(*) FROM cv_frames WHERE target_key=?", (tgt,))),
            source="cv_frames")
        add(f"{tag} nights", fmt_int(one(
            cv, "SELECT count(DISTINCT night) FROM cv_frames "
                "WHERE target_key=?", (tgt,))), source="cv_frames")
        add(f"{tag} points", fmt_int(one(
            cv, "SELECT count(*) FROM cv_lightcurve l JOIN cv_series s "
                "ON s.series_key=l.series_key WHERE s.target_key=? "
                "AND l.role='target' AND l.cal_mag IS NOT NULL", (tgt,))),
            source="cv_lightcurve",
            note="catalogue-tied target measurements")
        p = one(cv, "SELECT period_d FROM p3_ephemeris WHERE target_key=?",
                (tgt,))
        add(f"{tag} period d", fmt_float(p, 8), unit="d",
            source="p3_ephemeris", note="published VSX period, not measured")
        add(f"{tag} period min", fmt_float(None if p is None
                                           else float(p) * 1440.0, 1),
            unit="min", source="p3_ephemeris")

    # Three-filter full-orbit night census: the number that decides which
    # colour panels may exist at all, and therefore the one number in this
    # paper that must not be typed.  CV-S10's rule, applied here to every
    # target rather than only to AN UMa, reproduces p4_anuma's own count.
    mp = one(cv, "SELECT value FROM p4_meta WHERE key='full_orbit_min_points'")
    min_pts = int(mp) if mp is not None else FULL_ORBIT_MIN_POINTS_DEFAULT
    add("full orbit min points", fmt_int(min_pts), source="p4_meta",
        note="points a night must carry, over more than one orbit, to count "
             "as covering a full orbit")
    for tgt, tag in (("stlmi", "st lmi"), ("anuma", "an uma"),
                     ("vvpup", "vv pup"), ("euuma", "eu uma")):
        cen = coverage_census(cv, tgt, min_pts)
        add(f"{tag} three filter nights", fmt_int(len(cen["three"])),
            source="cv_lightcurve",
            note="nights with all three of g, r and i independently "
                 "covering a full orbit: "
                 + (", ".join(sorted(cen["three"])) or "none"))
        add(f"{tag} full orbit nights", fmt_int(len(cen["any"])),
            source="cv_lightcurve",
            note="nights covering a full orbit in at least one filter; "
                 "per band "
                 + ", ".join(f"{b}={len(cen['per_band'].get(b, ()))}"
                             for b in COLOUR_BANDS))

    # -- §2 Detector constants (the restored S2) ------------------------
    def dp(group, quantity):
        return one(man, "SELECT value FROM detector_params WHERE "
                        "era_group=? AND quantity=?", (group, quantity))

    add("hg ceiling adu", fmt_int(dp("High Gain", "ceiling_adu")), unit="ADU",
        source="detector_params", note="12-bit pileup clip")
    add("hg veto adu", fmt_int(dp("High Gain", "saturation_veto_adu")),
        unit="ADU", source="detector_params", note="0.92 of the ceiling")
    add("hg read noise e", fmt_float(dp("High Gain", "read_noise_e"), 2),
        unit="e-", source="detector_params")
    add("hg read noise e err", fmt_float(one(
        man, "SELECT uncertainty FROM detector_params WHERE era_group=? "
             "AND quantity=?", ("High Gain", "read_noise_e")), 2),
        unit="e-", source="detector_params",
        note="dominated by the gain bracket, not by statistics")
    glo = dp("High Gain", "gain_lower_bound_e_per_adu")
    ghi = dp("High Gain", "gain_upper_bound_e_per_adu")
    add("hg gain lo", fmt_float(glo, 2), unit="e-/ADU",
        source="detector_params", note="sky-level light-pair PTC slope")
    add("hg gain hi", fmt_float(ghi, 2), unit="e-/ADU",
        source="detector_params", note="dark-pair PTC slope")
    add("hg gain bracket", fmt_range(glo, ghi, 2), unit="e-/ADU",
        source="detector_params",
        note="a BRACKET, never a value: the two bounds are biased in "
             "opposite directions and the flat-field PTC that closes them "
             "awaits the October 2026 restart")
    add("mode zero ceiling adu", fmt_int(dp("Mode0", "ceiling_adu")),
        unit="ADU", source="detector_params")
    add("mode zero veto adu", fmt_int(dp("Mode0", "saturation_veto_adu")),
        unit="ADU", source="detector_params")
    add("stackpro ceiling adu", fmt_int(dp("High Gain StackPro",
                                           "ceiling_adu")),
        unit="ADU", source="detector_params")
    add("nsub", fmt_int(dp("High Gain StackPro", "nsub")),
        source="detector_params",
        note="sub-exposures coadded per StackPro frame before readout")
    # Table 1's own rows, as macros, so §2.1's sentence about the pileup
    # measurement cannot describe a table it no longer matches.
    cm = rows(man, "SELECT * FROM s2_ceiling_modes WHERE clip_adu IS NOT NULL")
    nf_lo, nf_hi = _minmax([r["n_frames"] for r in cm])
    add("ceiling frames range", fmt_range(nf_lo, nf_hi, 0),
        source="s2_ceiling_modes",
        note="science frames whose pixel histograms the pileup clip was "
             "measured in, per readout mode")
    vr = [r["veto_adu"] / r["clip_adu"] for r in cm
          if r["clip_adu"] and r["veto_adu"]]
    vlo, vhi = _minmax(vr)
    add("veto clip ratio range", fmt_range(vlo, vhi, 3),
        source="s2_ceiling_modes",
        note="the veto is 0.92 of the measured clip ROUNDED DOWN to the "
             "nearest 100 ADU, so the realised ratio is a little below 0.92 "
             "in every mode")
    add("veto clip ratio nominal", fmt_float(0.92, 2),
        source="ANALYSIS_STRATEGY §2",
        note="the nominal fraction of the measured clip the photometry "
             "vetoes at; an external choice, not a measurement")
    # THE DYNAMIC-RANGE RATIO, MEASURED, ONCE.  §2.1 called it "nearly
    # twenty" and Figure 3's caption "a factor of 16" -- the second being
    # the nominal bit-depth ratio 2^16/2^12 and not a measured quantity at
    # all.  In a paper whose Table 1 caption says no threshold was chosen
    # by eye, the same comparison must not carry two numbers, and neither
    # of them may be typed.  Both bounds are emitted, because "every 16-bit
    # mode" is not one ceiling: the 16-bit clips run 56,062 (High Gain
    # StackPro) to 65,535 (Mode0).
    _dpar = {(r["era_group"], r["quantity"]): float(r["value"])
             for r in rows(man, "SELECT era_group, quantity, value FROM "
                                "detector_params")}
    _rlo, _rhi = dynamic_range_ratios(_dpar)
    add("dynamic range ratio min", fmt_float(_rlo, 1),
        source="detector_params",
        note="the SMALLEST 16-bit measured ceiling divided by High Gain's "
             "(High Gain StackPro): every 16-bit mode is at LEAST this far "
             "above High Gain. It coincides with the nominal bit-depth "
             "ratio 2^16/2^12 = 16 and is not that ratio -- it is measured")
    add("dynamic range ratio max", fmt_float(_rhi, 1),
        source="detector_params",
        note="the LARGEST 16-bit measured ceiling divided by High Gain's, "
             "i.e. Mode0 against High Gain: the widest dynamic-range step "
             "across the 2024/2025 instrument seam")
    add("dynamic range ratio range", fmt_range(_rlo, _rhi, 1),
        source="detector_params",
        note="the measured ratio of each 16-bit mode's ceiling to High "
             "Gain's, low to high. §2.1 and Figure 3's caption both quote "
             "it through numbers_cv.dynamic_range_ratios so they cannot "
             "disagree, and neither may say 'a factor of X' where the "
             "modes span a range")

    # -- §3 Photometry: ensemble, tie, error model ----------------------
    solved = rows(cv, "SELECT * FROM cv_series WHERE status='solved' "
                      "AND check_rms_median IS NOT NULL")
    # The precision that matters is the one reached AT THE CV's OWN
    # MAGNITUDE, not the median over a field whose stars are mostly
    # brighter or fainter than the target.  ch_noise_series measures it
    # directly; cv_series' check-star median is a FIELD statistic and is
    # quoted under its own name so the two cannot be confused.
    #
    # WHAT SAMPLE THIS IS FITTED OVER, SAID HERE BECAUSE §3.1 MUST SAY IT.
    # ``prec_at_target_all`` interpolates the measured RMS-magnitude
    # relation at the target's magnitude over the constant stars NEAR that
    # magnitude, and that local sample is comparison AND check stars --
    # overwhelmingly comparison stars, whose residuals the ensemble solve
    # minimises by construction.  It is therefore NOT a held-out statistic,
    # and the paper may not claim it as one.  The held-out statistics are
    # ``check_rms_median`` (Table 2's sigma_chk), the tie accuracy and the
    # inflation factor; ``ch_check_bias`` measures how far a four-star
    # hold-out can be from a magnitude-matched field sample, and that
    # spread is emitted below because it bounds what the hold-out could
    # have said instead.
    prec = [r["prec_at_target_all"] for r in
            rows(ch, "SELECT prec_at_target_all FROM ch_noise_series")]
    lo, hi = _minmax(prec, 1000.0)
    add("precision lo mmag", fmt_float(lo, 0), unit="mmag",
        source="ch_noise_series",
        note="per-point precision at the CV's own magnitude, best series")
    add("precision hi mmag", fmt_float(hi, 0), unit="mmag",
        source="ch_noise_series")
    add("precision range mmag", fmt_range(lo, hi, 0), unit="mmag",
        source="ch_noise_series",
        note="the paper's headline per-point precision: the measured "
             "RMS-magnitude relation evaluated at each target's own "
             "magnitude, fitted over the constant ENSEMBLE AND CHECK stars "
             "local to that magnitude. It is a local fit, not a held-out "
             "statistic, and §3.1 says so")
    # Only the series that actually CONTRIBUTE a precision value: the five
    # with no target magnitude contribute none, and counting their empty
    # neighbourhoods would open the range at zero.
    nnt = [r["n_near_target"] for r in
           rows(ch, "SELECT n_near_target FROM ch_noise_series "
                    "WHERE prec_at_target_all IS NOT NULL "
                    "AND n_near_target IS NOT NULL")]
    nlo, nhi = _minmax(nnt)
    add("precision sample range", fmt_range(nlo, nhi, 0),
        source="ch_noise_series",
        note="constant stars near the target's magnitude that the local "
             "fit uses, per series: comparison and check stars together, "
             "the upper end being VV Pup's crowded low-latitude field")
    # The scope clause itself, as a macro, so §3.1 and Figure 2's caption
    # are the same sentence rather than two sentences that have to be kept
    # in step by hand.  See PRECISION_SCOPE_CLAUSE.
    add("precision scope clause", PRECISION_SCOPE_CLAUSE,
        source="ch_noise_series",
        note="the scope of the hold-out rule, stated once and used by both "
             "§3.1 and Figure 2's caption")
    cb = rows(ch, "SELECT * FROM ch_check_bias")
    blo2, bhi2 = _minmax([r["bias_ratio"] for r in cb])
    add("check bias ratio range", fmt_range(blo2, bhi2, 2),
        source="ch_check_bias",
        note="median check-star RMS divided by the median RMS of "
             "magnitude-matched FIELD stars in the same frames, per block: "
             "how far a four-star hold-out can sit from the field it is "
             "meant to represent")
    add("check bias blocks", fmt_int(len(cb)), source="ch_check_bias")
    clo, chi = _minmax([r["check_rms_median"] for r in solved], 1000.0)
    add("check rms range mmag", fmt_range(clo, chi, 0), unit="mmag",
        source="cv_series",
        note="median held-out check-star scatter per series: a field "
             "statistic, not the precision at the target")
    ilo, ihi = _minmax([r["chi2_inflation"] for r in solved])
    add("chisq inflation range", fmt_range(ilo, ihi, 2), source="cv_series",
        note="ratio of achieved scatter to the formal error bar")
    add("comp stars median", fmt_int(one(
        cv, "SELECT n_comp FROM cv_series WHERE status='solved' "
            "ORDER BY n_comp LIMIT 1 OFFSET (SELECT count(*)/2 FROM "
            "cv_series WHERE status='solved')")), source="cv_series")
    # Check stars per series, from the series that HAVE target photometry.
    # euuma|e78|g is the exception the paper names in §5.2: five comparison
    # stars, no check stars and no target detection at all.
    chk = sorted({r[0] for r in cv.execute(
        "SELECT DISTINCT n_check FROM cv_series WHERE status='solved' "
        "AND n_target_points > 0")})
    add("check stars per series",
        fmt_int(chk[0]) if len(chk) == 1 else fmt_range(chk[0], chk[-1], 0),
        source="cv_series",
        note="SOLVE check stars: held out of every ensemble solve that "
             "produced target photometry; the one solved series with none "
             "produced no target points either. These are NOT the tie check "
             "stars of §3.2, which are catalogue stars held out of the tie "
             "fit and number 15--513 per block; the paper names the two "
             "populations differently for exactly that reason")
    add("series no check", fmt_int(one(
        cv, "SELECT count(*) FROM cv_series WHERE status='solved' "
            "AND n_check = 0")), source="cv_series",
        note="solved series holding out no check stars")

    tie = rows(cv, "SELECT * FROM cv_cattie WHERE is_primary=1")
    add("tie series", fmt_int(len(tie)), source="cv_cattie")
    add("tie tied", fmt_int(sum(1 for r in tie
                                if str(r["verdict"]).startswith("TIED"))),
        source="cv_cattie", note="series with a usable catalogue tie")
    add("tie at goal", fmt_int(sum(1 for r in tie
                                   if r["verdict"] == "TIED-GOAL")),
        source="cv_cattie",
        note="blocks meeting the 0.01--0.02 mag accuracy goal on the "
             "sigma-clipped check-star RMS, which is the statistic the tie "
             "stage grades on")
    add("tie at goal unclipped", fmt_int(sum(
        1 for r in tie if str(r["verdict"]).startswith("TIED")
        and r["check_rms"] is not None and r["check_rms"] <= 0.020)),
        source="cv_cattie",
        note="blocks meeting the same goal with every held-out check star "
             "kept. Fewer, which is why the paper gives both counts rather "
             "than saying the goal is missed 'either way' without checking")
    add("tie untied", fmt_int(sum(1 for r in tie
                                  if r["verdict"] == "UNTIED")),
        source="cv_cattie")
    # NAME THE UNTIED BLOCK.  "The one block with no usable tie" was used
    # in §3.2 for the primary block cv_cattie grades UNTIED (EU UMa's 2026
    # Fast series) and in §7 for the ST LMi block that has no cv_cattie row
    # at all -- two different series under one description, in a paper that
    # already had to separate two senses of "full-orbit night".  Both are
    # named from the database here so neither sentence can float free.
    def _blk_label(series_key: str) -> str:
        tgt, era, filt = str(series_key).split("|")
        return (f"{TARGET_LABEL.get(tgt, tgt)}'s "
                f"{ERA_LABEL.get(int(era[1:]), era)} ${filt}$ block")
    _unt = next((r for r in tie if r["verdict"] == "UNTIED"), None)
    add("tie untied block",
        _blk_label(_unt["series_key"]) if _unt else "\\NumMissing",
        source="cv_cattie",
        note="THE primary block graded UNTIED, named: this is the block §3.2 "
             "means by 'does not carry a usable tie'. It is not the ST LMi "
             "block of §7, which never reached the tie stage and has no "
             "cv_cattie row at all")
    # The solved block whose target IS detected but which produced no
    # catalogue-tied point, and so has no tie row: §7's block.
    _nozp = rows(cv, """
        SELECT s.series_key, s.n_target_rows FROM cv_series s
        LEFT JOIN cv_cattie c ON c.series_key = s.series_key
                             AND c.is_primary = 1
        WHERE s.status='solved' AND s.n_target_points = 0
          AND s.n_target_rows > 0 AND c.series_key IS NULL""")
    _nz = _nozp[0] if len(_nozp) == 1 else None
    add("untied block no tie stage",
        _blk_label(_nz["series_key"]) if _nz else "\\NumMissing",
        source="cv_series",
        note="the solved block in which the target IS detected but the "
             "ensemble produced no zero point, so it carries no "
             "natural-system magnitude and never entered the tie stage: "
             "§7's 'block with no usable tie at all'")
    add("untied block no tie stage rows",
        fmt_int(_nz["n_target_rows"]) if _nz else None,
        source="cv_series",
        note="instrumental target detections in that block, with real "
             "instrumental magnitudes and errors and no zero point: the "
             "reason it may not be described as a non-detection")
    # THE TIE ACCURACY, BOTH WAYS.  ``check_rms_clip`` is the check-star
    # RMS after sigma-clipping; ``check_rms`` is the same statistic with
    # every held-out star in it.  They differ by a factor 2.6 in the median
    # and by up to 7x per block, which is far too large a choice to make
    # silently in the number the paper calls its largest systematic.  Both
    # are emitted, and the paper quotes both.
    tied = [r for r in tie if str(r["verdict"]).startswith("TIED")]
    add("tie median accuracy mmag",
        fmt_float(_median([r["check_rms_clip"] for r in tied], 1000.0), 0),
        unit="mmag", source="cv_cattie",
        note="median achieved accuracy on the held-out check stars AFTER "
             "sigma-clipping, over the tied primary blocks, against a "
             "10--20 mmag goal; the unclipped median is quoted beside it "
             "and is 2.6x larger")
    add("tie median accuracy unclipped mmag",
        fmt_float(_median([r["check_rms"] for r in tied], 1000.0), 0),
        unit="mmag", source="cv_cattie",
        note="the SAME statistic with every held-out check star kept: no "
             "outlier rejection at all. Neither is obviously the right one "
             "-- a clipped star may be a blend or a variable rather than a "
             "tie failure -- so the paper prints both and the calibration "
             "goal is missed either way")
    add("tie blocks clipped", fmt_int(sum(
        1 for r in tied if (r["n_check_outlier"] or 0) > 0)),
        source="cv_cattie",
        note="tied blocks in which at least one held-out check star was "
             "clipped from the accuracy statistic")
    olo, ohi = _minmax([r["n_check_outlier"] for r in tied
                        if (r["n_check_outlier"] or 0) > 0])
    add("tie clipped stars range", fmt_range(olo, ohi, 0),
        source="cv_cattie", note="stars clipped per affected block")
    klo, khi = _minmax([r["n_check"] for r in tied])
    add("tie check stars range", fmt_range(klo, khi, 0),
        source="cv_cattie",
        note="TIE check stars per tied block: catalogue stars withheld from "
             "the tie FIT, an entirely different population from the four "
             "SOLVE check stars withheld from the ensemble solution. Table 2 "
             "prints statistics of both, in adjacent columns, and its "
             "caption says which is which")
    # TWO DIFFERENT "WORST" BLOCKS, AND THEY ARE NOT THE SAME BLOCK.
    # §3.2 offered ST LMi's 2024 z block as "the worst case" for the
    # clipped/unclipped choice.  It is the LARGEST UNCLIPPED VALUE (320
    # mmag), which is not what "worst case" reads as: the largest
    # SENSITIVITY to the choice is AN UMa's Mode0 g block, which moves by a
    # factor 13 on removal of a few stars against the z block's 7.5.  Both
    # are emitted, each labelled with what it is the maximum of, and each
    # carries its OWN check-star count -- §3.2 used to print the
    # whole-survey range "15--513" where that block's 15 belonged.
    worst = max(tied, key=lambda r: (r["check_rms"] or 0.0)) if tied else None

    def _ratio(r):
        return ((r["check_rms"] or 0.0) / r["check_rms_clip"]
                if r["check_rms_clip"] else 0.0)
    worst_ratio = max(tied, key=_ratio) if tied else None
    for blk, tag, what in (
            (worst, "tie worst",
             "the block with the LARGEST UNCLIPPED tie error"),
            (worst_ratio, "tie most sensitive",
             "the block MOST SENSITIVE to the clipping choice, i.e. the "
             "largest ratio of unclipped to clipped tie error")):
        _lab = (f"{TARGET_LABEL.get(str(blk['series_key']).split('|')[0], '?')}"
                f" {ERA_LABEL.get(blk['era_id'], '?')} "
                f"${str(blk['series_key']).split('|')[-1]}$"
                if blk else "n/a")
        add(f"{tag} block", _lab, source="cv_cattie",
            note=what + f" ({blk['series_key'] if blk else 'n/a'})")
        add(f"{tag} unclipped mmag",
            fmt_float(None if blk is None else 1000.0 * blk["check_rms"], 0),
            unit="mmag", source="cv_cattie", note=what + ", unclipped")
        add(f"{tag} clipped mmag",
            fmt_float(None if blk is None
                      else 1000.0 * blk["check_rms_clip"], 0),
            unit="mmag", source="cv_cattie",
            note="the same block after clipping")
        add(f"{tag} check stars",
            fmt_int(None if blk is None else blk["n_check"]),
            source="cv_cattie",
            note="THIS block's own held-out tie check-star count, not the "
                 "whole-survey range: a sentence about one block that "
                 "prints '15--513' has printed the wrong number")
        add(f"{tag} clipped stars",
            fmt_int(None if blk is None else blk["n_check_outlier"]),
            source="cv_cattie",
            note="stars clipped from this block's accuracy statistic")
        add(f"{tag} ratio",
            fmt_float(None if blk is None else _ratio(blk), 1),
            source="cv_cattie",
            note="unclipped tie error divided by the clipped one for this "
                 "block: how far the outlier choice moves it")
    # THE BARS FIGURE 6 ACTUALLY DRAWS.  §3.3 introduces the tie bar on the
    # headline colour figure and used to leave its size to the picture.
    # Each panel's bar is the two contributing blocks' check-star residuals
    # added in quadrature, so it is larger than either block alone and
    # larger than the medians §3.2 quotes -- the Mode0 g-r panel's
    # unclipped bar is 200 mmag against a 10--20 mmag goal.  A reader who
    # meets that bar in the figure without having met it in the text has
    # been surprised by the paper's largest systematic.
    _by_sk = {r["series_key"]: r for r in tie}
    _bars_clip, _bars_raw = [], []
    for _era, _ba, _bb in COLOUR_PANEL_PAIRS:
        _a = _by_sk.get(f"stlmi|e{_era}|{_ba}")
        _b = _by_sk.get(f"stlmi|e{_era}|{_bb}")
        if not (_a and _b):
            continue
        _bars_clip.append(math.hypot(_a["check_rms_clip"] or 0.0,
                                     _b["check_rms_clip"] or 0.0))
        _bars_raw.append(math.hypot(_a["check_rms"] or 0.0,
                                    _b["check_rms"] or 0.0))
    add("tie bar range clipped mmag",
        fmt_range(*_minmax(_bars_clip, 1000.0), 0), unit="mmag",
        source="cv_cattie",
        note="the sigma-clipped tie bar drawn on Figure 6's four "
             "colour--phase panels, low to high: each is the two "
             "contributing blocks' check-star residuals in quadrature")
    add("tie bar range unclipped mmag",
        fmt_range(*_minmax(_bars_raw, 1000.0), 0), unit="mmag",
        source="cv_cattie",
        note="the same four bars with every held-out check star kept. The "
             "upper end is the Mode0 g-r panel and is an order of magnitude "
             "above the calibration goal, which is why §3.3 states it")
    add("tie goal lo mmag", fmt_int(10), unit="mmag",
        source="ANALYSIS_STRATEGY §5", note="the stated accuracy goal")
    add("tie goal hi mmag", fmt_int(20), unit="mmag",
        source="ANALYSIS_STRATEGY §5")
    add("tie extrapolated", fmt_int(sum(
        1 for r in tie if r["colour_position"] == "extrapolated")),
        source="cv_cattie",
        note="blocks that place the CV OUTSIDE the fitted colour range, so "
             "their colour term is an extrapolation")
    add("tie colour unknown", fmt_int(sum(
        1 for r in tie if r["colour_position"] in (None, "unknown"))),
        source="cv_cattie",
        note="blocks in which the CV's own colour could not be measured at "
             "all, so it cannot even be checked against the fitted range")
    add("tie colour unsafe", fmt_int(sum(
        1 for r in tie if r["colour_position"] not in
        ("inside-span", "inside-core"))),
        source="cv_cattie",
        note="extrapolated plus unknown: every block whose colour term is "
             "not demonstrably an interpolation")
    add("tie colour inside", fmt_int(sum(
        1 for r in tie if r["colour_position"] in
        ("inside-span", "inside-core"))), source="cv_cattie")

    veto = rows(cv, "SELECT * FROM p2_cloud_series")
    fl, fh = _minmax([r["frac_vetoed"] for r in veto], 100.0)
    add("cloud veto frac hi", fmt_float(fh, 1), unit="per cent",
        source="p2_cloud_series", note="worst series")
    add("cloud vetoed frames", fmt_int(one(
        cv, "SELECT count(*) FROM p2_cloud_frame WHERE vetoed=1")),
        source="p2_cloud_frame")
    add("cloud checked frames", fmt_int(one(
        cv, "SELECT count(*) FROM p2_cloud_frame")),
        source="p2_cloud_frame")

    ext = rows(cv, "SELECT * FROM p2_extinction WHERE kpp IS NOT NULL")
    add("extinction groups", fmt_int(len(ext)), source="p2_extinction")
    add("extinction significant", fmt_int(sum(1 for r in ext
                                              if r["significant"])),
        source="p2_extinction",
        note="groups whose colour term survives an honest error bar")
    sig_terms = [r["term_p95_mmag"] for r in ext if r["significant"]]
    add("extinction max effect mmag", fmt_float(max(sig_terms)
                                                if sig_terms else None, 1),
        unit="mmag", source="p2_extinction",
        note="largest 95th-percentile shift the significant terms would "
             "apply; no term is applied anywhere in this paper")
    blo, bhi = _minmax([r["bound_mmag"] for r in ext])
    add("extinction bound range mmag", fmt_range(blo, bhi, 1), unit="mmag",
        source="p2_extinction", note="3-sigma bounds carried per group")

    lim = rows(cv, "SELECT * FROM p2_limit_series WHERE "
                   "median_limit_cal_mag IS NOT NULL")
    # A limit is only informative if it goes DEEPER than the star's own
    # typical brightness.  The comparison is therefore against the MEDIAN
    # detection in that series, computed here rather than assumed.
    shallow = 0
    for r in lim:
        mags = [x[0] for x in cv.execute(
            "SELECT cal_mag FROM cv_lightcurve WHERE series_key=? AND "
            "role='target' AND cal_mag IS NOT NULL ORDER BY cal_mag",
            (r["series_key"],))]
        if mags and r["median_limit_cal_mag"] < mags[len(mags) // 2]:
            shallow += 1
    add("limit series", fmt_int(len(lim)), source="p2_limit_series",
        note="series with enough validated forced photometry to publish a "
             "limit at all")
    add("limit shallow series", fmt_int(shallow), source="p2_limit_series",
        note="series whose MEDIAN limit is shallower than their MEDIAN "
             "detection: for these the non-detections say nothing about "
             "how faint the star went, so every faint-state fraction is an "
             "upper bound on a low-state duty cycle and not a measurement")
    add("limit forced points", fmt_int(one(
        cv, "SELECT sum(n_limits) FROM p2_limit_series")),
        source="p2_limit_series")

    # -- §4 Time-series analysis ----------------------------------------
    per = rows(cv, "SELECT * FROM p3_period WHERE status='ok'")
    add("period series", fmt_int(len(per)), source="p3_period")
    add("period tight", fmt_int(sum(1 for r in per
                                    if r["constraint_class"] == "TIGHT")),
        source="p3_period")
    add("period weak", fmt_int(sum(1 for r in per
                                   if r["constraint_class"] == "WEAK")),
        source="p3_period")
    add("period uninformative", fmt_int(sum(
        1 for r in per if r["constraint_class"] == "UNINFORMATIVE")),
        source="p3_period")
    alo, ahi = _minmax([r["alias_frac_max"] for r in per])
    add("alias frac max", fmt_float(ahi, 2), source="p3_period",
        note="largest fraction of the peak carried by a +/-1 c/d alias; "
             "no multi-night period is claimed without naming its family")
    add("alias frac range", fmt_range(alo, ahi, 2), source="p3_period")
    add("period prior families", fmt_int(sum(1 for r in per
                                             if r["family_code"] == "PRIOR")),
        source="p3_period",
        note="series whose period is a CONFIRMATION of the catalogue value "
             "inside a named alias family, not a determination")

    # TWO timing numbers, and the paper must state both.  The ANALYTIC
    # budget (ch_timing) is what the cadence and the depth allow in
    # principle.  The INJECTION test (p3_sigmat) is what recovering a
    # synthetic edge from the real timestamps and the real noise actually
    # achieved.  They differ by a factor of a few, in the direction that
    # matters, and quoting only the first would be the single most
    # misleading number this paper could carry.
    ct = rows(ch, """SELECT * FROM ch_timing WHERE target_key='stlmi'
                     AND night_kind='richest'""")
    ideal = [r["sigma_t_s"] for r in ct
             if r["regime"] == "per-cycle" and r["ingress_req"] == 0.01]
    mism = [r["sigma_t_s"] for r in ct
            if r["regime"] == "per-cycle shape-mismatched"]
    add("sigma t ideal s", fmt_float(min(ideal) if ideal else None, 0),
        unit="s", source="ch_timing",
        note="analytic per-cycle timing budget with the edge shape known "
             "exactly, richest ST LMi night")
    add("sigma t mismatched range s",
        fmt_range(min(mism), max(mism), 0) if mism else None, unit="s",
        source="ch_timing",
        note="the same budget with the edge shape wrong by a factor of "
             "five and the depth by 20 per cent")
    add("sigma t mismatched hi s",
        fmt_float(max(mism) if mism else None, 0), unit="s",
        source="ch_timing")

    sig = rows(cv, "SELECT * FROM p3_sigmat")
    add("sigma t injected best s",
        fmt_float(min((r["sigma_t_s"] for r in sig
                       if r["sigma_t_s"]), default=None), 0),
        unit="s", source="p3_sigmat",
        note="best single cell of the injection-recovery test at the real "
             "timestamps")
    tot = sorted(r["total_error_s"] for r in sig if r["total_error_s"])
    add("sigma t injected median s", fmt_float(_median(tot), 0), unit="s",
        source="p3_sigmat",
        note="MEDIAN total error the injection test actually achieved: the "
             "number a per-cycle timing claim must be judged against, and "
             "it exceeds the 60 s threshold")
    add("sigma t injected range s",
        fmt_range(min(tot) if tot else None, max(tot) if tot else None, 0),
        unit="s", source="p3_sigmat")
    thr = one(cv, "SELECT value FROM p3_meta WHERE key='timing_threshold_s'")
    add("sigma t threshold s", fmt_float(thr, 0), unit="s",
        source="p3_meta",
        note="the strategy's own per-epoch timing threshold, as the timing "
             "stage recorded it")
    add("sigma t trials", fmt_int(one(cv, "SELECT max(n_try) FROM "
                                      "p3_sigmat")), source="p3_sigmat")
    add("sigma t recovered frac", fmt_percent(one(
        cv, "SELECT max(recovered_fraction) FROM p3_sigmat"), 0),
        unit="per cent", source="p3_sigmat",
        note="best recovery fraction over the injection grid")

    add("edges accepted", fmt_int(one(cv, "SELECT sum(accepted) FROM "
                                      "p3_edge")), source="p3_edge")
    add("edges attempted", fmt_int(one(cv, "SELECT count(*) FROM p3_edge")),
        source="p3_edge")
    add("st lmi edges accepted", fmt_int(one(
        cv, "SELECT sum(accepted) FROM p3_edge WHERE target_key='stlmi'")),
        source="p3_edge")
    # THE PUBLISHED EPOCHS.  p3_oc holds the per-cycle edge residuals; the
    # injection test above does not license a single cycle's edge as an
    # epoch, so §4.2 forbids publishing one and CV-S9 collapses them into
    # p3_oc_night: one epoch per night per band.  Every O-C number the
    # paper states comes from that table, and the per-cycle count is
    # emitted beside it under its own name so a reader can see the
    # reduction rather than have to infer it.
    add("st lmi edges timed", fmt_int(one(
        cv, "SELECT count(*) FROM p3_oc WHERE target_key='stlmi'")),
        source="p3_oc",
        note="accepted per-cycle edges: the INPUTS to the O-C, not epochs; "
             "none is published on its own")
    add("st lmi oc epochs", fmt_int(one(
        cv, "SELECT count(*) FROM p3_oc_night WHERE target_key='stlmi'")),
        source="p3_oc_night",
        note="published timing epochs: one per night per band, the mean of "
             "that night's accepted per-cycle edges")
    add("st lmi oc nights", fmt_int(one(
        cv, "SELECT count(DISTINCT night) FROM p3_oc_night WHERE "
            "target_key='stlmi'")), source="p3_oc_night")
    ocn = rows(cv, "SELECT * FROM p3_oc_night WHERE target_key='stlmi'")
    olo, ohi = _minmax([r["oc_s"] for r in ocn])
    # ``to`` rather than an en dash: the low end is negative, and
    # "-239--175" is a range a reader has to decode.
    add("st lmi oc range s", fmt_range(olo, ohi, 0, dash=" to "), unit="s",
        source="p3_oc_night")
    ncy_lo, ncy_hi = _minmax([r["n_cycles"] for r in ocn])
    add("st lmi oc cycles per epoch", fmt_range(ncy_lo, ncy_hi, 0),
        source="p3_oc_night",
        note="per-cycle edges averaged into each published epoch")
    # HOW OFTEN THE MEAN HAS ONE TERM IN IT.  The rule §4.2 states is about
    # the ERROR BAR, not the estimator: no per-cycle edge is published
    # carrying its own per-cycle sigma_t.  On a night with one accepted
    # edge the epoch IS that edge, with the injection budget attached
    # instead.  Stating the rule without stating this count -- an earlier
    # revision did -- makes a true rule read as a false one.
    add("st lmi single cycle epochs",
        fmt_int(sum(1 for r in ocn if r["n_cycles"] == 1)),
        source="p3_oc_night",
        note="published epochs built from exactly one timed cycle: the "
             "mean of one cycle is that cycle, republished with the "
             "injection-demonstrated budget rather than the edge fit's own "
             "error bar")
    add("st lmi epochs below threshold",
        fmt_int(sum(1 for r in ocn if not r["meets_threshold"])),
        source="p3_oc_night",
        note="published epochs whose error bar exceeds the 60 s threshold. "
             "The threshold is a per-epoch QUALITY LABEL, not a publication "
             "gate: an epoch outside it is published with an error bar that "
             "says so")
    cc = {r["target_key"]: r for r in rows(cv, "SELECT * FROM p3_cycle_count")}
    st_cc = cc.get("stlmi", {})
    add("st lmi oc rms s", fmt_float(st_cc.get("oc_night_rms_s"), 0),
        unit="s", source="p3_cycle_count",
        note="RMS of the published per-night epochs about the FITTED "
             "constant offset, not about the catalogue epoch's phase zero; "
             "see the offset macros below")
    # =====================================================================
    # THE CONSTANT THAT WAS SUBTRACTED, AND WHAT THAT MAKES THE O-C
    # ---------------------------------------------------------------------
    # CV-S9 subtracts the mean of the per-cycle O-C before writing p3_oc,
    # because the bright-phase edge is not at phase zero of the VSX
    # ephemeris and the raw residual therefore carries a constant that is a
    # property of the FEATURE.  The subtraction is standard and it is
    # right, but until this revision no word of it reached the manuscript,
    # which said in three places that the residuals were measured "against
    # the catalogue ephemeris" -- and a reader took that to mean the edge
    # falls at the catalogue epoch's phase.  It does not: it falls 0.157
    # cycles later, and that offset was fitted out of the same edges.
    #
    # Two consequences, and both are emitted here rather than left for a
    # reader to find in the products.  (1) The offset is itself a
    # measurement and is published as one.  (2) One parameter was estimated
    # from these edges, so the reduced chi-squared has nu = N - 1 degrees
    # of freedom and not N; the stored oc_night_chi2nu is chi-squared PER
    # EPOCH and is not the same number.
    # =====================================================================
    _Ps = (float(st_cc["period_d"]) * 86400.0
           if st_cc.get("period_d") else None)
    add("st lmi oc offset s", fmt_int(st_cc.get("oc_mean_s")), unit="s",
        source="p3_cycle_count",
        note="THE CONSTANT REMOVED FROM THE O-C. The mean of the per-cycle "
             "residuals against the catalogue ephemeris, subtracted before "
             "anything is plotted or fitted: the phase of the timed "
             "bright-phase edge relative to the VSX epoch's phase zero. It "
             "is a property of the feature and of the catalogue's choice of "
             "fiducial, not of the clock, and the O-C therefore tests "
             "whether this interval is CONSTANT, never whether it is zero")
    add("st lmi oc offset cycles",
        fmt_float(None if not (_Ps and st_cc.get("oc_mean_s"))
                  else st_cc["oc_mean_s"] / _Ps, 3),
        unit="cycles", source="p3_cycle_count",
        note="the same offset as a fraction of an orbit; it agrees with the "
             "phase of the steepest faintward gradient of the folded "
             "profile, measured independently, which is the check that it "
             "is the feature's own phase and not a clock error")
    _chi2_sum = sum((r["oc_s"] / r["oc_sigma_s"]) ** 2 for r in ocn
                    if r["oc_sigma_s"])
    _dof = len(ocn) - 1                # the one absorbed constant, above
    add("st lmi oc dof", fmt_int(_dof), source="p3_oc_night",
        note="degrees of freedom of the O-C: the published epoch count less "
             "the one constant absorbed from the same edges")
    add("st lmi oc chisq", fmt_float(_chi2_sum / _dof if _dof else None, 2),
        source="p3_oc_night",
        note="REDUCED chi-squared of the published O-C about zero on "
             f"{_dof} degrees of freedom, against the injection-licensed "
             "per-night errors. The denominator is the epoch count less "
             "one, because the constant offset above was estimated from "
             "these same edges; dividing by the epoch count instead -- "
             "which is what p3_cycle_count.oc_night_chi2nu stores -- "
             "understates it")
    add("st lmi sigma night median s",
        fmt_float(st_cc.get("sigma_night_median_s"), 0), unit="s",
        source="p3_cycle_count",
        note="median error bar on a published epoch: the injection random "
             "term divided by the square root of the cycles averaged, added "
             "in quadrature to the injection bias, which does not average "
             "down")
    add("st lmi sigma night range s",
        fmt_range(st_cc.get("sigma_night_lo_s"),
                  st_cc.get("sigma_night_hi_s"), 0), unit="s",
        source="p3_cycle_count")
    add("st lmi epochs at threshold",
        fmt_int(st_cc.get("n_night_at_threshold")), source="p3_cycle_count",
        note="published epochs whose error bar is inside the 60 s "
             "threshold; the rest are nights on which too few cycles were "
             "timed to average down to it")
    # TWO CYCLE COUNTS, AND THEY ARE NOT INTERCHANGEABLE.
    #   \NumStLmiCycles      -- catalogue epoch to our last timed edge.
    #                           The EXTRAPOLATION BASELINE: it is what the
    #                           sensitivity to a period ERROR scales with,
    #                           and it is what §4.3's cycle-count argument
    #                           is about.
    #   \NumStLmiEpochSpan*  -- first timed edge to last timed edge.  The
    #                           OBSERVED SPAN: it is what the leverage on a
    #                           period CHANGE scales with, and it is what
    #                           any sentence about "the epochs" is about.
    # They differ by a factor 2.5 here (21,869 against 8,688 cycles, 4.7
    # against 1.9 yr), and an earlier revision printed the first where the
    # second belonged and then glossed it as "some two years".
    add("st lmi cycles", fmt_int(st_cc.get("n_cycles_last")),
        source="p3_cycle_count",
        note="cycles between the CATALOGUE EPOCH and the last timed edge: "
             "the extrapolation baseline, not the span of the epochs")
    _P = st_cc.get("period_d")
    add("st lmi cycles years",
        fmt_float(None if not (_P and st_cc.get("n_cycles_last")) else
                  st_cc["n_cycles_last"] * float(_P) / 365.25, 1),
        unit="yr", source="p3_cycle_count",
        note="the same baseline in years")
    add("st lmi epoch span cycles", fmt_int(st_cc.get("n_cycles_span")),
        source="p3_cycle_count",
        note="cycles between the FIRST and LAST timed edge: the span the "
             "published epochs actually cover")
    add("st lmi epoch span days", fmt_int(st_cc.get("span_d")), unit="d",
        source="p3_cycle_count")
    add("st lmi epoch span years",
        fmt_float(None if st_cc.get("span_d") is None
                  else st_cc["span_d"] / 365.25, 1), unit="yr",
        source="p3_cycle_count")
    add("st lmi drift cycles", fmt_float(st_cc.get("drift_cycles"), 4),
        unit="cycles", source="p3_cycle_count",
        note="accumulated phase drift at the catalogue period's quoted "
             "precision; the integer cycle count survives it")
    add("st lmi fitted period d",
        fmt_float(st_cc.get("fitted_period_night_d"), 8), unit="d",
        source="p3_cycle_count",
        note="period refitted on the published per-night epochs; a "
             "consistency check, never a replacement for the catalogue "
             "value")
    add("st lmi fitted period sigma d",
        fmt_sci(st_cc.get("fitted_period_night_sigma_d"), 1), unit="d",
        source="p3_cycle_count")
    _fp, _fs, _cp = (st_cc.get("fitted_period_night_d"),
                     st_cc.get("fitted_period_night_sigma_d"),
                     st_cc.get("period_d"))
    add("st lmi period agreement sigma",
        fmt_float(abs(_fp - _cp) / _fs if _fp and _fs and _cp else None, 1),
        source="p3_cycle_count",
        note="separation between our refitted period and the catalogue "
             "value, in units of our own error bar")
    # =====================================================================
    # WHERE EVERY PUBLISHED ERROR BAR COMES FROM, AND THE CHECK ON IT
    # ---------------------------------------------------------------------
    # p3_sigmat -- the injection grid the whole timing budget rests on --
    # exists for ONE night of ONE target in ONE readout mode.  night_epochs
    # serves that budget out by band NAME, so a 2024 High Gain epoch in the
    # G slot silently receives the 2025 Mode0 g number, and a band the grid
    # never saw (z) falls through to the whole-grid median.  That transfer
    # crosses the very instrument seam §2.1 insists cannot be crossed for
    # saturation, so it has to be stated, and it has to be checked.  The
    # check is oc_night_chi2nu_edge: the same reduced chi-squared computed
    # from the edge fits' OWN Monte-Carlo errors, which are measured on
    # every night in its own era.
    # =====================================================================
    add("timing budget night", one(
        cv, "SELECT value FROM p3_meta WHERE key='sigmat_night'"),
        source="p3_meta",
        note="the single night the injection grid was run on")
    _bs = rows(cv, "SELECT DISTINCT series_key FROM p3_sigmat")
    add("timing budget series", fmt_int(len(_bs)), source="p3_sigmat",
        note="series the injection grid covers: "
             + ", ".join(sorted(r["series_key"] for r in _bs)))
    _rand = [r["sigma_random_s"] * math.sqrt(r["n_cycles"]) for r in ocn
             if r["sigma_random_s"] and r["n_cycles"]]
    add("sigma random range s", fmt_range(*_minmax(_rand), 0), unit="s",
        source="p3_oc_night",
        note="the injection test's random term per band before averaging; "
             "it falls as the square root of the cycles a night timed")
    add("sigma bias floor s",
        fmt_float(_median([r["sigma_floor_s"] for r in ocn]), 1), unit="s",
        source="p3_oc_night",
        note="the injection test's bias, which does NOT average down and "
             "is the irreducible floor under every published epoch")
    _budget_era = 76          # the era of the night p3_sigmat was run on
    add("st lmi epochs transferred budget",
        fmt_int(sum(1 for r in ocn if r["era_id"] != _budget_era)),
        source="p3_oc_night",
        note="published epochs whose error bar was measured in a DIFFERENT "
             "instrument era from the one they were observed in: the "
             "injection grid exists only for the 2025 Mode0 camera, and "
             "these carry it by band slot regardless")
    add("st lmi epochs fallback budget",
        fmt_int(sum(1 for r in ocn if str(r["budget_band"]) == "*")),
        source="p3_oc_night",
        note="published epochs in a band the injection grid never saw at "
             "all, which take the whole-grid median instead")
    _chi2_edge = (sum((r["oc_s"] / r["oc_sigma_edge_s"]) ** 2 for r in ocn)
                  if all(r["oc_sigma_edge_s"] for r in ocn) else None)
    add("st lmi oc chisq edge",
        fmt_float(None if not (_chi2_edge and _dof)
                  else _chi2_edge / _dof, 2),
        source="p3_oc_night",
        note="THE CHECK ON THE TRANSFER: the same reduced chi-squared, on "
             f"the same {_dof} degrees of freedom, computed from the edge "
             "fits' own Monte-Carlo errors instead of the transported "
             "injection budget. It agrees, so the null does not rest on "
             "the transfer")
    add("st lmi sigma night edge median s",
        fmt_float(st_cc.get("sigma_night_edge_median_s"), 0), unit="s",
        source="p3_cycle_count",
        note="median epoch error under the edge fits' own errors")
    # The same check split by era, because the transfer is exactly an
    # across-era assumption and a reader is entitled to see it tested
    # inside each era separately.
    for era, tag in ((7, "high gain"), (76, "mode zero")):
        sub = [r for r in ocn if r["era_id"] == era]
        add(f"st lmi epochs {tag}", fmt_int(len(sub)),
            source="p3_oc_night")
        # Chi-squared PER EPOCH here, not per degree of freedom: the one
        # absorbed constant was fitted across the whole edge set and
        # belongs to neither era, so charging it to one of them (or half
        # to each) would be arbitrary.  The prose names the denominator.
        add(f"st lmi oc chisq {tag}",
            fmt_float(sum((r["oc_s"] / r["oc_sigma_s"]) ** 2
                          for r in sub) / len(sub) if sub else None, 2),
            source="p3_oc_night",
            note=f"chi-squared PER EPOCH about zero over the {len(sub)} "
                 f"epochs of this instrument era alone. The denominator is "
                 f"the epoch count, not a degree-of-freedom count: the "
                 f"absorbed constant is global to the target and is not "
                 f"attributable to either era")
    # =====================================================================
    # THE PERIOD DERIVATIVE THE NULL BOUNDS
    # ---------------------------------------------------------------------
    # "We report no period change" is an absence until it carries a
    # number.  A weighted quadratic through the published epochs gives one.
    # =====================================================================
    add("st lmi pdot sigma", fmt_float(
        abs(st_cc["quad_coeff_s_per_cycle2"] /
            st_cc["quad_sigma_s_per_cycle2"])
        if st_cc.get("quad_sigma_s_per_cycle2") else None, 1),
        source="p3_cycle_count",
        note="significance of the quadratic term in the O-C: below 2, so "
             "no period change is detected")
    add("st lmi pdot limit", fmt_sci(st_cc.get("pdot_limit3"), 1),
        source="p3_cycle_count",
        note="3-sigma upper bound on |dP/dt| (dimensionless) from a "
             "weighted quadratic through the published epochs: what the "
             "null actually excludes")
    # WHICH 3-SIGMA CONVENTION.  ``pdot_limit3`` is |Pdot_fit| + 3 sigma,
    # not 3 sigma: the two differ by the fitted value, and a reader given
    # only the phrase "3-sigma bound" cannot tell which was meant. Both are
    # emitted so §5.1 can name the convention it uses.
    add("st lmi pdot three sigma",
        fmt_sci(None if st_cc.get("pdot_sigma") is None
                else 3.0 * st_cc["pdot_sigma"], 1),
        source="p3_cycle_count",
        note="three times the standard error on the fitted Pdot alone, "
             "WITHOUT the fitted value added: the other convention a reader "
             "might assume behind the phrase '3-sigma bound'")
    add("st lmi pdot limit s per yr", fmt_float(
        None if st_cc.get("pdot_limit3") is None
        else st_cc["pdot_limit3"] * 365.25 * 86400.0, 2),
        unit="s/yr", source="p3_cycle_count",
        note="the same bound as a change in the orbital period per year")
    # P / |dP/dt|, with the PERIOD in the numerator.  It was written as
    # 1 / (pdot_limit3 * 365.25), which drops P entirely — and since
    # pdot_limit3 is dimensionless (d/d), that expression is not a time at
    # all.  It published 7.6e5 yr where the paper's own P and Pdot give
    # 6.0e4 yr, a factor 12.7, in the abstract and a conclusion.
    add("st lmi pdot timescale yr", fmt_sci(
        pdot_timescale_yr(st_cc.get("period_d"), st_cc.get("pdot_limit3")),
        1), unit="yr",
        source="p3_cycle_count",
        note="P / |dP/dt| at the 3-sigma bound: the shortest period-change "
             "timescale these epochs are consistent with")
    # =====================================================================
    # THE INTER-BAND OFFSET: A NON-DETECTION, PUBLISHED AS ONE
    # ---------------------------------------------------------------------
    # An earlier revision put "the epoch of a bright-phase edge is band
    # dependent" in the abstract and the conclusions.  Every row of
    # p3_band_pair carries significant=0 and no pooled pair reaches 2
    # sigma.  These macros make the null quotable: the number of pairs, how
    # many are significant, the strongest one with its error bar and its
    # significance, and the tightest bound the pooled pairs place on any
    # such offset.  A claim that cannot be stated with these numbers beside
    # it is a claim this data set does not support.
    # =====================================================================
    bp = rows(cv, "SELECT * FROM p3_band_pair WHERE target_key='stlmi'")
    pooled = [r for r in bp if str(r["night"]).lower().startswith("(pooled")
              and r["sigma_s"]]
    add("band pairs", fmt_int(len(bp)), source="p3_band_pair",
        note="paired same-cycle edge-time differences: one row per night "
             "per band pair, plus one pooled row per band pair per era")
    add("band pairs pooled", fmt_int(len(pooled)), source="p3_band_pair",
        note="the publishable inter-band numbers; the per-night rows are "
             "their components, not independent results")
    add("band pairs significant",
        fmt_int(sum(1 for r in bp if r["significant"])),
        source="p3_band_pair",
        note="ZERO. No band pair, pooled or per night, reaches the "
             "3-sigma bar, and no pooled pair reaches even 2 sigma")
    add("band pair sigma bar", fmt_int(3), source="ANALYSIS_STRATEGY §4",
        note="the significance a band-to-band offset must reach to be "
             "called a detection")
    if pooled:
        top = max(pooled, key=lambda r: abs(r["delta_s"]) / r["sigma_s"])
        add("band offset top pair",
            f"${top['band_a']}-{top['band_b']}$", source="p3_band_pair",
            note="the MOST SIGNIFICANT pooled band pair, over "
                 f"{top['n_cycles']} paired cycles")
        add("band offset top s", fmt_float(top["delta_s"], 0), unit="s",
            source="p3_band_pair")
        add("band offset top err s", fmt_float(top["sigma_s"], 0), unit="s",
            source="p3_band_pair")
        add("band offset top sigma",
            fmt_float(abs(top["delta_s"]) / top["sigma_s"], 1),
            source="p3_band_pair",
            note="its significance: the largest any band pair here reaches")
        add("band offset top cycles", fmt_int(top["n_cycles"]),
            source="p3_band_pair")
        big = max(pooled, key=lambda r: abs(r["delta_s"]))
        add("band offset abs max s", fmt_float(abs(big["delta_s"]), 0),
            unit="s", source="p3_band_pair",
            note=f"largest pooled offset in absolute size "
                 f"(${big['band_a']}-{big['band_b']}$), at "
                 f"{abs(big['delta_s']) / big['sigma_s']:.1f} sigma")
        # THE BOUND: |delta| + 2 sigma, per pooled pair.  TWO of them are
        # emitted, and which one a sentence may use depends on its
        # quantifier.  "The g-r pair bounds an offset below X" is a
        # statement about g-r and takes the TIGHTEST.  "Any band-to-band
        # offset is below X" quantifies over the pairs and must take the
        # WEAKEST, or it excludes offsets these data do not exclude: the
        # High Gain G-I pair allows 313 s, 2.3x the g-r figure, and a
        # 250 s offset in G-I is entirely consistent with this data set.
        # Figure 9(b) prints all five honestly; the text used to print one.
        def _bound(r):
            return abs(r["delta_s"]) + 2.0 * r["sigma_s"]
        bnd = min(pooled, key=_bound)
        weak = max(pooled, key=_bound)
        _Pss = float(_P) * 86400.0 if _P else None
        for r_, tag_, note_ in (
                (bnd, "band offset bound",
                 "TIGHTEST 2-sigma upper bound over the pooled pairs. Use "
                 "it only where the sentence names this pair; a sentence "
                 "quantifying over ALL pairs must use the weakest bound "
                 "below"),
                (weak, "band offset weakest bound",
                 "WEAKEST 2-sigma upper bound over the pooled pairs, and "
                 "therefore the only one that may stand beside the word "
                 "'any': an offset smaller than this is consistent with "
                 "every pooled pair, which is what a limit on a "
                 "band-to-band offset has to mean")):
            b_ = _bound(r_)
            add(f"{tag_} pair", f"${r_['band_a']}-{r_['band_b']}$",
                source="p3_band_pair",
                note=f"pooled over {int(r_['n_cycles'])} paired cycles in "
                     f"{ERA_LABEL.get(r_['era_id'], 'era ' + str(r_['era_id']))}")
            add(f"{tag_} s", fmt_float(b_, 0), unit="s",
                source="p3_band_pair", note=note_)
            add(f"{tag_} cycles",
                fmt_float(b_ / _Pss if _Pss else None, 3), unit="cycles",
                source="p3_band_pair",
                note="the same bound as a fraction of an orbit")
        add("band offset bound era", ERA_LABEL.get(bnd["era_id"], "?"),
            source="p3_band_pair",
            note="the instrument era the tightest pooled pair belongs to")
        add("band offset weakest bound era",
            ERA_LABEL.get(weak["era_id"], "?"), source="p3_band_pair",
            note="the instrument era the weakest pooled pair belongs to")
        add("band offset pooled range s",
            fmt_range(min(r["delta_s"] for r in pooled),
                      max(r["delta_s"] for r in pooled), 0, dash=" to "),
            unit="s", source="p3_band_pair",
            note="every pooled offset; all are consistent with zero")
    add("st lmi phase spread", fmt_float(st_cc.get("phase_spread"), 3),
        unit="cycles", source="p3_cycle_count",
        note="circular scatter of the accepted edges; below the bar that "
             "says they time ONE feature")
    add("one feature bar", fmt_float(one(
        cv, "SELECT value FROM p3_meta WHERE "
            "key='one_feature_bar_cycles'"), 2), unit="cycles",
        source="p3_meta",
        note="phase scatter above which pooled edges are timing different "
             "features and no O-C may be built")
    for tgt, tag in (("anuma", "an uma"), ("vvpup", "vv pup")):
        add(f"{tag} phase spread", fmt_float(
            cc.get(tgt, {}).get("phase_spread"), 3), unit="cycles",
            source="p3_cycle_count",
            note="above the 0.05 bar: these epochs time different features "
                 "and CV-S9 refused to build an O-C from them")

    st = rows(cv, "SELECT * FROM p3_state_series WHERE separability "
                  "IS NOT NULL")
    add("state series graded", fmt_int(len(st)), source="p3_state_series")
    add("state series bimodal", fmt_int(sum(1 for r in st if r["bimodal"])),
        source="p3_state_series",
        note="series in which two accretion states separate")
    add("state separability bar", fmt_float(one(
        cv, "SELECT value FROM p3_meta WHERE key='state_bimodal_bar'"), 2),
        source="p3_meta", note="Otsu separability threshold")
    slo, shi = _minmax([r["separability"] for r in st])
    add("state separability range", fmt_range(slo, shi, 2),
        source="p3_state_series")
    # Figure 8 draws every class p3_state_night contains, so §4.4 states
    # how many nights fall in each rather than leaving three of the five
    # colours in the panel unexplained.
    for state in ("HIGH", "LOW", "INTERMEDIATE", "UNCLASSIFIED", "UNKNOWN"):
        add(f"state nights {state.lower()}", fmt_int(one(
            cv, "SELECT count(*) FROM p3_state_night WHERE state=?",
            (state,))), source="p3_state_night")
    add("state nights total", fmt_int(one(
        cv, "SELECT count(*) FROM p3_state_night")),
        source="p3_state_night",
        note="classified nights across all series and targets")

    det = rows(cv, "SELECT * FROM p3_detrend")
    dlo, dhi = _minmax([r["frac_detrend"] for r in det], 100.0)
    add("detrend loss range", fmt_range(dlo, dhi, 0), unit="per cent",
        source="p3_detrend",
        note="signal a naive detrending removes at these periods, measured "
             "by injecting a signal and filtering it; the joint fit is used "
             "instead everywhere in this paper. This is NOT the same "
             "measurement as the recovery contours below and the two point "
             "in opposite directions")
    # THE DIRECTION OF FIGURE 12(a)'s GAP, MEASURED RATHER THAN ASSUMED.
    # §4.6 and the caption both called the raw-versus-detrended gap "the
    # sensitivity cost of detrending".  Pairing ch_contour's 'season' and
    # 'season-dt' scopes period bin by period bin says otherwise: the
    # detrended contour is LOWER -- a smaller amplitude recovered 90 per
    # cent of the time, i.e. better sensitivity -- for four of the five
    # series, because the filter removes red noise along with signal. The
    # panel shows a trade, not a cost, and the text may not tell a reader
    # to read it as one.
    _dt_ratios: list[float] = []
    _dt_worse = 0
    for _s in sorted({str(r["scope"]).rsplit("|", 1)[0] for r in rows(
            ch, "SELECT scope FROM ch_contour WHERE scope LIKE '%|season%'")}):
        _raw = {r["period_d"]: r["amp90"] for r in rows(
            ch, "SELECT period_d, amp90 FROM ch_contour WHERE scope=? "
                "AND amp90 IS NOT NULL", (f"{_s}|season",))}
        _dtc = {r["period_d"]: r["amp90"] for r in rows(
            ch, "SELECT period_d, amp90 FROM ch_contour WHERE scope=? "
                "AND amp90 IS NOT NULL", (f"{_s}|season-dt",))}
        _shared = sorted(set(_raw) & set(_dtc))
        if not _shared:
            continue
        _med = _median([_dtc[p] / _raw[p] for p in _shared])
        _dt_ratios.append(_med)
        _dt_worse += int(_med > 1.0)
    add("detrend contour ratio range", fmt_range(*_minmax(_dt_ratios), 2),
        source="ch_contour",
        note="median over the shared injected periods of the DETRENDED 90 "
             "per cent recovery amplitude divided by the RAW one, per "
             "series. Below one means detrending RECOVERS a smaller "
             "amplitude, i.e. does better, not worse")
    add("detrend contour series worse", fmt_int(_dt_worse),
        source="ch_contour",
        note="series whose detrended contour is worse than their raw one")
    add("detrend contour series better", fmt_int(len(_dt_ratios) - _dt_worse),
        source="ch_contour",
        note="series whose detrended contour is BETTER than their raw one, "
             "which is the majority and the opposite of what §4.6 and "
             "Figure 12's caption used to assert")
    add("detrend contour series", fmt_int(len(_dt_ratios)),
        source="ch_contour",
        note="series carrying both a raw and a detrended season contour")

    # -- §5 Results: ST LMi's folded morphology -------------------------
    # The one QUALITATIVE claim §5.1 makes about the fold, made checkable.
    # An earlier draft said "a sharp rise and a slower decline", which is
    # the wrong way round and contradicted §4.3, Figure 5 and the edge
    # fitter, all three of which treat the FALLING edge as the sharp
    # feature.  These two macros are the gradients that settle it.
    fall_ph, ratios = [], []
    for sk, in cv.execute("SELECT series_key FROM cv_series WHERE "
                          "target_key='stlmi' AND status='solved' AND "
                          "lower(filter) IN ('g','r','i') "
                          "ORDER BY series_key"):
        prof = _folded_profile(cv, sk, "stlmi")
        if prof is None:
            continue
        centres, med = prof
        grad = [med[(i + 1) % len(med)] - med[i - 1] for i in range(len(med))]
        gv = [(g, c) for g, c in zip(grad, centres)
              if g is not None and math.isfinite(g)]
        if len(gv) < 8:
            continue
        fade = max(gv)                    # steepest FAINTWARD gradient
        rise = min(gv)                    # steepest BRIGHTWARD gradient
        fall_ph.append(fade[1])
        if rise[0] != 0:
            ratios.append(abs(fade[0] / rise[0]))
    flo, fhi = _minmax(fall_ph)
    add("st lmi fall phase range", fmt_range(flo, fhi, 2), unit="cycles",
        source="cv_lightcurve",
        note="phase of the steepest FAINTWARD gradient of the 40-bin "
             "folded profile, over every solved g/r/i ST LMi series in both "
             "eras: the falling edge §4.3 times")
    rlo2, rhi2 = _minmax(ratios)
    add("st lmi fall rise ratio range", fmt_range(rlo2, rhi2, 1),
        source="cv_lightcurve",
        note="steepest faintward gradient divided by the steepest "
             "brightward one, per series: the decline is the sharp feature "
             "and the rise the gradual one, in every band and both eras")

    # -- §5 Results: EU UMa's 2026 Fast-mode series ---------------------
    eu = rows(cv, "SELECT * FROM cv_series WHERE target_key='euuma' "
                  "AND era_id IN (78, 79)")
    add("eu uma fast frames", fmt_int(one(
        cv, "SELECT count(*) FROM cv_frames WHERE target_key='euuma' "
            "AND era_id IN (78,79)")), source="cv_frames",
        note="frames of the merged 2026 Fast-mode series: counted in the "
             "archive census, used in no measurement")
    add("eu uma fast comp stars",
        fmt_int(max((r["n_comp"] or 0) for r in eu) if eu else None),
        source="cv_series")
    add("eu uma fast check stars",
        fmt_int(max((r["n_check"] or 0) for r in eu) if eu else None),
        source="cv_series", note="ZERO: nothing held out, so there is no "
                                 "measurement of this series' accuracy")
    add("eu uma fast points",
        fmt_int(sum((r["n_target_points"] or 0) for r in eu)),
        source="cv_series",
        note="ZERO catalogue-tied target measurements: the star is not "
             "detected in this series at all, which is why excluding it "
             "removes no photometry from the paper")

    # -- §5 Results: YZ Cnc -------------------------------------------
    gate = rows(cv, "SELECT * FROM p4_gate")
    add("gate lines", fmt_int(len(gate)), source="p4_gate")
    add("gate passes", fmt_int(sum(1 for r in gate if r["passes"])),
        source="p4_gate", note="lines of the strategy's 4.19 S/N gate that "
                               "the photometry clears")
    # The two eras' MEASURED floors, per era, over the same statistic and
    # the same selection.  An earlier revision typed 27--46 for High Gain,
    # which was the floor at the shortest sampled lag only, and 5--12 for
    # the Sloan era, which was the whole range -- two different bases
    # printed side by side as though they were one comparison.
    def _floor_range(prefix):
        f = [r["sf_floor"] for r in rows(
            cv, "SELECT sf_floor FROM p4_flicker WHERE sf_floor IS NOT NULL "
                "AND upper(state)='QUIESCENT' AND series_key LIKE ?",
            (prefix,))]
        return _minmax(f, 1000.0)
    add("floor hg range mmag", fmt_range(*_floor_range("yzcnc|e7|%"), 0),
        unit="mmag", source="p4_flicker",
        note="measured noise floor from magnitude-matched field stars "
             "through the same frames, High Gain era, over every sampled "
             "timescale of the quiescent runs")
    add("floor sloan range mmag",
        fmt_range(*_floor_range("yzcnc|e72|%"), 0), unit="mmag",
        source="p4_flicker", note="same statistic, same selection, "
                                  "Sloan-era frames")
    for tag, prefix in (("hg", "yzcnc|e7|%"), ("sloan", "yzcnc|e72|%")):
        add(f"floor {tag} exp s", fmt_float(_median(
            [r["exptime"] for r in rows(
                cv, "SELECT exptime FROM cv_frames WHERE series_key LIKE ? "
                    "AND night IN (SELECT DISTINCT night FROM p4_flicker "
                    "WHERE upper(state)='QUIESCENT')", (prefix,))]), 0),
            unit="s", source="cv_frames",
            note="median exposure of the frames the floor was measured in")

    # TESTABLE QUIESCENT scopes only: those in a quiescent run for which
    # the red-noise contour could be measured.  Pooling the outburst runs
    # in -- an earlier version of this collector did -- reports an
    # outburst's slope as an orbital hump, and quotes amplitudes the
    # fallback never claimed.
    run = rows(cv, "SELECT * FROM p4_run WHERE hump_amp IS NOT NULL "
                   "AND amp90_self IS NOT NULL "
                   "AND upper(state)='QUIESCENT'")
    hlo, hhi = _minmax([r["hump_amp"] for r in run], 1000.0)
    add("hump amp range mmag", fmt_range(hlo, hhi, 0), unit="mmag",
        source="p4_run", note="fitted semi-amplitude on the quiescent runs")
    flo, fhi = _minmax([r["amp90_field"] for r in run], 1000.0)
    add("hump field contour mmag", fmt_range(flo, fhi, 0), unit="mmag",
        source="p4_run", note="90 per cent recovery against field stars")
    zlo, zhi = _minmax([r["amp90_self"] for r in run], 1000.0)
    add("hump self contour mmag", fmt_range(zlo, zhi, 0), unit="mmag",
        source="p4_run",
        note="90 per cent recovery against the star's own night-rolled "
             "residuals: the red-noise null the hump does not clear")
    add("hump scopes tested", fmt_int(len(run)), source="p4_run",
        note="quiescent scopes on which BOTH contours could be measured")
    add("hump scopes quiescent", fmt_int(one(
        cv, "SELECT count(*) FROM p4_run WHERE upper(state)='QUIESCENT'")),
        source="p4_run",
        note="all quiescent scopes, including those carrying only the "
             "instrumental contour")
    add("hump detections", fmt_int(sum(
        1 for r in run if str(r["detection"]).upper().startswith("DETECT"))),
        source="p4_run",
        note="ZERO: the hump clears the red-noise contour on no scope at "
             "all, and the instrumental contour on all but one")
    add("hump above field contour", fmt_int(sum(
        1 for r in run if r["hump_amp"] is not None
        and r["amp90_field"] is not None
        and r["hump_amp"] > r["amp90_field"])), source="p4_run",
        note="scopes on which the fitted hump exceeds the INSTRUMENTAL "
             "contour, i.e. on which the photometry could have seen a hump "
             "of the fitted size at all; on the remaining scope it could "
             "not, so that scope is uninformative rather than a "
             "non-detection")
    # Within-night filter-to-filter phase agreement and the night-to-night
    # shift, as CIRCULAR statistics on p4_run.hump_phase.  A phase is an
    # angle: a linear standard deviation of 0.98 and 0.09 would report a
    # disagreement of half a cycle where there is none.
    def _circ(ph):
        if not ph:
            return None, None
        z = sum(cmath.exp(2j * math.pi * float(p)) for p in ph) / len(ph)
        r = abs(z)
        sd = (math.sqrt(-2.0 * math.log(r)) / (2 * math.pi)
              if r > 1e-12 else 0.5)
        return (math.degrees(cmath.phase(z)) / 360.0) % 1.0, sd
    by_night: dict[str, list] = {}
    for r in rows(cv, "SELECT nights, hump_phase FROM p4_run WHERE "
                      "hump_phase IS NOT NULL AND upper(state)='QUIESCENT' "
                      "AND kind='run'"):
        by_night.setdefault(str(r["nights"]), []).append(r["hump_phase"])
    spreads = {n: _circ(v)[1] for n, v in by_night.items() if len(v) > 1}
    slo2, shi2 = _minmax(spreads.values())
    add("hump phase agreement lo", fmt_float(slo2, 3), unit="cycles",
        source="p4_run",
        note="best within-night filter-to-filter circular phase agreement, "
             "over " + ", ".join(sorted(spreads)))
    add("hump phase agreement hi", fmt_float(shi2, 3), unit="cycles",
        source="p4_run", note="worst of the same set")
    # The two consecutive May nights, whose ephemeris drift is below 0.01
    # cycles: the pair the coherence argument rests on.
    pair = sorted(n for n in by_night if n.startswith("2024-05"))
    if len(pair) == 2:
        m1, m2 = (_circ(by_night[pair[0]])[0], _circ(by_night[pair[1]])[0])
        d = abs(m1 - m2) % 1.0
        shift = min(d, 1.0 - d)
    else:
        shift = None
    add("hump phase night shift", fmt_float(shift, 3), unit="cycles",
        source="p4_run",
        note=(f"{pair[0]} to {pair[1]} " if len(pair) == 2 else "") +
             "shift of the circular mean phase, against <0.01 cycles of "
             "ephemeris drift: coherent within a night, not between nights")

    fk = rows(cv, "SELECT * FROM p4_flicker WHERE sf_excess IS NOT NULL "
                  "AND upper(state)='QUIESCENT'")
    add("flicker bins", fmt_int(len(fk)), source="p4_flicker",
        note="quiescent timescale bins in which the excess over the "
             "measured floor was MEASURABLE; sub-floor bins report not "
             "measured, never zero")
    add("flicker detected", fmt_int(sum(1 for r in fk if r["detected"])),
        source="p4_flicker",
        note="bins clearing 3 sigma on the variance excess")
    add("flicker bins all", fmt_int(one(
        cv, "SELECT count(*) FROM p4_flicker WHERE "
            "upper(state)='QUIESCENT'")), source="p4_flicker",
        note="every quiescent bin the timescale grid contains")
    xlo, xhi = _minmax([r["sf_excess"] for r in fk if r["detected"]], 1000.0)
    add("flicker amp range mmag", fmt_range(xlo, xhi, 0), unit="mmag",
        source="p4_flicker",
        note="over the measured floor, subtracted in quadrature")
    tlo, thi = _minmax([r["tau_s"] for r in fk])
    add("flicker tau range s", fmt_range(tlo, thi, 0), unit="s",
        source="p4_flicker",
        note="timescales the sampling actually populates")

    ob = rows(cv, "SELECT * FROM p4_outburst")
    add("outburst runs", fmt_int(len({r["night"] for r in ob})),
        source="p4_outburst")
    add("outburst run filters", fmt_int(len(ob)), source="p4_outburst")
    sig_ob = [r for r in ob if str(r["rate_verdict"]).upper()
              in ("BRIGHTENING", "FADING")]
    add("outburst rate significant", fmt_int(len(sig_ob)),
        source="p4_outburst",
        note="run-filters whose rate of change clears 3 sigma; the rest "
             "are graded flat and contribute no rate")
    rlo, rhi = _minmax([r["rate_mag_per_h"] for r in sig_ob])
    add("outburst rate range", fmt_range(rlo, rhi, 3), unit="mag/h",
        source="p4_outburst",
        note="the 3-sigma rates only; quoting the flat runs' fitted slopes "
             "beside them would advertise noise as measurement")
    add("outburst peak amp", fmt_float(max(
        (r["amp_above_quiescence"] for r in ob
         if r["amp_above_quiescence"] is not None), default=None), 2),
        unit="mag", source="p4_outburst",
        note="brightest dense run above quiescence")
    add("superoutburst amp", fmt_float(3.0, 1), unit="mag",
        source="ANALYSIS_STRATEGY §4",
        note="the amplitude a superoutburst reaches: the peak above is "
             "1.14 mag short of it")
    blo, bhi = _minmax([r["amp90_blind"] for r in ob], 1000.0)
    add("blind contour range mmag", fmt_range(blo, bhi, 0), unit="mmag",
        source="p4_outburst",
        note="blind-search 90 per cent contour on the outburst runs, "
             "against superhump semi-amplitudes from 50 mmag: this is what "
             "makes 'no superhump' a measurement")
    add("blind contour run filters", fmt_int(sum(
        1 for r in ob if r["amp90_blind"] is not None)),
        source="p4_outburst",
        note="run-filters long enough to place a blind periodogram maximum "
             "and therefore to carry a contour at all; the rest carry none")
    add("superhump floor mmag", fmt_float(_median(
        [r["superhump_floor"] for r in ob], 1000.0), 0), unit="mmag",
        source="p4_outburst", note="lower edge of published superhump "
                                   "semi-amplitudes, as the outburst stage "
                                   "recorded it")
    # Dense RUNS are nights, not run-filters: three filters through one
    # night are one run seen three ways.
    add("yz cnc dense runs", fmt_int(one(
        cv, "SELECT count(DISTINCT nights) FROM p4_run WHERE kind='run'")),
        source="p4_run",
        note="dense runs, counted as nights; none falls inside a "
             "superoutburst")
    add("yz cnc quiescent runs", fmt_int(one(
        cv, "SELECT count(DISTINCT nights) FROM p4_run WHERE kind='run' "
            "AND upper(state)='QUIESCENT'")), source="p4_run")
    add("yz cnc outburst runs", fmt_int(one(
        cv, "SELECT count(DISTINCT night) FROM p4_outburst")),
        source="p4_outburst")

    # -- §5 Results: AN UMa --------------------------------------------
    an = rows(cv, "SELECT * FROM p4_anuma")
    add("an uma capabilities", fmt_int(len(an)), source="p4_anuma")
    add("an uma supported", fmt_int(sum(
        1 for r in an if str(r["verdict"]).upper().startswith("SUPPORTED"))),
        source="p4_anuma")
    # AN UMa's timing precision.  There is NO injection test for this
    # target: p3_sigmat's grid was run on ST LMi's Mode0 night alone.  What
    # exists is (a) the analytic budget, which §4.2 forbids quoting as
    # achieved, and (b) the per-epoch errors of the edges actually fitted.
    # Both are emitted, under names that cannot be mistaken for each other.
    an_edge = rows(cv, "SELECT sigma_t_s, accepted FROM p3_edge "
                       "WHERE target_key='anuma'")
    elo, ehi = _minmax([r["sigma_t_s"] for r in an_edge])
    alo2, ahi2 = _minmax([r["sigma_t_s"] for r in an_edge if r["accepted"]])
    add("an uma edge sigma range s", fmt_range(elo, ehi, 0), unit="s",
        source="p3_edge",
        note="MEASURED per-epoch timing error of every bright-phase edge "
             "fitted for AN UMa; the whole range lies outside the 60 s "
             "threshold")
    add("an uma edge sigma accepted range s", fmt_range(alo2, ahi2, 0),
        unit="s", source="p3_edge",
        note="the same, over the accepted edges only")
    add("an uma edges accepted", fmt_int(sum(1 for r in an_edge
                                             if r["accepted"])),
        source="p3_edge")
    add("an uma edges fitted", fmt_int(len(an_edge)), source="p3_edge")
    # WHY THE EDGES ARE REJECTED, from the stored reasons.  §5.3 blamed a
    # slow filter cycle; the database says the opposite, and the remedy
    # that follows is different in kind.  A faster filter cycle samples
    # more edges; it does not make any one of them stand out of AN UMa's
    # flickering, and 16 of the 20 rejections are exactly that failure.
    rej = [str(r["reason"]) for r in rows(
        cv, "SELECT reason FROM p3_edge WHERE target_key='anuma' AND "
            "accepted=0")]
    snr_v = [float(m.group(1)) for m in
             (re.search(r"step SNR ([0-9.]+)", t) for t in rej) if m]
    add("an uma rejections", fmt_int(len(rej)), source="p3_edge")
    add("an uma rejections snr", fmt_int(len(snr_v)), source="p3_edge",
        note="edges rejected because the step signal-to-noise is below the "
             "bar: the edge is not distinguishable from this star's own "
             "flickering on that cycle")
    add("an uma rejections grid",
        fmt_int(sum(1 for t in rej if "search grid" in t)), source="p3_edge",
        note="edges whose best epoch landed on the boundary of the search "
             "grid, so the epoch is not bracketed")
    add("an uma rejections gap",
        fmt_int(sum(1 for t in rej if " gap," in t)), source="p3_edge",
        note="edges rejected for CADENCE -- the edge fell in a gap wider "
             "than 2.5 times the sampling. Exactly one of twenty, which is "
             "why a faster filter cycle is not the first remedy")
    add("an uma step snr range", fmt_range(*_minmax(snr_v), 1),
        source="p3_edge",
        note="measured step signal-to-noise of the rejected edges")
    add("an uma step snr bar", fmt_int(5), source="ANALYSIS_STRATEGY §4",
        note="step signal-to-noise an edge must reach to be accepted")
    an_ct = rows(ch, "SELECT * FROM ch_timing WHERE target_key='anuma' "
                     "AND night_kind='richest' AND "
                     "regime='per-cycle shape-mismatched'")
    blo2, bhi2 = _minmax([r["sigma_t_s"] for r in an_ct])
    add("an uma analytic budget range s", fmt_range(blo2, bhi2, 0), unit="s",
        source="ch_timing",
        note="ANALYTIC per-cycle budget for AN UMa's richest g night. This "
             "is the same estimator as the ST LMi budget the injection test "
             "overturned, it was never tested against injection for this "
             "target, and §4.2 forbids quoting it as an achieved precision")
    add("an uma injection cells", fmt_int(one(
        cv, "SELECT count(*) FROM p3_sigmat WHERE series_key LIKE "
            "'anuma|%'")), source="p3_sigmat",
        note="ZERO: no injection test was run for this target")
    duty = [r for r in an if r["capability"] == "absolute duty cycle"]
    add("an uma duty halfwidth pp",
        fmt_float(duty[0]["measured"] if duty else None, 0),
        unit="percentage points", source="p4_anuma",
        note="binomial half-width on a 50 per cent duty cycle from the "
             "independent nights below, against a 15 pp bar")
    add("an uma independent nights", fmt_int(one(
        cv, "SELECT max(n_used) FROM p3_state_series WHERE "
            "target_key='anuma'")), source="p3_state_series",
        note="ungated nights the state classification could use")

    # -- §7 Data products: what Table 2 actually shows -------------------
    # cv_series has more rows than Table 2 has lines, and §7 has to say so:
    # the table lists the blocks that produced target photometry, and the
    # unsolved series and the untied 2026 block are counted here instead.
    add("series in table", fmt_int(one(
        cv, "SELECT count(*) FROM cv_series WHERE n_target_points > 0")),
        source="cv_series",
        note="rows of Table 2: every series that produced at least one "
             "target measurement")
    add("series unsolved", fmt_int(one(
        cv, "SELECT count(*) FROM cv_series WHERE status <> 'solved'")),
        source="cv_series", note="series whose ensemble solve did not "
                                 "converge")
    add("series no target points", fmt_int(one(
        cv, "SELECT count(*) FROM cv_series WHERE n_target_points = 0 "
            "OR n_target_points IS NULL")), source="cv_series",
        note="series carrying no target measurement, solved or not")

    # -- Verdicts -------------------------------------------------------
    vd = rows(cv, "SELECT * FROM p4_verdict")
    add("verdicts", fmt_int(len(vd)), source="p4_verdict")

    return N


# ===========================================================================
# Tables.  Same law: a table of measurements is emitted, never typed.
# ===========================================================================
def _esc(text) -> str:
    """LaTeX-safe text for a table cell.

    Comparison-star counts and verdict strings arrive from the database
    containing underscores, per-cent signs and ampersands, every one of
    which ends a ``tectonic`` run with an error a hundred lines from the
    actual cause.
    """
    s = "" if text is None else str(text)
    for a, b in (("\\", "\\textbackslash{}"), ("&", "\\&"), ("%", "\\%"),
                 ("$", "\\$"), ("#", "\\#"), ("_", "\\_"), ("{", "\\{"),
                 ("}", "\\}"), ("~", "\\textasciitilde{}"),
                 ("^", "\\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def render_series_table(cv: sqlite3.Connection) -> str:
    """The per-series census: what was observed and how well it solved.

    One row per solved series, with the numbers a referee checks first --
    frames, points, comparison and check stars, achieved scatter, the
    error-bar inflation, and the catalogue-tie verdict.  Untied and
    unsolved series appear too, with their status, because a table that
    silently omits the failures is a different claim than the data support.
    """
    rr = rows(cv, """
        SELECT s.series_key, s.target_key, s.era_id, s.filter,
               s.n_frames_used, s.n_target_points, s.n_comp, s.n_check,
               s.check_rms_median, s.chi2_inflation, s.status,
               c.verdict AS tie_verdict, c.check_rms_clip
        FROM cv_series s
        LEFT JOIN cv_cattie c ON c.series_key = s.series_key
                             AND c.is_primary = 1
        WHERE s.n_target_points > 0
        ORDER BY s.target_key, s.era_id, s.filter
    """)
    n_all_series = one(cv, "SELECT count(*) FROM cv_series")
    n_unsolved = one(cv, "SELECT count(*) FROM cv_series WHERE "
                         "status <> 'solved'")
    # Both tie statistics, so the caption names which column this is.
    _tied = rows(cv, "SELECT check_rms, check_rms_clip FROM cv_cattie "
                     "WHERE is_primary=1 AND verdict LIKE 'TIED%'")
    _med_clip = _median([r["check_rms_clip"] for r in _tied]) or 0.0
    _med_raw = _median([r["check_rms"] for r in _tied]) or 0.0
    # THE TWO CHECK-STAR POPULATIONS, NAMED, IN THE TABLE THAT PRINTS BOTH.
    # $N_{\rm chk}$ is the four stars held out of the ENSEMBLE SOLVE.  The
    # "Tie acc." column is measured on a different and much larger set: the
    # catalogue stars held out of the TIE FIT, 15--513 per block.  Adjacent
    # columns, one name, two counts -- a reader cross-checking §3.2's
    # "15--513" against this table's "4" cannot reconcile them unless the
    # caption says they count different stars.
    _tk = rows(cv, "SELECT n_check FROM cv_cattie WHERE is_primary=1 "
                   "AND verdict LIKE 'TIED%'")
    _tk_lo, _tk_hi = _minmax([r["n_check"] for r in _tk])
    _sk_chk = sorted({r["n_check"] for r in rr if r["n_check"] is not None})
    # The solved-but-no-tied-points blocks, NAMED.  An earlier caption
    # described both as "the target is undetected", which is false of
    # ST LMi's 2024 y block: the star IS detected there on every frame,
    # and what is missing is the ensemble zero point and therefore the tie.
    _no_pts = rows(cv, """
        SELECT s.series_key, s.target_key, s.era_id, s.filter,
               s.n_target_rows, c.verdict AS tie_verdict
        FROM cv_series s
        LEFT JOIN cv_cattie c ON c.series_key = s.series_key
                             AND c.is_primary = 1
        WHERE s.n_target_points = 0 AND s.status = 'solved'
        ORDER BY s.target_key, s.era_id, s.filter""")

    def _blk(r):
        return (f"{TARGET_LABEL.get(r['target_key'], r['target_key'])}'s "
                f"{ERA_LABEL.get(r['era_id'], 'e' + str(r['era_id']))} "
                f"${r['filter']}$ block")
    _no_pts_txt = "; ".join(
        f"{_blk(r)} (" + (
            f"{r['n_target_rows']} instrumental target detections that the "
            f"solve could not place on a zero point, so the block never "
            f"reached the tie stage"
            if (r["n_target_rows"] or 0) > 0 and not r["tie_verdict"]
            else "the target is not detected in it at all") + ")"
        for r in _no_pts)
    out = [
        "\\begin{deluxetable*}{llccrrrrrrl}",
        "\\tablecaption{Per-series photometric census: the "
        f"{len(rr)} target--era--filter blocks that produced target "
        "photometry. Every column is a query against "
        "\\texttt{cv\\_series} and \\texttt{cv\\_cattie}; nothing in this "
        "table was typed. "
        "TWO COLUMNS OF THIS TABLE COUNT DIFFERENT HELD-OUT STARS, and the "
        "paper uses different names for them throughout. $N_{\\rm chk}$ is "
        "the number of SOLVE check stars, withheld from the ensemble "
        f"solution ({fmt_range(_sk_chk[0], _sk_chk[-1], 0) if _sk_chk else '?'}"
        " in every block here), and $\\sigma_{\\rm chk}$ is their median "
        "scatter. The `Tie acc.' column is not measured on those stars: it "
        "is the SIGMA-CLIPPED residual RMS of the TIE check stars "
        "(\\texttt{check\\_rms\\_clip}), the catalogue stars withheld from "
        f"the tie fit, of which each block holds out "
        f"{fmt_range(_tk_lo, _tk_hi, 0)}. Over the tied blocks its median is "
        f"{1000 * _med_clip:.0f}~mmag against {1000 * _med_raw:.0f}~mmag "
        "unclipped, and Section~\\ref{sec:tie} quotes both. $I$ is the "
        "ratio of achieved scatter to the formal error bar. "
        f"The release holds {n_all_series} series in all; the remaining "
        f"{n_all_series - len(rr)} produced no catalogue-tied target "
        f"measurement ({n_unsolved} ensemble solves that did not converge, "
        f"and {len(_no_pts)} solved series: {_no_pts_txt}; "
        "Section~\\ref{sec:vvpupeuuma}). "
        "\\label{tab:series}}",
        "\\tablehead{\\colhead{Target} & \\colhead{Camera} & "
        "\\colhead{Band} & \\colhead{Frames} & \\colhead{Points} & "
        "\\colhead{$N_{\\rm comp}$} & \\colhead{$N_{\\rm chk}$} & "
        "\\colhead{$\\sigma_{\\rm chk}$} & \\colhead{$I$} & "
        "\\colhead{Tie acc.} & \\colhead{Tie verdict}\\\\"
        "\\colhead{} & \\colhead{} & \\colhead{} & \\colhead{} & "
        "\\colhead{} & \\colhead{} & \\colhead{} & \\colhead{(mmag)} & "
        "\\colhead{} & \\colhead{(mmag)} & \\colhead{}}",
        "\\startdata",
    ]
    for r in rr:
        cells = [
            _esc(TARGET_LABEL.get(r["target_key"], r["target_key"])),
            _esc(ERA_LABEL.get(r["era_id"], f"e{r['era_id']}")),
            _esc(r["filter"]),
            fmt_int(r["n_frames_used"]) or "\\nodata",
            fmt_int(r["n_target_points"]) or "\\nodata",
            fmt_int(r["n_comp"]) or "\\nodata",
            fmt_int(r["n_check"]) or "\\nodata",
            fmt_float(None if r["check_rms_median"] is None
                      else 1000 * r["check_rms_median"], 0) or "\\nodata",
            fmt_float(r["chi2_inflation"], 2) or "\\nodata",
            fmt_float(None if r["check_rms_clip"] is None
                      else 1000 * r["check_rms_clip"], 0) or "\\nodata",
            _esc(r["tie_verdict"] or "no tie"),
        ]
        out.append(" & ".join(cells) + " \\\\")
    out += ["\\enddata", "\\end{deluxetable*}", ""]
    return "\n".join(out)


def render_anuma_table(cv: sqlite3.Connection) -> str:
    """AN UMa, capability by capability, per filter.

    One measured number against one stated bar per row: a reader who
    disagrees can disagree with the bar in one line, instead of with a
    composite score nobody can audit.
    """
    rr = rows(cv, "SELECT * FROM p4_anuma ORDER BY rank, filter")
    out = [
        "\\begin{deluxetable}{llrrll}",
        "\\tablecaption{AN UMa graded capability by capability and filter "
        "by filter. Each row is one measured value against one stated "
        "threshold, emitted from \\texttt{p4\\_anuma}. The morphology rows "
        "count \\emph{usable} full-orbit nights, which is a stricter "
        "criterion than Section~\\ref{sec:obs}'s census: the census counts "
        "nights whose SAMPLING covers an orbit, while a folded panel "
        "additionally requires the modulation to be detected in that band, "
        "which is why AN~UMa's six full-orbit $i$ nights yield none. "
        "The timing capability is graded on TWO rows per band, because it "
        "fails for two independent reasons and one row can carry only one "
        "number: the accepted-edge count, graded against a bar of accepted "
        "\\emph{per-cycle} edges even though Section~\\ref{sec:timing} "
        "publishes no per-cycle epoch for any target --- the per-cycle "
        "edges are the inputs a per-night epoch is averaged from, so a "
        "target with too few of them has no epoch to publish --- and the "
        "one-feature test of Section~\\ref{sec:oc}, the circular scatter of "
        "the accepted edges in orbital phase. "
        "\\label{tab:anuma}}",
        "\\tablehead{\\colhead{Capability} & \\colhead{Band} & "
        "\\colhead{Measured} & \\colhead{Bar} & \\colhead{Unit} & "
        "\\colhead{Verdict}}",
        "\\startdata",
    ]
    for r in rr:
        out.append(" & ".join([
            _esc(r["capability"]), _esc(r["filter"]),
            fmt_float(r["measured"], 2) or "\\nodata",
            fmt_float(r["bar"], 2) or "\\nodata",
            _esc(r["unit"]), _esc(r["verdict"])]) + " \\\\")
    out += ["\\enddata", "\\end{deluxetable}", ""]
    return "\n".join(out)


def render_verdict_table(cv: sqlite3.Connection) -> str:
    """The headline calls, each with the one number that decided it."""
    rr = rows(cv, "SELECT * FROM p4_verdict ORDER BY rank")
    out = [
        "\\begin{deluxetable*}{p{0.30\\textwidth}lp{0.42\\textwidth}}",
        "\\tablecaption{The closing decisions of this programme, each with "
        "the single measured number that decided it, from "
        "\\texttt{p4\\_verdict}. \\label{tab:verdicts}}",
        "\\tablehead{\\colhead{Question} & \\colhead{Verdict} & "
        "\\colhead{Deciding number}}",
        "\\startdata",
    ]
    for r in rr:
        out.append(" & ".join([
            _esc(r["question"]), _esc(r["verdict"]),
            _esc(r["deciding_number"])]) + " \\\\")
    out += ["\\enddata", "\\end{deluxetable*}", ""]
    return "\n".join(out)


def render_instrument_table(man: sqlite3.Connection) -> str:
    """The detector constants every veto and every error bar rests on."""
    rr = rows(man, """SELECT * FROM s2_ceiling_modes
                      WHERE clip_adu IS NOT NULL ORDER BY n_frames DESC""")
    ratios = [r["veto_adu"] / r["clip_adu"] for r in rr
              if r["clip_adu"] and r["veto_adu"]]
    vlo, vhi = (min(ratios), max(ratios)) if ratios else (0.0, 0.0)
    out = [
        "\\begin{deluxetable}{lrrrr}",
        "\\tablecaption{Measured detector limits per readout mode, from "
        "\\texttt{s2\\_ceiling\\_modes}. The clip is the pileup edge of "
        "the pixel histograms of the frames counted in column~2; the veto "
        "is $0.92$ of that clip rounded DOWN to the nearest 100~ADU, so "
        f"the realised ratio is {vlo:.3f}--{vhi:.3f}. No threshold in "
        "this paper was chosen by eye. \\label{tab:instrument}}",
        "\\tablehead{\\colhead{Readout mode} & \\colhead{Frames} & "
        "\\colhead{Hard max} & \\colhead{Clip} & \\colhead{Veto}\\\\"
        "\\colhead{} & \\colhead{} & \\colhead{(ADU)} & \\colhead{(ADU)} & "
        "\\colhead{(ADU)}}",
        "\\startdata",
    ]
    for r in rr:
        out.append(" & ".join([
            _esc(r["mode"]), fmt_int(r["n_frames"]) or "\\nodata",
            fmt_int(r["hard_max_adu"]) or "\\nodata",
            fmt_int(r["clip_adu"]) or "\\nodata",
            fmt_int(r["veto_adu"]) or "\\nodata"]) + " \\\\")
    out += ["\\enddata", "\\end{deluxetable}", ""]
    return "\n".join(out)


def render_tables(cv: sqlite3.Connection, man: sqlite3.Connection,
                  stamp: str = "") -> str:
    """``tables.tex``: every measured table in the paper, in one file."""
    head = [
        "%% tables.tex -- GENERATED FILE.  DO NOT EDIT.",
        "%% Emitted by pipeline/macro_phot/numbers_cv.py.  Each table is a "
        "query result;",
        "%% editing a cell here would put a number in the paper that the "
        "database does not hold.",
    ]
    if stamp:
        head.append(f"%% {stamp}")
    head.append("")
    return "\n".join(head) + "\n".join([
        render_instrument_table(man),
        render_series_table(cv),
        render_anuma_table(cv),
        render_verdict_table(cv),
    ])
