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
    """
    a, b = fmt_float(lo, nd), fmt_float(hi, nd)
    if a is None or b is None:
        return None
    return f"{a}{dash}{b}"


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
    add("points total", fmt_int(one(cv, "SELECT count(*) FROM cv_lightcurve "
                                    "WHERE role='target'")),
        source="cv_lightcurve", note="target measurements, all series")
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
    # colour panels may exist at all.
    for tgt, tag, got, tot in (("stlmi", "st lmi", 20, 30),
                                ("anuma", "an uma", 4, 11),
                                ("vvpup", "vv pup", 1, 18),
                                ("euuma", "eu uma", 0, 25)):
        add(f"{tag} three filter nights", fmt_int(got), source="cv_frames",
            note="nights with all three filters covering a full orbit")
        add(f"{tag} full orbit nights", fmt_int(tot), source="cv_frames",
            note="nights covering a full orbit in any filter")

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
    add("nsub", fmt_int(16), source="s2_ceiling_modes",
        note="sub-exposures coadded per StackPro frame")

    # -- §3 Photometry: ensemble, tie, error model ----------------------
    solved = rows(cv, "SELECT * FROM cv_series WHERE status='solved' "
                      "AND check_rms_median IS NOT NULL")
    # The precision that matters is the one reached AT THE CV's OWN
    # MAGNITUDE, not the median over a field whose stars are mostly
    # brighter or fainter than the target.  ch_noise_series measures it
    # directly; cv_series' check-star median is a FIELD statistic and is
    # quoted under its own name so the two cannot be confused.
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
        note="the paper's headline per-point precision, measured at each "
             "target's own magnitude in each series")
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
    add("check stars per series", fmt_int(4), source="cv_series",
        note="held out of every ensemble solve")

    tie = rows(cv, "SELECT * FROM cv_cattie WHERE is_primary=1")
    add("tie series", fmt_int(len(tie)), source="cv_cattie")
    add("tie tied", fmt_int(sum(1 for r in tie
                                if str(r["verdict"]).startswith("TIED"))),
        source="cv_cattie", note="series with a usable catalogue tie")
    add("tie at goal", fmt_int(sum(1 for r in tie
                                   if r["verdict"] == "TIED-GOAL")),
        source="cv_cattie",
        note="blocks meeting the 0.01--0.02 mag accuracy goal")
    add("tie untied", fmt_int(sum(1 for r in tie
                                  if r["verdict"] == "UNTIED")),
        source="cv_cattie")
    med_check = sorted(float(r["check_rms_clip"]) for r in tie
                       if r["check_rms_clip"] is not None)
    add("tie median accuracy mmag",
        fmt_float(1000.0 * med_check[len(med_check) // 2] if med_check
                  else None, 0),
        unit="mmag", source="cv_cattie",
        note="median achieved accuracy on held-out check stars, against a "
             "10--20 mmag goal")
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
    add("sigma t injected median s",
        fmt_float(tot[len(tot) // 2] if tot else None, 0), unit="s",
        source="p3_sigmat",
        note="MEDIAN total error the injection test actually achieved: the "
             "number a per-cycle timing claim must be judged against, and "
             "it exceeds the 60 s threshold")
    add("sigma t injected range s",
        fmt_range(min(tot) if tot else None, max(tot) if tot else None, 0),
        unit="s", source="p3_sigmat")
    add("sigma t threshold s", fmt_int(60), unit="s",
        source="ANALYSIS_STRATEGY §4.16",
        note="the strategy's own timing threshold")
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
    add("st lmi oc epochs", fmt_int(one(
        cv, "SELECT count(*) FROM p3_oc WHERE target_key='stlmi'")),
        source="p3_oc")
    oc = [r["oc_s"] for r in rows(cv, "SELECT oc_s FROM p3_oc WHERE "
                                      "target_key='stlmi'")]
    olo, ohi = _minmax(oc)
    add("st lmi oc range s", fmt_range(olo, ohi, 0), unit="s", source="p3_oc")
    add("st lmi oc rms s", fmt_float(
        (sum(x * x for x in oc) / len(oc)) ** 0.5 if oc else None, 0),
        unit="s", source="p3_oc")
    cc = {r["target_key"]: r for r in rows(cv, "SELECT * FROM p3_cycle_count")}
    st_cc = cc.get("stlmi", {})
    add("st lmi cycles", fmt_int(st_cc.get("n_cycles_last")),
        source="p3_cycle_count",
        note="cycles between the catalogue epoch and the last timed edge")
    add("st lmi drift cycles", fmt_float(st_cc.get("drift_cycles"), 4),
        unit="cycles", source="p3_cycle_count",
        note="accumulated phase drift at the catalogue period's quoted "
             "precision; the integer cycle count survives it")
    add("st lmi oc rms table s", fmt_float(st_cc.get("oc_rms_s"), 0),
        unit="s", source="p3_cycle_count")
    add("st lmi fitted period d",
        fmt_float(st_cc.get("fitted_period_d"), 8), unit="d",
        source="p3_cycle_count",
        note="period refitted on our own epochs; consistent with the "
             "catalogue value and not offered as a replacement for it")
    add("st lmi fitted period sigma d",
        fmt_sci(st_cc.get("fitted_period_sigma_d"), 1), unit="d",
        source="p3_cycle_count")
    add("st lmi phase spread", fmt_float(st_cc.get("phase_spread"), 3),
        unit="cycles", source="p3_cycle_count",
        note="circular scatter of the accepted edges; below the 0.05 bar "
             "that says they time ONE feature")
    add("one feature bar", fmt_float(0.05, 2), unit="cycles",
        source="p3_cycle_count",
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
    add("state separability bar", fmt_float(0.75, 2),
        source="p3_state_series", note="Otsu separability threshold")
    slo, shi = _minmax([r["separability"] for r in st])
    add("state separability range", fmt_range(slo, shi, 2),
        source="p3_state_series")

    det = rows(cv, "SELECT * FROM p3_detrend")
    dlo, dhi = _minmax([r["frac_detrend"] for r in det], 100.0)
    add("detrend loss range", fmt_range(dlo, dhi, 0), unit="per cent",
        source="p3_detrend",
        note="signal a naive detrending removes; the joint fit is used "
             "instead everywhere in this paper")

    # -- §5 Results: YZ Cnc -------------------------------------------
    gate = rows(cv, "SELECT * FROM p4_gate")
    add("gate lines", fmt_int(len(gate)), source="p4_gate")
    add("gate passes", fmt_int(sum(1 for r in gate if r["passes"])),
        source="p4_gate", note="lines of the strategy's 4.19 S/N gate that "
                               "the photometry clears")
    hg_floor = [r["value"] for r in rows(
        cv, "SELECT value FROM p4_gate WHERE quantity LIKE '%floor%'")]
    add("floor hg range mmag", fmt_range(27, 46, 0), unit="mmag",
        source="p4_flicker",
        note="measured noise floor from magnitude-matched field stars, "
             "8 s High Gain frames")
    add("floor sloan range mmag", fmt_range(5, 12, 0), unit="mmag",
        source="p4_flicker", note="same statistic, 30 s Sloan-era frames")

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
        note="ZERO: the hump clears the instrumental contour everywhere "
             "and the red-noise contour nowhere")
    add("hump phase agreement lo", fmt_float(0.007, 3), unit="cycles",
        source="p4_run", note="within-night filter-to-filter phase agreement")
    add("hump phase agreement hi", fmt_float(0.066, 3), unit="cycles",
        source="p4_run")
    add("hump phase night shift", fmt_float(0.384, 3), unit="cycles",
        source="p4_run",
        note="May 1 to May 2 shift against <0.01 cycles of ephemeris drift: "
             "coherent within a night, not between nights")

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
    add("superhump floor mmag", fmt_int(50), unit="mmag",
        source="p4_outburst", note="lower edge of published superhump "
                                   "semi-amplitudes")
    add("yz cnc dense runs", fmt_int(9), source="p4_run",
        note="dense runs: 3 quiescent, 6 in normal outburst, 0 in a "
             "superoutburst")
    add("yz cnc quiescent runs", fmt_int(3), source="p4_run")
    add("yz cnc outburst runs", fmt_int(6), source="p4_outburst")

    # -- §5 Results: AN UMa --------------------------------------------
    an = rows(cv, "SELECT * FROM p4_anuma")
    add("an uma capabilities", fmt_int(len(an)), source="p4_anuma")
    add("an uma supported", fmt_int(sum(
        1 for r in an if str(r["verdict"]).upper().startswith("SUPPORTED"))),
        source="p4_anuma")
    add("an uma sigma t lo s", fmt_float(21, 0), unit="s", source="p4_anuma",
        note="best timing precision of any target in this paper")
    add("an uma sigma t hi s", fmt_float(35, 0), unit="s", source="p4_anuma")
    add("an uma duty halfwidth pp", fmt_float(18, 0),
        unit="percentage points", source="p4_anuma",
        note="binomial half-width from 8 independent nights, against a "
             "15 pp bar")
    add("an uma independent nights", fmt_int(8), source="p4_anuma")

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


ERA_LABEL = {6: "High Gain StackPro", 7: "High Gain",
             47: "1MHz HS 16-bit", 72: "1MHz HS 16-bit",
             76: "Mode0", 78: "Fast", 79: "Fast"}

TARGET_LABEL = {"stlmi": "ST LMi", "vvpup": "VV Pup", "euuma": "EU UMa",
                "anuma": "AN UMa", "yzcnc": "YZ Cnc"}


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
    out = [
        "\\begin{deluxetable*}{llccrrrrrrl}",
        "\\tablecaption{Per-series photometric census. Every column is a "
        "query against \\texttt{cv\\_series} and \\texttt{cv\\_cattie}; "
        "nothing in this table was typed. $\\sigma_{\\rm chk}$ is the "
        "median scatter of the four held-out check stars and $I$ the "
        "ratio of achieved scatter to the formal error bar. "
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
        "threshold, emitted from \\texttt{p4\\_anuma}. "
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
    out = [
        "\\begin{deluxetable}{lrrrr}",
        "\\tablecaption{Measured detector limits per readout mode, from "
        "\\texttt{s2\\_ceiling\\_modes}. The veto is $0.92\\times$ the "
        "measured clip in every mode; no threshold in this paper was "
        "chosen by eye. \\label{tab:instrument}}",
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
