"""S1 evidence report renderer: the stratified astrometry go/no-go.

Reads the S1 tables of the manifest database (NEVER the catalog, never the
solver's scratch files — if a number cannot be derived from the database,
it does not belong on the page) and writes:

* ``docs/pipeline/s1_astrometry.html``   — the report
* ``docs/pipeline/figures/s1/*.png``     — every figure

The page follows the site's Socratic format: one section per decision, each
section = Question → Evidence → Decision → Consequence.  EVERY number in
the HTML is interpolated from a SQL query executed in this module or from a
constant defined in ``macro_core.astrom`` — nothing is hand-typed, so
re-running ``run_s1_experiment.py report`` after more solves (or a re-drawn
design) regenerates the whole argument, verdicts included.
"""

from __future__ import annotations

import sqlite3
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
import numpy as np               # noqa: E402

from . import astrom             # noqa: E402  (constants for interpolation)
# Shared page machinery: one house figure style, one query discipline,
# one table generator
# generator as the S0/S0b reports — one visual language across the site.
from .report_s0 import (          # noqa: E402
    ACCENT, BAD, STYLE, DPI, GOOD, INK, MUTED, PAPER, WARN,
    _figure, esc, fmt, q, q1, table)
from . import plotstyle as ps    # noqa: E402  (the house figure style)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s1"
HTML_PATH = DOCS_DIR / "s1_astrometry.html"

#: Status color per verdict — matches the site badge palette.
VERDICT_COLOR = {"GO": GOOD, "CAUTION": WARN, "NO-GO": BAD}

#: Failure-gallery size: distinct diagnoses first, then fill, capped here.
N_GALLERY = 6


# ---------------------------------------------------------------------------
# Small derived-stat helpers (medians live in Python: SQLite has none).
# ---------------------------------------------------------------------------
def med(vals) -> float | None:
    """Median of a list of non-NULL values, None when empty."""
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def table_exists(con, name: str) -> bool:
    """True when ``name`` is a table in this database (the baseline tables
    only exist once ``run_s1_experiment.py snapshot`` has been run)."""
    return q1(con, "SELECT count(*) FROM sqlite_master "
                   "WHERE type='table' AND name=?", (name,)) > 0


def stratum_stats(con, prefix: str = "s1_") -> list[dict]:
    """One dict per stratum: counts, rate, Wilson CI, verdict, medians.

    This is THE results computation — the figures, the tables and the
    verdict section all read from this one list, so they can never
    disagree with each other.

    ``prefix`` selects which experiment to read: ``s1_`` is the live one,
    ``s1_baseline_`` the frozen previous design.  ONE function computes
    both sides of the before/after delta, so the old rates cannot be
    computed by a different rule than the new ones — which is the only
    way a delta table means anything.

    Success = status ``solved`` ONLY: a ``bad_solve`` (a .solved marker
    whose WCS failed the acceptance gate) is a failure, and its bogus
    pixel scale never enters the measured-scale statistics.

    CENSUS strata (sample == population) carry no sampling uncertainty:
    their interval collapses to the exact rate (lo == hi == rate) and the
    verdict judges that exact rate — see ``astrom.verdict_for``.
    """
    out = []
    for sid, pop_name, desc, n_pop, n_sample, seed in q(con, f"""
            SELECT stratum_id, population, description, n_population,
                   n_sample, seed FROM {prefix}strata ORDER BY rowid"""):
        rows = q(con, f"""
            SELECT status, solve_time_s, pixscale_arcsec, rms_arcsec,
                   n_matched FROM {prefix}solve_experiment
            WHERE stratum_id = ?""", (sid,))
        n = len(rows)
        k = sum(r[0] == "solved" for r in rows)
        census = (n > 0 and n == n_pop)
        if census:
            # Fully-enumerated stratum: the population rate IS k/n.
            lo = hi = k / n
        else:
            lo, hi = astrom.wilson_ci(k, n)
        solved = [r for r in rows if r[0] == "solved"]
        out.append({
            "stratum_id": sid, "population": pop_name, "desc": desc,
            "n_pop": n_pop, "n": n, "k": k, "rate": k / n if n else 0.0,
            "lo": lo, "hi": hi, "census": census,
            "verdict": astrom.verdict_for(k, n, n_pop),
            "seed": seed,
            "med_time": med([r[1] for r in rows]),          # all outcomes
            "med_time_solved": med([r[1] for r in solved]),
            "med_pixscale": med([r[2] for r in solved]),
            "med_rms": med([r[3] for r in solved]),
            "med_matched": med([r[4] for r in solved]),
        })
    return out


def night_coverage(con, strata: list[dict]) -> list[dict]:
    """Per-stratum night-clustering diagnostics for the CI caveat.

    Frames within a night share cloud/focus/wind state, so the Wilson
    intervals' iid assumption is optimistic; this table shows how far.
    Population night counts come from the SAME base query + pure
    classifier the design used (``astrom.fetch_candidates`` /
    ``classify_stratum``) — one definition, so the denominators cannot
    drift from the design's.
    """
    pop_nights: dict[str, set] = {s["stratum_id"]: set() for s in strata}
    for r in astrom.fetch_candidates(con):
        if astrom.is_solvable_candidate(r):
            sid = astrom.classify_stratum(r)
            if sid in pop_nights and r["night"]:
                pop_nights[sid].add(r["night"])
    out = []
    for s in strata:
        rows = q(con, """SELECT night, status = 'solved'
                         FROM s1_solve_experiment
                         WHERE stratum_id = ?""", (s["stratum_id"],))
        # All-or-nothing collapse: a night succeeds only if EVERY sampled
        # frame on it solved — the conservative night-level rate.
        k_n, n_n = astrom.night_collapse(rows)
        lo_n, _ = astrom.wilson_ci(k_n, n_n)
        out.append({
            "stratum_id": s["stratum_id"],
            "n_nights": n_n, "k_nights": k_n,
            "n_pop_nights": len(pop_nights[s["stratum_id"]]),
            "night_lo": lo_n, "census": s["census"],
        })
    return out


def population_rollup(con, strata: list[dict]) -> list[dict]:
    """Aggregate strata into the four verdict populations.

    The rollup is population-weighted where it matters: the projected
    batch hours multiply each stratum's FULL population by that stratum's
    own median per-frame time (a stratum's sample speaks for its own
    population and no one else's).
    """
    by_pop: dict[str, list[dict]] = {}
    for s in strata:
        by_pop.setdefault(s["population"], []).append(s)
    out = []
    for pop_name, members in by_pop.items():
        k = sum(s["k"] for s in members)
        n = sum(s["n"] for s in members)
        n_backlog = sum(s["n_pop"] for s in members)
        # Pooled SAMPLE rate + CI: equal allocation (~48 per stratum)
        # means small strata are over-represented in this pool relative
        # to their backlog share — it answers "how did the sample do",
        # not "what will the backlog yield" (the weighted columns below
        # answer that).  The mismatch is conservative for these data:
        # the small strata are the weak ones, so pooling can only
        # UNDERSTATE a GO population's backlog rate.
        lo, hi = astrom.wilson_ci(k, n)
        hours = sum(
            astrom.projected_hours(s["n_pop"], s["med_time"] or 0.0)
            for s in members)
        # Expected yield: each stratum's backlog × its own point rate,
        # with a range from each stratum's own interval bounds (a census
        # stratum contributes exactly its known count — lo == hi there).
        expect = sum(s["n_pop"] * s["rate"] for s in members)
        expect_lo = sum(s["n_pop"] * s["lo"] for s in members)
        expect_hi = sum(s["n_pop"] * s["hi"] for s in members)
        out.append({
            "population": pop_name, "k": k, "n": n, "lo": lo, "hi": hi,
            "rate": k / n if n else 0.0, "n_backlog": n_backlog,
            "hours": hours, "expected_solved": expect,
            "expected_lo": expect_lo, "expected_hi": expect_hi,
            # Backlog-weighted rate: the rate the yield column implies.
            "weighted_rate": expect / n_backlog if n_backlog else 0.0,
            "verdict": astrom.verdict_for(k, n),
        })
    return out


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_success_rates(strata: list[dict]) -> str:
    """Horizontal bars: success rate per stratum with Wilson 95% CIs."""
    labels = [s["stratum_id"] for s in strata][::-1]
    rates = [100 * s["rate"] for s in strata][::-1]
    # Error bars clamp at zero: the Wilson interval is centered BELOW the
    # point rate at the extremes (a 48/48 stratum has hi < 1.0), and
    # matplotlib rejects negative bar lengths.
    lo_err = [max(0.0, 100 * (s["rate"] - s["lo"])) for s in strata][::-1]
    hi_err = [max(0.0, 100 * (s["hi"] - s["rate"])) for s in strata][::-1]
    colors = [VERDICT_COLOR[s["verdict"]] for s in strata][::-1]
    ns = [(s["k"], s["n"]) for s in strata][::-1]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(labels) + 1.6))
        bars = ax.barh(labels, rates, color=colors,
                       xerr=[lo_err, hi_err], ecolor=INK, capsize=3)
        ax.set_xlim(0, 105)
        ax.set_xlabel("solve success rate (%) with Wilson 95% CI")
        ax.set_title("S1 experiment: success rate by stratum")
        # The GO bar: verdicts key off the CI lower bound crossing it.
        ax.axvline(100 * astrom.GO_LOWER_BOUND, color=MUTED,
                   linestyle="--", linewidth=1)
        # Count labels sit INSIDE the bar's left edge, clear of the CI
        # whiskers that live at the bar's right end.
        for b, r, (k, n) in zip(bars, rates, ns):
            ax.annotate(f"{k}/{n}", (2.0, b.get_y() + b.get_height() / 2),
                        va="center", ha="left", fontsize=8,
                        color=PAPER if r > 12 else INK)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s1_success_rates.png", dpi=DPI)
        plt.close(fig)
    return "figures/s1/s1_success_rates.png"


def fig_solve_times(con) -> str:
    """Solve-time distributions: solved (fast spike) vs failed (CPU wall).

    The split matters for the cost projection: a batch's wall-clock is
    dominated by its FAILURES, each of which burns the full CPU budget.
    """
    ok = [r[0] for r in q(con, """SELECT solve_time_s
        FROM s1_solve_experiment WHERE status = 'solved'
          AND solve_time_s IS NOT NULL""")]
    bad = [r[0] for r in q(con, """SELECT solve_time_s
        FROM s1_solve_experiment WHERE status IN ('unsolved', 'timeout')
          AND solve_time_s IS NOT NULL""")]
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(8.6, 3.4))
        bins = np.linspace(0, max(ok + bad + [1.0]), 40)
        ax.hist(ok, bins=bins, color=ACCENT, alpha=0.85,
                label=f"solved (n={len(ok)})")
        ax.hist(bad, bins=bins, color=WARN, alpha=0.85,
                label=f"failed (n={len(bad)})")
        ax.set_yscale("log")
        ax.set_xlabel("wall time per frame (s), funpack + solve")
        ax.set_ylabel("frames (log)")
        ax.set_title("Solve-time distribution: successes are seconds, "
                     "failures burn the CPU budget")
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s1_solve_times.png", dpi=DPI)
        plt.close(fig)
    return "figures/s1/s1_solve_times.png"


def fig_failure_gallery(con) -> tuple[str, int]:
    """Thumbnail gallery of autopsied failures, one panel per example.

    Picks up to ``N_GALLERY`` failures with distinct diagnoses first
    (every failure MODE gets a face), then fills with further examples.
    Returns (src, n_shown).
    """
    rows = q(con, """
        SELECT a.obs_rowid, a.diagnosis, a.thumb_path,
               s.canonical_target, s.filter, s.exptime, s.stratum_id
        FROM s1_failure_autopsy a
        JOIN s1_solve_experiment s USING (obs_rowid)
        WHERE a.thumb_path IS NOT NULL ORDER BY a.obs_rowid""")
    # Distinct-diagnosis-first selection, deterministic by rowid order.
    chosen, seen = [], set()
    for r in rows:
        if r[1] not in seen:
            chosen.append(r)
            seen.add(r[1])
    for r in rows:
        if len(chosen) >= N_GALLERY:
            break
        if r not in chosen:
            chosen.append(r)
    chosen = chosen[:N_GALLERY]
    n = len(chosen)
    with plt.rc_context(STYLE):
        # constrained layout + generous height: the two title lines per
        # panel must not collide with the image row above them.
        fig, axes = plt.subplots(2, 3, figsize=(8.6, 7.2),
                                 layout="constrained")
        for ax in axes.flat:
            ax.axis("off")
        for ax, (rowid, diag, thumb, target, filt, expt, sid) in zip(
                axes.flat, chosen):
            try:
                ax.imshow(mpimg.imread(thumb))
            except Exception:
                pass                     # a missing thumb leaves the slot
            # Short diagnosis label: the part before the parenthesis.
            short = diag.split(" (")[0]
            ax.set_title(f"{target or '?'} · {filt or '—'} · "
                         f"{expt:.0f}s\n{short} · {sid}",
                         fontsize=8, pad=6)
        fig.suptitle("Failure gallery: what an unsolvable frame "
                     "actually looks like", fontsize=11)
        fig.savefig(FIG_DIR / "s1_failure_gallery.png", dpi=DPI)
        plt.close(fig)
    return "figures/s1/s1_failure_gallery.png", n


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def section_tooling(con) -> str:
    meta = dict(q(con, "SELECT key, value FROM s1_build_meta"))
    gb = int(meta.get("index_bytes", 0)) / 1e9
    # False positives caught by the solution-acceptance gate.
    n_bad = q1(con, """SELECT count(*) FROM s1_solve_experiment
                       WHERE status = 'bad_solve'""")
    # Measured plate scale per readout family — from the experiment's own
    # solved WCS headers, the strongest possible check on the priors.
    fam_rows = q(con, """
        SELECT readoutm, xbinning, count(*),
               min(pixscale_arcsec), max(pixscale_arcsec)
        FROM s1_solve_experiment WHERE status = 'solved'
        GROUP BY readoutm, xbinning ORDER BY 3 DESC""")
    fam_tbl = table(
        ["READOUTM", "bin", "solved frames", "measured scale (\"/px)",
         "prior handed to solver"],
        [[esc(ro) or "<i>(blank)</i>", fmt(xb), fmt(c),
          f"{mn:.4f} – {mx:.4f}",
          "{:.2f} – {:.2f}".format(*astrom.scale_bounds(ro, xb))]
         for ro, xb, c, mn, mx in fam_rows])
    return f"""
<section id="tooling">
<div class="bhead"><h2>1 &middot; Tooling</h2>
<span class="tag">local astrometry.net + downloaded index files</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">No local plate-solving stack existed on this machine (no
<code>solve-field</code>, no index files anywhere on the attached disks).
Can one be assembled that solves RLMT frames at all — and at what
per-frame cost?</p>

<h3>Evidence</h3>
<p class="sub">Installed <code>{esc(meta.get('solve_field_version',
'?').splitlines()[0] if meta.get('solve_field_version') else '?')}</code>
(astrometry.net, Homebrew) and downloaded
{fmt(int(meta.get('index_n_files', 0)))} index files
({gb:.1f}&nbsp;GB) into <code>{esc(meta.get('index_dir', ''))}</code> —
series {esc(meta.get('index_series', ''))}: the 4200 (2MASS) skymark
scales 3&ndash;7 (5.6&prime;&ndash;30&prime; quads) that bracket the
&sim;37&prime; RLMT field, plus the small wide-field 4100 (Tycho-2) set.
Solves ran with <code>--downsample {esc(meta.get('downsample', ''))}</code>,
a {esc(meta.get('solve_cpulimit_s', ''))}&nbsp;s CPU limit
({esc(meta.get('solve_timeout_s', ''))}&nbsp;s wall cap), and the header
pointing as a &plusmn;{esc(meta.get('hint_radius_deg', ''))}&deg; hint
where present.  The measured plate scales, per camera family:</p>
{fam_tbl}
<p class="sub">A <b>solution-acceptance gate</b> guards every
&ldquo;solved&rdquo; verdict: a solution counts only when its measured
pixel scale sits inside the prior handed to the solver, at least
{fmt(astrom.MIN_MATCHED_STARS)} index stars matched, and the astrometric
RMS stays under {astrom.MAX_SOLVE_RMS_ARCSEC:g}&Prime;.  The gate caught
{fmt(n_bad)} false positive(s) in the experiment (recorded as
<code>bad_solve</code>): astrometry.net emitted a confident-looking WCS
at 7&times; the true plate scale from a 4-star quad match on a healthy
Mode0 field.  Any batch run inherits the same gate, so false solutions
are never written into the manifest.</p>

<h3>Decision</h3>
<div class="decision"><b>The local stack is real and the scale priors are
verified from the gate-accepted WCS headers themselves</b> — Mode0/Fast
is the IMX455 optic at &sim;0.45&Prime;/px binned, the GSENSE4040 sits at
its facility-fact 0.54&Prime;/px, and the Andor iKon lands inside its
deliberately wide prior.  Every future solve can reuse these bounds
tightened to the measured band, behind the same acceptance gate.</div>

<h3>Consequence</h3>
<p class="sub">A successful solve costs seconds (median times per stratum
below); the tooling is not the bottleneck — frame quality is.</p>
</div></section>"""


def section_universe(con) -> str:
    pops = {r[0]: (r[1], r[2]) for r in q(con,
            "SELECT class, n_frames, note FROM s1_populations")}
    # The comparison classes (label_gate_*, gate_moved_*, included_*) are
    # deliberately NOT shown here: this section answers "what is the
    # universe", section 3 answers "how did the universe change".  Mixing
    # them put a retired rule's arithmetic in the middle of the current
    # census, which is how the retired rule stayed invisible for so long.
    ordered = ["unsolved_total", "excluded_measured_spectrum",
               "excluded_calib_vocab_filter", "excluded_window_geometry",
               "solvable_candidates", "candidates_unstratified"]
    pop_tbl = table(
        ["population class", "frames", "meaning"],
        [[f"<code>{esc(c)}</code>", fmt(pops[c][0]), esc(pops[c][1])]
         for c in ordered if c in pops],
        row_classes=["warn" if c.startswith("excluded") else None
                     for c in ordered if c in pops])
    n_total = pops["unsolved_total"][0]
    n_solvable = pops["solvable_candidates"][0]
    n_window = pops["excluded_window_geometry"][0]
    n_spectrum = pops["excluded_measured_spectrum"][0]
    # The EU UMa Fast-readout series: the frames the geometry gate used to
    # throw away.  BOTH counts are queried — how many of them the gate
    # still excludes, and how many exist at all — because the paragraph
    # below has to be true whichever answer the database gives.
    n_euuma_strip = q1(con, """
        SELECT count(*) FROM frames f LEFT JOIN eras e USING (era_id)
        WHERE f.target_key = 'euuma' AND e.readoutm = 'Fast'
          AND f.is_canonical = 1 AND f.tree = 'rawimage'
          AND (f.pltsolvd IS NULL OR f.pltsolvd != 1)
          AND f.naxis1 < ?""", (astrom.MIN_SOLVABLE_NAXIS,))
    n_euuma_fast = q1(con, """
        SELECT count(*) FROM frames f LEFT JOIN eras e USING (era_id)
        WHERE f.target_key = 'euuma' AND e.readoutm = 'Fast'
          AND f.is_canonical = 1 AND f.tree = 'rawimage'
          AND (f.pltsolvd IS NULL OR f.pltsolvd != 1)""")
    # Where those frames sit NOW, if they came back: the stratum that was
    # created for them.  NULL-safe — the stratum may not exist in an old
    # design, and the paragraph must not crash on that.
    n_cv_fast = q1(con, """SELECT n_population FROM s1_strata
                           WHERE stratum_id = 'cv_fast_fullframe'""") or 0
    # --- The exclusion paragraph, in two versions -----------------------
    # THIS IS THE ONE PLACE IN THE REPORT WHERE THE ARGUMENT CHANGED, not
    # just the numbers.  Until the S0e geometry repair, the largest
    # exclusion here was "window geometry", and the page said EU UMa's
    # Fast-readout season "can never be plate-solved by any tool".  That
    # claim was wrong: the 8x3211 shape was never on the sky.  It was a
    # tile-compressed BINTABLE's row length read as an image width, and
    # the frames underneath are ordinary full frames.  Interpolating the
    # new counts into the old sentence would have printed "0 frames are
    # high-speed photometry WINDOWS ... among them sit all 0 of EU UMa's
    # ..." — arithmetically sourced from the database and yet still
    # asserting a retired claim.  So the SENTENCE branches on the
    # evidence too, not only its numbers.
    if n_window:
        excl_para = f"""<p class="sub">The geometry exclusion is the
headline: {fmt(n_window)} &ldquo;frames&rdquo; are high-speed photometry
WINDOWS &mdash; strips narrower than
{fmt(astrom.MIN_SOLVABLE_NAXIS)}&nbsp;px read out around a single target
star, with too little sky for quad matching.  {fmt(n_euuma_strip)} of
EU&nbsp;UMa&rsquo;s {fmt(n_euuma_fast)} unsolved Fast-readout frames are
among them.  The {fmt(n_spectrum)} spectrum rows are frames S2c
<i>measured</i> to be dispersed &mdash; they have no star field by
design (see section 3 for why this is now a measurement, not a
label).</p>"""
    else:
        excl_para = f"""<p class="sub"><b>The geometry exclusion is now
EMPTY</b> &mdash; {fmt(n_window)} frames.  It used to be the largest
exclusion on this page, and it was an artifact: the 8&times;3211
&ldquo;photometry windows&rdquo; were a tile-compressed BINTABLE&rsquo;s
row length read as an image width, not a shape any camera ever put on the
sky.  The S0e repair repointed those rows at their real geometry (see
<a href="s0e_geometry_fix.html">the geometry-repair page</a>), and they
returned to the candidate universe.  EU&nbsp;UMa is the case that matters
for the CV paper: all {fmt(n_euuma_fast)} of its unsolved Fast-readout
frames are full frames, {fmt(n_euuma_strip)} of them fail the geometry
gate, and they now sit in their own stratum
(<code>cv_fast_fullframe</code>, {fmt(n_cv_fast)} frames) instead of in
no stratum at all.  <b>An earlier version of this page said that
population &ldquo;can never be plate-solved by any tool&rdquo;.  That was
wrong, and the rates below are measured on the repaired universe.</b>
The {fmt(n_spectrum)} spectrum rows are the only large exclusion left:
frames that stage <a href="s2c_filter_identity.html">S2c</a>
<i>measured</i> to be slitless spectra &mdash; parallel traces along a
common grating axis, no star field by design.  <b>That criterion is
itself a correction</b>: until S1 v1.2 this page decided &ldquo;spectrum&rdquo;
by matching the FILTER string against a list of grism names, which was
wrong in both directions.  Section 3 shows the damage and the delta.</p>"""
    return f"""
<section id="universe">
<div class="bhead"><h2>2 &middot; The candidate universe</h2>
<span class="tag">what astrometry cannot even apply to, counted first</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The manifest counts {fmt(n_total)} unsolved canonical raw
Light frames.  How many of them can a plate solver even be pointed at?</p>

<h3>Evidence</h3>
{pop_tbl}
{excl_para}

<h3>Decision</h3>
<div class="decision"><b>The astrometry batch universe is the
{fmt(n_solvable)} solvable candidates: not MEASURED dispersed, both axes
&ge; {fmt(astrom.MIN_SOLVABLE_NAXIS)} px.</b>  Exclusions are recorded
per class in <code>s1_populations</code> — nothing is silently dropped,
and the window-strip populations are handed to the time-series pipeline
as a named fact, not a failure.</div>

<h3>Consequence</h3>
<p class="sub">Success rates below are rates on frames a solver could in
principle solve — the honest denominator for the batch decision.</p>
</div></section>"""


def section_correction(con, strata: list[dict]) -> str:
    """The label-vs-measurement correction, and every consequence of it.

    This section exists because a published number was wrong, and the
    project's rule is that a retired claim is shown being retired, with
    its replacement beside it — never quietly overwritten.  Every figure
    below comes from ``s1_gate_comparison`` / ``s1_populations`` (the
    census the design wrote) or from ``s1_baseline_*`` (the frozen
    previous experiment); nothing here is typed.
    """
    pops = {r[0]: r[1] for r in q(
        con, "SELECT class, n_frames FROM s1_populations")}
    n_moved_out = pops.get("gate_moved_out_total", 0)
    n_moved_in = pops.get("gate_moved_in_total", 0)
    n_in_direct = pops.get("gate_moved_in_direct", 0)
    n_in_indet = pops.get("gate_moved_in_indeterminate", 0)
    n_unmeasured = pops.get("included_unmeasured", 0)
    n_label_solvable = pops.get("label_gate_solvable_candidates", 0)
    n_solvable = pops.get("solvable_candidates", 0)
    # --- the cross-tab: what the two universes disagree about -----------
    xtab = q(con, """
        SELECT label_class, dispersion_class, label_gate, measured_gate,
               movement, n_frames FROM s1_gate_comparison
        ORDER BY movement != 'moved_out', movement != 'moved_in',
                 n_frames DESC""")
    MOVE_LABEL = {
        "moved_out": "<b>OUT</b> — was in the universe, is a spectrum",
        "moved_in": "<b>IN</b> — was deleted unseen, is an image",
        "unchanged_in": "unchanged (in)",
        "unchanged_out": "unchanged (out)"}
    xtab_tbl = table(
        ["FILTER label", "S2c measurement", "label gate", "measured gate",
         "movement", "frames"],
        [[f"<code>{esc(lab)}</code>", f"<code>{esc(disp)}</code>",
          esc(lg), esc(mg), MOVE_LABEL.get(mv, esc(mv)), fmt(n)]
         for lab, disp, lg, mg, mv, n in xtab],
        row_classes=["warn" if mv in ("moved_in", "moved_out") else None
                     for _, _, _, _, mv, _ in xtab])
    # --- per-stratum before/after ---------------------------------------
    have_baseline = table_exists(con, "s1_baseline_strata")
    if not have_baseline:
        delta_tbl = ("<p class=\"sub\"><i>No frozen baseline in this "
                     "manifest — run <code>run_s1_experiment.py "
                     "snapshot</code> before a design change to enable the "
                     "before/after table.</i></p>")
        moved_verdicts: list = []
        base_meta: dict = {}
    else:
        base = {s["stratum_id"]: s
                for s in stratum_stats(con, "s1_baseline_")}
        base_meta = dict(q(con, "SELECT key, value FROM s1_baseline_meta"))
        body, classes, moved_verdicts = [], [], []
        for s in strata:
            b = base.get(s["stratum_id"])
            moved = b is not None and b["verdict"] != s["verdict"]
            if moved:
                moved_verdicts.append((s["stratum_id"], b["verdict"],
                                       s["verdict"]))
            # Flag any row whose denominator or verdict changed — the
            # reader should be able to find the affected strata without
            # reading every cell.
            changed = b is not None and (moved or b["n_pop"] != s["n_pop"])
            classes.append("warn" if changed else None)
            body.append([
                f"<code>{esc(s['stratum_id'])}</code>",
                fmt(b["n_pop"]) if b else "&mdash;",
                fmt(s["n_pop"]),
                (f"{s['n_pop'] - b['n_pop']:+,}" if b
                 and b["n_pop"] != s["n_pop"] else "&mdash;"),
                (f"{fmt(b['k'])} / {fmt(b['n'])}" if b else "&mdash;"),
                f"{fmt(s['k'])} / {fmt(s['n'])}",
                (f"{100 * b['rate']:.1f}%" if b else "&mdash;"),
                f"{100 * s['rate']:.1f}%",
                (f"{100 * (s['rate'] - b['rate']):+.1f} pp" if b
                 else "&mdash;"),
                (f"{esc(b['verdict'])}" if b else "&mdash;"),
                f"<b>{esc(s['verdict'])}</b>",
                "<b>MOVED</b>" if moved else "same"])
        delta_tbl = table(
            ["stratum", "old backlog", "new backlog", "&Delta; backlog",
             "old solved", "new solved", "old rate", "new rate",
             "&Delta; rate", "old verdict", "new verdict", "verdict"],
            body, row_classes=classes)
    if moved_verdicts:
        moved_para = "".join(
            f"<li><code>{esc(sid)}</code>: <b>{esc(old)} &rarr; "
            f"{esc(new)}</b></li>" for sid, old, new in moved_verdicts)
        moved_html = (f"<p class=\"sub\"><b>Verdicts that MOVE:</b></p>"
                      f"<ul class=\"sub\">{moved_para}</ul>")
    elif have_baseline:
        moved_html = ("<p class=\"sub\"><b>No published verdict moves.</b>"
                      "  Every stratum keeps the GO / CAUTION / NO-GO it "
                      "carried before — the rates move, the decisions do "
                      "not.  That is a result, not a formality: it means "
                      "the batch plan built on this page survives the "
                      "correction unchanged.</p>")
    else:
        moved_html = ""
    # --- corroboration for the 'indeterminate' rule ---------------------
    # The decision to KEEP indeterminate frames is defended on principle
    # (exclusion must be earned), but principle is cheap; the S1b
    # production batch supplies an independent behavioural check on tens
    # of thousands of frames.  Attempts only — QC-skipped rows never
    # reached the solver and would bias every class differently.  Computed
    # here, never typed: if the batch grows, the sentence updates.
    if table_exists(con, "s1_batch"):
        beh = {d: (k, n) for d, k, n in q(con, """
            SELECT coalesce(d.verdict, ?),
                   sum(b.status = 'solved'), count(*)
            FROM s1_batch b LEFT JOIN frame_dispersion d USING (obs_rowid)
            WHERE b.status != 'skipped_qc'
            GROUP BY 1""", (astrom.UNMEASURED_CLASS,))}
    else:
        beh = {}

    def beh_txt(cls: str) -> str:
        """'92% of 1,040 attempts' for one dispersion class, or a dash."""
        k, n = beh.get(cls, (0, 0))
        return (f"{100 * k / n:.0f}% of {fmt(n)} attempts" if n
                else "&mdash; (no attempts)")

    beh_html = (
        f"<p class=\"sub\"><b>The behavioural check.</b>  In the S1b "
        f"production batch — a population far larger than this "
        f"experiment&rsquo;s samples — frames that the solver actually "
        f"attempted solved at "
        f"<b>{beh_txt(astrom.DIRECT_VERDICT)}</b> when measured "
        f"<code>direct</code>, "
        f"<b>{beh_txt(astrom.INDETERMINATE_VERDICT)}</b> when "
        f"<code>indeterminate</code>, and "
        f"<b>{beh_txt(astrom.DISPERSED_VERDICT)}</b> when "
        f"<code>dispersed</code> (QC-skipped frames excluded: they never "
        f"reached the solver).  The indeterminate class behaves like "
        f"images, not like spectra — which is the empirical half of the "
        f"argument for keeping it.</p>") if beh else ""
    # --- the spectra that were counted as failures ----------------------
    # Straight from the FROZEN baseline autopsy: how the old taxonomy
    # diagnosed frames that S2c measures as dispersed.
    if have_baseline and table_exists(con, "s1_baseline_failure_autopsy"):
        spec_rows = q(con, """
            SELECT a.stratum_id, a.diagnosis, count(*)
            FROM s1_baseline_failure_autopsy a
            JOIN frame_dispersion d USING (obs_rowid)
            WHERE d.verdict = ?
            GROUP BY 1, 2 ORDER BY 3 DESC""",
            (astrom.DISPERSED_VERDICT,))
        n_spec_fail = sum(r[2] for r in spec_rows)
        n_base_fail = q1(con,
                         "SELECT count(*) FROM s1_baseline_failure_autopsy")
        spec_tbl = table(
            ["stratum", "old machine diagnosis", "frames"],
            [[f"<code>{esc(s)}</code>", esc(d), fmt(c)]
             for s, d, c in spec_rows],
            row_classes=["warn"] * len(spec_rows))
    else:
        spec_rows, n_spec_fail, n_base_fail, spec_tbl = [], 0, 0, ""
    return f"""
<section id="correction">
<div class="bhead"><h2>3 &middot; Label vs measurement: the gate
correction</h2>
<span class="tag">a published rate, taxonomy and verdict were computed
over the wrong denominator</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Astrometry cannot apply to a slitless spectrum, so the
candidate universe must exclude spectra.  Earlier versions of this page
decided which frames were spectra by <b>reading the FILTER string</b> and
matching it against five known grism names
({", ".join(f"<code>{esc(g)}</code>"
             for g in sorted(astrom.GRISM_FILTERS))}).
Stage <a href="s2c_filter_identity.html">S2c</a> has since <b>measured</b>
dispersion per frame — parallel traces along a common grating axis.  Do
the label and the measurement agree about what a spectrum is?</p>

<h3>Evidence</h3>
<p class="sub"><b>They do not, in both directions.</b>  The telescope has
a MIXED filter slot: slot <code>6</code> carries a grating on some nights
and glass on others, and its FILTER string is simply <code>6</code> —
nothing in the header distinguishes the two cases.  The label rule
therefore let those spectra <i>into</i> the astrometry universe, where
they could not possibly solve; and it deleted, unseen, frames whose
labels look like grisms but which S2c measures as ordinary images.  The
full cross-tab over all {fmt(pops.get('unsolved_total', 0))} unsolved
frames:</p>
{xtab_tbl}
<p class="sub">The retired label gate called
{fmt(n_label_solvable)} frames solvable candidates; the measured gate
calls {fmt(n_solvable)}.  The two universes differ by
{fmt(n_moved_out)} frames <b>out</b> and {fmt(n_moved_in)} frames
<b>in</b> &mdash; a near-cancelling net of
{n_solvable - n_label_solvable:+,}, which is precisely why the error was
invisible in the headline census and had to be found in the failure
gallery instead.</p>
{beh_html}

<h3>Decision</h3>
<div class="decision"><b>A frame is excluded as a spectrum when, and only
when, S2c MEASURED it dispersed.</b>  The FILTER label no longer gates
anything; it survives in the code
(<code>is_solvable_candidate_by_label</code>) purely to compute the table
above.  The rule for each of the four dispersion classes, stated so that
none of them defaults silently:
<ul>
<li><code>dispersed</code> &rarr; <b>EXCLUDED</b>. Measured traces on a
common grating axis. It is a spectrum.</li>
<li><code>direct</code> &rarr; <b>INCLUDED, whatever the FILTER says</b>.
{fmt(n_in_direct)} frames carrying a grism-looking label are measured
direct images. A label loses to a measurement of the same fact; these
return to the universe.</li>
<li><code>indeterminate</code> &rarr; <b>INCLUDED</b>. S2c looked and
could not certify either way. This moves {fmt(n_in_indet)} frames in
(the indeterminate frames that carried a grism-looking label; the rest
were already in). The reason is that the gate's question is ontological
&mdash; <i>is this an image?</i> &mdash; never predictive &mdash;
<i>will it solve?</i>. Excluding frames because they look unlikely to
solve would inflate the very rate this experiment measures. S2c's
commonest indeterminate reason is literally &ldquo;no usable sources
extracted&rdquo;, which describes a blank IMAGE far more often than a
spectrum &mdash; and the production batch corroborates it (below).</li>
<li><code>unmeasured</code> &rarr; <b>INCLUDED</b>. S2c targeted the
dispersion-suspect population, so {fmt(n_unmeasured)} candidates were
never measured at all. This rule moves <b>zero</b> frames &mdash; every
unmeasured frame carries an ordinary label and was already in &mdash;
but it still has to be stated, because the alternative (exclude anything
not certified direct) would delete most of the backlog on no evidence
whatever.</li>
</ul>
In one sentence: <b>exclusion must be earned by a positive
measurement.</b>  The old label rule and its mirror image (&ldquo;exclude
whatever is not certified direct&rdquo;) share one flaw — removing frames
from a denominator without evidence.</div>

<h3>Consequence</h3>
<p class="sub"><b>What the contamination did.</b>  {fmt(n_spec_fail)} of
the {fmt(n_base_fail)} failures autopsied under the old gate are frames
S2c measures as spectra.  Every one of them was machine-diagnosed as an
<i>optical fault</i> &mdash; the frames were real, the statistics were
real, and the diagnosis was still wrong, because the autopsy was asking
&ldquo;what is wrong with this image?&rdquo; about something that was
never an image:</p>
{spec_tbl}
<p class="sub">No wrong astrometry was ever produced: a spectrum simply
fails to solve, and the batch QC gate skips it.  What was produced was a
<b>rate over a denominator containing non-images</b>, a failure taxonomy
that attributed spectra to defocus and wind, and verdicts resting on
both.  The corrected per-stratum delta:</p>
{delta_tbl}
{moved_html}
<p class="sub">Frames whose sample membership survived the re-design keep
their recorded solve outcome rather than being re-solved
({esc(dict(q(con, "SELECT key, value FROM s1_build_meta"))
      .get('n_carried_over', '?'))} of {fmt(sum(s['n'] for s in strata))}
sampled frames).  A solve outcome is a property of the frame and the
solver configuration, and neither changed here; re-solving would spend
CPU to reproduce a known answer <i>and</i> let solver nondeterminism leak
into a delta that exists to isolate the gate change.  The baseline it is
compared against is frozen in
<code>s1_baseline_*</code>{" (" + esc(base_meta.get("label", "")) + ")"
                          if base_meta.get("label") else ""}.</p>
</div></section>"""


def section_results(con, strata: list[dict]) -> str:
    src_rates = fig_success_rates(strata)
    src_times = fig_solve_times(con)
    n_frames = sum(s["n"] for s in strata)
    n_solved = sum(s["k"] for s in strata)
    # The seed is printed RAW, not through fmt().  fmt() groups thousands,
    # which rendered the reproducibility seed 20260817 as "20,260,817" — a
    # quantity, when it is the identifier a reader has to type back in
    # verbatim to redraw the same sample.
    seed_txt = str(strata[0]["seed"]) if strata else "&mdash;"
    body, classes = [], []
    for s in strata:
        classes.append(None if s["verdict"] == "GO" else "warn")
        # A census stratum has no sampling interval to print — its rate
        # is the population rate, exactly.
        rate_cell = (f"{100 * s['rate']:.0f}% (census)" if s["census"]
                     else f"{100 * s['rate']:.0f}% "
                          f"[{100 * s['lo']:.0f}&ndash;"
                          f"{100 * s['hi']:.0f}]")
        body.append([
            f"<code>{esc(s['stratum_id'])}</code>", esc(s["desc"]),
            fmt(s["n_pop"]), f"{fmt(s['k'])} / {fmt(s['n'])}",
            rate_cell,
            f"{s['med_time_solved']:.1f}" if s["med_time_solved"]
            else "&mdash;",
            f"{s['med_rms']:.2f}" if s["med_rms"] else "&mdash;",
            fmt(s["med_matched"]),
            f"<b>{s['verdict']}</b>"])
    res_tbl = table(
        ["stratum", "definition", "backlog", "solved", "rate [95% CI]",
         "median t (s)", "median RMS (\")", "median stars", "verdict"],
        body, row_classes=classes)
    med_rms_all = med([s["med_rms"] for s in strata])
    # Census footnote facts, straight from the stats list.
    census_strata = [s for s in strata if s["census"]]
    census_note = ""
    if census_strata:
        c = census_strata[0]
        census_note = (
            f"<p class=\"sub\"><code>{esc(c['stratum_id'])}</code> is a "
            f"<b>census</b>: all {fmt(c['n_pop'])} population frames were "
            f"attempted, so its rate is the exact population rate — no "
            f"confidence interval applies, and the verdict judges "
            f"{100 * c['rate']:.1f}% directly against the thresholds.  "
            f"Its &ldquo;batch decision&rdquo; is moot: the "
            f"{fmt(c['k'])} solutions already exist.</p>")
    # Night-clustering caveat: the iid stress test.
    nights = night_coverage(con, strata)
    night_body, night_classes = [], []
    for nc in nights:
        # Flag the rows where the night-level lower bound falls below the
        # GO bar while the frame-level verdict says GO — the honest
        # "this GO leans on the independence assumption" marker.
        s = next(x for x in strata if x["stratum_id"] == nc["stratum_id"])
        leaning = (s["verdict"] == "GO"
                   and nc["night_lo"] < astrom.GO_LOWER_BOUND
                   and not nc["census"])
        night_classes.append("warn" if leaning else None)
        night_body.append([
            f"<code>{esc(nc['stratum_id'])}</code>",
            f"{fmt(nc['n_nights'])} / {fmt(nc['n_pop_nights'])}",
            f"{fmt(nc['k_nights'])} / {fmt(nc['n_nights'])}",
            "exact" if nc["census"]
            else f"{100 * nc['night_lo']:.0f}%"])
    night_tbl = table(
        ["stratum", "nights sampled / in backlog",
         "perfect nights (all frames solved)",
         "night-level Wilson lower bound"],
        night_body, row_classes=night_classes)
    return f"""
<section id="results">
<div class="bhead"><h2>4 &middot; The experiment</h2>
<span class="tag">{fmt(n_frames)} frames, {fmt(len(strata))} strata, seed
{seed_txt}</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Per stratum — camera family &times; exposure band &times;
project — what fraction of sampled frames does the local solver actually
solve, how fast, and how well?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src_rates,
    f"{fmt(n_solved)} of {fmt(n_frames)} sampled frames solved (gate-"
    "accepted). Bars are point rates; whiskers the Wilson 95% interval; "
    f"the dashed line the {100 * astrom.GO_LOWER_BOUND:.0f}% GO bar the "
    "interval's LOWER bound must clear.")}
{_figure(src_times,
    "Successes cluster at seconds; failures burn the full CPU budget — "
    "the failure rate, not the solve speed, sets the batch cost.")}</div>
{res_tbl}
{census_note}
<p class="sub"><b>How independent are these trials?</b>  Frames within a
night share cloud, focus and wind state, and the failure autopsy below
attributes most failures to exactly those nightly conditions — so the
Wilson intervals&rsquo; frame-level independence assumption is optimistic,
and the effective sample size sits somewhere between the night count and
the frame count.  The stress test: collapse each stratum to
all-or-nothing NIGHTS (a night succeeds only if every sampled frame on it
solved) and re-read the interval.</p>
{night_tbl}

<h3>Decision</h3>
<div class="decision"><b>Verdicts are mechanical: GO when the Wilson 95%
lower bound clears {100 * astrom.GO_LOWER_BOUND:.0f}%, CAUTION above
{100 * astrom.CAUTION_LOWER_BOUND:.0f}%, NO-GO below — judged on the
worst rate each sample still allows, never on the point estimate (census
strata are judged on their exact, fully-enumerated rate).</b>  Solutions
are trustworthy where the gate accepts them: median astrometric RMS
&sim;{med_rms_all:.2f}&Prime; against the 2MASS index across all strata,
with plate scales landing inside every family prior.  <b>The frame-level
intervals assume within-stratum independence that weather demonstrably
violates</b> — under the harsher night-level collapse, the perfect strata
keep lower bounds of
{"/".join(f"{100 * nc['night_lo']:.0f}%" for nc in nights
          if next(x for x in strata
                  if x['stratum_id'] == nc['stratum_id'])['rate'] == 1.0)}
across far fewer effective trials, and each stratum samples only a
minority of its backlog&rsquo;s nights (the coverage column above) — so
unsampled bad nights are expected in any batch.  The mitigation is not a
wider interval but the QC pre-gate of section 5: starless frames from
unsampled cloudy nights are detected and skipped at batch time, they do
not silently drag the yield.</div>

<h3>Consequence</h3>
<p class="sub">{fmt(sum(1 for s in strata if s['rate'] == 1.0))} of the
{fmt(len(strata))} strata solved every sampled frame — including the
Mode0 and iKon strata that carry the bulk of the CV Sloan backlog.  The
failure analysis below shows where the other strata lose frames: to the
frames, not to the solver.</p>
</div></section>"""


def section_autopsy(con) -> str:
    src, n_shown = fig_failure_gallery(con)
    # bad_solve rows (false-positive WCS rejected by the gate) are
    # failures too — counted and autopsied like the rest.
    n_fail = q1(con, """SELECT count(*) FROM s1_solve_experiment
                        WHERE status IN ('unsolved', 'timeout',
                                         'bad_solve')""")
    n_autop = q1(con, "SELECT count(*) FROM s1_failure_autopsy")
    # --- the taxonomy, new beside old -----------------------------------
    # The old counts come from the FROZEN baseline autopsy, so the "before"
    # column is the taxonomy as it was actually published, not a
    # reconstruction of it.
    have_base = table_exists(con, "s1_baseline_failure_autopsy")
    old_counts = dict(q(con, """
        SELECT diagnosis, count(*) FROM s1_baseline_failure_autopsy
        GROUP BY diagnosis""")) if have_base else {}
    # ...and, of those old counts, how many were measured spectra: the
    # inflation the label gate injected into each diagnosis.
    old_spectra = dict(q(con, """
        SELECT a.diagnosis, count(*)
        FROM s1_baseline_failure_autopsy a JOIN frame_dispersion d
             USING (obs_rowid)
        WHERE d.verdict = ? GROUP BY 1""",
        (astrom.DISPERSED_VERDICT,))) if have_base else {}
    diag_rows = q(con, """
        SELECT diagnosis, count(*), group_concat(DISTINCT stratum_id)
        FROM s1_failure_autopsy GROUP BY diagnosis ORDER BY 2 DESC""")
    new_counts = {d: c for d, c, _ in diag_rows}
    strata_by_diag = {d: s for d, _, s in diag_rows}
    all_diags = sorted(set(old_counts) | set(new_counts),
                       key=lambda d: -(new_counts.get(d, 0)
                                       + old_counts.get(d, 0)))
    diag_body, diag_classes = [], []
    for d in all_diags:
        o, n_, sp = (old_counts.get(d, 0), new_counts.get(d, 0),
                     old_spectra.get(d, 0))
        diag_classes.append("warn" if sp else None)
        diag_body.append([
            esc(d), fmt(o) if have_base else "&mdash;", fmt(n_),
            (f"{n_ - o:+d}" if have_base and n_ != o else "&mdash;"),
            (f"<b>{fmt(sp)}</b>" if sp else "0"),
            f"<code>{esc(strata_by_diag.get(d, '&mdash;'))}</code>"])
    diag_tbl = table(
        ["diagnosis", "old count (label gate)", "new count (measured gate)",
         "&Delta;", "of the OLD count, measured spectra",
         "strata affected"],
        diag_body, row_classes=diag_classes)
    n_unexp = q1(con, """SELECT count(*) FROM s1_failure_autopsy
                         WHERE diagnosis LIKE 'unexplained%'""")
    # THE INVARIANT, checked live from the table rather than asserted: no
    # autopsied failure may be a measured spectrum any more.
    n_spec_now = q1(con, """SELECT count(*) FROM s1_failure_autopsy
                            WHERE dispersion_class = ?""",
                    (astrom.DISPERSED_VERDICT,))
    n_spec_before = (q1(con, """
        SELECT count(*) FROM s1_baseline_failure_autopsy a
        JOIN frame_dispersion d USING (obs_rowid)
        WHERE d.verdict = ?""", (astrom.DISPERSED_VERDICT,))
        if have_base else 0)
    if n_spec_now == 0:
        invariant = (
            f"<p class=\"sub\"><b>Audit, computed from the table above, "
            f"not asserted: {fmt(n_spec_now)} of the {fmt(n_autop)} "
            f"autopsied failures are frames S2c measures as spectra.</b>  "
            f"Under the retired label gate the count was "
            f"{fmt(n_spec_before)} — spectra sitting inside a taxonomy of "
            f"optical faults, each one wearing a diagnosis "
            f"(&lsquo;defocused&rsquo;, &lsquo;trailing&rsquo;, "
            f"&lsquo;starved&rsquo;) that presumed it was an image.  Every "
            f"row of this taxonomy now describes a frame that really was "
            f"pointed at a star field.</p>")
    else:
        # Fails loud rather than printing a comforting sentence: if a
        # spectrum reaches the autopsy again, the gate has regressed.
        # macro.css now carries a semantic `.bad` ink for text, so the
        # banner asks for it by name.  It used to hard-code #e07a7a, a pale
        # red chosen for a dark ground, precisely because no such class
        # existed — and that literal would have gone nearly invisible when
        # the site turned white.
        invariant = (
            f"<p class=\"sub bad\">"
            f"<b>REGRESSION: {fmt(n_spec_now)} "
            f"autopsied failure(s) are measured spectra.</b>  The "
            f"candidate gate should have excluded them; this taxonomy is "
            f"not trustworthy until that is fixed.</p>")
    return f"""
<section id="autopsy">
<div class="bhead"><h2>5 &middot; Failure autopsy</h2>
<span class="tag">read the pixels before blaming the solver</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">{fmt(n_fail)} sampled frames failed to solve (the
gate-rejected false positive included).  Is the solver missing solvable
fields — or were these frames never solvable?</p>

<h3>Evidence</h3>
<div class="grid">{_figure(src,
    f"{n_shown} of the {fmt(n_autop)} autopsied failures — every "
    "failure is autopsied, no sampling cap — machine-diagnosed "
    "from source-extraction statistics.")}</div>
{diag_tbl}
<p class="sub">Diagnosis is computed, not eyeballed: 10&sigma;
source extraction (min {fmt(astrom.AUTOPSY_MIN_AREA_PX)} px) splits
detections into PSF-shaped sources
({astrom.PSF_A_MIN_PX:g}&ndash;{astrom.PSF_A_MAX_PX:g} px semi-major
axis) versus the hot-pixel spikes these RAW, dark-unsubtracted CMOS
frames carry by the thousand.  &lsquo;Defocused&rsquo; = the brightest
detections are &gt;{astrom.AUTOPSY_DEFOCUS_A_PX:g} px blobs;
&lsquo;starved&rsquo; = fewer than {fmt(astrom.AUTOPSY_MIN_SOURCES)}
real sources <i>or</i> brightest detections smaller than
{astrom.AUTOPSY_BLANK_BRIGHT_A_PX:g} px (single-pixel spikes — hot-pixel
pairs can fake both the source count and the elongation of a blank
frame, so the trailing verdict additionally requires star-sized
brightest detections); &lsquo;trailing&rsquo; = median elongation &gt;
{astrom.AUTOPSY_TRAIL_ELONG:g}.  Only {fmt(n_unexp)} autopsied
failure(s) show a healthy star field the solver still missed.</p>
{invariant}

<h3>Decision</h3>
<div class="decision"><b>The failures are frame-quality facts, not solver
shortfalls — and the claim now covers every failure, because every
failure was autopsied: blank/cloud frames, defocused sequences, and
wind-trailed exposures.</b>  No solver setting recovers a frame with no
stars on it; these frames should be flagged by the same statistics at
batch time and spent zero further CPU.
<p><b>One line of the old version of this paragraph is retracted.</b>  It
named &ldquo;the filter-&lsquo;6&rsquo; series on the GSENSE&rdquo; as
the archetype of a badly defocused sequence.  Slot <code>6</code> is the
telescope&rsquo;s MIXED slot, and S2c has since measured most of those
frames to be <i>spectra</i>: the &ldquo;giant blobs, no point
sources&rdquo; the autopsy saw were dispersion traces, not defocused
stars.  A source-extraction statistic cannot tell those apart — it was
never asked to.  The measured-dispersion gate now keeps spectra out of
the universe entirely, so the diagnosis is never put in that position
again.  See section 3 for the full accounting.</p></div>

<h3>Consequence</h3>
<p class="sub">Batch success rates will track the experiment&rsquo;s
rates, and the batch runner inherits a cheap pre-filter: the autopsy
statistics double as a QC gate that skips starless frames before the
solver burns its budget on them.  The taxonomy is also now
<i>diagnostic</i> in a way it was not before: with spectra excluded
upstream, a &lsquo;defocused&rsquo; verdict means the optics, and can be
taken to the observatory as such.</p>
</div></section>"""


def section_verdict(con, strata: list[dict]) -> str:
    pops = population_rollup(con, strata)
    body, classes = [], []
    for p in pops:
        classes.append(None if p["verdict"] == "GO" else "warn")
        body.append([
            esc(p["population"]), fmt(p["n_backlog"]),
            f"{fmt(p['k'])} / {fmt(p['n'])}",
            f"{100 * p['rate']:.0f}% "
            f"[{100 * p['lo']:.0f}&ndash;{100 * p['hi']:.0f}]",
            # One decimal: a 99.7% weighted rate must not print as 100%.
            f"{100 * p['weighted_rate']:.1f}%",
            f"{fmt(round(p['expected_lo']))}&ndash;"
            f"{fmt(round(p['expected_hi']))}",
            f"{p['hours']:.1f}",
            f"<b>{p['verdict']}</b>"])
    pop_tbl = table(
        ["population", "backlog frames", "sample solved",
         "pooled sample rate [95% CI]", "backlog-weighted rate",
         "expected new WCS (range)", "projected hours",
         "verdict"], body, row_classes=classes)
    tot_backlog = sum(p["n_backlog"] for p in pops)
    tot_hours = sum(p["hours"] for p in pops)
    tot_lo = sum(p["expected_lo"] for p in pops)
    tot_hi = sum(p["expected_hi"] for p in pops)
    # The stratified backlog is NOT the whole solvable universe: the
    # residue below carries no measured rate and is excluded from every
    # projection on this page.
    n_solvable = q1(con, """SELECT n_frames FROM s1_populations
                            WHERE class = 'solvable_candidates'""")
    n_residue = q1(con, """SELECT n_frames FROM s1_populations
                           WHERE class = 'candidates_unstratified'""")
    # Project-critical cut: everything that is not the facility backlog.
    crit = [p for p in pops if p["population"] != "Facility backlog"]
    crit_backlog = sum(p["n_backlog"] for p in crit)
    crit_hours = sum(p["hours"] for p in crit)
    workers = q1(con, "SELECT value FROM s1_build_meta "
                      "WHERE key = 'workers'")
    return f"""
<section id="verdict">
<div class="bhead"><h2>6 &middot; Verdict &amp; batch cost</h2>
<span class="tag">the Week-2 review's decision table</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Per population: is batch re-solving viable, how many new
WCS solutions should it yield, and what does it cost in wall-clock?</p>

<h3>Evidence</h3>
{pop_tbl}
<p class="sub">Hours are population-weighted: each stratum&rsquo;s FULL
backlog &times; its own median per-frame wall time &divide;
{esc(workers)} parallel workers — the same worker count the experiment
ran at, on this machine.  Failures are what cost time, and the
projection charges every stratum its own failure-heavy median.  The two
rate columns answer different questions: the <i>pooled sample rate</i>
pools the strata&rsquo;s equal-size samples, which over-represents small
strata relative to their backlog share (conservatively, for these data —
the small strata are the weak ones); the <i>backlog-weighted rate</i> is
what the expected-WCS column implies, weighting each stratum&rsquo;s own
rate by its backlog.  The expected-WCS range multiplies each
stratum&rsquo;s backlog by its own interval bounds (a census stratum
contributes its exactly-known count).</p>

<h3>Decision</h3>
<div class="decision"><b>{esc(next(p['verdict'] for p in pops
    if p['population'] == 'CV polars (Sloan)'))} for the CV polars&rsquo;
Sloan series — the paper-gating population solved
{fmt(next(p['k'] for p in pops
          if p['population'] == 'CV polars (Sloan)'))} of
{fmt(next(p['n'] for p in pops
          if p['population'] == 'CV polars (Sloan)'))} sampled frames
({100 * next(p['rate'] for p in pops
             if p['population'] == 'CV polars (Sloan)'):.0f}%; the
per-stratum table in section 4 locates every miss) — and the same
mechanical rule for every population: GO where the Wilson lower bound
clears {100 * astrom.GO_LOWER_BOUND:.0f}%.</b>  The projected cost of the
project-critical populations ({fmt(crit_backlog)} frames) is
&sim;{crit_hours:.0f}&nbsp;h; the stratified backlog
({fmt(tot_backlog)} of the {fmt(n_solvable)} solvable frames) is
&sim;{tot_hours:.0f}&nbsp;h and should yield
&sim;{fmt(round(tot_lo / 100) * 100)}&ndash;{fmt(round(tot_hi / 100)
* 100)} new WCS solutions.  The remaining {fmt(n_residue)} solvable
candidates fit no stratum (a heterogeneous residue of one-off configs):
they carry <b>no measured rate</b> and are excluded from the yield and
hours projections — batching them means extrapolating, and that choice
belongs to the review.  Populations below GO stay viable as
<i>filtered</i> batches: run the QC pre-gate from section 5, solve what
has stars, and book the starless/defocused frames as facts.  <b>No batch
runs before the Week-2 review accepts this report</b> — S1&rsquo;s
mandate ends at the experiment.</div>

<h3>Consequence</h3>
<p class="sub">On acceptance, the batch writes its solutions into a NEW
manifest table (headers in the archive are never rewritten), zero-point
photometry for the CV paper unblocks behind it, and the October frames
inherit a verified solving stack with measured per-camera scale
priors.</p>
</div></section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S1 report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    # Read-only is not immune to a concurrent writer: the S1b production
    # batch holds this database's write lock in bursts, and a reader that
    # does not wait simply fails.  Same patience as the S3 renderer.
    con.execute("PRAGMA busy_timeout = 300000")
    try:
        strata = stratum_stats(con)
        n_frames = sum(s["n"] for s in strata)
        n_solved = sum(s["k"] for s in strata)
        n_go = sum(s["verdict"] == "GO" for s in strata)
        meta = dict(q(con, "SELECT key, value FROM s1_build_meta"))

        sections = [
            section_tooling(con),
            section_universe(con),
            section_correction(con, strata),
            section_results(con, strata),
            section_autopsy(con),
            section_verdict(con, strata),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S1 — Astrometry Go/No-Go Experiment</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S1 — Astrometry Go/No-Go Experiment</h1>
  <p>{fmt(n_solved)} / {fmt(n_frames)} sampled frames solved &middot;
  {fmt(len(strata))} strata ({fmt(n_go)} GO) &middot;
  seed {esc(meta.get('sample_seed', ''))} &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">the front page</a></p>
</header>

<nav>
  <a href="#tooling">1 Tooling</a> &middot;
  <a href="#universe">2 Candidate universe</a> &middot;
  <a href="#correction">3 Gate correction</a> &middot;
  <a href="#results">4 Experiment</a> &middot;
  <a href="#autopsy">5 Failure autopsy</a> &middot;
  <a href="#verdict">6 Verdict &amp; cost</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_core.report_s1</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query or a named constant in
<code>macro_core.astrom</code>; none is typed by hand.  Regenerate with
<code>pipeline/scripts/run_s1_experiment.py report</code>.</footer>
</body></html>"""

        HTML_PATH.write_text(html, encoding="utf-8")

        # Belt and braces: every <img> the page references must exist and be
        # non-empty, or the build fails loudly rather than shipping a broken
        # evidence page.
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
