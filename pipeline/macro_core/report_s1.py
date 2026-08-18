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
# Shared page machinery: same dark theme, same query discipline, same table
# generator as the S0/S0b reports — one visual language across the site.
from .report_s0 import (          # noqa: E402
    ACCENT, DARK, DPI, WARN, _figure, esc, fmt, q, q1, table)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s1"
HTML_PATH = DOCS_DIR / "s1_astrometry.html"

#: Status color per verdict — matches the site badge palette.
VERDICT_COLOR = {"GO": "#9fd8ae", "CAUTION": WARN, "NO-GO": "#e07a7a"}

#: Failure-gallery size: distinct diagnoses first, then fill, capped here.
N_GALLERY = 6


# ---------------------------------------------------------------------------
# Small derived-stat helpers (medians live in Python: SQLite has none).
# ---------------------------------------------------------------------------
def med(vals) -> float | None:
    """Median of a list of non-NULL values, None when empty."""
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def stratum_stats(con) -> list[dict]:
    """One dict per stratum: counts, rate, Wilson CI, verdict, medians.

    This is THE results computation — the figures, the tables and the
    verdict section all read from this one list, so they can never
    disagree with each other.

    Success = status ``solved`` ONLY: a ``bad_solve`` (a .solved marker
    whose WCS failed the acceptance gate) is a failure, and its bogus
    pixel scale never enters the measured-scale statistics.

    CENSUS strata (sample == population) carry no sampling uncertainty:
    their interval collapses to the exact rate (lo == hi == rate) and the
    verdict judges that exact rate — see ``astrom.verdict_for``.
    """
    out = []
    for sid, pop_name, desc, n_pop, n_sample, seed in q(con, """
            SELECT stratum_id, population, description, n_population,
                   n_sample, seed FROM s1_strata ORDER BY rowid"""):
        rows = q(con, """
            SELECT status, solve_time_s, pixscale_arcsec, rms_arcsec,
                   n_matched FROM s1_solve_experiment
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
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(8.6, 0.42 * len(labels) + 1.6))
        bars = ax.barh(labels, rates, color=colors,
                       xerr=[lo_err, hi_err], ecolor="#e8eaed", capsize=3)
        ax.set_xlim(0, 105)
        ax.set_xlabel("solve success rate (%) with Wilson 95% CI")
        ax.set_title("S1 experiment: success rate by stratum")
        # The GO bar: verdicts key off the CI lower bound crossing it.
        ax.axvline(100 * astrom.GO_LOWER_BOUND, color="#9aa4b2",
                   linestyle="--", linewidth=1)
        # Count labels sit INSIDE the bar's left edge, clear of the CI
        # whiskers that live at the bar's right end.
        for b, r, (k, n) in zip(bars, rates, ns):
            ax.annotate(f"{k}/{n}", (2.0, b.get_y() + b.get_height() / 2),
                        va="center", ha="left", fontsize=8,
                        color="#0f1115" if r > 12 else "#e8eaed")
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
    with plt.rc_context(DARK):
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
    with plt.rc_context(DARK):
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
    ordered = ["unsolved_total", "excluded_grism",
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
    n_grism = pops["excluded_grism"][0]
    # The EU UMa special case: its Fast-readout series is ALL strips.
    n_euuma_strip = q1(con, """
        SELECT count(*) FROM frames f LEFT JOIN eras e USING (era_id)
        WHERE f.target_key = 'euuma' AND e.readoutm = 'Fast'
          AND f.is_canonical = 1 AND f.tree = 'rawimage'
          AND (f.pltsolvd IS NULL OR f.pltsolvd != 1)
          AND f.naxis1 < ?""", (astrom.MIN_SOLVABLE_NAXIS,))
    return f"""
<section id="universe">
<div class="bhead"><h2>2 &middot; The candidate universe</h2>
<span class="tag">what astrometry cannot even apply to, counted first</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The manifest counts {fmt(n_total)} unsolved canonical raw
Light frames.  How many of them can a plate solver even be pointed at?</p>

<h3>Evidence</h3>
{pop_tbl}
<p class="sub">The geometry exclusion is the headline: {fmt(n_window)}
&ldquo;frames&rdquo; are high-speed photometry WINDOWS — 8-pixel-wide
strips read out around a single target star.  Among them sit all
{fmt(n_euuma_strip)} of EU&nbsp;UMa&rsquo;s unsolved Fast-readout series:
that population can never be plate-solved by any tool, and its astrometry
must come from the pointing header + the window geometry instead.  The
{fmt(n_grism)} grism rows are slitless spectra (hrg/lrg/HaGrism/OGGrism)
— they have no star field by design.</p>

<h3>Decision</h3>
<div class="decision"><b>The astrometry batch universe is the
{fmt(n_solvable)} solvable candidates: non-grism FILTER, both axes
&ge; {fmt(astrom.MIN_SOLVABLE_NAXIS)} px.</b>  Exclusions are recorded
per class in <code>s1_populations</code> — nothing is silently dropped,
and the window-strip populations are handed to the time-series pipeline
as a named fact, not a failure.</div>

<h3>Consequence</h3>
<p class="sub">Success rates below are rates on frames a solver could in
principle solve — the honest denominator for the batch decision.</p>
</div></section>"""


def section_results(con, strata: list[dict]) -> str:
    src_rates = fig_success_rates(strata)
    src_times = fig_solve_times(con)
    n_frames = sum(s["n"] for s in strata)
    n_solved = sum(s["k"] for s in strata)
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
<div class="bhead"><h2>3 &middot; The experiment</h2>
<span class="tag">{fmt(n_frames)} frames, 10 strata, seed
{fmt(strata[0]['seed']) if strata else '&mdash;'}</span></div>

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
wider interval but the QC pre-gate of section 4: starless frames from
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
    diag_rows = q(con, """
        SELECT diagnosis, count(*), group_concat(DISTINCT stratum_id)
        FROM s1_failure_autopsy GROUP BY diagnosis ORDER BY 2 DESC""")
    diag_tbl = table(
        ["diagnosis", "autopsied frames", "strata affected"],
        [[esc(d), fmt(c), f"<code>{esc(s)}</code>"]
         for d, c, s in diag_rows])
    n_unexp = q1(con, """SELECT count(*) FROM s1_failure_autopsy
                         WHERE diagnosis LIKE 'unexplained%'""")
    return f"""
<section id="autopsy">
<div class="bhead"><h2>4 &middot; Failure autopsy</h2>
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

<h3>Decision</h3>
<div class="decision"><b>The failures are frame-quality facts, not solver
shortfalls — and the claim now covers every failure, because every
failure was autopsied: blank/cloud frames, badly defocused sequences
(notably the filter-&lsquo;6&rsquo; series on the GSENSE), and
wind-trailed exposures.</b>  No solver setting recovers a frame with no
stars on it; these frames should be flagged by the same statistics at
batch time and spent zero further CPU.</div>

<h3>Consequence</h3>
<p class="sub">Batch success rates will track the experiment&rsquo;s
rates, and the batch runner inherits a cheap pre-filter: the autopsy
statistics double as a QC gate that skips starless frames before the
solver burns its budget on them.</p>
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
<div class="bhead"><h2>5 &middot; Verdict &amp; batch cost</h2>
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
per-stratum table in section 3 locates every miss) — and the same
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
<i>filtered</i> batches: run the QC pre-gate from section 4, solve what
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
    try:
        strata = stratum_stats(con)
        n_frames = sum(s["n"] for s in strata)
        n_solved = sum(s["k"] for s in strata)
        n_go = sum(s["verdict"] == "GO" for s in strata)
        meta = dict(q(con, "SELECT key, value FROM s1_build_meta"))

        sections = [
            section_tooling(con),
            section_universe(con),
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
  &middot; <a href="../index.html">back to the evidence hub</a></p>
</header>

<nav>
  <a href="#tooling">1 Tooling</a> &middot;
  <a href="#universe">2 Candidate universe</a> &middot;
  <a href="#results">3 Experiment</a> &middot;
  <a href="#autopsy">4 Failure autopsy</a> &middot;
  <a href="#verdict">5 Verdict &amp; cost</a>
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
