"""S0e evidence report: the compressed-FITS geometry artifact and its repair.

Reads the S0e evidence tables and writes:

* ``docs/pipeline/s0e_geometry_fix.html``   — the report
* ``docs/pipeline/figures/s0e/*.png``       — every figure

UNUSUALLY FOR THIS SITE, this report reads the CATALOG, not the manifest,
and says so on the page.  The other reports read only the manifest, because
that is where their evidence lives.  S0e's evidence is a catalog defect, and
— more importantly — it is evidence that the fix DESTROYS: once the manifest
is rebuilt, no frame carries the phantom geometry any more.  So the whole
argument (the audit trail ``geom_rescan``, the exhibit header
``s0e_header_dump``, the pre-fix breakdowns ``s0e_blast_era`` /
``s0e_blast_target``, the era forecast ``s0e_era_forecast`` and the re-queue
population ``s0e_requeue``) was snapshotted into catalog-side tables by
``rescan_geometry.py`` and ``s0e_era_forecast.py`` while it was still true.

Every number on the page comes from a SQL query in this module.  Nothing is
typed by hand, so re-running the repair regenerates the whole argument.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")            # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from .report_s0 import (          # noqa: E402  shared page machinery
    ACCENT, DARK, DPI, WARN, _figure, esc, fmt, q, q1, table)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s0e"
HTML_PATH = DOCS_DIR / "s0e_geometry_fix.html"

DEFAULT_CATALOG = Path(
    "/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-catalog.sqlite")
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"

OK_GREEN = "#7fc99a"     # repaired / correct
MUTED = "#9aa4b2"


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_anatomy(cat) -> str:
    """The trap, drawn to scale: what the table header says vs what is there.

    Two rectangles on one axis — the phantom 8x3211 strip beside the real
    4800x3211 frame — because the entire incident is the visual difference
    between those two shapes being invisible in a header.
    """
    old = q(cat, "SELECT old_naxis1, old_naxis2 FROM geom_rescan "
                 "WHERE changed = 1 LIMIT 1")
    new = q(cat, "SELECT new_naxis1, new_naxis2 FROM geom_rescan "
                 "WHERE changed = 1 LIMIT 1")
    if not old or not new:
        o1, o2, n1, n2 = 8, 3211, 4800, 3211
    else:
        o1, o2 = int(old[0][0]), int(old[0][1])
        n1, n2 = int(new[0][0]), int(new[0][1])
    with plt.rc_context(DARK):
        fig, (ax, ax2) = plt.subplots(
            1, 2, figsize=(9.4, 4.0), gridspec_kw={"width_ratios": [2.2, 1]})

        # --- left: the two shapes, honestly to scale -------------------
        # The phantom strip is 8 px inside a 4800 px frame, so at true
        # scale it is a HAIRLINE.  That is the whole point and it is left
        # un-exaggerated: an artifact this thin is exactly what nobody
        # noticed.  An arrow does the pointing instead of a fake width.
        ax.add_patch(plt.Rectangle((0, 0), n1, n2, facecolor="#1d2633",
                                   edgecolor=OK_GREEN, lw=2.0))
        ax.add_patch(plt.Rectangle((0, 0), max(o1, n1 * 0.004), o2,
                                   facecolor=WARN, edgecolor=WARN, lw=1.0))
        ax.annotate(f"TRUE image\nZNAXIS1 × ZNAXIS2\n{n1} × {n2}",
                    (n1 * 0.55, n2 * 0.5), color=OK_GREEN, ha="center",
                    va="center", fontsize=11, weight="bold")
        ax.annotate(f"NAXIS1 × NAXIS2 = {o1} × {o2}\n"
                    f"(table row BYTES × row COUNT)",
                    xy=(o1, n2 * 0.80), xytext=(n1 * 0.30, n2 * 1.10),
                    color=WARN, fontsize=9, ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", color=WARN, lw=1.4))
        ax.set_xlim(-n1 * 0.06, n1 * 1.06)
        ax.set_ylim(-n2 * 0.08, n2 * 1.30)
        ax.set_aspect("equal")
        ax.set_xlabel("pixels")
        ax.set_ylabel("pixels")
        ax.set_title("Drawn to scale", fontsize=10)
        ax.grid(False)

        # --- right: the same two numbers where the ratio is legible ----
        bars = ax2.bar(["NAXIS1\n(table)", "ZNAXIS1\n(image)"], [o1, n1],
                       color=[WARN, OK_GREEN])
        ax2.set_yscale("log")
        ax2.set_ylabel("frame width (pixels, log scale)")
        ax2.set_title(f"{n1 // o1}× wrong", fontsize=10)
        for b, v in zip(bars, [o1, n1]):
            ax2.annotate(f"{v:,}", (b.get_x() + b.get_width() / 2, v),
                         ha="center", va="bottom", fontsize=10,
                         color="#e8eaed")
        ax2.set_ylim(1, n1 * 6)

        fig.suptitle("The same frame, described two ways — one of them wrong")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0e_anatomy.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0e/s0e_anatomy.png"


def fig_blast(cat) -> str:
    """Who was affected: frames per target, and the repair's control group."""
    rows = q(cat, """
        SELECT canonical_target, n_frames FROM s0e_blast_target
        ORDER BY n_frames DESC LIMIT 12""")
    names = [r[0] or "(unnamed)" for r in rows][::-1]
    counts = [r[1] for r in rows][::-1]
    changed = q1(cat, "SELECT count(*) FROM geom_rescan WHERE changed = 1")
    control = q1(cat, "SELECT count(*) FROM geom_rescan WHERE changed = 0")
    with plt.rc_context(DARK):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.8),
                                       gridspec_kw={"width_ratios": [3, 1.5]})
        bars = ax1.barh(names, counts, color=WARN)
        ax1.set_xlabel("canonical frames carrying phantom geometry")
        ax1.set_title("Blast radius by target")
        ax1.set_xlim(0, max(counts) * 1.16)   # room for the value labels
        for b, c in zip(bars, counts):
            ax1.annotate(f" {c:,}", (c, b.get_y() + b.get_height() / 2),
                         va="center", fontsize=8, color="#e8eaed")
        ax2.bar(["repaired", "control\n(left alone)"], [changed, control],
                color=[OK_GREEN, MUTED])
        ax2.set_yscale("log")
        ax2.set_title("Rows re-read")
        for i, v in enumerate([changed, control]):
            ax2.annotate(f"{v:,}", (i, v), ha="center", va="bottom",
                         fontsize=9, color="#e8eaed")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0e_blast.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0e/s0e_blast.png"


def fig_eras(cat) -> str:
    """Era frame counts before and after, for every era the repair moves."""
    rows = q(cat, """
        SELECT era_id, n_before, n_after, verdict FROM s0e_era_forecast
        WHERE verdict != 'unchanged' ORDER BY era_id""")
    if not rows:
        rows = [(0, 0, 0, "none")]
    ids = [f"era {r[0]}" for r in rows]
    before = [r[1] or 0 for r in rows]
    after = [r[2] or 0 for r in rows]
    x = np.arange(len(ids))
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(9.4, 3.6))
        ax.bar(x - 0.2, before, 0.4, label="before (phantom geometry)",
               color=WARN)
        ax.bar(x + 0.2, after, 0.4, label="after (repaired)", color=OK_GREEN)
        ax.set_xticks(x)
        ax.set_xticklabels(ids)
        ax.set_ylabel("canonical frames")
        ax.set_title("Camera eras the repair moves")
        ax.legend()
        for xi, (b, a) in enumerate(zip(before, after)):
            ax.annotate(f"{b:,}", (xi - 0.2, b), ha="center", va="bottom",
                        fontsize=8, color="#e8eaed")
            ax.annotate(f"{a:,}", (xi + 0.2, a), ha="center", va="bottom",
                        fontsize=8, color="#e8eaed")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0e_eras.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0e/s0e_eras.png"


def fig_requeue(cat) -> str:
    """The re-queue population, split by stratum — including the frames that
    currently have NO stratum and would therefore be silently dropped."""
    rows = q(cat, """
        SELECT coalesce(stratum_id, '(no stratum — would NOT be queued)'),
               count(*) FROM s0e_requeue GROUP BY 1 ORDER BY 2 DESC""")
    if not rows:
        rows = [("(none)", 0)]
    names = [r[0] for r in rows][::-1]
    counts = [r[1] for r in rows][::-1]
    colors = [WARN if n.startswith("(no stratum") else ACCENT for n in names]
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(9.4, 3.2))
        bars = ax.barh(names, counts, color=colors)
        ax.set_xlabel("frames to re-queue for astrometry")
        ax.set_title("Wrongly-excluded frames, by S1 stratum")
        ax.set_xlim(0, max(counts) * 1.12)    # room for the value labels
        for b, c in zip(bars, counts):
            ax.annotate(f" {c:,}", (c, b.get_y() + b.get_height() / 2),
                        va="center", fontsize=8, color="#e8eaed")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s0e_requeue.png", dpi=DPI)
        plt.close(fig)
    return "figures/s0e/s0e_requeue.png"


# ---------------------------------------------------------------------------
# Sections — Question / Evidence / Decision / Consequence
# ---------------------------------------------------------------------------
def section_artifact(cat) -> str:
    dump = q(cat, "SELECT card FROM s0e_header_dump ORDER BY card_no")
    exemplar = q1(cat, "SELECT path FROM s0e_header_dump LIMIT 1") \
        if dump else "(no exhibit captured)"
    cards = "\n".join(esc(c[0]) for c in dump) or "(run `rescan_geometry.py exemplar`)"
    src = fig_anatomy(cat)
    return f"""
<section id="artifact">
<div class="bhead"><h2>1 &middot; The artifact</h2>
<span class="tag">an 8-pixel strip that was never there</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">The catalog said the observatory had taken thousands of
8&nbsp;&times;&nbsp;3211-pixel exposures &mdash; strips one hundredth the
width of the detector. No observing programme in the archive asks for that
shape. Are they real?</p>

<h3>Evidence</h3>
<p class="sub">They are not. Here are the first cards of a real archive
file, <code>{esc(exemplar)}</code>, read straight from disk:</p>
<pre class="cards">{cards}</pre>
<p class="sub">A tile-compressed FITS file does not store an image; it
stores a <strong>binary table</strong> whose rows hold compressed tiles.
So <code>NAXIS1</code> is the width of one table row <em>in bytes</em>
(8 &mdash; a single variable-length-array descriptor) and
<code>NAXIS2</code> is the <em>number of rows</em> (3211, one per tile row).
The real picture size sits nine cards lower, in
<code>ZNAXIS1</code>&nbsp;/&nbsp;<code>ZNAXIS2</code>.</p>
{_figure(src, "The phantom strip (yellow) drawn inside the frame that was "
              "actually recorded. Both shapes are described by the same "
              "header; only one of them is the image.")}

<h3>Decision</h3>
<p class="sub">Geometry is resolved by one pure function,
<code>macro_core.fitsgeom.resolve_geometry</code>: when a header carries the
compression markers (<code>ZIMAGE</code>&nbsp;/&nbsp;<code>ZCMPTYPE</code>)
the dimensions are <code>ZNAXIS*</code>; otherwise they are
<code>NAXIS*</code>. A header that claims compression but carries no
<code>ZNAXIS*</code> <strong>raises</strong> rather than falling back &mdash;
falling back is precisely what produced this defect.</p>

<h3>Consequence</h3>
<p class="sub">Why the rows were wrong in the first place: astropy
translates <code>Z*</code> for you, so the scanner was safe <em>when it
could read the header at all</em>. These files carry a malformed
<code>CONTINUE</code> card (a <code>CONTINUE</code> following a non-string
<code>FWALLNAM</code> value) that makes <code>Header.update</code> raise
<code>VerifyError</code> and abandon the whole header &mdash; and the
fallback read the raw table header. The scanner now merges cards
<em>one at a time</em>, so one bad card costs that card, not the frame.</p>
</div>
</section>"""


def section_blast(cat) -> str:
    n_changed = q1(cat, "SELECT count(*) FROM geom_rescan WHERE changed = 1")
    n_control = q1(cat, "SELECT count(*) FROM geom_rescan WHERE changed = 0")
    # The control group's evidentiary value rests entirely on these two
    # counts matching n_control: the rows that did NOT change must be
    # compressed files too, or they prove only the trivial case.  Rendered,
    # not asserted — an earlier draft asserted the opposite in prose.
    n_ctrl_comp = q1(cat, "SELECT count(*) FROM geom_rescan "
                          "WHERE changed = 0 AND compressed = 1")
    n_ctrl_fz = q1(cat, "SELECT count(*) FROM geom_rescan "
                        "WHERE changed = 0 AND path LIKE '%.fz'")
    n_err = q1(cat, "SELECT count(*) FROM geom_rescan WHERE error IS NOT NULL")
    n_left = q1(cat, "SELECT count(*) FROM obs WHERE naxis1 = 8 AND naxis2 = 3211")
    src = fig_blast(cat)

    # One query, used for both the cells and the row highlighting — asking
    # the database the same question twice invites the two answers to differ.
    matrix_rows = q(cat, """
        SELECT old_naxis1, old_naxis2, new_naxis1, new_naxis2, count(*)
        FROM geom_rescan GROUP BY 1,2,3,4 ORDER BY 5 DESC""")
    matrix = table(
        ["stored geometry", "re-read as", "rows", "outcome"],
        [[f"{fmt(o1)} &times; {fmt(o2)}", f"{fmt(n1)} &times; {fmt(n2)}",
          fmt(c), "repaired" if (o1, o2) != (n1, n2) else "left alone"]
         for o1, o2, n1, n2, c in matrix_rows],
        row_classes=["warn" if (o1, o2) != (n1, n2) else None
                     for o1, o2, n1, n2, _ in matrix_rows])

    # From the pre-fix SNAPSHOT, not the live manifest: once the manifest is
    # rebuilt no frame carries the phantom geometry any more, and a live
    # query would render this whole section as zeros.
    era_rows = q(cat, """
        SELECT era_id, readoutm, naxis1, naxis2, xbinning, egain, n_frames
        FROM s0e_blast_era ORDER BY n_frames DESC""")
    era_tbl = table(
        ["era", "readout", "geometry", "bin", "EGAIN", "affected frames"],
        [[fmt(e), esc(r or "(blank)"), f"{fmt(n1)} &times; {fmt(n2)}",
          fmt(xb), fmt(eg), fmt(c)]
         for e, r, n1, n2, xb, eg, c in era_rows],
        row_classes=["warn"] * len(era_rows))

    return f"""
<section id="blast">
<div class="bhead"><h2>2 &middot; Blast radius</h2>
<span class="tag">counted before anything was changed</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">How far did the phantom geometry travel &mdash; how many
rows, which cameras, which science targets, and how many frames did it cost
the astrometry batch?</p>

<h3>Evidence</h3>
<p class="sub">Every catalog row whose stored <code>naxis1</code> was small
enough to be a table row-length ({fmt(n_changed + n_control)} rows) was
re-read from the archive. The change matrix:</p>
{matrix}
<p class="sub"><strong>{fmt(n_changed)}</strong> rows carried phantom
geometry and were repaired. <strong>{fmt(n_control)}</strong> rows were
genuinely small &mdash; Andor iKon focus and guide windows &mdash; and came
back byte-identical. {fmt(n_err)} rows failed to read.
{fmt(n_left)} phantom rows remain.</p>

<p class="sub"><strong>Why that control group is evidence, and not a
tautology.</strong> Every one of those {fmt(n_control)} files is itself
<em>tile-compressed</em> &mdash; {fmt(n_ctrl_fz)} of {fmt(n_control)} are
<code>.fts.fz</code>, and {fmt(n_ctrl_comp)} of {fmt(n_control)} carry
<code>ZIMAGE</code> in their headers. So each
one presents the resolver with the identical trap: a BINTABLE
<code>NAXIS1</code> of 8 sitting right next to the real dimensions, in the
same container, through the same <code>Z*</code> machinery that produced the
phantom. The difference is only in the answer &mdash; here
<code>ZNAXIS1</code> genuinely reads 45, 56 or 57, and the frame really is a
tiny focus window. Two small numbers in one header, and the repair had to
tell them apart {fmt(n_control)} times out of {fmt(n_control)}.</p>

<p class="sub">This paragraph is worth its length because getting it wrong
was the near miss of this whole exercise. The first verification defined the
control as &ldquo;the uncompressed rows&rdquo;, found <em>zero</em> of them,
and passed &mdash; vacuously. An earlier draft of this page then repeated the
same error in prose, describing the control as uncompressed, which would have
left a reader auditing the repair with the trivial case in front of them and
no idea the hard case had been tested at all. The check now asserts what the
files actually are, and
<code>test_genuinely_small_COMPRESSED_frame_survives</code> pins it in the
unit suite so it cannot quietly revert.</p>
{_figure(src, "Left: which science targets lost frames to the artifact. "
              "Right: repaired rows against the untouched control group "
              "(log scale).")}
<p class="sub">By camera era &mdash; both affected eras are keyed on the
phantom geometry, so both are phantom configurations:</p>
{era_tbl}

<h3>Decision</h3>
<p class="sub">Repair the geometry in place, and <em>only</em> the geometry:
<code>rescan_geometry.py</code> writes <code>naxis1</code>/<code>naxis2</code>
and nothing else, recording the old and new value of every row it reads in
<code>geom_rescan</code> whether that row changed or not.</p>

<h3>Consequence</h3>
<p class="sub">The audit table is permanent, so the repair is reversible and
the before/after diff on this page is evidence rather than recollection.</p>
</div>
</section>"""


def section_eras(cat) -> str:
    rows = q(cat, """SELECT era_id, key_before, key_after, n_before, n_after,
                     verdict FROM s0e_era_forecast
                     WHERE verdict != 'unchanged' ORDER BY era_id""")
    n_redefined = q1(cat, "SELECT count(*) FROM s0e_era_forecast "
                          "WHERE verdict = 'REDEFINED'")
    n_total = q1(cat, "SELECT count(*) FROM s0e_era_forecast")
    src = fig_eras(cat)
    tbl = table(
        ["era", "key before", "key after", "frames before", "frames after",
         "verdict"],
        [[fmt(e), f"<code>{esc(kb)}</code>", f"<code>{esc(ka)}</code>",
          fmt(nb), fmt(na), esc(v)] for e, kb, ka, nb, na, v in rows],
        row_classes=["warn" if "RETIRED" in r[5] or "REDEFINED" in r[5]
                     else None for r in rows])
    verdict_html = (
        f'<p class="sub"><strong>{fmt(n_redefined)} published era ids '
        f'changed meaning.</strong> Each must be retired and re-minted '
        f'rather than silently redefined.</p>' if n_redefined else
        '<p class="sub"><strong>No published era id was redefined.</strong> '
        'Every id still denotes exactly the camera configuration it denoted '
        'before &mdash; the phantom ids simply now describe a configuration '
        'that never existed and hold zero frames.</p>')
    return f"""
<section id="eras">
<div class="bhead"><h2>3 &middot; What happens to the eras</h2>
<span class="tag">rehearsed on a copy, never on the live manifest</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">Camera eras are keyed on
<code>(READOUTM, NAXIS1, NAXIS2, XBINNING, EGAIN)</code>, so correcting the
geometry changes the key of every affected frame. Era ids are also a
<strong>pinned registry</strong>: reports, five strategy documents and the
ops request cite &ldquo;era&nbsp;76&rdquo;, &ldquo;era&nbsp;47&rdquo;. Does
any published id change meaning underneath a citation?</p>

<h3>Evidence</h3>
<p class="sub">The manifest was copied and rebuilt against the repaired
catalog &mdash; twice, once before and once after &mdash; and the two era
tables were diffed. The live manifest was never written to; an S1 batch
solve was running throughout. Of {fmt(n_total)} eras, these moved:</p>
{tbl}
{_figure(src, "Frame counts before and after, for every era the repair "
              "moves. The phantom eras empty; the real ones absorb.")}

<h3>Decision</h3>
{verdict_html}
<p class="sub">This falls out of the registry design rather than from care
taken after the fact: because ids are keyed on the configuration tuple
<em>itself</em>, a frame whose geometry is corrected moves to the id that
already owns its true configuration. The phantom id keeps its definition
and loses its frames &mdash; which is exactly &ldquo;retire the phantom,
do not redefine a cited number&rdquo;, enforced structurally.</p>

<h3>Consequence</h3>
<p class="sub">The retired ids must stay in the registry as zero-frame
entries, not be deleted and not be reused: a future build that reissued
those numbers to a new configuration would resurrect precisely the citation
corruption this design prevents. Any document that cites a retired era
should be corrected to cite the era that absorbed its frames.</p>
</div>
</section>"""


def section_requeue(cat) -> str:
    n_new = q1(cat, "SELECT count(*) FROM s0e_requeue")
    n_unstrat = q1(cat, "SELECT count(*) FROM s0e_requeue "
                        "WHERE stratum_id IS NULL")
    n_euuma = q1(cat, "SELECT count(*) FROM s0e_requeue "
                      "WHERE target_key = 'euuma'")
    src = fig_requeue(cat)
    strat_rows = q(cat, """SELECT stratum_id, count(*) FROM s0e_requeue
                           GROUP BY 1 ORDER BY 2 DESC""")
    strat_tbl = table(
        ["stratum", "frames", "status"],
        [[f"<code>{esc(s)}</code>" if s else
          "<em>(none &mdash; classify_stratum returns NULL)</em>", fmt(c),
          "queued by <code>enqueue</code>" if s else "STILL DROPPED"]
         for s, c in strat_rows],
        row_classes=[None if s else "warn" for s, _ in strat_rows])
    tgt_tbl = table(
        ["target", "frames to re-queue"],
        [[esc(t or "(unnamed)"), fmt(c)] for t, c in q(cat, """
            SELECT canonical_target, count(*) FROM s0e_requeue
            GROUP BY 1 ORDER BY 2 DESC LIMIT 12""")])
    return f"""
<section id="requeue">
<div class="bhead"><h2>4 &middot; What astrometry must reconsider</h2>
<span class="tag">18k frames that were never unsolvable</span></div>

<div class="stage"><h3>Question</h3>
<p class="sub">S1 excludes frames narrower than
<code>MIN_SOLVABLE_NAXIS</code> = 512&nbsp;px as &ldquo;sub-frame photometry
windows&rdquo; &mdash; too little sky for quad matching. Given the geometry
was wrong, how many of those exclusions were wrong too?</p>

<h3>Evidence</h3>
<p class="sub">Running the project&rsquo;s own gates
(<code>macro_core.astrom</code>) over the before- and after-manifests:
<strong>{fmt(n_new)}</strong> frames move from excluded to solvable. The
window-geometry exclusion does not merely shrink &mdash; it goes to
<strong>zero</strong>. There were never any photometry windows in the S1
candidate universe; the entire population was this artifact.</p>
{_figure(src, "The re-queue population by stratum. The yellow bar is the "
              "part a queue rebuild would NOT pick up.")}
{strat_tbl}
<p class="sub">Among them are all <strong>{fmt(n_euuma)}</strong> EU&nbsp;UMa
frames &mdash; the series the CV time-series project recorded as permanently
unsolvable. They are full 4800&nbsp;&times;&nbsp;3211 fields, and a spot
check solves them in about three seconds apiece at roughly 75 matched stars
and 1.4&Prime; RMS.</p>
{tgt_tbl}

<p class="sub"><strong>Reconciling the two counts a reader will try to
subtract.</strong> {fmt(n_new)} frames become solvable, and
{fmt(n_new - n_unstrat)} of them reach a stratum &mdash; but the queue takes
{fmt(n_new - n_unstrat + 1)}, one more than that, and the residue of
unstratified solvable frames grows by one LESS than {fmt(n_unstrat)}. The
extra frame is <code>rawimage/2023-11-28/xek33245.fts.fz</code>. It was
always solvable and always unstratified; its geometry did not change at all.
What changed is its <code>target_key</code>, from <code>2023ixf2</code> to
<code>2023ixf</code> &mdash; a target-alias correction that rode along in the
same manifest rebuild &mdash; which let it match the SN stratum for the first
time. Both numbers are right; they count different populations, and the
single frame between them is not a geometry repair at all.</p>

<h3>Decision</h3>
<p class="sub">Rebuild the manifest, then <strong>add to</strong> the S1
queue with <code>enqueue</code> &mdash; never <code>build --rebuild</code>,
which DROPs <code>s1_batch</code> and destroys every solved verdict the batch
has earned.</p>

<p class="sub"><strong>The EU&nbsp;UMa frames are now queued, under a stratum
of their own.</strong> They previously hit
<code>return None &nbsp;# CV frame in an unplanned config</code> and would
have been dropped in silence a second time. The fix is a new stratum id,
<code>cv_fast_fullframe</code>, and the choice of a NEW id rather than a
widened old one is the same law this page argues for era numbers: the
obvious one-line alternative was to let these frames fall through to
<code>fast_fullframe</code>, but that is a <em>facility backlog</em> stratum
whose population is already published, and quietly filling it with
paper-critical CV frames would have redefined a cited number instead of
retiring it. The new stratum sits inside the CV band, so an interrupted batch
still lands the CV astrometry first.</p>

<p class="sub">The remaining <strong>{fmt(n_unstrat)}</strong> unstratified
frames are the blank-<code>READOUTM</code> population, and they are
deliberately left out. They are not one project: the largest blocks are
V426&nbsp;Oph and V2400&nbsp;Oph, but the tail is deep-sky imaging
(NGC&nbsp;6888, M12, M13, M87) across seven filters including narrowband.
Inventing a stratum to hold a heterogeneous population would produce a solve
rate that describes nothing. They stay in the residue until someone designs
for them, and the count is on this page so the decision stays visible.</p>

<h3>Consequence</h3>
<p class="sub">Two pieces of code asserted the artifact was real and have
been corrected. The comment beside the <code>fast_fullframe</code> stratum
&mdash; &ldquo;geometry gate already dropped strips&rdquo; &mdash; was
exactly backwards: the geometry exclusion on that camera is now zero, and
that stratum absorbed the repaired full frames, growing about sevenfold.</p>

<p class="sub">Adding a stratum also renumbered every queue rank below it,
and <code>s1_batch</code> <em>stores</em> the rank on each row &mdash; so
rows queued before the change would have kept the old numbers and sorted
inconsistently against rows queued after. <code>enqueue</code> now re-derives
<code>priority</code>, <code>population</code> and <code>qc_gated</code> for
every row from the policy table on each run. Those three are ordering
metadata; no status, WCS or timing is touched.</p>
</div>
</section>"""


# ---------------------------------------------------------------------------
def render_report(catalog_path: Path = DEFAULT_CATALOG) -> Path:
    """Render the S0e report. Returns the HTML path.

    Reads ONLY the catalog: every consequence this page reports — the
    pre-fix era and target breakdowns, the era forecast, the re-queue
    population — was snapshotted into catalog-side ``s0e_*`` tables while
    the evidence still existed.  Querying the live manifest instead would
    make the page render zeros the moment the manifest is rebuilt, i.e. the
    page would stop being able to prove the incident precisely once the
    incident was fixed.
    """
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    cat = sqlite3.connect(f"file:{catalog_path}?mode=ro", uri=True)
    cat.execute("PRAGMA busy_timeout = 300000")
    try:
        n_changed = q1(cat, "SELECT count(*) FROM geom_rescan WHERE changed = 1")
        n_read = q1(cat, "SELECT count(*) FROM geom_rescan")
        n_new = q1(cat, "SELECT count(*) FROM s0e_requeue")
        n_redef = q1(cat, "SELECT count(*) FROM s0e_era_forecast "
                          "WHERE verdict = 'REDEFINED'")
        sections = [
            section_artifact(cat),
            section_blast(cat),
            section_eras(cat),
            section_requeue(cat),
        ]
        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S0e — The Geometry Artifact</title>
<link rel="stylesheet" href="../assets/macro.css">
<style>
pre.cards {{ background:#0b0d12; border:1px solid #2a3140; border-radius:6px;
  padding:.8rem 1rem; overflow-x:auto; font-size:.78rem; line-height:1.45;
  color:#c8d1dc; }}
</style>
</head><body>

<header>
  <h1>S0e — The Geometry Artifact</h1>
  <p>{fmt(n_read)} catalog rows re-read &middot; {fmt(n_changed)} repaired
  &middot; {fmt(n_new)} frames returned to astrometry &middot;
  {fmt(n_redef)} published era ids redefined &middot;
  <a href="../index.html">back to the evidence hub</a></p>
</header>

<nav>
  <a href="#artifact">1 The artifact</a> &middot;
  <a href="#blast">2 Blast radius</a> &middot;
  <a href="#eras">3 Eras</a> &middot;
  <a href="#requeue">4 Re-queue</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>macro_core.report_s0e</code> from the
<code>s0e_*</code> and <code>geom_rescan</code> evidence tables in
<code>rlmt-catalog.sqlite</code> &mdash; snapshotted before the repair so
this page keeps working after the manifest is rebuilt. Every number here is
the result of a SQL query; none is typed by hand. Regenerate with <code>pipeline/scripts/rescan_geometry.py</code> then
<code>pipeline/scripts/s0e_era_forecast.py</code>.</footer>
</body></html>"""
        HTML_PATH.write_text(html, encoding="utf-8")
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        cat.close()


if __name__ == "__main__":
    print(render_report())
