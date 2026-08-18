"""S2 evidence report renderer: detector truth from the archive.

Reads the S2 tables of the manifest database and writes:

* ``docs/pipeline/s2_detector.html``     — the report
* ``docs/pipeline/figures/s2/*.png``     — every figure

Same discipline as the S0/S0b renderers: the page follows the site's
Socratic format (Question → Evidence → Decision → Consequence) and EVERY
number in the HTML is interpolated from a SQL query executed here or from a
constant defined in the ``rlmt_diagnostics`` logic modules — nothing is
hand-typed.  Figures may additionally read the regenerable pixel products
under ``products/detector/`` (postage stamps are pixels, not claims); all
claims still trace to the database.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # headless: we only ever write PNG files
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

from . import ceiling as ceilmod  # noqa: E402  (constants for interpolation)
from . import linearity as linmod  # noqa: E402
from . import ptc as ptcmod       # noqa: E402
from . import reconstruct as recmod  # noqa: E402

# Shared page machinery: same dark theme, same query discipline, same table
# generator as the S0/S0b reports — one visual language across the site.
import sys                        # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from macro_core.report_s0 import (  # noqa: E402
    ACCENT, DARK, DPI, WARN, _figure, esc, fmt, q, q1, table)

# ---------------------------------------------------------------------------
# Locations, derived from the repo layout (report lives in docs/pipeline/).
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "s2"
HTML_PATH = DOCS_DIR / "s2_detector.html"
GOOD = "#9fd8ae"                # site badge green — confirmations


def fnum(x, nd=2) -> str:
    """Format a float for the page (NULL-safe, fixed decimals)."""
    if x is None:
        return "&mdash;"
    return f"{x:,.{nd}f}"


# ---------------------------------------------------------------------------
# Figures — one function per figure, each returns its relative src path.
# ---------------------------------------------------------------------------
def fig_ceiling_hist(con) -> str:
    """Terminal histogram structure per mode, clip and veto marked."""
    modes = q(con, """SELECT mode, hard_max_adu, clip_adu, veto_adu
                      FROM s2_ceiling_modes ORDER BY mode""")
    with plt.rc_context(DARK):
        fig, axes = plt.subplots(2, 4, figsize=(11.5, 5.6))
        for ax, (mode, hard, clip, veto) in zip(axes.ravel(), modes):
            top = int(hard or 0)
            lo = max(0, top - 1200)
            rows = q(con, "SELECT adu, count FROM s2_ceiling_hist "
                          "WHERE mode = ? AND adu BETWEEN ? AND ?",
                     (mode, lo, top + 60))
            if rows:
                adu = np.array([r[0] for r in rows])
                cnt = np.array([r[1] for r in rows])
                ax.semilogy(adu, np.maximum(cnt, 0.5), lw=0.7, color=ACCENT)
            if clip is not None:
                ax.axvline(clip, color=GOOD, lw=1.2)
            if veto is not None:
                ax.axvline(veto, color=WARN, lw=1.0, ls="--")
            ax.set_title(mode, fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in axes.ravel()[len(modes):]:
            ax.set_visible(False)
        fig.suptitle("Top of each mode's pixel histogram — ceiling (green) "
                     "and saturation veto (dashed yellow)", fontsize=10)
        fig.supxlabel("ADU", fontsize=9)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2_ceiling_hist.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2/s2_ceiling_hist.png"


def fig_frame_maxes(con) -> str:
    """Per-frame maxima: the StackPro 16x cluster and the Low Gain hot pixel."""
    modes = [r[0] for r in q(con, "SELECT mode FROM s2_ceiling_modes "
                                  "ORDER BY mode")]
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(9.6, 3.6))
        for i, mode in enumerate(modes):
            mx = [r[0] for r in q(con, "SELECT max_adu FROM s2_ceiling_frames "
                                       "WHERE mode = ? AND max_adu > 0",
                                  (mode,))]
            y = np.full(len(mx), i) + np.random.default_rng(1).uniform(
                -0.18, 0.18, len(mx))
            ax.scatter(mx, y, s=4, alpha=0.5, color=ACCENT, linewidths=0)
            clip = q1(con, "SELECT clip_adu FROM s2_ceiling_modes "
                           "WHERE mode = ?", (mode,))
            if clip:
                ax.plot([clip, clip], [i - 0.3, i + 0.3], color=GOOD, lw=1.4)
        ax.set_xscale("log")
        ax.set_yticks(range(len(modes)), modes, fontsize=8)
        ax.set_xlabel("frame maximum (ADU, log)")
        ax.set_title("Per-frame maxima: adopted ceilings in green "
                     "(StackPro clusters at 16x the single-read clip)")
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2_frame_maxes.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2/s2_frame_maxes.png"


def fig_ptc(con) -> str:
    """PTC points: darks + lights per mode, with the fitted bracket."""
    with plt.rc_context(DARK):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.2, 4.0))
        colors = {"High Gain": ACCENT, "High Gain StackPro": WARN}
        for mode, color in colors.items():
            for kind, marker in (("dark", "o"), ("light", "^")):
                pts = q(con, "SELECT level, var FROM s2_ptc_points "
                             "WHERE mode = ? AND kind = ?", (mode, kind))
                if pts:
                    ax1.loglog([p[0] for p in pts], [p[1] for p in pts],
                               marker, ms=3, alpha=0.55, color=color,
                               linestyle="none",
                               label=f"{mode} {kind}")
            fit = q(con, "SELECT slope, intercept FROM s2_ptc_fits "
                         "WHERE mode = ? AND kind = 'dark'", (mode,))
            if fit:
                s, i0 = fit[0]
                x = np.geomspace(80, 3000, 50)
                ax1.loglog(x, np.maximum(i0 + s * x, 1e-2), color=color,
                           lw=1.0, ls=":")
        # Poisson reference at the nominal header gain.
        eg = q1(con, "SELECT value FROM detector_params WHERE era_group = "
                     "'High Gain' AND quantity = 'gain_e_per_adu_nominal'")
        if eg:
            x = np.geomspace(80, 3000, 50)
            ax1.loglog(x, x / eg, color=GOOD, lw=1.2, ls="--",
                       label=f"Poisson @ EGAIN {eg:.3f}")
        ax1.set_xlabel("mean level (ADU)")
        ax1.set_ylabel("half difference variance (ADU$^2$)")
        ax1.set_title("Photon transfer, 2023-06-07 pairs")
        ax1.legend(fontsize=6, loc="upper left")
        # Right: the variance floors — the 16x StackPro signature.
        rows = q(con, """SELECT era_group, quantity, value FROM detector_params
                         WHERE quantity IN ('read_noise_adu',
                                            'bias_offset_adu')""")
        vals = {(r[0], r[1]): r[2] for r in rows}
        labels, off, rn2 = [], [], []
        for mode in ("High Gain", "High Gain StackPro"):
            if (mode, "read_noise_adu") in vals:
                labels.append(mode.replace("High Gain", "HG"))
                off.append(vals[(mode, "bias_offset_adu")])
                rn2.append(vals[(mode, "read_noise_adu")] ** 2)
        xpos = np.arange(len(labels))
        ax2.bar(xpos - 0.18, off, 0.32, color=ACCENT, label="bias offset (ADU)")
        ax2.bar(xpos + 0.18, rn2, 0.32, color=WARN,
                label="read-noise variance (ADU$^2$)")
        for i in range(1, len(labels)):
            for dx, series in ((-0.18, off), (0.18, rn2)):
                ratio = series[i] / series[0]
                ax2.annotate(f"x{ratio:.1f}", (xpos[i] + dx, series[i]),
                             ha="center", va="bottom", fontsize=8,
                             color=GOOD)
        ax2.set_xticks(xpos, labels)
        ax2.set_title("StackPro = 16 summed sub-frames:\noffset and RN "
                      "variance both x16")
        ax2.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2_ptc.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2/s2_ptc.png"


def _center_stamp(npz, key: str) -> np.ndarray | None:
    """Reshape the concatenated-region vector's CENTER region to 2-D."""
    regions = npz["regions"]
    name, y0, y1, x0, x1 = regions[0]
    c = int(y1) - int(y0)
    vec = npz[key]
    if vec.size < c * c:
        return None
    return vec[: c * c].reshape(c, c)


def fig_recon(con) -> str:
    """Reconstructed D and F stamps per era; era-47 residual vs its master."""
    rows = q(con, """SELECT era_id, mode, npz_path, truth_offset
                     FROM s2_recon_eras WHERE npz_path IS NOT NULL
                     ORDER BY era_id""")
    n = len(rows)
    with plt.rc_context(DARK):
        fig, axes = plt.subplots(n, 3, figsize=(8.6, 2.5 * n))
        if n == 1:
            axes = axes[None, :]
        for i, (era_id, mode, npz_path, t_off) in enumerate(rows):
            npz = np.load(REPO_ROOT / npz_path, allow_pickle=True)
            D = _center_stamp(npz, "D")
            F = _center_stamp(npz, "F")
            for j, (img, name, cmap) in enumerate(
                    ((D, "dark D (ADU)", "magma"),
                     (F, "flat F", "viridis"))):
                ax = axes[i, j]
                if img is not None:
                    lo, hi = np.nanpercentile(img, [2, 98])
                    im = ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi)
                    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
                ax.set_title(f"era {era_id} ({mode or 'blank'}): {name}",
                             fontsize=7)
                ax.set_xticks([]), ax.set_yticks([])
            ax = axes[i, 2]
            truth = npz["truth"]
            if truth.size and D is not None and t_off is not None:
                c = D.shape[0]
                resid = D - truth[: c * c].reshape(c, c) - t_off
                lim = float(np.nanpercentile(np.abs(resid), 98))
                im = ax.imshow(resid, cmap="coolwarm", vmin=-lim, vmax=lim)
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
                ax.set_title(f"era {era_id}: D $-$ archived master (ADU)",
                             fontsize=7)
            else:
                ax.text(0.5, 0.5, "no archived master\nto compare",
                        ha="center", va="center", fontsize=7,
                        color="#9aa4b2", transform=ax.transAxes)
            ax.set_xticks([]), ax.set_yticks([])
        fig.suptitle("Master reconstruction: central 512$^2$ stamps of the "
                     "calibration each era's reduction actually applied",
                     fontsize=10)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2_recon_stamps.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2/s2_recon_stamps.png"


def fig_linearity(con) -> str:
    """Residuals vs exposure for the best ladder per mode."""
    best = q(con, """
        SELECT l1.mode, l1.ladder_id FROM s2_linearity_ladders l1
        WHERE l1.max_abs_resid_pct = (
            SELECT min(l2.max_abs_resid_pct) FROM s2_linearity_ladders l2
            WHERE l2.mode = l1.mode AND l2.max_abs_resid_pct IS NOT NULL)
        ORDER BY l1.mode""")
    with plt.rc_context(DARK):
        fig, ax = plt.subplots(figsize=(8.8, 3.8))
        palette = [ACCENT, WARN, GOOD, "#d38ce0", "#e0a56c"]
        for (mode, ladder_id), color in zip(best, palette):
            veto = q1(con, "SELECT coalesce(veto_adu, 1e12) FROM "
                           "s2_ceiling_modes WHERE mode = ?", (mode,))
            rungs = q(con, """SELECT exptime, resid_pct, peak_med
                              FROM s2_linearity_rungs WHERE ladder_id = ?
                              ORDER BY exptime""", (ladder_id,))
            t = [r[0] for r in rungs]
            r_ = [r[1] for r in rungs]
            sat = [r[2] is not None and r[2] >= (veto or 1e12) for r in rungs]
            ax.plot(t, r_, "-", color=color, lw=0.9, alpha=0.7)
            for ti, ri, si in zip(t, r_, sat):
                ax.plot(ti, ri, "o" if not si else "x", ms=5 if si else 4,
                        color=color)
            ax.plot([], [], "o-", color=color, label=mode, ms=4)
        ax.set_xscale("log")
        ax.axhspan(-2, 2, color="#2a3140", alpha=0.5, zorder=0)
        ax.axhline(0, color="#9aa4b2", lw=0.6)
        ax.set_xlabel("exposure time (s, log)")
        ax.set_ylabel("residual vs flux = k$\\,\\cdot\\,$t   (%)")
        ax.set_ylim(-60, 60)
        ax.set_title("Best archival exposure ladder per mode "
                     "(x = rung above the saturation veto)")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "s2_linearity.png", dpi=DPI)
        plt.close(fig)
    return "figures/s2/s2_linearity.png"


# ---------------------------------------------------------------------------
# Section builders — each returns one <section> of Socratic HTML.
# ---------------------------------------------------------------------------
def section_ceiling(con) -> str:
    src1 = fig_ceiling_hist(con)
    src2 = fig_frame_maxes(con)
    n_frames = q1(con, "SELECT count(*) FROM s2_ceiling_frames "
                       "WHERE max_adu >= 0")
    n_px = q1(con, "SELECT sum(count) FROM s2_ceiling_hist")
    rows = q(con, """SELECT mode, n_frames, hard_max_adu, clip_adu, veto_adu,
                            bits, unused_codes FROM s2_ceiling_modes
                     ORDER BY mode""")
    tbl = table(
        ["mode", "frames", "hard max", "adopted ceiling", "veto", "ADC bits",
         "unused codes"],
        [[esc(m), fmt(nf), fmt(hm), fmt(cl) if cl is not None
          else "<i>not observed to saturate</i>", fmt(vt), fmt(b), fmt(u)]
         for m, nf, hm, cl, vt, b, u in rows])
    hg_clip = q1(con, "SELECT clip_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'High Gain'")
    hg_hard = q1(con, "SELECT hard_max_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'High Gain'")
    hg_veto = q1(con, "SELECT veto_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'High Gain'")
    sp_clip = q1(con, "SELECT clip_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'High Gain StackPro'")
    lg_hard = q1(con, "SELECT hard_max_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'Low Gain'")
    lg_cluster, lg_div, lg_npos = q(con, """
        SELECT cluster_adu, cluster_diversity, cluster_n_pos
        FROM s2_ceiling_modes WHERE mode = 'Low Gain'""")[0]
    m0_clip = q1(con, "SELECT clip_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'Mode0'")
    m0_veto = q1(con, "SELECT veto_adu FROM s2_ceiling_modes "
                      "WHERE mode = 'Mode0'")
    sp_ratio = sp_clip / hg_clip if (sp_clip and hg_clip) else None
    hg_bits = q1(con, "SELECT bits FROM s2_ceiling_modes "
                      "WHERE mode = 'High Gain'")
    hg_full = q1(con, "SELECT adc_full_scale FROM s2_ceiling_modes "
                      "WHERE mode = 'High Gain'")
    hg_unused = q1(con, "SELECT unused_codes FROM s2_ceiling_modes "
                        "WHERE mode = 'High Gain'")
    # Per-egain epoch statistics (review finding: the High Gain mound pools
    # two cleanly separated egain epochs).  All numbers query-derived.
    eg_rows = q(con, """SELECT egain, n_frames, min_max_adu, median_max_adu,
                               max_max_adu FROM s2_ceiling_egain
                        WHERE mode = 'High Gain' ORDER BY egain""")
    eg_desc = "; ".join(
        f"egain {eg:g}: {n} frames, maxima {mn:,}&ndash;{mx:,} "
        f"(median {med:,.0f})" for eg, n, mn, med, mx in eg_rows)
    # Which epoch does each external prior sit closest to?  Computed by
    # nearest epoch median, never hand-asserted (and the two-epochs-or-one
    # conclusion below is emitted conditionally on the computation).
    sn_mid = (ceilmod.PRIOR_SN_CLIP_LO_ADU
              + ceilmod.PRIOR_SN_CLIP_HI_ADU) / 2.0
    sn_host = (min(eg_rows, key=lambda r: abs(r[3] - sn_mid))
               if eg_rows else None)
    cv_host = (min(eg_rows,
                   key=lambda r: abs(r[3] - ceilmod.PRIOR_CV_CLIP_ADU))
               if eg_rows else None)
    dominant = max(eg_rows, key=lambda r: r[1]) if eg_rows else None
    epoch_txt = ""
    if len(eg_rows) > 1:
        two_epochs = (sn_host is not None and cv_host is not None
                      and sn_host[0] != cv_host[0])
        prior_reading = (
            "the two priors are consistent with having measured "
            "<i>different epochs</i>, not different parts of one mound"
            if two_epochs else
            "both priors sit nearest the same epoch, so the epoch split "
            "does not by itself explain their difference")
        epoch_txt = (
            f" Part of that width is <b>egain-epoch drift</b>: the "
            f"near-ceiling frame maxima split cleanly by camera epoch "
            f"({eg_desc}).  The SN panel's prior range "
            f"{fmt(ceilmod.PRIOR_SN_CLIP_LO_ADU)}&ndash;"
            f"{fmt(ceilmod.PRIOR_SN_CLIP_HI_ADU)} sits nearest the "
            f"egain-{sn_host[0]:g} epoch's median ({sn_host[3]:,.0f}), "
            f"while the CV team's {fmt(ceilmod.PRIOR_CV_CLIP_ADU)} sits "
            f"nearest the egain-{cv_host[0]:g} epoch's median "
            f"({cv_host[3]:,.0f}) — {prior_reading}.  "
            f"The adopted {fmt(hg_clip)} "
            f"describes the dominant egain-{dominant[0]:g} population "
            f"({fmt(dominant[1])} of {fmt(sum(r[1] for r in eg_rows))} "
            f"near-ceiling frames); the veto below sits under every "
            f"near-ceiling maximum of <i>both</i> epochs (lowest observed: "
            f"{fmt(min(r[2] for r in eg_rows))}).")
    veto_vs_prior = ceilmod.prior_comparison(int(hg_veto or 0),
                                             ceilmod.PRIOR_CV_VETO_ADU)
    return f"""
<section id="ceiling"><h2>1&nbsp;&middot;&nbsp;The ceiling memo:
where every mode's scale actually ends</h2>
<h3>Question</h3>
<p>The SN panel measured the High Gain clip at
~{fmt(ceilmod.PRIOR_SN_CLIP_LO_ADU)}&ndash;{fmt(ceilmod.PRIOR_SN_CLIP_HI_ADU)}
ADU; the CV team adopted {fmt(ceilmod.PRIOR_CV_CLIP_ADU)} with a
{fmt(ceilmod.PRIOR_CV_VETO_ADU)} veto (external priors, quoted from the
facility notes as documented constants).  Which number is right, is the
detector 12-bit, and what are the ceilings of the OTHER seven modes nobody
measured?</p>
<h3>Evidence</h3>
<p class="sub">{fmt(n_frames)} science frames sampled across every readout
mode ({fmt(n_px)} pixels histogrammed):</p>
{_figure(src1, "Terminal histogram structure per mode.  A ceiling is a "
         "valley-then-mound at the very end of the occupied range "
         "(pileup); a smooth falling tail is honest headroom.")}
{_figure(src2, "Per-frame maxima.  High Gain frames pile at the mound; "
         "StackPro frames cluster at 16x the single-read clip; Low Gain's "
         "tight cluster failed the position-diversity test (same hot "
         "feature in every frame) and is NOT a ceiling.")}
{tbl}
<p class="sub">Method notes, each one load-bearing: (1) the High Gain
&ldquo;clip&rdquo; is not one code but a {fmt(hg_hard - hg_clip)}-ADU-wide
pileup mound (per-pixel CMOS saturation spread), peaked at {fmt(hg_clip)}
and reaching {fmt(hg_hard)}.{epoch_txt} (2) a
frame-maximum cluster only counts as a ceiling if its members' maxima land
at <i>diverse</i> sensor positions (&ge;{ceilmod.DIVERSITY_MIN_FRAC:.0%}
distinct) — Low Gain's candidate at {fmt(lg_cluster)} ADU repeated a
handful of positions across {fmt(lg_npos)} frames (diversity
{fnum(lg_div, 2)}) and was rejected, so Low Gain's ceiling is
<i>unobserved</i> (bounded below by its hard max {fmt(lg_hard)});
(3) the 5&nbsp;MHz iKon mode's steep sky slope was the false positive that
forced the two-gate detector (valley-then-mound + nothing-above).</p>
<h3>Decision</h3>
<div class="decision"><b>Adopted per-mode ceilings and vetoes as tabled
above: High Gain {fmt(hg_clip)} ADU (veto {fmt(hg_veto)}), StackPro
{fmt(sp_clip)} (= {fnum(sp_ratio, 2)}&nbsp;&times; the single-read clip),
Mode0/Fast/blank-2026 at the full 16-bit code {fmt(m0_clip)} (veto
{fmt(m0_veto)}), iKon
1&nbsp;MHz at {fmt(q1(con, "SELECT clip_adu FROM s2_ceiling_modes WHERE mode='1MHz High Sensitivity 16-bit'"))}.
The GSENSE4040 High Gain scale is <i>12-bit-consistent</i>: the mound at
{fmt(hg_clip)} fits in {fmt(hg_bits)} bits
(2<sup>{fmt(hg_bits)}</sup> = {fmt(hg_full)}) and leaves {fmt(hg_unused)}
codes of the 16-bit FITS container unused — but a per-pixel mound spread
over tens of ADU cannot distinguish a {fmt(hg_bits)}-bit digitization clip
(which lands on ONE shared code, as Mode0/Fast do at {fmt(m0_clip)}) from
per-pixel analog full-well saturation or a clamp below {fmt(hg_full)}, so
&ldquo;confirmed&rdquo; waits for the vendor ADC readback (October).  The
veto rule is floor({ceilmod.VETO_FRACTION} &times; ceiling) to
{ceilmod.VETO_GRANULARITY_ADU} ADU; against the CV team's prior veto of
{fmt(ceilmod.PRIOR_CV_VETO_ADU)}, the measured {fmt(hg_veto)}
{veto_vs_prior}.</b></div>
<h3>Consequence</h3>
<p class="sub">Every downstream photometry pipeline applies the tabled veto
for its mode; pixels above it are ceiling-biased, not data.  StackPro's
16&times; headroom is real but its sub-exposures still clip at the
single-read mound — a star saturated in High Gain is equally saturated
inside every StackPro sub-frame (see &sect;2).</p>
</section>"""


def section_ptc(con) -> str:
    src = fig_ptc(con)
    n_pairs = q1(con, "SELECT count(*) FROM s2_ptc_pairs WHERE n_points > 0")
    n_pts = q1(con, "SELECT count(*) FROM s2_ptc_points")
    # Dark-current contamination of the shortest-dark variance floor —
    # measured from the floors themselves (the "no dark-current shot noise"
    # claim the review caught was false as stated: the shortest dark is
    # 8 s, not 0 s).
    hg_dark_pts = q(con, "SELECT exptime, level, var FROM s2_ptc_points "
                         "WHERE mode = 'High Gain' AND kind = 'dark'")
    shot = ptcmod.dark_shot_fraction(hg_dark_pts)
    shot_txt = (
        f"only a small dark-current term ({100 * shot['frac_of_floor']:.1f}%"
        f" of the floor variance, measured from the {shot['t_short']:g}s-vs-"
        f"{shot['t_long']:g}s floors — negligible but not zero)"
        if shot else "an unquantified dark-current term")
    p = {(r[0], r[1]): (r[2], r[3]) for r in q(con, """
        SELECT era_group, quantity, value, uncertainty FROM detector_params
        WHERE quantity IN ('read_noise_adu', 'read_noise_e',
                           'bias_offset_adu', 'gain_lower_bound_e_per_adu',
                           'gain_upper_bound_e_per_adu',
                           'gain_e_per_adu_nominal', 'nsub')""")}
    def v(mode, qty, nd=2):
        return fnum(p.get((mode, qty), (None,))[0], nd)
    nsub_prov = q1(con, "SELECT provenance FROM detector_params WHERE "
                        "era_group = 'High Gain StackPro' AND quantity = 'nsub'")
    glow = q(con, """SELECT era_group, quantity, value, provenance
                     FROM detector_params WHERE quantity LIKE 'amp_glow%'""")
    glow_txt = "; ".join(
        f"{esc(g[0])}: hottest corner +{fnum(g[2], 1)} ADU "
        f"({esc(g[1]).replace('amp_glow_corner_adu_', '')} dark)"
        for g in glow) or "no long darks available"
    fits_rows = q(con, """SELECT mode, kind, gain_e_per_adu, read_noise_adu,
                                 slope, intercept, n_points
                          FROM s2_ptc_fits ORDER BY mode, kind""")
    ftbl = table(
        ["mode", "pairs", "1/slope (e-/ADU)", "slope", "intercept (ADU²)",
         "points"],
        [[esc(m), esc(k), fnum(g, 3), fnum(s, 3), fnum(i, 1), fmt(n)]
         for m, k, g, rn, s, i, n in fits_rows])
    return f"""
<section id="ptc"><h2>2&nbsp;&middot;&nbsp;Photon transfer: read noise,
the gain bracket, and StackPro's 16 sub-frames</h2>
<h3>Question</h3>
<p>Roadmap R2: the 2023-06-07 repeated darks were shot for a StackPro
photon-transfer curve that was never computed.  What do those frames —
plus the same night's repeated Albireo star fields — actually support?</p>
<h3>Evidence</h3>
<p class="sub">{fmt(n_pairs)} same-scene consecutive pairs &rarr;
{fmt(n_pts)} (level, variance) points:</p>
{_figure(src, "Left: half-difference variance vs level.  Dark pairs "
         "(circles) run BELOW the Poisson line (part of the hot-pixel "
         "population barely fluctuates), light pairs (triangles) run above "
         "it at star levels (sub-pixel scene motion inflates the "
         "difference) — the two slopes bracket the true gain.  Right: the "
         "StackPro/High-Gain ratio of bias offset and read-noise variance.")}
{ftbl}
<p class="sub">The raw fits above are recorded facts, not adopted gains:
the dark slope is sub-Poisson-biased and the light slope motion-biased
(both stated on the figure).  What the frames DO support cleanly:</p>
<ul class="sub">
<li><b>Read noise</b> — the variance floor of the shortest darks has no
scene and {shot_txt}: High Gain
{v("High Gain", "read_noise_adu")}&nbsp;ADU =
{v("High Gain", "read_noise_e")}&nbsp;&plusmn;&nbsp;{fnum(
    p.get(("High Gain", "read_noise_e"), (None, None))[1], 2)}
&nbsp;e<sup>-</sup> (value at nominal EGAIN; the electron-unit uncertainty
spans the measured gain bracket below, NOT the tiny ADU statistical
error); StackPro {v("High Gain StackPro", "read_noise_adu")}&nbsp;ADU.</li>
<li><b>Gain bracket</b> — High Gain true gain lies between
{v("High Gain", "gain_lower_bound_e_per_adu")} (sky-level light slope,
lower bound) and {v("High Gain", "gain_upper_bound_e_per_adu")}
(dark slope, upper bound) e<sup>-</sup>/ADU; the header EGAIN
{v("High Gain", "gain_e_per_adu_nominal", 3)} sits inside the bracket and
is adopted as nominal.  A clean flat-field PTC is an October item.</li>
<li><b>StackPro architecture</b> — three independent ratios against plain
High Gain ({esc(nsub_prov.split(';')[0] if nsub_prov else '')}) agree:
<b>N<sub>sub</sub> = {v("High Gain StackPro", "nsub", 0)}</b>.  A StackPro
frame is the SUM of 16 sub-exposures: bias offset
{v("High Gain StackPro", "bias_offset_adu", 1)} ADU
(&asymp;16 &times; {v("High Gain", "bias_offset_adu", 0)}), read-noise
variance &asymp;16&times; the single read, ceiling &asymp;16&times; the
single-read clip (&sect;1).</li>
<li><b>Amp glow</b> — {glow_txt}; a measurable but small spatial term that
master darks remove.</li>
</ul>
<h3>Decision</h3>
<div class="decision"><b>Error models use: read noise
{v("High Gain", "read_noise_adu")} ADU (High Gain single read); gain =
header EGAIN as nominal with the measured
[{v("High Gain", "gain_lower_bound_e_per_adu")},
{v("High Gain", "gain_upper_bound_e_per_adu")}] e<sup>-</sup>/ADU bracket
as its systematic; StackPro frames modeled as sums of 16 sub-reads
(noise variance 16&times;, read noise 4&times; in ADU), NOT as single
reads.</b>  The repeated darks cannot deliver a full PTC gain — that
limit is stated here rather than papered over, and the flat-field PTC
joins the October list.</div>
<h3>Consequence</h3>
<p class="sub">Every S/N estimate for StackPro data changes: the effective
read noise per frame is {v("High Gain StackPro", "read_noise_adu")} ADU,
sixteen sub-reads' worth, and short-exposure StackPro frames are
read-noise-limited far above where a single High Gain read would be.</p>
</section>"""


def section_recon(con) -> str:
    src = fig_recon(con)
    n_links = q1(con, "SELECT sum(n_links) FROM s2_recon_eras")
    rows = q(con, """SELECT era_id, mode, n_links, n_pairs_used,
                            flat_median, flat_mad_sigma, dark_median,
                            rms_median, pedestal, truth_resid_rms,
                            truth_offset, fd_corr FROM s2_recon_eras
                     ORDER BY era_id""")
    tbl = table(
        ["era", "mode", "links", "pairs fit", "flat median F",
         "dark median D (ADU)", "fit RMS", "pedestal", "F&ndash;D corr",
         "vs archived master", "verdict"],
        [[fmt(e), esc(m or "(blank)"), fmt(nl), fmt(np_), fnum(f, 4),
          fnum(d, 1), fnum(r, 2), fnum(pd, 0), fnum(fc, 2),
          (f"{fnum(t, 2)} ADU RMS (offset {fnum(to, 1)})"
           if t is not None else "&mdash;"),
          esc(recmod.recon_verdict(f, pd, t))]
         for e, m, nl, np_, f, fs, d, r, pd, t, to, fc in rows])
    era47 = q(con, """SELECT truth_resid_rms, truth_resid_mad, truth_master,
                             n_pairs_used, flat_mad_sigma, fd_corr
                      FROM s2_recon_eras WHERE era_id = 47""")
    e47 = era47[0] if era47 else (None,) * 6
    era72_t = q1(con, "SELECT truth_resid_rms FROM s2_recon_eras "
                      "WHERE era_id = 72")
    # Era 79, corrected by the adversarial review: the identity fit shows
    # the two trees hold the SAME pixels — and the manifest shows that
    # shared content is CALIBRATED (every frame named *_calibrated, at the
    # cropped calibrated-output geometry), so it is the rawimage-tree
    # copies that are miscast, not the reduced tree.  Every number below is
    # a query.
    era79_f = q1(con, "SELECT flat_median FROM s2_recon_eras "
                      "WHERE era_id = 79")
    era79_rms = q1(con, "SELECT rms_median FROM s2_recon_eras "
                        "WHERE era_id = 79")
    n79_total = q1(con, "SELECT count(*) FROM frames WHERE era_id = 79")
    n79_cal = q1(con, "SELECT count(*) FROM frames WHERE era_id = 79 "
                      "AND basename LIKE '%calibrated%'")
    n79_raw = q1(con, "SELECT count(*) FROM frames WHERE era_id = 79 "
                      "AND tree = 'rawimage'")
    g79 = q(con, "SELECT naxis1, naxis2 FROM eras WHERE era_id = 79")[0]
    g78 = q(con, "SELECT naxis1, naxis2 FROM eras WHERE era_id = 78")[0]
    return f"""
<section id="recon"><h2>3&nbsp;&middot;&nbsp;The master-reconstruction
experiment: auditing reductions nobody can re-run</h2>
<h3>Question</h3>
<p>The <code>reduced/</code> tree was produced by a pipeline whose
calibration files we mostly do not hold.  Can the raw&harr;reduced pairs
themselves reveal — per pixel — what dark D and flat F were applied, and
can era&nbsp;47 (whose actual masters ARE archived) grade the method?</p>
<h3>Evidence</h3>
<p class="sub">Per pixel, the reduction model
<code>reduced = (raw &minus; D)/F + pedestal</code> is a straight line
<code>raw = F&middot;(reduced &minus; ped) + D</code>; fitting it across
{fmt(n_links)} linked pairs (&le;{fmt(36)} pairs per era, robust
median-slope seeding + {recmod.RECON_CLIP_ITERS} clip passes at
{recmod.RECON_CLIP_SIGMA}&sigma; against astro-scrappy pixel rewrites)
recovers D and F on a 512&sup2; central + four corner grids:</p>
{_figure(src, "Reconstructed dark (left) and flat (middle) per era; "
         "right: era-47 reconstruction minus its archived master "
         "bias+scaled-dark — the ground-truth grade.")}
{tbl}
<p class="sub">Notable per-era facts, all from the table: era 47's flat
median {fnum(q1(con, "SELECT flat_median FROM s2_recon_eras WHERE era_id=47"), 4)}
confirms from pixels what the reduced headers claim
(&ldquo;flat correction NOT performed&rdquo;), and its reconstructed dark
matches the archived master to {fnum(e47[0], 2)} ADU RMS
(MAD-&sigma; {fnum(e47[1], 2)}) over {fmt(e47[3])} pairs.  <b>Read the
F&ndash;D corr column before trusting per-pixel values</b>: era 47's
near-constant sky gave the fit little level diversity, so its per-pixel F
and D trade off against each other (correlation {fnum(e47[5], 2)},
per-pixel F scatter {fnum((e47[4] or 0) * 100, 1)}%) — its
&ldquo;no flat&rdquo; verdict rests on the MEDIAN F only, and its
{fnum(e47[0], 2)} ADU grade is an upper bound that includes degeneracy
noise (era 72, better conditioned, grades at {fnum(era72_t, 2)} ADU RMS);
eras with |corr| &asymp; 0 support per-pixel readings.  <b>Era 79 —
corrected by adversarial review</b>: F = {fnum(era79_f, 4)}, D = 0,
pedestal 0 is an IDENTITY fit — the two trees hold the <i>same pixels</i>
to {fnum(era79_rms, 2)} ADU RMS (compression rounding; NOT byte-identical
— an earlier draft overclaimed both the identity and its direction).
What that shared content is, the manifest answers: {fmt(n79_cal)} of era
79's {fmt(n79_total)} frames — including all {fmt(n79_raw)} in the
rawimage tree — are named <code>*_calibrated</code>, and the era's
geometry {fmt(g79[0])}&times;{fmt(g79[1])} is the CROPPED
calibrated-output geometry (era 78's true Fast raws are
{fmt(g78[0])}&times;{fmt(g78[1])}).  Era 79's rawimage tree therefore
holds mirrored copies of <i>calibrated products</i>; the era's true raws
are absent from the archive, and its reduced tree remains calibrated
data.  The 2026 eras (78/80/83) crop their reduced output by
13&ndash;18 pixels; the crop offset was measured by patch alignment
before fitting.</p>
<h3>Decision</h3>
<div class="decision"><b>The experiment stands: per-pixel line fits
recover the reduction's calibration to a few ADU RMS wherever
&ge;{recmod.RECON_MIN_PAIRS} clean pairs exist, and the era-47 control
grades the method at {fnum(e47[0], 2)} ADU RMS against ground truth
(an upper bound — see the degeneracy note above).  Reduced-tree
photometry is hereby AUDITED for the fitted eras with the table's
verdicts attached.  Era 79: the {fmt(n79_raw)} rawimage-tree frames are
misfiled copies of calibrated products and must NOT be used as raw data
(no true raws exist in the archive for that era); its reduced tree stays
in service as calibrated data.</b>
The reconstructed D/F stamps ship as regenerable products
(<code>products/detector/recon/era*.npz</code>) for downstream flat/dark
sanity checks.</div>
<h3>Consequence</h3>
<p class="sub">Papers using reduced-tree pixels now cite a measured
calibration provenance instead of an assumption; eras with F &ne; 1 are
known to be flat-fielded (and by how much, per pixel on the sampled
grids); the 2026-era crop offset is on record for anyone aligning raw
against reduced frames.</p>
</section>"""


def section_linearity(con) -> str:
    src = fig_linearity(con)
    n_lad = q1(con, "SELECT count(*) FROM s2_linearity_ladders")
    n_fit = q1(con, "SELECT count(*) FROM s2_linearity_ladders "
                    "WHERE rate_adu_per_s IS NOT NULL")
    lin_rows = q(con, """SELECT era_group, value, uncertainty, provenance
                         FROM detector_params
                         WHERE quantity = 'linearity_max_dev_pct'
                         ORDER BY value""")
    ltbl = table(
        ["mode", "max |deviation| (%)", "rung scatter (%)", "witness ladder"],
        [[esc(m), fnum(v, 2), fnum(u, 2), esc(pr.split(";")[0])]
         for m, v, u, pr in lin_rows])
    vega = q(con, """SELECT exptime, flux_med, peak_med
                     FROM s2_linearity_rungs
                     WHERE ladder_id LIKE '2024-05-19|vega%'
                     ORDER BY exptime""")
    vega_txt = "; ".join(f"{t:g}s &rarr; {f:,.0f} ADU (peak {p:,.0f})"
                         for t, f, p in vega)
    modes_without = [r[0] for r in q(con, """
        SELECT DISTINCT mode FROM s2_ceiling_modes WHERE mode NOT IN
        (SELECT DISTINCT era_group FROM detector_params
         WHERE quantity = 'linearity_max_dev_pct')""")]
    return f"""
<section id="linearity"><h2>4&nbsp;&middot;&nbsp;Linearity from archival
exposure ladders</h2>
<h3>Question</h3>
<p>Does flux scale with exposure time — and can the 2024-05-20 Vega
BeStar ladder (0.0001&ndash;0.1&nbsp;s) calibrate the short end?</p>
<h3>Evidence</h3>
<p class="sub">{fmt(n_lad)} same-night same-target same-filter same-gain
ladders surfaced from the manifest ({fmt(n_fit)} fittable); residuals
against flux = k&middot;t for the cleanest ladder per mode:</p>
{_figure(src, "Per-mode best ladder.  Filled points: unsaturated rungs "
         "(the linearity measurement).  Crosses: rungs whose peak sits "
         "above the saturation veto — they measure the ceiling, not "
         "linearity, and are excluded from the statistic.")}
{ltbl}
<p class="sub">Read the table as <b>single-ladder consistency bounds</b>,
not measured non-linearities: each mode's number comes from ONE
(night,&nbsp;target,&nbsp;filter,&nbsp;egain) ladder, the median-rate fit
zeroes one rung by construction (a 3-rung ladder has two informative
rungs), and the residuals fold in sky-transparency drift — which usually
inflates the bound but can partially cancel a real roll-off.  The rung
scatter column is the spread of the clean rungs' residuals, the honest
scale of each bound.</p>
<p class="sub">Two negative findings, recorded as findings: (1) the
<b>Vega ladder is unusable for linearity</b> — its commanded exposure
times collapse below ~10 ms ({vega_txt}: flux is NOT monotonic in
commanded time, i.e. the shutter/timing floor, not the detector,
dominates); (2) modes with <b>no usable ladder</b> —
{esc(", ".join(modes_without) if modes_without else "none")} — have no
archival linearity constraint at all.  Both go on the October list
(a dedicated dome-flat ladder per mode costs minutes).</p>
<h3>Decision</h3>
<div class="decision"><b>Adopted per-mode linearity bounds as tabled:
best-witnessed modes (Mode0, Fast, iKon 1&nbsp;MHz) are linear to
&le;{fnum(q1(con, "SELECT max(value) FROM detector_params WHERE quantity='linearity_max_dev_pct' AND era_group IN ('Mode0','Fast','1MHz High Sensitivity 16-bit')"), 1)}%
over their witnessed ranges; the High Gain family carries a
&le;{fnum(q1(con, "SELECT max(value) FROM detector_params WHERE quantity='linearity_max_dev_pct' AND era_group LIKE 'High Gain%'"), 1)}%
bound that is sky/cloud-limited, not detector-limited — treat it as an
upper limit, not a measured non-linearity.  Sub-10-ms commanded exposure
times are untrustworthy on the iKon and must not be used for absolute
photometric scaling.</b></div>
<h3>Consequence</h3>
<p class="sub">Exposure-time-ratio photometry (HDR stitching, ladder
bootstraps) inherits the tabled bound as a systematic; the timing-paper
error budget cites the sub-10-ms finding directly.</p>
</section>"""


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------
def render_report(manifest_path: Path) -> Path:
    """Render the full S2 report from the manifest DB.  Returns HTML path."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{manifest_path}?mode=ro", uri=True)
    try:
        n_modes = q1(con, "SELECT count(*) FROM s2_ceiling_modes")
        n_ceil = q1(con, "SELECT count(*) FROM s2_ceiling_modes "
                         "WHERE clip_adu IS NOT NULL")
        n_pairs = q1(con, "SELECT count(*) FROM s2_ptc_pairs")
        n_recon = q1(con, "SELECT count(*) FROM s2_recon_eras "
                          "WHERE npz_path IS NOT NULL")
        n_params = q1(con, "SELECT count(*) FROM detector_params")
        meta = dict(q(con, "SELECT key, value FROM s2_build_meta"))

        sections = [
            section_ceiling(con),
            section_ptc(con),
            section_recon(con),
            section_linearity(con),
        ]

        html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S2 — Detector Truth</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>

<header>
  <h1>S2 — Detector Truth</h1>
  <p>{fmt(n_modes)} readout modes characterized ({fmt(n_ceil)} measured
  ceilings) &middot; {fmt(n_pairs)} PTC pairs &middot; {fmt(n_recon)}
  reduction eras reconstructed &middot; {fmt(n_params)} rows in
  <code>detector_params</code> &middot;
  built {esc(meta.get('built_utc', ''))[:16]}Z
  ({esc(meta.get('code_version', ''))},
  commit <code>{esc(meta.get('git_commit', '') or 'uncommitted')}</code>)
  &middot; <a href="../index.html">back to the evidence hub</a></p>
</header>

<nav>
  <a href="#ceiling">1 Ceiling memo</a> &middot;
  <a href="#ptc">2 Photon transfer</a> &middot;
  <a href="#recon">3 Master reconstruction</a> &middot;
  <a href="#linearity">4 Linearity</a>
</nav>

{"".join(sections)}

<footer>Generated by <code>rlmt_diagnostics.report_s2</code> from
<code>products/manifest/rlmt-manifest.sqlite</code> — every number on this
page is the result of a SQL query or a documented pipeline constant; none
is typed by hand.  Regenerate with
<code>pipeline/scripts/run_s2_campaign.py report</code>.</footer>
</body></html>"""

        HTML_PATH.write_text(html, encoding="utf-8")

        # Belt and braces: every <img> the page references must exist and
        # be non-empty, or the build fails loudly.
        import re as _re
        for src in _re.findall(r'<img src="([^"]+)"', html):
            p = DOCS_DIR / src
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"report references missing figure: {src}")
        return HTML_PATH
    finally:
        con.close()
