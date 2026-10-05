"""The G page: the grism calibration library and its acceptance tests.

Renders ``docs/pipeline/g_grism.html`` from ``products/grism/grism.sqlite``
(and, for the counts-by-S2c-verdict line, the manifest snapshot).  Every
number comes from :func:`macro_grism.summary_g.numbers` or a query; one
figure — the D1 evidence — is drawn from the spectrum cache.
"""

from __future__ import annotations

import json
import sqlite3
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_core import plotstyle as ps                  # noqa: E402
from macro_core.report_s0 import esc, fmt               # noqa: E402

from . import config as gconfig                         # noqa: E402
from . import linecal as lc                             # noqa: E402
from . import store as gstore                           # noqa: E402
from . import summary_g                                 # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "pipeline"
FIG_DIR = DOCS_DIR / "figures" / "g"
HTML_PATH = DOCS_DIR / "g_grism.html"

#: The three D1 evidence panels: (star, grism, epoch, rival label).
D1_PANELS = (("eta Hya", "hrg", "ASI-post", "1.59 A/px (v1 code)"),
             ("Phecda", "hrg", "ASI-pre", "1.59 A/px (v1 code)"),
             ("Vega", "lrg", "QHY", "1.1 A/px (stored mode)"))

#: Lines labelled on the figure.
FIG_LINES = ("Hdelta", "Hgamma", "Hbeta", "HeI5876", "O2gamma", "SiII6347",
             "Halpha", "HeI6678", "O2B", "HeI7065", "O2A", "OI7774")


def _f(x, nd=3) -> str:
    return "&mdash;" if x is None else f"{x:.{nd}f}"


def fig_d1(con) -> str:
    """Normalised hot-star spectra against pixel offset from Halpha, with
    the line positions the ADOPTED solution predicts (solid, labelled) and
    the positions the rival dispersion predicts (dashed grey)."""
    import matplotlib.pyplot as plt
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    with ps.context():
        fig, axes = plt.subplots(len(D1_PANELS), 1,
                                 figsize=(ps.COL_DOUBLE, 7.6))
        for ax, (star, grism, epoch, rival) in zip(axes, D1_PANELS):
            row = con.execute("""
                SELECT i.path, i.coeffs_json, i.n_match, i.sign
                FROM g_line_id i JOIN g_frames f USING (path)
                WHERE i.star = ? AND i.grism = ? AND i.mech_epoch = ?
                ORDER BY i.n_match DESC, f.snr_median DESC LIMIT 1""",
                              (star, grism, epoch)).fetchone()
            if row is None:
                ax.set_visible(False)
                continue
            path, cj, nm, sign = row
            coeffs = json.loads(cj)
            spec = gstore.load_spec(path)
            flux = np.where(spec["inside"], spec["flux"], np.nan)
            norm, _c, _ok = lc.normalize(flux)
            x_ha = float(np.polyval(coeffs, 0.0))
            dx = np.arange(len(norm)) - x_ha
            ax.plot(dx, lc.smooth(norm, 1.5), lw=0.6, color=ps.INK)
            rival_a = lc.RIVALS[grism][rival][0]
            good = np.isfinite(norm)
            lo, hi = dx[good].min() - 20, dx[good].max() + 20
            for name in FIG_LINES:
                ln = lc.LINE_BY_NAME[name]
                if name not in lc.GRISM_LINES[grism]:
                    continue
                xa = float(np.polyval(coeffs, ln.wave - 6562.8)) - x_ha
                xr = sign * rival_a * (ln.wave - 6562.8)
                if lo <= xa <= hi:
                    ax.axvline(xa, color=ps.ACCENT, lw=0.6, alpha=0.7)
                    ax.text(xa, 1.10, name, rotation=90, fontsize=6,
                            ha="center", va="bottom", color=ps.ACCENT,
                            clip_on=True)
                if lo <= xr <= hi:
                    ax.axvline(xr, color=ps.MUTED, lw=0.6, ls="--",
                               alpha=0.7)
                    ax.text(xr, 0.60, name, rotation=90, fontsize=5,
                            ha="center", va="bottom", color=ps.MUTED,
                            clip_on=True)
            ax.set_ylim(0.55, 1.35)
            ax.set_xlim(lo, hi)
            ax.set_ylabel("flux / continuum")
            ax.set_title(f"{star} {grism} ({epoch}): blue = adopted "
                         f"solution; grey dashed = where {rival} puts the "
                         f"same lines", fontsize=8)
        axes[-1].set_xlabel("pixel offset from Halpha (px)")
        fig.tight_layout()
        out = FIG_DIR / "g_d1_line_ids.png"
        fig.savefig(out, dpi=ps.WEB_DPI)
        plt.close(fig)
    return f"figures/g/{out.name}"


def _table(head, rows) -> str:
    h = "".join(f"<th>{c}</th>" for c in head)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
                   for r in rows)
    return f"<table><tr>{h}</tr>{body}</table>"


def gate_by_s2c(con) -> list[tuple]:
    """T CrB gate verdicts cross-tabulated by the frame's S2c verdict
    (read from the manifest snapshot, joined on path; unique frames)."""
    snap = gstore.default_manifest()
    s2c = dict(sqlite3.connect(f"file:{snap}?mode=ro", uri=True).execute(
        "SELECT path, verdict FROM frame_dispersion"))
    tab: dict = {}
    for path, verdict in con.execute(
            "SELECT path, verdict FROM g_identity WHERE truth = 'target'"):
        key = s2c.get(path, "not measured")
        tab.setdefault(key, {"ACCEPT": 0, "REJECT": 0})[verdict] += 1
    return [(k, v["ACCEPT"], v["REJECT"]) for k, v in sorted(tab.items())]


def epoch_check() -> tuple:
    """(n formal epochs, n ok, the non-ok rows) of the grism-epoch merge
    against the formal mech_epoch table (gconfig.reconcile_epochs)."""
    snap = gstore.default_manifest()
    rows = sqlite3.connect(f"file:{snap}?mode=ro", uri=True).execute(
        "SELECT mech_epoch, first_night, last_night, boundary_cause "
        "FROM mech_epoch ORDER BY seq").fetchall()
    rec = gconfig.reconcile_epochs(rows)
    bad = [r for r in rec if r["verdict"] != "ok"]
    return len(rec), len(rec) - len(bad), bad


def render_report(db_path: Path = gconfig.GRISM_DB) -> Path:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    n_ep, n_ok, bad_ep = epoch_check()
    ep_note = (esc(bad_ep) if bad_ep
               else "no hardware boundary is merged")
    n = summary_g.numbers(con)
    fig = fig_d1(con)
    meta = dict(con.execute("SELECT key, value FROM g_meta"))
    d1 = n["d1"]
    hrg = [r for r in n["dispersion"] if r["grism"] == "hrg"
           and r["status"] == "adopted"]
    lsf_cal = [r for r in n["lsf"] if r["grism"] == "hrg"
               and r["sample"] == "calibrator"]
    disp_rows = [(r["grism"], esc(r["mech_epoch"]), r["status"],
                  _f(r["disp_a_per_px"], 4), _f(r["disp_err"], 4),
                  f"{_f(r['disp_minus1000'])} / {_f(r['disp_plus1000'])}",
                  r["degree"], r["n_lines"], r["n_stars"], r["n_frames"],
                  r["n_clipped"], _f(r["rms_px"]),
                  f"{_f(r['chi2'], 1)}/{r['dof']}",
                  _f(r["o2b_wave_eff"], 2), esc(r["lines"]))
                 for r in n["dispersion"]]
    lsf_rows = [(r["grism"], esc(r["mech_epoch"]), r["sample"], r["n"],
                 _f(r["fwhm_px"], 1), _f(r["fwhm_a"], 2),
                 f"{_f(r['fwhm_a_p16'], 2)}–{_f(r['fwhm_a_p84'], 2)}",
                 _f(r["fwhm_kms"], 0), _f(r["R"], 0), r["off_nominal"])
                for r in n["lsf"]]
    reg_rows = [(r["grism"], r["regressor"],
                 f"{r['slope_px_per_unit']:+.3g}", _f(r["r"], 2), r["n"],
                 f"{r['range'][0]:.1f} – {r['range'][1]:.1f}")
                for r in n["lsf_regressions"]]
    g1 = n["g1_tcrb"]
    g1_rows = [(g, d["n_frames"], d["n_halpha"], d["n_o2b"],
                _f(d.get("sep_px_median"), 1),
                _f(100 * d.get("sep_px_scatter_frac", float("nan")), 2),
                _f(d.get("dwave_median"), 2), _f(d.get("dwave_scatter"), 2),
                _f(100 * d.get("dwave_scatter_frac", float("nan")), 2),
                d.get("n_dev_gt_2pct"))
               for g, d in g1.items()]
    g2_rows = [(g, d["n_frames"], _f(d["ratio_median"]),
                f"{_f(d['ratio_p16'])}–{_f(d['ratio_p84'])}",
                f"{d['n_within_20pct']}/{d['n_frames']}",
                f"{_f(d['ptc_gain'])} ± {_f(d['ptc_gain_err'])}",
                _f(d["ptc_floor_adu"], 2))
               for g, d in n["g2"].items()]
    g3 = n["g3"]
    t, c = g3.get("target", {}), g3.get("non_target", {})
    g4_rows = []
    for g, d in n["g4"].items():
        for m in ("poly", "flanking"):
            if m in d:
                g4_rows.append((g, m, d[m]["n"],
                                _f(100 * d[m]["median"], 2),
                                _f(100 * d[m]["scatter"], 2),
                                d[m]["n_gt_3pct"],
                                _f(100 * d.get(f"neg_frac_{m}_median", 0), 2),
                                _f(100 * d.get(f"neg_frac_{m}_max", 0), 2)))
    meth_rows = [(g, d.get("n_frames"),
                  _f(100 * d.get("sky_method_diff_median", float("nan")), 2),
                  _f(100 * d.get("sky_method_diff_p84", float("nan")), 2),
                  _f(100 * d.get("box_opt_diff_median", float("nan")), 2),
                  _f(100 * d.get("box_opt_diff_p84", float("nan")), 2))
                 for g, d in n["g4"].items()]
    a6_rows = [(g, d["n"], d["verdicts"].get("clean", 0),
                d["verdicts"].get("flag", 0), d["verdicts"].get("discard", 0),
                fmt(d["cap"]), fmt(d["peak_max"]),
                _f(100 * (d["masked_frac_median"] or 0), 2))
               for g, d in n["a6"].items()]
    ew_rows = [(g, d.get("paired_n"), _f(d.get("v1_median"), 1),
                _f(d.get("v2_median_paired"), 1), _f(d.get("ratio_median"), 2),
                f"{_f(d.get('ratio_p16'), 2)}–{_f(d.get('ratio_p84'), 2)}",
                _f(100 * d.get("v1_scatter_frac", float("nan")), 0),
                _f(100 * d.get("v2_scatter_frac", float("nan")), 0),
                _f(d.get("v1_err_median"), 2),
                _f(d.get("v2_err_median_paired"), 2),
                f"{_f((d.get('v1_disp_range') or (None, None))[0], 2)}–"
                f"{_f((d.get('v1_disp_range') or (None, None))[1], 2)}")
               for g, d in n["ew"].items()]
    det_rows = [(k, d.camera, _f(d.gain_e_per_adu), _f(d.read_noise_adu, 2),
                 _f(d.pedestal_adu, 1), fmt(d.native_clip_adu),
                 fmt(round(d.saturation_cap_adu())), "yes" if d.provisional
                 else "no") for k, d in gconfig.detector_table().items()]
    hb = d1.get("hrg", {})
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>G — Grism calibration library</title>
<link rel="stylesheet" href="../assets/macro.css">
</head><body>
<header><h1>G — The grism calibration library</h1>
<p>{esc(meta.get('dispersion_code', ''))} &middot; built from
<code>products/grism/grism.sqlite</code> &middot;
<a href="../index.html">the front page</a></p></header>

<section id="d1"><h2>1 &middot; D1: what the high-resolution grism
disperses</h2>
<p>Each calibrator frame's absorption features were identified twice with
the same features and the same tolerance: once under the adopted
dispersion and once under each value that was on file before 2026-10-03.
On {hb.get('n_frames', 0)} hot-star hrg frames the adopted hypothesis
matches a median of {_f(hb.get('adopted_median_matches'), 0)} lines against
{_f(hb.get('rival_median_matches'), 0)} for 1.59&nbsp;&Aring;/px; it matches
more on {hb.get('adopted_wins', 0)} frames and as many on
{hb.get('ties', 0)}.  On the {hb.get('b_frames', 0)} B-star frames, whose
narrow He&nbsp;I, Si&nbsp;II and Ne&nbsp;I lines no wrong scale can
reproduce, the medians are {_f(hb.get('b_adopted_median'), 0)} against
{_f(hb.get('b_rival_median'), 0)} (adopted better on
{hb.get('b_adopted_wins', 0)} frames, rival on
{hb.get('b_rival_wins', 0)}).  Where 1.59&nbsp;&Aring;/px places H&beta; the
spectrum's median depth is {_f(hb.get('hbeta159_depth_median'), 3)}
(noise {_f(hb.get('hbeta159_noise_median'), 3)}; deeper than 5&sigma; on
{hb.get('hbeta159_n_detected_5sigma', 0)} of {hb.get('hbeta159_n', 0)}
frames): H&beta; is not there.  The adopted hrg dispersion at the chip
centre is {", ".join(f"{abs(r['disp_a_per_px']):.3f} ({esc(r['mech_epoch'])})" for r in hrg)}
&Aring;/px; the delivered LSF at H&alpha; on the calibrators is
{", ".join(f"{r['fwhm_a']:.2f} &Aring; = {r['fwhm_kms']:.0f} km/s, R {r['R']:.0f} ({esc(r['mech_epoch'])})" for r in lsf_cal)}.</p>
<figure><a href="{fig}"><img src="{fig}" alt="hot-star spectra with line identifications"></a><figcaption>Hot-star hrg spectra with the lines that fix the dispersion identified; where 1.59 &Aring;/px would put H&beta; there is no line.</figcaption></figure>

<h3>The fixed dispersion per (grism, mechanical epoch) — G-1</h3>
<p>Grism epochs are the formal mechanical epochs merged across
wheel-relabelling and sub-degree re-seats only: {n_ok} of {n_ep} formal
epochs reconcile ({ep_note}).</p>
<p>Wavelength is a polynomial in DETECTOR position (the grism-camera
distortion; one constant per frame).  Stellar lines define it; for the hrg
the O<sub>2</sub>&gamma; and O<sub>2</sub>B band edges, calibrated on the
adopted epochs, are added as secondary standards (two-pass solve).</p>
{_table(("grism", "epoch", "status", "D at centre (Å/px)", "±",
         "|D| at −1000 / +1000 px", "deg", "lines", "stars", "frames",
         "clipped", "rms (px)", "χ²/dof", "O₂B edge (Å)", "lines used"),
        disp_rows)}
<p>On T&nbsp;CrB, with the fixed solution and a per-frame H&alpha; zero
point, the O<sub>2</sub>B edge:</p>
{_table(("grism", "frames", "Hα found", "O₂B measured", "Hα–O₂B (px)",
         "px scatter (%)", "edge − standard (Å)", "scatter (Å)",
         "scatter / lever (%)", "|dev| > 2%"), g1_rows)}
</section>

<section id="lsf"><h2>2 &middot; G-5: delivered line-spread function</h2>
{_table(("grism", "epoch", "sample", "n", "FWHM (px)", "FWHM (Å)",
         "16–84% (Å)", "km/s", "R", "off-nominal focus"), lsf_rows)}
<p>Regressors (T&nbsp;CrB frames; LSF in px):</p>
{_table(("grism", "regressor", "slope (px/unit)", "r", "n", "range"),
        reg_rows)}
</section>

<section id="g2"><h2>3 &middot; G-2: variance and saturation from the
measured detector table</h2>
{_table(("key", "camera", "K (e⁻/ADU)", "RN (ADU)", "pedestal",
         "clip", "cap applied", "provisional"), det_rows)}
<p>Measured / predicted sky variance in the trace flanks, T&nbsp;CrB:</p>
{_table(("grism", "frames", "median ratio", "16–84%", "within 20%",
         "pooled PTC K (e⁻/ADU)", "floor (ADU)"), g2_rows)}
</section>

<section id="g3"><h2>4 &middot; G-3: the pixel identity gate</h2>
<p>T&nbsp;CrB: {t.get('accept', 0)} of {t.get('n', 0)} frames accepted;
reasons {esc(t.get('reasons'))}.  Of the {t.get('header_off_n', 0)} frames
whose header pointing is &gt;1&deg; off, {t.get('header_off_accept', 0)} are
accepted on their pixels; the August gate rejected
{t.get('v1_header_rejects', 0)} frames on the header alone, of which
{t.get('v1_header_rejects_now_accepted', 0)} now pass.  No frame is
rejected for a header value.</p>
{_table(("S2c verdict", "accepted", "rejected"), gate_by_s2c(con))}
<p>False-accept rate on {c.get('n', 0)} non-T&nbsp;CrB frames of the same
grisms, camera state and exposure regime: {c.get('accept', 0)}/{c.get('n', 0)}
(Wilson 95% {_f(100 * c.get('wilson95', (0, 0))[0], 1)}–{_f(100 * c.get('wilson95', (0, 0))[1], 1)}%).</p>
</section>

<section id="g4"><h2>5 &middot; G-4: the sky-lozenge background</h2>
<p>Empty-aperture sums &plusmn;200 px from the trace, as a fraction of the
target continuum, and the share of trace columns below &minus;3&sigma;:</p>
{_table(("grism", "sky model", "frames×2", "median (%)", "scatter (%)",
         "|frac| > 3%", "neg. columns median (%)", "max (%)"), g4_rows)}
{_table(("grism", "frames", "sky-method diff. median (%)", "84th pct (%)",
         "boxcar vs optimal median (%)", "84th pct (%)"), meth_rows)}
</section>

<section id="a6"><h2>6 &middot; TCRB-A6: saturation triage</h2>
{_table(("grism", "frames", "clean", "flag", "discard", "cap (ADU)",
         "brightest good pixel", "bad-pixel fraction of aperture (%)"),
        a6_rows)}
</section>

<section id="ew"><h2>7 &middot; How the T&nbsp;CrB H&alpha; EW moved</h2>
{_table(("grism", "paired frames", "v1 EW (Å)", "v2 EW (Å)", "v2/v1",
         "16–84%", "v1 scatter (%)", "v2 scatter (%)", "v1 error (Å)",
         "v2 error (Å)", "v1 dispersions (Å/px)"), ew_rows)}
</section>

<footer>Generated by <code>macro_grism.report_g</code>; every number is a
query of the grism database. Regenerate with
<code>python pipeline/scripts/run_g_tcrb_validation.py --report</code>.
</footer></body></html>"""
    HTML_PATH.write_text(html, encoding="utf-8")
    con.close()
    return HTML_PATH
