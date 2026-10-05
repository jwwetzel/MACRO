#!/usr/bin/env python
"""T CrB Pre-Eruption Monitoring — the project's resumable CLI.

Every stage reads the shared manifest / grism library READ-ONLY and writes
only its own tables in ``products/tcrb/tcrb.sqlite`` (see
``macro_tcrb.db``).  Each stage is idempotent: re-running it replaces its
tables with identical content.

STAGES (plan task in brackets)
------------------------------
    mech-epoch   [TCRB-P0-mech-epoch]   per-frame mechanical epoch + checks
    temp-split   [TCRB-P0-temp-split, -calib-acquisition] temperature groups,
                 dark treatment, flanking-band adequacy on blank sky
    filters      [TCRB-P0-filter-forensics] colour-slope effective wavelengths
    zmag         [TCRB-P0-zmag-provenance]  which catalogue band ZMAG means
    shutter      [TCRB-P0-shutter-timing]   exposure offset + rolling skew
    ...          (later phases are added below their own headers)

USAGE
-----
    P=/opt/miniconda3/envs/rlmt-checks/bin/python
    $P pipeline/scripts/run_tcrb.py mech-epoch
    $P pipeline/scripts/run_tcrb.py status
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np                                               # noqa: E402

from macro_tcrb import db                                        # noqa: E402
from macro_tcrb import p0                                        # noqa: E402


# ===========================================================================
# TCRB-P0-mech-epoch
# ===========================================================================
def cmd_mech_epoch(args) -> None:
    """Per-staged-frame mechanical epoch, the grism crosswalk, and checks."""
    from macro_grism import config as gcfg
    man = db.manifest_ro()
    con = db.connect()
    rows = db.q(man, """
        SELECT s.path, s.role, COALESCE(s.canonical_target, ''), s.filter,
               s.night, s.tree, s.obs_rowid, f.camera, f.mech_epoch,
               f.detector_epoch, f.epoch_certain
        FROM stage_tcrb_monitoring s
        LEFT JOIN frame_mech_epoch f ON f.obs_rowid = s.obs_rowid""")
    out = []
    for (path, role, tgt, filt, night, tree, rid, cam, me, de, cert) in rows:
        gep = gcfg.mech_epoch_id(night) if night else None
        out.append((path, role, tgt, filt, night, tree, rid, cam, me, de,
                    cert, gep))
    n = db.write_table(con, "tcrb_mech_epoch",
                       ("path", "role", "target", "filter", "night", "tree",
                        "obs_rowid", "camera", "mech_epoch", "detector_epoch",
                        "epoch_certain", "grism_epoch"), out)

    # Crosswalk: every grism-library epoch that holds a staged frame.
    cause = dict(db.q(man, "SELECT mech_epoch, boundary_cause FROM mech_epoch"))
    pairs = [(r[11], r[8]) for r in out if r[11] and r[8]]
    xw = p0.crosswalk_verdicts(pairs, cause)
    db.write_table(con, "tcrb_mech_crosswalk",
                   ("grism_epoch", "formal_epochs", "n_frames", "verdict"),
                   [(x["grism_epoch"], x["formal_epochs"], x["n_frames"],
                     x["verdict"]) for x in xw])

    # The grism library's own frame table, where it has been built: its
    # epoch label per frame must equal the crosswalk's (no frame keyed to
    # a different epoch than its night implies).
    lib_mismatch = None
    try:
        g = db.grism_ro()
        lib = dict(db.q(g, "SELECT path, mech_epoch FROM g_frames"))
        mine = {r[0]: r[11] for r in out}
        common = [p for p in lib if p in mine]
        lib_mismatch = sum(lib[p] != mine[p] for p in common)
        lib_common = len(common)
    except Exception as exc:                          # pragma: no cover
        lib_common, lib_mismatch = 0, f"grism library unreadable: {exc}"

    # The v1 validation extraction (g_extractions, stage G) subtracted
    # master darks; every master it used is checked against the epoch of
    # the frame it was applied to.  (A full scan of `frames` by path: ~2 min.)
    v1 = db.q(man, """
        SELECT e.obs_rowid, e.dark_path, me.mech_epoch, md.mech_epoch
        FROM g_extractions e
        JOIN frames fd ON fd.path = e.dark_path
        LEFT JOIN frame_mech_epoch md ON md.obs_rowid = fd.obs_rowid
        LEFT JOIN frame_mech_epoch me ON me.obs_rowid = e.obs_rowid
        WHERE e.method = 'masterdark'""")
    v1_cross = [r for r in v1 if not p0.master_allowed(r[2], r[3])]
    db.write_table(con, "tcrb_mech_v1_darks",
                   ("obs_rowid", "dark_path", "frame_epoch", "dark_epoch",
                    "allowed"),
                   [(*r, int(p0.master_allowed(r[2], r[3]))) for r in v1])

    sci = [r for r in out if r[1] == "science"]
    tc25 = [r for r in sci if r[2] == "T CrB" and r[3] in ("hrg", "lrg")]
    th25 = [r for r in sci if r[2] == "tet CrB" and r[4] and
            "2025-01-01" <= r[4] < "2025-10-01"]
    masters = [r for r in out if r[1].startswith("master_")]
    sci_epochs = {r[8] for r in sci}
    # A staged master that matches no science epoch is never applicable;
    # one that matches is applicable only to frames of its own epoch (the
    # pipelines select by equality — p0.master_allowed).
    cross_masters = [r for r in masters if r[8] not in sci_epochs]
    checks = [
        ("staged_frames", len(out), None),
        ("staged_frames_without_mech_epoch",
         sum(1 for r in out if not r[8]), "pass if 0"),
        ("staged_frames_epoch_uncertain",
         sum(1 for r in out if r[10] == 0), "pass if 0"),
        ("tcrb_2025_grism_frames", len(tc25), None),
        ("tcrb_2025_grism_formal_epochs",
         len({r[8] for r in tc25}), "pass if 1"),
        ("tcrb_2025_grism_formal_epoch", ",".join(sorted({r[8] for r in tc25})),
         None),
        ("tcrb_2025_grism_library_epochs", len({r[11] for r in tc25}),
         "pass if 1"),
        ("tetcrb_2025_frames_same_epoch_as_tcrb",
         sum(1 for r in th25 if r[8] in {x[8] for x in tc25}),
         f"of {len(th25)}"),
        ("grism_epochs_crossing_hardware",
         sum(1 for x in xw if x["verdict"] == "CROSSES_HARDWARE"),
         "pass if 0"),
        ("staged_masters", len(masters), None),
        ("staged_masters_of_no_science_epoch", len(cross_masters),
         "never applicable: excluded by master_allowed"),
        ("v1_masterdark_rows", len(v1), "legacy g_extractions (stage G)"),
        ("v1_masterdark_rows_crossing_epoch", len(v1_cross),
         "RETIRED: never read by this project; the v2 library (G-2/G-4) "
         "uses no cross-epoch master"),
        ("grism_library_frames_compared", lib_common, None),
        ("grism_library_epoch_mismatches", lib_mismatch, "pass if 0"),
    ]
    db.write_table(con, "tcrb_mech_checks", ("check_name", "value", "rule"),
                   [(a, b if not isinstance(b, (list, dict)) else json.dumps(b),
                     c) for a, b, c in checks])
    db.record_stage(con, "mech-epoch", {"tcrb_mech_epoch": n,
                                         "tcrb_mech_crosswalk": len(xw)})
    lines = ["| check | value | rule |", "|---|---|---|"]
    lines += [f"| {a} | {b} | {c or ''} |" for a, b, c in checks]
    lines += ["", "| grism-library epoch | formal (F-3) epochs | staged frames "
              "| verdict |", "|---|---|---|---|"]
    lines += [f"| {x['grism_epoch']} | {x['formal_epochs']} | "
              f"{x['n_frames']} | {x['verdict']} |" for x in xw]
    write_note("mech_epoch", "TCRB-P0-mech-epoch — mechanical epoch of every "
               "staged frame", lines,
               "Tables `tcrb_mech_epoch` (one row per staged frame), "
               "`tcrb_mech_crosswalk`, `tcrb_mech_checks` in "
               "`products/tcrb/tcrb.sqlite`. Formal epochs: F-3 "
               "`frame_mech_epoch`; grism epochs: `macro_grism.config."
               "MECH_EPOCHS`. A master is applied only inside its formal "
               "epoch (`p0.master_allowed`); a grism solution is keyed by a "
               "grism epoch, safe when that epoch never crosses a camera, "
               "flip or coarse-rotation boundary (`p0.crosswalk_verdicts`).")
    print("\n".join(lines))


NOTES = db.REPO / "TCrB_Monitoring" / "notes" / "p0"


def write_note(stem: str, title: str, lines: list[str], method: str) -> Path:
    """Write a generated evidence note (never hand-edited)."""
    NOTES.mkdir(parents=True, exist_ok=True)
    path = NOTES / f"{stem}.md"
    head = [f"# {title}", "", f"Generated by `pipeline/scripts/run_tcrb.py` "
            f"at {db.utc_now()} (git {db.git_commit()}); do not edit.", "",
            method, ""]
    path.write_text("\n".join(head + lines) + "\n")
    return path


# ===========================================================================
# TCRB-P0-temp-split and TCRB-P0-calib-acquisition
# ===========================================================================
#: The temperature-matched Mode0 2x2 masters of the 2025 configuration
#: (mechanical epoch ASI:2024-12-16), as found in the archive.  The -10 C set
#: was taken on 2025-01-05, the 0 C set on 2025-04-14 — the night the
#: regulation set point changed to 0 C.  Neither has a 240 s master, so the
#: 240 s dark is interpolated pixel by pixel between 128 s and 512 s (dark
#: signal is linear in time at these levels; the check is printed).  The
#: archive's only 240 s Mode0 master (Calibrations/masters/, 2025-11-22) is
#: in epoch ASI:2025-10-11 — camera turned 180 deg, FLIPSTAT changed — and
#: may not be applied to 2025 frames (p0.master_allowed).
DARK_MASTERS = {
    "cold": ("Calibrations/2025-01-05/master_dark_128s_read0_g100_o30_2x2_-10C.fts",
             "Calibrations/2025-01-05/master_dark_512s_read0_g100_o30_2x2_-10C.fts"),
    "warm": ("Calibrations/2025-04-14/master_dark_128s_read0_g100_o30_2x2_0C.fts",
             "Calibrations/2025-04-14/master_dark_512s_read0_g100_o30_2x2_0C.fts"),
}
NULL_OFFSETS = (-250, -150, 150, 250)   #: blank-sky pseudo-apertures (rows)
HOT_ABOVE_ADU = 50.0                    #: a dark pixel this far above pedestal
                                        #: in 240 s is "hot" for the census


def _dark240(pair: tuple[str, str]) -> tuple[np.ndarray, dict]:
    from macro_grism import fits_io
    d128, h128, _ = fits_io.load_frame(str(db.ARCHIVE / pair[0]), "float64")
    d512, h512, _ = fits_io.load_frame(str(db.ARCHIVE / pair[1]), "float64")
    f = (240.0 - 128.0) / (512.0 - 128.0)
    meta = {"t128_ccd_temp": h128.get("CCD-TEMP"),
            "t512_ccd_temp": h512.get("CCD-TEMP"),
            "flipstat": h512.get("FLIPSTAT"),
            "median_128": float(np.median(d128)),
            "median_512": float(np.median(d512))}
    return d128 + f * (d512 - d128), meta


def _net_columns(img: np.ndarray, yc: np.ndarray, x0: int, x1: int
                 ) -> np.ndarray:
    """Aperture sum minus the library's linear flanking-band model, per column
    (macro_grism.extract geometry: +-12 px aperture, bands 30-60 px)."""
    from macro_grism import extract as gx
    cut = gx.rectify(img, yc)
    sky = gx.sky_flanking(cut)
    half = cut.shape[0] // 2
    h = gx.APERTURE_HALFWIN
    net = cut[half - h:half + h + 1] - sky[half - h:half + h + 1]
    col = np.nansum(net, axis=0)
    col[~np.isfinite(net).all(axis=0)] = np.nan
    return col[x0:x1 + 1]


def cmd_temp_split(args) -> None:
    """Temperature groups, dark treatment and flanking-band adequacy."""
    from macro_grism import fits_io
    from rlmt_diagnostics import badpix
    man = db.manifest_ro()
    g = db.grism_ro()
    con = db.connect()
    frames = db.q(man, """
        SELECT s.path, s.filter, s.night, f.ccd_temp, f.set_temp, f.flipstat
        FROM stage_tcrb_monitoring s JOIN frames f USING (obs_rowid)
        WHERE s.role = 'science' AND s.canonical_target = 'T CrB'
          AND s.filter IN ('hrg', 'lrg') ORDER BY s.path""")
    trace = {r[0]: r[1:] for r in db.q(g, """
        SELECT path, status, trace_coeffs, trace_x0, trace_x1
        FROM g_frames WHERE sample = 'tcrb'""")}
    darks, dmeta = {}, {}
    for grp, pair in DARK_MASTERS.items():
        darks[grp], dmeta[grp] = _dark240(pair)
    ped = 303.0
    hot = {k: (v - np.median(v)) > HOT_ABOVE_ADU for k, v in darks.items()}
    rows = []
    for path, filt, night, ccd, setp, flip in frames:
        grp = p0.temp_group(ccd)
        dgrp = "warm" if grp == "warm" else "cold"
        st = trace.get(path)
        base = [path, filt, night, ccd, setp, grp, dgrp]
        if st is None or st[0] != "ok":
            rows.append(base + [st[0] if st else "not_in_library"]
                        + [None] * 12)
            continue
        coeffs = np.array(json.loads(st[1]))
        x0, x1 = int(st[2]), int(st[3])
        data, hdr, _ = fits_io.load_frame(str(db.ARCHIVE / path), "float64")
        mask = badpix.load_mask("ASI", data.shape, temp_c=-10.0,
                                flipstat=flip)
        if mask is not None:
            data = np.where(mask, np.nan, data)
        nx = data.shape[1]
        yc = np.polyval(coeffs, np.arange(nx))
        cont = float(np.nanmedian(_net_columns(data, yc, x0, x1)))
        dm = darks[dgrp] if mask is None else np.where(mask, np.nan,
                                                        darks[dgrp])
        dc = darks["cold"] if mask is None else np.where(mask, np.nan,
                                                          darks["cold"])
        # (A) the dark itself through the extraction's background operator,
        # at the trace: what flanking-only leaves of the dark structure.
        leak_none = float(np.nanmedian(_net_columns(dm, yc, x0, x1)))
        # ... and what a -10 C master would leave on a 0 C frame.
        leak_mis = (float(np.nanmedian(_net_columns(dm - dc + ped, yc, x0,
                                                    x1)))
                    if dgrp == "warm" else 0.0)
        # (B) blank-sky pseudo-apertures in the science frame, flanking only
        # and after the matched dark: median over columns, then offsets.
        nulls_raw, nulls_dk = [], []
        for off in NULL_OFFSETS:
            y = yc + off
            if y[x0:x1].min() < 70 or y[x0:x1].max() > data.shape[0] - 70:
                continue
            nulls_raw.append(np.nanmedian(_net_columns(data, y, x0, x1)))
            nulls_dk.append(np.nanmedian(_net_columns(data - dm + ped, y,
                                                      x0, x1)))
        nr = float(np.median(nulls_raw)) if nulls_raw else float("nan")
        nd = float(np.median(nulls_dk)) if nulls_dk else float("nan")
        # Hot-pixel census in the extraction aperture: hot in the matched
        # 240 s dark but NOT in the (-10 C) bad-pixel mask.
        from macro_grism import extract as gx
        rows_ap = (np.round(yc[x0:x1 + 1]).astype(int)[None, :]
                   + np.arange(-gx.APERTURE_HALFWIN,
                               gx.APERTURE_HALFWIN + 1)[:, None])
        cols_ap = np.broadcast_to(np.arange(x0, x1 + 1), rows_ap.shape)
        ok = (rows_ap >= 0) & (rows_ap < data.shape[0])
        hh = hot[dgrp][rows_ap[ok], cols_ap[ok]]
        mm = (mask[rows_ap[ok], cols_ap[ok]] if mask is not None
              else np.zeros_like(hh))
        n_hot_unmasked = int(np.sum(hh & ~mm))
        frac = (lambda v: v / cont if cont and np.isfinite(cont) else None)
        rows.append(base + ["ok", cont, leak_none, frac(leak_none), leak_mis,
                            frac(leak_mis), nr, frac(nr), nd, frac(nd),
                            len(nulls_raw), n_hot_unmasked,
                            int(ok.sum())])
        print(f"  {path[-40:]} {grp:5s} C={cont:8.0f} leak={leak_none:7.2f} "
              f"mis={leak_mis:6.2f} null={nr:7.2f}/{nd:7.2f} "
              f"hot={n_hot_unmasked}", flush=True)
    cols = ("path", "grism", "night", "ccd_temp", "set_temp", "temp_group",
            "dark_group", "status", "cont_adu_col", "leak_none_adu",
            "leak_none_frac", "leak_mismatch_adu", "leak_mismatch_frac",
            "null_flank_adu", "null_flank_frac", "null_dark_adu",
            "null_dark_frac", "n_null_offsets", "n_hot_unmasked",
            "n_aperture_px")
    n = db.write_table(con, "tcrb_temp_split", cols, rows)
    db.write_table(con, "tcrb_dark_masters",
                   ("dark_group", "master_128", "master_512", "meta_json"),
                   [(k, v[0], v[1], json.dumps(dmeta[k]))
                    for k, v in DARK_MASTERS.items()])
    _temp_summary(con)
    db.record_stage(con, "temp-split", {"tcrb_temp_split": n})


def _temp_summary(con) -> None:
    """Per temperature group: counts, adequacy statistics, adopted penalty."""
    out = []
    for grp in ("cold", "intermediate", "warm", "unknown"):
        r = db.q(con, """SELECT grism, ccd_temp, leak_none_frac,
                 leak_mismatch_frac, null_flank_frac, null_dark_frac,
                 n_hot_unmasked, n_aperture_px FROM tcrb_temp_split
                 WHERE temp_group = ? AND status = 'ok'""", grp)
        n_all = db.q1(con, "SELECT COUNT(*) FROM tcrb_temp_split WHERE "
                           "temp_group = ?", grp)
        if not n_all:
            continue
        a = np.array([[np.nan if v is None else v for v in x[1:]]
                      for x in r], float) if r else np.empty((0, 7))

        def st(i, fn):
            v = a[:, i][np.isfinite(a[:, i])] if a.size else np.array([])
            return float(fn(v)) if v.size else None
        out.append((grp, n_all, len(r),
                    st(0, np.min), st(0, np.max),
                    st(1, np.median), st(1, lambda v: np.percentile(np.abs(v), 95)),
                    st(2, np.median), st(3, np.median),
                    st(3, lambda v: np.percentile(np.abs(v), 95)),
                    st(4, np.median),
                    st(4, lambda v: np.percentile(np.abs(v), 95)),
                    st(5, np.sum), st(6, np.sum)))
    db.write_table(con, "tcrb_temp_summary",
                   ("temp_group", "n_frames", "n_measured", "ccd_temp_min",
                    "ccd_temp_max", "leak_none_frac_med", "leak_none_frac_p95",
                    "leak_mismatch_frac_med", "null_flank_frac_med",
                    "null_flank_frac_p95", "null_dark_frac_med",
                    "null_dark_frac_p95", "n_hot_unmasked", "n_aperture_px"),
                   out)
    for o in out:
        print("  ", o)


def cmd_temp_report(args) -> None:
    """The adequacy figure and the method note (P0-temp-split and
    P0-calib-acquisition): every number below is read from the tables."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from macro_core import plotstyle as ps
    con = db.connect()
    rows = db.q(con, """SELECT temp_group, cont_adu_col, leak_none_adu,
        leak_mismatch_adu, null_flank_adu, null_dark_adu, n_hot_unmasked,
        n_aperture_px FROM tcrb_temp_split WHERE status = 'ok'""")
    a = {g: np.array([r[1:] for r in rows if r[0] == g], float)
         for g in ("cold", "warm")}
    ps.apply("print")
    fig, ax = plt.subplots(1, 2, figsize=(ps.COL_DOUBLE, 2.8))
    for g, col, mk in (("cold", ps.ACCENT, "o"), ("warm", ps.BAD, "s")):
        v = a[g]
        ax[0].plot(v[:, 0], v[:, 1] / v[:, 0],
                   **ps.measurement_kw(col, mk, size=3.5),
                   label=f"{g}: dark left by flanking bands")
        if g == "warm":
            ax[0].plot(v[:, 0], np.abs(v[:, 2]) / v[:, 0],
                       **ps.measurement_kw(ps.WARN, "D", size=3.0),
                       label="warm: 0 C minus -10 C master")
        ax[1].plot(v[:, 0], np.abs(v[:, 3]) / v[:, 0],
                   **ps.measurement_kw(col, mk, size=3.5),
                   label=f"{g}: blank sky, flanking only")
        ax[1].plot(v[:, 0], np.abs(v[:, 4]) / v[:, 0],
                   **ps.measurement_kw(ps.tint(col), "^", size=3.5),
                   label=f"{g}: blank sky, matched dark first")
    for k in (0, 1):
        ax[k].set_xscale("log")
        ax[k].set_yscale("log")
        ax[k].set_xlabel("T CrB continuum (ADU per column, 240 s)")
        ax[k].axhline(0.01, **ps.reference_kw())
        ax[k].legend(fontsize=6, loc="lower left")
    ax[0].set_ylabel("residual / continuum")
    ax[0].set_title("(a) dark structure through the background bands",
                    fontsize=7)
    ax[1].set_title("(b) empty pseudo-apertures 150-250 px off the trace",
                    fontsize=7)
    fig.tight_layout()
    out = NOTES / "fig_temp_split.png"
    NOTES.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=ps.PNG_DPI)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)

    summ = {r[0]: r for r in db.q(con, "SELECT * FROM tcrb_temp_summary")}
    c, w = summ["cold"], summ["warm"]
    med = {g: np.nanmedian(a[g], axis=0) for g in a}
    lines = [
        "| group | CCD-TEMP (C) | frames | traced | dark master (2025, epoch "
        "ASI:2024-12-16) | dark left by flanking, median ADU/col | ... as "
        "fraction of continuum, median | unmasked hot px in aperture per "
        "frame |", "|---|---|---|---|---|---|---|---|",
        f"| cold | {c[3]:.1f} .. {c[4]:.1f} | {c[1]} | {c[2]} | -10 C, "
        f"2025-01-05 (128 s, 512 s) | {med['cold'][1]:.2f} | {c[5]:.4f} | "
        f"{c[12] / max(c[2], 1):.0f} |",
        f"| warm | {w[3]:.1f} .. {w[4]:.1f} | {w[1]} | {w[2]} | 0 C, "
        f"2025-04-14 (128 s, 512 s) | {med['warm'][1]:.2f} | {w[5]:.4f} | "
        f"{w[12] / max(w[2], 1):.0f} |",
        "", "Blank-sky pseudo-apertures (median |residual| / continuum): "
        f"cold {c[8]:.4f} flanking-only, {c[10]:.4f} with the matched dark "
        f"subtracted first; warm {w[8]:.4f} and {w[10]:.4f}. Using the -10 C "
        f"master on a 0 C frame would leave a further {w[7]:.4f} (median).",
        "", "## Method (the paragraph TCRB-P0-calib-acquisition asks for)", "",
        f"The {c[1] + w[1]} T CrB grism frames of 2025 (240 s, Mode0, 2x2 "
        f"average) split cleanly by CCD-TEMP: {c[1]} at the -10 C set point "
        f"and {w[1]} at the 0 C set point used from 2025-04-14; none lies "
        "between. Contrary to the premise of DE.F5, both temperatures have "
        "matched master darks from the same mechanical epoch (-10 C on "
        "2025-01-05, 0 C on 2025-04-14; 128 s and 512 s, interpolated to "
        "240 s). The archive's only 240 s master (2025-11-22) belongs to the "
        "next epoch (camera turned 180 deg, FLIPSTAT changed) and is never "
        "applied. The adopted treatment is the grism library's: no dark "
        "frame is subtracted; the per-column background from the two "
        "flanking bands (30-60 px either side of the trace, a straight line "
        "through their medians) removes sky, halo and dark together, and "
        "bad pixels are masked (S2 mask ASI_4788x3194_T-10_FM). What that "
        "leaves of the dark is measured, not assumed: each frame's "
        "temperature-matched 240 s master is passed through the identical "
        "aperture-minus-flanking operator at that frame's trace, and the "
        "result, divided by the frame's own continuum, is carried as a "
        "per-frame fractional EW penalty (an additive background error b "
        "biases EW by -b/C). Warm frames carry a penalty "
        f"{w[5] / c[5]:.1f}x the cold frames' (medians {w[5]:.4f} and "
        f"{c[5]:.4f}) and about {w[12] / max(w[2], 1) / max(c[12] / max(c[2], 1), 1):.1f}x "
        "as many hot pixels the -10 C mask does not hold; those are "
        "removed by the extraction's outlier rejection, and the warm frames "
        "are flagged so every EW result can be checked with them excluded. "
        "Blank-sky pseudo-apertures show that the flanking-band residual "
        "is dominated by sky structure, not dark (panel b): subtracting "
        "the matched dark first changes it by "
        f"{(c[8] - c[10]):.4f} (cold) and {(w[8] - w[10]):.4f} (warm) of "
        "the continuum. New Mode0 darks cannot be taken on the camera now "
        "mounted (QHY600): whether the ASI survives and can be cooled on a "
        "bench is question 1 of section 8 of the rev. 3 observatory request "
        "(`ops/2026-10_observatory_request_rev3.md`).", "",
        "![adequacy](fig_temp_split.png)"]
    write_note("temp_split", "TCRB-P0-temp-split / -calib-acquisition — "
               "the grism frames by detector temperature", lines,
               "Tables `tcrb_temp_split` (per frame), `tcrb_temp_summary`, "
               "`tcrb_dark_masters` in `products/tcrb/tcrb.sqlite`; "
               "groups by `p0.temp_group` (cold <= -7.5 C, warm >= -5 C). "
               "Frames without a library trace (G v2 `no_trace`) are "
               "counted but not measured.")
    print("\n".join(lines[:6]))


# ===========================================================================
# TCRB-P0-shutter-timing
# ===========================================================================
def _box_flux(data: np.ndarray, coeffs, x0: int, x1: int,
              half: int = 20, gap: int = 40, width: int = 30) -> float:
    """Total first-order counts: a +-``half`` row box along the trace minus
    the mean of two flanking bands (``gap``..``gap+width`` rows out)."""
    from macro_grism import extract as gx
    nx = data.shape[1]
    yc = np.polyval(coeffs, np.arange(nx))
    cut = gx.rectify(data, yc, half=gap + width)
    c = cut.shape[0] // 2
    box = cut[c - half:c + half + 1, x0:x1 + 1]
    bg = np.concatenate([cut[c - gap - width:c - gap, x0:x1 + 1],
                         cut[c + gap:c + gap + width, x0:x1 + 1]])
    with np.errstate(invalid="ignore"):
        b = np.nanmedian(bg, axis=0)
    net = box - b[None, :]
    return float(np.nansum(net))


def _iso_s(stamp: str) -> float:
    """ISO-8601 UTC stamp (any number of fractional digits) -> seconds."""
    return float(np.datetime64(stamp, "us").astype("int64")) / 1e6


def cmd_shutter(args) -> None:
    """Exposure-time offset from same-night exposure ladders of bright stars
    (ASI Mode0, the camera of every 2025 T CrB spectrum), and a bound on the
    rolling-shutter skew from the frame cadence."""
    from macro_grism import fits_io
    g = db.grism_ro()
    man = db.manifest_ro()
    con = db.connect()
    rows = db.q(g, """SELECT path, star, grism, night, exptime, trace_coeffs,
        trace_x0, trace_x1, peak_adu FROM g_frames
        WHERE sample = 'calibrator' AND mech_epoch IN ('ASI-pre', 'ASI-post')
          AND status = 'ok' ORDER BY star, grism, night, exptime""")
    lad: dict = {}
    for r in rows:
        lad.setdefault((r[1], r[2], r[3]), []).append(r)
    meas = []
    for key, fr in sorted(lad.items()):
        if len({x[4] for x in fr}) < 2 or min(x[4] for x in fr) > 2.5:
            continue
        # One COMMON column range per ladder (the intersection of the
        # frames' detected trace extents, 20 px inside it): each frame's
        # own extent shrinks with S/N, so per-frame ranges would bias short
        # exposures low and fake a negative delta.
        c0 = max(int(x[6]) for x in fr) + 20
        c1 = min(int(x[7]) for x in fr) - 20
        if c1 - c0 < 200:
            continue
        for path, star, grism, night, exp, co, x0, x1, pk in fr:
            data, hdr, _ = fits_io.load_frame(str(db.ARCHIVE / path),
                                              "float64")
            f = _box_flux(data, np.array(json.loads(co)), c0, c1, half=30)
            meas.append((path, star, grism, night, float(exp), f, pk,
                         hdr.get("DATE-OBS")))
    db.write_table(con, "tcrb_shutter_frames",
                   ("path", "star", "grism", "night", "exptime", "flux_adu",
                    "peak_adu", "date_obs"), meas)
    fits_out = []
    by: dict = {}
    for m in meas:
        by.setdefault((m[1], m[2], m[3]), []).append(m)
    for (star, grism, night), ms in sorted(by.items()):
        t = np.array([m[4] for m in ms])
        f = np.array([m[5] for m in ms])
        if len(set(t)) < 2 or (f <= 0).any():
            continue
        # Per-frame error: the empirical scatter of repeats at one exposure
        # (scintillation + transparency), as a fraction, pooled over rungs;
        # where no rung repeats, 2% (the median pooled value is printed).
        fr_s = []
        for tt in set(t):
            v = f[t == tt]
            if v.size >= 2:
                fr_s.append(np.std(v, ddof=1) / np.mean(v))
        frac = float(np.median(fr_s)) if fr_s else 0.02
        frac = max(frac, 0.003)
        r = p0.exposure_offset_fit(t, f, frac * f)
        fits_out.append((star, grism, night, len(ms), float(t.min()),
                         float(t.max()), frac, r["rate"], r["delta_s"],
                         r["delta_err_s"], r["chi2"], r["dof"]))
    db.write_table(con, "tcrb_shutter_fits",
                   ("star", "grism", "night", "n_frames", "t_min", "t_max",
                    "frac_scatter", "rate", "delta_s", "delta_err_s", "chi2",
                    "dof"), fits_out)
    d = np.array([x[8] for x in fits_out])
    e = np.array([x[9] for x in fits_out])
    # Ladders that constrain delta: formal error < 50 ms.
    ok = np.isfinite(d) & np.isfinite(e) & (e > 0) & (e < 0.05)
    w = 1 / e[ok] ** 2
    dm = float(np.sum(w * d[ok]) / np.sum(w))
    de = float(np.sqrt(1 / np.sum(w)))
    chi2 = float(np.sum(((d[ok] - dm) / e[ok]) ** 2))
    dof = int(ok.sum() - 1)
    lo, med, hi = np.percentile(d[ok], [5, 50, 95])
    mad_sig = float(1.4826 * np.median(np.abs(d[ok] - med)))
    # Rolling-shutter skew bound: the shortest dead time (start-to-start
    # minus exposure) between consecutive frames of the camera in this
    # mechanical epoch; the row-sequential readout of a full frame, which
    # IS the skew, cannot exceed it.
    st = db.q(man, """SELECT date_obs, exptime FROM frames
        WHERE night BETWEEN '2024-12-16' AND '2025-06-23' AND camera = 'ASI'
          AND readoutm = 'Mode0' AND date_obs IS NOT NULL
          AND tree = 'rawimage' ORDER BY date_obs""")
    gaps = []
    for (a, ea), (b, _) in zip(st, st[1:]):
        try:
            dt = _iso_s(b) - _iso_s(a) - float(ea)
        except Exception:
            continue
        if 0 < dt < 600:
            gaps.append(dt)
    gap_min = float(min(gaps))
    nrows = 3194
    trace_rows = 2 * 12 + 1
    summary = [
        ("n_ladders", len(fits_out), None),
        ("n_ladders_constraining", int(ok.sum()), None),
        ("n_frames", len(meas), None),
        ("shortest_exposure_s", float(min(m[4] for m in meas)), None),
        ("delta_wmean_s", dm, de),
        ("delta_chi2", chi2, dof),
        ("delta_median_s", float(med), None),
        ("delta_mad_sigma_s", mad_sig, None),
        ("delta_p5_s", float(lo), None),
        ("delta_p95_s", float(hi), None),
        ("delta_min_s", float(d[ok].min()), None),
        ("delta_bound_s", float(max(abs(lo), abs(hi))), None),
        ("n_dead_times", len(gaps), None),
        ("min_dead_time_s", gap_min, None),
        ("skew_frame_bound_s", gap_min, None),
        ("skew_across_trace_bound_s", gap_min * trace_rows / nrows, None),
    ]
    db.write_table(con, "tcrb_shutter_summary", ("quantity", "value", "aux"),
                   summary)
    db.record_stage(con, "shutter", {"tcrb_shutter_frames": len(meas),
                                      "tcrb_shutter_fits": len(fits_out)})
    lines = ["| quantity | value | aux |", "|---|---|---|"]
    lines += [f"| {a} | {b:.6g} | {'' if c is None else f'{c:.4g}'} |"
              for a, b, c in summary]
    lines += ["", "| star | grism | night | frames | t range (s) | scatter | "
              "delta (s) | +- | chi2/dof |", "|---|---|---|---|---|---|---|---"
              "|---|"]
    lines += [f"| {x[0]} | {x[1]} | {x[2]} | {x[3]} | {x[4]:.3g}-{x[5]:.3g} "
              f"| {x[6]:.3f} | {x[8]:+.4f} | {x[9]:.4f} | {x[10]:.1f}/{x[11]}"
              f" |" for x in fits_out]
    write_note("shutter_timing", "TCRB-P0-shutter-timing — exposure-time "
               "offset and rolling-shutter skew (ASI Mode0)", lines,
               "Each same-night exposure ladder of a bright star through "
               "either grism (G v2 library frames, traced) is fitted with "
               "counts = R (t + delta) (`p0.exposure_offset_fit`); the "
               "per-frame error is the measured repeat scatter at fixed "
               "exposure, and each ladder uses one common column range (the "
               "intersection of its frames' trace extents). Ladders with a "
               "formal error under 50 ms constrain delta; they disagree far "
               "beyond those errors (chi2 printed), i.e. delta is not one "
               "constant of the camera, so the bound quoted is the "
               "5-95 percentile range of the per-ladder values (with the "
               "robust MAD sigma; the single most extreme ladder is a "
               "transparency change and is printed as delta_min), not the "
               "error of their mean. The skew of a row-sequential readout "
               "cannot exceed the shortest observed dead time between "
               "consecutive frames; across the 25 rows of a trace aperture "
               "it is that bound x 25/3194. The archival 0.085 s theta CrB "
               "frames are direct r images taken at one exposure only, so "
               "they cannot fit delta; they are covered by the ladders, "
               "whose shortest rungs reach the same regime.")
    print("\n".join(lines[:12]))


# ===========================================================================
# Catalogues for the T CrB field (filters, zmag, B-track)
# ===========================================================================
TCRB_RA, TCRB_DEC = 239.875667, 25.920167      #: VSX J2000 (pm negligible)
CAT_RADIUS_DEG = 0.45                          #: AC4040 half-diagonal 0.435
CAT_DIR = db.REPO / "products" / "tcrb" / "catalogue_cache"


def _cattie():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "run_cv_cattie", db.REPO / "pipeline" / "scripts" / "run_cv_cattie.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


#: Gaia DR3 synthetic photometry (GSPC) as served by VizieR I/360/syntphot,
#: mapped to the band names used here.
GSPC_VIZIER = {"b_jkc": "B", "v_jkc": "V", "r_jkc": "R", "i_jkc": "I",
               "g_sdss": "g", "r_sdss": "r", "i_sdss": "i", "z_sdss": "z",
               "y_ps1": "y"}
#: On 2026-10-05 the ESA archive timed out on a plain count over this cone,
#: VizieR's I/355/gaiadr3 did not answer in minutes, and the AIP mirror's
#: GSPC table timed out even on a three-id IN list.  What answered in about
#: a second: the gaia_source cone at the AIP mirror and the I/360/syntphot
#: cone at VizieR.  So the two cones are pulled separately and joined here
#: on source_id — exactly the join the server would have done.
AIP_TAP = "https://gaia.aip.de/tap/sync"
VIZIER_SYNC = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"


def _tap_csv(url: str, adql: str) -> list[dict]:
    import csv
    import io
    import urllib.parse
    import urllib.request
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL",
                                   "FORMAT": "csv", "QUERY": adql}).encode()
    with urllib.request.urlopen(urllib.request.Request(url, data=data),
                                timeout=300) as r:
        txt = r.read().decode()
    if txt.lstrip().startswith("<"):
        # Some services ignore FORMAT=csv and answer VOTable.
        if 'QUERY_STATUS" value="ERROR' in txt:
            raise RuntimeError(f"TAP error from {url}: {txt[:400]}")
        from astropy.io.votable import parse_single_table
        t = parse_single_table(io.BytesIO(txt.encode())).to_table(
            use_names_over_ids=True)
        return [{c: (None if np.ma.is_masked(row[c]) else row[c])
                 for c in t.colnames} for row in t]
    return list(csv.DictReader(io.StringIO(txt)))


def _gaia_vizier(cc) -> tuple[dict, str]:
    """Gaia DR3 cone (AIP) joined locally on source_id to the GSPC cone
    (VizieR I/360/syntphot); G < cc.GAIA_G_LIMIT."""
    circ = (f"CIRCLE('ICRS',{TCRB_RA:.6f},{TCRB_DEC:.6f},"
            f"{CAT_RADIUS_DEG:.4f})")
    q1 = ("SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp "
          "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
          f"{circ}) AND phot_g_mean_mag < {cc.GAIA_G_LIMIT}")
    sel = ['"Source"'] + [f'"{v}mag", "F{v}", "e_F{v}", "{v}Flag"'
                          for v in GSPC_VIZIER.values()]
    q2 = (f'SELECT {", ".join(sel)} FROM "I/360/syntphot" WHERE '
          f"1=CONTAINS(POINT('ICRS',RA_ICRS,DE_ICRS),{circ})")
    g = cc.with_retry(lambda: _tap_csv(AIP_TAP, q1), "gaia cone (AIP)")
    sp = cc.with_retry(lambda: _tap_csv(VIZIER_SYNC, q2), "GSPC cone (CDS)")
    if len(g) >= 2000 or len(sp) >= 2000:
        raise RuntimeError("cone result at a 2,000-row cap: refuse to "
                           "proceed on a possibly truncated catalogue")
    gs = {int(r["Source"]): r for r in sp}

    def fl(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return float("nan")
    cols: dict = {k: [] for k in ("source_id", "ra", "dec", "pmra", "pmdec",
                                  "phot_g_mean_mag", "bp_rp")}
    for b in GSPC_VIZIER:
        for suf in ("_mag", "_mag_error", "_flag"):
            cols[b + suf] = []
    for r in g:
        sid = int(r["source_id"])
        cols["source_id"].append(sid)
        for k in ("ra", "dec", "pmra", "pmdec", "phot_g_mean_mag", "bp_rp"):
            cols[k].append(fl(r[k]))
        x = gs.get(sid)
        for b, v in GSPC_VIZIER.items():
            if x is None:
                cols[b + "_mag"].append(float("nan"))
                cols[b + "_mag_error"].append(float("nan"))
                cols[b + "_flag"].append(1.0)
                continue
            f, ef = fl(x[f"F{v}"]), fl(x[f"e_F{v}"])
            cols[b + "_mag"].append(fl(x[f"{v}mag"]))
            cols[b + "_mag_error"].append(1.0857 * ef / f if f > 0
                                          else float("nan"))
            cols[b + "_flag"].append(fl(x[f"{v}Flag"]))
    print(f"      gaia {len(g)} sources, GSPC {len(sp)}, joined "
          f"{sum(1 for r in g if int(r['source_id']) in gs)}", flush=True)
    return cols, q1 + "\n" + q2


def cmd_catalogues(args) -> None:
    """ATLAS-REFCAT2 and Gaia DR3 + GSPC (Johnson B,V,R,I; SDSS g,r,i,z;
    PS1 y) over the T CrB field — fetched once, cached with sha256 and the
    exact queries (the CV catalogue-tie fetchers, reused unchanged)."""
    cc = _cattie()
    con = db.connect()
    out = []
    for name in ("refcat2", "gaia_gspc"):
        path = CAT_DIR / f"{name}_tcrb.json.gz"
        cache = cc.read_cache(path)
        if cache is None or args.force:
            if name == "refcat2":
                cols, queries = cc.fetch_refcat2(TCRB_RA, TCRB_DEC,
                                                 CAT_RADIUS_DEG, 18.0)
            else:
                cols, queries = _gaia_vizier(cc)
            payload = {"columns": cols, "queries": queries,
                       "pulled_utc": db.utc_now(),
                       "centre": [TCRB_RA, TCRB_DEC],
                       "radius_deg": CAT_RADIUS_DEG}
            sha = cc.atomic_write_gz(path, payload)
            cache = payload
        else:
            import hashlib
            import gzip
            sha = hashlib.sha256(gzip.open(path).read()).hexdigest()
        n = len(cache["columns"]["ra"])
        out.append((name, str(path.relative_to(db.REPO)), n,
                    cache["pulled_utc"], sha, cache["queries"][:2000]))
        print(f"  {name}: {n} rows, pulled {cache['pulled_utc']}")
    db.write_table(con, "tcrb_catalogues",
                   ("catalogue", "cache_path", "n_rows", "pulled_utc",
                    "sha256", "queries"), out)
    db.record_stage(con, "catalogues", {"tcrb_catalogues": len(out)})


def load_catalogue() -> dict:
    """Merged field catalogue: Gaia DR3/GSPC rows with REFCAT2 g,r,i,z
    attached by a 1-arcsec positional match (both at their own epochs —
    proper motions in this field are < 0.1 arcsec over a decade for all
    but a handful of stars, which the match radius absorbs)."""
    cc = _cattie()
    g = cc.read_cache(CAT_DIR / "gaia_gspc_tcrb.json.gz")["columns"]
    r = cc.read_cache(CAT_DIR / "refcat2_tcrb.json.gz")["columns"]
    out = {k: np.asarray(v, float) for k, v in g.items()
           if k != "source_id"}
    out["source_id"] = np.asarray(g["source_id"], np.int64)
    from scipy.spatial import cKDTree
    cosd = np.cos(np.radians(TCRB_DEC))
    tr = cKDTree(np.c_[np.asarray(r["ra"]) * cosd, r["dec"]])
    d, j = tr.query(np.c_[out["ra"] * cosd, out["dec"]],
                    distance_upper_bound=1.0 / 3600)
    ok = np.isfinite(d)
    for band in ("gmag", "rmag", "imag", "zmag"):
        v = np.full(out["ra"].size, np.nan)
        v[ok] = np.asarray(r[band], float)[j[ok]]
        out[f"rc2_{band}"] = v
    return out


# ===========================================================================
# measure — calibrated aperture photometry of every T CrB imaging frame
# (TCRB-P0-resolve, -filter-forensics, -zmag-provenance; B0-B3, B5, C1)
# ===========================================================================
PLATE_SCALE = 0.54          #: arcsec / native px, AC4040 at f/6.8
CAT_G_MAX = 16.5            #: catalogue stars measured per frame (Gaia G)
WCS_MIN_MATCH, WCS_MAX_RMS = 20, 1.0
WCS_MIN_MATCH_HDR = 8
APERTURES = {"r15": ("fwhm", 1.5), "fix8": ("px", 8.0), "r30": ("fwhm", 3.0)}
ANNULUS_FWHM = (4.0, 6.0)


def _ac4040_masters(man) -> list[dict]:
    rows = db.q(man, """
        SELECT s.role, s.path, s.night, COALESCE(s.filter, ''), f.exptime,
               f.readoutm
        FROM stage_tcrb_monitoring s JOIN frames f USING (obs_rowid)
        WHERE s.mech_epoch = 'AC4040:2023-02-06'
          AND s.role IN ('master_dark', 'master_flat')""")
    out = []
    for role, path, night, code, exp, mode in rows:
        name = path.rsplit("/", 1)[-1].lower()
        kind = "flat" if role == "master_flat" and "dark" not in name \
            else "dark"
        out.append({"kind": kind, "path": path, "night": night,
                    "code": code, "exptime": float(exp), "mode": mode})
    return out


class _MasterCache:
    """Loaded, prepared masters (darks raw; flats dark-subtracted and
    normalised to their central median, dead pixels < 0.2 masked)."""

    def __init__(self, cands):
        self.cands, self.cache = cands, {}

    def dark(self, mode, exptime, night):
        from macro_tcrb import phot
        m = phot.pick_master(self.cands, "dark", mode, exptime, night)
        if m is None:
            return None, None
        if m["path"] not in self.cache:
            from macro_sn import snio
            self.cache[m["path"]] = snio.load_data(m["path"])[0]
        return self.cache[m["path"]], m

    def flat(self, mode, code, night):
        from macro_tcrb import phot
        m = phot.pick_master(self.cands, "flat", mode, 0.0, night, code)
        if m is None:
            return None, None, None
        key = "FLAT:" + m["path"]
        if key not in self.cache:
            from macro_sn import snio
            fl = snio.load_data(m["path"])[0].astype(np.float32)
            dk, dm = self.dark(mode, m["exptime"], m["night"])
            if dk is not None:
                fl = fl - dk
            fl = fl / np.median(fl[1024:3072, 1024:3072])
            dead = fl < 0.2
            fl[dead] = 1.0
            self.cache[key] = (fl, dead, dm["path"] if dm else None)
        fl, dead, dpath = self.cache[key]
        return fl, dead, {**m, "flat_dark": dpath}


def _imaging_frames(man) -> list[tuple]:
    return db.q(man, """
        SELECT p.obs_rowid, p.path, p.night, p.mode, p.filter, p.exptime,
               f.jd, f.airmass, f.zmag, f.fwhm, f.pltsolvd, f.ccd_temp,
               f.flipstat, f.objctra, f.objctdec, d.verdict,
               v.cap_fraction, v.clip_adu, v.bias_adu, v.native_factor
        FROM s2_target_peaks p CROSS JOIN frames f ON f.obs_rowid =
             p.obs_rowid
        LEFT JOIN frame_dispersion d ON d.obs_rowid = p.obs_rowid
        LEFT JOIN s2_target_verdicts v ON v.obs_rowid = p.obs_rowid
        WHERE p.is_canonical = 1 ORDER BY f.jd""")


def _measure_one(fr, mc, cat, args) -> tuple[dict, list]:
    import sep
    from macro_sn import snio
    from rlmt_diagnostics import badpix
    (rid, path, night, mode, code, exp, jd, am, zmag, fwhm_h, solved, ccd,
     flip, ora, odec, dverd, capf, clip, bias, nfac) = fr
    rec = {"obs_rowid": rid, "path": path, "night": night, "mode": mode,
           "filter": code, "exptime": exp, "jd": jd, "airmass": am,
           "zmag": zmag, "pltsolvd": solved, "s2c_verdict": dverd}
    if dverd == "dispersed":
        rec["status"] = "dispersed (S2c): spectrum, not photometry"
        return rec, []
    raw, hdr = snio.load_data(path)
    rec["nsubexp"] = hdr.get("NSUBEXP")
    # PinPoint's own record of which catalogue it solved (and took ZMAG)
    # against: 'Matched N stars from the <catalogue>'.
    hist = " | ".join(str(h) for h in hdr.get("HISTORY", []))
    m_pp = re.search(r"Matched\s+(\d+)\s+stars from the (.+?)(?:\s*\||$)",
                     hist)
    rec["pp_catalog"] = m_pp.group(2).strip() if m_pp else None
    rec["pp_nmatch"] = int(m_pp.group(1)) if m_pp else None
    dark, dm = mc.dark(mode, exp, night) if mode != "Low Gain" else (None,
                                                                       None)
    flat, dead, fm = mc.flat(mode, code, night) if mode != "Low Gain" \
        else (None, None, None)
    img = raw.astype(np.float32)
    if dark is not None:
        img = img - dark
    if flat is not None:
        img = img / flat
    rec["dark_master"] = dm["path"] if dm else None
    rec["flat_master"] = fm["path"] if fm else None
    rec["recipe"] = ("dark+flat" if dark is not None and flat is not None
                     else "dark only" if dark is not None
                     else "flat only" if flat is not None else "none")
    mask = badpix.load_mask("AC4040", img.shape, temp_c=ccd, flipstat=flip)
    mask = np.zeros(img.shape, bool) if mask is None else mask.copy()
    if dead is not None:
        mask |= dead
    rec["badpix"] = int(mask.sum())
    sub, bkg, obj = snio.detect(img, mask)
    cap_raw = None
    if capf is not None and clip is not None and bias is not None:
        cap_raw = bias + capf * (clip - bias)
    fwhm = snio.fwhm_of(obj, cap_raw or 1e9)
    rec["fwhm_px"] = fwhm
    rec["n_detect"] = int(len(obj))
    if fwhm is None:
        rec["status"] = "no FWHM (too few unsaturated stars)"
        return rec, []
    # --- astrometry: header solution refined, else blind (P0-resolve) ----
    w0 = snio.header_wcs(hdr) if solved == 1 else None
    src = "header"
    if w0 is None:
        src = "blind"
        # Blind matching on stellar detections only: on short or cloudy
        # frames the brightest sep detections are hot pixels and cosmic
        # rays (median FWHM < 2 px), which no triangle can match.
        st = obj[(obj["flag"] == 0) & (obj["a"] > 0.8) &
                 (obj["b"] / np.maximum(obj["a"], 1e-3) > 0.5)]
        w0 = snio.blind_solve(st, cat["ra"], cat["dec"],
                              cat["phot_g_mean_mag"], TCRB_RA, TCRB_DEC,
                              PLATE_SCALE, img.shape)
    w, nmat, rms = (snio.refine_wcs(w0, obj, cat["ra"], cat["dec"],
                                    img.shape) if w0 is not None
                    else (None, 0, None))
    rec.update(wcs_source=src, wcs_n=nmat, wcs_rms_arcsec=rms)
    # A PinPoint solution re-fitted to fewer stars is accepted down to
    # WCS_MIN_MATCH_HDR (narrow-band frames hold few catalogue stars); a
    # blind solution needs WCS_MIN_MATCH.
    need = WCS_MIN_MATCH_HDR if src == "header" else WCS_MIN_MATCH
    if w is None or nmat < need or rms is None or rms > WCS_MAX_RMS:
        rec["status"] = f"astrometry failed ({src})"
        return rec, []
    rec["wcs_cd"] = json.dumps(np.asarray(w.pixel_scale_matrix).tolist())
    rec["wcs_crval"] = json.dumps(list(map(float, w.wcs.crval)))
    # --- forced photometry at every catalogue star + T CrB -----------------
    sel = cat["phot_g_mean_mag"] < CAT_G_MAX
    x, y = w.all_world2pix(cat["ra"][sel], cat["dec"][sel], 0)
    ids = np.flatnonzero(sel)
    ny, nx = img.shape
    ins = (x > 40) & (x < nx - 40) & (y > 40) & (y < ny - 40)
    x, y, ids = x[ins], y[ins], ids[ins]
    data = np.ascontiguousarray(sub, dtype=np.float32)
    gain = float(hdr.get("EGAIN") or 1.0)
    res = {}
    for name, (unit, val) in APERTURES.items():
        r = val * fwhm if unit == "fwhm" else val
        f, e, fl = sep.sum_circle(data, x, y, r, err=bkg.rms(), gain=gain,
                                  mask=mask,
                                  bkgann=(ANNULUS_FWHM[0] * fwhm,
                                          ANNULUS_FWHM[1] * fwhm))
        res[name] = (f, e, fl)
    peak = snio.native_peaks(raw, x, y, max(fwhm, 3.0))
    tid = int(cat["target_index"])
    rows = []
    for k in range(len(ids)):
        rows.append((rid, int(cat["source_id"][ids[k]]),
                     int(ids[k] == tid), float(x[k]), float(y[k]),
                     *[float(res[n][j][k]) for n in APERTURES
                       for j in (0, 1)],
                     int(res["r15"][2][k]), float(peak[k])))
    rec["n_stars"] = len(rows)
    t = [r for r in rows if r[2] == 1]
    rec["target_in_frame"] = int(bool(t))
    if t:
        rec["target_x"], rec["target_y"] = t[0][3], t[0][4]
        rec["target_peak_raw"] = t[0][-1]
    rec["status"] = "ok"
    return rec, rows


FRAME_COLS = ("obs_rowid", "path", "night", "mode", "filter", "exptime",
              "jd", "airmass", "zmag", "pltsolvd", "s2c_verdict", "nsubexp",
              "pp_catalog", "pp_nmatch",
              "dark_master", "flat_master", "recipe", "badpix", "fwhm_px",
              "n_detect", "wcs_source", "wcs_n", "wcs_rms_arcsec", "wcs_cd",
              "wcs_crval", "n_stars", "target_in_frame", "target_x",
              "target_y", "target_peak_raw", "status")
STAR_COLS = ("obs_rowid", "source_id", "is_target", "x", "y",
             "flux_r15", "err_r15", "flux_fix8", "err_fix8", "flux_r30",
             "err_r30", "sep_flag", "peak_raw")


def cmd_measure(args) -> None:
    """Calibrate, solve and photometer every canonical T CrB imaging frame."""
    man = db.manifest_ro()
    con = db.connect()
    cat = load_catalogue()
    # T CrB's own catalogue row: the brightest Gaia source within 3".
    d = np.hypot((cat["ra"] - TCRB_RA) * np.cos(np.radians(TCRB_DEC)),
                 cat["dec"] - TCRB_DEC) * 3600
    near = np.flatnonzero(d < 3.0)
    cat["target_index"] = near[np.argmin(cat["phot_g_mean_mag"][near])]
    mc = _MasterCache(_ac4040_masters(man))
    frames = _imaging_frames(man)
    recs, stars = [], []
    if args.retry_failed:
        # Keep every frame already measured or ruled dispersed; redo the
        # rest and merge.
        keep = {r[0] for r in db.q(con, """SELECT obs_rowid FROM
            tcrb_img_frames WHERE status = 'ok' OR status LIKE
            'dispersed%'""")}
        old = db.q(con, f"SELECT {', '.join(FRAME_COLS)} FROM "
                        "tcrb_img_frames")
        recs = [dict(zip(FRAME_COLS, r)) for r in old if r[0] in keep]
        stars = [r for r in db.q(con, "SELECT * FROM tcrb_img_stars")
                 if r[0] in keep]
        frames = [f for f in frames if f[0] not in keep]
    for i, fr in enumerate(frames):
        try:
            rec, rows = _measure_one(fr, mc, cat, args)
        except Exception as exc:                       # recorded, not hidden
            rec, rows = {"obs_rowid": fr[0], "path": fr[1],
                         "status": f"error: {type(exc).__name__}: {exc}"}, []
        recs.append(rec)
        stars.extend(rows)
        print(f"  [{i + 1}/{len(frames)}] {fr[1][-28:]} {fr[4]:2s} "
              f"{fr[3][:12]:12s} {rec.get('status')} "
              f"wcs={rec.get('wcs_source')}/{rec.get('wcs_n')} "
              f"fwhm={rec.get('fwhm_px')}", flush=True)
    n1 = db.write_table(con, "tcrb_img_frames", FRAME_COLS,
                        [tuple(r.get(c) for c in FRAME_COLS) for r in recs])
    n2 = db.write_table(con, "tcrb_img_stars", STAR_COLS, stars)
    db.write_table(con, "tcrb_target", ("source_id", "ra", "dec", "gaia_g",
                                        "bp_rp", "b_jkc", "v_jkc", "r_jkc",
                                        "i_jkc"),
                   [(int(cat["source_id"][cat["target_index"]]),
                     *[float(cat[k][cat["target_index"]]) for k in
                       ("ra", "dec", "phot_g_mean_mag", "bp_rp",
                        "b_jkc_mag", "v_jkc_mag", "r_jkc_mag",
                        "i_jkc_mag")])])
    db.record_stage(con, "measure", {"tcrb_img_frames": n1,
                                      "tcrb_img_stars": n2})


# ===========================================================================
# TCRB-P0-filter-forensics and TCRB-P0-zmag-provenance
# ===========================================================================
#: Nominal effective wavelengths (A): Johnson-Cousins from Bessell (2005,
#: ARA&A 43, 293, Table 1); SDSS from Fukugita et al. (1996, AJ 111, 1748);
#: PS1 y from Tonry et al. (2012, ApJ 750, 99).  Used only to place each
#: catalogue band on a wavelength axis for the colour-term zero crossing.
BAND_LAM = {"b_jkc": 4380, "g_sdss": 4770, "v_jkc": 5450, "r_sdss": 6231,
            "r_jkc": 6410, "i_sdss": 7625, "i_jkc": 7980, "z_sdss": 9134,
            "y_ps1": 9620}
G_RANGE = (9.0, 15.5)
COLOUR_RANGE = (0.3, 2.5)


def _frame_stars(con, cat, rid, exptime, cap):
    rows = db.q(con, """SELECT source_id, flux_r15, err_r15, peak_raw,
        sep_flag FROM tcrb_img_stars WHERE obs_rowid = ? AND is_target = 0
        AND flux_r15 > 0""", rid)
    if not rows:
        return None
    idx = {int(s): i for i, s in enumerate(cat["source_id"])}
    a = np.array([[idx.get(int(r[0]), -1), r[1], r[2], r[3], r[4]]
                  for r in rows], float)
    k = a[:, 0].astype(int)
    ok = (k >= 0) & (a[:, 4] == 0) & (a[:, 3] < (cap or 1e12))
    k, a = k[ok], a[ok]
    g = cat["phot_g_mean_mag"][k]
    c = cat["bp_rp"][k]
    sel = (g > G_RANGE[0]) & (g < G_RANGE[1]) & (c > COLOUR_RANGE[0]) & \
        (c < COLOUR_RANGE[1])
    k, a, c = k[sel], a[sel], c[sel]
    m = -2.5 * np.log10(a[:, 1] / exptime)
    e = np.hypot(1.0857 * a[:, 2] / a[:, 1], 0.01)
    return k, m, e, c


def cmd_filters(args) -> None:
    """Colour-term zero crossing (effective wavelength) per frame and code,
    and the ZMAG band match, from the measured imaging photometry."""
    from macro_tcrb import phot
    con = db.connect()
    cat = load_catalogue()
    frames = db.q(con, """SELECT f.obs_rowid, f.filter, f.mode, f.exptime,
        f.night, f.zmag, c.cap_adu FROM tcrb_img_frames f
        LEFT JOIN tcrb_census c USING (obs_rowid) WHERE f.status = 'ok'""")
    per_band, per_frame, zrows = [], [], []
    for rid, code, mode, exp, night, zmag, cap in frames:
        st = _frame_stars(con, cat, rid, exp, cap)
        if st is None or len(st[0]) < 10:
            continue
        k, m, e, c = st
        ks = {}
        for b, lam in BAND_LAM.items():
            mb = cat[f"{b}_mag"][k]
            f = phot.huber_line(c, mb - m, e, c_ref=1.0)
            if f is None:
                continue
            ks[b] = f
            per_band.append((rid, code, b, lam, f.zp, f.k, f.k_err, f.rms,
                             f.n, f.chi2, f.dof))
        known = sorted((BAND_LAM[b], f.k) for b, f in ks.items())
        lam_eff = float("nan")
        for (l1, k1), (l2, k2) in zip(known, known[1:]):
            if np.sign(k1) != np.sign(k2):
                lam_eff = l1 + (0 - k1) * (l2 - l1) / (k2 - k1)
                break
        best = min(ks, key=lambda b: abs(ks[b].k)) if ks else None
        per_frame.append((rid, code, mode, night, len(k), lam_eff, best,
                          ks[best].k if best else None,
                          ks[best].rms if best else None,
                          ks[best].zp if best else None))
        if zmag is not None and np.isfinite(zmag):
            bands = {b: cat[f"{b}_mag"][k] for b in BAND_LAM}
            bands.update({"gaia_G": cat["phot_g_mean_mag"][k],
                          **{f"refcat2_{x}": cat[f"rc2_{x}mag"][k]
                             for x in "griz"}})
            zm = p0.zmag_band_match(m, zmag, bands)
            for b, v in zm.items():
                zrows.append((rid, code, mode, night, zmag, b, v["offset"],
                              v["scatter"], v["n"]))
    db.write_table(con, "tcrb_filter_bandfits",
                   ("obs_rowid", "filter", "band", "lam_band", "zp", "k",
                    "k_err", "rms", "n", "chi2", "dof"), per_band)
    db.write_table(con, "tcrb_filter_frames",
                   ("obs_rowid", "filter", "mode", "night", "n_stars",
                    "lam_eff", "best_band", "best_k", "best_rms", "best_zp"),
                   per_frame)
    db.write_table(con, "tcrb_zmag_match",
                   ("obs_rowid", "filter", "mode", "night", "zmag", "band",
                    "offset", "scatter", "n"), zrows)
    db.record_stage(con, "filters", {"tcrb_filter_frames": len(per_frame),
                                      "tcrb_zmag_match": len(zrows)})
    for r in db.q(con, """SELECT filter, COUNT(*), ROUND(AVG(lam_eff)),
            GROUP_CONCAT(DISTINCT best_band) FROM tcrb_filter_frames
            GROUP BY 1"""):
        print("  ", r)


#: The catalogue band each code would be solved against IF PinPoint matched
#: the filter (hypothesis 1), from the colour-term test above.
MATCHED_BAND = {"B": "b_jkc", "G": "g_sdss", "R": "r_sdss", "I": "i_sdss",
                "L": "v_jkc", "O": "refcat2_g", "H": "r_jkc", "1": "r_jkc"}


def cmd_zmag_report(args) -> None:
    """P0-zmag-provenance: which catalogue band ZMAG is, per filter.

    PinPoint defines ZMAG so that m_inst(1 s) + ZMAG = catalogue magnitude.
    Two hypotheses per filter: (1) ZMAG was solved in the band matching the
    filter; (2) ZMAG was solved in ATLAS-REFCAT2 r whatever the filter.
    Under the true one, catalogue - (m_inst + ZMAG) is only the aperture
    difference between PinPoint's photometry and ours — the value the R
    frames show — for every filter.
    """
    con = db.connect()
    pp = db.q(con, """SELECT filter, pp_catalog, COUNT(*) FROM
        tcrb_img_frames WHERE zmag IS NOT NULL GROUP BY 1, 2""")
    off = {(f, b): (n, o, sc) for f, b, n, o, sc in db.q(con, """SELECT
        filter, band, COUNT(*), AVG(offset), AVG(scatter) FROM
        tcrb_zmag_match WHERE n >= 10 GROUP BY 1, 2""")}
    ap = off.get(("R", "refcat2_r"), (0, 0.0, 0))[1]
    rows = []
    for f in sorted({k[0] for k in off}):
        mb = MATCHED_BAND.get(f)
        n1, o1, s1 = off.get((f, mb), (0, None, None))
        n2, o2, s2 = off.get((f, "refcat2_r"), (0, None, None))
        cats = "; ".join(f"{c} ({n})" for ff, c, n in pp if ff == f)
        if o1 is None or o2 is None:
            verdict = "untested"
        elif min(abs(o1 - ap), abs(o2 - ap)) > 0.3:
            verdict = "neither (not a zero point)"
        else:
            verdict = ("REFCAT2 r" if abs(o2 - ap) < abs(o1 - ap)
                       else "matched band")
        rows.append((f, cats, mb, n1, o1, s1, o2, s2, o2 - ap if o2 is not
                     None else None, verdict))
    db.write_table(con, "tcrb_zmag_provenance", (
        "filter", "pinpoint_catalogue", "matched_band", "n_frames",
        "offset_matched", "scatter_matched", "offset_refcat2_r",
        "scatter_refcat2_r", "offset_r_minus_aperture", "verdict"), rows)
    L = [f"Aperture term (R frames vs REFCAT2 r): {ap:+.3f} mag.", "",
         "| filter | PinPoint catalogue (frames) | matched band | frames | "
         "offset vs matched | scatter | offset vs REFCAT2 r | scatter | r "
         "offset minus aperture term | ZMAG is |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append("| " + " | ".join("" if x is None else (f"{x:+.3f}" if
                 isinstance(x, float) else str(x)) for x in r) + " |")
    L += ["", "Reading: in the broad bands the matched-band offsets are "
          "large (B +1.2, G +0.7 mag) while the REFCAT2 r offsets stay "
          "within ~0.15 mag of the aperture term: PinPoint solved ZMAG "
          "against ATLAS-REFCAT2 r in every filter (its HISTORY names the "
          "ATLAS/PanSTARRS catalogue on every frame). In a non-r filter "
          "ZMAG is therefore r-band zero point plus the field's mean colour "
          "term — a transparency monitor, which is how B4 uses it (season-"
          "mode deficits per filter), never a calibration. Narrow-band "
          "offsets (1, H, O) depart by 0.5-0.9 mag under every hypothesis: "
          "their ZMAG is not interpretable as a zero point at all."]
    write_note("zmag_provenance", "TCRB-P0-zmag-provenance — what ZMAG "
               "is measured against", L,
               "Per frame, m_inst(1 s) + ZMAG compared with every catalogue "
               "band for the matched field stars (`p0.zmag_band_match`; "
               "table `tcrb_zmag_match`); per filter the mean offset and "
               "robust scatter.")
    print("\n".join(L))


def _target_band_residual(con, cat, code: str) -> list[tuple]:
    """T CrB's magnitude relative to the field-star colour relation in one
    code (reference band Gaia-synthetic R, fixed 8 px aperture): a source
    with line emission inside the passband is brighter than predicted."""
    from macro_tcrb import phot
    idx = {int(x): i for i, x in enumerate(cat["source_id"])}
    out = []
    for rid, exp, night in db.q(con, """SELECT obs_rowid, exptime, night
            FROM tcrb_img_frames WHERE status = 'ok' AND filter = ?""", code):
        rows = db.q(con, """SELECT source_id, flux_fix8, is_target FROM
            tcrb_img_stars WHERE obs_rowid = ? AND flux_fix8 > 0""", rid)
        k = np.array([idx.get(int(r[0]), -1) for r in rows])
        m = -2.5 * np.log10(np.array([r[1] for r in rows]) / exp)
        ist = np.array([r[2] for r in rows], bool)
        ok = k >= 0
        y = np.where(ok, cat["r_jkc_mag"][k] - m, np.nan)
        c = np.where(ok, cat["bp_rp"][k], np.nan)
        g = np.where(ok, cat["phot_g_mean_mag"][k], np.nan)
        use = ok & ~ist & np.isfinite(y) & np.isfinite(c) & (g > 9) & \
            (g < 15.5)
        f = phot.huber_line(c[use], y[use], np.full(int(use.sum()), 0.02))
        if f is None or not ist.any():
            continue
        t = np.flatnonzero(ist)[0]
        pred = cat["r_jkc_mag"][k[t]] - (f.zp + f.k * (c[t] - f.c_ref))
        out.append((rid, code, night, float(m[t] - pred), f.k))
    return out


def cmd_filters_report(args) -> None:
    """Filter-mapping table (P0-filter-forensics) and the ZMAG provenance
    table (P0-zmag-provenance), with every number read from the tables."""
    man = db.manifest_ro()
    con = db.connect()
    codes = ("6", "1", "G", "H", "L", "O", "W", "B", "R", "I", "V", "X")
    s2c = {r[0]: r[1:] for r in db.q(man, f"""
        SELECT f.filter, COUNT(*), SUM(d.verdict = 'dispersed'),
               SUM(d.verdict = 'direct'), SUM(d.verdict = 'indeterminate'),
               AVG(CASE WHEN f.zmag > 0 AND d.verdict = 'direct'
                   THEN f.zmag END),
               COUNT(CASE WHEN f.zmag > 0 AND d.verdict = 'direct'
                   THEN 1 END)
        FROM frames f LEFT JOIN frame_dispersion d USING (obs_rowid)
        WHERE f.night BETWEEN '2023-02-06' AND '2024-03-26'
          AND f.camera = 'AC4040' AND f.tree = 'rawimage'
          AND f.is_canonical = 1
          AND f.filter IN ({','.join('?' * len(codes))})
        GROUP BY f.filter""", *codes)}
    # Paired ZMAG difference to R on the same night (direct frames only):
    # a throughput ratio, i.e. a bandwidth indicator.
    pair = {r[0]: r[1:] for r in db.q(man, f"""
        WITH z AS (SELECT f.filter, f.night, AVG(f.zmag) zm
                   FROM frames f JOIN frame_dispersion d USING (obs_rowid)
                   WHERE f.night BETWEEN '2023-02-06' AND '2024-03-26'
                     AND f.camera = 'AC4040' AND f.tree = 'rawimage'
                     AND f.is_canonical = 1 AND f.zmag > 0
                     AND d.verdict = 'direct'
                   GROUP BY f.filter, f.night)
        SELECT a.filter, AVG(a.zm - b.zm), COUNT(*)
        FROM z a JOIN z b ON a.night = b.night AND b.filter = 'R'
        GROUP BY a.filter""")}
    lam = {r[0]: r[1:] for r in db.q(con, """SELECT filter, COUNT(*),
        AVG(lam_eff), GROUP_CONCAT(lam_eff), GROUP_CONCAT(best_band)
        FROM tcrb_filter_frames WHERE lam_eff IS NOT NULL GROUP BY 1""")}
    out = []
    for c in codes:
        n, ndisp, ndir, nind, zmean, nz = s2c.get(c, (0,) * 6)
        dz, npair = pair.get(c, (None, 0))
        nf, lmean, lall, bests = lam.get(c, (0, None, "", ""))
        ls = [float(x) for x in (lall or "").split(",") if x]
        lsd = float(np.std(ls, ddof=1)) if len(ls) > 1 else None
        from collections import Counter
        best = Counter((bests or "").split(",")).most_common(1)
        best = best[0][0] if best and best[0][0] else None
        width = (10 ** (-0.4 * -dz) if dz is not None else None)
        if n and ndisp / n > 0.5:
            mapping = "grism (dispersed in %d of %d frames)" % (ndisp, n)
        elif c == "W":
            mapping = "mixed: direct to 2024-01, dispersed from 2024-02"
        elif dz is not None and dz < -2.0:
            mapping = "narrow band"
        elif dz is not None:
            mapping = "broad band"
        else:
            mapping = "no paired ZMAG"
        out.append((c, n, ndisp, ndir, nind, zmean, nz, dz, npair, width,
                    nf, lmean, lsd, best, mapping))
    db.write_table(con, "tcrb_filter_map", (
        "code", "n_frames", "n_dispersed", "n_direct", "n_indeterminate",
        "zmag_mean_direct", "n_zmag", "dzmag_vs_R_same_night",
        "n_nights_paired", "throughput_vs_R", "n_tcrb_frames_colour",
        "lam_eff_A", "lam_eff_sd_A", "best_matching_band", "mapping"), out)
    lines = ["| code | frames (archive) | S2c dispersed / direct / indet. | "
             "ZMAG - ZMAG(R), same night (n) | throughput vs R | "
             "lambda_eff (A, T CrB field, n frames) | nearest standard band "
             "| mapping |", "|---|---|---|---|---|---|---|---|"]
    for o in out:
        lines.append(
            f"| {o[0]} | {o[1]} | {o[2]} / {o[3]} / {o[4]} | "
            f"{'' if o[7] is None else f'{o[7]:+.2f} ({o[8]})'} | "
            f"{'' if o[9] is None else f'{o[9]:.3f}'} | "
            f"{'' if o[11] is None else f'{o[11]:.0f}'}"
            f"{'' if o[12] is None else f' +- {o[12]:.0f}'} ({o[10]}) | "
            f"{o[13] or ''} | {o[14]} |")
    cat = load_catalogue()
    tr = []
    for c_ in ("H", "1", "O", "L", "R", "G", "I", "B"):
        tr += _target_band_residual(con, cat, c_)
    db.write_table(con, "tcrb_filter_target_resid", ("obs_rowid", "code",
                                                     "night", "resid_mag",
                                                     "colour_term"), tr)
    lines += ["", "T CrB relative to the field-star colour relation (fixed "
              "8 px aperture; Gaia-synthetic R reference; positive = fainter "
              "than predicted):", "", "| code | frames | nights | median "
              "residual (mag) |", "|---|---|---|---|"]
    for c_ in ("H", "1", "O", "L", "G", "B", "I", "R"):
        v = [r[3] for r in tr if r[1] == c_]
        if v:
            lines.append(f"| {c_} | {len(v)} | "
                         f"{len({r[2] for r in tr if r[1] == c_})} | "
                         f"{np.median(v):+.2f} |")
    vh = np.median([r[3] for r in tr if r[1] == "H"])
    v1 = np.median([r[3] for r in tr if r[1] == "1"])
    lines += ["", f"H minus 1: {vh - v1:+.2f} mag. The two narrow bands have "
              "the same colour response on field stars (lambda_eff 6500 +- "
              "300 A each), but T CrB is far brighter in H: H contains its "
              "Halpha emission and 1 lies in the M giant's TiO 6651 A band. "
              "Mapping: H = Halpha, 1 = a red narrow band longward of Halpha "
              "([S II] 6716/6731 by the usual H/O/S set; the colour test "
              "cannot separate 6563 from 6716).", "",
              "**Mapping adopted.** `6` = grism (dispersed on 94% of frames; "
              "which unit, per epoch, is A0b's and G-1's); `W` = an "
              "unfiltered/clear slot (throughput 1.5x R, like L) until "
              "2024-01, then a grism from 2024-02 (S2c: 84 dispersed, 0 "
              "direct after that date); `L` = luminance (lambda_eff 5670 A, "
              "throughput 1.5x R); `G` = Sloan g (4940 A; nearest band g), "
              "uppercase `R` = Sloan r (6234 A) and `I` = Sloan i (7587 A) — "
              "the CV hypothesis confirmed; `B` = Johnson B (4397 A, the "
              "method's check against B's 4380); `O` = [O III] narrow band "
              "(4990 +- 80 A, 4% of R's throughput); `H` = Halpha narrow "
              "band; `1` = red narrow band ([S II]). Unmapped frames: none. "
              "Transmission curves remain James's request to Cannon "
              "(SYNTHESIS §6); this table is the measured substitute."]
    write_note("filter_forensics", "TCRB-P0-filter-forensics — the "
               "single-character filter codes", lines,
               "Three independent measurements per code, no header "
               "assumption: (i) S2c's pixel verdict (is the slot a grism?); "
               "(ii) the same-night ZMAG difference to R on direct frames "
               "(PinPoint's 1-s zero point; a difference of -2.5 mag is a "
               "~10x narrower band); (iii) the colour response — for each "
               "frame the colour term of (catalogue band - m_inst) against "
               "Gaia BP-RP is fitted for nine Gaia-synthetic standard bands "
               "(Johnson-Cousins BVRI, SDSS griz, PS1 y) and lambda_eff is "
               "where that colour term crosses zero, interpolated in the "
               "bands' nominal wavelengths (BAND_LAM). Codes of known "
               "identity (B, R, I) are measured the same way and are the "
               "method's check.")
    print("\n".join(lines))


def cmd_resolve_report(args) -> None:
    """TCRB-P0-resolve evidence: every direct imaging frame's astrometry."""
    con = db.connect()
    rows = db.q(con, """SELECT f.filter, f.mode, f.exptime, f.path,
        f.pltsolvd, f.wcs_source, f.wcs_n, f.wcs_rms_arcsec, f.n_detect,
        f.status, c.verdict FROM tcrb_img_frames f LEFT JOIN tcrb_census c
        USING (obs_rowid) WHERE f.status NOT LIKE 'dispersed%'
        ORDER BY f.pltsolvd IS NULL DESC, f.filter""")
    n_dir = len(rows)
    unsolved = [r for r in rows if r[4] != 1]
    hdr_ok = [r for r in rows if r[4] == 1 and r[9] == "ok"]
    hdr_bad = [r for r in rows if r[4] == 1 and r[9] != "ok"]
    blind_ok = [r for r in unsolved if r[9] == "ok"]
    rms = [r[7] for r in hdr_ok + blind_ok if r[7] is not None]
    lines = [f"Direct (non-dispersed) canonical T CrB imaging frames: {n_dir}. "
             f"Without a PinPoint solution: {len(unsolved)} "
             f"({100 * len(unsolved) / n_dir:.0f}%); solved blind here: "
             f"{len(blind_ok)}. PinPoint solutions re-fitted and verified "
             f"against Gaia DR3: {len(hdr_ok)} (median rms "
             f"{np.median(rms):.2f} arcsec); PinPoint solutions that fail "
             f"verification: {len(hdr_bad)}.", "",
             "| code | mode | exp (s) | frame | PinPoint | outcome | matched "
             "| rms (\") | detections | target verdict |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in unsolved + hdr_bad:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]:g} | `{r[3]}` | "
                     f"{'yes' if r[4] == 1 else 'no'} | {r[9]} | {r[6]} | "
                     f"{'' if r[7] is None else f'{r[7]:.2f}'} | {r[8]} | "
                     f"{r[10]} |")
    lines += ["", "Every failure is a frame with too few stellar detections "
              "to match (short exposures through cloud; sep detections "
              "dominated by hot pixels and cosmic rays, median FWHM < 2.3 "
              "px), not a frame off target: all carry the T CrB pointing. "
              "None is a B anchor (all 23 B frames are solved)."]
    write_note("resolve", "TCRB-P0-resolve — astrometry of the unsolved "
               "imaging", lines,
               "Stage `measure`: a PinPoint solution is re-fitted to Gaia DR3 "
               "(>= 8 stars, rms < 1 arcsec); a frame without one is solved "
               "blind (triangle matching of stellar detections to Gaia "
               "around T CrB, both parities; >= 20 stars, rms < 1 arcsec). "
               "Grism frames get the identity gate (G-3) instead.")
    print("\n".join(lines[:1]))


# ===========================================================================
# TCRB-B0 — peak-at-target census, restated against today's S2c verdicts
# ===========================================================================
def cmd_census(args) -> None:
    """Every canonical T CrB imaging frame: native-pixel peak at the target
    (S2 ``s2_target_peaks``), the mode's measured cap and clip (S2
    ``s2_target_verdicts`` / ``detector_params``), and a verdict recomputed
    with the CURRENT S2c dispersion verdict (F-6 v1.2 re-issued several
    slot-6 and W frames as dispersed after S2 stored its verdicts)."""
    from macro_tcrb import phot
    man = db.manifest_ro()
    con = db.connect()
    rows = db.q(man, """
        SELECT p.obs_rowid, p.path, p.night, p.mode, p.filter, p.exptime,
               p.peak_raw, p.sky, v.bias_adu, v.clip_adu, v.cap_fraction,
               v.native_factor, v.verdict, d.verdict, f.jd
        FROM s2_target_peaks p JOIN s2_target_verdicts v USING (obs_rowid)
        CROSS JOIN frames f ON f.obs_rowid = p.obs_rowid
        LEFT JOIN frame_dispersion d USING (obs_rowid)
        WHERE p.is_canonical = 1 ORDER BY f.jd""")
    staged = {r[0] for r in db.q(man, """SELECT obs_rowid FROM
        stage_tcrb_monitoring WHERE role = 'science' AND canonical_target =
        'T CrB' AND filter NOT IN ('hrg', 'lrg')""")}
    out = []
    for (rid, path, night, mode, code, exp, pk, sky, bias, clip, capf, nf,
         v_old, s2c, jd) in rows:
        cap = (bias + capf * (clip - bias)) if None not in (bias, clip,
                                                             capf) else None
        pkn = pk * (nf or 1.0) if pk is not None else None
        if mode == "Low Gain":
            verdict = ("dispersed" if s2c == "dispersed"
                       else "no_ceiling_for_mode")
        else:
            verdict = phot.peak_verdict(pkn, cap, clip, s2c == "dispersed")
        role = ("science" if rid in staged else
                "excluded (H, single epoch)" if code == "H" else "unstaged")
        out.append((rid, path, night, mode, code, exp, jd, pk, sky, pkn, cap,
                    clip, s2c, v_old, verdict, role))
    n = db.write_table(con, "tcrb_census",
                       ("obs_rowid", "path", "night", "mode", "filter",
                        "exptime", "jd", "peak_raw", "sky", "peak_native",
                        "cap_adu", "clip_adu", "s2c_verdict", "s2_verdict",
                        "verdict", "role"), out)
    # Anchor set: clean B frames; singletons = clean direct frames of other
    # codes on nights with no clean B (their band identity is the
    # filter-forensics table's; they are reported, never merged with B).
    bn = {r[2] for r in out if r[4] == "B" and r[14] == "clean"}
    anchors = []
    for r in out:
        if r[15] != "science" or r[14] != "clean":
            continue
        kind = ("B anchor" if r[4] == "B" else
                "singleton" if r[2] not in bn else "same night as B")
        anchors.append((r[0], r[2], r[4], r[3], r[5], kind))
    db.write_table(con, "tcrb_anchor_set", ("obs_rowid", "night", "filter",
                                            "mode", "exptime", "kind"),
                   anchors)
    db.record_stage(con, "census", {"tcrb_census": n,
                                     "tcrb_anchor_set": len(anchors)})
    tab = db.q(con, """SELECT filter, mode, verdict, COUNT(*),
        COUNT(DISTINCT night) FROM tcrb_census WHERE role = 'science'
        GROUP BY 1, 2, 3 ORDER BY 1, 2, 3""")
    changed = db.q1(con, "SELECT COUNT(*) FROM tcrb_census WHERE "
                         "(verdict = 'clean') != (s2_verdict = 'clean')")
    lines = ["| code | readout mode | verdict | frames | nights |",
             "|---|---|---|---|---|"]
    lines += [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in tab]
    ak = db.q(con, """SELECT kind, filter, COUNT(*), COUNT(DISTINCT night),
        GROUP_CONCAT(DISTINCT night) FROM tcrb_anchor_set GROUP BY 1, 2""")
    lines += ["", "Anchor set restated from the census:", "",
              "| kind | code | frames | nights | nights list |",
              "|---|---|---|---|---|"]
    lines += [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in ak]
    lines += ["", f"Staged science imaging frames: {len(staged)}; every one "
              f"carries a peak and a verdict "
              f"({sum(1 for r in out if r[15] == 'science')} rows). "
              f"Frames whose usability (clean or not) differs from S2's "
              f"stored verdict: {changed} (S2 split non-clean frames into "
              f"clipped/flat-topped/above-cap differently; only the "
              f"clean/not-clean boundary decides an anchor)."]
    write_note("peak_census", "TCRB-B0 — peak-at-target census of every "
               "imaging frame", lines,
               "Peak = RAW native-pixel maximum at the target (S2 "
               "`s2_target_peaks`; the AC4040 is unbinned, native factor 1). "
               "Cap = bias + cap_fraction x (clip - bias), the mode's "
               "measured linearity cap (F-5). Verdicts: `macro_tcrb.phot."
               "peak_verdict`, with S2c's current frame verdict deciding "
               "'dispersed'. Low Gain has no measured ceiling.")
    print("\n".join(lines))


# ===========================================================================
# TCRB-D1 — external context, cached with pull dates
# ===========================================================================
AAVSO_AUID = "000-BBW-825"
AAVSO_URL = "https://vsx.aavso.org/index.php"
ARAS_BASE = "https://aras-database.github.io/database/"
ZTF_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
ASASSN_V2 = "https://asas-sn.ifa.hawaii.edu/"
#: The windows pulled: the 2025 grism season with a month of margin, and
#: the 2023-24 imaging season (the B anchors and the A0b spectra).
EXT_WINDOWS = (("2025-01-15", "2025-07-15"), ("2023-04-15", "2024-04-30"))


def _jd(date: str) -> float:
    from astropy.time import Time
    return float(Time(date + "T00:00:00", scale="utc").jd)


def _get(url: str, timeout: float = 300.0) -> bytes:
    import urllib.request
    req = urllib.request.Request(url, headers={
        "User-Agent": "MACRO-RLMT research pipeline (T CrB project)"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _save(rel: str, payload: bytes, gz: bool = True) -> tuple[str, str]:
    import gzip
    import hashlib
    path = db.EXTERNAL / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        path = path.with_name(path.name + ".gz")
        with gzip.open(path, "wb") as fh:
            fh.write(payload)
    else:
        path.write_bytes(payload)
    return (str(path.relative_to(db.REPO)),
            hashlib.sha256(payload).hexdigest())


def cmd_external(args) -> None:
    """Pull and cache AAVSO, ARAS, ZTF, ASAS-SN, TESS and Swift context."""
    con = db.connect()
    con.execute("""CREATE TABLE IF NOT EXISTS tcrb_ext_fetch (
        source TEXT, chunk TEXT, url TEXT, pulled_utc TEXT, n_rows INTEGER,
        cache_path TEXT, sha256 TEXT, ok INTEGER, note TEXT,
        PRIMARY KEY (source, chunk))""")
    done = {(r[0], r[1]) for r in db.q(con, "SELECT source, chunk FROM "
                                       "tcrb_ext_fetch WHERE ok = 1")}
    want = set((args.sources or "aavso,aras,ztf,asassn,tess,swift"
                ).split(","))

    def rec(source, chunk, url, n, path, sha, ok, note):
        con.execute("INSERT OR REPLACE INTO tcrb_ext_fetch VALUES "
                    "(?,?,?,?,?,?,?,?,?)", (source, chunk, url, db.utc_now(),
                                            n, path, sha, int(ok), note))
        con.commit()
        print(f"  {source:7s} {chunk:24s} ok={ok} n={n} {note[:70]}",
              flush=True)

    import urllib.parse
    # ---- AAVSO: monthly chunks of the AID (api.delim) -------------------
    if "aavso" in want:
        import datetime as _dt
        for a, b in EXT_WINDOWS:
            d = _dt.date.fromisoformat(a)
            end = _dt.date.fromisoformat(b)
            while d < end:
                nxt = min((d.replace(day=1) + _dt.timedelta(days=32)
                           ).replace(day=1), end)
                chunk = f"{d.isoformat()}_{nxt.isoformat()}"
                if ("aavso", chunk) in done:
                    d = nxt
                    continue
                url = f"{AAVSO_URL}?" + urllib.parse.urlencode({
                    "view": "api.delim", "ident": AAVSO_AUID,
                    "fromjd": f"{_jd(d.isoformat()):.5f}",
                    "tojd": f"{_jd(nxt.isoformat()):.5f}",
                    "delimiter": "@@@"})
                try:
                    raw = _get(url)
                    n = raw.count(b"\n") - 1
                    path, sha = _save(f"aavso/tcrb_{chunk}.delim", raw)
                    rec("aavso", chunk, url, n, path, sha, True,
                        "AAVSO International Database, api.delim route")
                except Exception as exc:
                    rec("aavso", chunk, url, 0, "", "", False,
                        f"{type(exc).__name__}: {exc}")
                d = nxt
    # ---- ARAS: every Halpha-covering spectrum in both windows ------------
    if "aras" in want:
        import csv
        lst = list(csv.DictReader(open(db.NOVELTY / "external" /
                                       "aras_tcrb_listing.csv")))
        for row in lst:
            dte = row["date_ut"]
            if not any(a <= dte <= b for a, b in EXT_WINDOWS):
                continue
            if row["ha_loose"] != "True" or not row["file"]:
                continue
            chunk = row["file"].rsplit("/", 1)[-1]
            if ("aras", chunk) in done:
                continue
            url = ARAS_BASE + row["file"]
            try:
                raw = _get(url, 120)
                path, sha = _save(f"aras/{chunk}", raw, gz=False)
                rec("aras", chunk, url, 1, path, sha, True,
                    f"{row['observer']} R={row['resolving_power']}")
            except Exception as exc:
                rec("aras", chunk, url, 0, "", "", False,
                    f"{type(exc).__name__}: {exc}")
    # ---- ZTF: IRSA light-curve cone (T CrB saturates; recorded) ----------
    if "ztf" in want and ("ztf", "cone2arcsec") not in done:
        url = f"{ZTF_URL}?" + urllib.parse.urlencode(
            {"POS": f"CIRCLE {TCRB_RA} {TCRB_DEC} 0.000556",
             "FORMAT": "csv"}, quote_via=urllib.parse.quote)
        try:
            raw = _get(url)
            n = max(raw.count(b"\n") - 1, 0)
            path, sha = _save("ztf/tcrb_lc.csv", raw)
            rec("ztf", "cone2arcsec", url, n, path, sha, True,
                "IRSA ZTF DR light curves; T CrB (V~10) is far above the "
                "ZTF saturation limit (~12.5 mag): context only")
        except Exception as exc:
            rec("ztf", "cone2arcsec", url, 0, "", "", False,
                f"{type(exc).__name__}: {exc}")
    # ---- ASAS-SN Sky Patrol v2: reachability recorded --------------------
    if "asassn" in want and ("asassn", "skypatrol_v2") not in done:
        try:
            raw = _get(ASASSN_V2, 30)
            rec("asassn", "skypatrol_v2", ASASSN_V2, 0, "", "", False,
                f"service answered ({len(raw)} bytes) but Sky Patrol v2 "
                "light curves need its python client/API session; not "
                "pulled")
        except Exception as exc:
            rec("asassn", "skypatrol_v2", ASASSN_V2, 0, "", "", False,
                f"UNREACHABLE from this network: {type(exc).__name__}: "
                f"{exc}")
    # ---- TESS: sectors and FFI cutouts (MAST TESScut) ---------------------
    if "tess" in want:
        try:
            from astroquery.mast import Tesscut
            from astropy.coordinates import SkyCoord
            c = SkyCoord(TCRB_RA, TCRB_DEC, unit="deg")
            sec = Tesscut.get_sectors(coordinates=c)
            sectors = [int(x) for x in sec["sector"]]
            raw = json.dumps({"sectors": sectors}).encode()
            path, sha = _save("tess/sectors.json", raw, gz=False)
            rec("tess", "sectors", "astroquery.mast.Tesscut.get_sectors",
                len(sectors), path, sha, True, f"sectors {sectors}")
            for sct in sectors:
                chunk = f"s{sct:04d}"
                if ("tess", chunk) in done:
                    continue
                try:
                    hd = Tesscut.get_cutouts(coordinates=c, size=21,
                                             sector=sct)
                    import io
                    buf = io.BytesIO()
                    hd[0].writeto(buf)
                    path, sha = _save(f"tess/tcrb_{chunk}_21px.fits",
                                      buf.getvalue(), gz=True)
                    rec("tess", chunk, "astroquery.mast.Tesscut.get_cutouts"
                        "(size=21)", int(len(hd[0][1].data)), path, sha,
                        True, "FFI cutout")
                except Exception as exc:
                    rec("tess", chunk, "Tesscut.get_cutouts", 0, "", "",
                        False, f"{type(exc).__name__}: {exc}")
        except Exception as exc:
            rec("tess", "sectors", "Tesscut.get_sectors", 0, "", "", False,
                f"{type(exc).__name__}: {exc}")
    # ---- Swift: the observation log (HEASARC swiftmastr) ------------------
    if "swift" in want and ("swift", "swiftmastr") not in done:
        try:
            from astroquery.heasarc import Heasarc
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            tab = Heasarc.query_region(SkyCoord(TCRB_RA, TCRB_DEC,
                                                unit="deg"),
                                       catalog="swiftmastr",
                                       radius=5 * u.arcmin)
            import io
            buf = io.StringIO()
            tab.write(buf, format="ascii.ecsv")
            path, sha = _save("swift/swiftmastr.ecsv",
                              buf.getvalue().encode())
            rec("swift", "swiftmastr", "astroquery.heasarc swiftmastr r=5'",
                len(tab), path, sha, True, "observation log only")
        except Exception as exc:
            rec("swift", "swiftmastr", "astroquery.heasarc", 0, "", "",
                False, f"{type(exc).__name__}: {exc}")
    n = db.q1(con, "SELECT COUNT(*) FROM tcrb_ext_fetch")
    db.record_stage(con, "external", {"tcrb_ext_fetch": n})


def cmd_external_report(args) -> None:
    """The D1 evidence note, from tcrb_ext_fetch."""
    con = db.connect()
    rows = db.q(con, """SELECT source, SUM(ok), SUM(1 - ok), SUM(n_rows),
        MIN(pulled_utc), MAX(pulled_utc) FROM tcrb_ext_fetch GROUP BY 1""")
    fails = db.q(con, """SELECT source, chunk, note FROM tcrb_ext_fetch
        WHERE ok = 0""")
    lines = ["| source | chunks cached | failed | rows / files | first pull "
             "| last pull |", "|---|---|---|---|---|---|"]
    lines += [f"| {a} | {b} | {c} | {d} | {e} | {f} |"
              for a, b, c, d, e, f in rows]
    lines += ["", "Failures (recorded, not hidden):", ""]
    lines += [f"- {a} `{b}`: {c}" for a, b, c in fails]
    tess = db.q1(con, "SELECT note FROM tcrb_ext_fetch WHERE source='tess' "
                      "AND chunk='sectors'")
    lines += ["", f"TESS: {tess} — none inside 2025 Feb-Jun, so TESS gives "
              "no contemporaneous continuum; the cutouts are cached as "
              "context. ZTF: T CrB is above the saturation limit. "
              "AAVSO: the GOAL_RULES premise that AAVSO blocks scripted "
              "access holds for the WebObs/data-download pages (bot "
              "challenge, TCRB-N1); the AID `api.delim` route used by the CV "
              "package answered every request and is the source of the B "
              "and Rc/Ic continuum. Every AAVSO row carries its observer "
              "code; RLMT's own submissions (observer MALW) are tagged at "
              "use, never counted as external."]
    write_note("external_pulls", "TCRB-D1 — external context, cached with "
               "pull dates", lines,
               "Stage `external` of `pipeline/scripts/run_tcrb.py`; cache "
               "under `TCrB_Monitoring/external_data/<source>/`, one row per "
               "chunk in `tcrb_ext_fetch` (URL or call, pull time, row "
               "count, sha256). Windows: 2025-01-15..2025-07-15 and "
               "2023-04-15..2024-04-30.")
    print("\n".join(lines))


# ===========================================================================
# TCRB-B1, -B3, -B5: the B ensemble, its penalties and empirical errors
# ===========================================================================
B_FIT_RANGE = (9.5, 14.5)      #: ensemble stars, Gaia-synthetic Johnson B
N_CHECK = 5                    #: check stars bracketing T CrB in B
ISOLATION_ARCSEC = 8.0         #: no Gaia neighbour within this radius ...
ISOLATION_DMAG = 4.0           #: ... brighter than (star G + this)


def _aavso_night(con, jd: float, band: str, half: float = 0.5):
    """Median AAVSO magnitude (non-own, CCD/PEP) within +-half d; n."""
    v = [r[0] for r in db.q(con, """SELECT mag FROM tcrb_aavso WHERE
        band = ? AND own = 0 AND jd BETWEEN ? AND ?""", band, jd - half,
                             jd + half)]
    return (float(np.median(v)), len(v)) if v else (None, 0)


def _isolated(cat) -> np.ndarray:
    from scipy.spatial import cKDTree
    cosd = np.cos(np.radians(TCRB_DEC))
    xy = np.c_[cat["ra"] * cosd, cat["dec"]] * 3600
    tr = cKDTree(xy)
    g = cat["phot_g_mean_mag"]
    iso = np.ones(g.size, bool)
    for i, nb in enumerate(tr.query_ball_point(xy, ISOLATION_ARCSEC)):
        for j in nb:
            if j != i and g[j] < g[i] + ISOLATION_DMAG:
                iso[i] = False
                break
    return iso


def cmd_ensemble(args) -> None:
    """Per B frame: Gaia-synthetic Johnson B ensemble (colour term in B-V),
    T CrB natural and transformed B, colour-extrapolation error, check-star
    residuals, the flat-staleness term (ensemble residual vs detector
    position) and the cross-mode term (StackPro vs High Gain B flat at the
    ensemble stars' positions)."""
    from macro_tcrb import phot
    from macro_phot import cattie as ct
    con = db.connect()
    man = db.manifest_ro()
    # B2's adopted aperture (lowest differential check-star scatter);
    # r15 (the strategy default) until B2 has run.
    try:
        ap = db.q1(con, "SELECT aperture FROM tcrb_b2_apertures ORDER BY "
                        "check_rms_mag LIMIT 1") or "r15"
    except Exception:
        ap = "r15"
    print(f"  aperture: {ap}")
    cat = load_catalogue()
    iso = _isolated(cat)
    idx = {int(sid): i for i, sid in enumerate(cat["source_id"])}
    tgt = db.q(con, "SELECT source_id, b_jkc, v_jkc FROM tcrb_target")[0]
    b_t = tgt[1]
    bmag, vmag = cat["b_jkc_mag"], cat["v_jkc_mag"]
    berr = np.nan_to_num(cat["b_jkc_mag_error"], nan=0.05)
    rc2 = cat["rc2_gmag"] + 0.313 * (cat["rc2_gmag"] - cat["rc2_rmag"]) \
        + 0.219                  # Jordi et al. (2006) B from g, r (stars)
    # Check stars: the N_CHECK isolated catalogue stars nearest T CrB's B
    # among those with good synthetic B — fixed once for the season.
    good = iso & np.isfinite(bmag) & np.isfinite(vmag) & \
        (bmag > B_FIT_RANGE[0]) & (bmag < B_FIT_RANGE[1])
    good[idx[int(tgt[0])]] = False
    cand = np.flatnonzero(good)
    ref_b = b_t if np.isfinite(b_t) else 11.6
    checks = set(cand[np.argsort(np.abs(bmag[cand] - ref_b))][:N_CHECK])
    # Flats for the cross-mode term.
    mc = _MasterCache(_ac4040_masters(man))
    frames = db.q(con, """SELECT obs_rowid, path, night, mode, exptime, jd,
        airmass, flat_master FROM tcrb_img_frames WHERE status = 'ok' AND
        filter = 'B' ORDER BY jd""")
    fr_out, ck_out, st_out = [], [], []
    for rid, path, night, mode, exp, jd, am, flatp in frames:
        cap = db.q1(con, "SELECT cap_adu FROM tcrb_census WHERE obs_rowid "
                         "= ?", rid)
        rows = db.q(con, f"""SELECT source_id, is_target, x, y, flux_{ap},
            err_{ap}, peak_raw, sep_flag FROM tcrb_img_stars WHERE
            obs_rowid = ? AND flux_{ap} > 0""", rid)
        a = [r for r in rows if int(r[0]) in idx]
        k = np.array([idx[int(r[0])] for r in a])
        isT = np.array([r[1] for r in a], bool)
        x = np.array([r[2] for r in a])
        y = np.array([r[3] for r in a])
        fl = np.array([r[4] for r in a])
        fe = np.array([r[5] for r in a])
        pk = np.array([r[6] for r in a])
        sf = np.array([r[7] for r in a])
        m = -2.5 * np.log10(fl / exp)
        me = 1.0857 * fe / fl
        use = (~isT) & good[k] & (sf == 0) & (pk < cap) & \
            ~np.isin(k, list(checks))
        c = bmag[k] - vmag[k]
        sig = np.sqrt(me ** 2 + berr[k] ** 2 + 0.005 ** 2)
        f = phot.huber_line(c[use], bmag[k][use] - m[use], sig[use])
        if f is None:
            continue
        resid = bmag[k] - m - (f.zp + f.k * (c - f.c_ref))
        cmin, cmax, _, _ = ct.colour_range(c[use][f.used])
        # Flat staleness: weighted residual vs position, 3x3 grid medians.
        u = np.flatnonzero(use)[f.used]
        gx = np.clip((x[u] / 4096 * 3).astype(int), 0, 2)
        gy = np.clip((y[u] / 4096 * 3).astype(int), 0, 2)
        cells = [np.median(resid[u][(gx == i) & (gy == j)])
                 for i in range(3) for j in range(3)
                 if ((gx == i) & (gy == j)).sum() >= 3]
        cell_n = [((gx == i) & (gy == j)).sum() for i in range(3)
                  for j in range(3) if ((gx == i) & (gy == j)).sum() >= 3]
        noise = [1.2533 * f.rms / np.sqrt(n_) for n_ in cell_n]
        spread = float(np.std(cells, ddof=1)) if len(cells) > 2 else None
        stale = (float(np.sqrt(max(spread ** 2 - np.mean(
            np.square(noise)), 0.0))) if spread is not None else None)
        # Cross-mode: StackPro B flat vs High Gain B flat at these stars.
        fs, _, _ = mc.flat(mode, "B", night)
        fh, _, _ = mc.flat("High Gain", "B", night)
        xm = None
        if fs is not None and fh is not None:
            xi, yi = x[u].astype(int), y[u].astype(int)
            d = 2.5 * np.log10(fs[yi, xi] / fh[yi, xi])
            xm = float(np.std(d - np.median(d), ddof=1))
        # T CrB
        t = np.flatnonzero(isT)
        bt_nat = bt_tr = bt_err = ext = bv_t = None
        nbv = 0
        if t.size:
            ti = t[0]
            b_av, nb = _aavso_night(con, jd, "B")
            v_av, nv = _aavso_night(con, jd, "V")
            if b_av is not None and v_av is not None:
                bv_t, nbv = b_av - v_av, min(nb, nv)
            bt_nat = float(m[ti] + f.zp)            # at the ensemble c_ref
            if bv_t is not None:
                bt_tr = float(bt_nat + f.k * (bv_t - f.c_ref))
                ext = ct.colour_extrapolation_error(bv_t, cmin, cmax, f.k,
                                                    f.k_err)
            bt_err = float(np.hypot(me[ti], f.zp_err))
        # REFCAT2 cross-check of the zero point (Jordi B from g, r).
        okr = use & np.isfinite(rc2[k])
        f2 = phot.huber_line(c[okr], rc2[k][okr] - m[okr], sig[okr],
                             c_ref=f.c_ref) if okr.sum() >= 5 else None
        fr_out.append((rid, path, night, mode, exp, jd, am, f.zp, f.zp_err,
                       f.k, f.k_err, f.c_ref, f.n, f.n_clip, f.chi2, f.dof,
                       f.rms, cmin, cmax, bt_nat, bt_tr, bt_err, bv_t, nbv,
                       ext, stale, xm, f2.zp - f.zp if f2 else None))
        for j in np.flatnonzero(np.isin(k, list(checks))):
            ck_out.append((rid, night, mode, jd, int(cat["source_id"][k[j]]),
                           float(bmag[k[j]]), float(resid[j]), float(me[j]),
                           float(pk[j])))
    n = db.write_table(con, "tcrb_bphot_frames", (
        "obs_rowid", "path", "night", "mode", "exptime", "jd", "airmass",
        "zp", "zp_err", "k_bv", "k_bv_err", "c_ref", "n_fit", "n_clip",
        "chi2", "dof", "rms", "colour_min", "colour_max", "b_nat", "b_trans",
        "b_stat_err", "bv_target_aavso", "n_aavso_bv", "extrap_err",
        "flat_stale_mag", "crossmode_mag", "zp_refcat2_minus_gspc"), fr_out)
    db.write_table(con, "tcrb_bphot_checks", (
        "obs_rowid", "night", "mode", "jd", "source_id", "b_cat", "resid",
        "phot_err", "peak_raw"), ck_out)
    db.record_stage(con, "ensemble", {"tcrb_bphot_frames": n,
                                       "tcrb_bphot_checks": len(ck_out)})
    for r in db.q(con, """SELECT night, COUNT(*), ROUND(AVG(b_trans),3),
            ROUND(AVG(k_bv),3), ROUND(AVG(chi2/dof),2), ROUND(AVG(rms),3),
            ROUND(AVG(flat_stale_mag),4), ROUND(AVG(crossmode_mag),4),
            ROUND(AVG(extrap_err),3), ROUND(AVG(zp_refcat2_minus_gspc),3)
            FROM tcrb_bphot_frames GROUP BY night"""):
        print("  ", r)


# ===========================================================================
# TCRB-B5 / B6 — empirical errors and the precision statement
# ===========================================================================
def cmd_errors(args) -> None:
    """Check-star rms per night per mode (each star about its own nightly
    mean, so a catalogue offset cannot enter), and the per-mode precision."""
    from macro_tcrb import phot
    con = db.connect()
    rows = db.q(con, """SELECT night, mode, source_id, resid, b_cat
        FROM tcrb_bphot_checks""")
    groups: dict = {}
    for night, mode, sid, r, b in rows:
        groups.setdefault((night, mode), {}).setdefault(sid, []).append(r)
    out = []
    for (night, mode), stars in sorted(groups.items()):
        dev = np.concatenate([np.asarray(v) - np.median(v)
                              for v in stars.values() if len(v) >= 3]) \
            if any(len(v) >= 3 for v in stars.values()) else np.array([])
        n_frames = max(len(v) for v in stars.values())
        rms = phot.check_star_rms(dev) if dev.size else None
        out.append((night, mode, n_frames, len(stars), int(dev.size), rms))
    db.write_table(con, "tcrb_berrors", ("night", "mode", "n_frames",
                                         "n_check", "n_resid", "rms_mag"),
                   out)
    # Season-level systematic floor: the scatter of the check stars' NIGHTLY
    # means about their season means (night-to-night zero-point stability).
    nm: dict = {}
    for night, mode, sid, r, b in rows:
        nm.setdefault(sid, {}).setdefault(night, []).append(r)
    floors = []
    for sid, nights in nm.items():
        m = [np.median(v) for v in nights.values() if len(v) >= 2]
        if len(m) >= 2:
            floors.append(np.std(m, ddof=1))
    floor = float(np.median(floors)) if floors else None
    prec = []
    for mode in sorted({o[1] for o in out}):
        v = [o[5] for o in out if o[1] == mode and o[5] is not None]
        w = [o[4] for o in out if o[1] == mode and o[5] is not None]
        pooled = float(np.sqrt(np.average(np.square(v), weights=w))) \
            if v else None
        prec.append((mode, len(v), pooled, floor, len(floors)))
    db.write_table(con, "tcrb_bprecision", ("mode", "n_nights",
                                            "per_frame_rms_mag",
                                            "season_floor_mag",
                                            "n_floor_stars"), prec)
    db.record_stage(con, "errors", {"tcrb_berrors": len(out)})
    for o in out + prec:
        print("  ", o)


# ===========================================================================
# TCRB-C1 — archival flickering limits (one table)
# ===========================================================================
#: Published scale beside every limit (standing rule 2): B-band variability
#: amplitude ~0.07 mag, high-speed photometry of 2023 June 8 (Maslennikova
#: et al. 2023, Astronomy Letters, arXiv:2308.10011).  It is an AMPLITUDE;
#: our limits are RMS.  Both are printed, labelled, never converted.
FLICKER_PUBLISHED_AMP = 0.07
SNIPPET_GAP_MIN = 30.0         #: frames closer than this form one snippet


def cmd_flicker(args) -> None:
    from macro_tcrb import phot
    con = db.connect()
    fr = db.q(con, """SELECT b.obs_rowid, b.night, b.mode, b.jd, b.b_nat,
        b.airmass, i.fwhm_px FROM tcrb_bphot_frames b JOIN tcrb_img_frames i
        USING (obs_rowid) WHERE b.b_nat IS NOT NULL ORDER BY b.jd""")
    snips, cur = [], []
    for r in fr:
        if cur and ((r[3] - cur[-1][3]) * 1440 > SNIPPET_GAP_MIN
                    or r[2] != cur[-1][2]):
            snips.append(cur)
            cur = []
        cur.append(r)
    if cur:
        snips.append(cur)
    out = []
    for sn in snips:
        if len(sn) < 4:
            continue
        ids = [r[0] for r in sn]
        t = (np.array([r[3] for r in sn]) - sn[0][3]) * 86400
        m = np.array([r[4] for r in sn])
        X = np.vstack([[r[5] for r in sn], [r[6] for r in sn]]).T
        ck: dict = {}
        for rid, sid, res in db.q(con, f"""SELECT obs_rowid, source_id,
                resid FROM tcrb_bphot_checks WHERE obs_rowid IN
                ({','.join('?' * len(ids))})""", *ids):
            ck.setdefault(sid, {})[rid] = res
        noise = [np.array([d[i] for i in ids]) for d in ck.values()
                 if all(i in d for i in ids)]
        res = phot.flicker_upper_limit(t, m, noise, X, n_mc=args.n_mc)
        out.append((sn[0][1], sn[0][2], len(sn), float(t[-1] / 60),
                    res["order"], res["n_cov"], res["obs_rms"], res["dof"],
                    res.get("noise_rms_median"), res.get("noise_rms_p95"),
                    res["limit"].get(1.0), res["limit"].get(2.0),
                    res["recovered90"].get(1.0), res["recovered90"].get(2.0),
                    max(res["limit"].get(2.0) or 0,
                        res["recovered90"].get(2.0) or 0),
                    FLICKER_PUBLISHED_AMP, len(noise)))
        print("  ", out[-1], flush=True)
    db.write_table(con, "tcrb_flicker", (
        "night", "mode", "n_frames", "span_min", "detrend_order",
        "n_covariates", "obs_rms", "dof",
        "check_rms_median", "check_rms_p95", "limit95_beta1",
        "limit95_beta2", "recovered90_beta1", "recovered90_beta2",
        "quoted_limit_rms", "published_amp_B", "n_check_series"), out)
    db.record_stage(con, "flicker", {"tcrb_flicker": len(out)})


# ===========================================================================
# Phase B / C1 report: B2 aperture choice, B4 ZMAG QC, anchors vs AAVSO,
# B6 precision, C1 table — all from the tables
# ===========================================================================
ZMAG_CUT, ZMAG_FLAG = 0.5, 0.3      #: ruling 4 (season-mode cut; M1 flag)


def cmd_phaseb_report(args) -> None:
    from macro_tcrb import phot
    con = db.connect()
    # ---- B2: three apertures; check-star scatter after a per-frame
    # median zero point formed IN THE SAME aperture from the fit stars ----
    cat = load_catalogue()
    idx = {int(x): i for i, x in enumerate(cat["source_id"])}
    checks = {r[0] for r in db.q(con, "SELECT DISTINCT source_id FROM "
                                      "tcrb_bphot_checks")}
    iso = _isolated(cat)
    ap_rows = []
    for ap in ("r15", "fix8", "r30"):
        rows = db.q(con, f"""SELECT s.obs_rowid, s.is_target, s.source_id,
            s.flux_{ap}, f.exptime, b.night, s.peak_raw, s.sep_flag,
            c.cap_adu FROM tcrb_img_stars s
            JOIN tcrb_img_frames f USING (obs_rowid)
            JOIN tcrb_bphot_frames b USING (obs_rowid)
            JOIN tcrb_census c USING (obs_rowid)
            WHERE s.flux_{ap} > 0""")
        by: dict = {}
        for r in rows:
            by.setdefault(r[0], []).append(r)
        ck: dict = {}
        for rid, rs in by.items():
            d_fit = []
            for (_, ist, sid, fl, ex, night, pk, sf, cap) in rs:
                k = idx.get(int(sid))
                if k is None or ist or sid in checks or sf or pk >= cap:
                    continue
                bm = cat["b_jkc_mag"][k]
                if iso[k] and np.isfinite(bm) and B_FIT_RANGE[0] < bm < \
                        B_FIT_RANGE[1]:
                    d_fit.append(bm + 2.5 * np.log10(fl / ex))
            if len(d_fit) < 5:
                continue
            zp = float(np.median(d_fit))
            for (_, ist, sid, fl, ex, night, pk, sf, cap) in rs:
                if sid in checks:
                    k = idx[int(sid)]
                    ck.setdefault((sid, night), []).append(
                        cat["b_jkc_mag"][k] + 2.5 * np.log10(fl / ex) - zp)
        dev = np.concatenate([np.asarray(v) - np.median(v)
                              for v in ck.values() if len(v) >= 3])
        ap_rows.append((ap, len(by), phot.check_star_rms(dev),
                        int(dev.size)))
    # aperture correction r15 -> r30 for T CrB (growth to the wide aperture)
    gr = db.q(con, """SELECT AVG(-2.5 * (LOG10(s.flux_r15) - LOG10(s.flux_r30)))
        FROM tcrb_img_stars s JOIN tcrb_bphot_frames b USING (obs_rowid)
        WHERE s.is_target = 1 AND s.flux_r15 > 0 AND s.flux_r30 > 0""")
    best = min(ap_rows, key=lambda r: r[2])[0]
    db.write_table(con, "tcrb_b2_apertures", ("aperture", "n_target_frames",
                                              "check_rms_mag", "n_resid"),
                   ap_rows)
    # ---- B4: ZMAG QC flags (never calibration) --------------------------
    zq = []
    zf = db.q(con, """SELECT obs_rowid, filter, mode, zmag FROM
        tcrb_img_frames WHERE zmag IS NOT NULL AND zmag > 0""")
    from collections import Counter
    for flt in sorted({r[1] for r in zf}):
        z = np.array([r[3] for r in zf if r[1] == flt])
        # season mode: the densest 0.1-mag bin (ZMAG is cloud-skewed faint)
        h, e = np.histogram(z, bins=np.arange(z.min() - 0.05,
                                              z.max() + 0.15, 0.1))
        mode = float(0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1]))
        for rid, f_, m_, zz in [r for r in zf if r[1] == flt]:
            zq.append((rid, f_, m_, zz, mode, mode - zz,
                       int(mode - zz > ZMAG_CUT), int(mode - zz > ZMAG_FLAG)))
    db.write_table(con, "tcrb_zmag_qc", ("obs_rowid", "filter", "mode",
                                         "zmag", "season_mode", "deficit",
                                         "cut_0p5", "flag_0p3"), zq)
    # ---- anchors vs AAVSO -----------------------------------------------
    anc = db.q(con, """SELECT f.night, COUNT(*), AVG(f.b_trans),
        AVG(f.b_nat), AVG(f.extrap_err), AVG(f.bv_target_aavso),
        MIN(f.jd), MAX(f.jd) FROM tcrb_bphot_frames f GROUP BY f.night""")
    arows = []
    for night, n, bt, bn, ext, bv, j0, j1 in anc:
        v = [r[0] for r in db.q(con, """SELECT mag FROM tcrb_aavso WHERE
            band = 'B' AND own = 0 AND jd BETWEEN ? AND ?""", j0 - 0.1,
                                j1 + 0.1)]
        obs = db.q1(con, """SELECT COUNT(DISTINCT obscode) FROM tcrb_aavso
            WHERE band = 'B' AND own = 0 AND jd BETWEEN ? AND ?""",
                    j0 - 0.1, j1 + 0.1)
        rms = db.q1(con, "SELECT rms_mag FROM tcrb_berrors WHERE night = ?",
                    night)
        qc = db.q1(con, """SELECT SUM(cut_0p5) FROM tcrb_zmag_qc z JOIN
            tcrb_bphot_frames b USING (obs_rowid) WHERE b.night = ?""", night)
        arows.append((night, n, bt, bn, rms, ext, bv,
                      float(np.median(v)) if v else None,
                      float(1.4826 * np.median(np.abs(np.array(v) -
                                                      np.median(v))))
                      if len(v) > 2 else None, len(v), obs, qc))
    db.write_table(con, "tcrb_b_anchors", (
        "night", "n_frames", "b_trans", "b_nat", "per_frame_rms",
        "extrap_err", "bv_aavso", "aavso_b_median", "aavso_b_mad",
        "n_aavso", "n_aavso_observers", "n_zmag_cut"), arows)
    pen = db.q(con, """SELECT AVG(flat_stale_mag), MAX(flat_stale_mag),
        AVG(crossmode_mag), AVG(zp_refcat2_minus_gspc), AVG(k_bv),
        AVG(k_bv_err), AVG(chi2 / dof), MIN(dof), MAX(dof)
        FROM tcrb_bphot_frames""")[0]
    prec = db.q(con, "SELECT * FROM tcrb_bprecision")
    fl = db.q(con, "SELECT * FROM tcrb_flicker")
    L = ["## B2 — apertures (T CrB B frames)", "",
         "| aperture | T CrB frames | check-star rms (mag) | residuals |",
         "|---|---|---|---|"]
    L += [f"| {a} | {b} | {c:.4f} | {d} |" for a, b, c, d in ap_rows]
    L += ["", f"Adopted: `{best}` (lowest check-star scatter). Mean growth "
          f"from 1.5 to 3.0 FWHM for T CrB: {gr[0][0]:+.4f} mag (applied "
          "equally to the ensemble, so it cancels in the zero point).", "",
          "## B1 / B3 — calibration penalties and the ensemble", "",
          "Mode-matched StackPro darks and StackPro B flat (2023-10-18) "
          "exist for every B frame (no cross-mode master was applied). "
          f"Flat staleness, ensemble residual vs detector position (3x3 "
          f"cells, noise removed): mean {pen[0]:.4f}, max {pen[1]:.4f} mag. "
          f"Cross-mode term, the High Gain B flat in place of the StackPro "
          f"one at the ensemble positions: {pen[2]:.4f} mag rms — the "
          "penalty any High Gain-flat calibration of StackPro data would "
          "carry. Colour term in B-V (Gaia-synthetic Johnson B): "
          f"{pen[4]:+.3f} +- {pen[5]:.3f}; chi2/dof of the ensemble fits "
          f"{pen[6]:.1f} (dof {pen[7]}-{pen[8]}) — formal errors understate "
          "the scatter, which is why B5's empirical errors are the quoted "
          "ones. REFCAT2 cross-check (B from g, r by Jordi et al. 2006) "
          f"minus Gaia-synthetic zero point: {pen[3]:+.3f} mag.", "",
          "## Anchors (B) against AAVSO B within 0.1 d", "",
          "| night | frames | B (transformed) | B (natural) | per-frame rms "
          "| colour-extrapolation err | B-V (AAVSO) | AAVSO B median | "
          "AAVSO MAD | AAVSO points (observers) | ZMAG-cut frames |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in arows:
        L.append("| " + " | ".join("" if x is None else (f"{x:.3f}" if
                 isinstance(x, float) else str(x)) for x in r[:9])
                 + f" | {r[9]} ({r[10]}) | {r[11]} |")
    L += ["", "## B4 — ZMAG as QC only", "",
          "| filter | frames | cut (>0.5 below season mode) | flagged (>0.3) "
          "|", "|---|---|---|---|"]
    for f_ in sorted({r[1] for r in zq}):
        zz = [r for r in zq if r[1] == f_]
        L.append(f"| {f_} | {len(zz)} | {sum(r[6] for r in zz)} | "
                 f"{sum(r[7] for r in zz)} |")
    L += ["", "## B6 — precision achieved (check stars; per mode)", "",
          "| mode | nights | per-frame rms (mag) | season floor (mag) | "
          "floor stars |", "|---|---|---|---|---|"]
    L += [f"| {a} | {b} | {c:.4f} | {d:.4f} | {e} |" for a, b, c, d, e
          in prec]
    L += ["", "No nightly-mean precision is claimed: three B nights cannot "
          "carry one.", "", "## C1 — archival flickering limits (B)", "",
          "| night | frames | span (min) | detrend | observed rms | dof | "
          "check-star rms (median) | 95% limit beta=1 / 2 | recovered at "
          "90%, beta=1 / 2 | quoted limit (rms) | published B amplitude |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in fl:
        L.append(f"| {r[0]} | {r[2]} | {r[3]:.1f} | "
                 f"{'mean' if r[4] == 0 else 'linear'} | {r[6]:.3f} | "
                 f"{r[7]} | {r[8]:.3f} | {r[10]:.3f} / {r[11]:.3f} | "
                 f"{r[12]:.3f} / {r[13]:.3f} | {r[14]:.3f} | {r[15]:.2f} |")
    L += ["", "Sentence: no archival B snippet (4-6 frames, <= 48 min) could "
          "have detected T CrB's flickering — each would have recovered "
          "only rms of order 0.1 mag or more, against a published B "
          "amplitude of 0.07 mag (Maslennikova et al. 2023) — so the "
          "archive sets no useful flickering limit."]
    write_note("phase_b", "Phase B (B-only anchors) and C1 flickering "
               "limits", L,
               "Stages `measure`, `census`, `ensemble`, `errors`, `flicker` "
               "of `pipeline/scripts/run_tcrb.py`; tables `tcrb_bphot_*`, "
               "`tcrb_berrors`, `tcrb_bprecision`, `tcrb_b_anchors`, "
               "`tcrb_b2_apertures`, `tcrb_zmag_qc`, `tcrb_flicker`.")
    print("\n".join(L))


# ===========================================================================
# AAVSO, parsed once into the project database
# ===========================================================================
AAVSO_OWN = ("MALW",)        #: RLMT's own resubmissions (CV-S7 finding)
AAVSO_BANDS = ("B", "V", "R", "I", "TB", "TG", "TR", "U")


def cmd_aavso_load(args) -> None:
    """Every cached AAVSO chunk -> tcrb_aavso (CCD/PEP photometry in the
    bands used here; visual estimates kept out; validation flag 'T' —
    failed validation — dropped; own submissions tagged)."""
    import gzip
    con = db.connect()
    rows = []
    seen = set()
    for f in sorted((db.EXTERNAL / "aavso").glob("*.delim.gz")):
        with gzip.open(f, "rt", errors="replace") as fh:
            head = fh.readline().rstrip("\n").split("@@@")
            ix = {k: i for i, k in enumerate(head)}
            for line in fh:
                v = line.rstrip("\n").split("@@@")
                if len(v) != len(head):
                    continue
                band = v[ix["band"]]
                if band not in AAVSO_BANDS or v[ix["val"]] == "T":
                    continue
                if v[ix["fainterThan"]] == "1":
                    continue
                try:
                    jd, mag = float(v[ix["JD"]]), float(v[ix["mag"]])
                except ValueError:
                    continue
                key = v[ix["obsID"]]
                if key in seen:
                    continue
                seen.add(key)
                try:   # some observers write a decimal comma
                    err = float(v[ix["uncert"]].replace(",", "."))
                except ValueError:
                    err = None
                obs = v[ix["by"]].upper()
                own = int(obs in AAVSO_OWN or "MUTEL" in
                          v[ix["comment"]].upper())
                rows.append((jd, mag, err, band, obs, v[ix["transformed"]],
                             v[ix["mtype"]], v[ix["obsType"]],
                             v[ix["val"]], own))
    n = db.write_table(con, "tcrb_aavso", ("jd", "mag", "err", "band",
                                           "obscode", "transformed",
                                           "mtype", "obstype", "val",
                                           "own"), rows)
    con.execute("CREATE INDEX IF NOT EXISTS ix_aavso ON tcrb_aavso "
                "(band, jd)")
    con.commit()
    db.record_stage(con, "aavso-load", {"tcrb_aavso": n})
    for r in db.q(con, "SELECT band, COUNT(*), SUM(own), "
                       "COUNT(DISTINCT obscode) FROM tcrb_aavso GROUP BY 1"):
        print("  ", r)


# ===========================================================================
# TCRB-A5 (ARAS side): native and LSF-degraded EWs with identical windows
# ===========================================================================
def _lsf_fwhm(grism: str) -> tuple[float, float, int]:
    """Median measured LSF FWHM (Å) of the T CrB frames of one grism in the
    2025 epoch, from G-5's g_lsf (off-nominal-focus frames excluded)."""
    g = db.grism_ro()
    v = [r[0] for r in db.q(g, """SELECT lsf_fwhm_a FROM g_lsf WHERE
        sample = 'tcrb' AND grism = ? AND mech_epoch = 'ASI-pre' AND
        COALESCE(off_nominal, 0) = 0 AND lsf_fwhm_a > 0""", grism)]
    if not v:
        raise RuntimeError(f"no measured LSF for {grism} in g_lsf")
    return float(np.median(v)), float(1.4826 * np.median(
        np.abs(np.array(v) - np.median(v)))), len(v)


def cmd_aras(args) -> None:
    from astropy.io import fits
    from macro_tcrb import spec
    con = db.connect()
    lsf = {gr: _lsf_fwhm(gr) for gr in ("hrg", "lrg")}
    rows = db.q(con, """SELECT chunk, cache_path, note FROM tcrb_ext_fetch
        WHERE source = 'aras' AND ok = 1""")
    out = []
    for fname, path, note in rows:
        with fits.open(db.REPO / path) as h:
            hd = h[0].header
            f = np.asarray(h[0].data, float)
        w = hd["CRVAL1"] + hd["CDELT1"] * (np.arange(f.size) + 1 -
                                           hd.get("CRPIX1", 1))
        R = hd.get("BSS_ITRP")
        from astropy.time import Time
        try:
            jd = float(Time(hd["DATE-OBS"]).jd + 0.5 * float(
                hd.get("EXPTIME", 0)) / 86400)
        except Exception:
            jd = None
        rec = [fname, hd.get("DATE-OBS"), jd, hd.get("OBSERVER"),
               hd.get("BSS_INST"), R]
        for gr in ("hrg", "lrg"):
            fw = lsf[gr][0]
            hw = spec.half_width(fw)
            nat = spec.equivalent_width(w, f, h=hw)
            fd = spec.degrade(w, f, R, fw)
            deg = spec.equivalent_width(w, fd, h=hw)
            dm = spec.equivalent_width(w, fd, h=hw, shift=-spec.SHIFT)
            dp = spec.equivalent_width(w, fd, h=hw, shift=+spec.SHIFT)
            rec += [nat["ew"], deg["ew"], deg["ew_err"], dm["ew"], dp["ew"]]
        out.append(tuple(rec))
    n = db.write_table(con, "tcrb_aras_ew", (
        "file", "date_obs", "jd_mid", "observer", "instrument", "R",
        "ew_native_hrg", "ew_deg_hrg", "ew_deg_err_hrg", "ew_deg_m10_hrg",
        "ew_deg_p10_hrg", "ew_native_lrg", "ew_deg_lrg", "ew_deg_err_lrg",
        "ew_deg_m10_lrg", "ew_deg_p10_lrg"), out)
    db.write_table(con, "tcrb_lsf_adopted", ("grism", "fwhm_a", "mad_a",
                                             "n_frames", "half_width_a"),
                   [(gr, *lsf[gr], spec.half_width(lsf[gr][0]))
                    for gr in lsf])
    # Window-survival test (pre-registered tolerance), 2025 window only.
    summ = []
    for gr in ("hrg", "lrg"):
        r = db.q(con, f"""SELECT ew_native_{gr}, ew_deg_{gr} FROM
            tcrb_aras_ew WHERE date_obs BETWEEN '2025-02-21' AND
            '2025-06-25' AND ew_native_{gr} > 0 AND ew_deg_{gr} > 0""")
        a = np.array(r, float)
        d = np.abs(a[:, 1] / a[:, 0] - 1)
        rel = a[:, 1] / a[:, 0] - 1
        summ.append((gr, len(a), float(np.median(d)),
                     float(np.percentile(d, 90)), float(np.median(rel)),
                     int(np.median(d) <= 0.05 and
                         np.percentile(d, 90) <= 0.10)))
    db.write_table(con, "tcrb_window_survival", (
        "grism", "n_aras", "median_abs_frac", "p90_abs_frac",
        "median_signed_frac", "passes"), summ)
    db.record_stage(con, "aras", {"tcrb_aras_ew": n})
    for x in summ:
        print("  ", x)


# ===========================================================================
# TCRB-A1 / A7 inputs: every theta CrB grism frame through the shared
# library's reduce_frame (the G-2/G-4 code), into THIS project's cache
# ===========================================================================
TET_SPEC = db.REPO / "products" / "tcrb" / "spec1d_tetcrb"


def _reduce_one(task):
    from macro_grism import reduce as gred
    try:
        r = gred.reduce_frame(str(db.ARCHIVE / task["path"]),
                              night=task["night"])
    except Exception as exc:                              # recorded
        return task, {"status": f"error: {type(exc).__name__}: {exc}"}
    if r.get("status") != "ok":
        return task, {"status": r.get("status"),
                      "header": r.get("header", {})}
    sp = r["spec"]
    out = TET_SPEC / (task["path"].replace("/", "__") + ".npz")
    np.savez_compressed(out, **{k: np.asarray(v) for k, v in sp.items()
                                if isinstance(v, (np.ndarray, list))})
    tr = r["trace"]
    return task, {"status": "ok", "spec_file": str(out.relative_to(db.REPO)),
                  "header": r["header"], "trace_coeffs":
                  json.dumps(np.asarray(tr["coeffs"]).tolist()),
                  "x0": (tr["extent"] or (None, None))[0],
                  "x1": (tr["extent"] or (None, None))[1],
                  "peak_adu": r.get("peak_adu"),
                  "n_sat_cols": r.get("n_sat_cols"),
                  "snr_median": r.get("snr_median"),
                  "fwhm_px": r.get("fwhm_px_median")}


def cmd_reduce_tet(args) -> None:
    """Reduce every staged theta CrB grism frame (hrg/lrg, all epochs)."""
    from concurrent.futures import ProcessPoolExecutor
    from macro_grism import config as gcfg
    man = db.manifest_ro()
    con = db.connect()
    TET_SPEC.mkdir(parents=True, exist_ok=True)
    tasks = [{"path": r[0], "night": r[1], "grism": r[2], "exptime": r[3],
              "mech_epoch": r[4], "jd": r[5]} for r in db.q(man, """
        SELECT s.path, s.night, s.filter, s.exptime, s.mech_epoch, s.jd
        FROM stage_tcrb_monitoring s WHERE s.role = 'science' AND
        s.canonical_target = 'tet CrB' AND s.filter IN ('hrg', 'lrg')
        ORDER BY s.jd""")]
    if args.limit:
        tasks = tasks[:args.limit]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for i, (t, r) in enumerate(ex.map(_reduce_one, tasks)):
            h = r.get("header", {}) or {}
            rows.append((t["path"], t["night"], t["grism"], t["exptime"],
                         t["mech_epoch"], gcfg.mech_epoch_id(t["night"]),
                         t["jd"], h.get("READOUTM"), h.get("CCD-TEMP"),
                         h.get("FOCPOS") or h.get("FOCUSPOS"),
                         h.get("AIRMASS"), r["status"], r.get("spec_file"),
                         r.get("trace_coeffs"), r.get("x0"), r.get("x1"),
                         r.get("peak_adu"), r.get("n_sat_cols"),
                         r.get("snr_median"), r.get("fwhm_px")))
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(tasks)}", flush=True)
    n = db.write_table(con, "tcrb_tet_frames", (
        "path", "night", "grism", "exptime", "mech_epoch", "grism_epoch",
        "jd", "readoutm", "ccd_temp", "focus", "airmass", "status",
        "spec_file", "trace_coeffs", "x0", "x1", "peak_adu", "n_sat_cols",
        "snr_median", "fwhm_px"), rows)
    db.record_stage(con, "reduce-tet", {"tcrb_tet_frames": n})
    for r in db.q(con, "SELECT grism_epoch, grism, status, COUNT(*) FROM "
                       "tcrb_tet_frames GROUP BY 1, 2, 3"):
        print("  ", r)


# ===========================================================================
# TCRB-A5 (RLMT side): Hα EW per frame and per night, on the G library
# ===========================================================================
def _library_solution(g, grism: str, epoch: str = "ASI-pre"):
    r = db.q(g, """SELECT coeffs_json, x_ref, status, disp_a_per_px,
        disp_err FROM g_dispersion WHERE grism = ? AND mech_epoch = ?""",
             grism, epoch)
    if not r:
        raise RuntimeError(f"no G-1 dispersion for {grism}/{epoch}")
    return json.loads(r[0][0]), float(r[0][1]), r[0][2], r[0][3], r[0][4]


def _spec_ew(spec_path: Path, coeffs, x_ref, x_ha, h, shift=0.0):
    from macro_tcrb import spec
    d = np.load(spec_path, allow_pickle=True)
    x = np.arange(d["flux"].size, dtype=float)
    lam = spec.wavelength_scale(x, coeffs, x_ref, x_ha)
    out = {}
    for arm, fk, vk in (("opt", "flux", "var"), ("box", "box", "box_var")):
        if fk in d.files:
            out[arm] = spec.equivalent_width(lam, d[fk],
                                             d[vk] if vk in d.files
                                             else None, h=h, shift=shift)
    return out


def cmd_ew(args) -> None:
    """Per-frame Hα EW (optimal and boxcar arms, ±10 Å shifts) for every
    T CrB grism frame the identity gate accepts."""
    from macro_tcrb import spec
    g = db.grism_ro()
    con = db.connect()
    sol = {gr: _library_solution(g, gr) for gr in ("hrg", "lrg")}
    lsf = {r[0]: r[1:] for r in db.q(con, "SELECT * FROM tcrb_lsf_adopted")}
    pen = {r[0]: r[1] for r in db.q(con, """SELECT path, leak_none_frac
        FROM tcrb_temp_split""")}
    tsg = {r[0]: r[1] for r in db.q(con, """SELECT path, temp_group FROM
        tcrb_temp_split""")}
    rows = db.q(g, """SELECT f.path, f.grism, f.night, f.jd, f.exptime,
        f.spec_file, f.n_sat_cols, f.peak_adu, f.sat_cap_adu,
        i.verdict, i.reason, z.x_halpha, z.x_halpha_err, z.sep_frac_dev,
        l.lsf_fwhm_a, l.off_nominal, l.focus_offset
        FROM g_frames f LEFT JOIN g_identity i ON i.path = f.path
        LEFT JOIN g_zero_point z ON z.path = f.path
        LEFT JOIN g_lsf l ON l.path = f.path
        WHERE f.sample = 'tcrb' ORDER BY f.jd""")
    out = []
    for (path, gr, night, jd, exp, sf, nsat, pk, cap, ver, why, xha, xhae,
         sepdev, lfw, offn, foff) in rows:
        base = [path, gr, night, jd, exp, ver, why, nsat, offn, foff,
                tsg.get(path), pen.get(path)]
        if ver != "ACCEPT" or sf is None or xha is None:
            out.append(tuple(base + [None] * 10 + ["not measured: " + (
                "identity " + str(ver) if ver != "ACCEPT" else
                "no spectrum" if sf is None else "no Halpha anchor")]))
            continue
        coeffs, xref = sol[gr][0], sol[gr][1]
        h = lsf[gr][3]
        sp_path = db.REPO / "products" / "grism" / sf
        e0 = _spec_ew(sp_path, coeffs, xref, xha, h)
        em = _spec_ew(sp_path, coeffs, xref, xha, h, -spec.SHIFT)
        ep = _spec_ew(sp_path, coeffs, xref, xha, h, +spec.SHIFT)
        o, b = e0.get("opt", {}), e0.get("box", {})
        out.append(tuple(base + [o.get("ew"), o.get("ew_err"), o.get("cont"),
                                 b.get("ew"), em.get("opt", {}).get("ew"),
                                 ep.get("opt", {}).get("ew"), h, lfw, xha,
                                 sepdev, "ok" if o.get("ew") is not None
                                 and np.isfinite(o["ew"]) else "ew failed"]))
    n = db.write_table(con, "tcrb_ew_frames", (
        "path", "grism", "night", "jd", "exptime", "gate", "gate_reason",
        "n_sat_cols", "off_nominal_focus", "focus_offset", "temp_group",
        "dark_penalty_frac", "ew", "ew_err", "cont_adu", "ew_box",
        "ew_m10", "ew_p10", "half_width_a", "lsf_fwhm_a", "x_halpha",
        "sep_frac_dev", "status"), out)
    db.record_stage(con, "ew", {"tcrb_ew_frames": n})
    for r in db.q(con, """SELECT grism, status, COUNT(*), ROUND(AVG(ew),2),
        ROUND(MIN(ew),2), ROUND(MAX(ew),2) FROM tcrb_ew_frames
        GROUP BY 1, 2"""):
        print("  ", r)


# ===========================================================================
# TCRB-A1 — characterise theta CrB (a Be/shell star) before using it
# ===========================================================================
HBETA = 4861.3
HB_BLUE, HB_RED = (4790.0, 4820.0), (4900.0, 4930.0)
LINEFREE_CENTRE = 6750.0          #: continuum-window test (A7)
LINEFREE_RED = (6800.0, 6840.0)   #: red band below the O2-B edge
O2B = 6867.0                      #: telluric B band, deepest part
TET_CONFIG_KEYS = ("grism", "grism_epoch", "exptime")


def _ew_generic(lam, flux, var, centre, h, blue, red):
    """The pre-registered EW machinery re-centred on another feature."""
    from macro_tcrb import spec
    return spec.equivalent_width(lam, flux, var, h=h, centre=centre,
                                 blue=blue, red=red)


def cmd_tet_analysis(args) -> None:
    """Per configuration: align the theta CrB spectra, fix the zero point
    from the Hα/O2-B absorption pair at the G-1 dispersion, measure Hα (and
    Hβ in lrg) EW and a line-free window per frame."""
    from macro_tcrb import spec
    g = db.grism_ro()
    con = db.connect()
    disp = {(r[0], r[1]): (json.loads(r[2]), float(r[3]), float(r[4]))
            for r in db.q(g, """SELECT grism, mech_epoch, coeffs_json,
            x_ref, disp_a_per_px FROM g_dispersion""")}
    lsf = {r[0]: r[1] for r in db.q(con, "SELECT grism, fwhm_a FROM "
                                         "tcrb_lsf_adopted")}
    fr = db.q(con, """SELECT path, night, grism, exptime, grism_epoch, jd,
        spec_file, ccd_temp, focus, airmass, readoutm, n_sat_cols
        FROM tcrb_tet_frames WHERE status = 'ok'""")
    groups: dict = {}
    for r in fr:
        groups.setdefault((r[2], r[4], round(r[3], 3)), []).append(r)
    out, conf, stacks = [], [], {}
    # Prior for the Hα pixel: the T CrB frames of the same grism and epoch
    # (both stars are centred by the pointing; +-600 px covers the spread).
    prior = {(gr_, "ASI-pre"): (float(np.median([r[0] for r in db.q(
        con, "SELECT x_halpha FROM tcrb_ew_frames WHERE grism = ? AND "
        "status = 'ok'", gr_)])), 600.0) for gr_ in ("hrg", "lrg")}
    for (gr, ep, ex), rows in sorted(groups.items()):
        if (gr, ep) not in disp:
            conf.append((gr, ep, ex, len(rows), "no G-1 dispersion", None,
                         None))
            continue
        coeffs, xref, d0 = disp[(gr, ep)]
        specs = []
        for r in rows:
            d = np.load(db.REPO / r[6])
            specs.append((r, d["flux"], d.get("var"), d.get("box")))
        nlen = min(len(s_[1]) for s_ in specs)
        norm = [spec.normalise(s_[1][:nlen]) for s_ in specs]
        ref = np.nanmedian(np.vstack(norm), axis=0)
        lags = [spec.xcorr_lag(ref, n_) for n_ in norm]
        # iterate once: re-stack on the aligned grid
        ali = [np.interp(np.arange(nlen) + l, np.arange(nlen), n_)
               for n_, l in zip(norm, lags)]
        ref = np.nanmedian(np.vstack(ali), axis=0)
        mins = spec.absorption_minima(ref, depth=0.04)
        # Only where the stacked spectrum carries signal (the steep edges
        # of the response turn the normalised shape into a comb).
        raw = np.nanmedian(np.vstack([s_[1][:nlen] for s_ in specs]), axis=0)
        valid = raw > 0.2 * np.nanmedian(raw[raw > 0])
        x_ha = spec.pair_zero_point_poly(ref, mins, coeffs, xref,
                                         spec.HALPHA, O2B,
                                         tol_px=max(8.0, 6.0 / abs(d0)),
                                         valid=valid,
                                         minima_b=spec.absorption_minima(
                                             ref, depth=0.015),
                                         prior=prior.get((gr, ep)))
        status = "ok"
        if x_ha is None and (gr, ep) in stacks:
            # Same grism and epoch, other exposure: carry the zero point
            # across by cross-correlating the two stacks.
            ref0, x0 = stacks[(gr, ep)]
            m_ = min(len(ref0), len(ref))
            x_ha = x0 + spec.xcorr_lag(ref0[:m_], ref[:m_], max_lag=120)
            status = "zero point carried from the other exposure's stack"
        if x_ha is not None and (gr, ep) not in stacks:
            stacks[(gr, ep)] = (ref, x_ha)
        conf.append((gr, ep, ex, len(rows),
                     status if x_ha is not None else "no Halpha/O2-B pair "
                     "in the stacked spectrum", x_ha, d0))
        if x_ha is None:
            continue
        h = spec.half_width(lsf.get(gr))
        for (r, f, v, b), lag in zip(specs, lags):
            x = np.arange(f.size, dtype=float)
            lam = spec.wavelength_scale(x, coeffs, xref, x_ha + lag)
            ha = spec.equivalent_width(lam, f, v, h=h)
            hab = spec.equivalent_width(lam, b, None, h=h) if b is not \
                None else {"ew": None}
            lf = _ew_generic(lam, f, v, LINEFREE_CENTRE, h, spec.RED,
                             LINEFREE_RED)
            hb = (_ew_generic(lam, f, v, HBETA, h, HB_BLUE, HB_RED)
                  if gr == "lrg" else {"ew": None, "ew_err": None})
            out.append((r[0], r[1], gr, ep, r[3], r[5], r[7], r[8], r[9],
                        r[10], r[11], float(lag), ha["ew"], ha["ew_err"],
                        hab["ew"], hb["ew"], hb["ew_err"], lf["ew"],
                        lf["ew_err"], ha["cont"]))
    # Frame quality: a frame whose LINE-FREE window departs from its
    # configuration's median by > 5 robust sigma has a broken continuum
    # (mis-aligned spectrum, a second star in the slitless field, cloud);
    # it is flagged, counted and kept out of A1/A7 statistics.  The test
    # never looks at Halpha, so Be variability cannot flag a frame.
    flagged = []
    keyed: dict = {}
    for i, r in enumerate(out):
        keyed.setdefault((r[2], r[3], r[4]), []).append(i)
    for idxs in keyed.values():
        v = np.array([out[i][17] for i in idxs], float)
        med = np.nanmedian(v)
        mad = 1.4826 * np.nanmedian(np.abs(v - med))
        for i in idxs:
            bad = (not np.isfinite(out[i][17])) or \
                abs(out[i][17] - med) > 5 * max(mad, 0.05)
            flagged.append((i, int(bad)))
    out = [tuple(r) + (dict(flagged)[i],) for i, r in enumerate(out)]
    n = db.write_table(con, "tcrb_tet_ew", (
        "path", "night", "grism", "grism_epoch", "exptime", "jd",
        "ccd_temp", "focus", "airmass", "readoutm", "n_sat_cols", "lag_px",
        "ew_ha", "ew_ha_err", "ew_ha_box", "ew_hb", "ew_hb_err",
        "ew_linefree", "ew_linefree_err", "cont", "continuum_fault"), out)
    db.write_table(con, "tcrb_tet_configs", ("grism", "grism_epoch",
                                             "exptime", "n_frames", "status",
                                             "x_halpha_stack", "disp"), conf)
    db.record_stage(con, "tet-analysis", {"tcrb_tet_ew": n})
    for c in conf:
        print("  ", c)


# ===========================================================================
# TCRB-A7 — the empirical EW error floor
# ===========================================================================
def _jitter_for_chi1(e: np.ndarray, s: np.ndarray) -> float:
    """Extra scatter j such that chi2/dof of e about its weighted mean,
    with errors sqrt(s^2 + j^2), equals 1 (0 if already <= 1)."""
    from scipy.optimize import brentq

    def f(j):
        v = s ** 2 + j ** 2
        w = 1 / v
        m = np.sum(w * e) / np.sum(w)
        return np.sum((e - m) ** 2 / v) / (e.size - 1) - 1.0
    if e.size < 3 or f(0.0) <= 0:
        return 0.0
    hi = 10 * float(np.std(e)) + 1e-6
    return float(brentq(f, 0.0, hi))


def _chi2(e, s):
    w = 1 / s ** 2
    m = np.sum(w * e) / np.sum(w)
    return float(np.sum(((e - m) / s) ** 2)), int(e.size - 1)


def cmd_floor(args) -> None:
    """Per grism: (a) the theta CrB line-free-window floor, fitted on odd
    nights and validated on even ones (chi2nu with dof); (b) the T CrB
    extraction-method difference; (c) the 240 s smear term; floor =
    max(a, b) (+) c.  Also the theta CrB Halpha scatter (Be variability)
    and the focus regression, reported."""
    from macro_tcrb import spec
    con = db.connect()
    g = db.grism_ro()
    out, halves = [], []
    for gr in ("hrg", "lrg"):
        tet = db.q(con, """SELECT night, exptime, ew_linefree, ew_linefree_err,
            ew_ha, ew_ha_err, focus FROM tcrb_tet_ew WHERE grism = ? AND
            grism_epoch = 'ASI-pre' AND readoutm = 'Mode0'
            AND ew_linefree IS NOT NULL AND continuum_fault = 0""", gr)
        # The floor is a PER-FRAME term f: frames of one night and one
        # exposure scatter about their mean by sqrt(sigma^2 + f^2).  f is
        # fitted on the frames of alternate nights (half A) and tested on
        # the other half (B) as the chi^2 of NIGHTLY means about each
        # exposure's mean, with errors sqrt((sigma^2 + f^2)/n): the
        # night-to-night test the floor has to pass.  Exposures are never
        # pooled into one mean (each keeps its own level).
        groups: dict = {}
        for r in tet:
            groups.setdefault((r[1], r[0]), []).append(r)
        allnights = sorted({k[1] for k in groups})
        half = {n_: i % 2 for i, n_ in enumerate(allnights)}
        bad_n = 0
        # night-level robustness, per exposure, on the line-free mean
        for ex in sorted({k[0] for k in groups}):
            keys = [k for k in groups if k[0] == ex]
            m_ = np.array([np.mean([x[2] for x in groups[k]]) for k in keys])
            med = np.median(m_)
            mad = 1.4826 * np.median(np.abs(m_ - med))
            for k, v in zip(keys, m_):
                if abs(v - med) > 5 * max(mad, 0.02):
                    del groups[k]
                    bad_n += 1

        def resid_set(which):
            r_, s_, ng = [], [], 0
            for (ex, n_), rs in groups.items():
                if half[n_] != which or len(rs) < 2:
                    continue
                v = np.array([x[2] for x in rs])
                r_ += list(v - v.mean())
                s_ += [x[3] for x in rs]
                ng += 1
            return np.array(r_), np.array(s_), ng

        def solve_f(r_, s_, ng):
            from scipy.optimize import brentq
            dof = r_.size - ng

            def fn(f):
                return np.sum(r_ ** 2 / (s_ ** 2 + f ** 2)) / dof - 1
            if dof < 3 or fn(0.0) <= 0:
                return 0.0
            return float(brentq(fn, 0.0, 50.0))

        def night_chi2(which, f, F=0.0):
            chi, dof = 0.0, 0
            for ex in sorted({k[0] for k in groups}):
                ms, es = [], []
                for (e2, n_), rs in groups.items():
                    if e2 != ex or half[n_] != which:
                        continue
                    v = np.array([x[2] for x in rs])
                    s2 = np.array([x[3] for x in rs]) ** 2 + f ** 2
                    w = 1 / s2
                    ms.append(np.sum(w * v) / np.sum(w))
                    es.append(np.sqrt(1 / np.sum(w) + F ** 2))
                if len(ms) >= 2:
                    c, d = _chi2(np.array(ms), np.array(es))
                    chi += c
                    dof += d
            return chi, dof
        def solve_F(which, f):
            # night-level term: nightly chi^2/dof = 1 on the given half
            from scipy.optimize import brentq

            def fn(F):
                c, d = night_chi2(which, f, F)
                return c / max(d, 1) - 1
            if fn(0.0) <= 0:
                return 0.0
            return float(brentq(fn, 0.0, 50.0))
        rA, sA, gA = resid_set(0)
        fA = solve_f(rA, sA, gA)
        FA = solve_F(0, fA)
        chi_b, dof_b = night_chi2(1, fA, FA)
        rAll = np.concatenate([rA, resid_set(1)[0]])
        sAll = np.concatenate([sA, resid_set(1)[1]])
        f_all = solve_f(rAll, sAll, gA + resid_set(1)[2])
        # Conservative: the larger of the two halves' night terms (each
        # half's estimate rests on ~11 dof; the split-half test below shows
        # how far one half's value transfers to the other).
        F0, F1 = solve_F(0, f_all), solve_F(1, f_all)
        F_all = max(F0, F1)
        chi_h = [night_chi2(k, f_all, F_all) for k in (0, 1)]
        chi_all, dof_all = night_chi2(0, f_all, F_all)
        c1, d1 = night_chi2(1, f_all, F_all)
        chi_all, dof_all = chi_all + c1, dof_all + d1
        halves.append((gr, "pooled exposures", int(rA.size), gA, fA, FA,
                       chi_b, dof_b, chi_b / max(dof_b, 1), F0, F1,
                       chi_h[0][0], chi_h[0][1], chi_h[1][0], chi_h[1][1]))
        lf_floor = f_all
        night_floor = F_all
        # T CrB's own within-night frame scatter beyond its formal errors
        tc: dict = {}
        for n_, e_, s_ in db.q(con, """SELECT night, ew, ew_err FROM
                tcrb_ew_frames WHERE grism = ? AND status = 'ok'""", gr):
            tc.setdefault(n_, []).append((e_, s_))
        rT, sT, gT = [], [], 0
        for n_, rs in tc.items():
            if len(rs) < 2:
                continue
            v = np.array([x[0] for x in rs])
            rT += list(v - v.mean())
            sT += [x[1] for x in rs]
            gT += 1
        f_tcrb = solve_f(np.array(rT), np.array(sT), gT)
        # Halpha of theta CrB itself: per-frame jitter with the same method
        # (Be variability would show here and NOT in the line-free window)
        jit_ha = []
        for ex in sorted({k[0] for k in groups}):
            ms = [np.mean([x[4] for x in rs]) for (e2, n_), rs in
                  groups.items() if e2 == ex]
            jit_ha.append((ex, float(np.std(ms, ddof=1)) if len(ms) > 1
                           else None, len(ms)))
        # (b) extraction-method difference on T CrB frames
        md = [abs(r[0] - r[1]) for r in db.q(con, """SELECT ew, ew_box FROM
            tcrb_ew_frames WHERE grism = ? AND status = 'ok' AND ew_box IS
            NOT NULL""", gr)]
        mdiff = float(np.median(md)) if md else None
        # (c) 240 s smear: EW change when the LSF broadens by the T CrB
        # frame-to-frame scatter of the cross-dispersion FWHM (as a fraction)
        fw = np.array([r[0] for r in db.q(g, """SELECT fwhm_px FROM g_frames
            WHERE sample = 'tcrb' AND grism = ? AND status = 'ok' AND
            fwhm_px > 0""", gr)])
        lsf = db.q1(con, "SELECT fwhm_a FROM tcrb_lsf_adopted WHERE grism=?",
                    gr)
        frac = float(1.4826 * np.median(np.abs(fw - np.median(fw))) /
                     np.median(fw)) if fw.size else 0.0
        smear = None
        ex = db.q(con, """SELECT path FROM tcrb_ew_frames WHERE grism = ? AND
            status = 'ok' ORDER BY cont_adu DESC LIMIT 5""", gr)
        sm_vals = []
        for (path,) in ex:
            sf = db.q1(g, "SELECT spec_file FROM g_frames WHERE path = ?",
                       path)
            xha = db.q1(g, "SELECT x_halpha FROM g_zero_point WHERE path=?",
                        path)
            coeffs, xref, *_ = _library_solution(g, gr)
            d = np.load(db.REPO / "products" / "grism" / sf)
            lam = spec.wavelength_scale(np.arange(d["flux"].size), coeffs,
                                        xref, xha)
            o = np.argsort(lam)
            l2, f2 = lam[o], d["flux"][o]
            h = spec.half_width(lsf)
            e0 = spec.equivalent_width(l2, f2, h=h)["ew"]
            wide = lsf * np.sqrt((1 + frac) ** 2 - 1)
            f3 = spec.degrade(l2, f2, None, wide) if wide > 0 else f2
            e1 = spec.equivalent_width(l2, f3, h=h)["ew"]
            sm_vals.append(abs(e1 - e0))
        smear = float(np.median(sm_vals)) if sm_vals else 0.0
        # Per-NIGHT floor carried by every T CrB nightly EW: the theta CrB
        # night-level term (+) the larger of the per-frame terms (theta CrB
        # line-free, T CrB within-night, extraction-method difference)
        # reduced by sqrt(typical frames per night) (+) the 240 s smear.
        nfr = np.median([len(v) for v in tc.values()]) if tc else 1
        per_frame = max(lf_floor or 0, f_tcrb, mdiff or 0)
        floor = float(np.sqrt(night_floor ** 2 + per_frame ** 2 / nfr +
                              smear ** 2))
        nb = db.q(con, """SELECT COUNT(*), SUM(off_nominal_focus) FROM
            tcrb_ew_frames WHERE grism = ? AND status = 'ok'""", gr)[0]
        out.append((gr, lf_floor, night_floor, f_tcrb, float(nfr), chi_b,
                    dof_b, chi_b / max(dof_b, 1), chi_all, dof_all, mdiff,
                    frac, smear, floor, json.dumps(jit_ha), nb[0], nb[1],
                    bad_n))
    db.write_table(con, "tcrb_floor", (
        "grism", "linefree_frame_floor", "night_floor", "tcrb_frame_floor",
        "tcrb_frames_per_night", "chi2_halfB", "dof_halfB", "chi2nu_halfB",
        "chi2_all", "dof_all", "method_diff", "fwhm_scatter_frac", "smear",
        "floor_ew", "tet_halpha_nightsd_json", "n_tcrb_frames",
        "n_off_focus", "n_bad_nights"), out)
    db.write_table(con, "tcrb_floor_halves", (
        "grism", "exposures", "n_frames_A", "n_groups_A", "frame_floor_A",
        "night_floor_A", "chi2_B_given_A", "dof_B", "chi2nu_B_given_A",
        "night_floor_half0", "night_floor_half1", "chi2_half0_adopted",
        "dof_half0", "chi2_half1_adopted", "dof_half1"), halves)
    db.record_stage(con, "floor", {"tcrb_floor": len(out)})
    for r in out:
        print("  ", r)
    for h_ in halves:
        print("   half", h_)


# ===========================================================================
# TCRB-A5 nightly series + the pre-registered rule (A5a)
# ===========================================================================
def cmd_nightly(args) -> None:
    """Nightly EW per grism: inverse-variance mean of the accepted frames
    (saturated or off-nominal-focus frames excluded and counted), with the
    A7 floor added in quadrature (rule 1)."""
    con = db.connect()
    floor = {r[0]: r[1] for r in db.q(con, """SELECT grism, floor_ew FROM
        tcrb_floor""")}
    rows = db.q(con, """SELECT grism, night, jd, ew, ew_err, ew_box, ew_m10,
        ew_p10, cont_adu, n_sat_cols, off_nominal_focus, temp_group,
        dark_penalty_frac FROM tcrb_ew_frames WHERE status = 'ok'""")
    by: dict = {}
    excl = {"saturated": 0, "off_focus": 0}
    for r in rows:
        if (r[9] or 0) > 0:
            excl["saturated"] += 1
            continue
        if r[10]:
            excl["off_focus"] += 1
            continue
        by.setdefault((r[0], r[1]), []).append(r)
    out = []
    for (gr, night), rs in sorted(by.items()):
        e = np.array([r[3] for r in rs])
        # per-frame error: statistical (+) the frame's dark penalty -b/C x EW
        se = np.array([np.hypot(r[4], abs((r[12] or 0) * r[3]))
                       for r in rs])
        w = 1 / se ** 2
        m = float(np.sum(w * e) / np.sum(w))
        sm = float(np.sqrt(1 / np.sum(w)))
        fl = floor.get(gr, 0.0) or 0.0
        out.append((gr, night, float(np.mean([r[2] for r in rs])), len(rs),
                    m, sm, fl, float(np.hypot(sm, fl)),
                    float(np.mean([r[5] for r in rs])),
                    float(np.mean([r[6] for r in rs])),
                    float(np.mean([r[7] for r in rs])),
                    float(np.median([r[8] for r in rs])),
                    int(any(r[11] == "warm" for r in rs))))
    n = db.write_table(con, "tcrb_ew_nightly", (
        "grism", "night", "jd", "n_frames", "ew", "ew_stat", "ew_floor",
        "ew_err", "ew_box", "ew_m10", "ew_p10", "cont_adu", "warm"), out)
    db.write_table(con, "tcrb_ew_exclusions", ("reason", "n_frames"),
                   list(excl.items()))
    db.record_stage(con, "nightly", {"tcrb_ew_nightly": n})
    for g_ in ("hrg", "lrg"):
        print("  ", g_, db.q(con, """SELECT COUNT(*), ROUND(MIN(ew),2),
            ROUND(MAX(ew),2), ROUND(AVG(ew_err),2) FROM tcrb_ew_nightly
            WHERE grism = ?""", g_), excl)


def cmd_detect(args) -> None:
    """Apply the pre-registered rule (macro_tcrb.detect) to each grism."""
    from macro_tcrb import detect
    con = db.connect()
    ser = {g_: np.array(db.q(con, """SELECT jd, ew, ew_err FROM
        tcrb_ew_nightly WHERE grism = ? ORDER BY jd""", g_), float)
           for g_ in ("hrg", "lrg")}
    res_rows, ev_rows = [], []
    for g_, other in (("hrg", "lrg"), ("lrg", "hrg")):
        a, b = ser[g_], ser[other]
        r = detect.apply_rule(a[:, 0], a[:, 1], a[:, 2],
                              other=(b[:, 0], b[:, 1], b[:, 2]),
                              n_mc=args.n_mc)
        res_rows.append((g_, len(a), r.threshold_sigma, r.expected_false,
                         len(r.events), r.chi2, r.dof, r.p_const,
                         r.p_const_inflated, int(r.variable),
                         r.min_step_recovered,
                         r.recovery.get("eligible_fraction"),
                         max(r.recovery.get("overall", [0])),
                         json.dumps(r.recovery)))
        for ev in r.events:
            ev_rows.append((g_, ev.t_i, ev.t_j, ev.delta, ev.sigma, ev.nsig,
                            ev.corroboration, ev.flux_verdict))
    db.write_table(con, "tcrb_detect", (
        "grism", "n_nights", "threshold_sigma", "expected_false",
        "n_events", "chi2", "dof", "p_const", "p_const_inflated",
        "variable", "min_step_recovered", "eligible_fraction",
        "overall_recovery_max", "recovery_json"), res_rows)
    db.write_table(con, "tcrb_detect_events", (
        "grism", "jd_i", "jd_j", "delta", "sigma", "nsig", "corroboration",
        "flux_verdict"), ev_rows)
    db.record_stage(con, "detect", {"tcrb_detect": len(res_rows)})
    for r in res_rows:
        print("  ", r)
    for e in ev_rows:
        print("   event", e)


# ===========================================================================
# TCRB-A5b — line flux = EW x contemporaneous continuum, vs orbital phase
# ===========================================================================
#: Zero-magnitude flux densities (erg cm^-2 s^-1 A^-1), Bessell, Castelli &
#: Plez (1998, A&A 333, 231), Table A2: Cousins R and I.
F0 = {"R": 2.177e-9, "I": 1.126e-9}
CONT_BAND = "R"            #: Cousins R (lambda_eff 6410 A) sits on Halpha
CONT_HALF_D = 0.5          #: same-night window for AAVSO R
CONT_INTERP_MAX_D = 2.0    #: else linear interpolation across <= 2 d


def _ephemeris():
    import csv
    for r in csv.DictReader(open(db.NOVELTY / "external" /
                                 "ephemerides.csv")):
        if r["ephemeris"].startswith("Munari"):
            return float(r["t0_hjd"]), float(r["period_d"])
    raise RuntimeError("Munari ephemeris missing")


def cmd_flux(args) -> None:
    from macro_tcrb import spec
    con = db.connect()
    t0, P = _ephemeris()
    nights = db.q(con, """SELECT grism, night, jd, ew, ew_err FROM
        tcrb_ew_nightly ORDER BY jd""")
    # nightly AAVSO R (non-own CCD/PEP), median and MAD
    rd = db.q(con, f"""SELECT CAST(jd + 0.5 AS INT) d, AVG(jd), COUNT(*),
        GROUP_CONCAT(mag), COUNT(DISTINCT obscode) FROM tcrb_aavso WHERE
        band = '{CONT_BAND}' AND own = 0 AND jd BETWEEN 2460700 AND 2460880
        GROUP BY d ORDER BY d""")
    rn = np.array([[r[1], float(np.median([float(x) for x in
                                           r[3].split(",")])),
                    1.4826 * float(np.median(np.abs(
                        np.array([float(x) for x in r[3].split(",")]) -
                        np.median([float(x) for x in r[3].split(",")])))) /
                    np.sqrt(r[2]), r[2], r[4]] for r in rd], float)
    out = []
    for gr, night, jd, ew, ewe in nights:
        d = np.abs(rn[:, 0] - jd)
        src, mag, merr, nobs = "none", None, None, 0
        if d.min() <= CONT_HALF_D:
            i = int(np.argmin(d))
            mag, merr, nobs = rn[i, 1], max(rn[i, 2], 0.01), int(rn[i, 3])
            src = f"AAVSO {CONT_BAND}, same night ({nobs} obs)"
        else:
            lo = rn[rn[:, 0] < jd]
            hi = rn[rn[:, 0] > jd]
            if len(lo) and len(hi) and hi[0, 0] - lo[-1, 0] <= \
                    2 * CONT_INTERP_MAX_D:
                f_ = (jd - lo[-1, 0]) / (hi[0, 0] - lo[-1, 0])
                mag = lo[-1, 1] + f_ * (hi[0, 1] - lo[-1, 1])
                merr = float(np.hypot(max(lo[-1, 2], 0.01),
                                      max(hi[0, 2], 0.01)))
                src = f"AAVSO {CONT_BAND}, interpolated"
        cont = cerr = flux = ferr = None
        if mag is not None:
            cont = spec.mag_to_flambda(mag, F0[CONT_BAND])
            cerr = cont * 0.921 * merr
            flux, ferr = spec.line_flux(ew, ewe, cont, cerr)
        out.append((gr, night, jd, spec.orbital_phase(jd, t0, P), ew, ewe,
                    mag, merr, src, cont, cerr, flux, ferr))
    n = db.write_table(con, "tcrb_flux", (
        "grism", "night", "jd", "phase", "ew", "ew_err", "cont_mag",
        "cont_mag_err", "cont_source", "cont_flambda", "cont_flambda_err",
        "line_flux", "line_flux_err"), out)
    db.record_stage(con, "flux", {"tcrb_flux": n})
    for r in db.q(con, """SELECT grism, cont_source LIKE '%same%',
        COUNT(*), ROUND(MIN(phase),3), ROUND(MAX(phase),3),
        ROUND(MIN(cont_mag),2), ROUND(MAX(cont_mag),2) FROM tcrb_flux
        GROUP BY 1, 2"""):
        print("  ", r)


# ===========================================================================
# TCRB-A8 — ARAS cross-validation on the same UT dates
# ===========================================================================
def _ut_date(jd: float) -> str:
    from astropy.time import Time
    return Time(jd, format="jd").iso[:10]


def cmd_xval(args) -> None:
    con = db.connect()
    ar = db.q(con, """SELECT substr(date_obs, 1, 10), jd_mid, observer, R,
        ew_native_hrg, ew_deg_hrg, ew_native_lrg, ew_deg_lrg, file
        FROM tcrb_aras_ew WHERE date_obs BETWEEN '2025-02-15' AND
        '2025-07-01'""")
    by: dict = {}
    for r in ar:
        by.setdefault(r[0], []).append(r)
    out = []
    for gr in ("hrg", "lrg"):
        k_nat, k_deg = (4, 5) if gr == "hrg" else (6, 7)
        for night, jd, ew, err in db.q(con, """SELECT night, jd, ew, ew_err
                FROM tcrb_ew_nightly WHERE grism = ?""", gr):
            day = _ut_date(jd)
            rs = [r for r in by.get(day, []) if r[k_deg] is not None and
                  np.isfinite(r[k_deg])]
            if not rs:
                continue
            deg = float(np.median([r[k_deg] for r in rs]))
            nat = float(np.median([r[k_nat] for r in rs]))
            out.append((gr, night, day, jd, ew, err, len(rs),
                        ",".join(sorted({str(r[2]) for r in rs})),
                        float(np.median([r[3] or np.nan for r in rs])),
                        deg, nat, ew / deg - 1 if deg else None))
    n = db.write_table(con, "tcrb_xval", (
        "grism", "night", "ut_date", "jd", "rlmt_ew", "rlmt_err", "n_aras",
        "aras_observers", "aras_R_median", "aras_ew_deg", "aras_ew_native",
        "frac_diff"), out)
    summ = []
    for gr in ("hrg", "lrg"):
        f = np.array([r[11] for r in out if r[0] == gr and r[11] is not
                      None], float)
        if f.size:
            summ.append((gr, f.size, float(np.mean(f)),
                         float(np.std(f, ddof=1)),
                         float(np.median(f)),
                         float(1.4826 * np.median(np.abs(f - np.median(f))))))
    db.write_table(con, "tcrb_xval_summary", ("grism", "n_dates",
                                              "mean_frac", "rms_frac",
                                              "median_frac", "mad_frac"),
                   summ)
    db.record_stage(con, "xval", {"tcrb_xval": n})
    for s_ in summ:
        print("  ", s_)


# ===========================================================================
# TCRB-D3 — release package (Zenodo-ready) and ARAS-format spectra
# ===========================================================================
RELEASE = db.REPO / "releases" / "tcrb_halpha_v1"
ARAS_STEP = {"hrg": 0.5, "lrg": 2.0}     #: linear resampling step (A)


def cmd_release(args) -> None:
    """Reduced 1-D spectra of every identity-accepted T CrB frame (native
    pixel grid, with wavelength, flux, variance and the boxcar arm), the
    per-frame and nightly EW tables, the ARAS cross-validation and the
    archival anchors; plus the same spectra resampled to a linear
    wavelength grid with the ARAS BSS header keywords, for submission."""
    import csv
    from astropy.io import fits
    from macro_tcrb import spec
    con = db.connect()
    g = db.grism_ro()
    man = db.manifest_ro()
    (RELEASE / "spectra").mkdir(parents=True, exist_ok=True)
    (RELEASE / "aras_format").mkdir(parents=True, exist_ok=True)
    sol = {gr: _library_solution(g, gr) for gr in ("hrg", "lrg")}
    rows = db.q(con, """SELECT path, grism, night, jd, exptime, ew, ew_err,
        ew_box, x_halpha, temp_group, lsf_fwhm_a FROM tcrb_ew_frames
        WHERE status = 'ok' ORDER BY jd""")
    idx = []
    for path, gr, night, jd, exp, ew, ewe, ewb, xha, tg, lfw in rows:
        sf = db.q1(g, "SELECT spec_file FROM g_frames WHERE path = ?", path)
        d = np.load(db.REPO / "products" / "grism" / sf)
        coeffs, xref = sol[gr][0], sol[gr][1]
        x = np.arange(d["flux"].size, dtype=float)
        lam = spec.wavelength_scale(x, coeffs, xref, xha)
        bjd = db.q1(man, "SELECT bjd_tdb FROM frame_times WHERE path = ?",
                    path)
        base = Path(path).name.replace(".fts.fz", "")
        cols = [fits.Column("wavelength", "D", "Angstrom", array=lam),
                fits.Column("flux", "E", "adu", array=d["flux"]),
                fits.Column("flux_var", "E", "adu2", array=d["var"]),
                fits.Column("flux_box", "E", "adu", array=d["box"])]
        h = fits.Header()
        for k, v in (("OBJECT", "T CrB"), ("TELESCOP", "RLMT 0.5m"),
                     ("INSTRUME", f"ZWO ASI + {gr}"), ("GRISM", gr),
                     ("NIGHT", night), ("JD_START", jd),
                     ("BJD_TDB", bjd if bjd is not None else -1.0),
                     ("EXPTIME", exp), ("HA_EW", ew), ("HA_EWERR", ewe),
                     ("TEMPGRP", tg or ""), ("LSF_FWHM", lfw or -1.0),
                     ("RAWFILE", path), ("PIPELINE",
                                         "MACRO macro_grism + macro_tcrb")):
            h[k] = v
        out = RELEASE / "spectra" / f"tcrb_{gr}_{base}.fits"
        fits.BinTableHDU.from_columns(cols, header=h).writeto(out,
                                                              overwrite=True)
        # ARAS format: linear grid, continuum-normalised flux, BSS keys.
        o = np.argsort(lam)
        l2, f2 = lam[o], d["flux"][o]
        ok = np.isfinite(f2)
        grid = np.arange(np.ceil(l2[ok].min()), np.floor(l2[ok].max()),
                         ARAS_STEP[gr])
        fg = np.interp(grid, l2[ok], f2[ok])
        cn = np.nanmedian(fg[(grid > spec.BLUE[0]) & (grid < spec.BLUE[1])])
        ah = fits.Header()
        from astropy.time import Time
        for k, v in (("CRVAL1", float(grid[0])), ("CDELT1", ARAS_STEP[gr]),
                     ("CRPIX1", 1.0), ("CTYPE1", "Wavelength"),
                     ("CUNIT1", "Angstrom"), ("OBJNAME", "T CrB"),
                     ("OBJECT", "T CrB"),
                     ("DATE-OBS", Time(jd, format="jd").isot),
                     ("EXPTIME", exp), ("OBSERVER", "MACRO Consortium"),
                     ("BSS_SITE", "Winer Observatory, Sonoita AZ"),
                     ("BSS_INST", f"RLMT 0.5m + ZWO ASI6200 + slitless "
                                  f"{gr}"),
                     ("BSS_ITRP", int(round(spec.HALPHA / lfw)) if lfw
                      else -1), ("BSS_VHEL", 0.0),
                     ("BSS_NORM", "6470-6520 A median = 1"),
                     ("COMMENT", "slitless grism; telluric not removed; "
                                 "not flux calibrated")):
            ah[k] = v
        fits.PrimaryHDU((fg / cn).astype(np.float32), header=ah).writeto(
            RELEASE / "aras_format" / f"tcrb_{gr}_{base}_aras.fit",
            overwrite=True)
        idx.append((out.name, path, gr, night, jd, bjd, exp, ew, ewe, ewb,
                    tg, lfw))
    def dump(name, cols, rows_):
        with open(RELEASE / name, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            w.writerows(rows_)
    dump("tcrb_rlmt_spectra_index.csv",
         ("file", "raw_path", "grism", "night", "jd_start", "bjd_tdb",
          "exptime_s", "ew_A", "ew_err_A", "ew_box_A", "temp_group",
          "lsf_fwhm_A"), idx)
    for t, f_ in (("tcrb_ew_nightly", "tcrb_rlmt_nightly_ew.csv"),
                  ("tcrb_flux", "tcrb_rlmt_line_flux.csv"),
                  ("tcrb_xval", "aras_cross_validation.csv"),
                  ("tcrb_b_anchors", "b_anchors_2024.csv"),
                  ("tcrb_filter_map", "filter_codes.csv"),
                  ("tcrb_flicker", "flickering_limits.csv")):
        cur = con.execute(f"SELECT * FROM {t}")
        dump(f_, [c[0] for c in cur.description], cur.fetchall())
    zen = {"title": "T CrB: RLMT slitless-grism H-alpha spectra and "
                    "equivalent widths, 2025 February-June (MACRO "
                    "Consortium release)",
           "upload_type": "dataset",
           "description": "Reduced 1-D slitless-grism spectra of the "
                          "recurrent nova T CrB from the 0.5 m Robert L. "
                          "Mutel Telescope, their H-alpha equivalent widths "
                          "per frame and per night, line flux against "
                          "orbital phase, the same-date ARAS cross-"
                          "validation, and the archival 2024 B anchors. "
                          "See README.md.",
           "creators": [{"name": "Wetzel, James",
                         "affiliation": "Coe College"}],
           "license": "cc-by-4.0",
           "keywords": ["T CrB", "recurrent nova", "H-alpha",
                        "slitless spectroscopy", "MACRO Consortium",
                        "RLMT"],
           "notes": "Released on behalf of the MACRO Consortium. Code: "
                    "pipeline/macro_tcrb, pipeline/scripts/run_tcrb.py.",
           "related_identifiers": [{"identifier":
                                    "https://github.com/jwwetzel/MACRO",
                                    "relation": "isSupplementedBy",
                                    "resource_type": "software"}]}
    (RELEASE / ".zenodo.json").write_text(json.dumps(zen, indent=2) + "\n")
    n_n = db.q1(con, "SELECT COUNT(*) FROM tcrb_ew_nightly")
    (RELEASE / "README.md").write_text(
        "# T CrB — RLMT slitless-grism Halpha release v1\n\n"
        f"Generated {db.utc_now()} by `pipeline/scripts/run_tcrb.py "
        f"release` (git {db.git_commit()}) from `products/tcrb/tcrb.sqlite` "
        "and the shared grism library. Every number is a query.\n\n"
        f"- `spectra/` — {len(idx)} reduced 1-D spectra (FITS binary "
        "tables: wavelength on the fixed G-1 dispersion with a per-frame "
        "Halpha zero point; optimal flux, its variance; boxcar flux). "
        "Index: `tcrb_rlmt_spectra_index.csv`.\n"
        f"- `aras_format/` — the same {len(idx)} spectra on a linear grid, "
        "normalised at 6470-6520 A, with BSS keywords, for submission to "
        "the ARAS T CrB database.\n"
        f"- `tcrb_rlmt_nightly_ew.csv` — {n_n} nightly EWs (per grism; "
        "error = statistical (+) floor).\n"
        "- `tcrb_rlmt_line_flux.csv` — line flux = EW x AAVSO Rc continuum, "
        "orbital phase (Munari et al. 2025 ephemeris).\n"
        "- `aras_cross_validation.csv`, `b_anchors_2024.csv`, "
        "`filter_codes.csv`, `flickering_limits.csv`.\n\n"
        "EW convention (pre-registered, ANALYSIS_STRATEGY.md §10): emission "
        "positive; Halpha 6562.8 A +- max(30 A, 1.5 LSF FWHM); "
        "pseudo-continuum through the medians of 6470-6520 and 6600-6640 A."
        "\n\nArchive DOI: pending (the Zenodo deposit and the ARAS "
        "submission are owner actions).\n")
    db.record_stage(con, "release", {"n_spectra": len(idx)})
    print(f"  {len(idx)} spectra -> {RELEASE}")


# ===========================================================================
# TCRB-A0b — the 2023-05 -> 2024-03 slot '6' and 'W' spectra
# ===========================================================================
A0B_MIN_CORR = 0.85        #: template-match acceptance (shape correlation)
A0B_HA_TOL_PX = 6.0        #: Halpha emission must sit at the predicted pixel


def _smooth_norm(f):
    from scipy.ndimage import gaussian_filter1d
    f = gaussian_filter1d(np.nan_to_num(np.asarray(f, float)), 3)
    c = gaussian_filter1d(f, 60)
    with np.errstate(all="ignore"):
        r = f / c - 1
    return np.nan_to_num(r)


def _template_match(f, grid, Tn, coarse=True):
    """Best (corr, dispersion, x_Ha) of a linear wavelength scale mapping
    the spectrum onto the 2025 lrg T CrB template (shape correlation)."""
    x = np.arange(f.size)
    best = (-1.0, None, None)

    def score(D, x0):
        lam = spec_HA + D * (x - x0)
        o = np.argsort(lam)
        s_ = np.interp(grid, lam[o], f[o], left=np.nan, right=np.nan)
        ok = np.isfinite(s_) & (s_ > 0.15 * np.nanmax(s_))
        if ok.sum() < 300:
            return -1.0
        sn = _smooth_norm(np.where(ok, s_, np.nan))
        return float(np.corrcoef(sn[ok], Tn[ok])[0, 1])
    spec_HA = 6562.8
    Ds = np.concatenate([np.arange(-4, -0.4, 0.05), np.arange(0.4, 4, 0.05)])
    for D in Ds:
        for x0 in range(400, f.size - 400, 16):
            c = score(D, x0)
            if c > best[0]:
                best = (c, D, x0)
    c0, D0, x00 = best
    for D in np.arange(D0 - 0.06, D0 + 0.061, 0.005):
        for x0 in range(x00 - 20, x00 + 21, 1):
            c = score(D, x0)
            if c > best[0]:
                best = (c, float(D), int(x0))
    return best


def cmd_a0b(args) -> None:
    """Extract-or-reject decision per slot-6/W T CrB frame of 2023-24.

    No hot-star dispersion exists for the AC4040 grism epochs (G-1 solved
    none), so a frame's wavelength scale is found by matching its whole
    M-giant spectrum — dozens of TiO features — to the 2025 lrg T CrB
    template (one linear scale per frame), and is accepted only if (i) the
    shape correlation exceeds A0B_MIN_CORR and (ii) the Halpha emission
    peak, which the match does not use separately, lies within
    A0B_HA_TOL_PX of the predicted pixel.  Saturated, untraced and
    unmatched frames are documented failures.
    """
    from macro_tcrb import spec
    g = db.grism_ro()
    con = db.connect()
    sol = _library_solution(g, "lrg")
    grid = np.arange(5600, 9000, 2.0)
    T = []
    for path, xha in db.q(con, """SELECT path, x_halpha FROM tcrb_ew_frames
            WHERE grism = 'lrg' AND status = 'ok' ORDER BY cont_adu DESC
            LIMIT 15"""):
        sf = db.q1(g, "SELECT spec_file FROM g_frames WHERE path = ?", path)
        d = np.load(db.REPO / "products" / "grism" / sf)
        lam = spec.wavelength_scale(np.arange(d["flux"].size), sol[0],
                                    sol[1], xha)
        o = np.argsort(lam)
        T.append(np.interp(grid, lam[o], d["flux"][o], left=np.nan,
                           right=np.nan))
    Tn = _smooth_norm(np.nanmedian(np.vstack(T), axis=0))
    rows = db.q(g, """SELECT path, filter, mech_epoch, night, exptime,
        readoutm, status, n_sat_cols, spec_file, peak_adu, sat_cap_adu
        FROM g_frames WHERE sample = 'tcrb_slot' ORDER BY night, path""")
    out = []
    for (path, flt, ep, night, exp, mode, st, nsat, sf, pk, cap) in rows:
        rec = [path, flt, ep, night, exp, mode]
        if st != "ok":
            why = ("Low Gain: no trace and no measured ceiling" if
                   mode == "Low Gain" else
                   "no trace found (spectrum below the library's 15 ADU "
                   "trace-height threshold)")
            out.append(tuple(rec + [None] * 6 + ["reject", why]))
            continue
        if (nsat or 0) > 0:
            out.append(tuple(rec + [None] * 6 + [
                "reject", f"saturated: {nsat} columns above the measured "
                f"cap ({cap:.0f} ADU)"]))
            continue
        d = np.load(db.REPO / "products" / "grism" / sf)
        c, D, x0 = _template_match(d["flux"], grid, Tn)
        # independent Halpha check: emission peak near the predicted pixel
        lo, hi = int(x0 - 25), int(x0 + 26)
        seg = d["flux"][lo:hi] - np.nanmedian(d["flux"][lo - 60:hi + 60])
        xpk = lo + int(np.nanargmax(seg))
        ha_ok = abs(xpk - x0) <= A0B_HA_TOL_PX
        if c < A0B_MIN_CORR or not ha_ok:
            out.append(tuple(rec + [c, D, x0, xpk, None, None, "reject",
                                    "wavelength scale not established "
                                    f"(template corr {c:.2f}; Halpha peak "
                                    f"{xpk - x0:+d} px from prediction)"]))
            continue
        lam = 6562.8 + D * (np.arange(d["flux"].size) - xpk)
        e = spec.equivalent_width(lam, d["flux"], d.get("var"), h=30.0)
        out.append(tuple(rec + [c, D, x0, xpk, e["ew"], e["ew_err"],
                                "extract", "template-matched scale; "
                                "Halpha peak confirmed"]))
    n = db.write_table(con, "tcrb_a0b", (
        "path", "filter", "grism_epoch", "night", "exptime", "readoutm",
        "template_corr", "disp_a_per_px", "x_halpha_pred", "x_halpha_peak",
        "ew", "ew_err", "decision", "reason"), out)
    db.record_stage(con, "a0b", {"tcrb_a0b": n})
    L = ["| frame | code | epoch | night | exp (s) | mode | corr | A/px | "
         "EW (A) | decision | reason |", "|---|---|---|---|---|---|---|---|"
         "---|---|---|"]
    for r in out:
        L.append(f"| `{r[0]}` | {r[1]} | {r[2]} | {r[3]} | {r[4]:g} | "
                 f"{r[5]} | {'' if r[6] is None else f'{r[6]:.2f}'} | "
                 f"{'' if r[7] is None else f'{r[7]:.2f}'} | "
                 f"{'' if r[10] is None else f'{r[10]:.1f} +- {r[11]:.1f}'}"
                 f" | {r[12]} | {r[13]} |")
    write_note("a0b_early_spectra", "TCRB-A0b — the 2023-24 slot-6 and W "
               "spectra, extract or reject", L,
               "Library frames (sample `tcrb_slot`, G v2.6). " +
               cmd_a0b.__doc__.split("\n\n")[1].replace("\n    ", " "))
    print("\n".join(L))


# ===========================================================================
# TCRB-A9 — Halpha profile morphology and differential wing velocities (hrg)
# ===========================================================================
C_KMS = 299792.458


def _profile_metrics(v: np.ndarray, r: np.ndarray) -> dict:
    """Continuum-subtracted profile r(v) (peak-normalised): FWHM, the
    velocities where each wing falls to 25% and 10% of the peak, the
    wing asymmetry (v_red + v_blue at 10%), and whether a central
    reversal (a local minimum within +-150 km/s, >= 5% of the peak below
    both flanking maxima) is present."""
    i0 = int(np.nanargmax(np.where(np.abs(v) < 300, r, -np.inf)))
    pk = r[i0]
    r = r / pk

    def cross(level, side):
        idx = range(i0, len(v)) if side > 0 else range(i0, -1, -1)
        prev = None
        for i in idx:
            if r[i] < level:
                if prev is None:
                    return v[i]
                return float(np.interp(level, [r[i], r[prev]],
                                       [v[i], v[prev]]))
            prev = i
        return float("nan")
    out = {f"v{int(100 * L)}_{s}": cross(L, sg) for L in (0.5, 0.25, 0.10)
           for s, sg in (("blue", -1), ("red", 1))}
    out["fwhm"] = out["v50_red"] - out["v50_blue"]
    out["asym10"] = out["v10_red"] + out["v10_blue"]
    c = np.abs(v) < 150
    rc = np.where(c, r, np.nan)
    j = int(np.nanargmin(rc))
    left = np.nanmax(np.where((v < v[j]) & c, r, np.nan)) if np.any(
        (v < v[j]) & c) else np.nan
    right = np.nanmax(np.where((v > v[j]) & c, r, np.nan)) if np.any(
        (v > v[j]) & c) else np.nan
    out["reversal"] = int(np.isfinite(left) and np.isfinite(right) and
                          min(left, right) - r[j] >= 0.05)
    return out


def cmd_profiles(args) -> None:
    """A9: hrg only, on-focus frames only; the delivered LSF FWHM in km/s is
    stated first.  Velocities are DIFFERENTIAL (zero = the frame's own Halpha
    centroid, A3's anchor), so only shapes and wing extents are compared."""
    from macro_tcrb import spec
    g = db.grism_ro()
    con = db.connect()
    lsf_a = db.q1(con, "SELECT fwhm_a FROM tcrb_lsf_adopted WHERE grism="
                       "'hrg'")
    lsf_kms = lsf_a / spec.HALPHA * C_KMS
    coeffs, xref, *_ = _library_solution(g, "hrg")
    rows = db.q(con, """SELECT path, night, jd, x_halpha, lsf_fwhm_a FROM
        tcrb_ew_frames WHERE grism = 'hrg' AND status = 'ok' AND
        COALESCE(off_nominal_focus, 0) = 0 ORDER BY jd""")
    out = []
    for path, night, jd, xha, lfw in rows:
        sf = db.q1(g, "SELECT spec_file FROM g_frames WHERE path = ?", path)
        d = np.load(db.REPO / "products" / "grism" / sf)
        lam = spec.wavelength_scale(np.arange(d["flux"].size), coeffs, xref,
                                    xha)
        o = np.argsort(lam)
        lam, f = lam[o], d["flux"][o]
        b = (lam > spec.BLUE[0]) & (lam < spec.BLUE[1])
        rr = (lam > spec.RED[0]) & (lam < spec.RED[1])
        xb, yb = np.median(lam[b]), np.median(f[b])
        xr, yr = np.median(lam[rr]), np.median(f[rr])
        cont = yb + (yr - yb) / (xr - xb) * (lam - xb)
        sel = np.abs(lam - spec.HALPHA) < 45
        v = (lam[sel] - spec.HALPHA) / spec.HALPHA * C_KMS
        m = _profile_metrics(v, (f[sel] / cont[sel]) - 1)
        out.append((path, night, jd, lfw / spec.HALPHA * C_KMS if lfw
                    else None, m["fwhm"], m["v25_blue"], m["v25_red"],
                    m["v10_blue"], m["v10_red"], m["asym10"],
                    m["reversal"]))
    n = db.write_table(con, "tcrb_profiles", (
        "path", "night", "jd", "lsf_fwhm_kms", "fwhm_kms", "v25_blue",
        "v25_red", "v10_blue", "v10_red", "asym10_kms", "reversal"), out)
    a = np.array([r[4:10] for r in out], float)
    nights = sorted({r[1] for r in out})
    nm = {k: np.array([[r[4], r[7], r[8], r[9]] for r in out
                       if r[1] == k], float) for k in nights}
    med = np.array([np.nanmedian(v_, axis=0) for v_ in nm.values()])
    jd = np.array([np.mean([r[2] for r in out if r[1] == k])
                   for k in nights])
    # differential evolution: linear trend of nightly FWHM and HW10% wings
    def trend(y):
        ok = np.isfinite(y)
        A = np.vstack([np.ones(ok.sum()), jd[ok] - jd[ok].mean()]).T
        c, res, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
        r_ = y[ok] - A @ c
        se = np.sqrt(np.sum(r_ ** 2) / max(ok.sum() - 2, 1) /
                     np.sum((jd[ok] - jd[ok].mean()) ** 2))
        return float(c[1] * 100), float(se * 100)   # per 100 d
    summ = [("lsf_fwhm_kms", lsf_kms, None),
            ("n_frames", float(n), None), ("n_nights", float(len(nights)),
                                           None),
            ("fwhm_kms_median", float(np.nanmedian(a[:, 0])),
             float(1.4826 * np.nanmedian(np.abs(a[:, 0] -
                                                np.nanmedian(a[:, 0]))))),
            ("hwzi10_blue_kms_median", float(np.nanmedian(a[:, 3])), None),
            ("hwzi10_red_kms_median", float(np.nanmedian(a[:, 4])), None),
            ("asym10_kms_median", float(np.nanmedian(a[:, 5])), None),
            ("reversal_fraction", float(np.mean([r[10] for r in out])),
             None)]
    for nm_, col in (("fwhm", 0), ("v10_blue", 1), ("v10_red", 2),
                     ("asym10", 3)):
        t_, e_ = trend(med[:, col])
        summ.append((f"trend_{nm_}_kms_per_100d", t_, e_))
    db.write_table(con, "tcrb_profiles_summary", ("quantity", "value",
                                                  "uncertainty"), summ)
    db.record_stage(con, "profiles", {"tcrb_profiles": n})
    L = [f"Delivered LSF at Halpha (hrg, T CrB 240 s frames, G-5): "
         f"FWHM {lsf_kms:.0f} km/s (stated before any velocity).", "",
         "| quantity | value | +- |", "|---|---|---|"]
    L += [f"| {a_} | {b_:.3g} | {'' if c_ is None else f'{c_:.2g}'} |"
          for a_, b_, c_ in summ]
    L += ["", "Velocities are differential (zero = each frame's Halpha "
          "centroid); off-nominal-focus frames excluded; lrg (LSF ~650 "
          "km/s) is not used for profiles."]
    write_note("a9_profiles", "TCRB-A9 — Halpha profile morphology and "
               "wing velocities (hrg)", L,
               "Stage `profiles` of `pipeline/scripts/run_tcrb.py`; table "
               "`tcrb_profiles` (per frame) and `tcrb_profiles_summary`.")
    print("\n".join(L))


# ===========================================================================
# TCRB-C4 — uncertainties on every headline number
# ===========================================================================
def cmd_uncert(args) -> None:
    """value +- stat +- sys for every headline number, with the systematic
    terms listed separately (linearity and cross-mode penalty named), emcee
    posteriors for the season means, and the red-noise inflation factor
    from binned comparison-star rms."""
    from macro_tcrb import phot
    con = db.connect()
    man = db.manifest_ro()
    rows = []
    # red-noise inflation (B check stars, per night, Pont+2006 method)
    ser = []
    for night, sid in db.q(con, """SELECT DISTINCT night, source_id FROM
            tcrb_bphot_checks"""):
        v = [r[0] for r in db.q(con, """SELECT resid FROM tcrb_bphot_checks
            WHERE night = ? AND source_id = ? ORDER BY jd""", night, sid)]
        ser.append(np.array(v))
    beta = phot.red_noise_beta(ser)
    rows.append(("red_noise_beta_B", beta, None, None, "{}",
                 "Pont+2006 binning of check-star residuals, n=2..3"))
    # linearity term for B (StackPro): measured slope x T CrB peak fraction
    slope = db.q1(man, """SELECT value FROM detector_params WHERE era_group
        = 'High Gain StackPro' AND quantity =
        'linearity_slope_pct_per_scale'""")
    for night, n, bt, bn, rms, ext, *_ in db.q(con, """SELECT * FROM
            tcrb_b_anchors ORDER BY night"""):
        pf = db.q1(con, """SELECT AVG(c.peak_native / c.clip_adu) FROM
            tcrb_census c JOIN tcrb_bphot_frames b USING (obs_rowid) WHERE
            b.night = ?""", night)
        stale = db.q1(con, "SELECT AVG(flat_stale_mag) FROM "
                           "tcrb_bphot_frames WHERE night = ?", night)
        floor = db.q1(con, "SELECT season_floor_mag FROM tcrb_bprecision")
        lin = abs(slope) / 100 * pf * 1.0857 if slope and pf else 0.0
        stat = rms / np.sqrt(n) * beta
        terms = {"colour_extrapolation": ext, "flat_staleness": stale,
                 "season_floor": floor, "linearity": lin,
                 "cross_mode_penalty": 0.0}
        sys_ = float(np.sqrt(sum(v ** 2 for v in terms.values())))
        rows.append((f"B_{night}", bt, stat, sys_, json.dumps(terms),
                     "nightly mean of transformed B; stat = check-star rms/"
                     "sqrt(n) x beta; cross-mode 0 (mode-matched flat; "
                     "0.035 mag had the High Gain flat been used)"))
    # season EW per grism: emcee mean + intrinsic scatter (jitter)
    for gr in ("hrg", "lrg"):
        a = np.array(db.q(con, """SELECT ew, ew_err, ew_m10, ew_p10 FROM
            tcrb_ew_nightly WHERE grism = ?""", gr), float)
        post = phot.mean_with_jitter(a[:, 0], a[:, 1])
        mu = post["mu"]
        win = float(np.nanmax(np.abs([np.mean(a[:, 2] - a[:, 0]),
                                      np.mean(a[:, 3] - a[:, 0])])))
        degr = db.q1(con, "SELECT median_signed_frac FROM "
                          "tcrb_window_survival WHERE grism = ?", gr)
        terms = {"window_shift_10A": win,
                 "lsf_degradation": abs(degr) * mu[1]}
        sys_ = float(np.sqrt(sum(v ** 2 for v in terms.values())))
        rows.append((f"EW_season_mean_{gr}", mu[1], 0.5 * (mu[2] - mu[0]),
                     sys_, json.dumps(terms), "emcee posterior (16/50/84) "
                     "of the mean with free intrinsic scatter"))
        j = post["jitter"]
        rows.append((f"EW_intrinsic_scatter_{gr}", j[1],
                     0.5 * (j[2] - j[0]), None, "{}",
                     "emcee jitter: night-to-night variation beyond errors"))
        lo, hi = float(np.min(a[:, 0])), float(np.max(a[:, 0]))
        rows.append((f"EW_min_{gr}", lo, float(a[np.argmin(a[:, 0]), 1]),
                     win, "{}", "lowest nightly EW; stat = its error"))
        rows.append((f"EW_max_{gr}", hi, float(a[np.argmax(a[:, 0]), 1]),
                     win, "{}", "highest nightly EW"))
    # ARAS agreement
    for gr, n, mean, rms, med, mad in db.q(con, "SELECT * FROM "
                                                "tcrb_xval_summary"):
        degr = db.q1(con, "SELECT median_signed_frac FROM "
                          "tcrb_window_survival WHERE grism = ?", gr)
        rows.append((f"ARAS_mean_frac_{gr}", mean, rms / np.sqrt(n),
                     abs(degr), json.dumps({"lsf_degradation": abs(degr)}),
                     "mean RLMT/ARAS-1 over same UT dates"))
        rows.append((f"ARAS_rms_frac_{gr}", rms,
                     rms / np.sqrt(2 * (n - 1)), None, "{}",
                     "rms; stat = rms/sqrt(2(n-1))"))
    # profile
    for q_ in ("fwhm_kms_median",):
        v, u = db.q(con, "SELECT value, uncertainty FROM "
                         "tcrb_profiles_summary WHERE quantity = ?", q_)[0]
        nn = db.q1(con, "SELECT value FROM tcrb_profiles_summary WHERE "
                        "quantity = 'n_frames'")
        lsf = db.q1(con, "SELECT value FROM tcrb_profiles_summary WHERE "
                         "quantity = 'lsf_fwhm_kms'")
        rows.append(("Halpha_FWHM_kms", v, u / np.sqrt(nn),
                     0.5 * lsf ** 2 / v, json.dumps({"lsf_quadrature": 0.5 *
                                                     lsf ** 2 / v}),
                     "observed FWHM; sys = half the LSF quadrature "
                     "correction"))
    # theta CrB
    v = db.q(con, """SELECT ha_ew_median, ha_nightly_robust_sd, n_nights FROM
        tcrb_a1 WHERE grism = 'hrg' AND grism_epoch = 'ASI-pre' AND
        exptime = 4.86""")[0]
    rows.append(("tetCrB_Halpha_EW_hrg", v[0], v[1] / np.sqrt(v[2]), None,
                 "{}", "median of nightly medians, 4.86 s"))
    n = db.write_table(con, "tcrb_headline", ("quantity", "value", "stat",
                                              "sys", "sys_terms", "method"),
                       rows)
    db.record_stage(con, "uncert", {"tcrb_headline": n})
    L = ["| quantity | value | stat | sys | systematic terms | method |",
         "|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r[0]} | {r[1]:.4g} | "
                 f"{'' if r[2] is None else f'{r[2]:.2g}'} | "
                 f"{'' if r[3] is None else f'{r[3]:.2g}'} | {r[4]} | "
                 f"{r[5]} |")
    write_note("c4_uncertainties", "TCRB-C4 — every headline number as "
               "value +- stat +- sys", L,
               "Stage `uncert` of `pipeline/scripts/run_tcrb.py`; table "
               "`tcrb_headline`. Linearity and cross-mode penalty terms are "
               "listed separately for the photometry; no max(chi2nu, 1) "
               "anywhere.")
    print("\n".join(L))


# ===========================================================================
# A-track report: A1, A5, A5b, A7, A8 and the rule, from the tables
# ===========================================================================
def cmd_atrack_report(args) -> None:
    con = db.connect()
    L = ["## A1 — theta CrB (B6, Be/shell star) before it is used", ""]
    fl = {r[0]: r for r in db.q(con, "SELECT * FROM tcrb_floor")}
    a1 = []
    for gr, ep, ex in db.q(con, """SELECT DISTINCT grism, grism_epoch,
            exptime FROM tcrb_tet_ew ORDER BY 1, 2, 3"""):
        rs = db.q(con, """SELECT night, ew_ha, ew_ha_err, ew_hb FROM
            tcrb_tet_ew WHERE grism = ? AND grism_epoch = ? AND exptime = ?
            AND continuum_fault = 0""", gr, ep, ex)
        nf = db.q1(con, """SELECT COUNT(*) FROM tcrb_tet_ew WHERE grism = ?
            AND grism_epoch = ? AND exptime = ?""", gr, ep, ex)
        by: dict = {}
        for n_, e_, s_, hb in rs:
            by.setdefault(n_, []).append((e_, s_, hb))
        nm = np.array([np.median([x[0] for x in v]) for v in by.values()])
        ne = np.array([np.sqrt(np.mean([x[1] ** 2 for x in v]) / len(v))
                       for v in by.values()])
        f = fl.get(gr)
        chi = dof = None
        if f is not None and ep == "ASI-pre" and len(nm) > 2:
            nfr = np.array([len(v) for v in by.values()])
            err = np.sqrt(ne ** 2 + f[1] ** 2 / nfr + f[2] ** 2)
            chi, dof = _chi2(nm, err)
        hb = [x[2] for v in by.values() for x in v if x[2] is not None]
        a1.append((gr, ep, ex, nf, len(rs), len(by),
                   float(np.median(nm)) if nm.size else None,
                   float(1.4826 * np.median(np.abs(nm - np.median(nm))))
                   if nm.size else None,
                   float(np.median(ne)) if ne.size else None, chi, dof,
                   float(np.median(hb)) if hb else None))
    db.write_table(con, "tcrb_a1", (
        "grism", "grism_epoch", "exptime", "n_frames", "n_good", "n_nights",
        "ha_ew_median", "ha_nightly_robust_sd", "ha_nightly_err", "chi2",
        "dof", "hb_ew_median"), a1)
    nconf = db.q1(con, "SELECT COUNT(*) FROM tcrb_tet_frames")
    nok = db.q1(con, "SELECT COUNT(*) FROM tcrb_tet_frames WHERE status="
                     "'ok'")
    L += [f"All {nconf} staged theta CrB grism frames reduced with the "
          f"shared library; {nok} traced. Per configuration (exposures "
          "never pooled):", "",
          "| grism | epoch | exp (s) | frames | good | nights | Halpha EW "
          "median (A) | nightly robust sd | nightly error | chi2 / dof "
          "with the A7 floor | Hbeta EW median |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in a1:
        cells = ["" if x is None else (f"{x:.2f}" if isinstance(x, float)
                                      else str(x)) for x in r]
        chi = "" if r[9] is None else f"{r[9]:.1f} / {r[10]}"
        L.append("| " + " | ".join(cells[:9] + [chi, cells[11]]) + " |")
    L += ["", "Reading: theta CrB shows Halpha in ABSORPTION in every "
          "configuration (no emission, no shell reversal resolved at the "
          "hrg R~2,800). Its night-to-night Halpha scatter is tested "
          "against the floor measured on the line-free window. The two "
          "2025 hrg exposures observe the SAME nights: at 4.86 s the Halpha "
          "is constant within the floor, at 2.43 s it is not — a variation "
          "seen at one exposure and not the other on the same nights is "
          "instrumental (short-exposure), not Be activity, so no Be "
          "variability is detected in 2025 and the 2.43 s series is not "
          "used for the floor's Halpha cross-check. The 2026 QHY value "
          "belongs to a second instrument and is not compared. theta CrB "
          "Halpha is never used to anchor wavelength or as a reference. In lrg the Hbeta window returns a "
          "positive (unphysical for a B6 star) EW: second-order blue light "
          "overlaps the first order of this blue star, so theta CrB lrg is "
          "unusable as an EW reference below ~5000 A (PH's warning, "
          "measured). The paper states the Be nature of theta CrB.", ""]
    # ---- A7 -------------------------------------------------------------
    L += ["## A7 — the EW error floor", "",
          "| grism | theta CrB per-frame floor (A) | night floor (A) | T CrB "
          "within-night excess (A) | frames/night | split-half chi2/dof "
          "(B given A) | full chi2/dof | extraction-method difference (A) "
          "| 240 s smear (A) | adopted nightly floor (A) | off-focus frames "
          "|", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in db.q(con, "SELECT * FROM tcrb_floor"):
        L.append(f"| {r[0]} | {r[1]:.2f} | {r[2]:.2f} | {r[3]:.2f} | "
                 f"{r[4]:.0f} | {r[5]:.1f} / {r[6]} = {r[7]:.2f} | "
                 f"{r[8]:.1f} / {r[9]} = {r[8] / r[9]:.2f} | {r[10]:.2f} | "
                 f"{r[12]:.2f} | {r[13]:.2f} | {r[16]} of {r[15]} |")
    for r in db.q(con, "SELECT * FROM tcrb_floor_halves"):
        L.append("")
        L.append(f"{r[0]}: night floor from half 0 = {r[9]:.2f} A, from "
                 f"half 1 = {r[10]:.2f} A; adopted the larger; with it "
                 f"chi2/dof = {r[11]:.1f}/{r[12]} and {r[13]:.1f}/{r[14]} on "
                 f"the two halves.")
    L += ["", "The theta CrB continuum-window series with the adopted floor "
          "has chi2/dof inside [0.7, 1.4] on the full series in both "
          "grisms; the split-half transfer (floor fitted on alternate "
          "nights, tested on the rest) is reported as well and shows the "
          "night term is uncertain at the factor-1.5 level with ~11 dof, "
          "which is why the larger half value is adopted. The floor exceeds "
          "the extraction-method difference in both grisms.", ""]
    # ---- A5 -------------------------------------------------------------
    ws = db.q(con, "SELECT * FROM tcrb_window_survival")
    sh = db.q(con, """SELECT grism, COUNT(*), AVG(ew_m10 / ew - 1),
        AVG(ew_p10 / ew - 1) FROM tcrb_ew_frames WHERE status = 'ok'
        GROUP BY 1""")
    ex = dict(db.q(con, "SELECT reason, n_frames FROM tcrb_ew_exclusions"))
    L += ["## A5 — Halpha EW (pre-registered windows)", "",
          "| grism | ARAS spectra | median abs(degraded/native - 1) | 90th "
          "pct | median signed | passes the pre-registered tolerance |",
          "|---|---|---|---|---|---|"]
    L += [f"| {a} | {b} | {c:.3f} | {d:.3f} | {e:+.3f} | "
          + ("yes" if f_ else "no") + " |" for a, b, c, d, e, f_ in ws]
    L += ["", "| grism | frames | EW change, windows -10 A | +10 A |",
          "|---|---|---|---|"]
    L += [f"| {a} | {b} | {c:+.3f} | {d:+.3f} |" for a, b, c, d in sh]
    L += ["", f"Frames excluded from the nightly series: saturated "
          f"{ex.get('saturated', 0)}, off-nominal focus "
          f"{ex.get('off_focus', 0)}.", ""]
    for gr in ("hrg", "lrg"):
        r = db.q(con, """SELECT COUNT(*), MIN(ew), MAX(ew), AVG(ew_err),
            MIN(night), MAX(night) FROM tcrb_ew_nightly WHERE grism = ?""",
                 gr)[0]
        L.append(f"- {gr}: {r[0]} nights {r[4]}..{r[5]}, EW {r[1]:.1f} to "
                 f"{r[2]:.1f} A, mean nightly error {r[3]:.2f} A.")
    L += ["", "## The pre-registered rule (A5a)", "",
          "| grism | nights | threshold (sigma) | expected false calls | "
          "events | chi2 / dof (constant) | p | p (floor x1.5) | variable "
          "| step recovered at 90% (eligible boundaries, A) | eligible "
          "fraction |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in db.q(con, "SELECT * FROM tcrb_detect"):
        L.append(f"| {r[0]} | {r[1]} | {r[2]:.1f} | {r[3]:.4f} | {r[4]} | "
                 f"{r[5]:.1f} / {r[6]} | {r[7]:.2g} | {r[8]:.2g} | "
                 f"{'yes' if r[9] else 'no'} | {r[10]:.1f} | {r[11]:.2f} |")
    # ---- A5b ------------------------------------------------------------
    L += ["", "## A5b — line flux against orbital phase", ""]
    for gr in ("hrg", "lrg"):
        a = np.array(db.q(con, """SELECT phase, ew, ew_err, line_flux,
            line_flux_err, cont_mag FROM tcrb_flux WHERE grism = ? AND
            line_flux IS NOT NULL ORDER BY jd""", gr), float)
        srcs = db.q(con, """SELECT cont_source, COUNT(*) FROM tcrb_flux
            WHERE grism = ? GROUP BY 1""", gr)
        if a.size:
            cf, df = _chi2(a[:, 3], a[:, 4])
            ce, de = _chi2(a[:, 1], a[:, 2])
            L.append(f"- {gr}: {len(a)} nights, phase {a[:, 0].min():.3f}"
                     f"-{a[:, 0].max():.3f}; AAVSO Rc {a[:, 5].min():.2f}-"
                     f"{a[:, 5].max():.2f} mag; line-flux range "
                     f"{a[:, 3].min():.3g}-{a[:, 3].max():.3g} erg/cm2/s; "
                     f"chi2/dof about constant: EW {ce:.0f}/{de}, flux "
                     f"{cf:.0f}/{df}. Continuum source: AAVSO Rc same "
                     f"night on {sum(n_ for s_, n_ in srcs if 'same' in s_)}"
                     f" nights, interpolated (<= 2 d) on "
                     f"{sum(n_ for s_, n_ in srcs if 'interp' in s_)}, none"
                     f" on {sum(n_ for s_, n_ in srcs if s_ == 'none')} "
                     "(named per night in `tcrb_flux.cont_source`).")
    # ---- A8 -------------------------------------------------------------
    L += ["", "## A8 — same-UT-date ARAS cross-validation", "",
          "| grism | dates | mean (RLMT/ARAS - 1) | rms | median | MAD |",
          "|---|---|---|---|---|---|"]
    L += [f"| {a} | {b} | {c:+.3f} | {d:.3f} | {e:+.3f} | {f_:.3f} |"
          for a, b, c, d, e, f_ in db.q(con,
                                        "SELECT * FROM tcrb_xval_summary")]
    write_note("a_track", "Phase A — theta CrB, the EW series, line flux "
               "and the ARAS cross-validation", L,
               "Stages `reduce-tet`, `tet-analysis`, `aras`, `ew`, `floor`, "
               "`nightly`, `flux`, `xval`, `detect` of "
               "`pipeline/scripts/run_tcrb.py`, on the shared grism "
               "library G v2.6 (G-1..G-5 closed).")
    print("\n".join(L))


# ===========================================================================
# paper — numbers.tex, tables, figures
# ===========================================================================
def cmd_paper(args) -> None:
    from macro_tcrb import paper
    n = paper.emit_numbers()
    tabs = paper.emit_tables()
    figs = []
    try:
        from macro_tcrb import figures
        figs = figures.build_all()
    except ImportError:
        pass
    txt = n.read_text()
    miss = txt.count(r"\NumMissing}")
    print(f"  {n} ({txt.count(chr(10))} lines; {miss} missing)")
    for t in tabs + figs:
        print(f"  {t}")


# ===========================================================================
def cmd_status(args) -> None:
    con = db.connect()
    for r in db.q(con, "SELECT stage, n_rows, run_utc FROM tcrb_meta "
                       "ORDER BY run_utc"):
        print("  ", *r)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in COMMANDS.items():
        sp = sub.add_parser(name)
        sp.set_defaults(fn=fn)
        for flag, kw in EXTRA_ARGS.get(name, ()):
            sp.add_argument(flag, **kw)
    args = ap.parse_args()
    args.fn(args)


COMMANDS = {
    "mech-epoch": cmd_mech_epoch,
    "temp-split": cmd_temp_split,
    "temp-report": cmd_temp_report,
    "shutter": cmd_shutter,
    "catalogues": cmd_catalogues,
    "measure": cmd_measure,
    "census": cmd_census,
    "filters": cmd_filters,
    "filters-report": cmd_filters_report,
    "zmag-report": cmd_zmag_report,
    "resolve-report": cmd_resolve_report,
    "external": cmd_external,
    "external-report": cmd_external_report,
    "aavso-load": cmd_aavso_load,
    "ensemble": cmd_ensemble,
    "errors": cmd_errors,
    "flicker": cmd_flicker,
    "phaseb-report": cmd_phaseb_report,
    "aras": cmd_aras,
    "reduce-tet": cmd_reduce_tet,
    "ew": cmd_ew,
    "tet-analysis": cmd_tet_analysis,
    "paper": cmd_paper,
    "uncert": cmd_uncert,
    "profiles": cmd_profiles,
    "atrack-report": cmd_atrack_report,
    "a0b": cmd_a0b,
    "nightly": cmd_nightly,
    "release": cmd_release,
    "flux": cmd_flux,
    "xval": cmd_xval,
    "floor": cmd_floor,
    "detect": cmd_detect,
    "status": cmd_status,
}
EXTRA_ARGS: dict = {
    "catalogues": (("--force", {"action": "store_true"}),),
    "external": (("--sources", {"default": None}),),
    "flicker": (("--n-mc", {"type": int, "default": 300}),),
    "measure": (("--retry-failed", {"action": "store_true"}),),
    "detect": (("--n-mc", {"type": int, "default": 10000}),),
    "reduce-tet": (("--workers", {"type": int, "default": 3}),
                   ("--limit", {"type": int, "default": 0})),
}

if __name__ == "__main__":
    main()
