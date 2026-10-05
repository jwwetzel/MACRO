#!/usr/bin/env python
"""Was the camera turned, or was the picture? — the hot-pixel flip test.

THE QUESTION (finding F-3 / TE.F1, plan review 2026-10-03)
----------------------------------------------------------
Across the 2025 monsoon shutdown the sky rotated by 179.3 degrees on the
ASI camera's frames, and the ``FLIPSTAT`` header card changed from
``'Flip/Mirror'`` to blank.  Two different things can do that:

* a SOFTWARE flip — MaxIm stopped (or started) mirroring the array it
  saves.  The sensor never moved; every file after the change is the same
  detector written out rotated by 180 degrees;
* a PHYSICAL rotation — somebody turned the camera on the telescope.

They need different handling.  After a software flip, a dark frame taken
before the change is still a perfect description of the sensor — once it is
rotated by 180 degrees — and is WRONG, pixel for pixel, if applied as
stored.  After a physical rotation, a dark is valid exactly as stored, and
only the sky-side calibrations (flats, dust, grism trace) are lost.

The sky cannot tell the two apart: both rotate it.  The sensor can.  A hot
pixel is a defect in the silicon.  It sits at fixed coordinates on the
chip, so in the saved files it stays put under a physical rotation and
moves to ``(W-1-x, H-1-y)`` under a software flip.  The telescope
engineer's memo asked for exactly this test: *"a hot-pixel map showing the
same pixel coordinates before and after 2025 monsoon would make the flip
physical-only and leave darks valid (flats still not)"*.

WHAT THIS SCRIPT DOES
---------------------
For one mechanical boundary (default: ``ASI:2025-10-11``) it builds a
hot-pixel map on each side and counts how many hot pixels coincide under
each of four array transforms — identity, 180-degree rotation, left-right
mirror, up-down mirror — and under a set of arbitrary shifts.

It does this TWICE, on independent evidence:

``dark``   long dark exposures and master darks (the direct measurement);
``light``  science frames of many different fields.  A star is in a
           different place in every field; a hot pixel is in the same place
           in all of them.  This arm shares no file, no exposure type and
           no reduction step with the dark arm, so agreement between the
           two is not an artefact of how the darks were made.

THE CONTROLS, without which the counts mean nothing:

* **containment, not coincidence** — the two sides are never taken under
  the same conditions (sensor temperature, sky brightness), so a fixed
  threshold selects a different subset of the same hot pixels on each.
  The statistic asks: of the pixels CERTAINLY hot on one side, how many are
  at least PLAUSIBLY hot on the other?  (See ``LOOSE_SIGMA``.)
* **shifted maps** — the same comparison after displacing one map by an
  arbitrary offset.  This is the chance-coincidence rate, measured rather
  than computed, because hot pixels cluster (columns, amplifier corners)
  and a uniform-density formula would understate it.
* **split halves** — each side's own frames divided in two and compared
  with each other.  This is the ceiling: how well a hot-pixel map
  reproduces when NOTHING changed.  A cross-boundary match is judged
  against this, not against 100%.

It also reads the bias level on both sides (DE.F5: "show the bias level on
both sides of the 2025 gap").

OUTPUTS
-------
* ``mech_flip_test``    one row per (boundary, arm, transform): hot pixels
                        on each side, coincidences, fraction.
* ``mech_flip_frames``  every frame read, with its side, arm and hot count.
* ``mech_flip_bias``    bias level per frame on each side.
* ``docs/pipeline/figures/s0b/s0b_mech_flip.png``

All three tables are written into the manifest in one transaction.  The
archive is opened READ-ONLY.

USAGE
-----
    /opt/miniconda3/envs/rlmt-checks/bin/python \\
        pipeline/scripts/check_mech_flip.py
    ... --boundary QHY600:2026-03-23     # any rotation boundary
    ... --n-frames 8                     # fewer frames per side (faster)
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
import warnings
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PIPELINE_ROOT))

from macro_core import staging as stg                        # noqa: E402

REPO_ROOT = PIPELINE_ROOT.parent
DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
FIG_PATH = (REPO_ROOT / "docs" / "pipeline" / "figures" / "s0b"
            / "s0b_mech_flip.png")

FLIP_CODE_VERSION = "mech-flip v1.0 (2026-10-03)"

#: The boundary TE.F1 asked about: the ASI camera's first epoch after the
#: 2025 monsoon shutdown.  It begins on 2025-10-11 — the night of the first
#: post-shutdown calibration frames, where FLIPSTAT first reads blank —
#: four nights before the first plate-solved science frame showed the
#: rotated sky.
DEFAULT_BOUNDARY = "ASI:2025-10-11"

#: Frames read per side and arm.  Twelve is enough that a pixel must be hot
#: in nine of them to count; reading more costs minutes on the archive
#: drive and changes nothing.
N_FRAMES = 12

#: A pixel is a hot-pixel CANDIDATE in one frame when it stands this many
#: robust sigma above its 3x3 median-filtered neighbourhood.  The median
#: filter removes isolated pixels and keeps everything smooth, so the
#: residual of a hot pixel is its full excess; a star, being wider than one
#: pixel, mostly survives the filter and leaves little.
CANDIDATE_SIGMA = 8.0

#: ... and it is a HOT PIXEL when it is a candidate in at least this
#: fraction of the frames of its side.  This is what removes stars and
#: cosmic rays from the light arm: they are never in the same place twice.
PERSISTENCE = 0.7

#: The LOOSE map: a lower threshold and a lower persistence.  The two sides
#: of a boundary are never taken under the same conditions — the sensor was
#: at -10 C before the 2025 shutdown and partly at 0 C after, the sky
#: brightness differs frame to frame — so "above 8 sigma in 70% of frames"
#: selects a different SUBSET of the same hot pixels on each side, and
#: comparing strict map with strict map understates the agreement (it gave
#: 0.40 on the first run for a match that is in fact near-complete).  The
#: statistic is therefore CONTAINMENT: of the pixels that are certainly hot
#: on one side (strict), what fraction is at least plausibly hot on the
#: other (loose)?  The shifted control is computed the same way, so the
#: looser map's higher density is paid for in the chance rate it is judged
#: against.
LOOSE_SIGMA = 5.0
LOOSE_PERSISTENCE = 0.3

#: Exposure-time window (s).  Dark current scales with time: below a minute
#: few hot pixels clear the threshold.
DARK_MIN_EXPTIME = 60.0
LIGHT_EXPTIME = (60.0, 400.0)

#: Number of arbitrary shifts for the chance-coincidence control, and the
#: seed that makes them the same on every run.
N_SHIFTS = 40
SHIFT_SEED = 20261003

#: A hot-pixel map with more pixels than this is drawn as a small patch;
#: a sparser one is drawn whole (see draw_figure).
PATCH_IF_MORE_THAN = 5000

#: The array transforms compared.  Each maps a boolean hot-pixel map of the
#: AFTER side into the BEFORE side's coordinates.
TRANSFORMS = {
    "identity": lambda m: m,
    "rot180": lambda m: m[::-1, ::-1],
    "mirror_lr": lambda m: m[:, ::-1],
    "mirror_ud": lambda m: m[::-1, :],
}


# ---------------------------------------------------------------------------
# Pure logic (unit-tested in pipeline/tests/test_inventory.py)
# ---------------------------------------------------------------------------
def residual_in_sigma(image: np.ndarray) -> np.ndarray | None:
    """One-pixel-wide residual of an image, in units of its own noise.

    ``image - median_filter(image, 3x3)`` isolates anything one pixel wide:
    the median filter removes isolated pixels and keeps everything smooth,
    so a hot pixel keeps its full excess while a star, being wider than one
    pixel, mostly cancels.  The noise of that residual is estimated
    robustly, so the hot pixels themselves do not inflate it.

    The scale ladder exists because MASTER darks are not raw frames: a
    median-combined, integer-valued master has a residual that is exactly
    zero in more than half its pixels, so 1.4826 x MAD is 0.  The next rung
    is half the 16th–84th percentile width; the last is the standard
    deviation of the residual with its top and bottom 0.1% removed.
    Returns ``None`` only for an image with no variation at all.
    """
    from scipy.ndimage import median_filter
    img = np.asarray(image, dtype=np.float32)
    resid = img - median_filter(img, size=3, mode="nearest")
    sigma = 1.4826 * float(np.median(np.abs(resid - np.median(resid))))
    if not (np.isfinite(sigma) and sigma > 0):
        p16, p84 = np.percentile(resid, [16.0, 84.0])
        sigma = 0.5 * float(p84 - p16)
    if not (np.isfinite(sigma) and sigma > 0):
        lo, hi = np.percentile(resid, [0.1, 99.9])
        core = resid[(resid >= lo) & (resid <= hi)]
        sigma = float(core.std()) if core.size else 0.0
    if not (np.isfinite(sigma) and sigma > 0):
        return None
    return resid / sigma


def candidate_mask(image: np.ndarray,
                   n_sigma: float = CANDIDATE_SIGMA) -> np.ndarray:
    """Boolean map of single-frame hot-pixel candidates: pixels more than
    ``n_sigma`` above their neighbourhood (see :func:`residual_in_sigma`).

    Only POSITIVE outliers: a dead pixel is a defect too, but in a sky
    frame every star's shoulder is a negative residual, and one sign is
    enough to answer the question.  An image with no noise scale has no
    outliers.
    """
    z = residual_in_sigma(image)
    if z is None:
        return np.zeros(np.asarray(image).shape, dtype=bool)
    return z > n_sigma


def persistent_mask(masks: list[np.ndarray],
                    persistence: float = PERSISTENCE) -> np.ndarray:
    """Pixels flagged in at least ``persistence`` of the frames.

    Raises ``ValueError`` on an empty list or mismatched shapes — a map
    built from frames of two different geometries would be nonsense, and
    silently cropping would hide that.
    """
    if not masks:
        raise ValueError("persistent_mask needs at least one frame")
    shape = masks[0].shape
    if any(m.shape != shape for m in masks):
        raise ValueError("frames of different geometry in one hot-pixel map")
    need = int(np.ceil(persistence * len(masks)))
    count = np.zeros(shape, dtype=np.int16)
    for m in masks:
        count += m
    return count >= max(need, 1)


def containment(strict_a: np.ndarray, loose_a: np.ndarray,
                strict_b: np.ndarray, loose_b: np.ndarray
                ) -> tuple[int, int]:
    """``(n_found, n_asked)`` for two sides already in ONE coordinate frame.

    ``n_asked`` is the number of certainly-hot pixels on both sides
    together; ``n_found`` is how many of them are at least plausibly hot on
    the OTHER side.  Symmetric by construction, so neither side's
    temperature or exposure decides the answer.
    """
    found = int(np.count_nonzero(strict_a & loose_b)) \
        + int(np.count_nonzero(strict_b & loose_a))
    asked = int(strict_a.sum()) + int(strict_b.sum())
    return found, asked


def compare_maps(strict_b: np.ndarray, loose_b: np.ndarray,
                 strict_a: np.ndarray, loose_a: np.ndarray,
                 n_shifts: int = N_SHIFTS, seed: int = SHIFT_SEED
                 ) -> list[dict]:
    """Containment under every transform, plus the shifted control.

    ``*_b`` are the BEFORE side's strict and loose maps, ``*_a`` the AFTER
    side's.  Each transform is applied to the after side.  The control
    rolls the after side by ``n_shifts`` arbitrary offsets of at least 16
    pixels on both axes (wrap-around keeps the pixel count and the column
    structure): whatever is "found" then is found by accident, and its
    mean and maximum are the chance rate — measured, not computed, because
    hot pixels cluster and a uniform-density formula would understate it.
    """
    def _row(name, found, asked):
        return {"transform": name, "n_before": int(strict_b.sum()),
                "n_after": int(strict_a.sum()), "n_match": found,
                "frac": found / asked if asked else 0.0}

    rows = []
    for name, fn in TRANSFORMS.items():
        rows.append(_row(name, *containment(strict_b, loose_b,
                                            fn(strict_a), fn(loose_a))))
    rng = np.random.default_rng(seed)
    ny, nx = strict_b.shape
    founds, asked = [], 0
    for _ in range(n_shifts):
        dy = int(rng.integers(16, max(17, ny - 16)))
        dx = int(rng.integers(16, max(17, nx - 16)))
        f, asked = containment(
            strict_b, loose_b,
            np.roll(strict_a, (dy, dx), (0, 1)),
            np.roll(loose_a, (dy, dx), (0, 1)))
        founds.append(f)
    rows.append(_row("shift_control_mean", float(np.mean(founds)), asked))
    rows.append(_row("shift_control_max", int(np.max(founds)), asked))
    return rows


def verdict_of(rows: list[dict], ceiling: float) -> str:
    """Read one arm's comparison as a sentence.  PURE.

    ``ceiling`` is the split-half containment of a hot-pixel map with
    nothing changed.  A transform "explains" the boundary when its
    containment fraction is at least half that ceiling AND at least ten
    times the worst shifted-control fraction.  Exactly one of identity / rot180 may
    do so; anything else is reported as inconclusive rather than rounded
    to an answer.
    """
    by = {r["transform"]: r["frac"] for r in rows}
    floor = max(10.0 * by["shift_control_max"], 0.5 * ceiling)
    winners = [t for t in ("identity", "rot180", "mirror_lr", "mirror_ud")
               if by[t] >= floor and by[t] > 0]
    if winners == ["rot180"]:
        return "software_flip"
    if winners == ["identity"]:
        return "physical_rotation"
    return "inconclusive"


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def read_image(path: str) -> np.ndarray:
    """The 2-D image of a FITS file (first HDU that holds one), as float32.

    astropy is imported here so that the pure functions above can be tested
    without it.  Header verification warnings are silenced: ~20k archive
    files carry a malformed CONTINUE card that has nothing to do with their
    pixels.
    """
    from astropy.io import fits
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with fits.open(path, memmap=False, ignore_missing_simple=True) as hdul:
            for hdu in hdul:
                data = getattr(hdu, "data", None)
                if data is not None and getattr(data, "ndim", 0) == 2:
                    return np.asarray(data, dtype=np.float32)
    raise ValueError(f"no 2-D image in {path}")


def epochs_around(con: sqlite3.Connection, boundary: str) -> tuple[str, str]:
    """``(before, after)`` mechanical-epoch ids for a boundary: ``after`` is
    the epoch named, ``before`` the same camera's previous epoch."""
    row = con.execute("SELECT camera, first_night FROM mech_epoch "
                      "WHERE mech_epoch = ?", (boundary,)).fetchone()
    if row is None:
        raise SystemExit(f"ERROR: no mechanical epoch {boundary!r}; run "
                         "build_s0b_inventory.py first")
    prev = con.execute(
        "SELECT mech_epoch FROM mech_epoch WHERE camera = ? "
        "AND first_night < ? ORDER BY first_night DESC LIMIT 1",
        row).fetchone()
    if prev is None:
        raise SystemExit(f"ERROR: {boundary} is the camera's first epoch — "
                         "there is no 'before' side")
    return prev[0], boundary


def _spread(df: pd.DataFrame, n: int, by: list[str]) -> pd.DataFrame:
    """Up to ``n`` rows spread evenly through ``df``, at most one per
    distinct value of ``by`` (so no night or target is sampled twice)."""
    one = df.sort_values(["night", "path"]).drop_duplicates(by)
    if len(one) <= n:
        return one
    idx = np.unique(np.linspace(0, len(one) - 1, n).round().astype(int))
    return one.iloc[idx]


def select_frames(con: sqlite3.Connection, epoch: str, arm: str,
                  n: int) -> pd.DataFrame:
    """The frames of one side and one arm, deterministically chosen.

    ``dark``: raw darks and master darks of the epoch's DETECTOR-side
    geometry at >= 60 s, one per (night, exposure).  ``light``: canonical
    raw science frames at 60–400 s in g/r/i, one per target, so that no
    star field repeats.  Both restricted to the modal image geometry of the
    epoch's science frames — a hot-pixel map cannot mix binnings — and to
    nights placed in the epoch WITH CERTAINTY: a frame from the gap before
    a rotation-dated boundary could belong to either side, and one such
    frame on the wrong side would dilute exactly the contrast being
    measured.
    """
    geom = con.execute(
        """SELECT CAST(f.naxis1 AS INT), CAST(f.naxis2 AS INT), COUNT(*)
           FROM frames f JOIN frame_mech_epoch m ON m.obs_rowid = f.obs_rowid
           WHERE m.mech_epoch = ? AND f.is_canonical = 1 AND f.error IS NULL
             AND f.tree = 'rawimage' GROUP BY 1, 2 ORDER BY 3 DESC LIMIT 1""",
        (epoch,)).fetchone()
    if geom is None:
        return pd.DataFrame()
    base = """SELECT f.obs_rowid, f.path, f.night, f.exptime, f.target_key,
                     f.imagetyp, f.set_temp
              FROM frames f JOIN frame_mech_epoch m
                ON m.obs_rowid = f.obs_rowid
              WHERE m.mech_epoch = ? AND m.epoch_certain = 1
                AND f.is_canonical = 1
                AND f.error IS NULL AND CAST(f.naxis1 AS INT) = ?
                AND CAST(f.naxis2 AS INT) = ? """
    if arm == "dark":
        df = pd.read_sql_query(
            base + "AND f.obs_rowid IN (SELECT obs_rowid FROM calib_frames "
                   "WHERE kind = 'dark') AND f.exptime >= ?",
            con, params=(epoch, geom[0], geom[1], DARK_MIN_EXPTIME))
        return _spread(df, n, ["night", "exptime"])
    df = pd.read_sql_query(
        base + "AND f.tree = 'rawimage' AND f.imagetyp LIKE 'Light%' "
               "AND f.filter IN ('g', 'r', 'i') AND f.target_key IS NOT NULL "
               "AND f.exptime BETWEEN ? AND ? "
               "AND f.obs_rowid NOT IN (SELECT obs_rowid FROM calib_frames)",
        con, params=(epoch, geom[0], geom[1], *LIGHT_EXPTIME))
    return _spread(df, n, ["target_key"])


def side_map(frames: pd.DataFrame, archive_root: str, label: str
             ) -> tuple[dict | None, list[dict]]:
    """Read one side's frames → ``(maps, per-frame log)``.

    ``maps`` holds a ``(strict, loose)`` pair of boolean hot-pixel maps for
    the whole side (``"full"``) and for its two halves (``"h1"``, ``"h2"``
    — alternate frames, so both halves span the epoch), or ``None`` when
    fewer than four frames could be read.
    """
    strict, loose, log = [], [], []
    for _, fr in frames.iterrows():
        path = stg.abs_archive_path(archive_root, fr["path"])
        entry = {"path": fr["path"], "night": fr["night"],
                 "exptime": fr["exptime"], "n_candidates": None,
                 "error": None}
        try:
            z = residual_in_sigma(read_image(path))
        except Exception as e:                   # noqa: BLE001 — recorded
            entry["error"] = f"{type(e).__name__}: {e}"[:200]
            log.append(entry)
            continue
        if z is None:
            entry["error"] = "image has no variation (no noise scale)"
            log.append(entry)
            continue
        if strict and z.shape != strict[0].shape:
            entry["error"] = f"geometry {z.shape} != {strict[0].shape}"
            log.append(entry)
            continue
        strict.append(z > CANDIDATE_SIGMA)
        loose.append(z > LOOSE_SIGMA)
        entry["n_candidates"] = int(strict[-1].sum())
        log.append(entry)
        print(f"[flip]   {label}: {fr['path']}  "
              f"{entry['n_candidates']:,} candidates", flush=True)
    if len(strict) < 4:
        return None, log

    def _pair(sl):
        return (persistent_mask(strict[sl], PERSISTENCE),
                persistent_mask(loose[sl], LOOSE_PERSISTENCE))
    return {"full": _pair(slice(None)), "h1": _pair(slice(0, None, 2)),
            "h2": _pair(slice(1, None, 2))}, log


def bias_levels(con: sqlite3.Connection, epoch: str, archive_root: str,
                n: int) -> list[dict]:
    """Median level of up to ``n`` bias frames (raw and master) of an epoch.

    The central 512x512 pixels only: a bias has no signal anywhere, the
    centre avoids amplifier-glow corners, and reading a crop of the array
    is all the question needs.
    """
    df = pd.read_sql_query(
        """SELECT c.path, c.night, c.is_master FROM calib_frames c
           WHERE c.mech_epoch = ? AND c.epoch_certain = 1
             AND c.kind = 'bias'
           ORDER BY c.is_master, c.night, c.path""", con, params=(epoch,))
    rows = []
    for is_master, grp in df.groupby("is_master"):
        for _, fr in _spread(grp, n, ["night"]).iterrows():
            try:
                img = read_image(stg.abs_archive_path(archive_root,
                                                      fr["path"]))
                cy, cx = img.shape[0] // 2, img.shape[1] // 2
                crop = img[cy - 256:cy + 256, cx - 256:cx + 256]
                rows.append({"mech_epoch": epoch, "path": fr["path"],
                             "night": fr["night"],
                             "is_master": int(is_master),
                             "median_adu": float(np.median(crop)),
                             "mad_adu": float(np.median(np.abs(
                                 crop - np.median(crop)))), "error": None})
            except Exception as e:               # noqa: BLE001 — recorded
                rows.append({"mech_epoch": epoch, "path": fr["path"],
                             "night": fr["night"],
                             "is_master": int(is_master), "median_adu": None,
                             "mad_adu": None,
                             "error": f"{type(e).__name__}: {e}"[:200]})
    return rows


# ---------------------------------------------------------------------------
# The figure
# ---------------------------------------------------------------------------
def draw_figure(tests: pd.DataFrame, maps: dict, boundary: str,
                out: Path) -> Path:
    """Two panels per arm: the match fractions, and a patch of the two
    hot-pixel maps laid over each other so the eye can check the bars."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from macro_core import plotstyle as ps

    arms = [a for a in ("dark", "light") if a in maps]
    out.parent.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(ps.STYLE):
        fig, axes = plt.subplots(len(arms), 2,
                                 figsize=(ps.COL_DOUBLE * 1.45,
                                          3.9 * len(arms)), squeeze=False)
        for i, arm in enumerate(arms):
            sub = tests[tests["arm"] == arm].set_index("transform")
            ax = axes[i][0]
            names = ["identity", "rot180", "mirror_lr", "mirror_ud",
                     "shift_control_max"]
            labels = ["as stored", "rotated 180", "mirror L-R",
                      "mirror U-D", "shifted\n(chance)"]
            vals = [sub.loc[n, "frac"] for n in names]
            colors = [ps.ACCENT, ps.BAD, ps.FAINT, ps.FAINT, ps.MUTED]
            bars = ax.bar(labels, vals, color=colors, edgecolor=ps.INK,
                          linewidth=0.5)
            for b, n in zip(bars, names):
                ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                        f"{sub.loc[n, 'n_match']:,.0f}", ha="center",
                        va="bottom", fontsize=7)
            ceil_b = sub.loc["split_half_before", "frac"]
            ceil_a = sub.loc["split_half_after", "frac"]
            ax.axhline(ceil_b, **ps.reference_kw(color=ps.GOOD))
            ax.axhline(ceil_a, **ps.reference_kw(color=ps.GOOD, style=":"))
            ax.text(0.99, max(ceil_a, ceil_b), "split-half ceiling "
                    f"(before {ceil_b:.2f}, after {ceil_a:.2f})",
                    transform=ax.get_yaxis_transform(), ha="right",
                    va="bottom", fontsize=7, color=ps.GOOD)
            ax.set_ylabel("hot pixels found on other side")
            ax.set_title(
                f"({'ac'[i]}) {arm} frames: {int(sub.loc['identity', 'n_before']):,}"
                f" hot pixels before, {int(sub.loc['identity', 'n_after']):,} "
                f"after")
            ax.tick_params(axis="x", labelsize=7)

            # The maps themselves, laid over each other.  A dense map (the
            # dark arm: tens of thousands of hot pixels) is shown as a
            # 160x160 pixel patch at the array centre; a sparse one (the
            # light arm: hundreds) as the whole array, or the patch would
            # be empty.
            before, after = maps[arm]
            ax2 = axes[i][1]
            ny, nx = before.shape
            if int(before.sum()) > PATCH_IF_MORE_THAN:
                h = 80
                sl = (slice(ny // 2 - h, ny // 2 + h),
                      slice(nx // 2 - h, nx // 2 + h))
                y0, x0, what = ny // 2 - h, nx // 2 - h, \
                    "a 160 x 160 px patch at the array centre"
                sizes = (38, 22, 34)
            else:
                sl = (slice(0, ny), slice(0, nx))
                y0, x0, what = 0, 0, "the whole array"
                sizes = (16, 9, 14)
            for mask, kw, lab in (
                    (before[sl], dict(marker="o", s=sizes[0],
                                      facecolors="none", edgecolors=ps.INK,
                                      linewidths=0.7), "before"),
                    (after[sl], dict(marker="x", s=sizes[1], color=ps.ACCENT,
                                     linewidths=0.8), "after, as stored"),
                    (after[::-1, ::-1][sl],
                     dict(marker="+", s=sizes[2], color=ps.BAD,
                          linewidths=0.8), "after, rotated 180")):
                py, px = np.nonzero(mask)
                ax2.scatter(px + x0, py + y0, label=f"{lab} ({len(px):,})",
                            **kw)
            ax2.set_xlim(x0, x0 + (sl[1].stop - sl[1].start))
            ax2.set_ylim(y0, y0 + (sl[0].stop - sl[0].start))
            ax2.set_aspect("equal")
            ax2.legend(fontsize=6.5, loc="upper center", ncol=3,
                       bbox_to_anchor=(0.5, -0.16), frameon=False)
            ax2.set_xlabel("array x (pixel)")
            ax2.set_ylabel("array y (pixel)")
            ax2.set_title(f"({'bd'[i]}) hot pixels in {what}")
        fig.suptitle(f"Hot-pixel test of the {boundary} boundary: "
                     "did the sensor move, or the saved array?", fontsize=10)
        fig.tight_layout()
        fig.savefig(out, dpi=ps.WEB_DPI)
        plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Atomic write
# ---------------------------------------------------------------------------
FLIP_TABLES = ("mech_flip_test", "mech_flip_frames", "mech_flip_bias",
               "mech_flip_meta")


def write_tables(manifest: Path, tables: dict[str, pd.DataFrame]) -> None:
    """Swap the result tables into the manifest in one transaction (the
    S0b discipline: build under temporary names, then DROP + RENAME)."""
    assert set(tables) == set(FLIP_TABLES)
    with closing(sqlite3.connect(manifest, timeout=600.0)) as con:
        con.execute("PRAGMA busy_timeout = 600000")
        for name, frame in tables.items():
            con.execute(f"DROP TABLE IF EXISTS {name}_flip_tmp")
            frame.to_sql(f"{name}_flip_tmp", con, index=False)
        con.execute("BEGIN")
        for name in tables:
            con.execute(f"DROP TABLE IF EXISTS {name}")
            con.execute(f"ALTER TABLE {name}_flip_tmp RENAME TO {name}")
        con.execute("COMMIT")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Hot-pixel test of a mechanical boundary: decide whether "
                    "a 180-degree change of sky orientation was a software "
                    "flip of the saved array or a physical rotation of the "
                    "camera, on dark frames and (independently) on science "
                    "frames.  Reads the archive read-only; writes three "
                    "tables to the manifest and one figure.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--archive-root", default=stg.DEFAULT_ARCHIVE_ROOT)
    p.add_argument("--boundary", default=DEFAULT_BOUNDARY,
                   help="mech_epoch id whose START is the boundary tested")
    p.add_argument("--n-frames", type=int, default=N_FRAMES,
                   help="frames read per side and arm")
    p.add_argument("--figure", type=Path, default=FIG_PATH)
    p.add_argument("--dry-run", action="store_true",
                   help="list the frames that would be read and stop")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    with closing(sqlite3.connect(f"file:{args.manifest}?mode=ro", uri=True,
                                 timeout=600.0)) as con:
        before_id, after_id = epochs_around(con, args.boundary)
        chosen = {(arm, side): select_frames(con, eid, arm, args.n_frames)
                  for arm in ("dark", "light")
                  for side, eid in (("before", before_id),
                                    ("after", after_id))}
        bias = (bias_levels(con, before_id, args.archive_root, 6)
                + bias_levels(con, after_id, args.archive_root, 6)) \
            if not args.dry_run else []
    print(f"[flip] boundary {args.boundary}: before = {before_id}, "
          f"after = {after_id}")
    for (arm, side), df in chosen.items():
        print(f"[flip]   {arm:5} {side:6}: {len(df)} frames")
    if args.dry_run:
        for (arm, side), df in chosen.items():
            for p in df["path"]:
                print(f"    {arm} {side} {p}")
        return 0

    test_rows, frame_rows, maps = [], [], {}
    for arm in ("dark", "light"):
        sides = {}
        for side in ("before", "after"):
            sides[side], log = side_map(chosen[(arm, side)],
                                        args.archive_root, f"{arm}/{side}")
            for r in log:
                frame_rows.append({"boundary": args.boundary, "arm": arm,
                                   "side": side, **r})
        if sides["before"] is None or sides["after"] is None:
            print(f"[flip] {arm}: fewer than 4 readable frames on one side "
                  "— arm skipped (recorded as such)")
            test_rows.append({"boundary": args.boundary, "arm": arm,
                              "transform": "(arm skipped: too few frames)",
                              "n_before": None, "n_after": None,
                              "n_match": None, "frac": None,
                              "verdict": "not_run"})
            continue
        (sb, lb), (sa, la) = sides["before"]["full"], sides["after"]["full"]
        rows = compare_maps(sb, lb, sa, la)
        # Split-half ceilings: each side against itself, same statistic.
        ceil = {}
        for side in ("before", "after"):
            (s1, l1), (s2, l2) = sides[side]["h1"], sides[side]["h2"]
            found, asked = containment(s1, l1, s2, l2)
            ceil[side] = found / asked if asked else 0.0
            rows.append({"transform": f"split_half_{side}",
                         "n_before": int(s1.sum()), "n_after": int(s2.sum()),
                         "n_match": found, "frac": ceil[side]})
        verdict = verdict_of(rows, min(ceil.values()))
        for r in rows:
            test_rows.append({"boundary": args.boundary, "arm": arm, **r,
                              "verdict": verdict})
        maps[arm] = (sb, sa)
        print(f"[flip] {arm}: " + "; ".join(
            f"{r['transform']} {r['frac']:.3f}" for r in rows)
              + f"  ->  {verdict}")

    tests = pd.DataFrame(test_rows, columns=[
        "boundary", "arm", "transform", "n_before", "n_after", "n_match",
        "frac", "verdict"])
    frames = pd.DataFrame(frame_rows, columns=[
        "boundary", "arm", "side", "path", "night", "exptime",
        "n_candidates", "error"])
    bias_df = pd.DataFrame(bias, columns=[
        "mech_epoch", "path", "night", "is_master", "median_adu", "mad_adu",
        "error"])
    meta = pd.DataFrame([
        {"key": "built_utc", "value": datetime.now(timezone.utc).isoformat()},
        {"key": "code_version", "value": FLIP_CODE_VERSION},
        {"key": "boundary", "value": args.boundary},
        {"key": "epoch_before", "value": before_id},
        {"key": "epoch_after", "value": after_id},
        {"key": "candidate_sigma", "value": str(CANDIDATE_SIGMA)},
        {"key": "persistence", "value": str(PERSISTENCE)},
        {"key": "loose_sigma", "value": str(LOOSE_SIGMA)},
        {"key": "loose_persistence", "value": str(LOOSE_PERSISTENCE)},
        {"key": "n_shifts", "value": str(N_SHIFTS)},
        {"key": "shift_seed", "value": str(SHIFT_SEED)},
    ])
    write_tables(args.manifest, {"mech_flip_test": tests,
                                 "mech_flip_frames": frames,
                                 "mech_flip_bias": bias_df,
                                 "mech_flip_meta": meta})
    print(f"[flip] wrote {', '.join(FLIP_TABLES)} -> {args.manifest}")
    if maps:
        print(f"[flip] figure -> "
              f"{draw_figure(tests, maps, args.boundary, args.figure)}")
    for _, b in bias_df.groupby("mech_epoch"):
        ok = b["median_adu"].dropna()
        if len(ok):
            print(f"[flip] bias level {b['mech_epoch'].iloc[0]}: median "
                  f"{np.median(ok):.1f} ADU over {len(ok)} frames "
                  f"(range {ok.min():.1f}–{ok.max():.1f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
