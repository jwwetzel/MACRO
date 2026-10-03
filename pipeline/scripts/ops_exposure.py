#!/usr/bin/env python
"""Exposure times for the QHY600 (a 16-bit camera), derived from pixels.

WHY THIS SCRIPT EXISTS
----------------------
Revision 2 of the observatory request asked for "3 x 1 s in r ... short to
stay below the High-Gain clip" and for 240 s grism frames "same mode, same
exposure" as 2025.  Both sentences describe cameras that are no longer on the
telescope (DE section 2, TE.F7): the High Gain clip at 3.5 kADU belonged to
the SBIG AC4040 (out of the beam since 2024-03), and the 2025 grism series is
ZWO ASI / Mode0.  October's camera is the QHY600M under pyscope, which clips
at 65,535 ADU (``detector_params``, era group "(blank 2026)").

T CrB has never been observed with the QHY (zero frames in eras 78-83), so
its exposures cannot be read off the archive.  They can be DERIVED from it,
and this script does so without a typed number in the chain:

GRISM (hrg, lrg)
    1. *Bridge.*  Stars observed through the same grism at the same
       exposure with BOTH cameras -- the Be programme's standards eta Hya,
       theta Vir, theta Crt, and theta CrB and kappa Dra -- give the ratio
       of QHY to ASI peak ADU per binned pixel, empirically, optics and
       gain setting included.  The QHY era is split into its MaxIm months
       (READOUTM 'Fast', era 78) and its pyscope-native days (blank
       READOUTM, era 81), because the pyscope headers carry
       ``FOCOFFCG = 0`` for the grisms: if the +650-count hrg focus offset
       of 2025 is no longer applied the trace is defocused and its peak
       drops, and that would show here as a second, lower ratio.
    2. *T CrB on the ASI.*  For every gate-accepted T CrB frame in
       ``g_extractions`` the along-trace peak pixel is measured twice:
       the continuum level (90th percentile of the per-column peak on the
       bright part of the trace) and the H-alpha emission peak (maximum
       within +/- ``HALPHA_HALF_PX`` of the fitted line position).
    3. *Prediction.*  QHY peak at 240 s = (2) x (1).  Headroom is quoted
       against the saturation veto in ``detector_params`` and against the
       NATIVE-pixel ceiling: the frames are 2x2 AVERAGE-binned, so a binned
       peak understates the brightest native pixel (DE.F9: 7-16 %; the
       conservative 16 % is applied).
    4. *Eruption ladder.*  The quiescent ADU rate scaled by
       10**(0.4 * (V_quiescent - V)) gives the exposure that puts the peak
       at ``TARGET_PEAK_ADU`` for V = 2 ... 10.  This is an order-of-
       magnitude guide -- the eruption spectrum is not the quiescent one --
       which is why the schedule block brackets it with a x4 ladder.

IMAGING (g, r, i)
    The pyscope reduction writes ``ZMAG`` (zero point for 1 ADU/s; it does
    not depend on exposure time -- checked, see ``zmag_vs_exptime``) and
    ``FWHM`` into every solved QHY frame; the manifest holds them.  The
    median of the as-found configuration (era 82, the reduced twins of the
    2026-06-28 -> 07-02 frames) converts an assumed T CrB magnitude into a
    total rate, and a Gaussian of the stated FWHM into a peak pixel:

        peak = rate * t * 4 ln 2 / (pi * FWHM_px**2)

    ROADMAP convention 8 says ZMAG is QC, never calibration.  Planning an
    exposure is QC.  No photometry is calibrated with this number.

    The table is a function of image FWHM because the request asks for a
    mild DEFOCUS: at 1 s the scintillation noise of a 0.5 m aperture is
    ~1 %, so the exposure must be several seconds, and in focus several
    seconds saturates a 9th-magnitude star.  Scintillation is Young's (1967)
    approximation as given by Dravins et al. (1998):

        sigma = 0.09 D^(-2/3) X^1.75 exp(-h / 8000 m) / sqrt(2 t)    [D in cm]

ASSUMPTIONS THAT ARE INPUTS, NOT RESULTS (all in ``ASSUMED``)
    * T CrB quiescent magnitudes V, g, r, i.  They come from the AAVSO
      record and Sloan transformations and are good to ~0.3 mag; the table
      therefore prints a bright and a faint bracket, and the request asks
      for the first-night frame to be checked against the 15-40 kADU
      acceptance window rather than trusting this prediction.
    * The archive is read-only and is only ever opened for reading.

OUTPUTS (``ops/generated/``)
    ``exposure_measurements.csv``   one row per frame measured
    ``exposure_bridge.md``          QHY/ASI ratio per star and grism
    ``exposure_grism.md``           T CrB grism prediction + headroom
    ``exposure_eruption.md``        eruption exposure ladder
    ``exposure_imaging.md``         g/r/i exposure vs FWHM
    ``calibration_census.md``       QHY exposures / filters / set-points in
                                    use, i.e. which darks and flats to take
    ``eruption_plan.md``            the eruption ladders per regime
and, under ``ops/eruption_block/``, one pyscope-style ``.sch`` per eruption
regime plus ``tcrb_eruption_exposures.csv`` (the authoritative content).
    ``fig_exposure.{png,pdf}``      the bridge and the T CrB peaks

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/ops_exposure.py measure     # ~10 min, reads frames
    $PY pipeline/scripts/ops_exposure.py report      # seconds, from the CSV
    $PY pipeline/scripts/ops_exposure.py report --inject ops/2026-10_observatory_request_rev3.md
"""

from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "pipeline"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from macro_core import plotstyle                     # noqa: E402
from ops_visibility import inject                    # noqa: E402

OPS_EXPOSURE_VERSION = "1.0 (2026-10-03)"

MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
ARCHIVE = Path("/Volumes/OWC StudioStack HDD/DATA/ASTRO/rlmt-archive")
OUT_DIR = REPO_ROOT / "ops" / "generated"
MEASUREMENTS = OUT_DIR / "exposure_measurements.csv"

#: Camera label -> manifest era ids (raw tree).  Era 76 is the ZWO ASI in
#: Mode0; 78 the QHY600 under MaxIm ('Fast'); 81 the QHY600 under pyscope
#: (blank READOUTM) -- the configuration October will find.
CAMERAS = {"ASI Mode0": (76,), "QHY MaxIm": (78,), "QHY pyscope": (81,)}
REFERENCE_CAMERA = "ASI Mode0"

#: Bridge groups: (manifest target, display name, filter, exposure [s]).
#: Chosen because each has >= 8 raw frames with BOTH cameras at one
#: exposure (query in the report).  Spica is excluded: 4 s / 0.5 s frames
#: are scintillation-dominated in the peak pixel.
BRIDGE = (
    ("HR 3454", "η Hya", "hrg", 60.0), ("HR 3454", "η Hya", "lrg", 15.0),
    ("HR 4963", "θ Vir", "hrg", 60.0), ("HR 4963", "θ Vir", "lrg", 15.0),
    ("HR 4468", "θ Crt", "hrg", 60.0), ("HR 4468", "θ Crt", "lrg", 15.0),
    ("tet CrB", "θ CrB", "hrg", 60.0), ("kap Dra", "κ Dra", "hrg", 60.0),
)
#: Frames sampled per (group, camera): evenly spaced in time, so the sample
#: is deterministic and spans the season.
N_PER_GROUP = 8

#: A column belongs to "the bright part of the trace" when its peak exceeds
#: this fraction of the ``TRACE_REF_PERCENTILE``-th percentile column peak.
#: The reference is the 95th percentile, not the maximum or the 99.5th: a
#: zero-order image (10-25 columns of 4,800) or T CrB's H-alpha line (~20
#: columns) would otherwise SET the threshold and push the continuum
#: columns below it.  A first-order trace covers 25-60 % of the columns,
#: so the 95th percentile always lies on it.
TRACE_FRACTION = 0.2
TRACE_REF_PERCENTILE = 95.0
#: Half-width [px] of the window searched for the H-alpha emission peak
#: around the fitted line position.  T CrB's line is 13-20 px wide (OA.E1).
HALPHA_HALF_PX = 25
#: Half-height [rows] of the band searched around the fitted trace.
TRACE_HALF_ROWS = 8
#: Inner and outer distance [rows] from the trace of the two flanking
#: bands whose median is the local sky.  The inner edge clears the wings
#: of a defocused trace (cross-dispersion FWHM up to ~15 px).
SKY_BAND_ROWS = (25, 60)
#: A bridge frame is rejected when its sky background exceeds this multiple
#: of its group's median (twilight / Moon) or when the trace is at the
#: clip (the statistic is then a lower limit, useless for a ratio).
BG_REJECT_FACTOR = 2.0
CLIP_REJECT_ADU = 55_000.0
#: ... or when its trace level is below this fraction of the group median:
#: the star is not in the frame, or is behind cloud (levels of 20-40 ADU
#: occur where the group median is 6,000).
FAINT_REJECT_FRACTION = 0.2
#: Average binning hides the brightest native pixel: DE.F9 measured the
#: binned peak to understate it by 7-16 % for FWHM 6-4 native px.  The
#: conservative end is applied to every headroom number.
NATIVE_PEAK_FACTOR = 1.16
#: QHY600 pedestal [ADU] in the as-found configuration (OFFSET 10): the
#: ``OSCNMEAN`` overscan mean of the reduced headers and the median sky
#: of every short QHY frame measured here agree at 172-174 ADU.
QHY_PEDESTAL_ADU = 174.0
#: Peak level [ADU, binned, above background] the exposures aim for:
#: mid-way in OA's 15-40 kADU acceptance window.
TARGET_PEAK_ADU = 25_000.0
ACCEPT_WINDOW_ADU = (15_000.0, 40_000.0)

#: Inputs that are assumptions, not measurements.  Quiescent T CrB: V from
#: the AAVSO record (strategy section 1: V ~ 10.2 in Aug 2026, ~10.0 at the
#: 2025 epochs); g, r, i from V, B-V ~ 1.4 and the Jester et al. (2005)
#: transformations for a red giant, good to ~0.3 mag.
ASSUMED = {
    "tcrb_V_quiescent": 10.0,
    "tcrb_mag": {"g": 10.7, "r": 9.5, "i": 8.3},
    "mag_bracket": 0.5,
    "aperture_cm": 50.8,          # OBSDIA card: 0.508 m
    "site_height_m": 1515.7,      # OBSELEV card
    "pixscale_arcsec": 0.449,     # SECPIX1 card, 2x2 binned
}
#: The manifest era whose reduced frames carry the as-found ZMAG / FWHM.
ASFOUND_REDUCED_ERA = 82
IMAGING_FILTERS = ("g", "r", "i")
#: FWHM grid [arcsec] of the imaging table: in focus -> deliberate defocus.
FWHM_GRID_ARCSEC = (2.0, 2.5, 4.0, 6.0, 8.0)
ERUPTION_V = (2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
#: Eruption regimes: (file tag, title, V bright, V faint, direct imaging?).
#: Direct imaging is not attempted while the star is brighter than V = 7:
#: it saturates at any exposure the camera can time (strategy risk 5), and
#: the grisms are the attenuator.
ERUPTION_REGIMES = (
    ("1_peak", "peak, V 2-4", 2.0, 4.0, False),
    ("2_early_decline", "early decline, V 4-7", 4.0, 7.0, False),
    ("3_late_decline", "late decline and secondary maximum, V 7-10",
     7.0, 10.0, True),
)
#: Ladder step between rungs, and the bracket applied beyond the regime's
#: bright and faint ends, both as exposure ratios.
LADDER_STEP = 4.0
#: Shortest and longest exposure a rung may take [s].  The floor is where
#: exposure-time accuracy and rolling-shutter skew stop being negligible
#: for an unverified driver (DE on TCRB-P0-shutter-timing); the ceiling is
#: the programme's standard grism exposure.
LADDER_MIN_S, LADDER_MAX_S = 0.05, 240.0
#: Defocused image size [arcsec] assumed for eruption imaging (as section
#: 6.A of the request).
ERUPTION_IMAGING_FWHM = 4.0
ERUPTION_DIR = REPO_ROOT / "ops" / "eruption_block"
#: T CrB, ICRS (SIMBAD), as the schedule block needs it.
TCRB_RA_HMS, TCRB_DEC_DMS = "15:59:30.16", "+25:55:12.6"

#: Raw-tree eras of the QHY600 (MaxIm months + pyscope days).
QHY_RAW_ERAS = (78, 81)
#: Dark exposures are requested for the smallest set of exposure times
#: that covers this fraction of all QHY light frames, plus every exposure
#: the rev. 3 programme itself asks for (``PROGRAMME_EXPOSURES_S``).
DARK_COVERAGE = 0.95
PROGRAMME_EXPOSURES_S = (5.0, 10.0, 15.0, 60.0, 120.0, 240.0)


# ---------------------------------------------------------------------------
# Pixel statistics (pure functions; tested on synthetic frames)
# ---------------------------------------------------------------------------
def despike(data: np.ndarray) -> np.ndarray:
    """3x3 median filter: removes isolated hot pixels (the IMX455's
    16.5 kADU 'rail' pixels, DE.F1).

    The price is a peak biased LOW: for a cross-dispersion profile of
    sigma = 2 px the filtered peak is ~6 % under the true one, and for an
    emission line that is also narrow along the trace ~12 %.  The bias is
    the same on both cameras, so it cancels in the bridge ratio; in the
    absolute T CrB level it is covered by ``NATIVE_PEAK_FACTOR`` being
    taken at the conservative end and by a headroom that is a factor of
    tens, not percent."""
    from scipy.ndimage import median_filter
    return median_filter(data, size=3)


def trace_peak_blind(data: np.ndarray) -> dict:
    """Peak statistics of the brightest spectrum in a frame, no trace model.

    For a bright standard the target trace owns the brightest pixel of
    almost every column it crosses.  Returns the sky background (frame
    median), the 90th percentile of the per-column peak over the bright
    part of the trace (``cont_p90``: the continuum level, insensitive to
    the ~10 columns of a zero-order image and to single emission lines),
    and the absolute maximum (``peak_max``).
    """
    bg = float(np.median(data))
    net = despike(data) - bg
    col = net.max(axis=0)
    on = col > TRACE_FRACTION * np.percentile(col, TRACE_REF_PERCENTILE)
    return {"bg": bg, "cont_p90": float(np.percentile(col[on], 90)),
            "peak_max": float(col.max()), "n_cols": int(on.sum())}


def trace_peak_along(data: np.ndarray, coeffs, x_line: float) -> dict:
    """Peak statistics ALONG a fitted trace (``np.polyval`` coefficients),
    above the LOCAL sky.

    Used for T CrB, which is faint enough at 240 s that field stars and
    their zero orders own the brightest pixels elsewhere in the frame, and
    whose sky is not flat: a slitless frame's background is the dispersed
    sky, a sharp-edged lozenge several hundred ADU above the frame median
    in a 240 s lrg exposure (OA.E4, E7).  A frame-median background would
    count that lozenge as starlight.  The sky is therefore taken per
    column as the median of two flanking bands ``SKY_BAND_ROWS`` either
    side of the trace.

    Returns, in ADU above local sky: ``cont_p90`` (90th percentile of the
    per-column peak over the bright part of the trace), ``line_peak`` (the
    maximum within ``HALPHA_HALF_PX`` of ``x_line``), ``peak_max`` (the
    brightest on-trace pixel anywhere -- the emission line if ``x_line``
    is right, and still the pixel that saturates first if it is not, which
    matters because OA.E1 found the lrg line positions in
    ``g_extractions`` unreliable); ``sky_local`` (the
    median flanking sky above the frame median, for the record); and
    ``contrast``, the continuum level over the same statistic measured in
    a flanking band where there is no star.  A polynomial that does not
    belong to this frame gives contrast ~ 1 and the caller drops the frame.
    """
    ny, nx = data.shape
    bg = float(np.median(data))
    net = despike(data) - bg
    lo_b, hi_b = SKY_BAND_ROWS
    h = TRACE_HALF_ROWS
    yc = np.polyval(coeffs, np.arange(nx))
    col = np.full(nx, np.nan)       # on-trace peak above local sky
    off = np.full(nx, np.nan)       # same statistic in an empty band
    sky = np.full(nx, np.nan)
    for xi in range(nx):
        y0 = int(round(yc[xi]))
        if y0 - hi_b < 0 or y0 + hi_b + 1 > ny:
            continue
        below = net[y0 - hi_b:y0 - lo_b, xi]
        above = net[y0 + lo_b + 1:y0 + hi_b + 1, xi]
        sky[xi] = np.median(np.concatenate([below, above]))
        col[xi] = net[y0 - h:y0 + h + 1, xi].max() - sky[xi]
        # An empty window of the same height, centred in the upper band.
        yb = y0 + (lo_b + hi_b) // 2
        off[xi] = net[yb - h:yb + h + 1, xi].max() - sky[xi]
    good = np.isfinite(col)
    on = good & (col > TRACE_FRACTION * np.nanpercentile(
        col, TRACE_REF_PERCENTILE))
    cont = float(np.percentile(col[on], 90))
    lo = max(int(x_line) - HALPHA_HALF_PX, 0)
    hi = min(int(x_line) + HALPHA_HALF_PX + 1, nx)
    line = float(np.nanmax(col[lo:hi])) if np.isfinite(col[lo:hi]).any() \
        else float("nan")
    off_level = float(np.percentile(off[on], 90))
    return {"bg": bg, "cont_p90": cont, "line_peak": line,
            "peak_max": float(col[on].max()),
            "sky_local": float(np.median(sky[on])),
            "contrast": cont / max(off_level, 1.0), "n_cols": int(on.sum())}


def gaussian_peak_fraction(fwhm_px: float) -> float:
    """Fraction of a 2-D Gaussian's flux in its central pixel (peak
    surface brightness x 1 px^2): 4 ln 2 / (pi FWHM^2)."""
    return 4.0 * np.log(2.0) / (np.pi * fwhm_px ** 2)


def rate_from_zmag(mag: float, zmag: float) -> float:
    """Total count rate [ADU/s] of a star of magnitude ``mag`` where
    ``zmag`` is the magnitude giving 1 ADU/s."""
    return 10.0 ** (0.4 * (zmag - mag))


def scintillation_mmag(exptime_s: float, airmass: float,
                       aperture_cm: float, height_m: float) -> float:
    """Young's scintillation approximation, in mmag (1.0857 x fractional)."""
    frac = (0.09 * aperture_cm ** (-2.0 / 3.0) * airmass ** 1.75
            * np.exp(-height_m / 8000.0) / np.sqrt(2.0 * exptime_s))
    return float(1085.7 * frac)


def exposure_for_peak(peak_rate_adu_s: float,
                      target_adu: float = TARGET_PEAK_ADU) -> float:
    """Exposure [s] that brings a peak pixel accumulating at
    ``peak_rate_adu_s`` to ``target_adu``."""
    return target_adu / peak_rate_adu_s


def robust(values) -> tuple:
    """``(median, upper quartile, 1.4826 x MAD, n)`` of a sample."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan"), float("nan"), 0
    med = float(np.median(v))
    return (med, float(np.percentile(v, 75)),
            float(1.4826 * np.median(np.abs(v - med))), int(v.size))


# ---------------------------------------------------------------------------
# Measurement (reads the archive, read-only)
# ---------------------------------------------------------------------------
def load_frame(path: Path) -> np.ndarray:
    """First 2-D image HDU of a (possibly fpack-ed) frame, as float32."""
    import warnings
    from astropy.io import fits
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # malformed CONTINUE cards
        with fits.open(path, mode="readonly") as hdul:
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim == 2:
                    return hdu.data.astype(np.float32)
    raise ValueError(f"no 2-D image in {path}")


def evenly(rows: list, n: int) -> list:
    """``n`` rows evenly spaced through a time-ordered list (all if fewer)."""
    if len(rows) <= n:
        return rows
    idx = np.linspace(0, len(rows) - 1, n).round().astype(int)
    return [rows[i] for i in idx]


def measure(manifest: Path = MANIFEST, archive: Path = ARCHIVE,
            out: Path = MEASUREMENTS, n_per_group: int = N_PER_GROUP,
            kinds: tuple = ("bridge", "tcrb")) -> int:
    """Measure every bridge frame and every accepted T CrB frame; write
    one CSV row per frame.  Returns the number of rows written.

    ``kinds`` restricts the run to one family of frames; rows of the other
    family already in the CSV are carried over unchanged, so re-measuring
    the 56 T CrB frames does not cost a re-read of the 142 bridge frames.
    """
    con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    jobs = []
    for key, name, filt, exptime in BRIDGE:
        for camera, eras in CAMERAS.items():
            marks = ",".join("?" * len(eras))
            rows = con.execute(
                f"SELECT path, night FROM frames WHERE tree='rawimage' "
                f"AND canonical_target=? AND filter=? "
                f"AND ABS(exptime-?) < 0.01 AND era_id IN ({marks}) "
                f"ORDER BY jd", (key, filt, exptime, *eras)).fetchall()
            for path, night in evenly(rows, n_per_group):
                jobs.append(dict(kind="bridge", target=name, filter=filt,
                                 exptime=exptime, camera=camera, night=night,
                                 path=path, coeffs=None, x_line=None))
    # T CrB: every gate-accepted frame with a fitted trace and line.
    for path, filt, night, exptime, c0, c1, c2, xh in con.execute(
            "SELECT path, filter, night, exptime, trace_c0, trace_c1, "
            "trace_c2, x_halpha FROM g_extractions WHERE target='T CrB' "
            "AND method='flanking' AND gate_verdict='ACCEPT' "
            "AND trace_c2 IS NOT NULL AND x_halpha IS NOT NULL "
            "ORDER BY jd"):
        jobs.append(dict(kind="tcrb", target="T CrB", filter=filt,
                         exptime=exptime, camera=REFERENCE_CAMERA,
                         night=night, path=path, coeffs=(c0, c1, c2),
                         x_line=xh))
    con.close()

    jobs = [j for j in jobs if j["kind"] in kinds]
    carried = []
    if out.exists() and set(kinds) != {"bridge", "tcrb"}:
        with out.open() as fh:
            carried = [r for r in csv.DictReader(fh)
                       if r["kind"] not in kinds]

    out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["kind", "target", "filter", "exptime", "camera", "night",
              "path", "bg", "cont_p90", "peak_max", "line_peak", "sky_local",
              "contrast", "n_cols", "error"]
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        for r in carried:
            w.writerow({k: r.get(k, "") for k in fields})
        for i, job in enumerate(jobs, 1):
            row = {k: job[k] for k in ("kind", "target", "filter", "exptime",
                                       "camera", "night", "path")}
            try:
                data = load_frame(archive / job["path"])
                if job["kind"] == "bridge":
                    row.update(trace_peak_blind(data))
                else:
                    row.update(trace_peak_along(data, job["coeffs"],
                                                job["x_line"]))
            except Exception as exc:            # recorded, never hidden
                row["error"] = f"{type(exc).__name__}: {exc}"
            w.writerow(row)
            fh.flush()
            print(f"[{i}/{len(jobs)}] {job['camera']:<12} {job['target']:<7}"
                  f" {job['filter']} {job['night']} "
                  f"{row.get('cont_p90', row.get('error'))}", flush=True)
    return len(jobs) + len(carried)


# ---------------------------------------------------------------------------
# Analysis of the measurements
# ---------------------------------------------------------------------------
def read_measurements(path: Path = MEASUREMENTS) -> list:
    """The CSV as a list of dicts with numeric fields converted."""
    rows = []
    with path.open() as fh:
        for r in csv.DictReader(fh):
            for k in ("exptime", "bg", "cont_p90", "peak_max", "line_peak",
                      "sky_local", "contrast", "n_cols"):
                r[k] = float(r[k]) if r[k] not in ("", None) else float("nan")
            rows.append(r)
    return rows


def bridge_table(rows: list) -> list:
    """Per (star, grism): robust level per camera and the QHY/ASI ratios.

    The level adopted per camera is the UPPER QUARTILE of the per-frame
    continuum peak: cloud and poor seeing only ever lower a peak, and the
    question being asked is how bright the trace gets.
    """
    out = []
    for _key, name, filt, exptime in BRIDGE:
        rec = {"target": name, "filter": filt, "exptime": exptime}
        for camera in CAMERAS:
            grp = [r for r in rows if r["kind"] == "bridge"
                   and r["target"] == name and r["filter"] == filt
                   and r["camera"] == camera and not r["error"]]
            bg_med = np.median([r["bg"] for r in grp]) if grp else np.nan
            lvl_med = np.median([r["cont_p90"] for r in grp]) if grp \
                else np.nan
            kept = [r["cont_p90"] for r in grp
                    if r["bg"] <= BG_REJECT_FACTOR * bg_med
                    and r["peak_max"] < CLIP_REJECT_ADU
                    and r["cont_p90"] >= FAINT_REJECT_FRACTION * lvl_med]
            med, q75, sig, n = robust(kept)
            rec[camera] = {"median": med, "q75": q75, "sigma": sig, "n": n,
                           "n_rejected": len(grp) - n}
        ref = rec[REFERENCE_CAMERA]["q75"]
        for camera in CAMERAS:
            rec[camera]["ratio"] = rec[camera]["q75"] / ref
        out.append(rec)
    return out


def adopted_ratio(bridge: list, filt: str, camera: str) -> tuple:
    """Median over stars of the QHY/ASI ratio for one grism and camera,
    with the full range over stars and the number of stars."""
    r = [b[camera]["ratio"] for b in bridge
         if b["filter"] == filt and b[camera]["n"] >= 3
         and np.isfinite(b[camera]["ratio"])]
    if not r:
        return float("nan"), float("nan"), float("nan"), 0
    return float(np.median(r)), float(min(r)), float(max(r)), len(r)


def tcrb_levels(rows: list, filt: str) -> dict:
    """T CrB on the ASI at its archival exposure: continuum level and
    brightest on-trace pixel [ADU above local sky], and the local sky
    itself, over the frames with a believable trace (contrast >= 3).
    ``line`` / ``line_max`` are the brightest on-trace pixel (median over
    frames / brightest frame): H-alpha in hrg; in lrg H-alpha or the
    zero-order image, whichever is brighter -- the conservative quantity
    for a saturation estimate either way."""
    grp = [r for r in rows if r["kind"] == "tcrb" and r["filter"] == filt
           and not r["error"] and r["contrast"] >= 3.0]
    n_all = sum(1 for r in rows if r["kind"] == "tcrb"
                and r["filter"] == filt)
    return {"cont": robust([r["cont_p90"] for r in grp]),
            "line": robust([r["peak_max"] for r in grp]),
            "line_max": max((r["peak_max"] for r in grp), default=np.nan),
            "sky": robust([r["sky_local"] for r in grp]),
            "exptime": grp[0]["exptime"] if grp else float("nan"),
            "n": len(grp), "n_all": n_all}


def detector_limits(manifest: Path = MANIFEST) -> dict:
    """Ceiling and saturation veto of the as-found QHY configuration, from
    ``detector_params`` (S2), not typed."""
    con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    try:
        d = dict(con.execute(
            "SELECT quantity, value FROM detector_params "
            "WHERE era_group='(blank 2026)'").fetchall())
    finally:
        con.close()
    return {"ceiling": d["ceiling_adu"], "veto": d["saturation_veto_adu"]}


def asfound_zeropoints(manifest: Path = MANIFEST) -> dict:
    """Median ZMAG and FWHM per filter in the as-found reduced era, plus
    the check that ZMAG does not depend on exposure time (the slope of
    ZMAG against 2.5 log10 t over ALL QHY reduced frames: 0 for a per-
    second zero point, 1 for a per-exposure one)."""
    con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    out = {}
    try:
        for filt in IMAGING_FILTERS:
            rows = con.execute(
                "SELECT zmag, fwhm FROM frames WHERE era_id=? AND filter=? "
                "AND zmag IS NOT NULL AND fwhm IS NOT NULL",
                (ASFOUND_REDUCED_ERA, filt)).fetchall()
            z = robust([r[0] for r in rows])
            f = robust([r[1] for r in rows])
            allq = np.array(con.execute(
                "SELECT zmag, exptime FROM frames WHERE era_id IN (79, 82) "
                "AND filter=? AND zmag IS NOT NULL AND exptime > 0",
                (filt,)).fetchall())
            slope = float(np.polyfit(2.5 * np.log10(allq[:, 1]),
                                     allq[:, 0], 1)[0])
            out[filt] = {"zmag": z[0], "zmag_sigma": z[2], "fwhm": f[0],
                         "fwhm_q25": float(np.percentile(
                             [r[1] for r in rows], 25)),
                         "n": z[3], "zmag_vs_exptime": slope,
                         "n_slope": int(len(allq))}
    finally:
        con.close()
    return out


def round_sig(x: float, sig: int = 2) -> float:
    """Round to ``sig`` significant figures (exposures are not typed to
    the millisecond)."""
    if x <= 0:
        return 0.0
    return float(f"{x:.{sig}g}")


def ladder(t_bright: float, t_faint: float, step: float = LADDER_STEP,
           t_min: float = LADDER_MIN_S, t_max: float = LADDER_MAX_S) -> list:
    """Geometric exposure ladder covering a brightness regime.

    ``t_bright`` and ``t_faint`` are the exposures that put the peak at
    the target level at the bright and faint ends of the regime.  The
    ladder starts one ``step`` SHORTER than ``t_bright`` and ends one
    ``step`` LONGER than ``t_faint`` -- the prediction is an order-of-
    magnitude one and the bracket is what makes at least one rung usable
    -- clipped to ``[t_min, t_max]`` and rounded to two figures.
    """
    lo = max(t_bright / step, t_min)
    hi = min(t_faint * step, t_max)
    rungs, t = [], lo
    while t < hi * (1 - 1e-9):
        rungs.append(round_sig(t))
        t *= step
    rungs.append(round_sig(hi))
    out = []
    for r in rungs:                       # drop duplicates after rounding
        if not out or r > out[-1]:
            out.append(r)
    return out


def eruption_rates(rows: list, bridge: list, zp: dict) -> dict:
    """Quiescent peak-pixel rate [ADU/s, brightest native pixel] per
    filter on the QHY: grisms from the bridge, imaging from ZMAG at the
    defocused FWHM.  Everything in the eruption ladders scales from this.
    """
    rate = {}
    for filt in ("hrg", "lrg"):
        t = tcrb_levels(rows, filt)
        ratio = adopted_ratio(bridge, filt, "QHY MaxIm")[0]
        rate[filt] = (max(t["line_max"], t["cont"][1]) * ratio
                      * NATIVE_PEAK_FACTOR / t["exptime"])
    for filt in IMAGING_FILTERS:
        rate[filt] = (rate_from_zmag(ASSUMED["tcrb_mag"][filt],
                                     zp[filt]["zmag"])
                      * gaussian_peak_fraction(
                          ERUPTION_IMAGING_FWHM / ASSUMED["pixscale_arcsec"])
                      * NATIVE_PEAK_FACTOR)
    return rate


def eruption_plan(rate: dict) -> list:
    """One dict per regime: its ladders per filter."""
    vq = ASSUMED["tcrb_V_quiescent"]
    plan = []
    for tag, title, v_bright, v_faint, imaging in ERUPTION_REGIMES:
        filters = ("lrg", "hrg") + (IMAGING_FILTERS if imaging else ())
        ladders = {}
        for filt in filters:
            t_b = exposure_for_peak(rate[filt] * 10 ** (0.4 * (vq - v_bright)))
            t_f = exposure_for_peak(rate[filt] * 10 ** (0.4 * (vq - v_faint)))
            ladders[filt] = ladder(t_b, t_f)
        plan.append({"tag": tag, "title": title, "v_bright": v_bright,
                     "v_faint": v_faint, "ladders": ladders})
    return plan


def write_eruption_block(plan: list, out_dir: Path = ERUPTION_DIR) -> list:
    """Write one pyscope-style ``.sch`` per regime and one CSV of every
    exposure.  Returns the paths written.

    The CSV is the authoritative statement of content.  The ``.sch``
    syntax follows the block fields pyscope writes into its own headers
    (``BLKNAME, BLKRA, BLKDEC, BLKPRI, BLKFILT, BLKEXP, BLKNEXP``); pyscope
    is not installed where this script runs, so the files have NOT been
    parsed by ``pyscope.telrun`` -- that is the dry run the request asks
    the site for.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    csv_path = out_dir / "tcrb_eruption_exposures.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["regime", "v_bright", "v_faint", "order", "target",
                    "filter", "exposure_s", "nexp"])
        for reg in plan:
            order = 0
            for filt, rungs in reg["ladders"].items():
                for t in rungs:
                    order += 1
                    w.writerow([reg["tag"], reg["v_bright"], reg["v_faint"],
                                order, "T CrB", filt, f"{t:g}", 1])
    written.append(csv_path)
    for reg in plan:
        lines = [
            "# DRAFT -- NOT LOADED, NOT VALIDATED BY pyscope.telrun.",
            "# T CrB eruption response: " + reg["title"],
            f"# Generated by pipeline/scripts/ops_exposure.py "
            f"v{OPS_EXPOSURE_VERSION}; do not edit by hand.",
            "# Exposures are a x4 ladder so that at least one rung is",
            "# well exposed whatever the star is doing; see README.md.",
            "# Site: fill in observer and code; check keyword spelling",
            "# against your pyscope version; then dry-run and send the log.",
            "",
            'title "T CrB eruption - ' + reg["title"] + '"',
            'observer "<observer e-mail>"',
            'code "<observing code>"',
            "",
            "block start",
            "    priority 1",
            "    do_not_interrupt true",
            '    source "T CrB"',
            f"    ra {TCRB_RA_HMS}",
            f"    dec {TCRB_DEC_DMS}",
            "    binning 2x2",
            "    repositioning true",
        ]
        for filt, rungs in reg["ladders"].items():
            for t in rungs:
                lines.append(f"    filter {filt} exposure {t:g} nexp 1")
        lines += ["block end", ""]
        path = out_dir / f"tcrb_eruption_{reg['tag']}.sch"
        path.write_text("\n".join(lines))
        written.append(path)
    return written


def eruption_plan_md(plan: list, rate: dict) -> str:
    """The ladders, as a table for the README."""
    lines = ["| Regime | lrg (s) | hrg (s) | g (s) | r (s) | i (s) |",
             "|---|---|---|---|---|---|"]
    for reg in plan:
        cells = [", ".join(f"{t:g}" for t in reg["ladders"][f])
                 if f in reg["ladders"] else "not attempted"
                 for f in ("lrg", "hrg") + IMAGING_FILTERS]
        lines.append(f"| {reg['title']} | " + " | ".join(cells) + " |")
    rates = ", ".join(f"{f} {rate[f]:.1f}" if rate[f] < 100
                      else f"{f} {rate[f]:,.0f}" for f in rate)
    lines += ["",
              f"Each ladder runs from ×{LADDER_STEP:g} shorter than the "
              f"exposure predicted for the bright end of the regime to "
              f"×{LADDER_STEP:g} longer than that for the faint end, in "
              f"steps of ×{LADDER_STEP:g}, clipped to "
              f"{LADDER_MIN_S:g}–{LADDER_MAX_S:g} s; target peak "
              f"{TARGET_PEAK_ADU:,.0f} ADU. Quiescent peak-pixel rates used "
              f"(ADU s⁻¹ at V = {ASSUMED['tcrb_V_quiescent']:g}): {rates}. "
              f"Imaging assumes the telescope defocused to "
              f"{ERUPTION_IMAGING_FWHM:g}″ FWHM.", "", _stamp()]
    return "\n".join(lines)


def calibration_census(manifest: Path = MANIFEST) -> dict:
    """What the QHY era needs calibrating, counted from the manifest.

    ``exposures``: light-frame count per exposure time (raw tree, QHY
    eras), most-used first, with the cumulative fraction and whether the
    exposure is in the requested dark set.  ``filters``: light frames per
    filter and the number of raw flats the archive holds for it
    (``calib_gaps.have_raw``).  ``temps``: frames per CCD temperature,
    rounded to the degree.  ``n_calib``: raw bias/dark/flat frames of any
    kind in these eras.
    """
    marks = ",".join("?" * len(QHY_RAW_ERAS))
    con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    try:
        exp = con.execute(
            f"SELECT ROUND(exptime, 3), COUNT(*) FROM frames WHERE "
            f"tree='rawimage' AND era_id IN ({marks}) AND exptime > 0 "
            f"GROUP BY 1 ORDER BY 2 DESC", QHY_RAW_ERAS).fetchall()
        filt = con.execute(
            f"SELECT filter, COUNT(*) FROM frames WHERE tree='rawimage' "
            f"AND era_id IN ({marks}) GROUP BY 1 ORDER BY 2 DESC",
            QHY_RAW_ERAS).fetchall()
        have = {}
        for spec, n in con.execute(
                f"SELECT spec, have_raw FROM calib_gaps WHERE "
                f"need_kind='flat' AND era_id IN ({marks})", QHY_RAW_ERAS):
            name = spec.split()[1]
            have[name] = have.get(name, 0) + (n or 0)
        temps = con.execute(
            f"SELECT ROUND(ccd_temp), COUNT(*) FROM frames WHERE "
            f"tree='rawimage' AND era_id IN ({marks}) GROUP BY 1 "
            f"ORDER BY 2 DESC", QHY_RAW_ERAS).fetchall()
        n_calib = con.execute(
            f"SELECT COUNT(*) FROM calib_frames WHERE era_id IN ({marks})",
            QHY_RAW_ERAS).fetchone()[0]
        span = con.execute(
            f"SELECT MIN(night), MAX(night), COUNT(*) FROM frames WHERE "
            f"tree='rawimage' AND era_id IN ({marks})",
            QHY_RAW_ERAS).fetchone()
    finally:
        con.close()
    total = sum(n for _, n in exp)
    rows, cum = [], 0
    for t, n in exp:
        need = cum / total < DARK_COVERAGE or t in PROGRAMME_EXPOSURES_S
        cum += n
        rows.append({"exptime": t, "n": n, "cum": cum / total, "dark": need})
    return {"exposures": rows, "filters": [(f, n, have.get(f, 0))
                                           for f, n in filt],
            "temps": temps, "n_calib": n_calib, "span": span,
            "total": total}


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def _stamp() -> str:
    return (f"*Generated by `pipeline/scripts/ops_exposure.py` "
            f"v{OPS_EXPOSURE_VERSION} from `ops/generated/"
            f"exposure_measurements.csv` and the manifest. "
            f"Do not edit by hand.*")


def _k(adu: float) -> str:
    """ADU with thousands separator, or an em dash."""
    return "—" if not np.isfinite(adu) else f"{adu:,.0f}"


def bridge_md(bridge: list) -> str:
    cams = [c for c in CAMERAS if c != REFERENCE_CAMERA]
    head = ("| Star | Grism | Exp (s) | " + REFERENCE_CAMERA
            + " peak (ADU), n | "
            + " | ".join(f"{c} peak (ADU), n | {c} / ASI" for c in cams)
            + " |")
    lines = [head, "|---|---|---:|---:|" + "---:|---:|" * len(cams)]
    for b in bridge:
        ref = b[REFERENCE_CAMERA]
        cells = [f"{_k(ref['q75'])} ± {_k(ref['sigma'])}, {ref['n']}"]
        for c in cams:
            m = b[c]
            if m["n"] == 0:
                cells += ["—", "—"]
            else:
                cells += [f"{_k(m['q75'])} ± {_k(m['sigma'])}, {m['n']}",
                          f"{m['ratio']:.2f}"]
        lines.append(f"| {b['target']} | {b['filter']} | {b['exptime']:g} | "
                     + " | ".join(cells) + " |")
    for filt in ("hrg", "lrg"):
        for c in cams:
            med, lo, hi, n = adopted_ratio(bridge, filt, c)
            if n:
                lines.append(f"| **adopted** | **{filt}** | | | **{c}: "
                             f"{med:.2f}** (range {lo:.2f}–{hi:.2f}, "
                             f"{n} stars) | |" + " |" * (2 * len(cams) - 2))
    lines += ["",
              "Peak = upper quartile over frames of the 90th-percentile "
              "per-column peak pixel along the trace, ADU above sky, 2×2 "
              "average-binned pixels; ± is 1.4826 × MAD over frames; n = "
              "frames kept (twilight, clipped and no-star frames rejected). Same "
              "star, same grism, same exposure on both cameras.", "",
              _stamp()]
    return "\n".join(lines)


def grism_md(rows: list, bridge: list, limits: dict) -> str:
    """T CrB 240 s on the QHY: predicted peak and headroom."""
    lines = ["| Grism | T CrB on ASI, frames used | Continuum peak (ADU) | "
             "Brightest on-trace pixel, median / max (ADU) | Sky under "
             "trace (ADU) | QHY / ASI | Predicted on QHY, median / max "
             "(ADU) | Brightest raw native pixel (ADU) | Fraction of "
             f"{limits['veto']:,.0f} ADU veto |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for filt in ("hrg", "lrg"):
        t = tcrb_levels(rows, filt)
        ratio, lo, hi, n = adopted_ratio(bridge, filt, "QHY MaxIm")
        med, mx = t["line"][0] * ratio, t["line_max"] * ratio
        # Raw value of the brightest native pixel: star (x native-peak
        # factor) + sky + pedestal, the last two carried at the ASI level
        # (the QHY pedestal is lower, so this is conservative).
        raw = (mx * NATIVE_PEAK_FACTOR + t["sky"][1] * ratio
               + QHY_PEDESTAL_ADU)
        lines.append(
            f"| {filt} {t['exptime']:g} s | {t['n']} of {t['n_all']} | "
            f"{_k(t['cont'][0])} ± {_k(t['cont'][2])} | "
            f"{_k(t['line'][0])} / {_k(t['line_max'])} | "
            f"{_k(t['sky'][0])} | {ratio:.2f} | "
            f"{_k(med)} / {_k(mx)} | {_k(raw)} | "
            f"{raw / limits['veto']:.1%} |")
    lines += ["",
              "T CrB frames are the gate-accepted rows of `g_extractions` "
              "(method `flanking`) whose fitted trace stands at least 3× "
              "above the same statistic in an empty flanking band. Levels "
              "are peak pixels above the local sky (median of two bands "
              f"{SKY_BAND_ROWS[0]}–{SKY_BAND_ROWS[1]} rows either side of "
              "the trace). The brightest on-trace pixel is the Hα line in "
              "hrg; in lrg it is the line or the zero-order image. The "
              "prediction uses the QHY/ASI ratio measured under MaxIm (the "
              "larger, hence conservative, of the two QHY ratios); the "
              "raw-pixel column adds the sky and the "
              f"{QHY_PEDESTAL_ADU:g} ADU pedestal and applies the ×"
              f"{NATIVE_PEAK_FACTOR:g} native-pixel factor to the brightest "
              "archival frame.", "", _stamp()]
    return "\n".join(lines)


def eruption_md(rows: list, bridge: list, limits: dict) -> str:
    """Exposure that puts the trace peak at TARGET_PEAK_ADU, vs V."""
    vq = ASSUMED["tcrb_V_quiescent"]
    rate = {}
    for filt in ("hrg", "lrg"):
        t = tcrb_levels(rows, filt)
        ratio = adopted_ratio(bridge, filt, "QHY MaxIm")[0]
        # Brightest of line and continuum, brightest archival frame: the
        # pixel that saturates first.
        rate[filt] = (max(t["line_max"], t["cont"][1]) * ratio
                      * NATIVE_PEAK_FACTOR / t["exptime"])
    lines = ["| V | Brightening (×) | hrg exposure (s) | lrg exposure (s) |",
             "|---:|---:|---:|---:|"]
    for v in ERUPTION_V:
        boost = 10.0 ** (0.4 * (vq - v))
        cells = []
        for filt in ("hrg", "lrg"):
            t_exp = round_sig(exposure_for_peak(rate[filt] * boost))
            cells.append(f"{t_exp:g}" if t_exp < 100 else f"{t_exp:,.0f}")
        lines.append(f"| {v:g} | {boost:,.0f} | {cells[0]} | {cells[1]} |")
    lines += ["",
              f"Exposure that brings the brightest native trace pixel to "
              f"{TARGET_PEAK_ADU:,.0f} ADU (ceiling {limits['ceiling']:,.0f}"
              f"), scaling the measured quiescent rate (hrg "
              f"{rate['hrg']:.1f}, lrg {rate['lrg']:.1f} ADU s⁻¹ at V = "
              f"{vq:g}, assumed) by the V-band brightening. Order of "
              "magnitude only: the eruption spectrum is hotter than the "
              "quiescent one and Hα will not scale with V. The schedule "
              "block therefore brackets each value by ×4 either side.", "",
              _stamp()]
    return "\n".join(lines)


def imaging_md(zp: dict) -> str:
    """g/r/i exposure for T CrB vs image FWHM, with scintillation."""
    scale = ASSUMED["pixscale_arcsec"]
    br = ASSUMED["mag_bracket"]
    lines = ["| Filter | ZMAG (as found), n | Median FWHM (″) | Assumed "
             "T CrB mag | " + " | ".join(
                 f"t at FWHM {f:g}″ (s)" for f in FWHM_GRID_ARCSEC) + " |",
             "|---|---:|---:|---:|" + "---:|" * len(FWHM_GRID_ARCSEC)]
    for filt in IMAGING_FILTERS:
        z = zp[filt]
        mag = ASSUMED["tcrb_mag"][filt]
        cells = []
        for fwhm in FWHM_GRID_ARCSEC:
            def t_for(m):
                peak = (rate_from_zmag(m, z["zmag"])
                        * gaussian_peak_fraction(fwhm / scale)
                        * NATIVE_PEAK_FACTOR)
                return exposure_for_peak(peak)
            cells.append(f"{t_for(mag):.1f} ({t_for(mag - br):.1f}–"
                         f"{t_for(mag + br):.1f})")
        lines.append(f"| {filt} | {z['zmag']:.2f} ± {z['zmag_sigma']:.2f}, "
                     f"{z['n']} | {z['fwhm']:.1f} | {mag:.1f} ± {br:g} | "
                     + " | ".join(cells) + " |")
    ap, h = ASSUMED["aperture_cm"], ASSUMED["site_height_m"]
    sc = " · ".join(
        f"{t:g} s: {scintillation_mmag(t, 1.2, ap, h):.1f} / "
        f"{scintillation_mmag(t, 1.8, ap, h):.1f}" for t in (1, 5, 10, 20))
    slopes = ", ".join(f"{f} {zp[f]['zmag_vs_exptime']:+.3f} "
                       f"(n = {zp[f]['n_slope']:,})" for f in IMAGING_FILTERS)
    lines += ["",
              f"Exposure that brings the brightest native pixel of T CrB to "
              f"{TARGET_PEAK_ADU:,.0f} ADU for a Gaussian image of the stated "
              f"FWHM; in parentheses the range for the star {br:g} mag "
              f"brighter – fainter than assumed. ZMAG and FWHM are medians "
              f"over the reduced frames of era {ASFOUND_REDUCED_ERA} "
              f"(2026-06-28 → 07-02, the as-found configuration).",
              "",
              f"Scintillation (Young; D = {ap:g} cm, h = {h:g} m), mmag at "
              f"airmass 1.2 / 1.8 — {sc}.",
              "",
              f"ZMAG is a per-second zero point: slope of ZMAG against "
              f"2.5 log₁₀ t over all QHY reduced frames — {slopes} "
              f"(0 = per second, 1 = per exposure).", "", _stamp()]
    return "\n".join(lines)


def calibration_md(census: dict) -> str:
    """The QHY calibration debt, as counts: which darks, which flats."""
    first, last, n_all = census["span"]
    darks = sorted(r["exptime"] for r in census["exposures"] if r["dark"])
    lines = ["| Exposure (s) | Light frames | Cumulative | Dark requested |",
             "|---:|---:|---:|---|"]
    for r in census["exposures"]:
        if not r["dark"] and r["n"] < 100:
            continue
        lines.append(f"| {r['exptime']:g} | {r['n']:,} | {r['cum']:.1%} | "
                     f"{'yes' if r['dark'] else 'no (scaled)'} |")
    lines += ["",
              f"Dark set requested ({len(darks)} exposures): "
              + ", ".join(f"{t:g}" for t in darks) + " s.", "",
              "| Filter | Light frames | Raw flats in archive |",
              "|---|---:|---:|"]
    for f, n, have in census["filters"]:
        if n >= 10:
            lines.append(f"| {f} | {n:,} | {have} |")
    temps = "; ".join(
        f"{'no card' if t is None else format(t, '+.0f') + ' °C'}: {n:,}"
        for t, n in census["temps"])
    lines += ["",
              f"QHY600 raw frames in the archive: {n_all:,} "
              f"({first} → {last}, manifest eras "
              f"{', '.join(map(str, QHY_RAW_ERAS))}). Raw bias, dark or "
              f"flat frames among them: {census['n_calib']}. "
              f"CCD-TEMP census — {temps}. Darks are requested at every "
              f"exposure needed to cover {DARK_COVERAGE:.0%} of the light "
              "frames plus the exposures this request itself asks for; "
              "rarer exposures are scaled from the nearest longer dark.",
              "", _stamp()]
    return "\n".join(lines)


def figure(path_stem: Path, rows: list, bridge: list) -> None:
    """Left: the camera bridge, star by star.  Right: T CrB peak pixels on
    the ASI against the 16-bit ceiling -- the picture of the headroom."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with plotstyle.context("web"):
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.6, 4.2),
                                       constrained_layout=True)
        cams = [c for c in CAMERAS if c != REFERENCE_CAMERA]
        labels = [f"{b['target']}\n{b['filter']}" for b in bridge]
        for i, c in enumerate(cams):
            y = [b[c]["ratio"] if b[c]["n"] >= 3 else np.nan for b in bridge]
            ax0.plot(np.arange(len(bridge)) + 0.12 * (i - 0.5), y,
                     label=c, **plotstyle.measurement_kw(
                         **plotstyle.series(i)))
        ax0.axhline(1.0, **plotstyle.reference_kw())
        ax0.set_xticks(range(len(bridge)))
        ax0.set_xticklabels(labels, fontsize="small")
        ax0.set_ylabel("trace peak, QHY600 / ASI (same star, exposure)")
        ax0.set_ylim(0, None)
        ax0.set_title("Camera bridge from stars both cameras observed")
        ax0.legend(frameon=False)

        for i, filt in enumerate(("hrg", "lrg")):
            grp = [r for r in rows if r["kind"] == "tcrb"
                   and r["filter"] == filt and not r["error"]
                   and r["contrast"] >= 3.0]
            st = plotstyle.series(i)
            ax1.plot([r["cont_p90"] for r in grp],
                     [r["peak_max"] for r in grp],
                     label=f"{filt}, 240 s (n = {len(grp)})",
                     **plotstyle.measurement_kw(**st))
        ax1.axhline(65_535, **plotstyle.reference_kw(color=plotstyle.BAD))
        ax1.text(0.02, 0.93, "16-bit ceiling, 65,535 ADU",
                 transform=ax1.transAxes, color=plotstyle.BAD,
                 fontsize="small")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_xlabel("continuum peak pixel (ADU above local sky)")
        ax1.set_ylabel("brightest on-trace pixel (ADU above local sky)")
        ax1.set_title("T CrB, 2025 grism frames (ASI): far below the clip")
        ax1.legend(frameon=False, loc="lower right")
        fig.savefig(path_stem.with_suffix(".png"), dpi=plotstyle.WEB_DPI)
        fig.savefig(path_stem.with_suffix(".pdf"))
        plt.close(fig)


def report(out_dir: Path = OUT_DIR, manifest: Path = MANIFEST) -> dict:
    """Render every fragment and the figure from the measurements CSV."""
    rows = read_measurements(out_dir / MEASUREMENTS.name)
    bridge = bridge_table(rows)
    limits = detector_limits(manifest)
    zp = asfound_zeropoints(manifest)
    fragments = {
        "exposure_bridge": bridge_md(bridge),
        "exposure_grism": grism_md(rows, bridge, limits),
        "exposure_eruption": eruption_md(rows, bridge, limits),
        "exposure_imaging": imaging_md(zp),
        "calibration_census": calibration_md(calibration_census(manifest)),
    }
    rate = eruption_rates(rows, bridge, zp)
    plan = eruption_plan(rate)
    write_eruption_block(plan)
    fragments["eruption_plan"] = eruption_plan_md(plan, rate)
    for name, body in fragments.items():
        (out_dir / f"{name}.md").write_text(body + "\n")
    figure(out_dir / "fig_exposure", rows, bridge)
    return fragments


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=("measure", "report"))
    ap.add_argument("--kind", choices=("bridge", "tcrb"), default=None,
                    help="measure: re-measure only this family of frames")
    ap.add_argument("--inject", type=Path, default=None,
                    help="Markdown file whose GENERATED blocks to rewrite")
    args = ap.parse_args()
    if args.command == "measure":
        n = measure(kinds=(args.kind,) if args.kind
                    else ("bridge", "tcrb"))
        print(f"wrote {n} rows to {MEASUREMENTS}")
        return 0
    fragments = report()
    for body in fragments.values():
        print(body, end="\n\n")
    if args.inject:
        done = inject(args.inject, fragments)
        print(f"injected into {args.inject}: {', '.join(done) or 'nothing'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
