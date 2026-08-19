"""S1 astrometry go/no-go experiment: strata logic + solve-field wrapper.

ROADMAP §1.1 S1 + §5: the CV paper's polar time series (ST LMi, VV Pup,
EU UMa, AN UMa in Sloan g/r/i) are 73–95% UNSOLVED (``pltsolvd`` = 0, no
zero point).  Before anyone batch-solves a ~59k-frame backlog we run a
*stratified experiment*: sample a few dozen frames from each homogeneous
(camera era family × exposure band × project) stratum, solve them with a
local astrometry.net, and report success rates with confidence intervals —
so the Week-2 review can make the batch decision on evidence, not vibes.

Everything in this module that decides something is a *pure function*
(no I/O, no globals mutated) unit-tested in ``pipeline/tests/test_astrom.py``:
stratum classification, deterministic sampling, Wilson intervals, solve
command construction, residual arithmetic, failure diagnosis.  The only
impure code is the thin ``solve_one_frame`` wrapper (funpack + solve-field
subprocesses) that the experiment script drives; it returns a plain dict
and writes nothing to the database itself.

The experiment writes three NEW tables into the manifest database
(existing tables are never modified — S0/S0b remain the source of truth):

* ``s1_strata``           — one row per stratum: definition, population,
                            sample size, the stored RNG seed.
* ``s1_solve_experiment`` — one row per sampled frame: solve outcome,
                            wall time, pixel scale, match count, RMS.
* ``s1_failure_autopsy``  — image-statistics post-mortem for EVERY
                            failure (source counts, trailing, saturation)
                            plus a machine diagnosis.
"""

from __future__ import annotations

import math
import zlib
from dataclasses import dataclass
from random import Random
from typing import Optional, Sequence

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these,
# so changing a value here changes both the pipeline and its documentation).
# --------------------------------------------------------------------------

#: Version string recorded into ``s1_build_meta``.
S1_CODE_VERSION = "S1 v1.1 (2026-08-18)"

#: Master RNG seed for the whole experiment.  Every stratum derives its own
#: child seed from (this, stratum_id) — see ``stable_seed`` — so the sample
#: is reproducible frame-for-frame from the manifest alone.
SAMPLE_SEED = 20260817

#: Frames sampled per stratum (fewer when the population is smaller).
#: 48 successes/failures give a Wilson 95% interval no wider than ±14% at
#: p = 0.5 and ±7% at p = 0.9 — tight enough to call go/no-go.
N_PER_STRATUM = 48

#: Per-frame solve budget.  solve-field gets ``SOLVE_CPU_LIMIT_S`` of CPU;
#: the subprocess is killed at ``SOLVE_TIMEOUT_S`` wall seconds regardless
#: (funpack + I/O overhead live in the gap between the two).
SOLVE_CPU_LIMIT_S = 55
SOLVE_TIMEOUT_S = 75

#: solve-field source-extraction downsampling (mission spec: 2).  Halves
#: the pixel grid before detection — big win on 4–5k-pixel frames.
SOLVE_DOWNSAMPLE = 2

#: When the header carries a pointing (frames.ra_deg/dec_deg parsed from
#: OBJCTRA/OBJCTDEC), restrict the index search to this radius.  Wide
#: enough to survive the worst mount pointing error in the S0 audit
#: (pointing_offset_deg is typically < 1 deg); narrow enough to cut the
#: search space by ~99%.
HINT_RADIUS_DEG = 15.0

#: FILTER strings that mean "this is a slitless spectrum, not an image".
#: A grism frame has no point sources for quad matching — it is excluded
#: from the astrometry candidate universe, not counted as a failure.
GRISM_FILTERS = frozenset({"hrg", "lrg", "hagrism", "oggrism", "grism"})

#: FILTER strings that collide with the calibration vocabulary (the era-76
#: filter-wheel glitch documented by S0b) — excluded the same way.
CALIB_VOCAB_FILTERS = frozenset({"bias", "dark", "flat"})

#: A frame narrower than this on either axis is a high-speed *photometry
#: window* (the archive holds 8×3211-pixel strips), not a field image:
#: too little sky for quad matching, excluded by geometry.  512 px at the
#: coarsest plate scale here (~1.1"/px) is still a ~9' field — the
#: smallest footprint worth attempting.
MIN_SOLVABLE_NAXIS = 512

#: The CV paper's four polars — the populations S1 exists to unblock.
CV_TARGET_KEYS = frozenset({"stlmi", "vvpup", "euuma", "anuma"})

#: Wilson z for 95% two-sided intervals.
WILSON_Z = 1.959963984540054

#: Solution-acceptance gate (see ``solution_sane``).  A ``.solved`` marker
#: alone is NOT success: astrometry.net can emit a false positive on a
#: sparse quad match (the experiment caught exactly one — frame 149276
#: recorded 3.38"/px from 4 matched stars on a 0.45"/px camera).  A
#: solution counts only when all three hold:
#:
#: * the measured pixel scale sits inside the prior handed to the solver
#:   (a solution outside its own search bounds indicts the match);
#: * at least ``MIN_MATCHED_STARS`` index stars matched — every genuine
#:   solve in the experiment matched >= 8, the false positive matched 4;
#: * the astrometric RMS stays under ``MAX_SOLVE_RMS_ARCSEC`` — genuine
#:   solves topped out at 4.83", the false positive scored 5.63".
MIN_MATCHED_STARS = 8
MAX_SOLVE_RMS_ARCSEC = 5.0

#: Verdict thresholds on the stratum success rate's Wilson LOWER bound —
#: pessimistic by construction: a stratum is GO only when even the
#: unluckiest reading of its sample clears the bar.
GO_LOWER_BOUND = 0.80
CAUTION_LOWER_BOUND = 0.50

#: Parallel solve workers for wall-clock projections (and the runner's
#: default).  Ten concurrent solve-field processes fit the Mac's cores
#: with the ``inparallel`` index config.
DEFAULT_WORKERS = 10


# --------------------------------------------------------------------------
# Candidate gates: which unsolved frames can astrometry even apply to?
# --------------------------------------------------------------------------

def norm_filter(filt) -> str:
    """Normalize a FILTER string for gate checks: lowercase, trimmed,
    NULL→''.  ('g', ' G ', None → 'g', 'g', '')."""
    return (filt or "").strip().lower()


def is_grism_filter(filt) -> bool:
    """True when the FILTER string names a dispersing element — the frame
    is a spectrum and can never be plate-solved."""
    return norm_filter(filt) in GRISM_FILTERS


def is_calib_vocab_filter(filt) -> bool:
    """True for the filter-wheel glitch strings ('dark' science frames)."""
    return norm_filter(filt) in CALIB_VOCAB_FILTERS


def is_window_geometry(naxis1, naxis2) -> bool:
    """True when the frame is a sub-frame photometry window (either axis
    below ``MIN_SOLVABLE_NAXIS``) — geometrically unsolvable.

    NULL axes count as windows too: a frame whose geometry is unknown
    cannot be promised a solvable field.
    """
    if naxis1 is None or naxis2 is None:
        return True
    return naxis1 < MIN_SOLVABLE_NAXIS or naxis2 < MIN_SOLVABLE_NAXIS


def is_solvable_candidate(row: dict) -> bool:
    """Full candidate gate over one frames-row dict (keys: filter, naxis1,
    naxis2).  The SQL base query already restricts to canonical rawimage
    Light frames with pltsolvd != 1; this function applies the pure gates.
    """
    return (not is_grism_filter(row.get("filter"))
            and not is_calib_vocab_filter(row.get("filter"))
            and not is_window_geometry(row.get("naxis1"), row.get("naxis2")))


# --------------------------------------------------------------------------
# The candidate-universe base query: every gate that lives naturally in SQL.
# ONE definition, shared by the experiment runner (to design the strata) and
# the report (to compute per-stratum night coverage) — two copies of this
# query would eventually disagree about what "the universe" is.
# The *pure* gates (grism filter, window geometry) stay in Python above so
# they remain unit-testable; SQL only narrows to "canonical raw Light frame
# that the headers call unsolved".
# --------------------------------------------------------------------------

BASE_SQL = """
SELECT f.obs_rowid, f.path, f.target_key, f.canonical_target, e.readoutm,
       f.xbinning, f.filter, f.exptime, f.naxis1, f.naxis2,
       f.ra_deg, f.dec_deg, f.night
FROM frames f LEFT JOIN eras e ON e.era_id = f.era_id
WHERE f.is_canonical = 1
  AND f.tree = 'rawimage'
  AND f.error IS NULL
  AND (f.imagetyp LIKE 'Light%' OR f.imagetyp IS NULL OR f.imagetyp = '')
  AND (f.pltsolvd IS NULL OR f.pltsolvd != 1)
"""

#: Column order of BASE_SQL, used to build row dicts.
BASE_COLS = ["obs_rowid", "path", "target_key", "canonical_target",
             "readoutm", "xbinning", "filter", "exptime", "naxis1",
             "naxis2", "ra_deg", "dec_deg", "night"]


def fetch_candidates(con) -> list[dict]:
    """All unsolved canonical raw Light frames, as row dicts (thin SQL
    wrapper — the decisions about these rows live in the pure gates)."""
    return [dict(zip(BASE_COLS, r)) for r in con.execute(BASE_SQL)]


# --------------------------------------------------------------------------
# Strata: one label per homogeneous (camera family × band × project) cell
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Stratum:
    """One experiment stratum: id, verdict population it feeds, and the
    human definition the report prints (the *code* definition is
    ``classify_stratum`` — the report quotes both)."""
    stratum_id: str
    population: str          # verdict grouping: which go/no-go it feeds
    description: str


#: The design, in priority order (CV-critical first, per the roadmap).
STRATA: tuple[Stratum, ...] = (
    Stratum("cv_mode0_sloan_short", "CV polars (Sloan)",
            "ST LMi / VV Pup / EU UMa / AN UMa — Mode0 (ASI/IMX455) bin2, "
            "Sloan g/r/i, exptime < 150 s"),
    Stratum("cv_mode0_sloan_long", "CV polars (Sloan)",
            "same targets/camera, Sloan g/r/i, exptime ≥ 150 s"),
    Stratum("cv_ikon_sloan", "CV polars (Sloan)",
            "VV Pup on the Andor iKon (1MHz 16-bit, 2048²), Sloan + blank "
            "filter"),
    Stratum("cv_gsense_misc", "CV polars (Sloan)",
            "CV polars on the GSENSE4040 (High Gain family) — small "
            "legacy set, mostly ST LMi"),
    Stratum("cv_fast_fullframe", "CV polars (Sloan)",
            "CV polars on the fast-readout camera, bin2 full frames — the "
            "EU UMa season, invisible until the 8x3211 metadata artifact "
            "was repaired (see docs/pipeline/s0e_geometry_fix.html)"),
    Stratum("sn_gsense_broadband", "SN 2023ixf",
            "SN 2023ixf on the GSENSE4040 (High Gain / StackPro bin1), "
            "all filters"),
    Stratum("dwarf_gsense_deep", "Dwarf/AGN",
            "Dw survey fields + NGC 5548 + NGC 5238 on the GSENSE4040, "
            "deep broadband/narrowband"),
    Stratum("mode0_backlog_short", "Facility backlog",
            "Mode0 bin2 non-CV imaging, exptime < 10 s (rho Oph campaign "
            "and friends)"),
    Stratum("mode0_backlog_long", "Facility backlog",
            "Mode0 bin2 non-CV imaging, exptime ≥ 10 s"),
    Stratum("fast_fullframe", "Facility backlog",
            "Fast-readout FULL frames (naxis ≥ 512; the 8-px photometry "
            "strips are excluded by geometry)"),
    Stratum("ikon_backlog", "Facility backlog",
            "Andor iKon imaging beyond VV Pup"),
)

#: Readout-mode strings of the Andor iKon CCD (era-47 family).
IKON_READOUTS = frozenset({"1mhz high sensitivity 16-bit",
                           "5mhz high sensitivity 16-bit"})

#: Readout-mode strings of the GSENSE4040 CMOS.
GSENSE_READOUTS = frozenset({"high gain", "low gain", "hdr",
                             "high gain stackpro", "low gain stackpro"})


def classify_stratum(row: dict) -> Optional[str]:
    """Assign one candidate frame to a stratum id (or None: not sampled).

    ``row`` keys used: target_key, canonical_target, readoutm, xbinning,
    filter, exptime, naxis1, naxis2.  The frame must already have passed
    ``is_solvable_candidate`` — classification assumes a solvable field.

    The rules mirror ``STRATA`` one-to-one, first match wins; they are
    deliberately explicit (no clever generality) so a student can check
    each stratum against its description by reading twenty lines.
    """
    if not is_solvable_candidate(row):
        return None
    tkey = (row.get("target_key") or "").strip().lower()
    canon = (row.get("canonical_target") or "").strip().lower()
    readout = (row.get("readoutm") or "").strip().lower()
    xbin = row.get("xbinning")
    filt = norm_filter(row.get("filter"))
    exptime = row.get("exptime")

    # --- CV polars ------------------------------------------------------
    if tkey in CV_TARGET_KEYS:
        if readout == "mode0" and xbin == 2 and filt in ("g", "r", "i"):
            # Two exposure bands: the 30–120 s survey cadence vs the
            # 240 s deep cadence — depth changes both star counts and
            # trailing risk, so they get separate verdicts.
            if exptime is not None and exptime < 150:
                return "cv_mode0_sloan_short"
            return "cv_mode0_sloan_long"
        if readout in IKON_READOUTS:
            # VV Pup's iKon season: Sloan filters plus the blank-FILTER
            # frames of the same series ('empty' string in the headers).
            if filt in ("g", "r", "i", "empty", ""):
                return "cv_ikon_sloan"
        if readout in GSENSE_READOUTS:
            return "cv_gsense_misc"
        if readout == "fast":
            # EU UMa's whole season.  This branch is NEW, and it exists
            # because of a metadata artifact rather than a design change.
            #
            # These 207 frames carried a phantom 8x3211 geometry in the
            # catalog (a tile-compressed BINTABLE's row length read as an
            # image width), so every one of them failed the geometry gate
            # and never reached any stratum at all.  The CV project recorded
            # them as "permanently unsolvable".  They are nothing of the
            # kind: the frames are ordinary 4800x3211 full frames, and a
            # spot check solves them in ~3 s at 74 matched stars and ~1.4"
            # RMS.
            #
            # Note what is NOT done here.  The fall-through below would have
            # swept them into 'fast_fullframe', a FACILITY BACKLOG stratum,
            # which would have quietly changed the population behind an id
            # that is already published.  A new id costs nothing and keeps
            # the rule this project already applies to era numbers: retire
            # or add, never redefine a number someone may have cited.
            return "cv_fast_fullframe"
        return None                     # CV frame in an unplanned config
    # --- SN 2023ixf -----------------------------------------------------
    if tkey == "2023ixf" and readout in GSENSE_READOUTS:
        return "sn_gsense_broadband"
    # --- Dwarf/AGN survey ----------------------------------------------
    if readout in GSENSE_READOUTS and (
            canon.startswith("dw") or tkey in ("ngc5548", "ngc5238")):
        return "dwarf_gsense_deep"
    # --- Facility backlog probes ---------------------------------------
    if readout == "mode0" and xbin == 2:
        if exptime is not None and exptime < 10:
            return "mode0_backlog_short"
        return "mode0_backlog_long"
    if readout == "fast":
        # The comment that used to sit here said "geometry gate already
        # dropped strips".  There were never any strips: the 8x3211 "sub-
        # frame photometry windows" were a tile-compressed BINTABLE's row
        # length misread as an image width, and repairing the catalog took
        # the geometry exclusion on this camera to ZERO.  This stratum
        # absorbed the repaired full frames and grew about sevenfold as a
        # result.  See macro_core/fitsgeom.py and
        # docs/pipeline/s0e_geometry_fix.html.
        return "fast_fullframe"
    if readout in IKON_READOUTS:
        return "ikon_backlog"
    return None                         # small residue: not worth a stratum


# --------------------------------------------------------------------------
# Deterministic sampling
# --------------------------------------------------------------------------

def stable_seed(base_seed: int, stratum_id: str) -> int:
    """Child seed for one stratum: crc32 of the id folded into the base.

    Python's ``hash`` is salted per process, so it cannot appear anywhere
    near a reproducibility contract; crc32 is stable across processes,
    platforms and Python versions.
    """
    return (base_seed * 2654435761 + zlib.crc32(stratum_id.encode())) \
        % (2 ** 31)


def sample_frames(rowids: Sequence[int], n: int, base_seed: int,
                  stratum_id: str) -> list[int]:
    """Reproducible sample of up to ``n`` rowids for one stratum.

    The candidate list is sorted first (query order is not a contract),
    then shuffled by a Random seeded from (base_seed, stratum_id) — so
    the sample depends only on the candidate SET and the stored seed.
    """
    ordered = sorted(rowids)
    Random(stable_seed(base_seed, stratum_id)).shuffle(ordered)
    return ordered[:n]


# --------------------------------------------------------------------------
# Success-rate arithmetic
# --------------------------------------------------------------------------

def wilson_ci(k: int, n: int, z: float = WILSON_Z) -> tuple[float, float]:
    """Wilson score 95% interval for k successes in n trials.

    Chosen over the normal approximation because strata sit near the
    boundaries (rates of ~0 or ~1 with n ≈ 48) where the naive interval
    collapses to zero width or leaks outside [0, 1].
    """
    if n == 0:
        return (0.0, 1.0)               # no data: total ignorance
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, center - half), min(1.0, center + half))


def verdict_for(k: int, n: int, n_population: Optional[int] = None) -> str:
    """GO / CAUTION / NO-GO from the Wilson lower bound (pessimistic:
    judge each stratum by the worst rate its sample still allows).

    CENSUS special case: when the sample IS the population
    (``n == n_population``), there is no sampling uncertainty to be
    pessimistic about — the population rate is known exactly, so the
    thresholds judge the exact rate, not a Wilson bound.  (Without this,
    a fully-enumerated 25/40 stratum would read NO-GO off an interval
    that models uncertainty which does not exist.)
    """
    if n_population is not None and n > 0 and n == n_population:
        rate = k / n
        if rate >= GO_LOWER_BOUND:
            return "GO"
        if rate >= CAUTION_LOWER_BOUND:
            return "CAUTION"
        return "NO-GO"
    lo, _ = wilson_ci(k, n)
    if lo >= GO_LOWER_BOUND:
        return "GO"
    if lo >= CAUTION_LOWER_BOUND:
        return "CAUTION"
    return "NO-GO"


def night_collapse(night_solved: Sequence[tuple]) -> tuple[int, int]:
    """Collapse per-frame outcomes to per-NIGHT all-or-nothing outcomes.

    Input: (night, solved_bool) pairs, one per sampled frame.  Returns
    (nights with every frame solved, nights sampled).  This is the
    stress test for the frame-level Wilson intervals: frames within a
    night share cloud/focus/wind state, so they are not independent
    trials — the number of NIGHTS is closer to the true effective sample
    size, and a night counts as a success only when nothing on it failed
    (the most conservative collapse).
    """
    fails_by_night: dict = {}
    for night, solved in night_solved:
        # A night is "perfect" until any of its frames fails.
        fails_by_night[night] = fails_by_night.get(night, False) \
            or (not solved)
    n_nights = len(fails_by_night)
    k_nights = sum(1 for failed in fails_by_night.values() if not failed)
    return (k_nights, n_nights)


def projected_hours(n_frames: int, median_solve_s: float,
                    workers: int = DEFAULT_WORKERS) -> float:
    """Wall-clock projection for a batch: n × median ÷ workers, in hours.

    The median (not the mean) is the honest scale here: the time
    distribution is a spike of fast successes plus a wall of
    timeout-capped failures, and the cap makes the mean a function of our
    own timeout setting rather than of the data.
    """
    if n_frames <= 0 or workers <= 0:
        return 0.0
    return n_frames * median_solve_s / workers / 3600.0


# --------------------------------------------------------------------------
# Plate-scale priors per camera family
# --------------------------------------------------------------------------

def scale_bounds(readoutm, xbinning) -> tuple[float, float]:
    """(arcsec/px lower, upper) bounds handed to solve-field.

    Derived from measured hardware, then widened ~±30%:

    * ASI/IMX455 (Mode0/Fast/blank readout): 3.76 µm pixels at 3454 mm
      focal length → 0.225"/px unbinned, 0.449"/px bin2 (the pilot solve
      measured 0.4508 — the prior is honest).
    * GSENSE4040 (High/Low Gain, HDR): 9 µm at ~0.54"/px unbinned.
    * Andor iKon (1MHz/5MHz 16-bit): pixel pitch unverified on this
      telescope — wide bounds, let the index decide.
    """
    readout = (readoutm or "").strip().lower()
    xbin = xbinning or 1
    if readout in GSENSE_READOUTS:
        lo, hi = 0.40, 0.75             # around 0.54"/px unbinned
    elif readout in IKON_READOUTS:
        # Unverified camera: stay wide.  The experiment measured
        # 0.809"/px unbinned (13.5 µm pixels at 3454 mm) — inside these
        # bounds; every sampled iKon frame was bin1, but the binning
        # multiplier below still applies so a hypothetical bin2 frame
        # gets an honest prior instead of falling outside it.
        lo, hi = 0.40, 1.30
    else:                               # Mode0 / Fast / blank → IMX455
        lo, hi = 0.15, 0.33             # around 0.225"/px unbinned
    # Binning scales the pixel on-sky linearly.
    return (lo * xbin, hi * xbin)


def solution_sane(pixscale_arcsec, n_matched, rms_arcsec,
                  scale_lo: float, scale_hi: float) -> bool:
    """Acceptance gate for a solve-field solution (see the constant block
    around ``MIN_MATCHED_STARS`` for the rationale and the caught case).

    A ``.solved`` marker whose WCS fails ANY check is a false positive:
    the experiment records it as ``bad_solve``, never as success, and a
    batch runner must never write its WCS into the manifest.  Missing
    values fail closed — a solution that cannot show its matched-star
    residuals has not earned acceptance.
    """
    if pixscale_arcsec is None or not \
            (scale_lo <= pixscale_arcsec <= scale_hi):
        return False                    # scale outside the solver's own prior
    if n_matched is None or n_matched < MIN_MATCHED_STARS:
        return False                    # too few index stars to trust a quad
    if rms_arcsec is None or rms_arcsec > MAX_SOLVE_RMS_ARCSEC:
        return False                    # residuals worse than any real solve
    return True


# --------------------------------------------------------------------------
# solve-field command construction (pure: a list of argv strings)
# --------------------------------------------------------------------------

def build_solve_command(fits_path: str, config_path: str, out_dir: str,
                        scale_lo: float, scale_hi: float,
                        ra_deg: Optional[float] = None,
                        dec_deg: Optional[float] = None,
                        downsample: int = SOLVE_DOWNSAMPLE,
                        cpulimit_s: int = SOLVE_CPU_LIMIT_S) -> list[str]:
    """argv for one solve-field run.

    Output policy: we keep only the three files the experiment reads —
    ``.solved`` (existence = success), ``.wcs`` (solution header) and
    ``.corr`` (matched-star table for RMS) — and suppress every other
    product (no plots, no new FITS, no xyls) to keep the scratch dir and
    the I/O footprint minimal at 10 workers.
    """
    cmd = [
        "solve-field", fits_path,
        "--config", config_path,
        "--dir", out_dir,               # all products land in scratch
        "--scale-units", "arcsecperpix",
        "--scale-low", f"{scale_lo:g}",
        "--scale-high", f"{scale_hi:g}",
        "--downsample", str(downsample),
        "--cpulimit", str(cpulimit_s),
        "--overwrite", "--no-plots",
        "--new-fits", "none",           # suppress the solved-FITS copy
        "--index-xyls", "none", "--match", "none", "--rdls", "none",
        "--crpix-center",               # report the center pixel's RA/Dec
    ]
    if ra_deg is not None and dec_deg is not None:
        # Header pointing available: restrict the sky search.  The radius
        # covers the worst S0 pointing offsets with a wide margin.
        cmd += ["--ra", f"{ra_deg:.5f}", "--dec", f"{dec_deg:.5f}",
                "--radius", f"{HINT_RADIUS_DEG:g}"]
    return cmd


# --------------------------------------------------------------------------
# Solution readout: residuals from the .corr matched-star table
# --------------------------------------------------------------------------

def sky_residuals_arcsec(field_ra, field_dec, index_ra, index_dec
                         ) -> list[float]:
    """Per-star residuals (arcsec) between solved field positions and the
    index catalog, small-angle spherical: ΔRA·cos(Dec) folded with ΔDec.

    Inputs are parallel sequences of degrees (the .corr table columns).
    """
    out = []
    for fr, fd, ir, id_ in zip(field_ra, field_dec, index_ra, index_dec):
        # RA wraps at 360°; residuals are sub-arcminute, so take the
        # short way around before scaling by cos(Dec).
        dra = ((fr - ir + 180.0) % 360.0) - 180.0
        dra *= math.cos(math.radians(fd))
        ddec = fd - id_
        out.append(math.hypot(dra, ddec) * 3600.0)
    return out


def rms(values: Sequence[float]) -> Optional[float]:
    """Root mean square, None for an empty sequence (no matched stars)."""
    if not values:
        return None
    return math.sqrt(sum(v * v for v in values) / len(values))


# --------------------------------------------------------------------------
# Failure autopsy: image statistics → machine diagnosis
# --------------------------------------------------------------------------

#: Diagnosis thresholds — one place, quoted by the report.
#: Extraction is deliberately CONSERVATIVE (10σ, ≥5 connected pixels): a
#: 3σ cut on the GSENSE CMOS counts tens of thousands of hot-pixel/pattern
#: peaks as "sources" and diagnosed blank frames as healthy fields (the
#: first autopsy round proved it); real stars survive 10σ easily on any
#: frame a solver could have used.
AUTOPSY_EXTRACT_SIGMA = 10.0   # detection threshold, in background σ
AUTOPSY_MIN_AREA_PX = 5        # min connected pixels — kills lone hot px
AUTOPSY_MIN_SOURCES = 15       # fewer PSF-shaped detections = starved
AUTOPSY_DEFOCUS_A_PX = 6.0     # bright-source size above = defocused
AUTOPSY_TRAIL_ELONG = 1.8      # median elongation above = trailing
AUTOPSY_SAT_FRAC = 0.02        # >2% saturated pixels = flooded frame

#: A detection is "PSF-shaped" (a plausible star) when its semi-major axis
#: sits between these bounds, in pixels.  Below: the single-pixel spikes
#: of the RAW frames' hot-pixel forest (no dark subtraction has happened —
#: the second autopsy round measured median a ≈ 0.7 px over thousands of
#: 10σ detections on visually BLANK frames).  Above: blobs, not stars.
PSF_A_MIN_PX = 1.2
PSF_A_MAX_PX = 8.0

#: Blank-frame veto on the BRIGHTEST detections: when the median
#: semi-major axis of the brightest N detections sits below this, the
#: brightest things on the frame are single-pixel spikes — no real star
#: survived, whatever the PSF-band counter says (adjacent hot pixels pair
#: into 1–2 px clumps that sneak past ``PSF_A_MIN_PX`` and fake both
#: "enough sources" and "elongated sources"; the fourth autopsy round
#: caught four visually blank 512 s frames diagnosed as trailing this
#: way).  Calibration from the autopsy table: the four blank frames all
#: measure bright_median_a ≈ 0.71 px; the faintest GENUINELY trailed
#: frame (68275, curved star trails on visual check) measures 1.12 px —
#: the threshold sits between the two clusters.
AUTOPSY_BLANK_BRIGHT_A_PX = 0.9

#: Bright-source size statistic: the median semi-major axis of this many
#: brightest detections — the defocus tell (defocused frames scatter huge
#: donuts among the hot-pixel spikes; the spikes hide them from a global
#: median, the brightest-N median finds them).
N_BRIGHTEST_FOR_SIZE = 20


def diagnose_failure(n_sources: Optional[int],
                     n_psf_sources: Optional[int],
                     median_elongation: Optional[float],
                     saturated_fraction: Optional[float],
                     bright_median_a_px: Optional[float] = None) -> str:
    """Name the most probable failure cause from image statistics.

    Ordered by evidential strength: a frame with almost no PSF-SHAPED
    sources failed for lack of stars whatever else its statistics say
    (hot-pixel spikes do not count — raw frames carry thousands); within
    the starless frames, giant bright blobs mean defocus, otherwise
    blank/clouds.  Once real stars exist: trailing, then saturation, then
    'unexplained' — a genuine solver shortfall worth human eyes, never
    silently binned.
    """
    if n_sources is None:
        return "unreadable"
    # Defocus first among the shape verdicts: a defocused frame fragments
    # its donuts into many elongated pieces, so it fakes both "enough PSF
    # sources" and "trailing" — but nothing except defocus makes the
    # BRIGHTEST detections giant (round-3 autopsy: the filter-'6' blob
    # frames scored median elongation 2.5 with bright-source size 10-27 px).
    if bright_median_a_px is not None and \
            bright_median_a_px > AUTOPSY_DEFOCUS_A_PX:
        return "defocused (giant blobs, no point sources)"
    # Starved has TWO tells, either sufficient: too few PSF-shaped
    # sources, or brightest detections that are sub-star-sized spikes
    # (a blank frame's hot-pixel pairs can inflate the PSF-band count
    # AND its elongation — the brightest-N size is the veto that cannot
    # be faked, because real stars always top the flux ranking when any
    # exist).  Checked before trailing so hot-pixel pairs can never
    # masquerade as wind-trailed stars again.
    if n_psf_sources is None or n_psf_sources < AUTOPSY_MIN_SOURCES \
            or (bright_median_a_px is not None
                and bright_median_a_px < AUTOPSY_BLANK_BRIGHT_A_PX):
        return "starved (blank/clouds/hot-pixel-only)"
    if median_elongation is not None and \
            median_elongation > AUTOPSY_TRAIL_ELONG:
        return "trailing (guiding/wind)"
    if saturated_fraction is not None and \
            saturated_fraction > AUTOPSY_SAT_FRAC:
        return "saturated/flooded"
    return "unexplained (stars present, solver still failed)"


def image_metrics(data, saturation_adu: float) -> dict:
    """Source-extraction statistics for one frame.

    Impure only in its dependency (sep); given the same pixel array it is
    deterministic, and the tests drive it with synthetic star fields.
    Returns: n_sources (all 10σ detections), n_psf_sources (the subset
    shaped like stars), median_elongation + median_a_px (over the PSF
    subset), bright_median_a_px (size of the brightest detections — the
    defocus tell), saturated_fraction, bkg_rms.
    """
    import numpy as np
    import sep
    arr = np.ascontiguousarray(data, dtype=np.float32)
    # Fraction of pixels at/above the saturation rail, measured BEFORE
    # background subtraction (saturation is an absolute ADU fact).
    sat_frac = float(np.mean(arr >= saturation_adu))
    bkg = sep.Background(arr)
    sub = arr - bkg.back()
    try:
        # Conservative extraction (see the threshold block above): only
        # unambiguous sources count toward "the field had stars".
        sources = sep.extract(sub, thresh=AUTOPSY_EXTRACT_SIGMA,
                              err=bkg.globalrms,
                              minarea=AUTOPSY_MIN_AREA_PX)
    except Exception:
        # Pathological frames (constant, NaN-ridden) — report as empty.
        sources = None
    if sources is None or len(sources) == 0:
        return {"n_sources": 0 if sources is not None else 0,
                "n_psf_sources": 0, "median_elongation": None,
                "median_a_px": None, "bright_median_a_px": None,
                "saturated_fraction": sat_frac,
                "bkg_rms": float(bkg.globalrms)}
    n = len(sources)
    a = np.asarray(sources["a"], dtype=float)
    b = np.asarray(sources["b"], dtype=float)
    flux = np.asarray(sources["flux"], dtype=float)
    # The PSF-shaped subset: plausible stars, not spikes and not blobs.
    psf = (a >= PSF_A_MIN_PX) & (a <= PSF_A_MAX_PX)
    n_psf = int(psf.sum())
    med_elong = med_a = None
    if n_psf:
        ap, bp = a[psf], b[psf]
        elong = ap[bp > 0] / bp[bp > 0]            # guard b = 0
        # np.median throughout (interpolating on even counts) — the same
        # median definition as median_a_px, so the two statistics never
        # disagree about what "median" means.
        med_elong = float(np.median(elong)) if len(elong) else None
        med_a = float(np.median(ap))
    # Size of the brightest N detections regardless of shape class: a
    # defocused frame's donuts dominate the flux ranking however many
    # hot-pixel spikes surround them.
    brightest = np.argsort(flux)[-N_BRIGHTEST_FOR_SIZE:]
    bright_a = float(np.median(a[brightest]))
    return {
        "n_sources": n,
        "n_psf_sources": n_psf,
        "median_elongation": med_elong,
        "median_a_px": med_a,
        "bright_median_a_px": bright_a,
        "saturated_fraction": sat_frac,
        "bkg_rms": float(bkg.globalrms),
    }


# --------------------------------------------------------------------------
# The impure wrapper: funpack + solve-field for ONE frame
# --------------------------------------------------------------------------

def solve_one_frame(archive_path: str, work_dir: str, config_path: str,
                    readoutm, xbinning,
                    ra_deg=None, dec_deg=None,
                    timeout_s: int = SOLVE_TIMEOUT_S) -> dict:
    """Solve one archive frame; return a result dict (no DB access here).

    Steps, each explicit so a failure names its stage:
    1. funpack the .fz into the private work dir (raw tree is
       fpack-compressed int16; solve-field's reader wants plain FITS),
    2. run solve-field with the family's scale prior and any pointing
       hint, under a hard wall-clock timeout,
    3. read back .solved/.wcs/.corr, compute match count + RMS,
    4. apply the ``solution_sane`` acceptance gate — a .solved marker
       whose WCS fails it is recorded as ``bad_solve``, never success,
    5. leave cleanup to the caller (the work dir is per-frame temp).
    """
    import shutil
    import subprocess
    import time
    from pathlib import Path

    src = Path(archive_path)
    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    result: dict = {"status": "error", "solve_time_s": None,
                    "solved_ra": None, "solved_dec": None,
                    "pixscale_arcsec": None, "rotation_deg": None,
                    "n_matched": None, "rms_arcsec": None,
                    "used_hint": int(ra_deg is not None
                                     and dec_deg is not None),
                    "log_tail": None}
    t0 = time.monotonic()
    try:
        # -- 1. decompress ------------------------------------------------
        if src.suffix == ".fz":
            plain = work / src.stem       # "x.fts.fz" -> "x.fts"
            proc = subprocess.run(
                ["funpack", "-O", str(plain), str(src)],
                capture_output=True, text=True, timeout=60)
            if proc.returncode != 0:
                result["log_tail"] = ("funpack: "
                                      + proc.stderr.strip()[-400:])
                result["solve_time_s"] = time.monotonic() - t0
                return result
        else:
            plain = work / src.name
            shutil.copyfile(src, plain)
        # -- 2. solve -----------------------------------------------------
        lo, hi = scale_bounds(readoutm, xbinning)
        cmd = build_solve_command(str(plain), config_path, str(work),
                                  lo, hi, ra_deg, dec_deg)
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=timeout_s)
        except subprocess.TimeoutExpired:
            result["status"] = "timeout"
            result["solve_time_s"] = time.monotonic() - t0
            return result
        result["solve_time_s"] = time.monotonic() - t0
        # -- 3. read the verdict back ------------------------------------
        stem = plain.name[:-len(plain.suffix)] if plain.suffix else plain.name
        solved_marker = work / f"{stem}.solved"
        if not solved_marker.exists():
            result["status"] = "unsolved"
            # Keep the tail of the solver log: the autopsy's first clue.
            result["log_tail"] = (proc.stdout.strip()[-400:]
                                  or proc.stderr.strip()[-400:])
            return result
        result["status"] = "solved"
        from astropy.io import fits as _fits
        from astropy.wcs import WCS as _WCS
        with _fits.open(work / f"{stem}.wcs") as hdul:
            w = _WCS(hdul[0].header)
            # Center RA/Dec: CRVAL is the center pixel by --crpix-center.
            result["solved_ra"] = float(w.wcs.crval[0])
            result["solved_dec"] = float(w.wcs.crval[1])
            # Pixel scale + rotation from the CD matrix determinant/angle.
            cd = w.pixel_scale_matrix
            result["pixscale_arcsec"] = float(
                3600.0 * math.sqrt(abs(cd[0, 0] * cd[1, 1]
                                       - cd[0, 1] * cd[1, 0])))
            result["rotation_deg"] = float(
                math.degrees(math.atan2(cd[1, 0], cd[0, 0])))
        corr_path = work / f"{stem}.corr"
        if corr_path.exists():
            with _fits.open(corr_path) as hdul:
                tab = hdul[1].data
                res = sky_residuals_arcsec(
                    tab["field_ra"], tab["field_dec"],
                    tab["index_ra"], tab["index_dec"])
                result["n_matched"] = len(res)
                result["rms_arcsec"] = rms(res)
        # -- 4. acceptance gate ------------------------------------------
        # The solver said "solved"; the gate decides whether we BELIEVE
        # it.  A rejected solution keeps its numbers (they are the
        # forensic evidence) but is statused ``bad_solve`` so it never
        # counts as a success and its WCS is never trusted downstream.
        if not solution_sane(result["pixscale_arcsec"],
                             result["n_matched"], result["rms_arcsec"],
                             lo, hi):
            result["status"] = "bad_solve"
            result["log_tail"] = (
                f"gate: pixscale {result['pixscale_arcsec']} vs "
                f"[{lo:.3f},{hi:.3f}], n_matched {result['n_matched']} "
                f"(min {MIN_MATCHED_STARS}), rms {result['rms_arcsec']} "
                f"(max {MAX_SOLVE_RMS_ARCSEC})")
        return result
    except Exception as exc:            # noqa: BLE001 — recorded, not hidden
        result["log_tail"] = f"{type(exc).__name__}: {exc}"[:400]
        if result["solve_time_s"] is None:
            result["solve_time_s"] = time.monotonic() - t0
        return result
