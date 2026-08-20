"""The external survey record: parsing it, and reading an accretion state off it.

WHY THIS MODULE EXISTS
----------------------
Two CV plan tasks need the public record of what these five stars were doing
while RLMT was pointed at them:

* ``CV-P0-aavso-yzcnc`` — the GATING one.  YZ Cnc is an SU UMa dwarf nova.
  Superhumps exist only in SUPEROUTBURST.  The strategy's Q3 therefore has
  two branches — a superhump period (and possibly dP_sh/dt), or an
  orbital-hump + flickering fallback — and which branch the paper takes is
  not a modelling choice, it is a fact about the accretion state on the
  nights we observed.  Nothing in our own pixels can settle it: our
  photometry is DIFFERENTIAL, so a light curve that is 1.4 mag brighter than
  last night looks exactly like a light curve that is not, until someone
  supplies an absolute reference.  The external record is that reference.

* ``CV-P0-survey-context`` — the long-baseline record for all five targets,
  against which our nights are a handful of dots.

Everything here is PURE: parsing, night arithmetic, aggregation and
classification, with no I/O and no network.  The fetching, caching and
reporting live in ``pipeline/scripts/run_cv_external.py``, so every judgement
this module makes can be tested on synthetic input whose right answer is
known by construction.

THE TWO TRAPS THIS MODULE IS BUILT AROUND
-----------------------------------------
**Trap 1 — the night label is not the UTC date.**  ``cv_frames.night`` is the
LOCAL evening date at Winer Observatory (UTC-7).  Every CV frame in the 2024
YZ Cnc season was taken after local midnight, so the UTC date is the night
label PLUS ONE DAY.  The strategy quotes the dense blocks as "Feb 21-24,
Mar 1-4, May 2-3" — those are UTC dates, and they are the SAME nights the
manifest calls 2024-02-20..23, 2024-02-29/03-02/03-03 and 2024-05-01/02.
Aligning an AAVSO record (Julian Date, absolute) against a night label
(local, ambiguous) off by one day would move every classification onto its
neighbour, and on a dwarf nova that rises 1.4 mag in 24 h that is the whole
answer.  So this module never compares date STRINGS: it converts everything
to JD and derives the UTC night label from the JD.

**Trap 2 — our own data came back to us wearing a stranger's coat.**  A large
share of the AAVSO rows for YZ Cnc in this window were submitted by observer
``MALW`` from THIS TELESCOPE, carrying the comment "TAKEN WITH MACRO
CONSORTIUM'S ROBERT L. MUTEL TELESCOPE AT WINER OBSERVATORY".  Classifying
our own nights from our own resubmitted photometry and calling the result an
independent confirmation would be circular — the exact failure the task
"pull the EXTERNAL record" exists to avoid.  :func:`is_own_observation` marks
those rows, every aggregate is computed twice (all rows, and independent-only
rows), and the classifier records WHICH of the two carried each verdict.

WHAT A STATE CLASSIFICATION RESTS ON
------------------------------------
Not a hand-typed magnitude.  :func:`quiescent_baseline` measures the star's
own faint level from the record itself (the faint-decile median), so the
ladder is anchored to the data in front of it rather than to a number
transcribed from a catalogue.  A night is then classified by its AMPLITUDE
above that baseline, and outburst EPISODES — contiguous runs of outburst
nights — are separated into normal and super by the two properties that
actually distinguish them in an SU UMa star: peak amplitude and duration.
That is the physics, not a threshold someone liked.
"""

from __future__ import annotations

import datetime as dt
import math
import re
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

# ===========================================================================
# 0.  Constants — every one of them a decision, named and defended
# ===========================================================================

#: Julian Date of MJD zero.  JD = MJD + this.
JD_MINUS_MJD = 2400000.5

#: Winer Observatory's UTC offset in hours (Arizona, no daylight saving).
#: Used ONLY to explain the night-label convention in prose and to derive
#: the local evening date from a UTC timestamp; no science number depends
#: on it, because all arithmetic is done in JD.
WINER_UTC_OFFSET_H = -7.0

#: Substrings that identify an AAVSO row as OUR OWN photometry resubmitted.
#: Matched case-insensitively against the observation's comment field.  Both
#: the consortium name and the telescope name are listed because a submitter
#: who trims the comment is likelier to keep one than to keep both.
OWN_COMMENT_MARKERS: tuple[str, ...] = ("MUTEL", "MACRO CONSORTIUM", "WINER")

#: AAVSO observer codes known to submit RLMT data.  Belt AND braces with the
#: comment markers above: a code can be re-used by a different program and a
#: comment can be dropped, so a row is ours if EITHER test fires.
OWN_OBSERVER_CODES: frozenset[str] = frozenset({"MALW"})

#: Bands that are a usable stand-in for Johnson V when classifying a state.
#: 'Vis.' is a visual estimate — coarse (~0.2-0.3 mag) but unbiased in the
#: mean, and on a 4-magnitude outburst amplitude that is ample.  'CV' is an
#: unfiltered CCD magnitude reduced to a V sequence, which is what the AAVSO
#: label means.  Sloan bands are NOT here: they are handled separately,
#: through a measured offset, because assuming SR == V is exactly the sort of
#: unexamined step that makes a wrong answer look tidy.
V_LIKE_BANDS: frozenset[str] = frozenset({"V", "Vis.", "CV", "Vis"})

#: The faint fraction of a season's points whose median defines "quiescence".
#: A dwarf nova spends most of its time near quiescence but not all of it, so
#: the plain median would sit in the outburst tail; the faint decile is deep
#: enough to be quiescence and wide enough not to be one noisy point.
QUIESCENT_DECILE = 0.10

#: Amplitude ladder above the quiescent baseline, in magnitudes.
#: Below LO: quiescence.  Above HI: outburst.  Between: 'elevated' — an
#: honest third answer for a star caught rising or decaying, which a
#: two-way cut would have to lie about.
AMP_QUIESCENT_MAX = 0.6
AMP_OUTBURST_MIN = 1.0

#: What separates a SUPEROUTBURST from a normal outburst in an SU UMa star.
#: Both criteria must hold.  Amplitude alone is not enough (a well-caught
#: normal outburst peak can be bright); duration alone is not enough (a run
#: of unrelated normal outbursts can look long through sparse sampling).
#: YZ Cnc's normal outbursts reach ~2-3 mag above quiescence and last 2-4
#: days; its superoutbursts reach ~4 mag and hold a plateau for 10-14 days.
SUPEROUTBURST_AMP_MIN = 3.0
SUPEROUTBURST_DAYS_MIN = 8.0

#: How far below the peak still counts as "on the plateau", in magnitudes.
#: The duration that grades an episode is the PLATEAU length, not the length
#: of the whole above-threshold excursion.  On a star like YZ Cnc, which is
#: rarely at true quiescence, a plain above-threshold run merges the
#: superoutburst with its own slow decline AND with the normal outbursts
#: that follow it -- the first build of this module graded one 2024 event as
#: a 38-day superoutburst, which is not a thing that exists.  Measuring the
#: time spent within 1.5 mag of peak recovers the physical quantity: the
#: same event's plateau is 12 days, which is what a superoutburst is.
PLATEAU_DEPTH_MAG = 1.5

#: A gap larger than this (days) ends an outburst episode.  Sparse amateur
#: coverage routinely skips a night; three days of silence is long enough
#: that joining across it would be invention rather than interpolation.
EPISODE_MAX_GAP_D = 3.0

STATE_QUIESCENT = "QUIESCENT"
STATE_ELEVATED = "ELEVATED"
STATE_OUTBURST = "OUTBURST"
STATE_UNKNOWN = "NO DATA"

EPISODE_SUPEROUTBURST = "SUPEROUTBURST"
EPISODE_NORMAL = "NORMAL OUTBURST"


# ===========================================================================
# 1.  Time — the one place an off-by-one day would cost the whole answer
# ===========================================================================

def jd_to_datetime(jd: float) -> dt.datetime:
    """Julian Date -> timezone-aware UTC datetime.

    Written out rather than delegated to astropy so the test suite can
    exercise it without importing an ephemeris: this is calendar
    arithmetic, not barycentric arithmetic.  (The BJD_TDB that the science
    actually uses is S3's job and is not computed here.)
    """
    return (dt.datetime(1858, 11, 17, tzinfo=dt.timezone.utc)
            + dt.timedelta(days=jd - JD_MINUS_MJD))


def jd_to_mjd(jd: float) -> float:
    """Julian Date -> Modified Julian Date."""
    return jd - JD_MINUS_MJD


def mjd_to_jd(mjd: float) -> float:
    """Modified Julian Date -> Julian Date."""
    return mjd + JD_MINUS_MJD


def utc_night(jd: float) -> str:
    """The UTC calendar date of an observation, as ``YYYY-MM-DD``.

    This is the label every external source is bucketed by, because it is
    the only one that means the same thing to AAVSO, to ZTF and to us.
    """
    return jd_to_datetime(jd).strftime("%Y-%m-%d")


def local_night_label(jd: float,
                      utc_offset_h: float = WINER_UTC_OFFSET_H) -> str:
    """The LOCAL evening date — the convention ``cv_frames.night`` uses.

    Subtracting the site's UTC offset walks an after-midnight UTC timestamp
    back into the evening it belongs to.  This function exists so the
    manifest's night labels can be mapped onto UTC dates BY COMPUTATION
    rather than by a reader assuming the two agree; they differ by a day for
    every frame in the YZ Cnc season.
    """
    return (jd_to_datetime(jd)
            + dt.timedelta(hours=utc_offset_h)).strftime("%Y-%m-%d")


def night_label_to_utc_date(night: str, shift_days: int = 1) -> str:
    """Map a local night label to its UTC date.

    ``shift_days=1`` is correct for observations taken after local midnight,
    which is every CV frame in the 2024 YZ Cnc season.  The parameter is
    explicit rather than hard-coded so a caller working with an
    evening-side data set cannot silently inherit the wrong convention.
    """
    d = dt.datetime.strptime(night, "%Y-%m-%d").date()
    return (d + dt.timedelta(days=shift_days)).strftime("%Y-%m-%d")


# ===========================================================================
# 2.  Parsers — one per source, each returning the same row shape
# ===========================================================================

@dataclass(frozen=True)
class Obs:
    """One external photometric point, normalised across all sources.

    ``is_own`` is carried on the row itself rather than recomputed later so
    that no aggregation can accidentally forget to ask.
    """

    jd: float
    mag: float
    err: Optional[float]
    band: str
    source: str
    observer: str = ""
    is_own: bool = False
    is_limit: bool = False          # 'fainter than' — a non-detection
    note: str = ""

    @property
    def utc_night(self) -> str:
        return utc_night(self.jd)


def is_own_observation(observer: str, comment: str) -> bool:
    """Is this AAVSO row OUR OWN photometry, resubmitted?

    EITHER test fires: the observer code, or a marker in the free-text
    comment.  Deliberately generous — a false positive costs one point of
    independent evidence, a false negative costs the independence of the
    whole verdict, and those two errors are not the same size.
    """
    if (observer or "").strip().upper() in OWN_OBSERVER_CODES:
        return True
    up = (comment or "").upper()
    return any(marker in up for marker in OWN_COMMENT_MARKERS)


def parse_aavso_delim(text: str, delimiter: str = "@@@",
                      source: str = "aavso") -> list[Obs]:
    """Parse the AAVSO VSX ``view=api.delim`` payload.

    The first line is the header; the column set is documented by AAVSO in
    the LCG v2 client and includes JD, mag, uncert, band, by, comment and
    fainterThan.  Rows whose magnitude will not parse are DROPPED rather
    than coerced: a magnitude that is not a number is not a faint
    magnitude, and substituting anything for it would put an invented point
    into a state classification.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split(delimiter)
    idx = {name.strip(): i for i, name in enumerate(header)}
    # A payload missing JD or mag is not a light curve; say so loudly rather
    # than returning an empty list that reads like "the star had no data".
    for required in ("JD", "mag"):
        if required not in idx:
            raise ValueError(
                f"AAVSO payload has no {required!r} column; got {header[:6]}")

    def cell(row: Sequence[str], name: str) -> str:
        i = idx.get(name, -1)
        return row[i].strip() if 0 <= i < len(row) else ""

    out: list[Obs] = []
    for line in lines[1:]:
        row = line.split(delimiter)
        try:
            jd = float(cell(row, "JD"))
            mag = float(cell(row, "mag"))
        except ValueError:
            continue
        if not (math.isfinite(jd) and math.isfinite(mag)):
            continue
        try:
            err = float(cell(row, "uncert"))
        except ValueError:
            err = None
        observer = cell(row, "by")
        comment = cell(row, "comment")
        out.append(Obs(
            jd=jd, mag=mag, err=err, band=cell(row, "band") or "?",
            source=source, observer=observer,
            is_own=is_own_observation(observer, comment),
            is_limit=cell(row, "fainterThan") == "1",
            note=cell(row, "obsType")))
    return out


def parse_ztf_csv(text: str, source: str = "ztf") -> list[Obs]:
    """Parse the IRSA ZTF ``nph_light_curves`` CSV.

    ``catflags != 0`` marks a point the ZTF pipeline itself flagged
    (contamination, bad pixels, saturation-adjacent).  They are KEPT but
    marked in ``note``, because for coverage accounting "ZTF looked and got
    a flagged measurement" is a different statement from "ZTF did not look",
    and the report needs to be able to tell them apart.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    header = [h.strip() for h in lines[0].split(",")]
    idx = {name: i for i, name in enumerate(header)}
    if "mjd" not in idx or "mag" not in idx:
        raise ValueError(f"ZTF payload has no mjd/mag column; got {header[:8]}")
    out: list[Obs] = []
    for line in lines[1:]:
        row = line.split(",")
        try:
            mjd = float(row[idx["mjd"]])
            mag = float(row[idx["mag"]])
        except (ValueError, IndexError):
            continue
        try:
            err = float(row[idx["magerr"]])
        except (ValueError, IndexError, KeyError):
            err = None
        flags = row[idx["catflags"]].strip() if "catflags" in idx else "0"
        band = (row[idx["filtercode"]].strip()
                if "filtercode" in idx else "?")
        out.append(Obs(jd=mjd_to_jd(mjd), mag=mag, err=err, band=band,
                       source=source, observer="ZTF",
                       note="" if flags in ("", "0") else f"catflags={flags}"))
    return out


def parse_asassn_search(html: str) -> list[dict]:
    """Parse the ASAS-SN Variable Stars Database cone-search results table.

    Returns one dict per candidate with ``uuid``, ``asassn_name``,
    ``other_names``, ``sep_arcsec``, ``mean_vmag``, ``amplitude`` and
    ``var_type``.  PURE, so a captured page can be replayed in a test.

    Why parse the table at all instead of grabbing the first UUID in the
    page: the search radius that actually finds these polars is 0.5 ARCMIN
    (the parameter is arcmin, not degrees -- a 0.01 that looked like a
    36 arcsec cone was really 0.6 arcsec and silently missed three of five
    targets).  A cone that wide can admit a neighbour, so the caller has to
    be able to check the separation and the cross-identification before
    accepting a match.  Taking whatever UUID appeared first would have
    turned a wrong star into a published light curve.
    """
    import re
    body = re.search(r"<tbody.*?</tbody>", html, re.S)
    if not body:
        return []
    out: list[dict] = []
    for row in re.findall(r"<tr.*?</tr>", body.group(0), re.S):
        uuid_m = re.search(r"/variables/([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}"
                           r"-[0-9a-f]{4}-[0-9a-f]{12})", row)
        cells = [re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", c)).strip()
                 for c in re.findall(r"<td.*?</td>", row, re.S)]
        if not uuid_m or len(cells) < 9:
            continue

        def num(i):
            try:
                return float(cells[i])
            except (ValueError, IndexError):
                return None
        out.append(dict(uuid=uuid_m.group(1), asassn_name=cells[0],
                        other_names=cells[1], ra=num(2), dec=num(3),
                        sep_arcsec=num(4), mean_vmag=num(5),
                        amplitude=num(6), var_type=cells[8]))
    return out


def pick_asassn_match(candidates: Sequence[dict], target_name: str,
                      max_sep_arcsec: float = 15.0
                      ) -> tuple[Optional[dict], str]:
    """Choose the ASAS-SN record that IS this star, or refuse to choose.

    Two independent checks, and the name check is the one that matters:
    ASAS-SN publishes cross-identifications, so a record that names our
    target in its "other names" column is the star itself rather than a
    plausible neighbour.  Separation alone would be a guess dressed as a
    match in an 8 arcsec-pixel survey.
    """
    if not candidates:
        return None, "cone search returned no variable-star record"
    want = re.sub(r"\s+", "", target_name).upper()
    named = [c for c in candidates
             if want in re.sub(r"\s+", "", c.get("other_names", "")).upper()]
    if named:
        best = min(named, key=lambda c: c.get("sep_arcsec") or 1e9)
        return best, (f"cross-identified as {target_name} in the ASAS-SN "
                      f"record's own name list, {best.get('sep_arcsec')} "
                      f"arcsec from the VSX position")
    near = [c for c in candidates
            if (c.get("sep_arcsec") or 1e9) <= max_sep_arcsec]
    if not near:
        return None, (f"nearest record is "
                      f"{candidates[0].get('sep_arcsec')} arcsec away and "
                      f"does not name {target_name} — not accepted")
    best = min(near, key=lambda c: c.get("sep_arcsec") or 1e9)
    return best, (f"POSITIONAL MATCH ONLY: {best['asassn_name']} at "
                  f"{best.get('sep_arcsec')} arcsec; the record does not "
                  f"name {target_name}, so this identification is weaker "
                  f"than the others on this page")


def parse_asassn_json(payload: dict, source: str = "asassn") -> list[Obs]:
    """Parse the ASAS-SN Variable Stars Database ``.json`` export.

    Rows are ``{hjd, camera, mag, mag_err, flux, flux_err}``.  The times are
    HELIOCENTRIC JD.  For state classification that is irrelevant — the
    heliocentric correction is at most ~8 minutes and a dwarf-nova outburst
    is a multi-day event — but it is recorded in ``note`` so nobody later
    mistakes these for the barycentric times S3 computes.
    """
    rows = payload.get("results", []) if isinstance(payload, dict) else []
    out: list[Obs] = []
    for r in rows:
        try:
            jd = float(r["hjd"])
            mag = float(r["mag"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (math.isfinite(jd) and math.isfinite(mag)):
            continue
        try:
            err = float(r.get("mag_err"))
        except (TypeError, ValueError):
            err = None
        out.append(Obs(jd=jd, mag=mag, err=err, band="V", source=source,
                       observer=str(r.get("camera", "")),
                       note="HJD (heliocentric), not barycentric"))
    return out


# ===========================================================================
# 3.  Aggregation — nights, and the independence split
# ===========================================================================

@dataclass(frozen=True)
class NightPoint:
    """What the external record says about one UTC night."""

    night: str
    mag: float                  # representative magnitude (median)
    n: int
    n_independent: int
    bands: tuple[str, ...]
    observers: tuple[str, ...]
    spread: float               # max - min, a crude variability/scatter flag
    independent: bool           # did INDEPENDENT observers carry this night?


def _median(xs: Sequence[float]) -> float:
    """Median without importing numpy — this module stays dependency-free."""
    s = sorted(xs)
    n = len(s)
    if n == 0:
        raise ValueError("median of empty sequence")
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def v_like(obs: Iterable[Obs]) -> list[Obs]:
    """Just the rows whose band is a usable V stand-in (see V_LIKE_BANDS)."""
    return [o for o in obs if o.band in V_LIKE_BANDS]


def band_offset(obs: Iterable[Obs], band: str,
                reference_bands: frozenset[str] = V_LIKE_BANDS
                ) -> tuple[Optional[float], int, float]:
    """Measure ``band`` minus V, on nights where BOTH were observed.

    Returns ``(median_offset, n_nights, scatter)``.  This is how a Sloan
    ``SR`` magnitude earns the right to be read on the same ladder as a
    visual V estimate — by measurement on shared nights, not by assumption.
    ``None`` when no night carries both, which is itself the answer: the
    band cannot be placed on the ladder and must not be used as if it could.
    """
    per_night_band: dict[str, list[float]] = {}
    per_night_ref: dict[str, list[float]] = {}
    for o in obs:
        if o.is_limit:
            continue
        if o.band == band:
            per_night_band.setdefault(o.utc_night, []).append(o.mag)
        elif o.band in reference_bands:
            per_night_ref.setdefault(o.utc_night, []).append(o.mag)
    diffs = [_median(per_night_band[n]) - _median(per_night_ref[n])
             for n in sorted(set(per_night_band) & set(per_night_ref))]
    if not diffs:
        return None, 0, float("nan")
    med = _median(diffs)
    scatter = _median([abs(d - med) for d in diffs]) * 1.4826
    return med, len(diffs), scatter


def nightly(obs: Iterable[Obs], bands: Optional[frozenset[str]] = None,
            independent_only: bool = False) -> list[NightPoint]:
    """Collapse observations into one representative point per UTC night.

    The median is used rather than the mean because a dwarf nova genuinely
    varies within a night (that is the whole point of our 8 s cadence) and
    one bright flicker should not drag a night's state tag upward.
    ``spread`` preserves the within-night range that the median discards.
    """
    buckets: dict[str, list[Obs]] = {}
    for o in obs:
        if o.is_limit:
            continue
        if bands is not None and o.band not in bands:
            continue
        if independent_only and o.is_own:
            continue
        buckets.setdefault(o.utc_night, []).append(o)
    out: list[NightPoint] = []
    for night in sorted(buckets):
        rows = buckets[night]
        mags = [r.mag for r in rows]
        out.append(NightPoint(
            night=night, mag=_median(mags), n=len(rows),
            n_independent=sum(1 for r in rows if not r.is_own),
            bands=tuple(sorted({r.band for r in rows})),
            observers=tuple(sorted({r.observer for r in rows if r.observer})),
            spread=max(mags) - min(mags),
            independent=any(not r.is_own for r in rows)))
    return out


# ===========================================================================
# 4.  The ladder — quiescence measured, not asserted
# ===========================================================================

def quiescent_baseline(points: Sequence[NightPoint],
                       decile: float = QUIESCENT_DECILE) -> Optional[float]:
    """The star's own faint level: the median of its faintest ``decile``.

    Measured from the nightly points rather than the raw observations so a
    single densely-sampled night cannot define the baseline by weight of
    numbers.  Returns ``None`` when there are too few nights to speak — a
    baseline from two points is a guess wearing a statistic's clothes.
    """
    mags = sorted((p.mag for p in points), reverse=True)   # faintest first
    if len(mags) < 5:
        return None
    k = max(1, int(round(len(mags) * decile)))
    return _median(mags[:k])


def classify_night(mag: Optional[float], baseline: float) -> str:
    """Amplitude above the measured baseline -> a state tag.

    Three outcomes, not two.  A star caught mid-rise is genuinely neither
    quiescent nor in outburst, and ELEVATED says so instead of rounding to
    whichever neighbour the threshold happens to favour.
    """
    if mag is None or not math.isfinite(mag):
        return STATE_UNKNOWN
    amp = baseline - mag                       # positive == brighter
    if amp >= AMP_OUTBURST_MIN:
        return STATE_OUTBURST
    if amp <= AMP_QUIESCENT_MAX:
        return STATE_QUIESCENT
    return STATE_ELEVATED


@dataclass(frozen=True)
class Episode:
    """One contiguous outburst, and the verdict on what kind it was.

    ``duration_d`` is the whole above-threshold excursion; ``plateau_d`` is
    the part of it spent within :data:`PLATEAU_DEPTH_MAG` of the peak.  Both
    are reported because they answer different questions -- the first is
    "how long was this star active", the second is "was this a
    superoutburst" -- and only the second grades the episode.
    """

    start: str
    end: str
    peak_night: str
    peak_amp: float
    duration_d: float
    plateau_start: str
    plateau_end: str
    plateau_d: float
    n_nights: int
    kind: str
    why: str


def plateau_span(points: Sequence[NightPoint], peak_night: str,
                 depth: float = PLATEAU_DEPTH_MAG
                 ) -> tuple[str, str, float]:
    """The contiguous stretch around the peak that stays within ``depth``.

    Walks outward from the peak night and stops at the first night that has
    faded by more than ``depth``.  Walking outward (rather than taking every
    night in the episode that happens to be bright enough) is what keeps a
    post-superoutburst rebrightening from being counted as part of the
    plateau it is separated from.

    ``points`` must be THE EPISODE'S OWN NIGHTS, not the whole season.  Given
    the season, the walk escapes the event entirely: a low-amplitude peak at
    V=13.5 has almost every night in a dwarf nova's record sitting within
    1.5 mag of it, and an early version of this function duly reported an
    887-day plateau for a one-night outburst.  A plateau is part of an
    event, so its span is bounded by that event.
    """
    by_night = {p.night: p for p in points}
    order = sorted(by_night)
    if peak_night not in by_night:
        return peak_night, peak_night, 1.0
    i = order.index(peak_night)
    peak_mag = by_night[peak_night].mag
    lo = hi = i
    while lo - 1 >= 0 and by_night[order[lo - 1]].mag <= peak_mag + depth:
        lo -= 1
    while hi + 1 < len(order) and by_night[order[hi + 1]].mag <= peak_mag + depth:
        hi += 1
    start, end = order[lo], order[hi]
    return start, end, float(days_between(start, end)) + 1.0


def find_episodes(points: Sequence[NightPoint], baseline: float,
                  max_gap_d: float = EPISODE_MAX_GAP_D) -> list[Episode]:
    """Group outburst nights into episodes and grade each one.

    An episode is a run of OUTBURST nights separated by no more than
    ``max_gap_d``.  ``duration_d`` is measured from the first to the last
    outburst night INCLUSIVE, which is a lower bound on the true event
    length: the rise before the first detection and the decay after the last
    are not counted, because they were not observed.

    The GRADE, however, uses the plateau (see :func:`plateau_span`), not the
    whole excursion.  On a star that is rarely at true quiescence the
    above-threshold run merges an event with its own decline and with
    whatever follows within ``max_gap_d`` — which is how the first build of
    this function produced a "38-day superoutburst".  Both numbers are
    carried on the Episode so a reader can see the difference rather than
    take the grade on trust.
    """
    ob = [p for p in points
          if classify_night(p.mag, baseline) == STATE_OUTBURST]
    if not ob:
        return []
    runs: list[list[NightPoint]] = [[ob[0]]]
    for prev, cur in zip(ob, ob[1:]):
        gap = (dt.datetime.strptime(cur.night, "%Y-%m-%d")
               - dt.datetime.strptime(prev.night, "%Y-%m-%d")).days
        if gap <= max_gap_d:
            runs[-1].append(cur)
        else:
            runs.append([cur])
    episodes: list[Episode] = []
    for run in runs:
        peak = min(run, key=lambda p: p.mag)
        peak_amp = baseline - peak.mag
        span = (dt.datetime.strptime(run[-1].night, "%Y-%m-%d")
                - dt.datetime.strptime(run[0].night, "%Y-%m-%d")).days + 1.0
        p_start, p_end, p_days = plateau_span(run, peak.night)
        bright = peak_amp >= SUPEROUTBURST_AMP_MIN
        long_ = p_days >= SUPEROUTBURST_DAYS_MIN
        if bright and long_:
            kind, why = (EPISODE_SUPEROUTBURST,
                         f"peak {peak_amp:.2f} mag above quiescence "
                         f"(>= {SUPEROUTBURST_AMP_MIN:.1f}) AND a plateau of "
                         f"{p_days:.0f} d within {PLATEAU_DEPTH_MAG:.1f} mag "
                         f"of peak, {p_start} to {p_end} "
                         f"(>= {SUPEROUTBURST_DAYS_MIN:.0f})")
        else:
            missing = []
            if not bright:
                missing.append(f"peak only {peak_amp:.2f} mag above "
                               f"quiescence (< {SUPEROUTBURST_AMP_MIN:.1f})")
            if not long_:
                missing.append(f"plateau only {p_days:.0f} d "
                               f"(< {SUPEROUTBURST_DAYS_MIN:.0f})")
            kind, why = EPISODE_NORMAL, "; ".join(missing)
        episodes.append(Episode(
            start=run[0].night, end=run[-1].night, peak_night=peak.night,
            peak_amp=peak_amp, duration_d=span, plateau_start=p_start,
            plateau_end=p_end, plateau_d=p_days, n_nights=len(run),
            kind=kind, why=why))
    return episodes


# ===========================================================================
# 5.  The gating question, answered in one place
# ===========================================================================

@dataclass(frozen=True)
class NightVerdict:
    """The state of ONE of our observing nights, and what backs it."""

    utc_night: str
    local_night: str
    n_frames: int
    state: str
    mag: Optional[float]
    amp: Optional[float]
    basis: str                  # 'independent', 'own', 'bracketed', 'none'
    evidence: str
    episode: str = ""


def bracket(points: Sequence[NightPoint], night: str
            ) -> tuple[Optional[NightPoint], Optional[NightPoint]]:
    """The nearest external nights before and after ``night``.

    What makes a verdict possible on a night nobody else observed: a star
    cannot be 4 magnitudes up between two points that are both near
    quiescence, because a superoutburst does not fit in the gap.
    """
    d = dt.datetime.strptime(night, "%Y-%m-%d")
    before = [p for p in points
              if dt.datetime.strptime(p.night, "%Y-%m-%d") < d]
    after = [p for p in points
             if dt.datetime.strptime(p.night, "%Y-%m-%d") > d]
    return (before[-1] if before else None, after[0] if after else None)


def days_between(a: str, b: str) -> int:
    """Whole days from date ``a`` to date ``b`` (both ``YYYY-MM-DD``)."""
    return (dt.datetime.strptime(b, "%Y-%m-%d")
            - dt.datetime.strptime(a, "%Y-%m-%d")).days


def superoutburst_excluded_by_brackets(
        before: Optional[NightPoint], after: Optional[NightPoint],
        baseline: float, gap_limit_d: float = SUPEROUTBURST_DAYS_MIN
) -> tuple[bool, str]:
    """Can a superoutburst be ruled out on an unobserved night, from its
    neighbours alone?

    The argument is physical and it does not need our own photometry.  An
    SU UMa superoutburst holds a plateau for at least ``gap_limit_d`` days.
    So if the nearest external points on BOTH sides are close enough in time
    that a plateau could not have fitted between them, and both are fainter
    than superoutburst level, then no superoutburst occurred in the gap.
    This is the reasoning that keeps the YZ Cnc verdict independent of the
    RLMT rows AAVSO happens to hold.
    """
    if before is None or after is None:
        return False, "no external point on one side — nothing to bracket with"
    span = days_between(before.night, after.night)
    if span > gap_limit_d:
        return False, (f"external points {span} d apart — a "
                       f"{gap_limit_d:.0f} d plateau could hide in that gap")
    amp_b = baseline - before.mag
    amp_a = baseline - after.mag
    if amp_b >= SUPEROUTBURST_AMP_MIN or amp_a >= SUPEROUTBURST_AMP_MIN:
        return False, ("a bracketing point is itself at superoutburst "
                       "level — the gap is inside an event, not outside one")
    return True, (
        f"bracketed by {before.night} (V={before.mag:.2f}, "
        f"{amp_b:+.2f} mag above quiescence) and {after.night} "
        f"(V={after.mag:.2f}, {amp_a:+.2f}) — only {span} d apart, so a "
        f">= {gap_limit_d:.0f} d superoutburst plateau cannot fit between "
        f"them, and neither endpoint is anywhere near superoutburst level")


def branch_recommendation(verdicts: Sequence[NightVerdict],
                          dense_min_frames: int) -> tuple[str, str]:
    """Which branch should the YZ Cnc superhump task take?

    Returns ``(branch, reasoning)``.  The rule is the physics: common
    superhumps are a superoutburst phenomenon in SU UMa systems, so the
    superhump branch is available if and only if a DENSE run sits inside an
    episode graded SUPEROUTBURST.  Anything else and the strategy's stated
    fallback is the honest branch — with the outburst nights, if there are
    any, called out separately, because a normal outburst caught at 8 s
    cadence is a real data set even though it is not a superhump data set.
    """
    dense = [v for v in verdicts if v.n_frames >= dense_min_frames]
    if not dense:
        return ("UNDECIDED",
                f"no run reaches {dense_min_frames} frames — nothing here "
                f"could carry a period analysis either way")
    in_super = [v for v in dense if v.episode == EPISODE_SUPEROUTBURST]
    if in_super:
        nights = ", ".join(v.utc_night for v in in_super)
        return ("SUPERHUMP",
                f"{len(in_super)} of {len(dense)} dense runs fall inside an "
                f"episode graded {EPISODE_SUPEROUTBURST} ({nights}); common "
                f"superhumps are expected and the period analysis is on.")

    # Count by STATE, which every night has, rather than by episode
    # membership, which only nights the INDEPENDENT record traced can have.
    # Counting by episode silently reported four outburst nights as
    # "elevated" simply because no outside observer covered them.
    def n(state):
        return sum(1 for v in dense if v.state == state)

    ob = [v for v in dense if v.state == STATE_OUTBURST]
    detail = (f"{len(dense)} dense runs: {n(STATE_QUIESCENT)} quiescent, "
              f"{len(ob)} in outburst, {n(STATE_ELEVATED)} elevated, "
              f"{n(STATE_UNKNOWN)} unclassifiable.")
    if ob:
        peak = min((v for v in ob if v.amp is not None),
                   key=lambda v: -v.amp, default=None)
        detail += (f" The brightest dense run reaches "
                   f"{peak.amp:.2f} mag above quiescence on {peak.utc_night}, "
                   f"short of the {SUPEROUTBURST_AMP_MIN:.1f} mag a "
                   f"superoutburst reaches." if peak else "")
    return ("FALLBACK",
            f"No dense run falls inside a superoutburst. {detail} Common "
            f"superhumps are a superoutburst phenomenon in SU UMa systems, "
            f"so the superhump branch has no data to stand on; the "
            f"strategy's orbital-hump + flickering fallback is the branch "
            f"this season supports — and the normal-outburst nights are a "
            f"distinct, separately publishable data set, not a consolation.")


# ===========================================================================
# 6.  Survey coverage accounting
# ===========================================================================

@dataclass(frozen=True)
class Coverage:
    """What one source delivered for one target."""

    target: str
    source: str
    n_points: int
    mjd_min: Optional[float]
    mjd_max: Optional[float]
    bands: tuple[str, ...]
    n_nights: int
    median_gap_d: Optional[float]
    notes: str

    @property
    def span_d(self) -> Optional[float]:
        if self.mjd_min is None or self.mjd_max is None:
            return None
        return self.mjd_max - self.mjd_min


def coverage(target: str, source: str, obs: Sequence[Obs],
             notes: str = "") -> Coverage:
    """Summarise a source's delivery for a target: span, cadence, bands.

    ``median_gap_d`` is the median interval between consecutive OBSERVED
    NIGHTS, which is the honest way to state a cadence for a record that is
    dense in some seasons and absent in others — a mean gap would be
    dominated by the annual solar conjunction gap and would describe no
    actual observing behaviour.
    """
    if not obs:
        return Coverage(target, source, 0, None, None, (), 0, None,
                        notes or "no rows returned")
    mjds = [jd_to_mjd(o.jd) for o in obs]
    nights = sorted({o.utc_night for o in obs})
    gaps = [float(days_between(a, b)) for a, b in zip(nights, nights[1:])]
    return Coverage(
        target=target, source=source, n_points=len(obs),
        mjd_min=min(mjds), mjd_max=max(mjds),
        bands=tuple(sorted({o.band for o in obs})),
        n_nights=len(nights),
        median_gap_d=_median(gaps) if gaps else None,
        notes=notes)


#: Bright-limit guidance per survey, in the band the survey works in.
#: These are the magnitudes at which each survey's photometry stops being
#: trustworthy from the BRIGHT side — the failure mode that matters here,
#: because several of these CVs are bright in outburst and one of them
#: (YZ Cnc at V~10.5) is far above two of these limits when it matters most.
#: Quoted from each survey's own documentation, and flagged as external
#: constants in the report footer rather than presented as measurements.
SURVEY_BRIGHT_LIMIT = {
    "ztf": (13.0, "ZTF photometry saturates around g,r ~ 12.5-13.5; brighter "
                  "epochs are flagged or dropped by the pipeline"),
    "asassn": (10.0, "ASAS-SN saturates near V,g ~ 10-11 in its 8 arcsec "
                     "pixels"),
    "aavso": (None, "no bright limit — amateur observers switch equipment "
                    "and exposure to suit the star"),
}


def saturation_risk(source: str, brightest_mag: Optional[float]
                    ) -> tuple[bool, str]:
    """Would this source have struggled at this star's bright end?

    Returns ``(at_risk, explanation)``.  Called with the target's brightest
    known magnitude so the report can state, per target and per survey,
    whether an absence of bright-epoch points is a real absence or a
    saturation artefact.  Getting this backwards would let a saturated
    survey's silence be read as "the star was not in outburst".
    """
    limit, why = SURVEY_BRIGHT_LIMIT.get(source, (None, "no limit recorded"))
    if limit is None or brightest_mag is None:
        return False, why
    if brightest_mag <= limit:
        return True, (f"target reaches {brightest_mag:.1f} mag, at or above "
                      f"the ~{limit:.1f} mag bright limit — {why}")
    return False, (f"target's bright end ({brightest_mag:.1f} mag) stays "
                   f"below the ~{limit:.1f} mag limit")
