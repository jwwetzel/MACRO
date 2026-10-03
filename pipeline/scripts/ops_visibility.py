#!/usr/bin/env python
"""Visibility of the observing-request targets from Winer, Oct 2026 -> Apr 2027.

WHY THIS SCRIPT EXISTS
----------------------
Revision 2 of the observatory request (``ops/2026-08_observatory_request.md``)
asked for a weekly >= 2 h T CrB run "resuming at re-opening" and an ST LMi
season "starting at re-opening".  Neither is possible: the committee's plan
review of 2026-10-03 (OA section 2, TE.F7) showed T CrB is above airmass 2 in
nautical darkness for about an hour in early October and not at all from late
October to mid-December, and ST LMi does not give a useful night before about
November 1.  A request that asks the site for the impossible is not read
twice.  Revision 3 therefore carries a visibility table, and the house rule is
that a published number is emitted by a script, never typed -- this is that
script.

WHAT IT COMPUTES
----------------
For every target and every night, the number of hours for which BOTH

* the target is above ``--airmass`` (default 2.0, i.e. altitude > 30 deg;
  plane-parallel sec z, which at z = 60 deg differs from Pickering's formula
  by 0.2 % -- irrelevant at the 0.1 h precision quoted), and
* the Sun is below each twilight limit in ``SUN_LIMITS_DEG``
  (-12 deg nautical, -15 deg the grism limit, -18 deg astronomical).

The -15 deg row exists because of OA's remark that twilight fills the
slitless sky lozenge: a grism exposure wants a darker sky than photometry of
a 10th-magnitude star does.

A "night" is labelled by the LOCAL date on which it begins: the window
scanned is local noon to the following local noon (Winer keeps MST = UTC-7
all year, so 19:00 UT to 19:00 UT), sampled every ``STEP_MIN`` minutes.  For
each (target, night) the script also records the longest CONTINUOUS window at
the nautical limit, its UT start and end, whether it is an evening or a
morning window, and the airmass range inside it -- "70 minutes" at airmass
1.7-2.0 in evening twilight is a different observation from 70 minutes at
airmass 1.1.

From the nightly scan it then derives the dates the request quotes: the first
and last night each target offers at least 0.25 h (a nightly snapshot block),
1.5 h (a whole polar orbit) and 2.0 h (a flickering run).

METHOD, AND WHY NOT ``AltAz``
----------------------------
212 nights x 720 samples x ~35 targets is 5.3 million altitude evaluations.
``SkyCoord.transform_to(AltAz)`` would do it, slowly; the spherical triangle
does it in a second:

    sin(alt) = sin(lat) sin(dec) + cos(lat) cos(dec) cos(LAST - ra)

with ``ra, dec`` in the TRUE equator and equinox of date (astropy ``TETE``;
precession since J2000 is 0.37 deg by 2027 and is not ignored) and LAST the
local apparent sidereal time.  ``pipeline/tests/test_ops_visibility.py`` pins
this against astropy's full ``AltAz`` transform to < 0.02 deg, and pins the
headline cells against the observational astronomer's independent run.
Refraction is ignored: it is < 0.03 deg at 30 deg altitude, and at the Sun
limits it is defined out (the limits are geometric by convention).

Earth orientation for 2027 is not yet tabulated; UT1-UTC is extrapolated by
astropy (|error| < 1 s, i.e. < 0.005 deg of hour angle).  The script never
touches the network.

WHERE THE COORDINATES COME FROM
-------------------------------
Targets the archive has already observed take their position from the
manifest: the median ``ra_deg, dec_deg`` over that target's frames, read
read-only.  (A median, because individual headers can be stale -- OA.E2.)
Targets the archive has never observed -- the new standards -- are listed in
``NEW_TARGETS`` with catalogue (SIMBAD, ICRS J2000) coordinates.  Every row of
the output says which it was.

OUTPUTS (all under ``ops/generated/``)
--------------------------------------
``visibility.csv``            one row per (target, tabulated night)
``visibility_table.md``       the hours table, nautical / astronomical
``visibility_windows.md``     longest nautical window per tabulated night
``visibility_dates.md``       first / last usable night per threshold
``fig_visibility.{png,pdf}``  hours vs date, one panel per programme

``--inject FILE`` rewrites the blocks delimited by
``<!-- BEGIN GENERATED: name -->`` / ``<!-- END GENERATED: name -->`` in FILE
with the freshly generated fragments, so the request document itself never
holds a hand-typed visibility number.

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    $PY pipeline/scripts/ops_visibility.py
    $PY pipeline/scripts/ops_visibility.py --inject ops/2026-10_observatory_request_rev3.md
    # the eruption block relaxes the airmass limit to 3:
    $PY pipeline/scripts/ops_visibility.py --airmass 3 \
        --out-dir ops/eruption_block/generated \
        --inject ops/eruption_block/README.md
"""

from __future__ import annotations

import argparse
import csv
import re
import sqlite3
import sys
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "pipeline"))

from macro_core import plotstyle                     # noqa: E402
from macro_core import timing                        # noqa: E402

OPS_VISIBILITY_VERSION = "1.0 (2026-10-03)"

MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"
OUT_DIR = REPO_ROOT / "ops" / "generated"

#: Winer keeps Mountain Standard Time all year (Arizona observes no DST).
UTC_OFFSET_H = -7.0
#: Sampling step of the nightly scan.  Two minutes resolves every quantity
#: quoted (hours to one decimal) with an order of magnitude to spare.
STEP_MIN = 2.0
#: Sun-altitude limits [deg]: nautical, grism (OA section 2), astronomical.
SUN_LIMITS_DEG = (-12.0, -15.0, -18.0)
#: First and last LOCAL evening date of the nightly scan.
SCAN_START = date(2026, 10, 1)
SCAN_END = date(2027, 4, 30)
#: The nights printed in the request's table.  The first five columns of
#: OA's table (Oct 5, Oct 15, Nov 1, Dec 15, Jan 1, Jan 15, Mar 1) are a
#: subset, so the two can be compared cell for cell.
TABLE_NIGHTS = (date(2026, 10, 5), date(2026, 10, 15), date(2026, 10, 25),
                date(2026, 11, 1), date(2026, 11, 15), date(2026, 12, 1),
                date(2026, 12, 15), date(2027, 1, 1), date(2027, 1, 15),
                date(2027, 2, 1), date(2027, 3, 1), date(2027, 4, 1))
#: Thresholds [h] for the "first / last usable night" table, and what each
#: one is the minimum for.
DATE_THRESHOLDS_H = ((0.25, "snapshot block"),
                     (1.5, "one polar orbit"),
                     (2.0, "flickering / whole-orbit run"))


@dataclass(frozen=True)
class Target:
    """One row of the request's target list."""
    name: str            #: display name
    programme: str       #: which paper / purpose asks for it
    role: str            #: science | standard | calibration
    manifest_key: str | None = None   #: ``frames.canonical_target``, if any
    ra_deg: float | None = None       #: ICRS, used when not in the manifest
    dec_deg: float | None = None


#: Targets the archive already holds: position from the manifest.
ARCHIVE_TARGETS = (
    Target("T CrB", "T CrB", "science", "T CrB"),
    Target("θ CrB", "T CrB", "standard", "tet CrB"),
    Target("Vega", "T CrB / Be", "standard", "Vega"),
    Target("ST LMi", "CV", "science", "ST LMi"),
    Target("VV Pup", "CV", "science", "VV Pup"),
    Target("EU UMa", "CV", "science", "EU UMa"),
    Target("AN UMa", "CV", "science", "AN UMa"),
    Target("YZ Cnc", "CV", "science", "YZ Cnc"),
    Target("λ Eri", "Be", "science", "Lam Eri"),
    Target("69 Ori", "Be", "science", "69 Ori"),
    Target("5 Cnc", "Be", "science", "5 Cnc"),
    Target("HD 70340", "Be", "science", "HD 70340"),
    Target("53 Boo", "Be", "science", "53 Boo"),
    Target("Phecda", "Be", "science", "PHECDA"),
    Target("φ Leo", "Be", "science", "phi Leo"),
    Target("Spica", "Be", "standard", "Spica"),
    Target("η Hya (HR 3454)", "Be", "standard", "HR 3454"),
    Target("θ Vir (HR 4963)", "Be", "standard", "HR 4963"),
    Target("θ Crt (HR 4468)", "Be", "standard", "HR 4468"),
    Target("M101", "SN 2023ixf", "science", "M101"),
    Target("NGC 5548", "Dwarf / AGN", "science", "NGC 5548"),
    Target("NGC 7027", "grism calibration", "calibration", "NGC 7027"),
)

#: Targets new to the archive: catalogue positions (SIMBAD, ICRS J2000).
NEW_TARGETS = (
    # CALSPEC flux standard 7 deg from T CrB, V = 10.8 (OA section 2).
    Target("BD+33°2642", "T CrB", "standard", None, 237.99954, 32.94842),
    # pi2 Ori, A1 Vn, the autumn nightly standard beside lambda Eri (OA).
    Target("π² Ori (HR 1544)", "Be", "standard", None, 72.65301, 8.90018),
    # Compact planetary nebulae: emission-line wavelength / LSF calibrators.
    Target("IC 418", "grism calibration", "calibration", None,
           81.86750, -12.69731),
    # Equatorial Landolt / Smith et al. u'g'r'i'z' fields inside the SDSS
    # footprint: one per half of the night, for colour terms.
    Target("SA 95 (SDSS field)", "photometric calibration", "calibration",
           None, 58.50000, 0.00000),
    Target("SA 98 (SDSS field)", "photometric calibration", "calibration",
           None, 103.02083, -0.32917),
)

#: Figure panels: (title, target names).  One programme per panel so the
#: reader compares targets that compete for the same hours.
PANELS = (
    ("T CrB and its standards", ("T CrB", "θ CrB", "BD+33°2642", "Vega")),
    ("Polars and YZ Cnc", ("ST LMi", "VV Pup", "EU UMa", "AN UMa", "YZ Cnc")),
    ("Be stars and autumn standards",
     ("λ Eri", "69 Ori", "5 Cnc", "π² Ori (HR 1544)", "η Hya (HR 3454)",
      "Phecda", "53 Boo")),
    ("Galaxies and calibrators",
     ("M101", "NGC 5548", "IC 418", "NGC 7027", "SA 95 (SDSS field)",
      "SA 98 (SDSS field)")),
)


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------
def manifest_position(con: sqlite3.Connection, key: str) -> tuple:
    """Median ``(ra_deg, dec_deg, n_frames)`` of one archive target.

    The median, not the mean and not the first row: OA.E2 found frames whose
    header position was left over from the previous slew (T CrB spectra
    labelled RA 05 29, Dec -17).  A handful of those cannot move a median.
    RA is median-ed on the unit circle so a target near 0h cannot average to
    12h.
    """
    rows = con.execute(
        "SELECT ra_deg, dec_deg FROM frames WHERE canonical_target = ? "
        "AND ra_deg IS NOT NULL AND dec_deg IS NOT NULL", (key,)).fetchall()
    if not rows:
        raise LookupError(f"manifest has no positioned frame for {key!r}")
    ra = np.radians([r[0] for r in rows])
    ra_med = np.degrees(np.arctan2(np.median(np.sin(ra)),
                                   np.median(np.cos(ra)))) % 360.0
    return float(ra_med), float(np.median([r[1] for r in rows])), len(rows)


def resolve_targets(manifest: Path = MANIFEST) -> list:
    """All targets with a position and a statement of where it came from.

    Returns a list of ``(Target, ra_deg, dec_deg, source)``.
    """
    con = sqlite3.connect(f"file:{manifest}?mode=ro", uri=True)
    try:
        out = []
        for t in ARCHIVE_TARGETS:
            ra, dec, n = manifest_position(con, t.manifest_key)
            out.append((t, ra, dec, f"manifest median of {n} frames"))
    finally:
        con.close()
    for t in NEW_TARGETS:
        out.append((t, t.ra_deg, t.dec_deg, "SIMBAD ICRS (new target)"))
    return out


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def _quiet_iers() -> None:
    """Forbid network access and accept extrapolated UT1 for 2027 dates."""
    import warnings
    from astropy.utils import iers
    from astropy.utils.exceptions import AstropyWarning
    # Polar motion is not tabulated for 2027; astropy falls back to the
    # 50-year mean and says so.  That is an arcsecond-level effect on a
    # quantity quoted to 0.1 h, so the warning is noise here.
    warnings.filterwarnings("ignore", category=AstropyWarning,
                            message=".*polar motions.*")
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    iers.conf.iers_degraded_accuracy = "ignore"


def night_grid(night: date, step_min: float = STEP_MIN):
    """UTC sample times (astropy ``Time``) for the night that BEGINS on the
    local date ``night``: local noon to the next local noon, exclusive."""
    from astropy.time import Time, TimeDelta
    import astropy.units as u
    t0 = Time(f"{night.isoformat()}T12:00:00", scale="utc") \
        - TimeDelta(UTC_OFFSET_H * u.hour)
    n = int(round(24 * 60 / step_min))
    return t0 + TimeDelta(np.arange(n) * step_min * u.min)


def altitude_deg(ra_deg, dec_deg, last_deg, lat_deg=timing.WINER_LAT_DEG):
    """Geometric altitude from the spherical triangle (no refraction).

    ``ra_deg, dec_deg`` must be in the true equator and equinox of date and
    ``last_deg`` the local APPARENT sidereal time, both in degrees.  Arrays
    broadcast.
    """
    lat = np.radians(lat_deg)
    dec = np.radians(dec_deg)
    ha = np.radians(np.asarray(last_deg) - np.asarray(ra_deg))
    sin_alt = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(ha)
    return np.degrees(np.arcsin(np.clip(sin_alt, -1.0, 1.0)))


def of_date(ra_deg, dec_deg, when):
    """ICRS -> true equator and equinox of date (``TETE``) at ``when``.

    Stars: precession + nutation + annual aberration.  One epoch per night is
    ample -- the frame moves 0.014 deg per year.
    """
    from astropy.coordinates import SkyCoord, TETE
    import astropy.units as u
    c = SkyCoord(np.atleast_1d(ra_deg) * u.deg, np.atleast_1d(dec_deg) * u.deg,
                 frame="icrs").transform_to(TETE(obstime=when))
    return c.ra.deg, c.dec.deg


def sun_and_last(times):
    """``(sun_alt_deg, last_deg)`` at Winer for an array of UTC times."""
    from astropy.coordinates import TETE, get_sun
    import astropy.units as u
    lon = timing.WINER_LON_DEG * u.deg
    last = times.sidereal_time("apparent", longitude=lon).deg
    sun = get_sun(times).transform_to(TETE(obstime=times))
    return altitude_deg(sun.ra.deg, sun.dec.deg, last), last


def airmass_to_alt_deg(airmass: float) -> float:
    """Altitude at which plane-parallel sec z equals ``airmass``."""
    return float(np.degrees(np.arcsin(1.0 / airmass)))


def longest_run(mask: np.ndarray) -> tuple:
    """``(start, stop)`` indices (stop exclusive) of the longest run of True
    in a boolean array; ``(0, 0)`` if there is none."""
    best = (0, 0)
    start = None
    for i, v in enumerate(np.append(mask, False)):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    return best


def night_visibility(night: date, ra_deg, dec_deg, airmass: float = 2.0,
                     step_min: float = STEP_MIN) -> list:
    """Visibility of N targets on one night.

    Returns one dict per target with ``hours[limit]`` for each Sun limit and
    the longest continuous nautical window (UT start/end, which half of
    the night, airmass range).
    """
    times = night_grid(night, step_min)
    sun_alt, last = sun_and_last(times)
    ra_d, dec_d = of_date(ra_deg, dec_deg, times[len(times) // 2])
    alt = altitude_deg(ra_d[:, None], dec_d[:, None], last[None, :])
    up = alt > airmass_to_alt_deg(airmass)
    dark = {lim: sun_alt < lim for lim in SUN_LIMITS_DEG}
    step_h = step_min / 60.0
    # Local midnight sits half-way through the grid: a window that ends
    # before it is an evening window, one that starts after it a morning
    # window, anything else straddles it.
    mid = len(times) // 2
    out = []
    for k in range(len(ra_d)):
        rec = {"hours": {lim: float(np.sum(up[k] & dark[lim]) * step_h)
                         for lim in SUN_LIMITS_DEG}}
        i0, i1 = longest_run(up[k] & dark[SUN_LIMITS_DEG[0]])
        if i1 > i0:
            z = 1.0 / np.sin(np.radians(alt[k, i0:i1]))
            rec.update(window_h=(i1 - i0) * step_h,
                       ut_start=times[i0].strftime("%H:%M"),
                       ut_end=times[i1 - 1].strftime("%H:%M"),
                       half=("evening" if i1 <= mid else
                             "morning" if i0 >= mid else "spans midnight"),
                       airmass_min=float(z.min()), airmass_max=float(z.max()))
        else:
            rec.update(window_h=0.0, ut_start="", ut_end="", half="",
                       airmass_min=float("nan"), airmass_max=float("nan"))
        out.append(rec)
    return out


def scan(targets: list, start: date = SCAN_START, end: date = SCAN_END,
         airmass: float = 2.0) -> dict:
    """Nightly scan.  Returns ``{"nights": [date...], "hours": {limit:
    array[n_targets, n_nights]}, "records": {night: [rec per target]}}``."""
    _quiet_iers()
    ra = np.array([t[1] for t in targets])
    dec = np.array([t[2] for t in targets])
    nights = [start + timedelta(days=i) for i in range((end - start).days + 1)]
    hours = {lim: np.zeros((len(targets), len(nights)))
             for lim in SUN_LIMITS_DEG}
    records = {}
    for j, night in enumerate(nights):
        recs = night_visibility(night, ra, dec, airmass)
        records[night] = recs
        for lim in SUN_LIMITS_DEG:
            hours[lim][:, j] = [r["hours"][lim] for r in recs]
    return {"nights": nights, "hours": hours, "records": records}


def usable_segments(nights: list, hours: np.ndarray, threshold_h: float,
                    max_gap: int = 7) -> list:
    """The apparitions: runs of nights with at least ``threshold_h`` hours.

    Returns ``[(first_night, last_night), ...]``.  Two qualifying nights
    separated by at most ``max_gap`` failing nights belong to one segment
    (a threshold crossed and re-crossed by rounding is not a new
    apparition); a longer gap -- solar conjunction -- splits the season in
    two, which is exactly what the request has to say about T CrB.
    """
    ok = np.flatnonzero(hours >= threshold_h - 1e-9)
    if ok.size == 0:
        return []
    segments = []
    first = prev = ok[0]
    for i in ok[1:]:
        if i - prev - 1 > max_gap:
            segments.append((nights[first], nights[prev]))
            first = i
        prev = i
    segments.append((nights[first], nights[prev]))
    return segments


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------
def _d(night: date) -> str:
    """'Oct 5' -- the year is carried by the table header."""
    return f"{night.strftime('%b')} {night.day}"


def _cell(h_naut: float, h_astro: float) -> str:
    """'1.1 / 0.6'; a bare '0' when the target is never usable."""
    if h_naut < 0.05:
        return "0"
    return f"{h_naut:.1f} / {h_astro:.1f}"


def _stamp(airmass: float) -> str:
    return (f"*Generated by `pipeline/scripts/ops_visibility.py` "
            f"v{OPS_VISIBILITY_VERSION}; Winer "
            f"{timing.WINER_LAT_DEG:+.4f}°, {timing.WINER_LON_DEG:+.4f}°; "
            f"airmass < {airmass:g}; {STEP_MIN:g}-min sampling; night = "
            f"local noon to noon (MST), labelled by its evening date. "
            f"Do not edit by hand.*")


def table_md(targets, result, airmass) -> str:
    """Hours per tabulated night, 'Sun < -12 / Sun < -18'."""
    nights = result["nights"]
    idx = [nights.index(n) for n in TABLE_NIGHTS]
    head = "| Target | Programme | " + " | ".join(_d(n) for n in TABLE_NIGHTS)
    lines = [head + " |",
             "|---|---|" + "---:|" * len(TABLE_NIGHTS)]
    for k, (t, *_rest) in enumerate(targets):
        cells = [_cell(result["hours"][-12.0][k, j],
                       result["hours"][-18.0][k, j]) for j in idx]
        lines.append(f"| {t.name} | {t.programme} | " + " | ".join(cells)
                     + " |")
    y0, y1 = TABLE_NIGHTS[0].year, TABLE_NIGHTS[-1].year
    lines += ["",
              f"Hours above airmass {airmass:g} with the Sun below −12° / "
              f"below −18°, nights of {y0}–{y1} (columns from Jan 1 onward "
              f"are {y1}). `0` = never above airmass {airmass:g} in nautical "
              f"darkness.", "", _stamp(airmass)]
    return "\n".join(lines)


def windows_md(targets, result, airmass, names) -> str:
    """Longest continuous nautical window per tabulated night, for the
    targets whose request text depends on WHEN in the night they are up."""
    lines = ["| Target | Night | Window (UT) | Length (h) | Half | Airmass "
             "in window | Sun < −15° (h) |",
             "|---|---|---|---:|---|---|---:|"]
    by_name = {t[0].name: k for k, t in enumerate(targets)}
    for name in names:
        k = by_name[name]
        for night in TABLE_NIGHTS:
            r = result["records"][night][k]
            if r["window_h"] < 0.05:
                lines.append(f"| {name} | {_d(night)} | — | 0 | — | — | 0 |")
                continue
            lines.append(
                f"| {name} | {_d(night)} | {r['ut_start']}–{r['ut_end']} | "
                f"{r['window_h']:.1f} | {r['half']} | "
                f"{r['airmass_min']:.2f}–{r['airmass_max']:.2f} | "
                f"{r['hours'][-15.0]:.1f} |")
    lines += ["", "Longest continuous interval above airmass "
              f"{airmass:g} with the Sun below −12°; the last column is the "
              "total with the Sun below −15°, the limit adopted for grism "
              "exposures.", "", _stamp(airmass)]
    return "\n".join(lines)


def dates_md(targets, result, airmass) -> str:
    """First / last usable night per duration threshold (Sun < -12 deg),
    plus the grism season (Sun < -15 deg, snapshot threshold)."""
    nights = result["nights"]

    def span(h, thr):
        segs = usable_segments(nights, h, thr)
        if not segs:
            return "never"
        parts = []
        for a, b in segs:
            left = "(open) " if a == nights[0] else ""
            right = " (open)" if b == nights[-1] else ""
            parts.append(f"{left}{a.isoformat()} → {b.isoformat()}{right}")
        return "; ".join(parts)

    head = "| Target | " + " | ".join(
        f"≥ {thr:g} h ({what})" for thr, what in DATE_THRESHOLDS_H) \
        + " | Grism: ≥ 0.25 h, Sun < −15° |"
    lines = [head, "|---|" + "---|" * (len(DATE_THRESHOLDS_H) + 1)]
    for k, (t, *_rest) in enumerate(targets):
        cells = [span(result["hours"][-12.0][k], thr)
                 for thr, _ in DATE_THRESHOLDS_H]
        cells.append(span(result["hours"][-15.0][k], 0.25))
        lines.append(f"| {t.name} | " + " | ".join(cells) + " |")
    lines += ["",
              f"First → last night, within the scan {nights[0].isoformat()} "
              f"→ {nights[-1].isoformat()}, offering at least the stated "
              f"hours above airmass {airmass:g} with the Sun below −12° "
              "(last column: below −15°). Two intervals = an evening and a "
              "morning apparition separated by solar conjunction. "
              "`(open)` = the interval runs past that end of the scan.",
              "", _stamp(airmass)]
    return "\n".join(lines)


def write_csv(path: Path, targets, result) -> None:
    """One row per (target, tabulated night); the machine-readable twin of
    the two Markdown tables."""
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["target", "programme", "role", "ra_icrs_deg",
                    "dec_icrs_deg", "position_source", "night_local",
                    "hours_sun_lt_m12", "hours_sun_lt_m15",
                    "hours_sun_lt_m18", "window_ut_start", "window_ut_end",
                    "window_h", "half", "airmass_min", "airmass_max"])
        for k, (t, ra, dec, src) in enumerate(targets):
            for night in TABLE_NIGHTS:
                r = result["records"][night][k]
                w.writerow([t.name, t.programme, t.role, f"{ra:.5f}",
                            f"{dec:.5f}", src, night.isoformat(),
                            f"{r['hours'][-12.0]:.2f}",
                            f"{r['hours'][-15.0]:.2f}",
                            f"{r['hours'][-18.0]:.2f}",
                            r["ut_start"], r["ut_end"],
                            f"{r['window_h']:.2f}", r["half"],
                            f"{r['airmass_min']:.3f}",
                            f"{r['airmass_max']:.3f}"])


def figure(path_stem: Path, targets, result, airmass) -> None:
    """Hours of nautical darkness above the airmass limit vs date.

    One panel per programme.  Colour AND line style per target (house
    rule: colour is never the only channel).  The two horizontal reference
    lines are the durations the request actually asks for.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    nights = result["nights"]
    x = mdates.date2num(nights)
    by_name = {t[0].name: k for k, t in enumerate(targets)}
    with plotstyle.context("web"):
        fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.4), sharex=True,
                                 sharey=True, constrained_layout=True)
        for ax, (title, names) in zip(axes.ravel(), PANELS):
            for i, name in enumerate(names):
                ax.plot(x, result["hours"][-12.0][by_name[name]],
                        label=name, lw=1.5, **plotstyle.line_series(i))
            for thr, what in DATE_THRESHOLDS_H[1:]:
                ax.axhline(thr, **plotstyle.reference_kw())
            ax.set_title(title)
            ax.legend(loc="upper left", ncol=2, fontsize="small",
                      frameon=False)
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
            ax.set_ylim(0, 12.5)
        for ax in axes[:, 0]:
            ax.set_ylabel(f"hours at airmass < {airmass:g},\nSun < −12° (h)")
        for ax in axes[-1, :]:
            ax.set_xlabel("night (local evening date)")
        fig.suptitle("What the sky allows from Winer after the October 2026 "
                     "re-opening (dashed: 1.5 h and 2 h)")
        fig.savefig(path_stem.with_suffix(".png"), dpi=plotstyle.WEB_DPI)
        fig.savefig(path_stem.with_suffix(".pdf"))
        plt.close(fig)


def inject(doc: Path, fragments: dict) -> list:
    """Replace each ``<!-- BEGIN GENERATED: name -->`` ... ``<!-- END
    GENERATED: name -->`` block in ``doc`` whose name is a key of
    ``fragments``.  Returns the names replaced.  Blocks with other names
    (owned by other generators) are left untouched."""
    text = doc.read_text()
    done = []
    for name, body in fragments.items():
        pat = re.compile(
            rf"(<!-- BEGIN GENERATED: {re.escape(name)} -->\n).*?"
            rf"(\n<!-- END GENERATED: {re.escape(name)} -->)", re.S)
        text, n = pat.subn(lambda m: m.group(1) + body + m.group(2), text)
        if n:
            done.append(name)
    doc.write_text(text)
    return done


#: Targets whose request text depends on the time of night they are up.
WINDOW_TARGETS = ("T CrB", "ST LMi", "λ Eri", "M101")


def build(airmass: float = 2.0, out_dir: Path = OUT_DIR,
          manifest: Path = MANIFEST) -> dict:
    """Run the scan and write every output.  Returns the Markdown fragments
    keyed by the name used in ``<!-- BEGIN GENERATED: ... -->`` markers."""
    targets = resolve_targets(manifest)
    result = scan(targets, airmass=airmass)
    out_dir.mkdir(parents=True, exist_ok=True)
    fragments = {
        "visibility_table": table_md(targets, result, airmass),
        "visibility_windows": windows_md(targets, result, airmass,
                                         WINDOW_TARGETS),
        "visibility_dates": dates_md(targets, result, airmass),
    }
    for name, body in fragments.items():
        (out_dir / f"{name}.md").write_text(body + "\n")
    write_csv(out_dir / "visibility.csv", targets, result)
    figure(out_dir / "fig_visibility", targets, result, airmass)
    return fragments


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--airmass", type=float, default=2.0,
                    help="airmass limit (default 2.0)")
    ap.add_argument("--inject", type=Path, default=None,
                    help="Markdown file whose GENERATED blocks to rewrite")
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR,
                    help="output directory (default ops/generated); the "
                         "eruption block keeps its airmass < 3 tables in "
                         "ops/eruption_block/generated")
    args = ap.parse_args()
    fragments = build(airmass=args.airmass, out_dir=args.out_dir)
    print(fragments["visibility_table"])
    print()
    print(fragments["visibility_dates"])
    if args.inject:
        done = inject(args.inject, fragments)
        print(f"\ninjected into {args.inject}: {', '.join(done) or 'nothing'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
