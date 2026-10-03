#!/usr/bin/env python
"""TCRB-N1 — is the RLMT T CrB grism series new?  The spectra-per-month table.

WHAT THIS SCRIPT ANSWERS
------------------------
The T CrB strategy sells the 2025 slitless-grism series as "the densest
homogeneous single-instrument H-alpha record" of the post-dip recovery, and
the journal editor (memo ``journal-editor.md``, TCRB-N1) ruled that the claim
must be *counted* before any pipeline work is spent on it:

    "spectra per month, RLMT vs ARAS vs Asiago/Munari, Feb-Jun 2025 ...
     if ARAS alone has denser, higher-resolution H-alpha coverage, the paper
     re-scopes to a validation/methods note (PASP) and says so."

This script does the counting.  It compares, on one calendar:

* **RLMT** — T CrB ``hrg`` / ``lrg`` frames (and the 2023-24 slot ``6`` / ``W``
  frames) read from the shared manifest, opened READ-ONLY;
* **ARAS** — the public ARAS Spectral Database listing for T CrB, one row per
  archived spectrum with its observer, site, resolving power and wavelength
  range (https://aras-database.github.io/database/tcrb.html).

Asiago/Munari, Tautenburg, Rozhen and the other professional campaigns do not
publish a per-spectrum log, so they cannot be counted by a script.  What they
*state* about their cadence is curated, with a URL per row, in
``external/literature_campaigns.csv`` and is carried into the generated
markdown verbatim — nothing about them is typed into this file.  The same
holds for the published H-alpha equivalent widths
(``external/published_halpha_points.csv``) and the three orbital ephemerides
(``external/ephemerides.csv``) used to place the RLMT window in orbital
phase.

EVERY NUMBER IS EMITTED, NONE IS TYPED
--------------------------------------
The only literals below are *definitions* (what counts as "covers H-alpha",
where the resolution classes break, the window) — each is a named constant
with its reason beside it, and each is echoed into the outputs.  All counts,
gaps and verdicts in ``NOVELTY_TABLE.md`` and the CSVs come from the two
databases.

THE THREE BIASES THIS COUNT COULD HIDE, AND WHAT IS DONE ABOUT EACH
-------------------------------------------------------------------
1. **Counting frames flatters RLMT.**  RLMT takes 1-4 back-to-back 240 s
   frames a night; ARAS observers upload one co-added spectrum per night.
   Frames are not independent epochs.  The unit of comparison is therefore
   the **UT calendar date with at least one H-alpha-covering spectrum**;
   frame/spectrum counts are reported beside it but never compared.
2. **Counting header labels flatters RLMT.**  A FILTER card that says ``hrg``
   does not prove a spectrum was recorded, or that the telescope was on
   T CrB (21 frames carry header pointings > 1 deg off).  RLMT is therefore
   counted in tiers: *labelled* (every canonical frame), and *pixel-verified
   dispersed* (S2c verdict ``dispersed`` in ``frame_dispersion``).  The
   identity gate (G-3) has not run on the full sample, so even the second
   tier is an UPPER bound on usable epochs; the table says so.
3. **Pooling ARAS flatters ARAS on density and RLMT on homogeneity.**  ARAS
   is ~20 observers with different spectrographs.  The fair test of a
   "homogeneous single-instrument" claim is against the best *single*
   ARAS series, so the script also ranks every (observer, site, resolution
   class) series inside the window.

The coverage definition is itself a choice, so the script recomputes the
headline under a stricter window (:data:`HA_STRICT`) and reports both.

OUTPUTS (all beside this file)
------------------------------
    external/aras_tcrb_listing.csv   parsed ARAS listing (one row / spectrum)
    external/aras_tcrb_pull.json     URL, pull time, sha256, row count
    rlmt_grism_dates.csv             RLMT UT dates x filter x tier
    novelty_monthly.csv              the spectra-per-month table, 2023-2026
    novelty_window.csv               the Feb-Jun 2025 window head-to-head
    aras_series_window.csv           every single-observer, single-class series
    aras_observers_window.csv        every ARAS observer, setups pooled
    novelty_verdict.json             the rule's clauses and the window facts
    NOVELTY_TABLE.md                 the tables above, rendered, with verdict
    fig_novelty_monthly.{png,pdf}    dates/month by source and resolution
    fig_novelty_timeline.{png,pdf}   who observed on which night, Feb-Jun 2025

USAGE
-----
    PY=/opt/miniconda3/envs/rlmt-checks/bin/python
    cd "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO"
    $PY TCrB_Monitoring/notes/novelty/build_novelty_table.py            # from cache
    $PY TCrB_Monitoring/notes/novelty/build_novelty_table.py --refresh  # re-pull ARAS

``--refresh`` performs one HTTP GET of a public page.  Nothing is written to
the manifest; nothing outside this directory is touched.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html as htmllib
import json
import re
import sqlite3
import sys
import urllib.request
from collections import defaultdict
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]                       # .../MACRO
sys.path.insert(0, str(REPO / "pipeline"))   # macro_core is not installed

from macro_core import plotstyle as ps       # noqa: E402  (house figure style)
import matplotlib.dates as mdates            # noqa: E402
import matplotlib.pyplot as plt              # noqa: E402  (after plotstyle sets Agg)

# ===========================================================================
# Definitions — the only literals in this file.  Each is echoed to the outputs.
# ===========================================================================

MANIFEST = REPO / "products" / "manifest" / "rlmt-manifest.sqlite"
EXT = HERE / "external"
ARAS_URL = "https://aras-database.github.io/database/tcrb.html"
ARAS_RAW = EXT / "aras_tcrb.html.gz"         # the page exactly as served
ARAS_CSV = EXT / "aras_tcrb_listing.csv"
ARAS_META = EXT / "aras_tcrb_pull.json"
LIT_CSV = EXT / "literature_campaigns.csv"   # curated, one URL per row
PUB_CSV = EXT / "published_halpha_points.csv"
EPH_CSV = EXT / "ephemerides.csv"            # curated, one URL per row

TARGET = "T CrB"

#: A spectrum "covers H-alpha" if it reaches both sides of 6563 A far enough
#: to hold the line (FWZI ~ 700 km/s = 15 A, ATel 17041) and a continuum
#: window on each side.  LOOSE is the headline; STRICT demands ~130 A a side,
#: enough for the TiO-free pseudo-continuum bands a late-M giant needs, and
#: is the sensitivity check.  Narrow echelle H-alpha windows (e.g. 6500-6700)
#: pass LOOSE and fail STRICT, which is exactly the population the check
#: is there to expose.
HA_LOOSE = (6500.0, 6620.0)
HA_STRICT = (6430.0, 6700.0)

#: Resolution classes.  ``low`` is the regime RLMT itself occupies (the D1
#: test will put the delivered RLMT resolving power somewhere between ~200
#: and ~1000, physicist.md P8 vs observational-astronomer.md E1); ``mid`` can
#: separate the H-alpha core from its wings; ``high`` resolves the
#: double-peaked disc profile (peak separation ~90 km/s, ATel 17075, needs
#: R >~ 7000).
R_CLASSES = (("low", 0.0, 2000.0), ("mid", 2000.0, 7000.0),
             ("high", 7000.0, np.inf))
#: RLMT's delivered resolving power is NOT known: it is disagreement D1 of
#: the chair's synthesis.  Both seats agree the T CrB H-alpha line is at
#: best 13 px wide on hrg; they disagree on the dispersion that converts
#: pixels to Angstrom.  The "higher-resolution" clause of the editor's rule
#: is therefore evaluated on BOTH branches, each at the narrowest (most
#: favourable) measured line width, and never on a single adopted value.
RLMT_HA_WIDTH_PX = 13.0                      # narrowest hrg H-alpha width, px
RLMT_R_BRANCH = {
    "D1 -> 0.47 A/px (observational-astronomer.md E1)":
        6563.0 / (0.47 * RLMT_HA_WIDTH_PX),
    "D1 -> 1.59 A/px (physicist.md P8, code value)":
        6563.0 / (1.59 * RLMT_HA_WIDTH_PX),
}

#: RLMT filter labels that are, or may be, grisms on T CrB.
RLMT_GRISM_2025 = ("hrg", "lrg")
RLMT_GRISM_EARLY = ("6", "W")                # TCRB-A0b: 2023-05 -> 2024-03

#: The head-to-head window is not typed: it is the first and last UT date of
#: the RLMT 2025 grism series, read from the manifest.  The monthly table
#: spans the whole context interval below.
CONTEXT = (date(2023, 1, 1), date(2026, 9, 30))

#: The editor's months (journal-editor.md, TCRB-N1).
EDITOR_MONTHS = ("2025-02", "2025-03", "2025-04", "2025-05", "2025-06")


# ===========================================================================
# Small pure helpers (unit-tested in test_novelty_table.py)
# ===========================================================================
def jd_to_utc_date(jd: float) -> date:
    """UT calendar date of a Julian Date.

    Both databases are put on this one clock.  RLMT's ``night`` column is the
    LOCAL evening date (Arizona, UT-7), one day behind the UT date of most
    of its frames; ARAS lists UT.  Binning one by local night and the other
    by UT date would shift RLMT epochs across month boundaries.
    """
    unix = (float(jd) - 2440587.5) * 86400.0
    return (datetime(1970, 1, 1, tzinfo=timezone.utc)
            + timedelta(seconds=unix)).date()


def covers(lmin: float, lmax: float, window: tuple[float, float]) -> bool:
    """True if the spectrum spans the whole of ``window`` (Angstrom)."""
    return lmin <= window[0] and lmax >= window[1]


def r_class(resolving_power: float) -> str:
    """Name of the resolution class a resolving power falls in."""
    for name, lo, hi in R_CLASSES:
        if lo <= resolving_power < hi:
            return name
    raise ValueError(f"resolving power {resolving_power!r} fits no class")


def gap_stats(dates) -> dict:
    """Cadence of a set of dates: count, median gap, longest gap (days).

    A density claim is a claim about gaps, not about counts: 60 dates in two
    clumps are not a 2-day cadence.  ``max_gap`` is what a referee will ask
    for; the median alone hides it.
    """
    d = sorted(set(dates))
    if len(d) < 2:
        return {"n_dates": len(d), "median_gap_d": float("nan"),
                "max_gap_d": float("nan")}
    gaps = np.diff([x.toordinal() for x in d])
    return {"n_dates": len(d), "median_gap_d": float(np.median(gaps)),
            "max_gap_d": float(gaps.max())}


def orbital_phase(jd: float, period: float, t0: float) -> float:
    """Orbital phase in [0, 1) of ``jd`` on a linear ephemeris.

    All three published ephemerides in ``external/ephemerides.csv`` share one
    convention: phase 0 is the red giant at maximum recession velocity
    (ascending quadrature).  The giant is then seen side-on at phases 0 and
    0.5 (ellipsoidal MAXIMA) and end-on at 0.25 and 0.75 (ellipsoidal
    MINIMA; 0.75 is the giant in front of the white dwarf).
    """
    return ((float(jd) - float(t0)) / float(period)) % 1.0


def date_to_jd(d: date) -> float:
    """JD at 00:00 UT of a calendar date."""
    return d.toordinal() + 1721424.5


def phase_table(facts: dict) -> list[dict]:
    """Orbital phase of the RLMT window under each published ephemeris.

    Why this is in a novelty script: the physicist's amendment TCRB-A5b
    rests on the window spanning about one ellipsoidal cycle, so that raw EW
    carries a guaranteed geometric oscillation.  That is arithmetic on a
    literature ephemeris and on the manifest's window, so it is emitted here
    rather than typed into a report — and done on three ephemerides so the
    answer is seen not to depend on which one is adopted.
    """
    out = []
    lo, hi = date_to_jd(facts["window_lo"]), date_to_jd(facts["window_hi"]) + 1.0
    for e in read_csv(EPH_CSV):
        per, t0 = float(e["period_d"]), float(e["t0_hjd"])
        perr = float(e["period_err_d"]) if e["period_err_d"] else float("nan")
        n_cyc = (lo - t0) / per
        out.append({"ephemeris": e["ephemeris"], "period_d": per, "t0_hjd": t0,
                    "phase_start": orbital_phase(lo, per, t0),
                    "phase_end": orbital_phase(hi, per, t0),
                    "ellipsoidal_cycles_spanned": (hi - lo) / (per / 2.0),
                    # how far the phase has drifted since T0 through the
                    # period uncertainty alone (heliocentric correction,
                    # < 0.006 d, is negligible beside it)
                    "phase_err_from_period": abs(n_cyc) * perr / per,
                    "source": e["source"]})
    return out


def month_key(d: date) -> str:
    return f"{d.year:04d}-{d.month:02d}"


def month_range(lo: date, hi: date) -> list[str]:
    """Every 'YYYY-MM' from ``lo`` to ``hi`` inclusive — empty months too.

    A month with no spectra is a row of zeros, not a missing row: the
    solar-conjunction gaps are part of the answer.
    """
    out, y, m = [], lo.year, lo.month
    while (y, m) <= (hi.year, hi.month):
        out.append(f"{y:04d}-{m:02d}")
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


# ===========================================================================
# ARAS: pull, parse, cache
# ===========================================================================
def parse_aras(page: str) -> list[dict]:
    """One dict per spectrum row of the ARAS T CrB page.

    The page is a single HTML table whose data rows have thirteen cells:
    date, UT time, JD-2400000, observer code, site code, resolving power,
    (a colour swatch), lambda_min, lambda_max, (a swatch), file link,
    preview link, comments.  Rows that do not have thirteen cells are the
    header; anything else that fails to parse raises, because a silently
    dropped row is a silently wrong count.
    """
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", page, flags=re.S):
        td = re.findall(r"<td>(.*?)</td>", tr, flags=re.S)
        if not td:
            continue                          # header row: <th> only
        if len(td) != 13:
            raise ValueError(f"ARAS row with {len(td)} cells: layout changed")
        link = re.search(r'href="(spectra/[^"]+)"', tr)
        rows.append({
            "date_ut": td[0].strip(),
            "time_ut": td[1].strip(),
            "jd": 2400000.0 + float(td[2]),
            # Codes are upper-cased: the listing carries the same observer
            # as e.g. "XDU" and "XDu", and a case variant must not be
            # counted as a second instrument in the homogeneity ranking.
            "observer": htmllib.unescape(td[3].strip()).upper(),
            "site": htmllib.unescape(td[4].strip()).upper(),
            "resolving_power": float(td[5]),
            "lambda_min": float(td[7]),
            "lambda_max": float(td[8]),
            "file": link.group(1) if link else "",
        })
    return rows


def declared_count(page: str) -> int | None:
    """The 'Number of spectra' the page itself advertises, if present."""
    m = re.search(r"Number of spectra:\s*(\d+)", page)
    return int(m.group(1)) if m else None


def refresh_aras() -> None:
    """GET the ARAS page and cache it, gzipped, exactly as served."""
    EXT.mkdir(exist_ok=True)
    req = urllib.request.Request(ARAS_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        raw = resp.read()
    with gzip.open(ARAS_RAW, "wb", compresslevel=9) as fh:
        fh.write(raw)
    meta = {"url": ARAS_URL,
            "pulled_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256_html": hashlib.sha256(raw).hexdigest(),
            "bytes_html": len(raw)}
    ARAS_META.write_text(json.dumps(meta, indent=2) + "\n")


def load_aras() -> tuple[list[dict], dict]:
    """Parse the cached page; verify it against its own declared row count."""
    if not ARAS_RAW.exists():
        raise SystemExit(f"no cached ARAS page at {ARAS_RAW}; run with --refresh")
    raw = gzip.open(ARAS_RAW, "rb").read()
    page = raw.decode("utf-8", errors="replace")
    rows = parse_aras(page)
    declared = declared_count(page)
    # A truncated download parses cleanly and undercounts.  The page states
    # how many spectra it lists; refuse to proceed if we did not get them all.
    if declared is not None and declared != len(rows):
        raise SystemExit(f"ARAS page declares {declared} spectra, parsed "
                         f"{len(rows)} — truncated or layout changed")
    meta = json.loads(ARAS_META.read_text()) if ARAS_META.exists() else {}
    meta.update({"url": ARAS_URL, "n_rows": len(rows),
                 "declared_n": declared,
                 "sha256_html": hashlib.sha256(raw).hexdigest(),
                 "first_date": min(r["date_ut"] for r in rows),
                 "last_date": max(r["date_ut"] for r in rows)})
    ARAS_META.write_text(json.dumps(meta, indent=2) + "\n")
    for r in rows:
        r["date"] = jd_to_utc_date(r["jd"])
        r["r_class"] = r_class(r["resolving_power"])
        r["ha_loose"] = covers(r["lambda_min"], r["lambda_max"], HA_LOOSE)
        r["ha_strict"] = covers(r["lambda_min"], r["lambda_max"], HA_STRICT)
    with ARAS_CSV.open("w", newline="") as fh:
        cols = ["date_ut", "time_ut", "jd", "observer", "site",
                "resolving_power", "lambda_min", "lambda_max", "r_class",
                "ha_loose", "ha_strict", "file"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["jd"]))
    return rows, meta


# ===========================================================================
# RLMT: read the manifest, read-only
# ===========================================================================
RLMT_SQL = """
SELECT f.obs_rowid, f.filter, f.night, f.jd, f.exptime,
       f.pointing_offset_deg, d.verdict AS s2c_verdict
FROM frames AS f
LEFT JOIN frame_dispersion AS d ON d.obs_rowid = f.obs_rowid
WHERE f.canonical_target = ?
  AND f.is_canonical = 1
  AND f.filter IN ({marks})
ORDER BY f.jd
"""


def load_rlmt() -> list[dict]:
    """Every canonical T CrB grism-labelled frame, with its S2c verdict.

    ``is_canonical = 1`` removes the reduced/external twins of the same
    exposure (the F-1 dedup); without it the 2025 series double-counts.
    The connection is opened ``mode=ro`` — the manifest is shared and this
    package owns no table in it.
    """
    filters = RLMT_GRISM_2025 + RLMT_GRISM_EARLY
    sql = RLMT_SQL.format(marks=",".join("?" * len(filters)))
    con = sqlite3.connect(f"file:{MANIFEST}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    try:
        rows = [dict(r) for r in con.execute(sql, (TARGET, *filters))]
    finally:
        con.close()
    for r in rows:
        r["date"] = jd_to_utc_date(r["jd"])
        r["dispersed"] = r["s2c_verdict"] == "dispersed"
        off = r["pointing_offset_deg"]
        r["header_off"] = off is not None and off > 1.0
    return rows


# ===========================================================================
# The tables
# ===========================================================================
def monthly_table(aras: list[dict], rlmt: list[dict]) -> list[dict]:
    """Spectra and dates per UT month, both sources, CONTEXT interval."""
    rows = []
    for mk in month_range(*CONTEXT):
        a = [r for r in aras if month_key(r["date"]) == mk and r["ha_loose"]]
        a_all = [r for r in aras if month_key(r["date"]) == mk]
        g = [r for r in rlmt if month_key(r["date"]) == mk]
        gd = [r for r in g if r["dispersed"]]
        rp = [r["resolving_power"] for r in a]
        row = {"month": mk,
               "aras_spectra_all": len(a_all),
               "aras_spectra_ha": len(a),
               "aras_dates_ha": len({r["date"] for r in a}),
               "aras_observers_ha": len({(r["observer"], r["site"]) for r in a}),
               "aras_R_median": float(np.median(rp)) if rp else float("nan")}
        for name, _, _ in R_CLASSES:
            sub = [r for r in a if r["r_class"] == name]
            row[f"aras_spectra_{name}R"] = len(sub)
            row[f"aras_dates_{name}R"] = len({r["date"] for r in sub})
        row["aras_dates_ha_strict"] = len(
            {r["date"] for r in a_all if r["ha_strict"]})
        row["rlmt_frames_labelled"] = len(g)
        row["rlmt_dates_labelled"] = len({r["date"] for r in g})
        row["rlmt_frames_dispersed"] = len(gd)
        row["rlmt_dates_dispersed"] = len({r["date"] for r in gd})
        for filt in RLMT_GRISM_2025 + RLMT_GRISM_EARLY:
            row[f"rlmt_dates_dispersed_{filt}"] = len(
                {r["date"] for r in gd if r["filter"] == filt})
        rows.append(row)
    return rows


def aras_series(aras_window: list[dict]) -> list[dict]:
    """Every single-observer, single-resolution-class ARAS series, ranked.

    This is the homogeneity test.  A "series" is one observer code at one
    site within one resolution class — the nearest proxy the listing offers
    for "one instrument in one configuration".  It is a proxy: an observer
    can change a grating within a class, so a series here is AT LEAST as
    homogeneous as stated only if its resolving-power spread is small, which
    is why min/median/max R are carried beside the count.
    """
    groups = defaultdict(list)
    for r in aras_window:
        groups[(r["observer"], r["site"], r["r_class"])].append(r)
    out = []
    for (obs, site, cls), rs in groups.items():
        rp = np.array([r["resolving_power"] for r in rs])
        out.append({"observer": obs, "site": site, "r_class": cls,
                    "n_spectra": len(rs),
                    **gap_stats(r["date"] for r in rs),
                    "R_min": float(rp.min()), "R_median": float(np.median(rp)),
                    "R_max": float(rp.max()),
                    "lambda_min": min(r["lambda_min"] for r in rs),
                    "lambda_max": max(r["lambda_max"] for r in rs)})
    out.sort(key=lambda s: (-s["n_dates"], -s["n_spectra"]))
    return out


def aras_observers(aras_window: list[dict]) -> list[dict]:
    """Every ARAS observer/site in the window, pooling that observer's setups.

    The weaker homogeneity test: one person at one telescope, but possibly
    alternating between a low-resolution and an echelle spectrograph.  An
    observer with many dates here and few in :func:`aras_series` is two
    instruments, and is reported as such (``n_classes``).
    """
    groups = defaultdict(list)
    for r in aras_window:
        groups[(r["observer"], r["site"])].append(r)
    out = []
    for (obs, site), rs in groups.items():
        out.append({"observer": obs, "site": site,
                    "n_classes": len({r["r_class"] for r in rs}),
                    "n_spectra": len(rs),
                    **gap_stats(r["date"] for r in rs)})
    out.sort(key=lambda s: (-s["n_dates"], -s["n_spectra"]))
    return out


def window_table(aras: list[dict], rlmt: list[dict]) -> tuple[list[dict], dict]:
    """Head-to-head inside the RLMT 2025 grism window.

    The window is RLMT's own first-to-last UT date, so RLMT is compared on
    the interval most favourable to it: ARAS coverage outside the window
    (which exists on both sides) is not counted against it here.
    """
    g25 = [r for r in rlmt if r["filter"] in RLMT_GRISM_2025]
    lo, hi = min(r["date"] for r in g25), max(r["date"] for r in g25)
    span = (hi - lo).days + 1
    inwin = lambda r: lo <= r["date"] <= hi                     # noqa: E731
    aw = [r for r in aras if inwin(r)]
    aw_ha = [r for r in aw if r["ha_loose"]]
    series = aras_series(aw_ha)
    observers = aras_observers(aw_ha)

    def line(label, rows, n_items, note):
        st = gap_stats(r["date"] for r in rows)
        rp = [r["resolving_power"] for r in rows if "resolving_power" in r]
        return {"series": label, "n_spectra_or_frames": n_items, **st,
                "date_fill_fraction": st["n_dates"] / span,
                "R_median": float(np.median(rp)) if rp else float("nan"),
                "note": note}

    rows = []
    for filt in RLMT_GRISM_2025:
        sub = [r for r in g25 if r["filter"] == filt]
        rows.append(line(f"RLMT {filt}, labelled", sub, len(sub),
                         "FILTER card only; identity ungated"))
        sub_d = [r for r in sub if r["dispersed"]]
        rows.append(line(f"RLMT {filt}, S2c dispersed", sub_d, len(sub_d),
                         "pixel-verified spectrum; identity ungated"))
    rows.append(line("RLMT hrg+lrg, labelled", g25, len(g25),
                     "upper bound on epochs"))
    g25d = [r for r in g25 if r["dispersed"]]
    rows.append(line("RLMT hrg+lrg, S2c dispersed", g25d, len(g25d),
                     "upper bound until G-3 identity gate runs"))
    rows.append(line("ARAS, all H-alpha (loose)", aw_ha, len(aw_ha),
                     f"{len({(r['observer'], r['site']) for r in aw_ha})} observers pooled"))
    aw_s = [r for r in aw if r["ha_strict"]]
    rows.append(line("ARAS, all H-alpha (strict)", aw_s, len(aw_s),
                     "sensitivity: wide continuum windows required"))
    for name, _, _ in R_CLASSES:
        sub = [r for r in aw_ha if r["r_class"] == name]
        rows.append(line(f"ARAS, {name}-R only", sub, len(sub), "pooled"))
    for s in series[:3]:
        sub = [r for r in aw_ha if (r["observer"], r["site"], r["r_class"])
               == (s["observer"], s["site"], s["r_class"])]
        rows.append(line(f"ARAS single series {s['observer']}/{s['site']} "
                         f"({s['r_class']}-R)", sub, len(sub),
                         "one observer, one site, one resolution class"))
    if observers:
        o = observers[0]
        sub = [r for r in aw_ha
               if (r["observer"], r["site"]) == (o["observer"], o["site"])]
        rows.append(line(f"ARAS single observer {o['observer']}/{o['site']} "
                         f"(all setups)", sub, len(sub),
                         f"one observer, {o['n_classes']} resolution classes"))

    # Dates on which the two archives overlap: the cross-validation sample
    # TCRB-A5 will actually have (same UT date; +-1 d is the looser match).
    rl_dates = {r["date"] for r in g25d}
    ar_dates = {r["date"] for r in aw_ha}
    near = {d for d in rl_dates
            if any((d + timedelta(days=k)) in ar_dates for k in (-1, 0, 1))}
    facts = {"window_lo": lo, "window_hi": hi, "window_days": span,
             "rlmt_header_off_frames": sum(r["header_off"] for r in g25),
             "rlmt_frames_labelled": len(g25),
             "rlmt_frames_dispersed": len(g25d),
             "rlmt_dates_labelled": len({r["date"] for r in g25}),
             "rlmt_dates_dispersed": len(rl_dates),
             "rlmt_local_nights_labelled": len({r["night"] for r in g25}),
             "rlmt_both_grisms_dates": len(
                 {r["date"] for r in g25d if r["filter"] == "hrg"}
                 & {r["date"] for r in g25d if r["filter"] == "lrg"}),
             "aras_spectra_window_all": len(aw),
             "aras_spectra_ha": len(aw_ha),
             "aras_dates_ha": len(ar_dates),
             "aras_observers_ha": len({(r["observer"], r["site"]) for r in aw_ha}),
             "aras_R_median": float(np.median([r["resolving_power"] for r in aw_ha])),
             "aras_vs_rlmt_R": [
                 {"branch": name, "rlmt_R": rl_r,
                  "aras_frac_spectra_above": float(np.mean(
                      [r["resolving_power"] > rl_r for r in aw_ha])),
                  "aras_dates_above": len({r["date"] for r in aw_ha
                                           if r["resolving_power"] > rl_r}),
                  "aras_median_above": bool(np.median(
                      [r["resolving_power"] for r in aw_ha]) > rl_r)}
                 for name, rl_r in RLMT_R_BRANCH.items()],
             "aras_dates_highR": len({r["date"] for r in aw_ha
                                      if r["r_class"] == "high"}),
             "overlap_same_date": len(rl_dates & ar_dates),
             "overlap_within_1d": len(near),
             "rlmt_dates_without_aras": len(rl_dates - ar_dates),
             "aras_dates_without_rlmt": len(ar_dates - rl_dates),
             "best_series": series[0] if series else None,
             "best_observer": observers[0] if observers else None,
             "series": series, "observers": observers}
    return rows, facts


def verdict(monthly: list[dict], facts: dict) -> dict:
    """Evaluate the editor's re-scope rule, clause by clause, from the counts.

    The rule (journal-editor.md, TCRB-N1): *if ARAS alone has denser,
    higher-resolution H-alpha coverage, the paper re-scopes to a
    validation/methods note (PASP)*.  Two clauses, both must hold.

    **Denser** is judged against RLMT's LABELLED dates — the most generous
    RLMT count — over the whole window and in every editor month, so the
    clause can only be passed against RLMT's best case.

    **Higher-resolution** cannot be judged against one RLMT number because
    D1 is open, so it is judged per D1 branch and in two ways: the pooled
    ARAS median (a weak statistic — ARAS is bimodal, a grism population
    near R ~ 500-1000 and an echelle population near R ~ 13 000, and many
    entries are a round 1000), and the more telling count of ARAS dates
    that individually exceed the RLMT value.  The rule is reported as
    firing on a branch when ARAS has MORE dates above the RLMT resolving
    power than half of RLMT's own labelled dates AND is denser overall:
    i.e. when the better-resolved ARAS subset is by itself a series of
    comparable density.  That threshold is a judgement, so the raw counts
    are printed beside it for a reader who would draw the line elsewhere.

    The last lines evaluate the strategy's fallback wording — "homogeneous
    single-instrument" (strategy A.8) — against the best single ARAS series
    and the best single ARAS observer.
    """
    ed = [m for m in monthly if m["month"] in EDITOR_MONTHS]
    months_aras_denser = [m["month"] for m in ed
                          if m["aras_dates_ha"] > m["rlmt_dates_dispersed"]]
    months_aras_denser_lab = [m["month"] for m in ed
                              if m["aras_dates_ha"] > m["rlmt_dates_labelled"]]
    denser = (facts["aras_dates_ha"] > facts["rlmt_dates_labelled"]
              and len(months_aras_denser_lab) == len(ed))
    branches = []
    for b in facts["aras_vs_rlmt_R"]:
        comparable = b["aras_dates_above"] > 0.5 * facts["rlmt_dates_labelled"]
        branches.append({**b, "aras_better_resolved_series_comparable": comparable,
                         "rule_fires": bool(denser and comparable)})
    best, best_obs = facts["best_series"], facts["best_observer"]
    return {"aras_denser": denser,
            "months_aras_denser_vs_dispersed": months_aras_denser,
            "months_aras_denser_vs_labelled": months_aras_denser_lab,
            "branches": branches,
            "rule_fires_on_every_branch": all(b["rule_fires"] for b in branches),
            "rule_fires_on_any_branch": any(b["rule_fires"] for b in branches),
            "best_single_series_dates": best["n_dates"] if best else 0,
            "best_single_observer_dates": best_obs["n_dates"] if best_obs else 0,
            "rlmt_margin_over_best_single_series":
                facts["rlmt_dates_dispersed"] - (best["n_dates"] if best else 0),
            "rlmt_margin_over_best_single_observer":
                facts["rlmt_dates_dispersed"]
                - (best_obs["n_dates"] if best_obs else 0)}


# ===========================================================================
# Writers
# ===========================================================================
def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    cols = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v)
                        for k, v in r.items()})


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def fmt(v) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if np.isnan(v):
            return "—"
        return f"{v:.0f}" if abs(v) >= 100 else f"{v:.2g}"
    return str(v)


def md_table(rows: list[dict], cols: list[tuple[str, str]]) -> str:
    head = "| " + " | ".join(h for _, h in cols) + " |"
    rule = "|" + "|".join("---" for _ in cols) + "|"
    body = ["| " + " | ".join(fmt(r[k]) for k, _ in cols) + " |" for r in rows]
    return "\n".join([head, rule, *body])


def write_markdown(monthly, window_rows, facts, v, meta, rlmt) -> None:
    """Render every table to NOVELTY_TABLE.md.  No number here is typed."""
    ed = [m for m in monthly if m["month"] in EDITOR_MONTHS]
    busy = [m for m in monthly
            if m["aras_spectra_all"] or m["rlmt_frames_labelled"]]
    early = [r for r in rlmt if r["filter"] in RLMT_GRISM_EARLY]
    early_d = [r for r in early if r["dispersed"]]
    best = facts["best_series"]
    L = []
    L.append("# TCRB-N1 — novelty table (generated; do not edit)\n")
    L.append(f"Generated by `build_novelty_table.py`. ARAS listing: {meta['url']} "
             f"pulled {meta.get('pulled_utc', 'unknown')} "
             f"({meta['n_rows']} spectra, page declares {meta['declared_n']}, "
             f"{meta['first_date']} → {meta['last_date']}, "
             f"sha256 `{meta['sha256_html'][:16]}…`). "
             f"RLMT: `products/manifest/rlmt-manifest.sqlite` opened read-only.\n")
    L.append("Definitions: *covers Hα* = λmin ≤ "
             f"{HA_LOOSE[0]:.0f} Å and λmax ≥ {HA_LOOSE[1]:.0f} Å "
             f"(strict check: {HA_STRICT[0]:.0f}–{HA_STRICT[1]:.0f} Å). "
             "Resolution classes: "
             + ", ".join(f"{n} = {lo:.0f} ≤ R < {hi:.0f}" if np.isfinite(hi)
                         else f"{n} = R ≥ {lo:.0f}" for n, lo, hi in R_CLASSES)
             + ". A *date* is a UT calendar date with ≥ 1 qualifying spectrum. "
             "RLMT *dispersed* = S2c pixel verdict `dispersed`; the identity "
             "gate (G-3) has not run, so RLMT dates are upper bounds.\n")

    L.append("## 1. The editor's table — Feb–Jun 2025, per UT month\n")
    L.append(md_table(ed, [
        ("month", "month"),
        ("rlmt_frames_labelled", "RLMT frames (labelled)"),
        ("rlmt_dates_labelled", "RLMT dates (labelled)"),
        ("rlmt_dates_dispersed", "RLMT dates (S2c dispersed)"),
        ("rlmt_dates_dispersed_lrg", "…lrg"),
        ("rlmt_dates_dispersed_hrg", "…hrg"),
        ("aras_spectra_ha", "ARAS Hα spectra"),
        ("aras_dates_ha", "ARAS dates"),
        ("aras_observers_ha", "ARAS observers"),
        ("aras_R_median", "ARAS median R"),
        ("aras_dates_lowR", "ARAS dates low-R"),
        ("aras_dates_midR", "ARAS dates mid-R"),
        ("aras_dates_highR", "ARAS dates high-R"),
        ("aras_dates_ha_strict", "ARAS dates (strict Hα)")]))
    L.append("")

    L.append(f"## 2. Head-to-head inside the RLMT window "
             f"({facts['window_lo']} → {facts['window_hi']} UT, "
             f"{facts['window_days']} d)\n")
    L.append(md_table(window_rows, [
        ("series", "series"), ("n_spectra_or_frames", "spectra / frames"),
        ("n_dates", "dates"), ("date_fill_fraction", "fraction of days"),
        ("median_gap_d", "median gap (d)"), ("max_gap_d", "max gap (d)"),
        ("R_median", "median R"), ("note", "note")]))
    L.append("")
    L.append(f"- RLMT frames with header pointing > 1° off target: "
             f"{facts['rlmt_header_off_frames']} of {facts['rlmt_frames_labelled']} "
             f"(identity must be established from pixels, G-3).")
    L.append(f"- RLMT local nights (labelled): {facts['rlmt_local_nights_labelled']}; "
             f"UT dates (labelled): {facts['rlmt_dates_labelled']}; "
             f"UT dates with both grisms dispersed: {facts['rlmt_both_grisms_dates']}.")
    L.append(f"- ARAS in the same window: {facts['aras_spectra_ha']} Hα-covering "
             f"spectra of {facts['aras_spectra_window_all']} archived, on "
             f"{facts['aras_dates_ha']} dates, by {facts['aras_observers_ha']} "
             f"observer/site pairs; median R = {facts['aras_R_median']:.0f}; "
             f"high-R (profile-resolving) spectra on "
             f"{facts['aras_dates_highR']} dates.")
    for b in facts["aras_vs_rlmt_R"]:
        L.append(f"- RLMT resolving power if {b['branch']}: R ≈ {b['rlmt_R']:.0f} "
                 f"at the narrowest measured Hα width ({RLMT_HA_WIDTH_PX:.0f} px). "
                 f"ARAS spectra above that: "
                 f"{100 * b['aras_frac_spectra_above']:.0f}% of spectra, on "
                 f"{b['aras_dates_above']} dates.")
    L.append(f"- Cross-validation sample for TCRB-A5: {facts['overlap_same_date']} "
             f"UT dates have both an RLMT dispersed frame and an ARAS Hα spectrum; "
             f"{facts['overlap_within_1d']} RLMT dates have an ARAS spectrum within "
             f"±1 d. RLMT dates with no same-date ARAS spectrum: "
             f"{facts['rlmt_dates_without_aras']}; ARAS dates with no RLMT frame: "
             f"{facts['aras_dates_without_rlmt']}.\n")

    L.append("## 3. Homogeneity — every single-observer ARAS series in the window\n")
    L.append(md_table(facts["series"], [
        ("observer", "observer"), ("site", "site"), ("r_class", "class"),
        ("n_spectra", "spectra"), ("n_dates", "dates"),
        ("median_gap_d", "median gap (d)"), ("max_gap_d", "max gap (d)"),
        ("R_min", "R min"), ("R_median", "R median"), ("R_max", "R max"),
        ("lambda_min", "λmin (Å)"), ("lambda_max", "λmax (Å)")]))
    L.append("")

    L.append("## 4. The re-scope rule, evaluated\n")
    L.append("| clause | value |\n|---|---|")
    L.append(f"| ARAS dates > RLMT *labelled* dates, whole window and every "
             f"editor month (**denser**) | {fmt(v['aras_denser'])} |")
    L.append(f"| months where ARAS dates > RLMT labelled dates | "
             f"{', '.join(v['months_aras_denser_vs_labelled']) or 'none'} |")
    L.append(f"| months where ARAS dates > RLMT S2c-dispersed dates | "
             f"{', '.join(v['months_aras_denser_vs_dispersed']) or 'none'} |")
    for b in v["branches"]:
        L.append(f"| *{b['branch']}*: RLMT R ≈ {b['rlmt_R']:.0f} | |")
        L.append(f"| … pooled ARAS median R exceeds it | "
                 f"{fmt(b['aras_median_above'])} |")
        L.append(f"| … ARAS dates individually above it | {b['aras_dates_above']} "
                 f"(RLMT labelled dates: {facts['rlmt_dates_labelled']}) |")
        L.append(f"| … better-resolved ARAS subset is itself a comparable series "
                 f"(> half of RLMT's dates) | "
                 f"{fmt(b['aras_better_resolved_series_comparable'])} |")
        L.append(f"| … **editor's rule fires on this branch** | "
                 f"**{fmt(b['rule_fires'])}** |")
    L.append(f"| **rule fires on every D1 branch** | "
             f"**{fmt(v['rule_fires_on_every_branch'])}** |")
    if best:
        L.append(f"| best single ARAS series (one observer, one class) | "
                 f"{best['observer']}/{best['site']} "
                 f"{best['r_class']}-R: {best['n_dates']} dates, median gap "
                 f"{fmt(best['median_gap_d'])} d, max gap {fmt(best['max_gap_d'])} d, "
                 f"R {best['R_min']:.0f}–{best['R_max']:.0f} |")
    bo = facts["best_observer"]
    if bo:
        L.append(f"| best single ARAS observer (all setups) | "
                 f"{bo['observer']}/{bo['site']}: {bo['n_dates']} dates in "
                 f"{bo['n_classes']} resolution classes, median gap "
                 f"{fmt(bo['median_gap_d'])} d, max gap {fmt(bo['max_gap_d'])} d |")
    L.append(f"| RLMT S2c-dispersed dates minus best single ARAS series | "
             f"{v['rlmt_margin_over_best_single_series']:+d} |")
    L.append(f"| RLMT S2c-dispersed dates minus best single ARAS observer | "
             f"{v['rlmt_margin_over_best_single_observer']:+d} |\n")

    L.append("## 5. RLMT 2023–24 slot '6' / 'W' frames (TCRB-A0b)\n")
    L.append(f"{len(early)} canonical frames labelled "
             f"{' / '.join(repr(f) for f in RLMT_GRISM_EARLY)} on "
             f"{len({r['date'] for r in early})} UT dates "
             f"({min(r['date'] for r in early)} → {max(r['date'] for r in early)}); "
             f"{len(early_d)} are S2c `dispersed`, on "
             f"{len({r['date'] for r in early_d})} dates. ARAS Hα-covering "
             "spectra in the same months are in the context table below.\n")

    L.append("## 6. Context — every month with data, "
             f"{CONTEXT[0]} → {CONTEXT[1]}\n")
    L.append(md_table(busy, [
        ("month", "month"), ("aras_spectra_ha", "ARAS Hα spectra"),
        ("aras_dates_ha", "ARAS dates"), ("aras_observers_ha", "observers"),
        ("aras_R_median", "median R"), ("aras_dates_highR", "high-R dates"),
        ("rlmt_frames_labelled", "RLMT frames"),
        ("rlmt_dates_labelled", "RLMT dates (labelled)"),
        ("rlmt_dates_dispersed", "RLMT dates (dispersed)")]))
    L.append("")

    ph = phase_table(facts)
    if ph:
        L.append("## 7. Orbital phase of the RLMT window "
                 "(phase 0 = red giant at maximum velocity; ellipsoidal minima "
                 "at 0.25 and 0.75)\n")
        L.append("| ephemeris | P (d) | T0 (HJD) | phase at window start | "
                 "phase at window end | ellipsoidal cycles spanned | "
                 "phase error from σ(P) | source |\n|---|---|---|---|---|---|---|---|")
        for r in ph:
            L.append(f"| {r['ephemeris']} | {r['period_d']:.4f} | {r['t0_hjd']:.2f} | "
                     f"{r['phase_start']:.3f} | {r['phase_end']:.3f} | "
                     f"{r['ellipsoidal_cycles_spanned']:.2f} | "
                     f"{fmt(r['phase_err_from_period'])} | {r['source']} |")
        L.append("")
    lit = read_csv(LIT_CSV)
    if lit:
        L.append("## 8. Campaigns that publish no per-spectrum log "
                 "(curated; one source per row)\n")
        L.append(md_table(lit, [
            ("campaign", "campaign"), ("instrument", "instrument"),
            ("resolving_power", "R"), ("coverage_stated", "coverage as stated"),
            ("halpha_product", "Hα product"), ("source", "source")]))
        L.append("")
    pub = read_csv(PUB_CSV)
    if pub:
        L.append("## 9. Published Hα equivalent widths of T CrB "
                 "(curated; one source per row)\n")
        L.append(md_table(pub, [
            ("date_ut", "UT date"), ("ew_abs_A", "EW, emission positive (Å)"),
            ("ew_err_A", "± (Å)"), ("instrument", "instrument"),
            ("source", "source")]))
        L.append("")
    (HERE / "NOVELTY_TABLE.md").write_text("\n".join(L) + "\n")


# ===========================================================================
# Figures — the plots do the talking
# ===========================================================================
def fig_monthly(monthly: list[dict], facts: dict) -> None:
    """Dates per month, ARAS stacked by resolution class, RLMT beside it.

    Bars are DATES, not spectra (bias 1 in the module docstring).  ARAS
    classes are stacked on 'dates in that class', which can exceed the
    pooled date count when two classes observe the same night; the pooled
    count is therefore drawn as its own marker so nothing is read off a
    stack that double-counts.
    """
    months = [m["month"] for m in monthly]
    x = np.arange(len(months))
    width = 0.42
    class_color = {"low": ps.SECOND, "mid": ps.WARN, "high": ps.BAD}
    class_hatch = {"low": "", "mid": "//", "high": "xx"}
    with ps.context("web"):
        fig, ax = plt.subplots(figsize=(ps.COL_DOUBLE * 1.35, 4.6))
        bottom = np.zeros(len(months))
        for name, lo, hi in R_CLASSES:
            y = np.array([m[f"aras_dates_{name}R"] for m in monthly], float)
            lab = (f"ARAS {name}-R ({lo:.0f} ≤ R < {hi:.0f})" if np.isfinite(hi)
                   else f"ARAS {name}-R (R ≥ {lo:.0f})")
            ax.bar(x - width / 2, y, width, bottom=bottom, label=lab,
                   color=class_color[name], hatch=class_hatch[name],
                   edgecolor=ps.INK, linewidth=0.3)
            bottom += y
        ax.plot(x - width / 2, [m["aras_dates_ha"] for m in monthly],
                **ps.measurement_kw(ps.INK, "D", size=3.2),
                label="ARAS pooled (distinct dates)")
        ax.bar(x + width / 2, [m["rlmt_dates_labelled"] for m in monthly], width,
               color=ps.tint(ps.GOOD), edgecolor=ps.INK, linewidth=0.3,
               label="RLMT grism, labelled")
        ax.bar(x + width / 2, [m["rlmt_dates_dispersed"] for m in monthly], width,
               color=ps.GOOD, edgecolor=ps.INK, linewidth=0.3,
               label="RLMT grism, S2c dispersed")
        # Shade the RLMT 2025 window so the eye finds the comparison.
        lo_m, hi_m = month_key(facts["window_lo"]), month_key(facts["window_hi"])
        ax.axvspan(months.index(lo_m) - 0.5, months.index(hi_m) + 0.5,
                   color=ps.WISP, alpha=0.35, zorder=0, lw=0)
        step = max(1, len(months) // 16)
        ax.set_xticks(x[::step])
        ax.set_xticklabels(months[::step], rotation=45, ha="right")
        ax.set_xlim(-0.8, len(months) - 0.2)
        ax.set_ylabel("dates with an H$\\alpha$ spectrum / month")
        ax.set_xlabel("UT month")
        ax.set_title("T CrB H$\\alpha$ spectroscopy: ARAS database vs RLMT grism "
                     "(shaded: RLMT 2025 window)")
        # Legend below the axes: inside, it would sit on the 2024 bars.
        ax.legend(ncol=3, fontsize="small", loc="upper center",
                  bbox_to_anchor=(0.5, -0.30), frameon=False)
        fig.tight_layout()
        fig.savefig(HERE / "fig_novelty_monthly.png", dpi=ps.WEB_DPI)
        fig.savefig(HERE / "fig_novelty_monthly.pdf")
        plt.close(fig)


def fig_timeline(aras: list[dict], rlmt: list[dict], facts: dict) -> None:
    """One row per series, one mark per night, across the RLMT window.

    The figure a referee would draw to test "densest homogeneous record":
    RLMT's two grisms at the top, every ARAS single-observer series beneath,
    ordered by number of dates, marker shape AND colour by resolution class.
    """
    lo, hi = facts["window_lo"], facts["window_hi"]
    pad = timedelta(days=12)
    inview = lambda d: lo - pad <= d <= hi + pad                 # noqa: E731
    class_style = {"low": (ps.SECOND, "o"), "mid": (ps.WARN, "s"),
                   "high": (ps.BAD, "^")}
    rows = []                                   # (label, [(date, colour, marker, filled)])
    for filt in RLMT_GRISM_2025:
        pts = []
        for r in rlmt:
            if r["filter"] == filt and inview(r["date"]):
                pts.append((r["date"], ps.GOOD, "D", r["dispersed"]))
        rows.append((f"RLMT {filt}", pts))
    groups = defaultdict(list)
    for r in aras:
        if r["ha_loose"] and inview(r["date"]):
            groups[(r["observer"], r["site"])].append(r)
    ranked = sorted(groups.items(),
                    key=lambda kv: -len({r["date"] for r in kv[1]
                                         if lo <= r["date"] <= hi}))
    for (obs, site), rs in ranked:
        pts = [(r["date"], *class_style[r["r_class"]], True) for r in rs]
        n = len({r["date"] for r in rs if lo <= r["date"] <= hi})
        if n == 0:
            continue                # observed only in the padding, not the window
        rows.append((f"ARAS {obs}/{site}  ({n})", pts))
    pooled = [(r["date"], *class_style[r["r_class"]], True)
              for r in aras if r["ha_loose"] and inview(r["date"])]
    rows.insert(2, (f"ARAS pooled  ({facts['aras_dates_ha']})", pooled))
    rows[0] = (f"RLMT hrg  ({len({p[0] for p in rows[0][1] if p[3] and lo <= p[0] <= hi})})",
               rows[0][1])
    rows[1] = (f"RLMT lrg  ({len({p[0] for p in rows[1][1] if p[3] and lo <= p[0] <= hi})})",
               rows[1][1])

    with ps.context("web"):
        fig, ax = plt.subplots(figsize=(ps.COL_DOUBLE * 1.35,
                                        0.24 * len(rows) + 1.3))
        for i, (label, pts) in enumerate(rows):
            y = len(rows) - 1 - i
            ax.axhline(y, color=ps.GRID, lw=0.6, zorder=0)
            for d, colour, marker, filled in pts:
                kw = (ps.measurement_kw(colour, marker, size=4.2) if filled
                      else ps.floor_kw(colour, marker, size=4.2))
                ax.plot([d], [y], **kw)
        ax.axvspan(lo, hi, color=ps.WISP, alpha=0.3, zorder=0, lw=0)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([r[0] for r in rows][::-1], fontsize="small")
        ax.set_ylim(-0.7, len(rows) - 0.3)
        ax.set_xlim(lo - pad, hi + pad)
        ax.xaxis.set_major_locator(mdates.DayLocator(bymonthday=(1, 15)))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
        ax.set_xlabel("UT date (2025)")
        ax.set_title("Who took an H$\\alpha$ spectrum of T CrB on which night\n"
                     "(in brackets: dates inside the shaded RLMT window)")
        handles = [ps.measurement_handle("RLMT, S2c dispersed", ps.GOOD, "D"),
                   plt.Line2D([], [], label="RLMT, labelled only",
                              **ps.floor_kw(ps.GOOD, "D", size=4.0))]
        for name, lo_r, hi_r in R_CLASSES:
            c, mk = class_style[name]
            lab = (f"ARAS {name}-R (< {hi_r:.0f})" if np.isfinite(hi_r)
                   else f"ARAS {name}-R (≥ {lo_r:.0f})")
            handles.append(ps.measurement_handle(lab, c, mk))
        ax.legend(handles=handles, ncol=3, fontsize="small", frameon=False,
                  loc="upper center", bbox_to_anchor=(0.45, -0.07))
        fig.tight_layout()
        fig.savefig(HERE / "fig_novelty_timeline.png", dpi=ps.WEB_DPI)
        fig.savefig(HERE / "fig_novelty_timeline.pdf")
        plt.close(fig)


# ===========================================================================
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refresh", action="store_true",
                    help="re-download the ARAS listing before counting")
    args = ap.parse_args(argv)
    if args.refresh:
        refresh_aras()

    aras, meta = load_aras()
    rlmt = load_rlmt()
    monthly = monthly_table(aras, rlmt)
    window_rows, facts = window_table(aras, rlmt)
    v = verdict(monthly, facts)

    write_csv(HERE / "novelty_monthly.csv", monthly)
    write_csv(HERE / "novelty_window.csv", window_rows)
    write_csv(HERE / "aras_series_window.csv", facts["series"])
    write_csv(HERE / "aras_observers_window.csv", facts["observers"])
    by_date = defaultdict(lambda: defaultdict(int))
    for r in rlmt:
        by_date[(r["date"], r["filter"])]["frames"] += 1
        by_date[(r["date"], r["filter"])]["dispersed"] += int(r["dispersed"])
        by_date[(r["date"], r["filter"])]["header_off"] += int(r["header_off"])
    write_csv(HERE / "rlmt_grism_dates.csv",
              [{"date_ut": d, "filter": f, **dict(c)}
               for (d, f), c in sorted(by_date.items())])
    write_markdown(monthly, window_rows, facts, v, meta, rlmt)
    fig_monthly(monthly, facts)
    fig_timeline(aras, rlmt, facts)

    # A terse console summary; the markdown is the record.
    print(f"ARAS: {meta['n_rows']} spectra ({meta['first_date']} → {meta['last_date']})")
    print(f"window {facts['window_lo']} → {facts['window_hi']} ({facts['window_days']} d)")
    for r in window_rows:
        print(f"  {r['series']:<52s} n={r['n_spectra_or_frames']:>4} "
              f"dates={r['n_dates']:>3} medgap={fmt(r['median_gap_d']):>4} "
              f"maxgap={fmt(r['max_gap_d']):>4} R~{fmt(r['R_median'])}")
    print("verdict:", json.dumps(v, default=str, indent=1))
    (HERE / "novelty_verdict.json").write_text(
        json.dumps({"verdict": v,
                    "facts": {k: val for k, val in facts.items()
                              if k not in ("series", "observers")}},
                   default=str, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
