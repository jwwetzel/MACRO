"""Pure arithmetic for the four remaining CV Phase-2 photometry tasks.

WHY THIS MODULE EXISTS
----------------------
Phase 2 produced 3.17 million light-curve rows and a catalogue tie.  Four
questions were left open, and each one is a place where a plausible-looking
light curve can be quietly wrong:

1.  **Cloud.**  The acquisition software writes a per-image photometric zero
    point (``ZMAG``) into every frame header, and the strategy's original
    cloud cut leaned on it.  For the Sloan-era polar frames that keyword is
    absent or zero on *every single frame* — VV Pup 0 of 1,353, EU UMa 0 of
    208.  So the primary cloud channel for the polars cannot be ZMAG.  It
    has to be the comparison ensemble's own summed flux, which is measured
    on every frame by construction.  :func:`ensemble_flux_ratio`,
    :func:`running_median` and :func:`veto_from_ratio` are that channel.
2.  **Second-order extinction.**  The Honeycutt ensemble absorbs whatever
    is common to a whole frame into ``ZP_j``, which removes first-order
    extinction exactly.  What it cannot remove is the part that depends on a
    star's COLOUR — ``k'' * colour * airmass`` — because that term is
    different for every star in the same frame.  :func:`two_way_center` and
    :func:`fit_kpp` measure it, with an uncertainty, so that "we ignored it"
    becomes a number in the error budget instead of an omission.
3.  **Cross-era transformation.**  ST LMi was observed in G/R/I in 2024 and
    in g/r/i in 2025-26, and the two seasons do not overlap in time.  The
    paper therefore runs two within-era analyses and never stitches them.
    That discipline is only worth anything if it is CHECKED, and the
    transformation between the two natural systems is only publishable if it
    was measured on comparison stars rather than assumed.
    :func:`fit_transform` does the measurement; :func:`ols_line` is the
    two-parameter regression underneath it.
4.  **Faint-phase limits.**  A polar in a low state drops below detection.
    Dropping those epochs is not neutral: it censors exactly the faint half
    of the distribution and biases every duty-cycle statistic upward.
    :func:`similarity_from_pairs` puts the target's position on a frame that
    never detected it, :func:`forced_aperture` measures there,
    :func:`limit_flux` turns the noise into an upper limit, and
    :func:`km_survival` re-derives the state statistics with the limits
    included instead of thrown away.

Everything here is deterministic, dependency-light (numpy only) and pure —
no database, no files, no network.  ``pipeline/scripts/run_cv_phase2.py``
does all of the I/O and calls in here for every number it stores.  Every
function is unit-tested in ``pipeline/tests/test_phase2.py``.

A NOTE ON THE NUMBERS THIS MODULE IS ALLOWED TO ASSUME
------------------------------------------------------
The CV characterization measured what this data set can actually do, and it
is less than the strategy hoped: per-point precision 9-77 mmag by series,
chi2 inflation 0.92-3.02.  Constants below are set against the MEASURED
numbers.  In particular the cloud threshold is not a round number someone
liked; it is calibrated in :func:`roc_table` against the frames that DO
carry an independent ZMAG, and the constant here is only the default the
calibration is expected to confirm or move.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

# ===========================================================================
# 0.  Tunable constants — the single source of truth.  The report reads
#     these rather than repeating them, so changing one changes the page.
# ===========================================================================

#: Half-width, in FRAMES, of the running-median window the cloud statistic
#: is measured against.  The window has to be long enough that a passing
#: cloud cannot become "the local normal" and short enough to track the
#: real, slow drift of airmass and seeing through a night.  The measured
#: per-filter cadence on the densest CV nights is 190-500 s, so +/-10
#: frames is roughly +/-1 hour of one filter's sampling — several times
#: longer than a cumulus transit and several times shorter than the
#: airmass swing of a night.
CLOUD_WINDOW_HALF = 10

#: Minimum number of frames a night must contribute before the running
#: median means anything.  Below this the "local normal" is one or two
#: frames, and a single bad frame would set the baseline it is judged
#: against.
CLOUD_MIN_FRAMES = 7

#: A comparison star joins the cloud ensemble only if it was measured on at
#: least this fraction of the night's frames.  Stars that come and go are
#: exactly the stars that a cloud removes, so admitting them would let the
#: ensemble's MEMBERSHIP absorb the signal the ensemble is meant to show.
CLOUD_CORE_MIN_FRAC = 0.8

#: Minimum core-ensemble size.  Below this the summed flux is dominated by
#: individual stars' own noise rather than by the sky's transparency.
CLOUD_MIN_CORE = 5

#: DEFAULT veto threshold on the flux ratio: a frame is vetoed when the
#: core ensemble delivers less than this fraction of its running-median
#: flux.  0.90 is 0.114 mag of extinction.  This is a starting point; the
#: production threshold is the one :func:`choose_threshold` returns from the
#: ZMAG-calibrated ROC, and the two are compared explicitly on the page.
CLOUD_THRESHOLD_DEFAULT = 0.90

#: The false-veto rate the calibration is allowed to spend.  The veto exists
#: to protect the light curve from cloud, and a veto that throws away 5% of
#: clear frames costs more signal than the cloud it removes — especially on
#: a 14-18 point per cycle cadence, where every frame is a phase bin.
CLOUD_MAX_FALSE_VETO = 0.01

#: Independent-evidence labels, in MAGNITUDES of ZMAG departure from the
#: night's running median.  A frame is called CLEAR when its independent
#: zero point sits within 0.02 mag of the local normal and ATTENUATED when
#: it has lost more than 0.15 mag.  The gap between them is deliberate: the
#: frames in between are neither, and forcing them into one class or the
#: other would manufacture the very agreement the calibration is testing.
ZMAG_CLEAR_MAG = 0.02
ZMAG_ATTEN_MAG = 0.15

#: Airmass window a measurement must fall inside to enter the
#: colour-extinction fit.  NOT a taste judgement: 62 matched CV frames carry
#: header AIRMASS values between 5 and 6,877, and VV Pup — the southernmost
#: target, at dec -19 — cannot exceed airmass 2.1 from this site at any hour
#: of any night.  Those cards are arithmetic wreckage, and because the fitted
#: design column is proportional to (X - Xref), a single X = 6,877 point
#: carries 40,000 times the leverage of a real one and sets the coefficient
#: by itself.  The first run of this stage returned k'' = -4e-5 +/- 1e-5 for
#: era 76 for exactly that reason, which is how the defect was found.
AIRMASS_MIN, AIRMASS_MAX = 1.0, 3.0

#: Residual clip for the colour-extinction fit, in robust sigmas.  Same
#: number the ensemble solver and the catalogue tie use, so that "an
#: outlier" means one thing across the whole project.
KPP_CLIP_SIGMA = 4.0

#: Significance bar for the second-order colour-extinction coefficient.
#: Three sigma, not two: this is one coefficient measured out of a family of
#: (era, filter) combinations, and a 2-sigma bar over ~10 fits produces a
#: "significant" result by arithmetic alone.
KPP_SIGNIFICANCE_T = 3.0

#: Confidence multiplier for an upper limit.  3 sigma of the measured
#: aperture noise; one-sided Gaussian, 99.87%.
LIMIT_SIGMA = 3.0

#: Forced photometry only reports a limit when the position it measured at
#: is trustworthy.  The frame-to-reference transform must have been fitted
#: from at least this many matched stars...
FORCED_MIN_STARS = 6

#: ...and must close on them to better than this RMS, in pixels.  The
#: aperture radius is 4.0 arcsec (9-18 px depending on era), so a transform
#: good to 1.5 px puts the target well inside its own aperture.
FORCED_MAX_RMS_PX = 1.5

#: A forced measurement whose signal-to-noise reaches this is a DETECTION
#: that source detection missed, not an upper limit, and is reported as
#: such.  Same number as the limit multiplier on purpose: "we can see it"
#: and "we can bound it" must be the same threshold, or the two categories
#: would overlap.
FORCED_DETECT_SNR = 3.0

#: A block may publish upper limits only if its forced POSITION has been
#: validated on frames where the target was actually detected: at least this
#: many such frames, closing to at most :data:`CLOSURE_MAX_MEDIAN_PX`.
#:
#: This gate exists because the first production run produced 66 "forced
#: detections" of EU UMa in the merged 2026 Fast block (era 78) — the block
#: the characterization already flags as carrying five comparison stars,
#: zero check stars and no error validation at all.  Its frame-to-reference
#: transforms closed to 645-1,650 pixels on 87 of 153 attempts, and where
#: they closed tightly the measured signal-to-noise alternated between 3 and
#: 55 frame to frame — an aperture landing on a bright neighbour on half the
#: frames.  That block has never detected its target, so there is no frame
#: on which the forced position can be checked, and a limit measured at an
#: unverifiable position is not a limit, it is a number.  The gate refuses
#: it by rule rather than leaving a reader to notice.
CLOSURE_MIN_FRAMES = 10
CLOSURE_MAX_MEDIAN_PX = 1.0

#: Instrumental-magnitude zero offset, copied from
#: ``macro_phot.photometry.INST_MAG_OFFSET`` so that a forced magnitude
#: lands on the same scale as a detected one.  Duplicated as a named
#: constant rather than imported because this module is meant to stay
#: importable with numpy alone.
INST_MAG_OFFSET = 25.0


# ===========================================================================
# 1.  Small robust statistics used by more than one task
# ===========================================================================
def median_abs_deviation(x: Sequence[float]) -> float:
    """Median absolute deviation, scaled to a Gaussian sigma.

    Returned as ``nan`` for fewer than two finite values: a "spread" of one
    number is not zero, it is unknown, and returning 0.0 would make a
    single-frame night look like the steadiest night of the campaign.
    """
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size < 2:
        return float("nan")
    med = float(np.median(a))
    return float(1.4826 * np.median(np.abs(a - med)))


def running_median(y: Sequence[float], half_width: int) -> np.ndarray:
    """Centred running median of ``y`` over ``+/- half_width`` samples.

    The window SHRINKS at the ends rather than padding or reflecting.  A
    reflected edge invents data, and the invented data would sit exactly
    where the first and last frames of a night live — the frames most likely
    to be taken through horizon murk, i.e. the ones the cloud veto most
    needs to judge honestly.

    NaNs are ignored inside each window; a window with no finite value
    yields NaN, which downstream code reads as "no local normal here" and
    declines to veto rather than vetoing blind.
    """
    a = np.asarray(y, dtype=float)
    n = a.size
    out = np.full(n, np.nan)
    if n == 0:
        return out
    h = max(0, int(half_width))
    for i in range(n):
        lo, hi = max(0, i - h), min(n, i + h + 1)
        w = a[lo:hi]
        w = w[np.isfinite(w)]
        if w.size:
            out[i] = float(np.median(w))
    return out


def mann_whitney_u(a: Sequence[float], b: Sequence[float]
                   ) -> tuple[float, float, float]:
    """Two-sided Mann-Whitney U test by the normal approximation.

    Returns ``(U, z, p)``.  Written out here rather than imported because
    this is the test the whole "does the veto sculpt the light curve?"
    argument rests on, and a reader auditing that argument should be able to
    read the test itself in the same repository as the claim.

    Ties are handled by mid-ranks and the tie correction to the variance,
    which matters here: magnitudes rounded through a float column produce
    few exact ties, but a series where most frames share a saturation flag
    produces many.

    A sample smaller than two on either side returns ``(nan, nan, nan)``:
    with one point there is no distribution to compare, and returning p=1
    would read as "tested and cleared".
    """
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    n1, n2 = x.size, y.size
    if n1 < 2 or n2 < 2:
        return float("nan"), float("nan"), float("nan")
    both = np.concatenate([x, y])
    order = np.argsort(both, kind="mergesort")
    ranks = np.empty(both.size, dtype=float)
    sorted_vals = both[order]
    i = 0
    tie_term = 0.0
    while i < sorted_vals.size:
        j = i
        while j + 1 < sorted_vals.size and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        mid = 0.5 * (i + j) + 1.0          # 1-based mid-rank
        ranks[order[i:j + 1]] = mid
        t = j - i + 1
        tie_term += t ** 3 - t
        i = j + 1
    r1 = float(ranks[:n1].sum())
    u1 = r1 - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    mu = n1 * n2 / 2.0
    var = n1 * n2 * (n + 1) / 12.0 - n1 * n2 * tie_term / (12.0 * n * (n - 1))
    if var <= 0:
        return float(u1), float("nan"), float("nan")
    z = (u1 - mu) / math.sqrt(var)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return float(u1), float(z), float(p)


def two_proportion_z(k1: int, n1: int, k2: int, n2: int
                     ) -> tuple[float, float, float]:
    """Two-sided z test on two proportions.  Returns ``(diff, z, p)``.

    Used to ask "is the veto rate in the target's faint quartile different
    from the veto rate in its bright quartile?" — the sharpest possible
    version of the sculpting question, because it compares the two ends of
    the light curve directly instead of averaging over the middle.
    """
    if n1 <= 0 or n2 <= 0:
        return float("nan"), float("nan"), float("nan")
    p1, p2 = k1 / n1, k2 / n2
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1.0 - p) * (1.0 / n1 + 1.0 / n2))
    if se <= 0:
        return p1 - p2, float("nan"), float("nan")
    z = (p1 - p2) / se
    return float(p1 - p2), float(z), float(math.erfc(abs(z) / math.sqrt(2.0)))


# ===========================================================================
# 2.  TASK 1 — the ensemble-flux-ratio cloud veto
# ===========================================================================
def core_ensemble(star_ids: Sequence[int], frame_ids: Sequence[int],
                  min_frac: float = CLOUD_CORE_MIN_FRAC) -> list[int]:
    """Which stars form the cloud ensemble on this night.

    A star qualifies when it was measured on at least ``min_frac`` of the
    night's frames.  This is the single most important rule in the whole
    veto, and it is worth saying why in full:

    the statistic below is a SUM of comparison flux.  If the membership of
    the sum is allowed to change from frame to frame, then a cloud that
    hides the faintest three stars reduces the sum twice over — once
    because the surviving stars are dimmer, and once because three
    contributors vanished.  The second effect is not transparency, it is
    bookkeeping, and it is unbounded.  Fixing the membership to stars that
    are essentially always there makes the sum a measurement of the sky
    rather than a measurement of the detector's detection threshold.
    """
    sid = np.asarray(star_ids)
    fid = np.asarray(frame_ids)
    n_frames = np.unique(fid).size
    if n_frames == 0:
        return []
    out = []
    for s in np.unique(sid):
        seen = np.unique(fid[sid == s]).size
        if seen >= min_frac * n_frames:
            out.append(int(s))
    return sorted(out)


def ensemble_flux_ratio(star_ids: Sequence[int], frame_index: Sequence[int],
                        flux_rate: Sequence[float], n_frames: int,
                        core: Sequence[int]) -> tuple[np.ndarray, np.ndarray]:
    """The ensemble's summed comparison flux per frame, on a common scale.

    Parameters
    ----------
    star_ids, frame_index, flux_rate
        One entry per (star, frame) measurement.  ``frame_index`` is a
        0-based index into the night's frame ordering; ``flux_rate`` is flux
        per second (flux divided by exposure time), so that a night mixing
        60 s and 240 s exposures — which every EU UMa night does — is not
        read as a night mixing clear and cloudy sky.
    n_frames
        Length of the night.
    core
        The star ids from :func:`core_ensemble`.

    Returns
    -------
    (ratio, n_used)
        ``ratio[j]`` is ``sum_s f_sj / sum_s fbar_s`` over the core stars
        PRESENT on frame j, where ``fbar_s`` is star s's median flux rate
        over the night.  Dividing by the reference sum of the same subset is
        what makes a frame that lost one core star comparable with a frame
        that kept them all: both are then "what fraction of its own normal
        did this ensemble deliver".  ``n_used[j]`` is that subset's size, so
        a reader can see when a ratio rests on three stars instead of forty.

    Frames with no core star present get ``nan`` and are never vetoed: an
    unmeasurable frame is not a cloudy frame, and conflating the two would
    veto every frame the extraction happened to fail on.
    """
    sid = np.asarray(star_ids)
    fj = np.asarray(frame_index, dtype=int)
    fr = np.asarray(flux_rate, dtype=float)
    core_set = {int(c) for c in core}
    keep = np.array([int(s) in core_set for s in sid], dtype=bool)
    keep &= np.isfinite(fr) & (fr > 0)
    sid, fj, fr = sid[keep], fj[keep], fr[keep]
    # Each core star's own reference level: the median of its flux rate over
    # the night.  The median, not the mean, because a night with two cloudy
    # frames must not have its baseline dragged down by them.
    fbar = {}
    for s in np.unique(sid):
        v = fr[sid == s]
        if v.size:
            fbar[int(s)] = float(np.median(v))
    num = np.zeros(n_frames)
    den = np.zeros(n_frames)
    cnt = np.zeros(n_frames, dtype=int)
    for s, j, f in zip(sid, fj, fr):
        b = fbar.get(int(s))
        if b is None or not math.isfinite(b) or b <= 0:
            continue
        num[j] += f
        den[j] += b
        cnt[j] += 1
    ratio = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return ratio, cnt


def veto_from_ratio(ratio: Sequence[float], threshold: float,
                    half_width: int = CLOUD_WINDOW_HALF
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Veto mask and the local normal each frame was judged against.

    Returns ``(vetoed, local)``.  ``local`` is the running median of the
    ratio; ``vetoed[j]`` is True when ``ratio[j] / local[j] < threshold``.

    Judging against a RUNNING median rather than the night's median is what
    makes this a cloud detector rather than an airmass detector.  Over a
    six-hour run the ensemble's summed flux falls smoothly by tens of
    percent as the field sets; that is extinction, it is already absorbed
    exactly by the ensemble's per-frame zero point, and vetoing it would
    delete the end of every night.  A cloud, by contrast, is a departure
    from the local normal on a timescale of minutes.

    A frame whose local normal is NaN — the ends of a night too short to
    supply a window — is NOT vetoed.  See :func:`running_median`.
    """
    r = np.asarray(ratio, dtype=float)
    local = running_median(r, half_width)
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.where(np.isfinite(local) & (local > 0), r / local, np.nan)
    vetoed = np.isfinite(rel) & (rel < float(threshold))
    return vetoed, local


def relative_ratio(ratio: Sequence[float],
                   half_width: int = CLOUD_WINDOW_HALF) -> np.ndarray:
    """``ratio / running_median(ratio)`` — the veto statistic itself."""
    r = np.asarray(ratio, dtype=float)
    local = running_median(r, half_width)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(np.isfinite(local) & (local > 0), r / local, np.nan)


def zmag_transmission(zmag: Sequence[float],
                      half_width: int = CLOUD_WINDOW_HALF) -> np.ndarray:
    """Independent transmission from the header zero point, as a fraction.

    ``ZMAG`` is the magnitude of a 1 ADU/s source as the acquisition
    software's own plate solve measured it: a frame that lost half its light
    to cloud has a ZMAG 0.75 mag SMALLER.  Referencing it to its own running
    median, exactly as the ensemble ratio is referenced, gives a number on
    the same footing that owes nothing to our photometry.

    This is the only genuinely independent cloud channel this archive has,
    and it exists for 2,506 CV frames and for none of the Sloan-era polar
    frames — which is precisely why the ensemble ratio has to be primary,
    and precisely why these 2,506 frames are worth calibrating on.
    """
    z = np.asarray(zmag, dtype=float)
    z = np.where(np.isfinite(z) & (z != 0.0), z, np.nan)
    local = running_median(z, half_width)
    with np.errstate(invalid="ignore"):
        return np.where(np.isfinite(local), 10.0 ** (0.4 * (z - local)),
                        np.nan)


def label_by_zmag(transmission: Sequence[float],
                  clear_mag: float = ZMAG_CLEAR_MAG,
                  atten_mag: float = ZMAG_ATTEN_MAG) -> np.ndarray:
    """Turn independent transmission into ``'clear'`` / ``'attenuated'`` /
    ``''`` labels.  The empty label is the intentional third class — see
    :data:`ZMAG_CLEAR_MAG`.
    """
    t = np.asarray(transmission, dtype=float)
    out = np.full(t.size, "", dtype=object)
    with np.errstate(invalid="ignore", divide="ignore"):
        dmag = np.where(t > 0, -2.5 * np.log10(t), np.nan)
    out[np.isfinite(dmag) & (np.abs(dmag) <= clear_mag)] = "clear"
    out[np.isfinite(dmag) & (dmag >= atten_mag)] = "attenuated"
    return out


def roc_table(rel_ratio: Sequence[float], labels: Sequence[str],
              thresholds: Sequence[float]) -> list[dict]:
    """False-veto and recall rate at each candidate threshold.

    One row per threshold: how many independently-CLEAR frames it would
    throw away (the cost) and how many independently-ATTENUATED frames it
    would catch (the benefit).  This is the whole calibration: the threshold
    is not chosen because 0.9 is a nice number, it is chosen because at that
    value the cost is bounded and the benefit is measured.
    """
    r = np.asarray(rel_ratio, dtype=float)
    lab = np.asarray(labels, dtype=object)
    clear = (lab == "clear") & np.isfinite(r)
    atten = (lab == "attenuated") & np.isfinite(r)
    n_clear, n_atten = int(clear.sum()), int(atten.sum())
    rows = []
    for t in thresholds:
        vc = int((r[clear] < t).sum())
        va = int((r[atten] < t).sum())
        rows.append({
            "threshold": float(t),
            "n_clear": n_clear, "n_clear_vetoed": vc,
            "false_veto_rate": (vc / n_clear) if n_clear else float("nan"),
            "n_attenuated": n_atten, "n_attenuated_vetoed": va,
            "recall": (va / n_atten) if n_atten else float("nan"),
        })
    return rows


def choose_threshold(rows: Sequence[dict],
                     max_false_veto: float = CLOUD_MAX_FALSE_VETO
                     ) -> tuple[Optional[float], str]:
    """The highest (most aggressive) threshold whose false-veto rate is
    still within budget, with the sentence that defends it.

    "Highest" is the right direction: a higher threshold vetoes more, so
    among the thresholds that stay inside the false-veto budget we want the
    one that catches the most cloud.  Returns ``(None, reason)`` when the
    calibration set cannot support any choice, which is an answer and not a
    failure — it means the default constant must be quoted as a default.
    """
    ok = [r for r in rows
          if math.isfinite(r["false_veto_rate"])
          and r["false_veto_rate"] <= max_false_veto]
    if not ok:
        return None, ("no threshold keeps the false-veto rate within "
                      f"{max_false_veto:.1%} on the independently-labelled "
                      "frames")
    best = max(ok, key=lambda r: r["threshold"])
    return best["threshold"], (
        f"highest threshold whose false-veto rate stays within "
        f"{max_false_veto:.1%}: {best['n_clear_vetoed']}/{best['n_clear']} "
        f"independently-clear frames vetoed "
        f"({best['false_veto_rate']:.2%}), catching "
        f"{best['n_attenuated_vetoed']}/{best['n_attenuated']} "
        f"independently-attenuated frames ({best['recall']:.1%})")


def sculpting_test(target_mag: Sequence[float], vetoed: Sequence[bool]
                   ) -> dict:
    """Does the veto preferentially remove the target's FAINT phases?

    This is the test that decides whether the veto is a cleaning step or a
    light-curve editor.  A veto that fires more often when the target is
    faint would carve the bottom out of every eclipse and every low state,
    and would do it invisibly, because the survivors would still look like a
    clean light curve.

    Two independent readings of the same question:

    * a Mann-Whitney U on the target magnitudes of vetoed vs kept frames —
      sensitive to any systematic shift of the whole distribution;
    * a two-proportion z on the veto rate in the target's FAINTEST quartile
      against its BRIGHTEST quartile — sensitive to an effect concentrated
      at the extremes, which is where an eclipse lives and where a
      whole-distribution test is weakest.

    The construction argues the answer should be null: the ratio is measured
    on comparison stars only, and the target's magnitude is differential, so
    cloud dims both and cancels.  But "should" is not "does", and the two
    can couple through the detection threshold — a cloudy frame loses the
    faint end of the ensemble AND the faint end of the target.  Hence the
    measurement.

    THE DIRECTION IS PART OF THE ANSWER.  Both tests are two-sided, because
    a one-sided test would be a decision to only notice the result we fear.
    But a significant asymmetry in which the BRIGHT quartile is vetoed more
    often is not the failure this test is looking for — it cannot carve a
    low state out of a light curve — so the verdict names the direction
    instead of collapsing both into one alarm.  Calling a bright-side excess
    "sculpting" would train the reader to ignore the word.
    """
    m = np.asarray(target_mag, dtype=float)
    v = np.asarray(vetoed, dtype=bool)
    good = np.isfinite(m)
    m, v = m[good], v[good]
    out = {"n": int(m.size), "n_vetoed": int(v.sum()),
           "median_mag_vetoed": float("nan"), "median_mag_kept": float("nan"),
           "u": float("nan"), "z": float("nan"), "p_mannwhitney": float("nan"),
           "faint_veto_rate": float("nan"), "bright_veto_rate": float("nan"),
           "prop_diff": float("nan"), "prop_z": float("nan"),
           "p_proportion": float("nan"), "verdict": "NOT TESTABLE"}
    if m.size < 8 or v.sum() == 0 or (~v).sum() == 0:
        out["verdict"] = ("NOT TESTABLE — fewer than 8 measured epochs, or "
                          "nothing vetoed, or everything vetoed")
        return out
    out["median_mag_vetoed"] = float(np.median(m[v]))
    out["median_mag_kept"] = float(np.median(m[~v]))
    u, z, p = mann_whitney_u(m[v], m[~v])
    out["u"], out["z"], out["p_mannwhitney"] = u, z, p
    # Quartiles in MAGNITUDE: larger magnitude = fainter.
    q1, q3 = float(np.percentile(m, 25)), float(np.percentile(m, 75))
    bright = m <= q1
    faint = m >= q3
    d, pz, pp = two_proportion_z(int(v[faint].sum()), int(faint.sum()),
                                 int(v[bright].sum()), int(bright.sum()))
    out["faint_veto_rate"] = (float(v[faint].sum()) / faint.sum()
                              if faint.sum() else float("nan"))
    out["bright_veto_rate"] = (float(v[bright].sum()) / bright.sum()
                               if bright.sum() else float("nan"))
    out["prop_diff"], out["prop_z"], out["p_proportion"] = d, pz, pp
    ps = [x for x in (p, pp) if math.isfinite(x)]
    if not ps:
        out["verdict"] = "NOT TESTABLE — no finite test statistic"
        return out
    if min(ps) >= 0.05 / max(1, len(ps)):
        out["verdict"] = "NO SCULPTING DETECTED"
        return out
    # Significant.  Which way?  "Faint-side" means the vetoed frames are the
    # FAINTER ones: the faint quartile is vetoed more often, or the vetoed
    # frames' median magnitude is larger (fainter) than the kept frames'.
    faint_side = 0
    if math.isfinite(d):
        faint_side += 1 if d > 0 else -1
    if (math.isfinite(out["median_mag_vetoed"])
            and math.isfinite(out["median_mag_kept"])):
        faint_side += (1 if out["median_mag_vetoed"]
                       > out["median_mag_kept"] else -1)
    out["verdict"] = ("FAINT-PHASE VETO EXCESS" if faint_side > 0
                      else "BRIGHT-PHASE VETO EXCESS")
    return out


# ===========================================================================
# 3.  TASK 2 — second-order colour-extinction terms
# ===========================================================================
def two_way_center(values: Sequence[float], row: Sequence[int],
                   col: Sequence[int], n_iter: int = 12) -> np.ndarray:
    """Remove the row means and the column means from a sparse table.

    WHY THIS EXISTS, AND WHY THE FIT IS WRONG WITHOUT IT
    ----------------------------------------------------
    The quantity being fitted is a residual that the Honeycutt solver has
    already reduced twice.  The solver's model is
    ``m_sj = M_s + ZP_j``: it has removed a free constant per STAR and a
    free constant per FRAME.  A design column that still contains star means
    and frame means is therefore partly in the space the solver already
    projected out, and the coefficient that comes back is diluted by an
    unknown amount.

    So the design column is projected the same way the data were.  With a
    complete table one pass of "subtract row means, subtract column means"
    would do it; the tables here are sparse (a star is missing from the
    frames it was not detected on) and the two projections do not commute,
    so the subtraction is iterated to convergence.  A dozen sweeps is far
    more than these tables need — they converge in three or four — and
    costs nothing.
    """
    v = np.asarray(values, dtype=float).copy()
    r = np.asarray(row, dtype=int)
    c = np.asarray(col, dtype=int)
    if v.size == 0:
        return v
    n_r, n_c = int(r.max()) + 1, int(c.max()) + 1
    for _ in range(int(n_iter)):
        for idx, n in ((r, n_r), (c, n_c)):
            tot = np.bincount(idx, weights=v, minlength=n)
            cnt = np.bincount(idx, minlength=n).astype(float)
            mean = np.where(cnt > 0, tot / np.where(cnt > 0, cnt, 1.0), 0.0)
            v = v - mean[idx]
    return v


@dataclass(frozen=True)
class KppFit:
    """One second-order colour-extinction coefficient and its context."""

    n_points: int
    n_stars: int
    n_frames: int
    n_clipped: int
    kpp: float
    kpp_err: float
    t_stat: float
    p_value: float
    significant: bool
    chi2nu: float
    rms_before: float
    rms_after: float
    colour_ref: float
    airmass_ref: float
    span_term_p95: float


def _empty_kpp(n: int, n_star: int, n_frame: int, c_ref=float("nan"),
               x_ref=float("nan")) -> KppFit:
    """A KppFit that says "not measurable" without pretending otherwise."""
    nan = float("nan")
    return KppFit(n, n_star, n_frame, 0, nan, nan, nan, nan, False, nan,
                  nan, nan, c_ref, x_ref, nan)


def fit_kpp(resid: Sequence[float], sigma: Sequence[float],
            colour: Sequence[float], airmass: Sequence[float],
            star_index: Sequence[int], frame_index: Sequence[int],
            t_bar: float = KPP_SIGNIFICANCE_T,
            clip_sigma: float = KPP_CLIP_SIGMA, passes: int = 2) -> KppFit:
    """Fit ``resid = k'' * (colour - Cref) * (airmass - Xref)``.

    One free parameter, weighted by ``1/sigma^2``, on a design column that
    has been two-way centred (see :func:`two_way_center`), with two robust
    clipping passes on the residuals.

    The uncertainty returned is the formal one INFLATED by ``sqrt(chi2nu)``
    when the fit is over-dispersed, which it always is here: the CV
    characterization measured chi2 inflation of 0.92-3.02 across these
    series, so the formal error understates the real one by up to a factor
    1.7.  A coefficient declared significant on an un-inflated error would
    be an artefact of that known under-estimate.

    ``span_term_p95`` is the 95th percentile of ``|k'' * z|`` over the data
    actually fitted: the size of the effect on a real measurement, which is
    the number the error budget wants, as opposed to the coefficient, which
    is the number a referee wants.
    """
    r = np.asarray(resid, dtype=float)
    s = np.asarray(sigma, dtype=float)
    c = np.asarray(colour, dtype=float)
    x = np.asarray(airmass, dtype=float)
    si = np.asarray(star_index, dtype=int)
    fi = np.asarray(frame_index, dtype=int)
    ok = (np.isfinite(r) & np.isfinite(c) & np.isfinite(x)
          & np.isfinite(s) & (s > 0))
    r, s, c, x, si, fi = r[ok], s[ok], c[ok], x[ok], si[ok], fi[ok]
    n_star, n_frame = int(np.unique(si).size), int(np.unique(fi).size)
    if r.size < 20 or n_star < 3 or n_frame < 3:
        return _empty_kpp(int(r.size), n_star, n_frame)
    # Re-index star/frame to dense 0..n-1 so bincount in two_way_center is
    # not asked to allocate an array the size of the largest frame_id.
    _, si_d = np.unique(si, return_inverse=True)
    _, fi_d = np.unique(fi, return_inverse=True)
    c_ref = float(np.median(c))
    x_ref = float(np.median(x))
    z = two_way_center((c - c_ref) * (x - x_ref), si_d, fi_d)
    w = 1.0 / (s * s)
    keep = np.ones(r.size, dtype=bool)
    k = err = chi2nu = float("nan")
    for _ in range(max(1, int(passes)) + 1):
        denom = float(np.sum(w[keep] * z[keep] * z[keep]))
        if not math.isfinite(denom) or denom <= 0 or int(keep.sum()) < 20:
            return _empty_kpp(int(r.size), n_star, n_frame, c_ref, x_ref)
        k = float(np.sum(w[keep] * z[keep] * r[keep]) / denom)
        dof = max(1, int(keep.sum()) - 1)
        chi2nu = float(np.sum(w[keep] * (r[keep] - k * z[keep]) ** 2) / dof)
        err = math.sqrt(1.0 / denom) * math.sqrt(max(1.0, chi2nu))
        full = r - k * z
        med = float(np.median(full))
        mad = median_abs_deviation(full)
        if not math.isfinite(mad) or mad <= 0:
            break
        new_keep = np.abs(full - med) <= clip_sigma * mad
        if int(new_keep.sum()) < 20 or bool((new_keep == keep).all()):
            break
        keep = new_keep
    t = k / err if err and err > 0 else float("nan")
    # Two-sided Gaussian tail.  A t distribution with >10^4 dof is Gaussian
    # to far more places than this measurement can distinguish.
    p = math.erfc(abs(t) / math.sqrt(2.0)) if math.isfinite(t) else float("nan")
    model = k * z[keep]
    return KppFit(
        n_points=int(keep.sum()), n_stars=n_star, n_frames=n_frame,
        n_clipped=int(r.size - keep.sum()),
        kpp=k, kpp_err=float(err), t_stat=float(t), p_value=float(p),
        significant=bool(math.isfinite(t) and abs(t) >= t_bar),
        chi2nu=chi2nu,
        rms_before=float(np.sqrt(np.mean(r[keep] ** 2))),
        rms_after=float(np.sqrt(np.mean((r[keep] - model) ** 2))),
        colour_ref=c_ref, airmass_ref=x_ref,
        span_term_p95=float(np.percentile(np.abs(model), 95)))


def bootstrap_kpp_error(resid: Sequence[float], sigma: Sequence[float],
                        colour: Sequence[float], airmass: Sequence[float],
                        star_index: Sequence[int],
                        frame_index: Sequence[int],
                        n_boot: int = 24, seed: int = 20260819) -> float:
    """Uncertainty on k'' from resampling STARS, not points.

    WHY THE FORMAL ERROR IS NOT ENOUGH
    ----------------------------------
    The formal weighted-least-squares error treats every (star, frame)
    measurement as an independent draw.  They are not.  A comparison star
    whose catalogue colour is 0.05 mag wrong contributes the SAME wrong
    colour to all 1,600 of its measurements, so its error enters the fit
    1,600 times with the same sign.  With half a million points and 1,800
    stars, the formal error is set by the point count and the real
    uncertainty is set by the star count — a factor of 17 apart in the
    square root alone.

    Resampling whole stars with replacement rebuilds the fit from a
    different draw of the thing that actually varies independently.  The
    spread of the resulting coefficients is the uncertainty a referee should
    be quoted, and it is the one this stage publishes whenever it exceeds
    the formal one.

    Returns NaN when there are too few stars for the resample to mean
    anything, in which case the formal error stands and the report says so.
    """
    si = np.asarray(star_index)
    uniq, inv = np.unique(si, return_inverse=True)
    if uniq.size < 12 or n_boot < 4:
        return float("nan")
    r = np.asarray(resid, dtype=float)
    s = np.asarray(sigma, dtype=float)
    c = np.asarray(colour, dtype=float)
    x = np.asarray(airmass, dtype=float)
    fi = np.asarray(frame_index)
    # Index lists per star, built once; the loop below only concatenates.
    order = np.argsort(inv, kind="mergesort")
    starts = np.searchsorted(inv[order], np.arange(uniq.size))
    ends = np.searchsorted(inv[order], np.arange(uniq.size), side="right")
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(int(n_boot)):
        pick = rng.integers(0, uniq.size, uniq.size)
        idx = np.concatenate([order[starts[p]:ends[p]] for p in pick])
        if idx.size < 20:
            continue
        # Each resampled COPY of a star must be its own star in the
        # centring, or the duplicates would be pooled and the resample
        # would not be a resample.
        rep = np.concatenate([np.full(ends[p] - starts[p], i)
                              for i, p in enumerate(pick)])
        f = fit_kpp(r[idx], s[idx], c[idx], x[idx], rep, fi[idx])
        if math.isfinite(f.kpp):
            vals.append(f.kpp)
    if len(vals) < 4:
        return float("nan")
    return float(np.std(vals, ddof=1))


def bound_mmag(fit: KppFit, colour_span: float, airmass_span: float,
               t_bar: float = KPP_SIGNIFICANCE_T) -> float:
    """The worst-case systematic, in mmag, of NOT applying the term.

    ``t_bar * kpp_err`` is the largest coefficient the data still allow;
    half the colour span times half the airmass span is the largest value of
    ``(C - Cref)(X - Xref)`` a real measurement reaches.  Their product is
    the honest ceiling on what the omission can cost, and it is that ceiling
    — not the fitted point value — that belongs in an error budget when the
    fit is consistent with zero.
    """
    if fit.kpp_err is None or not math.isfinite(fit.kpp_err):
        return float("nan")
    return float(1000.0 * t_bar * fit.kpp_err
                 * abs(colour_span) / 2.0 * abs(airmass_span) / 2.0)


def budget_sentence(fit: KppFit, colour_span: float, airmass_span: float,
                    t_bar: float = KPP_SIGNIFICANCE_T) -> str:
    """What a fitted (or unfitted) coefficient means for the error budget.

    The house rule this implements: an insignificant coefficient is NOT
    forced into the correction.  It is converted into a BOUND and the bound
    is carried in the budget.  Applying a term measured to be consistent
    with zero adds its own noise to every point and buys nothing — and on a
    data set whose per-point precision is 9-77 mmag, adding noise to chase a
    1 mmag effect is a net loss that would never show up in a residual plot.
    """
    if not math.isfinite(fit.kpp):
        return ("not fitted — too few points, stars or frames to identify a "
                "colour-airmass term at all")
    ceiling = bound_mmag(fit, colour_span, airmass_span, t_bar)
    if fit.significant:
        return (f"significant at {abs(fit.t_stat):.1f} sigma; the term moves "
                f"a real measurement by up to "
                f"{1000 * fit.span_term_p95:.1f} mmag (95th percentile of "
                f"|k''z| over the fitted points)")
    return (f"consistent with zero ({abs(fit.t_stat):.1f} sigma); "
            f"|k''| < {t_bar * fit.kpp_err:.4f} mag mag^-1 airmass^-1 at "
            f"{t_bar:.0f} sigma, which over the colour and airmass span "
            f"actually observed can cost at most {ceiling:.1f} mmag — "
            "carried in the error budget as a bound, not applied as a "
            "correction")


# ===========================================================================
# 4.  TASK 3 — cross-era transformation coefficients and the discipline
# ===========================================================================
@dataclass(frozen=True)
class LineFit:
    """``y = a + b * (x - x_ref)`` with uncertainties."""

    n: int
    a: float
    a_err: float
    b: float
    b_err: float
    x_ref: float
    rms: float
    chi2nu: float
    x_min: float
    x_max: float


def ols_line(x: Sequence[float], y: Sequence[float],
             sigma: Optional[Sequence[float]] = None,
             clip_sigma: float = 4.0, passes: int = 2) -> LineFit:
    """Weighted straight-line fit with sigma clipping, centred on ``median(x)``.

    Centring is not cosmetic.  Fitting ``y = a + b*x`` with x running from
    0.2 to 1.4 makes ``a`` an extrapolation to zero colour — a place no star
    in this campaign occupies — and correlates it strongly with ``b``, so
    the quoted ``a_err`` describes a number nobody wants.  Centring on the
    median colour makes ``a`` the offset AT THE COLOUR THE STARS ACTUALLY
    HAVE, decorrelated from the slope, which is the coefficient a data
    release should publish.

    Sigma clipping runs on the residuals with a MAD scale, twice, matching
    the convention the ensemble and the catalogue tie already use, so that
    "an outlier" means the same thing on every page of this project.
    """
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    sa = (np.ones_like(xa) if sigma is None
          else np.asarray(sigma, dtype=float))
    ok = np.isfinite(xa) & np.isfinite(ya) & np.isfinite(sa) & (sa > 0)
    xa, ya, sa = xa[ok], ya[ok], sa[ok]
    nan = float("nan")
    if xa.size < 3:
        return LineFit(int(xa.size), nan, nan, nan, nan, nan, nan, nan,
                       nan, nan)
    x_ref = float(np.median(xa))
    keep = np.ones(xa.size, dtype=bool)
    a = b = ae = be = rms = chi2nu = nan
    for _ in range(max(1, int(passes)) + 1):
        xk, yk, sk = xa[keep] - x_ref, ya[keep], sa[keep]
        if xk.size < 3:
            break
        w = 1.0 / (sk * sk)
        sw = float(w.sum())
        sx = float((w * xk).sum())
        sxx = float((w * xk * xk).sum())
        sy = float((w * yk).sum())
        sxy = float((w * xk * yk).sum())
        det = sw * sxx - sx * sx
        if det <= 0:
            break
        a = (sxx * sy - sx * sxy) / det
        b = (sw * sxy - sx * sy) / det
        resid = yk - (a + b * xk)
        dof = max(1, xk.size - 2)
        chi2nu = float((w * resid ** 2).sum() / dof)
        scale = math.sqrt(max(1.0, chi2nu))
        ae = math.sqrt(sxx / det) * scale
        be = math.sqrt(sw / det) * scale
        rms = float(np.sqrt(np.mean(resid ** 2)))
        # Re-clip on the full sample against the current model.
        full_resid = ya - (a + b * (xa - x_ref))
        med = float(np.median(full_resid))
        mad = median_abs_deviation(full_resid)
        if not math.isfinite(mad) or mad <= 0:
            break
        new_keep = np.abs(full_resid - med) <= clip_sigma * mad
        if new_keep.sum() < 3 or bool((new_keep == keep).all()):
            keep = new_keep if new_keep.sum() >= 3 else keep
            break
        keep = new_keep
    return LineFit(int(keep.sum()), float(a), float(ae), float(b), float(be),
                   x_ref, float(rms), float(chi2nu),
                   float(np.min(xa[keep])) if keep.any() else nan,
                   float(np.max(xa[keep])) if keep.any() else nan)


def fit_transform(mag_from: Sequence[float], mag_to: Sequence[float],
                  colour: Sequence[float],
                  sigma: Optional[Sequence[float]] = None) -> LineFit:
    """Coefficients of ``m_to - m_from = a + b * (colour - colour_ref)``.

    Fitted on COMPARISON STARS, which is the whole point.  Both magnitudes
    are natural-system magnitudes already zero-pointed to the same
    catalogue, so ``a`` is what is left of the zero point after that tie
    (it should be small, and if it is not, the tie is telling us something)
    and ``b`` is the bandpass difference — the difference of the two eras'
    colour terms, measured on stars rather than assumed from a filter
    curve.

    Published as METADATA.  It is never applied to a target magnitude; see
    :func:`verify_no_target_transform`, which checks that it was not.
    """
    mf = np.asarray(mag_from, dtype=float)
    mt = np.asarray(mag_to, dtype=float)
    return ols_line(colour, mt - mf, sigma)


def verify_no_era_mixing(series_eras: Sequence[tuple[str, int, int]]
                         ) -> tuple[int, list[str]]:
    """Check that no series contains a frame from another era.

    ``series_eras`` is ``(series_key, series_era, frame_era)`` triples, one
    per frame.  Returns ``(n_violations, examples)``.

    This is the assertion that carries the paper's cross-era discipline.
    The rule was designed into the series key — ``target|eNN|filter`` — but
    a key is a naming convention, and a naming convention is not a
    guarantee.  Checking the frames themselves turns it into one.
    """
    bad: list[str] = []
    for key, s_era, f_era in series_eras:
        if int(s_era) != int(f_era):
            bad.append(f"{key}: frame from era {f_era} in an era-{s_era} "
                       "series")
    return len(bad), bad[:20]


def verify_no_target_transform(rows: Sequence[tuple[str, float, float, float]],
                               tol: float = 1e-9) -> tuple[int, list[str]]:
    """Check that every target magnitude is ``mag - zp`` and nothing else.

    ``rows`` is ``(series_key, cal_mag, mag, zp_tie)``.  A calibrated target
    magnitude must equal ``mag - zp_tie`` EXACTLY (to float tolerance): a
    zero-point shift and no colour transformation.  Any row where the two
    disagree means a colour term was applied somewhere to a cataclysmic
    variable, which is the one operation this project has ruled out —
    CVs are blue, variable, and routinely outside the colour range any
    transformation was calibrated over.
    """
    bad: list[str] = []
    for key, cal, mag, zp in rows:
        if cal is None or mag is None or zp is None:
            continue
        if not math.isfinite(cal) or not math.isfinite(mag):
            continue
        d = float(cal) - (float(mag) - float(zp))
        if abs(d) > tol:
            bad.append(f"{key}: cal_mag - (mag - zp) = {d:+.3e}")
    return len(bad), bad[:20]


# ===========================================================================
# 5.  TASK 4 — forced photometry and upper limits
# ===========================================================================
@dataclass(frozen=True)
class Similarity:
    """A 4-parameter similarity: rotation+scale+translation, plus its fit
    quality.  Four parameters and not six because the frame-to-reference
    mapping across one night of one telescope IS a rotation and a scale; a
    full affine would happily absorb a genuine mismatch into a shear and
    report a small residual while placing the aperture in the wrong place.
    """

    n: int
    a: float      # x' = a*x - b*y + tx
    b: float
    tx: float
    ty: float
    rms_px: float

    @property
    def scale(self) -> float:
        return math.hypot(self.a, self.b)

    @property
    def rotation_deg(self) -> float:
        return math.degrees(math.atan2(self.b, self.a))


def similarity_from_pairs(src: np.ndarray, dst: np.ndarray) -> Similarity:
    """Least-squares similarity mapping ``src`` onto ``dst``.

    ``src`` and ``dst`` are ``(n, 2)`` arrays of matched positions.  Used to
    map the REFERENCE frame's pixel grid onto a given science frame's grid,
    so that the target's reference position can be evaluated on a frame
    where the target itself was never detected.

    The pairs come from the detections that WERE matched on that frame —
    typically 100-400 comparison stars — so this transform is derived from
    the same evidence the frame's photometry already rests on, and its
    residual RMS is a direct, per-frame statement of how well the forced
    position is known.
    """
    s = np.asarray(src, dtype=float)
    d = np.asarray(dst, dtype=float)
    ok = np.isfinite(s).all(axis=1) & np.isfinite(d).all(axis=1)
    s, d = s[ok], d[ok]
    n = s.shape[0]
    nan = float("nan")
    if n < 2:
        return Similarity(n, nan, nan, nan, nan, nan)
    # Linear model:  [x' ; y'] = [[a, -b], [b, a]] [x ; y] + [tx ; ty]
    # Stack it as a 2n x 4 design in (a, b, tx, ty).
    x, y = s[:, 0], s[:, 1]
    xp, yp = d[:, 0], d[:, 1]
    m = np.zeros((2 * n, 4))
    m[0::2, 0] = x
    m[0::2, 1] = -y
    m[0::2, 2] = 1.0
    m[1::2, 0] = y
    m[1::2, 1] = x
    m[1::2, 3] = 1.0
    rhs = np.empty(2 * n)
    rhs[0::2] = xp
    rhs[1::2] = yp
    sol, *_ = np.linalg.lstsq(m, rhs, rcond=None)
    a, b, tx, ty = (float(v) for v in sol)
    pred = m @ sol
    dx = pred[0::2] - xp
    dy = pred[1::2] - yp
    rms = float(np.sqrt(np.mean(dx ** 2 + dy ** 2)))
    return Similarity(n, a, b, tx, ty, rms)


def apply_similarity(t: Similarity, x: float, y: float) -> tuple[float, float]:
    """Map one position through a :class:`Similarity`."""
    return (t.a * x - t.b * y + t.tx, t.b * x + t.a * y + t.ty)


def gnomonic_project(ra_deg: Sequence[float], dec_deg: Sequence[float],
                     ra0_deg: float, dec0_deg: float
                     ) -> tuple[np.ndarray, np.ndarray]:
    """Tangent-plane (gnomonic) projection about ``(ra0, dec0)``, degrees.

    Needed for exactly one case, and it is the case that matters most: EU
    UMa's era-78 block, whose field tie failed with an HTTP error, carries
    no identified target star at all.  Its reference stars DO carry sky
    coordinates, so the target's pixel position can be recovered by fitting
    a plate model to those stars and evaluating it at the target's
    catalogued position.  Fitting in RA/Dec directly would fail near the
    pole and would not be linear anywhere; projecting first makes the plate
    model the linear thing it physically is.
    """
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    ra0, dec0 = math.radians(ra0_deg), math.radians(dec0_deg)
    cosc = (math.sin(dec0) * np.sin(dec)
            + math.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
    with np.errstate(invalid="ignore", divide="ignore"):
        xi = np.where(cosc != 0, np.cos(dec) * np.sin(ra - ra0) / cosc, np.nan)
        eta = np.where(cosc != 0,
                       (math.cos(dec0) * np.sin(dec)
                        - math.sin(dec0) * np.cos(dec) * np.cos(ra - ra0))
                       / cosc, np.nan)
    return np.degrees(xi), np.degrees(eta)


def affine_from_pairs(src: np.ndarray, dst: np.ndarray
                      ) -> tuple[np.ndarray, float]:
    """Least-squares 6-parameter affine mapping ``src`` onto ``dst``.

    Returns ``(matrix, rms)`` where ``matrix`` is ``(2, 3)`` such that
    ``dst = matrix @ [src_x, src_y, 1]``.  Used for the sky-to-pixel plate
    model above, where an affine IS the right model — a tangent plane maps
    to a detector through a scale, a rotation and, if the camera is not
    square to the sky, a shear.
    """
    s = np.asarray(src, dtype=float)
    d = np.asarray(dst, dtype=float)
    ok = np.isfinite(s).all(axis=1) & np.isfinite(d).all(axis=1)
    s, d = s[ok], d[ok]
    if s.shape[0] < 3:
        return np.full((2, 3), np.nan), float("nan")
    a = np.column_stack([s[:, 0], s[:, 1], np.ones(s.shape[0])])
    sol, *_ = np.linalg.lstsq(a, d, rcond=None)      # (3, 2)
    pred = a @ sol
    rms = float(np.sqrt(np.mean(np.sum((pred - d) ** 2, axis=1))))
    return sol.T, rms


@dataclass(frozen=True)
class ForcedMeasurement:
    """One aperture measurement at a position we chose rather than found."""

    flux: float          # background-subtracted aperture sum, ADU
    flux_err: float      # 1 sigma, ADU
    sky: float           # per-pixel sky level, ADU
    sky_rms: float       # per-pixel sky scatter, ADU
    n_pix: float         # pixels in the aperture
    n_sky: int           # pixels used for the sky estimate
    snr: float


def forced_aperture(image: np.ndarray, x: float, y: float, r_ap: float,
                    r_in: float, r_out: float,
                    gain_e_per_adu: Optional[float] = None
                    ) -> ForcedMeasurement:
    """Aperture photometry at a FIXED position, with a local sky annulus.

    Deliberately written out in numpy rather than delegated, because the
    quantity this whole task turns on is not the flux — it is the NOISE.
    An upper limit is a noise measurement wearing a magnitude's clothes, so
    the reader has to be able to see exactly which terms went into it:

    * ``n_pix * sky_rms^2``    — sky shot + read noise inside the aperture;
    * ``n_pix^2 / n_sky * sky_rms^2`` — the uncertainty of the SKY LEVEL
      itself, which is not negligible: the annulus holds a few hundred
      pixels, and subtracting a sky that is itself uncertain moves the whole
      aperture sum;
    * ``flux / gain``          — the source's own shot noise, dropped when
      the gain is unknown or zero (the Andor iKon writes ``EGAIN = 0``).
      For an undetected source this term is ~0 by construction, so dropping
      it costs nothing where it matters.

    The sky level is a MEDIAN over a sigma-clipped annulus.  A mean would be
    pulled by the wings of any neighbour, and in these fields — VV Pup sits
    at galactic latitude +2 — there is always a neighbour.

    Sub-pixel coverage is handled by a soft edge: each pixel contributes the
    fraction of its area inside the aperture, approximated by clipping the
    signed distance from the circle to +/-0.5 px.  Exact enough at r >= 8 px
    (the smallest aperture in this campaign is 8.9 px) and continuous, which
    a hard mask is not — a hard mask makes the measured flux jump as the
    forced position moves by a tenth of a pixel.
    """
    img = np.asarray(image, dtype=float)
    ny, nx = img.shape
    nan = float("nan")
    # The APERTURE must lie wholly inside the frame.  A partially clipped
    # aperture is the most dangerous possible failure for this task: it
    # under-counts the flux AND under-counts the pixels, so the noise it
    # reports is too small and the upper limit derived from it is too DEEP
    # — a fabricated constraint rather than a missing one.  Refusing is the
    # only safe answer.  (The sky annulus is allowed to clip: it is a
    # median over hundreds of pixels and loses only precision.)
    if not (x - r_ap >= 0 and y - r_ap >= 0
            and x + r_ap <= nx - 1 and y + r_ap <= ny - 1):
        return ForcedMeasurement(nan, nan, nan, nan, nan, 0, nan)
    lo_x = max(0, int(math.floor(x - r_out - 2)))
    hi_x = min(nx, int(math.ceil(x + r_out + 2)))
    lo_y = max(0, int(math.floor(y - r_out - 2)))
    hi_y = min(ny, int(math.ceil(y + r_out + 2)))
    nan = float("nan")
    if hi_x - lo_x < 3 or hi_y - lo_y < 3:
        return ForcedMeasurement(nan, nan, nan, nan, nan, 0, nan)
    cut = img[lo_y:hi_y, lo_x:hi_x]
    yy, xx = np.mgrid[lo_y:hi_y, lo_x:hi_x]
    d = np.hypot(xx - x, yy - y)
    finite = np.isfinite(cut)
    # --- local sky, sigma-clipped median of the annulus ---
    ann = finite & (d >= r_in) & (d <= r_out)
    if int(ann.sum()) < 20:
        return ForcedMeasurement(nan, nan, nan, nan, nan, int(ann.sum()), nan)
    vals = cut[ann]
    med = float(np.median(vals))
    mad = median_abs_deviation(vals)
    if math.isfinite(mad) and mad > 0:
        keep = np.abs(vals - med) <= 3.0 * mad
        if int(keep.sum()) >= 20:
            vals = vals[keep]
    sky = float(np.median(vals))
    sky_rms = float(np.std(vals))
    n_sky = int(vals.size)
    # --- aperture with a soft edge ---
    frac = np.clip(r_ap + 0.5 - d, 0.0, 1.0)
    frac = np.where(finite, frac, 0.0)
    n_pix = float(frac.sum())
    if n_pix <= 0:
        return ForcedMeasurement(nan, nan, nan, nan, nan, n_sky, nan)
    flux = float(np.sum(frac * (cut - sky)))
    var = n_pix * sky_rms ** 2 + (n_pix ** 2 / max(1, n_sky)) * sky_rms ** 2
    if gain_e_per_adu and gain_e_per_adu > 0 and flux > 0:
        var += flux / float(gain_e_per_adu)
    err = math.sqrt(var) if var > 0 else nan
    snr = flux / err if err and math.isfinite(err) and err > 0 else nan
    return ForcedMeasurement(flux, err, sky, sky_rms, n_pix, n_sky, snr)


def limit_flux(flux_err: float, k: float = LIMIT_SIGMA) -> float:
    """``k * sigma`` — the flux an undetected source is bounded below.

    Deliberately NOT ``flux + k*sigma``.  Both conventions are in use, and
    the difference matters at the few-tenths-of-a-magnitude level for a
    source sitting on a positive noise excursion.  ``k*sigma`` states the
    sensitivity of the measurement — "a source brighter than this would have
    been seen" — which is a property of the frame and is what a duty-cycle
    statistic needs.  ``flux + k*sigma`` states a Bayesian-flavoured bound
    on this particular realisation and is systematically fainter on half the
    frames by construction.  The chosen convention is stated on every row of
    the product so that nobody has to guess which one they are reading.
    """
    if flux_err is None or not math.isfinite(flux_err) or flux_err <= 0:
        return float("nan")
    return float(k) * float(flux_err)


def limit_magnitude(flux: float, exptime: float, zp: float,
                    offset: float = INST_MAG_OFFSET) -> float:
    """Turn a limiting flux into a magnitude on the series' own scale.

    ``m = -2.5 log10(flux / exptime) + offset - zp``, matching
    ``macro_phot.photometry.instrumental_mag`` and the ensemble convention
    ``mag = inst_mag - zp`` exactly, so a limit and a detection from the
    same frame sit on the same axis and can be plotted together.  If they
    did not, every figure in the faint-limit section would be comparing two
    different magnitude systems and calling it a light curve.
    """
    if (flux is None or exptime is None or zp is None
            or not math.isfinite(flux) or flux <= 0
            or not math.isfinite(exptime) or exptime <= 0
            or not math.isfinite(zp)):
        return float("nan")
    return float(-2.5 * math.log10(flux / exptime) + offset - zp)


def km_survival(values: Sequence[float], censored: Sequence[bool]
                ) -> tuple[np.ndarray, np.ndarray]:
    """Kaplan-Meier estimate of the magnitude distribution with limits.

    In magnitude, "fainter" is "larger", and an upper limit says the true
    magnitude is LARGER than the quoted value.  That is right-censoring in
    the usual survival sense with no relabelling needed, so the standard
    product-limit estimator applies directly.

    Returns ``(m, S)`` where ``S(m)`` is the estimated fraction of epochs at
    which the target was FAINTER than ``m``.  ``S`` is what makes the
    duty-cycle statistics honest: the naive version divides detections by
    detections and therefore reports the duty cycle of the epochs where the
    target was bright enough to be seen, which is a tautology, not a
    measurement.
    """
    v = np.asarray(values, dtype=float)
    c = np.asarray(censored, dtype=bool)
    ok = np.isfinite(v)
    v, c = v[ok], c[ok]
    if v.size == 0:
        return np.array([]), np.array([])
    order = np.argsort(v, kind="mergesort")
    v, c = v[order], c[order]
    n = v.size
    at_risk = n
    surv = 1.0
    ms, ss = [], []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and v[j + 1] == v[i]:
            j += 1
        d = int((~c[i:j + 1]).sum())          # events = actual detections
        if d > 0 and at_risk > 0:
            surv *= (1.0 - d / at_risk)
        ms.append(float(v[i]))
        ss.append(float(surv))
        at_risk -= (j - i + 1)
        i = j + 1
    return np.asarray(ms), np.asarray(ss)


def km_median(values: Sequence[float], censored: Sequence[bool]) -> float:
    """Median magnitude from the Kaplan-Meier curve, or NaN if the curve
    never reaches 0.5 — which is itself a result, and one this data set
    produces: a series where more than half the epochs are upper limits has
    no estimable median, and inventing one from the detections alone is the
    exact bias this task exists to remove."""
    m, s = km_survival(values, censored)
    if m.size == 0:
        return float("nan")
    below = np.nonzero(s <= 0.5)[0]
    if below.size == 0:
        return float("nan")
    return float(m[below[0]])


def state_statistics(detected_mag: Sequence[float],
                     limit_mag: Sequence[float]) -> dict:
    """Duty-cycle statistics computed twice: censored, and limit-aware.

    The pair is the deliverable.  A single "corrected" number would hide the
    size of the correction, and the size of the correction IS the finding —
    it is how much the published statistic would have been wrong by if the
    undetected epochs had simply been dropped.
    """
    det = np.asarray(detected_mag, dtype=float)
    det = det[np.isfinite(det)]
    lim = np.asarray(limit_mag, dtype=float)
    lim = lim[np.isfinite(lim)]
    n_det, n_lim = det.size, lim.size
    n_all = n_det + n_lim
    values = np.concatenate([det, lim]) if n_all else np.array([])
    censored = np.concatenate([np.zeros(n_det, bool), np.ones(n_lim, bool)])
    return {
        "n_detected": int(n_det),
        "n_limit": int(n_lim),
        "n_epochs": int(n_all),
        "detected_fraction_censored": 1.0 if n_det else float("nan"),
        "detected_fraction_true": (n_det / n_all) if n_all else float("nan"),
        "median_censored": float(np.median(det)) if n_det else float("nan"),
        "median_km": km_median(values, censored) if n_all else float("nan"),
        "faint_state_fraction": (n_lim / n_all) if n_all else float("nan"),
        "median_limit": float(np.median(lim)) if n_lim else float("nan"),
        "faintest_detection": float(np.max(det)) if n_det else float("nan"),
    }
