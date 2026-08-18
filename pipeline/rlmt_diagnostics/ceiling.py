"""Pure S2 ceiling logic: clip/pileup detection and the per-mode ceiling.

Everything here is a pure function on plain numpy arrays / scalars — no I/O,
no database, no globals.  The campaign script feeds these functions pixel
histograms accumulated from science frames; the unit tests feed them
hand-built histograms including the shapes that MUST NOT read as a clip
(a smooth sky histogram, an empty histogram, a single hot pixel).

THE PHYSICS BEING DETECTED
--------------------------
A detector that runs out of dynamic range does not fade out — it *piles up*:
photon arrivals above the ceiling all land at the top of the scale, so the
pixel histogram of a star field shows a falling bright tail, then a dip,
then a TERMINAL MOUND at the very end of the occupied range.  On a CCD-like
single ADC the mound collapses to one spike (Mode0/Fast pile onto exactly
65,535); on the GSENSE4040 CMOS every pixel clips at a slightly different
corrected code, so High Gain shows a mound a few tens of ADU wide (measured
peak near 3,496, frame maxima spread 3,487-3,584).  Both are the same
fingerprint: (1) the top of the occupied range holds a local maximum that
towers over the valley below it, and (2) essentially nothing lies above it.
A smooth falling tail (no saturation in the sampled frames), a lone hot
pixel, and a steep SKY slope (the 5 MHz iKon frames' sky peak — the false
positive this detector was redesigned to reject) all fail one of the two
gates; we then report the observed maximum but explicitly decline to call
it a ceiling (an honest "not seen" instead of a fabricated number).
"""

from __future__ import annotations

import math
from typing import Optional, Sequence

import numpy as np

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: The terminal mound's peak must stand at least this many times above the
#: valley below it.  The measured High Gain mound peaks ~15x over its
#: valley and the single-code clips (Mode0's 65,535) tower thousands-fold;
#: a monotonic falling tail scores ~1.  10x separates them with margin.
PILEUP_MIN_RATIO = 10.0

#: How far above the top of the DENSE range the mound peak may sit, i.e.
#: the window (in ADU bins) searched for the terminal local maximum.
PILEUP_WINDOW = 160

#: How far below the mound peak the valley is searched for.  The High Gain
#: mound rises out of its dip within ~180 ADU; 300 covers it with room
#: while staying local enough not to reach back into the sky histogram.
PILEUP_VALLEY_SPAN = 300

#: A bin only counts as part of the dense range (and a mound peak is only
#: credible) with at least this many pixels across the sampled frames.
#: Guards against a lone hot pixel masquerading as a clip.
PILEUP_MIN_PIXELS = 1000

#: Candidate clips below this ADU are ignored: a real ceiling sits in the
#: bright range, while bias/overscan artifacts can spike near zero.
PILEUP_MIN_ADU = 500

#: Saturation-veto policy: veto everything above this fraction of the clip,
#: rounded DOWN to :data:`VETO_GRANULARITY_ADU`.  0.92 encodes where the SN
#: panel saw the response leave linearity (approach shoulder starts a few
#: percent below the hard clip); on the measured High Gain clip it lands on
#: the same 3,200 ADU threshold the CV team adopted independently — the two
#: prior choices reconcile under one rule.
VETO_FRACTION = 0.92
VETO_GRANULARITY_ADU = 100

# --------------------------------------------------------------------------
# External priors (documented constants, NOT measurements).  The report may
# quote these as context and must COMPARE live measured numbers against them
# with :func:`prior_comparison` — never assert agreement in fixed prose,
# which would silently go stale if a re-run moved the measurement.
# --------------------------------------------------------------------------

#: The SN panel's High Gain clip range (ADU), from their 2023 supernova
#: photometry saturation study (facility notes; June-2023 detector-test era).
PRIOR_SN_CLIP_LO_ADU = 3530
PRIOR_SN_CLIP_HI_ADU = 3550

#: The CV team's adopted High Gain clip and saturation veto (ADU), from
#: their cataclysmic-variable pipeline configuration (Oct-2023-onward era).
PRIOR_CV_CLIP_ADU = 3500
PRIOR_CV_VETO_ADU = 3200


def prior_comparison(measured: int, prior: int) -> str:
    """Phrase comparing a measured number against an external prior.

    Returns ``"reproduces it exactly"`` only when the two integers are
    EQUAL; otherwise states the signed discrepancy.  The report interpolates
    this instead of hand-asserting agreement, so a re-run that shifts the
    measurement automatically rewrites the sentence instead of lying.
    """
    if measured == prior:
        return "reproduces it exactly"
    return f"differs from it by {measured - prior:+,d} ADU"


def find_clip(hist: np.ndarray) -> Optional[dict]:
    """Locate the clip/pileup spike in a full-pixel ADU histogram.

    Parameters
    ----------
    hist
        Integer array where ``hist[a]`` = number of sampled pixels with
        value exactly ``a`` ADU (index = ADU code, accumulated over frames).

    Returns
    -------
    dict or None
        ``None`` when the terminal structure fails the pileup test (mode
        not observed to saturate).  Otherwise a dict with:

        * ``clip_adu``     — the ADU code of the mound peak (the adopted
          per-mode ceiling; for a single-ADC clip this IS the clip code);
        * ``spike_count``  — pixels sitting exactly at that code;
        * ``tail_level``   — the valley level the mound towers over;
        * ``ratio``        — peak/valley ratio (>= :data:`PILEUP_MIN_RATIO`).

    Two gates, both required:

    1.  **Valley-then-mound**: the top :data:`PILEUP_WINDOW` bins of the
        dense range (bins holding >= :data:`PILEUP_MIN_PIXELS` pixels)
        contain a local maximum standing :data:`PILEUP_MIN_RATIO` times
        above the minimum within :data:`PILEUP_VALLEY_SPAN` bins below it.
        A monotonically falling tail has no such structure (its window
        maximum sits at the window edge with valley ratio ~1).
    2.  **Nothing above**: the total pixel count ABOVE the dense range must
        not exceed the mound peak's own bin count.  A sky slope truncated
        by the density threshold fails here (its distribution continues
        above); a real ceiling has only stray hot-pixel counts above.
    """
    h = np.asarray(hist, dtype=np.int64)
    if h.size == 0 or h.sum() == 0:
        return None
    dense = np.flatnonzero(h >= PILEUP_MIN_PIXELS)
    if dense.size == 0:
        return None
    a_top = int(dense[-1])                     # top of the dense range
    if a_top < PILEUP_MIN_ADU:
        return None                            # bright range never reached
    # Gate 1: valley-then-mound in the terminal window.
    win_lo = max(0, a_top - PILEUP_WINDOW + 1)
    p = win_lo + int(np.argmax(h[win_lo:a_top + 1]))     # mound peak
    val_lo = max(0, p - PILEUP_VALLEY_SPAN)
    valley = float(h[val_lo:p + 1].min())
    ratio = h[p] / max(valley, 1.0)
    if ratio < PILEUP_MIN_RATIO:
        return None
    # Gate 2: essentially nothing may live above the dense range.
    if int(h[a_top + 1:].sum()) > int(h[p]):
        return None
    return {"clip_adu": int(p), "spike_count": int(h[p]),
            "tail_level": valley, "ratio": float(ratio)}


def veto_threshold(clip_adu: int,
                   fraction: float = VETO_FRACTION,
                   granularity: int = VETO_GRANULARITY_ADU) -> int:
    """The saturation-veto threshold implied by a measured clip.

    ``floor(fraction * clip / granularity) * granularity`` — a round number
    a human can quote in an observing log.  3,526 ADU -> 3,200 ADU.
    """
    return int(math.floor(fraction * clip_adu / granularity) * granularity)


def bit_depth_reading(clip_adu: int) -> dict:
    """The is-it-N-bit CONSISTENCY reading for a measured clip.

    Returns the smallest ADC bit depth whose code range contains the clip
    (``bits``), that ADC's full code range (``adc_full_scale``), and the
    headroom the detector never uses (``unused_codes``).  A clip of 3,526
    needs 12 bits (2^12 = 4,096) and leaves 570 codes unused even though
    the FITS files store 16-bit integers.

    HONESTY LIMIT: this is ``ceil(log2(clip + 1))``, so ANY ceiling in
    (2^(b-1), 2^b] reads "b bits" — the histogram cannot distinguish a
    b-bit digitization clip from analog full-well saturation (or a firmware
    clamp) that happens to land below 2^b.  A pure digitization clip lands
    on ONE shared code (Mode0/Fast at 65,535 do); a per-pixel mound spread
    over tens of ADU (High Gain) points at per-pixel analog saturation
    instead.  Callers must therefore report "b-bit-consistent", not "b-bit
    confirmed" — the vendor ADC readback (October hardware item) owns the
    confirmed verdict.
    """
    bits = max(1, math.ceil(math.log2(clip_adu + 1)))
    full = 2 ** bits
    return {"bits": int(bits), "adc_full_scale": int(full),
            "unused_codes": int(full - 1 - clip_adu)}


def merge_hist(total: np.ndarray, frame_hist: np.ndarray) -> np.ndarray:
    """Accumulate one frame's histogram into a running total.

    Both arrays are indexed by ADU; the shorter one is treated as
    zero-padded.  Returns a NEW array (callers keep functional style; the
    campaign script owns the single mutable accumulator).
    """
    n = max(len(total), len(frame_hist))
    out = np.zeros(n, dtype=np.int64)
    out[:len(total)] += np.asarray(total, dtype=np.int64)
    out[:len(frame_hist)] += np.asarray(frame_hist, dtype=np.int64)
    return out


def frame_top_stats(pixels: np.ndarray, veto_adu: Optional[int]) -> dict:
    """Per-frame bright-end summary for the s2_ceiling_frames table.

    Pure on the pixel array: maximum, its (y, x) position (the hot-pixel-
    vs-ceiling diversity evidence), count at the maximum, the 99.9th
    percentile, and — when a veto threshold is already known — how many
    pixels that veto would remove.
    """
    p = np.asarray(pixels)
    mx = int(p.max())
    my, mxx = np.unravel_index(int(np.argmax(p)), p.shape)
    return {
        "max_adu": mx,
        "max_y": int(my), "max_x": int(mxx),
        "n_at_max": int((p == mx).sum()),
        "p999_adu": float(np.percentile(p, 99.9)),
        "n_ge_veto": int((p >= veto_adu).sum()) if veto_adu is not None else None,
    }


#: frame_max_cluster: a mode "clusters" when at least this fraction of its
#: frames share (to within :data:`CLUSTER_REL_WINDOW`) one maximum value.
CLUSTER_MIN_FRAC = 0.30
CLUSTER_REL_WINDOW = 0.02


def frame_max_cluster(maxes: Sequence[float] | np.ndarray) -> Optional[dict]:
    """Ceiling estimate from per-frame maxima, for sparse saturation.

    High Gain StackPro saturates so rarely per pixel that the full-pixel
    histogram's terminal mound stays below the density threshold — but the
    *frame maxima* still pile up at the ceiling (every frame containing one
    saturated star maxes out at the same place: measured cluster ~56,000
    ADU ~ 16 x the single-read clip).  If at least
    :data:`CLUSTER_MIN_FRAC` of frames have their maximum within
    +/-:data:`CLUSTER_REL_WINDOW` of the overall median maximum, that
    cluster's median is the ceiling estimate; otherwise None (frame maxima
    that scatter freely are unsaturated star peaks, not a ceiling).
    """
    m = np.asarray(list(maxes), dtype=np.float64)
    m = m[np.isfinite(m) & (m > 0)]
    if m.size < 10:
        return None
    # Find the MODAL value, not the median: unsaturated frames scatter
    # their maxima freely, so the cluster can sit anywhere in the sorted
    # order.  Log-spaced bins of +/-CLUSTER_REL_WINDOW relative width;
    # the fullest bin nominates the cluster, whose own median then defines
    # a refined +/-window that collects the full membership.
    edges = np.geomspace(m.min(), m.max() * (1 + CLUSTER_REL_WINDOW),
                         max(int(np.log(m.max() / m.min())
                                 / np.log(1 + 2 * CLUSTER_REL_WINDOW)) + 2, 2))
    counts, _ = np.histogram(m, bins=edges)
    k = int(np.argmax(counts))
    members = m[(m >= edges[k]) & (m < edges[k + 1])]
    if members.size == 0:
        return None
    center = float(np.median(members))
    inside = m[np.abs(m - center) <= CLUSTER_REL_WINDOW * center]
    frac = inside.size / m.size
    if frac < CLUSTER_MIN_FRAC:
        return None
    med = float(np.median(inside))
    # The cluster members' own spread (MAD-sigma) is the honest uncertainty
    # of a cluster-based ceiling — per-pixel clip variation plus any
    # era-to-era drift folded into the pooled sample.
    mad_adu = float(1.4826 * np.median(np.abs(inside - med)))
    return {"clip_adu": int(round(med)),
            "cluster_frac": float(frac), "n_frames": int(m.size),
            "mad_adu": mad_adu}


#: A frame-max cluster only counts as a ceiling when its members' argmax
#: POSITIONS are diverse: a true ceiling is hit by whatever bright star the
#: frame happens to contain (position varies with pointing), while a stable
#: hot pixel puts every frame's maximum at the same (y, x).  Measured on
#: the archive: Low Gain's fake cluster repeats one position; the StackPro
#: and iKon clusters scatter freely.
DIVERSITY_MIN_FRAC = 0.5

#: Two argmax positions closer than this many pixels are the same feature
#: (tracking jitter moves a star by a few pixels between frames).
DIVERSITY_TOL_PX = 6


def position_diversity(positions: Sequence[tuple[float, float]],
                       tol: float = DIVERSITY_TOL_PX) -> float:
    """Fraction of distinct argmax positions among cluster members.

    Greedy clustering: a position within ``tol`` pixels (Chebyshev) of an
    already-seen representative joins it; otherwise it founds a new one.
    Returns n_representatives / n_positions — 1.0 for fully diverse maxima,
    ~1/n for a single hot pixel repeated everywhere.  Empty input -> 0.0.
    """
    reps: list[tuple[float, float]] = []
    n = 0
    for y, x in positions:
        n += 1
        for ry, rx in reps:
            if abs(y - ry) <= tol and abs(x - rx) <= tol:
                break
        else:
            reps.append((float(y), float(x)))
    return len(reps) / n if n else 0.0


def mode_group(readoutm: Optional[str]) -> str:
    """Canonical mode label for grouping eras: the READOUTM string, with the
    blank 2026 headers folded onto one explicit label (the current camera's
    frames whose READOUTM card was simply not written — S0's finding)."""
    r = (readoutm or "").strip()
    return r if r else "(blank 2026)"
