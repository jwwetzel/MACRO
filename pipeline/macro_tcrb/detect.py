"""macro_tcrb.detect — the pre-registered EW-change detection rule (TCRB-A5a).

The rule is stated in ``TCrB_Monitoring/ANALYSIS_STRATEGY.md`` §10,
"Pre-registered EW-change detection rule (TCRB-A5a) — fixed
2026-10-05T02:40Z", and was written before any EW from the rebuilt grism
library existed.  This module is that text as code.  The constants below are
the rule's thresholds; changing one is a departure from the pre-registration
and must be reported as such in the paper, with both results.

Pure functions only — no file or database access — so every decision is
unit-tested on synthetic series (``pipeline/tests/test_tcrb.py``).

The nightly series of ONE grism is passed as three aligned arrays: night
time ``t`` (days), nightly EW ``e`` and its total error ``s`` (statistical
⊕ floor, rule 1).  The other grism, when given, is used only for the
corroboration clause 2(c).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np
from scipy import stats

# ---- the pre-registered thresholds (rule numbers from the strategy) -------
STEP_SIGMA = 3.0          #: 2(a) step significance
PERSIST_SIGMA_SAME = 2.0  #: 2(b) next night agrees with the new level within
PERSIST_SIGMA_DIFF = 3.0  #: 2(b) ... and still differs from the old level by
MAX_GAP_D = 4.0           #: 2 / 2(b) consecutive nights at most this far apart
CORROB_WINDOW_D = 1.0     #: 2(c) other grism within +-1 d of both nights
CORROB_SIGMA = 2.0        #: 2(c) same-sign difference in the other grism
VAR_P = 1e-3              #: 3 season-level variability p-value
VAR_P_INFLATED = 1e-2     #: 3 ... and still p < this with floor x 1.5
FLOOR_INFLATION = 1.5
MAX_FALSE_PER_SEASON = 0.05   #: 4 expected false step calls allowed
N_MC = 10_000                 #: 4 Monte Carlo realisations
RECOVERY_LEVEL = 0.90         #: 4 smallest step recovered in 90% of injections
FLUX_SIGMA = 2.0              #: 5 an event must survive in line flux at 2 sigma


@dataclass
class StepEvent:
    i: int                      # index of the "before" night
    j: int                      # index of the "after" night
    t_i: float
    t_j: float
    delta: float                # E_j - E_i
    sigma: float                # sqrt(s_i^2 + s_j^2)
    nsig: float
    persist_k: Optional[int]    # index of the night that held the new level
    corroboration: str          # 'corroborated' | 'contradicted' | 'single-grism'
    flux_verdict: str = "untested"   # 'accretion' | 'continuum-driven'


@dataclass
class RuleResult:
    threshold_sigma: float
    expected_false: float
    events: list = field(default_factory=list)
    chi2: float = float("nan")
    dof: int = 0
    p_const: float = float("nan")
    p_const_inflated: float = float("nan")
    variable: bool = False
    min_step_recovered: float = float("nan")
    recovery: dict = field(default_factory=dict)


def _sorted(t, e, s):
    t, e, s = (np.asarray(a, dtype=float) for a in (t, e, s))
    o = np.argsort(t)
    return t[o], e[o], s[o]


def step_candidates(t, e, s, threshold: float = STEP_SIGMA) -> list[StepEvent]:
    """Rule 2(a)+(b): significant, persistent steps between consecutive nights.

    Returns events with ``corroboration`` still unset ('single-grism').
    ``t, e, s`` must already be sorted by time.
    """
    out: list[StepEvent] = []
    n = len(t)
    for i in range(n - 1):
        j = i + 1
        if t[j] - t[i] > MAX_GAP_D:
            continue
        sig = float(np.hypot(s[i], s[j]))
        d = float(e[j] - e[i])
        if not (sig > 0 and abs(d) > threshold * sig):
            continue
        # 2(b) persistence, two-sided: the next observed night after j and
        # the previous observed night before i, each within 4 d.
        k, h = j + 1, i - 1
        if k >= n or t[k] - t[j] > MAX_GAP_D:
            continue
        if h < 0 or t[i] - t[h] > MAX_GAP_D:
            continue
        same = abs(e[k] - e[j]) <= PERSIST_SIGMA_SAME * np.hypot(s[k], s[j])
        still = abs(e[k] - e[i]) > PERSIST_SIGMA_DIFF * np.hypot(s[k], s[i])
        before = abs(e[h] - e[i]) <= PERSIST_SIGMA_SAME * np.hypot(s[h], s[i])
        if same and still and before and np.sign(e[k] - e[i]) == np.sign(d):
            out.append(StepEvent(i, j, float(t[i]), float(t[j]), d, sig,
                                 abs(d) / sig, k, "single-grism"))
    return out


def corroborate(ev: StepEvent, t2, e2, s2) -> str:
    """Rule 2(c): does the other grism show the same-sign step at >= 2 sigma?"""
    if t2 is None or len(t2) == 0:
        return "single-grism"
    t2, e2, s2 = (np.asarray(a, dtype=float) for a in (t2, e2, s2))
    a = np.flatnonzero(np.abs(t2 - ev.t_i) <= CORROB_WINDOW_D)
    b = np.flatnonzero(np.abs(t2 - ev.t_j) <= CORROB_WINDOW_D)
    if a.size == 0 or b.size == 0:
        return "single-grism"
    ia = a[np.argmin(np.abs(t2[a] - ev.t_i))]
    ib = b[np.argmin(np.abs(t2[b] - ev.t_j))]
    if ia == ib:
        return "single-grism"
    d2 = e2[ib] - e2[ia]
    s = np.hypot(s2[ia], s2[ib])
    if np.sign(d2) == np.sign(ev.delta) and abs(d2) >= CORROB_SIGMA * s:
        return "corroborated"
    return "contradicted"


def constancy(e, s) -> tuple[float, int, float]:
    """Rule 3: chi^2 about the inverse-variance mean, its dof and p-value."""
    e, s = np.asarray(e, float), np.asarray(s, float)
    if e.size < 2:
        return float("nan"), 0, float("nan")
    w = 1.0 / s ** 2
    mu = np.sum(w * e) / np.sum(w)
    chi2 = float(np.sum(((e - mu) / s) ** 2))
    dof = int(e.size - 1)
    return chi2, dof, float(stats.chi2.sf(chi2, dof))


def expected_false_steps(t, s, threshold: float, n_mc: int = N_MC,
                         seed: int = 20261005) -> float:
    """Rule 4: mean number of rule-2(a,b) calls on a constant series."""
    rng = np.random.default_rng(seed)
    t, s = np.asarray(t, float), np.asarray(s, float)
    tot = 0
    for _ in range(n_mc):
        e = rng.normal(0.0, s)
        tot += len(step_candidates(t, e, s, threshold))
    return tot / n_mc


def calibrated_threshold(t, s, n_mc: int = N_MC) -> tuple[float, float]:
    """Rule 4: 3 sigma, raised in 0.25 sigma steps until false calls <= 0.05."""
    thr = STEP_SIGMA
    while True:
        ef = expected_false_steps(t, s, thr, n_mc)
        if ef <= MAX_FALSE_PER_SEASON or thr >= 10:
            return thr, ef
        thr += 0.25


def min_recoverable_step(t, s, threshold: float, n_mc: int = 400,
                         seed: int = 7) -> float:
    """Rule 4: smallest step (EW units) recovered in >= 90% of injections.

    The step is injected at a uniformly drawn interior position of the real
    sampling (one with a persistence night on each side) on a constant series with the real
    errors; a recovery is a rule-2(a,b) call spanning the injected boundary.
    """
    rng = np.random.default_rng(seed)
    t, s = np.asarray(t, float), np.asarray(s, float)
    n = len(t)
    if n < 4:
        return float("nan")
    med = float(np.median(s))
    for amp in med * np.arange(1.0, 30.0, 0.25):
        hit = 0
        for r in range(n_mc):
            # step between b-1 and b; both sides need a persistence night
            b = int(rng.integers(2, n - 1))
            e = rng.normal(0.0, s) + np.where(np.arange(n) >= b, amp, 0.0)
            ev = step_candidates(t, e, s, threshold)
            hit += any(x.i == b - 1 and x.j == b for x in ev)
        if hit / n_mc >= RECOVERY_LEVEL:
            return float(amp)
    return float("nan")


def apply_rule(t, e, s, other: Optional[tuple] = None,
               n_mc: int = N_MC) -> RuleResult:
    """Apply rules 2-4 to one grism's nightly series.

    ``other`` = (t2, e2, s2) for the other grism (clause 2(c)) or None.
    Rule 5 (line flux) is applied afterwards by :func:`flux_verdict`.
    """
    t, e, s = _sorted(t, e, s)
    thr, ef = calibrated_threshold(t, s, n_mc)
    res = RuleResult(threshold_sigma=thr, expected_false=ef)
    res.events = step_candidates(t, e, s, thr)
    for ev in res.events:
        ev.corroboration = corroborate(ev, *(other or (None, None, None)))
    res.chi2, res.dof, res.p_const = constancy(e, s)
    _, _, res.p_const_inflated = constancy(e, s * FLOOR_INFLATION)
    res.variable = bool(res.p_const < VAR_P
                        and res.p_const_inflated < VAR_P_INFLATED)
    prof = recovery_profile(t, s, thr, np.median(s) * np.arange(
        1.0, 20.5, 0.5), n_mc=200)
    res.min_step_recovered = prof["amp90_eligible"]
    res.recovery = prof
    return res


def flux_verdict(f_i: float, sf_i: float, f_j: float, sf_j: float) -> str:
    """Rule 5: is the step still there in line flux at >= 2 sigma?"""
    sig = float(np.hypot(sf_i, sf_j))
    if not np.isfinite(sig) or sig <= 0:
        return "untested"
    return ("accretion" if abs(f_j - f_i) >= FLUX_SIGMA * sig
            else "continuum-driven")


def eligible_boundaries(t) -> np.ndarray:
    """Boundaries b (step between nights b-1 and b) at which rule 2 CAN
    fire: the pair and both persistence nights lie within MAX_GAP_D."""
    t = np.asarray(t, float)
    n = len(t)
    return np.array([b for b in range(2, n - 1)
                     if t[b] - t[b - 1] <= MAX_GAP_D
                     and t[b + 1] - t[b] <= MAX_GAP_D
                     and t[b - 1] - t[b - 2] <= MAX_GAP_D], int)


def recovery_profile(t, s, threshold: float, amps, n_mc: int = 300,
                     seed: int = 7) -> dict:
    """Rule-2 recovery of injected steps at the real sampling.

    Returns the overall recovery (step placed at any interior boundary)
    and the recovery at eligible boundaries, per amplitude, plus the
    smallest amplitude recovered in >= 90% of eligible-boundary
    injections.  The overall curve plateaus at the eligible fraction: a
    step across a gap > 4 d cannot be called by the pre-registered rule.
    """
    rng = np.random.default_rng(seed)
    t, s = np.asarray(t, float), np.asarray(s, float)
    n = len(t)
    elig = eligible_boundaries(t)
    out = {"eligible_fraction": len(elig) / max(n - 3, 1), "amps": [],
           "overall": [], "eligible": [], "amp90_eligible": float("nan")}
    for amp in amps:
        h_all = h_el = 0
        for _ in range(n_mc):
            b = int(rng.integers(2, n - 1))
            e = rng.normal(0.0, s) + np.where(np.arange(n) >= b, amp, 0.0)
            h_all += any(x.i == b - 1 and x.j == b
                         for x in step_candidates(t, e, s, threshold))
            if len(elig):
                b2 = int(rng.choice(elig))
                e2 = rng.normal(0.0, s) + np.where(np.arange(n) >= b2, amp,
                                                   0.0)
                h_el += any(x.i == b2 - 1 and x.j == b2
                            for x in step_candidates(t, e2, s, threshold))
        out["amps"].append(float(amp))
        out["overall"].append(h_all / n_mc)
        out["eligible"].append(h_el / n_mc)
        if np.isnan(out["amp90_eligible"]) and h_el / n_mc >= RECOVERY_LEVEL:
            out["amp90_eligible"] = float(amp)
    return out
