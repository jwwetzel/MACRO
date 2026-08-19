"""Pure arithmetic for the CATALOGUE TIE: turning a relative ensemble into
publishable natural-system magnitudes.

WHY THIS MODULE EXISTS
----------------------
A Honeycutt ensemble (``macro_phot.ensemble``) is invariant under
``M_i -> M_i + c, ZP_j -> ZP_j - c``.  The solver fixes that gauge freedom
by demanding ``mean(ZP) = 0``, which is arbitrary: it makes the light curves
internally consistent and externally meaningless.  Every magnitude in
``cv_lightcurve`` therefore sits on a gauge whose origin is "whatever the
comparison stars happened to average to".  The CV characterization graded
the strategy's calibration goal NOT SUPPORTED on exactly this: 7 of 14
(target, era) blocks carried no catalogue tie at all.

This module supplies the arithmetic that removes the arbitrariness, and
NOTHING ELSE — no I/O, no database, no network.  The staged script
``pipeline/scripts/run_cv_cattie.py`` does all of that and calls in here for
every number it stores.

THE MODEL, AND THE FOUR RULINGS IT OBEYS
----------------------------------------
For comparison star *s* in one (target, era, filter) block we regress

    delta_s  ==  m_ens(s) - m_cat(s)  =  ZP0 + k * (C_s - C_ref)

where ``m_ens`` is the star's ensemble mean magnitude (the arbitrary gauge),
``m_cat`` its catalogue magnitude in the band we believe the filter matches,
``C_s`` its catalogue colour, and ``C_ref`` the MEDIAN colour of the tie
stars.  Four committee rulings are built into that one line:

1.  **The natural system is the product.**  A published magnitude is
    ``m_nat = m_ens - ZP0``.  That is the telescope's own bandpass, moved
    onto the catalogue's zero point AT COLOUR ``C_ref`` — no colour
    transformation is ever applied to a light-curve point.  ``k`` is
    published as METADATA describing the bandpass mismatch, and anyone who
    wants standard-system magnitudes of a star whose colour they know can
    apply it themselves.  The science targets are CVs: blue, variable, and
    routinely outside the colour range over which any transformation was
    ever calibrated, so transforming them would replace a known-size
    bandpass error with an unknown-size extrapolation error.
2.  **The tie is solved on comparison stars.**  The target contributes
    nothing to the fit; it is excluded upstream from the comparison pool
    for the same reason (a polar's orbital modulation must not set the zero
    point it is measured against) and excluded here again by role.
3.  **Zero point AND a linear colour term, with uncertainties and a stated
    validity range.**  Centring the colour on ``C_ref`` is not cosmetic: it
    decorrelates ZP0 from k, so ``zp_err`` means "how well the zero point is
    known", not "how well it is known if the colour term were exactly
    right".  :func:`colour_range` reports the interval the fit actually
    interpolates over, and :func:`colour_position` says where a given
    target falls relative to it.
4.  **Saturation and blending veto catalogue stars too.**  A comparison star
    that saturates, or that shares its 4-arcsec aperture with a neighbour
    the catalogue resolves and this telescope does not, carries a
    magnitude error of arbitrary sign and size straight into ZP0.
    :func:`clean_mask` is the single gate every tie star passes through.

WHAT "ACHIEVED ACCURACY" MEANS HERE
-----------------------------------
Residual scatter about the fit is not accuracy — the fit was optimised on
those stars.  :func:`holdout_mask` deterministically withholds a quarter of
the eligible stars (plus every star the ensemble already designated a CHECK
star) from the solve; their ``m_nat - m_cat`` residual after the tie, with
the colour term applied to THEM (their colours are known and they are not
the science target, so rule 1 does not bind), is the number compared with
the strategy's 0.01-0.02 mag goal.

Every function here is deterministic and pure, and every one is unit-tested
in ``pipeline/tests/test_cattie.py``.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

# ===========================================================================
# Tunable constants -- the single source of truth.  The report interpolates
# these rather than repeating them, so a change here changes the page too.
# ===========================================================================

#: The catalogues this stage knows how to read.  ``refcat2`` is ATLAS
#: REFCAT2 (Tonry et al. 2018, ApJ 867, 105) via VizieR ``J/ApJ/867/105``;
#: ``gaia_gspc`` is Gaia DR3 standardised synthetic photometry
#: (``gaiadr3.synthetic_photometry_gspc``, Gaia Collaboration/Montegriffo
#: et al. 2023) via the ESA archive's TAP service.
CATALOGUES = ("refcat2", "gaia_gspc")

#: The PRIMARY catalogue -- the one whose tie is applied to the light
#: curves.  ATLAS-REFCAT2, for three reasons that are argued in the report
#: and summarised here:
#:
#: * it is the system the strategy itself names ("nightly REFCAT2 tie ->
#:   PS1 AB to 0.01-0.02 mag", ANALYSIS_STRATEGY.md section 5), so tying to
#:   it grades the goal that was actually set rather than a substitute;
#: * it reaches m ~ 19, about 1.5 mag deeper than Gaia's synthetic
#:   photometry (which needs a BP/RP spectrum, so it stops near G = 17.65).
#:   Depth is tie stars, and tie stars are what the colour term is fitted
#:   from;
#: * it ships its own blend metrology (the R1/R10 contamination radii),
#:   which is exactly the veto ruling 4 demands and which no other
#:   all-sky catalogue provides for free.
#:
#: Gaia synthetic photometry is fetched for EVERY block anyway, because two
#: catalogues on the same stars is the only honest measurement of the
#: systematic floor, and because it is the only one of the two that covers
#: the PS1 y band and the Johnson-Cousins system (see BAND_CANDIDATES).
PRIMARY_CATALOGUE = "refcat2"

#: Cone radius for a catalogue pull, degrees.  The measured reference-star
#: footprints run to 26.0 arcmin from their own field centre, and those
#: centres sit up to 0.7 arcmin from the target the cone is drawn about, so
#: the binding requirement is 26.7 arcmin.  0.55 deg = 33 arcmin covers it
#: with margin for the pointing scatter between eras.
CONE_RADIUS_DEG = 0.55

#: Faint limit of a catalogue pull.  The faintest comparison star in any
#: solved series sits near G = 18.6, so 19.5 keeps every possible match
#: while stopping the VV Pup cone (galactic latitude ~ +2 deg) from
#: returning a quarter of a million rows nobody will use.
CAT_MAG_MAX = 19.5

#: Sky-match tolerance between a reference star's position and a catalogue
#: position, arcsec.  Tighter than the 2.0 arcsec used for Gaia IDENTITY
#: matching in ``macro_phot.gaia``, deliberately: an identity may be
#: probable, a PHOTOMETRIC tie must be certain, and a 1.2 arcsec circle at
#: this plate scale (0.45-0.81 arcsec/px) is 1.5-2.7 pixels.
MATCH_TOL_ARCSEC = 1.2

#: A match is thrown away as AMBIGUOUS unless the second-nearest catalogue
#: source is at least this many times farther than the accepted one.  Two
#: catalogue entries at comparable distance mean we do not know which star
#: we measured, and a coin-flip identity is worse than no tie star.
AMBIGUITY_FACTOR = 2.5

#: Blend veto, part 1 -- the APERTURE.  Photometry uses a 4.0 arcsec radius
#: aperture (``macro_phot.photometry.APERTURE_RADIUS_ARCSEC``).  A catalogue
#: neighbour inside 1.5x that radius is inside the aperture for any seeing
#: this archive has, so its light is in our flux and not in the catalogue's
#: magnitude for the star we think we measured.
BLEND_APERTURE_ARCSEC = 6.0

#: ...but only if it is bright enough to matter.  2.5 mag fainter is a 10%
#: flux contribution -- about 0.10 mag of error, five times the accuracy
#: goal.  Fainter neighbours than that are tolerated and counted.
BLEND_DMAG = 2.5

#: Blend veto, part 2 -- the SKY ANNULUS (8-12 arcsec,
#: ``macro_phot.photometry.SKY_ANNULUS_ARCSEC``).  A neighbour brighter than
#: the star itself sitting in its own background annulus biases the
#: subtracted sky and therefore the star's flux, in the direction of making
#: it look faint.
BLEND_ANNULUS_ARCSEC = 12.0
BLEND_ANNULUS_DMAG = 0.0

#: Saturation veto.  ``cv_detections.saturated`` is already set per
#: measurement against the era's own readout-mode ceiling (High Gain's scale
#: ends at 3,496 ADU; the fraction of its detections that trip the flag is
#: MEASURED at render time from ``cv_detections`` rather than quoted here --
#: an earlier draft of this comment carried a hand-typed 7.2% that the
#: database disagreed with).  A tie star is
#: rejected if ANY of its measurements in the series tripped that flag:
#: tolerance zero, because a star that saturates on the best-seeing frames
#: has a mean magnitude biased faint by an amount nobody can estimate after
#: the fact, and there are always more candidate tie stars than the fit
#: needs.
SAT_FRAC_MAX = 0.0

#: ...and the shoulder below it.  A detector is non-linear before it clips.
#: A star whose MEDIAN peak exceeds this fraction of the applied veto is in
#: the shoulder even when never formally saturated, so it is rejected too.
NEAR_VETO_FRAC = 0.90

#: A block needs at least this many clean tie stars to be tied at all.  Two
#: free parameters (ZP0 and k) plus enough redundancy for a robust fit to
#: mean anything; below it the block is reported UNTIED rather than tied
#: badly.
MIN_TIE_STARS = 8

#: ...and at least this many held-out stars for the independent accuracy
#: check to be quotable.  With fewer, the achieved accuracy is reported but
#: flagged as an estimate from a small sample.
MIN_CHECK_STARS = 4

#: Fraction of eligible tie stars withheld from the solve.  A quarter is the
#: usual compromise: enough held-out stars for their scatter to mean
#: something, few enough that the fit keeps most of its colour baseline.
HOLDOUT_FRACTION = 0.25

#: Salt for the deterministic holdout hash.  Changing it reshuffles the
#: split; it is recorded in the build metadata so a published accuracy
#: number can be reproduced exactly.
HOLDOUT_SALT = "cattie-v1"

#: Robust fit: Huber tuning constant in units of the robust residual scale.
#: 1.345 is the classical choice (95% efficiency on Gaussian data).
HUBER_DELTA = 1.345

#: Iteratively-reweighted-least-squares cap for the Huber fit.
HUBER_MAX_ITER = 50
HUBER_TOL = 1e-9

#: Hard rejection after the Huber fit converges: a star this many robust
#: sigma from the line is not a mis-weighted point, it is a wrong
#: identification or an unrecognised variable, and it is dropped outright.
CLIP_SIGMA = 4.0

#: Floor on the robust residual scale, magnitudes.  One micro-magnitude:
#: no photometric residual below it is a measurement of anything.  Without
#: the floor, a perfectly-fitting sample drives the scale to machine
#: epsilon and the clip pass starts deleting stars for being 1e-16 mag off
#: the line -- a failure mode that only appears on synthetic data, which is
#: exactly where it would go unnoticed until it appeared on real data.
MIN_RESID_SCALE_MAG = 1e-6

#: The strategy's stated calibration goal (ANALYSIS_STRATEGY.md section 5),
#: in magnitudes.  Verdicts are graded against these two numbers and
#: nothing else.
ACCURACY_GOAL_MAG = 0.02
ACCURACY_STRETCH_MAG = 0.01

#: Published V-magnitude ranges of the five science targets, from the
#: AAVSO Variable Star Index (VSX), recorded 2026-08-19.
#:
#: These are the ONLY hand-entered numbers in this module, and they are
#: here for one purpose: an EXTERNAL check on the zero point.  Every other
#: validation in this stage compares the tie against the same catalogue the
#: tie was fitted to, so all of them would survive a sign error or a
#: gauge mistake applied uniformly to a whole block.  A calibrated target
#: that lands outside its own literature range would not.
#:
#: They are a SANITY bound, never a calibration: a CV's brightness depends
#: on its accretion state, the ranges are broad by construction, and V is
#: not any of the bands measured here.  Nothing is fitted to them and
#: nothing is corrected by them.  A target inside its range is evidence
#: that nothing is grossly wrong; a target outside it is a reason to stop.
TARGET_V_RANGE: dict[str, tuple[float, float, str]] = {
    "anuma": (14.9, 20.0, "AAVSO VSX, AN UMa (AM Her type)"),
    "vvpup": (14.5, 18.0, "AAVSO VSX, VV Pup (AM Her type)"),
    "stlmi": (15.0, 18.0, "AAVSO VSX, ST LMi (AM Her type)"),
    "euuma": (15.6, 19.5, "AAVSO VSX, EU UMa (AM Her type)"),
    "yzcnc": (10.3, 14.9, "AAVSO VSX, YZ Cnc (SU UMa type)"),
}

#: How far outside its literature range a calibrated target may sit before
#: the check is called a FAILURE, magnitudes.  Generous, because the bands
#: differ from V and a CV can be caught outside a catalogued range: this
#: test is meant to catch a wrong sign or a mis-scaled gauge (which would
#: be wrong by magnitudes), not a tenth.
LITERATURE_TOLERANCE_MAG = 1.5


def literature_check(target_key: str, median_mag: float) -> tuple[str, str]:
    """Does a calibrated target land where the literature says it lives?

    Returns ``(verdict, explanation)``.  ``verdict`` is one of ``inside``,
    ``near`` (within LITERATURE_TOLERANCE_MAG of the range), ``outside``,
    or ``unknown`` (no published range recorded for this target).
    """
    entry = TARGET_V_RANGE.get(target_key)
    if entry is None or not np.isfinite(median_mag):
        return "unknown", "no published range recorded for this target"
    lo, hi, src = entry
    if lo <= median_mag <= hi:
        return "inside", f"{lo:.1f}-{hi:.1f} ({src})"
    gap = (lo - median_mag) if median_mag < lo else (median_mag - hi)
    word = "near" if gap <= LITERATURE_TOLERANCE_MAG else "outside"
    return word, (f"{gap:.2f} mag {'below' if median_mag < lo else 'above'} "
                  f"the published {lo:.1f}-{hi:.1f} ({src})")


#: A trend of residual against magnitude or position is called SIGNIFICANT
#: at this many times its own standard error.  3 sigma on a single declared
#: test; the report states how many tests were run so the reader can apply
#: their own look-elsewhere penalty.
TREND_SIGMA = 3.0


# ===========================================================================
# Band mapping -- which catalogue column does each FILTER label mean?
# ===========================================================================

@dataclass(frozen=True)
class BandSpec:
    """One hypothesis for "what standard band is this filter label?".

    Attributes
    ----------
    catalogue
        ``'refcat2'`` or ``'gaia_gspc'``.
    mag_col, err_col
        Column names in the cached catalogue table.  ``err_col`` may be
        None for a catalogue that publishes no per-band error.
    system
        Human-readable photometric system, quoted verbatim in the report
        and stored in ``cv_cattie.band_system``.  This is the string a
        reader needs in order to know what an AB or Vega magnitude here
        means.
    colour_blue, colour_red
        The two catalogue columns whose difference is the colour index the
        colour term is fitted against.
    colour_label
        e.g. ``'g-r'`` -- for axis labels and for the stored metadata.
    hypothesis
        ``'primary'`` for the mapping we believe, ``'alternative'`` for one
        the data is asked to rule out (see the uppercase-label question
        below).
    """
    catalogue: str
    mag_col: str
    err_col: Optional[str]
    system: str
    colour_blue: str
    colour_red: str
    colour_label: str
    hypothesis: str = "primary"


# The RLMT filter wheel, as the archive's own FILTER headers describe it,
# holds g, r, i, z and y -- a Pan-STARRS-like set.  Two eras write those
# labels in UPPERCASE ('G', 'R', 'I' in eras 6 and 7) and the rest in
# lowercase ('g', 'r', 'i' in eras 72/76/78/80; 'r', 'y', 'z' in era 47).
# The obvious reading is that the case difference is a control-software
# convention and the glass is the same.  The obvious reading is not
# evidence, and mis-assigning a band puts a systematic straight into the
# zero point -- so for every uppercase label this table ALSO carries the
# Johnson-Cousins hypothesis, and ``solve`` fits both and lets the colour
# term and the residual scatter decide.  Gaia's synthetic photometry is
# what makes the test possible: it publishes both systems, for the same
# stars, from the same BP/RP spectra, so the two hypotheses differ ONLY in
# the band and not in the sample, the epoch or the calibration.
BAND_CANDIDATES: dict[str, tuple[BandSpec, ...]] = {
    # ---- lowercase: Sloan-like, both catalogues, no ambiguity ------------
    "g": (
        BandSpec("refcat2", "gmag", "e_gmag", "PS1 AB (ATLAS-REFCAT2)",
                 "gmag", "rmag", "g-r"),
        BandSpec("gaia_gspc", "g_sdss_mag", "g_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r"),
    ),
    "r": (
        BandSpec("refcat2", "rmag", "e_rmag", "PS1 AB (ATLAS-REFCAT2)",
                 "gmag", "rmag", "g-r"),
        BandSpec("gaia_gspc", "r_sdss_mag", "r_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r"),
    ),
    "i": (
        BandSpec("refcat2", "imag", "e_imag", "PS1 AB (ATLAS-REFCAT2)",
                 "rmag", "imag", "r-i"),
        BandSpec("gaia_gspc", "i_sdss_mag", "i_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "r_sdss_mag", "i_sdss_mag", "r-i"),
    ),
    "z": (
        BandSpec("refcat2", "zmag", "e_zmag", "PS1 AB (ATLAS-REFCAT2)",
                 "imag", "zmag", "i-z"),
        BandSpec("gaia_gspc", "z_sdss_mag", "z_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "i_sdss_mag", "z_sdss_mag", "i-z"),
    ),
    # ---- y: ONE catalogue can serve it, and it is not the primary one ----
    # ATLAS-REFCAT2 stops at z.  Gaia's synthetic photometry publishes
    # y_ps1, so the block is tieable -- through the SECONDARY catalogue,
    # which the report has to say out loud rather than quietly substituting.
    "y": (
        BandSpec("gaia_gspc", "y_ps1_mag", "y_ps1_mag_error",
                 "PS1 AB (Gaia DR3 synthetic)",
                 "i_sdss_mag", "z_sdss_mag", "i-z"),
    ),
    # ---- uppercase: the same glass, on the primary hypothesis ------------
    "G": (
        BandSpec("refcat2", "gmag", "e_gmag", "PS1 AB (ATLAS-REFCAT2)",
                 "gmag", "rmag", "g-r"),
        BandSpec("gaia_gspc", "g_sdss_mag", "g_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r"),
        BandSpec("gaia_gspc", "v_jkc_mag", "v_jkc_mag_error",
                 "Johnson V Vega (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r",
                 hypothesis="alternative"),
        BandSpec("gaia_gspc", "b_jkc_mag", "b_jkc_mag_error",
                 "Johnson B Vega (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r",
                 hypothesis="alternative"),
    ),
    "R": (
        BandSpec("refcat2", "rmag", "e_rmag", "PS1 AB (ATLAS-REFCAT2)",
                 "gmag", "rmag", "g-r"),
        BandSpec("gaia_gspc", "r_sdss_mag", "r_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r"),
        BandSpec("gaia_gspc", "r_jkc_mag", "r_jkc_mag_error",
                 "Cousins R Vega (Gaia DR3 synthetic)",
                 "g_sdss_mag", "r_sdss_mag", "g-r",
                 hypothesis="alternative"),
    ),
    "I": (
        BandSpec("refcat2", "imag", "e_imag", "PS1 AB (ATLAS-REFCAT2)",
                 "rmag", "imag", "r-i"),
        BandSpec("gaia_gspc", "i_sdss_mag", "i_sdss_mag_error",
                 "SDSS AB (Gaia DR3 synthetic)",
                 "r_sdss_mag", "i_sdss_mag", "r-i"),
        BandSpec("gaia_gspc", "i_jkc_mag", "i_jkc_mag_error",
                 "Cousins I Vega (Gaia DR3 synthetic)",
                 "r_sdss_mag", "i_sdss_mag", "r-i",
                 hypothesis="alternative"),
    ),
}


def band_candidates(filter_label: str,
                    catalogue: Optional[str] = None,
                    hypothesis: Optional[str] = None) -> tuple[BandSpec, ...]:
    """Every band hypothesis for one FILTER label, optionally filtered.

    Returns an empty tuple for a filter this stage cannot tie -- which is a
    RESULT, not an error: "no catalogue analogue" is one of the honest
    reasons a block is left relative.
    """
    out = BAND_CANDIDATES.get(filter_label, ())
    if catalogue is not None:
        out = tuple(b for b in out if b.catalogue == catalogue)
    if hypothesis is not None:
        out = tuple(b for b in out if b.hypothesis == hypothesis)
    return out


def primary_band(filter_label: str, catalogue: str) -> Optional[BandSpec]:
    """The believed band mapping for (filter, catalogue), or None."""
    cands = band_candidates(filter_label, catalogue, hypothesis="primary")
    return cands[0] if cands else None


# ===========================================================================
# Deterministic train / holdout split
# ===========================================================================

def holdout_mask(star_ids: Sequence[int], series_key: str,
                 fraction: float = HOLDOUT_FRACTION,
                 salt: str = HOLDOUT_SALT) -> np.ndarray:
    """True where a star is WITHHELD from the solve, deterministically.

    A random split would make the published accuracy number unreproducible;
    a split on ``star_id % 4`` would correlate with detection order, which
    correlates with brightness, which is the axis the check is meant to be
    blind to.  Hashing ``salt|series_key|star_id`` gives a split that is
    stable across runs, uncorrelated with anything physical, and different
    per series (so a star unlucky in one filter is not systematically
    unlucky in all three).
    """
    out = np.zeros(len(star_ids), dtype=bool)
    cut = int(round(float(fraction) * 2 ** 32))
    for n, sid in enumerate(star_ids):
        key = f"{salt}|{series_key}|{int(sid)}".encode()
        h = int.from_bytes(hashlib.sha256(key).digest()[:4], "big")
        out[n] = h < cut
    return out


# ===========================================================================
# The cleanliness gate -- ruling 4
# ===========================================================================

@dataclass(frozen=True)
class CleanCensus:
    """Why each candidate tie star was kept or thrown away.

    Every count here appears in the report; together they must add up to
    ``n_candidates``, which is the arithmetic that proves no star was
    silently lost.
    """
    n_candidates: int
    n_saturated: int
    n_near_veto: int
    n_blend_aperture: int
    n_blend_annulus: int
    n_ambiguous: int
    n_no_cat_mag: int
    n_flagged: int
    n_clean: int


def clean_mask(*, saturated_frac: np.ndarray, near_veto_frac: np.ndarray,
               blend_sep_arcsec: np.ndarray, blend_dmag: np.ndarray,
               annulus_sep_arcsec: np.ndarray, annulus_dmag: np.ndarray,
               second_sep_arcsec: np.ndarray, match_sep_arcsec: np.ndarray,
               cat_mag: np.ndarray, cat_colour: np.ndarray,
               cat_flag: np.ndarray,
               ) -> tuple[np.ndarray, CleanCensus]:
    """The single gate every candidate tie star passes through.

    All inputs are parallel arrays over candidate stars.  ``blend_*``
    describe the NEAREST catalogue neighbour (separation, and how much
    BRIGHTER the star is than it: positive dmag means the neighbour is
    fainter); ``annulus_*`` describe the nearest neighbour that is brighter
    than the star itself.  ``second_sep_arcsec`` is the distance to the
    second-nearest catalogue source from the star's own position, used for
    the ambiguity test.  ``cat_flag`` is non-zero where the catalogue
    itself says the photometry is suspect (REFCAT2's ``dupvar``, Gaia's
    per-band validity flag).

    Returns the boolean keep-mask and the census that explains it.
    Order matters only for the census (each star is charged to the FIRST
    reason it failed), never for the mask.
    """
    n = len(saturated_frac)
    ok = np.ones(n, dtype=bool)
    charged = np.zeros(n, dtype=bool)
    counts: dict[str, int] = {}

    def _reject(bad, name: str) -> None:
        # ``first`` is the subset not already charged to an earlier reason,
        # so the census partitions the candidates exactly once each; ``ok``
        # is cleared for the WHOLE failing set, so a star that fails two
        # tests is still rejected by both.
        bad = np.asarray(bad, dtype=bool)
        first = bad & ~charged
        counts[name] = int(first.sum())
        charged[first] = True
        ok[bad] = False

    # A missing catalogue magnitude or colour is not a veto, it is an
    # absence -- charged first so the veto counts below describe stars that
    # actually had photometry to lose.
    _reject(~np.isfinite(cat_mag) | ~np.isfinite(cat_colour), "no_cat_mag")
    _reject(np.asarray(cat_flag) != 0, "flagged")
    _reject(np.asarray(saturated_frac) > SAT_FRAC_MAX, "saturated")
    _reject(np.asarray(near_veto_frac) > 0.0, "near_veto")
    _reject((np.asarray(blend_sep_arcsec) < BLEND_APERTURE_ARCSEC)
            & (np.asarray(blend_dmag) < BLEND_DMAG), "blend_aperture")
    _reject((np.asarray(annulus_sep_arcsec) < BLEND_ANNULUS_ARCSEC)
            & (np.asarray(annulus_dmag) <= BLEND_ANNULUS_DMAG),
            "blend_annulus")
    _reject(np.asarray(second_sep_arcsec)
            < AMBIGUITY_FACTOR * np.maximum(np.asarray(match_sep_arcsec),
                                            1e-6), "ambiguous")

    census = CleanCensus(
        n_candidates=n,
        n_saturated=counts.get("saturated", 0),
        n_near_veto=counts.get("near_veto", 0),
        n_blend_aperture=counts.get("blend_aperture", 0),
        n_blend_annulus=counts.get("blend_annulus", 0),
        n_ambiguous=counts.get("ambiguous", 0),
        n_no_cat_mag=counts.get("no_cat_mag", 0),
        n_flagged=counts.get("flagged", 0),
        n_clean=int(ok.sum()))
    return ok, census


def neighbour_metrics(ra: np.ndarray, dec: np.ndarray, mag: np.ndarray,
                      cat_ra: np.ndarray, cat_dec: np.ndarray,
                      cat_mag: np.ndarray,
                      self_index: Optional[np.ndarray] = None,
                      aperture_arcsec: float = BLEND_APERTURE_ARCSEC) -> dict:
    """Blend metrology for each star against the full catalogue.

    Small-angle geometry (these fields are half a degree across, so the
    cos(dec) factor is the only spherical correction that matters).
    Returns arrays of

    ``nn_sep``     separation to the nearest OTHER catalogue source, arcsec
    ``nn_dmag``    that neighbour's magnitude minus the star's (positive =
                   the neighbour is fainter, i.e. harmless)
    ``aper_sep``   separation of the WORST contaminant inside
                   ``aperture_arcsec`` -- the one whose magnitude is closest
                   to (or brighter than) the star's -- or inf if the
                   aperture is empty
    ``aper_dmag``  that contaminant's dmag: the MINIMUM dmag inside the
                   aperture, which is the quantity ruling 4 is about
    ``bright_sep`` separation to the nearest neighbour BRIGHTER than the
                   star (inf if none)
    ``bright_dmag`` its dmag (<= 0 by construction; inf if none)

    WHY ``aper_*`` EXISTS AND ``nn_*`` IS NOT ENOUGH.  The first version of
    the veto tested the NEAREST neighbour, which is not the same star as the
    worst one: a faint source at 1 arcsec and an equal-brightness source at
    3 arcsec both sit inside a 4-arcsec aperture, the nearest test reports
    the faint one, and the equal-brightness contaminant -- the one that
    actually corrupts the magnitude -- is never looked at.  Writing the
    regression test for the BAND defect surfaced 55 such stars still inside
    the VV Pup fits after the band was fixed.  Aperture photometry does not
    care which contaminant is closest; it sums them all, and the brightest
    dominates.

    ``self_index`` gives, per star, the catalogue row that IS the star, so
    it can be excluded from its own neighbour search.  Pass None when the
    star positions are not catalogue rows.
    """
    ra = np.asarray(ra, dtype=float)
    dec = np.asarray(dec, dtype=float)
    mag = np.asarray(mag, dtype=float)
    cra = np.asarray(cat_ra, dtype=float)
    cdec = np.asarray(cat_dec, dtype=float)
    cmag = np.asarray(cat_mag, dtype=float)
    n, m = len(ra), len(cra)
    nn_sep = np.full(n, np.inf)
    nn_dmag = np.full(n, np.inf)
    ap_sep = np.full(n, np.inf)
    ap_dmag = np.full(n, np.inf)
    br_sep = np.full(n, np.inf)
    br_dmag = np.full(n, np.inf)
    if n == 0 or m == 0:
        return {"nn_sep": nn_sep, "nn_dmag": nn_dmag,
                "aper_sep": ap_sep, "aper_dmag": ap_dmag,
                "bright_sep": br_sep, "bright_dmag": br_dmag}
    cosd = np.cos(np.radians(np.median(cdec)))
    for k in range(n):
        d = np.hypot((cra - ra[k]) * cosd, cdec - dec[k]) * 3600.0
        if self_index is not None and self_index[k] >= 0:
            d[int(self_index[k])] = np.inf
        j = int(np.argmin(d))
        if np.isfinite(d[j]):
            nn_sep[k] = d[j]
            nn_dmag[k] = cmag[j] - mag[k]
        # The WORST contaminant inside the aperture: minimum dmag, not
        # minimum distance.  A NaN catalogue magnitude cannot be judged, so
        # it is excluded rather than treated as infinitely faint.
        inside = (d < aperture_arcsec) & np.isfinite(cmag)
        if inside.any():
            dm = np.where(inside, cmag - mag[k], np.inf)
            ja = int(np.argmin(dm))
            ap_sep[k] = float(d[ja])
            ap_dmag[k] = float(dm[ja])
        # The nearest neighbour that is BRIGHTER than this star.
        brighter = cmag < mag[k]
        if brighter.any():
            db = np.where(brighter, d, np.inf)
            jb = int(np.argmin(db))
            if np.isfinite(db[jb]):
                br_sep[k] = db[jb]
                br_dmag[k] = cmag[jb] - mag[k]
    return {"nn_sep": nn_sep, "nn_dmag": nn_dmag,
            "aper_sep": ap_sep, "aper_dmag": ap_dmag,
            "bright_sep": br_sep, "bright_dmag": br_dmag}


#: HEALPix level used to turn a set of Gaia source ids into indexed
#: ``source_id`` ranges (see :func:`hpx_ranges`).  Level 9 pixels are
#: 0.013 deg^2: fine enough that the ranges do not drag in a degree of
#: unwanted sky, coarse enough that a 0.55-deg cone collapses to a few
#: dozen runs rather than a few thousand.
GAIA_HPX_LEVEL = 9


def hpx_ranges(source_ids: Sequence[int],
               level: int = GAIA_HPX_LEVEL) -> list[tuple[int, int]]:
    """Turn a set of Gaia source ids into contiguous ``source_id`` ranges.

    A Gaia DR3 ``source_id`` encodes the source's level-12 HEALPix pixel in
    its high bits: ``pixel(level L) = source_id >> (35 + 2*(12 - L))``.  The
    tables are physically partitioned and ordered by that number, so a
    ``BETWEEN`` on ``source_id`` is a range scan over the archive's own
    partitioning.

    Why this exists at all.  Measured on the ESA archive during this build:
    ``source_id IN (400 primary keys)`` against
    ``gaidr3.synthetic_photometry_gspc`` was answered with "canceling
    statement due to statement timeout" -- the planner will not use the
    index for a large IN-list on a 220-million-row partitioned table, so it
    scans it.  Rewriting the same request as an OR of a few dozen
    ``BETWEEN`` ranges asks for exactly the rows, through exactly the
    structure the table was built around.

    Because source ids are HEALPix-ordered, sky-adjacent stars are usually
    id-adjacent and a cone collapses to a handful of runs.  Sources outside
    the cone that happen to fall inside a run come back too and are
    discarded by the caller -- cheap, and it cannot lose a star that should
    have been kept, which is the direction of error that matters.
    """
    shift = 35 + 2 * (12 - int(level))
    sid = np.asarray(sorted({int(s) for s in source_ids}), dtype=np.int64)
    if sid.size == 0:
        return []
    pix = np.unique(sid >> shift)
    runs: list[tuple[int, int]] = []
    start = prev = int(pix[0])
    for p in pix[1:]:
        p = int(p)
        if p == prev + 1:
            prev = p
        else:
            runs.append((start, prev))
            start = prev = p
    runs.append((start, prev))
    return [(a << shift, ((b + 1) << shift) - 1) for a, b in runs]


def match_by_sky(star_ra: np.ndarray, star_dec: np.ndarray,
                 cat_ra: np.ndarray, cat_dec: np.ndarray,
                 tol_arcsec: float = MATCH_TOL_ARCSEC
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nearest-catalogue-source match, with the ambiguity distance kept.

    Returns ``(index, sep_arcsec, second_sep_arcsec)`` per star; index is
    -1 where nothing lies within ``tol_arcsec``.  The second-nearest
    distance is returned rather than used, because the DECISION about
    ambiguity belongs to :func:`clean_mask` where every veto lives together
    and can be counted in one census.

    This is a nearest-neighbour match, not the greedy one-to-one match used
    for frame registration: two reference stars claiming the same catalogue
    source is itself evidence of a blend, and both are then thrown out by
    the ambiguity rule -- which is the correct outcome, and one that a
    one-to-one assignment would hide by silently giving one of them a
    different, wrong catalogue star.
    """
    sra = np.asarray(star_ra, dtype=float)
    sdec = np.asarray(star_dec, dtype=float)
    cra = np.asarray(cat_ra, dtype=float)
    cdec = np.asarray(cat_dec, dtype=float)
    n = len(sra)
    idx = np.full(n, -1, dtype=int)
    sep = np.full(n, np.inf)
    sep2 = np.full(n, np.inf)
    if n == 0 or len(cra) == 0:
        return idx, sep, sep2
    cosd = np.cos(np.radians(np.median(cdec)))
    for k in range(n):
        if not (np.isfinite(sra[k]) and np.isfinite(sdec[k])):
            continue
        d = np.hypot((cra - sra[k]) * cosd, cdec - sdec[k]) * 3600.0
        order = np.argsort(d)
        if d[order[0]] <= tol_arcsec:
            idx[k] = int(order[0])
            sep[k] = float(d[order[0]])
            sep2[k] = float(d[order[1]]) if len(order) > 1 else np.inf
    return idx, sep, sep2


# ===========================================================================
# The astrometric zero point of the block, measured against the catalogue
# ===========================================================================

#: Loose radius, arcsec, inside which a reference star is paired with a
#: catalogue source for the purpose of MEASURING a rigid offset (never for
#: photometry).  It has to be far wider than MATCH_TOL_ARCSEC or an offset
#: bigger than the photometric tolerance -- the only kind worth finding --
#: would be invisible to the very test meant to find it.  15 arcsec is about
#: 30 pixels at this plate scale: wide enough to catch a whole-frame
#: translation error, narrow enough that the paired source is still far more
#: often the right star than a random field star.
ASTROM_LOOSE_TOL_ARCSEC = 15.0

#: A measured offset is APPLIED only if it is at least this large.  Below it
#: the 1.2-arcsec photometric tolerance already absorbs the shift, and
#: "correcting" it would be fitting the match's own noise -- and would make
#: every block's astrometry depend on the catalogue, which is exactly the
#: dependence a photometric tie must not acquire silently.
ASTROM_REFINE_MIN_ARCSEC = 1.0

#: ...and only if the pairing is COHERENT: the scatter of the individual
#: (dRA, dDec) about their median must be small.  A genuine translation
#: error moves every star by the same vector, so the scatter stays at the
#: astrometric noise; a wrong plate solution, a wrong parity or a random
#: field-star pairing produces a scatter comparable to the search radius,
#: and there is then no single offset to remove.
ASTROM_REFINE_MAX_SCATTER_ARCSEC = 2.0

#: ...measured on at least this many stars.  A median of three pairings is
#: not a measurement of anything.
ASTROM_REFINE_MIN_STARS = 8


@dataclass(frozen=True)
class RigidOffset:
    """A whole-block translation between our sky positions and a catalogue.

    ``dra_arcsec`` is ``catalogue - ours`` in RA times cos(dec), so ADDING
    it to our positions moves them onto the catalogue.  ``scatter_arcsec``
    is the robust (MAD-based) scatter of the individual pairings about that
    median -- the number that decides whether a single vector describes the
    error at all.  ``applied`` records the DECISION, and ``reason`` records
    it in words, so a block that was left alone is as legible as one that
    was corrected.
    """
    dra_arcsec: float
    ddec_arcsec: float
    scatter_arcsec: float
    n: int
    applied: bool
    reason: str

    @property
    def size_arcsec(self) -> float:
        """Length of the offset vector, arcsec (NaN if not measured)."""
        return float(math.hypot(self.dra_arcsec, self.ddec_arcsec))


def rigid_offset(star_ra: Sequence[float], star_dec: Sequence[float],
                 cat_ra: Sequence[float], cat_dec: Sequence[float],
                 loose_tol_arcsec: float = ASTROM_LOOSE_TOL_ARCSEC,
                 min_stars: int = ASTROM_REFINE_MIN_STARS,
                 min_offset_arcsec: float = ASTROM_REFINE_MIN_ARCSEC,
                 max_scatter_arcsec: float = ASTROM_REFINE_MAX_SCATTER_ARCSEC,
                 ) -> RigidOffset:
    """Measure a whole-block translation between our positions and a catalogue.

    WHY THIS EXISTS.  A block's sky positions come either from an S1 plate
    solution or from a similarity fit of the reference frame onto a Gaia
    cone.  The similarity fit has four free parameters and two of them are
    translations, so a fit that locked onto a thin set of correspondences
    can be internally consistent -- correct scale, correct rotation, small
    residuals -- and still hand out every position shifted by the same
    vector.  Nothing upstream can see that, because the upstream fit has no
    absolute reference; the catalogue tie is the first stage that does.

    Adversarial review found exactly this: EU UMa era 78's comparison stars
    sat at a rigid -5.20 arcsec in RA*cos(dec) with 0.88 arcsec of scatter,
    so NOT ONE of them matched a catalogue source inside the 1.2-arcsec
    photometric tolerance, and the block was reported untieable for "bad
    astrometry".  The astrometry was not bad; it was displaced.

    The measurement is deliberately conservative.  The offset is APPLIED
    only when it is bigger than the photometric tolerance can absorb, when
    the individual pairings agree on it, and when enough stars contribute.
    Otherwise the positions are left exactly as the upstream stage wrote
    them and ``reason`` says why -- because a stage that quietly bends every
    block's astrometry onto the catalogue it is about to be tied to has
    destroyed the independence the tie depends on.
    """
    sra = np.asarray(star_ra, dtype=float)
    sdec = np.asarray(star_dec, dtype=float)
    cra = np.asarray(cat_ra, dtype=float)
    cdec = np.asarray(cat_dec, dtype=float)
    nan = float("nan")
    if sra.size == 0 or cra.size == 0:
        return RigidOffset(nan, nan, nan, 0, False,
                           "no stars or no catalogue sources to pair")
    # Pair each star with its NEAREST catalogue source inside the loose
    # radius.  This is the same nearest-neighbour rule the photometric match
    # uses, only with a radius chosen to see the thing being measured.
    idx, sep, _ = match_by_sky(sra, sdec, cra, cdec, loose_tol_arcsec)
    hit = idx >= 0
    n = int(hit.sum())
    if n < min_stars:
        return RigidOffset(nan, nan, nan, n, False,
                           f"only {n} star(s) paired within "
                           f"{loose_tol_arcsec:g}\" (need {min_stars})")
    cosd = np.cos(np.radians(sdec[hit]))
    dra = (cra[idx[hit]] - sra[hit]) * cosd * 3600.0
    ddec = (cdec[idx[hit]] - sdec[hit]) * 3600.0
    mra, mdec = float(np.median(dra)), float(np.median(ddec))
    # Robust scatter of the residual vectors about the median offset.  MAD
    # rather than standard deviation: a handful of field-star pairings must
    # not be able to talk this stage out of a real correction, nor into one.
    resid = np.hypot(dra - mra, ddec - mdec)
    scat = float(1.4826 * np.median(resid))
    size = float(math.hypot(mra, mdec))
    # COHERENCE is tested before SIZE, and the order carries meaning: a
    # block whose pairings scatter has no single translation to remove
    # whatever its median happens to be, and saying "the offset was small"
    # about such a block would describe a median of noise as a measurement.
    if not np.isfinite(scat) or scat > max_scatter_arcsec:
        return RigidOffset(mra, mdec, scat, n, False,
                           f"pairings scatter {scat:.2f}\" about the median, "
                           f"above the {max_scatter_arcsec:g}\" coherence "
                           f"limit — no single translation describes this")
    if size < min_offset_arcsec:
        return RigidOffset(mra, mdec, scat, n, False,
                           f"offset {size:.2f}\" is below the "
                           f"{min_offset_arcsec:g}\" threshold; the "
                           f"{MATCH_TOL_ARCSEC:g}\" match tolerance absorbs it")
    return RigidOffset(mra, mdec, scat, n, True,
                       f"rigid offset {size:.2f}\" "
                       f"(dRA*cos(dec) {mra:+.2f}\", dDec {mdec:+.2f}\") "
                       f"measured on {n} stars with {scat:.2f}\" scatter, "
                       f"removed before the photometric match")


def apply_offset(star_ra: Sequence[float], star_dec: Sequence[float],
                 dra_arcsec: float, ddec_arcsec: float
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Move sky positions by a rigid offset (the inverse is a sign flip).

    ``dra_arcsec`` is a great-circle offset, so it is divided by cos(dec) to
    become a coordinate offset in RA.  Small-angle throughout, which is
    exact enough for a few arcsec at these declinations.
    """
    ra = np.asarray(star_ra, dtype=float)
    dec = np.asarray(star_dec, dtype=float)
    cosd = np.cos(np.radians(dec))
    cosd = np.where(np.abs(cosd) < 1e-9, 1e-9, cosd)
    return (ra + float(dra_arcsec) / 3600.0 / cosd,
            dec + float(ddec_arcsec) / 3600.0)


# ===========================================================================
# The robust straight-line fit -- ruling 3
# ===========================================================================

@dataclass(frozen=True)
class LineFit:
    """One robust fit of ``y = zp + slope * (x - x_ref)``.

    ``zp`` is the intercept AT ``x_ref`` (not at x = 0), which is what makes
    ``zp_err`` and ``slope_err`` nearly independent.
    """
    zp: float
    zp_err: float
    slope: float
    slope_err: float
    x_ref: float
    resid_rms: float
    resid_mad: float
    n_used: int
    n_clipped: int
    chi2nu: float
    converged: bool


def _wls(x: np.ndarray, y: np.ndarray, w: np.ndarray
         ) -> tuple[float, float, np.ndarray]:
    """Weighted least squares for ``y = a + b*x``; returns (a, b, cov 2x2).

    Written out rather than delegated to numpy.polyfit so the covariance is
    the ANALYTIC one from the weights, and so a student can check it against
    any textbook: with S = sum(w), Sx = sum(w x), Sxx = sum(w x^2),
    Sy = sum(w y), Sxy = sum(w x y), the normal equations are
    [[S, Sx], [Sx, Sxx]] [a, b]^T = [Sy, Sxy]^T.
    """
    S = float(np.sum(w))
    Sx = float(np.sum(w * x))
    Sxx = float(np.sum(w * x * x))
    Sy = float(np.sum(w * y))
    Sxy = float(np.sum(w * x * y))
    det = S * Sxx - Sx * Sx
    if not np.isfinite(det) or abs(det) < 1e-300:
        return float("nan"), float("nan"), np.full((2, 2), np.nan)
    a = (Sxx * Sy - Sx * Sxy) / det
    b = (S * Sxy - Sx * Sy) / det
    cov = np.array([[Sxx, -Sx], [-Sx, S]], dtype=float) / det
    return a, b, cov


def robust_line_fit(x: np.ndarray, y: np.ndarray,
                    yerr: Optional[np.ndarray] = None,
                    x_ref: Optional[float] = None,
                    huber_delta: float = HUBER_DELTA,
                    clip_sigma: float = CLIP_SIGMA) -> LineFit:
    """Huber-robust weighted fit of ``y = zp + slope*(x - x_ref)``.

    Why robust at all: a comparison-star sample always contains a few stars
    that are variable, mis-identified, or blended with something the
    catalogue resolves and this telescope does not.  Ordinary least squares
    lets any one of them tilt the colour term, and the colour term is a
    published number.  Iteratively-reweighted least squares with Huber
    weights lets an outlier keep a vote proportional to 1/|residual| instead
    of |residual|, which is the difference between "downweighted" and
    "deleted" -- and only after convergence, where the identification is
    plainly wrong rather than merely noisy, does the CLIP_SIGMA pass delete.

    ``x_ref`` defaults to the median of x: centring is what decorrelates the
    intercept from the slope, so that ``zp_err`` answers "how well is the
    zero point known" rather than "how well would it be known if the colour
    term were exactly right".

    The reported ``zp_err``/``slope_err`` are the weighted-least-squares
    errors SCALED by sqrt(chi2nu) when chi2nu > 1.  That is the honest
    choice for photometry: catalogue errors are known to understate the
    real star-to-star scatter (bandpass mismatch beyond a linear colour
    term is a real, per-star effect), and quoting the formal error of an
    obviously-underdispersed model would understate the zero point's
    uncertainty by exactly the factor the fit itself measured.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    # UNIT WEIGHTS mean "we do not know the errors", not "the errors are
    # 1 magnitude".  The distinction decides how the covariance is scaled
    # at the end: with real errors the scaling only ever INFLATES (never
    # let a model claim to be better than its own scatter); with unit
    # weights the residual scatter IS the only error estimate available,
    # so the covariance must be scaled by it in both directions or every
    # trend test comes back insignificant no matter how large it is.
    unit_weights = yerr is None
    if unit_weights:
        sig = np.ones_like(x)
    else:
        sig = np.asarray(yerr, dtype=float).copy()
        sig[~np.isfinite(sig) | (sig <= 0)] = np.nan
        # A star with no usable error is not discarded; it is given the
        # median error, so that "we do not know this one's error" costs it
        # no more and no less than an average star's vote.
        med = np.nanmedian(sig[good]) if good.any() else 1.0
        sig[~np.isfinite(sig)] = med if np.isfinite(med) else 1.0
        good &= np.isfinite(sig)
    nan = float("nan")
    if int(good.sum()) < 3:
        return LineFit(nan, nan, nan, nan,
                       float(x_ref) if x_ref is not None else nan,
                       nan, nan, int(good.sum()), 0, nan, False)
    xr = float(np.median(x[good])) if x_ref is None else float(x_ref)
    xc = x - xr

    keep = good.copy()
    converged = False
    a = b = nan
    cov = np.full((2, 2), nan)
    scale = nan
    for _clip_pass in range(2):
        w = np.where(keep, 1.0 / np.square(np.where(keep, sig, 1.0)), 0.0)
        a, b, cov = _wls(xc[keep], y[keep], w[keep])
        for _ in range(HUBER_MAX_ITER):
            r = y - (a + b * xc)
            # Robust residual scale: MAD of the KEPT residuals, converted to
            # a Gaussian-equivalent sigma.  Using the plain RMS here would
            # let the very outliers we are protecting against inflate the
            # scale until nothing is downweighted.
            scale = 1.4826 * float(np.median(np.abs(r[keep] - np.median(r[keep]))))
            if not np.isfinite(scale) or scale <= 0:
                scale = float(np.std(r[keep]))
            scale = max(scale, MIN_RESID_SCALE_MAG)
            u = np.abs(r) / scale
            hub = np.where(u <= huber_delta, 1.0, huber_delta / np.maximum(u, 1e-12))
            w = np.where(keep, hub / np.square(np.where(keep, sig, 1.0)), 0.0)
            a2, b2, cov = _wls(xc[keep], y[keep], w[keep])
            moved = max(abs(a2 - a), abs(b2 - b))
            a, b = a2, b2
            if moved < HUBER_TOL:
                converged = True
                break
        r = y - (a + b * xc)
        bad = keep & (np.abs(r) > clip_sigma * scale)
        if not bad.any():
            break
        keep &= ~bad

    r = y - (a + b * xc)
    n_used = int(keep.sum())
    resid = r[keep]
    rms = float(np.sqrt(np.mean(np.square(resid)))) if n_used else nan
    mad = float(1.4826 * np.median(np.abs(resid - np.median(resid)))) \
        if n_used else nan
    dof = max(n_used - 2, 1)
    chi2 = float(np.sum(np.square(resid / sig[keep]))) if n_used else nan
    chi2nu = chi2 / dof if np.isfinite(chi2) else nan
    infl = (math.sqrt(chi2nu)
            if np.isfinite(chi2nu) and (unit_weights or chi2nu > 1.0)
            else 1.0)
    return LineFit(zp=float(a), zp_err=float(np.sqrt(cov[0, 0]) * infl),
                   slope=float(b), slope_err=float(np.sqrt(cov[1, 1]) * infl),
                   x_ref=xr, resid_rms=rms, resid_mad=mad,
                   n_used=n_used, n_clipped=int(good.sum()) - n_used,
                   chi2nu=float(chi2nu), converged=converged)


# ===========================================================================
# Colour range, and where a target falls in it
# ===========================================================================

def colour_range(colours: Sequence[float]) -> tuple[float, float, float, float]:
    """``(min, max, p05, p95)`` of the colours the fit interpolated over.

    Both are reported.  The full span says where the fit has ANY leverage;
    the 5th-95th percentile core says where it has enough stars for the
    linear term to be believed.  A target inside the core is calibrated; a
    target inside the span but outside the core is calibrated by a handful
    of stars; a target outside the span is EXTRAPOLATED and the report has
    to say so.
    """
    c = np.asarray([v for v in colours if np.isfinite(v)], dtype=float)
    if c.size == 0:
        nan = float("nan")
        return nan, nan, nan, nan
    return (float(c.min()), float(c.max()),
            float(np.percentile(c, 5)), float(np.percentile(c, 95)))


def colour_position(target_colour: float, cmin: float, cmax: float,
                    p05: float, p95: float) -> str:
    """One word for where a target's colour sits relative to the fit.

    ``'unknown'``      the target's colour was never measured here
    ``'inside-core'``  within the 5th-95th percentile of the tie stars
    ``'inside-span'``  within the full span but outside the core
    ``'extrapolated'`` outside the span entirely
    """
    if not np.isfinite(target_colour) or not np.isfinite(cmin):
        return "unknown"
    if p05 <= target_colour <= p95:
        return "inside-core"
    if cmin <= target_colour <= cmax:
        return "inside-span"
    return "extrapolated"


def colour_extrapolation_error(target_colour: float, cmin: float, cmax: float,
                               slope: float, slope_err: float) -> float:
    """Extra magnitude error a target incurs by sitting outside the range.

    Zero inside the span.  Outside it, the honest charge is the colour term
    times the distance beyond the nearest edge -- the size of the
    correction the linear model WOULD apply out there, which is also the
    size of the error if the model is wrong out there, plus the slope's own
    uncertainty over the same lever arm.  This number is not subtracted from
    anything; it is quoted beside the block so a reader can see what the
    extrapolation would cost if they insisted on it.
    """
    if not np.isfinite(target_colour) or not np.isfinite(cmin):
        return float("nan")
    if cmin <= target_colour <= cmax:
        return 0.0
    reach = (target_colour - cmax) if target_colour > cmax else (cmin - target_colour)
    return float(abs(slope) * reach + abs(slope_err) * reach)


#: Maximum time between the two filter measurements that may be differenced
#: into ONE colour, days.  0.005 d is 7.2 minutes.
#:
#: This constant exists because of a real error found in review.  The first
#: implementation formed the target's colour from its ENSEMBLE MEAN
#: magnitude in each filter -- the mean over every frame of that series --
#: and differenced the two.  For a constant star that is correct.  For a
#: polar it is not a colour at all: it is the difference between the
#: target's mean state during the blue campaign and its mean state during
#: the red campaign, and those campaigns need not overlap.  VV Pup era 76
#: was published at g-r = -1.73 by that recipe; its g points span 370 days
#: and its r points 55 of them, and on the epochs where both filters were
#: actually observed the colour is +0.04.  A 1.77-magnitude error, and it
#: moved the block's stated verdict from "inside-span" to "extrapolated".
#:
#: 7.2 minutes is chosen from the data rather than from taste: these series
#: interleave filters within a night, the answer is stable from 7 to 72
#: minutes of tolerance (it moves by under 20 mmag on every block), and the
#: tighter end costs almost no pairs.  It is still not simultaneity -- these
#: polars have ~100-minute orbits, so 7 minutes is 0.07 in phase -- which is
#: why :func:`paired_colour` also returns the pair-to-pair SCATTER, and why
#: the report quotes it beside every colour.
COLOUR_PAIR_TOL_DAYS = 0.005

#: Fewer paired epochs than this and the colour is reported as UNKNOWN
#: rather than as a number.  A CV colour from one or two pairs is one or two
#: orbital phases of a strongly modulated object, and calling that "the
#: target's colour" is the same mistake as the campaign-mean recipe, only
#: smaller.
MIN_COLOUR_PAIRS = 5


@dataclass(frozen=True)
class PairedColour:
    """A target colour formed only from near-simultaneous measurements."""
    colour: float             #: median of the paired differences
    scatter: float            #: robust pair-to-pair scatter (1.4826 x MAD)
    n_pairs: int
    dt_median_days: float     #: median time separation actually achieved
    note: str


def paired_colour(t_blue: Sequence[float], m_blue: Sequence[float],
                  t_red: Sequence[float], m_red: Sequence[float],
                  tol_days: float = COLOUR_PAIR_TOL_DAYS,
                  min_pairs: int = MIN_COLOUR_PAIRS) -> PairedColour:
    """Colour of a VARIABLE target, from epochs where both filters saw it.

    Each blue point is paired with the nearest red point in time; pairs
    farther apart than ``tol_days`` are dropped, and the colour is the
    MEDIAN of the surviving differences.  A red point may serve more than
    one blue point -- the alternative, a one-to-one assignment, would throw
    away good pairs for a tidiness that buys nothing here.

    Returns NaN with ``n_pairs`` below ``min_pairs`` when the two filters
    never sampled the same time: that is an honest "not measurable", and it
    is the correct answer for a block whose g and r campaigns share no
    night.  Magnitudes in, magnitude out -- the caller passes natural-system
    magnitudes (each filter's own mean already zero-pointed) and receives a
    natural-system colour.
    """
    tb = np.asarray(t_blue, dtype=float)
    mb = np.asarray(m_blue, dtype=float)
    tr = np.asarray(t_red, dtype=float)
    mr = np.asarray(m_red, dtype=float)
    nan = float("nan")
    okb = np.isfinite(tb) & np.isfinite(mb)
    okr = np.isfinite(tr) & np.isfinite(mr)
    tb, mb, tr, mr = tb[okb], mb[okb], tr[okr], mr[okr]
    if tb.size == 0 or tr.size == 0:
        return PairedColour(nan, nan, 0, nan,
                            "one of the two filters has no finite "
                            "measurement of the target in this era")
    order = np.argsort(tr)
    tr, mr = tr[order], mr[order]
    # Nearest red point to each blue point, found by insertion position:
    # the candidate is either the point just before or just after, so one
    # comparison settles it and the whole pairing stays O(n log n).
    if len(tr) == 1:
        j = np.zeros(len(tb), dtype=int)
    else:
        pos = np.clip(np.searchsorted(tr, tb), 1, len(tr) - 1)
        left, right = pos - 1, pos
        j = np.where(np.abs(tb - tr[left]) <= np.abs(tb - tr[right]),
                     left, right)
    dt = np.abs(tb - tr[j])
    keep = dt <= tol_days
    n = int(keep.sum())
    if n < min_pairs:
        return PairedColour(nan, nan, n, float(np.min(dt)) if dt.size else nan,
                            f"only {n} epoch pair(s) within "
                            f"{tol_days * 1440:.0f} min (need {min_pairs}) — "
                            f"the two filters did not sample the same time, "
                            f"so this target has no measured colour in this "
                            f"era")
    d = mb[keep] - mr[j[keep]]
    med = float(np.median(d))
    scat = float(1.4826 * np.median(np.abs(d - med)))
    return PairedColour(med, scat, n, float(np.median(dt[keep])),
                        f"median of {n} epoch pairs within "
                        f"{tol_days * 1440:.0f} min "
                        f"(median separation {np.median(dt[keep]) * 1440:.1f} "
                        f"min); pair-to-pair scatter {scat:.3f} mag")


def target_colour_solve(nat_colour: float, k_blue: float, k_red: float,
                        c_ref: float) -> float:
    """The target's CATALOGUE colour, inferred from its natural-system one.

    This does NOT transform the target's magnitude -- ruling 1 forbids that
    and nothing here writes a magnitude.  It answers a different question:
    "where on the colour axis of the fit does this target sit?", which has
    to be answered before anyone can say whether the tie applies to it.

    The algebra.  For each band, the natural-system magnitude relates to the
    catalogue one by ``m_nat = m_cat + k*(C - C_ref)`` (that is the fit,
    rearranged).  Subtracting the two bands, with ``C = m_cat,blue -
    m_cat,red`` the catalogue colour and ``N = m_nat,blue - m_nat,red`` the
    measured natural-system colour::

        N = C + (k_blue - k_red) * (C - C_ref)
          => C = (N + (k_blue - k_red) * C_ref) / (1 + k_blue - k_red)

    The denominator vanishes only if the two colour terms differ by exactly
    -1, i.e. if the instrument's colour baseline has collapsed to nothing;
    that returns NaN rather than a number, because in that case the
    measurement genuinely carries no colour information.

    NOTE the ``c_ref`` here must be the SAME reference colour for both
    bands, which the solver guarantees by computing it once per
    (target, era, catalogue, colour index) rather than per filter.
    """
    if not (np.isfinite(nat_colour) and np.isfinite(k_blue)
            and np.isfinite(k_red) and np.isfinite(c_ref)):
        return float("nan")
    den = 1.0 + k_blue - k_red
    if abs(den) < 1e-6:
        return float("nan")
    return float((nat_colour + (k_blue - k_red) * c_ref) / den)


# ===========================================================================
# Adversarial checks -- ruling 5
# ===========================================================================

@dataclass(frozen=True)
class Trend:
    """A regression of residual on one nuisance variable."""
    name: str
    slope: float
    slope_err: float
    n: int
    span: float
    significance: float          # |slope| / slope_err
    swing: float                 # slope * span: the trend's size in mag
    significant: bool


def residual_trend(name: str, x: Sequence[float], resid: Sequence[float],
                   sigma: float = TREND_SIGMA) -> Trend:
    """Regress tie residuals on one variable and say whether it is real.

    Used three ways, and each answers a specific accusation:

    * against catalogue MAGNITUDE -- detector non-linearity, or an aperture
      correction that varies with signal-to-noise, would show here;
    * against detector X and Y -- a flat-field residual would show here, and
      the CV characterization flagged a noise floor 28x above scintillation
      that nobody had yet looked at in (x, y);
    * against airmass or seeing, if a caller wants it.

    ``swing`` (slope times the observed span of x) is the number that
    matters for a reader: a formally significant slope of 1 mmag per
    magnitude across a 4-magnitude range is 4 mmag of systematic, which is
    a fifth of the accuracy goal, whereas the same significance across a
    0.2-magnitude range is nothing.
    """
    x = np.asarray(x, dtype=float)
    r = np.asarray(resid, dtype=float)
    ok = np.isfinite(x) & np.isfinite(r)
    nan = float("nan")
    if int(ok.sum()) < 5:
        return Trend(name, nan, nan, int(ok.sum()), nan, nan, nan, False)
    fit = robust_line_fit(x[ok], r[ok])
    span = float(np.ptp(x[ok]))
    sig = (abs(fit.slope) / fit.slope_err
           if np.isfinite(fit.slope_err) and fit.slope_err > 0 else nan)
    return Trend(name=name, slope=fit.slope, slope_err=fit.slope_err,
                 n=int(ok.sum()), span=span, significance=sig,
                 swing=fit.slope * span,
                 significant=bool(np.isfinite(sig) and sig >= sigma))


def plane_fit(x: Sequence[float], y: Sequence[float], z: Sequence[float]
              ) -> dict:
    """Least-squares plane ``z = c0 + cx*x + cy*y`` with errors.

    The flat-field question in its natural form: is there a TILT of the tie
    residual across the detector?  Reported as the two gradients and, more
    usefully, as ``swing`` -- the plane's total change corner to corner,
    which is the systematic a star suffers for being on the wrong side of
    the chip.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    nan = float("nan")
    if int(ok.sum()) < 6:
        return {"c0": nan, "cx": nan, "cy": nan, "cx_err": nan,
                "cy_err": nan, "n": int(ok.sum()), "swing": nan,
                "significance": nan, "significant": False, "rms": nan}
    xs, ys, zs = x[ok], y[ok], z[ok]
    x0, y0 = float(np.mean(xs)), float(np.mean(ys))
    A = np.column_stack([np.ones(len(xs)), xs - x0, ys - y0])
    coef, *_ = np.linalg.lstsq(A, zs, rcond=None)
    resid = zs - A @ coef
    dof = max(len(xs) - 3, 1)
    s2 = float(np.sum(resid ** 2)) / dof
    cov = s2 * np.linalg.inv(A.T @ A)
    cx_err, cy_err = math.sqrt(abs(cov[1, 1])), math.sqrt(abs(cov[2, 2]))
    swing = float(abs(coef[1]) * np.ptp(xs) + abs(coef[2]) * np.ptp(ys))
    sig = max(abs(coef[1]) / cx_err if cx_err > 0 else 0.0,
              abs(coef[2]) / cy_err if cy_err > 0 else 0.0)
    return {"c0": float(coef[0]), "cx": float(coef[1]), "cy": float(coef[2]),
            "cx_err": cx_err, "cy_err": cy_err, "n": int(len(xs)),
            "swing": swing, "significance": float(sig),
            "significant": bool(sig >= TREND_SIGMA),
            "rms": float(np.sqrt(np.mean(resid ** 2)))}


def check_accuracy(nat_mag: Sequence[float], cat_mag: Sequence[float],
                   colour: Sequence[float], zp: float, slope: float,
                   c_ref: float, clip_sigma: float = CLIP_SIGMA) -> dict:
    """Achieved ABSOLUTE accuracy on stars the fit never saw.

    For a HELD-OUT star the colour is known and the star is not the science
    target, so the full transformation may be applied to it -- that is the
    point: it tests the calibration as a user would apply it.  The residual

        (m_nat - slope*(C - C_ref)) - m_cat

    is what a published magnitude of that star would be wrong by.

    TWO numbers come back, and reporting only one of them would be a lie in
    a different direction each time.  Measured on this archive: 91% of
    held-out stars land inside 0.05 mag while two or three per block land
    0.4-1.3 mag out -- variables, blends the veto missed, mis-matches.

    * ``rms`` over ALL held-out stars is dominated by those few.  Quoting it
      as "the accuracy" would say the calibration is ten times worse than
      it is for a typical star.
    * ``rms_clip``, after removing stars beyond ``clip_sigma`` times the
      robust scale, is the accuracy of a typical star -- and it is the
      APPLES-TO-APPLES comparison, because the fit itself clipped at the
      same threshold, so an unclipped check against a clipped fit would be
      measuring the clipping rather than the calibration.
    * ``n_outlier`` / ``outlier_frac`` are the price of the clip, stated
      rather than absorbed: a user calibrating one arbitrary star has
      exactly this probability of drawing one of them.

    ``median`` is the residual zero-point bias -- the part that does NOT
    average away over stars, and therefore the more dangerous half.
    """
    m = np.asarray(nat_mag, dtype=float)
    c = np.asarray(cat_mag, dtype=float)
    col = np.asarray(colour, dtype=float)
    ok = np.isfinite(m) & np.isfinite(c) & np.isfinite(col)
    nan = float("nan")
    if not ok.any():
        return {"n": 0, "rms": nan, "rms_clip": nan, "median": nan,
                "mad": nan, "n_outlier": 0, "outlier_frac": nan,
                "resid": np.array([]), "keep": np.array([], dtype=bool)}
    r = (m[ok] - slope * (col[ok] - c_ref)) - c[ok]
    med = float(np.median(r))
    mad = float(1.4826 * np.median(np.abs(r - med)))
    scale = max(mad, MIN_RESID_SCALE_MAG)
    keep = np.abs(r - med) <= clip_sigma * scale
    n_out = int((~keep).sum())
    return {"n": int(ok.sum()),
            "rms": float(np.sqrt(np.mean(np.square(r)))),
            "rms_clip": (float(np.sqrt(np.mean(np.square(r[keep]))))
                         if keep.any() else nan),
            "median": med, "mad": mad,
            "n_outlier": n_out,
            "outlier_frac": n_out / float(len(r)),
            "resid": r, "keep": keep}


# ===========================================================================
# Verdicts -- one word per block, graded against the strategy's own numbers
# ===========================================================================

def block_verdict(n_tie: int, n_check: int, check_rms: float,
                  fit_converged: bool) -> str:
    """The block's grade, and the ONLY place the grading rule is written.

    ``check_rms`` is the OUTLIER-CLIPPED check RMS (``rms_clip`` from
    :func:`check_accuracy`), because the fit clipped at the same threshold
    and grading a clipped fit with an unclipped check would measure the
    clip.  The unclipped RMS and the outlier fraction are reported beside
    every verdict rather than folded into it.

    ``UNTIED``          fewer than MIN_TIE_STARS clean stars, or no fit.
    ``TIED-STRETCH``    achieved accuracy meets the 0.01 mag stretch goal.
    ``TIED-GOAL``       meets the 0.02 mag goal.
    ``TIED-ABOVE-GOAL`` tied, but the independent check misses the goal --
                        the magnitudes are on a standard system and their
                        accuracy is measured and worse than hoped.  That is
                        a usable result with a stated error, not a failure.
    ``TIED-UNVERIFIED`` tied, but too few held-out stars to measure the
                        accuracy at all.  The weakest passing grade: the
                        zero point exists and nothing independent confirms
                        it.
    """
    if n_tie < MIN_TIE_STARS or not fit_converged:
        return "UNTIED"
    if n_check < MIN_CHECK_STARS or not np.isfinite(check_rms):
        return "TIED-UNVERIFIED"
    if check_rms <= ACCURACY_STRETCH_MAG:
        return "TIED-STRETCH"
    if check_rms <= ACCURACY_GOAL_MAG:
        return "TIED-GOAL"
    return "TIED-ABOVE-GOAL"


def goal_verdict(verdicts: Sequence[str], n_series_total: int) -> tuple[str, str]:
    """Grade the strategy's calibration goal from the per-block verdicts.

    The goal is "every published magnitude on a standard system, good to
    0.01-0.02 mag".  It has two halves and both must hold:

    * COVERAGE -- every block that carries science must be tied;
    * ACCURACY -- the independent check must meet the stated band.

    Returns ``(verdict, deciding_number)``.  SUPPORTED requires every block
    tied AND the median achieved accuracy inside the goal.  Anything less
    is stated as what it is rather than rounded up.
    """
    tied = [v for v in verdicts if v.startswith("TIED")]
    met = [v for v in verdicts if v in ("TIED-GOAL", "TIED-STRETCH")]
    n_t, n_m = len(tied), len(met)
    deciding = (f"{n_t} of {n_series_total} solved series carry a catalogue "
                f"tie; {n_m} of {n_series_total} meet the "
                f"{ACCURACY_STRETCH_MAG:g}-{ACCURACY_GOAL_MAG:g} mag goal on "
                f"independent check stars")
    if n_t == n_series_total and n_m == n_series_total:
        return "SUPPORTED", deciding
    if n_t == n_series_total:
        return "SUPPORTED-WITH-CAVEATS", deciding
    if n_t == 0:
        return "NOT SUPPORTED", deciding
    return "PARTIALLY SUPPORTED", deciding
