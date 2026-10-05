"""macro_grism.config — the ONE place the grism track reads hardware facts.

WHY THIS MODULE EXISTS
----------------------
The committee review of 2026-10-03 found three hardware numbers typed into
``extract.py`` that were each wrong, and wrong in a way nobody could see
from the code that used them:

* the conversion gain (header ``EGAIN`` 0.2467 e-/ADU is the UNBINNED gain;
  the Mode0 frames are 2x2 *average*-binned, so one stored ADU stands for
  four native pixels and the effective gain is ~1.0 e-/ADU — DE.F1);
* the read noise (3.5 e- at the wrong gain made a 201 ADU^2 read term where
  bias pairs measure ~4.4 ADU^2);
* the "16.3 kADU rail" (there is no such rail: the 16.6 kADU pile-up is one
  saturated *native* hot pixel averaged with three normal ones — OA.E8).

The fix is not three better constants.  It is that nothing downstream may
hold a detector constant at all.  Every consumer in ``macro_grism`` asks
:func:`detector_for` and gets a :class:`Detector` record that says where
its numbers came from and whether they are still provisional.  When the
``detector`` work package lands its measured photon-transfer table the
numbers change HERE — either by editing :data:`DETECTOR_DEFAULTS` or,
without touching code, by dropping a JSON override file at
:data:`OVERRIDE_JSON` — and every variance, every saturation flag and every
error bar downstream follows.

The same reasoning applies to the *mechanical* history.  A grism bolted at
a fixed distance from a detector has one dispersion and one sign until
somebody moves the camera or the wheel (TE.F1, TE.F3, PH.P7).  The unit a
wavelength solution belongs to is therefore the **mechanical epoch**, not
the header "era" (era 76 alone spans a 180-degree camera rotation).  The
formal ``mech_epoch`` table is being built by the ``foundation-s0``
package; until it lands, :data:`MECH_EPOCHS` below is this track's working
copy of the same table, derived from the telescope engineer's evidence
(TE.F1) plus the grism-trace position angles in ``frame_dispersion``.  It
is keyed by NIGHT so the two can be reconciled by a join, and the report
lists every boundary so a disagreement is visible.

Nothing here touches pixels or disk except :func:`load_overrides`, which
reads one optional JSON file.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

# --------------------------------------------------------------------------
# Locations
# --------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PRODUCTS = REPO_ROOT / "products" / "grism"

#: The grism track's own database for this wave (the shared manifest is
#: read-only; the integrated rebuild lands these tables in the manifest).
GRISM_DB = PRODUCTS / "grism.sqlite"

#: Optional override file: ``{"ASI-Mode0-bin2avg": {"gain_e_per_adu": 1.03,
#: "read_noise_adu": 2.1, "provenance": "...", "provisional": false}}``.
#: This is the drop-in point for the ``detector`` package's final numbers.
OVERRIDE_JSON = PRODUCTS / "detector_config.json"


# --------------------------------------------------------------------------
# Detectors
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Detector:
    """Everything the extraction needs to know about one camera
    configuration.  All ADU quantities are in STORED (binned) ADU.

    ``gain_e_per_adu`` is the EFFECTIVE gain of a stored pixel — for an
    average-binned frame that is ``n_native x`` the native gain, which is
    what a flat-pair photon-transfer curve measures directly.

    ``native_clip_adu`` is the ADC ceiling of one NATIVE pixel.  For an
    average-binned frame a stored pixel can hide one saturated native
    pixel at a value far below the clip, so saturation is judged with
    :meth:`saturation_cap_adu`, not with the clip itself.
    """
    key: str
    camera: str
    gain_e_per_adu: float
    read_noise_adu: float
    pedestal_adu: float
    native_clip_adu: float
    n_native: int            # native pixels averaged into one stored pixel
    binning: str             # 'average' | 'sum' | 'none'
    pixel_um: float          # STORED pixel pitch (micron), binning included
    provenance: str
    provisional: bool

    def variance_adu2(self, counts_adu):
        """Per-pixel variance model in ADU^2 for raw stored counts:

            V = RN^2 + max(counts - pedestal, 0) / K

        read noise plus photon noise of everything above the electronic
        pedestal (source + sky + dark — their photons are all noise).  The
        pedestal itself carries no shot noise; the superseded model charged
        shot noise on it, which alone inflated the faint-pixel variance by
        ~300 ADU / K.
        """
        import numpy as np
        signal = np.clip(np.asarray(counts_adu, dtype=float)
                         - self.pedestal_adu, 0.0, None)
        return self.read_noise_adu ** 2 + signal / self.gain_e_per_adu

    def saturation_cap_adu(self, peaking: float = None) -> float:
        """Stored-ADU level above which a NATIVE pixel inside the stored
        pixel may have clipped (standing rule 4: saturation is judged in
        native pixels).

        For average binning a stored value S can hide a native peak of up
        to ``peaking x S``; DE.F9 measured 1.07-1.16 for seeing of 6-4
        native pixels.  The cap is therefore ``clip / peaking`` with the
        worst measured peaking, :data:`NATIVE_PEAKING`.  For unbinned or
        summed data the cap is the clip (less one ADU of headroom).
        """
        if self.binning != "average":
            return self.native_clip_adu - 1.0
        p = NATIVE_PEAKING if peaking is None else peaking
        return (self.native_clip_adu
                - self.pedestal_adu) / p + self.pedestal_adu


#: Stored pixel pitch the dispersion SEEDS are quoted for (the IMX455
#: cameras, 3.76 um native, 2x2 binned).  A seed is rescaled by pitch for
#: other cameras: the grism's dispersion is fixed in A per MICRON.
REFERENCE_PIXEL_UM = 7.52

#: Worst-case ratio of the brightest native pixel to the stored 2x2 average
#: around a stellar peak (DE.F9: 1.07-1.16 for FWHM 6-4 native px).  A
#: grism trace is at least as broad across the dispersion as a star, so
#: this bounds the trace too.
NATIVE_PEAKING = 1.16

#: The detector table.  PROVISIONAL rows carry the detector engineer's
#: quick-look flat-pair numbers (DE memo, 2026-10-03) and are replaced by
#: the ``detector`` package's measured table (override JSON or edit).
DETECTOR_DEFAULTS = {
    "ASI-Mode0-bin2avg": Detector(
        key="ASI-Mode0-bin2avg", camera="ZWO ASI (IMX455)",
        gain_e_per_adu=1.03, read_noise_adu=2.1, pedestal_adu=303.0,
        native_clip_adu=65535.0, n_native=4, binning="average", pixel_um=7.52,
        provenance="DE.F1 quick-look flat-pair PTC, Calibrations/"
                   "2024-12-26 (K 1.01-1.06 over nine pairs); bias-pair "
                   "RN 2.1 ADU; pedestal = master-dark median",
        provisional=True),
    "QHY600-bin2avg": Detector(
        key="QHY600-bin2avg", camera="QHY600M (IMX455)",
        gain_e_per_adu=1.0, read_noise_adu=2.1, pedestal_adu=172.0,
        native_clip_adu=65535.0, n_native=4, binning="average", pixel_um=7.52,
        provenance="PLACEHOLDER: no QHY flat pair exists on disk (DE: "
                   "'unknown until October flats'); pedestal from "
                   "s2 recon era 78; G-2 sky-PTC cross-check in "
                   "g_variance_check is the only measurement",
        provisional=True),
    "AC4040-HighGain": Detector(
        key="AC4040-HighGain", camera="SBIG AC4040M (GSENSE4040)",
        gain_e_per_adu=1.06, read_noise_adu=3.93, pedestal_adu=93.0,
        native_clip_adu=3496.0, n_native=1, binning="none", pixel_um=9.0,
        provenance="DE.F2 quick-look V flat pairs (K 1.02-1.06); S2 "
                   "read noise 3.93 ADU, pedestal 93, 12-bit clip 3496",
        provisional=True),
    "AC4040-StackPro": Detector(
        key="AC4040-StackPro", camera="SBIG AC4040M (GSENSE4040)",
        gain_e_per_adu=1.06, read_noise_adu=15.59, pedestal_adu=1484.5,
        native_clip_adu=56062.0, n_native=1, binning="sum", pixel_um=9.0,
        provenance="S2 v1.2: 16-sub-read sum, RN 15.59 ADU, pedestal "
                   "1484.5, ceiling 56062; gain as High Gain",
        provisional=True),
    "iKon-1MHz": Detector(
        key="iKon-1MHz", camera="Andor iKon CCD",
        gain_e_per_adu=0.977, read_noise_adu=7.3, pedestal_adu=305.0,
        native_clip_adu=64674.0, n_native=1, binning="none", pixel_um=13.5,
        provenance="DE.F2 quick-look flat pairs K 0.977 +/- 0.005, RN "
                   "7.3 ADU; S2 ceiling 64674; pedestal s2 recon era 47",
        provisional=True),
}


def load_overrides(path: Path = OVERRIDE_JSON) -> dict:
    """The optional override file, as ``{key: {field: value}}`` (empty
    when the file is absent — absence is the normal state this wave)."""
    if not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text())


def detector_table(overrides: Optional[dict] = None) -> dict:
    """The effective detector table: defaults with overrides applied
    field by field.  An override for an unknown key is an error — a typo
    must not silently leave a provisional number in force."""
    table = dict(DETECTOR_DEFAULTS)
    ov = load_overrides() if overrides is None else overrides
    for key, fields in ov.items():
        if key not in table:
            raise KeyError(f"detector override for unknown key {key!r}")
        table[key] = replace(table[key], **fields)
    return table


def detector_key(instrume: Optional[str], readoutm: Optional[str],
                 night: Optional[str] = None) -> str:
    """Which detector row a frame belongs to, from its INSTRUME and
    READOUTM cards (night is the tie-breaker for blank/pyscope headers).

    Camera identity is the header's ``INSTRUME`` (DE: READOUTM is only a
    driver label); READOUTM then separates the AC4040's two modes.
    """
    ins = (instrume or "").lower()
    rdm = (readoutm or "").lower()
    if "asi" in ins or rdm == "mode0":
        return "ASI-Mode0-bin2avg"
    if "qhy" in ins or rdm == "fast":
        return "QHY600-bin2avg"
    if "andor" in ins or "mhz" in rdm:
        return "iKon-1MHz"
    if "stackpro" in rdm:
        return "AC4040-StackPro"
    if "dl imaging" in ins or "gain" in rdm:
        return "AC4040-HighGain"
    # Blank headers: fall back on the mechanical epoch's camera.
    if night is not None:
        ep = mech_epoch_for(night)
        if ep is not None:
            return ep.default_detector
    raise KeyError(f"cannot identify detector: INSTRUME={instrume!r} "
                   f"READOUTM={readoutm!r} night={night!r}")


def detector_for(instrume: Optional[str], readoutm: Optional[str],
                 night: Optional[str] = None,
                 overrides: Optional[dict] = None) -> Detector:
    """The :class:`Detector` record for a frame's header cards."""
    return detector_table(overrides)[detector_key(instrume, readoutm, night)]


# --------------------------------------------------------------------------
# Mechanical epochs (working copy — see the module docstring)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class MechEpoch:
    """One mechanical state of camera + wheel, bounded by nights
    (inclusive, ``night`` = local observing night as in ``frames.night``).
    """
    epoch_id: str
    camera: str
    first_night: str
    last_night: str
    default_detector: str
    evidence: str


#: Working mechanical-epoch table for the grism track.  Sources: TE.F1's
#: state table (S1 rotation steps, FLIPSTAT, wheel renames) and the grism
#: trace position angles in ``frame_dispersion`` (slot '6' traces flip from
#: PA 177.7 to PA 2.0 across the 2023 monsoon: the element or the camera
#: was turned end for end, so the AC4040 era is two grism epochs).
MECH_EPOCHS = (
    MechEpoch("AC4040-a", "SBIG AC4040M", "2023-02-01", "2023-08-31",
              "AC4040-HighGain",
              "slot-6 trace PA 177.7 deg (2023-02..07); rot 179.0-179.4"),
    MechEpoch("AC4040-b", "SBIG AC4040M", "2023-09-01", "2024-03-31",
              "AC4040-HighGain",
              "slot-6 trace PA 2.0 deg (2023-10..2024-03): dispersion "
              "direction reversed relative to AC4040-a"),
    MechEpoch("ANDOR-1", "Andor iKon", "2024-04-01", "2024-10-31",
              "iKon-1MHz", "rot 180.8 deg; HaGrism/OGGrism wheel names"),
    MechEpoch("ANDOR-2", "Andor iKon", "2024-11-01", "2024-11-07",
              "iKon-1MHz", "re-seated, rot 180.0 deg"),
    MechEpoch("ANDOR-3", "Andor iKon", "2024-11-08", "2024-11-29",
              "iKon-1MHz", "re-seated, rot 180.7; wheel renamed hrg/lrg"),
    MechEpoch("ANDOR-4", "Andor iKon", "2024-11-30", "2024-12-13",
              "iKon-1MHz", "camera rotated, rot 204.3 deg"),
    MechEpoch("ASI-pre", "ZWO ASI (IMX455)", "2024-12-14", "2025-09-30",
              "ASI-Mode0-bin2avg",
              "rot +0.43 deg, FLIPSTAT 'Flip/Mirror', MaxIm 6.40"),
    MechEpoch("ASI-post", "ZWO ASI (IMX455)", "2025-10-01", "2026-03-20",
              "ASI-Mode0-bin2avg",
              "rot 179.75 deg, FLIPSTAT blank, MaxIm 6.30 (camera turned "
              "180 deg over the 2025 monsoon and re-seated)"),
    MechEpoch("QHY-n1", "QHY600M", "2026-03-21", "2026-03-21",
              "QHY600-bin2avg", "commissioning night, rot 3.3 deg"),
    MechEpoch("QHY", "QHY600M", "2026-03-22", "2099-12-31",
              "QHY600-bin2avg",
              "rot 179.74 deg; wheel reloaded, lrg re-seated (trace PA "
              "+7 deg); pyscope-native headers from 2026-06-28"),
)


def mech_epoch_for(night: str) -> Optional[MechEpoch]:
    """The mechanical epoch containing ``night`` ('YYYY-MM-DD'), or None
    for a night outside every known epoch (pre-2023)."""
    for ep in MECH_EPOCHS:
        if ep.first_night <= night <= ep.last_night:
            return ep
    return None


def mech_epoch_id(night: str) -> Optional[str]:
    ep = mech_epoch_for(night)
    return ep.epoch_id if ep else None


# --------------------------------------------------------------------------
# Grism identity
# --------------------------------------------------------------------------
#: FILTER-card spellings of the two units across the archive's three
#: vocabularies.  Slot '6' and 'W' are NOT here: what sat in those slots is
#: a measurement (S2c / F-6), never an assumption.
HIGH_RES_NAMES = ("hrg", "hagrism", "hag")
LOW_RES_NAMES = ("lrg", "oggrism")


def grism_unit(filter_card: Optional[str]) -> Optional[str]:
    """'hrg' / 'lrg' for a FILTER card that names a grism, else None."""
    f = (filter_card or "").strip().lower()
    if f in HIGH_RES_NAMES:
        return "hrg"
    if f in LOW_RES_NAMES:
        return "lrg"
    return None
