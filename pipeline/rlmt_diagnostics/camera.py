"""Camera / configuration identity for detector measurements (pure logic).

WHY THIS MODULE EXISTS
----------------------
S2 v1.2 keyed everything on the READOUTM header string and called the
result a "readout mode".  The 2026-10-03 review (DE, opening table) showed
that the archive holds FOUR physical cameras, that READOUTM is a driver
label and not an instrument, and that two of those labels (``Mode0`` and
``Fast``) belong to different sensors while one label (``High Gain``)
spans two EGAIN epochs with measurably different clip levels.  A gain or a
read noise is a property of a (camera, readout configuration) pair, so the
flat-pair photon-transfer stage needs a key that says which pair a frame
belongs to.  This module is that key, and nothing else:

* :func:`camera_of` — which sensor, from INSTRUME (and geometry as the
  fallback for the 2026 pyscope frames whose INSTRUME card is blank);
* :func:`config_key` — the (camera, readout, binning, gain epoch) label
  every flat-PTC table is keyed on;
* :func:`mode_label_of` — the S2 v1.2 mode label a configuration belongs
  to, so per-configuration results can still be joined to the existing
  per-mode tables (``s2_ceiling_modes``, ``detector_params``).

The EGAIN card is used ONLY to tell epochs apart (1.054 vs 1.057 on the
AC4040) — never as a gain.  What the gain IS is the measurement this key
exists to organise.

Pure functions on header values; no I/O.
"""

from __future__ import annotations

from typing import Optional

from .ceiling import mode_group

#: Sensor identity by INSTRUME substring (lower-cased match).  The strings
#: are what MaxIm DL / pyscope wrote; the camera names are the committee's
#: (detector-engineer memo, opening table).
_INSTRUME_CAMERA = (
    ("dl imaging", "AC4040"),          # SBIG AC4040M, GSENSE4040 sCMOS
    ("andor", "iKon"),                 # Andor iKon CCD, 2048^2
    ("asi camera", "ASI"),             # ZWO ASI (IMX455)
    ("qhy", "QHY600"),                 # QHY600M (IMX455)
)

#: Native (unbinned) sensor widths.  Used only when INSTRUME is blank —
#: the 2026 pyscope frames — to recognise the camera from its geometry:
#: 2x2-binned raw QHY frames are 4800 wide (reduced: 4787), ASI 4788.
_QHY_BINNED_WIDTHS = (4800, 4787)

#: Native pixels combined into one stored pixel, by camera and XBINNING.
#: Only the IMX455 cameras bin on-camera in this archive.
CAMERA_NATIVE_PIX_UM = {"AC4040": 9.0, "iKon": 13.5, "ASI": 3.76,
                        "QHY600": 3.76}


def camera_of(instrume: Optional[str], naxis1: Optional[float] = None,
              readoutm: Optional[str] = None) -> str:
    """The physical camera a frame came from.

    INSTRUME decides when present.  When it is blank (598 pyscope frames
    of 2026-06/07) the stored width decides: 4800/4787 columns is the
    2x2-binned QHY600 and nothing else in the archive.  Anything
    unrecognised returns ``"unknown"`` — an honest label that keeps such
    frames out of every per-camera fit instead of guessing them into one.
    """
    s = (instrume or "").strip().lower()
    for needle, cam in _INSTRUME_CAMERA:
        if needle in s:
            return cam
    if naxis1 is not None and int(naxis1) in _QHY_BINNED_WIDTHS:
        return "QHY600"
    r = (readoutm or "").strip()
    if r.startswith(("High Gain", "Low Gain")):
        return "AC4040"
    if "MHz" in r:
        return "iKon"
    if r == "Mode0":
        return "ASI"
    if r == "Fast":
        return "QHY600"
    return "unknown"


def config_key(instrume: Optional[str], readoutm: Optional[str],
               xbinning: Optional[float], naxis1: Optional[float] = None,
               egain: Optional[float] = None,
               preamp: Optional[str] = None) -> str:
    """The camera/configuration key of a frame.

    Examples (the configurations the archive actually holds)::

        AC4040 High Gain e1.054      AC4040 High Gain e1.057
        AC4040 StackPro e1.054       AC4040 StackPro e1.057
        AC4040 Low Gain
        iKon 1MHz   iKon 3MHz   iKon 5MHz   (iKon 1MHz 4x with a preamp)
        ASI Mode0 2x2      ASI Mode0 1x1      ASI Mode0 2x2 e0.780
        QHY600 Fast 2x2    QHY600 pyscope 2x2

    Rules: the AC4040's two EGAIN epochs are separate configurations (the
    clip level differs between them — DE.F9); StackPro is separate from
    plain High Gain (a 16-read sum has 16x the read variance); the ASI's
    binning is part of the key because 2x2 there is an on-camera AVERAGE
    with a different effective gain; the handful of ASI frames written
    with EGAIN 0.78 (a different GAIN setting) are kept apart; and the
    QHY600's MaxIm-era ``Fast`` frames and pyscope-era blank-READOUTM
    frames carry separate keys until a measurement shows them identical.
    """
    cam = camera_of(instrume, naxis1, readoutm)
    r = (readoutm or "").strip()
    xb = int(xbinning) if xbinning else 1
    if cam == "AC4040":
        if r.startswith("Low Gain"):
            return "AC4040 Low Gain"
        base = "StackPro" if "StackPro" in r else "High Gain"
        # Epoch from the EGAIN card: 1.0537 -> 1.054, 1.0570 -> 1.057.
        tag = f" e{egain:.3f}" if egain and egain > 0 else ""
        return f"AC4040 {base}{tag}"
    if cam == "iKon":
        speed = r.split("MHz")[0].strip() + "MHz" if "MHz" in r else "?"
        key = f"iKon {speed}" + (f" {xb}x{xb}" if xb != 1 else "")
        # The Andor pre-amplifier setting (header GAIN: '1x', '2x', '4x')
        # changes the e-/ADU by its factor: the 2024-04-09/10 engineering
        # flats were taken at 1x (K ~ 3.5), everything after at 4x
        # (K ~ 0.98).  The manifest does not carry this card yet, so the
        # key only gains the suffix when the caller has read the header.
        p_ = str(preamp or "").strip()
        return key + (f" {p_}" if p_ else "")
    if cam == "ASI":
        key = f"ASI Mode0 {xb}x{xb}"
        # The archive's standard ASI setting writes EGAIN 0.2467; a few
        # 2025 frames carry 0.78 (another GAIN setting) and must not be
        # pooled with it.
        if egain and egain > 0 and abs(egain - 0.2467) > 0.02:
            key += f" e{egain:.3f}"
        return key
    if cam == "QHY600":
        era = "Fast" if r == "Fast" else "pyscope"
        return f"QHY600 {era} {xb}x{xb}"
    return f"unknown {r or '(blank)'}"


def mode_label_of(config: str) -> Optional[str]:
    """The S2 mode label (``s2_ceiling_modes.mode``) a configuration maps to.

    The inverse of the pooling S2 v1.2 did: both AC4040 High Gain epochs
    belong to mode ``High Gain``; the ASI's binned and unbinned frames are
    both ``Mode0``; the pyscope QHY frames are ``(blank 2026)``.  None for
    an unknown configuration.
    """
    if config.startswith("AC4040 High Gain"):
        return "High Gain"
    if config.startswith("AC4040 StackPro"):
        return "High Gain StackPro"
    if config.startswith("AC4040 Low Gain"):
        return "Low Gain"
    if config.startswith("iKon "):
        speed = config.split()[1]
        return f"{speed} High Sensitivity 16-bit"
    if config.startswith("ASI Mode0"):
        return "Mode0"
    if config.startswith("QHY600 Fast"):
        return "Fast"
    if config.startswith("QHY600 pyscope"):
        return mode_group(None)
    return None


def n_native_per_pixel(config: str) -> int:
    """Native sensor pixels combined into one stored pixel (1 or 4)."""
    return 4 if config.rstrip().endswith("2x2") or " 2x2 " in config else 1


def camera_name(config: str) -> str:
    """The camera part of a configuration key (``"AC4040"``, ``"ASI"`` ...)."""
    return config.split()[0]
