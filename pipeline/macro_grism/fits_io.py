"""HDU-resolution FITS reader for the grism track.

The archive packages 2-D frames three ways (S0b established all three):

* ``plain``       — image in the primary HDU (the recovered master
                    calibrations: ``master_dark_240s_...fts`` is a float32
                    PrimaryHDU).
* ``fpack``       — fpack tile compression: empty primary + a
                    ``CompImageHDU`` (every ``rawimage`` .fts.fz frame).
* ``repackaged``  — era-C reprocessing: empty primary + an UNCOMPRESSED
                    ``ImageHDU`` extension.

Anything else is an unknown layout and must HARD-FAIL with a message that
names the file and describes every HDU — a silent guess about which HDU is
"the image" is how wrong pixels enter a paper.

``classify_hdus`` is pure (it sees only a summary list, so the tests drive
it without any FITS file); ``load_frame`` is the thin impure wrapper.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence


class GrismLayoutError(Exception):
    """Raised when a FITS file's HDU layout matches no known packaging."""


@dataclass(frozen=True)
class HduSummary:
    """What ``classify_hdus`` is allowed to know about one HDU: its class
    name, whether it holds 2-D image data, and its index in the file."""
    index: int
    kind: str          # astropy class name: PrimaryHDU / CompImageHDU / ...
    has_2d_data: bool


def summarize_hdus(hdulist) -> list[HduSummary]:
    """Build the pure-side summary from an open astropy HDUList.

    ``hdu.shape`` is used instead of touching ``.data`` so a compressed
    file is not decompressed twice; shape comes from the header alone.
    """
    out = []
    for i, hdu in enumerate(hdulist):
        shape = getattr(hdu, "shape", ()) or ()
        out.append(HduSummary(index=i, kind=type(hdu).__name__,
                              has_2d_data=len(shape) == 2))
    return out


def classify_hdus(summaries: Sequence[HduSummary]) -> tuple[str, int]:
    """Resolve the packaging: return (layout name, index of the data HDU).

    Rules, in order (first match wins — the orderings encode which layout
    each rule is *allowed* to claim):

    1. primary HDU itself holds a 2-D image        -> ``plain``
    2. exactly one CompImageHDU holds a 2-D image  -> ``fpack``
    3. empty primary + exactly one 2-D ImageHDU    -> ``repackaged``

    Two data-bearing HDUs, or none, is ambiguous — hard fail.  (A file with
    both a primary image AND extensions still resolves as ``plain``: the
    primary image is the canonical pixel source in every archive layout.)
    """
    if not summaries:
        raise GrismLayoutError("empty HDU list")
    if summaries[0].kind == "PrimaryHDU" and summaries[0].has_2d_data:
        return ("plain", 0)
    # Primary carries no image: the data must live in exactly one extension.
    comp = [s for s in summaries if s.kind == "CompImageHDU"
            and s.has_2d_data]
    img = [s for s in summaries if s.kind == "ImageHDU" and s.has_2d_data]
    if len(comp) == 1 and not img:
        return ("fpack", comp[0].index)
    if len(img) == 1 and not comp:
        return ("repackaged", img[0].index)
    # Anything else: refuse loudly, describing what was actually found.
    desc = ", ".join(f"[{s.index}] {s.kind}"
                     f"{' (2-D data)' if s.has_2d_data else ''}"
                     for s in summaries)
    raise GrismLayoutError(f"unknown FITS packaging: {desc}")


def load_frame(path: str, dtype: str = "float32"):
    """Open one archive frame, resolve its packaging by inspection, and
    return ``(data, header, layout)``.

    * ``data`` is a float array (the uint16 raws are promoted so that
      downstream arithmetic — background subtraction, dark subtraction —
      can go negative without wrapping around).
    * ``header`` is the header of the DATA HDU (the raws keep all their
      cards on the CompImageHDU, not the stub primary).
    * ``layout`` names which packaging rule fired, for the DB record.

    Errors carry the path: a batch log line must identify its file.
    """
    import numpy as np
    from astropy.io import fits

    with fits.open(path) as hdulist:
        try:
            layout, idx = classify_hdus(summarize_hdus(hdulist))
        except GrismLayoutError as exc:
            raise GrismLayoutError(f"{path}: {exc}") from exc
        data = np.asarray(hdulist[idx].data, dtype=dtype)
        header = hdulist[idx].header.copy()
    return data, header, layout
