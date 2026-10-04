"""S2 bad-pixel masks per camera: logic, product format and the lookup API.

WHY THIS MODULE EXISTS
----------------------
Three committee findings (DE cross-cutting 5, DE.F1, OA.E8) come down to
one missing product.  The IMX455 cameras deliver 2x2 ON-CAMERA AVERAGES,
so a native pixel that is saturated in every exposure (a "hot" pixel)
never reads 65,535: averaged with three ordinary neighbours it reads

    (65,535 + 3 x level) / 4   ~ 16.6 kADU      (one saturated native pixel)
    (2 x 65,535 + 2 x level)/4 ~ 32.9 kADU      (two)
    (3 x 65,535 + level) / 4   ~ 49.2 kADU      (three)

— the "rails".  S2 v1.2 saw the first of these piling up in the frame
maxima and the grism extractor read it as a 12-bit ADC ceiling
(``SATURATION_ADU = 16300``), masking valid spectrum above it.  There is
no such ceiling; there is a population of bad pixels, and the right
treatment is a MASK applied before detection, trace fitting and frame-max
statistics.  The GSENSE4040 has its own population (hot pixels and
random-telegraph-signal pixels whose level jumps between frames).

HOW A MASK IS BUILT (one method for every camera)
-------------------------------------------------
A bad pixel is one that misbehaves in the same place in frame after frame
whatever the telescope points at.  So the mask is built by PERSISTENCE
over a sample of science frames from many nights and fields:

1. each frame's pixel is compared with the median of its four nearest
   neighbours (:func:`neighbour_excess`) — a local high-pass that removes
   sky, vignetting and nebulosity and keeps single-pixel spikes;
2. a pixel "fires" in a frame when that excess is above
   :data:`HOT_SIGMA` robust sigmas of the frame's own excess distribution
   (and above :data:`HOT_MIN_EXCESS_ADU`), or when its value sits on one
   of the rails computed from ITS neighbours' level
   (:func:`rail_index`);
3. a pixel is masked when it fires in at least :data:`PERSIST_FRACTION`
   of the sample (:class:`PersistenceAccumulator`).

A star fires step 2 as well — but at a different pixel in every pointing,
so it never reaches step 3.  The method needs no darks, which matters: the
camera now on the telescope (QHY600) has none.  Where darks DO exist the
same accumulator is run on them as an independent check, and the campaign
records how well the two masks agree.

Noisy (RTS) pixels are found separately from a stack of same-exposure
darks (:func:`temporal_noise_mask`): their mean can look normal while
their frame-to-frame scatter is many times the read noise.

THE API OTHER PACKAGES CALL
---------------------------
:func:`load_mask` returns a boolean array (True = bad) for a camera and
frame shape; :func:`mask_for_frame` does the same from header values.
Both read the npz products written by the campaign under
``products/detector/badpix/`` and return None — never an invented empty
mask — when no product exists for that camera.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .camera import camera_of

# --------------------------------------------------------------------------
# Bit flags of the stored uint8 mask.
# --------------------------------------------------------------------------
#: Persistently hot: fires the neighbour-excess test in most frames.
BAD_HOT = 1
#: On a rail: at least one saturated native pixel inside an averaged
#: superpixel (IMX455 2x2 only).
BAD_RAIL = 2
#: Temporally noisy (RTS): frame-to-frame scatter far above read noise.
BAD_NOISY = 4

# --------------------------------------------------------------------------
# Tunable constants (single source of truth — the report interpolates these).
# --------------------------------------------------------------------------

#: A pixel fires when its excess over its neighbours' median exceeds this
#: many robust sigmas of the frame's excess distribution ...
HOT_SIGMA = 8.0

#: ... and this many ADU (a floor: on a very quiet frame 8 sigma can be a
#: few ADU, which is noise, not a defect).
HOT_MIN_EXCESS_ADU = 20.0

#: A pixel is masked when it fires in at least this fraction of the
#: sample.  0.6 keeps every pixel that is bad more often than not, and is
#: far above what any sky source can reach across different pointings.
PERSIST_FRACTION = 0.6

#: Rail tolerance: a value within this many ADU of a rail level counts as
#: on it.  The rail level is computed from the pixel's own neighbours, so
#: the tolerance only has to cover their noise and the saturated native
#: pixel's own clip variation.
RAIL_TOL_ADU = 120.0

#: ADC full scale of the IMX455 cameras' native pixels.
NATIVE_FULL_SCALE = 65535.0

#: Temporal-noise mask: a pixel is noisy when its scatter over a dark
#: stack exceeds this multiple of the stack's median per-pixel scatter.
NOISY_FACTOR = 5.0

#: ... and the stack must hold at least this many frames.
NOISY_MIN_FRAMES = 8

#: MAD -> sigma for a Gaussian.
MAD_TO_SIGMA = 1.4826

#: Where the campaign writes mask products (relative to the repo root).
BADPIX_SUBDIR = Path("products") / "detector" / "badpix"


def neighbour_excess(img: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Each pixel minus the median of its four nearest neighbours.

    Returns ``(excess, neighbour_median)`` as float32 arrays of the frame's
    shape.  Edges are handled by reflecting the frame (an edge pixel is
    compared with the neighbours it has).  The median of four values is
    the mean of the middle two — robust to ONE deviant neighbour, which is
    what lets two adjacent bad pixels still be found.
    """
    a = np.asarray(img, dtype=np.float32)
    p = np.pad(a, 1, mode="reflect")
    nb = np.stack([p[:-2, 1:-1], p[2:, 1:-1], p[1:-1, :-2], p[1:-1, 2:]])
    nb.sort(axis=0)
    med = 0.5 * (nb[1] + nb[2])
    return a - med, med


def robust_sigma(values: np.ndarray, sample: int = 2_000_000) -> float:
    """MAD-sigma of an array, from an evenly strided sample (speed)."""
    v = np.asarray(values).ravel()
    if v.size > sample:
        v = v[:: v.size // sample]
    med = np.median(v)
    return float(MAD_TO_SIGMA * np.median(np.abs(v - med)))


def rail_levels(level: float, n_native: int = 4,
                full_scale: float = NATIVE_FULL_SCALE) -> list[float]:
    """The rail values for k = 1 .. n_native-1 saturated native pixels.

    ``level`` is what an ordinary native pixel reads there (pedestal plus
    sky).  For level 303 and a 2x2 average: 16,611 / 32,919 / 49,227 ADU.
    """
    return [(k * full_scale + (n_native - k) * level) / n_native
            for k in range(1, n_native)]


def rail_index(img: np.ndarray, neighbour_median: np.ndarray,
               n_native: int = 4, tol: float = RAIL_TOL_ADU,
               full_scale: float = NATIVE_FULL_SCALE) -> np.ndarray:
    """Per pixel: how many saturated native pixels its value implies.

    0 where the value is not on a rail; k (1 .. n_native-1) where it lies
    within ``tol`` of ``(k * full + (n - k) * local_level) / n`` with the
    local level taken from the pixel's own neighbours.  A value at the
    full scale itself (all native pixels saturated) is NOT a rail — that
    is ordinary saturation, which a bright star produces legitimately.
    """
    a = np.asarray(img, dtype=np.float32)
    out = np.zeros(a.shape, dtype=np.uint8)
    for k in range(1, n_native):
        expect = (k * full_scale + (n_native - k) * neighbour_median) / n_native
        out[np.abs(a - expect) <= tol] = k
    return out


class PersistenceAccumulator:
    """Count, pixel by pixel, in how many frames a pixel fired.

    Feed frames one at a time with :meth:`add` (memory stays at a few
    frame-sized arrays however many frames go through); read the mask
    with :meth:`mask`.  ``n_native`` > 1 switches the rail test on.
    """

    def __init__(self, shape: tuple[int, int], n_native: int = 1):
        self.shape = tuple(shape)
        self.n_native = int(n_native)
        self.n_frames = 0
        self.hot = np.zeros(self.shape, dtype=np.uint16)
        self.rail = np.zeros(self.shape, dtype=np.uint16)
        #: Sum over frames of the rail index where on a rail (for the
        #: per-pixel typical k).
        self.rail_k_sum = np.zeros(self.shape, dtype=np.uint32)

    def add(self, img: np.ndarray) -> dict:
        """Accumulate one frame; returns that frame's summary numbers."""
        if tuple(img.shape) != self.shape:
            raise ValueError(f"frame shape {img.shape} != {self.shape}")
        excess, med = neighbour_excess(img)
        sigma = robust_sigma(excess)
        thr = max(HOT_SIGMA * sigma, HOT_MIN_EXCESS_ADU)
        fired = excess > thr
        self.hot += fired
        n_rail = 0
        if self.n_native > 1:
            k = rail_index(img, med, self.n_native)
            on = k > 0
            self.rail += on
            self.rail_k_sum += k
            n_rail = int(on.sum())
        self.n_frames += 1
        return {"sigma": sigma, "threshold": float(thr),
                "n_fired": int(fired.sum()), "n_rail": n_rail}

    def mask(self, persist_fraction: float = PERSIST_FRACTION) -> np.ndarray:
        """The uint8 bit mask (BAD_HOT | BAD_RAIL) from the counts so far."""
        need = max(2, int(np.ceil(persist_fraction * self.n_frames)))
        out = np.zeros(self.shape, dtype=np.uint8)
        out[self.hot >= need] |= BAD_HOT
        if self.n_native > 1:
            out[self.rail >= need] |= BAD_RAIL
        return out

    def rail_k(self) -> np.ndarray:
        """Typical rail index (rounded mean k over the frames it railed)."""
        with np.errstate(invalid="ignore", divide="ignore"):
            k = np.where(self.rail > 0, self.rail_k_sum / self.rail, 0.0)
        return np.rint(k).astype(np.uint8)


def temporal_noise_mask(stack: np.ndarray,
                        factor: float = NOISY_FACTOR) -> tuple[np.ndarray,
                                                               float]:
    """Noisy (RTS) pixels from a stack of same-exposure darks.

    ``stack`` is (n_frames, ny, nx).  A pixel's scatter is its MAD-sigma
    over the stack (robust to a single cosmic-ray hit); it is flagged when
    that exceeds ``factor`` times the median per-pixel scatter of the whole
    stack.  Returns ``(boolean mask, typical per-pixel sigma)``.  Fewer
    than :data:`NOISY_MIN_FRAMES` frames raises — a scatter from three
    frames is not a measurement.
    """
    s = np.asarray(stack, dtype=np.float32)
    if s.shape[0] < NOISY_MIN_FRAMES:
        raise ValueError(f"need >= {NOISY_MIN_FRAMES} frames, got {s.shape[0]}")
    med = np.median(s, axis=0)
    sig = MAD_TO_SIGMA * np.median(np.abs(s - med[None]), axis=0)
    typical = float(np.median(sig))
    return sig > factor * max(typical, 1e-6), typical


def mask_agreement(mask_a: np.ndarray, mask_b: np.ndarray) -> dict:
    """How well two boolean masks of one detector agree.

    Returns the two counts, the overlap, and the two conditional
    fractions (``a_in_b`` = fraction of a's pixels also in b, and the
    converse).  Used to grade the science-frame persistence mask against
    the dark-frame one where darks exist.
    """
    a = np.asarray(mask_a, dtype=bool)
    b = np.asarray(mask_b, dtype=bool)
    na, nb, both = int(a.sum()), int(b.sum()), int((a & b).sum())
    return {"n_a": na, "n_b": nb, "n_both": both,
            "a_in_b": (both / na) if na else None,
            "b_in_a": (both / nb) if nb else None}


# --------------------------------------------------------------------------
# Product format + lookup API
# --------------------------------------------------------------------------
def mask_key(camera: str, naxis1: int, naxis2: int,
             temp_group: Optional[int] = None,
             orientation: Optional[str] = None) -> str:
    """File stem of a mask product.

    ``ASI_4788x3194_T-10_FM``, ``QHY600_4800x3211_T-20_N``.  Temperature
    group: hot-pixel populations grow with temperature (OA.E8).
    Orientation (``FM`` = FLIPSTAT 'Flip/Mirror', ``N`` = none): a mask is
    a map in FILE coordinates, and the acquisition software's flip setting
    changed during the ASI era — the same sensor pixel is stored at a
    different place before and after.
    """
    key = f"{camera}_{int(naxis1)}x{int(naxis2)}"
    if temp_group is not None:
        key += f"_T{int(temp_group)}"
    if orientation is not None:
        key += f"_{orientation}"
    return key


def parse_mask_key(key: str) -> dict:
    """Inverse of :func:`mask_key` (missing parts come back as None)."""
    parts = key.split("_")
    nx, ny = (int(v) for v in parts[1].split("x"))
    tg, orient = None, None
    for part in parts[2:]:
        if part.startswith("T") and part[1:].lstrip("-").isdigit():
            tg = int(part[1:])
        else:
            orient = part
    return {"camera": parts[0], "naxis1": nx, "naxis2": ny,
            "temp_group": tg, "orientation": orient}


def orientation_of(flipstat: Optional[str]) -> str:
    """Orientation tag of a FLIPSTAT header value (``FM``, ``N``, ...)."""
    v = str(flipstat or "").strip()
    if not v:
        return "N"
    if v.lower().replace(" ", "") == "flip/mirror":
        return "FM"
    return "".join(ch for ch in v if ch.isalnum())[:8] or "N"


def save_mask(directory: Path, key: str, mask: np.ndarray,
              hot_count: np.ndarray, rail_k: Optional[np.ndarray],
              meta: dict) -> Path:
    """Write one mask product atomically; returns its path.

    The npz holds the uint8 bit mask, the per-pixel hot count (so a user
    can re-threshold without re-reading the archive), the per-pixel rail
    index, and a JSON string of metadata (frames used, thresholds, counts).
    """
    import os
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{key}.npz"
    tmp = path.with_suffix(f".tmp{os.getpid()}.npz")
    np.savez_compressed(
        tmp, mask=mask.astype(np.uint8),
        hot_count=hot_count.astype(np.uint16),
        rail_k=(rail_k if rail_k is not None
                else np.zeros((0, 0), np.uint8)),
        meta=np.array(json.dumps(meta, sort_keys=True)))
    tmp.replace(path)
    return path


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent.parent


def available_masks(directory: Optional[Path] = None) -> list[str]:
    """Keys of every mask product on disk."""
    d = Path(directory) if directory else _repo_root() / BADPIX_SUBDIR
    return sorted(p.stem for p in d.glob("*.npz")) if d.exists() else []


def load_mask(camera: str, shape: tuple[int, int],
              temp_c: Optional[float] = None,
              flipstat: Optional[str] = None,
              flags: int = BAD_HOT | BAD_RAIL | BAD_NOISY,
              crop: Optional[tuple[int, int]] = None,
              directory: Optional[Path] = None) -> Optional[np.ndarray]:
    """Boolean bad-pixel mask (True = bad) for a camera and frame shape.

    Parameters
    ----------
    camera
        ``"AC4040"``, ``"iKon"``, ``"ASI"`` or ``"QHY600"`` (see
        :func:`rlmt_diagnostics.camera.camera_of`).
    shape
        (ny, nx) of the frame the mask will be applied to.
    temp_c
        Sensor temperature; among the masks of the right orientation the
        nearest temperature group is used.  None takes the coldest.
    flipstat
        The frame's FLIPSTAT header value (None or '' = no flip).  Only a
        mask built from frames of the SAME orientation is returned.
    flags
        Which defect classes to include (bitwise OR of the BAD_* flags).
    crop
        ``(dy, dx)`` offset of a CROPPED frame inside the raw frame (the
        2026 reduction drops a few rows and columns —
        ``s2_recon_eras.crop_dy/crop_dx``).  The mask is built on raw
        geometry; with ``crop`` it is cut to ``shape`` at that offset.

    Returns None when no product matches — the caller must then proceed
    WITHOUT a mask and say so; this function never returns an all-False
    array for a camera, geometry or orientation it knows nothing about.
    """
    d = Path(directory) if directory else _repo_root() / BADPIX_SUBDIR
    ny, nx = int(shape[0]), int(shape[1])
    want = orientation_of(flipstat)
    cands = []
    for p in (d.glob(f"{camera}_*.npz") if d.exists() else []):
        k = parse_mask_key(p.stem)
        if k["orientation"] is not None and k["orientation"] != want:
            continue
        fits_raw = (k["naxis1"], k["naxis2"]) == (nx, ny)
        fits_crop = (crop is not None and k["naxis1"] >= nx
                     and k["naxis2"] >= ny)
        if fits_raw or fits_crop:
            cands.append((p, k["temp_group"], fits_raw))
    if not cands:
        return None
    if temp_c is None:
        cands.sort(key=lambda c: (not c[2], c[1] if c[1] is not None else 0))
    else:
        cands.sort(key=lambda c: (not c[2], abs((c[1] if c[1] is not None
                                                  else temp_c) - temp_c)))
    path, _tg, fits_raw = cands[0]
    with np.load(path, allow_pickle=False) as npz:
        m = (npz["mask"] & np.uint8(flags)) > 0
    if not fits_raw:
        dy, dx = crop
        m = m[dy:dy + ny, dx:dx + nx]
        if m.shape != (ny, nx):
            return None
    return m


def mask_for_frame(instrume: Optional[str], naxis1: int, naxis2: int,
                   readoutm: Optional[str] = None,
                   ccd_temp: Optional[float] = None,
                   flipstat: Optional[str] = None,
                   flags: int = BAD_HOT | BAD_RAIL | BAD_NOISY,
                   crop: Optional[tuple[int, int]] = None,
                   directory: Optional[Path] = None
                   ) -> Optional[np.ndarray]:
    """:func:`load_mask` addressed by a frame's header values."""
    cam = camera_of(instrume, naxis1, readoutm)
    return load_mask(cam, (int(naxis2), int(naxis1)), temp_c=ccd_temp,
                     flipstat=flipstat, flags=flags, crop=crop,
                     directory=directory)


def mask_metadata(key: str, directory: Optional[Path] = None) -> dict:
    """The metadata dict stored with a mask product."""
    d = Path(directory) if directory else _repo_root() / BADPIX_SUBDIR
    with np.load(d / f"{key}.npz", allow_pickle=False) as npz:
        return json.loads(str(npz["meta"]))


def sample_evenly(items: Sequence, n: int) -> list:
    """At most ``n`` items spread evenly through a sequence."""
    if len(items) <= n:
        return list(items)
    idx = np.unique(np.linspace(0, len(items) - 1, n).round().astype(int))
    return [items[i] for i in idx]
