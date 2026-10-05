"""Named frame samples of the grism track — one definition each.

Every G runner that needs "the T CrB series" or "the identity-gate
controls" gets it from here, so two scripts can never disagree about
which frames a number was computed from.  All functions take a READ-ONLY
manifest connection (the snapshot or the shared manifest) and return
task dicts ready for :func:`macro_grism.store.reduce_batch`.

Samples
-------
``tcrb``              every canonical raw T CrB frame taken through a
                      NAMED grism (hrg / lrg): the 2025 series.
``tcrb_slot``         T CrB frames through slot '6' or 'W' (2023-24),
                      whose identity as spectra is S2c's measurement.
``identity_controls`` non-T CrB frames through the same grisms in the
                      same mechanical epoch and exposure regime — the
                      population the pixel identity gate must REJECT.
``ngc5548_slot``      NGC 5548 slot-'6' frames (F-6 re-issue).
"""

from __future__ import annotations

import numpy as np

from . import config as gconfig

_COLS = ("obs_rowid, path, filter, night, jd, exptime, target_best, "
         "pointing_offset_deg")


def _rows(mcon, where: str, params=(), raw_only: bool = True) -> list[dict]:
    """Canonical frames matching ``where``.  ``raw_only`` restricts to
    the ``rawimage`` tree; it is switched off for campaigns whose
    canonical copies live in a partner tree (NGC 5548: ``macalester``)."""
    cols = [c.strip() for c in _COLS.split(",")]
    tree = "AND tree = 'rawimage' " if raw_only else ""
    cur = mcon.execute(
        f"SELECT {_COLS} FROM frames WHERE is_canonical = 1 "
        f"{tree}AND {where} ORDER BY night, path", params)
    return [dict(zip(cols, r)) for r in cur]


def tcrb(mcon) -> list[dict]:
    """The T CrB grism series (hrg + lrg), every frame, with the 2-D
    diagnostics switched on (variance, null apertures, legacy sky)."""
    out = _rows(mcon, "lower(target_best) = 't crb' "
                      "AND lower(filter) IN ('hrg','lrg')")
    for t in out:
        t["sample"], t["diagnostics"] = "tcrb", True
    return out


def tcrb_slot(mcon) -> list[dict]:
    """T CrB through slots '6' and 'W' (AC4040, 2023-05 -> 2024-03)."""
    out = _rows(mcon, "lower(target_best) = 't crb' "
                      "AND lower(filter) IN ('6','w')")
    for t in out:
        t["sample"], t["diagnostics"] = "tcrb_slot", True
    return out


def ngc5548_slot(mcon) -> list[dict]:
    """NGC 5548 through slot '6' (AC4040, 2023)."""
    out = _rows(mcon, "lower(target_best) IN ('ngc 5548','ngc5548') "
                      "AND filter = '6'", raw_only=False)
    for t in out:
        t["sample"] = "ngc5548_slot"
    return out


#: Identity-gate controls: exposure window (s) and sample size.  The
#: controls should look as much like a T CrB frame as a non-T CrB frame
#: can — same grisms, same camera state, long exposures of faint-ish
#: targets — because the gate's job is to tell T CrB from THOSE, not
#: from a 0.1 s frame of Spica.
CONTROL_EXPTIME = (60.0, 400.0)
CONTROL_N = 160
CONTROL_SEED = 20261003


def identity_controls(mcon) -> list[dict]:
    """Non-T CrB named-grism frames from the T CrB series' mechanical
    epoch, drawn at random (fixed seed) with at most 4 frames per
    target per grism so no single star dominates."""
    ep = gconfig.mech_epoch_for("2025-03-20")
    rows = _rows(mcon, "lower(filter) IN ('hrg','lrg') "
                       "AND lower(target_best) NOT LIKE 't crb%' "
                       "AND night BETWEEN ? AND ? "
                       "AND exptime BETWEEN ? AND ?",
                 (ep.first_night, ep.last_night, *CONTROL_EXPTIME))
    rng = np.random.default_rng(CONTROL_SEED)
    order = rng.permutation(len(rows))
    seen, out = {}, []
    for i in order:
        r = rows[i]
        key = ((r["target_best"] or "").lower(), r["filter"].lower())
        if seen.get(key, 0) >= 4:
            continue
        seen[key] = seen.get(key, 0) + 1
        r["sample"] = "identity_control"
        out.append(r)
        if len(out) >= CONTROL_N:
            break
    out.sort(key=lambda r: (r["night"], r["path"]))
    return out


SAMPLES = {"tcrb": tcrb, "tcrb_slot": tcrb_slot,
           "ngc5548_slot": ngc5548_slot,
           "identity_controls": identity_controls}
