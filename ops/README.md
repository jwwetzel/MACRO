# `ops/` — requests to the observatory

| File | Status |
|---|---|
| `2026-10_observatory_request_rev3.md` | **Current.** Revision 3 — DRAFT until James sends it. Written for the QHY600 under pyscope and for the sky actually available after the October 2026 re-opening. |
| `2026-08_observatory_request.md` | **SUPERSEDED by revision 3 (2026-10-03). Do not act on it.** Kept unedited as the record of what was asked in August; the plan review of 2026-10-03 (finding U10) found that it requests modes of cameras no longer mounted and observations the sky does not allow. |
| `generated/` | Script output embedded in revision 3: visibility tables and figure (`pipeline/scripts/ops_visibility.py`), exposure derivation, calibration census and figure (`pipeline/scripts/ops_exposure.py`). Never edited by hand. |
| `eruption_block/` | DRAFT pyscope-style schedule blocks for a T CrB eruption, with operating notes. Not loaded on the telescope. |

The superseded file is not edited in place because it is a declared output of the
pipeline's provenance graph (stage `OPS` in `pipeline/macro_core/provenance.py`); changing
its bytes would change its digest. Re-pointing that stage at revision 3 is a ledger change
for the chair (see `committee/work/ops/REPORT.md`).

## Regenerating revision 3's numbers

```
PY=/opt/miniconda3/envs/rlmt-checks/bin/python
$PY pipeline/scripts/ops_visibility.py --inject ops/2026-10_observatory_request_rev3.md
$PY pipeline/scripts/ops_exposure.py measure      # reads ~200 archive frames, read-only
$PY pipeline/scripts/ops_exposure.py report --inject ops/2026-10_observatory_request_rev3.md
$PY pipeline/scripts/ops_exposure.py report --inject ops/eruption_block/README.md
$PY pipeline/scripts/ops_visibility.py --airmass 3 --out-dir ops/eruption_block/generated \
    --inject ops/eruption_block/README.md
```
