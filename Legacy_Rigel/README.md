# Legacy archive (was "Legacy_Rigel") — census only

Committee ruling U8 (2026-10-03): this is a **census (L0–L1) with a pre-registered go/no-go (L2)**,
not a science project. The directory keeps its old name until the chair renames the ledger entry;
the premise behind the name is wrong — only the first months of 2015 are the Rigel/PL16803, the rest
is the 0.5 m with a succession of cameras (see the census page, section 3).

| What | Where |
|---|---|
| Go/no-go criteria, written before any header was read | `notes/GO_NOGO_PREREGISTERED.md` |
| Census database (every number on the page is a query of it) | `products/legacy/legacy.sqlite` |
| Evidence page | `docs/Legacy_Rigel/legacy_census.html` |
| Pure rules (unit-tested) | `pipeline/macro_legacy/census.py`, `scan.py` |
| Builds | `pipeline/scripts/build_legacy_scan.py` → `build_legacy_census.py` → `build_legacy_external.py` → `build_legacy_census.py` |
| Tests | `pipeline/tests/test_legacy_census.py` |

The archive itself (`/Volumes/OWC StudioStack HDD/DATA/ASTRO/legacy-archive`) is read-only to this
code. The census reads headers of `*.fz` files only; it opens no pixel, so it cannot speak to
saturation, comparison stars, the absolute clock or gain — the page says so wherever that matters.

Do not edit `notes/GO_NOGO_PREREGISTERED.md`: its sha256 is stored in the database at build time
(`census_meta.prereg_sha256`) and a unit test pins the gate constants to its text.
