# Committee brief — plan review of 2026-10-03

Repo: `/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO` (git: jwwetzel/MACRO).
Python: `/opt/miniconda3/envs/rlmt-checks/bin/python` (always use absolute paths; `cd` to the repo first).
Archives (read-only): `/Volumes/OWC StudioStack HDD/DATA/ASTRO/{rlmt-archive,legacy-archive,Calibrations}`; catalog `rlmt-catalog.sqlite` there.

## State of play
- Six projects: TCrB_Monitoring (2/34 tasks), CV_TimeSeries (34/34, 32-pp draft through six referee rounds),
  SN2023ixf_LightCurve (5/20), BeStar_Grism (4/24), DwarfGalaxy_AGN_Survey (1/31), Legacy_Rigel (0/11, undecided candidate).
- Shared pipeline: every stage currently reads STALE after the frame-dispersion (grism vs imaging) and
  compressed-FITS geometry reclassifications. S1 was recomputed on measured dispersion; SN Gate 0 has not been re-run.
- Winer Observatory was closed for monsoon; re-opening is expected this month (October 2026). No new frames yet.
- Legacy (pre-MACRO, 2015–2022) archive is on disk but not yet in the manifest.

## Where things are
- `ROADMAP.md`, `<Project>/ANALYSIS_STRATEGY.md` — the strategies (hand-authored, currently flagged STALE).
- Plan ledger: `pipeline/macro_core/project_plan.py`; print with `python pipeline/scripts/update_project_plan.py show`.
- Freshness: `python pipeline/scripts/check_pipeline_status.py`.
- Evidence pages: `docs/pipeline/*.html`, `docs/<Project>/*.html`. Code: `pipeline/`. Tests: `pipeline/tests/`.
- Databases: `products/manifest/rlmt-manifest.sqlite` (+ others under `products/`). Query with sqlite3 read-only.
- CV draft: `manuscripts/CV_TimeSeries/main.tex`, `main.pdf`. Ops request: `ops/2026-08_observatory_request.md`.

## What to deliver
Write ONE memo to `committee/reviews/2026-10-03/<your-seat>.md` (≤ 2,500 words), structured:

1. **Existing work** — findings on what is already built/claimed, each with severity, file:line evidence, and the closing test.
2. **Plan hardening, per project** (all six) — tasks to ADD, CHANGE, or DROP, each with a one-line acceptance criterion.
   Say explicitly which planned tasks are impossible or pointless with the data on disk, and what the honest reduced scope is.
3. **Cross-cutting** — your top five items for the shared pipeline or the portfolio.
4. **Verdict** — per project: ready to execute / execute with amendments / re-scope / stop.

Rules: read-only. Do not modify code, databases, manuscripts or docs. Do not launch long pipeline runs
(short read-only queries and opening a handful of FITS frames are fine). Do not write anywhere but your memo.
Your final message back should be a ≤ 200-word summary of your top findings.
