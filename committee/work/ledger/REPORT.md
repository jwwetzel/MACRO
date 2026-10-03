# Package `ledger` — report

**The plan ledger now says what the committee ruled.** Every ADD / CHANGE / DROP / DEFER / HOLD /
CLOSE / UNBLOCK of SYNTHESIS §3–§4 is in `pipeline/macro_core/project_plan.py` as data, each with
its ruling and finding ids, and the ledger validates, cites and renders. One thing is left red for
the chair: **`test_site.py::test_the_stylesheet_defines_every_status_class_the_site_emits` fails**,
because the two new statuses need two chip rules in `docs/assets/macro.css`, which is not my file.
The fix is in `renderer_and_stylesheet.patch` beside this report (`git apply --check` passes).

My own test file is green: `pytest pipeline/tests/test_project_plan.py` → **286 passed** (was 114).

## 1. What was done, against the brief

| Brief item | Result |
|---|---|
| 1. Closing statuses | `dropped` and `deferred` exist; both require a `Ruling` (document, section, finding ids, date). Completion is `done / in-scope`; dropped and the 2027 backlog are counted and shown apart. |
| 2. Apply §3–§4 | 124 tasks carry a ruling: 58 ADD, 45 CHANGE, 10 DROP, 4 DEFER, 3 CLOSE, 2 HOLD, 2 UNBLOCK. CV_TimeSeries is reopened with CV-R1…R15. F-1…F-10 and G-1…G-5 are a new *Shared foundation* group. |
| 3. Blockers (U4) | No blocker says a table was destroyed. The eight tasks that carried that text are unblocked. Gates are computed: task dependencies (now cross-group), stage verdicts, and database probes. |
| 4. Legacy premise and name | Displayed title is "Legacy Archive 2015–2022 (census only)"; premise corrected per TE.F10; the key `Legacy_Rigel` and the `RIG-` ids are unchanged, so no link breaks. |
| 5. Strategies, ROADMAP, titles | Each strategy has `## 10. Committee amendments 2026-10-03`; ROADMAP §0 has a regenerated portfolio table and `### 0.1 Committee amendments 2026-10-03`; four skeleton titles re-emitted from the ledger. |
| 6. Tests, render | 286 tests green; `sync` then `render` run (44 pages, WEB recorded 2026-10-03T07:23:45Z). |

## 2. Evidence

### 2.1 The ledger after `sync` (emitted by `update_project_plan.py show`; statuses from `project_plan_status`)

| Group | done / in-scope | in progress | blocked | pending | redo needed | dropped | 2027 backlog | total |
|---|---|---|---|---|---|---|---|---|
| Shared_Foundation | 0/15 | 0 | 0 | 15 | 0 | 0 | 0 | 15 |
| TCrB_Monitoring | 0/41 | 3 | 1 | 34 | 3 | 2 | 2 | 45 |
| CV_TimeSeries | 0/49 | 0 | 1 | 14 | 34 | 0 | 0 | 49 |
| SN2023ixf_LightCurve | 0/26 | 1 | 2 | 16 | 7 | 1 | 0 | 27 |
| BeStar_Grism | 0/28 | 1 | 2 | 21 | 4 | 0 | 2 | 30 |
| DwarfGalaxy_AGN_Survey | 0/28 | 2 | 2 | 23 | 1 | 7 | 0 | 35 |
| Legacy_Rigel | 0/15 | 3 | 0 | 12 | 0 | 0 | 0 | 15 |

Every "done" reads `redo needed` because every stage in the DAG is STALE or STALE_UPSTREAM today
(22 + 17). That is ruling U5 showing through the plan, not a ledger fault: when F-7 and F-10 land,
`sync` rule 2 turns them back to `done` with no hand edit. Before this package the same table read
CV 34/34, T CrB 2/34, SN 5/20, Be 4/24, Dwarf 1/31, Legacy 0/11.

### 2.2 Rulings by project (from `pp.ruling_action` over the ledger)

| Group | ADD | CHANGE | DROP | DEFER | HOLD | CLOSE | UNBLOCK |
|---|---|---|---|---|---|---|---|
| Shared_Foundation | 15 | | | | | | |
| TCrB_Monitoring | 10 | 11 | 2 | 2 | 1 | 1 | 2 |
| CV_TimeSeries | 15 | | | | | | |
| SN2023ixf_LightCurve | 7 | 7 | 1 | | | 2 | |
| BeStar_Grism | 3 | 13 | | 2 | 1 | | |
| DwarfGalaxy_AGN_Survey | 4 | 9 | 7 | | | | |
| Legacy_Rigel | 4 | 5 | | | | | |

The per-task tables (id, verb, ruling, acceptance criterion) are in each strategy's §10, emitted
by `python pipeline/scripts/update_project_plan.py amendments <Project>`.

### 2.3 U4 by query, not by sentence (`update_project_plan.py blockers`, read-only)

| Probe | Value now | Meaning |
|---|---|---|
| `SELECT count(*) FROM s2_linearity_ladders` | 66 | gate met — the table a blocker called destroyed |
| `SELECT count(*) FROM s2_ptc_fits` | 4 | gate met |
| `SELECT count(*) FROM s2_ceiling_modes` | 8 | gate met |
| canonical frames that are `_calibrated` twins (F-1) | 27,261 | acceptance probe: to go until F-1 lands |
| `frames` has a `focpos` column (F-2) | 0 | to go |
| `mech_epoch` table exists (F-3) | 0 | to go |
| NGC 5548 slot-'6' frames not measured dispersed (DW-P02-slot6) | 52 | to go until F-6 re-issues verdicts |

`blockers` exits non-zero if a blocked task's database gates are all met and nothing else holds it,
or if a `done` task's acceptance probe is unmet. Today it reports zero contradictions. Eight tasks
remain `blocked`, all external or held: TCRB-A9-profiles and BE-VR-hold (HOLD pending D1),
CV-R15 (James), SN-S1-narrowband-curves and SN-S6-halpha-curve (transmission curve, dated
2026-11-15), BE-S2-calibration (re-opening + rev. 3), DW-P02-filter-dossier and
DW-P2-band-transformations (Cannon). 32 project tasks now wait on a foundation task by dependency.

### 2.4 Citations

`pp.verify_citations(repo)` → `()` and `pp.verify_rulings(repo)` → `()`: every task source, every
ruling's section of the synthesis, and every finding id (resolved against the seat's memo) check
out. A ruling citing a finding nobody filed is a test failure.

### 2.5 Rendered pages

`docs/<Project>/index.html` × 6, `docs/evidence.html`, `index.html` and the site chrome were
regenerated. Mastheads now read, e.g., "Cataclysmic-Variable Time Series (reopened) … 0 of 49 plan
tasks complete (34 need redoing, 1 blocked) · AJ — major revision open". The hub reads "0 of 187
plan tasks complete across six projects". No page contains "tables were destroyed".

## 3. Findings

| Finding | Status | Why |
|---|---|---|
| SYNTHESIS §0 (three ways to close) | CLOSED | Statuses, mandatory ruling, in-scope count, separate phases; enforced by `validate()` and `set`. |
| U4 · DS.F6 · ED.E8 · TE.F11 · DE.F6 · RF cross-cutting 5 (stale blockers) | CLOSED in the ledger and CLI; **PARTIAL on the pages** | Dependency gates show on the pages today. Stage-verdict and probe gates, and the ruling line, show only after the chair applies `renderer_and_stylesheet.patch`. |
| DS.F5 / F-7, ledger half (`show` prints an unsynced count) | CLOSED | `show` states its sync age; `show --live` prints what `sync` would leave; `sync` was run. The provenance half is `foundation-s0`'s. |
| ED.E6 (titles, ROADMAP, venue lines) | CLOSED, one exception | CV skeleton title is not mine to edit. |
| ED.E7 (NGC 5548 still planned) | CLOSED | P41/42/44/45 dropped; P03 re-scoped; NGC 5548 out of title and claim. |
| TE.F10 (Legacy premise and name) | CLOSED in the ledger; **PARTIAL on the site** | `site.SHORT_LABEL` still says "Legacy Rigel" (6 occurrences on the landing page). |
| DE.F3, ledger note (CV-P15 "measured linearity curves") | CLOSED | Task text now says it is not a measurement near any veto. |
| RF minor 5 (injection night 02-28 vs 02-27) | CLOSED | Title corrected; the strategy's date is the same night by UTC. |
| U1, U2, U3, U6, U7, U8, U9 and D1–D3 as plan content | CLOSED (recorded) | The work itself belongs to other packages. |

## 4. Decisions the chair should know I made

1. **The shared foundation is a group beside `PROJECTS`, not a seventh project.** The site's
   navigation, `provenance._PROJECT_DIRS` and the "across six projects" counts are keyed on
   `PROJECTS`, and none of those files is mine. `pp.groups()` / `pp.ledger_tasks()` include it;
   `pp.all_tasks()` stays the six papers. It shows in `show` and `blockers`; it has no HTML page
   until the patch is applied (the patch adds a section to `docs/evidence.html`).
2. **`sync` gained a third rule.** A status recorded before a ruling's date gives way to the
   ledger status that ruling set. Without it the `in_progress` recorded on 2026-08-21 would have
   overridden the committee's closure of SN-G0c and SN-venue-decision. A status recorded on or
   after the ruling date is never touched. `PLAN_CODE_VERSION` is now v2.0.
3. **I ran `sync` for real** (49 rows at 2026-10-03T07:01:08Z). If `foundation-s0`'s S0 table swap
   dropped them, re-running `sync` reproduces them exactly.
4. **DW-P52 (transient search) is `dropped`**, with its re-opening condition in the text; the
   synthesis says "unless threshold/trials pre-declared with injection".
5. **BE-S2-calibration stays `blocked`, not `deferred`.** It needs new frames, but existing era-C
   science depends on it, so it is on the critical path. TCRB-P0-restart, by contrast, is deferred.
6. **Ids I had to mint** where the synthesis named none: TCRB-A5a-detection-rule, TCRB-A9-profiles,
   TCRB-B0-peak-census, BE-VR-hold, BE-S0-era-table, BE-S0-dispositions, BE-X1-dither-test,
   BE-X2-season2-observing, SN-G0-rerun, SN-G0d-s2c-broadband, SN-S5-template-table,
   SN-S5b-peak-epoch, SN-S6-0-excess-gate, SN-S6a-flash-colour, SN-S7b-residuals-table,
   RIG-L0-dedup-reconcile, RIG-L0-mech-epoch, RIG-L1-clock-audit, RIG-L2-prereg.
7. **SN-G0c and SN-venue-decision are re-bound to stage `SN-G0`** (from `G` and `S4`): both verdicts
   are rows of `sn_g0_verdict`.
8. **New tasks are `pending` even where a sibling package is working on them now.** I cannot see
   their state; see §6.
9. I ran `tests/test_site.py` once (read-only, 20 min) to check the blast radius of the new
   statuses. It is not my file. 48 passed, 1 failed as above.

## 5. Files changed

- `pipeline/macro_core/project_plan.py` — statuses, `Ruling`, `Probe`, sync rule 3, computed
  gates, `FOUNDATION`, all six plans, `paper_title`, title helpers.
- `pipeline/scripts/update_project_plan.py` — `show` (in-scope count, sync age, `--live`), `set
  --ruling`, `sync` over the whole ledger, new `blockers`, `amendments`, `titles`.
- `pipeline/tests/test_project_plan.py` — sections 13–19 (172 more test cases, 114 → 286); integrity sweep widened
  to the whole ledger.
- `{TCrB_Monitoring,CV_TimeSeries,SN2023ixf_LightCurve,BeStar_Grism,DwarfGalaxy_AGN_Survey}/ANALYSIS_STRATEGY.md`
  — banner, §10; CV's "Target journal" line.
- `ROADMAP.md` — §0 table and §0.1.
- `manuscripts/{TCrB_Monitoring,SN2023ixf_LightCurve,BeStar_Grism,DwarfGalaxy_AGN_Survey}/main.tex`
  — the `\title` line only (untracked directory; originals are in my scratchpad).
- Generated by `render`: `docs/**` pages and `index.html`.
- Manifest: 49 rows appended to `project_plan_status`; one `WEB` row in `stage_provenance`.
- This directory: `REPORT.md`, `renderer_and_stylesheet.patch`, `tools/` (the scripts that wrote
  §10 and ROADMAP §0).

## 6. For the chair

**Apply (not my files):**
- `git apply committee/work/ledger/renderer_and_stylesheet.patch` — adds `.chip.dropped` /
  `.chip.deferred` to `docs/assets/macro.css` (turns the failing site test green) and, in
  `report_projects.py`, the dropped/backlog counts beside the fraction, computed gates and the
  ruling in each task row, and a Shared-foundation table on the hub. Rendered and checked on a
  scratch copy; not run against the test suite.
- `site.py`: `SHORT_LABEL["Legacy_Rigel"]` → "Legacy Archive"; `Build.open_tasks` lists every task
  that is not `done`, so dropped and deferred tasks appear as open on the Case pages — it should
  exclude `pp.CLOSED_BY_RULING`.
- `provenance.py`: the `STRAT` stage is hand-authored and now stale (five strategies and ROADMAP
  edited); `committee/reviews/2026-10-03/SYNTHESIS.md` is cited by 124 tasks and is not a declared
  resource.
- `manuscripts/CV_TimeSeries/main.tex`: title and venue are the integrating package's.

**Record as packages report** (`update_project_plan.py set <id> <status> --evidence …`):
`foundation-s0` → F-1, F-2, F-3, F-7, F-9 · `detector` → F-4, F-5, TCRB-P0-gain-ptc,
TCRB-B0-peak-census, SN-S2-linearity · `grism` → G-1…G-5, F-6, then D1 decides TCRB-A9-profiles
and BE-VR-hold (open or drop) · `clock` → F-8, CV-R11 · `cv-stats` → CV-R1…R8 · `cv-literature` →
CV-R12 · `novelty-tcrb` → TCRB-N1 · `novelty-be` → BE-S-1a-bess, BE-N1-gate, the QQ Gem half of
BE-S0-dispositions · `novelty-dwarf` → DW-P31, DW-N1, inputs to DW-P36-0 · `legacy` → RIG-L2-prereg
and the L0/L1 tasks · `ops` → TCRB-P0-eruption-block. Then `sync` and `render`.

**Dated:** on 2026-11-15, if no transmission curve is on file, set SN-S1-narrowband-curves and
SN-S6-halpha-curve to `dropped` (the ledger already carries the ruling, so `set` will accept it).

## 7. Needs James

- D4 (instrument paper) — the foundation group's venue line says so.
- CV-R15: authorship policy, ORCIDs, Zenodo DOI, outside reader.
- Contact Cannon — DW-P02-filter-dossier and DW-P2-band-transformations are blocked on it.
- Send ops rev. 3 — BE-S2-calibration and the T CrB bench-dark question wait on it.
- The four re-emitted skeleton titles are my wording of the committee's scope; the Dwarf one is
  provisional by the strategy's own rule (fixed only after the literature cross-match).
