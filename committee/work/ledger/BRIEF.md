# Package `ledger` — make the plan ledger say what the committee ruled
Owns: `pipeline/macro_core/project_plan.py`, `pipeline/tests/test_project_plan.py`, `pipeline/scripts/update_project_plan.py`,
`*/ANALYSIS_STRATEGY.md`, `ROADMAP.md`, the title lines of the four skeleton manuscripts (NOT `manuscripts/CV_TimeSeries/`),
and `docs/Legacy_Rigel` naming via the ledger. You may run `update_project_plan.py show|sync|render`.
Do:
1. Add the two closing statuses of SYNTHESIS §0 (`dropped`, `deferred`) with mandatory ruling citation; completion count = done / (in-scope), with dropped shown separately and deferred listed as a "2027 backlog".
2. Apply every ADD / CHANGE / DROP / DEFER / HOLD / CLOSE / UNBLOCK in SYNTHESIS §3–§4 (and the detail behind each in the memos) to the ledger: new tasks get strategy citation = the committee finding, a one-line acceptance criterion, and honest status. CV_TimeSeries is REOPENED: CV-R1…R15 are open tasks; its headline must no longer read "complete". Add the shared-foundation tasks F-1…F-10, G-1…G-5 as a pipeline area if the ledger has one, else as a new "Shared foundation" group.
3. Replace prose blockers that contradict the database (U4) with blockers computed from a query or DAG verdict where feasible; at minimum remove every "S2's tables were destroyed" string and state the true blocker.
4. Correct the Legacy project's premise and name (TE.F10) — keep the directory id stable if renaming would break links; change the displayed title.
5. Amend each ANALYSIS_STRATEGY.md and ROADMAP §0 with a dated "Committee amendments 2026-10-03" section recording the rulings that bind it (venues, scopes, dropped analyses, figure caps), and regenerate skeleton titles (ED.E6).
6. Tests updated and green for your files. Then `render`.
