# Goal-mode rules (2026-10-04) — supersede GROUND_RULES.md where they differ

Goal (James): every in-scope task of CV_TimeSeries (49) and SN2023ixf_LightCurve (26) closed.
Committee reviews are BENCHED: do not spawn reviewers. Token budget is tight: no exploration
beyond what a task needs, no extra figures, terse reports.

- Read `committee/work/GROUND_RULES.md` (house standard, concurrency) and the SYNTHESIS §4 entry
  for your project. Python: `/opt/miniconda3/envs/rlmt-checks/bin/python`; `cd` to the repo each command.
- You MAY close your own tasks in the ledger as you finish them:
  `python pipeline/scripts/update_project_plan.py set <TASK> done --evidence <repo-relative path> --note '<one line, numbers>'`
  A task is done only when its acceptance criterion (in `update_project_plan.py show <PROJECT>`) is met.
  If a task is genuinely impossible on the data that exist, set it `blocked` with the true reason —
  never mark done what is not done. Never touch other projects' tasks.
- Every published number is emitted by a script from a database. Pipeline stages you re-run must be
  recorded: `python pipeline/scripts/check_pipeline_status.py record <STAGE> --run-utc <now ISO>`.
- Do not git commit/push; the chair commits. Do not run `build_site.py` or `update_project_plan.py render|sync`.
- Report = your final message (≤300 words): tasks closed with one-line evidence each, tasks left and why,
  files changed.

## Goal 2 (2026-10-05): T CrB, Be-Star, Dwarf, Legacy — all steps

- Same rules as above. Committee reviews stay benched EXCEPT standing rule 5: when your paper's abstract
  and figure set are ready, ask Seat 6 for approval by SendMessage to 'a211ee0ccc2f88799' (the journal
  editor seat; it writes its decision under committee/reviews/), then make its edits.
- The grism library (G-1…G-5, D1) is being rebuilt by the grism worker. Projects that need it (T CrB A-track,
  Be S3–S13) do everything else first and poll the ledger (`update_project_plan.py show <PROJECT>`, every
  ≤10 min, sparingly) for G-1..G-5 to read done before extracting spectra.
- External data: AAVSO blocks scripted access — use ARAS, ASAS-SN, ZTF, TESS (MAST), BeSS, VizieR.
  Never invent references. If a task truly needs something only James or an outside party can supply
  (e.g. a filter transmission curve, a dry-run log from the site), set it `blocked` with that reason.
- Archive frames are READ-ONLY. The manifest is shared: write only your project's tables.
