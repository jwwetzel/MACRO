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
