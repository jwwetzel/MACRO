# Ground rules for every work package (read before touching anything)

Repo: `/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO`. Python: `/opt/miniconda3/envs/rlmt-checks/bin/python`
(use absolute paths; always `cd` to the repo in each shell command). Archives are READ-ONLY:
`/Volumes/OWC StudioStack HDD/DATA/ASTRO/{rlmt-archive,legacy-archive,Calibrations}`.

Binding context: `committee/reviews/2026-10-03/SYNTHESIS.md` and the seven memos beside it. Your package
names the finding ids you own; read those findings in the memos in full.

House standard (James's words): no shortcuts, no statistical bias untested or unaddressed; every script
self-contained, documented and commented; clean code; plots do the talking; every published number is
emitted by a script from a database, never typed. Match the existing code's style, docstring density and
test conventions (`pipeline/tests/`). Figures use `pipeline/macro_core/plotstyle.py`.

Concurrency — several packages run at once in this one working tree:
- Edit ONLY the files your package owns (listed in your brief). If you need a change elsewhere, write it
  down in your report; do not make it.
- Do NOT `git commit`, `git add`, `git stash`, `git checkout` or otherwise touch git state. The chair commits.
- Do NOT edit `pipeline/macro_core/project_plan.py`, `provenance.py`, `site.py` or run `build_site.py` /
  `update_project_plan.py set|sync|render` unless your brief says you own them.
- `products/manifest/rlmt-manifest.sqlite` is shared. Open it read-only (`file:...?mode=ro`) unless your
  brief says you own named tables in it. Never run `build_s0_manifest.py` unless you own S0.
- Run only your own test files while working; the chair runs the full suite.
- Never delete data. Never modify the archives.

Standing statistical rules (SYNTHESIS §5) apply to everything you produce.

Report: write `committee/work/<package>/REPORT.md` — what you did, the evidence (tables, figure paths,
numbers with their queries), each finding id with status CLOSED / PARTIAL / REBUTTED / BLOCKED and why,
every file you changed, ledger changes you want the chair to make, and anything that needs James.
Your final message is a ≤ 250-word summary of that report.
