# INCIDENT — worker processes of other jobs killed by the `legacy` package (2026-10-03 02:08:28 CDT)

**What happened.** While restarting my own header scan I ran `pkill -f "multiprocessing.spawn"` to clear
my scan's orphaned pool workers. That pattern is not specific to my job: on macOS every Python
`multiprocessing` worker has that command line. The command therefore killed the pool workers of
**every** Python job running under this user at 02:08:28, not only mine. This was my error.

**Jobs affected (parents were NOT killed; their workers were):**

| PID | Job | Owner | Observed after |
|---|---|---|---|
| 80298 | `legacy_dedupe.py --workers 14` (in `/Volumes/OWC StudioStack HDD/DATA/ASTRO`) | James / chair (pre-existing job, started 00:51) | parent alive, **no workers left**; log stops at "20,000/202,961 checked" |
| 42324 | `pipeline/scripts/build_s3b_clock_transits.py --stage photometry --workers 14` | package `clock` | 14 new workers respawned at 02:08:28; tasks in flight at that instant were lost |
| 39944 | `pipeline/scripts/run_g_dispersion.py --extract --workers 4` | package `grism` | 4 new workers respawned at 02:08:28–31; tasks in flight were lost |

**What this means for each owner.**

- `legacy_dedupe.py`: it deletes an uncompressed `.fts` only in the parent, after a worker returns
  "redundant", so no file is half-deleted and nothing was removed wrongly. It is idempotent: re-running
  it re-walks and continues. **It needs to be restarted by its owner** (I have not restarted it — the
  archive is read-only to this package and the job deletes files). Until it is, the archive still holds
  ≈183k redundant uncompressed twins.
- `clock` and `grism`: if the job uses `multiprocessing.Pool.map/imap`, a task whose worker was killed
  never returns and the job can **hang at the end** instead of failing; if it uses
  `ProcessPoolExecutor`, it will have raised `BrokenProcessPool`. Either way the run that was in
  progress at 02:08:28 should be treated as interrupted and re-run (both scripts are expected to be
  resumable under the house standard).

No file, database or git state was modified by the kill itself.
