#!/usr/bin/env python
"""Re-run every stale pipeline stage, in DAG order, unattended.

Reads the ordered re-run plan from ``check_pipeline_status.py plan`` and
executes each stage's commands exactly as the plan prints them, then records
the stage.  Hand-authored stages (whose plan line is a parenthesised note
rather than a command) are recorded with an explicit note saying they were
re-recorded, not re-derived, so the provenance never claims work that was not
done.  Stops at the first failing command — a later stage must never be
recorded FRESH on top of an earlier one that failed.

Usage:  refresh_stale.py [--dry-run] [--only=STAGE,STAGE]      (log: products/refresh_stale.log)
"""
import re
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PY = sys.executable
STATUS = ["pipeline/scripts/check_pipeline_status.py"]
QUEUE_BUILDERS = ("run_s2c_dispersion.py build",)
HAND_NOTE = ("re-recorded after the 2026-10-04 S0 rebuild; hand-authored "
             "content NOT re-derived (committee reviews benched)")


def plan():
    """Return [(stage, [command, ...]), ...] from the status tool's plan."""
    out = subprocess.run([PY, *STATUS, "plan"], cwd=REPO, text=True,
                         capture_output=True).stdout
    stages, cur = [], None
    for line in out.splitlines():
        m = re.match(r"\s*\d+\.\s+(\S+)\s+\[", line)
        if m:
            cur = (m.group(1), [])
            stages.append(cur)
        elif cur and line.strip().startswith("$ "):
            cur[1].append(line.strip()[2:])
    return stages


def main() -> int:
    dry = "--dry-run" in sys.argv
    # --only S4,SN-G0,... runs just those stages (a second, parallel driver
    # for branches of the DAG that do not wait on the first driver's stage).
    only = next((a.split("=", 1)[1].split(",") for a in sys.argv
                 if a.startswith("--only=")), None)
    for stage, cmds in plan():
        if only and stage not in only:
            continue
        for cmd in cmds:
            if cmd.startswith("("):
                continue                      # hand-authored: no command
            if any(k in cmd for k in QUEUE_BUILDERS):
                # These rebuild a resumable measurement queue from scratch and
                # (rightly) refuse while measured rows exist.  A refresh keeps
                # the measurements and lets the following `run` finish
                # whatever is pending; a full re-measure is a deliberate act.
                print(f"[skip] {stage}: {cmd}  (queue exists; keep measurements)",
                      flush=True)
                continue
            if " record " in cmd and "<" in cmd:
                cmd = re.sub(r"--note '.*'", f"--note '{HAND_NOTE}'", cmd)
            if " record " in cmd:
                # A re-run is a NEW run: record it under its own timestamp
                # rather than overwriting the fingerprints of the earlier one.
                cmd += " --run-utc " + time.strftime("%Y-%m-%dT%H:%M:%S+00:00",
                                                     time.gmtime())
            cmd = re.sub(r"^python\b", f"'{PY}'", cmd)
            t0 = time.time()
            print(f"[{time.strftime('%H:%M:%S')}] {stage}: {cmd}", flush=True)
            if dry:
                continue
            rc = subprocess.run(cmd, shell=True, cwd=REPO).returncode
            print(f"    -> exit {rc} in {(time.time() - t0) / 60:.1f} min",
                  flush=True)
            if rc != 0:
                print(f"STOPPED at {stage}: a later stage must not be recorded "
                      "on top of a failed one.", flush=True)
                return rc
    print("ALL STALE STAGES REFRESHED", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
