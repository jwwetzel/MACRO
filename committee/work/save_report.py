#!/usr/bin/env python
"""Save a work package's final message as its REPORT.md.

Subagents are not permitted to write their own report files, so the chair
extracts each package's final assistant message from the task transcript and
files it verbatim.  Usage: save_report.py <transcript.jsonl> <package>
"""
import json, sys, pathlib
src, pkg = sys.argv[1], sys.argv[2]
last = None
for line in open(src):
    try:
        rec = json.loads(line)
    except ValueError:
        continue
    msg = rec.get("message", {})
    if msg.get("role") == "assistant":
        text = "".join(b.get("text", "") for b in msg.get("content", [])
                       if isinstance(b, dict) and b.get("type") == "text")
        if text.strip():
            last = text
out = pathlib.Path(__file__).parent / pkg / "REPORT.md"
out.write_text(last or "")
print(out, len(last or ""), "chars")
