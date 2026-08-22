#!/usr/bin/env python3
"""build_site.py — assemble the site's three navigation layers.

WHAT THIS DOES

    python pipeline/scripts/build_site.py

writes the landing page, and for the shared pipeline and each of the six
projects a Case, a Figures wall, a Draft Paper page where a compiled draft
exists, and an Evidence Detail index; then puts the same chrome — the project
row, the view tabs, the breadcrumb, the sticky question rail, previous/next —
around every page the thirteen report renderers already wrote.

It adds no facts.  Every question, figure, caption, progress fraction and
freshness verdict on the pages it writes is read out of the plan ledger, the
provenance DAG, the products databases, or the evidence pages' own section
headings.  See ``macro_core.site`` for where each one comes from.

WHEN TO RUN IT

* after any report renderer re-runs (a fresh page arrives without chrome);
* after ``update_project_plan.py set`` changes a status;
* ``update_project_plan.py render`` already calls it, so the normal working
  rhythm needs no extra step.

Running it twice produces the same bytes as running it once: each wrap
recovers the renderer's own markup from between the content markers before
re-wrapping, so chrome never nests inside chrome.

``--check`` writes nothing and exits non-zero if a rebuild would change any
page — the form a pre-commit hook or CI wants.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = SCRIPT_DIR.parent
REPO_ROOT = PIPELINE_DIR.parent
sys.path.insert(0, str(PIPELINE_DIR))

from macro_core import site                                  # noqa: E402

DEFAULT_MANIFEST = REPO_ROOT / "products" / "manifest" / "rlmt-manifest.sqlite"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST,
                    help="the manifest database, for live provenance "
                         "verdicts and recorded plan status")
    ap.add_argument("--docs", type=Path, default=REPO_ROOT / "docs")
    ap.add_argument("--check", action="store_true",
                    help="write nothing; exit 1 if a rebuild would change "
                         "any page")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    if args.check:
        return _check(args)

    written = site.build_site(manifest=args.manifest, docs_dir=args.docs,
                              repo_root=REPO_ROOT)
    if not args.quiet:
        for path in written:
            try:
                shown = path.relative_to(REPO_ROOT)
            except ValueError:                            # pragma: no cover
                shown = path
            print(f"wrote {shown}")
        print(f"\n{len(written)} page(s). Three layers on every one: the "
              f"project row, the view tabs, the question rail.")
    return 0


def _check(args) -> int:
    """Rebuild into memory and report any page whose bytes would move.

    The build writes files, so "check" means: snapshot, build, compare,
    restore.  Restoring is what makes this safe to run on a clean tree — a
    check that left the tree dirty would be a check nobody runs twice.
    """
    docs = Path(args.docs)
    before = {p: p.read_bytes() for p in docs.rglob("*.html")}
    written = site.build_site(manifest=args.manifest, docs_dir=docs,
                              repo_root=REPO_ROOT)
    after = {p: p.read_bytes() for p in docs.rglob("*.html")}
    changed = sorted(
        {p for p in set(before) | set(after)
         if before.get(p) != after.get(p)})
    for path, data in before.items():                 # restore the snapshot
        if after.get(path) != data:
            path.write_bytes(data)
    for path in set(after) - set(before):
        path.unlink()
    del written
    if changed:
        print(f"{len(changed)} page(s) would change:")
        for path in changed:
            print(f"  {path.relative_to(REPO_ROOT)}")
        print("\nrun: python pipeline/scripts/build_site.py")
        return 1
    print("site is current.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
