#!/usr/bin/env python
"""Reduce a named grism sample into the 1-D spectrum cache.

One job: pixels -> ``products/grism/spec1d/*.npz`` + a ``g_frames`` row
per frame (trace geometry, header regressors, quality numbers).  Every
other G runner reads the cache; none of them touches the archive.

    run_g_reduce.py --sample tcrb               # the 2025 T CrB series
    run_g_reduce.py --sample tcrb_slot          # slot '6'/'W' T CrB
    run_g_reduce.py --sample identity_controls  # non-target controls
    run_g_reduce.py --list                      # sample sizes, no work

Resumable: frames already reduced at the current code version are
skipped.  The archive and the manifest are opened read-only.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from macro_grism import config as gconfig                # noqa: E402
from macro_grism import db as gdb                        # noqa: E402
from macro_grism import samples as gsamples              # noqa: E402
from macro_grism import store as gstore                  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--manifest", default=str(gstore.default_manifest()))
    ap.add_argument("--archive", default=str(gstore.DEFAULT_ARCHIVE))
    ap.add_argument("--db", default=str(gconfig.GRISM_DB))
    ap.add_argument("--sample", action="append", default=[],
                    choices=sorted(gsamples.SAMPLES))
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args(argv)

    mcon = gdb.connect_manifest_ro(args.manifest)
    if args.list:
        for name, fn in sorted(gsamples.SAMPLES.items()):
            print(f"  {name:18s} {len(fn(mcon)):5d} frames")
        return 0
    con = gdb.connect_grism(args.db)
    for name in args.sample:
        tasks = gsamples.SAMPLES[name](mcon)
        if args.limit:
            tasks = tasks[:args.limit]
        n = gstore.reduce_batch(con, tasks, Path(args.archive),
                                workers=args.workers)
        print(f"{name}: {len(tasks)} frames in sample, {n} reduced now")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
