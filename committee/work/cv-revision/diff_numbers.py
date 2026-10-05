#!/usr/bin/env python
"""CV-R14 — macro-by-macro diff of the CV paper's numbers.tex.

Compares three snapshots of ``manuscripts/CV_TimeSeries/numbers.tex``:

  A  ``committee/work/cv-stats/numbers.tex.before-cv-stats``
     the macro file the first draft was built on (307 macros);
  B  ``committee/work/cv-revision/numbers.tex.before-cap``
     after the clean-tree rebuild (F-10) and the cv-stats emission, before
     the S2 linearity cap was applied to the photometry;
  C  ``manuscripts/CV_TimeSeries/numbers.tex``, the current file.

Every macro that changed value, appeared or disappeared between A and C is
listed with the step in which it moved (A->B or B->C) and the reason.  The
reason is assigned by rule from the step and the macro's provenance, and a
macro no rule explains is printed as UNEXPLAINED so a reader sees it.

Writes ``committee/work/cv-revision/numbers_macro_diff.md``.  Read-only
otherwise.  Usage:

    /opt/miniconda3/envs/rlmt-checks/bin/python committee/work/cv-revision/diff_numbers.py
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
A = REPO / "committee" / "work" / "cv-stats" / "numbers.tex.before-cv-stats"
B = HERE / "numbers.tex.before-cap"
C = REPO / "manuscripts" / "CV_TimeSeries" / "numbers.tex"
OUT = HERE / "numbers_macro_diff.md"

LINE = re.compile(r"\\newcommand\{\\(\w+)\}\{(.*)\}(?:\s*%\s*(.*))?$")


def load(path: Path) -> dict:
    """``macro -> (value, provenance comment)``."""
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = LINE.match(line.rstrip())
        if m:
            out[m.group(1)] = (m.group(2), (m.group(3) or "").strip())
    return out


#: Reasons for the A->B step (2026-10-03 cv-stats emission and the
#: 2026-10-04 clean-tree rebuild, F-10), by macro.
REASON_AB = {
    "NumCeilingFramesRange": "S2 re-run (detector package): the ceiling "
                             "histograms now cover up to 299 frames per mode",
    "NumTieSeries": "CV-S6 re-run in the clean rebuild: the 1 MHz y block "
                    "now reaches the tie stage and fails it (counted as a "
                    "primary block)",
    "NumTieUntied": "same: that block is the second untied primary block",
    "NumUntiedBlockNoTieStage": "same: no block now fails BEFORE the tie "
                                "stage, so the clause has no subject",
    "NumUntiedBlockNoTieStageRows": "same",
    "NumExtinctionBoundRangeMmag": "CV-S8 re-run in the clean rebuild on the "
                                   "re-tied series",
    "NumStLmiFittedPeriodSigmaD": "CV-S9 v1.1: the one-sided max(chi2_nu, 1) "
                                  "error rescaling was removed (standing "
                                  "rule 1)",
    "NumStLmiPeriodAgreementSigma": "same",
    "NumBandPairsSignificant": "CV-S9 v1.1: band pairs use scatter-based "
                               "errors (CV-R1); the g-i pair is now "
                               "significant",
    "NumSuperoutburstAmp": "literature correction (PH.P5): 4.0 mag "
                           "superoutburst amplitude of YZ Cnc, kato2002",
    "NumSuperoutburstAmpThreshold": "added when NumSuperoutburstAmp became "
                                    "the literature amplitude: the CV-S7 "
                                    "classification threshold (3.0 mag) is "
                                    "kept as its own constant",
    "NumMacrosTotal": "count of macros: the revision added the NumRv set",
    "NumMacrosPhot": "same",
    "NumMacrosExternal": "same (literature scales and method constants)",
}


def reason_bc(name: str, prov: str) -> str:
    """Why a macro moved between B and C (the linearity-cap re-run)."""
    if name.startswith("NumMacros"):
        return "count of macros: new checks stages (reduction, ramp, mech, " \
               "clock, cap, detector) added macros"
    if "rv_result" in prov or name.startswith("NumRv"):
        return ("revision statistics recomputed on the re-run chain (the "
                "photometry now withholds measurements above the S2 1% "
                "linearity cap; the rv stages had last run before the "
                "F-10 rebuild of CV-S9, so they also absorb that rebuild)")
    if any(t in prov for t in ("cv_", "p2_", "p3_", "p4_", "p5_", "ch_")):
        return ("photometry chain re-run with the S2 1% linearity cap "
                "(run_cv_photometry.py recap): comparison and target "
                "measurements above the cap are withheld, so ensembles, "
                "ties and every downstream product move")
    if "s2_" in prov or "detector_params" in prov:
        return "detector package re-emitted S2 tables"
    return ""


def main() -> None:
    a, b, c = load(A), load(B), load(C)
    rows = []
    for name in sorted(set(a) | set(c)):
        va = a.get(name, (None, ""))[0]
        vb = b.get(name, (None, ""))[0]
        vc, prov = c.get(name, (None, ""))
        if va == vc:
            continue
        if va is None:
            kind = "added"
        elif vc is None:
            kind = "removed"
        else:
            kind = "changed"
        steps = []
        if va != vb:
            steps.append(("A->B", REASON_AB.get(name) or
                          ("added by the cv-stats revision emission"
                           if va is None and name.startswith("NumRv") else
                           "")))
        if vb != vc:
            steps.append(("B->C", reason_bc(name, prov)))
        rows.append((name, kind, va, vb, vc, steps))
    n_un = sum(1 for r in rows for s in r[5] if not s[1])
    kinds = Counter(r[1] for r in rows)
    lines = [
        "# CV paper: macro-by-macro diff of numbers.tex (CV-R14)", "",
        "Emitted by `committee/work/cv-revision/diff_numbers.py`. Do not "
        "edit.", "",
        f"- A: `{A.relative_to(REPO)}` ({len(a)} macros) — the first draft",
        f"- B: `{B.relative_to(REPO)}` ({len(b)} macros) — after the "
        "F-10 rebuild and the cv-stats emission",
        f"- C: `{C.relative_to(REPO)}` ({len(c)} macros) — now", "",
        f"Changed {kinds['changed']}, added {kinds['added']}, removed "
        f"{kinds['removed']}; unexplained steps: {n_un}.", "",
        "## Pre-existing macros that changed or disappeared", "",
        "| macro | A | B | C | step: reason |", "|---|---|---|---|---|"]
    for name, kind, va, vb, vc, steps in rows:
        if kind == "added":
            continue
        why = "; ".join(f"{s}: {r or 'UNEXPLAINED'}" for s, r in steps)
        lines.append(f"| `{name}` | {va} | {vb} | {vc} | {why} |")
    added = [r for r in rows if r[1] == "added"]
    moved_bc = [r for r in added if r[3] is not None and r[3] != r[4]]
    new_bc = [r for r in added if r[3] is None]
    lines += ["", "## Macros added by the revision", "",
              f"{len(added)} macros did not exist in A. Of these, "
              f"{len(added) - len(new_bc)} existed in B; {len(moved_bc)} of "
              f"those changed value between B and C (listed below), and "
              f"{len(new_bc)} are new in C (the checks stages).", "",
              "### Changed between B and C", "",
              "| macro | B | C |", "|---|---|---|"]
    lines += [f"| `{r[0]}` | {r[3]} | {r[4]} |" for r in moved_bc]
    lines += ["", "### New in C", "", ", ".join(f"`{r[0]}`" for r in new_bc)]
    lines += ["", "Reason for every B->C change of an `rv_result` macro: "
              + reason_bc("NumRvX", "rv_result"), ""]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{kinds['changed']} changed, {kinds['added']} added, "
          f"{kinds['removed']} removed; {len(moved_bc)} revision macros "
          f"moved B->C; unexplained steps {n_un} -> {OUT}")


if __name__ == "__main__":
    main()
