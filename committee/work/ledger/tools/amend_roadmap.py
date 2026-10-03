"""ROADMAP §0: regenerate the portfolio table from the ledger's current
titles and venues, and add the dated committee-amendments subsection."""
import sys
sys.path.insert(0, ".")
sys.path.insert(0, "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO/pipeline")
from patchlib import Patcher
from macro_core import project_plan as pp

REPO = "/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO"
p = Patcher(f"{REPO}/ROADMAP.md")

# One line per project that is NOT in the ledger: the headline product and
# the committee's verdict.  Title and venue come from the ledger itself.
HEADLINE = {
    "TCrB_Monitoring": (
        "Hα equivalent-width series through the 2025 pre-eruption dip and "
        "recovery (slitless grism, Feb–Jun 2025), with ARAS cross-validation "
        "and B-only photometric anchors; line flux as EW × continuum",
        "Execute with amendments — highest priority. Novelty table first; "
        "profile/velocity work held pending D1; flickering is paper 2"),
    "CV_TimeSeries": (
        "ST LMi colour–phase curve in two instrument eras, a spin/spot-"
        "longitude timing bound, a YZ Cnc excluded-amplitude statement, "
        "a coverage audit of three more polars, and the release",
        "**Major revision — reopened.** Was 34/34; CV-R1 … CV-R15 block "
        "release; D3 (band offset: astrophysics or estimator bias) settled "
        "by test"),
    "SN2023ixf_LightCurve": (
        "gri validation and limits release, +5.4 → +50 d, as residuals "
        "against the published world dataset; flash-phase Hα excess from a "
        "differential colour only if it passes a predicted-excess gate",
        "Re-scope (U7). Grism NOT PROMOTED; model fit dropped; Gate 0 "
        "re-run on a fresh manifest"),
    "BeStar_Grism": (
        "~3-day-cadence Hα EW monitoring of bright Be stars, with detection "
        "claims only from the standards epoch (2025-12-05)",
        "Execute with amendments (U9). BeSS check first and alone; below "
        "two verified-active emitters the science paper stops; V/R held "
        "pending D1"),
    "DwarfGalaxy_AGN_Survey": (
        "Hα detections and non-detections for 13 dwarf-candidate fields "
        "plus NGC 5238 — vetting, not discovery",
        "Re-scope (U1). NGC 5548 broadband photometry is dead (slot '6' is "
        "a grism) and leaves the title; needs James to contact Cannon"),
    "Legacy_Rigel": (
        "Census of the pre-MACRO 2015–2022 archive: cameras first, then "
        "targets (directory `Legacy_Rigel/`; only ≈4% of it is the Rigel "
        "system)",
        "Census only, L0–L1 (U8). Go/no-go pre-registered; no science "
        "commitment"),
}

rows = ["| Project | Headline product | Venue posture | State after the "
        "2026-10-03 plan review |", "|---|---|---|---|"]
for project in pp.PROJECTS:
    headline, state = HEADLINE[project.key]
    rows.append(f"| **{project.key}** — {project.title} | {headline} | "
                f"{project.venue} | {state} |")
table = "\n".join(rows)

p.between("| Project | Headline product | Venue posture | State |",
          "Every panel's internal referee round traced its worst errors",
          table + '''

*Regenerated 2026-10-03 from the plan ledger's current titles and venues
(`pipeline/macro_core/project_plan.py`); the 2026-08-16 table it replaces sold
three polars and superhumps to ApJ, an "early" SN light curve with an ApJ upside,
and an NGC 5548 light curve — papers the committee has since ruled out (ED.E6).
Per-project task counts are deliberately not typed here: they are on each project
page and in `python pipeline/scripts/update_project_plan.py show`, which states
`done / in-scope` with the dropped count and the 2027 backlog beside it.*

### 0.1 Committee amendments 2026-10-03

**Binding.** Plan review of 2026-10-03 — `committee/reviews/2026-10-03/SYNTHESIS.md`
and the seven seat memos beside it. It amends this roadmap, the five strategies
(each now carries a §10 recording its own rulings) and the plan ledger. Where this
subsection and the rest of this document disagree, this subsection wins; the
sections below are kept as the 2026-08-16 record.

**What "complete" means (SYNTHESIS §0).** A plan task closes as **done** (criterion
met, evidence linked), **dropped** (ruled impossible or pointless on the data that
exist) or **deferred** (needs frames that do not exist yet → the project's *2027
backlog*, outside the paper's critical path and its completion count). A project is
complete when every in-scope task is done or dropped, its result has been through
the committee, and every seat records *satisfied*.

**Rulings that bind this roadmap.**

| # | Ruling | What it changes here |
|---|---|---|
| U1 | NGC 5548 broadband photometry is dead: slot '6' is a grism on every night | Settles §1.3 **C2** — the wheel did not change in May 2023; Dwarf Q4 is dead as photometry |
| U2 | Grism dispersion is hardware: one solution per (grism, mechanical epoch) from hot stars; per frame only a zero point | Replaces §1.1 **G**'s "per-frame zero point self-anchored … dispersion per filter/era" with a fixed dispersion per mechanical epoch |
| U3 | The CV paper is not finished — major revision | §0 table; Wave 2 no longer starts from a finished CV paper |
| U4 | Ledger blockers citing "S2 tables destroyed" are stale; blockers are computed from the database | The plan ledger; nothing in this document |
| U5 | The whole DAG reads STALE for a cosmetic reason, hiding real staleness | Convention 2 holds only at a clean-tree rebuild: no number is citable until F-10 |
| U6 | T CrB flickering and the period search leave paper 1 — no ≥2 h run is possible before mid-January 2027 | §2 Wave 0/1: the "weekly flickering runs at re-opening" are the 2027 backlog |
| U7 | SN 2023ixf is a +5.4 → +50 d validation/limits release (AJ/PASP); grism NOT PROMOTED; model fit dropped | §0 table; §1.2's "narrowband forward-modeling" survives only as the gated flash-phase colour |
| U8 | Legacy archive: census only (L0–L1), premise corrected | The "Rigel" name; a sixth project remains undecided |
| U9 | Be-star: BeSS novelty check first; pre-standards seasons are descriptive only | §2 Wave 2: no Be pipeline effort before the gate |
| U10 | The October request is rewritten (rev. 3) for the camera actually mounted (QHY600) and the sky actually available | §2 "Monsoon window": `ops/2026-08_observatory_request.md` is superseded |

**Disagreements resolved by test, not by vote.** **D1** — hrg dispersion 0.47 Å/px
or 1.59 Å/px (R ≈ 2,500 vs R ≈ 200): G-1's ≥3-line test on Vega and θ CrB; all
profile, velocity and V/R tasks are *held* until it returns. **D2** — Mode0 gain,
header 0.2467 or ≈1.0 e⁻/ADU: the flat-pair photon transfer (F-4). **D3** — the CV
band offset, astrophysical or estimator bias: per-band injection (CV-R2).
**D4** — a seventh, instrument/pipeline paper written first: **James's decision**;
until then the instrument material is built once as `docs/pipeline` evidence and
each science paper carries a condensed instrument section.

**Wave 0 is replaced by the shared foundation (SYNTHESIS §3).** F-1 QHY-era dedup ·
F-2 header re-scrape · F-3 `mech_epoch` table beneath the eras · F-4 flat-pair PTC
per camera → `detector_params` · F-5 one saturation/linearity policy · F-6 S2c
background-morphology test and a labelled truth set · F-7 data digests apart from
page digests · F-8 absolute clock from archived transits · F-9 evidence snapshot ·
F-10 clean-tree rebuild to all-FRESH · G-1 … G-5 the grism library. They are held
in the ledger as the *Shared foundation* group and gate project tasks by dependency.

**Conventions amended.**
- *Convention 1* (pixels come from `rawimage/`; `reduced/` is unaudited): 5,407 of
  the CV paper's 7,802 measured frames are server-reduced. Either the reduced tree
  is audited per era and this convention is amended, or those frames are re-reduced
  from raw (OA.E5, CV-R9). Undecided; it is not waived.
- *Convention 3* (global dedup): violated for every 2026 count until F-1 links the
  QHY-era `_calibrated` twins to their raw parents (TE.F6).
- *Convention 4* (era keys on READOUTM, geometry, EGAIN): necessary and not
  sufficient. The archive holds four cameras and at least eight mechanical states;
  a `mech_epoch` layer sits beneath the eras and bounds calibration validity
  (TE.F1, DE.F5, F-3).
- *§1.3 C1:* the "anchors on B/I" default is superseded — I is clipped at the target
  as well as R; T CrB's anchors are **B-only** (OA.E3). The one-afternoon hardware
  check is retired with the camera.
- *Standing statistical rules* (SYNTHESIS §5) are adopted portfolio-wide: χ²ν per
  stratum with dof and no `max(χ²ν, 1)`; every null with its recoverable effect size
  and predicted scale; signed matched-cell injection bias; saturation at the target
  in native pixels; a ≤ 250-word abstract approved before figures; a continuum light
  curve wherever EW is quoted.

**Sequencing as amended.** Shared foundation → T CrB (novelty table first) ∥ the CV
major revision → SN release → Be (only past its gate) → Dwarf (only past its
novelty gate). A novelty gate runs *before* pipeline work in T CrB, Be and Dwarf.
The legacy census runs at low priority and does not compete with T CrB.

**Only James can close (SYNTHESIS §6):** D4; the authorship policy, ORCIDs, Zenodo
DOI and an outside reader for the CV draft; contacting Cannon; sending the rev. 3
request to Winer and asking whether the ASI and AC4040 cameras still exist, for the
server `calibrations/` trees and a hardware change log; whether `manuscripts/` is
tracked.

''')
p.save()
print("roadmap ok")
