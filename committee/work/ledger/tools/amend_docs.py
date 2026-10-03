"""Append the dated "Committee amendments 2026-10-03" section to each
strategy.  The rulings are written here by hand from the chair's synthesis;
the task table under each is EMITTED from the plan ledger by
`update_project_plan.py amendments <Project>`, so the ids the ledger cites
are the ids the section holds."""
import subprocess
import sys
from pathlib import Path

REPO = Path("/Volumes/OWC StudioStack HDD/Dropbox/01_Research/MACRO")
PY = "/opt/miniconda3/envs/rlmt-checks/bin/python"
HEADING = "## 10. Committee amendments 2026-10-03"

PREAMBLE = """**Binding.** Plan review of 2026-10-03: the chair's synthesis
`committee/reviews/2026-10-03/SYNTHESIS.md` and the seven seat memos beside it.
**Where this section and anything above it disagree, this section wins**; the
text above is kept as the record of what was planned on 2026-08-16. Finding ids
are `<seat>.<id>` — DS data scientist, OA observational astronomer, PH physicist,
DE detector engineer, TE telescope engineer, ED editor, RF referee; `U#` and `D#`
are the synthesis's own rulings (§1) and test-resolved disagreements (§2).

A task now closes in exactly one of three ways (SYNTHESIS §0): **done**
(acceptance criterion met, evidence linked), **dropped** (ruled impossible or
pointless on the data that exist) or **deferred** (needs frames that do not exist
yet — the *2027 backlog*, outside this paper's critical path and its completion
count). The plan ledger (`pipeline/macro_core/project_plan.py`) carries each
ruling as data; this section is its citation."""

STANDING = """**Standing statistical rules (SYNTHESIS §5) — apply to everything this project
produces.** (1) Every χ²ν is reported per band/era/mode with its dof; χ²ν < 0.5 is
a defect equal to χ²ν > 2; no `max(χ²ν, 1)`. (2) Every null carries the effect size
it would have recovered and a "predicted scale" beside it; "k of N significant" is
banned. (3) Injection grids report the signed matched-cell bias, never a median of
absolute values. (4) Saturation is judged at the target, per frame, in native
pixels. (5) A ≤ 250-word abstract with number placeholders is approved by seat 6
before figures are built. (6) EW is a ratio: wherever the continuum varies, carry a
continuum light curve.

**Shared foundation.** This project stands on the Wave-0 tasks of SYNTHESIS §3
(F-1 … F-10, G-1 … G-5), held in the ledger as the *Shared foundation* group. Where
a task below waits on one of them, the ledger states it as a dependency and reads
its status from the database — no blocker here is prose about a table."""

BODY = {
"TCrB_Monitoring": """**Verdict:** execute with amendments — the highest-priority paper in the
portfolio, and the only one with a clock on it.

**Venue:** ApJ only if the novelty table (TCRB-N1) shows a density or homogeneity
advantage over ARAS; otherwise AJ or a PASP validation note. TCRB-N1 runs first.

**Scope as amended.**
- **The paper is an Hα equivalent-width series, Feb–Jun 2025, with ARAS
  cross-validation.** Counts are reported by S2c verdict on unique frames, not by
  filter name (DS; RF: 402 unique frames, not 471).
- **Grism dispersion is hardware (U2).** One solution per (grism, mechanical epoch)
  from hot stars; per frame only a zero point. The per-frame dispersions on disk
  (−1.8 to +2.0 Å/px) are unphysical and are retired. This replaces §4 Phase A
  step 3's "per-frame self-anchored" solution.
- **D1 — hrg is ≈0.47 Å/px (OA) or ≈1.59 Å/px (code; PH and TE reasoned from it).**
  Decided by test, not by vote (G-1: Hα plus ≥2 of telluric 6277/6867/7186/7594 and
  Hβ on Vega and θ CrB; ≥3 lines, residual < 1 px). **Profile morphology and wing
  velocities are HELD pending D1 — neither dropped nor kept.**
- **D2 — Mode0 gain.** The grism variance and the saturation threshold are read from
  the measured `detector_params` (F-4, G-2); the 16.3 kADU "rail" is a hot-pixel
  signature and is replaced by a hot-pixel mask and the true clip.
- **The identity gate is pixel-based (G-3).** The header is never the sole reason a
  frame is rejected.
- **§4 Phase A step 4 is dropped as written.** No absolute flux from a θ CrB
  response. It is replaced by TCRB-A5b: F(Hα) = EW × continuum from contemporaneous
  AAVSO or zero-order photometry, with EW and flux both shown against orbital phase
  (the giant's ellipsoidal continuum, ±10%, modulates any raw EW).
- **The EW-change detection rule is pre-registered before A5 runs** (TCRB-A5a), and
  ARAS spectra degraded to the measured LSF must reproduce their native EW.
- **Phase B anchors are B-only** (plus singletons): R and I are clipped at the
  target. A scripted peak-at-target census of every imaging frame decides it. Chair's
  ruling 2 ("anchors on B and I") is superseded; the nightly-mean precision promise
  of Phase B step 6 is deleted.
- **The five tasks blocked on "S2's tables were destroyed" are unblocked (U4).** The
  tables exist. P0-bitdepth is closed from the archive.
- **The mechanical epoch of every frame, the measured gain, and a temperature split
  of the grism frames are Phase-0 gates** (TCRB-P0-mech-epoch, -gain-ptc,
  -temp-split).
- **An eruption contingency exists as artefacts**: the pyscope block committed with
  a dry-run log (TCRB-P0-eruption-block) and an RNAAS/ATel skeleton (TCRB-N2).

**Dropped:** A4 (θ CrB response → absolute flux); C3 (per-season period search).
C1 shrinks to one sentence and one table.

**Deferred to the 2027 backlog (U6):** C2 (weekly ≥2 h flickering runs — none is
possible before mid-January 2027) and the October / January-2027 restart, which is
a *second instrument* (QHY600) and may be spliced to 2025 only through same-night
θ CrB + Vega on ≥3 nights. Flickering is paper 2.

**Figure cap:** six (was ten), each mapped to a sentence of the abstract.""",

"CV_TimeSeries": """**Verdict:** MAJOR REVISION (U3). The plan read 34 of 34; the paper is not
finished. CV_TimeSeries is **reopened**, with fifteen tasks (CV-R1 … CV-R15), all of
which block release. Every number the referee recomputed reproduces; four of the
inferences drawn from them do not survive.

**Venue:** AJ. The header's "Target journal: ApJ" is superseded.

**Scope as amended.**
- **The band-offset "null" is contradicted by the paper's own edges** (g−i ≈
  −110 ± 28 s; 9/9 nights or 11/12 cycles negative). The published ±51 s exists only
  because of `max(χ²ν, 1)` on a transported budget. Re-tested with paired,
  scatter-based statistics (CV-R1).
- **D3 — astrophysical or estimator bias? Not decided by the committee.** Decided by
  per-band injection with per-band ramp widths and the signed matched-cell bias
  carried (CV-R2). The paper reports whichever the test supports; the sentence "no
  offset" goes either way.
- **The O−C is refitted** with per-band constants, night-level epochs (N = 17) and
  an era-offset nuisance term; Ṗ is quoted each way and with 2024 dropped (CV-R3).
- **The superhump "measurement of absence" is logically inverted** and is restated
  as an excluded-amplitude range (CV-R4).
- **The Ṗ bound is a spin/spot-longitude bound**, not an orbital one, and is given
  three physical scales and a longitude stability in degrees by accretion state
  (CV-R5).
- **The colour–phase result leads.** It is quantified (amplitude, phase, two-era
  repeatability) and must survive a ≤ 120 s pairing window (CV-R6).
- **The edge fit is disclosed** (CV-R7); **the counts mean what they say** (CV-R8);
  **the reduction is described and tested** (CV-R9); **the instrument section uses
  measured gains** (CV-R10); **the clock residual is printed and the clock tested on
  transits** (CV-R11); **the paper meets the literature** (CV-R12).
- **Restructure:** abstract ≤ 250 words leading with the colour result; ≤ 18 pages;
  process prose to an appendix; internal language removed (CV-R13).
- **Nothing is citable until the chain is rebuilt FRESH at a clean commit** and a
  macro-by-macro diff of `numbers.tex` is filed (CV-R14, U5).
- **§9's "request one more g/r/i ST LMi season" is not a reason to wait.** Submit
  what exists (ED).

**Needs James (CV-R15):** authorship policy with the consortium, real ORCIDs, a
Zenodo DOI, availability statements, and an outside reader who works on polars.

**Not changed:** nulls reported as nulls, the tie treatment, the saturation vetoes,
the refusal to publish per-cycle error bars, the AN UMa grading and the 307/307
macro traceability were each recorded as *satisfied* by the seats that checked them.""",

"SN2023ixf_LightCurve": """**Verdict:** RE-SCOPE (U7).

**Venue:** AJ/PASP — decided, by the rule this strategy pre-registered (§2): ApJ only
if Gate 0 promoted the grism or recovered a narrowband bandpass. It did neither. The
"ApJ upside decided at week 3" is closed. Whether this release becomes the
validation section of an instrument paper is D4, which is James's decision; until
then it proceeds as its own short release.

**Scope as amended.**
- **A +5.4 → +50 d validation and limits release.** "Early" leaves the title. The
  first clean broadband night is +5.4 d; if that is at peak (PH: published optical
  maximum is nearer +5–7 d than this document's "+8–9 d"), the paper contains no
  rise and says so. The peak epoch is read from the literature.
- **The grism series is NOT PROMOTED** (SN-G0c closed): zero extracted spectra, no
  wavelength source, no contamination test. Its frames go in the release.
- **Gate 0 is re-run on a fresh manifest before anything cites it**: the 438/439
  arithmetic is fixed, the one unnamed exclusion rule is named, and all 632
  broadband frames are S2c-measured.
- **SN-S2-linearity is the key detector task of the portfolio**: add the
  twilight-flat ramp, the 2023-06-07 Albireo set, a gain per EGAIN epoch, and
  StackPro on its own curve.
- **Step 4:** the χ²ν ≈ 1 claim is made on held-out check stars. **Step 5:** the
  error model carries an explicit scintillation term. **Step 8:** nightly means and
  a whole-night bootstrap; **the intra-night periodogram is dropped.**
- **Step 6 is split.** *S6a, flash phase:* the (H − "[S II]") differential colour,
  with filter widths (~65 Å) from the zero-point ratios, behind a
  **predicted-excess gate**. *S6b, ejecta phase:* dropped unless a transmission
  curve arrives by **2026-11-15**.
- **A residuals-vs-published table precedes any variability limit.**
- **A template table** lists camera, mechanical epoch, FWHM and focus offset per
  template epoch.

**Dropped:** Step 7 (model consistency fits) — no clean data before +5.4 d and no
UV; one overlay on the published curves stands in for it.

**Figure cap:** five (was twelve).""",

"BeStar_Grism": """**Verdict:** execute with amendments (U9) — **BE-S-1a-bess first and alone.**

**Venue:** decided by the novelty gate, BE-N1-gate: ApJ if ≥4 BeSS-verified active
emitters, AJ/PASP for 2–3; **below two, the science paper stops** and the standards
and precision material moves to the instrument section.

**Scope as amended.**
- **Step −1(a) runs before any pipeline effort**, and BE-N1-gate records the
  decision. BE-figures and BE-draft do not start until it returns ≥2 verified
  emitters.
- **Step −1(b) is extended:** injection–recovery for every star, through
  detrending, with the contour published *before* any periodogram is opened.
- **The grism library G-1 … G-5 is inherited** (U2): fixed dispersion per (grism,
  mechanical epoch), measured gain and saturation, pixel-based identity, the
  sky-lozenge background, focus and CCD-TEMP regressors with a measured LSF. Step 4
  must demonstrate a stable dispersion per (grism, epoch) before any EW.
- **§3.4's era table is replaced:** the instrument has five mechanical states
  (Andor; ASI pre-monsoon; ASI post-monsoon flipped; QHY night 1; QHY, with the
  grisms re-seated at the swap), and its gains are measured, not read from headers.
- **V/R is HELD pending D1** — neither dropped nor kept. At ≈0.47 Å/px it is in
  scope; at ≈1.59 Å/px Step 5's own FWHM < 4 Å rule excludes it.
- **Step 6** carries a telluric H₂O regressor (7200 Å band depth). **Step 8:**
  season 1 has no standard, so season 1 is EW only. **Step 10:** a detection rule
  exists only from the standards epoch (2025-12-05); earlier seasons are
  descriptive, with no ΔEW claims.
- **Step 11 short tier:** the three qualifying stars, with a global false-alarm
  probability; **the λ Eri short search is dropped.**
- **A disposition column** for the 308 indeterminate and 24 direct staged frames;
  **QQ Gem is dispositioned.**
- **Step 0's header re-scrape** takes GAIN, OFFSET, SET-TEMP, COOLPOWR, CCD-TEMP,
  FOCPOS, FLIPSTAT, TELPIER, FWPOS and FWALLNAM (F-2).

**Deferred to the 2027 backlog:** the dither test in lieu of grism flats, and the
season-2 observing items (carried in the rev. 3 observatory request).

**Figure cap:** six (was twelve).""",

"DwarfGalaxy_AGN_Survey": """**Verdict:** RE-SCOPE (U1).

**Venue:** AJ at best — detection/non-detection with broadband depth. Calibrated Hα
fluxes only if Cannon's filter curve arrives. The "ApJ (NGC 5548 section
conditional)" posture is superseded.

**Scope as amended.**
- **NGC 5548 broadband photometry is dead (U1, unanimous).** Filter slot '6' is a
  grism on every night. NGC 5548 leaves the title. Phase 4.1, 4.2–4.3, 4.4 and 4.5
  are dropped; Phase 0.3's plate-solve no longer includes those frames (spectra do
  not plate-solve — that blocker could never clear). DW-P02-slot6 closes with a
  per-night verdict table from F-6's background-morphology test.
- **Two triage tasks stand in for Phase 4**, each of which earns a paragraph or is
  dropped: DW-P4x, a one-night broad-Hα extraction (two days; EW and nightly
  scatter; **no lags**), and DW-P4y, zero-order differential photometry feasibility.
- **Novelty before stacking:** DW-N1-novelty (the Phase 3.1 literature cross-match,
  as a gate) and DW-P36-0 (a depth-versus-expected-L(Hα) table) come before any
  coadd.
- **Phase 0.4 is unblocked** with the measured StackPro read noise (≈16.5 e⁻ per
  pixel per frame); the sum-versus-average contradiction about StackPro is
  resolved; depth promises are re-stated at the read-noise floor.
- **Phase 1.1:** no AC4040 flat will ever exist — the camera left in March 2024. The
  flat is a night-sky superflat from the June 2023 L frames.
- **Hα products:** a detection through a ~65 Å filter is a velocity statement;
  a non-detection carries the out-of-band caveat and is not an SFR limit.
- **W leaves NGC 5238's surface photometry.** The completeness map (Phase 5.4) is
  kept for NGC 5238 only.

**Dropped:** Phase 4.1, 4.2–4.3, 4.4, 4.5 (NGC 5548); Phase 5.3 (period search) and
5.6 (eclipse timing); Phase 5.2 (transient search) unless its threshold and trials
are pre-declared with injection.

**Needs James:** contact John Cannon — filter transmission curves, the candidate
dossier, and authorship. Four tasks and the lead result wait on a conversation that
no observatory closure ever blocked.

**Figure cap:** six (was fourteen).""",
}

BANNER = ("> **Amended 2026-10-03.** The plan review of 2026-10-03 binds this "
          "strategy; its rulings are in **§10 (Committee amendments "
          "2026-10-03)** at the end of this document and override anything "
          "above them that they contradict.\n\n")


def table(project: str) -> str:
    out = subprocess.run(
        [PY, str(REPO / "pipeline/scripts/update_project_plan.py"),
         "amendments", project], capture_output=True, text=True, check=True)
    return out.stdout.strip()


for project, body in BODY.items():
    path = REPO / project / "ANALYSIS_STRATEGY.md"
    text = path.read_text(encoding="utf-8")
    if HEADING in text:
        text = text[:text.index(HEADING)].rstrip()
        if text.endswith("---"):
            text = text[:-3].rstrip()
    lines = text.split("\n")
    if not any(l.startswith("> **Amended 2026-10-03.**") for l in lines[:6]):
        lines.insert(2, BANNER.rstrip("\n"))
        lines.insert(3, "")
        text = "\n".join(lines)
    section = "\n\n".join([
        HEADING, PREAMBLE, body, STANDING,
        "### Task amendments\n\nEmitted from the plan ledger by "
        f"`python pipeline/scripts/update_project_plan.py amendments {project}` "
        "— do not edit by hand; change the ledger and re-emit. *Action* is the "
        "committee's verb; *Ruling* is the place in the synthesis and the finding "
        "ids behind it.\n\n" + table(project)])
    path.write_text(text.rstrip() + "\n\n---\n\n" + section + "\n",
                    encoding="utf-8")
    print("amended", path.relative_to(REPO))
