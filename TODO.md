# MACRO — working to-do list

Maintained by Claude as work proceeds; edit freely. Last updated **2026-10-04** — committee reviews BENCHED by James (token budget); work is now: refresh every stale stage and republish the site.
Goal: all six projects complete — every in-scope task done or dropped by committee ruling, and every
committee seat satisfied. Detail behind each line: [`committee/reviews/2026-10-03/SYNTHESIS.md`](committee/reviews/2026-10-03/SYNTHESIS.md)
and the reports in `committee/work/<package>/REPORT.md`. Task-level status lives in the plan ledger
(`python pipeline/scripts/update_project_plan.py show`).

## 1. Waiting on James

- [ ] **Read and send observatory request rev. 3** to Winer — time-sensitive (as-found flats lose value nightly). `ops/2026-10_observatory_request_rev3.md`, cover note in `committee/work/ops/emails/01_…`
- [ ] **Decide: seventh instrument/pipeline paper (PASP)?** Editor recommends it, with SN (and possibly T CrB) as validation sections.
- [ ] **Contact Cannon** — filter transmission curves (needed by 2026-11-15 for SN Hα), candidate dossier, and whether a dwarf paper exists given his AJ 170, 23. Draft: `committee/work/ops/emails/02_…`
- [ ] **Take `AUTHORSHIP.md` to the consortium**; collect real ORCIDs; reserve a Zenodo DOI.
- [ ] **Be-star sample**: I am proceeding on the re-drawn sample (19 BeSS-verified active stars + QQ Gem). Object if you would rather stop the paper.
- [ ] **CV novelty wording** — prior simultaneous multiband ST LMi curves exist (1985, 1992); abstract/intro claim must narrow.
- [ ] **Email VSX (S. Otero)**: what the ST LMi / AN UMa epochs mark and how the ST LMi period was derived.
- [ ] **Download from AAVSO by browser**: T CrB Rc/Ic/B/V photometry 2005–2026 (site blocks scripted access).
- [ ] **Version manuscripts?** `manuscripts/` and `products/` are gitignored — a private repo/branch would give drafts a history.
- [ ] Choose an outside polar-expert reader for the CV draft (send only after the revision lands).
- [ ] Eruption block: observer email/code, who is called on an alert, pre-agreed eruption author list.
- [ ] Confirm T CrB quiescent magnitudes assumed in rev. 3 (V 10.0; g/r/i 10.7/9.5/8.3); ask site about a B filter.

## 2. In progress (running now)

- [ ] **Unattended refresh of every stale pipeline stage** (`pipeline/scripts/refresh_stale.py`, log `products/refresh_stale.log`), gated on the tests of the code the stopped packages left half-edited. Then site rebuild + push.
- Stopped mid-work, resumable later from their transcripts (partial code edits are in the working tree, uncommitted): grism library, detector, legacy census, CV statistics, clock check, foundation (S0 rebuilt: 169,107 canonical frames; S0b/S0c being refreshed by the driver).
- [x] Legacy dedupe finished: ~202.5k verified duplicates removed (~2.1 TB freed). 830 frames remain uncompressed — 828 fail fpack (truncated downloads, need re-fetch), 461 failed twin verification.

## 3. Next, in order

- [ ] File the six outstanding package reports; apply their requested ledger changes
- [ ] **Integrated clean rebuild** of the whole pipeline at a tagged commit → all stages FRESH; products snapshot
- [ ] Full test suite green; site rebuilt (new figures onto figure walls; "Legacy Rigel" label; dropped tasks not shown as open)
- [ ] **CV** — integrate statistics + literature into the manuscript (abstract ≤ 250 words, ≤ 18 pp, AJ); rebuild numbers; macro diff
- [ ] **T CrB** — re-extract all frames on the corrected grism library; EW series; ARAS same-date cross-validation; continuum correction; pre-registered homogeneity test (the only AJ path)
- [ ] **SN 2023ixf** — re-run Gate 0; linearity audit; +5.4 → +50 d light curve vs published; flash-phase Hα gate
- [ ] **Be stars** — re-drawn sample; injection–recovery per star; EW series on corrected library
- [ ] **Dwarf / NGC 5548** — reduced scope: Hα detections/limits for the informative fields + NGC 5238; one-night broad-Hα triage for NGC 5548
- [ ] **Legacy** — apply the pre-registered go/no-go; write up outcome
- [ ] ~~Committee result review per project~~ — BENCHED until James says otherwise

## 4. Done

- [x] Standing committee: seven profiles (`.claude/agents/`), charter (`AGENTS.md`)
- [x] First plan review: seven memos + binding synthesis
- [x] Plan ledger carries every ruling (124 tasks amended; CV reopened with CV-R1…R15)
- [x] Observatory request rev. 3, eruption block draft, `AUTHORSHIP.md` draft, three email drafts (none sent)
- [x] Novelty gates: T CrB (not densest → PASP note unless homogeneity test passes), Dwarf (stop as scoped → reduced scope), Be (planned sample fails; archive holds 19 verified-active stars)
- [x] CV literature: 113 verified references, comparison section draft, citation map

## 5. 2027 backlog (needs frames not yet taken; outside completion count)

- T CrB flickering runs (first ≥ 2 h night 2027-01-09) and the QHY-era splice to the 2025 series
- ST LMi new season (≥ 1.5 h from 2026-10-30)
- Be-star season 2 and the grism dither test
- QHY calibration set (bias, darks, flat pairs incl. grisms), per-mode PTC and absolute clock night at re-opening
