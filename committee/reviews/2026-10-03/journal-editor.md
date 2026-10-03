# Seat 6 — Journal Editor · plan review of 2026-10-03

Reviewed at commit `27c962f`. Read: ROADMAP, the five strategies, the plan ledger, `check_pipeline_status.py`, the
manifest database (read-only), all `manuscripts/*/main.tex`, the CV `numbers.tex`, `references.bib`, `main.log`, Figs. 6
and 9, the Rigel case page and the ops request. I opened no FITS frames; pixel questions belong to seats 2 and 4.

---

## 1. Existing work

### E1 — BLOCKER (for submission; not for the science) · The CV paper has no literature
`manuscripts/CV_TimeSeries/references.bib` holds **10 entries**; seven are methods/software. The science bibliography is
Cropper 1990, Ferrario 1993, Honeycutt 1992. Each is cited once. The text nonetheless:
- uses ZTF, ASAS-SN, AAVSO, VSX and ATLAS-REFCAT2 as data and as the calibration reference (abstract, `main.tex:72-139`;
  acknowledgments; `\facilities`, `:1565`) with no citation to any of them;
- asserts "published ZTF polar folds do" (`main.tex:172`) with no citation;
- adopts catalogue ephemerides for all five targets (`numbers.tex`, "VSX catalogue, external constant") without citing the
  original ephemeris papers the 1,071 s edge offset and the 21,869-cycle count are measured against;
- reports a YZ Cnc superhump null and "published superhump amplitudes" (Conclusions) without citing one superhump paper;
- bounds |Ṗ| without citing one prior ST LMi period or period-change measurement.

The strategy promised exactly this comparison (`CV_TimeSeries/ANALYSIS_STRATEGY.md:176`: Ferrario/Cropper, WD-binary timing
literature, Kato superhump surveys, Duffy et al. states). The Discussion (`main.tex:1156-1278`) contains none of it.
A paper whose claim is "what these data add" cannot go to referees when it never says what existed. Six internal
referee rounds did not catch this; they were all inward-facing.
**Closes when:** every external dataset, catalogue, ephemeris and comparison result has a citation; the Discussion has a
"comparison with previous work" subsection for ST LMi (colour-phase morphology, bright-phase longitude, period) and YZ Cnc;
bibliography ≥ ~40 entries and `main.log` still clean. **Would change my mind:** nothing — this one is mechanical.

### E2 — MAJOR · The abstract is 551 words and leads with the wrong result
`main.tex:72-139`, counted with `wc -w`. AAS journals cap at 250. More important than the length: the abstract spends
~300 words on a period-change null, band-offset upper limits and alias caveats, and **zero numbers** on the one positive,
visually unmistakable result the paper owns — `figures/fig06_colourphase.pdf`: ST LMi's g−r swings by ~0.75 mag through
the bright phase, with the same morphology in two instrument eras a year apart (212 and 557 quasi-simultaneous pairs),
state-tagged against the survey record. The text itself calls this "the headline physics figure" (`main.tex:825`). The
abstract calls the null "the strongest timing statement", which is true and is not the same as the paper's result.

The |Ṗ| < 3.6×10⁻⁹ bound (timescale > 6×10⁴ yr on a 1.9 yr baseline) is honest and physically uninformative; it
belongs in one sentence.
**Closes when:** abstract ≤ 250 words, structured as: what was observed → the colour-phase result with amplitudes and the
two-era repeatability → the timing null in one sentence with its bound → YZ Cnc in one sentence → the release. Seat 3 should
confirm the colour amplitude is quotable given the tie systematic (28–38 mmag clipped is small against 0.75 mag; I believe
it is). **Would change my mind:** seat 3 or 1 showing the colour swing is a calibration artefact.

### E3 — MAJOR · Length and structure are disproportionate to the result
32 pages (`main.log:516`), 13 figures, 245+ macros, for one target's colour curve, one null and a coverage audit of three
targets that support nothing. Process prose has leaked into the science: "A note on how this paper is written"
(`main.tex:196`), the release section's count of macros per database (`:1280-1300`), "the strongest reproducibility claim we know how to make" (§7), "first
in the queue for the next observing period" (`:1270-1271`). The integrity machinery is an appendix and a data-availability statement, not Section 7.
**Change:** target ~14–16 pages. VV Pup / EU UMa / AN UMa collapse to one table plus one figure ("coverage audit");
AN UMa's 4-of-N capability grading moves to an appendix or the release. The provenance narrative moves to an appendix.
**Closes when:** page count ≤ 18 and every main-text figure is cited in support of a claim in the abstract.

### E4 — MAJOR · Authorship, ORCID, contributions, data DOI
`main.tex:67`: a single author with ORCID `0000-0000-0000-0000`, for frames taken by a five-college consortium over three
years; the same placeholder is in all four skeletons. No author-contribution statement, no named observers, no data DOI
(grep for zenodo/doi: none), no licence. I find **no authorship policy anywhere in the repo**.
**Closes when:** a one-page `AUTHORSHIP.md` (criteria, order rule, student observers, opt-in window) exists and is agreed
by the consortium; real ORCIDs; a Zenodo DOI reserved for the three databases; AAS data-availability and
software-availability statements in place.

### E5 — MAJOR · The manuscript's numbers are built on a chain that now reads STALE
`check_pipeline_status.py`: every CV stage is `STALE_UPSTREAM`, and **CV-S5 is STALE on its own input**
(`cv_selection: rows 34 -> 33`). `numbers.tex` was built 2026-08-20 at commit `1380616`. 
**Closes when:** the CV chain is re-run to FRESH and a macro-by-macro diff of `numbers.tex` is filed with the response.
"No macro moved" is an acceptable and likely outcome; it has to be shown.

### E6 — MAJOR · Skeleton titles and the ROADMAP advertise papers the strategies already killed
- `manuscripts/TCrB_Monitoring/main.tex:4`: "…A Three-Year Quiescent Baseline" — the strategy's §1 says the data cannot
  support that paper and AAVSO owns it.
- `manuscripts/SN2023ixf_LightCurve/main.tex:4`: "Early Multi-Band Photometry" — the strategy calls this a desk reject;
  the first clean epoch is +5.4 d (`sn_g0_verdict`).
- `manuscripts/DwarfGalaxy_AGN_Survey/main.tex:4` promises "Photometric Monitoring of NGC 5548" (see E7).
- `ROADMAP.md:15` still sells CV as "3 polars + YZ Cnc superhumps … ApJ"; the delivered paper is one polar and a null.
- `CV_TimeSeries/ANALYSIS_STRATEGY.md:4` still says target journal ApJ.
**Closes when:** titles and the ROADMAP §0 table are regenerated from current verdicts. MINOR individually; MAJOR together
because the site now presents these to outsiders.

### E7 — MAJOR · NGC 5548: the plan still carries a photometry section the manifest says may not exist
`frame_dispersion` for NGC 5548 filter `6`: **91 dispersed / 17 direct / 35 indeterminate** of 143 frames; only 3 nights
have any direct frames. ROADMAP C2 warned of exactly this. Yet the ledger blocks `DW-P03-plate-solve` on "the S1b batch
reaches these frames" (`project_plan.py:1314ff`) and leaves `DW-P41…P46` waiting on it. The dispersed verdicts are all
`strength_class='ambiguous'` and source counts swing from 11 to 19,263 per frame, so the classification is not settled —
seat 2 must look at the frames. Either way: no title, abstract or site page may promise NGC 5548 monitoring until
`DW-P02-slot6-dispersion` closes with pixels.
**Would change my mind:** ≥ 10 nights of visually confirmed direct images.

### E8 — MINOR · Ledger blocker text contradicts the database
`TCRB-P0-bitdepth`, `TCRB-P0-ladders`, `DW-P04-noise-model` are blocked because "S2's tables were destroyed"
(`project_plan.py:675, 684, 1329`). The tables exist: `s2_ptc_fits` 4 rows, `s2_ceiling_modes` 8, `s2_linearity_ladders` 66,
`s2_build_meta.built_utc = 2026-08-20`. The blockers are stale; the true blocker is "S2 reads STALE_UPSTREAM".
**Closes when:** blocker strings are re-derived from the DAG.

### E9 — NOTE · What is good
Nulls are reported as nulls and limits as limits throughout the CV draft; the SN Gate 0 venue rule was pre-registered and
then obeyed (`sn_g0_verdict`: AJ/PASP, grism NOT PROMOTED). I am **satisfied** with the CV draft's claims language.

---

## 2. Plan hardening, per project

### CV_TimeSeries (34/34)
- **ADD CV-R1-literature** — closes E1. *Accept:* comparison subsection present; every external resource cited.
- **ADD CV-R2-abstract-restructure** — E2/E3. *Accept:* abstract ≤ 250 words leading with the colour-phase result; ≤ 18 pp.
- **ADD CV-R3-rebuild-on-fresh** — E5. *Accept:* chain FRESH; `numbers.tex` diff filed.
- **ADD CV-R4-authorship-release** — E4. *Accept:* author list agreed, ORCIDs, Zenodo DOI.
- **ADD CV-R5-external-read** — one reader outside the consortium who works on polars. *Accept:* their one-sentence summary
  of the result matches ours.
- **CHANGE venue:** AJ, not ApJ. A well-characterised single-target multi-colour dataset with a null and a release is an AJ
  paper. The YZ Cnc material (superhump null, flickering structure function) is separable; I would keep it as one short
  section rather than split to RNAAS, but it must earn its place in ≤ 1.5 pages.
- **Pointless:** waiting for "one more g/r/i season". Submit what exists.

### TCrB_Monitoring (2/34)
The only paper with a clock on it, and the strongest outside-reader case — *if* the novelty claim is proved first.
- **ADD TCRB-N1-novelty-table (do first, before A2–A8):** spectra per month, RLMT vs ARAS vs Asiago/Munari, Feb–Jun 2025.
  The strategy demands it; no task owns it. *Accept:* table on disk; if ARAS alone has denser, higher-resolution Hα
  coverage, the paper re-scopes to a validation/methods note (PASP) and says so.
- **ADD TCRB-N2-eruption-contingency (editorial):** pre-written RNAAS/ATel skeleton with the 2025 EW series as baseline.
  *Accept:* skeleton compiles with script-emitted numbers.
- **CHANGE TCRB-D4:** "full ten-figure set" → six. Centrepiece (EW vs time over AAVSO B), spectra montage, calibrator
  floor, cross-validation, flickering limits table as a table. *Accept:* each figure maps to an abstract sentence.
- **DROP from this paper:** `TCRB-C2-2026-runs` and `TCRB-C3-period-search`. C2's data do not exist and "accrue mainly in
  the 2027 season" by the ledger's own words; waiting forfeits the timeliness that is the paper's value. C3 on 23 d and 42 d clumps "would only advertise emptiness" (strategy ruling 5). Flickering
  becomes paper 2.
- **CHANGE Phase B** to the minimum that places anchors on the AAVSO curve; plan R-band loss as the default (ROADMAP C1).
- Honest scope: an EW time series, absolute flux only Mar–Apr 2025, from the 70 gate-accepted sample extractions extended
  to the 247. Target ApJ only if N1 shows a real density or homogeneity advantage; otherwise AJ/RNAAS.

### SN2023ixf_LightCurve (5/20)
- **CLOSE SN-venue-decision now.** `sn_g0_verdict` already answers it: AJ/PASP. Leaving it `in_progress` invites drift.
- **DROP SN-S7-model-consistency.** No clean data before +5.4 d, no UV; a consistency fit adds a referee target and no
  information. One overlay on published curves suffices.
- **DROP or demote SN-S6-halpha-curve** unless the transmission curve arrives by a hard date (I suggest 2026-11-15).
  Without it the ledger itself calls the product a "methods demonstration". *Accept:* dated go/no-go recorded.
- **CHANGE SN-G0c:** the timebox has lapsed with 0 extracted spectra, 0 contamination tests, no wavelength source. Either
  one more bounded attempt on the 3 flash-window nights, or close as NOT PROMOTED and publish the frames in the release.
- **CHANGE SN-figures:** twelve → five.
- **Honest reduced scope:** 438 clean broadband frames, +5.4 → +50 d, tied to PS1/REFCAT2, residuals against the published
  world dataset, night-to-night variability limits, saturated-frame inventory, release. That is a short **PASP or RNAAS + Zenodo**
  product — the photometric pipeline's external validation — or the validation section of the instrument paper (§3.1).
- **Re-run Gate 0** before anything else (brief: not re-run since reclassification).

### BeStar_Grism (4/24)
- **Do BE-S-1a-bess first and alone.** A database check costing days decides the venue and whether there is a science
  paper at all; it has been "pending" since August while pipeline tasks progressed. *Accept:* per-target
  active/inactive table with BeSS spectra dates bracketing our seasons.
- **ADD BE-N1-gate:** if < 2 verified-active emitters, the Be paper stops as a science paper and its standards/precision
  material moves into the instrument paper. *Accept:* decision recorded.
- **DROP BE-S11 short-tier search** on λ Eri (the strategy already expects the injection test to fail it); keep
  Phecda/Spica/φ Leo only as upper limits, in a table.
- **CHANGE BE-S8-flux-tiers:** season 1 has no standard; publish EW only for season 1.
- **CHANGE BE-figures:** twelve → six.
- **ADD:** disposition QQ Gem (ROADMAP C3 still open).

### DwarfGalaxy_AGN_Survey (1/31)
- **Contact Cannon this week.** Four tasks and the Hα lead result wait on a conversation the ledger says has not happened
  and that the monsoon never blocked. It is also an authorship question (E4). *Accept:* dated reply on file.
- **CHANGE DW-P41–P46 (NGC 5548):** gate behind `DW-P02-slot6-dispersion`; default outcome is DROP, with any surviving
  direct nights reported as an RNAAS at most. Remove NGC 5548 from the title now.
- **DROP DW-P53 period search, DW-P56 eclipse timing, DW-P52 transient search** from this paper. A variable-star census in
  dwarf fields is a different paper for a different reader.
- **Hα fluxes are impossible without the filter curve** (ledger says so). Without it the honest product is detection /
  non-detection morphology and broadband depth — AJ at best.
- **ADD DW-N1-novelty:** Phase 3.1 literature cross-match *before* any stacking. *Accept:* per-candidate table of existing
  Hα/spectroscopic status; if most are already classified, stop.
- **CHANGE DW-figures:** fourteen → six.

### Legacy_Rigel (0/11)
- The census plan is the right plan. **ADD RIG-N1:** the go/no-go must name a science question, an outside reader and a
  venue, or the outcome is "archive data release only" (a Zenodo record plus a short data note). Contact-binary O−C
  extensions are natural JAAVSO-class student papers; say that in advance.
- Do not let L0–L1 compete with T CrB for effort this autumn.

---

## 3. Cross-cutting — top five

1. **Write the instrument/data paper, and write it first (PASP).** Every strategy defers the same material to an appendix:
   High Gain ceiling and StackPro noise, filter-slot identity, timing audit, astrometric solvability, ensemble-photometry
   validation, grism identity gate. The CV draft is 32 pages largely because it carries all of that itself. One
   "RLMT: instrument, archive and pipeline" paper that the science papers cite would shorten each of them, is the
   one paper only this consortium can write, and is the natural home for student co-authors. The SN curve is its
   validation section. Six papers become: instrument (PASP), CV (AJ), T CrB (ApJ/AJ), Be (conditional), Dwarf
   (conditional), Rigel (undecided).
2. **Authorship and credit policy before anything circulates** (E4). One document, all papers.
3. **A novelty gate per project, executed before pipeline work**: T CrB N1, Be S-1a, Dwarf N1. Each is days of reading and
   decides venue. None is currently first in its "NEXT UP" list except Be's.
4. **Regenerate the public face from verdicts**: ROADMAP §0, skeleton titles, strategy venue lines, blocker strings (E6,
   E8). 
5. **Abstract-first discipline.** For each paper, a ≤ 250-word abstract with [NUMBER] placeholders is written and approved
   by this seat before figures are built, and figure counts are capped at what the abstract needs. The CV draft shows what
   happens otherwise.

---

## 4. Verdict

| Project | Verdict |
|---|---|
| CV_TimeSeries | **Execute with amendments.** Not releasable: E1 (blocker), E2–E5. Science and claims language: satisfied. Venue AJ. |
| TCrB_Monitoring | **Execute with amendments** — novelty table first; drop 2026 flickering and period search from paper 1. Highest priority. |
| SN2023ixf_LightCurve | **Re-scope** to a short validation/data paper (PASP/RNAAS) or a section of the instrument paper. Close the venue task. |
| BeStar_Grism | **Execute Step −1 only**, then decide. Not ready for pipeline effort. |
| DwarfGalaxy_AGN_Survey | **Re-scope.** Drop NGC 5548 from the headline and the variability census; nothing proceeds until Cannon is contacted. |
| Legacy_Rigel | **Execute the census** at low priority; go/no-go must name a reader and a venue. |

**Deferred to other seats:** whether NGC 5548's slot-6 frames are spectra (2); whether the ST LMi colour amplitude
survives the tie systematic (1, 3); whether the High Gain ceiling kills T CrB R (4).
