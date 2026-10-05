# Legacy archive — go/no-go criteria, pre-registered

**Written 2026-10-03, before the census existed.** No header scan of the legacy tree had been run, no
`legacy_manifest*.csv` had been opened, and `products/legacy/` did not exist when this file was saved.
Its sha256 and modification time are recorded in `committee/work/legacy/REPORT.md` and stored in
`products/legacy/legacy.sqlite` (`census_meta.prereg_sha256`) by the census build, so any later edit
of this file is detectable. The census script applies §4 mechanically; nobody chooses the outcome.

Authority: chair's synthesis of 2026-10-03, ruling U8 and §4 "Legacy archive"; findings DS
(RIG-L1-selection, RIG-L2-gonogo), ED (RIG-N1), PH (contact-binary criterion), OA, DE, TE.F10, RF.

## 1. What was already known when this was written (disclosure)

Pre-registration is only as honest as its statement of prior knowledge. Before writing this I had read:

- the seven committee memos and the synthesis, which state: ≈210k files, ~1,100 nights, 2015–2022;
  ≈7.8k files (≈4%) from the Rigel/PL16803 through 2015 day ≈120, the rest from the 0.5 m with six
  cameras (TE.F10);
- the ledger's "first scan" note: a programme "dominated by contact binaries and eclipsing systems
  (W UMa, XY Leo, TU Boo, RW Com, CC Com, RZ Com, AW Vir, TX Cnc), plus HAT-P-12, a comet and an asteroid";
- the top-level directory names of the archive (`2015` … `2022`, `archivar`, `old`).

I did **not** know, for any target: the number of nights, the filters, the run lengths, the seasons
covered, whether calibration frames exist, or whether any RLMT-era target is present at all. The
thresholds below were therefore set without sight of the quantities they test.

## 2. Definitions (fixed here so the census cannot tune them)

- **Canonical frame (dedup rule).** Only `*.fz` files are scanned. A frame's identity is its header key
  (camera identity, `DATE-OBS` to full recorded precision, exposure time, filter, image geometry).
  Files sharing a key are copies of one exposure; exactly one is canonical, chosen as: a raw file over a
  reduced/calibrated twin, then the lexicographically smallest path. Every threshold below counts
  canonical science frames only. Files that cannot be read, have no usable `DATE-OBS`, or are named in
  `legacy_manifest_BAD_collisions.csv` without a surviving file are *named exclusions*, listed by
  reason; the census is valid only if files-on-disk = rows + named exclusions.
- **Night.** Local noon to local noon at Winer (UTC−7): night = calendar date of (`DATE-OBS` − 19 h).
- **Target.** The alias-merged object name (case, spaces, underscores, hyphens folded; documented alias
  table). Calibration and pointing/focus frames are not targets.
- **Series.** One target × one filter × one camera. Frames from different cameras are never pooled.
- **Run.** Within one night and one series, a maximal sequence of frames with no gap longer than
  30 min. Run length = last exposure start − first exposure start. "Longest run" is per series.
- **Season.** A target's nights clustered by gaps > 120 d.
- **Calibrated night.** A night for which, on the same camera and binning, (a) a flat in the same
  filter exists within ±30 nights and (b) a bias or dark exists within ±30 nights. This is an
  *availability* statement only; it does not assert the calibrations are good.
- **Minimum-bearing night** (for an eclipsing/contact binary with a catalogue period P): a night whose
  series contains a run of length ≥ P/2 (which guarantees one minimum inside it whatever the
  ephemeris) **or**, where a catalogue ephemeris (epoch and period) is available, a run that contains
  a predicted primary or secondary minimum with ≥ 30 min of data on each side. Where no period is
  available the first clause is applied with run ≥ 3 h. Both counts are published separately; the gate
  uses their union.
- **Time-convention pass.** For a camera/software epoch, the header audit identifies what `DATE-OBS`
  means (start of exposure, UTC), an independent card (`JD`, `JD-OBS`, `MJD-OBS`, `TIME-OBS`,
  `DATE-AVG` or equivalent) agrees with it to < 2 s on ≥ 99% of frames that carry both, and there are
  no duplicate or non-monotonic timestamps within a run. A header audit can establish the
  *convention*; it cannot establish the *absolute* clock (that needs measured minima against
  literature, OA's 60 s criterion, which is photometry and outside a census).

## 3. The three candidate questions (named in advance, ED RIG-N1)

A "go" must be able to name a question, a reader and a venue. These are the only three that count;
they are written down now so that no question can be invented to fit what the census happens to show.

| | Question | Reader | Venue | Data gate |
|---|---|---|---|---|
| **Q1** | Do W UMa-type/eclipsing systems in the archive show secular period change over 2015–2022 (dP/dt ~ 10⁻⁷ d yr⁻¹ ⇒ ~7 min of O−C in ten years, detectable with ~1 min timings — PH)? | Observers and ephemeris maintainers of those systems (AAVSO eclipsing-binary section, BAV / O−C Gateway users) | JAAVSO (or IBVS-successor / OEJV) — a student-led timing paper | **G1** |
| **Q2** | Does the archive extend an RLMT-era baseline backwards — above all T CrB's 2015–2022 super-active phase in BVRI, or the polars' (ST LMi, VV Pup, EU UMa, AN UMa) accretion-state / timing history, or YZ Cnc? | Readers of the corresponding MACRO paper (T CrB monitoring; CV time series) | A section of that existing paper — **not** a sixth paper; RNAAS if it stands alone | **G2** |
| **Q3** | Is there one long, calibrated single-target photometric series that supports a stand-alone variability study (DS criterion)? | Named only if the target is known; otherwise the gate cannot pass on its own | AJ / JAAVSO according to target | **G3** |

## 4. Gates and the decision rule

Evaluated on canonical science frames, by script, in this order.

- **G0 (validity; must hold for any go).** Files-on-disk = rows + named exclusions; a camera/mechanical
  timeline exists with first/last night per camera; the time-convention pass (§2) holds for every
  camera epoch that contributes frames to a passing gate. If G0 fails, the outcome is **NOT DECIDABLE —
  fix the census**, never a go.
- **G1 (contact binaries, PH).** ≥ 3 eclipsing/contact systems each with minimum-bearing nights in
  ≥ 3 seasons. (1–2 systems: timings may be released as a table in the data note; not a go.)
- **G2 (overlap, DS/OA).** For at least one RLMT-era target: ≥ 10 nights in one series, spread over
  ≥ 2 seasons, in a filter the RLMT-era analysis of that target also uses or can be tied to (B, V, R, I,
  g, r, i or their Johnson/Sloan counterparts). For the timing targets (polars, YZ Cnc) at least one
  run ≥ 1 h is also required. Passing G2 is *conditional*: saturation at the target (standing rule 4)
  cannot be judged from headers and must be tested on pixels before any baseline is claimed.
- **G3 (single long series, DS).** ≥ 1 target with ≥ 30 calibrated nights in one series.

**Decision rule.**

1. G0 fails → **NOT DECIDABLE**.
2. G0 holds and none of G1–G3 passes → **NO-GO for a science project. Outcome = archive data release
   only**: the census database, a Zenodo record and a short data note (ED RIG-N1). L2-strategy is
   recorded as "not commissioned".
3. G0 holds and G2 passes (alone or with others) → **NO sixth project on that account**: the overlap
   frames are handed to the owning project (T CrB or CV) as a candidate baseline extension, subject to
   that project's own pixel gates. This is a transfer, not a go.
4. G0 holds and G1 or G3 passes → **GO-CANDIDATE**: the named question, reader and venue in §3 for the
   passing gate are stated, and a strategy *may* be commissioned. Under G3 the go additionally requires
   that the target's question and reader can be named; if they cannot, rule 2 applies. A GO-CANDIDATE
   is a recommendation; whether a sixth project is opened is James's decision, and per ED it must not
   compete with T CrB for effort in autumn 2026.

More than one of rules 3 and 4 can apply at once; each is reported.

## 5. What is forbidden after the census is seen

- Changing any threshold, definition or question in §§2–4. If a definition proves unworkable (for
  example a camera with no recorded filter card), the deviation is listed in the report as a
  deviation, with the result under both the written and the amended definition.
- Promoting a target that fails the gates because it "looks interesting". Such targets are listed in
  the census as *exploratory* and are explicitly not a reason for a go.
- Quoting frame counts as evidence of a usable series (DS RIG-L1-selection): only nights × filters ×
  longest run count.
- Treating header gain as measured (DE): calibration availability, including flat pairs, is reported;
  no gain is asserted by the census.
