# The "(verify)" items for the CV paper, answered from the literature

Package `cv-literature`. Owners: PH.P2, PH.P3, PH.P5, PH.P6, RF minor 1, plus the ephemeris-fiducial
question of the brief (which is PH.P4's open half).

**How to read the "read how" column.** *Direct* = I read the statement in the paper's own text (arXiv
full text, or the ADS scan with its OCR layer, or the scanned page image). *Quoted* = I read it in a later
paper that quotes it; the original was not opened. *Catalogue* = a database record. *Summary* = taken
from a fetched-page summary, the weakest grade; confirm at source before submission. Every number below
is carried, with its locator, in `literature_scales.py` and emitted to `literature_scales.md`; the
arithmetic is done there and nowhere else.

---

## Summary table

| Id | Committee statement to verify | Verdict | Cited value | Key | Read how |
|---|---|---|---|---|---|
| V1 | PH.P2: "an asynchronous polar relaxing like V1500 Cyg has Ṗ_spin ~ 4×10⁻⁸" | **Confirmed** | Ṗ = 3.86×10⁻⁸ (polarimetry, 1987–1992); synchronisation time ≈ 170 yr | `schmidt1995` | Quoted (in `pavlenko2018` §1) |
| V2 | PH.P2: secular orbital evolution below the gap is "\|Ṗ\| ~ 10⁻¹³–10⁻¹⁴ (τ ~ 10⁹ yr)" | **Confirmed as arithmetic** | P/τ for τ = 1–5 Gyr → 4×10⁻¹⁴ – 2×10⁻¹³ | `knigge2011` (for the Gyr timescale) | Derived; the timescale is the standard result of that paper, not a number I read off a table |
| V3 | PH.P3: the band-offset null needs an expected size; "chromatic shifts … plausibly tens of seconds" | **Confirmed, with a published measurement of the effect and its sign** | Bright-phase duration at half maximum: 0.29 P (white light), 0.31 P (J), 0.36 P (H), 0.39 P (K), simultaneous. Half-differences: 68 s (white→J), 273 s (J→K). Fall takes 0.06–0.07 P = 410–478 s. | `bailey1985` §7.3; `peacock1992`; `cropper1986` §4 | Direct (ADS scans) |
| V4 | PH.P4 / brief: what does each VSX epoch mark? | **Answered for 4 of 5; ST LMi is undocumented** | see §V4 table | `watson2006` + originals | Catalogue + Direct |
| V5 | PH.P5: "for YZ Cnc ~3 mag is roughly a normal outburst amplitude; superoutbursts are nearer 3.5–4" | **Confirmed** | Superoutburst: 15.0 → 11.0 = 4.0 mag (`kato2002` Table 3); peak 10.6 in 1988–89 (`vanparadijs1994`). Normal: 15.2 → 12.2 = 3.0 mag (`hakala2004`); peak ≈ 12 (`vanparadijs1994`). Supermaxima ≈ 1 mag brighter than normal maxima (`patterson1979`). | as listed | Direct |
| V6 | PH.P6: "ST LMi's ~12 MG field" | **Confirmed** | 12.1 ± 0.5 MG (`campbell2008c`, direct); 11.5 ± 0.5 / ≈12 MG (`ferrario1993`, quoted); second pole 30 ± 5 MG (`schmidt1983`, quoted) | as listed | Direct / Quoted |
| V7 | RF minor 1: cycle-count uniqueness rests on an *assumed* σ_P | **Confirmed, and now replaceable by a published σ_P** | VSX publishes no error (its ST LMi period and epoch are "from AAVSO data", rev. 2025-11-12). The only published error is Cropper's: P = 0.07908908 ± 0.00000008 d → drift over 21,869 cycles = **0.022 cycles**, not 0.0014. Still unique, by a factor of 23. | `cropper1986` §3 | Direct |

Two further items surfaced while verifying; they are not on the committee's list but change sentences in
the paper: **V8** (a more precise YZ Cnc period and an epoch exist; a prior detection of the orbital hump
exists) and **V9** (the 50 mmag superhump "floor" is a late-plateau value; the peak is 125–150 mmag).

---

## V1–V2. What the Ṗ bound constrains (PH.P2)

The paper's bound, `\NumStLmiPdotLimit` = 3.6×10⁻⁹, against three scales:

| Scale | Value | Source | Bound / scale |
|---|---|---|---|
| (i) secular orbital evolution, τ = 1–5 Gyr | \|Ṗ\| = 4×10⁻¹⁴ – 2×10⁻¹³ | arithmetic P/τ; `knigge2011` | bound is 10⁴–10⁵ × too weak |
| (ii) resynchronising white dwarf (V1500 Cyg) | Ṗ_spin = 3.86×10⁻⁸ | `schmidt1995` via `pavlenko2018`; context `stockman1988`, `schmidt1991`, `campbell1999` ("fastest about 50 years"), dispute `harrison2016`/`harrison2018` vs `pavlenko2018` | bound is **11× below** → excluded |
| (iii) libration of a synchronised white dwarf (DP Leo: 50° in 1979–2001, then stopped; ≈ 60 yr period, ≈ 25° amplitude) | peak curvature ≡ \|Ṗ\| ≲ 4×10⁻¹¹ at ST LMi's period | `beuermann2014` abstract (direct); arithmetic: A(2π/T)²/360 × P² | bound is 100× too weak |

So the physicist's reading is right on every count: the bound is a spin/spot-longitude limit, it excludes
only V1500 Cyg-like relaxation, and "the result the data actually support" is a longitude statement.
For that statement the literature supplies the comparison scale, which the memo did not have:

- our O−C rms `\NumStLmiOcRmsS` = 84 s = 0.0123 cycle = **4.4°** (arithmetic check of the memo: agrees);
- ST LMi's bright phase has changed length by 0.05 P between epochs **twice** in the record — 0.33 P (1982)
  → 0.38 P (1985) (`cropper1986` §4), and 0.24 P (1982 June) → 0.29 P (1983 Mar–May) (`bailey1985` §8).
  A symmetric change moves one edge by 0.025 P = **171 s = 9°**.

So an edge that holds to 4° over 1.9 yr is *tighter* than the star's own published history, and whether
it holds across accretion states (the P2 "close" condition) is a genuinely informative test, not a
formality. That number is CV-R5's to produce.

**Not done, flagged for CV-R3/R5:** Cropper's polarimetric phase zero is the *same feature* we time
(mid-decline), 43.7 yr earlier. If the cycle count across that gap is accepted, a single O−C point at
E ≈ −180,034 would tighten the curvature limit by orders of magnitude (½ Ṗ E² P against a phase
tolerance of ~0.2 cycle gives \|Ṗ\| ≲ 10⁻¹¹). It is **not** accepted here: under Cropper's own σ_P the
count across the gap is uncertain by 0.18 cycle at 1σ (unique at 1σ, not at 3σ), the two feature
definitions are not proven identical, and HJD-vs-BJD differs by up to ~1 min. It is recorded as a
candidate analysis with those three conditions, not as a result.

## V3. Expected size and sign of a band-dependent edge (PH.P3)

The memo argued "limb occultation of a compact region is achromatic to first order … plausibly tens of
seconds". The literature is stronger than that: the effect has been **measured in this star**.

- `bailey1985` §7.3, simultaneous data: bright-phase duration at half maximum **0.29 P white light, 0.31 P
  in J, 0.36 P in H, 0.39 P in K**; "the emitting regions at the different wavelengths must be at different
  heights", a difference of 0.01 R_wd sufficing.
- A longer bright phase ends later. Half-differences (one edge): white→J **68 s**; J→K **273 s**.
- The region is extended: ±15–18° of magnetic longitude (`potter2000`, *summary grade*), with structure
  (`cropper1994`, `stockman1996`); the fall itself takes 0.06 P (`peacock1992`) to 0.07 P (`cropper1986`),
  i.e. **410–478 s**.
- A log-linear interpolation of Bailey's durations to g and i gives **≈ 37 s, redder band later** — order
  of magnitude only: the infrared trend is an optical-depth effect and g, r, i sit at optically thin
  harmonics 18, 14, 12.

**Consequence for D3.** The predicted sign is *g edge earlier than i*, i.e. g − i **negative**, and the
predicted size is tens of seconds up to ~10² s. The chair's U3 quotes the paper's own edges as
g − i ≈ −110 ± 28 s with 9/9 or 11/12 pairs negative. That is the literature's sign and within a factor
of ~3 of the interpolated size. This does **not** decide D3 — an estimator bias with per-band ramp widths
could have the same sign (a wider, shallower g ramp is exactly what a smaller g amplitude produces) — but
it removes the option of calling a negative g − i "unexpected", and it means the sentence at
`main.tex:903–906` ("as it should") can cite a measurement, with its scale beside it (standing rule 2).

## V4. What each adopted epoch marks

Catalogue values retrieved 2026-10-03 from VizieR `B/vsx/vsx` and the VSX detail pages (`vsx.aavso.org`).

| Target | VSX P (d) | VSX T₀ (HJD) | T₀ marks | Origin of the VSX values | Best published ephemeris |
|---|---|---|---|---|---|
| ST LMi | 0.07908912 | 2459298.4236 | **not documented** | VSX rev. 5, 2025-11-12, S. Otero: "Period and epoch from AAVSO data". No paper, no error. | `cropper1986`: HJD 2445059.7024(3) + 0.07908908(8) E; zero = **linear-polarisation peak, "half way through the intensity's rapid decline"** (direct) |
| VV Pup | 0.0697468256 | 2427889.6474 | **maximum light** (VSX remark: "The epoch of light maximum is given"; identical in GCVS) | `walker1965` | same; `howell2006` §3 confirms it is HJD and "the time of high state maximum light" (direct) |
| EU UMa | 0.0626194 | 2458231.7799 | **minimum light** (inferred: VSX rev. 5, 2021-10-26, "Updated from Chen 2020 catalog"; Chen et al. define T₀ as "the minimum" of the fitted ZTF light curve) | `chen2020` | discovery/EUVE: `mittaz1992`, `howell1995` (values not re-read) |
| AN UMa | 0.07975274 | 2456725.4053 | **minimum light** (VSX remark: "The epoch of Min for periodic … light variations … is given") | epoch: VSX rev. 4, 2025-11-12, "Epoch from AAVSO data". Period: origin not documented — it is **not** the GCVS value (GCVS: 0.0797522 d, epoch 2442891.397) | `bonnetbidaud1996`: 2443190.9921(2) + 0.07975282(4) E, T₀ = linear-polarisation pulse of `liebert1982`; `ok2025`: BJD 2443190.9926(2) + 0.079752867(12) E |
| YZ Cnc | 0.0868 | none | — | `shafter1988` 0.0868(2) d via `downes2001` (VSX rev. 3) | `vanparadijs1994`: HJD 2447518.7255(61) + 0.086924(7) E, zero "probably" superior conjunction of the white dwarf (direct) |

**Where the ST LMi catalogue zero falls on the published ephemerides** (arithmetic on catalogue constants
only, in `literature_scales.py`):

- on Cropper's polarimetric ephemeris: phase **0.97 ± 0.18**; our edge (zero + 0.157) at **0.13 ± 0.18**,
  against Cropper's pulse at 0.0 and end-of-decline at 0.05;
- on the inferior-conjunction ephemeris of `kafka2007`/`robertson2008` (HJD 2450596.9006 + 0.07908908 E,
  *summary grade*): phase **0.80 ± 0.11**; our edge at **0.96 ± 0.11**, against their peak near 0.78 and
  sharp drop near 0.87. (Independent cross-check: Kafka et al. state conjunction zero = polarimetric phase
  0.17; 0.80 + 0.17 = 0.97. The two routes agree.)

Reading: the catalogue zero sits at the bright-phase maximum, and the 1,071 s offset is the distance from
maximum to the end of the bright phase. That is a *consistency*, with ±0.1–0.2 cycle of extrapolation
error; the fiducial itself is known only to the VSX moderator. **Action for James/chair:** one e-mail to
VSX (S. Otero) asking what the ST LMi and AN UMa epochs of 2025-11-12 mark and what data they used would
turn "not documented" into a citation.

Note also PH.P4's first half stands unchanged: the folded-profile gradient at phase 0.16–0.19 is the same
data on the same ephemeris and checks the edge fitter, not the offset.

## V5. YZ Cnc outburst amplitudes (PH.P5)

`\NumSuperoutburstAmp` = 3.0 ("ANALYSIS_STRATEGY §4") is a normal-outburst amplitude.

| Quantity | Value | Source |
|---|---|---|
| quiescence | V ≈ 15.0 (`kato2002`), 15.2 (`hakala2004`), ≈ 14.5 (`zhao2005`) | direct |
| normal maximum | ≈ 12 (`vanparadijs1994`), 12.2 (`hakala2004`); visual peaks 11.7–13 in `patterson1979` Fig. 2 (read off the scan) | direct |
| supermaximum | 11.0 (`kato2002`, `patterson1979` Table II: V = 11.0), 10.6 (`vanparadijs1994`), ≈ 10.5 (`zhao2005`) | direct |
| **normal amplitude** | **≈ 3.0 mag** | 15.2 − 12.2 |
| **superoutburst amplitude** | **≈ 4.0 mag** | 15.0 − 11.0 |
| super − normal | ≈ 1 mag (`patterson1979`), 1.4 mag (`vanparadijs1994`) | direct |
| recurrence | normal 10.3 ± 2.5 d, super 134 ± 19 d (`patterson1979` Table I); TESS: normal 8.1–9.5 d (`sun2026`) | direct |

The classification in the paper (brightest dense run 1.86 mag above quiescence ⇒ not a superoutburst)
survives, and with more room than the text claims: 1.86 mag is short even of a normal maximum. The
constant should become 4.0 mag with `kato2002` as its source.

## V6. Field strength and the physics of the colour curve (PH.P6)

- B = **12.1 ± 0.5 MG**, kT = 3.3 keV, i = 55°, β = 128°, cyclotron harmonics n = 4–7 seen at 2.25–1.30 µm
  (`campbell2008c`, direct). `ferrario1993`: 11.5 ± 0.5 MG as quoted by `kafka2007`, "≃ 12.0 MG" as quoted
  by `campbell2008c`. Second pole ≈ 30 MG (`schmidt1983`, quoted in `cropper1986` and `kafka2007`).
- Fundamental λ_c = 10710 Å × (10⁸ G / B) = **8.9 µm**; harmonic numbers at our bands: **i ≈ 12, r ≈ 14,
  g ≈ 18**. High, optically thin harmonics; intensity falls steeply with n, so the cyclotron component is
  much redder than the white dwarf + stream it sits on.
- Published sign and size of the colour change: bright phase "very red" (`cropper1986` §1, §4); orbital
  amplitude 2.0 mag in R_c (`peacock1992`, direct), ≈ 1.2 mag in V (`kafka2007`, direct), and B/R/I/J/H =
  0.6/1.4/2.0/1.5/1.3 mag (`peacock1992` as quoted by `campbell2008c`); "it becomes redder near phase 0.8
  … due to the contribution of the cyclotron continuum" (`kafka2007`).

So "why would a physicist care": the g − r swing is the ratio of harmonic ~18 to harmonic ~14 switching on
and off as the pole crosses the limb; its amplitude is set by how steeply the optically thin spectrum falls
(temperature and the size parameter Λ), and its *state dependence* is the new measurement.

## V7. The cycle count (RF minor 1)

| Quantity | Paper | With the published σ_P |
|---|---|---|
| σ_P | 0.5×10⁻⁸ d, assumed from the last digit of a VSX string | 8×10⁻⁸ d (`cropper1986`) |
| drift over 21,869 cycles | 0.0014 cycle | **0.022 cycle** |
| margin to half a cycle | ≈ 350× | **23×** |
| P_VSX − P_Cropper | — | +4×10⁻⁸ d = 0.5 σ_Cropper |
| P_refit − P_Cropper | — | −0.6 σ (combined) |

The count is unique either way; the text should quote the second column and say the first was an
assumption. Note what the VSX period *is*: a value fitted by a database moderator to AAVSO data in 2025,
with no error bar. Is the paper's agreement with it circular? The project's own record says no:
`cv_external` (in `products/phot/cv_timeseries.sqlite`, read-only query
`SELECT target, source, n_points, n_independent, notes FROM cv_external WHERE target='stlmi'`) holds 9,014
AAVSO rows for ST LMi of which **0** are flagged as RLMT data resubmitted to AAVSO (against 1,499 of 99,732
for YZ Cnc). So the catalogue period rests on other observers' photometry of the same years, and the
1.5σ agreement of our refit with it is a real, if errorless, check. The comparison that carries an error
bar is with Cropper's period, 40 years earlier, and it agrees at 0.6σ.

## V8. YZ Cnc: period, epoch, and a prior hump detection (new)

- A period 30× more precise than the one used exists: **0.086924 ± 0.000007 d** with epoch HJD 2447518.7255
  ± 0.0061 (`vanparadijs1994` §6, direct). VSX carries 0.0868 d (= `shafter1988`, ± 0.0002) and no epoch.
- `main.tex:1097` "which publishes no epoch": true of VSX, not of the literature. The published epoch
  propagates to 2024 with ± 12 cycles, so absolute phase is still unavailable — conclusion unchanged.
- `main.tex:1121` "ephemeris drift is below 0.01 cycles": under the published σ_P of the period actually
  used the one-day drift is **0.027 cycles**; the two periods differ by 0.016 cycles per day. Both ≪ the
  measured 0.384-cycle shift — conclusion unchanged, number wrong.
- **Prior detection.** `vanparadijs1994` §6–7: orbital modulation of **0.5 mag** full amplitude in
  quiescence, maximum at phase 0.8 (bright spot); ≈ 0.2 mag with maximum at phase 0.5 at late decline. Our
  fitted 30–70 mmag semi-amplitude is 4–8× below their 250 mmag. And the hump's phase *moved* between
  outburst stages in their data, which bears on RF.M7 ("coherent within a night and incoherent between
  nights is real — and is not orbital"): a phase change between nights at different outburst stages is
  what they saw in a modulation they call orbital.
- TESS: "no orbital-period signal was found" in the superoutburst sectors (`sun2026` §3.16); "quiescent and
  normal-outburst intervals lack a comparably persistent SH-band signal" (`dai2026` abstract).

## V9. Superhump amplitudes: what "floor" means (bears on PH.P1 / CV-R4)

| Stage | Full amplitude | Semi-amplitude | Source |
|---|---|---|---|
| maximum, low inclination | ≈ 0.25 mag | 125 mmag | `smak2010` abstract; `kato2012` §4.7 (sample mean 0.25 mag) |
| YZ Cnc, precursor / early plateau (TESS) | ≈ 0.3 mag | 150 mmag | `dai2026` §4.2 |
| YZ Cnc, post-plateau (TESS) | ≈ 0.5 mag | 250 mmag | `dai2026` §4.2 |
| classical "fully developed" | 0.3–0.4 mag | 150–200 mmag | Warner (1995), as quoted in `kato2012` |
| late plateau (decayed) | not read at source | the paper's 50 mmag | `kato2012` §4.7 states only that amplitudes "globally decrease with time" after the maximum; I did not read a late-plateau value off their Figs 85–88. The 50 mmag "floor" (`final_science.py:100–105`: "common superhumps run ~0.05–0.15 mag") therefore has **no citation yet**; its natural source is those figures. |

The paper's contour (82–237 mmag) therefore brackets the published *peak* semi-amplitudes (125–150 mmag)
and lies entirely above its own 50 mmag floor: it could have seen a fresh superhump on some run-filters and
a decayed one on none. One
further caution for P1's "it is also moot physically: the runs are in normal outbursts, where no superhump
is expected": not quite. In YZ Cnc itself late superhumps "survived during the next normal outburst" after
the 2011 superoutburst (`kato2014v` §3.5), and TESS shows persistent superhumps in normal outbursts in
three systems (`sun2026`). Whether any RLMT run is the first normal outburst after a superoutburst is
answerable from the AAVSO record the project already holds.
