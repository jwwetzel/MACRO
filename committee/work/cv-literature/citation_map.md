# Citation map — `manuscripts/CV_TimeSeries/main.tex` (as of 2026-08-20 15:06, 1,679 lines)

Package `cv-literature`, findings ED.E1 / RF.M8 / CV-R12. This package does **not** edit `main.tex`;
this file tells whoever does (CV-R13) which key goes where. Line numbers are those of the file as it
stands on disk today. Every key exists in `manuscripts/CV_TimeSeries/references.bib` and was verified
against its registry (`bib_verification.md`: 113 entries, 112 VERIFIED, 1 MANUAL).

Legend for the **Action** column: **CITE** = add the key(s), wording unchanged; **CITE+REWORD** = the
sentence is not correct as written once the literature is read, and the correction is stated;
**FIX** = an existing citation or credit is wrong.

AAS style keeps citations out of the abstract, so abstract uses (ll. 72–139) are cited at first use in
the body instead; they are listed here only where the body has no counterpart.

## A. Introduction (ll. 145–215)

| Line(s) | Uncited use | Key(s) | Action |
|---|---|---|---|
| 149–152 | what a polar is: field prevents a disc, locks the spin | `cropper1990`, `warner1995`, `wickramasinghe2000`, `ferrario2015`; discovery of the class: `tapia1977` | CITE |
| 152–157 | shock; cyclotron + bremsstrahlung; cyclotron emission "strongly beamed and strongly wavelength dependent"; colour encodes geometry and field | `chanmugam1981`, `wickramasinghe1985`, `wickramasinghe2000` | CITE |
| 157–160 | "Four decades of phase-resolved photometry and polarimetry — classically Cropper (1990) and Ferrario et al. (1993), and their successors — … ST LMi and VV Pup among them" | ST LMi: `stockman1983`, `schmidt1983`, `bailey1985`, `cropper1986`, `peacock1992`, `cropper1994`, `potter2000`, `campbell2008c`. VV Pup: `walker1965`, `visvanathan1979`, `schwope1997`, `howell2006`, `mason2008`, `campbell2008b` | CITE+REWORD. `ferrario1993` is not a general cyclotron paper: its real title is "Detection of cyclotron emission features in the infrared spectrum of ST LMi" (the bib entry carried an invented title; fixed). It is the ST LMi field-strength paper and should be cited as that. |
| 162–163 | "What those campaigns could not do is tag a cycle-resolved colour curve against the modern survey record" | `bailey1985`, `peacock1992` | CITE+REWORD. Simultaneous multiband single-night light curves of ST LMi exist (Bailey et al. 1985, optical+JHK; Peacock et al. 1992, "simultaneous multiband photometric and polarimetric light curves … on three occasions"). The defensible novelty is *number of nights × two seasons × an independent state tag*, not cycle-resolved colour as such. Same correction at ll. 173–176 and 1160–1163. |
| 164–165 | ZTF, ASAS-SN, AAVSO as data | `bellm2019`, `masci2019`; `shappee2014`, `kochanek2017`; `aavso_aid` | CITE |
| 165–168 | "enough to say whether a night was a high state or a low state" | `kafka2005`, `mason2017`, `duffy2022` | CITE (prior survey-based state work on these very stars; Duffy et al. tabulate ST LMi, AN UMa and EU UMa) |
| 172 | "published ZTF polar folds do" | `szkody2020`, `szkody2021`, `chen2020`, `rodriguez2023`, `vanroestel2025`; TESS period catalogues `dag2026`, `kepler2026` | CITE. `chen2020` is the sharpest example: VSX's EU UMa period and epoch *are* a ZTF fold. |
| 179–181 | "the SU UMa dwarf nova YZ Cnc … superhump timing" | `patterson1979`, `kato2009`, `warner1995` | CITE |
| 201–203 | external constants: "a catalogue period, an amplitude taken from the literature" | `watson2006`; `smak2010`, `kato2012` | CITE |

## B. Observations, times, photometry (ll. 217–535)

| Line(s) | Uncited use | Key(s) | Action |
|---|---|---|---|
| 223–224 | the telescope and site | none exists | No citable instrument paper (SYNTHESIS D4, James's decision). State that explicitly or cite the consortium's site. |
| 317–319 | mid-exposure BJD_TDB "computed from … the target coordinates" | `eastman2010`; implementation `astropy2013`, `astropy2018`, `astropy2022` | CITE |
| 331–336 | "eclipse timings of a detached eclipsing binary" | the binary's own ephemeris paper — **not in this bib**; owned by package `clock` (F-8, CV-R11) | CITE (key to be supplied by `clock`) |
| 344–345 | "aperture photometry on ensembles" | `honeycutt1992` (present). Aperture photometry and source detection are done with `sep`: `barbary2016`, `bertin1996`; frame registration with `astroalign`: `beroiz2020` | CITE |
| 87, 416 | ATLAS-REFCAT2 | `tonry2018` | CITE |
| 416–418 | (if the Gaia cross-identification / synthetic photometry of `macro_phot/gaia.py`, `cattie.py` is mentioned in CV-R9's reduction paragraph) | `gaia2016`, `gaia2023`, `montegriffo2023`; query layer `ginsburg2019` | CITE |
| 496, 1085 | "Sloan-era" $g,r,i$ | `tonry2018` (REFCAT2 is on the Pan-STARRS1 system; say so rather than "Sloan") | CITE+REWORD (minor) |

## C. Time-series analysis (ll. 537–795)

| Line(s) | Uncited use | Key(s) | Action |
|---|---|---|---|
| 541–542 | Lomb–Scargle | add `lomb1976`, `scargle1982` beside `vanderplas2018` | CITE |
| 557–559 | "We therefore adopt the published orbital ephemerides" | Table `tab:ephem` of `comparison_section.tex`: `watson2006` (all five, as retrieved), `cropper1986` (ST LMi), `walker1965` (VV Pup), `chen2020` (EU UMa), `bonnetbidaud1996`/`ok2025` (AN UMa), `shafter1988`+`downes2001`, `vanparadijs1994` (YZ Cnc) | CITE+REWORD. Two of the five are **not published ephemerides**: VSX's ST LMi period+epoch and AN UMa epoch were derived from AAVSO data by the VSX moderator (rev. 2025-11-12), with no error and no paper. Say "catalogue ephemerides (VSX)" and cite the published ones for comparison. |
| 564 | "Polar light curves have a sharp bright-phase edge" | `cropper1990`; for ST LMi `bailey1985`, `cropper1986`, `peacock1992` | CITE |
| 612–615, 905–906 | "if the cyclotron beaming is wavelength dependent then a bright-phase edge should sit at a slightly different phase in $i$ than in $g$" / "as it should" | `bailey1985` (bright-phase duration 0.29 P white, 0.31 J, 0.36 H, 0.39 K), `potter2000`, `cropper1994` | CITE+REWORD: give the predicted sign (redder later) and scale (PH.P3; `verify_items.md` V3) |
| 684–686 | "the epoch VSX publishes for a polar is a fiducial of its own" | `watson2006`; `cropper1986` | CITE+REWORD: for ST LMi the VSX fiducial is *undocumented*; Cropper's is the linear-polarisation pulse, mid-decline (PH.P4; `verify_items.md` V4) |
| 687–694 | "and independently, the steepest faintward gradient …" | — | REWORD (PH.P4): same data, same ephemeris — a consistency check of the edge fitter, not corroboration |
| 727–729 | "the accumulated phase drift at the catalogue period's quoted precision" | `cropper1986` | CITE+REWORD (RF minor 1): the catalogue quotes no precision; the only published σ_P is Cropper's 8×10⁻⁸ d, giving 0.022 cycles over 21,869 cycles, not 0.0014 |
| 738–741 | high/low accretion states | `kafka2005`, `kafka2007`, `mason2017`, `duffy2022` | CITE |

## D. Results (ll. 797–1150)

| Line(s) | Uncited use | Key(s) | Action |
|---|---|---|---|
| 801–803 | ST LMi, 113.9-min orbit | `stockman1983`, `cropper1986`, `watson2006` | CITE |
| 813–816 | bright phase in the same phase interval, "gradual rise and a sharp decline", amplitude "strongly band dependent" | `cropper1986`, `peacock1992`, `kafka2007`, `robertson2008` | CITE (all three properties are published; ours confirm them) |
| 824–831 | colour through the orbit | `kafka2007`, `cropper1986`, `peacock1992`, `ferrario1993`, `campbell2008c`, `wickramasinghe2000` | CITE (PH.P6; draft paragraph in `comparison_section.tex`) |
| 864–867 | "$|\dot P|<…$, or … s of orbital period per year" | `schmidt1995`, `campbell1999`, `beuermann2014`, `knigge2011` | CITE+REWORD (PH.P2): spin/spot-longitude period, not orbital; three scales |
| 903–906 | band-dependent edge "as it should" | `bailey1985` | CITE (see C) |
| 944–954 | VV Pup | `walker1965`, `howell2006`, `mason2008`, `campbell2008b`, `schwope1997` | CITE |
| 956–981 | EU UMa | `mittaz1992`, `howell1995`, `ramsay2004`, `chen2020`, `duffy2022` | CITE |
| 983–1057 | AN UMa | `krzeminski1977`, `liebert1982`, `bonnetbidaud1996`, `campbell2008b`, `ok2025`, `duffy2022` | CITE |
| 1017–1019 | "not distinguishable from AN UMa's own flickering" | `bonnetbidaud1996` (fast optical oscillations of AN UMa) | CITE |
| 1062–1064 | superhump timing in an SU UMa system | `patterson1979`, `kato2009`, `osaki1989`, `osaki1996` | CITE |
| 1068–1070 | "roughly `\NumSuperoutburstAmp` mag for a superoutburst" | `kato2002`, `patterson1979`, `vanparadijs1994`, `hakala2004` | FIX (PH.P5): 3.0 mag is a *normal*-outburst amplitude for YZ Cnc; the superoutburst is ≈4.0 mag (15.0→11.0). See `verify_items.md` V5; `numbers.tex:278`. |
| 1070–1074 | "independent AAVSO photometry" | `aavso_aid` | CITE |
| 1079–1081 | structure function | `simonetti1985` | CITE |
| 1097–1098 | "the catalogue period — which publishes no epoch" | `shafter1988`, `downes2001`, `watson2006`; `vanparadijs1994` | CITE+REWORD: a spectroscopic epoch exists (HJD 2447518.7255 ± 0.0061), but its error grows to ≈12 cycles by 2024, so the conclusion (phase zero arbitrary) stands |
| 1097–1125 | the quiescent orbital hump | `vanparadijs1994` (0.5 mag hump in quiescence, max at phase 0.8), `sun2026` | CITE (prior detection and TESS null; draft paragraph supplied) |
| 1120–1122 | "whose ephemeris drift is below 0.01 cycles" | `shafter1988`, `vanparadijs1994` | CITE+REWORD: under the published σ_P of the period actually used (0.0002 d) the one-day drift is 0.027 cycles; under van Paradijs et al. it is 0.0009. Either is ≪ 0.384, so the conclusion survives; the number does not. |
| 1123–1124 | "is real — and is not orbital" | `vanparadijs1994` | REWORD (RF.M7); note van Paradijs et al. saw the hump maximum move from phase 0.8 (quiescence) to 0.5 (late decline), so a phase change between nights is not by itself non-orbital |
| 1126–1133 | flickering | `moffett1974`, `vanparadijs1994`, `bruch2021` | CITE |
| 1140–1144 | "published superhump semi-amplitudes that start at `\NumSuperhumpFloorMmag` mmag … a measurement of absence" | `smak2010`, `kato2012`, `dai2026`, `kato2014v`, `sun2026` | CITE+REWORD (PH.P1/CV-R4): peak semi-amplitude ≈125–150 mmag; 50 mmag is the decayed late-plateau value; contour = smallest recoverable amplitude |

## E. Discussion, conclusions, back matter (ll. 1152–1572)

| Line(s) | Uncited use | Key(s) | Action |
|---|---|---|---|
| 1160–1167 | "what these data add … whether the accretion geometry is stable between states" | `bailey1985`, `cropper1986` (bright-phase duration changed 0.24→0.29 P and 0.33→0.38 P between epochs), `kafka2007`, `duffy2022` | CITE; insert `comparison_section.tex` after this subsection |
| 1169–1172 | "colour differences a change in the cyclotron beaming would produce" | `ferrario1993`, `campbell2008c` | CITE |
| 1220–1224 | "expected on cyclotron grounds" | `bailey1985` | CITE + predicted scale (standing rule 2) |
| 1272–1277 | "if one had, a superhump would have been visible" | `kato2009`, `dai2026` | REWORD (PH.P1) |
| 1503–1511 | conclusions item on superhumps | as ll. 1140–1144 | CITE+REWORD |
| 1553–1562 | acknowledgments: ZTF, ASAS-SN, ATLAS-REFCAT2, VSX named, none cited; Gaia not named at all although Gaia DR3 is queried (`pipeline/macro_phot/gaia.py`) | `bellm2019`, `masci2019`, `shappee2014`, `kochanek2017`, `tonry2018`, `watson2006`, `gaia2016`, `gaia2023` | CITE. The facilities' *mandatory acknowledgment wording* (ZTF, Gaia/DPAC, ASAS-SN, VSX, SIMBAD/VizieR) must be copied from each facility's own page at submission — not drafted here (needs James/chair). |
| 1565 | `\facilities{RLMT, PO:1.2m (ZTF), ASAS-SN, AAVSO}` | add `Gaia` | FIX |
| 1567–1571 | `\software{… photutils \citep{photutils} …}` | remove `photutils`; add `sep` (`barbary2016`, `bertin1996`), `astroalign` (`beroiz2020`), `astroquery` (`ginsburg2019`), `pandas` (`mckinney2010`); `astropy` → `astropy2013`, `astropy2018`, `astropy2022` | FIX. **photutils is not imported anywhere under `pipeline/` and is not installed in the `rlmt-checks` environment**; the photometry is `sep`. (`grep -rnE "^\s*(import\|from)\s+photutils" pipeline` → nothing; `import sep` in `macro_phot/extract.py`, `macro_core/astrom.py`, `scripts/run_cv_characterization.py`.) |

## F. Generated files (owned by the emitters, not by `main.tex`)

These are written by `pipeline/scripts/run_cv_paper.py` / `macro_phot/numbers_cv.py`; a citation there is
a change to an emitter, which this package does not own.

| File:line | Uncited use | Key(s) |
|---|---|---|
| `captions.tex:6` (Fig. 1) | "independent survey epochs (ZTF, ASAS-SN, AAVSO)" | `bellm2019`, `shappee2014`, `kochanek2017`, `aavso_aid` |
| `captions.tex:9` (Fig. 4) | "the published orbital frequency" | `watson2006` |
| `captions.tex:13` (Fig. 8) | "ZTF, ASAS-SN, AAVSO nightly means" | as Fig. 1 |
| `captions.tex:14` (Fig. 9) | "the catalogue PERIOD … the catalogue epoch" | `watson2006`, `cropper1986` |
| `captions.tex:15` (Fig. 10) | "the independent AAVSO nightly record" | `aavso_aid` |
| `tables.tex:81` (verdict row) | "NO — AND THAT IS A MEASUREMENT, NOT AN ABSENCE" | wording, PH.P1/CV-R4; `smak2010`, `kato2012` |
| `numbers.tex:25–47` | "[VSX catalogue, external constant]" | provenance note could carry `watson2006` + VSX revision date; ST LMi/AN UMa notes should say "from AAVSO data, no published error" |
| `numbers.tex:278` | `\NumSuperoutburstAmp{3.0}` "[ANALYSIS_STRATEGY §4, external constant]" | replace by a literature constant: 4.0 mag, `kato2002` (PH.P5) |
| `numbers.tex:281` | `\NumSuperhumpFloorMmag{50}` "[literature: …]" | name the literature: `smak2010`, `kato2012`; and add the peak value 125 mmag |

## G. Coverage check

Counted from the files, not by hand (`\cite` keys in `comparison_section.tex`; keys in this map; keys in
`main.tex`, each intersected with the bib):

- `main.tex` today cites **10** keys (one of them, `photutils`, wrongly).
- `comparison_section.tex` cites **62** distinct keys; it compiles against the bib with 0 undefined citations.
- This map places a further **≈ 40** keys that the section does not use (surveys, catalogues, software,
  time standard, the three audit polars, ZTF/TESS population papers).
- Applying both puts **more than 100 distinct references** in the paper, against the editor's floor of ≈ 40.
- Verified entries deliberately left unplaced, available if the text grows to need them: `ricker2015`
  (TESS mission), `ritter2003` (CV catalogue), `otulakowska2016` (outburst statistics), `zhao2005`,
  `verbunt1999`, `shahbaz1996` (ST LMi system parameters), `yates2026` (eclipse-time variations in
  white-dwarf binaries; the strategy's "WD-binary timing literature"). Unused entries cost nothing in BibTeX.
