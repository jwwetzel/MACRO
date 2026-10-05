# CV paper: macro-by-macro diff of numbers.tex (CV-R14)

Emitted by `committee/work/cv-revision/diff_numbers.py`. Do not edit.

- A: `committee/work/cv-stats/numbers.tex.before-cv-stats` (307 macros) — the first draft
- B: `committee/work/cv-revision/numbers.tex.before-cap` (1649 macros) — after the F-10 rebuild and the cv-stats emission
- C: `manuscripts/CV_TimeSeries/numbers.tex` (1934 macros) — now

Changed 48, added 1627, removed 0; unexplained steps: 0.

## Pre-existing macros that changed or disappeared

| macro | A | B | C | step: reason |
|---|---|---|---|---|
| `NumBandOffsetAbsMaxS` | 137 | 137 | 131 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumBandOffsetPooledRangeS` | -137 to 59 | -137 to 59 | -131 to 59 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumBandOffsetWeakestBoundCycles` | 0.046 | 0.046 | 0.045 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumBandOffsetWeakestBoundS` | 313 | 313 | 308 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumBandPairsSignificant` | 0 | 1 | 1 | A->B: CV-S9 v1.1: band pairs use scatter-based errors (CV-R1); the g-i pair is now significant |
| `NumBlindContourRangeMmag` | 82--237 | 82--237 | 82--226 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumCloudVetoedFrames` | 850 | 850 | 852 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumCompStarsMedian` | 228 | 228 | 225 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumDetrendContourSeries` | 5 | 5 | 6 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumDetrendContourSeriesBetter` | 4 | 4 | 5 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumExtinctionBoundRangeMmag` | 1.2--36.2 | 1.1--40.3 | 1.1--19.1 | A->B: CV-S8 re-run in the clean rebuild on the re-tied series; B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumExtinctionMaxEffectMmag` | 1.9 | 1.9 | 1.6 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumExtinctionSignificant` | 2 | 2 | 3 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumHgReadNoiseE` | 4.15 | 4.15 | 4.18 | B->C: detector package re-emitted S2 tables |
| `NumHgReadNoiseEErr` | 2.30 | 2.30 | 0.02 | B->C: detector package re-emitted S2 tables |
| `NumHumpAmpRangeMmag` | 30--70 | 30--70 | 30--71 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumHumpFieldContourMmag` | 9--52 | 9--52 | 9--43 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumLightcurveRows` | 3\,174\,618 | 3\,174\,618 | 3\,158\,147 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumMacrosExternal` | 25 | 59 | 71 | A->B: same (literature scales and method constants); B->C: count of macros: new checks stages (reduction, ramp, mech, clock, cap, detector) added macros |
| `NumMacrosPhot` | 251 | 1\,559 | 1\,832 | A->B: same; B->C: count of macros: new checks stages (reduction, ramp, mech, clock, cap, detector) added macros |
| `NumMacrosTotal` | 307 | 1\,649 | 1\,934 | A->B: count of macros: the revision added the NumRv set; B->C: count of macros: new checks stages (reduction, ramp, mech, clock, cap, detector) added macros |
| `NumOutburstRateSignificant` | 8 | 8 | 7 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumPrecisionSampleRange` | 20--815 | 20--815 | 20--824 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumPrecisionScopeClause` | The crosses at the CV's own magnitude and the noise floor annotated on each panel are both LOCAL FITS over the ensemble and check stars together, not held-out statistics: the held-out quantities in this paper are the check-star scatter $\sigma_{\rm chk}$, the catalogue-tie accuracy and the error-bar inflation factor | The crosses at the CV's own magnitude and the noise floor annotated on each panel are both LOCAL FITS over the ensemble and check stars together, not held-out statistics: the held-out quantities in this paper are the check-star scatter $\sigma_{\rm chk}$, the catalogue-tie accuracy and the error-bar inflation factor | The crosses at the CV's own magnitude and the noise floor annotated on each panel are both local fits over the ensemble and check stars together, not held-out statistics: the held-out quantities in this paper are the check-star scatter $\sigma_{\rm chk}$, the catalogue-tie accuracy and the error-bar inflation factor | B->C: text only: emphasis lowered to sentence case (seat 6, 2026-10-05), so the clause matches its caption |
| `NumStLmiEdgeFormalBarMinS` | 4.1 | 4.1 | 8.0 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiFallRiseRatioRange` | 1.4--3.0 | 1.4--3.0 | 1.3--3.0 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiFittedPeriodSigmaD` | 6.4 \times 10^{-8} | 6.0 \times 10^{-8} | 5.9 \times 10^{-8} | A->B: CV-S9 v1.1: the one-sided max(chi2_nu, 1) error rescaling was removed (standing rule 1); B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiOcChisq` | 0.92 | 0.92 | 0.91 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiOcChisqEdge` | 0.94 | 0.94 | 0.74 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiOcChisqHighGain` | 0.80 | 0.80 | 0.75 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiOcOffsetS` | 1\,071 | 1\,071 | 1\,072 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiOcRangeS` | -239 to 175 | -239 to 175 | -240 to 175 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiPeriodAgreementSigma` | 1.5 | 1.6 | 1.7 | A->B: same; B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStLmiSigmaNightEdgeMedianS` | 104 | 104 | 103 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStateNightsHigh` | 63 | 63 | 65 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStateNightsIntermediate` | 24 | 24 | 22 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumStateSeparabilityBimodalRange` | 0.76--0.99 | 0.76--0.99 | 0.75--0.99 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumSuperoutburstAmp` | 3.0 | 4.0 | 4.0 | A->B: literature correction (PH.P5): 4.0 mag superoutburst amplitude of YZ Cnc, kato2002 |
| `NumTieAtGoal` | 3 | 3 | 4 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieAtGoalUnclipped` | 1 | 1 | 2 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieBarRangeClippedMmag` | 28--38 | 28--38 | 24--38 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieBarRangeUnclippedMmag` | 56--200 | 56--200 | 59--200 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieBlocksClipped` | 19 | 19 | 18 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieClippedStarsRange` | 1--16 | 1--16 | 1--15 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieMedianRatioToGoal` | 1.2 | 1.2 | 1.3 | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |
| `NumTieSeries` | 26 | 27 | 27 | A->B: CV-S6 re-run in the clean rebuild: the 1 MHz y block now reaches the tie stage and fails it (counted as a primary block) |
| `NumTieUntied` | 1 | 2 | 2 | A->B: same: that block is the second untied primary block |
| `NumTieUntiedBlock` | EU~UMa's Fast $g$ block | EU~UMa's Fast $g$ block | EU~UMa's Fast $g$ block and ST~LMi's 1MHz HS 16-bit $y$ block | B->C: photometry chain re-run with the S2 1% linearity cap (run_cv_photometry.py recap): comparison and target measurements above the cap are withheld, so ensembles, ties and every downstream product move |

## Macros added by the revision

1627 macros did not exist in A. Of these, 1342 existed in B; 287 of those changed value between B and C (listed below), and 285 are new in C (the checks stages).

### Changed between B and C

| macro | B | C |
|---|---|---|
| `NumRvBiasVOneEraSevenGS` | -27 | -21 |
| `NumRvBiasVOneEraSevenRErrS` | 8 | 9 |
| `NumRvBiasVOneEraSevenRS` | 15 | 17 |
| `NumRvBiasVOneMinS` | -27 | -23 |
| `NumRvBiasVTwoEraSevenGS` | 11 | 16 |
| `NumRvBiasVTwoEraSevenIErrS` | 5 | 4 |
| `NumRvBiasVTwoEraSevenRErrS` | 6 | 7 |
| `NumRvBiasVTwoEraSevenRS` | 9 | 19 |
| `NumRvBiasVTwoMaxS` | 21 | 22 |
| `NumRvColGiEraSevenInterpN` | 137 | 135 |
| `NumRvColGiEraSevenInterpSeeingN` | 137 | 135 |
| `NumRvColGiEraSevenSixHundredN` | 214 | 212 |
| `NumRvColGiEraSevenSixHundredPhaseBlueErr` | 0.03 | 0.04 |
| `NumRvColGiRepeatInterpAmpDiffErrMag` | 0.09 | 0.08 |
| `NumRvColGiRepeatInterpChinu` | 2.80 | 2.73 |
| `NumRvColGiRepeatInterpSeeingChinu` | 3.21 | 2.97 |
| `NumRvColGiRepeatSixHundredChinu` | 2.91 | 2.65 |
| `NumRvColGrEraSevenInterpN` | 121 | 119 |
| `NumRvColGrEraSevenInterpSeeingN` | 121 | 119 |
| `NumRvColGrEraSevenSixHundredAmpErrMag` | 0.07 | 0.06 |
| `NumRvColGrEraSevenSixHundredAmpMag` | 0.74 | 0.75 |
| `NumRvColGrEraSevenSixHundredN` | 212 | 210 |
| `NumRvColGrRepeatInterpChinu` | 2.78 | 2.34 |
| `NumRvColGrRepeatInterpDof` | 15 | 14 |
| `NumRvColGrRepeatInterpR` | 0.96 | 0.97 |
| `NumRvColGrRepeatInterpScale` | 0.89 | 0.90 |
| `NumRvColGrRepeatInterpSeeingChinu` | 2.79 | 1.97 |
| `NumRvColGrRepeatInterpSeeingDof` | 15 | 14 |
| `NumRvColGrRepeatInterpSeeingR` | 0.96 | 0.97 |
| `NumRvColGrRepeatInterpSeeingScale` | 0.89 | 0.90 |
| `NumRvColGrRepeatSixHundredAmpDiffMag` | -0.01 | -0.02 |
| `NumRvColGrRepeatSixHundredChinu` | 3.62 | 3.39 |
| `NumRvColGrRepeatSixHundredScale` | 0.80 | 0.79 |
| `NumRvColRiEraSevenInterpAmpMag` | 0.40 | 0.39 |
| `NumRvColRiEraSevenInterpPhaseBlue` | 0.81 | 0.18 |
| `NumRvColRiEraSevenInterpPhaseBlueErr` | 0.28 | 0.43 |
| `NumRvColRiEraSevenInterpSeeingAmpMag` | 0.40 | 0.39 |
| `NumRvColRiEraSevenInterpSeeingPhaseBlue` | 0.81 | 0.18 |
| `NumRvColRiEraSevenInterpSeeingPhaseBlueErr` | 0.26 | 0.42 |
| `NumRvColRiEraSevenSixHundredN` | 252 | 251 |
| `NumRvColRiEraSevenSixHundredPhaseBlueErr` | 0.09 | 0.07 |
| `NumRvColRiRepeatInterpChinu` | 5.33 | 4.60 |
| `NumRvColRiRepeatInterpRmsMag` | 0.10 | 0.09 |
| `NumRvColRiRepeatInterpScale` | 0.58 | 0.59 |
| `NumRvColRiRepeatInterpSeeingChinu` | 3.75 | 3.44 |
| `NumRvColRiRepeatInterpSeeingRmsMag` | 0.10 | 0.09 |
| `NumRvColRiRepeatInterpSeeingScale` | 0.58 | 0.59 |
| `NumRvColRiRepeatSixHundredChinu` | 3.39 | 3.27 |
| `NumRvColRiRepeatSixHundredScale` | 0.17 | 0.16 |
| `NumRvColRiRepeatStrictAmpDiffMag` | -0.11 | -0.10 |
| `NumRvColRiRepeatStrictChinu` | 6.01 | 5.81 |
| `NumRvDbiasMagcGiEqualMaxS` | 2 | -6 |
| `NumRvDbiasMagcGiEqualMinS` | -42 | -43 |
| `NumRvDbiasMagcGiMinS` | -60 | -62 |
| `NumRvDbiasMagcGiSeMaxS` | 10 | 6 |
| `NumRvDbiasMagcGrEqualMinS` | -29 | -28 |
| `NumRvDbiasMagcGrMinS` | -52 | -51 |
| `NumRvDbiasMagcGrSeMaxS` | 5 | 6 |
| `NumRvDbiasMagcRiEqualMaxS` | 7 | 9 |
| `NumRvDbiasMagcRiEqualMinS` | -10 | -9 |
| `NumRvDbiasMagcRiMaxS` | 46 | 45 |
| `NumRvDbiasMagcRiMinS` | -60 | -58 |
| `NumRvDbiasMagcRiSeMaxS` | 6 | 7 |
| `NumRvDbiasVOneGiMaxS` | 25 | 26 |
| `NumRvDbiasVOneGrEqualMinS` | -36 | -35 |
| `NumRvDbiasVOneGrMinS` | -79 | -77 |
| `NumRvDbiasVOneRiMaxS` | 66 | 64 |
| `NumRvDbiasVTwoGiMaxS` | 4 | 0 |
| `NumRvDbiasVTwoGiSeMaxS` | 8 | 6 |
| `NumRvDbiasVTwoRiEqualMaxS` | 4 | 6 |
| `NumRvDbiasVTwoRiEqualMinS` | -11 | -8 |
| `NumRvDbiasVTwoRiMaxS` | 8 | 9 |
| `NumRvDbiasVTwoRiMinS` | -11 | -8 |
| `NumRvDbiasVTwoRiSeMaxS` | 6 | 7 |
| `NumRvDthreeMagcGiResidErrS` | 19 | 17 |
| `NumRvDthreeMagcGiResidS` | -34 | -31 |
| `NumRvDthreeMagcRiResidErrS` | 25 | 23 |
| `NumRvDthreeMagcRiResidNsigma` | 2.0 | 1.9 |
| `NumRvDthreeMagcRiResidS` | 49 | 44 |
| `NumRvDthreeVOneGiResidErrS` | 32 | 33 |
| `NumRvDthreeVOneGrResidErrS` | 29 | 30 |
| `NumRvDthreeVOneGrResidNsigma` | -1.2 | -1.1 |
| `NumRvDthreeVOneGrResidS` | -36 | -34 |
| `NumRvDthreeVTwoGiResidErrS` | 22 | 21 |
| `NumRvDthreeVTwoGiResidNsigma` | -3.7 | -3.9 |
| `NumRvDthreeVTwoRiResidErrS` | 26 | 24 |
| `NumRvDthreeVTwoRiResidNsigma` | 0.1 | -0.2 |
| `NumRvDthreeVTwoRiResidS` | 2 | -5 |
| `NumRvFitVOneChinuMin` | 1.7 | 1.6 |
| `NumRvLongHighMinusLowDeg` | 0.4 | 0.6 |
| `NumRvLongHighMinusLowP` | 0.81 | 0.70 |
| `NumRvLongRmsOverHistorical` | 0.43 | 0.42 |
| `NumRvLongRmsS` | 73 | 72 |
| `NumRvLongSlopeDegPerMag` | -0.7 | -0.9 |
| `NumRvLongSlopeP` | 0.79 | 0.75 |
| `NumRvLongStateDetectableDeg` | 4.3 | 4.2 |
| `NumRvLongStateHighErrDeg` | 0.8 | 0.7 |
| `NumRvLongStateHighMeanDeg` | -0.4 | -0.1 |
| `NumRvLongStateHighN` | 10 | 11 |
| `NumRvLongStateIntermediateErrDeg` | 2.0 | 2.6 |
| `NumRvLongStateIntermediateN` | 5 | 4 |
| `NumRvOcBudgetChiPerNG` | 1.06 | 1.04 |
| `NumRvOcBudgetChiPerNR` | 0.37 | 0.36 |
| `NumRvOcBudgetRmsRS` | 43 | 42 |
| `NumRvOcLimitRatioBand` | 0.90 | 0.89 |
| `NumRvOcLimitRatioBandNoTwoZeroTwoFour` | 4.57 | 4.56 |
| `NumRvOcLimitRatioBandPlusEra` | 2.26 | 2.24 |
| `NumRvOcLimitRatioNightNoTwoZeroTwoFour` | 3.81 | 3.79 |
| `NumRvOcMagcCycleScatterIS` | 46 | 44 |
| `NumRvOcMagcCycleScatterRS` | 47 | 46 |
| `NumRvOcMagcScatterBandConstIErrS` | 13 | 12 |
| `NumRvOcMagcScatterBandConstIS` | 5 | 7 |
| `NumRvOcMagcScatterBandConstRS` | -23 | -24 |
| `NumRvOcMagcScatterBandGMinusIS` | -54 | -55 |
| `NumRvOcMagcScatterBandLimitNightboot` | 2.7 \times 10^{-9} | 2.8 \times 10^{-9} |
| `NumRvOcMagcScatterBandPdotNsig` | 2.7 | 2.8 |
| `NumRvOcMagcScatterBandPeriodErrD` | 5.9 \times 10^{-8} | 5.8 \times 10^{-8} |
| `NumRvOcMagcScatterBandPlusEraEraOffsetErrS` | 165 | 161 |
| `NumRvOcMagcScatterBandPlusEraEraOffsetS` | 150 | 144 |
| `NumRvOcMagcScatterBandPlusEraLimitBudget` | 6.4 \times 10^{-9} | 6.3 \times 10^{-9} |
| `NumRvOcMagcScatterBandPlusEraLimitScatter` | 6.4 \times 10^{-9} | 6.3 \times 10^{-9} |
| `NumRvOcMagcScatterNightChinu` | 0.66 | 0.68 |
| `NumRvOcMagcScatterNightChisq` | 8.5 | 8.9 |
| `NumRvOcMagcScatterNightPlusEraChinu` | 0.67 | 0.69 |
| `NumRvOcMagcScatterNightPlusEraChisq` | 8.0 | 8.3 |
| `NumRvOcMagcScatterNightPlusEraEraOffsetErrS` | 136 | 134 |
| `NumRvOcMagcScatterNightPlusEraEraOffsetS` | 123 | 126 |
| `NumRvOcMagcScatterPooledLimitBudget` | 3.0 \times 10^{-9} | 2.9 \times 10^{-9} |
| `NumRvOcVOneCycleScatterGS` | 129 | 128 |
| `NumRvOcVOnePubBandChisq` | 20.5 | 20.3 |
| `NumRvOcVOnePubBandConstGS` | -66 | -65 |
| `NumRvOcVOnePubBandGMinusIS` | -88 | -86 |
| `NumRvOcVOnePubBandPdotNsig` | 2.0 | 2.1 |
| `NumRvOcVOnePubBandPeriodErrD` | 8.4 \times 10^{-8} | 8.3 \times 10^{-8} |
| `NumRvOcVOnePubBandPlusEraChisq` | 19.6 | 19.4 |
| `NumRvOcVOnePubBandPlusEraEraOffsetErrS` | 265 | 264 |
| `NumRvOcVOnePubBandPlusEraEraOffsetS` | 313 | 305 |
| `NumRvOcVOnePubBandPlusEraLimitScatter` | 10.0 \times 10^{-9} | 9.9 \times 10^{-9} |
| `NumRvOcVOnePubBandPlusEraRmsS` | 67 | 66 |
| `NumRvOcVOnePubNightPlusEraChinu` | 0.67 | 0.68 |
| `NumRvOcVOnePubNightPlusEraEraOffsetErrS` | 262 | 263 |
| `NumRvOcVOnePubNightPlusEraEraOffsetS` | 309 | 300 |
| `NumRvOcVOnePubNightPlusEraPdot` | 3.7 \times 10^{-9} | 3.6 \times 10^{-9} |
| `NumRvOcVOnePubPooledChinu` | 0.85 | 0.83 |
| `NumRvOcVOnePubPooledChisq` | 28.0 | 27.5 |
| `NumRvOcVOnePubPooledPdot` | 1.1 \times 10^{-9} | 1.2 \times 10^{-9} |
| `NumRvOcVOneScatterBandConstGErrS` | 33 | 32 |
| `NumRvOcVOneScatterBandConstGS` | -53 | -52 |
| `NumRvOcVOneScatterBandConstRErrS` | 18 | 17 |
| `NumRvOcVOneScatterBandGMinusIS` | -83 | -82 |
| `NumRvOcVOneScatterBandPdot` | 8.7 \times 10^{-10} | 8.8 \times 10^{-10} |
| `NumRvOcVOneScatterBandPlusEraEraOffsetS` | 221 | 218 |
| `NumRvOcVOneScatterBandPlusEraLimitBudget` | 8.2 \times 10^{-9} | 8.1 \times 10^{-9} |
| `NumRvOcVOneScatterBandPlusEraLimitScatter` | 8.2 \times 10^{-9} | 8.1 \times 10^{-9} |
| `NumRvOcVOneScatterBandPlusEraPdot` | 2.7 \times 10^{-9} | 2.6 \times 10^{-9} |
| `NumRvOcVOneScatterBandPlusEraPdotNsig` | 1.5 | 1.4 |
| `NumRvOcVOneScatterNightChinu` | 0.64 | 0.65 |
| `NumRvOcVOneScatterNightChisq` | 8.4 | 8.5 |
| `NumRvOcVOneScatterNightLimitScatter` | 2.7 \times 10^{-9} | 2.8 \times 10^{-9} |
| `NumRvOcVOneScatterNightPdot` | 8.7 \times 10^{-10} | 8.8 \times 10^{-10} |
| `NumRvOcVOneScatterNightPlusEraChinu` | 0.62 | 0.63 |
| `NumRvOcVOneScatterNightPlusEraChisq` | 7.4 | 7.5 |
| `NumRvOcVOneScatterNightPlusEraEraOffsetErrS` | 162 | 163 |
| `NumRvOcVOneScatterNightPlusEraEraOffsetS` | 201 | 199 |
| `NumRvOcVOneScatterPooledPdot` | 6.5 \times 10^{-10} | 6.7 \times 10^{-10} |
| `NumRvOcVOneScatterPooledRmsS` | 81 | 80 |
| `NumRvOcVOnecorrScatterBandConstGS` | -38 | -37 |
| `NumRvOcVOnecorrScatterBandLimitNightboot` | 4.3 \times 10^{-9} | 4.4 \times 10^{-9} |
| `NumRvOcVOnecorrScatterBandPdot` | 8.8 \times 10^{-10} | 8.6 \times 10^{-10} |
| `NumRvOcVOnecorrScatterBandPlusEraEraOffsetErrS` | 207 | 208 |
| `NumRvOcVOnecorrScatterBandPlusEraEraOffsetS` | 207 | 190 |
| `NumRvOcVOnecorrScatterBandPlusEraLimitBudget` | 8.0 \times 10^{-9} | 7.9 \times 10^{-9} |
| `NumRvOcVOnecorrScatterBandPlusEraLimitScatter` | 8.0 \times 10^{-9} | 7.9 \times 10^{-9} |
| `NumRvOcVOnecorrScatterBandPlusEraPdot` | 2.5 \times 10^{-9} | 2.4 \times 10^{-9} |
| `NumRvOcVOnecorrScatterBandPlusEraPdotNsig` | 1.4 | 1.3 |
| `NumRvOcVOnecorrScatterNightChinu` | 0.68 | 0.69 |
| `NumRvOcVOnecorrScatterNightChisq` | 8.8 | 8.9 |
| `NumRvOcVOnecorrScatterNightPdot` | 8.8 \times 10^{-10} | 8.6 \times 10^{-10} |
| `NumRvOcVOnecorrScatterNightPdotNsig` | 1.4 | 1.3 |
| `NumRvOcVOnecorrScatterNightPlusEraChinu` | 0.66 | 0.68 |
| `NumRvOcVOnecorrScatterNightPlusEraChisq` | 7.9 | 8.2 |
| `NumRvOcVOnecorrScatterNightPlusEraEraOffsetErrS` | 165 | 169 |
| `NumRvOcVOnecorrScatterNightPlusEraEraOffsetS` | 193 | 173 |
| `NumRvOcVOnecorrScatterNightPlusEraLimitBudget` | 7.8 \times 10^{-9} | 7.6 \times 10^{-9} |
| `NumRvOcVOnecorrScatterNightPlusEraLimitScatter` | 6.8 \times 10^{-9} | 6.7 \times 10^{-9} |
| `NumRvOcVOnecorrScatterNightPlusEraPdot` | 2.4 \times 10^{-9} | 2.2 \times 10^{-9} |
| `NumRvOcVOnecorrScatterNightPlusEraPdotNsig` | 1.7 | 1.5 |
| `NumRvOcVOnecorrScatterPooledPdot` | 7.3 \times 10^{-10} | 7.1 \times 10^{-10} |
| `NumRvOcVTwoCycleScatterRS` | 73 | 69 |
| `NumRvOcVTwoScatterBandConstGS` | -42 | -41 |
| `NumRvOcVTwoScatterBandConstIS` | 12 | 13 |
| `NumRvOcVTwoScatterBandConstRErrS` | 19 | 18 |
| `NumRvOcVTwoScatterBandConstRS` | 1 | -1 |
| `NumRvOcVTwoScatterBandPdot` | 9.2 \times 10^{-10} | 9.0 \times 10^{-10} |
| `NumRvOcVTwoScatterBandPeriodD` | 0.07908914 | 0.07908915 |
| `NumRvOcVTwoScatterBandPlusEraEraOffsetErrS` | 167 | 164 |
| `NumRvOcVTwoScatterBandPlusEraEraOffsetS` | 288 | 308 |
| `NumRvOcVTwoScatterBandPlusEraLimitBudget` | 7.1 \times 10^{-9} | 7.2 \times 10^{-9} |
| `NumRvOcVTwoScatterBandPlusEraLimitScatter` | 7.1 \times 10^{-9} | 7.2 \times 10^{-9} |
| `NumRvOcVTwoScatterBandPlusEraPdot` | 3.0 \times 10^{-9} | 3.2 \times 10^{-9} |
| `NumRvOcVTwoScatterBandPlusEraPdotNsig` | 2.2 | 2.4 |
| `NumRvOcVTwoScatterBandPlusEraRmsS` | 53 | 52 |
| `NumRvOcVTwoScatterBandRmsS` | 53 | 52 |
| `NumRvOcVTwoScatterNightChinu` | 0.99 | 1.05 |
| `NumRvOcVTwoScatterNightChisq` | 12.9 | 13.6 |
| `NumRvOcVTwoScatterNightLimitBudget` | 2.5 \times 10^{-9} | 2.4 \times 10^{-9} |
| `NumRvOcVTwoScatterNightLimitScatter` | 2.4 \times 10^{-9} | 2.5 \times 10^{-9} |
| `NumRvOcVTwoScatterNightPdot` | 9.2 \times 10^{-10} | 9.0 \times 10^{-10} |
| `NumRvOcVTwoScatterNightPdotNsig` | 1.8 | 1.7 |
| `NumRvOcVTwoScatterNightPlusEraChinu` | 0.92 | 0.95 |
| `NumRvOcVTwoScatterNightPlusEraChisq` | 11.1 | 11.4 |
| `NumRvOcVTwoScatterNightPlusEraEraOffsetErrS` | 175 | 176 |
| `NumRvOcVTwoScatterNightPlusEraEraOffsetS` | 243 | 268 |
| `NumRvOcVTwoScatterNightPlusEraLimitBudget` | 7.0 \times 10^{-9} | 7.1 \times 10^{-9} |
| `NumRvOcVTwoScatterNightPlusEraLimitScatter` | 6.8 \times 10^{-9} | 7.0 \times 10^{-9} |
| `NumRvOcVTwoScatterNightPlusEraPdot` | 2.7 \times 10^{-9} | 2.9 \times 10^{-9} |
| `NumRvOcVTwoScatterNightPlusEraPdotNsig` | 2.0 | 2.1 |
| `NumRvOcVTwoScatterNightRmsS` | 30 | 31 |
| `NumRvOcVTwoScatterPooledLimitBudget` | 2.5 \times 10^{-9} | 2.4 \times 10^{-9} |
| `NumRvOcVTwoScatterPooledPdot` | 6.4 \times 10^{-10} | 6.3 \times 10^{-10} |
| `NumRvOcVTwoScatterPooledRmsS` | 61 | 60 |
| `NumRvOffMagcGiEraAllChiBoot` | 0.83 | 1.00 |
| `NumRvOffMagcGiNightMeanS` | -88 | -89 |
| `NumRvOffMagcGrNightMeanS` | -25 | -24 |
| `NumRvOffMagcGrNightPPerm` | 0.25 | 0.27 |
| `NumRvOffMagcGrNightPPermBonf` | 0.74 | 0.81 |
| `NumRvOffMagcRiEraAllChiBoot` | 0.90 | 0.74 |
| `NumRvOffMagcRiEraAllChiPub` | 0.31 | 0.23 |
| `NumRvOffMagcRiEraAllMeanS` | -11 | -14 |
| `NumRvOffMagcRiEraAllPCluster` | 0.76 | 0.70 |
| `NumRvOffMagcRiEraAllPPerm` | 0.66 | 0.52 |
| `NumRvOffMagcRiEraAllSeS` | 24 | 22 |
| `NumRvOffMagcRiEraAllT` | -0.5 | -0.7 |
| `NumRvOffMagcRiEraSevenMeanS` | 71 | 55 |
| `NumRvOffMagcRiEraSevenSeS` | 65 | 50 |
| `NumRvOffMagcRiNightMeanS` | -25 | -27 |
| `NumRvOffMagcRiNightPPerm` | 0.23 | 0.20 |
| `NumRvOffMagcRiNightPPermBonf` | 0.70 | 0.59 |
| `NumRvOffMagcRiNightSeS` | 20 | 19 |
| `NumRvOffVOneGiEraAllChiBoot` | 2.50 | 2.49 |
| `NumRvOffVOneGiEraSevenMeanS` | -255 | -254 |
| `NumRvOffVOneGiEraSevenSeS` | 138 | 143 |
| `NumRvOffVOneGiNightMeanS` | -123 | -121 |
| `NumRvOffVOneGrEraAllChiBoot` | 1.14 | 1.17 |
| `NumRvOffVOneGrEraAllMeanS` | -18 | -15 |
| `NumRvOffVOneGrEraAllPCluster` | 0.49 | 0.56 |
| `NumRvOffVOneGrEraAllPPerm` | 0.54 | 0.60 |
| `NumRvOffVOneGrEraAllPWilcoxon` | 0.54 | 0.59 |
| `NumRvOffVOneGrEraAllSeS` | 28 | 29 |
| `NumRvOffVOneGrEraAllT` | -0.6 | -0.5 |
| `NumRvOffVOneGrEraSevenMeanS` | 26 | 54 |
| `NumRvOffVOneGrNightMeanS` | -39 | -37 |
| `NumRvOffVOneGrNightPPermBonf` | 0.86 | 0.88 |
| `NumRvOffVTwoGrNightMeanS` | -32 | -30 |
| `NumRvOffVTwoGrNightPPerm` | 0.027 | 0.031 |
| `NumRvOffVTwoGrNightPPermBonf` | 0.082 | 0.094 |
| `NumRvOffVTwoGrNightSeS` | 13 | 14 |
| `NumRvOffVTwoRiEraAllChiBoot` | 0.50 | 0.41 |
| `NumRvOffVTwoRiEraAllChiPub` | 0.42 | 0.30 |
| `NumRvOffVTwoRiEraAllMeanS` | -9 | -14 |
| `NumRvOffVTwoRiEraAllPCluster` | 0.77 | 0.61 |
| `NumRvOffVTwoRiEraAllPPerm` | 0.74 | 0.56 |
| `NumRvOffVTwoRiEraAllSeS` | 26 | 23 |
| `NumRvOffVTwoRiEraAllT` | -0.3 | -0.6 |
| `NumRvOffVTwoRiEraSevenMeanS` | 109 | 87 |
| `NumRvOffVTwoRiEraSevenSeS` | 54 | 33 |
| `NumRvOffVTwoRiNightMeanS` | -27 | -29 |
| `NumRvOffVTwoRiNightPPerm` | 0.29 | 0.24 |
| `NumRvOffVTwoRiNightPPermBonf` | 0.86 | 0.71 |
| `NumRvOffVTwoRiNightSeS` | 24 | 23 |
| `NumRvOffVTwoRiOverPredicted` | 0.5 | 0.8 |
| `NumRvPdotAsyncOverLimitEra` | 4.7 | 4.8 |
| `NumRvPdotLimitEraOverLibration` | 229 | 228 |
| `NumRvShContourMaxFullMag` | 0.47 | 0.45 |
| `NumRvShContourMaxMmag` | 237 | 226 |
| `NumRvShExcludingOneTwoFiveMmag` | 2 | 3 |
| `NumRvShExcludingOneZeroZeroMmag` | 2 | 3 |
| `NumRvShExcludingPeak` | 2 | 3 |
| `NumRvShNotExcludingPeak` | 16 | 15 |
| `NumRvSigmaInjVOneEraSevenGS` | 147 | 151 |
| `NumRvSigmaInjVOneEraSevenIS` | 46 | 47 |
| `NumRvSigmaInjVOneEraSevenRS` | 59 | 63 |
| `NumRvSigmaInjVTwoEraSevenGS` | 113 | 112 |
| `NumRvSigmaInjVTwoEraSevenIS` | 53 | 50 |
| `NumRvSigmaInjVTwoEraSevenRS` | 56 | 59 |
| `NumRvVOneGridStepEraSevenMaxS` | 8.1 | 8.0 |
| `NumRvVOneValleyMaxS` | 488 | 494 |

### New in C

`NumRvAnnulusInnerArcsec`, `NumRvAnnulusOuterArcsec`, `NumRvApertureRadiusArcsec`, `NumRvBandPairTrials`, `NumRvCapCheckCapped`, `NumRvCapCompAll`, `NumRvCapCompCapped`, `NumRvCapCompCappedPct`, `NumRvCapHgAdu`, `NumRvCapHgVetoAdu`, `NumRvCapSeries`, `NumRvCapStlmiFoursevenRComp`, `NumRvCapStlmiFoursevenRTarget`, `NumRvCapStlmiFoursevenYComp`, `NumRvCapStlmiFoursevenYTarget`, `NumRvCapStlmiFoursevenZComp`, `NumRvCapStlmiFoursevenZTarget`, `NumRvCapStlmiSevenGComp`, `NumRvCapStlmiSevenGTarget`, `NumRvCapStlmiSevenIComp`, `NumRvCapStlmiSevenITarget`, `NumRvCapStlmiSevenRComp`, `NumRvCapStlmiSevenRTarget`, `NumRvCapTargetAll`, `NumRvCapTargetCapped`, `NumRvClockBarS`, `NumRvClockEclipseCycles`, `NumRvClockEclipseEphErrS`, `NumRvClockEclipseEphS`, `NumRvClockEclipseErrS`, `NumRvClockEclipseGaiaSigma`, `NumRvClockEclipseOcS`, `NumRvClockEraAErrS`, `NumRvClockEraAOcS`, `NumRvClockEraASeries`, `NumRvClockEraBErrS`, `NumRvClockEraBOcS`, `NumRvClockEraBSeries`, `NumRvClockEraCErrS`, `NumRvClockEraCOcS`, `NumRvClockEraCSeries`, `NumRvClockEraDVOneNights`, `NumRvClockEraDVOneOffsetErrS`, `NumRvClockEraDVOneOffsetS`, `NumRvClockEraDVOneVerifiedNights`, `NumRvClockEraDVTwoNights`, `NumRvClockEraDVTwoOffsetErrS`, `NumRvClockEraDVTwoOffsetS`, `NumRvClockEraDVTwoVerifiedNights`, `NumRvClockEraEErrS`, `NumRvClockEraEOcS`, `NumRvClockEraESeries`, `NumRvClockEraFErrS`, `NumRvClockEraFOcS`, `NumRvClockEraFSeries`, `NumRvClockErasPass`, `NumRvClockErasTested`, `NumRvClockEuumaOffsetS`, `NumRvDetFastCapAdu`, `NumRvDetFastCapDevPct`, `NumRvDetFastGain`, `NumRvDetFastGainErr`, `NumRvDetFastLinSlopeErrPct`, `NumRvDetFastLinSlopePct`, `NumRvDetHgCapAdu`, `NumRvDetHgCapDevPct`, `NumRvDetHgGain`, `NumRvDetHgGainErr`, `NumRvDetHgInflFactor`, `NumRvDetHgInflMax`, `NumRvDetHgInflMeasMax`, `NumRvDetHgInflMeasMin`, `NumRvDetHgInflMin`, `NumRvDetHgLinSlopeErrPct`, `NumRvDetHgLinSlopePct`, `NumRvDetHgRnE`, `NumRvDetHgRnErrE`, `NumRvDetMeasuredGains`, `NumRvDetModezeroCapAdu`, `NumRvDetModezeroCapDevPct`, `NumRvDetModezeroGain`, `NumRvDetModezeroGainErr`, `NumRvDetModezeroInflFactor`, `NumRvDetModezeroInflMax`, `NumRvDetModezeroInflMeasMax`, `NumRvDetModezeroInflMeasMin`, `NumRvDetModezeroInflMin`, `NumRvDetModezeroLinSlopeErrPct`, `NumRvDetModezeroLinSlopePct`, `NumRvDetModezeroRnE`, `NumRvDetModezeroRnErrE`, `NumRvDetNominalGain`, `NumRvDetOnemhzCapAdu`, `NumRvDetOnemhzCapDevPct`, `NumRvDetOnemhzGain`, `NumRvDetOnemhzGainErr`, `NumRvDetOnemhzInflFactor`, `NumRvDetOnemhzInflMax`, `NumRvDetOnemhzInflMeasMax`, `NumRvDetOnemhzInflMeasMin`, `NumRvDetOnemhzInflMin`, `NumRvDetOnemhzLinSlopeErrPct`, `NumRvDetOnemhzLinSlopePct`, `NumRvDetOnemhzRnE`, `NumRvDetOnemhzRnErrE`, `NumRvDetStackproCapAdu`, `NumRvDetStackproCapDevPct`, `NumRvDetStackproGain`, `NumRvDetStackproGainErr`, `NumRvDetStackproLinSlopeErrPct`, `NumRvDetStackproLinSlopePct`, `NumRvDetStackproRnE`, `NumRvDetStackproRnErrE`, `NumRvDetStackproRnRatio`, `NumRvEdgeqVOneEraFoursevenZChinuMax`, `NumRvEdgeqVOneEraFoursevenZChinuMed`, `NumRvEdgeqVOneEraFoursevenZChinuMin`, `NumRvEdgeqVOneEraFoursevenZDofMax`, `NumRvEdgeqVOneEraFoursevenZDofMin`, `NumRvEdgeqVOneEraFoursevenZN`, `NumRvEdgeqVOneEraSevenGChinuMax`, `NumRvEdgeqVOneEraSevenGChinuMed`, `NumRvEdgeqVOneEraSevenGChinuMin`, `NumRvEdgeqVOneEraSevenGDofMax`, `NumRvEdgeqVOneEraSevenGDofMin`, `NumRvEdgeqVOneEraSevenGN`, `NumRvEdgeqVOneEraSevenIChinuMax`, `NumRvEdgeqVOneEraSevenIChinuMed`, `NumRvEdgeqVOneEraSevenIChinuMin`, `NumRvEdgeqVOneEraSevenIDofMax`, `NumRvEdgeqVOneEraSevenIDofMin`, `NumRvEdgeqVOneEraSevenIN`, `NumRvEdgeqVOneEraSevenRChinuMax`, `NumRvEdgeqVOneEraSevenRChinuMed`, `NumRvEdgeqVOneEraSevenRChinuMin`, `NumRvEdgeqVOneEraSevenRDofMax`, `NumRvEdgeqVOneEraSevenRDofMin`, `NumRvEdgeqVOneEraSevenRN`, `NumRvEdgeqVOneEraSevensixGChinuMax`, `NumRvEdgeqVOneEraSevensixGChinuMed`, `NumRvEdgeqVOneEraSevensixGChinuMin`, `NumRvEdgeqVOneEraSevensixGDofMax`, `NumRvEdgeqVOneEraSevensixGDofMin`, `NumRvEdgeqVOneEraSevensixGN`, `NumRvEdgeqVOneEraSevensixIChinuMax`, `NumRvEdgeqVOneEraSevensixIChinuMed`, `NumRvEdgeqVOneEraSevensixIChinuMin`, `NumRvEdgeqVOneEraSevensixIDofMax`, `NumRvEdgeqVOneEraSevensixIDofMin`, `NumRvEdgeqVOneEraSevensixIN`, `NumRvEdgeqVOneEraSevensixRChinuMax`, `NumRvEdgeqVOneEraSevensixRChinuMed`, `NumRvEdgeqVOneEraSevensixRChinuMin`, `NumRvEdgeqVOneEraSevensixRDofMax`, `NumRvEdgeqVOneEraSevensixRDofMin`, `NumRvEdgeqVOneEraSevensixRN`, `NumRvEdgeqVTwoEraFoursevenZChinuMax`, `NumRvEdgeqVTwoEraFoursevenZChinuMed`, `NumRvEdgeqVTwoEraFoursevenZChinuMin`, `NumRvEdgeqVTwoEraFoursevenZDofMax`, `NumRvEdgeqVTwoEraFoursevenZDofMin`, `NumRvEdgeqVTwoEraFoursevenZN`, `NumRvEdgeqVTwoEraSevenGChinuMax`, `NumRvEdgeqVTwoEraSevenGChinuMed`, `NumRvEdgeqVTwoEraSevenGChinuMin`, `NumRvEdgeqVTwoEraSevenGDofMax`, `NumRvEdgeqVTwoEraSevenGDofMin`, `NumRvEdgeqVTwoEraSevenGN`, `NumRvEdgeqVTwoEraSevenIChinuMax`, `NumRvEdgeqVTwoEraSevenIChinuMed`, `NumRvEdgeqVTwoEraSevenIChinuMin`, `NumRvEdgeqVTwoEraSevenIDofMax`, `NumRvEdgeqVTwoEraSevenIDofMin`, `NumRvEdgeqVTwoEraSevenIN`, `NumRvEdgeqVTwoEraSevenRChinuMax`, `NumRvEdgeqVTwoEraSevenRChinuMed`, `NumRvEdgeqVTwoEraSevenRChinuMin`, `NumRvEdgeqVTwoEraSevenRDofMax`, `NumRvEdgeqVTwoEraSevenRDofMin`, `NumRvEdgeqVTwoEraSevenRN`, `NumRvEdgeqVTwoEraSevensixGChinuMax`, `NumRvEdgeqVTwoEraSevensixGChinuMed`, `NumRvEdgeqVTwoEraSevensixGChinuMin`, `NumRvEdgeqVTwoEraSevensixGDofMax`, `NumRvEdgeqVTwoEraSevensixGDofMin`, `NumRvEdgeqVTwoEraSevensixGN`, `NumRvEdgeqVTwoEraSevensixIChinuMax`, `NumRvEdgeqVTwoEraSevensixIChinuMed`, `NumRvEdgeqVTwoEraSevensixIChinuMin`, `NumRvEdgeqVTwoEraSevensixIDofMax`, `NumRvEdgeqVTwoEraSevensixIDofMin`, `NumRvEdgeqVTwoEraSevensixIN`, `NumRvEdgeqVTwoEraSevensixRChinuMax`, `NumRvEdgeqVTwoEraSevensixRChinuMed`, `NumRvEdgeqVTwoEraSevensixRChinuMin`, `NumRvEdgeqVTwoEraSevensixRDofMax`, `NumRvEdgeqVTwoEraSevensixRDofMin`, `NumRvEdgeqVTwoEraSevensixRN`, `NumRvFitVOneSpanCadences`, `NumRvFramesLocalReduced`, `NumRvFramesServerReduced`, `NumRvLinCriterionPct`, `NumRvLocalFlatAgeMaxD`, `NumRvLocalFlatAgeMinD`, `NumRvMechMultiState`, `NumRvMechMultiStateList`, `NumRvMechSeries`, `NumRvMechStepMaxMmag`, `NumRvMechStepMedianMmag`, `NumRvMechStlmiEraSevenSixStates`, `NumRvMechStlmiGStepMmag`, `NumRvMechStlmiLateNights`, `NumRvRadialEraFoursevenMaxMmag`, `NumRvRadialEraSevenMaxMmag`, `NumRvRadialEraSevensixMaxMmag`, `NumRvRadialEraSeventwoMaxMmag`, `NumRvRadialMaxMmag`, `NumRvRadialMedianMmag`, `NumRvRadialSignificant`, `NumRvRampBelowBar`, `NumRvRampBlocks`, `NumRvRampClipSigma`, `NumRvRampConsistentWithBar`, `NumRvRampEraFoursevenBlocks`, `NumRvRampEraFoursevenMaxPct`, `NumRvRampEraFoursevenMinPct`, `NumRvRampEraSevenBlocks`, `NumRvRampEraSevenMaxPct`, `NumRvRampEraSevenMinPct`, `NumRvRampEraSevensixBlocks`, `NumRvRampEraSevensixMaxPct`, `NumRvRampEraSevensixMinPct`, `NumRvRampEraSeventwoBlocks`, `NumRvRampEraSeventwoMaxPct`, `NumRvRampEraSeventwoMinPct`, `NumRvRampFlatBarPct`, `NumRvRampMaxAbsErrPct`, `NumRvRampMaxAbsPct`, `NumRvRampMaxBlock`, `NumRvRampMedianAbsPct`, `NumRvRampSignificant`, `NumRvRampStlmiFourSevenZErrPct`, `NumRvRampStlmiFourSevenZPct`, `NumRvRampStlmiFourSevenZRadialMmag`, `NumRvRampStlmiSevenGErrPct`, `NumRvRampStlmiSevenGPct`, `NumRvRampStlmiSevenGRadialMmag`, `NumRvRampStlmiSevenIErrPct`, `NumRvRampStlmiSevenIPct`, `NumRvRampStlmiSevenIRadialMmag`, `NumRvRampStlmiSevenRErrPct`, `NumRvRampStlmiSevenRPct`, `NumRvRampStlmiSevenRRadialMmag`, `NumRvRampStlmiSevenSixGErrPct`, `NumRvRampStlmiSevenSixGPct`, `NumRvRampStlmiSevenSixGRadialMmag`, `NumRvRampStlmiSevenSixIErrPct`, `NumRvRampStlmiSevenSixIPct`, `NumRvRampStlmiSevenSixIRadialMmag`, `NumRvRampStlmiSevenSixRErrPct`, `NumRvRampStlmiSevenSixRPct`, `NumRvRampStlmiSevenSixRRadialMmag`, `NumRvRawredEraSevenCheckRmsMmag`, `NumRvRawredEraSevenFrames`, `NumRvRawredEraSevenOffsetMmag`, `NumRvRawredEraSevenSigmaChkMmag`, `NumRvRawredEraSevenTargetRmsMmag`, `NumRvRawredEraSevensixCheckRmsMmag`, `NumRvRawredEraSevensixFlatAgeD`, `NumRvRawredEraSevensixFrames`, `NumRvRawredEraSevensixOffsetMmag`, `NumRvRawredEraSevensixSigmaChkMmag`, `NumRvRawredEraSevensixTargetRmsMmag`, `NumRvRawredEraSeventwoCheckRmsMmag`, `NumRvRawredEraSeventwoFlatAgeD`, `NumRvRawredEraSeventwoFrames`, `NumRvRawredEraSeventwoOffsetMmag`, `NumRvRawredEraSeventwoSigmaChkMmag`, `NumRvRawredEraSeventwoTargetRmsMmag`, `NumRvRawredNights`, `NumRvRawredNightsAgree`, `NumRvRollingAcFourZeroFourZeroS`, `NumRvRollingAsiS`, `NumRvRollingQhyS`, `NumRvShutterIkonS`

Reason for every B->C change of an `rv_result` macro: revision statistics recomputed on the re-run chain (the photometry now withholds measurements above the S2 1% linearity cap; the rv stages had last run before the F-10 rebuild of CV-S9, so they also absorb that rebuild)

