# ADNI conditional pathology-information and prospective validation protocol

Status: protocol to freeze before participant-level execution

## Primary objective

The primary experiment asks whether accessible candidate measurements contain information about
continuous PET-defined Alzheimer pathology beyond contemporaneous cognitive/clinical information
and fixed inherited background. This is the direct analogue of testing whether wearable signals
contain virus information beyond activity.

- Candidate `X` blocks: plasma biomarkers, structural MRI, and (in a reciprocal analysis) genetics.
- Pathology `Y`: amyloid PET Centiloid or tau PET temporal meta-ROI SUVR.
- Context `Z`: cognition, demographics, time and technical variables; genetics is a time-invariant
  part of `Z` unless genetics is the tested `X`.

The prospective future-ADAS analysis below is downstream validation rather than the only target.

## Downstream prospective objective

The downstream prospective experiment tests whether a multimodal amyloid/tau/neurodegeneration
(A/T/N) panel contains information about repeated cognition over the following 12--36 months
beyond information already available at the origin visit. A larger repeated-MRI analysis is
retained as a secondary prospective experiment. A prespecified genetic extension evaluates
inherited risk without reducing the A/T/N cohort to the smaller genotype-complete intersection.

For candidate biomarker `j`, the target is

`H0,j: X(j) at origin is conditionally independent of future ADAS13 given origin-time Z`.

The experiment is designed to distinguish incremental marker information from a same-visit
association or a biologically plausible direction.

## Prospective data choice

- Candidate markers: amyloid PET Centiloid, flortaucipir temporal meta-ROI SUVR, and
  intracranial-volume-adjusted bilateral hippocampal volume.
- Origin: tau PET date, with amyloid PET within six months and MRI/baseline ADAS13 within
  90 days.
- Outcome: repeated future ADAS-Cog13 nearest 12, 24, and 36 months within frozen ±3-month
  windows.
- Origin-time reference set: current ADAS13, current diagnosis, permitted prior cognitive
  history, age, sex, education, APOE4, intracranial volume, site, ADNI phase/protocol, MRI
  field strength, time since baseline, MRI-to-origin date gap, and exact forecast interval.
- Independent unit: participant. Every split and uncertainty calculation operates at the
  participant level.

Within the downstream prospective validation, this A/T/N panel directly compares distinct
biomarker domains against
an established biomedical ordering. Tau PET is expected to provide the strongest unique
prognostic information, MRI may contribute complementary neurodegeneration information, and
amyloid may contribute less once tau and baseline cognition are known. The repeated-MRI
12-month experiment remains secondary because it has stronger longitudinal coverage.

## Genetic extension

- Genetic candidates: APOE ε4 dosage and one externally weighted non-APOE AD PRS/PHS.
- Data locations: APOE under `Biospecimen Results`; external polygenic-hazard outputs under
  `Genetic -> Genotype results`; GWAS inputs, if required, under `Downloads -> Genetic Files`.
- The score definition, allele alignment, QC, ancestry population, and APOE-region exclusion are
  frozen before outcome analysis. SNPs and score thresholds are not selected in ADNI.
- Genetic reference set: origin ADAS13 and diagnosis, age, sex, education, site/phase, and ancestry
  principal components. APOE is not included in `Z` when APOE is the tested candidate.

Two estimands are reported separately:

1. Genetic information beyond baseline clinical state, without A/T/N in the reference set.
2. Genetic information beyond baseline clinical state and the full A/T/N panel.

APOE or PRS attenuation in the second analysis is not a failed experiment: it can indicate that
current measured pathology mediates or absorbs the score's prognostic information. Because the
genetic candidates are time invariant, their participant score is their residualised genetic value
times the equal-weighted mean future-ADAS residual for that participant. All model fitting and
uncertainty calculations remain participant-level.

## Composite settings

Define `R` as demographics, site/phase, timing, and ancestry principal components; `C` as current
diagnosis, current ADAS13, and permitted prior cognition; `G` as APOE and the frozen non-APOE
PRS/PHS; and `B=(A,T,N)` as current amyloid, tau, and MRI measurements. Genetics has no time index.

The main MRI experiment uses an oriented hippocampal-volume, entorhinal-thickness, and ventricular-
volume panel as `X`, repeated future ADAS13 as `Y`, and nested reference sets:

1. `Z=R+C`;
2. `Z=R+C+G` (main incremental-information setting);
3. `Z=R+C+G+A+T` (unique MRI/neurodegeneration setting).

The reciprocal experiment uses `X=G` first with `Z=R+C`, then with `Z=R+C+B`. A joint global test
uses `X=(B,G)` with `Z=R+C`. The three resulting claims are kept distinct: joint information,
biomarker information beyond genetics, and genetic information beyond measured biomarkers.

Time-invariant genetics is stored once per participant. It may be broadcast internally when forming
visit-level predictions, but its inferential contribution is formed once as the residualised genetic
value times the equal-weighted participant mean of future-outcome residuals; duplicated genotype
rows never increase the effective sample size.

## Frozen pairing rules

1. An origin record must contain all candidate A/T/N markers, the reference variables, and an
   origin ADAS13 value.
2. Search only forward in time for future ADAS13 assessments.
3. Select assessments closest to 12, 24, and 36 months within ±3-month windows.
4. Resolve equal-distance ties by the earlier future date.
5. Record the exact forecast interval and do not use future variables in preprocessing,
   feature selection, imputation, or nuisance tuning.
6. Retain all eligible origin-future pairs, but aggregate inferential contributions to one
   equal-weighted score per participant.

## Analysis hierarchy

1. Prospective global cluster-GCM test across amyloid, tau, and MRI markers.
2. Holm-adjusted marker-wise tests.
3. Unique-marker tests: test each biological domain after adding the other two domains to the
   reference set.
4. Subject-held-out predictive comparison of nested clinical, amyloid, MRI, and tau models,
   reporting paired change in RMSE and MAE with participant-level uncertainty.
5. Separate genetic comparisons: clinical versus clinical-plus-genetics, and clinical-plus-A/T/N
   versus clinical-plus-A/T/N-plus-genetics.

The global test answers whether any candidate marker adds residual information. Marker-wise
tests locate detectable signals. Unique-marker tests ask whether a region contributes
complementary information. Prediction changes assess practical incremental value. These are
related but non-interchangeable claims.

## Candidate discrimination for PET pathology

For amyloid PET and tau PET separately, report:

1. a global test across all frozen candidate measurements;
2. multiplicity-adjusted block tests for plasma, MRI, and genetics;
3. unique-block tests with the other candidate blocks added to `Z`;
4. Holm-adjusted marker-wise and unique-marker tests within each block;
5. participant-held-out continuous-PET RMSE/MAE gains and secondary PET-positivity AUC/Brier gains.

The primary plasma panel is p-tau217, Aβ42/40, GFAP, and NfL. The primary MRI panel is oriented
hippocampal, entorhinal, and ventricular FreeSurfer 7 measures. Genetics consists of APOE ε4 dosage
and one frozen externally weighted non-APOE PRS/PHS. For tau PET, amyloid burden is included in `Z`
so that the test concerns tau-specific information beyond amyloid stage.

Results are classified as total information, unique conditional information, residual direction,
and held-out discrimination gain. A global rejection alone does not establish that every candidate
is useful, and a plausible direction alone does not establish incremental information.

Candidate p-values are not ranked. Direct comparison uses a common-cohort cross-fitted normalised
residual covariance for each oriented marker, participant-bootstrap intervals for marker and
pairwise-difference estimates, unique-marker/block tests, and paired differences in held-out loss.
The maximal available cohort for each block is a sensitivity analysis, not a basis for ranking
unequal-cohort p-values.

The final biomarker-style output is a candidate-by-target signature matrix: amyloid PET, tau PET
conditional on amyloid, and future cognition. Plasma p-tau217 is expected to show broad A/T signal;
Aβ42/40 to be more amyloid-focused; NfL and MRI to be more related to neurodegeneration/cognition;
and inherited risk to attenuate after measured pathology enters `Z`. Deviations are reported and
audited rather than deleted.

An ADNI4 remote/digital extension is reserved for feasibility screening: Novoic Storyteller and
ECog-12 form accessible `X`, standard in-clinic cognition and genetics form `Z`, and plasma or PET
pathology forms `Y`. It is not promoted until continuous-date linkage and pathology overlap are
large enough; remote timepoints are not equated with in-clinic visit codes.

## Sensitivity analyses

- Alternative forecast windows fixed before execution.
- Complete-case versus a pre-specified origin-time imputation strategy.
- Baseline diagnosis strata, reported as heterogeneity rather than separate discoveries.
- Influence analysis at the participant-score level.
- Random-slope stress calibration matched to the retained visit distribution.

## Reporting requirements

- Cohort flow from available origins to complete eligible pairs.
- Participants, origin-future pairs, visits per participant, and forecast-interval summaries.
- Global and Holm-adjusted marker results.
- Unique-marker results.
- Out-of-subject RMSE and MAE changes with uncertainty.
- Residual diagnostics, fold balance, and participant influence.
- Explicit distinction between contemporaneous association, prospective incremental
  information, prediction gain, and clinical utility.

## Data governance

Participant-level execution must occur in an institutionally approved environment consistent
with the current ADNI data-use agreement. Participant data, credentials, and restricted
derived rows must not enter the public repository. Only disclosure-checked aggregate outputs,
configuration, provenance, and source-table hashes may be released.
