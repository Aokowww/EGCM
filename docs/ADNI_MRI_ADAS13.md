# Exploratory ADNI MRI-ADAS13 analysis

## Question and cohort

This analysis asks whether residual structural MRI measures are associated
with concurrent ADAS13 after adjustment for baseline and technical
covariates. It does not estimate a causal effect or a rate of cognitive
change.

ADAS13 assessment date is the anchor. For each assessment, the preparation
script selects the nearest quality-controlled FreeSurfer 7 MRI scan within 30
days. Ties are resolved using completeness, quality information, field
strength and image identifier, without reference to ADAS13.

The A/B cohort requires at least five complete observations per participant.

| Quantity | Value |
|---|---:|
| Complete date-matched observations | 8,628 |
| Participants before the repeat restriction | 2,355 |
| Observations after the repeat restriction | 4,406 |
| Participants after the repeat restriction | 647 |
| Median visits | 6 |
| Visit range | 5 to 13 |
| Median absolute MRI-ADAS gap | 1 day |
| 90th percentile absolute gap | 19 days |

The adjustment set contains time since baseline, baseline age, sex, education,
APOE4, intracranial volume, study protocol, site, MRI field strength and the
MRI-ADAS date gap. Concurrent diagnosis and other cognitive measurements are
excluded.

## Results

The subject-level global multiplier p-value was 0.0001 for both calculations.
This is the resolution limit with 9,999 multiplier draws, not a zero p-value.

| Design | Marker | Statistic | Holm p | Normalized residual covariance |
|---|---|---:|---:|---:|
| A | Ventricles | 7.850 | 4.16e-15 | 0.307 |
| A | Hippocampus | -12.979 | 3.20e-38 | -0.466 |
| A | Entorhinal | -13.403 | 1.75e-40 | -0.419 |
| B | Ventricles | 9.007 | 6.35e-19 | 0.455 |
| B | Hippocampus | -8.462 | 5.25e-17 | -0.256 |
| B | Entorhinal | -8.367 | 5.94e-17 | -0.236 |

Higher residual ventricular volume was associated with higher ADAS13.
Residual hippocampal and entorhinal volumes had negative associations with
ADAS13. Removing any one participant preserved all six directions.

Design A subject folds contained 123 to 136 participants. Its maximum absolute
standardized subject scores ranged from 6.46 to 7.47. Design B's corresponding
range was 7.14 to 14.05; the ventricular result has a clear influence warning.

## Interpretation

Design A passed the core-null simulation and five of six stress settings. It
failed the complete stress gate because the random-slope setting did not meet
the frozen criterion. The ADNI result is therefore exploratory.

Design B can be computed in this denser repeated-visit cohort, but its earlier
null rejection rates were 0.29 and 0.39. Its p-values are reported as
sensitivity results, not calibrated evidence for a subject-conditional null.

The estimates concern same-occasion conditional association. They do not show
that MRI differences preceded ADAS13 differences.

## Reproduction

- Configuration: `configs/adni_mri_adas13_exploratory.yaml`
- Cohort builder: `code/adni/prepare_adni_mri_adas13.R`
- Analysis: `code/adni/run_adni_redesign.R`
- Diagnostic summary: `code/adni/summarize_adni_mri_adas13.R`
- Local analytic CSV SHA-256:
  `bb426c9b0df2922a310319e3e79b1288e1f0088ed3848be0d3683af3ad88a400`
- Local result RDS SHA-256:
  `56765ca2d903312fd136be39366e767c00f888b52afe3c887de1fab1fb62fa2e`

Neither local file is distributed.
