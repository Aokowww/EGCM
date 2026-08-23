# Exploratory ADNI plasma-Centiloid analysis

## Question

This analysis asks whether four plasma measurements contain information about
continuous amyloid PET Centiloids beyond cognition, demographics, APOE, and
acquisition timing. It is an incremental-information analysis, not a causal
effect estimate or a clinical validation study.

The plasma panel was fixed as p-tau217, A-beta42/40, NfL, and GFAP. Amyloid PET
was the anchor. The nearest complete plasma and cognitive measurements within
90 days were selected using deterministic date rules. All folds and uncertainty
calculations used participants, not rows, as the independent unit.

## Cohort

The complete analysis table contained 758 PET-anchored observations from 664
participants. Ninety-four participants had repeated observations. The median
absolute plasma-PET gap was seven days, the 90th percentile was 42 days, and
the maximum permitted gap was 90 days.

Participant rows and identifiers are not distributed. The public files contain
only aggregate counts, test summaries, and prediction metrics.

## Results

The four-marker global multiplier test returned `p = 0.0001`, the resolution
of 9,999 draws. Marker-unique tests added the other three plasma measurements
to the conditioning set:

| Marker | Normalised residual covariance | Holm p-value | Interpretation |
|---|---:|---:|---|
| p-tau217 | 0.5047 | 9.19e-16 | Unique incremental information detected |
| A-beta42/40 | -0.2640 | 8.25e-12 | Unique incremental information detected |
| NfL | -0.0979 | 0.1571 | No unique incremental information detected |
| GFAP | 0.0084 | 0.8559 | No unique incremental information detected |

In participant-held-out prediction, adding the plasma block reduced RMSE from
0.7942 to 0.6659, equivalent to a 29.69% reduction in mean squared error. The
participant bootstrap one-sided p-value was 0.0001.

A subsequent 20-assignment subject-fold sensitivity analysis preserved the
unique-marker distinction in every assignment: p-tau217 and A-beta42/40 were
Holm-detected in 20/20, while NfL and GFAP were detected in 0/20. The plasma
MSE reduction was positive in 20/20 assignments (25.4% to 29.3%). This checks
split stability within the same cohort and is not an external replication.

## Common-cohort MRI comparison

A separate comparison retained 641 observations from 570 participants with a
quality-controlled structural MRI within 90 days of PET. On this common cohort:

- plasma remained detectable after conditioning on the MRI block (`p = 0.0001`);
- MRI was not detectable after conditioning on plasma (`p = 0.2378`);
- plasma reduced held-out MSE by 28.78%;
- MRI reduced held-out MSE by 0.30%; and
- both blocks together reduced held-out MSE by 29.25%.

The comparison uses one cohort and one fold assignment. It does not rank
p-values obtained from different sample sizes.

## Interpretation boundary

The completed analysis distinguishes a candidate information signature:
p-tau217 and A-beta42/40 retained amyloid information under the specified
conditioning set, whereas NfL and GFAP did not. This conclusion is tied to the
selected ADNI cohort, date windows, measurements, and reference variables.
It does not establish temporal precedence, causal mechanisms, transportability,
or clinical utility.

The frozen configuration is `configs/adni_plasma_amyloid.yaml`. Public aggregate
results are stored under `results_public/adni/`. Reproduction requires separate
ADNI approval and locally downloaded source tables.
