# Prospective plasma information about future ADAS-Cog13

## Question

This exploratory analysis asks whether plasma biomarkers contain information
about future ADAS-Cog13 beyond current cognition, observed cognitive history,
demographics, observed APOE4 context and timing. It instantiates the
many-independent-participant branch of Case C. Fixed observed history is added
to the conditioning set; the analysis is not a single-series temporal test.

## Variables and alignment

The candidate block contains log-transformed p-tau217, A-beta42/40, NfL and
GFAP. A plasma assay date is the index date. The outcome is an ADAS-Cog13
assessment 270 to 1,278 days later, selecting the assessment closest to 730
days when several are available. Current ADAS-Cog13 is matched within 30 days
in the primary analysis and within 90 days in a sensitivity analysis. The
observed history is the most recent ADAS-Cog13 assessment 180 to 1,080 days
before the index.

The reference model contains age, sex, education, APOE4 dosage, current
ADAS-Cog13, prior ADAS-Cog13 and the alignment gaps. One earliest eligible
index is retained per participant. All nuisance predictions use the same five
participant-level folds.

## Cohort

The primary cohort contains 199 participants. The median absolute current
alignment gap is zero days. The median prior-assessment gap is 384 days. The
median future-assessment gap is 526 days, with an interquartile range of 373
to 740 days. The 90-day sensitivity cohort contains 203 participants.

## Conditional residual tests

The primary global four-marker multiplier test is 0.0001 using 9,999 draws.
Each unique-marker test conditions on the other three plasma markers.

| Marker | Primary normalized covariance | Primary Holm p | 90-day Holm p |
|---|---:|---:|---:|
| p-tau217 | 0.321 | 0.00123 | 0.000380 |
| A-beta42/40 | 0.096 | 0.3988 | 0.2078 |
| NfL | 0.057 | 0.7362 | 0.7131 |
| GFAP | -0.054 | 0.7362 | 0.8215 |

P-tau217 is the only marker retaining a detectable unique signal. Deleting
each participant in turn leaves its unique-marker p-value between
8.01e-5 and 5.59e-4.

## Held-out prediction

Adding p-tau217 to the history-conditioned model lowers primary-cohort RMSE
from 5.851 to 5.387, a 15.2% MSE reduction. The participant-bootstrap interval
for the mean MSE reduction is -0.496 to 11.006, with one-sided p = 0.0375. The
interval crosses zero, so this primary prediction comparison is suggestive
under a two-sided 95% criterion.

In the 90-day sensitivity cohort, RMSE falls from 5.886 to 5.351. The 17.3%
MSE reduction has an interval of 0.334 to 12.231 and one-sided p = 0.0186.
Adding observed history alone to the current-cognition model reduces MSE by
1.9% in the primary cohort. History remains in the reference set because it
defines the scientific question, not because it must increase predictive
accuracy.

## Interpretation boundary

The residual score is a conditional covariance moment, not an omnibus proof
of conditional independence. The prediction comparison uses one frozen fold
assignment and participant bootstrapping. ADNI is observational and the cohort
is selected by assay and follow-up availability. The analysis therefore
supports a candidate-information distinction, not a causal effect or a claim
of clinical biomarker validity.

Participant-level ADNI rows, identifiers, dates and predictions remain under
the restricted data directory. The repository contains only the runner,
configuration and reviewed aggregate summaries.
