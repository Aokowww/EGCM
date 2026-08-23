# Validation report: ADNI repeated subject-fold stability

Date: 24 August 2026

## Purpose

The primary ADNI analyses used pre-specified subject-level folds. This
supplementary sensitivity analysis repeats the fold assignment 20 times while
holding the cohort, variables, formulas, and scientific conditioning sets
fixed. It asks whether the candidate-information signature or held-out
prediction gain is an artefact of one split.

Repeated folds are not independent cohorts and are not pooled as replication
evidence. The originally frozen analyses remain primary.

## Plasma information about amyloid PET (Case A1)

The 758-row, 664-participant plasma--Centiloid cohort was unchanged. For each
fold assignment, all four unique-marker tests were recomputed by conditioning
on the other plasma markers, and the reference-versus-plasma prediction
comparison was refitted.

| Marker | Normalised covariance range | Same sign | Holm detection rate | Largest Holm p-value |
|---|---:|---:|---:|---:|
| p-tau217 | 0.456 to 0.497 | 20/20 | 20/20 | $6.36\times10^{-10}$ |
| A-beta42/40 | -0.276 to -0.244 | 20/20 | 20/20 | $1.78\times10^{-8}$ |
| NfL | -0.119 to -0.067 | 20/20 | 0/20 | 0.397 |
| GFAP | -0.016 to 0.041 | 19/20 relative to the median sign | 0/20 | 0.846 |

The plasma-block MSE reduction was positive in every split, with a median of
27.9% and a range of 25.4% to 29.3%. Thus the distinction between the two
amyloid-informative measurements and NfL/GFAP did not depend on the original
fold assignment.

## Plasma p-tau217 information about future ADAS-Cog13 (Case C1)

The primary 30-day cohort of 199 participants and its complete observed-history
conditioning set were unchanged. Each repeated split recomputed the p-tau217
unique-marker score given the other three plasma measurements and compared
history-only with history-plus-p-tau217 prediction.

The normalised residual covariance was positive in 20/20 splits, with a median
of 0.297 and a range of 0.230 to 0.351. A conservative four-marker Bonferroni
adjustment gave $p<0.05$ in 19/20 splits; the largest adjusted p-value was 0.098.
The MSE reduction was positive in 20/20 splits, with a median of 21.0% and a
range of 14.6% to 26.6%.

The prospective direction and predictive gain are therefore stable, while the
binary significance threshold shows mild fold sensitivity in the smaller C1
cohort. The paper reports both facts.

## Privacy and interpretation

Only aggregate split summaries are public. Participant identifiers, dates,
fold membership, residuals, and predictions remain in the restricted ADNI
workspace. These analyses support stability within the selected ADNI cohorts;
they do not supply external validation, causal effects, or clinical diagnostic
performance.

## Reproducibility hashes (SHA-256)

### Plasma--amyloid

- configuration: `1327be144ff8d62b2967465a4c874199e9a17c3130dbd6793d1d8fa140663303`
- runner: `68b842d657e40d752675e5fffd032b55a323702fb181a9b274a80ed9e990ba28`
- aggregate summary: `0c00c0e32d6d61cd22f2cc4e14ef6e0230e1e74589887ba10accc7a31d56891a`

### Prospective Case C1

- configuration: `70bca1df4fc242eb7f6fef3abda4189c38035cac1c495f503d2432ebceeae80b`
- runner: `b4e2be9b7035c2655437e6d78d83ad9c4493fe23b2a0616ff9049538f6310f75`
- aggregate summary: `28c5917892c518ab90d8e6d1ad7ffc9736a7468f3a93dd9d73ea5dd6bb5c731e`
