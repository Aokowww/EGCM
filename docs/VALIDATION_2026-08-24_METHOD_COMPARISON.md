# A1 inferential-unit comparison

Date: 2026-08-24

## Question

Does the A1 subject-score construction control false positives better than applying the same cross-fitted residual products as if all visits were independent?

## Frozen design

- 250 independent subjects and six visits per subject.
- Independent subject-specific effects in the two nuisance regressions.
- AR(1) within-subject innovations with correlation 0.55.
- Five subject-level folds.
- A base null, a random-slope null, and four contemporaneous residual-covariance alternatives.
- 1,000 repetitions per null and 500 per alternative.
- Primary method: subject-score GCM.
- Pre-specified calibration rule: the exact 95% upper confidence limit must not exceed 0.075.

The subject-score and visit-iid analyses use the same cross-fitted residuals. They differ only in whether the score and reference variance are constructed at the subject or visit level. The oracle benchmark uses the known conditional means and subject scores.

## Results

| Scenario | Subject-score GCM | Visit-iid GCM | Oracle subject score |
|---|---:|---:|---:|
| Base null | 60/1,000 (0.060) | 216/1,000 (0.216) | 57/1,000 (0.057) |
| Random-slope null | 55/1,000 (0.055) | 240/1,000 (0.240) | 57/1,000 (0.057) |

The exact 95% interval for the base-null subject-score row was 0.0461--0.0766. Its upper endpoint exceeded the frozen 0.075 ceiling by 0.0016, so this row did not pass the strict gate. The random-slope row passed (0.0417--0.0710). Both intervals include the nominal 0.05 level and closely track the oracle. The visit-iid analysis was severely anti-conservative in both null scenarios.

Under alternatives with residual-covariance coefficients 0.05, 0.10, 0.20, and 0.30, subject-score rejection rates were 0.176, 0.456, 0.942, and 1.000. Oracle rates were 0.182, 0.458, 0.942, and 1.000. Visit-iid rejection rates were larger, but cannot be interpreted as valid power because its null rejection rate was already inflated.

## Interpretation

This experiment isolates the independent-unit decision. It supports subject aggregation as the mechanism preventing large false-positive inflation. It does not erase the separate frozen Stage-2b random-slope failure, which used a different and more demanding generator and remains reported unchanged.

## Reproducibility

- Configuration SHA-256: `69a947f1cc251cfbdbd7ef7a87fac7ee40aff056546d51df826d1856b48ee908`
- Simulation script SHA-256: `1ce8fa18bbd2ddd05e1cd8cc1b866be8712ded3041ac8fbc88cc607ddca0cb62`
- Full summary SHA-256: `218e3c4b60980244e949c6c9625af129945034028d77084fd836a242bc18675a`

Replicate-level results remain in the non-public experiment directory. The public repository contains only the aggregate summary.
