# Validation report: Case-B proxy quality and B/C1 local power

Date: 24 August 2026

Configuration: `configs/non_iid_bc_gradients.yaml`

## Frozen questions

1. How accurately must a stable subject context be estimated from independent
   auxiliary measurements before treating it as observed becomes a reasonable
   approximation in this generator?
2. Do the observed-context Case-B and complete-history Case-C1 procedures show
   graded power under local residual-covariance alternatives?

The proxy experiment is diagnostic. The auxiliary measurements are independent
of the test-period innovations, so its results do not validate an arbitrary
BLUP, embedding, or generated covariate.

## Design

- 250 independent subjects;
- four visits in Case B and eight retained visits in Case C1;
- five subject-level folds;
- two-sided cluster-score test at 0.05;
- 1,000 repetitions per proxy-quality row;
- 500 repetitions per power row; and
- exact binomial 95% intervals.

The target null is true given the stable subject context in Case B and the
complete one-lag history in Case C1. Frozen null anchors for perfectly observed
context and complete history remain the earlier 1,000-repetition experiment;
they were not rerun or selected after seeing these gradients.

## Proxy-quality results

| Auxiliary measurements per subject | Rejections / 1,000 | Rate | Exact 95% interval |
|---:|---:|---:|---:|
| 4 | 999 | 0.999 | 0.9944--1.0000 |
| 20 | 292 | 0.292 | 0.2640--0.3213 |
| 100 | 72 | 0.072 | 0.0568--0.0898 |
| 500 | 59 | 0.059 | 0.0452--0.0754 |
| 2,000 | 58 | 0.058 | 0.0443--0.0743 |

The 500-measurement row missed the pre-existing strict upper-confidence-limit
criterion of 0.075 by 0.00045. The 2,000-measurement row passed. These numbers
describe this data-generating process rather than a universal measurement-count
threshold. They show why the manuscript states a rate condition on the generated
context instead of treating every estimated subject effect as observed.

## Local power

| Residual-covariance coefficient | Case B observed context | Case C1 complete history |
|---:|---:|---:|
| 0.05 | 172/500 (0.344) | 301/500 (0.602) |
| 0.10 | 432/500 (0.864) | 498/500 (0.996) |
| 0.20 | 500/500 (1.000) | 500/500 (1.000) |
| 0.30 | 500/500 (1.000) | 500/500 (1.000) |

Together with the frozen null anchors (0.052 for observed Case B and 0.030 for
complete-history C1), the new rows show graded sensitivity at small effects.
The different B and C1 power curves should not be used to rank the cases because
their generators and effective within-subject information differ.

## Reproducibility hashes (SHA-256)

- configuration: `0e8d54662bcfb1cb7681da72111b7211fd914790e56d86d1245fd06ef307fc38`
- simulation code: `cb83c321b8e71252d3c8dc2ede06f4514b82f392c08d5cdca96affe92de21fb7`
- aggregate summary: `7219a60e54cff6c67207e1be5865143c9b4c1b72fe5cb89e781007558224f22b`
- replicate results: `fc134d143ad06986e6bd108988c7f0bfd2f9bde5beca4c29353eda8d4175db64`

Replicate-level results remain outside the reviewed public artifact. The public
repository contains the configuration, runner, validation report, and aggregate
summary only.
