# Frozen Case B/C simulation plan

Date frozen: 23 August 2026

## Purpose

This experiment separates three questions that are easily conflated in non-i.i.d.
conditional-independence testing:

1. whether the target conditioning information is observed;
2. whether subject context can be estimated accurately enough to be treated as a
   generated conditioning variable; and
3. whether a finite observed history has been included completely.

The simulation is an extension study. It does not change the previously frozen Stage 2a
or Stage 2b results.

## Case B: subject context

The target null is

\[
X_{it}\perp Y_{it}\mid (Z_{it},\theta_i).
\]

The data-generating process has a shared subject context `theta_i`. Under the null, the
innovations of `X` and `Y` are independent after conditioning on `(Z, theta)`. Five rows are
reported:

- observed `theta` under the null (primary calibration row);
- a sparse auxiliary proxy based on four calibration measurements;
- a dense auxiliary proxy based on 2,000 calibration measurements;
- omitted `theta`, which intentionally tests a different, reduced conditioning set; and
- observed `theta` under a residual-covariance alternative.

The sparse and dense proxies are generated from measurements independent of the test-period
innovations. They isolate generated-regressor error from leakage. The dense proxy is a
favourable heuristic case, not a theorem for arbitrary latent-effect estimators.

## Case C: finite observed history

The target null is

\[
X_{it}\perp Y_{it}\mid (Z_{it},H_{it}),\qquad
H_{it}=(X_{i,t-1},Y_{i,t-1}).
\]

The process allows both series to depend on their past and allows `Y_t` to depend on
`X_{t-1}`. Under the null, the contemporaneous innovations are independent. Four rows are
reported:

- complete observed history under the null (primary calibration row);
- only `Y_{t-1}` included;
- all history omitted; and
- complete history under a contemporaneous innovation alternative.

The incomplete-history rows are scientific-conditioning diagnostics, not Type-I-error
assessments of the reduced nulls.

## Inference and frozen criterion

- Independent unit: subject.
- Nuisance fitting: five-fold subject-level cross-fitting using linear regressions.
- Subject score: mean residual product across analysed visits.
- Test: two-sided normal-reference cluster-score test at `alpha = 0.05`.
- Monte Carlo repetitions: 1,000 per row.
- Uncertainty: exact binomial 95% interval.
- Primary gate: for each primary calibration row, the exact 95% upper confidence limit must
  not exceed 0.075.

The simulation uses linear nuisance models because the conditional means are linear by
construction. Its purpose is to isolate information structure and inferential units, not to
compare regression algorithms.
