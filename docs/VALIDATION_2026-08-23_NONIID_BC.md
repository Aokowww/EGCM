# Validation report: non-i.i.d. CIT Cases B and C

Date: 23 August 2026

Configuration: `configs/non_iid_bc_simulation.yaml`

## Frozen design

- 1,000 Monte Carlo repetitions per scenario
- 250 independent subjects
- five subject-level folds
- two-sided cluster-score test at 0.05
- exact binomial 95% intervals
- primary gate: the exact upper interval endpoint is at most 0.075 for
  `B_observed_null` and `C_complete_history_null`

The primary gate passed.

## Results

| Scenario | Rejections / 1,000 | Rate | Exact 95% interval | Interpretation |
|---|---:|---:|---:|---|
| B: observed subject context, null | 52 | 0.052 | 0.0391–0.0676 | Primary calibration passed |
| B: sparse auxiliary proxy, target null | 1,000 | 1.000 | 0.9963–1.0000 | Proxy error leaves shared context in both residuals |
| B: dense auxiliary proxy, target null | 54 | 0.054 | 0.0408–0.0699 | Favourable generated-context heuristic |
| B: context omitted, target null | 1,000 | 1.000 | 0.9963–1.0000 | Tests a different reduced conditioning set |
| B: observed context, alternative | 1,000 | 1.000 | 0.9963–1.0000 | Power for the simulated residual-covariance alternative |
| C: complete finite history, null | 30 | 0.030 | 0.0203–0.0426 | Primary calibration passed; conservative in this DGP |
| C: partial history, target null | 1,000 | 1.000 | 0.9963–1.0000 | Omitted `X_(t-1)` changes the tested null |
| C: history omitted, target null | 1,000 | 1.000 | 0.9963–1.0000 | Tests a different reduced conditioning set |
| C: complete history, alternative | 1,000 | 1.000 | 0.9963–1.0000 | Power for the simulated contemporaneous innovation alternative |

## Claim boundary

The B observed-context and C complete-history rows are direct many-subject reductions to
Case A1. The dense-proxy row is deliberately favourable: the proxy is estimated from many
auxiliary measurements independent of test-period innovations. It is not evidence that an
arbitrary BLUP or latent embedding is valid. The sparse-proxy failure agrees with the earlier
Design B failure under bounded visits.

The partial- and omitted-history rows are not Type-I-error estimates for the reduced nulls.
The target null is true only after conditioning on the complete generated history. Their
rejection rates demonstrate that incomplete domain knowledge changes the hypothesis rather
than merely reducing efficiency.

## Reproducibility hashes (SHA-256)

- configuration: `aee6ba192f65adb55d1e9e1650aae92c5cea181b3a27a216d392066eed7c201a`
- simulation code: `6bab65160206c0c8c62a490b73947b1bacabe7d8353a650d24b9c6328bbf0b45`
- frozen plan: `544a7bae30e6158652da3fd72733ef717ee62801af5cb2994b0027664bfad4bc`
- aggregate summary: `beb66ef1df74b33bdc329cd28de6be799e0a1899d611e0e7f73b160f6cdc9896`
- replicate results: `6ab2f53b6165cd34aa607ec830f76cc81c031f9804b24ef8ebf7ba0b35d6f98f`
