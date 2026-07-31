# Stage 2b: stress validation

Stage 2b contains six null settings, each with 1,000 repetitions. The frozen
decision rule required the exact 95% upper confidence limit to be no greater
than 0.075 in every setting.

| Setting | Rejections | Rate | Exact 95% interval | Decision |
|---|---:|---:|---:|---|
| Combined stress | 53/1,000 | 0.053 | 0.0399 to 0.0688 | Pass |
| Informative visits | 48/1,000 | 0.048 | 0.0356 to 0.0631 | Pass |
| One-sided nuisance misspecification | 52/1,000 | 0.052 | 0.0391 to 0.0676 | Pass |
| Random slopes | 65/1,000 | 0.065 | 0.0505 to 0.0821 | Fail |
| Serial AR(1) | 46/1,000 | 0.046 | 0.0339 to 0.0609 | Pass |
| Unequal visits | 56/1,000 | 0.056 | 0.0426 to 0.0721 | Pass |

The random-slope upper limit was 0.0821. The overall Stage 2b decision is
therefore **fail**, although the other five settings passed. The settings are
not pooled.

## Integrity checks

- All 6,000 expected repetitions were present.
- Scenario and repetition keys were unique.
- No expected key, p-value or aggregate cell was missing.
- Both R test scripts passed.
- Reconstruction from the 6,000 checkpoints reproduced the output hashes.
- No simulation process remained active after aggregation.

## Frozen implementation

| File | SHA-256 |
|---|---|
| `configs/adni_simulation_stage2b_stress.yaml` | `3d84ab40025df336c9f51da2c259a54124de2af0afd14367f83b95c5e9e138fa` |
| `code/simulation/simulate_adni_design_a_stress.R` | `bc5ee92b608c86790bfa51d5ca545c14b8b0a63b1944ce253d5081376ef19337` |
| `code/simulation/adni_stress_dgp.R` | `5f30e55e78b9d928a10c5797529b3713b0b5d3d853d3921d4fb218bd19b93e61` |
| `code/adni/cluster_gcm_core.R` | `f1881bd59ef9a82fc1f4f6f5b4969625dcc1413d661183bb85becd4c806a2cda` |

## Final outputs

| File | SHA-256 |
|---|---|
| `replicate_results.csv` | `6754ddb02998595693543b4cdc11193ec782d55fb83b2258bf023261e6e04d96` |
| `summary.csv` | `8298aa3db167855e51a69711494af2cda073eac7e6e849275468d949155776a6` |
| `STATUS.txt` | `70ec2a66c29ab42a8db48a6491e95ecbabe87eec502125cf844ece2d84ea9c60` |

The current method must not be described as generally calibrated under random
slopes. A change aimed at that setting would be a new method version and would
need a newly frozen simulation study.
