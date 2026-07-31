# Stage 2a: core-null validation

The Stage 2a run evaluated Design A under the core population null. The
configuration was fixed before the 1,000 repetitions were aggregated.

| Repetitions | Rejections | Rate | Exact 95% interval |
|---:|---:|---:|---:|
| 1,000 | 47 | 0.047 | 0.0347 to 0.0620 |

The decision rule required 1,000 valid repetitions and an exact 95% upper
confidence limit no greater than 0.075. Stage 2a passed this rule.

Rejection was defined as a global multiplier p-value strictly below 0.05.
Three p-values were exactly 0.05 and were not counted.

## Integrity checks

- 1,000 unique repetition identifiers were present.
- No aggregate cells were missing.
- Every simulated data set contained 120 subjects.
- The subject-score unit test passed.
- Reconstructing the aggregate files from the checkpoints changed no results.

| Output | SHA-256 |
|---|---|
| `replicate_results.csv` | `822ce936671e974df054666a3f93af3b6a67d8412e986cc0a57fcd7adaa04225` |
| `summary.csv` | `18b3e3b8b1cdf0c99c2f0129fbc031c2281a1c63c455149cdbeb22819ede0459` |
| `STATUS.txt` | `0b2de562c37f758450abae10a77aecf837ff3fd5a57f6f7e2f5e740a50821242` |

Passing this core setting did not establish calibration under random slopes,
serial dependence, unequal visit counts or informative visit counts. Those
settings were evaluated separately in Stage 2b.
