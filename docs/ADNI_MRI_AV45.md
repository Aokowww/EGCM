# Exploratory ADNI MRI-AV45 analysis

This earlier analysis used structural MRI markers and continuous AV45 PET. The
primary adjustment set included concurrent cognition; a separate sensitivity
analysis omitted concurrent cognition.

The aligned data contained 672 observations from 467 participants. The median
participant contributed one complete observation and the maximum was three.
This sparse structure did not support a useful Design B calculation.

## Design A results

| Adjustment set | Global statistic | Multiplier p |
|---|---:|---:|
| Concurrent cognition and baseline covariates | 1.7663 | 0.2124 |
| Baseline covariates without concurrent cognition | 6.3355 | 0.0001 |

The cognition-adjusted analysis did not reject its population-level null. The
baseline-only sensitivity analysis did. These analyses answer different
questions and should not be substituted for one another.

| Adjustment set | Marker | Statistic | Holm p | Normalized covariance |
|---|---|---:|---:|---:|
| Concurrent cognition | Ventricles | 0.825 | 0.8185 | 0.041 |
| Concurrent cognition | Hippocampus | -1.766 | 0.2320 | -0.070 |
| Concurrent cognition | Entorhinal | -0.604 | 0.8185 | -0.027 |
| Baseline covariates | Ventricles | 2.391 | 0.0168 | 0.134 |
| Baseline covariates | Hippocampus | -6.335 | 7.10e-10 | -0.277 |
| Baseline covariates | Entorhinal | -4.373 | 2.45e-5 | -0.203 |

The multiplier value 0.0001 is the resolution limit from 9,999 draws. Failure
to reject in the primary analysis is not proof of conditional independence.
The sensitivity associations are not causal estimates.

The complete-case reduction, sparse visits and influential subject scores
limit the interpretation. The MRI-ADAS13 analysis uses a different outcome,
conditioning set and visit restriction; its results do not replace this
analysis.
