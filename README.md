# Conditional independence testing beyond i.i.d. data

This repository contains a taxonomy and residual-score framework for
conditional-independence testing (CIT) with non-i.i.d. data, together with the
code used for its simulation checks and exploratory ADNI analyses. The
theoretically developed branch is the many-independent-subject regime; the
method changes the unit of inference from visits to subjects.

The earlier thesis implementation is preserved in the
[`v1.0-thesis-legacy`](https://github.com/Aokowww/EGCM/tree/v1.0-thesis-legacy)
tag. It fitted mixed models and then applied a visit-level GCM calculation.
That implementation remains available for historical reproduction, but it is
not the method used in the current analysis.

The current full-length working article is available as a
[compiled PDF](article/cluster_gcm_working.pdf), with its
[LaTeX source and figures](article/README.md).

## Current status

The method is still under development. The evidence is mixed, and the
repository reports the unsuccessful checks as well as the successful ones.

| Check | Result |
|---|---|
| Design A, core population null | 47/1,000 rejections; passed the frozen criterion |
| Design A, six stress settings | Five passed; the random-slope setting failed |
| Earlier Design B, estimated BLUP context | 29/100 and 39/100 rejections; sensitivity only |
| Case B, observed subject context | 52/1,000 rejections; passed the extension criterion |
| Case B, dense independent context proxy | 54/1,000 rejections; favourable heuristic only |
| Case B, proxy-quality gradient | 0.999, 0.292, 0.072, 0.059 and 0.058 rejection at 4, 20, 100, 500 and 2,000 auxiliary measurements |
| Case C, complete finite history | 30/1,000 rejections; passed, conservatively |
| Case B/C1 local power | At effect 0.05: 0.344/0.602; at effect 0.10: 0.864/0.996 |
| A1 inferential-unit comparison | Subject-score: 60/1,000 and 55/1,000; visit-iid: 216/1,000 and 240/1,000 |
| ADNI MRI-ADAS13 analysis | Exploratory association analysis; 647 participants and 4,406 observations |
| ADNI plasma-Centiloid analysis | Exploratory incremental-information analysis; 664 participants and 758 observations |
| ADNI plasma-future ADAS13 Case C | Exploratory prospective analysis; 199 participants in the primary cohort; p-tau217 retained unique signal |
| ADNI repeated subject folds | Plasma marker distinction held in 20/20; prospective p-tau217 adjusted detection held in 19/20 and prediction gain in 20/20 |

The complete Stage 2b gate failed because the random-slope setting had 65
rejections in 1,000 repetitions and an exact 95% interval of 0.0505 to
0.0821. The frozen upper-limit criterion was 0.075. Design A must therefore
not be described as generally calibrated under random slopes.

Design B is included as a model-based sensitivity calculation. Its null
rejection rates were too high for confirmatory use.

## Taxonomy and statistical target

The article separates five targets before selecting an estimator:

- **A:** same-occasion CIT given observed visit-level context;
- **B:** same-occasion CIT additionally given subject context;
- **C:** same-occasion CIT given observed history, with bounded-history C1 and growing-history C2;
- **D:** trajectory-level CIT given an observed context trajectory; and
- **E:** trajectory-level CIT additionally given subject context.

Cases D and E are future work. Observed Case-B context and complete bounded
Case-C1 history can be appended to the conditioning set. Estimated or latent
context requires additional identification, convergence and leakage
conditions; a sparse subject-effect proxy is not automatically valid. C2
requires temporal weak-dependence or martingale theory rather than the current
subject-score theorem.

For independent subjects \(i=1,\ldots,N\), with repeated visits
\(t=1,\ldots,m_i\), Design A examines the population-level same-visit null

\[
H_0: X_{it} \perp\!\!\!\perp Y_{it}\mid Z_{it}.
\]

Latent subject effects are integrated out rather than conditioned upon. The
calculation:

1. assigns all visits from a subject to the same cross-fitting fold;
2. estimates population conditional means without a held-out-subject BLUP;
3. averages residual products within each subject;
4. studentizes across the \(N\) independent subject scores; and
5. uses subject-level multipliers for the global multi-marker test.

The scalar test targets a residual covariance. Conditional independence
implies a zero target, but a zero target alone does not establish conditional
independence.

The article proves the scalar A1 result under stated high-level
observation-process, nuisance-rate, cross-term, and variance conditions and
gives a fixed-marker Gaussian multiplier proposition. These results do not
automatically verify the assumptions for every nuisance learner; a concrete
GAMM, kernel, or neural implementation must satisfy the stated conditions.

## Simulation evidence

The public validation notes record the frozen criteria, exact confidence
intervals, integrity checks and file hashes:

- [cluster-GCM method](docs/METHOD_CLUSTER_GCM.md)
- [Stage 2a core-null validation](docs/VALIDATION_STAGE2A_CORE.md)
- [Stage 2b stress validation](docs/VALIDATION_STAGE2B_STRESS.md)
- [Case B/C frozen plan](docs/EXPERIMENT_PLAN_2026-08-23_NONIID_BC.md)
- [Case B/C validation](docs/VALIDATION_2026-08-23_NONIID_BC.md)
- [A1 inferential-unit comparison](docs/VALIDATION_2026-08-24_METHOD_COMPARISON.md)
- [Case-B proxy-quality and B/C1 local-power gradients](docs/VALIDATION_2026-08-24_NONIID_BC_GRADIENTS.md)

Short, machine-readable summaries are stored in
[`results_public/simulation`](results_public/simulation). The full checkpoint
directories are omitted from Git because they contain thousands of
intermediate files.

## Exploratory ADNI analyses

The current MRI-ADAS13 experiment anchors each observation on an ADAS13
assessment and selects the nearest quality-controlled FreeSurfer 7 MRI scan
within 30 days. The common A/B cohort requires at least five complete visits.
It contains 4,406 observations from 647 participants; the median number of
visits is six.

The conditioning set contains time since baseline, baseline age, sex,
education, APOE4, intracranial volume, protocol, site, MRI field strength and
the MRI-ADAS date gap. Concurrent diagnostic and cognitive variables are not
included.

For Design A, the global three-marker multiplier test reached the resolution
limit of 0.0001 with 9,999 draws. Residual ventricular volume was positively
associated with ADAS13, while residual hippocampal and entorhinal volumes were
negatively associated with ADAS13. These are contemporaneous conditional
associations. They do not show that structural change preceded cognitive
change, and they are not causal estimates.

The prospective Case-C analysis asks whether plasma p-tau217, A-beta42/40,
NfL and GFAP contain information about future ADAS-Cog13 beyond current and
prior cognition, demographics, observed APOE4 context and timing. The primary
30-day cohort contains 199 participants. The global four-marker test reached
0.0001, and only p-tau217 retained a unique marker signal after conditioning
on the other plasma measurements. A 90-day alignment sensitivity reproduced
the marker distinction. These results describe conditional residual
association and held-out prediction, not causality or clinical validation.

Repeating the subject-fold assignment 20 times preserved the plasma--amyloid
marker distinction in every assignment. In the prospective C1 cohort,
p-tau217 retained a positive unique score and a positive prediction gain in
all 20 assignments; a conservative four-marker adjusted test was below 0.05
in 19. Repeated folds quantify split sensitivity within the same cohorts and
are not independent external validation.

The cohort construction, diagnostics and bounded interpretation are reported
in:

- [MRI-ADAS13 experiment report](docs/ADNI_MRI_ADAS13.md)
- [earlier MRI-AV45 Design A analysis](docs/ADNI_MRI_AV45.md)
- [plasma-Centiloid and common-cohort MRI comparison](docs/ADNI_PLASMA_AMYLOID.md)
- [prospective plasma-future ADAS13 Case-C analysis](docs/ADNI_PROSPECTIVE_PTAU_CASE_C.md)
- [ADNI repeated subject-fold stability](docs/VALIDATION_2026-08-24_ADNI_SPLIT_STABILITY.md)
- [public aggregate tables](results_public/adni)

## Repository layout

```text
code/adni/            cluster-GCM implementation and ADNI preparation scripts
code/simulation/      frozen simulation drivers and data-generating processes
configs/              analysis and simulation settings
docs/                 design decisions and final validation reports
results_public/       reviewed aggregate results without participant-level data
tests/                checks for the score calculation and stress generators
```

The older thesis scripts remain under their original paths so that the legacy
tag and previous figures can still be traced.

## Running the checks

Run these commands from the repository root:

```bash
Rscript tests/test_cluster_gcm_core.R
Rscript tests/test_adni_stress_dgp.R
```

The core Stage 2a configuration is:

```bash
Rscript code/simulation/simulate_adni_designs_ab.R \
  configs/adni_simulation_stage2a_core_a_null.yaml
```

The six-setting Stage 2b stress run is:

```bash
Rscript code/simulation/simulate_adni_design_a_stress.R \
  configs/adni_simulation_stage2b_stress.yaml
```

The Case B/C extension is:

```bash
Rscript code/simulation/simulate_non_iid_bc.R \
  configs/non_iid_bc_simulation.yaml
```

The direct inferential-unit comparison is:

```bash
Rscript code/simulation/simulate_method_comparison.R \
  configs/method_comparison.yaml
```

The Case-B proxy-quality and B/C1 local-power gradients are:

```bash
Rscript code/simulation/simulate_non_iid_bc_gradients.R \
  configs/non_iid_bc_gradients.yaml
```

ADNI analyses require an approved ADNI account and locally prepared source
tables:

```bash
Rscript code/adni/prepare_adni_mri_adas13.R \
  data/adni_raw/YYYY-MM-DD 30 5

Rscript code/adni/run_adni_redesign.R \
  configs/adni_mri_adas13_exploratory.yaml

ADNI_RESTRICTED_ROOT=/path/to/private/adni \
  Rscript code/adni/run_adni_prospective_ptau_case_c.R \
  configs/adni_prospective_ptau_case_c.yaml

ADNI_RESTRICTED_ROOT=/path/to/private/adni \
  Rscript code/adni/run_adni_plasma_amyloid_split_stability.R \
  configs/adni_plasma_split_stability.yaml

ADNI_RESTRICTED_ROOT=/path/to/private/adni \
  ADNI_AGGREGATE_OUTPUT=/path/to/private/aggregate/prospective_ptau_split \
  Rscript code/adni/run_adni_prospective_ptau_case_c.R \
  configs/adni_prospective_ptau_split_stability.yaml
```

Package requirements and environment notes are in
[docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md).

## Data access

ADNI data are not distributed in this repository. Local analytic CSV files,
source tables, participant identifiers and serialized model objects are
excluded by `.gitignore`. Researchers must obtain access from ADNI and accept
the applicable data-use terms before running the ADNI scripts.

See [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md) for the expected local file
layout. Only reviewed aggregate tables are versioned here.

## Citation and licence

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). The code is
released under the MIT licence. ADNI data remain subject to ADNI's own access
and use conditions.
