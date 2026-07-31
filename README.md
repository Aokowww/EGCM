# Cluster-GCM for repeated-measures data

This repository contains the current cluster-GCM redesign of EGCM, together
with the code used for its simulation checks and exploratory ADNI analyses.
The redesign changes the unit of inference from visits to subjects.

The earlier thesis implementation is preserved in the
[`v1.0-thesis-legacy`](https://github.com/Aokowww/EGCM/tree/v1.0-thesis-legacy)
tag. It fitted mixed models and then applied a visit-level GCM calculation.
That implementation remains available for historical reproduction, but it is
not the method used in the current analysis.

## Current status

The method is still under development. The evidence is mixed, and the
repository reports the unsuccessful checks as well as the successful ones.

| Check | Result |
|---|---|
| Design A, core population null | 47/1,000 rejections; passed the frozen criterion |
| Design A, six stress settings | Five passed; the random-slope setting failed |
| Design B, two relevant null settings | 29/100 and 39/100 rejections |
| ADNI MRI-ADAS13 analysis | Exploratory association analysis; 647 participants and 4,406 observations |

The complete Stage 2b gate failed because the random-slope setting had 65
rejections in 1,000 repetitions and an exact 95% interval of 0.0505 to
0.0821. The frozen upper-limit criterion was 0.075. Design A must therefore
not be described as generally calibrated under random slopes.

Design B is included as a model-based sensitivity calculation. Its null
rejection rates were too high for confirmatory use.

## Statistical target

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

The theorem statements for the new cluster-GCM method remain proposed
results with a proof programme. They are not presented as completed theorems.

## Simulation evidence

The public validation notes record the frozen criteria, exact confidence
intervals, integrity checks and file hashes:

- [cluster-GCM method](docs/METHOD_CLUSTER_GCM.md)
- [Stage 2a core-null validation](docs/VALIDATION_STAGE2A_CORE.md)
- [Stage 2b stress validation](docs/VALIDATION_STAGE2B_STRESS.md)

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

The cohort construction, diagnostics and bounded interpretation are reported
in:

- [MRI-ADAS13 experiment report](docs/ADNI_MRI_ADAS13.md)
- [earlier MRI-AV45 Design A analysis](docs/ADNI_MRI_AV45.md)
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

ADNI analyses require an approved ADNI account and locally prepared source
tables:

```bash
Rscript code/adni/prepare_adni_mri_adas13.R \
  data/adni_raw/YYYY-MM-DD 30 5

Rscript code/adni/run_adni_redesign.R \
  configs/adni_mri_adas13_exploratory.yaml
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
