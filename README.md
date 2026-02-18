# EGCM: An Extension of the Generalized Covariance Measure (GCM)

This repository contains code and materials accompanying the master thesis:

> **An Extension of the Generalized Covariance Measure for Conditional Independence Testing**

The goal is to extend the Generalised Covariance Measure (GCM) framework for conditional independence (CI)
testing to **non-i.i.d. repeated-measurements data** by replacing the regression step with a
**Generalized Additive Mixed Model (GAMM)** (implemented via `gamm4`).

## Core idea

Given repeated measurements $(X_{it}, Y_{it}, Z_{it})$ for subject $i$ at time $t$, the EGCM workflow is:

1. Fit GAMMs for $\mathbb{E}[X_{it}\mid Z_{it}]$ and $\mathbb{E}[Y_{it}\mid Z_{it}]$ to account for
   within-subject dependence.
2. Compute residuals
   $$\varepsilon_{it}=X_{it}-\widehat{\mathbb{E}}[X_{it}\mid Z_{it}],\quad \xi_{it}=Y_{it}-\widehat{\mathbb{E}}[Y_{it}\mid Z_{it}].$$
3. Apply the GCM test to the residuals (or their products) to test
   $$H_0: X \perp\!\!\!\perp Y\mid Z.$$

## Simulation design (Chapter 6.1)

We consider $n\in\{20,50,100,200\}$ subjects and $h_i\in\{3,5,10,20,50,100\}$ measurements per subject.
Within-subject covariates are generated from a multivariate normal distribution with covariance
$\Sigma$ such that $\Sigma_{tt}=1$ and $\Sigma_{ts}=0.5$ for $t\neq s$.

### Scenario (a): simple correlation

$$
X_{it}=f(Z_{it})+0.3\,\varepsilon_{it},\qquad
Y_{it}=g(Z_{it})+bX_{it}+0.3\,\xi_{it}.
$$

### Scenario (b): random intercepts

$$
X_{it}=f(Z_{it})+\alpha_i+0.3\,\varepsilon_{it},\qquad
Y_{it}=g(Z_{it})+\beta_i+bX_{it}+0.3\,\xi_{it},
$$

with $(\alpha_i,\beta_i)$ subject-specific intercepts.

## ADNI case study (Chapter 6.2)

The ADNI analysis focuses on testing whether MRI markers
$X=(\texttt{Ventricles},\texttt{Hippocampus},\texttt{Entorhinal})$
are CI of amyloid PET signal $Y=\texttt{AV45}$ given cognitive tests
$Z=(\texttt{MOCA},\texttt{MMSE},\texttt{CDRSB})$.

**Note:** this repository does not distribute ADNI data. See `docs/DATA_ACCESS.md`.

## Repository layout

- `code/simulation/`: simulation scripts.
- `code/adni/`: ADNI analysis scripts.
- `code/notes/`: exploratory R Markdown notes.
- `configs/`: YAML configs documenting key settings.
- `thesis/`: the thesis PDF.

## Requirements

- R (>= 4.2 recommended)
- Key packages: `mgcv`, `gamm4`, `lme4`, `GeneralisedCovarianceMeasure`, `MASS`, `tidyverse`

For reproducibility, consider using `renv` (see `docs/REPRODUCIBILITY.md`).

## How to run

From the repository root:

```r
# Example: run one simulation script
source("code/simulation/sim_correlatedz.R")

# Example: run ADNI analysis
source("code/adni/ADNI_exp.R")
```

## License

MIT (see `LICENSE`).
