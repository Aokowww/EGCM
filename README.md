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

## Simulation studies (repeated measurements)

We validate the original GCM and the proposed EGCM on non-i.i.d. longitudinal/repeated-measurements data.
Across all simulation experiments we report **Type I error** (under the null) and **power** (under the alternative),
based on **100 Monte Carlo repetitions** per setting.

### Common notation

We generate repeated measurements
\[
\{(X_{it}, Y_{it}, Z_{it})\}_{i=1,\dots,n;\; t=1,\dots,h_i},
\]
where $i$ indexes subjects and $t$ indexes within-subject measurement occasions.
The conditional independence (CI) hypothesis of interest is
\[
H_0:\; X \perp\!\!\!\perp Y \mid Z.
\]

EGCM differs from GCM only in the regression step: it estimates conditional means using **GAMMs** to account
for within-subject dependence before applying a GCM-style residual covariance test.

---

### Simulation 1: Random intercept models (nonlinear mean functions)

This experiment corresponds to the main repeated-measurement toy setting.

**Design.** We vary the number of subjects and measurements as
\[
n\in\{20,50,100,200\},\qquad
h_i\in\{3,5,10,20,50,100\}.
\]
Within-subject covariates $Z_{it}$ are generated from a multivariate normal distribution with covariance matrix
$\Sigma$ such that $\Sigma_{tt}=1$ and $\Sigma_{ts}=0.5$ for $t\neq s$ (within-subject correlation).

We consider a collection of nonlinear mean functions (e.g. square/cubic/tanh/negative exponential) via
$f(\cdot)$ and $g(\cdot)$ to emulate diverse empirical relationships.

**Scenario (a): simple correlation structure**
\[
X_{it}=f(Z_{it})+0.3\,\varepsilon_{it},\qquad
Y_{it}=g(Z_{it})+bX_{it}+0.3\,\xi_{it}.
\]

**Scenario (b): random intercepts**
\[
X_{it}=f(Z_{it})+\alpha_i+0.3\,\varepsilon_{it},\qquad
Y_{it}=g(Z_{it})+\beta_i+bX_{it}+0.3\,\xi_{it},
\]
where $(\alpha_i,\beta_i)$ are subject-specific random intercepts (Gaussian).

**Null vs alternative.** We set $b=0$ (CI holds) to evaluate Type I error and $b=0.2$ to evaluate power.

**Implementation.**
- GCM uses GAM-based conditional mean estimation (i.i.d. assumption implicit).
- EGCM uses GAMMs with a random intercept for subject, e.g.
  `X ~ s(Z) + (1|subject)` and `Y ~ s(Z) + (1|subject)`.

**Scripts.**
- `code/simulation/sim_correlatedz.R`
- `code/simulation/sim_nonlinear_post.R`
- `code/notes/Note_random_intercept.Rmd`

---

### Simulation 2: Random slope models + missingness (robustness & misspecification)

This experiment stresses the methods under stronger subject-specific heterogeneity and missing data.

**Design.** We vary the number of subjects
\[
n\in\{20,50,100,200,300,400\},
\]
and fix the number of measurements per subject at $h=10$.
We generate $Z_{it}\sim \mathcal{N}(0,1)$, random intercepts $(\alpha_i,\beta_i)$, and random slopes
$(RE_{sx,i}, RE_{sy,i})$:
\[
\alpha_i,\beta_i \sim \mathcal{N}(0,5^2),\qquad
RE_{sx,i}, RE_{sy,i} \sim \mathcal{N}(0,1.5^2).
\]
We introduce missingness independently in $(X,Y,Z)$ with probability $0.2$.

**DGP.**
\[
X_{it}= \alpha_i + f(Z_{it}) + RE_{sx,i} Z_{it} + 0.3\,\varepsilon_{it},
\]
\[
Y_{it}= \beta_i + g(Z_{it}) + RE_{sy,i} Z_{it} + b\,f(Z_{it}) + 0.3\,\xi_{it}.
\]

**Model specification vs misspecification.**
We compare:
- correctly specified random effects: `(1 + Z | subject)`,
- misspecified random effects: `(1 | subject)` only.

**Scripts.**
- `code/simulation/random_slope.R`
- `code/simulation/Mimic_ADNI_random_slope.R`
- `code/simulation/Mimic_ADNI_zuni.R`
- `code/simulation/Mimic_ADNI_zuni - mis.R`

---

### Simulation 3: Synthetic ADNI-like data (irregular visit schedules)

This experiment mimics common properties of medical cohorts: irregular numbers of visits, incomplete records,
and heterogeneous subject trajectories.

**Key feature: variable number of measurements.**
We generate $h_i$ to emulate ADNI-like visit counts, sampling from a Gaussian distribution and truncating to a feasible range:
\[
h_i \sim \text{round}\big(\mathcal{N}(6.73, 4.75^2)\big),\qquad h_i \in \{1,\dots,24\}.
\]

We consider large-scale settings
\[
n\in\{200,600,1000,1200,1500\},
\]
and evaluate both a univariate and a bivariate $(X,Y,Z)$ configuration.

**Univariate configuration (CI under $b=0$).**
\[
X_{it}= f_1(Z_{1it}) + f_2(Z_{2it}) + \alpha_i + 0.3\,\varepsilon_{it},\qquad
Y_{it}= g_1(Z_{1it}) + g_2(Z_{2it}) + \beta_i + 0.3\,\xi_{it}.
\]

**Bivariate configuration (more complex dependence).**
\[
Z_{1it}, Z_{2it} \sim \mathcal{N}(0,1),\quad
\alpha_{1i},\alpha_{2i},\beta_{1i},\beta_{2i}\sim \mathcal{N}(0,5^2),
\]
\[
\begin{aligned}
X_{1it} &= f_1(Z_{1it}) + f_2(Z_{2it}) + \alpha_{1i} + 0.3\,\varepsilon_{1it}, \\
X_{2it} &= f_1(Z_{1it}) + f_2(Z_{2it}) + X_{1it} + \alpha_{2i} + 0.3\,\varepsilon_{2it}, \\
Y_{1it} &= g_1(Z_{1it}) + g_2(Z_{2it}) + \beta_{1i} + 0.3\,\xi_{1it}, \\
Y_{2it} &= g_1(Z_{1it}) + g_2(Z_{2it}) + Y_{1it} + \beta_{2i} + 0.3\,\xi_{2it}.
\end{aligned}
\]

**Script.**
- `code/simulation/Mimic_ADNI.R`


## ADNI case study (structural MRI vs amyloid PET, conditional on cognition)

### Scientific question (CI testing)

We test whether **structural MRI markers** provide information about early AD-related pathology beyond
standard cognitive assessments.

Let
- $X$ be a vector of structural MRI markers,
- $Z$ be cognitive test scores,
- $Y$ be an amyloid PET-derived continuous biomarker.

The statistical question is a conditional independence test:
\[
H_0:\; X \perp\!\!\!\perp Y \mid Z.
\]
Because ADNI is a longitudinal cohort with repeated observations per subject, the data are **non-i.i.d.**,
and EGCM is designed to account for within-subject dependence through mixed-effects smoothing.

### Variable specification

We follow the thesis variable choice:
\[
X=(\texttt{Ventricles},\texttt{Hippocampus},\texttt{Entorhinal}),\qquad
Z=(\texttt{MOCA},\texttt{MMSE},\texttt{CDRSB}),
\]
and
\[
Y=\texttt{AV45 ratio}=\frac{\text{Cortical Grey Matter}}{\text{Whole Cerebellum}}.
\]

### Why a continuous outcome ($Y$ = AV45 ratio)

A key practical choice is to use the **continuous** AV45 ratio rather than a **categorical** diagnosis label.
The standard GCM asymptotic justification relies on CLT-style arguments for residual products after conditional
mean estimation. When switching to classification likelihoods (e.g., log-likelihood-based fits),
additional theoretical work is typically needed to ensure analogous asymptotic validity.
Using a continuous PET biomarker keeps the analysis aligned with the standard GCM/EGCM setup.

### Modeling strategy (EGCM)

EGCM estimates conditional means using GAMMs with subject-level random effects, e.g.
- for a univariate marker: `X ~ s(Z1) + s(Z2) + s(Z3) + (1|subject)`,
- similarly for $Y$.

Residuals are formed as
\[
\varepsilon_{it}=X_{it}-\widehat{\mathbb{E}}[X_{it}\mid Z_{it}],\qquad
\xi_{it}=Y_{it}-\widehat{\mathbb{E}}[Y_{it}\mid Z_{it}],
\]
and CI is tested via a GCM-type statistic based on $\varepsilon_{it}\xi_{it}$.

### Multiple testing (optional per-marker analysis)

If analyzing each MRI marker separately,
\[
X_1 \perp\!\!\!\perp Y\mid Z,\quad
X_2 \perp\!\!\!\perp Y\mid Z,\quad
X_3 \perp\!\!\!\perp Y\mid Z,
\]
we recommend reporting FDR-adjusted $p$-values (e.g., Benjamini–Hochberg) to control false discoveries.

### Scripts

- Main analysis using continuous AV45 outcome: `code/adni/ADNI_exp.R`
- Alternative script (includes categorical endpoint exploration): `code/adni/ADNI.R`

### Data access & compliance

This repository **does not distribute ADNI data**.
To reproduce the analysis, you must obtain ADNI access independently and run the scripts locally.
See `docs/DATA_ACCESS.md` for the expected workflow and repository conventions.


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
