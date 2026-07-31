# Cluster-GCM method

## Target

Subjects are independent. Measurements within a subject may be dependent. For
subject \(i\) and visit \(t\), the population-level target is

\[
H_0: X_{it} \perp\!\!\!\perp Y_{it}\mid Z_{it}.
\]

Subject-specific latent effects are integrated out. This differs from a null
that conditions on an unobserved subject effect \(U_i\).

Let

\[
f_0(z)=E(X_{it}\mid Z_{it}=z), \qquad
g_0(z)=E(Y_{it}\mid Z_{it}=z).
\]

For pre-specified weights \(w_{it}\) satisfying
\(\sum_t w_{it}=1\), define the subject score

\[
S_i=\sum_{t=1}^{m_i}w_{it}
\{X_{it}-f_0(Z_{it})\}\{Y_{it}-g_0(Z_{it})\}.
\]

The current implementation uses \(w_{it}=1/m_i\). Each subject therefore has
the same total weight, regardless of visit count.

## Calculation

Subjects, not visits, are assigned to cross-fitting folds. Nuisance models are
fitted on the training subjects. Predictions for a held-out subject target the
population conditional mean and exclude that subject's fitted random effect.

The estimated subject scores are studentized across subjects:

\[
T_N =
\frac{\sqrt{N}\,\bar S_N}
{\{N^{-1}\sum_i(\widehat S_i-\bar S_N)^2\}^{1/2}}.
\]

For several markers, the global statistic is the largest absolute
marker-specific statistic. Gaussian multipliers are drawn at subject level and
shared across markers so that the estimated marker dependence is retained.

## Design A and Design B

Design A is the population-level procedure described above. Its independent
unit is the subject.

Design B uses estimated subject effects and within-subject splitting. It was
included to study the latent-effect interpretation discussed during the
redesign. With a bounded number of visits, a fitted subject effect is a
shrunken estimate of \(U_i\), not an observed conditioning variable. Splitting
visits prevents direct reuse of a held-out response but does not remove the
estimation error in \(U_i\).

Design B rejected in 29% and 39% of two simulated settings where its null was
true. It is therefore not used as a confirmatory test.

## Scope of the theoretical result

The working theory reduces the problem to a central limit theorem for
independent subject scores and remainder bounds for the cross-fitted nuisance
estimates. The required assumptions include:

- independent subjects;
- controlled cluster size or a cluster-level Lindeberg condition;
- population conditional-mean centering;
- a non-degenerate subject-score variance;
- subject-level cross-fitting;
- a product-rate condition for the two nuisance regressions; and
- consistency of the subject-score variance estimator.

These are proposed results with a proof programme. The repository does not
claim that every technical lemma has been completed.

The scalar score tests a residual-covariance moment. Conditional independence
implies that this moment is zero. The converse is not generally true.
