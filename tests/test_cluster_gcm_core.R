source("code/adni/cluster_gcm_core.R")

set.seed(42)
n_subjects <- 48L
n_visits <- 4L
n <- n_subjects * n_visits
subject <- rep(seq_len(n_subjects), each = n_visits)
z <- stats::rnorm(n)
u_x <- rep(stats::rnorm(n_subjects), each = n_visits)
u_y <- rep(stats::rnorm(n_subjects), each = n_visits)
synthetic <- data.frame(
  RID = subject,
  visit = rep(seq_len(n_visits), n_subjects),
  z = z,
  X1 = sin(z) + u_x + stats::rnorm(n),
  X2 = z^2 + 0.5 * u_x + stats::rnorm(n),
  Y = cos(z) + u_y + stats::rnorm(n)
)

subject_folds <- make_subject_folds(synthetic$RID, k = 4L, seed = 7L)
stopifnot(all(vapply(
  split(subject_folds, synthetic$RID),
  function(x) length(unique(x)) == 1L,
  logical(1)
)))
equal_subject_weights <- make_equal_subject_weights(synthetic$RID)
stopifnot(all(abs(
  tapply(equal_subject_weights, synthetic$RID, sum) - 1
) < 1e-12))

a <- run_design_a(
  synthetic,
  x_markers = c("X1", "X2"),
  y_outcome = "Y",
  rhs = "s(z, k = 5) + visit",
  id_col = "RID",
  subject_folds = 4L,
  bootstrap_reps = 99L,
  seed = 7L
)
stopifnot(
  a$design == "A",
  nrow(a$scores) == n_subjects,
  ncol(a$scores) == 2L,
  a$nuisance_weighting == "equal subject",
  all(is.finite(a$marker_tests$statistic)),
  a$global_test$p_value > 0,
  a$global_test$p_value <= 1
)

a_no_subject_re <- run_design_a(
  synthetic,
  x_markers = c("X1", "X2"),
  y_outcome = "Y",
  rhs = "s(z, k = 5) + visit",
  id_col = "RID",
  subject_folds = 4L,
  bootstrap_reps = 99L,
  seed = 7L,
  include_subject_re = FALSE
)
stopifnot(
  a_no_subject_re$design == "A",
  nrow(a_no_subject_re$scores) == n_subjects,
  all(is.finite(a_no_subject_re$marker_tests$statistic))
)

b <- run_design_b(
  synthetic,
  x_markers = c("X1", "X2"),
  y_outcome = "Y",
  rhs = "s(z, k = 5) + visit",
  id_col = "RID",
  visit_folds = 2L,
  min_visits = 3L,
  bootstrap_reps = 99L,
  seed = 7L
)
stopifnot(
  b$design == "B",
  nrow(b$scores) == n_subjects,
  grepl("sensitivity analysis only", b$validity, fixed = TRUE)
)

message("cluster_gcm_core.R tests passed")
