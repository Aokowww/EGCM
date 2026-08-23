#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else
  "configs/method_comparison_smoke.yaml"

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required.", call. = FALSE)
}

config <- yaml::read_yaml(config_path)
sim <- config$simulation
alpha <- as.numeric(sim$alpha)
output_dir <- config$output$directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_subject_folds <- function(id, k, seed) {
  ids <- unique(id)
  set.seed(seed)
  assignment <- sample(rep(seq_len(k), length.out = length(ids)))
  assignment[match(id, ids)]
}

ar1_noise <- function(n, rho) {
  innovation <- stats::rnorm(n)
  value <- numeric(n)
  value[[1L]] <- innovation[[1L]] / sqrt(1 - rho^2)
  if (n > 1L) {
    for (tt in 2:n) value[[tt]] <- rho * value[[tt - 1L]] + innovation[[tt]]
  }
  value
}

generate_data <- function(repetition, slope_sd, delta) {
  set.seed(as.integer(sim$seed) + 100000L * round(100 * delta) + repetition)
  n_subjects <- as.integer(sim$subjects)
  visits <- as.integer(sim$visits)
  rho <- as.numeric(sim$ar1)
  intercept_sd <- as.numeric(sim$random_intercept_sd)
  time <- rep(seq(-1, 1, length.out = visits), times = n_subjects)
  id <- rep(seq_len(n_subjects), each = visits)
  z <- stats::rnorm(length(id))

  bx0 <- stats::rnorm(n_subjects, sd = intercept_sd)
  by0 <- stats::rnorm(n_subjects, sd = intercept_sd)
  bx1 <- stats::rnorm(n_subjects, sd = slope_sd)
  by1 <- stats::rnorm(n_subjects, sd = slope_sd)
  ex <- unlist(lapply(seq_len(n_subjects), function(i) ar1_noise(visits, rho)))
  ey <- unlist(lapply(seq_len(n_subjects), function(i) ar1_noise(visits, rho)))

  mean_x <- 0.60 * z + 0.25 * time
  mean_y <- -0.50 * z + 0.15 * time
  x <- mean_x + bx0[id] + bx1[id] * time + ex
  y <- mean_y + by0[id] + by1[id] * time + delta * ex + ey
  data.frame(
    RID = id, time = time, Z = z, X = x, Y = y,
    oracle_residual_x = x - mean_x,
    oracle_residual_y = y - mean_y
  )
}

crossfit_residuals <- function(data, seed) {
  folds <- make_subject_folds(data$RID, as.integer(sim$subject_folds), seed)
  residual_x <- residual_y <- rep(NA_real_, nrow(data))
  for (fold in seq_len(max(folds))) {
    train <- folds != fold
    test <- folds == fold
    fit_x <- stats::lm(X ~ Z + time, data = data[train, , drop = FALSE])
    fit_y <- stats::lm(Y ~ Z + time, data = data[train, , drop = FALSE])
    residual_x[test] <- data$X[test] - stats::predict(
      fit_x, newdata = data[test, , drop = FALSE]
    )
    residual_y[test] <- data$Y[test] - stats::predict(
      fit_y, newdata = data[test, , drop = FALSE]
    )
  }
  list(x = residual_x, y = residual_y)
}

normal_score_test <- function(score) {
  score_mean <- mean(score)
  score_sd <- stats::sd(score)
  statistic <- sqrt(length(score)) * score_mean / score_sd
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  c(statistic = statistic, p_value = p_value, estimate = score_mean)
}

evaluate_methods <- function(repetition, slope_sd, delta, scenario) {
  data <- generate_data(repetition, slope_sd, delta)
  fitted <- crossfit_residuals(
    data, as.integer(sim$seed) + repetition + 1000L * round(100 * slope_sd)
  )
  product <- fitted$x * fitted$y
  oracle_product <- data$oracle_residual_x * data$oracle_residual_y

  subject_score <- as.numeric(tapply(product, data$RID, mean))
  oracle_subject_score <- as.numeric(tapply(oracle_product, data$RID, mean))
  values <- list(
    "subject-score GCM" = normal_score_test(subject_score),
    "visit-iid GCM" = normal_score_test(product),
    "oracle subject score" = normal_score_test(oracle_subject_score)
  )
  do.call(rbind, lapply(names(values), function(method) {
    value <- values[[method]]
    data.frame(
      scenario = scenario,
      delta = delta,
      repetition = repetition,
      method = method,
      statistic = unname(value[["statistic"]]),
      p_value = unname(value[["p_value"]]),
      estimate = unname(value[["estimate"]]),
      stringsAsFactors = FALSE
    )
  }))
}

scenario_spec <- rbind(
  data.frame(
    scenario = c("base_null", "random_slope_null"),
    slope_sd = c(
      as.numeric(sim$random_slope_sd_base),
      as.numeric(sim$random_slope_sd_stress)
    ),
    delta = 0,
    repetitions = as.integer(sim$repetitions_null)
  ),
  data.frame(
    scenario = paste0("alternative_delta_", unlist(sim$alternative_coefficients)),
    slope_sd = as.numeric(sim$random_slope_sd_base),
    delta = as.numeric(unlist(sim$alternative_coefficients)),
    repetitions = as.integer(sim$repetitions_power)
  )
)

message(
  "Running ", sum(scenario_spec$repetitions), " data sets across ",
  nrow(scenario_spec), " scenarios with ", sim$workers, " worker(s)."
)

run_one_scenario <- function(row_index) {
  spec <- scenario_spec[row_index, ]
  repetitions <- seq_len(spec$repetitions)
  worker <- function(rep) evaluate_methods(
    rep, spec$slope_sd, spec$delta, spec$scenario
  )
  pieces <- if (as.integer(sim$workers) > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      repetitions, worker, mc.cores = as.integer(sim$workers),
      mc.preschedule = TRUE
    )
  } else {
    lapply(repetitions, worker)
  }
  do.call(rbind, pieces)
}

results <- do.call(rbind, lapply(seq_len(nrow(scenario_spec)), run_one_scenario))
results$reject <- results$p_value < alpha

summaries <- do.call(
  rbind,
  lapply(split(results, list(results$scenario, results$method), drop = TRUE), function(part) {
    n <- nrow(part)
    rejected <- sum(part$reject)
    exact <- stats::binom.test(rejected, n)$conf.int
    data.frame(
      scenario = part$scenario[[1L]],
      delta = part$delta[[1L]],
      method = part$method[[1L]],
      repetitions = n,
      rejections = rejected,
      rejection_rate = rejected / n,
      monte_carlo_se = sqrt((rejected / n) * (1 - rejected / n) / n),
      exact_95_low = exact[[1L]],
      exact_95_high = exact[[2L]],
      mean_score = mean(part$estimate),
      stringsAsFactors = FALSE
    )
  })
)
rownames(summaries) <- NULL

primary <- config$gate$primary_method
summaries$primary_gate <- summaries$delta == 0 & summaries$method == primary
summaries$gate_evaluable <- summaries$primary_gate &
  summaries$repetitions >= as.integer(config$gate$minimum_null_repetitions)
summaries$gate_pass <- ifelse(
  summaries$gate_evaluable,
  summaries$exact_95_high <= as.numeric(config$gate$maximum_exact_95_high),
  NA
)

utils::write.csv(
  results, file.path(output_dir, "replicate_results.csv"), row.names = FALSE
)
utils::write.csv(
  summaries, file.path(output_dir, "summary.csv"), row.names = FALSE
)

message("Wrote results to ", output_dir)
