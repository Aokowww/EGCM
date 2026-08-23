#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else
  "configs/non_iid_bc_simulation_smoke.yaml"

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required.", call. = FALSE)
}

config <- yaml::read_yaml(config_path)
sim <- config$simulation
alpha <- as.numeric(sim$alpha)
output_dir <- config$output$directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

make_subject_folds_fast <- function(id, k, seed) {
  ids <- unique(id)
  set.seed(seed)
  assignment <- sample(rep(seq_len(k), length.out = length(ids)))
  assignment[match(id, ids)]
}

crossfit_cluster_test <- function(data, rhs_terms, seed) {
  folds <- make_subject_folds_fast(
    data$RID, as.integer(sim$subject_folds), seed
  )
  residual_x <- residual_y <- rep(NA_real_, nrow(data))
  rhs <- paste(rhs_terms, collapse = " + ")
  formula_x <- stats::as.formula(paste("X ~", rhs))
  formula_y <- stats::as.formula(paste("Y ~", rhs))

  for (fold in seq_len(max(folds))) {
    train <- folds != fold
    test <- folds == fold
    fit_x <- stats::lm(formula_x, data = data[train, , drop = FALSE])
    fit_y <- stats::lm(formula_y, data = data[train, , drop = FALSE])
    residual_x[test] <- data$X[test] - stats::predict(
      fit_x, newdata = data[test, , drop = FALSE]
    )
    residual_y[test] <- data$Y[test] - stats::predict(
      fit_y, newdata = data[test, , drop = FALSE]
    )
  }

  products <- residual_x * residual_y
  score <- as.numeric(tapply(products, data$RID, mean))
  score_mean <- mean(score)
  score_scale <- sqrt(mean((score - score_mean)^2))
  statistic <- sqrt(length(score)) * score_mean / score_scale
  p_value <- 2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  c(statistic = statistic, p_value = p_value, estimate = score_mean)
}

generate_case_b <- function(repetition, alternative = FALSE) {
  set.seed(as.integer(sim$seed) + 10000L + repetition)
  n_subjects <- as.integer(sim$subjects)
  visits <- as.integer(sim$visits_b)
  id <- rep(seq_len(n_subjects), each = visits)
  theta <- stats::rnorm(n_subjects)
  z <- 0.35 * theta[id] + stats::rnorm(length(id))
  error_x <- stats::rnorm(length(id))
  error_y <- stats::rnorm(length(id))
  delta <- if (alternative) as.numeric(sim$alternative_coefficient_b) else 0

  sparse_r <- as.numeric(sim$proxy_repeats_sparse)
  dense_r <- as.numeric(sim$proxy_repeats_dense)
  theta_sparse <- theta + stats::rnorm(n_subjects, sd = 1 / sqrt(sparse_r))
  theta_dense <- theta + stats::rnorm(n_subjects, sd = 1 / sqrt(dense_r))

  data.frame(
    RID = id,
    Z = z,
    theta = theta[id],
    theta_sparse = theta_sparse[id],
    theta_dense = theta_dense[id],
    X = 0.80 * z + theta[id] + error_x,
    Y = -0.60 * z + theta[id] + delta * error_x + error_y
  )
}

generate_case_c <- function(repetition, alternative = FALSE) {
  set.seed(as.integer(sim$seed) + 20000L + repetition)
  n_subjects <- as.integer(sim$subjects)
  kept_visits <- as.integer(sim$visits_c)
  total_visits <- kept_visits + as.integer(sim$burn_in_c)
  delta <- if (alternative) as.numeric(sim$alternative_coefficient_c) else 0

  rows <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    z <- stats::rnorm(total_visits)
    error_x <- stats::rnorm(total_visits)
    error_y <- stats::rnorm(total_visits)
    x <- y <- numeric(total_visits)
    for (tt in 2:total_visits) {
      x[tt] <- 0.65 * x[tt - 1L] + 0.45 * z[tt] + error_x[tt]
      y[tt] <- 0.55 * y[tt - 1L] + 0.50 * x[tt - 1L] -
        0.35 * z[tt] + delta * error_x[tt] + error_y[tt]
    }
    keep <- seq.int(total_visits - kept_visits + 1L, total_visits)
    rows[[i]] <- data.frame(
      RID = i,
      Z = z[keep],
      lagX = x[keep - 1L],
      lagY = y[keep - 1L],
      X = x[keep],
      Y = y[keep]
    )
  }
  do.call(rbind, rows)
}

scenario_spec <- data.frame(
  scenario = c(
    "B_observed_null",
    "B_proxy_sparse_null",
    "B_proxy_dense_null",
    "B_omit_context_null",
    "B_observed_alternative",
    "C_complete_history_null",
    "C_partial_history_null",
    "C_omit_history_null",
    "C_complete_history_alternative"
  ),
  case = c(rep("B", 5L), rep("C", 4L)),
  target_null = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
  conditioning_status = c(
    "observed context",
    "estimated context: sparse proxy",
    "estimated context: dense proxy",
    "target context omitted",
    "observed context",
    "complete finite history",
    "incomplete finite history",
    "target history omitted",
    "complete finite history"
  ),
  stringsAsFactors = FALSE
)

run_scenario <- function(scenario, repetition) {
  if (scenario == "B_observed_null") {
    data <- generate_case_b(repetition, FALSE)
    rhs <- c("Z", "theta")
  } else if (scenario == "B_proxy_sparse_null") {
    data <- generate_case_b(repetition, FALSE)
    rhs <- c("Z", "theta_sparse")
  } else if (scenario == "B_proxy_dense_null") {
    data <- generate_case_b(repetition, FALSE)
    rhs <- c("Z", "theta_dense")
  } else if (scenario == "B_omit_context_null") {
    data <- generate_case_b(repetition, FALSE)
    rhs <- "Z"
  } else if (scenario == "B_observed_alternative") {
    data <- generate_case_b(repetition, TRUE)
    rhs <- c("Z", "theta")
  } else if (scenario == "C_complete_history_null") {
    data <- generate_case_c(repetition, FALSE)
    rhs <- c("Z", "lagX", "lagY")
  } else if (scenario == "C_partial_history_null") {
    data <- generate_case_c(repetition, FALSE)
    rhs <- c("Z", "lagY")
  } else if (scenario == "C_omit_history_null") {
    data <- generate_case_c(repetition, FALSE)
    rhs <- "Z"
  } else if (scenario == "C_complete_history_alternative") {
    data <- generate_case_c(repetition, TRUE)
    rhs <- c("Z", "lagX", "lagY")
  } else {
    stop("Unknown scenario: ", scenario, call. = FALSE)
  }

  value <- crossfit_cluster_test(
    data, rhs, as.integer(sim$seed) + 100L * repetition
  )
  data.frame(
    scenario = scenario,
    repetition = repetition,
    statistic = unname(value[["statistic"]]),
    p_value = unname(value[["p_value"]]),
    estimate = unname(value[["estimate"]]),
    stringsAsFactors = FALSE
  )
}

message(
  "Running ", nrow(scenario_spec), " scenarios x ", sim$repetitions,
  " repetitions with ", sim$workers, " worker(s)."
)
results <- do.call(rbind, lapply(scenario_spec$scenario, function(scenario) {
  repetitions <- seq_len(as.integer(sim$repetitions))
  worker <- function(rep) run_scenario(scenario, rep)
  pieces <- if (as.integer(sim$workers) > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      repetitions, worker, mc.cores = as.integer(sim$workers),
      mc.preschedule = TRUE
    )
  } else {
    lapply(repetitions, worker)
  }
  do.call(rbind, pieces)
}))
results <- merge(results, scenario_spec, by = "scenario", sort = FALSE)
results$reject <- results$p_value < alpha

summaries <- do.call(rbind, lapply(split(results, results$scenario), function(part) {
  n <- nrow(part)
  rejected <- sum(part$reject)
  exact <- stats::binom.test(rejected, n)$conf.int
  data.frame(
    scenario = part$scenario[[1L]],
    case = part$case[[1L]],
    target_null = part$target_null[[1L]],
    conditioning_status = part$conditioning_status[[1L]],
    repetitions = n,
    rejections = rejected,
    rejection_rate = rejected / n,
    monte_carlo_se = sqrt((rejected / n) * (1 - rejected / n) / n),
    exact_95_low = exact[[1L]],
    exact_95_high = exact[[2L]],
    mean_score = mean(part$estimate),
    stringsAsFactors = FALSE
  )
}))
rownames(summaries) <- NULL
summaries <- summaries[match(scenario_spec$scenario, summaries$scenario), ]

primary <- unlist(config$gate$primary_scenarios, use.names = FALSE)
summaries$primary_gate <- summaries$scenario %in% primary
summaries$gate_evaluable <- summaries$primary_gate &
  summaries$repetitions >= as.integer(config$gate$minimum_repetitions)
summaries$gate_pass <- ifelse(
  summaries$gate_evaluable,
  summaries$exact_95_high <= as.numeric(config$gate$max_exact_95_high),
  NA
)

utils::write.csv(
  results, file.path(output_dir, "replicate_results.csv"), row.names = FALSE
)
utils::write.csv(
  summaries, file.path(output_dir, "summary.csv"), row.names = FALSE
)

plot_labels <- c(
  B_observed_null = "B: observed context (null)",
  B_proxy_sparse_null = "B: sparse proxy (null)",
  B_proxy_dense_null = "B: dense proxy (null)",
  B_omit_context_null = "B: context omitted (target null)",
  B_observed_alternative = "B: observed context (alternative)",
  C_complete_history_null = "C: complete history (null)",
  C_partial_history_null = "C: partial history (target null)",
  C_omit_history_null = "C: history omitted (target null)",
  C_complete_history_alternative = "C: complete history (alternative)"
)
plot_order <- rev(seq_len(nrow(summaries)))
grDevices::pdf(file.path(output_dir, "case_bc_rejection_rates.pdf"),
               width = 8.5, height = 5.8, useDingbats = FALSE)
old_par <- graphics::par(mar = c(4.5, 13.5, 1.0, 0.8))
graphics::plot(
  summaries$rejection_rate[plot_order], seq_along(plot_order),
  xlim = c(0, 1), ylim = c(0.5, length(plot_order) + 0.5),
  xlab = "Empirical rejection rate", ylab = "", yaxt = "n",
  pch = ifelse(summaries$target_null[plot_order], 16, 17),
  col = ifelse(summaries$case[plot_order] == "B", "#0072B2", "#D55E00")
)
graphics::segments(
  summaries$exact_95_low[plot_order], seq_along(plot_order),
  summaries$exact_95_high[plot_order], seq_along(plot_order),
  col = ifelse(summaries$case[plot_order] == "B", "#0072B2", "#D55E00"),
  lwd = 1.4
)
graphics::axis(
  2, at = seq_along(plot_order),
  labels = plot_labels[summaries$scenario[plot_order]], las = 1, cex.axis = 0.78
)
graphics::abline(v = alpha, lty = 2, col = "grey35")
graphics::box()
graphics::par(old_par)
grDevices::dev.off()

gate_rows <- summaries[summaries$primary_gate, , drop = FALSE]
overall_gate <- if (!all(gate_rows$gate_evaluable)) {
  "PENDING"
} else if (all(gate_rows$gate_pass)) {
  "PASS"
} else {
  "FAIL"
}
writeLines(c(
  config$stage$label,
  paste0("Configuration: ", normalizePath(config_path)),
  paste0("Overall primary calibration gate: ", overall_gate, "."),
  paste0(
    "The omitted and incomplete-context rows test reduced conditioning sets; ",
    "their rejection rates are not Type-I-error estimates for those reduced nulls."
  )
), file.path(output_dir, "STATUS.txt"))

message("Completed. Overall primary calibration gate: ", overall_gate)
