#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else
  "configs/non_iid_bc_gradients_smoke.yaml"

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

generate_case_b <- function(repetition, delta = 0, proxy_measurements = Inf) {
  set.seed(as.integer(sim$seed) + 10000L + repetition)
  n_subjects <- as.integer(sim$subjects)
  visits <- as.integer(sim$visits_b)
  id <- rep(seq_len(n_subjects), each = visits)
  theta <- stats::rnorm(n_subjects)
  proxy_noise <- stats::rnorm(n_subjects)
  z <- 0.35 * theta[id] + stats::rnorm(length(id))
  error_x <- stats::rnorm(length(id))
  error_y <- stats::rnorm(length(id))
  theta_proxy <- if (is.infinite(proxy_measurements)) {
    theta
  } else {
    theta + proxy_noise / sqrt(proxy_measurements)
  }

  data.frame(
    RID = id,
    Z = z,
    theta = theta[id],
    theta_proxy = theta_proxy[id],
    X = 0.80 * z + theta[id] + error_x,
    Y = -0.60 * z + theta[id] + delta * error_x + error_y
  )
}

generate_case_c <- function(repetition, delta = 0) {
  set.seed(as.integer(sim$seed) + 20000L + repetition)
  n_subjects <- as.integer(sim$subjects)
  kept_visits <- as.integer(sim$visits_c)
  total_visits <- kept_visits + as.integer(sim$burn_in_c)
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

proxy_grid <- as.numeric(unlist(sim$proxy_measurements, use.names = FALSE))
power_grid <- as.numeric(unlist(sim$power_coefficients, use.names = FALSE))

scenario_spec <- rbind(
  data.frame(
    scenario = paste0("B_proxy_m", proxy_grid, "_null"),
    experiment = "B proxy-quality null",
    case = "B",
    parameter = "auxiliary_measurements",
    parameter_value = proxy_grid,
    repetitions = as.integer(sim$proxy_repetitions),
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = paste0("B_observed_delta", power_grid),
    experiment = "B observed-context power",
    case = "B",
    parameter = "residual_covariance_coefficient",
    parameter_value = power_grid,
    repetitions = as.integer(sim$power_repetitions),
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = paste0("C_complete_delta", power_grid),
    experiment = "C1 complete-history power",
    case = "C1",
    parameter = "residual_covariance_coefficient",
    parameter_value = power_grid,
    repetitions = as.integer(sim$power_repetitions),
    stringsAsFactors = FALSE
  )
)

run_one <- function(spec, repetition) {
  if (spec$experiment == "B proxy-quality null") {
    data <- generate_case_b(
      repetition, delta = 0,
      proxy_measurements = as.numeric(spec$parameter_value)
    )
    rhs <- c("Z", "theta_proxy")
  } else if (spec$experiment == "B observed-context power") {
    data <- generate_case_b(
      repetition, delta = as.numeric(spec$parameter_value),
      proxy_measurements = Inf
    )
    rhs <- c("Z", "theta")
  } else if (spec$experiment == "C1 complete-history power") {
    data <- generate_case_c(
      repetition, delta = as.numeric(spec$parameter_value)
    )
    rhs <- c("Z", "lagX", "lagY")
  } else {
    stop("Unknown experiment: ", spec$experiment, call. = FALSE)
  }

  value <- crossfit_cluster_test(
    data, rhs,
    as.integer(sim$seed) + 100L * repetition +
      match(spec$scenario, scenario_spec$scenario)
  )
  data.frame(
    scenario = spec$scenario,
    repetition = repetition,
    statistic = unname(value[["statistic"]]),
    p_value = unname(value[["p_value"]]),
    estimate = unname(value[["estimate"]]),
    stringsAsFactors = FALSE
  )
}

message("Running ", nrow(scenario_spec), " frozen gradient scenarios.")
pieces <- lapply(seq_len(nrow(scenario_spec)), function(index) {
  spec <- scenario_spec[index, , drop = FALSE]
  repetitions <- seq_len(spec$repetitions[[1L]])
  worker <- function(rep) run_one(spec, rep)
  rows <- if (as.integer(sim$workers) > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      repetitions, worker, mc.cores = as.integer(sim$workers),
      mc.preschedule = TRUE
    )
  } else {
    lapply(repetitions, worker)
  }
  do.call(rbind, rows)
})
results <- do.call(rbind, pieces)
results <- merge(results, scenario_spec, by = "scenario", sort = FALSE)
results$reject <- results$p_value < alpha

summaries <- do.call(rbind, lapply(split(results, results$scenario), function(part) {
  n <- nrow(part)
  rejected <- sum(part$reject)
  exact <- stats::binom.test(rejected, n)$conf.int
  data.frame(
    scenario = part$scenario[[1L]],
    experiment = part$experiment[[1L]],
    case = part$case[[1L]],
    parameter = part$parameter[[1L]],
    parameter_value = part$parameter_value[[1L]],
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

utils::write.csv(
  results, file.path(output_dir, "replicate_results.csv"), row.names = FALSE
)
utils::write.csv(
  summaries, file.path(output_dir, "summary.csv"), row.names = FALSE
)
writeLines(c(
  config$stage$label,
  paste0("Configuration: ", normalizePath(config_path)),
  "Proxy-quality rows are diagnostic generated-context experiments.",
  "Power rows use the scientifically complete observed conditioning set."
), file.path(output_dir, "STATUS.txt"))

message("Completed frozen Case-B/C1 gradient experiment.")
