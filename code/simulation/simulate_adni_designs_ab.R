#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "configs/adni_simulation_smoke.yaml"
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required.", call. = FALSE)
}
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.", call. = FALSE)
}

source("code/adni/cluster_gcm_core.R")
config <- yaml::read_yaml(config_path)
sim <- config$simulation
canonical_regimes <- c(
  "A_null", "A_alt_B_null", "both_alt", "sparse_B_stress"
)
regimes <- if (!is.null(sim$regimes)) {
  unlist(sim$regimes, use.names = FALSE)
} else {
  canonical_regimes
}
designs <- if (!is.null(sim$designs)) {
  toupper(unlist(sim$designs, use.names = FALSE))
} else {
  c("A", "B")
}
if (length(regimes) == 0L || any(!regimes %in% canonical_regimes)) {
  stop("simulation.regimes contains an unsupported regime.", call. = FALSE)
}
if (length(designs) == 0L || any(!designs %in% c("A", "B"))) {
  stop("simulation.designs must contain A, B, or both.", call. = FALSE)
}
regimes <- unique(regimes)
designs <- unique(designs)
workers <- if (!is.null(sim$workers)) as.integer(sim$workers) else 1L
if (!is.finite(workers) || workers < 1L) {
  stop("simulation.workers must be a positive integer.", call. = FALSE)
}
max_tasks_per_run <- if (!is.null(sim$max_tasks_per_run)) {
  as.integer(sim$max_tasks_per_run)
} else {
  Inf
}
if (!is.infinite(max_tasks_per_run) &&
    (!is.finite(max_tasks_per_run) || max_tasks_per_run < 1L)) {
  stop("simulation.max_tasks_per_run must be positive.", call. = FALSE)
}
dir.create(config$output$directory, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(config$output$directory, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

generate_data <- function(regime, repetition) {
  set.seed(sim$seed + repetition + 1000L * match(
    regime,
    canonical_regimes
  ))

  n_subjects <- sim$subjects
  visits <- if (regime == "sparse_B_stress") {
    sample(
      seq.int(sim$visits_sparse_min, sim$visits_sparse_max),
      n_subjects,
      replace = TRUE
    )
  } else {
    rep(sim$visits_regular, n_subjects)
  }
  id <- rep(seq_len(n_subjects), visits)
  visit <- sequence(visits)
  time <- (visit - 1) / pmax(visits[id] - 1, 1)
  z_subject <- stats::rnorm(n_subjects)
  z <- 0.6 * z_subject[id] + stats::rnorm(length(id), sd = 0.8)

  re_correlation <- if (regime == "A_null") {
    0
  } else {
    sim$random_effect_correlation
  }
  re_sigma <- matrix(
    c(1, re_correlation, re_correlation, 1),
    nrow = 2L
  ) * sim$random_effect_sd^2
  re <- MASS::mvrnorm(n_subjects, mu = c(0, 0), Sigma = re_sigma)

  error_x <- matrix(stats::rnorm(length(id) * 3L), ncol = 3L)
  error_y <- stats::rnorm(length(id))
  if (regime == "both_alt") {
    error_y <- error_y + sim$alternative_error_coefficient *
      rowMeans(error_x)
  }

  data.frame(
    RID = id,
    visit = visit,
    time_years = time,
    z = z,
    X1 = sin(z) + 0.3 * time + re[id, 1L] + error_x[, 1L],
    X2 = z^2 - 1 + 0.2 * time + re[id, 1L] + error_x[, 2L],
    X3 = tanh(z) - 0.2 * time + re[id, 1L] + error_x[, 3L],
    Y = cos(z) + 0.4 * time + re[id, 2L] + error_y
  )
}

run_one <- function(regime, repetition) {
  data <- generate_data(regime, repetition)
  rhs <- "s(z, k = 4) + s(time_years, k = 4)"
  seed <- sim$seed + repetition

  rows <- list()
  if ("A" %in% designs) {
    result_a <- run_design_a(
      data = data,
      x_markers = c("X1", "X2", "X3"),
      y_outcome = "Y",
      rhs = rhs,
      id_col = "RID",
      subject_folds = sim$subject_folds,
      bootstrap_reps = sim$bootstrap_reps,
      seed = seed
    )
    rows[["A"]] <- data.frame(
      regime = regime,
      repetition = repetition,
      design = "A",
      null_true = regime == "A_null",
      subjects = nrow(result_a$scores),
      global_statistic = result_a$global_test$statistic,
      global_p_value = result_a$global_test$p_value
    )
  }
  if ("B" %in% designs) {
    result_b <- run_design_b(
      data = data,
      x_markers = c("X1", "X2", "X3"),
      y_outcome = "Y",
      rhs = rhs,
      id_col = "RID",
      visit_folds = sim$visit_folds,
      min_visits = sim$design_b_min_visits,
      bootstrap_reps = sim$bootstrap_reps,
      seed = seed
    )
    rows[["B"]] <- data.frame(
      regime = regime,
      repetition = repetition,
      design = "B",
      null_true = regime %in%
        c("A_null", "A_alt_B_null", "sparse_B_stress"),
      subjects = nrow(result_b$scores),
      global_statistic = result_b$global_test$statistic,
      global_p_value = result_b$global_test$p_value
    )
  }
  do.call(rbind, rows)
}

tasks <- expand.grid(
  regime = regimes,
  repetition = seq_len(sim$repetitions),
  stringsAsFactors = FALSE
)
tasks$checkpoint <- file.path(
  checkpoint_dir,
  sprintf("%s_%04d.csv", tasks$regime, tasks$repetition)
)
valid_checkpoint <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0L) return(FALSE)
  value <- try(utils::read.csv(path, stringsAsFactors = FALSE), silent = TRUE)
  !inherits(value, "try-error") &&
    nrow(value) == length(designs) &&
    setequal(as.character(value$design), designs) &&
    all(c(
      "regime", "repetition", "design", "null_true", "subjects",
      "global_statistic", "global_p_value"
    ) %in% names(value))
}
completed_before <- vapply(tasks$checkpoint, valid_checkpoint, logical(1))
pending <- which(!completed_before)
if (length(pending) > max_tasks_per_run) {
  pending <- pending[seq_len(max_tasks_per_run)]
}
message(
  "Starting ", length(pending), " pending regime-repetitions with ", workers,
  " worker(s); ", sum(completed_before), "/", nrow(tasks),
  " already checkpointed."
)
task_runner <- function(index) {
  result <- run_one(tasks$regime[[index]], tasks$repetition[[index]])
  temporary <- tempfile(
    pattern = "checkpoint_",
    tmpdir = checkpoint_dir,
    fileext = ".csv"
  )
  utils::write.csv(result, temporary, row.names = FALSE)
  if (!file.rename(temporary, tasks$checkpoint[[index]])) {
    unlink(temporary)
    stop("Could not atomically save checkpoint: ", tasks$checkpoint[[index]])
  }
  result
}
task_results <- if (length(pending) == 0L) {
  list()
} else if (workers == 1L) {
  lapply(pending, task_runner)
} else {
  parallel::mclapply(
    pending,
    task_runner,
    mc.cores = workers,
    mc.preschedule = FALSE
  )
}
failed <- vapply(task_results, inherits, logical(1), "try-error")
if (length(failed) > 0L && any(failed)) {
  stop(
    sum(failed), " simulation task(s) failed. First error: ",
    as.character(task_results[[which(failed)[[1L]]]]),
    call. = FALSE
  )
}
completed_after <- vapply(tasks$checkpoint, valid_checkpoint, logical(1))
checkpoint_results <- lapply(
  tasks$checkpoint[completed_after],
  utils::read.csv,
  stringsAsFactors = FALSE
)
results <- do.call(rbind, checkpoint_results)
if (is.null(results) || nrow(results) == 0L) {
  stop("No valid simulation checkpoints are available.", call. = FALSE)
}

summary_rows <- do.call(
  rbind,
  lapply(split(results, list(results$regime, results$design)), function(part) {
    rejected <- part$global_p_value < sim$alpha
    n <- nrow(part)
    rate <- mean(rejected)
    se <- sqrt(rate * (1 - rate) / n)
    exact <- stats::binom.test(sum(rejected), n)$conf.int
    gate_minimum <- if (!is.null(config$gate$minimum_repetitions)) {
      as.integer(config$gate$minimum_repetitions)
    } else {
      NA_integer_
    }
    gate_maximum <- if (!is.null(config$gate$max_exact_95_high)) {
      as.numeric(config$gate$max_exact_95_high)
    } else {
      NA_real_
    }
    gate_evaluable <- isTRUE(part$null_true[[1L]]) &&
      is.finite(gate_minimum) &&
      n >= gate_minimum &&
      is.finite(gate_maximum)
    data.frame(
      regime = part$regime[[1L]],
      design = part$design[[1L]],
      null_true = part$null_true[[1L]],
      repetitions = n,
      mean_subjects = mean(part$subjects),
      rejection_rate = rate,
      monte_carlo_se = se,
      exact_95_low = exact[[1L]],
      exact_95_high = exact[[2L]],
      gate_evaluable = gate_evaluable,
      gate_pass = if (gate_evaluable) exact[[2L]] <= gate_maximum else NA
    )
  })
)
rownames(summary_rows) <- NULL

utils::write.csv(
  results,
  file.path(config$output$directory, "replicate_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_rows,
  file.path(config$output$directory, "summary.csv"),
  row.names = FALSE
)
writeLines(
  c(
    if (!is.null(config$stage$label)) config$stage$label else
      "Monte Carlo calibration run.",
    paste0(
      "Each regime has ", sim$repetitions,
      " repetitions; exact binomial 95% intervals are reported."
    ),
    paste0(
      "Checkpoint progress: ", sum(completed_after), "/", nrow(tasks),
      " regime-repetitions complete."
    ),
    paste0("Designs: ", paste(designs, collapse = ", "), "."),
    if (!is.null(config$gate$minimum_repetitions)) paste0(
      "Configured gate minimum repetitions: ",
      config$gate$minimum_repetitions, "."
    ),
    if (!is.null(config$gate$max_exact_95_high)) paste0(
      "Configured gate maximum exact 95% upper bound: ",
      config$gate$max_exact_95_high, "."
    ),
    "Do not open or interpret ADNI outcome results until the formal gate passes."
  ),
  file.path(config$output$directory, "STATUS.txt")
)

if (all(completed_after)) {
  message("Simulation stage completed: ", config$output$directory)
} else {
  message(
    "Simulation batch completed: ", sum(completed_after), "/", nrow(tasks),
    ". Re-run the same command to continue from checkpoints."
  )
}
