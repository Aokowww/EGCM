#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "configs/adni_simulation_stage2b_stress_smoke.yaml"
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required.", call. = FALSE)
}

source("code/adni/cluster_gcm_core.R")
source("code/simulation/adni_stress_dgp.R")
config <- yaml::read_yaml(config_path)
sim <- config$simulation
scenarios <- sim$scenarios
if (is.null(scenarios) || length(scenarios) == 0L ||
    is.null(names(scenarios)) || any(names(scenarios) == "")) {
  stop("simulation.scenarios must be a named non-empty mapping.",
       call. = FALSE)
}

workers <- as.integer(stress_value(sim, "workers", 1L))
max_tasks_per_run <- stress_value(sim, "max_tasks_per_run", Inf)
if (!is.infinite(max_tasks_per_run)) {
  max_tasks_per_run <- as.integer(max_tasks_per_run)
}
if (!is.finite(workers) || workers < 1L) {
  stop("simulation.workers must be a positive integer.", call. = FALSE)
}
if (!is.infinite(max_tasks_per_run) &&
    (!is.finite(max_tasks_per_run) || max_tasks_per_run < 1L)) {
  stop("simulation.max_tasks_per_run must be positive.", call. = FALSE)
}

dir.create(config$output$directory, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(config$output$directory, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

run_one <- function(scenario_name, repetition) {
  scenario <- scenarios[[scenario_name]]
  data <- generate_adni_stress_data(
    scenario_name, scenario, sim, repetition
  )
  result <- run_design_a(
    data = data,
    x_markers = c("X1", "X2", "X3"),
    y_outcome = "Y",
    rhs = as.character(stress_value(
      scenario, "fitted_rhs", "s(z, k = 4) + s(time_years, k = 4)"
    )),
    id_col = "RID",
    subject_folds = as.integer(sim$subject_folds),
    bootstrap_reps = as.integer(sim$bootstrap_reps),
    seed = as.integer(sim$seed) + repetition,
    equal_subject_nuisance_weights = TRUE
  )
  visit_counts <- table(data$RID)
  data.frame(
    scenario = scenario_name,
    repetition = repetition,
    null_true = TRUE,
    subjects = nrow(result$scores),
    mean_visits = mean(visit_counts),
    min_visits = min(visit_counts),
    max_visits = max(visit_counts),
    global_statistic = result$global_test$statistic,
    global_p_value = result$global_test$p_value
  )
}

tasks <- expand.grid(
  scenario = names(scenarios),
  repetition = seq_len(as.integer(sim$repetitions)),
  stringsAsFactors = FALSE
)
tasks$checkpoint <- file.path(
  checkpoint_dir,
  sprintf("%s_%04d.csv", tasks$scenario, tasks$repetition)
)

required_columns <- c(
  "scenario", "repetition", "null_true", "subjects", "mean_visits",
  "min_visits", "max_visits", "global_statistic", "global_p_value"
)
valid_checkpoint <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0L) return(FALSE)
  value <- try(utils::read.csv(path, stringsAsFactors = FALSE), silent = TRUE)
  !inherits(value, "try-error") &&
    nrow(value) == 1L &&
    all(required_columns %in% names(value)) &&
    is.finite(value$global_p_value[[1L]]) &&
    value$global_p_value[[1L]] > 0 &&
    value$global_p_value[[1L]] <= 1
}

completed_before <- vapply(tasks$checkpoint, valid_checkpoint, logical(1))
pending <- which(!completed_before)
if (length(pending) > max_tasks_per_run) {
  pending <- pending[seq_len(max_tasks_per_run)]
}
message(
  "Starting ", length(pending), " pending stress repetitions with ", workers,
  " worker(s); ", sum(completed_before), "/", nrow(tasks),
  " already checkpointed."
)

task_runner <- function(index) {
  result <- run_one(tasks$scenario[[index]], tasks$repetition[[index]])
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
    sum(failed), " stress task(s) failed. First error: ",
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
  stop("No valid stress checkpoints are available.", call. = FALSE)
}

gate_minimum <- as.integer(config$gate$minimum_repetitions)
gate_maximum <- as.numeric(config$gate$max_exact_95_high)
summary_rows <- do.call(
  rbind,
  lapply(split(results, results$scenario), function(part) {
    rejected <- part$global_p_value < as.numeric(sim$alpha)
    n <- nrow(part)
    rate <- mean(rejected)
    exact <- stats::binom.test(sum(rejected), n)$conf.int
    gate_evaluable <- n >= gate_minimum
    data.frame(
      scenario = part$scenario[[1L]],
      repetitions = n,
      mean_subjects = mean(part$subjects),
      mean_visits = mean(part$mean_visits),
      rejection_rate = rate,
      monte_carlo_se = sqrt(rate * (1 - rate) / n),
      exact_95_low = exact[[1L]],
      exact_95_high = exact[[2L]],
      gate_evaluable = gate_evaluable,
      gate_pass = if (gate_evaluable) exact[[2L]] <= gate_maximum else NA
    )
  })
)
rownames(summary_rows) <- NULL
all_evaluable <- nrow(summary_rows) == length(scenarios) &&
  all(summary_rows$gate_evaluable)
overall_gate <- if (!all_evaluable) {
  "PENDING"
} else if (all(summary_rows$gate_pass)) {
  "PASS"
} else {
  "FAIL"
}

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
    config$stage$label,
    paste0("Scenarios: ", paste(names(scenarios), collapse = ", "), "."),
    paste0(
      "Each scenario has ", sim$repetitions,
      " repetitions; exact binomial 95% intervals are reported."
    ),
    paste0(
      "Checkpoint progress: ", sum(completed_after), "/", nrow(tasks),
      " stress repetitions complete."
    ),
    paste0("Overall stress gate: ", overall_gate, "."),
    paste0(
      "Gate requires at least ", gate_minimum,
      " repetitions per scenario and every exact 95% upper bound <= ",
      gate_maximum, "."
    ),
    "Do not open or interpret ADNI outcome results unless the overall stress gate passes."
  ),
  file.path(config$output$directory, "STATUS.txt")
)

if (all(completed_after)) {
  message("Stress simulation stage completed: ", config$output$directory)
} else {
  message(
    "Stress simulation batch completed: ", sum(completed_after), "/",
    nrow(tasks), ". Re-run the same command to continue from checkpoints."
  )
}
