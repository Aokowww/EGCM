#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
audit_only <- "--audit-only" %in% args
design_a_only <- "--design-a-only" %in% args
known_flags <- c("--audit-only", "--design-a-only")
positional_args <- args[!args %in% known_flags]
config_path <- if (length(positional_args) >= 1L) {
  positional_args[[1L]]
} else {
  "configs/adni_redesign.yaml"
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install the 'yaml' package before running this script.", call. = FALSE)
}

source("code/adni/cluster_gcm_core.R")
config <- yaml::read_yaml(config_path)

sha256_file <- function(path) {
  executable <- Sys.which("shasum")
  if (!nzchar(executable)) {
    warning("shasum is unavailable; input SHA-256 was not recorded.")
    return(NA_character_)
  }
  output <- system2(
    executable,
    args = c("-a", "256", shQuote(normalizePath(path))),
    stdout = TRUE
  )
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

data_path <- config$data$path
if (length(positional_args) >= 2L) {
  data_path <- positional_args[[2L]]
}
if (!file.exists(data_path)) {
  stop(
    "ADNI input was not found at ", data_path,
    ". Download it with your approved ADNI account and keep it under data/.",
    call. = FALSE
  )
}

adni <- utils::read.csv(data_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- unique(unlist(config$data$required_columns, use.names = FALSE))
assert_columns(adni, required)

for (column in names(config$data$expected_ranges)) {
  limits <- unlist(config$data$expected_ranges[[column]], use.names = FALSE)
  values <- adni[[column]]
  if (!is.numeric(values)) {
    stop(
      column, " must be numeric after import. Check the current ADNI dictionary ",
      "and missing-value encoding.",
      call. = FALSE
    )
  }
  invalid <- !is.na(values) & (
    !is.finite(values) | values < limits[[1L]] | values > limits[[2L]]
  )
  if (any(invalid)) {
    stop(
      column, " contains ", sum(invalid),
      " value(s) outside the pre-specified range [",
      limits[[1L]], ", ", limits[[2L]],
      "]. Resolve ADNI missing-value codes or units before analysis.",
      call. = FALSE
    )
  }
}

adni[[config$data$date_column]] <- as.Date(adni[[config$data$date_column]])
if (anyNA(adni[[config$data$date_column]])) {
  stop("EXAMDATE contains missing or unparseable dates.", call. = FALSE)
}

id_col <- config$data$id_column
adni <- adni[order(adni[[id_col]], adni[[config$data$date_column]]), , drop = FALSE]
baseline_date <- ave(
  as.numeric(adni[[config$data$date_column]]),
  adni[[id_col]],
  FUN = min
)
adni$time_years <- (
  as.numeric(adni[[config$data$date_column]]) - baseline_date
) / 365.25
if (!is.null(config$data$age_visit_column)) {
  adni$age_visit <- adni[[config$data$age_visit_column]]
} else {
  adni$age_visit <- adni[[config$data$baseline_age_column]] + adni$time_years
}

analysis_columns <- unique(c(
  required,
  "time_years",
  "age_visit",
  config$analysis$x_markers,
  config$analysis$y_outcome
))
complete <- stats::complete.cases(adni[, analysis_columns, drop = FALSE])
analysis_data <- adni[complete, , drop = FALSE]

if (nrow(analysis_data) == 0L) {
  stop("No complete rows remain after applying the pre-specified analysis set.",
       call. = FALSE)
}

gap_columns <- c(
  "MRI_gap_days", "CDR_gap_days", "MMSE_gap_days", "MOCA_gap_days"
)
gap_columns <- intersect(gap_columns, names(analysis_data))
if (length(gap_columns) > 0L) {
  observed_max_gap <- max(abs(as.matrix(analysis_data[, gap_columns, drop = FALSE])))
  configured_max_gap <- config$data$alignment$primary_max_abs_gap_days
  if (!is.numeric(configured_max_gap) ||
      length(configured_max_gap) != 1L ||
      observed_max_gap > configured_max_gap) {
    stop(
      "Observed source-table alignment exceeds the configured primary window.",
      call. = FALSE
    )
  }
}

collapse_rules <- config$data$factor_collapse_min_subjects
for (column in config$data$factor_columns) {
  value <- as.character(analysis_data[[column]])
  threshold <- if (!is.null(collapse_rules[[column]])) {
    as.integer(collapse_rules[[column]])
  } else {
    0L
  }
  if (threshold > 0L) {
    subjects_per_level <- tapply(
      as.character(analysis_data[[id_col]]),
      value,
      function(subject) length(unique(subject))
    )
    rare <- names(subjects_per_level[subjects_per_level < threshold])
    value[value %in% rare] <- "__OTHER__"
  }
  analysis_data[[column]] <- factor(value)
}

duplicate_key <- duplicated(
  analysis_data[, c(id_col, config$data$date_column), drop = FALSE]
)
if (any(duplicate_key)) {
  stop("Duplicate RID/EXAMDATE rows must be resolved before analysis.",
       call. = FALSE)
}

visit_counts <- table(as.character(analysis_data[[id_col]]))
audit <- data.frame(
  input_file = normalizePath(data_path),
  input_sha256 = sha256_file(data_path),
  input_rows = nrow(adni),
  complete_rows = nrow(analysis_data),
  excluded_rows = sum(!complete),
  n_subjects = length(visit_counts),
  subjects_one_visit = sum(visit_counts == 1L),
  median_visits = stats::median(as.numeric(visit_counts)),
  max_visits = max(visit_counts)
)

dir.create(config$output$directory, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  audit,
  file.path(config$output$directory, "cohort_audit.csv"),
  row.names = FALSE
)

if (audit_only) {
  message(
    "Cohort audit completed without fitting outcome models: ",
    file.path(config$output$directory, "cohort_audit.csv")
  )
  quit(save = "no", status = 0L)
}

if (!isTRUE(config$data$alignment$verified_against_source_tables)) {
  stop(
    "Cohort audit was written, but outcome models were not run. ",
    "Verify MRI/PET/cognitive timing against the current source tables, ",
    "record the chosen date window, then set ",
    "data.alignment.verified_against_source_tables to true in a frozen config.",
    call. = FALSE
  )
}
if (is.null(config$data$alignment$primary_max_abs_gap_days) ||
    !is.numeric(config$data$alignment$primary_max_abs_gap_days) ||
    config$data$alignment$primary_max_abs_gap_days <= 0) {
  stop(
    "Record a positive primary_max_abs_gap_days after verifying source-table ",
    "measurement dates and before outcome modeling.",
    call. = FALSE
  )
}

all_results <- list()
for (estimand_name in names(config$estimands)) {
  rhs <- config$estimands[[estimand_name]]$rhs
  design_a <- run_design_a(
    data = analysis_data,
    x_markers = config$analysis$x_markers,
    y_outcome = config$analysis$y_outcome,
    rhs = rhs,
    id_col = id_col,
    subject_folds = config$cross_fitting$design_a_subject_folds,
    bootstrap_reps = config$inference$bootstrap_reps,
    seed = config$reproducibility$seed
  )

  design_a$marker_tests$estimand <- estimand_name
  design_a$marker_tests$design <- "A"
  marker_results <- design_a$marker_tests

  if (!design_a_only) {
    design_b <- run_design_b(
      data = analysis_data,
      x_markers = config$analysis$x_markers,
      y_outcome = config$analysis$y_outcome,
      rhs = rhs,
      id_col = id_col,
      visit_folds = config$cross_fitting$design_b_visit_folds,
      min_visits = config$cross_fitting$design_b_min_visits,
      bootstrap_reps = config$inference$bootstrap_reps,
      seed = config$reproducibility$seed
    )
    design_b$marker_tests$estimand <- estimand_name
    design_b$marker_tests$design <- "B"
    marker_results <- rbind(marker_results, design_b$marker_tests)
  }

  utils::write.csv(
    marker_results,
    file.path(
      config$output$directory,
      paste0("marker_tests_", estimand_name, ".csv")
    ),
    row.names = FALSE
  )

  global_results <- data.frame(
    estimand = estimand_name,
    design = "A",
    statistic = design_a$global_test$statistic,
    p_value = design_a$global_test$p_value,
    interpretation = design_a$validity
  )
  if (!design_a_only) {
    global_results <- rbind(
      global_results,
      data.frame(
        estimand = estimand_name,
        design = "B",
        statistic = design_b$global_test$statistic,
        p_value = design_b$global_test$p_value,
        interpretation = design_b$validity
      )
    )
  }
  utils::write.csv(
    global_results,
    file.path(
      config$output$directory,
      paste0("global_tests_", estimand_name, ".csv")
    ),
    row.names = FALSE
  )
  utils::write.csv(
    transform(
      design_a$diagnostics$outcomes,
      estimand = estimand_name,
      design = "A"
    ),
    file.path(
      config$output$directory,
      paste0("nuisance_diagnostics_A_", estimand_name, ".csv")
    ),
    row.names = FALSE
  )
  if (!design_a_only) {
    utils::write.csv(
      transform(
        design_b$diagnostics$outcomes,
        estimand = estimand_name,
        design = "B"
      ),
      file.path(
        config$output$directory,
        paste0("nuisance_diagnostics_B_", estimand_name, ".csv")
      ),
      row.names = FALSE
    )
  }

  all_results[[estimand_name]] <- if (design_a_only) {
    list(A = design_a)
  } else {
    list(A = design_a, B = design_b)
  }
}

saveRDS(
  list(
    config = config,
    audit = audit,
    execution = list(design_a_only = design_a_only),
    results = all_results
  ),
  file.path(config$output$directory, "adni_redesign_full_results.rds")
)

message("ADNI redesign completed. Results: ", config$output$directory)
