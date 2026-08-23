#!/usr/bin/env Rscript

# Prospective ADNI Case-C analysis.
#
# Scientific question:
# Do plasma biomarkers contain information about future ADAS-Cog13 beyond
# contemporaneous cognition, demographics, observed APOE4 context, and a
# finite observed cognitive history?
#
# Privacy boundary:
# Participant-level ADNI data and derived rows remain under
# ADNI_RESTRICTED_ROOT. This script writes aggregate counts, tests, and
# cross-validated performance only. It never writes RID, dates, or row-level
# predictions.

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else
  "configs/adni_prospective_ptau_case_c.yaml"

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required.", call. = FALSE)
}
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required.", call. = FALSE)
}

config <- yaml::read_yaml(config_path)
spec <- config$analysis
restricted_root <- Sys.getenv("ADNI_RESTRICTED_ROOT", unset = "")
if (!nzchar(restricted_root) || !dir.exists(restricted_root)) {
  stop("Set ADNI_RESTRICTED_ROOT to the local restricted ADNI directory.",
       call. = FALSE)
}

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
default_output <- file.path(restricted_root, "aggregate", "prospective_ptau_case_c")
output_env <- config$outputs$environment_variable
output_dir <- Sys.getenv(output_env, unset = default_output)
output_parent <- dirname(output_dir)
dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)
output_abs <- normalizePath(output_parent, winslash = "/", mustWork = TRUE)
output_abs <- file.path(output_abs, basename(output_dir))
if (startsWith(output_abs, paste0(repo_root, "/")) || identical(output_abs, repo_root)) {
  stop("Aggregate output must first be written outside the public repository.",
       call. = FALSE)
}
dir.create(output_abs, recursive = TRUE, showWarnings = FALSE)

source("code/adni/cluster_gcm_core.R")

data_dir <- file.path(restricted_root, "ADNIMERGE2", "data")
required_objects <- c(
  "UPENN_PLASMA_FUJIREBIO_QUANTERIX", "ADAS", "ADSL", "APOERES"
)
required_files <- file.path(data_dir, paste0(required_objects, ".rda"))
if (any(!file.exists(required_files))) {
  stop("The restricted ADNI directory is missing a required ADNIMERGE2 table.",
       call. = FALSE)
}

load_object <- function(name) {
  env <- new.env(parent = emptyenv())
  load(file.path(data_dir, paste0(name, ".rda")), envir = env)
  env[[name]]
}

first_nonmissing <- function(values) {
  ok <- !is.na(values) & as.character(values) != ""
  if (any(ok)) values[which(ok)[1L]] else NA
}

count_apoe4 <- function(genotype) {
  genotype <- as.character(genotype)
  out <- nchar(gsub("[^4]", "", genotype))
  out[is.na(genotype) | genotype == ""] <- NA_integer_
  as.integer(out)
}

mean_same_day <- function(data, value) {
  key <- interaction(data$RID, data$VISDATE, drop = TRUE, lex.order = TRUE)
  pieces <- split(seq_len(nrow(data)), key)
  do.call(rbind, lapply(pieces, function(rows) {
    data.frame(
      RID = data$RID[rows[1L]],
      VISDATE = data$VISDATE[rows[1L]],
      TOTAL13 = mean(data[[value]][rows]),
      stringsAsFactors = FALSE
    )
  }))
}

make_demographics <- function(adsl, apoeres) {
  adsl$RID_JOIN <- suppressWarnings(as.integer(as.character(adsl$SUBJID)))
  adsl <- adsl[is.finite(adsl$RID_JOIN), , drop = FALSE]
  adsl <- adsl[order(adsl$RID_JOIN, adsl$ENRLDT), , drop = FALSE]
  pieces <- split(seq_len(nrow(adsl)), adsl$RID_JOIN)
  subject <- do.call(rbind, lapply(pieces, function(rows) {
    d <- adsl[rows, , drop = FALSE]
    data.frame(
      RID = d$RID_JOIN[1L],
      AGE = suppressWarnings(as.numeric(first_nonmissing(d$AGE))),
      SEX = as.character(first_nonmissing(d$SEX)),
      EDUC = suppressWarnings(as.numeric(first_nonmissing(d$EDUC))),
      ENRLDT = as.Date(first_nonmissing(d$ENRLDT)),
      APOE4 = count_apoe4(first_nonmissing(d$APOE)),
      stringsAsFactors = FALSE
    )
  }))

  apoeres <- apoeres[is.finite(apoeres$RID), , drop = FALSE]
  apoe_pieces <- split(seq_len(nrow(apoeres)), apoeres$RID)
  apoe <- do.call(rbind, lapply(apoe_pieces, function(rows) {
    data.frame(
      RID = apoeres$RID[rows[1L]],
      APOE4_RES = count_apoe4(first_nonmissing(apoeres$GENOTYPE[rows])),
      stringsAsFactors = FALSE
    )
  }))
  subject <- merge(subject, apoe, by = "RID", all.x = TRUE, sort = FALSE)
  missing_apoe <- is.na(subject$APOE4)
  subject$APOE4[missing_apoe] <- subject$APOE4_RES[missing_apoe]
  subject$APOE4_RES <- NULL
  subject
}

align_candidate <- function(index_row, cognition, current_window,
                            history_min, history_max,
                            future_min, future_max, target_future) {
  dates <- cognition$VISDATE
  delta <- as.numeric(dates - index_row$EXAMDATE)

  current_ok <- which(is.finite(delta) & abs(delta) <= current_window)
  history_ok <- which(is.finite(delta) & -delta >= history_min &
                        -delta <= history_max)
  future_ok <- which(is.finite(delta) & delta >= future_min &
                       delta <= future_max)
  if (!length(current_ok) || !length(history_ok) || !length(future_ok)) {
    return(NULL)
  }

  current <- current_ok[order(abs(delta[current_ok]), delta[current_ok])][1L]
  prior <- history_ok[order(-delta[history_ok])][1L]
  future <- future_ok[order(abs(delta[future_ok] - target_future),
                            delta[future_ok])][1L]

  history_dates <- as.numeric(dates[history_ok])
  history_values <- cognition$TOTAL13[history_ok]
  history_slope <- NA_real_
  if (length(unique(history_dates)) >= 2L) {
    history_slope <- stats::coef(stats::lm(
      history_values ~ I((history_dates - max(history_dates)) / 365.25)
    ))[[2L]]
  }

  data.frame(
    RID = index_row$RID,
    EXAMDATE = index_row$EXAMDATE,
    pT217_F = index_row$pT217_F,
    AB42_AB40_F = index_row$AB42_AB40_F,
    NfL_Q = index_row$NfL_Q,
    GFAP_Q = index_row$GFAP_Q,
    current_ADAS13 = cognition$TOTAL13[current],
    prior_ADAS13 = cognition$TOTAL13[prior],
    future_ADAS13 = cognition$TOTAL13[future],
    current_gap_days = delta[current],
    history_gap_days = -delta[prior],
    future_gap_days = delta[future],
    history_visits = length(history_ok),
    history_slope = history_slope,
    stringsAsFactors = FALSE
  )
}

prepare_cohort <- function(current_window, plasma, cognition, demographics) {
  cognition_by_id <- split(cognition, cognition$RID)
  rows <- vector("list", nrow(plasma))
  for (r in seq_len(nrow(plasma))) {
    subject_cognition <- cognition_by_id[[as.character(plasma$RID[r])]]
    if (is.null(subject_cognition)) next
    rows[[r]] <- align_candidate(
      plasma[r, , drop = FALSE], subject_cognition, current_window,
      as.numeric(spec$history_min_days), as.numeric(spec$history_max_days),
      as.numeric(spec$future_min_days), as.numeric(spec$future_max_days),
      as.numeric(spec$target_future_days)
    )
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) stop("No eligible prospective rows were found.", call. = FALSE)
  cohort <- do.call(rbind, rows)

  # One pre-specified index per subject: the earliest visit satisfying coverage
  # rules. Selection uses availability and dates, never the future value.
  cohort <- cohort[order(cohort$RID, cohort$EXAMDATE), , drop = FALSE]
  cohort <- cohort[!duplicated(cohort$RID), , drop = FALSE]
  cohort <- merge(cohort, demographics, by = "RID", all.x = TRUE, sort = FALSE)
  cohort$age_index <- cohort$AGE +
    as.numeric(cohort$EXAMDATE - cohort$ENRLDT) / 365.25
  cohort$female <- as.numeric(cohort$SEX == "Female")
  cohort$log_pT217 <- log(cohort$pT217_F)
  cohort$log_AB42_40 <- log(cohort$AB42_AB40_F)
  cohort$log_NfL <- log(cohort$NfL_Q)
  cohort$log_GFAP <- log(cohort$GFAP_Q)

  required <- c(
    "RID", "future_ADAS13", "current_ADAS13", "prior_ADAS13",
    "current_gap_days", "history_gap_days", "future_gap_days",
    "age_index", "female", "EDUC", "APOE4", unlist(spec$markers)
  )
  cohort <- cohort[complete.cases(cohort[, required]), , drop = FALSE]
  cohort <- cohort[is.finite(rowSums(as.matrix(cohort[, required[-1L]]))), , drop = FALSE]
  cohort
}

prediction_row <- function(data, reference_rhs, marker_terms, label,
                           folds, bootstrap_reps, seed) {
  reference <- crossfit_residuals(
    data, "future_ADAS13", reference_rhs, folds,
    prediction = "population", include_subject_re = FALSE
  )
  full_rhs <- paste(c(reference_rhs, marker_terms), collapse = " + ")
  full <- crossfit_residuals(
    data, "future_ADAS13", full_rhs, folds,
    prediction = "population", include_subject_re = FALSE
  )
  ref_loss <- reference$residual^2
  full_loss <- full$residual^2
  gain <- ref_loss - full_loss
  set.seed(seed)
  boot <- replicate(bootstrap_reps, mean(sample(gain, replace = TRUE)))
  data.frame(
    model = label,
    n_subjects = nrow(data),
    reference_rmse = sqrt(mean(ref_loss)),
    reference_plus_x_rmse = sqrt(mean(full_loss)),
    relative_mse_reduction = mean(gain) / mean(ref_loss),
    mean_paired_mse_reduction = mean(gain),
    bootstrap_ci_lower = unname(stats::quantile(boot, 0.025)),
    bootstrap_ci_upper = unname(stats::quantile(boot, 0.975)),
    bootstrap_p_one_sided = (1 + sum(boot <= 0)) / (bootstrap_reps + 1),
    stringsAsFactors = FALSE
  )
}

score_influence <- function(scores) {
  values <- as.numeric(scores[, 1L])
  loo_p <- vapply(seq_along(values), function(index) {
    loo <- matrix(values[-index], ncol = 1L,
                  dimnames = list(NULL, "score"))
    scalar_cluster_tests(loo)$p_value[[1L]]
  }, numeric(1L))
  data.frame(
    max_abs_score_share = max(abs(values)) / sum(abs(values)),
    leave_one_out_p_min = min(loo_p),
    leave_one_out_p_max = max(loo_p),
    stringsAsFactors = FALSE
  )
}

analyse_cohort <- function(data, current_window) {
  seed <- as.integer(spec$seed) + as.integer(current_window)
  folds <- make_subject_folds(data$RID, as.integer(spec$subject_folds), seed)
  marker_names <- unlist(spec$markers, use.names = FALSE)
  marker_terms <- paste0("s(", marker_names, ", k = 4)")
  base_rhs <- paste(c(
    "s(age_index, k = 4)", "female", "s(EDUC, k = 4)", "APOE4",
    "s(current_ADAS13, k = 5)", "current_gap_days", "future_gap_days"
  ), collapse = " + ")
  history_rhs <- paste(c(
    base_rhs, "s(prior_ADAS13, k = 5)", "history_gap_days"
  ), collapse = " + ")

  design <- run_design_a(
    data = data,
    x_markers = marker_names,
    y_outcome = "future_ADAS13",
    rhs = history_rhs,
    id_col = "RID",
    subject_folds = as.integer(spec$subject_folds),
    bootstrap_reps = as.integer(spec$bootstrap_reps),
    seed = seed,
    equal_subject_nuisance_weights = TRUE,
    include_subject_re = FALSE
  )
  marker_tests <- design$marker_tests
  marker_tests$current_window_days <- current_window
  marker_tests$test_type <- "history-conditioned marker-wise"
  marker_tests$normalized_residual_covariance <- marker_tests$estimate /
    (apply(design$residuals_x, 2L, stats::sd) * stats::sd(design$residuals_y))
  marker_influence <- do.call(rbind, lapply(
    seq_len(ncol(design$scores)),
    function(index) score_influence(design$scores[, index, drop = FALSE])
  ))
  marker_tests <- cbind(marker_tests, marker_influence)

  unique_rows <- lapply(seq_along(marker_names), function(index) {
    marker <- marker_names[index]
    other_terms <- marker_terms[-index]
    rhs <- paste(c(history_rhs, other_terms), collapse = " + ")
    y_fit <- crossfit_residuals(
      data, "future_ADAS13", rhs, folds,
      prediction = "population", include_subject_re = FALSE
    )
    x_fit <- crossfit_residuals(
      data, marker, rhs, folds,
      prediction = "population", include_subject_re = FALSE
    )
    scores <- subject_score_matrix(
      data$RID, matrix(x_fit$residual, ncol = 1L,
                       dimnames = list(NULL, marker)), y_fit$residual
    )
    out <- scalar_cluster_tests(scores)
    out$normalized_residual_covariance <- out$estimate /
      (stats::sd(x_fit$residual) * stats::sd(y_fit$residual))
    cbind(out, score_influence(scores))
  })
  unique_tests <- do.call(rbind, unique_rows)
  unique_tests$p_holm <- stats::p.adjust(unique_tests$p_value, method = "holm")
  unique_tests$current_window_days <- current_window
  unique_tests$test_type <- "unique marker given other plasma markers"

  prediction <- do.call(rbind, c(
    list(prediction_row(
      data, base_rhs,
      c("s(prior_ADAS13, k = 5)", "history_gap_days"),
      "current-only reference + observed history",
      folds, as.integer(spec$bootstrap_reps), seed + 50L
    )),
    list(prediction_row(
      data, base_rhs, marker_terms, "current-only reference + all markers",
      folds, as.integer(spec$bootstrap_reps), seed + 100L
    )),
    list(prediction_row(
      data, history_rhs, marker_terms, "history reference + all markers",
      folds, as.integer(spec$bootstrap_reps), seed + 200L
    )),
    lapply(seq_along(marker_names), function(index) prediction_row(
      data, history_rhs, marker_terms[index],
      paste0("history reference + ", marker_names[index]),
      folds, as.integer(spec$bootstrap_reps), seed + 300L + index
    ))
  ))
  prediction$current_window_days <- current_window

  global <- data.frame(
    current_window_days = current_window,
    n_subjects = nrow(data),
    statistic = design$global_test$statistic,
    p_value = design$global_test$p_value,
    bootstrap_reps = design$global_test$bootstrap_reps,
    stringsAsFactors = FALSE
  )

  cohort <- data.frame(
    current_window_days = current_window,
    n_subjects = nrow(data),
    history_slope_available = sum(is.finite(data$history_slope)),
    median_current_gap_days = stats::median(abs(data$current_gap_days)),
    median_history_gap_days = stats::median(data$history_gap_days),
    median_future_gap_days = stats::median(data$future_gap_days),
    p25_future_gap_days = unname(stats::quantile(data$future_gap_days, 0.25)),
    p75_future_gap_days = unname(stats::quantile(data$future_gap_days, 0.75)),
    median_current_ADAS13 = stats::median(data$current_ADAS13),
    median_prior_ADAS13 = stats::median(data$prior_ADAS13),
    median_future_ADAS13 = stats::median(data$future_ADAS13),
    stringsAsFactors = FALSE
  )

  list(
    cohort = cohort,
    global = global,
    marker_tests = marker_tests,
    unique_tests = unique_tests,
    prediction = prediction,
    nuisance = transform(design$diagnostics$outcomes,
                         current_window_days = current_window)
  )
}

analyse_primary_split <- function(data, current_window, split_index) {
  seed <- as.integer(spec$seed) + as.integer(current_window) +
    10000L + split_index
  folds <- make_subject_folds(
    data$RID, as.integer(spec$subject_folds), seed
  )
  marker_names <- unlist(spec$markers, use.names = FALSE)
  marker_terms <- paste0("s(", marker_names, ", k = 4)")
  base_rhs <- paste(c(
    "s(age_index, k = 4)", "female", "s(EDUC, k = 4)", "APOE4",
    "s(current_ADAS13, k = 5)", "current_gap_days", "future_gap_days"
  ), collapse = " + ")
  history_rhs <- paste(c(
    base_rhs, "s(prior_ADAS13, k = 5)", "history_gap_days"
  ), collapse = " + ")

  target_index <- match("log_pT217", marker_names)
  if (is.na(target_index)) {
    stop("The split-stability analysis requires log_pT217.", call. = FALSE)
  }
  unique_rhs <- paste(
    c(history_rhs, marker_terms[-target_index]), collapse = " + "
  )
  y_fit <- crossfit_residuals(
    data, "future_ADAS13", unique_rhs, folds,
    prediction = "population", include_subject_re = FALSE
  )
  x_fit <- crossfit_residuals(
    data, "log_pT217", unique_rhs, folds,
    prediction = "population", include_subject_re = FALSE
  )
  scores <- subject_score_matrix(
    data$RID,
    matrix(x_fit$residual, ncol = 1L,
           dimnames = list(NULL, "log_pT217")),
    y_fit$residual
  )
  test <- scalar_cluster_tests(scores)

  reference <- crossfit_residuals(
    data, "future_ADAS13", history_rhs, folds,
    prediction = "population", include_subject_re = FALSE
  )
  full_rhs <- paste(c(history_rhs, marker_terms[target_index]), collapse = " + ")
  full <- crossfit_residuals(
    data, "future_ADAS13", full_rhs, folds,
    prediction = "population", include_subject_re = FALSE
  )
  reference_mse <- mean(reference$residual^2)
  full_mse <- mean(full$residual^2)

  data.frame(
    split = split_index,
    seed = seed,
    current_window_days = current_window,
    normalized_residual_covariance = test$estimate /
      (stats::sd(x_fit$residual) * stats::sd(y_fit$residual)),
    p_value = test$p_value,
    p_bonferroni_four_markers = pmin(1, 4 * test$p_value),
    reference_rmse = sqrt(reference_mse),
    reference_plus_ptau217_rmse = sqrt(full_mse),
    relative_mse_reduction = (reference_mse - full_mse) / reference_mse,
    stringsAsFactors = FALSE
  )
}

plasma <- load_object("UPENN_PLASMA_FUJIREBIO_QUANTERIX")
adas <- load_object("ADAS")
adsl <- load_object("ADSL")
apoeres <- load_object("APOERES")

plasma$EXAMDATE <- as.Date(plasma$EXAMDATE)
marker_source <- c("pT217_F", "AB42_AB40_F", "NfL_Q", "GFAP_Q")
plasma <- plasma[complete.cases(plasma[, c("RID", "EXAMDATE", marker_source)]), ]
plasma <- plasma[apply(plasma[, marker_source], 1L, function(x) all(is.finite(x) & x > 0)), ]
plasma <- plasma[order(plasma$RID, plasma$EXAMDATE), ]
plasma <- plasma[!duplicated(plasma[, c("RID", "EXAMDATE")]), ]

adas$VISDATE <- as.Date(adas$VISDATE)
adas <- adas[is.finite(adas$RID) & !is.na(adas$VISDATE) &
               is.finite(adas$TOTAL13), ]
adas <- adas[is.na(adas$HAS_QC_ERROR) |
               adas$HAS_QC_ERROR != "Has QC error", ]
adas <- mean_same_day(adas, "TOTAL13")
adas <- adas[order(adas$RID, adas$VISDATE), ]

demographics <- make_demographics(adsl, apoeres)
windows <- as.integer(unlist(spec$current_windows_days, use.names = FALSE))
cohorts <- lapply(windows, function(window) {
  prepare_cohort(window, plasma, adas, demographics)
})
analyses <- lapply(seq_along(windows), function(index) {
  window <- windows[[index]]
  message("Preparing and analysing current window +/-", window, " days.")
  analyse_cohort(cohorts[[index]], window)
})

primary_window <- as.integer(spec$primary_current_window_days)
primary_index <- match(primary_window, windows)
if (is.na(primary_index)) {
  stop("The primary current window is absent from current_windows_days.",
       call. = FALSE)
}
split_repetitions <- if (is.null(spec$split_stability_repetitions)) {
  0L
} else {
  as.integer(spec$split_stability_repetitions)
}
split_stability <- split_stability_summary <- NULL
if (split_repetitions > 0L) {
  message(
    "Running ", split_repetitions,
    " repeated subject-fold assignments for the primary C1 analysis."
  )
  split_stability <- do.call(rbind, lapply(
    seq_len(split_repetitions),
    function(index) analyse_primary_split(
      cohorts[[primary_index]], primary_window, index
    )
  ))
  median_covariance <- stats::median(
    split_stability$normalized_residual_covariance
  )
  split_stability_summary <- rbind(
    data.frame(
      analysis = "prospective_ptau217_unique_information",
      quantity = "normalized_residual_covariance",
      split_repetitions = split_repetitions,
      median_value = median_covariance,
      minimum_value = min(split_stability$normalized_residual_covariance),
      maximum_value = max(split_stability$normalized_residual_covariance),
      sign_agreement_rate = mean(
        sign(split_stability$normalized_residual_covariance) ==
          sign(median_covariance)
      ),
      median_adjusted_p = stats::median(
        split_stability$p_bonferroni_four_markers
      ),
      maximum_adjusted_p = max(
        split_stability$p_bonferroni_four_markers
      ),
      detection_rate = mean(
        split_stability$p_bonferroni_four_markers < 0.05
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis = "prospective_ptau217_prediction",
      quantity = "relative_mse_reduction",
      split_repetitions = split_repetitions,
      median_value = stats::median(split_stability$relative_mse_reduction),
      minimum_value = min(split_stability$relative_mse_reduction),
      maximum_value = max(split_stability$relative_mse_reduction),
      sign_agreement_rate = mean(split_stability$relative_mse_reduction > 0),
      median_adjusted_p = NA_real_,
      maximum_adjusted_p = NA_real_,
      detection_rate = NA_real_,
      stringsAsFactors = FALSE
    )
  )
}

bind_result <- function(name) do.call(rbind, lapply(analyses, `[[`, name))
utils::write.csv(bind_result("cohort"), file.path(output_abs, "cohort_audit.csv"),
                 row.names = FALSE)
utils::write.csv(bind_result("global"), file.path(output_abs, "global_tests.csv"),
                 row.names = FALSE)
utils::write.csv(bind_result("marker_tests"),
                 file.path(output_abs, "marker_tests.csv"), row.names = FALSE)
utils::write.csv(bind_result("unique_tests"),
                 file.path(output_abs, "unique_marker_tests.csv"), row.names = FALSE)
utils::write.csv(bind_result("prediction"),
                 file.path(output_abs, "prediction_comparisons.csv"), row.names = FALSE)
utils::write.csv(bind_result("nuisance"),
                 file.path(output_abs, "nuisance_diagnostics.csv"), row.names = FALSE)
if (split_repetitions > 0L) {
  utils::write.csv(split_stability,
                   file.path(output_abs, "split_stability_results.csv"),
                   row.names = FALSE)
  utils::write.csv(split_stability_summary,
                   file.path(output_abs, "split_stability_summary.csv"),
                   row.names = FALSE)
}

writeLines(c(
  "Prospective plasma biomarker -> future ADAS-Cog13 aggregate analysis.",
  paste0("Configuration: ", normalizePath(config_path)),
  paste0("Primary current-alignment window: +/-",
         spec$primary_current_window_days, " days."),
  "Inference: subject-held-out, one eligible index visit per subject.",
  paste0("Repeated subject-fold sensitivity assignments: ",
         split_repetitions, "."),
  "Outputs contain aggregate statistics only; no participant rows or identifiers.",
  "Interpretation boundary: conditional residual association and held-out prediction, not causality or clinical validation."
), file.path(output_abs, "STATUS.txt"))

message("Completed aggregate-only prospective Case-C analysis: ", output_abs)
