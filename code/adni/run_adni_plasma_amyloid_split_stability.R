#!/usr/bin/env Rscript

# Repeated subject-level split analysis for the plasma--amyloid application.
# The script writes aggregate split-level statistics only. It never writes RID,
# dates, fold membership, residuals, or row-level predictions.

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else
  "configs/adni_plasma_split_stability.yaml"

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required.", call. = FALSE)
}
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required.", call. = FALSE)
}

config <- yaml::read_yaml(config_path)
base_config <- yaml::read_yaml(config$base_config)
spec <- config$analysis
source("code/adni/cluster_gcm_core.R")

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
default_parent <- Sys.getenv(
  "ADNI_RESTRICTED_ROOT", unset = file.path(repo_root, "restricted_adni")
)
default_output <- file.path(
  default_parent, "aggregate", config$outputs$default_subdirectory
)
output_dir <- Sys.getenv(
  config$outputs$environment_variable, unset = default_output
)
dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
output_abs <- file.path(
  normalizePath(dirname(output_dir), winslash = "/", mustWork = TRUE),
  basename(output_dir)
)
if (startsWith(output_abs, paste0(repo_root, "/")) ||
    identical(output_abs, repo_root)) {
  stop("Aggregate output must be written outside the public repository.",
       call. = FALSE)
}
dir.create(output_abs, recursive = TRUE, showWarnings = FALSE)

data_path <- base_config$data$path
if (!file.exists(data_path)) {
  stop(
    "The restricted analytic CSV is unavailable at the configured path. ",
    "Reconstruct it locally before running this analysis.", call. = FALSE
  )
}

data <- utils::read.csv(
  data_path, stringsAsFactors = FALSE, check.names = FALSE
)
data$EXAMDATE <- as.Date(data$EXAMDATE)
data <- data[order(data$RID, data$EXAMDATE), , drop = FALSE]
data$time_years <- ave(
  as.numeric(data$EXAMDATE), data$RID,
  FUN = function(x) (x - min(x)) / 365.25
)
data$age_visit <- data$AGE_AT_PET
data$PTGENDER <- factor(data$PTGENDER)
data$TRACER <- factor(data$TRACER)
markers <- unlist(base_config$analysis$x_markers, use.names = FALSE)
base_rhs <- base_config$estimands$cognition_genetics_conditioned$rhs
weights <- make_equal_subject_weights(data$RID)

standardize <- function(x) as.numeric(scale(x))
for (marker in markers) {
  data[[paste0("z_", marker)]] <- standardize(data[[marker]])
}
data$z_CENTILOIDS <- standardize(data$CENTILOIDS)

analyse_split <- function(split_index) {
  seed <- as.integer(spec$seed) + split_index
  folds <- make_subject_folds(
    data$RID, as.integer(spec$subject_folds), seed
  )

  marker_rows <- lapply(seq_along(markers), function(index) {
    marker <- markers[[index]]
    other <- setdiff(markers, marker)
    rhs <- paste(c(base_rhs, paste0("z_", other)), collapse = " + ")
    y_fit <- crossfit_residuals(
      data, "z_CENTILOIDS", rhs, folds,
      prediction = "population", case_weights = weights,
      include_subject_re = FALSE
    )
    x_fit <- crossfit_residuals(
      data, paste0("z_", marker), rhs, folds,
      prediction = "population", case_weights = weights,
      include_subject_re = FALSE
    )
    scores <- subject_score_matrix(
      data$RID, matrix(x_fit$residual, ncol = 1L), y_fit$residual
    )
    colnames(scores) <- marker
    out <- scalar_cluster_tests(scores)
    data.frame(
      split = split_index,
      seed = seed,
      marker = marker,
      normalized_residual_covariance = out$estimate /
        (stats::sd(x_fit$residual) * stats::sd(y_fit$residual)),
      p_value = out$p_value,
      stringsAsFactors = FALSE
    )
  })
  marker_rows <- do.call(rbind, marker_rows)
  marker_rows$p_holm <- stats::p.adjust(marker_rows$p_value, method = "holm")

  reference <- crossfit_residuals(
    data, "z_CENTILOIDS", base_rhs, folds,
    prediction = "population", case_weights = weights,
    include_subject_re = FALSE
  )
  full_rhs <- paste(
    c(base_rhs, paste0("z_", markers)), collapse = " + "
  )
  full <- crossfit_residuals(
    data, "z_CENTILOIDS", full_rhs, folds,
    prediction = "population", case_weights = weights,
    include_subject_re = FALSE
  )
  reference_mse <- mean(tapply(reference$residual^2, data$RID, mean))
  full_mse <- mean(tapply(full$residual^2, data$RID, mean))
  prediction <- data.frame(
    split = split_index,
    seed = seed,
    reference_rmse = sqrt(reference_mse),
    reference_plus_plasma_rmse = sqrt(full_mse),
    relative_mse_reduction = (reference_mse - full_mse) / reference_mse,
    stringsAsFactors = FALSE
  )
  list(markers = marker_rows, prediction = prediction)
}

repetitions <- as.integer(spec$split_repetitions)
message("Running ", repetitions, " repeated subject-fold assignments.")
analyses <- lapply(seq_len(repetitions), analyse_split)
marker_results <- do.call(rbind, lapply(analyses, `[[`, "markers"))
prediction_results <- do.call(rbind, lapply(analyses, `[[`, "prediction"))

marker_summary <- do.call(rbind, lapply(
  split(marker_results, marker_results$marker), function(part) {
    median_covariance <- stats::median(part$normalized_residual_covariance)
    data.frame(
      analysis = "plasma_to_amyloid_unique_marker",
      quantity = part$marker[[1L]],
      split_repetitions = nrow(part),
      median_value = median_covariance,
      minimum_value = min(part$normalized_residual_covariance),
      maximum_value = max(part$normalized_residual_covariance),
      sign_agreement_rate = mean(
        sign(part$normalized_residual_covariance) == sign(median_covariance)
      ),
      median_adjusted_p = stats::median(part$p_holm),
      maximum_adjusted_p = max(part$p_holm),
      detection_rate = mean(part$p_holm < 0.05),
      stringsAsFactors = FALSE
    )
  }
))
prediction_summary <- data.frame(
  analysis = "plasma_to_amyloid_prediction",
  quantity = "relative_mse_reduction",
  split_repetitions = nrow(prediction_results),
  median_value = stats::median(prediction_results$relative_mse_reduction),
  minimum_value = min(prediction_results$relative_mse_reduction),
  maximum_value = max(prediction_results$relative_mse_reduction),
  sign_agreement_rate = mean(prediction_results$relative_mse_reduction > 0),
  median_adjusted_p = NA_real_,
  maximum_adjusted_p = NA_real_,
  detection_rate = NA_real_,
  stringsAsFactors = FALSE
)
summary <- rbind(marker_summary, prediction_summary)

utils::write.csv(
  marker_results, file.path(output_abs, "split_marker_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  prediction_results, file.path(output_abs, "split_prediction_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary, file.path(output_abs, "summary.csv"), row.names = FALSE
)
writeLines(c(
  "Repeated subject-fold stability analysis for plasma--amyloid.",
  paste0("Configuration: ", normalizePath(config_path)),
  paste0("Repeated splits: ", repetitions, "."),
  "Outputs contain aggregate statistics only; no participant rows or fold memberships.",
  "Repeated folds are a sensitivity analysis, not independent replications."
), file.path(output_abs, "STATUS.txt"))

message("Completed aggregate-only split-stability analysis: ", output_abs)
