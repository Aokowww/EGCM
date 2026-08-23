#!/usr/bin/env Rscript

# Pre-specified extensions for the selected plasma--amyloid combination:
# marker-unique GCM tests and held-out prediction gain. The common cohort and
# subject folds are fixed by configs/adni_plasma_amyloid.yaml.

if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required.")
source("code/adni/cluster_gcm_core.R")
config <- yaml::read_yaml("configs/adni_plasma_amyloid.yaml")
data <- read.csv(config$data$path, stringsAsFactors = FALSE, check.names = FALSE)
data$EXAMDATE <- as.Date(data$EXAMDATE)
data <- data[order(data$RID, data$EXAMDATE), ]
data$time_years <- ave(as.numeric(data$EXAMDATE), data$RID,
                       FUN = function(x) (x - min(x)) / 365.25)
data$age_visit <- data$AGE_AT_PET
data$PTGENDER <- factor(data$PTGENDER)
data$TRACER <- factor(data$TRACER)
markers <- config$analysis$x_markers
base_rhs <- config$estimands$cognition_genetics_conditioned$rhs
folds <- make_subject_folds(data$RID, 5L, config$reproducibility$seed)
weights <- make_equal_subject_weights(data$RID)

standardize <- function(x) as.numeric(scale(x))
for (marker in markers) data[[paste0("z_", marker)]] <- standardize(data[[marker]])
data$z_CENTILOIDS <- standardize(data$CENTILOIDS)

unique_rows <- lapply(seq_along(markers), function(index) {
  marker <- markers[[index]]
  other <- setdiff(markers, marker)
  rhs <- paste(c(base_rhs, paste0("z_", other)), collapse = " + ")
  y_fit <- crossfit_residuals(data, "z_CENTILOIDS", rhs, folds,
                              prediction = "population", case_weights = weights,
                              include_subject_re = FALSE)
  x_fit <- crossfit_residuals(data, paste0("z_", marker), rhs, folds,
                              prediction = "population", case_weights = weights,
                              include_subject_re = FALSE)
  scores <- subject_score_matrix(data$RID, matrix(x_fit$residual, ncol = 1L),
                                 y_fit$residual)
  colnames(scores) <- marker
  out <- scalar_cluster_tests(scores)
  out$normalized_residual_covariance <- out$estimate /
    (sd(x_fit$residual) * sd(y_fit$residual))
  out
})
unique_tests <- do.call(rbind, unique_rows)
unique_tests$p_holm <- p.adjust(unique_tests$p_value, method = "holm")

base_y <- crossfit_residuals(data, "z_CENTILOIDS", base_rhs, folds,
                             prediction = "population", case_weights = weights,
                             include_subject_re = FALSE)
full_rhs <- paste(c(base_rhs, paste0("z_", markers)), collapse = " + ")
full_y <- crossfit_residuals(data, "z_CENTILOIDS", full_rhs, folds,
                             prediction = "population", case_weights = weights,
                             include_subject_re = FALSE)
subject_mse <- function(residual) tapply(residual^2, data$RID, mean)
base_mse <- subject_mse(base_y$residual)
full_mse <- subject_mse(full_y$residual)
gain <- base_mse - full_mse
set.seed(config$reproducibility$seed + 20L)
boot <- replicate(9999L, mean(sample(gain, replace = TRUE)))
prediction <- data.frame(
  n_subjects = length(gain),
  reference_rmse = sqrt(mean(base_mse)),
  reference_plus_plasma_rmse = sqrt(mean(full_mse)),
  relative_mse_reduction = mean(gain) / mean(base_mse),
  bootstrap_ci_lower = unname(quantile(boot, .025)),
  bootstrap_ci_upper = unname(quantile(boot, .975)),
  bootstrap_p_one_sided = (1 + sum(boot <= 0)) / (length(boot) + 1)
)

dir.create("results/adni_plasma_amyloid", recursive = TRUE, showWarnings = FALSE)
write.csv(unique_tests,
          "results/adni_plasma_amyloid/unique_marker_tests.csv", row.names = FALSE)
write.csv(prediction,
          "results/adni_plasma_amyloid/heldout_prediction_gain.csv", row.names = FALSE)
message("Plasma--amyloid unique-marker and prediction extensions completed.")
