#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
result_path <- if (length(args) >= 1L) args[[1L]] else {
  "results/adni_mri_adas13_exploratory/adni_redesign_full_results.rds"
}
output_path <- if (length(args) >= 2L) args[[2L]] else {
  "results/adni_mri_adas13_exploratory/effect_and_influence_diagnostics.csv"
}

result <- readRDS(result_path)
estimand_name <- "baseline_and_technical_covariates"
expected_sign <- c(Ventricles = 1, Hippocampus = -1, Entorhinal = -1)

rows <- lapply(c("A", "B"), function(design_name) {
  fit <- result$results[[estimand_name]][[design_name]]
  residual_sd_x <- apply(fit$residuals_x, 2L, stats::sd)
  residual_sd_y <- stats::sd(fit$residuals_y)
  estimate <- colMeans(fit$scores)
  leave_one_out <- vapply(
    seq_len(nrow(fit$scores)),
    function(index) colMeans(fit$scores[-index, , drop = FALSE]),
    numeric(ncol(fit$scores))
  )
  standardized_scores <- scale(fit$scores)

  data.frame(
    design = design_name,
    marker = colnames(fit$scores),
    normalized_residual_covariance = estimate / (residual_sd_x * residual_sd_y),
    leave_one_out_estimate_min = apply(leave_one_out, 1L, min),
    leave_one_out_estimate_max = apply(leave_one_out, 1L, max),
    leave_one_out_direction_stable = vapply(
      seq_along(estimate),
      function(index) all(sign(leave_one_out[index, ]) == sign(estimate[[index]])),
      logical(1L)
    ),
    max_abs_standardized_subject_score = apply(
      abs(standardized_scores), 2L, max
    ),
    matches_medical_direction =
      sign(estimate) == expected_sign[colnames(fit$scores)],
    stringsAsFactors = FALSE
  )
})

diagnostics <- do.call(rbind, rows)
utils::write.csv(diagnostics, output_path, row.names = FALSE)
message("Wrote effect and influence diagnostics: ", output_path)
