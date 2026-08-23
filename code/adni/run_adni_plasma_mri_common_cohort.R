#!/usr/bin/env Rscript

# Common-cohort comparison of plasma and structural MRI information about
# amyloid Centiloids. Uses the already frozen plasma-PET cohort and a nearest
# QC-eligible MRI match within the same 90-day window.

source("code/adni/cluster_gcm_core.R")
if (!requireNamespace("yaml", quietly = TRUE)) stop("Package yaml is required.")
config <- yaml::read_yaml("configs/adni_plasma_amyloid.yaml")
d <- read.csv(config$data$path, stringsAsFactors = FALSE, check.names = FALSE)
m <- read.csv("data/adni_raw/2026-07-26/UCSFFSX7_26Jul2026.csv",
              stringsAsFactors = FALSE, check.names = FALSE,
              na.strings = c("", "NA"))
d$EXAMDATE <- as.Date(d$EXAMDATE); m$EXAMDATE <- as.Date(m$EXAMDATE)
mri_cols <- c("ST10CV", "ST29SV", "ST88SV", "ST24CV", "ST83CV",
              "ST37SV", "ST96SV", "ST30SV", "ST89SV", "ST127SV", "ST9SV")
m <- m[complete.cases(m[, c("RID", "EXAMDATE", mri_cols)]), ]
m <- m[apply(m[, mri_cols], 1L, function(x) all(is.finite(x) & x > 0)), ]
m <- m[(is.na(m$OVERALLQC) | m$OVERALLQC != "Fail") &
       (is.na(m$VENTQC) | m$VENTQC != "Fail") &
       (is.na(m$HIPPOQC) | m$HIPPOQC != "Fail"), ]
m <- m[order(m$RID, m$EXAMDATE, m$IMAGEUID), ]
m <- m[!duplicated(m[, c("RID", "EXAMDATE")]), ]
by_id <- split(seq_len(nrow(m)), as.character(m$RID))
selected <- rep(NA_integer_, nrow(d)); gap <- rep(NA_real_, nrow(d))
for (r in seq_len(nrow(d))) {
  rows <- by_id[[as.character(d$RID[[r]])]]
  if (is.null(rows)) next
  g <- as.numeric(m$EXAMDATE[rows] - d$EXAMDATE[[r]])
  ok <- which(is.finite(g) & abs(g) <= 90)
  if (!length(ok)) next
  local <- ok[order(abs(g[ok]), g[ok])][[1L]]
  selected[[r]] <- rows[[local]]; gap[[r]] <- g[[local]]
}
keep <- !is.na(selected)
d <- d[keep, ]; mm <- m[selected[keep], ]
d$ICV <- mm$ST10CV
d$Hippocampus <- mm$ST29SV + mm$ST88SV
d$Entorhinal <- mm$ST24CV + mm$ST83CV
d$Ventricles <- rowSums(mm[, c("ST37SV", "ST96SV", "ST30SV", "ST89SV",
                                 "ST127SV", "ST9SV")])
d$MRI_gap_days <- gap[keep]
d <- d[complete.cases(d[, c("ICV", "Hippocampus", "Entorhinal", "Ventricles")]), ]
d <- d[order(d$RID, d$EXAMDATE), ]
d$time_years <- ave(as.numeric(d$EXAMDATE), d$RID,
                    FUN = function(x) (x - min(x)) / 365.25)
d$age_visit <- d$AGE_AT_PET
d$PTGENDER <- factor(d$PTGENDER); d$TRACER <- factor(d$TRACER)
plasma <- c("pTau217", "AB42_40", "NfL", "GFAP")
mri <- c("Ventricles", "Hippocampus", "Entorhinal")
all_x <- c(plasma, mri)
for (v in c(all_x, "CENTILOIDS", "ICV")) d[[paste0("z_", v)]] <- as.numeric(scale(d[[v]]))
base <- config$estimands$cognition_genetics_conditioned$rhs
base <- paste(base, "+ z_ICV + MRI_gap_days")
folds <- make_subject_folds(d$RID, 5L, config$reproducibility$seed)
w <- make_equal_subject_weights(d$RID)
fit_y <- function(rhs) crossfit_residuals(d, "z_CENTILOIDS", rhs, folds,
  prediction = "population", case_weights = w, include_subject_re = FALSE)
subject_mse <- function(residual) tapply(residual^2, d$RID, mean)
reference <- fit_y(base); plasma_fit <- fit_y(paste(c(base, paste0("z_", plasma)), collapse=" + "))
mri_fit <- fit_y(paste(c(base, paste0("z_", mri)), collapse=" + "))
both_fit <- fit_y(paste(c(base, paste0("z_", all_x)), collapse=" + "))
models <- list(reference=reference, plasma=plasma_fit, mri=mri_fit, both=both_fit)
mse <- lapply(models, function(x) subject_mse(x$residual))
pred <- data.frame(
  model = names(mse), n_subjects = vapply(mse, length, integer(1)),
  rmse = vapply(mse, function(x) sqrt(mean(x)), numeric(1)),
  relative_mse_reduction_vs_reference = vapply(mse, function(x) 1-mean(x)/mean(mse$reference), numeric(1))
)
unique_block <- function(block, other) {
  rhs <- paste(c(base, paste0("z_", other)), collapse=" + ")
  y <- fit_y(rhs)
  xfits <- lapply(block, function(v) crossfit_residuals(
    d, paste0("z_", v), rhs, folds, prediction="population", case_weights=w,
    include_subject_re=FALSE))
  xr <- vapply(xfits, function(x) x$residual, numeric(nrow(d))); colnames(xr)<-block
  scores <- subject_score_matrix(d$RID, xr, y$residual)
  global_multiplier_test(scores, 9999L, config$reproducibility$seed + length(block))
}
plasma_unique <- unique_block(plasma, mri); mri_unique <- unique_block(mri, plasma)
block <- data.frame(
  block=c("plasma_given_mri", "mri_given_plasma"), n_subjects=length(unique(d$RID)),
  statistic=c(plasma_unique$statistic, mri_unique$statistic),
  p_value=c(plasma_unique$p_value, mri_unique$p_value)
)
dir.create("results/adni_plasma_mri_common", recursive=TRUE, showWarnings=FALSE)
write.csv(pred, "results/adni_plasma_mri_common/prediction_models.csv", row.names=FALSE)
write.csv(block, "results/adni_plasma_mri_common/unique_block_tests.csv", row.names=FALSE)
write.csv(data.frame(rows=nrow(d), participants=length(unique(d$RID)),
                     median_abs_mri_gap=median(abs(d$MRI_gap_days))),
          "results/adni_plasma_mri_common/cohort.csv", row.names=FALSE)
message("Common-cohort plasma--MRI comparison completed.")
