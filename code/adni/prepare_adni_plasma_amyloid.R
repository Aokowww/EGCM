#!/usr/bin/env Rscript

# Prepare the PET-anchored plasma--amyloid analysis table. Restricted ADNI
# participant data remain under data/ and are excluded from version control.

args <- commandArgs(trailingOnly = TRUE)
snapshot_new <- if (length(args) >= 1L) args[[1L]] else "data/adni_raw/2026-08-13"
snapshot_context <- if (length(args) >= 2L) args[[2L]] else "data/adni_raw/2026-07-26"
window_days <- if (length(args) >= 3L) as.integer(args[[3L]]) else 90L

plasma_path <- file.path(snapshot_new, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_12Aug2026.csv")
pet_path <- file.path(snapshot_new, "UCBERKELEY_AMY_6MM_12Aug2026.csv")
cdr_path <- file.path(snapshot_context, "CDR_26Jul2026.csv")
mmse_path <- file.path(snapshot_context, "MMSE_26Jul2026.csv")
moca_path <- file.path(snapshot_context, "MOCA_26Jul2026.csv")
adsl_path <- file.path(snapshot_context, "extracted/ADNIMERGE2/data/ADSL.rda")
apoeres_path <- file.path(snapshot_context, "extracted/ADNIMERGE2/data/APOERES.rda")
inputs <- c(plasma_path, pet_path, cdr_path, mmse_path, moca_path, adsl_path, apoeres_path)
if (any(!file.exists(inputs))) stop("Missing ADNI source input.", call. = FALSE)

read_adni <- function(path) utils::read.csv(
  path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")
)
load_object <- function(path, name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env[[name]]
}
first_nonmissing <- function(values) {
  values <- values[!is.na(values) & as.character(values) != ""]
  if (length(values)) values[[1L]] else NA
}
count_apoe4 <- function(genotype) {
  genotype <- as.character(genotype)
  out <- nchar(gsub("[^4]", "", genotype))
  out[is.na(genotype) | genotype == ""] <- NA_integer_
  as.integer(out)
}

nearest <- function(anchor, candidate, candidate_date, value_columns, prefix) {
  by_id <- split(seq_len(nrow(candidate)), as.character(candidate$RID))
  out <- as.data.frame(matrix(NA, nrow(anchor), length(value_columns)))
  names(out) <- paste0(prefix, value_columns)
  out[[paste0(prefix, "date")]] <- as.Date(NA)
  out[[paste0(prefix, "gap_days")]] <- NA_real_
  for (r in seq_len(nrow(anchor))) {
    rows <- by_id[[as.character(anchor$RID[[r]])]]
    if (is.null(rows)) next
    gap <- as.numeric(candidate[[candidate_date]][rows] - anchor$EXAMDATE[[r]])
    ok <- which(is.finite(gap) & abs(gap) <= window_days)
    if (!length(ok)) next
    selected_local <- ok[order(abs(gap[ok]), gap[ok])][[1L]]
    selected <- rows[[selected_local]]
    for (column in value_columns) {
      out[[paste0(prefix, column)]][[r]] <- candidate[[column]][[selected]]
    }
    out[[paste0(prefix, "date")]][[r]] <- candidate[[candidate_date]][[selected]]
    out[[paste0(prefix, "gap_days")]][[r]] <- gap[[selected_local]]
  }
  out
}

plasma <- read_adni(plasma_path)
pet <- read_adni(pet_path)
cdr <- read_adni(cdr_path)
mmse <- read_adni(mmse_path)
moca <- read_adni(moca_path)
adsl <- load_object(adsl_path, "ADSL")
apoeres <- load_object(apoeres_path, "APOERES")

plasma$EXAMDATE <- as.Date(plasma$EXAMDATE)
pet$SCANDATE <- as.Date(pet$SCANDATE)
cdr$VISDATE <- as.Date(cdr$VISDATE)
mmse$VISDATE <- as.Date(mmse$VISDATE)
moca$VISDATE <- as.Date(moca$VISDATE)

# Anchor on PET and remove duplicates without using plasma or cognition.
pet <- pet[complete.cases(pet[, c("RID", "SCANDATE", "CENTILOIDS")]), ]
pet <- pet[is.finite(pet$CENTILOIDS) & (!is.na(pet$qc_flag) & pet$qc_flag >= 1), ]
pet <- pet[order(pet$RID, pet$SCANDATE, pet$PROCESSDATE, decreasing = FALSE), ]
pet <- pet[!duplicated(pet[, c("RID", "SCANDATE")]), ]
anchor <- data.frame(
  RID = as.integer(pet$RID), EXAMDATE = pet$SCANDATE,
  CENTILOIDS = as.numeric(pet$CENTILOIDS), SITE = as.character(pet$SITEID),
  TRACER = as.character(pet$TRACER), stringsAsFactors = FALSE
)

plasma_columns <- c("pT217_F", "AB42_AB40_F", "NfL_Q", "GFAP_Q")
for (column in plasma_columns) {
  plasma[[column]][plasma[[column]] < 0] <- NA_real_
}
plasma <- plasma[complete.cases(plasma[, c("RID", "EXAMDATE", plasma_columns)]), ]
plasma <- plasma[order(plasma$RID, plasma$EXAMDATE), ]
plasma <- plasma[!duplicated(plasma[, c("RID", "EXAMDATE")]), ]

for (object_name in c("cdr", "mmse", "moca")) {
  object <- get(object_name)
  date_column <- "VISDATE"
  object <- object[order(object$RID, object[[date_column]]), ]
  object <- object[!duplicated(object[, c("RID", date_column)]), ]
  assign(object_name, object)
}

aligned <- cbind(
  anchor,
  nearest(anchor, plasma, "EXAMDATE", plasma_columns, "PLASMA_"),
  nearest(anchor, cdr, "VISDATE", "CDRSB", "CDR_"),
  nearest(anchor, mmse, "VISDATE", "MMSCORE", "MMSE_"),
  nearest(anchor, moca, "VISDATE", "MOCA", "MOCA_")
)
aligned$pTau217 <- aligned$PLASMA_pT217_F
aligned$AB42_40 <- aligned$PLASMA_AB42_AB40_F
aligned$NfL <- aligned$PLASMA_NfL_Q
aligned$GFAP <- aligned$PLASMA_GFAP_Q
aligned$CDRSB <- aligned$CDR_CDRSB
aligned$MMSE <- aligned$MMSE_MMSCORE
aligned$MOCA <- aligned$MOCA_MOCA

adsl$RID <- suppressWarnings(as.integer(as.character(adsl$SUBJID)))
subject <- do.call(rbind, lapply(split(adsl, adsl$RID), function(d) data.frame(
  RID = d$RID[[1L]], AGE = as.numeric(first_nonmissing(d$AGE)),
  PTGENDER = as.character(first_nonmissing(d$SEX)),
  PTEDUCAT = as.numeric(first_nonmissing(d$EDUC)),
  ENRLDT = as.Date(first_nonmissing(d$ENRLDT)),
  APOE4 = count_apoe4(first_nonmissing(d$APOE)), stringsAsFactors = FALSE
)))
apoe <- do.call(rbind, lapply(split(apoeres, apoeres$RID), function(d) data.frame(
  RID = d$RID[[1L]], APOE4_RES = count_apoe4(first_nonmissing(d$GENOTYPE))
)))
subject <- merge(subject, apoe, by = "RID", all.x = TRUE)
subject$APOE4[is.na(subject$APOE4)] <- subject$APOE4_RES[is.na(subject$APOE4)]
subject$APOE4_RES <- NULL
aligned <- merge(aligned, subject, by = "RID", all.x = TRUE, sort = FALSE)
aligned$AGE_AT_PET <- aligned$AGE + as.numeric(aligned$EXAMDATE - aligned$ENRLDT) / 365.25
aligned$time_years <- ave(as.numeric(aligned$EXAMDATE), aligned$RID, FUN = function(x) (x - min(x)) / 365.25)

required <- c(
  "RID", "EXAMDATE", "CENTILOIDS", "pTau217", "AB42_40", "NfL", "GFAP",
  "CDRSB", "MMSE", "MOCA", "AGE_AT_PET", "PTGENDER", "PTEDUCAT", "APOE4",
  "SITE", "TRACER", "PLASMA_gap_days", "CDR_gap_days", "MMSE_gap_days", "MOCA_gap_days"
)
complete <- complete.cases(aligned[, required])
analysis <- aligned[complete, ]
analysis <- analysis[order(analysis$RID, analysis$EXAMDATE), ]

dir.create("data", showWarnings = FALSE)
dir.create("results/adni_plasma_amyloid", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(analysis, "data/ADNI_PLASMA_AMYLOID_90d.csv", row.names = FALSE, na = "")
counts <- table(analysis$RID)
audit <- data.frame(
  pet_anchor_rows = nrow(anchor), complete_rows = nrow(analysis),
  complete_subjects = length(counts), subjects_repeated = sum(counts > 1L),
  median_visits = if (length(counts)) median(as.numeric(counts)) else NA_real_,
  median_abs_plasma_gap_days = median(abs(analysis$PLASMA_gap_days)),
  p90_abs_plasma_gap_days = unname(quantile(abs(analysis$PLASMA_gap_days), .9)),
  max_abs_plasma_gap_days = max(abs(analysis$PLASMA_gap_days))
)
utils::write.csv(audit, "results/adni_plasma_amyloid/cohort_audit.csv", row.names = FALSE)
message("Prepared PET-anchored plasma--amyloid table: ", nrow(analysis),
        " rows / ", length(counts), " participants.")
