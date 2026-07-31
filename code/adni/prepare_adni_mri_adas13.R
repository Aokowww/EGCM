#!/usr/bin/env Rscript

# Build an exploratory, date-aligned structural MRI--ADAS13 cohort.
# ADAS13 is the anchor and the nearest QC-eligible FreeSurfer 7 scan is
# selected without using the ADAS13 score or the direction of marker values
# after the pre-specified MRI completeness, physiological-range, and QC checks.

args <- commandArgs(trailingOnly = TRUE)
snapshot_dir <- if (length(args) >= 1L) args[[1L]] else {
  "data/adni_raw/2026-07-26"
}
window_days <- if (length(args) >= 2L) as.integer(args[[2L]]) else 30L
min_visits <- if (length(args) >= 3L) as.integer(args[[3L]]) else 5L

if (!is.finite(window_days) || window_days <= 0L ||
    !is.finite(min_visits) || min_visits < 2L) {
  stop("Require a positive date window and min_visits >= 2.", call. = FALSE)
}

data_dir <- file.path(snapshot_dir, "extracted", "ADNIMERGE2", "data")
paths <- list(
  adas = file.path(data_dir, "ADAS.rda"),
  mri = file.path(data_dir, "UCSFFSX7.rda"),
  adsl = file.path(data_dir, "ADSL.rda"),
  ptdemog = file.path(data_dir, "PTDEMOG.rda"),
  apoeres = file.path(data_dir, "APOERES.rda")
)
missing_files <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]
if (length(missing_files) > 0L) {
  stop("Missing ADNI input(s): ", paste(missing_files, collapse = ", "),
       call. = FALSE)
}

load_rda_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!expected_name %in% loaded) {
    stop("Expected object ", expected_name, " was not found in ", path,
         call. = FALSE)
  }
  as.data.frame(environment[[expected_name]])
}

as_adni_date <- function(x) as.Date(as.character(x))

qc_error <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

first_nonmissing_by_rid <- function(data, columns, date_column = NULL) {
  if (!is.null(date_column)) {
    data <- data[order(data$RID, data[[date_column]], na.last = TRUE), ,
                 drop = FALSE]
  } else {
    data <- data[order(data$RID), , drop = FALSE]
  }
  rows <- split(seq_len(nrow(data)), as.character(data$RID))
  output <- data.frame(RID = as.integer(names(rows)))
  for (column in columns) {
    output[[column]] <- vapply(rows, function(index) {
      values <- data[[column]][index]
      values <- values[!is.na(values) & as.character(values) != ""]
      if (length(values) == 0L) NA_character_ else as.character(values[[1L]])
    }, character(1L))
  }
  output
}

count_apoe4 <- function(genotype) {
  genotype <- as.character(genotype)
  count <- nchar(gsub("[^4]", "", genotype))
  count[is.na(genotype) | genotype == ""] <- NA_integer_
  as.integer(count)
}

adas <- load_rda_object(paths$adas, "ADAS")
mri <- load_rda_object(paths$mri, "UCSFFSX7")
adsl <- load_rda_object(paths$adsl, "ADSL")
ptdemog <- load_rda_object(paths$ptdemog, "PTDEMOG")
apoeres <- load_rda_object(paths$apoeres, "APOERES")

adas$VISDATE <- as_adni_date(adas$VISDATE)
adas$TOTAL13 <- suppressWarnings(as.numeric(as.character(adas$TOTAL13)))
adas_valid <- stats::complete.cases(adas[, c("RID", "VISDATE", "TOTAL13")]) &
  is.finite(adas$TOTAL13) & adas$TOTAL13 >= 0 & adas$TOTAL13 <= 85 &
  (is.na(adas$HAS_QC_ERROR) | !qc_error(adas$HAS_QC_ERROR))
adas <- adas[adas_valid, , drop = FALSE]
adas <- adas[order(adas$RID, adas$VISDATE, adas$ID, na.last = TRUE), ,
             drop = FALSE]
adas <- adas[!duplicated(adas[, c("RID", "VISDATE")]), , drop = FALSE]

mri$EXAMDATE <- as_adni_date(mri$EXAMDATE)
mri_numeric <- c(
  "ST10CV", "ST29SV", "ST88SV", "ST24CV", "ST83CV",
  "ST37SV", "ST96SV", "ST30SV", "ST89SV", "ST127SV", "ST9SV"
)
for (column in mri_numeric) {
  mri[[column]] <- suppressWarnings(as.numeric(as.character(mri[[column]])))
}
mri_valid <- stats::complete.cases(
  mri[, c("RID", "EXAMDATE", mri_numeric), drop = FALSE]
) & apply(mri[, mri_numeric, drop = FALSE], 1L, function(value) {
  all(is.finite(value) & value > 0)
}) & mri$ST10CV <= 3000000 &
  (is.na(mri$OVERALLQC) | mri$OVERALLQC != "Fail") &
  (is.na(mri$VENTQC) | mri$VENTQC != "Fail") &
  (is.na(mri$HIPPOQC) | mri$HIPPOQC != "Fail")
mri <- mri[mri_valid, , drop = FALSE]

mri$.quality_rank <-
  4L * (mri$STATUS == "complete") +
  2L * (mri$OVERALLQC == "Pass") +
  1L * (mri$VENTQC == "Pass" & mri$HIPPOQC == "Pass")
mri$.quality_rank[is.na(mri$.quality_rank)] <- 0L
mri$.field_strength_rank <- suppressWarnings(
  as.numeric(gsub("[^0-9.]", "", mri$FIELD_STRENGTH))
)
mri <- mri[order(
  mri$RID, mri$EXAMDATE, -mri$.quality_rank, -mri$.field_strength_rank,
  mri$IMAGEUID, na.last = TRUE
), , drop = FALSE]
mri <- mri[!duplicated(mri[, c("RID", "EXAMDATE")]), , drop = FALSE]

mri_rows <- split(seq_len(nrow(mri)), as.character(mri$RID))
selected_mri <- rep(NA_integer_, nrow(adas))
mri_gap <- rep(NA_real_, nrow(adas))
for (row in seq_len(nrow(adas))) {
  candidates <- mri_rows[[as.character(adas$RID[[row]])]]
  if (is.null(candidates)) next
  gaps <- as.numeric(mri$EXAMDATE[candidates] - adas$VISDATE[[row]])
  eligible <- which(is.finite(gaps) & abs(gaps) <= window_days)
  if (length(eligible) == 0L) next
  selected_local <- eligible[order(abs(gaps[eligible]), gaps[eligible])][[1L]]
  selected_mri[[row]] <- candidates[[selected_local]]
  mri_gap[[row]] <- gaps[[selected_local]]
}

matched <- which(!is.na(selected_mri))
adas <- adas[matched, , drop = FALSE]
mri_gap <- mri_gap[matched]
mri <- mri[selected_mri[matched], , drop = FALSE]

analysis <- data.frame(
  RID = as.integer(adas$RID),
  EXAMDATE = adas$VISDATE,
  VISCODE = as.character(adas$VISCODE),
  VISCODE2 = as.character(adas$VISCODE2),
  ADAS13 = as.numeric(adas$TOTAL13),
  COLPROT = as.character(adas$COLPROT),
  SITE = as.character(adas$SITEID),
  MRI_DATE = mri$EXAMDATE,
  MRI_gap_days = mri_gap,
  MRI_FIELD_STRENGTH = as.character(mri$FIELD_STRENGTH),
  ICV = mri$ST10CV,
  Hippocampus = mri$ST29SV + mri$ST88SV,
  Entorhinal = mri$ST24CV + mri$ST83CV,
  Ventricles = rowSums(mri[, c(
    "ST37SV", "ST96SV", "ST30SV", "ST89SV", "ST127SV", "ST9SV"
  )], na.rm = FALSE),
  stringsAsFactors = FALSE
)

adsl$RID <- suppressWarnings(as.integer(as.character(adsl$SUBJID)))
adsl$ENRLDT <- as_adni_date(adsl$ENRLDT)
subject <- adsl[, c("RID", "AGE", "ENRLDT", "SEX", "EDUC", "APOE")]
names(subject)[names(subject) == "SEX"] <- "PTGENDER"
names(subject)[names(subject) == "EDUC"] <- "PTEDUCAT"
subject$APOE4 <- count_apoe4(subject$APOE)

ptdemog$VISDATE <- as_adni_date(ptdemog$VISDATE)
demog <- first_nonmissing_by_rid(
  ptdemog, c("PTGENDER", "PTEDUCAT"), date_column = "VISDATE"
)
subject <- merge(subject, demog, by = "RID", all.x = TRUE,
                 suffixes = c("", ".DEMOG"))
subject$PTGENDER <- ifelse(
  is.na(subject$PTGENDER) | subject$PTGENDER == "",
  subject$PTGENDER.DEMOG, subject$PTGENDER
)
subject$PTEDUCAT <- ifelse(
  is.na(subject$PTEDUCAT),
  suppressWarnings(as.numeric(subject$PTEDUCAT.DEMOG)), subject$PTEDUCAT
)

apoe <- first_nonmissing_by_rid(apoeres, "GENOTYPE")
apoe$APOE4.RES <- count_apoe4(apoe$GENOTYPE)
subject <- merge(subject, apoe[, c("RID", "APOE4.RES")], by = "RID",
                 all.x = TRUE)
subject$APOE4 <- ifelse(is.na(subject$APOE4), subject$APOE4.RES, subject$APOE4)

analysis$.row_order <- seq_len(nrow(analysis))
analysis <- merge(analysis, subject, by = "RID", all.x = TRUE, sort = FALSE)
analysis <- analysis[order(analysis$.row_order), , drop = FALSE]
analysis$.row_order <- NULL
analysis$AGE_AT_VISIT <- analysis$AGE +
  as.numeric(analysis$EXAMDATE - analysis$ENRLDT) / 365.25

required <- c(
  "RID", "EXAMDATE", "AGE", "AGE_AT_VISIT", "PTGENDER", "PTEDUCAT",
  "APOE4", "COLPROT", "SITE", "MRI_FIELD_STRENGTH", "MRI_gap_days",
  "ICV", "Hippocampus", "Entorhinal", "Ventricles", "ADAS13"
)
complete <- stats::complete.cases(analysis[, required, drop = FALSE])
complete_data <- analysis[complete, , drop = FALSE]
visit_counts <- table(as.character(complete_data$RID))
eligible_ids <- names(visit_counts[visit_counts >= min_visits])
repeated_data <- complete_data[
  as.character(complete_data$RID) %in% eligible_ids, , drop = FALSE
]

if (anyDuplicated(repeated_data[, c("RID", "EXAMDATE")])) {
  stop("Duplicate RID/EXAMDATE rows remain after matching.", call. = FALSE)
}

dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/adni_mri_adas13_exploratory", recursive = TRUE,
           showWarnings = FALSE)
full_path <- sprintf("data/ADNI_MRI_ADAS13_%dd_complete.csv", window_days)
repeated_path <- sprintf(
  "data/ADNI_MRI_ADAS13_%dd_min%d.csv", window_days, min_visits
)
utils::write.csv(complete_data, full_path, row.names = FALSE, na = "")
utils::write.csv(repeated_data, repeated_path, row.names = FALSE, na = "")

audit <- data.frame(
  date_window_days = window_days,
  date_matched_adas_rows = nrow(adas),
  complete_matched_rows = nrow(complete_data),
  complete_subjects = length(visit_counts),
  min_visits = min_visits,
  repeated_rows = nrow(repeated_data),
  repeated_subjects = length(eligible_ids),
  median_abs_mri_gap_days = stats::median(abs(complete_data$MRI_gap_days)),
  p90_abs_mri_gap_days = unname(stats::quantile(
    abs(complete_data$MRI_gap_days), 0.9
  )),
  stringsAsFactors = FALSE
)
utils::write.csv(
  audit,
  "results/adni_mri_adas13_exploratory/preparation_audit.csv",
  row.names = FALSE
)

message("Prepared MRI--ADAS13 tables:")
message("  ", full_path)
message("  ", repeated_path)
message("Repeated cohort: ", nrow(repeated_data), " rows from ",
        length(eligible_ids), " subjects")
