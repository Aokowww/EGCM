#!/usr/bin/env Rscript

# Build the visit-level ADNI analysis table used by Designs A and B.
#
# The PET examination is the anchor visit. MRI and cognitive assessments are
# matched within a pre-specified symmetric date window without inspecting any
# outcome-test result. Restricted ADNI data remain under data/ and are ignored
# by Git.

args <- commandArgs(trailingOnly = TRUE)
snapshot_dir <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "data/adni_raw/2026-07-26"
}
primary_window <- if (length(args) >= 2L) as.integer(args[[2L]]) else 90L
sensitivity_window <- if (length(args) >= 3L) as.integer(args[[3L]]) else 180L

if (!is.finite(primary_window) || primary_window <= 0L ||
    !is.finite(sensitivity_window) ||
    sensitivity_window < primary_window) {
  stop("Require 0 < primary_window <= sensitivity_window.", call. = FALSE)
}

paths <- list(
  pet = file.path(snapshot_dir, "UCBERKELEYAV45_04_26_22_26Jul2026.csv"),
  mri = file.path(snapshot_dir, "UCSFFSX7_26Jul2026.csv"),
  cdr = file.path(snapshot_dir, "CDR_26Jul2026.csv"),
  mmse = file.path(snapshot_dir, "MMSE_26Jul2026.csv"),
  moca = file.path(snapshot_dir, "MOCA_26Jul2026.csv"),
  adsl = file.path(snapshot_dir, "extracted", "ADNIMERGE2", "data", "ADSL.rda"),
  registry = file.path(
    snapshot_dir, "extracted", "ADNIMERGE2", "data", "REGISTRY.rda"
  ),
  ptdemog = file.path(
    snapshot_dir, "extracted", "ADNIMERGE2", "data", "PTDEMOG.rda"
  ),
  apoeres = file.path(
    snapshot_dir, "extracted", "ADNIMERGE2", "data", "APOERES.rda"
  )
)
missing_files <- names(paths)[!file.exists(unlist(paths, use.names = FALSE))]
if (length(missing_files) > 0L) {
  stop(
    "Missing ADNI input(s): ", paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

read_adni_csv <- function(path) {
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
}

load_rda_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!expected_name %in% loaded) {
    stop(
      "Expected object ", expected_name, " was not found in ", path,
      call. = FALSE
    )
  }
  environment[[expected_name]]
}

as_adni_date <- function(x) {
  value <- as.Date(as.character(x))
  value
}

deduplicate_by_date <- function(data, value_columns, date_column) {
  data <- data[stats::complete.cases(
    data[, c("RID", date_column, value_columns), drop = FALSE]
  ), , drop = FALSE]
  data <- data[order(
    data$RID,
    data[[date_column]],
    if ("HAS_QC_ERROR" %in% names(data)) {
      data$HAS_QC_ERROR
    } else {
      rep(FALSE, nrow(data))
    },
    if ("ID" %in% names(data)) data$ID else seq_len(nrow(data)),
    na.last = TRUE
  ), , drop = FALSE]
  data[!duplicated(data[, c("RID", date_column), drop = FALSE]), , drop = FALSE]
}

nearest_match <- function(anchor, candidate, candidate_date, value_columns,
                          prefix, window_days) {
  split_rows <- split(seq_len(nrow(candidate)), as.character(candidate$RID))
  result <- as.data.frame(
    matrix(NA, nrow = nrow(anchor), ncol = length(value_columns)),
    stringsAsFactors = FALSE
  )
  names(result) <- paste0(prefix, value_columns)
  result[[paste0(prefix, "date")]] <- as.Date(NA)
  result[[paste0(prefix, "gap_days")]] <- NA_real_

  for (row in seq_len(nrow(anchor))) {
    rows <- split_rows[[as.character(anchor$RID[[row]])]]
    if (is.null(rows)) next
    gaps <- as.numeric(candidate[[candidate_date]][rows] - anchor$EXAMDATE[[row]])
    eligible <- which(is.finite(gaps) & abs(gaps) <= window_days)
    if (length(eligible) == 0L) next

    # Deterministic tie-break: minimum absolute gap, then the earlier measure.
    selected_local <- eligible[order(abs(gaps[eligible]), gaps[eligible])][[1L]]
    selected <- rows[[selected_local]]
    for (column in value_columns) {
      result[[paste0(prefix, column)]][[row]] <- candidate[[column]][[selected]]
    }
    result[[paste0(prefix, "date")]][[row]] <-
      candidate[[candidate_date]][[selected]]
    result[[paste0(prefix, "gap_days")]][[row]] <- gaps[[selected_local]]
  }
  result
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
    output[[column]] <- vapply(
      rows,
      function(index) {
        values <- data[[column]][index]
        values <- values[!is.na(values) & as.character(values) != ""]
        if (length(values) == 0L) NA_character_ else as.character(values[[1L]])
      },
      character(1)
    )
  }
  output
}

count_apoe4 <- function(genotype) {
  genotype <- as.character(genotype)
  count <- nchar(gsub("[^4]", "", genotype))
  count[is.na(genotype) | genotype == ""] <- NA_integer_
  as.integer(count)
}

pet <- read_adni_csv(paths$pet)
mri <- read_adni_csv(paths$mri)
cdr <- read_adni_csv(paths$cdr)
mmse <- read_adni_csv(paths$mmse)
moca <- read_adni_csv(paths$moca)
adsl <- load_rda_object(paths$adsl, "ADSL")
registry <- load_rda_object(paths$registry, "REGISTRY")
ptdemog <- load_rda_object(paths$ptdemog, "PTDEMOG")
apoeres <- load_rda_object(paths$apoeres, "APOERES")

pet$EXAMDATE <- as_adni_date(pet$EXAMDATE)
mri$EXAMDATE <- as_adni_date(mri$EXAMDATE)
cdr$VISDATE <- as_adni_date(cdr$VISDATE)
mmse$VISDATE <- as_adni_date(mmse$VISDATE)
moca$VISDATE <- as_adni_date(moca$VISDATE)
registry$EXAMDATE <- as_adni_date(registry$EXAMDATE)
ptdemog$VISDATE <- as_adni_date(ptdemog$VISDATE)

pet_outcomes <- c(
  "SUMMARYSUVR_COMPOSITE_REFNORM",
  "SUMMARYSUVR_WHOLECEREBNORM"
)
pet <- deduplicate_by_date(pet, pet_outcomes, "EXAMDATE")
if (anyDuplicated(pet[, c("RID", "EXAMDATE"), drop = FALSE])) {
  stop("PET anchor has duplicate RID/EXAMDATE rows.", call. = FALSE)
}

mri_columns <- c(
  "ST10CV", "ST29SV", "ST88SV", "ST24CV", "ST83CV",
  "ST37SV", "ST96SV", "ST30SV", "ST89SV", "ST127SV", "ST9SV",
  "FIELD_STRENGTH", "STATUS", "OVERALLQC", "VENTQC", "HIPPOQC"
)
numeric_mri <- mri_columns[grepl("^ST[0-9]", mri_columns)]
mri_valid <- stats::complete.cases(
  mri[, c("RID", "EXAMDATE", numeric_mri), drop = FALSE]
)
mri_valid <- mri_valid & apply(
  mri[, numeric_mri, drop = FALSE],
  1L,
  function(value) all(is.finite(value) & value > 0)
)
# Preserve the pre-specified broad physiological/unit check used by the
# analysis configuration. Values above this bound were confined to un-QC'd
# partial rows in the downloaded snapshot.
mri_valid <- mri_valid & mri$ST10CV <= 3000000
mri_valid <- mri_valid &
  (is.na(mri$OVERALLQC) | mri$OVERALLQC != "Fail") &
  (is.na(mri$VENTQC) | mri$VENTQC != "Fail") &
  (is.na(mri$HIPPOQC) | mri$HIPPOQC != "Fail")
mri <- mri[mri_valid, , drop = FALSE]

# Prefer the most complete/QC-supported scan when more than one scan is
# available on the same date. No PET or cognitive value enters this choice.
mri$.quality_rank <- (
  4L * (mri$STATUS == "complete") +
  2L * (mri$OVERALLQC == "Pass") +
  1L * (mri$VENTQC == "Pass" & mri$HIPPOQC == "Pass")
)
mri$.quality_rank[is.na(mri$.quality_rank)] <- 0L
mri$.field_strength_rank <- suppressWarnings(
  as.numeric(gsub("[^0-9.]", "", mri$FIELD_STRENGTH))
)
mri <- mri[order(
  mri$RID, mri$EXAMDATE, -mri$.quality_rank, -mri$.field_strength_rank,
  mri$IMAGEUID, na.last = TRUE
), , drop = FALSE]
mri <- mri[!duplicated(mri[, c("RID", "EXAMDATE"), drop = FALSE]), ,
           drop = FALSE]
message(
  "Eligible MRI rows after value/QC checks and same-day deduplication: ",
  nrow(mri)
)

cdr <- deduplicate_by_date(cdr, "CDRSB", "VISDATE")
mmse <- deduplicate_by_date(mmse, "MMSCORE", "VISDATE")
moca <- deduplicate_by_date(moca, "MOCA", "VISDATE")
registry <- deduplicate_by_date(registry, c("COLPROT", "SITEID"), "EXAMDATE")

build_aligned <- function(window_days) {
  aligned <- data.frame(
    RID = as.integer(pet$RID),
    EXAMDATE = pet$EXAMDATE,
    VISCODE = pet$VISCODE,
    VISCODE2 = pet$VISCODE2,
    AV45 = pet$SUMMARYSUVR_COMPOSITE_REFNORM,
    AV45_WHOLE_CEREBELLUM = pet$SUMMARYSUVR_WHOLECEREBNORM,
    stringsAsFactors = FALSE
  )

  aligned <- cbind(
    aligned,
    nearest_match(
      aligned, mri, "EXAMDATE", mri_columns, "MRI_", window_days
    ),
    nearest_match(aligned, cdr, "VISDATE", "CDRSB", "CDR_", window_days),
    nearest_match(aligned, mmse, "VISDATE", "MMSCORE", "MMSE_", window_days),
    nearest_match(aligned, moca, "VISDATE", "MOCA", "MOCA_", window_days),
    nearest_match(
      aligned, registry, "EXAMDATE", c("COLPROT", "SITEID"),
      "REG_", window_days
    )
  )

  aligned$ICV <- aligned$MRI_ST10CV
  aligned$Hippocampus <- aligned$MRI_ST29SV + aligned$MRI_ST88SV
  aligned$Entorhinal <- aligned$MRI_ST24CV + aligned$MRI_ST83CV
  aligned$Ventricles <- rowSums(aligned[, c(
    "MRI_ST37SV", "MRI_ST96SV", "MRI_ST30SV", "MRI_ST89SV",
    "MRI_ST127SV", "MRI_ST9SV"
  )], na.rm = FALSE)
  aligned$CDRSB <- aligned$CDR_CDRSB
  aligned$MMSE <- aligned$MMSE_MMSCORE
  aligned$MOCA <- aligned$MOCA_MOCA
  aligned$COLPROT <- aligned$REG_COLPROT
  aligned$SITE <- aligned$REG_SITEID
  aligned
}

adsl$RID <- suppressWarnings(as.integer(as.character(adsl$SUBJID)))
adsl$ENRLDT <- as_adni_date(adsl$ENRLDT)
subject <- adsl[, c(
  "RID", "AGE", "ENRLDT", "SEX", "EDUC", "APOE", "SITEID", "ORIGPROT"
)]
names(subject)[names(subject) == "SEX"] <- "PTGENDER"
names(subject)[names(subject) == "EDUC"] <- "PTEDUCAT"
subject$APOE4 <- count_apoe4(subject$APOE)

demog <- first_nonmissing_by_rid(
  ptdemog,
  c("PTGENDER", "PTEDUCAT"),
  date_column = "VISDATE"
)
subject <- merge(subject, demog, by = "RID", all.x = TRUE, suffixes = c("", ".DEMOG"))
subject$PTGENDER <- ifelse(
  is.na(subject$PTGENDER) | subject$PTGENDER == "",
  subject$PTGENDER.DEMOG,
  subject$PTGENDER
)
subject$PTEDUCAT <- ifelse(
  is.na(subject$PTEDUCAT),
  suppressWarnings(as.numeric(subject$PTEDUCAT.DEMOG)),
  subject$PTEDUCAT
)

apoe <- first_nonmissing_by_rid(apoeres, "GENOTYPE")
apoe$APOE4.RES <- count_apoe4(apoe$GENOTYPE)
subject <- merge(
  subject,
  apoe[, c("RID", "APOE4.RES")],
  by = "RID",
  all.x = TRUE
)
subject$APOE4 <- ifelse(is.na(subject$APOE4), subject$APOE4.RES, subject$APOE4)

attach_subject <- function(aligned) {
  output <- merge(aligned, subject, by = "RID", all.x = TRUE, sort = FALSE)
  output <- output[match(
    paste(pet$RID, pet$EXAMDATE),
    paste(output$RID, output$EXAMDATE)
  ), , drop = FALSE]
  output$AGE_AT_VISIT <- output$AGE +
    as.numeric(output$EXAMDATE - output$ENRLDT) / 365.25
  output$PTGENDER <- factor(output$PTGENDER)
  output$COLPROT <- ifelse(
    is.na(output$COLPROT) | output$COLPROT == "",
    as.character(output$ORIGPROT),
    as.character(output$COLPROT)
  )
  output$SITE <- ifelse(is.na(output$SITE), output$SITEID, output$SITE)
  output
}

primary <- attach_subject(build_aligned(primary_window))
sensitivity <- attach_subject(build_aligned(sensitivity_window))

dir.create("data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/adni_preparation", recursive = TRUE, showWarnings = FALSE)
primary_path <- sprintf("data/ADNI_ANALYTIC_A_B_%dd.csv", primary_window)
sensitivity_path <- sprintf(
  "data/ADNI_ANALYTIC_A_B_%dd.csv", sensitivity_window
)
utils::write.csv(primary, primary_path, row.names = FALSE, na = "")
utils::write.csv(sensitivity, sensitivity_path, row.names = FALSE, na = "")

analysis_required <- c(
  "RID", "EXAMDATE", "AGE", "AGE_AT_VISIT", "PTGENDER", "PTEDUCAT",
  "APOE4", "COLPROT", "SITE", "ICV", "MOCA", "MMSE", "CDRSB",
  "Ventricles", "Hippocampus", "Entorhinal", "AV45"
)
cohort_row <- function(data, label, window_days) {
  complete <- stats::complete.cases(data[, analysis_required, drop = FALSE])
  counts <- table(data$RID[complete])
  data.frame(
    cohort = label,
    max_abs_gap_days = window_days,
    pet_anchor_rows = nrow(data),
    complete_rows = sum(complete),
    complete_subjects = length(counts),
    subjects_with_at_least_3_visits = sum(counts >= 3L),
    median_complete_visits = if (length(counts)) {
      stats::median(as.numeric(counts))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}
cohort_flow <- rbind(
  cohort_row(primary, "primary", primary_window),
  cohort_row(sensitivity, "date_window_sensitivity", sensitivity_window)
)
utils::write.csv(
  cohort_flow,
  "results/adni_preparation/cohort_flow.csv",
  row.names = FALSE
)

alignment_columns <- c(
  "MRI_gap_days", "CDR_gap_days", "MMSE_gap_days", "MOCA_gap_days",
  "REG_gap_days"
)
alignment_audit <- do.call(rbind, lapply(
  alignment_columns,
  function(column) {
    gap <- primary[[column]]
    observed <- abs(gap[is.finite(gap)])
    data.frame(
      source = sub("_gap_days$", "", column),
      matched_rows = length(observed),
      median_abs_gap_days = if (length(observed)) stats::median(observed) else NA,
      p90_abs_gap_days = if (length(observed)) {
        unname(stats::quantile(observed, 0.9))
      } else {
        NA
      },
      max_abs_gap_days = if (length(observed)) max(observed) else NA,
      stringsAsFactors = FALSE
    )
  }
))
utils::write.csv(
  alignment_audit,
  "results/adni_preparation/alignment_audit.csv",
  row.names = FALSE
)

variable_dictionary <- data.frame(
  analysis_name = c(
    "AV45", "AV45_WHOLE_CEREBELLUM", "Ventricles", "Hippocampus",
    "Entorhinal", "ICV", "CDRSB", "MMSE", "MOCA"
  ),
  source = c(
    "UCBERKELEYAV45", "UCBERKELEYAV45", rep("UCSFFSX7", 4L),
    "CDR", "MMSE", "MOCA"
  ),
  source_fields = c(
    "SUMMARYSUVR_COMPOSITE_REFNORM",
    "SUMMARYSUVR_WHOLECEREBNORM",
    "ST37SV + ST96SV + ST30SV + ST89SV + ST127SV + ST9SV",
    "ST29SV + ST88SV",
    "ST24CV + ST83CV",
    "ST10CV",
    "CDRSB", "MMSCORE", "MOCA"
  ),
  unit = c("ratio", "ratio", rep("mm3", 4L), "points", "points", "points"),
  role = c(
    "primary Y", "PET-reference sensitivity Y", rep("X", 4L),
    rep("conditioning variable", 3L)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  variable_dictionary,
  "results/adni_preparation/variable_dictionary.csv",
  row.names = FALSE
)

message("Prepared ADNI tables:")
message("  ", primary_path)
message("  ", sensitivity_path)
message("Cohort audit: results/adni_preparation/cohort_flow.csv")
