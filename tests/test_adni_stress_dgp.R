source("code/simulation/adni_stress_dgp.R")

simulation <- list(
  seed = 17L,
  subjects = 500L,
  scenarios = list(
    fixed = list(
      visit_type = "fixed",
      visits_fixed = 6L,
      random_slope_sd = 0.75
    ),
    serial = list(
      visit_type = "fixed",
      visits_fixed = 6L,
      ar1_rho = 0.6
    ),
    unequal = list(
      visit_type = "unequal",
      visits_min = 2L,
      visits_max = 8L
    ),
    informative = list(
      visit_type = "informative_x",
      visits_min = 2L,
      visits_max = 8L,
      visits_centre = 5,
      informative_strength = 2
    )
  )
)

fixed <- generate_adni_stress_data(
  "fixed", simulation$scenarios$fixed, simulation, 1L
)
stopifnot(
  length(unique(fixed$RID)) == simulation$subjects,
  all(table(fixed$RID) == 6L),
  stats::sd(unique(fixed[c("RID", "truth_u_x1")])$truth_u_x1) > 0.5,
  !anyNA(fixed)
)

serial <- generate_adni_stress_data(
  "serial", simulation$scenarios$serial, simulation, 1L
)
serial_pairs <- do.call(
  rbind,
  lapply(split(serial$truth_error_y, serial$RID), function(x) {
    cbind(x[-length(x)], x[-1L])
  })
)
stopifnot(stats::cor(serial_pairs[, 1L], serial_pairs[, 2L]) > 0.45)

unequal <- generate_adni_stress_data(
  "unequal", simulation$scenarios$unequal, simulation, 1L
)
unequal_counts <- table(unequal$RID)
stopifnot(
  min(unequal_counts) >= 2L,
  max(unequal_counts) <= 8L,
  length(unique(unequal_counts)) > 1L
)

informative <- generate_adni_stress_data(
  "informative", simulation$scenarios$informative, simulation, 1L
)
subject_truth <- informative[
  !duplicated(informative$RID),
  c("RID", "truth_u_x0", "truth_u_y0")
]
subject_truth$visits <- as.numeric(table(informative$RID)[
  as.character(subject_truth$RID)
])
stopifnot(
  stats::cor(subject_truth$truth_u_x0, subject_truth$visits) > 0.5,
  abs(stats::cor(subject_truth$truth_u_x0, subject_truth$truth_u_y0)) < 0.15
)

message("adni_stress_dgp.R tests passed")
