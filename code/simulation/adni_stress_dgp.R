# Synthetic null generators for the formal Design A stress calibration.
#
# Every scenario preserves the population-level null by keeping all X-side
# latent variables and errors independent of all Y-side latent variables and
# errors conditional on the observed Z variables. Dependence is allowed within
# a subject and within either side.

stress_value <- function(scenario, name, default) {
  value <- scenario[[name]]
  if (is.null(value)) default else value
}

simulate_subject_ar1 <- function(visits, rho = 0, sd = 1) {
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("AR(1) rho must be finite and strictly between -1 and 1.",
         call. = FALSE)
  }
  values <- numeric(sum(visits))
  start <- 1L
  innovation_sd <- sd * sqrt(1 - rho^2)
  for (m in visits) {
    rows <- start:(start + m - 1L)
    series <- numeric(m)
    series[[1L]] <- stats::rnorm(1L, sd = sd)
    if (m > 1L) {
      for (j in 2:m) {
        series[[j]] <- rho * series[[j - 1L]] +
          stats::rnorm(1L, sd = innovation_sd)
      }
    }
    values[rows] <- series
    start <- start + m
  }
  values
}

stress_visit_counts <- function(scenario, n_subjects, x_intercept) {
  visit_type <- as.character(stress_value(scenario, "visit_type", "fixed"))
  minimum <- as.integer(stress_value(scenario, "visits_min", 2L))
  maximum <- as.integer(stress_value(scenario, "visits_max", 8L))
  if (minimum < 2L || maximum < minimum) {
    stop("Stress visits_min/visits_max are invalid.", call. = FALSE)
  }

  if (visit_type == "fixed") {
    visits <- rep(as.integer(stress_value(scenario, "visits_fixed", 6L)),
                  n_subjects)
  } else if (visit_type == "unequal") {
    visits <- sample(seq.int(minimum, maximum), n_subjects, replace = TRUE)
  } else if (visit_type == "informative_x") {
    centre <- as.numeric(stress_value(scenario, "visits_centre", 5))
    strength <- as.numeric(stress_value(
      scenario, "informative_strength", 2
    ))
    visits <- round(
      centre + strength * tanh(x_intercept) + stats::rnorm(n_subjects, sd = 1)
    )
    visits <- pmin(maximum, pmax(minimum, visits))
  } else {
    stop("Unsupported stress visit_type: ", visit_type, call. = FALSE)
  }
  as.integer(visits)
}

generate_adni_stress_data <- function(scenario_name, scenario, simulation,
                                      repetition) {
  scenario_names <- names(simulation$scenarios)
  scenario_index <- match(scenario_name, scenario_names)
  if (is.na(scenario_index)) {
    stop("Unknown stress scenario: ", scenario_name, call. = FALSE)
  }
  set.seed(
    as.integer(simulation$seed) + as.integer(repetition) +
      10000L * scenario_index
  )

  n_subjects <- as.integer(simulation$subjects)
  random_intercept_sd <- as.numeric(stress_value(
    scenario, "random_intercept_sd", 1
  ))
  random_slope_sd <- as.numeric(stress_value(
    scenario, "random_slope_sd", 0
  ))
  ar1_rho <- as.numeric(stress_value(scenario, "ar1_rho", 0))
  interaction_x <- as.numeric(stress_value(
    scenario, "omitted_x_interaction", 0
  ))

  u_x0 <- stats::rnorm(n_subjects, sd = random_intercept_sd)
  u_y0 <- stats::rnorm(n_subjects, sd = random_intercept_sd)
  u_x1 <- stats::rnorm(n_subjects, sd = random_slope_sd)
  u_y1 <- stats::rnorm(n_subjects, sd = random_slope_sd)
  visits <- stress_visit_counts(scenario, n_subjects, u_x0)

  id <- rep(seq_len(n_subjects), visits)
  visit <- sequence(visits)
  time <- (visit - 1) / pmax(visits[id] - 1, 1)
  z_subject <- stats::rnorm(n_subjects)
  z <- 0.6 * z_subject[id] + stats::rnorm(length(id), sd = 0.8)

  error_x1 <- simulate_subject_ar1(visits, rho = ar1_rho)
  error_x2 <- simulate_subject_ar1(visits, rho = ar1_rho)
  error_x3 <- simulate_subject_ar1(visits, rho = ar1_rho)
  error_y <- simulate_subject_ar1(visits, rho = ar1_rho)
  x_interaction <- interaction_x * z * time

  data.frame(
    RID = id,
    visit = visit,
    time_years = time,
    z = z,
    X1 = sin(z) + 0.3 * time + x_interaction +
      u_x0[id] + u_x1[id] * time + error_x1,
    X2 = z^2 - 1 + 0.2 * time + x_interaction +
      0.8 * u_x0[id] + 0.8 * u_x1[id] * time + error_x2,
    X3 = tanh(z) - 0.2 * time + x_interaction +
      0.6 * u_x0[id] + 0.6 * u_x1[id] * time + error_x3,
    Y = cos(z) + 0.4 * time + u_y0[id] + u_y1[id] * time + error_y,
    truth_u_x0 = u_x0[id],
    truth_u_y0 = u_y0[id],
    truth_u_x1 = u_x1[id],
    truth_u_y1 = u_y1[id],
    truth_error_y = error_y
  )
}
