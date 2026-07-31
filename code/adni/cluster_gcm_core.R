# Cluster-robust GCM utilities for repeated-measures data.
#
# Design A targets population-level conditional independence and is the primary
# analysis. Design B conditions on an estimated subject effect and is explicitly
# a model-based sensitivity analysis for sparse longitudinal data.

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package '", package, "' is required.", call. = FALSE)
  }
}

assert_columns <- function(data, columns) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

make_subject_folds <- function(id, k = 5L, seed = 20260726L) {
  ids <- unique(as.character(id))
  if (length(ids) < 2L) {
    stop("At least two subjects are required.", call. = FALSE)
  }
  k <- min(as.integer(k), length(ids))
  counts <- table(as.character(id))

  set.seed(seed)
  ordered_ids <- names(sort(counts + stats::runif(length(counts), 0, 0.01),
                            decreasing = TRUE))
  fold_load <- numeric(k)
  fold_n <- integer(k)
  assignment <- integer(length(ordered_ids))
  names(assignment) <- ordered_ids

  for (subject in ordered_ids) {
    candidates <- which(fold_load == min(fold_load))
    if (length(candidates) > 1L) {
      candidates <- candidates[fold_n[candidates] == min(fold_n[candidates])]
    }
    chosen <- sample(candidates, 1L)
    assignment[subject] <- chosen
    fold_load[chosen] <- fold_load[chosen] + counts[[subject]]
    fold_n[chosen] <- fold_n[chosen] + 1L
  }

  unname(assignment[as.character(id)])
}

make_within_subject_folds <- function(id, k = 2L, seed = 20260726L) {
  k <- as.integer(k)
  if (k < 2L) {
    stop("Within-subject cross-fitting requires at least two folds.", call. = FALSE)
  }

  set.seed(seed)
  folds <- integer(length(id))
  subject_rows <- split(seq_along(id), as.character(id))
  for (rows in subject_rows) {
    if (length(rows) < k) {
      stop("Every Design B subject must have at least k visits.", call. = FALSE)
    }
    labels <- rep(seq_len(k), length.out = length(rows))
    folds[rows] <- sample(labels, length(labels), replace = FALSE)
  }
  folds
}

make_gam_formula <- function(outcome, rhs) {
  stats::as.formula(
    paste0("`", outcome, "` ~ ", rhs, " + s(.subject_re, bs = 're')")
  )
}

make_equal_subject_weights <- function(id) {
  counts <- table(as.character(id))
  as.numeric(1 / counts[as.character(id)])
}

fit_predict_gam <- function(train, test, outcome, rhs,
                            prediction = c("population", "subject"),
                            case_weights = NULL) {
  require_namespace("mgcv")
  prediction <- match.arg(prediction)
  if (is.null(case_weights)) {
    case_weights <- rep(1, nrow(train))
  }
  if (length(case_weights) != nrow(train) ||
      any(!is.finite(case_weights) | case_weights <= 0)) {
    stop("case_weights must be positive, finite, and match the training rows.",
         call. = FALSE)
  }
  train$.case_weight <- case_weights
  formula <- make_gam_formula(outcome, rhs)
  fit <- mgcv::gam(
    formula,
    data = train,
    weights = .case_weight,
    method = "REML",
    na.action = stats::na.fail,
    drop.unused.levels = TRUE
  )

  newdata <- test
  if (prediction == "population") {
    # Held-out subject levels are intentionally replaced before prediction.
    # The random-effect term is excluded, so this replacement has no effect on
    # the population prediction and avoids new-factor-level failures.
    reference_subject <- as.character(train$.subject_re[[1L]])
    newdata$.subject_re <- factor(
      rep(reference_subject, nrow(newdata)),
      levels = levels(train$.subject_re)
    )
    fitted <- stats::predict(
      fit,
      newdata = newdata,
      type = "response",
      exclude = "s(.subject_re)"
    )
  } else {
    fitted <- stats::predict(fit, newdata = newdata, type = "response")
  }

  list(prediction = as.numeric(fitted), model = fit)
}

crossfit_residuals <- function(data, outcome, rhs, folds,
                               prediction = c("population", "subject"),
                               case_weights = NULL) {
  prediction <- match.arg(prediction)
  if (is.null(case_weights)) {
    case_weights <- rep(1, nrow(data))
  }
  if (length(case_weights) != nrow(data)) {
    stop("case_weights must match the analysis rows.", call. = FALSE)
  }
  residual <- rep(NA_real_, nrow(data))
  prediction_value <- rep(NA_real_, nrow(data))
  models <- vector("list", max(folds))

  for (fold in seq_len(max(folds))) {
    test_rows <- which(folds == fold)
    train_rows <- which(folds != fold)
    fitted <- fit_predict_gam(
      train = data[train_rows, , drop = FALSE],
      test = data[test_rows, , drop = FALSE],
      outcome = outcome,
      rhs = rhs,
      prediction = prediction,
      case_weights = case_weights[train_rows]
    )
    prediction_value[test_rows] <- fitted$prediction
    residual[test_rows] <- data[[outcome]][test_rows] - fitted$prediction
    models[[fold]] <- fitted$model
  }

  list(residual = residual, prediction = prediction_value, models = models)
}

subject_score_matrix <- function(id, residual_x, residual_y) {
  residual_x <- as.matrix(residual_x)
  if (nrow(residual_x) != length(residual_y)) {
    stop("X and Y residuals have incompatible lengths.", call. = FALSE)
  }
  products <- residual_x * residual_y
  subjects <- unique(as.character(id))
  scores <- vapply(
    seq_len(ncol(products)),
    function(j) {
      as.numeric(tapply(products[, j], as.character(id), mean)[subjects])
    },
    numeric(length(subjects))
  )
  if (is.null(dim(scores))) {
    scores <- matrix(scores, ncol = 1L)
  }
  rownames(scores) <- subjects
  colnames(scores) <- colnames(residual_x)
  scores
}

scalar_cluster_tests <- function(scores) {
  scores <- as.matrix(scores)
  n_subjects <- nrow(scores)
  if (n_subjects < 2L) {
    stop("At least two subject scores are required.", call. = FALSE)
  }

  means <- colMeans(scores)
  variances <- colMeans(sweep(scores, 2L, means, "-")^2)
  statistics <- sqrt(n_subjects) * means / sqrt(variances)
  p_values <- 2 * stats::pnorm(abs(statistics), lower.tail = FALSE)

  data.frame(
    marker = colnames(scores),
    n_subjects = n_subjects,
    estimate = means,
    standard_error = sqrt(variances / n_subjects),
    statistic = statistics,
    p_value = p_values,
    p_holm = stats::p.adjust(p_values, method = "holm"),
    row.names = NULL
  )
}

global_multiplier_test <- function(scores, bootstrap_reps = 9999L,
                                   seed = 20260726L) {
  scores <- as.matrix(scores)
  n_subjects <- nrow(scores)
  centered <- sweep(scores, 2L, colMeans(scores), "-")
  scale_values <- sqrt(colMeans(centered^2))
  if (any(!is.finite(scale_values) | scale_values <= 0)) {
    stop("Every marker score must have positive finite variance.", call. = FALSE)
  }

  observed <- max(abs(sqrt(n_subjects) * colMeans(scores) / scale_values))
  set.seed(seed)
  bootstrap_statistics <- numeric(bootstrap_reps)
  for (b in seq_len(bootstrap_reps)) {
    multiplier <- stats::rnorm(n_subjects)
    bootstrap_z <- colSums(centered * multiplier) /
      sqrt(n_subjects) / scale_values
    bootstrap_statistics[b] <- max(abs(bootstrap_z))
  }

  list(
    statistic = observed,
    p_value = (1 + sum(bootstrap_statistics >= observed)) /
      (bootstrap_reps + 1),
    bootstrap_reps = bootstrap_reps
  )
}

nuisance_diagnostics <- function(data, x_markers, y_outcome,
                                 x_fits, y_fit, folds, id_col) {
  rows <- lapply(
    c(x_markers, y_outcome),
    function(outcome) {
      fit <- if (outcome == y_outcome) y_fit else x_fits[[outcome]]
      data.frame(
        outcome = outcome,
        rmse = sqrt(mean(fit$residual^2)),
        residual_mean = mean(fit$residual),
        residual_sd = stats::sd(fit$residual),
        stringsAsFactors = FALSE
      )
    }
  )
  list(
    outcomes = do.call(rbind, rows),
    fold_balance = data.frame(
      fold = sort(unique(folds)),
      subjects = as.integer(tapply(
        as.character(data[[id_col]]), folds,
        function(x) length(unique(x))
      )),
      visits = as.integer(table(folds))
    )
  )
}

run_design_a <- function(data, x_markers, y_outcome, rhs, id_col,
                         subject_folds = 5L, bootstrap_reps = 9999L,
                         seed = 20260726L,
                         equal_subject_nuisance_weights = TRUE) {
  assert_columns(data, c(id_col, x_markers, y_outcome))
  data$.subject_re <- factor(as.character(data[[id_col]]))
  folds <- make_subject_folds(data[[id_col]], k = subject_folds, seed = seed)
  case_weights <- if (isTRUE(equal_subject_nuisance_weights)) {
    make_equal_subject_weights(data[[id_col]])
  } else {
    rep(1, nrow(data))
  }

  y_fit <- crossfit_residuals(
    data, y_outcome, rhs, folds, prediction = "population",
    case_weights = case_weights
  )
  x_fits <- setNames(lapply(
    x_markers,
    function(marker) {
      crossfit_residuals(
        data, marker, rhs, folds, prediction = "population",
        case_weights = case_weights
      )
    }
  ), x_markers)
  x_residuals <- vapply(
    x_fits, function(fit) fit$residual, numeric(nrow(data))
  )
  colnames(x_residuals) <- x_markers
  scores <- subject_score_matrix(
    data[[id_col]], x_residuals, y_fit$residual
  )

  list(
    design = "A",
    target = "population-level X_it independent of Y_it given Z_it",
    inferential_unit = "subject",
    validity = "primary cluster-cross-fitted analysis",
    nuisance_weighting = if (isTRUE(equal_subject_nuisance_weights)) {
      "equal subject"
    } else {
      "equal visit"
    },
    folds = folds,
    scores = scores,
    marker_tests = scalar_cluster_tests(scores),
    global_test = global_multiplier_test(scores, bootstrap_reps, seed + 1L),
    diagnostics = nuisance_diagnostics(
      data, x_markers, y_outcome, x_fits, y_fit, folds, id_col
    ),
    residuals_x = x_residuals,
    residuals_y = y_fit$residual,
    predictions_x = vapply(
      x_fits, function(fit) fit$prediction, numeric(nrow(data))
    ),
    predictions_y = y_fit$prediction
  )
}

run_design_b <- function(data, x_markers, y_outcome, rhs, id_col,
                         visit_folds = 2L, min_visits = 3L,
                         bootstrap_reps = 9999L, seed = 20260726L) {
  assert_columns(data, c(id_col, x_markers, y_outcome))
  visit_counts <- table(as.character(data[[id_col]]))
  eligible <- names(visit_counts[visit_counts >= max(min_visits, visit_folds)])
  analysis_data <- data[as.character(data[[id_col]]) %in% eligible, , drop = FALSE]
  if (length(eligible) < 2L) {
    stop("Design B needs at least two subjects with sufficient repeated visits.",
         call. = FALSE)
  }
  analysis_data$.subject_re <- factor(as.character(analysis_data[[id_col]]))
  folds <- make_within_subject_folds(
    analysis_data[[id_col]], k = visit_folds, seed = seed
  )

  y_fit <- crossfit_residuals(
    analysis_data, y_outcome, rhs, folds, prediction = "subject"
  )
  x_fits <- setNames(lapply(
    x_markers,
    function(marker) {
      crossfit_residuals(
        analysis_data, marker, rhs, folds, prediction = "subject"
      )
    }
  ), x_markers)
  x_residuals <- vapply(
    x_fits, function(fit) fit$residual, numeric(nrow(analysis_data))
  )
  colnames(x_residuals) <- x_markers
  scores <- subject_score_matrix(
    analysis_data[[id_col]], x_residuals, y_fit$residual
  )

  list(
    design = "B",
    target = "model-based X_it independent of Y_it given Z_it and latent U_i",
    inferential_unit = "subject",
    validity = paste(
      "sensitivity analysis only:",
      "within-subject splitting removes self-influence but does not identify U_i",
      "when visit counts are bounded"
    ),
    excluded_subjects = setdiff(names(visit_counts), eligible),
    folds = folds,
    scores = scores,
    marker_tests = scalar_cluster_tests(scores),
    global_test = global_multiplier_test(scores, bootstrap_reps, seed + 1L),
    diagnostics = nuisance_diagnostics(
      analysis_data, x_markers, y_outcome, x_fits, y_fit, folds, id_col
    ),
    residuals_x = x_residuals,
    residuals_y = y_fit$residual,
    predictions_x = vapply(
      x_fits, function(fit) fit$prediction, numeric(nrow(analysis_data))
    ),
    predictions_y = y_fit$prediction
  )
}
