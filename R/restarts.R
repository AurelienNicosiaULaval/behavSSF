#' Diagnose sensitivity to HMM-SSF starting values
#'
#' Fit the same model from different starting values. Each run uses the same
#' steps and control sets, if supplied. The restart log shows which fits pass
#' the taxonomy checks and whether the optimiser agrees with the forward
#' likelihood calculation.
#'
#' The search groups fits that give the same Viterbi labels and have similar
#' likelihoods and parameters once their states are aligned by label. It works
#' from highest to lowest likelihood and compares each fit to the existing
#' representatives in that order. The parameter tolerance applies to the
#' largest absolute difference on the fitting scale, including log transition
#' odds relative to self-transition. Restart counts depend on the starting
#' values tried.
#'
#' @param data Observed `ssf_steps` or case-control data for [fit_ssf()].
#'   Case-control data must have a density measure declared.
#' @param taxonomy A specification from [ssf_taxonomy()] or a built-in name
#'   accepted by [assess_taxonomy()].
#' @param n_starts Number of random starts, used when `starts` is NULL.
#' @param seed Optional seed for the whole search. The previous RNG state is
#'   restored on exit.
#' @param spread,short_iter,long_iter,engine Fitting controls for [fit_ssf()].
#' @param normalization Defaults to `"exact"` for `ssf_steps` and
#'   `"importance"` for case-control data.
#' @param starts Optional list of start specifications that replaces the
#'   default random starts. Each entry is a list with optional `seed`,
#'   `spread` and `initial_values` arguments for [fit_ssf()]. An empty list
#'   entry uses the default random initialisation. Names of entries are
#'   optional group labels and may repeat.
#' @param formula Optional formula override. Defaults to the custom taxonomy
#'   formula or the built-in movement formula with a target term for HDT.
#' @param target Name of the target coefficient for the built-in HDT taxonomy.
#' @param assessment_args Optional named list of thresholds for
#'   [assess_taxonomy()], for example `list(min_occupancy = 0.10)`.
#' @param likelihood_tol,parameter_tol Absolute grouping tolerances for
#'   log likelihood and aligned working parameters.
#' @param objective_tol Maximum allowed difference between the optimiser
#'   and forward log likelihood.
#' @return An `ssf_restarts` object with `solutions` ranked by likelihood,
#'   `states` with movement and occupancy summaries for each solution,
#'   `restarts` with one outcome per start, `attempts` with optimiser logs,
#'   and `fits` with the representative `ssf_fit` objects. Solution numbers
#'   are scoped to this search. Use [compare_ssf_solutions()] for step labels.
#' @seealso [compare_ssf_solutions()], [plot_ssf_changes()]
#' @examples
#' gps <- readRDS(system.file("extdata", "deer-restart-steps.rds",
#'   package = "behavSSF"))
#' steps <- prepare_ssf_steps(gps, density_measure = "area")
#' labels <- ssf_taxonomy("step_lengths", c("Short steps", "Long steps"))
#' check <- diagnose_ssf_restarts(steps, labels, n_starts = 6, seed = 60200)
#' check$solutions
#' @export
diagnose_ssf_restarts <- function(data, taxonomy, n_starts = 30L, seed = NULL,
                                  spread = 0.5, short_iter = 100L, long_iter = 900L,
                                  engine = c("native", "R"),
                                  normalization = if (inherits(data, "ssf_steps")) "exact" else "importance",
                                  starts = NULL, formula = NULL,
                                  target = "cos_diff_theta_to_target", assessment_args = list(),
                                  likelihood_tol = 1e-5, parameter_tol = 1e-3,
                                  objective_tol = 1e-7) {
  engine <- match.arg(engine)
  normalization <- match.arg(normalization, c("exact", "importance"))
  custom <- inherits(taxonomy, "ssf_taxonomy")
  if (!custom) taxonomy <- match.arg(taxonomy, names(.reference_state_counts))
  k <- if (custom) length(taxonomy$labels) else .reference_state_counts[[taxonomy]]
  taxonomy_name <- if (custom) taxonomy$name else taxonomy
  .scalar(short_iter, "short_iter", 1, 1e6, TRUE)
  .scalar(long_iter, "long_iter", 1, 1e6, TRUE)
  .scalar(likelihood_tol, "likelihood_tol", 0)
  .scalar(parameter_tol, "parameter_tol", 0)
  .scalar(objective_tol, "objective_tol", 0)
  if (is.null(formula)) {
    terms <- c("step", "log_step", "cos_angle")
    formula <- if (custom) taxonomy$formula else stats::reformulate(
      c(terms, if (taxonomy == "HDT_F_R") target))
  }
  exact <- normalization == "exact"
  design <- if (exact) .design_exact(data, formula) else .design(data, formula)
  convention <- attr(design$data, "density_measure")
  if (length(convention) != 1L || !convention %in% c("area", "step_angle"))
    .stop("Declare the sampling density measure before diagnosing restarts.")
  if (!all(c("step", "log_step", "cos_angle") %in% colnames(design$X)))
    .stop("The fitted design must include step, log_step, and cos_angle.")
  if (!custom && taxonomy == "HDT_F_R" &&
      (length(target) != 1L || !target %in% colnames(design$X)))
    .stop("The HDT taxonomy requires a target coefficient in the fitted design.")
  thresholds <- c("min_occupancy", "persistence_bounds", "min_mean_gap", "score_tolerance")
  if (!is.list(assessment_args) || (length(assessment_args) &&
      (is.null(names(assessment_args)) || anyDuplicated(names(assessment_args)) ||
       !all(names(assessment_args) %in% thresholds))))
    .stop("assessment_args must be a named list of assessment thresholds.")
  if (is.null(starts)) {
    .scalar(n_starts, "n_starts", 1, 10000, TRUE)
    starts <- rep(list(list()), n_starts)
  }
  if (!is.list(starts) || !length(starts) || !all(vapply(starts, function(s)
    is.list(s) && (!length(s) || (!is.null(names(s)) && !anyDuplicated(names(s)) &&
      all(names(s) %in% c("seed", "spread", "initial_values")))), logical(1))))
    .stop("starts must be a nonempty list of start specifications.")
  group <- names(starts)
  if (is.null(group)) group <- rep("", length(starts))
  group[is.na(group) | !nzchar(group)] <- "Start"

  attempts <- .with_seed(seed, function() lapply(unname(starts), function(start) {
    fit <- tryCatch(
      fit_ssf(data, formula, n_states = k, n_starts = 1L,
        short_iter = short_iter, long_iter = long_iter,
        spread = if (is.null(start$spread)) spread else start$spread,
        seed = start$seed, engine = engine, initial_values = start$initial_values,
        normalization = normalization),
      error = identity
    )
    if (inherits(fit, "error"))
      return(list(reason = conditionMessage(fit), log = fit$attempts))
    assessment <- do.call(assess_taxonomy, c(
      list(fit = fit, taxonomy = taxonomy, convention = convention, target = target),
      assessment_args))
    reason <- assessment$reasons
    if (!isTRUE(fit$diagnostics$objective_agreement <= objective_tol))
      reason <- c(reason, "objective_agreement")
    attempt <- list(log_likelihood = fit$log_likelihood, reason = reason, log = fit$starts)
    if (length(reason)) return(attempt)

    labels <- assessment$labels
    label_order <- if (custom) match(taxonomy$labels, labels) else order(labels, method = "radix")
    coefficients <- c(fit$coefficients[, label_order, drop = FALSE])
    if (exact) coefficients <- .exact_transform(coefficients, design, k, inverse = TRUE)
    transition <- fit$transition[label_order, label_order, drop = FALSE]
    # Default admissibility bounds exclude degenerate self-transition. A
    # zero off-diagonal entry can still occur with three or more states.
    log_odds <- log(transition) - log(diag(transition))
    attempt$signature <- c(coefficients, log_odds[row(transition) != col(transition)])
    attempt$labels <- labels[fit$state]
    attempt$label_order <- label_order
    attempt$fit <- fit
    attempt$assessment <- assessment
    attempt
  }))

  log_likelihood <- vapply(attempts, function(a)
    if (is.null(a$log_likelihood)) NA_real_ else a$log_likelihood, numeric(1))
  eligible <- vapply(attempts, function(a) !length(a$reason), logical(1))
  solution <- rep(NA_integer_, length(attempts))
  representatives <- integer()
  for (i in which(eligible)[order(log_likelihood[eligible], decreasing = TRUE)]) {
    a <- attempts[[i]]
    same <- vapply(representatives, function(j) {
      b <- attempts[[j]]
      if (abs(log_likelihood[i] - log_likelihood[j]) > likelihood_tol ||
          !identical(a$labels, b$labels)) return(FALSE)
      difference <- abs(a$signature - b$signature)
      difference[a$signature == b$signature] <- 0
      max(difference) <= parameter_tol
    }, logical(1))
    if (!any(same)) representatives <- c(representatives, i)
    solution[i] <- if (any(same)) match(TRUE, same) else length(representatives)
  }
  rank <- seq_along(representatives)
  best <- representatives[1L]
  solutions <- data.frame(
    solution = rank, representative_start = representatives,
    log_likelihood = log_likelihood[representatives],
    delta_log_likelihood = log_likelihood[best] - log_likelihood[representatives],
    n_restarts = tabulate(solution, nbins = length(rank)),
    n_changed_vs_best = vapply(representatives, function(i)
      sum(attempts[[i]]$labels != attempts[[best]]$labels), integer(1))
  )
  states <- data.frame(solution = integer(), label = character(), n_steps = integer(),
    occupancy = numeric(), mean_step = numeric(), persistence = numeric(),
    self_transition = numeric())
  for (s in rank) {
    a <- attempts[[representatives[s]]]
    assessment <- a$assessment
    labels <- assessment$labels
    row_order <- a$label_order
    metrics <- assessment$metrics[row_order, , drop = FALSE]
    states <- rbind(states, data.frame(solution = s, label = labels[row_order],
      n_steps = tabulate(a$fit$state, nbins = k)[row_order],
      occupancy = assessment$occupancy[row_order], mean_step = metrics$mean_step,
      persistence = metrics$persistence, self_transition = diag(a$fit$transition)[row_order]))
  }
  structure(list(
    solutions = solutions, states = states,
    restarts = data.frame(start = seq_along(starts), group = group, solution = solution,
      status = ifelse(eligible, "admissible", ifelse(is.na(log_likelihood), "failed", "rejected")),
      log_likelihood = log_likelihood,
      reason = vapply(attempts, function(a) paste(a$reason, collapse = ", "), character(1))),
    attempts = lapply(attempts, `[[`, "log"),
    fits = lapply(attempts[representatives], `[[`, "fit"),
    assessments = lapply(attempts[representatives], `[[`, "assessment"),
    taxonomy = taxonomy_name
  ), class = "ssf_restarts")
}

#' @export
print.ssf_restarts <- function(x, ...) {
  cat("HMM-SSF restart diagnostic (", x$taxonomy, ")\n", sep = "")
  cat(nrow(x$restarts), " starts, ", nrow(x$solutions), " solutions, ",
      sum(x$restarts$status != "admissible"), " unsuccessful\n", sep = "")
  if (nrow(x$solutions)) {
    display <- x$solutions[c("solution", "n_restarts", "log_likelihood", "n_changed_vs_best")]
    names(display) <- c("Solution", "Starts", "Log likelihood", "Changes vs 1")
    print(display, row.names = FALSE)
  }
  invisible(x)
}

#' Compare step labels between two restart solutions
#'
#' @param diagnostic Output from [diagnose_ssf_restarts()].
#' @param from,to Distinct solution numbers from `diagnostic$solutions`.
#' @return An `ssf_solution_comparison` with a `steps` data frame in fitted
#'   observation order, a `transitions` count table, and `posterior_from`
#'   and `posterior_to` matrices with columns named by taxonomy label.
#'   The step table preserves observed columns, including coordinates and
#'   timestamps when supplied. Added columns are `.step_index`,
#'   `.segment_step`, `.state_from`, `.state_to` and `.changed`.
#'   `.source_row` identifies the row in the input as prepared by the
#'   fitting workflow. Keep original fix IDs in separate input columns.
#' @seealso [diagnose_ssf_restarts()], [plot_ssf_changes()]
#' @export
compare_ssf_solutions <- function(diagnostic, from = 1L, to = 2L) {
  if (!inherits(diagnostic, "ssf_restarts")) .stop("Use diagnose_ssf_restarts() first.")
  n <- nrow(diagnostic$solutions)
  if (n < 2L) .stop("At least two admissible solutions are required for comparison.")
  .scalar(from, "from", 1, n, TRUE)
  .scalar(to, "to", 1, n, TRUE)
  if (from == to) .stop("Choose two distinct solutions.")
  a <- diagnostic$fits[[from]]
  b <- diagnostic$fits[[to]]
  labels_from <- diagnostic$assessments[[from]]$labels
  labels_to <- diagnostic$assessments[[to]]$labels
  posterior_from <- a$posterior
  posterior_to <- b$posterior
  colnames(posterior_from) <- labels_from
  colnames(posterior_to) <- labels_to
  levels <- diagnostic$states$label[diagnostic$states$solution == from]
  steps <- as.data.frame(a$observations)
  added <- c(".step_index", ".segment_step", ".state_from", ".state_to", ".changed")
  if (any(added %in% names(steps))) .stop("Step data contain reserved comparison columns.")
  steps$.step_index <- seq_len(nrow(steps))
  steps$.segment_step <- sequence(rle(as.character(steps$ID))$lengths)
  steps$.state_from <- factor(labels_from[a$state], levels = levels)
  steps$.state_to <- factor(labels_to[b$state], levels = levels)
  steps$.changed <- steps$.state_from != steps$.state_to
  structure(list(from = from, to = to, steps = steps,
    transitions = table(from = steps$.state_from, to = steps$.state_to),
    posterior_from = posterior_from[, levels, drop = FALSE],
    posterior_to = posterior_to[, levels, drop = FALSE]), class = "ssf_solution_comparison")
}

#' Map steps whose state labels change between restart solutions
#'
#' @param comparison,x Output from [compare_ssf_solutions()].
#' @param ... Reserved for the generic.
#' @param coords Four column names for start x, start y, end x and end y,
#'   respectively. Use projected coordinates in common units.
#' @param segment Optional segment identifiers from the `ID` column.
#' @return A ggplot object. Every observed step is drawn as its own line.
#'   A coloured endpoint identifies each relabelled step, including steps
#'   whose line is too short to see at the map scale.
#' @seealso [compare_ssf_solutions()]
#' @export
plot_ssf_changes <- function(comparison, coords = c("x0", "y0", "x", "y"), segment = NULL) {
  if (!inherits(comparison, "ssf_solution_comparison"))
    .stop("Use compare_ssf_solutions() first.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    .stop("Install ggplot2 to draw the trajectory.")
  steps <- comparison$steps
  if (!is.null(segment)) {
    if (!length(segment) || !all(segment %in% steps$ID))
      .stop("segment must identify fitted trajectory segments.")
    steps <- steps[steps$ID %in% segment, , drop = FALSE]
  }
  if (!.has_ssf_coordinates(steps, coords))
    .stop("coords must identify four finite numeric columns for step starts and ends.")
  track <- stats::setNames(steps[coords], c("x0", "y0", "x", "y"))
  track$status <- factor(ifelse(steps$.changed, "Changed", "Unchanged"),
    levels = c("Changed", "Unchanged"))
  track <- track[order(steps$.changed), , drop = FALSE]
  ggplot2::ggplot(track) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$x0, y = .data$y0,
      xend = .data$x, yend = .data$y, colour = .data$status, linewidth = .data$status),
      lineend = "round") +
    ggplot2::geom_point(data = track[track$status == "Changed", , drop = FALSE],
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$status), size = 1.5) +
    ggplot2::scale_colour_manual(values = c(Changed = "#AE3E1B", Unchanged = "#859AA5"),
      name = NULL, drop = FALSE) +
    ggplot2::scale_linewidth_manual(values = c(Changed = .85, Unchanged = .48), guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = "Where step labels change",
      subtitle = sprintf("Solution %d to %d  |  %d of %d steps relabelled",
        comparison$from, comparison$to, sum(steps$.changed), nrow(steps)),
      x = "Easting", y = "Northing") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
      text = ggplot2::element_text(colour = "#183345"))
}

#' @rdname plot_ssf_changes
#' @export
plot.ssf_solution_comparison <- function(x, ..., coords = c("x0", "y0", "x", "y"),
                                         segment = NULL) {
  plot_ssf_changes(x, coords = coords, segment = segment)
}

utils::globalVariables(".data")
