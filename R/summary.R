#' Summarise an HMM-SSF fit
#'
#' For a guided fit, the summary lists the candidates that were tried and
#' gives the ecological labels of the selected fit. A fit from [fit_ssf()]
#' uses its numerical state names. The state table counts Viterbi assignments.
#' Kernel means and signed directional persistence are included when the
#' movement terms and their density convention are available.
#'
#' @param object An object from [fit_ssf()] or [fit_guided_ssf()].
#' @param x A fitted object or a summary returned by one of these methods.
#' @param assessment Optional accepted assessment of the same fit, used for
#'   ecological state labels.
#' @param term_labels Optional named character vector of readable names for
#'   additional covariates, including reference categories when known.
#'   Factor coding is not reconstructed from the coefficient names.
#' @param step_scale Positive multiplier for displayed mean step lengths.
#' @param step_unit Distance unit after applying `step_scale`. Units are not inferred.
#' @param ... Reserved for the generic.
#' @param digits Number of significant digits to display.
#' @return A list with a state table (`states`), step-selection estimates
#'   and standard errors (`coefficients`), additional covariate effects with
#'   conditional 95% Wald intervals (`ssf_effects`), transition probabilities
#'   (`transition`) and optimisation diagnostics. Kernel means are in the
#'   chosen step-length units and describe the baseline kernel when there are
#'   lag or hourly effects. Standard errors are available when the fit has
#'   a covariance matrix. Daily transition arrays retain the hour dimension.
#'   Guided summaries also contain `status`, `taxonomy` and `candidates`.
#'   The candidate table gives each code, its outcome and a description
#'   of its states when available.
#'   When no candidate is selected, only these three components are returned.
#' @name summary.ssf_fit
#' @export
summary.ssf_fit <- function(object, assessment = NULL, term_labels = NULL,
                            step_scale = 1, step_unit = "input units", ...) {
  resolved <- .state_display_fit(object, assessment)
  .summarise_ssf(resolved$fit, resolved$labels, term_labels, step_scale, step_unit)
}

#' @rdname summary.ssf_fit
#' @export
summary.guided_ssf <- function(object, assessment = NULL, term_labels = NULL,
                               step_scale = 1, step_unit = "input units", ...) {
  out <- if (is.null(object$fit)) list() else {
    resolved <- .state_display_fit(object, assessment)
    .summarise_ssf(resolved$fit, resolved$labels, term_labels, step_scale, step_unit)
  }
  out["taxonomy"] <- list(object$taxonomy)
  out$status <- object$status
  out$candidates <- .guided_candidates(object)
  class(out) <- "summary_guided_ssf"
  out
}

.guided_candidates <- function(object) {
  history <- object$history
  codes <- names(history)
  description <- vapply(codes, function(code) {
    if (code == "Other") return("Single state")
    labels <- .reference_state_labels[[code]]
    if (is.null(labels)) labels <- history[[code]]$assessment$labels
    if (code == "HDT_F_R") labels[1L] <- "Habitat-directed travelling"
    paste(labels[!is.na(labels)], collapse = " / ")
  }, character(1))
  status <- vapply(history, function(a) {
    if (!is.null(a$error)) return("Fit failed")
    if (isTRUE(a$assessment$accepted)) "Selected" else "Rejected"
  }, character(1))
  reason <- vapply(history, function(a) {
    if (!is.null(a$error)) return(a$error)
    reasons <- a$assessment$reasons
    if (identical(a$assessment$taxonomy, "HDT_F_R")) {
      reasons[reasons == "target_positive"] <-
        "The 95% confidence interval for the travelling state's target effect must lie entirely above zero."
    }
    paste(reasons, collapse = ", ")
  }, character(1))
  data.frame(taxonomy = codes, status = unname(status), reason = unname(reason),
             description = unname(description), stringsAsFactors = FALSE)
}

.print_guided_candidates <- function(candidates) {
  cat("Candidate taxonomies in the order tried\n")
  for (i in seq_len(nrow(candidates))) {
    description <- candidates$description[i]
    cat("  ", candidates$taxonomy[i],
        if (nzchar(description)) paste0(" (", description, ")"),
        " [", candidates$status[i], "]\n", sep = "")
    if (nzchar(candidates$reason[i])) {
      lines <- strwrap(candidates$reason[i], width = max(25L, getOption("width") - 4L))
      cat(paste0("    ", lines), sep = "\n")
    }
  }
}

.ssf_term_labels <- function(terms, term_labels = NULL) {
  pretty <- stats::setNames(gsub("_", " ", terms, fixed = TRUE), terms)
  pretty[terms == "cos_angle"] <- "Previous heading"
  pretty[terms == "cos_diff_theta_to_target"] <- "Direction to target"
  if (!is.null(term_labels)) {
    if (!is.character(term_labels) || is.null(names(term_labels)) || anyNA(term_labels) ||
        anyDuplicated(names(term_labels)) || !all(names(term_labels) %in% terms))
      .stop("term_labels must be named by displayed coefficient terms.")
    pretty[names(term_labels)] <- term_labels
  }
  pretty
}

.summarise_ssf <- function(fit, labels, term_labels = NULL, step_scale = 1, step_unit = "input units") {
  .scalar(step_scale, "step_scale", .Machine$double.eps)
  if (!is.character(step_unit) || length(step_unit) != 1L || is.na(step_unit) || !nzchar(step_unit))
    .stop("step_unit must be a nonempty unit label.")
  b <- fit$coefficients
  n_steps <- length(fit$state)
  counts <- tabulate(fit$state, nbins = fit$n_states)
  states <- data.frame(state = labels, n_steps = counts, fraction = counts / n_steps)
  if (isTRUE(fit$density_measure %in% c("area", "step_angle")) &&
      (!is.null(fit$step_family) || all(c("step", "log_step", "cos_angle") %in% rownames(b)))) {
    metrics <- movement_metrics(b, fit$density_measure, fit$step_family,
      if (is.null(fit$step_max)) Inf else fit$step_max)
    states$mean_step <- metrics$mean_step * step_scale
    states$persistence <- metrics$persistence
  }
  standard_error <- rep(NA_real_, length(b))
  covariance <- fit[["covariance"]]
  if (!is.null(covariance)) {
    variance <- diag(covariance)[seq_along(b)]
    valid <- is.finite(variance) & variance >= 0
    standard_error[valid] <- sqrt(variance[valid])
  }
  coefficients <- data.frame(
    state = rep(labels, each = nrow(b)), term = rep(rownames(b), times = ncol(b)),
    estimate = as.vector(b), std_error = standard_error
  )
  terms <- .state_effect_terms(fit)
  effects <- coefficients[coefficients$term %in% terms, , drop = FALSE]
  effects$label <- unname(.ssf_term_labels(terms, term_labels)[effects$term])
  effects$lower <- effects$estimate - stats::qnorm(.975) * effects$std_error
  effects$upper <- effects$estimate + stats::qnorm(.975) * effects$std_error
  transition <- fit$transition
  dimnames(transition)[1:2] <- list(labels, labels)
  structure(list(
    states = states, coefficients = coefficients, ssf_effects = effects, transition = transition,
    step_unit = step_unit,
    n_steps = n_steps, n_segments = length(unique(fit$observations$ID)),
    normalization = fit$normalization, density_measure = fit$density_measure,
    step_family = fit$step_family, log_likelihood = fit$log_likelihood, diagnostics = fit$diagnostics
  ), class = "summary_ssf_fit")
}

# Print a compact ecological reading first; the full summary retains numerical
# movement parameters and transitions for users who need to inspect the fit.
.print_ssf_overview <- function(x, digits = 3L) {
  cat("Fitted HMM-SSF\n")
  cat(nrow(x$states), "state(s) |", format(x$n_steps, big.mark = ","), "decoded steps |",
      x$n_segments, "segment(s)\n")
  cat("Step-length family:", if (is.null(x$step_family)) "gamma" else x$step_family,
      "| normalisation:", x$normalization, "| measure:", x$density_measure, "\n")
  cat("Log likelihood:", format(x$log_likelihood, digits = max(7L, digits)), "\n")
  d <- x$diagnostics
  cat("Converged:", d$converged, "| positive Hessian:", d$hessian_positive,
      "| scaled gradient:", format(if (is.null(d$max_abs_scaled_gradient)) NA_real_ else d$max_abs_scaled_gradient, digits = digits), "\n")
  cat("\nStates (decoded steps)\n")
  state_table <- data.frame(State = x$states$state, Steps = x$states$n_steps,
    Share = sprintf("%.1f%%", 100 * x$states$fraction), check.names = FALSE)
  if ("mean_step" %in% names(x$states)) {
    state_table[[paste0("Mean step (", x$step_unit, ")")]] <- signif(x$states$mean_step, digits)
  }
  print(state_table, row.names = FALSE, right = FALSE)
  cat("\nAdditional SSF covariates by state\n")
  if (!nrow(x$ssf_effects)) {
    cat("  No additional covariates in this model.\n")
  } else {
    cat("  Estimate [95% Wald CI]\n")
    for (state in x$states$state) {
      cat("\n  ", state, "\n", sep = "")
      z <- x$ssf_effects[x$ssf_effects$state == state, , drop = FALSE]
      values <- function(v) formatC(v, format = "fg", digits = digits, flag = "#")
      for (i in seq_len(nrow(z))) {
        interval <- if (is.finite(z$lower[i]) && is.finite(z$upper[i]))
          paste0("[", values(z$lower[i]), ", ", values(z$upper[i]), "]") else "[CI unavailable]"
        cat("    ", z$label[i], ": ", values(z$estimate[i]), " ", interval, "\n", sep = "")
      }
    }
    cat("\n  Effects condition on the other model terms.\n")
    cat(paste0("  ", strwrap("CIs condition on the fitted model and labels; they do not account for their selection.",
      width = max(25L, getOption("width") - 2L))), sep = "\n")
    cat("  Use ssf_state_contrasts() to compare coefficients between states.\n")
  }
  if ("mean_step" %in% names(x$states)) {
    cat("\n")
    cat(strwrap("Mean steps describe the intrinsic kernel; hourly and lag effects are held at baseline.",
      width = max(25L, getOption("width"))), sep = "\n")
  }
  invisible(x)
}

#' @rdname summary.ssf_fit
#' @export
print.summary_ssf_fit <- function(x, digits = 4L, ...) {
  .print_ssf_overview(x, digits)
  baseline <- x$coefficients[!x$coefficients$term %in% x$ssf_effects$term, , drop = FALSE]
  if (nrow(baseline)) {
    cat("\nBaseline movement coefficients (estimate and standard error)\n")
    print(baseline, digits = digits, row.names = FALSE)
  }
  if (length(dim(x$transition)) == 2L) {
    cat("\nTransition probabilities (rows are origin states)\n")
    print(x$transition, digits = digits)
  } else {
    cat("\nSelf-transition probabilities by state (range across hours)\n")
    stay <- vapply(seq_len(nrow(x$transition)), function(i)
      range(x$transition[i, i, ]), numeric(2))
    dimnames(stay) <- list(c("Min", "Max"), x$states$state)
    print(t(stay), digits = digits)
  }
  invisible(x)
}

#' @rdname summary.ssf_fit
#' @export
print.summary_guided_ssf <- function(x, digits = 4L, ...) {
  if (is.null(x$taxonomy)) cat("No candidate passed the guided assessment.\n")
  .print_guided_candidates(x$candidates)
  if (!is.null(x$taxonomy)) {
    cat("\n")
    print.summary_ssf_fit(x, digits = digits, ...)
  }
  invisible(x)
}
