#' Summarise an HMM-SSF fit
#'
#' For a guided fit, the summary lists the candidates that were tried and
#' gives the ecological labels of the selected fit. A fit from [fit_ssf()]
#' uses its numerical state names. The state table counts Viterbi assignments.
#' Kernel means and signed directional persistence are included when the
#' movement terms and their density convention are available.
#'
#' @param object An object from [fit_ssf()] or [fit_guided_ssf()].
#' @param x A summary returned by one of these methods.
#' @param ... Reserved for the generic.
#' @param digits Number of significant digits to display.
#' @return A list with a state table (`states`), step-selection estimates
#'   and standard errors (`coefficients`), transition probabilities
#'   (`transition`) and optimisation diagnostics. Kernel means are in the
#'   input step-length units and describe the baseline kernel when there are
#'   lag or hourly effects. Standard errors are available when the fit has
#'   a covariance matrix. Daily transition arrays retain the hour dimension.
#'   Guided summaries also contain `status`, `taxonomy` and `candidates`.
#'   The candidate table gives each code, its outcome and a description
#'   of its states when available.
#'   When no candidate is selected, only these three components are returned.
#' @export
summary.ssf_fit <- function(object, ...) {
  .summarise_ssf(object, colnames(object$coefficients))
}

#' @rdname summary.ssf_fit
#' @export
summary.guided_ssf <- function(object, ...) {
  out <- if (is.null(object$fit)) list() else
    .summarise_ssf(object$fit, object$assessment$labels)
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

.summarise_ssf <- function(fit, labels) {
  b <- fit$coefficients
  n_steps <- length(fit$state)
  counts <- tabulate(fit$state, nbins = fit$n_states)
  states <- data.frame(state = labels, n_steps = counts, fraction = counts / n_steps)
  if (isTRUE(fit$density_measure %in% c("area", "step_angle")) &&
      (!is.null(fit$step_family) || all(c("step", "log_step", "cos_angle") %in% rownames(b)))) {
    metrics <- movement_metrics(b, fit$density_measure, fit$step_family,
      if (is.null(fit$step_max)) Inf else fit$step_max)
    states$mean_step <- metrics$mean_step
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
  transition <- fit$transition
  dimnames(transition)[1:2] <- list(labels, labels)
  structure(list(
    states = states, coefficients = coefficients, transition = transition,
    n_steps = n_steps, n_segments = length(unique(fit$observations$ID)),
    normalization = fit$normalization, density_measure = fit$density_measure,
    step_family = fit$step_family, log_likelihood = fit$log_likelihood, diagnostics = fit$diagnostics
  ), class = "summary_ssf_fit")
}

#' @rdname summary.ssf_fit
#' @export
print.summary_ssf_fit <- function(x, digits = 4L, ...) {
  cat("HMM-SSF summary\n")
  cat(x$n_steps, "observed steps in", x$n_segments, "segment(s)\n")
  cat("Normalisation:", x$normalization, "| density measure:", x$density_measure, "\n")
  if (!is.null(x$step_family)) cat("Step-length family:", x$step_family, "\n")
  cat("Log likelihood:", format(x$log_likelihood, digits = digits), "\n")
  d <- x$diagnostics
  cat("Converged:", d$converged, "| positive Hessian:", d$hessian_positive,
      "| scaled gradient:", format(if (is.null(d$max_abs_scaled_gradient)) NA_real_ else d$max_abs_scaled_gradient, digits = digits), "\n")
  cat("\nStates (Viterbi assignments)\n")
  print(x$states, digits = digits, row.names = FALSE)
  cat("\nStep-selection coefficients\n")
  print(x$coefficients, digits = digits, row.names = FALSE)
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
