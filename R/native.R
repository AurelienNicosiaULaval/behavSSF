#' @useDynLib behavSSF, .registration = TRUE, .fixes = "C_"
NULL

.native_objective <- function(design, n_states) {
  x <- design$X
  storage.mode(x) <- "double"
  starts <- as.integer(vapply(design$groups, `[`, integer(1), 1L) - 1L)
  sizes <- as.integer(lengths(design$groups))
  segment <- as.integer(factor(design$segment, levels = unique(design$segment)))
  .cached_objective(function(theta) .Call(C_behavSSF_objective, as.double(theta), x,
    design$logw, starts, sizes, segment, as.integer(n_states)))
}

.native_exact_objective <- function(design, n_states) {
  x <- design$X
  storage.mode(x) <- "double"
  segment <- as.integer(factor(design$segment, levels = unique(design$segment)))
  context <- list(design$indices, as.double(design$relative_bearing),
                  as.double(design$target_active), as.double(design$shape_offset))
  if (!is.null(design$diel) && n_states == 2L) context <- c(context, design$diel)
  .cached_objective(function(theta) .Call(C_behavSSF_exact_objective, as.double(theta), x,
    segment, as.integer(n_states), context))
}

.cached_objective <- function(compute) {
  last_theta <- NULL
  last_value <- NULL
  invalid <- 0L
  errors <- character()
  evaluate <- function(theta) {
    if (!identical(theta, last_theta)) {
      last_value <<- tryCatch({
        value <- compute(theta)
        if (!is.finite(value$value) || any(!is.finite(value$gradient)))
          .stop("Nonfinite objective or score.")
        value
      }, error = function(e) {
        invalid <<- invalid + 1L
        errors <<- unique(c(errors, conditionMessage(e)))
        # Reject this trial point so the BFGS line search can shorten its step.
        # No emission or transition probability is replaced in a fitted model.
        list(value = 1e100, gradient = numeric(length(theta)))
      })
      last_theta <<- theta
    }
    last_value
  }
  list(fn = function(theta) evaluate(theta)$value,
       gr = function(theta) evaluate(theta)$gradient,
       diagnostics = function() list(invalid_evaluations = invalid, reasons = errors))
}
