.lse <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.stationary <- function(tpm) {
  k <- nrow(tpm)
  if (k == 1L) return(1)
  a <- t(tpm) - diag(k)
  a[k, ] <- 1
  b <- c(rep(0, k - 1L), 1)
  delta <- as.numeric(solve(a, b))
  if (any(!is.finite(delta)) || any(delta < -1e-10))
    .stop("Cannot determine a stationary initial distribution.")
  delta <- pmax(delta, 0)
  delta / sum(delta)
}

.check_hmm <- function(log_emissions, transition, segment, initial) {
  e <- as.matrix(log_emissions)
  p <- as.matrix(transition)
  if (!is.numeric(e) || !nrow(e) || !ncol(e) || anyNA(e) || any(e == Inf))
    .stop("log_emissions must be a nonempty numeric matrix without NA or +Inf.")
  k <- ncol(e)
  if (!is.numeric(p) || !identical(dim(p), c(k, k)) || any(!is.finite(p)) ||
      any(p < 0) || any(abs(rowSums(p) - 1) > 1e-10))
    .stop("transition must be a nonnegative row-stochastic K by K matrix.")
  if (length(segment) != nrow(e) || anyNA(segment) || anyDuplicated(rle(as.character(segment))$values))
    .stop("segment must identify contiguous independent sequences, one ID per row.")
  if (is.null(initial)) initial <- .stationary(p)
  if (!is.numeric(initial) || length(initial) != k || any(!is.finite(initial)) ||
      any(initial < 0) || abs(sum(initial) - 1) > 1e-10)
    .stop("initial must be a probability vector of length K.")
  list(e = e, p = p, segment = segment, initial = initial)
}

.forward <- function(e, p, segment, initial, decode = FALSE, hour = NULL) {
  k <- ncol(e)
  lp <- log(p)
  total <- 0
  if (decode) {
    posterior <- matrix(NA_real_, nrow(e), k)
    path <- integer(nrow(e))
  }
  for (ind in split(seq_len(nrow(e)), factor(segment, levels = unique(segment)))) {
    n <- length(ind)
    a <- matrix(0, n, k)
    scales <- numeric(n)
    start <- if (is.null(hour)) initial else initial[, hour[ind[1L]]]
    transition_at <- function(t) if (is.null(hour)) lp else lp[, , hour[ind[t]]]
    a[1, ] <- log(start) + e[ind[1], ]
    scales[1] <- .lse(a[1, ])
    if (!is.finite(scales[1])) .stop("The sequence has zero or nonfinite likelihood.")
    a[1, ] <- a[1, ] - scales[1]
    if (n > 1L) for (t in 2:n) {
      log_p <- transition_at(t)
      a[t, ] <- vapply(seq_len(k), function(j) .lse(a[t - 1L, ] + log_p[, j]), numeric(1)) + e[ind[t], ]
      scales[t] <- .lse(a[t, ])
      if (!is.finite(scales[t])) .stop("The sequence has zero or nonfinite likelihood.")
      a[t, ] <- a[t, ] - scales[t]
    }
    total <- total + sum(scales)
    if (decode) {
      b <- matrix(0, n, k)
      if (n > 1L) for (t in (n - 1L):1L) {
        log_p <- transition_at(t + 1L)
        b[t, ] <- vapply(seq_len(k), function(i)
          .lse(log_p[i, ] + e[ind[t + 1L], ] + b[t + 1L, ]), numeric(1)) - scales[t + 1L]
      }
      ab <- a + b
      posterior[ind, ] <- exp(ab - apply(ab, 1, .lse))
      score <- matrix(0, n, k)
      previous <- matrix(0L, n, k)
      score[1, ] <- log(start) + e[ind[1], ]
      if (n > 1L) for (t in 2:n) {
        log_p <- transition_at(t)
        for (j in seq_len(k)) {
          candidates <- score[t - 1L, ] + log_p[, j]
          previous[t, j] <- which.max(candidates)
          score[t, j] <- max(candidates) + e[ind[t], j]
        }
      }
      path[ind[n]] <- which.max(score[n, ])
      if (n > 1L) for (t in (n - 1L):1L)
        path[ind[t]] <- previous[t + 1L, path[ind[t + 1L]]]
    }
  }
  if (decode) list(log_likelihood = total, state = path, posterior = posterior)
  else total
}

#' Decode independent hidden Markov sequences in log space
#'
#' Uses a homogeneous transition matrix. Restarts both filtering and Viterbi
#' decoding at every segment boundary, including segments of length one.
#' Ties in Viterbi decoding are resolved by the first state index.
#'
#' @param log_emissions Numeric matrix with one row per step and one column per state.
#' @param transition Row-stochastic transition probability matrix.
#' @param segment Independent sequence identifiers, one per step.
#' @param initial Initial state probabilities, used at each boundary. If NULL,
#'   the stationary distribution is used; it must be uniquely computable.
#' @return A list with `log_likelihood`, global Viterbi `state`, and smoothed
#'   marginal state `posterior`. These probabilities condition on fixed parameters.
#' @export
#' @examples
#' decode_states(log(matrix(c(.8, .2, .3, .7), 2, byrow = TRUE)),
#'               matrix(c(.9, .1, .1, .9), 2, byrow = TRUE))
decode_states <- function(log_emissions, transition,
                          segment = rep(1L, nrow(as.matrix(log_emissions))), initial = NULL) {
  x <- .check_hmm(log_emissions, transition, segment, initial)
  .forward(x$e, x$p, x$segment, x$initial, decode = TRUE)
}

.emissions <- function(design, coefficients) {
  .importance_emissions(design, design$X %*% coefficients)
}

.importance_emissions <- function(design, lp) {
  out <- matrix(0, length(design$groups), ncol(lp))
  for (s in seq_along(design$groups)) {
    i <- design$groups[[s]]
    denominator <- apply(lp[i, , drop = FALSE] - design$logw[i], 2, .lse)
    # Published importance approximation: numerator is unweighted; the
    # denominator includes the case and J controls and is divided by J.
    out[s, ] <- lp[i[1], ] - denominator + log(length(i) - 1L)
  }
  out
}

#' Evaluate importance-sampling or exact SSF emission densities
#'
#' For J controls and the case in each choice set, returns
#' eta_case - log(sum(exp(eta) / w)) + log(J). The case is included in the
#' denominator and its numerator is unweighted, matching the published code.
#' This is an approximate density, not a probability over the finite choice set.
#' The measure underlying w must agree with the kernel being fitted.
#' With exact normalisation, controls and weights are absent and the returned
#' emissions are properly normalised canonical-kernel densities.
#'
#' @param data Case-control [prepare_ssf_data()] data, or observed
#'   [prepare_ssf_steps()] data for exact normalisation.
#' @param formula One-sided SSF formula without an identifiable intercept.
#' @param coefficients Numeric matrix, with design columns in rows and states
#'   in columns. Row names must exactly match the design column names.
#' @param normalization importance uses the published approximation; exact
#'   uses observed ssf_steps and the analytic canonical-kernel normaliser.
#' @param step_family Optional length family, as in [fit_ssf()]. Coefficients
#'   use the named natural parameters returned by that fit.
#' @return A matrix of log emission densities, one row per stratum.
#' @export
ssf_log_emissions <- function(data, formula, coefficients,
                               normalization = c("importance", "exact"), step_family = NULL) {
  normalization <- match.arg(normalization)
  family <- .step_family(step_family)
  exact <- normalization == "exact"
  d <- if (!is.null(family)) .family_design(data, formula, family, exact) else
    if (exact) .design_exact(data, formula) else .design(data, formula)
  b <- as.matrix(coefficients)
  if (!is.numeric(b) || nrow(b) != ncol(d$X) || !ncol(b) || any(!is.finite(b)) ||
      !identical(rownames(b), colnames(d$X)))
    .stop("coefficients must be finite and have row names matching the design columns in order.")
  if (!is.null(family)) return(.family_emissions(d, b, exact))
  if (exact) .exact_emissions(d, b) else .emissions(d, b)
}
