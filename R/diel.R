.diel_hour <- function(hour, n, segment = NULL) {
  if (!is.numeric(hour) || length(hour) != n ||
      any(!is.finite(hour) | hour < 0 | hour > 23 | hour != floor(hour)))
    .stop("Supply one integer step-origin hour between 0 and 23 per row.")
  if (!is.null(segment)) {
    within <- segment[-1L] == segment[-length(segment)]
    if (any(diff(hour)[within] %% 24 != 1L))
      .stop("Step-origin hours must advance by one within each segment.")
  }
  as.integer(hour)
}

.diel_basis <- function(harmonics) {
  hour <- 0:23
  wave <- outer(hour, seq_len(harmonics), function(h, j) h*j*2*pi/24)
  basis <- cbind(1, cos(wave), sin(wave))
  colnames(basis) <- c("Intercept", paste0("cos", seq_len(harmonics)),
    paste0("sin", seq_len(harmonics)))
  basis
}

.diel_schedule <- function(basis, coefficients) {
  eta <- basis %*% coefficients
  states <- paste0("S", 1:2)
  transition <- array(0, c(2L, 2L, 24L), list(states, states, as.character(0:23)))
  transition[1, 2, ] <- stats::plogis(eta[, 1L])
  transition[2, 1, ] <- stats::plogis(eta[, 2L])
  transition[1, 1, ] <- stats::plogis(-eta[, 1L])
  transition[2, 2, ] <- stats::plogis(-eta[, 2L])
  day <- Reduce(function(product, h) product %*% transition[, , h], seq_len(24L), init = diag(2))
  initial <- matrix(0, 2L, 24L, dimnames = list(states, as.character(0:23)))
  delta <- .stationary(day)
  for (h in seq_len(24L)) {
    delta <- drop(delta %*% transition[, , h])
    initial[, h] <- delta
  }
  list(transition = transition, initial = initial)
}

.diel_step_offsets <- function(coefficients) {
  if (is.null(coefficients)) return(NULL)
  .diel_basis(nrow(coefficients)/2)[, -1L, drop = FALSE] %*% coefficients
}

#' Score and decode a fitted HMM-SSF with daily transitions
#'
#' Applies a fitted hourly transition schedule to new observations. At the
#' start of each segment, the hidden state uses the corresponding hour of
#' the periodic stationary distribution. Step-length autoregression starts
#' afresh at each segment boundary. All daily lag terms use only
#' observations in that segment. Fitted hourly effects on mean step
#' length apply from the first step of each segment.
#'
#' @param data Observed hourly `ssf_steps`.
#' @param fit A fit from [fit_ssf()] with daily transitions.
#' @param transition_hour Integer hour of the step origin, from 0 to 23.
#'   Supply one value per row in a consistent time zone.
#' @return The log likelihood, Viterbi states, and smoothed state probabilities.
#' @export
decode_ssf_diel <- function(data, fit, transition_hour) {
  if (!inherits(fit, "ssf_fit") || is.null(fit$transition_coefficients))
    .stop("fit must be an HMM-SSF with daily transitions.")
  if (!identical(attr(data, "density_measure"), fit$density_measure))
    .stop("The fitted and new data must use the same density measure.")
  d <- .design_exact(data, fit$formula)
  hour <- .diel_hour(transition_hour, nrow(d$X), d$segment)
  offset <- .diel_step_offsets(fit$step_hour_coefficients)
  if (!is.null(offset)) offset <- offset[hour + 1L, , drop = FALSE]
  emissions <- .exact_emissions(d, fit$coefficients, offset)
  .forward(emissions, fit$transition, d$segment, fit$initial, decode = TRUE,
           hour = hour + 1L)
}
