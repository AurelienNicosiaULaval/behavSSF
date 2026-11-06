#' Simulate exact gamma-von Mises trajectories on independent segments
#'
#' Supports the movement terms step, log_step, cos_angle and an optional single
#' cos-direction-to-target term. No finite control-point approximation is used
#' for simulation. With a target effect, the nearest supplied fixed target is
#' recomputed at the current position before each step. Movement is unbounded;
#' neither the raster extent nor observed positions are imposed as barriers.
#' With `lag_log_step`, the first step has no AR adjustment. Later
#' steps use the previous simulated length within the same segment.
#' Daily lags from 22 to 26 hours use earlier simulated lengths within
#' the same segment.
#' Hourly effects on the gamma mean apply from the first step.
#'
#' @param fit An `ssf_fit` with an explicit density measure and valid kernels.
#' @param segments Data frame with unique ID, x, y, bearing and n_steps columns.
#'   These give starting positions, previous headings, and step counts.
#'   With daily transitions, include `hour` (0 to 23) for each segment's
#'   first step origin. The simulator advances one hour per step.
#' @param n_sims Number of independent simulated replicates.
#' @param target_points Optional finite two-column matrix of fixed target points.
#'   Required if a target coefficient is nonzero. Nearest-target queries use
#'   an exact internal two-dimensional tree. At zero target distance its
#'   directional effect is set to zero.
#' @param state_path Optional fixed integer state path for all segments. NULL
#'   simulates a new Markov path, restarting from stationarity at each segment.
#'   A fixed path gives a conditional diagnostic, not a full generative check.
#' @param target_term Name of the single directional target coefficient.
#' @param seed Optional seed; the previous RNG state is restored on exit.
#' @return An `ssf_simulations` list with tracks, segment metadata, method and
#'   seed. Track rows are generated endpoints, not the initial positions.
#'   Parameters remain fixed; parameter uncertainty is not propagated.
#' @export
simulate_ssf <- function(fit, segments, n_sims = 1L, target_points = NULL,
                          state_path = NULL, target_term = "cos_diff_theta_to_target",
                          seed = NULL) {
  if (!inherits(fit, "ssf_fit")) .stop("fit must be an ssf_fit object.")
  measure <- fit$density_measure
  if (length(measure) != 1L || !measure %in% c("area", "step_angle"))
    .stop("Simulation requires an explicit fitted density measure.")
  b <- fit$coefficients
  if (any(target_term %in% c(.step_ar_term, .step_day_terms))) .stop("The AR coefficient cannot be a target term.")
  allowed <- c("step", "log_step", "cos_angle", target_term, .step_ar_term, .step_day_terms)
  if (!all(rownames(b) %in% allowed))
    .stop("Exact simulation does not support additional endpoint covariates.")
  m <- movement_metrics(b, measure)
  if (!all(m$valid)) .stop("All simulated gamma kernels must have positive shape and rate.")
  bias <- if (target_term %in% rownames(b)) b[target_term, ] else rep(0, ncol(b))
  rho <- rep(0, ncol(b))
  if (.step_ar_term %in% rownames(b)) {
    if (!identical(fit$normalization, "exact")) .stop("AR simulation requires an exact fit.")
    rho <- b[.step_ar_term, ]
    if (any(abs(rho) >= 1)) .stop("AR coefficients must lie strictly between -1 and 1.")
  }
  if (any(abs(m$persistence)+abs(bias) > 1e4))
    .stop("The implemented rejection sampler requires angular concentration at most 10000.")
  required <- c("ID", "x", "y", "bearing", "n_steps")
  if (!is.data.frame(segments) || !nrow(segments) || !all(required %in% names(segments)) ||
      anyNA(segments$ID) || anyDuplicated(segments$ID) ||
      any(!vapply(segments[required[-1]], function(x) is.numeric(x) && all(is.finite(x)), logical(1))) ||
      any(segments$n_steps < 1 | segments$n_steps != floor(segments$n_steps)))
    .stop("segments must contain unique IDs, finite starting values and positive integer lengths.")
  .scalar(n_sims, "n_sims", 1, 100000, TRUE)
  total <- sum(segments$n_steps)
  if (total*n_sims > 5e7) .stop("The requested simulation exceeds 50 million generated rows.")
  if (is.null(target_points)) target_points <- matrix(numeric(), 0L, 2L)
  target_points <- as.matrix(target_points)
  if (!is.numeric(target_points) || ncol(target_points) != 2L || any(!is.finite(target_points)) ||
      (any(bias != 0) && !nrow(target_points)))
    .stop("Nonzero target effects require finite target coordinates.")
  storage.mode(target_points) <- "double"
  if (is.null(state_path)) state_path <- integer()
  if (length(state_path) && (length(state_path) != total || anyNA(state_path) ||
      any(!state_path %in% seq_len(fit$n_states)) || any(state_path != floor(state_path))))
    .stop("state_path must contain a valid state for every generated step.")
  diel <- !is.null(fit$transition_coefficients)
  day_terms <- intersect(.step_day_terms, rownames(b))
  rho_day <- matrix(0, ncol(b), length(.step_day_terms), dimnames = list(NULL, .step_day_terms))
  if (length(day_terms)) {
    if (!diel || !.step_ar_term %in% rownames(b))
      .stop("Daily lags require an hourly AR(1) fit.")
    rho_day[, day_terms] <- t(b[day_terms, , drop = FALSE])
    if (any(abs(rho) + rowSums(abs(rho_day)) >= 1))
      .stop("The absolute AR coefficients must sum to less than one in each state.")
  }
  initial <- if (diel) fit$initial else .stationary(fit$transition)
  if (!diel) .check_hmm(matrix(0, 1, fit$n_states), fit$transition, 1, initial)
  parameters <- as.matrix(data.frame(shape = m$shape, rate = m$rate,
                                      persistence = m$persistence, target = bias, rho = rho))
  parameters <- cbind(parameters, rho_day)
  step_offset <- .diel_step_offsets(fit$step_hour_coefficients)
  if (!is.null(step_offset)) parameters <- cbind(parameters, t(step_offset))
  starts <- as.matrix(segments[c("x", "y", "bearing", "n_steps")])
  if (diel) starts <- cbind(starts, hour = .diel_hour(segments$hour, nrow(segments)))
  storage.mode(starts) <- "double"
  .with_seed(seed, function() {
    matrices <- .Call(C_behavSSF_simulate, parameters, fit$transition, initial,
                       starts, target_points, as.integer(n_sims), as.integer(state_path))
    tracks <- lapply(matrices, function(x) data.frame(
      ID = rep(segments$ID, segments$n_steps), index = sequence(segments$n_steps),
      x = x[,1], y = x[,2], step = x[,3], angle = x[,4], state = as.integer(x[,5])))
    structure(list(tracks = tracks, segments = segments,
      method = if (length(state_path)) "fixed_state_path" else "markov",
      seed = seed, parameter_uncertainty = FALSE), class = "ssf_simulations")
  })
}
