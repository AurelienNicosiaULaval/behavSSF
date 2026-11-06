.rvm_exact <- function(n, mean, concentration) {
  if (concentration == 0) return(stats::runif(n, -pi, pi))
  out <- numeric(n)
  filled <- 0L
  acceptance <- besselI(concentration, 0, expon.scaled = TRUE)
  while (filled < n) {
    count <- min(5e6, max(32L, ceiling((n-filled) / acceptance * 1.05)))
    angle <- stats::runif(count, -pi, pi)
    accepted <- angle[log(stats::runif(count)) <= concentration * (cos(angle)-1)]
    take <- min(length(accepted), n-filled)
    if (take > 0L) out[filled + seq_len(take)] <- accepted[seq_len(take)]
    filled <- filled + take
  }
  atan2(sin(out + mean), cos(out + mean))
}

#' Generate gamma-von Mises control steps with explicit area densities
#'
#' All observed rows and their segment identifiers are retained. The target
#' bearing is fixed at the step origin and shared by the case and controls.
#' Controls are not filtered by a raster or spatial observation window.
#' This implementation is for an untruncated gamma step proposal and a von
#' Mises angular proposal. A pool can be nested by retaining alternative <= J.
#'
#' @param observed One row per observed step, with ID, stratum, step, angle,
#'   previous_bearing and target_bearing. Bearings and angles are in radians.
#' @param n_controls Positive integer number of controls per observed step.
#' @param proposal Named numeric vector with gamma `shape`, gamma `rate`,
#'   angular `mean_angle` and angular `concentration`.
#' @param seed Optional integer RNG seed, restored on exit.
#' @return Validated `ssf_data` with area-density weights and alternative zero
#'   for the case, 1:n_controls for controls. Source row order is preserved.
#' @export
generate_ssf_controls <- function(observed, n_controls, proposal, seed = NULL) {
  .scalar(n_controls, "n_controls", 1, 100000, TRUE)
  required <- c("ID", "stratum", "step", "angle", "previous_bearing", "target_bearing")
  if (!is.data.frame(observed) || !nrow(observed) || !all(required %in% names(observed)))
    .stop("observed must contain one row per step and the required bearings.")
  numeric <- c("step", "angle", "previous_bearing", "target_bearing")
  if (any(!vapply(observed[numeric], function(x) is.numeric(x) && all(is.finite(x)), logical(1))) ||
      any(observed$step <= 0)) .stop("Observed steps must be positive and bearings finite.")
  parameters <- c("shape", "rate", "mean_angle", "concentration")
  if (!is.numeric(proposal) || is.null(names(proposal)) ||
      !all(parameters %in% names(proposal)) || anyDuplicated(names(proposal)) ||
      any(!is.finite(proposal[parameters])) || proposal["shape"] <= 0 ||
      proposal["rate"] <= 0 || proposal["concentration"] < 0 ||
      proposal["concentration"] > 1e4)
    .stop("proposal must specify positive gamma parameters and finite angular parameters.")
  .with_seed(seed, function() {
    n <- nrow(observed)
    j <- as.integer(n_controls)
    steps <- matrix(stats::rgamma(n*j, shape = proposal["shape"], rate = proposal["rate"]), nrow = n)
    angles <- matrix(.rvm_exact(n*j, proposal["mean_angle"], proposal["concentration"]), nrow = n)
    step <- as.vector(t(cbind(observed$step, steps)))
    angle <- as.vector(t(cbind(observed$angle, angles)))
    row <- rep(seq_len(n), each = j+1L)
    log_angular <- proposal["concentration"] *
      (cos(angle-proposal["mean_angle"]) - 1) - log(2*pi) -
      log(besselI(proposal["concentration"], 0, expon.scaled = TRUE))
    log_density <- stats::dgamma(step, shape = proposal["shape"], rate = proposal["rate"], log = TRUE) -
      log(step) + log_angular
    out <- data.frame(ID = observed$ID[row], stratum = observed$stratum[row],
      obs = as.integer(rep(0:j, n) == 0L), w = exp(log_density),
      step = step, angle = angle, log_step = log(step), cos_angle = cos(angle),
      cos_diff_theta_to_target = cos(observed$previous_bearing[row] + angle -
                                       observed$target_bearing[row]),
      alternative = rep(0:j, n))
    prepare_ssf_data(out, density_measure = "area")
  })
}
