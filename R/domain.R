# Moments of the canonical radial kernel on a finite availability disk.
# For a negative linear coefficient, incomplete gamma ratios avoid quadrature
# at the origin. For a nonnegative coefficient, r = R * u^(1/shape) removes
# the integrable singularity. The constant exponential shift cancels.
.bounded_radial_moments <- function(shape, linear, radius) {
  mean <- sd <- rep(NA_real_, length(shape))
  valid <- is.finite(shape) & shape > 0 & is.finite(linear)
  for (i in which(valid)) {
    a <- shape[i]
    beta <- linear[i]
    moments <- tryCatch({
      if (beta < 0) {
        log_mass <- stats::pgamma(radius, a, rate = -beta, log.p = TRUE)
        log_first <- log(a / -beta) +
          stats::pgamma(radius, a + 1, rate = -beta, log.p = TRUE) - log_mass
        log_second <- log(a) + log(a + 1) - 2 * log(-beta) +
          stats::pgamma(radius, a + 2, rate = -beta, log.p = TRUE) - log_mass
        c(exp(log_first), exp(log_second))
      } else if (beta == 0) {
        c(radius * a / (a + 1), radius^2 * a / (a + 2))
      } else {
        z <- beta * radius
        integral <- vapply(0:2, function(power) stats::integrate(function(u)
          u^(power / a) * exp(z * (u^(1 / a) - 1)),
          0, 1, rel.tol = 1e-9, subdivisions = 1000L)$value, numeric(1))
        radius^(1:2) * integral[2:3] / integral[1]
      }
    }, error = function(e) c(NA_real_, NA_real_))
    variance <- moments[2] - moments[1]^2
    if (any(!is.finite(moments)) || moments[1] <= 0 ||
        variance < -1e-10 * moments[2]) {
      valid[i] <- FALSE
    } else {
      mean[i] <- moments[1]
      sd[i] <- sqrt(max(0, variance))
    }
  }
  list(mean = mean, sd = sd, valid = valid)
}
