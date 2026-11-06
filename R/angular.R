#' Combine directional persistence and a target-direction effect
#'
#' For turning angle phi and target bearing psi relative to the previous
#' heading, the angular predictor is a*cos(phi) + b*cos(phi-psi).
#' Its conditional von Mises resultant is (a+b*cos(psi), b*sin(psi)). Thus
#' abs(a) alone is the intrinsic persistence concentration, not the full
#' concentration when a target effect is present.
#'
#' @param persistence Signed coefficient a for cos(turning angle).
#' @param target_bias Signed coefficient b for direction toward the target.
#' @param relative_bearing Target bearing minus previous heading, in radians.
#' @return A data frame with `mean_angle` and nonnegative `concentration`.
#'   The mean is NA for a uniform angular distribution. Arguments must have
#'   a common length or length one.
#' @export
#' @examples
#' angular_parameters(2, 1, pi / 2)
angular_parameters <- function(persistence, target_bias = 0, relative_bearing = 0) {
  values <- list(persistence, target_bias, relative_bearing)
  n <- max(lengths(values))
  if (!n || any(!lengths(values) %in% c(1L, n)) ||
      any(!vapply(values, function(x) is.numeric(x) && all(is.finite(x)), logical(1))))
    .stop("Angular parameters must be finite numeric vectors of compatible lengths.")
  a <- rep(persistence, length.out = n)
  b <- rep(target_bias, length.out = n)
  psi <- rep(relative_bearing, length.out = n)
  x <- a + b * cos(psi)
  y <- b * sin(psi)
  concentration <- sqrt(x^2 + y^2)
  mean <- atan2(y, x)
  mean[concentration == 0] <- NA_real_
  data.frame(mean_angle = mean, concentration = concentration)
}
