.step_families <- list(
  gamma = c(step_shape = TRUE, step_rate = TRUE),
  lognormal = c(step_meanlog = FALSE, step_sdlog = TRUE),
  inverse_gaussian = c(step_mean = TRUE, step_shape = TRUE),
  halfnormal = c(step_sd = TRUE),
  weibull = c(step_shape = TRUE, step_scale = TRUE),
  exponential = c(step_rate = TRUE)
)

.step_family <- function(family) {
  if (is.null(family)) return(NULL)
  if (!is.character(family) || length(family) != 1L || is.na(family) ||
      !family %in% names(.step_families))
    .stop("step_family must be one of: ", paste(names(.step_families), collapse = ", "), ".")
  family
}

.family_design <- function(data, formula, family, exact) {
  if (!inherits(formula, "formula") || length(formula) != 2L)
    .stop("formula must be one-sided.")
  if (!is.null(attr(data, "step_max")) && is.finite(attr(data, "step_max")))
    .stop("Finite step_max currently requires the canonical importance kernel.")
  terms <- stats::terms(formula)
  radial <- .step_families[[family]]
  if (length(attr(terms, "offset")) ||
      any(all.vars(formula) %in% c(names(radial), "step", "log_step", .step_ar_term,
                                  .step_day_terms)))
    .stop("With step_family, the formula supplies angular or selection terms only.")
  if (exact) {
    if (!inherits(data, "ssf_steps"))
      .stop("Exact fitting requires observed data from prepare_ssf_steps().")
    allowed <- .exact_terms[3:4]
    if (!all(attr(terms, "term.labels") %in% allowed))
      .stop("Exact step-family fits support only cos_angle and a prospective target term.")
    data <- prepare_ssf_steps(as.data.frame(data), attr(data, "density_measure"))
    has_target <- allowed[2L] %in% attr(terms, "term.labels")
    if (has_target && !all(c("relative_target_bearing", "target_active") %in% names(data)))
      .stop("Exact target normalisation requires a bearing at the step origin.")
    x <- stats::model.matrix(formula, data)
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    d <- list(data = data, X = x, cases = seq_len(nrow(data)), segment = data$ID,
      relative_bearing = if (has_target) data$relative_target_bearing else rep(0, nrow(data)),
      target_active = if (has_target) as.double(data$target_active) else rep(0, nrow(data)))
  } else {
    d <- .design(data, formula, allow_empty = TRUE)
  }
  measure <- attr(d$data, "density_measure")
  if (length(measure) != 1L || !measure %in% c("area", "step_angle"))
    .stop("Step families require an explicit density measure.")
  r <- d$data[["step"]]
  if (!is.numeric(r) || length(r) != nrow(d$data) || any(!is.finite(r) | r <= 0))
    .stop("Step families require finite, strictly positive step lengths.")
  d$X <- cbind(matrix(0, nrow(d$data), length(radial),
                      dimnames = list(NULL, names(radial))), d$X)
  if (anyDuplicated(colnames(d$X))) .stop("The step-family parameter names are reserved.")
  d$family <- family
  d$positive <- which(radial)
  d$log_step <- log(r)
  d$area <- measure == "area"
  d
}

.family_positive_indices <- function(d, k) {
  as.vector(outer(d$positive, ncol(d$X) * (seq_len(k) - 1L), "+"))
}

.family_transform <- function(theta, d, k, inverse = FALSE) {
  i <- .family_positive_indices(d, k)
  if (inverse && any(theta[i] <= 0)) .stop("Step-family scale and shape parameters must be positive.")
  theta[i] <- if (inverse) log(theta[i]) else exp(theta[i])
  theta
}

.family_jacobian <- function(theta, d, k) {
  i <- .family_positive_indices(d, k)
  replace(rep(1, length(theta)), i, exp(theta[i]))
}

.family_start <- function(r, family) {
  m <- mean(r)
  spread <- stats::sd(r)
  cv2 <- if (is.finite(spread) && spread > 0) max(.01, (spread / m)^2) else 1
  log_sd <- stats::sd(log(r))
  switch(family,
    gamma = c(1 / cv2, 1 / (m * cv2)),
    lognormal = c(mean(log(r)), if (is.finite(log_sd) && log_sd > 0) log_sd else 1),
    inverse_gaussian = c(m, m / cv2),
    halfnormal = sqrt(mean(r^2)),
    weibull = c(1, m),
    exponential = 1 / m)
}

.family_log_density <- function(r, u, b, family) {
  switch(family,
    gamma = stats::dgamma(r, b[1L], rate = b[2L], log = TRUE),
    lognormal = stats::dlnorm(r, b[1L], b[2L], log = TRUE),
    inverse_gaussian = (log(b[2L]) - log(2 * pi) - 3 * u -
      b[2L] * (expm1(u - log(b[1L])) / sqrt(r))^2) / 2,
    halfnormal = log(2) + stats::dnorm(r, sd = b[1L], log = TRUE),
    weibull = stats::dweibull(r, b[1L], b[2L], log = TRUE),
    exponential = stats::dexp(r, b[1L], log = TRUE))
}

.family_emissions <- function(d, b, exact) {
  radial <- .step_families[[d$family]]
  if (any(!is.finite(b)) || any(b[d$positive, , drop = FALSE] <= 0))
    .stop("Step-family parameters must be finite and within their natural domains.")
  eta <- d$X %*% b
  for (i in seq_len(ncol(b))) {
    eta[, i] <- eta[, i] +
      .family_log_density(d$data$step, d$log_step, b[seq_along(radial), i], d$family)
    if (d$area) eta[, i] <- eta[, i] - d$log_step
  }
  if (!exact) return(.importance_emissions(d, eta))
  persistence <- if ("cos_angle" %in% rownames(b)) b["cos_angle", ] else rep(0, ncol(b))
  target <- if (.exact_terms[4L] %in% rownames(b)) b[.exact_terms[4L], ] else rep(0, ncol(b))
  for (i in seq_len(ncol(b)))
    eta[, i] <- eta[, i] -
      .angular_log_z(persistence[i], target[i], d$relative_bearing, d$target_active)
  eta
}

.family_metrics <- function(b, convention, family) {
  parameters <- .step_families[[family]]
  if (!is.numeric(b) || !ncol(b) || is.null(rownames(b)) || anyDuplicated(rownames(b)) ||
      !all(names(parameters) %in% rownames(b)))
    .stop("coefficients must contain the named parameters of step_family.")
  a <- b[names(parameters)[1L], ]
  z <- if (length(parameters) == 2L) b[names(parameters)[2L], ] else NA_real_
  valid <- is.finite(a)
  if (parameters[1L]) valid <- valid & a > 0
  if (length(parameters) == 2L) valid <- valid & is.finite(z) & z > 0
  a[!valid] <- NA_real_
  if (length(parameters) == 2L) z[!valid] <- NA_real_
  moments <- switch(family,
    gamma = list(a / z, sqrt(a) / z, a, z),
    lognormal = list(exp(a + z^2 / 2),
      exp(a + z^2 / 2) * sqrt(expm1(z^2)), NA_real_, NA_real_),
    inverse_gaussian = list(a, a * sqrt(a / z), z, NA_real_),
    halfnormal = list(a * sqrt(2 / pi), a * sqrt(1 - 2 / pi), NA_real_, NA_real_),
    weibull = {
      log_mean <- lgamma(1 + 1 / a)
      mu <- z * exp(log_mean)
      list(mu, mu * sqrt(pmax(0, expm1(lgamma(1 + 2 / a) - 2 * log_mean))), a, NA_real_)
    },
    exponential = list(1 / a, 1 / a, 1, a))
  persistence <- if ("cos_angle" %in% rownames(b)) b["cos_angle", ] else rep(0, ncol(b))
  valid <- valid & is.finite(moments[[1L]]) & is.finite(moments[[2L]]) &
    moments[[1L]] > 0 & is.finite(persistence)
  data.frame(state = seq_len(ncol(b)), shape = moments[[3L]], rate = moments[[4L]],
    mean_step = ifelse(valid, moments[[1L]], NA_real_),
    sd_step = ifelse(valid, moments[[2L]], NA_real_), persistence = persistence,
    concentration = abs(persistence), valid = valid, convention = convention,
    row.names = NULL)
}
