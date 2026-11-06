test_that("hourly step means have exact scores with and without autoregression", {
  set.seed(814)
  hour <- as.integer(c(22:23, 0:5, 13:19))
  steps <- prepare_ssf_steps(data.frame(ID = rep(1:2, c(8, 7)),
    step = rgamma(length(hour), 2), angle = runif(length(hour), -pi, pi)), "area")
  basis <- .diel_basis(2L)
  step_basis <- basis[, -1L]
  gamma <- cbind(c(-2, .3, -.4, .35, -.2), c(-1.3, -.45, .25, -.15, .15))
  clock <- cbind(c(.3, -.15, .11, -.08), c(-.4, .16, .12, .23))
  for (ar in c(FALSE, TRUE)) {
    formula <- if (ar) ~ step + log_step + cos_angle + lag_log_step else
      ~ step + log_step + cos_angle
    design <- .design_exact(steps, formula)
    design$diel <- list(hour = hour, basis = basis, step_basis = step_basis)
    b <- cbind(c(-3, -.8, .4), c(-2, .4, 1.3))
    if (ar) b <- rbind(b, c(.18, .22))
    rownames(b) <- colnames(design$X)
    n_movement <- length(b)
    n_transition <- length(gamma)
    transition_index <- n_movement + seq_len(n_transition)
    theta <- .exact_transform(c(b, gamma, clock), design, 2L, inverse = TRUE)
    native <- .native_exact_objective(design, 2L)
    reference <- function(par) {
      natural <- .exact_transform(par, design, 2L)
      movement <- matrix(natural[seq_len(n_movement)], nrow(b), 2L, dimnames = dimnames(b))
      schedule <- .diel_schedule(basis, matrix(natural[transition_index], 5L))
      offsets <- step_basis[hour + 1L, ] %*% matrix(tail(natural, length(clock)), 4L)
      e <- .exact_emissions(design, movement, offsets)
      .forward(e, schedule$transition, steps$ID, schedule$initial, hour = hour + 1L)
    }
    expect_equal(native$fn(theta), -reference(theta), tolerance = 1e-11)
    finite_diff <- vapply(seq_along(theta), function(j) {
      offset <- replace(numeric(length(theta)), j, 1e-5)
      -(reference(theta + offset) - reference(theta - offset))/2e-5
    }, numeric(1))
    expect_equal(native$gr(theta), finite_diff, tolerance = 1e-7)

    design$diel$step_basis <- NULL
    baseline <- .native_exact_objective(design, 2L)
    base_theta <- head(theta, n_movement + n_transition)
    expect_equal(native$fn(c(base_theta, rep(0, length(clock)))), baseline$fn(base_theta),
      tolerance = 1e-11)
  }
  expect_error(fit_ssf(steps, normalization = "exact", step_harmonics = 1),
    "require transition_hour")
  rownames(gamma) <- colnames(basis)
  warm <- fit_ssf(steps, formula, normalization = "exact", n_starts = 1L,
    transition_hour = hour, transition_harmonics = 2L, step_harmonics = 2L,
    short_iter = 1L, long_iter = 1L,
    initial_values = list(coefficients = b, transition_coefficients = gamma))
  expect_equal(dim(warm$step_hour_coefficients), c(4L, 2L))
})

test_that("hourly step offsets apply at segment starts and respect the state", {
  steps <- data.frame(ID = c("a", "b"), x = 0, y = 0, bearing = 0,
    n_steps = c(3L, 2L), hour = c(23L, 12L))
  basis <- .diel_basis(1L)
  schedule <- .diel_schedule(basis, rbind(c(-2, -1), c(0, 0), c(0, 0)))
  fit <- structure(c(list(n_states = 2L, normalization = "exact", density_measure = "area",
    coefficients = rbind(step = c(-3, -1), log_step = c(0, 0), cos_angle = c(.5, 1)),
    transition_coefficients = matrix(0, 3L, 2L)), schedule), class = "ssf_fit")
  states <- c(1L, 2L, 2L, 1L, 2L)
  before <- simulate_ssf(fit, steps, state_path = states, seed = 1)$tracks[[1L]]$step
  fit$step_hour_coefficients <- rbind(c(.3, -.2), c(-.15, .25))
  after <- simulate_ssf(fit, steps, state_path = states, seed = 1)$tracks[[1L]]$step
  hour <- c(23L, 0L, 1L, 12L, 13L)
  offsets <- basis[, -1L] %*% fit$step_hour_coefficients
  expect_equal(after/before, exp(offsets[cbind(hour + 1L, states)]), tolerance = 1e-12)
})
