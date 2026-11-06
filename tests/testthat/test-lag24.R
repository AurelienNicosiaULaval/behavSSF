test_that("the 24-hour score, stability transform and segment boundaries agree", {
  set.seed(2401)
  sizes <- c(28L, 24L, 29L)
  hour <- (sequence(sizes) + rep(c(20L, 5L, 23L), sizes) - 1L) %% 24L
  steps <- prepare_ssf_steps(data.frame(ID = rep(1:3, sizes),
    step = rgamma(sum(sizes), 2), angle = runif(sum(sizes), -pi, pi)), "area")
  formula <- ~ step + log_step + cos_angle + lag_log_step + lag24_log_step
  design <- .design_exact(steps, formula)
  expect_identical(design$lag_rows[["lag24_log_step"]], c(25:28, 77:81))
  i <- design$lag_rows[["lag24_log_step"]]
  expect_equal(unname(design$X[i, "lag24_log_step"]), log(steps$step[i-24L]))
  expect_error(fit_ssf(steps, formula, normalization = "exact"), "require transition_hour")
  expect_error(fit_ssf(steps[29:52, ], formula, normalization = "exact",
    transition_hour = hour[29:52]), "largest lag")

  basis <- .diel_basis(1L)
  step_basis <- basis[, -1L]
  design$diel <- list(hour = as.integer(hour), basis = basis, step_basis = step_basis)
  b <- cbind(c(-3, -.8, .4, .18, .23), c(-2, .4, 1.3, .22, -.12))
  rownames(b) <- colnames(design$X)
  gamma <- rbind(c(-2, -1), c(.3, -.4), c(.25, .2))
  clock <- rbind(c(.3, -.4), c(-.15, .23))
  theta <- .exact_transform(c(b, gamma, clock), design, 2L, inverse = TRUE)
  expect_equal(.exact_transform(theta, design, 2L), c(b, gamma, clock), tolerance = 1e-14)
  native <- .native_exact_objective(design, 2L)
  reference <- function(par) {
    z <- .exact_transform(par, design, 2L)
    movement <- matrix(z[1:10], 5L, 2L, dimnames = dimnames(b))
    schedule <- .diel_schedule(basis, matrix(z[11:16], 3L))
    offsets <- step_basis[hour + 1L, ] %*% matrix(z[17:20], 2L)
    e <- .exact_emissions(design, movement, offsets)
    .forward(e, schedule$transition, steps$ID, schedule$initial, hour = hour + 1L)
  }
  expect_equal(native$fn(theta), -reference(theta), tolerance = 1e-11)
  difference <- function(fn, j) {
    offset <- replace(numeric(length(theta)), j, 1e-5)
    (fn(theta + offset) - fn(theta - offset))/2e-5
  }
  finite_diff <- vapply(seq_along(theta), function(j) -difference(reference, j), numeric(1))
  expect_equal(native$gr(theta), finite_diff, tolerance = 1e-7)
  transform <- function(par) .exact_transform(par, design, 2L)
  J <- vapply(seq_along(theta), function(j) difference(transform, j), numeric(length(theta)))
  expect_equal(.exact_jacobian(transform(theta), design, 2L), J, tolerance = 1e-8)

  rownames(gamma) <- colnames(basis)
  rownames(clock) <- colnames(step_basis)
  fit <- fit_ssf(steps, formula, normalization = "exact", n_starts = 1L,
    transition_hour = hour, transition_harmonics = 1, step_harmonics = 1,
    short_iter = 2L, long_iter = 2L,
    initial_values = list(coefficients = b, transition_coefficients = gamma,
      step_hour_coefficients = clock))
  expect_equal(attr(logLik(fit), "df"), 20L)
  expect_equal(decode_ssf_diel(steps, fit, hour)$log_likelihood, fit$log_likelihood)
})

test_that("forward simulation uses its own 24-hour history and resets at gaps", {
  segments <- data.frame(ID = c("a", "b"), x = 0, y = 0, bearing = 0,
    n_steps = c(28L, 26L), hour = c(23L, 12L))
  basis <- .diel_basis(1L)
  gamma <- rbind(c(-2, -1), c(0, 0), c(0, 0))
  schedule <- .diel_schedule(basis, gamma)
  b <- rbind(step = c(-3, -1), log_step = c(0, 0), cos_angle = c(.5, 1),
    lag_log_step = c(.22, .19))
  fit <- structure(c(list(n_states = 2L, normalization = "exact", density_measure = "area",
    coefficients = b, transition_coefficients = gamma,
    step_hour_coefficients = rbind(c(.15, -.12), c(.25, .04))), schedule), class = "ssf_fit")
  states <- rep(c(1L, 2L, 2L), length.out = sum(segments$n_steps))
  original <- simulate_ssf(fit, segments, state_path = states, seed = 1)$tracks[[1L]]$step
  rho24 <- c(.23, -.13)
  fit$coefficients <- rbind(b, lag24_log_step = rho24)
  generated <- simulate_ssf(fit, segments, state_path = states, seed = 1)$tracks[[1L]]$step
  mu <- movement_metrics(b, "area")$mean_step
  age <- sequence(segments$n_steps)
  expected <- original
  for (i in which(age > 24L)) {
    j <- states[i]
    one <- b["lag_log_step", j]*log(expected[i-1L]/original[i-1L])
    day <- rho24[j]*log(expected[i-24L]/mu[j])
    expected[i] <- original[i]*exp(one + day)
  }
  expect_identical(generated[age <= 24L], original[age <= 24L])
  expect_equal(generated, expected, tolerance = 1e-13)
  fit$coefficients["lag24_log_step", 2] <- .95
  expect_error(simulate_ssf(fit, segments, seed = 1), "sum to less than one")
})
