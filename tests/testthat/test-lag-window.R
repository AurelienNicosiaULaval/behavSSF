test_that("joint daily lags have the right scores and stay inside the stability bound", {
  set.seed(22026)
  sizes <- c(29L, 26L, 28L)
  hour <- as.integer((sequence(sizes) + rep(c(20L, 5L, 23L), sizes) - 1L) %% 24L)
  steps <- prepare_ssf_steps(data.frame(ID = rep(1:3, sizes),
    step = rgamma(sum(sizes), 2), angle = runif(sum(sizes), -pi, pi)), "area")
  formula <- ~ step + log_step + cos_angle + lag_log_step + lag24_log_step +
    lag26_log_step + lag22_log_step + lag23_log_step + lag25_log_step
  design <- .design_exact(steps, formula)
  expect_identical(design$lag_rows[["lag26_log_step"]], c(27:29, 82:83))
  for (lag in 22:26) {
    term <- paste0("lag", lag, "_log_step")
    i <- design$lag_rows[[term]]
    expect_identical(i, which(sequence(sizes) > lag))
    expect_equal(unname(design$X[i, term]), log(steps$step[i-lag]))
  }

  basis <- .diel_basis(1L)
  step_basis <- basis[, -1L]
  design$diel <- list(hour = hour, basis = basis, step_basis = step_basis)
  b <- cbind(c(-3, -.8, .4, .18, .12, .04, -.06, .07, .02),
    c(-2, .4, 1.3, .22, -.06, .07, .12, .04, .05))
  rownames(b) <- colnames(design$X)
  gamma <- rbind(c(-2, -1), c(.3, -.4), c(.25, .2))
  clock <- rbind(c(.3, -.4), c(-.15, .23))
  natural <- c(b, gamma, clock)
  theta <- .exact_transform(natural, design, 2L, inverse = TRUE)
  transform <- function(par) .exact_transform(par, design, 2L)
  expect_equal(transform(theta), natural, tolerance = 1e-14)
  native <- .native_exact_objective(design, 2L)
  reference <- function(par) {
    z <- transform(par)
    movement <- matrix(z[1:18], 9L, 2L, dimnames = dimnames(b))
    schedule <- .diel_schedule(basis, matrix(z[19:24], 3L))
    offsets <- step_basis[hour + 1L, ] %*% matrix(z[25:28], 2L)
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
  J <- vapply(seq_along(theta), function(j) difference(transform, j), numeric(length(theta)))
  expect_equal(.exact_jacobian(natural, design, 2L), J, tolerance = 1e-8)
  too_large <- replace(natural, 4L, .9)
  expect_error(.exact_transform(too_large, design, 2L, inverse = TRUE), "sum to less than one")
})

test_that("simulation respects each daily lag and segment boundary", {
  segments <- data.frame(ID = c("a", "b"), x = 0, y = 0, bearing = 0,
    n_steps = c(28L, 27L), hour = c(23L, 12L))
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
  day <- rbind(lag22_log_step = c(.1, -.12), lag24_log_step = c(-.04, .15),
    lag26_log_step = c(.15, .06))
  fit$coefficients <- rbind(b, day)
  generated <- simulate_ssf(fit, segments, state_path = states, seed = 1)$tracks[[1L]]$step
  mu <- movement_metrics(b, "area")$mean_step
  age <- sequence(segments$n_steps)
  expected <- original
  for (i in which(age > 22L)) {
    j <- states[i]
    one <- b["lag_log_step", j]*log(expected[i-1L]/original[i-1L])
    lag <- c(22L, 24L, 26L)
    available <- lag < age[i]
    day_shift <- sum(day[available, j]*log(expected[i-lag[available]]/mu[j]))
    expected[i] <- original[i]*exp(one + day_shift)
  }
  expect_identical(generated[age <= 22L], original[age <= 22L])
  expect_equal(generated, expected, tolerance = 1e-13)
})
