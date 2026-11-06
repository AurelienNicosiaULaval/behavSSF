simulation_fixture <- function(target = FALSE) {
  b <- rbind(step = -3, log_step = 0, cos_angle = 0)
  if (target) b <- rbind(b, cos_diff_theta_to_target = 5)
  structure(list(coefficients = b, transition = matrix(1,1,1), n_states = 1L,
    density_measure = "area"), class = "ssf_fit")
}

test_that("exact forward simulation has gamma and uniform circular moments", {
  fit <- simulation_fixture()
  segments <- data.frame(ID = 1L, x = 0, y = 0, bearing = 0, n_steps = 10000)
  set.seed(202); previous <- .Random.seed
  x <- simulate_ssf(fit, segments, seed = 203)
  expect_identical(.Random.seed, previous)
  d <- x$tracks[[1]]
  expect_lt(abs(mean(d$step)-2/3), 6*sqrt(2)/3/sqrt(10000))
  expect_lt(abs(mean(cos(d$angle))), 6*sqrt(.5/10000))
  expect_true(all(d$state == 1))
  expect_identical(x, simulate_ssf(fit, segments, seed = 203))
})

test_that("target direction is queried at the current origin", {
  fit <- simulation_fixture(TRUE)
  segments <- data.frame(ID = 1L, x = 0, y = 0, bearing = 0, n_steps = 1L)
  targets <- rbind(c(1,0), c(-100,0), c(10,10))
  x <- simulate_ssf(fit, segments, n_sims = 10000, target_points = targets, seed = 204)
  angles <- vapply(x$tracks, function(d) d$angle[1], numeric(1))
  expected <- besselI(5,1)/besselI(5,0)
  variance <- (1+besselI(5,2)/besselI(5,0))/2-expected^2
  expect_lt(abs(mean(cos(angles))-expected), 6*sqrt(variance/10000))
  expect_lt(abs(mean(sin(angles))), .03)
  expect_error(simulate_ssf(fit, segments), "target coordinates")
})

test_that("segments and observable displacement calculations preserve boundaries", {
  fit <- simulation_fixture()
  segments <- data.frame(ID = c("a", "b"), x = c(0,100), y = c(0,200),
    bearing = c(0,1), n_steps = c(3,4))
  x <- simulate_ssf(fit, segments, seed = 205)$tracks[[1]]
  for (i in 1:2) {
    d <- x[x$ID == segments$ID[i], ]
    expect_equal(d$x[1], segments$x[i] + d$step[1]*cos(segments$bearing[i]+d$angle[1]))
    expect_equal(d$y[1], segments$y[i] + d$step[1]*sin(segments$bearing[i]+d$angle[1]))
  }
  diagnostics <- movement_diagnostics(x, segments, lags = c(1,2))
  expect_equal(unname(diagnostics["msd_lag_1"]), mean(x$step^2), tolerance = 1e-10)
  shifted <- x; shifted$x <- shifted$x+500; shifted$y <- shifted$y-800
  shifted_segments <- segments; shifted_segments$x <- shifted_segments$x+500
  shifted_segments$y <- shifted_segments$y-800
  expect_equal(movement_diagnostics(shifted, shifted_segments, lags = c(1,2)), diagnostics,
                tolerance = 1e-10)
  invalid <- x; invalid$x[2] <- invalid$x[2]+5
  expect_error(movement_diagnostics(invalid, segments), "do not match step lengths")
  conditional <- simulate_ssf(fit, segments, state_path = rep(1L,7), seed = 206)
  expect_identical(conditional$method, "fixed_state_path")
})

test_that("simulated Markov transitions use rows of the transition matrix", {
  fit <- simulation_fixture()
  fit$n_states <- 2L
  fit$coefficients <- cbind(fit$coefficients, fit$coefficients)
  fit$transition <- matrix(c(.8,.2,.4,.6), 2, byrow = TRUE)
  segments <- data.frame(ID = 1L, x = 0, y = 0, bearing = 0, n_steps = 30000L)
  states <- simulate_ssf(fit, segments, seed = 207)$tracks[[1]]$state
  empirical <- prop.table(table(states[-length(states)], states[-1]), 1)
  expect_lt(max(abs(unname(empirical)-fit$transition)), .03)
  expect_lt(abs(mean(states == 1)-2/3), .03)
})
