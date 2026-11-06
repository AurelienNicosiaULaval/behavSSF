test_that("forward, smoothing and Viterbi agree with exhaustive enumeration", {
  p <- matrix(c(.8, .2, .35, .65), 2, byrow = TRUE)
  initial <- c(.4, .6)
  emission <- matrix(c(.2, .8, .9, .1, .3, .7), 3, byrow = TRUE)
  paths <- as.matrix(expand.grid(1:2, 1:2, 1:2))
  mass <- apply(paths, 1, function(s) {
    initial[s[1]] * emission[1, s[1]] *
      p[s[1], s[2]] * emission[2, s[2]] *
      p[s[2], s[3]] * emission[3, s[3]]
  })
  expected_posterior <- sapply(1:2, function(k)
    vapply(1:3, function(t) sum(mass[paths[, t] == k]) / sum(mass), numeric(1)))
  result <- decode_states(log(emission), p, initial = initial)
  expect_equal(result$log_likelihood, log(sum(mass)), tolerance = 1e-12)
  expect_equal(result$posterior, expected_posterior, tolerance = 1e-12)
  expect_equal(result$state, unname(paths[which.max(mass), ]))
})

test_that("segment boundaries reset all recursions, including singletons", {
  p <- matrix(c(.999, .001, .001, .999), 2, byrow = TRUE)
  e <- log(rbind(c(.99, .01), c(.01, .99), c(.9, .1)))
  together <- decode_states(e, p, segment = c("a", "b", "b"))
  a <- decode_states(e[1, , drop = FALSE], p)
  b <- decode_states(e[2:3, ], p)
  expect_equal(together$log_likelihood, a$log_likelihood + b$log_likelihood)
  expect_equal(together$posterior, rbind(a$posterior, b$posterior))
  expect_equal(together$state, c(a$state, b$state))
})

test_that("extreme finite log densities do not underflow", {
  p <- matrix(c(.9, .1, .1, .9), 2, byrow = TRUE)
  e <- matrix(c(-10000, -10001, -10002, -10000), 2, byrow = TRUE)
  a <- decode_states(e, p)
  b <- decode_states(e + 10000, p)
  expect_equal(a$log_likelihood + 20000, b$log_likelihood, tolerance = 1e-10)
  expect_equal(a$posterior, b$posterior, tolerance = 1e-10)
  expect_equal(a$state, b$state)
  expect_error(decode_states(matrix(-Inf, 1, 2), p), "zero or nonfinite")
})

test_that("importance emissions match hand calculation and remain stable", {
  d <- data.frame(ID = 1, stratum = c(1, 1, 2, 2),
                  obs = c(1, 0, 1, 0), w = c(1, 2, 1, 2), x = 1:4)
  b <- matrix(1, dimnames = list("x", "S1"))
  e <- ssf_log_emissions(d, ~ x, b)
  expect_equal(as.numeric(e), rep(log(exp(1) / (exp(1) + exp(2) / 2)), 2))
  e_large <- ssf_log_emissions(d, ~ x, b * 10000)
  expect_true(all(is.finite(e_large)))
  # Reordering cases and controls inside a choice set leaves its density intact.
  expect_equal(ssf_log_emissions(d[c(2, 1, 4, 3), ], ~ x, b), e)
})

test_that("one-state maximum agrees with an exact binomial reference", {
  n <- 20
  chosen <- c(rep(1, 5), rep(0, 15))
  d <- data.frame(ID = 1, stratum = rep(seq_len(n), each = 2),
                  obs = rep(c(1, 0), n), w = 1,
                  x = as.vector(rbind(chosen, 1 - chosen)))
  fit <- fit_ssf(d, ~ x, n_states = 1, short_iter = 10, long_iter = 100, seed = 93)
  expect_equal(unname(coef(fit)[1, 1]), log(5 / 15), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(fit)), 5 * log(.25) + 15 * log(.75), tolerance = 1e-8)
  expect_true(fit$diagnostics$converged)
  expect_true(fit$diagnostics$hessian_positive)
  expect_equal(fit$diagnostics$objective_agreement, 0)
  expect_equal(fit$diagnostics$long_start, fit$diagnostics$short_parameters)
})

test_that("kernel conventions have analytically distinct gamma shapes", {
  b <- rbind(step = c(-2, -3), log_step = c(0, -1.5), cos_angle = c(1, -.5))
  area <- movement_metrics(b, "area")
  radial <- movement_metrics(b, "step_angle")
  expect_equal(area$mean_step, c(1, 1 / 6))
  expect_equal(area$sd_step, c(sqrt(2) / 2, sqrt(.5) / 3))
  expect_equal(radial$valid, c(TRUE, FALSE))
  expect_true(is.na(radial$mean_step[2]))
  expect_equal(radial$mean_step[1], .5)
})
