test_that("lagged exact densities nest the original kernel and reset by segment", {
  track <- prepare_ssf_steps(data.frame(ID = c("a", "a", "b", "b"),
    step = c(.2, .1, 2, .3), angle = c(0, .5, -.2, .4)), "area")
  base <- rbind(step = -5, log_step = 0, cos_angle = 0)
  f <- ~ step + log_step + cos_angle + lag_log_step
  beta <- rbind(base, lag_log_step = 0)
  original <- ssf_log_emissions(track, ~ step + log_step + cos_angle, base, "exact")
  expect_equal(ssf_log_emissions(track, f, beta, "exact"), original, tolerance = 1e-13)

  beta["lag_log_step", ] <- .6
  conditional_mean <- .4 * (c(.4, .2, .4, 2) / .4)^.6
  expected <- dgamma(track$step, shape = 2, rate = 2 / conditional_mean, log = TRUE) -
    log(track$step) - log(2 * pi)
  expect_equal(as.vector(ssf_log_emissions(track, f, beta, "exact")), expected, tolerance = 1e-12)
  polar <- convert_density_measure(track, "step_angle")
  b_polar <- convert_movement_coefficients(beta, "area", "step_angle")
  expect_equal(as.vector(ssf_log_emissions(polar, f, b_polar, "exact")),
    expected + log(track$step), tolerance = 1e-12)
})

test_that("native AR score agrees with the independently normalised R likelihood", {
  track <- prepare_ssf_steps(data.frame(ID = rep(c("a", "b"), each = 3L),
    step = c(.4, .2, .07, 1.8, .3, .6), angle = seq(-1, 1, length.out = 6),
    relative_target_bearing = seq(1, -.8, length.out = 6)), "area")
  f <- ~ step + log_step + lag_log_step + cos_angle + cos_diff_theta_to_target
  design <- behavSSF:::.design_exact(track, f)
  beta <- c(-5, .2, .4, .3, .2, -2, .5, -.2, 1, -.3, -1, -2)
  theta <- behavSSF:::.exact_transform(beta, design, 2L, inverse = TRUE)
  native <- behavSSF:::.native_exact_objective(design, 2L)
  reference <- function(x) {
    b <- behavSSF:::.exact_transform(x, design, 2L)
    model <- behavSSF:::.unpack(b, 5L, 2L, colnames(design$X))
    emissions <- ssf_log_emissions(track, f, model$coefficients, "exact")
    -decode_states(emissions, model$transition, track$ID)$log_likelihood
  }
  expect_equal(native$fn(theta), reference(theta), tolerance = 1e-11)
  h <- 1e-5
  score <- vapply(seq_along(theta), function(i) {
    plus <- minus <- theta
    plus[i] <- plus[i] + h
    minus[i] <- minus[i] - h
    (reference(plus) - reference(minus)) / (2 * h)
  }, numeric(1))
  expect_equal(native$gr(theta), score, tolerance = 1e-7)
})

test_that("AR simulation feeds back generated lengths and restarts the lag", {
  coefficients <- rbind(step = c(-5, -3), log_step = c(0, 1), cos_angle = c(0, 0),
    lag_log_step = c(.3, .7))
  fit <- structure(list(coefficients = coefficients, n_states = 2L,
    transition = matrix(.5, 2, 2), density_measure = "area", normalization = "exact"),
    class = "ssf_fit")
  segments <- data.frame(ID = c("a", "b"), x = 0, y = 0, bearing = 0, n_steps = 2L)
  path <- c(1L, 2L, 2L, 1L)
  simulated <- simulate_ssf(fit, segments, state_path = path, seed = 91)$tracks[[1]]
  expected <- behavSSF:::.with_seed(91, function() {
    r <- numeric(4)
    for (t in 1:4) {
      state <- path[t]
      shape <- coefficients["log_step", state] + 2
      rate <- -coefficients["step", state]
      anchor <- shape / rate
      adjustment <- if (t %% 2L == 1L) 1 else (r[t - 1L] / anchor)^coefficients["lag_log_step", state]
      stats::runif(2L)
      r[t] <- stats::rgamma(1L, shape, scale = adjustment / rate)
    }
    r
  })
  expect_equal(simulated$step, expected, tolerance = 1e-12)
})

test_that("the AR density remains numerically stable for a concentrated trial kernel", {
  track <- prepare_ssf_steps(data.frame(ID = "a", step = 1 + c(0, 1e-9, -1e-9),
    angle = 0), "area")
  formula <- ~ step + log_step + cos_angle + lag_log_step
  design <- behavSSF:::.design_exact(track, formula)
  theta <- c(log(1e18), log(1e18), 0, 0)
  beta <- matrix(behavSSF:::.exact_transform(theta, design, 1L), ncol = 1,
    dimnames = list(colnames(design$X), NULL))
  reference <- dgamma(track$step, shape = -beta["step", 1], rate = -beta["step", 1],
    log = TRUE) - log(track$step) - log(2*pi)
  expect_equal(as.vector(ssf_log_emissions(track, formula, beta, "exact")), reference, tolerance = 1e-10)
  native <- behavSSF:::.native_exact_objective(design, 1L)
  expect_equal(native$fn(theta), -sum(reference), tolerance = 1e-8)
})
