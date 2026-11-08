family_parameters <- list(
  gamma = rbind(step_shape = c(2, 3), step_rate = c(4, 2)),
  lognormal = rbind(step_meanlog = c(-.9, .1), step_sdlog = c(.7, .9)),
  inverse_gaussian = rbind(step_mean = c(.5, 1.5), step_shape = c(.6, 2)),
  halfnormal = rbind(step_sd = c(.5, 1.5)),
  weibull = rbind(step_shape = c(1.5, 2), step_scale = c(.5, 1.5)),
  exponential = rbind(step_rate = c(2, .7))
)

family_steps <- function(measure) {
  prepare_ssf_steps(data.frame(ID = rep(1:2, each = 6L), step = seq(.2, 1.9, length.out = 12),
    angle = sin(seq_len(12)), previous_bearing = 0,
    target_bearing = seq(-pi, pi, length.out = 12),
    target_active = rep(c(TRUE, TRUE, FALSE), 4)), density_measure = measure)
}

family_controls <- function(steps, n, seed) {
  proposal <- c(shape = 1, rate = 1, mean_angle = 0, concentration = .1)
  generate_ssf_controls(steps, n, proposal, seed = seed)
}

test_that("step-family objectives and scores agree with independent emission calculations", {
  steps <- family_steps("area")
  case_control <- family_controls(steps, 6, 41)
  case_control$cover <- cos(case_control$step) > 0.95
  p <- matrix(c(.8, .2, .3, .7), 2, byrow = TRUE)
  for (family in names(family_parameters)) {
    for (exact in c(FALSE, TRUE)) {
      data <- if (exact) steps else case_control
      formula <- if (exact) ~ cos_angle + cos_diff_theta_to_target else ~ cos_angle + cover
      d <- .family_design(data, formula, family, exact)
      b <- rbind(family_parameters[[family]], c(.8, .3), c(.25, -.2))
      rownames(b) <- colnames(d$X)
      work <- .family_transform(c(b, qlogis(p[1, 2]), qlogis(p[2, 1])), d, 2, inverse = TRUE)
      objective <- .native_family_objective(d, 2, exact)
      reference <- function(theta) {
        unpacked <- .unpack(.family_transform(theta, d, 2), ncol(d$X), 2, colnames(d$X))
        log_e <- ssf_log_emissions(data, formula, unpacked$coefficients,
          if (exact) "exact" else "importance", step_family = family)
        -decode_states(log_e, unpacked$transition, d$segment)$log_likelihood
      }
      expect_equal(objective$fn(work), reference(work), tolerance = 1e-11, info = family)
      numerical <- vapply(seq_along(work), function(i) {
        plus <- minus <- work
        plus[i] <- plus[i] + 1e-5
        minus[i] <- minus[i] - 1e-5
        (reference(plus) - reference(minus)) / 2e-5
      }, numeric(1))
      expect_equal(objective$gr(work), numerical, tolerance = 1e-6, info = family)
      other_measure <- if (exact) family_steps("step_angle") else
        convert_density_measure(case_control, "step_angle")
      norm <- if (exact) "exact" else "importance"
      e <- ssf_log_emissions(data, formula, b, norm, step_family = family)
      e_polar <- ssf_log_emissions(other_measure, formula, b, norm, step_family = family)
      expect_equal(e_polar, e + log(steps$step), tolerance = 1e-12, info = family)
    }
  }
})

test_that("explicit gamma represents the original model and exponential is its unit-shape case", {
  steps <- family_steps("area")
  b <- rbind(family_parameters$gamma, cos_angle = c(.8, .3))
  legacy <- rbind(step = -b["step_rate", ], log_step = b["step_shape", ] - 2,
                  cos_angle = b["cos_angle", ])
  for (exact in c(TRUE, FALSE)) {
    data <- if (exact) steps else family_controls(steps, 8, 82)
    norm <- if (exact) "exact" else "importance"
    expect_equal(ssf_log_emissions(data, ~ cos_angle, b, norm, "gamma"),
      ssf_log_emissions(data, ~ step + log_step + cos_angle, legacy, norm), tolerance = 1e-12)
    b["step_shape", ] <- 1
    expect_equal(ssf_log_emissions(data, ~ cos_angle, b[-1L, , drop = FALSE], norm, "exponential"),
      ssf_log_emissions(data, ~ cos_angle, b, norm, "gamma"), tolerance = 1e-12)
    b["step_shape", ] <- family_parameters$gamma["step_shape", ]
  }
})

test_that("all families simulate with their fitted mean and valid positive steps", {
  segment <- data.frame(ID = "a", x = 0, y = 0, bearing = 0, n_steps = 6000)
  for (family in names(family_parameters)) {
    b <- family_parameters[[family]][, 1L, drop = FALSE]
    metrics <- movement_metrics(b, "area", family)
    fit <- structure(list(coefficients = b, density_measure = "area", n_states = 1L,
      normalization = "exact", step_family = family, transition = matrix(1)), class = "ssf_fit")
    sim <- simulate_ssf(fit, segment, seed = 711)$tracks[[1L]]
    expect_true(all(sim$step > 0) && all(is.finite(sim$step)), info = family)
    expect_true(abs(mean(sim$step) - metrics$mean_step) <
      5 * metrics$sd_step / sqrt(nrow(sim)), info = family)
  }
})

test_that("one taxonomy works across families, including a model with uniform turns", {
  steps <- family_steps("area")
  tax <- ssf_taxonomy("one_state", "Movement", formula = ~ 1)
  guided <- fit_guided_ssf(steps, "area", taxonomies = list(tax),
    step_family = "weibull", normalization = "exact", n_starts = 1L)
  expect_equal(guided$status, "selected")
  expect_equal(summary(guided)$states$persistence, 0)
  expect_match(capture.output(summary(guided)), "Step-length family: weibull", all = FALSE)
  restart <- diagnose_ssf_restarts(steps, tax, n_starts = 2L, seed = 1,
    step_family = "halfnormal")
  expect_equal(nrow(restart$solutions), 1)
  expect_equal(restart$solutions$n_restarts, 2)
  expect_error(fit_ssf(steps, ~ step, step_family = "gamma", normalization = "exact"),
    "angular or selection")
  expect_error(fit_ssf(steps, ~ target_active, step_family = "gamma", normalization = "exact"),
    "prospective target")
  expect_error(fit_ssf(steps, step_family = "invalid", normalization = "exact"),
    "step_family must be one of")
})
