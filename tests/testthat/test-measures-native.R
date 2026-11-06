test_that("native objective and gradient agree with the R reference", {
  set.seed(261104)
  d <- data.frame(ID = rep(c("a", "b", "c"), c(3, 6, 12)),
    stratum = rep(1:7, each = 3), obs = rep(c(1, 0, 0), 7),
    w = exp(rnorm(21)), x = rnorm(21), z = rnorm(21))
  design <- behavSSF:::.design(d, ~ x + z)
  for (k in c(1L, 2L, 3L, 6L)) {
    theta <- rnorm(2*k + k*(k-1), sd = .4)
    native <- behavSSF:::.native_objective(design, k)
    reference <- function(v) {
      par <- behavSSF:::.unpack(v, 2L, k, colnames(design$X))
      -decode_states(behavSSF:::.emissions(design, par$coefficients),
        par$transition, design$segment)$log_likelihood
    }
    h <- 1e-6
    differences <- vapply(seq_along(theta), function(i) {
      plus <- minus <- theta
      plus[i] <- plus[i] + h; minus[i] <- minus[i] - h
      (reference(plus) - reference(minus)) / (2*h)
    }, numeric(1))
    expect_equal(native$fn(theta), reference(theta), tolerance = 1e-11)
    expect_equal(native$gr(theta), differences, tolerance = 2e-6)
    if (k > 1L) {
      bad <- theta; bad[-seq_len(2*k)] <- -1000
      expect_equal(native$fn(bad), 1e100)
      expect_equal(native$fn(theta), reference(theta), tolerance = 1e-11)
      expect_gt(native$diagnostics()$invalid_evaluations, 0)
    }
  }
})

test_that("combined angular parameters agree with numerical integration", {
  result <- angular_parameters(2, 1, pi / 2)
  expect_equal(result$concentration, sqrt(5))
  expect_equal(result$mean_angle, atan2(1, 2))
  density <- function(phi) exp(2*cos(phi) + cos(phi-pi/2))
  normalizer <- stats::integrate(density, -pi, pi)$value
  mx <- stats::integrate(function(phi) cos(phi)*density(phi), -pi, pi)$value / normalizer
  my <- stats::integrate(function(phi) sin(phi)*density(phi), -pi, pi)$value / normalizer
  resultant <- besselI(sqrt(5), 1) / besselI(sqrt(5), 0)
  expect_equal(c(mx, my), resultant*c(2, 1)/sqrt(5), tolerance = 1e-10)
  expect_true(is.na(angular_parameters(0)$mean_angle))
  expect_equal(angular_parameters(-2)$concentration, 2)
})

test_that("density and coefficient conversions preserve the spatial model", {
  d <- data.frame(ID = rep(1:2, each = 6), stratum = rep(1:4, each = 3),
    obs = rep(c(1, 0, 0), 4), w = seq(.3, 1.4, length.out = 12),
    step = c(.1, .4, .9, .3, .5, 1.1, .2, .7, 1.2, .6, .8, 1.5),
    cos_angle = cos(seq_len(12)))
  area <- prepare_ssf_data(d, density_measure = "area")
  radial <- convert_density_measure(area, "step_angle")
  beta <- rbind(step = c(-2, -4), log_step = c(-.2, -.7), cos_angle = c(1, -.5))
  radial_beta <- convert_movement_coefficients(beta, "area", "step_angle")
  e_area <- ssf_log_emissions(area, ~ step + log_step + cos_angle, beta)
  e_radial <- ssf_log_emissions(radial, ~ step + log_step + cos_angle, radial_beta)
  log_r <- log(area$step[area$obs == 1])
  expect_equal(e_radial, e_area + log_r, tolerance = 1e-12)
  transition <- matrix(c(.8, .2, .3, .7), 2, byrow = TRUE)
  a <- decode_states(e_area, transition, segment = rep(1:2, each = 2))
  b <- decode_states(e_radial, transition, segment = rep(1:2, each = 2))
  expect_equal(b$log_likelihood - a$log_likelihood, sum(log_r), tolerance = 1e-12)
  expect_equal(a$posterior, b$posterior, tolerance = 1e-12)
  expect_identical(a$state, b$state)
  expect_equal(movement_metrics(beta, "area")$mean_step,
                movement_metrics(radial_beta, "step_angle")$mean_step)
  expect_equal(convert_density_measure(radial, "area")$w, area$w)
  expect_error(convert_density_measure(prepare_ssf_data(d), "area"), "Annotate")
})

test_that("annotated measures prevent inconsistent ecological assessment", {
  d <- data.frame(ID = 1, stratum = rep(1:3, each = 2),
                  obs = rep(c(1, 0), 3), w = 1,
                  step = c(.1, .3, .4, .2, .7, .9), cos_angle = cos(1:6))
  expect_error(fit_guided_ssf(d, convention = "area"), "Annotate")
  f <- structure(list(density_measure = "area"), class = "ssf_fit")
  expect_error(assess_taxonomy(f, "Other", "step_angle"), "must match")
  f$density_measure <- "unspecified"
  expect_error(assess_taxonomy(f, "Other", "area"), "explicit fitted density measure")
})
