toy_data <- function() data.frame(ID = 1, stratum = rep(1:3, each = 2),
                                  obs = rep(c(1, 0), 3), w = 1, x = 1:6)

test_that("invalid choice sets and sampling densities are rejected", {
  d <- toy_data()
  d$obs[2] <- 1
  expect_error(prepare_ssf_data(d), "exactly one")
  d <- toy_data(); d$w[2] <- 0
  expect_error(prepare_ssf_data(d), "strictly positive")
  d <- toy_data(); d$obs[1] <- NA
  expect_error(prepare_ssf_data(d), "only 0 and 1")
  d <- toy_data(); d$x[3] <- Inf
  expect_error(fit_ssf(d, ~ x, 1), "finite covariates")
  d <- toy_data(); d$constant <- 1
  expect_error(fit_ssf(d, ~ constant, 1), "not identifiable")
})

test_that("time and ID contracts protect independent segments", {
  d <- toy_data(); d$ID <- rep(c("a", "b", "a"), each = 2)
  expect_error(prepare_ssf_data(d), "contiguous")
  d <- toy_data(); d$time <- rep(c(1, 3, 2), each = 2)
  expect_error(prepare_ssf_data(d, time = "time"), "strictly increasing")
  d <- toy_data(); d$time <- 1:6
  expect_error(prepare_ssf_data(d, time = "time"), "share its case time")
  d <- toy_data(); d$time <- rep(c(0, 1800, 9000), each = 2)
  expect_error(prepare_ssf_data(d, time = "time", expected_interval = 1800), "time gap")
  d <- toy_data(); d$ID <- c(1, 1, 2, 2, 2, 2); d$stratum <- c(1, 1, 1, 1, 2, 2)
  expect_equal(length(unique(prepare_ssf_data(d)$stratum)), 3L)
})

test_that("seeded fits restore caller RNG state and reproduce", {
  d <- toy_data(); d$x <- c(1, 0, 0, 1, 1, 0)
  set.seed(417)
  previous <- .Random.seed
  a <- fit_ssf(d, ~ x, 1, seed = 94)
  expect_identical(.Random.seed, previous)
  b <- fit_ssf(d, ~ x, 1, seed = 94)
  expect_equal(a$coefficients, b$coefficients)
  expect_equal(a$starts, b$starts)
  expect_error(fit_ssf(d, ~ x, 1, seed = -1), "seed")
})

fixture_fit <- function() {
  b <- rbind(step = c(-1, -2, -3), log_step = c(1, 0, -.5),
              cos_angle = c(3, 1, .2), cos_diff_theta_to_target = c(2, 0, 0))
  structure(list(coefficients = b, n_states = 3L, density_measure = "area",
    transition = matrix(c(.8, .1, .1, .1, .8, .1, .1, .1, .8), 3),
    covariance = diag(.01, 18), state = rep(1:3, each = 10),
    diagnostics = list(converged = TRUE, hessian_positive = TRUE,
                        max_abs_gradient = 0)), class = "ssf_fit")
}

test_that("taxonomy rules are invariant to state numbering", {
  f <- fixture_fit()
  a <- assess_taxonomy(f, "HDT_F_R", "area")
  expect_true(a$accepted)
  perm <- c(3, 1, 2)
  f$coefficients <- f$coefficients[, perm]
  f$transition <- f$transition[perm, perm]
  f$state <- match(f$state, perm)
  b <- assess_taxonomy(f, "HDT_F_R", "area")
  expect_true(b$accepted)
  expect_equal(b$labels, a$labels[perm])
})

test_that("occupancy, negative attraction and numerical failure are explicit", {
  f <- fixture_fit(); f$state <- c(rep(1L, 98), 2L, 3L)
  a <- assess_taxonomy(f, "HDT_F_R", "area")
  expect_false(a$accepted)
  expect_contains(a$reasons, "occupancy")
  f <- fixture_fit(); f$coefficients[4, 1] <- -2
  a <- assess_taxonomy(f, "HDT_F_R", "area")
  expect_contains(a$reasons, "target_positive")
  f <- fixture_fit(); f$diagnostics$converged <- FALSE
  expect_contains(assess_taxonomy(f, "T_F_R", "area")$reasons, "converged")
  f <- fixture_fit(); f$coefficients[2, ] <- -3
  expect_contains(assess_taxonomy(f, "T_F_R", "area")$reasons, "valid_kernel")
})

test_that("custom taxonomies rank states and expose consistently labelled rules", {
  f <- fixture_fit()
  candidate <- ssf_taxonomy("movement_modes", c("Local", "Intermediate", "Directed"),
    order_by = function(metrics, coefficients) coefficients["cos_angle", ],
    rules = function(s) {
      c(contrast = s$metrics["Directed", "mean_step"] > s$metrics["Local", "mean_step"],
        precise = s$standard_errors["step", "Directed"] < .2,
        transition = s$transition["Directed", "Directed"] > .5)
    })
  expect_error(assess_taxonomy(f, candidate, "area"), "distinct from the built-in checks")
  candidate$rules <- function(s) {
    c(contrast = s$metrics["Directed", "mean_step"] > s$metrics["Local", "mean_step"],
      precise = s$standard_errors["step", "Directed"] < .2,
      residence = s$transition["Directed", "Directed"] > .5,
      occupied = s$occupancy[["Local"]] > .2)
  }
  a <- assess_taxonomy(f, candidate, "area")
  expect_true(a$accepted)
  expect_identical(a$labels, c("Directed", "Intermediate", "Local"))
  expect_identical(a$taxonomy, "movement_modes")
  perm <- c(3L, 1L, 2L)
  f$coefficients <- f$coefficients[, perm]
  f$transition <- f$transition[perm, perm]
  f$state <- match(f$state, perm)
  expect_identical(assess_taxonomy(f, candidate, "area")$labels, a$labels[perm])
  f$covariance <- NULL
  expect_contains(assess_taxonomy(f, candidate, "area")$reasons, "precise")
  f$coefficients["cos_angle", ] <- 0
  b <- assess_taxonomy(f, candidate, "area")
  expect_identical(b$reasons, "state_order")
  expect_true(all(is.na(b$labels)))
})

test_that("the guided workflow accepts an ordered custom sequence and fallback", {
  set.seed(231)
  s <- prepare_ssf_steps(data.frame(ID = 1, step = rgamma(150, 2, 3),
    angle = runif(150, -pi, pi)), "area")
  reject <- ssf_taxonomy("declined", "Local", rules = function(s) c(contrast = FALSE))
  keep <- ssf_taxonomy("accepted", "Local")
  fit <- fit_guided_ssf(s, "area", normalization = "exact", step_ar = TRUE,
    taxonomies = list(reject, keep, "Other"), n_starts = 1, seed = 271)
  expect_identical(fit$taxonomy, "accepted")
  expect_identical(fit$assessment$labels, "Local")
  expect_identical(names(fit$history), c("declined", "accepted"))
  expect_contains(fit$history$declined$assessment$reasons, "contrast")
  expect_contains(rownames(coef(fit$fit)), "lag_log_step")
  fit <- fit_guided_ssf(s, "area", normalization = "exact",
    taxonomies = list(reject, "Other"), n_starts = 1, seed = 271)
  expect_identical(fit$taxonomy, "Other")
  expect_error(fit_guided_ssf(s, "area", normalization = "exact",
    taxonomies = list(keep, keep)), "names must be unique")
  expect_error(fit_guided_ssf(s, "area", normalization = "exact",
    taxonomies = list("HDT_F_R")), "requires a target")
  expect_error(ssf_taxonomy("Other", "Local"), "built-in")
  expect_error(ssf_taxonomy("local", c("A", " A ")), "distinct nonempty")
})

test_that("malformed HMM inputs fail without substituting emission densities", {
  expect_error(decode_states(matrix(NA_real_, 2, 2), diag(2)), "without NA")
  expect_error(decode_states(matrix(0, 2, 2), matrix(.6, 2, 2)), "row-stochastic")
  expect_error(decode_states(matrix(0, 3, 2), matrix(.5, 2, 2), c(1, 2, 1)), "contiguous")
  expect_error(decode_states(matrix(0, 2, 2), diag(2)), "singular")
})
