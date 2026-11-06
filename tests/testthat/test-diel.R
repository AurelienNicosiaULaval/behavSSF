test_that("the diel score and decoder agree with independent calculations", {
  set.seed(413)
  hour <- c(22:23, 0:6, 12:19, 4L)
  steps <- prepare_ssf_steps(data.frame(ID = rep(1:3, c(9, 8, 1)),
    step = rgamma(length(hour), 2), angle = runif(length(hour), -pi, pi)), "area")
  formula <- ~ step + log_step + cos_angle + lag_log_step
  design <- .design_exact(steps, formula)
  basis <- .diel_basis(2L)
  design$diel <- list(hour = as.integer(hour), basis = basis)
  b <- cbind(c(-3, -.8, .4, .18), c(-2, .4, 1.3, .22))
  rownames(b) <- colnames(design$X)
  gamma <- cbind(c(-2, .3, -.4, .35, -.2), c(-1.3, -.45, .25, -.15, .15))
  theta <- .exact_transform(c(b, gamma), design, 2L, inverse = TRUE)
  native <- .native_exact_objective(design, 2L)
  reference <- function(par, decode = FALSE) {
    pars <- .exact_transform(par, design, 2L)
    b <- matrix(pars[1:8], 4L)
    rownames(b) <- colnames(design$X)
    schedule <- .diel_schedule(basis, matrix(pars[-(1:8)], 5L))
    .forward(.exact_emissions(design, b), schedule$transition, steps$ID,
      schedule$initial, decode = decode, hour = hour + 1L)
  }
  expect_equal(native$fn(theta), -reference(theta), tolerance = 1e-11)
  finite_diff <- vapply(seq_along(theta), function(j) {
    offset <- replace(numeric(length(theta)), j, 1e-5)
    -(reference(theta + offset) - reference(theta - offset))/2e-5
  }, numeric(1))
  expect_equal(native$gr(theta), finite_diff, tolerance = 1e-7)

  schedule <- .diel_schedule(basis, gamma)
  for (h in seq_len(24L)) {
    preceding <- (h + 22L) %% 24L + 1L
    expect_equal(drop(schedule$initial[, preceding] %*% schedule$transition[, , h]),
      unname(schedule$initial[, h]), ignore_attr = TRUE, tolerance = 1e-12)
  }
  constant <- matrix(0, 5L, 2L)
  constant[1L, ] <- gamma[1L, ]
  homogeneous <- .native_exact_objective(.design_exact(steps, formula), 2L)
  expect_equal(native$fn(c(theta[1:8], constant)), homogeneous$fn(c(theta[1:8], gamma[1L, ])),
    tolerance = 1e-11)

  # Enumerate one short path to check the transition hour in both decoding passes.
  keep <- 1:3
  e <- .exact_emissions(design, b)[keep, ]
  h <- hour[keep]+1L
  d <- .forward(e, schedule$transition, rep("a", 3), schedule$initial, TRUE, h)
  paths <- expand.grid(rep(list(1:2), 3))
  weight <- apply(paths, 1L, function(path)
    schedule$initial[path[1L], h[1L]] *
      prod(exp(e[cbind(1:3, path)])) *
      prod(schedule$transition[cbind(path[-3L], path[-1L], h[-1L])]))
  expect_equal(d$log_likelihood, log(sum(weight)))
  expect_equal(d$state, as.integer(unlist(paths[which.max(weight), ])))
  expect_equal(d$posterior[, 2], colSums((paths == 2)*weight)/sum(weight), ignore_attr = TRUE)
})

test_that("diel simulation uses the segment start hour and rolls through midnight", {
  steps <- data.frame(ID = c("a", "b"), x = 0, y = 0, bearing = 0,
    n_steps = c(3L, 2L), hour = c(23L, 4L))
  schedule <- array(rep(diag(2), 24), c(2L, 2L, 24L))
  schedule[, , 1L] <- diag(2)[, 2:1]
  initial <- matrix(rep(c(1, 0), 24), 2L)
  initial[, 5L] <- c(0, 1)
  fit <- structure(list(n_states = 2L, normalization = "exact", density_measure = "area",
    coefficients = rbind(step = c(-3, -1), log_step = c(0, 0), cos_angle = c(.5, 1)),
    transition = schedule, initial = initial, transition_coefficients = matrix(0, 3L, 2L)),
    class = "ssf_fit")
  expect_identical(simulate_ssf(fit, steps, seed = 1)$tracks[[1L]]$state, c(1L, 2L, 2L, 2L, 2L))
  expect_error(simulate_ssf(fit, steps[setdiff(names(steps), "hour")]), "step-origin hour")
  expect_error(.diel_hour(c(23, 1), 2L, c("same", "same")), "advance by one")
})
