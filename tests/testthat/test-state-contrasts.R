test_that("state contrasts use covariance between arbitrary state effects", {
  fit <- display_fit()$fit
  b <- cbind(c(step = -2, log_step = 0, cos_angle = .2, other_variable = .3),
    c(step = -1, log_step = 0, cos_angle = 1.2, other_variable = .9))
  colnames(b) <- c("Local", "Directed")
  fit$coefficients <- b
  fit$n_states <- 2L
  fit$transition <- matrix(c(.8, .2, .2, .8), 2L)
  fit$covariance <- diag(.25, 8L)
  fit$covariance[4, 8] <- fit$covariance[8, 4] <- .15
  x <- ssf_state_contrasts(fit, terms = "other_variable")
  expect_equal(x$estimate, .6)
  expect_equal(x$std_error, sqrt(.25 + .25 - 2 * .15))
  expect_equal(x$lower, .6 - qnorm(.975) * sqrt(.2))
  reverse <- ssf_state_contrasts(fit, terms = "other_variable",
    comparisons = data.frame(from = "Directed", to = "Local"))
  expect_equal(reverse$estimate, -x$estimate)
  expect_equal(reverse$lower, -x$upper)
  expect_identical(ssf_state_contrasts(fit)$term, "other_variable")
  expect_error(ssf_state_contrasts(fit, "unknown"), "coefficient names")
  expect_error(ssf_state_contrasts(fit, comparisons = data.frame(from = "Local", to = "Local")), "from/to")
  expect_error(ssf_state_contrasts(fit, level = 1), "level")
  fit$covariance <- NULL
  expect_true(all(is.na(ssf_state_contrasts(fit)$lower)))
  skip_if_not_installed("ggplot2")
  p <- plot_ssf_states(fit, panel = "contrasts", terms = "other_variable")
  expect_equal(p$data$estimate, .6)
  expect_identical(as.character(p$data$comparison), "Directed - Local")
  expect_error(ssf_state_contrasts(display_fit()), "At least two")
})
