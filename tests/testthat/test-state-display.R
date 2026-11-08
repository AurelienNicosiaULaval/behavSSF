test_that("ecological summaries preserve decoded counts and movement units", {
  result <- display_fit()
  s <- ssf_state_summary(result, step_scale = 1000, step_unit = "m")
  expect_identical(s$state, "Other")
  expect_equal(s$n_steps, length(result$fit$state))
  expect_equal(s$fraction, 1)
  b <- coef(result$fit)
  expect_equal(s$mean_step, 1000 * (b["log_step", ] + 2) / -b["step", ])
  expect_equal(s$sd_step, 1000 * sqrt(b["log_step", ] + 2) / -b["step", ])
  expect_equal(s$stay_min, 1)
  expect_true(is.na(s$mean_run_steps))
  expect_identical(attr(s, "step_unit"), "m")
  expect_error(ssf_state_summary(result, step_scale = 0), "step_scale")
  expect_error(ssf_state_summary(result, step_unit = ""), "step_unit")
  expect_equal(ssf_state_summary(result$fit, assessment = result$assessment),
    ssf_state_summary(result))
  rejected <- result$assessment
  rejected$accepted <- FALSE
  expect_error(ssf_state_summary(result, assessment = rejected), "same fitted model")
  different <- result$assessment
  different$context$coefficients[1, 1] <- 0
  expect_error(ssf_state_summary(result, assessment = different), "same fitted model")
  result$fit <- NULL
  expect_error(ssf_state_summary(result), "No selected fit")
})

test_that("state panels retain the plotted observations and conditional intervals", {
  skip_if_not_installed("ggplot2")
  result <- display_fit()
  track <- plot_ssf_states(result, panel = "track", segment = "C", step_range = c(3, 5),
    colours = c(Other = "#21765A"), coord_unit = "m")
  expect_equal(track$data$x, result$fit$observations$x[24:26])
  expect_identical(track$labels$x, "Easting (m)")
  expect_s3_class(plot_ssf_states(result, panel = "table"), "ggplot")
  effects <- plot_ssf_states(result, panel = "coefficients", terms = "cos_angle")
  expect_equal(effects$data$term, factor("cos_angle", levels = "cos_angle"))
  expect_equal(effects$data$lower, effects$data$estimate - qnorm(.975) * effects$data$std_error)
  expect_error(plot_ssf_states(result, terms = "absent"), "coefficient names")
  expect_error(plot_ssf_states(result, colours = "invalid"), "Invalid state colour")
  expect_error(plot_ssf_states(result, colours = c(Another = "red")), "every fitted state")
  result$fit$covariance <- NULL
  expect_true(is.na(plot_ssf_states(result, panel = "coefficients", terms = "cos_angle")$data$lower))
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggcircular")
  display <- plot_ssf_states(result, panel = "distributions")
  length_panel <- display[[1]]
  angle_panel <- display[[2]]
  expect_equal(length_panel$data$log_distance, log10(result$fit$observations$step))
  expect_equal(angle_panel$data$angle, result$fit$observations$angle %% (2 * pi))
  built <- ggplot2::ggplot_build(angle_panel)$data[[1]]
  expect_equal(built$y[1], tail(built$y, 1), tolerance = 1e-12)
  expect_error(plot_ssf_states(result, angle_bandwidth = 0), "angle_bandwidth")
  expect_s3_class(plot(result, type = "states"), "patchwork")
  expect_s3_class(plot_ssf_states(result, panel = "distributions"), "patchwork")
})


test_that("default coefficient views show covariates without repeating baseline movement", {
  skip_if_not_installed("ggplot2")
  x <- readRDS(system.file("extdata", "condor-context-fit.rds", package = "behavSSF"))
  terms <- as.character(plot_ssf_states(x$context, panel = "coefficients")$data$term)
  expect_identical(unique(terms), "cos_angle:altitude_origin_100m")
  expect_identical(ssf_state_contrasts(x$context)$term, "cos_angle:altitude_origin_100m")
  expect_identical(ssf_state_contrasts(x$context, terms = "cos_angle")$term, "cos_angle")
  expect_equal(nrow(ssf_state_contrasts(x$exact)), 0L)
  explicit <- plot_ssf_states(x$exact, panel = "coefficients", terms = "cos_angle")
  expect_true(all(as.character(explicit$data$term) == "cos_angle"))
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggcircular")
  expect_length(plot_ssf_states(x$exact), 4L)
  expect_length(plot_ssf_states(x$context), 5L)
})
