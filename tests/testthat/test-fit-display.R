test_that("guided summaries use labels and report the candidate checks", {
  result <- display_fit()
  fit <- result$fit
  report <- summary(result)
  expect_s3_class(report, "summary_guided_ssf")
  expect_identical(report$states$state, "Other")
  expect_equal(report$states$n_steps, 41)
  expect_equal(report$states$fraction, 1)
  expect_equal(report$states$mean_step, result$assessment$metrics$mean_step)
  expect_identical(report$candidates$status, c("Rejected", "Selected"))
  expect_identical(report$candidates$reason, c("custom_check", ""))
  expect_identical(report$candidates$description, c("Moving", "Single state"))
  expect_identical(rownames(report$transition), "Other")
  expect_equal(report$coefficients$std_error,
               unname(sqrt(diag(fit$covariance))[seq_along(coef(fit))]))
  expect_output(print(result), "Other (Single state) [Selected]", fixed = TRUE)
  expect_output(print(report), "Other (Single state) [Selected]", fixed = TRUE)
  expect_identical(summary(fit)$states$state, colnames(coef(fit)))

  fit$covariance <- NULL
  fit$density_measure <- "unspecified"
  expect_true(all(is.na(summary(fit)$coefficients$std_error)))
  expect_null(summary(fit)$states$mean_step)
})

test_that("the reference codes are explained even when a fit failed", {
  failed <- list(error = "Optimisation failed.")
  rejected <- list(assessment = list(accepted = FALSE, taxonomy = "HDT_F_R",
    reasons = "target_positive"))
  result <- structure(list(history = list(HDT_F_R = rejected, T_F_R = failed,
    T_E = failed, Other = failed), status = "no_accepted_model"), class = "guided_ssf")
  report <- summary(result)
  expect_identical(report$candidates$description, c(
    "Habitat-directed travelling / Foraging / Resting",
    "Travelling / Foraging / Resting", "Travelling / Encamped", "Single state"
  ))
  expect_output(print(result), "HDT_F_R (Habitat-directed travelling / Foraging / Resting)", fixed = TRUE)
  expect_output(print(report), "95% confidence interval", fixed = TRUE)
  expect_output(print(report), "Optimisation failed.", fixed = TRUE)
})

test_that("an unsuccessful guided search remains inspectable", {
  result <- display_fit()
  # The first candidate was the complete history before the fallback.
  result$history <- result$history[1L]
  result$fit <- NULL
  result$assessment <- NULL
  result$taxonomy <- NULL
  result$status <- "no_accepted_model"
  expect_null(summary(result)$states)
  expect_identical(summary(result)$candidates$status, "Rejected")
  expect_error(plot(result), "No selected fit")
})

test_that("state plots preserve row alignment and gaps", {
  skip_if_not_installed("ggplot2")
  result <- display_fit()
  track <- plot(result, segment = "C", step_range = c(3, 5))
  expect_equal(track$data$x, result$fit$observations$x[24:26])
  expect_identical(as.character(track$data$state), rep("Other", 3))
  expect_equal(nrow(track$layers[[3L]]$data), 1)

  posterior <- plot(result, type = "posterior")
  layers <- ggplot2::ggplot_build(posterior)$data
  expect_length(unique(layers[[1L]]$group), 2)
  expect_equal(layers[[2L]]$x, 21)
  expect_equal(posterior$data$probability, as.vector(result$fit$posterior))

  result$fit$observations$x0 <- NULL
  expect_equal(plot(result)$data, posterior$data)
  expect_error(plot(result, type = "track"), "coords must identify")
  expect_error(plot(result, step_range = c(1, 42)), "step_range")
})

test_that("printed SSF covariates retain their intervals and ecological labels", {
  x <- readRDS(system.file("extdata", "condor-context-fit.rds", package = "behavSSF"))
  name <- c("cos_angle:altitude_origin_100m" = "Persistence per 100 m above reference altitude")
  report <- summary(x$context, assessment = x$context_assessment,
    term_labels = name, step_scale = 1000, step_unit = "m")
  expect_identical(report$states$state, x$context_assessment$labels)
  expect_identical(unique(report$ssf_effects$label), unname(name))
  z <- report$coefficients[report$coefficients$term == names(name), ]
  expect_equal(report$ssf_effects$lower, z$estimate - qnorm(.975) * z$std_error)
  expect_equal(report$ssf_effects$upper, z$estimate + qnorm(.975) * z$std_error)
  expect_equal(report$states$mean_step, summary(x$context)$states$mean_step * 1000)
  output <- capture.output(returned <- print(x$context, assessment = x$context_assessment,
    term_labels = name, step_scale = 1000, step_unit = "m"))
  expect_identical(returned, x$context)
  expect_true(any(grepl("Persistence per 100 m", output, fixed = TRUE)))
  expect_false(any(grepl("cos_angle:", output, fixed = TRUE)))
  expect_true(any(grepl("95% Wald CI", output, fixed = TRUE)))
  expect_output(print(report), "Baseline movement coefficients", fixed = TRUE)
  expect_output(print(x$exact), "No additional covariates", fixed = TRUE)
  x$context$covariance <- NULL
  expect_output(print(x$context), "CI unavailable", fixed = TRUE)
  expect_error(summary(x$context, term_labels = c(missing = "Unknown")), "term_labels")
  expect_error(summary(x$context, step_scale = 0), "step_scale")
})
