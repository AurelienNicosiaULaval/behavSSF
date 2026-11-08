test_that("restart solutions align state labels and keep failed starts visible", {
  gps <- readRDS(system.file("extdata", "deer-restart-steps.rds",
    package = "behavSSF", mustWork = TRUE))
  steps <- prepare_ssf_steps(gps, "area", time = "time_utc")
  taxonomy <- ssf_taxonomy("step_lengths", c("Short steps", "Long steps"))
  fit <- fit_ssf(steps, normalization = "exact", n_starts = 8, short_iter = 100,
    long_iter = 900, seed = 9871)
  reversed <- list(coefficients = fit$coefficients[, 2:1],
    transition = fit$transition[2:1, 2:1])
  set.seed(41)
  previous <- .Random.seed
  result <- diagnose_ssf_restarts(steps, taxonomy, seed = 60200, starts = list(
    default = list(seed = 60002),
    near_A = list(initial_values = list(coefficients = fit$coefficients, transition = fit$transition)),
    near_A = list(initial_values = reversed),
    invalid = list(initial_values = list(coefficients = matrix(0, 1, 2)))
  ))
  expect_identical(.Random.seed, previous)
  expect_identical(result$restarts$solution, c(2L, 1L, 1L, NA_integer_))
  expect_identical(result$solutions$n_restarts, c(2L, 1L))
  expect_equal(result$solutions$log_likelihood, c(59.65201, 52.38869), tolerance = 1e-6)
  expect_identical(result$restarts$status[4], "failed")
  expect_match(result$restarts$reason[4], "Initial coefficients")
  expect_true(all(vapply(result$attempts[1:3], function(log) "long" %in% log$phase,
    logical(1))))

  comparison <- compare_ssf_solutions(result)
  expect_identical(comparison$steps$end_fix, gps$end_fix)
  expect_identical(comparison$steps$time_utc, gps$time_utc)
  expect_identical(result$solutions$n_changed_vs_best, c(0L, 325L))
  expect_equal(comparison$transitions["Short steps", "Long steps"], 325)
  expect_equal(comparison$transitions["Long steps", "Short steps"], 0)
  expect_equal(unname(rowSums(comparison$posterior_from)), rep(1, nrow(gps)))
  expect_identical(colnames(comparison$posterior_to), taxonomy$labels)
  reverse <- compare_ssf_solutions(result, from = 2, to = 1)
  expect_identical(reverse$steps$.changed, comparison$steps$.changed)
  expect_equal(reverse$transitions["Long steps", "Short steps"], 325)
  expect_error(compare_ssf_solutions(result, from = 1, to = 1), "distinct")

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    map <- plot_ssf_changes(comparison, segment = "segment-14")
    layers <- ggplot2::ggplot_build(map)$data
    expect_equal(nrow(layers[[1]]), 94)
    expect_equal(nrow(layers[[2]]), 50)
    expect_error(plot_ssf_changes(comparison, coords = letters[1:4]), "coords")
  }
})

test_that("zero admissible solutions leave an inspectable restart log", {
  gps <- readRDS(system.file("extdata", "deer-restart-steps.rds",
    package = "behavSSF", mustWork = TRUE))
  steps <- prepare_ssf_steps(gps, "area")
  taxonomy <- ssf_taxonomy("rejected", "Any", rules = function(s) c(example = FALSE))
  result <- diagnose_ssf_restarts(steps, taxonomy, n_starts = 1, seed = 60200, keep_rejected = TRUE)
  expect_equal(nrow(result$solutions), 0)
  expect_identical(result$restarts$status, "rejected")
  expect_match(result$restarts$reason, "example")
  expect_error(compare_ssf_solutions(result), "At least two")
  expect_s3_class(result$all_fits[[1]], "ssf_fit")
  expect_false(result$all_assessments[[1]]$accepted)
  relaxed <- ssf_taxonomy("retained", "Any")
  reassessed <- assess_taxonomy(result$all_fits[[1]], relaxed, "area")
  expect_true(reassessed$accepted)

})
