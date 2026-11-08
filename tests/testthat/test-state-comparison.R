comparison_steps <- function() {
  set.seed(731)
  prepare_ssf_steps(data.frame(ID = rep(c("a", "b", "c"), each = 35),
    step = rgamma(105, 2, 5), angle = runif(105, -pi, pi), stratum = 1:105, previous_bearing = 0, target_bearing = 0), "area")
}

test_that("validation blocks preserve whole segments and constrain unbalanced cuts", {
  d <- comparison_steps()
  expect_identical(ssf_blocks(d), c(a = 1L, b = 2L, c = 3L))
  big <- d[c(rep(1, 100), 36, 71), ]
  expect_identical(ssf_blocks(big), c(a = 1L, b = 2L, c = 3L))
  expect_error(ssf_blocks(d, 4), "n_blocks")
  expect_error(compare_ssf_states(d, list(ssf_taxonomy("one", "Movement")),
    blocks = c(a = 1, b = 2)), "every segment")
})

test_that("held-out fits use training rows only and scores integrate the states", {
  d <- comparison_steps()
  tax <- ssf_taxonomy("one", "Movement", formula = ~ 1)
  seen <- list()
  initializer <- function(data, taxonomy, step_family, seed, block) {
    seen[[as.character(block)]] <<- unique(data$ID)
    list(anchor = list(coefficients = rbind(step_shape = 2, step_rate = 5),
      transition = matrix(1)))
  }
  set.seed(132)
  before <- .Random.seed
  x <- compare_ssf_states(d, list(tax), step_families = "gamma",
    normalization = "exact", initializers = initializer, seed = 9)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(x$start_profiles),4L)
  expect_equal(nrow(x$start_coefficients),8L)
  expect_true(all(x$start_profiles$kernel_valid))
  expect_equal(x$start_profiles$length_rank,rep(1L,4L))
  expect_setequal(seen[["0"]], c("a", "b", "c"))
  expect_setequal(seen[["1"]], c("b", "c"))
  expect_setequal(seen[["2"]], c("a", "c"))
  expect_setequal(seen[["3"]], c("a", "b"))
  # For uniform turns, Cartesian density is gamma(r) / (2*pi*r).
  expected <- vapply(1:3, function(block) {
    train <- d[d$ID != letters[block], ]
    test <- d[d$ID == letters[block], ]
    f <- fit_ssf(train, ~ 1, n_states = 1L, n_starts = 1,
      normalization = "exact", step_family = "gamma",
      initial_values = initializer(train, tax, "gamma", 1, 0)[[1L]])
    b <- coef(f)
    sum(dgamma(test$step, b["step_shape", ], b["step_rate", ], log = TRUE) -
      log(2 * pi * test$step))
  }, numeric(1))
  expect_equal(x$block_scores$log_score, expected, tolerance = 1e-7)
  expect_equal(x$candidates$mean_log_score, sum(expected) / nrow(d), tolerance = 1e-7)
  expect_equal(sum(x$segment_scores$steps), nrow(d))
  expect_output(print(x), "3 held-out blocks")
  expect_output(print(summary(x)), "Scores by reserved block")
  if (requireNamespace("ggplot2", quietly = TRUE)) expect_s3_class(plot(x), "ggplot")
  if (.Platform$OS.type == "unix") {
    y <- compare_ssf_states(d, list(tax), step_families = "gamma",
      normalization = "exact", initializers = initializer, seed = 9, workers = 2)
    expect_equal(y$candidates, x$candidates, tolerance = 1e-10)
    expect_identical(y$starts$seed, x$starts$seed)
    expect_identical(.Random.seed, before)
  }
})

test_that("incomplete coverage has no pooled score and ecology does not select starts", {
  d <- comparison_steps()
  good <- ssf_taxonomy("good", "Movement", formula = ~ 1)
  rejected <- ssf_taxonomy("rejected", "Movement", formula = ~ 1,
    rules = function(s) c(example_rule = FALSE))
  init <- function(data, taxonomy, step_family, seed, block) {
    if (block == 2L) stop("deliberate training failure")
    list(anchor = list(coefficients = rbind(step_shape = 2, step_rate = 5),
      transition = matrix(1)))
  }
  x <- compare_ssf_states(d, list(good, rejected), step_families = "gamma",
    normalization = "exact", initializers = init, score_failed_hessian = TRUE)
  expect_identical(x$candidates$blocks_scored, c(2L, 2L))
  expect_true(all(is.na(x$candidates$mean_log_score)))
  expect_true(all(is.na(x$candidates$point_mean_log_score)))
  expect_true(all(x$candidates$numerical_passed))
  expect_identical(x$candidates$taxonomy_passed, c(TRUE, FALSE))
  expect_match(x$candidates$checks_failed[2], "example_rule")
  expect_equal(coef(x$fits[[1]]), coef(x$fits[[2]]), tolerance = 1e-7)
  expect_true(all(grepl("deliberate", x$starts$error[x$starts$block == 2])))
})

test_that("conditional choice scores account for weights and unequal control counts", {
  steps <- comparison_steps()
  d <- generate_ssf_controls(steps, 3L,
    c(shape = 1, rate = 1, mean_angle = 0, concentration = 0), seed = 51)
  d <- d[-which(d$obs == 0L)[1L], ]
  fit <- fit_ssf(d, ~ 1, step_family = "gamma", n_states = 1,
    n_starts = 1, normalization = "importance", seed = 1)
  result <- .comparison_score(fit, d)
  # Independent categorical probabilities, with no /J remaining.
  b <- coef(fit)
  log_f <- dgamma(d$step, b["step_shape", ], b["step_rate", ], log = TRUE) -
    log(2 * pi * d$step)
  expected <- sum(vapply(split(seq_len(nrow(d)), d$stratum), function(rows) {
    z <- log_f[rows] - log(d$w[rows])
    used <- which(d$obs[rows] == 1L)
    z[used] - (max(z) + log(sum(exp(z - max(z)))))
  }, numeric(1)))
  expect_equal(sum(result$log_score), expected, tolerance = 1e-10)
})

test_that("a failed Hessian has a separate point score and cannot enter the strict pool", {
  d <- comparison_steps()
  original <- fit_ssf
  testthat::local_mocked_bindings(fit_ssf = function(data, ...) {
    f <- original(data, ...)
    if (!"b" %in% data$ID) f$diagnostics$hessian_positive <- FALSE
    f
  }, .package = "behavSSF")
  x <- compare_ssf_states(d, list(ssf_taxonomy("one", "Movement", formula = ~ 1)),
    step_families = "gamma", normalization = "exact", n_starts = 1,
    score_failed_hessian = TRUE)
  expect_identical(x$candidates$blocks_scored, 2L)
  expect_true(is.na(x$candidates$mean_log_score))
  expect_true(is.finite(x$candidates$point_mean_log_score))
  expect_false(x$block_scores$hessian_positive[x$block_scores$screen == "point" &
    x$block_scores$block == 2])
})

test_that("maps distinguish subdivision, changed parents and uncertain labels", {
  d <- as.data.frame(comparison_steps()[1:8, ])
  d$x0 <- 1:8; d$y0 <- 0; d$x <- 2:9; d$y <- 1
  a <- list(observations = d, state = c(1, 1, 1, 1, 2, 2, 2, 2))
  b <- list(observations = d, state = c(1, 1, 2, 3, 2, 3, 3, 3),
    posterior = rbind(c(.8,.1,.1), c(.8,.1,.1), c(.2,.6,.2), c(.1,.1,.8),
      c(.1,.8,.1), c(.1,.1,.8), c(.1,.1,.8), c(.1,.1,.8)))
  x <- structure(list(candidates = data.frame(candidate = c("a", "b"),
    taxonomy = c("two", "three")), fits = list(a = a, b = b),
    observations = d, assessments = list(a = list(labels = c("Short","Long"), accepted = TRUE),
      b = list(labels = c("Short","Middle","Long"), accepted = TRUE)),
    taxonomies = list(two = ssf_taxonomy("two", c("Short","Long")),
      three = ssf_taxonomy("three", c("Short","Middle","Long")))), class = "ssf_state_comparison")
  z <- compare_ssf_state_labels(x, "a", "b")
  expect_equal(unname(z$mapping), c("Short", "Short", "Long"))
  expect_equal(which(z$steps$.changed), c(4L, 5L))
  expect_equal(which(z$steps$.ambiguous), 3L)
  expect_equal(sum(z$crosswalk), 8)
  expect_error(compare_ssf_state_labels(x, "a", "a"), "distinct")
  if (requireNamespace("ggplot2", quietly = TRUE))
    expect_s3_class(plot(x, type = "track", from = "a", to = "b", step_index = 2:6), "ggplot")
})
