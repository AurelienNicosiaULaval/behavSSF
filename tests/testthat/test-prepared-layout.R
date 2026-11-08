test_that("prepared choice sets cannot be silently reordered", {
  raw <- data.frame(ID = 1L, stratum = rep(1:8, each = 4),
    obs = rep(c(1, 0, 0, 0), 8), w = 1,
    step = seq(.1, 3.2, length.out = 32), cos_angle = cos(1:32))
  d <- prepare_ssf_data(raw, density_measure = "area")
  f <- ~ step + log_step + cos_angle
  b <- rbind(step = c(-2, -1), log_step = c(0, .2), cos_angle = c(.2, 1))
  expected <- ssf_log_emissions(d, f, b)
  reverse <- unlist(lapply(split(seq_len(nrow(d)), d$stratum), rev), use.names = FALSE)
  broken <- d[reverse, ]
  expect_error(ssf_log_emissions(broken, f, b), "case must be the first row")
  for (engine in c("R", "native"))
    expect_error(fit_ssf(broken, f, engine = engine), "case must be the first row")
  interleaved <- d[c(which(d$obs == 1L), which(d$obs == 0L)), ]
  expect_error(ssf_log_emissions(interleaved, f, b), "must be contiguous")
  controls <- unlist(lapply(split(seq_len(nrow(d)), d$stratum),
    function(i) c(i[1L], rev(i[-1L]))), use.names = FALSE)
  expect_equal(ssf_log_emissions(d[controls, ], f, b), expected, tolerance = 1e-12)
  expect_equal(ssf_log_emissions(raw[reverse, ], f, b), expected, tolerance = 1e-12)
  expect_equal(ssf_log_emissions(d[d$stratum %in% c(2, 5), ], f, b), expected[c(2, 5), ])
  expect_error(ssf_log_emissions(d[-1L, ], f, b), "one case")
})

test_that("custom taxonomy formulas and rule failures explain missing terms", {
  b <- rbind(step = c(-4, -1), log_step = c(0, 0), cos_angle = c(.2, 1))
  fit <- structure(list(coefficients = b, n_states = 2L, density_measure = "area",
    formula = ~ step + log_step + cos_angle, state = rep(1:2, each = 20),
    transition = matrix(c(.8, .2, .2, .8), 2), covariance = diag(.04, 6),
    diagnostics = list(converged = TRUE, hessian_positive = TRUE,
      max_abs_gradient = 0)), class = "ssf_fit")
  missing <- ssf_taxonomy("elevation", c("Short", "Long"),
    ~ step + log_step + cos_angle + elevation)
  expect_error(assess_taxonomy(fit, missing, "area"), "absent.*elevation")
  valid <- ssf_taxonomy("lengths", c("Short", "Long"))
  expect_true(assess_taxonomy(fit, valid, "area")$accepted)
  fit$formula <- ~ step + log_step + cos_angle + habitat
  expect_true(assess_taxonomy(fit, valid, "area")$accepted)
  invalid_rule <- ssf_taxonomy("bad_rule", c("Short", "Long"),
    rules = function(s) c(elevation = s$coefficients["elevation", "Long"] > 0))
  expect_error(assess_taxonomy(fit, invalid_rule, "area"),
    "bad_rule.*rule failed:.*Available coefficients:.*state labels")
  tied <- ssf_taxonomy("tied", c("Short", "Long"), order_by = function(m, b) c(1, 1))
  expect_false(assess_taxonomy(fit, tied, "area")$accepted)
})
