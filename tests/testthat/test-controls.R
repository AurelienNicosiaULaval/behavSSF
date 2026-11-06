test_that("control generation uses area densities and preserves observed steps", {
  observed <- data.frame(ID = c("a", "b"), stratum = 1:2, step = c(.2, .5),
    angle = c(.1, -.3), previous_bearing = c(0, 1), target_bearing = c(.4, -.7))
  proposal <- c(shape = 2, rate = 3, mean_angle = .2, concentration = .8)
  set.seed(151); before <- .Random.seed
  pool <- generate_ssf_controls(observed, 20, proposal, seed = 197)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(pool), 42L)
  expect_equal(pool$step[pool$obs == 1], observed$step)
  expect_equal(pool$angle[pool$obs == 1], observed$angle)
  expected <- dgamma(pool$step, 2, rate = 3) / pool$step *
    exp(.8*cos(pool$angle-.2)) / (2*pi*besselI(.8, 0))
  expect_equal(pool$w, expected, tolerance = 1e-12)
  expect_identical(attr(pool, "density_measure"), "area")
  expect_identical(pool, generate_ssf_controls(observed, 20, proposal, seed = 197))
  nested <- pool[pool$alternative <= 5, ]
  expect_equal(as.integer(table(nested$stratum)), c(6L, 6L))
})

test_that("proposal samples have the specified gamma and circular moments", {
  observed <- data.frame(ID = 1, stratum = 1, step = .2, angle = .1,
    previous_bearing = 0, target_bearing = 0)
  proposal <- c(shape = 2, rate = 3, mean_angle = .2, concentration = .8)
  pool <- generate_ssf_controls(observed, 10000, proposal, seed = 198)
  controls <- pool[pool$obs == 0, ]
  expect_lt(abs(mean(controls$step) - 2/3), 6*sqrt(2)/3/sqrt(10000))
  first <- besselI(.8, 1) / besselI(.8, 0)
  second <- besselI(.8, 2) / besselI(.8, 0)
  se <- sqrt(((1+second)/2-first^2)/10000)
  expect_lt(abs(mean(cos(controls$angle-.2))-first), 6*se)
})
