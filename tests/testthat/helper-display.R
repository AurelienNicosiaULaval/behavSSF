display_fit <- function() {
  set.seed(62)
  n <- 41L
  step <- rgamma(n, shape = 2, rate = 3)
  data <- data.frame(ID = rep(c("A", "B", "C"), c(20, 1, 20)),
    step = step, angle = runif(n, -pi, pi),
    x0 = c(0, cumsum(step)[-n]), y0 = 0, x = cumsum(step), y = 0)
  steps <- prepare_ssf_steps(data, "area")
  fit_guided_ssf(steps, convention = "area", normalization = "exact",
    taxonomies = list(ssf_taxonomy("unmet_rule", "Moving", rules = function(x)
      c(custom_check = FALSE)), "Other"), seed = 63)
}

