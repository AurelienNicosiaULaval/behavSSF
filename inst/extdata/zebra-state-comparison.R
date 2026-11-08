# Reconstruct the common choices using the deposited vegetation grid.
# Positions and step lengths are in km. The target was fixed at the step origin.
zebra_state_data <- function(seed = 660301L, n_controls = 128L) {
  gps <- readRDS(system.file('extdata', 'hwange-zebra-steps.rds', package = 'behavSSF', mustWork = TRUE))
  gps <- gps[as.integer(table(gps$ID)[gps$ID]) > 1L, ]
  grid <- readRDS(system.file('extdata', 'hwange-zebra-habitat.rds', package = 'behavSSF', mustWork = TRUE))
  saved_seed <- if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) get('.Random.seed', .GlobalEnv) else NULL
  saved_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(saved_kind))
    if (is.null(saved_seed)) {
      if (exists('.Random.seed', .GlobalEnv, inherits = FALSE)) rm('.Random.seed', envir = .GlobalEnv)
    } else assign('.Random.seed', saved_seed, envir = .GlobalEnv)
  })
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  set.seed(seed)
  n <- nrow(gps)
  J <- as.integer(n_controls)
  stopifnot(length(J) == 1L, !is.na(J), J >= 1L, J == n_controls)
  use_gamma <- matrix(runif(n * J) < .5, n, J)
  gamma_draw <- matrix(rgamma(n * J, shape = .55, rate = 2), n, J)
  lognormal_draw <- matrix(rlnorm(n * J, meanlog = -1.8, sdlog = 1.5), n, J)
  r <- as.vector(t(cbind(gps$step, ifelse(use_gamma, gamma_draw, lognormal_draw))))
  turn <- as.vector(t(cbind(gps$angle, matrix(runif(n * J, -pi, pi), n, J))))
  row <- rep(seq_len(n), each = J + 1L)
  heading <- gps$previous_bearing[row] + turn
  x <- gps$x0[row] + r * cos(heading)
  y <- gps$y0[row] + r * sin(heading)
  column <- floor((x - grid$extent[1L]) / grid$resolution[1L]) + 1L
  raster_row <- floor((grid$extent[4L] - y) / grid$resolution[2L]) + 1L
  inside <- column >= 1L & column <= grid$ncol & raster_row >= 1L & raster_row <= grid$nrow
  code <- rep(NA_integer_, length(r))
  code[inside] <- grid$values[(raster_row[inside] - 1L) * grid$ncol + column[inside]]
  observed <- rep(seq.int(0L, J), n) == 0L
  stopifnot(all(!is.na(code[observed])),
    max(abs(x[observed] - gps$x)) < 1e-8, max(abs(y[observed] - gps$y)) < 1e-8)
  # After rejecting points outside the map, q is proportional to the
  # conditioned proposal. Its common per-step factor cancels from choices;
  # this convention does not estimate the absolute physical normaliser.
  q <- (.5 * dgamma(r, .55, rate = 2) + .5 * dlnorm(r, -1.8, 1.5)) / (2 * pi * r)
  choices <- data.frame(ID = gps$ID[row], stratum = row, obs = as.integer(observed),
    w = q, step = r, log_step = log(r), cos_angle = cos(turn),
    cos_diff_theta_to_target = cos(heading - gps$target_bearing[row]),
    habitat = relevel(factor(grid$labels[code]), ref = 'grassland'),
    angle = turn, previous_bearing = gps$previous_bearing[row],
    target_bearing = gps$target_bearing[row], target_active = gps$target_active[row],
    x0 = gps$x0[row], y0 = gps$y0[row], x = x, y = y, time = gps$time[row])
  behavSSF::prepare_ssf_data(choices[!is.na(code), ], density_measure = 'area')
}

# Two informative starts, calculated separately inside every training split.
# The movement anchor uses no endpoint habitat terms. The second start uses
# equal-sized groups of training step lengths. Neither sees the reserved block.
zebra_state_initializers <- function(data, taxonomy, step_family, seed, block) {
  stopifnot(step_family %in% c('gamma', 'weibull'))
  gps <- data[data$obs == 1L, c('ID', 'time', 'step', 'angle', 'previous_bearing',
    'target_bearing', 'target_active', 'x0', 'y0', 'x', 'y')]
  k <- length(taxonomy$labels)
  anchor <- tryCatch(behavSSF::fit_ssf(behavSSF::prepare_ssf_steps(gps, 'area'),
    ~ step + log_step + cos_angle + cos_diff_theta_to_target, n_states = k,
    n_starts = 5L, short_iter = 65L, long_iter = 350L, normalization = 'exact',
    seed = 61001L + 100L * block + k), error = identity)
  groups <- pmin(k, ceiling(rank(gps$step, ties.method = 'first') * k / nrow(gps)))
  b <- vapply(seq_len(k), function(i) {
    d <- gps[groups == i, ]
    mu <- mean(d$step)
    shape <- min(12, max(.4, mu^2 / var(d$step)))
    c(step = -shape / mu, log_step = shape - 2,
      cos_angle = max(-1, min(2, 2 * mean(cos(d$angle)))), cos_diff_theta_to_target = .1)
  }, numeric(4))
  transition <- matrix(.32 / max(1L, k - 1L), k, k)
  diag(transition) <- if (k == 1L) 1 else .68
  sources <- list(Length_groups = list(coefficients = b, transition = transition))
  if (!inherits(anchor, 'error')) sources <- c(list(Movement_fit = anchor), sources)
  if (k == 1L) sources <- sources[1L]
  lapply(sources, function(a) {
    b <- a$coefficients
    shape <- b['log_step', ] + 2
    rate <- -b['step', ]
    radial <- if (step_family == 'gamma') rbind(step_shape = shape, step_rate = rate) else {
      shape_w <- vapply(1 / shape, function(cv2) uniroot(function(v)
        expm1(lgamma(1 + 2 / v) - 2 * lgamma(1 + 1 / v)) - cv2, c(.10, 100))$root, numeric(1))
      rbind(step_shape = shape_w, step_scale = (shape / rate) / exp(lgamma(1 + 1 / shape_w)))
    }
    habitat_names <- grep('^habitat', colnames(model.matrix(taxonomy$formula, data[1L, ])), value = TRUE)
    coefficients <- rbind(radial, b[c('cos_angle', 'cos_diff_theta_to_target'), , drop = FALSE],
      matrix(0, length(habitat_names), k, dimnames = list(habitat_names, NULL)))
    list(coefficients = coefficients, transition = a$transition)
  })
}

# Independently polish additional K=4 starts: eight full-track starts or four starts per held-out split.
# All anchors and perturbations are constructed from training rows only.
zebra_state_initializers_extended <- function(data,taxonomy,step_family,seed,block) {
  sources <- zebra_state_initializers(data,taxonomy,step_family,seed,block)
  if (length(taxonomy$labels)!=4L) return(sources)
  set.seed(seed)
  base <- sources
  budget <- if (block == 0L) 8L else 4L
  for (i in seq_len(budget-length(sources))) {
    a <- base[[1L+(i-1L)%%length(base)]]
    radial <- grepl("^step_",rownames(a$coefficients))
    a$coefficients[radial,] <- a$coefficients[radial,]*exp(rnorm(sum(radial)*4L,0,.25))
    a$coefficients[!radial,] <- a$coefficients[!radial,]+rnorm(sum(!radial)*4L,0,.35)
    a$transition <- a$transition*exp(matrix(rnorm(16,0,.3),4))
    a$transition <- a$transition/rowSums(a$transition)
    sources[[paste0("Perturbed_",i)]] <- a
  }
  sources
}
