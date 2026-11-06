# Rscript --vanilla elephant-lag-window.R /path/to/elephant.csv results [workers]
# Source: Wall et al. (2014), https://doi.org/10.5441/001/1.f321pf80
library(behavSSF)

args <- commandArgs(TRUE)
stopifnot(length(args) %in% 2:3, as.character(packageVersion("behavSSF")) == "0.1.0.9003")
stopifnot(unname(tools::md5sum(args[1L])) == "31bd957da4e5c96a9b788411d009e798")
out <- args[2L]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
raw <- utils::read.csv(args[1L], check.names = FALSE)
# Use the first track, as in the momentuHMM example for this archive.
animal <- raw[["individual-local-identifier"]][1L]
gps <- raw[raw[["individual-local-identifier"]] == animal, ]
time <- as.POSIXct(gps$timestamp, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
order <- order(time)
gps <- gps[order, ]
time <- time[order]
xy <- as.matrix(gps[c("utm-easting", "utm-northing")])
stopifnot(animal == "Salif Keita", nrow(gps) == 21972L, !anyNA(time),
  !anyDuplicated(time), all(is.finite(xy)), all(gps[["utm-zone"]] == "30N"),
  all(gps[["sensor-type"]] == "gps"))
# Relative coordinates in km, with a common origin for all plotted steps.
xy <- sweep(xy, 2L, xy[1L, ]) / 1000
dxy <- diff(xy)
step <- sqrt(rowSums(dxy^2))
heading <- atan2(dxy[, 2L], dxy[, 1L])
interval <- diff(as.numeric(time))
valid <- interval == 3600 & step > 0
i <- which(valid & c(FALSE, valid[-length(valid)]))
turn <- heading[i] - heading[i - 1L]
observed <- data.frame(ID = paste0("segment-", 1L + cumsum(!valid)[i]),
  time = time[i + 1L], x0 = xy[i, 1L], y0 = xy[i, 2L], x = xy[i+1L, 1L], y = xy[i+1L, 2L],
  step = step[i], angle = atan2(sin(turn), cos(turn)), previous_bearing = heading[i - 1L])
steps <- prepare_ssf_steps(observed, density_measure = "area",
  time = "time", expected_interval = 3600)
segments <- rle(steps$ID)
checks <- data.frame(animal = animal, input_fixes = nrow(gps), fitted_steps = nrow(steps),
  gaps = sum(interval != 3600), zero_displacements = sum(step == 0),
  segments = length(segments$values), longest_steps = max(segments$lengths),
  first_fix = format(min(time), tz = "UTC", usetz = TRUE),
  last_fix = format(max(time), tz = "UTC", usetz = TRUE))
print(checks)


workers <- if (length(args) == 3L) as.integer(args[3L]) else 1L
stopifnot(length(workers) == 1L, !is.na(workers), workers >= 1L)
origin_hour <- function(time) as.integer(format(time - 3600, "%H", tz = "UTC"))
cut <- floor(.75*nrow(steps))
train <- steps[seq_len(cut), ]
test <- steps[seq.int(cut + 1L, nrow(steps)), ]
train_hour <- origin_hour(train$time)
test_hour <- origin_hour(test$time)
stopifnot(nrow(train) == 16434L, nrow(test) == 5478L)

formula <- ~ step + log_step + cos_angle + lag_log_step
taxonomy <- ssf_taxonomy("step_lengths", c("Short steps", "Long steps"),
  rules = function(s) c(length_ratio =
    s$metrics["Long steps", "mean_step"] > 2 * s$metrics["Short steps", "mean_step"]))
assessment <- list(persistence_bounds = c(.05, .995))
message("AR(1) only on the first three quarters")
baseline <- diagnose_ssf_restarts(train, taxonomy, n_starts = 12L, seed = 37103L,
  formula = formula, assessment_args = assessment)
print(baseline)
stopifnot(nrow(baseline$solutions) > 0L)
base <- baseline$fits[[1L]]
stopifnot(abs(base$log_likelihood - (-4928.94959172661)) < 1e-6)

search_hourly <- function(warm, step_harmonics, seed, fit_formula = formula, reference_start = NULL) {
  initials <- c(warm, rep(list(NULL), 6L))
  fit_one <- function(i) {
    fit <- tryCatch(fit_ssf(train, fit_formula, n_states = 2L, n_starts = 1L,
      initial_values = initials[[i]], normalization = "exact", transition_hour = train_hour,
      transition_harmonics = 2L, step_harmonics = step_harmonics,
      short_iter = 100L, long_iter = 900L, seed = seed + i), error = identity)
    group <- if (i <= length(warm)) "Previous solution" else "Random"
    if (inherits(fit, "error"))
      return(list(log = data.frame(start = i, group = group, log_likelihood = NA_real_,
        score = NA_real_, accepted = FALSE, reason = conditionMessage(fit))))
    a <- do.call(assess_taxonomy, c(list(fit, taxonomy, "area"), assessment))
    objective_agrees <- fit$diagnostics$objective_agreement <= 1e-7
    list(fit = fit, labels = a$labels, log = data.frame(start = i, group = group,
      log_likelihood = fit$log_likelihood, score = fit$diagnostics$max_abs_gradient,
      accepted = a$accepted && objective_agrees,
      reason = paste(c(a$reasons, if (!objective_agrees) "objective_agreement"), collapse = ", ")))
  }
  fits <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_along(initials), fit_one, mc.cores = workers, mc.set.seed = FALSE)
  } else lapply(seq_along(initials), fit_one)
  stopifnot(!any(vapply(fits, inherits, logical(1), "try-error")))
  log <- do.call(rbind, lapply(fits, `[[`, "log"))
  print(log)
  admissible <- which(log$accepted)
  stopifnot(length(admissible) > 0L)
  ranked <- admissible[order(log$log_likelihood[admissible], decreasing = TRUE)]
  distinct <- ranked[!duplicated(round(log$log_likelihood[ranked], 5L))]
  if (!is.null(reference_start)) {
    # Keep the same numerical representative as the preceding benchmark.
    stopifnot(reference_start %in% admissible,
      abs(log$log_likelihood[distinct[1L]] - log$log_likelihood[reference_start]) < 1e-7)
    distinct[1L] <- reference_start
  }
  list(log = log, fits = fits[distinct])
}

message("Hourly transitions only")
warm <- lapply(baseline$fits, function(fit)
  list(coefficients = stats::coef(fit), transition = fit$transition))
transitions <- search_hourly(warm, step_harmonics = 0L, seed = 41200L)
base <- transitions$fits[[1L]]$fit
stopifnot(abs(base$log_likelihood - (-4663.987008539)) < 1e-6)

message("Hourly transitions and hourly step length")
warm <- lapply(transitions$fits, function(entry)
  list(coefficients = stats::coef(entry$fit), transition_coefficients = entry$fit$transition_coefficients))
length_search <- search_hourly(warm, step_harmonics = 2L, seed = 42300L)
fit <- length_search$fits[[1L]]$fit
stopifnot(fit$log_likelihood > base$log_likelihood)

base <- fit
message("Hourly model with lags of one and 24 hours")
lag_formula <- stats::update(formula, ~ . + lag24_log_step)
warm <- lapply(length_search$fits, function(entry) {
  f <- entry$fit
  list(coefficients = rbind(stats::coef(f), lag24_log_step = c(0, 0)),
    transition_coefficients = f$transition_coefficients,
    step_hour_coefficients = f$step_hour_coefficients)
})
lag_search <- search_hourly(warm, step_harmonics = 2L, seed = 43600L,
  fit_formula = lag_formula, reference_start = 7L)
fit <- lag_search$fits[[1L]]$fit
stopifnot(fit$log_likelihood > base$log_likelihood)

base <- fit
message("Hourly model with AR(1) and daily lags 22 through 26")
window_formula <- stats::update(lag_formula,
  ~ . + lag22_log_step + lag23_log_step + lag25_log_step + lag26_log_step)
warm <- lapply(lag_search$fits, function(entry) {
  f <- entry$fit
  zeros <- matrix(0, 4L, 2L, dimnames = list(paste0("lag", c(22, 23, 25, 26), "_log_step"), NULL))
  list(coefficients = rbind(stats::coef(f), zeros),
    transition_coefficients = f$transition_coefficients,
    step_hour_coefficients = f$step_hour_coefficients)
})
window_search <- search_hourly(warm, step_harmonics = 2L, seed = 44900L,
  fit_formula = window_formula)
fit <- window_search$fits[[1L]]$fit
stopifnot(fit$log_likelihood > base$log_likelihood)

test_decode <- function(fit) decode_ssf_diel(test, fit, test_hour)
base_test <- test_decode(base)
window_test <- test_decode(fit)
stopifnot(abs(base_test$log_likelihood/nrow(test) - .06484210) < 1e-8)
fits <- list(current = base, window = fit)
fit_table <- data.frame(model = names(fits), training_steps = nrow(train),
  test_steps = nrow(test), train_ll = c(base$log_likelihood, fit$log_likelihood),
  test_ll = c(base_test$log_likelihood, window_test$log_likelihood),
  df = vapply(fits, function(x) attr(stats::logLik(x), "df"), numeric(1)),
  n_starts = c(nrow(lag_search$log), nrow(window_search$log)),
  n_accepted = c(sum(lag_search$log$accepted), sum(window_search$log$accepted)),
  n_solutions = c(length(lag_search$fits), length(window_search$fits)),
  best_start = c(lag_search$fits[[1L]]$log$start, window_search$fits[[1L]]$log$start))
fit_table$test_score_per_step <- fit_table$test_ll / nrow(test)
print(fit_table)

coefficients <- do.call(rbind, lapply(names(fits), function(model) {
  f <- fits[[model]]
  b <- stats::coef(f)
  state <- order(movement_metrics(b, "area")$mean_step)
  lag <- c(1L, 22:26)
  terms <- c("lag_log_step", paste0("lag", 22:26, "_log_step"))
  do.call(rbind, lapply(seq_along(state), function(i) {
    rho <- numeric(length(lag))
    included <- terms %in% rownames(b)
    rho[included] <- b[terms[included], state[i]]
    data.frame(model = model, state = taxonomy$labels[i], lag = lag,
      coefficient = rho, multiplier = 2^rho)
  }))
}))
print(coefficients)

hours <- 0:23
starts <- which(!duplicated(test$ID))
segments <- data.frame(ID = test$ID[starts], x = test$x0[starts], y = test$y0[starts],
  bearing = test$previous_bearing[starts], n_steps = rle(test$ID)$lengths, hour = test_hour[starts])
ends <- cumsum(segments$n_steps)
lags <- 1:48
pairs <- lapply(lags, function(h) unlist(Map(function(i, j) {
  if (j-i < h) integer() else seq.int(i, j-h)
}, starts, ends)))
stopifnot(all(vapply(lags, function(h) {
  i <- pairs[[h]]
  all(test$ID[i] == test$ID[i+h]) &&
    all(as.numeric(difftime(test$time[i+h], test$time[i], units = "hours")) == h)
}, logical(1))))
by_hour <- split(seq_along(test_hour), factor(test_hour, levels = hours))

# Width of the local maximum in 22--26 h, measured above the higher shoulder.
# The nearest half-height crossings are linearly interpolated on 18--30 h.
daily_peak <- function(acf) {
  apex <- 21L + which.max(acf[22:26])
  baseline <- max(acf[c(18L, 30L)])
  height <- acf[apex]
  if (height <= baseline)
    return(c(peak_hour = apex, peak_height = height, peak_width = NA_real_))
  half <- (height + baseline)/2
  left <- max(which(acf[18:(apex-1L)] <= half)) + 17L
  right <- min(which(acf[(apex+1L):30] <= half)) + apex
  crossing <- function(i) i + (half-acf[i])/(acf[i+1L]-acf[i])
  c(peak_hour = apex, peak_height = height, peak_width = crossing(right-1L)-crossing(left))
}

movement <- function(track) {
  acf <- vapply(lags, function(h) {
    i <- pairs[[h]]
    stats::cor(track$step[i], track$step[i+h])
  }, numeric(1))
  hourly <- vapply(by_hour, function(i) mean(track$step[i]), numeric(1))
  adjusted <- track$step - hourly[test_hour + 1L]
  i <- pairs[[24L]]
  c(stats::setNames(acf, paste0("acf_", lags)),
    stats::setNames(hourly, paste0("mean_", hours)),
    daily_peak(acf), around_peak = mean(acf[c(22:23, 25:26)]),
    peak_contrast = acf[24] - mean(acf[c(18, 30)]),
    centered_acf_24 = stats::cor(adjusted[i], adjusted[i+24L]))
}
obs <- movement(test)
n_sims <- 199L
batches <- split(seq_len(n_sims), ceiling(seq_len(n_sims)/10L))
simulate_metrics <- function(fit, model) {
  message("Forward simulations of the held-out segments: ", model)
  do.call(rbind, lapply(seq_along(batches), function(batch) {
    sim <- simulate_ssf(fit, segments, n_sims = length(batches[[batch]]),
      seed = 38105L + batch - 1L)
    values <- t(vapply(sim$tracks, movement, obs))
    data.frame(model = model, replicate = batches[[batch]], values, row.names = NULL)
  }))
}
draws <- do.call(rbind, Map(simulate_metrics, fits, names(fits)))
stopifnot(all(is.finite(as.matrix(draws[setdiff(names(obs), "peak_width")]))), nrow(draws) == 2L*n_sims)
bands <- do.call(rbind, lapply(names(fits), function(model) {
  values <- draws[draws$model == model, names(obs)]
  q <- vapply(values, stats::quantile, numeric(3), probs = c(.025, .5, .975), names = FALSE, na.rm = TRUE)
  data.frame(model = model, metric = names(obs), observed = unname(obs),
    lower = q[1L, ], median = q[2L, ], upper = q[3L, ],
    n_defined = colSums(is.finite(as.matrix(values))), row.names = NULL)
}))
print(subset(bands, metric %in% c("acf_1", "acf_24", "acf_48", "peak_contrast", "centered_acf_24", "peak_width", "peak_height", "around_peak")))
if (length(window_search$fits) > 1L) {
  alternate <- window_search$fits[[2L]]
  alt_test <- test_decode(alternate$fit)
  other_draws <- simulate_metrics(alternate$fit, "window alternative")
  best_labels <- window_search$fits[[1L]]$labels[window_test$state]
  other_labels <- alternate$labels[alt_test$state]
  sensitivity <- data.frame(start = alternate$log$start,
    train_ll_gap = fit$log_likelihood - alternate$fit$log_likelihood,
    test_score_per_step = alt_test$log_likelihood/nrow(test),
    changed_test_labels = sum(best_labels != other_labels),
    median_acf_24 = stats::median(other_draws$acf_24),
    median_peak_width = stats::median(other_draws$peak_width, na.rm = TRUE))
  utils::write.csv(sensitivity, file.path(out, "window-alternative.csv"), row.names = FALSE)
  utils::write.csv(other_draws, file.path(out, "window-alternative-draws.csv"), row.names = FALSE)
  print(sensitivity)
}

restarts <- rbind(data.frame(model = "current", lag_search$log),
  data.frame(model = "window", window_search$log))
tables <- list(coverage = checks, fits = fit_table, restarts = restarts,
  "baseline-restarts" = baseline$restarts, "lag-coefficients" = coefficients,
  draws = draws, intervals = bands, pairs = data.frame(hours = lags, count = lengths(pairs)))
for (name in names(tables))
  utils::write.csv(tables[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
