#' Compute observable movement summaries without crossing segment boundaries
#'
#' The function does not compare decoded states with latent simulated states.
#' Each segment includes its supplied starting point for displacement and MSD.
#' Short-step runs are measured within the available segment, including runs
#' censored by segment boundaries, identically for observed and simulated data.
#'
#' @param track Data frame with ID, x, y, step and angle columns, in order.
#' @param segments Starting metadata as supplied to [simulate_ssf()].
#' @param lags Positive integer displacement lags in numbers of steps.
#' @param short_step Positive step threshold in the same units as coordinates.
#' @return Named numeric vector. Distances use input units, MSD squared units,
#'   and short-run durations numbers of steps. Undefined correlations are NA.
#' @export
movement_diagnostics <- function(track, segments, lags = c(1L, 2L, 4L, 8L, 16L),
                                  short_step = 0.05) {
  .scalar(short_step, "short_step", 0)
  if (!is.numeric(lags) || !length(lags) || anyNA(lags) || any(lags < 1 | lags != floor(lags)))
    .stop("lags must be positive integers.")
  required <- c("ID", "x", "y", "step", "angle")
  if (!is.data.frame(track) || !nrow(track) || !all(required %in% names(track)) ||
      anyNA(track$ID) || any(!vapply(track[required[-1]], function(x)
        is.numeric(x) && all(is.finite(x)), logical(1))) || any(track$step < 0))
    .stop("track must contain finite observed movement quantities and segment IDs.")
  if (!is.data.frame(segments) || !all(c("ID", "x", "y", "n_steps") %in% names(segments)) ||
      anyDuplicated(segments$ID) || !setequal(track$ID, segments$ID) ||
      any(!vapply(segments[c("x", "y", "n_steps")], function(x)
        is.numeric(x) && all(is.finite(x)), logical(1))) ||
      any(segments$n_steps < 1 | segments$n_steps != floor(segments$n_steps)) ||
      anyDuplicated(rle(as.character(track$ID))$values))
    .stop("Segment metadata and contiguous track IDs must agree.")
  runs <- numeric(); endpoint <- excursion <- numeric(nrow(segments))
  before <- after <- numeric()
  sums <- counts <- numeric(length(lags))
  for (s in seq_len(nrow(segments))) {
    d <- track[track$ID == segments$ID[s], ]
    if (nrow(d) != segments$n_steps[s]) .stop("Segment lengths do not match track rows.")
    xy <- rbind(c(segments$x[s], segments$y[s]), as.matrix(d[c("x", "y")]))
    n <- nrow(d)
    coordinate_steps <- sqrt(rowSums((xy[-1, , drop = FALSE]-xy[-nrow(xy), , drop = FALSE])^2))
    if (max(abs(coordinate_steps-d$step)) > 1e-7*(1+max(d$step)))
      .stop("Segment coordinates do not match step lengths; check for missing observations.")
    displacement <- sweep(xy, 2, xy[1, ])
    radii <- sqrt(rowSums(displacement^2))
    endpoint[s] <- radii[n+1L]
    excursion[s] <- max(radii)
    run <- rle(d$step < short_step)
    runs <- c(runs, run$lengths[run$values])
    if (n > 1L) {
      before <- c(before, d$step[-n]); after <- c(after, d$step[-1L])
    }
    for (i in seq_along(lags)) if (n >= lags[i]) {
      h <- lags[i]
      differences <- xy[(h+1L):(n+1L), , drop = FALSE] - xy[1L:(n+1L-h), , drop = FALSE]
      sums[i] <- sums[i] + sum(differences^2)
      counts[i] <- counts[i] + nrow(differences)
    }
  }
  correlation <- if (length(before) > 1L && stats::sd(before) > 0 && stats::sd(after) > 0)
    stats::cor(before, after) else NA_real_
  msd <- ifelse(counts > 0, sums/counts, NA_real_)
  names(msd) <- paste0("msd_lag_", lags)
  c(mean_step = mean(track$step), sd_step = stats::sd(track$step),
    q95_step = unname(stats::quantile(track$step, .95)),
    short_step_fraction = mean(track$step < short_step), mean_cos_turn = mean(cos(track$angle)),
    lag1_step_correlation = correlation,
    short_run_q95 = if (length(runs)) unname(stats::quantile(runs, .95)) else 0,
    median_endpoint_displacement = stats::median(endpoint),
    q90_maximum_excursion = unname(stats::quantile(excursion, .9)), msd)
}
