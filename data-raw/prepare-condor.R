# Rscript --vanilla data-raw/prepare-condor.R Vultur_gryphus.csv inst/extdata/condor-flight-steps.csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the Zenodo CSV member and an output CSV path.")
stopifnot(unname(tools::md5sum(args[1L])) == "7001d0b11bb71f9371875c64ce7b8d98")
raw <- read.csv(args[1L], na.strings = "", check.names = FALSE)
raw <- raw[which(raw$individuals == "Kuzco"), ]

nearest_index <- function(grid, times) {
  before <- pmax(1L, findInterval(grid, times))
  after <- pmin(before + 1L, length(times))
  ifelse(abs(times[before] - grid) <= abs(times[after] - grid), before, after)
}

prepare_day <- function(d) {
  d <- d[order(d$date_in_unix), ]
  seconds <- d$date_in_unix
  stopifnot(!anyDuplicated(seconds), min(seconds) == 0)
  # Retain observed fixes near the 5 s clock and remove low-speed intervals.
  grid <- seq(0, max(seconds), by = 5)
  fix <- nearest_index(grid, seconds)
  original <- d$to_impute == "False" & d$to_discard == "False"
  use <- which(abs(seconds[fix] - grid) <= .15 & original[fix])
  pos <- fix[use]
  xy <- as.matrix(d[pos, c("CoordinatesX(m)", "CoordinatesY(m)")])
  delta <- diff(xy)
  heading <- atan2(delta[, 2L], delta[, 1L])
  distance <- sqrt(rowSums(delta^2))
  good <- diff(use) == 1L & is.finite(distance) & distance > 0 &
    d$horizontal_speed[pos][-length(pos)] >= 3 & d$horizontal_speed[pos][-1L] >= 3
  good[is.na(good)] <- FALSE
  step <- which(good[-1L] & good[-length(good)]) + 1L
  turn <- heading[step] - heading[step - 1L]
  id <- paste(d$only_date[1L], cumsum(!good)[step] + 1L, sep = "/")
  from <- pos[step]
  to <- pos[step + 1L]
  midpoint <- nearest_index((seconds[from] + seconds[to]) / 2, seconds)
  # Apply the source's pointwise curvature and vertical-speed cutoffs.
  proxy <- rep(NA_character_, length(step))
  k <- abs(d$K[midpoint])
  vz <- d$vertical_speed[midpoint]
  proxy[which(k > .01 & vz > 0)] <- "Thermalling proxy"
  proxy[which(k < .01 & vz < 0)] <- "Gliding proxy"
  out <- data.frame(ID = id, day = d$only_date[1L], seconds = grid[use[step + 1L]],
    x0_m = xy[step, 1L], y0_m = xy[step, 2L], x_m = xy[step + 1L, 1L],
    y_m = xy[step + 1L, 2L], altitude_m = d$ele[to],
    step = distance[step] / 1000, angle = atan2(sin(turn), cos(turn)),
    vertical_m_s = (d$ele[to] - d$ele[from]) / (seconds[to] - seconds[from]),
    source_proxy = proxy)
  out[ave(out$seconds, out$ID, FUN = length) >= 12L, ]
}
steps <- do.call(rbind, lapply(split(raw, raw$only_date), prepare_day))
stopifnot(nrow(steps) > 500, !anyDuplicated(steps[c("day", "seconds")]),
  all(vapply(split(steps$seconds, steps$ID), function(t) all(diff(t) == 5), logical(1))))
write.csv(steps, args[2L], row.names = FALSE)
