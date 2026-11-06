# Rscript --vanilla data-raw/prepare-deer.R /path/to/amt_0.3.1.0.tar.gz
# Source: amt (GPL-3), with GPS fixes credited to Verein für Wildtierforschung
# Dresden und Göttingen e.V. in amt::deer and forest cover credited to the JRC.
# Rebuilding this data file requires terra.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the amt 0.3.1.0 source archive.")
archive <- normalizePath(args[[1L]], mustWork = TRUE)
stopifnot(unname(tools::md5sum(archive)) == "b512f28cc913dc0858c896ee3852c6e8")
local({
  scratch <- tempfile()
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE))
  members <- c("amt/data/deer.rda", "amt/inst/external/sh_forest.rds")
  utils::untar(archive, files = members, exdir = scratch)
  fixes <- new.env()
  load(file.path(scratch, members[1L]), envir = fixes)
  gps <- as.data.frame(fixes$deer)
  stopifnot(nrow(gps) == 826L, all(stats::complete.cases(gps)))

  xy <- as.matrix(gps[c("x_", "y_")])
  delta <- diff(xy)
  interval <- diff(as.numeric(gps$t_))
  distance <- sqrt(rowSums(delta^2))
  heading <- atan2(delta[, 2L], delta[, 1L])
  regular <- diff(gps$burst_) == 0L & abs(interval - 6 * 3600) <= 10 * 60
  origin <- which(regular[-1L] & regular[-length(regular)]) + 1L
  end <- origin + 1L
  turn <- heading[origin] - heading[origin - 1L]

  # Read the cell at each original GPS endpoint before shifting the map coordinates.
  forest <- terra::unwrap(readRDS(file.path(scratch, members[2L])))
  endpoints <- terra::project(terra::vect(xy[end, , drop = FALSE], crs = "EPSG:3035"),
    terra::crs(forest))
  end_cover <- terra::extract(forest, endpoints, method = "simple", ID = FALSE)[[1L]]
  stopifnot(all(terra::res(forest) == 25), length(end_cover) == length(end),
    all(end_cover %in% 0:1))

  # Retain the source fix IDs and use a common origin for both map coordinates.
  km <- xy / 1000
  centre <- floor(colMeans(km[end, , drop = FALSE]))
  km <- sweep(km, 2L, centre)
  deer <- data.frame(ID = paste0("segment-", 1L + cumsum(!regular)[origin]),
    source_burst = gps$burst_[end], start_fix = origin, end_fix = end, time_utc = gps$t_[end],
    x0 = km[origin, 1L], y0 = km[origin, 2L], x = km[end, 1L], y = km[end, 2L],
    step = distance[origin] / 1000, angle = atan2(sin(turn), cos(turn)),
    end_in_forest = end_cover == 1)
  stopifnot(nrow(deer) == 759L, length(unique(deer$ID)) == 29L,
    max(table(deer$ID)) == 94L, all(deer$step > 0))
  saveRDS(deer, "inst/extdata/deer-restart-steps.rds", version = 2)
})
