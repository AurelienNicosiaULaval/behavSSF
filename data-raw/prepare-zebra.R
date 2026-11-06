# Run from the repository root with Codes.zip from Figshare version 2.
# Rscript data-raw/prepare-zebra.R /path/to/Codes.zip
library(RANN)
library(terra)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L,
          unname(tools::md5sum(args[1L])) == "9c5681ed705776999891bb53b62424db")

extract_dir <- tempfile("hwange-")
dir.create(extract_dir)
source_dir <- file.path("Codes", "Application", "Zebra")
files <- c("Zebra_Case_Control_review.RData", "vegetation2.grd", "vegetation2.gri")
utils::unzip(args[1L], files = file.path(source_dir, files), exdir = extract_dir)
source_dir <- file.path(extract_dir, source_dir)

choices <- as.data.frame(readRDS(file.path(source_dir, files[1L])))
# Retain the established 6,377-step cohort, then split any remaining gaps.
archive_segment <- cumsum(c(TRUE, !diff(choices$counter) %in% c(0, 1)))
keep <- ave(archive_segment, archive_segment, FUN = length) >= 21L * 20L &
  choices$obs == 1L
cases <- choices[keep, ]
n <- nrow(cases)
separated <- abs(cases$prev_x[-1L] - cases$x[-n]) > 1e-8 |
  abs(cases$prev_y[-1L] - cases$y[-n]) > 1e-8
boundary <- c(TRUE, diff(archive_segment[keep]) != 0L |
  diff(as.numeric(cases$time)) != 1800 | separated)

vegetation <- terra::rast(file.path(source_dir, files[2L]))
grass <- terra::xyFromCell(vegetation, which(terra::values(vegetation, mat = FALSE) == 1L))
origins <- as.matrix(cases[c("prev_x", "prev_y")])
nearest <- RANN::nn2(grass, origins, k = 1L)
offset <- grass[nearest$nn.idx[, 1L], , drop = FALSE] - origins

# Use relative coordinates for the map and elapsed seconds for the time check.
zebra <- data.frame(
  ID = paste0("segment_", cumsum(boundary)),
  time = as.numeric(cases$time) - as.numeric(cases$time[1L]),
  step = cases$step, angle = cases$angle,
  previous_bearing = cases$theta - cases$angle,
  target_bearing = atan2(offset[, 2L], offset[, 1L]),
  target_active = nearest$nn.dists[, 1L] > 1e-12,
  x = cases$x - cases$prev_x[1L], y = cases$y - cases$prev_y[1L],
  x0 = cases$prev_x - cases$prev_x[1L], y0 = cases$prev_y - cases$prev_y[1L]
)
stopifnot(nrow(zebra) == 6377L, sum(boundary) == 28L,
          max(abs(sqrt((zebra$x - zebra$x0)^2 + (zebra$y - zebra$y0)^2) -
                    zebra$step)) < 1e-10)
saveRDS(zebra, "inst/extdata/hwange-zebra-steps.rds", version = 2)
unlink(extract_dir, recursive = TRUE)
