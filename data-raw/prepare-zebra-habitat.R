# Run from the repository root, using the deposited Codes.zip (Figshare v2).
# Rscript --vanilla data-raw/prepare-zebra-habitat.R /path/to/Codes.zip
library(terra)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L,
  unname(tools::md5sum(args[1L])) == '9c5681ed705776999891bb53b62424db')
folder <- tempfile('zebra-habitat-')
dir.create(folder)
files <- c('vegetation2.grd', 'vegetation2.gri', 'Zebra_Case_Control_review.RData')
utils::unzip(args[1L], files = file.path('Codes/Application/Zebra', files), exdir = folder)
source_dir <- file.path(folder, 'Codes/Application/Zebra')
raster <- terra::rast(file.path(source_dir, files[1L]))
raw <- readRDS(file.path(source_dir, files[3L]))
grid <- list(values = as.integer(terra::values(raster, mat = FALSE)),
  nrow = nrow(raster), ncol = ncol(raster),
  extent = as.vector(terra::ext(raster)) -
    c(raw$prev_x[1L], raw$prev_x[1L], raw$prev_y[1L], raw$prev_y[1L]),
  resolution = terra::res(raster),
  labels = c('grassland', 'bushed grassland', 'bushland', 'woodland'))
saveRDS(grid, 'inst/extdata/hwange-zebra-habitat.rds', compress = 'xz')
unlink(folder, recursive = TRUE)
