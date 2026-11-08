# Reconstruct origin-altitude context and one common set of control steps.
# The altitude is observed before the step, not attached to an unobserved
# three-dimensional control endpoint. It interacts with turning persistence.
condor_context_data <- function(n_controls = 512L, seed = 901205L) {
  gps <- utils::read.csv(system.file("extdata", "condor-flight-steps.csv",
    package = "behavSSF", mustWork = TRUE))
  start <- c(TRUE, gps$ID[-1L] != gps$ID[-nrow(gps)])
  gps$altitude_origin_m <- c(NA_real_, gps$altitude_m[-nrow(gps)])
  gps <- gps[!start, ]
  centre <- stats::median(gps$altitude_origin_m)
  gps$altitude_origin_100m <- (gps$altitude_origin_m - centre) / 100
  gps$stratum <- seq_len(nrow(gps))
  heading <- atan2(gps$y_m - gps$y0_m, gps$x_m - gps$x0_m)
  gps$previous_bearing <- heading - gps$angle
  gps$target_bearing <- 0 # Unused: neither model includes a target term.
  mean <- mean(gps$step)
  shape <- min(10, max(1, mean^2 / stats::var(gps$step)))
  controls <- generate_ssf_controls(gps, n_controls,
    proposal = c(shape = shape, rate = shape / mean, mean_angle = 0, concentration = 0),
    seed = seed)
  row <- rep(seq_len(nrow(gps)), each = n_controls + 1L)
  for (name in c("altitude_origin_100m", "altitude_origin_m", "altitude_m", "vertical_m_s", "source_proxy"))
    controls[[name]] <- gps[[name]][row]
  controls$x0_m <- gps$x0_m[row]
  controls$y0_m <- gps$y0_m[row]
  controls$x_m <- controls$x0_m + 1000 * controls$step * cos(gps$previous_bearing[row] + controls$angle)
  controls$y_m <- controls$y0_m + 1000 * controls$step * sin(gps$previous_bearing[row] + controls$angle)
  # Endpoint altitude, vertical speed and proxy are recorded for the used
  # step only. They are never supplied as predictors for hypothetical endpoints.
  for (name in c("altitude_m", "vertical_m_s", "source_proxy")) controls[[name]][controls$obs == 0L] <- NA
  list(choices = controls, observed = gps, altitude_centre_m = centre,
    n_controls = n_controls, control_seed = seed)
}

condor_flight_taxonomy <- function(formula = ~ step + log_step + cos_angle) {
  ssf_taxonomy("flight_modes", c("Circling-like flight", "Directed flight"),
    formula = formula, order_by = "persistence",
    rules = function(s) {
      p <- stats::setNames(s$metrics$persistence, rownames(s$metrics))
      c(persistence_gap = p[["Directed flight"]] - p[["Circling-like flight"]] > 1,
        forward_persistence = p[["Directed flight"]] > 0)
    })
}
