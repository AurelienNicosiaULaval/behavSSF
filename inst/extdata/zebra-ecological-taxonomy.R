# An exploratory interpretation of the saved three-state habitat model.
# These definitions were suggested by this fitted example, not specified
# before examining the data. Apply them to new data before claiming confirmation.
zebra_ecological_taxonomy <- function() {
  ssf_taxonomy(
    "grassland_movement",
    c("Local movement", "Grassland movement", "Persistent grassland travel"),
    formula = ~ cos_angle + cos_diff_theta_to_target + habitat,
    order_by = "mean_step",
    rules = function(s) {
      local <- "Local movement"
      middle <- "Grassland movement"
      travel <- "Persistent grassland travel"
      habitats <- c("habitatbushed grassland", "habitatbushland", "habitatwoodland")
      c(
        separated_lengths = s$metrics[middle, "mean_step"] >
          2 * s$metrics[local, "mean_step"] &&
          s$metrics[travel, "mean_step"] > 2 * s$metrics[middle, "mean_step"],
        increasing_persistence = ssf_contrast(s, "cos_angle", local, middle)$lower > 0 &&
          ssf_contrast(s, "cos_angle", middle, travel)$lower > 0,
        middle_grassland_selection = all(ssf_contrast(s, habitats, middle)$upper < 0),
        travel_grassland_selection = all(ssf_contrast(s, habitats, travel)$upper < 0)
      )
    }
  )
}
