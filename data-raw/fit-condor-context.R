# Rscript --vanilla path/to/fit-condor.R output-directory
library(behavSSF)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, packageVersion("behavSSF") == "0.1.0.9014")
out <- args[1L]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
helper <- system.file("extdata", "condor-context.R", package = "behavSSF", mustWork = TRUE)
source(helper)
gps_path <- system.file("extdata", "condor-flight-steps.csv", package = "behavSSF", mustWork = TRUE)
gps <- utils::read.csv(gps_path)
steps <- prepare_ssf_steps(gps, "area", time = "seconds", expected_interval = 5)
exact <- fit_guided_ssf(steps, convention = "area", normalization = "exact",
  taxonomies = list(condor_flight_taxonomy(), "Other"),
  persistence_bounds = c(.05, .995), n_starts = 8L,
  short_iter = 100L, long_iter = 900L, seed = 1729L)
stopifnot(identical(exact$taxonomy, "flight_modes"))
input <- condor_context_data()
base_formula <- ~ step + log_step + cos_angle
context_formula <- ~ step + log_step + cos_angle + cos_angle:altitude_origin_100m
initial <- list(coefficients = stats::coef(exact$fit), transition = exact$fit$transition)
fit_base <- fit_ssf(input$choices, base_formula, n_states = 2L,
  initial_values = initial, n_starts = 8L, short_iter = 100L, long_iter = 900L,
  seed = 901206L)
initial$coefficients <- rbind(stats::coef(fit_base), "cos_angle:altitude_origin_100m" = 0)
initial$transition <- fit_base$transition
fit_context <- fit_ssf(input$choices, context_formula, n_states = 2L,
  initial_values = initial, n_starts = 8L, short_iter = 100L, long_iter = 900L,
  seed = 901207L)
assessment_base <- assess_taxonomy(fit_base, condor_flight_taxonomy(base_formula),
  convention = "area", persistence_bounds = c(.05, .995))
assessment_context <- assess_taxonomy(fit_context, condor_flight_taxonomy(context_formula),
  convention = "area", persistence_bounds = c(.05, .995))
print(assessment_base)
print(assessment_context)
x <- list(exact = exact, baseline = fit_base, context = fit_context,
  baseline_assessment = assessment_base, context_assessment = assessment_context,
  altitude_centre_m = input$altitude_centre_m, n_controls = input$n_controls,
  control_seed = input$control_seed, altitude_range_m = range(input$observed$altitude_origin_m),
  altitude_quantiles_m = stats::quantile(input$observed$altitude_origin_m, c(.1, .5, .9)),
  metadata = list(package = as.character(packageVersion("behavSSF")),
    gps_md5 = unname(tools::md5sum(gps_path)),
    helper_md5 = unname(tools::md5sum(helper)), script_md5 = unname(tools::md5sum(script))))
saveRDS(x, file.path(out, "condor-context-fit.rds"), compress = "xz")
write.csv(data.frame(model = c("Baseline", "Origin altitude"),
  steps = c(length(fit_base$state), length(fit_context$state)),
  log_likelihood = c(as.numeric(stats::logLik(fit_base)), as.numeric(stats::logLik(fit_context))),
  AIC = c(stats::AIC(fit_base), stats::AIC(fit_context)),
  accepted = c(assessment_base$accepted, assessment_context$accepted)),
  file.path(out, "model-comparison.csv"), row.names = FALSE)
write.csv(ssf_state_contrasts(fit_context, assessment = if (assessment_context$accepted)
  assessment_context else NULL), file.path(out, "state-contrasts.csv"), row.names = FALSE)
dir.create(file.path(out, "fitting-sources"), showWarnings = FALSE)
file.copy(script, file.path(out, "fitting-sources/fit-condor.R"), overwrite = TRUE)
file.copy(helper, file.path(out, "fitting-sources/condor-context.R"), overwrite = TRUE)
write.csv(data.frame(file = script, md5 = unname(tools::md5sum(script))),
  file.path(out, "source-manifest.csv"), row.names = FALSE)
writeLines(c(paste("Started:", started), paste("Finished:", Sys.time()),
  paste("Elapsed seconds:", as.numeric(difftime(Sys.time(), started, units = "secs"))),
  paste("Version:", packageVersion("behavSSF")), "Controls per step: 512",
  "The baseline and context models use the same 1287 observations and the same controls.",
  "This is an additional fitting pass, separate from the original four-study run."),
  file.path(out, "run-record.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "session-info.txt"))
