#' Recover movement summaries under an explicit kernel convention
#'
#' For a kernel proportional to exp(beta_step * r + beta_log * log(r)),
#' the gamma shape is beta_log + 2 for a density with respect to planar area,
#' and beta_log + 1 for a density with respect to step length and angle.
#' The rate is -beta_step. These are intrinsic, habitat-free kernel summaries;
#' they are not the marginal distribution under spatially varying habitat.
#' With `lag_log_step`, mean and standard deviation describe the baseline
#' gamma kernel at the first step of a segment. The mean also anchors the
#' autoregression for subsequent steps. Hourly step-length effects are
#' separate from these baseline summaries.
#' The reported concentration is the intrinsic persistence proxy abs(beta_cos).
#' With a target-direction effect, use [angular_parameters()] for the full
#' conditional angular concentration.
#'
#' @param coefficients A numeric coefficient matrix with rows `step`,
#'   `log_step`, and `cos_angle` and states in columns.
#' @param convention Required: `"area"` for the article's +2 convention or
#'   `"step_angle"` for a step-length/angle density. The caller must verify
#'   consistency with the measure of the sampling density w.
#' @return A data frame with state, gamma shape/rate, mean/sd step length,
#'   signed persistence, concentration, and kernel validity.
#'   Invalid gamma kernels have NA moments, never fabricated finite values.
#' @export
movement_metrics <- function(coefficients, convention) {
  convention <- match.arg(convention, c("area", "step_angle"))
  b <- as.matrix(coefficients)
  required <- c("step", "log_step", "cos_angle")
  if (!is.numeric(b) || !ncol(b) || any(!is.finite(b)) ||
      is.null(rownames(b)) || anyDuplicated(rownames(b)) || !all(required %in% rownames(b)))
    .stop("coefficients must be finite with unique step, log_step, and cos_angle rows.")
  shape <- b["log_step", ] + if (convention == "area") 2 else 1
  rate <- -b["step", ]
  valid <- shape > 0 & rate > 0
  mu <- sd <- rep(NA_real_, ncol(b))
  mu[valid] <- shape[valid] / rate[valid]
  sd[valid] <- sqrt(shape[valid]) / rate[valid]
  data.frame(state = seq_len(ncol(b)), shape = shape, rate = rate,
             mean_step = mu, sd_step = sd, persistence = b["cos_angle", ],
             concentration = abs(b["cos_angle", ]), valid = valid,
             convention = convention, row.names = NULL)
}

#' Assess a fitted ecological taxonomy
#'
#' Tests strict mean-step ordering and travelling concentration greater than
#' foraging (or encamped). For HDT, the unadjusted 95 percent Wald interval of
#' the target coefficient must be entirely positive. A significantly negative
#' target coefficient does not establish attraction. Numerical failure and
#' failed ecological rules remain distinct in the returned checks.
#'
#' @param fit An `ssf_fit` object.
#' @param taxonomy One of `HDT_F_R`, `T_F_R`, `T_E`, or `Other`, or a
#'   specification created with [ssf_taxonomy()].
#' @param convention Kernel convention passed to [movement_metrics()].
#' @param target Name of a target-direction coefficient for HDT.
#' @param min_occupancy Minimum Viterbi occupancy, between zero and one.
#' @param persistence_bounds Open lower and upper self-transition bounds.
#'   For daily transitions, the bounds apply at each of the 24 hours.
#' @param min_mean_gap Nonnegative minimum absolute gap between ordered means,
#'   in step-length units. This is not a significance test.
#' @param score_tolerance Maximum absolute objective gradient for accepting a
#'   numerically stationary fit. Missing or nonfinite scores fail this check.
#' @return A list with `accepted`, named logical `checks`, state `labels`,
#'   `metrics`, `occupancy`, and `reasons`. Wald intervals do not adjust for
#'   model selection and numerical convergence alone does not establish adequacy.
#' @export
assess_taxonomy <- function(fit, taxonomy, convention, target = "cos_diff_theta_to_target",
                            min_occupancy = 0.05, persistence_bounds = c(0.05, 0.95),
                            min_mean_gap = 0, score_tolerance = 0.01) {
  custom <- inherits(taxonomy, "ssf_taxonomy")
  if (!custom) taxonomy <- match.arg(taxonomy, names(.reference_state_counts))
  if (!inherits(fit, "ssf_fit")) .stop("fit must be an ssf_fit object.")
  measure <- fit$density_measure
  if (!is.character(measure) || length(measure) != 1L || is.na(measure) ||
      !measure %in% c("area", "step_angle"))
    .stop("Ecological assessment requires an explicit fitted density measure.")
  if (measure != convention)
    .stop("The metric convention must match the fitted density measure.")
  .scalar(min_occupancy, "min_occupancy", 0, 1)
  .scalar(min_mean_gap, "min_mean_gap", 0)
  .scalar(score_tolerance, "score_tolerance", 0)
  if (!is.numeric(persistence_bounds) || length(persistence_bounds) != 2L ||
      any(!is.finite(persistence_bounds)) || persistence_bounds[1] < 0 ||
      persistence_bounds[2] > 1 || persistence_bounds[1] >= persistence_bounds[2])
    .stop("persistence_bounds must be increasing bounds in [0, 1].")
  expected <- if (custom) length(taxonomy$labels) else .reference_state_counts[[taxonomy]]
  if (fit$n_states != expected) .stop("taxonomy and fitted state count disagree.")
  m <- movement_metrics(fit$coefficients, convention)
  k <- fit$n_states
  occ <- tabulate(fit$state, nbins = k) / length(fit$state)
  score <- fit$diagnostics$max_abs_gradient
  score_ok <- is.numeric(score) && length(score) == 1L && is.finite(score) &&
    score <= score_tolerance
  checks <- c(converged = isTRUE(fit$diagnostics$converged),
              hessian = isTRUE(fit$diagnostics$hessian_positive),
              score = score_ok,
              valid_kernel = all(m$valid))
  labels <- rep(NA_character_, k)
  if (k > 1L) {
    stay <- if (is.null(fit$transition_coefficients)) diag(fit$transition) else
      apply(fit$transition, 3L, diag)
    checks <- c(checks,
      transition = all(stay > persistence_bounds[1] & stay < persistence_bounds[2]),
      occupancy = all(occ >= min_occupancy & occ > 0))
  }
  if (custom) {
    assessment <- .check_custom_taxonomy(fit, taxonomy, m, occ, checks, min_mean_gap)
    checks <- assessment$checks
    labels <- assessment$labels
    taxonomy <- taxonomy$name
  } else if (k > 1L) {
    order_mu <- order(m$mean_step, decreasing = TRUE, na.last = TRUE)
    checks <- c(checks, mean_order = all(m$valid) &&
                  all(-diff(m$mean_step[order_mu]) > min_mean_gap))
    travelling <- order_mu[1]
    next_state <- order_mu[2]
    checks <- c(checks, direction_order =
      m$concentration[travelling] > m$concentration[next_state] &&
        m$persistence[travelling] >= 0)
    labels[order_mu] <- .reference_state_labels[[taxonomy]]
  } else labels <- .reference_state_labels[[taxonomy]]
  if (!custom && taxonomy == "HDT_F_R") {
    target_positive <- FALSE
    if (target %in% rownames(fit$coefficients) && !is.null(fit$covariance)) {
      travelling <- which(labels == "HDT")
      row <- match(target, rownames(fit$coefficients))
      index <- row + (travelling - 1L) * nrow(fit$coefficients)
      lower <- fit$coefficients[row, travelling] - stats::qnorm(.975) *
        sqrt(fit$covariance[index, index])
      target_positive <- is.finite(lower) && lower > 0
    }
    checks <- c(checks, target_positive = target_positive)
  }
  checks[is.na(checks)] <- FALSE
  list(accepted = all(checks), checks = checks, labels = labels, metrics = m,
       occupancy = occ, reasons = names(checks)[!checks], taxonomy = taxonomy)
}

#' Fit and assess the descending taxonomies of the reference workflow
#'
#' By default, tries HDT/F/R (if a target is supplied), then T/F/R, T/E, and
#' Other. Supply an ordered list to use custom taxonomies or mix them with
#' built-in candidates. The first candidate passing its checks is selected.
#' Returns every attempted candidate, including numerical failures. If no
#' candidate passes, status is `no_accepted_model`; it does not reinterpret
#' numerical failure as evidence of homogeneous behaviour.
#'
#' @param data Case-control step data.
#' @param convention Required kernel convention, see [movement_metrics()].
#' @param target Optional name of a target-direction column for the built-in
#'   HDT candidate. Custom candidates supply covariates in their formula.
#' @param n_starts,short_iter,long_iter,spread,seed Fitting controls, see [fit_ssf()].
#' @param min_occupancy,persistence_bounds,min_mean_gap,score_tolerance Assessment controls,
#'   see [assess_taxonomy()].
#' @param engine Fitting engine passed to [fit_ssf()].
#' @param normalization importance for case-control data, exact for observed
#'   ssf_steps and compatible analytically normalised kernels.
#' @param step_ar Logical. Include state-specific length autoregression in
#'   every candidate when using exact normalisation.
#' @param taxonomies NULL for the reference sequence, or a nonempty ordered list
#'   of [ssf_taxonomy()] specifications and/or built-in taxonomy names. Include
#'   `"Other"` explicitly to try the one-state fallback in a custom sequence.
#' @return A `guided_ssf` object with selection status, selected fit and
#'   assessment, and the complete candidate history. No model averaging,
#'   post-selection confidence calibration or HSMM is performed.
#' @export
fit_guided_ssf <- function(data, convention, target = NULL, n_starts = 5L,
                           short_iter = 30L, long_iter = 300L, spread = 0.5,
                           seed = NULL, min_occupancy = 0.05,
                           persistence_bounds = c(0.05, 0.95), min_mean_gap = 0,
                           engine = c("native", "R"), score_tolerance = 0.01,
                           normalization = c("importance", "exact"), step_ar = FALSE,
                           taxonomies = NULL) {
  convention <- match.arg(convention, c("area", "step_angle"))
  engine <- match.arg(engine)
  normalization <- match.arg(normalization)
  if (!is.logical(step_ar) || length(step_ar) != 1L || is.na(step_ar))
    .stop("step_ar must be TRUE or FALSE.")
  if (step_ar && normalization != "exact") .stop("step_ar requires exact normalisation.")
  if (identical(target, .step_ar_term)) .stop("The AR coefficient cannot be a target term.")
  if (!is.null(target) && (!is.character(target) || length(target) != 1L ||
                          is.na(target) || !target %in% names(data)))
    .stop("target must be NULL or the name of a column in data.")
  .scalar(min_occupancy, "min_occupancy", 0, 1)
  .scalar(min_mean_gap, "min_mean_gap", 0)
  .scalar(score_tolerance, "score_tolerance", 0)
  .scalar(n_starts, "n_starts", 1, 10000, TRUE)
  .scalar(short_iter, "short_iter", 1, 1e6, TRUE)
  .scalar(long_iter, "long_iter", 1, 1e6, TRUE)
  if (!is.numeric(persistence_bounds) || length(persistence_bounds) != 2L ||
      any(!is.finite(persistence_bounds)) || persistence_bounds[1] < 0 ||
      persistence_bounds[2] > 1 || persistence_bounds[1] >= persistence_bounds[2])
    .stop("persistence_bounds must be increasing bounds in [0, 1].")
  if (normalization == "exact") {
    if (!inherits(data,"ssf_steps")) .stop("Exact fitting requires prepare_ssf_steps().")
  } else if (!inherits(data, "ssf_data")) data <- prepare_ssf_data(data)
  data <- convert_density_measure(data, convention)
  required <- c("step", "log_step", "cos_angle")
  if (!all(required %in% names(data))) .stop("The guided workflow requires step, log_step, and cos_angle.")
  if (step_ar) required <- c(required, .step_ar_term)
  # Reject input errors before starting the candidate search.
  base_formula <- stats::reformulate(required)
  design <- if (normalization=="exact") .design_exact else .design
  design(data, base_formula)
  if (!is.null(target)) design(data, stats::reformulate(c(required, target)))
  if (is.null(taxonomies)) {
    taxonomies <- as.list(c(if (!is.null(target)) "HDT_F_R", "T_F_R", "T_E", "Other"))
  }
  if (!is.list(taxonomies) || !length(taxonomies))
    .stop("taxonomies must be NULL or a nonempty ordered list of candidate taxonomies.")
  candidates <- lapply(taxonomies, function(spec) {
    if (inherits(spec, "ssf_taxonomy")) {
      f <- if (step_ar) stats::update(spec$formula, ~ . + lag_log_step) else spec$formula
      design(data, f)
      return(list(name = spec$name, formula = f, n_states = length(spec$labels), spec = spec))
    }
    if (!is.character(spec) || length(spec) != 1L || !spec %in% names(.reference_state_counts))
      .stop("Each taxonomy must be a built-in name or a specification from ssf_taxonomy().")
    if (spec == "HDT_F_R" && is.null(target)) .stop("HDT_F_R requires a target column.")
    f <- if (spec == "HDT_F_R") stats::reformulate(c(required, target)) else base_formula
    list(name = spec, formula = f, n_states = .reference_state_counts[[spec]], spec = spec)
  })
  if (anyDuplicated(vapply(candidates, `[[`, character(1), "name")))
    .stop("Taxonomy names must be unique within a candidate sequence.")
  .with_seed(seed, function() {
    history <- list()
    chosen <- NULL
    for (candidate in candidates) {
      mode <- candidate$name
      attempt <- tryCatch({
        model <- fit_ssf(data, candidate$formula, candidate$n_states, n_starts, short_iter,
                         long_iter, spread, engine = engine, normalization = normalization)
        list(fit = model, assessment = NULL, error = NULL)
      }, error = function(e) list(fit = NULL, assessment = NULL,
                                  error = conditionMessage(e), starts = e$attempts))
      if (!is.null(attempt$fit)) {
        attempt$assessment <- assess_taxonomy(attempt$fit, candidate$spec, convention, target,
                                              min_occupancy, persistence_bounds, min_mean_gap,
                                              score_tolerance)
      }
      history[[mode]] <- attempt
      if (!is.null(attempt$assessment) && attempt$assessment$accepted) {
        chosen <- mode
        break
      }
    }
    out <- list(status = if (is.null(chosen)) "no_accepted_model" else "selected",
                taxonomy = chosen, fit = if (!is.null(chosen)) history[[chosen]]$fit else NULL,
                assessment = if (!is.null(chosen)) history[[chosen]]$assessment else NULL,
                history = history, convention = convention, seed = seed)
    class(out) <- "guided_ssf"
    out
  })
}

#' @export
print.guided_ssf <- function(x, ...) {
  if (is.null(x$fit)) cat("Guided HMM-SSF. No candidate selected.\n") else
    cat("Guided HMM-SSF with", nrow(x$fit$observations), "observed steps\n")
  .print_guided_candidates(.guided_candidates(x))
  invisible(x)
}
