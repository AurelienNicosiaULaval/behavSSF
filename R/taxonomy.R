#' Define a configurable state taxonomy
#'
#' Labels are assigned in increasing order of a fitted state signature.
#' A rule can then check separate movement, selection, or transition contrasts.
#' The same numerical and occupancy checks as in the guided workflow apply.
#'
#' @param name Unique identifier for the candidate. The four built-in names
#'   are reserved.
#' @param labels One to six unique state labels, in increasing signature order.
#' @param formula One-sided SSF formula. The default uses the canonical gamma
#'   terms, or only `cos_angle` if the fit specifies a length family. Supply
#'   angular and selection terms for a custom step-family formula.
#'   Assessment requires these terms in the fitted model; the fitted model
#'   may contain additional terms.
#' @param order_by Name of a numeric movement metric from [movement_metrics()]
#'   (`mean_step`, `sd_step`, `persistence`, `concentration`, `shape`, or
#'   `rate`), or a function taking the unlabelled metrics and coefficient
#'   matrix and returning one score per state. Tied or non-finite scores
#'   fail the label assignment.
#' @param rules Optional function of a single list containing `metrics`
#'   (rows named by label), `coefficients` and `standard_errors` (columns
#'   named by label), `transition` (rows and columns named by label), and
#'   `occupancy` (named by label). Its return value is a nonempty logical
#'   vector with unique names. Each FALSE or NA rejects the candidate.
#'   Standard errors are NA if the fit has no covariance matrix.
#'   The list also contains the joint natural-scale coefficient `covariance`
#'   and a `parameter_index` table identifying each row by `state` and `term`.
#'   Use [ssf_contrast()] to form a contrast without numeric parameter indices.
#' @return An `ssf_taxonomy` specification for [assess_taxonomy()] and
#'   [fit_guided_ssf()].
#' @examples
#' ssf_taxonomy("step_regimes", c("Short steps", "Long steps"),
#'   rules = function(s) {
#'     c(length_ratio = s$metrics["Long steps", "mean_step"] >
#'         2 * s$metrics["Short steps", "mean_step"])
#'   })
#' @export
ssf_taxonomy <- function(name, labels, formula = ~ step + log_step + cos_angle,
                         order_by = "mean_step", rules = NULL) {
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !grepl("^[A-Za-z][A-Za-z0-9_]*$", name) || name %in% names(.reference_state_counts))
    .stop("A custom taxonomy name must be an identifier distinct from the built-in names.")
  if (!is.character(labels) || !length(labels) || length(labels) > 6L ||
      anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(trimws(labels)))
    .stop("labels must contain one to six distinct nonempty strings.")
  if (!inherits(formula, "formula") || length(formula) != 2L)
    .stop("formula must be one-sided.")
  metrics <- c("mean_step", "sd_step", "persistence", "concentration", "shape", "rate")
  if (!is.function(order_by) && !(is.character(order_by) && length(order_by) == 1L &&
      order_by %in% metrics))
    .stop("order_by must be a numeric movement metric name or a function.")
  if (!is.null(rules) && !is.function(rules)) .stop("rules must be NULL or a function.")
  structure(list(name = name, labels = trimws(labels), formula = formula,
                 default_formula = missing(formula), order_by = order_by, rules = rules), class = "ssf_taxonomy")
}

.taxonomy_formula <- function(taxonomy, family) {
  if (!is.null(family) && isTRUE(taxonomy$default_formula)) ~ cos_angle else taxonomy$formula
}

.check_taxonomy_formula <- function(fit, taxonomy) {
  formula <- .taxonomy_formula(taxonomy, fit$step_family)
  required <- attr(stats::terms(formula), "term.labels")
  fitted <- if (inherits(fit$formula, "formula"))
    attr(stats::terms(fit$formula), "term.labels") else rownames(fit$coefficients)
  absent <- setdiff(required, fitted)
  if (length(absent)) .stop("Taxonomy '", taxonomy$name,
    "' requires terms absent from the fitted model: ", paste(absent, collapse = ", "), ".")
  invisible(NULL)
}

.taxonomy_design <- function(data, formula, family, exact) {
  if (!is.null(family)) return(.family_design(data, formula, family, exact))
  d <- if (exact) .design_exact(data, formula) else .design(data, formula)
  if (!all(c("step", "log_step", "cos_angle") %in% colnames(d$X)))
    .stop("The canonical taxonomy design must include step, log_step, and cos_angle.")
  d
}

.reference_state_labels <- list(
  HDT_F_R = c("HDT", "Foraging", "Resting"),
  T_F_R = c("Travelling", "Foraging", "Resting"),
  T_E = c("Travelling", "Encamped"),
  Other = "Other"
)
.reference_state_counts <- lengths(.reference_state_labels)

.check_custom_taxonomy <- function(fit, taxonomy, metrics, occupancy, checks, min_mean_gap) {
  k <- fit$n_states
  labels <- rep(NA_character_, k)
  if (!all(metrics$valid))
    return(list(labels = labels, checks = c(checks, state_order = FALSE)))

  order_by <- taxonomy$order_by
  score <- if (is.function(order_by)) order_by(metrics, fit$coefficients) else metrics[[order_by]]
  if (!is.numeric(score) || length(score) != k || !is.null(dim(score)))
    .stop("order_by must return one numeric score per fitted state.")
  checks <- c(checks, state_order = all(is.finite(score)) && !anyDuplicated(score))
  if (!checks[["state_order"]]) return(list(labels = labels, checks = checks))
  labels[order(score)] <- taxonomy$labels
  if (k > 1L && min_mean_gap > 0)
    checks <- c(checks, mean_gap = all(diff(sort(metrics$mean_step)) > min_mean_gap))
  if (is.null(taxonomy$rules)) return(list(labels = labels, checks = checks))

  rules <- tryCatch(taxonomy$rules(.taxonomy_context(fit, labels, metrics, occupancy)),
    error = function(e) .stop("Taxonomy '", taxonomy$name, "' rule failed: ",
      conditionMessage(e), ". Available coefficients: ",
      paste(rownames(fit$coefficients), collapse = ", "),
      "; state labels: ", paste(taxonomy$labels, collapse = ", "), "."))
  rule_names <- names(rules)
  if (!is.logical(rules) || !length(rules) || !is.null(dim(rules)) || is.null(rule_names) ||
      anyNA(rule_names) || any(!nzchar(trimws(rule_names))) || anyDuplicated(rule_names) ||
      any(rule_names %in% names(checks)))
    .stop("rules must return a nonempty named logical vector with unique names distinct from the built-in checks.")
  list(labels = labels, checks = c(checks, rules))
}

.taxonomy_context <- function(fit, labels, metrics, occupancy) {
  b <- fit$coefficients
  colnames(b) <- labels
  index <- data.frame(parameter = seq_along(b),
    state = rep(labels, each = nrow(b)), term = rep(rownames(b), ncol(b)),
    stringsAsFactors = FALSE)
  covariance <- if (is.null(fit$covariance)) NULL else
    fit$covariance[seq_along(b), seq_along(b), drop = FALSE]
  if (!is.null(covariance)) {
    keys <- make.unique(paste(index$state, index$term, sep = "::"))
    dimnames(covariance) <- list(keys, keys)
  }
  se <- b
  se[] <- if (is.null(covariance)) NA_real_ else sqrt(diag(covariance))
  rownames(metrics) <- labels
  transition <- fit$transition
  dimnames(transition)[1:2] <- list(labels, labels)
  list(metrics = metrics, coefficients = b, standard_errors = se,
    covariance = covariance, parameter_index = index,
    transition = transition, occupancy = stats::setNames(occupancy, labels))
}

#' Compute a coefficient contrast using labelled joint covariance
#'
#' Contrasts condition on the fitted model and the selected labels. Wald
#' intervals do not adjust for selection of a model, solution or taxonomy.
#' @param x A rule's input list, or the `context` returned by [assess_taxonomy()].
#' @param term One or more distinct coefficient row names.
#' @param from Source state label. With `to`, the contrast is `to - from`.
#' @param to Destination label, or NULL to report the coefficient in `from`.
#' @param level Confidence level strictly between zero and one.
#' @return A data frame with `term`, `from`, `to`, `estimate`, `std_error`,
#'   `lower` and `upper`. Uncertainty is NA if covariance is unavailable.
#' @export
ssf_contrast <- function(x, term, from, to = NULL, level = .95) {
  if (!is.list(x) || is.null(x$coefficients) || is.null(x$parameter_index))
    .stop("x must be a taxonomy rule input or assessment context.")
  if (!is.character(term) || !length(term) || anyNA(term) || anyDuplicated(term) ||
      !all(term %in% rownames(x$coefficients))) .stop("Unknown or duplicated contrast term.")
  labels <- colnames(x$coefficients)
  valid_label <- function(a) is.character(a) && length(a) == 1L && !is.na(a) && a %in% labels
  if (!valid_label(from) || (!is.null(to) && (!valid_label(to) || to == from)))
    .stop("from and to must identify fitted state labels.")
  .scalar(level, "level", .Machine$double.eps, 1 - .Machine$double.eps)
  index <- x$parameter_index
  do.call(rbind, lapply(term, function(name) {
    i <- index$parameter[index$state == from & index$term == name]
    j <- if (is.null(to)) integer() else index$parameter[index$state == to & index$term == name]
    weights <- numeric(nrow(index)); weights[i] <- if (is.null(to)) 1 else -1
    if (length(j)) weights[j] <- 1
    estimate <- sum(weights * c(x$coefficients))
    variance <- if (is.null(x$covariance)) NA_real_ else
      as.numeric(crossprod(weights, x$covariance %*% weights))
    se <- if (is.finite(variance) && variance >= 0) sqrt(variance) else NA_real_
    width <- stats::qnorm((1 + level) / 2) * se
    data.frame(term = name, from = from, to = if (is.null(to)) NA_character_ else to,
      estimate = estimate, std_error = se, lower = estimate - width, upper = estimate + width)
  }))
}
