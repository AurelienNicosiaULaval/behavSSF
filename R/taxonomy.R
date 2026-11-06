#' Define a configurable state taxonomy
#'
#' Labels are assigned in increasing order of a fitted state signature.
#' A rule can then check separate movement, selection, or transition contrasts.
#' The same numerical and occupancy checks as in the guided workflow apply.
#'
#' @param name Unique identifier for the candidate. The four built-in names
#'   are reserved.
#' @param labels One to six unique state labels, in increasing signature order.
#' @param formula One-sided SSF formula, including `step`, `log_step`, and
#'   `cos_angle`. It must be supported by the chosen normalisation.
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
  if (!inherits(formula, "formula") || length(formula) != 2L ||
      !all(c("step", "log_step", "cos_angle") %in% attr(stats::terms(formula), "term.labels")))
    .stop("formula must be one-sided and include step, log_step, and cos_angle.")
  metrics <- c("mean_step", "sd_step", "persistence", "concentration", "shape", "rate")
  if (!is.function(order_by) && !(is.character(order_by) && length(order_by) == 1L &&
      order_by %in% metrics))
    .stop("order_by must be a numeric movement metric name or a function.")
  if (!is.null(rules) && !is.function(rules)) .stop("rules must be NULL or a function.")
  structure(list(name = name, labels = trimws(labels), formula = formula,
                 order_by = order_by, rules = rules), class = "ssf_taxonomy")
}

.reference_state_counts <- c(HDT_F_R = 3L, T_F_R = 3L, T_E = 2L, Other = 1L)

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

  b <- fit$coefficients
  colnames(b) <- labels
  se <- b
  se[] <- if (is.null(fit$covariance)) NA_real_ else sqrt(diag(fit$covariance)[seq_along(b)])
  rownames(metrics) <- labels
  transition <- fit$transition
  dimnames(transition)[1:2] <- list(labels, labels)
  rules <- taxonomy$rules(list(metrics = metrics, coefficients = b, standard_errors = se,
                              transition = transition, occupancy = stats::setNames(occupancy, labels)))
  rule_names <- names(rules)
  if (!is.logical(rules) || !length(rules) || !is.null(dim(rules)) || is.null(rule_names) ||
      anyNA(rule_names) || any(!nzchar(trimws(rule_names))) || anyDuplicated(rule_names) ||
      any(rule_names %in% names(checks)))
    .stop("rules must return a nonempty named logical vector with unique names distinct from the built-in checks.")
  list(labels = labels, checks = c(checks, rules))
}
