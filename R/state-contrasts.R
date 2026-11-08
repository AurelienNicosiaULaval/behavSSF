#' Compare model effects between inferred states
#'
#' Returns differences between state-specific coefficients using their joint
#' covariance. This applies to any estimable numeric, factor or interaction
#' term in the fitted model. The intervals condition on that model and its
#' assigned labels; they do not adjust for choosing states, rules or contrasts.
#'
#' @inheritParams ssf_state_summary
#' @param terms Coefficient names. The default includes all non-radial terms.
#' @param comparisons Optional data frame with `from` and `to` state labels.
#'   Each contrast is `to - from`. The default gives every unordered pair once.
#' @param level Confidence level, strictly between zero and one.
#' @return A data frame with coefficient names, state pairs, differences,
#'   standard errors and confidence limits. Missing covariance gives NA limits.
#' @export
ssf_state_contrasts <- function(x, terms = NULL, comparisons = NULL,
                                assessment = NULL, level = .95) {
  resolved <- .state_display_fit(x, assessment)
  fit <- resolved$fit
  labels <- resolved$labels
  if (length(labels) < 2L) .stop("At least two fitted states are needed for a contrast.")
  radial <- c("step", "log_step", unique(unlist(lapply(.step_families, names))))
  b <- stats::coef(fit)
  if (is.null(terms)) terms <- setdiff(rownames(b), radial)
  if (!is.character(terms) || !length(terms) || anyNA(terms) || anyDuplicated(terms) ||
      !all(terms %in% rownames(b))) .stop("terms must contain distinct fitted coefficient names.")
  if (is.null(comparisons)) {
    pairs <- utils::combn(labels, 2L)
    comparisons <- data.frame(from = pairs[1L, ], to = pairs[2L, ])
  }
  if (!is.data.frame(comparisons) || !all(c("from", "to") %in% names(comparisons)) ||
      !nrow(comparisons) || anyNA(comparisons[c("from", "to")]) ||
      !all(comparisons$from %in% labels) || !all(comparisons$to %in% labels) ||
      any(comparisons$from == comparisons$to) || anyDuplicated(comparisons[c("from", "to")]))
    .stop("comparisons must give distinct from/to pairs of fitted state labels.")
  colnames(b) <- labels
  covariance <- fit[["covariance"]]
  context <- list(coefficients = b,
    parameter_index = data.frame(parameter = seq_along(b),
      state = rep(labels, each = nrow(b)), term = rep(rownames(b), ncol(b))),
    covariance = if (is.null(covariance)) NULL else
      covariance[seq_along(b), seq_along(b), drop = FALSE])
  do.call(rbind, lapply(seq_len(nrow(comparisons)), function(i)
    ssf_contrast(context, terms, as.character(comparisons$from[i]),
      as.character(comparisons$to[i]), level)))
}
