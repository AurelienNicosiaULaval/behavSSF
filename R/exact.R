#' Prepare observed steps for analytically normalised movement kernels
#'
#' No controls or proposal weights are required. Each row is one positive
#' observed step, and each ID is an uninterrupted independent segment.
#' Directional target effects require the target bearing relative to the
#' previous heading, chosen at the step origin before observing the endpoint.
#' The constructor can derive this from target_bearing and previous_bearing.
#' It cannot establish from numbers alone how a supplied target was chosen.
#'
#' @param data Data frame with ID, step and angle columns. Angles are radians.
#'   Optional relative_target_bearing, or target_bearing and previous_bearing,
#'   describes the prospective target direction. Optional target_active is
#'   logical; FALSE removes its directional effect when direction is undefined.
#' @param id Name of the independent segment column.
#' @param time Optional numeric or POSIXct time column.
#' @param expected_interval Optional positive interval between consecutive
#'   observations in a segment, in seconds for POSIXct times.
#' @param density_measure Required: area or step_angle. This determines the
#'   reference measure of the exact likelihood, not an importance proposal.
#' @return An ssf_steps data frame with canonical covariates and metadata.
#' @export
prepare_ssf_steps <- function(data, density_measure, id = "ID", time = NULL,
                              expected_interval = NULL) {
  measure <- match.arg(density_measure, c("area", "step_angle"))
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      !all(c(id, "step", "angle", time) %in% names(data)))
    .stop("Observed data must contain ID, positive step and finite angle columns.")
  out <- as.data.frame(data)
  if ("obs" %in% names(out) && (anyNA(out$obs) || any(out$obs!=1)))
    .stop("prepare_ssf_steps() accepts observed rows only, without controls.")
  ids <- as.character(out[[id]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(rle(ids)$values))
    .stop("Each independent segment must occupy one contiguous block of observed steps.")
  finite <- function(x) is.numeric(x) && all(is.finite(x))
  if (!finite(out$step) || any(out$step <= 0) || !finite(out$angle))
    .stop("Observed steps must be positive and angles finite radians.")
  if (!is.null(expected_interval)) {
    .scalar(expected_interval, "expected_interval", .Machine$double.eps)
    if (is.null(time)) .stop("expected_interval requires a time column.")
  }
  if (!is.null(time)) {
    tm <- out[[time]]
    if (!(is.numeric(tm) || inherits(tm, "POSIXct")) || any(!is.finite(as.numeric(tm))))
      .stop("time must contain finite numeric values or POSIXct times.")
    intervals <- unlist(lapply(split(as.numeric(tm), factor(ids, levels = unique(ids))), diff),
                        use.names = FALSE)
    if (any(intervals <= 0)) .stop("Case times must increase within every segment.")
    if (!is.null(expected_interval) && any(abs(intervals-expected_interval) > 1e-8))
      .stop("A time gap occurs within a segment; split at missing observations.")
  }
  verify <- function(name, value) {
    if (name %in% names(out) && (!finite(out[[name]]) ||
        any(abs(out[[name]]-value) > 1e-10*(1+abs(value)))))
      .stop("Supplied ", name, " disagrees with the observed step geometry.")
    value
  }
  out$log_step <- verify("log_step", log(out$step))
  out$cos_angle <- verify("cos_angle", cos(out$angle))
  if (!"relative_target_bearing" %in% names(out) &&
      all(c("target_bearing", "previous_bearing") %in% names(out))) {
    if (!finite(out$target_bearing) || !finite(out$previous_bearing))
      .stop("Target and previous bearings must be finite radians.")
    out$relative_target_bearing <- out$target_bearing-out$previous_bearing
  }
  if ("relative_target_bearing" %in% names(out)) {
    if (!finite(out$relative_target_bearing)) .stop("Relative target bearings must be finite.")
    out$relative_target_bearing <- atan2(sin(out$relative_target_bearing),
                                        cos(out$relative_target_bearing))
    if (!"target_active" %in% names(out)) out$target_active <- rep(TRUE, nrow(out))
    if (!is.logical(out$target_active) || anyNA(out$target_active))
      .stop("target_active must be a nonmissing logical column.")
    value <- as.numeric(out$target_active)*cos(out$angle-out$relative_target_bearing)
    out$cos_diff_theta_to_target <- verify("cos_diff_theta_to_target", value)
  } else if ("target_active" %in% names(out)) {
    .stop("target_active requires a relative target bearing.")
  }
  out$ID <- ids
  out$stratum <- seq_len(nrow(out))
  out$.source_row <- seq_len(nrow(out))
  rownames(out) <- NULL
  class(out) <- c("ssf_steps", "data.frame")
  attr(out, "density_measure") <- measure
  out
}

.exact_terms <- c("step", "log_step", "cos_angle", "cos_diff_theta_to_target")
.step_ar_term <- "lag_log_step"
.step_day_lags <- 22:26
.step_day_terms <- paste0("lag", .step_day_lags, "_log_step")

.angular_log_z <- function(persistence, bias, psi, active) {
  x <- persistence + active * bias * cos(psi)
  y <- active * bias * sin(psi)
  concentration <- sqrt(x * x + y * y)
  log(2 * pi) + log(besselI(concentration, 0, expon.scaled = TRUE)) + concentration
}

.design_exact <- function(data, formula) {
  if (!inherits(data, "ssf_steps"))
    .stop("Exact fitting requires observed data from prepare_ssf_steps().")
  if (!inherits(formula, "formula") || length(formula) != 2L)
    .stop("formula must be one-sided.")
  formula_terms <- stats::terms(formula)
  if (length(attr(formula_terms,"offset"))) .stop("Exact kernels do not support formula offsets.")
  terms <- attr(formula_terms, "term.labels")
  allowed <- c(.exact_terms, .step_ar_term, .step_day_terms)
  if (!all(.exact_terms[1:3] %in% terms) || any(!terms %in% allowed))
    .stop("Exact normalisation supports step + log_step + cos_angle, with optional target, AR(1) and lags 22 to 26.")
  has_target <- .exact_terms[4] %in% terms
  if (has_target && !all(c("relative_target_bearing", "target_active") %in% names(data)))
    .stop("Exact target normalisation requires the target bearing relative to the previous heading at the step origin.")
  measure <- attr(data, "density_measure")
  if (length(measure) != 1L || !measure %in% c("area", "step_angle"))
    .stop("Exact fitting requires an explicit density measure.")
  # Recheck canonical geometry after subsetting or user edits.
  checked <- prepare_ssf_steps(as.data.frame(data), density_measure = measure)
  if (.step_ar_term %in% terms) {
    first <- !duplicated(checked$ID)
    if (all(first)) .stop("lag_log_step requires consecutive steps in at least one segment.")
    lag <- c(0, checked$log_step[-nrow(checked)])
    lag[first] <- 0
    checked[[.step_ar_term]] <- lag
  }
  day_terms <- intersect(.step_day_terms, terms)
  lag_rows <- stats::setNames(vector("list", length(day_terms)), day_terms)
  if (length(day_terms)) {
    if (!.step_ar_term %in% terms) .stop("Daily lag terms require lag_log_step.")
    age <- sequence(rle(checked$ID)$lengths)
    for (term in day_terms) {
      delay <- .step_day_lags[match(term, .step_day_terms)]
      i <- which(age > delay)
      lag <- numeric(nrow(checked))
      lag[i] <- checked$log_step[i-delay]
      checked[[term]] <- lag
      lag_rows[[term]] <- i
    }
  }
  x <- stats::model.matrix(formula, checked)
  x <- x[,colnames(x) != "(Intercept)",drop=FALSE]
  if (any(!is.finite(x))) .stop("Exact model covariates must be finite.")
  list(data = checked, X = x, cases = seq_len(nrow(x)), segment = checked$ID,
    indices = as.integer(match(allowed,colnames(x),nomatch=0L)-1L),
    lag_rows = lag_rows,
    relative_bearing = if (has_target) checked$relative_target_bearing else rep(0,nrow(x)),
    target_active = if (has_target) as.double(checked$target_active) else rep(0,nrow(x)),
    shape_offset = if (measure == "area") 2 else 1)
}

#' Exact log normalisers for canonical gamma-von Mises kernels
#'
#' Integrates exp(beta_step*r + beta_log*log(r) + a*cos(phi) +
#' b*cos(phi-psi)) on the full unbounded plane or in step-angle coordinates.
#' The angular concentration is sqrt((a+b*cos(psi))^2+(b*sin(psi))^2).
#' No endpoint habitat covariates, truncation or spatial barriers are supported.
#'
#' @param coefficients Coefficient matrix with canonical movement rows and an
#'   optional cos_diff_theta_to_target row.
#' @param relative_bearing Prospective target bearings relative to previous
#'   headings. Scalar or vector of length equal to target_active.
#' @param target_active Logical scalar or vector; FALSE removes the target term.
#' @param convention Required density measure: area or step_angle.
#' @return Matrix of log normalisers, rows for bearings and columns for states.
#' @export
ssf_log_normalizer <- function(coefficients, convention, relative_bearing = 0,
                               target_active = TRUE) {
  convention <- match.arg(convention,c("area","step_angle"))
  b <- as.matrix(coefficients)
  if (any(!rownames(b) %in% .exact_terms)) .stop("Unsupported exact kernel coefficients.")
  m <- movement_metrics(b,convention)
  if (!all(m$valid)) .stop("Exact kernels require positive gamma shape and rate.")
  n <- max(length(relative_bearing),length(target_active))
  if (!n || !is.numeric(relative_bearing) || any(!is.finite(relative_bearing)) ||
      !is.logical(target_active) || anyNA(target_active) ||
      any(!c(length(relative_bearing),length(target_active)) %in% c(1L,n)))
    .stop("Target bearings and activation must have compatible finite lengths.")
  psi <- rep(relative_bearing,length.out=n)
  active <- as.numeric(rep(target_active,length.out=n))
  bias <- if (.exact_terms[4] %in% rownames(b)) b[.exact_terms[4],] else rep(0,ncol(b))
  answer <- vapply(seq_len(ncol(b)),function(k) {
    value <- lgamma(m$shape[k])-m$shape[k]*log(m$rate[k]) +
      .angular_log_z(m$persistence[k], bias[k], psi, active)
    if (any(!is.finite(value))) .stop("Nonfinite exact normaliser.")
    value
  },numeric(n))
  matrix(answer,nrow=n,ncol=ncol(b),dimnames=list(NULL,colnames(b)))
}

.exact_emissions <- function(design, coefficients, step_offset = NULL) {
  measure <- attr(design$data, "density_measure")
  ar <- .step_ar_term %in% rownames(coefficients)
  if (!ar && is.null(step_offset)) {
    return(design$X %*% coefficients - ssf_log_normalizer(coefficients,
      measure, design$relative_bearing, as.logical(design$target_active)))
  }
  rho <- if (ar) coefficients[.step_ar_term, ] else rep(0, ncol(coefficients))
  if (any(!is.finite(rho) | abs(rho) >= 1)) .stop("AR coefficients must lie strictly between -1 and 1.")
  day_terms <- names(design$lag_rows)
  if (length(day_terms) &&
      any(abs(rho) + colSums(abs(coefficients[day_terms, , drop = FALSE])) >= 1))
    .stop("The absolute AR coefficients must sum to less than one in each state.")
  within <- which(duplicated(design$segment))
  r <- design$data$step
  lag <- if (ar) design$X[within, .step_ar_term] else numeric()
  radial <- movement_metrics(coefficients, measure)
  if (!all(radial$valid)) .stop("Exact kernels require positive gamma shape and rate.")
  angular <- intersect(.exact_terms[3:4], rownames(coefficients))
  emissions <- design$X[, angular, drop = FALSE] %*% coefficients[angular, , drop = FALSE]
  bias <- if (.exact_terms[4] %in% angular) coefficients[.exact_terms[4], ] else rep(0, length(rho))
  jacobian <- if (measure == "area") log(r) else 0
  for (state in seq_along(rho)) {
    log_rate <- log(radial$rate[state])
    log_anchor <- log(radial$shape[state]) - log_rate
    rate <- rep(radial$rate[state], length(r))
    if (rho[state] != 0) rate[within] <- exp(log_rate - rho[state] * (lag - log_anchor))
    for (term in day_terms) {
      coefficient <- coefficients[term, state]
      if (coefficient == 0) next
      i <- design$lag_rows[[term]]
      rate[i] <- rate[i]*exp(-coefficient*(design$X[i, term] - log_anchor))
    }
    if (!is.null(step_offset)) rate <- rate*exp(-step_offset[, state])
    if (any(!is.finite(rate) | rate <= 0)) .stop("Nonfinite or nonpositive conditional gamma rate.")
    emissions[, state] <- emissions[, state] + stats::dgamma(r, radial$shape[state], rate = rate, log = TRUE) -
      jacobian - .angular_log_z(radial$persistence[state], bias[state],
        design$relative_bearing, design$target_active)
  }
  emissions
}

.ar_indices <- function(design, k) {
  cols <- match(c(.step_ar_term, .step_day_terms), colnames(design$X), nomatch = 0L)
  outer(cols[cols > 0L], ncol(design$X) * (seq_len(k) - 1L), "+")
}

# Invert the smooth map rho_i = sinh(w_i) / sum(cosh(w_j)).
.ar_inverse_scale <- function(rho) {
  stats::uniroot(function(u) sum(sqrt(rho^2 + u^2)) - 1,
    c(0, 2/length(rho)), tol = .Machine$double.eps)$root
}

.ar_block <- function(x, inverse) {
  if (!length(x)) return(x)
  if (inverse && sum(abs(x)) >= 1)
    .stop("The initial absolute AR coefficients must sum to less than one in each state.")
  if (length(x) == 1L) return(if (inverse) atanh(x) else tanh(x))
  if (length(x) == 2L) {
    transformed <- if (inverse) atanh(c(x[1]+x[2], x[1]-x[2])) else
      tanh(c(x[1]+x[2], x[1]-x[2]))
    return(c(sum(transformed), diff(rev(transformed)))/2)
  }
  if (inverse) return(asinh(x/.ar_inverse_scale(x)))
  scale <- max(abs(x))
  pos <- exp(x-scale)
  neg <- exp(-x-scale)
  (pos-neg)/sum(pos+neg)
}

.exact_transform <- function(theta, design, k, inverse = FALSE) {
  p <- ncol(design$X)
  step_index <- which(colnames(design$X) == "step") + p*(seq_len(k)-1L)
  log_index <- which(colnames(design$X) == "log_step") + p*(seq_len(k)-1L)
  if (inverse) {
    if (any(theta[step_index] >= 0) || any(theta[log_index]+design$shape_offset <= 0))
      .stop("Initial exact kernels require positive gamma shape and rate.")
    theta[step_index] <- log(-theta[step_index])
    theta[log_index] <- log(theta[log_index]+design$shape_offset)
  } else {
    theta[step_index] <- -exp(theta[step_index])
    theta[log_index] <- exp(theta[log_index])-design$shape_offset
  }
  ar <- .ar_indices(design, k)
  for (state in seq_len(k)) theta[ar[, state]] <- .ar_block(theta[ar[, state]], inverse)
  theta
}

.exact_jacobian <- function(natural, design, k) {
  p <- ncol(design$X)
  jac <- rep(1, length(natural))
  si <- which(colnames(design$X) == "step") + p*(seq_len(k)-1L)
  li <- which(colnames(design$X) == "log_step") + p*(seq_len(k)-1L)
  jac[si] <- natural[si]
  jac[li] <- natural[li]+design$shape_offset
  ar <- .ar_indices(design, k)
  if (nrow(ar) <= 1L) {
    jac[c(ar)] <- 1-natural[c(ar)]^2
    return(jac)
  }
  J <- diag(jac)
  for (state in seq_len(k)) {
    i <- ar[, state]
    rho <- natural[i]
    if (length(rho) == 2L) {
      plus <- 1 - (rho[1]+rho[2])^2
      minus <- 1 - (rho[1]-rho[2])^2
      J[i, i] <- matrix(c(plus+minus, plus-minus, plus-minus, plus+minus)/2, 2L)
    } else {
      u <- .ar_inverse_scale(rho)
      J[i, i] <- diag(sqrt(rho^2+u^2)) - outer(rho, rho)
    }
  }
  J
}
