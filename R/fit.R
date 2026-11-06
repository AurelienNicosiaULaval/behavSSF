.unpack <- function(theta, p, k, names) {
  beta <- matrix(theta[seq_len(p * k)], p, k,
                 dimnames = list(names, paste0("S", seq_len(k))))
  tpm <- matrix(1, k, k)
  if (k > 1L) {
    alpha <- matrix(theta[-seq_len(p * k)], k, k - 1L, byrow = TRUE)
    for (i in seq_len(k)) {
      logits <- numeric(k)
      logits[-i] <- alpha[i, ]
      tpm[i, ] <- exp(logits - .lse(logits))
    }
  }
  list(coefficients = beta, transition = tpm)
}

.optim_attempt <- function(theta, objective, maxit, phase, start, hessian = FALSE,
                           gradient = NULL) {
  warnings <- character()
  error <- NULL
  value <- tryCatch(withCallingHandlers(
    stats::optim(theta, objective, gr = gradient, method = "BFGS", hessian = hessian,
                 control = list(maxit = maxit, reltol = 1e-11)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e); NULL })
  valid <- !is.null(value) && is.finite(value$value) && value$value < 1e99 &&
    all(is.finite(value$par))
  log <- data.frame(phase = phase, start = start,
                    objective = if (valid) value$value else NA_real_,
                    convergence = if (is.null(value)) NA_integer_ else value$convergence,
                    error = if (is.null(error)) "" else error,
                    warnings = paste(unique(warnings), collapse = "; "))
  list(fit = if (valid) value else NULL, log = log)
}

.ssf_working_score <- function(par, objective) {
  if (!is.null(objective$gr)) return(objective$gr(par))
  vapply(seq_along(par), function(i) {
    h <- 1e-5 * max(1, abs(par[i]))
    plus <- minus <- par
    plus[i] <- plus[i] + h
    minus[i] <- minus[i] - h
    (objective$fn(plus) - objective$fn(minus)) / (2 * h)
  }, numeric(1))
}

.positive_ssf_hessian <- function(h) {
  if (is.null(h) || any(!is.finite(h))) return(FALSE)
  eig <- eigen((h + t(h)) / 2, symmetric = TRUE, only.values = TRUE)$values
  min(eig) > 1e-8 * max(1, max(eig))
}

.ssf_natural_score <- function(working, jacobian) {
  if (!is.matrix(jacobian)) return(working/jacobian)
  tryCatch(solve(t(jacobian), working), error = function(e) rep(NA_real_, length(working)))
}

.refine_ssf_score <- function(fit, objective, jacobian_at) {
  working <- .ssf_working_score(fit$par, objective)
  score <- .ssf_natural_score(working, jacobian_at(fit$par))
  before <- max(abs(score))
  steps <- 0L
  # BFGS stops on relative changes in the objective. A long sequence can
  # meet that criterion with a remaining score large enough to fail assessment.
  for (iteration in seq_len(3L)) {
    largest <- max(abs(score))
    if (fit$convergence != 0L || !is.finite(largest) || largest <= 1e-3 ||
        !.positive_ssf_hessian(fit$hessian)) break
    direction <- solve((fit$hessian + t(fit$hessian)) / 2, working)
    roundoff <- 64 * .Machine$double.eps * max(1, abs(fit$value))
    accepted <- FALSE
    for (fraction in c(1, .5, .25, .125)) {
      par <- fit$par - fraction * direction
      if (any(!is.finite(par))) next
      value <- objective$fn(par)
      if (!is.finite(value) || value >= 1e99 || value > fit$value + roundoff) next
      trial_working <- .ssf_working_score(par, objective)
      trial_score <- .ssf_natural_score(trial_working, jacobian_at(par))
      if (any(!is.finite(trial_score)) || max(abs(trial_score)) >= largest) next
      hessian <- tryCatch(stats::optimHess(par, objective$fn, gr = objective$gr),
        error = function(e) NULL)
      if (!.positive_ssf_hessian(hessian)) next
      fit$par <- par
      fit$value <- value
      fit$hessian <- hessian
      working <- trial_working
      score <- trial_score
      steps <- steps + 1L
      accepted <- TRUE
      break
    }
    if (!accepted) break
  }
  list(fit = fit, working = working, score = score, before = before, steps = steps)
}

#' Fit a hidden Markov step selection model
#'
#' Uses either the published SSF importance approximation or exact canonical
#' kernel normalisation, with the HMM forward algorithm,
#' through either a compiled analytic-gradient engine or an independent R reference. All fitted SSF and transition
#' parameters are passed from the best short run to the long run. Errors and
#' warnings are retained for every attempt. This implementation prioritises
#' verifiability, with cross-checks between the two engines.
#'
#' @param data Case-control data accepted by [prepare_ssf_data()], or observed
#'   [prepare_ssf_steps()] data when normalization is exact.
#' @param formula One-sided SSF formula. For exact fitting, add `lag_log_step`
#'   to estimate a state-specific AR(1) coefficient for the log gamma mean.
#'   The previous length is used within a segment, including across state
#'   changes. Add any of `lag22_log_step` through `lag26_log_step` with
#'   `lag_log_step` and `transition_hour` to include those hourly lags.
#'   A lag only contributes when its step is available in the same segment.
#'   The absolute lag coefficients sum to less than one in each state.
#' @param n_states Number of hidden states, between 1 and 6.
#' @param n_starts Number of short runs. A one-state model uses one run.
#' @param short_iter,long_iter Maximum BFGS iterations for short and long phases.
#' @param spread Positive standard deviation of initial perturbations. A scalar
#'   or named vector with one entry for each SSF design column.
#' @param seed Optional integer seed, restored to the previous RNG state on exit.
#' @param engine `"native"` uses the compiled objective and analytic gradient;
#'   `"R"` uses the independent R objective with numerical derivatives.
#' @param initial_values Optional list with a coefficient matrix and positive
#'   transition matrix, named `coefficients` and `transition`. Included as the
#'   first short-run start; remaining starts are random perturbations.
#'   For daily transitions, supply `transition_coefficients` from a previous
#'   fit with the same transition harmonics in place of `transition`.
#'   Supply `step_hour_coefficients` from a fit with the same length
#'   harmonics to reuse those effects. They start at zero by default.
#' @param normalization importance uses case-control data; exact requires
#'   observed data from prepare_ssf_steps() and supports only the canonical
#'   gamma-von Mises terms and optional prospective target direction.
#'   Exact fitting optimises log rate and log shape to preserve positive gamma
#'   parameters. AR coefficients use a hyperbolic tangent link. Spread applies
#'   on the working scale. Reported coefficients,
#'   covariance and score are on the natural coefficient scale.
#' @param transition_hour Optional integer hour (0 to 23) at the beginning
#'   of each step. Daily transitions require two states and exact fitting
#'   on hourly observations. Supply hours in a consistent time zone.
#' @param transition_harmonics Number of sine/cosine pairs, from 1 to 5,
#'   when `transition_hour` is supplied. Each switch direction has its own
#'   logit regression. Each segment starts in the periodic stationary
#'   distribution for its initial hour.
#' @param step_harmonics Number of sine/cosine pairs in the log mean of
#'   the gamma step length, from 0 to 5. A positive value requires
#'   `transition_hour`. Each state has its own curve. The curve adds to
#'   the log mean after any AR(1) adjustment. The gamma shape remains fixed
#'   within each state.
#' @return An `ssf_fit` object containing coefficients, transition matrix,
#'   log likelihood, covariance when the Hessian is positive definite, decoded
#'   states, smoothed probabilities, diagnostics, and the full start log.
#'   With daily transitions, `transition` is a 2 by 2 by 24 array indexed
#'   by step-origin hour. `transition_coefficients` contains the two
#'   switching regressions and `initial` gives state probabilities by hour.
#'   If `step_harmonics` is positive, `step_hour_coefficients` holds the
#'   hourly log-mean coefficients, one column per state. The movement
#'   coefficients describe the baseline before the hourly adjustment.
#'   Use [decode_ssf_diel()] to score new observations.
#'   After BFGS converges, a local Newton correction can reduce the score when
#'   the Hessian is positive definite. It accepts a correction when the
#'   largest absolute score falls and the objective improves or agrees to
#'   rounding precision. The Hessian is recomputed at the retained solution. The
#'   `score_refinement` diagnostic records the largest absolute score before
#'   this step and the number of accepted corrections. The convergence code
#'   and the original BFGS logs retain their optimiser meaning.
#'   For exact fits the stored Hessian uses working log-rate/log-shape parameters;
#'   covariance and score use the reported natural coefficient scale.
#' @export
fit_ssf <- function(data, formula = ~ step + log_step + cos_angle, n_states = 2L,
                    n_starts = 5L, short_iter = 30L, long_iter = 300L,
                    spread = 0.5, seed = NULL, engine = c("native", "R"),
                    initial_values = NULL, normalization = c("importance", "exact"),
                    transition_hour = NULL, transition_harmonics = 2L, step_harmonics = 0L) {
  user_call <- match.call()
  engine <- match.arg(engine)
  normalization <- match.arg(normalization)
  exact <- normalization == "exact"
  .scalar(n_states, "n_states", 1, 6, TRUE)
  .scalar(n_starts, "n_starts", 1, 10000, TRUE)
  .scalar(short_iter, "short_iter", 1, 1e6, TRUE)
  .scalar(long_iter, "long_iter", 1, 1e6, TRUE)
  d <- if (exact) .design_exact(data, formula) else .design(data, formula)
  p <- ncol(d$X)
  k <- as.integer(n_states)
  diel <- !is.null(transition_hour)
  daily_lags <- exact && length(d$lag_rows) > 0L
  if (daily_lags && !diel) .stop("Daily lag terms require transition_hour.")
  if (daily_lags && any(lengths(d$lag_rows) == 0L))
    .stop("Daily lags need a segment with more hourly steps than the largest lag.")
  .scalar(step_harmonics, "step_harmonics", 0, 5, TRUE)
  step_harmonics <- as.integer(step_harmonics)
  if (step_harmonics > 0L && !diel) .stop("Hourly step lengths require transition_hour.")
  if (diel) {
    if (!exact || k != 2L) .stop("Daily transitions require an exact two-state fit.")
    .scalar(transition_harmonics, "transition_harmonics", 1, 5, TRUE)
    d$diel <- list(hour = .diel_hour(transition_hour, nrow(d$X), d$segment),
                   basis = .diel_basis(transition_harmonics))
    if (length(unique(d$diel$hour[duplicated(d$segment)])) < ncol(d$diel$basis))
      .stop("Too few distinct transition hours to estimate the requested harmonics.")
    if (step_harmonics > 0L) {
      d$diel$step_basis <- .diel_basis(step_harmonics)[, -1L, drop = FALSE]
      hour <- unique(d$diel$hour) + 1L
      if (qr(cbind(1, d$diel$step_basis[hour, , drop = FALSE]))$rank < 2L*step_harmonics + 1L)
        .stop("Too few distinct step hours to estimate the requested length harmonics.")
    }
  }
  if (length(spread) > 1L) {
    if (is.null(names(spread)) || anyDuplicated(names(spread)) ||
        !setequal(names(spread), colnames(d$X)))
      .stop("A vector spread must be named for all SSF design columns.")
    spread <- spread[colnames(d$X)]
  }
  if (!is.numeric(spread) || !length(spread) || any(!is.finite(spread)) || any(spread <= 0))
    .stop("spread must contain strictly positive finite values.")
  spread <- rep(spread, length.out = p)
  unpack <- function(theta, kk) {
    if (exact) theta <- .exact_transform(theta,d,kk)
    if (!diel || kk == 1L) return(.unpack(theta, p, kk, colnames(d$X)))
    states <- paste0("S", 1:2)
    beta <- matrix(theta[seq_len(2L*p)], p, 2L, dimnames = list(colnames(d$X), states))
    transition_index <- 2L*p + seq_len(2L*ncol(d$diel$basis))
    gamma <- matrix(theta[transition_index], ncol(d$diel$basis), 2L,
      dimnames = list(colnames(d$diel$basis), states))
    z <- c(list(coefficients = beta, transition_coefficients = gamma), .diel_schedule(d$diel$basis, gamma))
    if (step_harmonics > 0L) z$step_hour_coefficients <- matrix(theta[-c(seq_len(2L*p), transition_index)],
      2L*step_harmonics, 2L, dimnames = list(colnames(d$diel$step_basis), states))
    z
  }
  emissions <- function(z) {
    if (!exact) return(.emissions(d, z$coefficients))
    offset <- if (!is.null(z$step_hour_coefficients))
      d$diel$step_basis[d$diel$hour + 1L, , drop = FALSE] %*% z$step_hour_coefficients
    .exact_emissions(d, z$coefficients, offset)
  }
  step_start <- numeric(2L*k*step_harmonics)
  objective_for <- function(kk) {
    if (engine == "native") return(if (exact) .native_exact_objective(d,kk) else
                                    .native_objective(d, kk))
    fn <- function(theta) {
      z <- unpack(theta,kk)
      initial <- if (is.null(z$initial)) .stationary(z$transition) else z$initial
      hour <- if (!is.null(z$initial)) d$diel$hour + 1L else NULL
      -.forward(emissions(z), z$transition, d$segment, initial, hour = hour)
    }
    safe_fn <- function(theta) {
      value <- tryCatch(fn(theta), error = function(e) Inf)
      if (is.finite(value)) value else 1e100
    }
    list(fn = safe_fn, gr = NULL)
  }
  supplied <- NULL
  if (!is.null(initial_values)) {
    b <- initial_values$coefficients
    transition <- initial_values[["transition"]]
    if (!is.matrix(b) || !is.numeric(b) || !identical(dim(b), c(p, k)) ||
        !identical(rownames(b), colnames(d$X)) || any(!is.finite(b)))
      .stop("Initial coefficients must match the design rows and state count.")
    gamma <- initial_values$transition_coefficients
    if (is.null(gamma)) {
      .check_hmm(matrix(0, 1, k), transition, 1, NULL)
      if (any(transition <= 0)) .stop("Initial transitions must be strictly positive.")
      alpha <- unlist(lapply(seq_len(k), function(i)
        log(transition[i, -i] / transition[i, i])), use.names = FALSE)
      if (diel) alpha <- c(rbind(alpha, matrix(0, ncol(d$diel$basis) - 1L, k)))
    } else {
      if (!diel || !is.null(transition) || !is.numeric(gamma) ||
          !identical(dim(gamma), c(ncol(d$diel$basis), k)) ||
          !identical(rownames(gamma), colnames(d$diel$basis)) || any(!is.finite(gamma)))
        .stop("Initial transition_coefficients must match the daily transition basis.")
      alpha <- c(gamma)
    }
    initial_step_hour <- initial_values[["step_hour_coefficients"]]
    if (!is.null(initial_step_hour)) {
      if (step_harmonics == 0L || !is.numeric(initial_step_hour) ||
          !identical(dim(initial_step_hour), c(2L*step_harmonics, k)) ||
          !identical(rownames(initial_step_hour), colnames(d$diel$step_basis)) ||
          any(!is.finite(initial_step_hour)))
        .stop("Initial step_hour_coefficients must match the daily length basis.")
    }
    supplied <- c(b, alpha, if (is.null(initial_step_hour)) step_start else c(initial_step_hour))
    if (exact) supplied <- .exact_transform(supplied,d,k,inverse=TRUE)
  }
  .with_seed(seed, function() {
    baseline_objective <- objective_for(1L)
    base_start <- rep(0,p)
    if (exact) {
      mu <- mean(d$data$step)
      variance <- stats::var(d$data$step)
      shape <- if (is.finite(variance) && variance>0) max(.1,mu^2/variance) else 1
      base_start[match("step",colnames(d$X))] <- -shape/mu
      base_start[match("log_step",colnames(d$X))] <- shape-d$shape_offset
      base_start <- .exact_transform(base_start,d,1L,inverse=TRUE)
    }
    baseline <- .optim_attempt(base_start, baseline_objective$fn, long_iter,
                                "baseline", 0L, gradient = baseline_objective$gr)
    if (is.null(baseline$fit)) stop(errorCondition(
      "The one-state initialisation failed.", class = "ssf_fitting_error",
      attempts = baseline$log))
    objective <- objective_for(k)
    if (k == 1L) n_starts <- 1L
    starts <- lapply(seq_len(n_starts), function(s) {
      if (s == 1L && !is.null(supplied)) return(supplied)
      beta <- matrix(baseline$fit$par, p, k)
      if (k > 1L) beta <- beta + matrix(stats::rnorm(p * k, sd = rep(spread, k)), p, k)
      alpha <- if (diel) c(rbind(rep(-2, k), matrix(0, ncol(d$diel$basis) - 1L, k))) else
        rep(-2, k * (k - 1L))
      c(beta, alpha, step_start)
    })
    attempts <- lapply(seq_along(starts), function(s)
      .optim_attempt(starts[[s]], objective$fn, short_iter, "short", s,
                      gradient = objective$gr))
    values <- vapply(attempts, function(x) if (is.null(x$fit)) Inf else x$fit$value, numeric(1))
    if (all(!is.finite(values))) stop(errorCondition(
      "Every short run failed.", class = "ssf_fitting_error",
      attempts = do.call(rbind, c(list(baseline$log), lapply(attempts, `[[`, "log")))))
    best <- which.min(values)
    short <- attempts[[best]]$fit
    # Continue from all fitted parameters, including transitions.
    long <- .optim_attempt(short$par, objective$fn, long_iter, "long", best, TRUE,
                            gradient = objective$gr)
    retained <- !is.null(long$fit) && long$fit$value <= short$value + 1e-8
    fit <- if (retained) long$fit else short
    if (!retained) fit$hessian <- tryCatch(stats::optimHess(fit$par, objective$fn, gr = objective$gr),
                                           error = function(e) NULL)
    jacobian_at <- function(par) {
      if (exact) .exact_jacobian(.exact_transform(par, d, k), d, k) else rep(1, length(par))
    }
    refinement <- .refine_ssf_score(fit, objective, jacobian_at)
    fit <- refinement$fit
    z <- unpack(fit$par,k)
    log_e <- emissions(z)
    decode <- if (diel) .forward(log_e, z$transition, d$segment, z$initial,
      decode = TRUE, hour = d$diel$hour + 1L) else decode_states(log_e, z$transition, d$segment)
    score <- refinement$score
    h <- fit$hessian
    positive <- .positive_ssf_hessian(h)
    covariance <- if (positive) solve(h) else NULL
    if (positive && exact) {
      jacobian <- jacobian_at(fit$par)
      covariance <- if (is.matrix(jacobian)) jacobian %*% covariance %*% t(jacobian) else
        covariance * outer(jacobian, jacobian)
    }
    out <- list(coefficients = z$coefficients, transition = z$transition,
                log_likelihood = decode$log_likelihood, n_states = k,
                state = decode$state, posterior = decode$posterior,
                observations = d$data[d$cases, , drop = FALSE],
                covariance = covariance, hessian = h,
                diagnostics = list(converged = fit$convergence == 0L,
                  convergence_code = fit$convergence, hessian_positive = positive,
                  long_run_retained = retained, best_start = best,
                  short_parameters = short$par, long_start = short$par,
                  objective_agreement = abs(fit$value + decode$log_likelihood),
                  score = score, max_abs_gradient = max(abs(score)),
                  working_score = refinement$working,
                  score_refinement = list(before = refinement$before, steps = refinement$steps),
                  objective = if (!is.null(objective$diagnostics)) objective$diagnostics() else NULL),
                starts = do.call(rbind, c(list(baseline$log),
                             lapply(attempts, `[[`, "log"), list(long$log))),
                formula = formula, seed = seed, engine = engine,
                density_measure = attr(d$data, "density_measure"),
                normalization = normalization,
                hessian_parameterization = if (daily_lags) "log_rate_log_shape_joint_ar" else
                  if (exact && .step_ar_term %in% colnames(d$X))
                  "log_rate_log_shape_atanh_ar" else if (exact) "log_rate_log_shape" else "natural",
                covariance_parameterization = "natural",
                approximation = if (exact) "analytic full-plane normalisation; no controls" else
                  "case included; unweighted numerator; denominator / J",
                call = user_call)
    if (diel) {
      out$transition_coefficients <- z$transition_coefficients
      out$initial <- z$initial
      if (step_harmonics > 0L) out$step_hour_coefficients <- z$step_hour_coefficients
    }
    class(out) <- "ssf_fit"
    out
  })
}

#' @export
print.ssf_fit <- function(x, ...) {
  heading <- if (is.null(x$transition_coefficients)) "Homogeneous HMM-SSF:" else "Daily-transition HMM-SSF:"
  cat(heading, x$n_states, "state(s),",
      nrow(x$observations), "observed steps\n")
  cat("Log likelihood:", format(x$log_likelihood),
      "| converged:", x$diagnostics$converged,
      "| positive Hessian:", x$diagnostics$hessian_positive, "\n")
  invisible(x)
}

#' @export
coef.ssf_fit <- function(object, ...) object$coefficients

#' @export
logLik.ssf_fit <- function(object, ...) {
  transition_df <- if (is.null(object$transition_coefficients)) object$n_states * (object$n_states - 1L) else
    length(object$transition_coefficients)
  structure(object$log_likelihood,
            df = length(object$coefficients) + transition_df + length(object$step_hour_coefficients),
            nobs = nrow(object$observations), class = "logLik")
}
