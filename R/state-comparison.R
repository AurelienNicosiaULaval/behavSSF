#' Make consecutive validation blocks without splitting segments
#'
#' Input order defines chronology. Cuts are placed near equally spaced
#' cumulative step counts, subject to leaving at least one segment per block.
#' @param data Prepared `ssf_data` or `ssf_steps`.
#' @param n_blocks Number of validation blocks, at least two.
#' @return An integer vector named by segment ID, for [compare_ssf_states()].
#' @export
ssf_blocks <- function(data, n_blocks = 3L) {
  obs <- .comparison_observations(data)
  ids <- unique(as.character(obs$ID))
  .scalar(n_blocks, "n_blocks", 2, length(ids), TRUE)
  counts <- tabulate(match(as.character(obs$ID), ids), length(ids))
  ends <- cumsum(counts)
  cuts <- integer(n_blocks - 1L)
  for (i in seq_along(cuts)) {
    lower <- if (i == 1L) 1L else cuts[i - 1L] + 1L
    candidates <- seq.int(lower, length(ids) - (n_blocks - i))
    cuts[i] <- candidates[which.min(abs(ends[candidates] - nrow(obs) * i / n_blocks))]
  }
  stats::setNames(findInterval(seq_along(ids) - 1L, cuts) + 1L, ids)
}

.comparison_observations <- function(data) {
  if (!inherits(data, "ssf_data") && !inherits(data, "ssf_steps"))
    .stop("data must be prepared with prepare_ssf_data() or prepare_ssf_steps().")
  if (inherits(data, "ssf_data")) data[data$obs == 1L, , drop = FALSE] else data
}

.comparison_screen <- function(fit, score_tolerance, objective_tolerance) {
  d <- fit$diagnostics
  point <- is.finite(as.numeric(stats::logLik(fit))) && isTRUE(d$converged) &&
    isTRUE(.assessment_score(fit) <= score_tolerance) &&
    isTRUE(d$objective_agreement <= objective_tolerance) &&
    all(movement_metrics(stats::coef(fit), fit$density_measure, fit$step_family,
      if (is.null(fit$step_max)) Inf else fit$step_max)$valid)
  c(point = point, strict = point && isTRUE(d$hessian_positive))
}

# Score entire held-out segments, never a decoded state sequence.
.comparison_score <- function(fit, data) {
  obs <- .comparison_observations(data)
  emissions <- ssf_log_emissions(data, fit$formula, stats::coef(fit),
    normalization = fit$normalization, step_family = fit$step_family)
  if (fit$normalization == "importance") {
    counts <- tabulate(match(data$stratum, obs$stratum), nrow(obs)) - 1L
    emissions <- emissions - log(counts * obs$w)
  }
  do.call(rbind, lapply(unique(as.character(obs$ID)), function(id) {
    use <- as.character(obs$ID) == id
    score <- decode_states(emissions[use, , drop = FALSE], fit$transition)$log_likelihood
    data.frame(segment = id, steps = sum(use), log_score = score)
  }))
}

#' Compare state taxonomies on held-out trajectory blocks
#'
#' Every candidate uses the same observations, controls and validation blocks.
#' All starts are fitted on training data only. The highest training likelihood
#' among numerically eligible starts is retained before applying the taxonomy.
#' The full-data fit supplies the labels. Held-out scores integrate over states
#' and start each independent segment at the fitted stationary distribution.
#'
#' @param data Prepared `ssf_data` (importance normalisation) or `ssf_steps`
#'   (exact normalisation). Segment IDs must be independent uninterrupted runs.
#' @param taxonomies Nonempty list of [ssf_taxonomy()] objects with distinct names.
#'   The number of labels determines each candidate's state count.
#' @param blocks Named integer vector assigning each segment to a block.
#'   Values must be consecutive integers starting at one. Default is three
#'   chronological blocks made by [ssf_blocks()]. No segment is split.
#' @param step_families Character vector of explicit positive length families,
#'   or NULL for the canonical gamma formula. Candidates cross families and
#'   taxonomies. Keep the same observation density measure throughout.
#' @param n_starts Number of independently polished starts per training fit
#'   when `initializers` is NULL.
#' @param initializers Optional function with arguments `data`, `taxonomy`,
#'   `step_family`, `seed` and `block`. It receives training rows only and
#'   returns a nonempty named list of `initial_values` for [fit_ssf()].
#'   Block zero denotes the full-data fit. It replaces `n_starts`.
#' @param short_iter,long_iter,engine,normalization Passed to [fit_ssf()].
#'   Only homogeneous transitions and lengths without lag terms are supported.
#' @param seed Nonnegative integer. Per-job seeds are recorded. The caller's
#'   random state is restored, including with parallel workers.
#' @param workers Number of processes on Unix. Other platforms run serially.
#' @param assessment_args Named arguments for [assess_taxonomy()], excluding
#'   `fit`, `taxonomy`, `convention` and `score_tolerance`.
#' @param score_tolerance,objective_tolerance Numerical eligibility thresholds.
#' @param score_failed_hessian Also report a separate point-prediction score
#'   that allows a failed Hessian check. Convergence, gradient, kernel and
#'   objective agreement checks still apply. It does not replace strict scores.
#' @return An `ssf_state_comparison` with `candidates`, `block_scores`,
#'   `segment_scores`, `starts`, full-data `fits` and `assessments`, and the
#'   supplied specifications. `start_profiles` and `start_coefficients` retain
#'   movement summaries and coefficients of every successful start, including
#'   rejected ones. `length_rank` aligns them descriptively, not biologically.
#'   A pooled score is NA unless every block was
#'   scored. `taxonomy_passed` refers to the full-data fit, not prediction.
#'   Scores are conditional log choice probabilities for importance models,
#'   or log densities for exact models. Higher is better on identical rows.
#'   This comparison does not identify a true biological number of states.
#' @export
compare_ssf_states <- function(data, taxonomies, blocks = ssf_blocks(data),
    step_families = NULL, n_starts = 5L, initializers = NULL,
    short_iter = 60L, long_iter = 400L, seed = 1L, workers = 1L,
    engine = c("native", "R"), normalization = c("importance", "exact"),
    assessment_args = list(), score_tolerance = .01, objective_tolerance = 1e-7,
    score_failed_hessian = FALSE) {
  obs <- .comparison_observations(data)
  engine <- match.arg(engine)
  normalization <- match.arg(normalization)
  if (!is.list(taxonomies) || !length(taxonomies) ||
      !all(vapply(taxonomies, inherits, logical(1), "ssf_taxonomy")))
    .stop("taxonomies must be a nonempty list of ssf_taxonomy objects.")
  names(taxonomies) <- vapply(taxonomies, `[[`, "", "name")
  if (anyDuplicated(names(taxonomies))) .stop("Taxonomy names must be distinct.")
  ids <- unique(as.character(obs$ID))
  if (!is.numeric(blocks) || is.null(names(blocks)) || anyNA(blocks) ||
      any(!is.finite(blocks)) || any(blocks < 1) || any(blocks != floor(blocks)) ||
      anyDuplicated(names(blocks)) || !setequal(names(blocks), ids) ||
      !identical(sort(unique(as.integer(blocks))), seq_len(max(blocks))))
    .stop("blocks must assign every segment once to integer blocks 1, ..., B.")
  blocks <- blocks[ids]
  n_blocks <- max(blocks)
  .scalar(n_blocks, "number of blocks", 2, length(ids), TRUE)
  for (nm in c("n_starts", "workers", "short_iter", "long_iter"))
    .scalar(get(nm), nm, 1, .Machine$integer.max, TRUE)
  .scalar(seed, "seed", 0, .Machine$integer.max, TRUE)
  .scalar(score_tolerance, "score_tolerance", 0)
  .scalar(objective_tolerance, "objective_tolerance", 0)
  if (!is.null(initializers) && !is.function(initializers))
    .stop("initializers must be NULL or a function.")
  if (!is.logical(score_failed_hessian) || length(score_failed_hessian) != 1L ||
      is.na(score_failed_hessian)) .stop("score_failed_hessian must be TRUE or FALSE.")
  if (!is.list(assessment_args) || (length(assessment_args) &&
      (is.null(names(assessment_args)) || anyDuplicated(names(assessment_args)) ||
       any(!nzchar(names(assessment_args))) || any(names(assessment_args) %in%
       c("fit", "taxonomy", "convention", "score_tolerance")))))
    .stop("assessment_args must be a named list without reserved arguments.")
  families <- if (is.null(step_families)) list(NULL) else {
    if (!is.character(step_families) || !length(step_families) ||
        anyNA(step_families) || anyDuplicated(step_families))
      .stop("step_families must be NULL or distinct family names.")
    lapply(step_families, .step_family)
  }
  grid <- expand.grid(taxonomy = names(taxonomies), family_index = seq_along(families),
    stringsAsFactors = FALSE)
  grid$family <- vapply(grid$family_index, function(i)
    if (is.null(families[[i]])) "canonical_gamma" else families[[i]], "")
  grid$candidate <- paste(grid$family, grid$taxonomy, sep = "::")
  grid$states <- vapply(taxonomies[grid$taxonomy], function(t) length(t$labels), integer(1))
  # Validate all designs before launching expensive jobs. In particular, lag
  # terms cannot be scored by the homogeneous forward calculation below.
  for (i in seq_len(nrow(grid))) {
    t <- taxonomies[[grid$taxonomy[i]]]
    f <- families[[grid$family_index[i]]]
    design <- .taxonomy_design(data, .taxonomy_formula(t, f), f, normalization == "exact")
    if (any(colnames(design$X) %in% c(.step_ar_term, .step_day_terms)))
      .stop("State comparison currently requires models without lag terms.")
  }
  jobs <- expand.grid(candidate_index = seq_len(nrow(grid)), block = 0:n_blocks)
  jobs$seed <- .with_seed(seed, function()
    sample.int(.Machine$integer.max, nrow(jobs)))
  row_block <- unname(blocks[as.character(data$ID)])
  run <- function(j) {
    job <- jobs[j, ]
    spec <- grid[job$candidate_index, ]
    tax <- taxonomies[[spec$taxonomy]]
    family <- families[[spec$family_index]]
    training <- data[row_block != job$block, , drop = FALSE]
    sources <- tryCatch(.with_seed(job$seed, function() {
      if (is.null(initializers)) stats::setNames(rep(list(NULL), n_starts),
        paste0("start_", seq_len(n_starts))) else
        initializers(data = training, taxonomy = tax, step_family = family,
          seed = job$seed, block = job$block)
    }), error = identity)
    source_error <- inherits(sources, "error")
    if (!source_error && (!is.list(sources) || !length(sources) ||
        is.null(names(sources)) || anyDuplicated(names(sources)) ||
        anyNA(names(sources)) || any(!nzchar(names(sources)))))
      sources <- simpleError("initializers must return a nonempty named list.")
    source_error <- inherits(sources, "error")
    if (source_error) sources <- list(initializers = sources)
    seeds <- .with_seed(job$seed, function() sample.int(.Machine$integer.max, length(sources)))
    fits <- lapply(seq_along(sources), function(i) {
      if (source_error) return(sources[[i]])
      tryCatch(fit_ssf(training, .taxonomy_formula(tax, family),
        n_states = spec$states, n_starts = 1L, short_iter = short_iter,
        long_iter = long_iter, seed = seeds[i], engine = engine,
        normalization = normalization, step_family = family,
        initial_values = sources[[i]]), error = identity)
    })
    audit <- do.call(rbind, lapply(seq_along(fits), function(i) {
      fit <- fits[[i]]
      failed <- inherits(fit, "error")
      d <- if (failed) list() else fit$diagnostics
      ok <- if (failed) c(point = FALSE, strict = FALSE) else
        .comparison_screen(fit, score_tolerance, objective_tolerance)
      data.frame(candidate = spec$candidate, block = job$block, start = names(sources)[i],
        seed = seeds[i], job_seed = job$seed,
        training_steps = nrow(.comparison_observations(training)),
        log_likelihood = if (failed) NA_real_ else as.numeric(stats::logLik(fit)),
        converged = isTRUE(d$converged), hessian_positive = isTRUE(d$hessian_positive),
        max_abs_gradient = if (failed) NA_real_ else d$max_abs_gradient,
        objective_agreement = if (failed) NA_real_ else d$objective_agreement,
        max_abs_scaled_gradient = if (failed) NA_real_ else .assessment_score(fit),
        curvature_relative_min = if (failed) NA_real_ else .ssf_curvature(fit$hessian)$relative_min,
        curvature_condition = if (failed) NA_real_ else .ssf_curvature(fit$hessian)$condition,
        min_occupancy = if (failed) NA_real_ else min(tabulate(fit$state, nbins = fit$n_states))/length(fit$state),
        strict = ok[["strict"]], point = ok[["point"]],
        error = if (failed) conditionMessage(fit) else "")
    }))
    profiles <- coefficients <- list()
    for (i in seq_along(fits)) {
      f <- fits[[i]]
      if (inherits(f, "error")) next
      m <- movement_metrics(stats::coef(f), f$density_measure, f$step_family,
        if (is.null(f$step_max)) Inf else f$step_max)
      rank <- match(seq_len(spec$states), order(m$mean_step, na.last = TRUE))
      rank[!is.finite(m$mean_step)] <- NA_integer_
      profiles[[length(profiles) + 1L]] <- data.frame(candidate = spec$candidate,
        block = job$block, start = names(sources)[i], state = seq_len(spec$states),
        length_rank = rank, kernel_valid = m$valid, mean_step = m$mean_step, persistence = m$persistence,
        occupancy = tabulate(f$state, nbins = spec$states)/length(f$state))
      b <- stats::coef(f)
      coefficients[[length(coefficients) + 1L]] <- data.frame(candidate = spec$candidate,
        block = job$block, start = names(sources)[i],
        state = rep(seq_len(spec$states), each = nrow(b)),
        length_rank = rep(rank, each = nrow(b)), term = rep(rownames(b), spec$states),
        estimate = as.vector(b))
    }
    choose <- function(flag) if (any(flag))
      which.max(ifelse(flag, audit$log_likelihood, -Inf)) else NA_integer_
    winner <- choose(audit$strict)
    full <- NULL
    assessment <- NULL
    if (job$block == 0L && !is.na(winner)) {
      full <- fits[[winner]]
      assessment <- do.call(assess_taxonomy, c(list(fit = full, taxonomy = tax,
        convention = full$density_measure, score_tolerance = score_tolerance), assessment_args))
    }
    scores <- list()
    if (job$block > 0L) {
      for (screen in c("strict", if (score_failed_hessian) "point")) {
        pick <- if (screen == "strict") winner else choose(audit$point)
        status <- data.frame(candidate = spec$candidate, block = job$block,
          screen = screen, start = if (is.na(pick)) NA_character_ else audit$start[pick],
          hessian_positive = if (is.na(pick)) NA else audit$hessian_positive[pick],
          steps = sum(obs$ID %in% names(blocks)[blocks == job$block]),
          log_score = NA_real_, error = "No eligible training fit")
        segments <- NULL
        if (!is.na(pick)) {
          segments <- tryCatch(.comparison_score(fits[[pick]],
            data[row_block == job$block, , drop = FALSE]), error = identity)
          if (inherits(segments, "error")) {
            status$error <- conditionMessage(segments)
            segments <- NULL
          } else if (all(is.finite(segments$log_score))) {
            status$log_score <- sum(segments$log_score)
            status$error <- ""
            segments$candidate <- spec$candidate
            segments$block <- job$block
            segments$screen <- screen
          } else {
            status$error <- "Non-finite held-out score"
            segments <- NULL
          }
        }
        scores[[screen]] <- list(status = status, segments = segments)
      }
    }
    list(audit = audit, fit = full, assessment = assessment, scores = scores,
      profiles = do.call(rbind, profiles), coefficients = do.call(rbind, coefficients))
  }
  result <- if (.Platform$OS.type == "unix" && workers > 1L)
    parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = workers,
      mc.preschedule = FALSE, mc.set.seed = FALSE) else lapply(seq_len(nrow(jobs)), run)
  if (any(vapply(result, function(x) !is.list(x) || is.null(x$audit), logical(1))))
    .stop("A comparison worker failed. Retry with workers = 1 to inspect the error.")
  starts <- do.call(rbind, lapply(result, `[[`, "audit"))
  scores <- unlist(lapply(result, `[[`, "scores"), recursive = FALSE)
  block_scores <- do.call(rbind, lapply(scores, `[[`, "status"))
  segment_scores <- do.call(rbind, lapply(scores, `[[`, "segments"))
  fits <- assessments <- stats::setNames(vector("list", nrow(grid)), grid$candidate)
  candidates <- grid[, c("candidate", "taxonomy", "family", "states")]
  candidates$numerical_passed <- candidates$taxonomy_passed <- FALSE
  candidates$AIC <- candidates$log_score <- candidates$mean_log_score <- NA_real_
  candidates$blocks_scored <- 0L
  candidates$checks_failed <- "No numerically eligible full-data fit"
  if (score_failed_hessian) candidates$point_mean_log_score <- NA_real_
  for (i in seq_len(nrow(grid))) {
    full <- result[[which(jobs$candidate_index == i & jobs$block == 0L)]]
    fits[i] <- list(full$fit)
    assessments[i] <- list(full$assessment)
    if (!is.null(full$fit)) {
      candidates$numerical_passed[i] <- TRUE
      candidates$taxonomy_passed[i] <- isTRUE(full$assessment$accepted)
      candidates$AIC[i] <- stats::AIC(full$fit)
      candidates$checks_failed[i] <- paste(full$assessment$reasons, collapse = ", ")
    }
    b <- block_scores[block_scores$candidate == grid$candidate[i] &
      block_scores$screen == "strict", ]
    candidates$blocks_scored[i] <- sum(is.finite(b$log_score))
    if (nrow(b) == n_blocks && all(is.finite(b$log_score))) {
      candidates$log_score[i] <- sum(b$log_score)
      candidates$mean_log_score[i] <- sum(b$log_score) / sum(b$steps)
    }
    if (score_failed_hessian) {
      b <- block_scores[block_scores$candidate == grid$candidate[i] &
        block_scores$screen == "point", ]
      if (nrow(b) == n_blocks && all(is.finite(b$log_score)))
        candidates$point_mean_log_score[i] <- sum(b$log_score) / sum(b$steps)
    }
  }
  structure(list(candidates = candidates, block_scores = block_scores,
    segment_scores = segment_scores, starts = starts,
    start_profiles = do.call(rbind, lapply(result, `[[`, "profiles")),
    start_coefficients = do.call(rbind, lapply(result, `[[`, "coefficients")), fits = fits,
    assessments = assessments, taxonomies = taxonomies, blocks = blocks,
    observations = obs, normalization = normalization, seed = seed,
    score_tolerance = score_tolerance, objective_tolerance = objective_tolerance,
    call = match.call()), class = "ssf_state_comparison")
}

#' @export
print.ssf_state_comparison <- function(x, ...) {
  cat("HMM-SSF state comparison\n")
  cat(nrow(x$observations), "steps |", length(unique(x$blocks)), "held-out blocks\n")
  cat("Higher log scores predict the held-out endpoints better.\n")
  print(x$candidates[, c("candidate", "states", "numerical_passed", "taxonomy_passed",
    "blocks_scored", "mean_log_score")], row.names = FALSE)
  cat("Pooled scores require every block. Labels use full-data fits.\n")
  invisible(x)
}

#' @export
summary.ssf_state_comparison <- function(object, ...) {
  structure(list(candidates = object$candidates, block_scores = object$block_scores,
    starts = object$starts, n_blocks = length(unique(object$blocks))),
    class = "summary_ssf_state_comparison")
}

#' @export
print.summary_ssf_state_comparison <- function(x, ...) {
  cat("Full-data checks and held-out scores\n")
  print(x$candidates, row.names = FALSE)
  cat("\nScores by reserved block\n")
  print(x$block_scores, row.names = FALSE)
  invisible(x)
}

#' Compare step labels between two state counts
#'
#' Each label in `to` is mapped to the label in `from` sharing the most
#' decoded steps. Several fine labels may share a parent. Ties use the label
#' order in the source taxonomy. This descriptive mapping is not a biological
#' validation or a measure of held-out classification accuracy.
#' @param comparison An [compare_ssf_states()] result.
#' @param from,to Candidate identifiers from `comparison$candidates$candidate`.
#'   Both need numerically eligible full-data fits with assignable labels.
#' @param threshold Flag steps whose largest smoothed probability in `to`
#'   is below this display threshold. It is not a calibrated error probability.
#' @return A list with `steps`, a `crosswalk` table, the majority `mapping`,
#'   and full-data taxonomy status. `.changed` compares the source label
#'   with the mapped parent. `.ambiguous` describes the destination fit only.
#' @export
compare_ssf_state_labels <- function(comparison, from, to, threshold = .7) {
  if (!inherits(comparison, "ssf_state_comparison")) .stop("Expected an ssf_state_comparison.")
  .scalar(threshold, "threshold", 0, 1)
  keys <- comparison$candidates$candidate
  if (!is.character(from) || length(from) != 1L || !from %in% keys ||
      !is.character(to) || length(to) != 1L || !to %in% keys || from == to)
    .stop("from and to must identify two distinct candidates.")
  a <- comparison$fits[[from]]
  b <- comparison$fits[[to]]
  aa <- comparison$assessments[[from]]
  ab <- comparison$assessments[[to]]
  if (is.null(a) || is.null(b) || anyNA(aa$labels) || anyNA(ab$labels))
    .stop("Both candidates need numerically eligible fits with assignable labels.")
  if (!identical(a$observations, b$observations))
    .stop("Candidates must have identical observations in the same order.")
  ta <- comparison$taxonomies[[comparison$candidates$taxonomy[match(from, keys)]]]
  tb <- comparison$taxonomies[[comparison$candidates$taxonomy[match(to, keys)]]]
  la <- factor(aa$labels[a$state], levels = ta$labels)
  lb <- factor(ab$labels[b$state], levels = tb$labels)
  crosswalk <- table(from = la, to = lb)
  mapping <- stats::setNames(ta$labels[apply(crosswalk, 2, which.max)], tb$labels)
  steps <- as.data.frame(comparison$observations)
  steps$.state_from <- la
  steps$.state_to <- lb
  steps$.parent <- unname(mapping[as.character(lb)])
  steps$.changed <- as.character(la) != steps$.parent
  steps$.max_posterior <- apply(b$posterior, 1, max)
  steps$.ambiguous <- steps$.max_posterior < threshold
  steps$.status <- factor(ifelse(steps$.changed,
    ifelse(steps$.ambiguous, "Changed and uncertain", "Changed"),
    ifelse(steps$.ambiguous, "Uncertain", "Same parent")),
    levels = c("Same parent", "Changed", "Uncertain", "Changed and uncertain"))
  list(steps = steps, crosswalk = crosswalk, mapping = mapping,
    from = from, to = to, threshold = threshold,
    taxonomy_passed = c(from = aa$accepted, to = ab$accepted))
}

#' Plot validation scores or changes along the trajectory
#' @param x An [compare_ssf_states()] result.
#' @param type Show held-out `scores` or a `track` comparison.
#' @param reference Optional candidate ID. Scores become paired gains per step
#'   against that candidate on each identical block, including its family.
#' @param from,to,threshold Passed to [compare_ssf_state_labels()] for a map.
#' @param coords Names of origin x/y and endpoint x/y columns.
#' @param segment Optional segment IDs to display.
#' @param step_index Optional observed row indices to display on a map.
#'   These refer to the complete trajectory before selecting `segment`.
#' @param ... Reserved for future use.
#' @return A ggplot object. Failed or missing fold scores are not plotted.
#'   Inspect `block_scores` and `candidates` for coverage and failure details.
#' @export
plot.ssf_state_comparison <- function(x, type = c("scores", "track"), reference = NULL,
    from = NULL, to = NULL, threshold = .7, coords = c("x0", "y0", "x", "y"),
    segment = NULL, step_index = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) .stop("Install ggplot2 to plot comparisons.")
  type <- match.arg(type)
  .data <- NULL
  if (type == "scores") {
    d <- merge(x$block_scores[x$block_scores$screen == "strict", ], x$candidates,
      by = "candidate", suffixes = c("", "_full"))
    d$.score <- d$log_score / d$steps
    if (!is.null(reference)) {
      if (length(reference) != 1L || !reference %in% x$candidates$candidate)
        .stop("reference must be one candidate identifier.")
      base <- d[d$candidate == reference, ]
      d$.score <- d$.score - base$.score[match(d$block, base$block)]
    }
    d$block <- factor(d$block)
    ggplot2::ggplot(d, ggplot2::aes(x = .data$block, y = .data$.score,
      colour = .data$taxonomy, group = .data$candidate)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
      ggplot2::geom_line(na.rm = TRUE) + ggplot2::geom_point(size = 2.8, na.rm = TRUE) +
      ggplot2::facet_wrap(~ family) +
      ggplot2::labs(x = "Reserved block", y = if (is.null(reference))
        "Log score per step" else "Log score gain per step", colour = "Taxonomy",
        subtitle = "Each block is scored using a fit trained on the other blocks",
        caption = "Only strict scores are shown. Check the candidate table for missing blocks and taxonomy failures.") +
      theme_ssf() + ggplot2::theme(legend.title = ggplot2::element_text())
  } else {
    z <- compare_ssf_state_labels(x, from, to, threshold)
    d <- z$steps
    if (!is.null(step_index)) {
      if (!is.numeric(step_index) || !length(step_index) || anyNA(step_index) ||
          any(!is.finite(step_index)) || any(step_index != floor(step_index)) ||
          any(step_index < 1 | step_index > nrow(d)) || anyDuplicated(step_index))
        .stop("step_index must contain distinct observed row indices.")
      d <- d[sort(step_index), ]
    }
    if (length(coords) != 4L || !all(coords %in% names(d)) ||
        !all(vapply(d[coords], is.numeric, logical(1))))
      .stop("coords must name four numeric origin and endpoint columns.")
    if (!is.null(segment)) {
      if (!length(segment) || !all(segment %in% d$ID)) .stop("Unknown segment.")
      d <- d[d$ID %in% segment, ]
    }
    d$.x0 <- d[[coords[1L]]]; d$.y0 <- d[[coords[2L]]]
    d$.x <- d[[coords[3L]]]; d$.y <- d[[coords[4L]]]
    origin <- d[1L, , drop = FALSE]
    d <- d[order(d$.status), , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(x = .data$.x0, y = .data$.y0,
      xend = .data$.x, yend = .data$.y, colour = .data$.status)) +
      ggplot2::geom_segment(linewidth = .7) +
      ggplot2::geom_point(ggplot2::aes(x = .data$.x, y = .data$.y, colour = .data$.status),
        inherit.aes = FALSE, size = .8) +
      ggplot2::geom_point(data = origin, ggplot2::aes(x = .data$.x0, y = .data$.y0),
        inherit.aes = FALSE, shape = 21, fill = "white", colour = "#193541", size = 2) +
      ggplot2::coord_equal() +
      ggplot2::scale_colour_manual(values = c("Same parent" = "#c5ccce",
        "Changed" = "#d87435", "Uncertain" = "#158c8c", "Changed and uncertain" = "#815e9f"),
        drop = FALSE) + ggplot2::labs(x = coords[3L], y = coords[4L], colour = NULL,
        subtitle = sprintf("%d steps change parent label; %d have destination support below %.2f",
          sum(d$.changed), sum(d$.ambiguous), threshold)) +
      theme_ssf()
  }
}
