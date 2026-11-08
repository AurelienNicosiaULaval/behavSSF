#' Read the movement signature of each inferred state
#'
#' Counts describe decoded steps. Length means and SDs describe the intrinsic
#' movement kernel, holding environmental terms constant. They are baseline
#' summaries when lengths have hourly or lag effects. Stay probabilities are
#' model estimates. A geometric mean run length is returned only for
#' homogeneous multi-state transitions; it is not an observed bout length.
#'
#' @param x A fitted `ssf_fit` or selected `guided_ssf` object.
#' @param step_scale Positive multiplier applied to step lengths and their SDs.
#'   For input in kilometres, use 1000 to display metres.
#' @param step_unit Label for the resulting distance unit. Units are not inferred.
#' @param assessment Optional accepted assessment from [assess_taxonomy()] for
#'   this fitted model. This makes a post-fit ecological interpretation visible
#'   without fitting the model again.
#' @return A data frame with state labels, decoded counts and fractions, fitted
#'   movement summaries, stay-probability ranges and expected run lengths in steps.
#' @export
ssf_state_summary <- function(x, step_scale = 1, step_unit = "input units", assessment = NULL) {
  resolved <- .state_display_fit(x, assessment)
  fit <- resolved$fit
  .scalar(step_scale, "step_scale", .Machine$double.eps)
  if (!is.character(step_unit) || length(step_unit) != 1L || is.na(step_unit) || !nzchar(step_unit))
    .stop("step_unit must be a nonempty unit label.")
  out <- .summarise_ssf(fit, resolved$labels)$states
  out$sd_step <- NA_real_
  if ("mean_step" %in% names(out)) {
    metrics <- movement_metrics(stats::coef(fit), fit$density_measure, fit$step_family,
      if (is.null(fit$step_max)) Inf else fit$step_max)
    out$mean_step <- metrics$mean_step * step_scale
    out$sd_step <- metrics$sd_step * step_scale
  } else {
    out$mean_step <- NA_real_
    out$persistence <- NA_real_
  }
  stay <- if (length(dim(fit$transition)) == 2L) {
    rbind(diag(fit$transition), diag(fit$transition))
  } else vapply(seq_len(fit$n_states), function(i) range(fit$transition[i, i, ]), numeric(2))
  out$stay_min <- stay[1L, ]
  out$stay_max <- stay[2L, ]
  out$mean_run_steps <- if (length(dim(fit$transition)) == 2L && fit$n_states > 1L)
    1 / (1 - out$stay_min) else NA_real_
  attr(out, "step_unit") <- step_unit
  out
}

.state_display_fit <- function(x, assessment = NULL) {
  fit <- if (inherits(x, "guided_ssf")) x$fit else x
  if (!is.null(assessment)) {
    if (!inherits(fit, "ssf_fit") || !inherits(assessment, "ssf_assessment") ||
        !isTRUE(assessment$accepted) || is.null(assessment$context) ||
        !isTRUE(all.equal(unname(assessment$context$coefficients), unname(stats::coef(fit)), tolerance = 0)) ||
        !isTRUE(all.equal(unname(assessment$context$transition), unname(fit$transition), tolerance = 0)) ||
        !isTRUE(all.equal(unname(assessment$occupancy), tabulate(fit$state, fit$n_states) / length(fit$state))))
      .stop("assessment must be an accepted assessment of this same fitted model.")
    return(list(fit = fit, labels = assessment$labels))
  }
  if (inherits(x, "guided_ssf")) {
    if (is.null(x$fit)) .stop("No selected fit to summarise. Use summary() to inspect the candidate decisions.")
    return(list(fit = x$fit, labels = x$assessment$labels))
  }
  if (!inherits(x, "ssf_fit")) .stop("x must be an ssf_fit or guided_ssf object.")
  list(fit = x, labels = colnames(x$coefficients))
}

#' Show how inferred states acquire an ecological interpretation
#'
#' Combines a labelled state table, a track, observed movement distributions
#' and coefficient intervals. The distributions group the observations by
#' their Viterbi labels; they are descriptive, not independent validation.
#' The same state colours are used throughout. Coefficient intervals condition
#' on the fitted model and labels and do not adjust for model or taxonomy selection.
#'
#' @inheritParams ssf_state_summary
#' @param panel `"overview"` combines the views. Use `"table"`, `"track"`,
#'   `"distributions"` or `"coefficients"` to draw them separately.
#' @param segment,step_range,coords Passed to [plot.ssf_fit()] for the map only.
#'   All other panels use the complete fitted observations.
#' @param terms Optional coefficient names for the coefficient panel. By default
#'   it shows angular, target, environmental and lag terms, excluding the
#'   distance-distribution parameters.
#' @param term_labels Optional named character vector of readable coefficient
#'   names, indexed by the names in `terms`.
#' @param angle_bandwidth Circular smoothing bandwidth in radians, shared by
#'   all states. The default is 0.35 (about 20 degrees).
#' @param coord_unit Unit label for the original map coordinates. Coordinates
#'   are not rescaled by `step_scale`.
#' @param colours Optional colour vector, named by state label or in fitted
#'   state order. This can keep the same ecological colours across models.
#' @return A ggplot object, or a patchwork for combined views. ggplot2 is required;
#'   combined views also require patchwork. The turning-angle panel requires
#'   ggcircular. Missing observations or covariance
#'   are identified in the relevant panel.
#' @export
plot_ssf_states <- function(x, panel = c("overview", "table", "track", "distributions", "coefficients"),
                            step_scale = 1, step_unit = "input units", segment = NULL,
                            step_range = NULL, coords = c("x0", "y0", "x", "y"), terms = NULL,
                            assessment = NULL, coord_unit = "input units", colours = NULL,
                            term_labels = NULL, angle_bandwidth = .35) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) .stop("Install ggplot2 to draw the state summary.")
  panel <- match.arg(panel)
  resolved <- .state_display_fit(x, assessment)
  fit <- resolved$fit
  labels <- resolved$labels
  table <- ssf_state_summary(x, step_scale, step_unit, assessment)
  palette <- stats::setNames(c("#166484", "#B6501C", "#784B8B", "#21765A", "#9E3D60", "#526372")[seq_along(labels)], labels)
  if (!is.null(colours)) {
    if (!is.character(colours) || anyNA(colours) ||
        (is.null(names(colours)) && length(colours) != length(labels)) ||
        (!is.null(names(colours)) && !all(labels %in% names(colours))))
      .stop("colours must supply a colour for every fitted state label.")
    if (is.null(names(colours))) names(colours) <- labels
    palette <- colours[labels]
    tryCatch(grDevices::col2rgb(palette), error = function(e) .stop("Invalid state colour."))
  }
  if (!is.character(coord_unit) || length(coord_unit) != 1L || is.na(coord_unit) || !nzchar(coord_unit))
    .stop("coord_unit must be a nonempty unit label.")
  ink <- "#193541"
  muted <- "#526773"
  number <- function(a) vapply(a, function(v) if (is.finite(v)) format(signif(v, 3), trim = TRUE,
    scientific = FALSE, big.mark = ",") else "Not estimated", "")
  table$.y <- rev(seq_len(nrow(table)))
  table$.label <- vapply(table$state, function(s) paste(strwrap(s, 23), collapse = "\n"), "")
  table$.count <- sprintf("%s (%.1f%%)", format(table$n_steps, big.mark = ",", trim = TRUE), 100 * table$fraction)
  table$.length <- number(table$mean_step)
  table$.sd <- paste("SD", number(table$sd_step))
  table$.stay <- if (length(dim(fit$transition)) == 2L) sprintf("%.1f%%", 100 * table$stay_min) else
    sprintf("%.0f%%-%.0f%%", 100 * table$stay_min, 100 * table$stay_max)
  table$.run <- if (length(dim(fit$transition)) > 2L) "Range across hours" else if (fit$n_states == 1L)
    "One-state model" else paste(number(table$mean_run_steps), "steps per run")
  headings <- data.frame(x = c(.2, 2.8, 5.5, 7.8, 10.3),
    label = c("State", "Observed steps", paste0("Mean step (", step_unit, ")"),
      "Turn coefficient", "Stay in state"))
  tab <- ggplot2::ggplot(table) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = .data$.y - .47), colour = "#E5E5E5", linewidth = .35) +
    ggplot2::geom_text(data = headings, ggplot2::aes(x = .data$x, label = .data$label),
      y = nrow(table) + .7, inherit.aes = FALSE, hjust = 0, size = 3.6, fontface = "bold", colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = .2, y = .data$.y, label = .data$.label, colour = .data$state),
      hjust = 0, size = 4, fontface = "bold", lineheight = .95) +
    ggplot2::geom_text(ggplot2::aes(x = 2.8, y = .data$.y + .12, label = .data$.count), hjust = 0, size = 3.7, colour = ink) +
    ggplot2::geom_segment(ggplot2::aes(x = 2.8, xend = 4.7, y = .data$.y - .17, yend = .data$.y - .17), colour = "#EEEEEE", linewidth = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = 2.8, xend = 2.8 + 1.9 * .data$fraction,
      y = .data$.y - .17, yend = .data$.y - .17, colour = .data$state), linewidth = 3) +
    ggplot2::geom_text(ggplot2::aes(x = 5.5, y = .data$.y + .18, label = .data$.length), hjust = 0, size = 4, colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = 5.5, y = .data$.y - .2, label = .data$.sd), hjust = 0, size = 3.2, colour = muted) +
    ggplot2::geom_text(ggplot2::aes(x = 7.8, y = .data$.y, label = number(.data$persistence)), hjust = 0, size = 4, colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = 10.3, y = .data$.y + .18, label = .data$.stay), hjust = 0, size = 4, colour = ink) +
    ggplot2::geom_text(ggplot2::aes(x = 10.3, y = .data$.y - .2, label = .data$.run), hjust = 0, size = 3.1, colour = muted) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::coord_cartesian(xlim = c(0, 12.2), ylim = c(.35, nrow(table) + 1), clip = "off") +
    theme_ssf() + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), legend.position = "none")
  if (panel == "table") return(tab)
  observations <- fit$observations
  observations$.state <- factor(labels[fit$state], levels = labels)
  empty <- function(title, text) ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = text, colour = muted, size = 4) +
    ggplot2::labs(title = title) + theme_ssf() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  step_plot <- empty("Observed step-length densities", "Positive step lengths are not recorded")
  if ("step" %in% names(observations) && all(is.finite(observations$step) & observations$step > 0)) {
    d <- data.frame(log_distance = log10(observations$step * step_scale), state = observations$.state)
    enough <- names(which(table(d$state) >= 2L))
    d <- d[d$state %in% enough, , drop = FALSE]
    if (nrow(d) >= 2L) {
      bandwidth <- stats::bw.nrd0(d$log_distance)
      step_plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data$log_distance,
        colour = .data$state, fill = .data$state)) +
        ggplot2::geom_density(bw = bandwidth, alpha = .12, linewidth = .85, n = 512) +
        ggplot2::scale_colour_manual(values = palette, name = NULL, drop = FALSE) +
        ggplot2::scale_fill_manual(values = palette, name = NULL, drop = FALSE) +
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = if (max(nchar(labels)) > 20) length(labels) else 1L),
          fill = ggplot2::guide_legend(nrow = if (max(nchar(labels)) > 20) length(labels) else 1L)) +
        ggplot2::scale_x_continuous(labels = function(z)
          vapply(10^z, function(v) format(signif(v, 3), scientific = FALSE, trim = TRUE), "")) +
        ggplot2::labs(title = "Observed step-length densities",
          x = paste0("Step length (", step_unit, "; log scale)"),
          y = "Density per log10 distance") + theme_ssf() +
        ggplot2::theme(legend.position = "bottom")
    }
  }
  turns <- NULL
  if (panel %in% c("overview", "distributions")) {
    if (!requireNamespace("ggcircular", quietly = TRUE))
      .stop("Install ggcircular for the circular turning-angle panel.")
    .scalar(angle_bandwidth, "angle_bandwidth", .Machine$double.eps)
    turns <- empty("Observed turning-angle densities", "Signed turning angles are not recorded")
    if ("angle" %in% names(observations) && all(is.finite(observations$angle))) {
      d <- data.frame(angle = observations$angle %% (2 * pi), state = observations$.state)
      turns <- ggplot2::ggplot(d, ggplot2::aes(x = .data$angle, colour = .data$state)) +
        ggcircular::geom_circular_density(bw = angle_bandwidth, linewidth = 1) +
        ggcircular::scale_x_circular_degrees(breaks = c(0, pi / 2, pi, 3 * pi / 2),
          labels = c("0\u00b0", "+90\u00b0", "180\u00b0", "-90\u00b0")) +
        ggplot2::scale_colour_manual(values = palette, name = NULL, drop = FALSE) +
        ggplot2::scale_y_continuous(limits = c(0, NA)) +
        ggcircular::coord_circular(zero = "north", direction = "counterclockwise", clip = "off") +
        ggplot2::labs(title = "Observed turning-angle densities",
          subtitle = "0\u00b0 = straight ahead; 180\u00b0 = reversal", x = NULL, y = NULL) +
        theme_ssf() + ggplot2::theme(legend.position = "none",
          plot.margin = ggplot2::margin(12, 24, 20, 24))
    }
  }
  if (panel == "distributions") {
    if (!requireNamespace("patchwork", quietly = TRUE)) .stop("Install patchwork for combined state views.")
    return(patchwork::wrap_plots(step_plot, turns))
  }
  coefficients <- .summarise_ssf(fit, labels)$coefficients
  radial <- c("step", "log_step", unique(unlist(lapply(.step_families, names))))
  if (is.null(terms)) terms <- setdiff(rownames(stats::coef(fit)), radial)
  if (!is.character(terms) || anyNA(terms) || anyDuplicated(terms) || !all(terms %in% rownames(stats::coef(fit))))
    .stop("terms must contain distinct fitted coefficient names.")
  effects <- coefficients[coefficients$term %in% terms, ]
  if (nrow(effects)) {
    effects$state <- factor(effects$state, levels = rev(labels))
    effects$lower <- effects$estimate - stats::qnorm(.975) * effects$std_error
    effects$upper <- effects$estimate + stats::qnorm(.975) * effects$std_error
    effects$term <- factor(effects$term, levels = terms)
    displayed_terms <- gsub("_", " ", terms, fixed = TRUE)
    displayed_terms[terms == "cos_angle"] <- "Previous heading"
    displayed_terms[terms == "cos_diff_theta_to_target"] <- "Direction to target"
    names(displayed_terms) <- terms
    if (!is.null(term_labels)) {
      if (!is.character(term_labels) || is.null(names(term_labels)) || anyNA(term_labels) ||
          anyDuplicated(names(term_labels)) || !all(names(term_labels) %in% terms))
        .stop("term_labels must be named by displayed coefficient terms.")
      displayed_terms[names(term_labels)] <- term_labels
    }
    estimates <- ggplot2::ggplot(effects, ggplot2::aes(x = .data$estimate, y = .data$state, colour = .data$state)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "#AAAAAA", linewidth = .4) +
      ggplot2::geom_linerange(ggplot2::aes(xmin = .data$lower, xmax = .data$upper), linewidth = .55, na.rm = TRUE) +
      ggplot2::geom_point(size = 2.5) + ggplot2::facet_wrap(~ term, scales = "free_x", ncol = 2,
        labeller = ggplot2::labeller(term = displayed_terms)) +
      ggplot2::scale_colour_manual(values = palette) +
      ggplot2::scale_y_discrete(labels = function(z)
        vapply(z, function(s) paste(strwrap(s, 24), collapse = "\n"), "")) +
      ggplot2::labs(title = "Estimated effects by state", x = "Estimate and 95% Wald interval", y = NULL,
        caption = if (anyNA(effects$std_error)) "Some coefficient intervals are unavailable." else NULL) +
      theme_ssf() + ggplot2::theme(legend.position = "none")
  } else estimates <- empty("Estimated effects", "No angular or selection coefficient in this model")
  if (panel == "coefficients") return(estimates)
  track <- .plot_ssf(fit, labels, "auto", segment, step_range, coords, palette = palette) +
    ggplot2::labs(title = if (.has_ssf_coordinates(observations, coords)) "Where the inferred states occur" else "Smoothed state probabilities") +
    ggplot2::theme(legend.position = "none")
  if (.has_ssf_coordinates(observations, coords)) track <- track + ggplot2::labs(
    x = paste0("Easting (", coord_unit, ")"), y = paste0("Northing (", coord_unit, ")"))
  if (panel == "track") return(track + ggplot2::theme(legend.position = "bottom"))
  if (!requireNamespace("patchwork", quietly = TRUE)) .stop("Install patchwork for the combined ecological state overview.")
  caption <- paste(
    if (inherits(x, "guided_ssf") || !is.null(assessment)) "Labels follow the chosen ecological taxonomy. Distributions group observations by their decoded state." else
      "State indices require an ecological interpretation. Distributions group observations by their decoded state.",
    "Overlaid densities use a shared bandwidth within each variable; distance densities are on a log10 scale.",
    "Means and SDs describe the intrinsic kernel. Circular angles are relative to the previous heading; radius represents density.",
    "Intervals condition on the fitted model and labels. The map shows the requested excerpt; other panels use the complete fit.", sep = "\n")
  patchwork::wrap_plots(A = tab, B = track, C = step_plot, D = turns, E = estimates,
    design = if (length(terms) > 2L) "AAA\nBCD\nEEE" else "AA\nBC\nDE",
    heights = if (length(terms) > 2L) c(1.6, 2.4, 3.2) else c(1.6, 2, 2)) +
    patchwork::plot_annotation(title = if (inherits(x, "guided_ssf") || !is.null(assessment)) "Ecological state overview" else "Fitted state overview",
      subtitle = sprintf("%d inferred states, %s observed steps%s", fit$n_states,
        format(length(fit$state), big.mark = ","), if (inherits(x, "guided_ssf") || !is.null(assessment))
          paste0(". Selected: ", paste(labels, collapse = " / ")) else ""),
      caption = caption, theme = theme_ssf())
}
