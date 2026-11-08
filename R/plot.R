#' Plot fitted behavioural states
#'
#' Draw a track coloured by its Viterbi labels or show the smoothed state
#' probabilities in step order. Guided fits use the ecological labels of
#' the selected taxonomy. The probability curves are conditional on the
#' fitted model and stop at segment boundaries.
#'
#' @param x An object from [fit_ssf()] or [fit_guided_ssf()].
#' @param ... For `type = "states"`, arguments passed to [plot_ssf_states()].
#'   Otherwise reserved for the generic.
#' @param type Use `"track"` for a map, `"posterior"` for state probabilities,
#'   or `"auto"` for a map when the chosen coordinate columns are available
#'   and a probability plot otherwise. `"states"` combines the ecological
#'   interpretation, observed movement distributions and estimated effects.
#' @param segment Optional segment identifiers from the fitted `ID` column.
#' @param step_range Optional pair giving the first and last step to draw,
#'   counted in fitted order after selecting segments.
#' @param coords Four numeric columns in the fitted observations for start x,
#'   start y, end x and end y, respectively. Use projected coordinates in
#'   common units.
#' @return A ggplot object.
#' @export
plot.ssf_fit <- function(x, ..., type = c("auto", "track", "posterior", "states"),
                         segment = NULL, step_range = NULL,
                         coords = c("x0", "y0", "x", "y")) {
  type <- match.arg(type)
  if (type == "states") return(plot_ssf_states(x, ..., segment = segment, step_range = step_range, coords = coords))
  .plot_ssf(x, colnames(x$coefficients), type, segment, step_range, coords)
}

#' @rdname plot.ssf_fit
#' @export
plot.guided_ssf <- function(x, ..., type = c("auto", "track", "posterior", "states"),
                            segment = NULL, step_range = NULL,
                            coords = c("x0", "y0", "x", "y")) {
  if (is.null(x$fit)) .stop("No selected fit to plot. Use summary() to inspect the candidates.")
  type <- match.arg(type)
  if (type == "states") return(plot_ssf_states(x, ..., segment = segment, step_range = step_range, coords = coords))
  .plot_ssf(x$fit, x$assessment$labels, type, segment, step_range, coords)
}

.has_ssf_coordinates <- function(steps, coords) {
  is.character(coords) && length(coords) == 4L && !anyDuplicated(coords) &&
    all(coords %in% names(steps)) &&
    all(vapply(steps[coords], function(z) is.numeric(z) && all(is.finite(z)), logical(1)))
}

.plot_ssf <- function(fit, labels, type, segment, step_range, coords, palette = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) .stop("Install ggplot2 to plot the fit.")
  observations <- fit$observations
  rows <- seq_len(nrow(observations))
  if (!is.null(segment)) {
    if (!length(segment) || !all(segment %in% observations$ID))
      .stop("segment must identify fitted trajectory segments.")
    rows <- rows[observations$ID %in% segment]
  }
  if (!is.null(step_range)) {
    if (!is.numeric(step_range) || length(step_range) != 2L ||
        any(!is.finite(step_range)) || any(step_range != floor(step_range)) ||
        step_range[1L] < 1L || step_range[2L] < step_range[1L] ||
        step_range[2L] > length(rows))
      .stop("step_range must give two ordered step numbers within the selected data.")
    rows <- rows[seq.int(step_range[1L], step_range[2L])]
  }

  steps <- observations[rows, , drop = FALSE]
  has_coords <- .has_ssf_coordinates(steps, coords)
  if (type == "auto") type <- if (has_coords) "track" else "posterior"
  if (type == "track") {
    if (!has_coords) .stop("coords must identify four finite numeric columns for step starts and ends.")
    track <- stats::setNames(as.data.frame(steps[coords]), c("x0", "y0", "x", "y"))
    track$state <- factor(labels[fit$state[rows]], levels = labels)
    start <- track[!duplicated(steps$ID), , drop = FALSE]
    p <- ggplot2::ggplot(track) +
      ggplot2::geom_segment(ggplot2::aes(x = .data$x0, y = .data$y0,
        xend = .data$x, yend = .data$y, colour = .data$state),
        linewidth = 0.8, lineend = "round") +
      ggplot2::geom_point(ggplot2::aes(x = .data$x, y = .data$y, colour = .data$state),
        size = 0.9, show.legend = FALSE) +
      ggplot2::geom_point(data = start, ggplot2::aes(x = .data$x0, y = .data$y0),
        shape = 21, fill = "white", size = 2.3) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = "Inferred states on the track", x = "Easting", y = "Northing")
  } else {
    id <- as.character(steps$ID)
    single <- !(duplicated(id) | duplicated(id, fromLast = TRUE))
    posterior <- data.frame(
      step = rep(rows, times = length(labels)),
      state = factor(rep(labels, each = length(rows)), levels = labels),
      probability = as.vector(fit$posterior[rows, , drop = FALSE]),
      segment = rep(id, times = length(labels)),
      single = rep(single, times = length(labels))
    )
    p <- ggplot2::ggplot(posterior,
      ggplot2::aes(x = .data$step, y = .data$probability, colour = .data$state)) +
      ggplot2::geom_line(data = posterior[!posterior$single, , drop = FALSE],
        ggplot2::aes(group = interaction(.data$segment, .data$state)), linewidth = 0.7) +
      ggplot2::geom_point(data = posterior[posterior$single, , drop = FALSE], size = 1.5) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      ggplot2::labs(title = "Smoothed state probabilities", x = "Step in fitted order",
                   y = "Probability")
  }
  if (is.null(palette)) palette <- c("#166484", "#B6501C", "#784B8B", "#21765A", "#9E3D60", "#526372")[seq_along(labels)]
  p + ggplot2::scale_colour_manual(values = palette, drop = FALSE,
                                   name = NULL) +
    theme_ssf() + ggplot2::theme(legend.position = "bottom")
}
#' A common white theme for movement plots
#'
#' Uses the same fonts, light grid lines and white panel, legend and facet
#' backgrounds as the package's worked examples. Add ordinary ggplot2 theme
#' settings afterwards to adjust a particular plot.
#' @return A complete ggplot2 theme. Requires the optional ggplot2 package.
#' @export
theme_ssf <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) .stop("Install ggplot2 to use theme_ssf().")
  ink <- "#193541"
  muted <- "#526773"
  ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = ink),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_rect(fill = "white", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E5E5E5", linewidth = .3),
      axis.text = ggplot2::element_text(colour = muted, size = 10),
      axis.title = ggplot2::element_text(colour = ink, size = 11),
      plot.title = ggplot2::element_text(face = "bold", size = 13,
        margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(colour = muted, size = 10,
        margin = ggplot2::margin(b = 8), lineheight = 1.1),
      plot.caption = ggplot2::element_text(colour = muted, size = 9,
        hjust = 0, margin = ggplot2::margin(t = 9)),
      plot.title.position = "plot", plot.caption.position = "plot",
      legend.position = "top", legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(9, 11, 8, 11))
}

