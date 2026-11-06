.stop <- function(...) stop(..., call. = FALSE)

.scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x))) {
    .stop(name, " must be a finite ", if (integer) "integer " else "number ",
          "in [", lower, ", ", upper, "].")
  }
}

.with_seed <- function(seed, fun) {
  if (is.null(seed)) return(fun())
  .scalar(seed, "seed", 0, .Machine$integer.max, TRUE)
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) previous <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  fun()
}

#' Validate and organise case-control step data
#'
#' Each ID must identify an independent, uninterrupted segment. Supply separate
#' IDs across observation gaps, even for the same animal. Strata follow their
#' first appearance in the input; rows are reordered within strata to put the
#' case first. If time is supplied, case times must increase within each segment.
#' No incomplete stratum or nonfinite value is silently removed.
#'
#' @param data A data frame of observed and control steps.
#' @param id,stratum,observed,weight Names of the segment, choice set,
#'   binary observation, and positive sampling density columns.
#' @param time Optional name of a numeric or POSIXct observation time column.
#' @param density_measure Measure of the supplied sampling density: `"area"`,
#'   `"step_angle"`, or `"unspecified"`. This annotation does not alter weights.
#' @param expected_interval Optional positive interval between consecutive
#'   cases within a segment. Requires time. Units are seconds for POSIXct,
#'   otherwise the units of the numeric time column. Gaps cause an error.
#' @return A validated data frame of class `ssf_data`, with standard columns
#'   `ID`, `stratum`, `obs`, and `w`, and a `.source_row` mapping.
#' @export
#' @examples
#' d <- data.frame(ID = 1, stratum = rep(1:3, each = 2),
#'                 obs = rep(c(1, 0), 3), w = 1, step = 1:6)
#' prepare_ssf_data(d)
prepare_ssf_data <- function(data, id = "ID", stratum = "stratum",
                             observed = "obs", weight = "w", time = NULL,
                             density_measure = c("unspecified", "area", "step_angle"),
                             expected_interval = NULL) {
  density_measure <- match.arg(density_measure)
  if (!is.null(expected_interval)) {
    .scalar(expected_interval, "expected_interval", .Machine$double.eps)
    if (is.null(time)) .stop("expected_interval requires a time column.")
  }
  if (!is.data.frame(data) || nrow(data) == 0L || anyDuplicated(names(data)))
    .stop("data must be a nonempty data frame with unique column names.")
  columns <- c(id, stratum, observed, weight, time)
  if (!is.character(columns) || anyNA(columns) || anyDuplicated(columns) ||
      !all(columns %in% names(data)))
    .stop("Column names must be distinct and present in data.")
  if (".source_row" %in% names(data))
    .stop(".source_row is reserved; supply the original data frame.")
  out <- as.data.frame(data)
  key_id <- as.character(out[[id]])
  key_set <- as.character(out[[stratum]])
  if (anyNA(key_id) || anyNA(key_set) || any(!nzchar(key_id)) || any(!nzchar(key_set)))
    .stop("Segment and stratum identifiers must be nonmissing and nonempty.")
  obs <- out[[observed]]
  w <- out[[weight]]
  if (!(is.numeric(obs) || is.logical(obs)) || anyNA(obs) || any(!obs %in% 0:1))
    .stop("observed must contain only 0 and 1.")
  if (!is.numeric(w) || any(!is.finite(w)) || any(w <= 0))
    .stop("Sampling densities must be finite and strictly positive.")
  keys <- interaction(factor(key_id), factor(key_set), drop = TRUE)
  set <- match(keys, unique(keys))
  groups <- split(seq_len(nrow(out)), set)
  counts <- vapply(groups, function(i) sum(obs[i]), numeric(1))
  sizes <- lengths(groups)
  if (any(counts != 1) || any(sizes < 2L))
    .stop("Every stratum must contain exactly one case and at least one control.")
  case_rows <- vapply(groups, function(i) i[obs[i] == 1], integer(1))
  segment <- match(key_id[case_rows], unique(key_id[case_rows]))
  if (anyDuplicated(rle(segment)$values))
    .stop("Each segment must occur in one contiguous block of strata.")
  if (!is.null(time)) {
    tm <- out[[time]]
    if (!(is.numeric(tm) || inherits(tm, "POSIXct")) || any(!is.finite(as.numeric(tm))))
      .stop("time must be finite numeric values or POSIXct times.")
    if (any(vapply(split(as.numeric(tm[case_rows]), segment),
                   function(x) any(diff(x) <= 0), logical(1))))
      .stop("Case times must be strictly increasing within every segment.")
    if (!is.null(expected_interval) && any(vapply(split(as.numeric(tm[case_rows]), segment),
      function(x) any(abs(diff(x)-expected_interval) > 1e-8), logical(1))))
      .stop("A time gap occurs within a segment; split the segment at missing observations.")
    if (any(vapply(groups, function(i) length(unique(tm[i])) != 1L, logical(1))))
      .stop("All rows in a stratum must share its case time.")
  }
  # Standardise only after reading all source columns, to avoid name collisions.
  order_rows <- unlist(lapply(groups, function(i) c(i[obs[i] == 1], i[obs[i] == 0])),
                       use.names = FALSE)
  out$ID <- key_id
  out$stratum <- set
  out$obs <- as.integer(obs)
  out$w <- w
  out$.source_row <- seq_len(nrow(out))
  out <- out[order_rows, , drop = FALSE]
  rownames(out) <- NULL
  if ("step" %in% names(out) && !"log_step" %in% names(out)) {
    if (!is.numeric(out$step) || any(!is.finite(out$step)) || any(out$step <= 0))
      .stop("step must be finite and positive to calculate log_step.")
    out$log_step <- log(out$step)
  }
  if ("angle" %in% names(out) && !"cos_angle" %in% names(out)) {
    if (!is.numeric(out$angle) || any(!is.finite(out$angle)))
      .stop("angle must be finite numeric radians.")
    out$cos_angle <- cos(out$angle)
  }
  class(out) <- c("ssf_data", "data.frame")
  attr(out, "density_measure") <- density_measure
  out
}

#' Convert a sampling density between planar area and polar step coordinates
#'
#' Since dx dy = r dr dtheta, q_area = q_step_angle / r. This function changes
#' weights and their metadata together. It does not change covariates or fit a
#' model. Use [convert_movement_coefficients()] for the corresponding change
#' of kernel coefficients.
#'
#' @param data Validated `ssf_data` with an explicit density measure and
#'   finite positive `step` values, or observed ssf_steps. For ssf_steps there
#'   are no proposal weights; only the reference-measure metadata changes.
#' @param to Target measure, `"area"` or `"step_angle"`.
#' @return An `ssf_data` object with converted w and updated density metadata.
#' @export
convert_density_measure <- function(data, to) {
  to <- match.arg(to, c("area", "step_angle"))
  if (!inherits(data, "ssf_data") && !inherits(data,"ssf_steps"))
    .stop("data must be validated ssf_data or ssf_steps.")
  from <- attr(data, "density_measure")
  if (is.null(from) || !from %in% c("area", "step_angle"))
    .stop("Annotate the sampling density measure before converting it.")
  if (from == to) return(data)
  if (inherits(data,"ssf_steps")) {
    attr(data,"density_measure") <- to
    return(data)
  }
  r <- data$step
  if (!is.numeric(r) || length(r) != nrow(data) || any(!is.finite(r)) || any(r <= 0))
    .stop("Measure conversion requires finite, strictly positive step lengths.")
  data$w <- if (to == "area") data$w / r else data$w * r
  if (any(!is.finite(data$w)) || any(data$w <= 0))
    .stop("Converted densities are not finite and positive.")
  attr(data, "density_measure") <- to
  data
}

#' Convert movement-kernel coefficients to another density measure
#'
#' A radial-angle density equals the planar density multiplied by r, so its
#' log_step coefficient is larger by one. All other coefficients are unchanged.
#' This conversion preserves the underlying spatial model when proposal
#' densities are also converted with [convert_density_measure()].
#'
#' @param coefficients Numeric coefficient matrix with a unique log_step row.
#' @param from,to Source and target measures, `"area"` or `"step_angle"`.
#' @return A coefficient matrix in the target measure.
#' @export
convert_movement_coefficients <- function(coefficients, from, to) {
  from <- match.arg(from, c("area", "step_angle"))
  to <- match.arg(to, c("area", "step_angle"))
  b <- as.matrix(coefficients)
  if (!is.numeric(b) || any(!is.finite(b)) || !ncol(b) ||
      is.null(rownames(b)) || sum(rownames(b) == "log_step") != 1L)
    .stop("A finite coefficient matrix with one log_step row is required.")
  if (from != to) b["log_step", ] <- b["log_step", ] + if (to == "step_angle") 1 else -1
  b
}

.design <- function(data, formula) {
  if (!inherits(data, "ssf_data")) data <- prepare_ssf_data(data)
  if (!inherits(formula, "formula") || length(formula) != 2L)
    .stop("formula must be one-sided, for example ~ step + log_step + cos_angle.")
  mm <- stats::model.matrix(formula, data = stats::model.frame(
    formula, data = data, na.action = stats::na.pass))
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (!ncol(mm) || nrow(mm) != nrow(data) || any(!is.finite(mm)))
    .stop("The SSF design must have finite covariates and at least one column.")
  groups <- split(seq_len(nrow(data)), data$stratum)
  centered <- mm
  for (i in groups) centered[i, ] <- sweep(mm[i, , drop = FALSE], 2,
                                            colMeans(mm[i, , drop = FALSE]))
  if (qr(centered)$rank < ncol(mm))
    .stop("SSF coefficients are not identifiable within choice sets.")
  cases <- which(data$obs == 1L)
  list(data = data, X = mm, groups = groups, cases = cases,
       logw = log(data$w), segment = data$ID[cases])
}
