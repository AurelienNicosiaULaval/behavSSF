# Compare state taxonomies on held-out trajectory blocks

Every candidate uses the same observations, controls and validation
blocks. All starts are fitted on training data only. The highest
training likelihood among numerically eligible starts is retained before
applying the taxonomy. The full-data fit supplies the labels. Held-out
scores integrate over states and start each independent segment at the
fitted stationary distribution.

## Usage

``` r
compare_ssf_states(
  data,
  taxonomies,
  blocks = ssf_blocks(data),
  step_families = NULL,
  n_starts = 5L,
  initializers = NULL,
  short_iter = 60L,
  long_iter = 400L,
  seed = 1L,
  workers = 1L,
  engine = c("native", "R"),
  normalization = c("importance", "exact"),
  assessment_args = list(),
  score_tolerance = 0.01,
  objective_tolerance = 1e-07,
  score_failed_hessian = FALSE
)
```

## Arguments

- data:

  Prepared `ssf_data` (importance normalisation) or `ssf_steps` (exact
  normalisation). Segment IDs must be independent uninterrupted runs.

- taxonomies:

  Nonempty list of
  [`ssf_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md)
  objects with distinct names. The number of labels determines each
  candidate's state count.

- blocks:

  Named integer vector assigning each segment to a block. Values must be
  consecutive integers starting at one. Default is three chronological
  blocks made by
  [`ssf_blocks()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_blocks.md).
  No segment is split.

- step_families:

  Character vector of explicit positive length families, or NULL for the
  canonical gamma formula. Candidates cross families and taxonomies.
  Keep the same observation density measure throughout.

- n_starts:

  Number of independently polished starts per training fit when
  `initializers` is NULL.

- initializers:

  Optional function with arguments `data`, `taxonomy`, `step_family`,
  `seed` and `block`. It receives training rows only and returns a
  nonempty named list of `initial_values` for
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).
  Block zero denotes the full-data fit. It replaces `n_starts`.

- short_iter, long_iter, engine, normalization:

  Passed to
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).
  Only homogeneous transitions and lengths without lag terms are
  supported.

- seed:

  Nonnegative integer. Per-job seeds are recorded. The caller's random
  state is restored, including with parallel workers.

- workers:

  Number of processes on Unix. Other platforms run serially.

- assessment_args:

  Named arguments for
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md),
  excluding `fit`, `taxonomy`, `convention` and `score_tolerance`.

- score_tolerance, objective_tolerance:

  Numerical eligibility thresholds.

- score_failed_hessian:

  Also report a separate point-prediction score that allows a failed
  Hessian check. Convergence, gradient, kernel and objective agreement
  checks still apply. It does not replace strict scores.

## Value

An `ssf_state_comparison` with `candidates`, `block_scores`,
`segment_scores`, `starts`, full-data `fits` and `assessments`, and the
supplied specifications. `start_profiles` and `start_coefficients`
retain movement summaries and coefficients of every successful start,
including rejected ones. `length_rank` aligns them descriptively, not
biologically. A pooled score is NA unless every block was scored.
`taxonomy_passed` refers to the full-data fit, not prediction. Scores
are conditional log choice probabilities for importance models, or log
densities for exact models. Higher is better on identical rows. This
comparison does not identify a true biological number of states.
