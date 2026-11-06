# Diagnose sensitivity to HMM-SSF starting values

Fit the same model from different starting values. Each run uses the
same steps and control sets, if supplied. The restart log shows which
fits pass the taxonomy checks and whether the optimiser agrees with the
forward likelihood calculation.

## Usage

``` r
diagnose_ssf_restarts(
  data,
  taxonomy,
  n_starts = 30L,
  seed = NULL,
  spread = 0.5,
  short_iter = 100L,
  long_iter = 900L,
  engine = c("native", "R"),
  normalization = if (inherits(data, "ssf_steps")) "exact" else "importance",
  starts = NULL,
  formula = NULL,
  target = "cos_diff_theta_to_target",
  assessment_args = list(),
  likelihood_tol = 1e-05,
  parameter_tol = 0.001,
  objective_tol = 1e-07
)
```

## Arguments

- data:

  Observed `ssf_steps` or case-control data for
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).
  Case-control data must have a density measure declared.

- taxonomy:

  A specification from
  [`ssf_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md)
  or a built-in name accepted by
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md).

- n_starts:

  Number of random starts, used when `starts` is NULL.

- seed:

  Optional seed for the whole search. The previous RNG state is restored
  on exit.

- spread, short_iter, long_iter, engine:

  Fitting controls for
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).

- normalization:

  Defaults to `"exact"` for `ssf_steps` and `"importance"` for
  case-control data.

- starts:

  Optional list of start specifications that replaces the default random
  starts. Each entry is a list with optional `seed`, `spread` and
  `initial_values` arguments for
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).
  An empty list entry uses the default random initialisation. Names of
  entries are optional group labels and may repeat.

- formula:

  Optional formula override. Defaults to the custom taxonomy formula or
  the built-in movement formula with a target term for HDT.

- target:

  Name of the target coefficient for the built-in HDT taxonomy.

- assessment_args:

  Optional named list of thresholds for
  [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md),
  for example `list(min_occupancy = 0.10)`.

- likelihood_tol, parameter_tol:

  Absolute grouping tolerances for log likelihood and aligned working
  parameters.

- objective_tol:

  Maximum allowed difference between the optimiser and forward log
  likelihood.

## Value

An `ssf_restarts` object with `solutions` ranked by likelihood, `states`
with movement and occupancy summaries for each solution, `restarts` with
one outcome per start, `attempts` with optimiser logs, and `fits` with
the representative `ssf_fit` objects. Solution numbers are scoped to
this search. Use
[`compare_ssf_solutions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_solutions.md)
for step labels.

## Details

The search groups fits that give the same Viterbi labels and have
similar likelihoods and parameters once their states are aligned by
label. It works from highest to lowest likelihood and compares each fit
to the existing representatives in that order. The parameter tolerance
applies to the largest absolute difference on the fitting scale,
including log transition odds relative to self-transition. Restart
counts depend on the starting values tried.

## See also

[`compare_ssf_solutions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_solutions.md),
[`plot_ssf_changes()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot_ssf_changes.md)

## Examples

``` r
gps <- readRDS(system.file("extdata", "deer-restart-steps.rds",
  package = "behavSSF"))
steps <- prepare_ssf_steps(gps, density_measure = "area")
labels <- ssf_taxonomy("step_lengths", c("Short steps", "Long steps"))
check <- diagnose_ssf_restarts(steps, labels, n_starts = 6, seed = 60200)
check$solutions
#>   solution representative_start log_likelihood delta_log_likelihood n_restarts
#> 1        1                    1       59.65201             0.000000          2
#> 2        2                    4       52.38869             7.263316          4
#>   n_changed_vs_best
#> 1                 0
#> 2               325
```
