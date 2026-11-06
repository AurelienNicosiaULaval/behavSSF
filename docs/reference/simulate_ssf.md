# Simulate exact gamma-von Mises trajectories on independent segments

Supports the movement terms step, log_step, cos_angle and an optional
single cos-direction-to-target term. No finite control-point
approximation is used for simulation. With a target effect, the nearest
supplied fixed target is recomputed at the current position before each
step. Movement is unbounded; neither the raster extent nor observed
positions are imposed as barriers. With `lag_log_step`, the first step
has no AR adjustment. Later steps use the previous simulated length
within the same segment. Daily lags from 22 to 26 hours use earlier
simulated lengths within the same segment. Hourly effects on the gamma
mean apply from the first step.

## Usage

``` r
simulate_ssf(
  fit,
  segments,
  n_sims = 1L,
  target_points = NULL,
  state_path = NULL,
  target_term = "cos_diff_theta_to_target",
  seed = NULL
)
```

## Arguments

- fit:

  An `ssf_fit` with an explicit density measure and valid kernels.

- segments:

  Data frame with unique ID, x, y, bearing and n_steps columns. These
  give starting positions, previous headings, and step counts. With
  daily transitions, include `hour` (0 to 23) for each segment's first
  step origin. The simulator advances one hour per step.

- n_sims:

  Number of independent simulated replicates.

- target_points:

  Optional finite two-column matrix of fixed target points. Required if
  a target coefficient is nonzero. Nearest-target queries use an exact
  internal two-dimensional tree. At zero target distance its directional
  effect is set to zero.

- state_path:

  Optional fixed integer state path for all segments. NULL simulates a
  new Markov path, restarting from stationarity at each segment. A fixed
  path gives a conditional diagnostic, not a full generative check.

- target_term:

  Name of the single directional target coefficient.

- seed:

  Optional seed; the previous RNG state is restored on exit.

## Value

An `ssf_simulations` list with tracks, segment metadata, method and
seed. Track rows are generated endpoints, not the initial positions.
Parameters remain fixed; parameter uncertainty is not propagated.
