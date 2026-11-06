# Fit a hidden Markov step selection model

Uses either the published SSF importance approximation or exact
canonical kernel normalisation, with the HMM forward algorithm, through
either a compiled analytic-gradient engine or an independent R
reference. All fitted SSF and transition parameters are passed from the
best short run to the long run. Errors and warnings are retained for
every attempt. This implementation prioritises verifiability, with
cross-checks between the two engines.

## Usage

``` r
fit_ssf(
  data,
  formula = ~step + log_step + cos_angle,
  n_states = 2L,
  n_starts = 5L,
  short_iter = 30L,
  long_iter = 300L,
  spread = 0.5,
  seed = NULL,
  engine = c("native", "R"),
  initial_values = NULL,
  normalization = c("importance", "exact"),
  transition_hour = NULL,
  transition_harmonics = 2L,
  step_harmonics = 0L
)
```

## Arguments

- data:

  Case-control data accepted by
  [`prepare_ssf_data()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_data.md),
  or observed
  [`prepare_ssf_steps()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_steps.md)
  data when normalization is exact.

- formula:

  One-sided SSF formula. For exact fitting, add `lag_log_step` to
  estimate a state-specific AR(1) coefficient for the log gamma mean.
  The previous length is used within a segment, including across state
  changes. Add any of `lag22_log_step` through `lag26_log_step` with
  `lag_log_step` and `transition_hour` to include those hourly lags. A
  lag only contributes when its step is available in the same segment.
  The absolute lag coefficients sum to less than one in each state.

- n_states:

  Number of hidden states, between 1 and 6.

- n_starts:

  Number of short runs. A one-state model uses one run.

- short_iter, long_iter:

  Maximum BFGS iterations for short and long phases.

- spread:

  Positive standard deviation of initial perturbations. A scalar or
  named vector with one entry for each SSF design column.

- seed:

  Optional integer seed, restored to the previous RNG state on exit.

- engine:

  `"native"` uses the compiled objective and analytic gradient; `"R"`
  uses the independent R objective with numerical derivatives.

- initial_values:

  Optional list with a coefficient matrix and positive transition
  matrix, named `coefficients` and `transition`. Included as the first
  short-run start; remaining starts are random perturbations. For daily
  transitions, supply `transition_coefficients` from a previous fit with
  the same transition harmonics in place of `transition`. Supply
  `step_hour_coefficients` from a fit with the same length harmonics to
  reuse those effects. They start at zero by default.

- normalization:

  importance uses case-control data; exact requires observed data from
  prepare_ssf_steps() and supports only the canonical gamma-von Mises
  terms and optional prospective target direction. Exact fitting
  optimises log rate and log shape to preserve positive gamma
  parameters. AR coefficients use a hyperbolic tangent link. Spread
  applies on the working scale. Reported coefficients, covariance and
  score are on the natural coefficient scale.

- transition_hour:

  Optional integer hour (0 to 23) at the beginning of each step. Daily
  transitions require two states and exact fitting on hourly
  observations. Supply hours in a consistent time zone.

- transition_harmonics:

  Number of sine/cosine pairs, from 1 to 5, when `transition_hour` is
  supplied. Each switch direction has its own logit regression. Each
  segment starts in the periodic stationary distribution for its initial
  hour.

- step_harmonics:

  Number of sine/cosine pairs in the log mean of the gamma step length,
  from 0 to 5. A positive value requires `transition_hour`. Each state
  has its own curve. The curve adds to the log mean after any AR(1)
  adjustment. The gamma shape remains fixed within each state.

## Value

An `ssf_fit` object containing coefficients, transition matrix, log
likelihood, covariance when the Hessian is positive definite, decoded
states, smoothed probabilities, diagnostics, and the full start log.
With daily transitions, `transition` is a 2 by 2 by 24 array indexed by
step-origin hour. `transition_coefficients` contains the two switching
regressions and `initial` gives state probabilities by hour. If
`step_harmonics` is positive, `step_hour_coefficients` holds the hourly
log-mean coefficients, one column per state. The movement coefficients
describe the baseline before the hourly adjustment. Use
[`decode_ssf_diel()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/decode_ssf_diel.md)
to score new observations. After BFGS converges, a local Newton
correction can reduce the score when the Hessian is positive definite.
It accepts a correction when the largest absolute score falls and the
objective improves or agrees to rounding precision. The Hessian is
recomputed at the retained solution. The `score_refinement` diagnostic
records the largest absolute score before this step and the number of
accepted corrections. The convergence code and the original BFGS logs
retain their optimiser meaning. For exact fits the stored Hessian uses
working log-rate/log-shape parameters; covariance and score use the
reported natural coefficient scale.
