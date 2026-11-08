# Reference calculations and explicit ecological assessment

Fitting estimates the model parameters. Decoding gives a state sequence
and state probabilities for each step. The ecological checks then assess
whether the fitted states meet the chosen rules. This vignette uses
small, constructed examples to show the calculations behind the workflow
of [Bouderbala, Nicosia and Fortin
(2026)](https://doi.org/10.1371/journal.pone.0360547).

## Exact small-model decoding

For two states and three steps there are only eight possible state
paths. The package’s tests enumerate all eight to verify the likelihood,
marginal probabilities and most probable path independently of the
recursion.

``` r

transition <- matrix(c(.8, .2, .35, .65), 2, byrow = TRUE)
emissions <- matrix(c(.2, .8, .9, .1, .3, .7), 3, byrow = TRUE)
decode_states(log(emissions), transition, initial = c(.4, .6))
#> $log_likelihood
#> [1] -2.325731
#> 
#> $state
#> [1] 2 1 1
#> 
#> $posterior
#>           [,1]      [,2]
#> [1,] 0.2331750 0.7668250
#> [2,] 0.8120190 0.1879810
#> [3,] 0.5481005 0.4518995
```

Every independent segment restarts from its initial distribution. A gap
in observations must therefore create a new segment ID before fitting.

## An exactly checkable one-state fit

Each artificial stratum contains two alternatives with covariate values
zero and one. Value one is selected five times and value zero fifteen
times. With equal sampling weights and one control, the maximum is the
binomial log odds log(5/15). This provides a fitting check without
relying on another optimizer.

``` r

n_steps <- 20L
chosen <- c(rep(1, 5), rep(0, 15))
data <- data.frame(ID = 1L, stratum = rep(seq_len(n_steps), each = 2),
                   obs = rep(c(1, 0), n_steps), w = 1,
                   x = as.vector(rbind(chosen, 1 - chosen)))
fit <- fit_ssf(data, ~ x, n_states = 1L, seed = 93)
print(fit)
#> Homogeneous HMM-SSF: 1 state(s), 20 observed steps in 1 segment(s)
#> Log likelihood: -11.2467 | converged: TRUE | positive Hessian: TRUE
c(estimate = coef(fit)[1, 1], reference = log(5 / 15))
#>  estimate reference 
#> -1.098612 -1.098612
fit$starts
#>      phase start objective convergence error warnings
#> 1 baseline     0   11.2467           0               
#> 2    short     1   11.2467           0               
#> 3     long     1   11.2467           0
```

## Movement signatures

The measure convention is required. For a kernel proportional to
exp(beta_step \* r + beta_log_step \* log(r)), the gamma shape is
beta_log_step + 2 under planar area and beta_log_step + 1 under step
length and angle. The sampling density must use the same measure.
Annotate it with `prepare_ssf_data(..., density_measure = "area")` or
`"step_angle"`. The guided workflow requires this metadata.
[`convert_density_measure()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/convert_density_measure.md)
and
[`convert_movement_coefficients()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/convert_movement_coefficients.md)
convert both parts consistently.

``` r

coefficients <- rbind(step = c(-1, -2, -3),
                      log_step = c(1, 0, -.5),
                      cos_angle = c(3, 1, .2))
movement_metrics(coefficients, convention = "area")
#>   state shape rate mean_step   sd_step persistence concentration valid
#> 1     1   3.0    1       3.0 1.7320508         3.0           3.0  TRUE
#> 2     2   2.0    2       1.0 0.7071068         1.0           1.0  TRUE
#> 3     3   1.5    3       0.5 0.4082483         0.2           0.2  TRUE
#>   convention
#> 1       area
#> 2       area
#> 3       area
movement_metrics(coefficients, convention = "step_angle")
#>   state shape rate mean_step   sd_step persistence concentration valid
#> 1     1   2.0    1 2.0000000 1.4142136         3.0           3.0  TRUE
#> 2     2   1.0    2 0.5000000 0.5000000         1.0           1.0  TRUE
#> 3     3   0.5    3 0.1666667 0.2357023         0.2           0.2  TRUE
#>   convention
#> 1 step_angle
#> 2 step_angle
#> 3 step_angle
```

These summaries concern the intrinsic movement kernel. Habitat-dependent
selection can change the realised step-length distribution. The complete
zebra generator uses a gamma step density divided by step length, then
multiplies it by the angular density. Its archived weights are therefore
area densities and are coherent with the article’s +2 shape. No
additional Jacobian correction should be applied to them.

## Guided fitting on a new dataset

After preparing complete chronological case-control data, call
`fit_guided_ssf(data, convention = ...)`, optionally specifying a target
direction column. The assessment requires numerical convergence, a
positive Hessian, valid kernel summaries, occupancy support, transition
bounds and explicit movement ordering. It records why each rejected
candidate failed. The result may legitimately be `no_accepted_model`.

The HDT rule requires the entire unadjusted Wald interval for the target
coefficient to be positive. Smoothed probabilities describe the
uncertainty about individual states with the fitted parameters held
fixed. Uncertainty about the selected taxonomy would require a separate
analysis.

## Combining angular effects

The intrinsic persistence coefficient and target coefficient jointly
define the conditional angular distribution. Its concentration depends
on the target bearing relative to the previous heading.

``` r

angular_parameters(persistence = 2, target_bias = 1, relative_bearing = pi / 2)
#>   mean_angle concentration
#> 1  0.4636476      2.236068
```

## Numerical engines

`fit_ssf(..., engine = "native")` uses the compiled objective and
analytic gradient, including the derivative of the stationary initial
distribution. `engine = "R"` retains the independent R likelihood. The
tests compare them on independent segments, including singleton
sequences and up to six states. Rejected line-search trial points are
counted in the numerical diagnostics; no emission density is replaced to
accept a final model.

## Controls and exact forward simulation

Here we generate independent lengths and angles for a small example. The
control sampler records the Jacobian used to express its density per
unit area.

``` r

set.seed(317)
n <- 80L
observed <- data.frame(ID = 1L, stratum = seq_len(n),
  step = rgamma(n, shape = 2, rate = 3), angle = runif(n, -pi, pi),
  previous_bearing = 0, target_bearing = 0)
proposal <- c(shape = 1, rate = 2, mean_angle = 0, concentration = 0)
choices <- generate_ssf_controls(observed, n_controls = 20L,
                                  proposal = proposal, seed = 318)
model <- fit_ssf(choices, n_states = 1L, seed = 319)
segments <- data.frame(ID = 1L, x = 0, y = 0, bearing = 0, n_steps = 100L)
simulations <- simulate_ssf(model, segments, n_sims = 10L, seed = 320)
movement_diagnostics(simulations$tracks[[1]], segments, lags = c(1, 2, 4))
#>                    mean_step                      sd_step 
#>                   0.70566757                   0.49694835 
#>                     q95_step          short_step_fraction 
#>                   1.57404925                   0.01000000 
#>                mean_cos_turn        lag1_step_correlation 
#>                  -0.08304845                   0.03210831 
#>                short_run_q95 median_endpoint_displacement 
#>                   1.00000000                   7.37675677 
#>        q90_maximum_excursion                    msd_lag_1 
#>                  10.11641562                   0.74245481 
#>                    msd_lag_2                    msd_lag_4 
#>                   1.33308951                   2.37394231
```

The default simulator draws a fresh Markov state path. If you supply a
fixed decoded path, the simulation instead answers a conditional
question about that particular sequence. Both require correct segment
starts. For POSIXct data sampled every 30 minutes,
`prepare_ssf_data(..., time = "time", expected_interval = 1800)` rejects
gaps left inside an ID. The zebra sensitivity workflow also checks step
origins and retains short isolated sequences instead of joining them.

## Analytic normalisation and fitting without controls

For the canonical full-plane kernel, the radial integral is
Gamma(shape)/rate^shape. The angular integral is 2*pi*I0(kappa), where
kappa is the resultant concentration from persistence and any
prospective target direction. This is an exact integral for the
supported kernel, without control sampling. Extra endpoint habitat
covariates, barriers and truncation do not have this normaliser and are
rejected.

``` r

steps <- prepare_ssf_steps(observed, density_measure = "area")
exact_model <- fit_ssf(steps, n_states = 1L, normalization = "exact", seed = 321)
print(exact_model)
#> Homogeneous HMM-SSF: 1 state(s), 80 observed steps in 1 segment(s)
#> Log likelihood: -140.145 | converged: TRUE | positive Hessian: TRUE
ssf_log_normalizer(coef(exact_model), convention = "area")
#>              S1
#> [1,] -0.2980994
ssf_log_emissions(steps, ~ step + log_step + cos_angle,
  coef(exact_model), normalization = "exact")[1:4, , drop = FALSE]
#>           S1
#> 1 -1.0225299
#> 2 -5.5995137
#> 3 -0.6871202
#> 4 -0.7188064
exact_simulations <- simulate_ssf(exact_model, segments, n_sims = 10L, seed = 322)
movement_diagnostics(exact_simulations$tracks[[1]], segments, lags = c(1, 2, 4))
#>                    mean_step                      sd_step 
#>                   0.74843606                   0.50813191 
#>                     q95_step          short_step_fraction 
#>                   1.65006853                   0.02000000 
#>                mean_cos_turn        lag1_step_correlation 
#>                  -0.10920458                   0.03515177 
#>                short_run_q95 median_endpoint_displacement 
#>                   1.00000000                   4.93519086 
#>        q90_maximum_excursion                    msd_lag_1 
#>                  10.27076117                   0.81577260 
#>                    msd_lag_2                    msd_lag_4 
#>                   1.60653289                   3.03652364
```

Target effects require `relative_target_bearing`, or `target_bearing`
and `previous_bearing`, selected at the start of each step. The endpoint
cosine alone is insufficient. Exact fits optimise log rate and log shape
to enforce positivity, then return coefficients, covariance and score on
the natural scale. The Hessian’s working scale is identified separately.
The tests check analytic normalisers against independent quadrature and
the compiled score against finite differences of the R likelihood,
including zero angular resultants.

Simulation bands from these examples describe repeated tracks at the
fitted parameters and the chosen model.

## Check the numerical fit in meaningful units

The raw gradient depends on the parameter units. Assessment uses the
working gradient divided by the square root of the corresponding
observed-Hessian diagonal. Its maximum must be at most 0.01. The native
Hessian uses numerical differences of the analytic gradient. Curvature
is checked after dividing each Hessian entry by the square root of the
two diagonal entries. The smallest eigenvalue must exceed 1e-8 times the
larger of one and the largest eigenvalue. These checks survive diagonal
linear changes of parameter units. They do not prove that the optimiser
found the global maximum.

`fit$diagnostics$max_abs_scaled_gradient` reports the scaled score.
`fit$diagnostics$max_abs_gradient` retains the raw natural-scale score.
`fit$diagnostics$curvature` reports the scaled curvature. A small
gradient with weak curvature calls for closer inspection, rather than an
ecological label.

## Declare a finite availability disk

For a canonical importance model, pass `step_max` when preparing the
original case and control rows. Use the same units as `step`:

``` r

choices <- prepare_ssf_data(raw_choices, density_measure = "area", step_max = 3.8)
fit <- fit_ssf(choices, ~ step + log_step + cos_angle + canopy)
movement_metrics(coef(fit), "area", step_max = 3.8)
```

[`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md)
and [`summary()`](https://rdrr.io/r/base/summary.html) read that support
from the fit. A positive linear distance coefficient can be integrable
on this disk. The package checks the behaviour near zero and uses the
truncated radial mean for ordering. It holds habitat covariates constant
when computing that mean. This is not a method for averaging movement
over an arbitrary habitat polygon. Finite support is currently
unavailable with explicit step families or exact normalisation.

## Compare state coefficients within a rule

The rule receives labelled joint covariance as well as marginal standard
errors. For example, this additional rule asks whether the usual 95%
Wald interval for the persistence difference lies above zero:

``` r

lengths <- ssf_taxonomy("lengths", c("Short", "Long"),
  rules = function(s) c(more_persistent =
    ssf_contrast(s, "cos_angle", from = "Short", to = "Long")$lower > 0))
```

The estimate is Long minus Short. Its variance includes the covariance
between the two fitted coefficients. A missing covariance gives a
missing interval, which fails the rule. This interval conditions on the
fitted model. It does not adjust for trying several taxonomies, state
counts or solutions.
