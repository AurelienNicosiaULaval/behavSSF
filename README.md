# behavSSF <a href="https://aureliennicosiaulaval.github.io/behavSSF/"><img src="man/figures/logo.svg" align="right" width="136" alt="behavSSF hex logo"></a>

[![R package checks](https://github.com/AurelienNicosiaULaval/behavSSF/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/AurelienNicosiaULaval/behavSSF/actions/workflows/R-CMD-check.yaml) [![GitHub release](https://img.shields.io/github/v/release/AurelienNicosiaULaval/behavSSF?include_prereleases&sort=semver&label=GitHub&color=146475)](https://github.com/AurelienNicosiaULaval/behavSSF/releases) [![Documentation](https://img.shields.io/badge/docs-pkgdown-146475)](https://aureliennicosiaulaval.github.io/behavSSF/) [![GPL-3 licence](https://img.shields.io/badge/licence-GPL--3-526372)](LICENSE)

Fit hidden Markov step-selection models and map the inferred behavioural states on an animal's GPS track. Use the state-labelling rules from [Bouderbala, Nicosia and Fortin (2026)](https://doi.org/10.1371/journal.pone.0360547) or define your own, then check how the labels depend on the fitted model.

## Installation

```r
install.packages(c("remotes", "ggplot2"))
remotes::install_github("AurelienNicosiaULaval/behavSSF", ref = "v0.1.0.9010")
```

Installation from source requires an R-compatible C++ toolchain. The [GitHub release](https://github.com/AurelienNicosiaULaval/behavSSF/releases/tag/v0.1.0.9010) also provides a source archive with the rendered vignettes.

## Try it on the Hwange zebra

The package includes the prepared GPS steps from the article. This example fits the guided model, summarises the states and maps the first 48 hours of the longest segment.

```r
library(behavSSF)

zebra <- readRDS(system.file("extdata", "hwange-zebra-steps.rds",
                            package = "behavSSF"))
zebra_steps <- prepare_ssf_steps(zebra, density_measure = "area",
                               time = "time", expected_interval = 1800)
zebra_fit <- fit_guided_ssf(
  zebra_steps, convention = "area", target = "cos_diff_theta_to_target",
  normalization = "exact", n_starts = 8, short_iter = 100, long_iter = 900,
  seed = 1729
)
print(zebra_fit)
summary(zebra_fit)$states
longest <- names(which.max(table(zebra$ID)))
plot(zebra_fit, segment = longest, step_range = c(1, 96))
```

## Worked examples

- [Hwange zebra](https://aureliennicosiaulaval.github.io/behavSSF/articles/zebra-hwange.html).
  Apply the article's rules and read the state labels and probabilities.
- [Red deer](https://aureliennicosiaulaval.github.io/behavSSF/articles/deer-restarts.html).
  Compare two fitted solutions and map the steps whose labels change.
- [Andean condor](https://aureliennicosiaulaval.github.io/behavSSF/articles/condor-taxonomies.html).
  Define your own taxonomy for flight states.
- [African elephant](https://aureliennicosiaulaval.github.io/behavSSF/articles/elephant-lag-window.html).
  Compare the observed daily movement cycle with simulated tracks.

## Exact normalisation without controls

Use `prepare_ssf_steps(..., density_measure = "area")` and `fit_ssf(..., normalization = "exact")` for the canonical terms `step + log_step + cos_angle`, optionally adding `cos_diff_theta_to_target`. `ssf_log_normalizer()` exposes the analytic integral and `ssf_log_emissions(..., normalization = "exact")` evaluates the corresponding densities. The guided workflow also accepts `normalization = "exact"`.

For area measure, shape = beta_log_step + 2 and rate = -beta_step. At relative target bearing psi, persistence a and target coefficient b give concentration kappa = sqrt((a + b*cos(psi))^2 + (b*sin(psi))^2). The full-plane normaliser is Gamma(shape) / rate^shape * 2*pi*I0(kappa). With step-angle measure the shape offset is +1 instead. Both gamma parameters must be positive.

```r
library(behavSSF)
library(stats)

# Mathematical observations, not an empirical animal trajectory.
set.seed(341)
observed <- data.frame(ID = 1L, step = rgamma(120, shape = 1.6, rate = 4),
  angle = runif(120, -pi, pi))
steps <- prepare_ssf_steps(observed, density_measure = "area")
exact_fit <- fit_ssf(steps, n_states = 1L, normalization = "exact", seed = 342)
print(exact_fit)
movement_metrics(coef(exact_fit), "area")
```

For a target effect, supply `relative_target_bearing` at each step origin, or both `target_bearing` and `previous_bearing`. A supplied endpoint cosine alone does not determine the normaliser and is rejected for exact target fitting. Targets must be chosen prospectively; the constructor cannot infer from numeric bearings how they were chosen. Optional logical `target_active = FALSE` removes an undefined direction, consistently with the forward simulator's zero-distance rule.

Exact fitting uses log rate and log shape internally; `spread` applies on that working scale. Returned coefficients, covariance and score use the natural coefficient scale. The stored Hessian retains its explicitly identified working scale. Unsupported endpoint habitat terms, truncation and barriers are rejected rather than approximated. Exact integration removes control-point approximation error for these kernels; it does not remove estimation uncertainty or establish generative adequacy.

## Length autoregression

Add `lag_log_step` to the exact formula to estimate one lag coefficient per state. If the current state is j, the gamma mean at a subsequent step is `mu[j] * (previous_step / mu[j])^rho[j]`, with fixed gamma shape in that state and `abs(rho[j]) < 1`. At the start of each segment the mean is `mu[j]`. The anchor `mu[j]` also defines the length ordering in the ecological assessment. The previous length enters after state changes within a segment.

```r
# Continue the exact example with its prepared `steps`.
ar_fit <- fit_ssf(steps, ~ step + log_step + cos_angle + lag_log_step,
                  n_states = 1L, normalization = "exact", seed = 342)
coef(ar_fit)["lag_log_step", ]
```

The exact likelihood normalises the gamma distribution at each step. The prospective target term can be added alongside length autoregression. Use `step_ar = TRUE` in `fit_guided_ssf()` to apply the extension to all candidate models. `simulate_ssf()` feeds back the previously generated length. At each segment boundary both the lag and the hidden-state process restart.

## Daily movement patterns

For exact two-state fits to hourly data, supply the step-origin hours through
`transition_hour` in `fit_ssf()`. Use `step_harmonics` for a daily cycle in
mean step length and add `lag22_log_step` through `lag26_log_step` alongside
`lag_log_step` to use earlier lengths within the same segment.
`decode_ssf_diel()` decodes new observations. For simulation, supply the
first step's origin hour in the segment table's `hour` column. The
[elephant example](https://aureliennicosiaulaval.github.io/behavSSF/articles/elephant-lag-window.html)
shows how to compare daily movement patterns.

## Case-control input

The package follows the article's convention for the case-control
likelihood. The denominator includes the case and divides the sum by the
number of controls. The numerator is unweighted.

Required columns are `ID` (independent segment), `stratum` (choice set), `obs` (one case coded 1, controls coded 0), and `w` (strictly positive sampling density at every endpoint, including the case). The guided workflow also requires `step`, `log_step` and `cos_angle`. The constructor calculates the last two when positive `step` and finite `angle` in radians are supplied.

Each ID must describe an uninterrupted segment. Use separate IDs across observation gaps. Strata are chronological within a segment. Supplying a time column to the constructor verifies this ordering. The constructor preserves all complete strata, puts the case first, and records original row indices. It never drops missing values to make a fit succeed. Raw GPS preparation, projection and raster extraction remain separate preprocessing steps. The supported gamma-von Mises control generator is described below.

```r
library(behavSSF)

# A small artificial example of complete choice sets, not a movement dataset.
set.seed(19)
n_steps <- 40L
n_controls <- 4L
step <- matrix(rexp(n_steps * (n_controls + 1L)), nrow = n_steps)
angle <- matrix(runif(length(step), -pi, pi), nrow = n_steps)
data <- do.call(rbind, lapply(seq_len(n_steps), function(i) {
  data.frame(ID = 1L, stratum = i, obs = c(1L, rep(0L, n_controls)),
             w = dexp(step[i, ]) / (2 * pi),
             step = step[i, ], angle = angle[i, ])
}))
data <- prepare_ssf_data(data, density_measure = "step_angle")
fit <- fit_ssf(data, n_states = 1L, seed = 42)
print(fit)
coef(fit)
```

## Kernel convention

There is no default convention for the movement summaries or guided workflow. Choose explicitly:

* `"area"`: planar endpoint density, gamma shape `beta_log_step + 2`, matching the article's movement-summary formula.
* `"step_angle"`: step-length/angle density, gamma shape `beta_log_step + 1`.

The convention must agree with the measure used for `w`. Annotate it using `prepare_ssf_data(..., density_measure = "area")` or `"step_angle"`. The guided workflow requires this metadata and converts densities to the requested convention before fitting. `convert_density_measure()` and `convert_movement_coefficients()` perform the paired Jacobian and coefficient conversions. Changing the convention only in `movement_metrics()` does not perform these conversions.

The full generator confirms that the archived zebra weights already contain `dgamma(r, ...)/r` multiplied by the angular density. All 144,753 saved weights match this area density exactly. The article's +2 gamma shape is coherent. A second division by r would be wrong. The initial audit's unresolved warning resulted from not reading this helper; it is corrected here. A nonpositive gamma shape or rate is still flagged for any newly fitted model.

With both persistence and target direction, the conditional angular concentration is the norm of their combined vectors. Use `angular_parameters()`; `abs(beta_cos_angle)` is only the intrinsic persistence proxy used by the reference taxonomy.

```r
library(behavSSF)
coefficients <- rbind(step = c(-1, -2, -3),
                      log_step = c(1, 0, -.5),
                      cos_angle = c(3, 1, .2))
movement_metrics(coefficients, convention = "area")

# For actual data, only after checking sampling-density and kernel conventions:
# selected <- fit_guided_ssf(data, convention = "area", seed = 42)
# print(selected)
```

## Assessment and failure reporting

The guided search tries HDT/F/R when a target column is supplied, then T/F/R, T/E, and a one-state fallback. It retains all attempted models, optimizer logs and rejection reasons. If even the fallback fails its numerical and kernel checks, the result is `no_accepted_model`. Numerical failure is not evidence of homogeneous behaviour.

Relative to the archived application, this implementation resumes the long run from the fitted short-run parameters, resets decoding at segment boundaries, applies the stated minimum occupancy, requires positive target attraction for HDT, checks strict mean ordering, and requires convergence and a positive Hessian. These are documented differences; exact equality of selected models after refitting is not promised.

Smoothed state probabilities condition on a fitted parameter vector. Wald intervals used for HDT do not account for candidate selection. Neither constitutes calibrated uncertainty about the final behavioural taxonomy.

## Verification and limits

The test suite checks the recursions against exhaustive path enumeration, a closed-form one-state likelihood, analytic gradients against independent finite differences, coordinate-measure equivalence, angular moments, segment boundaries, extreme finite predictors, invalid inputs, state permutations and ecological rule failures. See the separate audit for the frozen-parameter zebra comparison and executed checks.

Use homogeneous transitions by default. Exact two-state fits on hourly
data can also model daily transitions, daily mean step lengths and lags from
22 to 26 hours alongside the one-hour lag. The [elephant example](https://aureliennicosiaulaval.github.io/behavSSF/articles/elephant-lag-window.html)
compares the height and width of the daily peak on the same test quarter.
The simulated peak remains lower and narrower than the GPS peak.

The worked examples illustrate the package on specific tracks. Check model adequacy
when applying it to a new trajectory.

## Credits and provenance

Ilhem Bouderbala wrote the original research functions. Aurélien Nicosia
reworked them for the R package and added plots, ways to present fitted
results, and support for new state taxonomies. `behavSSF` carries over the short-run and long-run fitting strategy and the published state labelling rules. It has its own likelihood and decoding routines and adds exact normalisation for supported kernels, configurable taxonomies and diagnostics.

The methodological references are the article and the [Figshare archive, version 2](https://doi.org/10.6084/m9.figshare.34077912.v2), by Ilhem Bouderbala, Aurélien Nicosia and Daniel Fortin. The Figshare archive is distributed under CC BY 4.0. `citation("behavSSF")` gives separate entries for the package and the article.

The original research scripts call the `hmmSSF` package described by Klappstein, Thomas and Michelot (2023). [Movement Ecology 11:30](https://doi.org/10.1186/s40462-023-00392-3)

## Control-sampling sensitivity and forward simulation

`generate_ssf_controls()` samples untruncated gamma step lengths and von Mises turning angles, records area densities, and keeps alternative indices for nested control pools. The target bearing is supplied at each observed step origin. It does not filter endpoints by a spatial window.

`simulate_ssf()` samples directly from the fitted gamma and von Mises laws for the canonical movement model, including an optional direction-to-target effect. It restarts the state process and position at every supplied segment. It rejects unsupported endpoint covariates instead of approximating them silently. The default uses new Markov state paths. Supplying a fixed decoded path is a separate conditional diagnostic.

`movement_diagnostics()` compares observable quantities, preserving segment boundaries. It checks that coordinate increments agree with recorded step lengths. Simulated latent states are not treated as equivalent to observed decoded states. Simulation envelopes are conditional on fitted parameters, segment starts and segment lengths; they do not propagate parameter or model-selection uncertainty.

The temporal audit of the archived zebra reference found two extra gaps hidden by repeated counter values. Current sensitivity and generative analyses retain all 6,377 observed steps on 28 independent segments, including two singletons and two 19-step segments. The earlier 26-segment analysis is preserved as a historical reference. Use `expected_interval = 1800` with POSIXct times to reject unnoticed gaps in future 30-minute data.
