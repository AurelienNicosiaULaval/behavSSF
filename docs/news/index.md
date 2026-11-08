# Changelog

## behavSSF 0.1.0.9011

- Read a state’s decoded share, movement signature and persistence
  together with
  [`ssf_state_summary()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_state_summary.md).
- Use `plot(result, type = "states")` for an ecological overview linking
  the taxonomy labels, trajectory, observed movement distributions and
  coefficient intervals.
- Draw the overview’s table, map, distributions or coefficients
  separately with
  [`plot_ssf_states()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot_ssf_states.md).
- Compare overlaid step-length densities and circular turning-angle
  densities with ggcircular; display conditional 95% intervals for each
  state and model effect.
- Put the ecological interpretation at the end of the Hwange zebra
  worked example.

## behavSSF 0.1.0.9010

- Use a shared white theme for package plots and worked examples with
  [`theme_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/theme_ssf.md).
- Reject reordered or damaged prepared choice sets before evaluating a
  likelihood.
- Check that custom taxonomy formula terms occur in the fitted model and
  give context when a rule fails.
- Complete the provenance notes for the Hwange zebra and the red deer.
- Explain local optima, proposal-density conventions and limits of the
  numerical screen.

## behavSSF 0.1.0.9009

- Assess canonical importance kernels on a common availability disk with
  `step_max`. Positive linear distance coefficients can be valid on this
  bounded support; ordering and summaries use the corresponding
  truncated moments.
- Screen stationarity and curvature after diagonal scaling. Raw
  gradients remain available. The new diagnostics do not depend on
  diagonal linear changes of parameter units.
- Pass labelled joint coefficient covariance to taxonomy rules. Use
  [`ssf_contrast()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_contrast.md)
  for coefficient estimates and between-state contrasts with Wald
  intervals.
- Retain rejected fitted restart models for later reassessment with
  `keep_rejected = TRUE`. Grouping still concerns admissible solutions
  only.
- Keep movement profiles and coefficients for every successful
  state-comparison start, including rejected starts, for direct
  stability checks.
- Report a concise assessment and expose scaled gradients in fit
  summaries.

## behavSSF 0.1.0.9008

- Compare custom state taxonomies and length families on common held-out
  blocks with
  [`compare_ssf_states()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_states.md).
- Keep complete segments together with
  [`ssf_blocks()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_blocks.md).
  Each start is fitted on training rows only, and incomplete block
  coverage has no pooled score.
- Report numerical checks, full-data taxonomy checks and predictive
  scores separately, retaining every start and failure. Scores allowing
  Hessian failures are an explicit separate sensitivity analysis.
- Map label changes across different state counts with
  [`compare_ssf_state_labels()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_state_labels.md)
  and [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  separating changed parent labels from low smoothed assignment support.
- Add a zebra vignette with common controls, gamma and Weibull fits,
  held-out scores, trajectory maps and an environmental taxonomy rule.

## behavSSF 0.1.0.9007

- Choose gamma, log-normal, inverse Gaussian, half-normal, Weibull or
  exponential step lengths with `step_family`. Fit analytic movement
  kernels or use matched controls for endpoint selection. The same
  taxonomy can label states by their fitted mean length across families.

- Summarise fitted states and guided candidate decisions with
  [`summary()`](https://rdrr.io/r/base/summary.html). Use
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) to draw the
  track, inspect smoothed state probabilities, or map steps that change
  labels between two restart solutions.

- Add a first vignette that applies the published state rules to the
  Hwange zebra, with the grassland target chosen at the start of each
  step.

- Fit combinations of lags from 22 to 26 hours in exact hourly models.
  The likelihood and simulator use each lag only within its segment.

- Add `lag24_log_step` alongside `lag_log_step` to exact hourly fits,
  with segment-aware likelihoods and forward simulation.

- Use `step_harmonics` to model the hourly mean of the gamma step length
  separately in each state, alongside daily transitions and optional
  AR(1).

- Fit exact two-state models with daily transition probabilities using
  `transition_hour` in
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md).
  Apply the same hourly schedule in
  [`decode_ssf_diel()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/decode_ssf_diel.md)
  and forward simulations.

## behavSSF 0.1.0

- Refine the score after BFGS when the Hessian supports a local
  correction. Record the score before refinement and the number of
  accepted corrections. The absolute score checks still apply to the
  resulting fit.

## behavSSF 0.0.7

- Diagnose sensitivity to starting values for a fixed taxonomy with
  [`diagnose_ssf_restarts()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/diagnose_ssf_restarts.md).
  Inspect admissible solutions and rejected attempts, compare decoded
  labels with
  [`compare_ssf_solutions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_solutions.md),
  and map the relabelled steps with
  [`plot_ssf_changes()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot_ssf_changes.md).
- Add a worked example on the red deer trajectory from `amt`.

## behavSSF 0.0.6

- Define state labels, a ranking signature, and candidate-specific
  checks with
  [`ssf_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md).
  Pass an ordered sequence of custom or built-in candidates to
  [`fit_guided_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_guided_ssf.md).

## behavSSF 0.0.5

- Add state-specific first-order autoregression in the log mean of exact
  gamma step lengths. The lag resets at each segment start, and forward
  simulation feeds back generated lengths.
- Support `lag_log_step` in exact fitting and `step_ar = TRUE` in the
  guided workflow. Prospective target effects can be estimated jointly.

## behavSSF 0.0.4

- Add observed-step input and analytic normalisers for unbounded
  gamma-von Mises kernels, with an optional prospective directional
  target.
- Add exact fitting and emissions to both numerical engines, with
  compiled analytic scores and positive log-rate/log-shape
  parameterisation.
- Return exact-fit covariance and score on the natural coefficient
  scale, with the Hessian scale explicitly identified.
- Reject unsupported endpoint covariates and missing prospective target
  context before exact fitting.

## behavSSF 0.0.3

- Add a gamma-von Mises control sampler with area densities and nested
  alternative indices.
- Add exact forward simulation for canonical movement kernels, with
  segment resets and optional nearest-target direction.
- Add observable movement diagnostics and distinguish Markov simulations
  from fixed-path conditional checks.
- Add expected-interval and coordinate/step continuity checks.
- Correct the zebra sensitivity reference to 28 segments after detecting
  two hidden temporal gaps, retaining all observed steps.

## behavSSF 0.0.2

- Confirm that the archived zebra proposal weights already include the
  polar Jacobian and are coherent with area-based gamma summaries.
- Add explicit sampling-density metadata and paired measure/parameter
  conversions.
- Add a registered compiled objective with analytic gradients, checked
  against the independent R reference.
- Add numerical score screening and diagnostics for rejected optimizer
  trial points.
- Add the resultant angular parameters when persistence and target
  effects are combined.
- Refit the complete zebra choice sets under three seeds and compare
  archived and prospective target definitions.

## behavSSF 0.0.1

- Add a homogeneous HMM-SSF reference likelihood and fitting engine.
- Add log-space forward/backward and Viterbi recursions with segment
  resets.
- Add multistart logs and continuation from fitted short-run parameters.
- Make gamma-kernel conventions explicit and reject invalid moments.
- Add ecological taxonomy assessment with visible rejection reasons.
- Add independent numerical reference tests.
