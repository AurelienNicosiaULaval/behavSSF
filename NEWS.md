# behavSSF 0.1.0.9010

- Use a shared white theme for package plots and worked examples with `theme_ssf()`.
- Reject reordered or damaged prepared choice sets before evaluating a likelihood.
- Check that custom taxonomy formula terms occur in the fitted model and give context when a rule fails.
- Complete the provenance notes for the Hwange zebra and the red deer.
- Explain local optima, proposal-density conventions and limits of the numerical screen.

# behavSSF 0.1.0.9009

- Assess canonical importance kernels on a common availability disk with `step_max`. Positive linear distance coefficients can be valid on this bounded support; ordering and summaries use the corresponding truncated moments.
- Screen stationarity and curvature after diagonal scaling. Raw gradients remain available. The new diagnostics do not depend on diagonal linear changes of parameter units.
- Pass labelled joint coefficient covariance to taxonomy rules. Use `ssf_contrast()` for coefficient estimates and between-state contrasts with Wald intervals.
- Retain rejected fitted restart models for later reassessment with `keep_rejected = TRUE`. Grouping still concerns admissible solutions only.
- Keep movement profiles and coefficients for every successful state-comparison start, including rejected starts, for direct stability checks.
- Report a concise assessment and expose scaled gradients in fit summaries.

# behavSSF 0.1.0.9008

- Compare custom state taxonomies and length families on common held-out blocks with `compare_ssf_states()`.
- Keep complete segments together with `ssf_blocks()`. Each start is fitted on training rows only, and incomplete block coverage has no pooled score.
- Report numerical checks, full-data taxonomy checks and predictive scores separately, retaining every start and failure. Scores allowing Hessian failures are an explicit separate sensitivity analysis.
- Map label changes across different state counts with `compare_ssf_state_labels()` and `plot()`, separating changed parent labels from low smoothed assignment support.
- Add a zebra vignette with common controls, gamma and Weibull fits, held-out scores, trajectory maps and an environmental taxonomy rule.

# behavSSF 0.1.0.9007

* Choose gamma, log-normal, inverse Gaussian, half-normal, Weibull or
  exponential step lengths with `step_family`. Fit analytic movement kernels
  or use matched controls for endpoint selection. The same taxonomy can
  label states by their fitted mean length across families.

* Summarise fitted states and guided candidate decisions with `summary()`.
  Use `plot()` to draw the track, inspect smoothed state probabilities, or
  map steps that change labels between two restart solutions.

* Add a first vignette that applies the published state rules to the Hwange
  zebra, with the grassland target chosen at the start of each step.

* Fit combinations of lags from 22 to 26 hours in exact hourly models.
  The likelihood and simulator use each lag only within its segment.

* Add `lag24_log_step` alongside `lag_log_step` to exact hourly fits, with
  segment-aware likelihoods and forward simulation.
* Use `step_harmonics` to model the hourly mean of the gamma step length
  separately in each state, alongside daily transitions and optional AR(1).

* Fit exact two-state models with daily transition probabilities using
  `transition_hour` in `fit_ssf()`. Apply the same hourly schedule in
  `decode_ssf_diel()` and forward simulations.

# behavSSF 0.1.0

* Refine the score after BFGS when the Hessian supports a local correction.
  Record the score before refinement and the number of accepted corrections.
  The absolute score checks still apply to the resulting fit.

# behavSSF 0.0.7

* Diagnose sensitivity to starting values for a fixed taxonomy with
  `diagnose_ssf_restarts()`. Inspect admissible solutions and rejected attempts,
  compare decoded labels with `compare_ssf_solutions()`, and map the relabelled
  steps with `plot_ssf_changes()`.
* Add a worked example on the red deer trajectory from `amt`.

# behavSSF 0.0.6

* Define state labels, a ranking signature, and candidate-specific checks with
  `ssf_taxonomy()`. Pass an ordered sequence of custom or built-in candidates
  to `fit_guided_ssf()`.

# behavSSF 0.0.5

* Add state-specific first-order autoregression in the log mean of exact gamma step lengths. The lag resets at each segment start, and forward simulation feeds back generated lengths.
* Support `lag_log_step` in exact fitting and `step_ar = TRUE` in the guided workflow. Prospective target effects can be estimated jointly.

# behavSSF 0.0.4

* Add observed-step input and analytic normalisers for unbounded gamma-von Mises kernels, with an optional prospective directional target.
* Add exact fitting and emissions to both numerical engines, with compiled analytic scores and positive log-rate/log-shape parameterisation.
* Return exact-fit covariance and score on the natural coefficient scale, with the Hessian scale explicitly identified.
* Reject unsupported endpoint covariates and missing prospective target context before exact fitting.

# behavSSF 0.0.3

* Add a gamma-von Mises control sampler with area densities and nested alternative indices.
* Add exact forward simulation for canonical movement kernels, with segment resets and optional nearest-target direction.
* Add observable movement diagnostics and distinguish Markov simulations from fixed-path conditional checks.
* Add expected-interval and coordinate/step continuity checks.
* Correct the zebra sensitivity reference to 28 segments after detecting two hidden temporal gaps, retaining all observed steps.

# behavSSF 0.0.2

* Confirm that the archived zebra proposal weights already include the polar Jacobian and are coherent with area-based gamma summaries.
* Add explicit sampling-density metadata and paired measure/parameter conversions.
* Add a registered compiled objective with analytic gradients, checked against the independent R reference.
* Add numerical score screening and diagnostics for rejected optimizer trial points.
* Add the resultant angular parameters when persistence and target effects are combined.
* Refit the complete zebra choice sets under three seeds and compare archived and prospective target definitions.

# behavSSF 0.0.1

* Add a homogeneous HMM-SSF reference likelihood and fitting engine.
* Add log-space forward/backward and Viterbi recursions with segment resets.
* Add multistart logs and continuation from fitted short-run parameters.
* Make gamma-kernel conventions explicit and reject invalid moments.
* Add ecological taxonomy assessment with visible rejection reasons.
* Add independent numerical reference tests.
