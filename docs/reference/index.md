# Package index

## Fit and assess models

- [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md)
  : Fit a hidden Markov step selection model
- [`fit_guided_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_guided_ssf.md)
  : Fit and assess the descending taxonomies of the reference workflow
- [`ssf_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md)
  : Define a configurable state taxonomy
- [`assess_taxonomy()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/assess_taxonomy.md)
  : Assess a fitted ecological taxonomy
- [`decode_states()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/decode_states.md)
  : Decode independent hidden Markov sequences in log space
- [`decode_ssf_diel()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/decode_ssf_diel.md)
  : Score and decode a fitted HMM-SSF with daily transitions
- [`movement_metrics()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/movement_metrics.md)
  : Recover movement summaries under an explicit kernel convention
- [`summary(`*`<ssf_fit>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/summary.ssf_fit.md)
  [`summary(`*`<guided_ssf>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/summary.ssf_fit.md)
  [`print(`*`<summary_ssf_fit>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/summary.ssf_fit.md)
  [`print(`*`<summary_guided_ssf>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/summary.ssf_fit.md)
  : Summarise an HMM-SSF fit
- [`plot(`*`<ssf_fit>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot.ssf_fit.md)
  [`plot(`*`<guided_ssf>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot.ssf_fit.md)
  : Plot fitted behavioural states

## Prepare steps and evaluate kernels

- [`prepare_ssf_data()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_data.md)
  : Validate and organise case-control step data
- [`prepare_ssf_steps()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/prepare_ssf_steps.md)
  : Prepare observed steps for analytically normalised movement kernels
- [`generate_ssf_controls()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/generate_ssf_controls.md)
  : Generate gamma-von Mises control steps with explicit area densities
- [`ssf_log_normalizer()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_log_normalizer.md)
  : Exact log normalisers for canonical gamma-von Mises kernels
- [`ssf_log_emissions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_log_emissions.md)
  : Evaluate importance-sampling or exact SSF emission densities
- [`angular_parameters()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/angular_parameters.md)
  : Combine directional persistence and a target-direction effect
- [`convert_density_measure()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/convert_density_measure.md)
  : Convert a sampling density between planar area and polar step
  coordinates
- [`convert_movement_coefficients()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/convert_movement_coefficients.md)
  : Convert movement-kernel coefficients to another density measure

## Simulate and diagnose

- [`simulate_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/simulate_ssf.md)
  : Simulate exact gamma-von Mises trajectories on independent segments
- [`movement_diagnostics()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/movement_diagnostics.md)
  : Compute observable movement summaries without crossing segment
  boundaries
- [`diagnose_ssf_restarts()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/diagnose_ssf_restarts.md)
  : Diagnose sensitivity to HMM-SSF starting values
- [`compare_ssf_solutions()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/compare_ssf_solutions.md)
  : Compare step labels between two restart solutions
- [`plot_ssf_changes()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot_ssf_changes.md)
  [`plot(`*`<ssf_solution_comparison>`*`)`](https://aureliennicosiaulaval.github.io/behavSSF/reference/plot_ssf_changes.md)
  : Map steps whose state labels change between restart solutions
