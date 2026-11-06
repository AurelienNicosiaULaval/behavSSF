# Summarise an HMM-SSF fit

For a guided fit, the summary lists the candidates that were tried and
gives the ecological labels of the selected fit. A fit from
[`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md)
uses its numerical state names. The state table counts Viterbi
assignments. Kernel means and signed directional persistence are
included when the movement terms and their density convention are
available.

## Usage

``` r
# S3 method for class 'ssf_fit'
summary(object, ...)

# S3 method for class 'guided_ssf'
summary(object, ...)

# S3 method for class 'summary_ssf_fit'
print(x, digits = 4L, ...)

# S3 method for class 'summary_guided_ssf'
print(x, digits = 4L, ...)
```

## Arguments

- object:

  An object from
  [`fit_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_ssf.md)
  or
  [`fit_guided_ssf()`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_guided_ssf.md).

- ...:

  Reserved for the generic.

- x:

  A summary returned by one of these methods.

- digits:

  Number of significant digits to display.

## Value

A list with a state table (`states`), step-selection estimates and
standard errors (`coefficients`), transition probabilities
(`transition`) and optimisation diagnostics. Kernel means are in the
input step-length units and describe the baseline kernel when there are
lag or hourly effects. Standard errors are available when the fit has a
covariance matrix. Daily transition arrays retain the hour dimension.
Guided summaries also contain `status`, `taxonomy` and `candidates`. The
candidate table gives each code, its outcome and a description of its
states when available. When no candidate is selected, only these three
components are returned.
