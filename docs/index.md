# behavSSF

Fit hidden Markov step-selection models for animal movement. The package
lets you check how starting values affect the fitted states, give those
states names using your own rules, and compare simulated tracks with the
GPS data.

This is a development preview. [Source
code](https://github.com/AurelienNicosiaULaval/behavSSF) · [Download the
source
package](https://github.com/AurelienNicosiaULaval/behavSSF/releases/tag/v0.1.0.9006)

## Start with a worked example

- [Hwange
  zebra](https://aureliennicosiaulaval.github.io/behavSSF/articles/zebra-hwange.md).
  Apply the published rules to a GPS track and see which steps are
  labelled travelling, foraging or resting.
- [Red
  deer](https://aureliennicosiaulaval.github.io/behavSSF/articles/deer-restarts.md).
  Find the steps that change labels between two fitted solutions and see
  how that changes the share ending in forest.
- [Andean
  condor](https://aureliennicosiaulaval.github.io/behavSSF/articles/condor-taxonomies.md).
  Define a taxonomy for flight states and compare the fitted labels with
  a proxy from the GPS track.
- [African
  elephant](https://aureliennicosiaulaval.github.io/behavSSF/articles/elephant-lag-window.md).
  See how lagged step lengths change the daily peak in step-length
  correlation.

## Install from GitHub

``` r

install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/behavSSF")
library(behavSSF)
```

Installation requires a C++ toolchain. The [reference
workflow](https://aureliennicosiaulaval.github.io/behavSSF/articles/reference-workflow.md)
introduces the published rules and the functions that implement them.

## Credits

Ilhem Bouderbala wrote the original research functions. Aurélien Nicosia
reworked them for the R package and added plots, ways to present fitted
results, and support for new state taxonomies. The published method is
by Ilhem Bouderbala, Aurélien Nicosia and Daniel Fortin. [Article in
PLOS ONE](https://doi.org/10.1371/journal.pone.0360547)
