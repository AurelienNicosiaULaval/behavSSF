# behavSSF

Fit hidden Markov step-selection models for animal movement. The package lets
you check how starting values affect the fitted states, give those states names
using your own rules, and compare simulated tracks with the GPS data.

This is a development preview. [Source code](https://github.com/AurelienNicosiaULaval/behavSSF) ·
[Download the source package](https://github.com/AurelienNicosiaULaval/behavSSF/releases/tag/v0.1.0.9003)

## Start with a worked example

- [African elephant](articles/elephant-lag-window.html). See how lagged step
  lengths change the daily peak in step-length correlation.
- [Red deer](articles/deer-restarts.html). Find the steps that change labels
  between two fitted solutions and see how that changes the share ending in forest.
- [Andean condor](articles/condor-taxonomies.html). Define a taxonomy for
  flight states and compare the fitted labels with a proxy from the GPS track.

[![The elephant’s daily correlation peak, observed and simulated](articles/elephant-lag-window_files/figure-html/correlation-1.png)](articles/elephant-lag-window.html)

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github("AurelienNicosiaULaval/behavSSF")
library(behavSSF)
```

Installation requires a C++ toolchain. The [reference workflow](articles/reference-workflow.html)
introduces the published rules and the functions that implement them.

## Credits

Aurélien Nicosia and Ilhem Bouderbala are the package authors. Ilhem wrote the
original research functions. The published method is by Ilhem Bouderbala,
Aurélien Nicosia and Daniel Fortin. [Article in PLOS ONE](https://doi.org/10.1371/journal.pone.0360547)
