# Behavioural states of a Hwange zebra

When does a zebra travel, forage or rest? We use the GPS track of one
plains zebra in Hwange National Park, Zimbabwe, and the rules proposed
by [Bouderbala, Nicosia and Fortin
(2026)](https://doi.org/10.1371/journal.pone.0360547). The package
includes 6,377 observed steps, each 30 minutes long. Gaps separate the
track into 28 segments.

## Fit the models

The prepared data give each step’s length and turning angle, along with
the direction to the nearest grassland cell from the start of the step.
This makes the target available before the zebra moves.

``` r

library(behavSSF)

zebra <- readRDS(system.file("extdata", "hwange-zebra-steps.rds",
                            package = "behavSSF"))
steps <- prepare_ssf_steps(zebra, density_measure = "area",
                           time = "time", expected_interval = 1800)

result <- fit_guided_ssf(
  steps, convention = "area", target = "cos_diff_theta_to_target",
  normalization = "exact", n_starts = 8, short_iter = 100, long_iter = 900,
  seed = 1729
)
result
#> Guided SSF: selected (T_F_R)
#>   HDT_F_R: rejected; failed checks: target_positive
#>   T_F_R: accepted
```

With the default taxonomy, the function first tries habitat-directed
travelling (HDT), foraging and resting. The travelling state must have
the longest mean step and more concentrated turns than the foraging
state. For the HDT label, the 95% confidence interval for the travelling
state’s grassland coefficient must lie above zero. If that model fails
its checks, the function tries travelling, foraging and resting without
a grassland term, then travelling and encamped, and finally a single
state. It stops at the first model that passes.

## Read the labels

The article reported habitat-directed travelling, foraging and resting.
Its archived code chose the nearest grassland cell using the observed
step endpoint. Here the cell is chosen at the start of the step, and the
likelihood is integrated directly. In this fit the travelling state’s
grassland interval includes zero, so the package selects travelling,
foraging and resting.

``` r

fit <- result$fit
labels <- result$assessment$labels
zebra$state <- labels[fit$state]

state_summary <- data.frame(
  State = labels,
  Steps = tabulate(fit$state, nbins = fit$n_states),
  `Mean step (m)` = round(1000 * result$assessment$metrics$mean_step),
  check.names = FALSE
)
knitr::kable(state_summary)
```

| State      | Steps | Mean step (m) |
|:-----------|------:|--------------:|
| Travelling |   805 |           730 |
| Foraging   |  2662 |           172 |
| Resting    |  2910 |            37 |

The counts come from the Viterbi sequence, which assigns one label to
each step. The means come from the fitted movement distributions. The
labels are inferred from GPS movement. `fit$posterior` gives the
probability of each state at each step, conditional on the fitted model.
Its columns follow the order in `labels`.

## See the states on the track

``` r

library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.5.2

longest <- names(which.max(table(zebra$ID)))
view <- head(zebra[zebra$ID == longest, ], 96)
view$x <- view$x - view$x0[1]
view$y <- view$y - view$y0[1]
view$x0 <- view$x0 - view$x0[1]
view$y0 <- view$y0 - view$y0[1]

ggplot(view) +
  geom_segment(aes(x = x0, y = y0, xend = x, yend = y, colour = state),
               linewidth = 1, lineend = "round") +
  geom_point(data = view[1, ], aes(x0, y0), shape = 21, fill = "white",
             size = 3, stroke = 1) +
  scale_colour_manual(values = c(Travelling = "#166484", Foraging = "#B6501C",
                                 Resting = "#784B8B"),
                      breaks = c("Travelling", "Foraging", "Resting")) +
  coord_equal() +
  labs(title = "48 hours on the zebra's track",
       subtitle = "First two days of the longest segment. The circle marks the start.",
       x = "Relative easting (km)", y = "Relative northing (km)", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title.position = "plot")
```

![Map of 48 hours of the Hwange zebra track. Longer paths are coloured
blue for travelling, shorter paths orange for foraging and purple for
resting.](zebra-hwange_files/figure-html/track-1.png)

The track and habitat raster come from the article’s [Figshare
supplement](https://doi.org/10.6084/m9.figshare.34077912.v2). The
original GPS data are also available in the [archive of Klappstein,
Thomas and Michelot](https://doi.org/10.5281/zenodo.7872602).
