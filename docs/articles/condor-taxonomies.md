# Configurable taxonomies for condor flight

Can the turns in a condor’s GPS track distinguish curved flight from
directed flight? We give the two states our own names and fit the model
to horizontal steps and turns. After fitting, we compare the labels with
altitude and with the flight proxies from the study.

## GPS data

Kuzco’s GPS fixes come from the flight dataset of [Keskin et
al. (2024)](https://doi.org/10.5281/zenodo.12607008), distributed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The
associated study distinguishes thermalling and gliding using path
curvature and vertical speed. Our prepared subset contains 1,304 steps
in 17 segments across seven dates. It uses measured fixes near a
five-second clock, omits intervals whose endpoints have source ground
speed below 3 m/s, splits gaps and retains segments with at least 12
fitted steps. Coordinates are relative to a daily origin. [Keskin et
al. (2025)](https://doi.org/10.1016/j.isci.2025.112090)

``` r

library(behavSSF)
library(ggplot2)
theme_set(theme_ssf())
knitr::opts_chunk$set(dev.args = list(bg = "white"))
library(patchwork)
library(stats)

knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.align = "center",
                      fig.retina = 2, fig.width = 10.4, dev.args = list(bg = "white"))
path <- system.file("extdata", "condor-flight-steps.csv", package = "behavSSF",
                    mustWork = TRUE)
gps <- utils::read.csv(path)
steps <- prepare_ssf_steps(gps, density_measure = "area", time = "seconds",
                          expected_interval = 5)
```

## Define and fit the candidate

The two labels are assigned in increasing order of the signed angular
persistence coefficient. The candidate also requires a positive
coefficient for directed flight and a gap greater than 1 between the two
coefficients. These are example contrasts for this demonstration.
Numerical stability, Viterbi occupancy and self-transition checks apply
alongside the custom rules.

``` r

flight_modes <- ssf_taxonomy(
  name = "flight_modes",
  labels = c("Circling-like flight", "Directed flight"),
  order_by = "persistence",
  rules = function(states) {
    p <- setNames(states$metrics$persistence, rownames(states$metrics))
    c(persistence_gap = p[["Directed flight"]] - p[["Circling-like flight"]] > 1,
      forward_persistence = p[["Directed flight"]] > 0)
  }
)

result <- fit_guided_ssf(
  steps, convention = "area", normalization = "exact",
  taxonomies = list(flight_modes, "Other"),
  persistence_bounds = c(.05, .995),
  n_starts = 8, short_iter = 100, long_iter = 900, seed = 1729
)
print(result)
#> Guided HMM-SSF with 1304 observed steps
#> Candidate taxonomies in the order tried
#>   flight_modes (Circling-like flight / Directed flight) [Selected]
result$history$flight_modes$assessment$checks
#>           converged             hessian               score        valid_kernel 
#>                TRUE                TRUE                TRUE                TRUE 
#>          transition           occupancy         state_order     persistence_gap 
#>                TRUE                TRUE                TRUE                TRUE 
#> forward_persistence 
#>                TRUE
```

The search accepts the first candidate that passes. We set the upper
self-transition bound to 0.995 to allow a flight mode to persist across
many five-second steps. The turning kernel is symmetric, so left and
right turns of the same size have equal density within a state.

``` r

stopifnot(identical(result$taxonomy, "flight_modes"))
assessment <- result$assessment
fit <- result$fit
gps$state <- factor(assessment$labels[fit$state], levels = flight_modes$labels)
probabilities <- as.data.frame(fit$posterior)
names(probabilities) <- assessment$labels

state_summary <- data.frame(
  State = assessment$labels,
  Steps = tabulate(fit$state, nbins = fit$n_states),
  Share = round(100 * assessment$occupancy, 1),
  `Mean length (m)` = round(1000 * assessment$metrics$mean_step, 1),
  Persistence = round(assessment$metrics$persistence, 2),
  check.names = FALSE
)
names(state_summary)[3] <- "Share (%)"
knitr::kable(state_summary, row.names = FALSE)
```

| State                | Steps | Share (%) | Mean length (m) | Persistence |
|:---------------------|------:|----------:|----------------:|------------:|
| Circling-like flight |   936 |      71.8 |            43.5 |       -1.11 |
| Directed flight      |   368 |      28.2 |            66.7 |        4.87 |

## See the states along a continuous segment

The map shows the longest fitted segment. The altitude panel uses the
same steps and colours each endpoint by its Viterbi state. Altitude is
displayed after fitting. The probability strip shows the smoothed
probabilities conditional on the fitted model.

``` r

palette <- c("Circling-like flight" = "#13726B", "Directed flight" = "#A54C24")
ink <- "#183345"
flight_theme <- theme_ssf() + theme(legend.position = "bottom")
longest <- names(which.max(table(gps$ID)))
i <- which(gps$ID == longest)
track <- gps[i, ]
track$minute <- seq_len(nrow(track)) * 5 / 60
origin <- track[1, c("x0_m", "y0_m")]
for (name in c("x0_m", "x_m")) track[[name]] <- track[[name]] - origin$x0_m
for (name in c("y0_m", "y_m")) track[[name]] <- track[[name]] - origin$y0_m

route <- ggplot(track) +
  geom_segment(aes(x = x0_m, y = y0_m, xend = x_m, yend = y_m, colour = state),
               linewidth = .95, lineend = "round") +
  annotate("point", x = 0, y = 0, shape = 21, fill = "white", colour = ink, size = 2.6) +
  annotate("text", x = 0, y = -18, label = "Start", colour = ink, size = 3.1) +
  scale_colour_manual(values = palette, name = NULL) +
  coord_equal() +
  labs(title = "A   Horizontal route", x = "Relative east (m)", y = "Relative north (m)") +
  flight_theme

altitude <- ggplot(track, aes(minute, altitude_m)) +
  geom_line(colour = "#93A6AC", linewidth = .55) +
  geom_point(aes(colour = state), size = 1.65) +
  scale_colour_manual(values = palette, guide = "none") +
  labs(title = "B   Recorded altitude", x = NULL, y = "Altitude (m)") +
  flight_theme + theme(axis.text.x = element_blank())

strip <- do.call(rbind, lapply(names(probabilities), function(label) {
  data.frame(minute = track$minute, p = probabilities[i, label], state = label)
}))
strip$state <- factor(strip$state, levels = names(palette))
uncertainty <- ggplot(strip, aes(minute, p, fill = state)) +
  geom_col(width = 5 / 60) +
  scale_fill_manual(values = palette, guide = "none") +
  scale_y_continuous(breaks = c(0, .5, 1), labels = c("0%", "50%", "100%"),
                     expand = expansion(mult = c(0, .02))) +
  labs(title = "C   Smoothed probabilities", x = "Nominal minutes from segment start",
       y = "Probability") +
  flight_theme + theme(panel.grid.major = element_blank())

panels <- ((route | (altitude / uncertainty + plot_layout(heights = c(1.3, 1)))) +
  plot_layout(guides = "collect")) & theme(legend.position = "bottom")
panels + plot_annotation(
  title = "The same flight in space and altitude",
  subtitle = paste("Kuzco ·", nrow(track), "five-second steps ·", track$day[1]),
  caption = "GPS data from Keskin et al. (2024), Zenodo 12607008",
  theme = theme_ssf()
)
```

![GPS track coloured by two inferred flight states with aligned altitude
and smoothed state probabilities. The curved portions tend to coincide
with ascent.](condor-taxonomies_files/figure-html/flight-figure-1.png)

## Compare with the published flight proxies

At each step midpoint, we use the archived curvature and vertical speed
to apply the pointwise cutoffs from the flight study. The thermalling
proxy requires absolute curvature above 0.01 m⁻¹ and positive vertical
speed. The gliding proxy requires absolute curvature below 0.01 m⁻¹ and
negative vertical speed. The original study applies further bout and
wind processing. Here the proxy is a pointwise comparison, and curvature
shares information with the horizontal data in the fit. The comparison
is descriptive.

``` r

classified <- !is.na(gps$source_proxy)
agreement <- table(
  `Inferred state` = gps$state[classified],
  `Published proxy` = gps$source_proxy[classified]
)
knitr::kable(agreement)
```

|                      | Gliding proxy | Thermalling proxy |
|:---------------------|--------------:|------------------:|
| Circling-like flight |            16 |               541 |
| Directed flight      |           168 |                50 |

``` r


aligned <- agreement["Circling-like flight", "Thermalling proxy"] +
  agreement["Directed flight", "Gliding proxy"]
round(tapply(gps$vertical_m_s, gps$state, median), 2)
#> Circling-like flight      Directed flight 
#>                  0.4                 -0.4
```

Of the 775 steps with a proxy, 709 align with the corresponding flight
label. The other 529 steps have no proxy under these cutoffs. The median
measured altitude change per second is positive in the circling-like
state and negative in the directed state. These estimates describe a
single tracked bird and a single time scale.

``` r

matches <- as.data.frame(agreement)
names(matches) <- c("state", "proxy", "n")
matches$share <- ave(matches$n, matches$proxy, FUN = function(x) x / sum(x))
matches$label <- paste0(round(100 * matches$share), "% · ", matches$n)
ggplot(matches, aes(share, proxy, fill = state)) +
  geom_col(width = .55) +
  geom_text(aes(label = label), position = position_stack(vjust = .5),
            colour = "white", size = 3.2) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_x_continuous(breaks = seq(0, 1, .25), labels = paste0(seq(0, 100, 25), "%"),
                     expand = expansion(mult = 0)) +
  labs(x = "Share of steps within each proxy", y = NULL) +
  flight_theme + theme(panel.grid.major.y = element_blank())
```

![Stacked bars showing inferred flight labels within each proxy
category. Most gliding proxy steps are directed and most thermalling
proxy steps are
circling-like.](condor-taxonomies_files/figure-html/proxy-figure-1.png)

## Adapt the specification to another species

Choose the candidate labels and contrasts before inspecting the decoded
states. The constructor also accepts a ranking function of the
unlabelled movement metrics and coefficients. A rule can compare
labelled coefficients, standard errors, transition probabilities and
occupancy. If an assessment rule is false, the returned history records
its name. A custom sequence can include a simpler candidate or the
built-in `"Other"` fallback.

For a directional habitat hypothesis, include the compatible target term
in the candidate formula and define its bearing at the step origin. The
article’s four candidates remain available by leaving `taxonomies` at
its default. See
[`?ssf_taxonomy`](https://aureliennicosiaulaval.github.io/behavSSF/reference/ssf_taxonomy.md)
and
[`?fit_guided_ssf`](https://aureliennicosiaulaval.github.io/behavSSF/reference/fit_guided_ssf.md)
for the input contracts.
