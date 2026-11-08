# Zebra state comparison

The vegetation grid and track come from Bouderbala, Nicosia and Fortin (2026),
[Figshare archive v2](https://doi.org/10.6084/m9.figshare.34077912.v2),
[article](https://doi.org/10.1371/journal.pone.0360547). The input data are
redistributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Original Codes.zip MD5 is `9c5681ed705776999891bb53b62424db`.

`hwange-zebra-habitat.rds` stores the original 1,510 by 1,686 vegetation
cells as row-major integers, their extent in coordinates relative to the
first step origin, and their resolution in kilometres. Codes 1 to 4 denote
grassland, bushed grassland, bushland and woodland. Preparation is in
`data-raw/prepare-zebra-habitat.R` in the repository.

`zebra-state-comparison.R` defines example helpers. `zebra_state_data()`
reconstructs the matched choice sets, and `zebra_state_initializers()` makes
training-only starts. `zebra_state_initializers_extended()` adds perturbations
for four states: eight starts on the full track and four per training split.
They are example code, not an exported modelling API.
The targets were defined at the step origins in the prepared GPS data.
Two singleton segments are excluded for the established matched cohort,
leaving 6,375 steps in 26 segments. The chronological blocks contain 1,961,
2,085 and 2,329 steps. Seed 660301 generates 128 proposed controls per step,
with 565 outside-map controls discarded. The Cartesian sampling density
is the raw gamma/log-normal length mixture divided by `2*pi*step`.
After map rejection this is proportional to the conditioned proposal.
The common factor cancels from conditional choices; an absolute physical
normaliser would require the conditioned density or the original draw count.

`zebra-state-comparison.rds` is the saved result of the full code shown in
the state-count vignette, run on 7 November 2026 with behavSSF 0.1.0.9009 and
R 4.5.0 on macOS arm64 using six workers. It includes the starts audit, movement profiles and coefficients for each successful start,
block scores, full-data fits and assessments for one to four states under
gamma and Weibull lengths. The cache avoids repeating 32 fitting jobs when
rendering the vignette. The vignette gives code for an independent refit.
All candidates use the same controls and the origin-based grassland target.
