# Condor flight data

This CSV is an adapted subset of `Vultur_gryphus.csv` in the archive by
Keskin, Duriez, Lacerda, Flack and Nagy (2024), [Adaptive cross-country
optimisation strategies in thermal soaring birds](https://doi.org/10.5281/zenodo.12607008).
The dataset is licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The software license applies to the code.

The file contains 1,304 steps of Kuzco across seven tracking dates in 2018.
The derived fields include the nominal 5 s clock, segment IDs, horizontal
lengths and turns, altitude changes, and pointwise proxies computed from
archived curvature and vertical speed. The projected coordinates are in
metres relative to the origin used in the source file for each date. Step
lengths are in kilometres and turning angles are in radians. The extraction,
filtering and transformations are implemented in `data-raw/prepare-condor.R`
in the source project.
