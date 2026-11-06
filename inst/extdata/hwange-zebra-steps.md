# Hwange plains zebra steps

Prepared from Bouderbala, Nicosia and Fortin (2026), [Figshare archive,
version 2](https://doi.org/10.6084/m9.figshare.34077912.v2). The source
article is [PLOS ONE 21(10), e0360547](https://doi.org/10.1371/journal.pone.0360547).
The GPS data are also available in the archive of
[Klappstein, Thomas and Michelot (2023)](https://doi.org/10.5281/zenodo.7872602).
The prepared data are distributed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

The preparation script is `data-raw/prepare-zebra.R` in the source repository.
It retains 6,377 observed 30-minute steps and splits the series into 28
uninterrupted segments. It selects the nearest grassland cell from each step
origin using the archived vegetation raster.

`ID` identifies a segment. `time` is the end time in seconds since the first
step ended. `step` is in kilometres. `angle`, `previous_bearing` and
`target_bearing` are in radians. `target_active` marks defined target
bearings. `x0` and `y0` are step origins and `x` and `y` are endpoints,
in kilometres relative to the first step origin.
