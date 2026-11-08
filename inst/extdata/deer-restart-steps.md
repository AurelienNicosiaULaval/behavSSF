# Red deer steps

Prepared from `deer` and `sh_forest` in `amt` 0.3.1.0 (GPL-3).
The GPS fixes are credited to Verein für Wildtierforschung Dresden und
Göttingen e.V.; forest cover is credited to the European Commission's
Joint Research Centre. See the original `amt` data documentation for the
source acknowledgements: <https://cran.r-project.org/package=amt>.

The preparation script is `data-raw/prepare-deer.R`. It retains steps whose
interval and preceding interval are within ten minutes of six hours and
belong to the same original burst. The preceding interval supplies the
incoming heading. This gives 759 steps in 29 uninterrupted segments.
The matched comparison with `hmmSSF` omits three singleton segments,
leaving 756 steps in 26 segments.

Lengths and relative coordinates are in kilometres, bearings and turns in
radians, and `time_utc` records the original endpoint times. `end_in_forest`
records the forest raster value at each observed endpoint. It describes
habitat use, without an availability comparison. Source fix and burst
identifiers are retained for tracing the preprocessing.
