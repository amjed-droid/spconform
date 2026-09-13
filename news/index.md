# Changelog

## spconform 0.1.0

CRAN release: 2026-09-12

- Initial release of `spconform` on CRAN.
- Implemented locally weighted split conformal prediction for
  geostatistical (point-referenced) data
  ([`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md)).
- Implemented neighbourhood-weighted leave-one-out conformal prediction
  for areal (lattice) data
  ([`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md)).
- Added kernel weighting utilities for spatial and spatio-temporal
  predictions
  ([`spatial_kernel_weights()`](https://amjed-droid.github.io/spconform/reference/spatial_kernel_weights.md),
  [`areal_neighbor_weights()`](https://amjed-droid.github.io/spconform/reference/areal_neighbor_weights.md)).
- Implemented comprehensive spatial diagnostic tools and multi-panel
  visualization
  ([`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md)).
- Included S3 methods for printing, summarizing, and plotting conformal
  prediction intervals ([`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`coverage_report()`](https://amjed-droid.github.io/spconform/reference/coverage_report.md)).
