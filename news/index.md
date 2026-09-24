# Changelog

## spconform 0.1.1

CRAN release: 2026-09-24

- Enhanced S3 generics and methods for `spconform` objects:
  - Added
    [`predict.spconform()`](https://amjed-droid.github.io/spconform/reference/predict.spconform.md)
    to extract point predictions or prediction intervals (`fit`, `lwr`,
    `upr`).
  - Added
    [`residuals.spconform()`](https://amjed-droid.github.io/spconform/reference/residuals.spconform.md)
    to compute raw response residuals or absolute calibration scores.
  - Added
    [`as.data.frame.spconform()`](https://amjed-droid.github.io/spconform/reference/as.data.frame.spconform.md)
    to coerce conformal objects into tidy data frames.
- Enhanced spatial diagnostics:
  - [`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md)
    now returns a classed `spconform_diagnose` object with dedicated
    `print.spconform_diagnose()` and
    [`plot.spconform_diagnose()`](https://amjed-droid.github.io/spconform/reference/plot.spconform_diagnose.md)
    methods.
  - Added Moran’s $`I`$ test on prediction residuals and conformal
    hit/miss indicators.
- Documentation and code compliance:
  - Replaced all non-English code comments with English comments in
    `R/scp_areal.R` for full ASCII compliance.
  - Standardized all Rd help page titles to Title Case style.
  - Updated citation metadata to Mao, Martin, and Reich (JASA 2024).

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
