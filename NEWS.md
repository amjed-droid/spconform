

# spconform 0.1.1

* Enhanced S3 generics and methods for `spconform` objects:
  - Added `predict.spconform()` to extract point predictions or prediction intervals (`fit`, `lwr`, `upr`).
  - Added `residuals.spconform()` to compute raw response residuals or absolute calibration scores.
  - Added `as.data.frame.spconform()` to coerce conformal objects into tidy data frames.
* Enhanced spatial diagnostics:
  - `diagnose()` now returns a classed `spconform_diagnose` object with dedicated `print.spconform_diagnose()` and `plot.spconform_diagnose()` methods.
  - Added Moran's $I$ test on prediction residuals and conformal hit/miss indicators.
* Documentation and code compliance:
  - Replaced all non-English code comments with English comments in `R/scp_areal.R` for full ASCII compliance.
  - Standardized all Rd help page titles to Title Case style.
  - Updated citation metadata to Mao, Martin, and Reich (JASA 2024).

# spconform 0.1.0

* Initial release of `spconform` on CRAN.
* Implemented locally weighted split conformal prediction for geostatistical (point-referenced) data (`scp_geostatistical()`).
* Implemented neighbourhood-weighted leave-one-out conformal prediction for areal (lattice) data (`scp_areal()`).
* Added kernel weighting utilities for spatial and spatio-temporal predictions (`spatial_kernel_weights()`, `areal_neighbor_weights()`).
* Implemented comprehensive spatial diagnostic tools and multi-panel visualization (`diagnose()`).
* Included S3 methods for printing, summarizing, and plotting conformal prediction intervals (`print()`, `summary()`, `plot()`, `coverage_report()`).
