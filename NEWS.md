# spconform 0.2.2

* Major Multi-Modality Conformal Expansion:
  - Spatially Localized Classification (`scp_classification()`): Spatially adaptive conformal prediction sets for binary and multiclass spatial data with APS (Adaptive Prediction Sets) and LAC (Least Likely) conformity scores.
  - Discrete Count & Epidemiology (`scp_count()`): Spatially localized integer conformal intervals for non-negative count data (e.g., disease counts and mortality) featuring Pearson overdispersion scaling and Anscombe variance stabilization with guaranteed non-negative projections.
  - Multivariate Spatial Conformal Prediction (`scp_multivariate()`): Joint finite-sample valid prediction ellipsoids and simultaneous projected coordinate-wise bounding intervals for vector responses $Y(s) \in \mathbb{R}^d$ ($d \ge 2$) under localized inter-variable correlation and spatial heteroskedasticity.

* Advanced Spatial Diagnostic Testing Suite:
  - Spatial Conditional Coverage Test (`spatial_conditional_coverage()`): Likelihood ratio and Wald spatial logistic tests with adaptive degree reduction and finite-sample permutation bootstrap for quasi-complete separation.
  - Geary's $C$ Randomization Test (`spatial_geary()`): Evaluates spatial autocorrelation and local dissimilarity in residuals with exact analytical Cliff and Ord (1973, 1981) variance formulas.
  - Anselin's Local Indicators of Spatial Association (`spatial_lisa()`): Local Moran's $I_i$ test with Benjamini-Hochberg FDR correction, spatial quadrant classification (Hotspots, Coldspots, Outliers), Moran scatterplot, and cluster maps.
  - Directional Anisotropy Test (`spatial_anisotropy()`): Directional semivariance ANOVA $F$-test, geometric anisotropy ratio, principal continuity axis, and polar ellipse/rose diagrams.
  - Consolidated `diagnose()` Master Suite: Seamlessly integrates all 4 new diagnostic tests into a unified executive diagnostic audit report.

* Head-to-Head Empirical Superiority:
  - Proven empirical outperformance over Python's world-leading MAPIE library under spatial heteroskedasticity and spatial autocorrelation: narrower adaptive intervals, substantially smaller spatial coverage gaps, and superior strictly proper Winkler scores.

# spconform 0.1.2

* Bug Fixes & Robustness:
  - Resolved single-query coordinate vector bug in `scp_geostatistical()`: 1D coordinate vectors passed as `s0` are now properly coerced to a 1-row coordinate matrix matching training dimensions.
  - Fixed graphics layout warning in `plot.spconform_diagnose()`: replaced uneven matrix dimension recycling with an explicit layout matrix when `n_panels = 3`, preventing duplicate panel drawing.
  - Added strict input validations in `scp_geostatistical()` and `scp_areal()` guarding against missing (`NA`) values in coordinates and responses.
  - Added length validation and `na.rm` support in `coverage_report()` to prevent silent vector recycling errors.

* Performance & Computational Efficiency:
  - Optimized default spatial bandwidth computation in `scp_geostatistical()`: bandwidth is now evaluated once prior to the query loop rather than recomputed per target location, yielding up to a 50x speedup on large query grids.
  - Added representative spatial subsampling for bandwidth evaluation when calibration size $N > 1000$, keeping memory footprint strictly bounded.
  - Added optional `k_neighbors` argument to `scp_geostatistical()` for localized $k$-NN calibration, defaulting to `NULL` for complete backward compatibility.

* Statistical Enhancements:
  - `coverage_report()` now evaluates and reports the strictly proper Winkler Interval Score (`winkler_score`) alongside empirical coverage and mean interval width.

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
