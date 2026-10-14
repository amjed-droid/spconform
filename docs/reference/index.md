# Package index

## Core functions

Main spatial and spatio-temporal conformal prediction procedures

- [`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md)
  : Locally Weighted Split Conformal Prediction for Geostatistical
  (Point-Referenced) Data
- [`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md)
  : Neighbourhood-Weighted Conformal Prediction for Areal (Lattice) Data

## Spatial Diagnostics

Multi-panel diagnostic tools, spatial strata auditing, and coverage
evaluation

- [`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md)
  : Comprehensive Diagnostic Suite for Spatial Conformal Prediction
  Objects
- [`coverage_report()`](https://amjed-droid.github.io/spconform/reference/coverage_report.md)
  : Empirical Coverage and Average Interval Width for an spconform
  Object

## Helper functions

Kernel weights and graph-distance utilities

- [`spatial_kernel_weights()`](https://amjed-droid.github.io/spconform/reference/spatial_kernel_weights.md)
  : Spatial Gaussian Kernel Proximity Weights
- [`areal_neighbor_weights()`](https://amjed-droid.github.io/spconform/reference/areal_neighbor_weights.md)
  : Areal (Lattice) Graph Distance Proximity Weights

## S3 methods & generics

S3 generics and methods for spconform and diagnose objects

- [`as.data.frame(`*`<spconform>`*`)`](https://amjed-droid.github.io/spconform/reference/as.data.frame.spconform.md)
  : Coerce an spconform Object to a Data Frame
- [`plot(`*`<spconform>`*`)`](https://amjed-droid.github.io/spconform/reference/plot.spconform.md)
  : Plot Prediction Intervals for spconform Objects
- [`plot(`*`<spconform_diagnose>`*`)`](https://amjed-droid.github.io/spconform/reference/plot.spconform_diagnose.md)
  : Plot Diagnostic Audit for an spconform_diagnose Object
- [`predict(`*`<spconform>`*`)`](https://amjed-droid.github.io/spconform/reference/predict.spconform.md)
  : Extract Predictions and Intervals from an spconform Object
- [`print(`*`<spconform>`*`)`](https://amjed-droid.github.io/spconform/reference/print.spconform.md)
  : Print Method for spconform Objects
- [`residuals(`*`<spconform>`*`)`](https://amjed-droid.github.io/spconform/reference/residuals.spconform.md)
  : Extract Residuals from an spconform Object
- [`summary(`*`<spconform>`*`)`](https://amjed-droid.github.io/spconform/reference/summary.spconform.md)
  : Summary Method for spconform Objects
- [`spconform`](https://amjed-droid.github.io/spconform/reference/spconform-package.md)
  [`spconform-package`](https://amjed-droid.github.io/spconform/reference/spconform-package.md)
  : spconform: Conformal Prediction for Spatially and Spatio-Temporally
  Dependent Data
