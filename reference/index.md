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
  : Comprehensive Diagnostic Report for spconform Objects
- [`coverage_report()`](https://amjed-droid.github.io/spconform/reference/coverage_report.md)
  : Empirical coverage and average interval width for an spconform
  object

## Helper functions

Kernel weights and graph-distance utilities

- [`spatial_kernel_weights()`](https://amjed-droid.github.io/spconform/reference/spatial_kernel_weights.md)
  : Compute Gaussian-kernel spatial weights between a target location
  and a set of reference locations
- [`areal_neighbor_weights()`](https://amjed-droid.github.io/spconform/reference/areal_neighbor_weights.md)
  : Compute neighbourhood-based weights for areal (lattice) data

## S3 methods

Print, summary, and plot methods for spconform and diagnose objects
