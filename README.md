# spconform

[![R-CMD-check](https://github.com/amjed-droid/spconform/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/amjed-droid/spconform/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/spconform)](https://CRAN.R-project.org/package=spconform)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Conformal Prediction for Spatially and Spatio-Temporally Dependent Data in R**

`spconform` provides distribution-free, finite-sample prediction intervals
for spatial and spatio-temporal data by relaxing the exchangeability
assumption of standard conformal prediction using spatial-distance /
graph-neighbourhood kernel weights.

---

## Why `spconform`?

Existing R conformal-prediction packages target either:

- purely **temporal** dependence (`conformalForecast`, `AdaptiveConformal`), or
- i.i.d./exchangeable data (`conformalInference`, `conformalClassification`, `cfcausal`).

No published, documented R package on CRAN currently provides a unified
framework for **both geostatistical (point-referenced)** and **areal
(lattice)** spatial data, nor for the **spatio-temporal** case. A related
Python package (`geoconformal`) exists but is not part of the R ecosystem
and is not published in a statistical-software journal. A minimal,
undocumented GitHub prototype (`scp`) implements a single geostatistical
method without areal-data support, testing, or a CRAN release.

`spconform` fills this gap with:

1. `scp_geostatistical()` — locally weighted split conformal prediction for
   point-referenced spatial (or spatio-temporal) data, compatible with any
   user-supplied point predictor (kriging, GAM, random forest, ...).
2. `scp_areal()` — leave-one-out, neighbourhood-weighted conformal
   prediction for lattice/areal data using an adjacency structure.
3. `coverage_report()`, `print()`, `summary()`, `plot()` methods for
   diagnostics.

---

## Installation

```r
# Install from GitHub (development version)
remotes::install_github("amjed-droid/spconform")

# Alternatively, once on CRAN:
# install.packages("spconform")
