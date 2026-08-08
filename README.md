# spconform

**Conformal Prediction for Spatially and Spatio-Temporally Dependent Data in R**

`spconform` provides distribution-free, finite-sample prediction intervals
for spatial and spatio-temporal data by relaxing the exchangeability
assumption of standard conformal prediction using spatial-distance /
graph-neighbourhood kernel weights.

## Why this package

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

## Installation (development version)

```r
# install.packages("remotes")
remotes::install_github("yourusername/spconform")
```

## Quick example

```r
library(spconform)

set.seed(42)
n <- 400
s <- matrix(runif(2 * n), ncol = 2)
f <- function(s) sin(4 * s[, 1]) + cos(3 * s[, 2])
y <- f(s) + rnorm(n, sd = 0.3)

idx <- sample(n, 300)
s_train <- s[idx, ]; y_train <- y[idx]
s_test  <- s[-idx, ]; y_test  <- y[-idx]

pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2], s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

out <- scp_geostatistical(s_train, y_train, s_test, pred_fun, alpha = 0.1)
coverage_report(out, y_test)
plot(out, y_true = y_test)
```

## Methodological basis

The geostatistical procedure follows the localized split-conformal idea of
Mao, Martin and Reich (2020, *"Valid model-free spatial prediction"*),
extended here to optionally include a temporal kernel for spatio-temporal
data. The areal procedure applies an analogous graph-distance weighting
scheme over a user-supplied adjacency/contiguity matrix.

## License

GPL-3
