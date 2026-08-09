# spconform

[![R-CMD-check](https://github.com/amjed-droid/spconform/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/amjed-droid/spconform/actions/workflows/R-CMD-check.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21862025.svg)](https://doi.org/10.5281/zenodo.21862025)

**Conformal Prediction for Spatially and Spatio-Temporally Dependent Data in R**

`spconform` provides distribution-free prediction intervals with finite-sample coverage properties for spatial and spatio-temporal data. It relaxes the exchangeability assumption of standard conformal prediction by using spatial-distance and graph-neighbourhood kernel weights, offering a unified framework for both **geostatistical (point-referenced)** and **areal (lattice)** data structures.

---

## Why spconform?

Standard conformal prediction assumes exchangeable data — an assumption routinely violated in spatial settings where nearby observations are more similar than distant ones. Existing R packages address either purely **temporal** dependence (`conformalForecast`, `AdaptiveConformal`) or **i.i.d./exchangeable** data (`conformalInference`, `conformalClassification`, `cfcausal`), but none provides a documented, unit-tested, unified solution for spatial data on CRAN.

| Feature | `conformalInference` | `conformalForecast` | `scp` (GitHub) | `geoconformal` (Python) | **`spconform`** |
|:---|:---:|:---:|:---:|:---:|:---:|
| Language | R | R | R | Python | **R** |
| Geostatistical (point-ref.) | — | — | ✓ | ✓ | **✓** |
| Areal (lattice) | — | — | — | — | **✓** |
| Spatio-temporal | — | — (temp. only) | — | — | **✓ (opt.)** |
| Model-agnostic | ✓ | ✓ | ✓ | ✓ | **✓** |
| Unit-tested / CRAN-ready | ✓ | ✓ | — | — | **✓** |

`spconform` is, to our knowledge, the first R package to offer conformal prediction spanning both major spatial data structures with optional spatio-temporal extension.

> **Status:** `spconform` has passed `R CMD check --as-cran` with 0 errors,
> 0 warnings, 0 notes (tested locally and via win-builder on R-release and
> R-devel) and has been submitted to CRAN. A permanent, citable snapshot of
> version 0.1.0 is archived on Zenodo (DOI above). The accompanying
> manuscript is currently in preparation for submission to the *Journal of
> Statistical Software*.

---

## Installation

```r
# Install the development version from GitHub
remotes::install_github("amjed-droid/spconform")

# Once accepted on CRAN:
# install.packages("spconform")
```

**Dependencies:** The package imports only `stats` (base R). Suggested
packages (`sp`, `knitr`, `rmarkdown`) are used to build and run the
vignette; `mgcv`, `ranger`, and `bmstdr` are only needed to reproduce the
extended examples shown in the accompanying paper and are not required
for core package functionality.

---

## Quick start

### 1. Geostatistical data (point-referenced)

```r
library(spconform)
library(sp)

data(meuse)
s <- as.matrix(meuse[, c("x", "y")])
y <- log(meuse$zinc)

# Any user-supplied point predictor
pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
            I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2],
        s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

# 90% locally weighted conformal intervals
set.seed(123)
idx <- sample(nrow(s), floor(0.7 * nrow(s)))
out <- scp_geostatistical(s[idx, ], y[idx], s[-idx, ], pred_fun, alpha = 0.1)

print(out)
#> <spconform> geostatistical conformal prediction
#> Target coverage: 90.0%
#> Number of prediction points: 47

coverage_report(out, y[-idx])
#> $coverage
#> [1] 0.957
#> $mean_width
#> [1] 2.21
```

### 2. Areal / lattice data

```r
# Aggregate Meuse to a 6x6 grid (21 occupied cells)
xbreaks <- seq(min(meuse$x), max(meuse$x), length.out = 7)
ybreaks <- seq(min(meuse$y), max(meuse$y), length.out = 7)
meuse$cell_x  <- cut(meuse$x, xbreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_y  <- cut(meuse$y, ybreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_id <- (meuse$cell_y - 1) * 6 + meuse$cell_x

agg <- aggregate(log(zinc) ~ cell_id, data = meuse, FUN = mean)
names(agg) <- c("cell_id", "y")
cell_coords <- unique(meuse[, c("cell_id", "cell_x", "cell_y")])
agg <- merge(agg, cell_coords, by = "cell_id")
agg <- agg[order(agg$cell_id), ]

# Build adjacency matrix (Queen contiguity)
n_cells <- nrow(agg)
adj <- matrix(0, n_cells, n_cells)
for (i in 1:n_cells) {
  for (j in 1:n_cells) {
    if (i != j) {
      dx <- abs(agg$cell_x[i] - agg$cell_x[j])
      dy <- abs(agg$cell_y[i] - agg$cell_y[j])
      if (dx <= 1 && dy <= 1) adj[i, j] <- 1
    }
  }
}

# 80% neighbourhood-weighted conformal intervals
out_areal <- scp_areal(agg$y, adjacency = adj, alpha = 0.2, decay = 0.5)
coverage_report(out_areal, agg$y)
#> $coverage
#> [1] 0.81
#> $mean_width
#> [1] 1.79
```

### 3. Spatio-temporal data

```r
library(mgcv)
library(bmstdr)

data("nysptime")
df <- nysptime[complete.cases(nysptime[, c("utmx", "utmy", "y8hrmax", "Day", "Month")]), ]
df$day_idx <- ifelse(df$Month == 7, df$Day, 31 + df$Day)

s <- as.matrix(df[, c("utmx", "utmy")])
t <- df$day_idx
s_3d <- cbind(s, t)
y <- df$y8hrmax

# Spatio-temporal GAM predictor
pred_fun_st <- function(s_train, y_train, s_new) {
  train_df <- data.frame(x = s_train[, 1], y = s_train[, 2],
                         day = s_train[, 3], z = y_train)
  fit <- gam(z ~ te(x, y, day, k = c(8, 8, 4)), data = train_df)
  new_df <- data.frame(x = s_new[, 1], y = s_new[, 2], day = s_new[, 3])
  as.numeric(predict(fit, newdata = new_df))
}

set.seed(123)
n <- nrow(s_3d)
train_idx <- sample(n, floor(0.7 * n))
out_st <- scp_geostatistical(
  s_train = s_3d[train_idx, ],
  y_train = y[train_idx],
  s0 = s_3d[-train_idx, ],
  pred_fun = pred_fun_st,
  t_train = t[train_idx],
  t0 = t[-train_idx],
  temporal_bandwidth = 5,
  alpha = 0.1,
  split = 0.5
)

coverage_report(out_st, y[-train_idx])
#> $coverage
#> [1] 0.909
#> $mean_width
#> [1] 41.7
```

---

## Key features

- **Model-agnostic**: Works with any user-supplied point predictor (kriging, GAM, random forest, linear model, ...).
- **Finite-sample coverage**: Maintains coverage close to nominal level regardless of predictor misspecification (under local exchangeability).
- **Lightweight**: Imports only `stats`; no heavy spatial-modelling dependencies.
- **Fully documented**: S3 methods (`print`, `summary`, `plot`, `coverage_report`) included, plus a full introductory vignette.

---

## Citation

If you use `spconform` in your research, please cite:

> Jabbar, A. S. (2026). *spconform: Conformal Prediction for Spatially and
> Spatio-Temporally Dependent Data in R* (Version 0.1.0) [Computer software].
> Zenodo. https://doi.org/10.5281/zenodo.21862025
>
> *A companion manuscript describing the package methodology is currently
> submitted to the Journal of Statistical Software and will be cited here
> upon acceptance.*

```r
citation("spconform")
```

---

## Getting help

- **Bug reports & feature requests**: [GitHub Issues](https://github.com/amjed-droid/spconform/issues)
- **Documentation**: `?scp_geostatistical`, `?scp_areal`, `vignette("spconform-intro", package = "spconform")`
- **Reproducible scripts**: See `inst/scripts/` in the package source.

---

## License

This package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
