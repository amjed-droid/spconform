# Getting Help and Troubleshooting

Welcome to the **spconform** help and troubleshooting guide. Whether you
are fitting spatial conformal prediction intervals for the first time,
diagnosing empirical coverage, or interfacing custom machine learning
models, this guide provides answers to common questions, error
diagnoses, model integration recipes, and guidelines for reporting
issues.

------------------------------------------------------------------------

## 1. Quick Support Matrix

| Need | Recommended Channel | Details |
|:---|:---|:---|
| **Bug Reports & Crashes** | [GitHub Issues](https://github.com/amjed-droid/spconform/issues) | Provide a minimal reproducible example (`reprex`) and [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). |
| **Model Advice & Theory** | [GitHub Issues (Question)](https://github.com/amjed-droid/spconform/issues) | Questions on spatial exchangeability, choice of kernel, or sample size. |
| **Feature Requests** | [GitHub Issues (Enhancement)](https://github.com/amjed-droid/spconform/issues/new) | Propose new diagnostic metrics, weight kernels, or interface integrations. |
| **Academic Inquiries** | Direct Email | Reach maintainer Ahmed Sattar Jabbar at `ahmed.state.me@gmail.com`. |

------------------------------------------------------------------------

## 2. Frequently Asked Questions (FAQ)

### Q1: My empirical coverage is slightly below $`1 - \alpha`$ (e.g. 88.5% instead of 90%). Is my model invalid?

**Not necessarily.** Conformal prediction guarantees valid marginal
coverage *in expectation*:
``` math
P\left(Y_{n+1} \in \widehat{C}(S_{n+1})\right) \ge 1 - \alpha
```
For any finite test evaluation set of size $`m`$, the observed empirical
coverage follows a Binomial distribution $`\text{Binom}(m, 1 - \alpha)`$
with standard error:
``` math
\text{SE} = \sqrt{\frac{\alpha(1 - \alpha)}{m}}
```
For example, with $`m = 100`$ and nominal coverage $`90\%`$
($`\alpha = 0.10`$), the standard error is
$`\sqrt{0.1 \times 0.9 / 100} = 0.03`$ ($`3\%`$). An observed coverage
between $`84\%`$ and $`96\%`$ falls well within normal $`95\%`$ binomial
sampling variation. You can compute the exact binomial confidence
interval via:

``` r

report <- coverage_report(out, y_test)
binom.test(round(report$coverage * length(y_test)), length(y_test), p = 1 - out$alpha)
```

### Q2: Should I use `scp_geostatistical()` or `scp_areal()`?

- Use
  [`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md)
  when observations have continuous coordinates $`(x, y)`$ or
  $`(\text{lon}, \text{lat})`$, such as weather stations, air quality
  sensors, or soil core samples.
- Use
  [`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md)
  when data are aggregated over discrete spatial units (counties,
  municipalities, census tracts, or lattice grid cells) where
  relationships are defined by an adjacency matrix $`W`$.

### Q3: How do spatial bandwidth and $`k`$-NN affect interval widths?

- **Bandwidth ($`h`$):** In
  [`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md),
  weights decay as $`w_i = \exp(-\|s_0 - s_i\|^2 / (2h^2))`$.
  - A **smaller $`h`$** restricts calibration to immediately proximate
    neighbors, producing localized intervals that adapt to local
    heteroscedasticity, but may be noisy if local data are sparse.
  - A **larger $`h`$** borrows strength from more calibration points,
    yielding smoother, more stable intervals that approach standard
    split conformal prediction.
  - The default (`bandwidth = NULL`) employs the median pairwise
    distance heuristic, which provides a robust balance for isotropic
    spatial domains.
- **$`k`$-Nearest Neighbors (`k_neighbors`):** In `spconform` 0.1.2+,
  for large calibration datasets ($`N > 1000`$), setting
  e.g. `k_neighbors = 50` restricts kernel weighting to the 50 closest
  calibration stations, accelerating computation up to 50x while
  preserving asymptotic validity.

### Q4: Does `spconform` support spatio-temporal data?

**Yes.** Supply `t_train` and `t0` alongside coordinates in
[`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md):

``` r

out <- scp_geostatistical(
  s_train = s_train, y_train = y_train, s0 = s_test,
  pred_fun = my_st_predictor,
  t_train = t_train, t0 = t_test,
  temporal_bandwidth = 7,  # e.g. 7 days
  alpha = 0.1
)
```

------------------------------------------------------------------------

## 3. Interpreting the Diagnostic Suite (`diagnose()`)

The
[`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md)
function conducts five audits on spatial prediction intervals:

``` r

diag <- diagnose(out, y_true = y_test, s_test = s_test, plot = TRUE)
print(diag)
```

### 1. Marginal Coverage & Interval Width

- **`coverage`**: Fraction of true values falling within
  $`[L(s_0), U(s_0)]`$. Should be $`\ge 1 - \alpha`$.
- **`mean_width` vs `median_width`**: If the mean width is substantially
  larger than the median, spatial outliers or high-uncertainty boundary
  clusters are inflating interval widths.

### 2. Winkler Interval Score (WIS)

A strictly proper scoring rule balancing interval sharpness and
miscoverage penalties:
``` math
\text{WIS}_\alpha(l, u, y) = (u - l) + \frac{2}{\alpha}(l - y)\mathbf{1}(y < l) + \frac{2}{\alpha}(y - u)\mathbf{1}(y > u)
```
Lower scores indicate superior predictive performance (narrower
intervals with minimal violations).

### 3. Spatial Quadrant Stratification

Partitions the 2D spatial bounding box into $`n \times n`$ quadrants
(default $`4 \times 4 = 16`$). \* Inspect whether any specific quadrant
suffers severe localized under-coverage (e.g. coverage $`< 75\%`$). \*
Local deficits often highlight unobserved spatial covariates or sharp
topographical transitions.

### 4. Convex Hull Boundary Proximity

Examines whether intervals near the edge of the spatial domain
experience degraded coverage due to spatial extrapolation.
[`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md)
partitions test points into **Interior** vs. **Boundary** units.

### 5. Spatial Autocorrelation (Moran’s $`I`$)

Computes Moran’s $`I`$ on absolute conformal residuals: \*
**$`I \approx 0`$ ($`p > 0.05`$):** Nonconformity scores are spatially
independent; local exchangeability holds well. \* **$`I > 0`$
($`p < 0.05`$):** Prediction errors cluster spatially, indicating
unmodeled spatial autocorrelation in the underlying point predictor.

------------------------------------------------------------------------

## 4. Model Integration Recipes (`pred_fun`)

`spconform` is strictly **model-agnostic**. The `pred_fun` must have the
signature: `function(s_train, y_train, s_new)` and return a plain 1D
numeric vector of length `nrow(s_new)`.

### Recipe 1: Random Forest (`ranger`)

``` r

library(ranger)

pred_rf <- function(s_train, y_train, s_new) {
  df_tr <- data.frame(x = s_train[, 1], y = s_train[, 2], target = y_train)
  df_te <- data.frame(x = s_new[, 1],   y = s_new[, 2])
  
  fit <- ranger(target ~ x + y, data = df_tr, num.trees = 100)
  predict(fit, data = df_te)$predictions
}
```

### Recipe 2: Generalized Additive Model (`mgcv::gam`)

``` r

library(mgcv)

pred_gam <- function(s_train, y_train, s_new) {
  df_tr <- data.frame(x = s_train[, 1], y = s_train[, 2], target = y_train)
  df_te <- data.frame(x = s_new[, 1],   y = s_new[, 2])
  
  fit <- gam(target ~ s(x, y, k = 20), data = df_tr)
  as.numeric(predict(fit, newdata = df_te))
}
```

### Recipe 3: Ordinary Spatial Kriging (`gstat`)

``` r

library(sp)
library(gstat)

pred_kriging <- function(s_train, y_train, s_new) {
  sp_tr <- SpatialPointsDataFrame(coords = s_train, data = data.frame(y = y_train))
  sp_te <- SpatialPoints(coords = s_new)
  
  vgm_mod <- fit.variogram(variogram(y ~ 1, sp_tr), vgm("Exp"))
  krig_out <- krige(y ~ 1, sp_tr, sp_te, model = vgm_mod, debug.level = 0)
  krig_out$var1.pred
}
```

------------------------------------------------------------------------

## 5. Troubleshooting Common Errors & Warnings

This section provides an exhaustive guide to diagnosing and resolving
statistical, computational, and dimensional errors encountered across
[`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md),
[`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md),
and
[`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md).

------------------------------------------------------------------------

### 1. Coordinate & Dimension Mismatches

#### Error: `s_train and s0 must have the same number of coordinate columns`

- **Cause:** The number of spatial coordinate dimensions differs between
  training and query locations (for example, passing 2D coordinates
  `(x, y)` in `s_train` but 3D `(x, y, t)` in `s0`).

- **Fix:** Verify column counts prior to execution:

  ``` r

  stopifnot(ncol(s_train) == ncol(s0))
  ```

#### Issue: Single-point query vector treated as 1D vector instead of 1-row matrix

- **Cause:** Passing a simple coordinate vector `s0 = c(50.2, 14.3)`
  rather than a 1-row matrix or data frame.

- **Fix:** Coerce single query coordinates explicitly:

  ``` r

  s0 <- matrix(c(50.2, 14.3), nrow = 1)
  ```

------------------------------------------------------------------------

### 2. Missing & Non-Finite Values (`NA`, `NaN`, `Inf`)

#### Error: `s_train cannot contain NA or non-finite values` or `y_train contains missing values`

- **Cause:** Sensor dropouts, unrecorded responses, or missing spatial
  coordinates in training or target sets.

- **Fix:** Synchronize and filter rows simultaneously before passing
  data to `spconform`:

  ``` r

  valid_tr <- complete.cases(s_train, y_train)
  s_train  <- s_train[valid_tr, , drop = FALSE]
  y_train  <- y_train[valid_tr]

  valid_te <- complete.cases(s0)
  s0       <- s0[valid_te, , drop = FALSE]
  ```

------------------------------------------------------------------------

### 3. Custom Predictor (`pred_fun`) Return Types & Single-Row Drops

#### Error: `pred_fun must return a numeric vector of length nrow(s_new)`

- **Cause 1 (Type Mismatch):** Your model’s
  [`predict()`](https://rdrr.io/r/stats/predict.html) method returns a
  `matrix` with 1 column, a `data.frame`, or an object with custom
  attributes rather than an atomic numeric vector.

- **Fix 1:** Always wrap the return value in
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html):

  ``` r

  pred_fun <- function(s_train, y_train, s_new) {
    fit <- ranger(y ~ ., data = data.frame(s_train, y = y_train))
    as.numeric(predict(fit, data = data.frame(s_new))$predictions)
  }
  ```

- **Cause 2 (Single-Row Dimension Drop):** Subsetting `s_new[i, ]`
  inside `pred_fun` drops a 1-row matrix to a vector, causing downstream
  matrix multiplications or feature extractors to fail.

- **Fix 2:** Always include `drop = FALSE` in internal subsetting
  operations:

  ``` r

  s_sub <- s_new[idx, , drop = FALSE]
  ```

------------------------------------------------------------------------

### 4. Coordinate Reference System (CRS) & Distance Scale Discrepancies

#### Symptom: Prediction intervals are excessively wide or bandwidth estimation is degenerate

- **Cause:** Mixing unprojected geographic coordinates (degrees of
  longitude and latitude) with planar coordinates (UTM meters), or
  passing coordinates where axes have wildly disparate scales. In
  unprojected coordinates, Euclidean distance fails to account for earth
  curvature.

- **Fix:** Reproject coordinates to a conformal planar projection
  (e.g. UTM or EPSG:3857) prior to calibration:

  ``` r

  # Example using sf:
  # pts_projected <- sf::st_transform(pts_sf, crs = 32638)
  # s_train <- sf::st_coordinates(pts_projected)
  ```

------------------------------------------------------------------------

### 5. Areal Adjacency Matrix Structural Errors

#### Error: `adjacency must be a square numeric matrix`

- **Cause:** Passing a non-square matrix
  (`nrow(adjacency) != ncol(adjacency)`) or an adjacency list rather
  than a binary/weighted contiguity matrix.

- **Fix:** Validate matrix dimensions:

  ``` r

  stopifnot(is.matrix(adjacency), nrow(adjacency) == ncol(adjacency))
  ```

#### Warning: `Effective neighbourhood weight vanishingly small; falling back to unweighted quantile`

- **Cause:** An areal unit has no graph neighbors in the training set
  (e.g. an isolated island, enclaved district, or disconnected graph
  component with degree zero).

- **Fix:** Audit graph connectivity. Identify isolated units and bridge
  graph components or adjust contiguity criteria:

  ``` r

  isolated <- which(rowSums(adjacency) == 0)
  if (length(isolated) > 0) {
    message("Isolated units detected: ", paste(isolated, collapse = ", "))
  }
  ```

------------------------------------------------------------------------

### 6. Calibration Sample Size vs. Nominal Miscoverage ($`\alpha`$)

#### Symptom: Conformal quantile returns the maximum nonconformity score or upper bound is uninformative

- **Cause (Mathematical Degeneracy):** The held-out calibration set is
  too small for the requested miscoverage level $`\alpha`$. Conformal
  prediction requires:
  ``` math
  N_{\text{cal}} \ge \left\lceil \frac{1 - \alpha}{\alpha} \right\rceil
  ```
  For example, requesting a 99% interval ($`\alpha = 0.01`$) when
  $`N_{\text{cal}} = 40`$ mathematically cannot calibrate a 99%
  empirical quantile, as the highest empirical quantile available is
  $`1 - 1/(N_{\text{cal}} + 1) = 97.6\%`$.
- **Fix:** Either increase the calibration partition size via `split`
  (e.g. `split = 0.5` instead of `0.8`), or choose an $`\alpha`$
  commensurate with sample size (e.g. $`\alpha = 0.05`$ or $`0.10`$).

------------------------------------------------------------------------

### 7. Evaluation & Diagnostic Alignment

#### Error: `y_true must have the same length as object$pred`

- **Cause:** Passing the entire response vector `y` instead of the
  held-out evaluation subset corresponding to `s0`.

- **Fix:** Align evaluation indices:

  ``` r

  # Correct evaluation subset:
  coverage_report(out, y_true = y_test)
  diagnose(out, y_true = y_test, s_test = s_test)
  ```

#### Error: `s_test is required for spatial stratification and boundary audits`

- **Cause:** Calling `diagnose(out, y_true = y_test)` without providing
  query coordinates `s_test`. Without coordinates,
  [`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md)
  cannot compute spatial quadrant partitions, convex hull boundaries, or
  Moran’s $`I`$.

- **Fix:** Always provide `s_test`:

  ``` r

  diag <- diagnose(out, y_true = y_test, s_test = s_test, plot = TRUE)
  ```

------------------------------------------------------------------------

### 8. Graphics Layout & Multi-Panel Plot Devices

#### Warning: Plot margins too large or layout display clipped in RStudio

- **Cause:** The multi-panel diagnostic plot
  ([`diagnose()`](https://amjed-droid.github.io/spconform/reference/diagnose.md))
  renders 4 to 6 visual audits simultaneously, requiring adequate
  display geometry.

- **Fix:** Open an expanded graphic window before plotting, or extract
  numerical diagnostics without plotting:

  ``` r

  # Numerical diagnostics only (no graphical device):
  diag <- diagnose(out, y_true = y_test, s_test = s_test, plot = FALSE)
  print(diag)

  # For plotting, expand display device:
  dev.new(width = 11, height = 7)
  plot(diag)
  ```

------------------------------------------------------------------------

## 6. Submitting a Minimal Reproducible Example (`reprex`)

If you encounter unexpected behavior, please create a self-contained
reproducible example:

``` r

library(spconform)

# 1. Synthetic reproducible dataset
set.seed(42)
n_tr <- 80
n_te <- 20
s_train <- matrix(runif(n_tr * 2), ncol = 2)
y_train <- sin(3 * s_train[, 1]) + rnorm(n_tr, sd = 0.2)

s_test  <- matrix(runif(n_te * 2), ncol = 2)
y_test  <- sin(3 * s_test[, 1])  + rnorm(n_te, sd = 0.2)

# 2. Simple baseline predictor
pfun <- function(s_tr, y_tr, s_new) {
  rep(mean(y_tr), nrow(s_new))
}

# 3. Execution triggering the issue
out <- scp_geostatistical(s_train, y_train, s_test, pfun, alpha = 0.1)
cov <- coverage_report(out, y_test)
print(cov)

# 4. Environment details
sessionInfo()
```

------------------------------------------------------------------------

## 7. Citation & Academic Reference

If you use **spconform** in academic research, please cite:

``` bibtex
@Manual{spconform2026,
  title  = {spconform: Conformal Prediction for Spatially and Spatio-Temporally Dependent Data in R},
  author = {Ahmed Sattar Jabbar},
  year   = {2026},
  note   = {R package version 0.1.2},
  url    = {https://CRAN.R-project.org/package=spconform},
  doi    = {10.32614/CRAN.package.spconform}
}
```

------------------------------------------------------------------------

## 8. Maintainer Contact

- **Lead Author & Maintainer:** Ahmed Sattar Jabbar
- **Email:** <ahmed.state.me@gmail.com>
- **ORCID:**
  [0009-0006-4614-735X](https://orcid.org/0009-0006-4614-735X)
- **GitHub Repository:**
  [amjed-droid/spconform](https://github.com/amjed-droid/spconform)
