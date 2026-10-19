# Comprehensive Diagnostic Suite for Spatial Conformal Prediction Objects

Produces a comprehensive multi-panel diagnostic report evaluating
marginal coverage validity, Winkler Interval Score (WIS) sharpness,
conditional coverage across spatial strata, boundary proximity effects,
spatial residual autocorrelation (Moran's I), and nonconformity score
distributions.

## Usage

``` r
diagnose(object, y_true, s_test = NULL, n_bins = 4, plot = TRUE, ...)
```

## Arguments

- object:

  An object of class `"spconform"`, typically output from
  [`scp_geostatistical`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md)
  or
  [`scp_areal`](https://amjed-droid.github.io/spconform/reference/scp_areal.md).

- y_true:

  Numeric vector of true observed responses at target prediction
  locations. Must have the same length as `object$pred`.

- s_test:

  Optional numeric matrix of target prediction coordinates (\\m \times
  p\\). Required for spatial stratification, boundary effects, and
  spatial autocorrelation diagnostics.

- n_bins:

  Integer; number of spatial bins per dimension for conditional coverage
  (default 4).

- plot:

  Logical; if `TRUE` (default), generates a multi-panel diagnostic plot.

- ...:

  Additional graphical parameters passed to internal plotting methods.

## Value

An S3 object of class `"spconform_diagnose"` containing:

- marginal:

  List with coverage, mean width, median width, sd width, Winkler score,
  and sample size.

- conditional:

  Data frame of coverage and widths partitioned by spatial quadrant
  (class `"spconform_strata"`).

- boundary:

  Object of class `"spconform_boundary"` comparing edge vs interior
  units.

- spatial_autocorr:

  Object of class `"spconform_moran"` with Moran's I statistics.

- scores:

  Summary of nonconformity score moments and quantiles.

- raw_scores:

  Numeric vector of nonconformity scores.

- alpha:

  Miscoverage level of the evaluated object.

## Details

The diagnostic suite coordinates modular audits across the prediction
intervals:

1.  **Marginal Validity & Sharpness**: Evaluates empirical coverage and
    the strictly proper Winkler Interval Score (WIS; Winkler 1972).

2.  **Conditional Spatial Stratification**: Calls
    [`strata_diagnostics`](https://amjed-droid.github.io/spconform/reference/strata_diagnostics.md)
    to evaluate local coverage across quadrants.

3.  **Convex Hull Boundary Proximity**: Calls
    [`boundary_diagnostics`](https://amjed-droid.github.io/spconform/reference/boundary_diagnostics.md)
    to test whether edge-effect extrapolations degrade calibration.

4.  **Spatial Autocorrelation (Moran's I)**: Calls
    [`spatial_moran`](https://amjed-droid.github.io/spconform/reference/spatial_moran.md)
    to test for residual spatial dependence.

5.  **Score Distribution**: Analyzes empirical quantiles and
    normality/QQ structure.

## References

Winkler, R. L. (1972). A Decision-Theoretic Approach to Interval
Estimation. *Journal of the American Statistical Association*, 67(337),
187–191.

Mao, H., Martin, R., and Reich, B. J. (2024). Valid Conformal Prediction
for Dependent Data. *Journal of the American Statistical Association*,
119(546), 1269–1281.
[doi:10.1080/01621459.2022.2147531](https://doi.org/10.1080/01621459.2022.2147531)
.

Susmann, H. P. (2024). Adaptive Conformal Prediction for Spatial Data.
*Computo*, [doi:10.57750/edan-5f53](https://doi.org/10.57750/edan-5f53)
.

## See also

[`spatial_moran`](https://amjed-droid.github.io/spconform/reference/spatial_moran.md),
[`boundary_diagnostics`](https://amjed-droid.github.io/spconform/reference/boundary_diagnostics.md),
[`strata_diagnostics`](https://amjed-droid.github.io/spconform/reference/strata_diagnostics.md),
[`scp_geostatistical`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md),
[`scp_areal`](https://amjed-droid.github.io/spconform/reference/scp_areal.md)

## Examples

``` r
set.seed(42)
n <- 100
s_tr <- matrix(runif(n * 2), n, 2)
y_tr <- sin(3 * s_tr[, 1]) + rnorm(n, sd = 0.2)
s_te <- matrix(runif(30 * 2), 30, 2)
y_te <- sin(3 * s_te[, 1]) + rnorm(30, sd = 0.2)

pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
out <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1)

diag_res <- diagnose(out, y_true = y_te, s_test = s_te, plot = FALSE)
print(diag_res)
#> ======================================================================
#>   spconform Diagnostic Audit Report
#> ======================================================================
#> 
#> >> 1. Marginal Coverage Validity & Sharpness:
#>    - Target Nominal Level  :  90.0%
#>    - Empirical Coverage    :  90.0% (27 of 30 covered)
#>    - Interval Width        : Mean = 1.317 | Median = 1.356 | SD = 0.056
#>    - Winkler Score (WIS)   : Mean = 1.432 (Proper Interval Score)
#> 
#> >> 2. Spatial Autocorrelation Audit (Moran's I on Residuals):
#>    - Moran's I Statistic   : +0.1376 (Expected = -0.0345)
#>    - Standardized Z-Score  : +5.1604 (p-value = 0.0000)
#>    - Residual Pattern      : Warning (Significant spatial clustering detected, p < 0.05)
#> 
#> >> 3. Conditional Coverage across Spatial Quadrants:
#>    - Strata Q1-1   : Cov =  75.0% | Mean Width =  1.357 | WIS =  1.547 (n = 4)
#>    - Strata Q1-2   : Cov =  80.0% | Mean Width =  1.357 | WIS =  1.414 (n = 5)
#>    - Strata Q1-3   : Cov = 100.0% | Mean Width =  1.396 | WIS =  1.396 (n = 2)
#>    - Strata Q2-1   : Cov = 100.0% | Mean Width =  1.306 | WIS =  1.306 (n = 2)
#>    - Strata Q2-2   : Cov = 100.0% | Mean Width =  1.256 | WIS =  1.256 (n = 2)
#>    - Strata Q2-3   : Cov = 100.0% | Mean Width =  1.357 | WIS =  1.357 (n = 3)
#>    - Strata Q3-2   : Cov = 100.0% | Mean Width =  1.256 | WIS =  1.256 (n = 2)
#>    - Strata Q3-3   : Cov = 100.0% | Mean Width =  1.356 | WIS =  1.356 (n = 1)
#>    - Strata Q3-4   : Cov = 100.0% | Mean Width =  1.306 | WIS =  1.306 (n = 4)
#>    - Strata Q4-2   : Cov = 100.0% | Mean Width =  1.249 | WIS =  1.249 (n = 1)
#>    - Strata Q4-3   : Cov =  50.0% | Mean Width =  1.243 | WIS =  2.440 (n = 2)
#>    - Strata Q4-4   : Cov = 100.0% | Mean Width =  1.243 | WIS =  1.243 (n = 2)
#> 
#> >> 4. Domain Boundary Effect (Convex Hull):
#>    - Near Boundary (Edge) : Cov =  86.7% | Mean Width =  1.319 | WIS =  1.529 (n = 15)
#>    - Far Boundary (Core) : Cov =  93.3% | Mean Width =  1.315 | WIS =  1.334 (n = 15)
#> 
#> >> 5. Nonconformity Score Distribution Moments:
#>    * Mean = 0.6585 | Median = 0.6782 | SD = 0.0280 | Q90 = 0.6784 | Q95 = 0.6891
#> ======================================================================

# Programmatic access to modular components
as.data.frame(diag_res)
#>     bin n coverage mean_width mean_winkler
#> 1  Q1-1 4     0.75   1.356630     1.547029
#> 2  Q1-2 5     0.80   1.356775     1.414175
#> 3  Q1-3 2     1.00   1.395844     1.395844
#> 4  Q2-1 2     1.00   1.305993     1.305993
#> 5  Q2-2 2     1.00   1.255501     1.255501
#> 6  Q2-3 3     1.00   1.356581     1.356581
#> 7  Q3-2 2     1.00   1.255501     1.255501
#> 8  Q3-3 1     1.00   1.356484     1.356484
#> 9  Q3-4 4     1.00   1.305993     1.305993
#> 10 Q4-2 1     1.00   1.248659     1.248659
#> 11 Q4-3 2     0.50   1.243038     2.439982
#> 12 Q4-4 2     1.00   1.243038     1.243038
```
