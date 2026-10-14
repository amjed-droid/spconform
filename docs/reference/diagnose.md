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

  Data frame of coverage and widths partitioned by spatial quadrant.

- boundary:

  List comparing empirical metrics between interior and boundary units.

- spatial_autocorr:

  List containing Moran's I statistic, expected value, z-score, and
  p-value.

- scores:

  Summary of nonconformity score moments and quantiles.

- alpha:

  Miscoverage level of the evaluated object.

## Details

The diagnostic suite performs five rigorous audits on the prediction
intervals:

1.  **Marginal Validity & Sharpness**: Evaluates empirical coverage
    \\\frac{1}{m}\sum \mathbf{1}(Y_i \in C(s_i))\\, mean/median interval
    width, and the strictly proper **Winkler Interval Score (WIS)**:
    \$\$\text{WIS}\_\alpha(l, u, y) = (u - l) + \frac{2}{\alpha}(l -
    y)\mathbf{1}(y \< l) + \frac{2}{\alpha}(y - u)\mathbf{1}(y \> u)\$\$

2.  **Conditional Spatial Stratification**: Partitions the 2D spatial
    domain into \\\text{n\\bins} \times \text{n\\bins}\\ equal-area
    quadrants and assesses local coverage.

3.  **Convex Hull Boundary Proximity**: Assesses whether edge-effect
    extrapolations degrade coverage near the boundary of the spatial
    domain.

4.  **Spatial Autocorrelation (Moran's I)**: Evaluates whether
    prediction miscoverages or nonconformity scores exhibit residual
    spatial clustering via Moran's \\I\\.

5.  **Score Distribution**: Analyzes empirical quantiles and
    normality/QQ structure.

## References

Winkler, R. L. (1972). "A Decision-Theoretic Approach to Interval
Estimation." *Journal of the American Statistical Association*, 67(337),
187-191.

Mao, R., Martin, R., and Reich, B. J. (2023). "Valid Conformal
Prediction for Dependent Data." *Journal of the American Statistical
Association*,
[doi:10.1080/01621459.2022.2147531](https://doi.org/10.1080/01621459.2022.2147531)
.

## See also

[`scp_geostatistical`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md),
[`scp_areal`](https://amjed-droid.github.io/spconform/reference/scp_areal.md),
[`coverage_report`](https://amjed-droid.github.io/spconform/reference/coverage_report.md)

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
#>              spconform Comprehensive Diagnostic Audit Report          
#> ======================================================================
#> 
#> >> 1. Marginal Validity & Prediction Sharpness:
#>    * [PASS] Empirical Coverage : 0.900 (Target Nominal >= 0.900)
#>    * Mean Interval Width    : 1.3169 (Median = 1.3565, SD = 0.0559)
#>    * Winkler Interval Score : 1.4317 (Strictly Proper Loss)
#>    * Total Target Units     : n = 30 (Covered = 27, Miscovered = 3)
#> 
#> >> 2. Spatial Residual Autocorrelation (Moran's I Audit):
#>    * [NOTE] Moran's I Statistic: 0.1376 (Expected = -0.0345, z = 5.16, p-value = 0.0000)
#>    * Conclusion: Moderate spatial residual structure detected; localized weights active.
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
```
