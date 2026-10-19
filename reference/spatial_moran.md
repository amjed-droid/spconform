# Spatial Autocorrelation Diagnostics (Moran's I) for Conformal Residuals

Evaluates spatial autocorrelation and clustering among prediction
residuals and binary miscoverage indicators using Moran's \\I\\
statistic.

## Usage

``` r
spatial_moran(object, y_true, s_test)
```

## Arguments

- object:

  An object of class `"spconform"` or numeric vector of predictions.

- y_true:

  Numeric vector of true observed responses.

- s_test:

  Numeric matrix of spatial coordinates (\\m \times p\\).

## Value

An S3 object of class `"spconform_moran"` containing:

- moran_i:

  Observed Moran's I statistic for prediction residuals.

- expected:

  Theoretical expectation \\-1/(n-1)\\ under spatial independence.

- z_score:

  Standardized normal z-score under the randomization assumption.

- p_value:

  Two-sided p-value against the null hypothesis of spatial randomness.

- is_random:

  Logical; `TRUE` if p-value \>= 0.05.

- moran_i_miscoverage:

  Moran's I on binary hit/miss indicators if available.

## References

Cliff, A. D. and Ord, J. K. (1981). *Spatial Processes: Models &
Applications*. Pion, London.

Mao, H., Martin, R., and Reich, B. J. (2024). Valid Conformal Prediction
for Dependent Data. *Journal of the American Statistical Association*,
119(546), 1269–1281.
[doi:10.1080/01621459.2022.2147531](https://doi.org/10.1080/01621459.2022.2147531)
.

## Examples

``` r
set.seed(42)
n <- 50
s <- matrix(runif(n * 2), n, 2)
y <- sin(3 * s[, 1]) + rnorm(n, sd = 0.2)
pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
fit <- scp_geostatistical(s[1:35, ], y[1:35], s[36:50, ], pfun, alpha = 0.1)
m_res <- spatial_moran(fit, y[36:50], s[36:50, ])
print(m_res)
#> Spatial Autocorrelation Audit (Moran's I on Residuals)
#> ----------------------------------------------------
#> Observed Moran's I : +0.2147
#> Expected Value      : -0.0714
#> Standard Normal Z   : +4.2817
#> Two-sided p-value   : 0.0000 (Significant spatial clustering detected, p < 0.05)
```
