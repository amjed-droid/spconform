# Conditional Spatial Strata Diagnostics for Conformal Prediction

Partitions the spatial domain into coordinate grid strata/quadrants and
audits local empirical coverage, interval widths, and Winkler scores.

## Usage

``` r
strata_diagnostics(object, y_true, s_test, n_bins = 4)
```

## Arguments

- object:

  An object of class `"spconform"`.

- y_true:

  Numeric vector of true observed responses.

- s_test:

  Numeric matrix of spatial coordinates (\\m \times 2\\).

- n_bins:

  Integer; number of spatial bins per coordinate axis (default 4).

## Value

An S3 object of class `"spconform_strata"` containing a data frame of
strata evaluations and summary metrics.

## References

Winkler, R. L. (1972). A Decision-Theoretic Approach to Interval
Estimation. *Journal of the American Statistical Association*, 67(337),
187–191.

## Examples

``` r
set.seed(42)
n <- 50
s <- matrix(runif(n * 2), n, 2)
y <- sin(3 * s[, 1]) + rnorm(n, sd = 0.2)
pfun <- function(s_tr, y_tr, s_new) rep(mean(y_tr), nrow(s_new))
fit <- scp_geostatistical(s[1:35, ], y[1:35], s[36:50, ], pfun, alpha = 0.1)
s_res <- strata_diagnostics(fit, y[36:50], s[36:50, ], n_bins = 2)
print(s_res)
#> Conditional Spatial Strata Audit
#> --------------------------------
#>    - Strata Q1-1   : Cov = 100.0% | Mean Width =  1.401 | WIS =  1.401 (n = 4)
#>    - Strata Q1-2   : Cov = 100.0% | Mean Width =  1.659 | WIS =  1.659 (n = 2)
#>    - Strata Q2-1   : Cov = 100.0% | Mean Width =  1.142 | WIS =  1.142 (n = 3)
#>    - Strata Q2-2   : Cov = 100.0% | Mean Width =  1.142 | WIS =  1.142 (n = 6)
```
