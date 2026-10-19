# Domain Boundary Proximity Diagnostics for Conformal Prediction

Evaluates whether prediction intervals degrade near the convex-hull
boundary of the spatial prediction domain versus interior locations.

## Usage

``` r
boundary_diagnostics(object, y_true, s_test)
```

## Arguments

- object:

  An object of class `"spconform"`.

- y_true:

  Numeric vector of true observed responses.

- s_test:

  Numeric matrix of spatial coordinates (\\m \times 2\\).

## Value

An S3 object of class `"spconform_boundary"` comparing empirical
coverage, interval width, and Winkler scores across edge and interior
subdomains.

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
b_res <- boundary_diagnostics(fit, y[36:50], s[36:50, ])
print(b_res)
#> Domain Boundary Proximity Audit (Convex Hull)
#> ---------------------------------------------
#> Distance Median Threshold: 0.0000
#> Near Boundary (Edge) : n = 8 | Cov = 100.0% | Mean Width =  1.336 | WIS =  1.336
#> Far Boundary (Core) : n = 7 | Cov = 100.0% | Mean Width =  1.216 | WIS =  1.216
```
