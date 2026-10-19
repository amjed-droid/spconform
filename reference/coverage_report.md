# Empirical Coverage and Average Interval Width for an spconform Object

Evaluates empirical coverage, average prediction interval width, and the
strictly proper Winkler interval score for an `"spconform"` object.

## Usage

``` r
coverage_report(object, y_true, na.rm = TRUE)
```

## Arguments

- object:

  An object of class `"spconform"`.

- y_true:

  Numeric vector of true observed values at the prediction locations,
  same length/order as `object$pred`.

- na.rm:

  Logical; if `TRUE`, missing values in `y_true` are removed before
  computing evaluation metrics. Defaults to `TRUE`.

## Value

A named list with:

- coverage:

  Proportion of observations falling within `[lower, upper]`

.

- mean_width:

  Average prediction interval width

.

- winkler_score:

  Mean Winkler interval score across evaluation points

.

## References

Winkler, R. L. (1972). A Decision-Theoretic Approach to Interval
Estimation. *Journal of the American Statistical Association*, 67(337),
187–191.

## Examples

``` r
set.seed(1)
n <- 20
adj <- matrix(0, n, n)
for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
y <- cumsum(rnorm(n))
out <- scp_areal(y, adjacency = adj, alpha = 0.2)
coverage_report(out, y)
#> $coverage
#> [1] 0.75
#> 
#> $mean_width
#> [1] 1.697122
#> 
#> $winkler_score
#> [1] 2.849409
#> 
```
