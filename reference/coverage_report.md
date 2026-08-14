# Empirical coverage and average interval width for an spconform object

Empirical coverage and average interval width for an spconform object

## Usage

``` r
coverage_report(object, y_true)
```

## Arguments

- object:

  An object of class `"spconform"`.

- y_true:

  Numeric vector of true observed values at the prediction locations,
  same length/order as `object$pred`.

## Value

A named list with `coverage` (proportion of `y_true` falling within
`[lower, upper]`) and `mean_width` (average interval width).

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
```
