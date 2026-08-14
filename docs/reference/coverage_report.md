# Empirical Coverage and Average Interval Width for an spconform Object

Computes the empirical coverage rate and mean interval width by
comparing conformal prediction intervals against held-out ground truth
responses.

## Usage

``` r
coverage_report(object, y_true)
```

## Arguments

- object:

  An object of class `"spconform"`.

- y_true:

  Numeric vector of true observed values at the prediction locations.
  Must have the same length as `object$pred`.

## Value

A named list with:

- coverage:

  Proportion of true values falling within the prediction intervals.

- mean_width:

  Mean width of the prediction intervals.

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
