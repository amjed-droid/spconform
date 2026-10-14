# Comprehensive Diagnostic Report for spconform Objects

Produces a multi-panel diagnostic report assessing marginal coverage,
conditional coverage by spatial strata, boundary effects, and the
distribution of nonconformity scores.

## Usage

``` r
diagnose(object, y_true, s_test = NULL, n_bins = 4, plot = TRUE, ...)
```

## Arguments

- object:

  An object of class `"spconform"`, typically the output of
  [`scp_geostatistical()`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md)
  or
  [`scp_areal()`](https://amjed-droid.github.io/spconform/reference/scp_areal.md).

- y_true:

  Numeric vector of true response values at the prediction locations.
  Must have the same length as `object$pred`.

- s_test:

  Optional numeric matrix of prediction coordinates (one row per
  location). Required for spatial diagnostics.

- n_bins:

  Integer; number of spatial bins for conditional coverage (default 4).

- plot:

  Logical; if `TRUE` (default), generates a multi-panel diagnostic plot.

- ...:

  Additional arguments passed to plotting functions.

## Value

A list (invisibly) containing:

- marginal:

  Marginal coverage and mean width.

- conditional:

  Coverage and width by spatial bin.

- boundary:

  Coverage by distance from convex hull boundary.

- scores:

  Summary of nonconformity score distribution.

## Details

The conditional coverage analysis partitions the prediction locations
into `n_bins` equal-area spatial quadrants and reports coverage within
each. The boundary analysis classifies points by their distance to the
convex hull of the training data (for geostatistical output).

If empirical coverage exceeds the nominal target by more than 5
percentage points, a message is issued suggesting bandwidth reduction
for tighter intervals.

## Examples

``` r
# Minimal reproducible example (< 0.1s execution time)
set.seed(123)
s_tr <- matrix(runif(40), ncol = 2)
y_tr <- rnorm(20)
s_te <- matrix(runif(20), ncol = 2)
y_te <- rnorm(10)

pfun <- function(s_train, y_train, s_new) rep(mean(y_train), nrow(s_new))

out <- scp_geostatistical(s_tr, y_tr, s_te, pfun, alpha = 0.1)
diag_res <- diagnose(out, y_te, s_te, plot = FALSE)
#> Note: Empirical coverage (1) exceeds nominal (0.9) by >5%. Consider reducing 'bandwidth' for tighter intervals.
print(diag_res)
#> === spconform Diagnostic Report ===
#> 
#> Marginal coverage:
#>   Empirical: 1  (nominal: 0.9 )
#>   Mean width: 3.9064 
#>   n = 10 , covered = 10 
#> 
#> Conditional coverage by spatial bin:
#>   Q1-1: 1 (n=1, width=3.906)
#>   Q1-2: 1 (n=1, width=3.906)
#>   Q1-3: 1 (n=1, width=3.906)
#>   Q2-1: 1 (n=1, width=3.906)
#>   Q2-2: 1 (n=1, width=3.906)
#>   Q3-4: 1 (n=1, width=3.906)
#>   Q4-1: 1 (n=1, width=3.906)
#>   Q4-3: 1 (n=1, width=3.906)
#>   Q4-4: 1 (n=2, width=3.906)
#> 
#> Boundary effect:
#>   Near boundary:   1 (n=5)
#>   Far from boundary:1 (n=5)
#> 
#> Nonconformity scores:
#>   Mean: 1.9532 
#>   Median: 1.9532 
#>   SD: 0 
#>   90% quantile: 1.9532 
```
