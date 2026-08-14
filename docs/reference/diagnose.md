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
if (FALSE) { # \dontrun{
out <- scp_geostatistical(s_train, y_train, s_test, pred_fun, alpha = 0.1)
diag <- diagnose(out, y_test, s_test)
} # }
```
