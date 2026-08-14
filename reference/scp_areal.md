# Neighbourhood-Weighted Conformal Prediction for Areal (Lattice) Data

Constructs distribution-free prediction intervals for areal units (e.g.
counties, census tracts, grid cells) using a leave-one-unit-out
conformal procedure in which nonconformity scores are weighted by graph
(neighbourhood) distance to the held-out unit, via its adjacency
structure.

## Usage

``` r
scp_areal(y, X = NULL, adjacency, pred_fun = NULL, alpha = 0.1, decay = 1)
```

## Arguments

- y:

  Numeric vector of observed responses, one per areal unit.

- X:

  Optional numeric design matrix of covariates (one row per areal unit);
  passed to `pred_fun` if supplied. May be `NULL` for purely spatial
  (neighbourhood-based) prediction.

- adjacency:

  Square adjacency (contiguity) matrix describing the neighbourhood
  structure among the `length(y)` areal units. Non-zero entries are
  treated as neighbours (coerced to binary).

- pred_fun:

  A function with signature
  `function(y_train, X_train, idx_train, idx_target, adjacency)`
  returning a single numeric point prediction for the areal unit indexed
  by `idx_target`, fitted using the training units `idx_train`. If
  `NULL` (default), a simple neighbourhood-mean predictor is used (the
  mean of `y_train` over the target's graph neighbours, falling back to
  the global training mean if the target has no observed neighbours).

- alpha:

  Miscoverage level; intervals target \\1-\alpha\\ coverage. Default
  0.1.

- decay:

  Decay rate for neighbourhood weights; see
  [`areal_neighbor_weights`](https://amjed-droid.github.io/spconform/reference/areal_neighbor_weights.md).

## Value

An object of class `"spconform"` (areal variant) with the same structure
as
[`scp_geostatistical`](https://amjed-droid.github.io/spconform/reference/scp_geostatistical.md),
using unit indices in place of coordinates.

## Details

Uses a full leave-one-out (jackknife-style) conformal scheme: for each
areal unit \\i\\, a model is fit on all other units and used to predict
unit \\i\\; the resulting nonconformity scores across all units are
combined into a graph-distance-weighted quantile specific to each target
unit, so that the calibration set is dominated by
spatially/graph-proximate units rather than treating all units as
exchangeable.

If the effective neighbourhood weight for a target unit is vanishingly
small (e.g., an isolated unit with no neighbours), the procedure falls
back to an unweighted quantile over all other units and issues a
warning.

## Examples

``` r
set.seed(1)
n <- 25
adj <- matrix(0, n, n)
for (i in 1:(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }
y <- cumsum(rnorm(n)) + rnorm(n, sd = 0.2)

out <- scp_areal(y, adjacency = adj, alpha = 0.2)
print(out)
#> <spconform> areal conformal prediction
#> Target coverage: 80.0%
#> Number of prediction points: 25
#>     pred  lower upper
#> 1 -0.474 -1.920 0.972
#> 2 -1.105 -2.551 0.341
#> 3 -0.126 -0.769 0.516
#> 4 -0.421 -1.868 1.025
#> 5  0.159 -0.483 0.802
#> 6  0.511 -0.131 1.154
#> ... (19 more)
```
