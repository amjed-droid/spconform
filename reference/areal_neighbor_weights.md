# Compute neighbourhood-based weights for areal (lattice) data

Compute neighbourhood-based weights for areal (lattice) data

## Usage

``` r
areal_neighbor_weights(i0, adjacency, decay = 1)
```

## Arguments

- i0:

  Integer index of the target (unobserved / held-out) areal unit.

- adjacency:

  A square 0/1 (or weighted) adjacency matrix describing the
  neighbourhood structure of the areal units (e.g. a spatial contiguity
  matrix). Row/column `i0` corresponds to the target unit.

- decay:

  Numeric decay rate applied to graph distance (number of hops) from
  `i0`; larger values down-weight distant neighbours more aggressively.
  Defaults to 1.

## Value

A numeric vector of length `nrow(adjacency)` with weights based on graph
distance from `i0` (self-weight is 0, i.e. the target unit is excluded
from its own calibration set).
