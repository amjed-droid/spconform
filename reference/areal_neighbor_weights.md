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

## Examples

``` r
adj <- matrix(c(0, 1, 0,
                1, 0, 1,
                0, 1, 0), nrow = 3, byrow = TRUE)
w <- areal_neighbor_weights(1, adj, decay = 0.5)
print(w)
#> [1] 0.0000000 0.6065307 0.3678794
```
