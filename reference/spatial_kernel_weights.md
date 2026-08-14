# Compute Gaussian-kernel spatial weights between a target location and a set of reference locations

Compute Gaussian-kernel spatial weights between a target location and a
set of reference locations

## Usage

``` r
spatial_kernel_weights(
  s0,
  s,
  bandwidth = NULL,
  t0 = NULL,
  t = NULL,
  temporal_bandwidth = NULL
)
```

## Arguments

- s0:

  A numeric vector (or 1-row matrix) giving the coordinates of the
  target (prediction) location.

- s:

  A numeric matrix of coordinates (one row per observation) for the
  reference (calibration) set.

- bandwidth:

  Positive numeric bandwidth of the Gaussian kernel. If `NULL`
  (default), the bandwidth is chosen automatically as the median
  pairwise distance among `s` (Silverman-type rule).

- t0:

  Optional numeric scalar: time index of the target observation, for
  spatio-temporal weighting. If supplied, `t` must be supplied too.

- t:

  Optional numeric vector: time indices for the reference set, same
  length as `nrow(s)`.

- temporal_bandwidth:

  Positive numeric bandwidth for the temporal kernel; only used if
  `t0`/`t` are supplied. Defaults to the spatial bandwidth logic applied
  to `t`.

## Value

A numeric vector of weights (not necessarily summing to one), one per
row of `s`, giving higher weight to reference points that are spatially
(and, if requested, temporally) closer to `s0`.

## Details

Weights are computed as \$\$w_i = \exp(-\\s_0 - s_i\\^2 / (2 h^2))\$\$
for the spatial-only case, and multiplied by an analogous temporal
kernel term when `t0`/`t` are provided. This is the core device used to
relax the exchangeability assumption required by standard conformal
prediction: observations located near the point to be predicted
contribute more to the calibration of the prediction interval than
distant ones (Mao, Martin and Reich, 2020).

## Examples

``` r
set.seed(1)
s <- matrix(runif(20), ncol = 2)
s0 <- c(0.5, 0.5)
w <- spatial_kernel_weights(s0, s)
round(w, 3)
#>  [1] 0.805 0.831 0.940 0.759 0.780 0.784 0.687 0.663 0.953 0.662
```
