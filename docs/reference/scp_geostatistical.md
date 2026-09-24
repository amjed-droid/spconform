# Locally Weighted Split Conformal Prediction for Geostatistical (Point-Referenced) Data

Constructs distribution-free prediction intervals at new spatial (or
spatio-temporal) locations by combining split conformal prediction with
spatial-distance kernel weights, relaxing the exchangeability assumption
that standard conformal prediction relies on.

## Usage

``` r
scp_geostatistical(
  s_train,
  y_train,
  s0,
  pred_fun,
  alpha = 0.1,
  split = 0.5,
  bandwidth = NULL,
  t_train = NULL,
  t0 = NULL,
  temporal_bandwidth = NULL,
  seed = NULL
)
```

## Arguments

- s_train:

  Numeric matrix of training coordinates, one row per observation (e.g.
  longitude/latitude, or projected easting/northing).

- y_train:

  Numeric vector of training responses, same length as `nrow(s_train)`.

- s0:

  Numeric matrix of coordinates at which prediction intervals are
  desired, one row per target location.

- pred_fun:

  A function with signature `function(s_train, y_train, s_new)` that
  returns a numeric vector of point predictions for the locations in
  `s_new`, fitted using `s_train`/`y_train`. This lets the user plug in
  any spatial predictor (kriging, random forest, GAM, etc.); `spconform`
  handles only the conformal calibration layer.

- alpha:

  Miscoverage level; prediction intervals target \\1 - \alpha\\
  coverage. Default 0.1 (90 percent intervals).

- split:

  Proportion of the training data used for model fitting; the remainder
  is used as the calibration set for conformal scores. Default 0.5.

- bandwidth:

  Optional spatial kernel bandwidth passed to
  [`spatial_kernel_weights`](https://amjed-droid.github.io/spconform/reference/spatial_kernel_weights.md);
  if `NULL`, chosen automatically per target location.

- t_train, t0, temporal_bandwidth:

  Optional time indices (and temporal bandwidth) for spatio-temporal
  weighting; see
  [`spatial_kernel_weights`](https://amjed-droid.github.io/spconform/reference/spatial_kernel_weights.md).

- seed:

  Optional integer seed for the train/calibration split, for
  reproducibility.

## Value

An object of class `"spconform"`, a list with components:

- pred:

  Point predictions at `s0`.

- lower, upper:

  Lower/upper prediction interval bounds at `s0`.

- alpha:

  The requested miscoverage level.

- s0:

  The target coordinates.

- call:

  The matched call.

## Details

This implements a localized split-conformal procedure in the spirit of
Mao, Martin and Reich (2020, "Valid model-free spatial prediction") and
subsequent spatio-temporal extensions: nonconformity scores from a
held-out calibration set are combined into a weighted empirical quantile
where weights decay with distance (in space, and optionally time) from
the target location, so nearby, more relevant observations dominate the
calibration of the interval while validity is retained under a local
exchangeability argument.

## References

Mao, H., Martin, R., and Reich, B. J. (2020). Valid model-free spatial
prediction. *arXiv:2006.15640*.

## Examples

``` r
set.seed(42)
n <- 200
s <- matrix(runif(2 * n), ncol = 2)
f <- function(s) sin(4 * s[, 1]) + cos(3 * s[, 2])
y <- f(s) + rnorm(n, sd = 0.3)

idx <- sample(n, 150)
s_train <- s[idx, ]; y_train <- y[idx]
s0 <- s[-idx, ][1:10, ]

pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  nd <- data.frame(s_train = I(s_new))
  as.numeric(cbind(1, s_new[, 1], s_new[, 2], s_new[, 1]^2, s_new[, 2]^2) %*%
               coef(fit))
}

out <- scp_geostatistical(s_train, y_train, s0, pred_fun, alpha = 0.1, seed = 1)
print(out)
#> <spconform> geostatistical conformal prediction
#> Target coverage: 90.0%
#> Number of prediction points: 10
#>     pred  lower  upper
#> 1 -0.543 -1.075 -0.010
#> 2  1.060  0.495  1.625
#> 3 -0.816 -1.355 -0.277
#> 4  0.619  0.054  1.184
#> 5  1.091  0.526  1.656
#> 6 -1.373 -1.906 -0.841
#> ... (4 more)
```
