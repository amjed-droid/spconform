# Plot Prediction Intervals for spconform Objects

Visualizes the conformal prediction intervals and point predictions
across target locations or areal units. Optionally overlays true
observations if provided.

## Usage

``` r
# S3 method for class 'spconform'
plot(x, y_true = NULL, ...)
```

## Arguments

- x:

  An object of class `"spconform"`.

- y_true:

  Optional numeric vector of true observed responses matching `x$pred`.

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns the input object `x`.
