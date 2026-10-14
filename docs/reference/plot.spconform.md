# Plot Prediction Intervals for spconform Objects

Plots point predictions and conformal prediction intervals across
observation indices, optionally overlaying true target values for visual
evaluation.

## Usage

``` r
# S3 method for class 'spconform'
plot(x, y_true = NULL, ...)
```

## Arguments

- x:

  An object of class `"spconform"`.

- y_true:

  Optional numeric vector of true observed responses at prediction
  locations.

- ...:

  Further graphical parameters passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns the input object `x`.
