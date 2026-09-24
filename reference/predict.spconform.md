# Extract Predictions and Intervals from an spconform Object

Extracts point predictions or prediction interval bounds from a fitted
`spconform` object.

## Usage

``` r
# S3 method for class 'spconform'
predict(object, interval = c("none", "prediction"), ...)
```

## Arguments

- object:

  An object of class `"spconform"`.

- interval:

  Character string specifying the type of prediction intervals to
  extract: `"none"` (default, point predictions only) or `"prediction"`
  (matrix with fit, lwr, upr).

- ...:

  Further arguments passed to or from other methods.

## Value

If `interval = "none"`, a numeric vector of point predictions. If
`interval = "prediction"`, a numeric matrix with columns `fit`, `lwr`,
and `upr`.
