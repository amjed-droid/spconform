# Extract Residuals from an spconform Object

Calculates prediction residuals (observed minus fitted values) or
nonconformity scores for an `spconform` object given true responses.

## Usage

``` r
# S3 method for class 'spconform'
residuals(object, y_true, type = c("response", "abs"), ...)
```

## Arguments

- object:

  An object of class `"spconform"`.

- y_true:

  Numeric vector of true observed responses at prediction locations.

- type:

  Character string indicating residual type: `"response"` (raw residuals
  \\y - \hat{y}\\) or `"abs"` (absolute nonconformity scores \\\|y -
  \hat{y}\|\\).

- ...:

  Further arguments passed to or from other methods.

## Value

A numeric vector of residuals.
