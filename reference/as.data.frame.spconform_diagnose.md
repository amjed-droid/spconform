# Coerce an spconform_diagnose Object to a Data Frame

Extracts conditional spatial strata and diagnostic statistics as a tidy
data frame for programmatic inspection.

## Usage

``` r
# S3 method for class 'spconform_diagnose'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `"spconform_diagnose"`.

- row.names:

  Optional row names.

- optional:

  Logical; passed to `as.data.frame`.

- ...:

  Additional arguments (unused).

## Value

A tidy `data.frame` of spatial diagnostic strata.
