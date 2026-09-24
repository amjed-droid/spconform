# Coerce an spconform Object to a Data Frame

Converts a conformal prediction output object of class `"spconform"`
into a tidy `data.frame` containing coordinates (if present), point
predictions, lower/upper bounds, and interval widths.

## Usage

``` r
# S3 method for class 'spconform'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `"spconform"`.

- row.names:

  `NULL` or a character vector giving row names.

- optional:

  Logical; passed to `as.data.frame`.

- ...:

  Additional arguments (currently unused).

## Value

A `data.frame` with prediction columns `pred`, `lower`, `upper`,
`width`, and any coordinate columns.
