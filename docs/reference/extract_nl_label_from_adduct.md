# Extract a clean neutral-loss label from a pipeline adduct_match string

Strips trailing mass values and returns the NL portion only when the
adduct is a neutral-loss form. `"[M+H-H2O]+ 274.131742"` -\>
`"[M+H-H2O]+"` `"[M+H]+ 194.057909"` -\> `""` (plain adduct, not NL)

## Usage

``` r
extract_nl_label_from_adduct(x)
```

## Arguments

- x:

  Character scalar, possibly with trailing mass

## Value

Character scalar: the NL label (stripped) or `""`
