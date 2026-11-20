# Grouping table from a RAMClustR result

Grouping table from a RAMClustR result

## Usage

``` r
ramclustr_to_grouping(RC)
```

## Arguments

- RC:

  Output of
  [`RAMClustR::ramclustR()`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html),
  run on the same feature table as MS1Match (e.g. the xcms object).
  `RC$frt` then has the same RT unit as that table, i.e. seconds, and is
  used as it is.

## Value

See
[`prepare_grouping_table`](https://ruibings.github.io/MS1Match_Rpackage/reference/prepare_grouping_table.md).
Features with `featclus = 0` (not clustered) get `group = NA`.
