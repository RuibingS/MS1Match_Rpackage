# Perform primary-ion matching against a direct target list

Matches feature table m/z values against a compound target list using
`[M+H]+` (positive mode) or `[M-H]-` (negative mode) ion m/z.

The target list is automatically normalised via
[`normalize_target_list()`](https://ruibings.github.io/MS1Match_Rpackage/reference/normalize_target_list.md)
before matching, so input data frames from different sources (HMDB,
plate-based, PA14, etc.) are accepted without manual column renaming.

## Usage

``` r
PImatch_target(object)
```

## Arguments

- object:

  A `DirectTargetMatch` object

## Value

The same `DirectTargetMatch` object with `@results` populated
