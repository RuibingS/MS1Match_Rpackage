# Calculate monoisotopic exact mass from molecular formulae

Fast C++ implementation of monoisotopic mass calculation. No external
packages required — replaces the previous `MetaboCoreUtils` dependency.

## Usage

``` r
calculate_exact_mass(formulas)
```

## Arguments

- formulas:

  Character vector of molecular formulae (e.g.
  `c("C6H12O6", "C8H10N4O2")`).

## Value

Numeric vector of monoisotopic masses, in the same order as `formulas`.
Returns `NA` for unknown elements or unparseable formulae.

## Examples

``` r
calculate_exact_mass(c("C6H12O6", "C8H10N4O2", "C31H36N2O11"))
#> [1] 180.0634 194.0804 612.2319
# [1] 180.0634 194.0798 592.2375

calculate_exact_mass("C9H17NO6")
#> [1] 235.1056
# [1] 235.1056
```
