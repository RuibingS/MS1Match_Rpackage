# Import PAMD Metabolite Data and Calculate Adducts

Reads a CSV file, removes all-NA columns, computes exact mass from
formula, and calculates adduct m/z values using MetaboCoreUtils.

## Usage

``` r
PAMD_target_list_import(target_list_path, IonPolarity)
```

## Arguments

- target_list_path:

  Path to the target list CSV file containing metabolite data.

## Value

A data.frame with original data, exact masses, and adduct m/z.
