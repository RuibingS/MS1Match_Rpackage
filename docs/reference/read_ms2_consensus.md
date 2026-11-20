# Read MS2 consensus data into a named list

Converts MS2 consensus data into a named list keyed by `feature_name`,
ready for
[`run_isf_annotation()`](https://ruibings.github.io/MS1Match_Rpackage/reference/run_isf_annotation.md).

Three input types accepted:

- Option A:

  `Spectra` from `combineSpectra()` — empty spectra removed
  automatically.

- Option B:

  `Spectra` already filtered
  (`ms2_final <- ms2_consensus[lengths(mz(ms2_consensus)) > 0]`).

- Option C:

  Path to an MGF file with `feature_name=` headers.

When `feature_table` is provided, the MS2 scan RT (DDA ) is replaced by
the MS1 peak apex RT (`rtmed`) from the feature table. This unifies MS2
and MS1 RT at source, preventing spurious RT differences of up to ~20 s
caused by DDA firing on the rising edge of a chromatographic peak.

## Usage

``` r
read_ms2_consensus(x, feature_table = NULL, ppm_check = 20, verbose = TRUE)
```

## Arguments

- x:

  `Spectra` object or MGF file path.

- feature_table:

  `data.frame`. Optional MS1 feature table with columns `feature_name`
  and `mzmed`/`mz` and `rtmed`/`rt`. When provided:

  - MS2 RT is overridden with MS1 apex RT.

  - Precursor m/z is cross-checked within `ppm_check`.

- ppm_check:

  Numeric. ppm tolerance for precursor m/z check. Default `20`.

- verbose:

  Logical. Print summary. Default `TRUE`.

## Value

Named list keyed by `feature_name`. Each element:

- `feature_name`:

  Character.

- `precursor_mz`:

  Numeric.

- `rt`:

  Numeric. MS1 apex RT (seconds) if `feature_table` provided; otherwise
  MS2 scan RT.

- `fragment_mz`:

  Numeric vector.

- `fragment_int`:

  Numeric vector.

- `n_fragments`:

  Integer.

## Examples

``` r
if (FALSE) { # \dontrun{
# Option A — ms2_consensus directly (empty spectra removed inside)
ms2_list <- read_ms2_consensus(ms2_consensus, feature_table = pool_ft)

# Option B — pre-filtered ms2_final
ms2_final <- ms2_consensus[lengths(mz(ms2_consensus)) > 0]
ms2_list  <- read_ms2_consensus(ms2_final, feature_table = pool_ft)

# Option C — MGF file
ms2_list  <- read_ms2_consensus("xcms_ms2_consensus_final.mgf",
                                 feature_table = pool_ft)
} # }
```
