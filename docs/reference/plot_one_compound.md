# Plot MS1 and MS2 mirror panels stacked for one compound

Plot MS1 and MS2 mirror panels stacked for one compound

## Usage

``` r
plot_one_compound(
  ms1_entry,
  ms2_entry,
  annot_res_ms1,
  annot_res_ms2,
  scale_method = "log10",
  mz_digits = 1,
  padding = 0.05,
  precursor_extra_frac = 0.04,
  shared_xlim = TRUE,
  mz_label_mode = "all"
)
```

## Arguments

- ms1_entry:

  One element of `direct_match_res_MS1`

- ms2_entry:

  One element of `direct_match_res_MS2`

- annot_res_ms1:

  Full `StM14_annotate_res_ms1` data frame. Filtered internally to
  `annotation_phase == "Phase1"`.

- annot_res_ms2:

  Full `StM14_annotate_res_ms2` data frame. Filtered internally to
  `annotation_phase == "Phase1"`.

- scale_method:

  Intensity scale: `"none"`, `"sqrt"`, `"log10"`

- mz_digits:

  Decimal places for m/z labels

- padding:

  Fractional x-axis padding

- precursor_extra_frac:

  Extra right-side padding for precursor marker

- shared_xlim:

  Logical; share x-axis across both panels

- mz_label_mode:

  One of `"all"`, `"primary"`, `"adduct"`

## Value

A `patchwork` object
