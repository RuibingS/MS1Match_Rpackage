# Batch-export mirror plots for all compounds in a direct_match list

Batch-export mirror plots for all compounds in a direct_match list

## Usage

``` r
export_mirror_plots(
  input_list_ms1,
  input_list_ms2,
  annot_res_ms1,
  annot_res_ms2,
  out_dir = "mirror_plots",
  scale_method = "log10",
  mz_digits = 1,
  padding = 0.05,
  precursor_extra_frac = 0.04,
  shared_xlim = TRUE,
  mz_label_mode = "all",
  width = 6.5,
  height = 9,
  dpi = 300
)
```

## Arguments

- input_list_ms1:

  List; `direct_match_StM14_res_MS1`

- input_list_ms2:

  List; `direct_match_StM14_res_MS2`

- annot_res_ms1:

  Data frame; `StM14_annotate_res_ms1` (only
  `annotation_phase == "Phase1"` rows are used)

- annot_res_ms2:

  Data frame; `StM14_annotate_res_ms2` (only
  `annotation_phase == "Phase1"` rows are used)

- out_dir:

  Output directory (created if absent)

- scale_method:

  Intensity scale

- mz_digits:

  Decimal places for m/z labels

- padding:

  Fractional x-axis padding

- precursor_extra_frac:

  Extra right-side padding for precursor marker

- shared_xlim:

  Logical; share x-axis across MS1 and MS2 panels

- mz_label_mode:

  One of `"all"`, `"primary"`, `"adduct"`

- width:

  Plot width in inches

- height:

  Plot height in inches

- dpi:

  Resolution for PNG output

## Value

Invisibly, a list of per-compound output paths
