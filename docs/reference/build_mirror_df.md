# Build the mirror plot data frame for one spectrum pair

Query intensities are normalized to 0-100; reference intensities are
used as-is (already 0-100 from the MSP library).

## Usage

``` r
build_mirror_df(
  query_mz,
  query_int,
  query_feature_names,
  ref_mz_matched,
  ref_mz_full,
  ref_int_full,
  annot_lookup,
  group_ion_type,
  scale_method = "none",
  unmatched_label = "Unmatched",
  precursor_mz = NA_real_
)
```

## Arguments

- query_mz:

  Numeric vector of matched query feature m/z

- query_int:

  Numeric vector of matched query feature intensities (raw scale)

- query_feature_names:

  Character vector of feature IDs

- ref_mz_matched:

  Numeric vector of matched library peak m/z

- ref_mz_full:

  Numeric vector of all library peak m/z

- ref_int_full:

  Numeric vector of all library peak intensities (0-100 scale)

- annot_lookup:

  Named list from
  [`build_annot_lookup()`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_annot_lookup.md)

- group_ion_type:

  Character; group-level Ion_type from direct_match

- scale_method:

  One of `"none"`, `"sqrt"`, `"log10"`

- unmatched_label:

  Label for unmatched reference peaks

- precursor_mz:

  Numeric scalar; used to identify the precursor peak and restrict
  `[M+H]+` labelling to it. Pass `NA_real_` for MS1.

## Value

Data frame carrying attribute `scale_method`

## Details

NL labels come from the pipeline `adduct_match` field only (e.g.
`"[M+H-H2O]+ 274.131742"` -\> `"[M+H-H2O]+"`). The precursor peak is
never relabelled.
