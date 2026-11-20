# Annotate features as MS2 fragment ions (Tier 1c / 1c_exp)

Dispatched on a
[`PeakAnnotationISF-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF-class.md)
object. link to
[`run_isf_annotation()`](https://ruibings.github.io/MS1Match_Rpackage/reference/run_isf_annotation.md).

**Phase 1 — Intra-group MS2 fragment (Tier 1c)** For each Lib_group, the
direct library-match anchor (`reference_tier == "Tier1a_AnnoMatch"`) is
used as the ISF parent, provided it has an MS2 spectrum. The anchor may
carry any confirmed ESI adduct. Every other co-group feature is tested:
if its `mzmed` matches any fragment ion within `ppm_isf` and
`rt_isf_threshold` seconds, it receives `Tier1c_MS2Fragment`.

**Phase 2 — RT-gated MS2 fragment (Tier 1c_exp)** The same anchor
feature searches the full feature pool (features with no Phase 1
assignment) within `rt_isf_threshold` seconds. Matching features get
`Tier1c_MS2Fragment_exp`.

## Usage

``` r
# S4 method for class 'PeakAnnotationISF'
annotateISF(object, ...)
```

## Arguments

- object:

  A
  [`PeakAnnotationISF-class`](https://ruibings.github.io/MS1Match_Rpackage/reference/PeakAnnotationISF-class.md)
  object.

- ...:

  Currently unused.

## Value

`data.frame` with columns: `feature_name`, `mz`, `rt`, `parent_feature`,
`parent_mz`, `parent_rt`, `matched_fragment_mz`, `isf_ppm`, `isf_mda`,
`rt_diff`, `Lib_group`, `reference_tier` (`"1c"` or `"1c_exp"`),
`annotation_phase` (`"Phase1_MS2frag"` or `"Phase2_MS2frag"`).

## See also

[`createPeakAnnotationISF()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationISF.md),
[`run_isf_annotation()`](https://ruibings.github.io/MS1Match_Rpackage/reference/run_isf_annotation.md),
[`read_ms2_consensus()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_ms2_consensus.md)
