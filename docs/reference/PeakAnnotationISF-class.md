# PeakAnnotationISF S4 class

Container for MS2 fragment-based annotation of an MS1 feature table.

**Phase 1 — Intra-group MS2 fragment (Tier 1c)** For each confirmed
parent feature (identified by `parent_col` == `parent_val`) with a
consensus MS2 spectrum, every other feature in the same library group is
checked: if its `mzmed` matches any fragment ion within `ppm_isf` and
`rt_isf_threshold` seconds, it receives Tier 1c.

**Phase 2 — RT-gated MS2 fragment (Tier 1c_exp)** The full feature pool
(features not in any Phase 1 group) is searched within
`rt_isf_threshold` seconds of the parent. Matching features receive Tier
1c_exp.

## Slots

- `phase1_result`:

  data.frame. Output of `annotate()`.

- `ms2_list`:

  list. Named list from
  [`read_ms2_consensus()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_ms2_consensus.md).

- `FT`:

  data.frame. Full MS1 feature pool.

- `parent_col`:

  character. Column identifying parent features. Default `"PI_match"`.

- `parent_val`:

  ANY. Value marking a parent. `NULL` = any non-empty string. Default
  `NULL`.

- `ppm_isf`:

  numeric. ppm tolerance for fragment m/z matching. Default `10`.

- `mz_threshold_isf`:

  numeric. Absolute m/z floor (Da). Default `0.005`.

- `rt_isf_threshold`:

  numeric. RT window (seconds). Default `5`.
