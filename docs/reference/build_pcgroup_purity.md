# Assess CAMERA pcgroup purity (over-grouping)

For every pcgroup that contains at least one anchor feature, counts how
many distinct compounds (neutral masses) were merged into the same
pcgroup by CAMERA. A pure pcgroup contains exactly one compound.

## Usage

``` r
build_pcgroup_purity(CAMERA_ft, eval_full, anchors)
```

## Arguments

- CAMERA_ft:

  `data.frame`. CAMERA-annotated feature table. Must contain:
  `feature_name`, `pcgroup`.

- eval_full:

  `data.frame`. Output of
  [`build_eval_full`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_eval_full.md).

- anchors:

  `data.frame`. Output of
  [`build_anchors`](https://ruibings.github.io/MS1Match_Rpackage/reference/build_anchors.md).

## Value

A named `list` with two elements:

- `all_pcgroups`:

  Per-pcgroup summary: `pcgroup`, `pcgroup_n_features`,
  `pcgroup_n_compounds`, `pcgroup_outside_ref`, `median_rt`.

- `anchor_pcgroups`:

  Subset of `all_pcgroups` restricted to pcgroups that contain at least
  one anchor feature.
