# 3. Evaluate a feature-grouping tool against MS1Match

## Overview

Feature-grouping tools such as CAMERA, RAMClustR or MZmine ion identity
networking combine the features of one compound (adducts, isotopes,
in-source fragments) into groups. Untargeted metabolomics has no
absolute ground truth for these groups. MS1Match therefore provides a
**compound-centric reference feature set**: the features it assigned to
library compounds, stratified by evidence tier. This vignette evaluates
the groups of another tool against that reference set.

**This vignette does not run MS1Match.** It starts from a finished
MS1Match result:

- for an xcms feature table, see vignette *2a. MS1Match workflow:
  StM14*;
- for an MZmine or other CSV feature table, or data without MS2, see
  vignette *2c. MS1Match with other inputs*.

MS1Match is the reference here, not one of the tools being compared.
Evaluating MS1Match’s own groups against its own reference set would be
circular.

The examples use the StM14 data shipped with the package, with CAMERA
and RAMClustR as example tools. The last section collects the calls for
your own data.

``` r

library(MS1Match)
library(dplyr)
```

## What you need

| Input | Where it comes from |
|:---|:---|
| Feature table | The feature table MS1Match was run on |
| MS1Match result | anno_result from join_annotation_results() (vignette 2a or 2c) |
| Grouping result | Your grouping tool, run on the same feature table |

**The one requirement:** MS1Match and the grouping tool must work on the
**same feature table**, so that both refer to a feature by the same id
(for xcms `FT0001`, `FT0002`, …).

``` workflow
                 ┌─> MS1Match ──────> anno_result ──┐   (vignette 2a / 2c)
                 │                                  │
feature table ──>┤                                  ├─> evaluate_grouping()
                 │                                  │   (this vignette)
                 └─> grouping tool ─> grouping ─────┘
                     (CAMERA, RAMClustR, ...)
```

The evaluation has five steps:

1.  Load the feature table and the MS1Match result.
2.  Convert the grouping result into a standard table.
3.  Check that both used the same feature table.
4.  Evaluate one or several tools.
5.  Read the results.

## Step 1 — Load the feature table and the MS1Match result

For StM14, the feature table is the xcms result `StM14_out.rds`.
MS1Match, CAMERA and RAMClustR were all run on it. `anno_result` was
saved at the end of vignette 2a.

``` r

data_dir <- system.file("extdata", "RDS_files", "StM14", package = "MS1Match")

StM14_out     <- readRDS(file.path(data_dir, "StM14_out.rds"))
feature_table <- StM14_out$ft

anno_result   <- readRDS(file.path(data_dir, "anno_result_StM14.rds"))
```

The evaluation uses two tables of `anno_result`:

- `export_feature_compound_summary`: the reference features, one row per
  feature × compound;
- `export_reduced_anno_table`: used to find each compound’s Tier 1a
  `[M+H]+` anchor.

## Step 2 — Convert the grouping result into a standard table

Every tool is reduced to one row per feature with these columns:

| Column       | Meaning                                              |
|:-------------|:-----------------------------------------------------|
| feature_name | feature id, as in the MS1Match feature table         |
| group        | group / cluster id; NA = feature not grouped         |
| mz           | m/z of the feature, used to verify the feature table |
| rt           | retention time in seconds (optional)                 |

The table must list **every** feature the tool processed, including the
ones it did not group (`group = NA`). The converters below do this.

Pick the subsection that matches your tool.

### CAMERA

Use the peak list from
[`CAMERA::getPeaklist()`](https://rdrr.io/pkg/CAMERA/man/getPeaklist-methods.html).
If it has no `feature_name` column,
[`camera_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/camera_to_grouping.md)
creates xcms-style ids from the row order. To be explicit, pass
`feature_names = rownames(featureDefinitions(xdata))`.

``` r

CAMERA_peaklist <- read.csv(file.path(data_dir, "CAMERA_annotated_StM14_FT.csv"))

grouping_camera <- camera_to_grouping(CAMERA_peaklist)
head(grouping_camera)
```

    ##   feature_name group       mz         rt
    ## 1       FT0001   798 61.00759   44.89657
    ## 2       FT0002   799 70.01280   44.89657
    ## 3       FT0003   427 72.08098   55.67931
    ## 4       FT0004   358 72.93733 1370.31300
    ## 5       FT0005   100 73.53200   62.09441
    ## 6       FT0006   800 73.53203   44.71700

### RAMClustR

Pass the object returned by
[`RAMClustR::ramclustR()`](https://rdrr.io/pkg/RAMClustR/man/ramclustR.html).
Features RAMClustR did not cluster (`featclus = 0`) get `group = NA`.

``` r

RC <- readRDS(rc_path)

grouping_ramclustr <- ramclustr_to_grouping(RC)
head(grouping_ramclustr)
```

    ##   feature_name group       mz    rt
    ## 1       FT2639  C014 585.4481 3.050
    ## 2       FT2522  C014 550.6280 3.370
    ## 3       FT0975  C014 250.1779 3.397
    ## 4       FT1058  C014 261.1457 3.397
    ## 5       FT1385  C014 299.1613 3.397
    ## 6       FT1695  C014 331.2090 3.397

### MZmine

This applies when MS1Match was run on an MZmine feature table imported
with
[`read_mzmine_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md)
(vignette 2c). Take the grouping from the same MZmine full feature
table: `feature_group` (correlation grouping) or `ion_identities:iin_id`
(ion identity networks).

[`mzmine_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/mzmine_to_grouping.md)
applies the same conversions as the import:

- MZmine `id` becomes an `FT0001`-style id;
- RT is converted from minutes to seconds;
- features removed by the import’s detection filter are dropped.

MZmine numbers its groups from 0, so `0` is a real group. Set
`unassigned = ""` so that only empty cells count as “not grouped”.

``` r

mzmine_file <- system.file("extdata", "MZmine_feature_table", "StM14_MZmine_export",
                           "StM14_mzmine_full_feature_table.csv",
                           package = "MS1Match")

ft_mzmine  <- read_mzmine_feature_table(mzmine_file)   # table MS1Match ran on
mzmine_raw <- read.csv(mzmine_file, check.names = FALSE)

grouping_mzmine <- mzmine_to_grouping(
  mzmine_raw,
  group_col     = "ion_identities:iin_id",   # or "feature_group"
  keep_features = ft_mzmine$feature_name,
  unassigned    = ""                         # MZmine group ids start at 0
)
```

### Any other tool

This covers a tool that ran on the MS1Match feature table and exports
its result in its own format, e.g. a CSV. Name its feature-id and group
columns:

``` r

my_result <- read.csv("my_tool_groups.csv")

grouping_mytool <- prepare_grouping_table(
  my_result,
  feature_col = "feature_id",   # your tool's feature-id column (FT0001, ...)
  group_col   = "cluster",      # your tool's group column
  mz_col      = "mz",           # your tool's m/z column (required for the check)
  rt_col      = "rt",           # optional
  rt_unit     = "sec",          # or "min"; converted to seconds
  unassigned  = c("0", "")      # values meaning "not grouped"; check whether
                                #   your tool uses 0 as a real group id
)
```

#### The tool’s output has no feature ids

If the tool exported m/z and RT but no feature ids, give each row the id
of the feature with the same m/z and RT in the MS1Match feature table.
The tool ran on that table, so the values are identical.

``` r

my_result <- assign_feature_names(
  my_result, feature_table,
  mz_col  = "mz",
  rt_col  = "rt",
  rt_unit = "min"            # unit of my_result$rt; converted to seconds here
)

# RT is now in seconds: keep the default rt_unit = "sec"
grouping_mytool <- prepare_grouping_table(
  my_result, group_col = "cluster", mz_col = "mz", rt_col = "rt")
```

Do not create ids by sorting the grouping result. Sorted ids only agree
with the feature table if both hold exactly the same rows in the same
order. Matching on m/z and RT also works for subsets.

#### The groups are in the feature table itself

Some pre-processing software (e.g. MZmine, MS-DIAL) exports one table
that holds both the features and a group / cluster column. That table
serves both branches:

1.  Run MS1Match on it as described in vignette 2c. This gives `ft` and
    `anno_result`. (For MZmine tables, use
    [`mzmine_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/mzmine_to_grouping.md)
    above instead.)
2.  Take the grouping from the same `ft`:

``` r

grouping_mytool <- prepare_grouping_table(
  ft,
  group_col = "cluster",                    # the group column of your table
  mz_col    = "mzmed", rt_col = "rtmed"     # column names in ft; check names(ft)
)
```

3.  Continue with Step 3, using `feature_table = ft`.

### If the feature tables differ

Two different feature tables cannot be compared. A feature in one table
has no defined counterpart in the other, so “same group as the anchor”
has no meaning. Run both on one table instead:

- **Run the grouping tool on the MS1Match feature table.** For example,
  give CAMERA, RAMClustR or cliqueMS the same xcms object that MS1Match
  used.
- **Or run MS1Match on the grouping tool’s feature table.** For example,
  import an MZmine feature table with
  [`read_mzmine_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md),
  run MS1Match on it (vignette 2c), and evaluate the MZmine grouping
  with
  [`mzmine_to_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/mzmine_to_grouping.md).

## Step 3 — Check that both used the same feature table

The check answers yes or no. Every feature id shared by the MS1Match
result, the grouping table and the feature table must have **the same
m/z** (and RT, where available). Feature ids alone are not enough,
because separate xcms runs all name their features `FT0001`, `FT0002`, …

``` r

check_feature_match(anno_result, grouping_camera, feature_table)
```

    ## Same feature table: YES 
    ##   Shared features:    3633 (m/z differ: 0, RT differ: 0)
    ##   Only in grouping table: 0 | only in feature table: 0 (excluded from evaluation)
    ##   MS1Match reference: 678 checked, 0 not in feature table, 0 m/z differ

If the answer is NO, the printout lists the first mismatching features
and the reason: m/z differs, RT differs, or a reference feature is
missing from the feature table.
[`evaluate_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md)
runs the same check and stops on NO.

**Subsets are allowed.** One table may be a subset of the other, e.g.
MS1Match run on an RT-truncated feature table and the grouping tool on
the full one. The evaluation then uses only the features present in
both. The printout reports how many were excluded on each side.

A feature missing from the grouping table is treated as outside the
shared feature table, not as Ref-unassigned.

## Step 4 — Evaluate

### One tool

``` r

eval_camera <- evaluate_grouping(
  anno_result,
  grouping_camera,
  tool_name     = "CAMERA",
  feature_table = feature_table,
  rt_window     = c(60, 1200)
)
eval_camera
```

    ## Grouping evaluation: CAMERA 
    ##   Same feature table: YES
    ##   Shared features evaluated: 3633 | reference features excluded: 0
    ##   Compound anchors: 14 
    ## 
    ##        Tier n_ref pct_grouped pct_no_anchor pct_unassigned pct_confirmed
    ##      Tier1a   129         100             0              0          95.3
    ##  Tier1a_exp     9         100             0              0          66.7
    ##      Tier1b   169         100             0              0          92.3
    ##  Tier1b_exp   113         100             0              0          95.6
    ##      Tier1c    61         100             0              0          98.4
    ##  Tier1c_exp     5         100             0              0          80.0
    ##      Tier1d    49         100             0              0          91.8
    ##       Tier2   151         100             0              0          82.8
    ##       Total   686         100             0              0          91.4
    ##  pct_conflicting
    ##              4.7
    ##             33.3
    ##              7.7
    ##              4.4
    ##              1.6
    ##             20.0
    ##              8.2
    ##             17.2
    ##              8.6

`rt_window` only restricts features *outside* the reference set. It
keeps their counts comparable between tools: here, MS1Match was run on
features within 60–1200 s. Reference features are always evaluated.

### Several tools

Every tool is evaluated with the same call; only the grouping table and
the name change. Put the grouping tables in a named list:

``` r

groupings <- list(CAMERA = grouping_camera)
if (has_rc) groupings$RAMClustR <- grouping_ramclustr
# groupings$MyTool <- grouping_mytool     # any table from Step 2

evals <- lapply(names(groupings), function(tool)
  evaluate_grouping(anno_result, groupings[[tool]],
                    tool_name     = tool,
                    feature_table = feature_table,
                    rt_window     = c(60, 1200)))
names(evals) <- names(groupings)
```

## Step 5 — Results

Each evaluation contains three tables:

| Table | Content |
|:---|:---|
| benchmark | One row per tier: counts and % of reference features per status |
| status_counts | Status counts over all features, including those outside the reference set |
| full_table | One row per reference feature × compound, plus one row per feature outside the reference set |

### Tier benchmark

``` r

show_table(eval_camera$benchmark,
           caption = "CAMERA: tier-stratified benchmark (StM14)",
           highlight_last = TRUE, scroll = TRUE)   # last row = Total
```

| Tier | n_ref | n_grouped | pct_grouped | n_confirmed | pct_confirmed | n_conflicting | pct_conflicting | n_unassigned | pct_unassigned | n_no_anchor | pct_no_anchor | n_not_evaluable |
|----|----|----|----|----|----|----|----|----|----|----|----|----|
| Tier1a | 129 | 129 | 100 | 123 | 95.3 | 6 | 4.7 | 0 | 0 | 0 | 0 | 0 |
| Tier1a_exp | 9 | 9 | 100 | 6 | 66.7 | 3 | 33.3 | 0 | 0 | 0 | 0 | 0 |
| Tier1b | 169 | 169 | 100 | 156 | 92.3 | 13 | 7.7 | 0 | 0 | 0 | 0 | 0 |
| Tier1b_exp | 113 | 113 | 100 | 108 | 95.6 | 5 | 4.4 | 0 | 0 | 0 | 0 | 0 |
| Tier1c | 61 | 61 | 100 | 60 | 98.4 | 1 | 1.6 | 0 | 0 | 0 | 0 | 0 |
| Tier1c_exp | 5 | 5 | 100 | 4 | 80.0 | 1 | 20.0 | 0 | 0 | 0 | 0 | 0 |
| Tier1d | 49 | 49 | 100 | 45 | 91.8 | 4 | 8.2 | 0 | 0 | 0 | 0 | 0 |
| Tier2 | 151 | 151 | 100 | 125 | 82.8 | 26 | 17.2 | 0 | 0 | 0 | 0 | 0 |
| Total | 686 | 686 | 100 | 627 | 91.4 | 59 | 8.6 | 0 | 0 | 0 | 0 | 0 |

CAMERA: tier-stratified benchmark (StM14) {.table .table .table-striped
.table-hover .table-condensed style="width: auto !important; "}

### Status over all features

``` r

show_table(eval_camera$status_counts,
           caption = "CAMERA: reference-relative status of all features")
```

| ref_status                 | n    | pct  |
|----------------------------|------|------|
| Ref-unassigned (no-anchor) | 0    | 0.0  |
| Ref-unassigned             | 0    | 0.0  |
| Ref-confirmed              | 627  | 24.9 |
| Ref-conflicting            | 59   | 2.3  |
| Not evaluable              | 0    | 0.0  |
| Ref-unmatched, co-grouped  | 118  | 4.7  |
| Ref-unmatched, unrelated   | 1715 | 68.1 |

CAMERA: reference-relative status of all features {.table .table
.table-striped .table-hover .table-condensed
style="width: auto !important; "}

### Comparing tools

``` r

show_table(do.call(compare_grouping, unname(evals)),
           caption = "Status (%) per tier and tool",
           highlight_last = TRUE, scroll = TRUE)
```

| Tier | n_ref | CAMERA_pct_no_anchor | CAMERA_pct_unassigned | CAMERA_pct_confirmed | CAMERA_pct_conflicting | RAMClustR_pct_no_anchor | RAMClustR_pct_unassigned | RAMClustR_pct_confirmed | RAMClustR_pct_conflicting |
|----|----|----|----|----|----|----|----|----|----|
| Tier1a | 129 | 0 | 0 | 95.3 | 4.7 | 0 | 0.8 | 63.6 | 35.7 |
| Tier1a_exp | 9 | 0 | 0 | 66.7 | 33.3 | 0 | 11.1 | 44.4 | 44.4 |
| Tier1b | 169 | 0 | 0 | 92.3 | 7.7 | 0 | 0.0 | 73.4 | 26.6 |
| Tier1b_exp | 113 | 0 | 0 | 95.6 | 4.4 | 0 | 0.9 | 69.0 | 30.1 |
| Tier1c | 61 | 0 | 0 | 98.4 | 1.6 | 0 | 1.6 | 86.9 | 11.5 |
| Tier1c_exp | 5 | 0 | 0 | 80.0 | 20.0 | 0 | 20.0 | 40.0 | 40.0 |
| Tier1d | 49 | 0 | 0 | 91.8 | 8.2 | 0 | 0.0 | 65.3 | 34.7 |
| Tier2 | 151 | 0 | 0 | 82.8 | 17.2 | 0 | 0.0 | 52.3 | 47.7 |
| Total | 686 | 0 | 0 | 91.4 | 8.6 | 0 | 0.7 | 66.2 | 33.1 |

Status (%) per tier and tool {.table .table .table-striped .table-hover
.table-condensed style="width: auto !important; "}

### Single compounds

Use `full_table` to see how a tool grouped the features of one compound:

``` r

eval_camera$full_table %>%
  filter(Compound_name == Compound_name[1]) %>%
  select(feature_name, Compound_name, tier, group, anchor_group, ref_status)
```

    ##    feature_name         Compound_name       tier group anchor_group
    ## 1        FT0007 2-Methoxybenzoic acid     Tier1c     1            1
    ## 2        FT0043 2-Methoxybenzoic acid     Tier1c     1            1
    ## 3        FT0095 2-Methoxybenzoic acid     Tier1c     1            1
    ## 4        FT0110 2-Methoxybenzoic acid Tier1a_exp     1            1
    ## 5        FT0174 2-Methoxybenzoic acid     Tier1c     1            1
    ## 6        FT0281 2-Methoxybenzoic acid     Tier1c     1            1
    ## 7        FT0284 2-Methoxybenzoic acid     Tier1a     1            1
    ## 8        FT0292 2-Methoxybenzoic acid      Tier2     1            1
    ## 9        FT0294 2-Methoxybenzoic acid     Tier1c     1            1
    ## 10       FT0295 2-Methoxybenzoic acid     Tier1c     1            1
    ## 11       FT0296 2-Methoxybenzoic acid     Tier1c     1            1
    ## 12       FT0297 2-Methoxybenzoic acid      Tier2     1            1
    ## 13       FT0298 2-Methoxybenzoic acid     Tier1c     1            1
    ## 14       FT0301 2-Methoxybenzoic acid     Tier1b     1            1
    ## 15       FT0315 2-Methoxybenzoic acid     Tier1b     1            1
    ## 16       FT0321 2-Methoxybenzoic acid Tier1c_exp     1            1
    ## 17       FT0322 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 18       FT0451 2-Methoxybenzoic acid     Tier1a     1            1
    ## 19       FT0461 2-Methoxybenzoic acid     Tier1b     1            1
    ## 20       FT0465 2-Methoxybenzoic acid     Tier1b     1            1
    ## 21       FT0565 2-Methoxybenzoic acid     Tier1a     1            1
    ## 22       FT0571 2-Methoxybenzoic acid     Tier1a     1            1
    ## 23       FT0582 2-Methoxybenzoic acid     Tier1a     1            1
    ## 24       FT0591 2-Methoxybenzoic acid     Tier1b     1            1
    ## 25       FT0596 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 26       FT0678 2-Methoxybenzoic acid     Tier1a     1            1
    ## 27       FT0706 2-Methoxybenzoic acid     Tier1a     1            1
    ## 28       FT0853 2-Methoxybenzoic acid Tier1a_exp     1            1
    ## 29       FT0897 2-Methoxybenzoic acid     Tier1a     1            1
    ## 30       FT0900 2-Methoxybenzoic acid     Tier1b     1            1
    ## 31       FT0905 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 32       FT0961 2-Methoxybenzoic acid     Tier1a     1            1
    ## 33       FT0963 2-Methoxybenzoic acid      Tier2     1            1
    ## 34       FT0967 2-Methoxybenzoic acid     Tier1b     1            1
    ## 35       FT1000 2-Methoxybenzoic acid      Tier2     1            1
    ## 36       FT1003 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 37       FT1006 2-Methoxybenzoic acid      Tier2     1            1
    ## 38       FT1011 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 39       FT1442 2-Methoxybenzoic acid     Tier1a     1            1
    ## 40       FT1453 2-Methoxybenzoic acid     Tier1b     1            1
    ## 41       FT1564 2-Methoxybenzoic acid     Tier1a     1            1
    ## 42       FT1578 2-Methoxybenzoic acid     Tier1b     1            1
    ## 43       FT1581 2-Methoxybenzoic acid     Tier1b     1            1
    ## 44       FT1588 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 45       FT1624 2-Methoxybenzoic acid     Tier1a     1            1
    ## 46       FT1627 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 47       FT1635 2-Methoxybenzoic acid     Tier1b     1            1
    ## 48       FT1640 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 49       FT1642 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 50       FT1658 2-Methoxybenzoic acid     Tier1a     1            1
    ## 51       FT1671 2-Methoxybenzoic acid     Tier1b     1            1
    ## 52       FT1768 2-Methoxybenzoic acid     Tier1a     1            1
    ## 53       FT1771 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 54       FT1803 2-Methoxybenzoic acid     Tier1a     1            1
    ## 55       FT1809 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 56       FT1852 2-Methoxybenzoic acid      Tier2     1            1
    ## 57       FT1857 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 58       FT1863 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 59       FT2249 2-Methoxybenzoic acid     Tier1a     1            1
    ## 60       FT2252 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 61       FT2259 2-Methoxybenzoic acid      Tier2     1            1
    ## 62       FT2264 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 63       FT2267 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 64       FT2308 2-Methoxybenzoic acid     Tier1a     1            1
    ## 65       FT2309 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 66       FT2335 2-Methoxybenzoic acid Tier1a_exp     1            1
    ## 67       FT2371 2-Methoxybenzoic acid      Tier2     1            1
    ## 68       FT2375 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 69       FT2377 2-Methoxybenzoic acid Tier1b_exp     1            1
    ## 70       FT2379 2-Methoxybenzoic acid Tier1b_exp     1            1
    ##       ref_status
    ## 1  Ref-confirmed
    ## 2  Ref-confirmed
    ## 3  Ref-confirmed
    ## 4  Ref-confirmed
    ## 5  Ref-confirmed
    ## 6  Ref-confirmed
    ## 7  Ref-confirmed
    ## 8  Ref-confirmed
    ## 9  Ref-confirmed
    ## 10 Ref-confirmed
    ## 11 Ref-confirmed
    ## 12 Ref-confirmed
    ## 13 Ref-confirmed
    ## 14 Ref-confirmed
    ## 15 Ref-confirmed
    ## 16 Ref-confirmed
    ## 17 Ref-confirmed
    ## 18 Ref-confirmed
    ## 19 Ref-confirmed
    ## 20 Ref-confirmed
    ## 21 Ref-confirmed
    ## 22 Ref-confirmed
    ## 23 Ref-confirmed
    ## 24 Ref-confirmed
    ## 25 Ref-confirmed
    ## 26 Ref-confirmed
    ## 27 Ref-confirmed
    ## 28 Ref-confirmed
    ## 29 Ref-confirmed
    ## 30 Ref-confirmed
    ## 31 Ref-confirmed
    ## 32 Ref-confirmed
    ## 33 Ref-confirmed
    ## 34 Ref-confirmed
    ## 35 Ref-confirmed
    ## 36 Ref-confirmed
    ## 37 Ref-confirmed
    ## 38 Ref-confirmed
    ## 39 Ref-confirmed
    ## 40 Ref-confirmed
    ## 41 Ref-confirmed
    ## 42 Ref-confirmed
    ## 43 Ref-confirmed
    ## 44 Ref-confirmed
    ## 45 Ref-confirmed
    ## 46 Ref-confirmed
    ## 47 Ref-confirmed
    ## 48 Ref-confirmed
    ## 49 Ref-confirmed
    ## 50 Ref-confirmed
    ## 51 Ref-confirmed
    ## 52 Ref-confirmed
    ## 53 Ref-confirmed
    ## 54 Ref-confirmed
    ## 55 Ref-confirmed
    ## 56 Ref-confirmed
    ## 57 Ref-confirmed
    ## 58 Ref-confirmed
    ## 59 Ref-confirmed
    ## 60 Ref-confirmed
    ## 61 Ref-confirmed
    ## 62 Ref-confirmed
    ## 63 Ref-confirmed
    ## 64 Ref-confirmed
    ## 65 Ref-confirmed
    ## 66 Ref-confirmed
    ## 67 Ref-confirmed
    ## 68 Ref-confirmed
    ## 69 Ref-confirmed
    ## 70 Ref-confirmed

### Per-feature overview of all standards

One row per reference feature × compound (here: all features MS1Match
assigned to the 14 StM14 standards), with:

- the MS1Match evidence: tier and annotation;
- per tool: the feature’s group, the anchor’s group and the reference
  status (`NA` group = the tool did not group the feature);
- for CAMERA also its own isotope and adduct labels.

The statuses are taken from `evals` (Step 4), so they match the
benchmark.

``` r

# MS1Match evidence per reference row
evidence <- anno_result$export_feature_compound_summary %>%
  select(feature_name, Compound_name, Lib_group, reference_tier,
         any_of(c("grp_anno_match", "exp_anno_match", "iso_annotation"))) %>%
  mutate(MS1Match_anno = coalesce(na_if(grp_anno_match, ""),
                                  na_if(exp_anno_match, ""))) %>%
  select(-any_of(c("grp_anno_match", "exp_anno_match")))

# group, anchor group and status per tool, prefixed with the tool name
tool_cols <- function(tool) {
  evals[[tool]]$full_table %>%
    filter(!is.na(Compound_name)) %>%
    select(feature_name, Compound_name, Lib_group,
           !!paste0(tool, "_group")        := group,
           !!paste0(tool, "_anchor_group") := anchor_group,
           !!paste0(tool, "_ref_status")   := ref_status)
}

feature_overview <- Reduce(
  function(acc, tool) left_join(acc, tool_cols(tool),
                                by = c("feature_name", "Compound_name", "Lib_group")),
  names(evals),
  evidence %>% filter(feature_name %in% evals[[1]]$full_table$feature_name)
)

# CAMERA's own labels
camera_labels <- CAMERA_peaklist %>%
  select(feature_name, any_of(c("isotopes", "adduct"))) %>%
  rename_with(~ paste0("CAMERA_", .x), -feature_name) %>%
  mutate(across(-feature_name, ~ na_if(trimws(as.character(.x)), "")))

feature_overview <- feature_overview %>%
  left_join(feature_table %>% select(feature_name, mzmed, rtmed), by = "feature_name") %>%
  left_join(camera_labels, by = "feature_name") %>%
  relocate(Compound_name, feature_name, mzmed, rtmed, reference_tier,
           MS1Match_anno, any_of("iso_annotation")) %>%
  arrange(Compound_name, rtmed, mzmed)

feature_overview
```

    ## # A tibble: 686 × 16
    ##    Compound_name         feature_name mzmed rtmed reference_tier   MS1Match_anno
    ##    <chr>                 <chr>        <dbl> <dbl> <chr>            <chr>        
    ##  1 2-Methoxybenzoic acid FT2309        496.  507. Tier1b_exp_Isot… NA           
    ##  2 2-Methoxybenzoic acid FT0967        249.  507. Tier1b_Isotope   NA           
    ##  3 2-Methoxybenzoic acid FT1863        361.  507. Tier1b_exp_Isot… NA           
    ##  4 2-Methoxybenzoic acid FT1003        257.  507. Tier1b_exp_Isot… NA           
    ##  5 2-Methoxybenzoic acid FT1011        258.  507. Tier1b_exp_Isot… NA           
    ##  6 2-Methoxybenzoic acid FT0322        138.  507. Tier1b_exp_Isot… NA           
    ##  7 2-Methoxybenzoic acid FT0297        135.  507. Tier2_Unexplain… NA           
    ##  8 2-Methoxybenzoic acid FT0451        153.  507. Tier1a_AnnoMatch [M+H]+ 152.0…
    ##  9 2-Methoxybenzoic acid FT0461        154.  507. Tier1b_Isotope   NA           
    ## 10 2-Methoxybenzoic acid FT0963        249.  507. Tier2_Unexplain… NA           
    ## # ℹ 676 more rows
    ## # ℹ 10 more variables: iso_annotation <chr>, Lib_group <chr>,
    ## #   CAMERA_group <chr>, CAMERA_anchor_group <chr>, CAMERA_ref_status <chr>,
    ## #   RAMClustR_group <chr>, RAMClustR_anchor_group <chr>,
    ## #   RAMClustR_ref_status <chr>, CAMERA_isotopes <chr>, CAMERA_adduct <chr>

Per standard: how many of its reference features each tool grouped with
the anchor, split from it, or left ungrouped.

``` r

status_summary <- function(tool) {
  st <- feature_overview[[paste0(tool, "_ref_status")]]
  feature_overview %>%
    mutate(.st = st) %>%
    group_by(Compound_name) %>%
    summarise(
      !!paste0(tool, "_confirmed")   := sum(.st == "Ref-confirmed",   na.rm = TRUE),
      !!paste0(tool, "_conflicting") := sum(.st == "Ref-conflicting", na.rm = TRUE),
      !!paste0(tool, "_unassigned")  := sum(.st %in% c("Ref-unassigned",
                                                        "Ref-unassigned (no-anchor)")),
      .groups = "drop")
}

per_standard <- Reduce(
  function(acc, tool) left_join(acc, status_summary(tool), by = "Compound_name"),
  names(evals),
  feature_overview %>% count(Compound_name, name = "n_features")
) %>%
  arrange(desc(n_features))

show_table(per_standard,
           caption = "Reference features per standard and their status per tool",
           scroll = TRUE)
```

| Compound_name | n_features | CAMERA_confirmed | CAMERA_conflicting | CAMERA_unassigned | RAMClustR_confirmed | RAMClustR_conflicting | RAMClustR_unassigned |
|----|----|----|----|----|----|----|----|
| Naproxen | 113 | 107 | 6 | 0 | 102 | 10 | 1 |
| Indole-3-acetyl-L-valine | 86 | 82 | 4 | 0 | 56 | 30 | 0 |
| Trans-Ferulic acid | 75 | 71 | 4 | 0 | 14 | 61 | 0 |
| 2-Methoxybenzoic acid | 70 | 70 | 0 | 0 | 40 | 29 | 1 |
| p-Coumaric acid | 59 | 56 | 3 | 0 | 57 | 0 | 2 |
| Phlorizin | 56 | 38 | 18 | 0 | 36 | 20 | 0 |
| Indole-3-carboxaldehyde | 45 | 45 | 0 | 0 | 43 | 2 | 0 |
| Nortriptylin | 45 | 45 | 0 | 0 | 40 | 5 | 0 |
| Rutin | 44 | 21 | 23 | 0 | 14 | 29 | 1 |
| Novobiocin | 39 | 38 | 1 | 0 | 21 | 18 | 0 |
| 3-Indoleacetonitrile | 24 | 24 | 0 | 0 | 10 | 14 | 0 |
| Biochanin A | 18 | 18 | 0 | 0 | 9 | 9 | 0 |
| Trimetoprim | 7 | 7 | 0 | 0 | 7 | 0 | 0 |
| Kinetin | 5 | 5 | 0 | 0 | 5 | 0 | 0 |

Reference features per standard and their status per tool {.table .table
.table-striped .table-hover .table-condensed
style="width: auto !important; "}

``` r

write.csv(feature_overview, "StM14_feature_overview.csv", row.names = FALSE, na = "")
write.csv(per_standard,     "StM14_per_standard.csv",     row.names = FALSE)
```

## Reading the output

For each compound, the Tier 1a `[M+H]+` feature is the **anchor**; the
group the tool assigns to it is the compound’s expected group. Every
reference compound has an anchor by construction. Each reference feature
gets the first status that applies, in this order:

| Status | Meaning |
|:---|:---|
| Ref-unassigned (no-anchor) | The tool did not group the compound’s anchor; none of its features can be evaluated |
| Ref-unassigned | The anchor is grouped, this feature is not |
| Ref-confirmed | Same group as its compound’s anchor |
| Ref-conflicting | Different group than its compound’s anchor |
| Ref-unmatched, co-grouped | Feature outside the reference set, in a group with a compound anchor (possible over-grouping) |
| Ref-unmatched, unrelated | Feature outside the reference set, in a group without an anchor, or in no group |

The first four statuses apply to reference features and make up the tier
benchmark. The last two apply to features outside the reference set and
appear only in `status_counts` and `full_table`.

The benchmark also reports **coverage** (`pct_grouped`): the share of
reference features the tool grouped at all. Tools that group every
feature (e.g. CAMERA) have 100 % coverage and always 0 Ref-unassigned
and 0 Ref-unassigned (no-anchor).

A further status, *Not evaluable*, appears only if a compound’s anchor
is missing from the shared feature table (e.g. cut off by an RT window).
It should be 0;
[`evaluate_grouping()`](https://ruibings.github.io/MS1Match_Rpackage/reference/evaluate_grouping.md)
warns otherwise.

## Using your own data

The same calls, collected. Replace the three inputs:

``` r

# Inputs
feature_table <- my_xcms_result$ft                 # the table MS1Match ran on
anno_result   <- readRDS("my_anno_result.rds")     # from vignette 2a / 2c
grouping      <- camera_to_grouping(getPeaklist(my_xsa))
                 # or ramclustr_to_grouping(RC), mzmine_to_grouping(...),
                 # or prepare_grouping_table(...)

# Check, evaluate
check_feature_match(anno_result, grouping, feature_table)   # must YES

res <- evaluate_grouping(anno_result, grouping,
                         tool_name     = "MyTool",
                         feature_table = feature_table,
                         rt_window     = c(60, 1200))
res

# Export
write.csv(res$benchmark,     "MyTool_benchmark.csv",     row.names = FALSE)
write.csv(res$status_counts, "MyTool_status_counts.csv", row.names = FALSE)
write.csv(res$full_table,    "MyTool_full_table.csv",    row.names = FALSE)
```

To export all tools from Step 4 at once:

``` r

for (tool in names(evals)) {
  write.csv(evals[[tool]]$benchmark,
            paste0(tool, "_benchmark.csv"),  row.names = FALSE)
  write.csv(evals[[tool]]$full_table,
            paste0(tool, "_full_table.csv"), row.names = FALSE)
}
```

## Optional: CAMERA adduct annotation accuracy

Besides groups, CAMERA assigns adduct labels. For Tier 1a features,
their coverage and accuracy are computed with CAMERA-specific functions.
These join on `feature_name` only and do not repeat the feature-table
check. Run them only after Step 3 has said YES.

``` r

eval_full     <- build_eval_full(anno_result$export_feature_compound_summary,
                                 CAMERA_ft = CAMERA_peaklist)
anchors       <- build_anchors(eval_full,
                               full_anno = anno_result$export_reduced_anno_table)
anno_accuracy <- build_anno_accuracy(eval_full, anchors, ppm_tol = 10)
anno_accuracy$summary
```

The full CAMERA analysis (pcgroup purity, adduct confusion, case
studies) is in `analysis/StM14_CAMERA_evaluation.Rmd` in the GitHub
repository.

------------------------------------------------------------------------

## Session info

``` r

sessionInfo()
```

    ## R version 4.6.1 (2026-06-24 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26100)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.2.1         MS1Match_0.0.0.9000
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] DBI_1.3.0                   rlang_1.3.0                
    ##   [3] magrittr_2.0.5              clue_0.3-68                
    ##   [5] MassSpecWavelet_1.78.2      otel_0.2.0                 
    ##   [7] matrixStats_1.5.0           compiler_4.6.1             
    ##   [9] PTMods_1.0.0                systemfonts_1.3.2          
    ##  [11] vctrs_0.7.3                 reshape2_1.4.5             
    ##  [13] stringr_1.6.0               ProtGenerics_1.44.0        
    ##  [15] crayon_1.5.3                pkgconfig_2.0.3            
    ##  [17] MetaboCoreUtils_1.20.1      fastmap_1.2.0              
    ##  [19] XVector_0.52.0              utf8_1.2.6                 
    ##  [21] rmarkdown_2.32              preprocessCore_1.74.0      
    ##  [23] xcms_4.10.1                 ragg_1.5.2                 
    ##  [25] purrr_1.2.2                 xfun_0.62                  
    ##  [27] MultiAssayExperiment_1.38.0 cachem_1.1.0               
    ##  [29] jsonlite_2.0.0              progress_1.2.3             
    ##  [31] DelayedArray_0.38.2         BiocParallel_1.46.0        
    ##  [33] prettyunits_1.2.0           parallel_4.6.1             
    ##  [35] cluster_2.1.8.3             R6_2.6.1                   
    ##  [37] bslib_0.12.0                stringi_1.8.9              
    ##  [39] RColorBrewer_1.1-3          limma_3.68.5               
    ##  [41] GenomicRanges_1.64.0        jquerylib_0.1.4            
    ##  [43] iterators_1.0.14            Rcpp_1.1.2                 
    ##  [45] Seqinfo_1.2.0               SummarizedExperiment_1.42.0
    ##  [47] knitr_1.52                  IRanges_2.46.0             
    ##  [49] Matrix_1.7-6                igraph_2.3.4               
    ##  [51] tidyselect_1.2.1            rstudioapi_0.19.0          
    ##  [53] abind_1.4-8                 yaml_2.3.12                
    ##  [55] doParallel_1.0.17           codetools_0.2-20           
    ##  [57] affy_1.90.0                 lattice_0.23-1             
    ##  [59] tibble_3.3.1                plyr_1.8.9                 
    ##  [61] withr_3.0.3                 Biobase_2.72.0             
    ##  [63] S7_0.2.2                    evaluate_1.0.5             
    ##  [65] desc_1.4.3                  Spectra_1.22.2             
    ##  [67] xml2_1.5.2                  pillar_1.11.1              
    ##  [69] affyio_1.82.0               BiocManager_1.30.27        
    ##  [71] MatrixGenerics_1.24.0       foreach_1.5.2              
    ##  [73] stats4_4.6.1                MSnbase_2.38.1             
    ##  [75] MALDIquant_1.22.3           ncdf4_1.24                 
    ##  [77] generics_0.1.4              hms_1.1.4                  
    ##  [79] S4Vectors_0.50.3            ggplot2_4.0.3              
    ##  [81] scales_1.4.0                MsExperiment_1.14.0        
    ##  [83] glue_1.8.1                  MsFeatures_1.20.0          
    ##  [85] lazyeval_0.2.3              tools_4.6.1                
    ##  [87] mzID_1.50.0                 data.table_1.18.6.1        
    ##  [89] QFeatures_1.22.0            vsn_3.80.0                 
    ##  [91] mzR_2.46.0                  enviPat_2.8                
    ##  [93] fs_2.1.0                    XML_3.99-0.25              
    ##  [95] grid_4.6.1                  impute_1.86.0              
    ##  [97] tidyr_1.3.2                 MsCoreUtils_1.24.0         
    ##  [99] patchwork_1.3.2             PSMatch_1.16.0             
    ## [101] cli_3.6.6                   kableExtra_1.4.1           
    ## [103] textshaping_1.0.5           viridisLite_0.4.3          
    ## [105] S4Arrays_1.12.1             svglite_2.2.2              
    ## [107] AnnotationFilter_1.36.0     pcaMethods_2.4.0           
    ## [109] gtable_0.3.6                sass_0.4.10                
    ## [111] digest_0.6.39               progressr_1.1.0            
    ## [113] BiocGenerics_0.58.1         SparseArray_1.12.3         
    ## [115] ggrepel_0.9.8               htmlwidgets_1.6.4          
    ## [117] farver_2.1.2                htmltools_0.5.9            
    ## [119] pkgdown_2.2.1               lifecycle_1.0.5            
    ## [121] statmod_1.5.2               MASS_7.3-66
