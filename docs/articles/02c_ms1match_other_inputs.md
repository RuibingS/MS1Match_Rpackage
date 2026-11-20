# 2c. MS1Match with other inputs: MZmine, CSV, no MS2

## Overview

Vignettes 2a and 2b start from an xcms feature table. This vignette runs
MS1Match on a feature table from other software instead:

1.  **An MZmine feature table**, imported with
    [`read_mzmine_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md).
    This is the main example; it runs on data shipped with the package.
2.  **Any other CSV feature table**, imported with
    [`read_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_feature_table.md)
    (template only).
3.  **Data without MS2.** The example has MS1 files only, so the
    in-source fragment step is skipped.

The MS1Match steps themselves are the same as in vignette 2a; only the
import and the missing MS2 step differ. The resulting `anno_result` is
the input for vignette *3. Evaluate a feature-grouping tool against
MS1Match*.

``` r

library(MS1Match)
library(BiocParallel)
library(dplyr)
```

## 1. Import the MZmine feature table

The example is StM14 processed with MZmine: 9 MS1 files (4, 8 and 10 µM,
three replicates each), exported as MZmine’s **full feature table** CSV.

The full export is needed, not a filtered one (e.g. “only rows with an
ion identity”). MS1Match should see every feature, and vignette 3 can
only evaluate a grouping on features that are in the table.

``` r

mzmine_file <- system.file("extdata", "MZmine_feature_table", "StM14_MZmine_export",
                           "StM14_mzmine_full_feature_table.csv",
                           package = "MS1Match")

ft <- read_mzmine_feature_table(mzmine_file)

dim(ft)
```

    ## [1] 1172   14

``` r

head(ft[, 1:6])
```

    ## # A tibble: 6 × 6
    ##   feature_name mzmed rtmed intensity_mean intensity_median
    ##   <chr>        <dbl> <dbl>          <dbl>            <dbl>
    ## 1 FT0001       165.   62.2           79.3             78.4
    ## 2 FT0002       151.   62.2          110.             105. 
    ## 3 FT0003       128.   62.3           89.8             89.0
    ## 4 FT0004        85.0  62.3          231.             222. 
    ## 5 FT0005       126.   62.3          309.             288. 
    ## 6 FT0006       111.   62.4          137.             131. 
    ## # ℹ 1 more variable: `Pool-4umSi11_neu MS_GA3_01_50026` <dbl>

[`read_mzmine_feature_table()`](https://ruibings.github.io/MS1Match_Rpackage/reference/read_mzmine_feature_table.md)
converts the MZmine export to the layout the other MS1Match functions
expect:

- MZmine `id` becomes `feature_name` (`1` → `FT0001`);
- `mz` becomes `mzmed`, and `rt` is converted from minutes to seconds as
  `rtmed`;
- the per-file `datafile:<file>:area` columns become the sample
  intensities, plus `intensity_mean` and `intensity_median`;
- features found in fewer than `min_detect_fraction` of the samples are
  dropped (default 0.5, i.e. at least 5 of the 9 files here).

MZmine’s own groupings (`feature_group`, `ion_identities:iin_id`) are
not part of `ft`. Section 5 reads them from the same CSV.

The MZmine features in this example all elute between 62 and 1197 s, so
no extra RT window is needed. For other data, apply the same 60–1200 s
window as in vignette 2a before running MS1Match.

## 2. Library and lookup tables

The same library and tables as in vignette 2a:

``` r

adduct_table <- read.csv(
  system.file("extdata", "Adducts", "adducts_with_4M.csv", package = "MS1Match")
)

neutral_loss_table <- read.csv(
  system.file("extdata", "Neutral loss", "neutral loss database.csv", package = "MS1Match"),
  check.names = TRUE
)

Si14_library <- import_msp(
  file_path     = system.file("extdata", "StM14_library", "Si14_library.msp",
                              package = "MS1Match"),
  ionPolarity   = "pos",
  spectrum_type = "MS1"
)
```

    ## [parser] Parsed:          14
    ## [parser] Filtered:        0
    ## [parser] Pre-dedup total: 14
    ## [parser] Duplicates:      0
    ## [parser] Building result: 14 compounds...
    ## [parser] Done.

## 3. MS1Match pipeline without MS2

Library matching and the Phase 1 / Phase 2 annotation run as in vignette
2a. Without MS2 spectra, the in-source fragment step
([`createPeakAnnotationISF()`](https://ruibings.github.io/MS1Match_Rpackage/reference/createPeakAnnotationISF.md)
and `annotateISF()`) is skipped and
[`join_annotation_results()`](https://ruibings.github.io/MS1Match_Rpackage/reference/join_annotation_results.md)
gets `isf_result = NULL`.

``` r

direct_obj <- createDirectLibMatch(
  FT               = ft,
  compounds_list   = Si14_library,
  rt_thr_compound  = 30,
  rt_thr_coelution = 5,
  ppm              = 10,
  ms_level         = "MS1",
  polarity         = "pos"
)

direct_obj <- libAnnotate(direct_obj, compute_cosine = TRUE, ppm = 20,
                          BPPARAM = SerialParam())
direct_res <- results(direct_obj)

anno_obj <- createPeakAnnotationFromLibraryMatch(
  library_matches    = direct_res,
  adducts_df         = adduct_table,
  neutral_loss_table = neutral_loss_table,
  FT                 = ft,
  polarity           = "positive",
  ppm                = 30
)
lib_anno_res <- annotateLibMatch(anno_obj, cutoff_MH = TRUE)

anno_result <- join_annotation_results(
  FT              = ft,
  direct_lib_list = direct_res,
  lib_anno_list   = lib_anno_res,
  isf_result      = NULL,         # no MS2 data
  include_phase2  = TRUE
)
```

Without MS2, the tiers Tier 1c and Tier 1c_exp (in-source MS2 fragments)
are not assigned. Fragments that would otherwise be Tier 1c usually end
up in Tier 2 (unexplained, below `[M+H]+`). All other tiers work as
usual.

``` r

anno_result$export_feature_compound_summary %>%
  count(reference_tier, annotation_phase) %>%
  arrange(reference_tier)
```

    ## # A tibble: 6 × 3
    ##   reference_tier       annotation_phase     n
    ##   <chr>                <chr>            <int>
    ## 1 Tier1a_AnnoMatch     Phase1             110
    ## 2 Tier1a_exp_AnnoMatch Phase2               8
    ## 3 Tier1b_Isotope       Phase1              15
    ## 4 Tier1b_exp_Isotope   Phase2              11
    ## 5 Tier1d_NeutralLoss   Phase1              53
    ## 6 Tier2_Unexplained    Phase1             125

## 4. Save the result

``` r

saveRDS(ft,          "ft_StM14_mzmine.rds")
saveRDS(anno_result, "anno_result_StM14_mzmine.rds")
```

The tables of `anno_result` are exported exactly as in vignette 2a,
Section 5.

## 5. Next: evaluate MZmine’s own grouping

The same CSV also holds MZmine’s groupings. MS1Match ran on the full
feature table (Section 1); the groupings are read from that file and
restricted to the same features with `keep_features`. Vignette 3
explains the evaluation; the short version:

``` r

mzmine_raw <- read.csv(mzmine_file, check.names = FALSE)

groupings <- list(
  # correlation grouping (the counterpart of CAMERA pcgroups)
  MZmine_corr = mzmine_to_grouping(mzmine_raw, group_col = "feature_group",
                                   keep_features = ft$feature_name,
                                   unassigned = ""),
  # ion identity networks
  MZmine_IIN  = mzmine_to_grouping(mzmine_raw, group_col = "ion_identities:iin_id",
                                   keep_features = ft$feature_name,
                                   unassigned = "")
)

evals <- lapply(names(groupings), function(tool)
  evaluate_grouping(anno_result, groupings[[tool]], tool_name = tool,
                    feature_table = ft, rt_window = c(60, 1200)))
names(evals) <- names(groupings)

do.call(compare_grouping, unname(evals))
```

    ##         Tier n_ref MZmine_corr_pct_no_anchor MZmine_corr_pct_unassigned
    ## 1     Tier1a   110                       0.9                        2.7
    ## 2 Tier1a_exp     8                       0.0                       25.0
    ## 3     Tier1b    15                       0.0                       20.0
    ## 4 Tier1b_exp    11                       0.0                       18.2
    ## 5     Tier1d    53                       0.0                        1.9
    ## 6      Tier2   125                       0.0                        8.0
    ## 7      Total   322                       0.3                        6.5
    ##   MZmine_corr_pct_confirmed MZmine_corr_pct_conflicting
    ## 1                      92.7                         3.6
    ## 2                      12.5                        62.5
    ## 3                      66.7                        13.3
    ## 4                      54.5                        27.3
    ## 5                      98.1                         0.0
    ## 6                      84.0                         8.0
    ## 7                      85.7                         7.5
    ##   MZmine_IIN_pct_no_anchor MZmine_IIN_pct_unassigned MZmine_IIN_pct_confirmed
    ## 1                     10.9                      47.3                     38.2
    ## 2                     12.5                      87.5                      0.0
    ## 3                     40.0                      60.0                      0.0
    ## 4                     27.3                      63.6                      0.0
    ## 5                     34.0                      62.3                      0.0
    ## 6                     16.8                      78.4                      0.0
    ## 7                     18.9                      64.0                     13.0
    ##   MZmine_IIN_pct_conflicting
    ## 1                        3.6
    ## 2                        0.0
    ## 3                        0.0
    ## 4                        9.1
    ## 5                        3.8
    ## 6                        4.8
    ## 7                        4.0

Two things to keep in mind:

- `unassigned = ""` matters here. MZmine numbers its groups from 0, so
  group `0` is a real group, and only an empty cell means “not grouped”.
- Ion identity networks cover only the features MZmine could give an ion
  identity, so many reference features have no network. If a compound’s
  `[M+H]+` anchor has no network, all its features are Ref-unassigned
  (no-anchor); if the anchor has one, its features without a network are
  Ref-unassigned. Read both as IIN’s coverage (`pct_grouped` in the
  benchmark), not as errors. Ref-confirmed vs Ref-conflicting shows
  whether the networks it did form are right. `feature_group`
  (correlation grouping) covers most features and is the direct
  counterpart of CAMERA’s pcgroups.

### Per-feature overview

One row per feature that MS1Match assigned to a library compound
(feature × compound), with:

- the MS1Match evidence: tier and annotation;
- MZmine’s correlation group and ion identity network, plus the ion
  identity MZmine assigned (`NA` = MZmine did not group the feature);
- the anchor’s group and the reference status of the feature for each
  MZmine grouping.

The statuses are taken from `evals`, so they are the same as in the
summary table above.

``` r

# MZmine's own ion identity label per feature
mzmine_ion <- mzmine_raw %>%
  transmute(feature_name  = sprintf("FT%04d", as.integer(id)),
            MZmine_ion_id = na_if(trimws(`ion_identities:ion_identities`), ""))

# status per reference row, one block of columns per MZmine grouping
status_cols <- function(tool, prefix) {
  evals[[tool]]$full_table %>%
    filter(!is.na(Compound_name)) %>%
    select(feature_name, Compound_name, Lib_group,
           !!paste0(prefix, "_group")        := group,
           !!paste0(prefix, "_anchor_group") := anchor_group,
           !!paste0(prefix, "_ref_status")   := ref_status)
}

# MS1Match evidence per reference row
evidence <- anno_result$export_feature_compound_summary %>%
  select(feature_name, Compound_name, Lib_group, reference_tier,
         any_of(c("annotation_phase", "grp_anno_match", "exp_anno_match",
                  "NL_match", "iso_annotation")))

feature_overview <- evidence %>%
  inner_join(status_cols("MZmine_corr", "corr"),
             by = c("feature_name", "Compound_name", "Lib_group")) %>%
  left_join(status_cols("MZmine_IIN", "IIN"),
            by = c("feature_name", "Compound_name", "Lib_group")) %>%
  left_join(ft %>% select(feature_name, mzmed, rtmed), by = "feature_name") %>%
  left_join(mzmine_ion, by = "feature_name") %>%
  mutate(
    MS1Match_anno = coalesce(na_if(grp_anno_match, ""), na_if(exp_anno_match, "")),
    MZmine_grouped = case_when(
      !is.na(corr_group) & !is.na(IIN_group) ~ "corr + IIN",
      !is.na(corr_group)                     ~ "corr only",
      !is.na(IIN_group)                      ~ "IIN only",
      TRUE                                   ~ "not grouped"
    )
  ) %>%
  select(Compound_name, feature_name, mzmed, rtmed,
         reference_tier, MS1Match_anno, any_of("iso_annotation"),
         MZmine_grouped,
         corr_group, corr_anchor_group, corr_ref_status,
         IIN_group, MZmine_ion_id, IIN_anchor_group, IIN_ref_status) %>%
  arrange(Compound_name, rtmed, mzmed)

feature_overview
```

    ## # A tibble: 322 × 15
    ##    Compound_name         feature_name mzmed rtmed reference_tier   MS1Match_anno
    ##    <chr>                 <chr>        <dbl> <dbl> <chr>            <chr>        
    ##  1 2-Methoxybenzoic acid FT0243        495.  505. Tier1a_AnnoMatch [3M+K]+ 152.…
    ##  2 2-Methoxybenzoic acid FT0247        135.  506. Tier1a_AnnoMatch [M+H-H2O]+ 1…
    ##  3 2-Methoxybenzoic acid FT0248        359.  506. Tier2_Unexplain… NA           
    ##  4 2-Methoxybenzoic acid FT0250        172.  506. Tier1a_AnnoMatch [2M+H+K]2+ 1…
    ##  5 2-Methoxybenzoic acid FT0254        191.  506. Tier1a_AnnoMatch [M+K]+ 152.0…
    ##  6 2-Methoxybenzoic acid FT0255        327.  506. Tier1a_AnnoMatch [2M+Na]+ 152…
    ##  7 2-Methoxybenzoic acid FT0256        343.  506. Tier1a_AnnoMatch [2M+K]+ 152.…
    ##  8 2-Methoxybenzoic acid FT0258        135.  506. Tier2_Unexplain… NA           
    ##  9 2-Methoxybenzoic acid FT0261        135.  506. Tier1d_NeutralL… NA           
    ## 10 2-Methoxybenzoic acid FT0260        175.  506. Tier1a_AnnoMatch [M+Na]+ 152.…
    ## # ℹ 312 more rows
    ## # ℹ 9 more variables: iso_annotation <chr>, MZmine_grouped <chr>,
    ## #   corr_group <chr>, corr_anchor_group <chr>, corr_ref_status <chr>,
    ## #   IIN_group <chr>, MZmine_ion_id <chr>, IIN_anchor_group <chr>,
    ## #   IIN_ref_status <chr>

Which reference features did MZmine not group, per compound:

``` r

feature_overview %>%
  group_by(Compound_name) %>%
  summarise(
    n_features        = n(),
    n_corr_grouped    = sum(!is.na(corr_group)),
    n_IIN_grouped     = sum(!is.na(IIN_group)),
    n_not_grouped     = sum(MZmine_grouped == "not grouped"),
    corr_confirmed    = sum(corr_ref_status == "Ref-confirmed"),
    corr_conflicting  = sum(corr_ref_status == "Ref-conflicting"),
    IIN_confirmed     = sum(IIN_ref_status  == "Ref-confirmed"),
    IIN_conflicting   = sum(IIN_ref_status  == "Ref-conflicting"),
    .groups = "drop"
  ) %>%
  arrange(desc(n_features))
```

    ## # A tibble: 14 × 9
    ##    Compound_name           n_features n_corr_grouped n_IIN_grouped n_not_grouped
    ##    <chr>                        <int>          <int>         <int>         <int>
    ##  1 Indole-3-acetyl-L-vali…         59             55             9             4
    ##  2 Naproxen                        50             45            11             5
    ##  3 2-Methoxybenzoic acid           36             36            11             0
    ##  4 Trans-Ferulic acid              31             29             4             2
    ##  5 Phlorizin                       30             24             0             6
    ##  6 Nortriptylin                    28             28             3             0
    ##  7 p-Coumaric acid                 24             24             3             0
    ##  8 Indole-3-carboxaldehyde         20             20             5             0
    ##  9 Rutin                           15             13             3             2
    ## 10 Novobiocin                      11             10             4             1
    ## 11 3-Indoleacetonitrile             9              9             2             0
    ## 12 Biochanin A                      6              5             3             1
    ## 13 Kinetin                          2              2             0             0
    ## 14 Trimetoprim                      1              0             0             1
    ## # ℹ 4 more variables: corr_confirmed <int>, corr_conflicting <int>,
    ## #   IIN_confirmed <int>, IIN_conflicting <int>

``` r

# overall: how the reference features split across the two MZmine groupings
feature_overview %>% count(MZmine_grouped)
```

    ## # A tibble: 3 × 2
    ##   MZmine_grouped     n
    ##   <chr>          <int>
    ## 1 corr + IIN        58
    ## 2 corr only        242
    ## 3 not grouped       22

``` r

table(corr = feature_overview$corr_ref_status,
      IIN  = feature_overview$IIN_ref_status, useNA = "ifany")
```

    ##                             IIN
    ## corr                         Ref-confirmed Ref-conflicting Ref-unassigned
    ##   Ref-confirmed                         42              11            179
    ##   Ref-conflicting                        0               2             12
    ##   Ref-unassigned                         0               0             15
    ##   Ref-unassigned (no-anchor)             0               0              0
    ##                             IIN
    ## corr                         Ref-unassigned (no-anchor)
    ##   Ref-confirmed                                      44
    ##   Ref-conflicting                                    10
    ##   Ref-unassigned                                      6
    ##   Ref-unassigned (no-anchor)                          1

``` r

write.csv(feature_overview, "StM14_mzmine_feature_overview.csv",
          row.names = FALSE, na = "")
```

## 6. Any other CSV feature table

For other software (MS-DIAL, OpenMS, an Excel export, …), name the
columns of your table. All other columns, e.g. a group / cluster column,
are kept.

``` r

ft <- read_feature_table(
  "my_feature_table.csv",
  mz_col      = "mz",
  rt_col      = "rt",
  rt_unit     = "min",          # "sec" or "min"; converted to seconds
  sample_cols = "\\.mzML$",     # sample intensity columns: names or a regex
  id_col      = NULL            # existing id column, or NULL for FT0001, ...
)
```

Then continue with Sections 2 and 3. If your MS2 spectra are available,
run the in-source fragment step as in vignette 2a, Sections 2 and 4.

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
    ## [1] dplyr_1.2.1         BiocParallel_1.46.0 MS1Match_0.0.0.9000
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
    ##  [31] DelayedArray_0.38.2         prettyunits_1.2.0          
    ##  [33] parallel_4.6.1              cluster_2.1.8.3            
    ##  [35] R6_2.6.1                    bslib_0.12.0               
    ##  [37] stringi_1.8.9               RColorBrewer_1.1-3         
    ##  [39] limma_3.68.5                GenomicRanges_1.64.0       
    ##  [41] jquerylib_0.1.4             iterators_1.0.14           
    ##  [43] Rcpp_1.1.2                  Seqinfo_1.2.0              
    ##  [45] SummarizedExperiment_1.42.0 knitr_1.52                 
    ##  [47] IRanges_2.46.0              Matrix_1.7-6               
    ##  [49] igraph_2.3.4                tidyselect_1.2.1           
    ##  [51] rstudioapi_0.19.0           abind_1.4-8                
    ##  [53] yaml_2.3.12                 doParallel_1.0.17          
    ##  [55] codetools_0.2-20            affy_1.90.0                
    ##  [57] lattice_0.23-1              tibble_3.3.1               
    ##  [59] plyr_1.8.9                  withr_3.0.3                
    ##  [61] Biobase_2.72.0              S7_0.2.2                   
    ##  [63] evaluate_1.0.5              desc_1.4.3                 
    ##  [65] Spectra_1.22.2              pillar_1.11.1              
    ##  [67] affyio_1.82.0               BiocManager_1.30.27        
    ##  [69] MatrixGenerics_1.24.0       foreach_1.5.2              
    ##  [71] stats4_4.6.1                MSnbase_2.38.1             
    ##  [73] MALDIquant_1.22.3           ncdf4_1.24                 
    ##  [75] generics_0.1.4              hms_1.1.4                  
    ##  [77] S4Vectors_0.50.3            ggplot2_4.0.3              
    ##  [79] scales_1.4.0                MsExperiment_1.14.0        
    ##  [81] glue_1.8.1                  MsFeatures_1.20.0          
    ##  [83] lazyeval_0.2.3              tools_4.6.1                
    ##  [85] mzID_1.50.0                 data.table_1.18.6.1        
    ##  [87] QFeatures_1.22.0            vsn_3.80.0                 
    ##  [89] mzR_2.46.0                  enviPat_2.8                
    ##  [91] fs_2.1.0                    XML_3.99-0.25              
    ##  [93] grid_4.6.1                  impute_1.86.0              
    ##  [95] tidyr_1.3.2                 MsCoreUtils_1.24.0         
    ##  [97] patchwork_1.3.2             PSMatch_1.16.0             
    ##  [99] cli_3.6.6                   textshaping_1.0.5          
    ## [101] S4Arrays_1.12.1             AnnotationFilter_1.36.0    
    ## [103] pcaMethods_2.4.0            gtable_0.3.6               
    ## [105] sass_0.4.10                 digest_0.6.39              
    ## [107] progressr_1.1.0             BiocGenerics_0.58.1        
    ## [109] SparseArray_1.12.3          ggrepel_0.9.8              
    ## [111] htmlwidgets_1.6.4           farver_2.1.2               
    ## [113] htmltools_0.5.9             pkgdown_2.2.1              
    ## [115] lifecycle_1.0.5             statmod_1.5.2              
    ## [117] MASS_7.3-66
