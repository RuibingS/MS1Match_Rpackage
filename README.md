# MS1Match <img src="man/logo/MS1Match_Logo.png" align="right" width=120 height=120 alt="" />

**MS1Match** is an R package for library-guided, compound-centric feature
grouping and annotation of LC-MS data (MS1 and MS2). It offers:

1. Preprocessing of raw spectral data with the xcms framework
2. Direct library-based annotation of detected features
3. Annotation of adducts, isotopes, neutral losses and in-source MS2
   fragments, each with an evidence tier
4. Targeted annotation based on a user-defined compound list
5. Benchmarking of other feature-grouping tools (CAMERA, RAMClustR, MZmine,
   ...) against the tier-stratified MS1Match reference set

# Overview

<p align="center">
  <img src="man/figures/ms1match_overview.svg" width="100%"
       alt="MS1Match workflow: inputs, the four-step MS1Match pipeline producing anno_result, and the optional evaluation of feature-grouping tools" />
</p>

> **Access before publication.** The repository will be made public when the
> manuscript is published. Until then, it is shared with collaborators and
> reviewers by invitation. Installing or cloning it requires a GitHub
> personal access token from an invited account; see
> [Access token](#access-token-before-publication) below.

# Installation

## Requirements

- **R** `>= 4.5`, with the matching Bioconductor release
- **Rcpp**: required (the package contains compiled C++ code)
- **Build tools** to compile the C++ code:
  - **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
    matching your R version (R 4.5.x → **Rtools45**)
  - **macOS:** Apple Command Line Tools

    ```bash
    xcode-select --install
    ```
  - **Linux (Ubuntu/Debian):**

    ```bash
    sudo apt update
    sudo apt install build-essential r-base-dev
    ```

```r
install.packages("Rcpp")
```

## Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(
  c("xcms", "MSnbase", "Spectra", "MetaboCoreUtils",
    "MsExperiment", "MsBackendMgf", "BiocParallel", "S4Vectors"),
  ask = FALSE, update = TRUE
)
```

Optional, only for the tool comparisons in tutorial 3:

```r
BiocManager::install("CAMERA")
install.packages("RAMClustR")
```

## Install MS1Match

```r
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")

remotes::install_github("RuibingS/MS1Match_Rpackage",
                        build_vignettes = TRUE,
                        dependencies    = TRUE)
```

`build_vignettes = TRUE` knits the tutorials during installation, which
takes a few minutes. Leave it out for a faster install; the tutorials are
then not available offline with `vignette()`.

## Check the installation

```r
library(MS1Match)
packageVersion("MS1Match")
browseVignettes("MS1Match")
```

## Troubleshooting

- **`Permission denied` on a `.dll` (Windows):** the package is loaded in an
  R session and its DLL is locked. Restart R (*Session → Restart R* in
  RStudio), close other R sessions, and install again before loading any
  package.
- **"... dependencies are out of date":** a notice, not an error. To skip
  updating them, add `upgrade = "never"` to `install_github()`.
- **A Bioconductor package fails to install:** check that the R and
  Bioconductor versions match with `BiocManager::valid()`.

# Using Docker

A Dockerfile is provided to run the package in a reproducible environment
without installing the dependencies manually.

### 1. Clone the repository

```bash
git clone https://github.com/RuibingS/MS1Match_Rpackage.git
cd MS1Match_Rpackage
```

### 2. Build the Docker image

```bash
docker build -t ms1match .
```

### 3. Run an interactive R session

```bash
docker run -it --rm ms1match
```

- `-it` opens an interactive terminal for the R console
- `--rm` removes the container after you exit R (`q()`)

Inside the container:

```r
library(MS1Match)
```

### Optional: mount a local folder

To work with files from your computer:

```bash
docker run -it --rm -v "$(pwd)":/workspace ms1match
```

This mounts the current directory at `/workspace` inside the container. On
Windows (PowerShell), use `-v "${PWD}:/workspace"`.

# Tutorials

The tutorials follow the three steps of a typical analysis. Each one starts
from the output of the previous one, so you can enter wherever your data
already are.

```
raw files ──> [1] xcms ──> feature table ──> [2] MS1Match ──> anno_result ──> [3] evaluate a grouping tool
              (optional)                     (StM14 / SuperPool)              (CAMERA, RAMClustR, any tool)
```

| # | Tutorial | What it covers | Start here if you have |
|---|----------|----------------|------------------------|
| 1 | [Preprocessing with xcms](https://RuibingS.github.io/MS1Match_Rpackage/articles/01_preprocessing_xcms.html) | Peak picking, alignment and correspondence (`runxcmsFT()` or step by step); export of the feature table | Raw `.mzML` / `.mzXML` files |
| 2a | [MS1Match workflow: StM14](https://RuibingS.github.io/MS1Match_Rpackage/articles/02a_ms1match_StM14.html) | **Main tutorial.** 14 standards spiked into a *P. aeruginosa* extract, MS1 + MS2: consensus MS2 spectra, direct library match, Phase 1 / Phase 2 annotation, in-source fragments, tier summary, export | An xcms feature table |
| 2b | [MS1Match workflow: SuperPool](https://RuibingS.github.io/MS1Match_Rpackage/articles/02b_ms1match_SuperPool.html) | The same workflow on a pooled sample of MSMLS standards and a larger in-house library | An xcms feature table (larger example) |
| 2c | [Other inputs: MZmine, CSV, no MS2](https://RuibingS.github.io/MS1Match_Rpackage/articles/02c_ms1match_other_inputs.html) | MS1Match on an MZmine full feature table (`read_mzmine_feature_table()`) or any CSV (`read_feature_table()`), without MS2; evaluating MZmine's own groupings | A feature table from MZmine or other software |
| 3 | [Evaluate a grouping tool](https://RuibingS.github.io/MS1Match_Rpackage/articles/03_evaluate_grouping_tools.html) | Benchmarking CAMERA, RAMClustR, MZmine or any other tool against the MS1Match reference set: standard grouping table, same-feature-table check, tier-stratified results | An `anno_result` and a grouping result from the same feature table |

After installation, the same tutorials are available offline:

```r
browseVignettes("MS1Match")
vignette("02a_ms1match_StM14", package = "MS1Match")
```

### Which path do I follow?

- **Starting from raw files:** 1 → 2a → (3)
- **Feature table from xcms:** 2a → (3)
- **Feature table from MZmine, MS-DIAL or a CSV:** 2c → (3)
- **Benchmarking a grouping tool and you already have `anno_result`:** 3 only

Tutorial 3 never runs MS1Match itself. It reuses the `anno_result` and the
feature table from tutorial 2, so the grouping tool must have been run on
**the same feature table**. Tutorial 3 shows how to check this.

### Example data

All tutorials run on example data shipped in `inst/extdata`. Slow steps
(peak picking, MS2 consensus building) are loaded from cached `.rds` files;
the code to recompute them from raw files is included but not run.

| Dataset | Content | Library |
|---------|---------|---------|
| StM14 | 14 standards spiked into a *P. aeruginosa* metabolome extract at three concentrations, MS1 + DDA MS2 | Si14 (14 compounds) |
| StM14 (MZmine) | The same samples processed with MZmine 3: full feature table of the 9 MS1 files, with MZmine's correlation groups and ion identity networks | Si14 (14 compounds) |
| SuperPool | Pooled MSMLS standards, MS1 + DDA MS2 | CBIO in-house library, restricted to the pool's targets |

# Analyses for the manuscript

The `analysis/` folder holds the full benchmark analyses reported in the
manuscript: CAMERA and RAMClustR evaluations on StM14 and SuperPool. All
Ref-* counts come from `evaluate_grouping()`, the same function as in
tutorial 3; the CAMERA files add adduct-annotation accuracy, confusion
analysis and case studies. They are not built as vignettes. Render them with
`rmarkdown::render()` after installing the package; the knitted reports are
in `analysis/html/`.

| File | Content |
|------|---------|
| `analysis/StM14_CAMERA_evaluation.Rmd` | CAMERA pcgroups vs. the MS1Match reference set: tier benchmark, pcgroup purity, adduct annotation coverage / accuracy, failure-cause analysis, compound case studies |
| `analysis/StM14_RAMClustR_evaluation.Rmd` | The same grouping evaluation for RAMClustR, with coverage and cluster purity per compound |
| `analysis/SuperPool_CAMERA_evaluation.Rmd` | The CAMERA evaluation on SuperPool |
| `analysis/SuperPool_RAMClustR_evaluation.Rmd` | The RAMClustR evaluation on SuperPool |

# Citation

If you use MS1Match, please cite: *(manuscript in preparation)*.
