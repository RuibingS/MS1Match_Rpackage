# MS1Match <img src="man/logo/MS1Match_Logo.png" align="right" width=120 height=120 alt="" />

**MS1Match** is an R package for annotating mass spectrometry (MS1 and MS2) feature tables using compound libraries. It offers:

1. Preprocessing of raw spectral data using the xcms framework

2. Direct library-based annotation of detected features

3. Annotation of adducts, neutral losses, and isotopic patterns

4. Targeted annotation based on a user-defined compound list

```r
# Install required Bioconductor packages first
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("xcms", "MSnbase", "Spectra"))

# Install devtools if not already
install.packages("devtools")

# Install MS1Match from GitHub
devtools::install_github("RuibingS/MS1Match_Rpackage")

```
## 📘 Vignette

A full reproducible tutorial using the demo StM14 data set is available:

👉 [MS1Match Tutorial](vignettes/MS1Match_tutorial.html)


