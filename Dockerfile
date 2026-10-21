# We use rocker/r-ver:devel to get the development version of R, or the latest 
# which corresponds to R >= 4.5 as it becomes available. You can also specify 
# FROM rocker/r-ver:4.5.0 once it is released as the stable version.
FROM rocker/r-ver:latest

# Update system repositories and install build tools (for Rcpp) along with
# standard dependencies necessary to build bioconductor packages such as xcms/MSnbase.
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    build-essential \
    libuv1-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libnetcdf-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpng-dev \
    libglpk-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
 && rm -rf /var/lib/apt/lists/*

# Install the Rcpp package for compiling C++ code
RUN R -e "install.packages('Rcpp', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Install BiocManager and the specified Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "BiocManager::install(c('xcms', 'MSnbase', 'Spectra', 'MetaboCoreUtils', 'BiocParallel', 'IRanges', 'S4Vectors', 'progressr', 'purrr'), ask = FALSE, update = TRUE)"

# Copy the local package source to the container and install it
COPY . /MS1Match_Rpackage
RUN R -e "if (!requireNamespace('remotes', quietly = TRUE)) install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "options(warn = 2); remotes::install_local('/MS1Match_Rpackage', dependencies = c('Depends', 'Imports', 'LinkingTo'), repos = BiocManager::repositories())"

# Verify installations by pre-loading packages during build (will throw error if one fails)
RUN R -e " \
    library(MS1Match); \
    library(xcms); \
    library(Spectra); \
    library(BiocParallel); \
    library(tibble); \
    library(dplyr); \
    library(tidyr); \
    print('All packages successfully loaded!')"

# Setup your default command (e.g. interactive R session)
CMD ["R"]
