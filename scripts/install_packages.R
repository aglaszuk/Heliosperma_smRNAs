#!/usr/bin/env Rscript
# Install R packages required by the Heliosperma smRNA analysis scripts.
# Run once from the project root: Rscript scripts/install_packages.R

cran <- c(
  "RColorBrewer",
  "statmod",
  "gplots",
  "devtools",
  "ggpubr",
  "dplyr",
  "tidyverse",
  "ggplot2",
  "patchwork",
  "xtable",
  "pkgconfig",
  "drake",
  "circlize",
  "gridExtra"
)

bioc <- c(
  "RUVSeq",
  "edgeR",
  "HTSFilter",
  "topGO",
  "ComplexHeatmap"
)

# GOplot is on CRAN
cran <- c(cran, "GOplot")

message("Installing CRAN packages...")
install.packages(setdiff(cran, rownames(installed.packages())), repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

message("Installing Bioconductor packages...")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), ask = FALSE, update = FALSE)

message("Done. Open Heliosperma_smRNAs.Rproj in RStudio, then run scripts in order (see README.md).")
