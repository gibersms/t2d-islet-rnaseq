# Step 0 - Install the R packages this project needs (run once)
# ---------------------------------------------------------------------
# CRAN packages come from the main R repository; Bioconductor packages
# (DESeq2, GEOquery, apeglm) come from Bioconductor, the repository for
# biology/bioinformatics packages, via BiocManager.
#
# If R asks "Do you want to install from sources the packages which need
# compilation?", answer "no".

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("here", "ggplot2", "ggrepel", "pheatmap")
bioc_packages <- c("DESeq2", "GEOquery", "apeglm")

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) install.packages(missing_cran)

missing_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# Check everything loads, and show the versions (useful for the README)
for (pkg in c(cran_packages, bioc_packages)) {
  ok <- suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))
  cat(sprintf("%-10s %s\n", pkg,
              if (ok) as.character(packageVersion(pkg)) else "NOT INSTALLED"))
}
cat("\n", R.version.string, " | Bioconductor ",
    as.character(BiocManager::version()), "\n", sep = "")
