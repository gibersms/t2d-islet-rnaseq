# =====================================================================
# Project settings
# This is the only file you should need to edit.
# =====================================================================

# ---- Dataset ---------------------------------------------------------
# GEO series: RNA-seq of pancreatic islets from living donors who had
# part of their pancreas removed (Wigger et al., Nature Metabolism 2021).
GSE_ID <- "GSE164416"

# Genome build used in NCBI's count file names. Only change this if the
# GEO download page shows a different name (e.g. GRCh38.p14).
NCBI_BUILD <- "GRCh38.p13"

# ---- Groups to compare ----------------------------------------------
# Run scripts/01_download_data.R first. It prints every sample
# characteristic it finds. Copy the exact column name and values here.
#
# Example (yours may differ!):
#   GROUP_COLUMN  <- "diabetes status:ch1"
#   CONTROL_LABEL <- "ND"
#   CASE_LABEL    <- "T2D"
GROUP_COLUMN  <- "FILL_ME_IN"   # column that holds the diagnosis
CONTROL_LABEL <- "FILL_ME_IN"   # value meaning non-diabetic
CASE_LABEL    <- "FILL_ME_IN"   # value meaning type 2 diabetes

# ---- Covariates (optional) -------------------------------------------
# Other sample characteristics to adjust for in the model, written as
#   short_name = "exact GEO column name"
# Leave as character(0) to run without adjustment.
#
# Example:
#   COVARIATES <- c(sex = "gender:ch1", age = "age:ch1")
COVARIATES <- character(0)

# ---- Analysis thresholds ---------------------------------------------
FDR_CUTOFF <- 0.05   # adjusted p-value cut-off for "significant"
MIN_COUNT  <- 10     # keep a gene if it has at least MIN_COUNT reads in at
                     # least as many samples as the smallest group
