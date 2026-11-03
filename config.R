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
# Column in the sample information that holds the diagnosis, and the
# values for the two groups being compared.
GROUP_COLUMN  <- "diabetes status:ch1"
CONTROL_LABEL <- "ND"    # non-diabetic
CASE_LABEL    <- "T2D"   # type 2 diabetes

# ---- Optional sample filter -------------------------------------------
# Keep only samples with a given value in another column, e.g. the
# authors' filtered subset:
#   SUBSET_COLUMN <- "in_ins_filtered_data_subset:ch1"
#   SUBSET_VALUE  <- "TRUE"
# Leave both as NA to keep all samples in the two groups.
SUBSET_COLUMN <- NA
SUBSET_VALUE  <- NA

# ---- Covariates (optional) -------------------------------------------
# Other sample characteristics to adjust for in the model, written as
#   short_name = "exact GEO column name"
# GEO has no sex, age or BMI for this dataset, so we adjust for the date
# each library was prepared (a technical batch effect).
# Use character(0) to run without adjustment.
COVARIATES <- c(batch = "library_prep_date:ch1")

# ---- Analysis thresholds ---------------------------------------------
FDR_CUTOFF <- 0.05   # adjusted p-value cut-off for "significant"
MIN_COUNT  <- 10     # keep a gene if it has at least MIN_COUNT reads in at
                     # least as many samples as the smallest group
