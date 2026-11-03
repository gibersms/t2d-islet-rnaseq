# Step 2 - Prepare the data
# ---------------------------------------------------------------------
# What this does:
#   1. Keeps only the two groups you are comparing (e.g. ND vs T2D)
#   2. Matches each sample's information to its column in the count table
#   3. Builds the covariates (if any) and drops samples with missing values
#   4. Removes genes with very low counts
# Every exclusion is written to results/analysis_log.txt.

source(here::here("config.R"))
source(here::here("scripts", "helpers.R"))

start_log("Step 2 - Prepare data", fresh = TRUE)
check_filled(GROUP_COLUMN, "GROUP_COLUMN")
check_filled(CONTROL_LABEL, "CONTROL_LABEL")
check_filled(CASE_LABEL, "CASE_LABEL")

raw_dir <- here("data", "raw")
sample_info <- read.csv(file.path(raw_dir, "sample_info_all_columns.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
counts_all <- as.matrix(read.delim(file.path(raw_dir, counts_file_name()),
                                   row.names = 1, check.names = FALSE))
annotation <- read.delim(file.path(raw_dir, annotation_file_name()),
                         quote = "", check.names = FALSE,
                         stringsAsFactors = FALSE)

log_msg("Dataset: ", GSE_ID)
log_msg("Count table: ", nrow(counts_all), " genes x ", ncol(counts_all), " samples")

# ---- 1. Keep the two groups ------------------------------------------
if (!GROUP_COLUMN %in% colnames(sample_info)) {
  stop("GROUP_COLUMN '", GROUP_COLUMN, "' is not in the sample information.\n",
       "Available columns: ",
       paste(grep(":ch1$", colnames(sample_info), value = TRUE), collapse = ", "),
       call. = FALSE)
}
group <- trimws(as.character(sample_info[[GROUP_COLUMN]]))
for (label in c(CONTROL_LABEL, CASE_LABEL)) {
  if (!label %in% group) {
    stop("No samples have the value '", label, "' in ", GROUP_COLUMN, ".\n",
         "Values found: ", paste(unique(group), collapse = ", "), call. = FALSE)
  }
}

group_counts <- table(group, useNA = "ifany")
log_msg("Samples per group in the full dataset:")
for (i in seq_along(group_counts)) {
  log_msg("  ", names(group_counts)[i], ": ", group_counts[[i]])
}

keep <- group %in% c(CONTROL_LABEL, CASE_LABEL)
info <- sample_info[keep, , drop = FALSE]
log_msg("Kept ", sum(keep), " samples (", CONTROL_LABEL, " + ", CASE_LABEL,
        "); excluded ", sum(!keep), " from other groups.")

# Optional extra filter, e.g. keep only samples in the authors' QC subset
if (!is.na(SUBSET_COLUMN)) {
  if (!SUBSET_COLUMN %in% colnames(info)) {
    stop("SUBSET_COLUMN '", SUBSET_COLUMN, "' is not in the sample information.",
         call. = FALSE)
  }
  in_subset <- trimws(as.character(info[[SUBSET_COLUMN]])) == SUBSET_VALUE
  in_subset[is.na(in_subset)] <- FALSE
  log_msg("Subset filter: kept samples where '", SUBSET_COLUMN, "' = ",
          SUBSET_VALUE, "; excluded ", sum(!in_subset), " samples.")
  info <- info[in_subset, , drop = FALSE]
}

# ---- 2. Match sample information to count columns -----------------------
if (anyDuplicated(info$geo_accession)) stop("Duplicated sample IDs found.")
no_counts <- setdiff(info$geo_accession, colnames(counts_all))
if (length(no_counts) > 0) {
  log_msg("Excluded ", length(no_counts), " samples that have no NCBI counts: ",
          paste(no_counts, collapse = ", "))
  info <- info[info$geo_accession %in% colnames(counts_all), , drop = FALSE]
}

samples <- data.frame(
  sample = info$geo_accession,
  condition = factor(make.names(trimws(info[[GROUP_COLUMN]])),
                     levels = make.names(c(CONTROL_LABEL, CASE_LABEL))),
  row.names = info$geo_accession,
  stringsAsFactors = FALSE
)

# ---- 3. Covariates ------------------------------------------------------
if (length(COVARIATES) > 0) {
  if (is.null(names(COVARIATES)) || any(names(COVARIATES) == "") ||
      any(names(COVARIATES) != make.names(names(COVARIATES)))) {
    stop("Give each covariate a simple short name in config.R, ",
         "e.g. COVARIATES <- c(sex = \"gender:ch1\")", call. = FALSE)
  }
  for (cov_name in names(COVARIATES)) {
    column <- COVARIATES[[cov_name]]
    check_filled(column, paste0("COVARIATES['", cov_name, "']"))
    if (!column %in% colnames(info)) {
      stop("Covariate column '", column, "' is not in the sample information.",
           call. = FALSE)
    }
    samples[[cov_name]] <- as_covariate(info[[column]])
    log_msg("Covariate '", cov_name, "' (from '", column, "') treated as ",
            if (is.numeric(samples[[cov_name]])) "a number" else "a category")
  }

  incomplete <- !complete.cases(samples[, names(COVARIATES), drop = FALSE])
  if (any(incomplete)) {
    log_msg("Excluded ", sum(incomplete), " samples with missing covariate values: ",
            paste(samples$sample[incomplete], collapse = ", "))
    samples <- samples[!incomplete, , drop = FALSE]
  }

  for (cov_name in names(COVARIATES)) {
    if (is.factor(samples[[cov_name]])) {
      samples[[cov_name]] <- droplevels(samples[[cov_name]])
      if (nlevels(samples[[cov_name]]) < 2) {
        stop("Covariate '", cov_name, "' has only one value after filtering; ",
             "remove it from COVARIATES.", call. = FALSE)
      }
    }
  }
}

log_msg("Final samples per group:")
final_counts <- table(samples$condition)
for (i in seq_along(final_counts)) {
  log_msg("  ", names(final_counts)[i], ": ", final_counts[[i]])
}
if (any(final_counts < 3)) stop("Each group needs at least 3 samples.", call. = FALSE)

# Describe the groups (useful for the README)
for (cov_name in names(COVARIATES)) {
  v <- samples[[cov_name]]
  if (is.numeric(v)) {
    s <- tapply(v, samples$condition,
                function(z) sprintf("%.1f (SD %.1f)", mean(z), sd(z)))
    log_msg("  ", cov_name, ", mean by group: ",
            paste(names(s), s, sep = " = ", collapse = "; "))
  } else {
    log_msg("  ", cov_name, " by group:")
    log_msg(capture.output(print(table(samples$condition, v))))
  }
}

# ---- 4. Count matrix in the same order as the sample table --------------
counts <- counts_all[, rownames(samples), drop = FALSE]
stopifnot(identical(colnames(counts), rownames(samples)))
if (any(is.na(counts))) stop("The count table contains missing values.")
if (any(counts < 0) || any(counts != round(counts))) {
  stop("Counts must be non-negative whole numbers (raw counts).")
}
mode(counts) <- "integer"

samples$library_size <- colSums(counts)
log_msg("Reads counted per sample (millions): min ",
        round(min(samples$library_size) / 1e6, 1), ", median ",
        round(median(samples$library_size) / 1e6, 1), ", max ",
        round(max(samples$library_size) / 1e6, 1))

# ---- 5. Remove genes with very low counts -------------------------------
# Genes with almost no reads carry no information and add to the
# multiple-testing burden.
min_group <- min(final_counts)
keep_genes <- rowSums(counts >= MIN_COUNT) >= min_group
counts_filtered <- counts[keep_genes, , drop = FALSE]
log_msg("Genes kept after low-count filter (>= ", MIN_COUNT, " reads in >= ",
        min_group, " samples): ", sum(keep_genes), " of ", nrow(counts))

# ---- 6. Gene annotation -------------------------------------------------
wanted <- intersect(c("GeneID", "Symbol", "Description", "GeneType"),
                    colnames(annotation))
if (!all(c("GeneID", "Symbol") %in% wanted)) {
  stop("The annotation file doesn't have GeneID and Symbol columns.")
}
gene_annotation <- annotation[, wanted]
gene_annotation$GeneID <- as.character(gene_annotation$GeneID)

# ---- Save ---------------------------------------------------------------
# Raw (unscaled) sample table for the README and plots
dir.create(here("results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(samples, here("results", "tables", "samples_used.csv"), row.names = FALSE)

# For the model, numeric covariates are centred and scaled (mean 0, SD 1),
# as the DESeq2 authors recommend
model_samples <- samples
for (cov_name in names(COVARIATES)) {
  if (is.numeric(model_samples[[cov_name]])) {
    model_samples[[cov_name]] <- as.numeric(scale(model_samples[[cov_name]]))
  }
}

proc_dir <- here("data", "processed")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(counts_filtered, file.path(proc_dir, "counts_filtered.rds"))
saveRDS(model_samples, file.path(proc_dir, "sample_table.rds"))
saveRDS(gene_annotation, file.path(proc_dir, "gene_annotation.rds"))
log_msg("Saved prepared data to data/processed/")
