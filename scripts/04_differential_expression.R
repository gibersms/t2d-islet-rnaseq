# Step 4 - Differential expression with DESeq2
# ---------------------------------------------------------------------
# What this does:
#   1. Runs DESeq2: estimates size factors (sequencing-depth
#      normalisation), gene-wise dispersion, fits a negative binomial
#      model per gene and tests the condition effect (Wald test)
#   2. Adjusts p-values for testing thousands of genes
#      (Benjamini-Hochberg false discovery rate)
#   3. Shrinks log2 fold changes (apeglm) so noisy, low-count genes don't
#      show misleadingly large changes
#   4. Saves a full results table and a top-25 table

source(here::here("config.R"))
source(here::here("scripts", "helpers.R"))
suppressPackageStartupMessages({
  library(DESeq2)
})

start_log("Step 4 - Differential expression")
proc_dir <- here("data", "processed")
dds <- readRDS(file.path(proc_dir, "dds_unfitted.rds"))
gene_annotation <- readRDS(file.path(proc_dir, "gene_annotation.rds"))

control <- levels(dds$condition)[1]
case <- levels(dds$condition)[2]

# ---- 1. Fit the model ---------------------------------------------------
dds <- DESeq(dds)
sf <- sizeFactors(dds)
log_msg("Size factors range from ", round(min(sf), 2), " to ", round(max(sf), 2))

coef_name <- paste0("condition_", case, "_vs_", control)
if (!coef_name %in% resultsNames(dds)) {
  stop("Expected coefficient ", coef_name, " not found. Available: ",
       paste(resultsNames(dds), collapse = ", "))
}

# ---- 2. Test and adjust for multiple testing ----------------------------
res <- results(dds, name = coef_name, alpha = FDR_CUTOFF)
log_msg("DESeq2 summary for ", case, " vs ", control, ":")
log_msg(capture.output(summary(res)))

# ---- 3. Shrink fold changes ---------------------------------------------
res_shrunk <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm",
                        quiet = TRUE)
stopifnot(identical(rownames(res), rownames(res_shrunk)))

# ---- 4. Results tables --------------------------------------------------
results_table <- data.frame(
  GeneID = rownames(res),
  baseMean = res$baseMean,
  log2FC = res$log2FoldChange,
  lfcSE = res$lfcSE,
  log2FC_shrunk = res_shrunk$log2FoldChange,
  pvalue = res$pvalue,
  padj = res$padj,
  stringsAsFactors = FALSE
)
results_table <- merge(gene_annotation, results_table, by = "GeneID",
                       all.y = TRUE, sort = FALSE)
no_symbol <- is.na(results_table$Symbol) | results_table$Symbol == ""
results_table$Symbol[no_symbol] <- results_table$GeneID[no_symbol]
results_table <- results_table[order(results_table$padj, results_table$pvalue), ]
rownames(results_table) <- NULL

significant <- !is.na(results_table$padj) & results_table$padj < FDR_CUTOFF
n_up <- sum(significant & results_table$log2FC > 0)
n_down <- sum(significant & results_table$log2FC < 0)
log_msg("Genes tested: ", nrow(results_table),
        " (padj is NA for ", sum(is.na(results_table$padj)),
        " genes removed by DESeq2's independent filtering or outlier detection)")
log_msg("Genes with FDR < ", FDR_CUTOFF, ": ", sum(significant), " (",
        n_up, " higher in ", case, ", ", n_down, " lower in ", case, ")")
log_msg("Of these, genes with shrunken |log2FC| >= 1 (at least 2-fold): ",
        sum(significant & abs(results_table$log2FC_shrunk) >= 1))

tables_dir <- here("results", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
rounded <- results_table
for (col in c("baseMean", "log2FC", "lfcSE", "log2FC_shrunk")) {
  rounded[[col]] <- round(rounded[[col]], 3)
}
rounded$pvalue <- signif(rounded$pvalue, 3)
rounded$padj <- signif(rounded$padj, 3)
write.csv(rounded, file.path(tables_dir, "deseq2_results_all_genes.csv"),
          row.names = FALSE)

top <- head(rounded[!is.na(rounded$padj), ], 25)
top$direction <- ifelse(top$log2FC > 0, paste("Higher in", case),
                        paste("Lower in", case))
top$significant <- top$padj < FDR_CUTOFF
top_cols <- intersect(c("Symbol", "Description", "baseMean", "log2FC_shrunk",
                        "padj", "significant", "direction"), colnames(top))
write.csv(top[, top_cols], file.path(tables_dir, "top25_genes.csv"),
          row.names = FALSE)

saveRDS(dds, file.path(proc_dir, "dds_fitted.rds"))
saveRDS(res_shrunk, file.path(proc_dir, "res_shrunk.rds"))
saveRDS(results_table, file.path(proc_dir, "results_table.rds"))
log_msg("Saved results tables to results/tables/")
