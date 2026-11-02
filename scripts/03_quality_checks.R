# Step 3 - Quality checks
# ---------------------------------------------------------------------
# What this does:
#   1. Plots library size (total reads counted) for each sample
#   2. Applies a variance-stabilising transformation (VST) so samples can
#      be compared on one scale
#   3. Runs PCA on the 500 most variable genes and flags unusual samples
#      (flagged samples are NOT removed automatically - look at them)
#   4. Draws a sample-to-sample distance heatmap

source(here::here("config.R"))
source(here::here("scripts", "helpers.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

start_log("Step 3 - Quality checks")
proc_dir <- here("data", "processed")
counts <- readRDS(file.path(proc_dir, "counts_filtered.rds"))
samples <- readRDS(file.path(proc_dir, "sample_table.rds"))
samples_raw <- read.csv(here("results", "tables", "samples_used.csv"),
                        stringsAsFactors = FALSE)
samples_raw$condition <- factor(samples_raw$condition, levels = levels(samples$condition))

design_formula <- build_design()
log_msg("Model design: ", paste(deparse(design_formula), collapse = ""))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = samples,
                              design = design_formula)

# ---- 1. Library sizes ---------------------------------------------------
lib <- samples_raw
lib$millions <- lib$library_size / 1e6
lib$sample <- factor(lib$sample, levels = lib$sample[order(lib$condition, lib$millions)])
p <- ggplot(lib, aes(sample, millions, fill = condition)) +
  geom_col() +
  labs(x = "Samples (sorted)", y = "Reads counted (millions)", fill = NULL,
       title = "Library size per sample") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())
save_plot(p, "01_library_sizes.png")

low <- lib$sample[lib$millions < 0.25 * median(lib$millions)]
log_msg(if (length(low) == 0) {
  "No sample has fewer than 25% of the median number of reads."
} else {
  paste0("Samples with < 25% of the median reads (check them): ",
         paste(low, collapse = ", "))
})

# ---- 2. Variance-stabilising transformation -----------------------------
# blind = TRUE: the transformation ignores the groups, which is what you
# want for an unbiased quality check
vsd <- vst(dds, blind = TRUE)

# ---- 3. PCA -------------------------------------------------------------
mat <- assay(vsd)
top_genes <- head(order(apply(mat, 1, var), decreasing = TRUE), 500)
pca <- prcomp(t(mat[top_genes, ]))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

pca_df <- cbind(as.data.frame(pca$x[, 1:2]),
                samples_raw[match(rownames(pca$x), samples_raw$sample), ])

# Flag samples more than 3 standard deviations from the centre on PC1 or PC2
far <- abs(as.numeric(scale(pca_df$PC1))) > 3 | abs(as.numeric(scale(pca_df$PC2))) > 3
pca_df$flag <- ifelse(far, pca_df$sample, "")
log_msg(if (!any(far)) {
  "PCA: no sample is more than 3 SD from the centre on PC1 or PC2."
} else {
  paste0("PCA: samples more than 3 SD from the centre (kept, but look at them): ",
         paste(pca_df$sample[far], collapse = ", "))
})
log_msg("PC1 explains ", pct[1], "% and PC2 ", pct[2], "% of the variance ",
        "(500 most variable genes).")

pca_plot <- function(colour_by, title) {
  ggplot(pca_df, aes(PC1, PC2, colour = .data[[colour_by]])) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_text_repel(aes(label = flag), colour = "black", size = 3,
                    max.overlaps = Inf) +
    labs(x = paste0("PC1 (", pct[1], "% of variance)"),
         y = paste0("PC2 (", pct[2], "% of variance)"),
         colour = colour_by, title = title) +
    theme_minimal()
}
save_plot(pca_plot("condition", "PCA of the 500 most variable genes"),
          "02_pca_condition.png")

# Also colour the PCA by each covariate: if samples separate by sex or age
# rather than by diagnosis, that is worth knowing
for (cov_name in names(COVARIATES)) {
  save_plot(pca_plot(cov_name, paste("PCA coloured by", cov_name)),
            paste0("02_pca_", cov_name, ".png"))
}

# ---- 4. Sample-to-sample distances -------------------------------------
sample_dists <- dist(t(mat))
annotation_col <- data.frame(condition = samples$condition,
                             row.names = colnames(dds))
pheatmap(as.matrix(sample_dists),
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col,
         show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("#08306B", "#F7FBFF"))(100),
         main = "Sample-to-sample distances (darker = more similar)",
         filename = here("results", "figures", "03_sample_distances.png"),
         width = 7, height = 6)

saveRDS(dds, file.path(proc_dir, "dds_unfitted.rds"))
log_msg("Saved figures 01-03 to results/figures/")
