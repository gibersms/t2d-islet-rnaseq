# Step 5 - Figures
# ---------------------------------------------------------------------
# What this does:
#   04 P-value histogram (a diagnostic: is the test behaving sensibly?)
#   05 MA plot (fold change vs average expression)
#   06 Volcano plot (fold change vs significance, top genes labelled)
#   07 Heatmap of the top 30 significant genes
#   08 Normalised counts for the top 6 genes, per group

source(here::here("config.R"))
source(here::here("scripts", "helpers.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

set.seed(1)  # makes the jittered points in figure 08 reproducible
start_log("Step 5 - Figures")
proc_dir <- here("data", "processed")
dds <- readRDS(file.path(proc_dir, "dds_fitted.rds"))
res_shrunk <- readRDS(file.path(proc_dir, "res_shrunk.rds"))
results_table <- readRDS(file.path(proc_dir, "results_table.rds"))

control <- levels(dds$condition)[1]
case <- levels(dds$condition)[2]
tested <- results_table[!is.na(results_table$padj), ]
significant_ids <- tested$GeneID[tested$padj < FDR_CUTOFF]

# ---- 04 P-value histogram -----------------------------------------------
# Expect a flat spread with a peak near 0 if there are real differences.
p <- ggplot(results_table[!is.na(results_table$pvalue), ], aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "grey40") +
  labs(x = "Raw p-value", y = "Number of genes",
       title = "Distribution of p-values") +
  theme_minimal()
save_plot(p, "04_pvalue_histogram.png")

# ---- 05 MA plot -----------------------------------------------------------
png(here("results", "figures", "05_ma_plot.png"),
    width = 7, height = 5, units = "in", res = 300)
plotMA(res_shrunk, alpha = FDR_CUTOFF, ylim = c(-3, 3),
       main = paste0("MA plot: ", case, " vs ", control, " (shrunken log2FC)"))
dev.off()

# ---- 06 Volcano plot ------------------------------------------------------
higher <- paste("Higher in", case)
lower <- paste("Lower in", case)
tested$status <- "Not significant"
tested$status[tested$padj < FDR_CUTOFF & tested$log2FC > 0] <- higher
tested$status[tested$padj < FDR_CUTOFF & tested$log2FC < 0] <- lower
tested$status <- factor(tested$status, levels = c(higher, lower, "Not significant"))
tested$minus_log10_p <- -log10(pmax(tested$pvalue, .Machine$double.xmin))
labelled <- head(tested[tested$padj < FDR_CUTOFF, ], 15)

p <- ggplot(tested, aes(log2FC_shrunk, minus_log10_p, colour = status)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_text_repel(data = labelled, aes(label = Symbol), colour = "black",
                  size = 3, max.overlaps = Inf) +
  scale_colour_manual(values = setNames(c("#B2182B", "#2166AC", "grey75"),
                                        c(higher, lower, "Not significant")),
                      drop = FALSE) +
  labs(x = paste0("log2 fold change, ", case, " vs ", control, " (shrunken)"),
       y = "-log10(p-value)", colour = NULL,
       title = paste0("Differential expression: ", case, " vs ", control),
       subtitle = paste0(length(significant_ids), " genes with FDR < ", FDR_CUTOFF,
                         "; top 15 labelled")) +
  theme_minimal()
save_plot(p, "06_volcano.png", width = 8, height = 6)

# ---- 07 Heatmap of top significant genes ----------------------------------
top_ids <- head(significant_ids, 30)
if (length(top_ids) >= 2) {
  # blind = FALSE: here we want the transformation to use the model design
  vsd <- vst(dds, blind = FALSE)
  mat <- assay(vsd)[top_ids, , drop = FALSE]
  mat <- t(scale(t(mat)))          # z-score each gene across samples
  mat[mat > 3] <- 3
  mat[mat < -3] <- -3
  rownames(mat) <- make.unique(results_table$Symbol[match(top_ids, results_table$GeneID)])
  sample_order <- order(dds$condition)
  annotation_col <- data.frame(condition = dds$condition, row.names = colnames(dds))
  pheatmap(mat[, sample_order],
           cluster_cols = FALSE,
           annotation_col = annotation_col,
           show_colnames = FALSE,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(-3, 3, length.out = 101),
           main = paste0("Top ", length(top_ids), " genes (z-scored VST expression)"),
           filename = here("results", "figures", "07_heatmap_top_genes.png"),
           width = 8, height = 7)
} else {
  log_msg("Fewer than 2 significant genes, so the heatmap was skipped.")
}

# ---- 08 Counts for the top 6 genes ----------------------------------------
norm_counts <- counts(dds, normalized = TRUE)
top6 <- head(tested$GeneID, 6)
plot_df <- do.call(rbind, lapply(top6, function(g) {
  data.frame(gene = results_table$Symbol[match(g, results_table$GeneID)],
             condition = dds$condition,
             count = norm_counts[g, ])
}))
plot_df$gene <- factor(plot_df$gene, levels = unique(plot_df$gene))

p <- ggplot(plot_df, aes(condition, count + 0.5, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  scale_y_log10() +
  facet_wrap(~ gene, scales = "free_y") +
  labs(x = NULL, y = "Normalised count (log scale)",
       title = "Top 6 genes by adjusted p-value") +
  theme_minimal() +
  theme(legend.position = "none")
save_plot(p, "08_top_genes_counts.png", width = 8, height = 6)

log_msg("Saved figures 04-08 to results/figures/")
