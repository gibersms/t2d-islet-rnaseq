# Step 1 - Download the data from GEO
# ---------------------------------------------------------------------
# What this does:
#   1. Downloads the sample information (diagnosis, sex, age, ... for
#      each sample) and saves it as data/raw/sample_info_all_columns.csv
#   2. Downloads the raw gene counts that NCBI generated for this study
#      (one row per gene, one column per sample)
#   3. Downloads NCBI's gene annotation table (turns Gene IDs into gene
#      symbols such as INS or GCG)
#
# After it runs, read the printed list of sample characteristics and
# fill in GROUP_COLUMN, CONTROL_LABEL and CASE_LABEL in config.R.

source(here::here("config.R"))
source(here::here("scripts", "helpers.R"))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

options(timeout = 1200)  # downloads can be slow
raw_dir <- here("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Sample information --------------------------------------------
message("Downloading sample information for ", GSE_ID, " ...")
gse <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE, destdir = raw_dir)

# A GEO series can be split across several sequencing machines
# (platforms). Combine the sample tables from all of them.
sample_tables <- lapply(gse, function(es) {
  p <- pData(es)
  p[] <- lapply(p, as.character)
  p
})
all_columns <- unique(unlist(lapply(sample_tables, colnames)))
sample_tables <- lapply(sample_tables, function(p) {
  for (missing_col in setdiff(all_columns, colnames(p))) p[[missing_col]] <- NA
  p[, all_columns, drop = FALSE]
})
names(sample_tables) <- NULL
sample_info <- do.call(rbind, sample_tables)
rownames(sample_info) <- sample_info$geo_accession

write.csv(sample_info, file.path(raw_dir, "sample_info_all_columns.csv"),
          row.names = FALSE)
message("Saved information for ", nrow(sample_info), " samples.")

# Show the characteristics so you can choose the settings in config.R
char_cols <- grep(":ch1$", colnames(sample_info), value = TRUE)
message("\nSample characteristics (column name -> values):")
for (col in char_cols) {
  values <- sort(unique(sample_info[[col]]))
  more <- if (length(values) > 8) {
    paste0("  ... (", length(values), " different values)")
  } else ""
  message("  '", col, "' -> ", paste(head(values, 8), collapse = " | "), more)
}

# ---- 2. NCBI-generated raw counts -------------------------------------
counts_file <- file.path(raw_dir, counts_file_name())
counts_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
                     "?type=rnaseq_counts&acc=", GSE_ID,
                     "&format=file&file=", counts_file_name())

if (!file.exists(counts_file)) {
  message("\nDownloading NCBI raw counts ...")
  tryCatch(download.file(counts_url, counts_file, mode = "wb"),
           error = function(e) message("Download failed: ", conditionMessage(e)))
}

# Check we really got a count table (and not an error page)
first_line <- ""
if (file.exists(counts_file) && file.size(counts_file) > 1000) {
  con <- gzfile(counts_file)
  first_line <- readLines(con, n = 1, warn = FALSE)
  close(con)
  if (length(first_line) == 0) first_line <- ""
}
if (!startsWith(first_line, "GeneID")) {
  if (file.exists(counts_file)) file.remove(counts_file)
  stop("NCBI-generated raw counts were not found for ", GSE_ID, ".\n",
       "Open https://www.ncbi.nlm.nih.gov/geo/download/?acc=", GSE_ID, "\n",
       "and look for 'NCBI-generated data' > 'Series RNA-seq raw counts matrix'.\n",
       "If the file name there uses a different build (e.g. GRCh38.p14), change\n",
       "NCBI_BUILD in config.R. If there is no such file, see 'Plan B' in the guide.",
       call. = FALSE)
}

# ---- 3. Gene annotation ------------------------------------------------
annot_file <- file.path(raw_dir, annotation_file_name())
annot_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
                    "?format=file&type=rnaseq_counts&file=",
                    annotation_file_name())
if (!file.exists(annot_file)) {
  message("Downloading gene annotation ...")
  download.file(annot_url, annot_file, mode = "wb")
}

counts_preview <- read.delim(counts_file, nrows = 3, check.names = FALSE)
message("\nDone. The counts file has ", ncol(counts_preview) - 1, " samples.")
message("Next: fill in GROUP_COLUMN, CONTROL_LABEL and CASE_LABEL in config.R.")
