# Small helper functions shared by the analysis scripts.
# You don't need to edit this file.

suppressPackageStartupMessages(library(here))

LOG_FILE <- here("results", "analysis_log.txt")

# File names NCBI uses for its generated counts and gene annotation
counts_file_name <- function() {
  paste0(GSE_ID, "_raw_counts_", NCBI_BUILD, "_NCBI.tsv.gz")
}
annotation_file_name <- function() {
  paste0("Human.", NCBI_BUILD, ".annot.tsv.gz")
}

# Print a message AND save it to results/analysis_log.txt, so every
# decision (exclusions, filters, numbers) is written down.
log_msg <- function(...) {
  msg <- paste0(..., collapse = "\n")
  message(msg)
  dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)
  cat(msg, "\n", file = LOG_FILE, append = TRUE, sep = "")
}

start_log <- function(step, fresh = FALSE) {
  dir.create(dirname(LOG_FILE), recursive = TRUE, showWarnings = FALSE)
  if (fresh && file.exists(LOG_FILE)) file.remove(LOG_FILE)
  log_msg("\n==== ", step, " (", format(Sys.time(), "%Y-%m-%d %H:%M"), ") ====")
}

# Stop with a clear message if a setting in config.R hasn't been filled in
check_filled <- function(value, setting) {
  if (length(value) == 0 || any(is.na(value)) || any(value == "FILL_ME_IN")) {
    stop("Please fill in ", setting, " in config.R first ",
         "(see the output of scripts/01_download_data.R).", call. = FALSE)
  }
}

# Turn a column of GEO text into a covariate R can use in a model.
# Mostly-numeric columns (e.g. "62" or "62 years") become numbers;
# anything else (e.g. "Male"/"Female") becomes a factor.
as_covariate <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "n/a", "na", "unknown", "Unknown",
             "not available", "--")] <- NA
  present <- !is.na(x)

  if (sum(present) > 0 && mean(grepl("[0-9]", x[present])) >= 0.9) {
    number <- rep(NA_real_, length(x))
    for (i in which(present)) {
      hit <- regmatches(x[i], regexpr("-?[0-9]+(\\.[0-9]+)?", x[i]))
      if (length(hit) == 1) number[i] <- as.numeric(hit)
    }
    return(number)
  }

  out <- x
  out[present] <- make.names(tolower(x[present]))
  factor(out)
}

# Model formula, e.g. ~ sex + age + condition
# (condition goes last so the main result is the condition effect)
build_design <- function() {
  terms <- c(names(COVARIATES), "condition")
  as.formula(paste("~", paste(terms, collapse = " + ")))
}

save_plot <- function(plot, filename, width = 7, height = 5) {
  dir.create(here("results", "figures"), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(here("results", "figures", filename), plot,
                  width = width, height = height, dpi = 300, bg = "white")
}
