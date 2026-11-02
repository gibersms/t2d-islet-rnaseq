# Run the whole analysis from start to finish, then record the package
# versions used. Use this at the end to prove the workflow reproduces.

steps <- c("01_download_data.R",
           "02_prepare_data.R",
           "03_quality_checks.R",
           "04_differential_expression.R",
           "05_figures.R")

for (step in steps) {
  message("\n========== Running ", step, " ==========")
  source(here::here("scripts", step), local = new.env())
}

dir.create(here::here("docs"), showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), here::here("docs", "session_info.txt"))
message("\nAll steps finished. Package versions saved to docs/session_info.txt")
