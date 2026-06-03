#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

obj <- readRDS(file.path(rds_dir, "05_visualization_complete.rds"))

saveRDS(
  obj,
  file.path(rds_dir, "06_final_processed.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)
