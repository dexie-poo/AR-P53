#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

obj <- readRDS(file.path(out_dir, "A1281_epi_sub_reannotated_visualization_complete.rds"))

saveRDS(
  obj,
  file.path(out_dir, "A1281_final_processed.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)