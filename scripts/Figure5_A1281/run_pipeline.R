#!/usr/bin/env Rscript

message("===== Running Figure 5 A1281 pipeline =====")

source("scripts/Figure5_A1281/00_setup/00_config.R")
source("scripts/Figure5_A1281/01_raw_processing/01_process_raw_lognorm.R")
source("scripts/Figure5_A1281/01_full_annotation/01_full_annotation.R")
source("scripts/Figure5_A1281/02_epi_intermediate/02_epi_intermediate.R")
source("scripts/Figure5_A1281/03_final_epi_PIN/03_final_epi_PIN.R")
source("scripts/Figure5_A1281/04_monocle3/04_monocle3.R")
source("scripts/Figure5_A1281/05_GSEA/05_GSEA.R")

message("===== Figure 5 A1281 pipeline complete =====")
