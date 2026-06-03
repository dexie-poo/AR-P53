#!/usr/bin/env Rscript

.libPaths(c(
  "/gpfs/data/sunz04lab/leungd02/R/x86_64-pc-linux-gnu-library/4.4",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
  library(ggpubr)
  library(rstatix)
  library(pheatmap)
  library(grid)
  library(monocle3)
  library(viridis)
  library(fgsea)
  library(data.table)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(scales)
})

source("scripts/utils/helper_functions.R")

raw_dir  <- "/gpfs/data/sunz04lab/ST/20260324_Dexie/result_A1281_0427_24um/outs"

out_dir <- "/gpfs/scratch/leungd02/06032026/A1281_final"

rds_dir   <- file.path(out_dir, "rds")
table_dir <- file.path(out_dir, "tables")
plot_dir  <- file.path(out_dir, "plots")
gsea_dir  <- file.path(out_dir, "gsea")
mono_dir  <- file.path(out_dir, "monocle3")
log_dir   <- file.path(out_dir, "logs")

expr_dir <- file.path(plot_dir, "expression_plots")
vln_dir  <- file.path(plot_dir, "violin_plots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mono_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(vln_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.25

min_deg_for_gsea <- 10

gsea_min_size <- 5
gsea_max_size <- 500

gsea_nperm_simple <- 1000
gsea_param <- 0

MIN_FEATURES <- 50
MITO_THRESH <- 25
LIBSIZE_MIN <- 5

set.seed(1234)
