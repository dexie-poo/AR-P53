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

source("scripts/Figure5_A1281/utils/helper_functions.R")

project_dir <- "/gpfs/scratch/leungd02/ARP53"

raw_dir <- "/gpfs/data/sunz04lab/ST/20260324_Dexie/result_A1281_0427_24um/outs"

out_dir   <- file.path(project_dir, "outputs")
rds_dir   <- file.path(out_dir, "rds")
table_dir <- file.path(out_dir, "tables")
plot_dir  <- file.path(out_dir, "plots")
gsea_dir  <- file.path(out_dir, "gsea")
mono_dir  <- file.path(out_dir, "monocle3")
log_dir   <- file.path(project_dir, "logs")

qc_plot_dir        <- file.path(plot_dir, "qc")
raw_plot_dir       <- file.path(plot_dir, "raw")
annot_plot_dir     <- file.path(plot_dir, "annotation")
reclust_plot_dir   <- file.path(plot_dir, "recluster")
monocle_plot_dir   <- file.path(plot_dir, "monocle3")
gsea_plot_dir      <- file.path(plot_dir, "gsea")
replot_dir         <- file.path(plot_dir, "replots")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mono_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annot_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reclust_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(monocle_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(replot_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(mono_dir, "RDS"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mono_dir, "Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mono_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)

spatial_assay <- "Spatial.024um"

MIN_FEATURES <- 50
MITO_THRESH  <- 25
LIBSIZE_MIN  <- 5

padj_cutoff       <- 0.25
min_deg_for_gsea  <- 10
gsea_min_size     <- 5
gsea_max_size     <- 500
gsea_nperm_simple <- 1000
gsea_param        <- 0

custom_colors <- c(
  "red", "green", "blue", "yellow",
  "cyan", "dodgerblue", "lightblue", "gold",
  "orange", "purple", "lightgreen", "magenta",
  "pink", "brown", "grey", "darkgreen",
  "darkblue", "darkred", "turquoise", "salmon",
  "violet", "tan", "steelblue", "tomato",
  "black", "lightcoral", "darksalmon", "azure2",
  "aquamarine", "chartreuse3", "coral3", "cornflowerblue",
  "darkgoldenrod1", "darkolivegreen3", "deeppink", "firebrick2",
  "khaki3", "mediumorchid", "navy", "olivedrab3",
  "plum3"
)

set.seed(1234)
