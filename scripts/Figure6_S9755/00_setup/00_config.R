#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 modular pipeline configuration
# Project: ARP53 publication repository
# Purpose: Central paths, parameters, libraries, helpers
# =============================================================

.libPaths(c(
  "/gpfs/data/sunz04lab/leungd02/R/x86_64-pc-linux-gnu-library/4.4",
  .libPaths()
))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(tidyverse)
  library(Matrix)
})

set.seed(42)

# ---- Repository paths ----
PROJECT_DIR <- Sys.getenv("ARP53_PROJECT_DIR", unset = "/gpfs/scratch/leungd02/ARP53")
FIGURE_ID   <- "Figure6_S9755"
OUT_DIR     <- file.path(PROJECT_DIR, "outputs", FIGURE_ID)
RDS_DIR     <- file.path(OUT_DIR, "rds")
PLOT_DIR    <- file.path(OUT_DIR, "plots")
TABLE_DIR   <- file.path(OUT_DIR, "tables")
LOG_DIR     <- file.path(OUT_DIR, "logs")

for (d in c(OUT_DIR, RDS_DIR, PLOT_DIR, TABLE_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Raw data / historical source paths ----
RAW_10X_DIR <- "/gpfs/data/sunz04lab/ST/20260324_Dexie/result_S9755_0427_24um/outs"
BIN_SIZE    <- 24

# Historical processed object inspected by user:
# 5410 features x 25779 spots, all global_outliers FALSE, PCA/UMAP present.
HISTORICAL_PROCESSED_RDS <- "/gpfs/scratch/leungd02/04302026/S9755_24um/S9755_24um_04302026.rds"

# ---- QC parameters confirmed from raw-processing lineage ----
MIN_FEATURES <- 50
MITO_THRESH  <- 25
LIBSIZE_MIN  <- 5

# ---- Processing parameters ----
GENE_PCT_MIN        <- 1
GENE_MEAN_MIN       <- 0.05
N_HVG_RAW           <- 2000
RAW_NPCS            <- 50
RAW_CLUSTER_DIMS    <- 1:30
RAW_CLUSTER_RESES   <- c(0.2, 0.3, 0.5, 0.8, 1.0)
RAW_FINAL_RES       <- 0.5

# ---- Broad reclustering parameters ----
BROAD_NFEATURES     <- 1500
BROAD_DIMS          <- 1:30
BROAD_RESOLUTION    <- 1.8

# ---- Epi/PIN reclustering parameters ----
EPIPIN_NFEATURES    <- 3000
EPIPIN_DIMS         <- 1:30
EPIPIN_RESOLUTION   <- 2.0  # historical filename says res2.2, actual FindClusters uses resolution=2

SPATIAL_ASSAY <- "Spatial.024um"

# ---- Annotation levels ----
BROAD_LEVELS <- c("Epi", "PIN", "FB", "SM", "UreSM", "UreLE", "Imm")
EPIPIN_LEVELS <- c("Epi1", "Epi2", "Epi3", "Epi4", "PIN")

BROAD_COLORS <- c(
  "Epi"   = "red",
  "PIN"   = "black",
  "FB"    = "dodgerblue",
  "SM"    = "orange",
  "UreSM" = "purple",
  "UreLE" = "green",
  "Imm"   = "grey40"
)

EPIPIN_COLORS <- c(
  "Epi1" = "red",
  "Epi2" = "blue",
  "Epi3" = "green",
  "Epi4" = "orange",
  "PIN"  = "black"
)

CUSTOM_CLUSTER_COLORS <- c(
  "red", "green", "blue", "yellow", "cyan", "dodgerblue", "lightblue", "gold",
  "orange", "purple", "lightgreen", "magenta", "pink", "brown", "grey", "darkgreen",
  "darkblue", "darkred", "turquoise", "salmon", "violet", "tan", "steelblue", "tomato",
  "black", "lightcoral", "darksalmon", "azure2", "aquamarine", "chartreuse3", "coral3",
  "cornflowerblue", "darkgoldenrod1", "darkolivegreen3", "deeppink", "firebrick2",
  "khaki3", "mediumorchid", "navy", "olivedrab3", "plum3", "rosybrown", "seagreen3",
  "slateblue2", "springgreen3", "wheat3", "yellowgreen", "orchid3", "sienna3", "palevioletred3"
)

save_tiff <- function(plot, filename, width = 12, height = 10, dpi = 300) {
  ggsave(filename = filename, plot = plot, device = "tiff", width = width, height = height,
         dpi = dpi, compression = "lzw", limitsize = FALSE)
}

get_spatial_assay <- function(obj) {
  assays <- Assays(obj)
  raw_assays <- assays[!grepl("^SCT", assays)]
  if (SPATIAL_ASSAY %in% assays) return(SPATIAL_ASSAY)
  if ("Spatial" %in% raw_assays) return("Spatial")
  raw_assays[1]
}

get_counts <- function(obj, assay = get_spatial_assay(obj)) {
  tryCatch(GetAssayData(obj, assay = assay, layer = "counts"),
           error = function(e) GetAssayData(obj, assay = assay, slot = "counts"))
}

get_data <- function(obj, assay = get_spatial_assay(obj)) {
  tryCatch(GetAssayData(obj, assay = assay, layer = "data"),
           error = function(e) GetAssayData(obj, assay = assay, slot = "data"))
}

plot_basic_dim <- function(obj, group.by, cols = NULL, title = NULL, pt.size = 0.6) {
  DimPlot(obj, reduction = "umap", group.by = group.by, label = FALSE, pt.size = pt.size, cols = cols) +
    theme_classic() + labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
}
