#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 module 02
# Broad reclustering and annotation: Epi, PIN, FB, SM, UreSM, UreLE, Imm
# =============================================================
source("scripts/06_S9755/00_config.R")

input_rds <- file.path(RDS_DIR, "01_S9755_24um_lognorm_processed.rds")
if (!file.exists(input_rds) && file.exists(HISTORICAL_PROCESSED_RDS)) {
  message("Using historical processed RDS because module 01 output was not found: ", HISTORICAL_PROCESSED_RDS)
  input_rds <- HISTORICAL_PROCESSED_RDS
}

obj <- readRDS(input_rds)
assay <- get_spatial_assay(obj)
DefaultAssay(obj) <- assay

message("Broad reclustering using historical parameters: nfeatures=", BROAD_NFEATURES,
        ", dims=1:30, resolution=", BROAD_RESOLUTION)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = BROAD_NFEATURES, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = BROAD_DIMS, verbose = FALSE)
obj <- FindClusters(obj, resolution = BROAD_RESOLUTION, verbose = FALSE)
obj <- RunUMAP(obj, dims = BROAD_DIMS, reduction = "pca", verbose = FALSE)

cluster_vec <- as.character(obj$seurat_clusters)
cell_type_vec <- dplyr::recode(
  cluster_vec,
  "0"  = "Epi1",
  "1"  = "Epi2",
  "2"  = "SM1",
  "3"  = "Epi3",
  "4"  = "SM2",
  "5"  = "FB1",
  "6"  = "Epi4",
  "7"  = "Epi5",
  "8"  = "UreSM",
  "9"  = "Epi6",
  "10" = "FB1",
  "11" = "UreLE",
  "12" = "SM3",
  "13" = "UreSM",
  "14" = "Epi7",
  "15" = "UreSM",
  "16" = "Imm",
  "17" = "UreSM",
  "18" = "Epi8",
  "19" = "UreSM",
  "20" = "UreSM",
  "21" = "FB2",
  "22" = "FB3",
  "23" = "Epi9",
  "24" = "PIN",
  "25" = "FB4",
  .default = NA_character_
)

obj$cell_type <- factor(cell_type_vec, levels = c(
  "Epi1", "Epi2", "Epi3", "Epi4", "Epi5", "Epi6", "Epi7", "Epi8", "Epi9",
  "PIN", "FB1", "FB2", "FB3", "FB4", "SM1", "SM2", "SM3", "UreSM", "UreLE", "Imm"
))

obj$major_group <- dplyr::case_when(
  grepl("^Epi", as.character(obj$cell_type)) ~ "Epi",
  as.character(obj$cell_type) == "PIN" ~ "PIN",
  grepl("^FB", as.character(obj$cell_type)) ~ "FB",
  grepl("^SM", as.character(obj$cell_type)) ~ "SM",
  as.character(obj$cell_type) == "UreSM" ~ "UreSM",
  as.character(obj$cell_type) == "UreLE" ~ "UreLE",
  as.character(obj$cell_type) == "Imm" ~ "Imm",
  TRUE ~ NA_character_
)
obj$major_group <- factor(obj$major_group, levels = BROAD_LEVELS)
Idents(obj) <- "cell_type"

write.csv(as.data.frame(table(obj$seurat_clusters, obj$cell_type)),
          file.path(TABLE_DIR, "02_S9755_broad_cluster_to_cell_type_table.csv"), row.names = FALSE)
write.csv(as.data.frame(table(obj$major_group)),
          file.path(TABLE_DIR, "02_S9755_major_group_counts.csv"), row.names = FALSE)

save_tiff(plot_basic_dim(obj, "cell_type", title = "S9755 broad detailed annotation"),
          file.path(PLOT_DIR, "02_S9755_UMAP_broad_detailed_annotation.tiff"), width = 10, height = 8)
save_tiff(plot_basic_dim(obj, "major_group", cols = BROAD_COLORS, title = "S9755 broad annotation"),
          file.path(PLOT_DIR, "02_S9755_UMAP_broad_major_annotation.tiff"), width = 8, height = 7)
save_tiff(SpatialDimPlot(obj, group.by = "major_group", pt.size.factor = 2.5) +
            scale_fill_manual(values = BROAD_COLORS) + ggtitle("S9755 broad spatial annotation"),
          file.path(PLOT_DIR, "02_S9755_Spatial_broad_major_annotation.tiff"), width = 10, height = 10)

saveRDS(obj, file.path(RDS_DIR, "02_S9755_broad_annotated.rds"))
message("Saved: ", file.path(RDS_DIR, "02_S9755_broad_annotated.rds"))
