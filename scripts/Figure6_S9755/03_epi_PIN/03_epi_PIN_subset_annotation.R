#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 module 03
# Subset Epi/PIN, recluster, final annotation: Epi1-Epi4 and PIN
# =============================================================
source("scripts/06_S9755/00_config.R")

obj <- readRDS(file.path(RDS_DIR, "02_S9755_broad_annotated.rds"))
assay <- get_spatial_assay(obj)
DefaultAssay(obj) <- assay

message("Subsetting broad Epi + PIN compartments...")
epi_pin_sub <- subset(obj, subset = major_group %in% c("Epi", "PIN"))
epi_pin_sub$epi_pin_cell_type <- as.character(epi_pin_sub$cell_type)
Idents(epi_pin_sub) <- "epi_pin_cell_type"
saveRDS(epi_pin_sub, file.path(RDS_DIR, "03a_S9755_EpiPIN_subset_annotated.rds"))

message("Epi/PIN reclustering. Historical filename says res2.2; actual FindClusters resolution=2.0")
epi_pin_sub <- FindVariableFeatures(epi_pin_sub, selection.method = "vst", nfeatures = EPIPIN_NFEATURES, verbose = FALSE)
epi_pin_sub <- ScaleData(epi_pin_sub, features = rownames(epi_pin_sub), verbose = FALSE)
epi_pin_sub <- RunPCA(epi_pin_sub, verbose = FALSE)
epi_pin_sub <- FindNeighbors(epi_pin_sub, dims = EPIPIN_DIMS, verbose = FALSE)
epi_pin_sub <- FindClusters(epi_pin_sub, resolution = EPIPIN_RESOLUTION, verbose = FALSE)
epi_pin_sub <- RunUMAP(epi_pin_sub, dims = EPIPIN_DIMS, reduction = "pca", verbose = FALSE)

cluster_vec <- as.character(epi_pin_sub$seurat_clusters)
cell_type_v2 <- dplyr::recode(
  cluster_vec,
  "0"  = "Epi1", "2"  = "Epi1",
  "8"  = "Epi2", "6"  = "Epi2", "1" = "Epi2", "12" = "Epi2",
  "4"  = "Epi2", "14" = "Epi2", "10" = "Epi2", "5" = "Epi2",
  "11" = "Epi3", "3" = "Epi3", "9" = "Epi3",
  "16" = "Epi4", "13" = "Epi4",
  "15" = "PIN", "17" = "PIN", "7" = "PIN",
  .default = NA_character_
)

epi_pin_sub$epi_pin_cell_type_v2 <- factor(cell_type_v2, levels = EPIPIN_LEVELS)
Idents(epi_pin_sub) <- "epi_pin_cell_type_v2"

write.csv(as.data.frame(table(epi_pin_sub$seurat_clusters, epi_pin_sub$epi_pin_cell_type_v2)),
          file.path(TABLE_DIR, "03_S9755_EpiPIN_cluster_to_final_annotation.csv"), row.names = FALSE)
write.csv(as.data.frame(table(epi_pin_sub$epi_pin_cell_type_v2)),
          file.path(TABLE_DIR, "03_S9755_EpiPIN_final_counts.csv"), row.names = FALSE)

save_tiff(plot_basic_dim(epi_pin_sub, "seurat_clusters", title = "S9755 Epi/PIN reclustered clusters"),
          file.path(PLOT_DIR, "03_S9755_EpiPIN_UMAP_clusters.tiff"), width = 10, height = 8)
save_tiff(plot_basic_dim(epi_pin_sub, "epi_pin_cell_type_v2", cols = EPIPIN_COLORS,
                          title = "S9755 final Epi/PIN annotation"),
          file.path(PLOT_DIR, "03_S9755_EpiPIN_UMAP_annotated_v2.tiff"), width = 8, height = 7)
save_tiff(SpatialDimPlot(epi_pin_sub, group.by = "epi_pin_cell_type_v2", pt.size.factor = 2.5) +
            scale_fill_manual(values = EPIPIN_COLORS) + ggtitle("S9755 final Epi/PIN spatial annotation"),
          file.path(PLOT_DIR, "03_S9755_EpiPIN_Spatial_annotated_v2.tiff"), width = 10, height = 10)

saveRDS(epi_pin_sub, file.path(RDS_DIR, "03_S9755_EpiPIN_subset_annotated_v2.rds"))
message("Saved: ", file.path(RDS_DIR, "03_S9755_EpiPIN_subset_annotated_v2.rds"))
