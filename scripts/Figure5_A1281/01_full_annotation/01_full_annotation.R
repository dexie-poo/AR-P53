#!/usr/bin/env Rscript

source("scripts/Figure5_A1281/00_setup/00_config.R")

message("Loading raw processed RDS...")

integrated_p14 <- readRDS(
  "/gpfs/scratch/leungd02/04302026/A1281_24um/A1281_24um_04302026.rds"
)

spatial_assay <- "Spatial.024um"
DefaultAssay(integrated_p14) <- spatial_assay

message("Reclustering full object at resolution 0.8...")

integrated_p14 <- FindNeighbors(integrated_p14, verbose = FALSE)

integrated_p14 <- FindClusters(
  integrated_p14,
  resolution = 0.8,
  verbose = FALSE
)

if (!"umap" %in% names(integrated_p14@reductions)) {
  integrated_p14 <- RunUMAP(
    integrated_p14,
    dims = 1:30,
    reduction = "pca",
    verbose = FALSE
  )
}

cluster_vec <- as.character(integrated_p14$seurat_clusters)

cell_type_vec <- dplyr::recode(
  cluster_vec,
  "0"  = "FB1",
  "1"  = "FB2",
  "2"  = "FB3",
  "3"  = "PIN1",
  "4"  = "SM",
  "5"  = "FB4",
  "6"  = "SM",
  "7"  = "LumEpi1",
  "8"  = "SV",
  "9"  = "PIN2",
  "10" = "PIN3",
  "11" = "Myo",
  "12" = "TUMOR",
  "13" = "FB5",
  "14" = "Myo",
  "15" = "PIN4",
  "16" = "SV",
  "17" = "UrLE",
  "18" = "Immune",
  "19" = "SM",
  "20" = "LumEpi2",
  "21" = "BasalEpi"
)

integrated_p14 <- AddMetaData(
  integrated_p14,
  metadata = setNames(cell_type_vec, colnames(integrated_p14)),
  col.name = "cell_type"
)

cluster_order <- c(
  "LumEpi1", "LumEpi2",
  "BasalEpi",
  "PIN1", "PIN2", "PIN3", "PIN4",
  "TUMOR",
  "FB1", "FB2", "FB3", "FB4", "FB5",
  "SM",
  "Myo",
  "SV",
  "UrLE",
  "Immune"
)

integrated_p14$cell_type <- factor(
  integrated_p14$cell_type,
  levels = cluster_order
)

integrated_p14$major_group_EpiCombined <- dplyr::case_when(
  integrated_p14$cell_type %in% c(
    "LumEpi1", "LumEpi2", "BasalEpi",
    "PIN1", "PIN2", "PIN3", "PIN4",
    "TUMOR", "UrLE"
  ) ~ "Epi",
  integrated_p14$cell_type %in% c(
    "FB1", "FB2", "FB3", "FB4", "FB5"
  ) ~ "FB",
  integrated_p14$cell_type == "SM" ~ "SM",
  integrated_p14$cell_type == "Myo" ~ "Myo",
  integrated_p14$cell_type == "SV" ~ "SV",
  integrated_p14$cell_type == "Immune" ~ "Immune",
  TRUE ~ NA_character_
)

integrated_p14$major_group_EpiCombined <- factor(
  integrated_p14$major_group_EpiCombined,
  levels = c("Epi", "FB", "SM", "Myo", "SV", "Immune")
)

Idents(integrated_p14) <- "cell_type"

less_broad_colors <- setNames(
  custom_colors[seq_len(length(cluster_order))],
  cluster_order
)

major_colors <- c(
  "Epi" = "red",
  "FB" = "green",
  "SM" = "blue",
  "Myo" = "purple",
  "SV" = "orange",
  "Immune" = "cyan"
)

p_umap_celltype <- DimPlot(
  integrated_p14,
  reduction = "umap",
  group.by = "cell_type",
  cols = less_broad_colors,
  label = FALSE,
  pt.size = 0.8
) +
  theme_classic() +
  labs(title = "A1281 24um – less-broad cell type annotation")

save_tiff(
  p_umap_celltype,
  file.path(annot_plot_dir, "02_UMAP_less_broad_cell_type.tiff"),
  width = 12,
  height = 10
)

p_spatial_celltype <- SpatialDimPlot(
  integrated_p14,
  group.by = "cell_type",
  label = FALSE,
  pt.size.factor = 2.8
) +
  scale_fill_manual(values = less_broad_colors) +
  labs(title = "A1281 24um – less-broad cell type spatial")

save_tiff(
  p_spatial_celltype,
  file.path(annot_plot_dir, "02_Spatial_less_broad_cell_type.tiff"),
  width = 12,
  height = 10
)

p_umap_major <- DimPlot(
  integrated_p14,
  reduction = "umap",
  group.by = "major_group_EpiCombined",
  cols = major_colors,
  label = FALSE,
  pt.size = 0.8
) +
  theme_classic() +
  labs(title = "A1281 24um – major groups, Epi combined")

save_tiff(
  p_umap_major,
  file.path(annot_plot_dir, "02_UMAP_major_group_EpiCombined.tiff"),
  width = 12,
  height = 10
)

p_spatial_major <- SpatialDimPlot(
  integrated_p14,
  group.by = "major_group_EpiCombined",
  label = FALSE,
  pt.size.factor = 2.8
) +
  scale_fill_manual(values = major_colors) +
  labs(title = "A1281 24um – major groups, Epi combined spatial")

save_tiff(
  p_spatial_major,
  file.path(annot_plot_dir, "02_Spatial_major_group_EpiCombined.tiff"),
  width = 12,
  height = 10
)

write.csv(
  as.data.frame(table(integrated_p14$cell_type)),
  file.path(table_dir, "02_less_broad_cell_type_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(integrated_p14$major_group_EpiCombined)),
  file.path(table_dir, "02_major_group_EpiCombined_counts.csv"),
  row.names = FALSE
)

saveRDS(
  integrated_p14,
  file.path(rds_dir, "02_full_annotated_less_broad_and_major.rds")
)

message("Saved: 02_full_annotated_less_broad_and_major.rds")
