#!/usr/bin/env Rscript

source("scripts/Figure5_A1281/00_setup/00_config.R")

# =============================================================
# Module 4: Final epithelial reclustering → Epi/PIN1–PIN6
# Input:  outputs/rds/03_epi_intermediate_annotated.rds
# Output: outputs/rds/04_epi_PIN_reannotated.rds
# =============================================================

message("Loading Module 3 epithelial intermediate object...")

epi_sub <- readRDS(
  file.path(rds_dir, "03_epi_intermediate_annotated.rds")
)

message("Object loaded: ", ncol(epi_sub), " cells")

DefaultAssay(epi_sub) <- spatial_assay

# ── 1. Reclustering: exact retry-confirmed parameters ─────────

message("Running FindVariableFeatures ...")
epi_sub <- FindVariableFeatures(
  epi_sub,
  assay            = spatial_assay,
  selection.method = "vst",
  nfeatures        = 2000,
  verbose          = FALSE
)

message("Running ScaleData on variable features only ...")
epi_sub <- ScaleData(
  epi_sub,
  features = VariableFeatures(epi_sub),
  verbose  = FALSE
)

message("Running PCA ...")
epi_sub <- RunPCA(
  epi_sub,
  assay   = spatial_assay,
  npcs    = 30,
  verbose = FALSE
)

message("Running FindNeighbors ...")
epi_sub <- FindNeighbors(
  epi_sub,
  dims    = 1:15,
  verbose = FALSE
)

message("Running FindClusters ...")
epi_sub <- FindClusters(
  epi_sub,
  resolution = 0.6,
  verbose    = FALSE
)

message("Running UMAP ...")
epi_sub <- RunUMAP(
  epi_sub,
  dims        = 1:5,
  n.neighbors = 70,
  min.dist    = 0.3,
  seed.use    = 42,
  verbose     = FALSE
)

message("Cluster summary:")
print(table(epi_sub$seurat_clusters, useNA = "always"))

message("Total clusters: ",
        length(unique(as.character(epi_sub$seurat_clusters))))

# ── 2. Raw cluster UMAP / spatial plots ───────────────────────

reclust_order <- levels(Idents(epi_sub))
recluster_cols <- setNames(
  custom_colors[seq_len(length(reclust_order))],
  reclust_order
)

p_umap_raw <- DimPlot(
  object     = epi_sub,
  reduction  = "umap",
  cols       = recluster_cols,
  label      = TRUE,
  label.size = 5,
  repel      = TRUE,
  pt.size    = 0.8
) +
  labs(
    title = "A1281 – Final epithelial reclustering (UMAP)",
    x     = "UMAP 1",
    y     = "UMAP 2"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text       = element_text(size = 12),
    axis.title      = element_text(size = 13),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

save_tiff(
  p_umap_raw,
  file.path(reclust_plot_dir, "04_final_recluster_seurat_clusters_UMAP.tiff"),
  width  = 12,
  height = 10
)

p_spatial_raw <- SpatialDimPlot(
  object         = epi_sub,
  group.by       = "seurat_clusters",
  label          = FALSE,
  pt.size.factor = 2.8
) +
  scale_fill_manual(values = recluster_cols) +
  labs(title = "A1281 – Final epithelial reclustering spatial") +
  theme(
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

save_tiff(
  p_spatial_raw,
  file.path(reclust_plot_dir, "04_final_recluster_seurat_clusters_Spatial.tiff"),
  width  = 12,
  height = 10
)

# ── 3. Final Epi/PIN1–PIN6 annotation ────────────────────────

message("Applying final Epi/PIN annotation ...")

seurat_char <- as.character(epi_sub$seurat_clusters)

reannotate_labels <- dplyr::case_when(
  seurat_char == "0"                ~ "Epi",
  seurat_char %in% c("6", "10")     ~ "PIN1",
  seurat_char %in% c("1", "3")      ~ "PIN2",
  seurat_char == "7"                ~ "PIN3",
  seurat_char %in% c("2", "4", "8") ~ "PIN4",
  seurat_char == "9"                ~ "PIN5",
  seurat_char == "5"                ~ "PIN6",
  TRUE                              ~ NA_character_
)

epi_sub$reannotated <- reannotate_labels

reannotate_order <- c(
  "Epi",
  "PIN1", "PIN2", "PIN3",
  "PIN4", "PIN5", "PIN6"
)

epi_sub$reannotated <- factor(
  as.character(epi_sub$reannotated),
  levels = reannotate_order
)

message("Annotation summary before removing NA:")
print(table(epi_sub$reannotated, useNA = "always"))

message("Cross-table seurat_clusters vs reannotated:")
print(table(epi_sub$seurat_clusters, epi_sub$reannotated, useNA = "always"))

cells_keep <- colnames(epi_sub)[!is.na(epi_sub$reannotated)]
epi_sub <- subset(epi_sub, cells = cells_keep)

message("Cells retained after final annotation: ", ncol(epi_sub))

Idents(epi_sub) <- "reannotated"

reannotate_colors <- c(
  "Epi"  = "red",
  "PIN1" = "green",
  "PIN2" = "blue",
  "PIN3" = "yellow",
  "PIN4" = "purple",
  "PIN5" = "orange",
  "PIN6" = "cyan"
)

# ── 4. Final annotated UMAP / spatial plots ──────────────────

p_umap_annot <- DimPlot(
  object     = epi_sub,
  reduction  = "umap",
  group.by   = "reannotated",
  cols       = reannotate_colors,
  label      = TRUE,
  label.size = 5,
  repel      = TRUE,
  pt.size    = 0.8
) +
  labs(
    title = "A1281 – Final Epi/PIN annotation (UMAP)",
    x     = "UMAP 1",
    y     = "UMAP 2"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text       = element_text(size = 12),
    axis.title      = element_text(size = 13),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

save_tiff(
  p_umap_annot,
  file.path(annot_plot_dir, "04_final_Epi_PIN_UMAP.tiff"),
  width  = 12,
  height = 10
)

p_spatial_annot <- SpatialDimPlot(
  object         = epi_sub,
  group.by       = "reannotated",
  label          = FALSE,
  pt.size.factor = 2.8
) +
  scale_fill_manual(values = reannotate_colors) +
  labs(title = "A1281 – Final Epi/PIN annotation spatial") +
  theme(
    plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text     = element_text(size = 11),
    legend.title    = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

save_tiff(
  p_spatial_annot,
  file.path(annot_plot_dir, "04_final_Epi_PIN_Spatial.tiff"),
  width  = 12,
  height = 10
)

p_combined <- p_umap_annot | p_spatial_annot

save_tiff(
  p_combined,
  file.path(annot_plot_dir, "04_final_Epi_PIN_UMAP_Spatial_combined.tiff"),
  width  = 24,
  height = 10
)

# ── 5. Save tables and RDS ───────────────────────────────────

write.csv(
  as.data.frame(table(epi_sub$seurat_clusters)),
  file.path(table_dir, "04_final_seurat_cluster_counts.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(table(epi_sub$reannotated)),
  file.path(table_dir, "04_final_Epi_PIN_counts.csv"),
  row.names = FALSE
)

saveRDS(
  reannotate_order,
  file.path(rds_dir, "reannotate_order.rds")
)

saveRDS(
  reannotate_colors,
  file.path(rds_dir, "reannotate_colors.rds")
)

saveRDS(
  epi_sub,
  file.path(rds_dir, "04_epi_PIN_reannotated.rds")
)

message("Saved: outputs/rds/04_epi_PIN_reannotated.rds")
message("Module 4 complete.")
