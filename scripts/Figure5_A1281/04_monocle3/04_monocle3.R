#!/usr/bin/env Rscript

source("scripts/Figure5_A1281/00_setup/00_config.R")

epi_sub <- readRDS(file.path(rds_dir, "04_epi_PIN_reannotated.rds"))

DefaultAssay(epi_sub) <- spatial_assay
Idents(epi_sub) <- "reannotated"

mono_rds_dir <- file.path(mono_dir, "RDS")
mono_plot_dir <- file.path(mono_dir, "Plots")
mono_table_dir <- file.path(mono_dir, "Tables")

dir.create(mono_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mono_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mono_table_dir, recursive = TRUE, showWarnings = FALSE)

count_matrix <- tryCatch(
  GetAssayData(epi_sub, assay = spatial_assay, layer = "counts"),
  error = function(e) GetAssayData(epi_sub, assay = spatial_assay, slot = "counts")
)

cell_metadata <- epi_sub@meta.data

gene_metadata <- data.frame(
  gene_short_name = rownames(count_matrix),
  row.names = rownames(count_matrix)
)

cds <- new_cell_data_set(
  expression_data = count_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

reducedDims(cds)[["UMAP"]] <- Embeddings(epi_sub, reduction = "umap")
reducedDims(cds)[["PCA"]] <- Embeddings(epi_sub, reduction = "pca")

cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",
  resolution = 1e-3,
  verbose = FALSE
)

cds <- learn_graph(
  cds,
  use_partition = TRUE,
  learn_graph_control = list(
    minimal_branch_len = 10,
    orthogonal_proj_tip = FALSE
  ),
  verbose = FALSE
)

umap_coords <- reducedDims(cds)[["UMAP"]]
epi_cells <- colnames(cds)[colData(cds)$reannotated == "Epi"]
epi_centroid <- colMeans(umap_coords[epi_cells, , drop = FALSE])

graph_nodes <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)

dist_to_centroid <- sqrt(
  (graph_nodes[, 1] - epi_centroid[1])^2 +
    (graph_nodes[, 2] - epi_centroid[2])^2
)

closest_node <- names(which.min(dist_to_centroid))

cds <- order_cells(cds, root_pr_nodes = closest_node)

pt_df <- data.frame(
  cell_barcode = colnames(cds),
  pseudotime = pseudotime(cds),
  reannotated = as.character(colData(cds)$reannotated),
  monocle_cluster = as.character(clusters(cds)),
  monocle_partition = as.character(partitions(cds)),
  stringsAsFactors = FALSE
)

write.csv(
  pt_df,
  file.path(mono_table_dir, "05_pseudotime_reannotated.csv"),
  row.names = FALSE
)

saveRDS(
  cds,
  file.path(mono_rds_dir, "05_monocle3_cds_A1281_Epi_root.rds")
)

epi_sub$monocle3_pseudotime <- setNames(pseudotime(cds), colnames(cds))[colnames(epi_sub)]

reannotate_colors <- readRDS(file.path(rds_dir, "reannotate_colors.rds"))

p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  label_roots = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 1,
  graph_label_size = 3,
  cell_size = 0.8
) +
  scale_color_viridis_c(option = "C", name = "Pseudotime") +
  ggtitle("Pseudotime (Monocle3) – A1281 – Root: Epi") +
  theme_classic(base_size = 14)

save_tiff(
  p_pt,
  file.path(mono_plot_dir, "05_UMAP_pseudotime.tiff"),
  width = 10,
  height = 8
)

p_ct <- plot_cells(
  cds,
  color_cells_by = "reannotated",
  label_cell_groups = TRUE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 1,
  cell_size = 0.8
) +
  scale_color_manual(values = reannotate_colors, name = "Cell Type") +
  ggtitle("Cell Type – A1281 Reannotated") +
  theme_classic(base_size = 14)

save_tiff(
  p_ct,
  file.path(mono_plot_dir, "05_UMAP_CellType_trajectory.tiff"),
  width = 12,
  height = 8
)

p_spatial_pseudotime <- SpatialFeaturePlot(
  object = epi_sub,
  features = "monocle3_pseudotime",
  pt.size.factor = 2.8
) +
  scale_fill_viridis_c(option = "C", name = "Pseudotime", na.value = "lightgrey") +
  labs(title = "A1281 – Pseudotime (Spatial)") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

save_tiff(
  p_spatial_pseudotime,
  file.path(mono_plot_dir, "05_Spatial_pseudotime.tiff"),
  width = 12,
  height = 10
)

saveRDS(
  epi_sub,
  file.path(rds_dir, "05_epi_PIN_monocle_seurat.rds")
)

message("Module 5 complete.")
