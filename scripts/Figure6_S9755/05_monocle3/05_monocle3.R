#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 module 05
# Monocle3 trajectory on final Epi/PIN object; root selected from Epi2 centroid
# =============================================================
source("scripts/06_S9755/00_config.R")
suppressPackageStartupMessages({ library(monocle3); library(SingleCellExperiment) })

obj <- readRDS(file.path(RDS_DIR, "03_S9755_EpiPIN_subset_annotated_v2.rds"))
assay <- get_spatial_assay(obj)
DefaultAssay(obj) <- assay
Idents(obj) <- "epi_pin_cell_type_v2"

counts <- get_counts(obj, assay)
cell_metadata <- obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
reducedDims(cds)$PCA <- Embeddings(obj, "pca")[colnames(cds), , drop = FALSE]
reducedDims(cds)$UMAP <- Embeddings(obj, "umap")[colnames(cds), , drop = FALSE]

message("Clustering and learning graph in Monocle3...")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE, learn_graph_control = list(minimal_branch_len = 10))

message("Selecting root node from Epi2 UMAP centroid...")
umap_coords <- reducedDims(cds)$UMAP
epi2_cells <- rownames(colData(cds))[as.character(colData(cds)$epi_pin_cell_type_v2) == "Epi2"]
if (length(epi2_cells) == 0) stop("No Epi2 cells found for Monocle root selection.")
epi2_centroid <- colMeans(umap_coords[epi2_cells, , drop = FALSE])
principal_nodes <- cds@principal_graph_aux[["UMAP"]]$dp_mst
node_coords <- t(principal_nodes)
dist_to_centroid <- apply(node_coords, 1, function(x) sqrt(sum((x - epi2_centroid)^2)))
closest_node <- names(which.min(dist_to_centroid))
message("Root node: ", closest_node)
cds <- order_cells(cds, root_pr_nodes = closest_node)

pt <- pseudotime(cds)
pt_df <- data.frame(
  barcode = names(pt),
  pseudotime = as.numeric(pt),
  epi_pin_cell_type_v2 = as.character(colData(cds)$epi_pin_cell_type_v2),
  monocle_cluster = as.character(clusters(cds)),
  monocle_partition = as.character(partitions(cds))
)
write.csv(pt_df, file.path(TABLE_DIR, "05_S9755_EpiPIN_pseudotime.csv"), row.names = FALSE)

obj$monocle3_pseudotime <- pt[colnames(obj)]
obj$monocle3_partition <- as.character(partitions(cds)[colnames(obj)])
obj$monocle3_cluster <- as.character(clusters(cds)[colnames(obj)])
obj$monocle3_pseudotime[!is.finite(obj$monocle3_pseudotime)] <- NA

saveRDS(cds, file.path(RDS_DIR, "05_S9755_EpiPIN_monocle3_cds_Epi2_root.rds"))
saveRDS(obj, file.path(RDS_DIR, "05_S9755_EpiPIN_subset_annotated_v2_monocle.rds"))

message("Plotting Monocle outputs...")
save_tiff(plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE) + ggtitle("S9755 Epi/PIN pseudotime"),
          file.path(PLOT_DIR, "05_S9755_Monocle_UMAP_pseudotime.tiff"), width = 8, height = 7)
save_tiff(plot_cells(cds, color_cells_by = "epi_pin_cell_type_v2", label_groups_by_cluster = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE) + ggtitle("S9755 Epi/PIN trajectory by cell type"),
          file.path(PLOT_DIR, "05_S9755_Monocle_UMAP_cell_type.tiff"), width = 8, height = 7)
save_tiff(VlnPlot(obj, features = "monocle3_pseudotime", group.by = "epi_pin_cell_type_v2", pt.size = 0) +
            ggtitle("S9755 pseudotime by final Epi/PIN state"),
          file.path(PLOT_DIR, "05_S9755_Monocle_violin_pseudotime_by_cell_type.tiff"), width = 8, height = 6)
save_tiff(SpatialFeaturePlot(obj, features = "monocle3_pseudotime", pt.size.factor = 2.5) + ggtitle("S9755 spatial pseudotime"),
          file.path(PLOT_DIR, "05_S9755_Monocle_spatial_pseudotime.tiff"), width = 8, height = 8)

message("Graph test...")
gt <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
write.csv(as.data.frame(gt), file.path(TABLE_DIR, "05_S9755_Monocle_graph_test.csv"))
saveRDS(gt, file.path(RDS_DIR, "05_S9755_Monocle_graph_test.rds"))
