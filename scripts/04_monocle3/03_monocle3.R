#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

obj <- readRDS(file.path(rds_dir, "02_reannotated.rds"))
raw_assay <- readRDS(file.path(rds_dir, "raw_assay.rds"))

count_matrix <- tryCatch(
  GetAssayData(obj, assay = raw_assay, layer = "counts"),
  error = function(e) GetAssayData(obj, assay = raw_assay, slot = "counts")
)

cell_metadata <- obj@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(count_matrix),
  row.names = rownames(count_matrix)
)

cds <- new_cell_data_set(
  expression_data = count_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

reducedDims(cds)[["UMAP"]] <- Embeddings(obj, reduction = "umap")
reducedDims(cds)[["PCA"]] <- Embeddings(obj, reduction = "pca")

cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3, verbose = FALSE)
cds <- learn_graph(
  cds,
  use_partition = TRUE,
  learn_graph_control = list(minimal_branch_len = 10, orthogonal_proj_tip = FALSE),
  verbose = FALSE
)

umap_coords <- reducedDims(cds)[["UMAP"]]
epi_cells <- colnames(cds)[colData(cds)$reannotated == "Epi"]
epi_centroid <- colMeans(umap_coords[epi_cells, , drop = FALSE])

graph_nodes <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
dist_to_centroid <- sqrt((graph_nodes[, 1] - epi_centroid[1])^2 + (graph_nodes[, 2] - epi_centroid[2])^2)
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

write.csv(pt_df, file.path(mono_dir, "Tables", "pseudotime_reannotated.csv"), row.names = FALSE)


pt_named <- setNames(pseudotime(cds), colnames(cds))
obj$monocle3_pseudotime <- pt_named[colnames(obj)]
saveRDS(cds, file.path(rds_dir, "03_monocle_cds.rds"))
saveRDS(obj, file.path(rds_dir, "03_monocle_pseudotime_seurat.rds"))
