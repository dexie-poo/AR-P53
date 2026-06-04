#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

raw_obj <- Load10X_Spatial(
  data.dir = raw_dir,
  bin.size = 24
)

raw_assay <- get_spatial_assay(raw_obj)
DefaultAssay(raw_obj) <- raw_assay

raw_obj <- JoinLayers(raw_obj, assay = raw_assay)

counts_mat <- tryCatch(
  GetAssayData(raw_obj, assay = raw_assay, layer = "counts"),
  error = function(e) GetAssayData(raw_obj, assay = raw_assay, slot = "counts")
)

is_mito <- grepl("^mt-|^MT-", rownames(raw_obj))

raw_obj$sum <- Matrix::colSums(counts_mat)
raw_obj$detected <- Matrix::colSums(counts_mat > 0)

raw_obj$subsets_mito_percent <- ifelse(
  any(is_mito),
  Matrix::colSums(counts_mat[is_mito, , drop = FALSE]) / pmax(raw_obj$sum, 1) * 100,
  0
)

raw_obj$subsets_mito_percent[is.na(raw_obj$subsets_mito_percent)] <- 0

raw_obj$qc_lib_size  <- raw_obj$sum < LIBSIZE_MIN
raw_obj$qc_detected  <- raw_obj$detected < MIN_FEATURES
raw_obj$qc_mito_prop <- raw_obj$subsets_mito_percent > MITO_THRESH

raw_obj$global_outliers <- raw_obj$qc_lib_size | raw_obj$qc_detected | raw_obj$qc_mito_prop

obj <- subset(raw_obj, cells = colnames(raw_obj)[!raw_obj$global_outliers])

counts_post <- tryCatch(
  GetAssayData(obj, assay = raw_assay, layer = "counts"),
  error = function(e) GetAssayData(obj, assay = raw_assay, slot = "counts")
)

obj <- subset(obj, features = rownames(obj)[Matrix::rowSums(counts_post) > 0])

DefaultAssay(obj) <- raw_assay

obj <- NormalizeData(
  obj,
  assay = raw_assay,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

counts_mat <- tryCatch(
  GetAssayData(obj, assay = raw_assay, layer = "counts"),
  error = function(e) GetAssayData(obj, assay = raw_assay, slot = "counts")
)

data_mat <- tryCatch(
  GetAssayData(obj, assay = raw_assay, layer = "data"),
  error = function(e) GetAssayData(obj, assay = raw_assay, slot = "data")
)

pct_expr <- Matrix::rowSums(counts_mat > 0) / ncol(counts_mat) * 100
mean_expr <- Matrix::rowMeans(data_mat)

genes_pass <- names(pct_expr[pct_expr >= 1 & mean_expr >= 0.05])

obj <- obj[genes_pass, ]

obj <- FindVariableFeatures(
  obj,
  assay = raw_assay,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

obj <- ScaleData(
  obj,
  assay = raw_assay,
  verbose = FALSE
)

obj <- RunPCA(
  obj,
  assay = raw_assay,
  npcs = 50,
  verbose = FALSE
)

obj <- FindNeighbors(
  obj,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

obj <- FindClusters(
  obj,
  resolution = 0.5,
  verbose = FALSE
)

obj <- RunUMAP(
  obj,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

p_raw_umap <- DimPlot(
  obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.8
) +
  theme_classic() +
  labs(title = "A1281 24um raw processed clustered UMAP")

save_tiff(
  p_raw_umap,
  file.path(raw_plot_dir, "01_A1281_raw_processed_UMAP.tiff"),
  width = 12,
  height = 10
)

p_raw_spatial <- SpatialDimPlot(
  obj,
  label = FALSE,
  pt.size.factor = 2.8
) +
  labs(title = "A1281 24um raw processed spatial clusters")

save_tiff(
  p_raw_spatial,
  file.path(raw_plot_dir, "01_A1281_raw_processed_spatial.tiff"),
  width = 12,
  height = 10
)

saveRDS(obj, file.path(rds_dir, "01_raw_lognorm_processed.rds"))
saveRDS(raw_assay, file.path(rds_dir, "spatial_assay.rds"))

writeLines(
  capture.output(sessionInfo()),
  file.path(log_dir, "01_raw_processing_sessionInfo.txt")
)
