#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

raw_obj <- Load10X_Spatial(data.dir = raw_dir, bin.size = 24)
raw_assay <- get_raw_assay(raw_obj)
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
obj <- NormalizeData(obj, assay = raw_assay, verbose = FALSE)
obj <- FindVariableFeatures(obj, assay = raw_assay, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, assay = raw_assay, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, assay = raw_assay, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:15, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:5, n.neighbors = 50, min.dist = 0.3, seed.use = 42, verbose = FALSE)

saveRDS(obj, file.path(out_dir, "01_A1281_loaded_qc_clustered.rds"))
saveRDS(raw_assay, file.path(out_dir, "raw_assay.rds"))