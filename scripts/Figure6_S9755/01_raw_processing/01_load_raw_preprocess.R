#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 module 01
# Load raw 10x VisiumHD 24um data, QC-filter, LogNormalize, cluster
# =============================================================
source("scripts/06_S9755/00_config.R")

message("Loading raw 10x spatial data: ", RAW_10X_DIR)
obj <- Load10X_Spatial(data.dir = RAW_10X_DIR, bin.size = BIN_SIZE)
assay <- get_spatial_assay(obj)
DefaultAssay(obj) <- assay
obj <- JoinLayers(obj, assay = assay)

message("Computing QC metrics...")
counts <- get_counts(obj, assay)
is_mito <- grepl("^MT-|^mt-", rownames(obj))
obj$sum <- Matrix::colSums(counts)
obj$detected <- Matrix::colSums(counts > 0)
obj$subsets_mito_percent <- if (any(is_mito)) {
  Matrix::colSums(counts[is_mito, , drop = FALSE]) / pmax(obj$sum, 1) * 100
} else 0
obj$subsets_mito_percent[is.na(obj$subsets_mito_percent)] <- 0
if (!"cell_count" %in% colnames(obj@meta.data)) obj$cell_count <- 1

obj$qc_lib_size  <- obj$sum < LIBSIZE_MIN
obj$qc_detected  <- obj$detected < MIN_FEATURES
obj$qc_mito_prop <- obj$subsets_mito_percent > MITO_THRESH
obj$global_outliers <- obj$qc_lib_size | obj$qc_detected | obj$qc_mito_prop

qc_summary <- data.frame(
  metric = c("raw_spots", "qc_lib_size", "qc_detected", "qc_mito_prop", "global_outliers", "kept_spots"),
  value  = c(ncol(obj), sum(obj$qc_lib_size), sum(obj$qc_detected), sum(obj$qc_mito_prop),
             sum(obj$global_outliers), sum(!obj$global_outliers))
)
write.csv(qc_summary, file.path(TABLE_DIR, "01_S9755_QC_summary.csv"), row.names = FALSE)
print(qc_summary)

message("Saving QC plots...")
qc_df <- obj@meta.data
qc_df$spot <- rownames(qc_df)
make_qc_violin <- function(feature, threshold, flag, label, log_scale = FALSE) {
  p <- ggplot(qc_df, aes(x = "All spots", y = .data[[feature]])) +
    geom_violin(fill = "steelblue", alpha = 0.6, color = NA, scale = "width") +
    geom_jitter(aes(color = .data[[flag]]), size = 0.2, alpha = 0.3, width = 0.25) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"), name = NULL) +
    theme_classic() + labs(x = NULL, y = label, title = label) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
  if (log_scale) p <- p + scale_y_continuous(trans = "pseudo_log")
  p
}
save_tiff(
  make_qc_violin("sum", LIBSIZE_MIN, "qc_lib_size", "Total UMI counts", TRUE) |
    make_qc_violin("detected", MIN_FEATURES, "qc_detected", "Detected genes", TRUE) |
    make_qc_violin("subsets_mito_percent", MITO_THRESH, "qc_mito_prop", "Mitochondrial %", FALSE),
  file.path(PLOT_DIR, "01_S9755_QC_violins_pre_filter.tiff"), width = 14, height = 6
)

message("Filtering spots and zero-count genes...")
obj <- subset(obj, cells = colnames(obj)[!obj$global_outliers])
counts_post <- get_counts(obj, assay)
obj <- subset(obj, features = rownames(obj)[Matrix::rowSums(counts_post) > 0])

message("Post-QC dimensions: ", paste(dim(obj), collapse = " x "))

message("LogNormalize...")
obj <- NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

message("Gene prefiltering: pct >= ", GENE_PCT_MIN, "%, mean >= ", GENE_MEAN_MIN)
counts <- get_counts(obj, assay)
data <- get_data(obj, assay)
pct_expr <- Matrix::rowSums(counts > 0) / ncol(counts) * 100
mean_expr <- Matrix::rowMeans(data)
genes_pass <- names(pct_expr[pct_expr >= GENE_PCT_MIN & mean_expr >= GENE_MEAN_MIN])
obj <- obj[genes_pass, ]

message("HVG / ScaleData / PCA / clustering / UMAP...")
obj <- FindVariableFeatures(obj, assay = assay, selection.method = "vst", nfeatures = N_HVG_RAW, verbose = FALSE)
obj <- ScaleData(obj, assay = assay, verbose = FALSE)
obj <- RunPCA(obj, assay = assay, npcs = RAW_NPCS, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = RAW_CLUSTER_DIMS, verbose = FALSE)
for (res in RAW_CLUSTER_RESES) obj <- FindClusters(obj, resolution = res, verbose = FALSE)
obj <- FindClusters(obj, resolution = RAW_FINAL_RES, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = RAW_CLUSTER_DIMS, verbose = FALSE)

save_tiff(plot_basic_dim(obj, "seurat_clusters", title = "S9755 24um raw processed clusters"),
          file.path(PLOT_DIR, "01_S9755_UMAP_raw_processed.tiff"), width = 8, height = 7)
saveRDS(obj, file.path(RDS_DIR, "01_S9755_24um_lognorm_processed.rds"))
message("Saved: ", file.path(RDS_DIR, "01_S9755_24um_lognorm_processed.rds"))
