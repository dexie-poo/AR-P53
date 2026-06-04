#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 module 04
# Marker panel heatmaps, violin plots, UMAP FeaturePlots, spatial expression plots
# =============================================================
source("scripts/06_S9755/00_config.R")
suppressPackageStartupMessages({ library(pheatmap); library(grid) })

obj <- readRDS(file.path(RDS_DIR, "03_S9755_EpiPIN_subset_annotated_v2.rds"))
assay <- get_spatial_assay(obj)
DefaultAssay(obj) <- assay
Idents(obj) <- "epi_pin_cell_type_v2"
obj$epi_pin_cell_type_v2 <- factor(as.character(obj$epi_pin_cell_type_v2), levels = EPIPIN_LEVELS)

panel_genes <- list(
  AR  = c("Nkx3-1", "Azgp1", "Fkbp5", "Tmprss2", "Sgk1"),
  P53 = c("Cdkn1a", "Bbc3", "Pmaip1", "Gadd45a", "Bax", "Fas", "Rrm2b", "Sesn1", "Sesn2", "Noxa", "Puma", "Xpc", "Wig1", "Perp"),
  MET = c("Snai1", "Cdh2", "Cxcr4", "Hras", "Hgf", "Hpn", "Uba52", "Stat3", "Akt1", "Raf1", "Fak", "Foxm1", "Src"),
  WNT = c("Ctnnb1", "Cd44", "Dkk2", "Lgr4", "Ccnd1", "Sp5", "Dkk1", "Myc"),
  KRAS = c("Jun", "Dusp6", "Bcl2l1", "Mcl1", "Cdk4")
)
all_genes <- unique(unlist(panel_genes))
genes_found <- all_genes[all_genes %in% rownames(obj)]
genes_missing <- setdiff(all_genes, genes_found)
write.csv(data.frame(gene = all_genes, found = all_genes %in% genes_found),
          file.path(TABLE_DIR, "04_S9755_panel_gene_presence.csv"), row.names = FALSE)
if (length(genes_missing) > 0) message("Missing genes: ", paste(genes_missing, collapse = ", "))

message("Scaling panel genes if needed...")
obj <- ScaleData(obj, assay = assay, features = unique(c(VariableFeatures(obj), genes_found)), verbose = FALSE)

message("Building average-expression pheatmap...")
scale_mat <- as.matrix(tryCatch(GetAssayData(obj, assay = assay, layer = "scale.data"),
                                error = function(e) GetAssayData(obj, assay = assay, slot = "scale.data")))
scale_mat <- scale_mat[genes_found, , drop = FALSE]
clusters_vec <- as.character(obj$epi_pin_cell_type_v2)
avg_mat <- sapply(EPIPIN_LEVELS, function(cl) {
  cells <- colnames(obj)[clusters_vec == cl]
  if (length(cells) == 0) return(rep(NA_real_, nrow(scale_mat)))
  rowMeans(scale_mat[, cells, drop = FALSE], na.rm = TRUE)
})
rownames(avg_mat) <- genes_found
colnames(avg_mat) <- EPIPIN_LEVELS
clamp_val <- 0.3
avg_mat <- pmin(pmax(avg_mat, -clamp_val), clamp_val)

module_df <- bind_rows(lapply(names(panel_genes), function(module) {
  data.frame(gene = panel_genes[[module]], Module = module)
})) %>% filter(gene %in% genes_found)
module_df <- module_df[match(genes_found, module_df$gene), , drop = FALSE]
row_anno <- data.frame(Module = module_df$Module, row.names = module_df$gene)

ph <- pheatmap(avg_mat, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = row_anno,
               color = colorRampPalette(c("#5B0082", "black", "#FFD700"))(100),
               breaks = seq(-clamp_val, clamp_val, length.out = 101), fontsize_row = 11,
               fontsize_col = 12, border_color = NA,
               main = "S9755 Epi/PIN pathway panel genes", silent = TRUE)
tiff(file.path(PLOT_DIR, "04_S9755_EpiPIN_panel_gene_pheatmap.tiff"), width = 8, height = 10, units = "in", res = 300, compression = "lzw")
grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()

message("Violin plots...")
for (module in names(panel_genes)) {
  genes <- panel_genes[[module]][panel_genes[[module]] %in% rownames(obj)]
  if (length(genes) == 0) next
  p <- VlnPlot(obj, features = genes, group.by = "epi_pin_cell_type_v2", pt.size = 0, ncol = min(5, length(genes))) +
    plot_annotation(title = paste0("S9755 ", module, " genes"))
  save_tiff(p, file.path(PLOT_DIR, paste0("04_S9755_VlnPlot_", module, ".tiff")), width = 14, height = 8)
}

expr_dir <- file.path(PLOT_DIR, "04_expression_plots")
dir.create(expr_dir, recursive = TRUE, showWarnings = FALSE)
message("UMAP FeaturePlots and SpatialFeaturePlots...")
for (gene in genes_found) {
  save_tiff(FeaturePlot(obj, features = gene, reduction = "umap", pt.size = 0.4) + ggtitle(gene),
            file.path(expr_dir, paste0("UMAP_", gene, ".tiff")), width = 6, height = 5)
  save_tiff(SpatialFeaturePlot(obj, features = gene, pt.size.factor = 2.5) + ggtitle(gene),
            file.path(expr_dir, paste0("Spatial_", gene, ".tiff")), width = 7, height = 7)
}

saveRDS(obj, file.path(RDS_DIR, "04_S9755_EpiPIN_subset_annotated_v2_panel_scaled.rds"))
