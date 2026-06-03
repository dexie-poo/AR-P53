#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

obj <- readRDS(file.path(rds_dir, "04_gsea_complete.rds"))
raw_assay <- readRDS(file.path(rds_dir, "raw_assay.rds"))
reannotate_order <- readRDS(file.path(rds_dir, "reannotate_order.rds"))
reannotate_colors <- readRDS(file.path(rds_dir, "reannotate_colors.rds"))

genes_AR  <- c("Nkx3-1", "Azgp1", "Fkbp5", "Tmprss2", "Sgk1")
genes_P53 <- c("Cdkn1a", "Bbc3", "Pmaip1", "Gadd45a", "Bax", "Fas", "Rrm2b", "Sesn1", "Sesn2", "Noxa", "Puma", "Xpc", "Wig1", "Perp")
genes_MET <- c("Snai1", "Cdh2", "Cxcr4", "Hras", "Hgf", "Hpn", "Uba52", "Stat3", "Akt1", "Raf1", "Fak", "Foxm1", "Src")
genes_WNT <- c("Ctnnb1", "Cd44", "Dkk2", "Lgr4", "Ccnd1", "Sp5", "Dkk1", "Myc")
genes_ordered <- c(genes_AR, genes_P53, genes_MET, genes_WNT)

module_map <- data.frame(
  gene = genes_ordered,
  module = c(
    rep("AR Signalling", length(genes_AR)),
    rep("P53 Pathway", length(genes_P53)),
    rep("MET Signalling", length(genes_MET)),
    rep("WNT Signalling", length(genes_WNT))
  ),
  stringsAsFactors = FALSE
)

existing_scaled <- tryCatch(
  rownames(GetAssayData(obj, assay = raw_assay, layer = "scale.data")),
  error = function(e) rownames(GetAssayData(obj, assay = raw_assay, slot = "scale.data"))
)

genes_found <- genes_ordered[genes_ordered %in% rownames(obj)]
genes_to_rescale <- genes_found[!genes_found %in% existing_scaled]
if (length(genes_to_rescale) > 0) {
  obj <- ScaleData(obj, assay = raw_assay, features = unique(c(existing_scaled, genes_to_rescale)), verbose = FALSE)
}

scale_mat <- tryCatch(
  as.matrix(GetAssayData(obj, assay = raw_assay, layer = "scale.data")[genes_found, , drop = FALSE]),
  error = function(e) as.matrix(GetAssayData(obj, assay = raw_assay, slot = "scale.data")[genes_found, , drop = FALSE])
)

clusters_vec <- as.character(obj$reannotated)
avg_mat <- sapply(reannotate_order, function(cl) {
  cells <- colnames(obj)[clusters_vec == cl]
  if (length(cells) == 0) return(rep(NA, nrow(scale_mat)))
  rowMeans(scale_mat[, cells, drop = FALSE], na.rm = TRUE)
})

rownames(avg_mat) <- genes_found
colnames(avg_mat) <- reannotate_order

clamp_val <- 0.3
avg_mat_clamped <- pmin(pmax(avg_mat, -clamp_val), clamp_val)

module_found <- module_map[module_map$gene %in% genes_found, ]
module_found <- module_found[match(genes_found, module_found$gene), ]
row_anno <- data.frame(Module = module_found$module, row.names = genes_found)

module_colors <- c(
  "AR Signalling" = "#FF1493",
  "P53 Pathway" = "#E41A1C",
  "MET Signalling" = "#377EB8",
  "WNT Signalling" = "#4DAF4A"
)
anno_colors <- list(Module = module_colors)

n_colors <- 100
col_breaks <- seq(-clamp_val, clamp_val, length.out = n_colors + 1)
pal <- colorRampPalette(c("#5B0082", "#000000", "#FFD700"))(n_colors)

p_ph <- pheatmap(
  avg_mat_clamped,
  color = pal,
  breaks = col_breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno,
  annotation_colors = anno_colors,
  fontsize_row = 17,
  fontsize_col = 17,
  fontsize = 17,
  border_color = NA,
  main = "A1281 – AR / P53 / MET / WNT Modules",
  angle_col = 90,
  legend_breaks = c(-clamp_val, -clamp_val / 2, 0, clamp_val / 2, clamp_val),
  legend_labels = c(as.character(-clamp_val), as.character(-clamp_val / 2), "0", as.character(clamp_val / 2), as.character(clamp_val)),
  silent = TRUE
)

genes_vln <- c("Krt5", "Krt8", "Krt19", "Cdh1", "Cd44", "Ctnnb1")
genes_found_vln <- genes_vln[genes_vln %in% rownames(obj)]

expr_mat <- tryCatch(
  as.matrix(GetAssayData(obj, assay = raw_assay, layer = "data")[genes_found_vln, , drop = FALSE]),
  error = function(e) as.matrix(GetAssayData(obj, assay = raw_assay, slot = "data")[genes_found_vln, , drop = FALSE])
)

df_long <- as.data.frame(t(expr_mat)) %>%
  rownames_to_column("barcode") %>%
  mutate(cell_type = as.character(obj$reannotated)) %>%
  pivot_longer(cols = all_of(genes_found_vln), names_to = "gene", values_to = "expression") %>%
  mutate(
    cell_type = factor(cell_type, levels = reannotate_order),
    gene = factor(gene, levels = genes_found_vln)
  ) %>%
  dplyr::filter(!is.na(cell_type))

stat_results <- df_long %>%
  group_by(gene) %>%
  kruskal_test(expression ~ cell_type) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

write.csv(stat_results, file.path(out_dir, "A1281_reannotated_violin_KruskalWallis.csv"), row.names = FALSE)

dunn_results <- df_long %>%
  group_by(gene) %>%
  dunn_test(expression ~ cell_type, p.adjust.method = "BH") %>%
  add_significance("p.adj")

write.csv(dunn_results, file.path(out_dir, "A1281_reannotated_violin_Dunn_pairwise.csv"), row.names = FALSE)

genes_plot <- c("Krt5", "Krt8", "Pbsn", "Hoxb13")
genes_found_plot <- genes_plot[genes_plot %in% rownames(obj)]

p_umap_reannotated <- DimPlot(
  object = obj,
  reduction = "umap",
  group.by = "reannotated",
  cols = reannotate_colors,
  label = FALSE,
  pt.size = 0.8
) +
  labs(title = "A1281 – Reannotated Clusters (UMAP)", x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

p_spatial_reannotated <- SpatialDimPlot(
  object = obj,
  group.by = "reannotated",
  label = FALSE,
  pt.size.factor = 2.8
) +
  scale_fill_manual(values = reannotate_colors) +
  labs(title = "A1281 – Reannotated Clusters (Spatial)") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.size = unit(0.6, "cm")
  )

save_tiff(p_umap_reannotated, file.path(plot_dir, "A1281_reannotated_UMAP.tiff"), width = 12, height = 10)
save_tiff(p_spatial_reannotated, file.path(plot_dir, "A1281_reannotated_SpatialDimPlot.tiff"), width = 12, height = 10)
save_tiff(p_umap_reannotated | p_spatial_reannotated, file.path(plot_dir, "A1281_reannotated_UMAP_Spatial_combined.tiff"), width = 24, height = 10)

tiff(file.path(plot_dir, "A1281_reannotated_AR_P53_MET_WNT_pheatmap.tiff"),
     width = 7, height = 8, units = "in", res = 300, compression = "lzw")
grid::grid.newpage()
grid::grid.draw(p_ph$gtable)
dev.off()

plot_list <- list()
for (g in genes_found_vln) {
  df_gene <- df_long %>% dplyr::filter(gene == !!g)
  y_max <- max(df_gene$expression, na.rm = TRUE)
  y_pos <- y_max * 1.05
  if (!is.finite(y_pos) || y_pos == 0) y_pos <- 1

  stat_gene <- dunn_results %>%
    dplyr::filter(gene == !!g, p.adj < 0.05) %>%
    dplyr::select(group1, group2, p.adj, p.adj.signif) %>%
    as.data.frame()

  p <- ggplot(df_gene, aes(x = cell_type, y = expression, fill = cell_type)) +
    geom_violin(trim = FALSE, alpha = 0.8, scale = "width") +
    geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = reannotate_colors) +
    labs(title = g, x = NULL, y = "Expression") +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold.italic", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "none"
    )

  if (nrow(stat_gene) > 0) {
    p <- p + stat_pvalue_manual(
      stat_gene,
      label = "p.adj.signif",
      y.position = y_pos,
      step.increase = 0.09,
      tip.length = 0.01,
      size = 3
    )
  }

  save_tiff(p, file.path(vln_dir, paste0("A1281_reannotated_violin_", g, ".tiff")), width = 8, height = 7)
  plot_list[[g]] <- p
}

n_genes <- length(plot_list)
n_cols <- 3
n_rows <- ceiling(n_genes / n_cols)

p_combined_vln <- wrap_plots(plot_list, ncol = n_cols) +
  plot_annotation(
    title = "A1281 – Gene Expression by Cell Type (Reannotated)",
    subtitle = "Dunn post-hoc (BH): * p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

save_tiff(
  p_combined_vln,
  file.path(out_dir, "A1281_reannotated_violin_combined_panel.tiff"),
  width = 24,
  height = n_rows * 8
)

if (length(plot_list) == 0) {
  stop("plot_list is empty; no violin plots were generated.")
}

genes_found_plot <- genes_plot[genes_plot %in% rownames(obj)]

if (length(genes_found_plot) > 0) {

  p_feature_umap <- FeaturePlot(
    object = obj,
    reduction = "umap",
    features = genes_found_plot,
    ncol = 2,
    pt.size = 0.8,
    order = TRUE,
    cols = c("lightgrey", "red")
  ) &
    theme_classic() &
    theme(
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.5, "cm"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    )

  save_tiff(
    p_feature_umap,
    file.path(plot_dir, "A1281_epi_FeaturePlot_UMAP_Krt5_Krt8_Pbsn_Hoxb13.tiff"),
    width = 14,
    height = 12
  )

  p_feature_spatial <- SpatialFeaturePlot(
    object = obj,
    features = genes_found_plot,
    ncol = 2,
    pt.size.factor = 2.8,
    alpha = 0.8
  ) &
    theme(
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.5, "cm")
    )

  save_tiff(
    p_feature_spatial,
    file.path(plot_dir, "A1281_epi_SpatialFeaturePlot_Krt5_Krt8_Pbsn_Hoxb13.tiff"),
    width = 14,
    height = 12
  )

  for (g in genes_found_plot) {
    p_umap_gene <- FeaturePlot(
      object = obj,
      reduction = "umap",
      features = g,
      pt.size = 0.8,
      order = TRUE,
      cols = c("lightgrey", "red")
    ) +
      theme_classic() +
      labs(title = g) +
      theme(
        plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12)
      )

    save_tiff(
      p_umap_gene,
      file.path(expr_dir, paste0("A1281_epi_UMAP_", g, ".tiff")),
      width = 8,
      height = 7
    )

    p_spatial_gene <- SpatialFeaturePlot(
      object = obj,
      features = g,
      pt.size.factor = 2.8,
      alpha = 0.8
    ) +
      labs(title = g) +
      theme(
        plot.title = element_text(size = 16, face = "bold.italic", hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")
      )

    save_tiff(
      p_spatial_gene,
      file.path(expr_dir, paste0("A1281_epi_Spatial_", g, ".tiff")),
      width = 8,
      height = 7
    )

    p_combined_gene <- p_umap_gene | p_spatial_gene

    save_tiff(
      p_combined_gene,
      file.path(expr_dir, paste0("A1281_epi_combined_", g, ".tiff")),
      width = 18,
      height = 7
    )
  }

  p_full_combined <- p_feature_umap | p_feature_spatial

  save_tiff(
    p_full_combined,
    file.path(plot_dir, "A1281_epi_FeaturePlot_UMAP_Spatial_combined.tiff"),
    width = 28,
    height = 12
  )
}

saveRDS(obj, file.path(rds_dir, "05_visualization_complete.rds"))
