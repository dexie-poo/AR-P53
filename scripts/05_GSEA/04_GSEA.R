#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

obj <- readRDS(file.path(out_dir, "A1281_epi_sub_reannotated_monocle.rds"))
raw_assay <- readRDS(file.path(out_dir, "raw_assay.rds"))
reannotate_order <- readRDS(file.path(out_dir, "reannotate_order.rds"))

load("/gpfs/home/leungd02/RDS_Files/m_t2g.RData")

build_pathways <- function(pattern) {
  sub <- m_t2g[grep(pattern, m_t2g$gs_name), ]
  if (nrow(sub) == 0) return(list())
  pw <- split(sub$entrez_gene, sub$gs_name)
  lapply(pw, unique)
}

gmt_pathways <- list(
  Hallmark     = build_pathways("HALLMARK"),
  KEGG         = build_pathways("KEGG"),
  Reactome     = build_pathways("REACTOME"),
  Biocarta     = build_pathways("BIOCARTA"),
  WikiPathways = build_pathways("^WP_"),
  PID          = build_pathways("^PID_")
)
gmt_pathways <- gmt_pathways[sapply(gmt_pathways, length) > 0]

all_markers_gsea <- FindAllMarkers(
  obj,
  only.pos = FALSE,
  min.pct = 0,
  logfc.threshold = 0,
  verbose = FALSE
)
all_markers_gsea$BH_p.val <- p.adjust(all_markers_gsea$p_val, method = "BH")

wb_DE <- createWorkbook()
for (ct in unique(all_markers_gsea$cluster)) {
  markers <- all_markers_gsea %>% filter(cluster == ct)
  sn <- ensureUniqueSheetName(wb_DE, as.character(ct))
  addWorksheet(wb_DE, sn)
  writeData(wb_DE, sn, markers)
}
saveWorkbook(wb_DE, file.path(gsea_dir, "DE_AllCellType.xlsx"), overwrite = TRUE)

dir.create(file.path(gsea_dir, "RNK_Files"), recursive = TRUE, showWarnings = FALSE)
for (ct in unique(all_markers_gsea$cluster)) {
  markers <- all_markers_gsea %>% filter(cluster == ct)
  vec_l2fc <- make_rank_vector(markers, "log2FC")
  vec_stat <- make_rank_vector(markers, "signed_p")
  fwrite(data.table(gene = names(vec_l2fc), rank = as.numeric(vec_l2fc)),
         file.path(gsea_dir, "RNK_Files", paste0(ct, "_log2FC.rnk")),
         sep = "\t", col.names = FALSE)
  fwrite(data.table(gene = names(vec_stat), rank = as.numeric(vec_stat)),
         file.path(gsea_dir, "RNK_Files", paste0(ct, "_signedP.rnk")),
         sep = "\t", col.names = FALSE)
}

wb_GSEA <- createWorkbook()
gsea_results_list <- list()

for (ct in unique(all_markers_gsea$cluster)) {
  markers <- all_markers_gsea %>% filter(cluster == ct)
  if (nrow(markers) < min_deg_for_gsea) next
  vec_l2fc <- make_rank_vector(markers, "log2FC")
  vec_stat <- make_rank_vector(markers, "signed_p")
  if (length(vec_l2fc) < min_deg_for_gsea || length(vec_stat) < min_deg_for_gsea) next

  gsea_results_list[[as.character(ct)]] <- list()

  for (gmt_name in names(gmt_pathways)) {
    pathways <- gmt_pathways[[gmt_name]]
    pathways <- pathways[lengths(pathways) >= gsea_min_size]
    pathways <- pathways[lengths(pathways) <= gsea_max_size]
    if (length(pathways) == 0) next

    fg <- runFgseaByRank(
      pathways = pathways,
      vec_l2fc = vec_l2fc,
      vec_stat = vec_stat,
      gsea_min_size = gsea_min_size,
      gsea_max_size = gsea_max_size,
      gsea_nperm_simple = gsea_nperm_simple,
      gsea_param = gsea_param
    )

    gsea_results_list[[as.character(ct)]][[gmt_name]] <- list(
      log2FC = fg$log2FC,
      signedP = fg$signedP
    )

    sn1 <- ensureUniqueSheetName(wb_GSEA, paste(ct, gmt_name, "l2fc", sep = "_"))
    addWorksheet(wb_GSEA, sn1)
    writeData(wb_GSEA, sn1, fg$log2FC)
    applyCustomConditionalFormatting(wb_GSEA, sn1, fg$log2FC, padj_cutoff)

    sn2 <- ensureUniqueSheetName(wb_GSEA, paste(ct, gmt_name, "stat", sep = "_"))
    addWorksheet(wb_GSEA, sn2)
    writeData(wb_GSEA, sn2, fg$signedP)
    applyCustomConditionalFormatting(wb_GSEA, sn2, fg$signedP, padj_cutoff)
  }
}

saveWorkbook(wb_GSEA, file.path(gsea_dir, "GSEA_AllGenes_ranked_Mouse_m_t2g.xlsx"), overwrite = TRUE)
saveRDS(gsea_results_list, file.path(gsea_dir, "GSEA_AllGenes_ranked_Mouse_m_t2g.rds"))

sig_dir <- file.path(gsea_dir, "Significant_Pathways")
saveSignificantPathways(gsea_results_list, padj_threshold = padj_cutoff, sig_dir = sig_dir)

pin_vs_epi_dir <- file.path(gsea_dir, "PIN_vs_Epi")
dir.create(pin_vs_epi_dir, recursive = TRUE, showWarnings = FALSE)

pin_levels <- c("PIN1", "PIN2", "PIN3", "PIN4", "PIN5", "PIN6")
pin_vs_epi_results <- list()
wb_pin_epi <- createWorkbook()

for (pin in pin_levels) {
  cells_use <- colnames(obj)[obj$reannotated %in% c("Epi", pin)]
  sub_obj <- subset(obj, cells = cells_use)
  Idents(sub_obj) <- "reannotated"

  deg <- FindMarkers(
    sub_obj,
    ident.1 = pin,
    ident.2 = "Epi",
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    verbose = FALSE
  )
  deg$gene <- rownames(deg)
  deg$p_val <- deg$p_val %||% deg$p_val_adj
  deg$cluster <- pin

  vec_l2fc <- make_rank_vector(deg, "log2FC")
  vec_stat <- make_rank_vector(deg, "signed_p")

  pin_vs_epi_results[[pin]] <- list()

  for (gmt_name in names(gmt_pathways)) {
    pathways <- gmt_pathways[[gmt_name]]
    pathways <- pathways[lengths(pathways) >= gsea_min_size]
    pathways <- pathways[lengths(pathways) <= gsea_max_size]
    if (length(pathways) == 0) next

    fg <- runFgseaByRank(
      pathways = pathways,
      vec_l2fc = vec_l2fc,
      vec_stat = vec_stat,
      gsea_min_size = gsea_min_size,
      gsea_max_size = gsea_max_size,
      gsea_nperm_simple = gsea_nperm_simple,
      gsea_param = gsea_param
    )

    pin_vs_epi_results[[pin]][[gmt_name]] <- list(
      log2FC = fg$log2FC,
      signedP = fg$signedP
    )

    sn1 <- ensureUniqueSheetName(wb_pin_epi, paste(pin, gmt_name, "l2fc", sep = "_"))
    addWorksheet(wb_pin_epi, sn1)
    writeData(wb_pin_epi, sn1, fg$log2FC)
    applyCustomConditionalFormatting(wb_pin_epi, sn1, fg$log2FC, padj_cutoff)

    sn2 <- ensureUniqueSheetName(wb_pin_epi, paste(pin, gmt_name, "stat", sep = "_"))
    addWorksheet(wb_pin_epi, sn2)
    writeData(wb_pin_epi, sn2, fg$signedP)
    applyCustomConditionalFormatting(wb_pin_epi, sn2, fg$signedP, padj_cutoff)
  }

  write.csv(deg, file.path(pin_vs_epi_dir, paste0(pin, "_vs_Epi_DE.csv")), row.names = FALSE)
}

saveWorkbook(wb_pin_epi, file.path(pin_vs_epi_dir, "GSEA_PINvsEpi_Mouse_m_t2g.xlsx"), overwrite = TRUE)
saveRDS(pin_vs_epi_results, file.path(pin_vs_epi_dir, "GSEA_PINvsEpi_Mouse_m_t2g.rds"))

selected_pathways_cluster_vs_all <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

selected_pathways_pin_vs_epi <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)

classic_plot_cluster_vs_all <- selected_pathways_cluster_vs_all
classic_plot_pin_vs_epi <- selected_pathways_pin_vs_epi

all_gsea_flat_cluster <- list()
for (ct in names(gsea_results_list)) {
  for (db in names(gsea_results_list[[ct]])) {
    for (rank_method in c("log2FC", "signedP")) {
      dt <- gsea_results_list[[ct]][[db]][[rank_method]]
      if (is.null(dt) || nrow(dt) == 0) next
      dt_copy <- data.table::copy(as.data.table(dt))
      dt_copy[, cell_type := ct]
      dt_copy[, RankMethod := rank_method]
      dt_copy[, Database := db]
      all_gsea_flat_cluster[[paste(ct, db, rank_method, sep = "_")]] <- dt_copy
    }
  }
}
combined_all_cluster <- data.table::rbindlist(all_gsea_flat_cluster, use.names = TRUE, fill = TRUE)

plot_dt_cluster <- combined_all_cluster[pathway %in% selected_pathways_cluster_vs_all] %>%
  dplyr::filter(!is.na(padj), !is.na(NES)) %>%
  dplyr::mutate(
    cell_type = factor(as.character(cell_type), levels = reannotate_order),
    neg_log10_padj = -log10(padj),
    point_size = abs(NES),
    leadingEdge_n = sapply(leadingEdge, function(x) if (is.null(x) || all(is.na(x))) 0 else length(x))
  ) %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::group_by(cell_type, pathway, RankMethod) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = padj < padj_cutoff)

pathway_order_cluster <- selected_pathways_cluster_vs_all

selected_dir_cluster <- file.path(gsea_dir, "Selected_Pathway_Plots_ClusterVsAll")
dir.create(selected_dir_cluster, recursive = TRUE, showWarnings = FALSE)

build_selected_bubble(
  plot_dt = plot_dt_cluster %>% filter(RankMethod == "signedP"),
  cell_order = reannotate_order,
  pathway_order = pathway_order_cluster,
  title_text = "A1281 – Cluster vs All Selected Pathways (signedP-based GSEA)",
  out_file = file.path(selected_dir_cluster, "A1281_ClusterVsAll_SelectedPathways_signedP.tiff"),
  padj_thresh = padj_cutoff,
  width = 10,
  height = 4
)

all_gsea_flat_pin <- list()
for (ct in names(pin_vs_epi_results)) {
  for (db in names(pin_vs_epi_results[[ct]])) {
    for (rank_method in c("log2FC", "signedP")) {
      dt <- pin_vs_epi_results[[ct]][[db]][[rank_method]]
      if (is.null(dt) || nrow(dt) == 0) next
      dt_copy <- data.table::copy(as.data.table(dt))
      dt_copy[, cell_type := ct]
      dt_copy[, RankMethod := rank_method]
      dt_copy[, Database := db]
      all_gsea_flat_pin[[paste(ct, db, rank_method, sep = "_")]] <- dt_copy
    }
  }
}
combined_all_pin <- data.table::rbindlist(all_gsea_flat_pin, use.names = TRUE, fill = TRUE)

plot_dt_pin <- combined_all_pin[pathway %in% selected_pathways_pin_vs_epi] %>%
  dplyr::filter(!is.na(padj), !is.na(NES)) %>%
  dplyr::mutate(
    cell_type = factor(as.character(cell_type), levels = pin_levels),
    neg_log10_padj = -log10(padj),
    point_size = abs(NES),
    leadingEdge_n = sapply(leadingEdge, function(x) if (is.null(x) || all(is.na(x))) 0 else length(x))
  ) %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::group_by(cell_type, pathway, RankMethod) %>%
  dplyr::arrange(padj, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig = padj < padj_cutoff)

pathway_order_pin <- selected_pathways_pin_vs_epi

selected_dir_pin <- file.path(pin_vs_epi_dir, "Selected_Pathway_Plots")
dir.create(selected_dir_pin, recursive = TRUE, showWarnings = FALSE)

build_selected_bubble(
  plot_dt = plot_dt_pin %>% filter(RankMethod == "signedP"),
  cell_order = pin_levels,
  pathway_order = pathway_order_pin,
  title_text = "A1281 – PIN vs Epi Selected Pathways (signedP-based GSEA)",
  out_file = file.path(selected_dir_pin, "A1281_PINvsEpi_SelectedPathways_signedP.tiff"),
  padj_thresh = padj_cutoff,
  width = 10,
  height = 4
)

enrich_dir_cluster <- file.path(gsea_dir, "GSEA_Enrichment_Plots_ClusterVsAll")
dir.create(enrich_dir_cluster, recursive = TRUE, showWarnings = FALSE)

for (ct in unique(all_markers_gsea$cluster)) {
  markers <- all_markers_gsea %>% dplyr::filter(cluster == ct)
  vec_stat <- make_rank_vector(markers, "signed_p")
  if (length(vec_stat) < min_deg_for_gsea) next

  for (pw in classic_plot_cluster_vs_all) {
    if (!pw %in% names(gmt_pathways[["Hallmark"]])) next
    p <- make_gsea_desktop_plot(
      pathway_genes = gmt_pathways[["Hallmark"]][[pw]],
      stats = vec_stat,
      pathway_name = pw,
      phenotype_pos = paste0(ct, "_pos"),
      phenotype_neg = paste0(ct, "_neg"),
      gsea_param = gsea_param
    )
    if (is.null(p)) next
    tiff(
      filename = file.path(enrich_dir_cluster, paste0(ct, "_", pw, "_signedP_GSEA_desktop.tiff")),
      width = 5.2, height = 6.2, units = "in", res = 300, compression = "lzw"
    )
    print(p)
    dev.off()
  }
}

enrich_dir_pin <- file.path(pin_vs_epi_dir, "GSEA_Enrichment_Plots")
dir.create(enrich_dir_pin, recursive = TRUE, showWarnings = FALSE)

for (pin in pin_levels) {
  deg <- read.csv(file.path(pin_vs_epi_dir, paste0(pin, "_vs_Epi_DE.csv")))
  vec_stat <- make_rank_vector(deg, "signed_p")
  if (length(vec_stat) < min_deg_for_gsea) next

  for (pw in classic_plot_pin_vs_epi) {
    if (!pw %in% names(gmt_pathways[["Hallmark"]])) next
    p <- make_gsea_desktop_plot(
      pathway_genes = gmt_pathways[["Hallmark"]][[pw]],
      stats = vec_stat,
      pathway_name = pw,
      phenotype_pos = paste0(pin, "_pos"),
      phenotype_neg = paste0(pin, "_neg"),
      gsea_param = gsea_param
    )
    if (is.null(p)) next
    tiff(
      filename = file.path(enrich_dir_pin, paste0(pin, "_", pw, "_signedP_GSEA_desktop.tiff")),
      width = 5.2, height = 6.2, units = "in", res = 300, compression = "lzw"
    )
    print(p)
    dev.off()
  }
}

saveRDS(gmt_pathways, file.path(gsea_dir, "gmt_pathways.rds"))
saveRDS(pin_levels, file.path(out_dir, "pin_levels.rds"))
saveRDS(obj, file.path(out_dir, "A1281_epi_sub_reannotated_GSEA_complete.rds"))