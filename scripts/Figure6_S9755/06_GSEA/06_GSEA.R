#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 module 06
# Differential expression and fgsea pathway analysis for final Epi/PIN states
# =============================================================
source("scripts/06_S9755/00_config.R")
suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(scales)
})

GSEA_DIR <- file.path(OUT_DIR, "GSEA")
dir.create(GSEA_DIR, recursive = TRUE, showWarnings = FALSE)

obj_path <- file.path(RDS_DIR, "05_S9755_EpiPIN_subset_annotated_v2_monocle.rds")
if (!file.exists(obj_path)) obj_path <- file.path(RDS_DIR, "03_S9755_EpiPIN_subset_annotated_v2.rds")
obj <- readRDS(obj_path)
assay <- get_spatial_assay(obj)
DefaultAssay(obj) <- assay
Idents(obj) <- "epi_pin_cell_type_v2"
obj$epi_pin_cell_type_v2 <- factor(as.character(obj$epi_pin_cell_type_v2), levels = EPIPIN_LEVELS)

padj_cutoff       <- 0.25
min_deg_for_gsea  <- 10
gsea_min_size     <- 5
gsea_max_size     <- 500
gsea_nperm_simple <- 1000
gsea_param        <- 0
set.seed(1234)

force_include_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_ANDROGEN_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)

message("Loading m_t2g.RData...")
load("/gpfs/home/leungd02/RDS_Files/m_t2g.RData")
build_pathways <- function(pattern) {
  sub <- m_t2g[grep(pattern, m_t2g$gs_name), ]
  if (nrow(sub) == 0) return(list())
  lapply(split(sub$entrez_gene, sub$gs_name), unique)
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

ensureUniqueSheetName <- function(wb, base) {
  base <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", base)
  base <- substr(base, 1, 31)
  nm <- base; i <- 1
  while (nm %in% names(wb)) {
    suffix <- paste0("_", i)
    nm <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1
  }
  nm
}

addCategory <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  dt[, Category := fifelse(grepl("^REACTOME_", pathway), "REACTOME",
    fifelse(grepl("^WP_", pathway), "WikiPathway",
    fifelse(grepl("^(GOBP_|GO_)", pathway), "GO",
    fifelse(grepl("^HALLMARK_", pathway), "HALLMARK",
    fifelse(grepl("^KEGG_", pathway), "KEGG",
    fifelse(grepl("^BIOCARTA_", pathway), "BIOCARTA",
    fifelse(grepl("^PID_", pathway), "PID", "Other")))))))]
  if ("log2err" %in% names(dt)) dt[, log2err := NULL]
  setcolorder(dt, c("Category", setdiff(names(dt), "Category")))
  dt[order(padj)]
}

make_rank_vector <- function(markers, rank_type = c("log2FC", "signed_p")) {
  rank_type <- match.arg(rank_type)
  markers <- markers %>% filter(!is.na(gene), gene != "") %>%
    mutate(p_val_safe = pmax(p_val, .Machine$double.xmin),
           signed_logp = sign(avg_log2FC) * (-log10(p_val_safe)))
  if (rank_type == "log2FC") {
    rank_df <- markers %>% group_by(gene) %>% summarise(rank_value = mean(avg_log2FC, na.rm = TRUE), .groups = "drop")
  } else {
    rank_df <- markers %>% group_by(gene) %>% summarise(rank_value = signed_logp[which.max(abs(signed_logp))], .groups = "drop")
  }
  gene_ids <- suppressMessages(bitr(rank_df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE))
  rank_df <- inner_join(rank_df, gene_ids, by = c("gene" = "SYMBOL"))
  vec <- setNames(rank_df$rank_value, rank_df$ENTREZID)
  vec <- vec[!duplicated(names(vec))]
  vec <- vec[is.finite(vec)]
  sort(vec, decreasing = TRUE)
}

run_fgsea_pair <- function(pathways, vec_l2fc, vec_stat) {
  run_safe <- function(stats_vec) {
    tie_idx <- duplicated(stats_vec) | duplicated(stats_vec, fromLast = TRUE)
    if (any(tie_idx)) stats_vec[tie_idx] <- stats_vec[tie_idx] + rank(stats_vec[tie_idx], ties.method = "first") * 1e-12
    out <- tryCatch(fgseaMultilevel(pathways = pathways, stats = stats_vec, minSize = gsea_min_size,
                                    maxSize = gsea_max_size, gseaParam = gsea_param,
                                    nPermSimple = gsea_nperm_simple, eps = 0, nproc = 1),
                    error = function(e) data.table())
    out <- as.data.table(out)
    if (nrow(out) == 0) return(out)
    addCategory(out)
  }
  list(log2FC = run_safe(vec_l2fc), signedP = run_safe(vec_stat))
}

message("Running FindAllMarkers...")
all_markers <- FindAllMarkers(obj, only.pos = FALSE, min.pct = 0, logfc.threshold = 0, verbose = FALSE)
all_markers$BH_p.val <- p.adjust(all_markers$p_val, method = "BH")
write.csv(all_markers, file.path(TABLE_DIR, "06_S9755_EpiPIN_DE_all_markers.csv"), row.names = FALSE)

wb_de <- createWorkbook()
for (ct in unique(all_markers$cluster)) {
  sn <- ensureUniqueSheetName(wb_de, as.character(ct))
  addWorksheet(wb_de, sn)
  writeData(wb_de, sn, all_markers %>% filter(cluster == ct))
}
saveWorkbook(wb_de, file.path(GSEA_DIR, "S9755_EpiPIN_DE_AllCellType.xlsx"), overwrite = TRUE)

message("Writing RNK and running GSEA...")
dir.create(file.path(GSEA_DIR, "RNK_Files"), recursive = TRUE, showWarnings = FALSE)
wb_gsea <- createWorkbook()
gsea_results <- list()

for (ct in unique(all_markers$cluster)) {
  message("Processing ", ct)
  markers <- all_markers %>% filter(cluster == ct)
  if (nrow(markers) < min_deg_for_gsea) next
  vec_l2fc <- make_rank_vector(markers, "log2FC")
  vec_stat <- make_rank_vector(markers, "signed_p")
  fwrite(data.table(gene = names(vec_l2fc), rank = as.numeric(vec_l2fc)),
         file.path(GSEA_DIR, "RNK_Files", paste0(ct, "_log2FC.rnk")), sep = "\t", col.names = FALSE)
  fwrite(data.table(gene = names(vec_stat), rank = as.numeric(vec_stat)),
         file.path(GSEA_DIR, "RNK_Files", paste0(ct, "_signedP.rnk")), sep = "\t", col.names = FALSE)
  if (length(vec_l2fc) < min_deg_for_gsea || length(vec_stat) < min_deg_for_gsea) next
  gsea_results[[as.character(ct)]] <- list()
  for (db in names(gmt_pathways)) {
    pathways <- gmt_pathways[[db]]
    pathways <- pathways[lengths(pathways) >= gsea_min_size & lengths(pathways) <= gsea_max_size]
    if (length(pathways) == 0) next
    fg <- run_fgsea_pair(pathways, vec_l2fc, vec_stat)
    gsea_results[[as.character(ct)]][[db]] <- fg
    for (rank_method in names(fg)) {
      sn <- ensureUniqueSheetName(wb_gsea, paste(ct, db, rank_method, sep = "_"))
      addWorksheet(wb_gsea, sn)
      writeData(wb_gsea, sn, fg[[rank_method]])
    }
  }
  gc()
}

saveWorkbook(wb_gsea, file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_AllGenes_ranked_Mouse_m_t2g.xlsx"), overwrite = TRUE)
saveRDS(gsea_results, file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_AllGenes_ranked_Mouse_m_t2g.rds"))

message("Collecting significant and selected pathways...")
n_flat <- list()
for (ct in names(gsea_results)) {
  for (db in names(gsea_results[[ct]])) {
    for (rank_method in c("log2FC", "signedP")) {
      dt <- gsea_results[[ct]][[db]][[rank_method]]
      if (is.null(dt) || nrow(dt) == 0) next
      tmp <- copy(as.data.table(dt))
      tmp[, cell_type := ct]
      tmp[, Database := db]
      tmp[, RankMethod := rank_method]
      n_flat[[paste(ct, db, rank_method, sep = "_")]] <- tmp
    }
  }
}
combined <- rbindlist(n_flat, use.names = TRUE, fill = TRUE)
fwrite(combined, file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_all_flat_results.csv"))

sig <- combined[!is.na(padj) & padj <= padj_cutoff]
fwrite(sig, file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_significant_padj025.csv"))
selected <- combined[pathway %in% force_include_pathways]
fwrite(selected, file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_forced_selected_pathways.csv"))
saveRDS(selected, file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_forced_selected_pathways.rds"))

if (nrow(selected) > 0) {
  plot_dt <- selected %>% filter(!is.na(padj), !is.na(NES)) %>%
    mutate(cell_type = factor(cell_type, levels = EPIPIN_LEVELS), neg_log10_padj = -log10(padj), point_size = abs(NES)) %>%
    group_by(cell_type, pathway, RankMethod) %>% arrange(padj, desc(abs(NES)), .by_group = TRUE) %>% slice(1) %>% ungroup()
  for (rank_method in c("log2FC", "signedP")) {
    pd <- plot_dt %>% filter(RankMethod == rank_method)
    if (nrow(pd) == 0) next
    p <- ggplot(pd, aes(x = cell_type, y = pathway)) +
      geom_point(aes(size = point_size, color = NES), alpha = 0.9) +
      scale_color_gradientn(colours = c("blue", "purple", "red", "#b40426"), name = "NES") +
      theme_bw(base_size = 12) + labs(title = paste0("S9755 selected pathways: ", rank_method), x = NULL, y = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))
    save_tiff(p, file.path(GSEA_DIR, paste0("S9755_EpiPIN_SelectedPathways_BubblePlot_", rank_method, ".tiff")), width = 12, height = 6)
  }
}

write.csv(data.frame(
  Parameter = c("padj_cutoff", "min_gene_set_size", "max_gene_set_size", "nperm", "gseaParam", "object", "cell_type_column", "cell_types"),
  Value = c(padj_cutoff, gsea_min_size, gsea_max_size, gsea_nperm_simple, gsea_param, basename(obj_path), "epi_pin_cell_type_v2", paste(EPIPIN_LEVELS, collapse = ", "))
), file.path(GSEA_DIR, "S9755_EpiPIN_GSEA_run_parameters.csv"), row.names = FALSE)

message("GSEA complete: ", GSEA_DIR)


library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(data.table)

padj_cutoff <- 0.25

message("\n===== S9755 EpiPIN Final Selected Pathway Plot =====")

gsea_dir_s9755 <- "/gpfs/scratch/leungd02/05292026/GSEA"

selected_dir_s9755 <- file.path(
  gsea_dir_s9755,
  "Selected_Pathway_Plots",
  "SigFormat_ManualOrder"
)

dir.create(selected_dir_s9755, recursive = TRUE, showWarnings = FALSE)

cell_order_s9755 <- c("Epi1", "Epi2", "Epi3", "Epi4", "PIN")

pathway_fixed_order_s9755 <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_ANDROGEN_RESPONSE"
)

message("Loading S9755 GSEA RDS ...")

gsea_rds_s9755 <- readRDS(
  file.path(
    gsea_dir_s9755,
    "S9755_EpiPIN_GSEA_AllGenes_ranked_Mouse_m_t2g.rds"
  )
)

all_flat_s9755 <- list()

for (ct in names(gsea_rds_s9755)) {
  for (db in names(gsea_rds_s9755[[ct]])) {
    for (rank_method in c("log2FC", "signedP")) {
      dt <- gsea_rds_s9755[[ct]][[db]][[rank_method]]
      if (is.null(dt) || nrow(dt) == 0) next

      dt_copy <- data.table::copy(data.table::as.data.table(dt))
      dt_copy[, cell_type  := ct]
      dt_copy[, RankMethod := rank_method]
      dt_copy[, Database   := db]

      all_flat_s9755[[paste(ct, db, rank_method, sep = "_")]] <- dt_copy
    }
  }
}

combined_all_s9755 <- data.table::rbindlist(
  all_flat_s9755,
  use.names = TRUE,
  fill = TRUE
)

message("S9755 flattened rows: ", nrow(combined_all_s9755))

plot_dt_s9755 <- combined_all_s9755[
  pathway %in% pathway_fixed_order_s9755
] %>%
  dplyr::filter(!is.na(padj), !is.na(NES)) %>%
  dplyr::mutate(
    cell_type = factor(
      as.character(cell_type),
      levels = cell_order_s9755
    ),
    neg_log10_padj = -log10(padj),
    point_size     = abs(NES),
    leadingEdge_n  = sapply(leadingEdge, function(x) {
      if (is.null(x) || all(is.na(x))) return(0)
      length(x)
    })
  ) %>%
  dplyr::filter(!is.na(cell_type)) %>%
  dplyr::group_by(cell_type, pathway, RankMethod) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(abs(NES)),
    .by_group = TRUE
  ) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

plot_dt_s9755$pathway <- factor(
  as.character(plot_dt_s9755$pathway),
  levels = rev(pathway_fixed_order_s9755)
)

plot_dt_s9755$cell_type <- factor(
  as.character(plot_dt_s9755$cell_type),
  levels = cell_order_s9755
)

message(
  "S9755 Pathways found: ",
  paste(
    intersect(
      pathway_fixed_order_s9755,
      unique(as.character(plot_dt_s9755$pathway))
    ),
    collapse = ", "
  )
)

message(
  "S9755 Pathways NOT found: ",
  paste(
    setdiff(
      pathway_fixed_order_s9755,
      unique(as.character(plot_dt_s9755$pathway))
    ),
    collapse = ", "
  )
)

message(
  "S9755 Leading edge range: ",
  min(plot_dt_s9755$leadingEdge_n, na.rm = TRUE),
  " – ",
  max(plot_dt_s9755$leadingEdge_n, na.rm = TRUE)
)

nes_min_s9755 <- min(plot_dt_s9755$NES, na.rm = TRUE)
nes_max_s9755 <- max(plot_dt_s9755$NES, na.rm = TRUE)
nes_break     <- 1.5

if (
  !is.finite(nes_min_s9755) ||
  !is.finite(nes_max_s9755) ||
  nes_min_s9755 == nes_max_s9755
) {
  nes_vals_s9755 <- c(0, 0.5, 1)
  nes_cols_s9755 <- c("blue", "purple", "#b40426")
} else {
  bp <- max(
    0,
    min(
      1,
      (nes_break - nes_min_s9755) /
        (nes_max_s9755 - nes_min_s9755)
    )
  )

  nes_vals_s9755 <- sort(unique(c(0, bp * 0.9, bp, 1)))

  nes_cols_s9755 <- c(
    "blue",
    "purple",
    "red",
    "#b40426"
  )[seq_along(nes_vals_s9755)]
}

make_sig_format_plot_s9755 <- function(
  rank_label,
  plot_width  = 10,
  plot_height = 3
) {
  pd <- plot_dt_s9755 %>%
    dplyr::mutate(pathway = as.character(pathway)) %>%
    dplyr::filter(RankMethod == rank_label) %>%
    dplyr::mutate(
      sig = padj < padj_cutoff,
      pathway = factor(
        pathway,
        levels = rev(pathway_fixed_order_s9755)
      )
    )

  if (nrow(pd) == 0) {
    warning("No data for: ", rank_label)
    return(NULL)
  }

  p <- ggplot(
    pd,
    aes(x = neg_log10_padj, y = pathway)
  ) +
    geom_point(
      data = dplyr::filter(pd, sig == TRUE),
      aes(size = leadingEdge_n, color = NES),
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = nes_cols_s9755,
      values  = nes_vals_s9755,
      name    = "NES",
      limits  = c(nes_min_s9755, nes_max_s9755),
      oob     = scales::squish
    ) +
    scale_x_continuous(
      limits = c(-1.5, 10),
      breaks = c(0, 5),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_discrete(
      expand = expansion(add = c(0.7, 0.7)),
      labels = function(x) {
        vapply(
          x,
          function(s) {
            paste(
              strwrap(gsub("_", " ", s), width = 32),
              collapse = "\n"
            )
          },
          character(1)
        )
      }
    ) +
    scale_size_continuous(
      name   = "Leading edge\ngene count",
      range  = c(1.5, 5),
      breaks = c(10, 40, 80),
      limits = c(0, 120)
    ) +
    labs(
      title = paste0(
        "S9755 EpiPIN – Selected Pathways\n(",
        rank_label,
        "-based GSEA)"
      ),
      subtitle = paste0(
        "padj < ",
        padj_cutoff,
        "  |  Shown = significant only"
      ),
      x = expression(-log[10](padj)),
      y = NULL
    ) +
    facet_grid(
      . ~ cell_type,
      scales = "fixed",
      space  = "fixed"
    ) +
    theme_bw(base_size = 11) +
    theme(
      aspect.ratio     = 4,
      axis.text.x      = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text.x     = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x  = grid::unit(0.45, "lines"),
      plot.margin      = margin(12, 18, 12, 18)
    )

  out_tiff <- file.path(
    selected_dir_s9755,
    paste0(
      "S9755_EpiPIN_SigFormat_ManualOrder_",
      rank_label,
      ".tiff"
    )
  )

  ggsave(
    filename = out_tiff,
    plot = p,
    device = "tiff",
    width = plot_width,
    height = plot_height,
    dpi = 300,
    units = "in",
    compression = "lzw",
    limitsize = FALSE
  )

  message(
    "Saved (",
    plot_width,
    " x ",
    plot_height,
    " in): ",
    out_tiff
  )

  return(p)
}

plot_width_s9755  <- 10
plot_height_s9755 <- 3

p_s9755_signedP <- make_sig_format_plot_s9755(
  "signedP",
  plot_width_s9755,
  plot_height_s9755
)

p_s9755_log2FC <- make_sig_format_plot_s9755(
  "log2FC",
  plot_width_s9755,
  plot_height_s9755
)

message("\n✅ S9755 final selected pathway plots complete.")
message("Output: ", selected_dir_s9755)
