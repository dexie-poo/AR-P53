save_tiff <- function(plot, filename, width = 12, height = 10, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    device = "tiff",
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw"
  )
}

get_raw_assay <- function(obj) {
  assays <- Assays(obj)
  non_sct <- assays[!grepl("^SCT", assays, ignore.case = TRUE)]
  non_sct[1]
}

ensureUniqueSheetName <- function(wb, base) {
  base <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", base)
  base <- substr(base, 1, 31)
  nm <- base
  i <- 1
  while (nm %in% names(wb)) {
    suffix <- paste0("_", i)
    nm <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1
  }
  nm
}

applyCustomConditionalFormatting <- function(wb, sheet, df, padj_threshold = 0.25) {
  if (!all(c("padj", "NES") %in% names(df))) return()
  if (nrow(df) == 0) return()
  pcol <- int2col(which(names(df) == "padj"))
  ncol2 <- int2col(which(names(df) == "NES"))
  nr <- nrow(df)
  conditionalFormatting(
    wb, sheet, cols = 1:ncol(df), rows = 2:(nr + 1),
    rule = paste0("=$", pcol, "2>", padj_threshold),
    type = "expression",
    style = createStyle(fontColour = "lightgrey")
  )
  conditionalFormatting(
    wb, sheet, cols = 1:ncol(df), rows = 2:(nr + 1),
    rule = paste0("=AND($", pcol, "2<=", padj_threshold, ",$", ncol2, "2>0)"),
    type = "expression",
    style = createStyle(fontColour = "red")
  )
  conditionalFormatting(
    wb, sheet, cols = 1:ncol(df), rows = 2:(nr + 1),
    rule = paste0("=AND($", pcol, "2<=", padj_threshold, ",$", ncol2, "2<0)"),
    type = "expression",
    style = createStyle(fontColour = "blue")
  )
}

addCategory <- function(dt) {
  if (nrow(dt) == 0) return(dt)
  dt[, Category := fifelse(
    grepl("^REACTOME_", pathway), "REACTOME",
    fifelse(grepl("^WP_", pathway), "WikiPathway",
    fifelse(grepl("^(GOBP_|GO_)", pathway), "GO",
    fifelse(grepl("^HALLMARK_", pathway), "HALLMARK",
    fifelse(grepl("^KEGG_", pathway), "KEGG",
    fifelse(grepl("^BIOCARTA_", pathway), "BIOCARTA",
    fifelse(grepl("^PID_", pathway), "PID", "Other")))))))
  ]
  if ("log2err" %in% names(dt)) dt[, log2err := NULL]
  setcolorder(dt, c("Category", setdiff(names(dt), "Category")))
  dt[order(padj)]
}

make_rank_vector <- function(markers, rank_type = c("log2FC", "signed_p")) {
  rank_type <- match.arg(rank_type)
  markers <- markers %>%
    filter(!is.na(gene), gene != "") %>%
    mutate(
      p_val_safe = pmax(p_val, .Machine$double.xmin),
      signed_logp = sign(avg_log2FC) * (-log10(p_val_safe))
    )

  if (rank_type == "log2FC") {
    rank_df <- markers %>%
      group_by(gene) %>%
      summarise(rank_value = mean(avg_log2FC, na.rm = TRUE), .groups = "drop")
  } else {
    rank_df <- markers %>%
      group_by(gene) %>%
      summarise(rank_value = signed_logp[which.max(abs(signed_logp))], .groups = "drop")
  }

  gene_ids <- suppressMessages(
    bitr(
      rank_df$gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Mm.eg.db,
      drop = TRUE
    )
  )

  rank_df <- rank_df %>%
    inner_join(gene_ids, by = c("gene" = "SYMBOL"))

  vec <- setNames(rank_df$rank_value, rank_df$ENTREZID)
  vec <- vec[!duplicated(names(vec))]
  vec <- vec[is.finite(vec)]
  sort(vec, decreasing = TRUE)
}

runFgseaByRank <- function(pathways, vec_l2fc, vec_stat,
                           gsea_min_size, gsea_max_size,
                           gsea_nperm_simple, gsea_param = 0) {
  run_safe <- function(stats_vec) {
    tie_idx <- duplicated(stats_vec) | duplicated(stats_vec, fromLast = TRUE)
    if (any(tie_idx)) {
      stats_vec[tie_idx] <- stats_vec[tie_idx] +
        rank(stats_vec[tie_idx], ties.method = "first") * 1e-12
    }
    out <- tryCatch(
      fgseaMultilevel(
        pathways = pathways,
        stats = stats_vec,
        minSize = gsea_min_size,
        maxSize = gsea_max_size,
        gseaParam = gsea_param,
        nPermSimple = gsea_nperm_simple,
        eps = 0,
        nproc = 1
      ),
      error = function(e) data.table::data.table()
    )
    out <- as.data.table(out)
    if (nrow(out) == 0) return(out)
    addCategory(out)
  }

  list(
    log2FC = run_safe(vec_l2fc),
    signedP = run_safe(vec_stat)
  )
}

make_gsea_desktop_plot <- function(pathway_genes, stats, pathway_name,
                                   phenotype_pos = "pos",
                                   phenotype_neg = "neg",
                                   gsea_param = 0) {
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  pathway_genes <- intersect(pathway_genes, names(stats))
  if (length(pathway_genes) < 3) return(NULL)

  N <- length(stats)
  hits <- names(stats) %in% pathway_genes
  Nm <- N - sum(hits)
  rw <- abs(stats)^gsea_param
  NR <- sum(rw[hits])
  running_es <- cumsum(ifelse(hits, rw / NR, -1 / Nm))
  hit_positions <- which(hits)
  zero_cross <- which.min(abs(stats))

  df_es <- data.frame(rank = seq_along(running_es), ES = running_es)
  df_hits <- data.frame(rank = hit_positions)
  df_metric <- data.frame(rank = seq_along(stats), metric = as.numeric(stats))

  p1 <- ggplot(df_es, aes(x = rank, y = ES)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(color = "green", linewidth = 1.1) +
    labs(title = paste0("Enrichment plot:\n", pathway_name),
         y = "Enrichment score (ES)", x = NULL) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 0, 5)
    )

  p2 <- ggplot(df_hits, aes(x = rank)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1), color = "black", linewidth = 0.25) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))

  p3 <- ggplot(data.frame(x = seq_len(N), y = 1), aes(x = x, y = y, fill = x)) +
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = N / 2, guide = "none") +
    annotate("text", x = max(1, N * 0.02), y = 1.35,
             label = paste0("'", phenotype_pos, "' (positively correlated)"),
             color = "red", hjust = 0, size = 2.7) +
    annotate("text", x = N * 0.72, y = 0.65,
             label = paste0("'", phenotype_neg, "' (negatively correlated)"),
             color = "blue", hjust = 0, size = 2.7) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))

  p4 <- ggplot(df_metric, aes(x = rank, y = metric)) +
    geom_area(fill = "grey75") +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.25) +
    geom_vline(xintercept = zero_cross, color = "grey60", linetype = "dashed", linewidth = 0.3) +
    annotate("text", x = zero_cross, y = 0,
             label = paste0("Zero cross at ", zero_cross),
             vjust = -0.8, size = 2.5) +
    labs(x = "Rank in Ordered Dataset", y = "Ranked list metric") +
    scale_x_continuous(labels = scales::comma) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 5, 5, 5)
    )

  p1 / p2 / p3 / p4 + patchwork::plot_layout(heights = c(3.2, 0.7, 0.35, 2.1))
}

saveSignificantPathways <- function(gsea_results_list, padj_threshold = 0.25, sig_dir) {
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)
  categories <- sort(unique(unlist(lapply(gsea_results_list, names))))
  for (cat_name in categories) {
    cat_dir <- file.path(sig_dir, cat_name)
    dir.create(cat_dir, recursive = TRUE, showWarnings = FALSE)
    for (rank_method in c("log2FC", "signedP")) {
      collected <- list()
      for (ct in names(gsea_results_list)) {
        entry <- gsea_results_list[[ct]][[cat_name]][[rank_method]]
        if (is.null(entry) || nrow(entry) == 0 || !"padj" %in% names(entry)) next
        sig <- entry[padj <= padj_threshold]
        if (nrow(sig) == 0) next
        sig <- copy(sig)
        sig[, cell_type := ct]
        sig[, RankMethod := rank_method]
        setcolorder(sig, c("cell_type", "RankMethod", setdiff(names(sig), c("cell_type", "RankMethod"))))
        collected[[ct]] <- sig
      }
      if (length(collected) == 0) next
      out_dt <- rbindlist(collected, use.names = TRUE, fill = TRUE)
      fwrite(out_dt, file.path(cat_dir, paste0("significant_", rank_method, ".csv")))
      saveRDS(out_dt, file.path(cat_dir, paste0("significant_", rank_method, ".rds")))
    }
  }
}

build_selected_bubble <- function(plot_dt, cell_order, pathway_order, title_text, out_file,
                                  padj_thresh = 0.25, width = 10, height = 4) {
  plot_dt$cell_type <- factor(as.character(plot_dt$cell_type), levels = cell_order)
  plot_dt$pathway <- factor(as.character(plot_dt$pathway), levels = rev(pathway_order))
  plot_dt <- plot_dt %>% filter(!is.na(cell_type), !is.na(pathway))

  if (nrow(plot_dt) == 0) return(NULL)

  nes_min <- min(plot_dt$NES, na.rm = TRUE)
  nes_max <- max(plot_dt$NES, na.rm = TRUE)
  nes_break <- 1.5

  if (!is.finite(nes_min) || !is.finite(nes_max) || nes_min == nes_max) {
    nes_vals <- c(0, 0.5, 1)
    nes_cols <- c("blue", "purple", "#b40426")
  } else {
    bp <- max(0, min(1, (nes_break - nes_min) / (nes_max - nes_min)))
    nes_vals <- sort(unique(c(0, bp * 0.9, bp, 1)))
    nes_cols <- c("blue", "purple", "red", "#b40426")[seq_along(nes_vals)]
  }

  pathway_labels <- setNames(
    vapply(pathway_order, function(s) paste(strwrap(gsub("_", " ", s), width = 32), collapse = "\n"), character(1)),
    pathway_order
  )

  p <- ggplot(plot_dt, aes(x = neg_log10_padj, y = pathway)) +
    geom_point(
      data = plot_dt %>% filter(sig == TRUE),
      aes(size = leadingEdge_n, color = NES),
      alpha = 0.9
    ) +
    scale_color_gradientn(
      colours = nes_cols,
      values = nes_vals,
      name = "NES",
      limits = c(nes_min, nes_max),
      oob = scales::squish
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_y_discrete(
      expand = expansion(add = c(0.7, 0.7)),
      labels = function(x) pathway_labels[x]
    ) +
    scale_size_continuous(
      name = "Leading edge\ngene count",
      range = c(1.5, 5)
    ) +
    labs(
      title = title_text,
      subtitle = paste0("padj < ", padj_thresh, " | shown = significant only"),
      x = expression(-log[10](padj)),
      y = NULL
    ) +
    facet_grid(. ~ cell_type, scales = "fixed", space = "fixed") +
    theme_bw(base_size = 11) +
    theme(
      aspect.ratio = 4,
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text.x = element_blank(),
      strip.background = element_blank(),
      panel.spacing.x = grid::unit(0.45, "lines"),
      plot.margin = margin(12, 18, 12, 18)
    )

  ggsave(
    filename = out_file,
    plot = p,
    device = "tiff",
    width = width,
    height = height,
    dpi = 300,
    units = "in",
    compression = "lzw",
    limitsize = FALSE
  )

  p
}