#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

epi_sub <- readRDS(file.path(rds_dir, "05_epi_PIN_monocle_seurat.rds"))

DefaultAssay(epi_sub) <- spatial_assay
Idents(epi_sub) <- "reannotated"

gsea_pin_dir <- file.path(gsea_dir, "PIN_vs_Epi")
dir.create(gsea_pin_dir, recursive = TRUE, showWarnings = FALSE)

load("/gpfs/home/leungd02/RDS_Files/m_t2g.RData")

build_pathways <- function(pattern) {
  sub <- m_t2g[grep(pattern, m_t2g$gs_name), ]
  if (nrow(sub) == 0) return(list())
  pw <- split(sub$entrez_gene, sub$gs_name)
  lapply(pw, unique)
}

gmt_pathways <- list(
  Hallmark = build_pathways("HALLMARK"),
  KEGG = build_pathways("KEGG"),
  Reactome = build_pathways("REACTOME"),
  Biocarta = build_pathways("BIOCARTA"),
  WikiPathways = build_pathways("^WP_"),
  PID = build_pathways("^PID_")
)

gmt_pathways <- gmt_pathways[sapply(gmt_pathways, length) > 0]

addCategory <- function(dt) {
  if (nrow(dt) == 0) return(dt)

  dt[, Category := dplyr::case_when(
    grepl("^REACTOME_", pathway) ~ "REACTOME",
    grepl("^WP_", pathway) ~ "WikiPathway",
    grepl("^(GOBP_|GO_)", pathway) ~ "GO",
    grepl("^HALLMARK_", pathway) ~ "HALLMARK",
    grepl("^KEGG_", pathway) ~ "KEGG",
    grepl("^BIOCARTA_", pathway) ~ "BIOCARTA",
    grepl("^PID_", pathway) ~ "PID",
    TRUE ~ "Other"
  )]

  if ("log2err" %in% names(dt)) dt[, log2err := NULL]

  setcolorder(dt, c("Category", setdiff(names(dt), "Category")))

  dt[order(padj)]
}

make_rank_vector <- function(markers, rank_type = c("log2FC", "signed_p")) {
  rank_type <- match.arg(rank_type)

  markers <- markers %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::mutate(
      p_val_safe = pmax(p_val, .Machine$double.xmin),
      signed_logp = sign(avg_log2FC) * (-log10(p_val_safe))
    )

  if (rank_type == "log2FC") {
    rank_df <- markers %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(rank_value = mean(avg_log2FC, na.rm = TRUE), .groups = "drop")
  } else {
    rank_df <- markers %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(rank_value = signed_logp[which.max(abs(signed_logp))], .groups = "drop")
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
    dplyr::inner_join(gene_ids, by = c("gene" = "SYMBOL"))

  vec <- setNames(rank_df$rank_value, rank_df$ENTREZID)
  vec <- vec[!duplicated(names(vec))]
  vec <- vec[is.finite(vec)]
  sort(vec, decreasing = TRUE)
}

run_fgsea_pair <- function(pathways, vec_l2fc, vec_stat) {
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

pin_levels <- c("PIN1", "PIN2", "PIN3", "PIN4", "PIN5", "PIN6")

gsea_results_list <- list()
all_deg_list <- list()

for (pin in pin_levels) {
  message("Running DE and GSEA: ", pin, " vs Epi")

  cells_use <- colnames(epi_sub)[epi_sub$reannotated %in% c("Epi", pin)]
  obj_sub <- subset(epi_sub, cells = cells_use)

  Idents(obj_sub) <- "reannotated"

  markers <- FindMarkers(
    obj_sub,
    ident.1 = pin,
    ident.2 = "Epi",
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    verbose = FALSE
  )

  markers <- markers %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(comparison = paste0(pin, "_vs_Epi"))

  all_deg_list[[pin]] <- markers

  vec_l2fc <- make_rank_vector(markers, "log2FC")
  vec_stat <- make_rank_vector(markers, "signed_p")

  gsea_results_list[[pin]] <- list()

  for (gmt_name in names(gmt_pathways)) {
    pathways <- gmt_pathways[[gmt_name]]
    pathways <- pathways[lengths(pathways) >= gsea_min_size]
    pathways <- pathways[lengths(pathways) <= gsea_max_size]

    if (length(pathways) == 0) next

    fg <- run_fgsea_pair(pathways, vec_l2fc, vec_stat)

    gsea_results_list[[pin]][[gmt_name]] <- list(
      log2FC = fg$log2FC,
      signedP = fg$signedP
    )
  }
}

all_deg <- dplyr::bind_rows(all_deg_list)

write.csv(
  all_deg,
  file.path(gsea_pin_dir, "DE_PIN_vs_Epi_all.csv"),
  row.names = FALSE
)

saveRDS(
  gsea_results_list,
  file.path(gsea_pin_dir, "GSEA_PINvsEpi_Mouse_m_t2g.rds")
)

saveRDS(
  epi_sub,
  file.path(rds_dir, "06_gsea_complete.rds")
)

message("Module 6 complete.")
