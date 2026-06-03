#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

obj <- readRDS(file.path(out_dir, "01_A1281_loaded_qc_clustered.rds"))
raw_assay <- readRDS(file.path(out_dir, "raw_assay.rds"))

seurat_char <- as.character(obj$seurat_clusters)
reannotate_labels <- dplyr::case_when(
  seurat_char == "0" ~ "Epi",
  seurat_char %in% c("6") ~ "PIN1",
  seurat_char %in% c("1", "10", "3") ~ "PIN2",
  seurat_char == "7" ~ "PIN3",
  seurat_char %in% c("2", "4", "8") ~ "PIN4",
  seurat_char == "9" ~ "PIN5",
  seurat_char == "5" ~ "PIN6",
  TRUE ~ NA_character_
)

reannotate_order <- c("Epi", "PIN1", "PIN2", "PIN3", "PIN4", "PIN5", "PIN6")
reannotate_colors <- c(
  "Epi"  = "red",
  "PIN1" = "green",
  "PIN2" = "blue",
  "PIN3" = "yellow",
  "PIN4" = "purple",
  "PIN5" = "orange",
  "PIN6" = "cyan"
)

obj$reannotated <- factor(reannotate_labels, levels = reannotate_order)
obj <- subset(obj, cells = colnames(obj)[!is.na(obj$reannotated)])
Idents(obj) <- "reannotated"
DefaultAssay(obj) <- raw_assay

saveRDS(obj, file.path(out_dir, "A1281_epi_sub_reannotated_final.rds"))
saveRDS(obj, file.path(rds_dir, "01_qc_clustered.rds"))
saveRDS(raw_assay, file.path(rds_dir, "raw_assay.rds"))
