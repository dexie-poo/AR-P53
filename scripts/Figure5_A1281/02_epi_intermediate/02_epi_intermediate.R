#!/usr/bin/env Rscript

source("scripts/00_setup/00_config.R")

message("Loading confirmed intermediate epithelial object...")

epi_sub <- readRDS(
  "/gpfs/scratch/leungd02/05262026/A1281_epi_reannotated_reclustered.rds"
)

DefaultAssay(epi_sub) <- spatial_assay

message("Object loaded: ", ncol(epi_sub), " cells")
message("Metadata columns:")
print(colnames(epi_sub@meta.data))

if ("new_annot" %in% colnames(epi_sub@meta.data)) {
  message("new_annot:")
  print(table(epi_sub$new_annot, useNA = "always"))
}

if ("reannotated" %in% colnames(epi_sub@meta.data)) {
  message("existing reannotated:")
  print(table(epi_sub$reannotated, useNA = "always"))
}

saveRDS(
  epi_sub,
  file.path(rds_dir, "03_epi_intermediate_annotated.rds")
)

message("Saved: outputs/rds/03_epi_intermediate_annotated.rds")
