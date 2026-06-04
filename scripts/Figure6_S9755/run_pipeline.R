#!/usr/bin/env Rscript
# =============================================================
# Figure 6 / S9755 complete modular pipeline runner
# =============================================================
steps <- c(
  "scripts/06_S9755/01_load_raw_preprocess.R",
  "scripts/06_S9755/02_broad_clustering_annotation.R",
  "scripts/06_S9755/03_epi_PIN_subset_annotation.R",
  "scripts/06_S9755/04_marker_panels_plots.R",
  "scripts/06_S9755/05_monocle3.R",
  "scripts/06_S9755/06_GSEA.R"
)

for (s in steps) {
  message("\n=============================================================")
  message("Running: ", s)
  message("=============================================================")
  source(s)
  gc()
}
message("Figure 6 / S9755 pipeline complete.")
