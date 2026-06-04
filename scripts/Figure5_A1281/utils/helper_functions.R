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

get_spatial_assay <- function(obj) {
  assays <- names(obj@assays)
  raw_assays <- assays[!grepl("^SCT", assays, ignore.case = TRUE)]
  if (length(raw_assays) == 0) {
    stop("No non-SCT assay found.")
  }
  if ("Spatial.024um" %in% raw_assays) return("Spatial.024um")
  if ("Spatial" %in% raw_assays) return("Spatial")
  raw_assays[1]
}
