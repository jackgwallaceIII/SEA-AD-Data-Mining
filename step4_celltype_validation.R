# =============================================================================
# SEA-AD TRAPPC9 Analysis — Step 4: Cell Type Marker Validation
# =============================================================================
# Validates that "Immune" and "Astrocyte" subclasses express canonical markers.
# Prerequisites: step1 h5ad files must exist for each region.
# Output: results/<REGION>/celltype_marker_validation_dotplot.pdf
# =============================================================================

library(Seurat)
library(zellkonverter)
library(SingleCellExperiment)
library(ggplot2)

REGIONS <- c("HIP", "MTG", "PFC")

# Canonical cell type markers
ASTROCYTE_MARKERS <- c("GFAP", "AQP4", "SLC1A2", "ALDH1L1", "S100B")
MICROGLIA_MARKERS <- c("CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "TREM2")
ALL_MARKERS       <- c(ASTROCYTE_MARKERS, MICROGLIA_MARKERS)

for (region in REGIONS) {
  message("\n", strrep("-", 50))
  message("Validating cell types: ", region)

  h5ad_path <- file.path("seaad_data", region,
                          paste0("SEAAD_", region, "_subset_for_R.h5ad"))
  if (!file.exists(h5ad_path)) {
    message("  File not found — skipping: ", h5ad_path); next
  }

  sce <- readH5AD(h5ad_path, verbose = FALSE)
  seu <- as.Seurat(sce, counts = "X", data = NULL)
  seu@meta.data <- as.data.frame(
    lapply(colData(sce), function(x) if (is.factor(x)) as.character(x) else as.vector(x)),
    stringsAsFactors = FALSE, row.names = colnames(sce)
  )
  rm(sce); gc()

  DefaultAssay(seu) <- "originalexp"
  seu <- SetAssayData(seu, layer = "data",
                      new.data = GetAssayData(seu, layer = "counts"))
  seu@meta.data$Subclass <- as.character(seu@meta.data$Subclass)
  Idents(seu) <- "Subclass"

  markers_present <- ALL_MARKERS[ALL_MARKERS %in% rownames(seu)]
  message("  Markers found: ", paste(markers_present, collapse = ", "))

  if (length(markers_present) == 0) {
    message("  No markers found — skipping plot"); rm(seu); gc(); next
  }

  # Dot plot: rows = marker genes, columns = cell types
  p <- DotPlot(seu,
               features  = markers_present,
               cols      = c("lightgrey", "#2166AC"),
               dot.scale = 9) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
      axis.text.y  = element_text(face = "italic", size = 10),
      plot.title   = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey45", size = 10),
      panel.grid.major = element_line(color = "grey92"),
      strip.background = element_blank()
    ) +
    # Add a horizontal separator between astrocyte and microglial marker groups
    geom_hline(yintercept = length(ASTROCYTE_MARKERS[ASTROCYTE_MARKERS %in% markers_present]) + 0.5,
               linetype = "dashed", color = "grey60", linewidth = 0.6) +
    labs(
      title    = paste("Cell type marker validation —", region),
      subtitle = "Dot size = % expressing  |  Color = mean expression",
      x        = NULL, y        = NULL
    )

  out_dir <- file.path("results", region)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "celltype_marker_validation_dotplot.pdf")
  ggsave(out_file, p, width = 6, height = 7)
  message("  Saved: ", out_file)

  rm(seu); gc()
}

message("\nDone.")
