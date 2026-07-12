# =============================================================================
# SEA-AD TRAPPC9 Analysis Pipeline — Step 2: Seurat Analysis
# =============================================================================
# Dataset : Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)
# Gene    : TRAPPC9
# Question: How does TRAPPC9 expression differ across cell types between
#           healthy controls (ADNC = "Not AD") and AD donors (ADNC = "High")?
#
# Prerequisites:
#   1. Run step1_download_preprocess.py first to generate
#      seaad_data/<REGION>/SEAAD_<REGION>_subset_for_R.h5ad for each region
#   2. Install R packages listed below
#   3. Set H5AD_PATH below (or loop over regions — see REGIONS_TO_RUN)
#
# Outputs (saved to results/):
#   - TRAPPC9_dotplot_celltype_x_ADNC.pdf
#   - TRAPPC9_violin_by_subclass.pdf
#   - TRAPPC9_umap.pdf
#   - TRAPPC9_pct_expressing_barplot.pdf
#   - TRAPPC9_stats_per_celltype.csv
#   - TRAPPC9_pseudobulk_DE.csv
#   - TRAPPC9_summary_heatmap.pdf
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================

# Install Bioconductor if not present
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_pkgs <- c(
  # Core
  "Seurat",        # scRNA-seq analysis
  "SeuratObject",
  "zellkonverter", # read h5ad in R (Bioconductor)
  "SingleCellExperiment",
  # Data wrangling
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "stringr",
  # Visualization
  "ggplot2",
  "ggpubr",       # stat_compare_means
  "patchwork",    # combine ggplots
  "RColorBrewer",
  "viridis",
  "pheatmap",
  # Statistical testing
  "broom"
)

bioc_pkgs <- c("zellkonverter", "SingleCellExperiment", "BiocParallel",
               "DESeq2", "edgeR")  # for pseudobulk DE

install_if_missing <- function(pkgs, installer = install.packages) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "))
    installer(missing, ask = FALSE)
  }
}

install_if_missing(setdiff(required_pkgs, bioc_pkgs))
install_if_missing(bioc_pkgs, installer = function(p, ...) BiocManager::install(p, ...))

suppressPackageStartupMessages({
  library(Seurat)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(pheatmap)
  library(broom)
})


# =============================================================================
# 1. CONFIGURATION
# =============================================================================

# Regions to run. The loop below iterates over all of them automatically.
# Remove or comment out any region you don't want to process.
REGIONS_TO_RUN <- c("HIP", "MTG", "PFC")

# (H5AD_PATH is set inside the loop — do not set it here)
H5AD_PATH <- NULL   # placeholder; overwritten each iteration

GENE <- "TRAPPC9"

# Cell types to keep. Set to NULL to load all cell types.
CELL_TYPES_KEEP <- c("Immune", "Astrocyte")

# ADNC groups — controls vs AD.
CONTROL_LABEL <- "Not AD"
AD_LABEL      <- "High"

# Colors for ADNC groups (ordered from low to high pathology)
ADNC_COLORS <- c(
  "Not AD"       = "#2166AC",
  "Low"          = "#74ADD1",
  "Intermediate" = "#F4A582",
  "High"         = "#D6604D"
)

# Cell-type resolution to use for main plots
CELLTYPE_COL <- "Subclass"

# Minimum cells per group required for statistical testing
MIN_CELLS_FOR_STATS <- 10


# =============================================================================
# MAIN LOOP — runs once per region
# =============================================================================

for (region in REGIONS_TO_RUN) {

message("\n", strrep("=", 65))
message("PROCESSING REGION: ", region)
message(strrep("=", 65))

H5AD_PATH <- file.path("seaad_data", region,
                        paste0("SEAAD_", region, "_subset_for_R.h5ad"))

if (!file.exists(H5AD_PATH)) {
  message("  File not found — skipping: ", H5AD_PATH)
  next
}

# Auto-detect region label from filename (e.g. "SEAAD_HIP_subset..." → "HIP")
region_match <- regmatches(H5AD_PATH,
                           regexpr("SEAAD_([^_/]+)_subset", H5AD_PATH))
REGION_LABEL <- if (length(region_match) > 0)
  sub("SEAAD_(.+)_subset.*", "\\1", region_match) else region

# Per-region output folder (e.g. results/HIP/, results/MTG/)
OUTPUT_DIR <- file.path("results", REGION_LABEL)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

message("\n[1/7] Loading h5ad and converting to Seurat...")
message("  Region : ", REGION_LABEL)
message("  File   : ", H5AD_PATH)

# Read h5ad → SingleCellExperiment
sce <- readH5AD(H5AD_PATH, verbose = FALSE, hdf5_backed = TRUE)
message("  Full dataset: ", ncol(sce), " cells, ", nrow(sce), " genes")

# Subset to cell types of interest BEFORE converting to Seurat.
# For large regions (e.g. MEC ~970k cells), loading everything into Seurat
# can require 40-50GB RAM. Filtering at the SCE level is fast and cheap.
if (!is.null(CELL_TYPES_KEEP)) {
  subclass_col <- if ("Subclass" %in% names(colData(sce))) "Subclass" else NULL
  if (!is.null(subclass_col)) {
    keep_mask <- colData(sce)[[subclass_col]] %in% CELL_TYPES_KEEP
    sce <- sce[, keep_mask]
    message("  Subsetted to ", paste(CELL_TYPES_KEEP, collapse = ", "),
            ": ", ncol(sce), " cells retained")
  } else {
    message("  WARNING: 'Subclass' column not found — loading all cell types")
  }
}

# Convert to Seurat
seu <- as.Seurat(sce, counts = "X", data = NULL)

# Move obs metadata to Seurat meta.data, decoding factors → character strings.
# as.data.frame(colData()) alone converts h5ad categoricals to integer codes;
# lapply with as.character() preserves the string labels (e.g. "Astrocyte").
seu@meta.data <- as.data.frame(
  lapply(colData(sce), function(x) if (is.factor(x)) as.character(x) else as.vector(x)),
  stringsAsFactors = FALSE,
  row.names        = colnames(sce)
)

# If UMAP was stored, add it
if ("X_umap" %in% names(reducedDims(sce))) {
  umap_mat <- reducedDim(sce, "X_umap")
  colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
  rownames(umap_mat) <- colnames(seu)
  seu[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_")
  message("  UMAP coordinates imported.")
}

rm(sce); gc()

message("  Seurat object: ", ncol(seu), " cells, ", nrow(seu), " genes")
message("  Metadata columns: ", paste(colnames(seu@meta.data), collapse = ", "))


# =============================================================================
# 3. METADATA CLEANUP & GROUPING
# =============================================================================

message("\n[2/7] Cleaning metadata...")

# Ensure ADNC column exists
if (!"ADNC" %in% colnames(seu@meta.data)) {
  stop("ADNC column not found. Check columns: ", paste(colnames(seu@meta.data), collapse = ", "))
}

# Remove any remaining Reference donors (safety check)
seu <- subset(seu, ADNC %in% c("Not AD", "Low", "Intermediate", "High"))

# Ordered factor for ADNC
seu$ADNC <- factor(seu$ADNC,
                   levels = c("Not AD", "Low", "Intermediate", "High"),
                   ordered = TRUE)

# Binary grouping for simpler comparisons
seu$Disease_Group <- ifelse(seu$ADNC == "Not AD", "Control (Not AD)",
                     ifelse(seu$ADNC == "High",   "AD (High)",
                            as.character(seu$ADNC)))

# Make sure cell-type columns are character
# NOTE: seu[[col]] returns a data.frame in Seurat; use seu@meta.data[[col]] for the vector
if (CELLTYPE_COL %in% colnames(seu@meta.data)) {
  seu@meta.data[[CELLTYPE_COL]] <- as.character(seu@meta.data[[CELLTYPE_COL]])
}
if ("Class" %in% colnames(seu@meta.data)) {
  seu@meta.data[["Class"]] <- as.character(seu@meta.data[["Class"]])
}

Idents(seu) <- CELLTYPE_COL

message("  ADNC distribution:")
print(table(seu$ADNC))
message("\n  Cell type distribution (", CELLTYPE_COL, "):")
print(sort(table(seu@meta.data[[CELLTYPE_COL]]), decreasing = TRUE))


# =============================================================================
# 4. GENE EXPRESSION CHECK
# =============================================================================

message("\n[3/7] Checking ", GENE, " expression...")

if (!GENE %in% rownames(seu)) {
  # Try case-insensitive match
  partial <- grep(toupper(GENE), toupper(rownames(seu)), value = TRUE)
  if (length(partial) > 0) {
    message("  ", GENE, " not found. Possible matches: ", paste(partial, collapse = ", "))
    GENE <- partial[1]
    message("  Using: ", GENE)
  } else {
    stop(GENE, " not found in the dataset. Check gene names with head(rownames(seu), 20).")
  }
}

# Normalize data for visualization (log1p of per-10k UMI)
# Note: if .X in the h5ad was already normalized, skip NormalizeData
# (Step 1 script sets .X to raw UMI counts; if it wasn't available it left .X as normalized)
DefaultAssay(seu) <- "originalexp"   # zellkonverter default assay name

# Check if data is raw or normalized by looking at max value
max_val <- max(GetAssayData(seu, layer = "counts")[GENE, ])
message("  Max value for ", GENE, " in counts slot: ", round(max_val, 2))

if (max_val > 50) {
  message("  Data looks like raw counts. Running NormalizeData...")
  seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
} else {
  message("  Data appears already normalized. Using as-is for visualization.")
  # Copy counts to data slot for Seurat visualization functions
  seu <- SetAssayData(seu, layer = "data",
                      new.data = GetAssayData(seu, layer = "counts"))
}

# Expression summary
expr_vec <- GetAssayData(seu, layer = "data")[GENE, ]
n_expr   <- sum(expr_vec > 0)
pct_expr <- n_expr / ncol(seu) * 100

message("  ", GENE, " detected in ", n_expr, " / ", ncol(seu), " cells (",
        round(pct_expr, 2), "%)")
message("  Mean expression (in all cells): ", round(mean(expr_vec), 4))
message("  Mean expression (in expressing cells): ",
        round(mean(expr_vec[expr_vec > 0]), 4))


# =============================================================================
# 5. PER-CELL-TYPE EXPRESSION SUMMARY TABLE
# =============================================================================

message("\n[4/7] Computing per-cell-type expression statistics...")

meta <- seu@meta.data
meta$TRAPPC9_expr <- as.numeric(expr_vec[rownames(meta)])

# Function: compute summary stats for a grouping
summarise_expr <- function(df, group_col) {
  df %>%
    group_by(across(all_of(group_col))) %>%
    summarise(
      n_cells       = n(),
      n_expressing  = sum(TRAPPC9_expr > 0),
      pct_expressing = round(n_expressing / n_cells * 100, 2),
      mean_expr_all  = round(mean(TRAPPC9_expr), 4),
      mean_expr_pos  = round(mean(TRAPPC9_expr[TRAPPC9_expr > 0], na.rm = TRUE), 4),
      median_expr    = round(median(TRAPPC9_expr), 4),
      .groups = "drop"
    ) %>%
    arrange(desc(pct_expressing))
}

# By cell type only
stats_by_celltype <- summarise_expr(meta, CELLTYPE_COL)
message("  Expression by ", CELLTYPE_COL, ":")
print(stats_by_celltype)

# By cell type × ADNC (the main comparison)
stats_by_ct_adnc <- summarise_expr(meta, c(CELLTYPE_COL, "ADNC"))

# By class × ADNC (broad overview)
if ("Class" %in% colnames(meta)) {
  stats_by_class_adnc <- summarise_expr(meta, c("Class", "ADNC"))
}

# Wilcoxon test: Control (Not AD) vs AD (High) per cell type
message("  Running Wilcoxon tests (Control vs High AD) per cell type...")

cell_types <- unique(meta[[CELLTYPE_COL]])
wilcox_results <- lapply(cell_types, function(ct) {
  sub_ctrl <- meta %>% filter(.data[[CELLTYPE_COL]] == ct, ADNC == "Not AD")
  sub_ad   <- meta %>% filter(.data[[CELLTYPE_COL]] == ct, ADNC == "High")

  if (nrow(sub_ctrl) < MIN_CELLS_FOR_STATS || nrow(sub_ad) < MIN_CELLS_FOR_STATS) {
    return(tibble(
      Cell_Type       = ct,
      n_Control       = nrow(sub_ctrl),
      n_AD            = nrow(sub_ad),
      mean_Control    = NA_real_,
      mean_AD         = NA_real_,
      pct_Control     = NA_real_,
      pct_AD          = NA_real_,
      log2FC          = NA_real_,
      p_value         = NA_real_,
      note            = "Too few cells"
    ))
  }

  wt <- wilcox.test(sub_ad$TRAPPC9_expr, sub_ctrl$TRAPPC9_expr,
                    alternative = "two.sided", exact = FALSE)

  mean_ctrl <- mean(sub_ctrl$TRAPPC9_expr)
  mean_ad   <- mean(sub_ad$TRAPPC9_expr)
  pct_ctrl  <- mean(sub_ctrl$TRAPPC9_expr > 0) * 100
  pct_ad    <- mean(sub_ad$TRAPPC9_expr > 0) * 100

  tibble(
    Cell_Type    = ct,
    n_Control    = nrow(sub_ctrl),
    n_AD         = nrow(sub_ad),
    mean_Control = round(mean_ctrl, 4),
    mean_AD      = round(mean_ad, 4),
    pct_Control  = round(pct_ctrl, 2),
    pct_AD       = round(pct_ad, 2),
    log2FC       = round(log2(mean_ad + 1e-6) - log2(mean_ctrl + 1e-6), 4),
    p_value      = wt$p.value,
    note         = "OK"
  )
}) %>% bind_rows()

# Multiple-testing correction (Benjamini-Hochberg)
ok_idx <- which(wilcox_results$note == "OK")
if (length(ok_idx) > 0) {
  wilcox_results$p_adj_BH <- NA_real_
  wilcox_results$p_adj_BH[ok_idx] <- p.adjust(wilcox_results$p_value[ok_idx], method = "BH")
  wilcox_results$significant <- wilcox_results$p_adj_BH < 0.05
}

wilcox_results <- wilcox_results %>% arrange(p_value)
message("  Wilcoxon results (significant cell types):")
print(wilcox_results %>% filter(significant == TRUE | is.na(significant)))

# Save stats tables
write_csv(wilcox_results,
          file.path(OUTPUT_DIR, paste0(GENE, "_stats_per_celltype.csv")))
write_csv(stats_by_ct_adnc,
          file.path(OUTPUT_DIR, paste0(GENE, "_expression_by_celltype_ADNC.csv")))


# =============================================================================
# 6. VISUALIZATIONS
# =============================================================================

message("\n[5/7] Generating visualizations...")

# ----- Palette helpers -------------------------------------------------------
# Use at most 12 qualitative colors for cell types
n_ct <- length(unique(meta[[CELLTYPE_COL]]))
ct_colors <- setNames(
  colorRampPalette(brewer.pal(min(n_ct, 12), "Paired"))(n_ct),
  sort(unique(meta[[CELLTYPE_COL]]))
)


# ---- 6a. UMAP colored by TRAPPC9 expression (if UMAP available) ------------
if ("umap" %in% names(seu@reductions)) {
  message("  Plotting UMAP...")
  p_umap_gene <- FeaturePlot(
    seu,
    features   = GENE,
    reduction  = "umap",
    order      = TRUE,
    pt.size    = 0.3,
    raster     = FALSE,
    cols       = c("lightgrey", "#D6604D")
  ) +
    ggtitle(paste(GENE, "expression —", REGION_LABEL, "(SEA-AD)")) +
    theme_classic(base_size = 12)

  p_umap_ct <- DimPlot(
    seu,
    reduction  = "umap",
    group.by   = CELLTYPE_COL,
    pt.size    = 0.3,
    label      = TRUE,
    repel      = TRUE,
    raster     = FALSE,
    cols       = ct_colors
  ) + ggtitle("Cell types") + theme_classic(base_size = 12) +
    NoLegend()

  p_umap_adnc <- DimPlot(
    seu,
    reduction = "umap",
    group.by  = "ADNC",
    pt.size   = 0.3,
    raster    = FALSE,
    cols      = ADNC_COLORS
  ) + ggtitle("ADNC status") + theme_classic(base_size = 12)

  p_umap_combined <- (p_umap_gene | p_umap_ct | p_umap_adnc)

  ggsave(file.path(OUTPUT_DIR, paste0(GENE, "_umap.pdf")),
         p_umap_combined, width = 18, height = 6)
  message("  Saved UMAP plot.")
} else {
  message("  No UMAP reduction found — skipping UMAP plot.")
}


# ---- 6b. DotPlot: cell type × ADNC -----------------------------------------
message("  Plotting DotPlot...")

# Subset to Control and High AD only for cleaner contrast plot
seu_contrast <- subset(seu, ADNC %in% c("Not AD", "High"))

# Create a combined grouping identity: CellType|ADNC
seu_contrast$CT_ADNC <- paste0(
  seu_contrast@meta.data[[CELLTYPE_COL]], " | ",
  ifelse(seu_contrast$ADNC == "Not AD", "Control", "AD")
)
# Order: for each cell type, Control then AD
ct_levels   <- sort(unique(seu_contrast@meta.data[[CELLTYPE_COL]]))
dot_order   <- as.vector(rbind(
  paste0(ct_levels, " | Control"),
  paste0(ct_levels, " | AD")
))
dot_order   <- dot_order[dot_order %in% unique(seu_contrast$CT_ADNC)]
seu_contrast$CT_ADNC <- factor(seu_contrast$CT_ADNC, levels = rev(dot_order))
Idents(seu_contrast) <- "CT_ADNC"

p_dot <- DotPlot(
  seu_contrast,
  features     = GENE,
  cols         = c("lightgrey", "#D6604D"),
  dot.scale    = 8,
  scale        = TRUE
) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold")
  ) +
  ggtitle(paste(GENE, "— % expressing & mean expression\nby cell type and ADNC status")) +
  xlab("") + ylab("Cell type | Disease group")

ggsave(file.path(OUTPUT_DIR, paste0(GENE, "_dotplot_celltype_x_ADNC.pdf")),
       p_dot, width = 10, height = max(6, n_ct * 0.6 + 2))
message("  Saved DotPlot.")


# ---- 6c. Violin plot: expression by cell type, split by ADNC ---------------
message("  Plotting Violin plots...")

Idents(seu) <- CELLTYPE_COL

# Full 4-group violin (all ADNC levels)
p_vln_full <- VlnPlot(
  seu,
  features    = GENE,
  group.by    = CELLTYPE_COL,
  split.by    = "ADNC",
  split.plot  = TRUE,
  pt.size     = 0,
  cols        = ADNC_COLORS
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  ggtitle(paste(GENE, "expression by cell type and ADNC —", REGION_LABEL, "(SEA-AD)"))

ggsave(file.path(OUTPUT_DIR, paste0(GENE, "_violin_by_subclass.pdf")),
       p_vln_full, width = max(10, n_ct * 0.9 + 2), height = 6)

# Focused Control vs High AD violin, per cell type
p_vln_contrast <- VlnPlot(
  seu_contrast,
  features    = GENE,
  group.by    = CELLTYPE_COL,
  split.by    = "ADNC",
  split.plot  = TRUE,
  pt.size     = 0,
  cols        = ADNC_COLORS[c("Not AD", "High")]
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  ggtitle(paste(GENE, "— Control vs High AD"))

ggsave(file.path(OUTPUT_DIR, paste0(GENE, "_violin_control_vs_AD.pdf")),
       p_vln_contrast, width = max(10, n_ct * 0.9 + 2), height = 6)
message("  Saved Violin plots.")


# ---- 6d. Bar chart: % TRAPPC9-expressing cells by cell type × ADNC ---------
message("  Plotting % expressing bar chart...")

pct_df <- stats_by_ct_adnc %>%
  filter(!is.na(.data[[CELLTYPE_COL]])) %>%
  mutate(ADNC = factor(ADNC, levels = c("Not AD", "Low", "Intermediate", "High")))

colnames(pct_df)[colnames(pct_df) == CELLTYPE_COL] <- "Cell_Type"

p_bar <- ggplot(pct_df, aes(x = Cell_Type, y = pct_expressing, fill = ADNC)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = ADNC_COLORS, name = "ADNC") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "top"
  ) +
  labs(
    title = paste("% cells expressing", GENE, "by cell type and ADNC"),
    subtitle = paste0("SEA-AD — ", REGION_LABEL),
    x = "Cell type", y = "% TRAPPC9+ cells"
  )

ggsave(file.path(OUTPUT_DIR, paste0(GENE, "_pct_expressing_barplot.pdf")),
       p_bar, width = max(8, n_ct * 0.7 + 2), height = 6)
message("  Saved bar chart.")


# ---- 6e. Log2 fold-change plot (Control vs High AD per cell type) ----------
message("  Plotting log2FC plot...")

fc_df <- wilcox_results %>%
  filter(note == "OK", !is.na(log2FC)) %>%
  mutate(
    direction  = ifelse(log2FC > 0, "Up in AD", "Down in AD"),
    sig_label  = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE             ~ ""
    ),
    Cell_Type = factor(Cell_Type, levels = Cell_Type[order(log2FC)])
  )

p_fc <- ggplot(fc_df, aes(x = Cell_Type, y = log2FC, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sig_label,
                y = log2FC + ifelse(log2FC > 0, 0.05, -0.05),
                vjust = ifelse(log2FC > 0, 0, 1)),
            size = 5) +
  scale_fill_manual(values = c("Up in AD" = "#D6604D", "Down in AD" = "#2166AC"),
                    name = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(legend.position = "top") +
  labs(
    title    = paste(GENE, "log2 fold-change: AD (High) vs Control (Not AD)"),
    subtitle = "* p_adj<0.05  ** p_adj<0.01  *** p_adj<0.001  (Wilcoxon, BH-corrected)",
    x        = NULL,
    y        = "log2 FC (AD vs Control)"
  )

ggsave(file.path(OUTPUT_DIR, paste0(GENE, "_log2FC_per_celltype.pdf")),
       p_fc, width = 8, height = max(5, nrow(fc_df) * 0.4 + 2))
message("  Saved log2FC plot.")


# ---- 6f. Summary heatmap: mean expression by cell type × ADNC --------------
message("  Plotting summary heatmap...")

# Pivot to matrix
heat_df <- stats_by_ct_adnc %>%
  select(all_of(c(CELLTYPE_COL, "ADNC", "mean_expr_all"))) %>%
  pivot_wider(names_from = ADNC, values_from = mean_expr_all) %>%
  column_to_rownames(CELLTYPE_COL) %>%
  select(any_of(c("Not AD", "Low", "Intermediate", "High")))   # order cols

# Replace NA with 0 (cell types with no cells in that ADNC group)
heat_mat <- as.matrix(heat_df)
heat_mat[is.na(heat_mat)] <- 0

if (nrow(heat_mat) > 1) {
  pdf(file.path(OUTPUT_DIR, paste0(GENE, "_summary_heatmap.pdf")),
      width = 7, height = max(4, nrow(heat_mat) * 0.35 + 2))
  pheatmap(
    heat_mat,
    color            = colorRampPalette(c("white", "#FEE0D2", "#D6604D"))(50),
    cluster_cols     = FALSE,
    cluster_rows     = TRUE,
    border_color     = NA,
    main             = paste(GENE, "mean expression (log-normalized)\nby cell type and ADNC"),
    fontsize_row     = 9,
    fontsize_col     = 11,
    angle_col        = 45
  )
  dev.off()
  message("  Saved heatmap.")
} else {
  message("  Only one cell type — skipping heatmap.")
}


# =============================================================================
# 7. PSEUDOBULK DIFFERENTIAL EXPRESSION
# =============================================================================
# Pseudobulk: aggregate raw counts per donor × cell type, then run DESeq2.
# This is the most statistically rigorous approach for scRNA-seq DE.

message("\n[6/7] Pseudobulk DE analysis (Control vs High AD)...")

if (requireNamespace("DESeq2", quietly = TRUE)) {
  library(DESeq2)
  library(SummarizedExperiment)

  # Subset to Control and High AD
  seu_pb <- subset(seu, ADNC %in% c("Not AD", "High"))

  # Check we have a donor column
  donor_col <- intersect(c("Donor_ID", "Donor.ID", "donor_id"), colnames(seu_pb@meta.data))
  if (length(donor_col) == 0) {
    message("  WARNING: Donor ID column not found — skipping pseudobulk.")
  } else {
    donor_col <- donor_col[1]

    # Aggregate counts per donor × cell type
    pb_list <- lapply(unique(seu_pb@meta.data[[CELLTYPE_COL]]), function(ct) {
      ct_cells <- subset(seu_pb, cells = colnames(seu_pb)[seu_pb@meta.data[[CELLTYPE_COL]] == ct])
      if (ncol(ct_cells) < MIN_CELLS_FOR_STATS) return(NULL)

      donors <- unique(ct_cells@meta.data[[donor_col]])
      count_mat <- sapply(donors, function(d) {
        d_cells <- subset(ct_cells, cells = colnames(ct_cells)[ct_cells@meta.data[[donor_col]] == d])
        if (ncol(d_cells) == 0) return(NULL)
        rowSums(GetAssayData(d_cells, layer = "counts"))
      })
      if (is.null(count_mat) || length(count_mat) == 0) return(NULL)

      donor_adnc <- sapply(donors, function(d)
        as.character(unique(ct_cells@meta.data$ADNC[ct_cells@meta.data[[donor_col]] == d])[1]))
      meta_pb <- data.frame(
        donor     = donors,
        ADNC      = donor_adnc,
        row.names = donors
      )
      meta_pb$group <- ifelse(meta_pb$ADNC == "Not AD", "Control", "AD")

      list(counts = count_mat, meta = meta_pb, cell_type = ct)
    })
    pb_list <- Filter(Negate(is.null), pb_list)

    # Run DESeq2 for TRAPPC9 per cell type
    de_results <- lapply(pb_list, function(pb) {
      if (ncol(pb$counts) < 4) return(NULL)  # need enough samples
      # Check TRAPPC9 is in this subset
      if (!GENE %in% rownames(pb$counts)) return(NULL)

      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = round(pb$counts),
          colData   = pb$meta,
          design    = ~ group
        )
        dds <- dds[rowSums(counts(dds)) >= 10, ]
        if (!GENE %in% rownames(dds)) return(NULL)

        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds, contrast = c("group", "AD", "Control"),
                       alpha = 0.05)

        if (GENE %in% rownames(res)) {
          as.data.frame(res[GENE, ]) %>%
            mutate(Cell_Type = pb$cell_type, Gene = GENE,
                   n_donors_ctrl = sum(pb$meta$group == "Control"),
                   n_donors_AD   = sum(pb$meta$group == "AD"))
        } else {
          NULL
        }
      }, error = function(e) {
        message("  DESeq2 failed for ", pb$cell_type, ": ", conditionMessage(e))
        NULL
      })
    })

    de_results <- bind_rows(de_results) %>%
      select(Cell_Type, Gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
             n_donors_ctrl, n_donors_AD) %>%
      arrange(padj)

    if (nrow(de_results) > 0) {
      write_csv(de_results,
                file.path(OUTPUT_DIR, paste0(GENE, "_pseudobulk_DESeq2.csv")))
      message("  DESeq2 pseudobulk results:")
      print(de_results)
      message("  Saved pseudobulk DE results.")
    } else {
      message("  No DESeq2 results generated (possibly too few donors per cell type).")
    }
  }
} else {
  message("  DESeq2 not installed — skipping pseudobulk analysis.")
  message("  Install with: BiocManager::install('DESeq2')")
}


# =============================================================================
# 8. FINAL SUMMARY
# =============================================================================

message("\n[7/7] Summary")
message(paste0(strrep("=", 61)))
message("Gene: ", GENE)
message("Region: ", REGION_LABEL, " (SEA-AD)")
message("Total cells analyzed: ", ncol(seu))
message("Cells expressing ", GENE, ": ", n_expr, " (", round(pct_expr, 2), "%)")

sig_cts <- wilcox_results %>% filter(significant == TRUE)
if (nrow(sig_cts) > 0) {
  message("\nCell types with significant TRAPPC9 difference (AD vs Control):")
  print(sig_cts %>% select(Cell_Type, n_Control, n_AD, log2FC, p_adj_BH))
} else {
  message("\nNo cell types showed significant TRAPPC9 differences after BH correction.")
  message("This may reflect low expression levels — consider the pseudobulk results.")
}

message("\nOutput files in: ", OUTPUT_DIR)
message("  ", GENE, "_stats_per_celltype.csv       — Wilcoxon stats per cell type")
message("  ", GENE, "_expression_by_celltype_ADNC.csv — expression summary table")
message("  ", GENE, "_dotplot_celltype_x_ADNC.pdf  — DotPlot")
message("  ", GENE, "_violin_by_subclass.pdf        — Violin (all ADNC groups)")
message("  ", GENE, "_violin_control_vs_AD.pdf      — Violin (Control vs AD)")
message("  ", GENE, "_pct_expressing_barplot.pdf    — % expressing cells")
message("  ", GENE, "_log2FC_per_celltype.pdf       — log2FC barplot")
message("  ", GENE, "_summary_heatmap.pdf           — Mean expression heatmap")
message("  ", GENE, "_pseudobulk_DESeq2.csv         — Pseudobulk DE (if DESeq2 available)")
message("\nDone — ", region, "!")

# Clean up large objects before next region
rm(seu, meta); if (exists("seu_contrast")) rm(seu_contrast)
if (exists("seu_pb")) rm(seu_pb)
gc()

} # end for (region in REGIONS_TO_RUN)

message("\n", strrep("=", 65))
message("ALL REGIONS COMPLETE.")
message(strrep("=", 65))
