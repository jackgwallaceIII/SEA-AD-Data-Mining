# =============================================================================
# SEA-AD TRAPPC9 Analysis Pipeline — Step 3: MAGMA Cell Type Enrichment
# =============================================================================
# Dataset  : Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)
# GWAS     : Bellenguez et al. 2022 (Nature Genetics) — largest AD GWAS to date
# Question : Which cell types in each SEA-AD brain region show significant
#            genetic enrichment for Alzheimer's disease risk loci?
#
# Prerequisites:
#   1. Run step1_download_preprocess.py to generate per-region h5ad files
#      (seaad_data/<REGION>/SEAAD_<REGION>_subset_for_R.h5ad)
#   2. Install R packages listed below (script will attempt auto-install)
#
# What this script does:
#   1. Downloads the MAGMA binary (needed by MAGMA.Celltyping)
#   2. Downloads Bellenguez 2022 AD GWAS summary statistics
#   3. Maps GWAS SNPs to genes using MAGMA
#   4. Loads SEA-AD h5ad data and builds a CellTypeData (CTD) specificity object
#   5. Runs linear MAGMA cell type enrichment tests per cell type
#   6. Saves results tables and plots to results/MAGMA/
#
# Outputs (saved to results/MAGMA/):
#   - ctd/                    — CellTypeData objects (one per region)
#   - magma_genes.tsv         — gene-level MAGMA scores
#   - celltype_results_<REGION>.tsv — enrichment results per region
#   - celltype_enrichment_<REGION>.pdf — bar/dot plots
#
# MAGMA binary note:
#   MAGMA is a standalone executable (not an R package) from Radboud University.
#   The get_magma() call below downloads the correct binary for your OS.
#   It is saved to a local cache folder and reused on subsequent runs.
#   Windows users: if the download fails, manually download the Windows binary
#   from https://ctg.cncr.nl/software/magma and point MAGMA_PATH to it below.
# =============================================================================


# =============================================================================
# 0. CONFIGURATION
# =============================================================================

# Regions to analyse. Must match the keys used in step1_download_preprocess.py
# and have corresponding h5ad files in seaad_data/<REGION>/.
REGIONS_TO_ANALYSE <- c(
  "HIP", "MEC", "LEC", "ITG", "MTG", "STG",
  "FI", "PFC", "AnG", "V1C", "Caudate_Nucleus"
)

# Cell taxonomy level to use for enrichment testing.
# "Subclass"  — broad cell types (Astrocyte, Immune/Microglia, Oligodendrocyte …)
# "Supertype" — finer subtypes
CELLTYPE_LEVEL <- "Subclass"

# Genome build used by the SEA-AD dataset
GENOME_BUILD <- "GRCh38"

# Output directory
OUTPUT_DIR <- "results/MAGMA"

# Paths (populated automatically; override here if needed)
H5AD_DIR   <- "seaad_data"
MAGMA_PATH <- NULL   # NULL = auto-download via get_magma()


# =============================================================================
# 1. PACKAGES
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("dplyr", "ggplot2", "data.table", "stringr", "patchwork",
               "remotes")
bioc_pkgs <- c("zellkonverter", "SingleCellExperiment", "EWCE")

install_if_missing <- function(pkgs, installer = install.packages) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "))
    installer(missing, ask = FALSE)
  }
}

install_if_missing(cran_pkgs)
install_if_missing(bioc_pkgs, installer = function(p, ...) BiocManager::install(p, ...))

# MAGMA.Celltyping is on Bioconductor / GitHub
if (!"MAGMA.Celltyping" %in% installed.packages()[, "Package"]) {
  message("Installing MAGMA.Celltyping ...")
  BiocManager::install("MAGMA.Celltyping")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(patchwork)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(EWCE)
  library(MAGMA.Celltyping)
})

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "ctd"), recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# 2. DOWNLOAD MAGMA BINARY
# =============================================================================

message("\n=== Step 3a: Downloading MAGMA binary ===")

if (is.null(MAGMA_PATH)) {
  MAGMA_PATH <- MAGMA.Celltyping::get_magma()
}
message("  MAGMA binary: ", MAGMA_PATH)


# =============================================================================
# 3. DOWNLOAD BELLENGUEZ 2022 AD GWAS & MAP SNPS TO GENES
# =============================================================================
# Bellenguez et al. 2022, Nature Genetics
# "New insights into the genetic etiology of Alzheimer's disease and related
# dementias" — GWAS Catalog accession GCST90027158 / OpenGWAS ID ieu-b-5067
# =============================================================================

message("\n=== Step 3b: Obtaining Bellenguez 2022 GWAS summary statistics ===")

gwas_file <- file.path(OUTPUT_DIR, "bellenguez2022_AD.tsv.gz")

if (!file.exists(gwas_file)) {
  message("  Downloading Bellenguez 2022 GWAS from GWAS Catalog ...")
  # Harmonised summary statistics from GWAS Catalog
  gwas_url <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
    "GCST90027001-GCST90028000/GCST90027158/harmonised/",
    "GCST90027158.h.tsv.gz"
  )
  tryCatch(
    download.file(gwas_url, gwas_file, mode = "wb"),
    error = function(e) {
      stop(
        "GWAS download failed. You can manually download the Bellenguez 2022 ",
        "harmonised summary statistics from:\n",
        "  https://www.ebi.ac.uk/gwas/studies/GCST90027158\n",
        "(click 'Download' -> 'Harmonised summary statistics')\n",
        "and save the file to: ", gwas_file
      )
    }
  )
} else {
  message("  GWAS file already exists: ", gwas_file)
}

message("\n=== Step 3c: Mapping GWAS SNPs to genes with MAGMA ===")

magma_snp_loc_file <- file.path(OUTPUT_DIR, "magma_genes")

# get_genome_ref() downloads the 1000 Genomes EUR reference panel that MAGMA
# needs for LD-based gene scoring.
genome_ref_path <- MAGMA.Celltyping::get_genome_ref(
  genome_build = GENOME_BUILD
)
message("  Genome reference: ", genome_ref_path)

# map_snps_to_genes() runs MAGMA annotation + gene-level analysis in two steps.
# It produces <magma_genes>.genes.raw and <magma_genes>.genes.out files.
# This can take several minutes depending on your machine.
magma_out <- MAGMA.Celltyping::map_snps_to_genes(
  path_formatted   = gwas_file,
  genome_ref_path  = genome_ref_path,
  genome_build     = GENOME_BUILD,
  upstream_kb      = 10,
  downstream_kb    = 1.5,
  magma_path       = MAGMA_PATH
)
message("  MAGMA gene file: ", magma_out)


# =============================================================================
# 4. HELPER: BUILD CTD FROM A SEA-AD H5AD FILE
# =============================================================================
# CTD = CellTypeData object (from the EWCE package).
# It stores mean expression and specificity scores for each cell type,
# which MAGMA.Celltyping uses to run enrichment tests.
# =============================================================================

build_ctd <- function(region_key, celltype_col = CELLTYPE_LEVEL) {
  h5ad_path <- file.path(H5AD_DIR, region_key,
                          paste0("SEAAD_", region_key, "_subset_for_R.h5ad"))

  if (!file.exists(h5ad_path)) {
    warning("h5ad not found for ", region_key, ": ", h5ad_path,
            "\n  Run step1_download_preprocess.py first.")
    return(NULL)
  }

  message("  Loading ", region_key, " ...")
  sce <- zellkonverter::readH5AD(h5ad_path, use_hdf5 = TRUE, verbose = FALSE)

  if (!celltype_col %in% colnames(colData(sce))) {
    warning("Column '", celltype_col, "' not found in ", region_key,
            ". Available: ", paste(colnames(colData(sce)), collapse = ", "))
    return(NULL)
  }

  # Drop cells with NA cell type labels
  keep <- !is.na(colData(sce)[[celltype_col]])
  sce  <- sce[, keep]

  # EWCE requires a normalised expression matrix (log-normalised counts).
  # If .X contains raw UMIs (from step1), log-normalise here.
  mat <- assay(sce, "X")
  if (max(mat) > 100) {
    # Raw counts detected — log-normalise (library-size normalise then log1p)
    message("    Detected raw counts; applying library-size normalisation ...")
    lib_sizes <- Matrix::colSums(mat)
    mat       <- Matrix::t(Matrix::t(mat) / lib_sizes) * 1e4
    mat       <- log1p(mat)
  }

  ann_df <- data.frame(
    cell_id   = colnames(sce),
    cell_type = as.character(colData(sce)[[celltype_col]]),
    stringsAsFactors = FALSE
  )

  # generate_celltype_data() computes mean expression and specificity per type
  ctd <- EWCE::generate_celltype_data(
    exp          = mat,
    annotLevels  = list(level1 = ann_df$cell_type),
    groupName    = region_key,
    savePath     = file.path(OUTPUT_DIR, "ctd"),
    verbose      = FALSE
  )

  return(ctd)
}


# =============================================================================
# 5. BUILD CTD OBJECTS FOR ALL REGIONS
# =============================================================================

message("\n=== Step 3d: Building cell type specificity objects (CTDs) ===")

ctd_list <- list()
for (region in REGIONS_TO_ANALYSE) {
  message("\n  Processing region: ", region)
  ctd_list[[region]] <- tryCatch(
    build_ctd(region),
    error = function(e) {
      message("  ERROR in ", region, ": ", conditionMessage(e))
      NULL
    }
  )
}

# Remove any regions that failed to load
ctd_list <- Filter(Negate(is.null), ctd_list)
message("\n  Successfully built CTDs for: ", paste(names(ctd_list), collapse = ", "))


# =============================================================================
# 6. RUN MAGMA CELL TYPE ENRICHMENT TESTS
# =============================================================================

message("\n=== Step 3e: Running MAGMA cell type enrichment tests ===")

all_results <- list()

for (region in names(ctd_list)) {
  message("\n  Running enrichment for: ", region)
  ctd <- ctd_list[[region]]

  res <- tryCatch(
    MAGMA.Celltyping::celltype_associations_pipeline(
      magma_sumstats_path = magma_out,
      ctd                 = ctd,
      ctd_name            = region,
      ctd_species         = "human",
      ctd_lvl             = 1,
      magma_path          = MAGMA_PATH,
      output_path         = file.path(OUTPUT_DIR, region)
    ),
    error = function(e) {
      message("  ERROR: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(res)) {
    res$Region <- region
    all_results[[region]] <- res

    out_path <- file.path(OUTPUT_DIR, paste0("celltype_results_", region, ".tsv"))
    write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("  Results saved: ", out_path)
  }
}


# =============================================================================
# 7. COMBINE RESULTS AND PLOT
# =============================================================================

message("\n=== Step 3f: Combining results and generating plots ===")

if (length(all_results) == 0) {
  stop("No enrichment results produced. Check errors above.")
}

combined <- bind_rows(all_results)

# Save combined results table
combined_path <- file.path(OUTPUT_DIR, "celltype_enrichment_all_regions.tsv")
write.table(combined, combined_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("  Combined results: ", combined_path)

# ---- Plot: -log10(p) per cell type, faceted by region ----------------------

plot_df <- combined %>%
  filter(!is.na(P)) %>%
  mutate(
    neg_log10_p = -log10(P),
    Significant = P < 0.05 / n()                   # Bonferroni across all tests
  )

p_heatmap <- ggplot(plot_df,
    aes(x = Region, y = Celltype, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(Significant, "*", "")),
            size = 4, vjust = 0.8) +
  scale_fill_gradient(low = "white", high = "#b22222",
                      name = expression(-log[10](p))) +
  labs(title = paste0("MAGMA Cell Type Enrichment — Bellenguez 2022 AD GWAS"),
       subtitle = paste0("Cell type level: ", CELLTYPE_LEVEL,
                         "  |  * = Bonferroni-significant"),
       x = "Brain Region", y = "Cell Type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid   = element_blank())

heatmap_path <- file.path(OUTPUT_DIR, "celltype_enrichment_heatmap.pdf")
ggsave(heatmap_path, p_heatmap, width = max(8, length(all_results) * 1.2),
       height = 7)
message("  Heatmap saved: ", heatmap_path)

# ---- Per-region bar plots --------------------------------------------------
for (region in names(all_results)) {
  df <- all_results[[region]] %>%
    filter(!is.na(P)) %>%
    mutate(neg_log10_p  = -log10(P),
           Significant  = P < 0.05,
           Celltype     = reorder(Celltype, neg_log10_p))

  p_bar <- ggplot(df, aes(x = Celltype, y = neg_log10_p,
                           fill = Significant)) +
    geom_col() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey50") +
    scale_fill_manual(values = c("TRUE" = "#b22222", "FALSE" = "#aaaaaa"),
                      name = "p < 0.05") +
    coord_flip() +
    labs(title = paste0("MAGMA Enrichment — ", region),
         subtitle = "Bellenguez 2022 AD GWAS",
         x = NULL, y = expression(-log[10](p))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, paste0("celltype_enrichment_", region, ".pdf")),
         p_bar, width = 6, height = max(4, nrow(df) * 0.4))
}

message("\n", paste(rep("=", 70), collapse = ""))
message("Step 3 complete!")
message("  Output directory : ", OUTPUT_DIR)
message("  Combined results : ", combined_path)
message("  Heatmap          : ", heatmap_path)
message("  Next: inspect ", OUTPUT_DIR, " for per-region PDFs and TSV tables.")
message(paste(rep("=", 70), collapse = ""))
