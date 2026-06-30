# SEA-AD TRAPPC9 Analysis

Exploratory analysis of TRAPPC9 expression across Alzheimer's disease severity and brain regions in the SEA-AD single-nucleus RNA-seq dataset, with genetic validation via MAGMA cell type enrichment.

## Background

TRAPPC9 was discovered/characterized in part by Dr. Wenhui Hu's lab and is of interest for a potential role in Alzheimer's disease pathogenesis. This is an exploratory side project investigating whether TRAPPC9 expression varies with AD neuropathological severity (ADNC) in microglia and astrocytes across multiple brain regions.

Our analysis has two complementary aims:

1. **Expression analysis (Steps 1–2):** Does TRAPPC9 expression in microglia and astrocytes change with increasing AD neuropathological severity (ADNC), and is this consistent across brain regions?

2. **Genetic validation (Step 3):** Are the cell types showing TRAPPC9 expression changes also genetically implicated in AD risk? We address this using MAGMA cell type enrichment — a method that tests whether genes specifically expressed in a given cell type (e.g. microglia) are genome-wide enriched for variants associated with AD in a large GWAS. This is **not** a test of whether TRAPPC9 itself is a GWAS hit; rather, it asks whether the *cell types* we are focused on are the right biological context for AD genetics. If microglia and/or astrocytes show significant MAGMA enrichment for AD risk, it provides independent genetic evidence that those cell types are centrally involved in AD pathogenesis — strengthening the rationale for studying TRAPPC9 in those specific populations.

## Dataset

[Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)](https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html) snRNA-seq dataset (Gabitto et al., 2024, *Nature Neuroscience*).

- 40 donors spanning the full AD neuropathological spectrum
- 208,843 total nuclei in the hippocampus alone, subsetted to 38,026 immune (microglia) and astrocyte nuclei
- Donors classified by ADNC (Overall AD Neuropathological Change): Not AD / Low / Intermediate / High
- Originally hippocampus-only; the download pipeline now supports all 11 SEA-AD brain regions:

| Region code | Brain region |
|---|---|
| HIP | Hippocampus |
| MEC | Medial entorhinal cortex |
| LEC | Lateral entorhinal cortex |
| ITG | Inferior temporal gyrus |
| MTG | Middle temporal gyrus |
| STG | Superior temporal gyrus |
| FI | Frontal insula |
| PFC | Prefrontal cortex |
| AnG | Angular gyrus |
| V1C | Primary visual cortex |
| Caudate_Nucleus | Caudate nucleus (head of the caudate) |

## Pipeline

### Step 1: Download & preprocess (Python)
`step1_download_preprocess.py`
- Downloads the raw SEA-AD h5ad for each region in `REGIONS_TO_PROCESS` from AWS S3 (defaults to all 11 regions; trim the list to run a subset)
- Filters to ADNC groups of interest and Immune/Astrocyte subclasses
- Outputs a compact h5ad per region to `seaad_data/<region>/SEAAD_<region>_subset_for_R.h5ad` for downstream R analysis

### Step 2: Seurat analysis (R)
`step2_TRAPPC9_seurat_analysis.R`
- Loads h5ad into a Seurat object, decodes categorical metadata
- Normalizes TRAPPC9 expression
- Generates UMAP, dotplot, violin plots, bar charts, and summary heatmap
- Runs Wilcoxon tests (cell-level, Control vs High AD) and Jonckheere-Terpstra trend tests (across ADNC severity gradient)
- Runs donor-level pseudobulk DESeq2 differential expression

### Step 3: MAGMA cell type enrichment (R)
`step3_MAGMA_celltyping.R`

**What it tests:** Whether the cell types showing TRAPPC9 expression changes in AD (microglia and astrocytes) are also independently implicated in AD genetic risk. MAGMA cell type enrichment works by computing a *specificity score* for each gene in each cell type — how exclusively that gene is expressed in that cell type compared to all others — and then testing genome-wide whether the genes most specific to a given cell type are enriched for AD GWAS signal. A significant result for microglia, for example, means that AD risk variants cluster in genes that are preferentially expressed in microglia, which is genetic evidence that microglia are a key cell type in AD biology.

**Why this matters for TRAPPC9:** If microglia and astrocytes show significant MAGMA enrichment, it validates that our focus on TRAPPC9 in those cell types is biologically well-grounded. Conversely, finding enrichment in a cell type where TRAPPC9 expression does *not* change with AD severity would suggest we may be missing something.

**What it does:**
- Downloads the MAGMA binary automatically via `MAGMA.Celltyping::get_magma()`
- Downloads Bellenguez et al. 2022 AD GWAS summary statistics (GCST90027158, ~1 million variants, the largest AD GWAS to date) from the GWAS Catalog
- Maps GWAS SNPs to genes using MAGMA (1000 Genomes EUR reference panel, GRCh38)
- Builds a CellTypeData (CTD) specificity object for each SEA-AD brain region from the step 1 outputs
- Runs linear MAGMA enrichment tests per cell type across all 11 brain regions
- Outputs per-region results tables and bar chart PDFs, plus a combined cross-region heatmap in `results/MAGMA/`

### Additional analysis scripts (in `src/`)
- `step3_TRAPPC9_boxplot.R` — Control vs AD boxplots (pooled and by cell type)
- `step4_TRAPPC9_ADNC_gradient.R` — expression trends across the full ADNC severity gradient
- `step5_violin_with_points.R` — violin plots with overlaid jittered points and p-values
- `step6_pseudobulk_barplot.R` — visualization of donor-level pseudobulk log2FC results

## Key findings (preliminary)

- TRAPPC9 is broadly expressed: detected in ~64% of all immune/astrocyte nuclei
- Cell-level Wilcoxon tests show small but statistically significant differences between Control and High AD donors (astrocytes up, microglia down; both p_adj < 0.05)
- Jonckheere-Terpstra trend tests confirm a significant monotonic trend across ADNC severity in both cell types (astrocyte p = 1.7×10⁻⁸, immune p = 3.2×10⁻⁴)
- Donor-level pseudobulk DESeq2 (the more statistically rigorous test, n=4 control vs n=21 AD donors) does **not** reach significance after correction for either cell type, likely reflecting modest effect size and donor sample imbalance

These results should be interpreted as hypothesis-generating preliminary data rather than confirmed findings, given the donor-level non-significance.

## Requirements

**Python:** anndata, h5py, numpy, pandas, scipy (see `requirements.txt`)

**R (steps 1–2):** Seurat, zellkonverter, SingleCellExperiment, rhdf5, DESeq2, ggplot2, ggpubr, dplyr, clinfun

**R (step 3 — MAGMA):** MAGMA.Celltyping, EWCE, zellkonverter, SingleCellExperiment, ggplot2, dplyr, data.table, patchwork

**MAGMA binary:** downloaded automatically by `step3_MAGMA_celltyping.R` via `MAGMA.Celltyping::get_magma()`. Windows users can also manually download from https://ctg.cncr.nl/software/magma

## Usage

```bash
# Install Python dependencies
pip install -r requirements.txt

# Step 1 (Python)
python step1_download_preprocess.py

# Step 2 (R, in RStudio or terminal)
Rscript step2_TRAPPC9_seurat_analysis.R
```

Update `H5AD_PATH` and `OUTPUT_DIR` at the top of the R script to match your local file paths.

## Author

Jack Wallace, Hu Lab, Virginia Commonwealth University School of Medicine
