# SEA-AD TRAPPC9 Analysis

Exploratory analysis of TRAPPC9 expression across Alzheimer's disease severity in the SEA-AD hippocampus single-nucleus RNA-seq dataset.

## Background

TRAPPC9 was discovered/characterized in part by Dr. Wenhui Hu's lab and is of interest for a potential role in Alzheimer's disease pathogenesis. This is an exploratory side project investigating whether TRAPPC9 expression varies with AD neuropathological severity (ADNC) in microglia and astrocytes.

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

### Additional analysis scripts
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

**R:** Seurat, zellkonverter, SingleCellExperiment, rhdf5, DESeq2, ggplot2, ggpubr, dplyr, clinfun

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
