#!/usr/bin/env python3
"""
=============================================================================
SEA-AD TRAPPC9 Analysis Pipeline — Step 1: Download & Preprocess
=============================================================================
Dataset : Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)
Region  : Hippocampus
Source  : s3://sea-ad-single-cell-profiling/ (public, no AWS account needed)

What this script does:
  1. Downloads the SEA-AD Hippocampus snRNA-seq h5ad from AWS S3
  2. Removes "Reference" donors (5 young donors used only for cell type mapping)
  3. Subsets the AnnData to cells of interest and key metadata columns
  4. Saves a compact h5ad that R can load without the 2^32 integer overflow issue

IMPORTANT: The full dataset has >2^32 non-zero matrix elements, which exceeds
R's 32-bit integer limit. This step produces a smaller file that loads fine in R.
(See: https://community.brain-map.org/t/sea-ad-snrnaseq-dataset-format/3454)

Requirements (install with pip):
    pip install scanpy anndata numpy pandas scipy

Run:
    python step1_download_preprocess.py

Output:
    seaad_data/SEAAD_HIP_subset_for_R.h5ad
=============================================================================
"""

import os
import sys
import urllib.request
import urllib.error

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# CONFIGURATION — adjust if needed
# ---------------------------------------------------------------------------

REGION      = "Hippocampus"
REGION_CODE = "HIP"
OUTPUT_DIR  = "seaad_data"

# The SEA-AD team releases updated h5ad files periodically.
# The script tries these dates in order until it finds a valid file.
# Browse https://sea-ad-single-cell-profiling.s3.amazonaws.com/index.html#Hippocampus/RNAseq/
# to find the latest filename and add it here if needed.
CANDIDATE_DATES = [
    "2026-06-22",
    "2024-02-13",
    "2023-05-05",
    "2022-08-18",
]

BASE_URL = (
    f"https://sea-ad-single-cell-profiling.s3.amazonaws.com"
    f"/{REGION_CODE}/RNAseq"
)

# Metadata columns to keep (everything else is dropped to reduce file size)
# Rename map: {original_col_name: short_name_for_R}
KEEP_COLS = {
    "Overall AD neuropathological Change": "ADNC",
    "Supertype":                           "Supertype",
    "Subclass":                            "Subclass",
    "Class":                               "Class",
    "Braak":                               "Braak_Stage",
    "Thal":                                "Thal_Phase",
    "CERAD score":                         "CERAD_score",
    "Cognitive Status":                    "Cognitive_Status",
    "Donor ID":                            "Donor_ID",
    "Sex":                                 "Sex",
    "Age at Death":                        "Age_at_Death",
    "library_prep":                        "library_prep",
}

# ADNC groups to KEEP (exclude "Reference" — young mapping donors, not part of cohort)
ADNC_KEEP = ["Not AD", "Low", "Intermediate", "High"]

GENE_OF_INTEREST = "TRAPPC9"


# ---------------------------------------------------------------------------
# HELPER FUNCTIONS
# ---------------------------------------------------------------------------

def url_exists(url: str, timeout: int = 15) -> bool:
    """Return True if the URL responds with a 200 status."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout):
            return True
    except (urllib.error.URLError, urllib.error.HTTPError):
        return False


def find_download_url() -> str:
    """Try each candidate date and return the first working URL."""
    for date in CANDIDATE_DATES:
        fname = f"SEAAD_{REGION_CODE}_RNAseq_final-nuclei.{date}.h5ad"
        url   = f"{BASE_URL}/{fname}"
        print(f"  Checking {url} ...", end=" ", flush=True)
        if url_exists(url):
            print("FOUND")
            return url, fname
        print("not found")
    raise FileNotFoundError(
        "\n\nCould not find a valid Hippocampus h5ad on S3.\n"
        "Browse the bucket index to find the correct filename:\n"
        "  https://sea-ad-single-cell-profiling.s3.amazonaws.com"
        "/index.html#Hippocampus/RNAseq/\n"
        "Then add that date string to CANDIDATE_DATES at the top of this script."
    )


def download_with_progress(url: str, dest: str) -> None:
    """Download a file with a simple progress display."""
    def _hook(count, block_size, total_size):
        pct = count * block_size / total_size * 100
        downloaded_mb = count * block_size / 1e6
        total_mb = total_size / 1e6
        sys.stdout.write(
            f"\r  {downloaded_mb:.0f} / {total_mb:.0f} MB  ({pct:.1f}%)"
        )
        sys.stdout.flush()
    print(f"\nDownloading (this may take 10-30+ minutes depending on your connection):")
    urllib.request.urlretrieve(url, dest, reporthook=_hook)
    print("\n  Download complete.")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ---- 1. Find & download ------------------------------------------------
    print("=" * 70)
    print("Step 1a: Locating SEA-AD Hippocampus h5ad on AWS S3")
    print("=" * 70)

    url, fname = find_download_url()
    local_h5ad = os.path.join(OUTPUT_DIR, fname)

    if os.path.exists(local_h5ad):
        print(f"File already exists locally: {local_h5ad}")
        print("  (delete it to force re-download)")
    else:
        download_with_progress(url, local_h5ad)

    # ---- 2. Load in backed (disk-backed) mode to avoid RAM overflow --------
    print("\n" + "=" * 70)
    print("Step 1b: Loading AnnData in backed mode")
    print("=" * 70)
    adata = sc.read_h5ad(local_h5ad, backed="r")
    print(f"  Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # ---- 3. Inspect metadata -----------------------------------------------
    print("\n--- Metadata preview ---")
    print(f"obs columns: {list(adata.obs.columns)}")

    adnc_col = "Overall AD Neuropathological Change"
    if adnc_col in adata.obs.columns:
        print(f"\nADNC distribution:\n{adata.obs[adnc_col].value_counts()}")
    else:
        # Column may have a slightly different name — search for it
        candidates = [c for c in adata.obs.columns if "neuropath" in c.lower() or "adnc" in c.lower()]
        if candidates:
            adnc_col = candidates[0]
            print(f"  ADNC column found as: '{adnc_col}'")
            print(adata.obs[adnc_col].value_counts())
        else:
            print("  WARNING: ADNC column not found — check obs columns above and update adnc_col.")
            adnc_col = None

    if "Subclass" in adata.obs.columns:
        print(f"\nSubclass distribution:\n{adata.obs['Subclass'].value_counts()}")

    # ---- 4. Check TRAPPC9 --------------------------------------------------
    print(f"\n--- {GENE_OF_INTEREST} check ---")
    if GENE_OF_INTEREST in adata.var_names:
        gene_idx = list(adata.var_names).index(GENE_OF_INTEREST)
        # Read just this gene column (backed mode allows slicing)
        col = adata.X[:, gene_idx]
        if sp.issparse(col):
            n_expr = (col > 0).sum()
        else:
            n_expr = (np.array(col) > 0).sum()
        pct = n_expr / adata.n_obs * 100
        print(f"  {GENE_OF_INTEREST} detected in {n_expr:,} / {adata.n_obs:,} cells ({pct:.2f}%)")
    else:
        print(f"  WARNING: {GENE_OF_INTEREST} not in gene names.")
        partial = [g for g in adata.var_names if "TRAPPC9" in g.upper()]
        if partial:
            print(f"  Possible matches: {partial}")
        else:
            print("  No partial matches found. Gene may be absent from this dataset.")

    # ---- 5. Filter cells ---------------------------------------------------
    print("\n" + "=" * 70)
    print("Step 1c: Filtering cells")
    print("=" * 70)

    if adnc_col:
        mask = adata.obs[adnc_col].isin(ADNC_KEEP)
        print(f"  Keeping ADNC in {ADNC_KEEP}")
        print(f"  Cells retained: {mask.sum():,} / {len(mask):,}")
        # Backed AnnData supports boolean indexing
        adata_filt = adata[mask.values]
    else:
        print("  WARNING: skipping ADNC filter (column not found)")
        adata_filt = adata

    # ---- 6. Load subset into memory ----------------------------------------
    print("\n  Loading filtered subset into memory...")
    adata_sub = adata_filt.to_memory()
    print(f"  Subset shape: {adata_sub.n_obs:,} × {adata_sub.n_vars:,}")

    # ---- 7. Set counts matrix to raw UMIs (needed for DE in R) -------------
    if "UMIs" in adata_sub.layers:
        print("  Setting .X to raw UMI counts (from .layers['UMIs'])")
        adata_sub.layers["normalized"] = adata_sub.X.copy()   # save log1p normalized
        adata_sub.X = adata_sub.layers["UMIs"].copy()
    else:
        print("  .layers['UMIs'] not found — .X left as-is (may be log-normalized)")

    # ---- 8. Slim down metadata columns -------------------------------------
    print("\n  Selecting metadata columns...")
    present = {k: v for k, v in KEEP_COLS.items() if k in adata_sub.obs.columns}
    missing = [k for k in KEEP_COLS if k not in adata_sub.obs.columns]
    if missing:
        print(f"  NOTE: these columns were not found and will be skipped: {missing}")

    adata_sub.obs = adata_sub.obs[list(present.keys())].rename(columns=present)

    # Ensure ADNC has tidy name in renamed obs
    adnc_new = present.get(adnc_col, "ADNC") if adnc_col else "ADNC"
    if adnc_new in adata_sub.obs.columns:
        # Add binary grouping column for convenience in R
        adata_sub.obs["Disease_Group"] = adata_sub.obs[adnc_new].map(
            {"Not AD": "Control", "Low": "Low_AD",
             "Intermediate": "Intermediate_AD", "High": "High_AD"}
        ).astype(str).fillna("Unknown")

    # ---- 9. Drop obsm/obsp/varp to reduce file size (UMAP can be kept) ----
    # Keep UMAP if present
    keep_obsm = {k: v for k, v in adata_sub.obsm.items()
                 if "umap" in k.lower() or "X_umap" in k}
    adata_sub.obsm = keep_obsm

    # Remove obsp (connectivity graphs, large and not needed in R)
    adata_sub.obsp = {}

    # Remove varm, uns (not needed)
    adata_sub.varm = {}
    adata_sub.uns  = {}

    # ---- 10. Save ----------------------------------------------------------
    out_path = os.path.join(OUTPUT_DIR, "SEAAD_HIP_subset_for_R.h5ad")
    print(f"\n  Saving subset to: {out_path}")
    adata_sub.write_h5ad(out_path, compression="gzip")
    print(f"  Done! File size: {os.path.getsize(out_path) / 1e9:.2f} GB")

    print("\n" + "=" * 70)
    print("Step 1 complete!")
    print(f"  Output: {out_path}")
    print("  Next: open step2_TRAPPC9_seurat_analysis.R in RStudio")
    print("=" * 70)

    # Print final obs column summary for R user
    print("\nMetadata columns in output h5ad (obs):")
    for col in adata_sub.obs.columns:
        print(f"  {col}: {adata_sub.obs[col].dtype}  "
              f"  unique={adata_sub.obs[col].nunique()}")


if __name__ == "__main__":
    main()
