#!/usr/bin/env python3
"""
=============================================================================
SEA-AD TRAPPC9 Analysis Pipeline — Step 1: Download & Preprocess
=============================================================================
Dataset : Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)
Regions : Any combination of the 11 SEA-AD brain regions (see REGIONS below)
Source  : s3://sea-ad-single-cell-profiling/ (public, no AWS account needed)

What this script does, for each region in REGIONS_TO_PROCESS:
  1. Downloads the SEA-AD snRNA-seq h5ad for that region from AWS S3
  2. Removes "Reference" donors (5 young donors used only for cell type mapping)
  3. Subsets the AnnData to key metadata columns
  4. Saves a compact h5ad that R can load without the 2^32 integer overflow issue

IMPORTANT: The full datasets have >2^32 non-zero matrix elements, which exceeds
R's 32-bit integer limit. This step produces smaller files that load fine in R.
(See: https://community.brain-map.org/t/sea-ad-snrnaseq-dataset-format/3454)

Requirements (install with pip):
    pip install scanpy anndata numpy pandas scipy

Run:
    python step1_download_preprocess.py

Output:
    seaad_data/<REGION_CODE>/SEAAD_<REGION_CODE>_subset_for_R.h5ad
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

OUTPUT_DIR = "seaad_data"

# Each entry describes how to locate that region's h5ad on S3.
#   folder              : top-level S3 prefix for the region
#   subdir               : sub-prefix under the folder (most regions use "RNAseq")
#   filename_template    : filename pattern, with {code} and {date} placeholders
#   candidate_dates       : dates to try, newest first
#
# Most regions follow the standard SEA-AD naming convention:
#   <folder>/RNAseq/SEAAD_<code>_RNAseq_final-nuclei.<date>.h5ad
#
# Two regions are irregular on the S3 bucket as of this writing:
#   - PFC:  the bucket folder is "PFC" but the files inside still use the
#           legacy region code "DFC" (Dorsolateral Frontal Cortex).
#   - Caudate_Nucleus: has no "RNAseq" subdir and uses a different filename
#           pattern entirely ("CaH_all_nuclei.<date>.h5ad").
STANDARD_DATES = ["2026-06-22", "2024-02-13", "2023-05-05", "2022-08-18"]

REGIONS = {
    "HIP": {
        "name": "Hippocampus",
        "folder": "HIP", "code": "HIP", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "MEC": {
        "name": "Medial entorhinal cortex",
        "folder": "MEC", "code": "MEC", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "LEC": {
        "name": "Lateral entorhinal cortex",
        "folder": "LEC", "code": "LEC", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "ITG": {
        "name": "Inferior temporal gyrus",
        "folder": "ITG", "code": "ITG", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "MTG": {
        "name": "Middle temporal gyrus",
        "folder": "MTG", "code": "MTG", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES + ["2022-06-07"],
    },
    "STG": {
        "name": "Superior temporal gyrus",
        "folder": "STG", "code": "STG", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "FI": {
        "name": "Frontal insula",
        "folder": "FI", "code": "FI", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "PFC": {
        "name": "Prefrontal cortex",
        "folder": "PFC", "code": "DFC", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "AnG": {
        "name": "Angular gyrus",
        "folder": "AnG", "code": "AnG", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "V1C": {
        "name": "Primary visual cortex",
        "folder": "V1C", "code": "V1C", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "Caudate_Nucleus": {
        "name": "Caudate Nucleus (head of the caudate)",
        "folder": "Caudate_Nucleus", "code": "CaH", "subdir": "",
        "filename_template": "CaH_all_nuclei.{date}.h5ad",
        "candidate_dates": ["2025-12-29"],
    },
}

# Which regions to actually download & process. Defaults to all of them;
# trim this list to run a subset (e.g. ["HIP", "PFC"]).
REGIONS_TO_PROCESS = list(REGIONS.keys())

S3_BUCKET = "https://sea-ad-single-cell-profiling.s3.amazonaws.com"

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


def find_download_url(region_key: str, cfg: dict) -> tuple:
    """Try each candidate date for a region and return the first working URL."""
    base_url = "/".join(p for p in (S3_BUCKET, cfg["folder"], cfg["subdir"]) if p)
    for date in cfg["candidate_dates"]:
        fname = cfg["filename_template"].format(code=cfg["code"], date=date)
        url = f"{base_url}/{fname}"
        print(f"  Checking {url} ...", end=" ", flush=True)
        if url_exists(url):
            print("FOUND")
            return url, fname
        print("not found")
    raise FileNotFoundError(
        f"\n\nCould not find a valid {cfg['name']} ({region_key}) h5ad on S3.\n"
        "Browse the bucket index to find the correct filename:\n"
        f"  {S3_BUCKET}/index.html#{cfg['folder']}/\n"
        "Then add that date string to the region's candidate_dates in REGIONS."
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
# PER-REGION PIPELINE
# ---------------------------------------------------------------------------

def process_region(region_key: str) -> None:
    cfg = REGIONS[region_key]
    region_dir = os.path.join(OUTPUT_DIR, region_key)
    os.makedirs(region_dir, exist_ok=True)

    # ---- 1. Find & download ------------------------------------------------
    print("=" * 70)
    print(f"Step 1a [{region_key}]: Locating {cfg['name']} h5ad on AWS S3")
    print("=" * 70)

    url, fname = find_download_url(region_key, cfg)
    local_h5ad = os.path.join(region_dir, fname)

    if os.path.exists(local_h5ad):
        print(f"File already exists locally: {local_h5ad}")
        print("  (delete it to force re-download)")
    else:
        download_with_progress(url, local_h5ad)

    # ---- 2. Load in backed (disk-backed) mode to avoid RAM overflow --------
    print("\n" + "=" * 70)
    print(f"Step 1b [{region_key}]: Loading AnnData in backed mode")
    print("=" * 70)
    adata = sc.read_h5ad(local_h5ad, backed="r")
    print(f"  Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

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
    print(f"Step 1c [{region_key}]: Filtering cells")
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
    print(f"  Subset shape: {adata_sub.n_obs:,} x {adata_sub.n_vars:,}")

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
    adata_sub.obs["Region"] = region_key

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
    adata_sub.uns = {}

    # ---- 10. Save ----------------------------------------------------------
    out_path = os.path.join(region_dir, f"SEAAD_{region_key}_subset_for_R.h5ad")
    print(f"\n  Saving subset to: {out_path}")
    adata_sub.write_h5ad(out_path, compression="gzip")
    print(f"  Done! File size: {os.path.getsize(out_path) / 1e9:.2f} GB")

    print("\n" + "=" * 70)
    print(f"Step 1 complete for {region_key}!")
    print(f"  Output: {out_path}")
    print("=" * 70)

    # Print final obs column summary for R user
    print(f"\nMetadata columns in output h5ad (obs) for {region_key}:")
    for col in adata_sub.obs.columns:
        print(f"  {col}: {adata_sub.obs[col].dtype}  "
              f"  unique={adata_sub.obs[col].nunique()}")


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    unknown = [r for r in REGIONS_TO_PROCESS if r not in REGIONS]
    if unknown:
        raise ValueError(f"Unknown region(s) in REGIONS_TO_PROCESS: {unknown}. "
                          f"Valid options: {list(REGIONS.keys())}")

    for region_key in REGIONS_TO_PROCESS:
        process_region(region_key)

    print("\n" + "=" * 70)
    print(f"All regions processed: {REGIONS_TO_PROCESS}")
    print("  Next: open step2_TRAPPC9_seurat_analysis.R in RStudio")
    print("=" * 70)


if __name__ == "__main__":
    main()
