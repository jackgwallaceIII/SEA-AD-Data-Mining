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

MEMORY NOTE (added after MTG OOM crash):
Large regions (MTG: ~1.18M cells, PFC likely similar) retain close to 100% of
cells after ADNC filtering. Row-slicing a backed CSR matrix by boolean mask
the naive way (`X_backed[mask]`) forces anndata to concatenate index/data
slices for EVERY kept row into ONE array in RAM -- for MTG this needed 49+ GB
just for the indices array. Instead, this version streams the filtered CSR
matrix to the output file in row-chunks (default 50,000 rows), so peak memory
use is bounded by the chunk size, not the dataset size, no matter how big a
future region is.

Requirements (install with pip):
    pip install anndata numpy pandas scipy h5py

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

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

# read_elem / sparse_dataset / write_elem moved between anndata.experimental
# and anndata.io across versions -- try the modern locations first.
try:
    from anndata import read_elem
except ImportError:
    from anndata.experimental import read_elem
try:
    from anndata import sparse_dataset
except ImportError:
    try:
        from anndata.io import sparse_dataset
    except ImportError:
        from anndata.experimental import sparse_dataset
try:
    from anndata.io import write_elem
except ImportError:
    try:
        from anndata import write_elem
    except ImportError:
        from anndata.experimental import write_elem

# ---------------------------------------------------------------------------
# CONFIGURATION — adjust if needed
# ---------------------------------------------------------------------------

OUTPUT_DIR = "seaad_data"

# How many rows to process per chunk when streaming large sparse matrices.
# Lower this (e.g. 20_000) if you still see memory errors; raise it
# (e.g. 100_000) for speed if you have plenty of RAM.
ROW_CHUNK_SIZE = 50_000

STANDARD_DATES = ["2026-06-22", "2024-02-13", "2023-05-05", "2022-08-18"]

REGIONS = {
    "HIP": {
        "name": "Hippocampus",
        "folder": "HIP", "code": "HIP", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
    "MTG": {
        "name": "Middle temporal gyrus",
        "folder": "MTG", "code": "MTG", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES + ["2022-06-07"],
    },
    "PFC": {
        "name": "Prefrontal cortex",
        "folder": "PFC", "code": "DFC", "subdir": "RNAseq",
        "filename_template": "SEAAD_{code}_RNAseq_final-nuclei.{date}.h5ad",
        "candidate_dates": STANDARD_DATES,
    },
}

REGIONS_TO_PROCESS = list(REGIONS.keys())

S3_BUCKET = "https://sea-ad-single-cell-profiling.s3.amazonaws.com"

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

ADNC_KEEP = ["Not AD", "Low", "Intermediate", "High"]
GENE_OF_INTEREST = "TRAPPC9"


# ---------------------------------------------------------------------------
# HELPER FUNCTIONS — download / discovery (unchanged)
# ---------------------------------------------------------------------------

def url_exists(url: str, timeout: int = 15) -> bool:
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout):
            return True
    except (urllib.error.URLError, urllib.error.HTTPError):
        return False


def find_download_url(region_key: str, cfg: dict) -> tuple:
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
# NEW: memory-safe streaming CSR filter/copy
# ---------------------------------------------------------------------------

def stream_filtered_csr_to_group(h5in_group, mask_arr, h5out_parent, out_name,
                                  row_chunk=ROW_CHUNK_SIZE,
                                  indices_dtype=np.int32,
                                  gene_idx=None):
    """
    Read a disk-backed CSR matrix (h5in_group: h5py Group with data/indices/
    indptr datasets) in row-chunks, keep only rows where mask_arr is True,
    and write the result straight to h5out_parent[out_name] as a new
    anndata-compatible CSR group -- WITHOUT ever holding the full matrix
    in RAM. Peak memory is bounded by `row_chunk` rows at a time.

    If gene_idx is given, also returns the number of kept rows where that
    gene column has a nonzero value (for the TRAPPC9 detection count).
    """
    indptr_in = h5in_group["indptr"][:]          # small (n_obs+1 int64s), safe to load whole
    data_in = h5in_group["data"]
    indices_in = h5in_group["indices"]
    n_obs = indptr_in.shape[0] - 1
    n_vars = int(h5in_group.attrs["shape"][1])
    data_dtype = data_in.dtype

    n_kept = int(mask_arr.sum())
    out_group = h5out_parent.create_group(out_name)
    out_data = out_group.create_dataset("data", shape=(0,), maxshape=(None,),
                                         dtype=data_dtype, chunks=True)
    out_indices = out_group.create_dataset("indices", shape=(0,), maxshape=(None,),
                                            dtype=indices_dtype, chunks=True)
    out_indptr = np.zeros(n_kept + 1, dtype=np.int64)

    write_pos = 0
    kept_rows_so_far = 0
    gene_detected_count = 0

    for start in range(0, n_obs, row_chunk):
        end = min(start + row_chunk, n_obs)
        block_mask = mask_arr[start:end]
        if not block_mask.any():
            continue

        elem_start, elem_end = int(indptr_in[start]), int(indptr_in[end])
        block_data = data_in[elem_start:elem_end]
        block_indices = indices_in[elem_start:elem_end].astype(indices_dtype, copy=False)
        local_indptr = indptr_in[start:end + 1] - elem_start

        keep_rows = np.where(block_mask)[0]
        slices = [slice(int(local_indptr[i]), int(local_indptr[i + 1])) for i in keep_rows]

        if slices:
            block_data_kept = np.concatenate([block_data[s] for s in slices])
            block_indices_kept = np.concatenate([block_indices[s] for s in slices])
        else:
            block_data_kept = np.empty(0, dtype=data_dtype)
            block_indices_kept = np.empty(0, dtype=indices_dtype)

        if gene_idx is not None:
            for s in slices:
                if gene_idx in block_indices[s]:
                    gene_detected_count += 1

        n_new = block_data_kept.shape[0]
        new_size = write_pos + n_new
        out_data.resize((new_size,))
        out_indices.resize((new_size,))
        out_data[write_pos:new_size] = block_data_kept
        out_indices[write_pos:new_size] = block_indices_kept

        row_lengths = np.diff(local_indptr)[keep_rows]
        cum = write_pos + np.cumsum(row_lengths)
        out_indptr[kept_rows_so_far + 1: kept_rows_so_far + 1 + len(keep_rows)] = cum

        write_pos = new_size
        kept_rows_so_far += len(keep_rows)
        print(f"    ...streamed {end:,}/{n_obs:,} rows "
              f"({kept_rows_so_far:,} kept, {write_pos:,} nonzeros written)", end="\r")

    print()
    out_group.create_dataset("indptr", data=out_indptr)
    out_group.attrs["encoding-type"] = "csr_matrix"
    out_group.attrs["encoding-version"] = "0.1.0"
    out_group.attrs["shape"] = np.array([n_kept, n_vars])

    return n_kept, gene_detected_count


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

    # ---- 2. Open the h5ad without eagerly loading large groups -------------
    print("\n" + "=" * 70)
    print(f"Step 1b [{region_key}]: Opening h5ad (lazy, X/UMIs not yet loaded)")
    print("=" * 70)
    h5file = h5py.File(local_h5ad, "r")
    obs = read_elem(h5file["obs"])
    var = read_elem(h5file["var"])
    n_obs, n_vars = h5file["X"].attrs["shape"] if "shape" in h5file["X"].attrs else sparse_dataset(h5file["X"]).shape
    print(f"  Shape: {n_obs:,} cells x {n_vars:,} genes")

    # ---- 3. Inspect metadata -----------------------------------------------
    print("\n--- Metadata preview ---")
    print(f"obs columns: {list(obs.columns)}")

    adnc_col = "Overall AD Neuropathological Change"
    if adnc_col in obs.columns:
        print(f"\nADNC distribution:\n{obs[adnc_col].value_counts()}")
    else:
        candidates = [c for c in obs.columns if "neuropath" in c.lower() or "adnc" in c.lower()]
        if candidates:
            adnc_col = candidates[0]
            print(f"  ADNC column found as: '{adnc_col}'")
            print(obs[adnc_col].value_counts())
        else:
            print("  WARNING: ADNC column not found — check obs columns above and update adnc_col.")
            adnc_col = None

    if "Subclass" in obs.columns:
        print(f"\nSubclass distribution:\n{obs['Subclass'].value_counts()}")

    gene_idx = None
    print(f"\n--- {GENE_OF_INTEREST} check ---")
    if GENE_OF_INTEREST in var.index:
        gene_idx = var.index.get_loc(GENE_OF_INTEREST)
        print(f"  {GENE_OF_INTEREST} found in var_names at index {gene_idx} "
              f"(expression count will be reported after subsetting).")
    else:
        print(f"  WARNING: {GENE_OF_INTEREST} not in gene names.")
        partial = [g for g in var.index if "TRAPPC9" in str(g).upper()]
        if partial:
            print(f"  Possible matches: {partial}")
        else:
            print("  No partial matches found. Gene may be absent from this dataset.")

    # ---- 4. Filter cells (mask only, no data read yet) ---------------------
    print("\n" + "=" * 70)
    print(f"Step 1c [{region_key}]: Filtering cells")
    print("=" * 70)

    if adnc_col:
        mask = obs[adnc_col].isin(ADNC_KEEP)
        print(f"  Keeping ADNC in {ADNC_KEEP}")
        print(f"  Cells retained after ADNC filter: {mask.sum():,} / {len(mask):,}")
    else:
        print("  WARNING: skipping ADNC filter (column not found)")
        mask = pd.Series(True, index=obs.index)

    if SUBCLASS_KEEP is not None and "Subclass" in obs.columns:
        subclass_mask = obs["Subclass"].isin(SUBCLASS_KEEP)
        mask = mask & subclass_mask
        print(f"  Keeping Subclass in {SUBCLASS_KEEP}")
        print(f"  Cells retained after Subclass filter: {mask.sum():,} / {len(mask):,}")
    elif SUBCLASS_KEEP is not None:
        print("  WARNING: 'Subclass' column not found — skipping subclass filter")

    mask_arr = mask.to_numpy()
    obs_sub = obs[mask_arr].copy()

    # ---- 5. Slim metadata columns down BEFORE any big I/O ------------------
    print("\n  Selecting metadata columns...")
    present = {k: v for k, v in KEEP_COLS.items() if k in obs_sub.columns}
    missing = [k for k in KEEP_COLS if k not in obs_sub.columns]
    if missing:
        print(f"  NOTE: these columns were not found and will be skipped: {missing}")

    obs_sub = obs_sub[list(present.keys())].rename(columns=present)
    obs_sub["Region"] = region_key

    adnc_new = present.get(adnc_col, "ADNC") if adnc_col else "ADNC"
    if adnc_new in obs_sub.columns:
        obs_sub["Disease_Group"] = obs_sub[adnc_new].map(
            {"Not AD": "Control", "Low": "Low_AD",
             "Intermediate": "Intermediate_AD", "High": "High_AD"}
        ).astype(str).fillna("Unknown")

    # ---- 6. Stream-write the filtered matrices straight to the output file
    # This is the memory-safe replacement for the old `X_backed[mask_arr]`
    # line that OOM'd on MTG. Nothing here holds more than ROW_CHUNK_SIZE
    # rows' worth of data in RAM at once.
    out_path = os.path.join(region_dir, f"SEAAD_{region_key}_subset_for_R.h5ad")
    print("\n" + "=" * 70)
    print(f"Step 1d [{region_key}]: Streaming filtered matrix to {out_path}")
    print("  (processing in chunks of {:,} rows -- this avoids the OOM)".format(ROW_CHUNK_SIZE))
    print("=" * 70)

    with h5py.File(out_path, "w") as out_file:
        # Original X becomes the "normalized" layer in the output.
        print("  Writing layers/UMIs -> normalized ...")
        _, gene_detected_normalized = stream_filtered_csr_to_group(
            h5file["X"], mask_arr, out_file, "layers_tmp_normalized",
            gene_idx=gene_idx,
        )

        # Final .X is the raw UMI counts, if present; otherwise reuse the
        # normalized data we just streamed (rename below).
        if "layers" in h5file and "UMIs" in h5file["layers"]:
            print("  Writing raw UMI counts -> X ...")
            n_kept, gene_detected = stream_filtered_csr_to_group(
                h5file["layers"]["UMIs"], mask_arr, out_file, "X",
                gene_idx=gene_idx,
            )
            has_umis = True
        else:
            print("  .layers['UMIs'] not found in source file — X will be the "
                  "same data as the normalized layer.")
            n_kept, gene_detected = mask_arr.sum(), gene_detected_normalized
            has_umis = False

        # obsm (e.g. UMAP) is small and dense -- safe to read fully and mask
        # in memory before we close the source file.
        obsm_to_write = {}
        if "obsm" in h5file:
            for key in h5file["obsm"].keys():
                if "umap" in key.lower():
                    arr = read_elem(h5file["obsm"][key])
                    obsm_to_write[key] = np.asarray(arr)[mask_arr]

        h5file.close()

        # Rename/organize the layers group properly now that we know
        # whether UMIs existed.
        layers_group = out_file.create_group("layers")
        if has_umis:
            out_file.move("layers_tmp_normalized", "layers/normalized")
        else:
            out_file.move("layers_tmp_normalized", "X")

        if gene_idx is not None:
            pct = gene_detected / n_kept * 100 if n_kept else 0
            print(f"  {GENE_OF_INTEREST} detected in {gene_detected:,} / {n_kept:,} "
                  f"retained cells ({pct:.2f}%)")

        # ---- 7. Write obs/var (small — fine to write directly) -------------
        print("\n  Writing obs/var metadata...")
        write_elem(out_file, "obs", obs_sub)
        write_elem(out_file, "var", var)

        # ---- 8. obsm: keep UMAP only, if present (small, dense) ------------
        if obsm_to_write:
            obsm_group = out_file.create_group("obsm")
            for key, arr in obsm_to_write.items():
                write_elem(obsm_group, key, arr)

        out_file.attrs["encoding-type"] = "anndata"
        out_file.attrs["encoding-version"] = "0.1.0"

    print(f"\n  Done! File size: {os.path.getsize(out_path) / 1e9:.2f} GB")
    print("\n" + "=" * 70)
    print(f"Step 1 complete for {region_key}!")
    print(f"  Output: {out_path}")
    print("=" * 70)


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
