"""
step1_inventory.py — SEA-AD: load and inventory one region's matrix.

Purpose: confirm what is actually in the downloaded matrix before any
filtering, QC, or analysis. Nothing is modified and nothing is written.

Decisions made here, and why:

1. The file is opened in BACKED mode (backed='r'). This reads obs and var
   into memory but leaves the expression matrix X on disk. The matrix is
   ~15 GB; loading it would consume half the machine's RAM for no reason,
   since this step never touches expression values. Backed mode is also
   what every later step should use until the data has been subset down.

2. h5py is used first, before anndata, to read the raw file structure.
   anndata's reader makes assumptions about layout; h5py reports what is
   literally in the file. If the two disagree, that is worth knowing now
   rather than three steps later.

3. Nothing is filtered. Not doublets, not low-UMI nuclei, not any cell
   type. Filtering decisions come after seeing what the distributions
   look like, not before.

4. The region is a parameter. The ten regions in this release are
   structurally identical — same genes, same taxonomy, same CSR layout —
   so this step works on any of them, and the path is resolved from the
   pinned manifest rather than typed. MTG is the default because it is
   what this project analyses; see decision 1 in step0_discover.py.

Run from the project root:
    python scripts/step1_inventory.py
    python scripts/step1_inventory.py --region HIP
    python scripts/step1_inventory.py --region DFC-10X --normalization log2
"""

import argparse

import anndata as ad
import h5py
import pandas as pd

import seaad
from seaad import section

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 50)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    seaad.add_region_args(ap, normalization=True)
    args = ap.parse_args()

    section("RELEASE")
    _, manifest = seaad.open_cache()
    region = seaad.resolve_region(manifest, args.region)
    H5AD = seaad.h5ad_path(manifest, region, args.normalization)
    print(f"  region:        {region}")
    print(f"  normalization: {args.normalization}")

    section("FILE")
    print(f"path: {H5AD}")
    print(f"size: {H5AD.stat().st_size / 1e9:.2f} GB")

    # ---- Raw HDF5 structure -------------------------------------------
    # Reading the container directly. This shows how X is stored (dense
    # array vs sparse CSR/CSC group), which determines how expensive row
    # subsetting will be later. CSR is cheap to slice by row; CSC is not.
    section("HDF5 STRUCTURE (top level)")
    with h5py.File(H5AD, "r") as f:
        for key in f.keys():
            obj = f[key]
            kind = "group" if isinstance(obj, h5py.Group) else "dataset"
            print(f"  {key:16s} {kind}")

        if "X" in f:
            x = f["X"]
            if isinstance(x, h5py.Group):
                enc = x.attrs.get("encoding-type", "unknown")
                shape = x.attrs.get("shape", "unknown")
                print(f"\n  X encoding: {enc}")
                print(f"  X shape:    {shape}")
                print(f"  X members:  {list(x.keys())}")
                if "data" in x:
                    print(f"  X data dtype: {x['data'].dtype}")
                    print(f"  X nnz:        {x['data'].shape[0]:,}")
            else:
                print(f"\n  X is dense, shape {x.shape}, dtype {x.dtype}")

        if "layers" in f:
            print(f"\n  layers: {list(f['layers'].keys())}")
        if "obsm" in f:
            print(f"  obsm:   {list(f['obsm'].keys())}")

    # ---- AnnData view --------------------------------------------------
    section("ANNDATA (backed, X left on disk)")
    adata = ad.read_h5ad(H5AD, backed="r")
    print(f"shape: {adata.shape[0]:,} nuclei x {adata.shape[1]:,} genes")
    print(f"X in memory: {adata.isbacked is False}")

    section("obs COLUMNS")
    print(list(adata.obs.columns))

    section("obs HEAD")
    print(adata.obs.head(5).to_string())

    section("obs DTYPES")
    print(adata.obs.dtypes.to_string())

    # Categorical columns are the annotation-like ones. Printing their
    # level counts shows what grouping variables exist without assuming
    # any particular column name is present.
    section("CATEGORICAL COLUMNS — level counts")
    for col in adata.obs.columns:
        s = adata.obs[col]
        if str(s.dtype) in ("category", "object", "bool"):
            n = s.nunique(dropna=True)
            print(f"  {col:40s} {n:>8,} levels")
            if n <= 15:
                print(f"      {sorted(map(str, s.dropna().unique()))}")

    section("var COLUMNS")
    print(list(adata.var.columns))
    print(adata.var.head(5).to_string())

    section("MISSING VALUES in obs")
    na = adata.obs.isna().sum()
    na = na[na > 0]
    print(na.to_string() if len(na) else "  none")

    adata.file.close()
    print("\nDone. Nothing was modified.")


if __name__ == "__main__":
    main()
