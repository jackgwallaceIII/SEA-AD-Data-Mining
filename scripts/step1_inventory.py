"""
step1_inventory.py — SEA-AD MTG, load and inventory only.

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

Run from the project root:
    python scripts/step1_inventory.py
"""

import pathlib
import h5py
import anndata as ad
import pandas as pd

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 50)

H5AD = pathlib.Path(
    "abc_cache/expression_matrices/SEA-AD-Multiregion-10X/20260630/MTG-10X-raw.h5ad"
)


def section(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


def main():
    if not H5AD.exists():
        raise SystemExit(f"Not found: {H5AD}\nRun from the project root.")

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
