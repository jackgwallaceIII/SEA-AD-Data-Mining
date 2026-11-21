"""
step2_qc_metrics.py — SEA-AD: per-nucleus QC metrics for one region.

Computes metrics only. NOTHING is filtered and the h5ad is not modified.
Output is one CSV, one row per nucleus, written to results/.

Decisions made here, and why:

1. NO FILTERING IN THIS STEP.
   Thresholds cannot be justified before the distributions are seen, and
   they cannot be applied fairly before cell type is known. Microglia
   carry substantially less RNA than neurons, so a single global UMI
   floor tuned on the whole population would discard real microglia as
   debris — a bias directly against this project's cell type of interest.
   Thresholds therefore belong after the annotation join, chosen per
   cell type.

2. THE MATRIX IS STREAMED IN ROW CHUNKS VIA h5py, NOT LOADED.
   X holds ~6.61e9 nonzeros as int32. data + indices together need about
   53 GB in memory, which exceeds this machine's 32 GB. There is no chunk
   size or swap setting that makes a full load possible.
   CSR storage means rows (nuclei) are contiguous on disk, so reading a
   block of rows is cheap. indptr[i] gives the offset into data/indices
   where row i begins, so a block of rows maps to one contiguous slice of
   each array. This is why the layout discovered in step 1 matters.
   anndata's backed mode could also slice rows, but going through h5py
   directly makes the access pattern explicit and avoids anndata
   constructing index arrays over the whole matrix — which is what
   produced the 49.3 GiB allocation failure in the July pipeline.

3. MITOCHONDRIAL FRACTION IS COMPUTED, BUT MEANS SOMETHING DIFFERENT HERE.
   In whole-cell scRNA-seq, high mitochondrial content marks dying cells.
   This is single-NUCLEUS data: nuclei should contain almost no
   mitochondrial transcript, so an elevated fraction indicates ambient
   or cytoplasmic carryover rather than a stressed cell. Interpret it as
   a contamination signal, not a viability signal.

4. RIBOSOMAL FRACTION IS COMPUTED AS A SECOND CONTAMINATION PROXY.
   Ribosomal protein transcripts are also cytoplasm-enriched, so the same
   logic applies. Two independent proxies are more informative than one,
   and they should broadly agree; nuclei where they disagree are worth a
   look.

5. doublet_score AND umi_count ALREADY EXIST in cell_metadata.csv,
   precomputed by Allen. umi_count is recomputed here anyway, as an
   independent check: if this script's total_counts disagrees with
   Allen's umi_count for the same nucleus, something is wrong with the
   row alignment and it needs finding before anything downstream.

6. THE REGION IS A PARAMETER AND THE OUTPUT IS NAMED FOR IT.
   All ten regions in this release share the same gene set, taxonomy and
   CSR layout, so this computation is region-independent; only the matrix
   differs. The output carries the region in its filename, because a
   fixed name would be silently overwritten by the next region analysed
   and step3 would then join one region's QC metrics to another region's
   annotations. raw counts are required — log2 is pre-normalized, and
   summing it gives a meaningless total.

Runtime: scales with the matrix. MTG is 15.19 GB / 6.61e9 nonzeros and
takes about 5 minutes; the smallest region (LEC, 2.05 GB) is well under
one. Progress prints per chunk.

Run from the project root:
    python scripts/step2_qc_metrics.py
    python scripts/step2_qc_metrics.py --region HIP
"""

import argparse
import time

import anndata as ad
import h5py
import numpy as np
import pandas as pd

import seaad
from seaad import section

# Rows per chunk. At ~5,600 nonzeros per nucleus on average, 20,000 rows
# is roughly 112M nonzeros, or ~900 MB for data + indices together.
# Lower this if memory pressure shows up; raise it for speed.
CHUNK_ROWS = 20_000


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    seaad.add_region_args(ap)
    ap.add_argument("--chunk-rows", type=int, default=CHUNK_ROWS,
                    metavar="N", help=f"rows per read (default: {CHUNK_ROWS:,})")
    args = ap.parse_args()

    section("RELEASE")
    _, manifest = seaad.open_cache()
    region = seaad.resolve_region(manifest, args.region)
    # raw only: see decision 6. QC totals over log2 values are meaningless.
    H5AD = seaad.h5ad_path(manifest, region, "raw")
    OUT = seaad.qc_csv(region)
    chunk_rows = args.chunk_rows
    print(f"  region: {region}")
    print(f"  matrix: {H5AD.name} ({H5AD.stat().st_size / 1e9:.2f} GB)")
    print(f"  output: {OUT}")

    # ---- Gene annotation, and the two gene sets we need ----------------
    # var is small, so read it through anndata in backed mode. X stays on
    # disk; this only pulls obs and var.
    section("GENE SETS")
    adata = ad.read_h5ad(H5AD, backed="r")
    var = adata.var.copy()
    obs_names = adata.obs_names.to_numpy()
    adata.file.close()

    # Use pandas' string accessor rather than np.char. Under pandas 3.x
    # .astype(str).to_numpy() yields an object-dtype array, which numpy's
    # string ufuncs reject; .str.startswith works regardless of backing
    # dtype and returns a clean boolean array.
    sym = var["gene_symbol"].astype("string").fillna("")
    symbols = sym.to_numpy(dtype=object)

    # Mitochondrial genes: the MT- prefix is the GRCh38 convention for the
    # 37 mitochondrially encoded genes. Matching on symbol rather than a
    # curated list keeps this tied to the file's own annotation.
    mito_mask = sym.str.startswith("MT-").to_numpy(dtype=bool)

    # Ribosomal protein genes: RPS* (small subunit) and RPL* (large).
    # MRPS/MRPL are mitochondrial ribosomal proteins and are nuclear-
    # encoded, so they are deliberately NOT matched here — startswith
    # excludes them since they begin with M.
    ribo_mask = (
        sym.str.startswith("RPS") | sym.str.startswith("RPL")
    ).to_numpy(dtype=bool)

    print(f"  genes total:      {len(symbols):,}")
    print(f"  mitochondrial:    {mito_mask.sum()} "
          f"{sorted(map(str, symbols[mito_mask]))[:5]} ...")
    print(f"  ribosomal:        {ribo_mask.sum()}")

    if mito_mask.sum() == 0:
        print("  WARNING: no MT- genes found. Check the gene_symbol column.")

    # ---- Stream the matrix in row chunks -------------------------------
    section("STREAMING X")
    t0 = time.time()

    with h5py.File(H5AD, "r") as f:
        X = f["X"]
        n_obs, n_var = X.attrs["shape"]

        # indptr has n_obs + 1 entries and is small enough to hold whole.
        # indptr[i]:indptr[i+1] is the slice of data/indices belonging to
        # row i, which is what makes row-block reads contiguous.
        indptr = X["indptr"][:]
        data_ds = X["data"]
        indices_ds = X["indices"]

        total_counts = np.zeros(n_obs, dtype=np.float64)
        n_genes = np.zeros(n_obs, dtype=np.int32)
        mito_counts = np.zeros(n_obs, dtype=np.float64)
        ribo_counts = np.zeros(n_obs, dtype=np.float64)

        n_chunks = int(np.ceil(n_obs / chunk_rows))
        for c, r0 in enumerate(range(0, n_obs, chunk_rows)):
            r1 = min(r0 + chunk_rows, n_obs)

            # One contiguous read per array for this block of rows.
            lo, hi = int(indptr[r0]), int(indptr[r1])
            dat = data_ds[lo:hi].astype(np.float64)
            idx = indices_ds[lo:hi]

            # Rebuild which row each nonzero belongs to. row_lengths is
            # the number of detected genes per nucleus; repeating the
            # local row number by that length labels every nonzero.
            row_lengths = np.diff(indptr[r0 : r1 + 1]).astype(np.int64)
            row_ids = np.repeat(np.arange(r1 - r0, dtype=np.int64), row_lengths)
            size = r1 - r0

            # bincount with weights sums the values per row.
            total_counts[r0:r1] = np.bincount(row_ids, weights=dat, minlength=size)
            n_genes[r0:r1] = row_lengths

            m = mito_mask[idx]
            if m.any():
                mito_counts[r0:r1] = np.bincount(
                    row_ids[m], weights=dat[m], minlength=size
                )

            b = ribo_mask[idx]
            if b.any():
                ribo_counts[r0:r1] = np.bincount(
                    row_ids[b], weights=dat[b], minlength=size
                )

            if c % 10 == 0 or c == n_chunks - 1:
                el = time.time() - t0
                print(f"  chunk {c + 1:>3}/{n_chunks}  rows {r0:>9,}-{r1:<9,}  {el:6.1f}s")

    # ---- Assemble -------------------------------------------------------
    section("RESULT")
    qc = pd.DataFrame(
        {
            "cell_label": obs_names,
            "total_counts": total_counts,
            "n_genes_detected": n_genes,
            "mito_counts": mito_counts,
            "ribo_counts": ribo_counts,
        }
    )
    # Guard against divide-by-zero for any all-zero nucleus.
    denom = qc["total_counts"].replace(0, np.nan)
    qc["pct_mito"] = 100.0 * qc["mito_counts"] / denom
    qc["pct_ribo"] = 100.0 * qc["ribo_counts"] / denom

    print(f"  nuclei: {len(qc):,}")
    print(f"  all-zero nuclei: {(qc['total_counts'] == 0).sum():,}")
    print()
    print(
        qc[["total_counts", "n_genes_detected", "pct_mito", "pct_ribo"]]
        .describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99])
        .to_string()
    )

    qc.to_csv(OUT, index=False)
    print(f"\n  written: {OUT}  ({OUT.stat().st_size / 1e6:.1f} MB)")
    print(f"  elapsed: {(time.time() - t0) / 60:.1f} min")
    print("\nNothing was filtered. The h5ad was not modified.")


if __name__ == "__main__":
    main()
