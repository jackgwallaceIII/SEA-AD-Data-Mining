"""
step3_annotation.py — SEA-AD: per-nucleus annotation table for one region.

Joins every per-nucleus and per-donor annotation onto one table, one row
per nucleus. Touches NO expression data and writes no h5ad.

Output: results/step3_annotation_<REGION>.csv

Decisions made here, and why:

1. EVERY JOIN IS ON cell_label, NEVER cell_barcode.
   cell_barcode has 1,033,197 unique values across 1,178,694 nuclei —
   barcodes repeat across libraries. Joining on it would silently
   duplicate rows and inflate every downstream count. cell_label is
   barcode + sample suffix and is unique. The script asserts row count
   after each merge so a duplication cannot pass unnoticed.

2. cell_metadata.csv IS READ IN CHUNKS WITH usecols.
   The file is 1.64 GB and covers all ten regions. Only five columns are
   needed, and only rows for this region. Reading it whole would cost
   several GB of memory for data that is immediately discarded.
   The region filter is feature_matrix_label == the chosen region;
   those labels match the manifest's region names exactly.

3. COLUMN NAMES ARE RESOLVED, NOT ASSUMED.
   The taxonomy table schemas were not inspected before writing this, so
   each key column is looked up from a list of candidate names and the
   script fails with the actual column list if none match. This is
   cheaper than guessing wrong and re-running blind.

4. ALLEN'S umi_count IS CROSS-CHECKED AGAINST step 2's total_counts.
   These are independent computations of the same quantity. If they
   disagree, row alignment between the h5ad and the metadata CSV is
   broken, and that must be found now rather than inherited by every
   later step. The check is a correlation plus a max absolute difference.

5. THE TABLE IS RESTRICTED TO NUCLEI PRESENT IN THE MATRIX.
   cell_metadata.csv lists more nuclei per region than the released
   matrix contains. For MTG it lists 1,861,334 against the matrix's
   1,178,694; the other 682,640 (36.7%) are nuclei Allen sequenced and
   kept in the metadata but excluded from
   the distributed matrix and from the taxonomy: they carry umi_count and
   doublet_score and nothing else — no expression data, no cluster, no
   cell type.
   Those rows are dropped here, with the count logged, rather than
   carried as NaN. Carrying them would make every downstream row count
   wrong by 58%, and would push NaN handling into every later step for
   nuclei that can never contribute to any result.
   The authoritative list of in-matrix nuclei is step 2's output, whose
   cell_label column was read directly from the h5ad's obs_names. The
   reverse direction is asserted empty: every nucleus in the matrix must
   have a metadata row, or the region filter is wrong.

6. NO QC FILTERING IS APPLIED.
   doublet_score is carried through but not applied. Low-count nuclei are
   carried through. The step 2 finding was that Allen already filtered on
   mitochondrial fraction (<5%) and genes detected (>=500), so the
   remaining QC decisions are doublets and any tighter contamination
   threshold, and they belong in their own step where they can be made
   per cell type.

7. THE REGION IS A PARAMETER AND OUTPUTS ARE NAMED FOR IT.
   cell_metadata.csv and the taxonomy tables cover all ten regions, so
   only the filter value and the output name change between them. Both
   the QC input and the annotation output carry the region, so a second
   region cannot overwrite the first or be joined to it by accident.

Run from the project root:
    python scripts/step3_annotation.py
    python scripts/step3_annotation.py --region HIP
"""

import pathlib
import sys

import argparse

import numpy as np
import pandas as pd

import seaad
from seaad import DIRECTORY, TAXONOMY, section

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 60)

CHUNK = 500_000


def resolve_col(df, candidates, what):
    """Return the first candidate present in df, else fail loudly."""
    for c in candidates:
        if c in df.columns:
            return c
    raise SystemExit(
        f"\nCould not find {what}.\n"
        f"  tried:   {candidates}\n"
        f"  columns: {list(df.columns)}\n"
        f"Tell Claude this column list and the script will be corrected."
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    seaad.add_region_args(ap)
    args = ap.parse_args()

    section("RELEASE")
    cache, manifest = seaad.open_cache()
    REGION = seaad.resolve_region(manifest, args.region)
    QC_CSV = seaad.qc_csv(REGION)
    OUT = seaad.annotation_csv(REGION)
    print(f"  region: {REGION}")
    print(f"  output: {OUT}")

    if not QC_CSV.exists():
        raise SystemExit(
            f"\nMissing {QC_CSV}.\n"
            f"Run step2 for this region first:\n"
            f"  python scripts/step2_qc_metrics.py --region {REGION}"
        )

    # ---- Inspect schemas before joining anything ----------------------
    section("SCHEMAS")
    cm_path = cache.get_metadata_path(directory=DIRECTORY, file_name="cell_metadata")
    cm_head = pd.read_csv(cm_path, nrows=3)
    print(f"cell_metadata ({cm_path.name if hasattr(cm_path,'name') else cm_path}):")
    print(f"  {list(cm_head.columns)}")

    ccm_path = cache.get_metadata_path(
        directory=TAXONOMY, file_name="cell_to_cluster_membership"
    )
    ccm_head = pd.read_csv(ccm_path, nrows=3)
    print("\ncell_to_cluster_membership:")
    print(f"  {list(ccm_head.columns)}")
    print(ccm_head.to_string())

    cam = cache.get_metadata_dataframe(
        directory=TAXONOMY, file_name="cluster_to_cluster_annotation_membership"
    )
    print("\ncluster_to_cluster_annotation_membership:", cam.shape)
    print(f"  {list(cam.columns)}")
    print(cam.head(3).to_string())

    # ---- Resolve the key columns ---------------------------------------
    cm_cell = resolve_col(cm_head, ["cell_label"], "cell key in cell_metadata")
    cm_region = resolve_col(
        cm_head, ["feature_matrix_label"], "region column in cell_metadata"
    )
    ccm_cell = resolve_col(
        ccm_head, ["cell_label"], "cell key in cell_to_cluster_membership"
    )
    ccm_clust = resolve_col(
        ccm_head,
        ["cluster_alias", "cluster_label", "cluster"],
        "cluster key in cell_to_cluster_membership",
    )
    cam_clust = resolve_col(
        cam,
        ["cluster_alias", "cluster_label", "cluster"],
        "cluster key in cluster_to_cluster_annotation_membership",
    )
    cam_term = resolve_col(
        cam,
        ["cluster_annotation_term_name", "name", "cluster_annotation_term_label"],
        "term name in cluster_to_cluster_annotation_membership",
    )
    cam_set = resolve_col(
        cam,
        ["cluster_annotation_term_set_name", "cluster_annotation_term_set_label"],
        "term set in cluster_to_cluster_annotation_membership",
    )
    print(f"\nresolved: cm[{cm_cell}, {cm_region}] "
          f"ccm[{ccm_cell}, {ccm_clust}] cam[{cam_clust}, {cam_set}, {cam_term}]")

    # ---- cell_metadata, chunked, region-filtered -----------------------
    section(f"cell_metadata -> {REGION}")
    want = [cm_cell, cm_region, "donor_label", "doublet_score", "umi_count"]
    want = [c for c in want if c in cm_head.columns]

    parts, n_seen = [], 0
    for chunk in pd.read_csv(cm_path, usecols=want, chunksize=CHUNK):
        n_seen += len(chunk)
        parts.append(chunk[chunk[cm_region] == REGION])
    cm = pd.concat(parts, ignore_index=True)
    del parts

    print(f"  rows scanned (all regions): {n_seen:,}")
    print(f"  rows kept ({REGION}):       {len(cm):,}")
    print(f"  unique {cm_cell}:            {cm[cm_cell].nunique():,}")
    print(f"  donors:                     {cm['donor_label'].nunique()}")
    if cm[cm_cell].nunique() != len(cm):
        raise SystemExit("cell_label is not unique in cell_metadata — stop and inspect.")

    # ---- Restrict to nuclei that exist in the matrix --------------------
    # Done here, before any other join, so the restriction propagates to
    # the cluster-membership read below and nothing downstream ever sees
    # a nucleus without expression data. See decision 5 in the docstring.
    section("RESTRICT TO MATRIX")
    qc = pd.read_csv(QC_CSV)
    in_matrix = set(qc["cell_label"])
    before = len(cm)
    cm = cm[cm[cm_cell].isin(in_matrix)].reset_index(drop=True)
    n_dropped = before - len(cm)

    print(f"  nuclei in cell_metadata ({REGION}): {before:,}")
    print(f"  nuclei in released matrix:          {len(in_matrix):,}")
    print(f"  kept (present in both):             {len(cm):,}")
    print(f"  DROPPED (metadata only, no X):      {n_dropped:,} "
          f"({100.0 * n_dropped / before:.1f}% of metadata rows)")

    # Every nucleus in the matrix must have a metadata row. If any are
    # missing, the region filter or the cell_label convention is wrong
    # and every join below is unsafe.
    n_missing = len(in_matrix) - len(cm)
    if n_missing:
        raise SystemExit(
            f"\n{n_missing:,} nuclei are in the matrix but have no "
            f"cell_metadata row for {REGION}.\n"
            f"The region filter ({cm_region} == '{REGION}') or the "
            f"cell_label convention is wrong. Stop and inspect."
        )
    print(f"  nuclei in matrix but not metadata:  0  (as required)")

    # ---- cluster assignment per nucleus --------------------------------
    # Read chunked and keep only nuclei in our region, so the all-region
    # membership table never lands in memory whole.
    section("cell_to_cluster_membership")
    keep = set(cm[cm_cell])
    parts = []
    for chunk in pd.read_csv(ccm_path, usecols=[ccm_cell, ccm_clust], chunksize=CHUNK):
        parts.append(chunk[chunk[ccm_cell].isin(keep)])
    ccm = pd.concat(parts, ignore_index=True)
    del parts, keep

    print(f"  rows for {REGION}: {len(ccm):,}")
    print(f"  unique clusters:   {ccm[ccm_clust].nunique():,}")
    if ccm[ccm_cell].nunique() != len(ccm):
        raise SystemExit(
            "A nucleus maps to more than one cluster — inspect before joining."
        )

    # ---- cluster -> class / subclass / supertype -----------------------
    # One row per cluster per annotation level; pivot so each cluster has
    # one column per level.
    section("cluster -> annotation levels")
    print(f"  term sets present: {sorted(cam[cam_set].dropna().unique())}")
    levels = cam.pivot_table(
        index=cam_clust, columns=cam_set, values=cam_term, aggfunc="first"
    )
    levels.columns = [str(c).lower().replace(" ", "_") for c in levels.columns]
    levels = levels.reset_index()
    print(f"  clusters annotated: {len(levels):,}")
    print(levels.head(3).to_string())

    # ---- donor-level tables --------------------------------------------
    section("donor + disease")
    donor = cache.get_metadata_dataframe(directory=DIRECTORY, file_name="donor")
    disease = cache.get_metadata_dataframe(directory=DIRECTORY, file_name="disease")
    donor_keep = [
        c for c in ["donor_label", "donor_sex", "donor_age_value", "PMI", "Brain pH"]
        if c in donor.columns
    ]
    disease_keep = [
        c for c in [
            "donor_label", "Overall AD neuropathological Change", "Thal", "Braak",
            "CERAD score", "Cognitive Status", "APOE Genotype", "CPS_Global",
            "Severely Affected Donor", "LATE",
        ] if c in disease.columns
    ]
    print(f"  donor cols kept:   {donor_keep}")
    print(f"  disease cols kept: {disease_keep}")

    # ---- Merge, asserting row count at every step ----------------------
    section("MERGE")
    # n0 is now the matrix's nucleus count, not cell_metadata's. The
    # earlier version baselined on len(cm) before the restriction above,
    # which is why 682,640 annotation-free rows passed every assertion.
    n0 = len(cm)
    assert n0 == len(qc), f"baseline mismatch: cm={n0:,} qc={len(qc):,}"
    print(f"  baseline (nuclei in matrix): {n0:,}")
    ann = cm.drop(columns=[cm_region])

    ann = ann.merge(ccm, left_on=cm_cell, right_on=ccm_cell, how="left")
    if ccm_cell != cm_cell:
        ann = ann.drop(columns=[ccm_cell])
    print(f"  + cluster          {len(ann):>10,}  (expected {n0:,})")
    assert len(ann) == n0, "row count changed joining cluster membership"

    ann = ann.merge(levels, left_on=ccm_clust, right_on=cam_clust, how="left")
    if cam_clust != ccm_clust:
        ann = ann.drop(columns=[cam_clust])
    print(f"  + annotation       {len(ann):>10,}")
    assert len(ann) == n0, "row count changed joining annotation levels"

    ann = ann.merge(donor[donor_keep], on="donor_label", how="left")
    print(f"  + donor            {len(ann):>10,}")
    assert len(ann) == n0, "row count changed joining donor"

    ann = ann.merge(disease[disease_keep], on="donor_label", how="left")
    print(f"  + disease          {len(ann):>10,}")
    assert len(ann) == n0, "row count changed joining disease"

    ann = ann.merge(qc, on="cell_label", how="left")
    print(f"  + QC metrics       {len(ann):>10,}")
    assert len(ann) == n0, "row count changed joining QC metrics"
    # After the restriction this join must be total — no nucleus should
    # be missing its QC metrics.
    n_noqc = int(ann["total_counts"].isna().sum())
    assert n_noqc == 0, f"{n_noqc:,} nuclei have no QC metrics after restriction"

    # ---- Cross-check: Allen's umi_count vs our total_counts -------------
    section("CROSS-CHECK — umi_count vs total_counts")
    if "umi_count" in ann.columns and "total_counts" in ann.columns:
        both = ann[["umi_count", "total_counts"]].dropna()
        if len(both):
            r = np.corrcoef(both["umi_count"], both["total_counts"])[0, 1]
            diff = (both["umi_count"] - both["total_counts"]).abs()
            print(f"  nuclei compared:      {len(both):,}")
            print(f"  pearson r:            {r:.6f}")
            print(f"  max |difference|:     {diff.max():,.0f}")
            print(f"  nuclei differing >1:  {(diff > 1).sum():,}")
            if r < 0.999:
                print("  WARNING: these should be near-identical. Row alignment"
                      " between the h5ad and cell_metadata is suspect.")
        else:
            print("  no overlapping rows — join failed, inspect cell_label values")
    else:
        print("  one of the columns is missing; skipped")

    # ---- Completeness ----------------------------------------------------
    section("COMPLETENESS")
    na = ann.isna().sum()
    na = na[na > 0]
    print(na.to_string() if len(na) else "  no missing values")

    # ---- What is actually in here ----------------------------------------
    section("COMPOSITION")
    for col in ann.columns:
        if col.lower() in ("subclass", "class", "supertype"):
            vc = ann[col].value_counts(dropna=False)
            print(f"\n  {col} ({len(vc)} levels):")
            print(vc.head(35).to_string())

    if "CPS_Global" in ann.columns:
        print("\n  CPS_Global across donors:")
        per_donor = ann.groupby("donor_label")["CPS_Global"].first()
        print(per_donor.describe().to_string())

    section("WRITE")
    ann.to_csv(OUT, index=False)
    print(f"  {OUT}  ({OUT.stat().st_size / 1e6:.1f} MB)")
    print(f"  {len(ann):,} nuclei x {ann.shape[1]} columns")
    print(f"\n  {n_dropped:,} metadata-only nuclei were dropped (no expression")
    print("  data, no cell type). No QC filtering was applied to the rest.")
    print("  No expression data was read.")


if __name__ == "__main__":
    sys.exit(main())
