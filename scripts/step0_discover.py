"""
step0_discover.py — SEA-AD: provenance and discovery.

What this is: a record of how the dataset for this project was chosen and
located. It downloads no expression data and performs no analysis. Running
it from a clean checkout reproduces every fact the choice rests on, in
about a minute.

Why these choices were made, in the order they were made:

1. REGION: MTG is the default, and why.
   SEA-AD's multiregion release covers ten regions, and every step takes
   --region so any of them can be analysed. MTG is the default because
   it carries the deepest donor characterization (Luminex Abeta/tau protein
   measurements exist for MTG and no other region) and sits at mid-Braak
   involvement across the cohort, so there is dynamic range in pathology
   burden between donors rather than floor or ceiling effects.
   Cost accepted: no within-donor low-pathology comparator region, so
   disease-stage contrasts are between-donor and donor-level covariates
   must be modeled explicitly. Running a second region with --region is
   one way to recover a within-donor comparator; the donor and disease
   tables are shared across all ten, so they join without extra work.

2. RELEASE: SEA-AD-Multiregion-10X, not the 2024 flagship MTG atlas.
   The flagship MTG atlas is peer-reviewed and has more scrutiny behind
   it, but it is NOT distributed through the ABC Atlas — only the
   multiregion and caudate releases are (see the directory listing this
   script prints). Reaching the flagship object means hand-copying a URL
   from the SEA-AD S3 bucket, outside any manifest.
   The multiregion release was chosen because every file is pinned by a
   versioned manifest with an md5 per file, which makes "which data did I
   use" answerable later. That reproducibility was judged worth more than
   the flagship's additional review, particularly since this project's
   cell type of interest (microglia) sits at subclass level, which is the
   stable part of the taxonomy. Supertype-level calls in this release are
   preprint-stage and should be treated as provisional.

3. NORMALIZATION: raw, not log2.
   Each region ships as both. log2 is pre-normalized by choices that
   cannot be undone, and QC metrics require raw counts to be meaningful.
   log2 can be derived from raw; the reverse is not true.

4. CELL TYPES: all, not microglia only.
   Scope decision made to keep ligand-receptor analysis available later,
   which requires sender populations. Note this does not change the file
   downloaded — MTG-10X-raw.h5ad contains all cell types regardless.

5. ORDER: metadata before matrix.
   The donor and disease tables were inspected before committing to the
   15 GB download, to confirm that donor-level pathology measures join to
   this cohort. They do: both tables are 84 rows keyed on donor_label, and
   CPS_Global ships with the data, so no external join is needed for a
   continuous disease-stage variable.

Environment: conda env `seaad`, python 3.12, abc_atlas_access 1.4.2.
Exact resolved versions are in requirements-lock.txt.

Run from the project root:
    python scripts/step0_discover.py
    python scripts/step0_discover.py --region HIP
"""

import argparse

import pandas as pd

import seaad
from seaad import DIRECTORY, TAXONOMY, section

pd.set_option("display.width", 200)


# Subclass label for myeloid cells in the CCN20260630 taxonomy. Note the
# subclass is named "Immune", not "Microglia", and its supertypes are named
# Micro-PVM_*, which lumps microglia with perivascular macrophages. The
# taxonomy provides no label that separates them; any such separation has
# to be done downstream on markers.
IMMUNE_SUBCLASS = "CS20260630_SCLA_029"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    seaad.add_region_args(ap, normalization=True)
    args = ap.parse_args()

    # ---- Which release are we pinned to --------------------------------
    section("MANIFEST")
    cache, manifest = seaad.open_cache()
    REGION = seaad.resolve_region(manifest, args.region)
    NORMALIZATION = args.normalization
    print(f"  region:        {REGION}")
    print(f"  normalization: {NORMALIZATION}")

    # ---- What SEA-AD data exists in the ABC Atlas ----------------------
    # This is the evidence for decision 2 above: there is no MTG-2024
    # directory here, only caudate and multiregion.
    section("SEA-AD DIRECTORIES AVAILABLE")
    for d in cache.list_directories:
        if "SEA" in d.upper():
            print(f"  {d}")

    # ---- Available matrices, and the one we chose ----------------------
    section("REGIONS AVAILABLE IN " + DIRECTORY)
    em = manifest["file_listing"][DIRECTORY]["expression_matrices"]
    for r in seaad.available_regions(manifest):
        raw = em[r]["raw"]["files"]["h5ad"]["size"] / 1e9
        mark = "  <-- selected" if r == REGION else ""
        print(f"  {r:10s} raw {raw:6.2f} GB{mark}")

    # ---- Exact file, size and checksum ---------------------------------
    # Read straight from the manifest rather than from a helper, because
    # the helpers only report directory-level totals. The manifest is the
    # pinned one, not whichever file happens to sort last on disk.
    section("CHOSEN FILE — provenance")
    entry = seaad.matrix_entry(manifest, REGION, NORMALIZATION)
    print(f"  region:        {REGION}")
    print(f"  normalization: {NORMALIZATION}")
    print(f"  version:       {entry['version']}")
    print(f"  relative path: {entry['relative_path']}")
    print(f"  size:          {entry['size'] / 1e9:.2f} GB")
    print(f"  md5:           {entry['file_hash']}")
    print(f"  url:           {entry['url']}")

    # ---- Metadata inventory --------------------------------------------
    section("METADATA FILES AND SIZES")
    md = manifest["file_listing"][DIRECTORY]["metadata"]
    for name, val in md.items():
        size = val["files"].get("csv", {}).get("size", 0)
        print(f"  {name:28s} {size / 1e6:10.1f} MB")

    # ---- Donor-level tables: the join that had to work -----------------
    section("DONOR AND DISEASE TABLES")
    donor = cache.get_metadata_dataframe(directory=DIRECTORY, file_name="donor")
    disease = cache.get_metadata_dataframe(directory=DIRECTORY, file_name="disease")
    print(f"  donor:   {donor.shape[0]} rows x {donor.shape[1]} cols")
    print(f"  disease: {disease.shape[0]} rows x {disease.shape[1]} cols")

    shared = set(donor["donor_label"]) & set(disease["donor_label"])
    print(f"  donor_label overlap: {len(shared)} of {donor.shape[0]}")

    staging = [
        c for c in disease.columns
        if any(k in c for k in ("Thal", "Braak", "CERAD", "CPS", "neuropathological"))
    ]
    print(f"  staging columns present: {staging}")

    # ---- Taxonomy: where the cell type labels live ---------------------
    section("TAXONOMY — myeloid terms")
    terms = cache.get_metadata_dataframe(
        directory=TAXONOMY, file_name="cluster_annotation_term"
    )
    myeloid = terms[
        (terms["label"] == IMMUNE_SUBCLASS)
        | (terms["parent_term_label"] == IMMUNE_SUBCLASS)
    ]
    print(myeloid[["label", "name", "cluster_annotation_term_set_name"]].to_string())

    section("SUMMARY")
    print(f"  manifest release:  {cache.current_manifest}")
    print(f"  matrix:            {REGION}-{NORMALIZATION}.h5ad "
          f"({entry['size'] / 1e9:.2f} GB, version {entry['version']})")
    print(f"  md5:               {entry['file_hash']}")
    print(f"  donors:            {donor.shape[0]}")
    print("\n  Note: the matrix holds ~6.61e9 nonzeros as int32. A full")
    print("  in-memory load needs roughly 53 GB (data + indices) and is")
    print("  not possible on a 32 GB machine. All access must be backed")
    print("  or row-chunked. See step1_inventory.py.")


if __name__ == "__main__":
    main()
