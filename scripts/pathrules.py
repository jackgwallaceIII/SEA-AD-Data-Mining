"""
pathrules.py — classify an S3 object key: region, modality, product, format.

Not a step. Pure functions over (bucket, key) strings: no network, no
file contents. The catalog (step 4) runs every listed key through
classify() and groups the results into rows.

    python scripts/pathrules.py sea-ad-single-cell-profiling \\
        MTG/RNAseq/SEAAD_MTG_RNAseq_final-nuclei.2026-06-22.h5ad

Every rule below was written against keys seen in read-only listings of
the three open buckets on 2026-11-21 (see tests/test_pathrules.py, which
pins those keys). Rules are ordered; the first match wins.

Decisions made here, and why:

1. A "UNIT" IS NOT ALWAYS ONE S3 OBJECT.
   S3 has no folders, only keys. Some datasets are thousands of keys
   that only mean something together:
     - a Zarr store (combined_spatialdata_zarr/): one dataset, ~3,000
       chunk keys per tissue section;
     - a Xenium instrument output folder (output-XETG...): ~200 files;
     - deep-zoom image tiles (<slide>_files/<level>/<x>_<y>.jpeg):
       ~800 JPEG tiles per slide, which exist only for a web viewer;
     - MTG/RNAseq/pseudoprogression-plots/: 1,134,632 plot images
       (1.38 TB), one folder per gene;
     - Hive-partitioned Parquet (counts_by_gene/gene=A1BG/...): one
       table split into ~36,000 parts.
   Each of these collapses to a single unit (its folder). Without this
   the single-cell bucket alone would be ~1.9 million catalog rows, of
   which ~1,500 are the files anyone actually downloads.

2. REGION COMES FROM THE FIRST PATH FOLDER THAT IS A REGION, THEN FROM
   THE FILE NAME, OTHERWISE NOTHING.
   Folders are checked with regions.codes_for() — exact, never fuzzy.
   Project folders (DREAM/, Multiregion_2026/) are not regions, so the
   file name's underscore-separated tokens are tried next:
   SEAAD_MTG_RNAseq_DREAM -> MTG; ..._HIP_MEC_LEC_ITG_MTG_FI_STG_DFC_
   AnG_V1C_... -> all ten. A folder region always wins over the name.
   Multiregion_2026 model outputs name no region at all; they are
   cross-region by design, so they get region_from = 'project:
   multi-region release' with no codes, rather than being copied into
   all ten regions or flagged as missing. DREAM/SEAAD_A9_... resolves to
   DFC because 'A9' (Brodmann area 9) is a DFC alias in regions.csv —
   added by a human decision on the bucket README's evidence, not by a
   rule here. Name tokens only ever match vocabulary entries exactly.
   CAVEAT for interpretation: "filed under MTG" means "in the MTG folder".
   MTG/RNAseq/Supplementary Information/scCODA Results/PFC_Mathys_2023/
   is filed under MTG but is about a PFC dataset. region_from records
   where the answer came from so this can be checked.

3. MODALITY IS ONLY ASSERTED WHEN THE PATH SAYS SO.
   RNAseq/ATACseq folders or name tokens -> snRNA-seq / snATAC-seq.
   Whole buckets with one modality (spatial, neuropathology) use a
   bucket default. One project default exists: Multiregion_2026 ->
   snRNA-seq, because every data file there that names a modality says
   RNAseq and its models are scVI/scANVI (RNA) models; modality_from
   says 'project default' so it is never mistaken for path evidence.
   Caudate_Nucleus/CaH_all_nuclei.2025-12-29.h5ad names no modality
   and is left unclassified rather than assumed to be RNA.
   The spatial bucket is NOT split into MERFISH vs Xenium: only paths
   containing 'xenium' prove the platform, and the MTG sections
   (DAPI_Max.tif, cellpose-detected_transcripts.csv) do not say.

4. DEPRECATED IS A FLAG, NOT A DELETION.
   Anything under previous_objects/ is marked deprecated=True. It stays
   classifiable so the catalog can count it and say what it is; hiding
   it is the job of the show/fetch commands (default: hidden).

5. SKIPPED MEANS "NOT DATA", AND THE REASON IS KEPT.
   Folder-marker keys ending in '/', Jupyter checkpoint copies, and the
   neuropathology bucket's web-viewer (simple_viewer/) and scratch
   (temp/) folders are skipped with a reason, so the catalog can report
   how much was skipped instead of silently dropping it.

6. UNCLASSIFIED IS A RESULT, NOT AN ERROR.
   Anything the rules cannot place gets 'unclassified' in the field it
   is missing, and is listed in `missing`. These are for a human to
   review; a guessed label would be worse than no label.

7. TEST DATA GETS NO DEFAULTS.
   Any folder named TEST_* (seen: Multiregion_2026/previous_objects/
   parquet/TEST_dev_mouse/, a mouse development table) is labelled
   product 'test/development data' and does NOT inherit the project's
   defaults. Without this rule it was filed as human multi-region
   snRNA-seq, which is wrong on two counts; a default is only safe for
   files that really belong to the project's dataset.
"""

import dataclasses
import re
import sys

import regions as regions_mod

# ---- Per-bucket configuration ------------------------------------------

BUCKETS = {
    "sea-ad-single-cell-profiling": {
        "default_modality": None,
        "skip_top": {},
    },
    "sea-ad-spatial-transcriptomics": {
        "default_modality": "spatial transcriptomics",
        "skip_top": {},
    },
    "sea-ad-quantitative-neuropathology": {
        "default_modality": "neuropathology",
        "skip_top": {"simple_viewer": "web viewer assets",
                     "temp": "scratch folder"},
    },
}

# Project folders whose modality is inferred, not stated. Decision 3.
PROJECT_MODALITY = {"Multiregion_2026": "snRNA-seq"}

# Project folders whose files are cross-region by design. A file here
# that names no region in its path gets region_from = this label and is
# NOT reported missing. Its regions stay () rather than all ten codes:
# a model fitted across regions is not "MTG data", and copying it into
# ten region rows would make region totals sum to more than the bucket.
PROJECT_SCOPE = {"Multiregion_2026": "project: multi-region release"}

MODALITY_TOKENS = {"rnaseq": "snRNA-seq", "atacseq": "snATAC-seq",
                   "merfish": "spatial transcriptomics (MERFISH)",
                   "xenium": "spatial transcriptomics (Xenium)"}

JUNK_FILES = {".ds_store", "thumbs.db"}

DONOR = re.compile(r"H\d\d\.\d\d\.\d{3}")

# Extension -> format. Compound extensions are checked first.
FORMATS = [
    (".csv.gz", "csv.gz"), (".tsv.gz", "tsv.gz"), (".tar.gz", "tar.gz"),
    (".h5ad", "h5ad"), (".h5", "hdf5"), (".csv", "csv"), (".tsv", "tsv"),
    (".parquet", "parquet"), (".rds", "rds"), (".bw", "bigwig"),
    (".svs", "svs"), (".tif", "tiff"), (".tiff", "tiff"),
    (".jpeg", "jpeg"), (".jpg", "jpeg"), (".png", "png"), (".svg", "svg"),
    (".dzi", "dzi"), (".annotations", "halo annotations"),
    (".xlsx", "xlsx"), (".npy", "npy"), (".pt", "pytorch"),
    (".pkl", "pickle"),
    (".json", "json"), (".md", "markdown"), (".html", "html"),
    (".txt", "text"), (".xml", "xml"), (".zip", "zip"), (".gz", "gzip"),
]


@dataclasses.dataclass(frozen=True)
class Classification:
    bucket: str
    key: str
    unit: str            # the key itself, or the collapsed folder (ends '/')
    collapsed: bool
    regions: tuple       # canonical codes; () if none
    region_from: str     # 'folder' | 'filename' | 'project: ...' | ''
    modality: str
    modality_from: str   # 'path' | 'filename' | 'bucket default' | 'project default' | ''
    product: str
    format: str
    per_donor: bool
    deprecated: bool
    skip: str            # '' or the reason
    missing: tuple       # fields that came out 'unclassified'


# ---- helpers ------------------------------------------------------------

def file_format(name):
    low = name.lower()
    for ext, fmt in FORMATS:
        if low.endswith(ext):
            return fmt
    if "." in name.lstrip("."):
        return "unknown:" + low.rsplit(".", 1)[-1]
    return "unknown:no-extension"


def collapse_point(parts):
    """Index of the folder a key collapses into, plus its format/product.

    Returns (index, format, product) or None. The OUTERMOST match wins,
    so a Parquet file inside a Zarr store collapses into the store.
    """
    for i, seg in enumerate(parts[:-1]):
        low = seg.lower()
        fmt = None
        if low.endswith(".zarr") or low.endswith("_zarr"):
            fmt = "zarr"
        elif seg.startswith("output-XETG"):
            fmt = "xenium bundle"
        elif seg == "pseudoprogression-plots":
            fmt = "image folder"
        elif low.endswith("_files") and file_format(parts[-1]) in ("jpeg", "png"):
            fmt = "deep-zoom tiles"
        elif "=" in seg and i > 0:
            # Hive partition: collapse into the table folder above it.
            return i - 1, "partitioned parquet", COLLAPSE_PRODUCT["partitioned parquet"]
        if fmt:
            return i, fmt, COLLAPSE_PRODUCT[fmt]
    return None


# Collapsed-unit format -> product. One table, used by folder_rule() and
# by index.py reclassify (which only has the stored format to go on).
COLLAPSE_PRODUCT = {
    "zarr": "SpatialData store",
    "xenium bundle": "instrument output (Xenium)",
    "image folder": "plots",
    "deep-zoom tiles": "deep-zoom tiles",
    "partitioned parquet": "partitioned table",
}


@dataclasses.dataclass(frozen=True)
class FolderRule:
    action: str          # 'descend' | 'skip' | 'collapse'
    reason: str = ""     # skip reason
    format: str = ""     # collapse format
    product: str = ""    # collapse product


def folder_rule(bucket, prefix, child_folders=()):
    """What a folder walker should do with one folder (prefix ends '/').

    The folder-level twin of classify(): it lets the walker skip or
    collapse a folder WITHOUT listing every key inside it first.
    child_folders are the names of its immediate sub-folders, needed for
    the Hive-partition rule ('gene=A1BG' children make the parent a unit).
    tests/test_sources.py checks the two functions agree: every key
    inside a collapsed folder must classify to that folder as its unit.
    One known exception: a file sitting directly in a Hive table folder
    (no 'key=value' in its own path) cannot be recognised by classify(),
    which sees one key at a time. The walker can, so WHEN BOTH ARE
    AVAILABLE THE WALKER DECIDES UNIT MEMBERSHIP (the catalog does this).
    """
    parts = [p for p in prefix.split("/") if p]
    if not parts:
        return FolderRule("descend")
    cfg = BUCKETS.get(bucket, {"skip_top": {}})
    if parts[0] in cfg["skip_top"]:
        return FolderRule("skip", reason=cfg["skip_top"][parts[0]])
    if ".ipynb_checkpoints" in parts:
        return FolderRule("skip", reason="jupyter checkpoint copy")
    last = parts[-1]
    low = last.lower()
    fmt = None
    if low.endswith(".zarr") or low.endswith("_zarr"):
        fmt = "zarr"
    elif last.startswith("output-XETG"):
        fmt = "xenium bundle"
    elif last == "pseudoprogression-plots":
        fmt = "image folder"
    elif low.endswith("_files"):
        fmt = "deep-zoom tiles"
    elif any("=" in c for c in child_folders):
        fmt = "partitioned parquet"
    if fmt:
        return FolderRule("collapse", format=fmt, product=COLLAPSE_PRODUCT[fmt])
    return FolderRule("descend")


def region_for(vocab, parts):
    for seg in parts[:-1]:
        codes = vocab.codes_for(seg)
        if codes:
            return codes, "folder"
    stem = parts[-1].split(".", 1)[0]
    hits = []
    for tok in stem.split("_"):
        codes = vocab.codes_for(tok) if tok else None
        if not codes and "-" in tok:
            # 'Hip-MEC' -> HIP, MEC, but ONLY if every hyphen part is a
            # region. 'Higher-Order-Visual' must not become V1C via the
            # alias 'visual'.
            sub = [vocab.codes_for(p) for p in tok.split("-")]
            if all(s and len(s) == 1 for s in sub):
                codes = tuple(s[0] for s in sub)
        if codes:
            hits.extend(c for c in codes if c not in hits)
    if hits:
        return tuple(hits), "filename"
    return (), ""


def modality_for(bucket, parts, allow_defaults=True):
    for seg in parts[:-1]:
        # Whole folder name first ('RNAseq'), then its '_' tokens
        # ('xenium_output_and_alignment' -> 'xenium').
        for tok in [seg.lower()] + seg.lower().split("_"):
            m = MODALITY_TOKENS.get(tok)
            if m:
                return m, "path"
    for tok in re.split(r"[_.]", parts[-1].lower()):
        m = MODALITY_TOKENS.get(tok)
        if m:
            return m, "filename"
    if not allow_defaults:
        return "unclassified", ""
    default = BUCKETS.get(bucket, {}).get("default_modality")
    if default:
        return default, "bucket default"
    if parts[0] in PROJECT_MODALITY:
        return PROJECT_MODALITY[parts[0]], "project default"
    return "unclassified", ""


def product_for(parts, fmt):
    """What kind of thing this is. First match wins."""
    segs = {s.lower() for s in parts[:-1]}
    name = parts[-1]
    low = name.lower()

    if low in ("readme.md", "changelog.html", "index.html"):
        return "documentation"
    if low.startswith("trackdb"):
        return "track hub config"
    if "donors_objects" in segs or "donor_objects" in segs:
        return "per-donor object"
    if "supplementary information" in segs:
        return "supplementary results"
    if "model_outputs" in segs:
        return "model output"
    if "pseudobulk_objects" in segs:
        return "pseudobulk object"
    if "subclass_objects" in segs:
        return "per-subclass object"
    if "public_datasets" in segs:
        return "external dataset annotation"
    if fmt == "bigwig":
        return "coverage track"
    # neuropathology
    if fmt == "svs":
        return "whole-slide image"
    if segs & {"analysis_masks", "subregion_annotations", "halo_masked_subregions"}:
        return "analysis mask"
    if fmt == "dzi":
        return "deep-zoom tiles"
    if "downsampled" in segs and fmt == "svg":
        return "downsampled preview"
    if "downsampled" in segs and fmt == "csv":
        # Named for where it sits; contents not yet inspected.
        return "per-slide table"
    if low == "annotation.svg":
        return "annotation overlay"
    if fmt == "halo annotations" and len(parts) >= 2 and \
            low.rsplit(".", 1)[0] == parts[-2].lower():
        return "slide annotation"          # <slide>/<slide>.annotations
    if fmt == "json" and len(parts) >= 2 and DONOR.fullmatch(parts[-2]) and \
            low == parts[-2].lower() + ".json":
        return "per-donor json"            # <donor>/<donor>.json
    if "full-tabular-data" in segs:
        return "quantification table"
    # spatial
    if segs & {"combined_anndata_object", "all_donors-h5ad"}:
        return "cell-by-gene table"
    if "if_images" in segs:
        return "immunofluorescence image"
    if fmt == "tiff" and low.endswith("_max.tif"):
        return "image (max projection)"
    if "detected_transcripts" in low:
        return "transcript table"
    # single-cell datasets
    if "_metadata" in low or "cell-annotation" in low:
        return "cell metadata table"
    if re.search(r"(all|final)[-_]nuclei", low) and fmt in ("h5ad", "rds", "csv"):
        return "dataset"
    if "dream" in segs and fmt == "h5ad":
        return "dataset"
    return "unclassified"


# ---- the classifier -----------------------------------------------------

def classify(bucket, key, vocab=None):
    """Classify one S3 key. Pure function: strings in, record out."""
    vocab = vocab or regions_mod.default()
    parts = key.split("/")
    cfg = BUCKETS.get(bucket, {"default_modality": None, "skip_top": {}})

    skip = ""
    if key.endswith("/"):
        skip = "folder marker"
    elif ".ipynb_checkpoints" in parts:
        skip = "jupyter checkpoint copy"
    elif parts[-1].lower() in JUNK_FILES:
        skip = "operating-system junk file"
    elif parts[0] in cfg["skip_top"]:
        skip = cfg["skip_top"][parts[0]]

    cp = collapse_point(parts) if not skip else None
    if cp:
        i, fmt, product = cp
        unit = "/".join(parts[:i + 1]) + "/"
    else:
        unit, fmt = key, file_format(parts[-1])
        product = product_for(parts, fmt) if not skip else "skipped"

    return _describe(bucket, key, parts, unit, bool(cp), fmt, product, skip,
                     vocab)


def classify_unit(bucket, prefix, format, product, vocab=None):
    """Classify a collapsed folder unit found by the walker.

    The walker already decided format and product (from folder_rule);
    region, modality and the flags come from the folder path exactly as
    they would for a file in it.
    """
    vocab = vocab or regions_mod.default()
    parts = [p for p in prefix.split("/") if p] + [""]   # '' = no file name
    return _describe(bucket, prefix, parts, prefix, True, format, product,
                     "", vocab)


def _describe(bucket, key, parts, unit, collapsed, fmt, product, skip, vocab):
    """Region, modality, flags and 'missing' — shared by files and units."""
    # Decision 7: test/development data gets no project defaults.
    is_test = any(p.startswith("TEST_") for p in parts[:-1])
    regs, region_from = region_for(vocab, parts)
    if not regs and parts[0] in PROJECT_SCOPE and not is_test:
        region_from = PROJECT_SCOPE[parts[0]]
    modality, modality_from = modality_for(bucket, parts,
                                           allow_defaults=not is_test)
    if is_test:
        product = "test/development data"

    missing = tuple(f for f, v in (("region", regs or region_from),
                                   ("modality", modality),
                                   ("product", product), ("format", fmt))
                    if (not v) or v == "unclassified" or str(v).startswith("unknown"))
    # Documentation, test data and skipped keys legitimately have no
    # region or modality; that is 'none applies', not 'could not tell'.
    if product in ("documentation", "track hub config",
                   "test/development data") or skip:
        missing = tuple(m for m in missing if m not in ("region", "modality"))
        if modality == "unclassified":
            modality = ""     # 'none applies', not 'could not tell'

    return Classification(
        bucket=bucket, key=key, unit=unit, collapsed=collapsed,
        regions=regs, region_from=region_from,
        modality=modality, modality_from=modality_from,
        product=product, format=fmt,
        per_donor=bool(DONOR.search(key)),
        deprecated="previous_objects" in parts,
        skip=skip, missing=missing,
    )


def main(argv):
    if len(argv) < 2:
        print(__doc__.split("Every rule")[0].strip())
        print("\nusage: pathrules.py <bucket> <key> [<key> ...]")
        return 2
    bucket, keys = argv[0], argv[1:]
    for k in keys:
        c = classify(bucket, k)
        print(k)
        for f in dataclasses.fields(c):
            if f.name in ("bucket", "key"):
                continue
            print(f"    {f.name:14s} {getattr(c, f.name)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
