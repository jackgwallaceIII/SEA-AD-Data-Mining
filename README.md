# SEA-AD Data Mining

Find, preview, fetch and describe data from the
[Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)](https://brain-map.org/consortia/sea-ad/our-data):
single-nucleus RNA-seq and ATAC-seq, spatial transcriptomics (MERFISH,
Xenium), neuropathology slides and quantification, and the cohort's
donor, MRI and cognitive spreadsheets.

You pick a brain region and a data type. The tool shows what exists and
how big it is, downloads only what you approve, checks every file
against the server's checksum, and writes a short **structure report**
for the next person: shapes, column names, number types, and whether a
matrix holds raw counts. It never shows the data's contents.

The repository started as an RNA-seq QC and annotation pipeline on the
Allen Brain Cell (ABC) Atlas. That pipeline is still here, described
[further down](#the-rna-seq-qc-and-annotation-pipeline-abc-atlas).

---

## Quickstart

```bash
git clone https://github.com/jackgwallaceIII/SEA-AD-Data-Mining.git
cd SEA-AD-Data-Mining

conda create -n seaad python=3.12 -y
conda activate seaad
pip install -r requirements-lock.txt

python -m unittest discover -s tests -t tests
python scripts/index.py start
```

The tests need no data or network and take about 8 seconds; in a fresh
clone 19 skip, because they check RNA-seq pipeline files you have not
made or fetched yet. `start` checks your setup, shows what is already
here, and walks through the workflow with commands you can copy. **Run `start` first**; everything
below is also printed there.

Each new terminal starts in conda `base`, where the packages are not
installed. Run `conda activate seaad` before any command.

---

## The workflow: look, read, preview, download, describe

Every command is `python scripts/index.py <command>`. Shown here for LEC;
any region name or code works (see [Regions](#regions)).

**1. Look at what exists.** These read the local catalog and are instant.

```bash
python scripts/index.py show
python scripts/index.py show LEC
python scripts/index.py show LEC --files
python scripts/index.py search final-nuclei
python scripts/index.py sources
python scripts/index.py issues
```

**2. Read the region's structure report.**

```bash
python scripts/index.py report LEC
open results/reports/LEC.md
```

**3. Preview a download.** Nothing is transferred.

```bash
python scripts/index.py get LEC --product "cell metadata table" --dry-run
python scripts/index.py get --external web-donor-metadata --dry-run
```

**4. Download.** Use the same command without `--dry-run`. It prints the
plan and asks `y/N` before transferring anything.

```bash
python scripts/index.py get LEC --product "cell metadata table" --match final-nuclei
```

**5. Describe what you downloaded.**

```bash
python scripts/index.py inspect --downloaded
python scripts/index.py report --all
```

Per-file reports go to `results/inspect/<source>/<path>.md` (and `.json`);
region pages and an index go to `results/reports/`.

### Big folders and images, without downloading them

```bash
python scripts/index.py get --descriptors --match "caudate_nucleus/H20.33.001/1417642733/combined_spatialdata_zarr/" --dry-run
python scripts/index.py peek --match "1193169898/DAPI_Max.tif" --dry-run
```

- `get --descriptors` fetches only the small files that describe a Zarr
  store or a Xenium bundle: about 100 KB out of several GB.
- `peek` reads only an image's header from S3, a few hundred KB out of
  1 to 7 GB, and saves nothing.

Both preview first and ask before reading.

---

## What it covers

| Source | Access | What is in it |
|---|---|---|
| `sea-ad-single-cell-profiling` (S3) | open | snRNA-seq for 11 regions (2026 release), snATAC-seq for MTG and PFC (2024 release), per-donor objects, coverage tracks, model outputs |
| `sea-ad-spatial-transcriptomics` (S3) | open | MERFISH (MTG) and Xenium (caudate): images, SpatialData stores, transcript tables; plus one mapped cell-by-gene h5ad each for HIP and MEC |
| `sea-ad-quantitative-neuropathology` (S3) | open | whole-slide images, analysis masks, deep-zoom tile pyramids, quantification tables |
| SEA-AD website spreadsheets | open | donor metadata, MRI volumetrics, Luminex (MTG), quantitative neuropathology, pseudoprogression scores, cognitive scores, imaging availability |
| AD Knowledge Portal: `syn26223298` (raw data), `syn64410371` (Integrative Analysis Center) | controlled | linked only; this tool never downloads them |
| Allen Brain Cell Atlas | open | the RNA-seq matrices used by the pipeline below |

The three buckets hold about 9.2 TB of current data in 45,711 catalog
units (a unit is a file, or a folder treated as one item). There is
another 570 GB in deprecated `previous_objects/` folders, which are
hidden unless you ask for them.

- **The catalog** is `catalog/files/<bucket>.csv.gz`, built by listing the
  buckets. If your copy has no `catalog/files/`, build it; this lists
  the buckets and downloads nothing:
  `python scripts/index.py refresh --bucket sea-ad-single-cell-profiling`
  (repeat for the other two buckets). The single-cell bucket has about
  1.9 million objects, so the listing takes a while.
- **Website files** are listed by hand in `catalog/external_sources.csv`.
  `python scripts/index.py sources --check` asks each one for its size
  and date without downloading it.

---

## What it will and will not do

- **Nothing downloads without a preview and your yes.** That covers full
  files, header peeks and descriptor fetches. `--dry-run` only previews.
  `--max-gb` (default 50) refuses large transfers. With no terminal to
  ask, it refuses unless you pass `--yes`.
- **Controlled-access data is linked, never downloaded,** under any flags.
- **Every download is verified.** The file is written to a `.part` file,
  checked against the server's checksum and size, and only then renamed.
  A SHA-256 goes into `data_cache/_ledger.csv`, so you can always tell
  which bytes an analysis used. If a file changed on S3 after the catalog
  was built, the download stops.
- **Reports show structure, not contents.** Never shown:
  - donor IDs, barcodes and library names;
  - free text, such as clinical "describe" fields;
  - category values in files with fewer than 5 donors;
  - scanner notes in slide images;
  - Xenium run identifiers;
  - names and e-mail addresses inside panel files.
- **Read-only against S3.** Listing and header checks need no approval
  because no content moves.

Downloads go to `data_cache/`, which is gitignored.

---

## Reading the reports: what you can and cannot conclude

- **"Raw counts or normalised?"** For each matrix in an h5ad, the report
  judges this from a sample of up to a million stored values. Whole
  numbers suggest raw counts, decimals mean normalised or log values,
  negatives mean scaled data, and 0/1 means binarised. Treat it as a
  strong hint, not proof. For SEA-AD RNA objects checked so far, `X` is
  normalised and `layers["UMIs"]` holds the raw counts. Use the latter
  for QC and differential expression.
- **"Checked on 1 of 86"** means one example file was opened. It stands
  for the others only as far as they were made the same way. Region
  pages list the kinds of file nobody has opened yet.
- **Missing counts are rows, not donors.** In a per-nucleus table, one
  donor without a value adds all of their nuclei to the missing count.
- **Known upstream issues** (`python scripts/index.py issues`) are
  problems found in the data itself. Each comes with its evidence and the
  date it was found, and `show`, `get`, `inspect` and the region pages
  print them wherever they apply. Examples:
  - a track-hub file whose 96 links point to a folder that does not exist;
  - ATAC cell-type names that differ from the 2026 RNA release;
  - MERFISH images with no pixel size.

---

## Regions

Region codes, anatomical names and common abbreviations all work, e.g.
`DFC`, `PFC`, `DLPFC`, `prefrontal` and `A9`. Ambiguous words are refused
rather than guessed: `temporal` matches ITG, MTG and STG.

| code | region |
|---|---|
| AnG | angular gyrus |
| CaH | head of caudate |
| DFC | dorsolateral prefrontal cortex |
| ESOC | extrastriate occipital cortex |
| FI | agranular frontal insular cortex (area FI) |
| HIP | hippocampus |
| ITG | inferior temporal gyrus |
| LEC | lateral (rostral) entorhinal cortex |
| MEC | medial (caudal) entorhinal cortex |
| MTG | middle temporal gyrus |
| STG | superior temporal gyrus |
| V1C | primary visual cortex |

The full table, with Allen brain-atlas structure IDs, is in
`catalog/regions.csv`. Some folders cover two regions at once (HIP+MEC,
MTG+STG, V1C+ESOC, and the 10-region multi-region release). `show` and
the region pages list those separately, so a region's total counts its
own data only.

---

## Repository layout

```
scripts/
  index.py          the command line: start, show, search, get, peek, inspect,
                    report, sources, issues, refresh, reclassify, review
  start.py          the front door: setup check, local state, the workflow
  regions.py        region vocabulary (codes, names, aliases, atlas IDs)
  pathrules.py      reads region, data type and product from S3 paths
  sources.py        anonymous, read-only S3 listing
  browse.py         show / search over the catalog
  plan.py           download previews and their guards
  download.py       verified downloads and the ledger
  descriptors.py    descriptor-only fetches for Zarr stores and Xenium bundles
  peek.py           image headers over S3 range requests
  external.py       website files and controlled-access links
  inspect_data.py   structure reports: CSV, h5ad/HDF5, bigWig, TIFF/SVS,
                    deep zoom, Zarr, Xenium, Excel (images.py, stores.py,
                    xlsx.py hold the format readers)
  issues.py         known upstream issues
  report.py         one page per region
  seaad.py, step0b_fetch.py, step0-3_*.py   the RNA pipeline (below)
catalog/
  regions.csv, region_groups.csv   region vocabulary
  files/<bucket>.csv.gz            per-unit catalog, one per bucket
  catalog.csv, refresh_log.csv     summary and listing history
  known_issues.csv                 upstream problems, with evidence
  external_sources.csv             website files and Synapse studies
tests/              run without data or network
data_cache/         downloads (gitignored)
results/inspect/    per-file structure reports
results/reports/    per-region pages
abc_cache/          ABC Atlas downloads for the pipeline (gitignored)
```

Every script explains its design decisions in its opening docstring,
including what it deliberately does not do. Read it before changing a
script.

---

## The RNA-seq QC and annotation pipeline (ABC Atlas)

Per-nucleus QC and annotation for the SEA-AD multi-region snRNA-seq
release, served through the ABC Atlas. It was built around microglia in
MTG, but every step takes `--region`.

```bash
python scripts/step0b_fetch.py preset mtg
python scripts/step0_discover.py
python scripts/step1_inventory.py
python scripts/step2_qc_metrics.py
python scripts/step3_annotation.py
```

- `preset mtg` fetches 9 files, 17.1 GB, each md5-verified against a
  pinned release manifest.
- Steps 2 and 3 take about 5 and 4 minutes.
- The output, `results/step3_annotation_MTG-10X.csv`, has one row per
  nucleus (1,178,694 rows × 27 columns): cell type, donor pathology and
  QC metrics, with no expression values.

Smaller starts:

```bash
python scripts/step0b_fetch.py preset mtg-metadata
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X LEC-10X/raw
python scripts/step2_qc_metrics.py --region LEC
python scripts/step3_annotation.py --region LEC
```

- `preset mtg-metadata` fetches 1.9 GB of annotation tables only.
- The full pipeline on LEC needs 2.05 GB and under two minutes.

`step0b_fetch.py` browses and fetches anything in the ABC Atlas (`list`,
`search`, `regions`, `get`, `presets`), with the same habits as the
S3 tool: a size preview, `--dry-run`, `--max-gb`, and md5 checks against
the manifest.

**The release is pinned deliberately,** to `releases/20261101/manifest.json`.
A newer one usually exists; the scripts stay on the pin, which keeps
"which data did I use" answerable. The pin lives in `scripts/seaad.py`
and `scripts/step0b_fetch.py`; change it only on purpose, then re-run
every step.

### Things that will bite you

- **The matrix will not fit in memory.** 6.61e9 nonzeros need about
  53 GB. Open it with `backed='r'` or stream row blocks; the matrix is
  stored by row (CSR), so reading blocks of nuclei is cheap.
  `step2_qc_metrics.py` shows the pattern.
- **Use raw counts, not log2.** Both ship, but QC metrics need raw counts.
- **`cell_metadata.csv` has more rows than the matrix.** For MTG it lists
  1,861,334 nuclei against 1,178,694. The extra 36.7% were excluded from
  the released matrix and have no expression data or cell type. Step 3
  drops them and logs the count. The share differs by region (LEC 24.0%).
- **Never join on `cell_barcode`.** It repeats across libraries. Use
  `cell_label`, which is unique.
- **Allen already filtered** at `pct_mito < 5%` and `n_genes_detected >= 500`,
  so those distributions are truncated.
- **`pct_mito` is not a viability signal** in single-*nucleus* data. A
  high value means ambient or cytoplasmic carryover, not a dying cell.
- **Microglia are in the subclass `Immune`,** with supertypes
  `Micro-PVM_*` that lump microglia with perivascular macrophages.
- **Subclass lists differ by region.** MTG and LEC each have 24, of which
  22 are shared. Cross-region comparisons have to handle the rest.

---

## Tests

```bash
python -m unittest discover -s tests -t tests
```

There are about 400 tests, taking about 8 seconds, using stdlib
`unittest`. They use synthetic files and in-memory stand-ins for S3 and
the web, so they need no data or network. Tests that need the RNA-seq
pipeline's outputs (`results/*.csv`) or the cached ABC manifest
(`abc_cache/`) skip when those are absent: in a fresh clone, expect
`OK (skipped=19)`. The tests include
privacy checks that plant fake donor IDs, barcodes, free text and run
identifiers and assert they never appear in any report.
