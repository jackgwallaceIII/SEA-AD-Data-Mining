# SEA-AD snRNA-seq pipeline

Per-nucleus QC and annotation for the [SEA-AD](https://portal.brain-map.org/explore/seattle-alzheimers-disease)
multiregion single-nucleus RNA-seq release, served through the Allen
Brain Cell (ABC) Atlas. Built around microglia in middle temporal gyrus
(MTG), but every step takes `--region` and works on any of the ten
regions in the release.

**You do not need the data to get started.** It is not in this
repository and does not need to be sent to you — it is pinned to a
versioned manifest with an md5 per file, so the fetch script rebuilds a
byte-identical copy from the Allen S3 bucket.

---

## Quickstart

```bash
git clone https://github.com/jackgwallaceIII/SEA-AD-scRNA-analysis.git
cd SEA-AD-scRNA-analysis

conda create -n seaad python=3.12 -y
conda activate seaad
pip install -r requirements-lock.txt

# Check it works before downloading anything (should be ~45 passed, 15 skipped)
python -m unittest discover -s tests -t tests

# Get the data: 9 files, 17.1 GB, md5-verified
python scripts/step0b_fetch.py preset mtg
```

Then run the pipeline:

```bash
python scripts/step0_discover.py     # provenance; no download
python scripts/step1_inventory.py    # what is in the matrix
python scripts/step2_qc_metrics.py   # ~5 min, writes results/
python scripts/step3_annotation.py   # ~4 min, writes results/
```

Output is `results/step3_annotation_MTG-10X.csv` — 1,178,694 nuclei ×
27 columns, one row per nucleus, carrying cell type, donor pathology and
QC metrics. No expression values, so it is ~300 MB and loads in pandas.

### Not ready for 17 GB?

```bash
python scripts/step0b_fetch.py preset mtg-metadata   # 1.9 GB, annotation tables only
```

Or start with the smallest region — the full pipeline on LEC is 2.05 GB
and under two minutes of compute:

```bash
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X LEC-10X/raw
python scripts/step1_inventory.py  --region LEC
python scripts/step2_qc_metrics.py --region LEC
python scripts/step3_annotation.py --region LEC
```

---

## Pulling other data

`scripts/step0b_fetch.py` browses and fetches anything in the ABC Atlas
— 43 directories, 599 files, 0.92 TB — not just this project's slice.

```bash
# Browse
python scripts/step0b_fetch.py regions                        # brain regions by name
python scripts/step0b_fetch.py list                           # all directories
python scripts/step0b_fetch.py list SEA-AD-Multiregion-10X    # files in one
python scripts/step0b_fetch.py search PFC                     # find files anywhere

# Fetch
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X MTG-10X/raw
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X --kind metadata
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X 'HIP-10X/*' --dry-run

# Named reproducible file sets
python scripts/step0b_fetch.py presets
python scripts/step0b_fetch.py preset mtg
```

### Regions by name, not code

The region codes are opaque — the prefrontal cortex is filed as `DFC`,
and the frontal insula as `FI`. `regions` prints the whole table with
anatomical names (sourced from the release's own `value_sets.csv` and
its DHBA ontology ids), sizes, and what you already have cached:

| code | region | raw |
|---|---|---|
| `AnG-10X` | angular gyrus | 2.58 GB |
| `DFC-10X` | dorsolateral prefrontal cortex | 17.76 GB |
| `FI-10X` | agranular frontal insular cortex (area FI) | 2.48 GB |
| `HIP-10X` | hippocampus | 2.41 GB |
| `ITG-10X` | inferior temporal gyrus | 5.58 GB |
| `LEC-10X` | lateral (rostral) entorhinal cortex | 2.05 GB |
| `MEC-10X` | medial (caudal) entorhinal cortex | 12.51 GB |
| `MTG-10X` | middle temporal gyrus | 15.19 GB |
| `STG-10X` | superior temporal gyrus | 7.85 GB |
| `V1C-10X` | primary visual cortex | 8.05 GB |

76.46 GB for all ten raw matrices. Anatomical names and common
abbreviations work everywhere a region is accepted, so for the
prefrontal cortex all of these are equivalent:

```bash
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X PFC/raw
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X DLPFC/raw
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X prefrontal/raw
python scripts/step0b_fetch.py get SEA-AD-Multiregion-10X DFC-10X/raw

python scripts/step2_qc_metrics.py --region PFC
```

Ambiguous words are refused rather than guessed — `temporal` matches
ITG, MTG and STG, and `entorhinal` matches MEC and LEC, so both error
with the candidates listed.

It will not let you download something expensive by accident:

- Every fetch prints the file list with sizes and a total **before**
  transferring, and marks what is already cached
- `--dry-run` prints that plan and stops
- `--max-gb` (default 50) refuses larger transfers; `--confirm-gb`
  (default 1) prompts
- A bare `get <directory>` with no narrowing is refused outright — for
  `SEA-AD-Multiregion-10X` that would be 157 GB
- md5 is verified against the manifest on every file; re-running is a
  hash check, not a re-download

Add `--yes` for unattended runs, `--force` to re-download.

---

## The data

| | |
|---|---|
| Release | `releases/20261101/manifest.json` (pinned) |
| Directory | `SEA-AD-Multiregion-10X` + `SEA-AD-Multiregion-taxonomy` |
| Matrix | `MTG-10X-raw.h5ad`, 15.19 GB, md5 `affce583f5b2db357022d289a9e9892a` |
| Shape | 1,178,694 nuclei × 36,601 genes, CSR, 6.61e9 nonzeros |
| Donors | 84, with Braak / Thal / CERAD / CPS_Global / APOE |

**The release is pinned deliberately.** A newer one usually exists; the
scripts stay on the pin and say so. This is what makes "which data did I
use" answerable later. To move, change `MANIFEST` in `scripts/seaad.py`
and `scripts/step0b_fetch.py`, then re-run every step.

**Raw, not log2.** Both ship. QC metrics need raw counts — log2 is
pre-normalized by choices that cannot be undone, and summing it gives a
meaningless total.

`scripts/step0_discover.py` documents why this release, this region and
this normalization were chosen, and reproduces the evidence in about a
minute without downloading anything.

---

## Things that will bite you

**The matrix will not fit in memory.** 6.61e9 nonzeros as int32 needs
~53 GB for data + indices. There is no chunk size or swap setting that
makes a full load work on a normal machine. Open it with `backed='r'`,
or stream it in row blocks — CSR means nuclei are contiguous on disk, so
row-block reads are cheap. `step2_qc_metrics.py` shows the pattern.

**`cell_metadata.csv` has more rows than the matrix has nuclei.** It
lists 1,861,334 for MTG against the matrix's 1,178,694. The other
682,640 (36.7%) are nuclei Allen sequenced and kept in the metadata but
excluded from the released matrix — they have `umi_count` and
`doublet_score` and nothing else: no expression data, no cluster, no
cell type. `step3_annotation.py` drops them and logs the count. The
ratio differs by region (LEC is 24.0%).

**Never join on `cell_barcode`.** It repeats across libraries —
1,033,197 unique values for 1,178,694 nuclei. Use `cell_label`, which is
barcode + sample suffix and is unique. Every merge in step 3 asserts its
row count.

**Allen already filtered.** `pct_mito < 5%` and `n_genes_detected >= 500`
were applied before release, so those distributions are truncated. The
remaining QC decisions are doublets and any tighter threshold.

**`pct_mito` is not a viability signal here.** This is single-*nucleus*
data: nuclei should contain almost no mitochondrial transcript, so an
elevated fraction means ambient or cytoplasmic carryover, not a dying
cell. Same for `pct_ribo`.

**The microglia subclass is called `Immune`,** and its supertypes are
`Micro-PVM_*`, which lumps microglia with perivascular macrophages. The
taxonomy offers no label separating them; that has to be done on markers
downstream.

**Subclass vocabularies are region-specific.** MTG and LEC both have 24
subclasses, but not the same 24 — MTG has `L4 IT` and `L6 IT Car3`, LEC
has `EC IT` and `Sub-CA1`, with 22 shared. `Immune` is shared. Any
cross-region comparison has to handle this.

---

## Layout

```
scripts/
  seaad.py              shared config: the release pin, region resolution,
                        manifest-driven paths, output naming
  step0b_fetch.py       browse and fetch from S3 (the only downloader)
  step0_discover.py     provenance: why this release/region/normalization
  step1_inventory.py    read the matrix structure before trusting it
  step2_qc_metrics.py   per-nucleus QC by streaming CSR row blocks
  step3_annotation.py   join cell type + donor + disease + QC
tests/
  test_seaad.py         path and region resolution (no data needed)
  test_fetch.py         fetch selection and guards (no data needed)
  test_outputs.py       regression tests on results/ (skipped if absent)
results/                outputs, named per region; CSVs are gitignored
abc_cache/              downloaded data; gitignored, rebuilt by step0b_fetch
```

Every script runs from any working directory and carries its reasoning
in its own docstring — including the decisions that were *rejected* and
why. Read the docstring before changing a step.

## Tests

```bash
python -m unittest discover -s tests -t tests -v
```

60 tests, about 7 seconds with data present. Tests needing `results/` or
`abc_cache/` skip automatically, so a fresh clone passes with 15 skips
and no downloads. Stdlib `unittest` only — nothing to install beyond
`requirements-lock.txt`.
