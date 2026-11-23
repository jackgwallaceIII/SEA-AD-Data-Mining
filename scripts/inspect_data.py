"""
inspect_data.py — report the STRUCTURE of a downloaded file, never its contents.

Not a step. Run through index.py:

    python scripts/index.py inspect data_cache/sea-ad-single-cell-profiling/LEC/RNAseq/SEAAD_LEC_RNAseq_final-nuclei_metadata.2026-06-22.csv
    python scripts/index.py inspect --downloaded        # every file in the ledger

Writes results/inspect/<bucket>/<key>.json (for programs) and .md (for
the next person), and prints the .md.

(Named inspect_data.py, not inspect.py: a module called 'inspect' would
shadow Python's own standard-library module of that name.)

Decisions made here, and why:

1. THE FORMAT IS DETECTED FROM THE BYTES, NOT THE NAME.
   The first bytes are checked (Parquet 'PAR1', HDF5 signature, gzip,
   zip, pickle opcodes) before anything is parsed. A file called .csv
   that is really something else is reported as such instead of being
   fed to a CSV parser. The extension is reported alongside, and a
   disagreement is flagged.

2. WHAT IS REPORTED: SHAPE, COLUMN NAMES, TYPES, MISSING COUNTS, AND
   HOW MANY DISTINCT VALUES EACH COLUMN HAS.
   No cell values, no rows, no means/min/max. The one exception is
   decision 3.

3. CATEGORY LEVELS ARE SHOWN ONLY FOR SAFE, FEW-LEVEL TEXT COLUMNS.
   (Also withheld: free text, recognised by a prompt-like column name such
   as 'describe' / 'specify' / 'notes', or by any value over 40 characters.
   A few free-text answers in a donor table are a few people's records.)
   Per your rule: a text column with at most 30 distinct values (cell
   class, Braak stage, sex) has its levels listed, because those levels
   ARE the structure a next user needs. Withheld, with the reason given:
   - identifier-like column NAMES (donor, barcode, cell_label, sample,
     library, specimen, ... ID);
   - identifier-like VALUES even under an innocent name (SEA-AD donor
     IDs like H20.33.001; nucleotide barcodes);
   - columns where every row is unique (an identifier by definition);
   - numeric columns (their levels are measurements, not categories).
   Exception to the NAME rule: a column whose every level is a yes/no
   value (Y/N, True/False, Pass/Fail, Checked/Unchecked) is shown even
   if its name says 'donor' — 'Severely Affected Donor' is a Y/N flag,
   and a flag cannot identify anyone. The VALUE rule has no exception.
   Only the levels are listed, never how many rows have each.
   Entirely empty columns are listed at the top of the report.

4. THE FILE IS CHECKED AGAINST THE DOWNLOAD LEDGER FIRST.
   Its SHA-256 is recomputed and compared with data_cache/_ledger.csv.
   A mismatch means the file changed after it was downloaded and
   verified; the report says so prominently.

5. UNSAFE OR UNSUPPORTED FORMATS ARE REPORTED, NOT OPENED.
   Pickle files can run arbitrary code when loaded: size and SHA only,
   never opened. Parquet needs pyarrow, which is not installed (adding it
   is your decision). RDS needs R. Each says so.

6. LARGE FILES ARE STREAMED.
   CSVs are read in chunks of 200,000 rows, so a multi-GB metadata table
   never has to fit in memory. Distinct values are tracked only up to
   100,000 per column; beyond that the count is reported as "> 100,000".

7. HDF5 IS READ WITH h5py, NOT LOADED AS AN AnnData OBJECT.
   anndata.read_h5ad would pull the whole expression matrix into memory
   (tens of GB for a region). h5py reads only what is asked for: the
   matrix shapes and storage come from metadata, and obs/var are read one
   column at a time. obs and var columns go through exactly the same
   rules as CSV columns (one shared function), so a donor ID is withheld
   whichever format it arrives in. Files that are HDF5 but not AnnData
   (a 10x .h5, say) get a list of object names, shapes and types.

8. "RAW COUNTS OR NORMALISED?" IS ANSWERED FROM A SAMPLE, AS A VERDICT.
   This is the question that decides which matrix to use for QC or
   differential expression, and file names do not answer it. Up to a
   million stored values, taken from 20 places across the matrix, are
   checked for three things: negatives (scaled/centred), non-whole
   numbers (normalised or log-transformed), or whole numbers only (raw
   counts). Only the verdict is reported, never a value. The bucket README
   says X is log-normalised and raw counts are in layers; this lets the
   report confirm or contradict that for each file.

9. uns IS LISTED BY KEY NAME ONLY.
   uns holds anything the authors chose to save (colour palettes, model
   parameters, sometimes small tables). Its key names and shapes describe
   the file; its values are not read.

10. A FILE WITH FEW DONORS SHOWS NO DONOR-LEVEL LEVELS.
   Category levels are safe in a table of 34 donors: 'Female, Male' says
   nothing about anyone. In a per-donor file every donor-level column has
   ONE level, and that level is the person's own sex, APOE genotype or
   Braak stage. So when a file holds fewer than 5 donors (from its donor
   column, or a donor ID in its file name), levels are withheld for any
   column with no more distinct values than there are donors. Columns
   with more values than donors must vary from cell to cell (cell types),
   so they are still shown. Two related guards: numbers stored as text
   ('991.0') are treated as numbers (never listed), and the words 'nan'
   / 'None' are counted as missing rather than shown as a level.

11. ATAC IS NOT RNA: FEATURES MAY BE PEAKS, VALUES MAY BE 0/1.
   In an ATAC h5ad the columns (var) are usually stretches of genome
   ("chr1:10000-10500"), not genes; the report says which it found.
   Peak matrices are often binarised (open = 1), so the value verdict has
   a fourth answer: "only 0/1 values", which is not a count either.

12. bigWig COVERAGE TRACKS: HEADER ONLY.
   The header and the chromosome list sit at the start of the file and
   say which chromosomes are covered and, from chr1's length, which
   genome build the coordinates use. The signal blocks and the file's own
   summary (min/max/mean of the signal) are not read: they are data.
   Read with Python's struct module, so no new dependency is needed.

13. FOLDER FORMATS ARE DESCRIBED FROM THEIR DESCRIPTOR FILES.
   A Zarr store's consolidated metadata (zmetadata), a folder of Zarr
   metadata files, a Xenium experiment.xenium or gene_panel.json: see
   stores.py. Images (TIFF/SVS) and deep-zoom descriptors: see images.py.
"""

import csv
import hashlib
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data_cache"
OUT_DIR = ROOT / "results" / "inspect"
CATALOG_FILES = ROOT / "catalog" / "files"

LEVEL_CAP = 30           # show levels only at or below this many
TRACK_CAP = 100_000      # stop counting distinct values above this
CHUNK_ROWS = 200_000

ID_NAME = re.compile(
    r"(donor|barcode|cell[ _]?label|cell[ _]?id|sample|library|specimen|"
    r"exp_component|alignment|load_name|rna_amplification$|"
    r"(^|[ _])id($|[ _])|uuid|index$)", re.I)
DONOR_VALUE = re.compile(r"^H\d\d\.\d\d\.\d{3}")
BARCODE_VALUE = re.compile(r"^[ACGTN]{12,}(-\d+)?", re.I)
# A column whose every level is one of these cannot identify anyone,
# whatever it is called ('Severely Affected Donor' = Y/N).
FLAG_VALUES = {"y", "n", "yes", "no", "true", "false", "pass", "fail",
               "checked", "unchecked", "0", "1"}
# Free-text prompts (a clinician's "describe", "specify other") and long values:
# such text can describe one person, so it is never listed, however few values.
FREE_TEXT_NAME = re.compile(r"(describe|specify|comment|notes?\b|remark|explain|free.?text)",
                            re.I)
FREE_TEXT_LENGTH = 40
# Text that stands for "no value" when a file stored missing entries as
# words (this h5ad wrote missing as the string 'nan').
MISSING_TEXT = {"nan", "none", "null", "na", "n/a", ""}
# Below this many donors in one file, a column's levels could be a
# person's own values (a per-donor file: Sex, APOE, Braak of ONE donor).
MIN_DONORS_FOR_LEVELS = 5
# chr1:100-600, chr1-100-600, chr1_100_600, 1:100-600
PEAK_NAME = re.compile(r"^(chr)?[0-9A-Za-z_.]+[:\-_]\d+[\-_]\d+$")
BIGWIG_MAGIC_LE = (0x888FFC26).to_bytes(4, "little")
BIGBED_MAGIC_LE = (0x8789F2EB).to_bytes(4, "little")
CHROM_TREE_MAGIC = 0x78CA8C91
# chr1 length identifies the genome build without reading any signal.
CHR1_LENGTH = {248_956_422: "GRCh38 / hg38", 249_250_621: "GRCh37 / hg19",
               195_471_971: "GRCm38 / mm10", 195_154_279: "GRCm39 / mm39"}
DONOR_COLUMN = re.compile(r"^donor([ _]?id)?$", re.I)
DONOR_IN_NAME = re.compile(r"H\d\d\.\d\d\.\d{3}")


# ---- format detection -----------------------------------------------------

def sniff(path):
    """Format from the first and last bytes."""
    with open(path, "rb") as fh:
        head = fh.read(16)
        try:
            fh.seek(-4, 2)
            tail = fh.read(4)
        except OSError:
            tail = b""
    if head[:4] == b"PAR1" and tail == b"PAR1":
        return "parquet"
    if head[:8] == b"\x89HDF\r\n\x1a\n":
        return "hdf5"
    if head[:4] in (b"II*\x00", b"MM\x00*", b"II+\x00", b"MM\x00+"):
        return "tiff"                 # includes BigTIFF, Aperio .svs, OME-TIFF
    if head[:4] in (BIGWIG_MAGIC_LE, BIGWIG_MAGIC_LE[::-1]):
        return "bigwig"
    if head[:4] in (BIGBED_MAGIC_LE, BIGBED_MAGIC_LE[::-1]):
        return "bigbed"
    if head[:2] == b"\x1f\x8b":
        return "gzip"
    if head[:4] == b"PK\x03\x04":
        return "zip"
    if head[:1] == b"\x80" and head[1:2] in (b"\x02", b"\x03", b"\x04", b"\x05"):
        return "pickle"
    if head[:2] == b"X\n" or head[:5] == b"RDX2\n":
        return "rds"
    try:
        head.decode("utf-8")
        return "text"
    except UnicodeDecodeError:
        return "binary"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def ledger_entry(path, data_dir=DATA_DIR):
    import download
    want = pathlib.Path(path).resolve()
    hits = [r for r in download.read_ledger(data_dir) if r["file"].resolve() == want]
    return hits[-1] if hits else None


def shown_path(path, root=ROOT):
    """The 'file' field of a saved report: relative to the project, else the name.

    Reports are meant to be committed and shared, so they never carry an
    absolute path (which would expose a user name and home folder).
    """
    try:
        return pathlib.Path(path).resolve().relative_to(pathlib.Path(root).resolve()).as_posix()
    except ValueError:
        return pathlib.Path(path).name


def source_of(path, data_dir=DATA_DIR):
    """(bucket, key) from a data_cache path, or (None, filename)."""
    try:
        rel = pathlib.Path(path).resolve().relative_to(data_dir.resolve())
    except ValueError:
        return None, pathlib.Path(path).name
    parts = rel.parts
    return parts[0], "/".join(parts[1:])


# ---- one column, any format -----------------------------------------------

def looks_like_ids(values):
    """True if any of the first 200 values is a SEA-AD donor ID or a barcode."""
    for i, v in enumerate(values):
        if i >= 200:
            break
        v = str(v)
        if DONOR_VALUE.match(v) or BARCODE_VALUE.match(v):
            return True
    return False


def _all_numbers(values):
    """True if every value is a number written as text ('991.0', '12')."""
    seen = False
    for v in values:
        try:
            float(str(v))
        except ValueError:
            return False
        seen = True
    return seen


def describe_column(name, kind, n_rows, nulls, distinct, id_value, note=None):
    """The report entry for one column, from facts gathered by any reader.

    distinct: the distinct non-missing values, or None if there were more
    than TRACK_CAP of them. Only this function decides what is shown, so
    the withholding rules are identical for CSV, h5ad and anything later.
    """
    n_unique = None if distinct is None else len(distinct)
    entry = {"name": name, "type": kind or "empty", "missing": int(nulls),
             "distinct": n_unique if n_unique is not None else f"> {TRACK_CAP:,}",
             "levels": None, "levels_withheld": None}
    if note:
        entry["note"] = note
    entry["empty"] = n_rows > 0 and nulls == n_rows
    if entry["empty"]:
        # pandas stores an all-blank column as float64; that says
        # nothing about what the column was meant to hold. Nothing to list.
        entry["type"] = "empty"
        return entry
    is_flag = n_unique is not None and n_unique > 0 and \
        all(str(v).strip().lower() in FLAG_VALUES for v in distinct)
    reason = None
    if id_value:
        reason = "values look like donor IDs or barcodes"
    elif ID_NAME.search(name) and not is_flag:
        reason = "identifier-like column name"
    elif kind in ("text", "category") and not is_flag and (
            FREE_TEXT_NAME.search(name) or
            any(len(str(v)) > FREE_TEXT_LENGTH for v in (distinct or []))):
        reason = "free text (can describe individuals)"
    elif n_unique is not None and n_rows and n_unique == n_rows - nulls \
            and n_unique > LEVEL_CAP:
        reason = "every row unique (an identifier)"
    # A measurement stored as text ('Fresh Brain Weight' as a category of
    # '991.0') is still a measurement: treat it like a numeric column.
    if kind in ("text", "category") and not is_flag and distinct and \
            _all_numbers(distinct):
        entry["type"] = f"{kind} (numbers as text)"
        return entry
    # Only text-like columns can show levels, so only they get a
    # 'withheld' note; numeric levels are never shown at all.
    if kind in ("text", "boolean", "category"):
        if reason:
            entry["levels_withheld"] = reason
        elif n_unique is not None and n_unique <= LEVEL_CAP:
            entry["levels"] = sorted(str(v) for v in distinct)
    return entry


def donor_count(columns, key):
    """(number of donors in the file, how we know) or (None, None)."""
    for c in columns:
        if DONOR_COLUMN.match(c["name"].strip()) and isinstance(c.get("distinct"), int):
            return c["distinct"], f"column `{c['name']}`"
    if key and DONOR_IN_NAME.search(pathlib.PurePath(key).name):
        return 1, "the file name names one donor"
    return None, None


def apply_donor_floor(columns, key):
    """Withhold levels that could be the personal values of a few donors.

    In a file holding n donors, a column with at most n distinct values
    may be a donor attribute (sex, APOE, Braak), so its levels would be
    those donors' own values. A column with MORE than n values must vary
    within a donor, so it describes cells (Class, Subclass) and is safe.
    Applies only below MIN_DONORS_FOR_LEVELS donors.
    """
    n, source = donor_count(columns, key)
    if n is None or n >= MIN_DONORS_FOR_LEVELS:
        return None
    hidden = 0
    for c in columns:
        if c.get("levels") is not None and isinstance(c["distinct"], int) \
                and c["distinct"] <= n:
            c["levels"] = None
            c["levels_withheld"] = (f"only {n} donor(s) in this file: levels would be "
                                    f"personal values")
            hidden += 1
    return {"donors": n, "source": source, "withheld": hidden,
            "minimum": MIN_DONORS_FOR_LEVELS}


# ---- CSV ------------------------------------------------------------------

def _kind(series):
    import pandas as pd
    if pd.api.types.is_bool_dtype(series):
        return "boolean"
    if pd.api.types.is_integer_dtype(series):
        return "integer"
    if pd.api.types.is_float_dtype(series):
        return "float"
    return "text"


def _merge_kind(a, b):
    if a is None:
        return b
    if a == b:
        return a
    if {a, b} <= {"integer", "float"}:
        return "float"
    return "text"


def _open_text(path, compression=None):
    if compression:
        import gzip
        return gzip.open(path, "rt", newline="", encoding="utf-8", errors="replace")
    return open(path, newline="", encoding="utf-8", errors="replace")


def delimiter_for(path, compression=None):
    """',' or tab from the name (.csv / .tsv); for .txt, only if the first lines
    split consistently into 2+ fields. None means: plain text, not a table."""
    name = path.name.lower().removesuffix(".gz")
    if name.endswith(".csv"):
        return ","
    if name.endswith(".tsv"):
        return "\t"
    with _open_text(path, compression) as fh:
        lines = [ln for ln in (fh.readline() for _ in range(20)) if ln.strip()]
    for sep in ("\t", ",", ";", "|"):
        counts = {len(next(csv.reader([ln], delimiter=sep))) for ln in lines}
        if lines and len(counts) == 1 and counts.pop() >= 2:
            return sep
    return None


def count_lines(path, compression=None):
    with _open_text(path, compression) as fh:
        return sum(1 for _ in fh)


def inspect_csv(path, compression=None, sep=","):
    import pandas as pd
    with _open_text(path, compression) as fh:
        header = next(csv.reader(fh, delimiter=sep))
    dialect_note = ""
    cols = {name: {"kind": None, "nulls": 0, "distinct": set(), "overflow": False,
                   "id_value": False}
            for name in header}
    n_rows = 0
    reader = pd.read_csv(path, chunksize=CHUNK_ROWS, low_memory=False, sep=sep,
                         compression=compression or "infer")
    for chunk in reader:
        n_rows += len(chunk)
        for name in chunk.columns:
            c = cols.setdefault(name, {"kind": None, "nulls": 0, "distinct": set(),
                                       "overflow": False, "id_value": False})
            s = chunk[name]
            c["kind"] = _merge_kind(c["kind"], _kind(s))
            c["nulls"] += int(s.isna().sum())
            if not c["overflow"]:
                c["distinct"].update(s.dropna().unique().tolist())
                if len(c["distinct"]) > TRACK_CAP:
                    c["overflow"] = True
                    c["distinct"] = set()
            if not c["id_value"] and c["kind"] == "text":
                sample = s.dropna().astype(str).head(200)
                if sample.str.match(DONOR_VALUE).any() or \
                        sample.str.match(BARCODE_VALUE).any():
                    c["id_value"] = True

    columns = [describe_column(name, c["kind"], n_rows, c["nulls"],
                               None if c["overflow"] else c["distinct"], c["id_value"])
               for name, c in cols.items()]
    return {"rows": n_rows, "columns": columns, "note": dialect_note}


# ---- HDF5 / h5ad ----------------------------------------------------------

SAMPLE_VALUES = 1_000_000   # matrix values examined for 'what they look like'
SAMPLE_PIECES = 20          # ... taken from this many places across the matrix
TREE_CAP = 200              # objects listed for a non-AnnData HDF5 file
UNS_CAP = 200               # uns keys listed


def _attr(node, name, default=None):
    v = node.attrs.get(name, default)
    if isinstance(v, bytes):
        return v.decode()
    if hasattr(v, "tolist") and not isinstance(v, str):
        v = v.tolist()
        if isinstance(v, list):
            return [x.decode() if isinstance(x, bytes) else str(x) for x in v]
    return v


def _strings(ds):
    import h5py
    return ds.asstr()[()] if h5py.check_string_dtype(ds.dtype) else ds[()]


def _plain_values(node):
    """The non-missing values of a simple column, whether stored as one
    dataset (older files) or as values + mask (anndata >= 0.11 writes text
    this way: 'nullable-string-array')."""
    import h5py
    if isinstance(node, h5py.Group) and "values" in node:
        mask = node["mask"][()].astype(bool) if "mask" in node else None
        vals = _strings(node["values"])
        return vals if mask is None else vals[~mask]
    return _strings(node)


def _distinct(arr):
    """Distinct values, or None above TRACK_CAP."""
    import pandas as pd
    u = pd.unique(arr.ravel() if hasattr(arr, "ravel") else arr)
    return None if len(u) > TRACK_CAP else list(u)


def _text_missing(values):
    """Boolean mask: which text values are words for 'no value'."""
    import numpy as np
    return np.array([str(v).strip().lower() in MISSING_TEXT for v in values], dtype=bool)


NAN_NOTE = "missing values are stored as the text 'nan' (counted as missing here)"


def _read_column(node):
    """(n_rows, kind, nulls, distinct, id_value, note) for one obs/var column."""
    import h5py
    import numpy as np
    if isinstance(node, h5py.Group):
        enc = _attr(node, "encoding-type", "")
        if enc == "categorical":
            cats = node["categories"]
            codes = node["codes"][()]
            text = bool(h5py.check_string_dtype(cats.dtype))
            cat_values = _strings(cats)
            cat_missing = _text_missing(cat_values) if text else \
                np.zeros(len(cat_values), dtype=bool)
            valid = codes >= 0
            as_text = np.zeros(len(codes), dtype=bool)
            as_text[valid] = cat_missing[codes[valid]]
            used = np.unique(codes[valid & ~as_text])
            distinct = [cat_values[i] for i in used]
            return (len(codes), "category" if text else "category (numbers)",
                    int((~valid).sum() + as_text.sum()), distinct,
                    text and looks_like_ids(cat_values),
                    NAN_NOTE if as_text.any() else None)
        if enc.startswith("nullable-"):
            mask = node["mask"][()].astype(bool)
            vals = _strings(node["values"])[~mask]
            kind = {"nullable-integer": "integer", "nullable-boolean": "boolean"}.get(
                enc, "text")
            note = None
            if kind == "text":
                words = _text_missing(vals)
                note = NAN_NOTE if words.any() else None
                vals, extra = vals[~words], int(words.sum())
            else:
                extra = 0
            return (len(mask), kind, int(mask.sum()) + extra, _distinct(vals),
                    kind == "text" and looks_like_ids(vals), note)
        return (None, f"unsupported encoding ({enc or 'none'})", 0, [], False, None)
    arr = _strings(node)
    if arr.dtype.kind in "OUS":
        words = _text_missing(arr)
        vals = arr[~words]
        return (len(arr), "text", int(words.sum()), _distinct(vals), looks_like_ids(vals),
                NAN_NOTE if words.any() else None)
    if arr.dtype.kind == "b":
        return len(arr), "boolean", 0, _distinct(arr), False, None
    if arr.dtype.kind in "iu":
        return len(arr), "integer", 0, _distinct(arr), False, None
    if arr.dtype.kind == "f":
        nan = np.isnan(arr)
        return len(arr), "float", int(nan.sum()), _distinct(arr[~nan]), False, None
    return len(arr), f"other ({arr.dtype})", 0, None, False, None


def _dataframe(g, gene_axis=False):
    """obs or var: index facts plus one described entry per column."""
    import h5py
    if isinstance(g, h5py.Dataset):          # anndata < 0.7: one compound table
        return {"legacy": True,
                "columns": [{"name": n, "type": str(g.dtype[n])} for n in g.dtype.names]}
    index_name = _attr(g, "_index", "_index")
    order = _attr(g, "column-order", []) or []
    if isinstance(order, str):
        order = [order]
    out = {"legacy": False, "index": None, "columns": []}
    if index_name in g:
        n, kind, nulls, distinct, id_value, _ = _read_column(g[index_name])
        idx = {"name": index_name, "entries": n,
               "unique": distinct is not None and len(distinct) == n}
        if distinct is None:                 # too many to hold: re-check uniqueness
            import pandas as pd
            idx["unique"] = bool(pd.Index(_plain_values(g[index_name])).is_unique)
        if gene_axis:
            vals = _plain_values(g[index_name])[:1000]
            ens = sum(str(v).startswith(("ENSG", "ENSMUSG")) for v in vals)
            peaks = sum(bool(PEAK_NAME.match(str(v))) for v in vals)
            if peaks > len(vals) / 2:
                # ATAC: features are stretches of genome, not genes
                idx["looks_like"] = "genomic intervals (peaks), not genes"
            elif ens > len(vals) / 2:
                idx["looks_like"] = "Ensembl gene IDs"
            else:
                idx["looks_like"] = "gene symbols (not Ensembl IDs)"
        out["index"] = idx
    for name in order:
        n, kind, nulls, distinct, id_value, note = _read_column(g[name])
        out["columns"].append(describe_column(name, kind, n or 0, nulls, distinct,
                                              id_value, note))
    return out


def _value_sample(ds):
    """Up to SAMPLE_VALUES stored values, taken from SAMPLE_PIECES places."""
    import numpy as np
    if ds.size == 0:
        return np.array([])
    if ds.ndim == 1:
        total, per = ds.shape[0], max(1, SAMPLE_VALUES // SAMPLE_PIECES)
        starts = np.linspace(0, max(total - per, 0), SAMPLE_PIECES).astype(int)
        return np.concatenate([ds[s:s + per] for s in sorted(set(starts))])
    rows, cols = ds.shape[0], int(np.prod(ds.shape[1:]))
    per = max(1, SAMPLE_VALUES // SAMPLE_PIECES // max(cols, 1))
    starts = np.linspace(0, max(rows - per, 0), SAMPLE_PIECES).astype(int)
    return np.concatenate([ds[s:s + per].ravel() for s in sorted(set(starts))])


def _what_values_look_like(ds):
    """A verdict about the values, computed from a sample; no value is reported."""
    import numpy as np
    if ds.dtype.kind == "b":
        return "true/false"
    if ds.dtype.kind not in "iuf":
        return f"non-numeric ({ds.dtype})"
    x = _value_sample(ds)
    if x.dtype.kind == "f":
        x = x[np.isfinite(x)]
    if x.size == 0:
        return "no stored values"
    n = f"{x.size:,}"
    if (x < 0).any():
        return f"negative values present (sample of {n}): scaled or centred, not counts"
    whole = ds.dtype.kind in "iu" or bool(np.all(x == np.round(x)))
    if whole and x.max() <= 1:
        # ATAC peak matrices are often binarised: open (1) or not (0)
        return f"only 0/1 values (sample of {n}): binarised (present/absent), not counts"
    if ds.dtype.kind in "iu":
        return f"integer type, none negative (sample of {n}): counts"
    if whole:
        return f"whole numbers only (sample of {n}): looks like raw counts stored as decimals"
    return f"non-whole numbers (sample of {n}): transformed (e.g. normalised / log), not raw counts"


def _matrix(name, node, look=True):
    """Shape, storage and dtype of a dense or sparse matrix."""
    import h5py
    import numpy as np
    if isinstance(node, h5py.Dataset):
        shape = tuple(int(x) for x in node.shape)
        m = {"name": name, "storage": "dense", "shape": shape, "dtype": str(node.dtype),
             "stored": int(np.prod(shape)) if shape else 1, "density": None}
        if look and node.ndim == 2:
            m["values"] = _what_values_look_like(node)
        return m
    enc = _attr(node, "encoding-type", "")
    if enc in ("csr_matrix", "csc_matrix") and "data" in node:
        shape = tuple(int(x) for x in _attr(node, "shape"))
        nnz = int(node["data"].shape[0])
        cells = shape[0] * shape[1]
        m = {"name": name, "storage": f"{enc[:3]} sparse", "shape": shape,
             "dtype": str(node["data"].dtype), "stored": nnz,
             "density": nnz / cells if cells else None}
        if look:
            m["values"] = _what_values_look_like(node["data"])
        return m
    if enc == "dataframe":
        return {"name": name, "storage": "table", "shape": None,
                "dtype": f"{len(_attr(node, 'column-order', []) or [])} columns",
                "stored": None, "density": None}
    return {"name": name, "storage": f"unrecognised ({enc or 'group'})", "shape": None,
            "dtype": "", "stored": None, "density": None}


def _uns_keys(g, prefix="", out=None):
    """Key names under uns with their kind and shape; values never read."""
    import h5py
    out = [] if out is None else out
    for k in g:
        if len(out) >= UNS_CAP:
            break
        node, path = g[k], f"{prefix}{k}"
        if isinstance(node, h5py.Group) and _attr(node, "encoding-type", "dict") == "dict":
            out.append({"key": path, "kind": "group"})
            _uns_keys(node, path + "/", out)
        elif isinstance(node, h5py.Group):
            out.append({"key": path, "kind": _attr(node, "encoding-type", "group")})
        else:
            out.append({"key": path, "kind": "strings" if h5py.check_string_dtype(node.dtype)
                        else str(node.dtype), "shape": list(node.shape)})
    return out


def _h5_tree(f):
    """Generic listing for an HDF5 file that is not AnnData (e.g. a 10x .h5)."""
    import h5py
    out = []

    def visit(path, node):
        if len(out) >= TREE_CAP:
            return True
        if isinstance(node, h5py.Dataset):
            out.append({"path": path, "kind": "dataset", "shape": list(node.shape),
                        "dtype": str(node.dtype)})
        else:
            out.append({"path": path, "kind": "group"})
    f.visititems(visit)
    return out


def inspect_hdf5(path):
    import h5py
    with h5py.File(path, "r") as f:
        is_anndata = _attr(f, "encoding-type") == "anndata" or \
            ("X" in f and "obs" in f and "var" in f)
        if not is_anndata:
            return {"hdf5_tree": _h5_tree(f), "tree_truncated": None}, \
                "HDF5 but not AnnData: object names, shapes and types only"
        a = {"encoding": f"{_attr(f, 'encoding-type', 'anndata (pre-0.8, no tag)')} "
                         f"{_attr(f, 'encoding-version', '')}".strip(),
             "matrices": [], "obsm": [], "varm": [], "obsp": [], "varp": [], "uns": []}
        a["obs"] = _dataframe(f["obs"])
        a["var"] = _dataframe(f["var"], gene_axis=True)
        if "X" in f:
            a["matrices"].append(_matrix("X", f["X"]))
        for k in f.get("layers", {}):
            a["matrices"].append(_matrix(f"layers/{k}", f["layers"][k]))
        if "raw" in f and "X" in f["raw"]:
            a["matrices"].append(_matrix("raw/X", f["raw"]["X"]))
            if "var" in f["raw"]:
                a["raw_var"] = _dataframe(f["raw"]["var"], gene_axis=True)
        for slot in ("obsm", "varm", "obsp", "varp"):
            for k in f.get(slot, {}):
                a[slot].append(_matrix(f"{slot}/{k}", f[slot][k], look=False))
        if "uns" in f:
            a["uns"] = _uns_keys(f["uns"])
        x = next((m for m in a["matrices"] if m["name"] == "X"), None)
        n_obs = (a["obs"].get("index") or {}).get("entries")
        n_var = (a["var"].get("index") or {}).get("entries")
        if x and x["shape"]:
            n_obs, n_var = x["shape"][0], x["shape"][1]
        a["n_obs"], a["n_vars"] = n_obs, n_var
        legacy = a["obs"]["legacy"]
        return {"anndata": a}, ("inspected (AnnData, legacy layout: obs/var column names "
                                "and types only)" if legacy else "inspected (AnnData)")


# ---- bigWig / bigBed ------------------------------------------------------

def _read_chrom_tree(fh, offset, endian):
    """Chromosome names and lengths from the B+ tree; no signal is read."""
    import struct
    fh.seek(offset)
    magic, block, key_size, val_size, count, _ = struct.unpack(endian + "IIIIQQ", fh.read(32))
    if magic != CHROM_TREE_MAGIC:
        raise ValueError("chromosome tree signature not found")
    chroms = []

    def node(at):
        fh.seek(at)
        is_leaf, _, n = struct.unpack(endian + "BBH", fh.read(4))
        children = []
        for _ in range(n):
            key = fh.read(key_size).rstrip(b"\x00").decode("ascii", "replace")
            if is_leaf:
                _cid, size = struct.unpack(endian + "II", fh.read(8))
                chroms.append((key, size))
            else:
                children.append(struct.unpack(endian + "Q", fh.read(8))[0])
        for c in children:
            node(c)
    node(offset + 32)
    return chroms


def _chrom_summary(chroms):
    primary = [c for c, _ in chroms if re.fullmatch(r"(chr)?([0-9]+|X|Y|M|MT)", c)]
    return {"n": len(chroms), "primary": primary, "other": len(chroms) - len(primary)}


def inspect_bigwig(path, kind="bigwig"):
    """Header facts of a bigWig/bigBed: version, zoom levels, chromosomes,
    genome build (from chr1's length), compression. The signal itself,
    and the file's built-in summary statistics, are never read."""
    import struct
    with open(path, "rb") as fh:
        head = fh.read(64)
        endian = "<" if head[:4] in (BIGWIG_MAGIC_LE, BIGBED_MAGIC_LE) else ">"
        (_, version, zooms, chrom_tree, _data, _index, fields, defined, _autosql,
         _summary, buf, _) = struct.unpack(endian + "IHHQQQHHQQIQ", head)
        chroms = _read_chrom_tree(fh, chrom_tree, endian)
    lengths = dict(chroms)
    chr1 = lengths.get("chr1", lengths.get("1"))
    out = {"version": version, "zoom_levels": zooms, "compressed": buf > 0,
           "chromosomes": _chrom_summary(chroms),
           "genome_build": CHR1_LENGTH.get(chr1, "unknown (chr1 length not recognised)"
                                           if chr1 else "unknown (no chr1)")}
    if kind == "bigbed":
        out["bed_fields"] = {"total": fields, "standard": defined}
    return out


# ---- report ---------------------------------------------------------------

def inspect(path, data_dir=DATA_DIR, out_dir=OUT_DIR, write=True, issues_path=None):
    path = pathlib.Path(path)
    if path.is_dir():
        return inspect_zarr_folder(path, data_dir, out_dir, write, issues_path)
    if not path.is_file():
        raise SystemExit(f"\nNot a file: {path}")
    bucket, key = source_of(path, data_dir)
    fmt = sniff(path)
    ext = path.suffix.lower().lstrip(".")
    report = {"file": shown_path(path), "bucket": bucket, "key": key,
              "size": path.stat().st_size, "detected_format": fmt,
              "extension": ext, "warnings": []}

    if fmt == "pickle":
        report["status"] = "not opened: pickle files can execute code when loaded"
        report["sha256"] = sha256(path)
    else:
        report["sha256"] = sha256(path)
        led = ledger_entry(path, data_dir)
        if led is None:
            report["ledger"] = "not in the download ledger (downloaded some other way?)"
        elif led["sha256"] == report["sha256"]:
            report["ledger"] = f"matches ledger ({led['verified']}, {led['downloaded_at']})"
        else:
            report["ledger"] = "DIFFERS from the ledger"
            report["warnings"].append(
                "SHA-256 differs from the one recorded at download: the file was "
                "changed after it was verified. Re-download with --force.")

        import stores
        if fmt == "text" and path.name in ("zmetadata", ".zmetadata", "experiment.xenium",
                                           "gene_panel.json"):
            try:
                if path.name == "experiment.xenium":
                    report["xenium"] = stores.describe_xenium(path)
                elif path.name == "gene_panel.json":
                    report["panel"] = stores.describe_panel(path)
                else:
                    folder = key.rsplit("/", 1)[0] + "/" if "/" in key else ""
                    report["zarr"] = stores.describe_zarr(
                        stores.metadata_from_consolidated(path),
                        f"consolidated metadata of {folder or 'the store'}")
                report["status"] = "inspected (descriptor file; describes its folder)"
            except (ValueError, KeyError, TypeError) as exc:
                report["status"] = "descriptor file could not be interpreted"
                report["warnings"].append(f"{type(exc).__name__}: {exc}")
        elif fmt == "text" and ext in ("csv", "tsv", "txt") or fmt == "gzip" and \
                path.name.lower().endswith((".csv.gz", ".tsv.gz", ".txt.gz")):
            comp = "gzip" if fmt == "gzip" else None
            sep = delimiter_for(path, comp)
            if sep is None:
                # e.g. a genome-browser track config: text, but not a table.
                report["text_lines"] = count_lines(path, comp)
                report["status"] = (f"plain text, not a table ({report['text_lines']:,} "
                                    f"lines; contents not shown)")
            else:
                try:
                    report.update(inspect_csv(path, comp, sep))
                    report["delimiter"] = "tab" if sep == "\t" else sep
                    report["donor_floor"] = apply_donor_floor(report["columns"], key)
                    report["status"] = "inspected"
                except (ValueError, csv.Error, UnicodeDecodeError) as exc:
                    # pandas' ParserError is a ValueError. Its message gives line
                    # numbers and field counts, never cell contents.
                    report.pop("columns", None)
                    report["status"] = "could not be parsed as a table"
                    report["warnings"].append(f"{type(exc).__name__}: {exc}".strip())
        elif fmt == "zip" and ext in ("xlsx", "xlsm"):
            import xlsx
            try:
                wb = xlsx.inspect(path, describe_column, looks_like_ids, TRACK_CAP)
                for sh in wb["sheets"]:
                    sh["donor_floor"] = apply_donor_floor(sh["columns"], key)
                report["workbook"] = wb
                report["status"] = (f"inspected (Excel workbook, {len(wb['sheets'])} sheet(s); "
                                    f"read with the standard library)")
            except xlsx.WorkbookError as exc:
                report["status"] = "detected an Excel file but could not read it"
                report["warnings"].append(str(exc))
        elif fmt == "parquet":
            report["status"] = ("detected Parquet; not inspected: needs pyarrow, which "
                                "is not installed (adding it is your decision)")
        elif fmt == "hdf5":
            try:
                body, status = inspect_hdf5(path)
            except OSError as exc:
                # The signature says HDF5 but the library cannot read it:
                # truncated or damaged. Say so instead of crashing.
                body, status = {}, "detected HDF5 but it could not be opened"
                report["warnings"].append(
                    f"HDF5 signature present but the file is unreadable "
                    f"({type(exc).__name__}): truncated or damaged? Re-download "
                    f"with --force.")
            report.update(body)
            report["status"] = status
            if "anndata" in body and not body["anndata"]["obs"]["legacy"]:
                body["anndata"]["donor_floor"] = apply_donor_floor(
                    body["anndata"]["obs"]["columns"], key)
        elif fmt in ("bigwig", "bigbed"):
            try:
                report["track"] = inspect_bigwig(path, fmt)
                report["status"] = f"inspected ({fmt} header only; no signal read)"
            except Exception as exc:     # struct.error, ValueError, OSError: bad header
                report["status"] = f"detected {fmt} but its header could not be read"
                report["warnings"].append(
                    f"{fmt} header unreadable ({type(exc).__name__}): truncated or "
                    f"damaged? Re-download with --force.")
        elif fmt == "tiff":
            import images
            reader = images.LocalReader(path)
            try:
                report["image"] = images.tiff_header(reader)
                report["status"] = (f"inspected (image header only: {reader.bytes_read:,} "
                                    f"bytes read; no pixels)")
            except Exception as exc:     # HeaderError, struct.error, OSError: bad header
                report["status"] = "detected TIFF but its header could not be read"
                report["warnings"].append(f"TIFF header unreadable ({type(exc).__name__}: "
                                          f"{exc}): truncated or damaged?")
            finally:
                reader.close()
        elif fmt == "text" and ext == "dzi":
            import images
            try:
                d = images.dzi(path.read_text(encoding="utf-8", errors="replace"))
                if bucket:
                    folder = key[:-len(".dzi")] + "_files/"
                    d["tile_folder"] = folder.rsplit("/", 2)[-2] + "/"
                    n = catalog_objects(bucket, folder)
                    if n is not None:
                        d["catalog_tiles"] = n
                report["deep_zoom"] = d
                report["status"] = "inspected (deep-zoom descriptor)"
            except images.HeaderError as exc:
                report["status"] = "named .dzi but not a deep-zoom descriptor"
                report["warnings"].append(str(exc))
        elif fmt == "rds":
            report["status"] = "detected R serialized object (RDS); needs R to inspect"
        else:
            report["status"] = f"detected {fmt}; no structural reader for it yet"
        if ext == "csv" and fmt not in ("text",):
            report["warnings"].append(f"named .csv but the bytes say {fmt}")

    if bucket:
        import issues as issues_mod
        known = issues_mod.matching(issues_mod.load(issues_path),
                                    [{"bucket": bucket, "unit": key}])
        report["known_issues"] = [{k: i[k] for k in ("id", "severity", "summary",
                                                     "detail", "evidence", "found")}
                                  for i, _ in known]
    md = to_markdown(report)
    if write:
        write_report(report, md, bucket, key, out_dir)
    return report, md


def inspect_zarr_folder(path, data_dir=DATA_DIR, out_dir=OUT_DIR, write=True,
                        issues_path=None):
    """A folder of Zarr metadata files (no consolidated file)."""
    import stores
    bucket, key = source_of(path, data_dir)
    key = key.rstrip("/") + "/"
    meta = stores.metadata_from_folder(path)
    if not meta:
        raise SystemExit(f"\nNo Zarr metadata (.zgroup/.zarray/.zattrs/zarr.json) in {path}")
    report = {"file": shown_path(path), "bucket": bucket, "key": key,
              "size": sum(p.stat().st_size for p in path.rglob("*") if p.is_file()),
              "detected_format": "zarr (metadata files)", "extension": "", "warnings": [],
              "zarr": stores.describe_zarr(meta, f"{len(meta)} metadata files"),
              "status": "inspected (Zarr metadata only; no chunks read)"}
    md = to_markdown(report)
    if write:
        write_report(report, md, bucket, key.rstrip("/") + "/_store", out_dir)
    return report, md


def write_report(report, md, bucket, key, out_dir=OUT_DIR):
    """results/inspect/<bucket>/<key>.json and .md (shared with peek.py)."""
    stem = out_dir / (bucket or "local") / key
    stem.parent.mkdir(parents=True, exist_ok=True)
    stem.with_name(stem.name + ".json").write_text(json.dumps(report, indent=2))
    stem.with_name(stem.name + ".md").write_text(md)
    report["written"] = [str(stem.with_name(stem.name + ".json")),
                         str(stem.with_name(stem.name + ".md"))]


def catalog_objects(bucket, prefix, files_dir=None):
    """n_objects of a collapsed unit in the per-file catalog, or None."""
    import gzip
    f = (files_dir or CATALOG_FILES) / f"{bucket}.csv.gz"
    if not f.exists():
        return None
    with gzip.open(f, "rt", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["unit"] == prefix:
                return int(r["n_objects"] or 0)
    return None


def _human(n):
    for unit, div in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("kB", 1e3)):
        if n >= div:
            return f"{n / div:,.2f} {unit}"
    return f"{n} B"


READ_COLUMNS = [
    f"- *distinct* is the number of different values in the column; levels are "
    f"listed only for text/category columns with ≤ {LEVEL_CAP} of them.",
    "- Levels are the possible values only — how many rows carry each is not shown.",
    "- *missing* counts rows (cells), not donors: one donor with no value for a "
    "field adds all of their cells.",
    "- Identifier columns (donor, barcode, sample, library …) never show values, "
    "even when they have few."]


def _floor_md(floor):
    if not floor:
        return []
    return ["", f"**Donor-level levels withheld:** this file holds {floor['donors']} "
                f"donor(s) ({floor['source']}), fewer than {floor['minimum']}, so a column "
                f"with that few values could be a person's own data (sex, APOE, Braak …). "
                f"Levels hidden in {floor['withheld']} column(s); columns that vary between "
                f"cells (cell types) are still shown. Use the pooled region table for "
                f"levels across donors."]


def _column_table(cols):
    lines = []
    noted = [c["name"] for c in cols if c.get("note")]
    empty = [c["name"] for c in cols if c.get("empty")]
    withheld = [c["name"] for c in cols if c.get("levels_withheld")]
    if empty:
        lines += ["", f"**{len(empty)} column(s) are entirely empty** "
                      f"(every row missing): " + ", ".join(f"`{n}`" for n in empty)]
    if withheld:
        lines += ["", f"**{len(withheld)} column(s) with values withheld:** "
                      + ", ".join(f"`{n}`" for n in withheld)]
    if noted:
        lines += ["", f"**Missing values stored as the text 'nan'** in "
                      + ", ".join(f"`{n}`" for n in noted)
                      + " — counted as missing here, but pandas/scanpy will show them "
                        "as a real category called 'nan'."]
    lines += ["", "| column | type | missing | distinct | levels |",
              "|---|---|---:|---:|---|"]
    for c in cols:
        if "missing" not in c:                       # legacy table: name and type only
            lines.append(f"| {c['name']} | {c['type']} |  |  |  |")
            continue
        if c["levels"] is not None:
            lv = ", ".join(c["levels"])
        elif c["levels_withheld"]:
            lv = f"*withheld: {c['levels_withheld']}*"
        else:
            lv = ""
        d = c["distinct"] if isinstance(c["distinct"], str) else f"{c['distinct']:,}"
        lines.append(f"| {c['name']} | {c['type']} | {c['missing']:,} | {d} | {lv} |")
    return lines


def _shape(t):
    return " × ".join(f"{x:,}" for x in t) if t else ""


def _anndata_md(a):
    n_obs = f"{a['n_obs']:,}" if a["n_obs"] is not None else "?"
    n_var = f"{a['n_vars']:,}" if a["n_vars"] is not None else "?"
    lines = ["", f"## AnnData: {n_obs} cells (obs) × {n_var} features (var)",
             "", f"- **Encoding:** {a['encoding']}"]
    lines += ["", "### Matrices", "",
              "| matrix | storage | shape | dtype | stored values | density | "
              "what the values look like |", "|---|---|---|---|---:|---:|---|"]
    for m in a["matrices"]:
        dens = f"{100 * m['density']:.2f}%" if m.get("density") is not None else ""
        st = f"{m['stored']:,}" if m.get("stored") is not None else ""
        lines.append(f"| {m['name']} | {m['storage']} | {_shape(m['shape'])} | "
                     f"{m['dtype']} | {st} | {dens} | {m.get('values', '')} |")
    for label, key, gene in (("obs: one row per cell", "obs", False),
                             ("var: one row per feature (gene)", "var", True),
                             ("raw/var: features of raw/X", "raw_var", True)):
        df = a.get(key)
        if not df:
            continue
        lines += ["", f"### {label} ({len(df['columns'])} columns)"]
        if df.get("legacy"):
            lines.append("\n*Pre-0.7 AnnData layout: column names and stored types only.*")
        idx = df.get("index")
        if key == "obs":
            lines += _floor_md(a.get("donor_floor"))
        if idx:
            u = "all unique" if idx["unique"] else "NOT unique"
            extra = f"; {idx['looks_like']}" if gene and idx.get("looks_like") else ""
            lines += ["", f"Index `{idx['name']}`: {idx['entries']:,} entries, {u}"
                          f"{extra} (values not shown)."]
        if df["columns"]:
            lines += _column_table(df["columns"])
    extras = [m for slot in ("obsm", "varm", "obsp", "varp") for m in a[slot]]
    if extras:
        lines += ["", "### Embeddings and graphs (obsm, varm, obsp, varp)", "",
                  "| name | storage | shape | dtype |", "|---|---|---|---|"]
        for m in extras:
            lines.append(f"| {m['name']} | {m['storage']} | {_shape(m['shape'])} | "
                         f"{m['dtype']} |")
    if a["uns"]:
        lines += ["", f"### uns: unstructured extras ({len(a['uns'])} keys; names only)", ""]
        for u in a["uns"]:
            shape = f" {u['shape']}" if u.get("shape") else ""
            lines.append(f"- `{u['key']}` — {u['kind']}{shape}")
        if len(a["uns"]) >= UNS_CAP:
            lines.append(f"- *(listing stopped at {UNS_CAP} keys)*")
    lines += ["", "### How to read this",
              f"- *what the values look like* is judged from up to {SAMPLE_VALUES:,} "
              f"stored values taken from {SAMPLE_PIECES} places across the matrix; no "
              "value is reported. Whole numbers suggest raw counts; non-whole numbers "
              "mean a transformation (normalised, log); negatives mean scaled or "
              "centred. It is a strong hint, not proof.",
              "- Sparse matrices store only non-zero values; *density* is stored "
              "values ÷ (rows × columns).",
              "- *dtype* is the stored number type: float32 = decimals, int32/int64 = "
              "whole numbers, bool = true/false."] + READ_COLUMNS
    return lines


def to_markdown(r):
    lines = [f"# {pathlib.Path(r['file']).name}", ""]
    src = (f"https://{r['key']}" if r["bucket"] == "web" else
           f"s3://{r['bucket']}/{r['key']}" if r["bucket"] else r["file"])
    lines += [f"- **Source:** `{src}`",
              f"- **Size:** {_human(r['size'])} ({r['size']:,} bytes)",
              f"- **Format (from the bytes):** {r['detected_format']}"
              + (f" (extension .{r['extension']})" if r["extension"] else " (no extension)"),
              ] + ([f"- **SHA-256:** `{r['sha256']}`"] if r.get("sha256") else [])
    if "ledger" in r:
        lines.append(f"- **Integrity:** {r['ledger']}")
    lines.append(f"- **Status:** {r['status']}")
    for w in r["warnings"]:
        lines.append(f"- **WARNING:** {w}")
    if r.get("known_issues"):
        lines += ["", "## Known upstream issues for this file"]
        for i in r["known_issues"]:
            lines += [f"- **{i['severity'].upper()} {i['id']}:** {i['summary']}. {i['detail']}",
                      f"  *(evidence: {i['evidence']}; found {i['found']})*"]
    if "columns" in r:
        lines += ["", f"## Shape: {r['rows']:,} rows × {len(r['columns'])} columns"]
        lines += _floor_md(r.get("donor_floor"))
        lines += _column_table(r["columns"])
        lines += ["", "### How to read this"] + READ_COLUMNS + [
                  "- *type* is inferred from the text (CSV stores no types); 'text' "
                  "may hide numbers written inconsistently."]
    if "anndata" in r:
        lines += _anndata_md(r["anndata"])
    if "workbook" in r:
        for sh in r["workbook"]["sheets"]:
            above = sh.get("rows_above_header", 0)
            lines += ["", f"## Sheet '{sh['name']}'" + (" (hidden)" if sh["hidden"] else "")
                      + f": {sh['rows']:,} rows × {len(sh['columns'])} columns",
                      (f"Column names taken from row {sh['header_row']}"
                       + (f"; {above} title/group-label row(s) above it were skipped."
                          if above else ".")) if sh["header_row"] else "Empty sheet."]
            lines += _floor_md(sh.get("donor_floor"))
            if sh["columns"]:
                lines += _column_table(sh["columns"])
        lines += ["", "### How to read this"] + READ_COLUMNS + [
                  "- Excel stores dates as numbers with a display format, which is not decoded "
                  "here: a column of dates shows as float.",
                  "- Formulas are not evaluated; the value Excel saved with the file is used."]
    if "zarr" in r or "xenium" in r or "panel" in r:
        import stores
        if "zarr" in r:
            lines += stores.zarr_markdown(r["zarr"])
        if "xenium" in r:
            lines += stores.xenium_markdown(r["xenium"])
        if "panel" in r:
            lines += stores.panel_markdown(r["panel"])
    if "image" in r:
        import images
        lines += images.tiff_markdown(r["image"])
    if "deep_zoom" in r:
        import images
        lines += images.dzi_markdown(r["deep_zoom"])
    if "track" in r:
        t, c = r["track"], r["track"]["chromosomes"]
        lines += ["", "## Coverage track header",
                  f"- **Genome build (from chr1's length):** {t['genome_build']}",
                  f"- **Chromosomes:** {c['n']} — {', '.join(c['primary']) or 'none standard'}"
                  + (f" (+ {c['other']} other contigs)" if c["other"] else ""),
                  f"- **Zoom levels:** {t['zoom_levels']} (pre-computed summaries for "
                  f"genome-browser display)",
                  f"- **Compressed:** {'yes' if t['compressed'] else 'no'}; format "
                  f"version {t['version']}"]
        if "bed_fields" in t:
            lines.append(f"- **BED fields:** {t['bed_fields']['total']} "
                         f"({t['bed_fields']['standard']} standard)")
        lines += ["", "### How to read this",
                  "- A bigWig holds one number per stretch of genome (here: chromatin "
                  "accessibility summed over many cells). Only the header is read: which "
                  "chromosomes, which genome build. No signal values, and not the file's "
                  "own min/max/mean summary.",
                  "- The genome build matters when comparing with other data: coordinates "
                  "differ between hg19 and hg38, so peaks cannot be compared across builds "
                  "without conversion (liftOver)."]
    if "hdf5_tree" in r:
        lines += ["", "## HDF5 objects", "", "| path | kind | shape | dtype |",
                  "|---|---|---|---|"]
        for o in r["hdf5_tree"]:
            lines.append(f"| {o['path']} | {o['kind']} | {o.get('shape', '')} | "
                         f"{o.get('dtype', '')} |")
        if len(r["hdf5_tree"]) >= TREE_CAP:
            lines.append(f"\n*(listing stopped at {TREE_CAP} objects)*")
    return "\n".join(lines) + "\n"


def downloaded_files(data_dir=DATA_DIR):
    import download
    files = {str(r["file"]) for r in download.read_ledger(data_dir)}
    return sorted(p for p in files if pathlib.Path(p).is_file())
